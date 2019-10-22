/* First Reinforcement Learning experiment */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include "myblas.h"
#include "ranvar.h"
#include "pcninv.h"
#include "fileio.h"
#define GOODVAL -6

/* G is defined into the rust library */
void G(const double *x, int a, double *y, int b);

double higherThan(const double *x, int n){
        double tmp = 0;
        G(x, n, &tmp, 1);
//       printf("x: \n");
//        printVec(x, n);
//        printf("G(x) = %.3f\n", tmp);
//        getchar();
        return (tmp > GOODVAL) ? 1 : 0;
}

/* Produce toy-model data. More precisely, it is assumed to
 * have G from R^domain_dim to R^codomain_dim.
 * The array x is initialized with random data, uniformly
 * between -10 and 10 (arbitrarely choice, no meaning).
 * Then y is created by applying G on x and is then perturbed by a noise.
 * So the aim of the main script will is to re-compute x having only y.
 * Since the true values of x are known, true error can be computed.
 * Parameters:
 - noise: covariance of gaussian's noise;
 - x: array of dimension domain_dim, elements set randomly;
 - y: array of dimension codomain_dim; elements set as observations. */
void createToyData(double noise, double *x, int domain_dim,
                        double *y, int codomain_dim, int verbose)
{
        int i = 0;
        /* Randomize the parameters x */
        for (i = 0; i < domain_dim; ++i) {
                x[i] = rndmUniformIn(-1, 1, NULL);
        }
        if (verbose) {
                printf("Randomly-generated true parameters: \n");
                printVec(x, domain_dim);
        }
        
        /* Apply G to x */
        G((const double *) x, domain_dim, y, codomain_dim);
        if (verbose) {
                printf("\n** noise-free obs: \n");
                printVec(y, codomain_dim);
                printf("\n");
        }
        /* Put a noise on each value of y */
        for (i = 0; i < codomain_dim; ++i) {
                y[i] += rndmGaussian(0, noise, NULL);
        }
        if (verbose) {
                printf("Noised observation: \n");
                printVec(y, codomain_dim);
                printf("\n");
        }
}

int main(int argc, char *argv[]) {
        srand(time(NULL));
        /* Noise used to produce the toy models data;
         * Noise introduced in the MCMC algorithm */
        double data_noise = 1e-1; 
        double mcmc_noise = 1e-3;

        /* The algorithm is very sensitive to the number of
         * produced samples, and how many monte carlo cycles
         * are used to produce each of it.
         * Default values: 2^10, 2^12 (powers set later) */
        int n = 10;
        int mcmc = 12;

        /* Default value for domain and codomain of G */
        int domain_dim = 248;
        int num_observations = 1;

        /* The values above can be modified via command arguments */
        if (argc >= 3){
                n = atoi(argv[1]);
                mcmc = atoi(argv[2]);
        }
        n = (int) pow(2, n);
        mcmc = (int) pow(2, mcmc);
        
        /* In this example, we want to specifically recontruct -1,
           no noy model data */
        double *true_params = NULL;
        // malloc(sizeof(double) * domain_dim);
        // assert(true_params != NULL);
        double *observed = malloc(sizeof(double) * num_observations);
        assert(observed != NULL);
        
        /* We want to find a preimage of the value -1 */
        observed[0] = -1;

        /* -- no toy model, we want to recontruct -1
        createToyData(data_noise, true_params, domain_dim,
                        observed, num_observations, 1);
        */

        /* Now that the data are ready, set the bayes parameters */
        /* Output file where to write the posterior distribution */
        FILE *pfile = fopen("posterior.txt", "w");
        assert(pfile != NULL);
        FILE *ofile = fopen("Gposterior.txt", "w");
        assert(ofile != NULL);
        
        /* Covariance matrix for the gaussian */
        double *cov = malloc(sizeof(double) * domain_dim * domain_dim);
        /* Starting point where to start the chain */
        double *start = malloc(sizeof(double) * domain_dim);
        assert(cov != NULL && start != NULL);
        /* Set a random starting point, a small covariance matrix */
        for (int i = 0; i < domain_dim; ++i){
                start[i] = rndmUniformIn(-2., 2., NULL);
                for (int j = 0; j < domain_dim; ++j){
                        cov[i + j * domain_dim] = (i == j) ? 0.9 : 0.1;
                }
        }
        printf("** Starting Markov-chain point:\n");
        printVec(start, domain_dim);
        printf("\n%d samples, %d iterations per sample\n", n, mcmc);
        printf("--- press a key to continue ---\n");
        getchar();

        double integrals[2] = {0., 0.};
        /* Create the seed for the parallelization */
        unsigned int *seed_r = malloc(sizeof(unsigned int) * n);
        seed_r[0] = time(NULL);
        for (int i = 1; i < n; ++i){
                seed_r[i] = seed_r[i-1] + 1;
        }

        /* Proceed with the bayesian inversion:
         * n = number of posterior samples to produces;
         * mcmc = number of monte carlo iterations;
         * true_params = NULL = true known parameters (toy model data)
         * G = the linear interpolation defined above
         * observed = vector containing the y_i
         * domain_dim = domain's dimension
         * observed = codomain's dimension
         * noise during the mcmc chain = mcmc_noise
         * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
         * cov = my covariance matrix, prior gaussian
         * start = starting point for the chain
         * pfile = file to write posterior distribution (values, probabilities)
         * ofile = file where to write the posterior's image
         * higherThree is a function defined into heat.c,
           which values one iff the norm of the parameters exceedes three.
           I use it as Quantity of Interest since I'd like to compute
           the probability of having a "large" parameters;
         * intergral, here defined, just stores such a value;
         * seed_r = seeds for the parallel random number generation
         * 0 = no verbose/debug mode */
        bayInv(n, mcmc, true_params, G, observed,
               domain_dim, num_observations,
               mcmc_noise, 0.2, cov, start, pfile, ofile,
               higherThan, integrals, seed_r, 0);
        
        printf("The expected quantity of interest equals:\n");
        printf("(full samples)\t\t%.3f\n", integrals[0]);
        printf("(kmeans samples)\t%.3f\n", integrals[1]); 

        /* Free all the allocated memory */
//        free(true_params);
        free(observed);
        free(cov);
        free(start);
        fclose(pfile);
        fclose(ofile);
        return 0;
}
