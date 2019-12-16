/* This is a general interface for the Bayesian inverse Problem
 * on an operator G, notation as in the README.txt file in
 * the repository root. This main file is basically, in principle,
 * the same for *every* problem:
 * what really changes is the operatior G,
 * define accordingly in an external source file:


COMMENTA 

*/

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
#include "polin_g.c"

/* G is defined ad-hoc in polin_g.c */


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
                        double *y, int codomain_dim)
{
        int i = 0;
        /* Randomize the parameters x */
        for (i = 0; i < domain_dim; ++i) {
//                x[i] = i + 1;
                x[i] = rndmUniformIn(-10, 10, NULL);
        }
        /* Apply G to x */
        G((const double *) x, domain_dim, y, codomain_dim);
        printf("\n** noise-free obs: \n");
        printVec(y, codomain_dim);
        printf("\n");
        /* Put a noise on each value of y */
        for (i = 0; i < codomain_dim; ++i) {
                y[i] += rndmGaussian(0, noise, NULL);
        }
}

int main(int argc, char *argv[]) {
        srand(time(NULL));
        /* Noise used to produce the toy models data;
         * Noise introduced in the MCMC algorithm */
        double data_noise = 1e-1; 
        double mcmc_noise = 1e-1;

        /* The algorithm is very sensitive to the number of
         * produced samples, and how many monte carlo cycles
         * are used to produce each of it.
         * Default values: 2^10, 2^12 (powers set later) */
        int n = 10;
        int mcmc = 12;

        /* Default value for domain and codomain of G */
        int domain_dim = 15;
        int num_observations = 6;

        /* The values above can be modified via command arguments */
        if (argc >= 3){
                n = atoi(argv[1]);
                mcmc = atoi(argv[2]);
                if (argc == 5){
                /* Then also domain_dim and num_observations */
                        domain_dim = atoi(argv[3]);
                        num_observations = atoi(argv[4]);
                }
        }
        n = (int) pow(2, n);
        mcmc = (int) pow(2, mcmc);

        double *true_params = malloc(sizeof(double) * domain_dim);
        double *observed = malloc(sizeof(double) * num_observations);
        assert(true_params != NULL && observed != NULL);

        printf("Domain dimension: %d\n", domain_dim);
        printf("Codomain dimension: %d\n", num_observations);

        createToyData(data_noise, true_params, domain_dim,
                        observed, num_observations);
        printf("** true coeff: \n");
        printVec(true_params, domain_dim);
        printf("\n** noised obs: \n");
        printVec(observed, num_observations);   
        printf("\n");

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
                start[i] = rndmUniformIn(-10., 10., NULL);
                for (int j = 0; j < domain_dim; ++j){
                        cov[i + j * domain_dim] = (i == j) ? 0.9 : 0.1;
                }
        }
        printf("** Starting point:\n");
        printVec(start, domain_dim);
        printf("\n%d samples, %d iterations per sample\n", n, mcmc);
        printf("--- press a key to continue ---\n");
        getchar();

        /* Proceed with the bayesian inversion:
         * n = number of posterior samples to produces;
         * mcmc = number of monte carlo iterations;
         * true_params = true known parameters (toy model data)
         * G = the linear interpolation defined above
         * observed = vector containing the y_i
         * domain_dim = domain's dimension
         * observed = codomain's dimension
         * noise during the mcmc chain = mcmc_noise
         * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
         * cov = my covariance matrix, prior gaussian
         * start = starting point for the chain
         * pfile = file to write posterior distribution (values, probabilities)
         * ofile
         * intergral, here defined, just stores such a value;
         * NULL
         * 0 = no verbose/debug mode */
        double integral[1] = {0.};

        /* Create the seed for the parallelization */
        unsigned int *seed_r = malloc(sizeof(unsigned int) * n);
        seed_r[0] = time(NULL);
        for (int i = 1; i < n; ++i){
                seed_r[i] = seed_r[i-1] + 1;
        }


        bayInv(n, mcmc, true_params, G, observed,
               domain_dim, num_observations,
               mcmc_noise, 0.2, cov, start, pfile, ofile,
               NULL, integral, seed_r, 0);

        
        printf("The original data were:\n");
        printf("** True coefficients: \n");
        printVec(true_params, domain_dim);
        printf("\n** Their image under G: \n");
        printVec(observed, num_observations);   
        printf("\n");


//        printf("Expected quantity of interest: %.3f\n", integral[0]);

        /* Free all the allocated memory */
        free(true_params);
        free(observed);
        free(cov);
        free(start);
        free(seed_r);
        fclose(pfile);
        fclose(ofile);
        return 0;
}
