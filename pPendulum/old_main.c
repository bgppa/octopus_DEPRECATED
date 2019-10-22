/* This is a general interface for the Bayesian inverse Problem
 * on an operator G, notation as in the README.txt file in
 * the repository root. This main file is basically, in principle,
 * the same for *every* problem:
 * what really changes is the operatior G,
 * define accordingly in an external source file:
 * here heat_eq_g.c
 * Its definition in mathematical terms follows:
 * Straightforward case: simples heat equation on [0,1]
 * with zero boundary conditions. u = u(x, t)
 *     -d^2/dx u = du/dt on [0,1] for every time t>0
 *       u(x, 0) = u_D   on [0,1]
 *       u(0, t) = u(1, t) = 0 at every time t
 *
 * The G operator described in the README file is here so defined:
 *       - express u_0, starting condition (not necessarely known), 
 *         as basis expansion
 *         by using Fourier. Say that we stop at n = 3;
 *         It's our domain dimension;
 *       - by using the approximated u_0, solve the PDE and register
 *         the results at a fixed time, say 0.01
 *       - set, as output y, various space values of u_sol,
 *         say at x=0, 0.1, ..., x=1 (again, time has been fixed).
 *
 * Summing up we have the map:
 * G: R^domain_dim -> R^number_of_spatial_observations_at_time_0.01
 * 
 * and our aim will be to reverse it:
 * reconstruct an approximative initial condition by observing the y datas. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "myblas.h"
#include "ranvar.h"
#include "pcninv.h"
#include "fileio.h"

void G(const double *x, int a, double *y, int b);
/* G is defined in pendulum */

double higherThan(const double *x, int n){
        double tmp = 0;
        G(x, n, &tmp, 1);
        return tmp > -7 ? 1 : 0;
}


/* Produce toy-model data. More precisely, it is assumed to
 * have G from R^domain_dim to R^codomain_dim.
 * The array x is initialized with random data, uniformly
 * between -10 and 10 (arbitrarely choice, no meaning).
 * Then y is created by applying G on x.
 * Finally, the values on y are noised.
 * Consequently, the aim of the main script will be to re-compute x
 * by having only y. Since the true values of x are known,
 * in this case will also be possible to compute the true error. 
 * Parameters:
 - noise: covariance of gaussian's noise;
 - x: array of dimension domain_dim, elements set randomly;
 - y: array of dimension codomain_dim; elements set as observations. */
void createToyData(double noise, double *x, int domain_dim,
                        double *y, int codomain_dim)
{
        int i = 0;

        /* Randomize the parameters x */
        for (i = 0; i < domain_dim; ++i){
                x[i] = rndmUniformIn(-1, 1, NULL);
        }

        /* Apply G to x */
        printf("%d %d\n", domain_dim, codomain_dim);
        getchar();

        G((const double *) x, domain_dim, y, codomain_dim);
        printf("\n** noise-free obs: \n");
        printVec(y, codomain_dim);
        printf("\n");

        /* Put a noise on each value of y */
        for (i = 0; i < codomain_dim; ++i){
                y[i] += rndmGaussian(0, noise, NULL);
        }
}       
        

int main(int argc, char *argv[]){
        srand(time(NULL));

        /* Noise used to produce the toy models data;
         * Noise introduced in the MCMC algorithm */
        double data_noise = 1e-1; 
        double mcmc_noise = 1e-3;
        
        /* The algorithm is very sensitive to the number of
         * produced samples, and how many monte carlo cycles
         * are used to produce each of it.
         * Default values: 2^10, 2^12 (powers set later) */
        int n2 = 10;
        int mcmc2 = 12;

        /* Default value for the domain of G
         * and its codomain */
        int domain_dim = 248;
        int num_observations = 1;

        /* The values above can be modified via command arguments */
#if 1
        if (argc >= 3){
                n2 = atoi(argv[1]);
                mcmc2 = atoi(argv[2]);
                if (argc == 5){
                /* Then also domain_dim and num_observations */
                        domain_dim = atoi(argv[3]);
                        num_observations = atoi(argv[4]);
                }
        }
#endif

                
        double *true_params = malloc(sizeof(double) * domain_dim);
        double *observed = malloc(sizeof(double) * num_observations);
        assert(true_params != NULL && observed != NULL);

//      createToyData(data_noise, true_params, domain_dim,
//                       observed, num_observations);

//      printf("** true coeff: \n");
//      printVec(true_params, domain_dim);
//      printf("\n** noised obs: \n");
//      printVec(observed, num_observations);   
//      printf("\n");

        /* Now that the data are ready, set the bayes parameters */

        /* Output file where to write the posterior distribution */
        FILE *pfile = fopen("posterior_measure.txt", "w");
        assert(pfile != NULL);

        int n = (int) pow(2, n2);
        int mcmc = (int) pow(2, mcmc2);

        unsigned int *private_seed = malloc(sizeof(unsigned int) * n);
        private_seed[0] = time(NULL);
        for (int i = 1; i < n; ++i){
                private_seed[i] = private_seed[i-1] + 1;
        }

        /* Residual error produced by the bayesian inversion */
        double err = 0;
        int i, j;
        
        /* Estimated parameters */
        double *map = malloc(sizeof(double) * domain_dim);
        /* Covariance matrix for the gaussian */
        double *cov = malloc(sizeof(double) * domain_dim * domain_dim);
        /* Starting point where to start the chain */
        double *start = malloc(sizeof(double) * domain_dim);
        assert(map != NULL && cov != NULL && start != NULL);

        /* Reset map, set a random starting point, a small covariance matrix */
        for (i = 0; i < domain_dim; ++i) {
                map[i] = 0;
                start[i] = rndmUniformIn(-1., 1., NULL);
                for (j = 0; j < domain_dim; ++j) {
                        cov[i + j * domain_dim] = (i == j) ? 1.2 : 0.2;
                }
        }


        /* integral stores the Quantity of Interest computation.
         * It's a parameter to be given to the bayes routine */
        double integral[1] = {0.};

        /* Set y manually, no tou model example */
        observed[0] = -1.;
 

//      printf("** Starting point:\n");
//      printVec(start, domain_dim);
        printf("\n%d samples, %d iterations per sample\n", n, mcmc);
        printf("Press a key to continue\n");
        getchar();

        /* Proceed with the bayesian inversion:
         * n = number of posterior samples to produces;
         * mcmc = number of monte carlo iterations;
         * map = most frequent sample = solution = MAP
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
         * higherThree is a function defined into heat.c,
           which values one iff the norm of the parameters exceedes three.
           I use it as Quantity of Interest since I'd like to compute
           the probability of having a "large" parameters;
         * intergral, here defined, just stores such a value;
         * 0 = no verbose/debug mode */
       err = bayInv(n, mcmc, map, NULL, G, observed,
                        domain_dim, num_observations,
                        mcmc_noise, 0.2, cov, start, pfile,
                        higherThan, integral, private_seed, 0);
/*
        err = bayInv(n, mcmc, map, true_params, G, observed,
                        domain_dim, num_observations,
                        mcmc_noise, 0.2, cov, start, pfile,
                        NULL, integral, 0);
*/


        /* err contains the residual error, i.e. G(MAP) - observations */
        /* Print the results */
        printf("MAP: ");
        printVec(map, domain_dim);

        FILE *mapfile = fopen("mapfile.txt", "w");
        assert(mapfile != NULL);
        fprintVec(mapfile, map, domain_dim);        
        fclose(mapfile);

        printf("RES ERR: %.3f%%\n", err);
        printf("Given observed output:\n");
        printVec(observed, num_observations);
        printf("MAP output :\n");
        G(map, domain_dim, observed, num_observations);
        printVec(observed, num_observations);
       
        printf("Probability higher that 5: %.f%%\n", integral[0]*100);

        printf("Writing the output distribution...: \n");
//        getchar();
        /* Let's try now to write an output distribution */
        FILE *outf = fopen("output_distribution.txt", "w");
        assert(mapfile != NULL);
        FILE *inpf = fopen("posterior_measure.txt", "r");
        assert(inpf != NULL);
        double fff;
        double *xxx = malloc(sizeof(double) * domain_dim);
        double *yyy = malloc(sizeof(double) * num_observations);
        int ok = 1;
        double tot_prob = 0;
        assert(ok != EOF);
        while(ok != EOF){
                for(int i = 0; i < domain_dim; ++i){
                        ok = fscanf(inpf, "%lf", xxx+i);
                }
                ok = fscanf(inpf, "%lf", &fff);
                tot_prob += fff;
//                printf("Read:\n");
//                printVec(xxx, domain_dim);
//                getchar(); 
                G(xxx, domain_dim, yyy, num_observations);
                for(int i = 0; i < num_observations; ++i){
                        fprintf(outf, "%.3f ", yyy[i]);
                }
                fprintf(outf, "%.3f \n", fff);
        }
        printf("(total prob: %.3f)\n", tot_prob);
        free(xxx);
        free(yyy);
        fclose(outf);
        fclose(inpf);

        /* Free all the allocated memory */
        free(true_params);
        free(observed);
        free(map);
        free(cov);
        free(start);
        fclose(pfile);
        free(private_seed);
        return 0;
}
