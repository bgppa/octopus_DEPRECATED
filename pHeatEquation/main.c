/* Straightforward case: simples heat equation on [0,1]
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
 * G: R^basis_expansion -> R^number_of_spatial_observations_at_time_0.01
 * 
 * and our aim will be to reverse it:
 * reconstruct an approximative initial condition by observing the y datas. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "basics.h"
#include "ranvar.h"
#include "pcninv.h"
#include "fileio.h"

/* Eigenalues and eigenvectors of laplacian -dx/dx2 on [0,1]
 * with 0-boundary conditions
 * be careful with indeces: 0 in array, 1 in mathematics */
double alpha(int j)
{
        ++j;
        return pow((2. * M_PI * j), 2.);
}
double phi(int j, double x)
{
        ++j;
        return sin(2. * M_PI * j * x);
}


/* Solver is the G operator as described above (before the headers)
 * Parameters:
 - a is an array of dimension basis_expansion, which contains the first
   #basin_expansion coefficients of the initial condition given in input;
 - obs_number describes the number of space observations done at the solution
   at time 0.01. Must be strictly higher that 1,
   reflecing the idea that the boundary 0 is always observed
   (and should be 0; it would be a numerical confirm). 
 - y array of dimension obs_number that will contain the output just described.
*/
void solver(const double *a, int basis_expansion, double *y, int obs_number)
{
        assert(basis_expansion >= 1);
        assert(obs_number > 1);
        assert(a != NULL && y != NULL);

        double time_limit = 0.01;

        /* Set h, the space-step, in a way to perform a number of (equally
         * spaced) observations equal to obs_number.
         * Eg: if obs_number= 11, h is then 0.1, allowing a total of
         * 11 observations: 0, 0.1, 0.2, ... , 1. */
        double h = 1.0 / (obs_number - 1);

        double tmp_sum = 0;
        int i = 0;
        int j = 0;

        /* Having the starting condition Fourier's coefficient,
         * to compute the solution is straightforward:
         * balance them with the exponentian of eigenvalues */
        for (i = 0; i < obs_number; ++i){
                tmp_sum = 0;
                for (j = 0; j < basis_expansion; ++j){
                        /* Solution's formula */
                        tmp_sum +=
                                a[j] * exp(-alpha(j) * time_limit)
                                * phi(j, h * i);
//printf("%.3f\n", h*i);
//getchar();
                }
                y[i] = tmp_sum;
        }
}

/* Produce toy-model data to test the inversion's
 * effectiveness. It creates a random array of coefficients for
 * the starting condition. Solves the PDEs with them,
 * finally make the space observations and noise the results,
 * storing all in y.
 * Parameters:
 - noise: covariance of gaussian's noise;
 - a: array of dimension basis_expansion, elements set randomly;
 - y: array of dimension obs_number; elements set as observations. */
void createToyData(double noise, double *a, int basis_expansion,
                        double *y, int obs_number)
{
        int i=0;

        /* Randomize the parameters a */
        for (i = 0; i < basis_expansion; ++i){
                a[i] = rndmUniformIn(-10, 10);
        }

        /* Solve the PDE at fixed time, store the space observations in y */
        solver((const double *) a, basis_expansion, y, obs_number);
        printf("\n** noise-free obs: \n");
        printVec(y, obs_number);
        printf("\n");

        /* Put a noise on each value of y */
        for (i = 0; i < obs_number; ++i){
                y[i] += rndmGaussian(0, noise);
        }

        /* Note: once the toy data are constructed, the aim of the script
         * will be to reconstruct "a" by observing y.
         * Since the true "a" is available, the true error can be 
         * in this case estimated */
}       
        

int main(int argc, char *argv[]){
        srand(time(NULL));

        /* Noise used to produce the toy models data;
         * Noise introduced in the MCMC algorithm */
        double data_noise = 1e-1; 
        double mcmc_noise = 1e-1;
        
        /* The algorithm is very sensitive to the number of
         * produced samples, and how many monte carlo cycles
         * are used to produce each of it.
         * Default values: 2^10, 2^12 (powers set later) */
        int n2 = 10;
        int mcmc2 = 12;

        /* Default value for basis expansion (i.e. domain of G)
         * and number of observations (i.e. codomain of G) */
        int expansion = 3;
        int num_observations = 11;

        /* The values above can be modified via command arguments */
        if (argc >= 3){
                n2 = atoi(argv[1]);
                mcmc2 = atoi(argv[2]);
                if (argc == 5){
                /* Then also expansion and num_observations */
                        expansion = atoi(argv[3]);
                        num_observations = atoi(argv[4]);
                }
        }

                
        double *true_params = malloc(sizeof(double) * expansion);
        double *observed = malloc(sizeof(double) * num_observations);
        assert(true_params != NULL && observed != NULL);

        createToyData(data_noise, true_params, expansion,
                        observed, num_observations);

        printf("** true coeff: \n");
        printVec(true_params, expansion);
        printf("\n** noised obs: \n");
        printVec(observed, num_observations);   
        printf("\n");

        /* Now that the data are ready, set the bayes parameters */

        /* Output file where to write the posterior distribution */
        FILE *pfile = fopen("posterior_measure.txt", "w");
        assert(pfile != NULL);

        int n = (int) pow(2, n2);
        int mcmc = (int) pow(2, mcmc2);

        /* Residual error produced by the bayesian inversion */
        double err = 0;
        int i, j;
        
        /* Estimated parameters */
        double *map = malloc(sizeof(double) * expansion);
        /* Covariance matrix for the gaussian */
        double *cov = malloc(sizeof(double) * expansion * expansion);
        /* Starting point where to start the chain */
        double *start = malloc(sizeof(double) * expansion);
        assert(map != NULL && cov != NULL && start != NULL);

        /* Reset map, set a random starting point, a small covariance matrix */
        for (i = 0; i < expansion; ++i){
                map[i] = 0;
                start[i] = rndmUniformIn(-10., 10.);
                for (j = 0; j < expansion; ++j){
                        cov[i + j * expansion] = (i == j) ? 0.9 : 0.1;
                }
        }

        printf("** Starting point:\n");
        printVec(start, expansion);
        printf("\n%d samples, %d iterations per sample\n", n, mcmc);


        /* Proceed with the bayesian inversion:
         * n = number of posterior samples to produces;
         * mcmc = number of monte carlo iterations;
         * map = most frequent sample = solution = MAP
         * true_params = true known parameters (toy model data)
         * solver = the linear interpolation defined above
         * observed = vector containing the y_i
         * expansion = domain's dimension
         * observed = codomain's dimension
         * noise during the mcmc chain = mcmc_noise
         * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
         * cov = my covariance matrix, prior gaussian
         * start = starting point for the chain
         * pfile = file to write posterior distribution (values, probabilities)
         * 0 = no verbose/debug mode */
        err = bayInv(n, mcmc, map, true_params, solver, observed,
                        expansion, num_observations,
                        mcmc_noise, 0.2, cov, start, pfile, 0);

        /* err contains the residual error, i.e. G(MAP) - observations */
        /* Print the results */
        printf("MAP: ");
        printVec(map, expansion);
        printf("RES ERR: %.3f%%\n", err);
        printf("Observed output:\n");
        printVec(observed, num_observations);
        printf("MAP output :\n");
        solver(map, expansion, observed, num_observations);
        printVec(observed, num_observations);

        /* Free all the allocated memory */
        free(true_params);
        free(observed);
        free(map);
        free(cov);
        free(start);
        fclose(pfile);
        return 0;
}
