/* MCMC algorithm for the bayesian inversion problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ranvar.h"
#include "basics.h"
#include "kmeans.h"

/* Here I am trying to re-write the Bayesian Pcn inverse technique
 * in a more controlled way.
 * The pcnMcmc routine should not check for the validity of its parameters,
 * being evoked many times during the execution.
 * Rather is task is given to checkPcnParameters.
 * So now is user's responsability to use it once
 * before running the monte carlo chain.
 */

/* ??? Do I really need it? */

/* The following function checks the validiy of parameters
 * that are supposed to be used with pcnMcmc */
int checkPcnParameters(double *C, void (*G) (const double*, int, double*, int),
                       int ITER_NUM, double *y, double eta, double beta,
                       int dn, int dm, double **tmp, double *x0, double *x1)
{
        assert(C != NULL);
        assert(G != NULL);
        assert(ITER_NUM > 0);
        assert(y != NULL);
        assert(eta > 0);
        assert(beta > 0 && beta < 1);
        assert(dn > 0);
        assert(dm > 0);
        assert(tmp != NULL && *tmp != NULL);
        assert(x0 != NULL);
        assert(x1 != NULL);
        return 1;
}

/* This function performs a Monte Carlo Metropolis sampling by following
 * the pCN algorithm suitable for the Bayesian Inverse problem.
 - C : covariance matrix of the gaussian prior measure in R^dn
 - G : the operator on which we do the bayesian inversion. G : R^dn -> R^dm
 - ITER_NUM : number of steps for every MCMC
 - y : the observed points, array of dimension dm
 - eta : the noise variance
 - beta : the coefficient 0 < beta < 1 described into the pCN algorithm;
 - dn : dimension of the domain
 - dm : dimension of the codomain.
 - tmp is an array of 4 dm-dimensional arrays/pointers ALREADY initialized,
        used for some value
        schifting inside the steps (avoiding so to call malloc multiple times)
 - x0 : point in R^n on which we start the Markov Chain.
        Its value will be progressively modified until reaching
        the end of the chain;
 - x1 : point in R^n, already allocated in the memory. Used as "next point"
        in the chain. In principle is *not* an input value,
        but as tmp is more convenient to allocate one time instead
        or do it again at avery function call.
 - verbose : integer that enables a debug mode */
void pcnMcmc(const double *C, void (*G) (const double*, int, double*, int),
             int ITER_NUM, const double *y, double eta, double beta,
             int dn, int dm, double **tmp, double *x0, double *x1, int verbose)
{
        double log_alpha;
        int i = 0;
        int k = 0;
        /* I use: tmp[0], [1], [2], [3]. They
         * MUST be already initialized from main. CRUCIAL! */

        /* Perform ITER_NUM steps */
        for (i = 0; i < ITER_NUM; ++i){
                /* Key rule to keep in mind:
                 * x0 represents the previous point in every step,
                 * while x1 is the new proposal
                 * The preposed x1 follows the rules here discribed: 
                 * 1) start by sampling x1 as a 0-mean (so, NULL) 
                 * dn-dimensional gaussian
                 * with covariance matrix C. No verbose mode (0 last param) */
                rndmNdimGaussian(NULL, C, dn, x1, 0);

                /* 2) balance x1 w.r.t. previous x0  and weight beta */
                for (k = 0; k < dn; ++k){
                        x1[k] = beta * x1[k] + sqrt(1.0 - beta * beta) * x0[k];
                }

                /* Compute the potentials for Metropolis acceptance rate,
                 * whose results determine the acceptance of x1 */

                /* Put in tmp[0] the evaluation of G(x0) */
                G(x0, dn, tmp[0], dm);
                /* Copy tmp[0] into tmp[2] */
                copy(tmp[0], tmp[2], dm);
                /* Perform then: tmp[2] = y - tmp[2]
                 * i.e.: tmp[2] = y - tmp[0] = y - G(x0) as required */
                diff(y, tmp[2], dm);

                /* Put in tmp[1] the evaluation of G(x1) */
                G(x1, dn, tmp[1], dm);  
                /* Repeat the same reasoning as before, with tmp[3] */
                copy(tmp[1], tmp[3], dm);
                /* So now tmp[3] = y - G(x1) */
                diff(y, tmp[3], dm);

                /* So:  tmp[0] evaluation of G in x0
                 *      tmp[1] evaluation of G in x1
                 *      tmp[2] (componentwise) difference between y and G(x0)
                 *      tmp[3] (componentwise) difference between y and G(x1)*/

                /* Compute the logarithm of the acceptance rate alpha */
                log_alpha = (nrm2(tmp[2],dm) - nrm2(tmp[3],dm)) / (2.0*eta*eta);

                if (verbose){
                        printf("Verbose mode activated!\n");
                        printf("x0: ");
                        printVec(x0, dn);
                        printf("x1: ");
                        printVec(x1, dn);
                        printf("G(x0): ");
                        printVec(tmp[0], dm);
                        printf("G(x1): ");
                        printVec(tmp[1], dm);
                        printf("y : ");
                        printVec(y, dm);
                        printf("y - G(x0) : ");
                        printVec(tmp[2], dm);
                        printf("y - G(x1) : ");
                        printVec(tmp[3], dm);
                        printf("eta: %f\n", eta);
                        printf("|y - G(x0)|^2 : %f\n", nrm2(tmp[2],dm));
                        printf("|y - G(x1)|^2 : %f\n", nrm2(tmp[3],dm));
                        printf("log_alpha: %f\n", log_alpha);
                }

                /* Accept the new point if the rate is enough */
                if (log(rndmUniform()) <= log_alpha){
                        if (verbose){
                                printf("FROM ");
                                printVec(x0, dn);
                                printf("TO ");
                                printVec(x1, dn);
                                printf("accepted!\n");
                        }
                        /* The point is accepted: copy in x0 the value of x1,
                         * since it becomes now the new starting point.
                         * Start then the cycle again */
                        copy(x1, x0, dn);
                } /* end if log() <= log_alpha */
                else {
                        if (verbose){
                        /* Point refused: nothing to do, only verbose */
                                printf("FROM ");
                                printVec(x0, dn);
                                printf("TO ");
                                printVec(x1, dn);
                                printf("refused.\n");
                        }
                }

                if (verbose){ /* The verbose stops at every cycle */
                        getchar();
                }
        } /* End: now x1 contains a single sample from the target measure */
}

/* Try to implement an automatized Bayesian Inversion algorithm 
 * Key point: every pcnMcmc is now repeated multiple times,
 * the resulting distribution is the posterior distribution.
 * The most frequent point is returned as a solution, while, if a name
 * file is specified, this posterior distributon is written on a file.
 * Parameters:
 - SAMPLES: number of samples generated by the pcn algorithm
   described above. They will be processed with a k-means algorithm
   to determin the most frequent point. It will produce a posterior
   probability distribution with points = sqrt(SAMPLES)
 - MCMC_ITER: number of iteration for every monte carlo cycle;
 - MAP: an array of dimension domain_dim, which will contain the most
   frequent sampled point (i.e. the solution of the problem);
 - true_params: if the user known the true parameters, i.e. a toy
   problem is going to be studied, can insert them here.
   Otherwise, NULL. If non-null, the actual true error between
   true_params and MAP will be printed;
 - operator: is the R^domain_dim -> R^codomain_dim map
   to be inverted
 - observed_data: point in R^codomain_dim whose preimage is
   desired.
 - domain_dim and codomain_dim: see operator above;
 - noise_var: how intense is the noise on observed_data?
 - beta: a constant between 0 and 1 used for the pcn algorithm;
 - covariance_step: covariance of the prior gaussian distribution;
 - starting_point: where to start the search;
 - posterior_file: if not NULL, the posterior distribuion is written there */
double bayInv ( int SAMPLES, int MCMC_ITER, double *MAP,
                const double *true_params,
                void (*operator) (const double *, int, double *, int),
                const double *observed_data,
                int domain_dim,
                int codomain_dim,
                double noise_var,
                double beta,
                const double *covariance_step,
                double *starting_point,
                FILE *posterior_file,
                int verbose)
{

        int i = 0;
        
        /* Four tmp are used for technical switces during pcn monte carlo */
        double **tmp_for_mcmc = malloc(sizeof(double*) * 4);
        for (i = 0; i < 4; ++i){
                tmp_for_mcmc[i] = malloc(sizeof(double) * codomain_dim);
        }

        double *next_point = malloc(sizeof(double) * domain_dim);
        
        /* Here I'll write all the results of my pcn MCMC
         * total amount of points = SAMPLES
         * each of dimension domain_dim */
        double *posterior_points =malloc(sizeof(double) * SAMPLES * domain_dim);
        /* Omit now the check */

        /* During every iteration, the starting point will be
         * progressively modified until becoming the sampled one.
         * Consequently we need a temporar copy of its original value
         * in order to reset at the beginning of every iteration */
        double *copy_start = malloc(sizeof(double) * domain_dim);

        /* Check every allocated pointer */
        assert(posterior_points != NULL && copy_start != NULL &&
                next_point != NULL && tmp_for_mcmc != NULL);
        for (i = 0; i < 4; ++i){
                assert(tmp_for_mcmc[i] != NULL);
        }
        
        /* Sample from MCMC a number of times equal to SAMPLES */
        for (i = 0; i < SAMPLES; ++i){
                copy(starting_point, copy_start, domain_dim);
                pcnMcmc(covariance_step,
                        operator,
                        MCMC_ITER,
                        observed_data,
                        noise_var,
                        beta,
                        domain_dim,
                        codomain_dim,
                        tmp_for_mcmc,
                        copy_start,
                        next_point,     
                        verbose);

                /* Ok, now copy_start contains a single posterior sample:
                 * copy it in the right position of posterior_points */
                if (verbose){
                        printf("Sampled: \n");
                        printVec(copy_start, domain_dim);
                }
                copy(copy_start, posterior_points + i * domain_dim, domain_dim);
        }

        /* Now posterior_points contains all the sampled points,
         * i.e. is an empirical estimation of the posterior distribution */
        if (verbose){
                printf("--- %d samples have been generated --- \n", SAMPLES);
                getchar();
                printMat(posterior_points, SAMPLES, domain_dim);
                printf(" - - - - - - - - - - - - - - - - - - - - \n");
                getchar();      
        }

        /* Use now the k-means algorithm to elaborate data:
         * the actual posterior measure will contains sqrt(SAMPLES)
         * data (the square root is define into kmeans.c as a standard
         * choice). Think as kmeans ad a multidimensional way for
         * generating histograms. To every samples will be assigned its
         * frequency */

        int centroid_num = (int) sqrt(SAMPLES);
        int max_iteration_for_kmeans = 1000;
        kMean(posterior_points, SAMPLES, domain_dim, centroid_num,  
                posterior_file, max_iteration_for_kmeans, MAP);

        /* Ok, now the posterior with frequncies has been written to
         * the given file (directly with frequencies), and the MAP
         * estimator has been saved into MAP */
        if (verbose){
                printf("\nEstimated MAP: \n");
                printVec(MAP, domain_dim);
        }
        
        
        /* If the user known the true parameters, e.g. he is working with
         * toy-model data, the true relative error can be computed */
        if (true_params != NULL){
                double err = nrm2dist(MAP, true_params, domain_dim) * 100.;
                printf("ERR: %.3f%%\n", err / nrm2(true_params, domain_dim) );
        }
        
        double *MAP_output = malloc(sizeof(double) * codomain_dim);
        assert(MAP_output != NULL);
        operator(MAP, domain_dim, MAP_output, codomain_dim);

        /* So MAP_output contains the output generated by the operator when
         * the MAP, i.e. the most frequent parameters estimated with the
         * bayesian technique, are set as input. It can be used to compute
         * the residual error: (i.e. norm(MAP_output - observed_data)) */
        double res = nrm2dist(MAP_output, observed_data, codomain_dim) * 100.;
        res /= nrm2(observed_data, codomain_dim); 
        if (verbose){
                printf("RES: %.3f%%\n", res);
        }

        free(MAP_output);
        for (i = 0; i < 4; ++i){
                free(tmp_for_mcmc[i]);
        }
        free(tmp_for_mcmc);
        free(copy_start);
        free(next_point);
        free(posterior_points);
        return res; /* Return the residual error */
}
