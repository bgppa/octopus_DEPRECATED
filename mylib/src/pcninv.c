/* MCMC algorithm for the bayesian inversion problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include "ranvar.h"
#include "myblas.h"
#include "kmeans.h"

/* Given x (point of dimension n), y (point of dimension m),
 * return the norm of y - function(x)
 * Used later in the Monte Carlo walk */
double normGdiff(const double *x, int n, const double *y, int m,
               void (*function) (const double *, int, double *, int)){

        double *tmp = malloc(sizeof(double) * m);
        assert(tmp != NULL);
        /* tmp = function(x) */
        function(x, n, tmp, m);
        /* tmp = y - tmp  = y - function(x) */
        diff(y, tmp, m);
        double result = nrm2(tmp, m);
        free(tmp);
        return result;
}


/* This function performs a Monte Carlo Metropolis sampling by following
 * the pCN algorithm suitable for the Bayesian Inverse problem.
 * All the parameters are left untpuched except for x0.
 - cov          : covariance of the gaussian prior in R^dn
 - G            : operator to invert. G : R^dn -> R^dm
 - iter         : number of Monte Carlo steps
 - y            : observed output of G, array of dimension dm
 - eta          : the noise variance
 - beta         : the coefficient 0 < beta < 1 described into the pCN algorithm;
 - dn           : dimension of the domain
 - dm           : dimension of the codomain.
 - x0           : point in R^n on which we start the Markov Chain.
                  Its value will be progressively modified during the chain,
                  and at the end will contain a Monte Carlo sample.
 - private_seed : private seed for parallels random generation.
                  if NULL, no paralelization is done.
 - verbose      : integer that enables a debug mode */
void newPcnMcmc(const double *cov,
                void (*G)(const double*, int, double*, int),
                int iter,
                const double *y,
                double eta,
                double beta,
                int dn,
                int dm,
                double *x0,
                unsigned int *private_seed,
                int verbose)
{

        double log_alpha = 0;   /* log of the acceptance rate */
        double pot0 = 0;        /* potential: norm(y - G(x0)) */
        double pot1 = 0;        /* potential: norm(y - G(x1)) */
        double *x1 = malloc(sizeof(double) * dn);
        assert(x1 != NULL);
        int i = 0;
        int k = 0;

        for (i = 0; i < iter; ++i) {
                /* In every step, x0 is the previous point and 
                 * x1 the new proposal. Defined so:
                 * start by sampling x1 as a 0-mean (so, NULL) 
                 * dn-dimensional gaussian with cov matrix.
                 * private_seed is given (NULL = no parallelization)
                 * No verbose mode (0 last param). */
                rndmNdimGaussian(NULL, cov, dn, x1, private_seed, 0);

                /* Balance x1 with x0 and beta */
                for (k = 0; k < dn; ++k) {
                        x1[k] = beta * x1[k] + sqrt(1.0 - beta * beta) * x0[k];
                }
        
                /* Compute the potentials for Metropolis acceptance rate,
                 * whose results determine the acceptance of x1 */
                pot0 = normGdiff(x0, dn, y, dm, G);
                pot1 = normGdiff(x1, dn, y, dm, G);
                log_alpha = (pot0 - pot1) / (2.0 * eta * eta);

                if (verbose) {
                        printf("x0: ");
                        printVec(x0, dn);
                        printf("x1: ");
                        printVec(x1, dn);
                        printf("y : ");
                        printVec(y, dm);
                        printf("|y - G(x0)|^2 : %f\n", pot0);
                        printf("|y - G(x1)|^2 : %f\n", pot1);
                        printf("log_alpha: %f\n", log_alpha);
                }

                /* Accept the new point if the rate is enough */
                if (log(rndmUniform(private_seed)) <= log_alpha) {
                        if (verbose){
                                printf("From ");
                                printVec(x0, dn);
                                printf("to ");
                                printVec(x1, dn);
                                printf("accepted!\n");
                                getchar();
                        }
                        /* x1 accepted: copy in x0 the value of x1,
                         * becoming the the new starting point.*/
                        copy(x1, x0, dn);
                } else { /* If rejected, just verbose */
                        if (verbose){
                                printf("From ");
                                printVec(x0, dn);
                                printf("to ");
                                printVec(x1, dn);
                                printf("refused.\n");
                                getchar();
                        }
                }
                /* Cycle again, with the possibly new starting point */
        }
        /* x0 has been rewritten with a Monte Carlo sample */
        free(x1);
}

/* Assuming to already have a set of samples (x_i) from a measure MU,
 * compute integral(f dMU) as the sum f(x_i) successively divided by N,
 * the number of samples. Straightforward Monte Carlo method.
 * Parameters:
 - array of N samples, each of dimension dims
   (therefore, samples is a matrix N x dims)
 - pointer to a function f: R^dims -> R */
double trivialIntegration(double *samples, int N, int dims,
                          double (*f) (const double *, int))
{ 
        assert(samples != NULL);
        assert(f != NULL);
        assert(N > 0 && dims > 0);
        
        double sum = 0;
        for (int i = 0; i < N; ++i) {
                sum += f(samples + i * dims, dims);
        }
        sum /= N;
        return sum;
}

/* pcn produces a single sample.
 samplePOsterior produces many of them, including the possibility
 of parallelizing the extraction WRITE BETTER */
void samplePosterior(const int samples,
                     const int iter,
                     void (*operator) (const double *, int, double *, int),
                     const double *observed_data,
                     const int dom_dim,
                     const int cod_dim,
                     const double noise_var,
                     const double beta,
                     const double *cov_step,
                     const double *start_pnt,
                     unsigned int *private_seed,
                     double *post_pts,
                     const int verbose)
{
        /* Store here all the samples produced by using
         * the pcn algorithm above (POSTerior PoinTS).
         * Their number equals the parameter "samples",
         * so malloc has dimension samples * dom_dim.
         * They are initialized with the starting point's value */
        for (int i = 0; i < samples; ++i){
                copy(start_pnt, post_pts + i * dom_dim, dom_dim);
        }
        if (private_seed == NULL){
                /* No parallelization, use the algorithm as always */
                for (int i = 0; i < samples; ++i){
                        printf("...sampling %d of %d\n", i+1, samples);
                        newPcnMcmc(cov_step, operator,
                                   iter, observed_data,
                                   noise_var, beta,
                                   dom_dim,cod_dim,
                                   post_pts + i * dom_dim, NULL,
                                   verbose);
                }
        } else {
                #pragma omp parallel for
                for (int i = 0; i < samples; ++i) {
                        printf("...(thread %d): sampling %d of %d\n",
                                omp_get_thread_num(), i+1, samples);
                        newPcnMcmc(cov_step, operator,
                                   iter, observed_data,
                                   noise_var, beta,
                                   dom_dim, cod_dim,
                                   post_pts + i * dom_dim, private_seed + i,
                                   verbose);
                }
        }
}

/* Automatized Bayesian Inversion routine. 
 * It performs pcnMcmc multiple times in a way to produce many samples, so 
 * the resulting distribution is an approximation of the posterior distribution.
 * This collection of samples is used to compute the integral of a possible
 * quantity of interest. Then it is reduced in dimensionality via
 * multidimensional histrogram (k-means algorithm), and the most frequent
 * point (MAP) is printed; it's our "solution". 
 * Parameters:
 - samples      : number of samples generated by the pcn algorithm
 - iter         : number of iteration for every monte carlo cycle
 - true_params  : if the user knowns the true parameters, i.e. a toy
                  problem is going to be studied, can insert them here.
                  Otherwise, NULL. If non-null, the actual true error between
                  true_params and map will be printed
 - operator     : is the R^dom_dim -> R^cod_dim map to inverted
 - observed_data: output, point in R^cod_dim, whose preimage is desired.
 - dom_dim      : dimension of operator's domain
 - cod_dim      : dimension of operators's codomain (and observed_data)
 - noise_var    : how intense is the noise on observed_data?
 - beta         : a constant between 0 and 1 used for the pcn algorithm;
 - cov_step     : covariance of the prior gaussian distribution;
 - start_pnt    : where to start the Monte Carlo Chain (in pcnInv);
 - post_file    : if not NULL, FILE where the posterior distribuion is written
 - Gpost_file   : file when the posterior's dist IMAGE under G is written
 - qoi          : Quantity of interest: function R^dom_dim -> R to
                  integrate wrt the posterior measure. Can be set to NULL.
 - intgrted_qoi : pointer to a double that will contain the result of
                  the integration above. Can be set to NULL when not used.
 - private_seed : when non NULL, allows the sampling to be computed in parallel
                  by using rand_r its initialization whose rules are not
                  repeated here.
 - verbose      : when positive, print more messages */
void bayInv(const int samples,
            const int iter,
            const double *true_params,
            void (*operator) (const double *, int, double *, int),
            const double *observed_data,
            const int dom_dim,
            const int cod_dim,
            const double noise_var,
            const double beta,
            const double *cov_step,
            const double *start_pnt,
            FILE *post_file,
            FILE *Gpost_file,
            double (*qoi) (const double *x, int dim),
            double *intgrted_qoi,
            unsigned int *private_seed,
            const int verbose)
{
        /* 1. Sample the posterior distribution, storing into post_pts */
        double *post_pts = malloc(sizeof(double) * samples * dom_dim);
        assert(post_pts != NULL);
        samplePosterior(samples, iter, operator,
                        observed_data, dom_dim, cod_dim,
                        noise_var, beta, cov_step,
                        start_pnt, private_seed, post_pts, verbose);
        if (verbose) {
                printf("--- %d samples generated --- \n", samples);
                printMat(post_pts, samples, dom_dim);
                getchar();      
                printf("(press a key to continue)\n");
       }

        /* 2. Integrate the Quantity of Interest */
        if (qoi != NULL && intgrted_qoi != NULL) {
                *intgrted_qoi = trivialIntegration(post_pts, samples,
                                                     dom_dim, qoi);
        } else {
                printf("*Remark: no Quantity of Interest to integrate*.\n");
        }

        /* 3. Reduce the posterior measure by using the k-means algorithm
         * The reduced posterior distribution is stored into post_reduced */
        int clusters = (int) sqrt(samples);
        int max_iteration_for_kmeans = 2000;
        double *post_reduced = malloc(sizeof(double) * clusters * (dom_dim+1));
        assert(post_reduced != NULL);
        kMeans(post_pts, samples, dom_dim, clusters,  
               max_iteration_for_kmeans, post_reduced, verbose);
        /* Write the reduced posterior distribution to post_file file */
        if (post_file != NULL) {
                fprintMat(post_file, post_reduced, clusters, dom_dim + 1); 
        }

        /* 4. Compute the image of the reduced posterior distribution,
         * write it possibly to the file Gpost_file */
        double *post_output = malloc(sizeof(double) * clusters * (cod_dim + 1));        assert(post_output != NULL);
        for (int i = 0; i < clusters; ++i) {
                /* The first element of each row is the same as the first
                 * element of every post_reduced's row, i.e. the %frequency */
                post_output[i * (cod_dim + 1)]=post_reduced[i * (dom_dim + 1)];
                /* while the remaining are the post_reduced's images */
                operator(post_reduced + (i * dom_dim) + 1, dom_dim,
                         post_output + (i * (cod_dim + 1)) + 1, cod_dim);
        }
        if (Gpost_file != NULL) {
                fprintMat(Gpost_file, post_output, clusters ,cod_dim + 1);
        }

        /* 5. Error estimation and MAP computation */
        printf("MAP:\n");
        printVec(post_reduced + 1, dom_dim);
        printf("Its image under G:\n");
        printVec(post_output + 1, cod_dim);
        double err = 0;
        /* If the user knowns the true parameters, e.g. he is working with
         * toy-model data, the true relative error is computed */
        if (true_params != NULL) {
                err = nrm2dist(post_reduced + 1, true_params, dom_dim) * 100.;
                printf("ERR: %.3f%%\n", err / nrm2(true_params, dom_dim));
        }
        /* anyay, we compute the residuom (output's discrepance) */ 
        err = nrm2dist(post_output + 1, observed_data, cod_dim) * 100.;
        err /= nrm2(observed_data, cod_dim); 
        printf("RES: %.3f%%\n", err);

        free(post_reduced);
        free(post_pts);
        free(post_output);
}
