/* Header pcninv.h */
/* TO DO: transfer here the comments in pcninv.c */
#ifndef _PCNINV_H_
#define _PCNINV_H_

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
                int verbose);

/* Assuming to already have a set of samples (x_i) from a measure MU,
 * compute integral(f dMU) as the sum f(x_i) successively divided by N,
 * the number of sampler. Straightforward Monte Carlo method.
 * Parameters:
 - array of N samples, each of dimension dims
   (therefore, samples is a matrix N x dims)
 - pointer to a function f: R^dims -> R */
double trivialIntegration(double *samples, int N, int dims,
                          double (*f) (const double *, int));

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
            const int verbose);
#endif
