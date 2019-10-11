/* Header pcninv.h */
/* TO DO: transfer here the comments in pcninv.c */
#ifndef _PCNINV_H_
#define _PCNINV_H_

/* The following function check the validiy of parameters
 * that are supposed to be used with pcnMcmc */ 
int checkPcnParameters(double *C,
                       void (*G) (const double *, int, double *, int),
                       int ITER_NUM,
                       double *y,
                       double eta,
                       double beta,
                       int dn,
                       int dm,
                       double **tmp, 
                       double *x0,
                       double *x1);

void pcnMcmc(const double *C,
             void (*G) (const double *, int, double *, int),
             int ITER_NUM,
             const double *y,
             double eta,
             double beta,
             int dn,
             int dm,
             double **tmp,
             double *x0,
             double *x1,
             int verbose);
	 
double bayInv(int SAMPLES,
              int MCMC_ITER,
              double *MAP,
              const double *true_params,
              void (*operator) (const double *, int, double *,int),
              const double *observed_data,
              int domain_dim,
              int codomain_dim,
              double noise_var,
              double beta,
              const double *covariance_step,
              double *starting_point,
              FILE* posterior_file,
              int verbose);
#endif
