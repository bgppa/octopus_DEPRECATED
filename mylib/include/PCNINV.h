/* Header pcn_inversion.h */
#ifndef _PCNINV_H_
#define _PCNINV_H_



/* prior_mean and _variance are the parameters for the prior Gaussian distribution on R^dn
 * G : the operator on which we do the bayesian inversion. G : R^dn -> R^dm
 * ITER_NUM = number of MonteCarlo iterations
 * y : the observed point
 * eta : the noise variance
 * beta : the coefficient 0 < beta < 1 described into the pCN algorithm;
 * dn : dimension of the domain
 * dm : dimension of the codomain.
 * tmp is an array of pointers ALREADY initialized, used for operations in between (avoiding so to call malloc multiple times)
 * x_1 : point in R^n which contains the result of the monte carlo iterations
 * By repeating this procedure multiple times, the solution to the problem is given by the most frequent sample */

void multidim_pCN(double* prior_mean, double* prior_diag_variance, void(*G)(const double*,int,double*,int), int ITER_NUM, double* y, double eta, double beta, int dn, int dm, double**tmp, double* x_0, double* x_1);



/* ----------------- NEW PART ------------- */

/* The following function check the validiy of parameters
 * that are supposed to be used with pcnMcmc */ 
int checkPcnParameters(double* C,
                       void(*G)(const double*,int,double*,int),
                       int ITER_NUM,
                       double* y,
                       double eta,
                       double beta,
                       int dn,
                       int dm,
                       double** tmp, 
                       double* x0,
                       double* x1,
                       int verbose);


	 

#endif

