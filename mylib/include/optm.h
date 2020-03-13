#ifndef _OPTM_H_
#define _OPTM_H_

/* Find the minimum of U via gradient descent */
double gradDesc (double *x_n, int d, double (*U) (int, const double*),
                const double* lam, int iter_max, int verbose);

/* Find the minimum of U via a straighforward random search */
double rwMinimum (double* x_n, int d, double (*U) (int, const double*),
                        const double *cov, int iter_max);

double nwtnZero (double* x_n, int d, double (*U) (int, const double*),
                int iter_max, double tol);

#endif
