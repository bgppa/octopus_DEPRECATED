#ifndef _MPLS_H_
#define _MPLS_H

/* Single sample */
double uMpls (double (*U) (int, const double*), int dim,
                const double *x, int num_sampl, int num_iter, double *smpls,
		int (*okconstraint) (const double *, int));

/* Complete sample set, non parallel */
double uPcnSampler (double (*U) (int, const double*), int dim,
                const double *x, int num_sampl, int iter, 
		double *smpls, double beta,
		int (*okconstaint) (const double *, int));

/* Complete sample set, parallelized */
double prll_uPcnSampler (double (*U) (int, const double*), int dim,
			const double *x, int num_sampl,
			int iter, double *samples, double beta,
			unsigned int* seed_r,
			int (*okconstraint) (const double *, int));

double fMpls (double (*f) (int, const double*), int dim,
                const double *x, int num_sampl, int iter, double *smpls);

#endif
