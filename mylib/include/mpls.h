#ifndef _MPLS_H_
#define _MPLS_H

double uMpls (double (*U) (int, const double*), int dim,
                const double *x, int num_sampl, int num_iter, double *smpls);

double fMpls (double (*f) (int, const double*), int dim,
                const double *x, int num_sampl, int iter, double *smpls);

double uPcnSampler (double (*U) (int, const double*), int dim,
                const double *x, int num_sampl, int iter, 
		double *smpls, double beta);

double prll_uPcnSampler (double (*U) (int, const double*), int dim,
			const double *x, int num_sampl,
			int iter, double *samples, double beta,
			unsigned int* seed_r);
#endif
