/* Metropolis, a module implementig various naive Metropolis Markov
 * Chain Monte Carlo approaches */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include "myblas.h"
#include "ranvar.h"
#include "mpls.h"

/* Enable the parallel support for the Metropolis module.
 * When ENABLED, the only problem will be a false positive error
 * when using valgrind. Therefore disabling it can be useful
 * for temporary debugging processes */
#define MPLS_PARALL 1

double uMpls(double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, double *smpls)
{
	/* if U is a dim-dimensional R^dim -> R potential function,
	 * we aim at sampling the density exp(-U(x))dx via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */

	double alpha = 0;
	double accepted = 0; /* Counters for estimating the acceptance rate */
	double mean_accepted = 0;
	/* x0 and x1 are the previous and nex points in the chain */
	double *x0 = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
	double *cov = malloc(sizeof(double) * dim * dim);
	id(cov, dim);

	/* Produce N samples, each stored in smpls */
	for (int ii = 0; ii < num_sampl; ++ii) {
		copy(x, x0, dim);	/* Reset the starting value */
		accepted = 0.;		/* Reset the acceptance counter */
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* Gaussian proposal: step from x0 */
			rndmNdimGaussian(x0, cov, dim, x1, NULL, 0);
			alpha = exp(U(dim, x0) - U(dim, x1));
			if (alpha > 1.) { alpha = 1.; }
			/* Determine if accepting the new peoposed step */
			if (rndmUniform(NULL) <= alpha) {
				copy(x1, x0, dim);
				++accepted;
			}
		}
		/* A sample of the distribution is in x0: copy to smpls */
		copy(x0, smpls + (ii * dim), dim);
		accepted /= iter;
		mean_accepted += accepted;
	}
	free(x0);
	free(x1);
	free(cov);
	return mean_accepted / num_sampl * 100.;
}

/* Pcn-Metropolis algorithm for potentials U */
double prll_uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta, unsigned int *seed_r)
{
	/* if U is a dim-dimensional R^dim -> R potential function,
	 * we aim at sampling the density exp(-U(x))dx via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */
	setbuf(stdout, NULL);
#if! MPLS_PARALL
	printf("\n\n*warning: prll_uPcnSampler: calling a parallel function,"
		" but the support is DISABLED. Check the"
	        " constant MPLS_PARALL in mpls.c\n");
	printf("*It should work anyway, but with no benefits*\n");
#endif
	double *log_alpha = malloc(sizeof(double) * num_sampl);
	double *accepted = malloc(sizeof(double) * num_sampl);
	fillzero(accepted, num_sampl);
		 /* Counters for estimating the acceptance rate */
	double mean_accepted = 0;
	/* x0 and x1 are the previous and nex points in the chain */
	double *x1n = malloc(sizeof(double) * dim * num_sampl);
	fillzero(x1n, dim * num_sampl);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
	double *cov = malloc(sizeof(double) * dim * dim);
	id(cov, dim);
	double *allzero = malloc(sizeof(double) * dim);
	fillzero(allzero, dim);
	double *tmpn = malloc(sizeof(double) * dim * num_sampl);
	fillzero(tmpn, dim * num_sampl);
	/* Fill smpls with multiple copies of x */
	for (int i = 0; i < num_sampl; ++i) {
		copy(x, smpls + i * dim, dim);
	}
	/* Produce N samples, each stored in smpls */
#if MPLS_PARALL
	#pragma omp parallel for
#endif
	for (int n = 0; n < num_sampl; ++n) {
		if (n % 10 == 0) {
			printf(".");
		}
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* Sample the gaussian in tmp */
			rndmNdimGaussian(allzero, cov, dim,
				       	tmpn + n * dim, seed_r + n, 0);
			/* Propose x1 as the weightes sum between that gaussian
			 * and the previous x0, which is smpls[n*dim] */
			for (int j = 0; j < dim; ++j) {
				x1n[n * dim + j] = sqrt(1. - beta * beta) *
					smpls[n * dim + j] +
				       	tmpn[n * dim + j] * beta;
			}
			/* Determine if the new proposal is accepted */
			log_alpha[n]  = U(dim, smpls + n * dim) -
					U(dim, x1n + n * dim);
			if (log_alpha[n] > 0.) { log_alpha[n] = 0.; }
			if (log(rndmUniform(seed_r + n)) <= log_alpha[n]) {
				copy(x1n + n * dim, smpls + n * dim, dim);
				++accepted[n];
			}
		}
		accepted[n] /= iter;
		mean_accepted += accepted[n];
	}
	free(log_alpha);
	free(accepted);
	free(tmpn);
	free(allzero);
	free(x1n);
	free(cov);
	return mean_accepted / num_sampl * 100.;
}


/* Pcn-Metropolis algorithm for potentials U */
double uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta)
{
	/* if U is a dim-dimensional R^dim -> R potential function,
	 * we aim at sampling the density exp(-U(x))dx via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */
	setbuf(stdout, NULL);
	double log_alpha = 0;
	double accepted = 0; /* Counters for estimating the acceptance rate */
	double mean_accepted = 0;
	/* x0 and x1 are the previous and nex points in the chain */
	double *x0 = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
	double *cov = malloc(sizeof(double) * dim * dim);
	id(cov, dim);

	double *tmp1 = malloc(sizeof(double) * dim);
	double *tmp2 = malloc(sizeof(double) * dim);
	double *allzero = malloc(sizeof(double) * dim);
	fillzero(tmp1, dim);
	fillzero(tmp2, dim);
	fillzero(allzero, dim);

	/* Produce N samples, each stored in smpls */
	for (int ii = 0; ii < num_sampl; ++ii) {
		if (ii % 10 == 0) {
		printf(".");
		}
		copy(x, x0, dim);	/* Reset the starting value */
		accepted = 0.;		/* Reset the acceptance counter */
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* Weight x0 with the beta square root */
			for (int j = 0; j < dim; ++j) {
				tmp1[j] = sqrt(1. - beta*beta) * x0[j];
			}
			/* Sample the weighted gaussian */
			rndmNdimGaussian(allzero, cov, dim, tmp2, NULL, 0);
			for (int j = 0; j < dim; ++j) {
				tmp2[j] *= beta;
			}
			/* Propose x1 as the sum of the two */
			for (int j = 0; j < dim; ++j) {
				x1[j] = tmp1[j] + tmp2[j];
			}
			/* Determine if the new proposal is accepted */
			log_alpha = U(dim, x0) - U(dim, x1);
			if (log_alpha > 0.) { log_alpha = 0.; }
			if (log(rndmUniform(NULL)) <= log_alpha) {
				copy(x1, x0, dim);
				++accepted;
			}
		}
		/* A sample of the distribution is in x0: copy to smpls */
		copy(x0, smpls + (ii * dim), dim);
		accepted /= iter;
		mean_accepted += accepted;
	}
	free(tmp1);
	free(tmp2);
	free(allzero);
	free(x0);
	free(x1);
	free(cov);
	return mean_accepted / num_sampl * 100.;
}


	/* if U is a dim-dimensional R^dim -> R potential function,
	 * we aim at sampling the density exp(-U(x))dx via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */


double fMpls(double (*f) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, double *smpls)
{
	/* if f is a dim-dimensional R^dim -> R density
	 * we aim at sampling the density f via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */

	double alpha = 0;
	double accepted = 0; /* Counters for estimating the acceptance rate */
	double mean_accepted = 0;
	/* x0 and x1 are the previous and nex points in the chain */
	double *x0 = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
	double *cov = malloc(sizeof(double) * dim * dim);
	id(cov, dim);

	/* Produce N samples, each stored in smpls */
	for (int ii = 0; ii < num_sampl; ++ii) {
		copy(x, x0, dim);	/* Reset the starting value */
		accepted = 0.;		/* Reset the acceptance counter */
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* Gaussian proposal: step from x0 */
			rndmNdimGaussian(x0, cov, dim, x1, NULL, 0);
			alpha = f(dim, x1) / f(dim, x0);
			if (alpha > 1.) { alpha = 1.; }
			/* Determine if accepting the new peoposed step */
			if (rndmUniform(NULL) <= alpha) {
				copy(x1, x0, dim);
				++accepted;
			}
		}
		/* A sample of the distribution is in x0: copy to smpls */
		copy(x0, smpls + (ii * dim), dim);
		accepted /= iter;
		mean_accepted += accepted;
	}
	free(x0);
	free(x1);
	free(cov);
	return mean_accepted / num_sampl;
}
