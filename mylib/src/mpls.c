/* Metropolis, a module implementig various naive Metropolis Markov
 * Chain Monte Carlo approaches */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "omp.h"
#include "myblas.h"
#include "ranvar.h"
#include "mpls.h"

/* NOTE: for simplicity, the covariance matrices for the Metropolis
 * sampling are assumed to be diagonal, in order to spare time.
 * In the future you should extend this functionalisy, easy easy,
 * just replace rndmDiagGauss with rndmNdimGauss - for the moment I prefer
 * continuing using the first because faster and OK for my cases */

/* Enable the parallel support for the Metropolis module.
 * When ENABLED, the only problem will be a false positive error
 * when using valgrind. Therefore disabling it can be useful
 * for temporary debugging processes */
#define MPLS_PARALL 1

/* Pcn-Metropolis algorithm for potentials U, NON parallel */
/* if U is a dim-dimensional R^dim -> R potential function,
 * we aim at sampling the density exp(-U(x))dx via a simple random
 * walk metropolis algorithm. 
 * produce num_sampl samples of is. Storen in "samples" assumed to be
 * an array of dimension num_sampl * dim. x0 repr. the starting point */

/* IN THE NEW VERSION BELOW, x is an ARRAY OF DIMENSION dim * num_sampl,
 * IN A WAY TO CONTAIN A POSSIBLY DIFFERENT STARTING POINT FOR EVERY CHAIN */
 
double old_uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta, const double *cov,
	       	int (*okconstraint) (const double*, int), int verbose)
{
	setbuf(stdout, NULL);
	assert(okconstraint != NULL);
	double log_alpha = 0;
	/* Counters for estimating the acceptance rate */
	double accepted = 0; 
	double mean_accepted = 0;
	/* x0 and x1 are the previous and next points in the chain */
//	double *x0 = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	double *tmp = malloc(sizeof(double) * dim);
	/* Produce N samples, each stored in smpls */
	for (int n = 0; n < num_sampl; ++n) {
		if (n % 10 == 0) {
			printf(".");
		}
		copy (x, smpls + n * dim, dim);
//		copy(x, x0, dim);	/* Reset the starting value */
		accepted = 0.;		/* Reset the acceptance counter */
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* tmp contains centered gaussian covariance cov */
			fillzero(tmp, dim);
			rndmDiagGauss(tmp, cov, dim, NULL);
			/* Propose x1 as the sum of the two */
			for (int j = 0; j < dim; ++j) {
//				x1[j] = sqrt(1. - beta*beta) * x0[j]
				x1[j] = sqrt(1.-beta*beta) * (smpls+n*dim)[j]
					 + beta * tmp[j];
			}
			if (verbose) {
				printf("[step %d] Proposed: ", i);
				printVec(x1, dim);
			}
			if (okconstraint(x1, dim)) {
				/* Determine if the new proposal is accepted */
//				log_alpha = min(U(dim, x0) - U(dim, x1), 0.);
	log_alpha = min(U(dim, smpls+n*dim) - U(dim, x1), 0.);
				if (log(rndmUniform(NULL)) <= log_alpha) {
//					copy(x1, x0, dim);
				copy(x1, smpls+n*dim, dim);
					++accepted;
					if (verbose) {
						printf("\tAccepted!\n");
					}
				}
			}	/* Has been accepted or not */
		}	/* Single sample produced, is in x0 */
//		copy(x0, smpls + (n * dim), dim);
		accepted /= iter;
		mean_accepted += accepted;
	}	/* num_sampl samples produced */
	free(tmp);
//	free(x0);
	free(x1);
	return mean_accepted / num_sampl * 100.;
}

/* Pcn-Metropolis algorithm for potentials U, NON parallel */
/* if U is a dim-dimensional R^dim -> R potential function,
 * we aim at sampling the density exp(-U(x))dx via a simple random
 * walk metropolis algorithm. 
 * produce num_sampl samples of is. Storen in "samples" assumed to be
 * an array of dimension num_sampl * dim. x0 repr. the starting point */
double uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta, const double *cov,
	       	int (*okconstraint) (const double*, int), int verbose)
{
	setbuf(stdout, NULL);
	assert(okconstraint != NULL);
	double log_alpha = 0;
	/* Counters for estimating the acceptance rate */
	double accepted = 0; 
	double mean_accepted = 0;
	/* x0 and x1 are the previous and next points in the chain */
//	double *x0 = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	double *tmp = malloc(sizeof(double) * dim);
	/* Produce N samples, each stored in smpls */
	for (int n = 0; n < num_sampl; ++n) {
		if (n % 10 == 0) {
			printf(".");
		}
		copy (x, smpls, dim * num_sampl);
//		copy(x, x0, dim);	/* Reset the starting value */
		accepted = 0.;		/* Reset the acceptance counter */
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* tmp contains centered gaussian covariance cov */
			fillzero(tmp, dim);
			rndmDiagGauss(tmp, cov, dim, NULL);
			/* Propose x1 as the sum of the two */
			for (int j = 0; j < dim; ++j) {
//				x1[j] = sqrt(1. - beta*beta) * x0[j]
				x1[j] = sqrt(1.-beta*beta) * (smpls+n*dim)[j]
					 + beta * tmp[j];
			}
			if (verbose) {
				printf("[step %d]", i);
				printf(" current: ");
				printVec(smpls + n*dim, dim);
				printf(".......proposed: ");
				printVec(x1, dim);
			}
			if (okconstraint(x1, dim)) {
				/* Determine if the new proposal is accepted */
//				log_alpha = min(U(dim, x0) - U(dim, x1), 0.);
	log_alpha = min(U(dim, smpls+n*dim) - U(dim, x1), 0.);
	printf("\t\talpha = %f ", exp(log_alpha));
	getchar();
				if (log(rndmUniform(NULL)) <= log_alpha) {
//					copy(x1, x0, dim);
				copy(x1, smpls+n*dim, dim);
					++accepted;
					if (verbose) {
						printf("\tAccepted!\n");
					}
				}
			}	/* Has been accepted or not */
		}	/* Single sample produced, is in x0 */
//		copy(x0, smpls + (n * dim), dim);
		accepted /= iter;
		mean_accepted += accepted;
	}	/* num_sampl samples produced */
	free(tmp);
//	free(x0);
	free(x1);
	return mean_accepted / num_sampl * 100.;
}

/* if U is a dim-dimensional R^dim -> R potential function,
 * we aim at sampling the density exp(-U(x))dx via a simple random
 * walk metropolis algorithm. 
 * produce num_sampl samples of is. Storen in "samples" assumed to be
 * an array of dimension num_sampl * dim. x0 repr. the starting point */
/* Pcn-Metropolis algorithm for potentials U */
/* This NEW samplers requires the user to give start_pt as a list of 
 * starting points, x, so an array of dimension dim * num_sampl */
double prll_uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta, const double *cov,
	       	unsigned int *seed_r,
		int (*okconstraint) (const double *, int))
{
	setbuf(stdout, NULL);
	assert(okconstraint != NULL);
#if! MPLS_PARALL
	printf("\n\n*warning: prll_uPcnSampler: calling a parallel function,"
		" but the support is DISABLED. Check the"
	        " constant MPLS_PARALL in mpls.c\n");
	printf("*It should work anyway, but with no benefits*\n");
#endif
	double *log_alpha = malloc(sizeof(double) * num_sampl);
	/* Counters for estimating the acceptance rate */
	double *accepted = malloc(sizeof(double) * num_sampl);
	fillzero(accepted, num_sampl);
	double mean_accepted = 0;
	/* x1n contains the various proposals */
	double *x1n = malloc(sizeof(double) * dim * num_sampl);
	fillzero(x1n, dim * num_sampl);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
	double *allzero = malloc(sizeof(double) * dim);
	fillzero(allzero, dim);
	double *tmpn = malloc(sizeof(double) * dim * num_sampl);
	fillzero(tmpn, dim * num_sampl);
	/* Fill smpls with multiple copies of x */
	copy (x, smpls, dim * num_sampl);
	/*
	for (int i = 0; i < num_sampl; ++i) {
		copy(x + i * dim, smpls + i * dim, dim);
	}
	*/
	/* Produce N samples, each stored in smpls */
#if MPLS_PARALL
	#pragma omp parallel for
#endif
	for (int n = 0; n < num_sampl; ++n) {
		if (n % 100 == 0) {
			printf(".");
		}
		for (int i = 0; i < iter; ++i) { /* Produce a sample */
			/* Sample the gaussian in tmp */
			fillzero (tmpn + n * dim, dim);
			rndmDiagGauss (tmpn + n * dim, cov, dim, seed_r + n);
//			rndmNdimGaussian(allzero, cov, dim,
//				       	tmpn + n * dim, seed_r + n, 0);
			/* Propose x1 as the weightes sum between that gaussian
			 * and the previous x0, which is smpls[n*dim] */
			for (int j = 0; j < dim; ++j) {
				x1n[n * dim + j] = sqrt(1. - beta * beta) *
					smpls[n * dim + j] +
				       	tmpn[n * dim + j] * beta;
			}
			if (okconstraint(x1n + n * dim, dim)) {
			/* Determine if the new proposal is accepted */
				log_alpha[n]  =  min(
					U(dim, smpls + n * dim) -
					U(dim, x1n + n * dim), 0.);
				if (log(rndmUniform(seed_r + n))
						<= log_alpha[n]) {
				copy(x1n + n * dim, smpls + n * dim, dim);
				++accepted[n];
//				printf("n : %d. Current accepted: %f\n",n,
//						accepted[n]);
//				getchar();
				}
			}
		}
	}
	for (int n = 0; n < num_sampl; ++n) {
		mean_accepted += accepted[n] / iter;
	}
	free(log_alpha);
	free(accepted);
	free(tmpn);
	free(allzero);
	free(x1n);
	return mean_accepted / num_sampl * 100.;
}

/* if U is a dim-dimensional R^dim -> R potential function,
 * we aim at sampling the density exp(-U(x))dx via a simple random
 * walk metropolis algorithm. 
 * produce num_sampl samples of is. Storen in "samples" assumed to be
 * an array of dimension num_sampl * dim. x0 repr. the starting point */
/* Pcn-Metropolis algorithm for potentials U */
double old_prll_uPcnSampler (double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, 
		double *smpls, double beta, const double *cov,
	       	unsigned int *seed_r,
		int (*okconstraint) (const double *, int))
{
	setbuf(stdout, NULL);
	assert(okconstraint != NULL);
#if! MPLS_PARALL
	printf("\n\n*warning: prll_uPcnSampler: calling a parallel function,"
		" but the support is DISABLED. Check the"
	        " constant MPLS_PARALL in mpls.c\n");
	printf("*It should work anyway, but with no benefits*\n");
#endif
	double *log_alpha = malloc(sizeof(double) * num_sampl);
	/* Counters for estimating the acceptance rate */
	double *accepted = malloc(sizeof(double) * num_sampl);
	fillzero(accepted, num_sampl);
	double mean_accepted = 0;
	/* x1n contains the various proposals */
	double *x1n = malloc(sizeof(double) * dim * num_sampl);
	fillzero(x1n, dim * num_sampl);
	/* We'll use a simple gaussian rw with covariance cov, the identity */
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
			fillzero (tmpn + n * dim, dim);
			rndmDiagGauss (tmpn + n * dim, cov, dim, seed_r + n);
//			rndmNdimGaussian(allzero, cov, dim,
//				       	tmpn + n * dim, seed_r + n, 0);
			/* Propose x1 as the weightes sum between that gaussian
			 * and the previous x0, which is smpls[n*dim] */
			for (int j = 0; j < dim; ++j) {
				x1n[n * dim + j] = sqrt(1. - beta * beta) *
					smpls[n * dim + j] +
				       	tmpn[n * dim + j] * beta;
			}
			if (okconstraint(x1n + n * dim, dim)) {
			/* Determine if the new proposal is accepted */
				log_alpha[n]  =  min(
					U(dim, smpls + n * dim) -
					U(dim, x1n + n * dim), 0.);
				if (log(rndmUniform(seed_r + n))
						<= log_alpha[n]) {
				copy(x1n + n * dim, smpls + n * dim, dim);
				++accepted[n];
				}
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
	return mean_accepted / num_sampl * 100.;
}

/* TO TEST, COMPLETELY EXPERIMENTAL AND POSSIBLY WRONG */
double simple_pcn (double (*U) (int, const double*), 
		const double *start_pt, int dim,
		double *chain, int len,
		double beta, const double *cov,
		double burning_percentage,
		int (*okconstraint) (const double *, int))
{
	setbuf(stdout, NULL);
	assert(okconstraint != NULL);

	/* Set the burning time */
	int bt = len * burning_percentage;
	printf("Burning time: %d\n", bt);
	int iter = len + bt;

	double accepted = 0;
	double log_alpha = 0;
	double *x0 = malloc(sizeof(double) * dim);
	double *tmp = malloc(sizeof(double) * dim);
	double *x1 = malloc(sizeof(double) * dim);
	int chain_counter = 0;

	copy(start_pt, x0, dim);
	for (int k = 0; k < iter; ++k) {
		if (k % 100 == 0) {
			printf(".");
		}
		rndmDiagGauss (tmp, cov, dim, NULL);
		/* Propose x1 as the weightes sum between that gaussian
		 * and the previous x0, which is smpls[n*dim] */
		for (int j = 0; j < dim; ++j) {
			x1[j] = sqrt(1. - beta * beta) * x0[j] + tmp[j] * beta;
		}
		if (okconstraint(x1, dim)) {
//			printf("Current point:\n");
//			printVec(x0, dim);
			/* Determine if the new proposal is accepted */
//			printf("Proposed: \n");
//			printVec(x1, dim);
//			getchar();
			log_alpha = min(U(dim, x0) - U(dim, x1), 0.);
//			printf("log_alpha: %e\n", log_alpha);
			if (log(rndmUniform(NULL)) <= log_alpha) {
					copy(x1, x0, dim);
					++accepted;
//					printf("Accepted\n");
//				printf("n : %d. Current accepted: %f\n",n,
//						accepted[n]);
//				getchar();
			}
		}
		else{--k;};
		/* If we surpassed the burning time, copy the current state
		 * into the chain */
		if (k > bt || k == bt) {
			copy(x0, chain + chain_counter * dim, dim);
			++chain_counter;
		}
	}
	free(x0);
	free(x1);
	free(tmp);
//	printf("Produced chain:\n");
//	getchar();
//	printMat(chain, len, dim);
//	getchar();
	return accepted * 100. / iter;
}


#if 0
/* Single metropolis sampling, maybe not useful. Comment to spare space */
double uMpls(double (*U) (int, const double*), int dim, 
		const double *x, int num_sampl, int iter, double *smpls,
		int (*okconstraint) (const double *, int))
{
	/* if U is a dim-dimensional R^dim -> R potential function,
	 * we aim at sampling the density exp(-U(x))dx via a simple random
	 * walk metropolis algorithm. 
	 * produce num_sampl samples of is. Storen in "samples" assumed to be
	 * an array of dimension num_sampl * dim. x0 repr. the starting point */

	assert(okconstraint != NULL);
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
			/* If the proposed point satisfies some
			 * additional constraints */
			if (okconstraint != NULL && okconstraint(x1, dim)) {
				alpha = exp(U(dim, x0) - U(dim, x1));
				if (alpha > 1.) { alpha = 1.; }
				/* Determine accepting the new peoposed step */
				if (rndmUniform(NULL) <= alpha) {
					copy(x1, x0, dim);
					++accepted;
				}
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

/* Thi function need to be checked again. It samples given a density f,
 * while the functions above requires U and assume density of the form
 * exp(-U))*/
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
#endif
