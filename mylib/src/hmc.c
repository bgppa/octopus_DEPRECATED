/* Hamiltonian Monte Carlo Methods */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "myblas.h"
#include "ranvar.h"
#include "ode.h"
#include "hmc.h"

/* GENERAL REMARK: the following functions uses Verlet Integrator
 * in order to sample by using Hamiltonian Monte Carlo methods.
 * They ASSUME and RELY on the definition of the global variables
 * ham_U, ham_M1, as described in the header ode.h where the Verlet
 * integrator has been defined */

/* --------------------- FIRST PART: THE SINGLE STEP */
/* This is a single step of a hmc chain. Given a starting point x,
 * produce the following by integrating the Hamiltonian ODE until
 * time "time", N steps, by using Verlet. Verbose when v = 1.
 * The particle moves in dimension d, so Hamiltonian defined in R^2d */
int hmc_single_step(int d2, double *x, double time, int N, const double *M,
		const double *M1, double (*U) (int, const double*), int v,
		unsigned int *my_seed, 
		int (*okconstraint) (const double *, int))
{
	/* ADD INFO: but M1 and M are inverse of each other */
	/* M is assumed to have been given as the inverse of M1 */
	int accpt = 0;
	double *zeroes = malloc(sizeof(double) * d2 / 2);
	double *xi0 = malloc(sizeof(double) * d2 / 2);
	fillzero(zeroes, d2/2);
	rndmNdimGaussian(zeroes, M, d2 / 2, xi0, my_seed, 0);
	/* Set the starting point as the half x, half xi0 */
	double *x_prev = malloc(sizeof(double) * d2);
	for (int i = 0; i < d2; ++i) {
		if (i < d2 / 2) { x_prev[i] = x[i]; }
		else		{ x_prev[i] = xi0[i - d2/2]; }
	}
	/* Solve Verlet starting from this point */
	double delta_H = verlet(x_prev, d2, time, N, M1, U, 0);
	if (v) {printf("delta_H = %e\n", delta_H);}
	/* If the proposed point satisfied some given domain constraint,
	 * continue with the metropolis acceptance */
	if (okconstraint(x_prev, d2)) {
		double alpha = exp(-delta_H) < 1. ? exp(-delta_H) : 1.;
		/* If energy is sufficiently preserved, accept and update */
		if (rndmUniform(my_seed) <= alpha)
			{ accpt = 1, copy(x_prev, x, d2); }
	}
	free(x_prev), free(zeroes), free(xi0);
	return accpt;
}

/* ----------------------------- SECOND PART: CHAINS
 *
 * We now construct hmc chains based on the single step above.
 * Two possibilities, according to the time until the Verlet integrator
 * is run. The simplest is deterministic, pHmcChain, the more advanced
 * includes a RANDOMIZED time according to Nawaf's paper */

/* Preconditioned RANDOMIZED Hamiltonian Metropolis CHAIN */
/* Params: d2 is the dimension, x the starting point, h the fixed integration
 * step, lam is lambda, the "intensity": each step integrates the system
 * until time Geometric(h/lambda), ad described in the paper.
 * M is the precondition matrix, while chain_length says by itself.
 * Return the %acceptance rate of the samples */
double pRanHmcChain(int d2, double *x, double h, double lam, const double* M, 
		const double *M1, double (*U) (int, const double*),
		int chain_length, unsigned int *my_seed,
		int (*okconstraint) (const double *, int)) {
	/* Return a sample after chain_length hcm_iterations */
	double n_in_interval = 0;
	double time_interval = 0;
	double accpt = 0;
	for (int i = 0; i < chain_length; ++i) {
		/* By theoretical definition, the number of steps is a 
		 * geometrical rv with mean lam / h (so p = h/lam). */
		n_in_interval = rndmGeometricInt(h/lam, my_seed);
		time_interval = n_in_interval * h;
		accpt += 
			hmc_single_step(d2, x, time_interval, n_in_interval,
					M, M1, U, 1, my_seed, okconstraint);
getchar();
	}
	return accpt * 100. / (double) chain_length;
}

/* Preconditioned Hamiltonian Metropolis chain, deterministic time interval
 * and number of steps. Return the acceptance rate in percentage % */
double pHmcChain(int d2, double *x, double time_interval, int steps_interval, 
		const double* M, 
		const double* M1, 
		double (*U) (int, const double*),
		int chain_length, unsigned int *specific_seed,
		int (*okconstraint) (const double *, int)) {
	/* Return a sample after chain_length hcm_iterations */
	double accpt = 0;
	for (int i = 0; i < chain_length; ++i) {
		accpt += hmc_single_step(d2, x, time_interval, 
				steps_interval, M, M1, U, 0, specific_seed,
				okconstraint);
	}
	return accpt * 100. / (double) chain_length;
}


/* -------------------- THIRD PART: COMPLETE SAMPLERS
 * Each chain defined above produce ONE sample, which is just the bulding block
 * for statistical analysis where we need many of them.
 * That's the purpose of the two functions below, using the randomized or
 * the deterministic hmc respectively */

/* NOTE: if the potential U: R^d -> R, hmc produces samples on the doubled space
 * R^2d, since we add the momentum. If we take the *second* half of these
 * samples, we obtain the right samples from U */
/* Parameters:
 * d2 is the dimension 2 * d of the system;
 * x is the starting point
 * h the fixed Verlet integration step
 * lam the intensity (see the randomized chain above)
 * M the precondition matrix
 * chain_length
 * n_samples is the number of samples
 * raw_samples an array of dimension n_samples * d2/2 containing the results.
 * Remember that only the first half of the sample is needed, dimension d2/2,
 * being so the marginal distribution relative to U.
 * Return: the average acceptance rate */
double pRanHmcSampler(int d2, const double *x, double h, double lam,
	       	const double *M,  
		const double* M1, 
		double (*U) (int, const double*),
		int chain_length, 
		int n_samples, double *raw_samples,
		int (*okconstraint) (const double*, int))
{
	setbuf(stdout, NULL);
	/* Copy the starting condition */
	double *start_copy = malloc(sizeof(double) * d2);
	copy(x, start_copy, d2);
	double acceptance = 0;
	/* raw_samples contain all the position samples */
	fillzero(raw_samples, d2 / 2 * n_samples);
	printf("Sampling (%d in total)\n", n_samples);
	for (int i = 0; i < n_samples; ++i) {
		if (i % 10 == 0 ) {
			printf(".");
		}
		/* For each sample extraction, reset the initial condition */
		copy(x, start_copy, d2);
		/* Perform a single sample extraction from randomized hmc */
		acceptance += pRanHmcChain(d2, start_copy, h, lam,
			       M, M1, U, chain_length, NULL, okconstraint);
		/* The sample is contained in start_copy: copy its FIRST
		 * half into raw_samples: this part represents sample from U */
		copy(start_copy, raw_samples + i * (d2 / 2), d2/2);
	} printf("\n");
	free(start_copy);
	return acceptance / n_samples;
}

/* Parallel version of the algorithm above */
double prll_pRanHmcSampler(int d2, const double *x, double h, double lam,
	       	const double *M,  
		const double* M1, 
		double (*U) (int, const double*),
		int chain_length, 
		int n_samples, double *raw_samples, unsigned int *prll_seed,
		int (*okconstraint) (const double *, int))
{
	setbuf(stdout, NULL);
	double *accpt_i = malloc(sizeof(double) * n_samples);
	fillzero(accpt_i, n_samples);
	double acceptance = 0;
	/* raw_samples contain all the position samples */
	fillzero(raw_samples, (d2 / 2) * n_samples);
	printf("Sampling (%d in total)\n", n_samples);
	/* Array containing multiple copies of the starting condition,
	 * one for each n_sample. Better to have multiple copies in order
	 * not to have possible unexpected problems during parallelization */
	double *n_start_copy = malloc(sizeof(double) * d2 * n_samples);
	for (int i = 0; i < n_samples; ++i) {
		copy(x, n_start_copy + i * d2, d2);
	}
	#pragma omp parallel for 
	for (int i = 0; i < n_samples; ++i) {
		if (i % 10 == 0) {
			printf(".");
		}
		/* Remark: for each index i, n_start_copy + i*d2 is the
		 * position of the i-th copy of the starting point.
		 * This will we overwritten with the chain's result,
		 * and need finally to be partially copied (only the first
		 * half, being the margina w.r.t. our interested U */
		/* For each sample extraction, reset the initial condition */
		/* Perform a single sample extraction from randomized hmc */
		accpt_i[i] = pRanHmcChain(d2, n_start_copy + i * d2, h, lam,
			       M, M1, U, chain_length, prll_seed + i,
			       okconstraint);
		copy(n_start_copy + i * d2, raw_samples + i * (d2 / 2), d2/2);
	}
	for (int i = 0; i < n_samples; ++i) {
		acceptance += accpt_i[i];
	}
	free(accpt_i);
	free(n_start_copy);
	return acceptance / n_samples;
}



/* NOTE: if the potential U: R^d -> R, hmc produces samples on the doubled space
 * R^2d, since we add the momentum. If we take the *second* half of these
 * samples, we obtain the right samples from U */
/* Parameters:
 * d2 is the dimension 2 * d of the system;
 * x is the starting point
 * time_interval the fixed time interval for Verlet integration step
 * n_single_step the number of divisions of the time interval
 * M the precondition matrix
 * chain_length
 * n_samples is the number of samples
 * raw_samples an array of dimension n_samples * d2/2 containing the results.
 * Remember that only the first half of the sample is needed, dimension d2/2,
 * being so the marginal distribution relative to U.
 * Return: the average acceptance rate */
double pHmcSampler(int d2, const double *x, double time_interval, 
		int n_single_step, const double *M,  
		const double* M1, 
		double (*U) (int, const double*),
		int chain_length, 
		int n_samples, double *raw_samples,
		int (*okconstraint) (const double*, int))
{
	setbuf(stdout, NULL);
	/* Copy the starting condition */
	double *start_copy = malloc(sizeof(double) * d2);
	copy(x, start_copy, d2);
	double acceptance = 0;
	/* raw_samples contain all the position samples */
	fillzero(raw_samples, d2 / 2 * n_samples);
	printf("Sampling (%d in total)\n", n_samples);
	for (int i = 0; i < n_samples; ++i) {
		if (i % 10 == 0) {
			printf(".");
		}
		/* For each sample extraction, reset the initial condition */
		copy(x, start_copy, d2);
		/* Perform a single sample extraction */
		acceptance += pHmcChain(d2, start_copy, time_interval,
			       n_single_step, M, M1, U, chain_length,
			       NULL, okconstraint);
		/* The sample is contained in start_copy: copy its FIRST
		 * half into raw_samples: this part represents sample from U */
		copy(start_copy, raw_samples + i*(d2/2), d2/2);
	} printf("\n");
	free(start_copy);
	return acceptance / n_samples;
}

/* Parallel version of the algorithm above */
double prll_pHmcSampler(int d2, const double *x, double time_interval, 
		int n_single_step, const double *M,  
		const double* M1, 
		double (*U) (int, const double*),
		int chain_length, 
		int n_samples, double *raw_samples, unsigned int *prll_seed,
		int (*okconstraint) (const double*, int))
{
	setbuf(stdout, NULL);
	double *accpt_i = malloc(sizeof(double) * n_samples);
	fillzero(accpt_i, n_samples);
	double acceptance = 0;
	/* raw_samples contain all the position samples */
	fillzero(raw_samples, (d2 / 2) * n_samples);
	printf("Sampling (%d in total)\n", n_samples);
	/* Array containing multiple copies of the starting condition,
	 * one for each n_sample. Better to have multiple copies in order
	 * not to have possible unexpected problems during parallelization */
	double *n_start_copy = malloc(sizeof(double) * d2 * n_samples);
	for (int i = 0; i < n_samples; ++i) {
		copy(x, n_start_copy + i * d2, d2);
	}
	#pragma omp parallel for 
	for (int i = 0; i < n_samples; ++i) {
		if (i % 10 == 0) {
			printf(".");
		}
		/* Remark: for each index i, n_start_copy + i*d2 is the
		 * position of the i-th copy of the starting point.
		 * This will we overwritten with the chain's result,
		 * and need finally to be partially copied (only the first
		 * half, being the margina w.r.t. our interested U */
		/* For each sample extraction, reset the initial condition */
		/* Perform a single sample extraction from randomized hmc */
		accpt_i[i] = pHmcChain(d2, n_start_copy + i * d2, time_interval,
				n_single_step,
			       M, M1, U, chain_length, prll_seed + i,
			       okconstraint);
		copy(n_start_copy + i * d2, raw_samples + i * (d2 / 2), d2/2);
	}
	for (int i = 0; i < n_samples; ++i) {
		acceptance += accpt_i[i];
	}
	free(accpt_i);
	free(n_start_copy);
	return acceptance / n_samples;
}
