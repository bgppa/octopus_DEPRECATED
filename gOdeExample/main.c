/* General interface for generating automated test by using
 * pCN and / or hamiltonian Monte Carlo methods */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "myblas.h"
#include "mylapack.h"
#include "ranvar.h"
#include "ode.h"
#include "mpls.h"
#include "hmc.h"
#include "kmeans.h"
#include "g.h"

#define PERCENTAGE 3
#define PARALLEL 1
#define PCN 1
#define RAN_HMC 0
#define PRC_HMC 0

/* Only global variables required: eta, the noise, glob_y.
 * They will be initialized via the function randomInput,
 * then generally used in many function without altering their value */
double *glob_y;
double glob_eta;

/* The domain of G does not have any constraint in particular */
int alwaysOk (const double *x, int d)
{
	return 1;
}

/* Thanks to g.h, it is assumed to have access to the globals
 * glob_dDom, glob_dCod, const integer, and the function
 * G(int, const double*, int double) representing the operator
 * G:R^glob_dDom -> R^glob_dCod to invert. They authomatically imply the
 * potential's definition below: */
double potU (int dim, const double *u)
{
	double *Gu = malloc(sizeof(double) * glob_dCod);
	/* Put in Gu the value of G(u) */
	G(u, dim, Gu, glob_dCod);
	/* Return 0.5 * ((y - G(u) / eta)^2 */
	double r = 0.5 * pow(nrm2dist(glob_y, Gu, glob_dCod) / glob_eta, 2.);
	free(Gu);
	return r;
}
/* The probability density from which we want to sample is just
 * exp (-potU(x)) dx, task for which we use three different samples. */

/* Store in x, array of dimension glob_dDom, a random input for G.
 * Then set glob_y = G(u), and add some noise glob_eta depending on
 * the norm of y. */
int randomInput (double *x)
{
	assert(x != NULL);
	/* Generate x as random, gaussian with high covariance */
	double *high_cov = malloc(sizeof(double) * glob_dDom * glob_dDom);
	assert(high_cov != NULL);
	double *zeros = malloc(sizeof(double) * glob_dDom);
	assert(zeros != NULL);
	fillzero(zeros, glob_dDom);
	/* Set the covariance as the identity, then multiply its diagonal */
	id(high_cov, glob_dDom);
	for (int i = 0; i < glob_dDom; ++i) {
		high_cov[i * glob_dDom + i] *= 200.;
	}
	rndmNdimGaussian(zeros, high_cov, glob_dDom, x, NULL, 0);
	printf("X:\t");
	printVec(x, glob_dDom);
	/* Now produce y = G(x) + noise */
	G(x, glob_dDom, glob_y, glob_dCod);
	printf("G(X):\t");
	printVec(glob_y, glob_dCod);
	//glob_eta = nrm2(glob_y, glob_dCod) / 100. * (double) PERCENTAGE;
	glob_eta = 0.1;
	printf("ETA:\t%f\n", glob_eta);
	for (int i = 0; i < glob_dCod; ++i) {
		glob_y[i] += rndmGaussian(0., glob_eta, NULL);
	}
	printf("G(X) + NOISE: ");
	printVec(glob_y, glob_dCod);
	free(zeros);
	free(high_cov);
	return 1;
}

int main() {
	/* Let's start with the variables completely in common between all the
	 * available methods */
	int n_samples = 2000;	/* Number of samples to generate */
	int n_iterations = 1000;/* Depth of every chain generating 1 sample */

	/* Initialize the seed: srand() is used e.g. in randomSample above,
	 * while seed_r is used e.g. when parallelization is enables.
	 * Give a different seed to every thread */
	srand(time(NULL));
	unsigned int *seed_r = malloc(sizeof(unsigned int) * n_samples);
	assert(seed_r != NULL);
	seed_r[0] = time(NULL) + (unsigned) 1;
	for (int i = 1; i < n_samples; ++i) {
		seed_r[i] = seed_r[i - 1] + (unsigned) 1;
	}

	/* Variable for containing the true value of random input x */
	double *true_x = malloc(sizeof(double) * glob_dDom);
	assert(true_x != NULL);
	/* Initializing the global y, the observation whose input is
	 * going to be reconstructed */
	glob_y = malloc(sizeof(double) * glob_dCod);
	assert(glob_y != NULL);

	/* Container of all the samples and of the one refined thgought
	 * the use of kmeans with "centroid" number of centroids */
	double *raw_samples = malloc(sizeof(double) * glob_dDom * n_samples);
	assert(raw_samples != NULL);
	int centroids = 20;
	int max_kmeans = 500;	/* Maximum number of iteration in the kmeans */
	double *km_results = malloc(sizeof(double) * (glob_dDom+1) * centroids);
	assert(km_results != NULL);

	/* This script executes multiple tests in an automated way:
	(1) generate a random input;
 	(2) set y = G(u) + noise
	(3) recontruct x from y, estimating error and residuals.
	A test is considered "successful" if the highest probably error
	is less than tol_err%. It is suggested to execute multiple times a test
	with the same initial condition on x - as written in the for below -
	since sometimes a lower error is obtained with the 2nd or 3rd
	candidate x rather than with the first (probabilistic reasons */
	int tot_test = 10;
	double tol_err = 10;	
	int success = 0;
	double avrg_acceptance_rate = 0;
	double map_err = 0;

	/* Now we set variables according the technique utilized: if
	 * pcn, random Hmc or deterministic preconditioned hmc */
#if PCN
	printf("---- PCN MONTE CARLO METHOD ----\n");
	/* Define beta, parameter for the pcn metropolis monte carlo */
	double beta = 0.2;
	double *start_pt = malloc(sizeof(double) * glob_dDom);
		assert(start_pt != NULL);
	printf("beta: %f\n", beta);
	#else
	printf("---- HAMILTONIAN MONTE CARLO METHOD ----\n");
	/* Hamiltonian Monte Carlo activated */
	double *start_pt = malloc(sizeof(double) * glob_dDom * 2);
	assert(start_pt != NULL);
	/* Mass matrix for the preconditioned hamiltonian mc.
	 * As a default choice, the identity */
	double *M = malloc(sizeof(double) * glob_dDom * glob_dDom);
	assert(M != NULL);
	id(M, glob_dDom);
	/* Inverse of M, necessary for the Verlet integrator.
	 * Better to compute it now one times */
	double *M1 = malloc(sizeof(double) * glob_dDom * glob_dDom);
	assert(M1 != NULL);
	/* Since M has been left with the identity matrix, M1 is id too */
	id(M1, glob_dDom);
	printf("Mass matrix:\n");
	printMat(M, glob_dDom, glob_dDom);

	#if PRC_HMC
	/* In case of deterministic preconditioned hmc */
	/* Each step of the chain is a Verlet integration from time 0
	 * to time_interval, with n_steps in between */
	double time_interval = 3.0;
	int n_steps = 200;
	printf("PRECONDITIONED, DETERMINISTIC\n");
	printf("time_interval: %f, n_steps: %d\n", time_interval, n_steps);
	#else	
	printf("PRECONDITIONED, RANDOMIZED\n");
	/* Othewise the quantities above are randomized: each chain step
	 * is Verlet with fixed step h, but until time geom(h/lam) */
	double h = 0.01;
	double lam = 1.0;
	printf("h: %f, lambda: %f\n", h, lam);
	#endif
#endif

#if PARALLEL
	printf("PARALLELIZED\n");
#endif
//	getchar();
	/* Run a serie of multiple test */
	for (int a = 0; a < tot_test; ++a) {
		/* Reset the contenitive variables */
		fillzero(true_x, glob_dDom);
		fillzero(glob_y, glob_dDom);
		fillzero(raw_samples, glob_dDom * n_samples);
		fillzero(km_results, (glob_dDom + 1) * centroids);
#if PCN
		fillzero(start_pt, glob_dDom);
#else
		/* The Hamiltonian algorithm double the dimension */
		fillzero(start_pt, glob_dDom * 2);
#endif
		printf("###### TEST %d of %d ######\n", a + 1, tot_test);
		randomInput(true_x);
		avrg_acceptance_rate = 
#if PCN
	#if PARALLEL /* Parallel pcn */
        prll_uPcnSampler(potU, glob_dDom, start_pt,
                        n_samples, n_iterations, raw_samples, beta, seed_r,
			alwaysOk);
	#else	/* Ordinary pcn */
        uPcnSampler(potU, glob_dDom, start_pt,
                        n_samples, n_iterations, raw_samples, beta,
			alwaysOk);
	#endif
#elif RAN_HMC /* Randomized version */
        #if PARALLEL	/* Parallel randomized hmc */
	prll_pRanHmcSampler(glob_dDom * 2, start_pt, h, lam, 
                M, M1, potU, n_iterations, n_samples, raw_samples, seed_r,
		alwaysOk);
        #else	/* Ordinary randomized hmc */ 
        pRanHmcSampler(glob_dDom * 2, start_pt, h, lam,
                M1, M1, potU, n_iterations, n_samples, raw_samples,
		alwaysOk);
        #endif
#else /* Deterministic preconditioned hmc */
        #if PARALLEL	/* Parallel preconditioned */
        prll_pHmcSampler(glob_dDom * 2, start_pt,
                time_interval, n_steps, M1, M1,
                potU, n_iterations, n_samples, raw_samples, seed_r,
		alwaysOk);
        #else	/* Ordinary preconditioned hmc */
        pHmcSampler(glob_dDom * 2, start_pt,
                time_interval, n_steps, M1, M1,
                potU, n_iterations, n_samples, raw_samples, alwaysOk);
        #endif
#endif
	printf("\nAverage acceptance rate: %.2f%%\n", avrg_acceptance_rate);

	/* Now the samples have been stored into raw_samples. Clean them kmen*/
	kMeans(raw_samples, n_samples, glob_dDom, centroids, max_kmeans,
			km_results, 0);
	/* Visualize the data */
	kmnsVisual(km_results, centroids, glob_dDom);
	map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, true_x);
	printf("MAP error: %f%%\n", map_err);
	if (map_err < tol_err) { ++success; }
	}
	printf("%d OK out of %d\n", success, tot_test);

	free(seed_r);
	free(true_x);
	free(glob_y);
	free(start_pt);
	free(km_results);
	free(raw_samples);
	return 0;
}
