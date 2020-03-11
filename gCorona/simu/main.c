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

/* Only global variables required: eta, the noise, glob_y.
 * They will be initialized via the function randomInput,
 * then generally used in many function without altering their value */
double *glob_y;
double glob_eta;

void yFromFile (void) {
	FILE* file = fopen("italia.txt", "r");
	int count = 0;
	for (int i = 0; i < glob_dCod; ++i) {
		if (fscanf (file, "%lf", glob_y + i))
		       	++count;
	}
	printf("Read: %d numbers\n", count);
	printf("Read data: ");
	printVec(glob_y, glob_dCod);
//	getchar();
	fclose(file);
}

int positive (const double *x, int dim)
{
	if (x[0] <= 0 || x[1] <= 0 ) //|| x[1] < 10)
		return 0;
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
	G(u, glob_dDom, Gu, glob_dCod);
	/* Return 0.5 * ((y - G(u) / eta)^2 */
	double r = 0.5 * pow(nrm2dist(glob_y, Gu, glob_dCod) / glob_eta, 2.);
	free(Gu);
	return r;
}
/* The probability density from which we want to sample is just
 * exp (-potU(x)) dx, task for which we use three different samples. */

int main() {
	/* Let's start with the variables completely in common between all the
	 * available methods */
	int n_samples = 500;	/* Number of samples to generate */
	int n_iterations = 1000000;
/* Depth of every chain generating 1 sample */
	glob_eta = 0.1;

	/* Initialize the seed: srand() is used e.g. in randomSample above,
	 * while seed_r is used e.g. when parallelization is enables.
	 * Give a different seed to every thread */
	srand(time(NULL));
//	srand(1);
	unsigned int *seed_r = malloc(sizeof(unsigned int) * n_samples);
	assert(seed_r != NULL);
	seed_r[0] = time(NULL) + (unsigned) 1;
//	seed_r[0] = 2;
	for (int i = 1; i < n_samples; ++i) {
		seed_r[i] = seed_r[i - 1] + (unsigned) 1;
	}

	/* Variable for containing the true value of random input x */
	/* Initializing the global y, the observation whose input is
	 * going to be reconstructed */
	glob_y = malloc(sizeof(double) * glob_dCod);
	assert(glob_y != NULL);

	yFromFile();

	/* Container of all the samples and of the one refined thgought
	 * the use of kmeans with "centroid" number of centroids */
	double *raw_samples = malloc(sizeof(double) * glob_dDom * n_samples);
	assert(raw_samples != NULL);
	int centroids = 10;
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
	double avrg_acceptance_rate = 0;
	double map_err = 0;

	/* Now we set variables according the technique utilized: if
	 * pcn, random Hmc or deterministic preconditioned hmc */
	printf("---- PCN MONTE CARLO METHOD ----\n");
	printf("Domain dimension: %d\n", glob_dDom);
	printf("Codomain dimension: %d\n", glob_dCod);
	/* Define beta, parameter for the pcn metropolis monte carlo */
	double beta = 0.1;
	double *start_pt = malloc(sizeof(double) * glob_dDom);
	assert(start_pt != NULL);
	printf("beta: %f\n", beta);
	double *cov = malloc(sizeof(double) * glob_dDom * glob_dDom);
	id(cov, glob_dDom);
	cov[0] = 0.8;
	cov[3] = 1000000.;	/* Specific tuning for the Corona case */
	printf("Covariance matrix:\n");
	printMat(cov, glob_dDom, glob_dDom);

#if PARALLEL
	printf("PARALLELIZED\n");
#endif
//	getchar();
	/* Run a serie of multiple test */
	setbuf(stdout, NULL);
	/* Reset the contenitive variables */
	fillzero(raw_samples, glob_dDom * n_samples);
	fillzero(km_results, (glob_dDom + 1) * centroids);
	/* Starting point suitable for the Corona model */
	start_pt[0] = 0.17;//;fabs(rndmGaussian(0, 1., NULL));
       	start_pt[1] = 3000;//floor(rndmGaussian(100, 20, NULL));
	printf("Starting point: \n");
	printVec(start_pt, 2);
	avrg_acceptance_rate = 
	#if PARALLEL /* Parallel pcn */
        prll_uPcnSampler(potU, glob_dDom, start_pt, n_samples, n_iterations,
                        raw_samples, beta, cov, seed_r, positive);
	#else	/* Ordinary pcn */
        uPcnSampler(potU, glob_dDom, start_pt, n_samples, 
			n_iterations, raw_samples, beta, cov, positive, 0);
	#endif
	printf("\nAverage acceptance rate: %.2f%%\n", avrg_acceptance_rate);

	/* Now the samples have been stored into raw_samples. Clean them kmen*/
	kMeans(raw_samples, n_samples, glob_dDom, centroids, max_kmeans,
			km_results, 0);
	/* Visualize the data */
	kmnsVisual(km_results, centroids, glob_dDom);
	map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, NULL);
	printf("MAP error: %f%%\n", map_err);
	free(seed_r);
	free(glob_y);
	free(cov);
	free(start_pt);
	free(km_results);
	free(raw_samples);
	return 0;
}
