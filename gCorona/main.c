/* General interface for generating automated test by using
 * pCN and / or hamiltonian Monte Carlo methods */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "myblas.h"
#include "mylapack.h"
#include "ranvar.h"
#include "mpls.h"
#include "kmeans.h"
#include "g.c"

#define PERCENTAGE 3 /* Check better */
#define PARALLEL 1
#define TEST_MODE 1

/* Only global variables required: eta, the noise, glob_y.
 * They will be initialized via the function randomInput,
 * then generally used in many function without altering their value */
double *glob_y;
double *glob_eta;

/* Control function: ensures that the parameters for the Gompertz
 * model are always positive - to be given to the pcn samples */
int positive (const double *x, int dim)
{
	if (x[0] > 0 && x[0] < 1 && x[1] > 0)
		return 1;
	return 0;
/*
	for (int i = 0; i < glob_dDom; ++i) 
		if (x[i] < 0)
			return 0;
	return 1;
	*/
}
/* For the exponential interpolation, no restriction on parameters */
int allok(const double *x, int dim) {return 1;}

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
	/* New version, includes the possibility of a multidimensional noise.
	 * Gu will (y - G(u) / eta)_i, vectorized form: */
	for (int i = 0; i < glob_dCod; ++i) {
		Gu[i] = (glob_y[i] - Gu[i]) / sqrt(glob_eta[i]);
	}
	double r = 0.5 * pow(nrm2(Gu, glob_dCod), 2.);
	free(Gu);
	return r;
}


/* yFromFile and randomInput initialize glob_y and glob_eta, the observation
 * from where to reconstruct and the noise, depending on the case under
 * study. yFromFiles reads a dataset formatted for the corona, while
 * randonInput generated a random x, set y = G(x), and is therefore used
 * as a possibility to check toy models' data and the algothms effectinevess*/
int yFromFile (int ignore_n_days, const char *filename, int verbose) {
        FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf ("Error: file %s not found!\n", filename);
		return 0;
	}
        int count = 0;
        double discard = 0;
        /* Ignore n lines, each of them has two numbers */
        for (int i = 0; i < 2 * ignore_n_days; ++i) {
                fscanf (file, "%lf", &discard);
        }
        /* Start reading the dataset from the day after */
        for (int i = 0; i < glob_dCod; ++i) {
                /* Ignore the first digit, it's for the plot */
                fscanf (file, "%lf", &discard);
                if (fscanf (file, "%lf", glob_y + i)) {
                        ++count;
                }
        }
	/* Set the ODE initial condition as the first read value */
	glob_initCond = glob_y[0];
	if (verbose) {
	        printf("Ignored %d days\n", ignore_n_days);
	        printf("Read: %d numbers\n", count);
		printf("Starting ODE condition: %f\n", glob_initCond);
	        printf("Read data: ");
	        printVec(glob_y, glob_dCod);
		printf("\n");
	}
	/* Now set the noise accordingly to the magnitude of the generated y*/
	for (int i = 0; i < glob_dCod; ++i) {
		glob_eta[i] = fabs(glob_y[i]);
	}
	/* Alternative option for the noise-free cases */
//	fillWith(0.1, glob_eta, glob_dCod);
	if (verbose) {
		printf("Diaginal of the noise covariance matrix:\n");
		printVec(glob_eta, glob_dCod);
		printf("\n");
	}
        fclose(file);
	return count;
}
/* Store in x, array of dimension glob_dDom, a random input for G.
 * Then set glob_y = G(u), and add some noise glob_eta depending on
 * the norm of y. */
void randomInput (double *x)
{
	assert(x != NULL);
	x[0] = rndmUniformIn(0.01, 0.60, NULL);
	x[1] = rndmUniformIn(100000, 400000, NULL);
	printf("X:\t");
	for (int i = 0; i < glob_dDom; ++i) {
		printf("%e ", x[i]);
	} printf("\n");
	/* Random initial condition for the ODE */
	glob_initCond = abs(rndmGaussian(100, 50, NULL));
	printf("Set initial ODE condition to: %f\n", glob_initCond);
	/* Now produce y = G(x) */
	G(x, glob_dDom, glob_y, glob_dCod);
	/*
	printf("G(X):\t");
	printVec(glob_y, glob_dCod);
	printf("\n");
	*/
	/* Now set the noise, we propose two options: */
	for (int i = 0; i < glob_dCod; ++i) {
		/* According to the measure's value itself */
//		glob_eta[i] = fabs(glob_y[i]);
		/* The noise-free case, for parameter tuning */
		glob_eta[i] = 0.1;
	}
	printf("[noise-free]\n");
//	printf("Noise diagonal:\n");
//	printVec(glob_eta, glob_dCod);
	/* Re-set y by addding the noise, so to have y = G(X) + eta */
	for (int i = 0; i < glob_dCod; ++i) {
		glob_y[i] += rndmGaussian(0., glob_eta[i], NULL);
//		glob_y[i] += fabs(rndmGaussian(0., glob_eta[i], NULL));
	}
	printf("G(X) + NOISE: ");
	printVec(glob_y, glob_dCod);
	printf("\n");
	/* Alternative printing for file */
//	for (int i = 0; i < glob_dCod; ++i) {
//		printf("%d %.f\n", i+1, glob_y[i]);
//	}
}

int main(int argc, char **argv) {
	srand(time(NULL));
#if TEST_MODE
	glob_dCod = floor(rndmUniformIn(10, 50, NULL));
	printf("Number of observations: %d\n", glob_dCod);
#endif
/* When NOT in test mode, must specify the range of days to read
 * and the considered dataset to be searched into ../datasets */
#if !TEST_MODE
	if (argc != 4) {
		printf("syntax: fist_day last_day country.txt\n");
		return -1;
	}
	int from_day = atoi(argv[1]); /* Minimum is day 1 */
	int to_day = atoi(argv[2]);
	glob_dCod = to_day - from_day + 1;
	char *filename = malloc(sizeof(char) * 40);
	filename[0] = '\0';
	filename = strcat (filename, "datasets/");
	filename = strcat (filename, argv[3]);
	printf("%s [%d, %d]\n", argv[3], from_day, to_day);
//	printf("Codomain of dimension %d\n", glob_dCod);
#endif

	/* Let's start with the variables completely in common between all the
	 * available methods */
	int n_samples = 1000;	/* Number of samples to generate */
	int n_iterations = 20000;/* Depth of every chain generating 1 sample */
	int centroids = 2;
	/* Covariance matrix for the Gaussian walk during the chain */
	/* This is the covariance of the prior probability measure */
	double *cov = malloc(sizeof(double) * glob_dDom * glob_dDom);
	id(cov, glob_dDom);
	cov[0] = 0.01;
	cov[glob_dDom + 1] = pow(5000, 2);
	/* 0 < beta < 1. small = conservative; hight = explorative */
	double beta = 0.4;
	/* Every chain has a starting point calculated as
	 * startMin + abs(gaussian(startMean, startCov)
	 * in this way you are free to increase/decrease the randomization
	 * of the starting points, stored in start_pt defined later */
//	double startPtCov = pow(30000, 2);
//	double startPtMin = 150000;
//	double startPtMean = 150000;
	/* Multiple reconstruction with the same data are done, to check
	 * stochastic stability */
	int tot_test = 1;
	/* The error is set to be practically zero. In the future,
	 * make this variable higher to model measuremen errors */

	/* From now on, all the code should be untouched, since the parameters
	 * are all set above */
	double *mean = malloc (sizeof(double) * glob_dDom);
	double *var  = malloc (sizeof(double) * glob_dDom);
	/* Initialize the seed: srand() is used e.g. in randomSample above,
	 * while seed_r is used e.g. when parallelization is enables. */
	//srand(time(NULL));
//	srand(1);
	unsigned int *seed_r = malloc(sizeof(unsigned int) * n_samples);
	assert(seed_r != NULL);
	seed_r[0] = time(NULL) + (unsigned) 1;
//	seed_r[0] = 2;
	for (int i = 1; i < n_samples; ++i) {
		seed_r[i] = seed_r[i - 1] + (unsigned) 1;
	}

/* In test mode we generate a random true input x, and set y = G(x). Then we
 * noise y, and try to find x again. Here initialize true_x in memory */
#if TEST_MODE
	double *true_x = malloc(sizeof(double) * glob_dDom);
	assert(true_x != NULL);
#endif
	/* Memory for global y, observed values, its preimage is the goal */
	glob_y = malloc(sizeof(double) * glob_dCod);
	assert(glob_y != NULL);
	glob_eta = malloc(sizeof(double) * glob_dCod);
	assert(glob_eta != NULL);
	/* Container of all the Monte Carlo samples */
	double *raw_samples = malloc(sizeof(double) * glob_dDom * n_samples);
	assert(raw_samples != NULL);
	/* The samples will be organized by using k-means clustering.
	 * These values are then stored into km_results */
	int max_kmeans = 500;	/* Maximum number of iteration in the kmeans */
	double *km_results = malloc(sizeof(double) * (glob_dDom+1) * centroids);
	assert(km_results != NULL);

	/* Once y is fixed, we perform the recontruction multiple times
	 * to see if the values are the same, hint to understand it MCMC 
	 * reached convergence. 
	 * For simplicity we authomatically compute the error relative
	 * to the most frequent sampled point, called MAP, storing in map_err
	 * and defined to be successfull if lower than tol_err %.
	 * In case of TEST_MODE, since x is known, this is the true error
	 * (so x_found - x_true / x_true), otherwise is the residual
	 * (G(x_found) - y / y ). */
	double tol_err = 2;	
	int success = 0;
	double avrg_acceptance_rate = 0;
	double map_err = 0;

	/* Each Markov Chain needs a starting point. In theory all equivalent,
	 * in practice it afflicts the convergence speed. Stored in start_pt */
	double *start_pt = malloc(sizeof(double) * glob_dDom * n_samples);
	assert(start_pt != NULL);

#if TEST_MODE
	/* Set the data to process:
	* When TEST_MODE, we have a random generated true_x input,
	* whose evaluation under G initialize glob_y */
	randomInput(true_x);
	printf("--- press a key to continue ---\n");
	getchar();
//	return 0;
#else	
	/* Otherwise we are in SIMULATION_MODE, where there is NO known
	 * intput true_x, since it's the one we are looking to, 
	 * and the glob_y is read from a source .txt file, 
	 * whose argoment indicates the days/lines to ignore */
	yFromFile(from_day - 1, filename, 1);
	getchar();
#endif

	/* THE RECONSTRUCTION PROCESS STARTS NOW! */
	/* Run a serie of multiple test on the SAME DATA
	 * This is CRUCIAL to check if there is convergence in probability */
	for (int a = 0; a < tot_test; ++a) {
		setbuf(stdout, NULL);
		/* Reset the contenitive variables */
		fillzero(raw_samples, glob_dDom * n_samples);
		fillzero(km_results, (glob_dDom + 1) * centroids);
#if TEST_MODE
		printf("###### TEST %d of %d ######\n", a + 1, tot_test);
#endif
		/* Determine how the chain starting points are chosen */
		fillzero(start_pt, glob_dDom * n_samples);
	        for (int i = 0; i < n_samples * glob_dDom; i += glob_dDom) {
                	start_pt[i] = rndmUniformIn(0, 0.6, NULL);
	                start_pt[i + 1] = rndmUniformIn(1000, 900000, NULL);
//				500000;
				//startMin + 
//				abs(rndmGaussian(startMean, startCov, NULL));
//				rndmGaussian(startMean, startCov, NULL);
       		}
	//	printf("Starting points initialized\n");
/*
 * 		printMat(start_pt, n_samples, glob_dDom);
		getchar();
		*/
#if TEST_MODE
		printf("Err/residual distibutions at the BEGINNING:\n");
		kMeans(start_pt, n_samples, glob_dDom, centroids, max_kmeans,
			km_results, 0);
		/* Visualize the data */
		kmnsVisual(km_results, centroids, glob_dDom);
		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, true_x, 1);
		printf("INITIAL MAP _actual_ error: %f%%\n", map_err);
		meanAndVarG (start_pt, glob_dDom, n_samples, mean, var, NULL, 
				G, glob_dCod, true_x, glob_y);
#else
//		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
//			G, glob_dCod, glob_y, NULL, 1);
//		printf("INITIAL MAP _residual_ error: %f%%\n", map_err);
//		meanAndVarG (start_pt, glob_dDom, n_samples, mean, var, confi, 
//				G, glob_dCod, NULL, glob_y);
#endif
		/* Starting points organized. Perform MCMC */
		avrg_acceptance_rate = 
#if PARALLEL /* Parallel pcn */
	        prll_uPcnSampler(potU, glob_dDom, start_pt, n_samples,
				n_iterations, raw_samples, beta, cov, seed_r,
				positive);
		//		allok);
#else		/* Ordinary pcn */
	        uPcnSampler(potU, glob_dDom, start_pt, n_samples, 
			n_iterations, raw_samples, beta, cov, 
	//		allok, 1);
			positive, 1);
#endif
		printf("\nAverage acceptance rate: %.2f%%\n", 
				avrg_acceptance_rate);
		/* Now the samples have been stored into raw_samples. 
		 * Clean them by using kmean*/
		fillzero(km_results, (glob_dDom + 1) * centroids);
		kMeans(raw_samples, n_samples, glob_dDom, centroids, max_kmeans,
			km_results, 0);
		kmnsVisual(km_results, centroids, glob_dDom);
#if TEST_MODE
		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, true_x, 1);
		printf("FINAL MAP _actual_ error: %f%%\n", map_err);
		meanAndVarG (raw_samples,glob_dDom,n_samples, mean, var, NULL,
				G, glob_dCod, true_x, glob_y);
#else
		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, NULL, 1);
		printf("res_err : %f%%\n", map_err);
		meanAndVarG (raw_samples,glob_dDom,n_samples, mean, var, NULL, 
				G, glob_dCod, NULL, glob_y);
#endif
		if (map_err < tol_err) { ++success; }
#if TEST_MODE
		printf("***end of simulation %d***\n", a + 1);
#endif
	} /* The tests have been performed */
#if TEST_MODE
	printf("%d OK out of %d\n", success, tot_test);
#endif
	free(seed_r);
	free(glob_y);
	free(glob_eta);
	free(cov);
	free(start_pt);
	free(km_results);
	free(raw_samples);
	free(mean);
	free(var);
//	free(confi);
#if TEST_MODE
	free(true_x);
#else
	free(filename);
#endif
	return 0;
}
