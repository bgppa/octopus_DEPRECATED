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
#define TEST_MODE 0

/* Only global variables required: eta, the noise, glob_y.
 * They will be initialized via the function randomInput,
 * then generally used in many function without altering their value */
double *glob_y;
double *glob_eta;
double *glob_cov_diag;
double *glob_gauss_means;
double *glob_Am;

/* Thanks to g.h, it is assumed to have access to the globals
 * glob_dDom, glob_dCod, const integer, and the function
 * G(int, const double*, int double) representing the operator
 * G:R^glob_dDom -> R^glob_dCod to invert. They authomatically imply the
 * potential's definition below: .
 * NOTE THAT potU is defined in a way that:
 * - log e^{-potU} = PHI, i.e. phi = potU */
double phi (int dim, const double *u);


int allok(const double *v, int a) {
return 1;
}

int positive (const double *v, int a) {
        for (int i = 0; i < a; ++i) {
                if(v[i] < 0) return 0;
        }
        return 1;
}

/* yFromFile and randomInput initialize glob_y and glob_eta, the observation
 * from where to reconstruct and the noise, depending on the case under
 * study. yFromFiles reads a dataset formatted for the corona, while
 * randonInput generated a random x, set y = G(x), and is therefore used
 * as a possibility to check toy models' data and the algothms effectinevess*/
int yFromFile (int ignore_n_days, const char *filename, int verbose); 

/* Store in x, array of dimension glob_dDom, a random input for G.
 * Then set glob_y = G(u), and add some noise glob_eta depending on
 * the norm of y. */
void randomInput (double *x);


int main(int argc, char **argv) {
	srand(time(NULL));
	/* The variable glob_dDom is declared and initialized in g.c */
	glob_cov_diag = malloc(sizeof(double) * glob_dDom);
	glob_gauss_means = malloc(sizeof(double) * glob_dDom);
	fillzero(glob_gauss_means, glob_dDom);

	glob_Am = malloc(sizeof(double) * glob_dDom);
	char *name_file_posterior = malloc(sizeof(char) * 200);
	FILE *file_posterior = NULL;
#if TEST_MODE
	glob_dCod = floor(rndmUniformIn(10, 50, NULL));
	printf("Number of observations: %d\n", glob_dCod);
#endif
/* When NOT in test mode, must specify the range of days to read
 * and the considered dataset to be searched into ../datasets */
#if !TEST_MODE
	int N_samples = 12;
	int N_iter = 12;
	double expQ = 10000;
	/* 0 < beta < 1 */
	double beta;
	if (argc == 8) {
		N_samples = atoi(argv[4]);
		N_iter = atoi(argv[5]);
		expQ = atof(argv[6]);
		beta = atof(argv[7]);
	} else {
		printf("./main from to country.txt sampl iter expQ beta\n");
		return -1;
	}
	assert(beta <= 1 && beta >= 0);
	int from_day = atoi(argv[1]); /* Minimum is day 1 */
	int to_day = atoi(argv[2]);
	glob_dCod = to_day - from_day + 1;
	char *filename = malloc(sizeof(char) * 200);
	filename[0] = '\0';
//	filename = strcat (filename, "datasets/active/");
//	filename = strcat (filename, "datasets/total/");
	filename = strcat (filename, "datasets/deceased/");
	filename = strcat (filename, argv[3]);
	printf("%s [%d, %d]\n", argv[3], from_day, to_day);
	snprintf(name_file_posterior, 50, "posteriors/%d-%d-%s",
		       from_day, to_day, argv[3]);
	printf("Posterior on file %s\n", name_file_posterior);	
//	printf("Codomain of dimension %d\n", glob_dCod);
#endif


	/* Let's start with the variables completely in common between all the
	 * available methods */
	int n_samples = pow(2., N_samples);/* Number of samples to generate */
	int n_iterations = pow(2., N_iter);
	printf("%d samples each with %d iterations\n", n_samples, n_iterations);
	/* Depth of every chain generating 1 sample */
	int centroids = 10;
	/* Covariance matrix for the Gaussian walk during the chain */
	/* This is the covariance of the prior probability measure */
	/* For the moment we assume it to be diagonal, but we keep
	 * its matrix form for future compatibility. */
	double *cov = malloc(sizeof(double) * glob_dDom * glob_dDom);
	/* The covariance matrix will be initialized later */
	id(cov, glob_dDom);
/*	0 < beta < 1. small = conservative; hight = explorative */
	int tot_test = 1;


	/* From now on, all the code should be untouched, since the parameters
	 * are all set above */
	double *mean = malloc (sizeof(double) * glob_dDom);
	double *var  = malloc (sizeof(double) * glob_dDom);
	/* Initialize the seed: srand() is used e.g. in randomSample above,
	 * while seed_r is used e.g. when parallelization is enables. */
//	srand(time(NULL));
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
	//getchar();
#else	
	/* Otherwise we are in SIMULATION_MODE, where there is NO known
	 * intput true_x, since it's the one we are looking to, 
	 * and the glob_y is read from a source .txt file, 
	 * whose argoment indicates the days/lines to ignore */
	yFromFile(from_day - 1, filename, 1);
//	getchar();
#endif

	/* THE RECONSTRUCTION PROCESS STARTS NOW! */
	/* Run a serie of multiple test on the SAME DATA
	 * This is CRUCIAL to check if there is convergence in probability */

	/* Set the mean and covariance prior */
	/* Initialize the gaussian means MANUALLY */
//	glob_gauss_means[0] = 0.1;
//	glob_gauss_means[1] = glob_y[glob_dDom-1] + 
//		(expQ - glob_y[glob_dDom-1]) / 2.;
//	glob_gauss_means[2] = 1;

	/* Initialize the covariance matrix */
//	cov[0] = pow(0.02, 2.);
//	cov[4] = pow((expQ - glob_y[glob_dDom-1]) / 4., 2.);
//	cov[8] = pow(0.3, 2.);
	cov[0] = pow(0.05, 2);
	cov[4] = pow(expQ / 2., 2.);
	cov[8] = pow(0.5, 2.);

	/* Compute the product C^-1 * m, crucial for the true pCN
	 * potential in the case of non centered Gaussian */
	/* NOW YOU DO NOT USE THAT FUNCTION, TO BE CORRECTED */ 
	for (int i = 0; i < glob_dDom; ++i) {
		glob_cov_diag[i] = 1. / cov[glob_dDom * i + i];
		glob_Am[i] = glob_cov_diag[i] * glob_gauss_means[i];
	}
//	printf("Diagonal of A: \n");
//	printVec(glob_cov_diag, glob_dDom);
//	printf("Multiplied by mean: \n");
//	printVec(glob_gauss_means, glob_dDom);
	printf("beta : %f\n", beta);
	printf("expQ: %f\n", expQ);
	printf("Let's start!\n");
	getchar();

	for (int a = 0; a < tot_test; ++a) {
		setbuf(stdout, NULL);
		/* Reset the contenitive variables */
		fillzero(raw_samples, glob_dDom * n_samples);
		fillzero(km_results, (glob_dDom + 1) * centroids);
#if TEST_MODE
		printf("###### TEST %d of %d ######\n", a + 1, tot_test);
		snprintf(name_file_posterior, 80,
			       	"posteriors/toy_test%d.txt", a);		
#endif
		file_posterior = fopen(name_file_posterior, "w");
		if (file_posterior == NULL){
		printf("Unable to write on: %s\n", name_file_posterior);
		getchar();
		} else {
			printf("Writing to file: %s\n", name_file_posterior);	
		}
	/* Determine how the chain starting points are chosen */
		fillzero(start_pt, glob_dDom * n_samples);
	//	printf("Minimum Q: %f\n", glob_y[glob_dCod-1]);
		printf("Q random in %f, %f\n", glob_y[glob_dCod-1], expQ);
	        for (int i = 0; i < n_samples * glob_dDom; i += glob_dDom) {
                	start_pt[i] = rndmUniformIn(0.01, 1., NULL);
			start_pt[i + 1] = 
		       	rndmUniformIn(glob_y[glob_dCod - 1], expQ, NULL);
			/* In the Richard case */
			if (glob_dDom == 3) {
				/*
			start_pt[i + 2] = rndmUniformIn(0.01, 0.9, NULL);
			*/ /* One oder lower, using the lower cov */
			start_pt[i + 2] = rndmUniformIn(0.01, 1., NULL);
			}
       		}
	//	printf("Initial error distribution\n");
	//	meanAndVarRes (start_pt, glob_dDom, n_samples,
//			       G, glob_dCod, glob_y);


		/* Starting points...DONE. Perform MCMC */
		avrg_acceptance_rate = 
#if PARALLEL	/* Parallel pcn */
	        prll_uPcnSampler(phi, glob_dDom, start_pt, n_samples,
				n_iterations, raw_samples, beta, cov,
					positive, seed_r);
                                        //allok, seed_r);
#else		/* Ordinary pcn */
	        uPcnSampler(phi, glob_dDom, start_pt, n_samples, 
			n_iterations, raw_samples, beta, cov, 
		//		allok, 1);
				positive, 1);
#endif
	//	printf("\nFINAL ERROR DITRIBUTION: \n");
 	//	meanAndVarRes (raw_samples, glob_dDom, n_samples,
	//		       G, glob_dCod, glob_y);

		double *expectation = malloc(sizeof(double) * glob_dDom);
		double *variance = malloc(sizeof(double) * glob_dDom);
		meanAndVar (raw_samples, glob_dDom, n_samples, expectation,
				variance);
		printf("E[parameters] = ");
		printVec(expectation, glob_dDom);
		free(expectation);
		free(variance);
		
		/* Plot the posterior distribution in the file_posterior */
		fprintMat (file_posterior, raw_samples, n_samples, glob_dDom);
		printf("\nAverage acceptance rate: %.2f%%\n", 
				avrg_acceptance_rate);

		/* The samples have been stored into raw_samples. To increase
		 * output readability, order them by using kmean */
		fillzero(km_results, (glob_dDom + 1) * centroids);
		kMeans(raw_samples, n_samples, glob_dDom, centroids, 
				max_kmeans, km_results, 0);

		/* km_results contain the centroids in the format
		 * frequence,centroid. Let's create a vector containing
		 * only the centroids */
		double *km_clean = malloc(sizeof(double) * centroids *
				glob_dDom);
		for (int i = 0; i < centroids; ++i) {
                copy(km_results + i*(glob_dDom+1) + 1, 
				km_clean + i * glob_dDom, glob_dDom);
		}
		printf("List of centroids: \n");
		printMat(km_clean, centroids, glob_dDom);
		//getchar();
 		meanAndVarRes (km_clean, glob_dDom, centroids,
			       G, glob_dCod, glob_y);
		//getchar();


		kmnsVisual(km_results, centroids, glob_dDom);
#if TEST_MODE
		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, true_x, 1);
		printf("FINAL MAP _actual_ error: %f%%\n", map_err);
#else
		map_err = kmnsBayErr(km_results, centroids, glob_dDom,
			G, glob_dCod, glob_y, NULL, 1);
		printf("res_err : %f%%\n", map_err);
#endif
		if (map_err < tol_err) { ++success; }
		fclose(file_posterior);
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
	free(name_file_posterior);
#if TEST_MODE
	free(true_x);
#else
	free(filename);
#endif
	free(glob_cov_diag);
	free(glob_gauss_means);
	free(glob_Am);
	printf("[%d samples, %d iterations]\n", n_samples, n_iterations);
	return 0;
}

double phi (int dim, const double *u)
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

	/* Add now the perturbation coming from the non-centered variable */
//	r += dot(u, glob_Am, glob_dDom);
//	r -= 0.5 * dot(glob_gauss_means, glob_Am, glob_dDom); 
	free(Gu);
	return r;
}

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
		/* EXPERIMENTAL */
		glob_eta[i] =  fabs(glob_y[i]) / 20.;
//		glob_eta[i] = ten_power(glob_eta[i]) / 2.;
	}
	/* Alternative option for the noise-free cases */
//	fillWith(0.1, glob_eta, glob_dCod);
	if (verbose) {
		printf("Diagonal of the noise covariance matrix:\n");
		printVec(glob_eta, glob_dCod);
		printf("\n");
	}
        fclose(file);
	return count;
}

void randomInput (double *x)
{
	assert(x != NULL);
	/* Dimension of x = glob_dDom */
	x[0] = rndmUniformIn(0.04, 1.0, NULL);
	x[1] = rndmUniformIn(4000, 100000, NULL);
	/* If we are in the Richard case */
	if (glob_dDom == 3) {
		x[2] = rndmUniformIn(0.01, 4.0, NULL);
	}
	printf("X:\t");
	for (int i = 0; i < glob_dDom; ++i) {
		printf("%e ", x[i]);
	} printf("\n");
	/* Random initial condition for the ODE */
	glob_initCond = rndmUniformIn(20., 100., NULL);
	printf("Set initial ODE condition to: %f\n", glob_initCond);
	/* Now produce y = G(x) */
	G(x, glob_dDom, glob_y, glob_dCod);
	printf("G(X):\t");
	printVec(glob_y, glob_dCod);
	printf("\n");
	/* Now set the noise, we propose two options: */
	for (int i = 0; i < glob_dCod; ++i) {
		/* According to the measure's value itself */
		/* EXPERIMENAL */
		glob_eta[i] = fabs(glob_y[i]);
		//glob_eta[i] = pow(ten_power(glob_y[i]), 2.);
		/* The noise-free case, for parameter tuning */
//		glob_eta[i] = 0.1;
	}
//	printf("[noise-free]\n");
	printf("Noise diagonal:\n");
	printVec(glob_eta, glob_dCod);
	/* Re-set y by addding the noise, so to have y = G(X) + eta */
	for (int i = 0; i < glob_dCod; ++i) {
		glob_y[i] += rndmGaussian(0., glob_eta[i], NULL) / 20.;
//		glob_y[i] += fabs(rndmGaussian(0., glob_eta[i], NULL));
	}
	printf("G(X) + NOISE: ");
	printVec(glob_y, glob_dCod);
	printf("\n");
	/* Alternative printing for file */
//	for (int i = 0; i < glob_dCod; ++i) {
//		printf("%d %.f\n", i+1, glob_y[i]);
//	}
//	getchar();
}
