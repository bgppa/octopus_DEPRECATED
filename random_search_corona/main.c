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

#define PARALLEL 1
#define TEST_MODE 0

/* Only global variables required: eta, the noise, glob_y.
 * They will be initialized via the function randomInput,
 * then generally used in many function without altering their value */
double *glob_y;
double *glob_eta;


/* potU is the FUNCTION TO MAXIMIZE */
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
		/* EXPERIMENTAL */
		glob_eta[i] =  fabs(glob_y[i]);
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
/* Store in x, array of dimension glob_dDom, a random input for G.
 * Then set glob_y = G(u), and add some noise glob_eta depending on
 * the norm of y. */
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
//	getchar();
}

double exhaustive_search (double qMin, double qMax, double qN, double QMin,
				double QMax, double QN, double vMin,
				double vMax, double vN) {
	double res = 0.;
	double prms[3] = {0.};
	double sol[3] = {0.};
	double qh = (qMax - qMin) / qN;
	double Qh = (QMax - QMin) / QN;
	double vh = (vMax - vMin) / vN;
	for (double i = qMin; i < qMax; i += qh) {
		for (double j = QMin; j < QMax; j += qh) {
			for (double k = vMin; k < vMax; k += vh) {
				prms[0] = i;
				prms[1] = j;
				prms[2] = k;
				double ev = potU(glob_dDom, prms);
//				printf("%f %f %f: %f\n", i, j, k, ev);
				if(ev > res) {
					res = ev;
					copy(prms, sol, 3);
				}
//				getchar();
			}
		}
	}
	printf("maximum is %f in %f, %f, %f\n", res, sol[0], sol[1], sol[2]);
	return res;
}

double random_search (double qMin, double qMax,
			double QMin, double QMax, 
			double vMin, double vMax, int N, double *covDiag) {
	double res = 0;
	double prms[3] = {0.};
	double sol[3] = {0.};
	double ev = 0;
	/* First guess for the minimum */
	sol[0] = (qMax - qMin) / 2.;
	sol[1] = (QMax - QMin) / 2.;
	sol[2] = (vMin - vMax) / 2.;
	res = potU(glob_dDom, sol) + 0.5 * dot(covDiag, sol, 3);

	for (int i = 0; i < N; ++i) {
		prms[0] = rndmUniformIn(qMin, qMax, NULL);
		prms[1] = rndmUniformIn(QMin, QMax, NULL);
		prms[2] = rndmUniformIn(vMin, vMax, NULL);
		ev = potU(glob_dDom, prms) + 0.5 * dot(covDiag, sol, 3);
//		printf("%f %f %f: %f\n", i, j, k, ev);
		if(ev < res) {
			res = ev;
			copy(prms, sol, 3);
		}
//				getchar();
	}
	printf("minimum in %f in %f, %f, %f\n", res, sol[0], sol[1], sol[2]);
	return res;
}

int main(int argc, char **argv) {
	srand(time(NULL));
	char *name_file_posterior = malloc(sizeof(char) * 200);
	FILE *file_posterior = NULL;
	int tot_tests = 5;
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
	char *filename = malloc(sizeof(char) * 200);
	filename[0] = '\0';
	filename = strcat (filename, "../corona_modeling/datasets/deceased/");
	filename = strcat (filename, argv[3]);
	printf("%s [%d, %d]\n", argv[3], from_day, to_day);
	snprintf(name_file_posterior, 100, "posteriors-%d-%d-%s",
		       from_day, to_day, argv[3]);
	printf("Posterior on file %s\n", name_file_posterior);	
//	printf("Codomain of dimension %d\n", glob_dCod);
#endif
	/* Memory for global y, observed values, its preimage is the goal */
	glob_y = malloc(sizeof(double) * glob_dCod);
	assert(glob_y != NULL);
	glob_eta = malloc(sizeof(double) * glob_dCod);
	assert(glob_eta != NULL);
	double *true_x = malloc(sizeof(double) * glob_dDom);

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
	double cov1 = pow(4.2, 2);
	double cov2 = pow(9, 2);
	double cov3 = pow(10.2, 2);

	double cov_diag[3];
	cov_diag[0] = 1. / sqrt(cov1);
	cov_diag[1] = 1. / sqrt(cov2);
	cov_diag[2] = 1. / sqrt(cov3);

//exhaustive_search (0.01, 2.0, 1000, 7000, 500000, 500000, 0.01, 5.0, 10000);
for (int i = 0; i < tot_tests; ++i) {
random_search (0.01, 2.0, 1000, 5000, 0.01, 5.0, 10000000, cov_diag);
}

	/* THE RECONSTRUCTION PROCESS STARTS NOW! */
	return 0;
}
