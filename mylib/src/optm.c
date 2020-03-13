/* Basic optimization library for a target function U : R^d -> R */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "myblas.h"
#include "ranvar.h"
#include "ode.h"
#include "optm.h"

/* Using ode.h to have access to -gradient of a function */

/* Find the minimum of U via gradient descent */
double gradDesc (double *x_n, int d, double (*U) (int, const double*),
		const double* lam, int iter_max, int verbose)
{
	double *tmp = malloc(sizeof(double) * d);
	assert (tmp != NULL);
	int i;
	for (i = 0; i < iter_max; ++i) {
		/* tmp = - gradU(x_n) */
		minus_gradient_U (d, x_n, tmp, U);
		if (verbose) {
			printf("\nx_%d ", i);
			printVec(x_n, d);
			printf("Grad x_%d ", i);
			printVec(tmp, d);
		}
		for (int j = 0; j < d; ++j) {
			x_n[j] += lam[j] * tmp[j];
		}
		if (verbose) {
			printf("next: ");
			printVec(x_n, d);
			getchar();
		}
	}
	if (verbose) {
		printf("gradDesc: %d iterations\n", i);
	}
	free(tmp);
	return U(d, x_n);
}

/* Find the minimum of U via a straighforward random search */
double rwMinimum (double* x_n, int d, double (*U) (int, const double*),
			const double *cov, int iter_max)
{
	double curr_min = U(d, x_n);
	double *tmp = malloc(sizeof(double) * d);
	assert(tmp != NULL);
	for (int i = 0; i < iter_max; ++i) {
		rndmNdimGaussian (x_n, cov, d, tmp, NULL, 0);
	       if (U(d, tmp) < curr_min) {
	       		copy (tmp, x_n, d);
	 		curr_min = U(d, x_n);		
	       }
	}
	free(tmp);
	return U(d, x_n);
}

/* Find zeroes of U via Newton methods */
double nwtnZero (double* x_n, int d, double (*U) (int, const double*),
		int iter_max, double tol)
{
	double *tmp = malloc(sizeof(double) * d);
	double ux;
	double err = U(d, x_n);
	int i;
	for (i = 0; i < iter_max && err > tol; ++i) {
		ux = U(d, x_n);
	       	minus_gradient_U (d, x_n, tmp, U);
		for (int j = 0; j < d; ++j) {
			x_n[j] += ux / tmp[j];
		}
		err = fabs(U(d, x_n));
	 }
	printf("%d iterations\n", i);
	if (i == iter_max - 1) {
		printf("*warning* : newton: max iteration reached\n");
	}
	free(tmp);
	return U(d, x_n);
}	



/* Find the minimum of U via Monte Carlo methods */
/* To implement */

