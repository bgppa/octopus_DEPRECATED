#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int glob_dDom = 2; /* <- full general logistic case */
double glob_initCond = -1;
int glob_dCod = 10;

/* Generalized logistic map, the most esoteric one with 6 parameters! */
double general_6 (const double * cff, double t) 
{
	/* cff is an array containing the 5 parameters I need:
	 * A, K, B, ni, Q, C, so they are with array indeces:
	 * 0, 1, 2, 3 , 4, 5 */
	return cff[0] + (cff[1] - cff[0]) /
		pow( cff[5] + cff[4] * exp(- cff[2] * t), 1./cff[3]);
}

double general_5 (const double * cff, double t) 
{
	/* cff is an array containing the 5 parameters I need:
	 * K, B, ni, Q, C, so they are with array indeces:
	 * 0, 1, 2, 3 , 4, */
	double A = 571; /* STarting number of infected chinese */
	return A + (cff[0] - A) /
		pow( cff[4] + cff[3] * exp(- cff[1] * t), 1./cff[2]);
}

double gompertz (double alpha, double N, double N0, double t)
{
	return N * exp( log(N0 / N)*exp(-alpha * t) );
}

double logistic (double alpha, double N, double N0, double t)
{
	return (N0 * exp(alpha*t)) / (1. - N0/N * (1. - exp(alpha*t)));
}

/* G operator to be inverted: since the goal is recontructing the parameters 
 * starting from 20 time observations, G does the converse:
 * a_n is a vector whose components are what need for the chosen model
 * (a_n[0] is alpha, a_n[1] is N when needed), and obs a 21-dimensional
 * array which will contain the initial condition and the next 20 results */
void G(const double* a_n, int due, double *obs, int num_obs)
{
	(int) due;
	/* The solution at time zero must be: */
	if (glob_initCond == -1) {
		printf("Forgot to set the initial condition/days!");
		getchar();
	}
	obs[0] = glob_initCond;
	/* exp inizia da zero, gli altri da 1 */
	/* __ RICORDA DI MODIFICARE IL FOR */
	for (int i = 1; i < num_obs; ++i) {
//	for (int i = 0; i < num_obs; ++i) {
//		obs[i] = gompertz (a_n[0], a_n[1], obs[0], (double) i);
		obs[i] = logistic (a_n[0], a_n[1], obs[0], (double) i);
//		obs[i] = exp(a_n[0] * (double) i + a_n[1]);
//		obs[i] = general_5 (a_n, (double) i);
	}
}
