#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//int glob_dDom = 2; /* <- full general logistic case */
int glob_dDom = 2; /* Testing the general logistic law */
double glob_initCond = -1;
int glob_dCod = 10;


double richard (double alpha, double N, double N0, double ni, double t)
{
	double Q = - 1. + pow((N / N0), ni);
	return N / pow(1. + Q * exp(-alpha * ni * t), 1. / ni);
}

double gompertz (double alpha, double N, double N0, double t)
{
	return N * exp( log(N0 / N)*exp(-alpha * t) );
}

double logistic (double alpha, double N, double N0, double t)
{
	return (N0 * exp(alpha*t)) / (1. - N0/N * (1. - exp(alpha*t)));
}

#if 0
double exGomp (double alpha, double N, double N0, double a, double S, double t)
{
	double T = 1. / a * log(S / N0);
	if (t < T) {
		return N0 * exp (a * t);
	} else {
		return gompertz (alpha, N, S, (double) t);
	}
}
#endif

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
		obs[i] = gompertz (a_n[0], a_n[1], obs[0], (double) i);
	//	obs[i] = logistic (a_n[0], a_n[1], obs[0], (double) i);
//	obs[i] = exGomp (a_n[0], a_n[1], obs[0], a_n[2], a_n[3], (double) i);
//		obs[i] = exp(a_n[0] * (double) i + a_n[1]);
//	obs[i] = richard (a_n[0], a_n[1], obs[0], a_n[2], (double) i);
	}
}
