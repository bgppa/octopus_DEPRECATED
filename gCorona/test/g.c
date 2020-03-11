#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int glob_dDom = 2;
int glob_dCod = 21;

double gompertz (double alpha, double N, double N0, double t)
{
	return N * exp(log(N0 / N) * exp(-alpha * t));
}

double logistic (double alpha, double N, double N0, double t)
{
	return (N0 * exp(alpha*t)) / (1. - N0/N * (1. - exp(alpha*t)));
}

void G(const double* a_n, int due, double *obs, int num_obs)
{
	(int) due;
	obs[0] = 3;
	for (int i = 1; i < num_obs; ++i) {
		obs[i] = gompertz (a_n[0], a_n[1], obs[0], (double) i);
//		obs[i] = logistic (a_n[0], a_n[1], obs[0], (double) i);
	}
}
