#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int glob_dDom = 2;
int glob_dCod = 50;

double gompertz (double alpha, double N, double N0, double t)
{
	return N * exp( log(N0 / N)*exp(-alpha * t) );
}

double logistic (double alpha, double N, double N0, double t)
{
	return (N0 * exp(alpha*t)) / (1. - N0/N * (1. - exp(alpha*t)));
}

/* Log maps in order to avoid the appearence of big numbers */
double log_gompertz (double alpha, double N, double N0, double t)
{
	return log(N) + log(N0/N)*exp(-alpha * t);
}

double log_exponential (double alpha, double N0, double t)
{
	return log(N0) + alpha*t;
}

double log_logistic (double alpha, double N, double N0, double t)
{
	return log(N0) + alpha*t - log(1.-N0/N) - 
		log(1.- N0/N * (1. - exp(alpha*t)));
}

/* G operator to be inverted: since the goal is recontructing the parameters 
 * starting from 20 time observations, G does the converse:
 * a_n is a vector whose components are what need for the chosen model
 * (a_n[0] is alpha, a_n[1] is N when needed), and obs a 21-dimensional
 * array which will contain the initial condition and the next 20 results */
void G(const double* a_n, int due, double *obs, int num_obs)
{
	(int) due;
	obs[0] = 3.;
	for (int i = 1; i < num_obs; ++i) {
		obs[i] = gompertz (a_n[0], a_n[1], obs[0], (double) i);
//		obs[i] = logistic (a_n[0], a_n[1], obs[0], (double) i);
//		obs[i] = obs[0] * exp (sqrt(a_n[0] * (double) i));
//              obs[i] = obs[0] * exp (sqrt(a_n[0] * (double) i + a_n[1]));
//		obs[i] = log_gompertz (a_n[0], a_n[1], obs[0], (double) i);
	//	obs[i] = log_logistic (a_n[0], a_n[1], obs[0], (double) i);
        printf("%d %.f\n", i, obs[i]);
	}
}

int main() {
       double params[2] = {0.095, 51600};
        double t[glob_dCod];
        G(params, 2, t, glob_dCod);
        return 0;
} 
