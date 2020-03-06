#ifndef _HMC_H_
#define _HMC_H_

/* NOTE: if the potential U: R^d -> R, hmc produces samples on the doubled space
 * R^2d, since we add the momentum. If we take the *second* half of these
 * samples, we obtain the right samples from U */


/* Parameters:
 * d2 is the dimension 2 * d of the system;
 * x is the starting point
 * time_interval the fixed time interval for Verlet integration step
 * n_single_step the number of divisions of the time interval
 * M the precondition matrix
 * chain_length
 * n_samples is the number of samples
 * raw_samples an array of dimension n_samples * d2/2 containing the results.
 * Remember that only the first half of the sample is needed, dimension d2/2,
 * being so the marginal distribution relative to U.
 * Return: the average acceptance rate */
double pHmcSampler(int d2, const double *x, double time_interval,
                int n_single_step, const double *M, const double* M1,
	       double (*U) (int, const double*),int chain_length,
                int n_samples, double *raw_samples);

double prll_pHmcSampler(int d2, const double *x, double time_interval,
                int n_single_step, const double *M, const double* M1,
		double (*U) (int, const double*), int chain_length,
                int n_samples, double *raw_samples, unsigned int *prll_seed);
/* Parameters:
 * d2 is the dimension 2 * d of the system;
 * x is the starting point
 * h the fixed Verlet integration step
 * lam the intensity (see the randomized chain above)
 * M the precondition matrix
 * chain_length
 * n_samples is the number of samples
 * raw_samples an array of dimension n_samples * d2/2 containing the results.
 * Remember that only the first half of the sample is needed, dimension d2/2,
 * being so the marginal distribution relative to U.
 * Return: the average acceptance rate */
double pRanHmcSampler(int d2, const double *x, double h, double lam,
                const double *M, const double *M1,
	       double (*U) (int, const double*),	int chain_length,
                int n_samples, double *raw_samples);

double prll_pRanHmcSampler(int d2, const double *x, double h, double lam,
                const double *M, const double *M1,
		double (*U) (int, const double*), int chain_length,
                int n_samples, double *raw_samples, unsigned int *prll_seed);
#endif
