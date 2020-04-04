#include <math.h>
#include <stdlib.h>

#ifndef _RANVAR_
#define _RANVAR_

#if 0
/* --- PART 1: structures, initialization, I/O --- */

/* A structure used for representing a stochastic process */
struct Stproc{
    double *val; /* val[i] contains the value corresponding at time t[i] */
    double *t;	 /* t represents the time */
    int len;	 /* length of both val and t, i.e. how many values */
} Stproc;

/* Offer a default resetting for the Stproc structure.
 * It is not *needed*, but can help for avoiding null pointer
 * failures. For instance, if you declare a strproc structure
 * and then reset with that, you know that NULL values are
 * set as a sort of flags standing for "not-used-variable". */
void resetStproc(struct Stproc *myprocessPtr);

/* Free and reset a structure Stproc */
void freeStproc(struct Stproc *myprocessPtr);

/* Print a sampled stochastic process */
void printStproc(struct Stproc process);

/* Print a sampled stochastic process to a file names filename */
void fprintStproc(const char *filename, struct Stproc process);
#endif

/* ---- PART 2: random variables --- */

/* Random double in [0, 1) */
double rndmUniform(unsigned int *);

/* Return a random real in the interval [a, b) */
double rndmUniformIn(double a, double b, unsigned int *);

/* Exponential distribution with average lambda = lam */
double rndmExp(double lam, unsigned int *);

/* Sample according a geometrical distribution with parameter p */
int rndmGeometricInt(double p, unsigned int *);

/* One-dimensional gaussian with normal and variance.
 * The followed method is decribed in
 * <Stochastic Simulation> by Amussen-Glynn*/
double rndmGaussian(double mean, double variance, unsigned int *);

/* Simplified Gaussian sampling: sample a d-dimensional Gaussian having
 * DIAGONAL covariance matrix and mean x. The new values are directly
 * overwritten on x itself, in order to increase efficiency.
 * Note that it's the same of adding a zero-mean gaussian on the variable x */
void rndmDiagGauss (double *x, const double *diag_cov, int d,
		unsigned int *);

/* Given a d-dimensional mean m, a covariance matrix A 
 * of dimension d * d, produces a d-dim gaussian sample stored in res.
 * Passing a NULL m implies zero mean */ 
void rndmNdimGaussian(double *m, const double *A,
                        int d, double *res, unsigned int *, int verbose);

int meanAndVar (double *values, int dim, int n, double *mean, double *var,
		double *conf_int);

double meanAndVarG (double *values, int dim, int n, double *mean, double *var,
			double *conf_int,
                        void (*GG) (const double *, int, double *, int),
                        int codim, double *true_x, double *y);


/* ---- PART 3: Stochastic processes ---- */

/* Simulate a Pointed poisson process and store the results into container.
 * Parameters are: lambda (intensity per unit time)
 * max_time = time horizon considered;
 * time_steps: parts in which discretize the time, minimum 1;
 * max_jump: maximum amount of jumps that are expected to happen. Can be set
 * huge - no consequences - but smaller spares memory;
 * container: Stochastic Process variable where to store the results */ 
/*
int rndmPointPoisson(double lam, double max_time, int time_steps,
                        int max_jumps, struct  Stproc* container);
*/

/* Clearly, many other kind of variables must be added! */
#endif
