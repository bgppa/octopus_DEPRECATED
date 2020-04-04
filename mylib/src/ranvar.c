/* This file contains multiple functions for generating various 
 * kind of RANdom VARiables and stochastic processes. It will be
 * constantly uploaded accordingly to the author's practical
 * needs. */

/* COMMENTE ADEQUATELY the private seed, attempt to paralelizing!!! */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "myblas.h"
#include "ranvar.h"

/* Constants useful for the gaussian sampling */
#define P0 -0.322232431088
#define P1 -1
#define P2 -0.342242088547
#define P3 -0.0204231210245
#define P4 -0.0000453642210148
#define Q0 0.099348462606
#define Q1 0.588581570495
#define Q2 0.531103462366
#define Q3 0.10353775285
#define Q4 0.0038560700634

/* Uniform double in the range [0,1) */
double rndmUniform(unsigned int *private_seed)
{
        if (private_seed == NULL) {
		return (double) rand() / RAND_MAX;
        } else {
		return (double) rand_r(private_seed) / RAND_MAX;
	}
}

/* Return a random real in the interval [a, b) */
double rndmUniformIn(double a, double b, unsigned int *private_seed)
{
	assert (a < b);
	return a + rndmUniform(private_seed) * (b - a);
}

/* Exponential distribution with average lambda = lam */
double rndmExp(double lam, unsigned int *private_seed)
{
	assert(lam > 0.);
	return - log(rndmUniform(private_seed)) / lam;
}

/* Return a geometric discerete random variable in the naturals {1,..}
 * with parameter p, so having mean 1/p */
int rndmGeometricInt(double p, unsigned int *private_seed)
{
	assert(p > 0. && p < 1.);
	return floor(log(rndmUniform(private_seed)) / log(1. - p));
}

/* One-dimensional gaussian with normal and variance.
 * The followed method is decribed in
 * <Stochastic Simulation> by Amussen-Glynn */
double rndmGaussian(double mean, double variance, unsigned int *private_seed)
{
	assert(variance >= 0);
	/* Variance = 0 -> Dirac delta */
	if (variance == 0) {
		return mean;
	} else {	/* Sampler from the Knuth's textbook */
		/*
                double P0 = -0.322232431088;
                double P1 = -1;
                double P2 = -0.342242088547;
                double P3 = -0.0204231210245;
                double P4 = -0.0000453642210148;
                double Q0 = 0.099348462606;
                double Q1 = 0.588581570495;
                double Q2 = 0.531103462366;
                double Q3 = 0.10353775285;
                double Q4 = 0.0038560700634;
		*/
                double u = rndmUniform(private_seed);
                double y = sqrt(-2.0 * log(1.0 - u));
                double NUM = P0 + P1*y + P2*y*y + P3*y*y*y + P4*y*y*y*y;
                double DEN = Q0 + Q1*y + Q2*y*y + Q3*y*y*y + Q4*y*y*y*y;
                return mean + sqrt(variance) * (y + NUM/DEN);
        }
}

double old_rndmGaussian(double mean, double variance, unsigned int *private_seed)
{
	double res = -6.;
	for (int i = 0; i < 12; ++i) {
		res += rndmUniform(private_seed);
	}
	return mean + sqrt(variance) * res;
}



/* Simplified Gaussian sampling: sample a d-dimensional Gaussian having
 * DIAGONAL covariance matrix and mean x. The new values are directly
 * overwritten on x itself, in order to increase efficiency.
 * Note that it's the same of adding a zero-mean gaussian on the variable x */
void rndmDiagGauss (double* x, const double *diag_cov, int d,
		unsigned int* seed_r)
{
	assert(x != NULL && diag_cov != NULL && d > 0);
	for (int i = 0; i < d; ++i)
		x[i] += rndmGaussian(0, diag_cov[i * d + i], seed_r);
}

/* Given a d-dimensional mean, a d times d symmetric positive definite
 * matrix C, stores in res a d-dimensional gaussian sample.
 * The method followed is from the Leray's Monte Carlo book.
 * The verbose option is meant for debugging purposes.
 * Parameters:
 - an array m of dimension d representing the mean.
   If NULL, the choice is intended to be "zero mean";
 - C the covariance matrix, assumed symmetric and positive definite;
 - d is the dimension of m, and C (d x d)
 - res is a d-dimension array where the d-dimension gaussian sample
   will be written;
 - a non-zero value for verbose activates debug */
void rndmNdimGaussian(double *m, const double *C, int d,
                      double *res, unsigned int *private_seed, int verbose)
{
        assert(C != NULL);
        assert(res != NULL);
        assert(d > 0);

        int i, j, k;
        int to_free = 0;

        /* Passing a null pointer implies zero mean */
        if (m == NULL){
                m = malloc(sizeof(double) * d);
                for (i = 0; i < d; ++i) {
                        m[i] = 0;
                }
                to_free = 1;
        }
        assert(m != NULL);

        if (verbose) {
                printf("Covariance matrix:\n");
                printMat(C, d, d);
        }

        /* sig is an auxiliary dXd matrix defined according to Leray.
         * The entire full algorithm is NOT explained into details,
         * being a literal adaptation to Leray's proposition explained
         * well in the mathematical reference */
        double *sig = malloc(sizeof(double) * d * d);
        assert(sig != NULL);
        double tmp;

        /* In order to initialize sigma, set it all to -1.: just a signal; */
        for (i = 0; i < d; ++i) {
                for (j = 0; j < d; ++j) {
                        sig[i * d + j] = -1.;
                }
        }

        if (verbose) {
                printf("debug: understand if sigma is correctly initialized\n");
                printf("Set all to -1 just for watching the changes: \n");
                printMat(sig, d, d);
                getchar();
        }

        /* Start by defining the firts column */
        tmp = sqrt(C[0]); /* sqrt(C_11) in the notes */
        for (i = 0; i < d; ++i) {
                sig[i * d + 0] = C[i * d + 0] / tmp;
        }

        /* proceed defining the diagonal (i starts from 1, 0 already done) */
        for (i = 1; i < d; ++i) {
                tmp = 0;
                for(j = 0; j < i - 1; ++j) {
                        tmp += pow(C[i * d + j], 2.);
                }
                sig[i * d + i] = sqrt(C[i * d + i]) - tmp;
        }

        if (verbose) {
                printf("Defined: first column and diagonal\n");
                printMat(sig, d, d);
                getchar();
        }

        /* Set the upper tridiagonal part to zero */
        for (i = 0; i < d; ++i) {
                for (j = 0; j < d; ++j) {
                        if (j > i) {
                                sig[i * d + j] = 0.;
                        }
                }
        }

        if (verbose) {
                printf("Upper part setted to zero: \n");
                printMat(sig, d, d);
        }

        /* Contructing the sigma matrix: (explained in the book) */
        for (i = 0; i < d; ++i) {
                for (j = 0; j < d; ++j) {
                        if (j > 0 && j < i) {
                            tmp = 0;
                            for (k = 0; k < j - 1; ++k) {
                                    tmp += sig[i * d + k] * sig[j * d + k];
                            }
                            sig[i*d + j] = (C[i*d + j] - tmp) / sig[j * d + j];
                        }
                }
        }

        if (verbose) {
                printf("Final sigma matrix: \n");
                printMat(sig, d, d);
        }

        /* Construct y as a d-dimensional vector with each component normal
         * gaussian, needed for the algorithm*/
        double *y = malloc(sizeof(double) * d);
        assert(y != NULL);
        for (i = 0; i < d; ++i) {
                y[i] = rndmGaussian(0., 1., private_seed);
        }

        /* According to the algorithm, the final multidimensional
         * gaussian sample (stored in res), is given by:
         *                  res = sig * y + m                
         * (as always, we assume res to be already malloced)
         * I use routine from basics.h for linear algebra */
         axiny(sig, y, res, d);      /* now: res = sig*y */
         axpy(m, res, d, 1.);        /* res = m*1 + res */
        /* So at the end: res = m + sig*y, as requested */
        free(y);
        free(sig);
        if (to_free) { free(m); }
}


/* Simple mean and variance estimator */
int meanAndVar (double *values, int dim, int n, double *mean, double *var,
		double *conf_int)
{
	/* n values each array of length dim */
	fillzero (mean, dim);
	fillzero (var, dim);
	if (conf_int != NULL) {
	fillzero (conf_int, dim); }
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < dim; ++j) 
			{ mean[j] += values[j + i * dim]; }
	}
	for (int i = 0; i < dim; ++i) 
		{ mean[i] /= (double) n; }
	/* Ok, now the array of means is ready */
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < dim; ++j) 
			{ var[j] += pow( values[j + i * dim] - mean[j], 2.); }
	}
	for (int i = 0; i < dim; ++i) 
		{ var[i] /= (double) (n - 1); }
	if (conf_int != NULL ) {
	for (int i = 0; i < dim; ++i) 
		{ conf_int[i] = sqrt(pow(1.96, 2.) * var[i] / n); }
	}
return 1;
}

double meanAndVarG (double *values, int dim, int n, double *mean, double *var, 
			double *conf_int,
			void (*GG) (const double *, int, double *, int),
			int codim, double *true_x, double *y)
{
        /* n values each array of length dim */
        fillzero (mean, dim);
        fillzero (var, dim);
	if (conf_int != NULL) {
	fillzero (conf_int, dim); }
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j)
                        { mean[j] += values[j + i * dim]; }
        }
        for (int i = 0; i < dim; ++i)
                { mean[i] /= (double) n; }
        /* Ok, now the array of means is ready */
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j)
                        { var[j] += pow(values[j + i * dim] - mean[j], 2.); }
        }
        for (int i = 0; i < dim; ++i)
                { var[i] /= (double) (n - 1); }
	/* Now estimate the error of the mean w.r.t. y or the true_x */
	/* Estimate the 95% confidence interval for each mean point */
	if (conf_int != NULL) {
	for (int i = 0; i < dim; ++i) 
		{ conf_int[i] = sqrt(pow(1.96, 2.) * var[i] / n); }
	}
	/* Compute the G errors w.r.t. the bayesian inverse problem context */	
	double err = 0;
	if (true_x != NULL) {
		err = nrm2dist(mean, true_x, dim) / nrm2(true_x, dim) * 100.;
	} else {
		double *mean_img = malloc(sizeof(double) * codim);
		assert (mean_img != NULL);
		GG (mean, dim, mean_img, codim);
		err = nrm2dist(mean_img, y, codim) * 100. / nrm2(y, codim);
		free (mean_img);
	}
	printf("E(X): ");
	if (conf_int != NULL){
	for (int i = 0; i < dim; ++i) 
		{ printf("%e +- %e; ", mean[i], conf_int[i]); }
	}
	else printVec(mean, dim);
	printf("Var(X): ");
	printVec (var, dim);
//	printf("%s mean err: %.2f%%\n",true_x==NULL ? "_res_" : "_true_", err);
	return err;
}

/* FROM NOW ON THERE IS NO PARALELIZATION IMPLE<EMTED */
/* ----- PART 3: Stochastic processes ---- */


#if 0 /* TO IMPLEMENT AGAIN LATER */



/* The requirement for this library are given only by basic linear
 * algebra operations. It contains various sections, then described */

/* REMARK: it is responsability of the user to remember updating the
 * random seed, as commonly done with C random.c */

/* -------------------- PART 1 --------------------:

 * generical management functions, like e.g. for printing the values,
 * interacting with files, initializing data structures for representing data.
 * The crucially important structure Stproc is defined into
 * ranvar.h, but here re-written between comments in order not to
 * forget its mechanisms.
        struct Stproc{
                double* val; val[i] = value corresponding at time t[i]
                double* t;   time array
                int len;     length of both val and t, i.e. how many values
        } Stproc;

*/

/* Offer a default resetting for the Stproc structure.
 * It is not *needed*, but can help for avoiding null pointer
 * failures. For instance, if you declare a strproc structure
 * and then reset with that, you know that NULL values are
 * set as a sort of flags standing for "not-used-variable". */
void resetStproc(struct Stproc *myprocessPtr)
{
        assert(myprocessPtr != NULL);

        myprocessPtr->val = NULL;
        myprocessPtr->t = NULL;
        myprocessPtr->len = 0;
}

/* Free and reset a structure Stproc */
void freeStproc(struct Stproc *myprocessPtr)
{
        assert(myprocessPtr != NULL);

        free(myprocessPtr->val);
        free(myprocessPtr->t);
        resetStproc(myprocessPtr);
}

/* Print a sampled stochastic process */
void printStproc(struct Stproc process)
{
        for (int i = 0; i < process.len; ++i){
                printf("%e\t%e\n", process.t[i], process.val[i]);
        }
}

/* Print a sampled stochastic process to a file names filename */
void fprintStproc(const char *filename, struct Stproc process)
{
        assert(filename != NULL);

        FILE *file = fopen(filename, "w");
        if(file != NULL){
                for(int i = 0; i < process.len; ++i){
                        fprintf(file, "%e\t%e\n", process.t[i], process.val[i]);
                }
                fclose(file);
        }
        else{
                printf("*err*: fprintStproc: %s: file not found\n", filename);
        }
}





/* Simulate a Pointed poisson process and store the results into container.
 * Parameters are: 
 - lambda (intensity per unit time)
 - max_time = time horizon considered;
 - time_steps: parts in which discretize the time, minimum 1;
 - max_jump: maximum amount of jumps that are expected to happen. Can be set
   huge - no consequences - but smaller spares memory;
 - container: Stochastic Process variable where to store the results */
int rndmPointPoisson(double lam, double max_time, int time_steps,
                        int max_jumps, struct Stproc *container)
{
        if (lam <= 0){
                printf("Error: lambda must be positive\n");
                return 0;
        }
        if (time_steps < 1){
                printf("Error: at least 1 time step must be observed\n");
                return 0;
        }
        assert(container != NULL);

        /* S[0] = time for the first jump counting from t=0;
           S[i] = time between the i-1 and ith jump */
        double *S = malloc(sizeof(double) * max_jumps);
        /* T[i] = time at which happens the i-th jump */
        double *T = malloc(sizeof(double) * max_jumps);
        assert(S != NULL && T != NULL);

        /* S is known by definition to follow an exponential rule;
         * These two variables follows the relation:
         * T[i] = sum of all S[i] with smaller i.
         * so initialize them accordingly*/
        S[0] = T[0] = rndmExp(lam);
        int i;
        /* i starts from 1, being 0 the exponential above */
        for (i = 1; i < max_jumps; ++i){
                S[i] = rndmExp(lam);
                T[i] = T[i - 1] + S[i];
        }

        double h = max_time / time_steps;
        /* max_jump jumps simulated, time to save the poisson values */
        container->t = malloc(sizeof(double) * time_steps);
        container->val = malloc(sizeof(double) * time_steps);
        assert(container->t != NULL && container->val != NULL);
        container->len = time_steps;

        /* Simulate the Poisson process: */
        int n = 1;
        container->t[0] = 0;
        container->val[0] = 0;
        for (i = 1; i < time_steps; ++i){
                printf("initializing time %f\n", h * i);
                /* Time discretization is clear */
                container->t[i] = h * i;
                /* the proc values are in between the number of jumps:*/
                for (; n < max_jumps; ++n){
                        if (container->t[i] < T[0]){
                                /*
                                printf("debug:"
                                        "%f < %f\n", container->t[i], T[0]);
                                */
                                container->val[i] = 0;
                                break;
                        }
                        else {
                                if (container->t[i] >= T[n - 1] && 
                                    container->t[i] < T[n] ){
                                        /*
                                        printf("debug: %f < %f < %f ?\n",
                                         T[n-1], container->t[i], T[n]);
                                        printf("yes!\n");
                                        */
                                        container->val[i] = n;
                                        break;
                                }
                                else if (n == max_jumps - 1){
                                        printf("*ERR*: max_jumps"
                                                "%d too small.\n",  max_jumps);
                                        return 0;
                                }
                        }
                }
        }

        free(S);
        free(T);
        return 1;
}

/* Clearly, many other kind of variables must be added! */
#endif
