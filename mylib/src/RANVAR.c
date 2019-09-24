/* This file contains multiple functions for generating various 
 * kind of RANdom VARiables and stochastic processes. It will be
 * constantly uploaded accordingly to the author's practical
 * needed. */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include"BASIC.h"
#include"RANVAR.h"

/* The requirement for this library are given only by basic linear
 * algebra operations. It contains various sections, then described */

/* REMARK: it is responsability of the user to remember updating the
 * random seed, as commonly done with C random.c */




/* -------------------- PART 1 --------------------:
 * generical management functions, like e.g. for printing the values,
 * interacting with files, initializing data structures for representing data.
 * The crucially important structure Stproc is defined into
 * RANVAR.h, but here re-written between comments in order not to
 * forget its mechanisms.  */
#if 0 /* Keep this structure always commented! */
	typedef struct Stproc{
		double* val; /* val[i] = value corresponding at time t[i] */
		double* t;   /* time array */
		int len;    /* length of both val and t, i.e. how many values */
    } Stproc;
#endif

/* Offer a default resetting for the Stproc structure.
 * It is not *needed*, but can help for avoiding null pointer
 * failures. For instance, if you declare a strproc structure
 * and then reset with that, you know that NULL values are
 * set as a sort of flags standing for "not-used-variable". */
void resetStproc(Stproc* myprocessPtr){
	assert(myprocessPtr != NULL);

	myprocessPtr->val = NULL;
	myprocessPtr->t = NULL;
	myprocessPtr->len = 0;
}

/* Free and reset a structure Stproc */
void freeStproc(Stproc* myprocessPtr){
	assert(myprocessPtr != NULL);

	free(myprocessPtr->val);
	free(myprocessPtr->t);
	resetStproc(myprocessPtr);
}

/* Print a sampled stochastic process */
void printStproc(Stproc process){
	for(int i=0; i<process.len; ++i){
		printf("%e\t%e\n", process.t[i], process.val[i]);
	}
}

/* Print a sampled stochastic process to a file names filename */
void fprintStproc(const char* filename, Stproc process){
	assert(filename != NULL);

	FILE* file = fopen(filename, "w");
	if(file != NULL){
		for(int i=0; i<process.len; ++i)
		fprintf(file, "%e\t%e\n", process.t[i], process.val[i]);
		fclose(file);
	}
	else{
		printf("*ERROR*: fprintStproc: %s: file not found\n", filename);
	}
}


/* ----------- PART 2 ----------- 
 * various kind of Random Variables */

/* Uniform double in the range [0,1) */
double rndm_uniform(void){
		return (double) rand() / RAND_MAX;
}

/* Return a random real in the interval [a, b) */
double rndm_uniform_in(double a, double b){
    if(a < b)
        return a + rndm_uniform()*(b-a);
    else{
        printf("rndm_uniform_in: interval [%lf, %lf)" 
		" but %lf >= %lf!\n", a, b, a, b);
        return 0;
    }
}

/* Exponential distribution with average lambda = lam */
double rndm_exp(double lam){
	if(lam > 0){
		return - log(rndm_uniform()) / lam;
	}
	else{
		printf("Error: lambda must be positive!\n");
		return 0;
	}
}

/* One-dimensional gaussian with normal and variance.
* The followed method is decribed in <Stochastic Simulation> by Amussen-Glynn */
double rndm_gaussian(double mean, double variance){
	if(variance > 0){
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
		double u = rndm_uniform();
		double y = sqrt(-2.0 * log(1.0 - u));
		double NUM = P0 + P1*y + P2*y*y + P3*y*y*y + P4*y*y*y*y;
		double DEN = Q0 + Q1*y + Q2*y*y + Q3*y*y*y + Q4*y*y*y*y;
		return mean + sqrt(variance) * (y + NUM/DEN);
	}
	else{
		printf("Error: gaussian variance must be positive!\n");
		return 0;
	}
}

/* Given a d-dimensional mean, a d times d symmetric positive definite
 * matrix C, stores in res a d-dimensional gaussian sample.
 * The method followed is from the Leray's Monte Carlo book.
 * The verbose option is meant for debugging purposes */
void ndim_gaussian(double* m, double* C, int d, double* res, int verbose){
	assert(m != NULL);
	assert(C != NULL);
	assert(res != NULL);
	assert(d > 0);

	if(verbose){
		printf("Covariance matrix:\n");
		printMat(C, d, d);
	}

	/* sig is an auxiliary dXd matrix defined according to Leray.
	 * The entire full algorithm is NOT explained into details,
	 * being a literal adaptation to Leray's proposition explained
	 * well in the mathematical reference */
	double* sig = malloc(sizeof(double)*d*d);
	int i, j, k;
	double tmp;
	/* In order to initialize sigma, set it all to -1.: just a signal; */
	for(i=0; i<d; ++i){
		for(j=0; j<d; ++j){
			sig[i*d+j] = -1.;
		}
    	}

	if(verbose){
		printf("debug: understand if sigma is correctly initialized\n");
		printf("Set all to -1 just for watching the changes: \n");
		printMat(sig, d, d);
		getchar();
	}

	/* Start by defining the firts column */
	tmp = sqrt(C[0]); /* sqrt(C_11) in the notes */
	for(i=0; i<d; ++i){
		sig[i*d + 0] = C[i*d + 0] / tmp;
	}

	/* proceed defining the diagonal */
	for(i=1; i<d; ++i){
		tmp = 0;
		for(j=0; j<i-1; ++j){
			tmp += pow(C[i*d + j], 2.);
		}
		sig[i*d + i] = sqrt(C[i*d+i]) - tmp;
	}

	if(verbose){
		printf("Defined: first column and diagonal\n");
		printMat(sig, d, d);
		getchar();
	}

	for(i=0; i<d; ++i){
		for(j=0; j<d; ++j){
			if(j > i){
				sig[i*d+j] = 0.;
			}
		}
	}

	if(verbose){
		printf("Upper part setted to zero: \n");
		printMat(sig, d, d);
	}

	/* Contructing the sigma matrix: */
	for(i=0; i<d; ++i){
		for(j=0; j<d; ++j){
			if(j > 0 && j < i){
				tmp = 0;
				for(k=0; k<j-1; ++k){
					tmp += sig[i*d + k] * sig[j*d + k];
				}
			sig[i*d + j] = (C[i*d + j] - tmp) / sig[j*d + j];
			}
		}
	}

	if(verbose){
		printf("Final sigma matrix: \n");
		printMat(sig, d, d);
	}


	/* Construct y as a d-dimensional vector with each component normal
	 * gaussian, needed for the algorithm*/
	double*y = malloc(sizeof(double)*d);
	for(i=0; i<d; ++i){
		y[i] = rndm_gaussian(0., 1.);
	}

	/* According to the algorithm, the final multidimensional
	 * gaussian sample (stored in res), is given by:
	 *                  res = sig*y + m                
	 * (as always, we assume res to be already malloced)
	 * I use routine from BASIC.h for linear algebra */
	 axiny(sig, y, res, d);      /* now: res = sig*y */
	 axpy(m, res, d, 1.);        /* res = m*1 + res */
	/* So at the end: res = m + sig*y, as requested */
	free(y);
	free(sig);
}


/* ----- PART 3: Stochastic processes ---- */

/* Simulate a Pointed poisson process and store the results into container.
 * Parameters are: lambda (intensity per unit time)
 * max_time = time horizon considered;
 * time_steps: parts in which discretize the time, minimum 1;
 * max_jump: maximum amount of jumps that are expected to happen. Can be set
 * huge - no consequences - but smaller spares memory;
 * container: Stochastic Process variable where to store the results */
int rndm_point_poisson(double lam, double max_time, int time_steps, int max_jumps, Stproc* container){
	if(lam <= 0){
		printf("Error: lambda must be positive\n");
		return 0;
	}
	if(time_steps < 1){
		printf("Error: at least 1 time step must be observed\n");
		return 0;
	}
	assert(container != NULL);

	/* S[0] = time for the first jump counting from t=0;
 	S[i] = time between the i-1 and ith jump */
	double* S = malloc(sizeof(double) * max_jumps);
	/* T[i] = time at which happens the i-th jump */
	double* T = malloc(sizeof(double) * max_jumps);

	if(T == NULL || S == NULL){
		printf("Posson process: malloc failed\n");
		return 0;
	}
	else{
	/* S is known by definition to follow an exponential rule;
	 * These two variables follows the relation:
	 * T[i] = sum of all S[i] with smaller i.
         * so initialize them accordingly*/
	S[0] = T[0] = rndm_exp(lam);
	int i;
	for(i=1; i<max_jumps; ++i){
		S[i] = rndm_exp(lam);
		T[i] = T[i-1] + S[i];
	}

	double h = max_time / time_steps;
	/* max_jump jumps simulated, time to save the poisson values */
	container->t = malloc(sizeof(double)*time_steps);
	container->val = malloc(sizeof(double)*time_steps);
	container->len = time_steps;

	/* Simulate the Poisson process: */
	int n = 1;
	container->t[0] = 0;
	container->val[0] = 0;
	for(i=1; i<time_steps; ++i){
		printf("initializing time %f\n", h*i);
	        /* Time discretization is clear */
		container->t[i] = h*i;
		/* the proc values are in between the number of jumps:*/
		for(; n<max_jumps; ++n){
			if(container->t[i] < T[0]){
				printf("debug: %f < %f\n", \
					 container->t[i], T[0]);
				container->val[i] = 0;
				break;
			}
			else{
				if(container->t[i] >= T[n-1] && \
					container->t[i] < T[n] ){
					printf("debug: %f < %f < %f ?\n",\
					 T[n-1], container->t[i], T[n]);
					printf("yes!\n");
               				container->val[i] = n;// T[n-1];
					break;
               			}
				else if(n == max_jumps-1){
					printf("*ERROR: max_jumps %d too small.\n",  max_jumps);
					return 0;
               			}
           		}
       		}
    	}
	free(S);
	free(T);
	return 1;
	}
}


/* TO ADD: 1-Dimensional SDEs */

/* Clearly, many other kind of variables must be added! */
