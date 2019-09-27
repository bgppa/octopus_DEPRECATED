#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include "RANVAR.h"
#include "BASIC.h"

void multidim_pCN(double* prior_mean, double* prior_diag_variance, void(*G)(const double*,int,double*,int), int ITER_NUM, \
					double* y, double eta, double beta, int dn, int dm, double**tmp, double* x_0, double* x_1){
		/* prior_mean and _variance are the parameters for the prior Gaussian distribution on R^dn
		 * G : the operator on which we do the bayesian inversion. G : R^dn -> R^dm
		 * ITER_NUM = number of MonteCarlo iterations
		 * y : the observed point
		 * eta : the noise variance
		 * beta : the coefficient 0 < beta < 1 described into the pCN algorithm;
		 * dn : dimension of the domain
		 * dm : dimension of the codomain.
		 * tmp is an array of pointers ALREADY initialized, used for operations in between (avoiding so to call malloc multiple times)
		 * x_1 : point in R^n which contains the result of the monte carlo iterations
		 * By repeating this procedure multiple times, the solution to the problem is given by the most frequent sample */
		// x_0 is the starting point for the Markov Chain */
		//
		/* Set the starting point: the same as the prior mean.
		this function is temporarely desabled, assuming that the user given an initial point as he will */
		/*
		for(int i=0; i<dn; ++i)
			x_0[i] = prior_mean[i];
		*/

		/* x_1 is aleady given ad an argument to the function */
		double log_alpha;
		/* I am going to use: tmp[0],1,2,3, which MUST be already initialized from main. CRUCIAL!!!! */

		for(int i = 0; i < ITER_NUM; ++i){
				/* pCM method */
				/* Sample a new proposed x_1 balancing the prior gaussian and the previous point x_0 */ 
				for(int i=0; i<dn; ++i)
						x_1[i] = beta*rndm_gaussian(prior_mean[i], prior_diag_variance[i]) + sqrt(1.0 - beta*beta) * x_0[i];

				G(x_0, dn, tmp[0], dm); /* Put in tmp[0] the evaluation of G(x_0) */
				G(x_1, dn, tmp[1], dm); /* Put in tmp[1] the evaluation of G(x_1) */

				//NEW_diff_array(y, tmp[0], dm, tmp[2]);
				/* NEW_diff_array is an old function fron linalg.c, now replaced by the
				 * combination of copy and diff in BASIC.c
				 * NEW_diff_array does: tmp[2] = y - tmp[0]
				 * (vector pointwise difference, dimension dm)
				 * To achieve the same results we do: */
				 copy(tmp[0], tmp[2], dm);
				 /* So now tmp[2] contains the value of tmp[0] */
				 diff(y, tmp[2], dm);
				 /* Which is: tmp[2] = y - tmp[2]
				  * i.e.: tmp[2] = y - tmp[0] as required */

				 /* Repeat the same reasoning for tmp[3] */
				 //NEW_diff_array(y, tmp[1], dm, tmp[3]);
				 copy(tmp[1], tmp[3], dm);
				 diff(y, tmp[3], dm);
				 
				/* So:  tmp[0] the evaluation of G in x_0
				 * 		tmp[1] the evaluation of G in x_1
				 * 		tmp[2] the (componentwise) difference between y and G(x_0)
				 * 		tmp[3] the (componentwise) difference between y and G(x_1)*/
				/* sqnorm(vector, dim) computes the norm^2 of a dim-dimensional array */
				/* DEPRECATED: now replaced by nrm2(vector,dim) */
//				log_alpha = (sqnorm(tmp[2],dm) - sqnorm(tmp[3],dm)) / (2.0 * eta * eta);
				log_alpha = (nrm2(tmp[2],dm) - nrm2(tmp[3],dm)) / (2.0 * eta * eta);

#ifdef MYDEBUG 
				printf("DEBUG...\n");
				printf("x_0: ");
				printVec(x_0, dn);
				printf("x_1: ");
				printVec(x_1, dn);
				printf("G(x_0): ");
			   	printVec(tmp[0], dm);
				printf("G(x_1): ");
			   	printVec(tmp[1], dm);
				printf("y : ");
				printVec(y, dm);
				printf("y - G(x_0) : ");
			   	printVec(tmp[2], dm);	
				printf("y - G(x_1) : ");
				printVec(tmp[3], dm);
				printf("eta: %f\n", eta);
	//			printf("|y - G(x_0)|^2 : %f\n", sqnorm(tmp[2],dm) ); // sqnorm is deprecated and 
	//			printf("|y - G(x_1)|^2 : %f\n", sqnorm(tmp[3],dm) ); // replaced by nrm2
				printf("|y - G(x_0)|^2 : %f\n", nrm2(tmp[2],dm) ); 
				printf("|y - G(x_1)|^2 : %f\n", nrm2(tmp[3],dm) ); 
			printf("log_alpha: %f\n", log_alpha);
#endif
				if (log(rndm_uniform()) <= log_alpha){
#ifdef MYDEBUG
						printf("FROM ");
						printVec(x_0,dn);
						printf("TO ");
						printVec(x_1, dn);
						printf("accepted!\n");
						getchar();
#endif
						/* The point is accepted: copy in x_0 the value from x_1
						 * and start again the cycle */
						/* Maybe, can now use the function copy, from BASIC.h? */
						for(int i=0; i<dn; ++i)
								x_0[i] = x_1[i];

				} /* end if */
#ifdef MYDEBUG
				else {
						/* The point is refused: nothing needs to be done, except possible debug output */
						printf("FROM ");
						printVec(x_0,dn);
						printf("TO ");
						printVec(x_1, dn);
						printf("refused.\n");
						getchar();
				}
#endif
		} /* End of for: at this point the vector x_1 contain a sample from the target measure */
}








/* --------------------------- NEW PART --------------------------------------*/

/* Here I am trying to re-write the Bayesian Pcn inverse technique
 * in a more controlled way.
 * The pcnMcmc routine should not check for the validity of its parameters,
 * being evoked many times during the execution.
 * Rather is task is given to checkPcnParameters.
 * So now is user's responsability to use it once
 * before running the monte carlo chain.
 */


/* The following function check the validiy of parameters
 * that are supposed to be used with pcnMcmc */
int checkPcnParameters(double* C,
                       void(*G)(const double*,int,double*,int),
                       int ITER_NUM,
                       double* y,
                       double eta,
                       double beta,
                       int dn,
                       int dm,
                       double** tmp, 
                       double* x0,
                       double* x1,
                       int verbose){

	assert(P != NULL);
	assert(C != NULL);
	assert(G != NULL);
	assert(ITER_NUM > 0);
	assert(y != NULL);
	assert(eta > 0);
	assert(beta > 0 && beta < 1);
	assert(dn > 0);
	assert(dm > 0);
	assert(tmp != NULL && *tmp != NULL);
	assert(x0 != NULL);
	assert(x1 != NULL);
	return 1;
}

/* This function performs a Monte Carlo Metropolis sampling by following
 * the pCN algorithm suitable for the Bayesian Inverse problem.
 * C : covariance matrix of the gaussian prior measure in R^dn
 * G : the operator on which we do the bayesian inversion. G : R^dn -> R^dm
 * ITER_NUM : number of steps for every MCMC
 * y : the observed points, array of dimension dm
 * eta : the noise variance
 * beta : the coefficient 0 < beta < 1 described into the pCN algorithm;
 * dn : dimension of the domain
 * dm : dimension of the codomain.
 * tmp is an array of 4 dm-dimensional arrays/pointers ALREADY initialized,
	used for some value
	schifting insise the steps (avoiding so to call malloc multiple times)
 * x0 : point in R^n on which we start the Markov Chain
 * x1 : point in R^n; output value; at the end will contain the MC result
 * verbose : integer that enables a debug mode */

void pcnMcmc(double* C,
             void(*G)(const double*,int,double*,int),
             int ITER_NUM,
             double* y,
             double eta,
             double beta,
             int dn,
             int dm,
             double** tmp,
             double* x0,
             double* x1,
             int verbose){

	double log_alpha;
	int i=0;
	int k=0;
	/* I use: tmp[0],1,2,3, MUST be already initialized from main.CRUCIAL!*/

	/* Perform ITER_NUM steps */
	for(i=0; i < ITER_NUM; ++i){
		/* New step: balance between the prior gaussian and the previous point x0 */
		/* Start by sampling x1 as a 0-mean dn-dimensional gaussian
		 * with covariance matrix C */
		rndmNdimGaussian(NULL, C, dn, x1, 0);
		/* Now balance x1 according to the weight beta and x0 */
		for(k=0; k<dn; ++k){
			x1[i] = beta*x1[i] + sqrt(1.0 - beta*beta) * x0[i];
		}

		/* Compute the potentials used in the Metropolis acceptance rate */
		G(x0, dn, tmp[0], dm); /* Put in tmp[0] the evaluation of G(x0) */
		G(x1, dn, tmp[1], dm); /* Put in tmp[1] the evaluation of G(x1) */

		/* Copy tmp[0] into tmp[2] */
		copy(tmp[0], tmp[2], dm);
		/* perform then: tmp[2] = y - tmp[2]
		 * i.e.: tmp[2] = y - tmp[0] as required */
		diff(y, tmp[2], dm);

		/* Repeat the same reasoning for tmp[3] */
		copy(tmp[1], tmp[3], dm);
		diff(y, tmp[3], dm);

		/* So:  tmp[0] evaluation of G in x0
		 *	tmp[1] evaluation of G in x1
		 * 	tmp[2] (componentwise) difference between y and G(x0)
		 * 	tmp[3] (componentwise) difference between y and G(x1)*/
		/* Compute the logarithm of the acceptance rate alpha */
		log_alpha = (nrm2(tmp[2],dm) - nrm2(tmp[3],dm)) / (2.0*eta*eta);

		if(verbose){
			printf("Verbose mode activated!\n");
			printf("x0: ");
			printVec(x0, dn);
			printf("x1: ");
			printVec(x1, dn);
			printf("G(x0): ");
			printVec(tmp[0], dm);
			printf("G(x1): ");
			printVec(tmp[1], dm);
			printf("y : ");
			printVec(y, dm);
			printf("y - G(x0) : ");
			printVec(tmp[2], dm);	
			printf("y - G(x1) : ");
			printVec(tmp[3], dm);
			printf("eta: %f\n", eta);
			printf("|y - G(x0)|^2 : %f\n", nrm2(tmp[2],dm) ); 
			printf("|y - G(x1)|^2 : %f\n", nrm2(tmp[3],dm) ); 
			printf("log_alpha: %f\n", log_alpha);
		}

		/* Accept the new point if the rate is enough */
		if (log(rndm_uniform()) <= log_alpha){
			if(verbose){
				printf("FROM ");
				printVec(x0,dn);
				printf("TO ");
				printVec(x1, dn);
				printf("accepted!\n");
			}
			/* The point is accepted: copy in x0 the value from x1
			 * and start again the cycle */
			copy(x0, x1, dn);
			/* Maybe, can now use the function copy, from BASIC.h? */
			//			for(int i=0; i<dn; ++i)
			//					x0[i] = x1[i];

		} /* end if log() <= log_alpha */
		else {
			if(verbose){
			/* Point refused: nothing to do, only verbose */
				printf("FROM ");
				printVec(x0,dn);
				printf("TO ");
				printVec(x1, dn);
				printf("refused.\n");
			}
		}
		if(verbose){ /* The verbose stops at every cycle */
			getchar();
		}
	} /* End: now x1 contains a single sample from the target measure */
}
