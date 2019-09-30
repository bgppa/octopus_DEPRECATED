#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include "RANVAR.h"
#include "BASIC.h"
#include "KMEAN.h"

/* Let's try to completely remove the old pCN implementation
 * in favor of a more complete one */

#if 0
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

#endif







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
                       double* x1){

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

void pcnMcmc(const double* C,
             void(*G)(const double*,int,double*,int),
             int ITER_NUM,
             const double* y,
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
		/* Key rule to keep in mind:
		 * x0 represents the previous point in every step,
		 * while x1 is the new proposal
		 * The preposed x1 follows the rules here discribed: 
		 * 1) start by sampling x1 as a 0-mean (so, NULL) 
		 * dn-dimensional gaussian
		 * with covariance matrix C. No verbose mode (0 last param) */
		rndmNdimGaussian(NULL, C, dn, x1, 0);
		/* 2) balance x1 w.r.t. previous x0 according to the weight beta */
		for(k=0; k<dn; ++k){
			x1[i] = beta*x1[i] + sqrt(1.0 - beta*beta) * x0[i];
		}

		/* Compute the potentials used in the Metropolis acceptance rate,
		 * whose results determine the acceptance of x1 */

		/* Put in tmp[0] the evaluation of G(x0) */
		G(x0, dn, tmp[0], dm);
		/* Copy tmp[0] into tmp[2] */
		copy(tmp[0], tmp[2], dm);
		/* Perform then: tmp[2] = y - tmp[2]
		 * i.e.: tmp[2] = y - tmp[0] = y - G(x0) as required */
		diff(y, tmp[2], dm);

		/* Put in tmp[1] the evaluation of G(x1) */
		G(x1, dn, tmp[1], dm);	
		/* Repeat the same reasoning as before, with tmp[3] */
		copy(tmp[1], tmp[3], dm);
		/* So now tmp[3] = y - G(x1) */
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
		if (log(rndmUniform()) <= log_alpha){
			if(verbose){
				printf("FROM ");
				printVec(x0,dn);
				printf("TO ");
				printVec(x1, dn);
				printf("accepted!\n");
			}
			/* The point is accepted: copy in x0 the value of x1,
			 * since it becomes now the new starting point.
			 * Start then the cycle again */
			copy(x0, x1, dn);
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


/* Try to implement an automatized Bayesian Inversion algorithm */
/* Key point: every pcnMcmc is now repeated multiple times,
 * the resulting distribution is the posterior distribution.
 * The most frequent point is returned as a solution, while, if a name
 * file is specified, this posterior distributon is written on a file */

double bayesianInversion(int SAMPLES,
			int MCMC_ITER,
			double* MAP,
			const double* true_params,
			void(*operator)(const double*,int,double*,int),
			const double* observed_data,
			int domain_dim,
			int codomain_dim,
			double noise_var,
			double beta,
			const double* covariance_step,
			double* starting_point,
			FILE* posterior_file,
			int verbose){

	int i = 0;
	
	/* Four tmp are used for technical switces during pcn monte carlo */
	double** tmp_for_mcmc = malloc(sizeof(double*) * 4);
	for(i=0; i<4; ++i){
		tmp_for_mcmc[i] = malloc(sizeof(double) * codomain_dim);
	}
	
	/* Here I'll write all the results of my pcn MCMC
	 * total amount of points = SAMPLES
	 * each of dimension domain_dim */
	double* posterior_points = malloc(sizeof(double)*SAMPLES*domain_dim);
	/* Omit now the check */

	/* Sample from MCMC a number of times equal to SAMPLES */
	for(i=0; i<SAMPLES; ++i){
		pcnMcmc(covariance_step,
			operator,
			MCMC_ITER,
		        observed_data,
			noise_var,
			beta,
			domain_dim,
			codomain_dim,
			tmp_for_mcmc,
			starting_point,
			posterior_points+i*domain_dim,
			verbose);
	}

	/* Ok, now posterior_points should contain the list of all
	 * the sampled points! */
	
	if(verbose){
		printf("--- %d samples have been generated --- \n", SAMPLES);
		getchar();
		printMat(posterior_points, SAMPLES, domain_dim);
		printf(" - - - - - - - - - - - - - - - - - - - - \n");
		getchar();	
	}

	/* Method 1: use the python script */
	/* --- TO IMPLEMENT ---- */
	
	/* Method 2: use directly my C k-means function to elaborate data */

	int centroid_num = (int) sqrt(SAMPLES);
	int max_iteration_for_kmeans = 1000;
	kMean(posterior_points, SAMPLES, domain_dim, centroid_num,  
		posterior_file, max_iteration_for_kmeans, MAP);

	/* Ok, now the posterior with frequncies has been written to
	 * the given file (directly with frequencies), and the MAP
	 * estimator has been saved into MAP */
	if(verbose){
		printf("\nEstimated MAP: \n");
		printVec(MAP, domain_dim);
	}
	
	
	/* If the user known the true parameters, e.g. he is working with
	 * toy-model data, the true relative error can be computed */
	if(true_params != NULL){
		double err = nrm2dist(MAP, true_params, domain_dim) * 100.;
		printf("ERR: %.3f%%\n", err / nrm2(true_params, domain_dim) );
	}
	
	double* MAP_output = malloc(sizeof(double) * codomain_dim);
	operator(MAP, domain_dim, MAP_output, codomain_dim);

	/* So MAP_output contains the output generated by the operator when
	 * the MAP, i.e. the most frequent parameters estimated with the
	 * bayesian technique, are set as input. It can be used to compute
	 * the residual error: */
	double res = nrm2dist(MAP_output, observed_data, codomain_dim) * 100.;
	printf("RES: %.3f%%\n", res / nrm2(observed_data, codomain_dim) );

	free(MAP_output);
	for(i=0; i<4; ++i){
		free(tmp_for_mcmc[i]);
	}
	free(tmp_for_mcmc);

	return res;
}

	
#if 0
int readPoints(char* file_name, double* list_of_points, int dim_each, int lines){
        FILE* file = fopen(file_name, "r");
        if(file == NULL){
                printf("*error* unable to open %s\n", file_name);
                return 0;
        }
        else{
                int n=0;
                for(n=0; n < dim_each*lines; ++n){
                        fscanf(file, "%lf", list_of_points+n);
                }
                fclose(file);
                return 1;
        }
}
#endif
