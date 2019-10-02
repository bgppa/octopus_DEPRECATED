/* Trying to deduce the heat equation starting condition from the
 * solution at time 0.1 */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>
#include<math.h>
#include"BASIC.h"
#include"RANVAR.h"
#include"PCNINV.h"
#include"FILEIO.h"


double glb_time_observed = 0.01;

/* Eigenalues and eigenvectors of laplacian -dx/dx2 on [0,1]
 * with 0-boundary conditions
 * be careful with indeces: 0 in array, 1 in mathematics */
double alpha(int j){
	++j;
	return pow((2. * M_PI * j), 2.);
}

double phi(int j, double x){
	++j;
	return sin (2. * M_PI * j * x);
}

void solver(const double* a, int basis_expansion, double*y, int obs_number){
	/* Set h as the space-step in [0,1] in a way to perform a number
	 * of observations equal to obs_number.
	 * Eg: if obs_number= 11, h is then 0.1, allowing a total of
	 * 11 observations: 0, 0.1, 0.2, ... , 1. 
	 * Codomain must always be strictly higher that 1; reflecing the idea 
	 * that the boundary 0 is always observed (and should be 0; it would be
	 * a numerically confirm confirm). 
	 * Again, e.g. obs_number = 2 produces h = 1, then observing
	 * only 0 at 1 - boundary points */

	/* Basis expansion: >= 1;
	 * obs_number: > 1;
	 * y, a : non-NULL */
	
	double time_observed = glb_time_observed;	
	double h = time_observed / (obs_number - 1);
	double tmp_sum = 0;
	int i = 0;
	int j = 0;
	for(i=0; i < obs_number; ++i){
		tmp_sum = 0;
		for(j = 0; j < basis_expansion; ++j){
//			printf("a[%d] = %.3f\n", j, a[j]);
//			printf("alpha(%d) = %.3f\n", j, alpha(j));
//			printf("phi(%d, %.2f) = %.3f\n", j, h*i, phi(j, h*i));
			tmp_sum += a[j] * exp(-alpha(j)*time_observed) * phi(j, h * i);
//			printf("tmp_sum = %.3f\n", tmp_sum);
//			getchar();
		}
		y[i] = tmp_sum;
//		printf("Sol. point %.2f = %.3f\n", i*h, y[i]);
//		getchar();
	}
}

/* The user wants to produce toy-model data in order to test the inversion's
 * effectiveness. This function:
 * takes a noise intensity;
 * randomizes an array of basis coefficients a, that can later be re-read
 * in order to quantify the actual error;
 * solves to the euqation and does the observations by using "solver" and
 * storing the results into y;
 * finally, y is noised according to the previous given intensity */
void createToyData(double noise,
		double* a, int basis_expansion, double*y, int obs_number){
	int i=0;
	/* Randomize the parameters a */
	for(i=0; i<basis_expansion; ++i){
		a[i] = rndmUniformIn(-10, 10);
	}

	/* Solve the equation storing the results in y */
	solver( (const double*) a, basis_expansion, y, obs_number);

	printf("\n** noise-free obs: \n");
	printVec(y, obs_number);
	printf("\n");

	printf("(temporarely, we work in the noise-free case)\n");
	
#if 1
	/* Put a noise on each value of y */
	for(i=0; i<obs_number; ++i){
		y[i] += rndmGaussian(0, noise);
	}
#endif

	/* The "result" of this function is given by the coefficients
	 * now stored in a, and the noised solution in y which can be used
	 * as observations for reconstruction purposes */
}	
	


int main(int argc, char* argv[]){
	srand(time(NULL));
//	srand(2);

	double data_noise = 1e-4; /* Temporarely disabled in createToyData */
	double mcmc_noise = 1e-3;

	int n2 = 10;
	int mcmc2 = 12;
	int expansion = 3;
	int num_observations = 11;

	if(argc == 5 ){
		/* Then expansion and num_observations are given
		 * in input */
		expansion = atoi(argv[3]);
		num_observations = atoi(argv[4]);
	}

	if(argc >= 3){
		n2 = atoi(argv[1]);
		mcmc2 = atoi(argv[2]);
	}		
		
	/* Arguments? later */
	double* true_params = malloc(sizeof(double) * expansion);
	double* observed = malloc(sizeof(double) * num_observations);

	createToyData(data_noise, true_params, expansion,
			observed, num_observations);

	printf("** true coeff: \n");
	printVec(true_params, expansion);
	printf("\n** noised obs: \n");
	printVec(observed, num_observations);	
	printf("\n");

	/* Now that the data are ready, set the bayes parameters */
	/* Output file where to write the posterior distribution */
	FILE* pfile = fopen("posterior_measure.txt", "w");
	int n = (int) pow(2, n2);
	int mcmc = (int) pow(2, mcmc2);
	/* Residual error produced by the bayesian inversion */
	double err = 0;
	int i, j;
	
	/* Estimated parameters */
	double* map = malloc(sizeof(double)*expansion);
	/* Covariance matrix for the gaussian */
	double* cov = malloc(sizeof(double)*expansion*expansion);
	/* Starting point where to start the chain */
	double* start = malloc(sizeof(double)*expansion);
	if(map == NULL || cov == NULL || start == NULL){
		fprintf(stderr, "malloc failed\n");
		return -1;
	}

	/* Reset map, set a random starting point, a small covariance matrix */
	for(i=0; i<expansion; ++i){
		map[i] = 0;
		start[i] = rndmUniformIn(-10., 10.);
		for(j=0; j<expansion; ++j){
			cov[i + j*expansion] = (i==j)? 0.9: 0.1;
		}
	}

	printf("** Starting point:\n");
	printVec(start, expansion);
	printf("\n%d samples, %d iterations per sample\n", n, mcmc);


	/* Proceed with the bayesian inversion:
	 * n = number of posterior samples to produces;
	 * mcmc = number of monte carlo iterations;
	 * map = the vector which will contain the most frequent sample = solution = MAP
	 * NULL = the value of the true solution, a,b, is not known
	 * solver = the linear interpolation defined above
	 * observed = vector containing the y_i
	 * expansion = domain's dimension
	 * observed = codomain's dimension
	 * noise during the mcmc chain = mcmc_noise
	 * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
	 * cov = my covariance matrix, prior gaussian
	 * start = starting point for the chain
	 * pfile = file will contain the posterior distribution (values, probabilities)
	 * 0 = no verbose/debug mode */
	err=bayInv(n, mcmc, map, true_params, solver, observed,
		expansion, num_observations, mcmc_noise, 0.2, cov, start, pfile, 0);

	/* err contains the residual error */
	/* Print the results */
	printf("MAP: ");
	printVec(map, expansion);
	printf("RES ERR: %.3f%%\n", err);
	printf("Observed output:\n");
	printVec(observed, num_observations);
	printf("MAP output :\n");
	solver(map, expansion, observed, num_observations);
	printVec(observed, num_observations);

	/* Free all the allocated memory */
	free(true_params);
	free(observed);
	free(map);
	free(cov);
	free(start);
	fclose(pfile);
	return 0;
}
