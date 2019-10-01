#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>
#include"BASIC.h"
#include"RANVAR.h"
#include"PCNINV.h"
#include"FILEIO.h"

double* glb_eval_points;
int glb_d2;
int glb_d1;

/* eval_points will be an array of domension "codomain",
 * containing the x w.r.t the linear evaluation ax+b
 * will be done. I mean: the input file for this script
 * is a 2-column set of data. The first corresponding
 * to the x, the second to the y. I want to find the proper
 * parameter to fix f(x) = y, so the x have to be registered
 * in order to produce the G operator defined as the 
 * evaluation: a,b -> [a x_1 + b, a x_2 + b, ... , a x_codomain + b]
 * The parameter d1 in solver is fixed to be 2, but left because
 * of structural compatibility with the bayesin algorithm.
 * the codomain is determined by the number of lines read.
 * This value will be stored into glb_d2, the passed as argument (see main).
 * ab is a vector whose first value represents a, the second b.
 * The function is named solver just for tradition, coherence with the other
 * examples provided in the library */
void solver(const double* ab, int d1, double*y, int codomain){
	/* The variable d1 is actually not used, but left for
	 * pointer compatibility with bayInv */
	(void) d1;

	int i=0;
	for(i=0; i<codomain; ++i){
		y[i] = ab[0] * glb_eval_points[i] + ab[1];
	}
}


int main(int argc, char* argv[]){
	/* Takes a valid filename as argument */
	if(argc != 2){
		printf("Error: must specify a single argument: file dataset\n");
		return -1;
	}

	/* read_from_file is an array that will contain *all* the numbers
	 * read from the file (sequentially). It is firstly allocated
 	 * for containing one number (allocating now = do not forget to free
	 * later). By passing it to dataFromFile it will be constantly
	 * reallocated until becoming of the size equal to the total
	 * number of points. This value is stored into total_points.
	 * Therefore: read_from_file will have size total_points. */
	double* read_from_file = malloc(sizeof(double));
	int total_points = dataFromFile(argv[1], &read_from_file);
	if(total_points < 2){
		printf("error: too few points or reallocating problems\n");
		return -1;
	}

	/* Since *by assumption* this is a script for 1-dimensional
	 * linear regression, the number of columns is assumed to be 2
	 * (i.e. the file is seen as a series of couples x_i, y_i)
	 * and the number of lines computes accordingly */
	int columns = 2;
	int lines = total_points / columns;
	if(lines * columns != total_points){
		printf("Odd number of points - invalid dataset\n");
		return -1;
	}
	printf("%d lines and %d columns\n", lines, columns);
	printVec(read_from_file, total_points);

	/* Store now the columns into x_i and y_i, i.e.
	 * places of evaluations - the x_i -
	 * and observed outputs - the y_i */
	int i,j;
	double* eval_points = malloc(sizeof(double)*lines);
	double* observed = malloc(sizeof(double)*lines);
	if(eval_points == NULL || observed == NULL){
		fprintf(stderr, "malloc failed\n");
		return -1;
	}
	/* Set the first column as x_i, the second ad y_i */
	for(i=0; i<lines; ++i){
		eval_points[i] = read_from_file[i*columns];
		observed[i] = read_from_file[i*columns + 1];
	}

	printf("x_i : ");
	printVec(eval_points, lines);
	printf("y_i : ");
	printVec(observed, lines);

	/* By passing the parameters to the following global variables,
	 * we set the function "solver" ready to be used in the bayesian
	 * inverse function */
	glb_eval_points = eval_points;
	glb_d2 = lines;
	glb_d1 = 2;

	/* Now that the data are ready, set the bayes parameters */
	srand(time(NULL));
	/* Output file where to write the posterior distribution */
	FILE* pfile = fopen("posterior_measure.txt", "w");
	int n = 500;
	int mcmc = 2000;

	/* Residual error produced by the bayesian inversion */
	double err = 0;
	
	/* Estimated parameters */
	double* map = malloc(sizeof(double)*columns);
	/* Covariance matrix for the gaussian */
	double* cov = malloc(sizeof(double)*columns*columns);
	/* Starting point where to start the chain */
	double* start = malloc(sizeof(double)*columns);
	if(map == NULL || cov == NULL || start == NULL){
		fprintf(stderr, "malloc failed\n");
		return -1;
	}

	/* Reset map, set a random starting point, a small covariance matrix */
	for(i=0; i<columns; ++i){
		map[i] = 0;
		start[i] = rndmUniformIn(-10., 10.);
		for(j=0; j<columns; ++j){
			cov[i + j*columns] = (i==j)? 0.2 : 0.1;
		}
	}

	/* Proceed with the bayesian inversion:
	 * n = number of posterior samples to produces;
	 * mcmc = number of monte carlo iterations;
	 * map = the vector which will contain the most frequent sample = solution = MAP
	 * NULL = the value of the true solution, a,b, is not known
	 * solver = the linear interpolation defined above
	 * observed = vector containing the y_i
	 * glb_d1 = 2, the domain dimension, 2 since two parameters: a, b
	 * gld_d2 = the codomain dimension, i.e. the lines in the file, number of y_i
	 * 0.15 = noise
	 * 0.2 = beta, parameter for the pCN algorithm 0<beta<1
	 * cov = my covariance matrix, prior gaussian
	 * start = starting point for the chain
	 * pfile = output file where to write the posterior distribution (values, probabilities)
	 * 0 = no verbose/debug mode */
	err=bayInv(n, mcmc, map, NULL, solver, observed, glb_d1, glb_d2, 0.15, 0.2, cov, start, pfile, 0);

	/* err contains the residual error */
	/* Print the results */
	printf("MAP: ");
	printVec(map, glb_d1);
	printf("RES ERR: %.3f%%\n", err);

	/* Free all the allocated memory */
	glb_eval_points = NULL;
	free(eval_points);
	free(observed);
	free(read_from_file);
	free(map);
	free(cov);
	free(start);
	fclose(pfile);
	return 0;
}
