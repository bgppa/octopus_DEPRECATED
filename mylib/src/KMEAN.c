/* THIS LIBRARY HAS STILL TO BE TESTED AND COMMENTED */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"KMEAN.h"
#include "BASIC.h"
#define DEBUG 0
 
/* Replace DEBUG with verbose? */

/* Known problems:
 * remember the kMean algorithm is an heuristic attempt;
 * it fails if two centroids are computed to be equal,
 * anyway this status is NOT YET checked into the following
 * source code.
 * For the moment it is highly experimental and must be
 * properly checked, as well as commented, again */

#if 0
double distance(double* v1, double* v2, int n){
	double s=0;
	for(int i=0; i<n; ++i){
		s += pow((v1[i] - v2[i]), 2.);
	}
	return sqrt(s);
}

int areEqual(double* v1, double* v2, int n, double tol){
	if(distance(v1, v2, n) < tol){
		return 1;
	}
	else{
		return 0;
	}
}

/* replace "distance" with nrm2dist, "areEqual" with "isequal" */
#endif 

#if 0
void findCentroid(double* list_of_points, int dim_each, int lines, double* centroid){
	int i, j;
	for(i=0; i<dim_each; ++i){
		centroid[i] = 0;
		for(j=0; j<lines; ++j){
			centroid[i] += list_of_points[j*dim_each + i];
		}
		centroid[i] /= (double) lines;
	}
}
#endif

void findKthCentroid(double* list_of_points, int dim_each, int lines, double* centroid, int* labels, int K){
	int i, j;
	int actual_belonging = 0;
	for(i=0; i<dim_each; ++i){
		centroid[i] = 0;
		actual_belonging = 0;
		for(j=0; j<lines; ++j){
			if(labels[j] == K){
				++actual_belonging;
				centroid[i] += list_of_points[j*dim_each + i];
			}
		}
		centroid[i] /= (double) actual_belonging;
	}
}

/* This function is actually not used here - can be recycled to read observation */

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
			
#if 0
void printMat(double* A, int n, int m){
	int i;
	int j;
	for(i=0; i<n; ++i){
		for(j=0; j<m; ++j){
			printf("%f ", A[i*m+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printVec(double* v, int n){
	for(int i=0; i<n; ++i){
		printf("%f ", v[i]);
	}
	printf("\n");
}

void printIntVec(int* v, int n){
	for(int i=0; i<n; ++i){
		printf("%d ", v[i]);
	}
	printf("\n");
}
#endif


double kMean(double* dati, int l, int r, int N, FILE* file_output, int MAX_ITER, double* MAP){

	int i, j;
	int* labels = malloc(sizeof(int) * l);

	double* prev_centroids = malloc(sizeof(double)*r * N); /* List of N centroids */
	double* next_centroids = malloc(sizeof(double)*r * N);

	/* At each centroid is assigned its probability according to the frequencies */
	double* frequencies = malloc(sizeof(double) * N);

	/* La convergenza la implemento dopo usando la norma, etc...*/
	/* Step 0: initialize as centroid the first N available points */
	for(i=0; i<N; ++i){
		for(j=0; j<r; ++j){
			prev_centroids[i*r + j] = dati[i*r + j];
			next_centroids[i*r + j] = dati[i*r + j];
		}
	}

#if DEBUG
	printf("Ok, let's start with the following centroids\n");
	printMat(prev_centroids, N, r);
#endif

	/* Step 1: set the labels: for each point i-th, the i-ht label
	 * corresponds to the closest centroid */
	double dist=0;
	double tmp=0;

for(int iteration=0; iteration<MAX_ITER; ++iteration){
	dist=0;
	tmp=0;
	/* For each point... */
	for(i=0; i<l; ++i){
#if DEBUG
		printf("Finding the label for point %d\n", i);
#endif
		/* ...find the centroid where it belongs: */
		/* Set the label and distance from the centroid number 0 */
		labels[i] = 0;
		dist = nrm2dist(dati+(i*r), prev_centroids, r);
#if DEBUG
		printf("Distance from centroid 0: %f\n", dist);
#endif
		/* COmpare the distance with all the remaining centroids,
		 * then save the smallest */
		for(j=1; j<N; ++j){
			tmp = nrm2dist(dati+(i*r), prev_centroids+(j*r), r);
#if DEBUG
			printf("Distance from centroid %d: %f\n", j, tmp);
#endif
			if(tmp < dist){
				dist = tmp;
				labels[i] = j;
			}
		}
#if DEBUG
		printf("Chosen label: %d\n\n", labels[i]);
#endif
	}

#if DEBUG
	printf("Current labels:");
	printIntVec(labels, l);
#endif

	/*Ok, now the labels should be ready. The labels describe the clusters.
	 * from clusters comes the subdivision */

	for(i=0; i<N; ++i){
		findKthCentroid(dati, r, l, next_centroids+(i*r), labels, i);
	}

#if DEBUG
	printf("Iteration number %d\n", iteration);
	printf("Centroid from:\n");
	printMat(prev_centroids, N, r);
	printf("To:\n");
	printMat(next_centroids, N, r);
#endif
	/* It can happen to have enough number of iterations */
	if( isequal( prev_centroids, next_centroids, N*r ) ){
		break;
	}
	
	/*Copy next_centroids to prev_centroids */
	for(int h=0; h<N*r; ++h){
		prev_centroids[h] = next_centroids[h];
	}
#if DEBUG
	getchar();
#endif

}

#if DEBUG
printf("Centroids:\n");
printMat(next_centroids, N, r);
//getchar();

printf("Labels:\n");
printIntVec(labels, l);
//getchar();

#endif
	for(i=0; i<N; ++i){
		frequencies[i] = 0;
	}

	/* Let's compute the frequency of every centroid and save the highest */
	double highest_freq = 0;
	int its_index = 0;
	for(i=0; i<N; ++i){
		/* Computing the frequency of centroid number i */
		/* Sum one for each corresponding label */
		for(j=0; j<l; ++j){
			if(labels[j] == i){
				++frequencies[i];
			}
		}
		/* Divide for the total number and normalize to 100 */
		frequencies[i] = (frequencies[i] / (double) l) * 100.;
		/* Check if it's the most frequent */
		if(frequencies[i] > highest_freq){
			highest_freq = frequencies[i];
			its_index = i;
		}
	}

#ifdef DEBUG
printf("OK, frequenze salvate\n");
printVec(frequencies, N);

printf("Highest frequence: centroid number %d, frequency:  %f\n", its_index, highest_freq);
#endif

	/* Save the highest centroid into MAP */
	for(i=0; i<r; ++i){
		MAP[i] = next_centroids[r*its_index + i];
	}

#ifdef DEBUG
printf("MAP:");
printVec(MAP, r);
//getchar();
#endif
	
	/* Print that will be moved to file */
	/* if file_name != NULL... */
	for(i=0; i<N; ++i){
		for(j=0; j<r; ++j){
			fprintf(file_output, "%f ", next_centroids[i*r+j]);
		}
		fprintf(file_output, "%.2f\n", frequencies[i]);
	}
	printf("\n");

	double sss=0;
	for(i=0; i<N; ++i){
		sss+=frequencies[i];
	}
/*
	printf("Total frequencies sum: %.6f\n", sss);

	printf("MAP:\n");
	printVec(MAP, r);
	printf("\n");
*/
	free(labels);
	free(prev_centroids);
	free(next_centroids);
	free(frequencies);
	return highest_freq;
}

#if 0
/* Example of usage into a main function */

int main(int argc, char*argv[]){
	int l = atoi(argv[1]);
	int r = atoi(argv[2]);
	double* dati = malloc(sizeof(double) * l * r);
	double* MAP = malloc(sizeof(double) * r);
	int N = (int) sqrt(l);
	FILE* posterior = fopen("posterior.txt", "w");

	printf("Reading data...\n");
	readPoints("dati.txt", dati, r, l);
	double f = kMean(dati, l, r, N, posterior, 1000, MAP);

	printf("MAP:\n");
	printVec(MAP, r);
	printf("(with frequency %.3f%%)\n", f);

	return 0;
}

#endif
