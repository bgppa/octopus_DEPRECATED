/* STILL TO BE TESTED AND COMMENTED */
/* Split kMeans in smaller parts */
/* Attempt: adding the computation of entropy? */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kmeans.h"
#include "basics.h"
#define DEBUG 0
 
/* Replace DEBUG with verbose? */

/* Known problems:
 * remember the kMean algorithm is an heuristic attempt;
 * it fails if two centroids are computed to be equal,
 * anyway this status is NOT YET checked into the following
 * source code.
 * For the moment it is highly experimental and must be
 * properly checked, as well as commented, again */


/* Given a label K, find the centroid corresponding to
 * such a set of data - or something  similar; write better */

void findKthCentroid(double *list_of_points, int dim_each,
                int lines, double *centroid, int *labels, int K)
{
        int i, j;
        int actual_belonging = 0;
        for(i = 0; i < dim_each; ++i){
                centroid[i] = 0;
                actual_belonging = 0;
                for(j = 0; j < lines; ++j){
                        if(labels[j] == K){
                                ++actual_belonging;
                                centroid[i] += list_of_points[j * dim_each + i];
                        }
                }
                centroid[i] /= (double) actual_belonging;
        }
}


/* Split kMeans in smaller parts <- TO DO */

double kMeans(double *dati, int l, int r, int N,
         FILE *file_output, int MAX_ITER, double *MAP)
{
        int i, j;
        int *labels = malloc(sizeof(int) * l);

        /* List of N centroids */
        double *prev_centroids = malloc(sizeof(double) * r * N);
        double *next_centroids = malloc(sizeof(double) * r * N);

        /* At each centroid will be assigned probability
         * according to frequencies store them here */
        double *frequencies = malloc(sizeof(double) * N);

        /* The discrete entropy will be computed by summing the frequencies */
        double entropy = 0;

        /* Step 0: initialize as centroid the first N available points */
        for(i = 0; i < N; ++i){
                for(j = 0; j < r; ++j){
                        prev_centroids[i * r + j] = dati[i * r + j];
                        next_centroids[i * r + j] = dati[i * r + j];
                }
        }

#if DEBUG
        printf("Ok, let's start with the following centroids\n");
        printMat(prev_centroids, N, r);
#endif

        /* Step 1: set the labels: for each point i-th, the i-ht label
         * corresponds to the closest centroid */
        double dist = 0;
        double tmp = 0;

        for(int iteration = 0; iteration < MAX_ITER; ++iteration){
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
                        dist = nrm2dist(dati + (i * r), prev_centroids, r);
#if DEBUG
                        printf("Distance from centroid 0: %f\n", dist);
#endif
                /* Compare the distance with all the *remaining* centroids,
                 * then save the smallest (the index starts from 1). */
                        for(j = 1; j < N; ++j){
                        tmp = nrm2dist(dati+(i*r), prev_centroids+(j*r), r);
#if DEBUG
                        printf("Distance from centroid %d: %f\n", j, tmp);
#endif
                        if(tmp < dist){
                                dist = tmp;
                                labels[i] = j;
                        }
                        } /* End of for = comparison done */
#if DEBUG
                        printf("Chosen label: %d\n\n", labels[i]);
#endif
                } /* Each point has been assigned to a label */

#if DEBUG
                printf("Current labels:");
                printIntVec(labels, l);
#endif

                /* The labels are ready. Use them to create the clusters */
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
                if(isequal(prev_centroids, next_centroids, N * r)){
                        break;
                }
        
                /*Copy next_centroids to prev_centroids */
                for(int h = 0; h < N * r; ++h){
                        prev_centroids[h] = next_centroids[h];
                }
#if DEBUG
                getchar();
#endif

        } /* The centroid creation hs been repeated "MAX_ITER" times */

#if DEBUG
        printf("Centroids:\n");
        printMat(next_centroids, N, r);
        printf("Labels:\n");
        printIntVec(labels, l);
#endif
        for(i = 0; i < N; ++i){
                frequencies[i] = 0;
        }

        /* Compute the frequency of every centroid and save the highest */
        double highest_freq = 0;
        int its_index = 0;
        for(i = 0; i < N; ++i){
                /* Computing the frequency of centroid number i */
                /* Sum one for each corresponding label */
                for(j = 0; j < l; ++j){
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

        /* Compute the entropy */
        double ttt;
        for(i = 0; i < N; ++i){
                ttt = frequencies[i] / 100.;
                entropy += ttt * log(ttt);
        }
        entropy = -entropy;

        printf("\n - - ENTROPY: %.3f\n", entropy);

#if DEBUG
        printf("OK, frequenze salvate\n");
        printVec(frequencies, N);
        printf("Highest: centroid %d, freq %f\n", its_index, highest_freq);
#endif

        /* Save the highest centroid into MAP */
        for(i = 0; i < r; ++i){
                MAP[i] = next_centroids[r * its_index + i];
        }

#if DEBUG
        printf("MAP:");
        printVec(MAP, r);
#endif
        
        /* if file_name != NULL... <- TO CHECK */
        /* Write the results to a file;
         * in the bayesian setting, it is the discretized
         * posterior measure! */
        for(i = 0; i < N; ++i){
                for(j = 0; j < r; ++j){
                        fprintf(file_output, "%f ", next_centroids[i*r+j]);
                }
                fprintf(file_output, "%.2f\n", frequencies[i]);
        }
        printf("\n");

        free(labels);
        free(prev_centroids);
        free(next_centroids);
        free(frequencies);
        return highest_freq;
}
