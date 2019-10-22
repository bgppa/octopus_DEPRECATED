/* The purpose of this library is to include the kmeans algorithm.
 * Many multidimensional samples are produced by using monte carlo
 * techniques in other part of the library; here we want to reduce
 * this sample set to something smaller and more readable,
 * like the idea of producing an histogram from 1-dimensional sampled points */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "kmeans.h"
#include "myblas.h"
 
/* Known problems: in kMeans(), to add the checking that the first 
 * centroids are different points. Not-so-urgent,
 * since it is essentially impossible in the usage of
 * bayesian inverse problem */

/* Just print an array of integer - debug reasons */
void printIntVec(const int* v, int l){
        for(int i = 0; i < l; ++i){
                printf("%d ", v[i]);
        }
        printf("\n");
}

/* findKthCentroid is meant to be used ONLY by kMeans, no others.
 * So the checking are omitted since done there. Given:
 * "data" seen as a matrix with "len" rows and "col" columns;
 * an array "labels" of dimension len, seen as a marker for each data;
 * an integer k which determines at which marker we are looking to,
 * This function takes the average of all the points marked with k
 * and store it into "centroid", array of dimension "dim" */
void findKthCentroid(double *col, int len, int dim,
                     double *centroid, int *labels, int k)
{
        int i, j;
        int actual_belonging = 0;
        for (i = 0; i < dim; ++i) {
                centroid[i] = 0;
                actual_belonging = 0;
                for (j = 0; j < len; ++j) {
                        if (labels[j] == k) {
                                ++actual_belonging;
                                centroid[i] += col[j * dim + i];
                        }
                }
                centroid[i] /= (double) actual_belonging;
        }
}

/* highFreq is meant to be used ONLY by kMeans, no others.
 * So the checking are omitted since done there.
 * Comparison function for the qsort algorithm in kmeans.
 * First array considered higher according only the first element;
 * Interpretation: it will be the frequency, and we care about most freq; */ 
int highFreq(const void *a, const void *b){
        double a0 = *(const double*) a;
        double b0 = *(const double*) b; 
        return a0 < b0 ? 1 : 0;
}

/* "labels" is a sample of "n" integers each varying from 0 to "k" (not "k").
 * freq an array of dimension k. This function assigns:
 * freq[i] = %frequency of the integer i found in "labels" */
void computeFreq(double *freq, int k, const int *labels, int n){
        int i, j;
        for (i = 0; i < k; ++i) {
                freq[k] = 0;
                for (j = 0; j < n; ++j) {
                        if (labels[j] == i) {
                                ++freq[i];
                        }
                }
                /* Divide for the total number and normalize to 100 */
                freq[i] = (freq[i] / (double) n) * 100.;
        }
}



/* Given a set of numbers (1-dim samples), an hystogram store them into few
 * representative values with associated frequencies.
 * The kmeans algorithm performs the same goal in a multidimensional case.
 * The set of data is subdivided into "clusters", for each of them
 * a "centroid" (i.e. a representer) is computed.
 * The parameters are:
 - data         :       array of datas, understood to be a matrix with...
 - len          :       ...len lines...
 - dim          :       ...each containing dim values.
 - cent_num     :       number of centroid requested. Typical choice: sqrt(l)
 - max_iter     :       kmeans is an iterative-refinement algorithm.
                        This parameter set the maximum number of iterations;
 - done         :       array of dimension cent_num times dim+1.
                        It will contain the transformed data set.
                        Must have rows of lenght one more that the
                        original data, because the beginning of it
                        we'll store the centroid frequency.
 - verbose      :       when positive, produce more detailed step-by-step debug
*/ 

void kMeans(double *data, int len, int dim, int cent_num,
            int max_iter, double *done, int verbose)
{
        assert(data != NULL && len > 0 && dim > 0 && cent_num > 0);
        assert(max_iter > 0 && done != NULL);

        /* Each data point has a label to trace to which centroid is closer */
        int *labels = malloc(sizeof(int) * len);

        /* Previous and Next set of computed centroids (in each iteration) */
        double *prev_centroids = malloc(sizeof(double) * dim * cent_num);
        double *next_centroids = malloc(sizeof(double) * dim * cent_num);

        /* Each centroid will have a corresponding frequency, i.e.
         * how many data belongs to it. */
        double *frequencies = malloc(sizeof(double) * cent_num);

        assert(labels != NULL && prev_centroids != NULL);
        assert(next_centroids != NULL && frequencies != NULL);

        /* Step 0: choose as centroid the first N available points
           (then, each iteration will be a refinement)
           REMARK: ASSUMPTION: they are DIFFERENT, otherwise fail. */
        for (int i = 0; i < cent_num; ++i) {
                for (int j = 0; j < dim; ++j) {
                        prev_centroids[i * dim + j] = data[i * dim + j];
                        next_centroids[i * dim + j] = data[i * dim + j];
                }
        }

        if (verbose) {
                printf("--- kmeans() ---\n");
                printf("Starting centroids:\n");
                printMat(prev_centroids, cent_num, dim);
        }

        /* The kmeans algorithm performs varius iteration on the whole
         * set of data; each step indicates a further refinement.
         * From prev_centroids, a new list next_centroids is computed */
        for (int iteration = 0; iteration < max_iter; ++iteration) {
                double dist = 0;
                double tmp = 0;

                /* For each point... */
                for (int i = 0; i < len; ++i) {
                        if (verbose) { 
                                printf("Finding label for point %d\n", i);
                        }
                        /* ...find the centroid where it belongs: */
                        /* ...start by assigning centroid 0 to it */
                        labels[i] = 0;
                        dist = nrm2dist(data + (i * dim), prev_centroids, dim);
                        if (verbose) {
                                printf("Distance from centroid 0: %f\n", dist);
                        }
                        /* ..compare the distance with all the *remaining*
                         * and save the smallest (so the index from 1). */
                        for (int j = 1; j < cent_num; ++j) {
                                tmp = nrm2dist(data + (i * dim),
                                      prev_centroids + (j * dim), dim);
                                if (tmp < dist) {
                                        dist = tmp;
                                        labels[i] = j;
                                }
                        }
                        if (verbose) { 
                                printf("Chosen label: %d\n\n", labels[i]);
                        }
                } /* Each point has been assigned to a label */

                if (verbose) {
                        printf("Iteration %d, current labels:", iteration);
                        printIntVec(labels, len);
                }

                /* To each data corresponds now a label. Now, for each grup
                 * of data sharing the same label, find its centroid */
                for (int i = 0; i < cent_num; ++i) {
                        findKthCentroid(data, len, dim,
                                        next_centroids + (i * dim), labels, i);
                }
                if (verbose) {
                        printf("Centroid from:\n");
                        printMat(prev_centroids, cent_num, dim);
                        printf("To:\n");
                        printMat(next_centroids, cent_num, dim);
                }

                /* If there is no new centroid proposal, stop */
                if (isequal(prev_centroids, next_centroids, cent_num * dim)) {
                        break;
                }
        
                /* Otherwise copy next_centroids into prev_centroids and
                 * repeat the cycle  */
                copy(prev_centroids, next_centroids, cent_num * dim);
        } 

        if (verbose) {
                printf("Final centroids:\n");
                printMat(next_centroids, cent_num, dim);
                printf("Labels:\n");
                printIntVec(labels, len);
        }

        /* In the array frequencies, store the %frequency of each centroid.
         * (i.e. how many points are classified with its label */
        computeFreq(frequencies, cent_num, labels, len);
        if (verbose) {
                printf("%%frequencies: \n");
                printVec(frequencies, cent_num);
                printf("their sum = %f%%\n", nrm1(frequencies, cent_num));
        }

        /* Store the results into "done":
         * cent_num lines in total, each n-th line has the frequency
         * of the n-th centroid as first entry, then that centroid's value.
         * Therefore the dimension cent_num times dim + 1 */
        for (int i = 0; i < cent_num; ++i) {
                /* First entry = frequency */
                done[i * (dim + 1)] = frequencies[i];
                /* copy the i-th centroid into i-th line of done, shift by 1 */
                copy(next_centroids + i * dim, done + i * (1 + dim) + 1, dim);
        }

        /* Sort done in descreasing order of frequencies */
        qsort(done, cent_num, (dim + 1) * sizeof(double), highFreq);

        if (verbose){
                printf("Frequencies + centroids:\n");
                printMat(done, cent_num, dim + 1);
                printf("--- press a key to continue ---\n");
                getchar();
        }

        free(labels);
        free(prev_centroids);
        free(next_centroids);
        free(frequencies);
}
