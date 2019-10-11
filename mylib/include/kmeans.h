/* kmeans.h header file */
/* AGAIN, everything need to be properly checked, commented
 * and tested!!!! */
#ifndef _KMEANS_H_
#define _KMEANS_H_

void findKthCentroid(double *list_of_points, int dim_each,
                int lines, double *centroid, int *labels, int K);

double kMeans(double *dati, int l, int r, int N,
                 FILE *file_output, int MAX_ITER, double *MAP);

#endif
