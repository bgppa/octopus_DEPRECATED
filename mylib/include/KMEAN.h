/* KMEAN.h header file */
/* AGAIN, everything need to be properly checked, commented
 * and tested!!!! */
#ifndef _KMEAN_H_
#define _KMEAN_H_

#if 0 /* Obsolete functions here no more needed */
double distance(double* v1, double* v2, int n);
int areEqual(double* v1, double* v2, int n, double tol);
int readPoints(char* file_name, double* list_of_points, int dim_each, int lines);
#endif


void findKthCentroid(double* list_of_points, int dim_each, int lines, double* centroid, int* labels, int K);

double kMean(double* dati, int l, int r, int N, FILE* file_output, int MAX_ITER, double* MAP);

#endif
