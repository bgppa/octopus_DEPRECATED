/* kmeans.h header file */
/* AGAIN, everything need to be properly checked, commented
 * and tested!!!! */
#ifndef _KMEANS_H_
#define _KMEANS_H_

void kMeans(double *data, int len, int dim, int cent_num, int max_iter,
            double *done, int verbose);

#endif
