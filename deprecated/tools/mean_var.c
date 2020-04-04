/* Simple script for estimating mean and variance of real numbers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double mean (double *v, int dim)
{
        assert(v != NULL && dim > 0);
        double sum = 0;
        for (int i = 0; i < dim; ++i) {
                sum += v[i];
        }
        return sum / dim;
}

double var (double *v, int dim, double m)
{
        /* Given data, their lenghts and their mean,
         * compute the variance */
        double sum = 0;
        for (int i = 0; i < dim; ++i) {
                sum += pow((v[i] - m), 2.);
        }
        sum /= (dim - 1.);
        return sqrt(sum);
}

int main() {
        int maxlen = 5000;
        int actual = 0;
        double *vect = malloc(sizeof(double) * maxlen);
        assert(vect != NULL);
        double *tmp = NULL;
        while (scanf("%lf", vect+actual) != EOF && actual < maxlen) {
                printf("Read: %f [%d]\n", vect[actual], actual);
                getchar();
                ++actual;
                if (actual == maxlen) {
                        tmp  = realloc (vect, maxlen * 2);
                        if (tmp != NULL){
                                vect = tmp;
                                tmp = NULL;
                                maxlen *= 2;
                                printf("Realloc to %d...\n", maxlen);
                                getchar();
                        } else {
                                printf("Out of memory...\n");
                        }
                }
        }
        printf("Read %d data\n", actual);
        getchar();
        double mm = mean(vect, actual);
        printf("MEAN: %.3f\n", mm);
        printf("QVAR: %.3f\n", var(vect, actual, mm));
        free(vect);
        return 0;
}
