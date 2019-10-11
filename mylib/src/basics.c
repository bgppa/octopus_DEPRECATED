/* Code for implementing some basics vector/matrix operations.
 * The use of blas would have been more efficient.
 * The author is aware of that, he preferred to write them "again"
 * as a warming-up exercise, and to guarantee a complete 100% control */

/* General rules to keep in mind:
 1) usually, the comments are available in the headers, too;
 2) every function **assume** to receive an already allocated pointer;
 3) NULL pointers, negative/zero dimensions are not allowed
    and brutally handled via assert() */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "basics.h"

/* Print the content of a n times m matrix A */
void printMat(const double *A, int n, int m)
{
        assert(A != NULL);
        assert(n > 0);
        assert(m > 0);          

        int tot_dim = n*m;
        for(int i = 0; i < tot_dim; ++i){
                printf("%.3e%c", A[i], ((i + 1) % m) ? ' ' : '\n');
        }
}


/* Print the content of a n times m matrix A to file F */
void fprintMat(FILE *F, const double *A, int n, int m)
{
        assert(F != NULL);
        assert(A != NULL);
        assert(n > 0);
        assert(m > 0);          

        int tot_dim = n * m;
        for(int i = 0; i < tot_dim; ++i){
                fprintf (F, "%.3e%c", A[i], ((i + 1) % m) ? ' ' : '\n');
        }
}


/* Print the content of a d dimensional array */
void printVec(const double *v, int d)
{
        assert(v != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                printf("%.3e ", v[i]);
        printf("\n");
}

/* Print to a file the contento of a d dimensional array */
void fprintVec(FILE *F, const double *v, int d)
{
        assert(F != NULL);
        assert(v != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                fprintf(F, "%.3e ", v[i]);
        fprintf(F, "\n");
}


/* y = alpha*x + y, where y,x vectors of dimension d, alpha scalar */
void axpy(const double *x, double *y, int d, double alpha)
{
        assert(x != NULL);
        assert(y != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                y[i] += alpha * x[i];
}


/* y = alpha*x - y, where x,y vectors, alpha scalar */
void axmy(const double *x, double *y, int d, double alpha)
{
        assert (x != NULL);
        assert (y != NULL);
        assert (d > 0);

        for(int i = 0; i < d; ++i)
                y[i] -= alpha * x[i];
}


/* Dot product between two vectors of dimension d */
double dot(const double *x, const double *y, int d)
{
        assert(x != NULL);
        assert(y != NULL);
        assert(d > 0);

        double res = 0.;
        for(int i = 0; i < d; ++i)
                res += x[i] * y[i];
        return res;
}

/* l2 norm of a vector of dimension d */
double nrm2(const double *x, int d)
{
        assert (x != NULL);
        assert (d > 0);

        return sqrt(dot(x,x,d));
}

/* l2 distance between two vectors */
double nrm2dist(const double *v1, const double *v2, int d)
{
        assert(v1 != NULL);
        assert(v2 != NULL);
        assert(d > 0);
        
        double sum = 0;
        int i = 0;
        for(i = 0; i < d; ++i){
                sum += pow(v1[i] - v2[i], 2.);
        }
        return sqrt(sum);
}

/* l1 norm of a vector of dimension d */
double nrm1(const double *x, int d)
{
        assert(x != NULL);
        assert(d > 0);

        double sum = 0.;
        for(int i = 0; i < d; ++i)
                sum += fabs(x[i]);
        return sum;
}

/* sup norm of a vector of dimension d */
double nrmsup(const double *x, int d)
{
        assert(x != NULL);
        assert(d > 0);

        double sup = 0;
        for(int i = 0; i < d; ++i){
                if(fabs(x[i]) > sup)
                        sup = fabs(x[i]);
        }
        return sup;
}


/* Copy the vector in x into y; d is the dimension */
void copy(const double *x, double *y, int d)
{
        assert(x != NULL);
        assert(y != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                y[i] = x[i];
}


/* x = alpha*x, where alpha a scalar and x a vector. d the dimension */
void scal(double *x, int d, double alpha)
{
        assert(x != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                x[i] *= alpha;
}


/* y = beta*y + alpha*A*x
 * where A is a matrix of dimension dimy x dimx
 * x is a vector of dimension dimx
 * y is a vector of dimension dimy
 * alpha and beta are scalars */
void gemv(const double *A, double *y, const double *x, double alpha,
                double beta, int dimy, int dimx)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(x != NULL);
        assert(dimy > 0);
        assert(dimx > 0);

        double tmp = 0;
        int j;
        for(int i = 0; i < dimy; ++i, tmp = 0){
                for(j = 0; j < dimx; ++j) 
                        tmp += x[j] * A[dimx * i + j];
                y[i] += beta * y[i] + alpha * tmp;
        }
}


/* For a square matrix A of dimension d x d,
 * computes the product Ax and stores the result in the vector y */
void axiny(const double *A, const double *x, double *y, int d)
{
        assert(A != NULL);
        assert(x != NULL);
        assert(y != NULL);
        assert(d > 0);

        for(int i = 0, j = 0; i < d; ++i){
                y[i] = 0.;
                for(j = 0; j < d; ++j)
                        y[i] += x[j] * A[d*i+j];
        }
}


/* Define A as the identity matrix of dimension n */
void id(double *A, int n)
{
        assert(A != NULL);
        assert(n > 0);

        int j;
        for(int i = 0; i < n; ++i)
                for(j = 0; j < n; ++j){
                        if(i == j)
                                A[i * n + j] = 1.;
                        else
                                A[i * n + j] = 0.;
                }
}


/* Creates the i-basis vector of dimension d.
 * Eg: d = 3, i = 2, then e_2 = (0,1,0) */
void ei(double *v, int d, int i)
{
        assert(v != NULL);
        assert(d > 0);
        assert(i > 0);

        for(int k = 0; k < d; ++k)
                if(k == i-1)
                        v[k] = 1.;
                else
                        v[k] = 0.;
}


/* A is n x m, B is m x n,
 * put in B the transport of A */
void transp(const double *A, double *B, int n, int m)
{
        assert(A != NULL);
        assert(B != NULL);
        assert(n > 0);
        assert(m > 0);

        int i,j;
        for(i = 0; i < m; ++i){
                for(j = 0; j < n; ++j){
                        B[n * i + j] = A[m * j + i];
                }
        }
}


/* Copy in v the i-th column of a n x m matrix,
 * i going from 0 to m-1 (incl) */
void column(const double *A, double *v, int n, int m, int i)
{
        assert(A != NULL);
        assert(v != NULL);
        assert(n > 0);
        assert(m > 0);
        assert(i >= 0 && i <= m-1);

        for(int  j = 0; j < n; ++j){
                v[j] = A[i + m * j];
        }
}


/* Copy in v the i-th row of a n x m matrix,
 * i going from 0 to (n-1) incl. */
void row(const double *A, double *v, int n, int m, int i)
{
        assert(A != NULL);
        assert(v != 0);
        assert(n > 0);
        assert(m > 0);
        assert(i >= 0 && i <= n-1);

        for(int j = 0; j < m; ++j){
                v[j] = A[i * m + j];
        }
}

/* Assuming matrices:
 * A, n x k
 * B, k x m
 * stores the product AxB into C, matrix n x m */
void matmul(const double *A, const double *B, double *C, int n, int k, int m)
{
        assert(A != NULL);
        assert(B != NULL);
        assert(C != NULL);
        assert(n > 0);
        assert(m > 0);
        assert(k > 0);
        assert(m > 0);

        int i,j;
        double *ai = malloc(sizeof(double) * k);
        double *bj = malloc(sizeof(double) * k);
        for(i = 0; i < n; ++i){
                for(j = 0; j < m; ++j){
                        row(A, ai, n, k, i);
                        column(B, bj, k, m, j);
                        C[i * m + j] = dot(ai, bj, k);
                }
        }
        free(ai);
        free(bj);
}


/* A, B, square matrices of dimension n: store AxB into B */
void matmulsq(const double *A, double *B, int n)
{
        assert(A != NULL);
        assert(B != NULL);
        assert(n > 0);

        double *C = malloc(sizeof(double) * n * n);
        matmul(A, B, C, n, n, n);
        copy(C, B, n * n);
        free(C);
}

/* Put in w the difference v-w */
void diff(const double *v, double *w, int d)
{
        assert(v != NULL);
        assert(w != NULL);
        assert(d > 0);

        for(int i = 0; i < d; ++i)
                w[i] -= v[i];
}

/* Check if two vector are equal, giving a tolerance explicitly
 * you can give a matrix as a parameter, assuming d = n x m */
int isequaltol(const double *v, const double *w, int d, double tol)
{
        assert(v != NULL);
        assert(w != NULL);
        assert(d > 0);
        assert(tol > 0);

        double sum = 0.;
        for(int i = 0; i < d; ++i){
                sum += (v[i] - w[i]) * (v[i] - w[i]);
        }
        if(sqrt(sum) < tol)
                return 1;
        else
                return 0;
}

/* Check if two vector are equal by using the default EPS value
 * You can give a matrix a parameter, assuming d = n x m */
int isequal(const double *v, const double *w, int d)
{
        #ifdef EPS
                return(isequaltol(v, w, d, EPS));
        #else
                printf("EPS not defined!\n");
                return 0;
        #endif
}

