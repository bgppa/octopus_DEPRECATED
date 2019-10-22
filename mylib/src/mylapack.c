/* Elementary Linear Algebra opeations.
 * This library should in principle be replaced by LAPACK */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "mylapack.h"
#include "myblas.h"


/* Interpret A as an upper triangular matrix;
 * solve Ax = y and put the result into y */
void uptri(const double *A, double *y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        double sum = 0;
        int j;
        for (int k = n - 1; k >= 0; --k) {
                /* Compute the back-sostitution partial sum */
                for (j = k + 1, sum = 0; j <= n - 1; ++j) {
                        sum += y[j] * A[k * n + j];
                }
                y[k] = (y[k] - sum) / A[k * n + k];
        }
}


/* Interpret A as a lower triangular matrix;
 * solve Ax = y and put the result into y */
void lwtri(const double *A, double *y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        int j;
        double sum = 0;
        for (int k = 0; k < n; ++k) {
                /* Compute the forward-sostitution partial sum */
                for (j = 0, sum = 0; j < k; ++j) {
                        sum += y[j] * A[k * n + j];
                }
                y[k] = (y[k] - sum) / A[k * n + k];
        }
}


/* Given the linear system Ax = y, find its upper tringular
 * equivalent formulation. So the matrices A and y are **modified**
 * during the process.
 * To keep in mind: since the algorithm is via the Gauss elimination,
 * the determinant of A is invariant */
void touptri(double *A, double *y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        int m, k, j;
        double pvt;
        /* Multiply lines to achieve an upper triangular matrix */
        for (m = 0; m < n - 1; ++m) {
                for (j = m + 1; j < n; ++j) {
                        if (A[m * n + m] == 0){
                                printf("ZERO ON THE DIAGONAL. %d"
                                  "WILL FAIL! - feature to be fixed!\n", m+1);
                        }
                        pvt = A[j * n + m];
                        for (k = m; k < n; ++k) {
                            A[n * j + k] -= A[m * n + k] * (pvt / A[m * n + m]);
                        }
                        y[j] -= pvt * y[m] / A[m * n + m];
                }
        }
}


/* Solve Ax = y via gauss elimination */
int solgse(double *A, double *y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        printf("DEBUG: inverting the matrix: \n");
        printMat(A, n, n);
        touptri(A, y, n);
        printf("DEBUG: now the matrix is upper tridiagonal!");
        printMat(A, n, n);
        uptri(A, y, n);
        return 1;
}

/* Compute the determinant of a n x n matrix A,
 * using essentially the gauss elimination */
double det(double *A, int n)
{
        assert(A != NULL);
        assert(n > 0);

        int m, k, j;
        double pvt;
        /* Repeat similar steps as in the gauss elimination... */
        for (m = 0; m < n - 1; ++m) {
                for (j = m + 1; j < n; ++j) {
                        pvt = A[j * n + m];
                        for (k = m; k < n; ++k) {
                            A[n * j + k] -= A[m * n + k] * (pvt / A[m * n + m]);
                        }
                }
        }

        /* ...to obtain a matrix whose det is just the diagonal product */
        for (m = 0, pvt = 1.0; m < n; ++m) {
                pvt *= A[n * m + m]; 
        }
        return pvt;
}


/* Write in B the inverse matrix of A, dimension n x n.
 * Done via gauss elimination on Ax = e_i */
void invgse(const double *A, double *B, int n)
{
        assert(A != NULL);
        assert(B != NULL);
        assert(n > 0);

        int q;
        double *copyofA = malloc(sizeof(double) * n * n);
        /* Strategy: assign to B the solution of Ax = e_i */
        double *e_tmp = malloc(sizeof(double) * n);
        for (int i = 0; i < n; ++i) {
                copy(A, copyofA, n * n);
                /* e_tmp is now the i-th basis of the n-dim space
                 * Remember that the indeces are switched by one */
                ei(e_tmp, n, i + 1);
                solgse(copyofA, e_tmp, n);
                /* ^---< here a KEY STEP. Can indeed use a faster method.
                 * The result is now stored in e_tmp 
                 * so, e_tmp is the i-th column of the inverse matrix:
                 * let's copy it in B */
                for (q = 0; q < n; ++q) {
                        B[n * q + i] = e_tmp[q];
                }
        }
        free(e_tmp);
        free(copyofA);
}


/* Solve a linear system Ax=b (dimension d) by using the SOR iteration with
 * parameter omega, tolerance eps and maximum number of iteration nmax 
 * x is the starting point of the iteration;
 * the result is stored in b itself (so ***b is overwritten***!).
 * MEMENTO: key hypothesis for SOR: positive definiteness.
 * omega is the relaxation parameters;
 * eps the error bound requested,
 * and nmax the maximum number of allowed iterations */
int sor(const double *A, double *b, double *x, int d,
        double omega, double eps, int nmax)
{
        assert(A != NULL);
        assert(b != NULL);
        assert(x != NULL);
        assert(d > 0);
        assert(omega > 0 && omega < 2);
        assert(eps > 0);
        assert(nmax >= 0);
        
        double *tmp = malloc(sizeof(double) * d);
        double err = eps + 1.0; /* Set the error "large" */
        int ndone = 0;
        int i, j;
        double sum = 0;

        /* Repeat until reaching a small enough error */
        while (err > eps && ndone < nmax) {
                for (i = 0; i < d; ++i) {
                        for (j = 0, sum = 0; j < d; ++j) {
                                sum += A[i * d + j] * x[j];
                        }
                        x[i] -= (omega / A[i * d + i]) * (sum - b[i]);
                }
                /* Computes the error of the operation Ax - b */
                axiny(A, x, tmp, d);
                diff(b, tmp, d);
                err = nrm2(tmp, d);
                ++ndone;
        }
        
        copy(x, b, d);
        free(tmp);
        printf("Final SOR error: %.3f\n", err);
        return ndone; /* Return the number of iteration done */
}


/* Solve a linear equation Ax=b of dimension d by using the
 * Jacobi iterative method. x_next is the starting value,
 - eps the error bound requested,
 - nmax the maximum number of iterations.
 * The result is stored in b itself, so **b is overwritten** 
 * MEMENTO: key hypothesys for Jacobi: positive definiteness. */
int jacobi(const double *A, double *b, double *x_next,
           int d, double eps, int nmax)
{
        assert(A != NULL);
        assert(b != NULL);
        assert(x_next != NULL);
        assert(d > 0);
        assert(eps > 0);
        assert(nmax >= 0);

        double *x_prev = malloc(sizeof(double) * d);
        double *tmp = malloc(sizeof(double) * d);
        int i, j;
        double sum;
        double err = eps + 1.0;
        int ndone = 0;

        /* set the starting point to x_next */  
        copy(x_next, x_prev, d);
        while (err > eps && ndone < nmax) {
                for (i = 0; i < d; ++i) {
                        /*Compute the right-hand-side of the Jacobi algorithm*/
                        for (j = 0, sum = 0; j< d ; ++j) {
                                if (j != i) {
                                        sum += A[i * d + j] * x_prev[j];
                                }
                        }
                        x_next[i] = (-sum + b[i]) / A[i * d + i];
                        copy(x_next, x_prev, d);
                }

                /* x_next is now the candidate.
                 * Compute the error of Ax_next - b */
                axiny(A, x_next, tmp, d);
                diff(b, tmp, d);
                err = nrm2(tmp, d);
                ++ndone;
        }
        copy(x_next, b, d);
        free(x_prev);
        free(tmp);
        printf("Final jacobi error: %.3f\n", err);
        return ndone; /* Return the number of iteration done */
}

#ifdef TOL
#ifdef NMAX
#ifdef OMG
/* Solve Ax=y via a default-mode for the SOR iterator.
 * It's a simplified wrapping of the functions above.
 * omega = OMG, a #define constant specified in the header,
 * the starting value is 0 and the tolerance TOL, again
 * specified in the header file (in a way to be easily customized
 * by the user). Similarly, max_iteration = NMAX */
int solsor(double *A, double *y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        double *x = malloc(sizeof(double) * n);
        for (int i = 0; i < n; ++i) {
                x[i] = 0;
        }
        int val = sor(A, y, x, n, OMG, TOL, NMAX);
        free(x);
        return val;
}
#endif /* Check for OMG completed */

/* Solve Ax=y via a default-mode for the JACOBI iterator.
 * Starting value = 0, tolerance=TOL,
 * max_iteration = NMAX. See the description of solsor above
 * for further information on the constants */
int soljcb(double* A, double* y, int n)
{
        assert(A != NULL);
        assert(y != NULL);
        assert(n > 0);

        double* x = malloc(sizeof(double) * n);
        for (int i = 0; i < n; ++i) {
                x[i] = 0.;
        }
        int val = jacobi(A, y, x, n, TOL, NMAX);
        free(x);
        return val;
}

#endif /* If TOL and NMAX are available, soljdb is fine */
#endif 

/* Compute the inverse matrix of A and store the result in B.
 * The method solves Ax=e_i for each i.
 * The hser must specify the derised way to solve such a system,
 * solgse (Gauss), solsor (SOR), soljcb (JACOBI) */
void invmat(const double *A, double *B, int n,
            int (*solver) (double*, double*, int))
{
        assert(A != NULL);
        assert(B != NULL);
        assert(n > 0);
        assert(solver != NULL);

        int q;
        double *copyofA = malloc(sizeof(double) * n * n);
        /* Strategy: assign to B the solution of Ax = e_i */
        double *e_tmp = malloc(sizeof(double) * n);
        for (int i = 0; i < n; ++i) {
                copy(A, copyofA, n * n);
                /* e_tmp is now the i-th basis of the n-dim space.
                 * Remember that the indeces are switched by one */
                ei(e_tmp, n, i + 1);
                solver(copyofA, e_tmp, n);
                /* The result is now stored in e_tmp;
                 * so, e_tmp is the i-th column of the inverse matrix:
                 * let's copy it in B */
                for (q = 0; q < n; ++q) {
                        B[n * q + i] = e_tmp[q];
                }
        }
        free(e_tmp);
        free(copyofA);
}

/* LINear SYMmetrizer.
 * Think about the linear system Ax = b:
 * it can happen that A has some zeros on the diagonal,
 * making then impossible to use Gaussian elimination correctly.
 * Note that the system is equivalent to A^T A x = A^T b
 * The new matrix A^T A is indeed symmetric and with strictly positive
 * diagonal. Therefore can be easier to solve, maybe using a more efficient
 * iteration method instead of Gauss elimination.
 * This functions takes A, b and store A^T A in C, while A^T b in d.
 * Then you can use the data in C and d to solve the system equivalently.
 * the interger d is the dimension of b, which is also the matrix's dimension
 * (A is assumed to be square of dimension d * d).
 * The purpose of this function is prepare the system in a way not to suffer
 * from the zero-on-the diagonal problem which makes
 * gaussian elimination fails. */
void linsym(const double *A, const double *b, double *C, double *d, int dim)
{
        assert(A != NULL);
        assert(b != NULL);
        assert(C != NULL);
        assert(d != NULL);
        assert(dim > 0);

        double *AT = malloc(sizeof(double) * dim * dim);
        transp(A, AT, dim, dim);
        matmul(AT, A, C, dim, dim, dim);
        axiny(AT, b, d, dim);
        free(AT);
}
