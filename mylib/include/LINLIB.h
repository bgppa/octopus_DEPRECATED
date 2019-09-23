/* LINLIB.h */
#include"BASIC.h"

#ifndef _LINLIB_H_
#define _LINLIB_H_

/* OMG: the default parameter for the sor iteration.
 * Used by solsor; */
#define OMG 1.2
/* NMAX: default maximal number of iterations for a solver */
/* Used by solsor, soljcb */
#define NMAX 10000

/* TOL: default tolerance for iterative solvers */
#define TOL 1e-8

/* A is an upper triangular matrix;
 * solve Ax=y and put the result into y */
void uptri(const double* A, double* y, int n);

/* A is a lower triangular matrix;
 * solve Ax=y and put the result into y */
void lwtri(const double* A, double* y, int n);

/* Given the linear system Ax = y, find its upper tringular
 * equivalent formulation. So the matrices A and y are **modified**
 * during the process, but since the algorithm is via the Gauss elim,
 * the determinant of A is invariant */
void touptri(double* A, double* y, int n);

/* Solve Ax=y via gauss elimination */
int solgse(double* A, double* y, int n);

/* Compute the determinant of a n x n matrix A */
/* done basically via gauss elimination */
double det(double* A, int n);

/* Put in B the inverse matrix of A, dimension n x n;
 * done via gauss elimination on Ax = e_i */
void invgse(const double* A, double* B, int n);

/* Solve a linear system Ax=b (dimension dim) by using the SOR iteration with
 * parameter omega, tolerance eps and maximum number of iteration nmax 
 * x is the starting point of the iteration;
 * the result is stored in b itself (so ***b is overwritten***!)
 * Memento: positive definiteness. */
int sor(const double* A, double* b, double* x, int dim, double omega,\
		double eps, int nmax);

/* Solve a linear equation Ax=b of dimension d by using the
 * Jacobi iterative method. x_next is the starting value, eps the tolerance
 * and nmax the maximum number of iterations.
 * The result is stored in b itself, so **b is overwritten**.
 * Memento: positive definiteness. */
int jacobi(const double* A, double* b, double* x_next, int dim, \
			 double eps, int nmax);

/* Solve Ax=y via a default-mode for the SOR iterator.
 * omega = OMG, starting value = 0, tolerance=TOL
 * max_iteration = NMAX */
int solsor(double* A, double* y, int n);

/* Solve Ax=y via a default-mode for the JACOBI iterator.
 * Starting value = 0, tolerance=TOL,
 * max_iteration = NMAX */
int soljcb(double* A, double* y, int n);

/* Compute the inverse matrix of A and store the result in B.
 * The method solves Ax=e_i for each i.
 * The user must specify the derised way to solve such a system,
 * solgse (Gauss), solsor (SOR), soljcb (JACOBI) */
void invmat(const double* A, double* B, int n, int (*solver) (double*, double*, int));

/* LINear SYMmerizer.
 * Think about the linear system Ax = b
 * it can happen that A has some zeros on the diagonal,
 * making then impossible to use Gaussian elimination correctly.
 * Note that the system is equivalent to A^T A x = A^T b
 * The new matrix A^T A is indeed symmetric and with strictly positive
 * diagonal. Therefore can be easier to solve, maybe with iterative too.
 * This functions takes A, b and store A^T A in C, while A^T b in d.
 * the interger dim is the dimension of b, which is also the matrix dim
 * (A is assumed to be square of dimension dim * dim) */ 
void linsym(const double* A, const double*b, double* C, double* d, int dim);
#endif
