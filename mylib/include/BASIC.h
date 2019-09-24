/* This is: BASIC.h */
#ifndef _BASIC_H_
#define _BASIC_H_
/* EPS is a default tolerance value used to determine when
 * two double have to be considered equals */
#define EPS 1e-8

/* Print to the screen/file a matrix A of dimension n x m, then newline */
void printMat(const double* A, int n, int m);
void fprintMat (FILE* F, const double* A, int n, int m);

/* Pront to the screen/file a vector of dimension a, followed by newline */
void printVec(const double* v, int n);
void fprintVec(FILE* file, const double* v, int d);

/* y = alpha*x + y, where y,x vectors, alpha scalar */
void axpy(const double *x, double* y, int dim, double alpha);

/* y = alpha*x - y, where x,y vectors, alpha scalar */
void axmy(const double *x, double* y, int dim, double alpha);

/* Dot product between two vectors of dimension d */
double dot(const double* x, const double *y, int d);

/* l2 norm of a vector of dimension n */
double nrm2(const double* x, int d);

/* l1 norm of a vector of dimension n */
double nrm1(const double* x, int d);

/* sup norm of a vector of dimension n */
double nrmsup(const double* x, int d);

/* Copy the vector in x into y; d is the dimension */
void copy(const double* x, double* y, int d);

/* x = alpha*x, where alpha a scalar and x a vector. d the dimension */
void scal(double* x, int d, double alpha);

/* y = beta*y + alpha*A*x
 * where A is a matrix of dimension dimy x dimx
 * x is a vector of dimension dimx
 * y is a vector of dimension dimy
 * alpha and beta are scalars */
void gemv(const double* A, double* y, const double* x, double alpha,\
		 double beta, int dimy, int dimx);

/* For a square matrix A of dimension d x d,
 * computes the product Ax and stores the result in the vector y */
void axiny(const double* A, const double* x, double* y, int d);

/* Define A as the identity matrix of dimension n */
void id(double* A, int n);

/* Define v as e_i, the i-th basis of dimension d.
 * e.g., e_2 with d=3 is (0,1,0) */
void ei(double *v, int d, int i);

/* If A is n x m, B is m x n,
 * put in B the transport matrix of A */
void transp(const double *A, double* B, int n, int m);

/* Copy in v the i-th column of a n x m matrix,
 * i going from 0 to m-1 (incl) */
void column(const double* A, double* v, int n, int m, int i);

/* Copy in v the i-th row of a n x m matrix,
 * i going from 0 to (n-1) incl. */
void row(const double* A, double* v, int n, int m, int i);

/* Assuming matrices:
 * A, n x k
 * B, k x m
 * stores the product AxB into C, matrix n x m */
void matmul(const double *A, const double* B, double* C, int n, int k, int m);

/* A, B, square matrices of dimension n: store AxB into B */
void matmulsq(const double* A, double *B, int n);

/* Put in w the difference v-w */
void diff(const double *v, double *w, int d);

/* Check if two vector are equal by using the default EPS value
 * You can give a matrix a parameter, assuming d = n x m */
int isequal(const double *v, const double *w, int d);

/* Check if two vector are equal, giving a tolerance explicitly
 * you can give a matrix as a parameter, assuming d = n x m */
int isequaltol(const double* v, const double*w, int d, double tol);

#endif /* Header guard */
