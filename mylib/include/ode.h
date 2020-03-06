#ifndef _ODE_H_
#define _ODE_H_

double euler (double (*f) (double, double), double y0, double T, int N,
              double (*true_sol) (double));

double midpoint (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double));

double rkfourth (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double));

void euler_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N);

void midpoint_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N);

void rkfourth_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N);

/* Variables for the hamiltonian case */
//extern double* ham_M1;
//extern double (*ham_U) (int dim, const double *x);
//void ham_F (int dim, double t, const double *x, double *y);
double verlet (double *x, int d, double T, int N, 
		const double *M1, double(*U)(int, const double*), int verbose);
#endif
