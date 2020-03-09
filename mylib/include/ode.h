#ifndef _ODE_H_
#define _ODE_H_

double euler (double (*f) (double, double), double y0, double T, int N,
              double (*true_sol) (double), int verbose);

double midpoint (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double), int verbose);

double rkfourth (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double), int verbose);

void euler_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, int verbose);

void midpoint_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, int verbose);

void rkfourth_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, int verbose);

/* Variables for the hamiltonian case */
//extern double* ham_M1;
//extern double (*ham_U) (int dim, const double *x);
//void ham_F (int dim, double t, const double *x, double *y);
double verlet (double *x, int d, double T, int N, 
		const double *M1, double(*U)(int, const double*), int verbose);
#endif
