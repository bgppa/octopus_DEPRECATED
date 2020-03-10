#ifndef _ODE_H_
#define _ODE_H_

double euler (double (*f) (double, double), double y0, double T, int N,
              double (*true_sol) (double), int verbose);

double midpoint (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double), int verbose);

double rkfourth (double (*f) (double, double), double y0, double T, int N,
                 double (*true_sol) (double), int verbose);

double euler_d	(void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, double *full_dynamic, 
		 const double *M1, double (*U) (int, const double *), 
		 int verbose);

double midpoint_d(void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, double *full_dynamic,
		 const double *M1, double (*U) (int, const double *), 
		 int verbose);

double rkfourth_d (void (*F) (int, double, const double *, double*), int d,
                 double *y0, double T, int N, double *full_dynamic,
		 const double *M1, double (*U) (int, const double *), 
		 int verbose);

double verlet (double *x, int d, double T, int N, double *full_dynamic,
		const double *M1, double(*U)(int, const double*), int verbose);
#endif
