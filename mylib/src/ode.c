#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "myblas.h"
#include "ode.h"

/* ODE: 1D, ODE: nD, ODE: Hamiltonian Systems (Verlet) */
	
/* Here there are three 1-dimensional solver for the ODE
 * y'(x) = f(x, y(x)) on the interval [0, T]
 * y0 is the initial condition, y(0), while N the dicretization value:
 * h = T/N. A "true_solution" argument must be given: if NULL,
 * nothing happens, otherwise the local error's sum are computed.
 * The functions prints couples like (x, y(x)) and return
 * the numerical value of y at time T. */

/* First possible algorithm: Euler's method */
double euler (double (*f) (double, double), double y0, double T, int N,
	      double (*true_sol) (double), int verbose)
{
	double h = T / (double) N;
	double tmp_t = 0; /* Time step */
	double y_n = y0;
	double x_n = 0;
	double err = 0;
	printf("%f %f\n", tmp_t, y_n);
	for (int i = 0; i < N; ++i) {
		y_n += h * f(x_n, y_n);
		x_n += h;
		tmp_t += h;
		if (true_sol != NULL) {
			err += fabs(y_n - true_sol(x_n));
		}
		if (verbose) {
			printf("%f %f\n", tmp_t, y_n);
		}
	}
	if (true_sol != NULL) {
		printf("Euler (sum of local) errors: %e\n", err);
	}
	return y_n;
}

/* Midpoint's rule, i.e. Runge-Kutta of order 2 */
/* Return the value of y at time T */
double midpoint (double (*f) (double, double), double y0, double T, int N, 
		 double (*true_sol) (double), int verbose)
{
	double h = T / (double) N;
	double tmp_t = 0;
	double y_n = y0;
	double x_n = 0;
	double err = 0;
	double k1 = 0;
	double k2 = 0;
	if (verbose) {
		printf("%f %f\n", tmp_t, y_n);
	}
	for (int i = 0; i < N; ++i) {
		/* Update the value */
		k1 = h * f(x_n, y_n);
		k2 = h * f(x_n + h / 2., y_n + k1 / 2.);
		y_n += k2;
		x_n += h;
		tmp_t += h;
		if (true_sol != NULL) {
			err += fabs(y_n - true_sol(x_n));
		}
		if (verbose) {
			printf("%f %f\n", tmp_t, y_n);
		}
	}
	if (true_sol != NULL) {
		printf("Midpoint (sum of local) errors: %e\n", err);
	}
	return y_n;
}

/* Runge-Kutta algorithm of order 4. Return the value of y at time T */
double rkfourth (double (*f) (double, double), double y0, double T, int N, 
		 double (*true_sol) (double), int verbose)
{
	double h = T / (double) N;
	double tmp_t = 0;
	double y_n = y0;
	double x_n = 0;
	double err = 0;
	double k1 = 0;
	double k2 = 0;
	double k3 = 0;
	double k4 = 0;
	if (verbose) {
		printf("%f %f\n", tmp_t, y_n);
	}
	for (int i = 0; i < N; ++i) {
		/* Update the value */
		k1 = h * f(x_n, y_n);
		k2 = h * f(x_n + h / 2., y_n + k1 / 2.);
		k3 = h * f(x_n + h / 2., y_n + k2 / 2.);
		k4 = h * f(x_n + h, y_n + k3);
		y_n += k1 / 6. + k2 / 3. + k3 / 3. + k4 / 6.;
		x_n += h;
		tmp_t += h;
		if (true_sol != NULL) {
			err += fabs(y_n - true_sol(x_n));
		}
		if (verbose) {
			printf("%f %f\n", tmp_t, y_n);
		}
	}
	if (true_sol != NULL) {
		printf("rk4 (sum of local) errors: %e\n", err);
	}
	return y_n;
}

/* ------------------ END 1-DIM ODES ------------------------------ */

/* ----- AUXILIARY FUNCTIONS FOR HAMILTONIAN DYNAMICS -------------- */
/* Evaluate the -gradient of a function U : R^d -> R
 * d = domain's dimension, i.e...
 * ...dimension of the array q: point in which to evaluate;
 * ...dimension of the array grad, which will contain the results;
 * U is the function whose gradient is of interest */
void minus_gradient_U (int d, const double *q, double *grad,
		double (*U) (int d, const double *x))
{
	double h = 0.01;	/* Finite difference's step */
	double *q_plus_h = malloc(sizeof(double) * d);
	assert(q_plus_h != NULL);
	copy(q, q_plus_h, d);
	/* Do a finite difference evaluation on every i-th coordinate */
	for (int i = 0; i < d; ++i ) {
		q_plus_h[i] += h;			/* Add h... */
		grad[i] = (U(d, q_plus_h)-U(d, q)) / h;	/* FinDiff */
		grad[i] = - grad[i]; 			/* MINUS gradient */
		q_plus_h[i] -= h;			/* Restore the value */
	}
	free(q_plus_h);
}

/* An Hamiltonian system of the form of our interest, is completely
 * specified by two elements:
 * a MASS_MATRIX calles ham_M1, of dimension dim. It is usually
 * just the identity, otherwise we go into preconditined systems;
 * a POTENTIAL energy U : R^dim -> R.
 * Given these data, we move our particles in dimension 2*dim.
 * The MASS_MATRIX determine the KYNETIC energy via multiplication:
 * T (p) : R^dim -> R, defined as 1/2 p MASS p
 * Then, given the KYNETIC and the POTENTIAL, the Hamiltonian is their sum:
 * H(q, p) = T(p) + U(q) : R ^ 2*dim -> R
 * 
 * The solvers above (Euler, Midpoint, Runge-Kutta) are suitable for
 * ODE on the form dy(t) = F(t, y(t)).
 * So the question is: given the Hamiltonian, how can I build F
 * and then use the solvers above?
 * By *definition*, an hamiltinian system is giverned by F having a
 * specific gradiend form w.r.t. H, here not specified again.
 *
 * Summing up:
 * M1 + U = T + U = H = F = ODE_SOLVER(F)
 * In other words, given M1 and U we can authomatically solve an hamiltonian
 * system by expoiting the solvers above.
 * The only inconvenience is that we have to keep M1 and U as
 * global data, since they are need to build the subsequent functions.
 * What follow now is the C-translation of what just described. */

/* Compute the KYNETIC ENERGY T: 1/2 * p * M1 * p */
/* d : dimension of p;
 * p : array of length d;
 * M1 : matrix of dimension d * d. Mathematically it must be the INVERSE
 * of M, mass matrix defined in the main problem. The function *assumes*
 * to have M1 already as inverse of M, instead of inverting it every time */
double ham_T (int d, const double *p, const double *M1)
{
	assert(M1 != NULL);
	double *tmp = malloc(sizeof(double) * d);
	assert(tmp != NULL);
	fillzero(tmp, d);
	double result = 0.;
	/* tmp = M1 * p */
	for (int i = 0; i < d; ++i) {
		for (int j = 0; j < d; ++j) 
			{tmp[i] += p[j] * M1[i * d + j];}
	}
	/* result = p * tmp = p * M1 * p */
	for (int i = 0; i < d; ++i) 
		{result += p[i] * tmp[i];}
	free(tmp);
	return result / 2.;
}


/* Compute the Hamiltonian of a dynamical state x = (q, p).
 * d = total dimension, equal to 2 * dim (by construction of ham. system)
 * x is the current state, think it as (q, p), (position, momentum), dimension d
 * M1 the inverse of the mass matrix, dimension d * d
 * U is a function R^dim -> R describing the potential energy */
double ham_H (int d, const double *x, const double *M1,
		double (*pot_U) (int, const double *))
{
	assert(d % 2 == 0);
	assert(M1 != NULL);
	assert(pot_U != NULL);
	double result = 0.;
	/* Kynetic energy is computed using only p, so the second half of x */
	result += ham_T(d / 2, x + d / 2, M1);
	/* Potential energy is w.r.t. q, i.e. first half components of x */
	result += pot_U(d / 2, x);
	return result;
}


/* Called not from an external user, no need to check input validity.
 * d2 is the dimension on y, array splitted into (q, p) according to
 * the hamiltonian notation. res is an array of dimension d2, too,
 * while M1 the inverted mass matrix of dimension d.
 * U is a function going from R^d to R.
 * Remember: the ODE standard notation is y' = F(y).
 * Well, this function return F given U and M1, standard in the Hamiltonian
 * case. In this way we can solve hamiltonian system with Euler and RK too */
void hamiltonianF (int d2, const double *y, double *res, const double *M1,
		double (*U) (int, const double*) )
{
	int d = d2 / 2;
	/* Compute the product M1 * p, where recall p to be the second
	 * half of y, and store the result in the first half of res.
	 * Enough writing d, since res has total dimension d2 */
	axiny (M1, y + d, res, d);
	/* Compute -gradU(q), where q is the first half of the vector y,
	 * store then the result in the second half of res */
	minus_gradient_U (d, y, res + d, U);
}


/* Starting now the support for high-dimensional ODEs */
/* The value of y at time T is written on y0, array with initial conditions */
double euler_d (void (*F) (int, double, const double *, double*), int d,
		double *y0, double T, int N, double *full_dynamic,
		const double *M1, double (*U) (int, const double *), 
		int verbose )
{
	double delta_h = 1.;
	int ordinary_mode = 1; /* 1: ordinary, 0: hamiltonian mode */
	if (F == NULL) {
		/* Then the F must be given from U and M1, and we are in
		 * the hamiltonian case */
		assert(U != NULL);
		assert(M1 != NULL);
		assert(d / 2 == 0);
		ordinary_mode = 0;
		delta_h = - ham_H(d, y0, M1, U);
	}
	/* need: y_n, k1 tmp_eval, tmp_sum */
	double *y_n = malloc(sizeof(double) * d);
	double *k1 = malloc(sizeof(double) * d);
	double *tmp_eval = malloc(sizeof(double) * d);
	assert(y_n != NULL);
	assert(k1 != NULL);
	assert(tmp_eval != NULL);

	double t_n = 0;
	double h = T / (double) N;
	
	/* Assigning the initial conditions */
	for (int i = 0; i < d; ++i) {
		y_n[i] = y0[i];
		k1[i] = 0.;
	}	
	/* Printing the initial coditions */
	if (verbose) {
		printf("%f ", t_n);
		for(int i = 0; i < d; ++i) {
			printf("%e%c", y_n[i], i != (d - 1) ? ' ' : '\n');
		}
	}
	/* Storing the initial conditions */
	if (full_dynamic != NULL) {
		full_dynamic[0] = t_n;
		copy(y_n, full_dynamic + 1, d);
	}
	/* Computing y_n iteratively according to Euler */
	for(int n = 0; n < N; ++n) {
		/* Computing k1 */
		if (ordinary_mode) {
			F(d, t_n, y_n, tmp_eval);
		} else {
			hamiltonianF(d, y_n, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k1[i] = h * tmp_eval[i];
		}
		/* Now that I have k1compute y_(n+1) */
		for (int i = 0; i < d; ++i) {
			y_n[i] += k1[i];
		}
		t_n += h;
		if (verbose) {
			printf("%f ", t_n);
			for(int i = 0; i < d; ++i) {
				printf("%e%c", y_n[i], i != (d-1) ? ' ' : '\n');
			}	
		}
		if (full_dynamic != NULL) {
			full_dynamic[(n+1) * (d+1)] = t_n;
			copy(y_n, full_dynamic + (n+1) * (d + 1) + 1, d);
		}
//		getchar();
	}
	/* Copy the value into y0 */
	copy (y_n, y0, d);
	if (ordinary_mode == 0) {
		delta_h += ham_H (d, y0, M1, U);
	}
	free(k1);
	free(y_n);
	return delta_h;
}

/* The value of y at time T is written on y0, array with initial conditions */
double midpoint_d (void (*F) (int, double, const double *, double*), int d,
	         double *y0, double T, int N, double *full_dynamic, 
		 const double *M1, double (*U) (int, const double *), 
		 int verbose)
{
	double delta_h = 1.;
	int ordinary_mode = 1; /* 1: ordinary, 0: hamiltonian mode */
	if (F == NULL) {
		/* Then the F must be given from U and M1, and we are in
		 * the hamiltonian case */
		assert(U != NULL);
		assert(M1 != NULL);
		assert(d / 2 == 0);
		ordinary_mode = 0;
		delta_h = - ham_H (d, y0, M1, U);
	}

	/* need: y_n, k1, k2, k3, k4, F, tmp_eval, tmp_sum */
	double *y_n = malloc(sizeof(double) * d);
	double *k1 = malloc(sizeof(double) * d);
	double *k2 = malloc(sizeof(double) * d);
	double *tmp_eval = malloc(sizeof(double) * d);
	double *tmp_sum = malloc(sizeof(double) * d);
	assert(y_n != NULL);
	assert(k1 != NULL && k2 != NULL);
	assert(tmp_eval != NULL && tmp_sum != NULL);

	double t_n = 0;
	double h = T / (double) N;
	
	/* Assigning the initial conditions */
	for (int i = 0; i < d; ++i) {
		y_n[i] = y0[i];
		k1[i] = 0.;
		k2[i] = 0.;
	}	

	/* Printing the initial coditions */
	if (verbose) {
		printf("%f ", t_n);
		for(int i = 0; i < d; ++i) {
			printf("%e%c", y_n[i], i < (d - 1) ? ' ' : '\n');
		}
	} 
	if (full_dynamic != NULL) {
		full_dynamic[0] = t_n;
		copy(y_n, full_dynamic + 1, d);
	}
	/* Computing y_n iteratively according to Rounge-Kutta */
	for(int n = 0; n < N; ++n) {
		/* Computing k1 */
		if (ordinary_mode) {
			F(d, t_n, y_n, tmp_eval);
		} else {
			hamiltonianF(d, y_n, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k1[i] = h * tmp_eval[i];
		}
		/* Computing k2 */
		for (int i = 0; i < d; ++i) {
			tmp_sum[i] = y_n[i] + k1[i] / 2.;
		}
		F(d, t_n + h / 2., tmp_sum, tmp_eval);
		for (int i = 0; i < d; ++i) {
			k2[i] = h * tmp_eval[i];
		}
		/* Now that I have y_n, k1, k2, compute y_(n+1) */
		for (int i = 0; i < d; ++i) {
			y_n[i] += k2[i];
		}
		t_n += h;
		if (verbose) {
			printf("%f ", t_n);
			for(int i = 0; i < d; ++i) {
				printf("%e%c", y_n[i], i < (d - 1) ? ' ' :'\n');
			}
		}
		if (full_dynamic != NULL) {
			full_dynamic[(n+1) * (d+1)] = t_n;
			copy(y_n, full_dynamic + (n+1) * (d + 1) + 1, d);
		}

	}
	/* Copy the current value into y0 */
	copy (y_n, y0, d);
	if (ordinary_mode == 0) {
		delta_h += ham_H (d, y0, M1, U);
	}
	free(k1);
	free(k2);
	free(y_n);
	free(tmp_sum);
	free(tmp_eval);
	return delta_h;
}

/* The value of y at time T is written on y0, array with initial conditions */
double rkfourth_d (void (*F) (int, double, const double *, double*), int d,
		double *y0, double T, int N, double *full_dynamic,
		const double *M1, double (*U) (int, const double*), int verbose)
{
	double delta_h = 1.;
	int ordinary_mode = 1; /* 1: ordinary, 0: hamiltonian mode */
	if (F == NULL) {
		/* Then the F must be given from U and M1, and we are in
		 * the hamiltonian case */
		assert(U != NULL);
		assert(M1 != NULL);
		assert(d % 2 == 0);
		ordinary_mode = 0;
		delta_h = - ham_H (d, y0, M1, U);
	}

	/* need: y_n, k1, k2, k3, k4, F, tmp_eval, tmp_sum */
	double *y_n = malloc(sizeof(double) * d);
	double *k1 = malloc(sizeof(double) * d);
	double *k2 = malloc(sizeof(double) * d);
	double *k3 = malloc(sizeof(double) * d);
	double *k4 = malloc(sizeof(double) * d);
	double *tmp_eval = malloc(sizeof(double) * d);
	double *tmp_sum = malloc(sizeof(double) * d);
	assert(y_n != NULL);
	assert(k1 != NULL && k2 != NULL && k3 != NULL && k4 != NULL);
	assert(tmp_eval != NULL && tmp_sum != NULL);

	double t_n = 0;
	double h = T / (double) N;
	
	/* Assigning the initial conditions */
	copy (y0, y_n, d);
	fillzero (k1, d);
	fillzero (k2, d);
	fillzero (k3, d);
	fillzero (k4, d);

	/* Printing the initial coditions */
	if (verbose) {
		printf("%f ", t_n);
		for(int i = 0; i < d; ++i) {
			printf("%e%c", y_n[i], i < (d - 1) ? ' ' : '\n');
		}
	}
	if (full_dynamic != NULL) {
		full_dynamic[0] = t_n;
		copy(y_n, full_dynamic + 1, d);
	}
	/* Computing y_n iteratively according to Rounge-Kutta */
	for(int n = 0; n < N; ++n) {
		/* Computing k1 */
		if (ordinary_mode) {
			F(d, t_n, y_n, tmp_eval);
		} else {
			hamiltonianF(d, y_n, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k1[i] = h * tmp_eval[i];
		}
		/* Computing k2 */
		for (int i = 0; i < d; ++i) {
			tmp_sum[i] = y_n[i] + k1[i] / 2.;
		}
		if (ordinary_mode) {
			F(d, t_n + h / 2., tmp_sum, tmp_eval);
		} else {
			hamiltonianF(d, tmp_sum, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k2[i] = h * tmp_eval[i];
		}
		/* Cumputing k3 */
		for (int i = 0; i < d; ++i) {
			tmp_sum[i] = y_n[i] + k2[i] / 2.;
		}
		if (ordinary_mode) {
			F(d, t_n + h / 2., tmp_sum, tmp_eval);
		} else {
			hamiltonianF(d, tmp_sum, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k3[i] = h * tmp_eval[i];
		}
		/* Computing k4 */
		for (int i = 0; i < d; ++i) {
			tmp_sum[i] = y_n[i] + k3[i];
		}
		if (ordinary_mode) {
			F(d, t_n + h, tmp_sum, tmp_eval);
		} else {
			hamiltonianF(d, tmp_sum, tmp_eval, M1, U);
		}
		for (int i = 0; i < d; ++i) {
			k4[i] = h * tmp_eval[i];
		}
		/* Now that I have y_n, k1, k2, k3, k4, compute y_(n+1) */
		for (int i = 0; i < d; ++i) {
			y_n[i] += k1[i]/6. + k2[i]/3. + k3[i]/3. + k4[i]/6.;
		}
		t_n += h;
		if (verbose) {
			printf("%f ", t_n);
			for(int i = 0; i < d; ++i) {
				printf("%e%c",y_n[i], i < (d - 1) ? ' ' : '\n');
	 		}
		}
		if (full_dynamic != NULL) {
			full_dynamic[(n+1) * (d+1)] = t_n;
			copy(y_n, full_dynamic + (n+1) * (d + 1) + 1, d);
		}

	}
	/* Copy the current value into y0 */
	copy (y_n, y0, d);
	if (ordinary_mode == 0) {
		delta_h += ham_H (d, y0, M1, U);
	}
	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(y_n);
	free(tmp_sum);
	free(tmp_eval);
	return delta_h;
}


/* TO REWRITE */
/* Verlet integrator as described at page 129, Nawaf's paper.
 * Parameters: x of dimension d containing the starting
 * conditions, which will be overwritten with the system status
 * at time T. N is the number of steps. IMPORTANT:
 * defining the Hamiltonian system. Return the difference
 * of the hamiltonian between initial and final point */
double verlet (double *x, int d, double T, int N, double *full_dynamic,
		const double *M1, double (*U) (int, const double*), int verbose)
{
	int d_half = d / 2;
	double h = T / (double) N;
	/* decouple x into q and p, in a way that it's easier to work with */
	double *q_n = malloc(sizeof(double) * d_half);
	double *p_n = malloc(sizeof(double) * d_half);
	double *p_half = malloc(sizeof(double) * d_half);
	double *tmp_eval = malloc(sizeof(double) * d_half);
	assert(q_n != NULL && p_n != NULL);
	assert(p_half != NULL && tmp_eval != NULL);
	assert(M1 != NULL);
	assert(U != NULL);

	/* Set the initial conditions */
	for (int i = 0; i < d; ++i) {
		if (i < d_half) { q_n[i] = x[i]; }
			else	{ p_n[i - d_half] = x[i]; }
	}
	double dt = 0;


	/* Set the initial hamiltonian value */
	double initial_ham = ham_H(d, x, M1, U);

	/* If verbose : */
	if (verbose) {
		/* Print the initial conditions */
		printf("%f ", dt);
		for (int i = 0; i < d; ++i) {
			printf("%e%c", i < d_half ? q_n[i] : p_n[i - d_half],
					i < d - 1 ? ' ' : '\n');
		}
	}

	if (full_dynamic != NULL) {
		full_dynamic[0] = dt;
		for (int i = 0; i < d; ++i) {
			if (i < d_half) { full_dynamic[i +1] = q_n[i]; }
			else		{ full_dynamic[i+1] = p_n[i - d_half]; }
		}
	}
	/* Solve the system */
	for (int i = 0; i < N; ++i) {
		dt += h;
		/* Step 1: update p_half */
		minus_gradient_U(d/2, q_n, tmp_eval, U);
		for (int i = 0; i < d / 2; ++i) {
			p_half[i] = p_n[i] + h / 2. * tmp_eval[i];
		} 
		/* Step 2: update q_n */
		for (int i = 0; i < d / 2; ++i) {
			tmp_eval[i] = 0;
			for(int j = 0; j < d / 2; ++j) {
				tmp_eval[i] += M1[i*(d/2) + j] * p_half[j];
			}
		}
		for (int i = 0; i < d_half; ++i) {
			q_n[i] += h * tmp_eval[i];
		}
		/* Step 3: update p_n */
		minus_gradient_U(d/2, q_n, tmp_eval, U);
		for (int i = 0; i < d / 2; ++i) {
			p_n[i] = p_half[i] + h / 2. * tmp_eval[i];
		}
		/* If verbose, print the status */
		if (verbose) {
			printf("%f ", dt);
			for (int i = 0; i < d; ++i) {
				printf("%e%c", i < d_half ? q_n[i] : 
						p_n[i - d_half],
						i < d - 1 ? ' ' : '\n');
			}
		}
		
		if (full_dynamic != NULL) {
			full_dynamic[(i + 1) * (d + 1)] = dt;
			for (int k = 0; k < d; ++k) {
				if (k < d_half) { 
				full_dynamic[k + (i+1) * (d+1) + 1] = q_n[k]; 
				} else {
				full_dynamic[k + (i+1)*(d+1) + 1] = p_n[k - d_half];
				}
			}
		}	
	} /* System solved: copy the solution into x */
	for (int i = 0; i < d; ++i) {
		if (i < d_half ) { x[i] = q_n[i]; }
			else	 { x[i] = p_n[i - d_half]; }
	}
	free(q_n);
	free(p_n);
	free(p_half);
	free(tmp_eval);
	return ham_H(d, x, M1, U) - initial_ham;
}
