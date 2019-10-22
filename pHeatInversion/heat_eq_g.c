/* Source code for the G operator in the case of heat equation inversion.
 * Straightforward case: simples heat equation on [0,1]
 * with zero boundary conditions. u = u(x, t)
 *     -d^2/dx u = du/dt on [0,1] for every time t>0
 *       u(x, 0) = u_D   on [0,1]
 *       u(0, t) = u(1, t) = 0 at every time t
 *
 * The G operator (bayesian notation as in the README file) is here so defined:
 *       - express u_0, starting condition (not necessarely known), 
 *         as basis expansion
 *         by using Fourier. Say that we stop at n = 3;
 *         It's our domain dimension;
 *       - by using the approximated u_0, solve the PDE and register
 *         the results at a fixed time, say 0.01
 *       - set, as output y, various space values of u_sol,
 *         say at x=0, 0.1, ..., x=1 (again, time has been fixed).
 *
 * Summing up we have the map:
 * G: R^basis_expansion -> R^number_of_spatial_observations_at_time_0.01
 * 
 * and our aim will be to reverse it:
 * reconstruct an approximative initial condition by observing the y datas. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "myblas.h"

/* Eigenalues and eigenvectors of laplacian -dx/dx2 on [0,1]
 * with 0-boundary conditions
 * be careful with indeces: 0 in array, 1 in mathematics */
double alpha(int j)
{
        ++j;
        return pow((2. * M_PI * j), 2.);
}
double phi(int j, double x)
{
        ++j;
        return sin(2. * M_PI * j * x);
}


/* Solver is the G operator as described above (before the headers)
 * Parameters:
 - a is an array of dimension basis_expansion, which contains the first
   #basis_expansion coefficients of the initial condition given in input;
 - obs_number describes the number of space observations done at the solution
   at time 0.01. Must be strictly higher that 1,
   reflecing the idea that the boundary 0 is always observed
   (and should be 0; it would be a numerical confirm). 
 - y array of dimension obs_number that will contain the output just described.
*/
void G(const double *a, int basis_expansion, double *y, int obs_number)
{
        assert(basis_expansion >= 1);
        assert(obs_number > 1);
        assert(a != NULL && y != NULL);

        double time_limit = 0.01;

        /* Set h, the space-step, in a way to perform a number of (equally
         * spaced) observations equal to obs_number.
         * Eg: if obs_number= 11, h is then 0.1, allowing a total of
         * 11 observations: 0, 0.1, 0.2, ... , 1. */
        double h = 1.0 / (obs_number - 1);
        double tmp_sum = 0;
        int i = 0;
        int j = 0;

        /* Having the starting condition Fourier's coefficient,
         * to compute the solution is straightforward:
         * balance them with the exponentian of eigenvalues */
        for (i = 0; i < obs_number; ++i){
                tmp_sum = 0;
                for (j = 0; j < basis_expansion; ++j){
                        /* Solution's formula */
                        tmp_sum +=
                                a[j] * exp(-alpha(j) * time_limit)
                                * phi(j, h * i);
                }
                y[i] = tmp_sum;
        }
}

double higherThree(const double *x, int dim){
        (void) dim;
        return sin(x[0]) * cos(x[1]);
//        return(nrm2(x, dim) > 14. ? 1 : 0);
}
