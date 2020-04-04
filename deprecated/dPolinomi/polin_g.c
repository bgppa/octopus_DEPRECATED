/* Source code for the G operator for the polynomial case.
 * This is actually one of the basis toy problem on which I will
 * constantly perform my benchmarks.
 * What does G do?
 * So, given a certain fixed dimension N, we think about a polynomial
 * of degree N-1 (so, dimension 1 = only the termine noto, makes no sense)
 * G does: given such a polynomial, evaluate it on various points
 * starting from -10 and increasing accordingly.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "myblas.h"

double polin (const double *coeff, int deg, double x)
/* Given d + 1 coefficients, evaluate the polynomial \sum coeff x^i */
{
        /* Remark: deg + 1 = numero dei coefficienti */
        double sum = 0;
        for (int i = 0; i < deg; ++i) {
                sum += coeff[i] * pow(x, i);
        }
        return sum;
}

void G(const double *a, int degree, double *y, int obs_number)
/* Defining now the observation operator corresponding to see the
 * polynomian in certain points */
{
        for (int i = 0; i < obs_number; ++i) {
                y[i] = polin(a, degree, -3. + i);
        }
}
