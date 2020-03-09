#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myblas.h"
#include "ode.h"
#include "g.h"

int glob_dDom = 2;
int glob_dCod = 2;

void example1 (int d, double x, const double* Y, double *res)
{
	(void) d; /* Here the dimension is already assumed to be 2 */
	res[0] = Y[0] * x;
	res[1] = cos(Y[1]) * x * x;
}

void G(const double* in_cond, int d1, double* yy, int d2)
{
//	(int) d1;
//	(int) d2;
//	copy(in_cond, yy, glob_dDom);
	yy[0] = in_cond[0];
	yy[1] = in_cond[1];
	rkfourth_d (example1, glob_dDom, yy, 2.0, 200, 0);
//	euler_d(example1, d1, y, 2.0, 200);

}
