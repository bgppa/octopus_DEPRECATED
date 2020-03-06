#include <stdio.h>
#include <stdlib.h>
#include "g.h"
int glob_dDom = 1;      /* Domain's dimension */
int glob_dCod = 1;      /* Codomain's dimension */
void G (const double* x, int d1, double *y, int d2) {
        (int) d1;
        (int) d2;
        y[0] = x[0] * x[0] * x[0];
}


