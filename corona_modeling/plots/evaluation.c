#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "myblas.h"
//#include "mylapack.h"
//#include "ranvar.h"
#include "../g.c"

#define NOISED 0

int main(int argc, char** argv) {

	if (argc != 6) {
		printf("Syntax: ./richard q Q vi X0 tot_days\n");
		return -1;
	}
#if NOISED
	srand(time(NULL));
#endif
	int num_params = 3;
	double *params = malloc(sizeof(double) * num_params);

	params[0] = atof(argv[1]);
	params[1] = atof(argv[2]);
	params[2] = atof(argv[3]);
	glob_initCond = atof(argv[4]);
	int num_days = -1;
	num_days = atoi(argv[5]);
	double *t = malloc(sizeof(double) * (num_days+1));

	for (int i = 0; i <= num_days; ++i) {
//		t[i] = gompertz_der (params[0], params[1], glob_initCond, i);
		t[i] = richard (params[0], params[1], 
				glob_initCond, params[2], (double) i);
//		t[i] = logistic_der (params[0], params[1], glob_initCond, i);
//		t[i] = gompertz (params[0], params[1], glob_initCond, i);
//		t[i] = logistic (params[0], params[1], glob_initCond, i);
#if NOISED
		t[i] += rndmGaussian(0., t[i], NULL);
#endif
	}
	

	for (int i = 0; i <= num_days; ++i) {
		printf("%d %.f\n", i+9, t[i]);
//		printf("%d %.f\n", (i + 11) <= 31? i + 11 :
//			       	(i+11)	% 32 + 1, t[i]);
//		printf("%d %.f\n", i - 15, t[i]);
	}
	free(params);
	free(t);
        return 0;
} 
