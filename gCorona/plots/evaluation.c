#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../g.c"


int main(int argc, char** argv) {
	int num_params = 2;
	int num_days = 50;
	double *t = malloc(sizeof(double) * num_days);
	double *params = malloc(sizeof(double) * num_params);

	params[0] = atof(argv[1]);
	params[1] = atof(argv[2]);
	glob_initCond = atof(argv[3]);

        G(params, num_params, t, num_days);

	for (int i = 0; i < num_days; ++i) {
		printf("%d %.f\n", i+1, t[i]);
//		printf("%d %.f\n", (i + 11) <= 31? i + 11 :
//			       	(i+11)	% 32 + 1, t[i]);
//		printf("%d %.f\n", i - 15, t[i]);
	}
	free(params);
	free(t);
        return 0;
} 
