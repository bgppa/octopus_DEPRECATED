#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
	if (argc != 5) {
		printf("./exp alpha beta from_day to_day\n");
		return 0;
	}
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	int from = atoi(argv[3]);
	int to = atoi(argv[4]);
	for (int i = 0; i <= to - from; ++i) {
		printf("%d %.f\n", i + from, exp(a * (double) i + b));
	}
	return 0;
}
