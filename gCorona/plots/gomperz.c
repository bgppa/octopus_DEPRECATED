#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gompertz (double alpha, double N, double N0, double t)
{
	return N * exp(log(N0 / N) * exp(-alpha * t));
}


int main(){
        double param[2] = {1.108736e-01, 2.956381e+04};
        for (int i = 0; i < 30; ++i) {
                printf("%d %f\n", i, gompertz(param[0], param[1], 3, i));
        }
        return 0;
} 
        
