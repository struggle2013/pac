#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<omp.h>
void main(){
	const int nTrials = 5;
	for (int trial = 0; trial< nTrials;trial++){
		const double x_upper_bound=1.0+trial;
		const double x_lower_bound=0.0;
		const int nSteps = 1000000000;
		const double dx=(x_upper_bound-x_lower_bound)/nSteps;
		double integral = 0.0;
		double partial_num=0.0;
		const double t0 = omp_get_wtime();
#pragma	omp parallel firstprivate(partial_num) 
		{ 
#pragma omp for
			for(int i=0;i<nSteps;i++){
				const double x = x_lower_bound+dx*(double(i)+0.5);
				partial_num += 1.0/sqrt(x)*dx;
			}
#pragma omp atomic
			integral += partial_num;
		}
		const double t1 = omp_get_wtime();
		const double analytical_result=(sqrt(x_upper_bound)-sqrt(x_lower_bound))*2;
		const double numerical_result = integral;
		printf("Result = %.8f (should be %.8f; err = %.8f) Time =%.8f ms\n",
				numerical_result,analytical_result,analytical_result-numerical_result,(t1-t0)*1000.0);

	}
}
