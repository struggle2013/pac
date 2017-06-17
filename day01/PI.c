#include<stdio.h>
#include<omp.h>
static long num_steps=100000;
double step;
void main(){
	int i;
	double x,pi,sum=0.0,start_time,end_time;
	step=1.0/(double)num_steps;//把单位长度分成num_steps步，step代表每一步有多长
	start_time=omp_get_wtime();
	for(i=1;i<=num_steps;i++){//每一步要执行的内容
		x=(i-0.5)*step;
		sum=sum+4.0/(1.0+x*x);
	}
	pi = step*sum;
	end_time=omp_get_wtime();
	printf("PI = %lf Running time = %lf\n",pi,end_time-start_time);

}
