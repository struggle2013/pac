#include<stdio.h>
#include<omp.h>
static long num_steps=10000000;//设置总共多事步
double step;
#define NUM_THREADS 10
void main(){
	int i;
	double x,pi,sum[NUM_THREADS],start_time,end_time;
	step=1.0/(double)num_steps;
	omp_set_num_threads(NUM_THREADS);

	start_time=omp_get_wtime();
#pragma omp parallel
	{
		int id;
		double x;
		id=omp_get_thread_num();
		for(i=id,sum[id]=0.0;i<num_steps;i=i+NUM_THREADS){
			x=(i-0.5)*step;
			sum[id]+=4.0/(1.0+x*x);
		}
	}
	for(i=0,pi=0.0;i<NUM_THREADS;i++)
		pi += step*sum[i];
	end_time=omp_get_wtime();
	printf("PI = %lf Running time = %lf\n",pi,end_time-start_time);

}
