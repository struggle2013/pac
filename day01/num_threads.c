#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
void hello(void);//thread function
int main(int argc,char*argv[]){
	//Get number of threads from command line
	int thread_count = strtol(argv[1],NULL,10);
# pragma omp parallel num_threads(thread_count)
	hello();
	return 0;
}
void hello(){
	int my_rank = omp_get_thread_num();
	int thread_count = omp_get_num_threads();
	printf("Hello from thread %d of %d \n",my_rank,thread_count);
}
