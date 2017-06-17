#include<stdio.h>
#include<omp.h>
int main(int argv, char* argc[]){
    int tid=0;
    printf("start threads %d\n",tid);
#pragma omp parallel private(tid)
    {
        //omp_set_num_threads(10);
        tid = omp_get_num_threads();
        printf("omp threads num is %d \n",tid);
        sleep(10);
    }
    printf("end threads %d\n",tid);
    return 0;
}
