#include<stdio.h>
#include<omp.h>
int main(int argc, char*argv[]){
    int nn,tid;
#pragma omp parallel private(tid)
    {
        tid =omp_get_thread_num();
        printf("helilo world ! from openMp thread %d\n",tid);
    }
}
