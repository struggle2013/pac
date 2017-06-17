#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <strings.h>
#include <immintrin.h>

#define REALTYPE double

inline double dtime() // Return double style timestamp using gettimeofday().
{
    double tseconds;
    struct timeval mytime;

    gettimeofday(&mytime, NULL);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return (tseconds);
} /* dtime() */

void findcov_(int *NVAR, int *NROW, int *NV, REALTYPE *MType, REALTYPE *COV)
{
    // Readonly MType[*NRWO * *NVAR];
    // Output COV[*NV * *NVAR];

    long IVar1, IVar2, IVar3;

    int nv_i = *NV;
    int nvar_i = *NVAR;
    int nrow_i = *NROW;
    int nrow_l = ((nrow_i+15)>>4)<<4;
    int nvar_l = nvar_i + nv_i;
    printf("%d %d\n",nvar_l,nrow_l);
    REALTYPE *MT = (REALTYPE *)_mm_malloc((long)nvar_l*nrow_l*sizeof(REALTYPE), 64);
#pragma omp parallel private(IVar2, IVar3)
{
    REALTYPE t0 = dtime();
#pragma omp for  
    for (IVar1=0; IVar1<nvar_l; IVar1++) {
        long I2 = IVar1 % nvar_i;
        long VR = IVar1*nrow_l;
        for (IVar3=0; IVar3<nrow_i; IVar3++) {
            REALTYPE MTemp = MType[IVar3*nvar_i+I2];
            MT[VR+IVar3] = (MTemp!=MTemp || fabs(MTemp-65.0)<=0.001) ? 65.0 : MTemp;
        }
    }
    REALTYPE t1 = dtime();
#pragma omp for  
    for (IVar1=0; IVar1<nvar_i; IVar1++) {
        for (IVar2=0; IVar2<nv_i; IVar2++) {
            int NumMissing = 0;
            REALTYPE MeanX1 = 0.0f;
            REALTYPE MeanX2 = 0.0f;
            REALTYPE C = 0.0f;
            long VR = IVar1*nrow_l;
            long V1 = VR+(IVar2+1)*nrow_l;
#pragma vector aligned
            for (IVar3=0; IVar3<nrow_i; IVar3++) {
                REALTYPE M1 = MT[VR+IVar3];
                REALTYPE M2 = MT[V1+IVar3];
                if (M1!=65.0 && M2!=65.0) {
//                if (M1==M1 && fabs(M1-65.0f)>0.001f && M2==M2 && fabs(M2-65.0f)>0.001f) {
                    NumMissing++;
                    MeanX1+=M1;
                    MeanX2+=M2;

                    C+=M1*M2;
                }
            }
            REALTYPE IRep = NumMissing<=0 ? 0.0f : 1.0f/(REALTYPE)NumMissing;
            MeanX1*=IRep;
            MeanX2*=IRep;
            COV[IVar2*nvar_i+IVar1]=C*IRep-MeanX1*MeanX2;
        }
    }
    REALTYPE t2 = dtime();
//    printf("Loop1: %fs, Loop2: %fs.\n", t1-t0, t2-t1);
}
    _mm_free(MT);
}

