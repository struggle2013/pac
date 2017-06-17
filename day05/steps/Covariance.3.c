#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <immintrin.h>

#define REALTYPE double

void findcov_(int *NVAR, int *NROW, int *NV, double *MType, double *COV)
{
    // Readonly MType[*NRWO * *NVAR];
    // Output COV[*NV * *NVAR];

    int IVar1, IVar2, IVar3;

#pragma omp parallel private(IVar2, IVar3)
{
    int nv_i = *NV;
    int nvar_i = *NVAR;
    int nrow_i = *NROW;
    int *NumMissing = (int *)_mm_malloc(nv_i*sizeof(int), 64);
    double *MeanX1 = (double *)_mm_malloc(nv_i*sizeof(double), 64);
    double *MeanX2 = (double *)_mm_malloc(nv_i*sizeof(double), 64);
    double *C = (double *)_mm_malloc(nv_i*sizeof(double), 64);
#pragma omp for  
    for (IVar1=0; IVar1<nvar_i; IVar1++) {

        bzero((void*)NumMissing, nv_i*sizeof(int));
        bzero((void*)MeanX1, nv_i*sizeof(double));
        bzero((void*)MeanX2, nv_i*sizeof(double));
        bzero((void*)C, nv_i*sizeof(double));

        for (IVar3=0; IVar3<nrow_i; IVar3++) {
            int CT = 0;
            long VR = IVar3*nvar_i;
            double M1 = MType[VR+IVar1];
            if (M1==M1 && fabs(M1-65.0f)>0.001f) {
                if (IVar1+nv_i < nvar_i) {
                    for (IVar2=0; IVar2<nv_i; IVar2++) {
                        int I2 = IVar1+IVar2+1;
                        double M2 = MType[VR+I2];
                        if (M2==M2 && fabs(M2-65.0f)>0.001f) {
                            NumMissing[IVar2]++;
                            MeanX1[IVar2]+=M1;
                            MeanX2[IVar2]+=M2;

                            C[IVar2]+=M1*M2;
                        }
                    }
                }
                else {
                    for (IVar2=0; IVar2<nv_i; IVar2++) {
                        int I2 = (IVar1+IVar2+1)%nvar_i;
                        double M2 = MType[VR+I2];
                        if (M2==M2 && fabs(M2-65.0f)>0.001f) {
                            NumMissing[IVar2]++;
                            MeanX1[IVar2]+=M1;
                            MeanX2[IVar2]+=M2;

                            C[IVar2]+=M1*M2;
                        }
                    }
                }
            }
        }
        for (IVar2=0; IVar2<nv_i; IVar2++) {
            if (IVar1==0 && IVar2==2) printf("%d %f %f %f\n", NumMissing[2], MeanX1[2], MeanX2[2], C[2]);
            long VR=IVar2*nvar_i+IVar1;
            int I2 = NumMissing[IVar2];
            double IRep = I2<=0 ? 0.0f : 1.0f/(double)I2;
            MeanX1[IVar2]*=IRep;
            MeanX2[IVar2]*=IRep;
            COV[VR]=C[IVar2]*IRep-MeanX1[IVar2]*MeanX2[IVar2];
            if (IVar1==0 && IVar2==2) printf("%f %f %f %f\n", C[2]*IRep, MeanX1[2], MeanX2[2], COV[VR]);
        }
    }
    _mm_free(NumMissing);
    _mm_free(MeanX1);
    _mm_free(MeanX2);
    _mm_free(C);
}
}

