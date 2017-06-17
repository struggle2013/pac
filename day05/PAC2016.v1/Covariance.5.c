#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <strings.h>
#ifdef _USEHBW
#include <hbwmalloc.h>
#endif

#define REALTYPE double
#define RB 192

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

    int nv_i = *NV;
    int nvar_i = *NVAR;
    int nrow_i = *NROW;
    int nrow_l = ((nrow_i+15)>>4)<<4;
    int nvar_l = nvar_i + nv_i;
    printf("%d %d\n",nvar_l,nrow_l);
#pragma omp parallel
{
    int IVar1, IVar2, IVar3;

    int i,j,k;
    int nth = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int NB = nvar_i/nth;
    int nmod = nvar_i%nth;
    int svar = tid*NB+((tid>=nmod)?nmod:tid);
    int nvar_l = NB+((tid>=nmod)?0:1);
    int evar = svar+nvar_l;
#ifdef _TIMING
    printf("%d %d %d %d %d %d\n", nth, tid, nvar_l, nmod, svar, evar);
    REALTYPE t0 = dtime();
#endif
#ifdef _USEHBW
    REALTYPE *MT; hbw_posix_memalign((void**)&MT, 64, RB*(nvar_l+nv_i)*sizeof(REALTYPE));
    int *NM; hbw_posix_memalign((void**)&NM, 64, nvar_l*nv_i*sizeof(int));
    REALTYPE *MX1; hbw_posix_memalign((void**)&MX1, 64, nvar_l*nv_i*sizeof(REALTYPE));
    REALTYPE *MX2; hbw_posix_memalign((void**)&MX2, 64, nvar_l*nv_i*sizeof(REALTYPE));
    REALTYPE *CT; hbw_posix_memalign((void**)&CT, 64, nvar_l*nv_i*sizeof(REALTYPE));
#else
    REALTYPE *MT  = (REALTYPE *)_mm_malloc(RB*(nvar_l+nv_i)*sizeof(REALTYPE), 64);
    int *NM = (int *)_mm_malloc(nvar_l*nv_i*sizeof(int), 64);
    REALTYPE *MX1 = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
    REALTYPE *MX2 = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
    REALTYPE *CT = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
#endif
    bzero((void*)NM, nvar_l*nv_i*sizeof(int));
    bzero((void*)MX1, nvar_l*nv_i*sizeof(REALTYPE));
    bzero((void*)MX2, nvar_l*nv_i*sizeof(REALTYPE));
    bzero((void*)CT, nvar_l*nv_i*sizeof(REALTYPE));

#ifdef _TIMING
    REALTYPE t1 = dtime(), t2=0, t3=0;
    printf("%f memory allocation %d\n", t1-t0, tid);
#endif
    REALTYPE *MTT = MType, *MT1, *MT2;
    for (i=0; i<nrow_i; i+=RB) {
        int srow = i;
        int erow = srow+RB;
        if (erow>nrow_i) erow = nrow_i;
        MT2 = MT;

	// Data block transposing.
        for (IVar1=svar; IVar1<evar; IVar1++) {
		for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
			REALTYPE MTemp = MTT[j*nvar_i+IVar1];
			MT2[j] = MTemp;
		}
            MT2+=RB;
        }
        for (IVar1=evar; IVar1<evar+nv_i; IVar1++) {
            for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
                  REALTYPE MTemp = MTT[j*nvar_i+IVar1%nvar_i];
		  MT2[j] = MTemp;
            }
            MT2+=RB;
        }
#ifdef _TIMING
        t0 = dtime();
        t2 += t0-t1;
#endif

	// Block accumulation calculation.
	k=0;
        for (IVar1=svar; IVar1<evar; IVar1++) {
            MT1 = MT+(IVar1-svar)*RB;
			for (IVar2=0; IVar2<nv_i; IVar2++) {
				MT2 = MT1+(IVar2+1)*RB;
				int NumMissing = 0;
				REALTYPE MeanX1 = 0.0;
				REALTYPE MeanX2 = 0.0;
				REALTYPE C = 0.0;
#pragma vector aligned
				for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
					REALTYPE M1 = MT1[j];
					REALTYPE M2 = MT2[j];
					if (M1==M1 && fabs(M1-65.0) > 0.001 && M2==M2 && fabs(M2-65.0)>0.001) {
						NumMissing+=1;
						MeanX1+=M1;
						MeanX2+=M2;
						C+=M1*M2;
					}
				}
				NM[k]+=NumMissing;
				MX1[k]+=MeanX1;
				MX2[k]+=MeanX2;
				CT[k]+=C;
				k++;
            }
        }
#ifdef _TIMING
        t1 = dtime();
        t3 += t1-t0;
#endif
        MTT += nvar_i*RB;
    }

	// Final covarience caculation. 
	j=0;
    for (IVar1=svar; IVar1<evar; IVar1++) {
    	for (IVar2=0; IVar2<nv_i; IVar2++) {
            REALTYPE IRep = NM[j]<=0 ? 0.0 : 1.0/(REALTYPE)NM[j];
            COV[IVar2*nvar_i+IVar1]=(CT[j]-MX1[j]*MX2[j]*IRep)*IRep;
            j++;
        }
    }
#ifdef _TIMING
    printf("Th %d, trans %8.4fs, calc %8.4f, last %8.4f.\n", tid, t2, t3, dtime()-t1);
#endif
#ifdef _USEHBW
    hbw_free(MT);
    hbw_free(NM);
    hbw_free(MX1);
    hbw_free(MX2);
    hbw_free(CT);
#else
    _mm_free(MT);
    _mm_free(NM);
    _mm_free(MX1);
    _mm_free(MX2);
    _mm_free(CT);
#endif
}
}

