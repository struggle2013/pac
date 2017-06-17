#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <strings.h>
#ifdef _USEHBW
#include <hbwmalloc.h>
#endif
#include <immintrin.h>

#define REALTYPE double
#define RB 192

inline void Trans_manual(double *mtt, double *mt2, int nvar) {
	register int i;
	register __m512d m65 = _mm512_set1_pd(65.0);
	register __m512d mep = _mm512_set1_pd(0.001);
	register __mmask8 mask65;
	register __m512d mt;
	register __m512d m2;
	register __m256i vi = _mm256_set_epi32(7*nvar, 6*nvar, 5*nvar, 4*nvar, 3*nvar, 2*nvar, 1*nvar, 0);
#pragma unroll(8)
	for (i=0; i<RB ; i+=8) {
		mt = _mm512_i32gather_pd(vi, (void *)(mtt+i*nvar), 8);
		m2 = _mm512_sub_pd(mt, m65);
		m2 = _mm512_abs_pd(m2);
		mask65 = _mm512_cmp_pd_mask(m2, mep, 26);
		mt = _mm512_mask_mov_pd(mt, mask65, m65);
		_mm512_store_pd((void*)(mt2+i),mt);
	}
}

inline void COV_manual(double *mt1, double *mt2, double *mx1, double *mx2, double *covt, double *numt) {
	register int i;
	register __m512d m65 = _mm512_set1_pd(65.0);
	register __m512d mone = _mm512_set1_pd(1.0);
	register __m512d x1 = _mm512_setzero_pd();
	register __m512d x2 = _mm512_setzero_pd();
	register __m512d c = _mm512_setzero_pd();
	register __m512d nm = _mm512_setzero_pd();
	register __m512d m1;
	register __m512d m2;
	register __mmask8 maskm2;
	register __mmask8 maskm1;
	register __mmask16 mask;
#pragma unroll(8)
	for (i=0; i<RB ; i+=8) { 
		m1 = _mm512_load_pd(mt1 + i);
		maskm1 = _mm512_cmp_pd_mask(m1, m65, 4);
		m2 = _mm512_load_pd(mt2 + i);
		maskm2 = _mm512_cmp_pd_mask(m2, m65, 4);
		mask = _mm512_kand(maskm1, maskm2);
		x1 = _mm512_mask_add_pd(x1, mask, x1, m1);
		x2 = _mm512_mask_add_pd(x2, mask, x2, m2);
		c = _mm512_mask3_fmadd_pd(m1, m2, c, mask);
		nm = _mm512_mask_add_pd(nm, mask, nm, mone);
	}
	*mx1 += _mm512_reduce_add_pd(x1);
	*mx2 += _mm512_reduce_add_pd(x2);
	*covt += _mm512_reduce_add_pd(c);
	*numt += _mm512_reduce_add_pd(nm);
}

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
    REALTYPE *NM; hbw_posix_memalign((void**)&NM, 64, nvar_l*nv_i*sizeof(REALTYPE));
    REALTYPE *MX1; hbw_posix_memalign((void**)&MX1, 64, nvar_l*nv_i*sizeof(REALTYPE));
    REALTYPE *MX2; hbw_posix_memalign((void**)&MX2, 64, nvar_l*nv_i*sizeof(REALTYPE));
    REALTYPE *CT; hbw_posix_memalign((void**)&CT, 64, nvar_l*nv_i*sizeof(REALTYPE));
#else
    REALTYPE *MT  = (REALTYPE *)_mm_malloc(RB*(nvar_l+nv_i)*sizeof(REALTYPE), 64);
    REALTYPE *NM = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
    REALTYPE *MX1 = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
    REALTYPE *MX2 = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
    REALTYPE *CT = (REALTYPE *)_mm_malloc(nvar_l*nv_i*sizeof(REALTYPE), 64);
#endif
    bzero((void*)NM, nvar_l*nv_i*sizeof(REALTYPE));
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
			if (erow!=nrow_i) {
				Trans_manual(MTT+IVar1, MT2, nvar_i);
			}
			else {
				for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
					REALTYPE MTemp = MTT[j*nvar_i+IVar1];
					if (MTemp!=MTemp || fabs(MTemp-65.0)<=0.001) {
						MT2[j] = 65.0;
					}
					else {
						MT2[j] = MTemp;
					}
				}
			}
            MT2+=RB;
        }
        for (IVar1=evar; IVar1<evar+nv_i; IVar1++) {
            for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
                  REALTYPE MTemp = MTT[j*nvar_i+IVar1%nvar_i];
				  if (MTemp!=MTemp || fabs(MTemp-65.0)<=0.001) {
					  MT2[j] = 65.0;
				  }
				  else {
					  MT2[j] = MTemp;
				  }
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
				if (erow!=nrow_i) {
					COV_manual(MT1,MT2,MX1+k,MX2+k,CT+k,NM+k);
				}
				else {
					REALTYPE NumMissing = 0.0;
					REALTYPE MeanX1 = 0.0;
					REALTYPE MeanX2 = 0.0;
					REALTYPE C = 0.0;
#pragma vector aligned
					for (IVar3=srow,j=0; IVar3<erow; IVar3++,j++) {
						REALTYPE M1 = MT1[j];
						REALTYPE M2 = MT2[j];
						if (M1!=65.0 && M2!=65.0) {
							NumMissing+=1.0;
							MeanX1+=M1;
							MeanX2+=M2;
							C+=M1*M2;
						}
					}
					NM[k]+=NumMissing;
					MX1[k]+=MeanX1;
					MX2[k]+=MeanX2;
					CT[k]+=C;
				}
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
            REALTYPE IRep = NM[j]<=0.0 ? 0.0 : 1.0/NM[j];
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

