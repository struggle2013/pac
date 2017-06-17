    SUBROUTINE FindCOV(NVAR, NROW, NV, MType, COV)

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: NVAR              !Number of variables
    INTEGER, INTENT(IN)    :: NROW              !Number of rows in MType
    INTEGER, INTENT(IN)    :: NV                !Number of cals in MType
    REAL,    INTENT(IN)    :: MType(NVAR,NROW) !Data used to calculate covariance
    REAL,    INTENT(INOUT) :: COV(NVAR,NV)   !Covariance matrix, output

    INTEGER :: IVar1,IVar2,IVar3
    INTEGER :: I2,VR
    INTEGER :: NumMissing(NV)
    REAL :: M1(NV)
    REAL :: MeanX1(NV), MeanX2(NV), C(NV), IRep

!$OMP PARALLEL DO private(IVar1,IVar2,IVar3,VR,I2,NUmMissing,MeanX1,MeanX2,C,M1) shared(NVAR,NV,NROW,MType) 
    DO IVar1=1, NVAR

    NumMissing = 0
    MeanX1 = 0.0
    MeanX2 = 0.0
    C = 0.0

    DO IVar3=1, NROW
        M1 = 0.0
        VR=IVar3*NVAR-NVAR
        IF(MType(IVar1,IVar3).ne.MType(IVar1,IVar3) .OR. ABS(MType(IVar1,IVar3)-65.0).lt.1.0e-3) THEN
        ELSEIF (IVar1+NV .le. NVAR) THEN
!DIR$ IVDEP            
            DO IVar2=1, NV
                I2=IVar1+IVar2
                IF(MType(I2,IVar3).eq.MType(I2,IVar3) .AND. ABS(MType(I2,IVar3)-65.0).ge.1.0e-3) THEN
                    M1(IVar2) = MType(I2,IVar3)
                ENDIF
            ENDDO
            DO IVar2=1, NV
                NumMissing(IVar2)=NumMissing(IVar2)+1
                MeanX1(IVar2)=MeanX1(IVar2)+MType(IVar1,IVar3)
                MeanX2(IVar2)=MeanX2(IVar2)+M1(IVar2)

                C(IVar2)=C(IVar2)+MType(IVar1,IVar3)*M1(IVar2)
            ENDDO
        ELSE
            DO IVar2=1, NV
                I2=MOD(IVar1+IVar2-1,NVAR)+1
                IF(MType(I2,IVar3).eq.MType(I2,IVar3) .AND. ABS(MType(I2,IVar3)-65.0).ge.1.0e-3) THEN
                    NumMissing(IVar2)=NumMissing(IVar2)+1
                    MeanX1(IVar2)=MeanX1(IVar2)+MType(IVar1,IVar3)
                    MeanX2(IVar2)=MeanX2(IVar2)+MType(I2,IVar3)

                    C(IVar2)=C(IVar2)+MType(IVar1,IVar3)*MType(I2,IVar3)
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    DO IVar2=1, NV
        VR=(IVar2-1)*NV+IVar1
        I2 = NumMissing(IVar2)
        IF(I2 .LE. 0) THEN
            COV(IVar1,IVar2)=0.0
        ELSE
            IRep=1.0/I2
            MeanX1(IVar2)=MeanX1(IVar2)*IRep
            MeanX2(IVar2)=MeanX2(IVar2)*IRep

            COV(IVar1,IVar2)=C(IVar2)*IRep-MeanX1(IVar2)*MeanX2(IVar2)
        ENDIF
    ENDDO
    ENDDO
!$OMP END PARALLEL DO
    RETURN

    END SUBROUTINE FindCOV


