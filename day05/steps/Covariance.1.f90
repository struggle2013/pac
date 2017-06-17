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
    REAL :: MeanX1(NV), MeanX2(NV), C(NV), IRep
    LOGICAL,EXTERNAL :: IsMissingPheno

!$OMP PARALLEL DO private(IVar1,IVar2,IVar3,VR,I2,NUmMissing,MeanX1,MeanX2,C) shared(NVAR,NV,NROW,MType) 
    DO IVar1=1, NVAR

    NumMissing = NROW
    MeanX1 = 0.0
    MeanX2 = 0.0
    C = 0.0

    DO IVar3=1, NROW
        VR=IVar3*NVAR-NVAR
        IF(IsMissingPheno(MType(IVar1,IVar3))) THEN
            DO IVar2=1, NV
                NumMissing(IVar2)=NumMissing(IVar2)-1
            ENDDO
        ELSEIF (IVar1+NV .le. NVAR) THEN
            DO IVar2=1, NV
                I2=IVar1+IVar2
                IF(IsMissingPheno(MType(I2,IVar3))) THEN
                    NumMissing(IVar2)=NumMissing(IVar2)-1
                ELSE
                    MeanX1(IVar2)=MeanX1(IVar2)+MType(IVar1,IVar3)
                    MeanX2(IVar2)=MeanX2(IVar2)+MType(I2,IVar3)

                    C(IVar2)=C(IVar2)+MType(IVar1,IVar3)*MType(I2,IVar3)
                ENDIF
            ENDDO
        ELSE
            DO IVar2=1, NV
                I2=MOD(IVar1+IVar2-1,NVAR)+1
                IF(IsMissingPheno(MType(I2,IVar3))) THEN
                    NumMissing(IVar2)=NumMissing(IVar2)-1
                ELSE
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


