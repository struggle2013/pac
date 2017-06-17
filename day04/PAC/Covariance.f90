    REAL*8 FUNCTION Covariance(NObs, X1, X2)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: NObs
        REAL*8,    INTENT(IN) :: X1(NObs)
        REAL*8,    INTENT(IN) :: X2(NObs)


        INTEGER :: IObs, NumMissing
        REAL*8 :: MeanX1, MeanX2
        LOGICAL,EXTERNAL :: IsMissingPheno
        NumMissing=0
        MeanX1=0.0
        MeanX2=0.0
        Covariance=0.0 
        
        !$omp parallel do REDUCTION (+:NumMissing,MeanX2,MeanX1)
        DO IObs=1, NObs
            IF(IsMissingPheno(X1(IObs)) .OR. &
                IsMissingPheno(X2(IObs))) THEN
                NumMissing=NumMissing+1
       
            ELSE
                MeanX1=MeanX1+X1(IObs)
                MeanX2=MeanX2+X2(IObs)

                Covariance=Covariance+X1(IObs)*X2(IObs)
     
            ENDIF
            
        ENDDO
        !$omp end parallel do
        MeanX1=MeanX1/(NObs-NumMissing)
        MeanX2=MeanX2/(NObs-NumMissing)

        IF(NObs-NumMissing .LE. 0) THEN
            Covariance=0.0
        ELSE
            Covariance=Covariance/(NObs-NumMissing)-MeanX1*MeanX2
   
        ENDIF


        RETURN

    END FUNCTION Covariance
