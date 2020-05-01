C*==setkkr.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SETKKR(NQ,NKMQ,IND0Q,TAUKINV,MSSQ,NQMAX,NKKR,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SET UP THE KKR-MATRIX                                          *
C   *                                                                  *
C   *           M  =  [ m  - G ]                                       *
C   *           =       =    =                                         *
C   *                                                                  *
C   *   NOTE:  TAUKINV = -G  SETUP IN <STRSET>                         *
C   *          ADD ONLY  MSSQ                                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SETKKR16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKKR,NKMMAX,NQ,NQMAX
      INTEGER IND0Q(NQMAX),NKMQ(NQMAX)
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUKINV(NKKR,NKKR)
C
C Local variables
C
      INTEGER IKM1,IKM2,IQ,KM1,KM2
C
C*** End of declarations rewritten by SPAG
C
      DO IQ = 1,NQ
C
         IKM2 = IND0Q(IQ)
         DO KM2 = 1,NKMQ(IQ)
            IKM2 = IKM2 + 1
C
            IKM1 = IND0Q(IQ)
            DO KM1 = 1,NKMQ(IQ)
               IKM1 = IKM1 + 1
               TAUKINV(IKM1,IKM2) = TAUKINV(IKM1,IKM2)
     &                              + MSSQ(KM1,KM2,IQ)
            END DO
C
         END DO
      END DO
C
      END
