C*==getdelm.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GETDELM(DM,IT,IQ,NKMQ,MSSQ,MSS,NTMAX,NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the difference of the single site m - matrices       *
C   *                                                                  *
C   *                     DM = m[a]-m[c]                               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--GETDELM12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,IT,NKMMAX,NQMAX,NTMAX
      COMPLEX*16 DM(NKMMAX,NKMMAX),MSS(NKMMAX,NKMMAX,NTMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER NKMQ(NQMAX)
C
C Local variables
C
      INTEGER I1,I2
C
C*** End of declarations rewritten by SPAG
C
      DO I2 = 1,NKMQ(IQ)
         DO I1 = 1,NKMQ(IQ)
            DM(I1,I2) = MSS(I1,I2,IT) - MSSQ(I1,I2,IQ)
         END DO
      END DO
C
      END
