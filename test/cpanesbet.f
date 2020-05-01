C*==cpanesbet.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                     NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,KMROT,
     &                     DROTQ,NTMAX,NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   PERFORM  CPA-ITERATION ACCORDING    NEBETS'S  ALGORITHM        *
C   *                                                                  *
C   *   t[c, i+1] = SUM(a)  c(a) * {   t[a] +                          *
C   *     ( t[c,i] * ( t[a] - (t[a] - t[c,i])*B[c,i] )**(-1)   -1 )    *
C   *                                           * (t[a] - t[c,i])   }  *
C   *                                                                  *
C   *   THE CPA - ITERATION STEP FOR SITE IQ IS OMITTED IF             *
C   *   ICPA(IQ) = 0 ( SET IN <INITALL> )                              *
C   *   MSSQ contains  t[c]**(-1)                                      *
C   *   TAU    contains    BzInt  [ t**(-1) - G(k) ]**(-1)  dk         *
C   *                                                                  *
C   * 04/01/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CPANESBET23
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CPANESBET')
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPACORR,CPAERR
      INTEGER IPRINT,KMROT,NKMMAX,NQ,NQMAX,NTMAX
      REAL*8 CONC(NTMAX)
      COMPLEX*16 DROTQ(NKMMAX,NKMMAX,NQMAX),MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      INTEGER ICPA(NQMAX),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)
C
C Local variables
C
      COMPLEX*16 B(:,:),CSUM,DMAMC(:,:),ERR(:,:),W1(:,:),W2(:,:)
      INTEGER I,IA_ERR,IO,IQ,IT,J,M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE B,DMAMC,W1,W2,ERR
C
      ALLOCATE (ERR(NKMMAX,NKMMAX))
      ALLOCATE (B(NKMMAX,NKMMAX),DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TSSQ')
      CPAERR = 0.0D0
      CPACORR = 0.0D0
      CPACHNG = 0.0D0
C
C================================================================= IQ ==
      DO IQ = 1,NQ
C
         M = NKMMAX
         N = NKMQ(IQ)
C
C     --------------------------------------------------------
C      convert  TAU  to  B[c]  :   B[c]  =  t[c]**(-1) * TAU
C     --------------------------------------------------------
C
         CALL CMATMUL(N,M,MSSQ(1,1,IQ),TAUQ(1,1,IQ),B)
C
C     --------------------------------------------------------
C      invert MSSQ  :   TSSQ  =  MSSQ**(-1)
C     --------------------------------------------------------
C
         CALL CMATCOP(N,M,MSSQ(1,1,IQ),W1)
C
         CALL CMATINV(N,M,W1,TSSQ(1,1,IQ))
C
C     --------------------------------------------------------
C      invert MSST     :      TSST  =  MSST**(-1)
C     --------------------------------------------------------
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
C
C ------------------------- rotate the single site m-matrix if necessary
            IF ( KMROT.NE.0 ) THEN
C
               CALL ROTATE(MSST(1,1,IT),'L->G',W1,N,DROTQ(1,1,IQ),M)
C
            ELSE
C
               CALL CMATCOP(N,M,MSST(1,1,IT),W1)
C
            END IF
C
         END DO
C
C
         IF ( ICPA(IQ).NE.0 ) THEN
C
            DO J = 1,N
               DO I = 1,N
                  ERR(I,J) = C0
               END DO
            END DO
C
C
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
C     -------------------------------------------
C                   ( t[a] - t[c] )
C     -------------------------------------------
               CALL GETDELM(DMAMC,IT,IQ,NKMQ,TSSQ,TSST,NTMAX,NQMAX,M)
C
C     -------------------------------------------
C                   ( t[a] - t[c] ) * B[c]
C     -------------------------------------------
               CALL CMATMUL(N,M,DMAMC,B,W1)
C
C     -------------------------------------------
C               t[a]  -  ( t[a] - t[c] ) * B[c]
C     -------------------------------------------
               DO J = 1,N
                  DO I = 1,N
                     W1(I,J) = TSST(I,J,IT) - W1(I,J)
                  END DO
               END DO
C
C
C     -------------------------------------------
C      ( t[a]  -  ( t[a] - t[c] ) * B[c] )**(-1)
C     -------------------------------------------
               CALL CMATINV(N,M,W1,W2)
C
C     --------------------------------------------
C     t[c] * ( t[a] - ( t[a] - t[c] )*B[c] )**(-1)
C     --------------------------------------------
               CALL CMATMUL(N,M,TSSQ(1,1,IQ),W2,W1)
C
C     ----------------------------------------------------
C     t[c] * ( t[a] - ( t[a] - t[c] )*B[c] )**(-1)    - 1
C     ----------------------------------------------------
               DO I = 1,N
                  W1(I,I) = W1(I,I) - 1D0
               END DO
C
C   --------------------------------------------------------------------
C   ( t[c] * ( t[a] - (t[a] - t[c])*B[c] )**(-1)   - 1 ) * (t[a] - t[c])
C   --------------------------------------------------------------------
               CALL CMATMUL(N,M,W1,DMAMC,W2)
C
C     -------------------------------------------------------------
C     add  t[a] , then multiply with  c[a]  and at last sum over  a
C     -------------------------------------------------------------
               DO I = 1,N
                  DO J = 1,N
                     ERR(I,J) = ERR(I,J) + CONC(IT)
     &                          *(W2(I,J)+TSST(I,J,IT))
                  END DO
               END DO
C
C
            END DO
C================================================================= IT ==
C
C
            DO I = 1,N
               CPACHNG = DBLE
     &                   (MAX(CPACHNG,ABS((ERR(I,I)-TSSQ(I,I,IQ))/TSSQ
     &                   (I,I,IQ))))
            END DO
C
C
            DO I = 1,N
               DO J = 1,N
C-----------------------------------------------------------------------
C  t[c, i+1] :=  p * t[c, i+1]  +  ( 1 - p ) * t[c, i] ,  mit  p = 0...1
C-----------------------------------------------------------------------
                  TSSQ(I,J,IQ) = 1.0D0*ERR(I,J) + 0.0D0*TSSQ(I,J,IQ)
               END DO
            END DO
C
C-----------------------------------------------------------------------
C   t[c, i+1]  should be a symmetric matrix for all i,
C      so symmetrize  t[c, i+1]  explicitly
C-----------------------------------------------------------------------
            DO I = 1,N
               DO J = I + 1,N
                  TSSQ(I,J,IQ) = 0.5D0*(TSSQ(I,J,IQ)+TSSQ(J,I,IQ))
                  TSSQ(J,I,IQ) = TSSQ(I,J,IQ)
               END DO
            END DO
C
C     --------------------------------------------------------
C      invert TSSQ  :   MSSQ  =  TSSQ**(-1)
C     --------------------------------------------------------
            CALL CMATINV(N,M,TSSQ(1,1,IQ),MSSQ(1,1,IQ))
C
            CPAERR = CPACHNG
C
            IF ( IPRINT.GE.2 ) THEN
               CSUM = C0
               DO I = 1,N
                  CSUM = CSUM + MSSQ(I,I,IQ)
               END DO
               WRITE (6,99001) IQ,CPAERR,CPACORR,CSUM
            END IF
C
         END IF
C
      END DO
C================================================================= IQ ==
C
      DEALLOCATE (B,DMAMC,W1,W2,ERR,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (' CPA:  IQ',I3,'  ERR',F12.5,'  CORR',F13.5,'  M',
     &        18(1X,2(1PE14.6)))
      END
