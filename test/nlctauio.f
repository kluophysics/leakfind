C*==nlctauio.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NLCTAUIO(IFILTAU,WRTAUMQ,IPRINT,IE,ERYD,TAUQ,MSSQ,
     &                    OMEGAHAT,MUEHAT,NQ,NKMQ,CPACHNG,ICPAFLAG,
     &                    IQCPA,NDIMCLU,NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  read/write the   NLCPA   TAUQ, MSSQ and OMEGAHAT  matrices      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--NLCTAUIO11
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD
      INTEGER ICPAFLAG,IE,IFILTAU,IPRINT,IQCPA,NDIMCLU,NKMMAX,NQ,NQMAX
      LOGICAL WRTAUMQ
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MUEHAT(NDIMCLU,NDIMCLU),
     &           OMEGAHAT(NDIMCLU,NDIMCLU),TAUQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER NKMQ(NQMAX)
C
C Local variables
C
      CHARACTER*40 FMTE
      INTEGER I,IP,IQ,IQP,J,JP,M,N
      COMPLEX*16 MSSINP,RMSS,RTAU,TAUINP
      CHARACTER*5 STR5
C
C*** End of declarations rewritten by SPAG
C
      DATA FMTE/'(2F21.15,25X,I5,9X,I2,F15.6)'/
C
      IF ( WRTAUMQ ) THEN
C=======================================================================
C
         DO IQ = 1,NQ
C
            WRITE (IFILTAU,99001) ERYD,IQ,ICPAFLAG,CPACHNG
C
            N = NKMQ(IQ)
            DO I = 1,N
               DO J = 1,N
                  IF ( I.EQ.J ) THEN
                     WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),MSSQ(I,J,IQ)
                  ELSE
                     RTAU = TAUQ(I,J,IQ)/TAUQ(I,I,IQ)
                     RMSS = MSSQ(I,J,IQ)/MSSQ(I,I,IQ)
                     IF ( (CDABS(RTAU).GT.TOL) .OR. (CDABS(RMSS).GT.TOL)
     &                    ) WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),
     &                             MSSQ(I,J,IQ)
                  END IF
C
               END DO
            END DO
C
            IF ( IQ.EQ.IQCPA ) THEN
C
               WRITE (IFILTAU,99002) IQ
               DO J = 1,NDIMCLU
                  DO I = 1,NDIMCLU
                     WRITE (IFILTAU,99003) I,J,OMEGAHAT(I,J),MUEHAT(I,J)
                  END DO
               END DO
            END IF
C
         END DO
C
      ELSE
C=======================================================================
C
         M = NKMMAX*NKMMAX*NQMAX
         CALL CINIT(M,MSSQ)
         CALL CINIT(M,TAUQ)
C
         DO IQ = 1,NQ
C
            N = NKMQ(IQ)
C
 20         CONTINUE
            READ (IFILTAU,'(A5)',END=50) STR5
            IF ( STR5.NE.'*****' ) GOTO 20
C
            READ (IFILTAU,FMT=FMTE,END=50) ERYD,IQP
            IF ( IPRINT.GT.0 ) WRITE (6,99004) IE,ERYD,IQ,IQP
            IF ( IQ.NE.IQP ) GOTO 20
C
 40         CONTINUE
            READ (IFILTAU,99003) I,J,TAUINP,MSSINP
            IF ( (I+J).NE.0 ) THEN
               TAUQ(I,J,IQ) = TAUINP
               MSSQ(I,J,IQ) = MSSINP
               IF ( (I+J).LT.2*N ) GOTO 40
            END IF
C
            IF ( IQ.EQ.IQCPA ) THEN
               READ (IFILTAU,'(/)')
               DO J = 1,NDIMCLU
                  DO I = 1,NDIMCLU
                     READ (IFILTAU,99003) IP,JP,OMEGAHAT(I,J),
     &                     MUEHAT(I,J)
                     IF ( I.NE.IP .OR. J.NE.JP )
     &                     STOP '<NLCPATAUIO> >> OMEGA'
                  END DO
               END DO
            END IF
C
         END DO
         RETURN
C
 50      CONTINUE
         STOP 'EOF reached while reading NLCPA - TAU-file IFILTAU'
C
C=======================================================================
      END IF
C
99001 FORMAT (/,80('*')/,2F21.15,' RYD   TAU-C M-C  FOR IQ=',I5,:,
     &        '  CPA:',I2,F15.8)
99002 FORMAT (/,'OMEGA-HAT  and  MUE-HAT  FOR IQ=',I5)
99003 FORMAT (2I5,1P,4E22.14)
99004 FORMAT (/,10X,'energy   ',I10,2F10.5,'  IQ ',2I5)
      END
C*==nlcdumptau.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NLCDUMPTAU(IE,ERYD,IWRI,MUEHAT,OMEGAHAT,TAUHAT,IREL,N,
     &                      NDIMCLU)
C   ********************************************************************
C   *                                                                  *
C   *  dump the NLCPA matrices   MUEHAT, OMEGAHAT, TAUHAT              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--NLCDUMPTAU154
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IE,IREL,IWRI,N,NDIMCLU
      COMPLEX*16 MUEHAT(NDIMCLU,NDIMCLU),OMEGAHAT(NDIMCLU,NDIMCLU),
     &           TAUHAT(NDIMCLU,NDIMCLU)
C
C*** End of declarations rewritten by SPAG
C
      WRITE (IWRI,99001) IE,ERYD
C
      WRITE (IWRI,99002) 'NLCPA - matrices'
C
      CALL CMATSTRUCT('MUEHAT  ',MUEHAT,NDIMCLU,NDIMCLU,0,0,N,TOL,IWRI)
      CALL CMATSTRUCT('OMEGAHAT',OMEGAHAT,NDIMCLU,NDIMCLU,0,0,N,TOL,
     &                IWRI)
      CALL CMATSTRUCT('TAUHAT  ',TAUHAT,NDIMCLU,NDIMCLU,0,0,N,TOL,IWRI)
      IF ( IREL.LT.0 ) STOP 'IREL'
C
99001 FORMAT (/,10X,'IE =',I4,2X,' E =',2F12.6,'  RYD',/)
99002 FORMAT (//,1X,79('#'),/,10X,A,/,1X,79('#'),/)
      END
