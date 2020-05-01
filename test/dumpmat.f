C*==dumpmat.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DUMPMAT(IE,ERYD,IWRI,MSST,MSSQ,TAUT,TAUQ,KME,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *  dump the matrices   MSST, MSSQ, TAUT and TAUQ                   *
C   *                                                                  *
C   *  KME = .T.           in addition: MEZZ, MEZJ for DOS and SPIN    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NMEMAX
      USE MOD_SITES,ONLY:NQ,NQMAX
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--DUMPMAT16
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
      INTEGER IE,IWRI
      LOGICAL KME
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER IQ,IT,M,N
C
C*** End of declarations rewritten by SPAG
C
      WRITE (IWRI,99001) IE,ERYD
C
      N = NKM
      M = NKMMAX
C
      DO IT = ITBOT,ITTOP
         WRITE (IWRI,99002) 'atom type   ',IT
         CALL CMATSTRUCT('Mss[t]',MSST(1,1,IT),N,M,IREL,IREL,1,TOL,IWRI)
         CALL CMATSTRUCT('TAU[t]',TAUT(1,1,IT),N,M,IREL,IREL,1,TOL,IWRI)
         IF ( KME ) THEN
            CALL CMATSTRUCT('MEZZ[t] DOS ',MEZZ(1,1,IT,1),N,M,IREL,IREL,
     &                      1,TOL,IWRI)
            CALL CMATSTRUCT('MEZJ[t] DOS ',MEZJ(1,1,IT,1),N,M,IREL,IREL,
     &                      1,TOL,IWRI)
            CALL CMATSTRUCT('MEZZ[t] SPIN',MEZZ(1,1,IT,2),N,M,IREL,IREL,
     &                      1,TOL,IWRI)
            CALL CMATSTRUCT('MEZJ[t] SPIN',MEZJ(1,1,IT,2),N,M,IREL,IREL,
     &                      1,TOL,IWRI)
         END IF
      END DO
C
      DO IQ = 1,NQ
         WRITE (IWRI,99002) 'lattice site ',IQ
         CALL CMATSTRUCT('Mss[q]',MSSQ(1,1,IQ),N,M,IREL,IREL,1,TOL,IWRI)
         CALL CMATSTRUCT('TAU[q]',TAUQ(1,1,IQ),N,M,IREL,IREL,1,TOL,IWRI)
      END DO
C
99001 FORMAT (/,10X,'IE =',I4,2X,' E =',2F12.6,'  RYD',/)
99002 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),/)
      END
