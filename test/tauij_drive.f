C*==tauij_drive.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_DRIVE(IECURR,ERYD,P,TSSQ,MSSQ,TSST,MSST,TAUQ,
     &                       ICPAFLAG,CPACHNG,ITCPA,ICPACONV)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine to calculate the site off-diagonal elements      *
C   *  of the scattering path operator   TAU_ij                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NQMAX
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_FILES,ONLY:WRTAUIJ,IPRINT
      USE MOD_CALCMODE,ONLY:KKRMODE
      USE MOD_KSPACE,ONLY:IBZINT
      IMPLICIT NONE
C*--TAUIJ_DRIVE18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_DRIVE')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IECURR,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C***********************************************************************
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
C=======================================================================
C                        free CLUSTER calculation
C=======================================================================
C
         IF ( IBZINT.EQ.0 ) THEN
C
Ccc            CALL CLUSTER(KTAUIJ,IPRINTBAND,ERYD,MSSQ,TAUQ)
            STOP 'CLUSTER calculation'
C
C=======================================================================
C               BZ - integration using  SPECIAL POINTS
C=======================================================================
C
         ELSE IF ( IBZINT.EQ.2 ) THEN
C
            CALL TAUIJ_STD(ERYD,P,TAUQ,MSSQ)
C
         ELSE
            WRITE (6,99001) IBZINT
            STOP 'in <TAUIJ_DRIVE>: IBZINT <> 0, 2'
         END IF
C
         IF ( WRTAUIJ ) CALL TAUIJ_HOST_INVWRI(IECURR,ERYD,TAUQ,MSSQ)
C
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
         CALL TAUIJ_TB(ERYD,P,MSSQ)
C
         IF ( WRTAUIJ ) CALL TAUIJ_HOST_INVWRI(IECURR,ERYD,TAUQ,MSSQ)
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
C***********************************************************************
C
         CALL TAUIJ_EMBEDDED(IPRINT,ERYD,ICPAFLAG,CPACHNG,ITCPA,
     &                       ICPACONV,IECURR,TSST,MSST,TSSQ,MSSQ,TAUQ)
C
C***********************************************************************
      ELSE
C
         STOP '<TAUIJ_DRIVE>: KKRMODE incorrect'
C
      END IF
C***********************************************************************
C
99001 FORMAT (//,1X,79('#'),/,10X,'IBZINT=',i4,/)
      END
