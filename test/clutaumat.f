C*==clutaumat.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUTAUMAT(KIJ,IPRINT,IREL,NKKR,NLM,IQCLU1,IQCLU2,
     &                     NLMQCLU,NKMQCLU,CONC,WA,WB,WC,WD,IPIV,
     &                     ITCPAMAX,NOQCLU,ITOQCLU,CPATOL,NCPA,ICPACLU,
     &                     CPACHNG,ICPAFLAG,NSYM,NSYMACCEPTED,DROT,
     &                     IQCLUORG,SYMACCEPTED,SYMUNITARY,SYMMETRIZE,
     &                     NKMMAX,NQCLU,NTMAX,TSSTCLU,MSSTCLU,TSSQCLU,
     &                     MSSQCLU,TAUQCLU)
C   ********************************************************************
C   *                                                                  *
C   *   calculate cluster TAU-matrix                                   *
C   *                                                                  *
C   *   on entry:    WA = -G0                                          *
C   *                                                                  *
C   *   - invert TAUINV = m - G                                        *
C   *   - cut out diagonal elements   TAUQCLU                          *
C   *   - perform CPA - cycle if requested                             *
C   *                                                                  *
C   *   KIJ = 0:  ONLY diagonal elements  TAUQCLU  are calculated      *
C   *                                                                  *
C   *   KIJ = 1:  off-diagonal elements are calculated in addition     *
C   *                                                                  *
C   *             on exit:       for IREL =      1   2   3             *
C   *                                                                  *
C   *                       WA = TAUIJ(dn,dn)    *   *   *             *
C   *                       WB = TAUIJ(up,up)    -   *   *             *
C   *                       WC = TAUIJ(up,dn)    -   -   *             *
C   *                       WD = TAUIJ(dn,up)    -   -   *             *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--CLUTAUMAT35
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUTAUMAT')
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPATOL
      INTEGER ICPAFLAG,IPRINT,IQCLU1,IQCLU2,IREL,ITCPAMAX,KIJ,NCPA,NKKR,
     &        NKMMAX,NLM,NQCLU,NSYM,NSYMACCEPTED,NTMAX
      LOGICAL SYMMETRIZE
      REAL*8 CONC(NTMAX)
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),
     &           MSSQCLU(NKMMAX,NKMMAX,NQCLU),
     &           MSSTCLU(NKMMAX,NKMMAX,NTMAX),
     &           TAUQCLU(NKMMAX,NKMMAX,NQCLU),
     &           TSSQCLU(NKMMAX,NKMMAX,NQCLU),
     &           TSSTCLU(NKMMAX,NKMMAX,NTMAX),WA(NKKR,NKKR),
     &           WB(NKKR,NKKR),WC(NKKR,NKKR),WD(NKKR,NKKR)
      INTEGER ICPACLU(NQCLU),IPIV(NKKR),IQCLUORG(NSYMMAX,NQCLU),
     &        ITOQCLU(NTMAX,NQCLU),NKMQCLU(NQCLU),NLMQCLU(NQCLU),
     &        NOQCLU(NQCLU)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      REAL*8 CPACORR,CPAERR,CPAERRL
      INTEGER IA_ERR,IQCLU,ITCPA,NQCLU_CPA
      COMPLEX*16 MG0MAT(:,:),W1(:,:),WTAU(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MG0MAT,WTAU,W1
C
      SYMMETRIZE = .FALSE.
C
C----------------------------- store G-matrix in case of CPA calculation
C
      IF ( NCPA.NE.0 ) THEN
C
         ALLOCATE (MG0MAT(NKKR,NKKR),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MG0MAT')
C
         CALL ZCOPY(NKKR*NKKR,WA,1,MG0MAT,1)
      END IF
C
      IF ( SYMMETRIZE ) THEN
C
         ALLOCATE (WTAU(NKMMAX,NKMMAX,NQCLU),W1(NKMMAX,NKMMAX),
     &             STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WTAU')
C
      END IF
C
      CPAERRL = 1.0D+6
      ITCPA = 0
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      IF ( ITCPA.GT.1 ) CALL ZCOPY(NKKR*NKKR,MG0MAT,1,WA,1)
C
C=======================================================================
C         non-/scalar-relativistic AND non-spin-polarized calculation
C=======================================================================
      IF ( IREL.LE.1 ) THEN
C
         CALL CINVRIJ1(KIJ,WA,WC,MSSQCLU,TAUQCLU,IPIV,NQCLU,NLMQCLU,
     &                 NKKR,NKMMAX)
C
C=======================================================================
C         non-/scalar-relativistic AND spin-polarized calculation
C=======================================================================
      ELSE IF ( IREL.EQ.2 ) THEN
C
         CALL CINVRIJ2(KIJ,WA,WB,WC,MSSQCLU,TAUQCLU,IPIV,NQCLU,NLM,
     &                 NLMQCLU,NKKR,NKMMAX)
C
C=======================================================================
C             fully relativistic / spin-polarized calculation
C=======================================================================
      ELSE IF ( IREL.GE.3 ) THEN
C
         IF ( KIJ.EQ.0 ) THEN
            CALL CINVRII3(WA,WB,WC,WD,MSSQCLU,TAUQCLU,IPIV,IQCLU1,
     &                    IQCLU2,NQCLU,NLM,NLMQCLU,NKKR,NKMMAX)
         ELSE
            CALL CINVRIJ3(WA,WB,WC,WD,MSSQCLU,TAUQCLU,IPIV,NQCLU,NLM,
     &                    NLMQCLU,NKKR,NKMMAX)
         END IF
C
      END IF
C
      IF ( SYMMETRIZE ) CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,WTAU,
     &                                 TAUQCLU,W1,NQCLU,NKMQCLU,DROT,
     &                                 IQCLUORG,SYMUNITARY,SYMACCEPTED,
     &                                 NSYM,NSYMACCEPTED,NQCLU,NKMMAX)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc CPA
      IF ( NCPA.GT.0 ) THEN
C
C---------------------- perform CPA only for central site of the cluster
         NQCLU_CPA = 1
C
C         write(*,*) '**********************************  ITCPA', ITCPA
C         CALL  CMATSTR('TAU  ',TAUQCLU(1,1,1),18,18,2,2,1,1.D-8,6)
C         CALL  CMATSTR('TSST VOR ',TSSTCLU(1,1,1),18,18,2,2,1,1.D-8,6)
C         CALL  CMATSTR('MSSQ VOR ',MSSQCLU(1,1,1),18,18,2,2,1,1.D-8,6)
C         write(*,*) ' NKMQCLU', NKMQCLU
C         write(*,*) ' NOQCLU ', NOQCLU
C         write(*,*) ' ITOQCLU', ITOQCLU
C         write(*,*) ' CONC   ', CONC
C         write(*,*) ' NQCLU  ', NQCLU
C         write(*,*) ' NTMAX  ', NTMAX
C
         CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPACLU,NQCLU_CPA,
     &                 NKMQCLU,NOQCLU,ITOQCLU,CONC,TSSTCLU,MSSTCLU,
     &                 TSSQCLU,MSSQCLU,TAUQCLU,NTMAX,NQCLU,NKMMAX)
C         CALL  CMATSTR('MSSQ NACH',MSSQCLU(1,1,1),18,18,2,2,1,1.D-8,6)
C         stop
C
         DO IQCLU = 2,NQCLU
            MSSQCLU(:,:,IQCLU) = MSSQCLU(:,:,1)
         END DO
C
         IF ( SYMMETRIZE ) CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,
     &        WTAU,MSSQCLU,W1,NQCLU,NKMQCLU,DROT,IQCLUORG,SYMUNITARY,
     &        SYMACCEPTED,NSYM,NSYMACCEPTED,NQCLU,NKMMAX)
C
         IF ( IPRINT.EQ.1 ) WRITE (6,99001) CPAERR,CPACORR,CPACHNG
C
         IF ( CPAERR.LE.CPATOL ) THEN
            ICPAFLAG = 0
            WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
         ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
            WRITE (6,99003) ITCPA,CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 1
         ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
            WRITE (6,99004) ITCPA
            WRITE (6,99001) CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 2
         ELSE
            CPAERRL = CPAERR
            GOTO 100
         END IF
C
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc CPA
      IF ( ALLOCATED(MG0MAT) ) DEALLOCATE (MG0MAT)
      IF ( ALLOCATED(WTAU) ) DEALLOCATE (WTAU,W1)
C      stop
C
99001 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
99002 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99003 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99004 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        I4,' >>>> iteration stopped ',5X,10('!'))
C
      END
