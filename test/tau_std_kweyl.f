C*==tau_std_kweyl.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_STD_KWEYL(ERYD,P,IPRINT,ICPAFLAG,ITCPA,ICPACONV,
     &                         CPACHNG,TSST,MSST,TSSQ,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *    perform the loop over all k- vectors of a  Weyl-mesh          *
C   *                                                                  *
C   *    - set up:   G-matrix  (structure constants)                   *
C   *    - set up:   M          (KKR-MATRIX)                           *
C   *    - invert    M    to get  TAU(K)                               *
C   *    - update    TAU                                               *
C   *    - do the CPA if requested                                     *
C   *                                                                  *
C   * 12/08/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:DROT,IWEDGEROT,MROTK,NSYMACCEPTED,NSYM,
     &    SYMACCEPTED,SYMUNITARY,IQORGQP,NWEDGE
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,DROTQ
      USE MOD_KSPACE,ONLY:NKTAB,NPTMAX,NPTMIN
      USE MOD_ENERGY,ONLY:IMEMIN,IMEMAX
      USE MOD_LATTICE,ONLY:BBAS
      IMPLICIT NONE
C*--TAU_STD_KWEYL28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_STD_KWEYL')
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 A,CPACORR,CPAERR,CPAERRL,F1,F2,F3,KVEC(3),KVEC0(3),WKSUM,
     &       ZETAX,ZETAY,ZETAZ
      INTEGER I,I1,IA_ERR,IK,INFO,IPIV(:),IQ,IROT,IWEDGE,J,J1,N
      COMPLEX*16 MAUX(:,:),SUMQ(:,:,:),TAUK(:,:),W1(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,IPIV,MAUX,SUMQ,W1
C
      ALLOCATE (SUMQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (IPIV(NKKR),MAUX(NKKR,NKKR))
      ALLOCATE (TAUK(NKKR,NKKR),W1(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SUMQ')
C
C ------------------------------------- initialize set up of weyl-k-mesh
C
      A = 0.359D0
      ZETAX = DSQRT(3.0D0)*A
      ZETAY = DSQRT(5.0D0)*A
      ZETAZ = 2.0D0*DSQRT(13.0D0)*A
C
      IF ( (IMEMAX-IMEMIN).GT.0D0 ) THEN
         NKTAB = NPTMIN + 
     &           INT((NPTMAX-NPTMIN)*((IMEMAX-DIMAG(ERYD))/(IMEMAX-
     &           IMEMIN))**2)
      ELSE
         NKTAB = NPTMAX
      END IF
C
      IF ( IPRINT.GE.1 ) WRITE (6,*) '    NKTAB= ',NKTAB
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ICPACONV = 0
      CPAERRL = 1.0D+6
      ITCPA = 0
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
      CALL CINIT(NKMMAX*NKMMAX*NQ,SUMQ)
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
C --------------------------------------- create ->k in  units of  2PI/A
C -----------------------------------------------------  -0.5 < F  < 0.5
C
         F1 = IK*ZETAX - IDINT(IK*ZETAX) - 0.5D0
         F2 = IK*ZETAY - IDINT(IK*ZETAY) - 0.5D0
         F3 = IK*ZETAZ - IDINT(IK*ZETAZ) - 0.5D0
C
         DO I = 1,3
            KVEC0(I) = F1*BBAS(I,1) + F2*BBAS(I,2) + F3*BBAS(I,3)
         END DO
C
C WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
         DO IWEDGE = 1,NWEDGE
C
            IROT = IWEDGEROT(IWEDGE)
C
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KVEC0,1,0D0,KVEC,1)
C
            CALL STRSET(IK,KVEC,MAUX,TAUK,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,NKMMAX)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C------------------------------------------------------------ store TAUQ
            DO IQ = 1,NQ
               I1 = IND0Q(IQ) + 1
               N = NKMQ(IQ)
               DO J = 1,N
                  J1 = IND0Q(IQ) + J
                  CALL ZAXPY(N,C1,TAUK(I1,J1),1,SUMQ(1,J,IQ),1)
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END DO
C
C WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      END DO
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      WKSUM = NKTAB*NWEDGE
C
      CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,WKSUM,SUMQ,TAUQ,W1,NQ,NKMQ,
     &               DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM,
     &               NSYMACCEPTED,NQMAX,NKMMAX)
C
      IF ( NCPA.GT.0 ) THEN
C
         IF ( ICPAALG.EQ.1 ) THEN
            CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                    NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,NTMAX,
     &                    NQMAX,NKMMAX)
C
         ELSE
            CALL CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                     NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,KMROT,
     &                     DROTQ,NTMAX,NQMAX,NKMMAX)
         END IF
C
         CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,SUMQ,MSSQ,W1,NQ,NKMQ,
     &                  DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                  NSYMACCEPTED,NQMAX,NKMMAX)
C
         IF ( IPRINT.GE.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
         IF ( CPAERR.LE.CPATOL ) THEN
            ICPACONV = 1
            IF ( IPRINT.GE.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &                                CPACHNG
         ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
            WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 1
         ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
            WRITE (6,99003) ITCPA
            WRITE (6,99004) CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 2
         ELSE
            CPAERRL = CPAERR
            GOTO 100
         END IF
C
      END IF
C
      DEALLOCATE (TAUK,IPIV,MAUX,SUMQ,W1,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
      END
