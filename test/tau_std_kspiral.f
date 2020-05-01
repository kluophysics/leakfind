C*==tau_std_kspiral.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_STD_KSPIRAL(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                           ICPACONV,TSST,MSST,TSSQ,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *    this is the special version to deal with SPIN SPIRALS         *
C   *    for KMROT=3 or 4 and implies working with the                 *
C   *    spin-polarized scalar - relativistic mode   IREL = 2          *
C   *    in the (l,ml)-representation for REAL spherical harmonics     *
C   *                                                                  *
C   * HE 06/05/09                                                      *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:WKTAB,NKTAB,KTAB
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,IQORGQP
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_SITES,ONLY:QMVEC,NQMAX,ITOQ,NOQ,ICPA,DROTQ,NQ
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--TAU_STD_KSPIRAL36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_STD_KSPIRAL')
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
      REAL*8 CPACORR,CPAERR,CPAERRL,KVEC(3),WK,WKSUM
      COMPLEX*16 CSCL,CWK,MAUX(:,:),MAUX2(:,:),SUMQ(:,:,:),TAUK(:,:),
     &           TAUK2(:,:),W1(:,:)
      INTEGER I,I1DN,I1UP,IA_ERR,II1DN,II1UP,IK,IQ,IQP,ISYM,J,JDN,JJ0DN,
     &        JJ0UP,JJDN,JJUP,JUP,N,NLMQ(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUK,MAUX2,TAUK2,NLMQ,SUMQ,W1
C
      ALLOCATE (W1(NKMMAX,NKMMAX))
      ALLOCATE (NLMQ(NQMAX),SUMQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MAUX2(NKKR*2,NKKR*2))
      ALLOCATE (TAUK2(NKKR*2,NKKR*2))
      ALLOCATE (MAUX(NKKR,NKKR))
      ALLOCATE (TAUK(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      DO IQ = 1,NQ
         NLMQ(IQ) = NKMQ(IQ)/2
      END DO
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
      CALL CINIT(NKMMAX*NKMMAX*NQ,SUMQ)
C
      WKSUM = 0D0
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
C--------------------------------------------- create  -G(up)  and store
C
         CALL CINIT(NKKR*NKKR,MAUX)
         CALL CINIT(NKKR*NKKR,TAUK)
         CALL CINIT(NKKR*NKKR*4,MAUX2)
         CALL CINIT(NKKR*NKKR*4,TAUK2)
C
         DO I = 1,3
            KVEC(I) = KTAB(I,IK) + 0.5D0*QMVEC(I)
         END DO
C
         CALL STRSET(IK,KVEC,TAUK,MAUX,P)
C
C------------------------------------------- copy  -G(up)  to KKR-matrix
C
         DO J = 1,NKKR
            CALL ZCOPY(NKKR,MAUX(1,J),1,MAUX2(1,J),1)
         END DO
C
C-------------------------------------------------------- create  -G(dn)
C
         DO I = 1,3
            KVEC(I) = KTAB(I,IK) - 0.5D0*QMVEC(I)
         END DO
C
         CALL STRSET(IK,KVEC,TAUK,MAUX,P)
C
C------------------------------------------- copy  -G(dn)  to KKR-matrix
C
         DO J = 1,NKKR
            CALL ZCOPY(NKKR,MAUX(1,J),1,MAUX2(NKKR+1,NKKR+J),1)
         END DO
C
C
         DO IQ = 1,NQ
            N = NLMQ(IQ)
C
            I1DN = 1
            I1UP = N + 1
C
            II1DN = IND0Q(IQ) + 1
            II1UP = NKKR + II1DN
C
            JJ0DN = IND0Q(IQ)
            JJ0UP = NKKR + JJ0DN
C
            DO J = 1,N
               JDN = J
               JUP = N + J
               JJDN = JJ0DN + J
               JJUP = JJ0UP + J
               CALL ZAXPY(N,C1,MSSQ(I1DN,JDN,IQ),1,MAUX2(II1DN,JJDN),1)
               CALL ZAXPY(N,C1,MSSQ(I1UP,JUP,IQ),1,MAUX2(II1UP,JJUP),1)
               CALL ZCOPY(N,MSSQ(I1DN,JUP,IQ),1,MAUX2(II1DN,JJUP),1)
               CALL ZCOPY(N,MSSQ(I1UP,JDN,IQ),1,MAUX2(II1UP,JJDN),1)
            END DO
         END DO
C
         CALL CINVLU(MAUX2,TAUK2,NKKR*2,NKKR*2)
C
C------------------------------------------------------------ store SUMQ
         WK = WKTAB(IK)
         WKSUM = WKSUM + WK
         CWK = DCMPLX(WK,0D0)
         DO IQ = 1,NQ
            N = NLMQ(IQ)
C
            I1DN = 1
            I1UP = N + 1
C
            II1DN = IND0Q(IQ) + 1
            II1UP = NKKR + II1DN
C
            JJ0DN = IND0Q(IQ)
            JJ0UP = NKKR + JJ0DN
C
            DO J = 1,N
               JDN = J
               JUP = N + J
               JJDN = JJ0DN + J
               JJUP = JJ0UP + J
               CALL ZAXPY(N,CWK,TAUK2(II1DN,JJDN),1,SUMQ(I1DN,JDN,IQ),1)
               CALL ZAXPY(N,CWK,TAUK2(II1UP,JJUP),1,SUMQ(I1UP,JUP,IQ),1)
               CALL ZAXPY(N,CWK,TAUK2(II1DN,JJUP),1,SUMQ(I1DN,JUP,IQ),1)
               CALL ZAXPY(N,CWK,TAUK2(II1UP,JJDN),1,SUMQ(I1UP,JDN,IQ),1)
            END DO
C
C            CALL  CMATSTR('  MSSQ   run A   ',SUMQ(1,1,IQ),NKMMAX,
C     &                  NKMMAX,0,0,0,1.d-6,6)
C
         END DO
C-----------------------------------------------------------------------
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C ----------------------------- perform the sum over symmetry operations
C
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C-----------------------------------------------------------------------
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            DO IQ = 1,NQ
C
               N = NKMQ(IQ)
               IQP = IQORGQP(ISYM,IQ)
C-----------------------------------------------------------------------
               CALL ZGEMM('N','N',N,N,N,C1,DROT(1,1,ISYM),NKMMAX,
     &                    SUMQ(1,1,IQP),NKMMAX,C0,W1,NKMMAX)
C
               CALL ZGEMM('N','C',N,N,N,C1,W1,NKMMAX,DROT(1,1,ISYM),
     &                    NKMMAX,C1,TAUQ(1,1,IQ),NKMMAX)
C.......................................................................
            END DO
         END IF
      END DO
C
      CSCL = 1D0/DBLE(NSYMACCEPTED*WKSUM)
C
      DO IQ = 1,NQ
         DO J = 1,NKMQ(IQ)
            CALL ZSCAL(NKMQ(IQ),CSCL,TAUQ(1,J,IQ),1)
         END DO
      END DO
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
C
C
         IF ( IPRINT.EQ.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
         IF ( CPAERR.LE.CPATOL ) THEN
            ICPACONV = 1
            IF ( IPRINT.GT.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
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
      DEALLOCATE (MAUX,TAUK,NLMQ,SUMQ,W1)
C
      IF ( IPRINT.GE.2 .AND. NCPA.GT.0 )
     &      WRITE (12,'(''E '',2F10.5,'' CPA '',I5,3E12.5)') ERYD,
     &     ICPAFLAG,CPAERR,CPACORR,CPACHNG
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
