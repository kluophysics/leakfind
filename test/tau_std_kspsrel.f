C*==tau_std_kspsrel.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_STD_KSPSREL(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                           ICPACONV,PHASK,IE,TSST,MSST,TSSQ,MSSQ,
     &                           TAUQ)
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
C   *    this is the special version for the                           *
C   *    spin-polarized scalar - relativistic mode   IREL = 2          *
C   *                                                                  *
C   *    the multiple scattering problem is dealt with using           *
C   *    the (l,ml)-representation for REAL spherical harmonics !!     *
C   *    for this reason the array dimensions differ from the          *
C   *    standard settings                                             *
C   *                                                                  *
C   *    - LLOYD = .TRUE.: the k-dependent terms in the Lloyd formula  *
C   *                      are calculated                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:WKTAB,NGFEP,GFEP,NKTAB,KTAB
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP,NSYMMAX
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,DROTQ
      USE MOD_ANGMOM,ONLY:NKMMAX,NLMQ,IND0Q,NLM,NKKR
      USE MOD_CALCMODE,ONLY:ITEST,LLOYD,KMROT
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      IMPLICIT NONE
C*--TAU_STD_KSPSREL43
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_STD_KSPSREL')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IE,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(IE),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 CPACORR,CPAERR,CPAERRL,GX,GY,GZ,KN2,WK,WKSUM
      COMPLEX*16 CWK,DETK,DETKFE,DETKGT,DROT2(:,:,:),EDU,MAUX(:,:),
     &           MQS(:,:,:,:),MQSLU(:,:,:),MTS(:,:,:,:),PHASKS(2),
     &           SUMQ(:,:,:),TAUK(:,:),TAUQS(:,:,:,:),TQS(:,:,:,:),
     &           TTS(:,:,:,:),W1(:,:),WLM(:,:)
      INTEGER I,I1,IA_ERR,ICALL,IK,IQ,IS,ISYM,IT,J,J1,N
      SAVE DROT2
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE TAUK,MAUX,SUMQ,WLM,W1,MQSLU,DROT2
      ALLOCATABLE TQS,TTS,MQS,MTS,TAUQS
C
      ICALL = ICALL + 1
      IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KMROT <> 0')
C
      ALLOCATE (TAUK(NKKR,NKKR),W1(NLM,NLM),WLM(NLM,NLM))
      ALLOCATE (MAUX(NKKR,NKKR),SUMQ(NLM,NLM,NQ))
      ALLOCATE (TTS(NLM,NLM,NT,2),TQS(NLM,NLM,NQ,2))
      ALLOCATE (MTS(NLM,NLM,NT,2),MQS(NLM,NLM,NQ,2))
      ALLOCATE (TAUQS(NLM,NLM,NQ,2),MQSLU(NLM,NLM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MTS')
C
C --------- store rotation matrices DROT with appropriate dimensions NLM
C
      IF ( ICALL.EQ.1 ) THEN
         ALLOCATE (DROT2(NLM,NLM,NSYMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROT2')
C
         DO ISYM = 1,NSYM
            DO J = 1,NLM
               CALL ZCOPY(NLM,DROT(1,J,ISYM),1,DROT2(1,J,ISYM),1)
            END DO
         END DO
      END IF
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      EDU = ERYD/(2*PI/ALAT)**2
C
      DO IQ = 1,NQ
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,MSSQ(1,J,IQ),1,MQS(1,J,IQ,1),1)
            CALL ZCOPY(NLM,MSSQ(NLM+1,NLM+J,IQ),1,MQS(1,J,IQ,2),1)
         END DO
C
      END DO
C
      IF ( NCPA.NE.0 ) THEN
         DO IT = 1,NT
            DO J = 1,NLM
               CALL ZCOPY(NLM,MSST(1,J,IT),1,MTS(1,J,IT,1),1)
               CALL ZCOPY(NLM,MSST(NLM+1,NLM+J,IT),1,MTS(1,J,IT,2),1)
            END DO
         END DO
C
      END IF
C
      PHASK(IE) = C0
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO IS = 1,2
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         ICPACONV = 0
         CPAERRL = 1.0D+6
         ITCPA = 0
 50      CONTINUE
         ITCPA = ITCPA + 1
C
         CALL CINIT(NLM*NLM*NQ,SUMQ)
C
         IF ( LLOYD ) THEN
C
            CALL CINIT(NLM*NLM*NQ,MQSLU)
C
            DO IQ = 1,NQ
               N = NLMQ(IQ)
               DO J = 1,N
                  CALL ZCOPY(N,MQS(1,J,IQ,IS),1,MQSLU(1,J,IQ),1)
               END DO
               CALL CINVLU(MQSLU(1,1,IQ),WLM,N,NLM)
            END DO
C
         END IF
C
         DETKGT = C0
         DETKFE = C0
         WKSUM = 0D0
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
         DO IK = 1,NKTAB
C
            CALL STRSET(IK,KTAB(1,IK),TAUK,MAUX,P)
C
            DO IQ = 1,NQ
               I1 = IND0Q(IQ) + 1
               N = NLMQ(IQ)
               DO J = 1,N
                  J1 = IND0Q(IQ) + J
                  CALL ZAXPY(N,C1,MQS(1,J,IQ,IS),1,MAUX(I1,J1),1)
               END DO
            END DO
C
            CALL CINVLU(MAUX,TAUK,NKKR,NKKR)
C
C------------------------------------------------------------ store TAUQ
            WK = WKTAB(IK)
            WKSUM = WKSUM + WK
            CWK = DCMPLX(WK,0D0)
            DO IQ = 1,NQ
               I1 = IND0Q(IQ) + 1
               N = NLMQ(IQ)
               DO J = 1,N
                  J1 = IND0Q(IQ) + J
                  CALL ZAXPY(N,CWK,TAUK(I1,J1),1,SUMQ(1,J,IQ),1)
               END DO
            END DO
C
C--------------------------------------------------------- LLOYD FORMULA
            IF ( LLOYD ) THEN
               DETK = C0
               DO IQ = 1,NQ
                  I1 = IND0Q(IQ)
                  N = NLMQ(IQ)
                  DO J = 1,N
                     I = I1 + J
                     DETK = DETK + LOG(MAUX(I,I)/MQSLU(J,J,IQ))
                  END DO
               END DO
               DETKGT = DETKGT - WK*DETK
C
               DETK = C0
               DO I = 1,NGFEP
                  GX = GFEP(1,I)
                  GY = GFEP(2,I)
                  GZ = GFEP(3,I)
                  KN2 = (KTAB(1,IK)+GX)**2 + (KTAB(2,IK)+GY)
     &                  **2 + (KTAB(3,IK)+GZ)**2
                  DETK = DETK + LOG(KN2-EDU)
               END DO
C
               DETKFE = DETKFE - WK*DETK
            END IF
C-----------------------------------------------------------------------
         END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         DETKFE = DETKFE/DBLE(WKSUM)
         DETKGT = DETKGT/DBLE(WKSUM)
C
         PHASKS(IS) = DETKGT + DETKFE
C
         CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,WKSUM,SUMQ,
     &                  TAUQS(1,1,1,IS),W1,NQ,NLMQ,DROT2,IQORGQP,
     &                  SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,NQMAX,
     &                  NLM)
C
         IF ( NCPA.GT.0 ) THEN
C
            IF ( ICPAALG.EQ.1 ) THEN
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NLMQ,
     &                       NOQ,ITOQ,CONC,TTS(1,1,1,IS),MTS(1,1,1,IS),
     &                       TQS(1,1,1,IS),MQS(1,1,1,IS),TAUQS(1,1,1,IS)
     &                       ,NTMAX,NQMAX,NLM)
C
            ELSE
               CALL CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,
     &                        NLMQ,NOQ,ITOQ,CONC,TTS(1,1,1,IS),
     &                        MTS(1,1,1,IS),TQS(1,1,1,IS),MQS(1,1,1,IS),
     &                        TAUQS(1,1,1,IS),KMROT,DROTQ,NTMAX,NQMAX,
     &                        NLM)
            END IF
C
            CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,SUMQ,MQS(1,1,1,IS),
     &                     W1,NQ,NLMQ,DROT2,IQORGQP,SYMUNITARY,
     &                     SYMACCEPTED,NSYM,NSYMACCEPTED,NQMAX,NLM)
C
            IF ( IPRINT.EQ.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
            IF ( CPAERR.LE.CPATOL ) THEN
               ICPACONV = 1
               IF ( IPRINT.GT.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &              CPACHNG
            ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
               WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 1
            ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
               WRITE (6,99003) ITCPA
               WRITE (6,99004) CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 2
            ELSE
               CPAERRL = CPAERR
               GOTO 50
            END IF
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      PHASK(IE) = PHASKS(1) + PHASKS(2)
C
C ------------------------------------------ now copy to standard arrays
C
      DO IQ = 1,NQ
C
         IF ( NCPA.NE.0 ) THEN
            CALL CINIT(NKMMAX*NKMMAX,MSSQ(1,1,IQ))
            DO J = 1,NLM
               CALL ZCOPY(NLM,MQS(1,J,IQ,1),1,MSSQ(1,J,IQ),1)
               CALL ZCOPY(NLM,MQS(1,J,IQ,2),1,MSSQ(NLM+1,NLM+J,IQ),1)
            END DO
         END IF
C
         CALL CINIT(NKMMAX*NKMMAX,TAUQ(1,1,IQ))
         DO J = 1,NLM
            CALL ZCOPY(NLM,TAUQS(1,J,IQ,1),1,TAUQ(1,J,IQ),1)
            CALL ZCOPY(NLM,TAUQS(1,J,IQ,2),1,TAUQ(NLM+1,NLM+J,IQ),1)
         END DO
C
      END DO
C
      IF ( ITEST.EQ.5 ) THEN
         WRITE (77,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(PHASK(IE))
         WRITE (78,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKFE)
         WRITE (79,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKGT)
      END IF
C
      DEALLOCATE (TAUK,MAUX,SUMQ,TAUQS,WLM,W1,MQS,MTS,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      IF ( IPRINT.GE.2 .AND. NCPA.GT.0 )
     &      WRITE (12,'(''E '',2F10.5,'' CPA '',I5,3E12.5)') ERYD,
     &     ICPAFLAG,CPAERR,CPACORR,CPACHNG
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
