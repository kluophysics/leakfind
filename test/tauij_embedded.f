C*==tauij_embedded.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_EMBEDDED(IPRINT,ERYD,ICPAFLAG,CPACHNG,ITCPA,
     &                          ICPACONV,IECURR,TSST,MSST,TSSQ,MSSQ,
     &                          TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *    calculate the TAU-matrix for EMBEDDED-CLUSTER case            *
C   *                                                                  *
C   *    CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED:                             *
C   *                                                                  *
C   *    FALSE:  only diagonal elements are calculated                 *
C   *    TRUE:   all elements are calculated                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_CALCMODE,ONLY:IREL,KMROT,RELAX_CLU
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,NCPA
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,NQCLU,IQBOT_CLU,IQ_QCLU,
     &    NQHOST,IQBOT_HOST,IQTOP_HOST,QBAS_QCLU,QBAS0_QCLU
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP,NSYMMAX,SYMMETRIZE_TAU
      USE MOD_TYPES,ONLY:CONC,NTMAX,NT
      USE MOD_ANGMOM,ONLY:NKMMAX,IXM0_QCLU,NKKR_CLU,NSPIN,NXM,NXM_QCLU,
     &    NKM,NLM,NLMMAX,NKMQ,NLMQ,WKM1
      USE MOD_FILES,ONLY:IFILTAUIJ
      USE MOD_TAUIJ,ONLY:NTAUIJ,TAUIJ,CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED,
     &    TAUQQ_CLU
      IMPLICIT NONE
C*--TAUIJ_EMBEDDED31
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_EMBEDDED')
      REAL*8 WGTK
      PARAMETER (WGTK=999999D0)
      LOGICAL ADD
      PARAMETER (ADD=.FALSE.)
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD
      INTEGER ICPACONV,ICPAFLAG,IECURR,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 CPACORR,CPAERR,CPAERRL,SVEC(3)
      COMPLEX*16 DSYMX(:,:,:),ERYD_IN,MAUX(:,:),MSSQX(:,:,:),
     &           MSSTX(:,:,:),TAUINV(:,:),TAUQX(:,:,:),TSSQX(:,:,:),
     &           TSSTX(:,:,:),UDIA(:),UMAT(:,:,:),WLM(:,:),WXA(:,:),
     &           WXB(:,:)
      INTEGER I,I0,IA_ERR,INFO,IOFF,IPIV(:),IQ,IQBOT_HOST_IN,IQCLU,
     &        IQCLUORGQCLUP(:,:),IQCLUP,IQEXT,IQEXT1,IQEXT2,IQHOST,IQP,
     &        IQTOP_HOST_IN,IREC,ISPIN,ISPIN_IN,ISYM,J,JJ,NI,NKM_IN,
     &        NTAUIJ_IN,NXMMAX,NXM_IN
      LOGICAL INITIALIZE,RESTRICT_TO_DIAGONAL_TAU
      SAVE DSYMX
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE WXA,WXB,WLM
      ALLOCATABLE TSSQX,TAUQX,DSYMX,TSSTX,MSSTX,MSSQX
      ALLOCATABLE TAUINV,UDIA,UMAT,IPIV,MAUX
      ALLOCATABLE IQCLUORGQCLUP
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (IQCLUORGQCLUP(NSYMMAX,NQCLU))
      IQCLUORGQCLUP(:,:) = 0
      DO IQCLUP = 1,NQCLU
         IQP = NQHOST + IQCLUP
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) IQCLUORGQCLUP(ISYM,IQCLUP)
     &           = IQORGQP(ISYM,IQP) - NQHOST
         END DO
      END DO
C
      ALLOCATE (IPIV(NKKR_CLU),MAUX(NKKR_CLU,NKKR_CLU))
      ALLOCATE (TAUINV(NKKR_CLU,NKKR_CLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAU')
C
C ----------------------------------------------------------------------
      IF ( INITIALIZE ) THEN
         IF ( IREL.EQ.2 ) THEN
            ALLOCATE (DSYMX(NXM,NXM,NSYMMAX))
            DSYMX(1:NLM,1:NLM,1:NSYMMAX) = DROT(1:NLM,1:NLM,1:NSYMMAX)
         ELSE
            ALLOCATE (DSYMX(1,1,1))
         END IF
         INITIALIZE = .FALSE.
      END IF
C ----------------------------------------------------------------------
C
      ALLOCATE (TAUQX(NXM,NXM,NQ),UDIA(NKKR_CLU))
      IF ( IREL.EQ.2 ) THEN
         ALLOCATE (TSSTX(NLM,NLM,NT))
         ALLOCATE (MSSTX(NLM,NLM,NT),MSSQX(NLM,NLM,NQ))
      END IF
C
      ALLOCATE (TSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUQX')
      IF ( .NOT.ALLOCATED(TAUIJ) ) THEN
         ALLOCATE (TAUIJ(NXM,NXM,NSPIN,NTAUIJ),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUIJ')
      END IF
C
      NXMMAX = NKMMAX
C
      IF ( CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED .AND. 
     &     .NOT.ALLOCATED(TAUQQ_CLU) ) THEN
         ALLOCATE (TAUQQ_CLU(NXMMAX,NXMMAX,NQCLU,NQCLU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &        CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUQQ_CLU')
      END IF
C
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF ( RELAX_CLU ) THEN
         ALLOCATE (WLM(NLM,NLM),WXA(NXMMAX,NXMMAX),WXB(NXMMAX,NXMMAX))
         ALLOCATE (UMAT(NXMMAX,NXMMAX,NQCLU),STAT=IA_ERR)
      END IF
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
      IQEXT1 = NQHOST + 1
      IQEXT2 = NQHOST + NQCLU
C
      IF ( IREL.EQ.2 ) THEN
         TAUQ(1:NLM,NLM+1:NLM+NLM,IQEXT1:IQEXT2) = C0
         TAUQ(NLM+1:NLM+NLM,1:NLM,IQEXT1:IQEXT2) = C0
      END IF
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO ISPIN = 1,NSPIN
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( .NOT.(CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED) ) THEN
            RESTRICT_TO_DIAGONAL_TAU = .TRUE.
         ELSE IF ( NCPA.GT.0 ) THEN
            RESTRICT_TO_DIAGONAL_TAU = .TRUE.
         ELSE
            RESTRICT_TO_DIAGONAL_TAU = .FALSE.
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C-----------------------------------------------------------------------
C                      read  inverse TAU-matrix
C-----------------------------------------------------------------------
C
         IREC = 1 + NSPIN*(IECURR-1) + ISPIN
C
         READ (IFILTAUIJ,REC=IREC) ERYD_IN,IQBOT_HOST_IN,IQTOP_HOST_IN,
     &                             NTAUIJ_IN,NKM_IN,NXM_IN,ISPIN_IN,
     &                             (TAUQ(1:NKM,1:NKM,IQ),
     &                             MSSQ(1:NKM,1:NKM,IQ),IQ=IQBOT_HOST,
     &                             IQTOP_HOST),
     &                             ((TAUINV(I,J),I=1,NKKR_CLU),J=1,
     &                             NKKR_CLU)
C
         IF ( ISPIN_IN.NE.ISPIN .OR. NKM_IN.NE.NKM ) THEN
C
            WRITE (6,99001) ERYD_IN,IQBOT_HOST_IN,IQTOP_HOST_IN,
     &                      NTAUIJ_IN,NKM_IN,NXM_IN,ISPIN_IN
            CALL STOP_MESSAGE(ROUTINE,
     &                        'ISPIN_IN<>ISPIN .OR. NKM_IN<>NKM')
         END IF
C
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         IF ( RELAX_CLU ) THEN
C-----------------------------------------------------------------------
C                         set up  U-matrix
C-----------------------------------------------------------------------
C
            DO IQCLU = 1,NQCLU
C
               SVEC(1:3) = (QBAS_QCLU(1:3,IQCLU)-QBAS0_QCLU(1:3,IQCLU))
     &                     *ALAT
C
               IF ( IREL.LE.2 ) THEN
C
                  CALL INIT_UMAT(SVEC,WLM,ERYD)
C
                  UMAT(1:NLM,1:NLM,IQCLU) = WLM(1:NLM,1:NLM)
C
               ELSE
C
                  CALL INIT_UMAT(SVEC,WLM,ERYD)
C
                  WKM1(1:NKMMAX,1:NKMMAX) = C0
C
                  WKM1(1:NLM,1:NLM) = WLM(1:NLM,1:NLM)
                  WKM1(NLM+1:NKMMAX,NLM+1:NKMMAX) = WLM(1:NLM,1:NLM)
C
                  CALL CHANGEREP(NKM,NKMMAX,WKM1(1,1),'RLM>REL',
     &                           UMAT(1,1,IQCLU))
C
               END IF
C
            END DO
C
         END IF
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         ICPACONV = 0
         CPAERRL = 1.0D+6
         ITCPA = 0
 50      CONTINUE
         ITCPA = ITCPA + 1
C
C ----------------------------------------------------------------------
         IF ( IREL.EQ.2 ) THEN
            IOFF = NLM*(ISPIN-1)
            DO IQ = 1,NQ
               MSSQX(1:NLM,1:NLM,IQ) = MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+
     &                                 NLM,IQ)
            END DO
         END IF
C
C=======================================================================
C                           set up  KKR-matrix
C=======================================================================
C
         DO IQCLU = 1,NQCLU
            I0 = IXM0_QCLU(IQCLU)
            NI = NXM_QCLU(IQCLU)
            IQHOST = IQ_QCLU(IQCLU)
            IQEXT = NQHOST + IQCLU
C
C unrelaxed uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
            IF ( .NOT.RELAX_CLU ) THEN
C
               IF ( IREL.EQ.2 ) THEN
                  DO J = 1,NI
                     JJ = I0 + J
                     DO I = 1,NI
                        TAUINV(I0+I,JJ) = TAUINV(I0+I,JJ)
     &                     + MSSQX(I,J,IQEXT) - MSSQX(I,J,IQHOST)
                     END DO
                  END DO
               ELSE
                  DO J = 1,NI
                     JJ = I0 + J
                     DO I = 1,NI
                        TAUINV(I0+I,JJ) = TAUINV(I0+I,JJ)
     &                     + MSSQ(I,J,IQEXT) - MSSQ(I,J,IQHOST)
                     END DO
                  END DO
               END IF
C
C unrelaxed uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C-----------------------------------------------------------------------
C     transfer t-matrix of cluster atoms to global/unshifted frame
C-----------------------------------------------------------------------
C
            ELSE IF ( IREL.EQ.2 ) THEN
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,UMAT(1,1,IQCLU),NXMMAX,
     &                    MSSQX(1,1,IQEXT),NXMMAX,C0,WXA,NXMMAX)
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WXA,NXMMAX,
     &                    UMAT(1,1,IQCLU),NXMMAX,C0,WXB,NXMMAX)
C
               DO J = 1,NI
                  JJ = I0 + J
                  DO I = 1,NI
                     TAUINV(I0+I,JJ) = TAUINV(I0+I,JJ) + WXB(I,J)
     &                                 - MSSQX(I,J,IQHOST)
                  END DO
               END DO
C
            ELSE
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,UMAT(1,1,IQCLU),NXMMAX,
     &                    MSSQ(1,1,IQEXT),NXMMAX,C0,WXA,NXMMAX)
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WXA,NXMMAX,
     &                    UMAT(1,1,IQCLU),NXMMAX,C0,WXB,NXMMAX)
C
               DO J = 1,NI
                  JJ = I0 + J
                  DO I = 1,NI
                     TAUINV(I0+I,JJ) = TAUINV(I0+I,JJ) + WXB(I,J)
     &                                 - MSSQ(I,J,IQHOST)
                  END DO
               END DO
C
            END IF
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
         END DO
C
C=======================================================================
C                           invert KKR-matrix
C=======================================================================
C
C------------------------------- calculate only diagonal elements of TAU
C
         IF ( RESTRICT_TO_DIAGONAL_TAU ) THEN
C
            IF ( IREL.EQ.2 ) THEN
C
               CALL CINVDIABLK(NKKR_CLU,TAUINV,UDIA,TAUQX(1,1,IQBOT_CLU)
     &                         ,NXM,NQCLU,IXM0_QCLU,NXM_QCLU,ADD,WGTK)
C
            ELSE
C
               CALL CINVDIABLK(NKKR_CLU,TAUINV,UDIA,TAUQ(1,1,IQBOT_CLU),
     &                         NKMMAX,NQCLU,IXM0_QCLU,NXM_QCLU,ADD,WGTK)
C
            END IF
C
C----------------------------------------- calculate ALL elements of TAU
C                             NOTE: after inversion TAUINV contains TAU
         ELSE IF ( IREL.EQ.2 ) THEN
C
            CALL CINVDIABLK(NKKR_CLU,TAUINV,UDIA,TAUQX(1,1,IQBOT_CLU),
     &                      NXM,NQCLU,IXM0_QCLU,NXM_QCLU,ADD,WGTK)
C
            CALL STOP_MESSAGE(ROUTINE,'IREL.EQ.2 not allowed')
C
         ELSE
C
            CALL ZGETRF(NKKR_CLU,NKKR_CLU,TAUINV,NKKR_CLU,IPIV,INFO)
C
            CALL ZGETRI(NKKR_CLU,TAUINV,NKKR_CLU,IPIV,MAUX,
     &                  NKKR_CLU*NKKR_CLU,INFO)
C
            CALL TAU_SPLIT_QQ_CLU(NKMMAX,NQMAX,NQHOST,NKMQ,NQCLU,
     &                            NKKR_CLU,TAUINV,TAUQ,TAUQQ_CLU)
C
         END IF
C
C=======================================================================
C
         IF ( IPRINT.GE.1 ) WRITE (6,99006) IECURR,ERYD
C
C=======================================================================
C                       symmetrize TAU matrix
C=======================================================================
C
         SYMMETRIZE_TAU = .TRUE.
         WRITE (*,*) '***********************************'
         WRITE (*,*) ROUTINE
         WRITE (*,*) 'SYMMETRIZE_TAU ',SYMMETRIZE_TAU
         WRITE (*,*) '***********************************'
         WRITE (*,*) '***********************************'
         IF ( SYMMETRIZE_TAU ) THEN
            IF ( IREL.EQ.2 ) THEN
C
               CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,TAUINV,
     &                        TAUQX(1,1,IQBOT_CLU),WKM1,NQCLU,NXM_QCLU,
     &                        DSYMX,IQCLUORGQCLUP,SYMUNITARY,
     &                        SYMACCEPTED,NSYM,NSYMACCEPTED,NQCLU,NXM)
C
            ELSE
C
               CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,TAUINV,
     &                        TAUQ(1,1,IQBOT_CLU),WKM1,NQCLU,
     &                        NKMQ(IQBOT_CLU),DROT,IQCLUORGQCLUP,
     &                        SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                        NQCLU,NKMMAX)
C
            END IF
         END IF
C=======================================================================
C
         IF ( NCPA.GT.0 ) THEN
C
            IF ( IREL.NE.2 ) THEN
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,
     &                       ICPA(IQBOT_CLU),NQCLU,NKMQ(IQBOT_CLU),
     &                       NOQ(IQBOT_CLU),ITOQ(1,IQBOT_CLU),CONC,TSST,
     &                       MSST,TSSQ(1,1,IQBOT_CLU),
     &                       MSSQ(1,1,IQBOT_CLU),TAUQ(1,1,IQBOT_CLU),
     &                       NTMAX,NQCLU,NKMMAX)
C
Cc             CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,
Cc     &    1D0,TAUQX(1,1,IQBOT_CLU),MSSQ(1,1,
Cc   &              IQBOT_CLU),WA,NQCLU,NKMQ(IQBOT_CLU),DROT,
Cc     &              IQTBORGQTBP,
Cc   &              SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,NQCLU,
Cc   &              NKMMAX)
            ELSE
C
               IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KMROT <> 0')
C
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,
     &                       ICPA(IQBOT_CLU),NQCLU,NLMQ(IQBOT_CLU),
     &                       NOQ(IQBOT_CLU),ITOQ(1,IQBOT_CLU),CONC,
     &                       TSSTX,MSSTX,TSSQX(1,1,IQBOT_CLU),
     &                       MSSQX(1,1,IQBOT_CLU),TAUQX(1,1,IQBOT_CLU),
     &                       NTMAX,NQCLU,NLMMAX)
C
Cc             CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,
Cc     &    1D0,DSSQX(1,1,IQBOT_CLU),
Cc   &              MSSQX(1,1,IQBOT_CLU),WX,NQCLU,NLMQ(IQBOT_CLU),DSYMX,
Cc   &              IQTBORGQTBP,SYMUNITARY,SYMACCEPTED,NSYM,
Cc   &              NSYMACCEPTED,NQCLU,NLMMAX)
C
               MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,1:NQ)
     &            = MSSQX(1:NLM,1:NLM,1:NQ)
C
            END IF
C
            IF ( IPRINT.GE.1 ) WRITE (6,99005) CPAERR,CPACORR,CPACHNG
C
            IF ( CPAERR.LE.CPATOL ) THEN
               ICPACONV = 1
               IF ( IPRINT.GT.0 ) WRITE (6,99002) ITCPA,CPAERR,CPACORR,
     &              CPACHNG
            ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
               WRITE (6,99003) ITCPA,CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 1
            ELSE IF ( CPAERR.GT.20000*CPAERRL ) THEN
               WRITE (6,99004) ITCPA
               WRITE (6,99005) CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 2
            ELSE
               CPAERRL = CPAERR
               GOTO 50
            END IF
C
C--- CPA converged - do an extra run to get FULL TAU matrix if requested
            IF ( CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED .AND. 
     &           .NOT.RESTRICT_TO_DIAGONAL_TAU ) THEN
               RESTRICT_TO_DIAGONAL_TAU = .FALSE.
               GOTO 50
            END IF
C
         END IF
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C-----------------------------------------------------------------------
C   transfer TAU-matrix of cluster atoms back to local/shifted frame
C-----------------------------------------------------------------------
         IF ( RELAX_CLU ) THEN
C
            DO IQCLU = 1,NQCLU
               IQEXT = NQHOST + IQCLU
C
               IF ( IREL.EQ.2 ) THEN
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,UMAT(1,1,IQCLU),
     &                       NXMMAX,TAUQX(1,1,IQEXT),NXMMAX,C0,WXA,
     &                       NXMMAX)
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WXA,NXMMAX,
     &                       UMAT(1,1,IQCLU),NXMMAX,C0,TAUQX(1,1,IQEXT),
     &                       NXMMAX)
               ELSE
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,UMAT(1,1,IQCLU),
     &                       NXMMAX,TAUQ(1,1,IQEXT),NXMMAX,C0,WKM1,
     &                       NXMMAX)
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM1,NXMMAX,
     &                       UMAT(1,1,IQCLU),NXMMAX,C0,TAUQ(1,1,IQEXT),
     &                       NXMMAX)
C
               END IF
C
            END DO
C
         END IF
C relax xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
         IF ( IREL.EQ.2 ) THEN
C
            IQEXT1 = NQHOST + 1
            IQEXT2 = NQHOST + NQCLU
C
            TAUQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,IQEXT1:IQEXT2)
     &         = TAUQX(1:NLM,1:NLM,IQEXT1:IQEXT2)
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DEALLOCATE (TAUIJ)
C
99001 FORMAT (//,10X,'TROUBLE reading TAUIJ - file',/,10X,
     &        'read from file:',/,10X,'ERYD_IN        ',2E12.5,/,10X,
     &        'IQBOT_HOST_IN  ',I5,/,10X,'IQTOP_HOST_IN  ',I5,/,10X,
     &        'NTAUIJ_IN      ',I5,/,10X,'NKM_IN         ',I5,/,10X,
     &        'NXM_IN         ',I5,/,10X,'ISPIN_IN       ',I5,/)
99002 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99003 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99004 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99005 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
99006 FORMAT (/,10X,'<TAUIJ_EMBEDDED>:     IE=',I3,2X,'ERYD=',2F12.5)
      END
C*==tau_split_qq_clu.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_SPLIT_QQ_CLU(NXMMAX,NQMAX,NQHOST,NXMQ,NQCLU,NQXM,
     &                            TAU,TAUQ,TAUQQ_CLU)
C   ********************************************************************
C   *                                                                  *
C   *  map the TAU matrix of an embedded cluster calculation to        *
C   *  TAUQQ_CLU indexed individually by the cluster sites IQCLU       *
C   *                                                                  *
C   *  TAUQ contains the standard site diagonal TAU matrix             *
C   *       indexed by the common IQ index   IQ = IQHOST + IQCLU       *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--TAU_SPLIT_QQ_CLU533
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_SPLIT_QQ_CLU')
C
C Dummy arguments
C
      INTEGER NQCLU,NQHOST,NQMAX,NQXM,NXMMAX
      INTEGER NXMQ(NQMAX)
      COMPLEX*16 TAU(NQXM,NQXM),TAUQ(NXMMAX,NXMMAX,NQMAX),
     &           TAUQQ_CLU(NXMMAX,NXMMAX,NQCLU,NQCLU)
C
C Local variables
C
      INTEGER IKMQ1,IKMQ2,IQ,IQCLU,JKMQ1,JKMQ2,JQ,JQCLU,NXMI,NXMJ
C
C*** End of declarations rewritten by SPAG
C
      JKMQ1 = 1
      DO JQCLU = 1,NQCLU
         JQ = NQHOST + JQCLU
         NXMJ = NXMQ(JQ)
         JKMQ2 = JKMQ1 - 1 + NXMJ
C
         IKMQ1 = 1
         DO IQCLU = 1,NQCLU
            IQ = NQHOST + IQCLU
            NXMI = NXMQ(IQ)
            IKMQ2 = IKMQ1 - 1 + NXMI
C
            TAUQQ_CLU(1:NXMI,1:NXMJ,IQCLU,JQCLU)
     &         = TAU(IKMQ1:IKMQ2,JKMQ1:JKMQ2)
C
            IF ( IQCLU.EQ.JQCLU ) TAUQ(1:NXMI,1:NXMJ,IQ)
     &           = TAU(IKMQ1:IKMQ2,JKMQ1:JKMQ2)
C
            IKMQ1 = IKMQ1 + NXMI
         END DO
C
         JKMQ1 = JKMQ1 + NXMJ
      END DO
C
      IF ( IKMQ2.NE.NQXM ) CALL STOP_MESSAGE(ROUTINE,'IKMQ2 != NQXM')
      IF ( JKMQ2.NE.NQXM ) CALL STOP_MESSAGE(ROUTINE,'JKMQ2 != NQXM')
C
      END
