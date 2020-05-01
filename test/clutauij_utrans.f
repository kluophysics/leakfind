C*==clutauij_utrans.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUTAUIJ_UTRANS(MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine to supply the host TAUIJ data                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SITES,ONLY:NQMAX,NQCLU,NQHOST,QBAS_QCLU,QBAS0_QCLU,
     &    IQ_QCLU,IQBOT_HOST,IQTOP_HOST
      USE MOD_ANGMOM,ONLY:NKMMAX,IXM0_QCLU,NKKR_CLU,NSPIN,NXM,NXM_QCLU,
     &    NKM,NLM,NLMMAX
      USE MOD_ENERGY,ONLY:IGRID,EMIN,EMAX,EFERMI,NETAB,NEPATH,ETAB
      USE MOD_FILES,ONLY:IFILTAUIJ,RECLNG_TAUIJ,IFILTAUIJ_TMP
      USE MOD_TAUIJ,ONLY:NTAUIJ,IQ_QTAB1_TAUIJ,JQ_QTAB2_TAUIJ,
     &    NQTAB2_TAUIJ,JQTAB_TAUIJMAX,IQ_TAUIJ,JQ_TAUIJ,N5VEC_TAUIJ,
     &    N123TAUIJMAX,ITAUIJ_LOOP,NLOOP_TAUIJ,NQTAB1_TAUIJ,
     &    NTAUIJ_QTAB2_TAUIJ,TAUIJ
      IMPLICIT NONE
C*--CLUTAUIJ_UTRANS22
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 ERYD,ERYD_IN,TAUINV(:,:),UMAT(:,:,:),WLM(:,:),WXA(:,:),
     &           WXB(:,:)
      INTEGER I,I1,I2,IA_ERR,IE,IECURR,IEPATH,IQ,IQBOT_HOST_IN,IQCLU,
     &        IQEXT1,IQEXT2,IQTOP_HOST_IN,IREC,ISPIN,ISPIN_IN,J,J1,J2,
     &        JQCLU,NKM_IN,NTAUIJ_IN,NXMMAX,NXM_IN
      REAL*8 SVEC(3)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WXA,WXB,WLM
      ALLOCATABLE TAUINV,UMAT
C
      ALLOCATE (WLM(NLM,NLM),WXA(NXM,NXM),WXB(NXM,NXM))
      ALLOCATE (TAUINV(NKKR_CLU,NKKR_CLU),UMAT(NXM,NXM,NQCLU),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUTAUIJ_UTRANS->TAU'
C
C ----------------------------------------------------------------------
C
      IF ( .NOT.ALLOCATED(TAUIJ) ) THEN
         ALLOCATE (TAUIJ(NXM,NXM,NSPIN,NTAUIJ),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc <CLUTAUIJ_UTRANS>:->TAUIJ'
      END IF
C
      IF ( IREL.LE.1 ) THEN
         NXMMAX = NLMMAX
      ELSE
         NXMMAX = NKMMAX
      END IF
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IEPATH = 1,NEPATH
         DO IE = 1,NETAB(IEPATH)
            IREC = 1 + NSPIN*(IE-1) + 1
C
            READ (IFILTAUIJ,REC=IREC) ETAB(IE,1)
         END DO
      END DO
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C ----------------------------------------------------------------------
C       open the (SCRATCH) TAU-data file for the shifted positions
C             and write first record with general information
C ----------------------------------------------------------------------
C
      OPEN (UNIT=IFILTAUIJ_TMP,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNG_TAUIJ,ERR=100)
C
      WRITE (IFILTAUIJ_TMP,REC=1) RECLNG_TAUIJ,IGRID(1),EMIN,EMAX,
     &                            EFERMI,NETAB(1),NTAUIJ,N123TAUIJMAX,
     &                            NLOOP_TAUIJ,NQTAB1_TAUIJ,
     &                            JQTAB_TAUIJMAX,
     &                            (IQ_TAUIJ(I),JQ_TAUIJ(I),
     &                            (N5VEC_TAUIJ(J,I),J=1,5),I=1,NTAUIJ),
     &                            (ITAUIJ_LOOP(I),I=1,NLOOP_TAUIJ),
     &                            (IQ_QTAB1_TAUIJ(I),NQTAB2_TAUIJ(I),
     &                            (JQ_QTAB2_TAUIJ(I,J),
     &                            NTAUIJ_QTAB2_TAUIJ(I,J),J=1,
     &                            JQTAB_TAUIJMAX),I=1,NQTAB1_TAUIJ)
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IEPATH = 1,NEPATH
C
         DO IE = 1,NETAB(IEPATH)
            IECURR = IE
C
            ERYD = ETAB(IE,1)
C
            IF ( IREL.EQ.2 ) THEN
C
               IQEXT1 = NQHOST + 1
               IQEXT2 = NQHOST + NQCLU
C
               TAUQ(1:NLM,NLM+1:NLM+NLM,IQEXT1:IQEXT2) = C0
               TAUQ(NLM+1:NLM+NLM,1:NLM,IQEXT1:IQEXT2) = C0
C
            END IF
C
C-----------------------------------------------------------------------
C                         set up  U-matrix
C-----------------------------------------------------------------------
C
            DO IQCLU = 1,NQCLU
C
               SVEC(1:3) = QBAS_QCLU(1:3,IQCLU) - QBAS0_QCLU(1:3,IQCLU)
C
               IF ( IREL.LE.2 ) THEN
C
                  CALL INIT_UMAT(SVEC,UMAT(1,1,IQCLU),ERYD)
C
               ELSE
C
                  CALL INIT_UMAT(SVEC,WLM,ERYD)
C
                  WXA(1:NKMMAX,1:NKMMAX) = C0
C
                  WXA(1:NLM,1:NLM) = WLM(1:NLM,1:NLM)
                  WXA(NLM+1:NKMMAX,NLM+1:NKMMAX) = WLM(1:NLM,1:NLM)
C
                  CALL CHANGEREP(NKM,NKMMAX,WXA(1,1),'RLM>REL',
     &                           UMAT(1,1,IQCLU))
C
               END IF
C
            END DO
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            DO ISPIN = 1,NSPIN
C
C-----------------------------------------------------------------------
C                         read  matrices
C-----------------------------------------------------------------------
C
               IREC = 1 + NSPIN*(IECURR-1) + ISPIN
C
               READ (IFILTAUIJ,REC=IREC) ERYD_IN,IQBOT_HOST_IN,
     &               IQTOP_HOST_IN,NTAUIJ_IN,NKM_IN,NXM_IN,ISPIN_IN,
     &               (TAUQ(1:NKM,1:NKM,IQ),MSSQ(1:NKM,1:NKM,IQ),
     &               IQ=IQBOT_HOST,IQTOP_HOST),
     &               TAUINV(1:NKKR_CLU,1:NKKR_CLU)
C
               IF ( (ABS(ERYD_IN-ERYD).GT.1D-6) .OR. 
     &              (IQBOT_HOST_IN.NE.IQBOT_HOST) .OR. 
     &              (IQTOP_HOST_IN.NE.IQTOP_HOST) .OR. 
     &              (NTAUIJ_IN.NE.NTAUIJ_IN) .OR. (NKM_IN.NE.NKM) .OR. 
     &              (NXM_IN.NE.NXM) .OR. (ISPIN_IN.NE.ISPIN) ) THEN
                  WRITE (6,*) 'TAUIJ-file inconsistent'
                  STOP 'in <CLUTAUIJ_UTRANS>:  TAUIJ-file inconsistent'
               END IF
C
C=======================================================================
C                       transform matrices
C=======================================================================
C
               DO IQCLU = 1,NQCLU
C
C-----------------------------------------------------------------------
C                transform site-diagonal  m- and TAU-matrices
C-----------------------------------------------------------------------
C
                  IQ = IQ_QCLU(IQCLU)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXMMAX,C1,UMAT(1,1,IQCLU),
     &                       NXMMAX,MSSQ(1,1,IQCLU),NXMMAX,C0,WXA,
     &                       NXMMAX)
                  CALL ZGEMM('N','T',NXMMAX,NXMMAX,NXMMAX,C1,WXA,NXMMAX,
     &                       UMAT(1,1,IQCLU),NXMMAX,C0,MSSQ(1,1,IQCLU),
     &                       NXMMAX)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXMMAX,C1,UMAT(1,1,IQCLU),
     &                       NXMMAX,TAUQ(1,1,IQCLU),NXMMAX,C0,WXA,
     &                       NXMMAX)
                  CALL ZGEMM('N','T',NXMMAX,NXMMAX,NXMMAX,C1,WXA,NXMMAX,
     &                       UMAT(1,1,IQCLU),NXMMAX,C0,TAUQ(1,1,IQCLU),
     &                       NXMMAX)
C
C-----------------------------------------------------------------------
C                transform FULL inverse host TAU-matrix
C-----------------------------------------------------------------------
C
                  I1 = IXM0_QCLU(IQCLU) + 1
                  I2 = I1 + NXM_QCLU(IQCLU)
C
                  DO JQCLU = 1,NQCLU
                     J1 = IXM0_QCLU(JQCLU) + 1
                     J2 = J1 + NXM_QCLU(JQCLU)
C
                     WXA(1:NXM,1:NXM) = TAUINV(I1:I2,J1:J2)
C
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,UMAT(1,1,IQCLU),
     &                          NXM,WXA,NXM,C0,WXB,NXM)
                     CALL ZGEMM('N','T',NXM,NXM,NXM,C1,WXB,NXM,
     &                          UMAT(1,1,JQCLU),NXM,C0,WXA,NXM)
C
                     TAUINV(I1:I2,J1:J2) = WXA(1:NXM,1:NXM)
C
                  END DO
               END DO
C
C-----------------------------------------------------------------------
C                          WRITE  matrices
C-----------------------------------------------------------------------
C
               WRITE (IFILTAUIJ_TMP,REC=IREC) ERYD_IN,IQBOT_HOST_IN,
     &                IQTOP_HOST_IN,NTAUIJ_IN,NKM_IN,NXM_IN,ISPIN_IN,
     &                (TAUQ(1:NKM,1:NKM,IQ),MSSQ(1:NKM,1:NKM,IQ),
     &                IQ=IQBOT_HOST,IQTOP_HOST),
     &                TAUINV(1:NKKR_CLU,1:NKKR_CLU)
C
            END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         END DO
      END DO
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C-----------------------------------------------------------------------
C             from now on: use the TAU data from the scratch file
C-----------------------------------------------------------------------
C
      CLOSE (IFILTAUIJ)
C
      IFILTAUIJ = IFILTAUIJ_TMP
C
      WRITE (6,99001)
C
      RETURN
 100  CONTINUE
      STOP 'TAUIJ_tmp file could not be opened for writing'
C-----------------------------------------------------------------------
99001 FORMAT (2(/,1X,79('*')),/,10X,'INFO from   <CLUTAUIJ_UTRANS>',/,
     &        10X,'the TAUIJ - data for the host have been created',
     &        ' for the shifted positions',2(/,1X,79('*')),/)
      END
