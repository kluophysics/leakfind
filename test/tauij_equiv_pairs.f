C*==tauij_equiv_pairs.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE TAUIJ_EQUIV_PAIRS
C   ********************************************************************
C   *                                                                  *
C   *   determine set of equivalent cluster site pairs within TAU_ij   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQCLU,IQ_QCLU,NQHOST,QBAS_QCLU
      USE MOD_SYMMETRY,ONLY:NSYM,SYMACCEPTED,SYMUNITARY,ICL_IJQCLU,
     &    IQCLU_REP_IJQCLU,JQCLU_REP_IJQCLU,ISYM_GEN_IJQCLU,NCL_IJQCLU,
     &    NCOPY_IJQCLU,IQ0_IQCLU,NVEC0_IQCLU
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_LATTICE,ONLY:ABAS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_EQUIV_PAIRS')
C
C Local variables
C
      REAL*8 BR(3,3),BRINV(3,3),BRMAT(3,3),QVECP(3),QVECP_IQCLU1(3),
     &       QVECP_JQCLU1(3)
      REAL*8 DDOT
      INTEGER DEL_NVEC_IJ1(3),DEL_NVEC_IJ2(3),I,ILOOP_SYM_TYPE,
     &        IMAT_IJQCLU(:,:),IQ1,IQ2,IQCLU1,IQCLU2,IQEXT,IQEXT1,
     &        IQP_IQCLU1,ISYM,J,JQCLU1,JQCLU2,JQEXT1,JQP_JQCLU1,N1,
     &        NVEC_IQCLU1(3),NVEC_JQCLU1(3)
      LOGICAL KDONE(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KDONE,IMAT_IJQCLU
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (IMAT_IJQCLU(NQCLU,NQCLU),KDONE(NQCLU,NQCLU))
      KDONE(:,:) = .FALSE.
C
      IF ( ALLOCATED(ISYM_GEN_IJQCLU) )
     &     DEALLOCATE (ISYM_GEN_IJQCLU,IQ0_IQCLU,NVEC0_IQCLU,
     &     IQCLU_REP_IJQCLU,JQCLU_REP_IJQCLU,ICL_IJQCLU)
C
      ALLOCATE (ISYM_GEN_IJQCLU(NQCLU,NQCLU))
      ALLOCATE (IQ0_IQCLU(NQCLU),NVEC0_IQCLU(3,NQCLU))
      ALLOCATE (IQCLU_REP_IJQCLU(NQCLU,NQCLU))
      ALLOCATE (JQCLU_REP_IJQCLU(NQCLU,NQCLU))
      ALLOCATE (ICL_IJQCLU(NQCLU,NQCLU))
      ICL_IJQCLU(:,:) = 0
C
      WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
C
C=======================================================================
C                                      BR    basis vectors in real space
      BR(1:3,1:3) = ABAS(1:3,1:3)
C
      DO J = 1,3
         DO I = 1,3
            BRMAT(I,J) = DDOT(3,BR(1,I),1,BR(1,J),1)
         END DO
      END DO
C
      CALL RINVGJ(BRINV,BRMAT,3,3)
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (*,'(''BR''   ,/,3(3f10.5,/))') BR
         WRITE (*,'(''BRinv'',/,3(3f10.5,/))') BRINV
      END IF
C
C=======================================================================
C                  scan site-diagonal site-pairs
C           only IDENTICAL host sites are connected
C           symmetry EQUIVALENT host sites are not yet connected
C=======================================================================
C
      NCL_IJQCLU = 0
      ISYM = 1
C
      LOOP_IQCLU1_DIA:DO IQCLU1 = 1,NQCLU
C
         IQEXT = NQHOST + IQCLU1
C
         ICL_IJQCLU(IQCLU1,IQCLU1) = NCL_IJQCLU
C
         CALL SYM_GET_NVECP(ISYM,BR,BRINV,QBAS_QCLU(1,IQCLU1),QVECP,
     &                      IQ0_IQCLU(IQCLU1),NVEC0_IQCLU(1,IQCLU1))
C
         IQ1 = IQ_QCLU(IQCLU1)
         IF ( IQ0_IQCLU(IQCLU1).NE.IQ1 )
     &         CALL STOP_MESSAGE(ROUTINE,'IQ0_IQCLU(IQCLU1).NE.IQ1')
C
         IF ( IPRINT.GT.0 ) WRITE (6,99004) 'I',IQCLU1,IQEXT,ISYM,
     &                             IQ0_IQCLU(IQCLU1),
     &                             NVEC0_IQCLU(:,IQCLU1)
C
         IF ( KDONE(IQCLU1,IQCLU1) ) CYCLE LOOP_IQCLU1_DIA
C
C-----------------------------------------------------------------------
         LOOP_IQCLU2_DIA:DO IQCLU2 = 1,NQCLU
C
            IF ( .NOT.KDONE(IQCLU2,IQCLU2) ) CYCLE LOOP_IQCLU2_DIA
            IQ2 = IQ_QCLU(IQCLU2)
            IF ( IQ1.NE.IQ2 ) CYCLE LOOP_IQCLU2_DIA
C
C-------------------------------- sites IQCLU1 and IQCLU2 are equivalent
C
            ICL_IJQCLU(IQCLU1,IQCLU1) = ICL_IJQCLU(IQCLU2,IQCLU2)
            ISYM_GEN_IJQCLU(IQCLU1,IQCLU1) = ISYM
            IQCLU_REP_IJQCLU(IQCLU1,IQCLU1) = IQCLU2
            JQCLU_REP_IJQCLU(IQCLU1,IQCLU1) = IQCLU2
            KDONE(IQCLU1,IQCLU1) = .TRUE.
            GOTO 100
C
         END DO LOOP_IQCLU2_DIA
C-----------------------------------------------------------------------
C
         NCL_IJQCLU = NCL_IJQCLU + 1
         ICL_IJQCLU(IQCLU1,IQCLU1) = NCL_IJQCLU
         ISYM_GEN_IJQCLU(IQCLU1,IQCLU1) = ISYM
         IQCLU_REP_IJQCLU(IQCLU1,IQCLU1) = IQCLU1
         JQCLU_REP_IJQCLU(IQCLU1,IQCLU1) = IQCLU1
         KDONE(IQCLU1,IQCLU1) = .TRUE.
C
 100  END DO LOOP_IQCLU1_DIA
C=======================================================================
C
C
C=======================================================================
C                  scan site-OFF-diagonal site-pairs
C=======================================================================
      LOOP_IQCLU1_OFF:DO IQCLU1 = 1,NQCLU
         IQEXT1 = NQHOST + IQCLU1
C
         LOOP_JQCLU1_OFF:DO JQCLU1 = 1,NQCLU
            IF ( IQCLU1.EQ.JQCLU1 ) CYCLE LOOP_JQCLU1_OFF
            JQEXT1 = NQHOST + JQCLU1
C
            IF ( KDONE(IQCLU1,JQCLU1) ) CYCLE LOOP_JQCLU1_OFF
            NCL_IJQCLU = NCL_IJQCLU + 1
            ICL_IJQCLU(IQCLU1,JQCLU1) = NCL_IJQCLU
            ISYM_GEN_IJQCLU(IQCLU1,JQCLU1) = 1
            IQCLU_REP_IJQCLU(IQCLU1,JQCLU1) = IQCLU1
            JQCLU_REP_IJQCLU(IQCLU1,JQCLU1) = JQCLU1
            KDONE(IQCLU1,JQCLU1) = .TRUE.
C
C-----------------------------------------------------------------------
C          run twice over all symmetry operations
C          ILOOP_SYM_TYPE = 1: consider ONLY      UNITARY operations
C          ILOOP_SYM_TYPE = 2: consider ONLY ANTI-UNITARY operations
C-----------------------------------------------------------------------
C
            LOOP_SYM_TYPE:DO ILOOP_SYM_TYPE = 1,2
               LOOP_ISYM_OFF:DO ISYM = 1,NSYM
                  IF ( .NOT.SYMACCEPTED(ISYM) ) CYCLE LOOP_ISYM_OFF
                  IF ( ILOOP_SYM_TYPE.EQ.1 ) THEN
                     IF ( .NOT.SYMUNITARY(ISYM) ) CYCLE LOOP_ISYM_OFF
                  ELSE IF ( SYMUNITARY(ISYM) ) THEN
                     CYCLE LOOP_ISYM_OFF
                  END IF
C
                  CALL SYM_GET_NVECP(ISYM,BR,BRINV,QBAS_QCLU(1,IQCLU1),
     &                               QVECP_IQCLU1,IQP_IQCLU1,
     &                               NVEC_IQCLU1)
C
                  CALL SYM_GET_NVECP(ISYM,BR,BRINV,QBAS_QCLU(1,JQCLU1),
     &                               QVECP_JQCLU1,JQP_JQCLU1,
     &                               NVEC_JQCLU1)
C
                  DEL_NVEC_IJ1(:) = NVEC_IQCLU1(:) - NVEC_JQCLU1(:)
C
                  IF ( IPRINT.GT.0 ) THEN
                     WRITE (6,99004) 'I',IQCLU1,IQEXT1,ISYM,IQP_IQCLU1,
     &                               NVEC_IQCLU1
                     WRITE (6,99004) 'J',JQCLU1,JQEXT1,ISYM,JQP_JQCLU1,
     &                               NVEC_JQCLU1
                  END IF
C
C-----------------------------------------------------------------------
                  LOOP_IQCLU2_OFF:DO IQCLU2 = 1,NQCLU
                     LOOP_JQCLU2_OFF:DO JQCLU2 = 1,NQCLU
                        IF ( IQCLU2.EQ.JQCLU2 ) CYCLE LOOP_JQCLU2_OFF
                        IF ( KDONE(IQCLU2,JQCLU2) )
     &                       CYCLE LOOP_JQCLU2_OFF
C
                        IF ( IQP_IQCLU1.EQ.IQ0_IQCLU(IQCLU2) ) THEN
                           IF ( JQP_JQCLU1.EQ.IQ0_IQCLU(JQCLU2) ) THEN
                              DEL_NVEC_IJ2(:) = NVEC0_IQCLU(:,IQCLU2)
     &                           - NVEC0_IQCLU(:,JQCLU2)
                              DO I = 1,3
                                 IF ( DEL_NVEC_IJ1(I).NE.DEL_NVEC_IJ2(I)
     &                                ) GOTO 102
                              END DO
C
C------------------------- pair (1) and pair (2) of sites are equivalent
C
                              ICL_IJQCLU(IQCLU2,JQCLU2)
     &                           = ICL_IJQCLU(IQCLU1,JQCLU1)
                              ISYM_GEN_IJQCLU(IQCLU2,JQCLU2) = ISYM
                              IQCLU_REP_IJQCLU(IQCLU2,JQCLU2) = IQCLU1
                              JQCLU_REP_IJQCLU(IQCLU2,JQCLU2) = JQCLU1
                              KDONE(IQCLU2,JQCLU2) = .TRUE.
C
                           END IF
                        END IF
C
 102                 END DO LOOP_JQCLU2_OFF
                  END DO LOOP_IQCLU2_OFF
C-----------------------------------------------------------------------
C
               END DO LOOP_ISYM_OFF
            END DO LOOP_SYM_TYPE
C
         END DO LOOP_JQCLU1_OFF
      END DO LOOP_IQCLU1_OFF
C=======================================================================
C
      N1 = 0
      NCOPY_IJQCLU = 0
      DO IQCLU2 = 1,NQCLU
         DO JQCLU2 = 1,NQCLU
            IF ( IPRINT.GT.0 ) WRITE (6,99003) IQCLU2,JQCLU2,
     &                                ICL_IJQCLU(IQCLU2,JQCLU2),
     &                                ISYM_GEN_IJQCLU(IQCLU2,JQCLU2),
     &                                IQCLU_REP_IJQCLU(IQCLU2,JQCLU2),
     &                                JQCLU_REP_IJQCLU(IQCLU2,JQCLU2)
C
            IF ( ISYM_GEN_IJQCLU(IQCLU2,JQCLU2).EQ.1 ) THEN
C
               IF ( IQCLU_REP_IJQCLU(IQCLU2,JQCLU2).EQ.IQCLU2 .AND. 
     &              JQCLU_REP_IJQCLU(IQCLU2,JQCLU2).EQ.JQCLU2 ) THEN
C
                  N1 = N1 + 1
C
               ELSE
C
                  NCOPY_IJQCLU = NCOPY_IJQCLU + 1
C
               END IF
            END IF
         END DO
      END DO
C
      IF ( N1.NE.NCL_IJQCLU )
     &      CALL STOP_MESSAGE(ROUTINE,'IQ0_IQCLU(IQCLU1).NE.IQ1')
C
      WRITE (6,99002) NQCLU*NQCLU,NCL_IJQCLU,NCOPY_IJQCLU,
     &                NQCLU*NQCLU - NCL_IJQCLU - NCOPY_IJQCLU
C
      IMAT_IJQCLU(:,:) = 0
      DO IQCLU2 = 1,NQCLU
         DO JQCLU2 = 1,NQCLU
            IF ( IQCLU_REP_IJQCLU(IQCLU2,JQCLU2).EQ.IQCLU2 .AND. 
     &           JQCLU_REP_IJQCLU(IQCLU2,JQCLU2).EQ.JQCLU2 )
     &           IMAT_IJQCLU(IQCLU2,JQCLU2) = ICL_IJQCLU(IQCLU2,JQCLU2)
     &           + 1
         END DO
      END DO
C
      CALL IMATSTRUCT('generating site pairs',IMAT_IJQCLU,NQCLU,NQCLU,0,
     &                0,0,-6)
C
      IMAT_IJQCLU(:,:) = ICL_IJQCLU(:,:) + 1
      CALL IMATSTRUCT('classes ICL_IJQCLU of site pairs',IMAT_IJQCLU,
     &                NQCLU,NQCLU,0,0,0,-6)
C
      CALL IMATSTRUCT(
     &    'generating symmetry operation ISYM_GEN_IJQCLU for site pairs'
     &    ,ISYM_GEN_IJQCLU,NQCLU,NQCLU,0,0,0,6)
C
C=======================================================================
99001 FORMAT (//,1X,79('*'),/,35X,'<',A,'>',/,1X,79('*'),//,10X,
     &      'determine classes of symmetry connected cluster site pairs'
     &      ,/)
99002 FORMAT (/,10X,'total number of site pairs           NQCLU*NQCLU: '
     &        ,I6,/,10X,
     &        'inequivalent site pairs (classes)    NCL_IJQCLU:  ',I6,/,
     &        10X,'site pairs generated by copying      NCOPY_IJQCLU:',
     &        I6,/,10X,
     &        'site pairs generated by other symmetry operations:',I6,/)
99003 FORMAT (' I J : ',2I3,'  ICL:',I3,'   ISYM:',I3,' I0 J0:',2I3)
99004 FORMAT (10X,'SYM ',A,10I3)
      END
C*==sym_get_nvecp.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SYM_GET_NVECP(ISYM,BR,BRINV,QVEC,QVECP,IQP,NVECP)
C   ********************************************************************
C   *                                                                  *
C   *  general symmetry operation:                                     *
C   *                                                                  *
C   *        {D|->t} ->q = D ->q + ->t = R_q + ->q''                   *
C   *                                                                  *
C   * find  R_q  as NSFTSYMQ that ensures ->q'' to be within unit cell *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:MROTR,SYMTVEC
      USE MOD_SITES,ONLY:NQHOST,QBAS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYM_GET_NVECP')
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
C
C Dummy arguments
C
      INTEGER IQP,ISYM
      REAL*8 BR(3,3),BRINV(3,3),QVEC(3),QVECP(3)
      INTEGER NVECP(3)
C
C Local variables
C
      REAL*8 BV(3),CF(3),DVEC(3),QVECTST(3)
      REAL*8 DDOT
      LOGICAL FOUND
      INTEGER I,IQHOST
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      QVECP(:) = MATMUL(MROTR(:,:,ISYM),QVEC(:)) + SYMTVEC(:,ISYM)
C
      FOUND = .FALSE.
C
      DO IQHOST = 1,NQHOST
C
         DVEC(:) = QVECP(:) - QBAS(:,IQHOST)
C
         DO I = 1,3
            BV(I) = DDOT(3,BR(1,I),1,DVEC,1)
         END DO
C
         CALL DGEMV('N',3,3,1D0,BRINV,3,BV,1,0D0,CF,1)
C
         DVEC(:) = 0D0
         DO I = 1,3
            NVECP(I) = NINT(CF(I))
            IF ( ABS(CF(I)-NVECP(I)).GT.TOL ) GOTO 100
            DVEC(:) = DVEC(:) + CF(I)*BR(:,I)
         END DO
C
         FOUND = .TRUE.
         IQP = IQHOST
         QVECTST(:) = QBAS(:,IQP) + DVEC(:)
C
         EXIT
C
 100  END DO
C
      IF ( .NOT.FOUND ) CALL STOP_MESSAGE(ROUTINE,'IQHOST not found')
C
      IF ( .NOT.RVEC_SAME(3,QVECTST,QVECP,1D-8) ) THEN
         WRITE (6,99001) IQP,ISYM,QVECTST,QVECP
         CALL STOP_MESSAGE(ROUTINE,'TEST not passed for given TOL')
      END IF
C
99001 FORMAT ('##### IQ=',I4,' ISYM=',I3,' QVECTST = ',3F10.6,/22X,
     &        ' QVECP   = ',3F10.6,/)
C
      END
