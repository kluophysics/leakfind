C*==scfmad0d_back.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD0D_BACK
C   ********************************************************************
C   *                                                                  *
C   *   set up the Madelung potential parameters   VLMMAD_BACK         *
C   *   for the host for background an embedded cluster                *
C   *                                                                  *
C   *   RELAX_CLU=.T.   the cluster positions are shifted              *
C   *                   transform  VLMMAD_BACK  to new positions       *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      USE MOD_SITES,ONLY:AVMAD,NQHOST,NQCLU,IQ_QCLU,NLVMAD,NLMVMAD,
     &    NLQMAD,NLMQMAD,NLMAD,QBAS,QBAS_QCLU,QBAS0_QCLU,VLMMAD_BACK,
     &    CMNTQ,VLMMAD_HOST
      USE MOD_CALCMODE,ONLY:RELAX_CLU
      USE MOD_ANGMOM,ONLY:L_LM,NRGNT123TAB,DBLFACT_L
      IMPLICIT NONE
C*--SCFMAD0D_BACK20
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 DNRM2,GAUNT_RYLM
      INTEGER I,IA_ERR,IQCLU,IQ_HOST,IQ_IMP,IRGNTMAD(:,:),JQCLU,JQ_HOST,
     &        JQ_IMP,L,L0,L1,L2,L2MAX,LM0,LM1,LM2,LMQ,LMV,M0,M1,M2,
     &        NLM2MAX,NRGNTMAD,NRGNTMADMAX
      REAL*8 M1_PW_L(:),RGAUNT,RGNTMAD(:),RYLM(:),SABS,SHAT(3),SVEC(3),
     &       VLMMAD_BACK0(:,:),WGT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IRGNTMAD,RGNTMAD,M1_PW_L,RYLM,VLMMAD_BACK0
C
      IF ( NLMQMAD.LT.NLMVMAD ) STOP '<SCFMAD0D>: NLMQMAD < NLMVMAD'
C
      WRITE (6,99003) NLVMAD,NLMVMAD,NLQMAD,NLMQMAD
C
C-----------------------------------------------------------------------
C     set cluster site position to unrelaxed positions in QBAS and
C     calculate corresponding real space cluster Madelung matrix
C-----------------------------------------------------------------------
C
      QBAS(1:3,NQHOST+1:NQHOST+NQCLU) = QBAS0_QCLU(1:3,1:NQCLU)
C
      CALL SCFMAD0D
C
C-----------------------------------------------------------------------
C     reset cluster site positions in QBAS to original positions
C-----------------------------------------------------------------------
C
      QBAS(1:3,NQHOST+1:NQHOST+NQCLU) = QBAS_QCLU(1:3,1:NQCLU)
C
C-----------------------------------------------------------------------
C           calculate background Madelung potential
C                VLMMAD_BACK of the host system
C-----------------------------------------------------------------------
C
      ALLOCATE (VLMMAD_BACK(NLMVMAD,(NQHOST+1):(NQHOST+NQCLU)))
      VLMMAD_BACK(1:NLMVMAD,(NQHOST+1):(NQHOST+NQCLU)) = 0D0
C
      DO IQCLU = 1,NQCLU
C
         IQ_IMP = NQHOST + IQCLU
         IQ_HOST = IQ_QCLU(IQCLU)
C
         DO LMV = 1,NLMVMAD
C
            VLMMAD_BACK(LMV,IQ_IMP) = VLMMAD_HOST(LMV,IQ_HOST)
C
            DO JQCLU = 1,NQCLU
C
               JQ_IMP = NQHOST + JQCLU
               JQ_HOST = IQ_QCLU(JQCLU)
C
               DO LMQ = 1,NLMQMAD
C
                  VLMMAD_BACK(LMV,IQ_IMP) = VLMMAD_BACK(LMV,IQ_IMP)
     &               - AVMAD(IQ_IMP,JQ_IMP,LMV,LMQ)*CMNTQ(LMQ,JQ_HOST)
C
               END DO
            END DO
C
         END DO
      END DO
C
      IF ( IPRINT.GE.0 ) THEN
         LMV = 1
         WRITE (6,99002)
         DO IQCLU = 1,NQCLU
C
            IQ_IMP = NQHOST + IQCLU
            IQ_HOST = IQ_QCLU(IQCLU)
C
            WRITE (6,99001) IQCLU,IQ_HOST,CMNTQ(LMV,IQ_HOST),
     &                      VLMMAD_HOST(LMV,IQ_HOST),
     &                      VLMMAD_BACK(LMV,IQ_IMP)
         END DO
         WRITE (6,*) ' '
      END IF
C
C=======================================================================
C                         relaxed embedded cluster
C=======================================================================
C
      IF ( .NOT.RELAX_CLU ) RETURN
C
      ALLOCATE (M1_PW_L(0:4*NLMAD))
C
      M1_PW_L(0) = 1D0
      DO L = 1,4*NLMAD
         M1_PW_L(L) = -M1_PW_L(L-1)
      END DO
C
C-----------------------------------------------------------------------
C              calculation of the harmonic polynomials
C-----------------------------------------------------------------------
C
      L2MAX = 2*(NLMAD-1)
      NLM2MAX = (L2MAX+1)**2
C
      ALLOCATE (RYLM(NLM2MAX))
C
C-----------------------------------------------------------------------
C     set up of the gaunt coefficients with an index field
C     recognize that they are needed here only for L2=L0-L1
C-----------------------------------------------------------------------
C
      NRGNTMADMAX = 2*NRGNT123TAB(NLMAD)
      ALLOCATE (IRGNTMAD(NRGNTMADMAX,3))
      ALLOCATE (RGNTMAD(NRGNTMADMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD0D->DFAC'
C
      I = 1
      DO L0 = 0,(NLVMAD-1)
         DO L1 = 0,(NLVMAD-1)
            L2 = L0 - L1
            DO M0 = -L0,L0
               DO M1 = -L1,L1
                  DO M2 = -L2,L2
C
                     RGAUNT = GAUNT_RYLM(L0,M0,L1,M1,L2,M2)
C
                     IF ( ABS(RGAUNT).GT.1.D-10 ) THEN
                        IF ( I.GT.NRGNTMADMAX ) THEN
                           WRITE (6,99004) I,NRGNTMADMAX
                           STOP 'in <SCFMAD0D_BACK> '
                        END IF
C
                        RGNTMAD(I) = RGAUNT
                        IRGNTMAD(I,1) = L0*(L0+1) + M0 + 1
                        IRGNTMAD(I,2) = L1*(L1+1) + M1 + 1
                        IRGNTMAD(I,3) = L2*(L2+1) + M2 + 1
                        I = I + 1
                     END IF
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      NRGNTMAD = I - 1
C
      ALLOCATE (VLMMAD_BACK0(NLMVMAD,(NQHOST+1):(NQHOST+NQCLU)))
C
      VLMMAD_BACK0(1:NLMVMAD,(NQHOST+1):(NQHOST+NQCLU))
     &   = VLMMAD_BACK(1:NLMVMAD,(NQHOST+1):(NQHOST+NQCLU))
C
      VLMMAD_BACK(1:NLMVMAD,(NQHOST+1):(NQHOST+NQCLU)) = 0D0
C
      DO IQCLU = 1,NQCLU
C
         IQ_IMP = NQHOST + IQCLU
         IQ_HOST = IQ_QCLU(IQCLU)
C
         SVEC(1:3) = QBAS_QCLU(1:3,IQCLU) - QBAS0_QCLU(1:3,IQCLU)
C
         SABS = DNRM2(3,SVEC,1)
C
         IF ( SABS.LT.1D-6 ) THEN
C
            VLMMAD_BACK(1:NLMVMAD,IQ_IMP)
     &         = VLMMAD_BACK0(1:NLMVMAD,IQ_IMP)
C
         ELSE
C
C--------------------------- normalize vector to get spherical harmonics
C
            SHAT(1:3) = SVEC(1:3)/SABS
C
            CALL CALC_RHPLM(SHAT(1),SHAT(2),SHAT(3),RYLM,L2MAX,NLM2MAX)
C
C---> this loop has to be calculated only for   L2=L0-L1
C
            DO I = 1,NRGNTMAD
               LM0 = IRGNTMAD(I,1)
               LM1 = IRGNTMAD(I,2)
               LM2 = IRGNTMAD(I,3)
               L0 = L_LM(LM0)
               L1 = L_LM(LM1)
               L2 = L0 - L1
C
               WGT = CONST_4PI*M1_PW_L(L1+L0)*RGNTMAD(I)
C
               WGT = WGT*DBLFACT_L(L0)/(DBLFACT_L(L1)*DBLFACT_L(L2))
C
               WGT = WGT*SABS**L2*RYLM(LM2)
C
               VLMMAD_BACK(LM1,IQ_IMP) = VLMMAD_BACK(LM1,IQ_IMP)
     &            + WGT*VLMMAD_BACK0(LM0,IQ_IMP)
C
            END DO
C
         END IF
C
      END DO
C
      IF ( IPRINT.GE.0 ) THEN
         LMV = 1
         WRITE (6,99005)
         DO IQCLU = 1,NQCLU
C
            IQ_IMP = NQHOST + IQCLU
            IQ_HOST = IQ_QCLU(IQCLU)
C
            WRITE (6,99006) IQCLU,IQ_HOST,CMNTQ(LMV,IQ_HOST),
     &                      VLMMAD_BACK(LMV,IQ_IMP)
         END DO
         WRITE (6,*) ' '
      END IF
C
99001 FORMAT (10X,2I5,3F15.6)
99002 FORMAT (/,10X,'host Madelung parameters (L=(0,0)) ',
     &        'for unrelaxed posistions',/,/,12X,'IQCLU IQ_HOST',
     &        '    CMNTQ      VLMMAD_HOST    VLMMAD_BACK',/)
99003 FORMAT (//,1X,79('*'),/,32X,'<SCFMAD0D_BACK>',/,1X,79('*'),//,10X,
     &        'set up Madelung parameters for the host ',/,10X,
     &        'for embedded cluster calculations',//,10X,'NLVMAD =',I4,
     &        5X,'NLMVMAD =',I4,/,10X,'NLQMAD =',I4,5X,'NLMQMAD =',I4,/)
99004 FORMAT ('<SCFMAD0D_BACK>: I=',I5,' > NRGNTMADMAX =',I5)
99005 FORMAT (/,10X,'host Madelung parameters (L=(0,0)) ',
     &        'for relaxed posistions',/,/,12X,'IQCLU IQ_HOST',
     &        '    CMNTQ ',20X,'VLMMAD_BACK',/)
99006 FORMAT (10X,2I5,F15.6,15X,F15.6)
      END
