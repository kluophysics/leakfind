C*==fpscfnewpot.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPSCFNEWPOT(ETOT,ITRSCF,R2RHO)
C   ********************************************************************
C   *                                                                  *
C   * This subroutine calculates the effective potential and magnetic  *
C   * field for the next iteration of the scf-cycle.                   *
C   *                                                                  *
C   * FPVINTRA  calculates the intracell-part of the Hartree-potential *
C   * SCFMAD_POTcalculates the intercell-part of the Hartree-potential *
C   * FPVXCLM   calculates the spin-up and spin-down                   *
C   *           exchange-correlation-potential                         *
C   *                                                                  *
C   * R2RHO(...,1) means the charge-density ( spin-up  +  spin-down  ) *
C   * R2RHO(...,2) means the spin-density   ( spin-up  -  spin-down  ) *
C   *                                                                  *
C   * KTE=(0,1):       key for calculating total energies              *
C   * 0       no calculation                                           *
C   * 1       calculation in the each iteration                        *
C   *                                                                  *
C   * CALC_EFG         key for calculating electric field gradients    *
C   * F       no calculation (DEFAULT)                                 *
C   * T       calculation at the end of the SCF-cycle (RETURN)         *
C   *                                                                  *
C   * KVMAD :          only invoked in ECOUB                           *
C   *                                                                  *
C   * SCFSIM=0..1      parameter for treating charge-transfer-effects  *
C   * in the screened-impurity-model of the Madelung-potential         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYM,IQREPQ,ISYMGENQ,SYMEULANG,
     &    SYMMETRIZE_RHO,SYMMETRIZE_POT
      USE MOD_FILES,ONLY:IPRINT,IOTMP
      USE MOD_RMESH,ONLY:NM,NRMAX,NSFMAX,NRSFMAX,FLMSF,R,JRNS1,ISFLM,
     &    LMISF,KLMSF,NSF,JRCRI,JRCUT,NPAN,NMAX_LEBGRID,RHAT_LEBGRID,
     &    N_LEBGRID,W_LEBGRID,Y_LEBGRID
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      USE MOD_SITES,ONLY:RNNQ,NQ,NQMAX,ITOQ,NOQ,IQBOT,IQTOP,QMTET,QMPHI,
     &    QMGAM,MAGROT_Q,NQCLU,FORCE_QCLU,NQHOST,NLMQMAD,NLQMAD,CMNTQ
      USE MOD_TYPES,ONLY:NLMFPMAX,NLFPMAX,NTMAX,CONC,BNST,VNST,BT,VT,
     &    NLMFPT,KLMFP,NLFP,NLMFP,Z,IMT,ITBOT,ITTOP,VMTZ,BEXT
      USE MOD_CPA,ONLY:NCPA
      USE MOD_CALCMODE,ONLY:IREL,BREAKPOINT,CALC_KINETIC_ENERGY_DENSITY,
     &    CALC_EFG,CALC_FORCE
      USE MOD_CONSTANTS,ONLY:PI,SQRT_4PI
      USE MOD_SCF,ONLY:SCFSIM,NFCOR
      IMPLICIT NONE
C*--FPSCFNEWPOT49
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPSCFNEWPOT')
      LOGICAL USE_KLMFP,WRITE_IS_POT
      PARAMETER (USE_KLMFP=.TRUE.,WRITE_IS_POT=.FALSE.)
C
C Dummy arguments
C
      REAL*8 ETOT
      INTEGER ITRSCF
      REAL*8 R2RHO(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      REAL*8 BLMTMP(:,:),BNEW(:,:,:),CMNTISQ(:,:),CMNTIST(:,:),
     &       CMNTMTQ(:,:),CMNTMTT(:,:),DROT_QLM(:,:,:),DSYM_RLM(:,:,:),
     &       D_INV(:,:),D_LMLMP,D_ONE(:,:),D_WRK(:,:),ECOU(:,:),
     &       EPOTIN(:),EXC(:,:),FAUX(:,:,:),FHF_LMQ(:,:),FHF_LMT(:,:),
     &       F_LMQ(:,:),F_LMT(:,:),QNETT(:),RGNT_VSF(:),RWGT,SUMA,SUMB,
     &       VLMTMP(:,:),VNEW(:,:,:),VSIM,WRLM1(:,:),WRLM2(:,:),WTMPA(:)
     &       ,WTMPB(:),Y00,ZZOR
      INTEGER I,IA_ERR,IFLAG,IFLAG_NEUTRALITY,IINV,IM,IO,IQ,IQCLU,IQREP,
     &        IR,IRCRIT,IRMTIN,IRSF,ISF,ISPIN,ISYM,ISYMP,IT,J,KTE,KVMAD,
     &        L,LM,LM1,LM2,LM3,LMP,LMRGNT_VSF(:,:),M,NRGNT_VSF,
     &        NRGNT_VSF_LM(:),NSPIN,NSYMH
      LOGICAL Q1M_OCCURS
      CHARACTER*3 STRDNS(3)
      SAVE DROT_QLM,DSYM_RLM,LMRGNT_VSF,NRGNT_VSF,NRGNT_VSF_LM,
     &     Q1M_OCCURS,RGNT_VSF
C
C*** End of declarations rewritten by SPAG
C
      DATA STRDNS/'chr','spn','orb'/
C
      ALLOCATABLE LMRGNT_VSF,NRGNT_VSF_LM,WRLM1,WRLM2
      ALLOCATABLE CMNTIST,CMNTISQ,CMNTMTT,CMNTMTQ
      ALLOCATABLE DROT_QLM,D_INV,DSYM_RLM
      ALLOCATABLE ECOU,EPOTIN,EXC,RGNT_VSF,QNETT,D_ONE,D_WRK,BNEW,VNEW
      ALLOCATABLE VLMTMP,BLMTMP
      ALLOCATABLE WTMPB,WTMPA,FHF_LMT,F_LMT,FHF_LMQ,F_LMQ,FAUX
C
      CALL TRACK_INFO(ROUTINE)
C
      IFLAG_NEUTRALITY = 0
C
      KTE = 1
      KVMAD = 0
      IF ( IREL.GE.2 ) THEN
         NSPIN = 2
      ELSE
         NSPIN = 1
      END IF
C
C ------------------- NOTE: the same  LM-expansion is used for ALL sites
C
      IF ( NLFP.GT.NLFPMAX ) CALL STOP_MESSAGE(ROUTINE,'NLFP > NLFPMAX')
      IF ( NLFP.NE.NLQMAD ) CALL STOP_MESSAGE(ROUTINE,'NLFP.NE.NLQMAD')
C
      WRITE (6,99004)
C
      ALLOCATE (FAUX(NRMAX,NLMFPMAX,2))
      ALLOCATE (FHF_LMQ(2:4,NQMAX),F_LMQ(2:4,NQMAX))
      ALLOCATE (FHF_LMT(2:4,NTMAX),F_LMT(2:4,NTMAX))
      ALLOCATE (WTMPB(NLMFPMAX),WTMPA(NLMFPMAX))
      ALLOCATE (CMNTIST(NLMFPMAX,NTMAX),CMNTISQ(NLMFPMAX,NQMAX))
      ALLOCATE (CMNTMTT(NLMFPMAX,NTMAX),CMNTMTQ(NLMFPMAX,NQMAX))
      ALLOCATE (EXC(0:(NLFPMAX-1),NTMAX),QNETT(NTMAX))
      ALLOCATE (VNEW(NRMAX,NLMFPMAX,NTMAX),BNEW(NRMAX,NLMFPMAX,NTMAX))
      ALLOCATE (VLMTMP(NRSFMAX,NLMFPMAX),BLMTMP(NRSFMAX,NLMFPMAX))
      ALLOCATE (ECOU(0:(NLFPMAX-1),NTMAX),EPOTIN(NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate VNS')
C
      CALL RINIT(NRMAX*NLMFPMAX*NTMAX,BNEW)
C
C   ********************************************************************
C                       INITIALISATION - START
C   ********************************************************************
      IF ( ITRSCF.LE.1 ) THEN
C
         IF ( ALLOCATED(DSYM_RLM) ) THEN
            DEALLOCATE (DROT_QLM,RGNT_VSF,LMRGNT_VSF)
            DEALLOCATE (DSYM_RLM,NRGNT_VSF_LM)
         END IF
C
C ------- initialize local variables  RGNT_VSF, NRGNT_VSF_LM, LMRGNT_VSF
C
         ALLOCATE (NRGNT_VSF_LM(0:NLMFPMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'allocate RGNT_VSF')
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
         CALL RGNTSF_SETUP(NLFP-1,NLFP-1,NRGNT_VSF_LM,NRGNT_VSF)
C
         ALLOCATE (RGNT_VSF(NRGNT_VSF),LMRGNT_VSF(NRGNT_VSF,3))
C
         REWIND IOTMP
         DO I = 1,NRGNT_VSF
            READ (IOTMP) (LMRGNT_VSF(I,J),J=1,3),RGNT_VSF(I)
         END DO
         CLOSE (IOTMP)
C
C=======================================================================
C                   check consistency of shape functions
C=======================================================================
C
C-----------------------------------------------------------------------
C                    check consistency of KLMSF and LMISF
C-----------------------------------------------------------------------
C
         DO IM = 1,NM
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               IF ( KLMSF(LM,IM).NE.1 ) THEN
                  WRITE (6,*) '   KLMSF <> 1 for IM=',IM,' LM =',LM
                  CALL STOP_MESSAGE(ROUTINE,
     &                              'KLMSF and LMISF not consistent')
               END IF
            END DO
         END DO
C
C-----------------------------------------------------------------------
C        convolution of any function with the shape functions
C        must not lead to any new and symmetry forbidden LM-terms
C-----------------------------------------------------------------------
C
         IFLAG = 0
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
C
            DO J = 1,NRGNT_VSF_LM(NLMFP)
               LM1 = LMRGNT_VSF(J,1)
               IF ( LM1.GT.NLMFPMAX )
     &               CALL STOP_MESSAGE(ROUTINE,'LMRGNT_VSF > NLMFPMAX')
               LM2 = LMRGNT_VSF(J,2)
               LM3 = LMRGNT_VSF(J,3)
C-----------------------------------------------------------------------
               IF ( KLMSF(LM3,IM).EQ.1 .AND. KLMFP(LM2,IT).EQ.1 .AND. 
     &              KLMFP(LM2,IT).NE.1 ) THEN
                  IFLAG = 1
                  ISF = ISFLM(LM3,IM)
                  WRITE (6,99001) IT,ISF,LM3,LM1,LM2
C ----------------------------------------------------------------------
               END IF
            END DO
C ----------------------------------------------------------------------
C
         END DO
         IF ( IFLAG.EQ.1 )
     &         CALL STOP_MESSAGE(ROUTINE,'check for SFNs failed')
C
C=======================================================================
C              rotation matrices for REAL spherical harmonics
C=======================================================================
C
         ALLOCATE (DROT_QLM(NLMQMAD,NLMQMAD,IQBOT:IQTOP))
         ALLOCATE (D_ONE(NLMQMAD,NLMQMAD))
         ALLOCATE (D_WRK(NLMQMAD,NLMQMAD),D_INV(NLMQMAD,NLMQMAD))
C
C-----------------------------------------------------------------------
C                create 1-matrix and matrix for inversion
C-----------------------------------------------------------------------
C
         D_INV(1:NLMQMAD,1:NLMQMAD) = 0D0
         D_ONE(1:NLMQMAD,1:NLMQMAD) = 0D0
         I = 0
         DO L = 0,(NLQMAD-1)
            DO M = 1,2*L + 1
               I = I + 1
               D_INV(I,I) = (-1.0D0)**L
               D_ONE(I,I) = 1D0
            END DO
         END DO
C
         NSYMH = NSYM/2
C
         DO IQ = IQBOT,IQTOP
C
            IF ( IQREPQ(IQ).NE.IQ ) THEN
C
C---- connection of NON-representative to associated representative site
C---------------------------------------------------------- IQ <> IQREPQ
C
               IF ( IQREPQ(IQ).GE.IQ )
     &               CALL STOP_MESSAGE(ROUTINE,'IQREP >= IQ')
               ISYM = ISYMGENQ(IQ)
               ISYMP = ISYM - NSYMH
C
               IF ( ISYM.LE.NSYMH ) THEN
C
                  CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYM),
     &                             SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                             DROT_QLM(1,1,IQ))
C
               ELSE
C
                  CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYMP),
     &                             SYMEULANG(2,ISYMP),SYMEULANG(3,ISYMP)
     &                             ,D_WRK)
C
                  DROT_QLM(1:NLMQMAD,1:NLMQMAD,IQ)
     &               = MATMUL(D_INV(1:NLMQMAD,1:NLMQMAD),
     &               D_WRK(1:NLMQMAD,1:NLMQMAD))
C
               END IF
C--------------------------------------------- representative site IQREP
C----------------------------------------------------------- IQ = IQREPQ
C
            ELSE IF ( MAGROT_Q(IQ) ) THEN
C
               CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,QMPHI(IQ),QMTET(IQ),
     &                          QMGAM(IQ),DROT_QLM(1,1,IQ))
C
            ELSE
C
               DROT_QLM(1:NLMQMAD,1:NLMQMAD,IQ)
     &            = D_ONE(1:NLMQMAD,1:NLMQMAD)
C
            END IF
C---------------------------------------------------------------- IQREPQ
C
         END DO
C
C=======================================================================
C
         ALLOCATE (DSYM_RLM(NLMQMAD,NLMQMAD,NSYM))
C
         DSYM_RLM(:,:,:) = 0.0D0
C
         NSYMH = NSYM/2
C
         IINV = NSYMH + 1
C
         DSYM_RLM(1:NLMQMAD,1:NLMQMAD,1) = D_ONE(1:NLMQMAD,1:NLMQMAD)
C
         DSYM_RLM(1:NLMQMAD,1:NLMQMAD,IINV) = D_INV(1:NLMQMAD,1:NLMQMAD)
C
         DO ISYM = 2,NSYMH
C
            CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYM),
     &                       SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                       DSYM_RLM(1,1,ISYM))
C
C-----------------------------------------------------------------------
C                          include inversion
C-----------------------------------------------------------------------
C
            ISYMP = NSYMH + ISYM
C
            CALL DGEMM('N','N',NLMQMAD,NLMQMAD,NLMQMAD,1D0,
     &                 DSYM_RLM(1,1,ISYM),NLMQMAD,DSYM_RLM(1,1,IINV),
     &                 NLMQMAD,0D0,DSYM_RLM(1,1,ISYMP),NLMQMAD)
C
         END DO
C
         Q1M_OCCURS = .FALSE.
         DO IT = ITBOT,ITTOP
            DO LM = 2,4
               IF ( KLMFP(LM,IT).NE.0 ) Q1M_OCCURS = .TRUE.
            END DO
         END DO
C
         DEALLOCATE (D_ONE,D_WRK,D_INV)
C
      END IF
C
C **********************************************************************
C                       INITALISATION - END
C **********************************************************************
C
      ALLOCATE (WRLM1(NLMQMAD,NLMQMAD),WRLM2(NLMQMAD,NLMQMAD))
C
C-----------------------------------------------------------------------
C                         symmetrize RHO
C-----------------------------------------------------------------------
C
      IF ( SYMMETRIZE_RHO ) THEN
         DO ISPIN = 1,NSPIN
            CALL SYMMETRIZE_FLMT_REP(STRDNS(ISPIN),R2RHO(1,1,1,ISPIN),
     &                               FAUX,DSYM_RLM,DROT_QLM,WRLM1,WRLM2)
         END DO
      END IF
C
C  ---------------------------------------------------------------------
C   calculate charge moments   CMNTMTT  w.r.t. muffin tin sphere
C                              CMNTIST  w.r.t. interstitial regime
C                              for each atomic type IT
C   and intra site Coulomb potential  VNEW
C   using the charge density          R2RHO
C   USE_KLMFP=.T.:   restrict calculation to symmetry allowed terms
C  ---------------------------------------------------------------------
C
      CALL FPVINTRA(USE_KLMFP,CMNTMTT,CMNTIST,R2RHO,VNEW,LMRGNT_VSF,
     &              NRGNT_VSF_LM,RGNT_VSF,NRGNT_VSF)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      IF ( SYMMETRIZE_POT ) THEN
CC
C         CALL SYMMETRIZE_FLMT_REP('V  ',VNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C         CALL SYMMETRIZE_FLMT_REP('B  ',BNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C      END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  ---------------------------------------------------------------------
C  Calculate charge-moments  CMNTMTQ  w.r.t. muffin tin sphere
C                            CMNTISQ  w.r.t. interstitial regime
C                            CMNTQ    total
C                            for each site IQ
C  In case of  CPA  this is the concentration-average of the
C  atom-type-dependent charge-moments  CMNTMTT, CMNTIST
C  otherwise:  CMNTMTQ(IQ) = CMNTMTT(IT) and CMNTISQ(IQ) = CMNTIST(IT)
C
C  NOTE:  CMNTQ includes the nuclear contribution -Z/SQRT(4*pi)
C               used to set the up Madelung potential in <SCFMAD_POT>
C  ---------------------------------------------------------------------
C
      CALL SCFCMNT(IPRINT,Q1M_OCCURS,IFLAG_NEUTRALITY,CMNTIST,CMNTISQ,
     &             CMNTMTT,CMNTMTQ,CMNTQ,QNETT,DROT_QLM,DSYM_RLM)
C
C-----------------------------------------------------------------------
C       perform the correction of the Madelung-potential due to
C      charge-transfer according to the screened-impurity-model
C-----------------------------------------------------------------------
C
      SCFSIM = 0D0
      IF ( NCPA.NE.0 .AND. ABS(SCFSIM).GT.1D-6 ) THEN
C
         DO IQ = IQBOT,IQTOP
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               IM = IMT(IT)
               VSIM = -SCFSIM*2.0D0*SQRT(4.0D0*PI)*QNETT(IT)/RNNQ(IQ)
               DO IR = 1,JRCRI(IM)
                  VNEW(IR,1,IT) = VNEW(IR,1,IT) + VSIM
               END DO
            END DO
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C   add the Madelung potential to the Coulomb potential in  VNEW
C   using the Madelung coefficients    AVMAD
C   together with the charge moments   CMNTMTQ  and  CMNTISQ
C   USE_KLMFP=.T.:   restrict calculation to symmetry allowed terms
C-----------------------------------------------------------------------
C
      CALL SCFMAD_POT(USE_KLMFP,NLFP,NLMFP,VNEW,DROT_QLM)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      IF ( SYMMETRIZE_POT ) THEN
CC
C         CALL SYMMETRIZE_FLMT_REP('V  ',VNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C         CALL SYMMETRIZE_FLMT_REP('B  ',BNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C      END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF ( NFCOR ) CALL FPVNFCOR(USE_KLMFP,NLFP,NLMFP,R2RHO,
     &                           RHAT_LEBGRID,W_LEBGRID,N_LEBGRID,
     &                           NMAX_LEBGRID,DSYM_RLM,VNEW)
C
C ----------------------------------------------------------------------
C                             total energy
C ----------------------------------------------------------------------
C
      IF ( KTE.NE.0 ) THEN
C
         CALL FPEPOTINB(EPOTIN,NSPIN,R2RHO,VT,BT,VNST,BNST)
C
         CALL FPECOUB(CMNTMTT,ECOU,NLFP,NSPIN,R2RHO,VNEW,KVMAD,
     &                NRGNT_VSF_LM,LMRGNT_VSF,RGNT_VSF,NRGNT_VSF)
C
      END IF
C
C ----------------------------------------------------------------------
C                             force
C ----------------------------------------------------------------------
C
      IF ( CALC_FORCE .AND. Q1M_OCCURS ) THEN
C
         CALL FORCE(CMNTMTT,VNEW,BNEW,F_LMT,FHF_LMT)
C
         IF ( NQ.GT.1 ) WRITE (6,*) ' '
C
         F_LMQ(2:4,1:NQ) = 0D0
         FHF_LMQ(2:4,1:NQ) = 0D0
C
         DO IQ = IQBOT,IQTOP
C
            IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
                  DO LM = 2,4
                     F_LMQ(LM,IQ) = F_LMQ(LM,IQ) + CONC(IT)*F_LMT(LM,IT)
                     FHF_LMQ(LM,IQ) = FHF_LMQ(LM,IQ) + CONC(IT)
     &                                *FHF_LMT(LM,IT)
                  END DO
               END DO
C
C--------------------------- rotate Q-moments from local to global frame
C                 use inverse of DROT_QLM  with  1/DROT_QLM = DROT_QLM^T
C
               IF ( MAGROT_Q(IQ) ) THEN
C
                  WTMPA(2:4) = F_LMQ(2:4,IQ)
                  WTMPB(2:4) = FHF_LMQ(2:4,IQ)
C
                  DO LM = 2,4
                     SUMA = 0D0
                     SUMB = 0D0
                     DO LMP = 2,4
                        D_LMLMP = DROT_QLM(LM,LMP,IQ)
                        SUMA = SUMA + WTMPA(LMP)*D_LMLMP
                        SUMB = SUMB + WTMPB(LMP)*D_LMLMP
                     END DO
                     F_LMQ(LM,IQ) = SUMA
                     FHF_LMQ(LM,IQ) = SUMB
                  END DO
C
               END IF
C
            ELSE
C
               IQREP = IQREPQ(IQ)
               ISYM = ISYMGENQ(IQ)
               DO LM = 2,4
                  SUMA = 0D0
                  SUMB = 0D0
                  DO LMP = 2,4
                     D_LMLMP = DSYM_RLM(LM,LMP,ISYM)
                     SUMA = SUMA + F_LMQ(LMP,IQREP)*D_LMLMP
                     SUMB = SUMB + FHF_LMQ(LMP,IQREP)*D_LMLMP
                  END DO
                  F_LMQ(LM,IQ) = SUMA
                  FHF_LMQ(LM,IQ) = SUMB
               END DO
            END IF
C
         END DO
C
         WRITE (6,99002)
         DO IQ = IQBOT,IQTOP
            WRITE (6,99003) IQ,FHF_LMQ(4,IQ),FHF_LMQ(2,IQ),FHF_LMQ(3,IQ)
     &                      ,F_LMQ(4,IQ),F_LMQ(2,IQ),F_LMQ(3,IQ)
         END DO
         WRITE (6,*) ' '
C
C ----------------------  store forces for embedded cluster calculations
         IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
            DO IQ = IQBOT,IQTOP
               IQCLU = IQ - NQHOST
               IF ( IQCLU.LT.1 .OR. IQCLU.GT.NQCLU )
     &              CALL STOP_MESSAGE(ROUTINE,'IQCLU out of range')
               FORCE_QCLU(1,IQCLU) = F_LMQ(4,IQ)
               FORCE_QCLU(2,IQCLU) = F_LMQ(2,IQ)
               FORCE_QCLU(3,IQCLU) = F_LMQ(3,IQ)
            END DO
         END IF
C
      END IF
C
C ----------------------------------------------------------------------
C           determine  ONLY  electric field gradient  EFG
C             using electrostatic part of the potential  V
C             the nucelar potential does not contribute
C ----------------------------------------------------------------------
C
      IF ( CALC_EFG ) THEN
         CALL FPEFG(CMNTMTT,R2RHO,VNEW)
         RETURN
      END IF
C
C ----------------------------------------------------------------------
C           calculate kinetic energy denisty if requested
C ----------------------------------------------------------------------
C
      IF ( CALC_KINETIC_ENERGY_DENSITY )
     &     CALL FPNRCORE_TAU(IPRINT,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &     NMAX_LEBGRID)
C
C ----------------------------------------------------------------------
C           add the exchange correlation potential to  V  and  B
C ----------------------------------------------------------------------
C
      CALL FPVXCLM(EXC,KTE,NSPIN,R2RHO,VNEW,BNEW,RGNT_VSF,LMRGNT_VSF,
     &             NRGNT_VSF_LM,NRGNT_VSF,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &             NMAX_LEBGRID)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      IF ( SYMMETRIZE_POT ) THEN
CC
C         CALL SYMMETRIZE_FLMT_REP('V  ',VNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C         CALL SYMMETRIZE_FLMT_REP('B  ',BNEW,FAUX,DSYM_RLM,DROT_QLM,
C     &                            WRLM1,WRLM2)
CC
C      END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      IF ( KTE.EQ.1 ) CALL FPETOT(ECOU,EFERMI,EPOTIN,EXC,ETOT)
C
C-----------------------------------------------------------------------
C            add contribution Coulomb potential of the nucleus to  V
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRCRIT = JRCUT(NPAN(IM),IM)
         ZZOR = 2.0D0*Z(IT)*SQRT_4PI
         DO I = 1,IRCRIT
            VNEW(I,1,IT) = VNEW(I,1,IT) - ZZOR/R(I,IM)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                  symmetrize potential VNEW and BNEW
C-----------------------------------------------------------------------
C
      IF ( SYMMETRIZE_POT ) THEN
C
         CALL SYMMETRIZE_FLMT_REP('V  ',VNEW,FAUX,DSYM_RLM,DROT_QLM,
     &                            WRLM1,WRLM2)
C
         CALL SYMMETRIZE_FLMT_REP('B  ',BNEW,FAUX,DSYM_RLM,DROT_QLM,
     &                            WRLM1,WRLM2)
C
      END IF
C
C-----------------------------------------------------------------------
C     determine muffin tin zero and shift potential  V  accordingly
C-----------------------------------------------------------------------
C
      CALL FPVMUFTIN(NLMFP,VNEW,VMTZ)
C
C-----------------------------------------------------------------------
C             write unconvoluted potential data to file
C              with shape function for next iteration
C-----------------------------------------------------------------------
C
      IF ( WRITE_IS_POT ) THEN
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'inter_stitial_pot.dat')
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IRMTIN = JRCUT(1,IM)
            IRCRIT = JRCRI(IM)
C
            WRITE (IOTMP,99006) IT,IM,IRMTIN,IRCRIT
            DO LM = 1,NLMFP
               DO IR = IRMTIN + 1,IRCRIT
                  WRITE (IOTMP,99006) IT,IM,LM,IR,VNEW(IR,LM,IT),
     &                                BNEW(IR,LM,IT)
               END DO
            END DO
C
         END DO
C
         WRITE (6,99007)
C
         CLOSE (IOTMP)
      END IF
C
C-----------------------------------------------------------------------
C             convolute potential functions  B  and  V
C              with shape function for next iteration
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCRI(IM)
C
         DO LM = 1,NLMFP
            DO IR = 1,IRCRIT - IRMTIN
               VLMTMP(IR,LM) = 0.0D0
               BLMTMP(IR,LM) = 0.0D0
            END DO
         END DO
C
         DO J = 1,NRGNT_VSF_LM(NLMFP)
            LM1 = LMRGNT_VSF(J,1)
            IF ( LM1.GT.NLMFPMAX )
     &            CALL STOP_MESSAGE(ROUTINE,'LMRGNT_VSF > NLMFPMAX')
            LM2 = LMRGNT_VSF(J,2)
            LM3 = LMRGNT_VSF(J,3)
            IF ( KLMSF(LM3,IM).EQ.1 ) THEN
               ISF = ISFLM(LM3,IM)
C ----------------------------------------------------------------------
               IF ( ISF.LE.NSFMAX .AND. ISF.GE.1 ) THEN
C ----------------------------------------------------------------------
                  DO IR = IRMTIN + 1,IRCRIT
                     IRSF = IR - IRMTIN
                     RWGT = RGNT_VSF(J)*FLMSF(IRSF,ISF,IM)
C
                     VLMTMP(IRSF,LM1) = VLMTMP(IRSF,LM1)
     &                                  + RWGT*VNEW(IR,LM2,IT)
                     BLMTMP(IRSF,LM1) = BLMTMP(IRSF,LM1)
     &                                  + RWGT*BNEW(IR,LM2,IT)
                  END DO
C ----------------------------------------------------------------------
               ELSE
                  WRITE (6,*) ' ISF: ',ISF,J,LM1,LM2,LM3,IM,NSFMAX
                  CALL STOP_MESSAGE(ROUTINE,'ISF out of range  !!')
               END IF
C ----------------------------------------------------------------------
            END IF
         END DO
C ----------------------------------------------------------------------
C
         DO LM = 1,NLMFP
            DO IR = IRMTIN + 1,IRCRIT
               IRSF = IR - IRMTIN
               VNEW(IR,LM,IT) = VLMTMP(IRSF,LM)
               BNEW(IR,LM,IT) = BLMTMP(IRSF,LM)
            END DO
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C    convert to the standard SPR-KKR potential functions
C    VT    BT     spherical part INCLUDING Y00
C    VNST  BNST   NON-spherical parts for LM>1
C-----------------------------------------------------------------------
C
      Y00 = 1D0/SQRT_4PI
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         DO IR = 1,JRCRI(IM)
            VT(IR,IT) = VNEW(IR,1,IT)*Y00
            BT(IR,IT) = BNEW(IR,1,IT)*Y00 - BEXT
         END DO
C
         DO LM = 2,NLMFPT(IT)
            IF ( KLMFP(LM,IT).NE.0 .AND. LM.LE.NLMFPMAX ) THEN
               DO IR = JRNS1(IM),JRCRI(IM)
                  VNST(IR,LM,IT) = VNEW(IR,LM,IT)
                  BNST(IR,LM,IT) = BNEW(IR,LM,IT)
               END DO
            END IF
         END DO
C
      END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 5,6
      IF ( BREAKPOINT.EQ.5 ) CALL BREAKPOINT_5(CMNTMTT,CMNTIST,DROT_QLM)
C
      IF ( BREAKPOINT.EQ.6 ) CALL STOP_BREAKPOINT(ROUTINE)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 5,6
C
      DEALLOCATE (CMNTIST,CMNTISQ,CMNTMTT,CMNTMTQ)
      DEALLOCATE (ECOU,EPOTIN,EXC,QNETT,VLMTMP,BLMTMP)
C
      IF ( IFLAG_NEUTRALITY.GT.0 .AND. SYSTEM_TYPE(1:16)
     &     .NE.'EMBEDDED-CLUSTER' .AND. SYSTEM_TYPE(1:10)
     &     .NE.'LIV       ' .AND. SYSTEM_TYPE(1:10).NE.'LIR       ' )
     &     THEN
         WRITE (6,99005)
         CALL STOP_MESSAGE(ROUTINE,'IFLAG_NEUTRALITY > 0')
      END IF
C
CC-----------------------------------------------------------------------
99001 FORMAT (10X,'ERROR: for type IT=',I3,'  shape function ISF=',I3,
     &        'with LM=',I3,/,10X,'couples with LM1=',I3,
     &        ' to forbidden',' LM2=',I3)
99002 FORMAT (/,5X,'site dependent forces  (a.u.)',//18X,
     &        'Hellmann-Feynman',23X,' total'//,6X,
     &        'IQ      F_x       F_y       F_z',
     &        '           F_x       F_y       F_z')
99003 FORMAT (I8,2(1X,3F10.4,:,3X))
99004 FORMAT (2(/,1X,79('*')),/,19X,'setting up new potential',
     &        ' in <FPSCFNEWPOT>',2(/,1X,79('*')),/)
99005 FORMAT (2(/,1X,79('#')),/,10X,
     &        'inconsitencies detected by <FPSCFNEWPOT>',/,10X,
     &        'see output for information - execution stopped',
     &        2(/,1X,79('#')),/)
99006 FORMAT (4I4,2E25.14)
99007 FORMAT (/,10X,'interstitial potential written to ',
     &        'file    inter_stitial_pot.dat',/)
      END
C*==symmetrize_flmt_rep.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE SYMMETRIZE_FLMT_REP(STR3,FLMT,FAUX,DSYM_RLM,DROT_QLM,
     &                               WRLM1,WRLM2)
C   ********************************************************************
C   *                                                                  *
C   *   symmetrize the function  FLMT  expanded in real spherical      *
C   *   harmonics for the actual types ITBOT ... ITTOP                 *
C   *                                                                  *
C   *         f_L = 1/h SUM_{SYM} SUM_{L'}  f_L' D{SYM}_L'L            *
C   *                                                                  *
C   *   It is assumed that only FLMT for the representative type       *
C   *   is available. For that reasons only symmetry operations SYM    *
C   *   that do not transfer the site q_rep to q'' are accounted for   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYM,IQREPQ,SYMACCEPTED,IQPSYMQ
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_SITES,ONLY:IQBOT,IQTOP,NLMQMAD,IQAT
      USE MOD_TYPES,ONLY:NLMFPMAX,NTMAX,NLMFPT,KLMFP,IMT,ITBOT,ITTOP
      IMPLICIT NONE
C*--SYMMETRIZE_FLMT_REP779
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPSCFNEWPOT')
      LOGICAL CHECK_SYMMETRIZE
      PARAMETER (CHECK_SYMMETRIZE=.TRUE.)
C
C Dummy arguments
C
      CHARACTER*3 STR3
      REAL*8 DROT_QLM(NLMQMAD,NLMQMAD,IQBOT:IQTOP),
     &       DSYM_RLM(NLMQMAD,NLMQMAD,NSYM),FAUX(NRMAX,NLMFPMAX),
     &       FLMT(NRMAX,NLMFPMAX,NTMAX),WRLM1(NLMQMAD,NLMQMAD),
     &       WRLM2(NLMQMAD,NLMQMAD)
C
C Local variables
C
      REAL*8 D_LMLMP,RELDIFF,RELDIFFMAX,SORG,SSYM,XNORM,XORG,XSYM
      INTEGER IFLAG_SIGN,IM,IQ,IR,ISYM,IT,IT_SELECT,LM,LMP,LM_SIGN,
     &        NSYMLOC
C
C*** End of declarations rewritten by SPAG
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      IT_SELECT = 17
      IT_SELECT = 0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C-----------------------------------------------------------------------
C                         symmetrize FLMT
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( IT.EQ.IT_SELECT ) THEN
            IM = IMT(IT)
            WRITE (101,*) '######################################',IT,
     &                    NLMFPT(IT)
            DO IR = 100,JRCRI(IM)
               WRITE (101,'(40f15.8)') (FLMT(IR,LM,IT),LM=2,NLMFPT(IT))
            END DO
         END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
         IM = IMT(IT)
         IQ = IQAT(1,IT)
         IF ( IQREPQ(IQ).NE.IQ )
     &        CALL STOP_MESSAGE(ROUTINE,'IQREPQ - Q ')
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( IT.EQ.IT_SELECT ) THEN
            WRITE (*,*) '############ '
            WRITE (*,*) '############ IQ IT',IQ,IT
Ccc         CALL  RMATSTR('DROT_QLM     ',DROT_QLM(1,1,IQ),NLMQMAD,
Ccc     &                NLMQMAD,1,1,0,1D-8,6)
         END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
         FAUX(:,:) = 0D0
         NSYMLOC = 0
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               IF ( IQ.EQ.IQPSYMQ(ISYM,IQ) ) THEN
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF ( IT.EQ.IT_SELECT ) THEN
                     WRITE (*,*) '############ '
                     WRITE (*,*) '############ ISYM',ISYM
Ccc                  CALL  RMATSTR('DSYM_RLM     ',DSYM_RLM(1,1,ISYM),
Ccc     &                         NLMQMAD,NLMQMAD,1,1,0,1D-8,6)
                  END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  NSYMLOC = NSYMLOC + 1
C
C--------- if needed: transfer the symmetry operation to the local frame
C
                  CALL DGEMM('N','N',NLMQMAD,NLMQMAD,NLMQMAD,1D0,
     &                       DSYM_RLM(1,1,ISYM),NLMQMAD,DROT_QLM(1,1,IQ)
     &                       ,NLMQMAD,0D0,WRLM1,NLMQMAD)
C
                  CALL DGEMM('T','N',NLMQMAD,NLMQMAD,NLMQMAD,1D0,
     &                       DROT_QLM(1,1,IQ),NLMQMAD,WRLM1,NLMQMAD,0D0,
     &                       WRLM2,NLMQMAD)
C
                  DO LM = 2,NLMFPT(IT)
                     IF ( KLMFP(LM,IT).NE.0 ) THEN
C
                        DO LMP = 2,NLMFPT(IT)
                           IF ( KLMFP(LMP,IT).NE.0 ) THEN
C
C                             D_LMLMP = DSYM_RLM(LM,LMP,ISYM)
                              D_LMLMP = WRLM2(LM,LMP)
C
                              CALL DAXPY(JRCRI(IM),D_LMLMP,
     &                           FLMT(1,LMP,IT),1,FAUX(1,LM),1)
C
                           END IF
                        END DO
C
                     END IF
                  END DO
C
               END IF
            END IF
         END DO
C
         XNORM = 1D0/DBLE(NSYMLOC)
         RELDIFFMAX = 0D0
         IFLAG_SIGN = 0
         DO LM = 2,NLMFPT(IT)
            IF ( KLMFP(LM,IT).NE.0 ) THEN
C
C-----------------------------------------------------------------------
               IF ( CHECK_SYMMETRIZE ) THEN
C
                  SORG = 0D0
                  SSYM = 0D0
                  RELDIFF = 0D0
                  DO IR = 1,JRCRI(IM)
                     XORG = FLMT(IR,LM,IT)
                     XSYM = FAUX(IR,LM)*XNORM
                     SORG = SORG + XORG
                     SSYM = SSYM + XSYM
                     IF ( ABS(XORG).GT.1D-10 ) THEN
                        RELDIFF = ABS(1D0-XSYM/XORG)
                     ELSE IF ( ABS(XSYM).GT.1D-10 ) THEN
                        RELDIFF = ABS(1D0-XORG/XSYM)
                     END IF
                     RELDIFFMAX = MAX(RELDIFF,RELDIFFMAX)
                  END DO
                  IF ( ABS(SORG).GT.1D-8 .OR. ABS(SSYM).GT.1D-8 ) THEN
                     IF ( ABS(SORG).LT.1D-8 .OR. ABS(SSYM).LT.1D-8 )
     &                    THEN
                        IFLAG_SIGN = 1
                        LM_SIGN = LM
                     ELSE IF ( SORG*SSYM.LT.0D0 ) THEN
                        IFLAG_SIGN = 1
                        LM_SIGN = LM
                     END IF
                     IF ( IFLAG_SIGN.NE.0 ) THEN
                        IFLAG_SIGN = 0
                        WRITE (6,99002) STR3,IT,LM_SIGN,SORG,SSYM
                     END IF
                  END IF
C
               END IF
C-----------------------------------------------------------------------
C
               FLMT(1:JRCRI(IM),LM,IT) = FAUX(1:JRCRI(IM),LM)*XNORM
C
            END IF
         END DO
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( IT.EQ.IT_SELECT ) THEN
            WRITE (*,*) 'AAAAAA'
            WRITE (102,*) '######################################',IT,
     &                    NLMFPT(IT)
            DO IR = 100,JRCRI(IM)
               WRITE (*,*) 'AAAAAA',IR
               WRITE (102,'(40f15.8)') (FLMT(IR,LM,IT),LM=2,NLMFPT(IT))
            END DO
            WRITE (*,*) 'BBBBBB'
         END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
         IF ( CHECK_SYMMETRIZE ) WRITE (6,99001) STR3,IT,RELDIFFMAX
      END DO
C
99001 FORMAT (5X,'<SYMMETRIZE_FLMT_REP>  ',A,':  IT=',I3,3X,
     &        'max rel diff = ',F20.10)
99002 FORMAT (5X,'<SYMMETRIZE_FLMT_REP>  ',A,': IT=',I3,
     &        ' sign problem for LM =',I3,' SUM org - sym: ',2F20.12)
C
      END
