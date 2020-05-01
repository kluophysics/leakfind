C*==fpvnfcor.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE FPVNFCOR(USE_KLMFP,NLFP,NLMFP,R2RHO,RHAT_LEBGRID,
     &                    W_LEBGRID,N_LEBGRID,NMAX_LEBGRID,DROTRLM_V,
     &                    VNEW)
C   ********************************************************************
C   *                                                                  *
C   *  set up the near-field correction to the Coulomb potential       *
C   *                                                                  *
C   *  calculate the Coulomb potential due to the nearest neighbors    *
C   *  and correct the Madelung potential accordingly                  *
C   *                                                                  *
C   * 2011 HE + MKO                                                    *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_LATTICE,ONLY:ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,ADAINV_I,
     &    ADAINV_R,SYSTEM_DIMENSION,ALAT
      USE MOD_SITES,ONLY:QBAS,NQ,NQMAX,NQ_L,NQ_R,IMQ,IQBOT,IQTOP,
     &    NLMQMAD,CMNTQ,KFP_LMQ,NOQ,ITOQ
      USE MOD_SYMMETRY,ONLY:NSYM,IQREPQ,ISYMGENQ,SYMEULANG
      USE MOD_RMESH,ONLY:R,DRDI,JRCUT,NPAN,FLMSF,KLMSF,NRMAX,NSFMAX,
     &    ISFLM,NM,NRSFMAX,NLMSFMAX
      USE MOD_TYPES,ONLY:KLMFP,NTMAX,NLFPMAX,NLMFPMAX,Z,CONC
      USE MOD_CONSTANTS,ONLY:PI,SQRT_4PI
      USE MOD_FILES,ONLY:IPRINT,IOTMP,IDUMMY,RDUMMY
      USE MOD_ANGMOM,ONLY:NL,L_LM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPVNFCOR')
      REAL*8 CONST_4PISQ
      PARAMETER (CONST_4PISQ=4D0*PI*PI)
      INTEGER NQNBRMAX
      PARAMETER (NQNBRMAX=27)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER NLFP,NLMFP,NMAX_LEBGRID,N_LEBGRID
      LOGICAL USE_KLMFP
      REAL*8 DROTRLM_V(NLMQMAD,NLMQMAD,NSYM),
     &       R2RHO(NRMAX,NLMFPMAX,NTMAX,3),RHAT_LEBGRID(3,NMAX_LEBGRID),
     &       VNEW(NRMAX,NLMFPMAX,NTMAX),W_LEBGRID(NMAX_LEBGRID,NLMFPMAX)
C
C Local variables
C
      REAL*8 AC,AFRLN(:),ALN(:),AVMADNBR_Q(:,:,:,:),CLURAD_NFC,DFAC(:,:)
     &       ,DQCLU_NFC(:),DROT_RLM(:,:,:),DVNF,DYMAX,D_LMLMP,FAC,FAC0,
     &       FAC1,FAC2,FAC3,FLMSFP(:,:),FLN(:),INT0(:),INT1R(:,:),
     &       INT2R(:,:),INT3R(:,:),M2SUM(:,:,:,:),M3SUM(:,:,:),MR13PWL1,
     &       MR13PWML1,R1,R13,R13LN,R23(:),R23PW(:,:),R23SQ(:,:),
     &       R2RHOQ(:,:),R2RHOQP(:,:),R2RHOTILQ(:,:,:),R3,RABS3,RABS_J,
     &       RAT,RCRITMAX,RCRIT_I,RCRIT_J,RGNT_FPTILLAM(:),
     &       RGNT_TILFPSF(:),RHAT3(3),RHAT_J(3),RL,RPWL(:),RPWL_M(:,:,:)
     &       ,RQCLU_NFC(:,:),RQNBR_Q(:,:,:),RSQ,RVEC3(3),RVEC_I(3),
     &       RVEC_J(3),RYLM(:),TIME,TIME0,USUM(:),V1(:),V2(:),VIC(:),
     &       VINT1(:),VINT2(:),VINTRA(:),VINTRAQ(:,:,:),VLMJ,
     &       VMADQ(:,:,:),VNFCORQ(:,:,:),VNFCORQ1(:,:,:),VR23(:),
     &       VSPHERE(:),WINT(:),WJA,WJB,WR,WR23PW(:),X_DRDI(:),Y,ZZOR
      REAL*8 DDOT,DNRM2
      CHARACTER*15 FILNAM
      INTEGER I,IA_ERR,IM,IO,IPRINT_LOW,IQ,IQCNTR,IQREP,IQ_QCLU_NFC(:),
     &        IQ_QNBR_Q(:,:),IR,IRCRIT,IRMTIN,IRSF,IR_START,
     &        IR_START_MIN,ISF,ISYM,IT,I_LEBGRID,J,JM,JQ,JQCLU,JQNBR,JR,
     &        JRA,JRB,JRSTART_IR(:),JR_START,JR_START_MIN,KFP_L1(:),
     &        KRHOTIL_LMQ(:,:),L,L1,LM,LM1,LM2,LM3,LMAX_FP,LMAX_RHOTIL,
     &        LMAX_TMP,LMP,LMRGNT_FPTILLAM(:,:),LMRGNT_TILFPSF(:,:),M,
     &        M1,N,NLAM,NLM_RHOTIL,NLM_TMP,NL_RHOTIL,NL_TMP,NQCLU_I_NFC,
     &        NQCLU_L_NFC,NQCLU_NFC,NQCLU_R_NFC,NQNBR_Q(:),
     &        NRGNT_FPTILLAM,NRGNT_FPTILLAM_LM1(:),NRGNT_TILFPSF,
     &        NRGNT_TILFPSF_LM1(:),NSHLCLU_NFC,NSYMH,SLN(:),V
      LOGICAL INITIALIZE,MOL
      SAVE ALN,AVMADNBR_Q,DFAC,DROT_RLM,FLN,IQ_QNBR_Q,LMRGNT_FPTILLAM,
     &     LMRGNT_TILFPSF,NLAM,NQNBR_Q,NRGNT_FPTILLAM,
     &     NRGNT_FPTILLAM_LM1,NRGNT_TILFPSF,NRGNT_TILFPSF_LM1,
     &     RGNT_FPTILLAM,RGNT_TILFPSF,RQNBR_Q,SLN
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE RGNT_FPTILLAM,LMRGNT_FPTILLAM,NRGNT_FPTILLAM_LM1
      ALLOCATABLE RGNT_TILFPSF,LMRGNT_TILFPSF,NRGNT_TILFPSF_LM1
      ALLOCATABLE DFAC,ALN,FLN,SLN,R23SQ,M2SUM,M3SUM,RPWL_M,X_DRDI
      ALLOCATABLE INT0,INT1R,INT2R,INT3R,R23,JRSTART_IR
      ALLOCATABLE VNFCORQ,VMADQ,VINTRA,R2RHOTILQ,R2RHOQ,R23PW
      ALLOCATABLE VSPHERE,VIC,RYLM,VNFCORQ1,R2RHOQP,FLMSFP
      ALLOCATABLE RQCLU_NFC,DQCLU_NFC,IQ_QCLU_NFC,AVMADNBR_Q
      ALLOCATABLE RPWL,V1,V2,VINT1,VINT2
      ALLOCATABLE NQNBR_Q,RQNBR_Q,IQ_QNBR_Q,VINTRAQ
      ALLOCATABLE DROT_RLM,KRHOTIL_LMQ,KFP_L1
      ALLOCATABLE WINT,AFRLN,USUM,VR23,WR23PW
C
      ALLOCATE (V1(NRMAX),V2(NRMAX),RPWL(NRMAX),X_DRDI(NRMAX))
      ALLOCATE (VINT1(NRMAX),VINT2(NRMAX))
C
      WRITE (6,99004)
C
C=======================================================================
C     angular momentum parameters for convoluted charge density
C=======================================================================
C
      IF ( NLFP.NE.2*NL-1 ) STOP '<FPVNFCOR>  NLFP .NE.  2*NL - 1'
      LMAX_FP = NLFP - 1
C
C      LMAX_RHOTIL = 2*(NL-1) + 4*(NL-1)
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      LMAX_RHOTIL = LMAX_FP
Cc      LMAX_RHOTIL = 8
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      NL_RHOTIL = LMAX_RHOTIL + 1
      NLM_RHOTIL = NL_RHOTIL*NL_RHOTIL
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         MOL = .FALSE.
         IF ( IPRINT.LE.0 ) THEN
            IPRINT_LOW = -1
         ELSE
            IPRINT_LOW = IPRINT
         END IF
C
         ALLOCATE (NQNBR_Q(NQMAX),RQNBR_Q(3,NQNBRMAX,NQMAX))
         ALLOCATE (IQ_QNBR_Q(NQNBRMAX,NQMAX))
         NQNBR_Q(1:NQMAX) = 0
C
         RCRITMAX = 0D0
         DO IM = 1,NM
            RCRITMAX = MAX(RCRITMAX,R(JRCUT(NPAN(IM),IM),IM))
         END DO
         WRITE (6,99001) RCRITMAX
C
         DO IQ = IQBOT,IQTOP
C
            IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
               IM = IMQ(IQ)
               RCRIT_I = R(JRCUT(NPAN(IM),IM),IM)
C
               IQCNTR = IQ
               NSHLCLU_NFC = 0
               CLURAD_NFC = (2*RCRITMAX+1D-6)/ALAT
C
               CALL CLUSSITES(IOTMP,IPRINT_LOW,MOL,SYSTEM_DIMENSION,
     &                        ABAS,ABAS_L,ABAS_I,ABAS_R,ADAINV_L,
     &                        ADAINV_I,ADAINV_R,QBAS,CLURAD_NFC,IQCNTR,
     &                        NQCLU_NFC,NQCLU_L_NFC,NQCLU_I_NFC,
     &                        NQCLU_R_NFC,NSHLCLU_NFC,NQ,NQ_L,NQ_R,
     &                        NQMAX)
C
               IF ( ALLOCATED(RQCLU_NFC) )
     &              DEALLOCATE (RQCLU_NFC,IQ_QCLU_NFC,DQCLU_NFC)
               ALLOCATE (RQCLU_NFC(3,NQCLU_NFC),DQCLU_NFC(NQCLU_NFC))
               ALLOCATE (IQ_QCLU_NFC(NQCLU_NFC),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 ) STOP 'allocation in <FPVNFCOR>'
C
               READ (IOTMP) ((RQCLU_NFC(J,I),J=1,3),DQCLU_NFC(I),
     &                      IQ_QCLU_NFC(I),I=1,NQCLU_NFC),
     &                      (IDUMMY,I=1,NSHLCLU_NFC)
               CLOSE (IOTMP)
C
               WRITE (6,99005)
C
               NQNBR_Q(IQ) = 0
               DO JQCLU = 1,NQCLU_NFC
                  JQ = IQ_QCLU_NFC(JQCLU)
                  JM = IMQ(JQ)
                  RCRIT_J = R(JRCUT(NPAN(JM),JM),JM)
                  IF ( (DQCLU_NFC(JQCLU).LE.(RCRIT_I+RCRIT_J+1D-6)/ALAT)
     &                 .AND. (DQCLU_NFC(JQCLU).GT.0.2D0) ) THEN
                     NQNBR_Q(IQ) = NQNBR_Q(IQ) + 1
                     IF ( NQNBR_Q(IQ).GT.NQNBRMAX ) STOP 'nqnbrmax'
                     JQNBR = NQNBR_Q(IQ)
                     RQNBR_Q(1:3,JQNBR,IQ) = RQCLU_NFC(1:3,JQCLU)
                     IQ_QNBR_Q(JQNBR,IQ) = JQ
                     WRITE (6,99003) IQ,IM,RCRIT_I,NQNBR_Q(IQ),JQ,JM,
     &                               RCRIT_J,DQCLU_NFC(JQCLU)*ALAT,
     &                               (RCRIT_I+RCRIT_J)
                  END IF
               END DO
               WRITE (6,99006) IQ,NQNBR_Q(IQ)
C
            END IF
C
         END DO
C
C-----------------------------------------------------------------------
C
         NLAM = LMAX_FP + LMAX_RHOTIL + 1
C
         ALLOCATE (ALN(-NLAM:NLAM),FLN(0:NLAM),SLN(-NLAM:NLAM))
         ALLOCATE (DFAC(0:NL_RHOTIL,0:NL_RHOTIL))
C
         CALL SACK_FACUL(DFAC,NL_RHOTIL)
C
         CALL SACK_FLN(SLN,ALN,FLN,NLAM)
C
C-----------------------------------------------------------------------
C   Gaunt coefficients for convolution of charge denisty n_L
C   with shape functions f_L'':    nbar_L = n_L' * f_L''
C-----------------------------------------------------------------------
C
C               RGNT_TILFPSF, NRGNT_TILFPSF_LM1, LMRGNT_TILFPSF
C
         ALLOCATE (NRGNT_TILFPSF_LM1(0:NLMFPMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: FPSCFNEWPOT -> RGNT_TILFPSF'
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
         CALL RGNTSF_SETUP(LMAX_FP,LMAX_RHOTIL,NRGNT_TILFPSF_LM1,
     &                     NRGNT_TILFPSF)
C
         ALLOCATE (RGNT_TILFPSF(NRGNT_TILFPSF))
         ALLOCATE (LMRGNT_TILFPSF(NRGNT_TILFPSF,3))
C
         REWIND IOTMP
         DO I = 1,NRGNT_TILFPSF
            READ (IOTMP) (LMRGNT_TILFPSF(I,J),J=1,3),RGNT_TILFPSF(I)
         END DO
         CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C   Gaunt coefficients for expansion of potential V_L in terms of
C   convoluted charge density:        V_L = nbar_L' * Y_L''
C-----------------------------------------------------------------------
C
C               RGNT_FPTILLAM, NRGNT_FPTILLAM_LM1, LMRGNT_FPTILLAM
C
         ALLOCATE (NRGNT_FPTILLAM_LM1(0:NLMFPMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: FPSCFNEWPOT -> RGNT_FPTILLAM'
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
         CALL RGNTSF_SETUP(LMAX_FP,LMAX_RHOTIL,NRGNT_FPTILLAM_LM1,
     &                     NRGNT_FPTILLAM)
C
         ALLOCATE (RGNT_FPTILLAM(NRGNT_FPTILLAM))
         ALLOCATE (LMRGNT_FPTILLAM(NRGNT_FPTILLAM,3))
C
         REWIND IOTMP
         DO I = 1,NRGNT_FPTILLAM
            READ (IOTMP) (LMRGNT_FPTILLAM(I,J),J=1,3),RGNT_FPTILLAM(I)
         END DO
         CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C           get rotation matrices for REAL spherical harmonics
C     allow to deal with full potential arrays  AND  shape functions
C     with              NLFP = 2 * NL - 1       and  NLSF = 4 * NL - 3
C ----------------------------------------------------------------------
C
         LMAX_TMP = 4*(NL-1)
         NL_TMP = LMAX_TMP + 1
         NLM_TMP = NL_TMP**2
C
         ALLOCATE (DROT_RLM(NLM_TMP,NLM_TMP,NSYM))
C
         NSYMH = NSYM/2
C
C-----------------------------------------------------------------------
         DO ISYM = 1,NSYM
C
            IF ( ISYM.LE.NSYMH ) THEN
C
               CALL ROTMAT_RYLM(NL_TMP,NLM_TMP,SYMEULANG(1,ISYM),
     &                          SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                          DROT_RLM(1,1,ISYM))
C
            ELSE
C
               LM2 = 0
               DO L = 0,LMAX_TMP
                  LM1 = LM2 + 1
                  LM2 = LM2 + (2*L+1)
                  IF ( MOD(L,2).EQ.1 ) DROT_RLM(LM1:LM2,LM1:LM2,ISYM)
     &                 = -DROT_RLM(LM1:LM2,LM1:LM2,ISYM-NSYMH)
               END DO
C
            END IF
C
         END DO
C-----------------------------------------------------------------------
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      ALLOCATE (R2RHOTILQ(NRMAX,NLM_RHOTIL,NQMAX))
      ALLOCATE (KRHOTIL_LMQ(NLM_RHOTIL,NQMAX))
      ALLOCATE (R2RHOQ(NRMAX,NLMFP),R2RHOQP(NRMAX,NLMFP))
      ALLOCATE (FLMSFP(NRSFMAX,NLMSFMAX))
      R2RHOTILQ(:,:,:) = 0D0
      KRHOTIL_LMQ(:,:) = 0
C
C=======================================================================
C           convolute type averaged charge density for each site IQ
C              with shape function and multiply with DRDI
C=======================================================================
C
      DO IQ = IQBOT,IQTOP
C
         IM = IMQ(IQ)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCUT(NPAN(IM),IM)
C
C=======================================================================
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            R2RHOQ(:,:) = 0D0
C-------------------------------------------------------------------- IO
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               X_DRDI(1:IRCRIT) = CONC(IT)*DRDI(1:IRCRIT,IM)
C
C-------------------------------------------------------------------  LM
               DO LM = 1,NLMFP
                  IF ( KLMFP(LM,IT).NE.KFP_LMQ(LM,IQ) )
     &                  STOP 'KLMFP(LM,IT) .NE. KFP_LMQ(LM,IQ)'
C----------------------------------------------------------------- KLMFP
                  IF ( .NOT.USE_KLMFP .OR. KLMFP(LM,IT).NE.0 ) THEN
C
                     KRHOTIL_LMQ(LM,IQ) = 1
                     DO IR = 1,IRCRIT
                        R2RHOQ(IR,LM) = R2RHOQ(IR,LM) + X_DRDI(IR)
     &                                  *R2RHO(IR,LM,IT,1)
                     END DO
C
                  END IF
C----------------------------------------------------------------- KLMFP
               END DO
C-------------------------------------------------------------------  LM
            END DO
C-------------------------------------------------------------------- IO
C
C
C-------------------------------------------------------------------  LM
            DO LM = 1,NLMFP
C----------------------------------------------------------------- KLMFP
               IF ( .NOT.USE_KLMFP .OR. KFP_LMQ(LM,IQ).NE.0 ) THEN
C
                  DO IR = 1,IRMTIN
                     R2RHOTILQ(IR,LM,IQ) = R2RHOQ(IR,LM)
                  END DO
C
                  DO J = NRGNT_TILFPSF_LM1(LM-1) + 1,
     &               NRGNT_TILFPSF_LM1(LM)
                     LM2 = LMRGNT_TILFPSF(J,2)
                     LM3 = LMRGNT_TILFPSF(J,3)
                     ISF = ISFLM(LM3,IM)
                     IF ( KLMSF(LM3,IM).EQ.1 .AND. ISF.LE.NSFMAX ) THEN
                        KRHOTIL_LMQ(LM2,IQ) = 1
                        DO IR = IRMTIN + 1,IRCRIT
                           R2RHOTILQ(IR,LM2,IQ) = R2RHOTILQ(IR,LM2,IQ)
     &                        + RGNT_TILFPSF(J)*R2RHOQ(IR,LM)
     &                        *FLMSF(IR-IRMTIN,ISF,IM)
                        END DO
                     END IF
                  END DO
C
               END IF
C----------------------------------------------------------------- KLMFP
            END DO
C-------------------------------------------------------------------  LM
C
C=======================================================================
C
C
C=======================================================================
C               deal with all equivalent sites
C=======================================================================
C
            DO JQ = IQBOT,IQTOP
               IQREP = IQREPQ(JQ)
               IF ( IQREP.EQ.IQ .AND. IQ.NE.JQ ) THEN
C
                  ISYM = ISYMGENQ(JQ)
C
                  DO LM = 1,NLM_RHOTIL
                     DO LMP = 1,NLM_RHOTIL
                        D_LMLMP = DROT_RLM(LM,LMP,ISYM)
                        IF ( ABS(D_LMLMP).GT.1D-8 ) THEN
                           KRHOTIL_LMQ(LM,JQ) = 1
                           DO IR = 1,IRCRIT
                              R2RHOTILQ(IR,LM,JQ) = R2RHOTILQ(IR,LM,JQ)
     &                           + R2RHOTILQ(IR,LMP,IQREP)*D_LMLMP
                           END DO
                        END IF
                     END DO
                  END DO
C
C
C
                  R2RHOQP(:,:) = 0D0
                  DO LM = 1,NLMFP
                     DO LMP = 1,NLMFP
                        D_LMLMP = DROT_RLM(LM,LMP,ISYM)
                        IF ( ABS(D_LMLMP).GT.1D-8 ) THEN
                           DO IR = 1,IRCRIT
                              R2RHOQP(IR,LM) = R2RHOQP(IR,LM)
     &                           + R2RHOQ(IR,LMP)*D_LMLMP
                           END DO
                        END IF
                     END DO
                  END DO
C
                  FLMSFP(:,:) = 0D0
                  DO LM = 1,NLMSFMAX
                     DO LMP = 1,NLMSFMAX
                        D_LMLMP = DROT_RLM(LM,LMP,ISYM)
                        IF ( ABS(D_LMLMP).GT.1D-8 .AND. KLMSF(LMP,IM)
     &                       .NE.0 ) THEN
                           ISF = ISFLM(LMP,IM)
                           DO IRSF = 1,NRSFMAX
                              FLMSFP(IRSF,LM) = FLMSFP(IRSF,LM)
     &                           + FLMSF(IRSF,ISF,IM)*D_LMLMP
                           END DO
                        END IF
                     END DO
                  END DO
C
C-------------------------------------------------------------------  LM
                  DO LM = 1, - NLMFP
C----------------------------------------------------------------- KLMFP
Cc               IF ( .NOT.USE_KLMFP .OR. KFP_LMQ(LM,JQ).NE.0 ) THEN
C
                     KRHOTIL_LMQ(LM,JQ) = 1
                     DO IR = 1,IRMTIN
                        R2RHOTILQ(IR,LM,JQ) = R2RHOQP(IR,LM)
                     END DO
C
                     DO J = NRGNT_TILFPSF_LM1(LM-1) + 1,
     &                  NRGNT_TILFPSF_LM1(LM)
                        LM2 = LMRGNT_TILFPSF(J,2)
                        LM3 = LMRGNT_TILFPSF(J,3)
Cc                     IF ( KLMSF(LM3,IM).EQ.1  ) THEN
                        KRHOTIL_LMQ(LM2,JQ) = 1
                        DO IR = IRMTIN + 1,IRCRIT
                           R2RHOTILQ(IR,LM2,JQ) = R2RHOTILQ(IR,LM2,JQ)
     &                        + RGNT_TILFPSF(J)*R2RHOQP(IR,LM)
     &                        *FLMSFP(IR-IRMTIN,LM3)
                        END DO
Cc                     END IF
                     END DO
C
Cc               END IF
C----------------------------------------------------------------- KLMFP
                  END DO
C-------------------------------------------------------------------  LM
C
               END IF
            END DO
C
         END IF
C
      END DO
C=======================================================================
C
C=======================================================================
C         near-field corrections based on  SACK's method
C=======================================================================
C    correction of Coulomb-potential of Wigner-Seitz-cell  IQ
C    due to the overlap with the Wigner-Seitz-sphere of cell  JQ
C    cf.  Eq.(38) of  R.A. Sack  J. Math. Phys. 5, 260 (1964)
C=======================================================================
C
      CALL CPU_TIME(TIME0)
C
C  =====================================================================
C                 find lower limits for IRSTART and JRSTART
C  =====================================================================
C
      ALLOCATE (JRSTART_IR(NRMAX))
C
      IR_START_MIN = NRMAX
      DO IQ = IQBOT,IQTOP
C
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            IM = IMQ(IQ)
            IRCRIT = JRCUT(NPAN(IM),IM)
C
C  ----------------------------------------------------------- neighbors
            DO JQNBR = 1,NQNBR_Q(IQ)
C
               JQ = IQ_QNBR_Q(JQNBR,IQ)
               JM = IMQ(JQ)
               IRCRIT = JRCUT(NPAN(JM),JM)
C
C -------------------------------------- r_3 = R_i - R_j = - (R_j - R_i)
               R3 = ALAT*DNRM2(3,RQNBR_Q(1,JQNBR,IQ),1)
C
C  evaluate start mesh point - overlap of critical radius of cell J
C                              with point r_i in cell I possible
               DO IR = 2,IRCRIT
                  IF ( R(IR,IM)+R(IRCRIT,JM).GT.R3 ) THEN
                     IR_START = IR
                     EXIT
                  END IF
               END DO
C
               IR_START_MIN = MIN(IR_START_MIN,IR_START)
C
            END DO
         END IF
      END DO
C
      WRITE (6,*) '############# IRSTART_MIN',IR_START_MIN
C
      ALLOCATE (VNFCORQ(NRMAX,NLMFPMAX,NQ),R23SQ(0:NLAM,NRMAX))
      ALLOCATE (INT0(NLMFPMAX))
      ALLOCATE (INT1R(IR_START_MIN:NRMAX,NLMFPMAX))
      ALLOCATE (INT2R(IR_START_MIN:NRMAX,NLMFPMAX))
      ALLOCATE (INT3R(IR_START_MIN:NRMAX,NLMFPMAX))
      ALLOCATE (M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,
     &          NLMFPMAX))
      ALLOCATE (M3SUM(0:(NLAM-1),NLM_RHOTIL,NLMFPMAX))
      ALLOCATE (R23(NRMAX))
      ALLOCATE (RPWL_M(IR_START_MIN:NRMAX,-NLAM:+NLAM,NM))
      ALLOCATE (R23PW(IR_START_MIN:NRMAX,-(LMAX_RHOTIL+1):+(LMAX_RHOTIL+
     &          1)))
      ALLOCATE (AFRLN(0:NLAM),USUM(0:NLAM),VR23(NRMAX))
      ALLOCATE (WR23PW(NRMAX),WINT(NRMAX),KFP_L1(0:LMAX_FP))
C
      VNFCORQ(:,:,:) = 0D0
C
      DO IM = 1,NM
         IRCRIT = JRCUT(NPAN(IM),IM)
         DO L = 0,LMAX_RHOTIL
            IF ( L.EQ.0 ) THEN
               RPWL_M(IR_START_MIN:IRCRIT,0,IM) = 1D0
            ELSE
               DO IR = IR_START_MIN,IRCRIT
                  RPWL_M(IR,+L,IM) = RPWL_M(IR,+L-1,IM)*R(IR,IM)
                  RPWL_M(IR,-L,IM) = RPWL_M(IR,-L+1,IM)/R(IR,IM)
               END DO
            END IF
         END DO
      END DO
C
      JR_START_MIN = NRMAX
C
C  =====================================================================
      DO IQ = IQBOT,IQTOP
C
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            KFP_L1(:) = 0
            LM1 = 0
            DO L1 = 0,LMAX_FP
               DO M1 = -L1,L1
                  LM1 = LM1 + 1
                  IF ( KFP_LMQ(LM1,IQ).NE.0 ) KFP_L1(L1) = 1
               END DO
            END DO
C
            IM = IMQ(IQ)
            IRCRIT = JRCUT(NPAN(IM),IM)
C
C  ----------------------------------------------------------- neighbors
            DO JQNBR = 1,NQNBR_Q(IQ)
C
               JQ = IQ_QNBR_Q(JQNBR,IQ)
               JM = IMQ(JQ)
               IRCRIT = JRCUT(NPAN(JM),JM)
C
C -------------------------------------- r_3 = R_i - R_j = - (R_j - R_i)
               RVEC3(1:3) = -RQNBR_Q(1:3,JQNBR,IQ)
               RABS3 = DNRM2(3,RVEC3,1)
               RHAT3(1:3) = RVEC3(1:3)/RABS3
               R3 = ALAT*RABS3
C
C  evaluate start mesh point - overlap of critical radius of cell J
C                              with point r_i in cell I possible
               DO IR = 2,IRCRIT
                  IF ( R(IR,IM)+R(IRCRIT,JM).GT.R3 ) THEN
                     IR_START = IR
                     EXIT
                  END IF
               END DO
C
C find first mesh point of cell J that overlaps with point r_i in cell I
C
               V = IRCRIT - 1
               DO IR = IR_START,IRCRIT
                  R1 = R(IR,IM)
                  JRSTART_IR(IR) = 2
                  DO JR = V,1, - 1
                     IF ( R1+R(JR,JM).LT.R3 ) THEN
                        JRSTART_IR(IR) = JR + 1
                        JR_START_MIN = MIN(JR_START_MIN,JRSTART_IR(IR))
                        V = JR
                        EXIT
                     END IF
                  END DO
               END DO
C
               IF ( IR_START_MIN.GT.JR_START_MIN )
     &               STOP '<FPVNFCOR>: IR_START_MIN > JR_START_MIN'
C
               JR_START = JRSTART_IR(IRCRIT)
C
C  calculate powers (r'/r'')^i up to i=NLAM
C
               DO JR = JR_START,IRCRIT
                  R23(JR) = R(JR,JM)/R3
                  RSQ = R23(JR)*R23(JR)
                  R23SQ(0,JR) = 1.0D0
                  DO V = 1,NLAM
                     R23SQ(V,JR) = R23SQ(V-1,JR)*RSQ
                  END DO
               END DO
C
C-----------------------------------------------------------------------
C         tabulate   R23^N  N = -(LMAX_FP+1), ... ,+(LMAX_FP+1)
C-----------------------------------------------------------------------
C
               R23PW(JR_START:IRCRIT,0) = 1D0
               DO N = 1,LMAX_RHOTIL + 1
                  R23PW(JR_START:IRCRIT,+N) = R23PW(JR_START:IRCRIT,N-1)
     &               *R23(JR_START:IRCRIT)
                  R23PW(JR_START:IRCRIT,-N)
     &               = 1D0/R23PW(JR_START:IRCRIT,+N)
               END DO
C
C-----------------------------------------------------------------------
C
               CALL SACK_CYSUM(M3SUM,RHAT3,LMRGNT_FPTILLAM,
     &                         NRGNT_FPTILLAM_LM1,RGNT_FPTILLAM,
     &                         NRGNT_FPTILLAM,NLMFP,NLM_RHOTIL,NLAM)
C
               CALL SACK_MMSUM(M2SUM,M3SUM,R2RHOTILQ(1,1,JQ),
     &                         KRHOTIL_LMQ(1,JQ),IR_START_MIN,LMAX_FP,
     &                         LMAX_RHOTIL,NLM_RHOTIL,NLAM)
C
               CALL SACK1(IQ,INT1R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                    JR_START,IRCRIT,JM,LMAX_FP,LMAX_RHOTIL,
     &                    NL_RHOTIL,NLAM)
C
               CALL SACK2(IQ,INT2R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                    JR_START,IRCRIT,JM,LMAX_FP,LMAX_RHOTIL,
     &                    NL_RHOTIL,NLAM)
C
               CALL SACK3(IQ,INT3R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                    JR_START,IRCRIT,JM,LMAX_FP,LMAX_RHOTIL,
     &                    NL_RHOTIL,NLAM)
C
C----------------------------------------------------------- loop over r
               DO IR = IR_START,IRCRIT
C
                  R13 = R(IR,IM)/R3
                  R13LN = 2.0D0*LOG(R13)
C
                  JR_START = JRSTART_IR(IR)
C
                  CALL SACK0(IQ,INT0,SLN,ALN,FLN,R13LN,R23SQ,M2SUM,
     &                       R23PW,KFP_L1,WINT,AFRLN,USUM,VR23,WR23PW,
     &                       IR_START_MIN,JR_START,IRCRIT,JM,LMAX_FP,
     &                       LMAX_RHOTIL,NLAM)
C
                  LM1 = 0
                  MR13PWML1 = 1D0
                  MR13PWL1 = -1D0/R13
                  DO L1 = 0,LMAX_FP
C
                     MR13PWML1 = -MR13PWML1/R13
                     MR13PWL1 = -MR13PWL1*R13
C
                     FAC0 = -CONST_4PISQ/(2.0D0*R3)*MR13PWML1
                     FAC1 = -CONST_4PISQ/R3*MR13PWML1
                     FAC2 = CONST_4PISQ/R3*ABS(MR13PWL1)
                     FAC3 = CONST_4PISQ/R3*MR13PWL1
C
C                    FAC1 = -CONST_4PISQ / R3 * (-R13)**ML1  ML1 = -(L1+1)
C                    FAC2 =  CONST_4PISQ / R3 *   R13 **L1
C                    FAC3 =  CONST_4PISQ / R3 * (-R13)**L1
C
                     DO M1 = -L1,L1
                        LM1 = LM1 + 1
                        IF ( KFP_LMQ(LM1,IQ).NE.0 ) THEN
C
C----------------------------------------------------multiply by 2 = e^2
                           DVNF = FAC1*INT1R(JR_START,LM1)
     &                            + FAC2*INT2R(JR_START,LM1)
     &                            - FAC3*INT3R(JR_START,LM1)
     &                            - 2D0*FAC0*INT0(LM1)
C
                           VNFCORQ(IR,LM1,IQ) = VNFCORQ(IR,LM1,IQ)
     &                        + DVNF
C
                        END IF
                     END DO
C
                  END DO
C
               END DO
C
            END DO
C  ----------------------------------------------------------- neighbors
C
C-----------------------------------------------------------------------
C         apply near field corrections for each type on site IQ
C-----------------------------------------------------------------------
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               DO LM = 1,NLMFP
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
C
                     DO IR = IR_START,IRCRIT
                        VNEW(IR,LM,IT) = VNEW(IR,LM,IT)
     &                     + VNFCORQ(IR,LM,IQ)
                     END DO
C
                  END IF
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END IF
      END DO
C  =====================================================================
C
      IF ( .NOT.CHECK ) RETURN
C
      CALL CPU_TIME(TIME)
      WRITE (6,*) 'Time used for Nearfield Correction',TIME - TIME0
C
C-----------------------------------------------------------------------
C                write results for comparison
C-----------------------------------------------------------------------
      IF ( CHECK .OR. IQBOT.GT.0 ) THEN
         DO IQ = IQBOT,IQTOP
            IF ( IQ.EQ.IQREPQ(IQ) ) THEN
C
               IM = IMQ(IQ)
               IRCRIT = JRCUT(NPAN(IM),IM)
C
               DO LM = 1,NLMFP
                  IF ( KFP_LMQ(LM,IQ).NE.0 ) THEN
C                    WRITE (*,*) '############### L LM',L_LM(LM),LM
C
                     WRITE (FILNAM(1:10),'("sack",I3.3,I3.3)') IQ,LM
                     CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:10))
C
                     DO IR = 100,IRCRIT
C
                        WRITE (IOTMP,99002) R(IR,IM),VNFCORQ(IR,LM,IQ)
                     END DO
                     CLOSE (IOTMP)
                  END IF
               END DO
            END IF
         END DO
      END IF
C
C#######################################################################
      DYMAX = 0D0
      DO IQ = IQBOT,IQTOP
         IF ( IQ.EQ.IQREPQ(IQ) ) THEN
C
            IM = IMQ(IQ)
            IRCRIT = JRCUT(NPAN(IM),IM)
C
            DO LM = 1,NLMFP
               IF ( KFP_LMQ(LM,IQ).NE.0 ) THEN
C                    WRITE (*,*) '############### L LM',L_LM(LM),LM
C
                  WRITE (FILNAM(1:10),'("sack",I3.3,I3.3)') IQ,LM
                  CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'A1/'//FILNAM(1:10)
     &                                 )
C
                  DO IR = 100,IRCRIT
                     READ (IOTMP,99002) RDUMMY,Y
                     DYMAX = MAX(DYMAX,ABS(Y-VNFCORQ(IR,LM,IQ)))
                  END DO
                  CLOSE (IOTMP)
               END IF
            END DO
         END IF
      END DO
      WRITE (*,*) '############ DELTA_max :',DYMAX
      IF ( DYMAX.GT.6D-13 ) THEN
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
         WRITE (*,*) '############ DELTA_max :',DYMAX
      END IF
C#######################################################################
C
C
C
C
C=======================================================================
C
      IF ( .NOT.CHECK ) RETURN
C
C=======================================================================
C              Zabloudil's method -- only for checking
C=======================================================================
C
      CALL CPU_TIME(TIME0)
C
      ALLOCATE (VSPHERE(NMAX_LEBGRID),RYLM(NLMFPMAX))
      ALLOCATE (VNFCORQ1(NRMAX,NLMFPMAX,NQ))
      ALLOCATE (VMADQ(NRMAX,NLMFPMAX,NQ),VIC(NRMAX))
      ALLOCATE (VINTRAQ(NRMAX,NLMFPMAX,NQ),VINTRA(NRMAX))
C
      VINTRAQ(:,:,:) = 0D0
      VNFCORQ1(:,:,:) = 0D0
      VMADQ(:,:,:) = 0D0
C
C-----------------------------------------------------------------------
      ALLOCATE (AVMADNBR_Q(NQNBRMAX,NQMAX,NLMQMAD,NLMQMAD))
      CALL SCFMADNBR(AVMADNBR_Q,NQNBR_Q,RQNBR_Q,NQNBRMAX)
C-----------------------------------------------------------------------
C
      DO IQ = IQBOT,IQTOP
C
         IM = IMQ(IQ)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCUT(NPAN(IM),IM)
C
C=======================================================================
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
C-------------------------------------------------------------------- IO
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
C--------------------------------------------------------------------  L
               DO L = 0,LMAX_FP
C
                  IF ( L.EQ.0 ) THEN
                     RPWL(:) = 1D0
                  ELSE
                     DO IR = 1,IRCRIT
                        RPWL(IR) = RPWL(IR)*R(IR,IM)
                     END DO
                  END IF
C
                  FAC = 8.0D0*PI/DBLE(2*L+1)
C
C--------------------------------------------------------------------  M
                  DO M = -L,L
                     LM = L*L + L + M + 1
C
C----------------------------------------------------------------- KLMFP
                     IF ( .NOT.USE_KLMFP .OR. KLMFP(LM,IT).NE.0 ) THEN
C
C-----------------------------------------------------------------------
C     set up of the integrands V1 and V2
C-----------------------------------------------------------------------
C
                        DO IR = 1,IRCRIT
                           RL = RPWL(IR)
                           V1(IR) = R2RHOTILQ(IR,LM,IQ)*RL
                           V2(IR) = R2RHOTILQ(IR,LM,IQ)/RL/R(IR,IM)
                        END DO
C
C-----------------------------------------------------------------------
C     now integrate V1 and V2
C-----------------------------------------------------------------------
C
                        CALL RRADINT_R(IM,V1,VINT1)
                        CALL RRADINT_INW_R(IM,V2,VINT2)
C
C-----------------------------------------------------------------------
C     collect all parts
C-----------------------------------------------------------------------
C
                        DO IR = 1,IRCRIT
                           RL = RPWL(IR)
                           VINTRA(IR) = FAC*(VINT1(IR)/R(IR,IM)/RL+VINT2
     &                                  (IR)*RL)
                        END DO
C
                        IF ( LM.EQ.1 ) THEN
                           ZZOR = 2.0D0*Z(IT)*SQRT_4PI
                           VINTRA(1:IRCRIT) = VINTRA(1:IRCRIT)
     &                        - ZZOR/R(1:IRCRIT,IM)
                        END IF
C
                        VINTRAQ(1:IRCRIT,LM,IQ)
     &                     = VINTRAQ(1:IRCRIT,LM,IQ) + CONC(IT)
     &                     *VINTRA(1:IRCRIT)
C
                     END IF
C----------------------------------------------------------------- KLMFP
C
                  END DO
C--------------------------------------------------------------------  M
C
               END DO
C--------------------------------------------------------------------  L
C
            END DO
C-------------------------------------------------------------------- IO
C
C=======================================================================
         ELSE
C
            IQREP = IQREPQ(IQ)
            ISYM = ISYMGENQ(IQ)
            DO LM = 1,NLMQMAD
               DO LMP = 1,NLMQMAD
                  D_LMLMP = DROTRLM_V(LM,LMP,ISYM)
                  DO IR = 1,IRCRIT
                     VINTRAQ(IR,LM,IQ) = VINTRAQ(IR,LM,IQ)
     &                  + VINTRAQ(IR,LMP,IQREP)*D_LMLMP
                  END DO
               END DO
            END DO
C
         END IF
C
      END DO
C=======================================================================
C
C
C=======================================================================
C                           METHOD 1
C=======================================================================
C
      DO IQ = IQBOT,IQTOP
C
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            IM = IMQ(IQ)
            IRCRIT = JRCUT(NPAN(IM),IM)
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO L = 0,LMAX_FP
C
               IF ( L.EQ.0 ) THEN
                  RPWL(:) = 1D0
               ELSE
                  DO IR = 1,IRCRIT
                     RPWL(IR) = -RPWL(IR)*R(IR,IM)
                  END DO
               END IF
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
               DO M = -L,L
                  LM = L*L + L + M + 1
C
C=======================================================================
C
                  IF ( .NOT.USE_KLMFP .OR. KFP_LMQ(LM,IQ).NE.0 ) THEN
C
                     AC = 0.0D0
C
                     DO JQNBR = 1,NQNBR_Q(IQ)
                        JQ = IQ_QNBR_Q(JQNBR,IQ)
                        DO LM2 = 1,NLMFP
                           IF ( KFP_LMQ(LM,JQ).NE.0 ) AC = AC + 
     &                          AVMADNBR_Q(JQNBR,IQ,LM,LM2)
     &                          *CMNTQ(LM2,JQ)
                        END DO
                     END DO
C
C-----------------------------------------------------------------------
C  add the type independent intercell-potential VIC (see NOTE above)
C-----------------------------------------------------------------------
C
                     VIC(1:IRCRIT) = RPWL(1:IRCRIT)*AC
C
                     DO IR = 1,IRCRIT
                        VMADQ(IR,LM,IQ) = VMADQ(IR,LM,IQ) + VIC(IR)
                     END DO
C
                  END IF
C=======================================================================
C
               END DO
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
            END DO
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
            DO IR = 1,IRCRIT
C
               VSPHERE(1:N_LEBGRID) = 0.0D0
               V2(:) = 0D0
C
C---------------------------------------------------------------- sphere
               DO I_LEBGRID = 1,N_LEBGRID
C
                  RVEC_I(1:3) = RHAT_LEBGRID(1:3,I_LEBGRID)*R(IR,IM)
C
                  IF ( I_LEBGRID.EQ.1 ) THEN
                     RVEC_J(1:3) = RQNBR_Q(1:3,1,IQ)
                     RABS_J = DNRM2(3,RVEC_J,1)
                     RHAT_J(1:3) = RVEC_J(1:3)/RABS_J
                     RVEC_I(1:3) = RHAT_J(1:3)*R(IR,IM)
                  END IF
C
C------------------------------------------------------------- neighbors
                  DO JQNBR = 1,NQNBR_Q(IQ)
                     JQ = IQ_QNBR_Q(JQNBR,IQ)
                     JM = IMQ(JQ)
                     IRCRIT = JRCUT(NPAN(JM),JM)
                     RCRIT_J = R(IRCRIT,JM)
C
                     RVEC_J(1:3) = -RQNBR_Q(1:3,JQNBR,IQ)
     &                             *ALAT + RVEC_I(1:3)
                     RABS_J = DNRM2(3,RVEC_J,1)
                     RHAT_J(1:3) = RVEC_J(1:3)/RABS_J
C
                     CALL CALC_RHPLM(RHAT_J(1),RHAT_J(2),RHAT_J(3),RYLM,
     &                               NLFPMAX-1,NLMFPMAX)
C
                     IF ( RABS_J.LT.RCRIT_J ) THEN
C
                        JRB = 2
                        DO WHILE ( R(JRB,JM).LT.(RABS_J-1D-8) )
                           JRB = JRB + 1
                           IF ( JRB.GT.IRCRIT ) STOP 'in <FPVNFCOR> JRB'
                        END DO
                        JRA = JRB - 1
                        WR = 1D0/(R(JRB,JM)-R(JRA,JM))
                        WJB = (RABS_J-R(JRA,JM))*WR
                        WJA = (R(JRB,JM)-RABS_J)*WR
                        IF ( ABS(1-WJA-WJB).GT.1D-8 ) WRITE (*,*) WJA,
     &                       WJB
                        IF ( RABS_J-1D-8.LT.R(JRA,JM) ) WRITE (*,*) 'A',
     &                       RABS_J,R(JRA,JM)
                        IF ( RABS_J-1D-8.GT.R(JRB,JM) ) WRITE (*,*) 'B',
     &                       RABS_J,R(JRB,JM)
C
                        DO LM = 1,NLMFP
                           IF ( KFP_LMQ(LM,JQ).NE.0 ) THEN
C
                              VLMJ = WJA*VINTRAQ(JRA,LM,JQ)
     &                               + WJB*VINTRAQ(JRB,LM,JQ)
C
                              VSPHERE(I_LEBGRID) = VSPHERE(I_LEBGRID)
     &                           + VLMJ*RYLM(LM)
C
                           END IF
                        END DO
C
                     ELSE
C
                        DO LM = 1,NLMFP
                           IF ( KFP_LMQ(LM,JQ).NE.0 ) THEN
C
                              VLMJ = VINTRAQ(IRCRIT,LM,JQ)
     &                               *(RCRIT_J/RABS_J)**(L_LM(LM)+1)
C
                              VSPHERE(I_LEBGRID) = VSPHERE(I_LEBGRID)
     &                           + VLMJ*RYLM(LM)
C
                           END IF
                        END DO
C
                     END IF
C
                     DO LM = 1,NLMFP
                        IF ( KFP_LMQ(LM,JQ).NE.0 ) THEN
C
                           VLMJ = VINTRAQ(IRCRIT,LM,JQ)*(RCRIT_J/RABS_J)
     &                            **(L_LM(LM)+1)
C
                           V2(I_LEBGRID) = V2(I_LEBGRID) + VLMJ*RYLM(LM)
C
                        END IF
                     END DO
C
                  END DO
C------------------------------------------------------------- neighbors
               END DO
C---------------------------------------------------------------- sphere
C
               DO LM = 1,NLMFP
                  IF ( KFP_LMQ(LM,IQ).NE.0 ) VNFCORQ1(IR,LM,IQ)
     &                 = DDOT(I_LEBGRID,VSPHERE,1,W_LEBGRID(1,LM),1)
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C                write results for comparison
C-----------------------------------------------------------------------
            DO LM = 1,NLMFP
               IF ( KFP_LMQ(LM,IQ).NE.0 ) THEN
                  WRITE (*,*) '############### L LM',L_LM(LM),LM
C
                  WRITE (FILNAM(1:11),'("zablo",I3.3,I3.3)') IQ,LM
                  CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:11))
C
                  DO IR = 100,IRCRIT
C
                     RAT = VNFCORQ1(IR,LM,IQ)/VMADQ(IR,LM,IQ)
                     RAT = MIN(+1.1D0,RAT)
                     RAT = MAX(-1.1D0,RAT)
C
                     WRITE (IOTMP,99002) R(IR,IM),VMADQ(IR,LM,IQ),
     &                      VNFCORQ1(IR,LM,IQ),VNFCORQ1(IR,LM,IQ)
     &                      - VMADQ(IR,LM,IQ),RAT
                  END DO
                  CLOSE (IOTMP)
               END IF
            END DO
C-----------------------------------------------------------------------
C
         END IF
      END DO
C=======================================================================
C
      CALL CPU_TIME(TIME)
      WRITE (6,*) 'Time used for Nearfield Correction ',TIME - TIME0
C
      STOP '<FPVNFCOR>:   check completed '
C
99001 FORMAT (10X,'maximum for r_crit = ',F10.6)
99002 FORMAT (10E20.12)
99003 FORMAT (4X,2I3,F10.6,3X,3I3,F10.6,2F12.6)
99004 FORMAT (//,1X,79('*'),/,33X,'<FPVNFCOR>',/,1X,79('*'),//,10X,
     &        'set up near field corrections to Coulomb potential',/)
99005 FORMAT (/,'     IQ IM    rci     JNBR JQ JM    rcj',
     &        '         d_ij     rci + rcj')
99006 FORMAT (/,I7,3X,'number of neighbors ',I5)
      END
