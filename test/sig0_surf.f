C*==sig0_surf.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG0_SURF(S10AQAB,S10BQAB,S10DQAB,S2AQAB,S2BQAB,S2DQAB,
     &                     S3AQAB,S3BQAB,S3DQAB,S4AQAB,S4BQAB,S4DQAB,
     &                     MREG_TAB,MREG_TBA,UMAT_VTA,MSSTA,TSSTA,MSSQA,
     &                     TAUQA_TMP,MSSTB,TSSTB,UMAT_VTB,MSSQB,
     &                     TAUQB_TMP,SIG0Q,NSPINPROJ,SIG_MODE)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *                J_m Im TAU00(t) J_n Im TAU00(t)                   *
C   *                                                                  *
C   * the prefactor 1/4 is included here in the S-terms (Bastin)       *
C   * constants are added in <SIG_SUM> when converting to SI units     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_FILES,ONLY:IPRINT,DATSET,LDATSET,IOTMP
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,WKM2,IPIVKM
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI,ITOQ,NOQ
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SOTSI,EESI
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_THERMAL,ONLY:NFT,FMAT_FT,MFMAT_FT,NVIBFLU,NVIBRA,NFLUCT,
     &    NVFTMAX,X_VFT,NVFO_Q,IVFT_VFOQ,NVT
      USE MOD_CALCMODE,ONLY:PUBLIC_VERSION
      IMPLICIT NONE
C*--SIG0_SURF28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG0_SURF')
      REAL*8 TOL
      PARAMETER (TOL=1D-12)
C
C Dummy arguments
C
      INTEGER NSPINPROJ
      CHARACTER*10 SIG_MODE
      COMPLEX*16 MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX),
     &           MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX),
     &           MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           S10AQAB(3,3,NSPINPROJ,NQMAX),
     &           S10BQAB(3,3,NSPINPROJ,NQMAX),
     &           S10DQAB(3,3,NSPINPROJ,NQMAX),
     &           S2AQAB(3,3,NSPINPROJ,NQMAX),S2BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S2DQAB(3,3,NSPINPROJ,NQMAX),
     &           S3AQAB(3,3,NSPINPROJ,NQMAX),S3BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S3DQAB(3,3,NSPINPROJ,NQMAX),
     &           S4AQAB(3,3,NSPINPROJ,NQMAX),S4BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4DQAB(3,3,NSPINPROJ,NQMAX),SIG0Q(3,3,NSPINPROJ,NQMAX)
     &           ,TAUQA_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TAUQB_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TSSTA(NKMMAX,NKMMAX,NTMAX),TSSTB(NKMMAX,NKMMAX,NTMAX),
     &           UMAT_VTA(NKMMAX,NKMMAX,NVT),UMAT_VTB(NKMMAX,NKMMAX,NVT)
C
C Local variables
C
      REAL*8 DUM1,PRESI,PRE_SIGOFF,X
      CHARACTER*80 FILNAM,SIUNITS
      INTEGER I,IA_ERR,IC,IFLUCT,IFT,IO,IQ,ISPR,ISPRX,IT,IVFT,IVIBRA,
     &        IVT,IXX,IXXX,J,M,MUE,N,NUE
      COMPLEX*16 IMTAU_VFTA(:,:),IMTAU_VFTB(:,:),IMTSS_VFTA(:,:),
     &           MIRR_DUMM(:,:,:,:,:),MIRR_FTAB(:,:,:,:,:),
     &           MIRR_FTBA(:,:,:,:,:),MIRR_TAB(:,:,:,:,:),
     &           MIRR_TBA(:,:,:,:,:),MIRR_VFTAB(:,:,:,:,:),
     &           MIRR_VFTBA(:,:,:,:,:),MREG_FTAB(:,:,:,:),
     &           MREG_FTBA(:,:,:,:),MREG_VFTAB(:,:,:,:),
     &           MREG_VFTBA(:,:,:,:),MSS_VFTA(:,:),MSS_VFTB(:,:),
     &           REBS_VFTB(:,:),RETAU_VFTB(:,:),RETSS_VFTB(:,:),
     &           SIG0II(3,3),SIG0II_VFT(:,:,:),SIG0IR(3,3),
     &           SIG0IRSS_SS(:,:),SIG0IRSS_SS_VFT(:,:,:,:),
     &           SIG0IR_BS(:,:),SIG0IR_BS_VFT(:,:,:),SIG0IR_IRR(3,3),
     &           SIG0IR_IRR_VFT(:,:,:),SIG0IR_SS(:,:),
     &           SIG0IR_SS_VFT(:,:,:,:),SIG0IR_VFT(:,:,:),SUM_A(3,3),
     &           SUM_B(3,3),SUM_D(3,3),SUM_TOT(3,3),TAU_VFTA(:,:),
     &           TAU_VFTB(:,:),TRACE,TSS_FTA(:,:),TSS_FTB(:,:),
     &           TSS_VFTA(:,:),TSS_VFTB(:,:)
      LOGICAL LPRSIG
      CHARACTER*1 MSTR,NSTR
      CHARACTER*5 TXTVAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IMTAU_VFTA,IMTAU_VFTB,RETAU_VFTB
      ALLOCATABLE SIG0II_VFT,SIG0IR_VFT,SIG0IR_IRR_VFT
      ALLOCATABLE REBS_VFTB,RETSS_VFTB
      ALLOCATABLE SIG0IR_BS,SIG0IR_SS,SIG0IR_BS_VFT,SIG0IR_SS_VFT
      ALLOCATABLE SIG0IRSS_SS,SIG0IRSS_SS_VFT
      ALLOCATABLE IMTSS_VFTA,TAU_VFTA,TAU_VFTB
      ALLOCATABLE MSS_VFTA,MSS_VFTB,TSS_FTA,TSS_FTB,TSS_VFTA,TSS_VFTB
      ALLOCATABLE MREG_FTAB,MREG_FTBA,MREG_VFTAB,MREG_VFTBA
      ALLOCATABLE MIRR_TAB,MIRR_FTAB,MIRR_VFTAB,MIRR_DUMM
      ALLOCATABLE MIRR_TBA,MIRR_FTBA,MIRR_VFTBA
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      ALLOCATE (MIRR_TAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_TBA(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_DUMM(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_FTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_FTBA(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_VFTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_VFTBA(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MREG_FTAB(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MREG_FTBA(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MREG_VFTAB(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MREG_VFTBA(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MSS_VFTA(M,M),MSS_VFTB(M,M))
      ALLOCATE (TSS_FTA(M,M),TSS_FTB(M,M))
      ALLOCATE (TSS_VFTA(M,M),TSS_VFTB(M,M))
      ALLOCATE (TAU_VFTA(M,M),TAU_VFTB(M,M))
      ALLOCATE (SIG0IR_BS(3,3),SIG0IR_SS(3,3))
      ALLOCATE (SIG0IR_BS_VFT(3,3,NVFTMAX),SIG0IR_SS_VFT(3,3,NVFTMAX,2))
      ALLOCATE (SIG0IRSS_SS(3,3),SIG0IRSS_SS_VFT(3,3,NVFTMAX,2))
      ALLOCATE (SIG0II_VFT(3,3,NVFTMAX),SIG0IR_VFT(3,3,NVFTMAX))
      ALLOCATE (SIG0IR_IRR_VFT(3,3,NVFTMAX))
      ALLOCATE (REBS_VFTB(M,M),RETSS_VFTB(M,M))
      ALLOCATE (IMTAU_VFTA(M,M),IMTAU_VFTB(M,M))
      ALLOCATE (IMTSS_VFTA(M,M))
      ALLOCATE (RETAU_VFTB(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc IMTAUT')
C
      LPRSIG = .FALSE.
C      LPRSIG = .TRUE.
C
      TXTVAR = 'sigma'
C
      IF ( SIG_MODE.EQ.'INTEGRAL  ' ) THEN
         PRE_SIGOFF = 1D0
      ELSE
         PRE_SIGOFF = 0.5D0
      END IF
C
      WRITE (6,99002)
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO ISPR = 1,NSPR
         IXX = LIST_ISPR(ISPR)
C
         IF ( IXX.EQ.8 ) THEN
C---------- torkance prefactor for SI output (C * m)
            PRESI = SOTSI
            SIUNITS = '10**-30 C*m'
         ELSE IF ( IXX.EQ.10 ) THEN
C---------- Edelstein prefactor for SI output (m/V)
            PRESI = EESI
            SIUNITS = 'm/V'
         ELSE
C---------- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
            PRESI = 1D-8*CONSI
            SIUNITS = '1/(muOhm*cm)'
         END IF
C
         WRITE (6,99014)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(IXX)
         WRITE (6,99014)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
C
            N = NKM
C
            WRITE (6,99003) IQ
C
            SUM_A(1:3,1:3) = S10AQAB(1:3,1:3,IXX,IQ)
     &                       + S2AQAB(1:3,1:3,IXX,IQ)
     &                       + S3AQAB(1:3,1:3,IXX,IQ)
     &                       + S4AQAB(1:3,1:3,IXX,IQ)
C
            SUM_B(1:3,1:3) = S10BQAB(1:3,1:3,IXX,IQ)
     &                       + S2BQAB(1:3,1:3,IXX,IQ)
     &                       + S3BQAB(1:3,1:3,IXX,IQ)
     &                       + S4BQAB(1:3,1:3,IXX,IQ)
C
            SUM_D(1:3,1:3) = S10DQAB(1:3,1:3,IXX,IQ)
     &                       + S2DQAB(1:3,1:3,IXX,IQ)
     &                       + S3DQAB(1:3,1:3,IXX,IQ)
     &                       + S4DQAB(1:3,1:3,IXX,IQ)
C
            SUM_A(1:3,1:3) = 0.25D0*SUM_A(1:3,1:3)
            SUM_B(1:3,1:3) = 0.25D0*SUM_B(1:3,1:3)
            SUM_D(1:3,1:3) = 0.25D0*SUM_D(1:3,1:3)
C
            SUM_TOT(1:3,1:3) = SUM_B(1:3,1:3) - SUM_D(1:3,1:3)
     &                         - SUM_A(1:3,1:3) + SUM_B(1:3,1:3)
C
            WRITE (6,99017)
            WRITE (6,99008) ((SUM_TOT(MUE,NUE),NUE=1,3),MUE=1,3)
C
            IF ( IPRINT.GT.4 ) THEN
C
               X = 0.25D0
               WRITE (6,99018) 'S10AQAB',IQ
               CALL PR_COND_TENSOR(S10AQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S2AQAB',IQ
               CALL PR_COND_TENSOR(S2AQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S3AQAB',IQ
               CALL PR_COND_TENSOR(S3AQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S4AQAB',IQ
               CALL PR_COND_TENSOR(S4AQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S10BQAB',IQ
               CALL PR_COND_TENSOR(S10BQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S2BQAB',IQ
               CALL PR_COND_TENSOR(S2BQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S3BQAB',IQ
               CALL PR_COND_TENSOR(S3BQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S4BQAB',IQ
               CALL PR_COND_TENSOR(S4BQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S10DQAB',IQ
               CALL PR_COND_TENSOR(S10DQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S2DQAB',IQ
               CALL PR_COND_TENSOR(S2DQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S3DQAB',IQ
               CALL PR_COND_TENSOR(S3DQAB(1:3,1:3,IXX,IQ),X)
               WRITE (6,99018) 'S4DQAB',IQ
               CALL PR_COND_TENSOR(S4DQAB(1:3,1:3,IXX,IQ),X)
C
               WRITE (6,*) ' '
               WRITE (6,99018) 'S10',IQ
               SUM_D = S10BQAB(1:3,1:3,IXX,IQ) - S10DQAB(1:3,1:3,IXX,IQ)
     &                 - S10AQAB(1:3,1:3,IXX,IQ)
     &                 + S10BQAB(1:3,1:3,IXX,IQ)
               CALL PR_COND_TENSOR(SUM_D,X)
C
               WRITE (6,99018) 'S2',IQ
               SUM_D = S2BQAB(1:3,1:3,IXX,IQ) - S2DQAB(1:3,1:3,IXX,IQ)
     &                 - S2AQAB(1:3,1:3,IXX,IQ) + S2BQAB(1:3,1:3,IXX,IQ)
               CALL PR_COND_TENSOR(SUM_D,X)
C
               WRITE (6,99018) 'S3',IQ
               SUM_D = S3BQAB(1:3,1:3,IXX,IQ) - S3DQAB(1:3,1:3,IXX,IQ)
     &                 - S3AQAB(1:3,1:3,IXX,IQ) + S3BQAB(1:3,1:3,IXX,IQ)
               CALL PR_COND_TENSOR(SUM_D,X)
C
               WRITE (6,99018) 'S4',IQ
               SUM_D = S4BQAB(1:3,1:3,IXX,IQ) - S4DQAB(1:3,1:3,IXX,IQ)
     &                 - S4AQAB(1:3,1:3,IXX,IQ) + S4BQAB(1:3,1:3,IXX,IQ)
               CALL PR_COND_TENSOR(SUM_D,X)
            END IF
C
C--------------------------------- check if imaginary part of SIG0Q is 0
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(DIMAG(SUM_TOT(MUE,NUE))).GT.1D-5 )
     &                 WRITE (6,99007) DIMAG(SUM_TOT(MUE,NUE)),MUE,NUE,
     &                                 IQ
               END DO
            END DO
C-----------------------------------------------------------------------
C
            CALL CINIT(3*3,SIG0II)
            CALL CINIT(3*3,SIG0IR)
            CALL CINIT(3*3,SIG0IR_IRR)
C
            CALL CINIT(3*3,SIG0IR_BS)
            CALL CINIT(3*3,SIG0IR_SS)
            CALL CINIT(3*3*NVFTMAX*2,SIG0IR_SS_VFT)
            CALL CINIT(3*3,SIG0IRSS_SS)
            CALL CINIT(3*3*NVFTMAX*2,SIG0IRSS_SS_VFT)
            CALL CINIT(3*3*NVFTMAX,SIG0IR_BS_VFT)
C
            CALL CINIT(3*3*NVFTMAX,SIG0II_VFT)
            CALL CINIT(3*3*NVFTMAX,SIG0IR_VFT)
            CALL CINIT(3*3*NVFTMAX,SIG0IR_IRR_VFT)
C
C================================================================= IO ==
            LOOP_IO:DO IO = 1,NOQ(IQ)
C
               IT = ITOQ(IO,IQ)
C
               MIRR_TAB(:,:,:,:,:) = 0D0
               MIRR_TBA(:,:,:,:,:) = 0D0
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
               IFT = (IT-1)*NFLUCT
               LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
                  IF ( NFLUCT.EQ.1 ) THEN
C
                     TSS_FTA(:,:) = TSSTA(:,:,IT)
                     TSS_FTB(:,:) = TSSTB(:,:,IT)
C
                     DO ISPRX = 1,NSPR
                        IXXX = LIST_ISPR(ISPRX)
                        MREG_FTAB(:,:,:,IXXX) = MREG_TAB(:,:,:,IXXX,IT)
                        MREG_FTBA(:,:,:,IXXX) = MREG_TBA(:,:,:,IXXX,IT)
                        MIRR_FTAB(:,:,:,:,IXXX) = MIRR_TAB(:,:,:,:,IXXX)
                        MIRR_FTBA(:,:,:,:,IXXX) = MIRR_TBA(:,:,:,:,IXXX)
                     END DO
C
                  ELSE
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSSTA(1,1,IT),'SAS+',TSS_FTA)
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSSTB(1,1,IT),'SAS+',TSS_FTB)
C
                     DO ISPRX = 1,NSPR
                        IXXX = LIST_ISPR(ISPRX)
                        CALL ME_ROTATE_REG(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                     MREG_TAB(1,1,1,IXXX,IT),
     &                     MREG_TBA(1,1,1,IXXX,IT),MREG_FTAB(1,1,1,IXXX)
     &                     ,MREG_FTBA(1,1,1,IXXX))
C
                        CALL ME_ROTATE_IRR(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                     MIRR_TAB(1,1,1,1,IXXX),MIRR_TBA(1,1,1,1,IXXX)
     &                     ,MIRR_TBA(1,1,1,1,IXXX),
     &                     MIRR_FTAB(1,1,1,1,IXXX),
     &                     MIRR_FTBA(1,1,1,1,IXXX),
     &                     MIRR_DUMM(1,1,1,1,IXXX))
                     END DO
C
                  END IF
C
C============================================================= IVIBRA ==
C                                             perform local displacement
                  IVT = (IT-1)*NVIBRA
                  LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
C
                     IVT = IVT + 1
C
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                     IF ( NVIBRA.EQ.1 ) THEN
C
                        TSS_VFTA(:,:) = TSS_FTA(:,:)
                        TSS_VFTB(:,:) = TSS_FTB(:,:)
C
                        DO ISPRX = 1,NSPR
                           IXXX = LIST_ISPR(ISPRX)
                           MREG_VFTAB(:,:,:,IXXX)
     &                        = MREG_FTAB(:,:,:,IXXX)
                           MREG_VFTBA(:,:,:,IXXX)
     &                        = MREG_FTBA(:,:,:,IXXX)
                           MIRR_VFTAB(:,:,:,:,IXXX)
     &                        = MIRR_FTAB(:,:,:,:,IXXX)
                           MIRR_VFTBA(:,:,:,:,IXXX)
     &                        = MIRR_FTBA(:,:,:,:,IXXX)
                        END DO
C
                     ELSE
C
                        CALL CMAT_U_TRANS(UMAT_VTA(1,1,IVT),TSS_FTA,
     &                     'UAUT',TSS_VFTA)
C
                        CALL CMAT_U_TRANS(UMAT_VTB(1,1,IVT),TSS_FTA,
     &                     'UAUT',TSS_VFTB)
C
                        DO ISPRX = 1,NSPR
                           IXXX = LIST_ISPR(ISPRX)
                           CALL ME_UTRANS_REG(UMAT_VTA(1,1,IVT),
     &                        UMAT_VTB(1,1,IVT),MREG_FTAB(1,1,1,IXXX),
     &                        MREG_FTBA(1,1,1,IXXX),
     &                        MREG_VFTAB(1,1,1,IXXX),
     &                        MREG_VFTBA(1,1,1,IXXX))
C
                           CALL ME_UTRANS_IRR(UMAT_VTA(1,1,IVT),
     &                        UMAT_VTB(1,1,IVT),MIRR_FTAB(1,1,1,1,IXXX),
     &                        MIRR_FTBA(1,1,1,1,IXXX),
     &                        MIRR_FTBA(1,1,1,1,IXXX),
     &                        MIRR_VFTAB(1,1,1,1,IXXX),
     &                        MIRR_VFTBA(1,1,1,1,IXXX),
     &                        MIRR_DUMM(1,1,1,1,IXXX))
                        END DO
C
                     END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                     IF ( NVIBFLU.EQ.1 ) THEN
                        MSS_VFTA(:,:) = MSSTA(:,:,IT)
                        MSS_VFTB(:,:) = MSSTB(:,:,IT)
                     ELSE
                        CALL CMATINV3(N,M,IPIVKM,TSS_VFTA,WKM1,MSS_VFTA)
                        CALL CMATINV3(N,M,IPIVKM,TSS_VFTB,WKM1,MSS_VFTB)
                     END IF
C
C-----------------------------------------------------------------------
C
                     CALL GET_TAUT(MSS_VFTA,MSSQA(1,1,IQ),
     &                             TAUQA_TMP(1,1,IQ),TAU_VFTA)
C
                     CALL GET_TAUT(MSS_VFTB,MSSQB(1,1,IQ),
     &                             TAUQB_TMP(1,1,IQ),TAU_VFTB)
C
C-----------------------------------------------------------------------
C
                     DO J = 1,N
                        DO I = 1,N
                           IMTAU_VFTA(I,J)
     &                        = (TAU_VFTA(I,J)-DCONJG(TAU_VFTA(J,I)))
     &                        /(2D0*CI)
                           IMTAU_VFTB(I,J)
     &                        = (TAU_VFTB(I,J)-DCONJG(TAU_VFTB(J,I)))
     &                        /(2D0*CI)
                           RETAU_VFTB(I,J)
     &                        = (TAU_VFTB(I,J)+DCONJG(TAU_VFTB(J,I)))
     &                        /2D0
C
                           REBS_VFTB(I,J)
     &                        = ((TAU_VFTB(I,J)-TSS_VFTB(I,J))
     &                        +DCONJG(TAU_VFTB(J,I)-TSS_VFTB(I,J)))/2D0
                           RETSS_VFTB(I,J)
     &                        = (TSS_VFTB(I,J)+DCONJG(TSS_VFTB(J,I)))
     &                        /2D0
C
                           IMTSS_VFTA(I,J)
     &                        = (TSS_VFTA(I,J)-DCONJG(TSS_VFTA(J,I)))
     &                        /(2D0*CI)
                        END DO
                     END DO
C
                     LOOP_MUE:DO MUE = 1,3
                        LOOP_NUE:DO NUE = 1,3
C
C----------------------------------------------------- j_m Im G j_n Im G
C
                           CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                                  IMTAU_VFTB,WKM1)
C
                           CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                           CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                                  WKM2,WKM1)
C
                           DO I = 1,NKM
                              SIG0II(MUE,NUE) = SIG0II(MUE,NUE)
     &                           + X_VFT(IVFT)*WKM1(I,I)
C
                              SIG0II_VFT(MUE,NUE,IVFT)
     &                           = SIG0II_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*WKM1(I,I)
                           END DO
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
                           IF ( .NOT.PUBLIC_VERSION ) THEN
C
                              IF ( MUE.LT.0 )
     &                             CALL SIG0_SURF_NONPUB(MIRR_VFTAB,
     &                             MIRR_VFTBA,MREG_VFTAB,MREG_VFTBA,
     &                             REBS_VFTB,IMTAU_VFTA,IMTSS_VFTA,
     &                             RETAU_VFTB,RETSS_VFTB,SIG0IR_BS,
     &                             SIG0IR_SS,SIG0IR_BS_VFT,
     &                             SIG0IR_SS_VFT,SIG0IRSS_SS,SIG0IR_IRR,
     &                             SIG0IRSS_SS_VFT,SIG0IR_VFT,
     &                             SIG0IR_IRR_VFT,SIG0IR,N,M,MUE,NUE,
     &                             IVFT,IXX,NSPINPROJ)
C
C
C-------------------------------- i ( J_m Im G j_n - j_n Im G J_m ) Re G
C-----------------------------------------------------------regular part
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           RETAU_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           WKM2,WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           RETAU_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           WKM2,WKM1)
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IR(MUE,NUE) = SIG0IR(MUE,NUE)
     &                           + X_VFT(IVFT)*CI*TRACE
C
                              SIG0IR_VFT(MUE,NUE,IVFT)
     &                           = SIG0IR_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*CI*TRACE
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C---------------------------------------------------------irregular part
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,
     &                           MIRR_VFTBA(1,1,NUE,MUE,IXX),WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,
     &                           MIRR_VFTAB(1,1,MUE,NUE,IXX),WKM1)
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IR_IRR(MUE,NUE) = SIG0IR_IRR(MUE,NUE)
     &                           + X_VFT(IVFT)*TRACE*CI
C
                              SIG0IR_IRR_VFT(MUE,NUE,IVFT)
     &                           = SIG0IR_IRR_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*TRACE*CI
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-------------------------------- here some tools for analysis follow
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C-----------------------------------------------------------------------
C-------------------------------------------- single site contribution 1
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-----------------------------------------------------------------------
C----------------------------------------------------Z tau Z * ZJ   term
C
                              SIG0IR_SS(MUE,NUE) = SIG0IR_SS(MUE,NUE)
     &                           - X_VFT(IVFT)*TRACE*CI
                              SIG0IR_SS_VFT(MUE,NUE,IVFT,1)
     &                           = SIG0IR_SS_VFT(MUE,NUE,IVFT,1)
     &                           - X_VFT(IVFT)*TRACE*CI
C----------------------------------------------------Z tau Z * ZtZ  term
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           RETSS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           WKM2,WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           RETSS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           WKM2,WKM1)
C
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IR_SS(MUE,NUE) = SIG0IR_SS(MUE,NUE)
     &                           + X_VFT(IVFT)*CI*TRACE
                              SIG0IR_SS_VFT(MUE,NUE,IVFT,2)
     &                           = SIG0IR_SS_VFT(MUE,NUE,IVFT,2)
     &                           + X_VFT(IVFT)*TRACE*CI
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C------------------------------------------ back scattering contribution
C-------------------------------------------- Z tau Z * Z(tau-t)Z   term
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           REBS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           WKM2,WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           REBS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           WKM2,WKM1)
C
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IR_BS(MUE,NUE) = SIG0IR_BS(MUE,NUE)
     &                           + X_VFT(IVFT)*CI*TRACE
                              SIG0IR_BS_VFT(MUE,NUE,IVFT)
     &                           = SIG0IR_BS_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*CI*TRACE
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-------------------------------------------- single site contribution 2
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C----------------------------------------------------Z t   Z * ZtZ  term
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),
     &                           RETSS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTSS_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           WKM2,WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,MREG_VFTAB(1,1,MUE,IXX),
     &                           RETSS_VFTB,WKM1)
C
                              CALL CMATMUL(N,M,IMTSS_VFTA,WKM1,WKM2)
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,NUE,1),
     &                           WKM2,WKM1)
C
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IRSS_SS(MUE,NUE)
     &                           = SIG0IRSS_SS(MUE,NUE) + X_VFT(IVFT)
     &                           *CI*TRACE
                              SIG0IRSS_SS_VFT(MUE,NUE,IVFT,2)
     &                           = SIG0IRSS_SS_VFT(MUE,NUE,IVFT,2)
     &                           + X_VFT(IVFT)*TRACE*CI
C----------------------------------------------------Z t   Z * ZJ   term
C
                              CALL CMATMUL(N,M,IMTSS_VFTA,
     &                           MIRR_VFTBA(1,1,NUE,MUE,IXX),WKM1)
C
                              TRACE = C0
                              DO I = 1,NKM
                                 TRACE = TRACE + WKM1(I,I)
                              END DO
C
                              CALL CMATMUL(N,M,IMTSS_VFTA,
     &                           MIRR_VFTAB(1,1,MUE,NUE,IXX),WKM1)
C
                              DO I = 1,NKM
                                 TRACE = TRACE - WKM1(I,I)
                              END DO
C
                              SIG0IRSS_SS(MUE,NUE)
     &                           = SIG0IRSS_SS(MUE,NUE) - X_VFT(IVFT)
     &                           *TRACE*CI
                              SIG0IRSS_SS_VFT(MUE,NUE,IVFT,1)
     &                           = SIG0IRSS_SS_VFT(MUE,NUE,IVFT,1)
     &                           - X_VFT(IVFT)*TRACE*CI
C
                           END IF
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C-----------------------------------------------------------------------
C
                        END DO LOOP_NUE
C-----------------------------------------------------------------------
C
                     END DO LOOP_MUE
C-----------------------------------------------------------------------
C
                  END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
               END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
            END DO LOOP_IO
C================================================================= IO ==
C
C
            IF ( IPRINT.GT.4 ) THEN
               WRITE (6,99018) 'SIG0II',IQ
               CALL PR_COND_TENSOR(SIG0II(1:3,1:3),1D0)
               WRITE (6,99001) 'SIG0IR',IQ
               CALL PR_COND_TENSOR(SIG0IR(1:3,1:3),1D0)
               WRITE (6,99001) 'SIG0IR_IRR',IQ
               CALL PR_COND_TENSOR(SIG0IR_IRR(1:3,1:3),1D0)
            END IF
C
C ----------------------------------- suppress small elements and sum up
C
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(SIG0II(MUE,NUE)).LT.TOL ) SIG0II(MUE,NUE)
     &                 = C0
C
                  IF ( ABS(DIMAG(SIG0II(MUE,NUE))).LT.TOL )
     &                 SIG0II(MUE,NUE) = DREAL(SIG0II(MUE,NUE))
C
                  IF ( ABS(SIG0IR(MUE,NUE)).LT.TOL ) SIG0IR(MUE,NUE)
     &                 = C0
C
                  IF ( ABS(DIMAG(SIG0IR(MUE,NUE))).LT.TOL )
     &                 SIG0IR(MUE,NUE) = DREAL(SIG0IR(MUE,NUE))
C
                  IF ( ABS(SIG0IR_IRR(MUE,NUE)).LT.TOL )
     &                 SIG0IR_IRR(MUE,NUE) = C0
C
                  IF ( ABS(DIMAG(SIG0IR_IRR(MUE,NUE))).LT.TOL )
     &                 SIG0IR_IRR(MUE,NUE) = DREAL(SIG0IR_IRR(MUE,NUE))
               END DO ! NUE
            END DO ! MUE
C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
C---------------anti-symmetric contribution neglected (i.e. no AHE, SHE)
C - only 1/2 of SIG0IR and SIG0IR_IRR -> Crepieux formula PRB 64, 014416
            SIG0Q(1:3,1:3,IXX,IQ) = SIG0II(1:3,1:3)
     &                              + PRE_SIGOFF*SIG0IR(1:3,1:3)
     &                              - PRE_SIGOFF*SIG0IR_IRR(1:3,1:3)
C
C
            WRITE (6,99001) 'SIG0Q ALT',IQ
            CALL PR_COND_TENSOR(SIG0Q(1:3,1:3,IXX,IQ),1D0)
C
            SIG0Q(1:3,1:3,IXX,IQ) = SUM_TOT(1:3,1:3)
C
            WRITE (6,99001) 'SIG0Q NEU',IQ
            CALL PR_COND_TENSOR(SIG0Q(1:3,1:3,IXX,IQ),1D0)
C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
C --------------------------------- write out site-diagonal term sigma_0
C
            IF ( IPRINT.GE.3 ) THEN
C
               WRITE (6,99004) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (SIG0II(MUE,NUE),NUE=1,3)
               END DO
               WRITE (6,99005) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (SIG0IR(MUE,NUE),NUE=1,3)
               END DO
               WRITE (6,99010) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (SIG0IR_IRR(MUE,NUE),NUE=1,3)
               END DO
C-------------------------------------- printing BS and SS contributions
               IF ( IPRINT.GT.4 ) THEN
                  WRITE (6,99012) TXTVAR
                  DO MUE = 1,3
                     WRITE (6,99008) (SIG0IR_SS(MUE,NUE),NUE=1,3)
                  END DO
                  WRITE (6,99013) TXTVAR
                  DO MUE = 1,3
                     WRITE (6,99008) (SIG0IR_BS(MUE,NUE),NUE=1,3)
                  END DO
               END IF
C
C------------------------------------------------ printing type resolved
               IF ( IPRINT.GT.4 ) THEN
C
                  WRITE (6,'(/,10x,"Type resolved")')
C
                  DO IO = 1,NVFO_Q(IQ)
                     IVFT = IVFT_VFOQ(IO,IQ)
C
                     WRITE (6,'(/,10x,"Type IVFT = ",i5)') IVFT
C
                     WRITE (6,99004) TXTVAR
                     DO MUE = 1,3
                        WRITE (6,99008) (SIG0II_VFT(MUE,NUE,IVFT),NUE=1,
     &                                  3)
                     END DO
                     WRITE (6,99005) TXTVAR
                     DO MUE = 1,3
                        WRITE (6,99008) (SIG0IR_VFT(MUE,NUE,IVFT),NUE=1,
     &                                  3)
                     END DO
                     WRITE (6,99010) TXTVAR
                     DO MUE = 1,3
                        WRITE (6,99008) (SIG0IR_IRR_VFT(MUE,NUE,IVFT),
     &                                  NUE=1,3)
                     END DO
                  END DO ! IO
               END IF
C------------------------------------------------------------- SI  units
               WRITE (6,99011) TRIM(SIUNITS)
               WRITE (6,99004) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (PRESI*SIG0II(MUE,NUE),NUE=1,3)
               END DO
               WRITE (6,99005) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (PRESI*SIG0IR(MUE,NUE),NUE=1,3)
               END DO
               WRITE (6,99010) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (PRESI*SIG0IR_IRR(MUE,NUE),NUE=1,3)
               END DO
C
               WRITE (6,99012) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (PRESI*SIG0IR_SS(MUE,NUE),NUE=1,3)
               END DO
               WRITE (6,99013) TXTVAR
               DO MUE = 1,3
                  WRITE (6,99008) (PRESI*SIG0IR_BS(MUE,NUE),NUE=1,3)
               END DO
C
               IF ( IPRINT.GE.2 ) THEN
                  WRITE (6,99006) TXTVAR
                  DO MUE = 1,3
                     WRITE (6,99009) (DREAL(SIG0Q(MUE,NUE,IXX,IQ)),NUE=1
     &                               ,3)
                  END DO
               END IF
C
               WRITE (6,99011) TRIM(SIUNITS)
C
               DO MUE = 1,3
                  WRITE (6,99009) (PRESI*DREAL(SIG0Q(MUE,NUE,IXX,IQ)),
     &                            NUE=1,3)
               END DO
C
            END IF
C
C ----------------------------- check if imaginary part of SIG0Q is zero
C
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(DIMAG(SIG0Q(MUE,NUE,IXX,IQ))).GT.1D-5 )
     &                 WRITE (6,99007) DIMAG(SIG0Q(MUE,NUE,IXX,IQ)),MUE,
     &                                 NUE,IQ
               END DO
            END DO
C
C------------------------------------------ write conductivities to file
            IF ( IXX.EQ.1 .AND. LPRSIG .AND. NQ.EQ.1 .AND. IVFT.GT.1 )
     &           THEN
               DO MUE = 1,3
                  DO NUE = 1,3
                     WRITE (MSTR(1:1),'(i1)') MUE
                     WRITE (NSTR(1:1),'(i1)') NUE
                     FILNAM = DATSET(1:LDATSET)//MSTR//NSTR//'_SIG0.dat'
                     OPEN (IOTMP,FILE=FILNAM,STATUS='replace')
                     WRITE (IOTMP,99016)
C
                     DUM1 = PRESI*DREAL(SIG0II(MUE,NUE)
     &                      +PRE_SIGOFF*(SIG0IR(MUE,NUE)
     &                      -SIG0IR_IRR(MUE,NUE)))
C
                     WRITE (IOTMP,99015) X_VFT(1),DUM1,
     &                      (DREAL(PRESI*SIG0IR_BS_VFT(MUE,NUE,IVFT)),
     &                      (DREAL(PRESI*SIG0IR_SS_VFT(MUE,NUE,IVFT,IC))
     &                      ,IC=1,2),IVFT=1,2),
     &                      ((DREAL(PRESI*SIG0IRSS_SS_VFT(MUE,NUE,IVFT,
     &                      IC)),IC=1,2),IVFT=1,2)
                     CLOSE (IOTMP)
                  END DO ! NUE
               END DO ! MUE
            END IF
C-----------------------------------------------------------------------
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DEALLOCATE (SIG0IR_IRR_VFT,SIG0IR_VFT,SIG0II_VFT)
      DEALLOCATE (SIG0IR_BS,SIG0IR_SS)
      DEALLOCATE (IMTAU_VFTA,IMTAU_VFTB,RETAU_VFTB,IMTSS_VFTA,
     &            STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC IMTAU_VFTA')
C
99001 FORMAT (/,'#### <SIG0_SURF> ',A,' IQ:',I3,/,3(2X,2E12.5))
99002 FORMAT (//,1X,79('*'),/,34X,'<SIG0_SURF>',/,1X,79('*'),/)
99003 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99004 FORMAT (/,10X,'site-diagonal term ',A,'_0    j_m Im G j_n Im G',/)
99005 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',/)
99006 FORMAT (/,10X,'site-diagonal term ',A,'_0',/)
99007 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99008 FORMAT (3(F16.10,F16.10))
99009 FORMAT (10X,3F14.6)
99010 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',/,
     &        13X,'(contribution from irregular solution)',/)
99011 FORMAT (/,'    SI units [',A,']',/)
99012 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',
     &        ' - single-site contribution',/)
99013 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',
     &        ' - back-scattering contribution',/)
99014 FORMAT (/,37('*'),' SIG0_SURF ',37('*'))
99015 FORMAT (99E17.8)
99016 FORMAT ('#',5x,'conc       | ',1x,'0sum |',5x,
     &        'bs1_Z*tau*Z*Z(tau-t)Z | ss1_Z*tau*Z*ZJ|',
     &        'ss1_Z*tau*Z*ZtZ|bs2_Z*tau*Z*Z(tau-t)Z|',
     &        'ss2_Z*tau*Z*ZJ|ss2_Z*tau*Z*ZtZ|',
     &        'ss1_ZtZ*ZJ    |    ss1_ZtZ*ZtZ  |  ',
     &        'ss2_ZtZ*ZJ    |    ss2_ZtZ*ZtZ')
99017 FORMAT (/,10X,'site-diagonal term SIGMA_0   ',
     &        '(Fermi-surface contribution)',/)
99018 FORMAT (/,'#### <SIG0_SURF> ',A,' IQ:',I3)
      END
