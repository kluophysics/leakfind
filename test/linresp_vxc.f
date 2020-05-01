C*==linresp_vxc.f    processed by SPAG 6.70Rc at 21:33 on 19 Dec 2016
      SUBROUTINE LINRESP_VXC(RHO2NS,RGNT_VSF,LMRGNT_VSF,NRGNT_VSF_LM,
     &                       NRGNT_VSF,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &                       NMAX_LEBGRID)
C   ********************************************************************
C   *                                                                  *
C   *  this is an adapted version of the standard SCF routine  FPVXCLM *
C   *                                                                  *
C   *   to get the CHANGE in the XC-potential                          *
C   *   due to a CHANGE in the densities                               *
C   *   to be used within linear response calculations                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:RHO2NS_GG,NOBSE,NPERT,K_PERT_LMTP,
     &    HX_PERT_LMTP,GNTTAB
      USE MOD_SCF,ONLY:SCFVXC
      USE MOD_FILES,ONLY:
      USE MOD_RMESH,ONLY:NRSFTOT,NSF,LMISF,JRCRI,FLMSF,R,JRCUT,NPAN,
     &    NRMAX,R2DRDI,NSFMAX,KLMSF,ISFLM,NLSF
      USE MOD_TYPES,ONLY:IMT,NTMAX,NLFP,NLMFP,NLMFPT,NLMFPMAX,ITBOT,
     &    ITTOP
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,CONST_4PI,RY_EV
      USE MOD_ANGMOM,ONLY:L_LM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_VXC')
C
C Dummy arguments
C
      INTEGER NMAX_LEBGRID,NRGNT_VSF,N_LEBGRID
      INTEGER LMRGNT_VSF(NRGNT_VSF,3),NRGNT_VSF_LM(0:NLMFPMAX)
      REAL*8 RGNT_VSF(NRGNT_VSF),RHO2NS(NRMAX,NLMFPMAX,NTMAX,3),
     &       W_LEBGRID(NMAX_LEBGRID,NLMFPMAX),
     &       Y_LEBGRID(NMAX_LEBGRID,NLMFPMAX)
C
C Local variables
C
      REAL*8 AGRDRHO(:,:),AGRDRHOD(:,:),AGRDRHOU(:,:),D2FDR2(:),DFDR(:),
     &       EMM(:),ENM(:),ENN(:),FLM(:,:),FPIPR2,G123,GAMMA_LEBO(:,:),
     &       GAMMA_O(:,:,:),GAMMA_X_FLM(:,:),GDGAG(:,:),GDGAGD(:,:),
     &       GDGAGU(:,:),GFANG(:,:,:),GFLM(:,:,:),HANG(:,:),HLM(:,:),
     &       HX_LM(:,:),ICHR,I_STONER,LAPRHO(:,:),LAPRHOD(:,:),
     &       LAPRHOU(:,:),RHO(:,:),RHO4PI(:,:),RHOD(:,:),RHOU(:,:),
     &       RINT1(:),RPWM2(:),RWGT,VXC_LEBOP(:,:,:),
     &       VXC_LMTOP(:,:,:,:,:),X
      REAL*8 DDOT
      LOGICAL GGA
      INTEGER IM,IOBSE,IPAN1,IPERT,IR,IRCRIT,IRMTIN,IRSF,IRTOP,ISF,
     &        ISPIN,IT,J,LDV,LH,LM,LM1,LM2,LM3,LMDV,LMH,LMSF,LSF,MH,
     &        NANG,NLG,NR
      EXTERNAL DAXPY,DDOT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RHO4PI,GAMMA_O,GAMMA_LEBO
      ALLOCATABLE RHO,RHOD,RHOU
      ALLOCATABLE AGRDRHO,AGRDRHOD,AGRDRHOU
      ALLOCATABLE LAPRHO,LAPRHOD,LAPRHOU
      ALLOCATABLE GDGAG,GDGAGD,GDGAGU
      ALLOCATABLE GFLM,GFANG,FLM,HLM,HANG
      ALLOCATABLE RPWM2,DFDR,D2FDR2,EMM,ENM,ENN,RINT1
      ALLOCATABLE VXC_LEBOP,GAMMA_X_FLM,VXC_LMTOP,HX_LM
C
      ALLOCATE (EMM(NMAX_LEBGRID),ENM(NMAX_LEBGRID),ENN(NMAX_LEBGRID))
      ALLOCATE (RHO4PI(NMAX_LEBGRID,2))
      ALLOCATE (GAMMA_O(NRMAX,NLMFPMAX,NOBSE),RINT1(NRMAX))
      ALLOCATE (GAMMA_LEBO(NMAX_LEBGRID,NOBSE))
      ALLOCATE (VXC_LEBOP(NMAX_LEBGRID,NOBSE,NPERT))
      ALLOCATE (GAMMA_X_FLM(NRMAX,NLMFPMAX))
      ALLOCATE (VXC_LMTOP(NRMAX,NLMFPMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (HX_LM(NRMAX,NLMFPMAX))
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C                      GGA - parametrisations
C=======================================================================
      IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
         GGA = .TRUE.
C
         NANG = NMAX_LEBGRID
         NR = NRMAX
         ALLOCATE (RHO(NANG,NR),RHOD(NANG,NR),RHOU(NANG,NR))
         ALLOCATE (AGRDRHO(NANG,NR),AGRDRHOD(NANG,NR),AGRDRHOU(NANG,NR))
         ALLOCATE (LAPRHO(NANG,NR),LAPRHOD(NANG,NR),LAPRHOU(NANG,NR))
         ALLOCATE (GDGAG(NANG,NR),GDGAGD(NANG,NR),GDGAGU(NANG,NR))
         ALLOCATE (GFLM(NR,NLMFP,3),GFANG(NANG,NR,3))
         ALLOCATE (FLM(NR,NLMFP),HLM(NR,NLMFP),HANG(NANG,NR))
         ALLOCATE (RPWM2(NRMAX),DFDR(NRMAX),D2FDR2(NRMAX))
C
      ELSE
C
         GGA = .FALSE.
C
      END IF
C=======================================================================
C
      NLG = N_LEBGRID
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IPAN1 = NPAN(IM)
         IRTOP = JRCUT(IPAN1,IM)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCRI(IM)
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C                  get normalized induced densities GAMMA_O
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
         DO IOBSE = 1,NOBSE
            IPERT = IOBSE
C
C----- convolute RHO2NS_GG with shape functions to get spherical density
C
            DO IR = 1,IRMTIN
               RINT1(IR) = RHO2NS_GG(IR,1,IT,IOBSE,IPERT)
     &                     *SQRT_4PI*R2DRDI(IR,IM)
            END DO
C
            RINT1((IRMTIN+1):IRCRIT) = 0D0
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               IF ( LM.LE.NLMFPMAX ) THEN
C
                  DO IRSF = 1,NRSFTOT(IM)
                     IR = IRSF + IRMTIN
C
                     RINT1(IR) = RINT1(IR)
     &                           + RHO2NS_GG(IR,LM,IT,IOBSE,IPERT)
     &                           *FLMSF(IRSF,ISF,IM)*R2DRDI(IR,IM)
                  END DO
C
               END IF
            END DO
C
C ---------------------------------------------------- integrate density
C
            CALL RRADINT(IM,RINT1,ICHR)
C
C ------------------------------------------------------- normalize to 1
C
            GAMMA_O(1:IRCRIT,1:NLMFPMAX,IOBSE)
     &         = RHO2NS_GG(1:IRCRIT,1:NLMFPMAX,IT,IOBSE,IPERT)/ICHR
C
            WRITE (6,99003) IT,IOBSE,ICHR
C
         END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
         IF ( .NOT.GGA ) THEN
C
            DO IR = 1,IRTOP
C
C----------------------------- generate the densities on an angular mesh
C
               RHO4PI(1:N_LEBGRID,1:2) = 0.0D0
C
               FPIPR2 = CONST_4PI/R(IR,IM)**2
               DO ISPIN = 1,2
                  DO LM = 1,NLMFP
                     CALL DAXPY(N_LEBGRID,RHO2NS(IR,LM,IT,ISPIN)*FPIPR2,
     &                          Y_LEBGRID(1,LM),1,RHO4PI(1,ISPIN),1)
                  END DO
               END DO
C
               GAMMA_LEBO(:,:) = 0.0D0
C
               DO IOBSE = 1,NOBSE
                  DO LM = 1,NLMFP
                     CALL DAXPY(N_LEBGRID,GAMMA_O(IR,LM,IOBSE),
     &                          Y_LEBGRID(1,LM),1,GAMMA_LEBO(1,IOBSE),1)
                  END DO
               END DO
C
C----------- interchange UP and DOWN to be consistent with <EXCVWN> etc.
C
               RHO4PI(1:N_LEBGRID,2) = -RHO4PI(1:N_LEBGRID,2)
C
C-------------------------------------- calculate the ex.-cor. potential
C
               IF ( SCFVXC.EQ.'NONE      ' ) THEN
C
               ELSE IF ( SCFVXC(1:3).EQ.'VWN' ) THEN
C
                  CALL CHIVWN(RHO4PI(1,1),RHO4PI(1,2),EMM,ENM,ENN,IT,
     &                        N_LEBGRID)
C
               ELSE
C
                  WRITE (6,99001) SCFVXC
                  CALL STOP_MESSAGE(ROUTINE,
     &                              'parametrisation not available')
C
               END IF
C
C?????????????????????????????????????????? NORMALISAZION 1/4PI ?
               VXC_LEBOP(:,:,:) = 0D0
               VXC_LEBOP(1:NLG,2,2) = GAMMA_LEBO(1:NLG,2)*EMM(1:NLG)
               VXC_LEBOP(1:NLG,2,1) = GAMMA_LEBO(1:NLG,2)*ENM(1:NLG)
               VXC_LEBOP(1:NLG,1,2) = GAMMA_LEBO(1:NLG,1)*ENM(1:NLG)
               VXC_LEBOP(1:NLG,1,1) = GAMMA_LEBO(1:NLG,1)*ENN(1:NLG)
C
C---------------- expand the ex.-cor. potential into spherical harmonics
C
               DO IPERT = 1,NPERT
                  DO IOBSE = 1,NOBSE
                     DO LM = 1,NLMFP
C
                        VXC_LMTOP(IR,LM,IT,IOBSE,IPERT)
     &                     = DDOT(N_LEBGRID,VXC_LEBOP(1,IOBSE,IPERT),1,
     &                     W_LEBGRID(1,LM),1)
C
                     END DO
                  END DO
               END DO
C
            END DO
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
C
         ELSE
C
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
C
            DO IR = 1,IRTOP
               RPWM2(IR) = 1D0/R(IR,IM)**2
            END DO
C
C----------- interchange UP and DOWN to be consistent with <EXCVWN> etc.
C
C=======================================================================
C                           spin UP density
C=======================================================================
C
C----------------------------- density  FLM  in real spherical harmonics
C
            DO LM = 1,NLMFP
               DO IR = 1,IRTOP
                  X = 0.5D0*RPWM2(IR)
C                 FLM(IR,LM) = (RHO2NS(IR,LM,IT,1)+RHO2NS(IR,LM,IT,2))*X
                  FLM(IR,LM) = (RHO2NS(IR,LM,IT,1)-RHO2NS(IR,LM,IT,2))*X
               END DO
            END DO
C
            CALL FPVXCLM_DERIVE(IM,IRTOP,RHOU,AGRDRHOU,LAPRHOU,GDGAGU,
     &                          GFLM,GFANG,FLM,HLM,HANG,RPWM2,DFDR,
     &                          D2FDR2,NLMFP,Y_LEBGRID,W_LEBGRID,
     &                          N_LEBGRID,NMAX_LEBGRID)
C
C=======================================================================
C                           spin DOWN density
C=======================================================================
C
            DO LM = 1,NLMFP
               DO IR = 1,IRTOP
                  X = 0.5D0*RPWM2(IR)
C                 FLM(IR,LM) = (RHO2NS(IR,LM,IT,1)-RHO2NS(IR,LM,IT,2))*X
                  FLM(IR,LM) = (RHO2NS(IR,LM,IT,1)+RHO2NS(IR,LM,IT,2))*X
               END DO
            END DO
C
            CALL FPVXCLM_DERIVE(IM,IRTOP,RHOD,AGRDRHOD,LAPRHOD,GDGAGD,
     &                          GFLM,GFANG,FLM,HLM,HANG,RPWM2,DFDR,
     &                          D2FDR2,NLMFP,Y_LEBGRID,W_LEBGRID,
     &                          N_LEBGRID,NMAX_LEBGRID)
C
C=======================================================================
C                           total density
C=======================================================================
C
            DO LM = 1,NLMFP
               DO IR = 1,IRTOP
                  X = RPWM2(IR)
                  FLM(IR,LM) = RHO2NS(IR,LM,IT,1)*X
               END DO
            END DO
C
            CALL FPVXCLM_DERIVE(IM,IRTOP,RHO,AGRDRHO,LAPRHO,GDGAG,GFLM,
     &                          GFANG,FLM,HLM,HANG,RPWM2,DFDR,D2FDR2,
     &                          NLMFP,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &                          NMAX_LEBGRID)
C
C=======================================================================
C
            IF ( SCFVXC.EQ.'PBEXXXXXXX' ) THEN
C
               DO IR = 1,IRTOP
C
C                  CALL EXCPBE(VXC,EXCIJGGA(1,IR),WEXCIJ,N_LEBGRID,
C     &                        NMAX_LEBGRID,RHO(1,IR),AGRDRHO(1,IR),
C     &                        LAPRHO(1,IR),GDGAG(1,IR),RHOU(1,IR),
C     &                        AGRDRHOU(1,IR),LAPRHOU(1,IR),GDGAGU(1,IR),
C     &                        RHOD(1,IR),AGRDRHOD(1,IR),LAPRHOD(1,IR),
C     &                        GDGAGD(1,IR))
C
C---------------- expand the ex.-cor. potential into spherical harmonics
C
                  VXC_LMTOP(:,:,:,:,:) = 0D0
C
               END DO
C
            ELSE
C
               CALL STOP_MESSAGE(ROUTINE,
     &                           'GGA parametrisation not available')
C
            END IF
C
         END IF
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
C
C ----------------------------------------------------------------------
C             integrate  Delta Vxc x  GAMMA(SPN) x shape functions
C                    to get an estimate for I_Stoner
C ----------------------------------------------------------------------
C
         GAMMA_X_FLM(:,:) = 0D0
C
         GAMMA_X_FLM(1:IRMTIN,1:NLMFPMAX)
     &      = GAMMA_O(1:IRMTIN,1:NLMFPMAX,2)
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
                     GAMMA_X_FLM(IR,LM1) = GAMMA_X_FLM(IR,LM1)
     &                  + RWGT*GAMMA_O(IR,LM2,2)
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
         RINT1(:) = 0D0
         DO LM = 1,NLMFP
            DO IR = 1,IRCRIT
               RINT1(IR) = RINT1(IR) + VXC_LMTOP(IR,LM,IT,2,2)
     &                     *GAMMA_X_FLM(IR,LM)*R2DRDI(IR,IM)
            END DO
         END DO
C
         CALL RRADINT(IM,RINT1,I_STONER)
C
         WRITE (6,99002) IT,I_STONER*RY_EV
C
C ----------------------------------------------------------------------
C         copy perturbation term and convolute with shape functions
C             restrict to terms diagonal w.r.t. operator
C ----------------------------------------------------------------------
C
         DO IPERT = 1,NPERT
            IOBSE = IPERT
C
            HX_LM(:,:) = 0D0
C
C---------------------------- convolute HX_LM with shape functions FLMSF
C------------------------------ and weight R2DRDI for radial integration
C
C------------------------------------- perturbation HLM in the MT regime
C
            DO LMDV = 1,NLMFPT(IT)
               IF ( K_PERT_LMTP(LMDV,IT,IPERT) ) THEN
                  DO IR = 1,IRMTIN
                     HX_LM(IR,LMDV) = VXC_LMTOP(IR,LMDV,IT,IPERT,IPERT)
     &                                *R2DRDI(IR,IM)
                  END DO
               END IF
            END DO
C
C------------------------------------- perturbation HLM in the IS regime
C
            DO LMDV = 1,NLMFPT(IT)
               LDV = L_LM(LMDV)
               IF ( K_PERT_LMTP(LMDV,IT,IPERT) ) THEN
C
                  DO ISF = 1,NSF(IM)
                     LMSF = LMISF(ISF,IM)
                     LSF = L_LM(LMSF)
                     IF ( LSF.LE.(NLSF-1) ) THEN
C
                        DO LH = ABS(LSF-LDV),MIN((LSF+LDV),(NLFP-1)),2
                           LMH = LH*LH
                           DO MH = -LH,LH
                              LMH = LMH + 1
C
                              IF ( LMH.NE.(LH*LH+LH+MH+1) )
     &                             CALL STOP_MESSAGE(ROUTINE,'# LMH #')
C
                              G123 = GNTTAB(LMDV,LMH,LMSF)
C
                              IF ( ABS(G123).GT.1D-10 ) THEN
C
                                 K_PERT_LMTP(LMH,IT,IPERT) = .TRUE.
C
                                 DO IR = IRMTIN + 1,IRCRIT
                                    IRSF = IR - IRMTIN
C
                                    HX_LM(IR,LMH) = HX_LM(IR,LMH)
     &                                 + VXC_LMTOP(IR,LMDV,IT,IPERT,
     &                                 IPERT)*R2DRDI(IR,IM)
     &                                 *FLMSF(IRSF,ISF,IM)*G123
C
                                 END DO
C
                              END IF
C
                           END DO
                        END DO
C
                     END IF
                  END DO
C
               END IF
C
            END DO
C
            HX_PERT_LMTP(1:IRCRIT,1:NLMFPT(IT),IT,IPERT)
     &         = HX_LM(1:IRCRIT,1:NLMFPT(IT))
C
         END DO
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (//,10X,'SCFVXC = ',A,'  not known ',//)
99002 FORMAT (10X,' IT= ',I2,'   Stoner I',F13.7,' eV')
99003 FORMAT (/,' IT ',I3,2X,'OBS:',I3,2X,'INT',F20.10)
      END
