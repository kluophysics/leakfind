C*==thermal_properties.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE THERMAL_PROPERTIES(IT,IECURR,IEPATH,WE,TSST,MSST,MSSQ,
     &                              TAUQ,TMAT_LOC_AVG,MREG,MIRR,
     &                              DOBS_TX_GLO,OBS_TX_GLO,NOBS)
C   ********************************************************************
C   *                                                                  *
C   * subroutine to calculate thermal averages                         *
C   *                                                                  *
C   *  - the Weiss field within RDLM                                   *
C   *  - the magnetic torque                                           *
C   *  - observables                                                   *
C   *                                                                  *
C   *  with all vectorial quantities referring to the GLOBAL frame     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_CONSTANTS,ONLY:C0,C1,CI,PI,CPRE,KB_SI,EV_J,RY_EV
      USE MOD_SITES,ONLY:IQAT,NQMAX
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,ITRQ,IWFD,NOBSMAX,NKMMAX,NKMQ,
     &    NKM,IPIVKM,WKM1,WKM2,WKM3,WKM4,NPOL,TXT_OBS,AME_G
      USE MOD_TYPES,ONLY:NTMAX,TXT_T,W_WEISS_T,DOBS_TEX_GLO
      USE MOD_CALCMODE,ONLY:IREL,GF_CONV_RH,PROGNAME
      USE MOD_ENERGY,ONLY:NETAB
      USE MOD_THERMAL,ONLY:NFT,FMAT_FT,MFMAT_FT,NVIBFLU,IQ_AVFT,NVIBRA,
     &    NFLUCT,UMAT_VT,X_VFT,X_FT,EHAT_FLUCT,NHAT_MAG,W0_FLUCT,
     &    TEMP_LAT
      IMPLICIT NONE
C*--THERMAL_PROPERTIES29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_PROPERTIES')
      REAL*8 HPRE
      PARAMETER (HPRE=-(1D0/PI)*3.D0)
      INTEGER NOBS_THERMAL
      PARAMETER (NOBS_THERMAL=3)
C
C Dummy arguments
C
      INTEGER IECURR,IEPATH,IT,NOBS
      COMPLEX*16 WE
      COMPLEX*16 DOBS_TX_GLO(0:3,NOBSMAX,NTMAX),
     &           MIRR(NKMMAX,NKMMAX,3,NOBS),MREG(NKMMAX,NKMMAX,3,NOBS),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TMAT_LOC_AVG(NKMMAX,NKMMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
      REAL*8 OBS_TX_GLO(0:3,NOBSMAX,NTMAX)
C
C Local variables
C
      REAL*8 BETA_W_WEISS,EHAT_AVG(3),E_DOT_N,OBS_ABS(:),SUM_X_FT,
     &       SUM_X_VFT
      COMPLEX*16 CAUX,DMAMC(:,:),DTIL_VFT(:,:),JMAT(:,:,:),LOGDROT,
     &           MIRR_FT(:,:,:,:),MIRR_VFT(:,:,:,:),MREG_FT(:,:,:,:),
     &           MREG_VFT(:,:,:,:),MSS_VFT(:,:),MTIL_VFT(:,:),
     &           TAU_VFT(:,:),TMAT(:,:),TSS_FT(:,:),TSS_VFT(:,:)
      COMPLEX*16 CMATDET,CMATTRC
      REAL*8 DDOT,DNRM2
      INTEGER I,IA_ERR,IC,IFLUCT,IFT,IOBS,IPOL,IQ,IVFT,IVIBRA,IVT,M,N
      LOGICAL INITIALIZE
      SAVE JMAT
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE TMAT,OBS_ABS
      ALLOCATABLE MREG_FT,MIRR_FT,MREG_VFT,MIRR_VFT
      ALLOCATABLE TSS_FT,TSS_VFT,MSS_VFT,TAU_VFT,DTIL_VFT,MTIL_VFT
      ALLOCATABLE DMAMC,JMAT
C
      ALLOCATE (OBS_ABS(NOBSMAX))
      ALLOCATE (TMAT(NKMMAX,NKMMAX))
C
      ALLOCATE (MREG_FT(NKMMAX,NKMMAX,3,NOBS))
      ALLOCATE (MIRR_FT(NKMMAX,NKMMAX,3,NOBS))
      ALLOCATE (MREG_VFT(NKMMAX,NKMMAX,3,NOBS))
      ALLOCATE (MIRR_VFT(NKMMAX,NKMMAX,3,NOBS))
      ALLOCATE (TSS_FT(NKMMAX,NKMMAX),TSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (DTIL_VFT(NKMMAX,NKMMAX))
      ALLOCATE (MTIL_VFT(NKMMAX,NKMMAX))
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (TAU_VFT(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MIRR')
C
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C
         ALLOCATE (JMAT(NKMMAX,NKMMAX,3))
C
C-- calculate J(IPOL)  angular matrix elements of TOTAL angular momentum
C--------------------------------- spherical coordinates   (-), (0), (+)
C
         JMAT(1:NKM,1:NKM,:) = AME_G(1:NKM,1:NKM,:,IOMT)
     &                         + 0.5D0*AME_G(1:NKM,1:NKM,:,ISMT)
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
         CALL CMAT_CONVERT_POLAR(JMAT,'S>C')
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,'not adapted to FULLPOT')
      IF ( IREL.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'IREL <> 3 ')
      IF ( NOBS_THERMAL.GT.NOBSMAX )
     &      CALL STOP_MESSAGE(ROUTINE,'NOBS_THERMAL > NOBSMAX')
      IF ( NOBS_THERMAL.GT.NOBS )
     &      CALL STOP_MESSAGE(ROUTINE,'NOBS_THERMAL > NOBS')
      IF ( GF_CONV_RH )
     &      CALL STOP_MESSAGE(ROUTINE,'not adapted to RH convention')
C
      DOBS_TX_GLO(:,:,IT) = 0.0D0
C
      M = NKMMAX
      N = NKMQ(IQAT(1,IT))
C
      TMAT_LOC_AVG(:,:) = C0
C
C---------------------------------------------- <e_i> = sum(f) w_f * e_i
      IFT = (IT-1)*NFLUCT
      EHAT_AVG(:) = 0.D0
      DO IFLUCT = 1,NFLUCT
         IFT = IFT + 1
         EHAT_AVG(:) = EHAT_AVG(:) + EHAT_FLUCT(:,IFLUCT,IT)*X_FT(IFT)
      END DO
C
      BETA_W_WEISS = DABS(W_WEISS_T(IT))*(EV_J*RY_EV/KB_SI)/TEMP_LAT
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
      SUM_X_VFT = 0D0
      SUM_X_FT = 0D0
C
C============================================================= IFLUCT ==
C                                             account for local rotation
      IFT = (IT-1)*NFLUCT
      LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
         IFT = IFT + 1
         SUM_X_FT = SUM_X_FT + X_FT(IFT)
C
C--------------------------------------------- compensate local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
         IF ( NFLUCT.EQ.1 ) THEN
C
            TSS_FT(:,:) = TSST(:,:,IT)
C
            MREG_FT(:,:,:,:) = MREG(:,:,:,:)
            MIRR_FT(:,:,:,:) = MIRR(:,:,:,:)
C
         ELSE
C
            CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,TSST(1,1,IT)
     &                     ,'SAS+',TSS_FT)
C
            DO IOBS = 1,NOBS_THERMAL
C
               IF ( IOBS.EQ.IDOS ) THEN
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           MREG(1,1,3,IOBS),'SAS+',
     &                           MREG_FT(1,1,3,IOBS))
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           MIRR(1,1,3,IOBS),'SAS+',
     &                           MIRR_FT(1,1,3,IOBS))
C
               ELSE
C
                  CALL ME_ROTATE(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                           MREG(1,1,1,IOBS),MREG_FT(1,1,1,IOBS))
C
                  CALL ME_ROTATE(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                           MIRR(1,1,1,IOBS),MIRR_FT(1,1,1,IOBS))
C
               END IF
C
            END DO
C
         END IF
C
C============================================================= IVIBRA ==
C                                         account for local displacement
         IVT = (IT-1)*NVIBRA
         LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
            IVT = IVT + 1
C
            IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)*NFLUCT + IFLUCT
C
C----------------------------------------- compensate local displacement
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
            IF ( NVIBRA.EQ.1 ) THEN
C
               TSS_VFT(:,:) = TSS_FT(:,:)
C
               MREG_VFT(:,:,:,:) = MREG_FT(:,:,:,:)
               MIRR_VFT(:,:,:,:) = MIRR_FT(:,:,:,:)
C
            ELSE
C
               CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,'UAUT',TSS_VFT)
C
               DO IOBS = 1,NOBS_THERMAL
                  DO IPOL = 1,3
C
                     CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),
     &                                 MREG_FT(1,1,IPOL,IOBS),'UAUT',
     &                                 MREG_VFT(1,1,IPOL,IOBS))
C
                     CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),
     &                                 MIRR_FT(1,1,IPOL,IOBS),'UAUT',
     &                                 MIRR_VFT(1,1,IPOL,IOBS))
                  END DO
               END DO
C
            END IF
C
C-----------------------------------------------------------------------
C                 all matrices now in GLOBAL frame
C-----------------------------------------------------------------------
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
            IF ( NVIBFLU.EQ.1 ) THEN
               MSS_VFT(:,:) = MSST(:,:,IT)
            ELSE
               CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WKM1,MSS_VFT)
            END IF
C
C-------------------------------------- get projected TAU matrix TAU_vft
C
            IQ = IQ_AVFT(1,IVFT)
C
C----------------------------------------------- Delta m = m_vft - m_CPA
C
            DMAMC(1:N,1:N) = MSS_VFT(1:N,1:N) - MSSQ(1:N,1:N,IQ)
C
C--------------------------------------------- (m_vft - m_CPA) * TAU_CPA
C
            CALL CMATMUL(N,M,DMAMC,TAUQ(1,1,IQ),MTIL_VFT)
C
C---------------------------- M~_vft = ( 1 + (m_vft - m_CPA) * TAU_CPA )
C
            DO I = 1,N
               MTIL_VFT(I,I) = C1 + MTIL_VFT(I,I)
            END DO
C
C--------------------------------------------------- D~_vft = 1 / M~_vft
C
            CALL CMATINV3(N,M,IPIVKM,MTIL_VFT,WKM1,DTIL_VFT)
C
C------------------------------------------- TAU~_vft = TAU_CPA * D~_vft
C
            CALL CMATMUL(N,M,TAUQ(1,1,IQ),DTIL_VFT,TAU_VFT)
C
C=======================================================================
C                calculate the Weiss field and torque
C=======================================================================
C
C-----------------------------LOGDROT = log | 1 + [m(t) - m(c)] tau_c |
C
            CAUX = CMATDET(N,M,MTIL_VFT)
C
            LOGDROT = LOG(CAUX)
C
C-----------------------------------------------------------------------
C      RDLM Weiss field         S1 =  Int de e.n LOGDROT
C      PRB,74,144411,(2006) Eq. (34)
C-----------------------------------------------------------------------
C
            E_DOT_N = DDOT(3,EHAT_FLUCT(1,IFLUCT,IT),1,NHAT_MAG,1)
C
            CAUX = HPRE*WE*E_DOT_N*W0_FLUCT(IFLUCT)*LOGDROT
C
            OBS_TX_GLO(0,IWFD,IT) = OBS_TX_GLO(0,IWFD,IT) + DIMAG(CAUX)
C
C-----------------------------------------------------------------------
C       Torque         TORQUE =  Int de dP(e)/dn logdrot  for  T > 0 K
C       PRB,74,144411,(2006) Eq. (39)
C-----------------------------------------------------------------------
C
            IF ( BETA_W_WEISS.LT.150.0D0 ) THEN
C
               DO IC = 1,3
C
                  CAUX = -CPRE*WE*BETA_W_WEISS*(EHAT_FLUCT(IC,IFLUCT,IT)
     &                   -EHAT_AVG(IC))*X_FT(IFT)*LOGDROT
C
                  OBS_TX_GLO(IC,ITRQ,IT) = OBS_TX_GLO(IC,ITRQ,IT)
     &               + DIMAG(CAUX)
C
               END DO
C
            ELSE
C
C-----------------------------------------------------------------------
C         TORQUE = Trace  i TAU [(l+s)t^-1 - t^-1(l+s)]   for  T = O K
C-----------------------------------------------------------------------
C
               DO IC = 1,3
                  CALL CMATMUL(N,M,JMAT(1,1,IC),MSS_VFT,WKM2)
                  CALL CMATMUL(N,M,MSS_VFT,JMAT(1,1,IC),WKM3)
                  WKM4(1:N,1:N) = WKM2(1:N,1:N) - WKM3(1:N,1:N)
                  CALL CMATMUL(N,M,TAU_VFT,WKM4,WKM2)
C
                  CAUX = CPRE*CI*WE*CMATTRC(N,M,WKM2)
C
                  OBS_TX_GLO(IC,ITRQ,IT) = OBS_TX_GLO(IC,ITRQ,IT)
     &               + DIMAG(CAUX)
C
               END DO
C
            END IF
C=======================================================================
C            calculate magnetic moments in GLOBAL frame
C=======================================================================
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
            CALL TAUGFCONV(MSS_VFT,TAU_VFT,TMAT)
C
            DO IOBS = 1,NOBS_THERMAL
               DO IPOL = 1,NPOL
C
                  WKM1(1:N,1:N) = MATMUL(MREG_VFT(1:N,1:N,IPOL,IOBS),
     &                            TMAT(1:N,1:N))
C
                  WKM2(1:N,1:N) = WKM1(1:N,1:N)
     &                            - MIRR_VFT(1:N,1:N,IPOL,IOBS)
C
                  DOBS_TX_GLO(IPOL,IOBS,IT) = DOBS_TX_GLO(IPOL,IOBS,IT)
     &               - X_VFT(IVFT)*CMATTRC(N,M,WKM2)/PI
C
               END DO
            END DO
C
C=======================================================================
C            bring TAU_VFT to LOCAL frame and take the average
C=======================================================================
C
            IF ( NVIBRA.EQ.1 ) THEN
C
               WKM2(:,:) = TAU_VFT(:,:)
C
            ELSE
C
               CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TAU_VFT,'UTAU',WKM2)
C
            END IF
C
            IF ( NFLUCT.EQ.1 ) THEN
C
               TMAT_LOC_AVG(:,:) = TMAT_LOC_AVG(:,:) + X_VFT(IVFT)
     &                             *WKM2(:,:)
C
            ELSE
C
               CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,WKM2,
     &                        'S+AS',WKM3)
C
               TMAT_LOC_AVG(:,:) = TMAT_LOC_AVG(:,:) + X_VFT(IVFT)
     &                             *WKM3(:,:)
C
            END IF
C
C-----------------------------------------------------------------------
C
            SUM_X_VFT = SUM_X_VFT + X_VFT(IVFT)
C
         END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
      END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
C      IF ( ABS(SUM_X_FT-1D0).GT.1D-6 )
C     &      CALL STOP_MESSAGE(ROUTINE,'SUM_X_FT != 1')
C
      TMAT_LOC_AVG(:,:) = TMAT_LOC_AVG(:,:)/SUM_X_VFT
C
      DOBS_TX_GLO(1:NPOL,:,IT) = DOBS_TX_GLO(1:NPOL,:,IT)/SUM_X_VFT
C
      OBS_TX_GLO(1:NPOL,:,IT) = OBS_TX_GLO(1:NPOL,:,IT)
     &                          + DIMAG(WE*DOBS_TX_GLO(1:NPOL,:,IT))
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C         store average spin-projected DOS for later write out
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      IF ( PROGNAME(4:6).EQ.'GEN' ) DOBS_TEX_GLO(:,:,IT,IECURR)
     &     = DIMAG(DOBS_TX_GLO(:,:,IT))
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      IF ( IECURR.LE.NETAB(IEPATH) ) RETURN
C
C=======================================================================
C============================================ EF reached for path NEPATH
C=======================================================================
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C                     write magnetic moments
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
      WRITE (6,99001) IT,TXT_T(IT),OBS_TX_GLO(3,IDOS,IT)
C
      DO IOBS = 2,NOBS_THERMAL
         OBS_ABS(IOBS) = DNRM2(3,OBS_TX_GLO(1,IOBS,IT),1)
C
         WRITE (6,99002) TXT_OBS(IOBS),OBS_ABS(IOBS),
     &                   (OBS_TX_GLO(IPOL,IOBS,IT),IPOL=1,3)
      END DO
C
C=======================================================================
99001 FORMAT (10X,'IT = ',I4,4X,A,/,10X,'CHR',F10.6)
99002 FORMAT (10X,A,F10.6,2X,3F10.6)
C
      END
