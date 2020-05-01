C*==sig.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the conductivity tensor                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,CHIZ,
c modified by XJQ; not_store_chiz and use_jtauk_square, for big system, don't store huge matrix chiz
c                  prt_sig_realz, print sig(E) for kappa program using Kubo-Greenwood
     &    not_store_chiz, use_jtauk_square, 
     &    prt_sig_realz, ie0_sig, ie1_sig
c end-mod-xjq
      USE MOD_THERMAL,ONLY:NVT,X_VFT,NVFO_Q,IVFT_VFOQ,NVFTMAX,UMAT_VT,
     &    I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT
      USE MOD_CPA,ONLY:ITCPAMAX,NCPA
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,MEZJ,MEZZ,MSST,SSST,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI,IQBOT,IQTOP,ITOQ,
     &    NOQ
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,ABAS,ABAS_I,ALAT
      USE MOD_TYPES,ONLY:NT,NTMAX,ITBOT,ITTOP,CONC
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,CONSI2,
     &    LCHECKME_SPIN,NSPINPROJ,SIG_PREFAC_AU,LSP_NO_ALPHA,
     &    LSP_NO_NABLA,FLAG_ISPINPROJ,SIG_MODE,SOTSI,SOTSI2,
     &    SOT_PREFAC_AU,EESI,EESI2,EE_PREFAC_AU
      USE MOD_CALCMODE,ONLY:ORBPOL,SIGMA_PROJECT,THERMAL_VIBRA_FLUCT,
     &    DMFT,LDAU,KKRMODE,GF_CONV_RH,SOLVER_FP,LHS_SOL_EQ_RHS_SOL,
     &    PUBLIC_VERSION
      USE MOD_ENERGY,ONLY:ETAB,WETAB,NEPATH,NETAB,EFERMI,NEMAX,EIMAG,
     &    EMIN,EILOW,IGRID
      USE MOD_FILES,ONLY:IFILCBWF0,IFILCBWF_INCREMENT_ENERGY,
     &    IFILCBWF_INCREMENT_CC,RECLNGWF,RECLNGWF_SPH,IPRINT,WRTAUMQ,
     &    WRMAT,IOTMP,IOTMP2,FOUND_SECTION,FOUND_REAL,FOUND_INTEGER
c      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,MPI_ELOOP,MPI_ID,NPROCS
      use mod_mpi
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY_LINRESP
      USE MOD_CONSTANTS,ONLY:A0_SI,E0_SI,HBAR_SI,PI,CI,C0,RY_EV
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      IMPLICIT NONE
C*--SIG37
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG')
C
C Local variables
C
      LOGICAL CALCINT,CHANGE_FRAME,FOUND,GETIRRSOL,OPTICS,
     &        SIGMA_INTEGRATION,SINGLE_SITE_TEST,VERTEX
      REAL*8 COND_AU_TO_SI,CPACHNG,DELTA_E,DMFTMIX,EE_AU_TO_SI,
     &       EE_PREFAC_AU_STREDA,EMINTMP,EPS_OPT,ERYDTOP,GAM_OPT,
     &       MROT_FRAME(3,3),OMEGABOT,OMEGATAB(:),OMEGATOP,RHOAU_NV,
     &       RHOAU_VC,RHO_NV(3,3),RHO_VC(3,3),SIGMAAU_NV,SIGMAAU_VC,
     &       SIGMA_NV(3,3),SIGMA_VC(3,3),SIG_PREFAC_AU_STREDA,
     &       SOT_AU_TO_SI,SOT_PREFAC_AU_STREDA,TIME,TIME0,VUC_AU
      COMPLEX*16 CWORK(:),DERYD,ERYDA,ERYDACC,ERYDB,ERYDBCC,
     &           MAQAB(:,:,:,:,:),MAQBA(:,:,:,:,:),MAQBADUM(:,:,:,:,:),
     &           MBQAB(:,:,:,:,:),MBQBA(:,:,:,:,:),MCQAB(:,:,:,:,:),
     &           MCQBA(:,:,:,:,:),MDQAB(:,:,:,:,:),MDQBA(:,:,:,:,:),
     &           MDQBADUM(:,:,:,:,:),MEOFF_AHE_T(:,:,:,:,:),
     &           MEOFF_BARG_T(:,:,:,:,:),ME_SZ_T(:,:,:,:,:),
     &           MREG_TAB(:,:,:,:,:,:,:),MREG_TBA(:,:,:,:,:,:,:),
     &           MSSQA(:,:,:),MSSQB(:,:,:),MSSTA(:,:,:),MSSTB(:,:,:),PA,
     &           PACC,PB,PBCC,PHASK_DUMMY(:),S10AQAB(:,:,:,:),
     &           S10AQBA(:,:,:,:),S10BQAB(:,:,:,:),S10BQBA(:,:,:,:),
     &           S10CQAB(:,:,:,:),S10CQBA(:,:,:,:),S10DQAB(:,:,:,:),
     &           S10DQBA(:,:,:,:),S10_VFTAB(:,:,:,:,:,:),
     &           S10_VFTBA(:,:,:,:,:,:),S2AQAB(:,:,:,:),S2AQBA(:,:,:,:),
     &           S2BQAB(:,:,:,:),S2BQBA(:,:,:,:),S2CQAB(:,:,:,:),
     &           S2CQBA(:,:,:,:),S2DQAB(:,:,:,:),S2DQBA(:,:,:,:),
     &           S2_VFTAB(:,:,:,:,:,:),S2_VFTBA(:,:,:,:,:,:),
     &           S3AQAB(:,:,:,:),S3AQBA(:,:,:,:),S3BQAB(:,:,:,:),
     &           S3BQBA(:,:,:,:),S3CQAB(:,:,:,:),S3CQBA(:,:,:,:),
     &           S3DQAB(:,:,:,:),S3DQBA(:,:,:,:),S3_VFTAB(:,:,:,:,:,:),
     &           S3_VFTBA(:,:,:,:,:,:),S4AQAB(:,:,:,:),S4AQBA(:,:,:,:),
     &           S4BQAB(:,:,:,:),S4BQBA(:,:,:,:),S4CQAB(:,:,:,:),
     &           S4CQBA(:,:,:,:),S4DQAB(:,:,:,:),S4DQBA(:,:,:,:),
     &           S4_VFTAB(:,:,:,:,:,:),S4_VFTBA(:,:,:,:,:,:),
     &           SIG0Q(:,:,:,:),SIG0Q_FSEA(:,:,:,:),SIG0Q_OPT(:,:,:,:,:)
     &           ,SIG0Q_SURF(:,:,:,:),SIG0Q_SURF_PROJ(:,:,:,:),
     &           SIG1Q_FSEA_NV(:,:,:,:,:),SIG1Q_FSEA_VC(:,:,:,:,:),
     &           SIG1Q_NV(:,:,:,:,:),SIG1Q_OPT_NV(:,:,:,:,:,:),
     &           SIG1Q_OPT_VC(:,:,:,:,:,:)
      INTEGER IA_ERR,ICPACONV,ICPAFLAG,IDIMS,IE,IEA,IEB,IEPATH,
     &        IFIL_DUMP_TAU,IFIL_LHSA,IFIL_LHSA_CC,IFIL_LHSB,
     &        IFIL_LHSB_CC,IFIL_RHSA,IFIL_RHSA_CC,IFIL_RHSB,
     &        IFIL_RHSB_CC,IFIL_SPHA,IFIL_SPHA_CC,IFIL_SPHB,
     &        IFIL_SPHB_CC,ILOOP_SIGMA_PROJECT,IOMEGA,IPRINTL,IPROCE(:),
     &        IPROJ,IQ,ISP,ISPR,ITCPA,IWRIRRWF,IWRREGWF,M,NEA,NEB,
     &        NEINTEG,NOMEGA,NSP
      CHARACTER*10 INTEGER_TO_STRING
      COMPLEX*16 SIG1Q_SURF_NV(:,:,:,:,:),SIG1Q_SURF_NV_PROJ(:,:,:,:,:),
     &           SIG1Q_SURF_VC(:,:,:,:,:),SIG1Q_SURF_VC_PROJ(:,:,:,:,:),
     &           SIG1Q_VC(:,:,:,:,:),SIGOFFQ(:,:,:,:),TAUQA(:,:,:),
     &           TAUQB(:,:,:),TAUQZ(:,:,:,:),TSSQA(:,:,:),TSSQB(:,:,:),
     &           TSSTA(:,:,:),TSSTB(:,:,:),UMAT_VTA(:,:,:),
     &           UMAT_VTB(:,:,:),W_SIGMA_INTEGRATION
c modified by XJQ: prt_sig_realz
      logical lopen, lexist
      character*6 ctemp_lat
      integer ie_sig, errorcode, mpierr
      real*8 rho_diag_si(3)
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      DATA CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/,IFIL_DUMP_TAU/6/,VERTEX/.TRUE./
      DATA SINGLE_SITE_TEST/.FALSE./
C
      ALLOCATABLE OMEGATAB,SIG0Q_OPT,IPROCE,CWORK
      ALLOCATABLE SIG1Q_OPT_NV,SIG1Q_OPT_VC
      ALLOCATABLE S10AQAB,S2AQAB,S3AQAB,S4AQAB
      ALLOCATABLE S10BQAB,S2BQAB,S3BQAB,S4BQAB
      ALLOCATABLE S10CQAB,S2CQAB,S3CQAB,S4CQAB
      ALLOCATABLE S10DQAB,S2DQAB,S3DQAB,S4DQAB
      ALLOCATABLE S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB
      ALLOCATABLE S10AQBA,S2AQBA,S3AQBA,S4AQBA
      ALLOCATABLE S10BQBA,S2BQBA,S3BQBA,S4BQBA
      ALLOCATABLE S10CQBA,S2CQBA,S3CQBA,S4CQBA
      ALLOCATABLE S10DQBA,S2DQBA,S3DQBA,S4DQBA
      ALLOCATABLE S10_VFTBA,S2_VFTBA,S3_VFTBA,S4_VFTBA
C
      ALLOCATABLE TAUQZ
      ALLOCATABLE MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA
      ALLOCATABLE SIGOFFQ
      ALLOCATABLE TSSTA,TSSTB,UMAT_VTA,UMAT_VTB
      ALLOCATABLE MEOFF_AHE_T,MEOFF_BARG_T,ME_SZ_T
      ALLOCATABLE PHASK_DUMMY
      ALLOCATABLE TAUQA,TAUQB,MSSQA,MSSQB,TSSQA,TSSQB
      ALLOCATABLE MSSTA,MSSTB
      ALLOCATABLE MREG_TAB
      ALLOCATABLE MREG_TBA
      ALLOCATABLE SIG0Q,SIG0Q_FSEA,SIG0Q_SURF,SIG0Q_SURF_PROJ
      ALLOCATABLE SIG1Q_NV,SIG1Q_FSEA_NV,SIG1Q_SURF_NV
      ALLOCATABLE SIG1Q_VC,SIG1Q_FSEA_VC,SIG1Q_SURF_VC
      ALLOCATABLE SIG1Q_SURF_NV_PROJ,SIG1Q_SURF_VC_PROJ
      ALLOCATABLE MAQBADUM,MDQBADUM
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      ALLOCATE (TAUQA(M,M,NQMAX),TAUQB(M,M,NQMAX))
      ALLOCATE (MSSQA(M,M,NQMAX),MSSQB(M,M,NQMAX))
      ALLOCATE (TSSQA(M,M,NQMAX),TSSQB(M,M,NQMAX))
C
      ALLOCATE (SIG0Q_FSEA(3,3,NSPINPROJ,NQMAX))
C
      NSP = NSPINPROJ
C
      ALLOCATE (MAQAB(M,M,3,NSP,NQMAX),MAQBA(M,M,3,NSP,NQMAX))
      ALLOCATE (MBQAB(M,M,3,NSP,NQMAX),MBQBA(M,M,3,NSP,NQMAX))
      ALLOCATE (MCQAB(M,M,3,NSP,NQMAX),MCQBA(M,M,3,NSP,NQMAX))
      ALLOCATE (MDQAB(M,M,3,NSP,NQMAX),MDQBA(M,M,3,NSP,NQMAX))
      ALLOCATE (MAQBADUM(M,M,3,NSP,NQMAX),MDQBADUM(M,M,3,NSP,NQMAX))
C
      ALLOCATE (S10AQAB(3,3,NSP,NQMAX),S10AQBA(3,3,NSP,NQMAX))
      ALLOCATE (S10BQAB(3,3,NSP,NQMAX),S10BQBA(3,3,NSP,NQMAX))
      ALLOCATE (S10CQAB(3,3,NSP,NQMAX),S10CQBA(3,3,NSP,NQMAX))
      ALLOCATE (S10DQAB(3,3,NSP,NQMAX),S10DQBA(3,3,NSP,NQMAX))
      ALLOCATE (S2AQAB(3,3,NSP,NQMAX),S2AQBA(3,3,NSP,NQMAX))
      ALLOCATE (S2BQAB(3,3,NSP,NQMAX),S2BQBA(3,3,NSP,NQMAX))
      ALLOCATE (S2CQAB(3,3,NSP,NQMAX),S2CQBA(3,3,NSP,NQMAX))
      ALLOCATE (S2DQAB(3,3,NSP,NQMAX),S2DQBA(3,3,NSP,NQMAX))
      ALLOCATE (S3AQAB(3,3,NSP,NQMAX),S3AQBA(3,3,NSP,NQMAX))
      ALLOCATE (S3BQAB(3,3,NSP,NQMAX),S3BQBA(3,3,NSP,NQMAX))
      ALLOCATE (S3CQAB(3,3,NSP,NQMAX),S3CQBA(3,3,NSP,NQMAX))
      ALLOCATE (S3DQAB(3,3,NSP,NQMAX),S3DQBA(3,3,NSP,NQMAX))
      ALLOCATE (S4AQAB(3,3,NSP,NQMAX),S4AQBA(3,3,NSP,NQMAX))
      ALLOCATE (S4BQAB(3,3,NSP,NQMAX),S4BQBA(3,3,NSP,NQMAX))
      ALLOCATE (S4CQAB(3,3,NSP,NQMAX),S4CQBA(3,3,NSP,NQMAX))
      ALLOCATE (S4DQAB(3,3,NSP,NQMAX),S4DQBA(3,3,NSP,NQMAX))
C
      ALLOCATE (SIG0Q(3,3,NSPINPROJ,NQMAX))
      ALLOCATE (SIG0Q_SURF(3,3,NSPINPROJ,NQMAX),
     &          SIG1Q_SURF_NV(3,3,NSPINPROJ,NQMAX,NQMAX))
      ALLOCATE (SIG1Q_SURF_VC(3,3,NSPINPROJ,NQMAX,NQMAX))
      ALLOCATE (SIG1Q_NV(3,3,NSPINPROJ,NQMAX,NQMAX),
     &          SIG1Q_FSEA_NV(3,3,NSPINPROJ,NQMAX,NQMAX))
      ALLOCATE (SIG1Q_VC(3,3,NSPINPROJ,NQMAX,NQMAX),
     &          SIG1Q_FSEA_VC(3,3,NSPINPROJ,NQMAX,NQMAX))
      ALLOCATE (TAUQZ(NKMMAX,NKMMAX,NQMAX,2),STAT=IA_ERR)
      ALLOCATE (SIGOFFQ(3,3,NSPINPROJ,NQMAX))
      SIGOFFQ(:,:,:,:) = C0
C----------------------------------------------------------- dummy array
      ALLOCATE (PHASK_DUMMY(1))
C
C ======================================================================
      WRITE (6,99002)
C ======================================================================
      CALL CPU_TIME(TIME0)
C ======================================================================
C
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc MAQAB')
C
      VERTEX = .TRUE.
C
C ======================================================================
C   set the range of sites IQ and types IT to maximum
C   to allow for TB-calculations on 2D systems
C ======================================================================
C
      ITBOT = 1
      ITTOP = NT
      IQBOT = 1
      IQTOP = NQ
C
C ======================================================================
C     FP - mode: default settings are overwritten here
C                set back to ZJ convention
C
      IF ( GF_CONV_RH ) THEN
C
         GF_CONV_RH = .FALSE.
         SOLVER_FP = 'BS        '
C
         CALL INFO_MESSAGE(ROUTINE,'WF/GF convention set to ZJ')
C
      END IF
C
C ======================================================================
C
      IF ( NO_SYMMETRY_LINRESP ) WRITE (6,*) 'NO_SYMMETRY'
C
C=======================================================================
C=======================================================================
C                             INITIALIZE
C=======================================================================
C=======================================================================
C
C=======================================================================
C                        finite TEMPERATURE
C=======================================================================
C
      TEMP_LAT = 0D0
      N_TEMP_LAT = 1
C---------------------------------------- to get appropriate array sizes
      CALL THERMAL_INIT(0,N_TEMP_LAT,TEMP_LAT)
C
C=======================================================================
C
      M = NKMMAX
C
      ALLOCATE (TSSTA(M,M,NTMAX),TSSTB(M,M,NTMAX))
      ALLOCATE (MSSTA(M,M,NTMAX),MSSTB(M,M,NTMAX))
      ALLOCATE (UMAT_VTA(M,M,NVT),UMAT_VTB(M,M,NVT))
      ALLOCATE (MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2))
      ALLOCATE (MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2))
      ALLOCATE (MEOFF_AHE_T(NKMMAX,NKMMAX,3,3,NTMAX))
      ALLOCATE (MEOFF_BARG_T(NKMMAX,NKMMAX,3,3,NTMAX))
      ALLOCATE (ME_SZ_T(NKMMAX,NKMMAX,3,3,NTMAX),STAT=IA_ERR)
C
      NSP = NSPINPROJ
      ALLOCATE (S10_VFTAB(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S2_VFTAB(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S3_VFTAB(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S4_VFTAB(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S10_VFTBA(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S2_VFTBA(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S3_VFTBA(3,3,NSP,NVFTMAX,2,2))
      ALLOCATE (S4_VFTBA(3,3,NSP,NVFTMAX,2,2))
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( DMFT .OR. LDAU ) THEN
         CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),DMFTMIX,NEMAX)
         CALL CPU_TIME(TIME)
         if(mpi_id==0) write(*,*) 'init_mod_dmft_ldau takes ',
     &                            time-time0,' seconds'
         CALL CPU_TIME(TIME0)
      ENDIF
C
C ======================================================================
C ------------------------- conversion of resistivity result to SI-units
C ======================================================================
C
C-------------------------------------- calculate volume of atomic cells
C
      CALL CALC_VOL_AC
C
C ----------------------------------------------------------------------
C -------------------------------  get volume of unit cell VUC in  [m^3]
C -----------------------------  in case of layered systems: I-ZONE only
C
      IF ( SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV' .OR. 
     &     SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
         CALL RVECSPAT(ABAS_I(1,1),ABAS_I(1,2),ABAS_I(1,3),VUC_AU,1)
      ELSE
         CALL RVECSPAT(ABAS(1,1),ABAS(1,2),ABAS(1,3),VUC_AU,1)
      END IF
C
      VUC_AU = ABS(VUC_AU)*ALAT**3
C
C#######################################################################
C------ conversion factor: conductivity from a.u. (Rydberg!) to SI units
      COND_AU_TO_SI = ((E0_SI**2)/2D0)*(1D0/(HBAR_SI*A0_SI))
C
C------------------------------------ prefactor Bastin formula - in a.u.
C     SIG_PREFAC_AU = 1D0/(2D0*PI*VUC_AU)
C--- revised Nov 2015: normalisation according to volume of atomic cells
      SIG_PREFAC_AU = 1D0/(2D0*PI)
C
C-------- take care of fact that MEs to momentum operator   (m*c*\alpha)
C-------- and not to (c*\alpha) are calculated:
C-------- multiply by ( 1 / m0_au )**2 = ( 1 / 0.5 )**2 = 4
      SIG_PREFAC_AU = SIG_PREFAC_AU*4D0
C
C-------- NB: a factor 1/4 is already contained in sig0_surf, . . .
C---------    compensate for that
      SIG_PREFAC_AU = SIG_PREFAC_AU*4D0
C
C--------------------------------------------------- convert to SI units
      CONSI = COND_AU_TO_SI*SIG_PREFAC_AU
C
C----------------------------------------- prefactor Streda-term in a.u.
C-------------------------- NB: a factor 1/(2*i) is contained in
C--------------------------     sigoff where  ( G^(+) - G^(-) ) / (2i)
C--------------------------     is rewritten as Im (G^(+))
C     SIG_PREFAC_AU_STREDA = 1D0/(PI*VUC_AU)
C--- revised Nov 2015: normalisation according to volume of atomic cells
      SIG_PREFAC_AU_STREDA = 1D0/PI
C
CDK      CONSI2 = E0_SI**2*A0_SI**2/(2D0*HBAR_SI*PI*VUC_SI)
      CONSI2 = COND_AU_TO_SI*SIG_PREFAC_AU_STREDA
C
C#######################################################################
C
C#######################################################################
C------ conversion factor: torkance (SOT) from a.u. (Rydberg!) to SI units
C------ result in units of [10^-30 C * m]
      SOT_AU_TO_SI = E0_SI*A0_SI*1D30
C
C------------------------------------ prefactor Bastin formula - in a.u.
      SOT_PREFAC_AU = 1D0/(2D0*PI)
C
C-------- take care of fact that ME to momentum operator   (m*c*\alpha)
C-------- and not to (c*\alpha) is calculated for the perturbation:
C-------- multiply by 1 / m0_au =  1 / 0.5 = 2
      SOT_PREFAC_AU = SOT_PREFAC_AU*2D0
C
C-------- NB: a factor 1/4 is already contained in sig0_surf, . . .
C---------    compensate for that
      SOT_PREFAC_AU = SOT_PREFAC_AU*4D0
C
C--------------------------------------------------- convert to SI units
      SOTSI = SOT_AU_TO_SI*SOT_PREFAC_AU
C
C----------------------------------------- prefactor Streda-term in a.u.
C-------------------------- NB: a factor 1/(2*i) is contained in
C--------------------------     sigoff where  ( G^(+) - G^(-) ) / (2i)
C--------------------------     is rewritten as Im (G^(+))
      SOT_PREFAC_AU_STREDA = 1D0/PI
C
      SOTSI2 = SOT_AU_TO_SI*SOT_PREFAC_AU_STREDA
C
C#######################################################################
C
C#######################################################################
C------ conversion factor: Edelstein conductivity (EE) from a.u. (Rydberg!) to SI units
C------ result in units of [m/V = C/N = (A * s**3)/(kg * m)]
      EE_AU_TO_SI = A0_SI**2/E0_SI
C
C------------------------------------ prefactor Bastin formula - in a.u.
      EE_PREFAC_AU = 1D0/(2D0*PI)
C
C-------- take care of fact that ME to momentum operator   (m*c*\alpha)
C-------- and not to (c*\alpha) is calculated for the perturbation:
C-------- multiply by 1 / m0_au =  1 / 0.5 = 2
      EE_PREFAC_AU = EE_PREFAC_AU*2D0
C
C-------- NB: a factor 1/4 is already contained in sig0_surf, . . .
C---------    compensate for that
      EE_PREFAC_AU = EE_PREFAC_AU*4D0
C
C--------------------------------------------------- convert to SI units
      EESI = EE_AU_TO_SI*EE_PREFAC_AU
C
C----------------------------------------- prefactor Streda-term in a.u.
C-------------------------- NB: a factor 1/(2*i) is contained in
C--------------------------     sigoff where  ( G^(+) - G^(-) ) / (2i)
C--------------------------     is rewritten as Im (G^(+))
      EE_PREFAC_AU_STREDA = 1D0/PI
C
      EESI2 = EE_AU_TO_SI*EE_PREFAC_AU_STREDA
C
C#######################################################################
C
      ALLOCATE (LIST_ISPR(NSPINPROJ))
C
      NSPR = NSPINPROJ
      LIST_ISPR(1:11) = (/1,2,3,4,5,6,7,8,9,10,11/)
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_INTEGER_ARRAY('SIGSPROJ',LIST_ISPR,NSPR,
     &                                  NSPINPROJ,0,9999,0)
C
         CALL SECTION_FIND_KEYWORD('CHECKME',LCHECKME_SPIN)
         IF ( LCHECKME_SPIN ) WRITE (6,'(10X,"<SIG> checking MEs")')
C
C
C -------------- allow for switching off alpha part in spin-pol operator
         CALL SECTION_FIND_KEYWORD('SP_NO_ALPHA',FOUND)
         IF ( FOUND ) THEN
            LSP_NO_ALPHA = .TRUE.
            WRITE (6,99008) ROUTINE
         END IF
C
C -------------- allow for switching off nabla part in spin-pol operator
         CALL SECTION_FIND_KEYWORD('SP_NO_NABLA',FOUND)
         IF ( FOUND ) THEN
            LSP_NO_NABLA = .TRUE.
            WRITE (6,99009) ROUTINE
         END IF
C
C ----------------------------- run code in single-site mode for testing
         CALL SECTION_FIND_KEYWORD('SINGLESITE',SINGLE_SITE_TEST)
C
         IF ( SINGLE_SITE_TEST ) KKRMODE = 'SINGLE-SITE'
C
      END IF
C
      ALLOCATE (FLAG_ISPINPROJ(NSPINPROJ))
      FLAG_ISPINPROJ(:) = .FALSE.
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      IF ( PUBLIC_VERSION ) THEN
         NSPR = 1
         LIST_ISPR(1) = 1
         CALL INFO_MESSAGE(ROUTINE,
     &                 'PUBLIC VERSION gives only diagonal conductivity'
     &                 )
      END IF
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      WRITE (6,'(10x,"<SIG> number of chosen operators: ",i5,/)') NSPR
C
      DO ISPR = 1,NSPR
         WRITE (6,'(10X,"<SIG> operators: ",I3,1X,A80)') ISPR,
     &          STR_ISP_PROJ(LIST_ISPR(ISPR))
         FLAG_ISPINPROJ(LIST_ISPR(ISPR)) = .TRUE.
      END DO
C
C -------------- check whether current (alpha) - related term is present
C
      IF ( LIST_ISPR(1).NE.1 ) CALL STOP_REGULAR(ROUTINE,
     &     'current-operator alpha (1) has to be part of list')
C
C ======================================================================
C
      IF ( IBZINT.NE.2 .AND. IBZINT.NE.6 .AND. IBZINT.NE.4 )
     &     CALL STOP_MESSAGE(ROUTINE,'IBZINT   <>   2 | 4 | 6')
C
      IF ( IBZINT.EQ.4 ) THEN
C
         ERYDTOP = 1.4D0*EFERMI
C
         CALL FRELPOLE(ERYDTOP,IPRINT)
      END IF
C
      ITCPAMAX = 30
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      MPI_ELOOP = .FALSE.
      MPI_KLOOP = .FALSE.
C
      IF ( MPI ) MPI_KLOOP = .TRUE.
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C ======================================================================
C                  set SIG_MODE and parameters
C ======================================================================
C
      CALL INIT_MROT_FRAME(CHANGE_FRAME,MROT_FRAME)
C
C      SIG_MODE = 'INTEGRAL  '
      SIG_MODE = 'STREDA    '
C
      OPTICS = .FALSE.
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      EPS_OPT = 0.01D0
      GAM_OPT = 0.5D0/RY_EV
C
      NEA = 30
      NEB = 30
      NOMEGA = 100
      OMEGABOT = 0.1D0
      OMEGATOP = 15.0D0
      IGRID(1) = 9
      IGRID(2) = 9
C
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_FIND_KEYWORD('OPTICS',OPTICS)
C
         CALL SECTION_SET_INTEGER('NEA',NEA,9999,0)
         CALL SECTION_SET_INTEGER('NEB',NEB,9999,0)
         CALL SECTION_SET_INTEGER('IGRID',IGRID(1),9999,0)
         CALL SECTION_SET_INTEGER('NOMEGA',NOMEGA,9999,0)
         OPTICS = OPTICS .OR. FOUND_INTEGER
         CALL SECTION_SET_REAL('OMEGA1',OMEGABOT,9999D0,0)
         OPTICS = OPTICS .OR. FOUND_REAL
         CALL SECTION_SET_REAL('OMEGA2',OMEGATOP,9999D0,0)
         OPTICS = OPTICS .OR. FOUND_REAL
         CALL SECTION_SET_REAL('EPSOPT',EPS_OPT,9999D0,0)
         OPTICS = OPTICS .OR. FOUND_REAL
         CALL SECTION_SET_REAL('GAMOPT',GAM_OPT,9999D0,0)
         OPTICS = OPTICS .OR. FOUND_REAL
      END IF
C
      IF ( OPTICS ) SIG_MODE = 'OPTICS    '
C
      IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
         IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &        'PUBLIC VERSION does not supply optics')
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C-------------------------------------------------------------------- AC
         LIST_ISPR(1:11) = (/1,0,0,0,0,0,0,0,0,0,0/)
         NSPR = 1
C
         OMEGABOT = OMEGABOT/RY_EV
         OMEGATOP = OMEGATOP/RY_EV
         ALLOCATE (OMEGATAB(NOMEGA))
C
         NEA = MAX(3,NEA)
         NEB = MAX(3,NEB)
         NEMAX = MAX(NEA,NEB)
         NEPATH = 2
C
         IF ( ALLOCATED(ETAB) ) DEALLOCATE (ETAB,WETAB)
         ALLOCATE (ETAB(NEMAX,2),WETAB(NEMAX,2))
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C          GRID = 3 implies calculation of absorptive part only
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         IF ( IGRID(1).EQ.3 ) THEN
C
            IF ( ABS(EIMAG).LT.1D-8 ) EIMAG = 0.001D0
C
            IGRID(2) = 3
            SIG_MODE = 'OPTICS-ABS'
C
            DELTA_E = (EFERMI-EMIN)/DBLE(NEB)
C
            NEA = NINT(OMEGATOP/DELTA_E)
            NEMAX = MAX(NEA,NEB)
            DEALLOCATE (ETAB,WETAB)
            ALLOCATE (ETAB(NEMAX,2),WETAB(NEMAX,2))
C
            NOMEGA = NEA
C
C----------------------------------------------------- unoccupied states
C                                   straight path:  EF --- EF + OMEGATOP
C
            NETAB(1) = NEA
C
            DO IE = 1,NETAB(1)
               ETAB(IE,1) = EFERMI + (IE-0.5D0)*DELTA_E + CI*EIMAG
               WETAB(IE,1) = 1D0
            END DO
C
C------------------------------------------------------- occupied states
C                                            straight path:  EMIN --- EF
            NETAB(2) = NEB
C
            DO IE = 1,NETAB(2)
               ETAB(IE,2) = EMIN + (IE-0.5D0)*DELTA_E + CI*EIMAG
               WETAB(IE,2) = DELTA_E
            END DO
C
            OMEGABOT = DELTA_E
            OMEGATOP = DELTA_E*NEA
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         ELSE IF ( IGRID(1).EQ.9 ) THEN
C
            EIMAG = GAM_OPT*0.45D0
C
C----------------------------------------------------- unoccupied states
C                             elliptic path:  EF --- EF + OMEGATOP + 0.2
C
            EMINTMP = EFERMI - (OMEGATOP+0.2D0)
            NETAB(1) = NEA
            CALL EPATH(9,EMINTMP,EFERMI,EIMAG,NETAB(1),ETAB(1,1),
     &                 WETAB(1,1),EILOW,IPRINT,NEMAX)
C
            DO IE = 1,NETAB(1)
               ETAB(IE,1) = ETAB(IE,1) + 2D0*DREAL((EFERMI-ETAB(IE,1)))
            END DO
C
C------------------------------------------------------- occupied states
C                                            elliptic path:  EMIN --- EF
            NETAB(2) = NEB
            CALL EPATH(9,EMIN,EFERMI,EIMAG,NETAB(2),ETAB(1,2),WETAB(1,2)
     &                 ,EILOW,IPRINT,NEMAX)
C
         ELSE
            CALL STOP_MESSAGE(ROUTINE,
     &                        'IGRID='//INTEGER_TO_STRING(IGRID(1))
     &                        //' <> 3, 9')
         END IF
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         ALLOCATE (SIG0Q_OPT(3,3,NSPINPROJ,NQMAX,NOMEGA))
         ALLOCATE (SIG1Q_OPT_NV(3,3,NSPINPROJ,NQMAX,NQMAX,NOMEGA))
         ALLOCATE (SIG1Q_OPT_VC(3,3,NSPINPROJ,NQMAX,NQMAX,NOMEGA))
         SIG0Q_OPT(:,:,:,:,:) = C0
         SIG1Q_OPT_NV(:,:,:,:,:,:) = C0
         SIG1Q_OPT_VC(:,:,:,:,:,:) = C0
         DO IOMEGA = 1,NOMEGA
            OMEGATAB(IOMEGA) = OMEGABOT + (OMEGATOP-OMEGABOT)
     &                         *DBLE(IOMEGA-1)/DBLE(NOMEGA-1)
         END DO
C
         SIGMA_PROJECT = .FALSE.
C
         ALLOCATE (IPROCE(NEA))
         IPROCE(1:NEA) = 0
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) THEN
            MPI_ELOOP = .TRUE.
            MPI_KLOOP = .FALSE.
            CALL MPI_DISTRIBUTE(IPROCE,NEA,MPI_ELOOP,'E')
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
         WRITE (6,99010) NEB,DREAL(ETAB(1,2)),DREAL(ETAB(NEB,2)),NEA,
     &                   DREAL(ETAB(1,1)),DREAL(ETAB(NEA,1)),NOMEGA,
     &                   OMEGABOT,OMEGATOP,OMEGABOT*RY_EV,
     &                   OMEGATOP*RY_EV,EIMAG,EIMAG*RY_EV,GAM_OPT,
     &                   GAM_OPT*RY_EV,EPS_OPT,EPS_OPT*RY_EV,NSPR,
     &                   LIST_ISPR(1:NSPR)
         WRITE (6,99011) MPI,MPI_ELOOP,MPI_KLOOP,SINGLE_SITE_TEST
C
C-------------------------------------------------------------------- AC
C
      ELSE
C
C-------------------------------------------------------------------- DC
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
         IF ( PUBLIC_VERSION ) NETAB(1) = 1
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
         IF ( NETAB(1).GT.1 .and. .not.prt_sig_realz ) THEN
            NEINTEG = NETAB(1)
         ELSE
            NEINTEG = 0
         END IF
         NEA = NEINTEG + 1
         NEB = 0
C
         ALLOCATE (IPROCE(NEA))
         IPROCE(1:NEA) = 0
C
C        DERYD = 0.001D0
C--------------------------- new value, checked by KC on various systems
         DERYD = 0.0001D0
C-------------------------------------------------------------------- DC
      END IF
C
C-----------------------------------------------------------------------
C                 see MOD_FILES for the setting scheme
C-----------------------------------------------------------------------
C
      IFIL_RHSA = IFILCBWF0 + 1
      IFIL_RHSB = IFIL_RHSA + IFILCBWF_INCREMENT_ENERGY
C
      IFIL_RHSA_CC = IFIL_RHSA + IFILCBWF_INCREMENT_CC
      IFIL_RHSB_CC = IFIL_RHSB + IFILCBWF_INCREMENT_CC
C
      OPEN (IFIL_RHSA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      OPEN (IFIL_RHSB,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      OPEN (IFIL_RHSA_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
      OPEN (IFIL_RHSB_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      IF ( FULLPOT ) THEN
C
         CALL SET_IFIL_SPH(IFIL_RHSA,IFIL_SPHA)
C
         CALL SET_IFIL_SPH(IFIL_RHSB,IFIL_SPHB)
C
         CALL SET_IFIL_SPH(IFIL_RHSA_CC,IFIL_SPHA_CC)
C
         CALL SET_IFIL_SPH(IFIL_RHSB_CC,IFIL_SPHB_CC)
C
         OPEN (IFIL_SPHA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
         OPEN (IFIL_SPHB,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
         OPEN (IFIL_SPHA_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
         OPEN (IFIL_SPHB_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
      END IF
C
C
      IF ( .NOT.LHS_SOL_EQ_RHS_SOL ) THEN
C
         CALL SET_IFIL_LHS(IFIL_RHSA,IFIL_LHSA)
C
         CALL SET_IFIL_LHS(IFIL_RHSB,IFIL_LHSB)
C
         CALL SET_IFIL_LHS(IFIL_RHSA_CC,IFIL_LHSA_CC)
C
         CALL SET_IFIL_LHS(IFIL_RHSB_CC,IFIL_LHSB_CC)
C
         OPEN (IFIL_LHSA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
         OPEN (IFIL_LHSB,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
         OPEN (IFIL_LHSA_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
         OPEN (IFIL_LHSB_CC,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
      END IF
C
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C              loop to get the SPIN-projected conductivities       START
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      ILOOP_SIGMA_PROJECT = 0
C
 100  CONTINUE
      IF ( SIGMA_PROJECT ) THEN
C
         IF ( ILOOP_SIGMA_PROJECT.EQ.0 ) THEN
            ALLOCATE (SIG0Q_SURF_PROJ(3,3,NQMAX,3))
            ALLOCATE (SIG1Q_SURF_NV_PROJ(3,3,NQMAX,NQMAX,3))
            ALLOCATE (SIG1Q_SURF_VC_PROJ(3,3,NQMAX,NQMAX,3))
         END IF
C
         ILOOP_SIGMA_PROJECT = ILOOP_SIGMA_PROJECT + 1
C
      END IF
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C***********************************************************************
C***********************************************************************
C****************   START OF SIGMA CALCULATION   ***********************
C***********************************************************************
C***********************************************************************
C
C=======================================================================
C===============         LOOP OVER TEMPERATURE          ========== START
C=======================================================================
C
      LOOP_TEMP_LAT:DO I_TEMP_LAT = 1,N_TEMP_LAT
C
C ======================================================================
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C            and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ
C
         CALL THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
C=======================================================================
C
         if(mpi_id==0 .and. prt_sig_realz) then
           ctemp_lat=''
           if(temp_lat>0) write(ctemp_lat,'(i6)') int(temp_lat)
           inquire(iotmp2,opened=lopen)
           if(lopen) then
             write(*,*) 'iotmp2 is already opened'
             call mpi_abort(mpi_comm_world,errorcode,mpierr)
           endif
           if(temp_lat>0) then
             open(iotmp2,file='sig_realz_t'//trim(adjustl(ctemp_lat))//
     &'.dat',action='write',status='replace')
           else
             open(iotmp2,file='sig_realz.dat',
     &         action='write',status='replace')
           endif
           write(iotmp2,'(a,i4,a,i4)') 'ie0= ',ie0_sig,'ie1= ',ie1_sig
           write(iotmp2,*) 'energy, rho_xx(yy/zz), rho_avg, inv_cond'
         endif
         do ie_sig=ie0_sig,ie1_sig
C=======================================================================
         SIG0Q(:,:,:,:) = C0
         SIG0Q_FSEA(:,:,:,:) = C0
         SIG0Q_SURF(:,:,:,:) = C0
C
         SIG1Q_NV(:,:,:,:,:) = C0
         SIG1Q_FSEA_NV(:,:,:,:,:) = C0
         SIG1Q_SURF_NV(:,:,:,:,:) = C0
C
         SIG1Q_VC(:,:,:,:,:) = C0
         SIG1Q_FSEA_VC(:,:,:,:,:) = C0
         SIG1Q_SURF_VC(:,:,:,:,:) = C0
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop  START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         WRITE (6,99006)
C
         WRTAUMQ = .FALSE.
C
         NEPATH = 1
         IEPATH = 1
C
         ICPAFLAG = 0
         CPACHNG = 0.0D0
C
         ERYDA_EQ_ERYDB = .TRUE.
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
         IE = 0
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
         LOOP_IEA:DO IEA = 1,NEA
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ELOOP .AND. MPI_ID.NE.IPROCE(IEA) ) CYCLE LOOP_IEA
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IE = IE + 1
C
C-------------------------------------------------------------------- AC
            IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
C-------------------------------------------------------------------- AC
               SIGMA_INTEGRATION = .TRUE.
               ERYDA = ETAB(IEA,1)
C-------------------------------------------------------------------- DC
            ELSE IF ( IEA.LT.NEA ) THEN
               SIGMA_INTEGRATION = .TRUE.
               ERYDA = ETAB(IEA,IEPATH)
               NEB = 2
C
            ELSE
               SIGMA_INTEGRATION = .FALSE.
               if(prt_sig_realz) then
                 eryda = etab(ie_sig,1)
               else
                 ERYDA = EFERMI + CI*EIMAG
               endif
               NEB = 1
            END IF
C-------------------------------------------------------------------- DC
C
            IF ( LCHECKME_SPIN ) ERYDA = (0.6D0,0.02D0)
C
            IF ( IPRINT.GE.1 ) WRITE (6,99001) EFERMI,IEA,ERYDA
C
C--------------------------------------------- adjust number of k-points
            CALL SIG_DRV_KMESH(IEA,NEA,SIG_MODE)
C
            ERYDA_EQ_ERYDB = .NOT.SIGMA_INTEGRATION
C ===================================== solve SS - differential equation
C
            IF ( DMFT .OR. LDAU ) then
              if(prt_sig_realz) then
                DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &            = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,ie_sig)
              else
                DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &            = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
              endif
            endif
C
C------------------------------------------------------------------ Z(-)
C
            ERYDACC = DCONJG(ERYDA)
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSA_CC,
     &                    GETIRRSOL,ERYDACC,PACC,IPRINT,TSST,MSST,SSST,
     &                    MEZZ,MEZJ,ORBPOL)
C
C------------------------------------------------------------------ Z(+)
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSA,GETIRRSOL,
     &                    ERYDA,PA,IPRINT,TSSTA,MSSTA,SSST,MEZZ,MEZJ,
     &                    ORBPOL)
C
C ======================================================================
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
CDK## debug UMAT Begin
            IF ( 1.EQ.0 .AND. THERMAL_VIBRA_FLUCT )
     &           CALL CHECK_UMAT(ERYDA)
CDK## debug UMAT End
C
            CALL THERMAL_INIT_UFMAT(ERYDA)
C
            UMAT_VTA(:,:,:) = UMAT_VT(:,:,:)
C
C ======================================================================
C
            CALL MSSINIT(TSSTA,MSSTA,TSSQA,MSSQA)
C
C***********************************************************************
            IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
               CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYDA,PA,IPRINT,
     &                              ITCPA,ICPACONV,CONC,NOQ,ITOQ,
     &                              PHASK_DUMMY,1,NTMAX,TSSTA,MSSTA,
     &                              TSSQA,MSSQA,TAUQA)
C
C***********************************************************************
            ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
               CALL TAU_TB(ERYDA,PA,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,
     &                     TSSTA,MSSTA,TSSQA,MSSQA,TAUQA)
C
C***********************************************************************
            ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
               TAUQA(:,:,:) = TSSQA(:,:,:)
C
C***********************************************************************
            ELSE
C
               CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
            END IF
C***********************************************************************
C
C ??????????????????????????????????????????????????????????????????????
            IF ( IPRINT.GE.5 .OR. WRMAT ) THEN
C
               CALL PROJTAU(ICPAFLAG,CPACHNG,ERYDA,MSST,MSSQA,TAUQA,
     &                      TAUT)
C
               CALL DUMPTAU(IE,ERYDA,IFIL_DUMP_TAU,MSST,MSSQA,TAUT,
     &                      TAUQA)
C
            END IF
C ??????????????????????????????????????????????????????????????????????
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
            LOOP_IEB:DO IEB = 1,NEB
C
C-------------------------------------------------------------------- AC
               IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
                  ERYDB = ETAB(IEB,2)
C-------------------------------------------------------------------- AC
C-------------------------------------------------------------------- DC
               ELSE IF ( SIGMA_INTEGRATION ) THEN
                  IE = IE + 1
                  ERYDB = ERYDA + (IEB-1.5D0)*DERYD
               ELSE
c                  ERYDB = EFERMI + CI*EIMAG
                  erydb = eryda
               END IF
C-------------------------------------------------------------------- AC
C
               IF ( LCHECKME_SPIN ) ERYDB = (0.8D0,0.01D0)
C
               WRITE (6,99005) EFERMI,IEB,ERYDB,IEA,ERYDA,IE
C
C ======================================================================
C ===================================== solve SS - differential equation
C
               IF ( DMFT .OR. LDAU ) then
                 if(prt_sig_realz) then
                   DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &               = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,ie_sig)
                 else
                   DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &               = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
                 endif
               endif
C
C------------------------------------------------------------------ Z(-)
C
               ERYDBCC = DCONJG(ERYDB)
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSB_CC,
     &                       GETIRRSOL,ERYDBCC,PBCC,IPRINT,TSST,MSST,
     &                       SSST,MEZZ,MEZJ,ORBPOL)
C
C------------------------------------------------------------------ Z(+)
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSB,
     &                       GETIRRSOL,ERYDB,PB,IPRINT,TSSTB,MSSTB,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
C ======================================================================
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
               CALL THERMAL_INIT_UFMAT(ERYDB)
C
               UMAT_VTB(:,:,:) = UMAT_VT(:,:,:)
C
C ======================================================================
C
               CALL MSSINIT(TSSTB,MSSTB,TSSQB,MSSQB)
C
C ======================================================================
C
C
               IF ( .NOT.SIGMA_INTEGRATION ) THEN
C
                  DO IQ = IQBOT_CHI,IQTOP_CHI
                     TAUQB(:,:,IQ) = TAUQA(:,:,IQ)
                     MSSQB(:,:,IQ) = MSSQA(:,:,IQ)
                  END DO
C
C***********************************************************************
               ELSE IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
                  CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYDB,PB,IPRINT,
     &                                 ITCPA,ICPACONV,CONC,NOQ,ITOQ,
     &                                 PHASK_DUMMY,1,NTMAX,TSSTB,MSSTB,
     &                                 TSSQB,MSSQB,TAUQB)
C
C***********************************************************************
               ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
                  CALL TAU_TB(ERYDB,PB,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,
     &                        TSSTB,MSSTB,TSSQB,MSSQB,TAUQB)
C
C***********************************************************************
               ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
                  TAUQB(:,:,:) = TSSQB(:,:,:)
C
C***********************************************************************
               ELSE
C
                  CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
C***********************************************************************
C
               END IF
C
C ======================================================================
C
C***********************************************************************
               IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
                  if( not_store_chiz .or. use_jtauk_square ) then
                    if( SIGMA_INTEGRATION )
     &                stop 'SIGMA_INTEGRATION when not_store_chiz 
     &                      or use_jtauk_square'
                    if( SIG_MODE(1:6).EQ.'OPTICS' )
     &                stop 'OPTICS when not_store_chiz or
     &                      use_jtauk_square'
                    if( .not.NO_SYMMETRY_LINRESP )
     &                stop 'SYMMETRY_LINRESP when not store_chiz or
     &                      use_jtauk_square'
                    tauqz(:,:,:,1) = tauqb(:,:,:)
                    do iq = iqbot_chi, iqtop_chi
                      tauqz(:,:,iq,2) = transpose(dconjg(tauqb(:,:,iq)))
                    enddo
c                  IF ( IBZINT.EQ.4 ) THEN
                  ELSEIF ( IBZINT.EQ.4 ) THEN
c end-mod-xjq
C
                     CALL STOP_MESSAGE(ROUTINE,
     &                   'only <SIGKLOOPS> available for BZ integration'
     &                   )
C
                  ELSE IF ( NO_SYMMETRY_LINRESP ) THEN
C
                     CALL SIGKLOOPS_NOSYM(ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,
     &                  TAUQZ,MSSQA,MSSQB)
C
                  ELSE
C
                     CALL SIGKLOOPS(ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQZ,
     &                              MSSQA,MSSQB)
C
                  END IF
C
C***********************************************************************
               ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
                  CALL SIGKLOOPSDRV_TB(ERYDA,PA,MSSQA,TAUQZ)
C
C***********************************************************************
               ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
                  CHIZ(:,:,:) = C0
C
C***********************************************************************
               ELSE
C
                  CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
               END IF
C***********************************************************************
C
C ??????????????????????????????????????????????????????????????????????
               IF ( WRMAT ) THEN
C
                  OPEN (IOTMP,FILE='CHIZ_sig.dump')
                  WRITE (IOTMP,'(2D25.12)') CHIZ
C
                  CALL STOP_REGULAR(ROUTINE,
     &                              'CHIZ dumped to CHIZ_sig.dump')
C
               END IF
C ??????????????????????????????????????????????????????????????????????
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
c               IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
               IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP .or. not_store_chiz .or.
     &              use_jtauk_square ) THEN
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C-----------------------------------------------------------------------
C             matrix elements for energies    E_a   and    E_b
C-----------------------------------------------------------------------
C
                  CALL SIGME(ILOOP_SIGMA_PROJECT,ERYDA,ERYDACC,ERYDB,
     &                       ERYDBCC,IFIL_RHSA,IFIL_RHSA_CC,IFIL_RHSB,
     &                       IFIL_RHSB_CC,MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,
     &                       MCQBA,MDQAB,MDQBA,MREG_TAB,MREG_TBA,
     &                       S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,
     &                       S10AQAB,S10BQAB,S10CQAB,S10DQAB,S2AQAB,
     &                       S2BQAB,S2CQAB,S2DQAB,S3AQAB,S3BQAB,S3CQAB,
     &                       S3DQAB,S4AQAB,S4BQAB,S4CQAB,S4DQAB,
     &                       UMAT_VTA,MSSTA,TSSTA,MSSQA,TAUQA,UMAT_VTB,
     &                       MSSTB,TSSTB,MSSQB,TAUQB,NSPINPROJ,
     &                       MEOFF_BARG_T,MEOFF_AHE_T,ME_SZ_T)
C
C ======================================================================
C                          Fermi sea integration
C ======================================================================
                  IF ( SIGMA_INTEGRATION ) THEN
C
C-------------------------------------------------------------------- AC
                     IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
                        W_SIGMA_INTEGRATION = WETAB(IEA,1)*WETAB(IEB,2)
C-------------------------------------------------------------------- AC
C-------------------------------------------------------------------- DC
                     ELSE IF ( IEB.EQ.1 ) THEN
                        W_SIGMA_INTEGRATION = -WETAB(IEA,1)/DERYD
                     ELSE
                        W_SIGMA_INTEGRATION = +WETAB(IEA,1)/DERYD
C-------------------------------------------------------------------- DC
                     END IF
C
C-----------------------------------------------------------------------
C             matrix elements for energies    E_b   and    E_a
C-----------------------------------------------------------------------
C NB: to avoid overwriting regular MEs MAQBA,MDQBA (needed in SIG1_FSEA)
C     dummy arrays are introduced
C
                     CALL SIGME(ILOOP_SIGMA_PROJECT,ERYDB,ERYDBCC,ERYDA,
     &                          ERYDACC,IFIL_RHSB,IFIL_RHSB_CC,
     &                          IFIL_RHSA,IFIL_RHSA_CC,MAQBADUM,MAQAB,
     &                          MBQBA,MBQAB,MCQBA,MCQAB,MDQBADUM,MDQAB,
     &                          MREG_TBA,MREG_TAB,S10_VFTBA,S2_VFTBA,
     &                          S3_VFTBA,S4_VFTBA,S10AQBA,S10BQBA,
     &                          S10CQBA,S10DQBA,S2AQBA,S2BQBA,S2CQBA,
     &                          S2DQBA,S3AQBA,S3BQBA,S3CQBA,S3DQBA,
     &                          S4AQBA,S4BQBA,S4CQBA,S4DQBA,UMAT_VTB,
     &                          MSSTB,TSSTB,MSSQB,TAUQB,UMAT_VTA,MSSTA,
     &                          TSSTA,MSSQA,TAUQA,NSPINPROJ,
     &                          MEOFF_BARG_T,MEOFF_AHE_T,ME_SZ_T)
C
                     IPRINTL = IPRINT
                     IPRINT = 5
C-----------------------------------------------------------------------
C             SIG 0   --   site diagonal contributions
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------------- AC
                     IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
C-------------------------------------------------------------------- AC
                        CALL SIG0_OPTICS(S10AQAB,S10BQAB,S10CQAB,
     &                     S10DQAB,S2AQAB,S2BQAB,S2CQAB,S2DQAB,S3AQAB,
     &                     S3BQAB,S3CQAB,S3DQAB,S4AQAB,S4BQAB,S4CQAB,
     &                     S4DQAB,S10AQBA,S10BQBA,S10CQBA,S10DQBA,
     &                     S2AQBA,S2BQBA,S2CQBA,S2DQBA,S3AQBA,S3BQBA,
     &                     S3CQBA,S3DQBA,S4AQBA,S4BQBA,S4CQBA,S4DQBA,
     &                     SIG0Q_OPT,NSPINPROJ,NOMEGA,OMEGATAB,ERYDA,
     &                     ERYDB,EPS_OPT,GAM_OPT,IEA,IEB)
C-------------------------------------------------------------------- DC
                     ELSE
                        CALL SIG0_FSEA(W_SIGMA_INTEGRATION,S10AQAB,
     &                                 S10DQAB,S2AQAB,S2DQAB,S3AQAB,
     &                                 S3DQAB,S4AQAB,S4DQAB,S10AQBA,
     &                                 S10DQBA,S2AQBA,S2DQBA,S3AQBA,
     &                                 S3DQBA,S4AQBA,S4DQBA,NSPINPROJ,
     &                                 SIG0Q_FSEA)
C-------------------------------------------------------------------- DC
                     END IF
C
C-----------------------------------------------------------------------
C             SIG 1   --   without vertex corrections
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------------- AC
                     IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
                        CALL SIG1_OPTICS(MAQAB,MBQAB,MCQAB,MDQAB,MAQBA,
     &                     MCQBA,MBQBA,MDQBA,NSPINPROJ,NOMEGA,OMEGATAB,
     &                     ERYDA,ERYDB,EPS_OPT,GAM_OPT,IEA,IEB,
     &                     SIG1Q_OPT_NV,'N')
C-------------------------------------------------------------------- AC
C-------------------------------------------------------------------- DC
                     ELSE
                        CALL SIG1_FSEA(W_SIGMA_INTEGRATION,MAQAB,MAQBA,
     &                                 MDQAB,MDQBA,SIG1Q_FSEA_NV,'N',
     &                                 NSPINPROJ)
C-------------------------------------------------------------------- DC
                     END IF
C
C-----------------------------------------------------------------------
C                        vertex corrections
C-----------------------------------------------------------------------
C
                     IF ( VERTEX ) THEN
C
                        CALL LINRESP_VERTEX(TAUQA,TAUQB,TSSTA,TSSTB,
     &                     UMAT_VTA,UMAT_VTB,MSSQA,MSSQB)
C
C-----------------------------------------------------------------------
C             SIG 1   --   with vertex corrections
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------------- AC
                        IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
                           CALL SIG1_OPTICS(MAQAB,MBQAB,MCQAB,MDQAB,
     &                        MAQBA,MCQBA,MBQBA,MDQBA,NSPINPROJ,NOMEGA,
     &                        OMEGATAB,ERYDA,ERYDB,EPS_OPT,GAM_OPT,IEA,
     &                        IEB,SIG1Q_OPT_VC,'V')
C-------------------------------------------------------------------- AC
C-------------------------------------------------------------------- DC
                        ELSE
                           CALL SIG1_FSEA(W_SIGMA_INTEGRATION,MAQAB,
     &                        MAQBA,MDQAB,MDQBA,SIG1Q_FSEA_VC,'V',
     &                        NSPINPROJ)
C-------------------------------------------------------------------- DC
                        END IF
C
C
                     ELSE
C
                        SIG1Q_FSEA_VC(:,:,:,:,:)
     &                     = SIG1Q_FSEA_NV(:,:,:,:,:)
C
                     END IF
C
                  ELSE
C
C ======================================================================
C                           Fermi surface
C ======================================================================
C
                     IPRINTL = IPRINT
                     IPRINT = 5
C
C-----------------------------------------------------------------------
C             SIG 0   --   site diagonal contributions
C-----------------------------------------------------------------------
C
                     CALL SIG0_SURF(S10AQAB,S10BQAB,S10DQAB,S2AQAB,
     &                              S2BQAB,S2DQAB,S3AQAB,S3BQAB,S3DQAB,
     &                              S4AQAB,S4BQAB,S4DQAB,MREG_TAB,
     &                              MREG_TBA,UMAT_VTA,MSSTA,TSSTA,MSSQA,
     &                              TAUQA,MSSTB,TSSTB,UMAT_VTB,MSSQB,
     &                              TAUQB,SIG0Q_SURF,NSPINPROJ,SIG_MODE)
C
C-----------------------------------------------------------------------
C             SIG 1   --   without vertex corrections
C-----------------------------------------------------------------------
C
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
                     if( use_jtauk_square ) then
                       CALL SIG1_SURF_JTAUKSQUARE(
     &                        ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQZ,
     &                        MSSQA,MSSQB,
     &                        MAQAB,MBQAB,MCQAB,MDQAB,
     &                        SIG1Q_SURF_NV,'N',MAQBA,MCQBA,MBQBA,
     &                        MDQBA,NSPINPROJ)
                     elseif( not_store_chiz ) then
                       CALL SIG1_SURF_NOTSTORECHIZ(
     &                        ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQZ,
     &                        MSSQA,MSSQB,
     &                        MAQAB,MBQAB,MCQAB,MDQAB,
     &                        SIG1Q_SURF_NV,'N',MAQBA,MCQBA,MBQBA,
     &                        MDQBA,NSPINPROJ)
                     else
                       CALL SIG1_SURF(MAQAB,MBQAB,MCQAB,MDQAB,
     &                        SIG1Q_SURF_NV,'N',MAQBA,MCQBA,MBQBA,
     &                        MDQBA,NSPINPROJ)
                     endif
c end-mod-xjq
C-----------------------------------------------------------------------
C                              Streda term
C-----------------------------------------------------------------------
C
                     IF ( SIG_MODE.EQ.'STREDA    ' .AND. (1.EQ.1) .AND. 
     &                    (.NOT.PUBLIC_VERSION) )
     &                    CALL SIGOFF(SIGMA_INTEGRATION,MSSTB,TSSTB,
     &                    MSSQB,TAUQB,SIGOFFQ,IFIL_RHSB,NSPINPROJ,
     &                    MEOFF_BARG_T,MEOFF_AHE_T,ME_SZ_T,X_VFT,NVFO_Q,
     &                    IVFT_VFOQ,NVFTMAX)
C
C-----------------------------------------------------------------------
C                        vertex corrections
C-----------------------------------------------------------------------
C
c modified by XJQ: debug
                     if ( VERTEX .and. ncpa > 0 ) then
c end-mod-xjq
C
                        CALL LINRESP_VERTEX(TAUQA,TAUQB,TSSTA,TSSTB,
     &                     UMAT_VTA,UMAT_VTB,MSSQA,MSSQB)
C
C-----------------------------------------------------------------------
C             SIG 1   --   with vertex corrections
C-----------------------------------------------------------------------
C
                        CALL SIG1_SURF(MAQAB,MBQAB,MCQAB,MDQAB,
     &                                 SIG1Q_SURF_VC,'V',MAQBA,MCQBA,
     &                                 MBQBA,MDQBA,NSPINPROJ)
C
                     ELSE
C
                        SIG1Q_SURF_VC(:,:,:,:,:)
     &                     = SIG1Q_SURF_NV(:,:,:,:,:)
C
                     END IF
C
C
                  END IF
C
C ======================================================================
C
                  IPRINT = IPRINTL
C
C=======================================================================
C
               END IF ! MPI
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            END DO LOOP_IEB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
         END DO LOOP_IEA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C=======================================================================
C                         collect results
C=======================================================================
C
         SIG0Q(:,:,:,:) = SIG0Q_SURF(:,:,:,:) + SIG0Q_FSEA(:,:,:,:)
C
         SIG1Q_NV(:,:,:,:,:) = SIG1Q_SURF_NV(:,:,:,:,:)
     &                         + SIG1Q_FSEA_NV(:,:,:,:,:)
C
         SIG1Q_VC(:,:,:,:,:) = SIG1Q_SURF_VC(:,:,:,:,:)
     &                         + SIG1Q_FSEA_VC(:,:,:,:,:)
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
            IF ( SIGMA_PROJECT ) THEN
C
               IF ( ILOOP_SIGMA_PROJECT.EQ.1 ) THEN
                  WRITE (6,99004) 
     &                       'total conductivity SIG for J(+-)=J(-+)=0 '
               ELSE IF ( ILOOP_SIGMA_PROJECT.EQ.2 ) THEN
                  WRITE (6,99004) 'projected conductivity SIG(++) '
               ELSE IF ( ILOOP_SIGMA_PROJECT.EQ.3 ) THEN
                  WRITE (6,99004) 'projected conductivity SIG(--) '
               END IF
C
            END IF
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
            LOOP_SPR:DO ISPR = 1,NSPR
C
               ISP = LIST_ISPR(ISPR)
               WRITE (6,99007)
               WRITE (6,'(/,A80)') STR_ISP_PROJ(ISP)
               WRITE (6,99007)
C
               IF ( NEA.GT.1 ) THEN
C
                  WRITE (6,'(/,"<SIG> printing SIG0Q_FSEA")')
                  CALL PR_COND_TENSOR(SIG0Q_FSEA(1,1,ISP,1),1D0)
C
                  WRITE (6,'(/,"<SIG> printing SIG1Q_FSEA_NV")')
                  CALL PR_COND_TENSOR(SIG1Q_FSEA_NV(1,1,ISP,1,1),1D0)
C
                  WRITE (6,'(/,"<SIG> printing SIG1Q_FSEA_VC")')
                  CALL PR_COND_TENSOR(SIG1Q_FSEA_VC(1,1,ISP,1,1),1D0)
               END IF
C
               WRITE (6,'(40("+-"))')
               WRITE (6,'(/,"<SIG> printing FSURF contribs")')
               WRITE (6,'(/,"<SIG> printing SIG0Q_SURF")')
               CALL PR_COND_TENSOR(SIG0Q_SURF(1,1,ISP,1),1D0)
C
               WRITE (6,'(/,"<SIG> printing SIG1Q_SURF_NV")')
               CALL PR_COND_TENSOR(SIG1Q_SURF_NV(1,1,ISP,1,1),1D0)
C
               WRITE (6,'(/,"<SIG> printing SIG1Q_SURF_VC")')
               CALL PR_COND_TENSOR(SIG1Q_SURF_VC(1,1,ISP,1,1),1D0)
            END DO LOOP_SPR
C
C
            IPRINTL = IPRINT
            IPRINT = 5
            WRITE (6,*) 'Complete conductivity'
C
C-------------------------------------------------------------------- AC
            IF ( SIG_MODE(1:6).EQ.'OPTICS' ) THEN
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C             collect results for energy points of a E-path
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
               IF ( MPI_ELOOP ) THEN
C
                  CALL DRV_MPI_BARRIER
C
                  IDIMS = 3*3*NSPINPROJ*NQMAX*NOMEGA
                  ALLOCATE (CWORK(IDIMS))
                  CALL DRV_MPI_REDUCE_C(SIG0Q_OPT(1,1,1,1,1),CWORK(1),
     &                                  IDIMS)
C
                  DEALLOCATE (CWORK)
                  IDIMS = 3*3*NSPINPROJ*NQMAX*NQMAX*NOMEGA
                  ALLOCATE (CWORK(IDIMS))
                  CALL DRV_MPI_REDUCE_C(SIG1Q_OPT_NV(1,1,1,1,1,1),
     &                                  CWORK(1),IDIMS)
                  CALL DRV_MPI_REDUCE_C(SIG1Q_OPT_VC(1,1,1,1,1,1),
     &                                  CWORK(1),IDIMS)
C
               END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
               IF ( MPI_ID.EQ.0 )
     &              CALL SIGSUM_OPTICS(SIG0Q_OPT,SIG1Q_OPT_NV,
     &              SIG1Q_OPT_VC,NSPINPROJ,NOMEGA,OMEGATAB,I_TEMP_LAT)
C
C-------------------------------------------------------------------- AC
C-------------------------------------------------------------------- DC
            ELSE
C
               CALL SIGSUM(SIG1Q_FSEA_NV,SIG1Q_FSEA_VC,SIG0Q_FSEA,
     &                     SIG1Q_SURF_NV,SIG1Q_SURF_VC,SIG0Q_SURF,
     &                     SIGMAAU_NV,SIGMAAU_VC,RHOAU_NV,RHOAU_VC,
     &                     SIGMA_NV,SIGMA_VC,RHO_NV,RHO_VC,NSPINPROJ,
     &                     SIGOFFQ,SIG_MODE,NEA,I_TEMP_LAT,CHANGE_FRAME,
     &                     MROT_FRAME)
C
C-------------------------------------------------------------------- DC
            END IF
C
C
            IPRINT = IPRINTL
C
         END IF
C
         if(mpi_id==0 .and. prt_sig_realz) then
           rho_diag_si(1) = 1e8*rho_vc(1,1)/consi
           rho_diag_si(2) = 1e8*rho_vc(2,2)/consi
           rho_diag_si(3) = 1e8*rho_vc(3,3)/consi
           write(iotmp2,'(2f9.6,5f11.6)') etab(ie_sig,1),
     &       rho_diag_si(1),rho_diag_si(2),rho_diag_si(3),
     &       sum(rho_diag_si(1:3))/3d0,
     &       3d0/sum((1/rho_diag_si(1:3)))
         endif
c
         enddo ! ie_sig
c
         if(mpi_id==0 .and. prt_sig_realz) close(iotmp2)
C
         IF( MPI_ID .EQ. 0 ) THEN
            CALL CPU_TIME(TIME)
            WRITE (6,99003) TIME - TIME0
         END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO LOOP_TEMP_LAT
C
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== END =
C=======================================================================
C=======================================================================
C
      DEALLOCATE (MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA)
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      IF ( SIGMA_PROJECT ) THEN
C
         IPROJ = ILOOP_SIGMA_PROJECT
C
         SIG0Q_SURF_PROJ(1:3,1:3,1:NQ,IPROJ) = SIG0Q(1:3,1:3,1,1:NQ)
         SIG1Q_SURF_NV_PROJ(1:3,1:3,1:NQ,1:NQ,IPROJ)
     &      = SIG1Q_NV(1:3,1:3,1,1:NQ,1:NQ)
         SIG1Q_SURF_VC_PROJ(1:3,1:3,1:NQ,1:NQ,IPROJ)
     &      = SIG1Q_VC(1:3,1:3,1,1:NQ,1:NQ)
C
         SIG1Q_FSEA_NV(1:3,1:3,1:NSPINPROJ,1:NQ,1:NQ) = 0D0
         SIG1Q_FSEA_VC(1:3,1:3,1:NSPINPROJ,1:NQ,1:NQ) = 0D0
C
         IF ( ILOOP_SIGMA_PROJECT.LT.3 ) GOTO 100
C
C-----------------------------------------------------------------------
C                     print spin flip results
C-----------------------------------------------------------------------
C
         SIG0Q(1:3,1:3,1,1:NQ) = SIG0Q_SURF_PROJ(1:3,1:3,1:NQ,1)
     &                           - SIG0Q_SURF_PROJ(1:3,1:3,1:NQ,2)
     &                           - SIG0Q_SURF_PROJ(1:3,1:3,1:NQ,3)
C
         SIG1Q_NV(1:3,1:3,1,1:NQ,1:NQ)
     &      = SIG1Q_SURF_NV_PROJ(1:3,1:3,1:NQ,1:NQ,1)
     &      - SIG1Q_SURF_NV_PROJ(1:3,1:3,1:NQ,1:NQ,2)
     &      - SIG1Q_SURF_NV_PROJ(1:3,1:3,1:NQ,1:NQ,3)
C
         SIG1Q_VC(1:3,1:3,1,1:NQ,1:NQ)
     &      = SIG1Q_SURF_VC_PROJ(1:3,1:3,1:NQ,1:NQ,1)
     &      - SIG1Q_SURF_VC_PROJ(1:3,1:3,1:NQ,1:NQ,2)
     &      - SIG1Q_SURF_VC_PROJ(1:3,1:3,1:NQ,1:NQ,3)
C
         WRITE (6,99004) 'spin-fip contributions SIG(+-) + SIG (-+)'
C
         CALL SIGSUM(SIG1Q_FSEA_NV,SIG1Q_FSEA_VC,SIG0Q,SIG1Q_SURF_NV,
     &               SIG1Q_SURF_VC,SIG0Q_SURF,SIGMAAU_NV,SIGMAAU_VC,
     &               RHOAU_NV,RHOAU_VC,SIGMA_NV,SIGMA_VC,RHO_NV,RHO_VC,
     &               NSPINPROJ,SIGOFFQ,SIG_MODE,NEA,I_TEMP_LAT,
     &               CHANGE_FRAME,MROT_FRAME)
C
      END IF
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C              loop to get the SPIN-projected conductivities         END
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
C   ====================================================================
99001 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,25X,'EFERMI =',F10.5,/,10X,'IEB  =',I4,5X,
     &        'ERYDB  =',2F10.5)
99002 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*            ****   ***   ****   *     *    ***              *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  **   *    *   *             *'
     &  ,/,10X,
     &  '*           *        *   *       * * * *  *     *            *'
     &  ,/,10X,
     &  '*            ****    *   *  ***  *  *  *  *******            *'
     &  ,/,10X,
     &  '*                *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*            ****   ***   ****   *     *  *     *            *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99003 FORMAT (/,5X,'execution time for <SIG>:',F14.3,' secs',/)
99004 FORMAT (//,2(/,1X,79('P')),/,10X,A,2(/,1X,79('P')),/)
99005 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,25X,'EFERMI =',F10.5,/,10X,'IEB  =',I4,5X,
     &        'ERYDB  =',F10.5,F20.15,/,10X,'IEA  =',I4,5X,'ERYDA  =',
     &        F10.5,F20.15/,10X,'IE   =',I4,5X,:,'ERYD   =',2F10.5)
99006 FORMAT (/,79('*'),/,79('*'),/,30X,'ENERGY LOOP SIGMA',/,79('*'),/,
     &        79('*'),/)
99007 FORMAT (/,40('*'),/,40('*'))
99008 FORMAT (10X,A6,'for spinpol:  alpha-part will be switched off ')
99009 FORMAT (10X,A6,'for spinpol:  nabla-part will be switched off ')
99010 FORMAT (//,10X,'optical conductivity calculated for ',//,10X,
     &        'initial  NEB =',I3,'  E =',F10.5,' -- ',F10.5,' Ry',/,
     &        10X,'final    NEA =',I3,'  E =',F10.5,' -- ',F10.5,' Ry',
     &        /,10X,'omega    NOM =',I3,'  E =',F10.5,' -- ',F10.5,
     &        ' Ry',/,10X,'                   E =',F10.5,' -- ',F10.5,
     &        ' eV',//,10X,'EIMAG    ',F10.5,' Ry  ',F10.5,' eV',/,10X,
     &        'GAMMA    ',F10.5,' Ry  ',F10.5,' eV',/,10X,'EPSILON  ',
     &        F10.5,' Ry  ',F10.5,' eV',//,10X,'NSPR     ',I3,/,10X,
     &        'LIST_ISPR ',20I2)
99011 FORMAT (/,10X,'MPI = ',L1,5X,'MPI_ELOOP = ',L1,5X,'MPI_KLOOP = ',
     &        L1,//,10X,'SINGLE_SITE_TEST = ',L1,//)
      END
C*==pr_cond_tensor.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE PR_COND_TENSOR(TENSOR,SCALEFAC)
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--PR_COND_TENSOR1530
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      REAL*8 SCALEFAC
      COMPLEX*16 TENSOR(3,3)
C
C Local variables
C
      CHARACTER*100 FMT1
      INTEGER I,J
      COMPLEX*16 MAT(3,3),RD
C
C*** End of declarations rewritten by SPAG
C
      FMT1 = '(E14.6,'' '',E10.4)'
C
      DO I = 1,3
         DO J = 1,3
            IF ( ABS(TENSOR(I,J)).LT.TOL ) THEN
               MAT(I,J) = C0
            ELSE
               RD = TENSOR(I,J)
               IF ( ABS(DIMAG(RD)).LT.TOL ) RD = DCMPLX(DREAL(RD),0D0)
               IF ( ABS(DREAL(RD)).LT.TOL ) RD = DCMPLX(0D0,DIMAG(RD))
               MAT(I,J) = RD
            END IF
         END DO
      END DO
C
C      WRITE (6,'(3X,84("-"))')
      DO I = 1,3
         WRITE (6,'(A3)',ADVANCE='NO') '  ('
         IF ( ABS(MAT(I,1)).LT.TOL ) THEN
            WRITE (6,'(11X,A1,13X)',ADVANCE='NO') '0'
         ELSE
            WRITE (6,FMT=FMT1,ADVANCE='NO') MAT(I,1)*SCALEFAC
         END IF
         WRITE (6,'(A4)',ADVANCE='NO') ' | '
         IF ( ABS(MAT(I,2)).LT.TOL ) THEN
            WRITE (6,'(11X,A1,13X)',ADVANCE='NO') '0'
         ELSE
            WRITE (6,FMT=FMT1,ADVANCE='NO') MAT(I,2)*SCALEFAC
         END IF
         WRITE (6,'(A4)',ADVANCE='NO') ' | '
         IF ( ABS(MAT(I,3)).LT.TOL ) THEN
            WRITE (6,'(11X,A1,13X)',ADVANCE='NO') '0'
         ELSE
            WRITE (6,FMT=FMT1,ADVANCE='NO') MAT(I,3)*SCALEFAC
         END IF
         WRITE (6,'(A2)') ' )'
      END DO
C      WRITE (6,'(3X,84("-"))')
C
      END
C*==check_umat.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE CHECK_UMAT(ERYDA)
      USE MOD_THERMAL,ONLY:UMAT_VT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,LMAT3,WKM1,WKM2,IPIVKM,WKM3,LMAT2,
     &    RREL
      IMPLICIT NONE
C*--CHECK_UMAT1612
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYDA
C
C Local variables
C
      COMPLEX*16 ERYD,LMAT(:,:),XM(:,:),XMD(:,:)
      INTEGER IVT,IVTM,K,M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE LMAT,XM,XMD
C
      ALLOCATE (LMAT(NKMMAX,NKMMAX))
      ALLOCATE (XM(NKMMAX,NKMMAX))
      ALLOCATE (XMD(NKMMAX,NKMMAX))
C
      XM = MATMUL(TRANSPOSE(RREL),RREL)
      XMD = TRANSPOSE(DCONJG(XM))
C
      N = NKM
      M = NKMMAX
C
      IF ( 1.EQ.0 ) THEN
C---------------------------------------------------------- non-rel repr
         LMAT = LMAT2
         K = 2
      ELSE
C----------------------------------------------------------     rel repr
         LMAT = LMAT3
         K = 3
      END IF
      CALL CMATSTRUCT('LMAT',LMAT,N,M,K,K,0,1D-8,1)
C
      ERYD = ERYDA
      WRITE (6,*) ERYD
C
      CALL THERMAL_INIT_UFMAT(ERYD)
C
C      IVTM = 7
      IVTM = 1
      WKM1 = UMAT_VT(:,:,IVTM)
      CALL CMATSTRUCT('u(-s)',WKM1,N,M,K,K,0,1D-8,1)
C----------------------------------------------------------- change svec
C      IVT = 14
      IVT = 4
      WKM1 = UMAT_VT(:,:,IVT)
      CALL CMATSTRUCT('u(s)',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM2 = UMAT_VT(:,:,IVTM)
      WKM3 = MATMUL(WKM1,WKM2)
      CALL CMATSTRUCT('u(s) u(-s)',WKM3,N,M,K,K,0,1D-8,1)
C
      WKM2 = DCONJG(TRANSPOSE(UMAT_VT(:,:,IVT)))
      WKM3 = MATMUL(WKM1,WKM2)
      CALL CMATSTRUCT('u(s) u(s)^+',WKM3,N,M,K,K,0,1D-8,1)
C
      WKM3 = UMAT_VT(:,:,IVT)
      CALL CMATINV3(N,M,IPIVKM,WKM3,WKM1,WKM2)
      CALL CMATSTRUCT('u(s)^(-1)',WKM2,N,M,K,K,0,1D-8,1)
C
C
      WKM1 = TRANSPOSE(UMAT_VT(:,:,IVT))
      CALL CMATSTRUCT('u(s)^T',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM1 = TRANSPOSE(DCONJG(UMAT_VT(:,:,IVT)))
      CALL CMATSTRUCT('u(s)^+',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM2 = TRANSPOSE(UMAT_VT(:,:,IVT))
      WKM1 = MATMUL(LMAT,(MATMUL(WKM2,LMAT)))
      CALL CMATSTRUCT('L u(s)^T L',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM1 = MATMUL(LMAT,
     &       (MATMUL(DCONJG(TRANSPOSE(UMAT_VT(:,:,IVT))),LMAT)))
      CALL CMATSTRUCT('L u(s)^+ L',WKM1,N,M,K,K,0,1D-8,1)
C
C-------------------------------------------------------- considering z
      ERYD = DCMPLX(DREAL(ERYDA),DIMAG(ERYDA))
      WRITE (6,*) ERYD
C
      CALL THERMAL_INIT_UFMAT(ERYD)
C
      WKM1 = UMAT_VT(:,:,IVT)
      CALL CMATSTRUCT('u(z)',WKM1,N,M,K,K,0,1D-8,1)
C
      IF ( K.EQ.2 ) THEN
C
         WKM2 = DCONJG(UMAT_VT(:,:,IVT))
         WKM1 = MATMUL(LMAT,MATMUL(WKM2,LMAT))
         CALL CMATSTRUCT('L (u(z))^* L',WKM1,N,M,K,K,0,1D-8,1)
C
      ELSE IF ( K.EQ.3 ) THEN
C
         WKM1 = DCONJG(UMAT_VT(:,:,IVT))
C
         WKM1 = MATMUL(XMD,MATMUL(WKM1,XM))
         WKM1 = MATMUL(LMAT,MATMUL(WKM1,LMAT))
C
         CALL CMATSTRUCT('L (A^T A)^+ U(z)^* (A^T A) L',WKM1,N,M,K,K,0,
     &                   1D-8,1)
      END IF
C
C-------------------------------------------------------- considering z*
      ERYD = DCMPLX(DREAL(ERYDA),-DIMAG(ERYDA))
      WRITE (6,*) ERYD
C
      CALL THERMAL_INIT_UFMAT(ERYD)
C
      WKM1 = UMAT_VT(:,:,IVT)
      CALL CMATSTRUCT('u(z*)',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM1 = MATMUL(LMAT,MATMUL(TRANSPOSE(WKM1),LMAT))
      CALL CMATSTRUCT('L u(z*)^T L',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM2 = TRANSPOSE(DCONJG(UMAT_VT(:,:,IVT)))
      WKM1 = MATMUL(LMAT,MATMUL(WKM2,LMAT))
      CALL CMATSTRUCT('L u(z*)^+ L',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM2 = UMAT_VT(:,:,IVT)
      WKM1 = MATMUL(LMAT,MATMUL(WKM2,LMAT))
      CALL CMATSTRUCT('L u(z*) L',WKM1,N,M,K,K,0,1D-8,1)
C
C------------------------------------------------------------ check RRLM
      CALL CMATSTRUCT('RREL',RREL,N,M,K,K,0,1D-8,1)
C
      WKM1 = MATMUL(TRANSPOSE(DCONJG(RREL)),RREL)
      CALL CMATSTRUCT('RREL^+ RREL',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM1 = MATMUL(RREL,TRANSPOSE(RREL))
      CALL CMATSTRUCT('RREL RREL^T',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM1 = MATMUL(DCONJG(RREL),TRANSPOSE(DCONJG(RREL)))
      CALL CMATSTRUCT('RREL^* RREL^-1',WKM1,N,M,K,K,0,1D-8,1)
C
      WKM2 = MATMUL(TRANSPOSE(RREL),RREL)
      CALL CMATSTRUCT('RREL^T RREL',WKM2,N,M,K,K,0,1D-8,1)
C
      WKM3 = MATMUL(TRANSPOSE(DCONJG(RREL)),DCONJG(RREL))
      CALL CMATSTRUCT('RREL^+ RREL^*',WKM3,N,M,K,K,0,1D-8,1)
C
      WKM3 = MATMUL(RREL,LMAT)
      CALL CMATSTRUCT('RREL LMAT',WKM3,N,M,K,K,0,1D-8,1)
C
      WKM3 = MATMUL(LMAT,RREL)
      CALL CMATSTRUCT('LMAT RREL',WKM3,N,M,K,K,0,1D-8,1)
C
      ERYD = ERYDA
      WRITE (6,*) ERYD
      CALL THERMAL_INIT_UFMAT(ERYD)
C
      WKM1 = TRANSPOSE(UMAT_VT(:,:,IVT))
C
      IF ( K.EQ.3 ) THEN
         WKM1 = MATMUL(XMD,MATMUL(WKM1,XM))
         CALL CMATSTRUCT('(A^T A)^+ U^T (A^T A)',WKM1,N,M,K,K,0,1D-8,1)
      END IF
C
      END
