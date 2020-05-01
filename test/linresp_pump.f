C*==linresp_pump.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_PUMP
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the conductivity tensor                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,CHIZ
      USE MOD_CPA,ONLY:ITCPAMAX
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,MSST,SSST,TAUT,WKM1,WKM2,NMEMAX,
     &    NCPLWF,IDOS,ISMT
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT,IQTOP,ITOQ,NOQ
      USE MOD_TYPES,ONLY:NT,NTMAX,ITBOT,ITTOP,CONC,NCPLWFMAX,IKMCPLWF,
     &    IMT
      USE MOD_CALCMODE,ONLY:ORBPOL,DMFT,LDAU,KKRMODE,GF_CONV_RH,
     &    SOLVER_FP,LHS_SOL_EQ_RHS_SOL
      USE MOD_ENERGY,ONLY:ETAB,WETAB,NETAB,EFERMI,NEMAX
      USE MOD_FILES,ONLY:IFILCBWF0,IFILCBWF_INCREMENT_ENERGY,RECLNGWF,
     &    RECLNGWF_SPH,IPRINT,WRTAUMQ
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,MPI_ELOOP,MPI_ID
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY_LINRESP
      USE MOD_CONSTANTS,ONLY:PI,CI,C1,C0
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      IMPLICIT NONE
C*--LINRESP_PUMP28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_PUMP')
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
      INTEGER NOBS
      PARAMETER (NOBS=2)
C
C Local variables
C
      LOGICAL CALCINT,CHANGE_FRAME,CHECK,CHECK_OBS_ME,GETIRRSOL
      COMPLEX*16 CAUXA,CAUXB,CDOBS_T(:,:,:),COBS_T(:,:,:),DERYD,EC_BOT,
     &           EC_TOP,ERYD_ADV,ERYD_RET,ETAB_ADV(:),ETAB_C(:),
     &           ETAB_RET(:),E_EXTEND_BOT,E_EXTEND_TOP,E_REFERENCE,
     &           JFL_ADV(:,:,:),JFR_RET(:,:,:),JGL_ADV(:,:,:),
     &           JGR_RET(:,:,:),MEZJ_ADV(:,:,:,:),MEZJ_AR(:,:,:,:,:,:),
     &           MEZJ_RA(:,:,:,:,:,:),MEZJ_RET(:,:,:,:),
     &           MEZZ_ADV(:,:,:,:),MEZZ_AR(:,:,:,:,:,:),
     &           MEZZ_RA(:,:,:,:,:,:),MEZZ_RET(:,:,:,:),MSSQ_ADV(:,:,:),
     &           MSSQ_RET(:,:,:),MSST_ADV(:,:,:),MSST_RET(:,:,:),
     &           OTIL_0(:,:,:),PHASK_DUMMY(:),P_ADV,P_RET,TAUQZ(:,:,:,:)
     &           ,TAUQ_ADV(:,:,:),TAUQ_RET(:,:,:),TAUT_ADV(:,:,:,:),
     &           TAUT_RET(:,:,:,:),TSSQ_ADV(:,:,:),TSSQ_RET(:,:,:),
     &           TSST_ADV(:,:,:),TSST_RET(:,:,:),WKMOBS(:,:,:),
     &           ZFL_ADV(:,:,:),ZFR_RET(:,:,:),ZGL_ADV(:,:,:),
     &           ZGR_RET(:,:,:)
      COMPLEX*16 CMATTRC
      REAL*8 CPACHNG,DMFTMIX,E_EPSILON,MROT_FRAME(3,3),TIME,TIME0
      INTEGER IA_ERR,ICPACONV,ICPAFLAG,IE,IE_0,IE_0_BOT,IE_0_TOP,IE_ADV,
     &        IE_RET,IFIL_LHSA,IFIL_LHSB,IFIL_RHSA,IFIL_RHSB,IFIL_SPHA,
     &        IFIL_SPHB,IM,IOBS,IPRINTL,IPROCE(:),IQ,IRTOP,IT,ITCPA,
     &        IT_ADV,IT_RET,IWRI,IWRIRRWF,IWRREGWF,M,N,NEC,NE_ADV,
     &        NE_RET,NE_RET_BOT,NE_RET_TOP
C
C*** End of declarations rewritten by SPAG
C
C      PARAMETER (E_EPSILON=1.0D-3)
C     PARAMETER (CHECK_OBS_ME=.FALSE.)
C
      DATA CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/
C
      ALLOCATABLE IPROCE
C
      ALLOCATABLE TAUQZ
      ALLOCATABLE ETAB_RET,ETAB_ADV,ETAB_C,WKMOBS
      ALLOCATABLE PHASK_DUMMY
      ALLOCATABLE TAUQ_RET,MSSQ_RET,TSSQ_RET,TAUT_RET,TSST_RET,MSST_RET
      ALLOCATABLE TAUQ_ADV,MSSQ_ADV,TSSQ_ADV,TAUT_ADV,TSST_ADV,MSST_ADV
      ALLOCATABLE CDOBS_T,COBS_T,OTIL_0
      ALLOCATABLE MEZJ_RET,MEZZ_RET,MEZJ_ADV,MEZZ_ADV
      ALLOCATABLE ZGR_RET,ZFR_RET,JGR_RET,JFR_RET
      ALLOCATABLE ZGL_ADV,ZFL_ADV,JGL_ADV,JFL_ADV
      ALLOCATABLE MEZJ_RA,MEZJ_AR,MEZZ_RA,MEZZ_AR
C
      CALL UNDER_CONSTRUCTION(ROUTINE)
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      ALLOCATE (TAUQ_RET(M,M,NQMAX),TAUQ_ADV(M,M,NQMAX))
      ALLOCATE (MSSQ_RET(M,M,NQMAX),MSSQ_ADV(M,M,NQMAX))
      ALLOCATE (TSSQ_RET(M,M,NQMAX),TSSQ_ADV(M,M,NQMAX))
      ALLOCATE (MEZZ_RET(M,M,NTMAX,NMEMAX),MEZZ_ADV(M,M,NTMAX,NMEMAX))
      ALLOCATE (MEZJ_RET(M,M,NTMAX,NMEMAX),MEZJ_ADV(M,M,NTMAX,NMEMAX))
C
      ALLOCATE (ZGR_RET(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFR_RET(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGR_RET(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFR_RET(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGRA')
C
      ALLOCATE (ZGL_ADV(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFL_ADV(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGL_ADV(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFL_ADV(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGLB')
C
      ALLOCATE (TAUQZ(NKMMAX,NKMMAX,NQMAX,2),STAT=IA_ERR)
C----------------------------------------------------------- dummy array
      ALLOCATE (PHASK_DUMMY(1))
C
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc MAQAB')
C
      CHECK = .FALSE.
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
      N = NKM
      M = NKMMAX
C
      ALLOCATE (TSST_RET(M,M,NTMAX),TSST_ADV(M,M,NTMAX))
      ALLOCATE (MSST_RET(M,M,NTMAX),MSST_ADV(M,M,NTMAX),STAT=IA_ERR)
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( DMFT .OR. LDAU ) CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),
     &     DMFTMIX,NEMAX)
C
C ======================================================================
      WRITE (6,99003)
C ======================================================================
C
      KKRMODE = 'SINGLE-SITE'
C
C ======================================================================
C
      CALL CPU_TIME(TIME0)
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
C                          set parameters
C ======================================================================
C
      CALL INIT_MROT_FRAME(CHANGE_FRAME,MROT_FRAME)
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      CALL SECTION_SET_INTEGER('NEA',NE_RET,9999,0)
C
      EC_BOT = 0.0D0
      EC_TOP = EFERMI
C
      E_EXTEND_BOT = 0.2D0
      E_EXTEND_TOP = 1.0D0
C
      NEC = 50
      NEC = 3
      DERYD = (EC_TOP-EC_BOT)/DBLE(NEC)
C
      NE_RET_BOT = INT(E_EXTEND_BOT/DERYD) + 1
      NE_RET_TOP = INT(E_EXTEND_TOP/DERYD) + 1
C
      NE_RET = NE_RET_BOT + NEC + NE_RET_TOP
      NE_ADV = NE_RET
C
      IE_0_BOT = NE_RET_BOT + 1
      IE_0_TOP = NE_RET_BOT + NEC
C
      E_REFERENCE = EC_BOT - DERYD*(NE_RET_BOT+0.5D0)
C
      NEMAX = MAX(NE_RET,NE_ADV)
C
      IF ( ALLOCATED(ETAB) ) DEALLOCATE (ETAB,WETAB)
      ALLOCATE (ETAB(NEMAX,2),WETAB(NEMAX,2))
C
      ALLOCATE (OTIL_0(M,M,NE_RET),WKMOBS(M,M,NOBS))
      ALLOCATE (TAUT_RET(M,M,NTMAX,NE_RET),TAUT_ADV(M,M,NTMAX,NE_RET))
      ALLOCATE (MEZZ_RA(NKMMAX,NKMMAX,3,NOBS,NE_RET,NE_RET))
      ALLOCATE (MEZJ_RA(NKMMAX,NKMMAX,3,NOBS,NE_RET,NE_RET))
      ALLOCATE (MEZZ_AR(NKMMAX,NKMMAX,3,NOBS,NE_RET,NE_RET))
      ALLOCATE (MEZJ_AR(NKMMAX,NKMMAX,3,NOBS,NE_RET,NE_RET))
      ALLOCATE (CDOBS_T(NE_RET,NOBS,NTMAX),COBS_T(NE_RET,NOBS,NTMAX))
      CDOBS_T(:,:,:) = C0
      COBS_T(:,:,:) = C0
C
      ALLOCATE (ETAB_RET(NE_RET),ETAB_ADV(NE_RET),ETAB_C(NE_RET))
      ALLOCATE (IPROCE(NE_RET))
      IPROCE(1:NE_RET) = 0
C
      E_EPSILON = 1.0D-6
C###      E_EPSILON = 0D0
C
      DO IE = 1,NE_RET
C
         ETAB_RET(IE) = E_REFERENCE + DERYD*IE + CI*E_EPSILON
         ETAB_ADV(IE) = E_REFERENCE + DERYD*IE - CI*E_EPSILON
         ETAB_C(IE) = E_REFERENCE + DERYD*IE
C
      END DO
C
      WRITE (6,99008) 'DERYD           ',DREAL(DERYD)
      WRITE (6,99008) 'E_REFERENCE     ',E_REFERENCE
      WRITE (6,99009) 'nea_bot         ',NE_RET_BOT
      WRITE (6,99008) 'E(1)            ',DREAL(E_REFERENCE+DERYD*1)
      WRITE (6,99008) 'E(nea_bot)      ',E_REFERENCE + DERYD*NE_RET_BOT
      WRITE (6,99009) 'iec_bot         ',IE_0_BOT
      WRITE (6,99008) 'E(iec_bot)      ',E_REFERENCE + DERYD*IE_0_BOT
      WRITE (6,99008) 'E(iec_bot)      ',ETAB_C(IE_0_BOT)
      WRITE (6,99009) 'iec_top         ',IE_0_TOP
      WRITE (6,99008) 'E(iec_top)      ',E_REFERENCE + DERYD*IE_0_TOP
      WRITE (6,99008) 'E(iec_top)      ',ETAB_C(IE_0_TOP)
      WRITE (6,99008) 'EFERMI          ',EFERMI
      WRITE (6,99008) 'E(iec_top+1)    ',
     &                E_REFERENCE + DERYD*(IE_0_TOP+1)
      WRITE (6,99008) 'E(iec_top+1)    ',ETAB_C(IE_0_TOP+1)
      WRITE (6,99008) 'E(nea)          ',E_REFERENCE + DERYD*NE_RET
      WRITE (6,99009) 'nea             ',NE_RET
C
C-----------------------------------------------------------------------
C                 see MOD_FILES for the setting scheme
C-----------------------------------------------------------------------
C
      IFIL_RHSA = IFILCBWF0 + 1
      IFIL_RHSB = IFIL_RHSA + IFILCBWF_INCREMENT_ENERGY
C
      OPEN (IFIL_RHSA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      OPEN (IFIL_RHSB,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      IF ( FULLPOT ) THEN
C
         CALL SET_IFIL_SPH(IFIL_RHSA,IFIL_SPHA)
C
         CALL SET_IFIL_SPH(IFIL_RHSB,IFIL_SPHB)
C
         OPEN (IFIL_SPHA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
         OPEN (IFIL_SPHB,STATUS='SCRATCH',FORM='UNFORMATTED',
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
         OPEN (IFIL_LHSA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
         OPEN (IFIL_LHSB,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=RECLNGWF)
C
      END IF
C
      IWRI = 100
C
C***********************************************************************
C***********************************************************************
C****************       START CALCULATION        ***********************
C***********************************************************************
C***********************************************************************
C
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop  START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      WRITE (6,99006)
C
      WRTAUMQ = .FALSE.
C
      ICPAFLAG = 0
      CPACHNG = 0.0D0
C
      ERYDA_EQ_ERYDB = .TRUE.
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      LOOP_IEA:DO IE_RET = 1,NE_RET
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ELOOP .AND. MPI_ID.NE.IPROCE(IE_RET) ) CYCLE LOOP_IEA
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
         ERYD_RET = ETAB_RET(IE_RET)
C
C--------------------------------------------- adjust number of k-points
C        CALL SIG_DRV_KMESH(IEA,NEA,SIG_MODE)
C
C ===================================== solve SS - differential equation
C
         IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &        = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE_RET)
C
C---------------------------------------------------------------- Z(RET)
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSA,GETIRRSOL,
     &                 ERYD_RET,P_RET,IPRINT,TSST_RET,MSST_RET,SSST,
     &                 MEZZ_RET,MEZJ_RET,ORBPOL)
C
C -------------------------------------- read in wavefunctions for R-RET
C
         IT_RET = 1
C
         IM = IMT(IT_RET)
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,IT_RET,1,ZGR_RET,ZFR_RET,
     &                        JGR_RET,JFR_RET,IRTOP,NCPLWF,IKMCPLWF)
C
C ======================================================================
C
         CALL MSSINIT(TSST_RET,MSST_RET,TSSQ_RET,MSSQ_RET)
C
C***********************************************************************
         IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
            CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD_RET,P_RET,IPRINT,
     &                           ITCPA,ICPACONV,CONC,NOQ,ITOQ,
     &                           PHASK_DUMMY,1,NTMAX,TSST_RET,MSST_RET,
     &                           TSSQ_RET,MSSQ_RET,TAUQ_RET)
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
            CALL TAU_TB(ERYD_RET,P_RET,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,
     &                  TSST_RET,MSST_RET,TSSQ_RET,MSSQ_RET,TAUQ_RET)
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
            TAUQ_RET(:,:,:) = TSSQ_RET(:,:,:)
C
C***********************************************************************
         ELSE
C
            CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
         END IF
C***********************************************************************
C
C ======================================================================
C
C
C ======================================================================
         IF ( IPRINT.GE.3 .OR. CHECK )
     &        CALL DUMPTAU(IE,ERYD_RET,IWRI,MSST,MSSQ_RET,TAUT,TAUQ_RET)
C ======================================================================
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
         LOOP_IEB:DO IE_ADV = 1,NE_ADV
C
Cccccc            IF ( IE_ADV.NE.IE_RET ) CYCLE LOOP_IEB
C
            ERYD_ADV = ETAB_ADV(IE_ADV)
C
            WRITE (6,99005) EFERMI,IE_RET,ERYD_RET,IE_ADV,ERYD_ADV
C
C ======================================================================
C ===================================== solve SS - differential equation
C
            IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &           = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
C---------------------------------------------------------------- Z(ADV)
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_RHSB,GETIRRSOL,
     &                    ERYD_ADV,P_ADV,IPRINT,TSST_ADV,MSST_ADV,SSST,
     &                    MEZZ_ADV,MEZJ_ADV,ORBPOL)
C
C -------------------------------------- read in wavefunctions for L-ADV
C
            IT_ADV = IT_RET
C
            CALL WAVFUN_READ_REL(IFIL_RHSB,IT_ADV,1,ZGL_ADV,ZFL_ADV,
     &                           JGL_ADV,JFL_ADV,IRTOP,NCPLWF,IKMCPLWF)
C
C ======================================================================
C
            CALL MSSINIT(TSST_ADV,MSST_ADV,TSSQ_ADV,MSSQ_ADV)
C
C ======================================================================
C
C
C***********************************************************************
            IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
               CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD_ADV,P_ADV,
     &                              IPRINT,ITCPA,ICPACONV,CONC,NOQ,ITOQ,
     &                              PHASK_DUMMY,1,NTMAX,TSST_ADV,
     &                              MSST_ADV,TSSQ_ADV,MSSQ_ADV,TAUQ_ADV)
C
C***********************************************************************
            ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
               CALL TAU_TB(ERYD_ADV,P_ADV,ICPAFLAG,CPACHNG,ITCPA,
     &                     ICPACONV,TSST_ADV,MSST_ADV,TSSQ_ADV,MSSQ_ADV,
     &                     TAUQ_ADV)
C
C***********************************************************************
            ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
               TAUQ_ADV(:,:,:) = TSSQ_ADV(:,:,:)
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
               IF ( IBZINT.EQ.4 ) THEN
C
                  CALL STOP_MESSAGE(ROUTINE,
     &                   'only <SIGKLOOPS> available for BZ integration'
     &                   )
C
               ELSE IF ( NO_SYMMETRY_LINRESP ) THEN
C
                  CALL SIGKLOOPS_NOSYM(ERYD_RET,P_RET,ERYD_ADV,P_ADV,
     &                                 TAUQ_RET,TAUQ_ADV,TAUQZ,MSSQ_RET,
     &                                 MSSQ_ADV)
C
               ELSE
C
                  CALL SIGKLOOPS(ERYD_RET,P_RET,ERYD_ADV,P_ADV,TAUQ_RET,
     &                           TAUQ_ADV,TAUQZ,MSSQ_RET,MSSQ_ADV)
C
               END IF
C
C***********************************************************************
            ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
               CALL SIGKLOOPSDRV_TB(ERYD_RET,P_RET,MSSQ_RET,TAUQZ)
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
C
            IT_ADV = IT_RET
C
            IF ( IE_ADV.EQ.IE_RET ) THEN
               CHECK_OBS_ME = .TRUE.
            ELSE
               CHECK_OBS_ME = .FALSE.
            END IF
            CHECK_OBS_ME = .FALSE.
C
            CALL CALC_OBS_ME(IT_RET,ZGR_RET,ZFR_RET,JGR_RET,JFR_RET,
     &                       NCPLWF,IKMCPLWF,ZGL_ADV,ZFL_ADV,JGL_ADV,
     &                       JFL_ADV,NCPLWF,IKMCPLWF,MEZZ_RET,MEZJ_RET,
     &                       MEZZ_RA(1,1,1,1,IE_RET,IE_ADV),
     &                       MEZJ_RA(1,1,1,1,IE_RET,IE_ADV),NOBS,'XXX',
     &                       CHECK_OBS_ME)
C
            CALL CALC_OBS_ME(IT_RET,ZGL_ADV,ZFL_ADV,JGL_ADV,JFL_ADV,
     &                       NCPLWF,IKMCPLWF,ZGR_RET,ZFR_RET,JGR_RET,
     &                       JFR_RET,NCPLWF,IKMCPLWF,MEZZ_ADV,MEZJ_ADV,
     &                       MEZZ_AR(1,1,1,1,IE_ADV,IE_RET),
     &                       MEZJ_AR(1,1,1,1,IE_ADV,IE_RET),NOBS,'XXX',
     &                       CHECK_OBS_ME)
C
C MEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEME
            IF ( IE_ADV.EQ.IE_RET ) THEN
C
               IT = 1
               IQ = 1
               TAUT_RET(:,:,IT,IE_RET) = TAUQ_RET(:,:,IQ)
               TAUT_ADV(:,:,IT,IE_ADV) = TAUQ_ADV(:,:,IQ)
C
C
               DO IOBS = 1,NOBS
C
                  CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ_RET(1,1,IT,IOBS),M,
     &                       TAUT_RET(1,1,IT,IE_RET),M,C0,WKM1,M)
C
                  CAUXA = CMATTRC(N,M,WKM1)
C
                  CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ_ADV(1,1,IT,IOBS),M,
     &                       TAUT_ADV(1,1,IT,IE_ADV),M,C0,WKM1,M)
C
                  CAUXB = CMATTRC(N,M,WKM1)
C
                  CDOBS_T(IE_RET,IOBS,IT) = (CAUXA-CAUXB)/(2*CI)
C
                  IF ( IE_RET.GT.1 ) COBS_T(IE_RET,IOBS,IT)
     &                 = COBS_T(IE_RET-1,IOBS,IT)
     &                 + CDOBS_T(IE_RET,IOBS,IT)*DERYD
C
               END DO
            END IF
C MEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEME
C
C
C
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C-----------------------------------------------------------------------
C             matrix elements for energies    E_a   and    E_b
C-----------------------------------------------------------------------
C
Cc               CALL SIGME()
C
C-----------------------------------------------------------------------
C             matrix elements for energies    E_b   and    E_a
C-----------------------------------------------------------------------
C
Cc                  CALL SIGME()
C
               IPRINTL = IPRINT
               IPRINT = 5
C
C-----------------------------------------------------------------------
C             SIG 0   --   site diagonal contributions
C-----------------------------------------------------------------------
C
CC                  CALL SIG0_FSEA()
C
C ======================================================================
C
               IPRINT = IPRINTL
C
C=======================================================================
C
            END IF    ! MPI
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
      DO IE_RET = 1,NE_RET
C
         WRITE (100,99002) DREAL(ETAB_RET(IE_RET)),
     &                     (CDOBS_T(IE_RET,IOBS,IT),
     &                     COBS_T(IE_RET,IOBS,IT),IOBS=1,NOBS)
C
      END DO
C
      CDOBS_T(:,:,:) = C0
      COBS_T(:,:,:) = C0
C
      IT = 1
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      DO IE_RET = 1,NE_RET
C
C00000000000000000000000000000000000000000000000000000000000000000000000
C------------------------------------------------------ O^RET0 * Im TAU0
         DO IE_0 = IE_0_BOT,IE_0_TOP
            WKM1(1:N,1:N) = TAUT_RET(1:N,1:N,IT,IE_0)
            WKM2(1:N,1:N) = WKM1(1:N,1:N) - DCONJG(WKM1(1:N,1:N))
            CALL ZGEMM('N','N',N,N,N,C1,
     &                 MEZZ_RA(1:N,1:N,1,IDOS,IE_RET,IE_0),M,WKM2,M,C0,
     &                 OTIL_0(1:N,1:N,IE_0),M)
         END DO
C00000000000000000000000000000000000000000000000000000000000000000000000
C
         WKMOBS(1:N,1:N,:) = C0
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         DO IE_ADV = 1,NE_RET
C
            WKM1(1:N,1:N) = C0
C00000000000000000000000000000000000000000000000000000000000000000000000
            DO IE_0 = IE_0_BOT,IE_0_TOP
               CALL ZGEMM('N','N',N,N,N,C1,OTIL_0(1,1,IE_0),M,
     &                    MEZZ_RA(1,1,1,IDOS,IE_0,IE_ADV),M,C1,WKM1,M)
            END DO
C00000000000000000000000000000000000000000000000000000000000000000000000
C
            CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,TAUT_ADV(1,1,IT,IE_ADV),
     &                 M,C0,WKM2,M)
C
            CALL ZGEMM('N','N',N,N,N,C1,WKM2,M,
     &                 MEZZ_RA(1,1,1,IDOS,IE_ADV,IE_RET),M,C1,
     &                 WKMOBS(1,1,1),M)
C
            CALL ZGEMM('N','N',N,N,N,C1,WKM2,M,
     &                 MEZZ_RA(1,1,2,ISMT,IE_ADV,IE_RET),M,C1,
     &                 WKMOBS(1,1,2),M)
C
         END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
         DO IOBS = 1,NOBS
            CALL ZGEMM('N','N',N,N,N,C1,TAUT_RET(1,1,IT,IE_RET),M,
     &                 WKMOBS(1,1,IOBS),M,C0,WKM1,M)
C
            CDOBS_T(IE_RET,IOBS,IT) = CMATTRC(N,M,WKM1)
C
            IF ( IE_RET.GT.1 ) COBS_T(IE_RET,IOBS,IT)
     &           = COBS_T(IE_RET-1,IOBS,IT) + CDOBS_T(IE_RET,IOBS,IT)
C
         END DO
C
      END DO
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      CAUXA = (-CI/2D0)*(DERYD/(2D0*PI))**2*(1D0/(2D0*CI))
      CAUXB = (-CI/2D0)*(DERYD/(2D0*PI))**3*(1D0/(2D0*CI))
C
      CDOBS_T(:,:,:) = CDOBS_T(:,:,:)*CAUXA
      COBS_T(:,:,:) = COBS_T(:,:,:)*CAUXB
C
      DO IE_RET = 1,NE_RET
C
         WRITE (200,99002) DREAL(ETAB_RET(IE_RET)),
     &                     (CDOBS_T(IE_RET,IOBS,IT),
     &                     COBS_T(IE_RET,IOBS,IT),IOBS=1,NOBS)
C
      END DO
C
C=======================================================================
C                         collect results
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C
         IPRINTL = IPRINT
         IPRINT = 5
         WRITE (6,*) 'Complete conductivity'
C
C        CALL SIGSUM()
C
         IPRINT = IPRINTL
C
         CALL CPU_TIME(TIME)
         WRITE (6,99004) TIME - TIME0
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      CALL UNDER_CONSTRUCTION(ROUTINE)
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
C   ====================================================================
99001 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,25X,'EFERMI =',F10.5,/,10X,'IEB  =',I4,5X,
     &        'ERYDB  =',2F10.5)
99002 FORMAT (F10.4,8E15.5)
99003 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*              *****   *    *  *     *  *****                *'
     &  ,/,10X,
     &  '*              *    *  *    *  **   **  *    *               *'
     &  ,/,10X,
     &  '*              *    *  *    *  * * * *  *    *               *'
     &  ,/,10X,
     &  '*              *****   *    *  *  *  *  *****                *'
     &  ,/,10X,
     &  '*              *       *    *  *     *  *                    *'
     &  ,/,10X,
     &  '*              *       *    *  *     *  *                    *'
     &  ,/,10X,
     &  '*              *        ****   *     *  *                    *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99004 FORMAT (/,5X,'execution time for <SIG>:',F14.3,' secs',/)
99005 FORMAT (/,79('E'),//,25X,'EFERMI =',F10.5,/,10X,'IEA  =',I4,5X,
     &        'ERYDA  =',F10.5,F15.10,/,10X,'IEB  =',I4,5X,'ERYDB  =',
     &        F10.5,F15.10)
99006 FORMAT (/,79('*'),/,79('*'),/,30X,'ENERGY LOOP SIGMA',/,79('*'),/,
     &        79('*'),/)
99007 FORMAT (/,10X,'MPI = ',L1,5X,'MPI_ELOOP = ',L1,5X,'MPI_KLOOP = ',
     &        L1,//)
99008 FORMAT (10x,A,2F20.12)
99009 FORMAT (10x,A,2I10)
      END
