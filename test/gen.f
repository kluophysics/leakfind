C*==gen.f    processed by SPAG 6.70Rc at 17:10 on 19 Apr 2017
      SUBROUTINE GEN(RUNELOOP,TASK,DOSREP,DOSCORSYS)
C
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine program for       KKRGEN                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CPA,ONLY:CPATOL,NLCPAWRCFG,USENLCPA,NCPA
      USE MOD_RMESH,ONLY:NMMAX,NRMAX,FULLPOT,JRNSMIN,JRNS1,JRCRI,JRWS
      USE MOD_KSPACE,ONLY:IBZINT,NKPTS0
      USE MOD_SITES,ONLY:NQ,IQAT,IQBOT,IQTOP,NLQMAD
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,SEARCHEF,NEPATH,SPLITSS,IGRID,
     &    NETAB,EILOW,EIMAG,EMAX,EMIN,EFERMI,WETAB,PHASA,PHASK,PHAST
      USE MOD_FILES,ONLY:IFILCBWF,LSYSTEM,SYSTEM,LDATSET,DATSET,PLOTPRS,
     &    PLOT2DPR,IPRINT,WRMAT,WRTAU,RDTAU,WRTAUMQ,RDTAUMQ,IFILTAU,
     &    IFILLOG,TAUFIL,WRNOS,CALC_JALF,CALC_JORB,FOUND_SECTION
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX,IMT,VT,BT,TXT_T,NT,ABIT,NLT,
     &    ITBOT,ITTOP,RHOSPN,RHOCHR,RHOORB,MUEORB,MUESPN,NAT,CONC,
     &    NVALTOT,RHO2NS,JDNST,OBS_T,OBS_LT,OBS_TX,OBS_LTX,DOBS_TX,
     &    DOBS_LTX,OBS_T_GLO,OBS_TX_GLO,DOBS_TX_GLO,DOBS_LTEX,
     &    DOBS_TEX_GLO
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,NOBSMAX,NKMMAX,NLMAX,NKMQ,NKM,
     &    NMEHFMAX,NLABIMAX,NL,MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,
     &    TSST,WKM1
      USE MOD_CALCMODE,ONLY:PROGNAME,DMFT,ORBPOL,IREL,ICORE,ISMQHFI,
     &    ITEST,LLOYD,BREITINT,KKRMODE,CLUTYPE,UPDATE_EFERMI,POTFIX,
     &    USE_CONST_POTENTIAL,LDAU,ME_CC_BRA_RWF,THERMAL_VIBRA_FLUCT,
     &    WAVE_FUNCTIONS_AVAILABLE,GF_CONV_RH
      USE MOD_THERMAL,ONLY:I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,MPI_ELOOP,NPROCS,MPI_ID
      USE MOD_DMFT_LDAU,ONLY:SEVT,SEVNST,SEBT,DMFTSIG,DMFTSIGMA,SEBNST,
     &    KSELF
      USE MOD_TAUIJ,ONLY:CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GEN')
      INTEGER NPOLMAX
      PARAMETER (NPOLMAX=3)
C
C Dummy arguments
C
      CHARACTER*3 DOSCORSYS,DOSREP
      LOGICAL RUNELOOP
      CHARACTER*10 TASK
C
C Local variables
C
      REAL*8 ABITNEW(:,:,:,:),AXCN(:,:,:,:),BCOR(:),BCORS(:),BXCNMM(:,:)
     &       ,BXCNMN(:,:),BXCNNM(:,:),BXCNNN(:,:),CLURAD,CMNTTX(:,:),
     &       CPACHNG,CPACHTAB(:),CPAERR,CPU_TIME_LOCAL_END,
     &       CPU_TIME_LOCAL_START,DMFTMIX,EBANDLD,ECOR(:),EFERMI0,
     &       EFERMILD,EMM(:,:),ENM(:,:),ENN(:,:),ERYDTOP,FCOR(:,:,:),
     &       GAMMAM(:,:),GAMMAN(:,:),GCOR(:,:,:),NEFTL(:,:),NLCPAMAD,
     &       OBS_POS_LT(:,:,:,:),OBS_POS_LTX(:,:,:,:),OBS_POS_T(:,:,:),
     &       OBS_POS_TX(:,:,:),ORBSQINT(:,:),RDUMARG,RHO2NSX(:,:,:,:),
     &       RHO2NS_EL(:,:,:,:),RHOCHRX(:,:),RHOCHR_EL(:,:),RHOORBX(:,:)
     &       ,RHOSPNX(:,:),SCFMIX,SCLNOS,SHFTEF,SZCOR(:),TOTDOS,TOTDOSX,
     &       TOTNOS,TOTNOSX
      LOGICAL CALCFMAG,CALCGFUNMAT,CALCINT,CALCJXC,GETIRRSOL,GETTAUIJ,
     &        INTEFERMI,MPI_OUTPUT_FLAG(:),POSANI,POSANIPREP,POSANITAU,
     &        POSITRON,RENORMALIZE,RHOPOS_ON_DISK,WRLOG
      COMPLEX*16 CTOTDOS(:),EREF(:),ERYD,ERYD9,GFMAT(:,:,:,:),
     &           HINT(:,:,:),P
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IE,IECPAFAIL(:),IECURR,IEPANEL,
     &        IEPATH,IERR,IKMCOR(:,:),IL,ILA,IM,IN,INC,IP,IPRINTBAND,
     &        IPROCE(:),IQ,IQ9,IQBOTSAV,IQTOPSAV,IR,IT,ITBOTSAV,ITCPA,
     &        ITEREF,ITTOPSAV,ITXRAY,ITXRAYDUM,IWRIRRWF,IWRREGWF,
     &        IZERO(:),J,KAPCOR(:),LMAXYLM,M,MM05COR(:),N,N1,N2,
     &        NCPAFAIL,NCSTMAX,NE,NKPCOR(:),NLCPAOUT,NN,NWRLOG
      CHARACTER*80 LINE
      CHARACTER*5 STR5
C
C*** End of declarations rewritten by SPAG
C
C
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
      DATA SHFTEF/0.0D0/,TOTNOS/0D0/,TOTDOS/0D0/
      DATA ITXRAYDUM/0/,SCFMIX/0D0/,GETTAUIJ/.FALSE./
      DATA RDUMARG/999999D0/,RENORMALIZE/.FALSE./
      DATA EBANDLD/0D0/
      DATA SCLNOS/0D0/
C ======================================================================
C
      ALLOCATABLE GFMAT,HINT
      ALLOCATABLE CPACHTAB,CTOTDOS,IECPAFAIL,EREF,IPROCE
C
C------------------------------------------ variables depending on NTMAX
C
      ALLOCATABLE BCOR,BCORS,ABITNEW
C
C------------------------------------------- variables depending on NCST
C
      ALLOCATABLE FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR
      ALLOCATABLE NKPCOR
C
      ALLOCATABLE RHO2NSX,RHOCHRX,RHOORBX,RHOSPNX
      ALLOCATABLE RHOCHR_EL,RHO2NS_EL
      ALLOCATABLE CMNTTX
      ALLOCATABLE OBS_POS_LT,OBS_POS_LTX,OBS_POS_T,OBS_POS_TX
C
C-------------------------------------------- variables depending on MPI
C
      ALLOCATABLE MPI_OUTPUT_FLAG
C
C----------------------------------- variables used for STONER parameter
C
      ALLOCATABLE AXCN,BXCNMM,BXCNMN,BXCNNM,BXCNNN,NEFTL,EMM,ENM,ENN
      ALLOCATABLE GAMMAM,GAMMAN,ORBSQINT
C
C=======================================================================
C
C----------------------------------------- variables depending on NMEMAX
C
      ALLOCATE (HINT(NLMAX,NTMAX,NMEHFMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: HINT')
C
C------------------------------------------ variables depending on NEMAX
C
      ALLOCATE (IECPAFAIL(NEMAX),CPACHTAB(NEMAX),IPROCE(NEMAX))
      ALLOCATE (CTOTDOS(NEMAX),EREF(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IECPAFAIL')
C
C------------------------------------------ variables depending on NTMAX
C
      ALLOCATE (BCOR(NTMAX),BCORS(NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: BCOR')
C
C------------------------------------------ variables depending on NLMAX
C
      IF ( ALLOCATED(DOBS_LTEX) ) DEALLOCATE (DOBS_LTEX,DOBS_TEX_GLO)
      ALLOCATE (DOBS_LTEX(0:3,NOBSMAX,NLMAX,NTMAX,NEMAX))
      ALLOCATE (DOBS_TEX_GLO(0:3,NOBSMAX,NTMAX,NEMAX))
C
C-----------------------------------------------------------------------
C         use a constant potential for testing purposes
C-----------------------------------------------------------------------
C
      IF ( USE_CONST_POTENTIAL ) THEN
         VT(:,:) = POTFIX
         BT(:,:) = 0D0
      END IF
C
C-----------------------------------------------------------------------
C  these arrays are needed only if information on a specific core shell
C   has to be obtained from <CORE> e.g. for core level spectroscopies
C-----------------------------------------------------------------------
C
      IF ( TASK(1:5).EQ.'APS  ' .OR. TASK(1:5).EQ.'AES  ' .OR. TASK(1:5)
     &     .EQ.'NRAES' .OR. TASK(1:7).EQ.'MECHECK' .OR. TASK(1:6)
     &     .EQ.'WFPLOT' .OR. TASK(1:6).EQ.'MEPLOT' .OR. 
     &     TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &     TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' .OR. 
     &     TASK.EQ.'CLXPS     ' ) THEN
         NCSTMAX = 14
      ELSE
         NCSTMAX = 1
      END IF
C
      ALLOCATE (IZERO(NCSTMAX),SZCOR(NCSTMAX),ECOR(NCSTMAX))
      ALLOCATE (KAPCOR(NCSTMAX),MM05COR(NCSTMAX))
      ALLOCATE (IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX))
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GCOR')
C
C-----------------------------------------------------------------------
C
      WRLOG = .FALSE.
C     CALCGFUNMAT = .TRUE.
      CALCGFUNMAT = .FALSE.
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( .NOT.FULLPOT ) THEN
         DO IM = 1,NMMAX
            JRNS1(IM) = JRNSMIN
            JRCRI(IM) = JRWS(IM)
         END DO
      END IF
      IF ( DMFT .OR. LDAU ) THEN
C
         CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),DMFTMIX,NEMAX)
C
         IF ( TASK(1:3).EQ.'AES' .OR. TASK(1:5).EQ.'NRAES' .OR. 
     &        TASK(1:5).EQ.'APS  ' .OR. TASK(1:5).EQ.'FSCAT' .OR. 
     &        TASK(1:6).EQ.'SOCPAR' .OR. TASK(1:6).EQ.'PSHIFT' .OR. 
     &        TASK(1:5).EQ.'CLXPS' .OR. TASK(1:7).EQ.'WFPLOT' .OR. 
     &        TASK(1:7).EQ.'MECHECK' .OR. TASK(1:7).EQ.'PLOTPRS' .OR. 
     &        TASK(1:8).EQ.'PLOT2DPR' .OR. TASK(1:2).EQ.'T1' )
     &        CALL STOP_MESSAGE(ROUTINE,
     &             'the routine to be called does work for DMFT or LDAU'
     &             )
      END IF
C=======================================================================
C
C=======================================================================
C                        Green''s function matrix
C=======================================================================
C
      IF ( CALCGFUNMAT ) THEN
         ALLOCATE (GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GFMAT')
      END IF
C
C--------------------------------- some consitency checks for parameters
      IF ( .NOT.RDTAUMQ .AND. .NOT.RDTAU ) THEN
C
      END IF
C
C-----------------------------------------------------------------------
C                            BREIT INTERACTION
C-----------------------------------------------------------------------
      IF ( BREITINT .AND. SEARCHEF ) THEN
C
         ALLOCATE (ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX))
C
         ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0D0
C
         CALL SCFBIMAD
C
      END IF
C
C=======================================================================
C                        Lloyd formula
C=======================================================================
C
      ERYDTOP = 1.4D0*EFERMI
C
      IF ( ITEST.EQ.5 ) CALL LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,MSST,SSST,
     &                                C0,0)
C
C--------------------------------------------------- free electron poles
C
      IF ( PROGNAME.EQ.'KKRGEN    ' ) THEN
         IF ( LLOYD .OR. IBZINT.EQ.4 ) THEN
C
            ERYDTOP = 1.4D0*EFERMI
            ERYDTOP = 1.4D0*EMAX
C
            CALL FRELPOLE(ERYDTOP,IPRINT)
         END IF
      END IF
C
C-----------------------------------------------------------------------
C
      CALL CINIT(NEMAX,EREF)
C
      IWRREGWF = 0
      IWRIRRWF = 0
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
C     IF( SPLITSS ) THEN
C        GETIRRSOL = .FALSE.
C     END IF
C
C----------- NO calculation of matrix elements, irregular wave function,
C--------------------------- and NO output of wave functions in <SSITE>,
C------------------ nor of TAU-matrix in <PROJTAU> for CLUSTER ELOOP-run
      IF ( IBZINT.EQ.0 ) THEN
         IWRREGWF = 0
         IWRIRRWF = 0
C        CALCINT   = .FALSE.
C        GETIRRSOL = .FALSE.
C
         CALCINT = .TRUE.
         GETIRRSOL = .TRUE.
C
         WRTAU = .FALSE.
         WRTAUMQ = .FALSE.
      END IF
C
C-------------------- allow tasks that need ELOOP to be done in parallel
C
      IF ( TASK(1:2).EQ.'TL' ) THEN
         CALCFMAG = TASK(3:3).EQ.'1'
         CALCJXC = TASK(4:4).EQ.'1'
         IWRREGWF = 1
         IWRIRRWF = 1
         CALCINT = .TRUE.
         GETIRRSOL = .TRUE.
         ICORE = 1
         NEPATH = 1
      ELSE
         CALCFMAG = .FALSE.
         CALCJXC = .FALSE.
      END IF
C
      IF ( (CALCJXC .AND. IREL.EQ.3) .AND. GF_CONV_RH ) THEN
         GF_CONV_RH = .FALSE.
         CALL INFO_MESSAGE(ROUTINE,
     &         'GF_CONV_RH set .FALSE. for calculation of J_ij - tensor'
     &         )
      END IF
C
      IF ( CALCFMAG .OR. ISMQHFI.EQ.1 ) THEN
C
         OBS_T_GLO(:,:,:) = 0D0
         OBS_T(:,:,:) = 0D0
         OBS_LT(:,:,:,:) = 0D0
C
      END IF
C
C--------------------- initialize arrays for current density calculation
C
      IF ( CALC_JALF ) THEN
C
         CALL AMEBI
C
         ALLOCATE (JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX))
         JDNST(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0D0
C
      END IF
C
      IF ( TASK(1:5).EQ.'CLXPS' .AND. .NOT.RDTAU ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         LMAXYLM = 2*(NLMAX-1)
C
         CALL CLXPS(1,CALCINT,GETIRRSOL,IPRINT,TAUT,TSST,MSST,MEZZ,MEZJ,
     &              GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,IKMCOR,
     &              IZERO,ITXRAY,NCSTMAX,NPOLMAX,LMAXYLM)
         RUNELOOP = .TRUE.
      END IF
C
C=======================================================================
C               prepare calculation of site off-diagonal
C=======================================================================
C
C
      IF ( KKRMODE(1:16).EQ.'EMBEDDED-CLUSTER' .AND. 
     &     CLUTYPE.EQ.'embedded' .AND. CALCJXC )
     &     CALC_OFF_DIAGONAL_TAUIJ_EMBEDDED = .TRUE.
C
      IF ( ((CALCJXC .AND. KKRMODE(1:16).NE.'EMBEDDED-CLUSTER') .OR. 
     &     TASK(1:8).EQ.'NUCSSCPL' .OR. TASK(1:8).EQ.'FORCETEN' .OR. 
     &     (KKRMODE(1:16).NE.'EMBEDDED-CLUSTER' .AND. 
     &     CLUTYPE.EQ.'embedded')) .AND. IBZINT.NE.0 ) THEN
C
         GETTAUIJ = .TRUE.
C
C-----------------------------------------------------------------------
         IF ( KKRMODE(1:16).EQ.'EMBEDDED-CLUSTER' .OR. 
     &        CLUTYPE.NE.'embedded' ) THEN
C
            CALL INPUT_FIND_SECTION('TASK',0)
C
            CALL SECTION_SET_REAL('CLURAD',CLURAD,1.5D0,0)
C
            CALL INIT_MOD_TAUIJ_STAR(CLURAD)
C
         END IF
C
C
         IF ( IBZINT.EQ.2 ) THEN
C
            CALL INIT_MOD_TAUIJ_KMESH
C
         ELSE
            WRITE (6,*) 'IBZINT should be 2  but not',IBZINT
            WRITE (6,*) 'for calculation of TAU(i,j) and TASK = JXC'
            STOP
         END IF
C
      ELSE
         GETTAUIJ = .FALSE.
      END IF
C
      ICORE = 1
C
C=======================================================================
C                        deal with positron states
C                   note hyrarchy of tasks for input
C                  POSANI > POSANIPREP > POSANITAU >  POSITRON
C=======================================================================
C     POSANIPREP   prepare calculation of positron anihilation spactra
C                  i.a. calculate posit. structure for posit. potential
C     POSANI       calculation of positron anihilation spactra
C                  i.a. use electronic potential >> no call of POSANIPOT
C     POSANITAU    calculate positronic lifetime  for posit. potential
C     POSITRON     calculate DOS, etc. for posit. potential
C=======================================================================
C
      POSANI = .FALSE.
      POSANIPREP = .FALSE.
      POSANITAU = .FALSE.
      POSITRON = .FALSE.
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('POSANI',POSANI)
         IF ( .NOT.POSANI ) THEN
            CALL SECTION_FIND_KEYWORD('POSANIPREP',POSANIPREP)
            IF ( .NOT.POSANIPREP ) THEN
               CALL SECTION_FIND_KEYWORD('POSANITAU',POSANITAU)
               IF ( .NOT.POSANITAU ) THEN
                  CALL INPUT_FIND_SECTION('MODE',0)
                  IF ( FOUND_SECTION )
     &                  CALL SECTION_FIND_KEYWORD('POSITRON',POSITRON)
               END IF
            END IF
         END IF
      END IF
C
      IF ( POSITRON .OR. POSANIPREP .OR. POSANITAU ) THEN
C Allocate variables needed for Positron calculations
         ALLOCATE (OBS_POS_T(0:3,NOBSMAX,NTMAX))
         ALLOCATE (OBS_POS_TX(0:3,NOBSMAX,NTMAX))
         ALLOCATE (OBS_POS_LT(0:3,NOBSMAX,NLMAX,NTMAX))
         ALLOCATE (OBS_POS_LTX(0:3,NOBSMAX,NLMAX,NTMAX))
         ALLOCATE (CMNTTX(NLMFPMAX,NTMAX))
         OBS_POS_T(:,:,:) = 0D0
         OBS_POS_TX(:,:,:) = 0D0
         OBS_POS_LT(:,:,:,:) = 0D0
         OBS_POS_LTX(:,:,:,:) = 0D0
         CMNTTX = 0D0
C
         UPDATE_EFERMI = .FALSE.
         SCLNOS = 1D0
         ICORE = 0
         NLQMAD = 1
C
         CALL POSANIPOT
C
         IF ( POSANIPREP .OR. POSANITAU ) THEN
            IWRREGWF = 1
            IWRIRRWF = 1
            CALCINT = .TRUE.
            GETIRRSOL = .TRUE.
         END IF
C
         IF ( POSANITAU ) THEN
            IF ( IGRID(1).NE.5 )
     &            CALL STOP_MESSAGE(ROUTINE,'POSANITAU .AND. IGRID.NE.5'
     &           )
C
            IF ( FULLPOT ) THEN
               ALLOCATE (RHO2NS_EL(NRMAX,NLMFPMAX,NTMAX,3))
               ALLOCATE (RHOCHR_EL(1,1))
               RHO2NS_EL(1:NRMAX,1:NLMFPMAX,1:NTMAX,1:3)
     &            = RHO2NS(1:NRMAX,1:NLMFPMAX,1:NTMAX,1:3)
            ELSE
               ALLOCATE (RHO2NS_EL(1,1,1,1))
               ALLOCATE (RHOCHR_EL(NRMAX,NTMAX))
               RHOCHR_EL(1:NRMAX,1:NTMAX) = RHOCHR(1:NRMAX,1:NTMAX)
            END IF
C
            OBS_T_GLO(:,:,:) = 0D0
            OBS_T(:,:,:) = 0D0
            OBS_LT(:,:,:,:) = 0D0
C
            IF ( IREL.LE.1 ) THEN
               NVALTOT = 1
            ELSE
               NVALTOT = 2
            END IF
C
         END IF
C
         ALLOCATE (RHOORBX(NRMAX,NTMAX))
         ALLOCATE (RHOCHRX(NRMAX,NTMAX),RHOSPNX(NRMAX,NTMAX))
         RHOCHRX(:,:) = 0.0D0
         RHOSPNX(:,:) = 0.0D0
         RHOORBX(:,:) = 0.0D0
C
         IF ( FULLPOT ) THEN
            ALLOCATE (RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3))
         ELSE
C                            rho2nsx needed as dummy in ASA calculations
            ALLOCATE (RHO2NSX(1,1,1,1))
            RHO2NSX = 0D0
         END IF
C
         CALL POSANIRHOIO('READ ',RHOPOS_ON_DISK)
C
         IF ( RHOPOS_ON_DISK ) THEN
C
            RUNELOOP = .FALSE.
C
         ELSE IF ( .NOT.GETTAUIJ ) THEN
C
            RHOCHR(:,:) = 0.0D0
            RHOSPN(:,:) = 0.0D0
            RHOORB(:,:) = 0.0D0
C
            IF ( FULLPOT ) THEN
               RHO2NS(:,:,:,:) = 0.0D0
               RHO2NSX(:,:,:,:) = 0.0D0
            END IF
C
         END IF
C
      END IF
C
      IF ( POSANI ) TASK = 'POSANI    '
C
C=======================================================================
      IF ( RUNELOOP ) THEN
C=======================================================================
C
         IF ( ICORE.NE.0 ) CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,
     &                               MM05COR,NKPCOR,IKMCOR,IZERO,
     &                               ITXRAYDUM,BCOR,BCORS,NCSTMAX)
C
         IF ( ABS(EMAX-EFERMI).LT.1D-3 ) THEN
            INTEFERMI = .TRUE.
         ELSE
            INTEFERMI = .FALSE.
         END IF
         IF ( TASK(1:5).EQ.'APS  ' .OR. TASK(1:5).EQ.'AES  ' .OR. 
     &        TASK(1:5).EQ.'NRAES' .OR. TASK(1:7).EQ.'MECHECK' .OR. 
     &        TASK(1:6).EQ.'WFPLOT' .OR. TASK(1:6).EQ.'MEPLOT' .OR. 
     &        TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &        TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' .OR. 
     &        TASK.EQ.'CLXPS     ' .OR. TASK.EQ.'VBPES     ' .OR. 
     &        TASK.EQ.'ARPES     ' .OR. TASK.EQ.'COMPTON   ' .OR. 
     &        TASK.EQ.'POSANI    ' ) INTEFERMI = .FALSE.
C
C ======================================================================
C                   mode for matrix element calculation
C
C   apply Complex Conjugation to radial wave function of <BRA| state
C   for matrix elements specific to spectroscopy
C ======================================================================
C
         IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &        TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' .OR. 
     &        TASK.EQ.'CLXPS     ' .OR. TASK.EQ.'VBPES     ' .OR. 
     &        TASK.EQ.'ARPES     ' ) THEN
C
            ME_CC_BRA_RWF = .TRUE.
C
         ELSE
C
            ME_CC_BRA_RWF = .FALSE.
C
         END IF
C
C
         ITEREF = 0
C=======================================================================
 50      CONTINUE
         EFERMI0 = EFERMI
         ITEREF = ITEREF + 1
C
         IF ( WRTAU .OR. (IBZINT.EQ.0) ) THEN
            REWIND 9
            WRITE (9,99002) NT,NQ
            WRITE (9,'(5X,'' IT ='',I3,'' IQ ='',I3,'': '',A)')
     &             (IT,IQAT(1,IT),TXT_T(IT),IT=1,NT)
            DO IQ = 1,NQ
               WRITE (9,'(''site'',I3,'' of '',A)') IQ,SYSTEM(1:LSYSTEM)
            END DO
         ELSE IF ( WRTAUMQ ) THEN
            REWIND 9
            WRITE (9,99003) NQ
            DO IQ = 1,NQ
               WRITE (9,'(''site'',I3,'' of '',A)') IQ,SYSTEM(1:LSYSTEM)
            END DO
         END IF
C
         NCPAFAIL = 0
C
         CALL CPU_TIME(CPU_TIME_LOCAL_START)
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
         IPROCE(1:NETAB(1)) = 0
C
         IF ( MPI ) THEN
C
            MPI_ELOOP = .TRUE.
            MPI_KLOOP = .FALSE.
C
            ALLOCATE (MPI_OUTPUT_FLAG(0:(NPROCS-1)))
C
            N1 = NETAB(1)/NPROCS
            N2 = MOD(NETAB(1),NPROCS)
C
            IE = NETAB(1) + 1
            DO IP = 1,NPROCS
               IF ( IP.LE.N2 ) THEN
                  NN = N1 + 1
               ELSE
                  NN = N1
               END IF
C
               DO INC = 1,NN
                  IE = IE - 1
                  IPROCE(IE) = IP - 1
               END DO
               MPI_OUTPUT_FLAG(IP-1) = .TRUE.
            END DO
C
            IF ( WRTAU .OR. WRTAUMQ )
     &           CALL MPI_MERGE_FILES(IFILTAU,TAUFIL,MPI_OUTPUT_FLAG,0)
C
         ELSE
C
            ALLOCATE (MPI_OUTPUT_FLAG(0:0))
            MPI_OUTPUT_FLAG(0) = .FALSE.
C
            MPI_ELOOP = .FALSE.
            MPI_KLOOP = .FALSE.
C
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C=======================================================================
C===============         LOOP OVER TEMPERATURE          ========== START
C=======================================================================
         TEMP_LAT = 0D0
         N_TEMP_LAT = 1
C
         CALL THERMAL_INIT(0,N_TEMP_LAT,TEMP_LAT)
C
         IF ( THERMAL_VIBRA_FLUCT ) THEN
            IWRREGWF = 1
            IWRIRRWF = 1
            CALCINT = .TRUE.
            GETIRRSOL = .TRUE.
            ICORE = 1
            NEPATH = 1
         END IF
C
         WAVE_FUNCTIONS_AVAILABLE = GETIRRSOL .AND. (IWRREGWF.EQ.1)
     &                              .AND. (IWRIRRWF.EQ.1)
C
C-----------------------------------------------------------------------
C
         I_TEMP_LAT = 0
 100     CONTINUE
         I_TEMP_LAT = I_TEMP_LAT + 1
C
C-----------------------------------------------------------------------
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C            and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ
C
         CALL THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
C=======================================================================
C
         IEPANEL = 1
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         LOOP_EPATH:DO IEPATH = 1,NEPATH
C
            TOTNOS = 1D0
            CALL CINIT(NETAB(IEPATH),CTOTDOS)
C
            IF ( DMFT .OR. LDAU )
     &           CALL DMFT_READSIG(NETAB(IEPATH),ETAB(1,IEPATH),IPRINT)
C
C--------------------------------- moments, hyperfine field, band energy
C
            OBS_TX_GLO(:,:,:) = 0D0
            OBS_TX(:,:,:) = 0D0
            OBS_LTX(:,:,:,:) = 0D0
C
            DOBS_LTEX(:,:,:,:,:) = 0D0
            DOBS_TEX_GLO(:,:,:,:) = 0D0
C
            LOOP_E:DO IE = 1,NETAB(IEPATH)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
               IF ( MPI_ELOOP .AND. MPI_ID.NE.IPROCE(IE) ) CYCLE LOOP_E
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
               IECURR = IE
               ICPAFLAG = 0
               CPACHNG = 0.0D0
C
C======= read site-diagonal TAU- and inverse eff. single site t-matrix m
               IF ( RDTAUMQ .AND. IEPATH.EQ.1 ) THEN
C
                  CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C
                  DO IQ = 1,NQ
C
 102                 CONTINUE
                     READ (9,'(A5)',END=200) STR5
                     IF ( STR5.NE.'*****' ) GOTO 102
C
                     READ (9,'(2F21.15,25X,I2,9X,I2,F15.6)',END=200)
     &                     ERYD9,IQ9,ICPAFLAG,CPACHNG
                     IF ( IQ.NE.IQ9 ) CALL STOP_MESSAGE(ROUTINE,
     &                    'reading TAU:  IQ<>IQ9')
C
 104                 CONTINUE
                     READ (9,'(2I5,1P,4E22.14)') I,J,TAUQ(I,J,IQ),
     &                     MSSQ(I,J,IQ)
C
                     IF ( (I+J).LT.2*NKMQ(IQ) ) GOTO 104
C
                  END DO
C
                  ERYD = ERYD9
                  ETAB(IE,IEPATH) = ERYD9
C
               ELSE
C
                  ERYD = ETAB(IE,IEPATH)
C
               END IF
C
C======================================solve SS - differential equation
C
               IF ( ORBPOL(1:6).EQ.'SIGMA ' .OR. DMFT ) THEN
                  DO IT = 1,NTMAX
                     IM = IMT(IT)
                     DO IR = 1,JRCRI(IM)
                        SEVT(IR,IT) = C0
                        SEBT(IR,IT) = C0
                     END DO
                     DO IR = JRNS1(IM),JRCRI(IM)
                        DO IN = 1,NLMFPMAX
                           SEVNST(IR,IN,IT) = C0
                           SEBNST(IR,IN,IT) = C0
                        END DO
                     END DO
                  END DO
               END IF
C
               IF ( ORBPOL(1:6).EQ.'SIGMA ' ) THEN
                  CALL DMFT_READSELFENE(ERYD,EFERMI,KSELF,SEVT,SEBT,
     &                                  SEVNST,SEBNST,IMT,JRNS1,JRCRI,
     &                                  JRNSMIN,NTMAX,NLMFPMAX,NRMAX,
     &                                  FULLPOT)
                  IF ( IE.EQ.1 ) THEN
                     WRITE (6,'(1X,79(''*''),/)')
                     DO IT = 1,NT
                        IF ( KSELF(IT).EQ.1 ) WRITE (6,
     &                       '(19X,''SELFENERGY applied for IT:'',I3)')
     &                       IT
                     END DO
                     WRITE (6,'(1X,79(''*''),/)')
                  END IF
               END IF
C
               IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &              = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
               IWRREGWF = 1
               IF ( FULLPOT ) IWRIRRWF = 1
C
C
               IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6)
     &              .EQ.'I-ZONE' ) THEN
                  ITBOTSAV = ITBOT
                  ITTOPSAV = ITTOP
                  IQBOTSAV = IQBOT
                  IQTOPSAV = IQTOP
                  ITBOT = 1
                  ITTOP = NT
                  IQBOT = 1
                  IQTOP = NQ
               END IF
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,
     &                       GETIRRSOL,ERYD,P,IPRINT,TSST,MSST,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
               IF ( THERMAL_VIBRA_FLUCT ) CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
               IF ( .NOT.RDTAUMQ .AND. IEPATH.EQ.1 )
     &              CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
               IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6)
     &              .EQ.'I-ZONE' ) THEN
                  ITBOT = ITBOTSAV
                  ITTOP = ITTOPSAV
                  IQBOT = IQBOTSAV
                  IQTOP = IQTOPSAV
               END IF
C
C ======================================================================
C                            CALCULATE TAU
C ======================================================================
C
               IF ( ITEST.EQ.2 ) THEN
C
                  IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
                     CALL SYMCHECK_STD(ERYD,P,ICPAFLAG,CPACHNG,IPRINT,
     &                                 ITCPA,ICPACONV,TSST,MSST,TSSQ,
     &                                 MSSQ,TAUQ)
C
                     CALL SYMCHECKIJ_STD(ERYD,P,ICPAFLAG,CPACHNG,IPRINT,
     &                  ITCPA,ICPACONV,TSST,MSST,TSSQ,MSSQ,TAUQ)
C
                  ELSE
C
                     CALL SYMCHECK_TB(ERYD,P,ICPAFLAG,CPACHNG,IPRINT,
     &                                ITCPA,ICPACONV,TSST,MSST,TSSQ,
     &                                MSSQ,TAUQ)
C
C                  CALL SYMCHECKIJ(MSSQ,ERYD,P,TAUQ,ICPAFLAG,CPACHNG,
C     &                            IPRINT,ITCPA,ICPACONV,MSST)
C
                  END IF
C
               END IF
C
C ----------------------------------------------------------- IEPATH = 1
               IF ( IEPATH.EQ.1 ) THEN
C
                  IF ( .NOT.RDTAUMQ ) THEN
C
                     IPRINTBAND = IPRINT
C
                     CALL TAU_DRIVE(IECURR,IPRINTBAND,ERYD,P,TSSQ,MSSQ,
     &                              TSST,MSST,TAUQ,ICPAFLAG,CPACHNG,
     &                              ITCPA,ICPACONV,PHASK)
C
                     IF ( ICPAFLAG.NE.0 ) THEN
                        NCPAFAIL = NCPAFAIL + 1
                        CPACHTAB(NCPAFAIL) = CPACHNG
                        IECPAFAIL(NCPAFAIL) = IECURR
                     END IF
C
                  END IF
C
C ----------------------------------------------- project TAU on type IT
C
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
C
C ---------------------------------------------------------------- NLCPA
                  IF ( USENLCPA ) THEN
                     NLCPAOUT = 0
                     IF ( NEPATH.NE.1 )
     &                     CALL STOP_MESSAGE(ROUTINE,'NEPATH <> 1')
                     IF ( NLCPAWRCFG ) THEN
                        IF ( IGRID(1).EQ.3 ) THEN
                           NLCPAOUT = 2
                        ELSE
                           NLCPAOUT = 3
                        END IF
                     END IF
C
                     CALL NLCPA(NLCPAOUT,NLCPAMAD,ICPAFLAG,CPAERR,
     &                          IPRINT,WRMAT,WRTAUMQ,IFILTAU,TAUQ,ITCPA,
     &                          ICPACONV,MSSQ,MSST,TAUT,IECURR,ERYD,P,
     &                          MEZZ,MEZJ,NKPTS0)
                  END IF
C ---------------------------------------------------------------- NLCPA
C
               END IF
C ----------------------------------------------------------- IEPATH = 1
C
C ======================================================================
C                     deal with Lloyd's formula
C ======================================================================
C
               IF ( LLOYD ) CALL LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,MSST,
     &              SSST,ERYD,IE)
C
C ======================================================================
C  if the epath is split into BS- and SS-part manipulate TAU accordingly
C
               IF ( SPLITSS ) THEN
                  M = NKMMAX
                  IF ( IEPATH.EQ.1 ) THEN
                     TAUT(1:M,1:M,ITBOT:ITTOP)
     &                  = TAUT(1:M,1:M,ITBOT:ITTOP)
     &                  - TSST(1:M,1:M,ITBOT:ITTOP)
                  ELSE
                     TAUT(1:M,1:M,ITBOT:ITTOP)
     &                  = TSST(1:M,1:M,ITBOT:ITTOP)
                  END IF
               END IF
C
C ======================================================================
C
               IF ( WRLOG .AND. IE.EQ.NETAB(NEPATH) ) THEN
                  LINE = DATSET(1:LDATSET)//'.log'
                  OPEN (UNIT=IFILLOG,FILE=LINE)
                  NWRLOG = 2
               ELSE
                  NWRLOG = 1
               END IF
C
               IF ( (MPI_ID.EQ.0 .OR. .NOT.MPI_KLOOP) .AND. 
     &              .NOT.CALCJXC )
     &              CALL CALCDOS(NWRLOG,INTEFERMI,SPLITSS,IEPATH,
     &              NCPAFAIL,IPRINT,ERYD,MEZZ,MEZJ,TSST,MSST,TAUT,MSSQ,
     &              TAUQ,IECURR,WETAB(IECURR,IEPATH),BCOR,BCORS,
     &              DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX,DOBS_TX_GLO,
     &              OBS_TX_GLO,CTOTDOS)
C
C ----------------------------------------------------------------------
C      calculate charge densities for positron liftime calculations
C ----------------------------------------------------------------------
C                  only implemented for IEPANEL = 1
C
               IF ( POSANITAU ) THEN
C
                  CALL CHRDNS(.FALSE.,CTOTDOS,SHFTEF,TOTDOSX,TOTDOS,
     &                        TOTNOSX,ERYD,WETAB(IECURR,IEPATH),RDUMARG,
     &                        IEPANEL,IEPATH,IECURR,OBS_POS_LTX,
     &                        OBS_POS_TX,BCOR,BCORS,MEZZ,MEZJ,TSST,MSST,
     &                        TAUT,MSSQ,TAUQ,RHOCHRX,RHOSPNX,RHOORBX,
     &                        RHO2NSX,CMNTTX)
C
                  CALL CHRDNS_NORM(.FALSE.,SCLNOS,TOTDOS,TOTNOS,TOTNOSX,
     &                             IEPANEL,IEPATH,IECURR,OBS_POS_LT,
     &                             OBS_POS_T,OBS_POS_LTX,OBS_POS_TX,
     &                             RHOCHRX,RHOSPNX,RHOORBX,RHO2NSX,
     &                             CMNTTX)
C
               END IF
C
C
C=======================================================================
C                     Green's function matrix
C=======================================================================
C
               IF ( CALCGFUNMAT ) CALL GFUNMAT_DRIVE(ERYD,IE,TAUT,MSST,
     &              GFMAT)
C
C=======================================================================
C                     calculate l,m_l-resolved DOS
C=======================================================================
C
               IF ( DOSREP.NE.'XXX' )
     &              CALL WRDOSREP(DOSREP,DOSCORSYS,MEZZ,MEZJ,TAUT,
     &              IECURR)
C
C
C=======================================================================
C                 calculate site off-diagonal TAU_ij
C=======================================================================
C
C
C
               IF ( GETTAUIJ ) CALL TAUIJ_DRIVE(IECURR,ERYD,P,TSSQ,MSSQ,
     &              TSST,MSST,TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV)
C
C
C=======================================================================
C                   calculate XC-coupling constants
C=======================================================================
C
               IF ( (MPI_ID.EQ.0 .OR. .NOT.MPI_KLOOP) .AND. CALCJXC )
     &              THEN
C
                  IF ( IREL.EQ.2 ) THEN
C
C-----------------------------------------------------------------------
C               site diagonal XC-coupling constant  J_0
C-----------------------------------------------------------------------
C
                     CALL XCPLJ0(TAUT,MSST,IECURR)
C
C-----------------------------------------------------------------------
C               site off-diagonal XC-coupling constants  J_ij
C-----------------------------------------------------------------------
C
                     IF ( IBZINT.NE.0 ) THEN
C
                        CALL XCPLJIJ(IECURR,TAUQ,MSSQ,MSST)
C
                     ELSE
C
                        CALL CLUXCPLJIJ(IECURR,MSSQ,MSST)
C
                     END IF
C
                  ELSE IF ( IREL.EQ.3 ) THEN
C
                     IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
                        CALL XCPLTENSOR(IECURR,TAUQ,MSSQ,MSST)
C
                     ELSE
C
                        CALL CLUXCPLTENSOR(IECURR,TAUQ,MSSQ,MSST)
C
                     END IF
C
                  END IF
C
Ccc                  CALL STOP_REGULAR(ROUTINE,
Ccc     &                         'XC-coupling constants calculated')
C
               END IF
C
C=======================================================================
C                       analyze bonding
C=======================================================================
               IF ( TASK.EQ.'BONDING   ' ) THEN
C
                  IF ( IBZINT.NE.0 ) THEN
C
                     CALL STOP_MESSAGE(ROUTINE,
     &                             'BOND ANALYSIS only for cluster mode'
     &                             )
C
                  ELSE
C
                     CALL CLUBONDING(IECURR,ERYD,MEZZ,MEZJ,MSSQ,MSST,
     &                               TSST)
C
                  END IF
C
               END IF
C
C=======================================================================
C           calculate nuclear spin-spin coupling constants
C=======================================================================
               IF ( TASK(1:8).EQ.'NUCSSCPL' ) THEN
                  IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &                 'TASK does not YET work for FULL POTENTIAL mode')
C
                  CALL NUCSSCPL(IECURR,MEZZ,TAUQ,MSSQ,MSST)
C
               END IF
C
C=======================================================================
C                          calculate force tensor
C=======================================================================
C
               IF ( TASK(1:8).EQ.'FORCETEN' )
     &              CALL FORCETENSOR(IECURR,P,TAUQ,MSSQ,MSST)
C
C=======================================================================
C           determine decomposition of relativistic hyperfine field
C=======================================================================
C
               IF ( (ISMQHFI.EQ.1) .AND. (IREL.GT.1) ) THEN
                  IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &                 'TASK does not YET work for FULL POTENTIAL mode')
C
                  CALL HFF_SPLIT(NWRLOG,IEPATH,EFERMI-EFERMI0,IPRINT,
     &                           MSST,TAUT,IECURR,OBS_LTX,HINT)
C
               END IF
C
               IF ( CALCFMAG ) THEN
C
                  IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &                 'TASK does not YET work for FULL POTENTIAL mode')
C
                  CALL DENSITY(IECURR,TAUT,OBS_TX)
C
               END IF
C
C=======================================================================
C            calculate current density and magnetic moment
C=======================================================================
C
               IF ( CALC_JALF ) CALL CURDNS(.FALSE.,RENORMALIZE,SHFTEF,
     &              SCLNOS,JDNST,WETAB(IE,IEPATH),TAUT)
C
C=======================================================================
C            calculate orbital current density and magnetic moment
C=======================================================================
C
               IF ( CALC_JORB ) CALL ORBCURDNS(.FALSE.,RENORMALIZE,
     &              SHFTEF,SCLNOS,WETAB(IE,IEPATH),MSST,TAUT)
C
C=======================================================================
               IF ( ITEST.EQ.6 ) THEN
                  N = NKM
                  M = NKMMAX
                  DO IT = 1,NT
                     WRITE (6,99001) 'atom type   ',IT
C
                     CALL CHANGEREP(N,M,TSST(1,1,IT),'REL>RLM',WKM1)
C
                     CALL CHANGEREP(N,M,MSST(1,1,IT),'REL>RLM',WKM1)
C
                     CALL CHANGEREP(N,M,TAUT(1,1,IT),'REL>RLM',WKM1)
C
                     CALL CHANGEREP(N,M,MEZZ(1,1,IT,1),'REL>RLM',WKM1)
C
                     CALL CHANGEREP(N,M,MEZJ(1,1,IT,1),'REL>RLM',WKM1)
C
                  END DO
C
                  DO IQ = 1,NQ
                     WRITE (6,99001) 'lattice site ',IQ
C
                     CALL CHANGEREP(N,M,MSSQ(1,1,IQ),'REL>RLM',WKM1)
C
                     CALL CHANGEREP(N,M,TAUQ(1,1,IQ),'REL>RLM',WKM1)
C
                  END DO
               END IF
C
               IF ( (IPRINT.GE.3 .OR. WRMAT) .AND. IEPATH.EQ.1 )
     &              CALL DUMPMAT(IE,ERYD,6,MSST,MSSQ,TAUT,TAUQ,.TRUE.,
     &                           MEZZ,MEZJ)
C
            END DO LOOP_E
C
C ======================================================================
C                              write DOS
C ======================================================================
C
            CALL CALCDOS_WRITE(IEPATH)
C
         END DO LOOP_EPATH
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ELOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IECURR = NETAB(1) + 1
C
C=======================================================================
C                   calculate XC-coupling constants
C=======================================================================
C
            IF ( CALCJXC ) THEN
C
               IF ( IREL.EQ.2 ) THEN
C
C-----------------------------------------------------------------------
C               site diagonal XC-coupling constant  J_0
C-----------------------------------------------------------------------
C
                  CALL XCPLJ0(TAUT,MSST,IECURR)
C
C-----------------------------------------------------------------------
C               site off-diagonal XC-coupling constants  J_ij
C-----------------------------------------------------------------------
C
                  IF ( IBZINT.NE.0 ) THEN
C
                     CALL XCPLJIJ(IECURR,TAUQ,MSSQ,MSST)
C
                  ELSE
C
                     CALL CLUXCPLJIJ(IECURR,MSSQ,MSST)
C
                  END IF
C
               ELSE IF ( IREL.EQ.3 ) THEN
C
                  IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
                     CALL XCPLTENSOR(IECURR,TAUQ,MSSQ,MSST)
C
                  ELSE
C
                     CALL CLUXCPLTENSOR(IECURR,TAUQ,MSSQ,MSST)
C
                  END IF
C
               END IF
C
            END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C=======================================================================
C
         IF ( THERMAL_VIBRA_FLUCT .AND. (I_TEMP_LAT.LT.N_TEMP_LAT) )
     &        GOTO 100
C
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== END =
C=======================================================================
C=======================================================================
C
         IF ( WRTAU .OR. WRTAUMQ ) CALL MPI_MERGE_FILES(IFILTAU,TAUFIL,
     &        MPI_OUTPUT_FLAG,1)
C
         IF ( LLOYD ) THEN
C
            CALL LLOYDPT2(TOTNOS,CTOTDOS,SCLNOS,EBANDLD,EFERMILD,SCFMIX,
     &                    NCPAFAIL,IECPAFAIL,PHASK,PHAST,PHASA)
C
            IF ( WRNOS ) CALL CALCNOS(CTOTDOS,PHASK,PHAST,PHASA)
C
         END IF
C
         CALL CPU_TIME(CPU_TIME_LOCAL_END)
C
         WRITE (6,99004) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                   CPU_TIME_LOCAL_END - CPU_TIME_LOCAL_START
C
         IF ( BREITINT .AND. SEARCHEF ) THEN
C
            DO IT = 1,NT
               OBS_T(0,ISMT,IT) = 0.0D0
               OBS_T(0,IOMT,IT) = 0.0D0
               DO IL = 1,NLT(IT)
                  OBS_T(0,ISMT,IT) = OBS_T(0,ISMT,IT)
     &                               + OBS_LT(0,ISMT,IL,IT)
                  OBS_T(0,IOMT,IT) = OBS_T(0,IOMT,IT)
     &                               + OBS_LT(0,IOMT,IL,IT)
               END DO
            END DO
C
            CALL SCFBIPOTMAD(ABITNEW)
C
            DO IT = 1,NT
               DO M = -1, + 1
                  DO ILA = 1,NLABIMAX
                     DO I = 1,NRMAX
                        ABIT(I,ILA,M,IT) = ABITNEW(I,ILA,M,IT)
                     END DO
                  END DO
               END DO
            END DO
C
         END IF
C
         IF ( IEPATH.EQ.1 ) THEN
            IF ( NCPAFAIL.NE.0 ) THEN
               WRITE (6,99005) CPATOL,NCPAFAIL,
     &                         (IECPAFAIL(I),DREAL(ETAB(IECPAFAIL(I),1))
     &                         ,CPACHTAB(I),I=1,NCPAFAIL)
               WRITE (6,'(1X,79(''*''),/)')
            ELSE IF ( NCPA.NE.0 ) THEN
               WRITE (6,99006)
            END IF
         END IF
C
         IF ( SEARCHEF .AND. (ABS(EFERMI-EFERMI0).GT.1.D-8) ) THEN
            IF ( ITEREF.LT.40 ) THEN
               WRITE (6,99007)
               DO IEPATH = 1,NEPATH
                  CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,
     &                       NETAB(IEPATH),ETAB(1,IEPATH),
     &                       WETAB(1,IEPATH),EILOW,IPRINT,NEMAX)
               END DO
               GOTO 50
            ELSE
               WRITE (6,99008)
            END IF
         END IF
C
C                    search E_Fermi with 1mRy precission -- if requested
C=======================================================================
C
C=======================================================================
      END IF
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      IF ( MPI ) THEN
         IF ( TASK.NE.'BSF       ' .AND. TASK.NE.'COMPTON   ' .AND. 
     &        TASK.NE.'POSANI    ' .AND. .NOT.CALCJXC ) THEN
            WRITE (6,99009)
            CALL MPI_FINALIZE(IERR)
            MPI = .FALSE.
         END IF
      END IF
C
C
C
C=======================================================================
C                         execute special TASK
C=======================================================================
C
C
C================================================= XAS - XMO - XES - XRS
      IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &     TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
C
         CALL XRAYSPEC(CALCINT,GETIRRSOL,EREF,GCOR,FCOR,ECOR,SZCOR,
     &                 KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,
     &                 BCORS,NCSTMAX,NPOLMAX)
C
C=================================================================== AES
      ELSE IF ( TASK(1:3).EQ.'AES' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL AES(IPRINT,TSST,MSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,
     &            KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,BCORS,
     &            NCSTMAX)
C
C================================================================= NRAES
      ELSE IF ( TASK(1:5).EQ.'NRAES' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL NRAES(TSST,MSST,SSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,
     &              KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,
     &              BCORS,NCSTMAX)
C
C=================================================================== APS
      ELSE IF ( TASK(1:5).EQ.'APS  ' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL NRAPS(TSST,MSST,SSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,
     &              KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,
     &              BCORS,NCSTMAX)
C
C================================================================= FSCAT
      ELSE IF ( TASK(1:5).EQ.'FSCAT' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL FSCAT(IPRINT,TSST,MSST)
C
C================================================================ SOCPAR
      ELSE IF ( TASK(1:6).EQ.'SOCPAR' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL SOCPAR(IPRINT,TSST,MSST,MEZZ,MEZJ)
C
C================================================================ PSHIFT
      ELSE IF ( TASK(1:6).EQ.'PSHIFT' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         CALL PSHIFT(IPRINT,TSST,MSST,SSST,MEZZ,MEZJ)
C
C================================================================= CLXPS
      ELSE IF ( TASK(1:5).EQ.'CLXPS' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         LMAXYLM = 2*(NLMAX-1)
C
         CALL CLXPS(2,CALCINT,GETIRRSOL,IPRINT,TAUT,TSST,MSST,MEZZ,MEZJ,
     &              GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,IKMCOR,
     &              IZERO,ITXRAY,NCSTMAX,NPOLMAX,LMAXYLM)
C
C=================================================== VBPES - ARPES - BIS
      ELSE IF ( TASK(1:5).EQ.'VBPES' .OR. TASK(1:5).EQ.'ARPES' .OR. 
     &          TASK(1:3).EQ.'BIS' ) THEN
C
         CALL VBPES(TASK,TSST,MSST,NE,MEZZ,MEZJ,IPRINT,PHASK,TAUQ,TAUT,
     &              CALCINT,GETIRRSOL,NPOLMAX)
C
C================================================================ WFPLOT
      ELSE IF ( TASK(1:7).EQ.'WFPLOT' ) THEN
C
         CALL WFPLOT(TSST,MSST,IPRINT,SSST,MEZZ,MEZJ,NCSTMAX)
C
C=============================================================== MECHECK
      ELSE IF ( TASK(1:7).EQ.'MECHECK' ) THEN
C
         CALL MECHECK(TSST,MSST,IPRINT,SSST,MEZZ,MEZJ,NCSTMAX)
C
C=============================================================== PLOTPRS
      ELSE IF ( TASK(1:7).EQ.'PLOTPRS' .AND. .NOT.POSANITAU ) THEN
C
         IF ( PLOTPRS(0) ) CALL FPPLOTPRS
C
C============================================================== PLOT2DPR
      ELSE IF ( TASK(1:8).EQ.'PLOT2DPR' .AND. .NOT.POSANITAU ) THEN
C
         IF ( PLOTPRS(0) ) CALL FPPLOTPRS
C
         CALL FPPLOT2DPR
C
C================================================================ MEPLOT
      ELSE IF ( TASK(1:6).EQ.'MEPLOT' ) THEN
C
         CALL ME_PLOT(NCSTMAX)
C
C==================================================================== T1
      ELSE IF ( TASK(1:2).EQ.'T1' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK does not YET work for FULL POTENTIAL mode')
C
         IF ( IREL.LE.1 ) THEN
            CALL NONRELT1(MEZJ,MEZZ,MSST,SSST,TAUT,ERYD,P,BCORS)
C
         ELSE
            CALL RELT1(TAUT,CALCINT,GETIRRSOL,IPRINT,TSST,MSST,MEZZ,
     &                 MEZJ)
         END IF
C
C================================================================= EKREL
      ELSE IF ( TASK(1:5).EQ.'EKREL' ) THEN
C
         CALL EKREL
C
C=================================================================== BSF
      ELSE IF ( TASK(1:3).EQ.'BSF' ) THEN
C
         CALL BLOCHSF
C
C============================================================= POSANITAU
      ELSE IF ( POSANITAU ) THEN
C
         IF ( .NOT.RHOPOS_ON_DISK .AND. MPI_ID.EQ.0 ) THEN
C
C------------------------------------------ perform corrections for E=EF
C                 needs from run for IE=NETAB:  TOTDOS, TAU- and WF-data
C                                               TOTNOS is recalculated
C
            TOTNOS = 0.0D0
            MUESPN = 0.0D0
            MUEORB = 0.0D0
C
            DO IT = ITBOT,ITTOP
               OBS_T(0,IDOS,IT) = 0.0D0
               OBS_T(0,ISMT,IT) = 0.0D0
               OBS_T(0,IOMT,IT) = 0.0D0
               DO IL = 1,NL
                  OBS_T(0,IDOS,IT) = OBS_T(0,IDOS,IT)
     &                               + OBS_LT(0,IDOS,IL,IT)
                  OBS_T(0,ISMT,IT) = OBS_T(0,ISMT,IT)
     &                               + OBS_LT(0,ISMT,IL,IT)
                  OBS_T(0,IOMT,IT) = OBS_T(0,IOMT,IT)
     &                               + OBS_LT(0,IOMT,IL,IT)
               END DO
               TOTNOS = TOTNOS + OBS_T(0,IDOS,IT)*CONC(IT)*NAT(IT)
               MUESPN = MUESPN + OBS_T(0,ISMT,IT)*CONC(IT)*NAT(IT)
               MUEORB = MUEORB + OBS_T(0,IOMT,IT)*CONC(IT)*NAT(IT)
            END DO
C
C ----------------------------------------------------------------------
C
            CALL CHRDNS(.TRUE.,CTOTDOS,SHFTEF,TOTDOSX,TOTDOS,TOTNOS,
     &                  ERYD,C0,EFERMILD,IEPANEL,IEPATH,IECURR,
     &                  OBS_POS_LT,OBS_POS_T,BCOR,BCORS,MEZZ,MEZJ,TSST,
     &                  MSST,TAUT,MSSQ,TAUQ,RHOCHRX,RHOSPNX,RHOORBX,
     &                  RHO2NSX,CMNTTX)
C
C
            CALL POSANIRHOIO('WRITE',RHOPOS_ON_DISK)
C
         END IF
C
         PLOTPRS(0:3) = .TRUE.
C
         CALL FPPLOTPRS
C
         IF ( FULLPOT ) THEN
C
            PLOT2DPR(1) = .FALSE.
            PLOT2DPR(2) = .TRUE.
C
            CALL FPPLOT2DPR
C
         END IF
C
C
         IF ( .NOT.GETTAUIJ ) CALL POSANILAMBDA(RHOCHR_EL,RHO2NS_EL)
C
C============================================================ POSANIPREP
      ELSE IF ( POSANIPREP ) THEN
C
         CALL POSANIPREPARE(ERYD,P,IWRREGWF,IWRIRRWF,IPRINT,MSSQ,MSST,
     &                      CALCINT,GETIRRSOL,TAUQ,TAUT,MEZZ,MEZJ)
C
C====================================================== COMPTON - POSANI
      ELSE IF ( TASK(1:7).EQ.'COMPTON' .OR. TASK(1:6).EQ.'POSANI' ) THEN
C
         CALL COMPTON(TASK,MSSQ,MSST,CALCINT,GETIRRSOL,TSST,TAUQ,MEZZ,
     &                MEZJ)
C
C================================================================ TORQUE
      ELSE IF ( TASK(1:6).EQ.'TORQUE' ) THEN
C
         CALL TORQUE
C
C================================================================== JALF
C
      ELSE IF ( CALC_JALF ) THEN
C
         CALL CURDNS(.TRUE.,RENORMALIZE,SHFTEF,SCLNOS,JDNST,C0,TAUT)
C
C================================================================== JORB
C
      ELSE IF ( CALC_JORB ) THEN
C
         CALL ORBCURDNS(.TRUE.,RENORMALIZE,SHFTEF,SCLNOS,C0,MSST,TAUT)
C
C================================================================== JORB
C
      ELSE IF ( TASK(1:6).EQ.'STONER' ) THEN
C
         ALLOCATE (AXCN(NRMAX,2,NLMAX,NTMAX),BXCNMM(NRMAX,NTMAX))
         ALLOCATE (BXCNMN(NRMAX,NTMAX),BXCNNM(NRMAX,NTMAX))
         ALLOCATE (BXCNNN(NRMAX,NTMAX),NEFTL(NTMAX,NLMAX))
         ALLOCATE (EMM(NRMAX,NTMAX),ENM(NRMAX,NTMAX),ENN(NRMAX,NTMAX))
         ALLOCATE (GAMMAM(NRMAX,NTMAX),GAMMAN(NRMAX,NTMAX))
         ALLOCATE (ORBSQINT(0:NLMAX,2))
C
         DO IT = 1,NT
            DO IL = 1,NL
               NEFTL(IT,IL) = DIMAG(DOBS_LTX(0,IDOS,IL,IT))
            END DO
         END DO
C
         CALL CHIGAMMA(IFILCBWF,GAMMAM,GAMMAN,BXCNMM,BXCNNM,BXCNMN,
     &                 BXCNNN,AXCN,NEFTL,ORBSQINT,EMM,ENM,ENN)
C
C================================================================== TEST
      ELSE IF ( ITEST.EQ.7 ) THEN
         CALL CHECK_GFUN
C         CALL CHECK_GFUN(MSSQ,MSST,TSST,SSST,TAUQ,MEZZ,MEZJ)
C
      END IF
C
      IF ( MPI ) THEN
         CALL DRV_MPI_BARRIER
         CALL MPI_FINALIZE(IERR)
      END IF
C
      IF ( CALCGFUNMAT ) DEALLOCATE (GFMAT)
      DEALLOCATE (HINT)
      DEALLOCATE (CPACHTAB,CTOTDOS,IECPAFAIL,EREF)
      DEALLOCATE (BCOR,BCORS)
      DEALLOCATE (FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR)
      DEALLOCATE (NKPCOR)
C
      RETURN
C
C=======================================================================
 200  CONTINUE
      WRITE (6,*) '      end of file 9 reached for IECURR= ',IECURR
      CALL STOP_MESSAGE(ROUTINE,'check TAU-file')
C=======================================================================
99001 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),/)
99002 FORMAT (/,79('*'),/,10X,' TAU(LAM,LAM'')',/,79('*'),/,5X,' NT =',
     &        I3,' NQ =',I3,' FMT = 3')
99003 FORMAT (/,79('*'),/,10X,' TAU(LAM,LAM'') and M(LAM,LAM'')',/,
     &        79('*'),/,13X,' NQ =',I3,' FMT = 3')
99004 FORMAT (/,' CPU - time used in  <',A,'>: ',F10.3,' sec',/)
99005 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99006 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99007 FORMAT (/,1X,79('*'),/,10X,'correction to',
     &        '  E_Fermi > 1mRy   NEWSTART !',/,1X,79('*'),/)
99008 FORMAT (/,1X,79('*'),/,10X,
     &        'E_Fermi not found within 4 iterations !',/,1X,79('*'),/)
99009 FORMAT (/,1X,79('*'),/,10X,'THIS TASK DOES NOT RUN IN MPI MODE !',
     &        /,10X,'CHANGING TO SEQUENTIAL MODE         ',/,1X,79('*'),
     &        /)
      END
C*==mpi_merge_files.f    processed by SPAG 6.70Rc at 17:10 on 19 Apr 2017
      SUBROUTINE MPI_MERGE_FILES(IFIL,FILNAM0,MPI_OUTPUT_FLAG,KEY)
C   ********************************************************************
C   *                                                                  *
C   *  deal with DOS and TAU files in case of MPI run of GEN           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI,MPI_ID,NPROCS
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*80 FILNAM0
      INTEGER IFIL,KEY
      LOGICAL MPI_OUTPUT_FLAG(0:(NPROCS-1))
C
C Local variables
C
      CHARACTER*80 FILNAM_ID
      INTEGER IPROC
      CHARACTER*120 LINE
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.MPI ) RETURN
C
      FILNAM_ID = FILNAM0(1:LEN_TRIM(FILNAM0))//'_'
C
      CALL STRING_ADD_N(FILNAM_ID,MPI_ID)
C
C=======================================================================
C                            initialize
C=======================================================================
      IF ( KEY.EQ.0 ) THEN
C
         WRITE (6,*) 'PROC ',MPI_ID,' opening ',IFIL,FILNAM_ID
C
         IF ( MPI_OUTPUT_FLAG(MPI_ID) ) THEN
C
            CLOSE (IFIL)
C
            OPEN (IFIL,FILE=FILNAM_ID)
C
         END IF
C
C=======================================================================
C               close file and collect in master file
C         go through list of writing processes in reverse order
C=======================================================================
      ELSE
C
         CLOSE (IFIL)
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
         IF ( MPI_ID.EQ.0 ) THEN
C
            OPEN (IFIL,FILE=FILNAM0,POSITION='APPEND')
C
            DO IPROC = NPROCS - 1,0, - 1
C
               IF ( MPI_OUTPUT_FLAG(IPROC) ) THEN
C
                  FILNAM_ID = FILNAM0(1:LEN_TRIM(FILNAM0))//'_'
C
                  CALL STRING_ADD_N(FILNAM_ID,IPROC)
C
                  OPEN (IOTMP,FILE=FILNAM_ID,STATUS='OLD')
 5                CONTINUE
                  READ (IOTMP,99001,END=10) LINE
                  WRITE (IFIL,99001) LINE
                  GOTO 5
C
 10               CONTINUE
                  CLOSE (IOTMP,STATUS='DELETE')
C
               END IF
C
            END DO
C
         END IF
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
         CLOSE (IFIL)
         OPEN (IFIL,FILE=FILNAM0,STATUS='OLD')
C
      END IF
99001 FORMAT (A120)
      END
