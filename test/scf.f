C*==scf.f    processed by SPAG 6.70Rc at 08:44 on  8 Mar 2017
      SUBROUTINE SCF(KBZI,KMOL,FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine program for       KKRSCF                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_DIMENSION,SYSTEM_TYPE,ALAT
      USE MOD_CPA,ONLY:CPATOL,NCPA,ITCPAMAX
      USE MOD_RMESH,ONLY:NMMAX,NRMAX,FULLPOT,JRNSMIN,JRNS1,JRCRI,JRWS,
     &    R2DRDI,R
      USE MOD_SYMMETRY,ONLY:NELMTMAX
      USE MOD_KSPACE,ONLY:NELMT,IBZINT
      USE MOD_SITES,ONLY:NQMAX,DROTQ,ITOQ,NOQ,IQAT,QMPHI,QMTET,QMVEC,
     &    NQ_L,NQ_I,IQBOT,IQTOP,NQHOST,NQ,CMNTQ,VLMMAD_HOST,NLMQMAD
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NEPATH,SPLITSS,IGRID,NETAB,EILOW,
     &    EIMAG,EMIN,EFERMI,SEARCHEF,WETAB,NEPANEL,PHASA,PHASK,PHAST
      USE MOD_FILES,ONLY:LSYSTEM,SYSTEM,STRINP,LDATSET,DATSET,IFILCBWF,
     &    WRLOG,IFILLOG,IPRINT,WRMAT,WRTAU,WRTAUMQ,PLOTPRS,POTFIL,
     &    LPOTFIL,POTFIL_CLU,LPOTFIL_CLU,IFILSCLALAT,DEBUG,IFILBREAK,
     &    FOUND_SECTION,FOUND_STRING
      USE MOD_TYPES,ONLY:NLMFPMAX,BEXT,NTMAX,NLMFPT,BNST,VNST,NAT,Z,IMT,
     &    BT,VT,LTXT_T,TXT_T,SEMCORSHLT,NSEMCORSHLT,ECORTAB,NCORT,
     &    RHOSPN,RHOCHR,CONC,NT,RHO2NS,RHOORB,VMTZ,QEL,ABIT,AOPT,ITBOT,
     &    ITTOP,CMNTT,ETOT,MUEORB,MUESPN,JDNST,OBS_T,OBS_LT,OBS_TX,
     &    OBS_LTX,OBS_X,W_WEISS_T
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,NMVECMAX,NKMMAX,NLMAX,NKM,NK,NL,
     &    AMEBI1,AMEBI2,GBIG,GBIL,NLABIMAX,MEZJ,MEZZ,TSSQ,MSSQ,MSST,
     &    SSST,TAUQ,TAUT,TSST
c modified by XJQ: parallel on energies
      USE MOD_MPI,ONLY:NPROCS,MPI_ID,MPI_ELOOP
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq
      USE MOD_CALCMODE,ONLY:PROGNAME,DMFT,NONMAG,KMROT,BLCOUPL,ORBPOL,
     &    BREITINT,IREL,ITEST,LLOYD,SCL_ALAT,FULLCHARGE,LDAU,BREAKPOINT,
     &    SOLVER_FP,BORN_USE_SPH,THERMAL_VIBRA_FLUCT,CALC_EFG,
     &    WAVE_FUNCTIONS_AVAILABLE
      USE MOD_THERMAL,ONLY:I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT,
     &    FLUCT_DIR_SETTING
      USE MOD_TB,ONLY:L_R_RELATION,IQ_L_IQ_R
      USE MOD_SCF,ONLY:SCFSTATUS,SCFSTATUS_CLU,SCFALG,SCFVXC,RMSAVA,
     &    RMSAVB,RMSAVV,SCFMIX,SCFMIXOP,SCFSIM,SCFTOL,ITRSCF,NSCFITER,
     &    SCF_CHECK_SPLITSS,NFCOR,W_WEISS_ERROR,SCF_THETA_DEPENDENT_POT
      USE MOD_DMFT_LDAU,ONLY:SEVT,SEVNST,SEBT,DMFTSIG,DMFTSIGMA,
     &    EREFLDAU,EREFDMFT,SEBNST,KSELF,ELDAU,DMFT_FIX_DYN_SE,DMFTSCF,
     &    DMFTDBLC
      USE MOD_CONSTANTS,ONLY:C0
c modified by XJQ: scf of vibrations
      use mod_scfvb_cpa_sigma,only:lscfvb
c end-mod-xjq
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF')
      LOGICAL COPY_POT_LB_TO_IZ
      PARAMETER (COPY_POT_LB_TO_IZ=.TRUE.)
      REAL*8 SHFTEFLIM
      PARAMETER (SHFTEFLIM=1.0D0)  !REC changed 0.1 to 1.0
C
C Dummy arguments
C
      LOGICAL FULLPOT5
      INTEGER KBZI,KMOL
C
C Local variables
C
      REAL*8 ABITNEW(:,:,:,:),AOPT0(:,:,:),BCOR(:),BCORS(:),BINT1,BINTT,
     &       BSCAL,BSUM_IZ,BSUM_LB,CMNTTX(:,:),CPACHNG,CPACHTAB(:),
     &       CPU_TIME_SCF_ITER_END,CPU_TIME_SCF_ITER_START,DEL,DMFTMIX,
     &       DQMVEC,EBANDLD,ECOR(:),ECTOP,EFERMI0,EFERMILD,EMAXSEMCOR,
     &       EMINSEMCOR,ERRAVANG,ERYDTOP,FCOR(:,:,:),GCOR(:,:,:),
     &       LOG10RMSB,LOG10RMSV,MVGAM(:,:),MVPHI(:,:),MVTET(:,:),
     &       NSEMCOR,QLZ(:,:),QMDIR(3),QMGAM(:),QMVEC0,RAUX,RDUMARG,
     &       RHO2NSX(:,:,:,:),RHOCHRX(:,:),RHOOPC(:,:,:),RHOOPO(:,:,:),
     &       RHOORBX(:,:),RHOSPNX(:,:),RINT(:),RMSERA,RMSERAS,RWK(:),
     &       RWKLMT(:,:),RWORK(:,:),SCLNOS,SHFTEF,SZCOR(:),TOTDOS,
     &       TOTDOSEF,TOTDOSEFX,TOTDOSX,TOTNOS,TOTNOSX,W1MIX(:,:),
     &       W1NSMIX(:,:,:),W2MIX(:,:),W2NSMIX(:,:,:),WW1,WW2
      LOGICAL CALCINT,COREHOLE,FOUND,GETIRRSOL,ITERMVDIR,MIXRHO,
     &        NOCHKVAL,OPDONE,RENORMALIZE,SCFDONE,SCFRESTART,SEMICORE,
     &        WRSSDOS
      COMPLEX*16 CTOTDOS(:),CWKE(:),CWKKMTE(:,:,:),CWKT(:),EGFNEW(:),
     &           EGFOLD(:),ERYD,ETABSS_BI(:),GFMAT(:,:,:,:),P,WEGFOLD(:)
     &           ,WETABSS_BI(:)
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IDIMS,IE,IECPAFAIL(:),IECURR,
     &        IEPANEL,IEPATH,IERR,IKMCOR(:,:),IL,ILA,IM,IMIX,IMV,INC,IO,
     &        IO_IZ,IO_LB,IPRINTBAND,IPRINTCORE,IPROC,IPROCE(:),IQ,
     &        IQBOTSAV,IQMVEC,IQMVEC1,IQMVEC2,IQTOPSAV,IQ_IZ,IQ_LB,
     &        IQ_RBULK,IR,IS,ISTBRY,IT,ITBOTSAV,ITCPA,ITDEPT,ITRSCF0,
     &        ITTOPSAV,ITXRAYDUM,IT_IZ,IT_LB,IT_RBULK,IVEC,IWRIRRWF,
     &        IWRREGWF,IZERO(:),JRNS1TMP(:),JRNSMINTMP,JTOP,KAPCOR(:),
     &        LL,LM,LWK,LWKE,LWKKMTE,LWKLMT,M,MM05COR(:),N,NCORT0(:),
     &        NCPAFAIL,NCSTMAX,NESEMCOR,NESS,NETAB1VB,NKPCOR(:),NVEC
      CHARACTER*2 MAGCOUPL
c modified by XJQ: parallel on energies
      logical lparalleled, lopen
      integer ie0, ie1, ne_mpi, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      integer, dimension(:), allocatable :: disp_proc, ne_proc
      real*8, allocatable, dimension(:,:,:) :: array3d_real
      real*8, allocatable, dimension(:,:,:,:) :: array4d_real
      complex*16, allocatable, dimension(:,:,:) :: gfmat_tmp
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
      DATA SHFTEF/0.0D0/,RENORMALIZE/.FALSE./,ITXRAYDUM/0/
      DATA RDUMARG/999999D0/,EBANDLD/0D0/,WRSSDOS/.FALSE./
      DATA SEMICORE/.FALSE./,COREHOLE/.FALSE./,SCLNOS/0D0/
C
      ALLOCATABLE FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR
      ALLOCATABLE NKPCOR
      ALLOCATABLE GFMAT,EGFOLD,EGFNEW
      ALLOCATABLE AOPT0,ABITNEW,RHO2NSX,RHOCHRX,RHOOPC,RHOOPO
      ALLOCATABLE ETABSS_BI,WETABSS_BI
      ALLOCATABLE RHOORBX,RHOSPNX,RWORK
      ALLOCATABLE CPACHTAB,CTOTDOS,CWKE,CWKKMTE,IECPAFAIL,IPROCE
      ALLOCATABLE MVTET,MVGAM,MVPHI
      ALLOCATABLE RWK,RWKLMT
      ALLOCATABLE CMNTTX
C
C------------------------------------------ variables depending on NQMAX
C
      ALLOCATABLE QMGAM
C
C------------------------------------------ variables depending on NTMAX
C
      ALLOCATABLE NCORT0
      ALLOCATABLE QLZ,BCOR,BCORS
      ALLOCATABLE CWKT
      ALLOCATABLE WEGFOLD
      ALLOCATABLE W1MIX,W2MIX,W1NSMIX,W2NSMIX,JRNS1TMP
C
C------------------------------------------ variables depending on NRMAX
C
      ALLOCATABLE RINT
C
      CALL TRACK_INFO(ROUTINE)
C
C------------------------------------------- variables depending on NCST
C
C  these arrays are needed only if information on a specific core shell
C   has to be obtained from <CORE> e.g. for core level spectroscopies
C-----------------------------------------------------------------------
C
      NCSTMAX = 14
C
      ALLOCATE (IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX))
      ALLOCATE (KAPCOR(NCSTMAX),MM05COR(NCSTMAX))
      ALLOCATE (IZERO(NCSTMAX),SZCOR(NCSTMAX),ECOR(NCSTMAX))
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SZCOR')
C
C----------------------------------------- variables depending on NMEMAX
C
      ALLOCATE (MVTET(NTMAX,NMVECMAX))
      ALLOCATE (MVGAM(NTMAX,NMVECMAX))
      ALLOCATE (MVPHI(NTMAX,NMVECMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MVPHI')
C
C------------------------------------------ variables depending on NEMAX
C
      ALLOCATE (IECPAFAIL(NEMAX))
      ALLOCATE (CPACHTAB(NEMAX),CTOTDOS(NEMAX),CWKE(NEMAX))
      ALLOCATE (CWKKMTE(NKMMAX,NTMAX,NEMAX),IPROCE(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CWKKMTE')
C
C-----------------------------------------------------------------------
C
      IF ( FULLPOT .OR. FULLPOT5 ) THEN
         ALLOCATE (RHOORBX(1,1))
         ALLOCATE (RHOCHRX(1,1),RHOSPNX(1,1))
         ALLOCATE (RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3))
      ELSE
         ALLOCATE (RHOORBX(NRMAX,NTMAX))
         ALLOCATE (RHOCHRX(NRMAX,NTMAX),RHOSPNX(NRMAX,NTMAX))
         ALLOCATE (RHO2NSX(1,1,1,1))
      END IF
      ALLOCATE (RWORK(NRMAX,NTMAX))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RWORK')
C
C-------------------------------- variables depending on NLMAX and NTMAX
C
      ALLOCATE (CMNTTX(NLMFPMAX,NTMAX))
      ALLOCATE (W1MIX(NRMAX,NTMAX),W2MIX(NRMAX,NTMAX))
C
C------------------------------------------ variables depending on NQMAX
C
      ALLOCATE (QMGAM(NQMAX))
C
C------------------------------------------ variables depending on NTMAX
C
      ALLOCATE (NCORT0(NTMAX))
      ALLOCATE (QLZ(NTMAX,2))
      ALLOCATE (BCOR(NTMAX),BCORS(NTMAX))
      ALLOCATE (CWKT(NTMAX))
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NCORT0')
C
C------------------------------------------ variables depending on NRMAX
C
      ALLOCATE (RINT(NRMAX))
C
C-----------------------------------------------------------------------
C
      IF ( IBZINT.EQ.0 .AND. KBZI.EQ.1 ) CALL STOP_MESSAGE(ROUTINE,
     &     'TAU  via  cluster  calculation BUT KBZI=1')
      IF ( IBZINT.NE.0 .AND. KBZI.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'TAU  via  BZ-integration BUT KBZI=0')
      IF ( ITEST.EQ.5 .AND. NPROCS>1 )
     &      CALL STOP_MESSAGE(ROUTINE,'NO MPI when checking Lloyd')
c modified by XJQ: parallel on energies using my own subroutines
c      IF ( NPROCS>1 .AND. NETAB(1).LT.NPROCS ) THEN
c         WRITE (6,99024) NETAB(1),NPROCS
c         CALL STOP_MESSAGE(ROUTINE,
c     &        'parallel E-loop requires: num. of procs. NPROCS <= NETAB'
c     &        )
c      END IF
c end-mod-xjq
C
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'parallel E-loop requires: num. of procs. NPROCS <= NETAB'
     &        )
         IF ( SPLITSS ) CALL STOP_MESSAGE(ROUTINE,
     &        'SPLISS not allowed for calculation mode')
      ELSE IF ( BREITINT ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'NO full potential mode allowed for calculation mode')
         IF ( SPLITSS ) CALL STOP_MESSAGE(ROUTINE,
     &        'SPLISS not allowed for calculation mode')
      END IF
C
      P = C0
C
      IMV = 1
      DO IQ = IQBOT,IQTOP
         QMGAM(IQ) = 0D0
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            MVPHI(IT,IMV) = QMPHI(IQ)
            MVTET(IT,IMV) = QMTET(IQ)
            MVGAM(IT,IMV) = QMGAM(IQ)
         END DO
      END DO
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
      IF(NPROCS>1) THEN
        MPI_ELOOP=.true.
      else
        mpi_eloop=.false.
      endif
C
      IF ( MPI_ID.NE.0 ) THEN
C
         WRLOG = .FALSE.
         WRTAU = .FALSE.
         WRTAUMQ = .FALSE.
C
      END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C=======================================================================
C                   open and position log-file
C=======================================================================
C
      IF ( WRLOG ) THEN
         OPEN (UNIT=IFILLOG,FILE=DATSET(1:LDATSET)//'.log')
C
         N = 0
         ITRSCF0 = 0
         DO I = 1,10000
            READ (IFILLOG,'(I5)',END=50,ERR=20) ITRSCF0
 20         CONTINUE
            N = N + 1
         END DO
C
 50      CONTINUE
         REWIND (IFILLOG)
         DO I = 1,N
            READ (IFILLOG,*)
         END DO
C
         WRITE (IFILLOG,99001)
      END IF
C
C=======================================================================
C
C      ITERMVDIR = .TRUE.
      ITERMVDIR = .FALSE.
C
      DO IT = ITBOT,ITTOP
         NCORT0(IT) = NCORT(IT)
      END DO
      NETAB1VB = NETAB(1)
C
C---------------------------------------------------- k-integration mesh
C
      IF ( IBZINT.EQ.3 ) THEN
         IF ( NELMT.GT.NELMTMAX ) THEN
            WRITE (6,99004) 'number of elements =',NELMT,
     &                      ' > array size=',NELMTMAX
            CALL STOP_MESSAGE(ROUTINE,'NELMT > NELMTMAX')
         END IF
C
      END IF
C
C ----------------------------------------------------------------------
C                       test Lloyd formula
C ----------------------------------------------------------------------
C
      IF ( ITEST.EQ.5 ) CALL LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,MSST,SSST,
     &                                C0,0)
C
      ERYDTOP = 1.4D0*EFERMI
C
      IF ( LLOYD .AND. PROGNAME.EQ.'KKRSCF    ' )
     &     CALL FRELPOLE(ERYDTOP,IPRINT)
C
C ----------------------------------------------------------------------
C
      IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
C
         POTFIL_CLU = POTFIL_CLU(1:LPOTFIL_CLU)//'_new'
         LPOTFIL_CLU = LPOTFIL_CLU + 4
C
      END IF
C
C ======================================================================
C                                            prepare SCF - cycle  STEP 1
C
      IF ( SPLITSS ) THEN
         GETIRRSOL = .FALSE.
      ELSE
         GETIRRSOL = .TRUE.
      END IF
      CALCINT = .TRUE.
      NSCFITER = 200
      SCFMIX = 0.20D0
      SCFMIXOP = SCFMIX
      SCFTOL = 1D-5
      ISTBRY = 1
      ITDEPT = 40
      SCFSIM = 0.0D0
      MAGCOUPL = 'NO'
      SCFRESTART = .FALSE.
      MIXRHO = .FALSE.
      NOCHKVAL = .FALSE.
      SCFDONE = .FALSE.
C
      CALL INPUT_FIND_SECTION('SCF',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_STRING('VXC',SCFVXC,'9999',0)
         CALL SECTION_SET_STRING('ALG',SCFALG,'9999',0)
         IF ( .NOT.FOUND_STRING )
     &        CALL SECTION_SET_STRING('SCFALG',SCFALG,'9999',0)
         CALL SECTION_SET_INTEGER('NITER',NSCFITER,9999,0)
         CALL SECTION_SET_REAL('MIX',SCFMIX,9999D0,0)
         CALL SECTION_SET_REAL('MIXOP',SCFMIXOP,9999D0,0)
         CALL SECTION_SET_REAL('TOL',SCFTOL,9999D0,0)
         CALL SECTION_SET_REAL('SIM',SCFSIM,9999D0,0)
         CALL SECTION_SET_INTEGER('ISTBRY',ISTBRY,9999,0)
         CALL SECTION_SET_INTEGER('ITDEPT',ITDEPT,9999,0)
         CALL SECTION_FIND_KEYWORD('SEMICORE',SEMICORE)
         CALL SECTION_FIND_KEYWORD('SCFRESTART',SCFRESTART)
         CALL SECTION_SET_STRING('MAGCOUPL',MAGCOUPL,'9999',0)
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' )
     &         CALL SECTION_FIND_KEYWORD('NFCOR',NFCOR)
         CALL SECTION_FIND_KEYWORD('MIXRHO',MIXRHO)
         CALL SECTION_FIND_KEYWORD('NOCHKVAL',NOCHKVAL)
         IF ( MAGCOUPL.NE.'FM' .AND. MAGCOUPL.NE.'AF' .AND. 
     &        MAGCOUPL.NE.'NO' ) THEN
            WRITE (6,*) 'MAGCOUPL = ',MAGCOUPL
            CALL STOP_MESSAGE(ROUTINE,'should be  FM  or  AF')
         END IF
         CALL SECTION_FIND_KEYWORD('TETDEPPOT',SCF_THETA_DEPENDENT_POT)
      END IF
      IF ( SCFSIM.LT.0D0 .OR. SCFSIM.GT.1D0 ) SCFSIM = 0.5D0
C
      IF ( SCFALG(1:8).EQ.'BROYDEN2' ) THEN
         IMIX = 4
      ELSE IF ( SCFALG(1:8).EQ.'ANDERSON' ) THEN
         IMIX = 7
      ELSE IF ( SCFALG(1:6).EQ.'TCHEBY' ) THEN
         IMIX = 4
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'TCHEBY algorithm not available for FULLP')
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'SCFALG not found')
      END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
      IF ( MPI_ID.EQ.0 ) THEN
         WRITE (6,99022)
         WRITE (6,99012) 'pot-file  :  ',POTFIL(1:LPOTFIL)
         WRITE (6,99012) 'XC-pot    :  ',SCFVXC
         WRITE (6,99012) 'SCF-alg   :  ',SCFALG
         WRITE (6,99013) 'SCF-niter :  ',NSCFITER
         WRITE (6,99014) 'SCF-mix   :  ',SCFMIX
         IF ( NCPA.NE.0 .AND. ABS(SCFSIM).GT.1D-6 ) WRITE (6,99014)
     &         'SCF-SIM   :  ',SCFSIM
         IF ( ORBPOL(1:4).NE.'NONE' ) WRITE (6,99014) 'SCF-mix-OP:  ',
     &        SCFMIXOP
         WRITE (6,99014) 'SCF-tol   :  ',SCFTOL
         WRITE (6,99015) 'BREIT-int :  ',BREITINT
         WRITE (6,99012) 'OP-scheme :  ',ORBPOL
         WRITE (6,99015) 'FULLPOT   :  ',(FULLPOT .OR. FULLPOT5)
         WRITE (6,99015) 'LLOYD     :  ',LLOYD
         WRITE (6,99015) 'NFCOR     :  ',NFCOR
         WRITE (6,99015) 'TETDEPPOT :  ',SCF_THETA_DEPENDENT_POT
         IF ( (SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV')
     &        ) WRITE (6,99015) 'POT LB->IZ:  ',COPY_POT_LB_TO_IZ
C
      END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C=======================================================================
C                  set up the Madelung matrices
C=======================================================================
C
      CALL SCFMAD_SETUP(FULLPOT5)
C
      ALLOCATE (RWKLMT(NLMFPMAX,NTMAX))
C
      IF ( NCPA.NE.0 .AND. ABS(SCFSIM).GT.1D-6 ) CALL GETRNN
C
C ======================================================================
C                      write start - potential for cluster
C ======================================================================
C
      IF ( KMOL.EQ.1 ) THEN
         CALL INPUT_FIND_SECTION('SCF',0)
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('POTUPD',FOUND)
            IF ( FOUND ) THEN
C
               CALL POTWRSPR
C
               STOP
C
            END IF
         END IF
      END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 1,5,6
      IF ( BREAKPOINT.EQ.1 ) CALL SFNDUMP
C
C      VNST=VNST
C      BNST=0d0
C      BT=0d0
C      igrid=0
C      IBZINT=1
C      NPTMAX=0
C
      IF ( BREAKPOINT.EQ.5 ) THEN
         IF ( IREL.NE.3 )
     &         CALL STOP_MESSAGE(ROUTINE,'IREL <> 3 for BREAK = 5')
         NOCHKVAL = .TRUE.
         OPEN (UNIT=IFILBREAK,STATUS='SCRATCH',FORM='UNFORMATTED')
         WRITE (IFILBREAK) EMIN,EFERMI,VT,BT,VNST,BNST
      END IF
C
      IF ( BREAKPOINT.EQ.6 ) THEN
         IF ( IREL.NE.3 )
     &         CALL STOP_MESSAGE(ROUTINE,'IREL <> 3 for BREAK = 6')
         NOCHKVAL = .TRUE.
         OPEN (UNIT=IFILBREAK,FILE='BREAK_V_rotated.dat',
     &         FORM='UNFORMATTED')
         READ (IFILBREAK) EMIN,EFERMI,VT,BT,VNST,BNST
      END IF
      IF ( BREAKPOINT.EQ.5 .OR. BREAKPOINT.EQ.6 ) THEN
         WRITE (6,99029) NINT(QMTET(1)),NINT(QMPHI(1)),SOLVER_FP,
     &                   BREAKPOINT,EMIN,EFERMI
         IF ( SOLVER_FP(1:4).EQ.'BORN' ) WRITE (6,99030) NINT(QMTET(1)),
     &        NINT(QMPHI(1)),SOLVER_FP,BORN_USE_SPH
         WRITE (6,*) 'BREAK    BORN_USE_SPH',BORN_USE_SPH
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 1,5,6
C
C ======================================================================
C                                            prepare SCF - cycle  STEP 2
C
C=======================================================================
C                            SCF - START
C            initialize potential and find Fermi energy
C=======================================================================
C
      IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
      CALL FLUSH(6)
C
      IF ( SCFSTATUS(1:5).EQ.'START' .OR. SCFSTATUS_CLU(1:5)
     &     .EQ.'START' .OR. SCFRESTART ) THEN
C
         CALL INPUT_FIND_SECTION('CONTROL',0)
         IF ( FOUND_SECTION ) CALL SECTION_FIND_KEYWORD('WRSSDOS',
     &        WRSSDOS)
C
CccccccccIF ( FULLPOT5 ) CALL SCFMAD_SETUP(FULLPOT5)
C
         IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
C
            CALL CLUSCFINITPOT
C
         ELSE IF ( SUB_SYSTEM(1:6).NE.'I-ZONE' .OR. SYSTEM_TYPE(1:3)
     &             .NE.'LIR' ) THEN
C
            CALL SCFINITPOT(WRSSDOS,QEL,TAUQ,TAUT,NOCHKVAL)
C
         END IF
C
      ELSE IF ( .NOT.SEMICORE ) THEN
C
         IF ( .NOT.NOCHKVAL ) CALL SCFCHKNVAL(EMIN,ECTOP)
C
         DO IEPATH = 1,NEPATH
            CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,NETAB(IEPATH),
     &                 ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                 NEMAX)
         END DO
C
C
      END IF
C
C
      IF ( KMOL.EQ.1 ) THEN
         CALL INPUT_FIND_SECTION('SCF',0)
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('POTGEN',FOUND)
            IF ( FOUND ) CALL CLUPOTGEN
         END IF
      END IF
C
C ======================================================================
C                  prepare SCF - cycle                            STEP 3
C                     only for master process
C ----------------------------------------------------------------------
C
      IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
C=======================================================================
C     switch from ASA to FULL POTENTIAL calculation if requested
C      - update the Fermi energy according to additional VMTZ shift
C      - write out the new potential for testing purposes
C      - set up the radial mesh for full potential mode
C      - create the energy mesh according to new Fermi energy
C-----------------------------------------------------------------------
C
      IF ( .NOT.FULLPOT .AND. FULLPOT5 ) THEN
C
         CALL FPASASTART
C
         IF ( PLOTPRS(0) ) CALL FPPLOTPRS
C
         IF ( MPI_ID.EQ.0 ) CALL POTWRSPR
C
         DO IEPATH = 1,NEPATH
            CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,NETAB(IEPATH),
     &                 ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                 NEMAX)
         END DO
C
      END IF
C
C=======================================================================
C
      IF ( .NOT.ALLOCATED(SEVT) )
     &     ALLOCATE (SEVT(NRMAX,NTMAX),SEBT(NRMAX,NTMAX))
      CALL CINIT(NRMAX*NTMAX,SEBT)
      CALL CINIT(NRMAX*NTMAX,SEVT)
C
C=======================================================================
C                         DMFT-calculational mode
C=======================================================================
C
      IF ( DMFT .OR. LDAU ) THEN
C
         IF ( IREL.LT.2 ) CALL STOP_MESSAGE(ROUTINE,
     &        'LDAU/DMFT not implemented in MODE SREL')
C
         CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),DMFTMIX,NEMAX)
C
         ALLOCATE (GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX))
C
         GFMAT(1:NKMMAX,1:NKMMAX,1:NEMAX,1:NTMAX) = C0
C
C
         ALLOCATE (EGFOLD(NEMAX),EGFNEW(NEMAX))
         ALLOCATE (WEGFOLD(NEMAX))
C
         CALL ZCOPY(NETAB(1),ETAB(1,1),1,EGFOLD,1)
         CALL ZCOPY(NETAB(1),WETAB(1,1),1,WEGFOLD,1)
C
         CALL CINIT(NRMAX*NTMAX,SEBT)
         CALL CINIT(NRMAX*NTMAX,SEVT)
         CALL CINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,SEBNST)
         CALL CINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,SEVNST)
C
         IF ( .NOT.FULLPOT ) THEN
            DO IM = 1,NMMAX
               JRNS1(IM) = JRNSMIN
               JRCRI(IM) = JRWS(IM)
            END DO
         END IF
C
      END IF
C=======================================================================
C
C-----------------------------------------------------------------------
C                    print present potential - if requested
C-----------------------------------------------------------------------
C
      IF ( IPRINT.GT.0 .AND. MPI_ID .EQ. 0) THEN
         DO IT = ITBOT,ITTOP
            STRINP = 'V0_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &               (1:LTXT_T(IT))
            LL = 3 + LDATSET + 1 + LTXT_T(IT)
            OPEN (80,FILE=STRINP(1:LL))
            IM = IMT(IT)
            JTOP = JRWS(IM)
            DO IR = 1,JTOP
               WRITE (80,'(2ES15.6)') R(IR,IM),VT(IR,IT)*R(IR,IM)
            END DO
            CLOSE (80)
         END DO
      END IF
C
C-----------------------------------------------------------------------
C                      initialize SCF - mixing procedures
C-----------------------------------------------------------------------
C
      ITRSCF = 0
C
      IF ( FULLPOT ) THEN
         ALLOCATE (JRNS1TMP(NMMAX))
         IF ( MIXRHO ) THEN
            JRNSMINTMP = 1
            JRNS1TMP(1:NMMAX) = 1
         ELSE
            JRNSMINTMP = JRNSMIN
            JRNS1TMP(1:NMMAX) = JRNS1(1:NMMAX)
         END IF
         ALLOCATE (W1NSMIX(JRNSMINTMP:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (W2NSMIX(JRNSMINTMP:NRMAX,NLMFPMAX,NTMAX))
      END IF
C
      IF ( .NOT.(MIXRHO) ) THEN
         W1MIX(1:NRMAX,1:NTMAX) = VT(1:NRMAX,1:NTMAX)
         W2MIX(1:NRMAX,1:NTMAX) = BT(1:NRMAX,1:NTMAX)
         IF ( FULLPOT ) THEN
            W1NSMIX(JRNSMINTMP:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = VNST(JRNSMINTMP:NRMAX,1:NLMFPMAX,1:NTMAX)
            W2NSMIX(JRNSMINTMP:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = BNST(JRNSMINTMP:NRMAX,1:NLMFPMAX,1:NTMAX)
C
         END IF
      ELSE IF ( FULLPOT ) THEN
C
         W1MIX(1:NRMAX,1:NTMAX) = RHO2NS(1:NRMAX,1,1:NTMAX,1)
         W2MIX(1:NRMAX,1:NTMAX) = RHO2NS(1:NRMAX,1,1:NTMAX,2)
C
         W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &      = RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1)
         W2NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &      = RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,2)
      ELSE
C
         W1MIX(1:NRMAX,1:NTMAX) = RHOCHR(1:NRMAX,1:NTMAX)
         W2MIX(1:NRMAX,1:NTMAX) = RHOSPN(1:NRMAX,1:NTMAX)
C
      END IF
C
      IF ( FULLPOT ) THEN
C JM: Only after first iteration rho2ns is consistent
C               ----- no initialisation here ---
C
         IF ( .NOT.MIXRHO ) CALL SCFBROYPT1(IPRINT,ITRSCF,IMIX,ISTBRY,
     &        ITDEPT,SCFMIX,W1MIX,W2MIX,W1NSMIX,W2NSMIX,RMSAVV,RMSAVB,
     &        JRNS1TMP,JRNSMINTMP)
C
C
      ELSE IF ( SCFALG(1:8).EQ.'BROYDEN2' .OR. SCFALG(1:8)
     &          .EQ.'ANDERSON' ) THEN
         CALL SCFBROYDEN(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,SCFMIX,R,
     &                   R2DRDI,JRWS,IMT,W1MIX,W2MIX,RMSAVV,RMSAVB,
     &                   NMMAX,NRMAX,NTMAX)
      ELSE IF ( SCFALG(1:6).EQ.'TCHEBY' ) THEN
         CALL SCFTCHEBY(ITRSCF,VT,BT,SCFMIX,RMSAVV,RMSAVB)
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'SCFALG not found')
C
      END IF
C
      IF ( .NOT.(MIXRHO) ) THEN
         VT(1:NRMAX,1:NTMAX) = W1MIX(1:NRMAX,1:NTMAX)
         BT(1:NRMAX,1:NTMAX) = W2MIX(1:NRMAX,1:NTMAX)
         IF ( FULLPOT ) THEN
            VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = W1NSMIX(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
            BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = W2NSMIX(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
         END IF
      ELSE IF ( FULLPOT ) THEN
         RHO2NS(1:NRMAX,1,1:NTMAX,1) = W1MIX(1:NRMAX,1:NTMAX)
C
         RHO2NS(1:NRMAX,1,1:NTMAX,2) = W2MIX(1:NRMAX,1:NTMAX)
C
         RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1)
     &      = W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
C
         RHO2NS(JRNSMIN:NRMAX,2:NLMFPMAX,1:NTMAX,2)
     &      = W2NSMIX(JRNSMIN:NRMAX,2:NLMFPMAX,1:NTMAX)
      ELSE
         RHOCHR(1:NRMAX,1:NTMAX) = W1MIX(1:NRMAX,1:NTMAX)
         RHOSPN(1:NRMAX,1:NTMAX) = W2MIX(1:NRMAX,1:NTMAX)
      END IF
C
C------------------------------------- update of magnetisation direction
C
      IF ( ITERMVDIR ) CALL SCFITERANG(ITRSCF,ITOQ,DROTQ,MVPHI,MVTET,
     &                                 MVGAM,QMPHI,QMTET,QMGAM,NQHOST,
     &                                 NK,ERRAVANG,NQMAX,NTMAX,NMVECMAX,
     &                                 NKMMAX)
C
C----------------------------------------------------------------- BREIT
      IF ( BREITINT ) THEN
C
         ALLOCATE (ABITNEW(NRMAX,NLABIMAX,-1:+1,NTMAX))
         ALLOCATE (JDNST(NRMAX,NLABIMAX,-1:+1,NTMAX))
         JDNST(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0D0
C
         CALL SCFBIMAD
C
      END IF
C
C-----------------------------------------------------------------------
C                 initialize treatment of core holes
C-----------------------------------------------------------------------
C
      CALL SCFRHOCORHOL(COREHOLE,ITRSCF,IKMCOR,NKPCOR,KAPCOR,MM05COR,
     &                  IZERO,ECOR,SZCOR,GCOR,FCOR,BCOR,BCORS,NCSTMAX)
C
C-----------------------------------------------------------------------
C                      ORBITAL POLARISATION
C-----------------------------------------------------------------------
      IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL ) THEN
C
         ALLOCATE (RHOOPC(NRMAX,NTMAX,2),RHOOPO(NRMAX,NTMAX,2))
         ALLOCATE (AOPT0(NRMAX,2,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: AOPT0')
C
         AOPT(1:NRMAX,1:2,1:NTMAX) = 0D0
         AOPT0(1:NRMAX,1:2,1:NTMAX) = 0D0
C
      END IF
C                                            prepare SCF - cycle  STEP 3
C ======================================================================
C
C
C ======================================================================
C                   outer loop for spin spirals
C ======================================================================
C
      IQMVEC1 = 0
      IQMVEC2 = 0
      QMVEC0 = 0.0D0
      DQMVEC = 0.0D0
C
      CALL INPUT_FIND_SECTION('MODE',0)
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_INTEGER('IQMVEC1',IQMVEC1,9999,0)
         CALL SECTION_SET_INTEGER('IQMVEC2',IQMVEC2,9999,0)
         CALL SECTION_SET_REAL('QMVEC0',QMVEC0,9999D0,0)
         CALL SECTION_SET_REAL('DQMVEC',DQMVEC,9999D0,0)
      END IF
      IF ( KMROT.GE.3 ) THEN
         WRITE (6,99005) IQMVEC1,IQMVEC2,QMVEC0,DQMVEC
         WRITE (6,99006) (QMTET(IQ),IQ=IQBOT,IQTOP)
         WRITE (6,'(/)')
      END IF
C
      QMDIR(1:3) = QMVEC(1:3)
C
CQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQS
C                           SPIN SPIRAL LOOP
CQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQS
C
      DO IQMVEC = IQMVEC1,IQMVEC2
C
         IF ( KMROT.GE.3 ) THEN
            DO IE = 1,3
               QMVEC(IE) = QMVEC0 + DQMVEC*IQMVEC*QMDIR(IE)
            END DO
            WRITE (6,99007) IQMVEC,QMVEC
         ELSE
            QMVEC(1:3) = 0.0D0
         END IF
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
            NEPATH = 1
         END IF
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
C
C=======================================================================
C                         SCF - cycle  START
C=======================================================================
         CALL FLUSH(6)
C
         IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
         ITRSCF = 0
 150     CONTINUE
         ITRSCF = ITRSCF + 1
C
         EFERMI0 = EFERMI
C
C-----------------------------------------------------------------------
C   suppress aspherical potential in case of FULLCHARGE calculations
C-----------------------------------------------------------------------
C
         IF ( FULLCHARGE ) THEN
            VNST(:,:,:) = 0D0
            BNST(:,:,:) = 0D0
            WRITE (6,99027)
         END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C                                                   broadcast potentials
         IF ( NPROCS>1 ) THEN
C
            IDIMS = NRMAX*NTMAX
            CALL DRV_MPI_BCAST_R(0,VT(1,1),IDIMS)
            CALL DRV_MPI_BCAST_R(0,BT(1,1),IDIMS)
C
            IF ( FULLPOT ) THEN
               IDIMS = (NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX
               CALL DRV_MPI_BCAST_R(0,VNST(1,1,1),IDIMS)
               CALL DRV_MPI_BCAST_R(0,BNST(1,1,1),IDIMS)
            END IF
C
            IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
               IDIMS = NRMAX*NLABIMAX*2*NTMAX
               CALL DRV_MPI_BCAST_R(0,AOPT(1,1,1),IDIMS)
            END IF
C
            CALL DRV_MPI_BCAST_R(0,EFERMI,1)
            CALL DRV_MPI_BCAST_C(0,ETAB(1,1),NEMAX*2)
            CALL DRV_MPI_BCAST_C(0,WETAB(1,1),NEMAX*2)
C
         END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C-----------------------------------------------------------------------
C    RDLM calculations - broadcast the new Weiss field in case of MPI
C                      - update the probabilities
C-----------------------------------------------------------------------
         IF ( FLUCT_DIR_SETTING.EQ.'RDLM      ' ) THEN
            IF ( NPROCS>1 ) CALL DRV_MPI_BCAST_R(0,W_WEISS_T(1),NTMAX)
            CALL THERMAL_CONT
         END IF
C-----------------------------------------------------------------------
C
         SCFRESTART = .FALSE.
C
         CALL CPU_TIME(CPU_TIME_SCF_ITER_START)
C
         IF ( ITRSCF.NE.1 ) THEN
            IF ( WRTAU ) THEN
               REWIND 9
               WRITE (9,99020) NT,NQHOST
               WRITE (9,'(5X,'' IT ='',I3,'' IQ ='',I3,'': '',A)')
     &                (IT,IQAT(1,IT),TXT_T(IT),IT=1,NT)
               DO IQ = IQBOT,IQTOP
                  WRITE (9,'(''site'',I3,'' of '',A)') IQ,
     &                   SYSTEM(1:LSYSTEM)
               END DO
            ELSE IF ( WRTAUMQ ) THEN
               REWIND 9
               WRITE (9,99021) NQHOST
               DO IQ = IQBOT,IQTOP
                  WRITE (9,'(''site'',I3,'' of '',A)') IQ,
     &                   SYSTEM(1:LSYSTEM)
               END DO
            END IF
         END IF
C
C--------------------------------------------------- free electron poles
C
         IF ( LLOYD .AND. (ITRSCF.LE.1 .OR. EFERMI.GT.ERYDTOP) ) THEN
C
            ERYDTOP = 1.4D0*EFERMI
C
            IF ( PROGNAME.EQ.'KKRSCF    ' )
     &           CALL FRELPOLE(ERYDTOP,IPRINT)
C
            SEARCHEF = .TRUE.
C
         END IF
C
C=======================================================================
C            initialize quantities determined by E-integration
C=======================================================================
C
         TOTNOS = 0D0
         TOTDOSEF = 0D0
C
C------------------------------------------------------- total densities
         IF ( FULLPOT .OR. FULLPOT5 ) THEN
            RHO2NS(:,:,ITBOT:ITTOP,:) = 0.0D0
         ELSE
            RHO2NS(:,:,:,:) = 0.0D0
         END IF
         RHOCHR(:,ITBOT:ITTOP) = 0.0D0
         RHOSPN(:,ITBOT:ITTOP) = 0.0D0
         RHOORB(:,ITBOT:ITTOP) = 0.0D0
C
C--------------------------------- moments, hyperfine field, band energy
         OBS_T(:,:,:) = 0D0
         OBS_LT(:,:,:,:) = 0D0
         CMNTT(:,:) = 0D0
C
C---------- spin projected densities connected with orbital polarisation
C
         IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL ) THEN
            RHOOPC(:,:,:) = 0D0
            RHOOPO(:,:,:) = 0D0
         END IF
C
C----------------------------------------------------------------- BREIT
         IF ( BREITINT ) ABITNEW(:,:,:,:) = 0.0D0
C
         IF ( IPRINT.LE.0 ) THEN
            IPRINTCORE = -1
            IPRINTBAND = -1
         ELSE
            IPRINTCORE = IPRINT
            IPRINTBAND = IPRINT
         END IF
C
C ======================================================================
C                             CORE STATES
C                      only for master process
C ======================================================================
C
C            DO IT = ITBOT,ITTOP
C               NCORT(IT) = NCORT0(IT)
C            END DO
C
C
         IF ( MPI_ID.EQ.0 ) CALL CORE(IPRINTCORE,GCOR,FCOR,ECOR,SZCOR,
     &                                KAPCOR,MM05COR,NKPCOR,IKMCOR,
     &                                IZERO,ITXRAYDUM,BCOR,BCORS,
     &                                NCSTMAX)
C
C-----------------------------------------------------------------------
C                                                  treat semicore states
C
C        IF ( SEMICORE ) THEN
         IF ( NT.LT.0 ) THEN
C
            CALL SCFSEMINIT(SEMICORE,NT,NL,NCORT0,NCORT,EMIN,EMINSEMCOR,
     &                      EMAXSEMCOR,NESEMCOR,ECORTAB,NSEMCORSHLT,
     &                      SEMCORSHLT,NTMAX)
C
            NSEMCOR = 0D0
            DO IT = ITBOT,ITTOP
               NSEMCOR = NSEMCOR + (NCORT0(IT)-NCORT(IT))*NAT(IT)
     &                   *CONC(IT)
            END DO
C
            CALL CORE(IPRINTCORE,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,
     &                NKPCOR,IKMCOR,IZERO,ITXRAYDUM,BCOR,BCORS,NCSTMAX)
C
C--------------------------------------------------- rewrite energy mesh
            NETAB(1) = NESEMCOR + NETAB1VB
C
            IEPATH = 1
C
            CALL EPATH(5,EMINSEMCOR,EMAXSEMCOR,EIMAG,NESEMCOR,
     &                 ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                 NEMAX)
C
            CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,NETAB1VB,
     &                 ETAB(NESEMCOR+1,IEPATH),WETAB(NESEMCOR+1,IEPATH),
     &                 EILOW,IPRINT,NEMAX)
C
            IF ( NEPATH.EQ.2 )
     &            CALL STOP_MESSAGE(ROUTINE,'NEPATH = 2 for SEMICORE')
C
         END IF
C                                                               SEMICORE
C-----------------------------------------------------------------------
C
C ======================================================================
C                       account for CORE HOLE
C ======================================================================
C
         IF ( COREHOLE ) CALL SCFRHOCORHOL(COREHOLE,ITRSCF,IKMCOR,
     &        NKPCOR,KAPCOR,MM05COR,IZERO,ECOR,SZCOR,GCOR,FCOR,BCOR,
     &        BCORS,NCSTMAX)
C
         NCPAFAIL = 0
C
         IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
         IF ( .NOT.(DMFTSCF) ) THEN
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP EPANEL
            DO IEPANEL = 1,NEPANEL
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP EPATH
               DO IEPATH = 1,NEPATH
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                  IPROCE(1:NEMAX) = 0
C
                  IF ( NPROCS.GT.1 ) THEN
                     IPROC = 0
                     INC = 1
                     DO IE = 1,NETAB(IEPATH) - 1
                        IPROC = IPROC + INC
                        IF ( IPROC.EQ.NPROCS ) THEN
                           IPROC = NPROCS - 1
                           INC = -1
                        ELSE IF ( IPROC.EQ.-1 ) THEN
                           IPROC = 0
                           INC = 1
                        END IF
                        IPROCE(IE) = IPROC
                     END DO
                  END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute E-points IE = 1, ..., NETAB(IEPATH)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C      CALL MPI_DISTRIBUTE(IPROCE,NETAB(IEPATH),MPI_ELOOP,'E')
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
                  PHASA(:,:,:) = C0
                  PHAST(:,:,:) = C0
                  PHASK(:) = C0
                  CTOTDOS(:) = C0
C
C=======================================================================
C            initialize arrays used for an energy path
C=======================================================================
C
C------------------------------------------------------- total densities
                  RHO2NSX(:,:,:,:) = 0.0D0
                  RHOCHRX(:,:) = 0.0D0
                  RHOSPNX(:,:) = 0.0D0
                  RHOORBX(:,:) = 0.0D0
C
C--------------------------------- moments, hyperfine field, band energy
C
                  OBS_TX(:,:,:) = 0D0
                  OBS_LTX(:,:,:,:) = 0D0
                  CMNTTX(:,:) = 0D0
C
c modified by XJQ: parallel on energies using my own subroutines
                  call get_comm_level(1,'parent',lparalleled,
     &              parent_comm,nprocs_parent,parent_rank)
                  call mpi_barrier(parent_comm,mpierr)
                  call mpi_multilevel_distribute(routine,1,
     &              netab(iepath),lparalleled,ie0,ie1)
c end-mod-xjq
                  LOOP_IE:DO IE = 1,NETAB(IEPATH)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ: parallel on energies using my own subroutines
c                  collecting to root is safe, since inter_comm having root
c                  always has all values of energy
c                     IF ( MPI_ID.NE.IPROCE(IE) ) CYCLE LOOP_IE
                     if(ie<ie0 .or. ie>ie1) cycle loop_ie
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                     IF ( DEBUG ) WRITE (6,99026) IEPANEL,IEPATH,IE,ERYD
C
                     IECURR = IE
                     ERYD = ETAB(IE,IEPATH)
C
                     ICPAFLAG = 0
                     CPACHNG = 0.0D0
C
C ======================================================================
C                       solve SS - differential equation
C ======================================================================
                     IWRREGWF = 1
                     IWRIRRWF = 1
                     GETIRRSOL = .TRUE.
C
                     WAVE_FUNCTIONS_AVAILABLE = GETIRRSOL .AND. 
     &                  (IWRREGWF.EQ.1) .AND. (IWRIRRWF.EQ.1)
C
                     IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6)
     &                    .EQ.'I-ZONE' ) THEN
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
                     IF ( DMFT .OR. LDAU )
     &                    DMFTSIG(1:NKM,1:NKM,ITBOT:ITTOP)
     &                    = DMFTSIGMA(1:NKM,1:NKM,ITBOT:ITTOP,IE)
C
                     CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,
     &                             GETIRRSOL,ERYD,P,IPRINT,TSST,MSST,
     &                             SSST,MEZZ,MEZJ,ORBPOL)
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
                     IF ( THERMAL_VIBRA_FLUCT )
     &                    CALL THERMAL_INIT_UFMAT(ERYD)
c modified by XJQ: scf of vibrations
                     if(lscfvb) then
                       thermal_vibra_fluct = .true.
                       CALL THERMAL_INIT_UFMAT(ERYD)
                       thermal_vibra_fluct = .false.
                     endif
c end-mod-xjq
C
C ======================================================================
C
C ------------------------------------------------------- initialize m_q
C
                     CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
                     IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6)
     &                    .EQ.'I-ZONE' ) THEN
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
C ----------------------------------------------------------- IEPATH = 1
                     IF ( IEPATH.EQ.1 ) THEN
C
C -------------------------------------------------- do BZ - integration
C
                        CALL TAU_DRIVE(IECURR,IPRINTBAND,ERYD,P,TSSQ,
     &                                 MSSQ,TSST,MSST,TAUQ,ICPAFLAG,
     &                                 CPACHNG,ITCPA,ICPACONV,PHASK)
C
C ----------------------------------------------------------------------
C
                        IF ( ICPAFLAG.NE.0 ) THEN
                           NCPAFAIL = NCPAFAIL + 1
                           CPACHTAB(NCPAFAIL) = CPACHNG
                           IECPAFAIL(NCPAFAIL) = IECURR
                        END IF
C
C -------------------- alloys dealt by ATA (ITCPAMAX=0) : reset m-matrix
C
                        IF ( NCPA.NE.0 .AND. ITCPAMAX.EQ.0 )
     &                       CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ----------------------------------------------- project TAU on type IT
C
                        CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,
     &                               TAUQ,TAUT)
C
                        IF ( IPRINT.GE.5 .OR. WRMAT )
     &                       CALL DUMPTAU(IECURR,ERYD,6,MSST,MSSQ,TAUT,
     &                       TAUQ)
C
                     END IF
C ----------------------------------------------------------- IEPATH = 1
C
C ======================================================================
C                     deal with Lloyd's formula
C ======================================================================
C
                     IF ( LLOYD ) CALL LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,
     &                    MSST,SSST,ERYD,IE)
C
C ======================================================================
C  if the epath is split into BS- and SS-part manipulate TAU accordingly
C
                     IF ( SPLITSS ) THEN
                        M = NKMMAX
                        IF ( IEPATH.EQ.1 ) THEN
                           TAUT(1:M,1:M,ITBOT:ITTOP)
     &                        = TAUT(1:M,1:M,ITBOT:ITTOP)
     &                        - TSST(1:M,1:M,ITBOT:ITTOP)
                        ELSE
                           TAUT(1:M,1:M,ITBOT:ITTOP)
     &                        = TSST(1:M,1:M,ITBOT:ITTOP)
                        END IF
                     END IF
C
C ======================================================================
C                         set up charge density
C ======================================================================
C
                     CALL CHRDNS(.FALSE.,CTOTDOS,SHFTEF,TOTDOSX,TOTDOS,
     &                           TOTNOSX,ERYD,WETAB(IE,IEPATH),RDUMARG,
     &                           IEPANEL,IEPATH,IECURR,OBS_LTX,OBS_TX,
     &                           BCOR,BCORS,MEZZ,MEZJ,TSST,MSST,TAUT,
     &                           MSSQ,TAUQ,RHOCHRX,RHOSPNX,RHOORBX,
     &                           RHO2NSX,CMNTTX)
C
C---------------------------------------------------------------- BROOKS
C
                     IF ( ORBPOL(1:6).EQ.'BROOKS' )
     &                    CALL SCFOPPOT(.FALSE.,SHFTEF,IECURR,OBS_T,
     &                    TAUT,RHOOPC,RHOOPO,QLZ,AOPT)
C
                     IF ( IEPATH.EQ.NEPATH .AND. IE.EQ.NETAB(IEPATH)
     &                    .AND. BLCOUPL ) THEN
C
                        IF ( ORBPOL(1:6).EQ.'NONE' )
     &                       CALL RINIT(NRMAX*2*NTMAX,AOPT)
C
                        DO IT = ITBOT,ITTOP
                           JTOP = JRWS(IMT(IT))
                           DO IS = 1,2
                              DO I = 1,JTOP
                                 AOPT(I,IS,IT) = AOPT(I,IS,IT) - BEXT
                              END DO
                           END DO
                        END DO
                     END IF
C
C----------------------------------------------------------------- BREIT
                     IF ( BREITINT ) THEN
C
                        CALL CURDNS(.FALSE.,RENORMALIZE,SHFTEF,SCLNOS,
     &                              JDNST,WETAB(IECURR,IEPATH),TAUT)
C
                        CALL SCFBIPOT(.FALSE.,1,SHFTEF,ABITNEW,
     &                                WETAB(IECURR,IEPATH),TAUT,TSST)
C
                     END IF
C
C=======================================================================
C                     Green's function matrix
C=======================================================================
C
                     IF ( DMFT .OR. LDAU )
     &                    CALL GFUNMAT_DRIVE(ERYD,IE,TAUT,MSST,GFMAT)
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
                  END DO LOOP_IE
C-------------------------------------------------------------------- IE
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C             collect results for energy points of a E-path
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c                  IF ( NPROCS>1 ) THEN
c modified by XJQ: parallel on energies using my own subroutines
                  if(lparalleled) then
                     call mpi_barrier(parent_comm,mpierr)
c using inter_comm will allow nprocs > ne
                     call get_comm_level(1,'inter ',lparalleled,
     &                 inter_comm,nprocs_inter,inter_rank)
c colloct ie0, ie1 and ne_mpi for further collecting
                     allocate(disp_proc(0:nprocs_inter-1))
                     allocate(ne_proc(0:nprocs_inter-1))
                     ne_mpi = ie1-ie0+1
                     call collect_varstartnvar(1,inter_comm,nprocs_inter
     &                 ,ie0,ne_mpi,disp_proc(0:),ne_proc(0:))
C
                     call mpi_barrier(parent_comm,mpierr)
                     IF ( FULLPOT ) THEN
                        NVEC = 3
                        IDIMS = NRMAX
                        DO IVEC = 1,NVEC
                           DO IT = ITBOT,ITTOP
                              DO LM = 1,NLMFPT(IT)
C
c                                 CALL DRV_MPI_REDUCE_R
c     &                              (RHO2NSX(1,LM,IT,IVEC),RWORK(1,1),
c     &                              IDIMS)
c XJQ: maybe later I should write my own mpi_reduce interface with communicator input
                                 call mpi_reduce(RHO2NSX(1,LM,IT,IVEC),
     &                             rwork(1,1),idims,mpi_double_precision
     &                             ,mpi_sum,0,inter_comm,mpierr)
                                 call mpi_barrier(parent_comm,mpierr)
                                 if(inter_rank==0) 
     &                             RHO2NSX(1:idims,LM,IT,IVEC)=
     &                               rwork(1:idims,1)
C
                              END DO
                           END DO
                        END DO
                     ELSE
                        IDIMS = NRMAX*NTMAX
c                        CALL DRV_MPI_REDUCE_R(RHOCHRX(1,1),RWORK(1,1),
c     &                     IDIMS)
c                        CALL DRV_MPI_REDUCE_R(RHOSPNX(1,1),RWORK(1,1),
c     &                     IDIMS)
c                        CALL DRV_MPI_REDUCE_R(RHOORBX(1,1),RWORK(1,1),
c     &                     IDIMS)
                        call mpi_reduce(rhochrx(1,1),rwork(1,1),idims,
     &                    mpi_double_precision,mpi_sum,
     &                    0,inter_comm,mpierr)
                        call mpi_barrier(parent_comm,mpierr)
                        if(inter_rank==0) rhochrx(1:nrmax,1:ntmax)=
     &                    rwork(1:nrmax,1:ntmax)
                        call mpi_reduce(rhospnx(1,1),rwork(1,1),idims,
     &                    mpi_double_precision,mpi_sum,
     &                    0,inter_comm,mpierr)
                        call mpi_barrier(parent_comm,mpierr)
                        if(inter_rank==0) rhospnx(1:nrmax,1:ntmax)=
     &                    rwork(1:nrmax,1:ntmax)
                        call mpi_reduce(rhoorbx(1,1),rwork(1,1),idims,
     &                    mpi_double_precision,mpi_sum,
     &                    0,inter_comm,mpierr)
                        call mpi_barrier(parent_comm,mpierr)
                        if(inter_rank==0) rhoorbx(1:nrmax,1:ntmax)=
     &                    rwork(1:nrmax,1:ntmax)
                     END IF
C
C---------------------------------------------------------------- BROOKS
                     IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
                        DO IS = 1,2
c                           CALL DRV_MPI_REDUCE_R(RHOOPC(1,1,IS),
c     &                        RWORK(1,1),IDIMS)
c                           CALL DRV_MPI_REDUCE_R(RHOOPO(1,1,IS),
c     &                        RWORK(1,1),IDIMS)
                           call mpi_reduce(rhoopc(1,1,is),rwork(1,1),
     &                       nrmax*ntmax,mpi_double_precision,mpi_sum,
     &                       0,inter_comm,mpierr)
                           call mpi_barrier(parent_comm,mpierr)
                           if(inter_rank==0) rhoopo(1:nrmax,1:ntmax,is)=
     &                       rwork(1:nrmax,1:ntmax)
                           call mpi_reduce(rhoopc(1,1,is),rwork(1,1),
     &                       nrmax*ntmax,mpi_double_precision,mpi_sum,
     &                       0,inter_comm,mpierr)
                           call mpi_barrier(parent_comm,mpierr)
                           if(inter_rank==0) rhoopo(1:nrmax,1:ntmax,is)=
     &                       rwork(1:nrmax,1:ntmax)
                        END DO
                     END IF
C
C----------------------------------------------------------------- BREIT
c XJQ: change mpi_reduce later
                     IF ( BREITINT ) THEN
                        IDIMS = NRMAX
                        DO ILA = 1,NLABIMAX
                           DO IVEC = -1,1
                              DO IT = ITBOT,ITTOP
C
                                 CALL DRV_MPI_REDUCE_R
     &                              (ABITNEW(1,ILA,IVEC,IT),RWORK(1,1),
     &                              IDIMS)
C
                                 CALL DRV_MPI_REDUCE_R
     &                              (JDNST(1,ILA,IVEC,IT),RWORK(1,1),
     &                              IDIMS)
C
                              END DO
                           END DO
                        END DO
                     END IF
C
C-----------------------------------------------------------------------
C
                     LWK = 4*NOBSMAX*NLMAX*NTMAX
c                     ALLOCATE (RWK(LWK))
                     allocate(array4d_real(0:3,1:NOBSMAX,
     &                 1:NLMAX,1:ntmax))
c                     CALL DRV_MPI_REDUCE_R(OBS_LTX(0,1,1,1),RWK(1),LWK)
                     call mpi_reduce(OBS_LTX(0,1,1,1),
     &                 array4d_real(0,1,1,1),LWK,
     &                 mpi_double_precision,mpi_sum,0,inter_comm,mpierr)
                     call mpi_barrier(parent_comm,mpierr)
                     if(inter_rank==0) obs_ltx(0:,1:,1:,1:)=
     &                 array4d_real(0:,1:,1:,1:)
                     deallocate(array4d_real)
                     
c                     DEALLOCATE (RWK)
C
                     LWK = 4*NOBSMAX*NTMAX
c                     ALLOCATE (RWK(LWK))
                     allocate(array3d_real(0:3,1:NOBSMAX,
     &                 1:ntmax))
c                     CALL DRV_MPI_REDUCE_R(OBS_TX(0,1,1),RWK(1),LWK)
                     call mpi_reduce(OBS_TX(0,1,1),
     &                 array3d_real(0,1,1),LWK,
     &                 mpi_double_precision,mpi_sum,0,inter_comm,mpierr)
                     call mpi_barrier(parent_comm,mpierr)
                     if(inter_rank==0) obs_tx(0:,1:,1:)=
     &                 array3d_real(0:,1:,1:)
                     deallocate(array3d_real)
c                     DEALLOCATE (RWK)
C
                     LWKLMT = NLMFPMAX*NTMAX
c                     CALL DRV_MPI_REDUCE_R(CMNTTX(1,1),RWKLMT(1,1),
c     &                  LWKLMT)
                     call mpi_reduce(CMNTTX(1,1),RWKLMT(1,1),LWKLMT,
     &                 mpi_double_precision,mpi_sum,0,inter_comm,mpierr)
                     call mpi_barrier(parent_comm,mpierr)
                     if(inter_rank==0) cmnttx(1:,1:)=rwklmt(1:,1:)
C
c XJQ: change mpi_reduce later
                     IF ( LLOYD ) THEN
                        LWKE = NETAB(IEPATH)
C
                        CALL DRV_MPI_REDUCE_C(CTOTDOS(1),CWKE(1),LWKE)
                        CALL DRV_MPI_REDUCE_C(PHASK(1),CWKE(1),LWKE)
C
                        LWKKMTE = NKMMAX*NTMAX*NETAB(IEPATH)
C
                        CALL DRV_MPI_REDUCE_C(PHASA(1,1,1),
     &                     CWKKMTE(1,1,1),LWKKMTE)
                        CALL DRV_MPI_REDUCE_C(PHAST(1,1,1),
     &                     CWKKMTE(1,1,1),LWKKMTE)
                     END IF
C
C------------------------------------------------------------------ DMFT
                     IF ( DMFT .OR. LDAU ) THEN
c mpi_reduce is very slow when ntmax is large, replace it by mpi_gatherv
c                        DEALLOCATE (CWKT)
c                        ALLOCATE (CWKT(NKMMAX))
c                        LWKKMTE = NKMMAX
c                        DO IVEC = 1,NKMMAX
c                           DO IE = 1,NEMAX
c                              DO IT = 1,NTMAX
c                                 CALL DRV_MPI_REDUCE_C
c     &                              (GFMAT(1,IVEC,IE,IT),CWKT(1),
c     &                              LWKKMTE)
c                              END DO
c                           END DO
c                        END DO
c                        DEALLOCATE (CWKT)
c                        ALLOCATE (CWKT(NTMAX))
                        allocate(gfmat_tmp(nkmmax,nkmmax,nemax))
                        call mpi_type_contiguous(nkmmax*nkmmax,
     &                    mpi_double_complex,mpi_row,mpierr)
                        call mpi_type_commit(mpi_row,mpierr)
                        do it=1,ntmax
                          call mpi_gatherv(gfmat(1,1,ie0,it),ne_mpi,
     &                      mpi_row,gfmat_tmp(1,1,1),
     &                      ne_proc(0:),disp_proc(0:),
     &                      mpi_row,0,inter_comm,mpierr)
                          call mpi_barrier(parent_comm,mpierr)
                          if(inter_rank==0) 
     &                      gfmat(1:nkmmax,1:nkmmax,1:nemax,it)=
     &                      gfmat_tmp(1:nkmmax,1:nkmmax,1:nemax)
                        enddo
                        call mpi_type_free(mpi_row,mpierr)
                        call mpi_barrier(parent_comm,mpierr)
                        deallocate(gfmat_tmp)
                     END IF
C
C----------------------------------------------------------------- BREIT
c XJQ: change mpi_reduce later
                     IF ( BREITINT ) THEN
                        IDIMS = NRMAX
                        DO ILA = 1,NLABIMAX
                           DO IVEC = -1,1
                              DO IT = ITBOT,ITTOP
                                 CALL DRV_MPI_REDUCE_R
     &                              (ABITNEW(1,ILA,IVEC,IT),RWORK(1,1),
     &                              IDIMS)
                              END DO
                           END DO
                        END DO
                     END IF
C
                     call mpi_barrier(parent_comm,mpierr)
c free mpi_distribution on iqqp
                     call mpi_multilevel_free(1)
                     deallocate(disp_proc)
                     deallocate(ne_proc)
c end-mod-xjq
                  END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C ======================================================================
C  renormalize GF-related quantities like charge density if necessary
C ======================================================================
                  IF ( MPI_ID.EQ.0 ) THEN
C
C------------------------------------------ perform corrections for E=EF
C                 needs from run for IE=NETAB:  TOTDOS, TAU- and WF-data
C                                               TOTNOS is recalculated
                     OBS_X(:,:) = 0.0D0
                     TOTNOSX = 0.0D0
                     DO IT = ITBOT,ITTOP
                        RAUX = 0D0
                        DO IL = 1,NL
                           RAUX = RAUX + OBS_LTX(0,IDOS,IL,IT)
                        END DO
                        OBS_X(:,:) = OBS_X(:,:) + OBS_TX(:,:,IT)
     &                               *CONC(IT)*NAT(IT)
                        TOTNOSX = TOTNOSX + RAUX*CONC(IT)*NAT(IT)
                     END DO
C
                     IF ( IEPANEL.EQ.NEPANEL ) THEN
                        TOTDOSEFX = DIMAG(CTOTDOS(NETAB(IEPATH)))
                        TOTDOSEF = TOTDOSEF + TOTDOSEFX
                     END IF
                     IF ( SCF_CHECK_SPLITSS ) WRITE (6,99025) TOTNOS,
     &                    TOTNOSX,TOTDOSEF,TOTDOSEFX
C
                     RENORMALIZE = .FALSE.
C
C ---------------------------------------------------------------- LLOYD
                     IF ( LLOYD ) THEN
C
                        CALL LLOYDPT2(TOTNOSX,CTOTDOS,SCLNOS,EBANDLD,
     &                                EFERMILD,SCFMIX,NCPAFAIL,
     &                                IECPAFAIL,PHASK,PHAST,PHASA)
C
                        RENORMALIZE = .TRUE.
C
                     END IF
C
                     CALL CHRDNS_NORM(RENORMALIZE,SCLNOS,TOTDOSEF,
     &                                TOTNOS,TOTNOSX,IEPANEL,IEPATH,
     &                                IECURR,OBS_LT,OBS_T,OBS_LTX,
     &                                OBS_TX,RHOCHRX,RHOSPNX,RHOORBX,
     &                                RHO2NSX,CMNTTX)
C
                  END IF
C ======================================================================
C
               END DO
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP EPATH
C
            END DO
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP EPANEL
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C
C ======================================================================
C                apply corrections for the Fermi energy
C ======================================================================
C
            IF ( MPI_ID.EQ.0 ) THEN
C
C--------------------------------- correct at E_F and add core densities
C
               CALL CHRDNS(.TRUE.,CTOTDOS,SHFTEF,TOTDOSEFX,TOTDOSEF,
     &                     TOTNOS,ERYD,C0,EFERMILD,IEPANEL,IEPATH,
     &                     IECURR,OBS_LT,OBS_T,BCOR,BCORS,MEZZ,MEZJ,
     &                     TSST,MSST,TAUT,MSSQ,TAUQ,RHOCHRX,RHOSPNX,
     &                     RHOORBX,RHO2NSX,CMNTTX)
C
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
C  Makes sure that density is always positive
C     For full potential, deal with the spherical part only
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
               IF ( FULLPOT ) THEN
                  DO IT = ITBOT,ITTOP
                     DO IR = 1,JRWS(IMT(IT))
                        IF ( RHO2NS(IR,1,IT,1).LT.0.0D0 ) THEN
                           IF ( Z(IT).NE.0 .AND. IPRINT.GE.0 )
     &                           WRITE (6,'(a,a,i5,i7,f15.8,1pe17.5)')
     &                           'Warning: Negative density for non-',
     &                          'vacuum atom:   IT, IR, R, RHOCHR =',IT,
     &                          IR,R(IR,IMT(IT)),RHO2NS(IR,1,IT,1)
                           RHO2NS(IR,1,IT,1) = 0.0D0
                           RHO2NS(IR,1,IT,2) = 0.0D0
                        ELSE IF ( ABS(RHO2NS(IR,1,IT,2))
     &                            .GE.RHO2NS(IR,1,IT,1) ) THEN
                           IF ( Z(IT).NE.0 .AND. IPRINT.GE.0 )
     &                           WRITE (6,'(a,a,i5,i7,f15.8,2(1pe17.5))'
     &                          ) 'Warning:  |SPN| > CHR for non-vacuum'
     &                            ,
     &                            ' atom:   IT, IR, R, RHOCHR, RHOSPN ='
     &                            ,IT,IR,R(IR,IMT(IT)),RHO2NS(IR,1,IT,1)
     &                            ,RHO2NS(IR,1,IT,2)
                           RHO2NS(IR,1,IT,2)
     &                        = SIGN(RHO2NS(IR,1,IT,1),RHO2NS(IR,1,IT,2)
     &                        )
                        END IF
                     END DO
                  END DO
               ELSE
                  DO IT = ITBOT,ITTOP
                     DO IR = 1,JRWS(IMT(IT))
                        IF ( RHOCHR(IR,IT).LT.0.0D0 ) THEN
                           IF ( Z(IT).NE.0 .AND. IPRINT.GE.0 )
     &                           WRITE (6,'(a,a,i5,i7,f15.8,1pe17.5)')
     &                           'Warning: Negative density for non-',
     &                          'vacuum atom:   IT, IR, R, RHOCHR =',IT,
     &                          IR,R(IR,IMT(IT)),RHOCHR(IR,IT)
                           RHOCHR(IR,IT) = 0.0D0
                           RHOSPN(IR,IT) = 0.0D0
                           RHOORB(IR,IT) = 0.0D0
                        ELSE IF ( ABS(RHOSPN(IR,IT)).GT.RHOCHR(IR,IT) )
     &                            THEN
                           IF ( Z(IT).NE.0 .AND. IPRINT.GE.0 )
     &                           WRITE (6,'(a,a,i5,i7,f15.8,2(1pe17.5))'
     &                          ) 'Warning:  |SPN| > CHR for non-vacuum'
     &                            ,
     &                            ' atom:   IT, IR, R, RHOCHR, RHOSPN ='
     &                            ,IT,IR,R(IR,IMT(IT)),RHOCHR(IR,IT),
     &                            RHOSPN(IR,IT)
                           RHOSPN(IR,IT)
     &                        = SIGN(RHOCHR(IR,IT),RHOSPN(IR,IT))
                        END IF
                     END DO
                  END DO
               END IF
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
C
C --------------------------------------------------------------- BROOKS
               IF ( ORBPOL(1:6).EQ.'BROOKS' )
     &              CALL SCFOPPOT(.TRUE.,SHFTEF,IECURR,OBS_T,TAUT,
     &              RHOOPC,RHOOPO,QLZ,AOPT)
C
C----------------------------------------------------------------- BREIT
               IF ( BREITINT ) CALL CURDNS(.TRUE.,RENORMALIZE,SHFTEF,
     &              SCLNOS,JDNST,C0,TAUT)
C
C-------------------------- use Femri energy obtained from Lloyd formula
C
               IF ( LLOYD ) THEN
C               EFERMI = EFERMILD
C               SHFTEF = EFERMILD - EFERMI0
C
C              SHFTEF = 0.5D0*(EFERMILD - EFERMI0)
C
C              SHFTEF = EFERMILD - EFERMI0
C              EFERMI = EFERMI0 + SHFTEF
C              EFERMILD = EFERMI
               END IF
C
               IF ( ABS(SHFTEF).GT.SHFTEFLIM ) THEN
                  EFERMI = EFERMI0 + SIGN(SHFTEFLIM,SHFTEF)
                  SHFTEF = SIGN(SHFTEFLIM,SHFTEF)
                  WRITE (6,99019) SHFTEF,EFERMI
               END IF
C
C-------------------------------------- use Fermi energy fixed by L-BULK
C
C
               IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SYSTEM_TYPE(1:3)
     &              .NE.'VIV' ) THEN
                  IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) EFERMI = EFERMI0
               END IF
C ----------------------------------------------------------------------
C                        Breit interaction
C ----------------------------------------------------------------------
C
               IF ( BREITINT ) THEN
C
                  CALL INPUT_FIND_SECTION('MODE',0)
                  NESS = 100
                  IF ( FOUND_SECTION )
     &                 CALL SECTION_SET_INTEGER('BRNESS',NESS,9999,0)
C
                  IF ( .NOT.ALLOCATED(ETABSS_BI) )
     &                 ALLOCATE (ETABSS_BI(NESS),WETABSS_BI(NESS))
C
                  CALL EPATH(3,EMIN,EFERMI0,0.0D0,NESS,ETABSS_BI,
     &                       WETABSS_BI,EILOW,IPRINT,NESS)
C
                  DO IE = 1,NESS
C
                     ERYD = ETABSS_BI(IE)
C
                     WRITE (6,*) ' SS-BIPOT  IE,E ',IE,ERYD
C
                     CALL SSITE(1,1,IFILCBWF,CALCINT,.TRUE.,ERYD,P,
     &                          IPRINT,NKM,TSST,MSST,SSST,MEZZ,MEZJ,
     &                          ORBPOL)
C
                     CALL SCFBIPOT(.FALSE.,2,0.0D0,ABITNEW,
     &                             WETABSS_BI(IE),TAUT,TSST)
C
                  END DO
C
                  CALL SCFBIPOT(.TRUE.,2,SHFTEF,ABITNEW,WETAB(1,1),TAUT,
     &                          TSST)
C
C
                  ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX)
     &               = 2.0D0*ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX)
C
C              ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0.0D0
                  WRITE (*,*) 'on-site j*A'
                  CALL ETOTAVEC(ABITNEW)
C
C              ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0.0D0
                  CALL SCFBIPOTMAD(ABITNEW)
C
                  WRITE (*,*) 'on+off-site j*A'
                  CALL ETOTAVEC(ABITNEW)
C
                  DO IT = ITBOT,ITTOP
                     DO M = -1, + 1
                        DO ILA = 1,NLABIMAX
                           DO I = 1,NRMAX
                              ABIT(I,ILA,M,IT) = ABITNEW(I,ILA,M,IT)
                           END DO
                        END DO
                     END DO
                  END DO
C
                  ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0.0D0
                  CALL SCFBIPOTMAD(ABITNEW)
                  WRITE (*,*) 'offsite j*A'
                  CALL ETOTAVEC(ABITNEW)
C
                  JDNST(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0D0
                  ABITNEW(1:NRMAX,1:NLABIMAX,-1:+1,1:NTMAX) = 0.0D0
               END IF
C
C-----------------------------------------------------------------------
               IF ( MIXRHO ) THEN
                  IF ( FULLPOT ) THEN
C
                     W1MIX(1:NRMAX,1:NTMAX)
     &                  = RHO2NS(1:NRMAX,1,1:NTMAX,1)
                     W2MIX(1:NRMAX,1:NTMAX)
     &                  = RHO2NS(1:NRMAX,1,1:NTMAX,2)
C
                     W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &                  = RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1)
                     W2NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &                  = RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,2)
                  ELSE
C
                     W1MIX(1:NRMAX,1:NTMAX) = RHOCHR(1:NRMAX,1:NTMAX)
                     W2MIX(1:NRMAX,1:NTMAX) = RHOSPN(1:NRMAX,1:NTMAX)
C
                  END IF
                  IF ( FULLPOT ) THEN
C JM: Only after first iteration rho2ns is consistent
C                --- solution is to write RHO2NS into potential file ? ---
                     IF ( ITRSCF.EQ.1 )
     &                    CALL SCFBROYPT1(IPRINT,0,IMIX,ISTBRY,ITDEPT,
     &                    SCFMIX,W1MIX,W2MIX,W1NSMIX,W2NSMIX,RMSAVV,
     &                    RMSAVB,JRNS1TMP,JRNSMINTMP)
C
                     CALL SCFBROYPT1(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,
     &                               SCFMIX,W1MIX,W2MIX,W1NSMIX,W2NSMIX,
     &                               RMSAVV,RMSAVB,JRNS1TMP,JRNSMINTMP)
C
                     IF ( ITRSCF.EQ.1 ) THEN
                        WRITE (6,*) 
     &                          'WARNING: IF FIRST ITERATION AND MIXRHO'
                        WRITE (6,*) '         RMS ERROR SET TO 1E08'
                        RMSAVV = 9999999.9D0
                     END IF
C
C
                  ELSE IF ( SCFALG(1:8).EQ.'BROYDEN2' .OR. SCFALG(1:8)
     &                      .EQ.'ANDERSON' ) THEN
                     CALL SCFBROYDEN(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,
     &                               SCFMIX,R,R2DRDI,JRWS,IMT,W1MIX,
     &                               W2MIX,RMSAVV,RMSAVB,NMMAX,NRMAX,
     &                               NTMAX)
                  ELSE IF ( SCFALG(1:6).EQ.'TCHEBY' ) THEN
                     CALL SCFTCHEBY(ITRSCF,VT,BT,SCFMIX,RMSAVV,RMSAVB)
                  ELSE
                     CALL STOP_MESSAGE(ROUTINE,'SCFALG not found')
C
                  END IF
C
                  IF ( FULLPOT ) THEN
                     RHO2NS(1:NRMAX,1,1:NTMAX,1)
     &                  = W1MIX(1:NRMAX,1:NTMAX)
C
                     RHO2NS(1:NRMAX,1,1:NTMAX,2)
     &                  = W2MIX(1:NRMAX,1:NTMAX)
C
                     RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1)
     &                  = W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
C
                     RHO2NS(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,2)
     &                  = W2NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
                  ELSE
                     RHOCHR(1:NRMAX,1:NTMAX) = W1MIX(1:NRMAX,1:NTMAX)
                     RHOSPN(1:NRMAX,1:NTMAX) = W2MIX(1:NRMAX,1:NTMAX)
C
                     DO IT = ITBOT,ITTOP
                        IM = IMT(IT)
C
                        IF ( IREL.LE.1 ) THEN
                           DO IR = 1,JRWS(IM)
                              RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
                           END DO
                        ELSE
                           DO IR = 1,JRWS(IM)
                              RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
                           END DO
                        END IF
C
                        CALL RRADINT(IM,RINT,QEL(IT))
                     END DO
                  END IF
               END IF
C-----------------------------------------------------------------------
C
               IF ( FULLPOT ) THEN
C
                  CALL FPSCFNEWPOT(ETOT,ITRSCF,RHO2NS)
C
                  IF ( NONMAG ) THEN
                     CALL RINIT(NRMAX*NTMAX,BT)
                     CALL RINIT((NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX,BNST)
                     IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL )
     &                    CALL RINIT(NRMAX*2*NTMAX,AOPT)
                  END IF
C
C
C-----------------------------------------------------------------------
C   suppress aspherical potential in case of FULLCHARGE calculations
C-----------------------------------------------------------------------
C
                  IF ( FULLCHARGE ) THEN
                     VNST(:,:,:) = 0D0
                     BNST(:,:,:) = 0D0
                     WRITE (6,99027)
                  END IF
C
               ELSE
C
                  CALL SCFNEWPOT(.TRUE.,ETOT,ECORTAB)
C
                  IF ( NONMAG ) THEN
                     CALL RINIT(NRMAX*NTMAX,BT)
                     IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL )
     &                    CALL RINIT(NRMAX*2*NTMAX,AOPT)
                  END IF
C
                  IF ( MIXRHO .AND. KMOL.EQ.0 ) CALL VMUFTIN(VMTZ)
               END IF
C
C-----------------------------------------------------------------------
C
               IF ( PLOTPRS(0) ) CALL FPPLOTPRS
C
C ----------------------------------------------------------------------
               IF ( MAGCOUPL.NE.'NO' ) THEN
C
                  DO IT = ITBOT,ITTOP
                     IM = IMT(IT)
                     JTOP = JRWS(IM)
                     DO I = 1,JTOP
                        RINT(I) = R2DRDI(I,IM)*BT(I,IT)
                     END DO
                     CALL RRADINT(IM,RINT,BINTT)
                     IF ( IT.EQ.1 ) THEN
                        BINT1 = BINTT
                     ELSE IF ( (MAGCOUPL.EQ.'FM' .AND. (BINT1*BINTT.LT.
     &                         0D0)) .OR. 
     &                         (MAGCOUPL.EQ.'AF' .AND. (BINT1*BINTT.GT.
     &                         0D0)) ) THEN
                        WRITE (6,99008) MAGCOUPL,IT
                        BSCAL = -1.0D0
                        CALL DSCAL(JTOP,BSCAL,BT(1,IT),1)
                     END IF
C
                  END DO
               END IF
C
C------------------------------------------------------ mix AT -- BROOKS
               RMSAVA = 0D0
               IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
                  WW1 = SCFMIXOP
                  WW2 = 1.D0 - SCFMIXOP
                  DO IT = ITBOT,ITTOP
                     IM = IMT(IT)
                     JTOP = JRWS(IM)
                     RMSERA = 0D0
                     DO IS = 1,2
                        DO I = 1,JTOP
                           DEL = AOPT(I,IS,IT) - AOPT0(I,IS,IT)
                           RINT(I) = DEL*DEL*R2DRDI(I,IM)
                           AOPT(I,IS,IT) = WW1*AOPT(I,IS,IT)
     &                        + WW2*AOPT0(I,IS,IT)
                           AOPT0(I,IS,IT) = AOPT(I,IS,IT)
                        END DO
                        CALL RRADINT(IM,RINT,RMSERAS)
                        RMSERA = MAX(RMSERA,SQRT(RMSERAS))
                     END DO
                     RMSAVA = MAX(RMSAVA,RMSERA)
                     WRITE (6,99016) IT,RMSERA
                  END DO
               END IF
C
               OPDONE = .TRUE.
C
C-----------------------------------------------------------------------
               IF ( .NOT.MIXRHO ) THEN
                  IF ( .NOT.FULLPOT ) THEN
C
                     IF ( KMOL.EQ.0 ) CALL VMUFTIN(VMTZ)
C
                     IF ( SCFALG(1:8).EQ.'BROYDEN2' .OR. SCFALG(1:8)
     &                    .EQ.'ANDERSON' ) THEN
                        CALL SCFBROYDEN(IPRINT,ITRSCF,IMIX,ISTBRY,
     &                                  ITDEPT,SCFMIX,R,R2DRDI,JRWS,IMT,
     &                                  VT,BT,RMSAVV,RMSAVB,NMMAX,NRMAX,
     &                                  NTMAX)
                     ELSE IF ( SCFALG(1:6).EQ.'TCHEBY' ) THEN
                        CALL SCFTCHEBY(ITRSCF,VT,BT,SCFMIX,RMSAVV,
     &                                 RMSAVB)
                     ELSE
                        CALL STOP_MESSAGE(ROUTINE,'SCFALG not found')
                     END IF
C
                  ELSE
C
                     CALL SCFBROYPT1(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,
     &                               SCFMIX,VT,BT,VNST,BNST,RMSAVV,
     &                               RMSAVB,JRNS1,JRNSMIN)
C
                  END IF
C
               END IF
C
               IF ( IPRINT.GT.0 ) THEN
                  DO IT = ITBOT,ITTOP
                     STRINP = 'V_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &                        (1:LTXT_T(IT))//'_I'
                     CALL STRING_ADD_N(STRINP,ITRSCF)
                     LL = LEN_TRIM(STRINP)
                     OPEN (80,FILE=STRINP(1:LL))
                     WRITE (80,'(''#  ITERATION '',I3,5X,A)') ITRSCF,
     &                      TXT_T(IT)
                     IM = IMT(IT)
                     JTOP = JRWS(IM)
                     DO IR = 1,JTOP
                        WRITE (80,'(2E15.6)') R(IR,IM),VT(IR,IT)
     &                         *R(IR,IM)
                     END DO
                     CLOSE (80)
                  END DO
               END IF
C-----------------------------------------------------------------------
C
C            CALL SCFBROYANG(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,SCFMIX,
C     &                      ITOQ,DROTQ,MVPHI,MVTET,MVGAM,QMPHI,
C     &                      QMTET,QMGAM,ANGFM,ANGSM,ANGFM1,ANGSM1,
C     &                      ANGINP,ANGOUT,NQHOST,NK,ERRAVANG,NQMAX,
C     &                      NTMAX,NMVECMAX,NKMMAX)
C
               IF ( ITERMVDIR ) CALL SCFITERANG(ITRSCF,ITOQ,DROTQ,MVPHI,
     &              MVTET,MVGAM,QMPHI,QMTET,QMGAM,NQHOST,NK,ERRAVANG,
     &              NQMAX,NTMAX,NMVECMAX,NKMMAX)
C
               CALL CPU_TIME(CPU_TIME_SCF_ITER_END)
               WRITE (6,99010) CPU_TIME_SCF_ITER_END - 
     &                         CPU_TIME_SCF_ITER_START
C
               IF ( ABS(RMSAVV).GT.1D-12 ) THEN
                  LOG10RMSV = LOG10(RMSAVV)
               ELSE
                  LOG10RMSV = 0D0
               END IF
               IF ( ABS(RMSAVB).GT.1D-12 ) THEN
                  LOG10RMSB = LOG10(RMSAVB)
               ELSE
                  LOG10RMSB = 0D0
               END IF
C
               IF ( (RMSAVV.LE.SCFTOL) .AND. (RMSAVB.LE.SCFTOL) .AND. 
     &              (RMSAVA.LE.SCFTOL) .AND. OPDONE ) THEN
C
                  SCFDONE = .TRUE.
                  IF ( (DMFT .OR. LDAU) .AND. ITRSCF.EQ.1 )
     &                 SCFDONE = .FALSE.
C
               ELSE
C
                  SCFDONE = .FALSE.
C
               END IF
C
C-----------------------------------------------------------------------
C              RDLM calculations - iterate the Weiss field
C-----------------------------------------------------------------------
               IF ( FLUCT_DIR_SETTING.EQ.'RDLM      ' ) THEN
                  CALL SCF_BROYDEN_W_WEISS(IPRINT,ITRSCF,SCFMIX,OBS_T,
     &               W_WEISS_ERROR)
C
                  IF ( W_WEISS_ERROR.GT.1D-5 ) SCFDONE = .FALSE.
C
               END IF
C
C ======================================================================
C        SCF cycle done - present type of SCF calculation converged
C ======================================================================
C
               IF ( SCFDONE ) THEN
C
C-----------------------------------------------------------------------
C                             3D systems
C-----------------------------------------------------------------------
C
                  IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C
                     SCFSTATUS = 'CONVERGED'
C
C-----------------------------------------------------------------------
C                             2D systems
C-----------------------------------------------------------------------
C
                  ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
                     IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
C
                        IF ( L_R_RELATION.EQ.'DIFFERENT ' ) THEN
C
                           SCFSTATUS = 'ITR-R-BULK'
C
C-----------------------------------------------------------------------
C                          L-BULK converged                           2D
C-----------------------------------------------------------------------
                        ELSE
C
                           SCFSTATUS = 'ITR-I-ZONE'
C
C-----------------------------------------------------------------------
C                          L-BULK converged                LIR or LIV 2D
C-----------------------------------------------------------------------
                           IF ( (SYSTEM_TYPE(1:3).EQ.'LIR' .OR. 
     &                          SYSTEM_TYPE(1:3).EQ.'LIV') .AND. 
     &                          COPY_POT_LB_TO_IZ ) THEN
C
                              IQ_LB = 0
                              DO IQ_IZ = NQ_L + 1,NQ_L + NQ_I
                                 IQ_LB = IQ_LB + 1
                                 IF ( IQ_LB.GT.NQ_L ) IQ_LB = 1
C
                                 DO IO_IZ = 1,NOQ(IQ_IZ)
                                    IT_IZ = ITOQ(IO_IZ,IQ_IZ)
C
                                    DO IO_LB = 1,NOQ(IQ_LB)
                                       IT_LB = ITOQ(IO_LB,IQ_LB)
C
                                       IF ( Z(IT_IZ).EQ.Z(IT_LB) ) THEN
C
                                         VT(1:NRMAX,IT_IZ)
     &                                      = VT(1:NRMAX,IT_LB)
C
                                         IF ( FULLPOT )
     &                                      VNST(JRNSMIN:NRMAX,
     &                                      1:NLMFPMAX,IT_IZ)
     &                                      = VNST(JRNSMIN:NRMAX,
     &                                      1:NLMFPMAX,IT_LB)
C
C---------------- copy B-field only if B is non-zero for LBULK and IZONE
C                account for sign of B-field, i.e. orientation of moment
                                         BSUM_LB = SUM(BT(1:NRMAX,IT_LB)
     &                                      )
                                         BSUM_IZ = SUM(BT(1:NRMAX,IT_IZ)
     &                                      )
                                         IF ( ABS(BSUM_LB).GT.1D-8 .AND. 
     &                                      ABS(BSUM_IZ).GT.1D-8 ) THEN
C
                                         IF ( BSUM_LB*BSUM_IZ.GT.0D0 )
     &                                      THEN
                                         BT(1:NRMAX,IT_IZ)
     &                                      = BT(1:NRMAX,IT_LB)
                                         ELSE
                                         BT(1:NRMAX,IT_IZ)
     &                                      = -BT(1:NRMAX,IT_LB)
                                         END IF
                                         IF ( FULLPOT ) THEN
                                         IF ( BSUM_LB*BSUM_IZ.GT.0D0 )
     &                                      THEN
                                         BNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                      IT_IZ)
     &                                      = BNST(JRNSMIN:NRMAX,1:
     &                                      NLMFPMAX,IT_LB)
                                         ELSE
                                         BNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                      IT_IZ)
     &                                      = -BNST(JRNSMIN:NRMAX,
     &                                      1:NLMFPMAX,IT_LB)
                                         END IF
                                         END IF
C
                                         END IF
C
                                         QEL(IT_IZ) = QEL(IT_LB)
C
                                         EXIT
                                       END IF
C
                                    END DO
                                 END DO
                              END DO
                           END IF
C
C-----------------------------------------------------------------------
C                          L-BULK converged                       LIR 2D
C-----------------------------------------------------------------------
                           IF ( SYSTEM_TYPE(1:3).EQ.'LIR' ) THEN
C
                              DO IQ_RBULK = NQ_L + NQ_I + 1,NQHOST
C
                                 IQ_LB = IQ_L_IQ_R(IQ_RBULK)
C
                                 CMNTQ(1:NLMQMAD,IQ_RBULK)
     &                              = CMNTQ(1:NLMQMAD,IQ_LB)
                                 VLMMAD_HOST(1:NLMFPMAX,IQ_RBULK)
     &                              = VLMMAD_HOST(1:NLMFPMAX,IQ_LB)
C
                                 DO IO = 1,NOQ(IQ_LB)
                                    IT_LB = ITOQ(IO,IQ_LB)
                                    IT_RBULK = ITOQ(IO,IQ_RBULK)
C
                                    VT(1:NRMAX,IT_RBULK)
     &                                 = VT(1:NRMAX,IT_LB)
                                    BT(1:NRMAX,IT_RBULK)
     &                                 = BT(1:NRMAX,IT_LB)
C
                                    IF ( FULLPOT ) THEN
                                       VNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                    IT_RBULK)
     &                                    = VNST(JRNSMIN:NRMAX,
     &                                    1:NLMFPMAX,IT_LB)
                                       BNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                    IT_RBULK)
     &                                    = BNST(JRNSMIN:NRMAX,
     &                                    1:NLMFPMAX,IT_LB)
                                    END IF
                                    QEL(IT_RBULK) = QEL(IT_LB)
                                 END DO
                              END DO
C
C-----------------------------------------------------------------------
C                          L-BULK converged                       LIV 2D
C-----------------------------------------------------------------------
                           ELSE IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
C
                              DO IQ_RBULK = NQ_L + NQ_I + 1,NQHOST
                                 DO IO = 1,NOQ(IQ_RBULK)
C
                                    IT_RBULK = ITOQ(IO,IQ_RBULK)
C
                                    VT(1:NRMAX,IT_RBULK) = 0.0D0
                                    BT(1:NRMAX,IT_RBULK) = 0.0D0
C
                                    IF ( FULLPOT ) THEN
                                       VNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                    IT_RBULK) = 0.0D0
                                       BNST(JRNSMIN:NRMAX,1:NLMFPMAX,
     &                                    IT_RBULK) = 0.0D0
                                    END IF
                                    QEL(IT_RBULK) = 0.0D0
C
                                 END DO
                              END DO
C
                           END IF
C
                        END IF
C
                     ELSE IF ( SUB_SYSTEM(1:6).EQ.'R-BULK' ) THEN
C
                        SCFSTATUS = 'ITR-I-ZONE'
C
                     ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
C
                        SCFSTATUS = 'CONVERGED'
C
                     ELSE
                        CALL STOP_MESSAGE(ROUTINE,'SUB_SYSTEM ?')
                     END IF
                  END IF
               END IF
C
               CALL POTWRSPR
C
               IF ( NCPAFAIL.NE.0 ) THEN
                  WRITE (6,99017) CPATOL,NCPAFAIL,
     &                            (IECPAFAIL(I),DREAL(ETAB(IECPAFAIL(I),
     &                            1)),CPACHTAB(I),I=1,NCPAFAIL)
                  WRITE (6,'(1X,79(''*''),/)')
               ELSE IF ( NCPA.NE.0 ) THEN
                  WRITE (6,99018)
               END IF
            END IF
         END IF
C==================================================================
C
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C
C=======================================================================
C                  calculate new self energy DMFTSIGMA via DMFT - method
C
         IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
c         IF ( DMFT .AND. .NOT.SCFDONE ) THEN
         IF ( DMFT ) THEN
C
c            IF ( ITRSCF.EQ.1 ) EFERMI = EFERMI0
            if(nprocs>1) call drv_mpi_bcast_r(0,efermi,1)
C
            DO IEPATH = 1,NEPATH
               CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,NETAB(IEPATH),
     &                    ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                    NEMAX)
            END DO
C
            CALL ZCOPY(NETAB(1),ETAB(1,1),1,EGFNEW,1)
            OPDONE = .TRUE.
            IF ( ITRSCF.EQ.1 ) OPDONE = .FALSE.
C
            IF ( NPROCS>1 ) THEN
               CALL DRV_MPI_BCAST_C(0,GFMAT(1,1,1,1),
     &                              NKMMAX*NKMMAX*NTMAX*NETAB(1))
               CALL DRV_MPI_BCAST_R(0,OBS_LT(0,1,1,1),
     &                              4*NOBSMAX*NLMAX*NTMAX)
            END IF
C
            IF ( DMFT_FIX_DYN_SE ) THEN
               WRITE (6,*) 'Dynamical part of sigma fixed.'
               WRITE (6,*) 'Iterating only LDA+U part.'
               CALL DMFT_DRV(ITRSCF,GFMAT,DMFTSIGMA,WEGFOLD,EGFOLD,
     &                       EGFNEW,DMFTMIX,IPRINT,EREFDMFT,OBS_LT,
     &                       SCLNOS)
            ELSE
               CALL DMFT_DRV(ITRSCF,GFMAT,DMFTSIGMA,WEGFOLD,EGFOLD,
     &                       EGFNEW,DMFTMIX,IPRINT,EREFDMFT,OBS_LT,
     &                       SCLNOS)
            END IF
            IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
C
            CALL ZCOPY(NETAB(1),EGFNEW,1,EGFOLD,1)
            CALL ZCOPY(NETAB(1),WETAB(1,1),1,WEGFOLD,1)
C
         END IF
C
C=======================================================================
C                  calculate new self energy DMFTSIGMA via LDAU - method
C
         IF ( MPI_ID.EQ.0 ) THEN
            IF ( LDAU ) THEN
C
               CALL LDAU_DRV(GFMAT,DMFTSIGMA,WETAB(1,1),EREFLDAU,ELDAU,
     &                       OBS_LT,SCLNOS,.TRUE.)
               DO IT = ITBOT,ITTOP
                  IF ( KSELF(IT).EQ.1 ) WRITE (6,99028) IT,NAT(IT)
     &                 *CONC(IT)*ELDAU(IT)
                  ETOT = ETOT + NAT(IT)*CONC(IT)*ELDAU(IT)
               END DO
               WRITE (6,'(/)')
C
            END IF
C
            IF ( .NOT.DMFT_FIX_DYN_SE ) THEN
               IF ( LDAU .OR. DMFT )
     &              CALL DMFT_BROYD_SIG(ITRSCF,IMIX,ISTBRY,ITDEPT,
     &              SCFMIXOP,DMFTSIGMA,SCFDONE,WETAB(1,1),SHFTEF,RMSAVV)
            END IF
         END IF
C=======================================================================
C
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
         IF ( NPROCS>1 ) THEN
            CALL DRV_MPI_BARRIER
C
            CALL DRV_MPI_BCAST_R(0,EMIN,1)
            CALL DRV_MPI_BCAST_R(0,EFERMI,1)
            CALL DRV_MPI_BCAST_R(0,RMSAVV,1)
            CALL DRV_MPI_BCAST_R(0,RMSAVB,1)
            CALL DRV_MPI_BCAST_R(0,RMSAVA,1)
            CALL DRV_MPI_BCAST_R(0,SCFTOL,1)
            CALL DRV_MPI_BCAST_R(0,TOTDOSEF,1)
            CALL DRV_MPI_BCAST_R(0,ETOT,1)
            CALL DRV_MPI_BCAST_R(0,LOG10RMSV,1)
            CALL DRV_MPI_BCAST_R(0,LOG10RMSB,1)
            CALL DRV_MPI_BCAST_R(0,SCLNOS,1)
            CALL DRV_MPI_BCAST_R(0,MUESPN,1)
            CALL DRV_MPI_BCAST_R(0,MUEORB,1)
C
            CALL DRV_MPI_BCAST_I(0,ITRSCF,1)
            CALL DRV_MPI_BCAST_L(0,OPDONE,1)
            CALL DRV_MPI_BCAST_L(0,SCFDONE,1)
C
            IF ( DMFT ) THEN
               CALL DRV_MPI_BCAST_C(0,DMFTSIGMA(1,1,1,1),
     &                              NKMMAX*NKMMAX*NTMAX*NETAB(1))
               CALL DRV_MPI_BCAST_C(0,EREFDMFT(1),NTMAX)
               CALL CINIT(NKMMAX*NKMMAX*NEMAX*NT,GFMAT)
               CALL DRV_MPI_BCAST_L(0,DMFT_FIX_DYN_SE,1)
            END IF
C
            IF ( LDAU ) THEN
               CALL DRV_MPI_BCAST_C(0,DMFTSIGMA(1,1,1,1),
     &                              NKMMAX*NKMMAX*NTMAX*NETAB(1))
               CALL DRV_MPI_BCAST_C(0,EREFLDAU(1),NTMAX)
               CALL CINIT(NKMMAX*NKMMAX*NEMAX*NT,GFMAT)
            END IF
C
            IDIMS = NRMAX*NTMAX
            CALL DRV_MPI_BCAST_R(0,VT(1,1),IDIMS)
            CALL DRV_MPI_BCAST_R(0,BT(1,1),IDIMS)
            IF ( FULLPOT ) THEN
               IDIMS = (NRMAX-JRNSMIN+1)*NLMFPMAX*NTMAX
               CALL DRV_MPI_BCAST_R(0,VNST(1,1,1),IDIMS)
               CALL DRV_MPI_BCAST_R(0,BNST(1,1,1),IDIMS)
            END IF
C
C----------------------------------------------------------------- BREIT
            IF ( BREITINT ) THEN
               IDIMS = NRMAX*NLABIMAX*3*NTMAX
               CALL DRV_MPI_BCAST_R(0,ABIT(1,1,-1,1),IDIMS)
               CALL DRV_MPI_BCAST_R(0,JDNST(1,1,-1,1),IDIMS)
            END IF
C
            CALL DRV_MPI_BARRIER
C
         END IF
c
         if(dmft) call dmft_write_sigepath(dmftsigma,netab(1),egfnew)
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C     Warning for total energies within LDA+U or DMFT
         IF ( DMFT .OR. (LDAU .AND. DMFTDBLC(1:3).EQ.'AMF') ) THEN
            WRITE (6,99023)
            ETOT = 0.0D0
         END IF
C
C
C
C==================================================================
C                             Finalise or prepare new scf iteration
C
         IF ( SCFDONE ) THEN
C
            WRITE (6,99011) ITRSCF,RMSAVV,RMSAVB,EFERMI,SHFTEF,TOTDOSEF,
     &                      MUESPN,MUEORB,ETOT,
     &                      'SCF - cycle converged !!!!!!!!! '
C
            CALL SCF_UCALC_DNDV(1,OBS_LT)
C
            CALL FLUSH(6)
C
            IF ( WRLOG ) WRITE (IFILLOG,99002) ITRSCF0 + ITRSCF,ITRSCF,
     &                          LOG10RMSV,LOG10RMSB,EFERMI,SHFTEF,
     &                          TOTDOSEF,MUESPN,MUEORB,ETOT,SCLNOS
C
            IF ( KMROT.GE.3 ) WRITE (6,99009) IQMVEC,QMVEC,ETOT,MUESPN
            IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
C
         ELSE IF ( ITRSCF.LT.NSCFITER ) THEN
C
            CALL SCF_UCALC_DNDV(0,OBS_LT)
C
            WRITE (6,99011) ITRSCF,RMSAVV,RMSAVB,EFERMI,SHFTEF,TOTDOSEF,
     &                      MUESPN,MUEORB,ETOT,
     &                      'SCF - cycle not converged  -  continue'
C
            IF ( WRLOG ) WRITE (IFILLOG,99002) ITRSCF0 + ITRSCF,ITRSCF,
     &                          LOG10RMSV,LOG10RMSB,EFERMI,SHFTEF,
     &                          TOTDOSEF,MUESPN,MUEORB,ETOT,SCLNOS
C
            IF ( WRLOG ) CALL FLUSH(IFILLOG)
C
            IF ( .NOT.NOCHKVAL ) CALL SCFCHKNVAL(EMIN,ECTOP)
C
            DO IEPATH = 1,NEPATH
               CALL EPATH(IGRID(IEPATH),EMIN,EFERMI,EIMAG,NETAB(IEPATH),
     &                    ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                    NEMAX)
            END DO
            IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
            GOTO 150
C
         ELSE
C
            WRITE (6,99011) ITRSCF,RMSAVV,RMSAVB,EFERMI,SHFTEF,TOTDOSEF,
     &                      MUESPN,MUEORB,ETOT,
     &                 'SCF - cycle not converged - NSCFITER exhausted '
C
            IF ( WRLOG ) WRITE (IFILLOG,99002) ITRSCF0 + ITRSCF,ITRSCF,
     &                          LOG10RMSV,LOG10RMSB,EFERMI,SHFTEF,
     &                          TOTDOSEF,MUESPN,MUEORB,ETOT,SCLNOS
C
         END IF
C==================================================================
C
C=======================================================================
C                         SCF - cycle   END
C=======================================================================
         ITRSCF = 0
C
         IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
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
      END DO
CQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQS
C                       SPIN SPIRAL LOOP END
CQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQSPQS
C
C=======================================================================
C                write log-file if the lattice parameter is scaled
C                return to start the next SCF run
C=======================================================================
      IF ( SCL_ALAT ) THEN
         WRITE (IFILSCLALAT,99003) ALAT,EFERMI,TOTDOSEF,MUESPN,MUEORB,
     &                             ETOT,SCLNOS,LOG10RMSV,LOG10RMSB
C
         CALL FLUSH(IFILSCLALAT)
C
         RETURN
      END IF
C
C=======================================================================
C                 calculate electric field gradient
C=======================================================================
      IF ( FULLPOT .AND. CALC_EFG ) CALL FPSCFNEWPOT(ETOT,ITRSCF,RHO2NS)
C
      IF ( NPROCS>1 ) CALL MPI_FINALIZE(IERR)
C
      DEALLOCATE (FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR)
      IF ( BREITINT ) DEALLOCATE (AMEBI1,AMEBI2,GBIG,GBIL,ABITNEW)
      IF ( ALLOCATED(RHOOPC) ) DEALLOCATE (AOPT0,RHOOPC,RHOOPO)
      DEALLOCATE (RHO2NSX,RHOCHRX,RHOORBX,RHOSPNX,RWORK)
      DEALLOCATE (CPACHTAB,CTOTDOS,CWKE,CWKKMTE,IECPAFAIL,IPROCE)
      DEALLOCATE (MVTET,MVGAM,MVPHI)
      DEALLOCATE (OBS_LT)
      DEALLOCATE (QMGAM,NCORT0,QLZ,BCOR,BCORS)
      DEALLOCATE (CWKT)
      DEALLOCATE (NKPCOR,RINT)
      IF ( DMFT ) THEN
         DEALLOCATE (GFMAT,EGFOLD,EGFNEW)
         DEALLOCATE (DMFTSIGMA)
      END IF
C
      STOP
C=======================================================================
99001 FORMAT ('# SCF ITER    RMSAVV    RMSAVB    EFERMI    SHFTEF',
     &     '  TOTDOSEF MUESPN   MUEORB           ETOT     Lloyd scaling'
     &     )
c99002 FORMAT (2I5,2F10.4,2F10.5,F8.3,2F9.4,F20.8,F12.8)
99002 FORMAT (2I5,2F10.4,2F10.5,F15.3,2F9.4,F20.8,F12.8)
99003 FORMAT (5F10.5,F20.8,3F10.5)
99004 FORMAT (/,A,/,2(A,I4))
99005 FORMAT (/,10X,'*************************************************',
     &        /,10X,'**************   QMVEC   LOOP   *****************',
     &        /,10X,'*************************************************',
     &        /,10X,'IQMVEC1   =',I3,7X,'IQMVEC2   =',I3,/,10X,
     &        'QMVEC0    =',F10.5,/,10X,'DMVEC0    =',F10.5)
99006 FORMAT (//,49('*'),/,
     &        '**************   QMVEC   START  *****************',/,
     &        49('*'),/,10X,'QMTET(IQ) =',9F7.2)
99007 FORMAT (I5,'  QMVEC   ',5X,3F10.5)
99008 FORMAT (10X,'flipping magnetisation for MAGCOUPL = ',A,5X,
     &        'and atom type IT =',I3)
99009 FORMAT (/,49('*'),/,10X,I5,'  QMVEC   ',5X,3F10.5,F15.6,F10.5,/,
     &        49('*'))
99010 FORMAT (/,5X,'execution time for last iteration',F14.3,' sec',/)
99011 FORMAT (/,1X,79('*'),/,I4,' ERR',2E10.3,' EF',2F9.5,' D',F8.3,
     &        ' M',2F9.4,/,5X,'ETOT',F20.8,4X,A,/,1X,79('*'),/)
99012 FORMAT (10X,A,A)
99013 FORMAT (10X,A,I10)
99014 FORMAT (10X,A,5F20.10)
99015 FORMAT (10X,A,L3)
99016 FORMAT (' rms-error for type',i3,':  A = ',1p,d11.4)
99017 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99018 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99019 FORMAT (' ',79('*'),/,' maximum  E_F CORRECTION  exceeded',/,
     &        ' NEW E_F CORRECTION',F9.5,/,' NEW FERMI ENERGY  ',F9.5,/)
99020 FORMAT (/,79('*'),/,10X,' TAU(LAM,LAM'')',/,79('*'),/,5X,' NT =',
     &        I3,' NQ =',I3,' FMT = 2')
99021 FORMAT (/,79('*'),/,10X,' TAU(LAM,LAM'') and M(LAM,LAM'')',/,
     &        79('*'),/,13X,' NQ =',I3,' FMT = 2')
99022 FORMAT (//,1X,79('*'),/,28X,'parameters for SCF-cycle',/,1X,
     &        79('*'))
99023 FORMAT (//,1X,79('*'),/,10X,'WARNING !!!: LDA+U (AMF) or DMFT ',
     &        'Total energies not yet implemented ',/,1X,79('*'))
99024 FORMAT (///,10X,'NETAB(1) = ',I5,/,10X,'NPROCS1 = ',I5,//)
99025 FORMAT (/,' <SCF>:     TOTNOS    ',F15.10,'  TOTNOSX   ',F15.10,/,
     &        ' <SCF>:     TOTDOSEF  ',F15.10,'  TOTDOSEFX ',F15.10)
99026 FORMAT (/,' <<DEBUG>>',/,' <SCF>:  IEPANEL =',I3,' IEPATH =',I3,
     &        ' IE =',I3,' ERYD =',2F15.10,/)
99027 FORMAT (//,1X,79('i'),/,10X,
     &        'INFO:  aspherical potential set to 0',
     &        ' for FULLCHARGE calculations',/,1X,79('i'))
99028 FORMAT (10X,'LDAU energy correction to IT',I3,' = ',F16.12)
99029 FORMAT ('BREAK ',2I3,2X,A,' BP=',I3,' EMIN,EFERMI:',2F16.12)
99030 FORMAT ('BREAK ',2I3,2X,A,' BORN_USE_SPH=',L1)
      END
