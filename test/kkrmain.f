C*==kkrxxx.f    processed by SPAG 6.70Rc at 13:06 on 21 Apr 2017
      PROGRAM KKRXXX
C   ********************************************************************
C   *                                                                  *
C   *  KKR main program                                                *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *  IREL                                                            *
C   *                                                                  *
C   *  0: non-relativistic         para-magnetic                       *
C   *     the non-relativistic (l,ml)-representation for               *
C   *     REAL spherical harmonics is used throughout                  *
C   *     there is NO spin-index                                       *
C   *                                                                  *
C   *  1: scalar relativistic      para-magnetic                       *
C   *     see:   IREL=0                                                *
C   *                                                                  *
C   *  2: scalar relativistic     (spin-polarized)                     *
C   *     the non-relativistic (l,ml)-representation for               *
C   *     REAL spherical harmonics is used ONLY for the K-space        *
C   *     integration; i.e. to deal with multiple scattering           *
C   *     the t-matrix and overlap integrals DZZ, DZJ, SZZ, and SZJ    *
C   *     are set up in the relativistic (kappa,mue)-representation.   *
C   *     This implies that all calculations can be done as for IREL=3 *
C   *                                                                  *
C   *  3: relativistic            (spin-polarized)                     *
C   *     the relativistic (kappa,mue)-representation used throughout  *
C   *     exceptions:                                                  *
C   *     - the k-dependent structure constants G(E,k) are set up in   *
C   *     <STRSET> using the non-relativistic (l,ml)-representation    *
C   *     for REAL spherical harmonics and then transformed            *
C   *     - in the cluster approach the multiple scattering part is    *
C   *     done using the non-relativistic (l,ml)-representation for    *
C   *     COMPLEX spherical harmonics                                  *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *  KMROT                                                           *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site IQ   *
C   *  2: global COMMON rotation of the magnetisation                  *
C   *  3: spin spiral    Theta =  90 }    QMVEC <> 0-vector            *
C   *  4: spin spiral    Theta <> 90 }    IREL = 2                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TABLES,ONLY:TAB_CHSYM
      USE MOD_SCF,ONLY:SCFSTATUS,SCFSTATUS_HOST
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_CALCMODE,ONLY:ORBPOL,BREITINT,EXTFIELD,PROGNAME,IREL,
     &    SELFENERGY,ISMQHFI,IHYPER,ICORE,IXRAY,ITEST,MOL,NONMAG,KMROT,
     &    TXTKMROT,LLOYD,KKRMODE,IREL_HOST,TASK,GF_CONV_RH,FULLCHARGE,
     &    LHS_SOL_EQ_RHS_SOL,MOMENTS_ROTATED,SETPOT,PUBLIC_VERSION
      USE MOD_SYMMETRY,ONLY:IQREPQ,ISYMGENQ,IQEQ,NQEQ,SNTAUUVMAX,
     &    NELMTMAX,NO_SYMMETRY,NO_SYMMETRY_CLU,NO_SYMMETRY_LINRESP,
     &    SYMMETRIZE_MSS,SYMMETRIZE_TAU,SYMMETRIZE_RHO,SYMMETRIZE_POT
      USE MOD_SITES,ONLY:QMGAM,QMTET,QMPHI,MDIRQ,NOMAX,NOQ,ICPA,IMQ,NQ,
     &    NQMAX,IQAT,ITOQ,QMVEC,NQ_L,NQ_R,QBAS,ALFDEG,BETDEG,GAMDEG,
     &    NQHOST,NQCLU,IQBOT,IQTOP,MAGROT_Q
      USE MOD_TYPES,ONLY:CONC,NLT,NLMFPT,NCORT,NVALT,TXT_T,NAT,Z,RHOSPN,
     &    RHOSPNC,VT,BT,VAMEF,VAMEG,NT,NTMAX,NPOTMAX,NCPLWFMAX,NLAFPMAX,
     &    NVALTOT,NLSHELLMAX,BEXT,ITBOT,ITTOP,NLIN_T,NKM_T,NTCLU,NTHOST,
     &    IMT,NLFP,NLMFP,NLFPMAX,NLMFPMAX
      USE MOD_ANGMOM,ONLY:NLQ,NKMQ,NLINQ,IND0Q,NL,NLMAX,NK,NKMAX,NLM,
     &    NKM,NKMMAX,NKMPMAX,NLIN,NMUEMAX,LINMAX,NKKR,NLMMAX,NLABIMAX,
     &    WKM1,WKM2,WKM3
      USE MOD_RMESH,ONLY:RWS,RMT,LMISF,NM,NMMAX,NRMAX,SPHERCELL,FULLPOT,
     &    FINITE_NUCLEUS,NMCLU,NMHOST,NLSFMAX,NLMSFMAX,NLSF,NLMSF
      USE MOD_CPA,ONLY:CPAFIL,WRCPA,NCPA,USENLCPA
      USE MOD_ENERGY,ONLY:EMIN,EMAX,SPLITSS,NETAB
      USE MOD_KSPACE,ONLY:IBZINT,NKTABMAX
      USE MOD_FILES,ONLY:LSYSTEM,LDATSET,LDATSET0,LINFO,RECLNGWF,
     &    RECLNGWF_SPH,DATSET,DATSET0,SYSTEM,INFO,IOTMP,LRECREAL8,
     &    IFILLDAU,IFILCORWF,TAUFIL,IFILGFWF,IFILCBWF,IFILCBWF_SPH,
     &    IFILCBWF_LHS,IFILMEZZL,TITLE,LTITLE,POTFMTOUT,IFILINP,IPRINT,
     &    WRLOG,WRKAPDOS,WRMAT,WRPOLAR,WRTAU,WRTAUMQ,NOWRDOS,RDTAU,
     &    RDTAUMQ,PLOT2DPR,PLOTPRS,POTFIL,DOSFIL,NOSFIL,LPOTFIL,IFILTAU,
     &    RDDOS,IFILDOS,IFILNOS,WRNOS,CPU_TIME_PROGRAM_START,
     &    CPU_TIME_LAST_CALL,WALL_TIME_PROGRAM_START,
     &    WALL_TIME_LAST_CALL,IFILBUILDBOT,WRBUILDBOT,FOUND_SECTION,
     &    FOUND_INTEGER,N_FOUND
      USE MOD_LATTICE,ONLY:SWS,ALAT,BRAVAIS,SYSTEM_DIMENSION,
     &    SYSTEM_TYPE,SUB_SYSTEM,ABAS,ABAS_L,ABAS_R
      USE MOD_MPI,ONLY:MPI_ID,MPI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ITESTMAX
      PARAMETER (ITESTMAX=8)
      LOGICAL DUMP_MODULES
      PARAMETER (DUMP_MODULES=.FALSE.)
C
C Local variables
C
      CHARACTER*80 ADSI,FILNAM
      LOGICAL BREITINT5,FOUND,FULLPOT5,INITELOOP,POSANIPREP,RUNELOOP,
     &        SPHERCELL5,USELPLS1,USETAUNN
      CHARACTER*3 DOSCORSYS,DOSREP
      COMPLEX*16 ERYD
      REAL*8 ETOP,RASRAD,RASSCL,RMTRED0(:),RSUM,VUC
      INTEGER I,IA,IA_ERR,IBLK,IBLKTOP,IDIMS,IEQ,IFLAG,IM_QAUX(:),IO,
     &        IPOT,IQ,IQ_ATAUX(:,:),IT,IT0_TAUX(:),IT_OQAUX(:,:),IZ,J,
     &        J1,J1TOP,JQ,KBZI,KMOL,KMROT_COMMON,LADSI,LF,LPROGNAME,MC,
     &        NA_TAUX(:),NETAU,NFPLIM,NLES,NLQ0(:),NSFLIM,NT0,NTAUX,
     &        NTLIM,POTFMTINP,Z_TAUX(:)
      CHARACTER*40 ROUTINE,TXTTEST(0:ITESTMAX)
      CHARACTER*10 STR10
      INTEGER TABNVAL
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
      DATA TXTTEST/'         NO TEST OPTIONS SET           ',
     &     'SET TAU=T-SS IN <ELOOP>                ',
     &     '<SYMCHECK>  check symmetry of  TAU     ',
     &     '<CHITKTKCHK>  TAU(k)*TAU(k) integration',
     &     'CHECK <SETTMAT> + <DIRAC>     CHECK2=T ',
     &     'test of Lloyd-formula                  ',
     &     'test relativistic routines             ',
     &     'test analytical property of G-function ',
     &     'test structure constants               '/
C ======================================================================
      DATA KBZI/1/,RUNELOOP/.TRUE./,INITELOOP/.TRUE./
      DATA LADSI/0/
      DATA FULLPOT5/.FALSE./
      DATA BREITINT5/.FALSE./
      DATA NETAU/0/
      DATA DOSREP,DOSCORSYS/'XXX','LOC'/
      DATA KMROT_COMMON/0/
      DATA USETAUNN/.FALSE./
      DATA POSANIPREP/.FALSE./
      DATA ERYD/(999999D0,999999D0)/
C
C ======================================================================
C
C-------------------------------------- variables depending on NQ and NT
      ALLOCATABLE IQ_ATAUX,IT_OQAUX,NLQ0,IM_QAUX
C---------------------------------------- variables depending only on NT
      ALLOCATABLE NA_TAUX,Z_TAUX,IT0_TAUX
C---------------------------------------- variables depending only on NM
      ALLOCATABLE RMTRED0
C
      CALL CPU_TIME(CPU_TIME_PROGRAM_START)
      CALL CPU_TIME(CPU_TIME_LAST_CALL)
      CALL SYSTEM_CLOCK(WALL_TIME_PROGRAM_START)
      WALL_TIME_LAST_CALL = WALL_TIME_PROGRAM_START
C
C=======================================================================
C                initialize some parameters for mod_files
C=======================================================================
C
      CALL INIT_MOD_FILES
C
C ======================================================================
C                   initialize MPI - parameters
C ======================================================================
C
      CALL INIT_MOD_MPI
C
C ======================================================================
C
      PROGNAME = 'KKRXXX    '
C
      IF ( PROGNAME(4:6).EQ.'SCF' ) TASK = 'SCF       '
C
      LPROGNAME = LEN_TRIM(PROGNAME)
C
      WRITE (6,99013) PROGNAME(1:LPROGNAME)
      IF ( PUBLIC_VERSION ) WRITE (6,99014)
C
      ROUTINE = PROGNAME(1:LPROGNAME)
C
      CALL WRDATE('programm execution',6)
      CALL GITSTAMP
C
      REWIND (IFILINP,IOSTAT=IFLAG)
      IF ( IFLAG.NE.0 ) THEN
         STR10 = PROGNAME
         CALL STRING_CONVERT_TO_LC(STR10)
         WRITE (6,99003) PROGNAME,STR10
         STOP
      END IF
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C                        file-info 1
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
      CALL INPUT_FIND_SECTION('CONTROL',1)
C
      CALL SECTION_SET_STRING('DATASET',DATSET,'9999',1)
      CALL SECTION_SET_STRING('ADSI',ADSI,'9999',0)
      CALL SECTION_SET_STRING('POTFIL',POTFIL,'9999',1)
      CALL SECTION_FIND_KEYWORD('NOWRDOS',NOWRDOS)
      CALL SECTION_SET_INTEGER('PRINT',IPRINT,9999,0)
C
      LPOTFIL = LEN_TRIM(POTFIL)
      LDATSET = LEN_TRIM(DATSET)
      DATSET0 = DATSET
      LDATSET0 = LDATSET
C
      LADSI = LEN_TRIM(ADSI)
      IF ( LADSI.GT.0 ) THEN
         DATSET = DATSET(1:LDATSET)//'_'//ADSI(1:LADSI)
         LDATSET = LDATSET + 1 + LADSI
      END IF
C
      OPEN (UNIT=4,STATUS='OLD',ERR=100,FILE=POTFIL)
      DOSFIL = DATSET(1:LDATSET)//'.dos'
      TAUFIL = DATSET(1:LDATSET)//'.tau'
C
C=======================================================================
C                   determine the record length used for a real variable
C=======================================================================
C
      CALL GETLRECREAL(IOTMP,IPRINT,LRECREAL8)
C
C=======================================================================
C     read potfile to fix NQ and NM and to get an upper limit for NT
C              allocate all arrays depending ONLY on NQ
C=======================================================================
C
      CALL POTRD1(IPRINT,POTFMTINP,SCFSTATUS,NQ,NTLIM,NM)
C
      SCFSTATUS_HOST = SCFSTATUS
C
      NQHOST = NQ
      NMHOST = NM
      IQBOT = 1
      IQTOP = NQ
C
C=======================================================================
C  check for additional sites and meshes in case of cluster calculation
C=======================================================================
C
      CALL CLUGETLIMITS(NQCLU,NMCLU)
C
      NQMAX = NQ + NQCLU
      NMMAX = NM + NMCLU
C
      IF ( NQCLU.GT.0 .AND. PROGNAME(1:3).EQ.'KKR' )
     &     CALL STOP_MESSAGE(ROUTINE,
     &     'the input requires a  CLU*  program but  KKR*  was called')
C
C
C=======================================================================
C                allocate site dependent arrays
C=======================================================================
C
      CALL INIT_MOD_SITES
C
C=======================================================================
C        read all variables from potfile that might influence NT
C        set temporary  NT  and allocate all arrays depending on NT
C=======================================================================
C
      NTMAX = NTLIM
C
      ALLOCATE (NA_TAUX(NTLIM),Z_TAUX(NTLIM),IT0_TAUX(NTLIM))
      ALLOCATE (IM_QAUX(NQMAX))
      ALLOCATE (RMTRED0(NMMAX))
      ALLOCATE (IQ_ATAUX(NQMAX,NTLIM),IT_OQAUX(NTLIM,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQ_ATAUX')
      IQ_ATAUX = 9999999
      IT_OQAUX = 9999999
C
      CALL POTRD2(POTFMTINP,IREL,BRAVAIS,QBAS,SYSTEM,LSYSTEM,NQ,NT0,NM,
     &            Z_TAUX,NA_TAUX,IM_QAUX,IQ_ATAUX,NOQ,IT_OQAUX,RMTRED0,
     &            SCFSTATUS,NONMAG,BREITINT,ORBPOL,EXTFIELD,KMROT,QMTET,
     &            QMPHI,QMVEC,FULLPOT,SPHERCELL,NQMAX,NTLIM,NMMAX,
     &            SYSTEM_DIMENSION,SYSTEM_TYPE,ALAT,ABAS,ABAS_L,ABAS_R,
     &            NQ_L,NQ_R)
C
C-----------------------------------------------------------------------
C optimize basis geometry - ONLY for SCF start of 2D or 3D (host)system
C-----------------------------------------------------------------------
C
      CALL OPTIMIZE_BASIS_DRIVE(ABAS,QBAS,NQ,NQMAX)
C
      IREL_HOST = IREL
C
      IF ( (SCFSTATUS(1:5).EQ.'START' .AND. PROGNAME(4:6).NE.'SCF') .OR. 
     &     (SCFSTATUS_HOST(1:5).EQ.'START' .AND. 
     &     (PROGNAME(1:3).EQ.'CLU' .OR. PROGNAME.EQ.'KKRGEN    ' .OR. 
     &     PROGNAME(4:6).EQ.'CHI')) ) THEN
C
         WRITE (6,99020) PROGNAME,SCFSTATUS,SCFSTATUS_HOST
         CALL STOP_MESSAGE(ROUTINE,
     &       'SCFSTART for   PROG <> KKRSCF  >>>>  call  KKRSCF instead'
     &       )
      END IF
C
      NT = NT0
C
      CALL INIT_MOD_LATTICE(IPRINT,SYSTEM_DIMENSION)
C
C=======================================================================
C         rotate lattice vectors if   ROTLAT   set in input file
C=======================================================================
C
      CALL SYMROTLAT
C
C=======================================================================
C             read all missing information from input file
C          eventually the potfile settings will be overwritten
C=======================================================================
C
      CALL INITVAR(NQ,TASK,POTFMTINP,POTFMTOUT,ITEST,TXTTEST,IREL,ICORE,
     &             IHYPER,ISMQHFI,IXRAY,RDTAU,RDTAUMQ,RDDOS,WRLOG,WRMAT,
     &             WRTAU,WRTAUMQ,USETAUNN,USELPLS1,RUNELOOP,INITELOOP,
     &             SPLITSS,NOWRDOS,IBZINT,KMROT,KMROT_COMMON,QMVEC,
     &             MDIRQ,ALFDEG,BETDEG,GAMDEG,LLOYD,NONMAG,PLOTPRS,
     &             PLOT2DPR,BREITINT5,NQMAX,ITESTMAX)
C
C=======================================================================
C             prepare spin spiral calculations
C=======================================================================
C
      IF ( KMROT.GE.3 ) CALL INIT_SPIRAL(KMROT,QMVEC,MDIRQ,QMTET,QMPHI,
     &     NQ,NQMAX)
C
C=======================================================================
C                 allocate arrays for TB calculations
C=======================================================================
C
C      IF( KKRMODE(1:6).EQ.'TB-KKR' )
      CALL INIT_MOD_TB(NM,NQMAX)
C
C=======================================================================
C         parameter for angular momentum expansion  NL = l_max + 1
C=======================================================================
C
      ALLOCATE (NLQ0(NQMAX))
C--------------------------------------- NL according to site occupation
      DO IQ = 1,NQ
C
         NLQ0(IQ) = 3
C
         DO IO = 1,NOQ(IQ)
            IZ = Z_TAUX(IT_OQAUX(IO,IQ))
            IF ( IZ.GT.56 .AND. IZ.LT.71 ) NLQ0(IQ) = MAX(4,NLQ0(IQ))
            IF ( IZ.GT.88 ) NLQ0(IQ) = MAX(4,NLQ0(IQ))
         END DO
C
C--------- in case of VBPES going beyond single scattering approximation
         IF ( USETAUNN .OR. USELPLS1 ) NLQ0(IQ) = NLQ0(IQ) + 1
C
      END DO
C
C-------------------------------------------------- NL set by input file
      CALL INPUT_FIND_SECTION('SITES',0)
      IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER_ARRAY('NL',NLQ0,
     &     N_FOUND,NQMAX,0,9999,0)
C
C------------------------------------------------------ NL for vacancies
      CALL SECTION_SET_INTEGER('NLES',NLES,9999,0)
C
      IF ( FOUND_INTEGER ) THEN
         DO IQ = 1,NQ
            IF ( NOQ(IQ).EQ.1 ) THEN
               IF ( Z_TAUX(IT_OQAUX(1,IQ)).EQ.0 ) NLQ0(IQ) = NLES
            END IF
         END DO
      END IF
C
      NL = 0
      DO IQ = 1,NQ
         NL = MAX(NL,NLQ0(IQ))
      END DO
      NL = MAX(2,NL)
C
C=======================================================================
C                      NL-derived parameters
C=======================================================================
C
      NLMAX = NL
      NLMMAX = NLMAX**2
      NKMMAX = 2*NLMMAX
      NKMAX = 2*NLMAX - 1
      NKMPMAX = NKMMAX + 2*NLMAX
      NMUEMAX = 2*NLMAX
      LINMAX = 2*NLMAX*(2*NLMAX-1)
      NLABIMAX = NLMAX
C
      NLFP = 2*NL - 1
      NLMFP = NLFP**2
      NLFPMAX = NLFP
      NLMFPMAX = NLMFP
C
      NLSF = 4*(NL-1) + 1
      NLMSF = NLSF**2
      NLSFMAX = NLSF
      NLMSFMAX = NLMSF
C
C=======================================================================
C             initialize general tables dependent on  NL
C=======================================================================
C
      CALL INIT_MOD_ANGMOM_GEN(NL)
C
C=======================================================================
C          determine symmetry of the system and fix NT this way
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TAU',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_FIND_KEYWORD('MOL',MOL)
C
      CALL INIT_MOD_SYMMETRY(TASK,INITELOOP,IREL,NONMAG,MOL,KMROT,QMVEC,
     &                       IPRINT,NT,NT0,NA_TAUX,IQ_ATAUX,IT_OQAUX,
     &                       IT0_TAUX,NL,BRAVAIS,NKMMAX,NTMAX,NTLIM)
C
      NTAUX = NT
C
      NTHOST = NT
C     NTMAX = NT
      ITBOT = 1
      ITTOP = NT
C
C=======================================================================
C               allocate and supply  IQAT  and  ITOQ
C=======================================================================
C
C=======================================================================
C         check for additional types in case of cluster calculation
C=======================================================================
C
      IF ( NQCLU.NE.0 ) CALL CLUSYMMETRY(NTCLU,NMCLU)
C
      IF ( NTCLU.NE.0 ) NTMAX = NTMAX + NTCLU
C
      ALLOCATE (IQAT(NQMAX,NTMAX),ITOQ(NTMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQAT')
      IQAT = 9999999
      ITOQ = 9999999
C
      DO IT = 1,NT
         DO IQ = 1,NQ
            ITOQ(IT,IQ) = IT_OQAUX(IT,IQ)
            IQAT(IQ,IT) = IQ_ATAUX(IQ,IT)
         END DO
      END DO
C
      IMQ(1:NQMAX) = IM_QAUX(1:NQMAX)
C
C=======================================================================
C              dimensions for full potential calculations
C=======================================================================
C FULLPOT = .F. AND FULLPOT5 = .T. --> switch from ASA to FULL POTENTIAL
C
      SPHERCELL5 = .FALSE.
      CALL INPUT_FIND_SECTION('SCF',0)
      IF ( .NOT.FULLPOT ) THEN
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('FULLPOT',FULLPOT5)
            IF ( FULLPOT5 ) CALL SECTION_FIND_KEYWORD('SPHERCELL',
     &           SPHERCELL5)
            IF ( SPHERCELL5 ) SPHERCELL = .TRUE.
         END IF
      END IF
C
      CALL SECTION_FIND_KEYWORD('FULLCHARGE',FULLCHARGE)
C
      IF ( FULLCHARGE ) FULLPOT5 = .TRUE.
C
      IF ( FULLPOT .OR. FULLPOT5 ) THEN
C
         CALL FPDIMENSIONS(IQAT,NA_TAUX,NTMAX,NFPLIM,NSFLIM)
C
         NPOTMAX = 2*NTMAX
         NLAFPMAX = 2*NLMAX
C
      END IF
C
C=======================================================================
C         allocate all arrays depending on radial meshes
C=======================================================================
C
      CALL INIT_MOD_RMESH
C
C=======================================================================
C            construct LEBEDEV mesh for angular integrations
C=======================================================================
C
      IF ( FULLPOT .OR. FULLPOT5 ) CALL FPSPHERE
C
C=======================================================================
C               initialize radial mesh parameters for
C             FULLPOT calculations and shape functions
C=======================================================================
C
      CALL INIT_MOD_RMESH_FULLPOT(NTLIM,NT0,NA_TAUX,RMTRED0,NSFLIM,
     &                            FULLPOT,FULLPOT5)
C
C=======================================================================
C                      calculation mode  STEP 1
C=======================================================================
C
      CALL INIT_MOD_CALCMODE1(BEXT)
C
C=======================================================================
C              allocate all arrays depending on NT
C=======================================================================
C
      IF ( BREITINT5 ) BREITINT = .TRUE.
C
      CALL INIT_MOD_TYPES(NLSHELLMAX,FULLPOT,FULLPOT5,SPHERCELL)
C
      LMISF = 0
C
      CALL POTRDSPR(POTFMTINP)
C
C-----------------------------------------------------------------------
C optimize basis geometry - ONLY for SCF start of 2D or 3D (host)system
C once more as <POTRDSPR> rereads orginal settings of  QBAS
C-----------------------------------------------------------------------
C
      CALL OPTIMIZE_BASIS_DRIVE(ABAS,QBAS,NQ,NQMAX)
C
      IF ( BREITINT5 ) BREITINT = .TRUE.
C
      NT = NTAUX
C
      DO IT = 1,NT
         Z(IT) = Z_TAUX(IT)
         NAT(IT) = NA_TAUX(IT)
         DO IQ = 1,NQ
            ITOQ(IT,IQ) = IT_OQAUX(IT,IQ)
            IQAT(IQ,IT) = IQ_ATAUX(IQ,IT)
         END DO
      END DO
C
      IF ( FULLPOT .OR. FULLPOT5 .OR. SPHERCELL ) THEN
         DO IT = 1,NT
            NLMFPT(IT) = (2*(NLQ0(IQAT(1,IT))-1)+1)**2
         END DO
      END IF
C
      IF ( .NOT.SPHERCELL .AND. SPHERCELL5 ) SPHERCELL = .TRUE.
C
      DEALLOCATE (Z_TAUX,NA_TAUX,IT_OQAUX,IQ_ATAUX)
C
C=======================================================================
C ONCE MORE:  rotate lattice vectors if   ROTLAT   set in input file
C=======================================================================
C
      CALL SYMROTLAT
C
C=======================================================================
C                 manipulate potential if requested
C=======================================================================
C
      SELECT CASE (SETPOT)
C
      CASE ('UP')
C
         DO IT = 1,NT
            VT(:,IT) = VT(:,IT) + BT(:,IT)
         END DO
         BT(:,:) = 0D0
C
         WRITE (6,'(/,19X,''paramagnetic calculation for V=V(up)'')')
C
      CASE ('DOWN')
C
         DO IT = 1,NT
            VT(:,IT) = VT(:,IT) - BT(:,IT)
         END DO
         BT(:,:) = 0D0
C
         WRITE (6,'(/,19X,''paramagnetic calculation for V=V(down)'')')
C
      CASE ('AVG')
C
         BT(:,:) = 0D0
C
         WRITE (6,'(/,19X,''paramagnetic calculation for V=<V>'')')
C
      CASE DEFAULT
C
C
      END SELECT
C
      IF ( NONMAG ) THEN
         CALL RINIT(NRMAX*NTMAX,BT)
         CALL RINIT(NRMAX*NTMAX,RHOSPN)
         CALL RINIT(NRMAX*NTMAX,RHOSPNC)
         WRITE (6,'(/,19X,''paramagnetic calculation  B set to 0'')')
      ELSE
         WRITE (6,'(/,19X,''spin-polarized calculation'')')
      END IF
C
      WRITE (6,'(/,19X,''SELFENERGY applied:'',A)') SELFENERGY
C
C=======================================================================
C                      finite nucleus calculation
C=======================================================================
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_FIND_KEYWORD('FINNUC',
     &     FINITE_NUCLEUS)
C
      IF ( FINITE_NUCLEUS ) THEN
         CALL POTFINNUC(VT)
         WRITE (6,'(/,19X,''a finite nucleus will be assumed '')')
      END IF
C
C=======================================================================
C                        type IT dependent settings
C=======================================================================
C
C-----------------------------------------------------------------------
C       update all type - dependent variables if necessary   (NT>NT0)
C-----------------------------------------------------------------------
C
C
      CALL INPUT_FIND_SECTION('SITES',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER_ARRAY('NVAL',NVALT,
     &     NT,NTMAX,1,9999,0)
C
      DO IT = 1,NT0
         IF ( ABS(NVALT(IT)).LT.1D-5 ) NVALT(IT) = TABNVAL(Z(IT))
      END DO
C
      IF ( NT0.NE.NT ) THEN
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &        'reduced symmetry + FULLPOT:  still to be done')
C
         CALL TYPES_SWAP(NT0,IT0_TAUX)
C
      END IF
C
      NVALTOT = 0.0D0
      DO IT = 1,NT
         TXT_T(IT) = TAB_CHSYM(Z(IT))
         NCORT(IT) = Z(IT) - NVALT(IT)
         NVALTOT = NVALTOT + CONC(IT)*NAT(IT)*NVALT(IT)
      END DO
C
C=======================================================================
C                  site IQ dependent settings
C=======================================================================
C
      NCPA = 0
      NOMAX = 0
      DO IQ = 1,NQ
C
         NOQ(IQ) = 0
         DO IT = 1,NT
            DO IA = 1,NAT(IT)
               IF ( IQAT(IA,IT).EQ.IQ ) THEN
                  NOQ(IQ) = NOQ(IQ) + 1
                  ITOQ(NOQ(IQ),IQ) = IT
               END IF
            END DO
         END DO
         NOMAX = MAX(NOMAX,NOQ(IQ))
C
C-------------------------- NO CPA if number of occupant NOQ = 1 on site
         IF ( NOQ(IQ).EQ.1 ) THEN
            ICPA(IQ) = 0
         ELSE
            ICPA(IQ) = 1
         END IF
C
         RSUM = 0.0D0
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            RSUM = RSUM + CONC(IT)
         END DO
C
         IF ( DABS(1.0D0-RSUM).GT.1D-6 )
     &         CALL STOP_MESSAGE(ROUTINE,'SUM OF CONC <> 1')
C
         NCPA = NCPA + ICPA(IQ)
C
         NQEQ(IQ) = 0
         DO JQ = 1,NQ
            IF ( NOQ(IQ).EQ.NOQ(JQ) ) THEN
               IFLAG = 1
               DO IO = 1,NOQ(IQ)
                  IF ( ITOQ(IO,IQ).NE.ITOQ(IO,JQ) ) IFLAG = 0
               END DO
               IF ( IFLAG.EQ.1 ) THEN
                  NQEQ(IQ) = NQEQ(IQ) + 1
                  IQEQ(NQEQ(IQ),IQ) = JQ
               END IF
            END IF
         END DO
      END DO
C
C---------------------------------------------- check consistency of NLQ
C
      DO IQ = 1,NQ
         DO I = 1,NQEQ(IQ)
            JQ = IQEQ(I,IQ)
            IF ( NLQ0(IQ).NE.NLQ0(JQ) ) WRITE (6,99002) IQ,NLQ0(IQ),JQ,
     &           NLQ0(JQ)
         END DO
      END DO
C
C=======================================================================
C
      CALL EXTEND_TXT_T
C
      WRITE (6,99010) NQ
C
      DO IQ = 1,NQ
C
         IF ( ABS(QMTET(IQ)).LT.1D-6 .AND. ABS(QMPHI(IQ)).LT.1D-6 .AND. 
     &        ABS(QMGAM(IQ)).LT.1D-6 ) THEN
            MAGROT_Q(IQ) = .FALSE.
         ELSE
            MAGROT_Q(IQ) = .TRUE.
         END IF
C
         WRITE (6,99001) IQ,NLQ0(IQ),IMQ(IQ),NINT(QMGAM(IQ)),
     &                   NINT(QMTET(IQ)),NINT(QMPHI(IQ)),RWS(IMQ(IQ)),
     &                   ICPA(IQ),NOQ(IQ),IQREPQ(IQ),ISYMGENQ(IQ),
     &                   (IQEQ(IEQ,IQ),IEQ=1,NQEQ(IQ))
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            WRITE (6,'(43X,I5,1X,A,2X,I4,F5.2)') IT,TAB_CHSYM(Z(IT)),
     &             NVALT(IT),CONC(IT)
         END DO
      END DO
C
      WRITE (6,99011) SWS,NVALTOT
      IF ( KMROT.EQ.1 .OR. KMROT.EQ.2 ) MOMENTS_ROTATED = .TRUE.
      WRITE (6,99015) KMROT,TXTKMROT(KMROT),MOMENTS_ROTATED,QMVEC
C
      WRITE (6,99016) NT
      DO IT = 1,NT
         WRITE (6,99017) IT,TXT_T(IT),NLQ0(IQAT(1,IT)),IMT(IT),
     &                   RMT(IMT(IT)),RWS(IMT(IT)),NAT(IT),CONC(IT),
     &                   (IQAT(IA,IT),IA=1,NAT(IT))
      END DO
      WRITE (6,'(//)')
C
C=======================================================================
C
      IF ( .NOT.(FULLPOT .OR. FULLPOT5) ) THEN
         VUC = ABAS(1,1)*(ABAS(2,2)*ABAS(3,3)-ABAS(2,3)*ABAS(3,2))
     &         + ABAS(1,2)*(ABAS(2,3)*ABAS(3,1)-ABAS(2,1)*ABAS(3,3))
     &         + ABAS(1,3)*(ABAS(2,1)*ABAS(3,2)-ABAS(2,2)*ABAS(3,1))
         VUC = ABS(VUC)*ALAT**3
C
         IF ( ABS(VUC-NQ*SWS**3*4D0*PI/3D0).GT.1D-5 ) WRITE (6,99012)
     &        PROGNAME,VUC,NQ*SWS**3*4D0*PI/3D0,
     &        ABS(VUC-NQ*SWS**3*4D0*PI/3D0)
      END IF
C
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TAU',0)
C
C --------------- allow to use the old section name KMESH instead of TAU
      IF ( .NOT.FOUND_SECTION ) CALL INPUT_FIND_SECTION('KMESH',0)
C
      IF ( FOUND_SECTION ) THEN
C ------------ allow to use the old format BZINT=CLUSTER
         CALL SECTION_SET_STRING('BZINT',STR10,'9999',0)
         IF ( STR10(1:7).EQ.'CLUSTER' ) IBZINT = 0
C ------------
         CALL SECTION_FIND_KEYWORD('CLUSTER',FOUND)
         IF ( FOUND ) IBZINT = 0
         CALL SECTION_FIND_KEYWORD('MOL',MOL)
         IF ( MOL ) IBZINT = 0
      END IF
C
      IF ( IBZINT.EQ.0 ) KBZI = 0
C
      KMOL = 1 - KBZI
C
C=======================================================================
C
      IF ( LSYSTEM.EQ.0 ) CALL GETSYSTEM(80)
C
      IF ( IREL.EQ.0 ) THEN
         INFO = 'KKR'
         LINFO = 3
      ELSE IF ( IREL.LE.2 ) THEN
         INFO = 'SRA-KKR'
         LINFO = 7
      ELSE
         INFO = 'SPR-KKR'
         LINFO = 7
      END IF
      IF ( .NOT.(FULLPOT .OR. FULLPOT5) ) THEN
         INFO = INFO(1:LINFO)//'-ASA'
         LINFO = LINFO + 4
      ELSE IF ( SPHERCELL ) THEN
         INFO = 'SPH-FP-'//INFO(1:LINFO)
         LINFO = 7 + LINFO
      ELSE
         INFO = 'FP-'//INFO(1:LINFO)
         LINFO = 3 + LINFO
      END IF
      TITLE = INFO(1:LINFO)//' calculation for '//
     &        SYSTEM(1:MAX((80-LINFO-17),LSYSTEM))
      LTITLE = LEN_TRIM(TITLE)
C
C=======================================================================
C                      calculation mode  STEP 2
C=======================================================================
C
      CALL INIT_MOD_CALCMODE2(ITEST,KMROT_COMMON,ALFDEG,BETDEG,GAMDEG,
     &                        BEXT,FULLPOT5)
C
C=======================================================================
C            initialize arrays and tables dependent on  NL
C=======================================================================
C
      CALL INIT_MOD_ANGMOM(IREL,ISMQHFI,IPRINT,NLQ0,NKKR,NLFPMAX)
C
C=======================================================================
C                open file TAUFIL for the TAU - matrices
C=======================================================================
      IF ( TASK.EQ.'PHONONS   ' .OR. TASK.EQ.'CHI       ' .OR. 
     &     (TASK(1:5).EQ.'SIGMA' .AND. .NOT.USENLCPA) .OR. TASK(1:6)
     &     .EQ.'MAGNET' ) THEN
         RDTAU = .FALSE.
         RDTAUMQ = .FALSE.
         WRTAU = .FALSE.
         WRTAUMQ = .FALSE.
      END IF
C
      IF ( TASK(1:3).EQ.'BSF' .AND. NCPA.EQ.0 ) THEN
         RDTAUMQ = .FALSE.
         INITELOOP = .TRUE.
      END IF
      IF ( TASK(1:5).EQ.'ARPES' .AND. PROGNAME(1:7).NE.'KKRSPEC' ) THEN
         RUNELOOP = .FALSE.
         IF ( NCPA.NE.0 ) RDTAUMQ = .TRUE.
      END IF
C
      IF ( RDTAU .OR. RDTAUMQ ) THEN
C
         OPEN (UNIT=IFILTAU,FILE=TAUFIL)
C
         WRITE (6,99004)
C
         CALL READTAU(IFILTAU,ERYD,0,NETAU,WKM1,0,NTHOST+NTCLU,RDTAU,
     &                WKM2,WKM3,0,NQHOST+NQCLU,0,NKMMAX,NKMMAX,IPRINT)
C
         IF ( NETAU.EQ.0 ) THEN
C
            IF ( RDTAU ) THEN
               WRTAU = .TRUE.
               RDTAU = .FALSE.
            ELSE
               WRTAUMQ = .TRUE.
               RDTAUMQ = .FALSE.
            END IF
            RUNELOOP = .TRUE.
            INITELOOP = .TRUE.
C
            IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &           TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' )
     &           NOWRDOS = .FALSE.
            WRITE (6,99005)
            REWIND 9
C
         ELSE
C
            RUNELOOP = .FALSE.
            WRITE (6,'(/,1X,79(''*''))')
         END IF
      END IF
C
      IF ( (WRTAU .OR. WRTAUMQ) .AND. (MPI_ID.EQ.0) ) THEN
C
         OPEN (UNIT=IFILTAU,FILE=TAUFIL)
C
         IF ( WRTAU ) THEN
            WRITE (IFILTAU,99006) NTHOST + NTCLU,NQHOST + NQCLU
            WRITE (IFILTAU,99007) (IT,IQAT(1,IT),TXT_T(IT),IT=1,NTHOST+
     &                            NTCLU)
            DO IQ = 1,NQHOST + NQCLU
               WRITE (IFILTAU,'(''site'',I3,'' of '',A)') IQ,
     &                SYSTEM(1:LSYSTEM)
            END DO
         ELSE
            WRITE (IFILTAU,99008) NQHOST + NQCLU
            DO IQ = 1,NQHOST + NQCLU
               WRITE (IFILTAU,'(''site'',I3,'' of '',A)') IQ,
     &                SYSTEM(1:LSYSTEM)
            END DO
         END IF
      END IF
C
      WRITE (6,99018) SWS,ALAT,NQ,NL,NLM,NQ*NLM,NK,NKM,NKKR,NLIN,
     &                POTFMTINP,POTFMTOUT,RDTAU,RDTAUMQ,WRTAU,WRTAUMQ,
     &                (PLOTPRS(I),I=1,3),(PLOT2DPR(I),I=1,2),IPRINT,
     &                ITEST,WRMAT,INITELOOP,RUNELOOP,SPLITSS,
     &                NO_SYMMETRY,NO_SYMMETRY_CLU,NO_SYMMETRY_LINRESP,
     &                SYMMETRIZE_MSS,SYMMETRIZE_TAU,SYMMETRIZE_RHO,
     &                SYMMETRIZE_POT,LLOYD,GF_CONV_RH,NONMAG,
     &                (FULLPOT .OR. FULLPOT5),FULLCHARGE,TASK,
     &                SYSTEM_DIMENSION,SYSTEM_TYPE,SUB_SYSTEM,KKRMODE
C
      IF ( ITEST.NE.0 ) WRITE (6,99009) TXTTEST(ITEST)
C
      RASRAD = 0.4D0
      RASSCL = 3
C
      IF ( IPRINT.GT.0 ) CALL SYMPLOTUC(RASRAD,RASSCL)
C
C ======================================================================
C                 parameter for  energy - path
C ======================================================================
C
      CALL INIT_MOD_ENERGY(IPRINT,TASK,ITEST,NCPA,INITELOOP,SCFSTATUS,
     &                     RDTAU,RDTAUMQ,NOWRDOS,NETAU)
C
Csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C                         structure constants
Csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      IF ( (INITELOOP .OR. (TASK(1:5).EQ.'EKREL') .OR. 
     &     (TASK(1:3).EQ.'BSF') .OR. (TASK(1:7).EQ.'COMPTON')) .AND. 
     &     TASK(1:6).NE.'SOCPAR' .AND. TASK(1:6).NE.'PSHIFT' .AND. 
     &     TASK(1:5).NE.'FSCAT ' .AND. KBZI.NE.0 ) THEN
C
         IF ( KBZI.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'structure constants requested but KBZI=0')
C
         IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
            ETOP = MAX(ABS(EMIN),ABS(EMAX))
C
            CALL STRINIT(ETOP)
C
            IF ( ITEST.EQ.8 ) CALL STRTEST(ETOP)
C
         ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C
            CALL INIT_MOD_TBCLU(IPRINT)
C
            CALL INIT_TBCALC(IPRINT)
C
         ELSE
C
            CALL STOP_MESSAGE(ROUTINE,'KKRMODE not found')
C
         END IF
C
      END IF
C
C=======================================================================
C                            k-space integration
C=======================================================================
C
      CALL INIT_MOD_KSPACE(INITELOOP,MOL,KMROT,ITEST,NKTABMAX,NELMTMAX,
     &                     SNTAUUVMAX)
C
C=======================================================================
C                     parameter for CPA - iteration
C=======================================================================
C
      IF ( INITELOOP ) CALL INIT_MOD_CPA(NCPA)
C
C=======================================================================
C      calculate rotation matrices  MROTQ  and  DROTQ
C      connecting LOCAL and GLOBAL frame individually for
C      all sites IQ according to the flag  KMROT
C=======================================================================
C
      CALL INIT_MOD_SITES_ROTMAG
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C                        file-info 2
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      OPEN (UNIT=13,STATUS='SCRATCH')
      OPEN (UNIT=14,STATUS='SCRATCH')
      OPEN (UNIT=15,STATUS='SCRATCH')
C
      WRITE (6,FMT='(/,1X,79(''f''),//,10X,''FILES:'',/)')
      WRITE (6,'(10X,A,A)') 'potential :  ( 4) ',POTFIL(1:LPOTFIL)
C
C ======================================================================
C                           open DOSFIL
C              open as scratch - file if DOS is not involved
C ======================================================================
C
      IF ( NOWRDOS .AND. .NOT.RDDOS ) THEN
         OPEN (UNIT=IFILDOS,STATUS='SCRATCH')
      ELSE IF ( RDDOS ) THEN
         OPEN (UNIT=IFILDOS,FILE=DOSFIL)
      ELSE
         IF ( MPI_ID.EQ.0 ) THEN
            OPEN (UNIT=IFILDOS,FILE=DOSFIL)
         ELSE
            OPEN (UNIT=IFILDOS,STATUS='SCRATCH')
         END IF
C
         CALL WRHEAD(IFILDOS,DOSFIL,'DOS       ',NETAB(1))
         WRITE (IFILDOS,99019) 'DOS-FMT:  ','OLD-SPRKKR'
C
      END IF
C
      IF ( WRNOS ) THEN
         NOSFIL = DATSET(1:LDATSET)//'.nos'
         OPEN (UNIT=IFILNOS,FILE=NOSFIL)
      END IF
C
      IF ( TASK(1:5).NE.'EKREL' ) THEN
         IF ( .NOT.NOWRDOS .OR. TASK(1:5).EQ.'APS  ' .OR. TASK(1:5)
     &        .EQ.'AES  ' .OR. TASK(1:5).EQ.'NRAES' .OR. TASK(1:5)
     &        .EQ.'RELAX' ) WRITE (6,'(10X,A,A)') 'DOS       :  (10) ',
     &                             DOSFIL
C
         IF ( RDTAU .OR. WRTAU .OR. RDTAUMQ .OR. WRTAUMQ ) THEN
            IF ( (TASK(1:3).NE.'T1 ') ) WRITE (6,'(10X,A,A)')
     &            'TAU       :  ( 9) ',TAUFIL
         END IF
C
         IF ( IREL.EQ.3 ) THEN
            CALL INPUT_FIND_SECTION('CONTROL',0)
            IF ( FOUND_SECTION ) THEN
               LF = LDATSET + 4
               CALL SECTION_FIND_KEYWORD('WRKAPDOS',WRKAPDOS)
               IF ( WRKAPDOS ) THEN
                  FILNAM = DATSET(1:LDATSET)//'.kap'
                  CALL WRHEAD(13,FILNAM,'DOS-KAPPA ',NETAB(1))
                  WRITE (6,'(10X,A,A)') 'DOS-kappa :  (13) ',
     &                                  FILNAM(1:LF)
               END IF
               CALL SECTION_FIND_KEYWORD('WRPOLAR',WRPOLAR)
               IF ( WRPOLAR ) THEN
                  FILNAM = DATSET(1:LDATSET)//'.pol'
                  CALL WRHEAD(14,FILNAM,'POLAR     ',NETAB(1))
                  WRITE (6,'(10X,A,A)') 'polaris.  :  (14) ',
     &                                  FILNAM(1:LF)
               END IF
               CALL SECTION_SET_STRING('DOSREP',DOSREP,'9999',0)
               IF ( DOSREP.NE.'XXX' ) WRITE (6,'(10X,A,A)')
     &               'DOS-REP   :       ',DOSREP
               CALL SECTION_SET_STRING('DOSCORSYS',DOSCORSYS,'9999',0)
               IF ( DOSCORSYS.EQ.'GLO' ) WRITE (6,'(10X,A,A)')
     &               'DOS CORSYS:       ',DOSCORSYS
               IF ( DOSCORSYS.NE.'GLO' .AND. DOSCORSYS.NE.'LOC' ) THEN
                  WRITE (6,*) 'Illegal DOSCORSYS value:',DOSCORSYS
                  CALL STOP_MESSAGE(ROUTINE,'DOSCORSYS')
               END IF
            END IF
         END IF
C
      ELSE
C
         IF ( NCPA.NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'TASK=EKREL not possible for disordered alloys  (NCPA<>0)'
     &        )
C
         FILNAM = DATSET(1:LDATSET)//'.bnd'
         CALL WRHEAD(IFILDOS,FILNAM,'DISPERSION',NETAB(1))
C
         WRITE (6,'(10X,A,A)') 'E(k)      :  (10) ',FILNAM
C
      END IF
C
      IF ( TASK(1:3).EQ.'BSF' ) THEN
         FILNAM = DATSET(1:LDATSET)//'.bsf'
         CALL WRHEAD(11,FILNAM,'BSF       ',NETAB(1))
C
         WRITE (6,'(10X,A,A)') 'Bloch-SF  :  (11) ',FILNAM
      END IF
C
      IF ( WRCPA ) WRITE (6,'(10X,A,A)') 'CPA       :  ',CPAFIL
C
      WRITE (6,*) ' '
C
      CLOSE (UNIT=4)
C
C=======================================================================
C
      IF ( IBZINT.NE.0 .AND. KBZI.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'TAU via BZ-integration but KBZI=0')
C
C=======================================================================
C                initialize FULL POTENTIAL type calculations
C=======================================================================
C
      CALL INIT_MOD_TYPES_FULLPOT(FULLPOT5)
C
C=======================================================================
C
C--------- in case of VBPES going beyond single scattering approximation
C------------------- reduce NL temporarily to deal with the valence band
C
      IF ( USETAUNN .OR. USELPLS1 ) THEN
         NL = 0
         NKKR = 0
         DO IQ = 1,NQ
            NLQ(IQ) = NLQ(IQ) - 1
            NL = MAX(NL,NLQ(IQ))
            NKMQ(IQ) = 2*NLQ(IQ)**2
            NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
            NKKR = NKKR + 2*NLQ(IQ)**2
            IF ( IQ.EQ.1 ) THEN
               IND0Q(IQ) = 0
            ELSE
               IND0Q(IQ) = IND0Q(IQ-1) + 2*NLQ(IQ-1)**2
            END IF
         END DO
C
         NLM = NL**2
         NKM = 2*NLM
C
         IF ( FULLPOT ) THEN
C
            IF ( MPI_ID.EQ.0 ) THEN
C
               CALL FPCOUPL(IPRINT,NKM,NL)
C
               REWIND (IOTMP)
               READ (IOTMP) MC,IBLKTOP,J1TOP
C
               READ (IOTMP) ((((((VAMEG(I,J,IPOT,J1,IBLK,IT),I=1,MC),J=1
     &                      ,MC),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,
     &                      IBLKTOP),IT=1,NT)
               IF ( IREL.GE.3 ) READ (IOTMP)
     &                                ((((((VAMEF(I,J,IPOT,J1,IBLK,IT),
     &                                I=1,MC),J=1,MC),IPOT=1,NLMFPMAX),
     &                                J1=1,J1TOP),IBLK=1,IBLKTOP),IT=1,
     &                                NT)
               CLOSE (IOTMP)
C
            END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
            IF ( MPI ) THEN
C
               CALL DRV_MPI_BARRIER
C
               IDIMS = NCPLWFMAX*NCPLWFMAX*NLMFPMAX*2*NKMMAX*NTMAX
               CALL DRV_MPI_BCAST_C(0,VAMEG(1,1,1,1,1,1),IDIMS)
               IF ( IREL.GE.3 )
     &              CALL DRV_MPI_BCAST_C(0,VAMEF(1,1,1,1,1,1),IDIMS)
C
               CALL DRV_MPI_BCAST_L(0,LHS_SOL_EQ_RHS_SOL,1)
C
            END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
         END IF
      END IF
C
C
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         NLT(IT) = NLQ(IQAT(1,IT))
         NLIN_T(IT) = NLINQ(IQ)
         NKM_T(IT) = NKMQ(IQ)
      END DO
C
C=======================================================================
C
C=======================================================================
C                 open scratch file for wave functions
C=======================================================================
C
C --------------------------------------- record length for CBWF - files
      IF ( FULLPOT .OR. FULLPOT5 .OR. SPHERCELL ) THEN
C
         RECLNGWF = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
         RECLNGWF_SPH = (2*NKMMAX**2+2*NRMAX*NKMAX)*2*LRECREAL8
C
         IF ( IREL.EQ.2 ) THEN
C
            RECLNGWF = RECLNGWF*2
            RECLNGWF_SPH = RECLNGWF_SPH*2
C
         ELSE IF ( IREL.EQ.3 ) THEN
C
            RECLNGWF_SPH = (4*NKMMAX**2+NKMMAX+2*NKMMAX+
     &                     4*NRMAX*2*NKMMAX)*2*LRECREAL8
C
         END IF
C
      ELSE
C
         RECLNGWF = (5+2*(2+2*NRMAX*2))*LRECREAL8
         RECLNGWF_SPH = 1
C
      END IF
C
C--------------------------------------- file numbers set in <MOD_FILES>
      OPEN (IFILCORWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      OPEN (IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      CALL SET_IFIL_LHS(IFILCBWF,IFILCBWF_LHS)
      OPEN (IFILCBWF_LHS,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      CALL SET_IFIL_SPH(IFILCBWF,IFILCBWF_SPH)
      OPEN (IFILCBWF_SPH,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
      OPEN (IFILGFWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
      OPEN (IFILLDAU,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
      OPEN (IFILMEZZL,STATUS='SCRATCH',FORM='UNFORMATTED')
C
      IF ( DUMP_MODULES ) THEN
         CALL OPEN_IOTMP_FILE(PROGNAME,IOTMP,'dump_5.0')
C      INCLUDE 'dump_modules_insert'
         CLOSE (IOTMP)
         STOP
      END IF
C
C=======================================================================
C                 check input for compatibility
C=======================================================================
C
      CALL CHECK_INPUT((FULLPOT .OR. FULLPOT5),POSANIPREP)
C
C=======================================================================
C                      write builbot info
C=======================================================================
C
      IF ( WRBUILDBOT ) THEN
         OPEN (IFILBUILDBOT,FILE='buildbot.log')
         CALL WRDATE(PROGNAME(1:LPROGNAME),IFILBUILDBOT)
      END IF
C
C=======================================================================
C
      IF ( PROGNAME(1:6).EQ.'KKRSCF' ) THEN
C
         CALL SCF_DRIVE(KBZI,KMOL,FULLPOT5)
C
C=======================================================================
C
      ELSE IF ( PROGNAME(1:6).EQ.'KKRGEN' ) THEN
C
         IF ( TASK.EQ.'TUTORIAL   ' ) THEN
C
            CALL TUTORIAL(RUNELOOP,TASK,DOSREP)
C
         ELSE
C
            CALL GEN(RUNELOOP,TASK,DOSREP,DOSCORSYS)
C
         END IF
C
C=======================================================================
C
      ELSE IF ( PROGNAME(1:6).EQ.'KKRCHI' ) THEN
C
         CALL CHIDRIVE(KBZI,IPRINT,WRMAT,TASK)
C
C=======================================================================
C
      ELSE IF ( PROGNAME(1:3).EQ.'CLU' ) THEN
C
         CALL CLUDRIVE(DOSREP,DOSCORSYS)
C
C=======================================================================
C
      ELSE IF ( PROGNAME(1:6).EQ.'KKROPM' ) THEN
C
         CLOSE (IFILCORWF)
C
         CALL OPM(KBZI)
C
C=======================================================================
C
      ELSE IF ( PROGNAME(1:7).EQ.'KKRSPEC' ) THEN
C
         CALL SPECMAIN
C
C=======================================================================
C
      END IF
C
      WRITE (6,*) 'KKR-run finished successfully'
      STOP
C
C=======================================================================
C
 100  CONTINUE
      WRITE (6,99021) POTFIL(1:LPOTFIL)
      CALL STOP_MESSAGE(ROUTINE,
     &       'the requested potential file could not be opened  - check'
     &       )
C
C=======================================================================
99001 FORMAT (I6,I4,I4,1X,3I5,F7.4,I3,I4,18X,I4,I5,1X,3I3,:,/,(67X,3I3))
99002 FORMAT (/,1X,79('#'),/,10X,'TROUBLE:  NLQ set inconsistently',
     &        ' for equivalent sites',/,10X,'IQ =',I3,'    NLQ =',I3,/,
     &        10X,'JQ =',I3,'    NLQ =',I3,/)
99003 FORMAT (/,1X,79('#'),/,10X,'STOP in program  ',A,/,10X,
     &        'no or empty input file supplied',/,10X,
     &        'proper program call:',6X,A,'<  inputfile',/,1X,79('#'),/)
99004 FORMAT (/,1X,79('*'),/,10X,'TAU will be read from file ',/)
99005 FORMAT (/,1X,79('*'),/,10X,'TAU - file is empty !!!!   ==>>   ',
     &        'TAU will be created !!!!',/,1X,79('*'),/)
99006 FORMAT (/,70('*'),/,10X,' TAU(LAM,LAM'')',/,70('*'),/,5X,' NT =',
     &        I5,' NQ =',I5,' FMT = 3')
99007 FORMAT (5X,' IT =',I5,' IQ =',I5,': ',A)
99008 FORMAT (/,70('*'),/,10X,' TAU(LAM,LAM'') and M(LAM,LAM'')',/,
     &        70('*'),/,13X,' NQ =',I5,' FMT = 3')
99009 FORMAT (2(/,1X,79('*')),/,10X,'test option selected:',/,10X,A40,
     &        2(/,1X,79('*')),/)
99010 FORMAT (//,'   number of lattice sites (NQ):',I4,//,
     &        '   site NL mesh MGAM MTET MPHI   RWS  CPA NOQ',
     &        ' IT ELMT VAL CONC QREP SYMG eq.sit.')
99011 FORMAT (3X,'average',19X,F8.4,10X,'total',F8.2,/)
99012 FORMAT (2(/,1X,79('#')),/,10X,'WARNING from   ',A,/,10X,
     &        'radial mesh parameters inconsistent '/,10X,
     &        'a1*(a2xa3) * a^3     ',F15.8,/,10X,
     &        'NQ * SWS^3 * 4*PI/3  ',F15.8,/,10X,
     &        'DEVIATION            ',F15.8,/,10X,
     &        'radial mesh parameters will be modified',2(/,1X,79('#')),
     &        /)
99013 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*        ****   *****   *****   *    *  *    *  *****        *'
     &  ,/,10X,
     &  '*       *    *  *    *  *    *  *   *   *   *   *    *       *'
     &  ,/,10X,
     &  '*       *       *    *  *    *  *  *    *  *    *    *       *'
     &  ,/,10X,
     &  '*        ****   *****   *****   * *     * *     *****        *'
     &  ,/,10X,
     &  '*            *  *       *  *    ** *    ** *    *  *         *'
     &  ,/,10X,
     &  '*       *    *  *       *   *   *   *   *   *   *   *        *'
     &  ,/,10X,
     &  '*        ****   *       *    *  *    *  *    *  *    *       *'
     &  ,/,10X,'*',60X,'*',/,10X,'*',60('-'),'*',/,10X,'*   ',A10,
     &  '  VERSION  7.7.0          (C) 2017  H. Ebert   *',/,10X,62('*')
     &  ,//)
99014 FORMAT (10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*                  PUBLIC VERSION                            *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99015 FORMAT (3X,'KMROT: ',I5,3X,A,3X,'MOMENTS_ROTATED: ',L2,/,3X,
     &        'QMVEC: ',3F7.3,/)
99016 FORMAT (//,'   number of atomic types  (NT):',I4,//
     &  '   type TXTT     NL mesh    RMT     RWS   NAT  CONC   on sites'
     &  )
99017 FORMAT (I6,2X,A,I3,I4,F9.4,F8.4,I4,F7.3,3X,8I3,:,/,(50X,8I3))
99018 FORMAT (//,10X,'SWS =',F10.6,5X,'ALAT =',F10.6,//,10X,'NQ =',I3,
     &        5X,'NL =',I3,5X,'NLM =',I4,3X,'NQ*NLM =',I4,/,22X,'NK =',
     &        I3,5X,'NKM =',I4,5X,'NKKR =',I4,5X,'NLIN =',I4,//,10X,
     &        'POTFMTINP =',I2,5X,'POTFMTOUT =',I2,/,10X,'RDTAU     =',
     &        L2,5X,'RDTAUMQ   =',L2,/,10X,'WRTAU     =',L2,5X,
     &        'WRTAUMQ   =',L2,/,10X,'PLOTPOT   =',L2,5X,'PLOTRHO   =',
     &        L2,5X,'PLOTSFN   =',L2,/,10X,'PLOT2DPOT =',L2,5X,
     &        'PLOT2DRHO =',L2,/,10X,'IPRINT    =',I2,5X,'ITEST     =',
     &        I2,5X,'WRMAT     =',L2,/,10X,'INITELOOP =',L2,5X,
     &        'RUNELOOP  =',L2,5X,'SPLITSS   =',L2,/,10X,
     &        'NO_SYM set for:',/,10X,'HOST      =',L2,5X,'CLU       =',
     &        L2,5X,'LINRESP   =',L2,/,10X,'symmetrize:',/,10X,
     &        'MSS       =',L2,5X,'TAU       =',L2,5X,'RHO       =',L2,
     &        5X,'POT       =',L2,/,10X,'LLOYD     =',L2,/,10X,
     &        'GF_CONV_RH=',L2,/,10X,'NONMAG    =',L2,/,10X,
     &        'FULLPOT   =',L2,/,10X,'FULLCHARGE=',L2,//,10X,
     &        'TASK      = ',A,//,10X,'SYSDIM    = ',A,/,10X,
     &        'SYSTYPE   = ',A,/,10X,'SUBSYSTEM = ',A,/,10X,
     &        'KKRMODE   = ',A)
99019 FORMAT (A10,A,A)
99020 FORMAT (//,10X,'PROGNAME       : ',A,/,10X,'SCFSTATUS      : ',A,
     &        /,10X,'SCFSTATUS_HOST : ',A,/)
99021 FORMAT (//,10X,'POTFIL : ',A,/)
      END
