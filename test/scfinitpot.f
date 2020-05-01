C*==scfinitpot.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
C
      SUBROUTINE SCFINITPOT(WRSSDOS,QEL,TAUQ,TAUT,NOCHKVAL)
C   ********************************************************************
C   *                                                                  *
C   *   prepare the SCF-cycle starting from  scratch                   *
C   *                                                                  *
C   *   - set up a guess for the charge and spin density               *
C   *     RHOCHR  and  RHOSPN, respectively, using atomic data         *
C   *     and the Mattheiss construction                               *
C   *   - set up the potential functions  VT and BT  using <SCFNEWPOT> *
C   *     default:    use the ionicity from the Mattheiss construction *
C   *                 QION may be specified in the input instead       *
C   *     USEVMATT:   get VT as superposition of atomic potentials     *
C   *                                                                  *
C   *   - calculate the DOS along a straight path along the real       *
C   *     axis to get a reasonable guess for the Fermi energy  EFERMI  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_FILES,ONLY:TITLE,LDATSET,DATSET,IPRINT,WRMAT,IFILDOS,
     &    IFILCBWF,FOUND_SECTION,FOUND_REAL,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_ENERGY,ONLY:NETAB,NEMAX,WETAB,ETAB,EFERMI,EMIN
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,RWS
      USE MOD_LATTICE,ONLY:ALAT,ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,
     &    ADAINV_R,ADAINV_I,SYSTEM_TYPE,SYSTEM_DIMENSION
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_CALCMODE,ONLY:IREL,ORBPOL,LLOYD,SEMICORE,MOL,NONMAG,SOLVER
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,NL,NMEMAX,NKMMAX,NLMAX,NKMQ,NKM
      USE MOD_TYPES,ONLY:NT,NTMAX,VT,CONC,Z,NAT,IMT,NLT,NT_L,ITBOT,
     &    ITTOP,NVALTOT,BT,CMNTT,NLMFPMAX,DOBS_LTX,DOBS_TX,OBS_LTX,
     &    OBS_TX,DOBS_TX_GLO,OBS_TX_GLO
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,IQAT,QBAS,NQ_L,NQ_R,NQ_I,
     &    IQBOT,IQTOP,AVMAD
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SCF,ONLY:SCFSTATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFINITPOT')
      REAL*8 DPAS
      PARAMETER (DPAS=0.05D0)
      INTEGER NRAT,NELIM
      PARAMETER (NRAT=251,NELIM=800)
      REAL*8 ESTEPTOL,EIMAG,EWINSSDOS,NOSTOL,DEFTOL
      PARAMETER (ESTEPTOL=1D-8,EIMAG=0.1D0,EWINSSDOS=2.0D0,NOSTOL=1D-4,
     &           DEFTOL=0.02D0)
      LOGICAL USE_EFERMILD
      PARAMETER (USE_EFERMILD=.FALSE.)
C
C Dummy arguments
C
      LOGICAL NOCHKVAL,WRSSDOS
      REAL*8 QEL(NTMAX)
      COMPLEX*16 TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BCOR(:),BCORS(:),CPACHNG,CPACHTAB(:),DNSINI(:,:,:,:),DOS1,
     &       DOS2,DQ,DQ0,ECORTAB(:,:),ECTOP,EFERMIDOS,EFERMILD,
     &       EFERMI_GUESS,EILOW,EINCR,EMAX,ESTEP,ETOT,MSPIN(:),NOSLD,
     &       QION(:),QIONINP(:),QIONSCL,R1AT(:),RWK1(:,:,:,:),RWK2(:),
     &       SUMAT(:,:),SUMDOS,TDOSINI(:),TOTDOS,TOTNOS,VMTZ,VR2AT(:,:),
     &       WA(:),WSPIN,WT,X
      LOGICAL CALCINT,CHECK,EFERMI_GUESS_SUPPLIED,GETIRRSOL,USEIONMATT,
     &        USEQION,USEVMATT
      COMPLEX*16 CSUMA,CSUMSS,CTOTDOS(:),CTOTNOS(:),CWKE(:),
     &           CWKKMTE(:,:,:),EBOT,ERYD,ETAB0(:),ETOP,ETRY,
     &           MEZJ(:,:,:,:),MEZZ(:,:,:,:),MSSQ(:,:,:),MSST(:,:,:),
     &           NOSBOT,NOSTOP,NOSTRY,P,PCOEF(:,:),PHASA(:,:,:),PHASK(:)
     &           ,PHAST(:,:,:),SEBT(:,:),SEVT(:,:),SSST(:,:,:),
     &           TSSQ(:,:,:),TSST(:,:,:),WETAB0(:)
      CHARACTER*80 DOSFIL
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IE,IE1,IECPAFAIL(:),IECURR,
     &        IELD,IELOOP,IETOP,IL,IM,IO,IQ,IQBOTTMP,IQTOPTMP,IR,IRTOP,
     &        IT,ITBOTTMP,ITCPA,ITRY,ITTOPTMP,IWRIRRWF,IWRREGWF,JRAT(:),
     &        LWK1,LWK2,LWKE,LWKKMTE,N,NCPAFAIL,NELD,NEMAX0,NESTART,
     &        NETAB0,NOS0,NPAD,NWRLOG
      CHARACTER*10 SOLVER_SELECTED
      REAL*8 TABESTMOM,YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MSPIN,QION,BCOR,CPACHTAB
      ALLOCATABLE TSSQ,MSSQ,MSST,TSST,JRAT,VR2AT
      ALLOCATABLE IECPAFAIL,BCORS,SUMAT
      ALLOCATABLE WA,ECORTAB
      ALLOCATABLE QIONINP,R1AT,SEBT,SEVT
      ALLOCATABLE CTOTDOS,PHASK,ETAB0,WETAB0,PHASA,PHAST,DNSINI,TDOSINI
      ALLOCATABLE RWK1,RWK2,CTOTNOS,PCOEF
      ALLOCATABLE CWKE,CWKKMTE,MEZJ,MEZZ,SSST
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (CWKE(NPROCS))
      ALLOCATE (CWKKMTE(NKMMAX,NTMAX,NPROCS))
      ALLOCATE (MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (SSST(NKMMAX,NKMMAX,NTMAX))
C
      ALLOCATE (RWK1(NLMAX,NTMAX,2,NPROCS),RWK2(NPROCS))
      ALLOCATE (SEBT(NRMAX,NTMAX),SEVT(NRMAX,NTMAX))
      ALLOCATE (PHAST(NKMMAX,NTMAX,NELIM),CTOTNOS(NELIM))
      ALLOCATE (PHASA(NKMMAX,NTMAX,NELIM),PHASK(NELIM),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ETAB0')
C
      ALLOCATE (MSPIN(NTMAX),QION(NTMAX))
      ALLOCATE (BCOR(NTMAX),CPACHTAB(NELIM),CTOTDOS(NELIM))
      ALLOCATE (ETAB0(NELIM),WETAB0(NELIM))
      ALLOCATE (DNSINI(NLMAX,NTMAX,2,NELIM),TDOSINI(NELIM))
      ALLOCATE (TSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSSQ(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MSSQ')
C
      ALLOCATE (TSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (MSST(NKMMAX,NKMMAX,NTMAX),WA(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MSST')
C
      ALLOCATE (JRAT(NTMAX),BCORS(NTMAX),R1AT(NTMAX))
      ALLOCATE (SUMAT(NRAT,NTMAX),ECORTAB(120,NTMAX),QIONINP(NTMAX))
      ALLOCATE (IECPAFAIL(NELIM),VR2AT(NRAT,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: BCORS')
C
      WRITE (6,99023)
C
C-----------------------------------------------------------------------
C  use fast (may be less accurate) radial differential equation SOLVER
C  to speed up the Fermi energy search
C-----------------------------------------------------------------------
      SOLVER_SELECTED = SOLVER
      SOLVER = 'BS'
      WRITE (6,99001) SOLVER
C-----------------------------------------------------------------------
C
      CHECK = .FALSE.
C     CHECK = .TRUE.
      IF ( IREL.GT.1 ) THEN
         WSPIN = 1D0
      ELSE
         WSPIN = 2D0
      END IF
C
      CALL CINIT(NRMAX*NTMAX,SEBT)
      CALL CINIT(NRMAX*NTMAX,SEVT)
      CALL RINIT(120*NTMAX,ECORTAB)
C
      DO IT = 1,NT
         QIONINP(IT) = 0D0
         MSPIN(IT) = TABESTMOM(Z(IT))
      END DO
C
      USEVMATT = .FALSE.
      USEQION = .FALSE.
      USEIONMATT = .TRUE.
      QIONSCL = 1.0D0
C
      CALL INPUT_FIND_SECTION('SCF',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL_ARRAY('QION',QIONINP,N_FOUND,NTMAX,0,
     &                               9999D0,0)
         USEQION = FOUND_REAL_ARRAY
         IF ( USEQION ) USEIONMATT = .FALSE.
         CALL SECTION_SET_REAL_ARRAY('MSPIN',MSPIN,N_FOUND,NTMAX,0,
     &                               9999D0,0)
         CALL SECTION_FIND_KEYWORD('USEVMATT',USEVMATT)
         IF ( USEVMATT ) THEN
            USEQION = .FALSE.
            USEIONMATT = .FALSE.
         END IF
C
         CALL SECTION_SET_REAL('QIONSCL',QIONSCL,9999D0,0)
      END IF
      IF ( NONMAG ) MSPIN(1:NTMAX) = 0D0
C
      IF ( USEVMATT .AND. SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
         WRITE (6,99003)
         USEVMATT = .FALSE.
      END IF
C
C ======================================================================
C               calculate the starting charge and spin density
C ======================================================================
C
C-----------------------------------------------------------------------
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
         NT_L = 0
         DO IQ = 1,NQ_L
            DO IO = 1,NOQ(IQ)
               NT_L = MAX(NT_L,ITOQ(IO,IQ))
            END DO
         END DO
C
         DO IQ = NQ_L + 1,NQ
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               IF ( IT.LT.NT_L ) CALL STOP_MESSAGE(ROUTINE,'IT.LT.NT_L')
            END DO
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      CALL SCFINITRHO(USEIONMATT,USEQION,IPRINT,QEL,QION,QIONINP,MSPIN,
     &                R1AT,VR2AT,QIONSCL)
C
      CMNTT(1:NLMFPMAX,ITBOT:ITTOP) = 0D0
C
C ======================================================================
C                     calculate the starting potential
C ======================================================================
C
C HE NOTE: for LIR and LIV the potential for ALL types is set
C HE NOTE: however the proper 2D Madelung matrix AVMAD for the I-ZONE
C HE NOTE: is not available, but only that of the 3D L-BULK only.
C HE NOTE: this should have 0-s in the LI-, II- and RI-blocks
C HE NOTE: Although the Madelung term for the I-Zone is still missing
C HE NOTE: the guess potential for the I-Zone seems to be acceptable
C
      ITBOTTMP = ITBOT
      ITTOPTMP = ITTOP
      IQBOTTMP = IQBOT
      IQTOPTMP = IQTOP
      ITBOT = 1
      ITTOP = NT
      IQBOT = 1
      IQTOP = NQ
C
      CALL SCFNEWPOT(.FALSE.,ETOT,ECORTAB)
C
      IF ( NONMAG ) BT(1:NRMAX,ITBOT:ITTOP) = 0D0
C
      ITBOT = ITBOTTMP
      ITTOP = ITTOPTMP
      IQBOT = IQBOTTMP
      IQTOP = IQTOPTMP
C
C-----------------------------------------------------------------------
C  for  SYSTEM_TYPE = VIV i.e. slab  set VT = 0 for left and right host
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
C
         DO IQ = 1,NQ_L + NQ_I + NQ_R
            IF ( IQ.GT.NQ_L .AND. IQ.LE.(NQ_L+NQ_I) ) CYCLE
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               VT(:,IT) = 0.0D0
               BT(:,IT) = 0.0D0
            END DO
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C     for  2D LIR- (bulk) and  LIV (surface) calculations
C      -  reset the Madelung matrix to that of the left host
C      -  modify the left host accordingly   IQ=1,NQ_L
C-----------------------------------------------------------------------
C
C HE NOTE: as long as the 2D Madelung matrix AVMAD is still missing
C HE NOTE: (see note above) this reset of the L-BULK potential
C HE NOTE: could be omitted
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &     (SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV') )
     &     THEN
         IF ( IQTOP.NE.NQ_L )
     &         CALL STOP_MESSAGE(ROUTINE,'after reset IQTOP <> NQ_L')
C
         CALL SCFMAD3D(IPRINT,ALAT,ABAS_L,QBAS,AVMAD,NQ_L,NQMAX)
C
         CALL SCFNEWPOT(.FALSE.,ETOT,ECORTAB)
C
         IF ( NONMAG ) BT(1:NRMAX,ITBOT:ITTOP) = 0D0
C
      END IF
C-----------------------------------------------------------------------
C
C--------------------- create potential using the Mattheiss construction
C----------------------- overwrite the potential based on RHO(Mattheiss)
C------------------------- keep the B-field BT set up from atomic charge
C
      IF ( USEVMATT ) THEN
C
         CALL SCF0MATT(MOL,NQ,NQ_L,NQ_R,NQMAX,QBAS,ABAS,ABAS_L,ABAS_R,
     &                 ABAS_I,ADAINV_L,ADAINV_R,ADAINV_I,
     &                 SYSTEM_DIMENSION,ALAT,DPAS,SUMAT,VR2AT,R1AT,JRAT,
     &                 RWS,CONC,NOQ,IMT,IQAT,ITOQ,NT,NTMAX,NRAT,WA,
     &                 NRMAX)
C
         DO IT = 1,NT
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            DO IR = 1,JRAT(IT)
               WA(IR) = R1AT(IT)*EXP(DBLE(IR-1)*DPAS)
            END DO
C
            DO IR = 1,IRTOP
               VT(IR,IT) = YLAG(R(IR,IM),WA,SUMAT(1,IT),0,3,JRAT(IT))
     &                     /R(IR,IM)**2
            END DO
         END DO
C
         WRITE (6,99011)
      ELSE IF ( USEIONMATT ) THEN
         WRITE (6,99012)
      ELSE
         WRITE (6,99013)
      END IF
C
      CALL VMUFTIN(VMTZ)
C
C-----------------------------------------------------------------------
C     for  2D LIR- (bulk) and  LIV (surface) calculations
C     the muffin tin shift VMTZ is determined by the left bulk
C     perform the SAME muffin tin shift VMTZ of the potential also
C     for the atoms in the interaction zone
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &     (SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV') )
     &     THEN
C
         DO IQ = NQ_L + 1,NQ_L + NQ_I
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
               VT(1:IRTOP,IT) = VT(1:IRTOP,IT) - VMTZ
               VT((IRTOP+1):NRMAX,IT) = 0.0D0
            END DO
         END DO
C
      END IF
C
C ======================================================================
C         write out the starting potential - EFERMI not known yet
C ======================================================================
C
      CALL POTWRSPR
C
C ======================================================================
C                    check the number of valence electrons
C ======================================================================
C
      IF ( .NOT.SEMICORE .AND. .NOT.NOCHKVAL )
     &     CALL SCFCHKNVAL(EMIN,ECTOP)
C
C ======================================================================
C                      find the Fermi energy EFERMI
C                 or just write the SS-DOS for checking
C ======================================================================
C
      IF ( IPRINT.GE.1 ) THEN
         NWRLOG = 1
      ELSE
         NWRLOG = 0
      END IF
C
      DOSFIL = DATSET(1:LDATSET)//'_SS-DOS.dos'
C
      IF ( WRSSDOS .AND. MPI_ID.EQ.0 ) THEN
C
         OPEN (UNIT=IFILDOS,FILE=DOSFIL)
C
         TITLE = 'DOS for SCF-start'
         CALL WRHEAD(IFILDOS,DOSFIL,'DOS       ',NELIM)
         WRITE (IFILDOS,99014) 'DOS-FMT:  ','OLD-SPRKKR'
C
         EINCR = 0.02D0
         NPROCS = INT(EWINSSDOS/EINCR) + 1
C
      ELSE
C
         OPEN (IFILDOS,STATUS='SCRATCH')
C
      END IF
C
C ----------------------------------------------------------------------
C
      DO IT = 1,NTMAX
         BCOR(IT) = 0D0
         BCORS(IT) = 0D0
      END DO
C
      IWRREGWF = 0
      IWRIRRWF = 0
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      NCPAFAIL = 0
C
C=======================================================================
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C=======================================================================
C
C--------------------------------------------------------- set up E-path
C            odd value for EINCR to recover settings in previous version
      EINCR = 1.6D0/DBLE(80-1)
      IF ( LLOYD ) EINCR = EINCR*1.5D0
      EMAX = EMIN + EINCR*(NELIM-1)
C
      CALL EPATH(3,EMIN,EMAX,EIMAG,NELIM,ETAB0,WETAB0,EILOW,IPRINT,
     &           NELIM)
C
      CALL RINIT(NLMAX*NTMAX*2*NELIM,DNSINI)
      CALL RINIT(NELIM,TDOSINI)
      CALL CINIT(NKMMAX*NTMAX*NELIM,PHASA)
      CALL CINIT(NKMMAX*NTMAX*NELIM,PHAST)
      CALL CINIT(NELIM,PHASK)
C
C ======================================================================
C      read  EFERMI_GUESS  guess for the Fermi energy
C      this will speed up the initialisation if supplied
C      in this case the energy loop and writing of DOS can be skipped
C ======================================================================
C
      CALL SECTION_SET_REAL('EFGUESS',EFERMI_GUESS,9999D0,0)
      EFERMI_GUESS_SUPPLIED = FOUND_REAL
C
      IF ( .NOT.(EFERMI_GUESS_SUPPLIED) ) THEN
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C     there are 3 cases for the follwing E-loop strcuture
C
C 1: search EFERMI, NPROCS=1 --> IELOOP is a dummy index
C                                the E-loop is realised by back-jumps
C                                to label 100
C 2: search EFERMI, NPROCS>1 --> IELOOP runs from 1 .. NPROCS
C                                each proc deals with 1 energy within
C                                this E-segment - the full E-loop is
C                                realised by back-jumps to label 100
C 3: WRSSDOS -> NPROCS=NELIM --> the IELOOP - loop is done only by the
C                                major process printing the SS-DOS
C
         IECURR = 0
 50      CONTINUE
         IE1 = IECURR + 1
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         DO IELOOP = 1,NPROCS
            IECURR = IECURR + 1
            ERYD = ETAB0(IECURR)
C
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
C--------------------------------------------------- MPI_ID = (IELOOP-1)
            IF ( (MPI_ID.EQ.(IELOOP-1)) .OR. (MPI_ID.EQ.0 .AND. WRSSDOS)
     &           ) THEN
C
C ======================================================================
C                       solve SS - differential equation
C ======================================================================
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,
     &                       GETIRRSOL,ERYD,P,IPRINT,TSST,MSST,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
C ======================================================================
               IF ( .NOT.WRSSDOS ) THEN
C
C ======================================================================
C                            CALCULATE TAU
C ======================================================================
C
                  CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
                  CALL TAU_DRIVE(IECURR,IPRINT,ERYD,P,TSSQ,MSSQ,TSST,
     &                           MSST,TAUQ,ICPAFLAG,CPACHNG,ITCPA,
     &                           ICPACONV,PHASK)
C
                  IF ( ICPAFLAG.NE.0 ) THEN
                     NCPAFAIL = NCPAFAIL + 1
                     CPACHTAB(NCPAFAIL) = CPACHNG
                     IECPAFAIL(NCPAFAIL) = IECURR
                  END IF
C                                                    END:  CALCULATE TAU
C ======================================================================
C
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
C
                  IF ( IPRINT.GE.5 .OR. WRMAT )
     &                 CALL DUMPMAT(IECURR,ERYD,6,MSST,MSSQ,TAUT,TAUQ,
     &                 .TRUE.,MEZZ,MEZJ)
C
C ======================================================================
C                     deal with Lloyd's formula
C ======================================================================
C
                  IF ( LLOYD ) CALL LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,MSST,
     &                 SSST,ERYD,IECURR)
C
C ======================================================================
C                use single site t-matrix to calculate SS-DOS
C ======================================================================
C
               ELSE
C
                  TAUT(:,:,1:NT) = TSST(:,:,1:NT)
C
               END IF
C
C ======================================================================
C
               NETAB0 = NETAB(1)
               NETAB(1) = NELIM
               NEMAX0 = NEMAX
               NEMAX = NELIM
C
               CALL CALCDOS(NWRLOG,.FALSE.,.FALSE.,1,NCPAFAIL,-1,ERYD,
     &                      MEZZ,MEZJ,TSST,MSST,TAUT,MSSQ,TAUQ,IECURR,
     &                      WETAB0(IECURR),BCOR,BCORS,DOBS_LTX,DOBS_TX,
     &                      OBS_LTX,OBS_TX,DOBS_TX_GLO,OBS_TX_GLO,
     &                      CTOTDOS)
C
               NETAB(1) = NETAB0
               NEMAX = NEMAX0
C
               IF ( IDOS.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'IDOS <> 1')
               IF ( ISMT.NE.2 ) CALL STOP_MESSAGE(ROUTINE,'ISMT <> 2')
C
               TOTDOS = 0.0D0
               DO IT = 1,NT
                  DO IL = 1,NLT(IT)
                     DNSINI(IL,IT,IDOS,IECURR)
     &                  = DIMAG(DOBS_LTX(0,IDOS,IL,IT))
                     DNSINI(IL,IT,ISMT,IECURR)
     &                  = DIMAG(DOBS_LTX(0,ISMT,IL,IT))
C
                     TOTDOS = TOTDOS + DNSINI(IL,IT,IDOS,IECURR)
     &                        *CONC(IT)*NAT(IT)
                  END DO
               END DO
               TDOSINI(IECURR) = TOTDOS
C
            END IF
C--------------------------------------------------- MPI_ID = (IELOOP-1)
         END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         IF ( WRSSDOS ) THEN
            IF ( MPI_ID.EQ.0 ) WRITE (6,99007) DOSFIL(1:(LDATSET+11)),
     &                                EMIN,EMIN + EWINSSDOS
            STOP
         END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) THEN
C
            LWK1 = NLMAX*NTMAX*2*NPROCS
C
            CALL DRV_MPI_REDUCE_R(DNSINI(1,1,1,IE1),RWK1(1,1,1,1),LWK1)
C
            LWK2 = NPROCS
C
            CALL DRV_MPI_REDUCE_R(TDOSINI(IE1),RWK2(1),LWK2)
C
            IF ( LLOYD ) THEN
C
               LWKE = NPROCS
C
               CALL DRV_MPI_REDUCE_C(CTOTDOS(IE1),CWKE(1),LWKE)
               CALL DRV_MPI_REDUCE_C(PHASK(IE1),CWKE(1),LWKE)
C
               LWKKMTE = NKMMAX*NTMAX*NPROCS
C
               CALL DRV_MPI_REDUCE_C(PHASA(1,1,IE1),CWKKMTE(1,1,1),
     &                               LWKKMTE)
               CALL DRV_MPI_REDUCE_C(PHAST(1,1,IE1),CWKKMTE(1,1,1),
     &                               LWKKMTE)
C
            END IF
C
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
         NESTART = IECURR
         TOTNOS = 0.0D0
         DO IE = 1,NESTART
            TOTNOS = TOTNOS + DREAL(WETAB0(IE))*TDOSINI(IE)
            DQ = NVALTOT/WSPIN - TOTNOS
            IETOP = IE
            IF ( MPI ) CALL DRV_MPI_BCAST_R(0,DQ,1)
            IF ( DQ.LT.0D0 ) GOTO 100
         END DO
C
         IF ( IPRINT.GT.-10 ) THEN
            IF ( IECURR.EQ.1 ) WRITE (6,99008)
            WRITE (6,99009) IECURR,ERYD,TOTNOS,TDOSINI(IETOP)
         END IF
         IF ( IECURR.LT.NELIM ) GOTO 50
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      ENERGY - LOOP    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C=======================================================================
C  the number of valence electrons TOTNOS could not be matched - STOP
C=======================================================================
         WRITE (6,99018)
         STOP
      END IF
C
C=======================================================================
C  the number of valence electrons TOTNOS has been matched - CONTINUE
C=======================================================================
C
C-----------------------------------------------------------------------
C                reset radial differential equation SOLVER
C-----------------------------------------------------------------------
 100  CONTINUE
      SOLVER = SOLVER_SELECTED
      WRITE (6,99002) SOLVER
C-----------------------------------------------------------------------
C
C=======================================================================
C       write DOS - improve EFERMI via Lloyd - new EPATH     MPI_ID.EQ.0
C
      IF ( MPI_ID.EQ.0 ) THEN
C
C=======================================================================
C -----------------------------------------------  EFERMI_GUESS_SUPPLIED
         IF ( EFERMI_GUESS_SUPPLIED ) THEN
C
            EFERMI = EFERMI_GUESS
C
C=======================================================================
C                                                      set up new E-path
            EMAX = EFERMI
C
            CALL EPATH(5,EMIN,EMAX,EIMAG,NETAB(1),ETAB(1,1),WETAB(1,1),
     &                 EILOW,IPRINT,NEMAX)
C
            WRITE (6,99020) EFERMI
C
C ------------------------------------------ .NOT. EFERMI_GUESS_SUPPLIED
         ELSE
C
            EFERMI = DREAL(ETAB0(IETOP)) + DQ/TDOSINI(IETOP)
C
            IF ( NCPAFAIL.NE.0 ) THEN
               WRITE (6,99015) CPATOL,NCPAFAIL,
     &                         (IECPAFAIL(I),DREAL(ETAB0(IECPAFAIL(I))),
     &                         CPACHTAB(I),I=1,NCPAFAIL)
               WRITE (6,99017)
            ELSE IF ( NCPA.NE.0 ) THEN
               WRITE (6,99016)
            END IF
C
C ======================================================================
C                             write  DOS
C ======================================================================
C
            DOSFIL = DATSET(1:LDATSET)//'_SCFSTART.dos'
            OPEN (UNIT=IFILDOS,FILE=DOSFIL)
C
            TITLE = 'DOS for SCF-start'
            CALL WRHEAD(IFILDOS,DOSFIL,'DOS       ',IETOP)
            WRITE (IFILDOS,99014) 'DOS-FMT:  ','OLD-SPRKKR'
C
            DO IE = 1,IETOP
               WRITE (IFILDOS,99010) ETAB0(IE),
     &                               ((((DNSINI(IL,IT,1,IE)+DNSINI(IL,
     &                               IT,2,IE))/WSPIN),IL=1,NL),
     &                               (((DNSINI(IL,IT,1,IE)
     &                               -DNSINI(IL,IT,2,IE))/WSPIN),IL=1,
     &                               NL),IT=1,NT)
            END DO
C
C ======================================================================
C            improve guess for EFERMI using the Lloyd formula
C ======================================================================
C
            IF ( LLOYD .AND. USE_EFERMILD ) THEN
C
               NELD = IETOP
C
               NPAD = NELD
               ALLOCATE (PCOEF(NPAD,NPAD),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: PCOEF')
C
C------------------------------------- remove jumps from PHAST and PHASA
C
               CALL SBRANCH(PHAST,NKM,NT,NKMMAX,NTMAX,NELD)
               CALL SBRANCH(PHASA,NKM,NT,NKMMAX,NTMAX,NELD)
C
C------------------------------------------------- sum up to get CTOTNOS
C
               SUMDOS = 0D0
C
               DO IELD = 1,NELD
                  CSUMA = C0
                  CSUMSS = C0
                  DO IT = 1,NT
                     N = NKMQ(IQAT(1,IT))
                     WT = NAT(IT)*CONC(IT)
                     DO I = 1,N
                        CSUMSS = CSUMSS - 
     &                           WT*(PHAST(I,IT,IELD)-PHASA(I,IT,IELD))
                        CSUMA = CSUMA + PHASA(I,IT,IELD)
                     END DO
                  END DO
C
                  CTOTNOS(IELD) = (PHASK(IELD)+CSUMSS)/PI
C
                  IF ( CHECK ) THEN
                     X = DREAL(ETAB0(IELD))
                     WRITE (62,'(5F12.5)') X,DIMAG(-CSUMA)
                     WRITE (63,'(2F12.5)') X,DIMAG(CSUMSS)
                     IF ( IELD.EQ.1 ) THEN
                        DOS1 = DIMAG(CTOTDOS(IELD))
                        WRITE (64,'(9F12.5)') X,0D0,0D0,0D0,NVALTOT,DOS1
                     ELSE
                        ESTEP = DREAL(ETAB0(IELD)-ETAB0(IELD-1))
                        DOS2 = DIMAG(CTOTDOS(IELD))
                        SUMDOS = SUMDOS + (DOS1+DOS2)*ESTEP/2D0
                        DOS1 = DOS2
                        NOSLD = DIMAG(CTOTNOS(IELD)-CTOTNOS(1))
                        WRITE (64,'(9F12.5)') X,NOSLD,SUMDOS,
     &                         NOSLD - SUMDOS,NVALTOT,DOS1
                     END IF
                  END IF
C
               END DO
C
C ----------------------------------------------------------------------
C                get coefficients for Pade approximation
C ----------------------------------------------------------------------
C
               CALL PADECOEFF(CTOTNOS,ETAB0,PCOEF,NPAD,NPAD)
C
C ----------------------------------------------------------------------
C                    extrapolate NOS to real axis
C ----------------------------------------------------------------------
C
               EBOT = EMIN
               CALL PADEAPPROX(NOSBOT,EBOT,ETAB0,PCOEF,NPAD,NPAD)
C
               ETOP = EFERMI
               CALL PADEAPPROX(NOSTOP,ETOP,ETAB0,PCOEF,NPAD,NPAD)
C
C ----------------------------------------------------------------------
C                   find new Fermi energy EFERMILD
C ----------------------------------------------------------------------
C
               ETRY = EFERMI
C
               DQ = NVALTOT/WSPIN - DIMAG(NOSTOP-NOSBOT)
               DQ0 = DQ
               ESTEP = SIGN(0.01D0,DQ)
C
               ITRY = 0
 110           CONTINUE
               DO I = 1,200
                  ETRY = ETRY + ESTEP
                  ITRY = ITRY + 1
C
                  CALL PADEAPPROX(NOSTRY,ETRY,ETAB0,PCOEF,NPAD,NPAD)
C
                  DQ = NVALTOT/WSPIN - DIMAG(NOSTRY-NOSBOT)
                  IF ( CHECK ) WRITE (6,99004) ITRY,DREAL(ETRY),ESTEP,
     &                                DQ,DIMAG(NOSTRY-NOSBOT)
C
                  IF ( DQ*DQ0.LT.0D0 ) THEN
                     ESTEP = -ESTEP/2D0
                     DQ0 = DQ
                     IF ( ABS(ESTEP).GT.ESTEPTOL ) GOTO 110
                     ETRY = ETRY + ESTEP
                     GOTO 120
                  END IF
               END DO
               CALL STOP_MESSAGE(ROUTINE,'no Fermi energy found')
C
 120           CONTINUE
               EFERMILD = DREAL(ETRY)
               EFERMIDOS = EFERMI
               EFERMI = EFERMILD
Cc          efermi = efermidos
C
               NOS0 = NINT(DIMAG(NOSBOT))
               IF ( ABS(NOS0-DIMAG(NOSBOT)).GT.NOSTOL .AND. NCPA.NE.0 )
     &              WRITE (6,99005) NOSTOL,NOSBOT
C
               IF ( ABS(EFERMIDOS-EFERMILD).GT.DEFTOL ) WRITE (6,99006)
     &              DEFTOL
C
            END IF
C
C=======================================================================
C                                                      set up new E-path
            EMAX = EFERMI
C
            CALL EPATH(5,EMIN,EMAX,EIMAG,NETAB(1),ETAB(1,1),WETAB(1,1),
     &                 EILOW,IPRINT,NEMAX)
C
            IF ( .NOT.(LLOYD) ) THEN
               WRITE (6,99019) EFERMI,DOSFIL(1:(LDATSET+13))
            ELSE IF ( USE_EFERMILD ) THEN
               WRITE (6,99021) EFERMI,EFERMIDOS,DOSFIL(1:(LDATSET+13))
            ELSE
               WRITE (6,99022) EFERMI,DOSFIL(1:(LDATSET+13))
            END IF
C
         END IF
C ------------------------------------------ .NOT. EFERMI_GUESS_SUPPLIED
C ======================================================================
C
C ======================================================================
C    write out the starting potential once more - EFERMI is fixed now
C ======================================================================
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
            SCFSTATUS = 'ITR-BULK  '
         ELSE IF ( SYSTEM_TYPE(1:3).EQ.'LIV' .OR. SYSTEM_TYPE(1:3)
     &             .EQ.'LIR' ) THEN
            SCFSTATUS = 'ITR-L-BULK'
         ELSE IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
            SCFSTATUS = 'ITR-I-ZONE'
         END IF
C
         CALL POTWRSPR
C
C ======================================================================
      END IF
C
C       write DOS - improve EFERMI via Lloyd - new EPATH     MPI_ID.EQ.0
C=======================================================================
C
      CLOSE (10)
C
      DEALLOCATE (BCOR,CPACHTAB,MSSQ,MSST,TSST,JRAT,VR2AT)
      DEALLOCATE (IECPAFAIL,BCORS,SUMAT,WA,ECORTAB)
      DEALLOCATE (QIONINP,R1AT,SEBT,SEVT,RWK1,RWK2,CTOTNOS)
      DEALLOCATE (CTOTDOS,PHASK,ETAB0,WETAB0,PHASA,PHAST,DNSINI,TDOSINI)
      IF ( ALLOCATED(PCOEF) ) DEALLOCATE (PCOEF)
C
C ======================================================================
99001 FORMAT (/,5X,'INFO: SOLVER temporarily set to: ',A,/)
99002 FORMAT (/,5X,'INFO: SOLVER reset to: ',A,/)
99003 FORMAT (2(/,1X,79('#')),/,28X,'WARNING from <SCFINITPOT>',//,10X,
     &        'USEVMATT  was set in input file for  2D-system ',/,10X,
     &        'this option is not available --- the starting potential '
     &        ,/,10X,'will be derived from RHO (Mattheiss)',//,
     &        2(/,1X,79('#')),/)
99004 FORMAT (I5,8F12.6)
99005 FORMAT (10X,'Im(NOSBOT) deviates from next integer by more than',
     &        F10.6,/,10X,'NOSBOT = ',2F14.10)
99006 FORMAT (10X,'Fermi energies from DOS and Lloyd formula ',/,10X,
     &        'deviate by more than',F10.6,/)
99007 FORMAT (//,1X,79('*'),/,34X,'<SCFINITPOT>'/,1X,79('*'),//,10X,
     &        'SS-DOS written to file  ',A,/,10X,
     &        'for the energy range  E = ',F6.3,' ... ',F6.3,' Ry ',/)
99008 FORMAT (/,10X,'searching the Fermi energy ',/)
99009 FORMAT (I14,4X,'E',2F10.4,4X,'NOS',F10.4,4X,'DOS',F10.4)
99010 FORMAT (8E10.4,:,/,(10X,7E10.4))
99011 FORMAT (/,1X,79('*'),/,10X,'start potential set up using',
     &        ' Mattheiss construction for  V',/,1X,79('*'),/)
99012 FORMAT (/,1X,79('*'),/,10X,'start potential set up using',
     &        ' Mattheiss construction for  RHO',/,20X,
     &        'the corresponding ionicity has been used',/,1X,79('*'),/)
99013 FORMAT (/,1X,79('*'),/,10X,'start potential set up using',
     &        ' Mattheiss construction for  RHO',/,1X,79('*'),/)
99014 FORMAT (A10,A,A)
99015 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99016 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99017 FORMAT (1X,79('*'),/)
99018 FORMAT (//,1X,79('*'),/,10X,'from <SCFINITPOT>  ',/,10X,
     &        'no Fermi-energy could be found ',/,1X,79('*'),/)
99019 FORMAT (//,1X,79('*'),/,10X,'Fermi energy fixed by <SCFINITPOT>',
     &        /,10X,'EFERMI = ',F10.6,' Ry    using DOS',/,10X,
     &        'DOS written to file  ',A,/,1X,79('*'),/)
99020 FORMAT (//,1X,79('*'),/,10X,'Fermi energy fixed by <SCFINITPOT>',
     &        /,10X,'EFERMI = ',F10.6,' Ry    using guess from input',/,
     &        1X,79('*'),/)
99021 FORMAT (//,1X,79('*'),/,10X,'Fermi energy fixed by <SCFINITPOT>',
     &        /,10X,'EFERMI = ',F10.6,' Ry    using Lloyd''s formula',/,
     &        10X,'EFERMI = ',F10.6,' Ry    using DOS (for comparison)',
     &        /,10X,'DOS written to file  ',A,/,1X,79('*'),/)
99022 FORMAT (//,1X,79('*'),/,10X,'Fermi energy fixed by <SCFINITPOT>',
     &        /,10X,'EFERMI = ',F10.6,' Ry    using Lloyd''s formula',/,
     &        10X,'DOS written to file  ',A,/,1X,79('*'),/)
99023 FORMAT (//,1X,79('*'),/,34X,'<SCFINITPOT>',/,11X,
     &     'setting up a starting potential and fixing the Fermi energy'
     &     ,/,1X,79('*'),/)
      END
C*==scfinitrho.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCFINITRHO(USEIONMATT,USEQION,IPRINT,QEL,QION,QIONINP,
     &                      MSPIN,R1AT,VR2AT,QIONSCL)
C   ********************************************************************
C   *                                                                  *
C   *   prepare the SCF-cycle starting from  scratch                   *
C   *                                                                  *
C   *   - set up a guess for the charge and spin density               *
C   *     RHOCHR  and  RHOSPN, respectively, using atomic data         *
C   *     and the Mattheiss construction                               *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TYPES,ONLY:NT,NTMAX,TXT_T,Z,CONC,NAT,IMT,NT_L,RHOCHR,
     &    RHOSPN
      USE MOD_RMESH,ONLY:R,DRDI,R2DRDI,RWS,JRWS,NRMAX
      USE MOD_CALCMODE,ONLY:MOL
      USE MOD_LATTICE,ONLY:ALAT,ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,
     &    ADAINV_R,ADAINV_I,SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_SITES,ONLY:NQ,NQMAX,QBAS,NQ_L,NQ_R,IQAT,NOQ,ITOQ,NQ_I
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFINITRHO')
      REAL*8 DPAS
      PARAMETER (DPAS=0.05D0)
      INTEGER NRAT
      PARAMETER (NRAT=251)
C
C Dummy arguments
C
      INTEGER IPRINT
      REAL*8 QIONSCL
      LOGICAL USEIONMATT,USEQION
      REAL*8 MSPIN(NTMAX),QEL(NTMAX),QION(NTMAX),QIONINP(NTMAX),
     &       R1AT(NTMAX),VR2AT(NRAT,NTMAX)
C
C Local variables
C
      REAL*8 AUX,DELQ,NQCORR,QEL0(:),QES,QION0,QMAG,RHOAT(:,:),
     &       RHOSAT(:,:),RHOVAT(:,:),RINT(:),SUMAT(:,:),SUMMSPIN,SUMQEL,
     &       SUMQION,SUMQIONINP,SUMQIONMC,SUMZ,WA(:),WB(:),WC(:),XAUX(2)
     &       ,YAUX(2),ZZ(:)
      LOGICAL ESPRESENT
      INTEGER IA_ERR,IM,IO,IQ,IR,IRTOP,IT,JRAT(:)
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RINT,JRAT,ZZ,SUMAT,QEL0
      ALLOCATABLE WA,WB,WC,RHOAT,RHOSAT,RHOVAT
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (JRAT(NTMAX),QEL0(NTMAX),RINT(NRMAX))
      ALLOCATE (RHOVAT(NRAT,NTMAX),WA(NRMAX),ZZ(NTMAX))
      ALLOCATE (RHOAT(NRAT,NTMAX),WC(NRMAX),WB(NRMAX))
      ALLOCATE (RHOSAT(NRAT,NTMAX),SUMAT(NRAT,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RHOVAT')
C
      CALL RINIT(NRAT*NTMAX,RHOAT)
      CALL RINIT(NRAT*NTMAX,RHOSAT)
      CALL RINIT(NRAT*NTMAX,RHOVAT)
C
      ESPRESENT = .FALSE.
C
C ======================================================================
C               calculate the starting charge and spin density
C ======================================================================
C
      DO IT = 1,NT
         ZZ(IT) = DBLE(Z(IT))
      END DO
C
      CALL SCF0ATOM(IPRINT,NT,ZZ,R1AT,TXT_T,RHOAT,RHOSAT,RHOVAT,VR2AT,
     &              NRAT)
C
C=======================================================================
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         WRITE (6,99001) ' '
      ELSE
         WRITE (6,99001) 'for I-regime'
      END IF
C=======================================================================
C
      CALL SCF0MATT(MOL,NQ,NQ_L,NQ_R,NQMAX,QBAS,ABAS,ABAS_L,ABAS_R,
     &              ABAS_I,ADAINV_L,ADAINV_R,ADAINV_I,SYSTEM_DIMENSION,
     &              ALAT,DPAS,SUMAT,RHOAT,R1AT,JRAT,RWS,CONC,NOQ,IMT,
     &              IQAT,ITOQ,NT,NTMAX,NRAT,WA,NRMAX)
C
      SUMQION = 0D0
      SUMQIONINP = 0D0
      SUMQIONMC = 0D0
      SUMMSPIN = 0D0
      QES = 0D0
C
      IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
C
         DO IT = 1,NT_L
            SUMAT(1:NRAT,IT) = 0.0D0
            RHOSAT(1:NRAT,IT) = 0.0D0
         END DO
C
         DO IQ = NQ_L + NQ_I + 1,NQ
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               SUMAT(1:NRAT,IT) = 0.0D0
               RHOSAT(1:NRAT,IT) = 0.0D0
            END DO
         END DO
C
      END IF
C
      NQCORR = 0.0D0
C
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
C---------------------------------------------------- set up atomic mesh
C
         DO IR = 1,JRAT(IT)
            WA(IR) = R1AT(IT)*EXP(DBLE(IR-1)*DPAS)
         END DO
C
C---------------------------------------------- renormalize spin density
C
         DO IR = 1,IRTOP
            WB(IR) = YLAG(R(IR,IM),WA,RHOSAT(1,IT),0,3,JRAT(IT))
     &               *DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,AUX)
C
         IF ( AUX.LE.0.0D0 ) THEN
            AUX = 1D0
         ELSE
            AUX = MSPIN(IT)/AUX
         END IF
C
         DO IR = 1,IRTOP
            WB(IR) = AUX*WB(IR)
            RHOSPN(IR,IT) = WB(IR)/R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,QMAG)
C
C-----------------------------------------------------------------------
C                use RHOVAT to set up magnetisation density
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG).LT.1D-6 .AND. ABS(MSPIN(IT)).GT.1D-6 ) THEN
C
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
            DO IR = 1,IRTOP
               IF ( R(IR,IM).LT.WA(1) ) THEN
                  XAUX(1) = 0.0D0
                  XAUX(2) = WA(1)
                  YAUX(1) = 0.0D0
                  YAUX(2) = RHOVAT(1,IT)
                  WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
                  XAUX(1) = WA(JRAT(IT))
                  XAUX(2) = 5.0D0*R(IRTOP,IM)
                  YAUX(1) = RHOVAT(JRAT(IT),IT)
                  YAUX(2) = 0.0D0
                  WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE
                  WB(IR) = YLAG(R(IR,IM),WA,RHOVAT(1,IT),0,3,JRAT(IT))
               END IF
               WB(IR) = WB(IR)*DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,AUX)
C
            IF ( AUX.LE.0.0D0 ) THEN
               AUX = 1D0
            ELSE
               AUX = MSPIN(IT)/AUX
            END IF
C
            DO IR = 1,IRTOP
               WB(IR) = AUX*WB(IR)
               RHOSPN(IR,IT) = WB(IR)/R2DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,QMAG)
C
         END IF
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG-MSPIN(IT)).GT.1D-6 )
     &         CALL STOP_MESSAGE(ROUTINE,'QMAG <> MSPIN')
C
C---------------------------------------------- deal with charge density
C
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
         DO IR = 1,IRTOP
            IF ( R(IR,IM).LT.WA(1) ) THEN
               XAUX(1) = 0.0D0
               XAUX(2) = WA(1)
               YAUX(1) = 0.0D0
               YAUX(2) = SUMAT(1,IT)
               WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
            ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
               XAUX(1) = WA(JRAT(IT))
               XAUX(2) = 5.0D0*R(IRTOP,IM)
               YAUX(1) = SUMAT(JRAT(IT),IT)
               YAUX(2) = 0.0D0
               WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
            ELSE
               WB(IR) = YLAG(R(IR,IM),WA,SUMAT(1,IT),0,3,JRAT(IT))
            END IF
C
            WB(IR) = WB(IR)*DRDI(IR,IM)
            RHOCHR(IR,IT) = WB(IR)/R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,QEL0(IT))
C
         IF ( Z(IT).LE.0 ) THEN
            QES = QES + QEL0(IT)*CONC(IT)*NAT(IT)
            ESPRESENT = .TRUE.
         END IF
C
         IF ( USEIONMATT ) THEN
            QION(IT) = Z(IT) - QEL0(IT)
         ELSE
            QION(IT) = QIONINP(IT)
         END IF
C
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
         SUMQIONINP = SUMQIONINP + QIONINP(IT)*CONC(IT)*NAT(IT)
         SUMQIONMC = SUMQIONMC + (Z(IT)-QEL0(IT))*CONC(IT)*NAT(IT)
         SUMMSPIN = SUMMSPIN + MSPIN(IT)*CONC(IT)*NAT(IT)
         IF ( (Z(IT).GT.0 .AND. CONC(IT).GT.0.001D0) .OR. 
     &        (Z(IT).EQ.0 .AND. CONC(IT).LT.0.999D0) ) NQCORR = NQCORR + 
     &        CONC(IT)*NAT(IT)
      END DO
      NQCORR = NINT(NQCORR)
C
      DELQ = SUMQION/NQCORR
      SUMQION = 0.0D0
      DO IT = 1,NT
         IF ( (Z(IT).GT.0 .AND. CONC(IT).GT.0.001D0) .OR. 
     &        (Z(IT).EQ.0 .AND. CONC(IT).LT.0.999D0) ) QION(IT)
     &        = QION(IT) - DELQ
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
      END DO
C
      QION(1:NT) = QION(1:NT)*QIONSCL
C
      IF ( USEQION ) THEN
         WRITE (6,99003) ' Q (inp)     '
      ELSE
         WRITE (6,99003) ' '
      END IF
      DO IT = 1,NT
         IF ( USEQION ) THEN
            WRITE (6,99004) IT,TXT_T(IT),QIONINP(IT),Z(IT) - QEL0(IT),
     &                      QION(IT),MSPIN(IT)
         ELSE
            WRITE (6,99004) IT,TXT_T(IT),Z(IT) - QEL0(IT),QION(IT),
     &                      MSPIN(IT)
         END IF
      END DO
      IF ( USEQION ) THEN
         WRITE (6,99005) SUMQIONINP,SUMQIONMC,SUMQION,SUMMSPIN
         IF ( ESPRESENT ) WRITE (6,99002) '             ',QES
      ELSE
         WRITE (6,99005) SUMQIONMC,SUMQION,SUMMSPIN
         WRITE (6,99007) QIONSCL
         IF ( ESPRESENT ) WRITE (6,99002) ' ',QES
      END IF
C
C-----------  Renormalization of the charge density according to SUMQION
C
      SUMQEL = 0D0
      SUMZ = 0D0
C
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
C---------------------------------------------------- set up atomic mesh
C
         DO IR = 1,JRAT(IT)
            WA(IR) = R1AT(IT)*EXP(DBLE(IR-1)*DPAS)
         END DO
C
         QION0 = Z(IT) - QEL0(IT)
C
         IF ( Z(IT).GT.0 ) THEN
C----------- normalize valence orbital density and add to charge density
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
            DO IR = 1,IRTOP
               IF ( R(IR,IM).LT.WA(1) ) THEN
                  XAUX(1) = 0.0D0
                  XAUX(2) = WA(1)
                  YAUX(1) = 0.0D0
                  YAUX(2) = RHOVAT(1,IT)
                  AUX = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
                  XAUX(1) = WA(JRAT(IT))
                  XAUX(2) = 5.0D0*R(IRTOP,IM)
                  YAUX(1) = RHOVAT(JRAT(IT),IT)
                  YAUX(2) = 0.0D0
                  AUX = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE
                  AUX = YLAG(R(IR,IM),WA,RHOVAT(1,IT),0,3,JRAT(IT))
               END IF
C
               WC(IR) = AUX/R(IR,IM)**2
               WB(IR) = AUX*DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,AUX)
C
            AUX = (QION(IT)-QION0)/AUX
            RHOCHR(1:IRTOP,IT) = RHOCHR(1:IRTOP,IT) - AUX*WC(1:IRTOP)
C----------------------------- normalize charge density for empty sphere
         ELSE IF ( ABS(QION0).GT.1D-10 ) THEN
            AUX = QION(IT)/QION0
            RHOCHR(1:IRTOP,IT) = AUX*RHOCHR(1:IRTOP,IT)
         ELSE
            RHOCHR(1:IRTOP,IT) = 0.0D0
         END IF
C
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,RINT,QEL(IT))
C
         SUMQEL = SUMQEL + QEL(IT)*CONC(IT)*NAT(IT)
         SUMZ = SUMZ + Z(IT)*CONC(IT)*NAT(IT)
C
      END DO
C
      IF ( ABS(SUMQEL-SUMZ).GT.1D-5 ) THEN
         WRITE (6,99006) SUMQEL,SUMZ
         CALL STOP_MESSAGE(ROUTINE,
     &       'no charge neutrality    setting up a guess charge density'
     &       )
      END IF
C
C=======================================================================
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) RETURN
C=======================================================================
C
C***********************************************************************
C deal now with left side and overwrite charge density for IT=1,..,NT_L
C***********************************************************************
C
C=======================================================================
C             left side is vacuum, i.e. slab calculation
C=======================================================================
      IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
C
         DO IT = 1,NT_L
            IM = IMT(IT)
            IRTOP = JRWS(IM)
            QION(IT) = 0D0
            QEL(IT) = 0D0
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = 0D0
               RHOSPN(IR,IT) = 0D0
            END DO
         END DO
C
         DO IQ = NQ_L + NQ_I + 1,NQ
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               IM = IMT(IT)
               IRTOP = JRWS(IM)
               QION(IT) = 0D0
               QEL(IT) = 0D0
               DO IR = 1,IRTOP
                  RHOCHR(IR,IT) = 0D0
                  RHOSPN(IR,IT) = 0D0
               END DO
            END DO
         END DO
C
         RETURN
C
      END IF
C=======================================================================
C
      WRITE (6,99001) 'for L-regime'
C
C
      CALL SCF0MATT(MOL,NQ_L,NQ_L,NQ_R,NQMAX,QBAS,ABAS_L,ABAS_L,ABAS_R,
     &              ABAS_I,ADAINV_L,ADAINV_R,ADAINV_I,'3D        ',ALAT,
     &              DPAS,SUMAT,RHOAT,R1AT,JRAT,RWS,CONC,NOQ,IMT,IQAT,
     &              ITOQ,NT_L,NTMAX,NRAT,WA,NRMAX)
C
      SUMQION = 0D0
      SUMQIONINP = 0D0
      SUMQIONMC = 0D0
      SUMMSPIN = 0D0
      QES = 0D0
C
      DO IT = 1,NT_L
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
C---------------------------------------------------- set up atomic mesh
C
         DO IR = 1,JRAT(IT)
            WA(IR) = R1AT(IT)*EXP(DBLE(IR-1)*DPAS)
         END DO
C
C---------------------------------------------- renormalize spin density
C
         DO IR = 1,IRTOP
            WB(IR) = YLAG(R(IR,IM),WA,RHOSAT(1,IT),0,3,JRAT(IT))
     &               *DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,AUX)
C
         IF ( AUX.LE.0.0D0 ) THEN
            AUX = 1D0
         ELSE
            AUX = MSPIN(IT)/AUX
         END IF
C
         DO IR = 1,IRTOP
            WB(IR) = AUX*WB(IR)
            RHOSPN(IR,IT) = WB(IR)/R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,QMAG)
C
C-----------------------------------------------------------------------
C                use RHOVAT to set up magnetisation density
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG).LT.1D-6 .AND. ABS(MSPIN(IT)).GT.1D-6 ) THEN
C
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
            DO IR = 1,IRTOP
               IF ( R(IR,IM).LT.WA(1) ) THEN
                  XAUX(1) = 0.0D0
                  XAUX(2) = WA(1)
                  YAUX(1) = 0.0D0
                  YAUX(2) = RHOVAT(1,IT)
                  WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
                  XAUX(1) = WA(JRAT(IT))
                  XAUX(2) = 5.0D0*R(IRTOP,IM)
                  YAUX(1) = RHOVAT(JRAT(IT),IT)
                  YAUX(2) = 0.0D0
                  WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE
                  WB(IR) = YLAG(R(IR,IM),WA,RHOVAT(1,IT),0,3,JRAT(IT))
               END IF
               WB(IR) = WB(IR)*DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,AUX)
C
            IF ( AUX.LE.0.0D0 ) THEN
               AUX = 1D0
            ELSE
               AUX = MSPIN(IT)/AUX
            END IF
C
            DO IR = 1,IRTOP
               WB(IR) = AUX*WB(IR)
               RHOSPN(IR,IT) = WB(IR)/R2DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,QMAG)
C
         END IF
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG-MSPIN(IT)).GT.1D-6 )
     &         CALL STOP_MESSAGE(ROUTINE,'QMAG <> MSPIN')
C
C---------------------------------------------- deal with charge density
C
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
         DO IR = 1,IRTOP
            IF ( R(IR,IM).LT.WA(1) ) THEN
               XAUX(1) = 0.0D0
               XAUX(2) = WA(1)
               YAUX(1) = 0.0D0
               YAUX(2) = SUMAT(1,IT)
               WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
            ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
               XAUX(1) = WA(JRAT(IT))
               XAUX(2) = 5.0D0*R(IRTOP,IM)
               YAUX(1) = SUMAT(JRAT(IT),IT)
               YAUX(2) = 0.0D0
               WB(IR) = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
            ELSE
               WB(IR) = YLAG(R(IR,IM),WA,SUMAT(1,IT),0,3,JRAT(IT))
            END IF
C
            WB(IR) = WB(IR)*DRDI(IR,IM)
            RHOCHR(IR,IT) = WB(IR)/R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,WB,QEL0(IT))
C
         IF ( Z(IT).LE.0 ) THEN
            QES = QES + QEL0(IT)*CONC(IT)*NAT(IT)
            ESPRESENT = .TRUE.
         END IF
C
         IF ( USEIONMATT ) THEN
            QION(IT) = Z(IT) - QEL0(IT)
         ELSE
            QION(IT) = QIONINP(IT)
         END IF
C
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
         SUMQIONINP = SUMQIONINP + QIONINP(IT)*CONC(IT)*NAT(IT)
         SUMQIONMC = SUMQIONMC + (Z(IT)-QEL0(IT))*CONC(IT)*NAT(IT)
         SUMMSPIN = SUMMSPIN + MSPIN(IT)*CONC(IT)*NAT(IT)
      END DO
C
      DELQ = SUMQION/DBLE(NQ_L)
      SUMQION = 0.0D0
      DO IT = 1,NT_L
         QION(IT) = QION(IT) - DELQ
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
      END DO
C
      IF ( USEQION ) THEN
         WRITE (6,99003) ' Q (inp)     '
      ELSE
         WRITE (6,99003) ' '
      END IF
      DO IT = 1,NT_L
         IF ( USEQION ) THEN
            WRITE (6,99004) IT,TXT_T(IT),QIONINP(IT),Z(IT) - QEL0(IT),
     &                      QION(IT),MSPIN(IT)
         ELSE
            WRITE (6,99004) IT,TXT_T(IT),Z(IT) - QEL0(IT),QION(IT),
     &                      MSPIN(IT)
         END IF
      END DO
      IF ( USEQION ) THEN
         WRITE (6,99005) SUMQIONINP,SUMQIONMC,SUMQION,SUMMSPIN
         IF ( ESPRESENT ) WRITE (6,99002) '             ',QES
      ELSE
         WRITE (6,99005) SUMQIONMC,SUMQION,SUMMSPIN
         IF ( ESPRESENT ) WRITE (6,99002) ' ',QES
      END IF
C
C-----------  Renormalization of the charge density according to SUMQION
C
      SUMQEL = 0D0
      SUMZ = 0D0
C
      DO IT = 1,NT_L
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
C---------------------------------------------------- set up atomic mesh
C
         DO IR = 1,JRAT(IT)
            WA(IR) = R1AT(IT)*EXP(DBLE(IR-1)*DPAS)
         END DO
C
         QION0 = Z(IT) - QEL0(IT)
C
         IF ( Z(IT).GT.0 ) THEN
C----------- normalize valence orbital density and add to charge density
C
C OS: Make sure that density is not negative due to extrapolation
C    Linear interpolation between zero and first atomic data point
C    and between last real point and zero density far away.
C
            DO IR = 1,IRTOP
               IF ( R(IR,IM).LT.WA(1) ) THEN
                  XAUX(1) = 0.0D0
                  XAUX(2) = WA(1)
                  YAUX(1) = 0.0D0
                  YAUX(2) = RHOVAT(1,IT)
                  AUX = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE IF ( R(IR,IM).GT.WA(JRAT(IT)) ) THEN
                  XAUX(1) = WA(JRAT(IT))
                  XAUX(2) = 5.0D0*R(IRTOP,IM)
                  YAUX(1) = RHOVAT(JRAT(IT),IT)
                  YAUX(2) = 0.0D0
                  AUX = YLAG(R(IR,IM),XAUX,YAUX,0,2,2)
               ELSE
                  AUX = YLAG(R(IR,IM),WA,RHOVAT(1,IT),0,3,JRAT(IT))
               END IF
C
               WC(IR) = AUX/R(IR,IM)**2
               WB(IR) = AUX*DRDI(IR,IM)
            END DO
C
            CALL RRADINT(IM,WB,AUX)
C
            AUX = (QION(IT)-QION0)/AUX
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = RHOCHR(IR,IT) - AUX*WC(IR)
            END DO
         ELSE
C----------------------------- normalize charge density for empty sphere
            AUX = QION(IT)/QION0
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = AUX*RHOCHR(IR,IT)
            END DO
         END IF
C
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,RINT,QEL(IT))
C
         SUMQEL = SUMQEL + QEL(IT)*CONC(IT)*NAT(IT)
         SUMZ = SUMZ + Z(IT)*CONC(IT)*NAT(IT)
C
      END DO
C
      IF ( ABS(SUMQEL-SUMZ).GT.1D-5 ) THEN
         WRITE (6,99006) SUMQEL,SUMZ
         CALL STOP_MESSAGE(ROUTINE,
     &          'no charge neutrality setting up a guess charge density'
     &          )
      END IF
C
C=======================================================================
C
      DEALLOCATE (RINT,JRAT,RHOAT,SUMAT)
      DEALLOCATE (WA,WB,WC,RHOSAT,ZZ,RHOVAT)
      DEALLOCATE (QEL0)
C
C=======================================================================
99001 FORMAT (/,1X,79('*'),/,34X,'<SCFINITRHO>'/,1X,79('*'),//,10X,
     &        'initializing charge density ',A,/)
99002 FORMAT (/,5X,'charge in empty spheres QES',A,F10.4)
99003 FORMAT (5X,'charge and spin moment from input and Mattheiss ',
     &        'construction',//,23X,A,4X,
     &        ' Z-Q_el     Q(corr)      m_spin')
99004 FORMAT (5X,'type',I3,' ',A,' :',4F12.4)
99005 FORMAT (5X,'sum ',12x,' :',4F12.4)
99006 FORMAT (/,5X,'number of electrons    ',F12.6,/,5X,
     &        'weighted sum over  Z   ',F12.6,/)
99007 FORMAT (/,5X,'QION will be scaled by     ',F10.4)
      END
C*==scf0matt.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0MATT(MOL,NQ,NQ_L,NQ_R,NQMAX,QBAS,ABAS,ABAS_L,
     &                    ABAS_R,ABAS_I,ADAINV_L,ADAINV_R,ADAINV_I,
     &                    SYSTEM_DIMENSION,ALAT,DPAS,SUMYAT,YAT,R1AT,
     &                    JRAT,RWS,CONC,NOQ,IMT,IQAT,ITOQ,NT,NTMAX,NRAT,
     &                    WA,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   use the Mattheis prescription to sum up contributions          *
C   *   of neighboring sites for a central site                        *
C   *   summation may be made for   charge   or   potential            *
C   *   depending on supplied input  YAT                               *
C   *   store results in variable    SUMYAT                            *
C   *                                                                  *
C   *   MOL = .FALSE.  -->  solid state calculation                    *
C   *                       contributions from neighboring unit cells  *
C   *         .TRUE.   -->  cluster calculation                        *
C   *                       NO neighbor contributions                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,IOTMP
      USE MOD_LATTICE,ONLY:SWS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0MATT')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
      INTEGER NTR,NTRM
      PARAMETER (NTR=4,NTRM=(2*NTR+1)*(2*NTR+1)*(2*NTR+1))
C
C Dummy arguments
C
      REAL*8 ALAT,DPAS
      LOGICAL MOL
      INTEGER NQ,NQMAX,NQ_L,NQ_R,NRAT,NRMAX,NT,NTMAX
      CHARACTER*10 SYSTEM_DIMENSION
      REAL*8 ABAS(3,3),ABAS_I(3,3),ABAS_L(3,3),ABAS_R(3,3),ADAINV_I(3,3)
     &       ,ADAINV_L(3,3),ADAINV_R(3,3),CONC(NTMAX),QBAS(3,NQMAX),
     &       R1AT(*),RWS(*),SUMYAT(NRAT,*),WA(NRMAX),YAT(NRAT,NTMAX)
      INTEGER IMT(NTMAX),IQAT(NQMAX,*),ITOQ(NTMAX,NQMAX),JRAT(NTMAX),
     &        NOQ(NQMAX)
C
C Local variables
C
      REAL*8 BBRX(:),BBRY(:),BBRZ(:),DIST2,DISTMAX2,DQMATT(:),MATTRAD,
     &       RAD,RQMATT(:,:)
      LOGICAL DONE,KDONET(NT)
      INTEGER I,IA_ERR,IFIL,IJKLIM,IO,IPRINT_LOW,IQ,IQCNTR,IQMATT,IQP,
     &        IQ_IQMATT(:),IT,ITP,IVEC,J,JO,JQ,JT,K,NQMATT,NQMATT_I,
     &        NQMATT_L,NQMATT_R,NSHLMATT,NVEC
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BBRX,BBRY,BBRZ,DQMATT,IQ_IQMATT,RQMATT
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( NTRM.LE.0 ) STOP
      ALLOCATE (BBRX(NTRM),BBRY(NTRM),BBRZ(NTRM))
C
      IF ( IPRINT.LE.0 ) THEN
         IPRINT_LOW = -1
      ELSE
         IPRINT_LOW = IPRINT
      END IF
C
      IF ( IPRINT_LOW.GE.0 ) WRITE (6,99001) SWS/ALAT,SWS
      DISTMAX2 = (5.0D0*SWS)**2
C
      IF ( MOL ) THEN
         IJKLIM = 0
      ELSE
         IJKLIM = NTR
      END IF
C
      IVEC = 0
      DO I = -IJKLIM,IJKLIM
         DO J = -IJKLIM,IJKLIM
            DO K = -IJKLIM,IJKLIM
               IVEC = IVEC + 1
               BBRX(IVEC) = I*ABAS(1,1) + J*ABAS(1,2) + K*ABAS(1,3)
               BBRY(IVEC) = I*ABAS(2,1) + J*ABAS(2,2) + K*ABAS(2,3)
               BBRZ(IVEC) = I*ABAS(3,1) + J*ABAS(3,2) + K*ABAS(3,3)
            END DO
         END DO
      END DO
      NVEC = IVEC
C
C-----------------------------------------------------------------------
C                     add on-site contribution
C-----------------------------------------------------------------------
      DO IT = 1,NT
         JRAT(IT) = MIN(NRAT,NINT(LOG(RWS(IMT(IT))/R1AT(IT))/DPAS)+2)
         DO I = 1,JRAT(IT)
            SUMYAT(I,IT) = YAT(I,IT)
         END DO
         KDONET(IT) = .FALSE.
      END DO
C
C=======================================================================
C
      IF ( NQ.GT.00 ) THEN
         DO IQ = 1,NQ
            DONE = .TRUE.
            DO IO = 1,NOQ(IQ)
               DONE = DONE .AND. KDONET(ITOQ(IO,IQ))
            END DO
C
            IF ( .NOT.DONE ) THEN
C
C-----------------------------------------------------------------------
C                 generate cluster around site IQ
C-----------------------------------------------------------------------
               IQCNTR = IQ
               NSHLMATT = 0
               MATTRAD = 5.0D0*(SWS/ALAT)
C
               CALL CLUSSITES(IOTMP,IPRINT_LOW,MOL,SYSTEM_DIMENSION,
     &                        ABAS,ABAS_L,ABAS_I,ABAS_R,ADAINV_L,
     &                        ADAINV_I,ADAINV_R,QBAS,MATTRAD,IQCNTR,
     &                        NQMATT,NQMATT_L,NQMATT_I,NQMATT_R,
     &                        NSHLMATT,NQ,NQ_L,NQ_R,NQMAX)
C
               IF ( IQ.GT.1 ) DEALLOCATE (RQMATT,DQMATT,IQ_IQMATT)
               ALLOCATE (RQMATT(3,NQMATT),DQMATT(NQMATT))
               ALLOCATE (IQ_IQMATT(NQMATT),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQMATTS')
C
               READ (IOTMP) ((RQMATT(J,I),J=1,3),DQMATT(I),IQ_IQMATT(I),
     &                      I=1,NQMATT)
               CLOSE (IOTMP)
               IF ( IPRINT.GE.5 ) WRITE (6,'(I5,4f10.6,I5)')
     &              (I,(RQMATT(J,I),J=1,3),DQMATT(I),IQ_IQMATT(I),I=1,
     &              NQMATT)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
                  KDONET(IT) = .TRUE.
C
C-----------------------------------------------------------------------
C       add contributions from surrounding ---> IQMATT = 2,NQMATT
C-----------------------------------------------------------------------
C
                  DO IQMATT = 2,NQMATT
                     JQ = IQ_IQMATT(IQMATT)
                     RAD = DQMATT(IQMATT)*ALAT
C
                     DO JO = 1,NOQ(JQ)
                        JT = ITOQ(JO,JQ)
                        CALL SCF0SUMUP(RAD,R1AT(JT),R1AT(IT),DPAS,
     &                                 YAT(1,JT),SUMYAT(1,IT),JRAT(IT),
     &                                 NRAT,CONC(JT),WA,NRMAX)
                     END DO
                  END DO
C
               END DO
C
            END IF
C
         END DO
C
      END IF
C
C=======================================================================
      IF ( .NOT.CHECK ) RETURN
C=======================================================================
C
      IFIL = 100
      WRITE (IFIL,'(f17.12)') ((SUMYAT(I,IT),I=1,JRAT(IT)),IT=1,NT)
C
C=======================================================================
C
      DO IT = 1,NT
         DO I = 1,JRAT(IT)
            SUMYAT(I,IT) = YAT(I,IT)
         END DO
      END DO
C
C--------------------------------------------------  loop over the types
      DO IT = 1,NT
         IQ = IQAT(1,IT)
C
         DO IVEC = 1,NVEC
            DO IQP = 1,NQ
               DIST2 = (QBAS(1,IQ)-QBAS(1,IQP)+BBRX(IVEC))
     &                 **2 + (QBAS(2,IQ)-QBAS(2,IQP)+BBRY(IVEC))
     &                 **2 + (QBAS(3,IQ)-QBAS(3,IQP)+BBRZ(IVEC))**2
C
               IF ( (DIST2.LT.DISTMAX2) .AND. (DIST2.GT.0.001D0) ) THEN
C------------------------------------------------ distance between sites
                  RAD = SQRT(DIST2)*ALAT
C
                  DO IO = 1,NOQ(IQP)
                     ITP = ITOQ(IO,IQP)
C
                     CALL SCF0SUMUP(RAD,R1AT(ITP),R1AT(IT),DPAS,
     &                              YAT(1,ITP),SUMYAT(1,IT),JRAT(IT),
     &                              NRAT,CONC(ITP),WA,NRMAX)
                  END DO
               END IF
            END DO
         END DO
C
      END DO
C
      IFIL = 200
      WRITE (IFIL,'(f17.12)') ((SUMYAT(I,IT),I=1,JRAT(IT)),IT=1,NT)
      DEALLOCATE (BBRX,BBRY,BBRZ)
C
      WRITE (6,99002)
      CALL STOP_MESSAGE(ROUTINE,'testing CHARGE SET UP')
C----------------------------------------------------------------------
99001 FORMAT (/,10X,'average Wigner-Seitz radius ',F12.6,' a  = ',F12.6,
     &        ' a.u.',/)
99002 FORMAT (/,2(/,1X,79('*')),/,10X,
     &        'test of charge setuo in <SCF0MATT> completed ',/,10X,
     &        'charge according to new scheme     in  fort.100  ',/,10X,
     &        'charge according to Sasha''s scheme in  fort.200  ',/,
     &        2(/,1X,79('*')),/)
      END
C*==scf0sumup.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0SUMUP(R,R1,R2,DPAS,RO1,RO2,JRAT,NRAT,AK,WA,NRMAX)
C **********************************************************************
C *                                                                    *
C *   calculation of the spherically averaged                          *
C *   charge density (potential) from atom at distance  R              *
C *   the    density (potential) has to be scaled by   4*PI*r^2  (r^2) *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 AK,DPAS,R,R1,R2
      INTEGER JRAT,NRAT,NRMAX
      REAL*8 RO1(NRAT),RO2(NRAT),WA(NRMAX)
C
C Local variables
C
      INTEGER J2,JMAX,JMIN,JP,JR
      REAL*8 REZ,RMAX,RMIN,RP,XJMAX,XJMIN
C
C*** End of declarations rewritten by SPAG
C
      WA(1:NRMAX) = 0D0
C
      IF ( R.LT.1D-2 ) RETURN
C
      RO1(NRAT) = 0
C
      DO JR = JRAT,1, - 1
         RP = R2*EXP(DPAS*(JR-1))
         RMIN = R - RP
         RMAX = R + RP
         XJMIN = LOG(RMIN/R1)/DPAS + 1.D0
         JMIN = INT(XJMIN)
         IF ( JMIN.GE.NRAT-2 ) RETURN
         XJMAX = LOG(RMAX/R1)/DPAS + 1.D0
         IF ( XJMAX.GE.NRAT ) XJMAX = NRAT - 0.001D0
         JMAX = INT(XJMAX) + 1
         J2 = JMAX - JMIN + 1
         JP = J2
         IF ( J2/2*2.EQ.J2 ) JP = J2 + 1
         IF ( JMIN+(JP-1).GT.NRAT ) THEN
            JP = 2*((NRAT-JMIN)/2) + 1
            IF ( JMIN+(JP-1).GT.NRAT ) JP = JP - 2
         END IF
C
         IF ( JP.GE.2 ) THEN
C
            CALL SCF0QD(DPAS,RO1(JMIN),WA,JP)
C
            REZ = -(XJMIN-JMIN)*WA(2) + (XJMAX-JMAX)*(WA(J2)-WA(J2-1))
     &            + WA(J2)
C
            RO2(JR) = RO2(JR) + AK*REZ*RP/(2*R)
         END IF
      END DO
      END
C*==scf0atom.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0ATOM(IPRINT,NT,Z,R1AT,TXT_T,RHOAT,RHOSAT,RHOVAT,
     &                    VR2AT,NRAT)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0ATOM')
C
C Dummy arguments
C
      INTEGER IPRINT,NRAT,NT
      REAL*8 R1AT(*),RHOAT(NRAT,*),RHOSAT(NRAT,*),RHOVAT(NRAT,*),
     &       VR2AT(NRAT,*),Z(*)
      CHARACTER*8 TXT_T(*)
C
C Local variables
C
      INTEGER I,ISR,IT,IZN
C
C*** End of declarations rewritten by SPAG
C
      DO IT = 1,NT
         IF ( Z(IT).LT.0 ) CALL STOP_MESSAGE(ROUTINE,' Z(IT).LT.0')
         IF ( Z(IT).LT.0.3D0 ) THEN
C empty sphere
            DO I = 1,NRAT
               VR2AT(I,IT) = 0.0D0
               RHOAT(I,IT) = 0.D0
               RHOSAT(I,IT) = 0.D0
               RHOVAT(I,IT) = 0.D0
            END DO
            R1AT(IT) = 1D-5
         ELSE
            DO ISR = 1,IT - 1
               IF ( NINT(Z(ISR)).EQ.NINT(Z(IT)) ) THEN
                  DO I = 1,NRAT
                     VR2AT(I,IT) = VR2AT(I,ISR)
                     RHOAT(I,IT) = RHOAT(I,ISR)
                     RHOSAT(I,IT) = RHOSAT(I,ISR)
                     RHOVAT(I,IT) = RHOVAT(I,ISR)
                  END DO
                  R1AT(IT) = R1AT(ISR)
                  GOTO 100
               END IF
            END DO
C
            IZN = NINT(Z(IT))
            WRITE (6,'(/,5X,''for atom '',A,''    Z ='',I4/)') TXT_T(IT)
     &             ,IZN
C
            CALL SCF0RHFDS(IPRINT,IZN,RHOAT(1,IT),RHOSAT(1,IT),
     &                     RHOVAT(1,IT),VR2AT(1,IT),R1AT(IT))
         END IF
 100  END DO
      END
C*==scf0rhfds.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0RHFDS(IPRINT,IZ,RHOAT,RHOSAT,RHOVAT,VR2AT,DR1AT)
C   ********************************************************************
C   *                                                                  *
C   *   HARTREE FOCK DIRAC SLATER   J P DESCLAUX     CEA PARIS 1969    *
C   *                                                                  *
C   ********************************************************************
C   03-03-94 - PERLOV : NDP cannot calculate    (0.0)**(1.D0/3.D0)
C                           changed to (0.0 + 1.d-42)**(1.D0/3.D0)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0RHFDS')
C
C Dummy arguments
C
      REAL*8 DR1AT
      INTEGER IPRINT,IZ
      REAL*8 RHOAT(251),RHOSAT(251),RHOVAT(251),VR2AT(251)
C
C Local variables
C
      CHARACTER*4 BAR(10),TITRE(30)
      REAL*8 D(251),DC(251),DCOP,DE,DEN(30),DEXE,DEXV,DFL(30),
     &       DGC(251,30),DP(251),DPAS,DPC(251,30),DQ(251),DQ1(30),
     &       DR(251),DV(251),DVAL,DVF(251),EMAX,TEST,TESTE,TESTV,TESTY,
     &       TETS,VAL,VMAX,Y,YMAX,YN,Z
      INTEGER I,ICEL,IM,IMAX,ION,ITER,J,JJ,JVAL,K,L,LVAL,MAG(30),N,
     &        NCORB(30),NEL(30),NES,NITER,NK(30),NMAX(30),NORB,NP,
     &        NQL(30),NQN(30),NSTOP,NUC
C
C*** End of declarations rewritten by SPAG
C
      DATA YN/0D0/,Y/0D0/
C
      NSTOP = 1
C
 100  CONTINUE
      CALL SCF0INSLD(IPRINT,IZ,DEN,DQ1,DFL,NQN,NQL,NK,NMAX,NEL,NORB,MAG,
     &               DV,DR,DPAS,Z,NSTOP,NES,NP,NUC,DEXV,DEXE,DCOP,TITRE,
     &               BAR,TEST,TESTE,TESTY,TESTV,NITER,ION,NCORB)
C
      ITER = 1
      DO I = 1,NP
         DO J = 1,NORB
            DGC(I,J) = 0.D0
            DPC(I,J) = 0.D0
         END DO
      END DO
      IF ( IPRINT.GT.2 ) THEN
         WRITE (6,99001) BAR
         WRITE (6,99002)
      END IF
      N = -(ION+1)
 200  CONTINUE
      DO I = 1,NP
         D(I) = 0.D0
      END DO
      TETS = TEST
      YMAX = 0.D0
      VMAX = 0.D0
      EMAX = 0.D0
C
C ======================================================================
C  solution of the Dirac equation for each orbital
C ======================================================================
      DO J = 1,NORB
         DE = DEN(J)
 250     CONTINUE
         CALL SCF0RESLD(NQN(J),NQL(J),NK(J),IMAX,DEN(J),DFL(J),DQ1(J),J,
     &                  DV,DR,DP,DQ,DPAS,Z,NSTOP,NES,TETS,NP,NUC)
         IF ( NSTOP.EQ.0 ) THEN
            VAL = ABS((DEN(J)-DE)/DE)
            IF ( VAL.GT.EMAX ) EMAX = VAL
            NMAX(J) = IMAX
            DO I = 1,NP
               VAL = DGC(I,J) - DP(I)
               IF ( ABS(DP(I)).GT.1.D0 ) VAL = VAL/DP(I)
               IF ( ABS(VAL).GE.ABS(YMAX) ) THEN
                  YMAX = VAL
                  Y = DP(I)
                  YN = DGC(I,J)
               END IF
               VAL = DPC(I,J) - DQ(I)
               IF ( ABS(DQ(I)).GT.1.D0 ) VAL = VAL/DQ(I)
               IF ( ABS(VAL).GE.ABS(YMAX) ) THEN
                  YMAX = VAL
                  Y = DQ(I)
                  YN = DPC(I,J)
               END IF
               DGC(I,J) = DP(I)
               DPC(I,J) = DQ(I)
               D(I) = D(I) + NEL(J)*(DP(I)*DP(I)+DQ(I)*DQ(I))
            END DO
         ELSE IF ( NSTOP.NE.362 .OR. ITER.GE.10 .OR. TETS.GT.TEST ) THEN
            IF ( IPRINT.GT.2 ) WRITE (6,99008) NSTOP,NQN(J),TITRE(J)
            GOTO 100
         ELSE
            TETS = TESTV
            GOTO 250
         END IF
      END DO
C
      CALL SCF0POTSL(DC,D,DP,DR,DPAS,Z,NP,ION)
C
      IF ( NUC.GT.0 ) THEN
         DO I = 1,NUC
            DC(I) = DC(I) + Z/DR(I) + Z*((DR(I)/DR(NUC))**2-3.D0)
     &              /(DR(NUC)+DR(NUC))
         END DO
      END IF
      DO I = 1,NP
         DVAL = ABS(DC(I)-DV(I))
         IF ( (DR(I)*DC(I)).LE.N ) DVAL = -DVAL/DC(I)
         IF ( DVAL.GT.VMAX ) THEN
            VMAX = DVAL
            J = I
         END IF
      END DO
      IF ( IPRINT.GT.2 ) WRITE (6,99003) ITER,VMAX,DR(J),DV(J),DC(J),
     &                          EMAX,YMAX,YN,Y
      IF ( TETS.GT.TEST .OR. EMAX.GT.TESTE .OR. VMAX.GT.TESTV .OR. 
     &     YMAX.GT.TESTY ) THEN
         ITER = ITER + 1
         IF ( ITER.LE.NITER ) THEN
C ======================================================================
C  potential for the next iteration
C ======================================================================
C
            DVAL = 1.D0 - DCOP
            DO I = 1,NP
               DVF(I) = DC(I)
               DV(I) = DVAL*DV(I) + DCOP*DC(I)
            END DO
            GOTO 200
         ELSE
            IF ( IPRINT.GT.2 ) WRITE (6,99004) NITER
            NSTOP = 2
         END IF
      END IF
      WRITE (6,99005)
C ======================================================================
C  the average values of R
C ======================================================================
      DO I = 1,NP
         DVF(I) = DC(I)
         DQ(I) = 0.D0
      END DO
      DVAL = 0.D0
      DO I = 1,NORB
         IM = NMAX(I)
         DVAL = DVAL + NEL(I)*DEN(I)
         DO J = 1,IM
            DC(J) = DGC(J,I)*DGC(J,I) + DPC(J,I)*DPC(J,I)
         END DO
         L = 5
         IF ( ABS(NK(I)).EQ.1 ) L = L - 1
         DO J = 1,L
            DP(J) = DFL(I) + DFL(I)
            IF ( J.EQ.1 ) THEN
               N = 4
            ELSE IF ( J.EQ.2 ) THEN
               N = 2
            ELSE IF ( J.EQ.3 ) THEN
               N = 1
            ELSE IF ( J.EQ.4 ) THEN
               N = -1
            ELSE
               N = -3
            END IF
            CALL SCF0SOMM(DR,DC,DQ,DPAS,DP(J),N,IM)
         END DO
         WRITE (6,99006) NQN(I),TITRE(I),DEN(I)*2,NEL(I),NCORB(I),
     &                   (DP(J),J=2,3)
      END DO
C ======================================================================
C  total energy averaged over sphere
C ======================================================================
C
      DC(1) = 1.D0
      DO I = 1,NP
         DP(I) = D(I)/DR(I)
      END DO
      IF ( NUC.GT.0 ) THEN
         DO I = 1,NUC
            DP(I) = D(I)*(3.D0-DR(I)*DR(I)/(DR(NUC)*DR(NUC)))
     &              /(DR(NUC)+DR(NUC))
         END DO
         DC(1) = 4.D0
      END IF
      CALL SCF0SOMM(DR,DP,DQ,DPAS,DC(1),0,NP)
      DO I = 1,NP
         DP(I) = D(I)*DVF(I)
         D(I) = D(I)*((D(I)*DR(I)+1.D-42)**(1.D0/3.D0))
      END DO
      DC(2) = 3.D0
      DC(3) = 1.D0
      IF ( NUC.NE.0 ) DC(3) = 4.D0
      CALL SCF0SOMM(DR,DP,DQ,DPAS,DC(3),0,NP)
      CALL SCF0SOMM(DR,D,DQ,DPAS,DC(2),-1,NP)
      DC(2) = -3.D0*DC(2)/(105.27578D0**(1.D0/3.D0))
      DC(1) = -Z*DC(1)
      DC(4) = DVAL - DC(3)
      DVAL = DVAL + (DC(1)-DC(3)+(DEXE-DEXV)*DC(2))/2.D0
      DC(3) = (DC(3)-DC(1)-DEXV*DC(2))/2.D0
      DC(2) = DC(2)*DEXE/2.D0
      IF ( NORB.NE.1 ) THEN
C ======================================================================
C  overlapping integrals
C ======================================================================
C
         DO I = 2,NORB
            K = I - 1
            DO J = 1,K
               IF ( NQL(I).EQ.NQL(J) .AND. NK(I).EQ.NK(J) ) THEN
                  IM = NMAX(J)
                  IF ( NMAX(I).LT.IM ) IM = NMAX(I)
                  DO L = 1,IM
                     DQ(L) = DPC(L,I)*DPC(L,J)
                     DC(L) = DGC(L,I)*DGC(L,J)
                  END DO
                  DVAL = DFL(I) + DFL(J)
                  CALL SCF0SOMM(DR,DC,DQ,DPAS,DVAL,0,IM)
                  IF ( DVAL.GT.1.D-3 .AND. IPRINT.GT.2 ) WRITE (6,99007)
     &                 NQN(I),TITRE(I),NQN(J),TITRE(J),DVAL
               END IF
            END DO
         END DO
      END IF
C ======================================================================
C  calculation od the core and total densities
C ======================================================================
      ICEL = 0
      DO J = 1,NORB
         IF ( NCORB(J).NE.0 ) ICEL = ICEL + NEL(J)
      END DO
      DR1AT = DR(1)
      DO I = 1,251
C        ROC(I)=0.D0
         RHOAT(I) = 0.D0
         RHOSAT(I) = 0.D0
      END DO
      DO J = 1,NORB
C       IF(NCORB(J).EQ.1) GO TO 400
         DO I = 1,251
            RHOAT(I) = RHOAT(I) + NEL(J)
     &                 *(DGC(I,J)*DGC(I,J)+DPC(I,J)*DPC(I,J))
            RHOSAT(I) = RHOSAT(I) + MAG(J)
     &                  *(DGC(I,J)*DGC(I,J)+DPC(I,J)*DPC(I,J))
         END DO
      END DO
C
C------------------------------------- store density for valence orbital
C
      LVAL = -1
      IF ( (IZ.GE.21 .AND. IZ.LE.29) .OR. (IZ.GE.39 .AND. IZ.LE.47) .OR. 
     &     (IZ.GE.57 .AND. IZ.LE.80) .OR. (IZ.GE.89 .AND. IZ.LE.111) )
     &     LVAL = 2
C
      J = NORB + 1
      DO JJ = 1,NORB
         J = J - 1
         IF ( NEL(J).GT.0 ) THEN
            IF ( LVAL.LT.0 ) THEN
               JVAL = J
               GOTO 300
            ELSE IF ( NQL(J).EQ.LVAL ) THEN
               JVAL = J
               GOTO 300
            END IF
         END IF
      END DO
C
      CALL STOP_MESSAGE(ROUTINE,'quantum number inconsitent')
C
 300  CONTINUE
      DO I = 1,251
         RHOVAT(I) = DGC(I,JVAL)*DGC(I,JVAL) + DPC(I,JVAL)*DPC(I,JVAL)
         VR2AT(I) = 2D0*DV(I)*DR(I)*DR(I)
      END DO
C
99001 FORMAT (5X,10A4/)
99002 FORMAT (5X,'ITER',4X,'DVMAX',10X,'R',14X,'VN-1',13X,'VN',10X,
     &        'DEMAX',6X,'DPMAX',9X,'PN-1',13X,'PN')
99003 FORMAT (5X,I5,1P,E11.2,3(1P,E16.6),2(1P,E11.2),2(1P,E16.6))
99004 FORMAT (5X,'Number of iterations is exceeded',i4)
99005 FORMAT (5X,'Level  Energy        Occ. Val.      R**2          R')
99006 FORMAT (5X,I2,A2,1X,1P,E14.7,2X,I2,3X,I2,2X,1P,E14.7,2X,1P,E14.7)
99007 FORMAT (34X,I1,A2,I3,A2,F19.7)
99008 FORMAT ('  NSTOP=',I4,'  dlq OPbiTAli   ',I3,A2)
      END
C*==scf0insld.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0INSLD(IPRINT,IZ,DEN,DQ1,DFL,NQN,NQL,NK,NMAX,NEL,
     &                     NORB,MAG,DV,DR,DPAS,Z,NSTOP,NES,NP,NUC,DEXV,
     &                     DEXE,DCOP,TITRE,BAR,TEST,TESTE,TESTY,TESTV,
     &                     NITER,ION,NCORB)
C **********************************************************************
C *                                                                    *
C *  Reading of start potential                                        *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0INSLD')
C
C Dummy arguments
C
      REAL*8 DCOP,DEXE,DEXV,DPAS,TEST,TESTE,TESTV,TESTY,Z
      INTEGER ION,IPRINT,IZ,NES,NITER,NORB,NP,NSTOP,NUC
      CHARACTER*4 BAR(10),TITRE(30)
      REAL*8 DEN(30),DFL(30),DQ1(30),DR(251),DV(251)
      INTEGER MAG(30),NCORB(30),NEL(30),NK(30),NMAX(30),NQL(30),NQN(30)
C
C Local variables
C
      REAL*8 D1,DR1,DVAL,DVC,R,VAL
      CHARACTER*4 ENDA,ITXCH(3,2),TTIRE(9)
      REAL*8 FPOT
      INTEGER I,ICH,IMCH,ION1,IZ1,J,K,L,NKAI,NUC1
C
C*** End of declarations rewritten by SPAG
C
      DATA TTIRE/'S   ','P*  ','P   ','D*  ','D   ','F*  ','F   ',
     &     'G*  ','G   '/,ENDA/'****'/
      DATA ITXCH/'BART','H-HE','DIN ','X-AL','PHA ','    '/
C
      IF ( NSTOP.EQ.1 ) THEN
         DPAS = 0.05D0
         DR1 = 1D-2
C
         NES = 15
         NITER = 50
         TESTE = 5.D-06
         TESTY = 1.D-05
         TESTV = 1.D-05
         TEST = 1.D-07
         NP = 251
         NSTOP = 30
         DEXV = 1.D0
         DEXE = 1.5D0
         DCOP = 0.3D0
         DVC = 137.0359895D0            !137.0373D0
         IZ1 = 0
         ION1 = 0
         NUC1 = -1
         ION = 0
         I = 0
         J = 0
         L = 0
         K = 0
         NUC = 0
         DVAL = 0.D0
         DO I = 1,7
            BAR(I) = '----'
         END DO
C
C ======================================================================
C BAR TITRE contains  24 symbols
C ======================================================================
         IF ( BAR(1).EQ.ENDA )
     &         CALL STOP_MESSAGE(ROUTINE,'BAR(1).EQ.ENDA')
         I = 0
         CALL SCF0DEFORB(NQN,NK,NEL,NCORB,NORB,IZ,MAG)
C
C ======================================================================
C IZ    - atomic number
C ION   - number of electrons (ION=IZ)
C NORB  - number of orbitals
C ICUT = 0  --> UL ON CORRIGE LE POTENTIEL EN -(ION+1)/R
C IPRAT = 0 --> PRATT procedure is performed
C I - number of points for integration (I = 251 if I taken as I=0)
C J - number of iterations to fit the energy (J = 15 if J taken as J=0)
C K - number of iterations (K = 50 if K taken as K=0)
C L=0  --> standard tolerance
C NUC > 0 --> the finite nucleus size
C ======================================================================
C
         DO ICH = 8,10
            IMCH = ICH - 7
            BAR(ICH) = ITXCH(IMCH,1)
         END DO
         IF ( NORB.LE.NSTOP ) THEN
            IF ( I.GT.0 ) THEN
               I = 2*(I/2) + 1
               IF ( I.LE.NP ) THEN
                  NP = I
               ELSE
                  IF ( IPRINT.GT.2 ) WRITE (6,99011) I
                  GOTO 100
               END IF
            END IF
            IF ( J.GT.0 ) NES = J
            IF ( K.GT.0 ) NITER = K
C
C ======================================================================
C DEXV - coefficient for the exchange potential  (DEXV = 1.0 for  SLATER)
C DEXE - coefficient for the exchange energy
C DEXV should be: DEXV = 2.*DEXE/3. for virial theorem
C
C TEST  - energy tolerance for SCF0RESLD
C TESTE - tolerance for selfconsitancy of single-electron energies
C TESTE - tolerance for selfconsitancy of wave functions
C TESTE - tolerance for selfconsitancy of potential
C ======================================================================
            DCOP = 0.3D0
C
C VI(N+1)=(1.-DCOP)*VI(N)+DCOP*VF(N)
C
            Z = DFLOAT(IZ)
            IF ( NUC.GT.0 ) THEN
C
C ======================================================================
C DVAL - atmic masse when NUC > 0
C ======================================================================
C
               DVAL = Z*(DVAL**(1.D0/3.D0))*2.2677D-05/EXP(4.D0*DPAS)
               IF ( DVAL.LE.DR1 ) THEN
                  DR1 = DVAL
                  NUC = 5
                  GOTO 20
               ELSE
                  DVAL = DVAL*EXP(4.D0*DPAS)
                  DO I = 6,NP
                     D1 = DR1*EXP(DFLOAT(I-1)*DPAS)
                     IF ( D1.GE.DVAL ) GOTO 10
                  END DO
                  IF ( IPRINT.GT.2 ) WRITE (6,99015)
                  GOTO 100
               END IF
 10            CONTINUE
               NUC = I
               DR1 = DR1*DVAL/D1
            END IF
 20         CONTINUE
            IF ( IPRINT.GT.2 ) THEN
               WRITE (6,99002) IZ,ION,NITER,TESTE,TESTY,TESTV
               WRITE (6,99009) NP,DR1,IZ,DPAS
               WRITE (6,99010) TEST,NES
               WRITE (6,99001)
               WRITE (6,99012) DEXV,DEXE
            END IF
            K = 0
            DVAL = Z*Z/(DVC*DVC)
            IF ( NUC.GT.0 .AND. IPRINT.GT.2 ) WRITE (6,99014)
C
C ------------------------------------------------ orbital resolved data
            DO I = 1,NORB
C
C ======================================================================
C DEN   - orbital energy in Landau units (<0)
C NQN   - main quantum number
C NK    - KAPPA quantum number
C NEL   - occupation of the orbital
C NCORB - valency indication (NCORB = 1 for BAl.)
C ======================================================================
C
               K = K + NEL(I)
C     IF (DEN(I)) 19,18,18
               DEN(I) = -Z*Z/(4.D0*NQN(I)*NQN(I))
               NQL(I) = ABS(NK(I))
               IF ( NK(I).LT.0 ) NQL(I) = NQL(I) - 1
               IF ( NUC.GT.0 ) THEN
                  NKAI = ABS(NK(I))
                  DFL(I) = DFLOAT(NKAI)
               ELSE
                  DFL(I) = NK(I)*NK(I)
                  DFL(I) = SQRT(DFL(I)-DVAL)
               END IF
               L = 2*ABS(NK(I))
               IF ( NQL(I).LT.NQN(I) .AND. NEL(I).LE.L .AND. NQN(I)
     &              .GT.0 .AND. NQL(I).LE.4 ) THEN
                  J = NQL(I) + ABS(NK(I))
                  TITRE(I) = TTIRE(J)
                  IF ( IPRINT.GT.2 ) WRITE (6,99006) NQN(I),TITRE(I),
     &                 NCORB(I),NEL(I),DEN(I)
               ELSE
                  IF ( IPRINT.GT.2 ) WRITE (6,99004) DEN(I),NQN(I),
     &                 NQL(I),J,NEL(I)
                  GOTO 100
               END IF
            END DO
            IF ( K.EQ.(IZ-ION) ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99013) DCOP
               IF ( NUC.EQ.NUC1 ) THEN
                  IF ( IZ.EQ.IZ1 .AND. ION.EQ.ION1 ) GOTO 50
                  IF ( IZ.EQ.IZ1 ) GOTO 40
               END IF
               DR(1) = DR1/Z
               DO I = 2,NP
                  DR(I) = DR(1)*EXP(DFLOAT(I-1)*DPAS)
               END DO
            ELSE
               IF ( IPRINT.GT.2 ) WRITE (6,99005)
               GOTO 100
            END IF
C
C --------------------------------------------------- starting potential
C
 40         CONTINUE
            VAL = -ION - 1
            IF ( IZ.NE.IZ1 .OR. ION.LE.ION1 .OR. NUC.NE.NUC1 ) THEN
               DO I = 1,NP
                  R = DR(I)
                  DV(I) = FPOT(R,Z,VAL)
               END DO
               IF ( NUC.GT.0 ) THEN
                  DO I = 1,NUC
                     DV(I) = DV(I) + Z/DR(I)
     &                       + Z*((DR(I)/DR(NUC))**2-3.D0)
     &                       /(DR(NUC)+DR(NUC))
                  END DO
               END IF
               GOTO 50
            END IF
C
            DO I = 1,NP
               IF ( (DR(I)*DV(I)).GT.VAL ) DV(I) = VAL/DR(I)
            END DO
            VAL = Z + DV(1)*DR(1)
            IF ( NUC.GT.0 ) VAL = Z + DV(NUC)*DR(NUC)
            IF ( ABS(VAL).GE.0.1D0 ) THEN
               IF ( IPRINT.GT.2 ) WRITE (6,99007)
               GOTO 100
            END IF
         ELSE
            IF ( IPRINT.GT.2 ) WRITE (6,99003) NORB
            GOTO 100
         END IF
 50      CONTINUE
         IF ( NORB.NE.1 ) THEN
            DO I = 2,NORB
               K = I - 1
               DO J = 1,K
                  IF ( NQN(I).EQ.NQN(J) .AND. NK(I).EQ.NK(J) ) THEN
                     IF ( IPRINT.GT.2 ) WRITE (6,99008)
                     GOTO 100
                  END IF
               END DO
            END DO
         END IF
         IZ1 = IZ
         ION1 = ION
         NUC1 = NUC
         DO I = 1,NORB
            NMAX(I) = NP
            L = 1
            J = NQN(I) - NQL(I)
            IF ( (J-2*(J/2)).EQ.0 ) L = -L
C
C MKO       NKAI = ABS(NK(I))
C MKO       DQ1(I) = DFLOAT(L*NK(I)/NKAI)
C MKO following line more stable for vectorisation of loop
C
            DQ1(I) = DFLOAT(L*SIGN(1,NK(I)))
            IF ( NUC.NE.0 ) THEN
               IF ( NK(I).LT.0 ) DQ1(I) = DQ1(I)*(NK(I)-DFL(I))*DVC/Z
            END IF
         END DO
         RETURN
      END IF
 100  CONTINUE
      IF ( IPRINT.GT.2 ) WRITE (6,99016)
      IF ( BAR(1).EQ.ENDA ) THEN
      END IF
      NSTOP = 1
99001 FORMAT (' Exchange: Barth-Hedin   ')
99002 FORMAT (' Atomic number  ',I3,'   Ionicity  ',I2/1X,
     &        'Maximal number of iterations ',
     &        I4/' Precision in energy ',1PE9.2/14x,'wave function ',
     &        1PE9.2/14X,'potential',1PE9.2/)
99003 FORMAT (' NORB=',I3,'too large ************** ')
99004 FORMAT (' Input error       ',E15.8,I1,2I2)
99005 FORMAT (' Erroneous number of electrons  **************')
99006 FORMAT (7X,I1,A2,2I8,1PE23.7)
99007 FORMAT (' Error in potential   ')
99008 FORMAT (' Bad configuration ')
99009 FORMAT (' Integration will be in ',i3,
     &        ' points'/' the first point is ',f7.4,'/',i2,'  step is ',
     &        f7.4,/)
99010 FORMAT (' In <SCF0RESLD> the relative precision in energy',' is ',
     &        1pe9.2,/' the number of attempts ',i3,/)
99011 FORMAT ('  NP=',I3,' is too large *******')
99012 FORMAT ('  Exchange: Slater X-alpha, DEXV=',F8.4,'     DEXE=',
     &        F8.4,/)
99013 FORMAT (' Potential is mixed with  ADMIX=',1PE14.7,/)
99014 FORMAT (10X,' Finite nucleus'/)
99015 FORMAT (' Error in atomic weight  **************** ')
99016 FORMAT (' The next case ')
      END
C*==scf0somm.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0SOMM(DR,DP,DQ,DPAS,DA,M,NP)
C **********************************************************************
C *                                                                    *
C *  integration using Simpson method                                  *
C *  of (DP+DQ)*DR**M  from 0 to R=DR(NP)                              *
C *  DPAS - exponential step for R->0   (DP+DQ)=CTE*R**DA              *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DA,DPAS
      INTEGER M,NP
      REAL*8 DP(*),DQ(*),DR(*)
C
C Local variables
C
      REAL*8 D1,DB,DC,DL
      INTEGER I,MM
C
C*** End of declarations rewritten by SPAG
C
      MM = M + 1
      D1 = DA + DFLOAT(MM)
      DA = 0.D0
      DB = 0.D0
      DO I = 1,NP
         DL = DR(I)**MM
         IF ( I.NE.1 .AND. I.NE.NP ) THEN
            DL = DL + DL
            IF ( (I-2*(I/2)).EQ.0 ) DL = DL + DL
         END IF
         DC = DP(I)*DL
         IF ( DC.LT.0 ) THEN
            DB = DB + DC
         ELSE IF ( ABS(DC).GT.1D-12 ) THEN
            DA = DA + DC
         END IF
         DC = DQ(I)*DL
         IF ( DC.LT.0 ) THEN
            DB = DB + DC
         ELSE IF ( ABS(DC).GT.1D-12 ) THEN
            DA = DA + DC
         END IF
      END DO
      DA = DPAS*(DA+DB)/3.D0
      DC = EXP(DPAS) - 1.D0
      DB = D1*(D1+1.D0)*DC*EXP((D1-1.D0)*DPAS)
      DB = DR(1)*(DR(2)**M)/DB
      DC = (DR(1)**MM)*(1.D0+1.D0/(DC*(D1+1.D0)))/D1
      DA = DA + DC*(DP(1)+DQ(1)) - DB*(DP(2)+DQ(2))
      END
C*==fpot.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      FUNCTION FPOT(R,Z,WA)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 R,WA,Z
      REAL*8 FPOT
C
C Local variables
C
      REAL*8 S,WC,WD,WE
C
C*** End of declarations rewritten by SPAG
C
C
C ======================================================================
C     Tomas-Fermi potential in the R point
C Z  - atomic number
C WA - number of electrons -Z-1
C ======================================================================
C
      S = Z
      WC = SQRT((R*(S+WA)**(1.D0/3.D0))/0.8853D0)
      WD = WC*(0.60112D0*WC+1.81061D0) + 1.D0
      WE = WC*(WC*(WC*(WC*(0.04793D0*WC+0.21465D0)+0.77112D0)+1.39515D0)
     &     +1.81061D0) + 1.D0
      WC = (Z+WA)*(WD/WE)**2 - WA
      FPOT = -WC/R
      END
C*==scf0potsl.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0POTSL(DV,D,DP,DR,DPAS,Z,NP,ION)
C **********************************************************************
C *                                                                    *
C *       4-point integration of potential                             *
C *   DV  - potential                                                  *
C *   D   - density                                                    *
C *   DP  - BLOC DE TRAVAIL                                            *
C *   DR  - radial mesh                                                *
C *   DPAS - exponential step                                          *
C *   DEXV - multiliing factor for the exchange                        *
C *   Z    - atomic number                                             *
C *   NP   - number of points                                          *
C *   ION  - number of electrons (ION = Z)                             *
C *   SI ICUT EST NUL ON CORRIGE EVENTUELLEMENT                        *
C *        LE POTENTIEL EN -(ION+1)/R                                  *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DPAS,Z
      INTEGER ION,NP
      REAL*8 D(*),DP(*),DR(*),DV(*)
C
C Local variables
C
      REAL*8 DARS,DAS,DBRS,DCF,DCP,DLO,DLO2,DNY,DR2,DS,DS1,DSF,DSF2,
     &       DSF3,DSP,DSP2,DSP3,DVXC
      INTEGER I,J,K
C
C*** End of declarations rewritten by SPAG
C
      DAS = DPAS/24.D0
      DO I = 1,NP
         DV(I) = D(I)*DR(I)
      END DO
      DLO = EXP(DPAS)
      DLO2 = DLO*DLO
      DP(2) = DR(1)*(D(2)-D(1)*DLO2)/(12.D0*(DLO-1.D0))
      DP(1) = DV(1)/3.D0 - DP(2)/DLO2
      DP(2) = DV(2)/3.D0 - DP(2)*DLO2
      J = NP - 1
      DO I = 3,J
         DP(I) = DP(I-1) + DAS*(13.D0*(DV(I)+DV(I-1))-(DV(I-2)+DV(I+1)))
      END DO
      DP(NP) = DP(J)
      DV(J) = DP(J)
      DV(NP) = DP(J)
      DO I = 3,J
         K = NP + 1 - I
         DV(K) = DV(K+1)
     &           /DLO + DAS*(13.D0*(DP(K+1)/DLO+DP(K))-(DP(K+2)/DLO2+
     &           DP(K-1)*DLO))
      END DO
      DV(1) = DV(3)/DLO2 + DPAS*(DP(1)+4.D0*DP(2)/DLO+DP(3)/DLO2)/3.D0
      DLO = -DFLOAT(ION+1)
      DO I = 1,NP
C ------------------------------------------------ BARTH-HEDIN potential
C
         DR2 = DR(I)**2
         DS1 = (D(I)/(3.D0*DR2)+1.D-42)**(1.D0/3.D0)
C
         IF ( ABS(DS1).LT.1.D-10 ) THEN
            DVXC = 0.D0
         ELSE
            DS = 1.D0/DS1
            DSF = DS/75.D0
            DSF2 = DSF*DSF
            DSF3 = DSF2*DSF
            DSP = DS/30.D0
            DSP2 = DSP*DSP
            DSP3 = DSP2*DSP
            DCF = (1.D0+DSF3)*LOG(1.D0+1.D0/DSF) + 0.5D0*DSF - DSF2 - 
     &            0.3333333333D0
            DCP = (1.D0+DSP3)*LOG(1.D0+1.D0/DSP) + 0.5D0*DSP - DSP2 - 
     &            0.3333333333D0
            DNY = 5.1297628D0*(0.0504D0*DCP-0.0254D0*DCF)
            DARS = -1.22177412D0/DS + DNY
            DBRS = -0.0504D0*LOG(1.D0+30.D0/DS) - DNY
            DVXC = DARS + DBRS
         END IF
         DV(I) = DV(I) - (Z-0.5D0*DR(I)*DVXC)
C
         IF ( DV(I).GT.DLO ) DV(I) = DLO
         DV(I) = DV(I)/DR(I)
      END DO
      END
C*==scf0deforb.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0DEFORB(XN,NK,XZ,IVAL,NORB,IZN,MAG)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0DEFORB')
C
C Dummy arguments
C
      INTEGER IZN,NORB
      INTEGER IVAL(30),MAG(30),NK(30),XN(30),XZ(30)
C
C Local variables
C
      INTEGER I,IXT,IXV,IXZ1(105),IXZ2(105),IXZ3(105),IXZ4(105),
     &        IXZ5(105),IZ,IZV,J,JF1,JFUL,JFULR(105),JREL(105),NMAG
      REAL*8 XJ(35),XL(35),XMAG,XML,XXJ(31),XXL(31),XXN(31),XXZ(31),ZN
C
C*** End of declarations rewritten by SPAG
C
C=====================================================================
C
C     STORED DATA FOR FULL ORBITALS (AS THEY ARE NUMBERED IN 'ATOMGL'):
C
C     STORED DATA FOR RELATIVISTIC ORBITALS:
      DATA XXN/1D0,2D0,2D0,2D0,3D0,3D0,3D0,4D0,3D0,3D0,4D0,4D0,5D0,4D0,
     &     4D0,5D0,5D0,6D0,4D0,4D0,5D0,5D0,6D0,6D0,7D0,5D0,5D0,6D0,6D0,
     &     7D0,7D0/
      DATA XXL/0D0,0D0,1D0,1D0,0D0,1D0,1D0,0D0,2D0,2D0,1D0,1D0,0D0,2D0,
     &     2D0,1D0,1D0,0D0,3D0,3D0,2D0,2D0,1D0,1D0,0D0,3D0,3D0,2D0,2D0,
     &     1D0,1D0/
      DATA XXJ/0.5D0,0.5D0,0.5D0,1.5D0,0.5D0,0.5D0,1.5D0,0.5D0,1.5D0,
     &     2.5D0,0.5D0,1.5D0,0.5D0,1.5D0,2.5D0,0.5D0,1.5D0,0.5D0,2.5D0,
     &     3.5D0,1.5D0,2.5D0,0.5D0,1.5D0,0.5D0,2.5D0,3.5D0,1.5D0,2.5D0,
     &     0.5D0,1.5D0/
      DATA XXZ/2D0,2D0,2D0,4D0,2D0,2D0,4D0,2D0,4D0,6D0,2D0,4D0,2D0,4D0,
     &     6D0,2D0,4D0,2D0,6D0,8D0,4D0,6D0,2D0,4D0,2D0,6D0,8D0,4D0,6D0,
     &     2D0,4D0/
      DATA JFULR/0,1,1,2,2,3,3,3,3,4,4,5,5,6,6,6,6,7,7,8,8,8,8,7,9,9,9,
     &     9,7,10,10,11,11,11,11,12,12,13,13,13,12,12,14,12,12,12,12,15,
     &     15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,19,19,19,
     &     19,20,20,20,20,21,21,21,21,17,17,22,22,23,23,23,23,24,24,25,
     &     25,25,25,25,25,26,26,26,26,26,26,26,26,27,27,27,27/
      DATA JREL/1,1,2,2,3,3,4,4,4,4,5,5,6,6,7,7,7,7,8,8,9,9,9,10,10,10,
     &     10,10,10,10,11,11,12,12,12,12,13,13,14,14,14,15,15,15,15,15,
     &     15,15,16,16,17,17,17,17,18,18,21,21,19,19,19,19,20,21,20,20,
     &     20,20,20,20,21,21,21,21,22,22,22,22,22,22,23,23,24,24,24,24,
     &     25,25,28,28,28,28,28,26,27,28,28,27,27,27,27,27,28,28,28/
      DATA IXZ1/1,0,1,0,1,0,1,2,3,0,1,0,1,0,1,2,3,0,1,0,1,2,3,1,1,2,3,4,
     &     1,0,1,0,1,2,3,0,1,0,1,2,1,1,1,1,1,0,1,0,1,0,1,2,3,0,1,0,0,1,
     &     3,4,5,0,1,1,3,4,5,6,7,0,1,2,3,0,1,2,3,1,1,0,1,0,1,2,3,0,1,0,
     &     0,0,2,3,4,0,1,1,2,4,5,6,7,0,1,2,3/
      DATA IXZ2/23*0,4,4*0,4,11*0,4,4,0,4,4,4,4,16*0,1,13*0,6,6,16*0,1,
     &     1,8*0/
      DATA IXZ3/23*0,1,4*0,6,12*0,1,0,3,4,6,6,9*0,1,1,19*0,8,8,9*0,1,2,
     &     1,1,1,12*0/
      DATA IXZ4/77*0,4,4,26*0/
      DATA IXZ5/77*0,5,6,26*0/
C=====================================================================
      ZN = IZN
      IZ = INT(ZN+0.000001D0)
      IZV = IZN
C------ J-REPRESENTATION ------
      JFUL = JFULR(IZ)
C
      J = JREL(IZ)
      DO I = 1,JFUL
         XN(I) = NINT(XXN(I))
         XL(I) = XXL(I)
         XJ(I) = XXJ(I)
         XZ(I) = NINT(XXZ(I))
      END DO
C
      XZ(JFUL+1) = IXZ1(IZ)
      XZ(JFUL+2) = IXZ2(IZ)
      XZ(JFUL+3) = IXZ3(IZ)
      XZ(JFUL+4) = IXZ4(IZ)
      XZ(JFUL+5) = IXZ5(IZ)
C
C--------------------------------
      IF ( JFUL.NE.J ) THEN
         JF1 = JFUL + 1
         DO I = JF1,J
            XN(I) = NINT(XXN(I))
            XL(I) = XXL(I)
            XJ(I) = XXJ(I)
C--------------------------------
         END DO
      END IF
C=====================================================================
      IXV = 0
      IXT = 0
      DO I = J,1, - 1
         IXT = IXT + XZ(I)
         NK(I) = NINT(XL(I))
         MAG(I) = 0
         XMAG = 0
         IF ( XJ(I).GT.XL(I) ) NK(I) = NINT(-XL(I)-1)
         IF ( IXV+XZ(I).LE.IZV ) THEN
            IVAL(I) = 1
            IXV = IXV + XZ(I)
            IF ( NK(I).GT.1 ) THEN
               XML = XZ(I) + XZ(I+1)
               NMAG = NK(I)*2 + 1
               XMAG = NMAG - ABS(XML-NMAG)
               MAG(I) = NINT(XMAG)
            END IF
         ELSE
            IVAL(I) = 0
         END IF
      END DO
      IF ( IXV.NE.IZV ) CALL STOP_MESSAGE(ROUTINE,'IXV.NE.IZV')
      NORB = J
      END
C*==scf0resld.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0RESLD(NQN,NQL,NK,IMAX,DE,DFL,DQ1,JC1,DV,DR,DP,DQ,
     &                     DPAS,Z,NSTOP,NES,TEST,NP,NUC)
C **********************************************************************
C *                                                                    *
C *   Solution of the Dirac equation                                   *
C *   NQN   - main quantum number                                      *
C *   NQL   - orbital quantum number                                   *
C *   NK    - KAPPA quantum number                                     *
C *   IMAX  - the last point for wave function tabulation              *
C *   DE    - energy                                                   *
C *   DFL   - the power factor in the wave function expansion          *
C *   DQ1   - the DP/DQ ratio in the origin of coordinates             *
C *                                                                    *
C *   DV   - potential in a.u.(<0)                                     *
C *   DR   - radial mesh                                               *
C *   DP   - large component                                           *
C *   DQ   - small component                                           *
C *   DPAS - exponential step                                          *
C *   Z    - atomic number                                             *
C *   NSTOP - numerical integration check                              *
C *   NES   - maximal number of iterations to find the energy          *
C *   TEST  - tolerance for the energy                                 *
C *   NP    - maximal number of points                                 *
C *   when NUC .neq. 0 --> finite nucleous size                        *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DE,DFL,DPAS,DQ1,TEST,Z
      INTEGER IMAX,JC1,NES,NK,NP,NQL,NQN,NSTOP,NUC
      REAL*8 DP(251),DQ(251),DR(251),DV(251)
C
C Local variables
C
      REAL*8 DB,DBE,DD,DEP(5),DEQ(5),DK,DKOEF,DM,DPM,DPNO(4,30),DPQ,DQM,
     &       DQNO(4,30),DSAL,DSUM,DVAL,DVC,ELIM,VAL
      INTEGER I,IES,IMAT,IMM,J,JC,K,LLL,M,ND,NOEUD
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
C  DEP,DEQ  - derivatives of DP and DQ
C  DB = ENERGIE/DVC
C  DVC - speed of light in a.u.
C  DSAL=2.*DVC
C  DK  - KAPPA quantum number
C  DM  - (exponential step)/720., DKOEF=1./720.
C ======================================================================
C
      DATA DKOEF/.1388888888888888D-2/
      JC = 1
      IF ( JC1.GT.0 ) JC = JC1
      NSTOP = 0
      DVC = 137.0373D0
      DSAL = DVC + DVC
      IMM = 0
      IES = 0
      DK = NK
      LLL = (NQL*(NQL+1))/2
      ND = 0
      NOEUD = NQN - NQL
      IF ( LLL.NE.0 ) THEN
         ELIM = DV(1) + LLL/(DR(1)*DR(1))
         DO I = 2,NP
            VAL = DV(I) + LLL/(DR(I)*DR(I))
            IF ( VAL.LE.ELIM ) ELIM = VAL
         END DO
         IF ( ELIM.LT.0 ) THEN
            IF ( DE.LE.ELIM ) DE = ELIM*0.5D0
         ELSE
            NSTOP = 17
C
C  2*V+L*(L+1)/R**2  positive everyware
C
            RETURN
         END IF
      ELSE
         ELIM = -Z*Z/(1.5D0*NQN*NQN)
         IF ( DE.LE.ELIM ) DE = ELIM*0.5D0
      END IF
 100  CONTINUE
      IF ( IMM.NE.1 ) THEN
         DO I = 7,NP,2
            IMAT = NP + 1 - I
            IF ( (DV(IMAT)+DFLOAT(LLL)/(DR(IMAT)*DR(IMAT))-DE).LE.0.D0 )
     &           EXIT
         END DO
         IF ( IMAT.LE.5 ) THEN
            DE = DE*0.5D0
            IF ( DE.LT.-TEST .AND. ND.LE.NOEUD ) GOTO 100
            NSTOP = 28
C
C  2*V+L*(L+1)/R**2-2*E  positive everyware
C
            RETURN
         END IF
      END IF
C ======================================================================
C initial values for the integration over the inner space
C ======================================================================
      DB = DE/DVC
      CALL SCF0INOUH(DP,DQ,DR,DQ1,DFL,DV(1),Z,TEST,NUC,NSTOP,JC,DPNO,
     &               DQNO,DEP,DEQ,DB,DVC,DSAL,DK,DM)
      IF ( NSTOP.NE.0 ) GOTO 99999
C     NSTOP=45
C ======================================================================
C no convergency in the origin of coordinates
C ======================================================================
      ND = 1
      DO I = 1,5
         DVAL = DR(I)**DFL
         IF ( I.NE.1 ) THEN
            IF ( ABS(DP(I-1)).GT.1D-12 ) THEN
               IF ( (DP(I)/DP(I-1)).LE.0.D0 ) ND = ND + 1
            END IF
         END IF
         DP(I) = DP(I)*DVAL
         DQ(I) = DQ(I)*DVAL
         DEP(I) = DEP(I)*DVAL
         DEQ(I) = DEQ(I)*DVAL
      END DO
      K = -1 + 2*(NOEUD-2*(NOEUD/2))
      IF ( (DP(1)*DFLOAT(K)).GT.0.D0 ) THEN
         IF ( (DFLOAT(K)*DFLOAT(NK)*DQ(1)).GE.0.D0 ) THEN
            DM = DPAS*DKOEF
C ======================================================================
C  integration over the inner space
C ======================================================================
C
            DO I = 6,IMAT
               DP(I) = DP(I-1)
               DQ(I) = DQ(I-1)
               CALL SCF0INTH(DP(I),DQ(I),DV(I),DR(I),DEP,DEQ,DB,DVC,
     &                       DSAL,DK,DM)
               IF ( ABS(DP(I-1)).GT.1D-12 ) THEN
                  IF ( (DP(I)/DP(I-1)).LE.0.D0 ) THEN
                     ND = ND + 1
                     IF ( ND.GT.NOEUD ) GOTO 200
                  END IF
               END IF
            END DO
            IF ( ND.EQ.NOEUD ) THEN
C ======================================================================
C initial values for the integration over the outer space
C ======================================================================
C
               DQM = DQ(IMAT)
               DPM = DP(IMAT)
               IF ( IMM.NE.1 ) THEN
                  DO I = 1,NP,2
                     IMAX = NP + 1 - I
                     IF ( ((DV(IMAX)-DE)*DR(IMAX)*DR(IMAX)).LE.300.D0 )
     &                    EXIT
                  END DO
               END IF
               DD = SQRT(-DE*(2.D0+DB/DVC))
               DPQ = -DD/(DSAL+DB)
               DM = -DM
               DO I = 1,5
                  J = IMAX + 1 - I
                  DP(J) = EXP(-DD*DR(J))
                  DEP(I) = -DD*DP(J)*DR(J)
                  DQ(J) = DPQ*DP(J)
                  DEQ(I) = DPQ*DEP(I)
               END DO
               M = IMAX - 5
C ======================================================================
C  integration over the outer space
C ======================================================================
C
               DO I = IMAT,M
                  J = M + IMAT - I
                  DP(J) = DP(J+1)
                  DQ(J) = DQ(J+1)
                  CALL SCF0INTH(DP(J),DQ(J),DV(J),DR(J),DEP,DEQ,DB,DVC,
     &                          DSAL,DK,DM)
               END DO
C ======================================================================
C  matching of the large component
C ======================================================================
C
               DVAL = DPM/DP(IMAT)
               IF ( DVAL.GT.0.D0 ) THEN
                  DO I = IMAT,IMAX
                     DP(I) = DP(I)*DVAL
                     DQ(I) = DQ(I)*DVAL
                  END DO
C  norm factor calculation
C **********************************************************************
                  DSUM = 3.D0*DR(1)*(DP(1)**2+DQ(1)**2)
     &                   /(DPAS*(DFL+DFL+1.D0))
                  DO I = 3,IMAX,2
                     DSUM = DSUM + DR(I)*(DP(I)**2+DQ(I)**2)
     &                      + 4.D0*DR(I-1)*(DP(I-1)**2+DQ(I-1)**2)
     &                      + DR(I-2)*(DP(I-2)**2+DQ(I-2)**2)
                  END DO
                  DSUM = DPAS*(DSUM+DR(IMAT)*(DQM*DQM-DQ(IMAT)*DQ(IMAT))
     &                   )*0.3333333333333333D0
C  energy calculation
C **********************************************************************
                  DBE = DP(IMAT)*(DQM-DQ(IMAT))*DVC/DSUM
                  IMM = 0
                  VAL = ABS(DBE/DE)
                  IF ( VAL.LE.TEST ) THEN
                     IF ( JC1.LE.0 .AND. IMAX.GT.(-JC1) ) THEN
                        IMAX = -JC1
                        DSUM = 3.D0*DR(1)*(DP(1)**2+DQ(1)**2)
     &                         /(DPAS*(DFL+DFL+1.D0))
                        DO I = 3,IMAX,2
                           DSUM = DSUM + DR(I)*(DP(I)**2+DQ(I)**2)
     &                            + 4.D0*DR(I-1)*(DP(I-1)**2+DQ(I-1)**2)
     &                            + DR(I-2)*(DP(I-2)**2+DQ(I-2)**2)
                        END DO
                        DSUM = DPAS*(DSUM+DR(IMAT)
     &                         *(DQM*DQM-DQ(IMAT)*DQ(IMAT)))
     &                         *0.3333333333333333D0
                     END IF
                     DSUM = SQRT(DSUM)
                     DQ1 = DQ1/DSUM
                     DO I = 1,IMAX
                        DP(I) = DP(I)/DSUM
                        DQ(I) = DQ(I)/DSUM
                     END DO
                     DO I = 1,4
                        DPNO(I,JC) = DPNO(I,JC)/DSUM
                        DQNO(I,JC) = DQNO(I,JC)/DSUM
                     END DO
                     IF ( IMAX.NE.NP ) THEN
                        J = IMAX + 1
                        DO I = J,NP
                           DP(I) = 0.D0
                           DQ(I) = 0.D0
                        END DO
                     END IF
                     NSTOP = 0
                     GOTO 99999
                  ELSE
 102                 CONTINUE
                     DVAL = DE + DBE
                     IF ( DVAL.LT.0.D0 ) THEN
                        DE = DVAL
                        IF ( VAL.LE.0.1D0 ) IMM = 1
                        IES = IES + 1
                        IF ( IES.LE.NES ) GOTO 100
                        NSTOP = 362
C  the number of iterations exceeds the maximum
C **********************************************************************
                        RETURN
                     ELSE
                        DBE = DBE*0.5D0
                        VAL = VAL*0.5D0
                        IF ( VAL.GT.TEST ) GOTO 102
                        NSTOP = 345
C  zero energy
C **********************************************************************
                        RETURN
                     END IF
                  END IF
               ELSE
                  NSTOP = 312
C  error in the sing of large component
C **********************************************************************
                  RETURN
               END IF
            ELSE
               DE = 0.8D0*DE
               IF ( DE.LT.-TEST ) GOTO 100
               NSTOP = 206
C  number of nodes is too small
C **********************************************************************
               RETURN
            END IF
         END IF
      END IF
      NSTOP = 53
C  error in the expantion in the origin of coordinates
C **********************************************************************
      RETURN
 200  CONTINUE
      DE = 1.2D0*DE
      IF ( DE.GT.ELIM ) GOTO 100
      NSTOP = 210
C  number of nodes is too big
C **********************************************************************
      RETURN
99999 CONTINUE
      END
C*==scf0inouh.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0INOUH(DP,DQ,DR,DQ1,DFL,DV,Z,TEST,NUC,NSTOP,JC,DPNO,
     &                     DQNO,DEP,DEQ,DD,DVC,DSAL,DK,DM)
C **********************************************************************
C *                                                                    *
C * Initial values for the integration over the inner space            *
C *   DP   - large component                                           *
C *   DQ   - small component                                           *
C *   DR   - radial mesh                                               *
C *   DQ1   - the DP/DQ ratio in the origin of coordinates             *
C *   DFL   - the power factor of the main term of the expansion       *
C *           in the origin                                            *
C *   DV   - potential in the first point                              *
C *   Z    - atomic number                                             *
C *   TEST  - tolerance                                                *
C *   when NUC .neq. 0 --> finite nucleous size                        *
C *   NSTOP - check of convergency at the power expansion              *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DD,DFL,DK,DM,DQ1,DSAL,DV,DVC,TEST,Z
      INTEGER JC,NSTOP,NUC
      REAL*8 DEP(5),DEQ(5),DP(*),DPNO(4,30),DQ(*),DQNO(4,30),DR(*)
C
C Local variables
C
      REAL*8 DBE,DEVA1,DEVA2,DEVA3,DPR,DQR,DSUM,DVAL
      INTEGER I,J,M
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
C  DEP,DEQ  - derivatives of DP and DQ
C  DD = ENERGIE/DVC
C  DVC - speed of light in a.u.
C  DSAL=2.*DVC
C  DK  - KAPPA quantum number
C  DM  - (exponential step)/720.
C ======================================================================
C
      DO I = 1,10
         DP(I) = 0.D0
         DQ(I) = 0.D0
      END DO
      IF ( NUC.LE.0 ) THEN
         DVAL = Z/DVC
         DEVA1 = -DVAL
         DEVA2 = DV/DVC + DVAL/DR(1) - DD
         DEVA3 = 0.D0
         IF ( DK.LE.0 ) THEN
            DBE = (DK-DFL)/DVAL
         ELSE
            DBE = DVAL/(DK+DFL)
         END IF
         DQ(10) = DQ1
         DP(10) = DBE*DQ1
      ELSE
         DVAL = DV + Z*(3.D0-DR(1)*DR(1)/(DR(NUC)*DR(NUC)))
     &          /(DR(NUC)+DR(NUC))
         DEVA1 = 0.D0
         DEVA2 = (DVAL-3.D0*Z/(DR(NUC)+DR(NUC)))/DVC - DD
         DEVA3 = Z/(DR(NUC)*DR(NUC)*DR(NUC)*DSAL)
         IF ( DK.LE.0 ) THEN
            DP(10) = DQ1
         ELSE
            DQ(10) = DQ1
         END IF
      END IF
      DO I = 1,5
         DP(I) = DP(10)
         DQ(I) = DQ(10)
         DEP(I) = DP(I)*DFL
         DEQ(I) = DQ(I)*DFL
      END DO
      M = 1
 100  CONTINUE
      DM = M + DFL
      DSUM = DM*DM - DK*DK + DEVA1*DEVA1
      DQR = (DSAL-DEVA2)*DQ(M+9) - DEVA3*DQ(M+7)
      DPR = DEVA2*DP(M+9) + DEVA3*DP(M+7)
      DVAL = ((DM-DK)*DQR-DEVA1*DPR)/DSUM
      DSUM = ((DM+DK)*DPR+DEVA1*DQR)/DSUM
      J = -1
      DO I = 1,5
         DPR = DR(I)**M
         DQR = DSUM*DPR
         DPR = DVAL*DPR
         IF ( M.NE.1 ) THEN
            IF ( ABS(DPR/DP(I)).LE.TEST .AND. ABS(DQR/DQ(I)).LE.TEST )
     &           J = 1
         END IF
         DP(I) = DP(I) + DPR
         DQ(I) = DQ(I) + DQR
         DEP(I) = DEP(I) + DPR*DM
         DEQ(I) = DEQ(I) + DQR*DM
      END DO
      IF ( J.NE.1 ) THEN
         DP(M+10) = DVAL
         DQ(M+10) = DSUM
         M = M + 1
         IF ( M.LE.20 ) GOTO 100
         NSTOP = 45
      END IF
      DO I = 1,4
         DPNO(I,JC) = DP(I+9)
         DQNO(I,JC) = DQ(I+9)
      END DO
      END
C*==scf0inth.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0INTH(DP,DQ,DV,DR,DEP,DEQ,DB,DVC,DSAL,DK,DM)
C **********************************************************************
C *                                                                    *
C * 5-point integration using the Adams method of large componet DP   *
C * and small componet DQ in the point DR                              *
C * DV - potential in this point                                      *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DB,DK,DM,DP,DQ,DR,DSAL,DV,DVC
      REAL*8 DEP(5),DEQ(5)
C
C Local variables
C
      REAL*8 DKOEF1,DKOEF2,DPR,DQR,DSUM
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
C  DEP,DEQ  - derivatives of DP and DQ
C  DD = ENERGIE/DVC
C  DVC - speed of light in a.u.
C  DSAL=2.*DVC
C  DK  - KAPPA quantum number
C  DM  - (exponential step)/720.
C  DKOEF1=405./502., DKOEF2=27./502.
C ======================================================================
C
      DATA DKOEF1/.9462151394422310D0/,DKOEF2/.5378486055776890D-1/
      DPR = DP + DM*((251.D0*DEP(1)+2616.D0*DEP(3)+1901.D0*DEP(5))
     &      -(1274.D0*DEP(2)+2774.D0*DEP(4)))
      DQR = DQ + DM*((251.D0*DEQ(1)+2616.D0*DEQ(3)+1901.D0*DEQ(5))
     &      -(1274.D0*DEQ(2)+2774.D0*DEQ(4)))
      DO I = 2,5
         DEP(I-1) = DEP(I)
         DEQ(I-1) = DEQ(I)
      END DO
      DSUM = (DB-DV/DVC)*DR
      DEP(5) = -DK*DPR + (DSAL*DR+DSUM)*DQR
      DEQ(5) = DK*DQR - DSUM*DPR
      DP = DP + DM*((106.D0*DEP(2)+646.D0*DEP(4)+251.D0*DEP(5))
     &     -(19.D0*DEP(1)+264.D0*DEP(3)))
      DQ = DQ + DM*((106.D0*DEQ(2)+646.D0*DEQ(4)+251.D0*DEQ(5))
     &     -(19.D0*DEQ(1)+264.D0*DEQ(3)))
      DP = DKOEF1*DP + DKOEF2*DPR
      DQ = DKOEF1*DQ + DKOEF2*DQR
      DEP(5) = -DK*DP + (DSAL*DR+DSUM)*DQ
      DEQ(5) = DK*DQ - DSUM*DP
      END
C*==scf0qd.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SCF0QD(H,Y,Z,NDIM)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATES NDIM INTEGRALS                                      *
C   *   BETWEEN X1 AND XN OF FUNCTION Z ( N.LE.NDIM                    *
C   *    XN=X1+(N-1)*H ) FOR THE SPECIAL CASE:                         *
C   *   NDIM=2*ND+1 Y(1)=0 . Z(2*I+1) IS CALCULATED                    *
C   *   BY THE SIMPSON'S METHOD AND Z(2*I+2)=Z(2*I+1)+                 *
C   *   DELT(2*I+2) . TO INCREASE A ACCURACY IT IS                     *
C   *   USED DELT(2*I+2)=DELT4(2*I+2)*(Z(2*I+3)-Z(2*I+1))/             *
C   *   (DELT4(2*I+2)+DELT4(2*I+3)) WHERE DELT4(2*I+2)                 *
C   *   IS CALCULATED BY USING A CUBIC INTERPOLATION                   *
C   *   BETWEEN Y(2*I),Y(2*I+1),Y(2*I+2),Y(2*I+3).                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF0QD')
C
C Dummy arguments
C
      REAL*8 H
      INTEGER NDIM
      REAL*8 Y(NDIM),Z(NDIM)
C
C Local variables
C
      REAL*8 C24,C3,DELT,DELT3,DELT4,DELT5,SUM2,Y1,Y2,Y3,Y4,Y5
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      IF ( (NDIM/2+NDIM/2).EQ.NDIM )
     &      CALL STOP_MESSAGE(ROUTINE,'works only for NDIM=2*ND+1')
C
      C3 = H/3.D0
      C24 = H/24.D0
      Z(1) = 0.D0
      IF ( NDIM.NE.1 ) THEN
         Y2 = Y(1)
         Y3 = Y(2)
         Y4 = Y(3)
         Z(2) = C24*(10*Y2+16*Y3-2*Y4)
         SUM2 = Y2 + 4*Y3 + Y4
         Z(3) = C3*SUM2
         IF ( NDIM.NE.3 ) THEN
            Y5 = Y(4)
C
C the main loop of a integration
C
            DO I = 4,NDIM - 1,2
               Y1 = Y3
               Y2 = Y4
               Y3 = Y5
               Y4 = Y(I+1)
               DELT3 = Y2 + 4*Y3 + Y4
               SUM2 = SUM2 + DELT3
               Z(I+1) = C3*SUM2
               DELT4 = -Y1 + 13*Y2 + 13*Y3 - Y4
               IF ( I.LT.NDIM-1 ) THEN
                  Y5 = Y(I+2)
                  DELT5 = -Y1 + 12*Y2 + 26*Y3 + 12*Y4 - Y5
               ELSE
                  DELT5 = DELT3*8
               END IF
               IF ( DELT5.GT.1.D-20 ) THEN
                  DELT = 8*DELT4*DELT3/DELT5
                  Z(I) = Z(I-1) + C24*DELT
               ELSE
                  Z(I) = Z(I-1) + C24*DELT4
               END IF
            END DO
         END IF
      END IF
      END
