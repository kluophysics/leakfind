C*==chisusc.f    processed by SPAG 6.70Rc at 16:03 on 16 Jan 2017
      SUBROUTINE CHISUSC(IPROCE,CHIPRINT)
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine for susceptibility calculations                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:DDTAUTAUT,TKTKTT,NTKTKMAX
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R
      USE MOD_KSPACE,ONLY:IBZINT,NKTAB,WKTAB,KTAB,NKTABMAX
      USE MOD_SYMMETRY,ONLY:SYMUNITARY,IQORGQP,DROT,NSYMACCEPTED,NSYM,
     &    SYMACCEPTED
      USE MOD_ANGMOM,ONLY:IDOS,NL,LINMAX,NKMMAX,NLMAX,NKKR,TXT_L,IND0Q,
     &    NKMQ,NKM,MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,NQMAX,DROTQ,ITOQ,NOQ,IQAT,ICPA
      USE MOD_TYPES,ONLY:NT,NTMAX,NLT,NAT,IMT,BT,VT,LTXT_T,TXT_T,RHOSPN,
     &    RHOCHR,CONC,DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX,DOBS_TX_GLO,
     &    OBS_TX_GLO
      USE MOD_CALCMODE,ONLY:ITEST,NONMAG,KMROT,ORBPOL,TASK
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NEPATH,IGRID,NETAB,EILOW,EIMAG,
     &    EMAX,EMIN,EFERMI,WETAB,PHASK
      USE MOD_FILES,ONLY:IFILCBWF,LSYSTEM,SYSTEM,LDATSET,DATSET,IOTMP,
     &    IPRINT,WRTAU,RDTAU,CMDUC,FOUND_SECTION,FOUND_REAL,
     &    FOUND_STRING,FOUND_REAL_ARRAY
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:CHI_AU2CGS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NOP,NOBS,ISPN,IORB,IHFI,IHVV
      PARAMETER (NOP=3,NOBS=NOP+2,ISPN=2,IORB=3,IHFI=4,IHVV=5)
      REAL*8 CU
      PARAMETER (CU=CHI_AU2CGS*1D+6)
C
C Dummy arguments
C
      INTEGER CHIPRINT
      INTEGER IPROCE(NEMAX)
C
C Local variables
C
      CHARACTER*80 ABXCFIL,FILNAM,LINE,TAUFIL
      REAL*8 AUXL(NLMAX),AXCN(NRMAX,2,NLMAX,NTMAX),AXCNP(NRMAX,2,NTMAX),
     &       BCOR(NTMAX),BCORS(NTMAX),BCP(NTMAX),BCPS(NTMAX),BCPT(:,:),
     &       BEXTCP,BT0(:,:),BXCMNV(NRMAX,NTMAX,4),BXCNDUM(NRMAX,NTMAX),
     &       BXCNMM(NRMAX,NTMAX),BXCNMN(NRMAX,NTMAX),BXCNNM(NRMAX,NTMAX)
     &       ,BXCNNN(NRMAX,NTMAX),BXCNP(NRMAX,NTMAX),
     &       CHILANLT(0:NLMAX,NTMAX,2),CHIMIX,CHITOL,CPACHNG,
     &       CPACHTAB(NEMAX),ECOR(:),EMM(NRMAX,NTMAX),ENM(NRMAX,NTMAX),
     &       ENN(NRMAX,NTMAX),EXC(:),FCOR(:,:,:),GAMMA(NRMAX,NTMAX),
     &       GAMMAN(NRMAX,NTMAX),GAMMAP(NRMAX,NTMAX),GCOR(:,:,:),
     &       NEFTL(NTMAX,NLMAX),ORBSQINT(0:NLMAX,2),RHO4PI(NRMAX,2),
     &       RHOCHRC(NRMAX,NTMAX),RHORINTL(NRMAX,NLMAX,NTMAX),
     &       RHOSPNC(NRMAX,NTMAX),RMSMN(2,4),RPWOP(NRMAX,NLMAX*2),
     &       RWORK(0:NLMAX,NTMAX,2),SIGNAF(NTMAX),SWP,SZCOR(:),TIME1,
     &       TIME2,VAUX(NRMAX,2),VMTZOFF,WEXC(NRMAX,2)
      LOGICAL BXCCORR,CALCINT,CALCNEF,GETIRRSOL,INSIDE,RDABXC,SLAB,
     &        SURFACE,TREATAF
      COMPLEX*16 CHILAN,CTOTDOS(1),DMATT(NKMMAX,NKMMAX,NTMAX),
     &           DOS1(NTMAX),DOS1L(NLMAX,NTMAX),DOSL1(NLMAX,NTMAX),
     &           DTILT(NKMMAX,NKMMAX,NTMAX),ERYD,P,
     &           RHORTL1(NRMAX,NLMAX,NTMAX,2,NOBS),
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,2),TOTDOS
      CHARACTER*10 CHIORBPOL,STR10
      INTEGER I,IA_ERR,ICOMP,ICPACONV,ICPAFLAG,IDIMS,IE,IECPAFAIL(NEMAX)
     &        ,IECURR,IFMT,IKMCOR(:,:),IL,ILOP,IM,INC,IOBS,IPROC,IQ,IR,
     &        IT,ITCPA,ITERM,ITP,ITRCHI,ITRCHI0,ITXRAYDUM,IWRI,IWRIRRWF,
     &        IWRREGWF,IZERO(:),JTOP,KAPCOR(:),LABXCFIL,LFN,LL,LTAUFIL,
     &        MM05COR(:),NCHIITER,NCPAFAIL,NCSTMAX,NIN,NKPCOR(:),NQ9,
     &        NSPINOBS(NOBS),NSPINOP(NOP),NT9,NTKMAX
      CHARACTER*40 PATH(0:7)
      CHARACTER*15 STR15
C
C*** End of declarations rewritten by SPAG
C
      DATA PATH/'ONLY REAL ENERGIES   Gauss-Legendre    ',
     &     'ONLY REAL ENERGIES   Trapez rule       ',
     &     'RECTANGULAR COMPLEX PATH               ',
     &     'STRAIGHT CPLX. PATH PARA. TO REAL AXIS ',
     &     'RECT. CPLX. GRID, RETURN ON LOG SCALE  ',
     &     'ARC IN COMPLEX PLANE                   ',
     &     'STANDARD X-RAY E-MESH                  ',
     &     'X-RAY E-MESH for integration           '/
      DATA ITXRAYDUM/0/
C ======================================================================
C
C------------------------------------------- variables depending on NCST
      ALLOCATABLE FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR
      ALLOCATABLE NKPCOR
C
      ALLOCATABLE BT0,BCPT,EXC
C
C-----------------------------------------------------------------------
C  these arrays are needed only if information on a specific core shell
C   has to be obtained from <CORE> e.g. for core level spectroscopies
C-----------------------------------------------------------------------
C
      NCSTMAX = 1
C
      ALLOCATE (EXC(NRMAX))
      ALLOCATE (IZERO(NCSTMAX),SZCOR(NCSTMAX),ECOR(NCSTMAX))
      ALLOCATE (KAPCOR(NCSTMAX),MM05COR(NCSTMAX))
      ALLOCATE (IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX))
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: CHI -> GCOR'
C
      ALLOCATE (BT0(NRMAX,NTMAX),BCPT(NRMAX,NTMAX))
C
      CALL RVECCOP(NRMAX*NT,BT,BT0)
C
C-----------------------------------------------------------------------
C
      IWRI = 6
C
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      IWRREGWF = 1
      IWRIRRWF = 1
      WRTAU = .TRUE.
C
C-----------------------------------------------------------------------
C
      IF ( .NOT.(MPI) ) THEN
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .FALSE.
      ELSE
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .TRUE.
         WRITE (6,*) '###############################'
         WRITE (6,*) '########## <CHISUSC> ##########'
         WRITE (6,*) '##########    MPI    ##########'
         WRITE (6,*) '###############################'
         WRITE (6,*) '##### UNDER CONSTRUCTION  #####'
         WRITE (6,*) '###############################'
         WRITE (6,*) '###############################'
         STOP
      END IF
C
C-----------------------------------------------------------------------
C
      NSPINOP(1) = 1
      NSPINOP(2) = 1
      NSPINOP(3) = 2
C
      DO IOBS = 1,NOBS
         IF ( IOBS.LE.NOP ) THEN
            NSPINOBS(IOBS) = NSPINOP(IOBS)
         ELSE
            NSPINOBS(IOBS) = 1
         END IF
      END DO
C
      NTKMAX = LINMAX
C
      WRITE (6,99005)
C
C ======================================================================
C === orbital polarisation =============================================
C ======================================================================
C
      CHIPRINT = IPRINT
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_INTEGER('PRINT',CHIPRINT,9999,0)
C
         CALL SECTION_SET_STRING('OP',CHIORBPOL,'9999',0)
         IF ( .NOT.FOUND_STRING ) CHIORBPOL = 'BROOKS    '
         ILOP = 3
         IF ( ORBPOL(1:8).EQ.'BROOKS-F' ) ILOP = 4
      END IF
C
C=======================================================================
C
      IPRINT = 0
      WRITE (6,99006) CHIPRINT,CHIORBPOL
C
C ======================================= RESUME AXC-BXC-SCF-CYCLE ? ===
      ITRCHI0 = 0
      ITRCHI = 0
      ABXCFIL = DATSET(1:LDATSET)//'.abxc'
      LABXCFIL = LDATSET + 5
      RDABXC = .FALSE.
      CALL SECTION_FIND_KEYWORD('RDABXC',RDABXC)
      IF ( RDABXC ) THEN
         CALL SECTION_SET_INTEGER('ITRCHI0',ITRCHI0,9999,0)
         WRITE (6,99007) 'resuming',ABXCFIL(1:LABXCFIL),ITRCHI0
      ELSE
         WRITE (6,99007) 'saving',ABXCFIL(1:LABXCFIL)
      END IF
C
C ================================ READ VARIABLES FOR ABXC-SCF-CYCLE ===
C
      CALL SECTION_SET_REAL('MIX',CHIMIX,0.5D0,0)
C
      CALL SECTION_SET_REAL('TOL',CHITOL,1.0D-5,0)
C
      CALL SECTION_SET_INTEGER('NCHIITER',NCHIITER,100,0)
C
      IF ( NCHIITER.GT.1 ) THEN
         WRITE (6,99008) CHITOL,CHIMIX,NCHIITER
      ELSE
         WRITE (6,99010)
      END IF
C
C ======================================================================
C
C ---------------------------- CHECK FOR SURFACE OR SLAB CALCULATION ---
      VMTZOFF = 0.0D0
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( FOUND_SECTION ) THEN
         IF ( IBZINT.EQ.0 ) THEN
            CALL SECTION_FIND_KEYWORD('SURFACE',SURFACE)
            CALL SECTION_FIND_KEYWORD('SLAB',SLAB)
            IF ( SURFACE .OR. SLAB ) THEN
               CALL SECTION_SET_REAL('VMTZOFF',VMTZOFF,9999D0,0)
               IF ( FOUND_REAL ) THEN
                  WRITE (6,*) 'MUFFIN-TIN ZERO SHIFTED BY VMTZOFF=',
     &                        VMTZOFF,' Ry'
                  DO IT = 1,NT
                     DO IR = 1,NRMAX
                        VT(IR,IT) = VT(IR,IT) - VMTZOFF
                     END DO
                  END DO
               END IF
            END IF
         END IF
C
         NONMAG = .FALSE.
         CALL SECTION_FIND_KEYWORD('NONMAG',NONMAG)
      END IF
C
C ======================================================================
C                  DOS at the Fermi level
C ======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( RDTAU ) WRTAU = .FALSE.
C
      CALCNEF = .FALSE.
C
      DO IT = 1,NT
         STR10 = 'NEF'
         CALL STRING_ADD_N(STR10,IT)
C
         CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,NIN,NLMAX,0,9999D0,0)
C
         IF ( FOUND_REAL_ARRAY ) THEN
            DO IL = 1,NLMAX
               NEFTL(IT,IL) = AUXL(IL)
            END DO
         ELSE
            CALCNEF = .TRUE.
            WRITE (6,99011)
            EXIT
         END IF
      END DO
C
      IF ( CALCNEF ) THEN
C
         ERYD = DCMPLX(EFERMI,0.001D0)
         IECURR = 1
         ICPAFLAG = 0
C
         CALL SSITE(0,0,IFILCBWF,.TRUE.,.TRUE.,ERYD,P,IPRINT,NKM,TSST,
     &              MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
         CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                        ICPACONV,CONC,NOQ,ITOQ,PHASK,IECURR,NTMAX,
     &                        TSST,MSST,TSSQ,MSSQ,TAUQ)
C
         CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
         CALL CALCDOS(0,.FALSE.,.FALSE.,1,0,-1,ERYD,MEZZ,MEZJ,TSST,MSST,
     &                TAUT,MSSQ,TAUQ,IECURR,WETAB(1,1),BCOR,BCORS,
     &                DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX,DOBS_TX_GLO,
     &                OBS_TX_GLO,CTOTDOS)
C
         TOTDOS = 0.0D0
         WRITE (6,99012)
         DO IT = 1,NT
            DO IL = 1,NLT(IT)
               NEFTL(IT,IL) = DIMAG(DOBS_LTX(0,IDOS,IL,IT))
               TOTDOS = TOTDOS + NEFTL(IT,IL)*CONC(IT)*NAT(IT)
            END DO
C
            WRITE (6,99013) IT,(NEFTL(IT,IL),IL=1,NLT(IT))
         END DO
C
         IF ( ICPAFLAG.NE.0 ) STOP ' CPA - cycle not converged'
C
      END IF
C
C ======================================================================
      CALL SECTION_SET_REAL_ARRAY('SIGNAF',SIGNAF,NIN,NTMAX,0,9999D0,0)
      TREATAF = FOUND_REAL_ARRAY
C
      IF ( TREATAF ) THEN
         DO IT = 1,NT
            SIGNAF(IT) = SIGN(1.0D0,SIGNAF(IT))
         END DO
         WRITE (6,*) '>>>>>>>>>>>>>>>>>> SIGNAF',(SIGNAF(IT),IT=1,NT)
      END IF
C
C -----------------------------------------------------------------
C *** FORMAT OF TAU x TAU-FILE:
C    1)  0=FORMATTED   -> 1 file for tau and tautau formatted
C or 2)  1=UNFORMATTED -> 2 files: one for tau      formatted
C                                  one for tautau   UNformatted
C -----------------------------------------------------------------
C
      TAUFIL = DATSET(1:LDATSET)//'.tau'
      LTAUFIL = LDATSET + 4
C
      FOUND_SECTION = .FALSE.
      IF ( RDTAU ) THEN
         OPEN (9,FILE=TAUFIL(1:LTAUFIL),STATUS='old',FORM='formatted')
         REWIND 9
      END IF
C
      NEPATH = 1
      WRITE (6,99016) PATH(IGRID(1))
      WRITE (6,99014) 'energies  :  ',NETAB(1)
      WRITE (6,99015) 'E_min     :  ',EMIN
      WRITE (6,99015) 'E_max     :  ',EMAX
      IF ( IGRID(1).NE.5 ) WRITE (6,99015) 'Im(E)     :  ',EIMAG
      WRITE (6,99015) 'E_Fermi   :  ',EFERMI
C
      CALL EPATH(IGRID(1),EMIN,EMAX,EIMAG,NETAB(1),ETAB(1,1),WETAB(1,1),
     &           EILOW,IPRINT,NEMAX)
C
C ======================================================================
C                     get core charge density
C ======================================================================
C
      WRITE (IWRI,'(1X,79(''-''))')
      WRITE (IWRI,'(A,/)') 'running CORE to get core charge density'
C
      CALL RVECCOP(NRMAX*NT,RHOCHR,RHOCHRC)
      CALL RVECCOP(NRMAX*NT,RHOSPN,RHOSPNC)
C
      CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &          IKMCOR,IZERO,ITXRAYDUM,BCOR,BCORS,NCSTMAX)
      WRITE (IWRI,'(1X,79(''-''),//)')
C
      DO IT = 1,NT
         IM = IMT(IT)
         DO I = 1,JRWS(IM)
            SWP = RHOCHRC(I,IT)
            RHOCHRC(I,IT) = RHOCHR(I,IT)
            RHOCHR(I,IT) = SWP
            SWP = RHOSPNC(I,IT)
            RHOSPNC(I,IT) = RHOSPN(I,IT)
            RHOSPN(I,IT) = SWP
         END DO
      END DO
C
C ======================================================================
C                      run ssite for E_F
C ======================================================================
      ERYD = DCMPLX(EFERMI,0.0D0)
C
      CALL SSITE(IWRREGWF,IWRIRRWF,IFILCBWF,CALCINT,GETIRRSOL,ERYD,P,
     &           IPRINT,NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C ======================================================================
C                    calculate GAMMA(R,IT)
C ======================================================================
C
      CALL CHIGAMMA(IFILCBWF,GAMMA,GAMMAN,BXCNMM,BXCNNM,BXCNMN,BXCNNN,
     &              AXCN,NEFTL,ORBSQINT,EMM,ENM,ENN)
C
      IF ( RDABXC ) THEN
         OPEN (29,FILE=ABXCFIL(1:LABXCFIL),STATUS='OLD',
     &         FORM='FORMATTED')
         REWIND 29
         DO IT = 1,NT
            IM = IMT(IT)
            JTOP = JRWS(IM)
C
            READ (29,'(''IT = '',I2)',END=300,ERR=200) ITP
            IF ( ITP.NE.IT ) THEN
               WRITE (6,*) ' INCORRECT FORMAT OF AXC-BXC-FILE ',
     &                     ABXCFIL(1:LABXCFIL)
               WRITE (6,*) 'TRY WITHOUT RDABXC IN INPUT-FILE'
               STOP
            ELSE
               READ (29,'(16X,5(1E15.7,1X))',END=300,ERR=200)
     &               (BXCNP(I,IT),AXCNP(I,1,IT),AXCNP(I,1,IT),
     &               GAMMAP(I,IT),I=1,JTOP)
C
               DO I = 1,JTOP
                  BXCNMM(I,IT) = BXCNP(I,IT)
                  AXCN(I,1,ILOP,IT) = AXCNP(I,1,IT)
                  AXCN(I,2,ILOP,IT) = AXCNP(I,2,IT)
                  GAMMA(I,IT) = GAMMAP(I,IT)
               END DO
            END IF
         END DO
      ELSE
         OPEN (29,FILE=ABXCFIL(1:LABXCFIL),FORM='FORMATTED')
      END IF
C
C ======================================================================
C           run core for getting core polarization hyperfine field
C ======================================================================
C *** choose a small external field to stay in linear response regime
C *** e.g. 10^-4 Ryd/mu_B  ( see Diss. HF -> Eq. (4.59) )
C
      BEXTCP = 1D-4
C
      DO IT = 1,NT
         IM = IMT(IT)
         JTOP = JRWS(IM)
         DO I = 1,JTOP
            VAUX(I,1) = 0.0D0
            VAUX(I,2) = 0.0D0
            RHO4PI(I,1) = RHOCHR(I,IT)
            RHO4PI(I,2) = GAMMA(I,IT)*BEXTCP
            WEXC(I,1) = 0.0D0
            WEXC(I,2) = 0.0D0
         END DO
         CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,JTOP,NRMAX)
         DO I = 1,JTOP
            BCPT(I,IT) = (VAUX(I,1)-VAUX(I,2))/2.0D0
         END DO
      END DO
C
      WRITE (IWRI,99009)
C
      CALL RVECCOP(NRMAX*NT,RHOCHR,RHOCHRC)
      CALL RVECCOP(NRMAX*NT,RHOSPN,RHOSPNC)
      CALL RVECCOP(NRMAX*NT,BCPT,BT)
C
      CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &          IKMCOR,IZERO,ITXRAYDUM,BCP,BCPS,NCSTMAX)
C
      CALL RVECCOP(NRMAX*NT,RHOCHRC,RHOCHR)
      CALL RVECCOP(NRMAX*NT,RHOSPNC,RHOSPN)
      CALL RVECCOP(NRMAX*NT,BT0,BT)
C
      DO IT = 1,NT
         BCP(IT) = BCP(IT)/BEXTCP
         BCPS(IT) = BCPS(IT)/BEXTCP
      END DO
      WRITE (IWRI,'(1X,79(''-''),/)')
C
C ======================================================================
C                initialize CHI - iteration - cycle
C ======================================================================
C  ITRCHI0: initial value set to 0 or to ITRCHI0 from input file
C
      DO IT = 1,NT
         IM = IMT(IT)
         DO I = 1,JRWS(IM)
            BXCMNV(I,IT,1) = BXCNMM(I,IT)
            BXCMNV(I,IT,2) = BXCNNM(I,IT)
            BXCMNV(I,IT,3) = BXCNMN(I,IT)
            BXCMNV(I,IT,4) = BXCNNN(I,IT)
         END DO
      END DO
C
      IF ( NCHIITER.GT.1000 ) THEN
         IF ( .NOT.RDABXC ) THEN
            DO IT = 1,NT
C
               CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_ABXC.dat',
     &                       9,FILNAM,LFN,2,IOTMP,'ABXC-file ',10,NTMAX)
C
               WRITE (IOTMP,99019) IT,TXT_T(IT)(1:LTXT_T(IT)),NL
               WRITE (IOTMP,'(5(1E15.7,1X))')
     &                (R(I,IM),BXCNMM(I,IT),AXCN(I,1,ILOP,IT),
     &                AXCN(I,2,ILOP,IT),GAMMA(I,IT),I=1,JTOP)
               CLOSE (IOTMP)
            END DO
         END IF
C
      END IF
C
C=======================================================================
C                    XC - SCF - cycle  START
C=======================================================================
C
 100  CONTINUE
      ITRCHI = ITRCHI + 1
C
      IF ( RDTAU ) THEN
         REWIND 9
C
         READ (9,'(//,10X,I3,5X,I3,6X,I2)',END=500,ERR=400) NT9,NQ9,IFMT
         WRITE (IWRI,*) ' IFMT=',IFMT
         IF ( NQ9.NE.NQ .OR. NT9.GT.NT ) THEN
            WRITE (IWRI,*) 'TAU-FILE inconsistent '
            WRITE (IWRI,*) ' NQ NQ9 ',NQ,NQ9,' NT NT9 ',NT,NT9
            STOP 'STOP on reading the head of taufile!'
         END IF
         DO I = 1,(NT9+NQ9)
            READ (9,'(A80)',END=500,ERR=400) LINE
            WRITE (IWRI,'(A80)') LINE
         END DO
      END IF
C
      IF ( WRTAU .AND. (MPI_ID.EQ.0) ) THEN
         CALL INPUT_FIND_SECTION('ENERGY',0)
         OPEN (9,FILE=TAUFIL(1:LTAUFIL),FORM='formatted')
C
         REWIND 9
C
         WRITE (9,99018) CMDUC(3:LEN_TRIM(CMDUC)),NT,NQ
         WRITE (9,'(5X,'' IT ='',I3,'' IQ ='',I3,'': '',A)')
     &          (IT,IQAT(1,IT),TXT_T(IT),IT=1,NT)
         DO IQ = 1,NQ
            WRITE (9,'(''site'',I3,'' of '',A)') IQ,SYSTEM(1:LSYSTEM)
         END DO
C
      END IF
C
      NCPAFAIL = 0
C
      CALL CPU_TIME(TIME1)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IPROCE(1:NEMAX) = 0
C
      IF ( NPROCS.GT.1 ) THEN
         IPROC = 0
         INC = 1
         DO IE = 1,NETAB(1) - 1
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
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      WRITE (IWRI,99017) RDTAU,WRTAU
C
      DO IE = 1,NETAB(1)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCE(IE) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IECURR = IE
C
            ERYD = ETAB(IE,1)
C
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
            CALL CINIT(NTKTKMAX*NTMAX,DDTAUTAUT)
            CALL CINIT(NTKTKMAX*NTMAX*NTMAX,TKTKTT)
            CALL CINIT(NKMMAX*NKMMAX*NTMAX,TAUT)
C
C ===================================== solve SS - differential equation
C
            CALL SSITE(IWRREGWF,IWRIRRWF,IFILCBWF,CALCINT,GETIRRSOL,
     &                 ERYD,P,IPRINT,NKM,TSST,MSST,SSST,MEZZ,MEZJ,
     &                 ORBPOL)
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C                                                  START:  CALCULATE TAU
            IF ( .NOT.RDTAU ) THEN
C
               IF ( ITEST.EQ.3 ) CALL CHITKTKCHK(.TRUE.,ERYD,P,IPRINT,
     &              ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,MSST,TSSQ,MSSQ,
     &              TAUQ,TAUQZ)
C
               CALL CHIKLOOPS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,NQ,NKM,
     &                        NKKR,NKMQ,IND0Q,KMROT,DROTQ,CPATOL,NCPA,
     &                        ICPA,ICPAALG,ITCPA,ITCPAMAX,ICPACONV,DROT,
     &                        IQORGQP,SYMUNITARY,SYMACCEPTED,NOQ,ITOQ,
     &                        CONC,WKTAB,KTAB,NKTAB,NSYM,NSYMACCEPTED,
     &                        NTMAX,NQMAX,NKMMAX,NKTABMAX,TSST,MSST,
     &                        TSSQ,MSSQ,TAUQ,TAUQZ)
C
               IF ( ICPAFLAG.NE.0 ) THEN
                  NCPAFAIL = NCPAFAIL + 1
                  CPACHTAB(NCPAFAIL) = CPACHNG
                  IECPAFAIL(NCPAFAIL) = IE
               END IF
C
               CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
C                                                     end: projected TAU
C ======================================================================
C
               CALL CHITKTKQ2T(MSSQ,MSST,TAUQ,DMATT,DTILT)
C
C                                              end: TAUTAUQQ -> TAUTAUTT
C ======================================================================
C
            END IF
C                                                    END:  CALCULATE TAU
C=======================================================================
C
            CALL CHICALC(SIGNAF,TREATAF,ILOP,RPWOP,MEZZ,MEZJ,BCP,BCPS,
     &                   TAUT,TKTKTT,DDTAUTAUT,.FALSE.,CHILANLT,
     &                   RHORINTL,RHORTL1,EMM,ENM,ENN,GAMMA,GAMMAN,
     &                   BXCNMM,BXCNNM,BXCNMN,BXCNNN,AXCN,ORBSQINT,DOS1,
     &                   DOS1L,DOSL1,CHIPRINT,ERYD,IE,IECURR,NTKMAX,
     &                   NTKTKMAX,IDOS,ISPN,IORB,IHFI,IHVV,NSPINOP,
     &                   NSPINOBS)
C
C ======================================================================
            IF ( IPRINT.GE.3 ) CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,
     &           TAUQ)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      END DO
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      CALL CPU_TIME(TIME2)
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C ======================================================================
C                        print Landau susceptibility
C ======================================================================
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         IDIMS = (1+NLMAX)*NTMAX*2
C
         CALL DRV_MPI_REDUCE_R(CHILANLT(0,1,1),RWORK(0,1,1),IDIMS)
C
         IF ( MPI_ID.EQ.0 .AND. TASK.EQ.'CHILANDAU ' ) THEN
C
            CHILAN = 0D0
            DO ITERM = 1,2
               DO IT = 1,NT
                  CHILANLT(0,IT,ITERM) = 0D0
C
                  DO IL = 1,NLT(IT)
                     CHILANLT(0,IT,ITERM) = CHILANLT(0,IT,ITERM)
     &                  + CHILANLT(IL,IT,ITERM)
                  END DO
C
                  CHILAN = CHILAN + NAT(IT)*CONC(IT)
     &                     *CHILANLT(0,IT,ITERM)
               END DO
C
               WRITE (6,99001) (TXT_L(IL),IL=1,NL)
               DO IT = 1,NT
                  WRITE (6,99002) IT,TXT_T(IT),
     &                            (CHILANLT(IL,IT,ITERM)*CU,IL=0,NLT(IT)
     &                            )
               END DO
               IF ( NT.GT.1 ) WRITE (6,99003) CHILAN*CU
               WRITE (6,99004)
            END DO
C
         END IF
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      WRITE (IWRI,99021) TIME2 - TIME1
C
      IF ( NCPAFAIL.NE.0 ) THEN
         WRITE (IWRI,99022) CPATOL,NCPAFAIL,
     &                      (IECPAFAIL(I),DREAL(ETAB(IECPAFAIL(I),1)),
     &                      CPACHTAB(I),I=1,NCPAFAIL)
         WRITE (IWRI,'(1X,79(''*''),/)')
      ELSE IF ( NCPA.NE.0 ) THEN
         WRITE (IWRI,99023)
      END IF
C
C ======================================================================
C                 CHI - iteration - cycle
C ======================================================================
C
C ======================================================================
C---------------------- suppress the iteration of the spin magnetisation
      NCHIITER = 1
C ======================================================================
C
      IF ( NCHIITER.GT.1 ) THEN
C
C----------------------------------- mix AXCN/BXCN and check convergence
C
         CALL RINIT(NRMAX*NTMAX,BXCNDUM)
C
         DO IT = 1,NT
            IM = IMT(IT)
            DO I = 1,JRWS(IM)
               BXCMNV(I,IT,1) = BXCNMM(I,IT)
               BXCMNV(I,IT,2) = BXCNNM(I,IT)
               BXCMNV(I,IT,3) = BXCNMN(I,IT)
               BXCMNV(I,IT,4) = BXCNNN(I,IT)
            END DO
         END DO
C
         DO ICOMP = 1,4
C
            DO IT = 1,NT
               CALL RINIT(NRMAX*NTMAX,BXCNDUM)
            END DO
            CALL SCFTCHEBY(ITRCHI,BXCMNV(1,1,ICOMP),BXCNDUM,CHIMIX,
     &                     RMSMN(1,ICOMP),RMSMN(2,ICOMP))
C
         END DO
C
         DO IT = 1,NT
            IM = IMT(IT)
C ******* TAKE CARE OF NUMERICAL INACCURACY FOR SMALL R ***
            INSIDE = .TRUE.
            BXCCORR = .FALSE.
            DO I = 1,JRWS(IM)
               IF ( INSIDE .AND. (BXCNMM(I,IT).LT.0.0D0) ) THEN
                  BXCNMM(I,IT) = 0.0D0
                  BXCCORR = .TRUE.
               ELSE
                  INSIDE = .FALSE.
               END IF
            END DO
            IF ( BXCCORR ) WRITE (6,*) 
     &                       'WARNING IN CHI !!! BXCN CORRECTED FOR IT='
     &                       ,IT
         END DO
C
         REWIND 29
         STR15 = '_ABXC'
         CALL STRING_ADD_N(STR15,ITRCHI0+ITRCHI)
         LL = LEN_TRIM(STR15)
         STR15 = STR15(1:LL)//'.dat'
         LL = LL + 4
C
         DO IT = 1,NT
            IM = IMT(IT)
            CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,STR15(1:LL),LL,
     &                    FILNAM,LFN,2,IOTMP,'ABXC-file ',10,NTMAX)
C
            WRITE (IOTMP,99020) IT,TXT_T(IT)(1:LTXT_T(IT)),NL
            WRITE (IOTMP,'(5(1E15.7,1X))')
     &             (R(I,IM),BXCNMM(I,IT),AXCN(I,1,ILOP,IT),
     &             AXCN(I,2,ILOP,IT),GAMMA(I,IT),I=1,JTOP)
            CLOSE (IOTMP)
         END DO
C
         WRITE (6,'(/,79(''-''))')
         WRITE (6,*) 'CHI - iteration - cycle : ITRCHI = ',
     &               ITRCHI0 + ITRCHI
         WRITE (6,*) 'ERR: BXCMM = ',RMSMN(1,1),RMSMN(2,1)
         WRITE (6,*) 'ERR: BXCNM = ',RMSMN(1,2),RMSMN(2,2)
         WRITE (6,*) 'ERR: BXCMN = ',RMSMN(1,3),RMSMN(2,3)
         WRITE (6,*) 'ERR: BXCNN = ',RMSMN(1,4),RMSMN(2,4)
C
         IF ( DSQRT(RMSMN(1,1)**2+RMSMN(1,2)**2+RMSMN(1,3)**2+RMSMN(1,4)
     &        **2).LE.CHITOL ) THEN
C
            WRITE (6,*) 'CHI - iteration - cycle CONVERGED !!!'
C
         ELSE IF ( ITRCHI.LT.NCHIITER ) THEN
C
            WRITE (6,*) 'CHI - iteration - cycle NOT CONVERGED ...'
            WRITE (6,*) 'NEW ITER NEW ITER NEW ITER NEW ITER NEW ITER'
            WRITE (6,*) '... STARTING NEXT ITERATION',ITRCHI0 + ITRCHI + 
     &                  1
            WRITE (6,'(79(''-''),/)')
C
C *** REOPEN FILE 9 ***
            OPEN (9,FILE=TAUFIL(1:LTAUFIL),STATUS='OLD',
     &            FORM='FORMATTED')
C *** GOTO NEXT ITERATION ***
            GOTO 100
         ELSE
            WRITE (6,*) 'CHI - iteration - cycle : NCHIITER EXHAUSTED'
         END IF
C
      END IF
C                                end mix AXCN/BXCN and check convergence
C ======================================================================
C
      WRITE (6,*) '          <CHISUSZ> - run completed'
      STOP
C
 200  CONTINUE
      WRITE (IWRI,*) ' stop in <CHISUSZ>'
      WRITE (IWRI,*) ' ERROR WHILE READING FILE 29'
      WRITE (IWRI,*) ' IT, ITP = ',IT,ITP
      STOP
 300  CONTINUE
      WRITE (IWRI,*) ' stop in <CHISUSZ>'
      WRITE (IWRI,*) ' END OF FILE REACHED WHILE READING FILE 29 '
      STOP
 400  CONTINUE
      WRITE (IWRI,*) ' stop in <CHISUSZ>'
      WRITE (IWRI,*) ' ERROR WHILE READING FILE 9: IECURR= ',IECURR,
     &               ' (IE= ',IE,')'
      STOP
 500  CONTINUE
      WRITE (IWRI,*) ' stop in <CHISUSZ>'
      WRITE (IWRI,*) ' END OF FILE 9 REACHED FOR IECURR= ',IECURR,
     &               ' (IE= ',IE,')'
      STOP
C
99001 FORMAT (/,1X,79('='),//,10X,
     &        'Landau magnetic susceptibility  in [10^(-6) cm^3/mol]',
     &        //,1X,79('='),//,10X,'type',9X,'sum         ',5(A,:,10X))
99002 FORMAT (10X,I2,1X,A,6F20.10)
99003 FORMAT (10X,'total  ',6F11.4)
99004 FORMAT (/,1X,79('='),/)
99005 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*     ****   *    *  ***   ****   *    *   ****    ****      *'
     &  ,/,10X,
     &  '*    *    *  *    *   *   *    *  *    *  *    *  *    *     *'
     &  ,/,10X,
     &  '*    *       *    *   *   *       *    *  *       *          *'
     &  ,/,10X,
     &  '*    *       ******   *    ****   *    *   ****   *          *'
     &  ,/,10X,
     &  '*    *       *    *   *        *  *    *       *  *          *'
     &  ,/,10X,
     &  '*    *    *  *    *   *   *    *  *    *  *    *  *    *     *'
     &  ,/,10X,
     &  '*     ****   *    *  ***   ****    ****    ****    ****      *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99006 FORMAT (1X,79('*'),/,24X,'settings for CHI - calculation ',/,1X,
     &        79('*'),//,5X,'print level            ',I3,/,5X,
     &        'orbital polarisation   ',A)
99007 FORMAT (5X,A,' AXC-BXC-SCF-cycle ',/,5X,
     &        'using the file           ',A,:,/,5X,
     &        'last iteration           ',I12)
99008 FORMAT (5X,'BXC/AXC-SCF-tolerance    ',F12.6,/,5X,
     &        'BXC/AXC-mixing-parameter ',F12.6,/,5X,
     &        'max. number of iterations',I12)
99009 FORMAT (1X,79('-'),/,5X,'running CORE to get ',
     &        'core polarisation hyperfine field')
99010 FORMAT (/,5X,'no iteration of induced normalized densities',/)
99011 FORMAT (/,5X,'DOS  n(E_F)  not supplied in input file  ....',
     &        '  will be calculated ',/)
99012 FORMAT (/,5X,'l-resolved DOS  n(E_F) [1/Ry] at the Fermi level',/)
99013 FORMAT (5X,'IT = ',I2,2X,10(F14.4))
99014 FORMAT (10X,A,5I10)
99015 FORMAT (10X,A,5F10.6)
99016 FORMAT (/,5X,'energy path:  ',A)
99017 FORMAT (/,10X,'RDTAU     =',L2,'  WRTAU     =',L2)
99018 FORMAT ('ENERGY ',A,/,1X,79('*'),/,5X,' NT =',I3,' NQ =',I3,
     &        ' FMT = 3')
99019 FORMAT ('# spin and orbital exchange correlation fields',/,
     &        '# for type  IT=',I2,4X,A,'   and   NL=',I3,/,
     &        '# initial values ',/,
     &        '#   R               BXC             AXC(1)          ',
     &        'AXC(2)          GAMMA')
99020 FORMAT ('# spin and orbital exchange correlation fields',/,
     &        '# for type  IT=',I2,4X,A,'   and   NL=',I3,/,
     &        '#   R               BXC             AXC(1)          ',
     &        'AXC(2)          GAMMA')
99021 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99022 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99023 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
      END
