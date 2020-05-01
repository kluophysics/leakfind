C*==negfdrive.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE NEGFDRIVE
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to run calculations in the                           *
C   *  NON-equilibrium Green's function mode                           *
C   *                                                                  *
C   *  NOTE:                                                           *
C   *      *  surface Green's function approach by M. Ogura is used    *
C   *      *  the Green's function convention is   RH                  *
C   *      *  the conductance is calculated using the                  *
C   *         (L,ML,MS) representation - even for IREL = 3             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,WKM1,WKM2,WKM3,MEZJ,MEZZ,
     &    TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_CALCMODE,ONLY:ORBPOL,THERMAL_VIBRA_FLUCT,DMFT,LDAU,
     &    GF_CONV_RH,SOLVER_FP
      USE MOD_THERMAL,ONLY:I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT
      USE MOD_CONSTANTS,ONLY:PI,CI,C0,C1
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,NCPA
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      USE MOD_ENERGY,ONLY:ETAB,WETAB,NEPATH,NETAB,EFERMI,NEMAX,EIMAG
      USE MOD_FILES,ONLY:IFILCBWF,RECLNGWF,IPRINT,WRTAUMQ,IFILTAU,
     &    TAUFIL,IFILBUILDBOT,WRBUILDBOT
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,NPROCS,MPI_ID,MPI_ELOOP
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI,NQTB,IQBOT_TB,
     &    IQTOP_TB
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY
      USE MOD_TB,ONLY:NSLAY_PER_PLAY
      USE MOD_TYPES,ONLY:NT,NTMAX,DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX,
     &    DOBS_TX_GLO,OBS_TX_GLO
      IMPLICIT NONE
C*--NEGFDRIVE36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CHECK_DOS
      PARAMETER (CHECK_DOS=.FALSE.)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NEGFDRIVE')
C
C Local variables
C
      LOGICAL BAR_FLAG_Q(:),CALCINT,CCZ,CHECK,CHECK_ZPM,DOS,GETIRRSOL,
     &        LBAR_FLAG_Q(:),RBAR_FLAG_Q(:),STT,TRANSMISSION
      REAL*8 BCOR(:),BCORS(:),CPACHNG,CPACHNGMAX,DMFTMIX,ERYDTOP,
     &       MEDIFF(:,:,:),TIME,TIME0
      COMPLEX*16 CDOS,CTOTDOS(:),ERYD,ERYD9,GLESQ(:,:,:),LMAT3(NKM,NKM),
     &           LMEDIFF(:,:),LMEZZMP(:,:),
     &           MEZZMM(NKMMAX,NKMMAX,NTMAX,NMEMAX),MEZZMP(:,:,:,:),P,
     &           PMIN,PPLS,ZMIN,ZPLS
      COMPLEX*16 CMATTRC
      INTEGER IA_ERR,ICPACONV,ICPAFLAG,IE,IEPATH,IESORT(:),IFIL_ZMIN,
     &        IFIL_ZPLS,INC,IPROC,IPROCE(:),IQ,IQBOT_LBAR,IQBOT_RBAR,
     &        IQTOP_LBAR,IQTOP_RBAR,ITCPA,IWRI,IWRIRRWF,IWRREGWF,JE,M,
     &        NCPAFAIL,NEB,NETAU,NQBAR,OUT
C
C*** End of declarations rewritten by SPAG
C
      DATA OUT/6/,CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/
      DATA ERYD9/(999999D0,999999D0)/
C
      ALLOCATABLE MEZZMP,GLESQ,MEDIFF,LMEZZMP,LMEDIFF
      ALLOCATABLE IPROCE,IESORT
      ALLOCATABLE BAR_FLAG_Q,LBAR_FLAG_Q,RBAR_FLAG_Q
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      ALLOCATABLE BCOR,BCORS,CTOTDOS
C
      ALLOCATE (MEDIFF(NKMMAX,NKMMAX,2),LMEZZMP(NKMMAX,NKMMAX))
      ALLOCATE (LMEDIFF(NKMMAX,NKMMAX))
      ALLOCATE (BCOR(NTMAX),BCORS(NTMAX),GLESQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (CTOTDOS(NEMAX))
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      M = NKMMAX
C
      ALLOCATE (MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),STAT=IA_ERR)
      ALLOCATE (BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ))
C
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M
      IF ( IA_ERR.NE.0 ) STOP 'alloc:<NEGFDRIVE_MO> -> MAQAB'
C
C      CHECK = .TRUE.
      CHECK = .FALSE.
C
      IQTOP_TB = IQBOT_TB + NQTB - 1
CSW   IQTOP_TB coming from mod_sites is incorrect! (NQ + 1)
C
C ======================================================================
C     FP - mode: default settings are overwritten here
C                set back to ZJ convention
C
Csw      GF_CONV_RH = .FALSE.
      GF_CONV_RH = .TRUE.
      SOLVER_FP = 'BS        '
C
      IF ( GF_CONV_RH ) THEN
         WRITE (6,99010)
      ELSE
         WRITE (6,99011)
      END IF
C
C ======================================================================
      WRITE (6,99003) 'NEGFDRIVE_MO'
C ======================================================================
C
      IF ( NO_SYMMETRY ) WRITE (6,*) 'NO_SYMMETRY'
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      TEMP_LAT = 0D0
      N_TEMP_LAT = 1
C
      CALL THERMAL_INIT(0,N_TEMP_LAT,TEMP_LAT)
C
C=======================================================================
C                  set barrier specific parameters
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
C     CALL SECTION_SET_INTEGER('NQBAR',NQBAR,9999,0)
C     IF ( .NOT.FOUND_INTEGER ) NQBAR = 6
      NQBAR = NSLAY_PER_PLAY
C
      IF ( NQBAR.GT.NQTB/2 ) STOP '<NEGFDRIVE_MO>: NQBAR.GT.NQTB/2'
C
      IQBOT_LBAR = IQBOT_TB
      IQTOP_LBAR = IQBOT_LBAR + NQBAR - 1
      IQTOP_RBAR = IQTOP_TB
      IQBOT_RBAR = IQTOP_RBAR - NQBAR + 1
C
      BAR_FLAG_Q(1:NQ) = .FALSE.
      LBAR_FLAG_Q(1:NQ) = .FALSE.
      RBAR_FLAG_Q(1:NQ) = .FALSE.
      DO IQ = IQBOT_LBAR,IQTOP_LBAR
         BAR_FLAG_Q(IQ) = .TRUE.
         LBAR_FLAG_Q(IQ) = .TRUE.
      END DO
      DO IQ = IQBOT_RBAR,IQTOP_RBAR
         BAR_FLAG_Q(IQ) = .TRUE.
         RBAR_FLAG_Q(IQ) = .TRUE.
      END DO
C
      WRITE (6,99002) NQBAR,IQBOT_LBAR,IQTOP_LBAR,IQBOT_RBAR,IQTOP_RBAR
C
      CALL SECTION_FIND_KEYWORD('CHECK_ZPM',CHECK_ZPM)
      IF ( CHECK_ZPM ) WRITE (6,*) 
     & 'SWSWSW    writing ssite & ms quantities to check z/z*    SWSWSW'
      CALL SECTION_FIND_KEYWORD('CCZ',CCZ)
      IF ( CCZ ) WRITE (6,*) 
     & 'SWSWSWSWSWSWSWSWSWSW     ZPLS <-> ZMIN     SWSWSWSWSWSWSWSWSWSW'
C
      CALL SECTION_FIND_KEYWORD('DOS',DOS)
      CALL SECTION_FIND_KEYWORD('STT',STT)
      CALL SECTION_FIND_KEYWORD('TRANSMISSION',TRANSMISSION)
      IF ( DOS .AND. TRANSMISSION ) STOP 
     &     'Currently NEGF DOS and TRANSMISSION are mutually exclusive'
      IF ( DOS .AND. STT ) STOP 
     &               'Currently NEGF DOS and STT are mutually exclusive'
      IF ( .NOT.STT .AND. .NOT.TRANSMISSION ) DOS = .TRUE.
C
      IF ( DOS ) WRITE (6,*) 
     & 'SWSWSWSW    Calculating the local DOS via G< and G+    SWSWSWSW'
      IF ( STT ) WRITE (6,*) 
     & 'SW      Calculating the Spin Transfer Torque (E_F,L->R)      SW'
      IF ( TRANSMISSION ) WRITE (6,*) 
     & 'SW    Calculating the Transmission from L to R: T(E,k_II)    SW'
      IF ( STT ) TRANSMISSION = .TRUE.
C
C=======================================================================
C
      M = NKMMAX
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( DMFT .OR. LDAU ) CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),
     &     DMFTMIX,NEMAX)
C
C ======================================================================
      WRITE (6,99004)
C ======================================================================
C
      CALL CPU_TIME(TIME0)
C
      IF ( IBZINT.NE.2 .AND. IBZINT.NE.6 .AND. IBZINT.NE.4 ) THEN
         WRITE (6,'("<NEGFDRIVE_MO>:  IBZINT   <>   2 | 4 | 6")')
         STOP
      END IF
C
      IF ( IBZINT.EQ.4 ) THEN
C
         ERYDTOP = 1.4D0*EFERMI
C
         CALL FRELPOLE(ERYDTOP,IPRINT)
      END IF
C
      NCPAFAIL = 0
      ITCPAMAX = 30
C
      CALL READTAU(9,ERYD9,0,NETAU,TAUT,0,NT,.FALSE.,TAUQ(1,1,1),
     &             MSSQ(1,1,1),1,NQ,0,NKM,NKMMAX,IPRINT)
C
      NEB = NETAB(1)
C
      IF ( .NOT.(MPI) ) THEN
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .FALSE.
      ELSE IF ( NEB.EQ.1 ) THEN
         MPI_KLOOP = .TRUE.
         MPI_ELOOP = .FALSE.
      ELSE
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .TRUE.
      END IF
C
      IFIL_ZPLS = IFILCBWF
      IFIL_ZMIN = IFILCBWF + 1
C
C     IFIL_ZPLS = IFILCBWF is already open
      OPEN (IFIL_ZMIN,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%    TAUQ       MSSQ                                   BLOCK BEGIN
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=======================================================================
C   calculate first TAUQ and MSSQ and store them for E-loops B and A
C   NB:   for sequential run (MPI_ELOOP=.FALSE.):
C=======================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      IWRI = 100
C
      IF ( NCPA.NE.0 .AND. NETAU.EQ.0 .AND. N_TEMP_LAT.EQ.1 ) THEN
C
         CALL THERMAL_INIT(1,N_TEMP_LAT,TEMP_LAT)
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop                               START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
         NEPATH = 1
         IEPATH = 1
         CPACHNGMAX = 0D0
C
         DO IE = 1,NEB
C
            IF ( NEB.EQ.1 ) THEN
               ERYD = EFERMI + CI*EIMAG
CSW     if EMIN isn''t supplied in input this is set to (-0.2,Im(E))!
            ELSE
               ERYD = ETAB(IE,IEPATH)
            END IF
C
            WRITE (6,99008) EFERMI,IE,ERYD
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
C ===================================== solve SS - differential equation
C
            IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &           = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_ZPLS,GETIRRSOL,
     &                    ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
C ------------- for tau-matrix calculation,  here we make it temporarily
C ------------- parallel in k (even if we are parallel in E in sigma
C ------------  calculation)
            IF ( MPI_ELOOP ) MPI_KLOOP = .TRUE.
C
            CALL TAU_TB(ERYD,P,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,
     &                  MSST,TSSQ,MSSQ,TAUQ)
C
            IF ( MPI_ELOOP ) MPI_KLOOP = .FALSE.
C
            IF ( ICPAFLAG.NE.0 ) THEN
               NCPAFAIL = NCPAFAIL + 1
               CPACHNGMAX = MAX(CPACHNGMAX,ABS(CPACHNG))
            END IF
C
C------------------------ to avoid asynchronisations using NFS sytems
C------------------------ we use this somewhat awkward scheme
            IF ( MPI ) THEN
               CALL DRV_MPI_BARRIER
               IF ( MPI_ID.EQ.0 ) THEN
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
                  CALL FLUSH(IFILTAU)
               ELSE
                  CLOSE (IFILTAU)
               END IF
               CALL DRV_MPI_BARRIER
               IF ( MPI_ID.NE.0 ) OPEN (UNIT=IFILTAU,FILE=TAUFIL)
               CALL DRV_MPI_BARRIER
            ELSE
               CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
            END IF
C
C ======================================================================
C
            IF ( IPRINT.GE.3 .OR. CHECK )
     &           CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,TAUQ)
C
C-----------------------------------------------------------------------
            IF ( IE.EQ.NETAB(IEPATH) ) THEN
C
               IF ( NCPAFAIL.NE.0 ) THEN
                  WRITE (OUT,99005) NCPAFAIL,CPATOL,CPACHNGMAX
               ELSE IF ( NCPA.NE.0 ) THEN
                  WRITE (OUT,99006)
               END IF
C
            END IF
C
            IF ( MPI ) CALL DRV_MPI_BARRIER
C
         END DO ! NEB
C
         CALL READTAU(9,ERYD9,0,NETAU,TAUT,0,NT,.FALSE.,TAUQ(1,1,1),
     &                MSSQ(1,1,1),1,NQ,0,NKM,NKMMAX,IPRINT)
      END IF
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop                                 END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      E-Punkte IE = 1...,(NEB-1) fuer E-Schleife
C       auf Prozessoren verteilen wie in SCF
C      IE=NEB auf IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCE(NEB),IESORT(NEB))
      DO IE = 1,NEB
         IESORT(IE) = IE
         IPROCE(IE) = 0
      END DO
C
      IF ( MPI_ELOOP ) THEN
         IPROC = 0
         INC = 1
         DO JE = 1,NEB - 1
            IE = IESORT(JE)
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
C
         CALL DRV_MPI_BARRIER
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C=======================================================================
C=======================================================================
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== START
C
CDK ToDo for Fermi sea calculation one needs to properly take care of
C        sequential run -- vibrations haveto be taken into account
C        --------------- sanity check -----------
      IF ( THERMAL_VIBRA_FLUCT .AND. MPI_ELOOP ) THEN
         WRITE (6,99009)
         STOP
      END IF
C
      DO I_TEMP_LAT = 1,N_TEMP_LAT
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop 2   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         WRTAUMQ = .FALSE.
C
         NEPATH = 1
         IEPATH = 1
C
         ICPAFLAG = 0
         CPACHNG = 0.0D0
C
C ======================================================================
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C
         CALL THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
C_EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         DO IE = 1,NEB
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ID.EQ.IPROCE(IE) .OR. MPI_KLOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
               IF ( NEB.EQ.1 ) THEN
                  ZPLS = EFERMI + CI*EIMAG
CSW     if EMIN isn''t supplied in input this is set to (-0.2,Im(E))!
               ELSE
                  ZPLS = ETAB(IE,IEPATH)
               END IF
C
CSW testing symmetry of tau(+/-) in negfkloop_gles!!
               IF ( CCZ ) ZPLS = DCONJG(ZPLS)
CSW testing symmetry of tau(+/-) in negfkloop_gles!!
               ZMIN = DCONJG(ZPLS)
C
               WRITE (6,99008) EFERMI,IE,ZPLS,ZMIN
C
C ===================================== solve SS - differential equation
C
               IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &              = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
               IF ( DMFT .OR. LDAU ) STOP 'DMFT to be done'
C
C ======================================================================
C
C------------------------------------------------------------------ Z(-)
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_ZMIN,
     &                       GETIRRSOL,ZMIN,PMIN,IPRINT,TSST,MSST,SSST,
     &                       MEZZMM,MEZJ,ORBPOL)
C
               IF ( CHECK_ZPM ) THEN
                  DO IQ = 1,NQ
                     CALL CMATSTRUCT('TSST_MIN',TSST(:,:,IQ),18,18,3,3,
     &                               1,1.D-8,557)
                  END DO
               END IF
C
C------------------------------------------------------------------ Z(+)
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL_ZPLS,
     &                       GETIRRSOL,ZPLS,PPLS,IPRINT,TSST,MSST,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
               IF ( CHECK_ZPM ) THEN
                  DO IQ = 1,NQ
                     CALL CMATSTRUCT('TSST_PLS',TSST(:,:,IQ),18,18,3,3,
     &                               1,1.D-8,558)
                  END DO
               END IF
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
               IF ( THERMAL_VIBRA_FLUCT ) CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
               CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
               IF ( MPI_ELOOP .OR. THERMAL_VIBRA_FLUCT ) THEN
C
                  CALL TAU_TB(ZPLS,PPLS,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,
     &                        TSST,MSST,TSSQ,MSSQ,TAUQ)
C
               ELSE IF ( NCPA.NE.0 ) THEN
C
                  DO IQ = IQBOT_CHI,IQTOP_CHI
                     CALL READTAU(9,ERYD9,IE,NETAU,TAUT,0,NT,.FALSE.,
     &                            TAUQ(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,1,NKM,
     &                            NKMMAX,IPRINT)
                     IF ( ABS(ZPLS-ERYD9).GT.1D-8 ) STOP 
     &                    'in <NEGFDRIVE_MO>: >>> FILE 9 ERYDB <> ERYD9'
                  END DO
C
               END IF
C
               IF ( IPRINT.GE.1 ) THEN
                  DO IQ = 1,NQ
                     WRITE (44,*) 'IQ = ',IQ
                     CALL CMATSTRUCT('TSST',TSST(:,:,IQ),18,18,3,3,1,
     &                               1.D-8,44)
                  END DO
               END IF
C ======================================================================
C
               CALL NEGFKLOOPSDRV_TB(IFIL_ZPLS,IFIL_ZMIN,ZPLS,PPLS,MSSQ,
     &                               TAUQ,GLESQ,IQBOT_LBAR,IQTOP_LBAR,
     &                               IQBOT_RBAR,IQTOP_RBAR,BAR_FLAG_Q,
     &                               LBAR_FLAG_Q,RBAR_FLAG_Q,MEZZ,
     &                               MEZZMM,MEZZMP,CCZ,CHECK_ZPM,LMAT3,
     &                               STT,TRANSMISSION)
C
               IF ( IPRINT.GE.1 ) CALL CMATSTRUCT('MSSQ',MSSQ(:,:,2),18,
     &              18,3,3,1,1.D-8,66)
C
C ======================================================================
C                         project on types
C ======================================================================
C
               CALL PROJTAU(ICPAFLAG,CPACHNG,ZPLS,MSST,MSSQ,TAUQ,TAUT)
C
C***********************************************************************
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
               IF ( .NOT.MPI_KLOOP .OR. MPI_ID.EQ.0 ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                  IF ( IPRINT.GE.1 ) CALL CMATSTRUCT('MSSQ',MSSQ(:,:,2),
     &                 18,18,3,3,1,1.D-8,66)
C
C ======================================================================
C
               END IF    ! MPI
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            END IF ! MPI
C
            IF ( DOS .OR. STT ) THEN
               OPEN (10,FILE='coupled.dos')
C
               IF ( IPRINT.GE.1 ) THEN
                  WRITE (6,*) 'epls',ZPLS
                  CALL CMATSTRUCT('mezz',MEZZ(:,:,1,1),18,18,3,3,1,
     &                            1.D-8,95)
                  CALL CMATSTRUCT('mezj',MEZJ(:,:,1,1),18,18,3,3,1,
     &                            1.D-8,95)
                  CALL CMATSTRUCT('taut',TAUT(:,:,2),18,18,3,3,1,1.D-8,
     &                            95)
                  CALL CMATSTRUCT('tsst',TSST(:,:,2),18,18,3,3,1,1.D-8,
     &                            95)
                  CALL CMATSTRUCT('msst',MSST(:,:,2),18,18,3,3,1,1.D-8,
     &                            95)
               END IF
C
               CALL CALCDOS(0,.FALSE.,.FALSE.,1,NCPAFAIL,IPRINT,ZPLS,
     &                      MEZZ,MEZJ,TSST,MSST,TAUT,MSSQ,TAUQ,IE,
     &                      WETAB(1,1),BCOR,BCORS,DOBS_LTX,DOBS_TX,
     &                      OBS_LTX,OBS_TX,DOBS_TX_GLO,OBS_TX_GLO,
     &                      CTOTDOS)
C
               CLOSE (10)
            END IF
C
C=======================================================================
C                       use GLESQ to calculate DOS
C=======================================================================
            IF ( MPI_ID.EQ.0 .AND. .NOT.TRANSMISSION .OR. STT ) THEN
C
               WRITE (6,*) ''
               WRITE (6,FMT='(A,3X,2F10.6)') 'ZPLS =',ZPLS
C
               IF ( IE.EQ.1 ) WRITE (123,FMT=
     &                               '(A,1X,A,2X,A,3X,A,7X,A,12X,A)')
     &                                '#','Re(ZPLS)','Im(ZPLS)','IQ',
     &                               'Re(CDOS)','Im(CDOS)'
C
               LOOP_GLES_DOS:DO IQ = 1,NQ
C
                  IF ( IQ.GE.IQBOT_LBAR .AND. IQ.LE.IQTOP_RBAR ) THEN
C
                     CALL ZGEMM('N','N',NKM,NKM,NKM,C1,MEZZMP(1,1,IQ,1),
     &                          NKMMAX,GLESQ(1,1,IQ),NKMMAX,C0,WKM1,
     &                          NKMMAX)
C
                     CDOS = -CMATTRC(NKM,NKMMAX,WKM1)/(2*PI)
C
                     WRITE (6,99001) IQ,CDOS
                     WRITE (123,FMT='(2F10.6,2X,I2,2F20.10)') ZPLS,IQ,
     &                      CDOS
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                     IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99012)
     &                    ROUTINE(1:LEN_TRIM(ROUTINE)),ZPLS,IQ,
     &                    DREAL(CDOS),DIMAG(CDOS)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
                     IF ( CHECK_DOS .AND. IE.EQ.1 ) THEN
                        CALL ZGEMM('N','N',NKM,NKM,NKM,C1,LMAT3,NKMMAX,
     &                             MEZZMP(1,1,IQ,1),NKMMAX,C0,LMEZZMP,
     &                             NKMMAX)
C
                        WKM1(:,:) = SQRT((DREAL(MEZZMP(:,:,IQ,1)))**2+(
     &                              DIMAG(MEZZMP(:,:,IQ,1)))**2)
                        WKM2(:,:) = SQRT((DREAL(LMEZZMP(:,:)))**2+(DIMAG
     &                              (LMEZZMP(:,:)))**2)
                        WKM3(:,:) = SQRT((DREAL(MEZZ(:,:,IQ,1)))**2+(
     &                              DIMAG(MEZZ(:,:,IQ,1)))**2)
C
                        MEDIFF(:,:,1) = DREAL(WKM1-WKM3)
                        MEDIFF(:,:,2) = DREAL(WKM2-WKM3)
                        LMEDIFF(:,:) = MEZZMP(:,:,IQ,1) - LMEZZMP(:,:)
C
                        WRITE (999,*) 'IQ =',IQ
                        WRITE (9991,*) 'IQ =',IQ
                        WRITE (9992,*) 'IQ =',IQ
                        WRITE (9993,*) 'IQ =',IQ
                        CALL CMATSTRUCT('MEZZ(IQ,1)',MEZZ(:,:,IQ,1),18,
     &                                  18,3,3,1,1.D-8,9991)
                        CALL CMATSTRUCT('MEZZMP(IQ,1)',MEZZMP(:,:,IQ,1),
     &                                  18,18,3,3,1,1.D-8,9991)
                        CALL RMATSTRUCT('MEDIFF(IQ,1)',MEDIFF(:,:,1),18,
     &                                  18,3,3,1,1.D-8,9991)
                        CALL CMATSTRUCT('MEZZ(IQ,1)',MEZZ(:,:,IQ,1),18,
     &                                  18,3,3,1,1.D-8,9992)
                        CALL CMATSTRUCT('L*MEZZMP(IQ,1)',LMEZZMP(:,:),
     &                                  18,18,3,3,1,1.D-8,9992)
                        CALL RMATSTRUCT('MEDIFF(IQ,2)',MEDIFF(:,:,2),18,
     &                                  18,3,3,1,1.D-8,9992)
                        CALL CMATSTRUCT('LMEDIFF(IQ)',LMEDIFF(:,:),18,
     &                                  18,3,3,1,1.D-8,9993)
                        CALL CMATSTRUCT('GLESQ(IQ)',GLESQ(:,:,IQ),18,18,
     &                                  3,3,1,1.D-8,999)
                     END IF
                  END IF
C
               END DO LOOP_GLES_DOS
C
            END IF
C======================================================================
C
C
         END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO ! I_TEMP_LAT
C
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== END =
C=======================================================================
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      CALL CPU_TIME(TIME)
      WRITE (6,99007) TIME - TIME0
      CALL STOP_REGULAR(ROUTINE,'NEGF-job done')
C
C   ====================================================================
C
99001 FORMAT ('CDOS for site IQ ',I3,': ',2F20.10)
99002 FORMAT (/,10X,'barrier parameters  ',/,10X,
     &        'number of barrier layers  NQBAR =',I3,/,10X,
     &        'left  barrier range   ',I3,'  ... ',I3,/,10X,
     &        'right barrier range   ',I3,'  ... ',I3,/,10X)
99003 FORMAT (/,79('*'),/,10X,A,/,79('*'),/)
99004 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*              *     *  ******   *****   ******              *'
     &  ,/,10X,
     &  '*              **    *  *       *     *  *                   *'
     &  ,/,10X,
     &  '*              * *   *  *       *        *                   *'
     &  ,/,10X,
     &  '*              *  *  *  ****    *   ***  ****                *'
     &  ,/,10X,
     &  '*              *   * *  *       *     *  *                   *'
     &  ,/,10X,
     &  '*              *    **  *       *     *  *                   *'
     &  ,/,10X,
     &  '*              *     *  ******   *****   *                   *'
     &  ,/,10X,'*',60X,'*',/,10X,'*  VERSION   MO-T3',43X,'*',/,10X,'*',
     &  60X,'*',/,10X,62('*'),//)
99005 FORMAT (/,1X,79('*'),/,10X,'CPA not converged for',I3,
     &        ' energies:',/,10X,'tolerance for CPA-cycle:',F15.7,/,10X,
     &        'maximum deviation:      ',F15.7,/,1X,79('*'),/)
99006 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99007 FORMAT (/,5X,'execution time for <NEGFDRIVE_MO>:',F14.3,' secs',/)
99008 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,24X,'EFERMI =',F10.6,/,10X,'IE  =',I4,5X,
     &        'ZPLS   =',2F10.6,:,/,24X,'ZMIN   =',2F10.6)
99009 FORMAT ('<NEGFDRIVE> sequential mode and finite temp not yet ',
     &        'available')
99010 FORMAT ('<NEGFDRIVE> INFO: setting WF/GF convention to RH')
99011 FORMAT ('<NEGFDRIVE> INFO: setting WF/GF convention to ZJ')
99012 FORMAT ('# BUILDBOT: ',A,':  CDOS for ZPLS =',2F10.6,', IQ =',I2,
     &        /,(1PE22.14))
      END
