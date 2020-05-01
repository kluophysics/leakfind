C*==chrdns.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS(EFCORRECT,CTOTDOS,SHFTEF,TOTDOSX,TOTDOS,TOTNOS,
     &                  ERYD,WEINP,EFERMILD,IEPANEL,IEPATH,IECURR,
     &                  OBS_LTX,OBS_TX,BCOR,BCORS,MEZZ,MEZJ,TSST,MSST,
     &                  TAUT,MSSQ,TAUQ,RHOCHRX,RHOSPNX,RHOORBX,RHO2NSX,
     &                  CMNTTX)
C   ********************************************************************
C   *                                                                  *
C   * subroutine to calculate the dos and magnetic spin and orbital    *
C   * moment densities within an atomic cell                           *
C   *                                                                  *
C   *  n(r) = -1/PI * IM    INT    TRACE G(R,R,E) DE                   *
C   *                     E=0..EF        =                             *
C   *                                                                  *
C   *  m(r) = -1/PI * Im    Int    TRACE BET*SIG*G(R,R,E) dE           *
C   *                     E=0..EF         =   =  =                     *
C   *                                                                  *
C   * radial integrals  DZZ, DZJ, SZZ,  .... calculated in <FPSSITE>   *
C   *                                                                  *
C   * the backscattering contribution and the full crystal             *
C   * contributions are calculated for all quantities                  *
C   * crystal:          G = SUM Z*TAU*Z - SUM Z*J                      *
C   * backscattering:   G = SUM Z*(TAU-T)*Z                            *
C   *                                                                  *
C   *  NOTE: for IREL <= 2 (non- and scalar relativistic) the          *
C   *        minor component has no meaning and is not written/read    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  NOTE: all matrix element related quantities like DOS are        *
C   *        calculated in the LOCAL frame using the corresponding     *
C   *        matrices TMAT_LOC, MEZZ and MEZJ in the LOCAL frame       *
C   *                                                                  *
C   *        all densities are evaluated in the LOCAL frame using the  *
C   *        corresponding projected multiple scattering matrix        *
C   *        TAU_LOC rotated to the LOCAL frame                        *
C   *                                                                  *
C   ********************************************************************
      USE MOD_SYMMETRY,ONLY:MREP_Q,IQREPQ
      USE MOD_SITES,ONLY:IQAT,DROTQ,MROTQ,IQBOT,IQTOP,ITOQ,NOQ,NLQMAD,
     &    NQMAX
      USE MOD_THERMAL,ONLY:X_FT,FTET_FT,NFLUCT,TEMP_LAT,MNT_TEMP_MCS
      USE MOD_ENERGY,ONLY:EFERMI,NEMAX,NETAB,SPLITSS,NEPANEL
      USE MOD_ANGMOM,ONLY:NL,NLMAX,NLMMAX,NMEMAX,NOBSMAX,NKMMAX,NKM,
     &    NCPLWF,TXT_OBS,IDOS,ISMT,IOMT,IHFF,IBND,ITRQ,IWFD,NSPIN,
     &    NMUEMAX,WKM1
      USE MOD_CALCMODE,ONLY:PROGNAME,UPDATE_EFERMI,IREL,LLOYD,
     &    LHS_SOL_EQ_RHS_SOL,L_PROJECTED_ME,IHYPER,KMROT,
     &    THERMAL_VIBRA_FLUCT,NONMAG,MOMENTS_ROTATED
      USE MOD_RMESH,ONLY:NRMAX,NRSFTOT,NSF,LMISF,FLMSF,ISFLM,JRWS,JRMT,
     &    JRCRI,R,FULLPOT,R2DRDI_W_RADINT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NLMFPMAX,NCPLWFMAX,NTMAX,RHOORB,
     &    RHOSPN,RHOCHR,RHOSPNC,RHOCHRC,RHO2NS,NAT,CONC,TXT_T,LTXT_T,
     &    NVALTOT,NT,IMT,MUEORB,MUESPN,NVALT,Z,IKMCPLWF,JFRA,JGRA,ZGRA,
     &    ZFRA,JGLA,JFLA,ZFLA,ZGLA,DOBS_LTX,DOBS_TX,DOBS_BS_EF_TX,
     &    DOBS_BS_EF_LTX,TOTDOS_BS_EF,NLT,QEL
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_CONSTANTS,ONLY:C0,C1,PI,SQRT_4PI,RY_EV
      USE MOD_SCF,ONLY:SCF_CHECK_SPLITSS
      USE MOD_FILES,ONLY:SYSTEM,LSYSTEM,IOTMP,IOTMP1,DATSET,LDATSET,
     &    IPRINT,IFILCBWF,IFILMEZZL,IDUMMY,IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--CHRDNS63
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHRDNS')
      INTEGER NOBS_THERMAL,NOBS
      PARAMETER (NOBS_THERMAL=3,NOBS=3)
      LOGICAL DAMP_SHFTEF,CHECK_OBS_ME,CHECK_CORE_DENSITY
      PARAMETER (DAMP_SHFTEF=.FALSE.,CHECK_OBS_ME=.FALSE.,
     &           CHECK_CORE_DENSITY=.TRUE.)
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
      INTEGER NME,NME_RED
      PARAMETER (NME=7,NME_RED=4)
C
C Dummy arguments
C
      LOGICAL EFCORRECT
      REAL*8 EFERMILD,SHFTEF,TOTDOS,TOTDOSX,TOTNOS
      COMPLEX*16 ERYD,WEINP
      INTEGER IECURR,IEPANEL,IEPATH
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),CMNTTX(NLMFPMAX,NTMAX),
     &       OBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX),OBS_TX(0:3,NOBSMAX,NTMAX),
     &       RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3),RHOCHRX(NRMAX,NTMAX),
     &       RHOORBX(NRMAX,NTMAX),RHOSPNX(NRMAX,NTMAX)
      COMPLEX*16 CTOTDOS(NEMAX),MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BDUM(3),DEL_DOS,DEL_OMT,DEL_SMT,DQ,DRHOCHR(:),DRHOORB(:),
     &       DRHOSPN(:),EVIL(:,:,:),E_BAND,INTEG_CHR,INTEG_ORB,
     &       INTEG_SPN,MIX_DL(:),MIX_IL(:),MROT(3,3),OBS_ABS(:),
     &       OBS_LOC_TX(:,:),OBS_SUM_L_TX(:,:,:),OBS_TOT(:,:),
     &       OBS_TPX(:,:),RSQ,SHFTEFDMP,SHFTEFPREV,WANG_FT,WR,WSPIN,WSUM
      COMPLEX*16 CTOTALDOS,DOBS_LMX(:,:,:,:),DOBS_SUM_L_TX(:,:,:),
     &           EVDL(:,:,:),METAU(:,:),METAU_L(:,:),MEZJL(:,:,:,:),
     &           MEZZL(:,:,:,:),MZBJA(:,:,:,:),MZBZA(:,:,:,:),TMAT(:,:),
     &           TMAT_LOC(:,:),WE
      REAL*8 DDOT,DNRM2
      CHARACTER*80 FILTORQUE(NT),MNTFIL
      INTEGER I,IA_ERR,IFIL_LHSB,IFIL_RHSA,IFLAG,IFLUCT,IFT,IL,IM,IME,
     &        IO,IOBS,IPOL,IQ,IQREP,IR,IRSF,IRTOP,ISF,ISPIN,IT,ITP,
     &        ITREP,J,L,LM,M,MX,MX0,MXTOP,N,NDENS
      LOGICAL INITIALIZE
      CHARACTER*20 STR20
      SAVE DOBS_SUM_L_TX,EVIL,MNTFIL,OBS_SUM_L_TX
C
C*** End of declarations rewritten by SPAG
C
      DATA SHFTEFPREV/1D0/,INITIALIZE/.TRUE./
C
      ALLOCATABLE MZBJA,MZBZA
      ALLOCATABLE DOBS_SUM_L_TX,OBS_SUM_L_TX,DOBS_LMX
      ALLOCATABLE OBS_LOC_TX,OBS_TPX,OBS_TOT,OBS_ABS
      ALLOCATABLE METAU
      ALLOCATABLE MEZJL,MEZZL,METAU_L,MIX_DL,MIX_IL
      ALLOCATABLE TMAT,TMAT_LOC
      ALLOCATABLE EVDL,EVIL
      ALLOCATABLE DRHOCHR,DRHOORB,DRHOSPN
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (DRHOCHR(NRMAX),DRHOORB(NRMAX),DRHOSPN(NRMAX))
      ALLOCATE (OBS_LOC_TX(0:3,NOBSMAX),OBS_ABS(NOBSMAX))
      ALLOCATE (OBS_TPX(0:3,NOBSMAX),OBS_TOT(0:3,NOBSMAX))
      ALLOCATE (TMAT(NKMMAX,NKMMAX),TMAT_LOC(NKMMAX,NKMMAX))
C
      ALLOCATE (DOBS_LMX(0:3,NOBSMAX,NLMAX,NMUEMAX))
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (OBS_SUM_L_TX(0:3,NOBSMAX,NTMAX))
         ALLOCATE (DOBS_SUM_L_TX(0:3,NOBSMAX,NTMAX))
C
C----------------------------------------------------------- TEMPERATURE
         IF ( THERMAL_VIBRA_FLUCT ) THEN
C
            MNTFIL = DATSET(1:LDATSET)//'_mnt_T.dat'
C
            OPEN (IOTMP,FILE=MNTFIL)
            WRITE (IOTMP,'(A)') '# data file: '//
     &                          MNTFIL(1:LEN_TRIM(MNTFIL))
            CLOSE (IOTMP)
C
C--------------------------------------- files for temporary test output
            DO IT = 1,NT
               FILTORQUE(IT) = DATSET(1:LDATSET)//'_Torque_T_'//
     &                         TXT_T(IT)(1:LTXT_T(IT))//'.dat'
            END DO
C
         END IF
C----------------------------------------------------------- TEMPERATURE
C
C-------------------------------------------------- check index settings
         IF ( IREL.LE.1 .AND. NKM.NE.NL**2 )
     &         CALL STOP_MESSAGE(ROUTINE,'IREL <= 1 AND NKM <> NL**2')
         IF ( IREL.LE.2 .AND. NKMMAX.NE.2*NLMMAX )
     &        CALL STOP_MESSAGE(ROUTINE,
     &        'IREL <= 2 AND NKMMAX.NE.2*NLMMAX')
         IF ( NME.LT.NME_RED )
     &        CALL STOP_MESSAGE(ROUTINE,'NME < NME_RED')
         IF ( NME.GT.NMEMAX ) CALL STOP_MESSAGE(ROUTINE,'NME > NMEMAX')
         IF ( IDOS.GT.NME_RED )
     &         CALL STOP_MESSAGE(ROUTINE,'IDOS > NME_RED')
         IF ( ISMT.GT.NME_RED )
     &         CALL STOP_MESSAGE(ROUTINE,'ISMT > NME_RED')
         IF ( IOMT.GT.NME_RED )
     &         CALL STOP_MESSAGE(ROUTINE,'IOMT > NME_RED')
         IF ( IHFF.GT.NME_RED )
     &         CALL STOP_MESSAGE(ROUTINE,'IHFF > NME_RED')
         IF ( IBND.GT.NME ) CALL STOP_MESSAGE(ROUTINE,'IBND > NME')
C
         INITIALIZE = .FALSE.
C
      END IF
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      IFIL_RHSA = IFILCBWF
C
      CALL SET_IFIL_LHS(IFIL_RHSA,IFIL_LHSB)
C
C=======================================================================
C
      DOBS_LMX(:,:,:,:) = 0D0
      OBS_SUM_L_TX(:,:,:) = 0D0
      DOBS_SUM_L_TX(:,:,:) = 0D0
C
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGRA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGRA')
      IF ( IREL.EQ.3 ) THEN
         ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JGLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFLA(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGLA')
      END IF
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      IF ( FULLPOT ) THEN
         CMNTTX(1,1) = 0D0
         RHOCHRX(1,1) = 0D0
         RHOORBX(1,1) = 0D0
         RHOSPNX(1,1) = 0D0
         DRHOCHR(1) = 0D0
         DRHOORB(1) = 0D0
         DRHOSPN(1) = 0D0
C
         ALLOCATE (METAU(NKMMAX,NMEMAX))
C
C-----------------------------------------------------------------------
C         read l-projected matrix elements MEZZL from IFILMEZZL
C-----------------------------------------------------------------------
         IF ( L_PROJECTED_ME ) THEN
C
            ALLOCATE (MIX_DL(NMEMAX))
            ALLOCATE (MIX_IL(NMEMAX))
            ALLOCATE (EVDL(NLMAX,NTMAX,NMEMAX))
            EVDL(:,:,:) = C0
            IF ( .NOT.ALLOCATED(EVIL) )
     &           ALLOCATE (EVIL(NLMAX,NTMAX,NMEMAX))
            IF ( IECURR.EQ.1 .AND. IEPATH.EQ.1 ) EVIL(:,:,:) = 0D0
C
            ALLOCATE (METAU_L(NKMMAX,NMEMAX))
            ALLOCATE (MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
            ALLOCATE (MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEZZL')
C
            REWIND IFILMEZZL
C
C--------------------------------------------------------- position file
            DO ITP = 1,ITBOT - 1
C
               READ (IFILMEZZL) IDUMMY,MEZZL,MEZJL
C
               IF ( IDUMMY.EQ.ITBOT-1 ) EXIT
               IF ( IDUMMY.EQ.ITBOT .AND. ITP.EQ.1 ) THEN
                  REWIND IFILMEZZL
                  EXIT
               END IF
               IF ( IDUMMY.GE.ITBOT ) CALL STOP_MESSAGE(ROUTINE,
     &              'reading IFILMEZZL: IDUMMY >= ITBOT')
            END DO
C
         END IF
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      ELSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      NLQMAD = 1
C ----------------------------------------------------------------------
C
         RHO2NSX(1,1,1,1) = 0D0
         L_PROJECTED_ME = .FALSE.
         LHS_SOL_EQ_RHS_SOL = .TRUE.
C
      END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C         calculate full vector in case of rotated magnetic moment
C          or thermal spin fluctuations and/or lattice vibrations
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
      IF ( THERMAL_VIBRA_FLUCT .OR. (MOMENTS_ROTATED .AND. IREL.EQ.3) )
     &     THEN
         ALLOCATE (MZBJA(NKMMAX,NKMMAX,3,NOBS))
         ALLOCATE (MZBZA(NKMMAX,NKMMAX,3,NOBS),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBZA')
      END IF
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C
      IF ( IREL.GT.1 ) THEN
         WSPIN = 1D0
      ELSE
         WSPIN = 2D0
      END IF
      WE = WEINP
C
C ======================================================================
C
      IF ( IREL.LT.3 ) THEN
         NDENS = NSPIN
      ELSE IF ( PROGNAME(4:6).EQ.'SCF' ) THEN
         NDENS = 2
      ELSE
         NDENS = 3
      END IF
C
C ======================================================================
C ======================================================================
C              correct Fermi energy to ensure charge neutrality
C              and correct charge density and related quantities
C ======================================================================
C ======================================================================
C
      IF ( EFCORRECT ) THEN
C
         DQ = TOTNOS - NVALTOT/WSPIN
C
         IF ( (SUB_SYSTEM(1:6).EQ.'I-ZONE' .AND. SYSTEM_TYPE(1:3)
     &        .NE.'VIV') .OR. SUB_SYSTEM(1:6).EQ.'R-BULK' .OR. 
     &        .NOT.UPDATE_EFERMI ) THEN
C
            SHFTEF = 0.0D0
            WE = 0D0
C
         ELSE
C
            SHFTEF = -DQ/TOTDOS
            WE = -DQ/TOTDOSX
C
         END IF
C
         IF ( SCF_CHECK_SPLITSS ) WRITE (6,99001) TOTNOS,DQ,TOTDOS,
     &        TOTDOSX
C
C------------------------------------------------- damp shift for EFERMI
         IF ( DAMP_SHFTEF ) THEN
            IF ( SHFTEF*SHFTEFPREV.GE.0D0 ) THEN
               SHFTEFDMP = 0.6D0
            ELSE
               SHFTEFDMP = 0.1D0
            END IF
            SHFTEFPREV = SHFTEF
            SHFTEF = SHFTEF*SHFTEFDMP
         END IF
C-----------------------------------------------------------------------
C
         EFERMI = EFERMI + SHFTEF
C
         WRITE (6,'(/)')
         IF ( UPDATE_EFERMI ) THEN
            WRITE (6,99002) SYSTEM(1:LSYSTEM),
     &                ' results extrapolated to corrected Fermi energy:'
         ELSE
            WRITE (6,99002) SYSTEM(1:LSYSTEM),
     &                   ' results integrated up to fixed Fermi energy:'
         END IF
         IF ( LLOYD ) THEN
            IF ( UPDATE_EFERMI ) THEN
               WRITE (6,99003) DQ,SHFTEF,EFERMI,EFERMILD
            ELSE
               WRITE (6,99004) DQ,EFERMI
            END IF
            IF ( FULLPOT ) EFERMI = EFERMILD
         ELSE IF ( UPDATE_EFERMI ) THEN
            WRITE (6,99003) DQ,SHFTEF,EFERMI
         ELSE
            WRITE (6,99004) DQ,EFERMI
         END IF
C
      END IF
C
C ======================================================================
C
      IF ( .NOT.SPLITSS ) THEN
         STR20 = 'CRYSTAL TERMS       '
      ELSE IF ( IEPATH.EQ.1 ) THEN
         STR20 = 'BACKSCATTERING TERMS'
      ELSE
         STR20 = 'SINGLE SITE TERMS   '
      END IF
C
      CTOTALDOS = 0.0D0
      TOTDOS = 0.0D0
      TOTNOS = 0.0D0
      MUESPN = 0.0D0
      MUEORB = 0.0D0
      E_BAND = 0.0D0
C
      DOBS_LTX(:,:,:,:) = C0
      DOBS_TX(:,:,:) = C0
      DOBS_SUM_L_TX(:,:,:) = C0
C
      OBS_SUM_L_TX(:,:,:) = 0.0D0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT A
      LOOP_IT_A:DO IT = ITBOT,ITTOP
C
         N = NKM
         M = NKMMAX
C
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
C=======================================================================
C                   read valence band wave functions
C                   NCPLWF  =   NCPLWF_RA =   NCPLWF_LB
C                 IKMCPLWF  = IKMCPLWF_RA = IKMCPLWF_LB
C=======================================================================
C
         IF ( IREL.LE.2 .AND. .NOT.FULLPOT )
     &        CALL WAVFUN_READ_SRA(IFILCBWF,IT,1,ZGRA,JGRA,IRTOP,NCPLWF,
     &        IKMCPLWF)
C
         IF ( IREL.EQ.3 ) THEN
C
C ----------------------------------------- read in wavefunctions for RA
C
            CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,
     &                           IRTOP,NCPLWF,IKMCPLWF)
C
            IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
               ZGLA(:,:,:) = ZGRA(:,:,:)
               ZFLA(:,:,:) = ZFRA(:,:,:)
               JGLA(:,:,:) = JGRA(:,:,:)
               JFLA(:,:,:) = JFRA(:,:,:)
C
            ELSE
C ----------------------------------------- read in wavefunctions for LA
C
               CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGLA,ZFLA,JGLA,JFLA,
     &                              IRTOP,NCPLWF,IKMCPLWF)
C
            END IF
         END IF
C=======================================================================
C
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C      calculate matrix elements  MZBZA and MZBJA in the LOCAL frame
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
         IF ( THERMAL_VIBRA_FLUCT .OR. (MOMENTS_ROTATED .AND. IREL.EQ.3)
     &        ) THEN
C-------------------------------------------------------------- IREL = 3
C                                                              LHS = RHS
            IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
               CALL CALC_OBS_ME(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA,MZBJA,NOBS,'S>C',
     &                          CHECK_OBS_ME)
C
C-------------------------------------------------------------- IREL = 3
C                                                             LHS != RHS
            ELSE
C
               CALL CALC_OBS_ME(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA,MZBJA,NOBS,'S>C',
     &                          CHECK_OBS_ME)
C
            END IF
C
         END IF
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C
C----------------------------------------------------------- TEMPERATURE
         IF ( .NOT.THERMAL_VIBRA_FLUCT ) THEN
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
            CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT)
C
C--- provide site diagonal multiple scattering matrix in the LOCAL frame
C
            IF ( KMROT.EQ.0 ) THEN
C
               TMAT_LOC(:,:) = TMAT(:,:)
C
            ELSE
C
               IQ = IQAT(1,IT)
C
               CALL ROTATE(TMAT,'G->L',TMAT_LOC,N,DROTQ(1,1,IQ),M)
C
            END IF
C
         ELSE
C---------------------------------------------------- finite temperature
C
            CALL THERMAL_PROPERTIES(IT,IECURR,IEPATH,WE,TSST,MSST,MSSQ,
     &                              TAUQ,TMAT,MZBZA,MZBJA,DOBS_TX,
     &                              OBS_TX,NOBS)
C
            CALL TAUGFCONV(MSST(1,1,IT),TMAT,TMAT_LOC)
C
         END IF
C----------------------------------------------------------- TEMPERATURE
C
         IM = IMT(IT)
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         IF ( FULLPOT ) THEN
C
            DO IME = 1,NME
C
               CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ(1,1,IT,IME),M,
     &                    TMAT_LOC,M,C0,WKM1,M)
C
               DO I = 1,N
                  METAU(I,IME) = WKM1(I,I) - CPRE*MEZJ(I,I,IT,IME)
               END DO
            END DO
C
C-----------------------------------------------------------------------
C         read l-projected matrix elements MEZZL from IFILMEZZL
C-----------------------------------------------------------------------
            IF ( L_PROJECTED_ME ) THEN
               READ (IFILMEZZL) IDUMMY,MEZZL,MEZJL
C
               IF ( IDUMMY.NE.IT )
     &               CALL STOP_MESSAGE(ROUTINE,'reading IFILMEZZL')
            END IF
C-----------------------------------------------------------------------
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         ELSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
            DRHOCHR(1:NRMAX) = 0.0D0
            DRHOORB(1:NRMAX) = 0.0D0
            DRHOSPN(1:NRMAX) = 0.0D0
C
         END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         LOOP_IL:DO IL = 1,NLT(IT)
C
            L = IL - 1
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
            IF ( FULLPOT ) THEN
C
               DO ISPIN = 1,NSPIN
                  IF ( IREL.LE.2 ) THEN
                     MX0 = (L**2) + (NL**2)*(ISPIN-1)
                     MXTOP = 2*L + 1
                  ELSE
                     MX0 = 2*(L-1)*L + L + L
                     MXTOP = 4*L + 2
                  END IF
                  DO IME = 1,NME
                     DO MX = 1,MXTOP
                        DOBS_LTX(0,IME,IL,IT) = DOBS_LTX(0,IME,IL,IT)
     &                     + METAU(MX0+MX,IME)
                     END DO
                  END DO
               END DO
C
C-----------------------------------------------------------------------
C                  calculate l-projected properties
C-----------------------------------------------------------------------
               IF ( L_PROJECTED_ME ) THEN
C
                  DO IME = 1,NME
                     CALL ZGEMM('N','N',N,N,N,CPRE,MEZZL(1,1,IL,IME),M,
     &                          TMAT_LOC,M,C0,WKM1,M)
C
                     DO I = 1,N
                        METAU_L(I,IME) = WKM1(I,I)
     &                     - CPRE*MEZJL(I,I,IL,IME)
                     END DO
                  END DO
C
                  DO ISPIN = 1,NSPIN
                     IF ( IREL.LE.2 ) THEN
                        MX0 = (L**2) + (NL**2)*(ISPIN-1)
                        MXTOP = 2*L + 1
                     ELSE
                        MX0 = 2*(L-1)*L + L + L
                        MXTOP = 4*L + 2
                     END IF
                     DO IME = 1,NME
                        DO MX = 1,MXTOP
                           EVDL(IL,IT,IME) = EVDL(IL,IT,IME)
     &                        + METAU_L(MX0+MX,IME)
                        END DO
                     END DO
                  END DO
C
                  EVIL(IL,IT,1:NME) = EVIL(IL,IT,1:NME)
     &                                + DIMAG(WE*EVDL(IL,IT,1:NME))
C
               END IF
C-----------------------------------------------------------------------
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
C ----------------------------------------------------------------------
C                                           supplies DOBS_LTX via module
C
            ELSE IF ( IREL.LE.2 ) THEN
C
               CALL CHRDNS_SRA(IT,IL,WE,IRTOP,DOBS_LMX,NME_RED,DRHOCHR,
     &                         DRHOSPN,MEZZ,MEZJ,TMAT_LOC,ZGRA,JGRA)
C
            ELSE
C
               CALL CHRDNS_REL(IT,IL,WE,IRTOP,DOBS_LMX,NME_RED,DRHOCHR,
     &                         DRHOORB,DRHOSPN,MEZZ,MEZJ,TMAT_LOC,ZGRA,
     &                         ZFRA,JGRA,JFRA)
C
C
            END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
            DOBS_LTX(0,IBND,IL,IT) = WSPIN*DOBS_LTX(0,IDOS,IL,IT)*ERYD
C
            DOBS_TX(0,:,IT) = DOBS_TX(0,:,IT) + DOBS_LTX(0,:,IL,IT)
C
            OBS_LTX(0,:,IL,IT) = OBS_LTX(0,:,IL,IT)
     &                           + DIMAG(WE*DOBS_LTX(0,:,IL,IT))
C
            OBS_TX(0,:,IT) = OBS_TX(0,:,IT)
     &                       + DIMAG(WE*DOBS_LTX(0,:,IL,IT))
C ......................................................................
C
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
            IF ( (IPRINT.GT.0) .AND. EFCORRECT ) WRITE (6,99005) IECURR,
     &           ERYD,L,IT,TXT_T(IT),STR20,OBS_LTX(0,IDOS,IL,IT),
     &           OBS_LTX(0,ISMT,IL,IT),OBS_LTX(0,IOMT,IL,IT),
     &           (OBS_LTX(0,IHFF,IL,IT)*1D-3),DOBS_LTX(0,IDOS,IL,IT),
     &           DOBS_LTX(0,ISMT,IL,IT),DOBS_LTX(0,IOMT,IL,IT),
     &           (DOBS_LTX(0,IHFF,IL,IT)*1D-6)
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
         END DO LOOP_IL
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         CTOTALDOS = CTOTALDOS + DOBS_TX(0,IDOS,IT)*CONC(IT)*NAT(IT)
C
         TOTDOS = TOTDOS + DIMAG(DOBS_TX(0,IDOS,IT))*CONC(IT)*NAT(IT)
         TOTNOS = TOTNOS + OBS_TX(0,IDOS,IT)*CONC(IT)*NAT(IT)
         MUESPN = MUESPN + OBS_TX(0,ISMT,IT)*CONC(IT)*NAT(IT)
         MUEORB = MUEORB + OBS_TX(0,IOMT,IT)*CONC(IT)*NAT(IT)
         E_BAND = E_BAND + OBS_TX(0,IBND,IT)*CONC(IT)*NAT(IT)
C
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
         IF ( IPRINT.GT.0 ) THEN
C
            WRITE (6,99006) IECURR,ERYD,IT,TXT_T(IT),OBS_TX(0,IDOS,IT),
     &                      OBS_TX(0,ISMT,IT),OBS_TX(0,IOMT,IT),
     &                      (OBS_TX(0,IHFF,IT)*1D-3),
     &                      DIMAG(DOBS_TX(0,IDOS,IT)),
     &                      DIMAG(DOBS_TX(0,ISMT,IT)),
     &                      DIMAG(DOBS_TX(0,IOMT,IT)),
     &                      (DIMAG(DOBS_TX(0,IHFF,IT))*1D-3)
C
            IF ( IT.LT.NT ) THEN
               WRITE (6,'(1X,79(''-''))')
            ELSE IF ( EFCORRECT ) THEN
               IF ( MOMENTS_ROTATED ) THEN
                  WRITE (6,99010) TOTDOS,TOTNOS
               ELSE
                  WRITE (6,99010) TOTDOS,TOTNOS,MUESPN,MUEORB
               END IF
            ELSE
               WRITE (6,'('' '',79(''=''))')
            END IF
         END IF
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
         IF ( .NOT.(FULLPOT) ) THEN
C
C ----------------------------------------------------------------------
C                   calculate electric dipole moments
C ----------------------------------------------------------------------
C
            IF ( NLQMAD.NE.1 ) THEN
C
               IF ( IREL.NE.3 ) THEN
C
                  CALL CHRDNS_DIP_MNT_SRA(IT,IM,WE,IRTOP,TMAT_LOC,ZGRA,
     &               JGRA,WSPIN,CMNTTX)
C
               ELSE
C
                  CALL CHRDNS_DIP_MNT_REL(IT,IM,WE,IRTOP,TMAT_LOC,ZGRA,
     &               ZFRA,JGRA,JFRA,CMNTTX)
C
               END IF
C
            END IF
C
C ======================================================================
C                                                          add densities
            DO IR = 1,IRTOP
               RHOCHRX(IR,IT) = RHOCHRX(IR,IT) + DRHOCHR(IR)
               RHOSPNX(IR,IT) = RHOSPNX(IR,IT) + DRHOSPN(IR)
               RHOORBX(IR,IT) = RHOORBX(IR,IT) + DRHOORB(IR)
            END DO
C------------------------------------------------------------- IREL <= 2
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
         ELSE IF ( IREL.LE.2 ) THEN
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
            CALL CHRDNS_SRA_FP(IT,IRTOP,WE,IFILCBWF,TMAT_LOC,RHO2NSX,
     &                         ZGRA,JGRA)
C
C-------------------------------------------------------------- IREL = 3
C                                                              LHS = RHS
         ELSE IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
            CALL CHRDNS_REL_FP(NDENS,IT,IRTOP,WE,IFILCBWF,TMAT_LOC,
     &                         RHO2NSX,ZGRA,ZFRA,JGRA,JFRA,ZGRA,ZFRA,
     &                         JGRA,JFRA)
C
C-------------------------------------------------------------- IREL = 3
C                                                             LHS != RHS
         ELSE
C
            CALL CHRDNS_REL_FP(NDENS,IT,IRTOP,WE,IFILCBWF,TMAT_LOC,
     &                         RHO2NSX,ZGLA,ZFLA,JGLA,JFLA,ZGRA,ZFRA,
     &                         JGRA,JFRA)
C
         END IF
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C       calculate full vector in case of rotated magnetic moment
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
         IF ( MOMENTS_ROTATED .AND. .NOT.THERMAL_VIBRA_FLUCT .AND. 
     &        IREL.EQ.3 ) CALL CALC_OBS(IT,TMAT_LOC,MZBZA,MZBJA,WE,
     &                                  OBS_TX,NOBS)
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      END DO LOOP_IT_A
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT A
C
C
C ======================================================================
C
      IF ( .NOT.EFCORRECT ) THEN
C
C------------------------------------- store backscattering contribution
C
         IF ( SPLITSS .AND. (IEPANEL.EQ.NEPANEL) .AND. (IEPATH.EQ.1)
     &        .AND. (IECURR.EQ.NETAB(IEPATH)) ) THEN
C
            TOTDOS_BS_EF = TOTDOS
C
            DOBS_BS_EF_TX(:,:,ITBOT:ITTOP) = DOBS_TX(:,:,ITBOT:ITTOP)
            DOBS_BS_EF_LTX(:,:,:,ITBOT:ITTOP)
     &         = DOBS_LTX(:,:,:,ITBOT:ITTOP)
C
         END IF
C
         CTOTDOS(IECURR) = CTOTALDOS
C
C???????????????????????????????????????????????????????????????????????
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         IF ( FULLPOT ) THEN
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         ELSE IF ( SCF_CHECK_SPLITSS .AND. (IECURR.EQ.NETAB(IEPATH)) )
     &             THEN
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
            WRITE (6,99022)
            DO IT = ITBOT,ITTOP
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
               INTEG_CHR = DDOT(IRTOP,R2DRDI_W_RADINT(1,IM),1,
     &                     RHOCHRX(1,IT),1)
C
               WRITE (6,99021) IT,'charge ',OBS_TX(0,IDOS,IT),INTEG_CHR,
     &                         OBS_TX(0,IDOS,IT)/INTEG_CHR
C
               IF ( .NOT.NONMAG ) THEN
                  INTEG_SPN = DDOT(IRTOP,R2DRDI_W_RADINT(1,IM),1,
     &                        RHOSPNX(1,IT),1)
                  INTEG_ORB = DDOT(IRTOP,R2DRDI_W_RADINT(1,IM),1,
     &                        RHOORBX(1,IT),1)
C
                  IF ( IREL.GE.2 .AND. ABS(INTEG_SPN).GT.TOL .AND. 
     &                 ABS(OBS_TX(0,ISMT,IT)).GT.TOL ) WRITE (6,99021)
     &                 IT,'spin   ',OBS_TX(0,ISMT,IT),INTEG_SPN,
     &                 OBS_TX(0,ISMT,IT)/INTEG_SPN
                  IF ( IREL.GE.3 .AND. ABS(INTEG_ORB).GT.TOL .AND. 
     &                 ABS(OBS_TX(0,IOMT,IT)).GT.TOL ) WRITE (6,99021)
     &                 IT,'orbital',OBS_TX(0,IOMT,IT),INTEG_ORB,
     &                 OBS_TX(0,IOMT,IT)/INTEG_ORB
               END IF
            END DO
C
         END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
         CALL CHRDNS_CHECK_OBS(NME,DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX)
C???????????????????????????????????????????????????????????????????????
C
         DEALLOCATE (ZGRA,ZFRA,JGRA,JFRA)
         IF ( ALLOCATED(JFLA) ) DEALLOCATE (JGLA,JFLA,ZFLA,ZGLA)
C
         RETURN
C
      END IF
C
C
C ======================================================================
C ======================================================================
C                   correction at the Fermi energy
C ======================================================================
C ======================================================================
C
C===========================================================TEMPERATURE
C             project m_z in local frame onto the z-axis of global frame
C                                   in case of thermal spin fluctuations
C
      IF ( THERMAL_VIBRA_FLUCT .AND. NFLUCT.GT.1 ) THEN
C
         OPEN (IOTMP,FILE=MNTFIL,POSITION='APPEND')
C
         DO IT = ITBOT,ITTOP
C
            IFT = (IT-1)*NFLUCT
            WSUM = 0D0
            DO IFLUCT = 1,NFLUCT
               IFT = IFT + 1
               WANG_FT = COS(FTET_FT(IFT)*PI/180D0)
               WSUM = WSUM + X_FT(IFT)*WANG_FT
            END DO
C
            DOBS_LTX(:,ISMT,:,IT) = DOBS_LTX(:,ISMT,:,IT)*WSUM
            DOBS_LTX(:,IOMT,:,IT) = DOBS_LTX(:,IOMT,:,IT)*WSUM
C
            DOBS_TX(0,ISMT,IT) = DOBS_TX(0,ISMT,IT)*WSUM
            DOBS_TX(0,IOMT,IT) = DOBS_TX(0,IOMT,IT)*WSUM
C
            OBS_TX(0,ISMT,IT) = OBS_TX(0,ISMT,IT)*WSUM
            OBS_TX(0,IOMT,IT) = OBS_TX(0,IOMT,IT)*WSUM
C
            DO IOBS = 2,NOBS_THERMAL
               OBS_ABS(IOBS) = DNRM2(3,OBS_TX(1,IOBS,IT),1)
            END DO
C
C------------------------------------------------- temporary test output
            WRITE (6,'(10X,A,F8.2,5X,A,I3,5X,2(A,F20.10),A)') 'T = ',
     &             TEMP_LAT,'IT = ',IT,'TORQUE: ',OBS_TX(1,ITRQ,IT),
     &             ' Ry; ',OBS_TX(1,ITRQ,IT)*RY_EV*10**3,' meV '
            OPEN (IOTMP1,FILE=FILTORQUE(IT),POSITION='APPEND')
            WRITE (IOTMP1,99029) TEMP_LAT,OBS_TX(1,ITRQ,IT)
            CLOSE (IOTMP1)
C------------------------------------------------- temporary test output
C
            WRITE (IOTMP,99027) TEMP_LAT,MNT_TEMP_MCS,OBS_TX(3,IDOS,IT),
     &                          (OBS_ABS(IOBS),
     &                          (OBS_TX(IPOL,IOBS,IT),IPOL=1,3),IOBS=2,
     &                          NOBS_THERMAL),OBS_TX(0,ISMT,IT)/WSUM,
     &                          OBS_TX(0,ISMT,IT),OBS_TX(0,IOMT,IT)
     &                          /WSUM,OBS_TX(0,IOMT,IT),WSUM
C
         END DO
C
         CLOSE (IOTMP)
C
      END IF
C ======================================================================
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      IF ( FULLPOT ) THEN
         RHOCHR(:,:) = 0D0
         RHOSPN(:,:) = 0D0
         RHOORB(:,:) = 0D0
      END IF
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
      IFLAG = 0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT_B:DO IT = ITBOT,ITTOP
C
         WRITE (6,99006) (IECURR+1),EFERMI,0.0D0,IT,TXT_T(IT)
C
         BDUM(1) = BCORS(IT)*1D-3
         BDUM(2) = (BCOR(IT)-BCORS(IT))*1D-3
         BDUM(3) = BCOR(IT)*1D-3
C
         IF ( IHYPER.EQ.0 ) THEN
            BDUM(:) = 0D0
            OBS_LTX(:,IHFF,:,IT) = 0D0
            OBS_TX(:,IHFF,IT) = 0D0
            BCOR(IT) = 0D0
         END IF
C
C-----------------------------------------------------------------------
C                  l-projected matrix quantities
C-----------------------------------------------------------------------
         IF ( L_PROJECTED_ME ) THEN
C
            WRITE (6,99007) (DIMAG(EVDL(IL,IT,IDOS)),EVIL(IL,IT,IDOS),
     &                      DIMAG(EVDL(IL,IT,ISMT)),EVIL(IL,IT,ISMT),
     &                      DIMAG(EVDL(IL,IT,IOMT)),EVIL(IL,IT,IOMT),
     &                      EVIL(IL,IT,IHFF)*1D-3,BDUM(IL),IL=1,
     &                      MIN(3,NLT(IT)))
            IF ( NLT(IT).GT.3 ) WRITE (6,99008)
     &                                 (DIMAG(EVDL(IL,IT,IDOS)),
     &                                 EVIL(IL,IT,IDOS),
     &                                 DIMAG(EVDL(IL,IT,ISMT)),
     &                                 EVIL(IL,IT,ISMT),
     &                                 DIMAG(EVDL(IL,IT,IOMT)),
     &                                 EVIL(IL,IT,IOMT),EVIL(IL,IT,IHFF)
     &                                 *1D-3,IL=4,NLT(IT))
C
            DOBS_SUM_L_TX(0,IME,IT) = 0D0
            OBS_SUM_L_TX(0,IME,IT) = 0D0
            DO IME = 1,NME
               DO IL = 1,NLT(IT)
                  DOBS_SUM_L_TX(0,IME,IT) = DOBS_SUM_L_TX(0,IME,IT)
     &               + EVDL(IL,IT,IME)
                  OBS_SUM_L_TX(0,IME,IT) = OBS_SUM_L_TX(0,IME,IT)
     &               + EVIL(IL,IT,IME)
               END DO
            END DO
C
            DO IME = 1,NME
               MIX_DL(IME) = DIMAG(DOBS_TX(0,IME,IT)-DOBS_SUM_L_TX(0,IME
     &                       ,IT))
               MIX_IL(IME) = OBS_TX(0,IME,IT) - OBS_SUM_L_TX(0,IME,IT)
            END DO
C
            WRITE (6,99009) 'mix',MIX_DL(IDOS),MIX_IL(IDOS),MIX_DL(ISMT)
     &                      ,MIX_IL(ISMT),MIX_DL(IOMT),MIX_IL(IOMT),
     &                      MIX_IL(IHFF)*1D-3
C
C-----------------------------------------------------------------------
         ELSE
            WRITE (6,99007) (DIMAG(DOBS_LTX(0,IDOS,IL,IT)),OBS_LTX(0,
     &                      IDOS,IL,IT),DIMAG(DOBS_LTX(0,ISMT,IL,IT)),
     &                      OBS_LTX(0,ISMT,IL,IT),
     &                      DIMAG(DOBS_LTX(0,IOMT,IL,IT)),
     &                      OBS_LTX(0,IOMT,IL,IT),OBS_LTX(0,IHFF,IL,IT)
     &                      *1D-3,BDUM(IL),IL=1,MIN(3,NLT(IT)))
            IF ( NLT(IT).GT.3 ) WRITE (6,99008)
     &                                 (DIMAG(DOBS_LTX(0,IDOS,IL,IT)),
     &                                 OBS_LTX(0,IDOS,IL,IT),
     &                                 DIMAG(DOBS_LTX(0,ISMT,IL,IT)),
     &                                 OBS_LTX(0,ISMT,IL,IT),
     &                                 DIMAG(DOBS_LTX(0,IOMT,IL,IT)),
     &                                 OBS_LTX(0,IOMT,IL,IT),
     &                                 OBS_LTX(0,IHFF,IL,IT)*1D-3,IL=4,
     &                                 NLT(IT))
C
         END IF
C-----------------------------------------------------------------------
C
         WRITE (6,99009) 'sum',DIMAG(DOBS_TX(0,IDOS,IT)),
     &                   OBS_TX(0,IDOS,IT),DIMAG(DOBS_TX(0,ISMT,IT)),
     &                   OBS_TX(0,ISMT,IT),DIMAG(DOBS_TX(0,IOMT,IT)),
     &                   OBS_TX(0,IOMT,IT),(OBS_TX(0,IHFF,IT)*1D-3),
     &                   ((OBS_TX(0,IHFF,IT)+BCOR(IT))*1D-3)
C
         IF ( .NOT.FULLPOT ) THEN
            IF ( NT.GT.1 ) WRITE (6,99023) OBS_TX(0,IBND,IT)
C
C            WRITE (*,'(A,i3,3f10.5)') 'dipole moment ',IT,CMNTTX(2:4,IT)
            WRITE (*,'(A,i3,3f24.16)') 'dipole moment ',IT,
     &                                 CMNTTX(2:4,IT)
         END IF
C
         IF ( IT.LT.ITTOP ) THEN
            WRITE (6,'(1X,79(''-''))')
         ELSE IF ( MOMENTS_ROTATED .AND. IREL.EQ.3 ) THEN
            WRITE (6,99011) TOTDOS,TOTNOS,E_BAND
         ELSE
            WRITE (6,99012) TOTDOS,TOTNOS,MUESPN,MUEORB,E_BAND
         END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C NOTE: HFF taken out for the moment as it is too sensitive
         IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99028)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                            DIMAG(DOBS_TX(0,IDOS,IT)),
     &                            OBS_TX(0,IDOS,IT),OBS_TX(0,ISMT,IT),
     &                            OBS_TX(0,IOMT,IT),OBS_TX(0,IBND,IT)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO LOOP_IT_B
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C         calculate full vector in case of rotated magnetic moment
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      IF ( MOMENTS_ROTATED .AND. IREL.EQ.3 ) THEN
C
         WRITE (6,99013)
         OBS_TOT(:,:) = 0D0
C
         DO IQ = IQBOT,IQTOP
C
            IQREP = IQREPQ(IQ)
C
            DO IO = 1,NOQ(IQ)
C
               IT = ITOQ(IO,IQ)
               ITREP = ITOQ(IO,IQREP)
C
C-----------------------------------------------------------------------
               IF ( IQREP.EQ.IQ ) THEN
C
                  MROT(1:3,1:3) = TRANSPOSE(MROTQ(1:3,1:3,IQ))
C
                  OBS_LOC_TX(:,1:NOBS) = OBS_TX(:,1:NOBS,IT)
                  DO IOBS = 2,NOBS
                     OBS_TX(1:3,IOBS,IT)
     &                  = MATMUL(MROT(1:3,1:3),OBS_LOC_TX(1:3,IOBS))
                     OBS_ABS(IOBS) = DNRM2(3,OBS_LOC_TX(1,IOBS),1)
                  END DO
                  OBS_TOT(:,:) = OBS_TOT(:,:) + CONC(IT)*OBS_TX(:,:,IT)
C
                  WRITE (6,99016) IQ,IT,CONC(IT),
     &                            (TXT_OBS(IOBS),OBS_ABS(IOBS),
     &                            OBS_TX(1:3,IOBS,IT),
     &                            OBS_LOC_TX(1:3,IOBS),IOBS=2,3)
C
C....................................................check consistency
C                           <CALC_OBS>             <CHRDNS>  z-component
                  DEL_DOS = OBS_LOC_TX(3,IDOS) - OBS_TX(0,IDOS,IT)
                  DEL_SMT = OBS_LOC_TX(3,ISMT) - OBS_TX(0,ISMT,IT)
                  DEL_OMT = OBS_LOC_TX(3,IOMT) - OBS_TX(0,IOMT,IT)
                  IF ( ABS(DEL_DOS).GT.TOL .OR. ABS(DEL_SMT).GT.TOL .OR. 
     &                 ABS(DEL_OMT).GT.TOL ) WRITE (6,99026) DEL_DOS,
     &                 DEL_SMT,DEL_OMT
C-----------------------------------------------------------------------
               ELSE
C
                  MROT(1:3,1:3) = TRANSPOSE(MREP_Q(1:3,1:3,IQ))
C
                  OBS_TPX(:,1:NOBS) = OBS_TX(:,1:NOBS,IT)
                  DO IOBS = 2,NOBS
                     OBS_TPX(1:3,IOBS)
     &                  = MATMUL(MROT(1:3,1:3),OBS_TX(1:3,IOBS,ITREP))
                     OBS_ABS(IOBS) = DNRM2(3,OBS_TPX(1,IOBS),1)
                  END DO
C
                  OBS_TOT(:,:) = OBS_TOT(:,:) + CONC(IT)*OBS_TPX(:,:)
C
                  WRITE (6,99017) IQ,IT,CONC(IT),
     &                            (TXT_OBS(IOBS),OBS_ABS(IOBS),
     &                            OBS_TPX(1:3,IOBS),IOBS=2,3)
C
               END IF
C-----------------------------------------------------------------------
               WRITE (6,*) ' '
C
            END DO
C
         END DO
C
         DO IOBS = 2,NOBS
            OBS_ABS(IOBS) = DNRM2(3,OBS_TOT(1,IOBS),1)
         END DO
C
         WRITE (6,99018) TXT_OBS(1),OBS_TOT(0,IDOS),
     &                   (TXT_OBS(IOBS),OBS_ABS(IOBS),OBS_TOT(1:3,IOBS),
     &                   IOBS=2,3)
         WRITE (6,*) ' '
C
      END IF
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTORQUE
C                         torque calculations
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTORQUE
      IF ( THERMAL_VIBRA_FLUCT .AND. IREL.EQ.3 ) THEN
C
         WRITE (6,99014)
C
         DO IT = ITBOT,ITTOP
C
            IQ = IQAT(1,IT)
C
            WRITE (6,99015) IQ,IT,OBS_TX(0,IWFD,IT),
     &                      (OBS_TX(J,ITRQ,IT),J=1,3)
C
         END DO
C
      END IF
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTORQUE
C
C ---------------------------------add BAND densities to TOTAL densities
C ----------------------------- account for spin degeneracy for IREL <=1
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT_C:DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         IF ( FULLPOT ) THEN
C
            IRTOP = JRCRI(IM)
C
            DO J = 1,NDENS
               DO LM = 1,NLMFPMAX
                  DO IR = 1,IRTOP
                     RSQ = R(IR,IM)*R(IR,IM)
                     RHO2NS(IR,LM,IT,J)
     &                  = (RHO2NS(IR,LM,IT,J)+RHO2NSX(IR,LM,IT,J))*RSQ
                  END DO
               END DO
            END DO
C
C- convolute RHO2NS with shape functions to get spherical charge density
C
            DO IR = 1,JRMT(IM)
               RHOCHR(IR,IT) = RHO2NS(IR,1,IT,1)*SQRT_4PI
               RHOSPN(IR,IT) = RHO2NS(IR,1,IT,2)*SQRT_4PI
               RHOORB(IR,IT) = RHO2NS(IR,1,IT,3)*SQRT_4PI
            END DO
C
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               IF ( LM.LE.NLMFPMAX ) THEN
C
                  DO IRSF = 1,NRSFTOT(IM)
                     IR = IRSF + JRMT(IM)
C
                     RHOCHR(IR,IT) = RHOCHR(IR,IT) + FLMSF(IRSF,ISF,IM)
     &                               *RHO2NS(IR,LM,IT,1)
                     RHOSPN(IR,IT) = RHOSPN(IR,IT) + FLMSF(IRSF,ISF,IM)
     &                               *RHO2NS(IR,LM,IT,2)
                     RHOORB(IR,IT) = RHOORB(IR,IT) + FLMSF(IRSF,ISF,IM)
     &                               *RHO2NS(IR,LM,IT,3)
C
                  END DO
C
               END IF
            END DO
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         ELSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
            IRTOP = JRWS(IM)
C
            IF ( IREL.LE.1 ) THEN
               DO IR = 1,IRTOP
                  RHOCHR(IR,IT) = RHOCHR(IR,IT) + 2*RHOCHRX(IR,IT)
                  RHOSPN(IR,IT) = 0.0D0
                  RHOORB(IR,IT) = 0.0D0
               END DO
            ELSE
               DO IR = 1,IRTOP
                  RHOCHR(IR,IT) = RHOCHR(IR,IT) + RHOCHRX(IR,IT)
                  RHOSPN(IR,IT) = RHOSPN(IR,IT) + RHOSPNX(IR,IT)
                  RHOORB(IR,IT) = RHOORB(IR,IT) + RHOORBX(IR,IT)
               END DO
            END IF
C
         END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
C ------------------------------------------------- check BAND densities
C
         CALL CHRDNS_CHECK_RHO_BAND(IT,RHOCHR,OBS_TX,IDOS,IFLAG)
         IF ( NDENS.GE.2 ) CALL CHRDNS_CHECK_RHO_BAND(IT,RHOSPN,OBS_TX,
     &        ISMT,IFLAG)
         IF ( NDENS.EQ.3 ) CALL CHRDNS_CHECK_RHO_BAND(IT,RHOORB,OBS_TX,
     &        IOMT,IFLAG)
C
C -------------------------------- add CORE densities to TOTAL densities
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         IF ( FULLPOT ) THEN
C
            DO IR = 1,IRTOP
               RSQ = R(IR,IM)*R(IR,IM)/SQRT_4PI
               RHO2NS(IR,1,IT,1) = RHO2NS(IR,1,IT,1) + RHOCHRC(IR,IT)
     &                             *RSQ
               RHO2NS(IR,1,IT,2) = RHO2NS(IR,1,IT,2) + RHOSPNC(IR,IT)
     &                             *RSQ
            END DO
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         ELSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
            RHOCHR(1:IRTOP,IT) = RHOCHR(1:IRTOP,IT)
     &                           + RHOCHRC(1:IRTOP,IT)
            RHOSPN(1:IRTOP,IT) = RHOSPN(1:IRTOP,IT)
     &                           + RHOSPNC(1:IRTOP,IT)
C
            QEL(IT) = DDOT(IRTOP,RHOCHR(1,IT),1,R2DRDI_W_RADINT(1,IM),1)
C
            IF ( SCF_CHECK_SPLITSS ) WRITE (*,*) '<SCF>:   QEL ',QEL(IT)
C
         END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
      END DO LOOP_IT_C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( IFLAG.EQ.0 ) THEN
         WRITE (6,99019) ROUTINE(1:LEN_TRIM(ROUTINE)),TOL
      ELSE
         WRITE (6,99020) ROUTINE(1:LEN_TRIM(ROUTINE))
      END IF
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      IF ( FULLPOT ) THEN
C =================================================== CHECK_CORE_DENSITY
C               integrals for RHO2NS from core solver
C
         IF ( CHECK_CORE_DENSITY ) THEN
            IFLAG = 0
C
            WRITE (6,99024)
C
            DO IT = ITBOT,ITTOP
               IM = IMT(IT)
               IRTOP = JRCRI(IM)
C
               ISF = ISFLM(1,IM)
C
               INTEG_CHR = 0.0D0
               INTEG_SPN = 0.0D0
               DO IR = 1,JRMT(IM)
                  WR = R2DRDI_W_RADINT(IR,IM)
                  INTEG_CHR = INTEG_CHR + WR*RHOCHRC(IR,IT)
                  INTEG_SPN = INTEG_SPN + WR*RHOSPNC(IR,IT)
               END DO
C
               DO IRSF = 1,NRSFTOT(IM)
                  IR = IRSF + JRMT(IM)
                  WR = R2DRDI_W_RADINT(IR,IM)
                  WR = WR*FLMSF(IRSF,ISF,IM)/SQRT_4PI
                  INTEG_CHR = INTEG_CHR + WR*RHOCHRC(IR,IT)
                  INTEG_SPN = INTEG_SPN + WR*RHOSPNC(IR,IT)
               END DO
C
               WRITE (6,99025) IT,'CORE',' CHR:',INTEG_CHR
               WRITE (6,99025) IT,'CORE',' SPN:',INTEG_SPN
               IF ( ABS(Z(IT)-NVALT(IT)-INTEG_CHR).GT.1D-8 ) IFLAG = 1
C
            END DO
            IF ( IFLAG.EQ.1 )
     &            CALL STOP_MESSAGE(ROUTINE,'RHO core is not OK')
         END IF
C =================================================== CHECK_CORE_DENSITY
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
      END IF
C
      CALL CHRDNS_CHECK_OBS(NME,DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX)
C
C ======================================================================
      DEALLOCATE (ZGRA,ZFRA,JGRA,JFRA)
      IF ( ALLOCATED(JFLA) ) DEALLOCATE (JGLA,JFLA,ZFLA,ZGLA)
C ======================================================================
C
99001 FORMAT (/,' <CHRDNS>:  TOTNOS    ',F15.10,'  DQ        ',F15.10,/,
     &        ' <CHRDNS>:  TOTDOS    ',F15.10,'  TOTDOSX   ',F15.10)
99002 FORMAT ((' ',79('*'),/),/,' SPRKKR-run for: ',A,//,A)
99003 FORMAT (' CHARGE MISFIT     ',F12.8,' els.',/,
     &        ' E_F CORRECTION    ',F12.8,/,' NEW FERMI ENERGY  ',F12.8,
     &        '   extrapolated     ',/,:,' NEW FERMI ENERGY  ',F12.8,
     &        '   via Lloyd formula',/)
99004 FORMAT (/,' CHARGE MISFIT     ',F12.8,' els.',/,
     &        ' FERMI ENERGY      ',F12.8,/)
99005 FORMAT (/,I4,' E=',2F7.4,3X,'L=',I2,3X,'IT=',I4,2X,A,2X,A20,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE)  ',8X,F8.3,10X,F8.3,10X,
     &        F8.3,10X,F8.1,/,10X,2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1)
99006 FORMAT (/,I4,' E=',2F7.4,10X,'IT=',I4,2X,A,:,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE) crystal  ',F8.3,4X,F14.3,4X,
     &        F14.3,4X,F14.1,/,' TOTAL   crystal  ',F8.3,4X,F14.3,4X,
     &        F14.3,4X,F14.1)
C for testing:
C99010 FORMAT ('         DOS      NOS     P_spin   m_spin',
C     &        '    P_orb    m_orb    B_val      B_core',/,'  s ',2F9.4,
C     &        F10.4,F9.4,F10.5,F9.5,F16.6,' s  ',F8.2,:,/,'  p ',2F9.4,
C     &        F10.4,F9.4,F10.5,F9.5,F16.6,' ns ',F8.2,:,/,'  d ',2F9.4,
C     &        F10.4,F9.4,F10.5,F9.5,F16.6,' cor',F8.2)
99007 FORMAT ('         DOS      NOS     P_spin   m_spin',
     &        '    P_orb    m_orb    B_val      B_core',/,'  s ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' s  ',F8.2,:,/,'  p ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' ns ',F8.2,:,/,'  d ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' cor',F8.2)
99008 FORMAT ('  f ',2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,:,/,'  g ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2)
C for testing:
C99013 FORMAT (1X,A3,2F9.4,F10.4,F9.4,F10.5,F9.5,F16.6,:,' v+c',F8.2)
99009 FORMAT (1X,A3,2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,' v+c',F8.2)
99010 FORMAT (' ',79('-'),/,' TDOS/NOS ',2F8.3,:,' MUE-SPIN:',F8.3,
     &        '  MUE-ORB:',F8.3)
99011 FORMAT (' ',79('-'),/,' TOT',2F9.4,/,' E_band',F19.8,' [Ry]',/,
     &        ' ',79('='))
99012 FORMAT (' ',79('-'),/,' TOT',2F9.4,10X,F9.4,10X,F9.5,/,' E_band',
     &        F19.8,' [Ry]',/,' ',79('='))
99013 FORMAT (/,5X,'vector observables',//,4X,'IQ  IT CONC OBS',17X,
     &        'global frame',14X,'local frame',/)
99014 FORMAT (/,5X,'torque and Weiss field calculations',//,4X,
     &        'IQ  IT  WEISS TORQUE',/)
99015 FORMAT (I6,I5,F15.6,2X,3F15.6)
99016 FORMAT (I6,I4,F5.2,1X,A,F8.4,2X,3F8.4,2X,3F8.4,:,
     &        (/,16X,A,F8.4,2X,3F8.4,2X,3F8.4))
99017 FORMAT (I6,I4,F5.2,1X,A,F8.4,2X,3F8.4,:,(/,16X,A,F8.4,2X,3F8.4))
99018 FORMAT (3X,'TOT',10X,A,F8.4,:,/,(16X,A,F8.4,2X,3F8.4))
99019 FORMAT (/,10X,'integrals in <',A,'> agree within ',1PE9.1,/)
99020 FORMAT (/,1X,79('#'),/,25X,'integrals in <',A,'>  NOT OK',/,1X,
     &        79('#'),/)
99021 FORMAT (' IT ',I3,2X,A,2X,'DNS',F20.10,/,18X,'INT',F20.10,5X,
     &        'RAT',F15.10)
99022 FORMAT (//,' <CHRDNS>:     CHECK temporary densities ',/,15X,
     &        'for .NOT.EFCORRECT and IECURR = NETAB(IEPATH) ',/)
99023 FORMAT (' E_band',F19.8,' [Ry]')
99024 FORMAT (/,' <FPCHRDNS>:  check integrals for RHO2NS from core ',
     &        'solver',/)
99025 FORMAT (10X,'IT =',I3,2(2X,A),3X,2F16.12)
99026 FORMAT (/,1X,79('!'),/,10X,' delta DOS ',F10.12,/,10X,
     &        ' delta SMT ',F10.12,/,10X,' delta OMT ',F10.12,/,1X,
     &        79('!'),/)
99027 FORMAT (20F12.5)
99028 FORMAT ('# BUILDBOT: ',A,':  DOS NOS SMT OMT BND for IT =',I5,/,
     &        (1PE22.14))
99029 FORMAT (F10.2,F16.10)
C
      END
C*==chrdns_rel.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_REL(IT,IL,WE,IRTOP,DOBS_LMX,NME_RED,DRHOCHR,
     &                      DRHOORB,DRHOSPN,MEZZ,MEZJ,TMAT_LOC,ZGRA,
     &                      ZFRA,JGRA,JFRA)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <FPCHRDNS>                              *
C   *                                                                  *
C   *   dealing with FULLY RELATIVISTIC case                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NTMAX,DOBS_LTX,NCPLWFMAX,IKMCPLWF
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,NLMAX,NMUEMAX,NKMMAX,NMEMAX,
     &    NOBSMAX,NCPLWF,AME_G,IMKM_IKM
      USE MOD_CONSTANTS,ONLY:PI,C0
      IMPLICIT NONE
C*--CHRDNS_REL1372
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IL,IRTOP,IT,NME_RED
      COMPLEX*16 WE
      COMPLEX*16 DOBS_LMX(0:3,NOBSMAX,NLMAX,NMUEMAX),
     &           JFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TMAT_LOC(NKMMAX,NKMMAX),ZFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
      REAL*8 DRHOCHR(NRMAX),DRHOORB(NRMAX),DRHOSPN(NRMAX)
C
C Local variables
C
      INTEGER IA,IB,IKMA,IKMB,IKMCB(2),IME,IMKMA,IMKMB,IR,JA,JB,KAP1,
     &        KAP2,L,LAMA,LAMB,MJM05,MUE,NSOL
      INTEGER IKAPMUE
      REAL*8 MJ
      COMPLEX*16 WDS,WOF,WOG,WSF,WSG,WT_LOC,ZFJF,ZFZF,ZGJG,ZGZG
C
C*** End of declarations rewritten by SPAG
C
      L = IL - 1
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      KAP1 = -L - 1
      KAP2 = L
      IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      DOBS_LMX(:,:,:,:) = C0
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      MUE = 0
      DO MJM05 = -L - 1, + L
         MJ = DBLE(MJM05) + 0.5D0
         MUE = MUE + 1
C
C-----------------------------------------------------------------------
C           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
         IF ( ABS(MJ).GT.DBLE(L) ) THEN
            NSOL = 1
         ELSE
            NSOL = 2
            IKMCB(2) = IKAPMUE(KAP2,MJM05)
         END IF
         IKMCB(1) = IKAPMUE(KAP1,MJM05)
C-----------------------------------------------------------------------
C
         DO JB = 1,NSOL
C
            LAMB = IKMCPLWF(JB,IKMCB(JB))
C
            DO JA = 1,NSOL
C
               LAMA = IKMCPLWF(JA,IKMCB(JA))
C
               WT_LOC = -TMAT_LOC(LAMB,LAMA)/PI
C
               DO IME = 1,NME_RED
                  DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &               + WT_LOC*MEZZ(LAMB,LAMA,IT,IME)
               END DO
C
               LOOP_IB:DO IB = 1,NCPLWF(LAMB)
                  IKMB = IKMCPLWF(IB,LAMB)
                  IMKMB = IMKM_IKM(IKMB)
C
                  LOOP_IA:DO IA = 1,NCPLWF(LAMA)
                     IKMA = IKMCPLWF(IA,LAMA)
                     IMKMA = IMKM_IKM(IKMA)
C
                     WDS = WE*WT_LOC*AME_G(IKMB,IKMA,1,IDOS)
                     WSG = WE*WT_LOC*AME_G(IKMB,IKMA,2,ISMT)
                     WOG = WE*WT_LOC*AME_G(IKMB,IKMA,2,IOMT)
                     IF ( IKMA.EQ.IKMB ) THEN
                        WSF = WE*WT_LOC*AME_G(IMKMB,IMKMA,2,ISMT)
                        WOF = WE*WT_LOC*AME_G(IMKMB,IMKMA,2,IOMT)
                     ELSE
                        WSF = C0
                        WOF = C0
                     END IF
C
                     DO IR = 1,IRTOP
                        ZGZG = ZGRA(IR,IB,LAMB)*ZGRA(IR,IA,LAMA)
                        ZFZF = ZFRA(IR,IB,LAMB)*ZFRA(IR,IA,LAMA)
C
                        DRHOCHR(IR) = DRHOCHR(IR)
     &                                + DIMAG(WDS*ZGZG+WDS*ZFZF)
                        DRHOSPN(IR) = DRHOSPN(IR)
     &                                + DIMAG(WSG*ZGZG-WSF*ZFZF)
                        DRHOORB(IR) = DRHOORB(IR)
     &                                + DIMAG(WOG*ZGZG-WOF*ZFZF)
                     END DO
                  END DO LOOP_IA
               END DO LOOP_IB
C
C---------------- no irregular contributions to the backscattering terms
C
               IF ( LAMB.EQ.LAMA .AND. .NOT.SPLITSS ) THEN
C
                  DO IME = 1,NME_RED
                     DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &                  + MEZJ(LAMB,LAMA,IT,IME)/PI
                  END DO
C
                  DO IB = 1,NSOL
                     IKMB = IKMCPLWF(IB,LAMB)
                     IMKMB = IMKM_IKM(IKMB)
C
                     DO IA = 1,NSOL
                        IKMA = IKMCPLWF(IA,LAMA)
                        IMKMA = IMKM_IKM(IKMA)
C
                        WDS = WE*AME_G(IKMB,IKMA,1,IDOS)/PI
                        WSG = WE*AME_G(IKMB,IKMA,2,ISMT)/PI
                        WOG = WE*AME_G(IKMB,IKMA,2,IOMT)/PI
C
                        IF ( IKMA.EQ.IKMB ) THEN
                           WSF = WE*AME_G(IMKMB,IMKMA,2,ISMT)/PI
                           WOF = WE*AME_G(IMKMB,IMKMA,2,IOMT)/PI
                        ELSE
                           WSF = C0
                           WOF = C0
                        END IF
C
                        DO IR = 1,IRTOP
                           ZGJG = ZGRA(IR,IB,LAMB)*JGRA(IR,IA,LAMA)
                           ZFJF = ZFRA(IR,IB,LAMB)*JFRA(IR,IA,LAMA)
C
                           DRHOCHR(IR) = DRHOCHR(IR)
     &                        + DIMAG(WDS*ZGJG+WDS*ZFJF)
                           DRHOSPN(IR) = DRHOSPN(IR)
     &                        + DIMAG(WSG*ZGJG-WSF*ZFJF)
                           DRHOORB(IR) = DRHOORB(IR)
     &                        + DIMAG(WOG*ZGJG-WOF*ZFJF)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
         DO IME = 1,NME_RED
            DOBS_LTX(0,IME,IL,IT) = DOBS_LTX(0,IME,IL,IT)
     &                              + DOBS_LMX(0,IME,IL,MUE)
         END DO
C
      END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
      END
C*==chrdns_sra.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_SRA(IT,IL,WE,IRTOP,DOBS_LMX,NME_RED,DRHOCHR,
     &                      DRHOSPN,MEZZ,MEZJ,TMAT,ZGRA,JGRA)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <CHRDNS>                                *
C   *                                                                  *
C   *   dealing with SCALAR RELATIVISTIC case                          *
C   *                                                                  *
C   *  NOTE: for IREL <= 2 (non- and scalar relativistic) the          *
C   *        minor component has no meaning and is not written/read    *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NTMAX,DOBS_LTX,NCPLWFMAX
      USE MOD_ANGMOM,ONLY:IHFF,NLM,NSPIN,NKMMAX,NL,NLMAX,NMEMAX,NMUEMAX,
     &    NOBSMAX
      USE MOD_CONSTANTS,ONLY:PI,C0
      IMPLICIT NONE
C*--CHRDNS_SRA1561
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IL,IRTOP,IT,NME_RED
      COMPLEX*16 WE
      COMPLEX*16 DOBS_LMX(0:3,NOBSMAX,NLMAX,NMUEMAX),
     &           JGRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),TMAT(NKMMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
      REAL*8 DRHOCHR(NRMAX),DRHOSPN(NRMAX)
C
C Local variables
C
      REAL*8 CHR,SPNWGT
      INTEGER ILS,IME,IR,IS,L,LMS,LMSOFF,ML,MUE
      COMPLEX*16 WT,WZJ,WZZ,ZGJG,ZGZG
C
C*** End of declarations rewritten by SPAG
C
      DOBS_LMX(:,:,:,:) = C0
C
      L = IL - 1
      MUE = 0
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      DO ML = -L,L
         MUE = MUE + 1
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         DO IS = 1,NSPIN
C
            ILS = NL*(IS-1) + IL
C
            LMSOFF = NLM*(IS-1)
            LMS = LMSOFF + L*L + L + ML + 1
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
C
            WT = -TMAT(LMS,LMS)/PI
C
            DO IME = 1,NME_RED
               IF ( IME.NE.IHFF ) THEN
                  DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &               + WT*MEZZ(LMS,LMS,IT,IME)
               ELSE
                  DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &               + WT*MEZZ(LMS,LMS,IT,IME)*SPNWGT
               END IF
            END DO
C
            IF ( SPLITSS ) THEN
C
               WZZ = WE*WT
               DO IR = 1,IRTOP
                  ZGZG = ZGRA(IR,1,ILS)*ZGRA(IR,1,ILS)
C
                  CHR = DIMAG(WZZ*ZGZG)
C
                  DRHOCHR(IR) = DRHOCHR(IR) + CHR
                  DRHOSPN(IR) = DRHOSPN(IR) + CHR*SPNWGT
               END DO
C
            ELSE
C
               DO IME = 1,NME_RED
                  IF ( IME.NE.IHFF ) THEN
                     DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &                  + MEZJ(LMS,LMS,IT,IME)/PI
                  ELSE
                     DOBS_LMX(0,IME,IL,MUE) = DOBS_LMX(0,IME,IL,MUE)
     &                  + MEZJ(LMS,LMS,IT,IME)/PI*SPNWGT
                  END IF
               END DO
C
               WZZ = WE*WT
               WZJ = WE/PI
               DO IR = 1,IRTOP
                  ZGZG = ZGRA(IR,1,ILS)*ZGRA(IR,1,ILS)
                  ZGJG = ZGRA(IR,1,ILS)*JGRA(IR,1,ILS)
C
                  CHR = DIMAG(WZZ*ZGZG+WZJ*ZGJG)
C
                  DRHOCHR(IR) = DRHOCHR(IR) + CHR
                  DRHOSPN(IR) = DRHOSPN(IR) + CHR*SPNWGT
               END DO
C
            END IF
C
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         DO IME = 1,NME_RED
            DOBS_LTX(0,IME,IL,IT) = DOBS_LTX(0,IME,IL,IT)
     &                              + DOBS_LMX(0,IME,IL,MUE)
         END DO
C
      END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
      END
C*==chrdns_dip_mnt_rel.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_DIP_MNT_REL(IT,IM,WE,IRTOP,TMAT_LOC,ZGRA,ZFRA,
     &                              JGRA,JFRA,CMNTTX)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <CHRDNS>                                *
C   *                                                                  *
C   *   dealing with FULLY RELATIVISTIC case                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_RMESH,ONLY:R,R2DRDI_W_RADINT,NRMAX
      USE MOD_ANGMOM,ONLY:A_Y1M,NKMMAX,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_TYPES,ONLY:NKM_T,NTMAX,NLMFPMAX,NCPLWFMAX,IKMCPLWF
      IMPLICIT NONE
C*--CHRDNS_DIP_MNT_REL1695
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,IT
      COMPLEX*16 WE
      REAL*8 CMNTTX(NLMFPMAX,NTMAX)
      COMPLEX*16 JFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGRA(NRMAX,NCPLWFMAX,NKMMAX),TMAT_LOC(NKMMAX,NKMMAX),
     &           ZFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      LOGICAL CNON0
      INTEGER ICWF,IKM,IR,ISOL,JKM,JSOL,JWF,LMQ,MQ
      COMPLEX*16 JFZF,JGZG,MIRR(-1:+1),MREG(-1:+1),RIRR,RREG,WIRR,WREG,
     &           ZFZF,ZGZG
      REAL*8 WR(NRMAX)
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C
      DO IR = 1,IRTOP
         WR(IR) = R(IR,IM)*R2DRDI_W_RADINT(IR,IM)
      END DO
C
C=======================================================================
      WIRR = WE/PI
C
      DO ISOL = 1,NKM_T(IT)
         DO JSOL = 1,NKM_T(IT)
C
            WREG = -WE*TMAT_LOC(JSOL,ISOL)/PI
C
            MREG(-1:1) = 0D0
            MIRR(-1:1) = 0D0
C
C-----------------------------------------------------------------------
            DO ICWF = 1,NCPLWF(ISOL)
               IKM = IKMCPLWF(ICWF,ISOL)
C
               DO JWF = 1,NCPLWF(JSOL)
                  JKM = IKMCPLWF(JWF,JSOL)
C
                  DO MQ = -1, + 1
                     IF ( CNON0(A_Y1M(IKM,JKM,MQ)) ) GOTO 5
                  END DO
                  CYCLE
C-----------------------------------------------------------------------
C
 5                CONTINUE
                  RREG = 0D0
                  RIRR = 0D0
C
                  IF ( ISOL.EQ.JSOL .AND. .NOT.SPLITSS ) THEN
C
                     DO IR = 1,IRTOP
                        ZGZG = ZGRA(IR,ICWF,ISOL)*ZGRA(IR,JWF,JSOL)
                        ZFZF = ZFRA(IR,ICWF,ISOL)*ZFRA(IR,JWF,JSOL)
                        RREG = RREG + (ZGZG+ZFZF)*WR(IR)
C
                        JGZG = JGRA(IR,ICWF,ISOL)*ZGRA(IR,JWF,JSOL)
                        JFZF = JFRA(IR,ICWF,ISOL)*ZFRA(IR,JWF,JSOL)
                        RIRR = RIRR + (JGZG+JFZF)*WR(IR)
                     END DO
C
                     DO MQ = -1, + 1
                        MREG(MQ) = MREG(MQ) + RREG*A_Y1M(IKM,JKM,MQ)
                        MIRR(MQ) = MIRR(MQ) + RIRR*A_Y1M(IKM,JKM,MQ)
                     END DO
C
                  ELSE
C
                     DO IR = 1,IRTOP
                        ZGZG = ZGRA(IR,ICWF,ISOL)*ZGRA(IR,JWF,JSOL)
                        ZFZF = ZFRA(IR,ICWF,ISOL)*ZFRA(IR,JWF,JSOL)
                        RREG = RREG + (ZGZG+ZFZF)*WR(IR)
                     END DO
C
                     DO MQ = -1, + 1
                        MREG(MQ) = MREG(MQ) + RREG*A_Y1M(IKM,JKM,MQ)
                     END DO
C
                  END IF
C
               END DO
            END DO
C-----------------------------------------------------------------------
C
            IF ( ISOL.EQ.JSOL .AND. .NOT.SPLITSS ) THEN
C
               DO MQ = -1, + 1
                  LMQ = 3 + MQ
                  CMNTTX(LMQ,IT) = CMNTTX(LMQ,IT)
     &                             + DIMAG(MREG(MQ)*WREG+MIRR(MQ)*WIRR)
               END DO
C
            ELSE
C
               DO MQ = -1, + 1
                  LMQ = 3 + MQ
                  CMNTTX(LMQ,IT) = CMNTTX(LMQ,IT) + DIMAG(MREG(MQ)*WREG)
               END DO
C
            END IF
C
         END DO
      END DO
C=======================================================================
C
      END
C*==chrdns_dip_mnt_sra.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_DIP_MNT_SRA(IT,IM,WE,IRTOP,TMAT,ZGRA,JGRA,WSPIN,
     &                              CMNTTX)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <CHRDNS>                                *
C   *                                                                  *
C   *   dealing with SCALAR RELATIVISTIC case                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_RMESH,ONLY:NRMAX,R,R2DRDI_W_RADINT
      USE MOD_TYPES,ONLY:NLT,NTMAX,NLMFPMAX,NCPLWFMAX
      USE MOD_ANGMOM,ONLY:A_Y1M,NLM,L_LM,NKMMAX,NL,NSPIN
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--CHRDNS_DIP_MNT_SRA1837
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,IT
      COMPLEX*16 WE
      REAL*8 WSPIN
      REAL*8 CMNTTX(NLMFPMAX,NTMAX)
      COMPLEX*16 JGRA(NRMAX,NCPLWFMAX,NKMMAX),TMAT(NKMMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      LOGICAL CNON0
      INTEGER IL1,IL2,ILS1,ILS2,IR,IS,LM1,LM2,LMQ,LMS1,LMS2,LMSOFF,MQ,
     &        NLM_LOC
      COMPLEX*16 JGZG,MIRR(-1:+1),MREG(-1:+1),RIRR,RREG,WIRR,WREG,ZGZG
      REAL*8 WR(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WR
C
      ALLOCATE (WR(NRMAX))
C
C=======================================================================
C
      DO IR = 1,IRTOP
         WR(IR) = R(IR,IM)*R2DRDI_W_RADINT(IR,IM)
      END DO
C
C=======================================================================
C
      NLM_LOC = NLT(IT)*NLT(IT)
C
      WIRR = WSPIN*WE/PI
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO IS = 1,NSPIN
         LMSOFF = NLM*(IS-1)
C
C=======================================================================
         DO LM1 = 1,NLM_LOC
            IL1 = L_LM(LM1) + 1
            ILS1 = NL*(IS-1) + IL1
            LMS1 = LMSOFF + LM1
C
            DO LM2 = 1,NLM_LOC
               IL2 = L_LM(LM2) + 1
               ILS2 = NL*(IS-1) + IL2
               LMS2 = LMSOFF + LM2
C
               WREG = -WSPIN*WE*TMAT(LMS1,LMS2)/PI
C
               MREG(-1:1) = 0D0
               MIRR(-1:1) = 0D0
C
               DO MQ = -1, + 1
                  IF ( CNON0(A_Y1M(LM1,LM2,MQ)) ) GOTO 10
               END DO
               CYCLE
C-----------------------------------------------------------------------
C
 10            CONTINUE
               RREG = 0D0
               RIRR = 0D0
C
               IF ( LM1.EQ.LM2 .AND. .NOT.SPLITSS ) THEN
C
                  DO IR = 1,IRTOP
                     ZGZG = ZGRA(IR,1,ILS1)*ZGRA(IR,1,ILS2)
                     RREG = RREG + ZGZG*WR(IR)
C
                     JGZG = JGRA(IR,1,ILS1)*ZGRA(IR,1,ILS2)
                     RIRR = RIRR + JGZG*WR(IR)
                  END DO
C
                  DO MQ = -1, + 1
                     MREG(MQ) = MREG(MQ) + RREG*A_Y1M(LM1,LM2,MQ)
                     MIRR(MQ) = MIRR(MQ) + RIRR*A_Y1M(LM1,LM2,MQ)
                  END DO
C
               ELSE
C
                  DO IR = 1,IRTOP
                     ZGZG = ZGRA(IR,1,ILS1)*ZGRA(IR,1,ILS2)
                     RREG = RREG + ZGZG*WR(IR)
                  END DO
C
                  DO MQ = -1, + 1
                     MREG(MQ) = MREG(MQ) + RREG*A_Y1M(LM1,LM2,MQ)
                  END DO
C
               END IF
C
C-----------------------------------------------------------------------
C
               IF ( LM1.EQ.LM2 ) THEN
C
                  DO MQ = -1, + 1
                     LMQ = 3 + MQ
                     CMNTTX(LMQ,IT) = CMNTTX(LMQ,IT)
     &                                + DIMAG(MREG(MQ)*WREG+MIRR(MQ)
     &                                *WIRR)
                  END DO
C
               ELSE
C
                  DO MQ = -1, + 1
                     LMQ = 3 + MQ
                     CMNTTX(LMQ,IT) = CMNTTX(LMQ,IT)
     &                                + DIMAG(MREG(MQ)*WREG)
                  END DO
C
               END IF
C
            END DO
         END DO
C
C=======================================================================
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      END
C*==chrdns_check_rho_band.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_CHECK_RHO_BAND(IT,RHO,OBS_TX,IKEY,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,ISMT,IOMT
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS,JRCRI,R2DRDI_W_RADINT,
     &    DRDI_W_RADINT
      USE MOD_TYPES,ONLY:IMT,NTMAX
      IMPLICIT NONE
C*--CHRDNS_CHECK_RHO_BAND1987
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
      INTEGER IFLAG,IKEY,IT
      REAL*8 OBS_TX(0:3,NOBSMAX,NTMAX),RHO(NRMAX,NTMAX)
C
C Local variables
C
      REAL*8 DDOT
      INTEGER IM,IRTOP
      REAL*8 INTEGRAL
      CHARACTER*7 TXT7
C
C*** End of declarations rewritten by SPAG
C
      IF ( IKEY.EQ.IDOS ) THEN
         TXT7 = 'charge '
      ELSE IF ( IKEY.EQ.ISMT ) THEN
         IF ( IREL.LT.2 ) RETURN
         TXT7 = 'spin   '
      ELSE IF ( IKEY.EQ.IOMT ) THEN
         IF ( IREL.LT.3 ) RETURN
         TXT7 = 'orbital'
      ELSE
         RETURN
      END IF
C
      IM = IMT(IT)
C
      IF ( FULLPOT ) THEN
C
         IRTOP = JRCRI(IM)
C
         INTEGRAL = DDOT(IRTOP,RHO(1,IT),1,DRDI_W_RADINT(1,IM),1)
C
      ELSE
C
         IRTOP = JRWS(IM)
C
         INTEGRAL = DDOT(IRTOP,RHO(1,IT),1,R2DRDI_W_RADINT(1,IM),1)
C
      END IF
C
      IF ( IREL.LE.1 ) INTEGRAL = INTEGRAL/2D0
C
      IF ( ABS(INTEGRAL-OBS_TX(0,IKEY,IT)).GT.TOL ) THEN
         IFLAG = IFLAG + 1
         IF ( ABS(INTEGRAL).GT.TOL ) THEN
            WRITE (6,99001) IT,TXT7,OBS_TX(0,IKEY,IT),INTEGRAL,
     &                      OBS_TX(0,IKEY,IT)/INTEGRAL
         ELSE
            WRITE (6,99001) IT,TXT7,OBS_TX(0,IKEY,IT),INTEGRAL
         END IF
      END IF
C
99001 FORMAT (' IT ',I3,2X,A,2X,'DNS',F20.10,/,18X,'INT',F20.10,5X,
     &        'RAT',F15.10)
C
      END
C*==chrdns_check_obs.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_CHECK_OBS(NME,DOBS_LT,DOBS_T,OBS_LT,OBS_T)
C   ********************************************************************
C   *                                                                  *
C   *  check consistency of OBS-arrays                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NLMAX,NPOL,NOBSMAX
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP,NLT
      IMPLICIT NONE
C*--CHRDNS_CHECK_OBS2079
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHRDNS_CHECK_OBS')
C
C Dummy arguments
C
      INTEGER NME
      COMPLEX*16 DOBS_LT(0:3,NOBSMAX,NLMAX,NTMAX),
     &           DOBS_T(0:3,NOBSMAX,NTMAX)
      REAL*8 OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX),OBS_T(0:3,NOBSMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER IFLAG,IL,IME,IPOL,IT
      REAL*8 RSUM
C
C*** End of declarations rewritten by SPAG
C
      IFLAG = 0
C
      IF ( NME.GT.1 ) RETURN
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         DO IPOL = 0,NPOL
            IF ( IPOL.GT.0 ) CYCLE
            DO IME = 1,NME
C
               RSUM = 0D0
               CSUM = 0D0
               DO IL = 1,NLT(IT)
                  RSUM = RSUM + OBS_LT(IPOL,IME,IL,IT)
                  CSUM = CSUM + DOBS_LT(IPOL,IME,IL,IT)
               END DO
C
               IF ( ABS(RSUM-OBS_T(IPOL,IME,IT)).GT.1D-10 ) THEN
                  IFLAG = 1
                  WRITE (6,99001) IT,IPOL,IME,RSUM,OBS_T(IPOL,IME,IT)
               END IF
               IF ( ABS(CSUM-DOBS_T(IPOL,IME,IT)).GT.1D-10 ) THEN
                  IFLAG = 1
                  WRITE (6,99002) IT,IPOL,IME,CSUM,DOBS_T(IPOL,IME,IT)
               END IF
C
            END DO
         END DO
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( IFLAG.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'l-sum inconsistent')
C
99001 FORMAT (/,5X,'l-sum for  OBS and IT=',I3,'  IPOL=',I2,'  IME=',I2,
     &        2E20.12)
99002 FORMAT (/,5X,'l-sum for DOBS and IT=',I3,'  IPOL=',I2,'  IME=',I2,
     &        2E20.12,2X,2E20.12)
      END
