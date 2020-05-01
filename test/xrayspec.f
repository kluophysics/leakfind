C*==xrayspec.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XRAYSPEC(CALCINT,GETIRRSOL,ETAB0,GCOR,FCOR,ECOR,SZCOR,
     &                    KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,
     &                    BCOR,BCORS,NCSTMAX,NPOLMAX)
C   ********************************************************************
C   *                                                                  *
C   *   TASK = RXAS/RXES                                               *
C   *                                                                  *
C   *       calculation of  x-ray-absorption/emission spectra for      *
C   *       circularly and linearly polarized radiation                *
C   *                                                                  *
C   *   TASK = XMO                                                     *
C   *                                                                  *
C   *       calculation of full optical conductivity tensor SIGMA      *
C   *                                                                  *
C   *   TASK = XRS                                                     *
C   *                                                                  *
C   *       calculation of X-ray resonant scattering spectrum          *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *                 INPUT-FILES                                      *
C   *                 ***********                                      *
C   *                                                                  *
C   *    FILE   SUBROUTINE  CONTENT                                    *
C   *                                                                  *
C   *      9    <XRAYSPEC>  E, TAU(KAP,MUE)                            *
C   *      5    <XRAYSPEC>  parameter: IPRINT, EFERMI,     ...         *
C   *      4    <POTFIT>    potential                                  *
C   *                                                                  *
C   *                 OUTPUT-CONTROL                                   *
C   *                 **************                                   *
C   *                                                                  *
C   *   IPRINT   SUBROUTINE  FILE  ADDITIONAL OUTPUT                   *
C   *                                                                  *
C   *     0      <XRAYSPEC>   7    abs-rate                            *
C   *            <XRAYSPEC>   6    some control output                 *
C   *   >=1      <XRAYME>     8    radial-matrix elements              *
C   *   >=2      <XRAYME>     6    control-output                      *
C   *            <XRAYME>     6    radial-matrix elements              *
C   *   >=4      <AMEADA>     6    matrices  A+, A-, A1+, A1- ....     *
C   *   >=5      <XRAYSPEC>   6    matrices  TAULMS, TAU (KAP,MUE)     *
C   *                                                                  *
C   *                                                                  *
C   *  E0      electron charge in esu                                  *
C   *  CCGS    speed of light in cm/s                                  *
C   *  M0      electron mass   in g                                    *
C   *  HBAR    Planck - const. in erg*s                                *
C   *  EVORY   Ry to eV                                                *
C   *  ERGOEV  eV to erg                                               *
C   *  ERGORY  Ry to erg                                               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EMAX,WETAB
      USE MOD_RMESH,ONLY:R,NRMAX,FULLPOT,JRWS,RWS,R2DRDI_W_RADINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,MEZJ,MEZZ,MSSQ,MSST,SSST,TAUQ,TAUT,
     &    TSST,WKM1,MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4,IPIVKM
      USE MOD_SITES,ONLY:NQ,IQAT,ITOQ,MROTQ,DROTQ
      USE MOD_CALCMODE,ONLY:SELFENERGY,ORBPOL,DMFT,GF_CONV_RH,SOLVER_FP,
     &    MOMENTS_ROTATED,THERMAL_VIBRA_FLUCT
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,LCXRAY,NCXRAY,Z,BT,VT,IMT,NT,
     &    LTXT_T,TXT_T,CTL
      USE MOD_FILES,ONLY:LRECREAL8,LSYSTEM,SYSTEM,LDATSET,DATSET,WRTAU,
     &    WRTAUMQ,IPRINT,IFILCBWF,IFILCORWF,IFILBUILDBOT,WRBUILDBOT,
     &    FOUND_SECTION,N_FOUND,S10DUMMY
      USE MOD_CONSTANTS,ONLY:E0_CGS,A0_CGS,HBAR_CGS,M0_CGS,C_CGS,C0,CI,
     &    C1,PI,EV_ERG,RY_EV,SQRT_2
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA,KSELF
      USE MOD_THERMAL,ONLY:NFT,FMAT_FT,MFMAT_FT,NVIBFLU,IQ_AVFT,NVIBRA,
     &    NFLUCT,I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT,UMAT_VT
      IMPLICIT NONE
C*--XRAYSPEC73
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='XRAYSPEC')
      REAL*8 DELTA
      PARAMETER (DELTA=1D-8)
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER ITXRAY,NCSTMAX,NPOLMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      COMPLEX*16 ETAB0(NEMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
C
C Local variables
C
      COMPLEX*16 ABSRATE(:,:),ABSRTA(:,:),ABSRTB(:,:),C1MAT(:,:),CECOR,
     &           CSCL,CSUM1,CSUM2,ECOR_KM(:),ERYD,ERYDCHK,ETERM1,ETERM2,
     &           MCV(3),MCVX(3),MECV(:,:,:),MECVX(:,:,:),MEVC(:,:,:),
     &           MEVXC(:,:,:),MIRR_2_FT(:,:,:,:),
     &           MIRR_2_NONDIP(:,:,:,:,:),MIRR_2_TMP(:,:,:,:),
     &           MIRR_3_FT(:,:,:,:),MIRR_3_NONDIP(:,:,:,:,:),
     &           MIRR_3_NONDIPX(:,:,:,:,:),MIRR_3_TMP(:,:,:,:),
     &           MIRR_4_FT(:,:,:,:),MIRR_4_NONDIP(:,:,:,:,:),
     &           MIRR_4_TMP(:,:,:,:),MSS_FT(:,:),MSS_VQ(:,:),MVC(3),
     &           MVXC(3),MZAZB_FT(:,:,:),MZAZB_NONDIP(:,:,:,:),
     &           MZAZB_NONDIPX(:,:,:,:),MZAZB_TMP(:,:,:),MZBZA_FT(:,:,:)
     &           ,MZBZA_NONDIP(:,:,:,:),MZBZA_NONDIPX(:,:,:,:),
     &           MZBZA_TMP(:,:,:),P,RATA,RATB,SIG(:,:,:,:),TAU12,
     &           TAU_VFT(:,:),TAU_VQ(:,:),TMAT_LOC(:,:),TMIN,TMINX,TPLS,
     &           TSS_FT(:,:),WE
      REAL*8 C,CPACHNG,DBDR(:),DELE,DVDR(:),EILOW,EIMAG,EMAX3,EPHOT,
     &       ESIG(:),MROT(3,3),MROT_FRAME(3,3),NIRR,NREG,PHASE_FACTOR,
     &       PHAS_POL(3),SCLRATE,SCLSIG,SCLSIG1,SCLSIGALT,SCLXRS,SIGTOL,
     &       TAUCORE(:),TAUCOREK(2),VUC,WTERM1,WTERM2,XNORM(2),XSCL
      LOGICAL CHANGE_FRAME,FOUND,TAUCOREWARNING,XES,XMO_DYNAMICAL
      CHARACTER*1 CHPOL1(5),SHELL(5)
      CHARACTER*2 CHSIG(3,3),CL
      CHARACTER*80 FILNAM,SPEC
      INTEGER I,IA_ERR,ICPAFLAG,ICST,IDIPOL,IE,IELOOP,IENERG,IESIG,IFIL,
     &        IFLUCT,IFT,IKM,IKM_CST,IM,IOL,IPOL,IPOLREV(3),IPOLXC,
     &        IPOL_NEW_OLD(3),IQ,ISIG,IT,IVFT,IVIBRA,IVT,IX1,IX2,J,JT,K,
     &        KTYP,KUNIT,LAM1,LAM2,LFN,LSP,LSUBSH(0:4),M,N,NCST,NE,NE3,
     &        NELOOP,NESIG,NPOL,NQ9,NT9,NTXRAY,NTXRSGRP,TSELECT
      CHARACTER*3 MEFORM,SUBSH(0:4),SUBSHP(0:4),TAU_FILE_FMT
      CHARACTER*10 OUTPUT,TASK
      CHARACTER*7 STRWARN
      REAL*4 WCOREHOLE
C
C*** End of declarations rewritten by SPAG
C
      DATA SHELL/'K','L','M','N','O'/,TAU_FILE_FMT/'NEW'/
      DATA LSUBSH/1,3,3,3,3/
      DATA SUBSH/'1  ','2,3','4,5','6,7','8,9'/
      DATA SUBSHP/'1  ','23 ','45 ','67 ','89 '/
      DATA CHPOL1/'+','-','z','x','y'/
      DATA CHSIG/'xx','xy','xz','yx','yy','yz','zx','zy','zz'/
C
      DATA TAUCOREK/0D0,0D0/,TAUCOREWARNING/.FALSE./
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE DBDR,MECV,MEVC,DVDR,ESIG,ECOR_KM
      ALLOCATABLE MECVX,MEVXC,ABSRTA,ABSRTB
      ALLOCATABLE SIG,ABSRATE,TAUCORE
      ALLOCATABLE C1MAT,TMAT_LOC
      ALLOCATABLE MZAZB_TMP,MZBZA_TMP
      ALLOCATABLE MIRR_2_TMP,MIRR_3_TMP,MIRR_4_TMP
      ALLOCATABLE MZAZB_NONDIP,MZBZA_NONDIP
      ALLOCATABLE MIRR_2_NONDIP,MIRR_3_NONDIP,MIRR_4_NONDIP
      ALLOCATABLE MZAZB_NONDIPX,MZBZA_NONDIPX
      ALLOCATABLE MIRR_3_NONDIPX
      ALLOCATABLE TSS_FT,MZAZB_FT,MZBZA_FT,TAU_VQ,MSS_VQ
      ALLOCATABLE MSS_FT,TAU_VFT
      ALLOCATABLE MIRR_2_FT,MIRR_3_FT,MIRR_4_FT
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (MECV(NKMMAX,NCSTMAX,3),ECOR_KM(NKMMAX))
      ALLOCATE (MEVC(NKMMAX,NCSTMAX,3),DBDR(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEVC')
C
      ALLOCATE (MZAZB_NONDIPX(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MZBZA_NONDIPX(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MIRR_3_NONDIPX(NKMMAX,NKMMAX,3,3,3))
C##########################################################################
      ALLOCATE (ABSRTA(NCSTMAX,NPOLMAX))
      ALLOCATE (ABSRTB(NCSTMAX,NPOLMAX))
      ALLOCATE (ABSRATE(NCSTMAX,NPOLMAX))
      ALLOCATE (MZAZB(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MZBZA(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MIRR_2(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_3(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_4(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MZAZB_NONDIP(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MZBZA_NONDIP(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MIRR_2_NONDIP(NKMMAX,NKMMAX,3,3,3))
      ALLOCATE (MIRR_3_NONDIP(NKMMAX,NKMMAX,3,3,3))
      ALLOCATE (MIRR_4_NONDIP(NKMMAX,NKMMAX,3,3,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZAZB')
C
      ALLOCATE (MSS_VQ(NKMMAX,NKMMAX),TAU_VQ(NKMMAX,NKMMAX))
      ALLOCATE (MSS_FT(NKMMAX,NKMMAX),TAU_VFT(NKMMAX,NKMMAX))
      ALLOCATE (TSS_FT(NKMMAX,NKMMAX))
      ALLOCATE (MZAZB_FT(NKMMAX,NKMMAX,3))
      ALLOCATE (MZBZA_FT(NKMMAX,NKMMAX,3))
C
      ALLOCATE (MIRR_2_FT(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MIRR_3_FT(NKMMAX,NKMMAX,3,3))
      ALLOCATE (MIRR_4_FT(NKMMAX,NKMMAX,3,3))
C
      ALLOCATE (TMAT_LOC(NKMMAX,NKMMAX))
      ALLOCATE (MECVX(NKMMAX,NCSTMAX,3))
      ALLOCATE (MEVXC(NKMMAX,NCSTMAX,3),DVDR(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEVXC')
C
      ALLOCATE (TAUCORE(NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUCORE')
C
      IPOLREV(1:3) = (/3,2,1/)
C---- IPOL_NEW_OLD:        print out spectra in OLD sequence (+),(-),(z)
      IPOL_NEW_OLD(1:3) = (/3,1,2/)
      PHAS_POL(1:3) = (/-1D0,+1D0,+1D0/)
C
      WRITE (6,99004)
      IF ( NPOLMAX.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'NPOLMAX <> 3')
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
      SIGTOL = 1D-10
C
      WRTAU = .FALSE.
      WRTAUMQ = .FALSE.
C
      OUTPUT = 'BARN      '
      MEFORM = 'ADA'
C
      NE3 = 40
      EMAX3 = 8D0
      TSELECT = 0
      WTERM1 = 1D0
      WTERM2 = 1D0
C
      CALL INIT_MROT_FRAME(CHANGE_FRAME,MROT_FRAME)
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      IF ( CHANGE_FRAME .OR. MOMENTS_ROTATED ) THEN
C
         ALLOCATE (C1MAT(NKMMAX,NKMMAX))
         C1MAT(:,:) = C0
         DO IKM = 1,NKM
            C1MAT(IKM,IKM) = C1
         END DO
C
         ALLOCATE (MZAZB_TMP(NKMMAX,NKMMAX,3))
         ALLOCATE (MZBZA_TMP(NKMMAX,NKMMAX,3))
         ALLOCATE (MIRR_2_TMP(NKMMAX,NKMMAX,3,3))
         ALLOCATE (MIRR_3_TMP(NKMMAX,NKMMAX,3,3))
         ALLOCATE (MIRR_4_TMP(NKMMAX,NKMMAX,3,3))
C
      END IF
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      NPOL = 3
C
      IF ( FOUND_SECTION ) THEN
C
         TASK = 'RXAS      '
C
         CALL SECTION_FIND_KEYWORD('XES',XES)
         IF ( XES ) THEN
            NPOL = 3
            TASK = 'RXES      '
         ELSE
            NPOL = 2
         END IF
         NPOL = 3
C
         CALL SECTION_SET_STRING('ME',MEFORM,'9999',0)
         CALL SECTION_SET_STRING('OUTPUT',OUTPUT,'9999',0)
C
         CALL SECTION_FIND_KEYWORD('XMO',FOUND)
         IF ( FOUND ) THEN
            TASK = 'XMO       '
         ELSE
            CALL SECTION_FIND_KEYWORD('XRS',FOUND)
            IF ( FOUND ) TASK = 'XRS       '
         END IF
C
         IF ( TASK.EQ.'RXAS      ' .OR. TASK.EQ.'XMO       ' ) THEN
            CALL SECTION_FIND_KEYWORD('DYNAMICAL',XMO_DYNAMICAL)
            IF ( XMO_DYNAMICAL ) THEN
               DEALLOCATE (DVDR,DBDR,MECV,MEVC)
               DEALLOCATE (MECVX,MEVXC,ABSRTA,ABSRTB)
               DEALLOCATE (ABSRATE,TAUCORE)
C
               CALL XRAYSPECDYN(CALCINT,GETIRRSOL,ETAB0,TAUQ,TAUT,MSSQ,
     &                          TSST,MSST,SSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,
     &                          SZCOR,KAPCOR,MM05COR,NKPCOR,IKMCOR,
     &                          IZERO,ITXRAY,BCOR,BCORS,NCSTMAX,NPOLMAX)
C
            END IF
         END IF
C
         IF ( TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
            OUTPUT = 'SIGMA     '
            CALL SECTION_SET_INTEGER('NE3',NE3,9999,0)
            CALL SECTION_SET_REAL('EMAX3',EMAX3,9999D0,0)
            CALL SECTION_SET_INTEGER('TSELECT',TSELECT,9999,0)
            CALL SECTION_SET_REAL_ARRAY('TAUCORE',TAUCOREK,N_FOUND,2,0,
     &                                  9999D0,0)
         END IF
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'task xas  not found in input')
C
      END IF
C
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
      NE = 0
C
      IF ( THERMAL_VIBRA_FLUCT ) CALL TAU_READ(NE,ERYD,MSSQ,TAUQ)
C
      IF ( NE.EQ.0 ) THEN
C
         TAU_FILE_FMT = 'OLD'
C
         CALL READTAU(9,ERYD,0,NE,TAUT,0,NT,.TRUE.,TAUQ,MSSQ,0,NQ,0,
     &                NKMMAX,NKMMAX,IPRINT)
C
      END IF
C
      C = 2.0D0*1.370373D+02
      NREG = -1D0/(0.5D0*C)
      NIRR = NREG*NREG
C
      VUC = 0D0
      DO IQ = 1,NQ
         JT = ITOQ(1,IQ)
         IM = IMT(JT)
         VUC = VUC + (4D0*PI/3D0)*RWS(IM)**3
      END DO
C
      XSCL = 1.0D0
C
      SCLSIG = (E0_CGS**2*HBAR_CGS/(0.5D0*M0_CGS*A0_CGS**3*EV_ERG*RY_EV)
     &         )*(0.5D0*C)**2*(1D0/(2D0*PI))*XSCL*1D-15
C
      SCLSIGALT = (E0_CGS**2*HBAR_CGS*C_CGS**2/A0_CGS**3)
     &            *(1D0/(EV_ERG*RY_EV)**2)*(1D0/(2D0*PI))*XSCL*1D-15
C
      SCLSIG1 = (PI*E0_CGS**2*HBAR_CGS/(0.5D0*M0_CGS*A0_CGS**3*EV_ERG*
     &          RY_EV))*(0.5D0*C)**2*(1.0D0/PI)*XSCL*1D-15
C
      SCLRATE = (PI*E0_CGS**2*HBAR_CGS/(0.5D0*M0_CGS*EV_ERG*RY_EV))
     &          *(4*PI/C_CGS)*(0.5D0*C)**2*(1.0D0/PI)*XSCL*1D+18
C
      SCLXRS = (M0_CGS*C_CGS**2/E0_CGS**2)*(1D0/(4*PI*HBAR_CGS*C_CGS))
     &         *(EV_ERG*RY_EV)*SCLRATE*1D-18
C
      WRITE (6,99001) SCLSIG,SCLSIGALT,SCLSIG1,SCLRATE,SCLXRS,VUC
C
      SCLSIG = SCLSIG/VUC
      SCLSIG1 = SCLSIG1/VUC
      IF ( TASK.EQ.'XRS       ' ) SCLSIG = SCLXRS
C
      KUNIT = 2
      IF ( TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
         IF ( TASK.EQ.'XMO       ' ) WRITE (6,99024)
         IF ( TASK.EQ.'XRS       ' ) WRITE (6,99025)
         IF ( TSELECT.EQ.1 ) THEN
            WTERM2 = 0D0
            WRITE (6,99026) 'only the resonant term'
         ELSE IF ( TSELECT.EQ.2 ) THEN
            WTERM1 = 0D0
            WRITE (6,99026) 'only the non-resonant term'
         ELSE
            WRITE (6,99026) 'both terms'
         END IF
      ELSE IF ( OUTPUT(1:5).EQ.'SIGMA' ) THEN
         WRITE (6,99023) 'units of  [10^15 1/s]'
         KUNIT = 2
         SCLRATE = SCLSIG1
      ELSE
         WRITE (6,99023) 'units of  [Mb]'
         KUNIT = 1
      END IF
C
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C                               xray
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
      OPEN (UNIT=IFILCORWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
      NTXRAY = 1
      IDIPOL = 0
C
      CALL SECTION_SET_INTEGER('DIPOL',IDIPOL,9999,0)
      IF ( IDIPOL.LT.0 .OR. IDIPOL.GT.1 ) IDIPOL = 0
C
C ------------------------------------ calculate angular matrix elements
C
      CALL AME_INIT(MEFORM,6)
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO JT = 1,NTXRAY
C
         CALL SECTION_SET_INTEGER('IT',IT,1,0)
C
         CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
         CALL SECTION_SET_REAL_ARRAY('TAUCORE',TAUCOREK,N_FOUND,2,0,
     &                               9999D0,0)
C
         IM = IMT(IT)
         C = CTL(IT,1)
C
         IF ( LDATSET.NE.0 ) THEN
            FILNAM = DATSET(1:LDATSET)//'_'//TXT_T(IT)(1:LTXT_T(IT))
     &               //'_'
            LFN = LDATSET + LTXT_T(IT) + 2
         ELSE
            FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'_'
            LFN = LTXT_T(IT) + 1
         END IF
C
         FILNAM = FILNAM(1:LFN)//SHELL(NCXRAY(IT))
         LFN = LFN + 1
         SPEC = '  '//SHELL(NCXRAY(IT))
C
         IF ( NCXRAY(IT).NE.1 ) THEN
            SPEC = SPEC(1:3)//SUBSH(LCXRAY(IT))
            LSP = 3 + LSUBSH(LCXRAY(IT))
            FILNAM = FILNAM(1:LFN)//SUBSHP(LCXRAY(IT))
            LFN = LFN + MAX(1,LSUBSH(LCXRAY(IT))-1)
         ELSE
            LSP = 3
         END IF
C
         FILNAM = FILNAM(1:LFN)//'.rat'
         SPEC = SPEC(1:LSP)//' - XAS spectrum of '//TXT_T(IT)
     &          (1:LTXT_T(IT))//' in  '//SYSTEM(1:LSYSTEM)
C
         OPEN (UNIT=7,FILE=FILNAM(1:LFN+4))
         WRITE (6,'(10X,A,A,/)') 'RATE-FILE :  ( 7) ',FILNAM(1:LFN+4)
C
         ITXRAY = IT
C
         WRITE (6,'(3X,A,/)') SPEC
         WRITE (6,99019) IDIPOL,MEFORM
         IF ( TASK.NE.'XMO       ' .AND. TASK.NE.'XRS       ' ) THEN
            WRITE (6,99020) 'NPOL',NPOL,(CHPOL1(IPOL),IPOL=1,NPOL)
         ELSE
            WRITE (6,99020) 'NSIG',9,((CHSIG(I,J),I=1,3),J=1,3)
         END IF
C
         WRITE (6,99021) IT,NCXRAY(IT),LCXRAY(IT)
C
         CALL WRHEAD(7,FILNAM,TASK,NE)
C
         NTXRSGRP = 1
         WRITE (7,99016) 'NTXRSGRP  ',NTXRSGRP
         WRITE (7,99016) 'IDIPOL    ',IDIPOL
         IF ( TASK.NE.'XMO       ' ) THEN
            WRITE (7,99017) 'NPOL      ',NPOL,(CHPOL1(IPOL),IPOL=1,NPOL)
         ELSE
            WRITE (7,99017) 'NSIG      ',9,((CHSIG(I,J),I=1,3),J=1,3)
         END IF
         WRITE (7,99018) 'SPECTRUM  ',SPEC(1:LSP)
         WRITE (7,99016) 'IT        ',IT
         WRITE (7,99016) 'NCXRAY    ',NCXRAY(IT)
         WRITE (7,99016) 'LCXRAY    ',LCXRAY(IT)
C
         NCST = 4*LCXRAY(IT) + 2
C
         IF ( NCST.GT.NCSTMAX ) THEN
            WRITE (6,99002) LCXRAY(IT),NCSTMAX
            CALL STOP_MESSAGE(ROUTINE,' LCXRAY, NCSTMAX ???')
         END IF
C
C=========================================================== CORE STATES
C               calculate core states and store in file
C=========================================================== CORE STATES
C
         CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &             IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
         DO IFIL = 6,7
            WRITE (IFIL,99005) NCST,(NKPCOR(ICST),ICST=1,NCST)
            WRITE (IFIL,99006)
         END DO
C
         ECOR_KM(:) = 0D0
C
         DO ICST = 1,NCST
C
            DO K = 1,NKPCOR(ICST)
               XNORM(K) = 0.0D0
               DO N = 1,JRWS(IM)
                  XNORM(K) = XNORM(K) + R2DRDI_W_RADINT(N,IM)
     &                       *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
               END DO
            END DO
C
            DO IFIL = 6,7
               WRITE (IFIL,99007) ICST,NCXRAY(IT),LCXRAY(IT),
     &                            KAPCOR(ICST),(2*MM05COR(ICST)+1),
     &                            IKMCOR(ICST,1),XNORM(1),ECOR(ICST),
     &                            ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                            IZERO(ICST)
               IF ( NKPCOR(ICST).EQ.2 ) WRITE (IFIL,99009)
     &              IKMCOR(ICST,2),XNORM(2)
            END DO
C
            IF ( KAPCOR(ICST).EQ.LCXRAY(IT) ) THEN
               KTYP = 1
            ELSE
               KTYP = 2
            END IF
            IF ( TAUCOREK(KTYP).LT.1D-6 ) THEN
               TAUCORE(ICST) = DBLE(WCOREHOLE(Z(IT),NCXRAY(IT),LCXRAY(IT
     &                         ),KTYP))/RY_EV
            ELSE
               TAUCORE(ICST) = TAUCOREK(KTYP)/RY_EV
            END IF
            WRITE (6,99008) TAUCORE(ICST),TAUCORE(ICST)*RY_EV
C
            IKM_CST = IKMCOR(ICST,1)
            ECOR_KM(IKM_CST) = ECOR(ICST)
C
         END DO
C
         WRITE (7,99010) NE,2,KUNIT,'VUC       ',VUC
C
C ======================================================================
C                          MEFORM  =  GRV
C ======================================================================
C
         IF ( MEFORM.EQ.'GRV' ) THEN
            CALL DVDRSPLINE(VT(1,IT),R(1,IM),DVDR,JRWS(IM))
C
            CALL DVDRSPLINE(BT(1,IT),R(1,IM),DBDR,JRWS(IM))
         END IF
C
C ======================================================================
C ============================================ conduction band = START =
C
C   ********************************************************************
C        the following section:
C             -reads in the energy and tau-matrix,
C             -calls the routine that calc. the radial matrix elements
C             -sums the various contributions to the total rate
C             -prints rates on file 7
C   ********************************************************************
         WRITE (6,99003)
C
         REWIND 9
C
         READ (9,'(////,10X,I3,5X,I3,6X,I2)') NT9,NQ9
         IF ( NT.NE.NT9 .OR. NQ.NE.NQ9 )
     &         CALL STOP_MESSAGE(ROUTINE,'NT9,NQ9 ???')
C
         READ (9,'(A)') (S10DUMMY,I=1,(NT9+NQ))
C
C ================================================================
C      calculation of FULL SIGMA TENSOR or SCATTERING AMPLITUDE
C ======================================================================
C
         IF ( TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
C
            IF ( NTXRAY.GT.1 ) WRITE (6,*) 'WARNING NTXRAY set to 1'
            NTXRAY = 1
C
            NESIG = 50
            DELE = 4.0D0
C
            ALLOCATE (ESIG(NEMAX+NESIG))
            ALLOCATE (SIG(NEMAX+NESIG,NCSTMAX,3,3),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SIG')
C
            CALL CINIT(NEMAX*NCSTMAX*3*3,SIG)
C
            DO IE = 1,NESIG
               ESIG(IE) = DREAL(ETAB(1,1))
     &                    - DELE*(DBLE(NESIG-IE-1)/DBLE(NESIG))**3
            END DO
            DO IE = 1,NE
               NESIG = NESIG + 1
               ESIG(NESIG) = DREAL(ETAB(IE,1))
            END DO
C
            NELOOP = 3
C
         ELSE
C
            NELOOP = 1
C
         END IF
C ================================================================
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
C ================================================================
C =================================== ENERGY LOOP === START ======
         DO IELOOP = 1,NELOOP
C
C
            IF ( IELOOP.EQ.2 ) THEN
               DO IE = 1,NE
                  ETAB(IE,1) = DREAL(ETAB(IE,1))
               END DO
            END IF
C
            IF ( IELOOP.EQ.3 ) THEN
               EIMAG = 0D0
               EILOW = 0D0
               CALL EPATH(7,EMAX,EMAX3,EIMAG,NE3,ETAB(1,1),WETAB(1,1),
     &                    EILOW,IPRINT,NEMAX)
               NE = NE3
            END IF
C
            DO IENERG = 1,NE
C
               TAUT(:,:,IT) = 0D0
C
               IF ( IELOOP.NE.1 ) THEN
C
                  ERYD = ETAB(IENERG,1)
C
                  TAUT(:,:,IT) = TSST(:,:,IT)
C
               ELSE IF ( TAU_FILE_FMT.EQ.'NEW' ) THEN
C
                  CALL TAU_READ(NE,ERYD,MSSQ,TAUQ)
C
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
C
               ELSE
C
                  CALL READTAU(9,ERYD,IENERG,NE,TAUT(1,1,IT),IT,NT,
     &                         .TRUE.,TAUQ,MSSQ,0,NQ,1,NKM,NKMMAX,
     &                         IPRINT)
C
               END IF
C
C
               IF ( SELFENERGY(1:4).EQ.'NONE' .OR. FULLPOT ) THEN
C
                  IF ( DMFT ) THEN
C
                     IF ( IELOOP.EQ.1 ) THEN
                        DMFTSIG(1:NKM,1:NKM,1:NT)
     &                     = DMFTSIGMA(1:NKM,1:NKM,1:NT,IENERG)
                     ELSE IF ( IELOOP.EQ.2 ) THEN
                        DMFTSIG(1:NKM,1:NKM,1:NT)
     &                     = DREAL(DMFTSIGMA(1:NKM,1:NKM,1:NT,IENERG))
                     ELSE
                        DMFTSIG(1:NKM,1:NKM,1:NT) = C0
                     END IF
                     DO I = 1,NT
                        IF ( KSELF(I).EQ.1 ) WRITE (6,99027) I
                     END DO
                  END IF
C
                  CALL RUNSSITE(CALCINT,1,1,IFILCBWF,GETIRRSOL,ERYD,P,
     &                          IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
                  ETAB0(IENERG) = ERYD
C
                  IF ( ABS(ETAB(IENERG,1)-ERYD).GT.1D-6 ) THEN
                     WRITE (6,*) ' TROUBLE with energies for IE',IENERG
                     WRITE (6,*) ' ETAB        ',ETAB(IENERG,1)
                     WRITE (6,*) ' ERYD (in)   ',ERYD
                  END IF
C
               ELSE
C
                  ERYDCHK = ERYD
C
                  IF ( ABS(ERYDCHK-ERYD).GT.1D-6 ) THEN
                     WRITE (6,*) ' TROUBLE with <CVSSITE> for IE',IENERG
                     WRITE (6,*) ' ETAB0       ',ETAB0(IENERG)
                     WRITE (6,*) ' ERYD (in)   ',ERYDCHK
                     WRITE (6,*) ' ERYD (new)  ',ERYD
                  END IF
C
                  CALL STOP_MESSAGE(ROUTINE,'<CVSSITE>')
C
               END IF
C
C ======================================================================
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
               CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
C --------------------------------------------- PRINT ------ START -----
               IF ( IPRINT.GT.0 ) WRITE (6,*) '  ENERGY : ',ERYD
C
               STRWARN = '       '
               IF ( TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' )
     &              THEN
                  DO ICST = 1,NCST
                     IF ( 2D0*DIMAG(ERYD).GE.TAUCORE(ICST) ) THEN
                        TAUCOREWARNING = .TRUE.
                        STRWARN = 'Im(E) ?'
                     END IF
                  END DO
               END IF
C
               WRITE (6,99011) IELOOP,IENERG,ERYD,IQAT(1,IT),IT,STRWARN
C
               IF ( IPRINT.GE.5 ) THEN
C
                  CALL CMATSTRUCT('TAUT(KAP,MUE)',TAUT(1,1,IT),NKM,
     &                            NKMMAX,3,3,0,1D-8,6)
C
                  WRITE (6,*) ' '
               END IF
C --------------------------------------------- PRINT ------- END ------
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C NOTE: the routine uses TMAT = TAU_t or GMAT_t in the LOCAL frame
C       TMAT will be rotated to the LOCAL frame below if necessary
C
               CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT_LOC)
C
C=======================================================================
C=======================================================================
C     calculate matrix elements in the LOCAL frame for ALL polarisations
C     using SPHERICAL coordinates     ipol = 1,2,3  ==  (-),(z),(+)
C=======================================================================
C=======================================================================
C
C---------------------------------- suppress surface terms in ME_ALF_ALF
               CECOR = 0D0
C
               CALL ME_ALF_ALF(1,NKM,IFILCORWF,CECOR,1,NKM,IFILCBWF,
     &                         ERYD,.FALSE.,IT,MZAZB,MZBZA,MIRR_2,
     &                         MIRR_3,MIRR_4,CTL(IT,1),K_CALC_ME)
C
C use the old scheme to calculate the non-dipolar (B) matrix elements MEB
C   NOTE: <XRAYME> gives only  MEBIRR(+,-), MEBIRR(-,+), MEBIRR(z,z)
C
               IF ( .NOT.FULLPOT .AND. IDIPOL.EQ.0 ) THEN
C
                  CALL ME_ALF_ALF_NONDIP(1,NKM,IFILCORWF,ECOR_KM,1,NKM,
     &               IFILCBWF,ERYD,.FALSE.,IT,MZAZB_NONDIPX,
     &               MZBZA_NONDIPX,MIRR_2_NONDIP,MIRR_3_NONDIPX,
     &               MIRR_4_NONDIP,CTL(IT,1),K_CALC_ME)
C
                  MIRR_3_NONDIPX(:,:,:,:,:) = -MIRR_3_NONDIPX(:,:,:,:,:)
                  DO I = 1,NKM
                     DO J = 1,NKM
                        IF ( I.NE.J ) MIRR_3_NONDIPX(I,J,:,:,2) = 0
                     END DO
                  END DO
C
                  MZAZB_NONDIP(:,:,:,1) = MZAZB_NONDIPX(:,:,:,2)
                  MZBZA_NONDIP(:,:,:,1) = MZBZA_NONDIPX(:,:,:,2)
                  MIRR_3_NONDIP(:,:,:,:,:) = 0D0
                  MIRR_3_NONDIP(:,:,1,3,1) = MIRR_3_NONDIPX(:,:,1,3,2)
                  MIRR_3_NONDIP(:,:,2,2,1) = MIRR_3_NONDIPX(:,:,2,2,2)
                  MIRR_3_NONDIP(:,:,3,1,1) = MIRR_3_NONDIPX(:,:,3,1,2)
C
               END IF
C
C
C----------------------------------------------------------------------
C           transformation from spherical to CARTESIAN coordinates
C              from now on        ipol = 1,2,3  ==  (x),(y),(z)
C----------------------------------------------------------------------
C
               CALL CMAT_CONVERT_POLAR(MZAZB,'S>C')
               CALL CMAT_CONVERT_POLAR(MZBZA,'S>C')
C
               CALL CMAT_CONVERT_POLAR2(MIRR_3,'S>C')
C
               IQ = IQAT(1,IT)
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C   deal with rotated magnetisation and/or changed frame of reference
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
               IF ( CHANGE_FRAME .OR. MOMENTS_ROTATED ) THEN
C
C--- provide site diagonal multiple scattering matrix in the LOCAL frame
C
                  IF ( MOMENTS_ROTATED ) THEN
C
                     WKM1(:,:) = TMAT_LOC(:,:)
C
                     CALL ROTATE(WKM1,'G->L',TMAT_LOC,NKM,DROTQ(1,1,IQ),
     &                           NKMMAX)
C
                  END IF
C
                  WRITE (6,*) '####### CHANGE_FRAME    ',CHANGE_FRAME
                  WRITE (6,*) '####### MOMENTS_ROTATED ',MOMENTS_ROTATED
C
C----------------------------------------------------------------------
C  transformation of polarisation from LOCAL to GLOBAL frame
C  in case of rotated magnetic moment on site IQ
C  NOTE:  rotation resticted to polarisation by use of  C1MAT = 1
C----------------------------------------------------------------------
C
                  IF ( MOMENTS_ROTATED ) THEN
C
                     MZAZB_TMP(:,:,:) = MZAZB(:,:,:)
                     MZBZA_TMP(:,:,:) = MZBZA(:,:,:)
C
                     MIRR_3_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
C
                     MIRR_2_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
                     MIRR_4_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
C
                     CALL ME_ROTATE_REG(1,1,C1MAT,MROTQ(1,1,IQ),
     &                                  MZAZB_TMP,MZBZA_TMP,MZAZB,MZBZA)
C
                     CALL ME_ROTATE_IRR(1,1,C1MAT,MROTQ(1,1,IQ),
     &                                  MIRR_2_TMP,MIRR_3_TMP,
     &                                  MIRR_4_TMP,MIRR_2,MIRR_3,MIRR_4)
C
                  END IF
C
                  MZAZB_NONDIP(:,:,:,:) = C0
                  MZBZA_NONDIP(:,:,:,:) = C0
                  MIRR_3_NONDIP(:,:,:,:,:) = C0
C
               END IF
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
               ABSRTA(:,:) = C0
               ABSRTB(:,:) = C0
               ABSRATE(:,:) = C0
C
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
               IFT = (IT-1)*NFLUCT
               LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
                  IF ( NFLUCT.EQ.1 ) THEN
C
                     TSS_FT(:,:) = TSST(:,:,IT)
C
                     MZAZB_FT(:,:,:) = MZAZB(:,:,:)
                     MZBZA_FT(:,:,:) = MZBZA(:,:,:)
                     MIRR_2_FT(:,:,:,:) = MIRR_2(:,:,:,:)
                     MIRR_3_FT(:,:,:,:) = MIRR_3(:,:,:,:)
                     MIRR_4_FT(:,:,:,:) = MIRR_4(:,:,:,:)
C
                  ELSE
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSST(1,1,IT),'SAS+',TSS_FT)
C
                     CALL ME_ROTATE_REG(IFT,NFT,FMAT_FT,MFMAT_FT,MZAZB,
     &                                  MZBZA,MZAZB_FT,MZBZA_FT)
C
                     CALL ME_ROTATE_IRR(IFT,NFT,FMAT_FT,MFMAT_FT,MIRR_2,
     &                                  MIRR_3,MIRR_4,MIRR_2_FT,
     &                                  MIRR_3_FT,MIRR_4_FT)
C
                  END IF
C
C============================================================= IVIBRA ==
C                                             perform local displacement
                  IVT = (IT-1)*NVIBRA
                  LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                     IVT = IVT + 1
C
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
C
                     IF ( NVIBRA.EQ.1 ) THEN
C
                        MSS_VQ(:,:) = MSSQ(:,:,IQ)
                        TAU_VQ(:,:) = TAUQ(:,:,IQ)
C
                     ELSE
C
C------------------------------- TAU_VQ: TAU_vq = U_v^(-1) * TAU_q * U_v
C
                        CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TAUQ(1,1,IQ),
     &                     'UTAU',TAU_VQ)
C
                        CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),MSSQ(1,1,IQ),
     &                     'UTAU',MSS_VQ)
C
                     END IF
C
C----------------------------------------------- MSS_FT: m_ft = 1 / t_ft
C
                     IF ( NFLUCT.EQ.1 ) THEN
                        MSS_FT(:,:) = MSST(:,:,IT)
                     ELSE
                        CALL CMATINV3(N,M,IPIVKM,TSS_FT,WKM1,MSS_FT)
                     END IF
C
C-------------------------------------- get projected TAU matrix TAU_vft
C-------------------------------------------------- TAU(t) = TAU * D~(t)
C
                     IF ( NVIBFLU.EQ.1 ) THEN
C
                        TAU_VFT(:,:) = TMAT_LOC(:,:)
C
                     ELSE
C
                        CALL GET_TAUT(MSS_FT,MSS_VQ,TAU_VQ,TAU_VFT)
C
                     END IF
C
C ======================================================================
C                            X-ray absorption
C ======================================================================
                     IF ( TASK.NE.'XMO       ' .AND. 
     &                    TASK.NE.'XRS       ' ) THEN
C
C----------------------------------------------------------------------
C  transformation of polarisation from GLOBAL to SPECIAL frame
C  aligned to photon beam
C  NOTE:  rotation resticted to polarisation by use of  C1MAT = 1
C----------------------------------------------------------------------
C
                        IF ( CHANGE_FRAME ) THEN
C
                           MROT = MROT_FRAME
C
                           MROT = TRANSPOSE(MROTQ(1:3,1:3,IQ))
                           MROT = 0D0
                           MROT(1,1) = 1
                           MROT(2,2) = 1
                           MROT(3,3) = 1
C
                           MROT = TRANSPOSE(MROT_FRAME)
C
                           WRITE (*,*) ' MROT_FRAME '
                           WRITE (*,'(3f10.3)') MROT_FRAME
                           WRITE (*,*) ' MROTQ '
                           WRITE (*,'(3f10.3)') MROTQ(1:3,1:3,IQ)
C
                           MZAZB_TMP(:,:,:) = MZAZB_FT(:,:,:)
                           MZBZA_TMP(:,:,:) = MZBZA_FT(:,:,:)
C
                           MIRR_3_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
C
                           MIRR_2_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
                           MIRR_4_TMP(:,:,:,:) = MIRR_3(:,:,:,:)
C
                           CALL ME_ROTATE_REG(1,1,C1MAT,MROT,MZAZB_TMP,
     &                        MZBZA_TMP,MZAZB_FT,MZBZA_FT)
C
                           CALL ME_ROTATE_IRR(1,1,C1MAT,MROT,MIRR_2_TMP,
     &                        MIRR_3_TMP,MIRR_4_TMP,MIRR_2,MIRR_3,
     &                        MIRR_4)
C
                        END IF
C
C----------------------------------------------------------------------
C           transformation from cartesian to SPHERICAL coordinates
C              from now on        ipol = 1,2,3  ==  (-),(z),(+)
C----------------------------------------------------------------------
C
                        CALL CMAT_CONVERT_POLAR(MZAZB_FT,'C>S')
                        CALL CMAT_CONVERT_POLAR(MZBZA_FT,'C>S')
C
                        CALL CMAT_CONVERT_POLAR2(MIRR_3,'C>S')
C
                        MZAZB_FT(:,:,:) = MZAZB_FT(:,:,:)*NREG
                        MZBZA_FT(:,:,:) = MZBZA_FT(:,:,:)*NREG
                        MIRR_3(:,:,:,:) = MIRR_3(:,:,:,:)*NIRR
C
C ######################################################################
C              calculate the resulting absorption rate
C     using SPHERICAL coordinates     ipol = 1,2,3  ==  (-),(z),(+)
C ######################################################################
C
C    - IM [ SUM(LAM,LAM')  M(I,LAM) # TAU(LAM,LAM') # M(LAM',I)# ]
C    + IM [ SUM(LAM)       I(I,LAM) ]
C
C ######################################################################
C
                        DO IPOL = 1,NPOL
C
                           IPOLXC = IPOLREV(IPOL)
C
                           PHASE_FACTOR = PHAS_POL(IPOL)
     &                        *PHAS_POL(IPOLXC)
C
                           DO ICST = 1,NCST
C
                              IKM_CST = IKMCOR(ICST,1)
C
                              CSCL = SCLRATE/(ETAB0(IENERG)-ECOR(ICST))
C
                              RATA = PHASE_FACTOR*MIRR_3(IKM_CST,
     &                               IKM_CST,IPOLXC,IPOL)
C
                              RATB = PHASE_FACTOR*MIRR_3_NONDIP(IKM_CST,
     &                               IKM_CST,IPOLXC,IPOL,1)
C
                              DO LAM1 = 1,NKM
                                 DO LAM2 = 1,NKM
C
                                    TAU12 = TAU_VFT(LAM1,LAM2)
C
                                    RATA = RATA - 
     &                                 PHASE_FACTOR*MZAZB_FT(IKM_CST,
     &                                 LAM1,IPOL)
     &                                 *TAU12*MZBZA_FT(LAM2,IKM_CST,
     &                                 IPOLXC)
C
                                    RATB = RATB - 
     &                                 PHASE_FACTOR*MZAZB_NONDIP
     &                                 (IKM_CST,LAM1,IPOL,1)
     &                                 *TAU12*MZBZA_NONDIP(LAM2,IKM_CST,
     &                                 IPOLXC,1)
C
                                 END DO
                              END DO
C
                              ABSRTA(ICST,IPOL) = ABSRTA(ICST,IPOL)
     &                           + DIMAG(RATA)
                              ABSRTB(ICST,IPOL) = ABSRTB(ICST,IPOL)
     &                           + DIMAG(RATB)
                              ABSRATE(ICST,IPOL) = ABSRATE(ICST,IPOL)
     &                           + DIMAG(RATA+RATB)
C
                           END DO
                        END DO
C
                        DO ICST = 1,NCST
                           CSCL = SCLRATE/(ETAB0(IENERG)-ECOR(ICST))
C----------------------------------------- avoid imaginary part for RATE
C                       CSCL = DREAL(CSCL)
                           ABSRTA(ICST,:) = ABSRTA(ICST,:)*CSCL
                           ABSRTB(ICST,:) = ABSRTB(ICST,:)*CSCL
                           ABSRATE(ICST,:) = ABSRATE(ICST,:)*CSCL
                        END DO
C
C---- IPOL_NEW_OLD:        print out spectra in OLD sequence (+),(-),(z)
C
                        WRITE (7,99013) ETAB0(IENERG),
     &                                  DREAL(ETAB0(IENERG))*RY_EV
                        DO ICST = 1,NCST
                           DO IPOL = 1,NPOL
                              WRITE (7,99014) ICST,IPOL,
     &                               ABSRTA(ICST,IPOL_NEW_OLD(IPOL)),
     &                               ABSRTB(ICST,IPOL_NEW_OLD(IPOL)),
     &                               ABSRATE(ICST,IPOL_NEW_OLD(IPOL)),
     &                               CHPOL1(IPOL)
                           END DO
                        END DO
                        WRITE (7,*) ' '
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                        IF ( WRBUILDBOT .AND. IENERG.LE.3 )
     &                       WRITE (IFILBUILDBOT,99028)
     &                       ROUTINE(1:LEN_TRIM(ROUTINE)),IT,IENERG,
     &                       ETAB0(IENERG),
     &                       ((DREAL(ABSRTA(ICST,IPOL_NEW_OLD(IPOL))),
     &                       DREAL(ABSRTB(ICST,IPOL_NEW_OLD(IPOL))),
     &                       IPOL=1,NPOL),ICST=1,NCST)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
                     ELSE
C
C ======================================================================
C      calculation of FULL SIGMA TENSOR or SCATTERING AMPLITUDE
C ======================================================================
C
C
C ======================================================================
C                            matrix elements
C ======================================================================
C                index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C                 M(X) =   [  M(+) + M(-) ] / SQRT(2)
C                 M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
C   MEA is  PURE IMAGINARY  for REAL  energies   !!
C
                        DO J = 1,NKM
                           DO ICST = 1,NCST
C
                              IKM_CST = IKMCOR(ICST,1)
C
                              DO IPOL = 1,3
                                 MCV(IPOL) = MZAZB_FT(IKM_CST,J,IPOL)
                                 MVC(IPOL)
     &                              = -MZAZB_FT(IKM_CST,J,IPOLREV(IPOL))
                                 MCVX(IPOL) = -DCMPLX(MCV(IPOL))
                                 MVXC(IPOL) = -DCMPLX(MVC(IPOL))
                              END DO
C
                              MECV(J,ICST,1) = (MCV(1)+MCV(2))/SQRT_2
                              MECV(J,ICST,2) = CI*(-MCV(1)+MCV(2))
     &                           /SQRT_2
                              MECV(J,ICST,3) = MCV(3)
C
                              MEVC(J,ICST,1) = (MVC(1)+MVC(2))/SQRT_2
                              MEVC(J,ICST,2) = CI*(-MVC(1)+MVC(2))
     &                           /SQRT_2
                              MEVC(J,ICST,3) = MVC(3)
C
                              MECVX(J,ICST,1) = (MCVX(1)+MCVX(2))/SQRT_2
                              MECVX(J,ICST,2) = CI*(-MCVX(1)+MCVX(2))
     &                           /SQRT_2
                              MECVX(J,ICST,3) = MCVX(3)
C
                              MEVXC(J,ICST,1) = (MVXC(1)+MVXC(2))/SQRT_2
                              MEVXC(J,ICST,2) = CI*(-MVXC(1)+MVXC(2))
     &                           /SQRT_2
                              MEVXC(J,ICST,3) = MVXC(3)
C
                           END DO
                        END DO
C from now on >>>>>>>>>>>>>>>>>   index 3:  ipol= 1,2,3  ==  (x),(y),(z)
C
C
                        DO ICST = 1,NCST
C
                           DO IX1 = 1,3
                              DO IX2 = 1,3
C
                                 CSUM1 = C0
                                 CSUM2 = C0
C
                                 DO LAM1 = 1,NKM
                                    DO LAM2 = 1,NKM
C
                                       IF ( IELOOP.EQ.1 ) THEN
                                         TPLS = TMAT_LOC(LAM1,LAM2)
     &                                      - TSST(LAM1,LAM2,IT)
                                         TMIN = TMAT_LOC(LAM2,LAM1)
     &                                      - TSST(LAM2,LAM1,IT)
                                       ELSE
                                         TPLS = TSST(LAM1,LAM2,IT)
                                         TMIN = TSST(LAM2,LAM1,IT)
C
                                         IF ( ABS(DREAL(MEVC(LAM2,ICST,1
     &                                      ))).GT.1D-6 ) WRITE (6,*)
     &                                       'WARNING Re(M)<>0 ',LAM2,
     &                                      ICST,' X '
                                         IF ( ABS(DIMAG(MEVC(LAM2,ICST,2
     &                                      ))).GT.1D-6 ) WRITE (6,*)
     &                                       'WARNING Re(M)<>0 ',LAM2,
     &                                      ICST,' Y '
                                         IF ( ABS(DREAL(MEVC(LAM2,ICST,3
     &                                      ))).GT.1D-6 ) WRITE (6,*)
     &                                       'WARNING Re(M)<>0 ',LAM2,
     &                                      ICST,' Z '
                                       END IF
                                       TMINX = DCONJG(TMIN)
C
                                       CSUM1 = CSUM1 + 
     &                                    MECV(LAM1,ICST,IX1)
     &                                    *TPLS*MEVC(LAM2,ICST,IX2)
     &                                    - MECVX(LAM1,ICST,IX1)
     &                                    *TMINX*MEVXC(LAM2,ICST,IX2)
C
                                       CSUM2 = CSUM2 + 
     &                                    MECV(LAM1,ICST,IX2)
     &                                    *TPLS*MEVC(LAM2,ICST,IX1)
     &                                    - MECVX(LAM1,ICST,IX2)
     &                                    *TMINX*MEVXC(LAM2,ICST,IX1)
C
                                    END DO
                                 END DO
C
                                 CSUM1 = WTERM1*SCLSIG*CSUM1
                                 CSUM2 = WTERM2*SCLSIG*CSUM2
C
                                 DO IESIG = 1,NESIG
C
                                    EPHOT = ESIG(IESIG) - ECOR(ICST)
C
                                    IF ( TASK.EQ.'XRS       ' ) THEN
                                       WE = WETAB(IENERG,1)*EPHOT
                                    ELSE
                                       WE = WETAB(IENERG,1)
                                    END IF
C
                                    ETERM1 = (ETAB(IENERG,1)-ECOR(ICST)
     &                                 -CI*DELTA)
     &                                 *(EPHOT+ECOR(ICST)-ETAB(IENERG,1)
     &                                 +CI*TAUCORE(ICST))
C
                                    ETERM2 = (ETAB(IENERG,1)-ECOR(ICST)
     &                                 +CI*DELTA)
     &                                 *(EPHOT-ECOR(ICST)+ETAB(IENERG,1)
     &                                 +CI*TAUCORE(ICST))
C
                                    SIG(IESIG,ICST,IX1,IX2)
     &                                 = SIG(IESIG,ICST,IX1,IX2)
     &                                 - WE*CSUM1/ETERM1 - 
     &                                 WE*CSUM2/ETERM2
C
                                 END DO
C
                              END DO
C
                           END DO
C
C
                           DO IESIG = 1,NESIG
                              IF ( ABS(SIG(IESIG,ICST,1,1)-SIG(IESIG,
     &                             ICST,2,2)).GT.SIGTOL ) WRITE (6,*)
     &                              'WARNING SIG(xx) <> SIG(yy) '
                              IF ( ABS(SIG(IESIG,ICST,1,2)+SIG(IESIG,
     &                             ICST,2,1)).GT.SIGTOL ) WRITE (6,*)
     &                              'WARNING SIG(xy) <> -SIG(yx) '
                           END DO
C
                        END DO
C
C ======================================================================
                     END IF
C
C
C
C
                  END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
               END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
C
            END DO
C
         END DO
C ========================================= energy loop ==== END =======
C ======================================================================
C
C ================================================================
         IF ( TASK.EQ.'XMO       ' ) THEN
C
            DO IESIG = 1,NESIG
C
               WRITE (7,99013) ESIG(IESIG),0D0,ESIG(IESIG)*RY_EV
               DO ICST = 1,NCST
                  ISIG = 0
                  DO IX2 = 1,3
                     DO IX1 = 1,3
                        ISIG = ISIG + 1
C
                        WRITE (7,99015) ICST,ISIG,
     &                                  SIG(IESIG,ICST,IX1,IX2)
                     END DO
                  END DO
               END DO
               WRITE (7,*) ' '
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT .AND. IESIG.LE.3 )
     &              WRITE (IFILBUILDBOT,99029)
     &              ROUTINE(1:LEN_TRIM(ROUTINE)),IT,IESIG,ESIG(IESIG),
     &              (((DREAL(SIG(IESIG,ICST,IX1,IX2)),
     &              DIMAG(SIG(IESIG,ICST,IX1,IX2)),ICST=1,NCST),IX1=1,3)
     &              ,IX2=1,3)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            END DO
         END IF
C
         WRITE (6,99012) NE
         IF ( TAUCOREWARNING ) WRITE (6,99022)
C   ********************************************************************
C
C
C ============================================ conduction band == END ==
C ======================================================================
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      STOP
C
99001 FORMAT (5X,'scaling factors:',/,5X,'SCLSIG                 ',
     &        E15.6,'  >> SIGMA  in [10^15 1/s]',/,5X,
     &        'SCLSIG  (alternative)  ',E15.6,/,5X,
     &        'SCLSIG1 (absorption)   ',E15.6,
     &        '  >> SIGMA1 in [10^15 1/s]',/,5X,
     &        'SCLRATE (abs. coeff.)  ',E15.6,
     &        '  >> MUE    in [Mb]      ',/,5X,
     &        'SCLXRS  (scat.coeff.)  ',E15.6,
     &        '  >> f         (dim.less)',/,5X,
     &        'V (unit cell)  [a.u.]  ',E15.6,/)
99002 FORMAT (//,60('*'),/,10X,' STOP in <XRAYSPEC>',/,10X,'LCXRAY=',I3,
     &        ' too large for NCSTMAX=',I3)
99003 FORMAT ('       *************************',/,
     &        '       calculating      SPECTRUM',/,
     &        '       *************************',//)
99004 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*          *     *       *****     **     *     *            *'
     &  ,/,10X,
     &  '*           *   *        *    *   *  *     *   *             *'
     &  ,/,10X,
     &  '*            * *         *    *  *    *     * *              *'
     &  ,/,10X,
     &  '*             *     ***  *****   ******      *               *'
     &  ,/,10X,
     &  '*            * *         *  *    *    *      *               *'
     &  ,/,10X,
     &  '*           *   *        *   *   *    *      *               *'
     &  ,/,10X,
     &  '*          *     *       *    *  *    *      *               *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
C
99005 FORMAT (//,' CORE STATES :',//,' NCST:  ',I4,/,' NKPCOR:',20I4)
99006 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99007 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99008 FORMAT (32X,'width:',F12.4,F12.3,/)
99009 FORMAT (22X,I4,F12.6)
99010 FORMAT (/,' number of energies',I5,/,' output format IFMT',2I5,/,
     &        A10,F15.8,/)
99011 FORMAT (' LOOP',I2,I5,' E=',2F6.3,'  IQ=',I2,'  IT=',I2,3X,A)
99012 FORMAT (//,' X-RAY ABSORPTION-RATE CALCULATED FOR ',I4,
     &        ' ENERGY-POINTS   >>> FEIERABEND')
99013 FORMAT (2X,2F11.8,' RYD',F8.4,' EV (REL EF)',2E14.6)
99014 FORMAT (2I2,1X,2E15.8,1X,2E15.8,2X,2E15.8,:,'(',A,')')
99015 FORMAT (2I2,1X,2E15.8,1X,:,'(',A,')')
99016 FORMAT (A10,2I10)
99017 FORMAT (A10,I10,9(2X,A))
99018 FORMAT (A10,A)
99019 FORMAT (5X,'IDIPOL:     ',I5,'  0: RATE=A+B',/,24X,'1: RATE=A',/,
     &        5X,'MODE:',11X,'  2: IM[ SUM M*TAU*M - SUM IRR ]',/,5X,
     &        'MEFORM:',7X,A)
99020 FORMAT (5X,A,':',I12,3X,9(2X,A))
99021 FORMAT (/,4X,' CORE QUANTUM-NUMBERS  FOR  IT=',I2,':   N=',I2,
     &        '  L=',I2,/)
99022 FORMAT (//,5X,70('*'),/,15X,
     &        'WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',/,15X,
     &        '>>>>>>>   2 * Im(E) > TAU(core)   <<<<<<<<',/,15X,
     &        '2 * imaginary part of  E  exceeds the core level width ',
     &        /,15X,'this may lead to arte facts !!!!!!!!!!!!!!',/,5X,
     &        70('*'),/)
99023 FORMAT (5X,'absorption coefficient will be given in ',A,/)
99024 FORMAT (5X,'optical conductivity tensor will be given in ',
     &        'units of  [10^15 1/s]',/)
99025 FORMAT (5X,'   scattering amplitude  f  will be given',
     &        ' dimensionless',/)
99026 FORMAT (5X,A,' will be included in the calculation',/)
99027 FORMAT (5X,'SELFENERGY applied for IT=',I3,/)
99028 FORMAT ('# BUILDBOT: ',A,': ((ABSRTA, ABSRTB),ipol=1,npol)',
     &        ',icst=1,ncst)  for IT =',I5,/,'#',10X,'energy  IE =',I5,
     &        ' E = ',2F10.6,/,(1PE22.14))
99029 FORMAT ('# BUILDBOT: ',A,': ((ABSRTA, ABSRTB),ipol=1,npol)',
     &        ',icst=1,ncst)  for IT =',I5,/,'#',10X,'SIGMA  IESIG =',
     &        I5,' ESIG = ',F10.6,/,(1PE22.14))
      END
