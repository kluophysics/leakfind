C*==compton.f    processed by SPAG 6.70Rc at 15:06 on  6 Apr 2017
      SUBROUTINE COMPTON(TASK,MSSQ,MSST,CALCINT,GETIRRSOL,TSST,TAUQ,
     &                   MEZZ,MEZJ)
C    *******************************************************************
C    *                                                                 *
C    *   TASK:                                                         *
C    *                                                                 *
C    *   COMPTON  -->    Magnetic Compton Scattering                   *
C    *                or Magnetic Momentum Densities                   *
C    *                ( the latter computes 2D cut of mom-density      *
C    *                  at PNMAX if specified, otherwise at PNMAX = 0) *
C    *   POSANI   -->    Positron Annihilation                         *
C    *                                                                 *
C    *   it is assumed that the magnetisation points along the         *
C    *   the global z-axis                                             *
C    *   the spin polarisation refers to that common quantisation axis *
C    *                                                                 *
C    *   revised  Dec. 2016 Misha                                      *
C    *                                                                 *
C    *******************************************************************
C
      USE MOD_RMESH,ONLY:NMMAX,NRMAX,JRWS,R,FULLPOT,JRCRI
      USE MOD_CALCMODE,ONLY:IREL,ORBPOL,DMFT,LDAU,GF_CONV_RH,SOLVER_FP
      USE MOD_ENERGY,ONLY:EFERMI,NEMAX,EIMAG,EILOW,EMAX,EMIN,SPLITSS,
     &    NETAB,IGRID,WETAB,ETAB
      USE MOD_LATTICE,ONLY:ALAT,BDBINV,BBAS
      USE MOD_ANGMOM,ONLY:NL,NMEMAX,NKM,NKMMAX,IND0Q,NKKR,NKMQ
      USE MOD_SITES,ONLY:NQ,QBAS,NQMAX,QMPHI,QMTET,ITOQ,IQAT,NOQ
      USE MOD_TYPES,ONLY:NT,NTMAX,IMT,CTL,CONC,NCPLWFMAX,NCORMAX
      USE MOD_FILES,ONLY:IPRINT,DATSET,LDATSET,IFILCBWF,IOTMP,IOTMP1,
     &    IOTMP2,IOTMP3,IFILBUILDBOT,WRBUILDBOT,FOUND_REAL,
     &    FOUND_REAL_ARRAY,N_FOUND
      USE MOD_CONSTANTS,ONLY:PI,CI,C1,C0
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      USE MOD_MPI,ONLY:MPI,MPI_ID,MPI_ELOOP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTON')
      INTEGER PPOINTS,NPMAX
      PARAMETER (PPOINTS=75,NPMAX=300)
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      CHARACTER*10 TASK
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 A,ADD,AK(3),ALPH(:,:),AP(3),AUX,AUX_D_R,AUX_D_R1,AUX_D_R2,
     &       B,BET(:,:),BP(3),C,CCJJ(:,:),CSQR,DELT(:,:),DELTA,DPN,DPP1,
     &       DPP2,DQVEC(3),EPFE,EPSP,GAM(:,:),KVEC(3),MAXR,MD(:,:),
     &       MD2D(:,:,:),MD2DBND1,MD2DBND2,MD2DCOR1,MD2DCOR2,
     &       MDCOR2D(:,:,:),MDCOR_S(:,:),MDCOR_SQ(:,:,:),MDCOR_ST(:,:),
     &       MD_1QQ(:,:),MD_2QQ(:,:),MD_SQQ(:,:,:,:),PFE,PFEDIR(3),
     &       PFEDIRL(3),PFEMAX,PHI,PNMAX,PNVEC(3),PNVEC0(3),PPMAX,
     &       PPVEC1(3),PPVEC10(3),PPVEC1L,PPVEC2(3),PPVEC20(3),PPVEC2L,
     &       PROTLG(3,3),PRQQ,PVEC(3),PVECDU(3),RWORK(:),TET,TIME1,
     &       TIME2,WGTCOR,WPOS,YVALPA(:),YVALPB(:),ZECH(PPOINTS+1)
      INTEGER API(3),AUX_I_R,I,I1,IA_ERR,IDIMS,IE,IEPATH,
     &        IKMCPLWF_CORT(:,:,:),IKMCPLWF_T(:,:,:),IM,INFO,IOA,IOB,
     &        IPIV(:),IPN,IPNBOT,IPNTOP,IPP1,IPP1BOT,IPP1TOP,IPP2,
     &        IPP2BOT,IPP2TOP,IPROCE(:),IQ,IQA,IQB,IRTOP,IS,IT,ITA,ITB,
     &        I_R,J,J1,K,L,M,MEMODE,N,NCPLWF_CORT(:,:),NCPLWF_T(:,:),NE,
     &        NEPATH,NIN,NP,NPN,NPP,NPP1,NPP2
      COMPLEX*16 CADD,CSUM1(2),CSUM2(2),CSUM3(2),DMAMC(:,:),DMATT(:,:,:)
     &           ,DMATTPOS(NTMAX),DTILT(:,:,:),DTILTPOS(NTMAX),EPRQQ,
     &           ERYD,F_JL_PW(:),F_PW_ZR(:),F_ZL_PW(:),I_JL_PW(:),
     &           I_PW_ZR(:),I_ZL_PW(:),JF(:,:,:),JFL(:,:,:),JFR(:,:,:),
     &           JG(:,:,:),JGL(:,:,:),JGR(:,:,:),MAUX(:,:),ME_IRR(:,:),
     &           ME_PWZR(:,:,:),ME_ZLPW(:,:,:),MSSQL(:,:),P,
     &           RME_COR_P(:,:,:,:),RME_IRR_P(:,:,:,:),
     &           RME_PW_ZR_P(:,:,:,:),RME_ZL_PW_P(:,:,:,:),
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAUA(:,:,:),TAUAB(:,:,:,:),
     &           TAUABPOS(NTMAX,NTMAX),TAUAPOS(NTMAX),TAUK(:,:),
     &           TAUKAB(:,:),TAUKQQL(:,:),TAUQL(:,:),
     &           TAUQQPOS(NQMAX,NQMAX),W1(:,:),WGT,ZF(:,:,:),ZFL(:,:,:),
     &           ZFR(:,:,:),ZG(:,:,:),ZGL(:,:,:),ZGPOS(NRMAX,NTMAX),
     &           ZGR(:,:,:),ZZ
      COMPLEX*16 CJLZ
      REAL*8 DDOT,DNRM2
      CHARACTER*4 EXTENSION
      CHARACTER*10 EXTENSION2D,HEADER_TEXT
      CHARACTER*80 FILNAM
      LOGICAL FOUND,POSANI,SETPNMAX,SHFT,SPEC_DIM_EQ_2D,
     &        SPEC_DIM_EQ_2D_MD2D,USE_ZJOVZZ
      CHARACTER*21 SPEC_TYPE
      CHARACTER*27 TITLE
C
C*** End of declarations rewritten by SPAG
C
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE CCJJ,ALPH,DELT,TAUA,MD,MD_SQQ,IPIV,MAUX,TAUK
      ALLOCATABLE MD_1QQ,MD_2QQ,DMAMC,TAUAB,DMATT,DTILT,ME_PWZR,ME_ZLPW
      ALLOCATABLE MSSQL,TAUQL,W1,TAUKAB,YVALPA,YVALPB,GAM,BET
      ALLOCATABLE RME_PW_ZR_P,TAUKQQL,MD2D,MDCOR2D
      ALLOCATABLE RME_ZL_PW_P,RME_IRR_P,ME_IRR
      ALLOCATABLE NCPLWF_T,IKMCPLWF_T,IPROCE,RWORK
C
      ALLOCATABLE I_PW_ZR,I_ZL_PW
      ALLOCATABLE F_PW_ZR,F_ZL_PW
      ALLOCATABLE I_JL_PW,F_JL_PW
C
      ALLOCATABLE JFL,JGL,ZGL,ZFL,ZG,ZF
      ALLOCATABLE JFR,JGR,ZGR,ZFR,JG,JF
C
      ALLOCATABLE IKMCPLWF_CORT,NCPLWF_CORT,RME_COR_P
      ALLOCATABLE MDCOR_ST,MDCOR_SQ,MDCOR_S
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (IPROCE(NEMAX))
      ALLOCATE (JFL(NRMAX,NCPLWFMAX,NKM),JGL(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGL(NRMAX,NCPLWFMAX,NKM),ZFL(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JFR(NRMAX,NCPLWFMAX,NKM),JGR(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGR(NRMAX,NCPLWFMAX,NKM),ZFR(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (I_PW_ZR(NRMAX))
      ALLOCATE (F_PW_ZR(NRMAX))
      ALLOCATE (I_ZL_PW(NRMAX))
      ALLOCATE (F_ZL_PW(NRMAX))
      ALLOCATE (I_JL_PW(NRMAX))
      ALLOCATE (F_JL_PW(NRMAX))
C
      ALLOCATE (NCPLWF_T(NKMMAX,NTMAX),MDCOR_ST(2,NTMAX))
      ALLOCATE (IKMCPLWF_T(NCPLWFMAX,NKMMAX,NTMAX))
C
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZF')
C
      ALLOCATE (RME_COR_P(NPMAX,NCPLWFMAX,NCORMAX,NTMAX))
      ALLOCATE (IKMCPLWF_CORT(NCPLWFMAX,NCORMAX,NTMAX))
      ALLOCATE (NCPLWF_CORT(NCORMAX,NTMAX))
      RME_COR_P(:,:,:,:) = 0D0
      IKMCPLWF_CORT(:,:,:) = 0
      NCPLWF_CORT(:,:) = 0
C
      ALLOCATE (CCJJ(PPOINTS+1,0:NL),ALPH(PPOINTS,0:NL))
      ALLOCATE (DELT(PPOINTS,0:NL),TAUA(NKMMAX,NKMMAX,NTMAX))
C
      ALLOCATE (IPIV(NKKR),MAUX(NKKR,NKKR))
      ALLOCATE (TAUK(NKKR,NKKR))
      ALLOCATE (MD_1QQ(NQMAX,NQMAX),MD_2QQ(NQMAX,NQMAX))
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (ME_PWZR(NKMMAX,2,NTMAX),ME_ZLPW(NKMMAX,2,NTMAX))
      ALLOCATE (ME_IRR(2,NTMAX))
      ALLOCATE (MSSQL(NKMMAX,NKMMAX),TAUQL(NKMMAX,NKMMAX))
      ALLOCATE (W1(NKMMAX,NKMMAX),TAUKAB(NKMMAX,NKMMAX))
      ALLOCATE (YVALPA(0:NL),YVALPB(0:NL))
      ALLOCATE (GAM(PPOINTS,0:NL),BET(PPOINTS,0:NL))
      ALLOCATE (TAUKQQL(NKMMAX,NKMMAX))
      ALLOCATE (TAUAB(NKMMAX,NKMMAX,NTMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUKQQL')
C
      IF ( IREL.LT.3 ) THEN
         WRITE (6,*) 'SORRY:  <COMPTON> does not work'
         WRITE (6,*) 'in the non/scalar relativistic mode'
         WRITE (6,*) 'use fully relativistic mode instead'
         STOP
      END IF
C
      CALL CPU_TIME(TIME1)
C
      IF ( TASK(1:7).EQ.'COMPTON' ) THEN
         POSANI = .FALSE.
         TITLE = 'Magnetic Compton Scattering'
      ELSE
         POSANI = .TRUE.
         TITLE = 'Positron Annihilation'
      END IF
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
C=======================================================================
C                   set USE_ZJOVZZ
C                   DEFAULT: .FALSE.
C  USE_ZJOVZZ = T:  avoid irregular solution by using MEZJ/MEZZ
C  USE_ZJOVZZ = F:  standard approach
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
      CALL SECTION_FIND_KEYWORD('USE_ZJOVZZ',USE_ZJOVZZ)
C
C=======================================================================
C                   set SPLITSS
C                   DEFAULT: .FALSE.
C=======================================================================
C
      CALL SECTION_FIND_KEYWORD('SPLITSS',SPLITSS)
C
      IF ( USE_ZJOVZZ ) SPLITSS = .FALSE.
C
C=======================================================================
C  MEMODE = 1:  interpolate p-dependent part of radial matrix elements
C  MEMODE = 2:  properly calculate radial matrix elements for each p
C=======================================================================
C
      MEMODE = 1
      CALL SECTION_SET_INTEGER('MEMODE',MEMODE,9999,0)
C
      NP = MIN(250,NPMAX)
C
      ALLOCATE (RME_IRR_P(NPMAX,NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (RME_PW_ZR_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX))
      ALLOCATE (RME_ZL_PW_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZF')
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C     Set up rotation matrix to rotate momentum vector from global
C     to local frame
C
      A = QMPHI(1)*PI/180.D0
      B = QMTET(1)*PI/180.D0
      PROTLG(1,1) = COS(A)*COS(B)
      PROTLG(1,2) = SIN(A)*COS(B)
      PROTLG(1,3) = -SIN(B)
      PROTLG(2,1) = -SIN(A)
      PROTLG(2,2) = COS(A)
      PROTLG(2,3) = 0.D0
      PROTLG(3,1) = COS(A)*SIN(B)
      PROTLG(3,2) = SIN(A)*SIN(B)
      PROTLG(3,3) = COS(B)
C
C======================================================================
C
      C = CTL(1,1)
      CSQR = C**2
C
      CALL READTAU(9,ERYD,0,NE,W1,0,NT,.FALSE.,TAUQ,MSSQ,0,NQ,0,NKMMAX,
     &             NKMMAX,IPRINT)
C
      IF ( NE.GT.NEMAX ) CALL STOP_MESSAGE(ROUTINE,'NE > NEMAX')
C
      CLOSE (IFILCBWF)
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=(4+2*(2+4*NRMAX))*2*8)
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( SPLITSS ) THEN
         NEPATH = 2
         DO IEPATH = 1,NEPATH
            CALL EPATH(IGRID(IEPATH),EMIN,EMAX,EIMAG,NETAB(IEPATH),
     &                 ETAB(1,IEPATH),WETAB(1,IEPATH),EILOW,IPRINT,
     &                 NEMAX)
         END DO
      ELSE
         NEPATH = 1
      END IF
C
C-----------------------------------------------------------------------
C         Set up the normalized right-hand triple of vectors
C         normal PNVEC and perpendicular (in plane) ones PPVEC1 & PPVEC2
C         according to the input file:
C         if PNVEC is specified explicitly, read and normalize
C            check whether PPVEC1 is also specified: if so
C               read it, normalize, construct PPVEC2 from PNVEC & PPVEC1
C            otherwise
C               guess PPVEC1 & PPVEC2 according to the given PNVEC.
C         otherwise
C            check if PPVEC1 and PPVEC2 are specified in the input file
C               if both found, read and normalize
C               if non is found use the default triple
C               if only one is specified, then quit with error status
C-----------------------------------------------------------------------
C
C-------------------------------------------------------- look for PNVEC
C
      CALL SECTION_SET_REAL_ARRAY('PNVEC',PNVEC,N_FOUND,3,0,9999D0,0)
C
      IF ( FOUND_REAL_ARRAY ) THEN
C
         AUX = 1D0/DNRM2(3,PNVEC,1)
         CALL DSCAL(3,AUX,PNVEC,1)
C
C-------------------------------------- PNVEC specified. look for PPVEC1
C
         CALL SECTION_SET_REAL_ARRAY('PPVEC1',PPVEC1,NIN,3,0,9999D0,0)
C
         IF ( FOUND_REAL_ARRAY ) THEN
C------------------------------------------------------ PPVEC1 specified
C
            AUX = 1D0/DNRM2(3,PPVEC1,1)
            CALL DSCAL(3,AUX,PPVEC1,1)
C
C----------- constracr PPVEC2. ignore PPVEC2 on input even if it exists!
C----------------------------------------------- PPVEC2 = PNVEC x PPVEC1
            CALL RVECXPRO(PNVEC,PPVEC1,PPVEC2)
C
            AUX = 1D0/DNRM2(3,PPVEC2,1)
            CALL DSCAL(3,AUX,PPVEC2,1)
C
         ELSE
C-------------- PPVEC1 not specified. ignore PPVEC2 even if it specified
C---------------------------- construct/guess PPVEC1 & PPVEC2 from PNVEC
C
            PPVEC1(1:3) = 0D0
            PPVEC2(1:3) = 0D0
C
            TET = ACOS(PNVEC(3))
C
C--------------------------------- treat TET = 0 and TET = PI separately
C
            IF ( ABS(TET).LT.1D-14 ) THEN
               PPVEC1(1) = 1D0
               PPVEC2(2) = 1D0
            ELSE IF ( ABS(ABS(TET)-PI).LT.1D-14 ) THEN
               PPVEC1(1) = -1D0
               PPVEC2(2) = 1D0
            ELSE
               PHI = ACOS(PNVEC(1)/SQRT(PNVEC(1)**2+PNVEC(2)**2))
               IF ( PNVEC(2).LT.0 ) PHI = -PHI
C
               PPVEC1(1) = SIN(TET+PI/2D0)*COS(PHI)
               PPVEC1(2) = SIN(TET+PI/2D0)*SIN(PHI)
               PPVEC1(3) = COS(TET+PI/2D0)
               PPVEC2(1) = COS(PHI+PI/2D0)
               PPVEC2(2) = SIN(PHI+PI/2D0)
               PPVEC2(3) = 0D0
            END IF
         END IF
C
      ELSE
C---------------------- PNVEC is not specified. look for PPVEC1 & PPVEC2
C
         CALL SECTION_SET_REAL_ARRAY('PPVEC1',PPVEC1,NIN,3,0,9999D0,0)
         FOUND = FOUND_REAL_ARRAY
C
         CALL SECTION_SET_REAL_ARRAY('PPVEC2',PPVEC2,NIN,3,0,9999D0,0)
C
         IF ( FOUND .AND. FOUND_REAL_ARRAY ) THEN
C
            AUX = 1D0/DNRM2(3,PPVEC1,1)
            CALL DSCAL(3,AUX,PPVEC1,1)
            AUX = 1D0/DNRM2(3,PPVEC2,1)
            CALL DSCAL(3,AUX,PPVEC2,1)
C
C----------------------------- construct PNVEC:  PNVEC = PPVEC1 x PPVEC2
            CALL RVECXPRO(PPVEC1,PPVEC2,PNVEC)
C
            AUX = 1D0/DNRM2(3,PNVEC,1)
            CALL DSCAL(3,AUX,PNVEC,1)
C
         ELSE IF ( .NOT.FOUND .AND. .NOT.FOUND_REAL_ARRAY ) THEN
C------------------------------------- non of the vectors are specified!
C------------------------------------------- take the default DIRECTIONS
            PNVEC(1) = 0D0
            PNVEC(2) = 0D0
            PNVEC(3) = 1D0
C
            PPVEC1(1) = 1D0
            PPVEC1(2) = 0D0
            PPVEC1(3) = 0D0
C
            PPVEC2(1) = 0D0
            PPVEC2(2) = 1D0
            PPVEC2(3) = 0D0
C
            FOUND = .TRUE.
         ELSE
            WRITE (6,*) 'SORRY: incomplete input for DIRECTIONS!!!'
            WRITE (6,*) 'NONE, PNVEC, PNVEC \& PPVEC1, PPVEC1 \& PPVEC2'
            WRITE (6,*) 'have to be specified in the input file!!!'
            STOP
         END IF
      END IF
C
      PPVEC10(:) = PPVEC1(:)
      PPVEC20(:) = PPVEC2(:)
      PNVEC0(:) = PNVEC(:)
C
      CALL SECTION_FIND_KEYWORD('MD2D',SPEC_DIM_EQ_2D_MD2D)
C
      IF ( SPEC_DIM_EQ_2D_MD2D ) THEN
         SPEC_DIM_EQ_2D = .TRUE.
      ELSE
         CALL SECTION_FIND_KEYWORD('CP2D',SPEC_DIM_EQ_2D)
      END IF
C
      CALL SECTION_SET_REAL('PNMAX',PNMAX,5D0,0)
      SETPNMAX = FOUND_REAL
C
      CALL SECTION_SET_INTEGER('NPN',NPN,50,0)
C
      IF ( NPN.GT.1 ) THEN
         DPN = PNMAX/DBLE(NPN-1)
      ELSE
         DPN = PNMAX
         NPN = 1
      END IF
C
C--------------------------------------------- perpendicular directions:
C
      CALL SECTION_SET_REAL('PPMAX',PPMAX,10D0,0)
C
      CALL SECTION_SET_INTEGER('NPP',NPP,50,0)
C
      NPP1 = NPP
      NPP2 = NPP
C
      IF ( NPP.GT.1 ) THEN
         DPP1 = PPMAX/DBLE(NPP1-1)
         DPP2 = PPMAX/DBLE(NPP2-1)
      ELSE
         DPP1 = PPMAX
         DPP2 = PPMAX
         NPP = 1
         NPP1 = 1
         NPP2 = 1
      END IF
C
C-----------------------------------------------------------------------
C
      CALL DSCAL(3,DPN,PNVEC,1)
      CALL DSCAL(3,DPP1,PPVEC1,1)
      CALL DSCAL(3,DPP2,PPVEC2,1)
C
      PPVEC1L = 0.0D0
      PPVEC2L = 0.0D0
C
C   special treatment of the integration limits if SHFT (shift) is given
C
      CALL SECTION_FIND_KEYWORD('SHFT',SHFT)
C
C------------------------------------- fix limits for K-loops PN,PP1,PP2
C---------------------------------------- allocations will be done later
C
      IF ( SHFT ) THEN
         IPP1TOP = INT((NPP1-1)/2)
         IPP1BOT = -IPP1TOP
         IPP2TOP = INT((NPP2-1)/2)
         IPP2BOT = -IPP2TOP
         IPNTOP = INT((NPN-1)/2)
         IPNBOT = -IPNTOP
      ELSE
         IPP1TOP = NPP1 - 1
         IPP1BOT = 0
         IPP2TOP = NPP2 - 1
         IPP2BOT = 0
         IPNTOP = NPN - 1
         IPNBOT = 0
      END IF
C
C------------------------------------- fix limits in case of CP2D & MD2D
      IF ( SPEC_DIM_EQ_2D ) THEN
C
C----------------------------------- distinguish   COMPTON - CP2D & MD2D
C----------------------------------------------------- and   POSANI MD2D
C
         IF ( TASK(1:7).EQ.'COMPTON' .AND. SPEC_DIM_EQ_2D_MD2D ) THEN
C
C-------------------------------- we take 2D mom-density cut at PNMAX if
C----------------------- latter was specified explicitly, zero otherwise
            NPN = 1
            IPNTOP = 1
            IPNBOT = 1
            IF ( .NOT.SETPNMAX ) PNMAX = 1D-6
            CALL DSCAL(3,PNMAX,PNVEC,1)
         END IF
C
         PPVEC1L = DNRM2(3,PPVEC1,1)
         PPVEC2L = DNRM2(3,PPVEC2,1)
C
      END IF
C
C-----------------------------------------------------------------------
C
      CALL COMPTON_CHECKVEC(PPVEC1,PPVEC2,PNVEC)
C
      WRITE (6,99001) TITLE,EFERMI,MEMODE,GF_CONV_RH,SPLITSS,USE_ZJOVZZ,
     &                NE,PPVEC10,PPVEC20,PNVEC0,NPP1,PPVEC1,PPMAX,DPP1,
     &                NPP2,PPVEC2,PPMAX,DPP2,NPN,PNVEC,PNMAX,DPN
C
      IF ( SPEC_DIM_EQ_2D ) WRITE (6,99003) SHFT,PNVEC,PPVEC1L,PPVEC2L
C
C=======================================================================
      IF ( POSANI ) THEN
         CALL POSANIREAD(DATSET,LDATSET,NT,IPRINT,JRWS,IMT,TAUAPOS,
     &                   DMATTPOS,DTILTPOS,TAUABPOS,TAUQQPOS,ZGPOS,
     &                   NQMAX,NTMAX,NRMAX,NMMAX)
C=======================================================================
C             Compton scattering: initialize variables with dummy values
      ELSE
         CALL CINIT(NTMAX,TAUAPOS)
         CALL CINIT(NTMAX,DMATTPOS)
         CALL CINIT(NTMAX,DTILTPOS)
         CALL CINIT(NTMAX*NTMAX,TAUABPOS)
         CALL CINIT(NQMAX*NQMAX,TAUQQPOS)
         CALL CINIT(NRMAX*NTMAX,ZGPOS)
      END IF
C=======================================================================
C
      IF ( IPRINT.GE.0 ) THEN
         DO J = 1,3
            WRITE (6,99002) J,(BBAS(I,J),I=1,3)
         END DO
C
      END IF
C
      IF ( SHFT ) THEN
         PFEMAX = SQRT((PNMAX/2)**2+0.5D0*PPMAX**2)
      ELSE
         PFEMAX = SQRT(PNMAX**2+2.D0*PPMAX**2)
      END IF
C
C
      MAXR = 0.D0
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
         MAXR = MAX(MAXR,R(IRTOP,IM))
      END DO
      DELTA = PFEMAX*MAXR/DBLE(PPOINTS)
C
      WRITE (6,99009) MAXR,PFEMAX,DELTA
C
      DO K = 1,PPOINTS + 1
         ZECH(K) = DELTA*DBLE(K-1)
         ZZ = DCMPLX(ZECH(K))
         DO L = 0,NL
            CCJJ(K,L) = DBLE(CJLZ(L,ZZ))
         END DO
      END DO
C
      DO L = 0,NL
         CALL SPLINEL(CCJJ(1,L),ZECH,PPOINTS+1,YVALPA,YVALPB,.FALSE.,L,
     &                ALPH,BET,GAM,DELT,NL)
      END DO
C
C-------------- Now one can allocate arrays with the correct dimensions
C
      ALLOCATE (MD(IPNBOT:IPNTOP,2))
      ALLOCATE (MD_SQQ(IPNBOT:IPNTOP,2,NQMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MD_SQQ')
C
      MD(-IPNBOT:IPNTOP,1:2) = 0D0
      MD_SQQ(IPNBOT:IPNTOP,1:2,1:NQMAX,1:NQMAX) = 0D0
C
      ALLOCATE (MDCOR_S(IPNBOT:IPNTOP,2))
      ALLOCATE (MDCOR_SQ(IPNBOT:IPNTOP,2,NQMAX))
C
      MDCOR_S(:,:) = 0D0
      MDCOR_SQ(:,:,:) = 0D0
C
      IF ( SPEC_DIM_EQ_2D ) THEN
         ALLOCATE (MDCOR2D(IPP1BOT:IPP1TOP,IPP2BOT:IPP2TOP,2))
         ALLOCATE (MD2D(IPP1BOT:IPP1TOP,IPP2BOT:IPP2TOP,2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MD2D')
C
         MDCOR2D(:,:,:) = 0D0
         MD2D(:,:,:) = 0D0
      END IF
C
C
C=======================================================================
C                               CORE STATES
C
C      Step 0: calculate |p|-dependent part for grid NP, PFEMAX
C=======================================================================
C
      CALL COMPTON_CORE0(C,NCPLWF_CORT,IKMCPLWF_CORT,TASK,NP,PFEMAX,
     &                   ZGPOS,RME_COR_P,ZECH,ALPH,BET,GAM,DELT,NPMAX,
     &                   PPOINTS)
C
C=======================================================================
C
C            PREPARATIONS DONE - now run loops over:
C
C             E-paths, E-loop and k-loops PN,PP1,PP2
C
C=======================================================================
C
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP PATH
      LOOP_IEPATH:DO IEPATH = 1,NEPATH
C
         IF ( DMFT .OR. LDAU )
     &        CALL DMFT_READSIG(NETAB(IEPATH),ETAB(1,IEPATH),IPRINT)
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
         MPI_ELOOP = MPI
         CALL MPI_DISTRIBUTE(IPROCE,NETAB(IEPATH),MPI_ELOOP,'E')
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         LOOP_IE:DO IE = 1,NETAB(IEPATH)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ID.NE.IPROCE(IE) ) CYCLE LOOP_IE
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IF ( IEPATH.EQ.1 ) THEN
C========read site-diagonal TAU- and inverse eff. single site t-matrix m
C                                        for present energy ERYD (KEY=2)
               DO IQ = 1,NQ
C
                  ERYD = ETAB(IE,1)
C
                  CALL READTAU(9,ERYD,IE,NE,W1,0,NT,.FALSE.,TAUQ(1,1,IQ)
     &                         ,MSSQ(1,1,IQ),IQ,NQ,2,NKMQ(IQ),NKMMAX,
     &                         IPRINT)
C
               END DO
               ERYD = ETAB(IE,1)
               WRITE (6,'(10X,A,I10,2F10.5)') 'energy   ',IE,ERYD
            ELSE
               ERYD = ETAB(IE,2)
               WRITE (6,'(10X,A,I10,2F10.5)') 'energy   ',IE,ERYD
            END IF
C
            IF ( DIMAG(ERYD).GT.0D0 ) GETIRRSOL = .TRUE.
C
C--------------------------------calculate WF and DOS - radial integrals
C
            IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &           = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
            IF ( (DMFT .OR. LDAU) .AND. IEPATH.EQ.2 )
     &           DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &           = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE))
     &           ,0.0D0)
C
            CALL RUNSSITE(CALCINT,1,1,IFILCBWF,GETIRRSOL,ERYD,P,IPRINT,
     &                    TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            WRITE (6,'(10X,A)') '<SSITE> - run completed '
C
C11111111111111111111111111111111111111111111111111111111111111111111111
            IF ( IEPATH.EQ.1 ) THEN
C----------------------------------- set up the matrices TAUQL and MSSQL
C----------------------------------- that are referred to the local axis
C----- the resulting projection matrices DMATT and DTILT, as well as the
C------------------ matrices TAUA and TAUAB also refer to the local axis
C-------------------------  TAUQ(IQ) for equivalent sites are the same !
C
               M = NKMMAX
               N = NKMQ(1)
C
               DO ITA = 1,NT
                  IQ = IQAT(1,ITA)
C
                  IF ( NKMQ(IQ).NE.N ) THEN
                     WRITE (6,*) 'NKMQ(1) = N: ',N
                     WRITE (6,*) 'NKMQ(IQ):    ',NKMQ(IQ)
                     WRITE (6,*) 'IQ  IT:      ',IT,IQ
                     CALL STOP_MESSAGE(ROUTINE,'NKMQ assumed the same'//
     &                                 ' for all sites  IQ')
                  END IF
C
                  DO J = 1,N
                     CALL ZCOPY(N,MSSQ(1,J,IQ),1,MSSQL(1,J),1)
                     CALL ZCOPY(N,TAUQ(1,J,IQ),1,TAUQL(1,J),1)
                  END DO
C
                  CALL GETDMAT(TAUQL,DMATT(1,1,ITA),DTILT(1,1,ITA),
     &                         DMAMC,N,MSSQL,MSST(1,1,ITA),M)
C
                  CALL ZGEMM('N','N',N,N,N,C1,DMATT(1,1,ITA),M,TAUQL,M,
     &                       C0,TAUA(1,1,ITA),M)
C
               END DO
C
               DO ITA = 1,NT
                  DO ITB = 1,NT
C
                     CALL ZGEMM('N','N',N,N,N,C1,TAUA(1,1,ITA),M,
     &                          DTILT(1,1,ITB),M,C0,TAUAB(1,1,ITA,ITB),
     &                          M)
C
                  END DO
               END DO
C
C-----------------------------------------------------------------------
C     subtract single site t-matrix in case of split energy path
C-----------------------------------------------------------------------
C
               IF ( NEPATH.EQ.2 ) THEN
                  DO IT = 1,NT
                     IQ = IQAT(1,IT)
                     DO J = 1,NKMQ(IQ)
                        DO I = 1,NKMQ(IQ)
                           TAUA(I,J,IT) = TAUA(I,J,IT) - TSST(I,J,IT)
                        END DO
                     END DO
                  END DO
               END IF
C
C-----------------------------------------------------------------------
C     account for irregular solution using      MEZJ/MEZZ
C-----------------------------------------------------------------------
C
               IF ( USE_ZJOVZZ ) THEN
                  DO IT = 1,NT
                     IQ = IQAT(1,IT)
                     DO I = 1,NKMQ(IQ)
                        TAUA(I,I,IT) = TAUA(I,I,IT) - MEZJ(I,I,IT,1)
     &                                 /MEZZ(I,I,IT,1)
                     END DO
                  END DO
               END IF
C
C11111111111111111111111111111111111111111111111111111111111111111111111
            ELSE
C22222222222222222222222222222222222222222222222222222222222222222222222
C
               DO IT = 1,NT
                  IQ = IQAT(1,IT)
                  DO J = 1,NKMQ(IQ)
                     DO I = 1,NKMQ(IQ)
                        TAUA(I,J,IT) = TSST(I,J,IT)
                     END DO
                  END DO
               END DO
            END IF
C22222222222222222222222222222222222222222222222222222222222222222222222
C
            WRITE (6,'(10X,A)') 'all matrices set up  -  enter k-loop '
C
C     calculate relativistic momentum
C
            P = SQRT(ERYD*(1.0D0+ERYD/CTL(1,1)**2))
C
C     -------------- calculate energy - dependent terms of str.constants
C
            CALL STRCC(ERYD,.FALSE.)
C
C=======================================================================
C  MEMODE = 1:  interpolate |p|-dependent part of radial matrix elements
C               Step 0: calculate |p|-dependent part for grid NP, PFEMAX
C=======================================================================
C
            IF ( MEMODE.EQ.1 ) CALL COMPTONME0(C,ZG,ZF,JG,JF,I_PW_ZR,
     &           I_ZL_PW,I_JL_PW,F_PW_ZR,F_ZL_PW,F_JL_PW,ZGR,ZFR,JGR,
     &           JFR,ZGL,ZFL,JGL,JFL,NCPLWF_T,IKMCPLWF_T,TASK,NP,PFEMAX,
     &           ZGPOS,RME_PW_ZR_P,RME_ZL_PW_P,RME_IRR_P,ZECH,ALPH,BET,
     &           GAM,DELT,NPMAX,PPOINTS)
C
C=======================================================================
C
C     KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C     KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C   to get numerically more stable results, it is in general recommended
C   to add values in an ascending order (smalls first)
C
C---------------------------------------------- determine the max radius
            AUX_D_R = SQRT(DBLE(IPNTOP**2+IPP1TOP**2+IPP2TOP**2))
            AUX_I_R = NINT(AUX_D_R)
C
            WRITE (6,'(10X,A,I5)') 'radius = ',AUX_I_R
C
            LOOP_I_R:DO I_R = AUX_I_R,0, - 1
C
               IF ( MOD(I_R,5).EQ.0 ) WRITE (6,'(10X,A,I5,A)')
     &               'outer p-loop: ',I_R,' until p-points done'
C
               LOOP_IPN:DO IPN = IPNBOT,IPNTOP
C
C----------------------------------------- outside the integration shell
C
                  IF ( IPN.GT.I_R ) CYCLE LOOP_IPN
C
                  AUX_D_R = DBLE(IPN**2)
C
                  LOOP_IPP2:DO IPP2 = IPP2BOT,IPP2TOP
C
                     AUX_D_R2 = AUX_D_R + IPP2**2
C
C----------------------------------------- outside the integration shell
C
                     IF ( SQRT(AUX_D_R2).GT.I_R ) CYCLE LOOP_IPP2
C
                     LOOP_IPP1:DO IPP1 = IPP1BOT,IPP1TOP
C
                        AUX_D_R1 = SQRT(AUX_D_R2+IPP1**2)
C
C----------------------------------------- outside the integration shell
C
                        IF ( (AUX_D_R1.GT.I_R) .OR. 
     &                       (AUX_D_R1.LE.(I_R-1)) ) CYCLE LOOP_IPP1
C
                        IF ( SPEC_DIM_EQ_2D ) MD(:,:) = 0D0
C
                        DO I = 1,3
                           PVEC(I) = PPVEC1(I)*IPP1 + PPVEC2(I)
     &                               *IPP2 + PNVEC(I)*IPN
                        END DO
C
                        PFE = DNRM2(3,PVEC,1)
C
                        IF ( DABS(PFE).LT.1.D-7 ) THEN
                           PVEC(:) = 1.0D-7
                           PFE = DNRM2(3,PVEC,1)
                        END IF
C
                        PFEDIR(1:3) = PVEC(1:3)/PFE
C
C ----------------------------- Rotate PFEDIR from global to local frame
C
                        DO I = 1,3
                           PFEDIRL(I) = 0.D0
                           DO J = 1,3
                              PFEDIRL(I) = PFEDIRL(I) + PROTLG(I,J)
     &                           *PFEDIR(J)
                           END DO
                        END DO
C
                        EPFE = (CSQR/2D0)
     &                         *(SQRT(1D0+4D0*PFE**2/CSQR)-1D0)
C
                        EPSP = 4D0*PI*SQRT((EPFE+CSQR)/(2D0*EPFE+CSQR))
C
C=======================================================================
C                     calculate matrix elments
C=======================================================================
C
                        IF ( MEMODE.EQ.1 ) THEN
C
C-----------------------------------------------------------------------
C  MEMODE = 1:  interpolate |p|-dependent part of radial matrix elements
C               Step 1: calculate remaining ->p - dependent part
C-----------------------------------------------------------------------
C
                           CALL COMPTONME1(NCPLWF_T,IKMCPLWF_T,NP,PFE,
     &                        PFEMAX,ME_PWZR,ME_ZLPW,ME_IRR,RME_PW_ZR_P,
     &                        RME_ZL_PW_P,RME_IRR_P,PFEDIR,EPSP,NPMAX)
C
C-----------------------------------------------------------------------
C  MEMODE = 2:  calculate matrix elements without using interpolation
C-----------------------------------------------------------------------
C
                        ELSE
C
                           CALL STOP_MESSAGE(ROUTINE,
     &                        'MEMODE = 2 not available yet')
C
                        END IF
C
C=======================================================================
C                          CORE CONTRIBUTION
C=======================================================================
C
                        IF ( IEPATH.EQ.1 .AND. IE.EQ.1 )
     &                       CALL COMPTON_CORE1(NCPLWF_CORT,
     &                       IKMCPLWF_CORT,NP,PFE,PFEMAX,RME_COR_P,
     &                       MDCOR_ST,PFEDIR,EPSP,NPMAX)
C
C=======================================================================
C
                        DO I = 1,3
                           PVECDU(I) = PVEC(I)*ALAT/(2*PI)
                        END DO
C
                        DO I = 1,3
                           BP(I) = DDOT(3,BBAS(1,I),1,PVECDU,1)
                        END DO
C
                        DO I = 1,3
                           AP(I) = 0D0
                           DO J = 1,3
                              AP(I) = AP(I) + BDBINV(I,J)*BP(J)
                           END DO
                        END DO
C
                        DO I = 1,3
                           API(I) = NINT(AP(I))
                           AK(I) = AP(I) - API(I)
                        END DO
C
                        DO I = 1,3
                           KVEC(I) = 0D0
                           DO J = 1,3
                              KVEC(I) = KVEC(I) + BBAS(I,J)*AK(J)
                           END DO
                        END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
                        IF ( IEPATH.EQ.1 ) THEN
C
                           CALL STRSET(IPN,KVEC,MAUX,TAUK,P)
C
C
                           CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,
     &                                 NKKR,NKMMAX)
C
                           CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
                           CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,
     &                                 NKKR*NKKR,INFO)
C
                        END IF
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ SITE A
                        LOOP_IQA:DO IQA = 1,NQ
                           LOOP_IOA:DO IOA = 1,NOQ(IQA)
                              ITA = ITOQ(IOA,IQA)
C
                              N = NKMQ(IQA)
                              M = NKMMAX
C
C
C     ------------------------------------------------------------------
C     TERM 1
                              CALL COMPTON_SUM(ME_PWZR,ITA,TAUA(1,1,ITA)
     &                           ,ME_ZLPW,ITA,N,M,CSUM1)
C
                              IF ( POSANI ) THEN
                                 WPOS = -DIMAG(TAUAPOS(ITA))
                              ELSE
                                 WPOS = 1D0
                              END IF
C
                              WGT = CONC(ITA)*WETAB(IE,IEPATH)*WPOS
C
                              DO IS = 1,2
                                 CADD = WGT*CSUM1(IS)
                                 ADD = -DIMAG(CADD)/PI
C
                                 MD(IPN,IS) = MD(IPN,IS) + ADD
                                 MD_SQQ(IPN,IS,IQA,IQA)
     &                              = MD_SQQ(IPN,IS,IQA,IQA) + ADD
                              END DO
C
C     ------------------------------------------------------------------
C                                core contribution
C     ------------------------------------------------------------------
                              IF ( IEPATH.EQ.1 .AND. IE.EQ.1 ) THEN
C
                                 WGTCOR = CONC(ITA)*WPOS
                                 DO IS = 1,2
                                    ADD = WGTCOR*MDCOR_ST(IS,ITA)
C
                                    MDCOR_S(IPN,IS) = MDCOR_S(IPN,IS)
     &                                 + ADD
                                    MDCOR_SQ(IPN,IS,IQA)
     &                                 = MDCOR_SQ(IPN,IS,IQA) + ADD
                                 END DO
C
                              END IF
C
C     ------------------------------------------------------------------
C                     add site-diagonal irregular contribution
C     ------------------------------------------------------------------
                              IF ( .NOT.SPLITSS .AND. .NOT.USE_ZJOVZZ )
     &                             THEN
C
                                 DO IS = 1,2
                                    CADD = WGT*ME_IRR(IS,ITA)
                                    ADD = -DIMAG(CADD)/PI
C
                                    MD(IPN,IS) = MD(IPN,IS) + ADD
                                    MD_SQQ(IPN,IS,IQA,IQA)
     &                                 = MD_SQQ(IPN,IS,IQA,IQA) + ADD
                                 END DO
C
                              END IF
C
C     ------------------------------------------------------------------
C
C=======================================================================
                              IF ( IEPATH.EQ.1 ) THEN
                                 DO IOB = 1,NOQ(IQA)
                                    ITB = ITOQ(IOB,IQA)
C
C     ---------------------------------------------------------------
C     TERM 2
                                    CALL COMPTON_SUM(ME_PWZR,ITA,
     &                                 TAUAB(1,1,ITA,ITB),ME_ZLPW,ITB,N,
     &                                 M,CSUM2)
C
                                    IF ( POSANI ) THEN
                                       WPOS = -DIMAG(TAUABPOS(ITA,ITB))
                                    ELSE
                                       WPOS = 1D0
                                    END IF
C
                                    WGT = CONC(ITA)*CONC(ITB)
     &                                 *WETAB(IE,IEPATH)*WPOS
C
                                    DO IS = 1,2
                                       CADD = -WGT*CSUM2(IS)
                                       ADD = -DIMAG(CADD)/PI
C
                                       MD(IPN,IS) = MD(IPN,IS) + ADD
                                       MD_SQQ(IPN,IS,IQA,IQA)
     &                                    = MD_SQQ(IPN,IS,IQA,IQA) + ADD
                                    END DO
C     ------------------------------------------------------------------
C
                                 END DO
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ SITE B
                                 LOOP_IQB:DO IQB = 1,NQ
C
                                    DQVEC(1:3) = QBAS(1:3,IQA)
     &                                 - QBAS(1:3,IQB)
                                    PRQQ = DDOT(3,PVEC,1,DQVEC,1)
                                    EPRQQ = EXP(-CI*PRQQ*ALAT)
C
                                    I1 = IND0Q(IQA) + 1
C
C ------------------------------ get the matrix TAUKQQ = TAU(k)[IQA,IQB]
C
                                    DO J = 1,N
                                       J1 = IND0Q(IQB) + J
                                       CALL ZCOPY(N,TAUK(I1,J1),1,
     &                                    TAUKQQL(1,J),1)
                                    END DO
C
                                    LOOP_IOB:DO IOB = 1,NOQ(IQB)
C
                                       ITB = ITOQ(IOB,IQB)
C
C     ------------------------------------------------------------------
C     TERM 3
                                       CALL ZGEMM('N','N',N,N,N,C1,
     &                                    DMATT(1,1,ITA),M,TAUKQQL,M,C0,
     &                                    W1,M)
C
                                       CALL ZGEMM('N','N',N,N,N,C1,W1,M,
     &                                    DTILT(1,1,ITB),M,C0,TAUKAB,M)
C
C
                                       CALL COMPTON_SUM(ME_PWZR,ITA,
     &                                    TAUKAB,ME_ZLPW,ITB,N,M,CSUM3)
C
C
                                       IF ( POSANI ) THEN
                                         WPOS = -
     &                                      DIMAG(DMATTPOS(ITA)*TAUQQPOS
     &                                      (IQA,IQB)*DTILTPOS(ITB))
                                       ELSE
                                         WPOS = 1D0
                                       END IF
C
                                       WGT = EPRQQ*CONC(ITA)*CONC(ITB)
     &                                    *WETAB(IE,IEPATH)*WPOS
C
                                       DO IS = 1,2
                                         CADD = WGT*CSUM3(IS)
                                         ADD = -DIMAG(CADD)/PI
C
                                         MD(IPN,IS) = MD(IPN,IS) + ADD
                                         MD_SQQ(IPN,IS,IQA,IQB)
     &                                      = MD_SQQ(IPN,IS,IQA,IQB)
     &                                      + ADD
                                       END DO
C     ------------------------------------------------------------------
C
                                    END DO LOOP_IOB
                                 END DO LOOP_IQB
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ SITE B
C
                              END IF
C=======================================================================
C
                           END DO LOOP_IOA
                        END DO LOOP_IQA
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ SITE A
C
C
                        IF ( SPEC_DIM_EQ_2D ) THEN
                           DO IS = 1,2
                              MDCOR2D(IPP1,IPP2,IS)
     &                           = MDCOR2D(IPP1,IPP2,IS)
     &                           + MDCOR_S(IPN,IS)
                              MD2D(IPP1,IPP2,IS) = MD2D(IPP1,IPP2,IS)
     &                           + MD(IPN,IS)
                           END DO
                        END IF
C
                     END DO LOOP_IPP1
                  END DO LOOP_IPP2
               END DO LOOP_IPN
            END DO LOOP_I_R
C     KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C     KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         END DO LOOP_IE
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      END DO LOOP_IEPATH
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP PATH
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C             collect results for energy points of a E-path
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         IF ( .NOT.SPEC_DIM_EQ_2D ) THEN
C
            IDIMS = (IPNTOP-IPNBOT+1)*2*NQMAX*NQMAX
            ALLOCATE (RWORK(IDIMS))
C
            RWORK(:) = 0D0
            IDIMS = (IPNTOP-IPNBOT+1)*2
            CALL DRV_MPI_REDUCE_R(MD(IPNBOT,1),RWORK,IDIMS)
C
            RWORK(:) = 0D0
            IDIMS = (IPNTOP-IPNBOT+1)*2*NQMAX*NQMAX
            CALL DRV_MPI_REDUCE_R(MD_SQQ(IPNBOT,1,1,1),RWORK,IDIMS)
C
            RWORK(:) = 0D0
            IDIMS = (IPNTOP-IPNBOT+1)*2*NQMAX
            CALL DRV_MPI_REDUCE_R(MDCOR_SQ(IPNBOT,1,1),RWORK,IDIMS)
C
         ELSE
C
            IDIMS = (IPP1TOP-IPP1BOT+1)*(IPP2TOP-IPP2BOT+1)*2
            ALLOCATE (RWORK(IDIMS))
            CALL DRV_MPI_REDUCE_R(MD2D(IPP1BOT,IPP2BOT,1),RWORK,IDIMS)
C
            CALL DRV_MPI_REDUCE_R(MDCOR2D(IPP1BOT,IPP2BOT,1),RWORK,
     &                            IDIMS)
C
         END IF
C
         CALL DRV_MPI_BARRIER
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C=======================================================================
C                         write spectra
C=======================================================================
C
C1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D
         IF ( .NOT.SPEC_DIM_EQ_2D ) THEN
C
C--------------------- distinguish   COMPTON - CP2D & MD2D
C----------------------------- and   POSANI MD2D
C
            IF ( TASK(1:7).EQ.'COMPTON' ) THEN
               EXTENSION = '.mcp'
               SPEC_TYPE = 'Compton Scattering'
               HEADER_TEXT = 'COMPTON   '
            ELSE IF ( TASK(1:6).EQ.'POSANI' ) THEN
               EXTENSION = '.mpa'
               SPEC_TYPE = 'Positron Annihilation'
               HEADER_TEXT = 'POSITRON  '
            END IF
C
C------------------------------------------------------------------- BND
            FILNAM = DATSET(1:LDATSET)//'_BND'//EXTENSION
            OPEN (UNIT=IOTMP1,FILE=FILNAM)
            WRITE (6,99011) SPEC_TYPE,IOTMP1,FILNAM
            CALL WRHEAD(IOTMP1,FILNAM,HEADER_TEXT,NETAB(1))
C
            WRITE (IOTMP1,99006) PNVEC0
            WRITE (IOTMP1,99004) PNMAX
            WRITE (IOTMP1,99005) NPN
C
C------------------------------------------------------------------- COR
            FILNAM = DATSET(1:LDATSET)//'_COR'//EXTENSION
            OPEN (UNIT=IOTMP2,FILE=FILNAM)
            WRITE (6,99011) SPEC_TYPE,IOTMP2,FILNAM
            CALL WRHEAD(IOTMP2,FILNAM,HEADER_TEXT,NETAB(1))
C
            WRITE (IOTMP2,99006) PNVEC0
            WRITE (IOTMP2,99004) PNMAX
            WRITE (IOTMP2,99005) NPN
C------------------------------------------------------------------- TOT
            FILNAM = DATSET(1:LDATSET)//'_TOT'//EXTENSION
            OPEN (UNIT=IOTMP3,FILE=FILNAM)
            WRITE (6,99011) SPEC_TYPE,IOTMP3,FILNAM
            CALL WRHEAD(IOTMP3,FILNAM,HEADER_TEXT,NETAB(1))
C
            WRITE (IOTMP3,99006) PNVEC0
            WRITE (IOTMP3,99004) PNMAX
            WRITE (IOTMP3,99005) NPN
C-----------------------------------------------------------------------
            IF ( IPRINT.GT.0 ) OPEN (IOTMP,FILE=DATSET(1:LDATSET)//
     &                               '_data'//EXTENSION,
     &                               STATUS='UNKNOWN')
C
            DO I = IPNBOT,IPNTOP
               IF ( IPRINT.GT.0 ) WRITE (IOTMP,'(6E15.5)') I*DPN,MD(I,1)
     &              ,MD(I,2),MD(I,1) + MD(I,2),MD(I,1) - MD(I,2)
C
C------------------------------------------------------------------- BND
               DO IQA = 1,NQ
                  DO IQB = 1,NQ
                     MD_1QQ(IQA,IQB) = MD_SQQ(I,1,IQA,IQB)
                     MD_2QQ(IQA,IQB) = MD_SQQ(I,2,IQA,IQB)
                  END DO
               END DO
               WRITE (IOTMP1,99008) I*DPN,
     &                              ((MD_1QQ(IQA,IQB)+MD_2QQ(IQA,IQB),
     &                              MD_1QQ(IQA,IQB)-MD_2QQ(IQA,IQB),
     &                              IQA=1,NQ),IQB=1,NQ)
C
C------------------------------------------------------------------- COR
               MD_1QQ(:,:) = 0D0
               MD_2QQ(:,:) = 0D0
               DO IQA = 1,NQ
                  MD_1QQ(IQA,IQA) = MDCOR_SQ(I,1,IQA)
                  MD_2QQ(IQA,IQA) = MDCOR_SQ(I,2,IQA)
               END DO
               WRITE (IOTMP2,99008) I*DPN,
     &                              ((MD_1QQ(IQA,IQB)+MD_2QQ(IQA,IQB),
     &                              MD_1QQ(IQA,IQB)-MD_2QQ(IQA,IQB),
     &                              IQA=1,NQ),IQB=1,NQ)
C
C------------------------------------------------------------------- TOT
               DO IQA = 1,NQ
                  DO IQB = 1,NQ
                     MD_1QQ(IQA,IQB) = MD_SQQ(I,1,IQA,IQB)
                     MD_2QQ(IQA,IQB) = MD_SQQ(I,2,IQA,IQB)
                  END DO
                  MD_1QQ(IQA,IQA) = MD_1QQ(IQA,IQA) + MDCOR_SQ(I,1,IQA)
                  MD_2QQ(IQA,IQA) = MD_2QQ(IQA,IQA) + MDCOR_SQ(I,2,IQA)
               END DO
               WRITE (IOTMP3,99008) I*DPN,
     &                              ((MD_1QQ(IQA,IQB)+MD_2QQ(IQA,IQB),
     &                              MD_1QQ(IQA,IQB)-MD_2QQ(IQA,IQB),
     &                              IQA=1,NQ),IQB=1,NQ)
C
C-----------------------------------------------------------------------
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT .AND. I.LE.MIN(3,NPN) ) THEN
                  WRITE (IFILBUILDBOT,99012) 'band IQA,IQB',
     &                   ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                   ((MD_SQQ(I,1,IQA,IQB),MD_SQQ(I,2,IQA,IQB),
     &                   IQA=1,NQ),IQB=1,NQ)
                  WRITE (IFILBUILDBOT,99012) 'core IQA    ',
     &                   ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                   (MDCOR_SQ(I,1,IQA),MDCOR_SQ(I,2,IQA),IQA=1,NQ)
               END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            END DO
C
C1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D-1D
C
         ELSE
C
C2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D
C
C------------- outputs four column version sum and difference (Magnetic)
C------------------------- if we do not print with full double precision
C------------------- then printing sum and difference is more preferable
C----------------------------------------------- than spin_up spind_dwn!
C-----------------------------------------------------------------------
C------------------------------------ distinguish  COMPTON:  CP2D & MD2D
C-------------------------------------------- and   POSANI:  MD2D
C
            IF ( TASK(1:7).EQ.'COMPTON' ) THEN
               IF ( .NOT.SPEC_DIM_EQ_2D_MD2D ) THEN
                  EXTENSION2D = '_mcp2d.dat'
               ELSE
C --------------- output a cut of mom-density at NPMAX (2D)
                  EXTENSION2D = '_mmd2d.dat'
               END IF
            ELSE IF ( TASK(1:6).EQ.'POSANI' ) THEN
               EXTENSION2D = '_pan2d.dat'
            END IF
C
C------------------------------------------------------------------- BND
            FILNAM = DATSET(1:LDATSET)//'_BND'//EXTENSION2D
            OPEN (UNIT=IOTMP1,FILE=FILNAM)
C
C------------------------------------------------------------------- COR
            FILNAM = DATSET(1:LDATSET)//'_COR'//EXTENSION2D
            OPEN (UNIT=IOTMP2,FILE=FILNAM)
C
C------------------------------------------------------------------- TOT
            FILNAM = DATSET(1:LDATSET)//'_TOT'//EXTENSION2D
            OPEN (UNIT=IOTMP3,FILE=FILNAM)
C
C
            DO IPP1 = IPP1BOT,IPP1TOP
               DO IPP2 = IPP2BOT,IPP2TOP
                  MD2DBND1 = MD2D(IPP1,IPP2,1)
                  MD2DBND2 = MD2D(IPP1,IPP2,2)
                  MD2DCOR1 = MDCOR2D(IPP1,IPP2,1)
                  MD2DCOR2 = MDCOR2D(IPP1,IPP2,2)
                  WRITE (IOTMP1,99007) PPVEC1L*IPP1,PPVEC2L*IPP2,
     &                                 MD2DBND1 + MD2DBND2,
     &                                 MD2DBND1 - MD2DBND2
                  WRITE (IOTMP2,99007) PPVEC1L*IPP1,PPVEC2L*IPP2,
     &                                 MD2DCOR1 + MD2DCOR2,
     &                                 MD2DCOR1 - MD2DCOR2
                  WRITE (IOTMP3,99007) PPVEC1L*IPP1,PPVEC2L*IPP2,
     &                                 (MD2DBND1+MD2DBND2)
     &                                 + (MD2DCOR1+MD2DCOR2),
     &                                 (MD2DBND1-MD2DBND2)
     &                                 + (MD2DCOR1-MD2DCOR2)
                  IF ( IPP2.EQ.IPP2TOP ) THEN
                     WRITE (IOTMP1,'(a)')
                     WRITE (IOTMP2,'(a)')
                     WRITE (IOTMP3,'(a)')
                  END IF
               END DO
            END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT .AND. I.LE.MIN(3,NPN) ) THEN
               WRITE (IFILBUILDBOT,99013) 'band',
     &                ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                ((MD2D(IPP1,IPP2,1),MD2D(IPP1,IPP2,2),
     &                IPP1=IPP1BOT,(IPP1BOT+2)),IPP2=IPP2BOT,(IPP2BOT+2)
     &                )
               WRITE (IFILBUILDBOT,99013) 'core',
     &                ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                ((MDCOR2D(IPP1,IPP2,1),MDCOR2D(IPP1,IPP2,2),
     &                IPP1=IPP1BOT,(IPP1BOT+2)),IPP2=IPP2BOT,(IPP2BOT+2)
     &                )
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         END IF
C
C2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D-2D
C
         CALL CPU_TIME(TIME2)
         WRITE (6,99010) TIME2 - TIME1
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C-----------------------------------------------------------------------
      DEALLOCATE (CCJJ,ALPH,DELT,TAUA,MD_SQQ,IPIV,MAUX,TAUK)
      DEALLOCATE (MD_1QQ,MD_2QQ,DMAMC,TAUAB,DMATT,DTILT,ME_PWZR,ME_ZLPW)
      DEALLOCATE (MSSQL,TAUQL,W1,TAUKAB,YVALPA,YVALPB,GAM,BET)
      DEALLOCATE (TAUKQQL,ZG,ZF,RME_PW_ZR_P)
C-----------------------------------------------------------------------
C
      CALL STOP_REGULAR(ROUTINE,' ')
C=======================================================================
C
99001 FORMAT (/,79('*'),/,24X,A,/,79('*'),//,10X,
     &        'Fermi energy            ',F10.6,' Ry',//,10X,
     &        'MEMODE     =',I4,10X,
     &        '1: interpolate M(p), 2: calculate M(p)',//,10X,
     &        'GF_CONV_RH =   ',L1,//,10X,'SLITSS     =   ',L1,//,10X,
     &        'USE_ZJOVZZ =   ',L1,//,10X,'NE         =',I4,//,10X,
     &        'p-mesh settings to start (default or input)',//,10X,
     &        'pp1:',13X,'V:',3F8.4,/,10X,'pp2:',13X,'V:',3F8.4,/,10X,
     &        'pn :',13X,'V:',3F8.4,//,10X,'p-mesh settings (adjusted)',
     &        //,10X,'pp1:  NPP1=',I4,'  V:',3F8.4,'  max:',F7.3,
     &        '  del:',F7.4,/,10X,'pp2:  NPP2=',I4,'  V:',3F8.4,
     &        '  max:',F7.3,'  del:',F7.4,/,10X,'pn :  NPN =',I4,'  V:',
     &        3F8.4,'  max:',F7.3,'  del:',F7.4,/)
99002 FORMAT (10X,'B-vector u',I3,3X,3F7.3,'  2*pi/a')
99003 FORMAT (/,10X,'2-dimensional momentum distribution   shift=',L3,/,
     &        10X,'pnvec   :',3F12.9,/,10X,'ppvec1L :',F12.9,/,10X,
     &        'ppvec2L :',F12.9)
99004 FORMAT (10x,'Pzmax(a.u.) ',F8.4)
99005 FORMAT (10x,'ngridpts.   ',I3)
99006 FORMAT (10x,'PNVEC       ',3F8.4)
99007 FORMAT (4E24.12)
99008 FORMAT (30E20.12)
99009 FORMAT (/,10X,'parameters for spline interpolation',//,10X,
     &        'max(r) =',F10.6,4X,'PFEMAX = ',F10.4,4X,'DELTA =',F10.6,
     &        /)
99010 FORMAT (/,10X,'CPU time used :',F10.4,' sec',/)
99011 FORMAT (10X,A,'   :  (',I2,') ',A)
99012 FORMAT ('# BUILDBOT: ',A,':  1D Compton ',A,/,(1PE22.14))
99013 FORMAT ('# BUILDBOT: ',A,':  2D Compton ',A,/,(1PE22.14))
      END
C*==compton_checkvec.f    processed by SPAG 6.70Rc at 15:06 on  6 Apr 2017
      SUBROUTINE COMPTON_CHECKVEC(A,B,C)
C    *******************************************************************
C    *                                                                 *
C    *   check whether vectors A and B are perpendicular to C          *
C    *                                                                 *
C    *                  -> A x -> B || -> C                            *
C    *                                                                 *
C    *                                                                 *
C    *******************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTON_CHECKVEC')
C
C Dummy arguments
C
      REAL*8 A(3),B(3),C(3)
C
C Local variables
C
      REAL*8 CXV,C_CROSS_V(3),C_DOT_V,V(3)
      REAL*8 DDOT,DNRM2
C
C*** End of declarations rewritten by SPAG
C
      CALL RVECXPRO(A,B,V)
C
      CALL RVECXPRO(C,V,C_CROSS_V)
C
      C_DOT_V = DDOT(3,C,1,V,1)
C
      CXV = DNRM2(3,C_CROSS_V,1)
C
      IF ( C_DOT_V.LT.0D0 .OR. CXV.GT.1D-8 ) THEN
C
         WRITE (6,99001) A,'a',DNRM2(3,A,1)
         WRITE (6,99001) B,'b',DNRM2(3,B,1)
         WRITE (6,99001) C,'c',DNRM2(3,C,1)
         WRITE (6,99001) V,'v',DNRM2(3,V,1)
         WRITE (6,99001) C_CROSS_V,'c x v',DNRM2(3,C_CROSS_V,1)
C
         CALL STOP_MESSAGE(ROUTINE,'A x -> B || -> C is FALSE')
C
      END IF
C
99001 FORMAT (3F12.8,A4,F12.8)
C
      END
C*==compton_sum.f    processed by SPAG 6.70Rc at 15:06 on  6 Apr 2017
      SUBROUTINE COMPTON_SUM(ME_PWZR,ITA,TAUAB,ME_ZLPW,ITB,N,M,CSUMS)
C     ******************************************************************
C     *                                                                *
C     *    evaluate the terms   M(A) * TAU(A,B) * M(B)                 *
C     *                                                                *
C     ******************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ITA,ITB,M,N
      COMPLEX*16 CSUMS(2),ME_PWZR(M,2,*),ME_ZLPW(M,2,*),TAUAB(M,M)
C
C Local variables
C
      COMPLEX*16 CSUM1,CSUM2
      INTEGER IKM1,IKM2,IS
C
C*** End of declarations rewritten by SPAG
C
      DO IS = 1,2
         CSUM1 = 0D0
         DO IKM1 = 1,N
            CSUM2 = 0D0
            DO IKM2 = 1,N
               CSUM2 = CSUM2 + TAUAB(IKM1,IKM2)*ME_ZLPW(IKM2,IS,ITB)
C
            END DO
C
            CSUM1 = CSUM1 + ME_PWZR(IKM1,IS,ITA)*CSUM2
         END DO
C
         CSUMS(IS) = CSUM1
      END DO
C
      END
