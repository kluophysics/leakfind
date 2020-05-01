C*==chicalc.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHICALC(SIGNAF,TREATAF,ILOP,RPWOP,MEZZ,MEZJ,BCP,BCPS,
     &                   TAUT,TKTKTT,DDTAUTAUT,CHILANDAU,CHILANLT,
     &                   RHORINTL,RHORTL1,EMM,ENM,ENN,GAMMAS,GAMMAN,
     &                   BXCNMM,BXCNNM,BXCNMN,BXCNNN,AXCN,ORBSQINT,DOS1,
     &                   DOS1L,DOSL1,CHIPRINT,ERYD,IE,IECURR,NTKMAX,
     &                   NTKTKMAX,IDOS,ISPN,IORB,IHFI,IHVV,NSPINOP,
     &                   NSPINOBS)
C   ********************************************************************
C   *                                                                  *
C   *   main routine for the calculation of magn. susceptibility       *
C   *   - call <CHIRADINT> to calculate the radial integrals           *
C   *   - perform the energy integration                               *
C   *   - print results                                                *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1        1 IDOS                           *
C   *  spin moment           b s_z    2 ISPN                           *
C   *  orbital moment        b l_z    3 IORB   ____ NOP                *
C   *  hyperfine interaction H_hf,z   4 IHFI                           *
C   *  VanVlecK hyperfine    H_VV,z   5 IHVV   ____ NOBS = NOP+2       *
C   *  diamagnetic suscept   r^2      6 ICDIA                          *
C   *  diamagnetic hyperfine r^-1     7 IKDIA                          *
C   *  hyperfine r-weight    r^-3     8 IRM3   ____ NOBSDNS = NOBS+3   *
C   *                                                                  *
C   *  MD, HF, HE 1996 - 2000                                          *
C   *  SM, HE  2003                                                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,WETAB,IGRID,NETAB
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_RMESH,ONLY:JRWS,NMMAX,NRMAX,R2DRDI,DRDI,R
      USE MOD_SITES,ONLY:IQAT
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,LINMAX,NLMAX,NMUEMAX,TXT_L,
     &    IKM2LIN,IKM1LIN,NLINQ,NLQ
      USE MOD_TYPES,ONLY:NT,NTMAX,LTXT_T,TXT_T,CONC,NAT,RHOCHRC,RHOCHR,
     &    Z,IMT,NLT
      USE MOD_FILES,ONLY:LDATSET,DATSET,IOTMP
      USE MOD_CONSTANTS,ONLY:PI,MB_CGS,C0,RY_ERG,CHI_AU2CGS
      IMPLICIT NONE
C*--CHICALC43
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NOP,NOBS,NOBSDNS,ICDIA,IKDIA,IRM3
      PARAMETER (NOP=3,NOBS=NOP+2,NOBSDNS=NOBS+3,ICDIA=6,IKDIA=7,IRM3=8)
      REAL*8 CU,TOLRADINT,KU
      PARAMETER (CU=CHI_AU2CGS*1D+6,TOLRADINT=1D-8,
     &           KU=100D0*MB_CGS/RY_ERG)
C
C Dummy arguments
C
      LOGICAL CHILANDAU,TREATAF
      INTEGER CHIPRINT,IDOS,IE,IECURR,IHFI,IHVV,ILOP,IORB,ISPN,NTKMAX,
     &        NTKTKMAX
      COMPLEX*16 ERYD
      REAL*8 AXCN(NRMAX,2,NLMAX,NTMAX),BCP(NTMAX),BCPS(NTMAX),
     &       BXCNMM(NRMAX,NTMAX),BXCNMN(NRMAX,NTMAX),BXCNNM(NRMAX,NTMAX)
     &       ,BXCNNN(NRMAX,NTMAX),CHILANLT(0:NLMAX,NTMAX,2),
     &       EMM(NRMAX,NTMAX),ENM(NRMAX,NTMAX),ENN(NRMAX,NTMAX),
     &       GAMMAN(NRMAX,NTMAX),GAMMAS(NRMAX,NTMAX),ORBSQINT(0:NLMAX,2)
     &       ,RHORINTL(NRMAX,NLMAX,NTMAX),RPWOP(NRMAX,NLMAX*2),
     &       SIGNAF(NTMAX)
      COMPLEX*16 DDTAUTAUT(NTKTKMAX,NTMAX),DOS1(NTMAX),
     &           DOS1L(NLMAX,NTMAX),DOSL1(NLMAX,NTMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           RHORTL1(NRMAX,NLMAX,NTMAX,2,NOBS),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TKTKTT(NTKTKMAX,NTMAX,NTMAX)
      INTEGER NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      REAL*8 AA(:,:),AAHF(:,:),BB(:,:),BBHFN(:),BBHFO(:),BBHFS(:),
     &       CHIAFSS(:),CHIAFSS0(:),CHIN(:),CHINN(:),CHINN0(:),
     &       CHINNL(:,:),CHINR(:,:),CHIO(:,:),CHIOO(:,:),CHIOO0(:,:),
     &       CHIOOL(:,:,:),CHIOOLTR(:,:,:,:),CHIOOTR(:,:,:),CHIOR(:,:,:)
     &       ,CHIOS(:,:),CHIOS0(:,:),CHIS(:),CHISO(:),CHISO0(:),
     &       CHISR(:,:),CHISS(:),CHISS0(:),CHISSL(:,:),CHISSLTR(:,:,:),
     &       CHISSRI,CHISSTR(:,:),CHISUM,DELTAEF,DOBSEF(:,:,:,:),DOSEF,
     &       GAMMAOP(:,:),GAMMAOPINT(2),ICHIAFSS0(:),ICHIO(:),ICHIS(:),
     &       MJ,MJMAX,MJMIN,NORM,POBS(:,:,:,:),POBSR(:,:,:,:,:),RAT,
     &       RHOLOP(:,:),RHOLOPS(:,:),RHORTLEF(:,:,:,:,:),RINT(:),
     &       RINTCDIA(:),RINTCDIAC(:),RINTGAMMAOP(:,:),RINTKDIA(:),
     &       RINTKDIAC(:),RINTX(:,:,:,:,:,:),RINTZ(:,:,:,:,:,:),RSQU,
     &       RSSZLOP(:,:),SFTO(:),SFTO0(:),SFTOL(:,:),SFTOL0(:,:),
     &       SFTOLTR(:,:,:),SFTOTR(:,:),SFTS(:),SFTS0(:),SFTSL(:,:),
     &       SFTSL0(:,:),SFTSLTR(:,:,:),SFTSTR(:,:),SUMEF,SUMWTZNN,
     &       SUMWTZNO(2),SUMWTZNS,T0X(:,:,:,:,:),T0XL(:,:,:,:,:,:),
     &       T0XR(:,:,:,:,:,:)
      LOGICAL COUPLI,INSIDE,NONMAG,NTERM
      COMPLEX*16 D0XL(:,:,:,:,:,:),D0XL1(:,:,:,:,:,:),D0XR(:,:,:,:,:,:),
     &           D0XR1(:,:,:,:,:,:),D0ZL(:,:,:,:,:,:),D0ZL1(:,:,:,:,:,:)
     &           ,D1ZL(:,:,:,:,:,:),D1ZL1(:,:,:,:,:,:),D2,DELE,
     &           DIJXL(:,:,:,:,:,:,:),DIJXL1(:,:,:,:,:,:,:),
     &           DIJXLTR(:,:,:,:,:,:,:),DIJXLTR1(:,:,:,:,:,:,:),
     &           DIJXR(:,:,:,:,:,:,:),DIJXR1(:,:,:,:,:,:,:),
     &           DIJZL(:,:,:,:,:,:,:),DIJZL1(:,:,:,:,:,:,:),
     &           DIJZLTR(:,:,:,:,:,:,:),DIJZLTR1(:,:,:,:,:,:,:),
     &           DOBS(:,:,:,:),DOBS1(:,:,:,:),DOBSE(:,:,:,:,:),DOS(:),
     &           DOSI(:),DOSL(:,:),DOSM(:),DZJ(:,:),DZL(:,:,:,:,:,:),
     &           DZL1(:,:,:,:,:,:),DZLR(:,:,:,:,:,:,:),
     &           DZLR1(:,:,:,:,:,:,:),DZZ(:,:),GGINT(:),JFA(:,:,:),
     &           JGA(:,:,:),RHORTL(:,:,:,:,:),SMT(:),SMTL(:,:),SMTM(:),
     &           SZJ(:,:),SZZ(:,:),TAUTLIN(:,:),W1,W2,ZFA(:,:,:),
     &           ZGA(:,:,:)
      CHARACTER*80 FILNAM
      INTEGER I,I1,I2,ICALL,IEF,IL,IM,INFO,IOBS,IOFF1,IOFF2,IOP,IPIV(:),
     &        IQ,IR,IRTOP,IS,IS1,IS2,ISKIP,ISOCOUPLING,IT,J,JT,JTP,K1,
     &        K2,KMAX,L,LFN,LIN,LMAX,MUE,MUETOP,NDIM,NL,NSPIN,NTK
      REAL*8 T0Z(:,:,:,:,:),T0ZL(:,:,:,:,:,:),T1X(:,:,:,:,:),
     &       T1Z(:,:,:,:,:),T1ZL(:,:,:,:,:,:),TCDIA(NTMAX),
     &       TCDIACRI(NTMAX),TCDIARI(NTMAX),TIJX(:,:,:,:,:,:),
     &       TIJXL(:,:,:,:,:,:,:),TIJXLTR(:,:,:,:,:,:,:),
     &       TIJXR(:,:,:,:,:,:,:),TIJXTR(:,:,:,:,:,:),TIJZ(:,:,:,:,:,:),
     &       TIJZL(:,:,:,:,:,:,:),TIJZLTR(:,:,:,:,:,:,:),
     &       TIJZTR(:,:,:,:,:,:),TKDIA(NTMAX),TKDIACRI(NTMAX),
     &       TKDIARI(NTMAX),TOBS(:,:,:,:),TOL,TRADD,TRM3(NTMAX),TRTEST,
     &       TSSZRLOP(:,:),TZ(:,:,:,:,:),TZL(:,:,:,:,:,:),
     &       TZLR(:,:,:,:,:,:,:),TZR(:,:,:,:,:,:),VH(:),VH_A(:),VH_B(:),
     &       WIHF,WIN,WIO(2),WIS,WJ,WJP,WT0XNN(NTMAX),WT0XNO(NTMAX,2),
     &       WT0XNS(NTMAX),WTIJXNN(NTMAX),WTIJXNO(NTMAX,2),
     &       WTIJXNS(NTMAX),WTXNN(NTMAX),WTXNO(NTMAX,2),WTXNS(NTMAX),
     &       WTZNN(NTMAX),WTZNO(NTMAX,2),WTZNS(NTMAX)
      CHARACTER*2 TXTS(2)
      SAVE D0XL,D0XL1,D0XR,D0XR1,D0ZL,D0ZL1,D1ZL,D1ZL1,DIJXL,DIJXL1,
     &     DIJXLTR,DIJXLTR1,DIJXR,DIJXR1,DIJZL,DIJZL1,DIJZLTR,DIJZLTR1,
     &     DOBS,DOBS1,DOBSE,DOS,DOSI,DOSL,DOSM,DZJ,DZL,DZL1,DZLR,DZLR1,
     &     DZZ,GGINT,JFA,JGA,RHORTL,RHORTLEF,SMT,SMTL,SMTM,SZJ,SZZ,T0XL,
     &     T0XR,T0Z,T0ZL,T1X,T1ZL,TAUTLIN,TIJX,TIJXL,TIJXLTR,TIJXR,
     &     TIJXTR,TIJZ,TIJZL,TIJZLTR,TIJZTR,TOBS,TSSZRLOP,TZ,TZL,TZLR,
     &     TZR,ZFA,ZGA
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/,TXTS/'dn','up'/
C
C
C
      ALLOCATABLE RINTCDIA,RINTKDIA,CHIOOLTR,CHISSLTR,ICHIAFSS0
      ALLOCATABLE RINTCDIAC,RINTKDIAC,AA,BB,VH,T0X,T1Z
      ALLOCATABLE RINTGAMMAOP,AAHF,CHIN,CHIO,CHIS,POBS,IPIV,RINT,SFTO
      ALLOCATABLE SFTS,SFTO0,SFTS0,BBHFN,BBHFO,BBHFS,ICHIO,CHINN
      ALLOCATABLE CHIOO,CHINR,CHIOR,CHIOS,CHISO,CHISR,CHISS,ICHIS
      ALLOCATABLE POBSR,SFTOL,SFTSL,RINTX,RINTZ,CHINN0,CHIOO0
      ALLOCATABLE CHIOS0,CHISO0,CHISS0,SFTOL0,SFTSL0,DOBSEF,CHINNL
      ALLOCATABLE CHIOOL,CHISSL,RHOLOP,SFTOTR,SFTSTR,CHIAFSS
      ALLOCATABLE GAMMAOP,CHIOOTR,CHISSTR,RHOLOPS,SFTOLTR,CHIAFSS0
      ALLOCATABLE SFTSLTR,RSSZLOP
      ALLOCATABLE DZL1,TZL,D0ZL1,D0XL1,T0ZL,T0XL,D1ZL1,T1ZL
      ALLOCATABLE T0Z,T1X,D0XL,D0XR,D0ZL,D1ZL,RHORTLEF
      ALLOCATABLE DIJXL,DIJXLTR,DIJXR,DIJZL,DIJZLTR
      ALLOCATABLE DIJZL1,DIJXL1,TIJZL,TIJXL,TZLR,DZLR1
      ALLOCATABLE TIJX,TIJXTR,TIJZ,TIJZTR,TSSZRLOP,TZ,TZR
      ALLOCATABLE TOBS,DOBS1,T0XR,D0XR1,TIJXR,DIJXR1
      ALLOCATABLE DOBS,DOS,DOSI,DOSL,DOSM,DZJ,DZL,DZLR,DZZ,GGINT
      ALLOCATABLE JFA,JGA,RHORTL,SMT,SMTL,SMTM,SZJ,SZZ
      ALLOCATABLE TAUTLIN,ZFA,ZGA,VH_A,VH_B
      ALLOCATABLE DIJZLTR1,DIJXLTR1,TIJZLTR,TIJXLTR,DOBSE
C
      ALLOCATE (VH_A(NRMAX),VH_B(NRMAX))
      ALLOCATE (RINTCDIA(NRMAX),RINTKDIA(NRMAX))
      ALLOCATE (CHIOOLTR(NLMAX,NTMAX,NTMAX,0:2))
      ALLOCATE (CHISSLTR(NLMAX,NTMAX,NTMAX),ICHIAFSS0(NTMAX))
      ALLOCATE (RINTCDIAC(NRMAX),RINTKDIAC(NRMAX))
      ALLOCATE (AA(NTMAX*4,NTMAX*4),BB(NTMAX*4,1),VH(NRMAX))
      ALLOCATE (T0X(NTMAX,2,NOBS,2,NOP),T1Z(NTMAX,2,NOBS,2,NOP))
      ALLOCATE (RINTGAMMAOP(NRMAX,2),AAHF(NTMAX*4,NTMAX*4))
      ALLOCATE (CHIN(NTMAX),CHIO(NTMAX,0:2),CHIS(NTMAX))
      ALLOCATE (POBS(0:NLMAX,NTMAX,2,NOBS),IPIV(4*NTMAX))
      ALLOCATE (RINT(NRMAX),SFTO(NTMAX),SFTS(NTMAX),SFTO0(NTMAX))
      ALLOCATE (SFTS0(NTMAX),BBHFN(4*NTMAX),BBHFO(NTMAX*4))
      ALLOCATE (BBHFS(NTMAX*4),ICHIO(NTMAX),CHINN(NTMAX))
      ALLOCATE (CHIOO(NTMAX,0:2),CHINR(NRMAX,NTMAX))
      ALLOCATE (CHIOR(NRMAX,NTMAX,0:2),CHIOS(NTMAX,0:2))
      ALLOCATE (CHISO(NTMAX),CHISR(NRMAX,NTMAX),CHISS(NTMAX))
      ALLOCATE (ICHIS(NTMAX),POBSR(NRMAX,0:NLMAX,NTMAX,2,NOBS))
      ALLOCATE (SFTOL(NLMAX,NTMAX),SFTSL(NLMAX,NTMAX))
      ALLOCATE (RINTX(NRMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (RINTZ(NRMAX,NTMAX,2,NOBS,2,NOP),CHINN0(NTMAX))
      ALLOCATE (CHIOO0(NTMAX,0:2),CHIOS0(NTMAX,0:2),CHISO0(NTMAX))
      ALLOCATE (CHISS0(NTMAX),SFTOL0(NLMAX,NTMAX))
      ALLOCATE (SFTSL0(NLMAX,NTMAX),DOBSEF(NLMAX,NTMAX,2,NOBSDNS))
      ALLOCATE (CHINNL(NLMAX,NTMAX),CHIOOL(NLMAX,NTMAX,0:2))
      ALLOCATE (CHISSL(NLMAX,NTMAX))
      ALLOCATE (RHOLOP(NRMAX,NTMAX),SFTOTR(NTMAX,NTMAX))
      ALLOCATE (SFTSTR(NTMAX,NTMAX),CHIAFSS(NTMAX),GAMMAOP(NRMAX,2))
      ALLOCATE (CHIOOTR(NTMAX,NTMAX,0:2),CHISSTR(NTMAX,NTMAX))
      ALLOCATE (RHOLOPS(NRMAX,2),SFTOLTR(NLMAX,NTMAX,NTMAX))
      ALLOCATE (CHIAFSS0(NTMAX),SFTSLTR(NLMAX,NTMAX,NTMAX))
      ALLOCATE (RSSZLOP(NRMAX,NTMAX))
C
C*** End of declarations rewritten by SPAG
C
      NTERM = .TRUE.
      NONMAG = .FALSE.
C
      ICALL = ICALL + 1
C
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( IORB.NE.NOP ) STOP 'in <CHICALC>  IORB <> NOP'
         IF ( IHVV.NE.NOBS ) STOP 'in <CHICALC>  IHVV <> NOBS'
         IF ( IRM3.NE.NOBSDNS ) STOP 'in <CHICALC>  IRM3 <> NOBSDNS'
C
         ALLOCATE (T1X(NTMAX,2,NOBS,2,NOP),T0Z(NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DZL1(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DZLR1(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TZLR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D1ZL1(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (T1ZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0ZL1(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0XR1(NRMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0XL1(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0XL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0XR(NRMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D0ZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (D1ZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (T0ZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (T0XL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (RHORTLEF(NRMAX,NLMAX,NTMAX,2,NOBS))
         ALLOCATE (T0XR(NRMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJZL1(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXR(NRMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJZL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJZLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXL1(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXR1(NRMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJZL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJXL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJXR(NRMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJZLTR1(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DIJXLTR1(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJZLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJXLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJX(NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJXTR(NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJZ(NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TIJZTR(NTMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TSSZRLOP(NRMAX,NTMAX))
         ALLOCATE (TZ(NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TZR(NRMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (TOBS(NLMAX,NTMAX,2,NOBSDNS))
         ALLOCATE (DOBS1(NLMAX,NTMAX,2,NOBSDNS))
         ALLOCATE (DOBSE(NLMAX,NTMAX,2,NOBSDNS,NEMAX))
         ALLOCATE (DOBS(NLMAX,NTMAX,2,NOBSDNS))
         ALLOCATE (DOS(NTMAX))
         ALLOCATE (DOSI(NTMAX))
         ALLOCATE (DOSL(NLMAX,NTMAX))
         ALLOCATE (DOSM(NMUEMAX))
         ALLOCATE (DZJ(LINMAX,NTMAX))
         ALLOCATE (DZL(NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DZLR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
         ALLOCATE (DZZ(LINMAX,NTMAX))
         ALLOCATE (GGINT(NTMAX))
         ALLOCATE (JFA(NRMAX,2,2))
         ALLOCATE (JGA(NRMAX,2,2))
         ALLOCATE (RHORTL(NRMAX,NLMAX,NTMAX,2,NOBS))
         ALLOCATE (SMT(NTMAX))
         ALLOCATE (SMTL(NLMAX,NTMAX))
         ALLOCATE (SMTM(NMUEMAX))
         ALLOCATE (SZJ(LINMAX,NTMAX))
         ALLOCATE (SZZ(LINMAX,NTMAX))
         ALLOCATE (TAUTLIN(NTKMAX,NTMAX))
         ALLOCATE (ZFA(NRMAX,2,2))
         ALLOCATE (ZGA(NRMAX,2,2))
C
      END IF
C
      IF ( IE.EQ.1 ) THEN
C
         CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,DZL1)
         CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D1ZL1)
         CALL CINIT(NRMAX*NTMAX*2*NOBS*2*NOP,D0XR1)
         CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D0ZL1)
         CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D0XL1)
         CALL CINIT(NLMAX*NTMAX*2*NOBSDNS,DOBS1)
C
         CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJZL1)
         CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXL1)
         CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJZLTR1)
         CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXLTR1)
         CALL CINIT(NRMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXR1)
         CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,DZLR1)
         CALL RINIT(NTMAX*2*NOBS*2*NOP,TZ)
         CALL RINIT(NLMAX*NTMAX*2*NOBS*2*NOP,TZL)
         CALL RINIT(NLMAX*NTMAX*2*NOBS*2*NOP,T1ZL)
         CALL RINIT(NLMAX*NTMAX*2*NOBS*2*NOP,T0ZL)
         CALL RINIT(NLMAX*NTMAX*2*NOBS*2*NOP,T0XL)
         CALL RINIT(NRMAX*NTMAX*2*NOBS*2*NOP,T0XR)
         CALL RINIT(NLMAX*NTMAX*2*NOBSDNS,TOBS)
         CALL RINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,TIJZL)
         CALL RINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,TIJXL)
         CALL RINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,TIJZLTR)
         CALL RINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,TIJXLTR)
         CALL RINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,TZLR)
         CALL RINIT(NRMAX*NTMAX*NTMAX*2*NOBS*2*NOP,TIJXR)
C
      END IF
      CALL RINIT(NRMAX*(NLMAX+1)*NTMAX*2*NOBS,POBSR)
C
C------------------------------------------------- recover old variables
C
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         IQ = IQAT(1,IT)
         DO LIN = 1,NLINQ(IQ)
            I1 = IKM1LIN(LIN)
            I2 = IKM2LIN(LIN)
            DZZ(LIN,IT) = MEZZ(I1,I2,IT,1)
            DZJ(LIN,IT) = MEZJ(I1,I2,IT,1)
            SZZ(LIN,IT) = MEZZ(I1,I2,IT,2)
            SZJ(LIN,IT) = MEZJ(I1,I2,IT,2)
            TAUTLIN(LIN,IT) = TAUT(I1,I2,IT)
         END DO
      END DO
C
C=======================================================================
      IF ( IECURR.EQ.1 ) THEN
C
         DO IT = 1,NT
            DOSI(IT) = C0
            GGINT(IT) = C0
            DO IL = 1,NLMAX
               DOSL1(IL,IT) = C0
            END DO
         END DO
         CALL RINIT(NRMAX*NLMAX*NTMAX,RHORINTL)
         CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS,RHORTL1)
      END IF
C=======================================================================
C
      NTK = NLINQ(1)
      CALL CHIRADINT(ZGA,ZFA,JGA,JFA,NOBS,NOBSDNS,NOP,IDOS,ISPN,IORB,
     &               IHFI,IHVV,ICDIA,IKDIA,IRM3,NSPINOBS,NSPINOP,DZL,
     &               DZLR,D0ZL,D0XL,D0XR,D1ZL,DIJZL,DIJXL,DIJZLTR,
     &               DIJXLTR,DIJXR,DOBS,RHORTL,TAUTLIN,TKTKTT,DDTAUTAUT,
     &               BXCNMM,BXCNNM,BXCNMN,BXCNNN,AXCN,CHIPRINT,NTK,
     &               NTKMAX,NTKTKMAX)
C
      CALL RINIT(NTMAX*2*NOBS*2*NOP,TZ)
      CALL RINIT(NRMAX*NTMAX*2*NOBS*2*NOP,TZR)
      CALL RINIT(NTMAX*2*NOBS*2*NOP,T1Z)
      CALL RINIT(NTMAX*2*NOBS*2*NOP,T1X)
      CALL RINIT(NTMAX*2*NOBS*2*NOP,T0Z)
      CALL RINIT(NTMAX*2*NOBS*2*NOP,T0X)
      CALL RINIT(NTMAX*NTMAX*2*NOBS*2*NOP,TIJZ)
      CALL RINIT(NTMAX*NTMAX*2*NOBS*2*NOP,TIJX)
      CALL RINIT(NTMAX*NTMAX*2*NOBS*2*NOP,TIJZTR)
      CALL RINIT(NTMAX*NTMAX*2*NOBS*2*NOP,TIJXTR)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         NL = NLQ(IQ)
         LMAX = NL - 1
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         LIN = 0
         DOS(IT) = C0
         SMT(IT) = C0
         DOSI(IT) = C0
         GGINT(IT) = C0
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,LMAX
            IL = L + 1
C
            DOSL(IL,IT) = C0
            SMTL(IL,IT) = C0
C
            IF ( IREL.GT.1 ) THEN
               MJMAX = DBLE(L) + 0.5D0
               KMAX = 2
               MUETOP = 2*L + 2
            ELSE
               MJMAX = DBLE(L)
               KMAX = 1
               MUETOP = 2*L + 1
            END IF
            MJMIN = -MJMAX
            MJ = MJMIN - 1D0
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MUE = 1,MUETOP
               MJ = MJ + 1D0
C
               DOSM(MUE) = C0
               SMTM(MUE) = C0
C
               IF ( IREL.LE.1 ) THEN
                  COUPLI = .FALSE.
C           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
               ELSE IF ( ABS(MJ).GT.DBLE(L) ) THEN
                  COUPLI = .FALSE.
               ELSE
                  COUPLI = .TRUE.
               END IF
C
               DO K1 = 1,KMAX
                  DO K2 = 1,KMAX
                     LIN = LIN + 1
C
                     DOSM(MUE) = DOSM(MUE) - TAUTLIN(LIN,IT)*DZZ(LIN,IT)
     &                           /PI
                     SMTM(MUE) = SMTM(MUE) - TAUTLIN(LIN,IT)*SZZ(LIN,IT)
     &                           /PI
C
C    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
                     IF ( K1.EQ.K2 ) THEN
                        DOSM(MUE) = DOSM(MUE) + DZJ(LIN,IT)/PI
                        SMTM(MUE) = SMTM(MUE) + SZJ(LIN,IT)/PI
                        IF ( .NOT.COUPLI ) GOTO 10
                     END IF
                  END DO
               END DO
C
 10            CONTINUE
               DOSL(IL,IT) = DOSL(IL,IT) + DOSM(MUE)
               SMTL(IL,IT) = SMTL(IL,IT) + SMTM(MUE)
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
            DOS(IT) = DOS(IT) + DOSL(IL,IT)
            SMT(IT) = SMT(IT) + SMTL(IL,IT)
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         IF ( IECURR.EQ.1 ) THEN
            DO IL = 1,NL
               DOS1L(IL,IT) = DOSL(IL,IT)
            END DO
            DOS1(IT) = DOS(IT)
         END IF
C
C ======================================================================
C               perform the energy - integration for all variables
C ======================================================================
C
C -------------- set the weight factors W1 and W2 to deal with any IGRID
         IF ( IGRID(1).EQ.5 ) THEN
            W1 = 0D0
            W2 = WETAB(IE,1)
         ELSE IF ( IECURR.EQ.1 ) THEN
            W1 = 0D0
            W2 = W1
         ELSE
            DELE = ETAB(IE,1) - ETAB(IE-1,1)
            W1 = 0.5D0*DELE
            W2 = W1
         END IF
C
         DO IOP = 1,NOP
            DO IS2 = 1,NSPINOP(IOP)
               DO IOBS = 1,NOBS
                  DO IS1 = 1,NSPINOBS(IOBS)
                     DO IL = 1,NL
C
                        D2 = DZL(IL,IT,IS1,IOBS,IS2,IOP)/PI
                        TZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = TZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + DIMAG(W2*D2+W1*DZL1(IL,IT,IS1,IOBS,IS2,IOP)
     &                     )
                        DZL1(IL,IT,IS1,IOBS,IS2,IOP) = D2
C
                        DO I = 1,IRTOP
                           D2 = DZLR(I,IL,IT,IS1,IOBS,IS2,IOP)/PI
                           TZLR(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                        = TZLR(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG(W2*D2+W1*DZLR1(I,IL,IT,IS1,IOBS,
     &                        IS2,IOP))
                           DZLR1(I,IL,IT,IS1,IOBS,IS2,IOP) = D2
C
                        END DO
C
                        D2 = D0ZL(IL,IT,IS1,IOBS,IS2,IOP)/PI
                        T0ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = T0ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + DIMAG(W2*D2+W1*D0ZL1(IL,IT,IS1,IOBS,IS2,
     &                     IOP))
                        D0ZL1(IL,IT,IS1,IOBS,IS2,IOP) = D2
C
                        D2 = D1ZL(IL,IT,IS1,IOBS,IS2,IOP)/PI
                        T1ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = T1ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + DIMAG(W2*D2+W1*D1ZL1(IL,IT,IS1,IOBS,IS2,
     &                     IOP))
                        D1ZL1(IL,IT,IS1,IOBS,IS2,IOP) = D2
C
                        D2 = D0XL(IL,IT,IS1,IOBS,IS2,IOP)/PI
                        T0XL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = T0XL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + DIMAG(W2*D2+W1*D0XL1(IL,IT,IS1,IOBS,IS2,
     &                     IOP))
                        D0XL1(IL,IT,IS1,IOBS,IS2,IOP) = D2
C
                     END DO
C
                     DO I = 1,IRTOP
                        D2 = D0XR(I,IT,IS1,IOBS,IS2,IOP)/PI
                        T0XR(I,IT,IS1,IOBS,IS2,IOP)
     &                     = T0XR(I,IT,IS1,IOBS,IS2,IOP)
     &                     + DIMAG(W2*D2+W1*D0XR1(I,IT,IS1,IOBS,IS2,IOP)
     &                     )
                        D0XR1(I,IT,IS1,IOBS,IS2,IOP) = D2
                     END DO
C
                     DO JT = 1,NT
                        DO IL = 1,NL
C
                           D2 = DIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)/PI
                           TIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        = TIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG(W2*D2+W1*DIJZL1(IL,IT,JT,IS1,IOBS,
     &                        IS2,IOP))
                           DIJZL1(IL,IT,JT,IS1,IOBS,IS2,IOP) = D2
C
                           D2 = DIJXL(IL,IT,JT,IS1,IOBS,IS2,IOP)/PI
                           TIJXL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        = TIJXL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG(W2*D2+W1*DIJXL1(IL,IT,JT,IS1,IOBS,
     &                        IS2,IOP))
                           DIJXL1(IL,IT,JT,IS1,IOBS,IS2,IOP) = D2
C
                           D2 = DIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)/PI
                           TIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        = TIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG
     &                        (W2*D2+W1*DIJZLTR1(IL,IT,JT,IS1,IOBS,IS2,
     &                        IOP))
                           DIJZLTR1(IL,IT,JT,IS1,IOBS,IS2,IOP) = D2
C
                           D2 = DIJXLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)/PI
                           TIJXLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        = TIJXLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG
     &                        (W2*D2+W1*DIJXLTR1(IL,IT,JT,IS1,IOBS,IS2,
     &                        IOP))
                           DIJXLTR1(IL,IT,JT,IS1,IOBS,IS2,IOP) = D2
C
                        END DO
C
                        DO I = 1,IRTOP
                           D2 = DIJXR(I,IT,JT,IS1,IOBS,IS2,IOP)/PI
                           TIJXR(I,IT,JT,IS1,IOBS,IS2,IOP)
     &                        = TIJXR(I,IT,JT,IS1,IOBS,IS2,IOP)
     &                        + DIMAG(W2*D2+W1*DIJXR1(I,IT,JT,IS1,IOBS,
     &                        IS2,IOP))
                           DIJXR1(I,IT,JT,IS1,IOBS,IS2,IOP) = D2
                        END DO
                     END DO
C
                  END DO
               END DO
            END DO
         END DO
C
         DO IL = 1,NL
            DO IOBS = 1,NOBSDNS
               IF ( IOBS.EQ.IORB ) THEN
                  NSPIN = 2
               ELSE
                  NSPIN = 1
               END IF
               DO IS1 = 1,NSPIN
C
                  D2 = DOBS(IL,IT,IS1,IOBS)
                  TOBS(IL,IT,IS1,IOBS) = TOBS(IL,IT,IS1,IOBS)
     &               + DIMAG(W2*D2+W1*DOBS1(IL,IT,IS1,IOBS))
                  DOBS1(IL,IT,IS1,IOBS) = D2
                  DOBSE(IL,IT,IS1,IOBS,IECURR) = D2
C
               END DO
            END DO
C
            IF ( ABS(DOBS1(IL,IT,1,IDOS)-DOSL(IL,IT)).GT.1D-8 )
     &           WRITE (6,99001) 'L-DOS',TXT_L(IL-1),IT,TXT_T(IT),
     &                           DOBS1(IL,IT,1,IDOS),DOSL(IL,IT),
     &                           (1.D0-DOBS1(IL,IT,1,IDOS)/DOSL(IL,IT))
C
            GGINT(IT) = GGINT(IT)
     &                  + DCMPLX(0.D0,TZL(IL,IT,1,IDOS,1,IDOS))
            DOSI(IT) = DOSI(IT) + DCMPLX(0.D0,TOBS(IL,IT,1,IDOS))
C
            DO I = 1,IRTOP
               D2 = RHORTL(I,IL,IT,1,ISPN)
               RHORINTL(I,IL,IT) = RHORINTL(I,IL,IT)
     &                             + DIMAG(W2*D2+W1*RHORTL1(I,IL,IT,1,
     &                             ISPN))
               RHORTL1(I,IL,IT,1,ISPN) = D2
            END DO
C
         END DO
C
C=======================================================================
C              test of the sum rule     d G / dE = - G * G
C=======================================================================
C
         IF ( (CHIPRINT.GE.1) .OR. (IE.EQ.NETAB(1)) ) THEN
C
            IF ( CHIPRINT.GE.1 ) THEN
               CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_GG.dat',7,
     &                       FILNAM,LFN,2,IOTMP,'G*G   file',10,NTMAX)
               IF ( IE.EQ.1 ) THEN
                  WRITE (IOTMP,99004) IT,TXT_T(IT)(1:LTXT_T(IT)),NL
                  REWIND IOTMP
               END IF
               DO ISKIP = 1,4 + IE - 1
                  READ (IOTMP,*)
               END DO
C
               WRITE (IOTMP,'(f8.4,15(1x,f8.4))') DREAL(ERYD),
     &                DIMAG(DOS(IT)),DIMAG(GGINT(IT)+DOS1(IT)),
     &                (DIMAG(DOBS(IL,IT,1,IDOS)),
     &                (DCMPLX(0.0D0,TZL(IL,IT,1,IDOS,1,IDOS))
     &                +DOS1L(IL,IT)),IL=1,NL)
               CLOSE (IOTMP)
C
            END IF
C
            WRITE (6,'(/,1X,79(''-''))')
            WRITE (6,99016) IE,ERYD,IT,TXT_T(IT),'CRYSTAL TERMS       ',
     &                      DOS(IT),GGINT(IT) + DOS1(IT),DOSI(IT),
     &                      (IL-1,DOBS(IL,IT,1,IDOS),
     &                      DCMPLX(0.0D0,TZL(IL,IT,1,IDOS,1,IDOS))
     &                      +DOS1L(IL,IT),
     &                      DCMPLX(0.0D0,TOBS(IL,IT,1,IDOS)),IL=1,NL)
            WRITE (6,'(1X,79(''-''))')
         END IF
C
C=======================================================================
C    sum over angular momentum and store energy integrated variables
C=======================================================================
C
         IF ( IE.EQ.NETAB(1) ) THEN
C
            DO IOP = 1,NOP
               DO IS2 = 1,NSPINOP(IOP)
                  DO IOBS = 1,NOBS
                     DO IS1 = 1,NSPINOBS(IOBS)
                        DO IL = 1,NL
                           TZ(IT,IS1,IOBS,IS2,IOP)
     &                        = TZ(IT,IS1,IOBS,IS2,IOP)
     &                        + TZL(IL,IT,IS1,IOBS,IS2,IOP)
C
                           T0Z(IT,IS1,IOBS,IS2,IOP)
     &                        = T0Z(IT,IS1,IOBS,IS2,IOP)
     &                        + T0ZL(IL,IT,IS1,IOBS,IS2,IOP)
C
                           T0X(IT,IS1,IOBS,IS2,IOP)
     &                        = T0X(IT,IS1,IOBS,IS2,IOP)
     &                        + T0XL(IL,IT,IS1,IOBS,IS2,IOP)
C
                           T1Z(IT,IS1,IOBS,IS2,IOP)
     &                        = T1Z(IT,IS1,IOBS,IS2,IOP)
     &                        + T1ZL(IL,IT,IS1,IOBS,IS2,IOP)
C
                           DO JT = 1,NT
                              TIJZ(IT,JT,IS1,IOBS,IS2,IOP)
     &                           = TIJZ(IT,JT,IS1,IOBS,IS2,IOP)
     &                           + TIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)
C
                              TIJX(IT,JT,IS1,IOBS,IS2,IOP)
     &                           = TIJX(IT,JT,IS1,IOBS,IS2,IOP)
     &                           + TIJXL(IL,IT,JT,IS1,IOBS,IS2,IOP)
C
                              TIJZTR(IT,JT,IS1,IOBS,IS2,IOP)
     &                           = TIJZTR(IT,JT,IS1,IOBS,IS2,IOP)
     &                           + TIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
C
                              TIJXTR(IT,JT,IS1,IOBS,IS2,IOP)
     &                           = TIJXTR(IT,JT,IS1,IOBS,IS2,IOP)
     &                           + TIJXLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
                           END DO
C
                           DO I = 1,IRTOP
                              TZR(I,IT,IS1,IOBS,IS2,IOP)
     &                           = TZR(I,IT,IS1,IOBS,IS2,IOP)
     &                           + TZLR(I,IL,IT,IS1,IOBS,IS2,IOP)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            DO IL = 1,NL
               DO IOBS = 1,NOBS
                  DO IS1 = 1,NSPINOBS(IOBS)
C
                     DO I = 1,IRTOP
                        RHORTLEF(I,IL,IT,IS1,IOBS)
     &                     = DIMAG(RHORTL(I,IL,IT,IS1,IOBS))
                     END DO
                  END DO
               END DO
            END DO
C
            DO IL = 1,NL
               CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,
     &                       '_'//TXT_L(IL-1)//'_I-spin.dat',13,FILNAM,
     &                       LFN,2,IOTMP,'I-spn file',10,NTMAX)
               WRITE (IOTMP,99005) IT,TXT_T(IT)(1:LTXT_T(IT)),IL,
     &                             TXT_L(IL-1)
               WRITE (IOTMP,'(1E16.7,2X,6E16.7)')
     &                (R(I,IM),RHORINTL(I,IL,IT),
     &                TZLR(I,IL,IT,1,ISPN,1,ISPN),
     &                TZLR(I,IL,IT,1,IORB,1,IORB)
     &                +TZLR(I,IL,IT,1,IORB,2,IORB)
     &                +TZLR(I,IL,IT,2,IORB,1,IORB)
     &                +TZLR(I,IL,IT,2,IORB,2,IORB),
     &                TZLR(I,IL,IT,1,ISPN,1,ISPN),
     &                TZLR(I,IL,IT,1,IHFI,1,ISPN),I=1,IRTOP)
               CLOSE (IOTMP)
            END DO
C
C=======================================================================
C                 Test for integration consistency
C     !!!       (was created only for the values with non-precize
C     !!!                calculation of 4-th term)
C=======================================================================
C
            DO I = 1,IRTOP
               RINTZ(I,IT,1,ISPN,1,ISPN) = TZR(I,IT,1,ISPN,1,ISPN)
     &            *R2DRDI(I,IM)
               RINTX(I,IT,1,ISPN,1,ISPN) = T0XR(I,IT,1,ISPN,1,ISPN)
     &            *R2DRDI(I,IM)
               RINTZ(I,IT,1,IHFI,1,ISPN) = TZR(I,IT,1,IHFI,1,ISPN)
     &            *DRDI(I,IM)
               RINTX(I,IT,1,IHFI,1,ISPN) = T0XR(I,IT,1,IHFI,1,ISPN)
     &            *DRDI(I,IM)
C
               DO IS = 1,2
                  DO IS2 = 1,2
                     RINTZ(I,IT,IS,IORB,IS2,IORB)
     &                  = TZR(I,IT,IS,IORB,IS2,IORB)*R2DRDI(I,IM)
                     RINTX(I,IT,IS,IORB,IS2,IORB)
     &                  = T0XR(I,IT,IS,IORB,IS2,IORB)*R2DRDI(I,IM)
                  END DO
C
               END DO
            END DO
C
            CALL RRADINT(IM,RINTZ(1,IT,1,ISPN,1,ISPN),TRTEST)
            IF ( ABS(1D0-TRTEST/TZ(IT,1,ISPN,1,ISPN)).GT.TOLRADINT )
     &           WRITE (6,99003) 'TSSZ SS',IT,
     &                           TRTEST/TZ(IT,1,ISPN,1,ISPN)
C
C
            CALL RRADINT(IM,RINTX(1,IT,1,ISPN,1,ISPN),TRTEST)
            IF ( ABS(1D0-TRTEST/T0X(IT,1,ISPN,1,ISPN)).GT.TOLRADINT )
     &           WRITE (6,99003) 'TSSX SS',IT,
     &                           TRTEST/T0X(IT,1,ISPN,1,ISPN)
C
C
            CALL RRADINT_HFF(IM,IRTOP,RINTZ(1,IT,1,IHFI,1,ISPN),TRTEST)
C
            IF ( ABS(1D0-TRTEST/TZ(IT,1,IHFI,1,ISPN)).GT.TOLRADINT )
     &           WRITE (6,99003) 'TSSZ HF-S',IT,
     &                           TRTEST/TZ(IT,1,IHFI,1,ISPN)
C
            CALL RRADINT_HFF(IM,IRTOP,RINTX(1,IT,1,IHFI,1,ISPN),TRTEST)
C
            IF ( ABS(1D0-TRTEST/T0X(IT,1,IHFI,1,ISPN)).GT.TOLRADINT )
     &           WRITE (6,99003) 'TSSX HF-S',IT,
     &                           TRTEST/T0X(IT,1,IHFI,1,ISPN)
C
C
            TRTEST = 0.D0
            DO IS = 1,2
               DO IS2 = 1,2
                  CALL RRADINT(IM,RINTZ(1,IT,IS,IORB,IS2,IORB),TRADD)
                  TRTEST = TRTEST + TRADD
               END DO
            END DO
            IF ( ABS(1D0-TRTEST/(TZ(IT,1,IORB,1,IORB)+TZ(IT,1,IORB,2,
     &           IORB)+TZ(IT,2,IORB,1,IORB)+TZ(IT,2,IORB,2,IORB)))
     &           .GT.TOLRADINT ) WRITE (6,99003) 'TLLZ OO',IT,
     &                                  TRTEST/(TZ(IT,1,IORB,1,IORB)
     &                                  +TZ(IT,1,IORB,2,IORB)
     &                                  +TZ(IT,2,IORB,1,IORB)
     &                                  +TZ(IT,2,IORB,2,IORB))
C
            TRTEST = 0.0D0
            DO IS = 1,2
               DO IS2 = 1,2
                  CALL RRADINT(IM,RINTX(1,IT,IS,IORB,IS2,IORB),TRADD)
                  TRTEST = TRTEST + TRADD
               END DO
            END DO
            IF ( ABS(1D0-TRTEST/(T0X(IT,1,IORB,1,IORB)+T0X(IT,1,IORB,2,
     &           IORB)+T0X(IT,2,IORB,1,IORB)+T0X(IT,2,IORB,2,IORB)))
     &           .GT.TOLRADINT ) WRITE (6,99003) 'TLLX OO',IT,
     &                                  TRTEST/(T0X(IT,1,IORB,1,IORB)
     &                                  +T0X(IT,1,IORB,2,IORB)
     &                                  +T0X(IT,2,IORB,1,IORB)
     &                                  +T0X(IT,2,IORB,2,IORB))
C
C=======================================================================
C                 calculate diamagnetic contributions
C=======================================================================
C
            DO I = 1,IRTOP
               RSQU = R(I,IM)*R(I,IM)
               RINTCDIA(I) = -RHOCHR(I,IT)*R2DRDI(I,IM)*RSQU
               RINTKDIA(I) = -RHOCHR(I,IT)*DRDI(I,IM)*R(I,IM)
               RINTCDIAC(I) = -RHOCHRC(I,IT)*R2DRDI(I,IM)*RSQU
               RINTKDIAC(I) = -RHOCHRC(I,IT)*DRDI(I,IM)*R(I,IM)
            END DO
C
            CALL RRADINT(IM,RINTCDIA,TCDIARI(IT))
            CALL RRADINT(IM,RINTKDIA,TKDIARI(IT))
            CALL RRADINT(IM,RINTCDIAC,TCDIACRI(IT))
            CALL RRADINT(IM,RINTKDIAC,TKDIACRI(IT))
C
            TCDIA(IT) = 0.0D0
            TKDIA(IT) = 0.0D0
            TRM3(IT) = 0.0D0
C
            DO IL = 1,NL
               TCDIA(IT) = TCDIA(IT) + TOBS(IL,IT,1,ICDIA)
               TKDIA(IT) = TKDIA(IT) + TOBS(IL,IT,1,IKDIA)
               TRM3(IT) = TRM3(IT) + TOBS(IL,IT,1,IRM3)
            END DO
C
            RAT = TCDIA(IT)/(TCDIARI(IT)-TCDIACRI(IT))
            IF ( ABS(1D0-RAT).GT.TOLRADINT ) WRITE (6,99003) 'TCDIA',IT,
     &           RAT
            RAT = TKDIA(IT)/(TKDIARI(IT)-TKDIACRI(IT))
            IF ( ABS(1D0-RAT).GT.TOLRADINT ) WRITE (6,99003) 'TKDIA',IT,
     &           RAT
C
C=======================================================================
C               variables connected with the OP-term
C=======================================================================
C
            CALL RINIT(IRTOP,RHOLOP(1,IT))
            CALL RINIT(IRTOP,TSSZRLOP(1,IT))
C
            DO IL = 1,NL
               DO I = 1,IRTOP
C ------------------------ take care of numerical inaccuracy for small r
                  IF ( RHORINTL(I,IL,IT).LT.0.0D0 ) RHORINTL(I,IL,IT)
     &                 = 0.0D0
               END DO
            END DO
C
            IF ( ILOP.GT.0 ) THEN
               DO I = 1,IRTOP
                  RHOLOP(I,IT) = RHORINTL(I,ILOP,IT)
                  TSSZRLOP(I,IT) = TZLR(I,ILOP,IT,1,ISPN,1,ISPN)
                  RINT(I) = RHOLOP(I,IT)*R2DRDI(I,IM)
               END DO
C
               CALL RRADINT(IM,RINT,NORM)
C
               DO I = 1,IRTOP
                  RHOLOP(I,IT) = RHOLOP(I,IT)/NORM
               END DO
C
            END IF
C
         END IF
C === IE .EQ. NETAB ====================================================
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
      IF ( IGRID(1).EQ.3 ) THEN
         WRITE (10,99018) ETAB(IE,1),
     &                    (((0.5D0*DIMAG(DOBS(IL,IT,1,IDOS)+SMTL(IL,IT))
     &                    ),IL=1,NL),
     &                    ((0.5D0*DIMAG(DOBS(IL,IT,1,IDOS)-SMTL(IL,IT)))
     &                    ,IL=1,NL),IT=1,NT)
C
         WRITE (20,99018) ETAB(IE,1),
     &                    (((0.5D0*DIMAG(DOS1L(IL,IT)+TZL(IL,IT,1,IDOS,
     &                    1,IDOS)+SMTL(IL,IT))),IL=1,NL),
     &                    ((0.5D0*DIMAG(DOS1L(IL,IT)
     &                    +TZL(IL,IT,1,IDOS,1,IDOS)-SMTL(IL,IT))),IL=1,
     &                    NL),IT=1,NT)
      END IF
C
C
      IF ( IE.LT.NETAB(1) ) RETURN
C
C
C ======================================================================
C              terms connected with high-field susceptibility
C ======================================================================
C
C
C--------------------------------------------- use data for Fermi energy
C
      IEF = NETAB(1)
C
      CALL RINIT(NTMAX*4*NTMAX*4,AAHF)
C
      DO IT = 1,NT
         DO IL = 1,NLT(IT)
            DO IOBS = 1,NOBS
               DO IS1 = 1,NSPINOBS(IOBS)
                  DOBSEF(IL,IT,IS1,IOBS)
     &               = DIMAG(DOBSE(IL,IT,IS1,IOBS,IEF))
               END DO
            END DO
         END DO
      END DO
C
      DOSEF = 0D0
      DO IT = 1,NT
         DO IL = 1,NLT(IT)
            DOSEF = DOSEF + NAT(IT)*CONC(IT)*DOBSEF(IL,IT,1,IDOS)
         END DO
      END DO
C
      DO IOBS = 1,NOBS
         DO IS1 = 1,NSPINOBS(IOBS)
            DO IT = 1,NT
               POBS(0,IT,IS1,IOBS) = 0D0
               DO I = 1,IRTOP
                  POBSR(I,0,IT,IS1,IOBS) = 0D0
               END DO
               DO IL = 1,NLT(IT)
                  POBS(IL,IT,IS1,IOBS) = DOBSEF(IL,IT,IS1,IOBS)/DOSEF
                  POBS(0,IT,IS1,IOBS) = POBS(0,IT,IS1,IOBS)
     &                                  + POBS(IL,IT,IS1,IOBS)
                  DO I = 1,IRTOP
                     POBSR(I,IL,IT,IS1,IOBS)
     &                  = RHORTLEF(I,IL,IT,IS1,IOBS)/DOSEF
                     POBSR(I,0,IT,IS1,IOBS) = POBSR(I,0,IT,IS1,IOBS)
     &                  + POBSR(I,IL,IT,IS1,IOBS)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      SUMWTZNS = 0D0
      DO IS2 = 1,2
         SUMWTZNO(IS2) = 0D0
      END DO
      SUMWTZNN = 0D0
      DO JT = 1,NT
         WJ = NAT(JT)*CONC(JT)
C
         WTZNS(JT) = WJ*TZ(JT,1,1,1,ISPN)
         WT0XNS(JT) = WJ*T0X(JT,1,1,1,ISPN)
         DO IS2 = 1,2
            WTZNO(JT,IS2) = WJ*TZ(JT,1,1,IS2,IORB)
            WT0XNO(JT,IS2) = WJ*T0X(JT,1,1,IS2,IORB)
         END DO
C
C---------- There is no charge coupling with external magnetic field
C---------- Therefore WTZNN(JT) = WJ*TZ(JT,1,1,1,IDOS) = 0
C
         WTZNN(JT) = 0.D0
         WT0XNN(JT) = WJ*T0X(JT,1,1,1,IDOS)
C
         WTIJXNS(JT) = 0D0
         DO IS2 = 1,2
            WTIJXNO(JT,IS2) = 0D0
         END DO
         WTIJXNN(JT) = 0D0
C
         DO JTP = 1,NT
C
            WJP = NAT(JTP)*CONC(JTP)
            WTIJXNS(JT) = WTIJXNS(JT) + WJP*TIJX(JTP,JT,1,1,1,ISPN)
            DO IS2 = 1,2
               WTIJXNO(JT,IS2) = WTIJXNO(JT,IS2)
     &                           + WJP*TIJX(JTP,JT,1,1,IS2,IORB)
            END DO
            WTIJXNN(JT) = WTIJXNN(JT) + WJP*TIJX(JTP,JT,1,1,1,IDOS)
C
         END DO
C
         SUMWTZNS = SUMWTZNS + WTZNS(JT)
         WTXNS(JT) = WT0XNS(JT) + WTIJXNS(JT)
         DO IS2 = 1,2
            SUMWTZNO(IS2) = SUMWTZNO(IS2) + WTZNO(JT,IS2)
            WTXNO(JT,IS2) = WT0XNO(JT,IS2) + WTIJXNO(JT,IS2)
         END DO
         SUMWTZNN = SUMWTZNN + WTZNN(JT)
         WTXNN(JT) = WT0XNN(JT) + WTIJXNN(JT)
C
      END DO
C
C---------- suppress terms connected with particle number N if requested
C
      IF ( .NOT.NTERM ) THEN
         SUMWTZNN = 0D0
         DO JT = 1,NT
            WTXNN(JT) = 0D0
         END DO
C
         NDIM = 3*NT
C
      ELSE IF ( NT.EQ.1 ) THEN
C
         NDIM = 3*NT
C
      ELSE
C
         NDIM = 4*NT
C
      END IF
C
      DO IT = 1,NT
C
C----contribution to high-field suscept due to shift of the Fermi energy
C
         WIS = POBS(0,IT,1,ISPN)
         WIO(1) = POBS(0,IT,1,IORB)
         WIO(2) = POBS(0,IT,2,IORB)
         WIN = POBS(0,IT,1,IDOS)
C
         BBHFS(IT) = WIS*SUMWTZNS
         BBHFO(IT) = WIS*(SUMWTZNO(1)+SUMWTZNO(2))
         BBHFN(IT) = WIS*SUMWTZNN
         DO JT = 1,NT
            AAHF(IT,JT) = -WIS*WTXNS(JT)
            DO IS2 = 1,2
               IOFF2 = IS2*NT
               AAHF(IT,IOFF2+JT) = -WIS*WTXNO(JT,IS2)
            END DO
C
            IF ( NT.EQ.1 ) THEN
               BBHFN(IT) = BBHFN(IT) + WIS*WTXNN(JT)
            ELSE
               AAHF(IT,3*NT+JT) = -WIS*WTXNN(JT)
            END IF
         END DO
C
         DO IS1 = 1,2
            IOFF1 = IS1*NT
            BBHFS(IOFF1+IT) = WIO(IS1)*SUMWTZNS
            BBHFO(IOFF1+IT) = WIO(IS1)*(SUMWTZNO(1)+SUMWTZNO(2))
            BBHFN(IOFF1+IT) = WIO(IS1)*SUMWTZNN
            DO JT = 1,NT
               AAHF(IOFF1+IT,JT) = -WIO(IS1)*WTXNS(JT)
               DO IS2 = 1,2
                  IOFF2 = IS2*NT
                  AAHF(IOFF1+IT,IOFF2+JT) = -WIO(IS1)*WTXNO(JT,IS2)
               END DO
C
               IF ( NT.EQ.1 ) THEN
                  BBHFN(IOFF1+IT) = BBHFN(IOFF1+IT) + WIO(IS1)*WTXNN(JT)
               ELSE
                  AAHF(IOFF1+IT,3*NT+JT) = -WIO(IS1)*WTXNN(JT)
               END IF
            END DO
         END DO
C
         IOFF1 = 3*NT
         BBHFS(IOFF1+IT) = WIN*SUMWTZNS
         BBHFO(IOFF1+IT) = WIN*(SUMWTZNO(1)+SUMWTZNO(2))
         BBHFN(IOFF1+IT) = WIN*SUMWTZNN
         DO JT = 1,NT
            AAHF(IOFF1+IT,JT) = -WIN*WTXNS(JT)
            DO IS2 = 1,2
               IOFF2 = IS2*NT
               AAHF(IOFF1+IT,IOFF2+JT) = -WIN*WTXNO(JT,IS2)
            END DO
C
            IF ( NT.EQ.1 ) THEN
               BBHFN(IOFF1+IT) = BBHFN(IOFF1+IT) + WIN*WTXNN(JT)
            ELSE
               AAHF(IOFF1+IT,3*NT+JT) = -WIN*WTXNN(JT)
            END IF
C            END DO
         END DO
C
C-------------- the terms due to spin-charge and orbital-charge coupling
C
         IOFF1 = 3*NT
         BBHFS(IOFF1+IT) = BBHFS(IOFF1+IT) - TZ(IT,1,IDOS,1,ISPN)
         BBHFO(IOFF1+IT) = BBHFO(IOFF1+IT)
     &                     - (TZ(IT,1,IDOS,1,IORB)+TZ(IT,1,IDOS,2,IORB))
         AAHF(3*NT+IT,3*NT+IT) = AAHF(3*NT+IT,3*NT+IT)
     &                           - (1D0-T0X(IT,1,IDOS,1,IDOS))
         DO JT = 1,NT
            AAHF(3*NT+IT,3*NT+JT) = AAHF(3*NT+IT,3*NT+JT)
     &                              + TIJX(IT,JT,1,IDOS,1,IDOS)
         END DO
C
         IF ( NT.EQ.1 ) THEN
            BBHFN(IT) = BBHFN(IT) - T0X(IT,1,ISPN,1,IDOS)
            DO JT = 1,NT
               BBHFN(IT) = BBHFN(IT) - TIJX(IT,JT,1,ISPN,1,IDOS)
            END DO
            DO IS1 = 1,2
               IOFF1 = IS1*NT
               BBHFN(IOFF1+IT) = BBHFN(IOFF1+IT)
     &                           - T0X(IT,IS1,IORB,1,IDOS)
C
               DO JT = 1,NT
                  BBHFN(IOFF1+IT) = BBHFN(IOFF1+IT)
     &                              - TIJX(IT,JT,IS1,IORB,1,IDOS)
               END DO
            END DO
C
         ELSE
C
            AAHF(IT,3*NT+IT) = AAHF(IT,3*NT+IT) + T0X(IT,1,ISPN,1,IDOS)
            AAHF(3*NT+IT,IT) = AAHF(3*NT+IT,IT) + T0X(IT,1,IDOS,1,ISPN)
            DO JT = 1,NT
               AAHF(IT,3*NT+JT) = AAHF(IT,3*NT+JT)
     &                            + TIJX(IT,JT,1,ISPN,1,IDOS)
               AAHF(3*NT+IT,JT) = AAHF(3*NT+IT,JT)
     &                            + TIJX(IT,JT,1,IDOS,1,ISPN)
            END DO
C
            DO IS1 = 1,2
               IOFF1 = IS1*NT
               AAHF(IOFF1+IT,3*NT+IT) = AAHF(IOFF1+IT,3*NT+IT)
     &                                  + T0X(IT,IS1,IORB,1,IDOS)
               AAHF(3*NT+IT,IOFF1+IT) = AAHF(3*NT+IT,IOFF1+IT)
     &                                  + T0X(IT,1,IDOS,IS1,IORB)
C
               DO JT = 1,NT
                  AAHF(IOFF1+IT,3*NT+JT) = AAHF(IOFF1+IT,3*NT+JT)
     &               + TIJX(IT,JT,IS1,IORB,1,IDOS)
                  AAHF(3*NT+IT,IOFF1+JT) = AAHF(3*NT+IT,IOFF1+JT)
     &               + TIJX(IT,JT,1,IDOS,IS1,IORB)
               END DO
            END DO
C
         END IF
C
      END DO
C
C ======================================================================
C ======================================================================
C     calculate atom type-resolved  SPIN  and  ORBITAL  susceptibility
C ======================================================================
C ======================================================================
C                                                               S-O-LOOP
      DO ISOCOUPLING = 1,2
C
C--------------------------- set up and solve linear system of equations
C
         CALL RINIT(NTMAX*4*NTMAX*4,AA)
         CALL RINIT(NTMAX*4,BB)
C
         DO IT = 1,NT
C
            BB(IT,1) = TZ(IT,1,ISPN,1,ISPN)
            AA(IT,IT) = 1D0 - T0X(IT,1,ISPN,1,ISPN)
            DO JT = 1,NT
               AA(IT,JT) = AA(IT,JT) - TIJX(IT,JT,1,ISPN,1,ISPN)
            END DO
C
            DO IS1 = 1,2
               IOFF1 = IS1*NT
               BB(IOFF1+IT,1) = TZ(IT,IS1,IORB,1,IORB)
     &                          + TZ(IT,IS1,IORB,2,IORB)
               AA(IOFF1+IT,IOFF1+IT) = 1D0
               DO IS2 = 1,2
                  IOFF2 = IS2*NT
                  AA(IOFF1+IT,IOFF2+IT) = AA(IOFF1+IT,IOFF2+IT)
     &               - T0X(IT,IS1,IORB,IS2,IORB)
C
                  DO JT = 1,NT
                     AA(IOFF1+IT,IOFF2+JT) = AA(IOFF1+IT,IOFF2+JT)
     &                  - TIJX(IT,JT,IS1,IORB,IS2,IORB)
                  END DO
               END DO
C
            END DO
C
         END DO
C
C---------------- exchange coupling of spin and orbital susceptibilities
C-------------- if ignored:  the coefficient matrix AA is block diagonal
C
         IF ( ISOCOUPLING.EQ.1 ) THEN
            WRITE (6,99020) 'OFF'
         ELSE
            WRITE (6,99020) 'ON'
C
            DO IT = 1,NT
               DO IS2 = 1,2
                  IOFF2 = IS2*NT
                  BB(IT,1) = BB(IT,1) + TZ(IT,1,ISPN,IS2,IORB)
                  AA(IT,IOFF2+IT) = AA(IT,IOFF2+IT)
     &                              - T0X(IT,1,ISPN,IS2,IORB)
                  DO JT = 1,NT
                     AA(IT,IOFF2+JT) = AA(IT,IOFF2+JT)
     &                                 - TIJX(IT,JT,1,ISPN,IS2,IORB)
                  END DO
               END DO
C
               DO IS1 = 1,2
                  IOFF1 = IS1*NT
                  BB(IOFF1+IT,1) = BB(IOFF1+IT,1)
     &                             + TZ(IT,IS1,IORB,1,ISPN)
                  AA(IOFF1+IT,IT) = AA(IOFF1+IT,IT)
     &                              - T0X(IT,IS1,IORB,1,ISPN)
                  DO JT = 1,NT
                     AA(IOFF1+IT,JT) = AA(IOFF1+IT,JT)
     &                                 - TIJX(IT,JT,IS1,IORB,1,ISPN)
                  END DO
C
               END DO
            END DO
C
         END IF
C
C------------------- add  terms connected with high-field susceptibility
C
         IF ( .NOT.NONMAG ) THEN
C
C------------------------------- the terms created by Fermi energy shift
            IF ( ISOCOUPLING.EQ.1 ) THEN
C
               DO I = 1,NT
                  BB(I,1) = BB(I,1) - BBHFS(I) - BBHFO(I) - BBHFN(I)
                  DO J = 1,NT
                     AA(I,J) = AA(I,J) - AAHF(I,J)
                  END DO
               END DO
C
               DO I = NT + 1,NDIM
                  BB(I,1) = BB(I,1) - BBHFS(I) - BBHFO(I) - BBHFN(I)
                  DO J = NT + 1,NDIM
                     AA(I,J) = AA(I,J) - AAHF(I,J)
                  END DO
               END DO
C
            ELSE
C
               DO I = 1,NDIM
                  BB(I,1) = BB(I,1) - BBHFS(I) - BBHFO(I) - BBHFN(I)
                  DO J = 1,NDIM
                     AA(I,J) = AA(I,J) - AAHF(I,J)
                  END DO
               END DO
C
            END IF
C
         END IF
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
         CALL DGESV(NDIM,1,AA,NTMAX*4,IPIV,BB,NTMAX*4,INFO)
C
         IF ( INFO.NE.0 ) STOP 'in <CHICALC> after ZGESV'
C
C-----------------------------------------------------------------------
C----------------------------------- determine enhanced susceptibilities
C
         DO IT = 1,NT
C
            CHIS(IT) = BB(IT,1)
C
            CHIO(IT,0) = 0D0
            DO IS1 = 1,2
               IOFF1 = IS1*NT
               CHIO(IT,IS1) = BB(IOFF1+IT,1)
               CHIO(IT,0) = CHIO(IT,0) + CHIO(IT,IS1)
            END DO
C
            CHIN(IT) = BB(3*NT+IT,1)
C
         END DO
C
C--------------------------------------- determine shift in Fermi energy
C
         IF ( .NOT.NONMAG ) THEN
C
            SUMEF = SUMWTZNS
            DO IS2 = 1,2
               SUMEF = SUMEF + SUMWTZNO(IS2)
            END DO
            SUMEF = SUMEF + SUMWTZNN
C
            DO JT = 1,NT
C
               SUMEF = SUMEF + WTXNS(JT)*CHIS(JT)
C
               DO IS2 = 1,2
                  SUMEF = SUMEF + WTXNO(JT,IS2)*CHIO(JT,IS2)
               END DO
               IF ( NT.EQ.1 ) THEN
                  SUMEF = SUMEF + WTXNN(JT)
               ELSE
                  SUMEF = SUMEF + WTXNN(JT)*CHIN(JT)
               END IF
C             SUMEF = SUMEF + WTXNN(JT)*CHIN(JT)
            END DO
C
            DELTAEF = -SUMEF/DOSEF
C
         END IF
C
         IF ( .NOT.NONMAG ) WRITE (6,99021) ' WITHOUT '
C
C ======================================================================
C             determine remaining susceptibility terms
C ======================================================================
C
         DO IT = 1,NT
C
C-------------------------------------------------------------------- SS
C
            CHISS0(IT) = TZ(IT,1,ISPN,1,ISPN)
            CHISS(IT) = CHISS0(IT) + T0X(IT,1,ISPN,1,ISPN)*CHIS(IT)
C
            DO JT = 1,NT
               CHISS(IT) = CHISS(IT) + TIJX(IT,JT,1,ISPN,1,ISPN)
     &                     *CHIS(JT)
               CHISSTR(IT,JT) = TIJXTR(IT,JT,1,ISPN,1,ISPN)*CHIS(JT)
            END DO
            ICHIS(IT) = 1.D0 - CHISS0(IT)/CHISS(IT)
C
            DO IL = 1,NL
               CHISSL(IL,IT) = TZL(IL,IT,1,ISPN,1,ISPN)
     &                         + T0XL(IL,IT,1,ISPN,1,ISPN)*CHIS(IT)
               DO JT = 1,NT
                  CHISSL(IL,IT) = CHISSL(IL,IT)
     &                            + TIJXL(IL,IT,JT,1,ISPN,1,ISPN)
     &                            *CHIS(JT)
                  CHISSLTR(IL,IT,JT) = TIJXLTR(IL,IT,JT,1,ISPN,1,ISPN)
     &                                 *CHIS(JT)
               END DO
            END DO
C
C-------------------------------------------------------------------- OO
C
            CHIOO0(IT,0) = 0D0
            CHIOO(IT,0) = 0D0
            DO JT = 1,NT
               DO IS1 = 1,2
                  CHIOOTR(IT,JT,IS1) = 0D0
                  DO IL = 1,NLT(IT)
                     CHIOOLTR(IL,IT,JT,IS1) = 0D0
                  END DO
               END DO
            END DO
C
            DO IS1 = 1,2
               CHIOO0(IT,IS1) = TZ(IT,IS1,IORB,1,IORB)
     &                          + TZ(IT,IS1,IORB,2,IORB)
C
               CHIOO(IT,IS1) = CHIOO0(IT,IS1)
C
               DO IS2 = 1,2
                  CHIOO(IT,IS1) = CHIOO(IT,IS1)
     &                            + T0X(IT,IS1,IORB,IS2,IORB)
     &                            *CHIO(IT,IS2)
                  DO JT = 1,NT
                     CHIOO(IT,IS1) = CHIOO(IT,IS1)
     &                               + TIJX(IT,JT,IS1,IORB,IS2,IORB)
     &                               *CHIO(JT,IS2)
                     CHIOOTR(IT,JT,IS1) = CHIOOTR(IT,JT,IS1)
     &                  + TIJXTR(IT,JT,IS1,IORB,IS2,IORB)*CHIO(JT,IS2)
                  END DO
               END DO
C
               CHIOO0(IT,0) = CHIOO0(IT,0) + CHIOO0(IT,IS1)
               CHIOO(IT,0) = CHIOO(IT,0) + CHIOO(IT,IS1)
               DO JT = 1,NT
                  CHIOOTR(IT,JT,0) = CHIOOTR(IT,JT,0)
     &                               + CHIOOTR(IT,JT,IS1)
               END DO
            END DO
C
            ICHIO(IT) = 1.D0 - CHIOO0(IT,0)/CHIOO(IT,0)
C
            DO IL = 1,NLT(IT)
               CHIOOL(IL,IT,0) = 0D0
               DO IS1 = 1,2
                  CHIOOL(IL,IT,IS1) = 0D0
C
                  DO IS2 = 1,2
                     CHIOOL(IL,IT,IS1) = CHIOOL(IL,IT,IS1)
     &                  + TZL(IL,IT,IS1,IORB,IS2,IORB)
     &                  + T0XL(IL,IT,IS1,IORB,IS2,IORB)*CHIO(IT,IS2)
                     DO JT = 1,NT
                        CHIOOL(IL,IT,IS1) = CHIOOL(IL,IT,IS1)
     &                     + TIJXL(IL,IT,JT,IS1,IORB,IS2,IORB)
     &                     *CHIO(JT,IS2)
                        CHIOOLTR(IL,IT,JT,IS1) = CHIOOLTR(IL,IT,JT,IS1)
     &                     + TIJXLTR(IL,IT,JT,IS1,IORB,IS2,IORB)
     &                     *CHIO(JT,IS2)
                     END DO
                  END DO
C
                  CHIOOL(IL,IT,0) = CHIOOL(IL,IT,0) + CHIOOL(IL,IT,IS1)
                  DO JT = 1,NT
                     CHIOOLTR(IL,IT,JT,0) = CHIOOLTR(IL,IT,JT,0)
     &                  + CHIOOLTR(IL,IT,JT,IS1)
                  END DO
               END DO
            END DO
C
C-------------------------------------------------------------------- SO
C
            CHISO0(IT) = 0D0
            CHISO(IT) = 0D0
            DO IS2 = 1,2
C
               CHISO0(IT) = CHISO0(IT) + TZ(IT,1,ISPN,IS2,IORB)
C
               CHISO(IT) = CHISO(IT) + TZ(IT,1,ISPN,IS2,IORB)
     &                     + T0X(IT,1,ISPN,IS2,IORB)*CHIO(IT,IS2)
               DO JT = 1,NT
C
                  CHISO(IT) = CHISO(IT) + TIJX(IT,JT,1,ISPN,IS2,IORB)
     &                        *CHIO(JT,IS2)
               END DO
            END DO
C
C-------------------------------------------------------------------- OS
C
            CHIOS0(IT,0) = 0D0
            CHIOS(IT,0) = 0D0
            DO IS1 = 1,2
               CHIOS0(IT,IS1) = TZ(IT,IS1,IORB,1,ISPN)
               CHIOS(IT,IS1) = TZ(IT,IS1,IORB,1,ISPN)
     &                         + T0X(IT,IS1,IORB,1,ISPN)*CHIS(IT)
               DO JT = 1,NT
                  CHIOS(IT,IS1) = CHIOS(IT,IS1)
     &                            + TIJX(IT,JT,IS1,IORB,1,ISPN)*CHIS(JT)
               END DO
               CHIOS(IT,0) = CHIOS(IT,0) + CHIOS(IT,IS1)
               CHIOS0(IT,0) = CHIOS0(IT,0) + CHIOS0(IT,IS1)
            END DO
C
C------------------------------------------------- radial magnetisations
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
C--------------------------------------- radial spin magnetisation CHISR
C--- SS contribution
C
            CALL DCOPY(IRTOP,TZR(1,IT,1,ISPN,1,ISPN),1,CHISR(1,IT),1)
            CALL DAXPY(IRTOP,CHIS(IT),T0XR(1,IT,1,ISPN,1,ISPN),1,
     &                 CHISR(1,IT),1)
C
            DO JT = 1,NT
               CALL DAXPY(IRTOP,CHIS(JT),TIJXR(1,IT,JT,1,ISPN,1,ISPN),1,
     &                    CHISR(1,IT),1)
            END DO
C
C--- SO contribution
C
            DO IS2 = 1,2
               CALL DAXPY(IRTOP,1D0,TZR(1,IT,1,ISPN,IS2,IORB),1,
     &                    CHISR(1,IT),1)
               CALL DAXPY(IRTOP,CHIO(IT,IS2),T0XR(1,IT,1,ISPN,IS2,IORB),
     &                    1,CHISR(1,IT),1)
               DO JT = 1,NT
                  CALL DAXPY(IRTOP,CHIO(JT,IS2),
     &                       TIJXR(1,IT,JT,1,ISPN,IS2,IORB),1,
     &                       CHISR(1,IT),1)
               END DO
            END DO
C
            DO IS1 = 1,2
C
C------------------------------------ radial orbital magnetisation CHIOR
C--- OS contribution
C
               CALL DCOPY(IRTOP,TZR(1,IT,IS1,IORB,1,ISPN),1,
     &                    CHIOR(1,IT,IS1),1)
               CALL DAXPY(IRTOP,CHIS(IT),T0XR(1,IT,IS1,IORB,1,ISPN),1,
     &                    CHIOR(1,IT,IS1),1)
               DO JT = 1,NT
                  CALL DAXPY(IRTOP,CHIS(JT),
     &                       TIJXR(1,IT,JT,IS1,IORB,1,ISPN),1,
     &                       CHIOR(1,IT,IS1),1)
               END DO
C
C--- OO contribution
C
               DO IS2 = 1,2
                  CALL DAXPY(IRTOP,1D0,TZR(1,IT,IS1,IORB,IS2,IORB),1,
     &                       CHIOR(1,IT,IS1),1)
                  CALL DAXPY(IRTOP,CHIO(IT,IS2),
     &                       T0XR(1,IT,IS1,IORB,IS2,IORB),1,
     &                       CHIOR(1,IT,IS1),1)
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIO(JT,IS2),
     &                          TIJXR(1,IT,JT,IS1,IORB,IS2,IORB),1,
     &                          CHIOR(1,IT,IS1),1)
                  END DO
               END DO
C
            END DO
C
            CALL DCOPY(IRTOP,CHIOR(1,IT,1),1,CHIOR(1,IT,0),1)
            CALL DAXPY(IRTOP,1D0,CHIOR(1,IT,2),1,CHIOR(1,IT,0),1)
C
C------------------------------------ orbital polarisation related terms
C
            DO I = 1,IRTOP
               RINT(I) = TSSZRLOP(I,IT)*R2DRDI(I,IM)
            END DO
C
            CALL RRADINT(IM,RINT,NORM)
C
            DO I = 1,IRTOP
               RSSZLOP(I,IT) = TSSZRLOP(I,IT)/NORM
            END DO
C
            TOL = 1D-8
C
C----------------------------------------------- consistency checks spin
C
            IF ( CHIPRINT.GE.1 ) THEN
               IF ( ISOCOUPLING.EQ.1 ) THEN
                  CHISUM = CHISS(IT)
               ELSE
                  CHISUM = CHISS(IT) + CHISO(IT)
               END IF
C
               IF ( ABS(CHIS(IT)-CHISUM).GT.TOL ) WRITE (6,99001) 'SUM',
     &              'spin',IT,TXT_T(IT),CHIS(IT),CHISUM,
     &              (1D0-CHIS(IT)/CHISUM)
C
               CHISUM = 0D0
               DO IL = 1,NL
                  CHISUM = CHISUM + CHISSL(IL,IT)
               END DO
               IF ( ABS(CHISS(IT)-CHISUM).GT.TOL ) WRITE (6,99001)
     &               'L-sum','spin',IT,TXT_T(IT),CHISS(IT),CHISUM,
     &              (1D0-CHISS(IT)/CHISUM)
C
               DO I = 1,IRTOP
                  RINT(I) = CHISR(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RRADINT(IM,RINT,CHISUM)
               IF ( ISOCOUPLING.EQ.1 ) CHISUM = CHISUM - CHISO(IT)
               IF ( ABS(CHIS(IT)-CHISUM).GT.TOL ) WRITE (6,99001)
     &               'R-int','spin',IT,TXT_T(IT),CHIS(IT),CHISUM,
     &              (1D0-CHIS(IT)/CHISUM)
C
C-------------------------------------------- consistency checks orbital
C
               DO IS = 1,2
                  IF ( ISOCOUPLING.EQ.1 ) THEN
                     CHISUM = CHIOO(IT,IS)
                  ELSE
                     CHISUM = CHIOO(IT,IS) + CHIOS(IT,IS)
                  END IF
C
                  IF ( ABS(CHIO(IT,IS)-CHISUM).GT.TOL ) WRITE (6,99002)
     &                  'SUM','orbital',IT,TXT_T(IT),IS,CHIO(IT,IS),
     &                 CHISUM,(1D0-CHIO(IT,IS)/CHISUM)
C
                  CHISUM = 0D0
                  DO IL = 1,NL
                     CHISUM = CHISUM + CHIOOL(IL,IT,IS)
                  END DO
                  IF ( ABS(CHIOO(IT,IS)-CHISUM).GT.TOL ) WRITE (6,99002)
     &                  'L-sum','orbital  IT=',IT,TXT_T(IT),IS,
     &                 CHIOO(IT,IS),CHISUM,(1D0-CHIOO(IT,IS)/CHISUM)
C
                  DO I = 1,IRTOP
                     RINT(I) = CHIOR(I,IT,IS)*R2DRDI(I,IM)
                  END DO
                  CALL RRADINT(IM,RINT,CHISUM)
                  IF ( ISOCOUPLING.EQ.1 ) CHISUM = CHISUM - CHIOS(IT,IS)
                  IF ( ABS(CHIO(IT,IS)-CHISUM).GT.TOL ) WRITE (6,99002)
     &                  'R-int','orbital',IT,TXT_T(IT),IS,CHIO(IT,IS),
     &                 CHISUM,(1D0-CHIO(IT,IS)/CHISUM)
               END DO
            END IF
C
         END DO
C
C-----------------------------------------------------------------------
         IF ( CHIPRINT.GE.1 ) THEN
            WRITE (6,99008) (IT,T1Z(IT,1,ISPN,1,ISPN),T0Z(IT,1,ISPN,1,
     &                      ISPN),IT=1,NT)
            WRITE (6,99009) (IT,T0X(IT,1,ISPN,1,ISPN),IT=1,NT)
C
            DO IT = 1,NT
               WRITE (6,'('' IT='',I2)') IT
               WRITE (6,99010) (JT,'Z',TIJZ(IT,JT,1,ISPN,1,ISPN),JT=1,
     &                         NT)
               WRITE (6,99010) (JT,'X',TIJX(IT,JT,1,ISPN,1,ISPN),JT=1,
     &                         NT)
            END DO
C
            WRITE (6,99011) ((IT,IS1,T1Z(IT,IS1,IORB,IS1,IORB),T0Z(IT,
     &                      IS1,IORB,IS1,IORB),IS1=1,2),IT=1,NT)
            WRITE (6,99012) (((IT,IS1,IS2,T0X(IT,IS1,IORB,IS2,IORB),IS1=
     &                      1,2),IS2=1,2),IT=1,NT)
C
            DO IT = 1,NT
               WRITE (6,'('' IT='',I2)') IT
               WRITE (6,99013) ((JT,IS1,TIJZ(IT,JT,IS1,IORB,IS1,IORB),
     &                         IS1=1,2),JT=1,NT)
               WRITE (6,99014) (((JT,IS1,IS2,TIJX(IT,JT,IS1,IORB,IS2,
     &                         IORB),IS1=1,2),IS2=1,2),JT=1,NT)
            END DO
         END IF
C
C
C ======================================================================
C                    calculate Knight shift
C ======================================================================
C
         DO IT = 1,NT
C
C------------------------------------------------------------------ HF-S
C
            SFTS0(IT) = TZ(IT,1,IHFI,1,ISPN)
            SFTS(IT) = TZ(IT,1,IHFI,1,ISPN) + T0X(IT,1,IHFI,1,ISPN)
     &                 *CHIS(IT)
            DO JT = 1,NT
               SFTS(IT) = SFTS(IT) + TIJX(IT,JT,1,IHFI,1,ISPN)*CHIS(JT)
               SFTSTR(IT,JT) = TIJXTR(IT,JT,1,IHFI,1,ISPN)*CHIS(JT)
            END DO
C
            DO IL = 1,NLT(IT)
               SFTSL0(IL,IT) = TZL(IL,IT,1,IHFI,1,ISPN)
               SFTSL(IL,IT) = TZL(IL,IT,1,IHFI,1,ISPN)
     &                        + T0XL(IL,IT,1,IHFI,1,ISPN)*CHIS(IT)
               DO JT = 1,NT
                  SFTSL(IL,IT) = SFTSL(IL,IT)
     &                           + TIJXL(IL,IT,JT,1,IHFI,1,ISPN)
     &                           *CHIS(JT)
                  SFTSLTR(IL,IT,JT) = +TIJXLTR(IL,IT,JT,1,IHFI,1,ISPN)
     &                                *CHIS(JT)
               END DO
            END DO
C
C------------------------------------------------------------------ HF-O
C
            SFTO0(IT) = 0D0
            SFTO(IT) = 0D0
            DO JT = 1,NT
               SFTOTR(IT,JT) = 0D0
               DO IL = 1,NLT(IT)
                  SFTOLTR(IL,IT,JT) = 0D0
               END DO
            END DO
            DO IS2 = 1,2
               SFTO0(IT) = SFTO0(IT) + TZ(IT,1,IHFI,IS2,IORB)
               SFTO(IT) = SFTO(IT) + TZ(IT,1,IHFI,IS2,IORB)
     &                    + T0X(IT,1,IHFI,IS2,IORB)*CHIO(IT,IS2)
               DO JT = 1,NT
                  SFTO(IT) = SFTO(IT) + TIJX(IT,JT,1,IHFI,IS2,IORB)
     &                       *CHIO(JT,IS2)
                  SFTOTR(IT,JT) = SFTOTR(IT,JT)
     &                            + TIJXTR(IT,JT,1,IHFI,IS2,IORB)
     &                            *CHIO(JT,IS2)
               END DO
            END DO
C
            DO IL = 1,NLT(IT)
               SFTOL0(IL,IT) = 0D0
               SFTOL(IL,IT) = 0D0
               DO IS2 = 1,2
                  SFTOL0(IL,IT) = SFTOL0(IL,IT)
     &                            + TZL(IL,IT,1,IHFI,IS2,IORB)
                  SFTOL(IL,IT) = SFTOL(IL,IT)
     &                           + TZL(IL,IT,1,IHFI,IS2,IORB)
     &                           + T0XL(IL,IT,1,IHFI,IS2,IORB)
     &                           *CHIO(IT,IS2)
                  DO JT = 1,NT
                     SFTOL(IL,IT) = SFTOL(IL,IT)
     &                              + TIJXL(IL,IT,JT,1,IHFI,IS2,IORB)
     &                              *CHIO(JT,IS2)
                     SFTOLTR(IL,IT,JT) = SFTOLTR(IL,IT,JT)
     &                  + TIJXLTR(IL,IT,JT,1,IHFI,IS2,IORB)*CHIO(JT,IS2)
                  END DO
               END DO
            END DO
C
         END DO
C
C ======================================================================
C            print out of SUSCEPTIBILITY and KNIGHT SHIFT
C                         WITHOUT HF-terms
C ======================================================================
C
C
         CALL CHIOUTPUT(LTXT_T,TXT_T,DATSET,LDATSET,BCP,BCPS,CHIO,CHIOO,
     &                  CHIOO0,CHIOOL,CHIOOLTR,CHIOOTR,CHIOR,CHIOS,
     &                  CHIOS0,CHIPRINT,CHIS,CHISO,CHISO0,CHISR,CHISS,
     &                  CHISS0,CHISSL,CHIN,CHINN,CHINN0,CHINNL,
     &                  CHILANDAU,CHILANLT,CHISSLTR,CHISSTR,CONC,DOBSEF,
     &                  ICHIO,ICHIS,IDOS,IHFI,IHVV,IMT,IORB,ISPN,JRWS,
     &                  NAT,NL,NLT,NSPINOBS,NT,NONMAG,R,SFTO,SFTO0,
     &                  SFTOL,SFTOL0,SFTOLTR,SFTOTR,SFTS,SFTS0,SFTSL,
     &                  SFTSL0,SFTSLTR,SFTSTR,TCDIACRI,TCDIARI,TKDIACRI,
     &                  TKDIARI,TOBS,TXTS,TZ,TZL,TZR,ICDIA,IKDIA,IRM3,
     &                  NOP,NOBS,NOBSDNS,NLMAX,NMMAX,NRMAX,NTMAX,IOTMP)
C
C
C ======================================================================
C                                      for FERRO add additional HF-terms
         IF ( .NOT.NONMAG ) THEN
C
            WRITE (6,99021) 'INCLUDING'
C
C ======================================================================
C             determine remaining susceptibility terms
C ======================================================================
C
            DO IT = 1,NT
C
C--------------------------------------------------------- high field SS
C
               WIS = POBS(0,IT,1,ISPN)
               CHISS0(IT) = CHISS0(IT) - WIS*SUMWTZNS
               CHISS(IT) = CHISS(IT) - WIS*SUMWTZNS
C
               DO JT = 1,NT
                  CHISS(IT) = CHISS(IT) - WIS*WTXNS(JT)*CHIS(JT)
               END DO
C
               DO IL = 1,NL
                  WIS = POBS(IL,IT,1,ISPN)
                  CHISSL(IL,IT) = CHISSL(IL,IT) - WIS*SUMWTZNS
                  DO JT = 1,NT
                     CHISSL(IL,IT) = CHISSL(IL,IT) - WIS*WTXNS(JT)
     &                               *CHIS(JT)
                  END DO
               END DO
C
C--------------------------------------------------------- high field OO
C
               CHIOO0(IT,0) = 0D0
               CHIOO(IT,0) = 0D0
C
               DO IS1 = 1,2
                  WIO(IS1) = POBS(0,IT,IS1,IORB)
                  DO IS2 = 1,2
                     CHIOO0(IT,IS1) = CHIOO0(IT,IS1) - WIO(IS1)
     &                                *SUMWTZNO(IS2)
                     CHIOO(IT,IS1) = CHIOO(IT,IS1) - WIO(IS1)
     &                               *SUMWTZNO(IS2)
                  END DO
                  DO JT = 1,NT
                     DO IS2 = 1,2
                        CHIOO(IT,IS1) = CHIOO(IT,IS1) - WIO(IS1)
     &                                  *WTXNO(JT,IS2)*CHIO(JT,IS2)
                     END DO
                  END DO
C
                  CHIOO0(IT,0) = CHIOO0(IT,0) + CHIOO0(IT,IS1)
                  CHIOO(IT,0) = CHIOO(IT,0) + CHIOO(IT,IS1)
               END DO
C
               DO IL = 1,NLT(IT)
                  CHIOOL(IL,IT,0) = 0D0
                  DO IS1 = 1,2
                     WIO(IS1) = POBS(IL,IT,IS1,IORB)
                     DO IS2 = 1,2
                        CHIOOL(IL,IT,IS1) = CHIOOL(IL,IT,IS1) - WIO(IS1)
     &                     *SUMWTZNO(IS2)
                     END DO
                     DO JT = 1,NT
                        DO IS2 = 1,2
                           CHIOOL(IL,IT,IS1) = CHIOOL(IL,IT,IS1)
     &                        - WIO(IS1)*WTXNO(JT,IS2)*CHIO(JT,IS2)
                        END DO
                     END DO
                     CHIOOL(IL,IT,0) = CHIOOL(IL,IT,0)
     &                                 + CHIOOL(IL,IT,IS1)
                  END DO
               END DO
C
C--------------------------------------------------------- high field SO
C
               WIS = POBS(0,IT,1,ISPN)
C
               DO IS2 = 1,2
                  CHISO0(IT) = CHISO0(IT) - WIS*SUMWTZNO(IS2)
                  CHISO(IT) = CHISO(IT) - WIS*SUMWTZNO(IS2)
C
                  DO JT = 1,NT
                     CHISO(IT) = CHISO(IT) - WIS*WTXNO(JT,IS2)
     &                           *CHIO(JT,IS2)
                  END DO
               END DO
C
C--------------------------------------------------------- high field OS
C
               CHIOS(IT,0) = 0D0
               CHIOS0(IT,0) = 0D0
C
               DO IS1 = 1,2
                  WIO(IS1) = POBS(0,IT,IS1,IORB)
C
                  CHIOS0(IT,IS1) = CHIOS0(IT,IS1) - WIO(IS1)*SUMWTZNS
                  CHIOS(IT,IS1) = CHIOS(IT,IS1) - WIO(IS1)*SUMWTZNS
C
                  DO JT = 1,NT
                     CHIOS(IT,IS1) = CHIOS(IT,IS1) - WIO(IS1)*WTXNS(JT)
     &                               *CHIS(JT)
                  END DO
C
                  CHIOS(IT,0) = CHIOS(IT,0) + CHIOS(IT,IS1)
                  CHIOS0(IT,0) = CHIOS0(IT,0) + CHIOS0(IT,IS1)
               END DO
C
C--------------------------------------------------------- high field NN
C
               WIN = POBS(0,IT,1,IDOS)
               CHINN0(IT) = -WIN*SUMWTZNN
               CHINN(IT) = CHINN0(IT)
               IF ( NT.EQ.1 ) THEN
                  CHINN(IT) = CHINN(IT) + T0X(IT,1,IDOS,1,IDOS)
                  DO JT = 1,NT
                     CHINN(IT) = CHINN(IT) - WIN*WTXNN(JT)
                     CHINN(IT) = CHINN(IT) + TIJX(IT,JT,1,IDOS,1,IDOS)
                  END DO
               ELSE
                  CHINN(IT) = CHINN(IT) + T0X(IT,1,IDOS,1,IDOS)*CHIN(IT)
                  DO JT = 1,NT
                     CHINN(IT) = CHINN(IT) - WIN*WTXNN(JT)*CHIN(JT)
                     CHINN(IT) = CHINN(IT) + TIJX(IT,JT,1,IDOS,1,IDOS)
     &                           *CHIN(JT)
                  END DO
               END IF
C
C   angular moment dependent suscept.
C
               DO IL = 1,NL
                  WIN = POBS(IL,IT,1,IDOS)
                  CHINNL(IL,IT) = -WIN*SUMWTZNN
                  IF ( NT.EQ.1 ) THEN
                     CHINNL(IL,IT) = CHINNL(IL,IT)
     &                               + T0XL(IL,IT,1,IDOS,1,IDOS)
                     DO JT = 1,NT
                        CHINNL(IL,IT) = CHINNL(IL,IT)
     &                                  + TIJXL(IL,IT,JT,1,IDOS,1,IDOS)
                        CHINNL(IL,IT) = CHINNL(IL,IT) - WIN*WTXNN(JT)
                     END DO
                  ELSE
                     CHINNL(IL,IT) = CHINNL(IL,IT)
     &                               + T0XL(IL,IT,1,IDOS,1,IDOS)
     &                               *CHIN(IT)
                     DO JT = 1,NT
                        CHINNL(IL,IT) = CHINNL(IL,IT)
     &                                  + TIJXL(IL,IT,JT,1,IDOS,1,IDOS)
     &                                  *CHIN(JT)
                        CHINNL(IL,IT) = CHINNL(IL,IT) - WIN*WTXNN(JT)
     &                                  *CHIN(JT)
                     END DO
                  END IF
               END DO
C
C------------------------------------------------- radial magnetisations
C
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
C--------------------------------------- radial spin magnetisation CHISR
C--- SS contribution
C
               CALL DCOPY(IRTOP,TZR(1,IT,1,ISPN,1,ISPN),1,CHISR(1,IT),1)
               CALL DAXPY(IRTOP,CHIS(IT),T0XR(1,IT,1,ISPN,1,ISPN),1,
     &                    CHISR(1,IT),1)
C
               DO JT = 1,NT
                  CALL DAXPY(IRTOP,CHIS(JT),TIJXR(1,IT,JT,1,ISPN,1,ISPN)
     &                       ,1,CHISR(1,IT),1)
               END DO
C
C--- SO contribution
C
               DO IS2 = 1,2
                  CALL DAXPY(IRTOP,1D0,TZR(1,IT,1,ISPN,IS2,IORB),1,
     &                       CHISR(1,IT),1)
                  CALL DAXPY(IRTOP,CHIO(IT,IS2),
     &                       T0XR(1,IT,1,ISPN,IS2,IORB),1,CHISR(1,IT),1)
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIO(JT,IS2),
     &                          TIJXR(1,IT,JT,1,ISPN,IS2,IORB),1,
     &                          CHISR(1,IT),1)
                  END DO
               END DO
C
C--- SN contribution (different for the cases NT=1 and NT>1)
C
               IF ( NT.EQ.1 ) THEN
                  CALL DAXPY(IRTOP,1.D0,T0XR(1,IT,1,ISPN,1,IDOS),1,
     &                       CHISR(1,IT),1)
C
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,1.D0,TIJXR(1,IT,JT,1,ISPN,1,IDOS),
     &                          1,CHISR(1,IT),1)
                  END DO
               ELSE
                  CALL DAXPY(IRTOP,CHIN(IT),T0XR(1,IT,1,ISPN,1,IDOS),1,
     &                       CHISR(1,IT),1)
C
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIN(JT),
     &                          TIJXR(1,IT,JT,1,ISPN,1,IDOS),1,
     &                          CHISR(1,IT),1)
                  END DO
               END IF
C
C--- contributions due to \Delta E_FERMI
C
               DO I = 1,IRTOP
                  WIS = POBSR(I,0,IT,1,ISPN)
                  CHISR(I,IT) = CHISR(I,IT) - WIS*SUMWTZNS - 
     &                          WIS*(SUMWTZNO(1)+SUMWTZNO(2))
                  DO JT = 1,NT
                     CHISR(I,IT) = CHISR(I,IT) - WIS*WTXNO(JT,1)
     &                             *CHIO(JT,1) - WIS*WTXNO(JT,2)
     &                             *CHIO(JT,2) - WIS*WTXNS(JT)*CHIS(JT)
                     IF ( NT.EQ.1 ) THEN
                        CHISR(I,IT) = CHISR(I,IT) - WIS*WTXNN(JT)
                     ELSE
                        CHISR(I,IT) = CHISR(I,IT) - WIS*WTXNN(JT)
     &                                *CHIN(JT)
                     END IF
                  END DO
C
               END DO
C
               DO IS1 = 1,2
C
C------------------------------------ radial orbital magnetisation CHIOR
C--- OS contribution
C
                  CALL DCOPY(IRTOP,TZR(1,IT,IS1,IORB,1,ISPN),1,
     &                       CHIOR(1,IT,IS1),1)
                  CALL DAXPY(IRTOP,CHIS(IT),T0XR(1,IT,IS1,IORB,1,ISPN),
     &                       1,CHIOR(1,IT,IS1),1)
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIS(JT),
     &                          TIJXR(1,IT,JT,IS1,IORB,1,ISPN),1,
     &                          CHIOR(1,IT,IS1),1)
                  END DO
C
C--- OO contribution
C
                  DO IS2 = 1,2
                     CALL DAXPY(IRTOP,1D0,TZR(1,IT,IS1,IORB,IS2,IORB),1,
     &                          CHIOR(1,IT,IS1),1)
                     CALL DAXPY(IRTOP,CHIO(IT,IS2),
     &                          T0XR(1,IT,IS1,IORB,IS2,IORB),1,
     &                          CHIOR(1,IT,IS1),1)
                     DO JT = 1,NT
                        CALL DAXPY(IRTOP,CHIO(JT,IS2),
     &                             TIJXR(1,IT,JT,IS1,IORB,IS2,IORB),1,
     &                             CHIOR(1,IT,IS1),1)
                     END DO
                  END DO
C
C--- ON contribution (different for the cases NT=1 and NT>1)
C
                  IF ( NT.EQ.1 ) THEN
                     CALL DAXPY(IRTOP,1.D0,T0XR(1,IT,IS1,IORB,1,IDOS),1,
     &                          CHIOR(1,IT,IS1),1)
C
                     DO JT = 1,NT
                        CALL DAXPY(IRTOP,1.D0,
     &                             TIJXR(1,IT,JT,IS1,IORB,1,IDOS),1,
     &                             CHIOR(1,IT,IS1),1)
                     END DO
                  ELSE
                     CALL DAXPY(IRTOP,CHIN(IT),
     &                          T0XR(1,IT,IS1,IORB,1,IDOS),1,
     &                          CHIOR(1,IT,IS1),1)
C
                     DO JT = 1,NT
                        CALL DAXPY(IRTOP,CHIN(JT),
     &                             TIJXR(1,IT,JT,IS1,IORB,1,IDOS),1,
     &                             CHIOR(1,IT,IS1),1)
                     END DO
                  END IF
C
C--- contributions due to \Delta E_FERMI
C
                  DO I = 1,IRTOP
                     WIO(IS1) = POBSR(I,0,IT,IS1,IORB)
                     CHIOR(I,IT,IS1) = CHIOR(I,IT,IS1) - WIO(IS1)
     &                                 *SUMWTZNS - WIO(IS1)
     &                                 *(SUMWTZNO(1)+SUMWTZNO(2))
                     DO JT = 1,NT
                        CHIOR(I,IT,IS1) = CHIOR(I,IT,IS1) - WIO(IS1)
     &                     *WTXNO(JT,1)*CHIO(JT,1) - WIO(IS1)
     &                     *WTXNO(JT,2)*CHIO(JT,2) - WIO(IS1)*WTXNS(JT)
     &                     *CHIS(JT)
C
                        IF ( NT.EQ.1 ) THEN
                           CHIOR(I,IT,IS1) = CHIOR(I,IT,IS1) - WIO(IS1)
     &                        *WTXNN(JT)
                        ELSE
                           CHIOR(I,IT,IS1) = CHIOR(I,IT,IS1) - WIO(IS1)
     &                        *WTXNN(JT)*CHIN(JT)
                        END IF
                     END DO
                  END DO
C
               END DO
C
C------------------------------------ orbital polarisation related terms
C
               CALL DCOPY(IRTOP,CHIOR(1,IT,1),1,CHIOR(1,IT,0),1)
               CALL DAXPY(IRTOP,1D0,CHIOR(1,IT,2),1,CHIOR(1,IT,0),1)
C
               DO I = 1,IRTOP
                  RINT(I) = TSSZRLOP(I,IT)*R2DRDI(I,IM)
               END DO
C
               CALL RRADINT(IM,RINT,NORM)
C
               DO I = 1,IRTOP
                  RSSZLOP(I,IT) = TSSZRLOP(I,IT)/NORM
               END DO
C
C-------------------------------------- radial charge perturbation CHINR
C--- NS contribution
C
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
               CALL DCOPY(IRTOP,TZR(1,IT,1,IDOS,1,ISPN),1,CHINR(1,IT),1)
               CALL DAXPY(IRTOP,CHIS(IT),T0XR(1,IT,1,IDOS,1,ISPN),1,
     &                    CHINR(1,IT),1)
C
               DO JT = 1,NT
                  CALL DAXPY(IRTOP,CHIS(JT),TIJXR(1,IT,JT,1,IDOS,1,ISPN)
     &                       ,1,CHINR(1,IT),1)
               END DO
C
C--- NO contribution
C
               DO IS = 1,2
C
                  CALL DAXPY(IRTOP,1.D0,TZR(1,IT,1,IDOS,IS,IORB),1,
     &                       CHINR(1,IT),1)
                  CALL DAXPY(IRTOP,CHIO(IT,IS),T0XR(1,IT,1,IDOS,IS,IORB)
     &                       ,1,CHINR(1,IT),1)
C
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIO(JT,IS),
     &                          TIJXR(1,IT,JT,1,IDOS,IS,IORB),1,
     &                          CHINR(1,IT),1)
                  END DO
C
               END DO
C
C--- NN contribution
C
               IF ( NT.EQ.1 ) THEN
                  CALL DAXPY(IRTOP,1.D0,T0XR(1,IT,1,IDOS,1,IDOS),1,
     &                       CHINR(1,IT),1)
C
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,1.D0,TIJXR(1,IT,JT,1,IDOS,1,IDOS),
     &                          1,CHINR(1,IT),1)
                  END DO
               ELSE
                  CALL DAXPY(IRTOP,CHIN(IT),T0XR(1,IT,1,IDOS,1,IDOS),1,
     &                       CHINR(1,IT),1)
C
                  DO JT = 1,NT
                     CALL DAXPY(IRTOP,CHIN(JT),
     &                          TIJXR(1,IT,JT,1,IDOS,1,IDOS),1,
     &                          CHINR(1,IT),1)
                  END DO
               END IF
C
C--- contributions due to \Delta E_FERMI
C
               DO I = 1,IRTOP
                  WIN = POBSR(I,0,IT,1,IDOS)
                  CHINR(I,IT) = CHINR(I,IT)
     &                          - WIN*(SUMWTZNO(1)+SUMWTZNO(2))
     &                          - WIN*SUMWTZNS
C
                  DO JT = 1,NT
                     CHINR(I,IT) = CHINR(I,IT) - WIN*WTXNO(JT,1)
     &                             *CHIO(JT,1) - WIN*WTXNO(JT,2)
     &                             *CHIO(JT,2) - WIN*WTXNS(JT)*CHIS(JT)
C
                     IF ( NT.EQ.1 ) THEN
                        CHINR(I,IT) = CHINR(I,IT) - WIN*WTXNN(JT)
                     ELSE
                        CHINR(I,IT) = CHINR(I,IT) - WIN*WTXNN(JT)
     &                                *CHIN(JT)
                     END IF
                  END DO
               END DO
C
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C----------------------- Calculation of the Hartree contribution
C
               DO IR = 1,IRTOP
                  RINT(IR) = CHINR(IR,IT)*R2DRDI(IR,IM)
               END DO
               CALL RRADINT_R(IM,RINT,VH_A)
C
               DO IR = 1,IRTOP
                  RINT(IR) = CHINR(IR,IT)*R(IR,IM)*DRDI(IR,IM)
               END DO
C
               CALL RRADINT_R(IM,RINT,VH_B)
C
               DO IR = 1,IRTOP
                  VH(IR) = VH_A(IR)/R(IR,IM) + VH_B(IRTOP) - VH_B(IR)
               END DO
C
C-----------------------------------------------------------------------
C
               TOL = 1D-8
C
C----------------------------------------------- consistency checks spin
C
               IF ( CHIPRINT.GE.1 ) THEN
                  IF ( ISOCOUPLING.EQ.1 ) THEN
                     CHISUM = CHISS(IT)
                  ELSE
                     CHISUM = CHISS(IT) + CHISO(IT)
                  END IF
C
                  IF ( ABS(CHIS(IT)-CHISUM).GT.TOL ) WRITE (6,99001)
     &                  'SUM','spin',IT,TXT_T(IT),CHIS(IT),CHISUM,
     &                 (1D0-CHIS(IT)/CHISUM)
C
                  CHISUM = 0D0
                  DO IL = 1,NL
                     CHISUM = CHISUM + CHISSL(IL,IT)
                  END DO
                  IF ( ABS(CHISS(IT)-CHISUM).GT.TOL ) WRITE (6,99001)
     &                  'L-sum','spin',IT,TXT_T(IT),CHISS(IT),CHISUM,
     &                 (1D0-CHISS(IT)/CHISUM)
C
                  DO I = 1,IRTOP
                     RINT(I) = CHISR(I,IT)*R2DRDI(I,IM)
                  END DO
                  CALL RRADINT(IM,RINT,CHISUM)
                  IF ( ISOCOUPLING.EQ.1 ) CHISUM = CHISUM - CHISO(IT)
                  IF ( ABS(CHIS(IT)-CHISUM).GT.TOL ) WRITE (6,99001)
     &                  'R-int','spin',IT,TXT_T(IT),CHIS(IT),CHISUM,
     &                 (1D0-CHIS(IT)/CHISUM)
C
C-------------------------------------------- consistency checks orbital
C
                  DO IS = 1,2
                     IF ( ISOCOUPLING.EQ.1 ) THEN
                        CHISUM = CHIOO(IT,IS)
                     ELSE
                        CHISUM = CHIOO(IT,IS) + CHIOS(IT,IS)
                     END IF
C
                     IF ( ABS(CHIO(IT,IS)-CHISUM).GT.TOL )
     &                    WRITE (6,99002) 'SUM','orbital',IT,TXT_T(IT),
     &                    IS,CHIO(IT,IS),CHISUM,(1D0-CHIO(IT,IS)/CHISUM)
C
                     CHISUM = 0D0
                     DO IL = 1,NL
                        CHISUM = CHISUM + CHIOOL(IL,IT,IS)
                     END DO
                     IF ( ABS(CHIOO(IT,IS)-CHISUM).GT.TOL )
     &                    WRITE (6,99002) 'L-sum','orbital  IT=',IT,
     &                    TXT_T(IT),IS,CHIOO(IT,IS),CHISUM,
     &                    (1D0-CHIOO(IT,IS)/CHISUM)
C
                     DO I = 1,IRTOP
                        RINT(I) = CHIOR(I,IT,IS)*R2DRDI(I,IM)
                     END DO
                     CALL RRADINT(IM,RINT,CHISUM)
                     IF ( ISOCOUPLING.EQ.1 ) CHISUM = CHISUM - 
     &                    CHIOS(IT,IS)
                     IF ( ABS(CHIO(IT,IS)-CHISUM).GT.TOL )
     &                    WRITE (6,99002) 'R-int','orbital',IT,TXT_T(IT)
     &                    ,IS,CHIO(IT,IS),CHISUM,
     &                    (1D0-CHIO(IT,IS)/CHISUM)
                  END DO
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
            IF ( CHIPRINT.GE.1 ) THEN
               WRITE (6,99008) (IT,T1Z(IT,1,ISPN,1,ISPN),T0Z(IT,1,ISPN,1
     &                         ,ISPN),IT=1,NT)
               WRITE (6,99009) (IT,T0X(IT,1,ISPN,1,ISPN),IT=1,NT)
C
               DO IT = 1,NT
                  WRITE (6,'('' IT='',I2)') IT
                  WRITE (6,99010) (JT,'Z',TIJZ(IT,JT,1,ISPN,1,ISPN),JT=1
     &                            ,NT)
                  WRITE (6,99010) (JT,'X',TIJX(IT,JT,1,ISPN,1,ISPN),JT=1
     &                            ,NT)
               END DO
C
               WRITE (6,99011) ((IT,IS1,T1Z(IT,IS1,IORB,IS1,IORB),T0Z(IT
     &                         ,IS1,IORB,IS1,IORB),IS1=1,2),IT=1,NT)
               WRITE (6,99012) (((IT,IS1,IS2,T0X(IT,IS1,IORB,IS2,IORB),
     &                         IS1=1,2),IS2=1,2),IT=1,NT)
C
               DO IT = 1,NT
                  WRITE (6,'('' IT='',I2)') IT
                  WRITE (6,99013) ((JT,IS1,TIJZ(IT,JT,IS1,IORB,IS1,IORB)
     &                            ,IS1=1,2),JT=1,NT)
                  WRITE (6,99014) (((JT,IS1,IS2,TIJX(IT,JT,IS1,IORB,IS2,
     &                            IORB),IS1=1,2),IS2=1,2),JT=1,NT)
               END DO
            END IF
C
C
C ======================================================================
C                    calculate Knight shift
C ======================================================================
C
            DO IT = 1,NT
C
C------------------------------------------------------ high field  HF-S
C
C
               WIHF = POBS(0,IT,1,IHFI)
               SFTS0(IT) = SFTS0(IT) + WIHF*SUMWTZNS
               SFTS(IT) = SFTS(IT) + WIHF*SUMWTZNS
               DO JT = 1,NT
                  SFTS(IT) = SFTS(IT) + WIHF*WTXNS(JT)*CHIS(JT)
               END DO
C
               DO IL = 1,NLT(IT)
                  WIHF = POBS(IL,IT,1,IHFI)
                  SFTSL0(IL,IT) = SFTSL0(IL,IT) + WIHF*SUMWTZNS
                  SFTSL(IL,IT) = SFTSL(IL,IT) + WIHF*SUMWTZNS
                  DO JT = 1,NT
                     SFTSL(IL,IT) = SFTSL(IL,IT) + WIHF*WTXNS(JT)
     &                              *CHIS(JT)
                  END DO
               END DO
C
C------------------------------------------------------ high field  HF-O
C
               WIHF = POBS(0,IT,1,IHFI)
               DO IS2 = 1,2
                  SFTO0(IT) = SFTO0(IT) + WIHF*SUMWTZNO(IS2)
                  SFTO(IT) = SFTO(IT) + WIHF*SUMWTZNO(IS2)
                  DO JT = 1,NT
                     SFTO(IT) = SFTO(IT) + WIHF*WTXNO(JT,IS2)
     &                          *CHIO(JT,IS2)
                  END DO
               END DO
C
               DO IL = 1,NLT(IT)
                  WIHF = POBS(IL,IT,1,IHFI)
C
                  DO IS2 = 1,2
                     SFTOL0(IL,IT) = SFTOL0(IL,IT) + WIHF*SUMWTZNO(IS2)
                     SFTOL(IL,IT) = SFTOL(IL,IT) + WIHF*SUMWTZNO(IS2)
                     DO JT = 1,NT
                        SFTOL(IL,IT) = SFTOL(IL,IT) + WIHF*WTXNO(JT,IS2)
     &                                 *CHIO(JT,IS2)
                     END DO
                  END DO
               END DO
C
            END DO
C
C ======================================================================
C            print out of SUSCEPTIBILITY and KNIGHT SHIFT
C                         including  HF-terms
C ======================================================================
C
            CALL CHIOUTPUT(LTXT_T,TXT_T,DATSET,LDATSET,BCP,BCPS,CHIO,
     &                     CHIOO,CHIOO0,CHIOOL,CHIOOLTR,CHIOOTR,CHIOR,
     &                     CHIOS,CHIOS0,CHIPRINT,CHIS,CHISO,CHISO0,
     &                     CHISR,CHISS,CHISS0,CHISSL,CHIN,CHINN,CHINN0,
     &                     CHINNL,CHILANDAU,CHILANLT,CHISSLTR,CHISSTR,
     &                     CONC,DOBSEF,ICHIO,ICHIS,IDOS,IHFI,IHVV,IMT,
     &                     IORB,ISPN,JRWS,NAT,NL,NLT,NSPINOBS,NT,NONMAG,
     &                     R,SFTO,SFTO0,SFTOL,SFTOL0,SFTOLTR,SFTOTR,
     &                     SFTS,SFTS0,SFTSL,SFTSL0,SFTSLTR,SFTSTR,
     &                     TCDIACRI,TCDIARI,TKDIACRI,TKDIARI,TOBS,TXTS,
     &                     TZ,TZL,TZR,ICDIA,IKDIA,IRM3,NOP,NOBS,NOBSDNS,
     &                     NLMAX,NMMAX,NRMAX,NTMAX,IOTMP)
C
         END IF
C                                      for FERRO all terms sumed up now
C ======================================================================
C
C
         IF ( .NOT.NONMAG ) WRITE (6,99022) DELTAEF
C
C ======================================================================
C                 ANTI-FERRO MAGNETIC COUPLING
C ======================================================================
C -------------------------- calculate atom type-resolved susceptibility
C
         IF ( TREATAF ) THEN
C
            CALL RINIT(NTMAX*4*NTMAX*4,AA)
            CALL RINIT(NTMAX*4,BB)
            CALL RINIT(NTMAX,CHIAFSS)
            CALL RINIT(NTMAX,CHIAFSS0)
            DO IT = 1,NT
               AA(IT,IT) = 1D0 - T0X(IT,1,ISPN,1,ISPN)*SIGNAF(IT)
               BB(IT,1) = T0Z(IT,1,ISPN,1,ISPN)*SIGNAF(IT)
               DO JT = 1,NT
                  AA(IT,JT) = AA(IT,JT) - TIJX(IT,JT,1,ISPN,1,ISPN)
     &                        *SIGNAF(JT)
                  BB(IT,1) = BB(IT,1) + TIJZ(IT,JT,1,ISPN,1,ISPN)
     &                       *SIGNAF(JT)
               END DO
            END DO
            DO IT = 1,NT
               CHIAFSS0(IT) = BB(IT,1)
            END DO
            CALL DGESV(NT,1,AA,NTMAX,IPIV,BB,NTMAX,INFO)
            IF ( INFO.NE.0 ) THEN
               WRITE (6,*) ' STOP IN < CHICALC > AFTER ZGESV FOR SPN'
               WRITE (6,*) ' INFO = ',INFO
               STOP
            END IF
            DO IT = 1,NT
               IM = IMT(IT)
               IRTOP = JRWS(IM)
               CHIAFSS(IT) = BB(IT,1)
               ICHIAFSS0(IT) = (CHIAFSS(IT)-CHIAFSS0(IT))/CHIAFSS(IT)
            END DO
            WRITE (6,99006)
            WRITE (6,99007)
            WRITE (6,99017) (IT,NINT(SIGNAF(IT)),CHIAFSS0(IT)*CU,
     &                      ICHIAFSS0(IT),CHIAFSS(IT)*CU,IT=1,NT)
C
         END IF
C
      END DO
C                                                               S-O-LOOP
C=======================================================================
C
      IF ( NT.GT.100 ) STOP
C
C ======================================================================
C                          update AXCN and BXCN
C ======================================================================
C
C >> take care of numerical inaccuracy for small r
C
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         INSIDE = .TRUE.
         DO I = 1,IRTOP
            IF ( INSIDE .AND. (CHISR(I,IT).LT.0.0D0) ) THEN
               CHISR(I,IT) = 0.0D0
            ELSE
               INSIDE = .FALSE.
            END IF
            RINT(I) = CHISR(I,IT)*R2DRDI(I,IM)
         END DO
C
         CALL RRADINT(IM,RINT,CHISSRI)
C
         IF ( ABS(1D0-CHISSRI/CHIS(IT)).GT.TOLRADINT ) WRITE (6,99003)
     &         'CHISS',IT,CHISSRI/CHIS(IT)
C
         DO I = 1,IRTOP
C
            GAMMAS(I,IT) = CHISR(I,IT)/CHIS(IT)
            IF ( NT.EQ.1 ) THEN
               GAMMAN(I,IT) = CHINR(I,IT)
            ELSE
               GAMMAN(I,IT) = CHINR(I,IT)/CHIN(IT)
            END IF
C
            BXCNMM(I,IT) = GAMMAS(I,IT)*EMM(I,IT)/(4D0*PI)
            BXCNMN(I,IT) = GAMMAS(I,IT)*ENM(I,IT)/(4D0*PI)
            BXCNNM(I,IT) = GAMMAN(I,IT)*ENM(I,IT)/(4D0*PI)
            BXCNNN(I,IT) = GAMMAN(I,IT)*ENN(I,IT)/(4D0*PI)
C
            DO IS = 1,2
               GAMMAOP(I,IS) = RSSZLOP(I,IT)
               RINTGAMMAOP(I,IS) = GAMMAOP(I,IS)*R2DRDI(I,IM)
               RHOLOPS(I,IS) = RHOLOP(I,IT)
            END DO
C
         END DO
C
         CALL RRADINT(IM,RINTGAMMAOP(1,1),GAMMAOPINT(1))
         CALL RRADINT(IM,RINTGAMMAOP(1,2),GAMMAOPINT(2))
C
         CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_GAMMA-OP.dat',
     &                 13,FILNAM,LFN,2,IOTMP,'gammas(OP) ',10,NTMAX)
         WRITE (IOTMP,99015)
         WRITE (IOTMP,'(4(1E15.7,1X))') (R(I,IM),GAMMAOP(I,1),GAMMAOP(I,
     &                                  2),RHOLOP(I,IT),I=1,IRTOP)
         CLOSE (IOTMP)
C
         CALL CHIORBPOL(AXCN(1,1,1,IT),RHOLOPS,ORBSQINT,R2DRDI(1,IM),
     &                  RPWOP,Z(IT),ORBPOL,TXT_T(IT),IT,IRTOP,NLMAX,
     &                  NRMAX)
      END DO
C
      WRITE (6,99019) CU,KU
C
C --------------------------------------- calculate magnetic form factor
C
      CALL CHIFMAGIND(RHOCHR,CHISR,CHIOR)
C
C ======================================================================
C
99001 FORMAT (/,' ####### WARNING <CHICALC>: ',A,' inconsistency for ',
     &        A,'  IT=',I2,2X,A,:,/,' ####### WARNING ',3F12.6)
99002 FORMAT (/,' ####### WARNING <CHICALC>: ',A,' inconsistency for ',
     &        A,'  IT=',I2,2X,A,' IS=',I2,:,/,' ####### WARNING ',
     &        3F12.6)
99003 FORMAT (/,' ####### WARNING <CHICALC>:  inconsistency for ',A,
     &        '  (rad.integ.)  IT=',I2,F15.10)
99004 FORMAT ('# check the sum rule dG/dE = - G * G ',/,
     &        '# for type  IT=',I2,4X,A,'   and   NL=',I3,/,
     &        '#            DOS -INT(G*G) dE      .........',/,
     &        '#   E        tot      tot      s        s        p      '
     &        ,'  p        d        d')
99005 FORMAT ('# radial functions ',/,'# for type  IT=',I2,4X,A,
     &        '   IL=',I2,' == ',A,/,
     &        '#    R                 RHORINTL        ',
     &        'TSSZRL          TLLZRL          TSSXRL          ',
     &        'THFSZRL         THFSXR')
99006 FORMAT (/,10X,'atom type resolved spin susceptibility',
     &        '  [10^(-6) cm^3/mol]')
99007 FORMAT (16X,'for an anti-ferromagnetic spin configuration',/)
99008 FORMAT (/,8(:,' IT=',I2,' T1SSZ =',1E13.5,' T0SSZ =',1E13.5,/))
99009 FORMAT (8(:,' IT=',I2,' T0SSX  =',1E13.5,/))
99010 FORMAT (8(:,' JT=',I2,' TIJSS',A,' =',1E13.5,/))
99011 FORMAT (/,16(:,' IT=',I2,' IS=',I2,' T1LLZS =',1E13.5,' T0LLZS =',
     &        1E13.5,/))
99012 FORMAT (32(:,' IT=',I2,' IS1=',I2,' IS2=',I2,' T0LLXS  =',1E13.5,
     &        /))
99013 FORMAT (32(:,' JT=',I2,' IS1=',I2,' TIJLLZS =',1E13.5,/))
99014 FORMAT (32(:,' JT=',I2,' IS1=',I2,' IS2=',I2,' TIJLLXS =',1E13.5,
     &        /))
99015 FORMAT ('# FILE fort.(700+IT)',/,
     &        '# LAST: R, GAMMAOP(1), GAMMAOP(2), RHOLOPS')
99016 FORMAT (I4,' E=',2F7.4,3X,'IT=',I2,2X,A2,2X,A20,/,16X,
     &        'DOS  [1/Ry]     |','     GGINT  [1/Ry]     |',
     &        '     DOSINT',/,' SUM(L)=',2E12.4,2E12.4,2E12.4,
     &        5(:,/,'  L=',I2,': ',2E12.4,2E12.4,2E12.4))
99017 FORMAT (14X,'SIGN     CHIS0          ICHIS              CHIS ',
     &        20(:,/,5X,'IT=',I2,4X,I3,2E16.7,2X,1E16.7))
99018 FORMAT (8E10.4,:,/,(10X,7E10.4))
99019 FORMAT (/,5X,'conversion factors   a.u. -->  cgs ',/5X,
     &        'susceptibility ',E15.6,'  10^(-6) cm^3/mol / a.u.',/,5X,
     &        'Knight shift   ',E15.6,'  % / a.u.',/)
99020 FORMAT (//,1X,79('*'),/,1X,79('*'),/,10X,
     &       'exchange coupling of spin and orbital susceptibilities:  '
     &       ,A,/,1X,79('*'),/,1X,79('*'),/)
99021 FORMAT (//,1X,79('*'),/,1X,16X,'results  ',A,
     &        '  additional high-field terms',/,1X,79('*'),/)
99022 FORMAT (10X,'norm. change of Fermi energy ',F15.6,' a.u.',//,1X,
     &        79('='),/)
      END
C*==chioutput.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIOUTPUT(LTXT_T,TXT_T,DATSET,LDATSET,BCP,BCPS,CHIO,
     &                     CHIOO,CHIOO0,CHIOOL,CHIOOLTR,CHIOOTR,CHIOR,
     &                     CHIOS,CHIOS0,CHIPRINT,CHIS,CHISO,CHISO0,
     &                     CHISR,CHISS,CHISS0,CHISSL,CHIN,CHINN,CHINN0,
     &                     CHINNL,CHILANDAU,CHILANLT,CHISSLTR,CHISSTR,
     &                     CONC,DOBSEF,ICHIO,ICHIS,IDOS,IHFI,IHVV,IMT,
     &                     IORB,ISPN,JRWS,NAT,NL,NLT,NSPINOBS,NT,NONMAG,
     &                     R,SFTO,SFTO0,SFTOL,SFTOL0,SFTOLTR,SFTOTR,
     &                     SFTS,SFTS0,SFTSL,SFTSL0,SFTSLTR,SFTSTR,
     &                     TCDIACRI,TCDIARI,TKDIACRI,TKDIARI,TOBS,TXTS,
     &                     TZ,TZL,TZR,ICDIA,IKDIA,IRM3,NOP,NOBS,NOBSDNS,
     &                     NLMAX,NMMAX,NRMAX,NTMAX,IOTMP)
C   ********************************************************************
C   *                                                                  *
C   *   print results for susceptibility and Knight shift              *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1        1 IDOS                           *
C   *  spin moment           b s_z    2 ISPN                           *
C   *  orbital moment        b l_z    3 IORB   ____ NOP                *
C   *  hyperfine interaction H_hf,z   4 IHFI                           *
C   *  VanVlecK hyperfine    H_VV,z   5 IHVV   ____ NOBS = NOP+2       *
C   *  diamagnetic suscept   r^2      6 ICDIA                          *
C   *  diamagnetic hyperfine r^-1     7 IKDIA                          *
C   *  hyperfine r-weight    r^-3     8 IRM3   ____ NOBSDNS = NOBS+3   *
C   *                                                                  *
C   *  SM, HE  2003                                                    *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:RY_ERG,MB_CGS,ALPHA_FS,CHI_AU2CGS
      USE MOD_ANGMOM,ONLY:TXT_L
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--CHIOUTPUT2515
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHIOUTPUT')
      REAL*8 CU,KU,FKVV,FKDIA
      PARAMETER (CU=CHI_AU2CGS*1D+6,KU=100D0*MB_CGS/RY_ERG,
     &           FKVV=100D0*ALPHA_FS*ALPHA_FS,FKDIA=FKVV/3D0)
C
C Dummy arguments
C
      LOGICAL CHILANDAU,NONMAG
      INTEGER CHIPRINT,ICDIA,IDOS,IHFI,IHVV,IKDIA,IORB,IOTMP,IRM3,ISPN,
     &        LDATSET,NL,NLMAX,NMMAX,NOBS,NOBSDNS,NOP,NRMAX,NT,NTMAX
      CHARACTER*80 DATSET
      REAL*8 BCP(NTMAX),BCPS(NTMAX),CHILANLT(0:NLMAX,NTMAX,2),
     &       CHIN(NTMAX),CHINN(NTMAX),CHINN0(NTMAX),CHINNL(NLMAX,NTMAX),
     &       CHIO(NTMAX,0:2),CHIOO(NTMAX,0:2),CHIOO0(NTMAX,0:2),
     &       CHIOOL(NLMAX,NTMAX,0:2),CHIOOLTR(NLMAX,NTMAX,NTMAX,0:2),
     &       CHIOOTR(NTMAX,NTMAX,0:2),CHIOR(NRMAX,NTMAX,0:2),
     &       CHIOS(NTMAX,0:2),CHIOS0(NTMAX,0:2),CHIS(NTMAX),CHISO(NTMAX)
     &       ,CHISO0(NTMAX),CHISR(NRMAX,NTMAX),CHISS(NTMAX),
     &       CHISS0(NTMAX),CHISSL(NLMAX,NTMAX),
     &       CHISSLTR(NLMAX,NTMAX,NTMAX),CHISSTR(NTMAX,NTMAX),
     &       CONC(NTMAX),DOBSEF(NLMAX,NTMAX,2,NOBSDNS),ICHIO(NTMAX),
     &       ICHIS(NTMAX),R(NRMAX,NMMAX),SFTO(NTMAX),SFTO0(NTMAX),
     &       SFTOL(NLMAX,NTMAX),SFTOL0(NLMAX,NTMAX),
     &       SFTOLTR(NLMAX,NTMAX,NTMAX),SFTOTR(NTMAX,NTMAX),SFTS(NTMAX),
     &       SFTS0(NTMAX),SFTSL(NLMAX,NTMAX),SFTSL0(NLMAX,NTMAX),
     &       SFTSLTR(NLMAX,NTMAX,NTMAX),SFTSTR(NTMAX,NTMAX),
     &       TCDIACRI(NTMAX),TCDIARI(NTMAX),TKDIACRI(NTMAX),
     &       TKDIARI(NTMAX),TOBS(NLMAX,NTMAX,2,NOBSDNS),
     &       TZ(NTMAX,2,NOBS,2,NOP),TZL(NLMAX,NTMAX,2,NOBS,2,NOP),
     &       TZR(NRMAX,NTMAX,2,NOBS,2,NOP)
      INTEGER IMT(NTMAX),JRWS(NMMAX),LTXT_T(NTMAX),NAT(NTMAX),NLT(NTMAX)
     &        ,NSPINOBS(NOBS)
      CHARACTER*2 TXTS(2)
      CHARACTER*8 TXT_T(NTMAX)
C
C Local variables
C
      REAL*8 CHISUM,CHITOT,DSUM(2,NOBS),RSUM,SFT0
      CHARACTER*80 FILNAM
      INTEGER I,IL,IM,IOBS,IRTOP,IS,IS1,IT,JT,LFN
C
C*** End of declarations rewritten by SPAG
C
      CHITOT = 0D0
C
      DO IT = 1,NT
C
         WRITE (6,99007) IT,TXT_T(IT)
C
C----------------------------------------------- polarisation at E_Fermi
C
         IF ( .NOT.NONMAG .OR. CHIPRINT.GT.1 ) THEN
C
            DO IOBS = 1,IORB
               DO IS1 = 1,NSPINOBS(IOBS)
                  DSUM(IS1,IOBS) = 0D0
                  DO IL = 1,NLT(IT)
                     DSUM(IS1,IOBS) = DSUM(IS1,IOBS)
     &                                + DOBSEF(IL,IT,IS1,IOBS)
                  END DO
               END DO
            END DO
C
            WRITE (6,99008) ((DSUM(IS1,IOBS),IS1=1,NSPINOBS(IOBS)),
     &                      IOBS=1,IORB),
     &                      (TXT_L(IL-1),((DOBSEF(IL,IT,IS1,IOBS),IS1=1,
     &                      NSPINOBS(IOBS)),IOBS=1,IORB),IL=1,NLT(IT))
         END IF
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CHISUM = CHIS(IT) + CHIO(IT,0) + TCDIARI(IT)/3D0
C
         IF ( CHILANDAU ) CHISUM = CHISUM + CHILANLT(0,IT,1)
     &                             + CHILANLT(0,IT,2)
C
         CHITOT = CHITOT + NAT(IT)*CONC(IT)*CHISUM
C
         WRITE (6,99009) CHISUM*CU
C
         WRITE (6,99011) 'CHI S      :',CHIS(IT)*CU
C
         WRITE (6,99011) 'CHI SS     :',CHISS(IT)*CU,CHISS0(IT)*CU,
     &                   (CHISS(IT)-CHISS0(IT))*CU
C
         WRITE (6,99012) (TXT_L(IL-1),CHISSL(IL,IT)*CU,TZL(IL,IT,1,ISPN,
     &                   1,ISPN)*CU,
     &                   (CHISSL(IL,IT)-TZL(IL,IT,1,ISPN,1,ISPN))*CU,
     &                   IL=1,NLT(IT))
C
         WRITE (6,99011) 'ICHI(S)    :',ICHIS(IT)
C
         WRITE (6,99011) 'CHI SO     :',CHISO(IT)*CU,CHISO0(IT)*CU,
     &                   (CHISO(IT)-CHISO0(IT))*CU
C
         WRITE (6,99011) 'CHI O      :',CHIO(IT,0)*CU
C
         WRITE (6,99011) 'CHI OO     :',CHIOO(IT,0)*CU,CHIOO0(IT,0)*CU,
     &                   (CHIOO(IT,0)-CHIOO0(IT,0))*CU
         WRITE (6,99012) (TXT_L(IL-1),CHIOOL(IL,IT,0)*CU,(TZL(IL,IT,1,
     &                   IORB,1,IORB)+TZL(IL,IT,1,IORB,2,IORB)
     &                   +TZL(IL,IT,2,IORB,1,IORB)
     &                   +TZL(IL,IT,2,IORB,2,IORB))*CU,
     &                   (CHIOOL(IL,IT,0)-TZL(IL,IT,1,IORB,1,IORB)
     &                   -TZL(IL,IT,1,IORB,2,IORB)
     &                   -TZL(IL,IT,2,IORB,1,IORB)
     &                   -TZL(IL,IT,2,IORB,2,IORB))*CU,IL=1,NLT(IT))
C
         IF ( .NOT.NONMAG .OR. CHIPRINT.GT.1 ) THEN
            DO IS = 1,2
               WRITE (6,99013) TXTS(IS),
     &                         (TXT_L(IL-1),CHIOOL(IL,IT,IS)*CU,
     &                         (TZL(IL,IT,IS,IORB,1,IORB)
     &                         +TZL(IL,IT,IS,IORB,2,IORB))*CU,
     &                         (CHIOOL(IL,IT,IS)
     &                         -TZL(IL,IT,IS,IORB,1,IORB)
     &                         -TZL(IL,IT,IS,IORB,2,IORB))*CU,IL=1,
     &                         NLT(IT))
            END DO
         END IF
C
         WRITE (6,99011) 'ICHI(O)    :',ICHIO(IT)
C
         WRITE (6,99011) 'CHI OS     :',CHIOS(IT,0)*CU,CHIOS0(IT,0)*CU,
     &                   (CHIOS(IT,0)-CHIOS0(IT,0))*CU
C
         IF ( CHILANDAU ) THEN
            WRITE (6,99017) 'CHI LAN    :',
     &                      (CHILANLT(0,IT,1)+CHILANLT(0,IT,2))*CU
            WRITE (6,99014) (TXT_L(IL-1),(CHILANLT(0,IT,1)+CHILANLT(0,IT
     &                      ,2))*CU,IL=1,NLT(IT))
         END IF
C
         WRITE (6,99016) 'CHI DIA    :',TCDIARI(IT)*CU/3D0,'core:',
     &                   TCDIACRI(IT)*CU/3D0,'  cb:',
     &                   (TCDIARI(IT)-TCDIACRI(IT))*CU/3D0
         WRITE (6,99014) (TXT_L(IL-1),TOBS(IL,IT,1,ICDIA)*CU/3D0,IL=1,
     &                   NLT(IT))
C
         WRITE (6,99011) 'CHI N     :',CHIN(IT)*CU
         WRITE (6,99011) 'CHI NN     :',CHINN(IT)*CU,CHINN0(IT)*CU,
     &                   (CHINN(IT)-CHINN0(IT))*CU
C
C    WRITE (6,99012) (TXT_L(IL-1),CHINNL(IL,IT)*CU,TZL(IL,IT,1,IDOS,1,
C    &                   IDOS)*CU,
C    &                   (CHINNL(IL,IT)-TZL(IL,IT,1,IDOS,1,IDOS))*CU,
C    &                   IL=1,NLT(IT))
Cc   WRITE (6,99012) (TXT_L(IL-1),CHINNL(IL,IT)*CU,TZL(IL,IT,1,IDOS,1,
Cc    &                   ISPN)*CU,
Cc    &                   (CHINNL(IL,IT)-TZL(IL,IT,1,IDOS,1,ISPN))*CU,
Cc    &                   IL=1,NLT(IT))
C
         WRITE (6,99012) (TXT_L(IL-1),CHINNL(IL,IT)*CU,0.,CHINNL(IL,IT)
     &                   *CU,IL=1,NLT(IT))
C
C
C -------------------------------------------------- OUTPUT KNIGHT SHIFT
C
         SFT0 = SFTS(IT) + SFTO(IT) + BCP(IT)*CHIS(IT)
C
         WRITE (6,99010) SFT0*KU + TKDIARI(IT)*FKDIA,SFT0*KU
C
         WRITE (6,99011) 'K HF-S  cb :',SFTS(IT)*KU,SFTS0(IT)*KU,
     &                   (SFTS(IT)-SFTS0(IT))*KU
C
         WRITE (6,99012) (TXT_L(IL-1),SFTSL(IL,IT)*KU,SFTSL0(IL,IT)*KU,
     &                   (SFTSL(IL,IT)-SFTSL0(IL,IT))*KU,IL=1,NLT(IT))
C
         WRITE (6,99016) 'K HF-S  CP :',BCP(IT)*CHIS(IT)*KU,'s-el',
     &                   BCPS(IT)*CHIS(IT)*KU,'non-s',(BCP(IT)-BCPS(IT))
     &                   *CHIS(IT)*KU
C
         WRITE (6,99011) 'K HF-O  cb :',SFTO(IT)*KU,SFTO0(IT)*KU,
     &                   (SFTO(IT)-SFTO0(IT))*KU
C
         WRITE (6,99012) (TXT_L(IL-1),SFTOL(IL,IT)*KU,SFTOL0(IL,IT)*KU,
     &                   (SFTOL(IL,IT)-SFTOL0(IL,IT))*KU,IL=1,NLT(IT))
C
         WRITE (6,99016) 'K DIA      :',TKDIARI(IT)*FKDIA,' core',
     &                   TKDIACRI(IT)*FKDIA,'   cb',
     &                   (TKDIARI(IT)-TKDIACRI(IT))*FKDIA
         WRITE (6,99014) (TXT_L(IL-1),TOBS(IL,IT,1,IKDIA)*FKDIA,IL=1,
     &                   NLT(IT))
C
         WRITE (6,*)
         WRITE (6,99016) 'BCP [MG]   :',BCP(IT)/1D6,' s-el',BCPS(IT)
     &                   /1D6,'non-s',(BCP(IT)-BCPS(IT))/1D6
C
         RSUM = 0D0
         DO IL = 1,NLT(IT)
            RSUM = RSUM + CHIOOL(IL,IT,0)*TOBS(IL,IT,1,IRM3)
     &             /TOBS(IL,IT,1,IDOS)
         END DO
C
         WRITE (6,99016) 'K VV       :',
     &                   (TZ(IT,1,IHVV,1,IORB)+TZ(IT,1,IHVV,2,IORB))
     &                   *FKVV,'   AP',RSUM*FKVV,'1/r^3'
         WRITE (6,99015) (TXT_L(IL-1),(TZL(IL,IT,1,IHVV,1,IORB)+TZL(IL,
     &                   IT,1,IHVV,2,IORB))*FKVV,FKVV*CHIOOL(IL,IT,0)
     &                   *TOBS(IL,IT,1,IRM3)/TOBS(IL,IT,1,IDOS),
     &                   TOBS(IL,IT,1,IRM3)/TOBS(IL,IT,1,IDOS),IL=1,
     &                   NLT(IT))
C
C-----------------------------------------------------------------------
C                                                 transferred quantities
C
         IF ( NT.GT.1 ) THEN
C
            WRITE (6,'(/,1X,79(''-''))')
C
            WRITE (6,99002) 'spin susceptibility CHISS',
     &                      (TXT_L(IL-1),IL=1,NLT(IT))
            WRITE (6,99003) CHISS(IT)*CU,(CHISSL(IL,IT)*CU,IL=1,NLT(IT))
            DO JT = 1,NT
               WRITE (6,99004) JT,CHISSTR(IT,JT)*CU,
     &                         (CHISSLTR(IL,IT,JT)*CU,IL=1,NLT(IT))
            END DO
C
            WRITE (6,99002) 'orbital susceptibility CHIOO',
     &                      (TXT_L(IL-1),IL=1,NLT(IT))
            WRITE (6,99003) CHIOO(IT,0)*CU,
     &                      (CHIOOL(IL,IT,0)*CU,IL=1,NLT(IT))
            DO JT = 1,NT
               WRITE (6,99004) JT,CHIOOTR(IT,JT,0)*CU,
     &                         (CHIOOLTR(IL,IT,JT,0)*CU,IL=1,NLT(IT))
               IF ( .NOT.NONMAG ) THEN
                  DO IS = 1,2
                     WRITE (6,99005) TXTS(IS),CHIOOTR(IT,JT,IS)*CU,
     &                               (CHIOOLTR(IL,IT,JT,IS)*CU,IL=1,
     &                               NLT(IT))
                  END DO
               END IF
            END DO
C
            WRITE (6,99002) 'Knight shift K HF-S',
     &                      (TXT_L(IL-1),IL=1,NLT(IT))
            WRITE (6,99003) SFTS(IT)*KU,(SFTSL(IL,IT)*KU,IL=1,NLT(IT))
            DO JT = 1,NT
               WRITE (6,99004) JT,SFTSTR(IT,JT)*KU,
     &                         (SFTSLTR(IL,IT,JT)*KU,IL=1,NLT(IT))
            END DO
C
            WRITE (6,99002) 'Knight shift K HF-O',
     &                      (TXT_L(IL-1),IL=1,NLT(IT))
            WRITE (6,99003) SFTO(IT)*KU,(SFTOL(IL,IT)*KU,IL=1,NLT(IT))
            DO JT = 1,NT
               WRITE (6,99004) JT,SFTOTR(IT,JT)*KU,
     &                         (SFTOLTR(IL,IT,JT)*KU,IL=1,NLT(IT))
            END DO
C
         END IF
C
         WRITE (6,'(/,1X,79(''-''))')
C
         CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_CHI-r.dat',10,
     &                 FILNAM,LFN,2,IOTMP,'CHI(r)file',10,NTMAX)
         WRITE (IOTMP,99006) IT,TXT_T(IT)(1:LTXT_T(IT)),NL
C
C         do I=1,IRTOP
C             WRITE (34,'(4F14.6)') R(I,IM),CHISR(I,IT)/CHISS(IT)
C         end do
         WRITE (IOTMP,'(13E14.6)') (R(I,IM),CHISR(I,IT),TZR(I,IT,1,ISPN,
     &                             1,ISPN),CHIOR(I,IT,1),CHIOR(I,IT,2),
     &                             (TZR(I,IT,1,IORB,1,IORB)
     &                             +TZR(I,IT,1,IORB,2,IORB)),
     &                             (TZR(I,IT,2,IORB,1,IORB)
     &                             +TZR(I,IT,2,IORB,2,IORB)),
     &                             TZR(I,IT,1,IHFI,1,ISPN),I=1,IRTOP)
         CLOSE (IOTMP)
C
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) THEN
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI S',IT,CHIS(IT)*CU,CHISS(IT)
     &                                 *CU,CHISS0(IT)*CU,
     &                                 (CHISSL(IL,IT)*CU,
     &                                 TZL(IL,IT,1,ISPN,1,ISPN)*CU,IL=1,
     &                                 NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'ICHI(S)',IT,ICHIS(IT)
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI SO',IT,CHISO(IT)*CU,
     &                                 CHISO0(IT)*CU
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI O',IT,CHIO(IT,0)*CU,
     &                                 CHIOO(IT,0)*CU,CHIOO0(IT,0)*CU,
     &                                 (CHIOOL(IL,IT,0)*CU,
     &                                 (TZL(IL,IT,1,IORB,1,IORB)
     &                                 +TZL(IL,IT,1,IORB,2,IORB)
     &                                 +TZL(IL,IT,2,IORB,1,IORB)
     &                                 +TZL(IL,IT,2,IORB,2,IORB))*CU,
     &                                 IL=1,NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'ICHI(O)',IT,ICHIO(IT)
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI OS',IT,CHIOS(IT,0)*CU,
     &                                 CHIOS0(IT,0)*CU
C
            IF ( CHILANDAU ) WRITE (IFILBUILDBOT,99018)
     &                              ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'CHI LAN',IT,
     &                              (CHILANLT(0,IT,1)+CHILANLT(0,IT,2))
     &                              *CU,
     &                              ((CHILANLT(0,IT,1)+CHILANLT(0,IT,2))
     &                              *CU,IL=1,NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI DIA',IT,TCDIARI(IT)*CU/3D0,
     &                                 TCDIACRI(IT)*CU/3D0,
     &                                 (TOBS(IL,IT,1,ICDIA)*CU/3D0,IL=1,
     &                                 NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'CHI N',IT,CHIN(IT)*CU,CHINN(IT)
     &                                 *CU,CHINN0(IT)*CU
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'K HF-S  cb',IT,SFTS(IT)*KU,
     &                                 SFTS0(IT)*KU,
     &                                 (SFTSL(IL,IT)*KU,SFTSL0(IL,IT)
     &                                 *KU,(SFTSL(IL,IT)-SFTSL0(IL,IT))
     &                                 *KU,IL=1,NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'K HF-S  CP',IT,BCP(IT)*CHIS(IT)
     &                                 *KU,BCPS(IT)*CHIS(IT)*KU
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'K HF-O  cb',IT,SFTO(IT)*KU,
     &                                 SFTO0(IT)*KU,
     &                                 (SFTOL(IL,IT)*KU,SFTOL0(IL,IT)
     &                                 *KU,IL=1,NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'K DIA',IT,TKDIARI(IT)*FKDIA,
     &                                 TKDIACRI(IT)*FKDIA,
     &                                 (TOBS(IL,IT,1,IKDIA)*FKDIA,IL=1,
     &                                 NLT(IT))
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'BCP [MG]',IT,BCP(IT)/1D6,
     &                                 BCPS(IT)/1D6
C
            WRITE (IFILBUILDBOT,99018) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'K VV',IT,
     &                                 (TZ(IT,1,IHVV,1,IORB)+TZ(IT,1,
     &                                 IHVV,2,IORB))*FKVV,RSUM*FKVV,
     &                                 ((TZL(IL,IT,1,IHVV,1,IORB)
     &                                 +TZL(IL,IT,1,IHVV,2,IORB))*FKVV,
     &                                 FKVV*CHIOOL(IL,IT,0)
     &                                 *TOBS(IL,IT,1,IRM3)
     &                                 /TOBS(IL,IT,1,IDOS),
     &                                 TOBS(IL,IT,1,IRM3)
     &                                 /TOBS(IL,IT,1,IDOS),IL=1,NLT(IT))
C
         END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C-----------------------------------------------------------------------
      END DO
C
      IF ( NT.EQ.1 ) THEN
         WRITE (6,'(/,1X,79(''=''))')
      ELSE
         WRITE (6,99001) CHITOT*CU
      END IF
C
99001 FORMAT (/,1X,79('='),//,10X,'total magnetic susceptibility',F15.6,
     &        ' x 10^(-6) cm^3/mol',//,1X,79('='),/)
99002 FORMAT (/,10X,'transferred ',A,//,23X,'sum         ',5(A,:,10X))
99003 FORMAT (10X,'total  ',6F11.4)
99004 FORMAT (10X,'JT =',I3,6F15.4)
99005 FORMAT (15X,A,6F11.4)
99006 FORMAT ('# for type  IT=',I2,4X,A,'   and   NL=',I3,/,
     &        '# R             CHISSR        TSSZR         TSSXR ',
     &        ' TLLZR         CHILLSR(1)    CHILLSR(2)    TLLZSR(1)   ',
     &        '  TLLXSR(2)     TLLZSR(1)*SCLOH2S(1) ..(2)  THFSZRL    ',
     &        '   THFSXR*CHISS')
99007 FORMAT (/,1X,79('='),//,10X,'results for atom type  IT=',I2,2X,A)
99008 FORMAT (//,10X,'polarisations at the Fermi level',//,30X,
     &        'DOS         P_spin    P_orb(dn)   P_orb(up)',/,22X,F15.6,
     &        3F12.6/,(21X,A,F15.6,3F12.6))
99009 FORMAT (//,10X,'magnetic susceptibility',//,10X,'CHI tot    :',
     &        F15.6,' x 10^(-6) cm^3/mol')
99010 FORMAT (//,10X,'Knight shift',//,10X,'K   tot    :',F15.6,' % '/,
     &        10X,'K   (-dia) :',F15.6,' % ')
99011 FORMAT (/,10X,A12,F15.6,5X,:,' 0:',F12.6,5X,'XC:',F12.6)
C
99012 FORMAT (21X,A,F15.6,2(8X,F12.6),/,(21X,A,F15.6,2(8X,F12.6)))
C
99013 FORMAT (18X,A,1X,A,F15.6,2(8X,F12.6),/,(21X,A,F15.6,2(8X,F12.6)))
99014 FORMAT (21X,A,43X,F12.6,/,(21X,A,43X,F12.6))
99015 FORMAT (21X,A,F15.6,2(8X,F12.6),/,(21X,A,F15.6,2(8X,F12.6)))
99016 FORMAT (/,10X,A12,F15.6,2X,:,A6,F12.6,3X,A5,F12.6)
99017 FORMAT (/,10X,A12,43X,F12.6)
99018 FORMAT ('# BUILDBOT: ',A,2X,A12,'  for IT =',I5,/,(1PE22.14))
      END
