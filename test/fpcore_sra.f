C*==fpcoresra.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPCORESRA(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,
     &                     NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,BCORS,
     &                     ECOR_LT,NLCORE,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SUBROUTINE TO CALCULATE THE RELATIVISTIC CORE WAVE             *
C   *   FUNCTIONS FOR A SPIN-DEPENDENT POTENTIAL                       *
C   *                                                                  *
C   *   FOR A GIVEN POTENTIAL THE NUMBER OF CORE AND VALENCE           *
C   *   ELECTRONS IS DETERMINED AND ALL CORE STATES THEN CALCULATED    *
C   *   > THE ROUTINE IS ORGANIZED AS DESCLAUXS ROUTINE <RESLD>        *
C   *     BUT FINDS THE CORRECTION TO THE E-EIGENVALUE AND THE         *
C   *     MATCHING PARAMETERS BY A NEWTON RAPHSON ALGORITHM            *
C   *     THIS IS IN VARIANCE TO THE METHOD SUGGESTED BY CORTONA       *
C   *   > SET THE SWITCH 'CHECK'  TO COPARE E-EIGENVALUES WITH         *
C   *     RESULTS OBTAINED WITH THE CONVENTIONAL E-CORRECTION          *
C   *     ALGORITHM, WHICH WORKS ONLY IF NO COUPLING IS PRESENT !      *
C   *   > THE FUNCTIONS  {GC,FC}(I,J) J=1,NSOL ARE THE LINEAR          *
C   *     INDEPENDENT SOLUTIONS TO THE DIFFERENTIAL EQUATIONS WITH     *
C   *     KAPPA-CHARACTER I=1,NSOL;   FOR OUTWARD AND INWARD           *
C   *     INTEGRATION THE SAME ARRAYS ARE USED !                       *
C   *   > THE PROPER SOLUTIONS SATISFYING THE BOUNDARY CONDITIONS      *
C   *     AT R=0 AND(!) R=INFINITY ARE STORED IN {GCK,FCK}(K,S)        *
C   *     S,K=1,NSOL   SOLUTION S=1 FOR  KAPPA = - L - 1               *
C   *                           S=2 FOR  KAPPA = + L (IF EXISTENT)     *
C   *   > THE SWITCH FINITE_NUCLEUS SELECTS WHETHER                    *
C   *     A FINITE NUCLEUS SHOULD BE USED                              *
C   *                                                                  *
C   *   HYPERFINE FIELD SPLITTING INTRODUCED IF ICORE=1 MB JUN. 1995   *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_CALCMODE,ONLY:ISMQHFI,NONMAG
      USE MOD_RMESH,ONLY:NMMAX,NRMAX,NPANMAX,FLMSF,ISFLM,NRSFTOT,NPAN,
     &    JRCUT,JRWS,DRDI,R2DRDI,R,NRCMAX,FINITE_NUCLEUS
      USE MOD_TYPES,ONLY:NT,NTMAX,LCXRAY,NCXRAY,IMT,Z,BT,VT,CTL,NCORT,
     &    ITBOT,ITTOP,RHOSPNC,RHOCHRC
      USE MOD_ANGMOM,ONLY:NMEHFMAX,NLMAX,NKMMAX,QMOFF,QOFF,QMDIA,QDIA,
     &    SMOFF,SOFF,SMDIA,SDIA
      IMPLICIT NONE
C*--FPCORESRA42
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 UNEND,TOLVAR,DVSTEP,TOL_RELCHANGE
      PARAMETER (UNEND=600.0D0,TOLVAR=1.0D-8,DVSTEP=0.01D0,
     &           TOL_RELCHANGE=0.1D0)
      INTEGER ITERMAX,NLSHELLMAX,NVAR
      PARAMETER (ITERMAX=400,NLSHELLMAX=15,NVAR=2)
C
C Dummy arguments
C
      INTEGER IPRINT,ITXRAY,NCSTMAX,NLCORE
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       ECOR_LT(NLCORE,NTMAX),FCOR(NRMAX,2,NCSTMAX),
     &       GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
C
C Local variables
C
      REAL*8 AUX,AUX1(:),BB(:),BEX,BHF(2,2),BHF1(2,2),BHF2(2,2),BSH,
     &       BSOL,CGD(2),CGMD(2),CGO,DEC,DEDV(4,4),DOVRC(:),DP(:,:,:),
     &       DQ(:,:,:),DRDIC(:),DROVRN(:),DV(4),DVDE(4,4),EC,ECC,
     &       ECORDUM(:),ECORTAB(:,:),ELIM,ERR(4),ERRNEW(4),FC(:,:,:),
     &       FCK(:,:,:),GC(:,:,:),GCK(:,:,:),MAX_RELCHANGE,MJ,NORM,
     &       PIW(2,2),POW(2,2),QIW(2,2),QOW(2,2),R2DRDIC(:),RAT,RC(:),
     &       RELCHANGE,RINT(:),RINTCHR(:),RINTSPN(:),RMAT4(4,4),RNUC,
     &       RSQ_F00,SCL,SHF(:,:,:),SIMP,SPLIT(:),SPLIT1(:),SPLIT2(:,:),
     &       SPLIT3(:,:),SZ,SZ1,SZCORDUM(:),SZD,VAL,VAR(4),VARNEW(4),
     &       VARTAB(4,20),VV(:),VZ,W,WP(:,:,:),WQ(:,:,:)
      LOGICAL CHECK
      INTEGER I,IC,IC1,IC2,ICST,ICSTDUM,IE,IFLAG,II,IL,ILC,ILSHELL,IM,
     &        IMIN,IN,IPIV4(4),IRMTIN,ISF1,ISH,ISTART,IT,ITER,IV,J,JLIM,
     &        JTOP,JV,K,KAP(2),KAP1,KAP2,KC,L,LLL,LOOP,LQNTAB(15),
     &        MUEM05,N,NLSHELL,NMATCH,NN,NODE,NQN,NQNTAB(15),NRC,NSH,
     &        NSOL,NZERO,NZEROTRIED,S,T
      INTEGER IKAPMUE
      LOGICAL RNON0
      REAL*8 RNUCTAB
      CHARACTER*10 TXTB(1:5)
      CHARACTER*3 TXTK(4)
      CHARACTER*1 TXT_L(0:3)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BB,FC,GC,DP,DQ,RC,FCK,GCK,AUX1,DRDIC,DOVRC,DROVRN
      ALLOCATABLE R2DRDIC,ECORTAB
      ALLOCATABLE SHF,SPLIT,SPLIT1,SPLIT2,SPLIT3,VV,WP,WQ
      ALLOCATABLE RINT,RINTCHR,RINTSPN
      ALLOCATABLE SZCORDUM,ECORDUM
      ALLOCATE (BB(NRMAX*3),FC(2,2,NRMAX*2),GC(2,2,NRMAX*2))
      ALLOCATE (DP(2,2,NRMAX*2),DQ(2,2,NRMAX*2),RC(NRMAX*3))
      ALLOCATE (FCK(2,2,NRMAX*2),GCK(2,2,NRMAX*2),AUX1(2*NRMAX))
      ALLOCATE (DRDIC(NRMAX*3),DOVRC(NRMAX*3),DROVRN(2*NRMAX))
      ALLOCATE (R2DRDIC(NRMAX*3),ECORTAB(120,NTMAX))
      ALLOCATE (SHF(2,2,NMEHFMAX),SPLIT(NMEHFMAX),SPLIT1(NMEHFMAX),
     &          SPLIT2(NMEHFMAX,NTMAX),SPLIT3(NMEHFMAX,NTMAX),
     &          VV(NRMAX*3),WP(2,2,NRMAX*2),WQ(2,2,NRMAX*2))
      ALLOCATE (SZCORDUM(200),ECORDUM(200))
C
C
C ----------------------------------------------------------------------
C  the core states are assumed to have only  s,p,d,f
C  but no higher l-components, so  NLCORE=4
C
C
      DATA TXTB/'B_nses','B_nseo','B_ssc ','B_sum ','B_tot '/
      DATA TXT_L/'s','p','d','f'/
      DATA TXTK/'1/2','3/2','5/2','7/2'/
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
      DATA CHECK/.FALSE./
C
      ALLOCATE (RINT(2*NRMAX),RINTCHR(NRMAX),RINTSPN(NRMAX))
C
      CALL RINIT(120*NTMAX,ECORTAB)
      CALL RINIT(4*20,VARTAB)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
C------------------------------------- cycle loop in case of empty spere
         IF ( Z(IT).EQ.0 ) CYCLE
C
         SCL = 1D0
C
         IM = IMT(IT)
         JTOP = JRWS(IM)
C
C------------------------------ use expanded radial mesh for core states
C
         NRC = 2*JTOP
         IF ( NRC.GT.NRCMAX ) STOP 'in <CORE>: NRC > NRCMAX'
C
C------------------------------------------------------------------- OLD
         RAT = R(JTOP,IM)/R(JTOP-1,IM)
C
         DO N = 1,JTOP
            RC(N) = R(N,IM)
            DRDIC(N) = DRDI(N,IM)
            R2DRDIC(N) = R2DRDI(N,IM)
            DOVRC(N) = DRDIC(N)/RC(N)
         END DO
C
         DO N = (JTOP+1),NRC
            RC(N) = RAT*RC(N-1)
            DRDIC(N) = (RAT-1.0D0)*RC(N-1)
            R2DRDIC(N) = RC(N)*RC(N)*DRDIC(N)
            DOVRC(N) = DRDIC(N)/RC(N)
         END DO
C------------------------------------------------------------------- OLD
C
C     integration boundary for hyperfine fields for finite nucleus
C     2 mesh points more for executing appropriate interpolation
C     to real nuclear radius RNUC
C
         IF ( FINITE_NUCLEUS ) THEN
            RNUC = RNUCTAB(Z(IT))
            IN = 1
            DO WHILE ( RC(IN).LE.RNUC )
               IN = IN + 1
            END DO
            JLIM = IN + 2
            IF ( MOD(JLIM,2).EQ.0 ) JLIM = JLIM - 1
C
            DO I = 1,NRC
               DROVRN(I) = (RC(I)/RNUC)**3*DRDIC(I)
            END DO
         END IF
C
C=======================================================================
C
         LOOP = 1
         BCOR(IT) = 0.0D0
         BCORS(IT) = 0.0D0
         DO I = 1,NMEHFMAX
            SPLIT2(I,IT) = 0.0D0
            SPLIT3(I,IT) = 0.0D0
         END DO
C
         DO N = 1,NRMAX
            RHOCHRC(N,IT) = 0.0D0
            RHOSPNC(N,IT) = 0.0D0
         END DO
C
         DO IL = 1,NLCORE
            ECOR_LT(IL,IT) = 0.0D0
         END DO
C
         DO N = 1,JRWS(IM)
            VV(N) = VT(N,IT)
            BB(N) = BT(N,IT)
         END DO
C
         DO N = (JRWS(IM)+1),NRC
            VV(N) = 0.0D0
            BB(N) = 0.0D0
         END DO
C
         NLSHELL = 0
         IF ( Z(IT).GT.2 ) NLSHELL = 1
         IF ( Z(IT).GT.10 ) NLSHELL = 3
         IF ( Z(IT).GT.18 ) NLSHELL = 5
         IF ( Z(IT).GT.30 ) NLSHELL = 6
         IF ( Z(IT).GT.36 ) NLSHELL = 8
         IF ( Z(IT).GT.48 ) NLSHELL = 9
         IF ( Z(IT).GT.54 ) NLSHELL = 11
         IF ( Z(IT).GT.70 ) NLSHELL = 12
         IF ( Z(IT).GT.80 ) NLSHELL = 13
         IF ( Z(IT).GT.86 ) NLSHELL = 15
C
         IF ( NCORT(IT).NE.0 ) THEN
            NLSHELL = 0
            N = 0
            DO ILSHELL = 1,NLSHELLMAX
               L = LQNTAB(ILSHELL)
               N = N + 2*(2*L+1)
               IF ( N.EQ.NCORT(IT) ) NLSHELL = ILSHELL
            END DO
            IF ( NLSHELL.EQ.0 ) THEN
               WRITE (6,*) 'NLSHELL not found for IT=',IT,' NCORT=',
     &                     NCORT(IT)
               STOP ' in <CORE>'
            END IF
         END IF
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) WRITE (6,99001) IT,Z(IT)
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C                   ---------------------------------------
C                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
C                   ---------------------------------------
         ICSTDUM = 0
C
         IC = 0
         DO ILSHELL = 1,NLSHELL
            NQN = NQNTAB(ILSHELL)
            L = LQNTAB(ILSHELL)
            IL = L + 1
            ILC = MIN(NLMAX,IL)
            NSH = 2*(2*L+1)
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C     SKIP SHELL IF NOT NEEDED IN A  XRAY - CALCULATION
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            IF ( ITXRAY.NE.0 ) THEN
               IF ( IT.NE.ITXRAY ) CYCLE
               IF ( (NQN.NE.NCXRAY(IT)) .OR. (L.NE.LCXRAY(IT)) ) CYCLE
               DO ICST = 1,NCSTMAX
                  DO KC = 1,2
                     DO N = 1,NRMAX
                        GCOR(N,KC,ICST) = 0.0D0
                        FCOR(N,KC,ICST) = 0.0D0
                     END DO
                  END DO
               END DO
            END IF
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
            ISH = 0
            BSH = 0.0D0
            DO I = 1,NMEHFMAX
               SPLIT1(I) = 0.0D0
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MUEM05 = -L - 1, + L
               MJ = MUEM05 + 0.5D0
C
C
               KAP1 = -L - 1
               KAP2 = L
               KAP(1) = KAP1
               KAP(2) = KAP2
C
               LLL = L*(L+1)
               IF ( ABS(MJ).GT.L ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
               DO S = 1,NSOL
                  NZEROTRIED = 0
                  IC = IC + 1
                  ISH = ISH + 1
C
                  T = 3 - S
C                  ----------------------------------------
C                   USE EC OF PREVIOUS RUNS AS START-VALUE
C                   TAKE SPIN-ORBIT SPLITTING INTO ACCOUNT
C                  ----------------------------------------
                  IF ( ISH.GT.1 ) THEN
                     EC = ECORTAB(IC-1,IT)
                     IF ( S.EQ.2 ) EC = ECORTAB(IC-1,IT)*1.1D0
                     IF ( ISH.GE.4 ) EC = ECORTAB(IC-2,IT)
                     GOTO 15
                  END IF
C
C
C                                      --------------------
C                                         FIND  E-LIMIT
C                                      --------------------
                  IF ( LLL.EQ.0 ) THEN
                     ELIM = -2*DBLE(Z(IT)**2)/(1.5D0*NQN*NQN)
                  ELSE
                     ELIM = VV(1) + LLL/RC(1)**2
                     DO N = 2,NRC
                        VAL = VV(N) + LLL/RC(N)**2
                        IF ( VAL.LE.ELIM ) ELIM = VAL
                     END DO
                  END IF
C
                  EC = -DBLE(Z(IT)**2)/(2.0D0*NQN*NQN)
C
                  ISTART = 1
 5                CONTINUE
                  IF ( EC.LE.ELIM ) EC = ELIM*0.7D0
C
C                                      --------------------
C                                         FIND    NZERO
C                                      --------------------
                  DO N = 1,(NRC-1)
                     IF ( (VV(N)-EC)*RC(N)**2.GT.UNEND ) THEN
                        IF ( MOD(N,2).EQ.0 ) THEN
                           NZERO = N + 1
                        ELSE
                           NZERO = N
                        END IF
                        GOTO 10
                     END IF
                  END DO
                  NZERO = NRC - 1
                  NZEROTRIED = NZEROTRIED + 1
                  IF ( NZEROTRIED.GT.20 ) THEN
                     WRITE (6,99002) IT,NQN,L,EC,(NRC-1)
                     STOP 'in <FPCORESRA>    NZERO not found'
                  END IF
C                                      --------------------
C                                         FIND    NMATCH
C                                      --------------------
 10               CONTINUE
                  N = NZERO + 1
                  DO NN = 1,NZERO
                     N = N - 1
                     IF ( (VV(N)+LLL/RC(N)**2-EC).LT.0D0 ) THEN
                        NMATCH = N
                        GOTO 15
                     END IF
                  END DO
                  WRITE (6,99003) IT,NQN,L,EC
C
 15               CONTINUE
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),EC,L,MJ,'OUT',VV,BB,
     &                              RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,
     &                              DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,Z(IT),IM,JRCUT,NPAN,NRC,
     &                              NMMAX,NPANMAX)
C
                  NODE = 0
                  DO N = 2,NMATCH
                     IF ( GC(S,S,N)*GC(S,S,N-1).LT.0D0 ) NODE = NODE + 1
                  END DO
C
C
                  IF ( IPRINT.GE.2 ) WRITE (6,99013) IT,NQN,L,KAP(S),
     &                 (2*MUEM05+1),IC,ISH,0,EC,NMATCH,RC(NMATCH),NZERO,
     &                 RC(NZERO),NODE,(GC(S,S,NMATCH)/GC(S,S,NMATCH-1))
C
                  IF ( NODE.NE.(NQN-L-1) ) THEN
                     IF ( NODE.GT.(NQN-L-1) ) THEN
                        EC = 1.2D0*EC
                     ELSE
                        EC = 0.8D0*EC
                     END IF
                     GOTO 5
                  ELSE IF ( (GC(S,S,NMATCH)/GC(S,S,NMATCH-1).LE.0D0)
     &                      .OR. 
     &                      (GC(S,S,NMATCH)/GC(S,S,NMATCH-1).GE.1D0) )
     &                      THEN
                     EC = 0.9D0*EC
                     GOTO 5
                  END IF
C
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),EC,L,MJ,'INW',VV,BB,
     &                              RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,
     &                              DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,Z(IT),IM,JRCUT,NPAN,NRC,
     &                              NMMAX,NPANMAX)
C
C                                      --------------------
C                                       START VALUES FOR
C                                           PARAMETERS
C                                      --------------------
C
                  VAR(1) = EC
                  VAR(2) = POW(S,S)/PIW(S,S)
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  DV(:) = VAR(:)
C
                  ITER = 0
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 20               CONTINUE
                  ITER = ITER + 1
C
C                         ----------------------------------
C                         CHECK WHETHER NUMBER OF NODES O.K.
C                         ----------------------------------
                  IF ( ITER.GT.1 ) THEN
                     NODE = 0
                     DO N = 2,(NMATCH-1)
                        IF ( GC(S,S,N)*GC(S,S,N-1).LT.0D0 )
     &                       NODE = NODE + 1
                     END DO
                     IF ( IPRINT.GE.3 ) WRITE (6,99013) IT,NQN,L,KAP(S),
     &                    (2*MUEM05+1),IC,ISH,ITER,EC,NMATCH,RC(NMATCH),
     &                    NZERO,RC(NZERO),NODE,
     &                    (GC(S,S,NMATCH)/GC(S,S,NMATCH-1))
C
                     IF ( NODE.NE.(NQN-L-1) ) THEN
                        IF ( NODE.GT.(NQN-L-1) ) THEN
                           EC = 1.2D0*EC
                        ELSE
                           EC = 0.8D0*EC
                        END IF
                        ISTART = ISTART + 1
                        IF ( ISTART.LT.20 ) GOTO 5
                     END IF
                  END IF
C
                  DO IV = 2,NVAR
                     DO JV = 1,NVAR
                        VARNEW(JV) = VAR(JV)
                     END DO
                     VARNEW(IV) = VAR(IV) + DV(IV)*DVSTEP
C
                     IF ( RNON0(VAR(IV)) ) THEN
                        IF ( ABS(DV(IV)/VAR(IV)).LT.TOLVAR ) VARNEW(IV)
     &                       = VAR(IV)
     &                         *(1.0D0+DSIGN(DVSTEP*TOLVAR,DV(IV)))
                     ELSE IF ( .NOT.NONMAG ) THEN
                        WRITE (6,99004) ' VAR(',IV,
     &                                  ') = 0 for (T,N,L,K,M;S,NSOL) ',
     &                                  IT,NQN,L,KAP(S),(2*MUEM05+1),
     &                                  '/2  ',S,NSOL
                     END IF
C
                     CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                     DO IE = 1,NVAR
                        IF ( ABS(ERRNEW(IE)-ERR(IE)).LT.1D-16 ) THEN
CJMJMJMJM                           DEDV(IE,IV) = 0.0D0
                           DEDV(IE,IV) = 1.0D-18
                           IF ( (IE.EQ.IV) .AND. NONMAG ) DEDV(IE,IV)
     &                          = 1.0D0
                        ELSE
                           DEDV(IE,IV) = (ERRNEW(IE)-ERR(IE))
     &                        /(VARNEW(IV)-VAR(IV))
                        END IF
                     END DO
                  END DO
C
                  DO JV = 1,NVAR
                     VARNEW(JV) = VAR(JV)
                  END DO
                  VARNEW(1) = VAR(1) + DV(1)*DVSTEP
                  IF ( ABS(DV(1)/VAR(1)).LT.TOLVAR ) VARNEW(1) = VAR(1)
     &                 *(1.0D0+DSIGN(DVSTEP*TOLVAR,DV(1)))
C
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'OUT',
     &                              VV,BB,RC,DRDIC,DOVRC,NMATCH,NZERO,
     &                              GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,
     &                              CGD,CGMD,CGO,Z(IT),IM,JRCUT,NPAN,
     &                              NRC,NMMAX,NPANMAX)
C
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'INW',
     &                              VV,BB,RC,DRDIC,DOVRC,NMATCH,NZERO,
     &                              GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,
     &                              CGD,CGMD,CGO,Z(IT),IM,JRCUT,NPAN,
     &                              NRC,NMMAX,NPANMAX)
C
                  CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                  DO IE = 1,NVAR
                     DEDV(IE,1) = (ERRNEW(IE)-ERR(IE))
     &                            /(VARNEW(1)-VAR(1))
                  END DO
C
                  IF ( NVAR.NE.4 ) THEN
C
C                    CALL RINVGJ(DVDE,DEDV,4,NVAR)
                     CALL RMATINV3(NVAR,4,IPIV4,DEDV,RMAT4,DVDE)
C
                  ELSE IF ( ABS(DEDV(3,3)).LT.1D-15 .OR. ABS(DEDV(4,4))
     &                      .LT.1D-15 ) THEN
C
C                    CALL RINVGJ(DVDE,DEDV,4,2)
                     CALL RMATINV3(2,4,IPIV4,DEDV,RMAT4,DVDE)
C
                     DVDE(3,3) = 0D0
                     DVDE(4,4) = 0D0
C
                  ELSE
C
C                    CALL RINVGJ(DVDE,DEDV,4,NVAR)
                     CALL RMATINV3(NVAR,4,IPIV4,DEDV,RMAT4,DVDE)
C
                  END IF
C
C--------- damp correction if relative change greater than TOL_RELCHANGE
C
                  IFLAG = 0
                  MAX_RELCHANGE = 0D0
                  DO IV = 1,NVAR
                     DV(IV) = 0.0D0
                     DO IE = 1,NVAR
                        DV(IV) = DV(IV) + DVDE(IV,IE)*ERR(IE)
                     END DO
                     IF ( IV.EQ.1 ) THEN
                        RELCHANGE = ABS(DV(IV)/VAR(IV))
                        IF ( RELCHANGE.GT.TOL_RELCHANGE ) THEN
                           IFLAG = 1
                           MAX_RELCHANGE = MAX(MAX_RELCHANGE,RELCHANGE)
                        END IF
                     END IF
                  END DO
C
                  IF ( IFLAG.EQ.1 ) DV(1:NVAR) = DV(1:NVAR)
     &                 *TOL_RELCHANGE/MAX_RELCHANGE
C
                  DO IV = 1,NVAR
                     VAR(IV) = VAR(IV) - DV(IV)
                  END DO
C
                  IF ( VAR(1).GT.0.0D0 ) THEN
                     WRITE (6,*) ' WARNING FROM <CORE> E=',VAR(1)
                     VAR(1) = -0.2D0
                  END IF
C
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),VAR(1),L,MJ,'OUT',VV,
     &                              BB,RC,DRDIC,DOVRC,NMATCH,NZERO,GC,
     &                              FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,Z(IT),IM,JRCUT,NPAN,NRC,
     &                              NMMAX,NPANMAX)
C
                  CALL FPCORESRADIR(IT,CTL(IT,ILC),VAR(1),L,MJ,'INW',VV,
     &                              BB,RC,DRDIC,DOVRC,NMATCH,NZERO,GC,
     &                              FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,Z(IT),IM,JRCUT,NPAN,NRC,
     &                              NMMAX,NPANMAX)
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  EC = VAR(1)
C
                  IF ( IPRINT.GE.2 ) WRITE (6,99005) LOOP,SCL,VAR(1),
     &                 (VAR(IV),IV=1,4),(DV(IV),IV=1,4),(ERR(IE),IE=1,4)
C
C----------------------------------  CHECK RELATIVE CHANGE IN PARAMETERS
C ----------------------- PARAMETERS 3 AND 4 = 0 FOR PARAMAGNETIC CASE !
                  IF ( ITER.LT.ITERMAX ) THEN
                     DO IV = 1,NVAR
                        VARTAB(IV,ISH) = VAR(IV)
                        IF ( (ABS(VAR(IV))+ABS(VAR(IV))).LT.1.0D-12 ) ! changed d-20 to D-12
     &                       THEN
                        ELSE IF ( ABS(DV(IV)/VAR(IV)).GT.TOLVAR ) THEN
                           GOTO 20
                        END IF
                     END DO
                  ELSE
                     WRITE (6,99006) ITERMAX,(VAR(IV),IV=1,4),
     &                               (DV(IV),IV=1,4),(ERR(IE),IE=1,4)
                  END IF
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C                         ---------------------------------
C                         NORMALIZE WAVEFUNCTIONS ACCORDING
C                               TO MATCHING CONDITIONS
C                         ---------------------------------
C
C                                    INWARD - SOLUTION
                  DO N = NMATCH,NZERO
                     DO J = 1,NSOL
                        DO I = 1,NSOL
                           GC(I,J,N) = GC(I,J,N)*VAR(2)
                           FC(I,J,N) = FC(I,J,N)*VAR(2)
                        END DO
                     END DO
                  END DO
C
                  IF ( NSOL.EQ.2 ) THEN
C                                   OUTWARD - SOLUTION
                     DO N = 1,(NMATCH-1)
                        DO I = 1,NSOL
                           GC(I,T,N) = GC(I,T,N)*VAR(3)
                           FC(I,T,N) = FC(I,T,N)*VAR(3)
                        END DO
                     END DO
C                                    INWARD - SOLUTION
                     DO N = NMATCH,NZERO
                        DO I = 1,NSOL
                           GC(I,T,N) = GC(I,T,N)*VAR(4)
                           FC(I,T,N) = FC(I,T,N)*VAR(4)
                        END DO
                     END DO
                  END IF
C
C                                    SUM FOR EACH KAPPA
                  DO N = 1,NZERO
                     DO K = 1,NSOL
                        GCK(K,S,N) = 0.0D0
                        FCK(K,S,N) = 0.0D0
                        DO J = 1,NSOL
                           GCK(K,S,N) = GCK(K,S,N) + GC(K,J,N)
                           FCK(K,S,N) = FCK(K,S,N) + FC(K,J,N)
                        END DO
                     END DO
                  END DO
C
C                       -----------------------------------
C                       CALCULATE  NORM  AND NORMALIZE TO 1
C                       -----------------------------------
C
                  DO N = (NZERO+1),JTOP
                     DO K = 1,NSOL
                        GCK(K,S,N) = 0.0D0
                        FCK(K,S,N) = 0.0D0
                     END DO
                  END DO
C
                  CALL RINIT(2*NRMAX,RINT)
C
                  IRMTIN = JRCUT(1,IM)
                  ISF1 = ISFLM(1,IM)
C
                  DO N = 1,IRMTIN
                     DO K = 1,NSOL
                        RINT(N) = RINT(N) + R2DRDI(N,IM)
     &                            *(GCK(K,S,N)*GCK(K,S,N)+FCK(K,S,N)
     &                            *FCK(K,S,N))
                     END DO
                  END DO
C
                  DO J = 1,NRSFTOT(IM)
                     N = IRMTIN + J
                     DO K = 1,NSOL
                        RINT(N) = RINT(N) + R2DRDI(N,IM)
     &                            *(GCK(K,S,N)*GCK(K,S,N)+FCK(K,S,N)
     &                            *FCK(K,S,N))*FLMSF(J,ISF1,IM)/SQRT_4PI
                     END DO
                  END DO
C
                  CALL RRADINT(IM,RINT,AUX)
C
C ------------------------------ OMIT NORMALIZATION FOR XRAY CALCULATION
C ------------------------------ TO RECOVER OLD (ERROUNOUS DATA) -------
                  IF ( ITXRAY.GT.0 ) THEN
                     NORM = 1.0D0
                  ELSE
                     NORM = 1.0D0/SQRT(AUX)
                  END IF
C
                  DO N = 1,MAX(NZERO,JTOP)
                     DO K = 1,NSOL
                        GCK(K,S,N) = GCK(K,S,N)*NORM
                        FCK(K,S,N) = FCK(K,S,N)*NORM
                     END DO
                  END DO
C
C                       -----------------------------------
C                            CALCULATE  SPIN CHARACTER
C                       -----------------------------------
C
                  W = R2DRDIC(1)
                  SZ = 0.0D0
                  DO K = 1,NSOL
                     SZ = SZ + W*(GCK(K,S,1)**2*CGD(K)+FCK(K,S,1)
     &                    **2*CGMD(K))
                  END DO
C
                  SIMP = -1.0D0
                  DO N = 2,NZERO
                     SIMP = -SIMP
                     W = (3.0D0+SIMP)*R2DRDIC(N)
                     DO K = 1,NSOL
                        SZ = SZ + W*(GCK(K,S,N)**2*CGD(K)+FCK(K,S,N)
     &                       **2*CGMD(K))
                     END DO
                  END DO
C
                  N = NZERO
                  W = R2DRDIC(N)
                  DO K = 1,NSOL
                     SZ = SZ - W*(GCK(K,S,N)**2*CGD(K)-FCK(K,S,N)
     &                    **2*CGMD(K))
                  END DO
C
C
                  IF ( NSOL.GT.1 ) THEN
C
                     W = R2DRDIC(1)
                     SZ = SZ + W*GCK(1,S,1)*GCK(2,S,1)*CGO*2
C
                     SIMP = -1.0D0
                     DO N = 2,NZERO
                        SIMP = -SIMP
                        W = (3.0D0+SIMP)*R2DRDIC(N)
                        SZ = SZ + W*GCK(1,S,N)*GCK(2,S,N)*CGO*2
                     END DO
C
                     N = NZERO
                     W = R2DRDIC(N)
                     SZ = SZ - W*GCK(1,S,N)*GCK(2,S,N)*CGO*2
C
                  END IF
C
                  SZ = SZ/3.0D0
C
C--- BEGIN --------- Consistency check on 2nd solution -----------------
C
C  +  in the large B regime the solver does not find the second solution
C     belonging to same mu and l and falls back to first sol
C  +  here we test whether sol2 = sol1, if so, we try to drive the
C     solver to the second solution, by setting appropriate initial
C     values for starting the search algorithm anew
C-----------------------------------------------------------------------
C               --------------------------------------------------
C               first: + estimate x-splitting parameter BEX
C                      + use state j = l + 1/2, mue = - l - 1
C               --------------------------------------------------
                  IF ( MUEM05.EQ.-L-1 .AND. L.GT.0 ) THEN
                     SZD = SZ
C
                     BEX = 0D0
                     BEX = R2DRDIC(1)
     &                     *(CGD(1)*GCK(1,S,1)**2+CGMD(1)*FCK(1,S,1)**2)
     &                     *BB(1)
C
                     SIMP = -1.0D0
                     DO N = 2,NZERO
                        SIMP = -SIMP
                        W = (3.0D0+SIMP)*R2DRDIC(N)
                        BEX = BEX + W*(CGD(1)*GCK(1,S,N)**2+CGMD(1)
     &                        *FCK(1,S,N)**2)*BB(N)
                     END DO
C
                     N = NZERO
                     BEX = BEX - R2DRDIC(N)
     &                     *(CGD(1)*GCK(1,S,N)**2+CGMD(1)*FCK(1,S,N)**2)
     &                     *BB(N)
                     BEX = BEX/3.0D0
C
                     IF ( IPRINT.GT.0 ) WRITE (6,'("BEX = ",e16.7)') BEX
                  END IF
C
                  IF ( NSOL.EQ.2 ) THEN
C
                     IF ( S.EQ.1 ) SZ1 = SZ
C
C
                  END IF
C
C              --------------------------------------------------------
C              check consistency of sol2 and reinitialise if necessary
C              --------------------------------------------------------
C              NB: ICST here still contains index of 1st sol
                  IF ( S.EQ.2 ) THEN
                     IF ( (SZ*SZCORDUM(ICSTDUM)).GT.0D0 ) THEN
                        VAR(3) = -VAR(3)
                        VAR(4) = -VAR(4)
                        IF ( IPRINT.GT.0 ) WRITE (6,'(9f15.7)') EC,
     &                       ECORDUM(ICSTDUM),SZ,SZCORDUM(ICSTDUM)
                        IF ( IPRINT.GT.0 ) WRITE (6,*) ITER
                        IF ( BEX.GT.0D0 ) THEN
C                  state mue=-l-1/2 is the higher of the j=l+1/2 states
                           IF ( SZD*SZ.GT.0D0 ) THEN
                              VAR(1) = ECORDUM(ICSTDUM) - 2D0*BEX
                           ELSE
                              VAR(1) = ECORDUM(ICSTDUM) + 2D0*BEX
                           END IF
C                  state mue=-l-1/2 is the lower of the j=l+1/2 states
                        ELSE IF ( SZD*SZ.GT.0D0 ) THEN
                           VAR(1) = ECORDUM(ICSTDUM) - 2D0*BEX
                        ELSE
                           VAR(1) = ECORDUM(ICSTDUM) + 2D0*BEX
                        END IF
                        IF ( ITER.LT.ITERMAX ) GOTO 20
                     END IF
C
C       check 3 - check whether the two states with same mu
C                 have opposite sz
                     IF ( ABS(SZ+SZ1).GT.0.1D0 .AND. ABS(BEX)
     &                    .GT.1D-6 .AND. ITER.LT.ITERMAX ) GOTO 20
                  END IF
C-----------------------------------------------------------------------
C--- END ----------- Consistency check on 2nd solution -----------------
C-----------------------------------------------------------------------
Cdiemo end
C
C
C
C                       -----------------------------------
C                       CALCULATE  CHARGE AND SPIN DENSITY
C                       -----------------------------------
C
                  DO N = 1,JRWS(IM)
                     DO K = 1,NSOL
                        RHOCHRC(N,IT) = RHOCHRC(N,IT)
     &                                  + (GCK(K,S,N)**2+FCK(K,S,N)**2)
                        RHOSPNC(N,IT) = RHOSPNC(N,IT)
     &                                  + (GCK(K,S,N)**2*CGD(K)
     &                                  -FCK(K,S,N)**2*CGMD(K))
                     END DO
                  END DO
C
                  IF ( NSOL.GT.1 ) THEN
                     DO N = 1,JRWS(IM)
                        RHOSPNC(N,IT) = RHOSPNC(N,IT) + GCK(1,S,N)
     &                                  *GCK(2,S,N)*CGO*2
                     END DO
                  END IF
C
C                         ------------------------------
C                         CALCULATE   HYPERFINE - FIELD
C                         ------------------------------
                  CALL FPCORESRAHFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,
     &                              RC,DRDIC,0.0D0,NZERO,NRC,JRCUT(0,IM)
     &                              ,NPAN(IM),NPANMAX)
                  IF ( FINITE_NUCLEUS ) THEN
                     CALL FPCORESRAHFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF1,GCK,
     &                                 FCK,RC,DRDIC,RNUC,JLIM,NRC,
     &                                 JRCUT(0,IM),NPAN(IM),NPANMAX)
                     CALL FPCORESRAHFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF2,GCK,
     &                                 FCK,RC,DROVRN,RNUC,JLIM,NRC,
     &                                 JRCUT(0,IM),NPAN(IM),NPANMAX)
C
                     DO I = 1,NSOL
                        DO J = 1,NSOL
                           BHF(I,J) = BHF(I,J) - BHF1(I,J) + BHF2(I,J)
                        END DO
                     END DO
                  END IF
C
                  BSOL = 0.0D0
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        BSOL = BSOL + BHF(I,J)
                        BSH = BSH + BHF(I,J)
                        BCOR(IT) = BCOR(IT) + BHF(I,J)
                     END DO
                  END DO
                  IF ( KAP1.EQ.-1 ) BCORS(IT) = BCORS(IT) + BHF(1,1)
                  ECORTAB(IC,IT) = EC
                  ECOR_LT(IL,IT) = ECOR_LT(IL,IT) + EC
C
C     ------------------
C     SPLIT HFF-FIELD
C     ------------------
                  IF ( ISMQHFI.EQ.1 ) THEN
                     CALL CORE_HFF_SPLIT(IM,RNUC,NZERO,KAP1,KAP2,NSOL,
     &                  MJ,GCK,FCK,NRC,SHF,S,NMEHFMAX,NKMMAX,RC,DRDIC,
     &                  SDIA,SMDIA,SOFF,SMOFF,QDIA,QMDIA,QOFF,QMOFF,
     &                  JLIM)
C
                     DO K = 1,NMEHFMAX
                        SPLIT(K) = 0.0D0
                        DO J = 1,NSOL
                           DO I = 1,NSOL
                              SPLIT(K) = SPLIT(K) + SHF(I,J,K)
                              SPLIT1(K) = SPLIT1(K) + SHF(I,J,K)
                              SPLIT2(K,IT) = SPLIT2(K,IT) + SHF(I,J,K)
                           END DO
                        END DO
                     END DO
                     DO K = 1,NMEHFMAX
                        IF ( KAP1.EQ.-1 ) SPLIT3(K,IT) = SPLIT3(K,IT)
     &                       + SHF(1,1,K)
                     END DO
                  END IF
C
Ccc         write(*,*) 'ECORE ',IC,EC
C
                  ICSTDUM = ICSTDUM + 1
                  SZCORDUM(ICSTDUM) = SZ
                  ECORDUM(ICSTDUM) = EC
C
C
C---------------------------------------------------- L-SHELL UNCOMPLETE
                  IF ( ISH.GE.NSH ) THEN
C----------------------------------------------------- L-SHELL COMPLETED
                     IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) THEN
                        WRITE (6,99009) IT,NQN,TXT_L(L),
     &                                  TXTK(IABS(KAP(S))),(2*MUEM05+1),
     &                                  KAP(S),ITER,EC,BSOL*.001D0,
     &                                  BSH*.001D0
                        IF ( ISMQHFI.EQ.1 ) THEN
                           DO K = 1,NMEHFMAX
                              WRITE (6,99010) TXTB(K),SPLIT(K)*.001D0,
     &                               SPLIT1(K)*.001D0
                           END DO
                           WRITE (6,99011) 'total error in %',
     &                            100.0D0*(1.0D0-SPLIT(4)/SPLIT(5))
                        END IF
                     END IF
C                              ----------------------------
C                              CHECK CONSISTENCY OF RESULTS
C                              ----------------------------
                     IF ( L.NE.0 ) THEN
                        IC1 = IC - NSH + 1
                        IC2 = IC
                        IF ( ECORTAB(IC2,IT).GE.ECORTAB(IC1,IT) ) THEN
                           IMIN = IC1
                           VZ = +1.0D0
                        ELSE
                           IMIN = IC2
                           VZ = -1.0D0
                        END IF
                        IFLAG = 0
                        II = 1
                        DO I = IC1 + 1,IC2,2
                           IF ( VZ*(ECORTAB(I,IT)-ECORTAB(I-II,IT))
     &                          .LT.0D0 ) IFLAG = 1
                           II = 2
                        END DO
                        IF ( ECORTAB(IC1+2,IT).GT.ECORTAB(IMIN,IT) )
     &                       IFLAG = 1
                        DO I = IC1 + 4,IC2 - 1,2
                           IF ( ECORTAB(I,IT).GT.ECORTAB(IMIN,IT) )
     &                          IFLAG = 1
                           IF ( VZ*(ECORTAB(I,IT)-ECORTAB(I-II,IT))
     &                          .GT.0D0 ) IFLAG = 1
                        END DO
C
                     END IF
                  ELSE IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) THEN
                     WRITE (6,99009) IT,NQN,TXT_L(L),TXTK(IABS(KAP(S))),
     &                               (2*MUEM05+1),KAP(S),ITER,EC,
     &                               BSOL*.001D0
                     IF ( ISMQHFI.EQ.1 ) THEN
                        DO K = 1,NMEHFMAX
                           WRITE (6,99010) TXTB(K),SPLIT(K)*.001D0
                        END DO
                        WRITE (6,99011) 'total error in %',
     &                                  100.0D0*(1.0D0-SPLIT(4)/SPLIT(5)
     &                                  )
                     END IF
                  END IF
C-----------------------------------------------------------------------
C
                  IF ( IPRINT.GE.1 ) WRITE (6,99012)
     &                 ((BHF(I,J)*.001D0,I=1,NSOL),J=1,NSOL)
C
C                          --------------------------------
C                            IF THE SWITCH CHECK IS SET:
C                            RECALCULATE THE EIGENVALUE
C                          USING THE CONVENTIONAL ALGORITHM
C                          --------------------------------
                  IF ( CHECK ) THEN
                     ECC = 0.95D0*EC
C
 22                  CONTINUE
                     CALL FPCORESRADIR(IT,CTL(IT,ILC),ECC,L,MJ,'OUT',VV,
     &                                 BB,RC,DRDIC,DOVRC,NMATCH,NZERO,
     &                                 GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,
     &                                 QIW,CGD,CGMD,CGO,Z(IT),IM,JRCUT,
     &                                 NPAN,NRC,NMMAX,NPANMAX)
C
                     CALL FPCORESRADIR(IT,CTL(IT,ILC),ECC,L,MJ,'INW',VV,
     &                                 BB,RC,DRDIC,DOVRC,NMATCH,NZERO,
     &                                 GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,
     &                                 QIW,CGD,CGMD,CGO,Z(IT),IM,JRCUT,
     &                                 NPAN,NRC,NMMAX,NPANMAX)
C
                     NORM = POW(S,S)/PIW(S,S)
                     DO N = NMATCH,NZERO
                        GC(S,S,N) = GC(S,S,N)*NORM
                        FC(S,S,N) = FC(S,S,N)*NORM
                     END DO
C
                     NORM = 0.0D0
C
                     CALL RINIT(2*NRMAX,RINT)
C
                     DO N = 1,NZERO
                        DO K = 1,NSOL
                           RINT(N) = RINT(N) + R2DRDIC(N)
     &                               *(GCK(K,S,N)*GCK(K,S,N)+FCK(K,S,N)
     &                               *FCK(K,S,N))
                        END DO
                     END DO
C
                     CALL RRADINT_R(IM,RINT,AUX1)
C
                     NORM = AUX1(MIN(NZERO,JRCUT(NPAN(IM),IM)))
C
                     DEC = POW(S,S)
     &                     *(QOW(S,S)-RC(NMATCH)*CTL(IT,ILC)*FC(S,S,
     &                     NMATCH))/NORM
                     ECC = ECC + DEC
                     IF ( ABS(DEC/ECC).GT.TOLVAR ) GOTO 22
                     WRITE (6,'(7X,''CHECK-E:'',10X,F12.5,/)') ECC
                  END IF
C
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C           STORE CORE WAVE FUNCTIONS IF REQUIRED
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                  IF ( ITXRAY.NE.0 ) THEN
                     IF ( NSOL.EQ.2 ) THEN
                        IF ( S.EQ.2 ) THEN
                           ICST = L + 1 + MUEM05
                        ELSE
                           ICST = L + 1 + MUEM05 + 2*L + 1
                        END IF
                     ELSE IF ( MJ.LT.0 ) THEN
                        ICST = 2*L + 1
                     ELSE
                        ICST = 4*L + 2
                     END IF
C
                     MM05COR(ICST) = MUEM05
                     NKPCOR(ICST) = NSOL
                     KAPCOR(ICST) = KAP(S)
                     IKMCOR(ICST,1) = IKAPMUE(KAP(S),MUEM05)
                     IZERO(ICST) = MIN(NZERO,JRWS(IM))
                     SZCOR(ICST) = SZ
                     ECOR(ICST) = EC
C
                     DO N = 1,IZERO(ICST)
                        GCOR(N,1,ICST) = GCK(S,S,N)
                        FCOR(N,1,ICST) = FCK(S,S,N)
                     END DO
                     IF ( NSOL.EQ.2 ) THEN
                        DO N = 1,IZERO(ICST)
                           GCOR(N,2,ICST) = GCK(T,S,N)
                           FCOR(N,2,ICST) = FCK(T,S,N)
                        END DO
                        IKMCOR(ICST,2) = IKAPMUE(KAP(T),MUEM05)
                     END IF
                  END IF
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
               END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) WRITE (6,99007) BCOR(IT)
     &        *.001D0,BCORS(IT)*.001D0
         IF ( ISMQHFI.EQ.1 ) THEN
            DO N = 1,NMEHFMAX
               WRITE (6,99008) SPLIT2(N,IT)*.001D0,SPLIT3(N,IT)*.001D0
            END DO
            WRITE (6,99011) 'total error',
     &                      100.0D0*(1.0D0-SPLIT2(4,IT)/SPLIT2(5,IT))
            WRITE (6,99011) 'total error',
     &                      100.0D0*(1.0D0-SPLIT3(4,IT)/SPLIT3(5,IT))
         END IF
C
         IF ( (ITXRAY.EQ.0) .OR. (IT.EQ.ITXRAY) ) THEN
C
            IRMTIN = JRCUT(1,IM)
            ISF1 = ISFLM(1,IM)
C
            DO N = 1,IRMTIN
               RINTCHR(N) = RHOCHRC(N,IT)*R2DRDI(N,IM)
               RINTSPN(N) = RHOSPNC(N,IT)*R2DRDI(N,IM)
            END DO
C
            DO J = 1,NRSFTOT(IM)
               N = IRMTIN + J
               RSQ_F00 = R2DRDI(N,IM)*FLMSF(J,ISF1,IM)/SQRT_4PI
C
               RINTCHR(N) = RHOCHRC(N,IT)*RSQ_F00
               RINTSPN(N) = RHOSPNC(N,IT)*RSQ_F00
            END DO
C
            CALL RRADINT(IM,RINTCHR,AUX)
C
            WRITE (6,99014) 'charge',IT,AUX
C
            CALL RRADINT(IM,RINTSPN,AUX)
C
            WRITE (6,99014) ' spin ',IT,AUX
C
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( NT.EQ.0 ) THEN
         OPEN (92,FILE='RHOCOR.dat')
         IM = 1
         DO I = 1,JRWS(IM)
C
C------------ NOTE: RHOCHR (RHOSPN) and RHO2NS differ by factor SQRT_4PI
C
            W = R(I,IM)**2
C
            WRITE (92,'(12E14.5)') R(I,IM),
     &                             (RHOCHRC(N,IT)*W,RHOSPNC(N,IT)*W,
     &                             IT=ITBOT,ITTOP)
         END DO
         CLOSE (92)
         WRITE (6,*) ' core charge density written to   RHOCOR.dat '
      END IF
C
      DEALLOCATE (RINT,RINTCHR,RINTSPN)
C
99001 FORMAT (/,'  ATOM   IT      : ',I5,/,'  ATOMIC NUMBER  : ',I5,//,
     &        '  IT',10X,'MUE  KAP ITER    ENERGY       B  (k-Gauss)  ')
99002 FORMAT ('  IT=',I4,' NQN=',I2,' L=',I2,' EC=',F10.3,
     &        '  NZERO set to  (NRC-1) =',I4)
99003 FORMAT (//,'  STOP IN <<CORE>>',/,'  IT=',I4,' NQN=',I2,' L=',I2,
     &        /,'  no matching-radius found for  EC=',F10.3)
99004 FORMAT (A,I1,A,5I4,A,2I2,A)
99005 FORMAT (' LOOP    =  ',I3,' BSCL=',F10.5,/,' E=',F14.7,' VAR  ',
     &        4E11.4,/,17X,' CORR ',4E11.4,/,17X,' ERR  ',4E11.4)
99006 FORMAT (' iteration not converged after',I3,' steps !',/,
     &        ' parameters:',4E18.10,/,' last corr.:',4E18.10,/,
     &        ' last error:',4E18.10)
99007 FORMAT (2X,57('-'),/,42X,F17.3,/,38X,'(S) ',F17.3,/,2X,57('*'),//)
99008 FORMAT (2X,57('-'),/,37X,F17.3,/,33X,'(S) ',F17.3,/,2X,57('*'),//)
99009 FORMAT (2I4,A1,A3,I3,'/2',2I4,2X,F15.8,F17.3,:,F17.3,/)
99010 FORMAT (A,:,32X,F17.3,:,F17.3,/)
99011 FORMAT (A,:,37X,F6.3)
99012 FORMAT (37X,F17.3)
99013 FORMAT (/,' IT=',I4,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,
     &        '/2    IC=',I3,'  ISH=',I2,/,' E(',I2,')   =',F15.5,/,
     &        ' NMATCH  =',I5,'    R=',F10.5,/,' NZERO   =',I5,'    R=',
     &        F10.5,/,' NODES   =',I5,'  RAT=',E10.4)
99014 FORMAT (10X,'integrated core ',A,' density for atom type ',I4,':',
     &        F12.8)
C
      END
C*==fpcoresradir.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPCORESRADIR(IT,C,E,L,MJ,WAY,VV,BB,RC,DRDIC,DOVRC,
     &                        NMATCH,NZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                        PIW,QIW,CGD,CGMD,CGO,Z,IM,JRCUT,NPAN,NRC,
     &                        NMMAX,NPANMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE RADIAL SPIN-POLARISED DIRAC EQUATIONS         *
C   *   FOR THE CORE WAVE FUNCTIONS                                    *
C   *                                                                  *
C   *   SIMILAR TO LOUCKS' METHOD TO SOLVE THE COUPLED SET OF          *
C   *   DIFFERENTIAL EQUATIONS                                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      IMPLICIT NONE
C*--FPCORESRADIR1178
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,INVMAX
      PARAMETER (MPSMAX=20,NPEMAX=20,INVMAX=3)
      REAL*8 TOL
      PARAMETER (TOL=1.0D-9)
      INTEGER ITMAX,NABM
      PARAMETER (ITMAX=50,NABM=4)
C
C Dummy arguments
C
      REAL*8 C,CGO,E,MJ
      INTEGER IM,IT,L,NMATCH,NMMAX,NPANMAX,NRC,NZERO,Z
      CHARACTER*3 WAY
      REAL*8 BB(NRC),CGD(2),CGMD(2),DOVRC(NRC),DP(2,2,NRC),DQ(2,2,NRC),
     &       DRDIC(NRC),FC(2,2,NRC),GC(2,2,NRC),PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),RC(NRC),VV(NRC),WP(2,2,NRC),WQ(2,2,NRC)
      INTEGER JRCUT(0:NPANMAX,NMMAX),NPAN(NMMAX)
C
C Local variables
C
      REAL*8 A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,B14,BB1,BB2,
     &       BC(0:NPEMAX),BETA,BH,BHLP(NABM+4),BOVA,BPP,BQQ,CG1,CG2,CG4,
     &       CG5,CG8,CSQR,D14,DET,DH,DHLP(NABM+4),DIFFA,DIFFB,DMUE,DVC,
     &       EMVPP,EMVQQ,GAM(2),GPM,H24,HLP(NABM+4),HLP1,KAP(2),MP1(2,2)
     &       ,MP2(2,2),MP3(2,2),MP4(2,2),MQ1(2,2),MQ2(2,2),MQ3(2,2),
     &       MQ4(2,2),P1(2,2),P2(2,2),P3(2,2),P4(2,2),PC(2,2,0:MPSMAX),
     &       PNEW(2,2),POLD(2,2),Q1(2,2),Q2(2,2),Q3(2,2),Q4(2,2),
     &       QC(2,2,0:MPSMAX),QNEW(2,2),QOLD(2,2),R14,RH,RHLP(NABM+4),
     &       RPWGPM,RR,SO2,SO6,SRK,TZ,V14,VC(0:NPEMAX),VH,VHLP(NABM+4),
     &       W1,W2,W3,W4
      INTEGER I,IPAN,IPANMATCH,IPANNZERO,IRK,IV,J,JCORR,K,KAP1,KAP2,M,
     &        MPS,N,NDIV,NEND,NHLP,NMESHLC,NN,NSOL,NSTART
      LOGICAL RNON0
      REAL*8 W5,W6,W7,X14,XH
      REAL*8 YLAG
      SAVE A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,B14,BB1,BB2,BC,
     &     BETA,BH,BHLP,BOVA,BPP,BQQ,CG1,CG2,CG4,CG5,CG8,CSQR,D14,DET,
     &     DH,DHLP,DIFFA,DIFFB,DMUE,DVC,EMVPP,EMVQQ,GAM,GPM,H24,HLP,
     &     HLP1,I,IPAN,IPANMATCH,IPANNZERO,IRK,IV,J,JCORR,K,KAP,KAP1,
     &     KAP2,M,MP1,MP2,MP3,MP4,MPS,MQ1,MQ2,MQ3,MQ4,N,NDIV,NEND,NHLP,
     &     NMESHLC,NN,NSOL,NSTART,P1,P2,P3,P4,PC,PNEW,POLD,Q1,Q2,Q3,Q4,
     &     QC,QNEW,QOLD,R14,RH,RHLP,RPWGPM,RR,SO2,SO6,SRK,TZ,V14,VC,VH,
     &     VHLP,W1,W2,W3,W4,W5,W6,W7,X14,XH
C
C*** End of declarations rewritten by SPAG
C
      H24 = 1.0D0/24.0D0
      DVC = C
      CSQR = DVC*DVC
C
C refinment factor for Runge-Kutta-starter
C
      NDIV = 60
C
C determine nearest cut point lower than nmatch,nzero
C
      IPANMATCH = NPAN(IM)
      DO WHILE ( JRCUT(IPANMATCH,IM).GE.NMATCH )
         IPANMATCH = IPANMATCH - 1
      END DO
C
      IPANNZERO = NPAN(IM)
      DO WHILE ( JRCUT(IPANNZERO,IM).GE.NZERO )
         IPANNZERO = IPANNZERO - 1
      END DO
C
C     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
C
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         TZ = DBLE(NINT(-VV(1)*RC(1)))
         VC(0) = VV(1) - (-TZ)/RC(1)
      ELSE
         TZ = 2.0D0*DBLE(Z)
         VC(0) = VV(1)
      END IF
      DO I = 1,2
         DO J = 1,2
            DO K = 1,NPEMAX
               PC(I,J,K) = 0.0D0
               QC(I,J,K) = 0.0D0
            END DO
         END DO
      END DO
C
      BC(0) = BB(1)
C
C    CALCULATE G-COEFFICIENTS OF B-FIELD
C
      KAP1 = -L - 1
      KAP2 = +L
C
      CG1 = -MJ/(KAP1+0.5D0)
      CG5 = -MJ/(-KAP1+0.5D0)
      CGD(1) = CG1
      CGMD(1) = CG5
      KAP(1) = DBLE(KAP1)
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         GAM(1) = DSQRT(KAP(1)**2-(TZ/DVC)**2)
      ELSE
         GAM(1) = DABS(KAP(1))
      END IF
C
      IF ( DABS(MJ).GT.L ) THEN
         CG2 = 0.0D0
         CG4 = 0.0D0
         CG8 = 0.0D0
         NSOL = 1
         CGD(2) = 0.0D0
         CGO = 0.0D0
         CGMD(2) = 0.0D0
         GAM(2) = 0.0D0
         KAP(2) = 0.0D0
      ELSE
         CG2 = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CG4 = -MJ/(KAP2+0.5D0)
         CG8 = -MJ/(-KAP2+0.5D0)
         NSOL = 2
         CGD(2) = CG4
         CGO = CG2
         CGMD(2) = CG8
         KAP(2) = DBLE(KAP2)
         IF ( .NOT.FINITE_NUCLEUS ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/DVC)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
      END IF
C
C
      IF ( WAY.EQ.'INW' ) THEN
C
C #####################################################################
C
C             INWARD INTEGRATION
C
         DMUE = SQRT(-E-E*E/CSQR)
         BOVA = -DMUE/(1.0D0+E/CSQR)
C
         DO N = (NZERO-3),NZERO
C
            RR = RC(N)
C
            DO J = 1,NSOL
               I = 3 - J
               WP(J,J,N) = DEXP(-DMUE*RR)
               DP(J,J,N) = -DMUE*DRDIC(N)*WP(J,J,N)
               WQ(J,J,N) = BOVA*WP(J,J,N)
               DQ(J,J,N) = BOVA*DP(J,J,N)
C
               WP(I,J,N) = 0.0D0
               WQ(I,J,N) = 0.0D0
               DP(I,J,N) = 0.0D0
               DQ(I,J,N) = 0.0D0
            END DO
         END DO
C
C =============================================================== N ====
C
C        ------------------------------------------------------------
C              initialize inward integration with runge - kutta
C        ------------------------------------------------------------
C
         DO IPAN = IPANNZERO + 1,IPANMATCH + 1, - 1
C
            IF ( IPAN.LT.IPANNZERO+1 ) THEN
C
               SRK = 1.0D0/DBLE(NDIV)
               SO2 = SRK/2.0D0
               SO6 = SRK/6.0D0
C
               NMESHLC = JRCUT(IPAN,IM)
               N = NMESHLC
C
               IF ( IPAN.LT.NPAN(IM) ) THEN
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        WP(I,J,N) = WP(I,J,N+1)
                        WQ(I,J,N) = WQ(I,J,N+1)
                     END DO
                  END DO
               END IF
C
               EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
               EMVPP = -(E-VV(N))*DRDIC(N)
               BQQ = BB(N)*DRDIC(N)/CSQR
               BPP = BB(N)*DRDIC(N)
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     DP(I,J,N) = -KAP(I)*WP(I,J,N)*DOVRC(N)
     &                           + (EMVQQ+BQQ*CGMD(I))*WQ(I,J,N)
C
                     DQ(I,J,N) = KAP(I)*WQ(I,J,N)*DOVRC(N)
     &                           + (EMVPP+BPP*CGD(I))*WP(I,J,N)
     &                           + BPP*CGO*WP(3-I,J,N)
C
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P1(I,J) = WP(I,J,N)
                     Q1(I,J) = WQ(I,J,N)
                     MP1(I,J) = DP(I,J,N)
                     MQ1(I,J) = DQ(I,J,N)
                  END DO
               END DO
C
               X14 = DBLE(N)
               NHLP = NABM + 4
               HLP1 = DBLE(NMESHLC-NHLP)
               DO I = 1,NHLP
                  HLP(I) = DBLE(I)
                  VHLP(I) = VV(NMESHLC-NHLP+I)
                  BHLP(I) = BB(NMESHLC-NHLP+I)
                  DHLP(I) = DRDIC(NMESHLC-NHLP+I)
                  RHLP(I) = RC(NMESHLC-NHLP+I)
               END DO
               DO IRK = 1,(NABM-1)*NDIV
C
                  XH = X14 - SO2
                  VH = YLAG(XH-HLP1,HLP,VHLP,0,3,NHLP)
                  BH = YLAG(XH-HLP1,HLP,BHLP,0,3,NHLP)
                  RH = YLAG(XH-HLP1,HLP,RHLP,0,3,NHLP)
                  DH = YLAG(XH-HLP1,HLP,DHLP,0,3,NHLP)
                  EMVQQ = (E-VH+CSQR)*DH/CSQR
                  EMVPP = -(E-VH)*DH
                  BQQ = BH*DH/CSQR
                  BPP = BH*DH
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        P2(I,J) = P1(I,J) - SO2*MP1(I,J)
                        Q2(I,J) = Q1(I,J) - SO2*MQ1(I,J)
                     END DO
                  END DO
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        MP2(I,J) = -KAP(I)*P2(I,J)
     &                             *DH/RH + (EMVQQ+BQQ*CGMD(I))*Q2(I,J)
C
                        MQ2(I,J) = KAP(I)*Q2(I,J)
     &                             *DH/RH + (EMVPP+BPP*CGD(I))*P2(I,J)
     &                             + BPP*CGO*P2(3-I,J)
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        P3(I,J) = P1(I,J) - SO2*MP2(I,J)
                        Q3(I,J) = Q1(I,J) - SO2*MQ2(I,J)
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        MP3(I,J) = -KAP(I)*P3(I,J)
     &                             *DH/RH + (EMVQQ+BQQ*CGMD(I))*Q3(I,J)
C
                        MQ3(I,J) = KAP(I)*Q3(I,J)
     &                             *DH/RH + (EMVPP+BPP*CGD(I))*P3(I,J)
     &                             + BPP*CGO*P3(3-I,J)
                     END DO
                  END DO
C
                  X14 = X14 - SRK
                  V14 = YLAG(X14-HLP1,HLP,VHLP,0,3,NHLP)
                  B14 = YLAG(X14-HLP1,HLP,BHLP,0,3,NHLP)
                  R14 = YLAG(X14-HLP1,HLP,RHLP,0,3,NHLP)
                  D14 = YLAG(X14-HLP1,HLP,DHLP,0,3,NHLP)
C
                  EMVQQ = (E-V14+CSQR)*D14/CSQR
                  EMVPP = -(E-V14)*D14
                  BQQ = B14*D14/CSQR
                  BPP = B14*D14
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        P4(I,J) = P1(I,J) - SRK*MP3(I,J)
                        Q4(I,J) = Q1(I,J) - SRK*MQ3(I,J)
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        MP4(I,J) = -KAP(I)*P4(I,J)*D14/R14 + 
     &                             (EMVQQ+BQQ*CGMD(I))*Q4(I,J)
C
                        MQ4(I,J) = KAP(I)*Q4(I,J)*D14/R14 + 
     &                             (EMVPP+BPP*CGD(I))*P4(I,J)
     &                             + BPP*CGO*P4(3-I,J)
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        P1(I,J) = P1(I,J)
     &                            - SO6*(MP1(I,J)+2*(MP2(I,J)+MP3(I,J))
     &                            +MP4(I,J))
                        Q1(I,J) = Q1(I,J)
     &                            - SO6*(MQ1(I,J)+2*(MQ2(I,J)+MQ3(I,J))
     &                            +MQ4(I,J))
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        MP1(I,J) = -KAP(I)*P1(I,J)*D14/R14 + 
     &                             (EMVQQ+BQQ*CGMD(I))*Q1(I,J)
C
                        MQ1(I,J) = KAP(I)*Q1(I,J)*D14/R14 + 
     &                             (EMVPP+BPP*CGD(I))*P1(I,J)
     &                             + BPP*CGO*P1(3-I,J)
                     END DO
                  END DO
C
                  IF ( MOD(IRK,NDIV).EQ.0 ) THEN
                     N = NMESHLC - IRK/NDIV
                     IF ( ABS(X14-DBLE(N)).GT.1.0D-5 ) THEN
                        WRITE (*,*) ' <dirac> runge-kutta: ',IRK,NDIV,N,
     &                              X14
                        STOP
                     END IF
                     DO J = 1,NSOL
                        DO I = 1,NSOL
                           WP(I,J,N) = P1(I,J)
                           WQ(I,J,N) = Q1(I,J)
                           DP(I,J,N) = MP1(I,J)
                           DQ(I,J,N) = MQ1(I,J)
                        END DO
                     END DO
C
                  END IF
C
               END DO
C
            END IF
C
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
C
C      DO 540 NN = 1,(NZERO-3 -NMATCH)
C
            IF ( IPAN.EQ.NPAN(IM)+1 ) THEN
               NSTART = NZERO
               NEND = MAX(NMATCH,JRCUT(NPAN(IM),IM))
            ELSE
               NSTART = MIN(NZERO,JRCUT(IPAN,IM))
               NEND = MAX(NMATCH,JRCUT(IPAN-1,IM)+1)
            END IF
C
C
            DO NN = 1,(NSTART-3-NEND)
C
               N = NSTART - 3 - NN
C
C    EVALUATE PREDICTOR
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     PNEW(I,J) = WP(I,J,N+1)
     &                           - H24*(55.0D0*DP(I,J,N+1)-59.0D0*DP(I,
     &                           J,N+2)+37.0D0*DP(I,J,N+3)
     &                           -9.0D0*DP(I,J,N+4))
                     QNEW(I,J) = WQ(I,J,N+1)
     &                           - H24*(55.0D0*DQ(I,J,N+1)-59.0D0*DQ(I,
     &                           J,N+2)+37.0D0*DQ(I,J,N+3)
     &                           -9.0D0*DQ(I,J,N+4))
                  END DO
               END DO
C
               EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
               EMVPP = -(E-VV(N))*DRDIC(N)
               BQQ = BB(N)*DRDIC(N)/CSQR
               BPP = BB(N)*DRDIC(N)
C
C    EVALUATE CORRECTOR
C
               DO JCORR = 1,ITMAX
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        POLD(I,J) = PNEW(I,J)
                        QOLD(I,J) = QNEW(I,J)
                        DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                              + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                        DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                              + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                              + BPP*CGO*PNEW(3-I,J)
C
                        PNEW(I,J) = WP(I,J,N+1)
     &                              - H24*(9.0D0*DP(I,J,N)+19.0D0*DP(I,
     &                              J,N+1)-5.0D0*DP(I,J,N+2)+DP(I,J,N+3)
     &                              )
                        QNEW(I,J) = WQ(I,J,N+1)
     &                              - H24*(9.0D0*DQ(I,J,N)+19.0D0*DQ(I,
     &                              J,N+1)-5.0D0*DQ(I,J,N+2)+DQ(I,J,N+3)
     &                              )
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        DIFFA = POLD(I,J) - PNEW(I,J)
                        IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) )
     &                       GOTO 10
                        IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) )
     &                       GOTO 10
C
                        DIFFB = QOLD(I,J) - QNEW(I,J)
                        IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) )
     &                       GOTO 10
                        IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) )
     &                       GOTO 10
                     END DO
                  END DO
                  GOTO 20
C
 10            END DO
               WRITE (6,99001) KAP1,N,RC(N),DIFFA,DIFFB,IT,L,INT(2*MJ),
     &                         'OUT'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
C
 20            CONTINUE
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     WP(I,J,N) = PNEW(I,J)
                     WQ(I,J,N) = QNEW(I,J)
                     DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                           + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                     DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                           + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                           + BPP*CGO*PNEW(3-I,J)
                  END DO
               END DO
C
            END DO
C
         END DO
C              ! ipan - Schleife
C
C =============================================================== N ====
C
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
         DO N = NMATCH,NZERO
            DO J = 1,NSOL
               DO I = 1,NSOL
                  GC(I,J,N) = WP(I,J,N)/RC(N)
                  FC(I,J,N) = WQ(I,J,N)/(RC(N)*DVC)
               END DO
            END DO
         END DO
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PIW(I,J) = WP(I,J,NMATCH)
               QIW(I,J) = WQ(I,J,NMATCH)
            END DO
         END DO
         GOTO 99999
      END IF
C
C #####################################################################
C
C             OUTWARD INTEGRATION
C
C  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
C
      MPS = 20
C
      AA12 = -TZ/CSQR
      AA21 = TZ
      EMVQQ = (E-VC(0)+CSQR)/CSQR
      EMVPP = -E + VC(0)
      BQQ = BC(0)/CSQR
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KAP(J))-GAM(J))
            QC(J,J,0) = (KAP(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = 0.0D0
            QC(I,J,0) = 0.0D0
         END DO
C
         DO J = 1,NSOL
C
            DO M = 1,MPS
               DO I = 1,NSOL
                  BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                  BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                  *CGO*PC(3-I,J,M-1)
                  AA11 = GAM(J) + M + KAP(I)
                  AA22 = GAM(J) + M - KAP(I)
                  DET = AA11*AA22 - AA12*AA21
                  PC(I,J,M) = (BB1*AA22-AA12*BB2)/DET
                  QC(I,J,M) = (AA11*BB2-BB1*AA21)/DET
               END DO
            END DO
C
         END DO
C
      ELSE
C EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
C EXPANSION OF POTENTIAL actually UP TO zeroth ORDER
C
C       DO IV=1,INVMAX
C        DO N=1,INVMAX
C         CM(N,IV)=RC(N)**(IV-1)
C        ENDDO
C       ENDDO
C
C       CALL RINVGJ(CMI,CM,INVMAX,INVMAX)
C
         DO IV = 1,INVMAX
            VC(IV-1) = 0.0D0
C        DO N=1,INVMAX
C         VC(IV-1)=VC(IV-1)+CMI(IV,N)*VV(N)
C        ENDDO
         END DO
         DO J = 1,NSOL
            I = 3 - J
            IF ( KAP(J).GT.0 ) THEN
C ARBITRARY STARTING VALUES
               ALPHA = 0.0D0
               BETA = 0.174D0
            ELSE
               BETA = 0.0D0
               ALPHA = 0.174D0
            END IF
            PC(J,J,0) = ALPHA
            QC(J,J,0) = BETA
            PC(I,J,0) = 0.0D0
            QC(I,J,0) = 0.0D0
         END DO
C
         W4 = BC(0)*CGO
         W2 = VC(1)/CSQR
         W5 = VC(1)
         W6 = VC(2)/CSQR
         W7 = VC(2)
         DO J = 1,NSOL
            DO I = 1,NSOL
               W1 = EMVQQ + BQQ*CGMD(I)
               W3 = -EMVPP + BC(0)*CGD(I)
               A11 = GAM(J) + KAP(I) + 1D0
               A12 = GAM(J) - KAP(I) + 1D0
               IF ( RNON0(A11) ) PC(I,J,1) = W1/A11*QC(I,J,0)
               IF ( RNON0(A12) ) QC(I,J,1)
     &              = (-W3*PC(I,J,0)+W4*PC(3-I,J,0))/A12
C
            END DO
         END DO
         DO J = 1,NSOL
            DO I = 1,NSOL
               W1 = EMVQQ + BQQ*CGMD(I)
               W3 = -EMVPP + BC(0)*CGD(I)
               A11 = GAM(J) + KAP(I) + 2D0
               A12 = GAM(J) - KAP(I) + 2D0
               IF ( RNON0(A11) ) PC(I,J,2) = (W1*QC(I,J,1)-W2*QC(I,J,0))
     &              /A11
               IF ( RNON0(A12) ) QC(I,J,2)
     &              = (-W3*PC(I,J,1)+W4*PC(3-I,J,1)+W5*PC(I,J,0))/A12
            END DO
         END DO
         DO J = 1,NSOL
            DO M = 3,MPS
               DO I = 1,NSOL
                  W1 = EMVQQ + BQQ*CGMD(I)
                  W3 = -EMVPP + BC(0)*CGD(I)
                  A21 = GAM(J) + KAP(I) + DBLE(M)
                  A22 = GAM(J) - KAP(I) + DBLE(M)
                  IF ( RNON0(A21) ) PC(I,J,M)
     &                 = (W1*QC(I,J,M-1)-W2*QC(I,J,M-2)-W6*QC(I,J,M-3))
     &                 /A21
                  IF ( RNON0(A22) ) QC(I,J,M)
     &                 = (-W3*PC(I,J,M-1)+W4*PC(3-I,J,M-1)
     &                 +W5*PC(I,J,M-2)+W7*PC(I,J,M-3))/A22
               END DO
            END DO
         END DO
      END IF
C
C  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
C  FOR THE FIRST 4 R - MESH - POINTS
C
      DO N = 1,4
         RR = RC(N)
C
         DO J = 1,NSOL
            RPWGPM = RR**GAM(J)
C
            DO I = 1,NSOL
               WP(I,J,N) = PC(I,J,0)*RPWGPM
               WQ(I,J,N) = QC(I,J,0)*RPWGPM
               DP(I,J,N) = PC(I,J,0)*RPWGPM*GAM(J)*DOVRC(N)
               DQ(I,J,N) = QC(I,J,0)*RPWGPM*GAM(J)*DOVRC(N)
            END DO
C
            DO M = 1,MPS
               RPWGPM = RPWGPM*RR
               GPM = GAM(J) + M
C
               DO I = 1,NSOL
                  WP(I,J,N) = WP(I,J,N) + PC(I,J,M)*RPWGPM
                  WQ(I,J,N) = WQ(I,J,N) + QC(I,J,M)*RPWGPM
                  DP(I,J,N) = DP(I,J,N) + PC(I,J,M)*RPWGPM*GPM*DOVRC(N)
                  DQ(I,J,N) = DQ(I,J,N) + QC(I,J,M)*RPWGPM*GPM*DOVRC(N)
               END DO
C
            END DO
         END DO
      END DO
C
C        ------------------------------------------------------------
C              initialize outward integration with runge - kutta
C        ------------------------------------------------------------
C
C
      DO IPAN = 0,IPANMATCH
C
         IF ( IPAN.GT.0 ) THEN
C
            SRK = 1.0D0/DBLE(NDIV)
            SO2 = SRK/2.0D0
            SO6 = SRK/6.0D0
C
            IF ( IPAN.LT.NPAN(IM) ) THEN
               N = JRCUT(IPAN,IM) + 1
            ELSE
               N = JRCUT(IPAN,IM)
            END IF
C
            NSTART = N
C
            IF ( IPAN.LT.NPAN(IM) ) THEN
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     WP(I,J,N) = WP(I,J,N-1)
                     WQ(I,J,N) = WQ(I,J,N-1)
                  END DO
               END DO
            END IF
C
            EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
            EMVPP = -(E-VV(N))*DRDIC(N)
            BQQ = BB(N)*DRDIC(N)/CSQR
            BPP = BB(N)*DRDIC(N)
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  DP(I,J,N) = -KAP(I)*WP(I,J,N)*DOVRC(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*WQ(I,J,N)
C
                  DQ(I,J,N) = KAP(I)*WQ(I,J,N)*DOVRC(N)
     &                        + (EMVPP+BPP*CGD(I))*WP(I,J,N)
     &                        + BPP*CGO*WP(3-I,J,N)
C
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P1(I,J) = WP(I,J,N)
                  Q1(I,J) = WQ(I,J,N)
                  MP1(I,J) = DP(I,J,N)
                  MQ1(I,J) = DQ(I,J,N)
               END DO
            END DO
C
            X14 = DBLE(N)
            NHLP = NABM + 4
            HLP1 = DBLE(N) - 1
            DO I = 1,NHLP
               HLP(I) = DBLE(I)
               VHLP(I) = VV(N+I-1)
               BHLP(I) = BB(N+I-1)
               DHLP(I) = DRDIC(N+I-1)
               RHLP(I) = RC(N+I-1)
            END DO
            DO IRK = 1,(NABM-1)*NDIV
C
               XH = X14 + SO2
               VH = YLAG(XH-HLP1,HLP,VHLP,0,3,NHLP)
               BH = YLAG(XH-HLP1,HLP,BHLP,0,3,NHLP)
               RH = YLAG(XH-HLP1,HLP,RHLP,0,3,NHLP)
               DH = YLAG(XH-HLP1,HLP,DHLP,0,3,NHLP)
               EMVQQ = (E-VH+CSQR)*DH/CSQR
               EMVPP = -(E-VH)*DH
               BQQ = BH*DH/CSQR
               BPP = BH*DH
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P2(I,J) = P1(I,J) + SO2*MP1(I,J)
                     Q2(I,J) = Q1(I,J) + SO2*MQ1(I,J)
                  END DO
               END DO
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     MP2(I,J) = -KAP(I)*P2(I,J)
     &                          *DH/RH + (EMVQQ+BQQ*CGMD(I))*Q2(I,J)
C
                     MQ2(I,J) = KAP(I)*Q2(I,J)
     &                          *DH/RH + (EMVPP+BPP*CGD(I))*P2(I,J)
     &                          + BPP*CGO*P2(3-I,J)
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P3(I,J) = P1(I,J) + SO2*MP2(I,J)
                     Q3(I,J) = Q1(I,J) + SO2*MQ2(I,J)
                  END DO
               END DO
C
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     MP3(I,J) = -KAP(I)*P3(I,J)
     &                          *DH/RH + (EMVQQ+BQQ*CGMD(I))*Q3(I,J)
C
                     MQ3(I,J) = KAP(I)*Q3(I,J)
     &                          *DH/RH + (EMVPP+BPP*CGD(I))*P3(I,J)
     &                          + BPP*CGO*P3(3-I,J)
                  END DO
               END DO
C
               X14 = X14 + SRK
               V14 = YLAG(X14-HLP1,HLP,VHLP,0,3,NHLP)
               B14 = YLAG(X14-HLP1,HLP,BHLP,0,3,NHLP)
               R14 = YLAG(X14-HLP1,HLP,RHLP,0,3,NHLP)
               D14 = YLAG(X14-HLP1,HLP,DHLP,0,3,NHLP)
C
               EMVQQ = (E-V14+CSQR)*D14/CSQR
               EMVPP = -(E-V14)*D14
               BQQ = B14*D14/CSQR
               BPP = B14*D14
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P4(I,J) = P1(I,J) + SRK*MP3(I,J)
                     Q4(I,J) = Q1(I,J) + SRK*MQ3(I,J)
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     MP4(I,J) = -KAP(I)*P4(I,J)*D14/R14 + 
     &                          (EMVQQ+BQQ*CGMD(I))*Q4(I,J)
C
                     MQ4(I,J) = KAP(I)*Q4(I,J)*D14/R14 + 
     &                          (EMVPP+BPP*CGD(I))*P4(I,J)
     &                          + BPP*CGO*P4(3-I,J)
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P1(I,J) = P1(I,J)
     &                         + SO6*(MP1(I,J)+2*(MP2(I,J)+MP3(I,J))
     &                         +MP4(I,J))
                     Q1(I,J) = Q1(I,J)
     &                         + SO6*(MQ1(I,J)+2*(MQ2(I,J)+MQ3(I,J))
     &                         +MQ4(I,J))
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     MP1(I,J) = -KAP(I)*P1(I,J)*D14/R14 + 
     &                          (EMVQQ+BQQ*CGMD(I))*Q1(I,J)
C
                     MQ1(I,J) = KAP(I)*Q1(I,J)*D14/R14 + 
     &                          (EMVPP+BPP*CGD(I))*P1(I,J)
     &                          + BPP*CGO*P1(3-I,J)
                  END DO
               END DO
C
               IF ( MOD(IRK,NDIV).EQ.0 ) THEN
                  N = NSTART + IRK/NDIV
                  IF ( ABS(X14-DBLE(N)).GT.1.0D-5 ) THEN
                     WRITE (*,*) ' <dirac> runge-kutta: ',IRK,NDIV,N,X14
                     STOP
                  END IF
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        WP(I,J,N) = P1(I,J)
                        WQ(I,J,N) = Q1(I,J)
                        DP(I,J,N) = MP1(I,J)
                        DQ(I,J,N) = MQ1(I,J)
                     END DO
                  END DO
C
               END IF
C
            END DO
C
         END IF
C
C
C =============================================================== N ====
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
C
         NSTART = JRCUT(IPAN,IM) + 1 + NABM
         NEND = MIN(NMATCH,JRCUT(IPAN+1,IM))
C
         DO N = NSTART,NEND
C
C
C    EVALUATE PREDICTOR
C
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PNEW(I,J) = WP(I,J,N-1)
     &                        + H24*(55.0D0*DP(I,J,N-1)-59.0D0*DP(I,J,
     &                        N-2)+37.0D0*DP(I,J,N-3)-9.0D0*DP(I,J,N-4))
                  QNEW(I,J) = WQ(I,J,N-1)
     &                        + H24*(55.0D0*DQ(I,J,N-1)-59.0D0*DQ(I,J,
     &                        N-2)+37.0D0*DQ(I,J,N-3)-9.0D0*DQ(I,J,N-4))
               END DO
            END DO
C
            EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
            EMVPP = -(E-VV(N))*DRDIC(N)
            BQQ = BB(N)*DRDIC(N)/CSQR
            BPP = BB(N)*DRDIC(N)
C
C    EVALUATE CORRECTOR
C
C
            DO JCORR = 1,ITMAX
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     POLD(I,J) = PNEW(I,J)
                     QOLD(I,J) = QNEW(I,J)
                     DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                           + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                     DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                           + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                           + BPP*CGO*PNEW(3-I,J)
C
                     PNEW(I,J) = WP(I,J,N-1)
     &                           + H24*(9.0D0*DP(I,J,N)+19.0D0*DP(I,J,
     &                           N-1)-5.0D0*DP(I,J,N-2)+DP(I,J,N-3))
                     QNEW(I,J) = WQ(I,J,N-1)
     &                           + H24*(9.0D0*DQ(I,J,N)+19.0D0*DQ(I,J,
     &                           N-1)-5.0D0*DQ(I,J,N-2)+DQ(I,J,N-3))
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     DIFFA = POLD(I,J) - PNEW(I,J)
                     IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 40
                     IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 40
C
                     DIFFB = QOLD(I,J) - QNEW(I,J)
                     IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 40
                     IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 40
                  END DO
               END DO
               GOTO 60
C
 40         END DO
            WRITE (6,99001) KAP1,N,RC(N),DIFFA,DIFFB,IT,L,INT(2*MJ),
     &                      'OUT'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
 60         CONTINUE
            DO J = 1,NSOL
               DO I = 1,NSOL
                  WP(I,J,N) = PNEW(I,J)
                  WQ(I,J,N) = QNEW(I,J)
                  DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                        + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                        + BPP*CGO*PNEW(3-I,J)
               END DO
            END DO
C
         END DO
C
      END DO
C                                                            ipan - loop
C =============================================================== N ====
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
C
      DO N = 1,NMATCH
         DO J = 1,NSOL
            DO I = 1,NSOL
               GC(I,J,N) = WP(I,J,N)/RC(N)
               FC(I,J,N) = WQ(I,J,N)/(RC(N)*DVC)
            END DO
         END DO
      END DO
C
      DO J = 1,NSOL
         DO I = 1,NSOL
            POW(I,J) = WP(I,J,NMATCH)
            QOW(I,J) = WQ(I,J,NMATCH)
         END DO
      END DO
C
      RETURN
C
99001 FORMAT (' P/C NOT CONV. IN <DIRAC> ',2I4,2X,F10.7,2X,2E12.4,3I2,
     &        '/2 ',A3)
99999 CONTINUE
      END
C*==fpcoresrahff.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPCORESRAHFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,RC,
     &                        DRDIC,RNUC,NZEROIN,NRC,JRCUT,NPAN,NPANMAX)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
C   *                  CURRENT  CORE STATE S                           *
C   *                                                                  *
C   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:B_AU2CGS
      IMPLICIT NONE
C*--FPCORESRAHFF2130
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,KAP1,KAP2,NPAN,NPANMAX,NRC,NSOL,NZEROIN,S
      REAL*8 MJ,RNUC
      REAL*8 BHF(2,2),DRDIC(NRC),FCK(2,2,NRC),GCK(2,2,NRC),RC(NRC)
      INTEGER JRCUT(0:NPANMAX)
C
C Local variables
C
      REAL*8 AMEHF(2,2),XX(5),YI(:),YY(5),ZI(:)
      INTEGER I,K1,K2,N,NZERO
      LOGICAL RNON0
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE YI,ZI
      ALLOCATE (YI(NRC),ZI(NRC))
C
C CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
C                                      ELECTRON CHARGE     IN ESU
C                                      BOHR-RADIUS         IN CM
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
      AMEHF(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
      IF ( NSOL.EQ.2 ) THEN
         AMEHF(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
         AMEHF(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
         AMEHF(2,1) = AMEHF(1,2)
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NZERO = NZEROIN
C      nzero=min(nzeroin,JRCUT(NPAN))
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
            DO N = 1,NZERO
               YI(N) = DRDIC(N)
     &                 *(GCK(K1,S,N)*FCK(K2,S,N)+FCK(K1,S,N)*GCK(K2,S,N)
     &                 )
            END DO
C
            DO N = NZERO + 1,JRCUT(NPAN)
               YI(N) = 0D0
            END DO
C
            CALL RRADINT_R(IM,YI,ZI)
C
            IF ( RNON0(RNUC) ) THEN
               DO I = 1,5
                  XX(I) = RC(NZERO-5+I)
                  YY(I) = ZI(NZERO-5+I)
               END DO
               ZI(NZERO) = YLAG(RNUC,XX,YY,0,4,5)
            END IF
            XX(1) = 1.0D0 !RC( 1)
            XX(2) = 6.0D0 !RC( 6)
            XX(3) = 11.0D0
C                         !RC(11)
            YY(1) = ZI(NZERO) - ZI(1)
            YY(2) = ZI(NZERO) - ZI(6)
            YY(3) = ZI(NZERO) - ZI(11)
            BHF(K1,K2) = B_AU2CGS*AMEHF(K1,K2)*YLAG(0.0D0,XX,YY,0,2,3)
         END DO
      END DO
C
      END
