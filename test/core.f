C*==core.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SUBROUTINE TO CALCULATE THE RELATIVISTIC CORE WAVE             *
C   *   FUNCTIONS FOR A SPIN-DEPENDENT POTENTIAL                       *
C   *                                                                  *
C   *   FOR A GIVEN POTENTIAL THE NUMBER OF CORE AND VALENCE           *
C   *   ELECTRONS IS DETERMINED AND ALL CORE STATES THEN CALCULATED    *
C   *   > THE ROUTINE IS ORGANIZED AS DESCLAUX'S ROUTINE <RESLD>       *
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
C   *   HYPERFINE FIELD SPLITTING introduced if icore=1 MB JUN. 1995   *
C   *                                                                  *
C   *   SCALEB:                                                        *
C   *   if the B-field is quite high it might happen that the routine  *
C   *   fails to find both 'spin-orbit-split' solutions.               *
C   *   in that case the whole l-shell is rerun with the B-field       *
C   *   gradually switched on, i.e. scaled with a parameter that       *
C   *   increases from 0 to 1 during the iteration loop  HE Nov. 95    *
C   *                                                                  *
C   *                                                                  *
C   *   ITXRAY controls mode of operation !                            *
C   *                                                                  *
C   *   ITXRAY =  0  run over all core states to get charge density    *
C   *   ITXRAY >  0  deal exclusively with state  NCXRAY,LCXRAY        *
C   *   ITXRAY <  0  state  NCXRAY,LCXRAY  is checked to be a          *
C   *                bound state or not. on return:                    *
C   *                ITXRAY = |ITXRAY| indicates bound state           *
C   *                ITXRAY = -1       indicates NO bound state found  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_CALCMODE,ONLY:ISMQHFI,NONMAG,TASK
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,NRCMAX,FINITE_NUCLEUS,
     &    R2DRDI_W_RADINT,FULLPOT,NRSFTOT,FLMSF,JRCUT,ISFLM,R_COR,
     &    DRDI_COR,DRDIOVR_COR
      USE MOD_TYPES,ONLY:NTMAX,LCXRAY,NCXRAY,ECORTAB,ECOR_LT,RHOSPNC,
     &    RHOCHRC,IMT,Z,BT,VT,CTL,NCORT,ITBOT,ITTOP
      USE MOD_ANGMOM,ONLY:NMEHFMAX,NLMAX,NKMMAX,QMOFF,QOFF,QMDIA,QDIA,
     &    SMOFF,SOFF,SMDIA,SDIA,NLCORE
      IMPLICIT NONE
C*--CORE60
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CORE')
      REAL*8 UNEND,TRYMIX,DVSTEP,TOL_RELCHANGE,TOLVAR
      PARAMETER (UNEND=600.0D0,TRYMIX=0.01D0,DVSTEP=0.01D0,
     &           TOL_RELCHANGE=0.1D0,TOLVAR=1.0D-8)
      INTEGER ITERMAX,NLSHELLMAX
      PARAMETER (ITERMAX=400,NLSHELLMAX=15)
      LOGICAL CHECK_AGAINST_OLD_ITERATION_SCHEME
      PARAMETER (CHECK_AGAINST_OLD_ITERATION_SCHEME=.FALSE.)
C
C Dummy arguments
C
      INTEGER IPRINT,ITXRAY,NCSTMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
C
C Local variables
C
      REAL*8 AUX,BB(:),BEX,BHF(2,2),BHF1(2,2),BHF2(2,2),BSH,BSOL,BSUM,
     &       CGD(2),CGMD(2),CGO,DEC,DEDV(4,4),DOVRC(:),DP(:,:,:),
     &       DQ(:,:,:),DRDIC(:),DROVRN(:),DV(4),DVDE(4,4),EC,ECC,
     &       ECORDUM(:),ELIM,ERR(4),ERRNEW(4),FC(:,:,:),FCK(:,:,:),
     &       GC(:,:,:),GCK(:,:,:),INTEG_CHR,INTEG_SPN,MAX_RELCHANGE,MJ,
     &       NIW(2),NORM,NOW(2),PIW(2,2),POW(2,2),QIW(2,2),QOW(2,2),
     &       RATT,RC(:),RELCHANGE,RMAT4(4,4),RNUC,RR,SCLFAC,SHF(:,:,:),
     &       SPLIT(:),SPLIT1(:),SPLIT2(:,:),SPLIT3(:,:),SZ,SZ1,
     &       SZCORDUM(:),SZD,VAL,VAR(4),VARNEW(4),VARTAB(4,20),VV(:),VZ,
     &       W,WP(:,:,:),WQ(:,:,:),WR
      LOGICAL BNDSTACHK,FERRO,SCALEB,SUPPRESSB
      INTEGER I,IA_ERR,IC,IC1,IC2,ICOR,ICST,ICSTDUM,IE,IFLAG,II,IL,ILC,
     &        ILSHELL,IM,IMEHF,IMIN,IN,IPIV4(4),IR,IRSF,IRTOP,
     &        IRTOP_SPHERE,IR_MATCH,IR_ZERO,IR_ZERO_TRIED,ISF_00,ISH,
     &        ISTART,IT,ITER,IV,J,JLIM,JV,K,KAP(2),KAP1,KAP2,KC,L,LLL,
     &        LOOP,LQNTAB(NLSHELLMAX),MUEM05,NLSHELL,NN,NODE,NQN,
     &        NQNTAB(NLSHELLMAX),NRC,NSH,NSOL,NVAR,S,T
      INTEGER IKAPMUE
      REAL*8 RNUCTAB
      CHARACTER*10 TXTB(1:5)
      CHARACTER*3 TXTK(4)
      CHARACTER*1 TXT_L(0:3)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DRDIC,DOVRC,SPLIT,BB,FC,GC,SPLIT1,SPLIT2,SPLIT3
      ALLOCATABLE DP,DQ,RC,WP,WQ,VV,DROVRN,FCK,GCK,SHF
      ALLOCATABLE SZCORDUM,ECORDUM
C
C ----------------------------------------------------------------------
C  the core states are assumed to have only  s,p,d,f
C  but no higher l-components, so  NLCORE=4
C
      DATA TXTB/'B_nses','B_nseo','B_ssc ','B_sum ','B_tot '/
      DATA TXT_L/'s','p','d','f'/
      DATA TXTK/'1/2','3/2','5/2','7/2'/
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
C
      ALLOCATE (RC(NRCMAX),DRDIC(NRMAX*3))
      ALLOCATE (DOVRC(NRCMAX))
      ALLOCATE (FCK(2,2,NRMAX*2),GCK(2,2,NRMAX*2))
      ALLOCATE (DP(2,2,NRCMAX),DQ(2,2,NRCMAX))
      ALLOCATE (WP(2,2,NRMAX*2),WQ(2,2,NRMAX*2))
      ALLOCATE (DROVRN(2*NRMAX),VV(NRMAX*3),BB(NRMAX*3))
      ALLOCATE (GC(2,2,NRCMAX),FC(2,2,NRCMAX))
      ALLOCATE (SPLIT(NMEHFMAX),SPLIT1(NMEHFMAX),SPLIT2(NMEHFMAX,NTMAX))
      ALLOCATE (SPLIT3(NMEHFMAX,NTMAX),STAT=IA_ERR)
      ALLOCATE (SZCORDUM(200),ECORDUM(200))
      IF ( .NOT.ALLOCATED(ECORTAB) ) ALLOCATE (ECORTAB(120,NTMAX))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RC')
C
      ECOR_LT(:,:) = 0.0D0
      ECORTAB(:,:) = 0.0D0
      VARTAB(:,:) = 0.0D0
C
      IF ( ITXRAY.LT.0 ) THEN
         ITXRAY = ABS(ITXRAY)
         BNDSTACHK = .TRUE.
      ELSE
         BNDSTACHK = .FALSE.
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         SUPPRESSB = .FALSE.
         SCALEB = .FALSE.
         SCLFAC = 1D0
C
C=======================================================================
C               create expanded radial mesh for core states
C=======================================================================
C
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
            IRTOP_SPHERE = JRCUT(1,IM)
            ISF_00 = ISFLM(1,IM)
         ELSE
            IRTOP = JRWS(IM)
            IRTOP_SPHERE = JRWS(IM)
            ISF_00 = ISFLM(1,IM)
         END IF
C
         NRC = 2*IRTOP
         IF ( NRC.GT.NRCMAX ) CALL STOP_MESSAGE(ROUTINE,'NRC > NRCMAX')
C
         DO IR = 1,NRC
            RC(IR) = R_COR(IR,IM)
            DRDIC(IR) = DRDI_COR(IR,IM)
            DOVRC(IR) = DRDIOVR_COR(IR,IM)
         END DO
C
C     integration boundary for hyperfine fields for finite nucleus
C     2 mesh points more for executing appropriate interpolation
C     to real nuclear radius RNUC
C
         IF ( FINITE_NUCLEUS .AND. Z(IT).NE.0 ) THEN
            RNUC = RNUCTAB(Z(IT))
            IF ( ABS(RNUC).LT.1.D-18 ) THEN
               WRITE (6,*) 'FINITE NUCLEUS: IT, Z(IT) =',IT,Z(IT)
               WRITE (6,*) 'RNUCTAB(Z(IT)) =',RNUCTAB(Z(IT))
               STOP 'RNUC is zero !'
            END IF
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
 50      CONTINUE
         BCOR(IT) = 0.0D0
         BCORS(IT) = 0.0D0
         DO I = 1,NMEHFMAX
            SPLIT2(I,IT) = 0.0D0
            SPLIT3(I,IT) = 0.0D0
         END DO
C
         DO IR = 1,NRMAX
            RHOCHRC(IR,IT) = 0.0D0
            RHOSPNC(IR,IT) = 0.0D0
         END DO
C
         DO IL = 1,NLCORE
            ECOR_LT(IL,IT) = 0.0D0
         END DO
C
C------------------------------------- cycle loop in case of empty spere
         IF ( Z(IT).LE.0 ) CYCLE LOOP_IT
C------------------------------------- cycle loop in case of empty spere
C
         BSUM = 0.0D0
         DO IR = 1,JRWS(IM)
            VV(IR) = VT(IR,IT)
            BB(IR) = BT(IR,IT)
            BSUM = BSUM + ABS(BB(IR))
         END DO
C
         IF ( SUPPRESSB ) THEN
            DO IR = 1,JRWS(IM)
               BB(IR) = 0.0D0
            END DO
            BSUM = 0.0D0
         END IF
C
         DO IR = (JRWS(IM)+1),NRC
            VV(IR) = 0.0D0
            BB(IR) = 0.0D0
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
            ICOR = 0
            DO ILSHELL = 1,NLSHELLMAX
               L = LQNTAB(ILSHELL)
               ICOR = ICOR + 2*(2*L+1)
               IF ( ICOR.EQ.NCORT(IT) ) NLSHELL = ILSHELL
            END DO
            IF ( NLSHELL.EQ.0 ) THEN
               WRITE (6,*) 'NLSHELL not found for IT=',IT,' NCORT=',
     &                     NCORT(IT)
               STOP ' in <CORE>'
            END IF
         END IF
C
         IF ( BSUM.GT.1.0D-5 ) THEN
            FERRO = .TRUE.
            IF ( NONMAG ) FERRO = .FALSE.
         ELSE
            FERRO = .FALSE.
            IF ( IPRINT.GE.1 ) WRITE (6,99001)
         END IF
C---------------------- activate next line to suppress spin-polarisation
C        FERRO = .FALSE.
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) WRITE (6,99002) IT,Z(IT)
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C                   ---------------------------------------
C                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
C                   ---------------------------------------
         ICSTDUM = 0
         IC = 0
C
         DO ILSHELL = 1,NLSHELL
            NQN = NQNTAB(ILSHELL)
            L = LQNTAB(ILSHELL)
            IL = L + 1
            ILC = MIN(NLMAX,IL)
            NSH = 2*(2*L+1)
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C          SKIP SHELL IF NOT NEEDED IN A  XRAY - CALCULATION
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            IF ( ITXRAY.NE.0 ) THEN
               IF ( IT.NE.ITXRAY ) CYCLE
               IF ( (NQN.NE.NCXRAY(IT)) .OR. (L.NE.LCXRAY(IT)) ) CYCLE
               DO ICST = 1,NCSTMAX
                  DO KC = 1,2
                     DO IR = 1,NRMAX
                        GCOR(IR,KC,ICST) = 0.0D0
                        FCOR(IR,KC,ICST) = 0.0D0
                     END DO
                  END DO
               END DO
            END IF
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
               IF ( FERRO ) THEN
                  NVAR = 2*NSOL
               ELSE
                  NVAR = 2
               END IF
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
               DO S = 1,NSOL
                  IF ( IPRINT.GE.2 ) WRITE (6,'(/,70(''-''))')
                  IR_ZERO_TRIED = 0
                  IC = IC + 1
                  ISH = ISH + 1
                  T = 3 - S
C
C                  ----------------------------------------
C                   USE EC OF PREVIOUS RUNS AS START-VALUE
C                   TAKE SPIN-ORBIT SPLITTING INTO ACCOUNT
C                  ----------------------------------------
                  IF ( ISH.GT.1 ) THEN
                     EC = ECORTAB(IC-1,IT)
                     IF ( S.EQ.2 ) EC = ECORTAB(IC-1,IT)*1.1D0
                     IF ( ISH.GE.4 ) EC = ECORTAB(IC-2,IT)
                     GOTO 65
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
                     DO IR = 2,NRC
                        VAL = VV(IR) + LLL/RC(IR)**2
                        IF ( VAL.LE.ELIM ) ELIM = VAL
                     END DO
                  END IF
C
                  EC = -DBLE(Z(IT)**2)/(2.0D0*NQN*NQN)
C
                  ISTART = 1
 55               CONTINUE
                  IF ( EC.LE.ELIM ) EC = ELIM*0.7D0
C
C                                      --------------------
C                                         FIND    NZERO
C                                      --------------------
                  DO IR = 1,(NRC-1)
                     IF ( (VV(IR)-EC)*RC(IR)**2.GT.UNEND ) THEN
                        IF ( MOD(IR,2).EQ.0 ) THEN
                           IR_ZERO = IR + 1
                        ELSE
                           IR_ZERO = IR
                        END IF
                        GOTO 60
                     END IF
                  END DO
                  IR_ZERO = NRC - 1
                  IR_ZERO_TRIED = IR_ZERO_TRIED + 1
                  IF ( IR_ZERO_TRIED.GT.20 ) THEN
                     IF ( BNDSTACHK ) THEN
                        ITXRAY = -1
                        RETURN
                     END IF
                     WRITE (6,99003) IT,NQN,L,EC,(NRC-1)
                     STOP 'in <CORE>    NZERO not found'
                  END IF
C                                      --------------------
C                                         FIND  IR_MATCH
C                                      --------------------
 60               CONTINUE
                  IR = IR_ZERO + 1
                  DO NN = 1,IR_ZERO
                     IR = IR - 1
                     IF ( (VV(IR)+LLL/RC(IR)**2-EC).LT.0D0 ) THEN
                        IR_MATCH = IR
                        GOTO 65
                     END IF
                  END DO
                  WRITE (6,99004) IT,NQN,L,EC
C
 65               CONTINUE
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),EC,L,MJ,'OUT',VV,BB,
     &                              RC,DRDIC,DOVRC,IR_MATCH,IR_ZERO,GC,
     &                              FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,NRC)
C
                  NODE = 0
                  DO IR = 2,IR_MATCH
                     IF ( GC(S,S,IR)*GC(S,S,IR-1).LT.0D0 ) NODE = NODE + 
     &                    1
                  END DO
C
C
                  IF ( IPRINT.GE.2 ) WRITE (6,99016) IT,NQN,L,KAP(S),
     &                 (2*MUEM05+1),IC,ISH,0,EC,IR_MATCH,RC(IR_MATCH),
     &                 IR_ZERO,RC(IR_ZERO),NODE,
     &                 (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1))
C
                  IF ( NODE.NE.(NQN-L-1) ) THEN
                     IF ( NODE.GT.(NQN-L-1) ) THEN
                        EC = 1.2D0*EC
                     ELSE
                        EC = 0.8D0*EC
                     END IF
                     GOTO 55
                  ELSE IF ( (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1).LE.0D0)
     &                      .OR. 
     &                      (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1).GE.1D0)
     &                      ) THEN
                     EC = 0.9D0*EC
                     GOTO 55
                  END IF
C
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),EC,L,MJ,'INW',VV,BB,
     &                              RC,DRDIC,DOVRC,IR_MATCH,IR_ZERO,GC,
     &                              FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                              CGMD,CGO,NRC)
C
C                                      --------------------
C                                       START VALUES FOR
C                                           PARAMETERS
C                                      --------------------
C
                  VAR(1) = EC
                  VAR(2) = POW(S,S)/PIW(S,S)
C
                  IF ( NSOL.NE.2 .OR. NVAR.EQ.2 ) THEN
                     DO IV = 3,4
                        ERR(IV) = 0.0D0
                        ERRNEW(IV) = 0.0D0
                        VAR(IV) = 0.0D0
                        VARNEW(IV) = 0.0D0
                        DV(IV) = 0.0D0
                     END DO
                  ELSE IF ( ISH.GE.4 ) THEN
                     DO IV = 1,4
                        VAR(IV) = VARTAB(IV,ISH-2)
                     END DO
                  ELSE
C
                     DO J = 1,NSOL
                        NOW(J) = 0.0D0
                     END DO
                     DO IR = 1,IR_MATCH - 1
                        RR = RC(IR)**3
                        DO J = 1,NSOL
                           NOW(J) = NOW(J) + GC(J,J,IR)**2*RR
                        END DO
                     END DO
C
                     DO J = 1,NSOL
                        NIW(J) = 0.0D0
                     END DO
                     DO IR = IR_MATCH,IR_ZERO - 1
                        RR = RC(IR)**3
                        DO J = 1,NSOL
                           NIW(J) = NIW(J) + GC(J,J,IR)**2*RR
                        END DO
                     END DO
C
                     RATT = POW(T,T)/PIW(T,T)
                     VAR(3) = TRYMIX*(NOW(S)+NIW(S)*VAR(2))
     &                        /(NOW(T)+NIW(T)*RATT)
                     VAR(4) = RATT*VAR(3)/VAR(2)
                  END IF
C
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  DO IV = 1,NVAR
                     DV(IV) = VAR(IV)
                  END DO
C
                  ITER = 0
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 70               CONTINUE
                  ITER = ITER + 1
C
                  IF ( SCALEB ) THEN
                     SCLFAC = MIN(1.0D0,0.1D0*ITER)
                     IF ( NSOL.EQ.1 ) SCLFAC = 1.0D0
                     DO IR = 1,JRWS(IM)
                        BB(IR) = BT(IR,IT)*SCLFAC
                     END DO
                  END IF
C                         ----------------------------------
C                         CHECK WHETHER NUMBER OF NODES O.K.
C                         ----------------------------------
                  IF ( ITER.GT.1 ) THEN
                     NODE = 0
                     DO IR = 2,(IR_MATCH-1)
                        IF ( GC(S,S,IR)*GC(S,S,IR-1).LT.0D0 )
     &                       NODE = NODE + 1
                     END DO
                     IF ( IPRINT.GE.3 ) WRITE (6,99016) IT,NQN,L,KAP(S),
     &                    (2*MUEM05+1),IC,ISH,ITER,EC,IR_MATCH,
     &                    RC(IR_MATCH),IR_ZERO,RC(IR_ZERO),NODE,
     &                    (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1))
C
                     IF ( NODE.NE.(NQN-L-1) ) THEN
                        IF ( NODE.GT.(NQN-L-1) ) THEN
                           EC = 1.2D0*EC
                        ELSE
                           EC = 0.8D0*EC
                        END IF
                        ISTART = ISTART + 1
                        IF ( ISTART.LT.20 ) GOTO 55
                     END IF
                  END IF
C
                  DO IV = 2,NVAR
                     DO JV = 1,NVAR
                        VARNEW(JV) = VAR(JV)
                     END DO
                     VARNEW(IV) = VAR(IV) + DV(IV)*DVSTEP
C
                     IF ( ABS(VAR(IV)).GT.1D-16 ) THEN
                        IF ( ABS(DV(IV)/VAR(IV)).LT.TOLVAR ) VARNEW(IV)
     &                       = VAR(IV)
     &                         *(1.0D0+DSIGN(DVSTEP*TOLVAR,DV(IV)))
                     ELSE IF ( FERRO ) THEN
                        WRITE (6,99011) ' VAR(',IV,
     &                                  ') = 0 for (T,N,L,K,M;S,NSOL) ',
     &                                  IT,NQN,L,KAP(S),(2*MUEM05+1),
     &                                  '/2  ',S,NSOL,'  --- suppress B'
                        LOOP = 2
                        SUPPRESSB = .TRUE.
                        GOTO 50
                     ELSE IF ( SUPPRESSB ) THEN
                        WRITE (6,*) 'suppressing B did not help !!'
                        STOP 'in <CORE>'
                     END IF
C
                     CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                     DO IE = 1,NVAR
                        IF ( ABS(ERRNEW(IE)-ERR(IE)).LT.1D-16 ) THEN
                           DEDV(IE,IV) = 1.0D-18
                           IF ( (IE.EQ.IV) .AND. .NOT.FERRO )
     &                          DEDV(IE,IV) = 1.0D0
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
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'OUT',
     &                              VV,BB,RC,DRDIC,DOVRC,IR_MATCH,
     &                              IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                              PIW,QIW,CGD,CGMD,CGO,NRC)
C
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'INW',
     &                              VV,BB,RC,DRDIC,DOVRC,IR_MATCH,
     &                              IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                              PIW,QIW,CGD,CGMD,CGO,NRC)
C
                  CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                  DO IE = 1,NVAR
                     DEDV(IE,1) = (ERRNEW(IE)-ERR(IE))
     &                            /(VARNEW(1)-VAR(1))
                  END DO
C
                  IF ( ABS(DEDV(3,3)).LT.1D-15 .OR. ABS(DEDV(4,4))
     &                 .LT.1D-15 ) THEN
C
                     CALL RMATINV3(2,4,IPIV4,DEDV,RMAT4,DVDE)
C
                     DVDE(3,3) = 0D0
                     DVDE(4,4) = 0D0
C
                  ELSE
C
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
                     IF ( IPRINT.GE.1 ) WRITE (6,*)
     &                     ' warning from <CORE> E=',VAR(1),IT,NQN,L
                     VAR(1) = -0.2D0
                  END IF
C
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),VAR(1),L,MJ,'OUT',VV,
     &                              BB,RC,DRDIC,DOVRC,IR_MATCH,IR_ZERO,
     &                              GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,
     &                              CGD,CGMD,CGO,NRC)
C
                  CALL CORE_DIR_ABM(IT,CTL(IT,ILC),VAR(1),L,MJ,'INW',VV,
     &                              BB,RC,DRDIC,DOVRC,IR_MATCH,IR_ZERO,
     &                              GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,
     &                              CGD,CGMD,CGO,NRC)
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  EC = VAR(1)
C
                  IF ( IPRINT.GE.2 ) WRITE (6,99005) LOOP,SCLFAC,VAR(1),
     &                 (VAR(IV),IV=1,4),(DV(IV),IV=1,4),(ERR(IE),IE=1,4)
C
C----------------------------------  check relative change in parameters
C ----------------------- parameters 3 and 4 = 0 for paramagnetic case !
                  IF ( ITER.LT.ITERMAX ) THEN
                     DO IV = 1,NVAR
                        VARTAB(IV,ISH) = VAR(IV)
                        IF ( (ABS(VAR(IV))+ABS(VAR(IV))).LT.1.0D-30 )  
     &                       THEN
                           IF ( FERRO ) WRITE (6,'(A,I3,A)') ' VAR ',IV,
     &                          ' = 0 ??????!!!!!'
                        ELSE IF ( ABS(DV(IV)/VAR(IV)).GT.TOLVAR ) THEN
                           GOTO 70
                        END IF
                     END DO
                  ELSE
                     IF ( BNDSTACHK ) THEN
                        ITXRAY = -1
                        RETURN
                     END IF
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
                  DO IR = IR_MATCH,IR_ZERO
                     DO J = 1,NSOL
                        DO I = 1,NSOL
                           GC(I,J,IR) = GC(I,J,IR)*VAR(2)
                           FC(I,J,IR) = FC(I,J,IR)*VAR(2)
                        END DO
                     END DO
                  END DO
C
                  IF ( NSOL.EQ.2 ) THEN
C                                   OUTWARD - SOLUTION
                     DO IR = 1,(IR_MATCH-1)
                        DO I = 1,NSOL
                           GC(I,T,IR) = GC(I,T,IR)*VAR(3)
                           FC(I,T,IR) = FC(I,T,IR)*VAR(3)
                        END DO
                     END DO
C                                    INWARD - SOLUTION
                     DO IR = IR_MATCH,IR_ZERO
                        DO I = 1,NSOL
                           GC(I,T,IR) = GC(I,T,IR)*VAR(4)
                           FC(I,T,IR) = FC(I,T,IR)*VAR(4)
                        END DO
                     END DO
                  END IF
C
C                                    SUM FOR EACH KAPPA
                  DO IR = 1,IR_ZERO
                     DO K = 1,NSOL
                        GCK(K,S,IR) = 0.0D0
                        FCK(K,S,IR) = 0.0D0
                        DO J = 1,NSOL
                           GCK(K,S,IR) = GCK(K,S,IR) + GC(K,J,IR)
                           FCK(K,S,IR) = FCK(K,S,IR) + FC(K,J,IR)
                        END DO
                     END DO
                  END DO
C
C                       -----------------------------------
C                       CALCULATE  NORM  AND NORMALIZE TO 1
C                       -----------------------------------
C
                  DO IR = (IR_ZERO+1),NRC
                     DO K = 1,NSOL
                        GCK(K,S,IR) = 0.0D0
                        FCK(K,S,IR) = 0.0D0
                     END DO
                  END DO
C
                  AUX = 0D0
                  DO IR = 1,IRTOP_SPHERE
                     DO K = 1,NSOL
                        AUX = AUX + R2DRDI_W_RADINT(IR,IM)
     &                        *(GCK(K,S,IR)*GCK(K,S,IR)+FCK(K,S,IR)
     &                        *FCK(K,S,IR))
                     END DO
                  END DO
C
                  IF ( FULLPOT ) THEN
                     DO IRSF = 1,NRSFTOT(IM)
                        IR = IRTOP_SPHERE + IRSF
                        DO K = 1,NSOL
                           AUX = AUX + R2DRDI_W_RADINT(IR,IM)
     &                           *(GCK(K,S,IR)*GCK(K,S,IR)+FCK(K,S,IR)
     &                           *FCK(K,S,IR))*FLMSF(IRSF,ISF_00,IM)
     &                           /SQRT_4PI
                        END DO
                     END DO
                  END IF
C
                  NORM = 1.0D0/SQRT(AUX)
C
                  DO IR = 1,IR_ZERO
                     DO K = 1,NSOL
                        GCK(K,S,IR) = GCK(K,S,IR)*NORM
                        FCK(K,S,IR) = FCK(K,S,IR)*NORM
                     END DO
                  END DO
C
C                       -----------------------------------
C                            CALCULATE  SPIN CHARACTER
C                       -----------------------------------
C
                  SZ = 0.0D0
                  DO IR = 1,IRTOP
                     W = R2DRDI_W_RADINT(IR,IM)
                     DO K = 1,NSOL
                        SZ = SZ + W*(GCK(K,S,IR)**2*CGD(K)+FCK(K,S,IR)
     &                       **2*CGMD(K))
                     END DO
                  END DO
C
                  IF ( NSOL.GT.1 ) THEN
                     DO IR = 1,IRTOP
                        W = R2DRDI_W_RADINT(IR,IM)
                        SZ = SZ + W*GCK(1,S,IR)*GCK(2,S,IR)*CGO*2
                     END DO
                  END IF
C
C-----------------------------------------------------------------------
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
                     DO IR = 1,IRTOP
                        W = R2DRDI_W_RADINT(IR,IM)
                        BEX = BEX + W*(CGD(1)*GCK(1,S,IR)**2+CGMD(1)
     &                        *FCK(1,S,IR)**2)*BB(IR)
                     END DO
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
                        IF ( ITER.LT.ITERMAX ) GOTO 70
                     END IF
C
C--- check 3 - check whether the two states with same MU but opposite SZ
C
                     IF ( ABS(SZ+SZ1).GT.0.1D0 .AND. ABS(BEX)
     &                    .GT.1D-6 .AND. ITER.LT.ITERMAX ) GOTO 70
                  END IF
C-----------------------------------------------------------------------
C--- END ----------- Consistency check on 2nd solution -----------------
C-----------------------------------------------------------------------
C
C                       -----------------------------------
C                       CALCULATE  CHARGE AND SPIN DENSITY
C                       -----------------------------------
C
                  DO IR = 1,JRWS(IM)
                     DO K = 1,NSOL
                        RHOCHRC(IR,IT) = RHOCHRC(IR,IT)
     &                     + (GCK(K,S,IR)**2+FCK(K,S,IR)**2)
                        RHOSPNC(IR,IT) = RHOSPNC(IR,IT)
     &                     + (GCK(K,S,IR)**2*CGD(K)-FCK(K,S,IR)
     &                     **2*CGMD(K))
                     END DO
                  END DO
C
                  IF ( NSOL.GT.1 ) THEN
                     DO IR = 1,JRWS(IM)
                        RHOSPNC(IR,IT) = RHOSPNC(IR,IT) + GCK(1,S,IR)
     &                     *GCK(2,S,IR)*CGO*2
                     END DO
                  END IF
C
C                         ------------------------------
C                         CALCULATE   HYPERFINE - FIELD
C                         ------------------------------
C
                  IF ( Z(IT).NE.0 ) THEN
C
                     CALL CORE_HFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,
     &                             RC,DRDIC,0.0D0,IR_ZERO,IRTOP,NRC)
C
                     IF ( FINITE_NUCLEUS ) THEN
                        CALL CORE_HFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF1,GCK,
     &                                FCK,RC,DRDIC,RNUC,JLIM,IRTOP,NRC)
                        CALL CORE_HFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF2,GCK,
     &                                FCK,RC,DROVRN,RNUC,JLIM,IRTOP,NRC)
C
                        DO I = 1,NSOL
                           DO J = 1,NSOL
                              BHF(I,J) = BHF(I,J) - BHF1(I,J)
     &                           + BHF2(I,J)
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
C
                  ELSE
C
                     BHF(:,:) = 0D0
                     BCORS(IT) = 0D0
C
                  END IF
C
                  ECORTAB(IC,IT) = EC
                  ECOR_LT(IL,IT) = ECOR_LT(IL,IT) + EC
C
C     ------------------
C     SPLIT HFF-FIELD
C     ------------------
                  IF ( ISMQHFI.EQ.10 .AND. Z(IT).NE.0 ) THEN
                     CALL CORE_HFF_SPLIT(IM,RNUC,IRTOP,KAP1,KAP2,NSOL,
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
                  ICSTDUM = ICSTDUM + 1
                  SZCORDUM(ICSTDUM) = SZ
                  ECORDUM(ICSTDUM) = EC
C
C---------------------------------------------------- l-shell UNCOMPLETE
                  IF ( ISH.GE.NSH ) THEN
C----------------------------------------------------- l-shell COMPLETED
                     IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) THEN
                        WRITE (6,99012) IT,NQN,TXT_L(L),
     &                                  TXTK(IABS(KAP(S))),(2*MUEM05+1),
     &                                  KAP(S),ITER,EC,BSOL*.001D0,SZ,
     &                                  BSH*.001D0
                        IF ( ISMQHFI.EQ.10 ) THEN
                           DO K = 1,NMEHFMAX
                              WRITE (6,99013) TXTB(K),SPLIT(K)*.001D0,
     &                               SPLIT1(K)*.001D0
                           END DO
                           WRITE (6,99014) 'total error in %',
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
                        IF ( FERRO .AND. (IFLAG.EQ.1) ) THEN
                           WRITE (6,99007)
                           SCALEB = .TRUE.
                           IF ( LOOP.EQ.1 ) THEN
                              LOOP = 2
                              WRITE (6,99008) IT
                              GOTO 50
                           END IF
                        END IF
C
                     END IF
                  ELSE IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) THEN
                     WRITE (6,99012) IT,NQN,TXT_L(L),TXTK(IABS(KAP(S))),
     &                               (2*MUEM05+1),KAP(S),ITER,EC,
     &                               BSOL*.001D0,SZ
                     IF ( ISMQHFI.EQ.10 ) THEN
                        DO K = 1,NMEHFMAX
                           WRITE (6,99013) TXTB(K),SPLIT(K)*.001D0
                        END DO
                        WRITE (6,99014) 'total error in %',
     &                                  100.0D0*(1.0D0-SPLIT(4)/SPLIT(5)
     &                                  )
                     END IF
                  END IF
C-----------------------------------------------------------------------
C
                  IF ( IPRINT.GE.1 ) WRITE (6,99015)
     &                 ((BHF(I,J)*.001D0,I=1,NSOL),J=1,NSOL)
C
C???????????????????????????????????????????????????????????????????????
                  IF ( CHECK_AGAINST_OLD_ITERATION_SCHEME ) THEN
C
C                          --------------------------------
C                            IF THE SWITCH CHECK IS SET:
C                            RECALCULATE THE EIGENVALUE
C                          USING THE CONVENTIONAL ALGORITHM
C                          --------------------------------
C
                     ECC = 0.95D0*EC
C
 72                  CONTINUE
                     CALL CORE_DIR_ABM(IT,CTL(IT,ILC),ECC,L,MJ,'OUT',VV,
     &                                 BB,RC,DRDIC,DOVRC,IR_MATCH,
     &                                 IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,
     &                                 QOW,PIW,QIW,CGD,CGMD,CGO,NRC)
C
                     CALL CORE_DIR_ABM(IT,CTL(IT,ILC),ECC,L,MJ,'INW',VV,
     &                                 BB,RC,DRDIC,DOVRC,IR_MATCH,
     &                                 IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,
     &                                 QOW,PIW,QIW,CGD,CGMD,CGO,NRC)
C
                     NORM = POW(S,S)/PIW(S,S)
                     DO IR = IR_MATCH,IR_ZERO
                        GC(S,S,IR) = GC(S,S,IR)*NORM
                        FC(S,S,IR) = FC(S,S,IR)*NORM
                     END DO
C
                     NORM = 0.0D0
                     DO IR = 1,IRTOP_SPHERE
                        WR = R2DRDI_W_RADINT(IR,IM)
                        NORM = NORM + WR*(GC(S,S,IR)**2+FC(S,S,IR)**2)
                     END DO
C
                     IF ( FULLPOT ) THEN
                        DO IRSF = 1,NRSFTOT(IM)
                           IR = IRTOP_SPHERE + IRSF
                           WR = R2DRDI_W_RADINT(IR,IM)
     &                          *FLMSF(IRSF,ISF_00,IM)/SQRT_4PI
                           NORM = NORM + 
     &                            WR*(GC(S,S,IR)**2+FC(S,S,IR)**2)
                        END DO
                     END IF
C
                     DEC = POW(S,S)
     &                     *(QOW(S,S)-RC(IR_MATCH)*CTL(IT,ILC)*FC(S,S,
     &                     IR_MATCH))/NORM
                     ECC = ECC + DEC
C
                     IF ( ABS(DEC/ECC).GT.TOLVAR ) GOTO 72
C
                     WRITE (6,'(7X,''CHECK-E:'',10X,F12.5,/)') ECC
                  END IF
C???????????????????????????????????????????????????????????????????????
C
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C                STORE CORE WAVE FUNCTIONS IF REQUESTED
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
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
                     IZERO(ICST) = MIN(IR_ZERO,JRWS(IM))
                     SZCOR(ICST) = SZ
                     ECOR(ICST) = EC
C
                     DO IR = 1,IZERO(ICST)
                        GCOR(IR,1,ICST) = GCK(S,S,IR)
                        FCOR(IR,1,ICST) = FCK(S,S,IR)
                     END DO
                     IF ( NSOL.EQ.2 ) THEN
                        DO IR = 1,IZERO(ICST)
                           GCOR(IR,2,ICST) = GCK(T,S,IR)
                           FCOR(IR,2,ICST) = FCK(T,S,IR)
                        END DO
                        IKMCOR(ICST,2) = IKAPMUE(KAP(T),MUEM05)
                     END IF
                  END IF
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
               END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GE.0 ) WRITE (6,99009) BCOR(IT)
     &        *.001D0,BCORS(IT)*.001D0
         IF ( ISMQHFI.EQ.10 ) THEN
            DO IMEHF = 1,NMEHFMAX
               WRITE (6,99010) SPLIT2(IMEHF,IT)*.001D0,SPLIT3(IMEHF,IT)
     &                         *.001D0
            END DO
            WRITE (6,99014) 'total error',
     &                      100.0D0*(1.0D0-SPLIT2(4,IT)/SPLIT2(5,IT))
            WRITE (6,99014) 'total error',
     &                      100.0D0*(1.0D0-SPLIT3(4,IT)/SPLIT3(5,IT))
         END IF
C
         IF ( (ITXRAY.EQ.0 .OR. IT.EQ.ITXRAY) .AND. TASK(1:6)
     &        .EQ.'POSANI' .AND. TASK(1:7).NE.'COMPTON' ) THEN
            IF ( IPRINT.GT.-2 ) THEN
C
               INTEG_CHR = 0D0
               INTEG_SPN = 0D0
               DO IR = 1,IRTOP_SPHERE
                  INTEG_CHR = INTEG_CHR + RHOCHRC(IR,IT)
     &                        *R2DRDI_W_RADINT(IR,IM)
                  INTEG_SPN = INTEG_SPN + RHOSPNC(IR,IT)
     &                        *R2DRDI_W_RADINT(IR,IM)
               END DO
C
               IF ( FULLPOT ) THEN
                  DO IRSF = 1,NRSFTOT(IM)
                     IR = IRTOP_SPHERE + IRSF
                     WR = R2DRDI_W_RADINT(IR,IM)*FLMSF(IRSF,ISF_00,IM)
     &                    /SQRT_4PI
C
                     INTEG_CHR = INTEG_CHR + RHOCHRC(IR,IT)*WR
                     INTEG_SPN = INTEG_SPN + RHOSPNC(IR,IT)*WR
                  END DO
               END IF
C
               WRITE (6,99017) 'charge',IT,INTEG_CHR
               WRITE (6,99017) ' spin ',IT,INTEG_SPN
C
            END IF
         END IF
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C                STORE CORE WAVE FUNCTIONS IF REQUESTED
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
         IF ( ITXRAY.NE.0 .AND. IT.EQ.ITXRAY )
     &        CALL CORE_WRITE_WAVFUN(GCOR,FCOR,MM05COR,IKMCOR,ITXRAY,
     &        NCSTMAX)
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
      END DO LOOP_IT
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (/,10X,'potential is not exchange split ')
99002 FORMAT (/,'  ATOM   IT      : ',I5,/,'  ATOMIC NUMBER  : ',I5,//,
     &        '  IT',10X,'MUE  KAP ITER      ENERGY',
     &        '        B (k-Gauss)      <S_z>  ')
99003 FORMAT ('  IT=',I4,' NQN=',I2,' L=',I2,' EC=',F10.3,
     &        '  NZERO set to  (NRC-1) =',I4)
99004 FORMAT (//,'  STOP IN <<CORE>>',/,'  IT=',I4,' NQN=',I2,' L=',I2,
     &        /,'  no matching-radius found for  EC=',F10.3)
99005 FORMAT (' LOOP    =  ',I3,' BSCL=',F10.5,/,' E=',F14.7,' VAR  ',
     &        4E11.4,/,17X,' CORR ',4E11.4,/,17X,' ERR  ',4E11.4)
99006 FORMAT (' iteration not converged after',I3,' steps !',/,
     &        ' parameters:',4E18.10,/,' last corr.:',4E18.10,/,
     &        ' last error:',4E18.10)
99007 FORMAT (' >>> check data E(KAP,MJ) should be monotonous ',
     &        ' and  E(+L,MJ) < E(-L-1,MJ) ',//)
99008 FORMAT (' all core states for atom type IT=',I4,
     &        ' will be recalculated ',/,
     &   ' with the spin dependent exchange field gradually switched on'
     &   )
99009 FORMAT (2X,57('-'),/,42X,F17.3,/,38X,'(S) ',F17.3,/,2X,57('*'),//)
99010 FORMAT (2X,57('-'),/,37X,F17.3,/,33X,'(S) ',F17.3,/,2X,57('*'),//)
99011 FORMAT (A,I1,A,5I4,A,2I2,A)
99012 FORMAT (2I4,A1,A3,I3,'/2',2I4,2X,F15.8,F17.3,F10.3,:,/,42X,F17.3,
     &        /)
99013 FORMAT (A,:,32X,F17.3,:,F17.3,/)
99014 FORMAT (A,:,37X,F6.3)
99015 FORMAT (37X,F17.3)
99016 FORMAT (/,' IT=',I4,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,
     &        '/2    IC=',I3,'  ISH=',I2,/,' E(',I2,')   =',F15.5,/,
     &        ' IR_MATCH=',I5,'    R=',F10.5,/,' IR_ZERO =',I5,'    R=',
     &        F10.5,/,' NODES   =',I5,'  RAT=',E11.4)
99017 FORMAT (10X,'integrated core ',A,' density for atom type ',I4,':',
     &        F12.8)
      END
C*==core_write_wavfun.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE_WRITE_WAVFUN(GCOR,FCOR,MM05COR,IKMCOR,ITXRAY,
     &                             NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  write core wave functions to file   IFILCORWF                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,FULLPOT
      USE MOD_TYPES,ONLY:LCXRAY,IMT,NCPLWFMAX,NT
      USE MOD_ANGMOM,ONLY:NKM
      USE MOD_FILES,ONLY:IFILCORWF
      IMPLICIT NONE
C*--CORE_WRITE_WAVFUN1261
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CORE_WRITE_WAVFUN')
C
C Dummy arguments
C
      INTEGER ITXRAY,NCSTMAX
      REAL*8 FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),MM05COR(NCSTMAX)
C
C Local variables
C
      INTEGER I,ICST,IKM,IKMCPLWFCOR(NCPLWFMAX),IM,IRTOP,IT,J,K,NCST,
     &        NSOL,NSOLDUM
C
C*** End of declarations rewritten by SPAG
C
      IT = ITXRAY
      IM = IMT(IT)
C
      NCST = 2*(2*LCXRAY(IT)+1)
C
      IF ( NCST.GT.NCSTMAX ) CALL STOP_MESSAGE(ROUTINE,'NCST > NCSTMAX')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
C-----------------------------------------------------------------------
C fill file for core wave function with dummy data for ALL l=1,..(NL-1)
C-----------------------------------------------------------------------
C
      NSOLDUM = 1
C
      DO IKM = 1,NKM
C
         WRITE (IFILCORWF,REC=IKM+(IT-1)*NKM) IT,'COR',IKM,IRTOP,
     &          NSOLDUM,(IKM,K=1,NSOLDUM),((C0,I=1,IRTOP),K=1,NSOLDUM),
     &          ((C0,I=1,IRTOP),K=1,NSOLDUM)
C
         WRITE (IFILCORWF,REC=IKM+(IT-1+NT)*NKM) IT,'IRR',IKM,IRTOP,
     &          ((C0,I=1,IRTOP),K=1,NSOLDUM),
     &          ((C0,I=1,IRTOP),K=1,NSOLDUM)
C
      END DO
C
C=======================================================================
C              write core wave function for present core shell
C                  including dummy irregular solution
C=======================================================================
C
      DO ICST = 1,NCST
C
         IF ( ABS(MM05COR(ICST)+0.5D0).GT.LCXRAY(IT) ) THEN
            NSOL = 1
         ELSE
            NSOL = 2
         END IF
         DO J = 1,NSOL
            IKMCPLWFCOR(J) = IKMCOR(ICST,J)
         END DO
C
         WRITE (IFILCORWF,REC=IKMCOR(ICST,1)+(IT-1)*NKM) IT,'COR',
     &          IKMCOR(ICST,1),IRTOP,NSOL,(IKMCPLWFCOR(K),K=1,NSOL),
     &          ((DCMPLX(GCOR(I,K,ICST),0.0D0),I=1,IRTOP),K=1,NSOL),
     &          ((DCMPLX(FCOR(I,K,ICST),0.0D0),I=1,IRTOP),K=1,NSOL)
C
         WRITE (IFILCORWF,REC=IKMCOR(ICST,1)+(IT-1+NT)*NKM) IT,'IRR',
     &          IRTOP,((C0,I=1,IRTOP),K=1,NSOL),
     &          ((C0,I=1,IRTOP),K=1,NSOL)
C
      END DO
C
      END
