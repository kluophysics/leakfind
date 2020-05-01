C*==core_dir_abm.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE_DIR_ABM(IT,C,E,L,MJ,WAY,VV,BB,RC,DRDIC,DOVRC,
     &                        IR_MATCH,IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,
     &                        QOW,PIW,QIW,CGD,CGMD,CGO,NRC)
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
C
      USE MOD_TYPES,ONLY:IMT,Z
      USE MOD_RMESH,ONLY:FULLPOT,FINITE_NUCLEUS,JRCUT,NPAN,JRWS
      IMPLICIT NONE
C*--CORE_DIR_ABM19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,INVMAX
      PARAMETER (MPSMAX=20,NPEMAX=20,INVMAX=3)
      REAL*8 TOL
      PARAMETER (TOL=1.0D-8) !REC changed from 1d-9
      INTEGER ITMAX,NABM
      PARAMETER (ITMAX=100,NABM=4) !REC changed from 50
C
C Dummy arguments
C
      REAL*8 C,CGO,E,MJ
      INTEGER IR_MATCH,IR_ZERO,IT,L,NRC
      CHARACTER*3 WAY
      REAL*8 BB(NRC),CGD(2),CGMD(2),DOVRC(NRC),DP(2,2,NRC),DQ(2,2,NRC),
     &       DRDIC(NRC),FC(2,2,NRC),GC(2,2,NRC),PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),RC(NRC),VV(NRC),WP(2,2,NRC),WQ(2,2,NRC)
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
      INTEGER I,IM,IPAN,IPANMATCH,IPANNZERO,IR,IRK,IV,J,JCORR,K,KAP1,
     &        KAP2,M,MPS,NDIV,NEND,NHLP,NMESHLC,NN,NSOL,NSTART
      LOGICAL RNON0
      REAL*8 W5,W6,W7,X14,XH
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      DIFFA = 0.0D0
      DIFFB = 0.0D0
C
      H24 = 1.0D0/24.0D0
      DVC = C
      CSQR = DVC*DVC
C
      IM = IMT(IT)
C
C refinement factor for Runge-Kutta-starter
C
      NDIV = 60
C
C determine nearest cut point lower than nmatch,nzero
C
      IF ( FULLPOT ) THEN
C
         IPANMATCH = NPAN(IM)
         DO WHILE ( JRCUT(IPANMATCH,IM).GE.IR_MATCH )
            IPANMATCH = IPANMATCH - 1
         END DO
C
         IPANNZERO = NPAN(IM)
         DO WHILE ( JRCUT(IPANNZERO,IM).GE.IR_ZERO )
            IPANNZERO = IPANNZERO - 1
         END DO
C
      ELSE
C
         IPANMATCH = 1
         IPANNZERO = 1
C
         NPAN(IM) = 1
         JRCUT(0,IM) = 0
         JRCUT(1,IM) = JRWS(IM)
C
      END IF
C
C     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
C
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         TZ = DBLE(NINT(-VV(1)*RC(1)))
         VC(0) = VV(1) - (-TZ)/RC(1)
      ELSE
         TZ = 2.0D0*DBLE(Z(IT))
         VC(0) = VV(1)
      END IF
C
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
C
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
C
         IF ( .NOT.FINITE_NUCLEUS ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/DVC)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
C
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
         DO IR = (IR_ZERO-3),IR_ZERO
C
            RR = RC(IR)
C
            DO J = 1,NSOL
               I = 3 - J
               WP(J,J,IR) = DEXP(-DMUE*RR)
               DP(J,J,IR) = -DMUE*DRDIC(IR)*WP(J,J,IR)
               WQ(J,J,IR) = BOVA*WP(J,J,IR)
               DQ(J,J,IR) = BOVA*DP(J,J,IR)
C
               WP(I,J,IR) = 0.0D0
               WQ(I,J,IR) = 0.0D0
               DP(I,J,IR) = 0.0D0
               DQ(I,J,IR) = 0.0D0
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
               IR = NMESHLC
C
               IF ( IPAN.LT.NPAN(IM) ) THEN
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        WP(I,J,IR) = WP(I,J,IR+1)
                        WQ(I,J,IR) = WQ(I,J,IR+1)
                     END DO
                  END DO
               END IF
C
               EMVQQ = (E-VV(IR)+CSQR)*DRDIC(IR)/CSQR
               EMVPP = -(E-VV(IR))*DRDIC(IR)
               BQQ = BB(IR)*DRDIC(IR)/CSQR
               BPP = BB(IR)*DRDIC(IR)
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     DP(I,J,IR) = -KAP(I)*WP(I,J,IR)*DOVRC(IR)
     &                            + (EMVQQ+BQQ*CGMD(I))*WQ(I,J,IR)
C
                     DQ(I,J,IR) = KAP(I)*WQ(I,J,IR)*DOVRC(IR)
     &                            + (EMVPP+BPP*CGD(I))*WP(I,J,IR)
     &                            + BPP*CGO*WP(3-I,J,IR)
C
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     P1(I,J) = WP(I,J,IR)
                     Q1(I,J) = WQ(I,J,IR)
                     MP1(I,J) = DP(I,J,IR)
                     MQ1(I,J) = DQ(I,J,IR)
                  END DO
               END DO
C
               X14 = DBLE(IR)
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
                     IR = NMESHLC - IRK/NDIV
                     IF ( ABS(X14-DBLE(IR)).GT.1.0D-5 ) THEN
                        WRITE (*,*) ' <dirac> runge-kutta: ',IRK,NDIV,
     &                              IR,X14
                        STOP
                     END IF
                     DO J = 1,NSOL
                        DO I = 1,NSOL
                           WP(I,J,IR) = P1(I,J)
                           WQ(I,J,IR) = Q1(I,J)
                           DP(I,J,IR) = MP1(I,J)
                           DQ(I,J,IR) = MQ1(I,J)
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
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
C
            IF ( .NOT.(FULLPOT) ) THEN
               NSTART = IR_ZERO
               NEND = IR_MATCH
            ELSE IF ( IPAN.EQ.NPAN(IM)+1 ) THEN
               NSTART = IR_ZERO
               NEND = MAX(IR_MATCH,JRCUT(NPAN(IM),IM))
            ELSE
               NSTART = MIN(IR_ZERO,JRCUT(IPAN,IM))
               NEND = MAX(IR_MATCH,JRCUT(IPAN-1,IM)+1)
            END IF
C
            DO NN = 1,(NSTART-3-NEND)
C
               IR = NSTART - 3 - NN
C
C    EVALUATE PREDICTOR
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     PNEW(I,J) = WP(I,J,IR+1)
     &                           - H24*(55.0D0*DP(I,J,IR+1)-
     &                           59.0D0*DP(I,J,IR+2)+37.0D0*DP(I,J,IR+3)
     &                           -9.0D0*DP(I,J,IR+4))
                     QNEW(I,J) = WQ(I,J,IR+1)
     &                           - H24*(55.0D0*DQ(I,J,IR+1)-
     &                           59.0D0*DQ(I,J,IR+2)+37.0D0*DQ(I,J,IR+3)
     &                           -9.0D0*DQ(I,J,IR+4))
                  END DO
               END DO
C
               EMVQQ = (E-VV(IR)+CSQR)*DRDIC(IR)/CSQR
               EMVPP = -(E-VV(IR))*DRDIC(IR)
               BQQ = BB(IR)*DRDIC(IR)/CSQR
               BPP = BB(IR)*DRDIC(IR)
C
C    EVALUATE CORRECTOR
C
               DO JCORR = 1,ITMAX
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        POLD(I,J) = PNEW(I,J)
                        QOLD(I,J) = QNEW(I,J)
                        DP(I,J,IR) = -KAP(I)*PNEW(I,J)*DOVRC(IR)
     &                               + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                        DQ(I,J,IR) = KAP(I)*QNEW(I,J)*DOVRC(IR)
     &                               + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                               + BPP*CGO*PNEW(3-I,J)
C
                        PNEW(I,J) = WP(I,J,IR+1)
     &                              - H24*(9.0D0*DP(I,J,IR)+
     &                              19.0D0*DP(I,J,IR+1)
     &                              -5.0D0*DP(I,J,IR+2)+DP(I,J,IR+3))
                        QNEW(I,J) = WQ(I,J,IR+1)
     &                              - H24*(9.0D0*DQ(I,J,IR)+
     &                              19.0D0*DQ(I,J,IR+1)
     &                              -5.0D0*DQ(I,J,IR+2)+DQ(I,J,IR+3))
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
               WRITE (6,99001) KAP1,IR,RC(IR),DIFFA,DIFFB,IT,L,INT(2*MJ)
     &                         ,' IN'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
C
 20            CONTINUE
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     WP(I,J,IR) = PNEW(I,J)
                     WQ(I,J,IR) = QNEW(I,J)
                     DP(I,J,IR) = -KAP(I)*PNEW(I,J)*DOVRC(IR)
     &                            + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                     DQ(I,J,IR) = KAP(I)*QNEW(I,J)*DOVRC(IR)
     &                            + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                            + BPP*CGO*PNEW(3-I,J)
                  END DO
               END DO
C
            END DO
C
         END DO
C              ! ipan - loop
C
C =============================================================== N ====
C
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
         DO IR = IR_MATCH,IR_ZERO
            DO J = 1,NSOL
               DO I = 1,NSOL
                  GC(I,J,IR) = WP(I,J,IR)/RC(IR)
                  FC(I,J,IR) = WQ(I,J,IR)/(RC(IR)*DVC)
               END DO
            END DO
         END DO
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PIW(I,J) = WP(I,J,IR_MATCH)
               QIW(I,J) = WQ(I,J,IR_MATCH)
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
C
      IF ( .NOT.FINITE_NUCLEUS ) THEN
C
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
      DO IR = 1,4
         RR = RC(IR)
C
         DO J = 1,NSOL
            RPWGPM = RR**GAM(J)
C
            DO I = 1,NSOL
               WP(I,J,IR) = PC(I,J,0)*RPWGPM
               WQ(I,J,IR) = QC(I,J,0)*RPWGPM
               DP(I,J,IR) = PC(I,J,0)*RPWGPM*GAM(J)*DOVRC(IR)
               DQ(I,J,IR) = QC(I,J,0)*RPWGPM*GAM(J)*DOVRC(IR)
            END DO
C
            DO M = 1,MPS
               RPWGPM = RPWGPM*RR
               GPM = GAM(J) + M
C
               DO I = 1,NSOL
                  WP(I,J,IR) = WP(I,J,IR) + PC(I,J,M)*RPWGPM
                  WQ(I,J,IR) = WQ(I,J,IR) + QC(I,J,M)*RPWGPM
                  DP(I,J,IR) = DP(I,J,IR) + PC(I,J,M)
     &                         *RPWGPM*GPM*DOVRC(IR)
                  DQ(I,J,IR) = DQ(I,J,IR) + QC(I,J,M)
     &                         *RPWGPM*GPM*DOVRC(IR)
               END DO
C
            END DO
         END DO
      END DO
C
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
               IR = JRCUT(IPAN,IM) + 1
            ELSE
               IR = JRCUT(IPAN,IM)
            END IF
C
            NSTART = IR
C
            IF ( IPAN.LT.NPAN(IM) ) THEN
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     WP(I,J,IR) = WP(I,J,IR-1)
                     WQ(I,J,IR) = WQ(I,J,IR-1)
                  END DO
               END DO
            END IF
C
            EMVQQ = (E-VV(IR)+CSQR)*DRDIC(IR)/CSQR
            EMVPP = -(E-VV(IR))*DRDIC(IR)
            BQQ = BB(IR)*DRDIC(IR)/CSQR
            BPP = BB(IR)*DRDIC(IR)
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  DP(I,J,IR) = -KAP(I)*WP(I,J,IR)*DOVRC(IR)
     &                         + (EMVQQ+BQQ*CGMD(I))*WQ(I,J,IR)
C
                  DQ(I,J,IR) = KAP(I)*WQ(I,J,IR)*DOVRC(IR)
     &                         + (EMVPP+BPP*CGD(I))*WP(I,J,IR)
     &                         + BPP*CGO*WP(3-I,J,IR)
C
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P1(I,J) = WP(I,J,IR)
                  Q1(I,J) = WQ(I,J,IR)
                  MP1(I,J) = DP(I,J,IR)
                  MQ1(I,J) = DQ(I,J,IR)
               END DO
            END DO
C
            X14 = DBLE(IR)
            NHLP = NABM + 4
            HLP1 = DBLE(IR) - 1
            DO I = 1,NHLP
               HLP(I) = DBLE(I)
               VHLP(I) = VV(IR+I-1)
               BHLP(I) = BB(IR+I-1)
               DHLP(I) = DRDIC(IR+I-1)
               RHLP(I) = RC(IR+I-1)
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
                  IR = NSTART + IRK/NDIV
                  IF ( ABS(X14-DBLE(IR)).GT.1.0D-5 ) THEN
                     WRITE (*,*) ' <dirac> runge-kutta: ',IRK,NDIV,IR,
     &                           X14
                     STOP
                  END IF
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        WP(I,J,IR) = P1(I,J)
                        WQ(I,J,IR) = Q1(I,J)
                        DP(I,J,IR) = MP1(I,J)
                        DQ(I,J,IR) = MQ1(I,J)
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
         IF ( FULLPOT ) THEN
            NEND = MIN(IR_MATCH,JRCUT(IPAN+1,IM))
         ELSE
            NEND = IR_MATCH
         END IF
C
         DO IR = NSTART,NEND
C
C    EVALUATE PREDICTOR
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PNEW(I,J) = WP(I,J,IR-1)
     &                        + H24*(55.0D0*DP(I,J,IR-1)-59.0D0*DP(I,J,
     &                        IR-2)+37.0D0*DP(I,J,IR-3)
     &                        -9.0D0*DP(I,J,IR-4))
                  QNEW(I,J) = WQ(I,J,IR-1)
     &                        + H24*(55.0D0*DQ(I,J,IR-1)-59.0D0*DQ(I,J,
     &                        IR-2)+37.0D0*DQ(I,J,IR-3)
     &                        -9.0D0*DQ(I,J,IR-4))
               END DO
            END DO
C
            EMVQQ = (E-VV(IR)+CSQR)*DRDIC(IR)/CSQR
            EMVPP = -(E-VV(IR))*DRDIC(IR)
            BQQ = BB(IR)*DRDIC(IR)/CSQR
            BPP = BB(IR)*DRDIC(IR)
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
                     DP(I,J,IR) = -KAP(I)*PNEW(I,J)*DOVRC(IR)
     &                            + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                     DQ(I,J,IR) = KAP(I)*QNEW(I,J)*DOVRC(IR)
     &                            + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                            + BPP*CGO*PNEW(3-I,J)
C
                     PNEW(I,J) = WP(I,J,IR-1)
     &                           + H24*(9.0D0*DP(I,J,IR)+19.0D0*DP(I,J,
     &                           IR-1)-5.0D0*DP(I,J,IR-2)+DP(I,J,IR-3))
                     QNEW(I,J) = WQ(I,J,IR-1)
     &                           + H24*(9.0D0*DQ(I,J,IR)+19.0D0*DQ(I,J,
     &                           IR-1)-5.0D0*DQ(I,J,IR-2)+DQ(I,J,IR-3))
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
            WRITE (6,99001) KAP1,IR,RC(IR),DIFFA,DIFFB,IT,L,INT(2*MJ),
     &                      'OUT'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
 60         CONTINUE
            DO J = 1,NSOL
               DO I = 1,NSOL
                  WP(I,J,IR) = PNEW(I,J)
                  WQ(I,J,IR) = QNEW(I,J)
                  DP(I,J,IR) = -KAP(I)*PNEW(I,J)*DOVRC(IR)
     &                         + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,IR) = KAP(I)*QNEW(I,J)*DOVRC(IR)
     &                         + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                         + BPP*CGO*PNEW(3-I,J)
               END DO
C
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
      DO IR = 1,IR_MATCH
         DO J = 1,NSOL
            DO I = 1,NSOL
               GC(I,J,IR) = WP(I,J,IR)/RC(IR)
               FC(I,J,IR) = WQ(I,J,IR)/(RC(IR)*DVC)
            END DO
         END DO
      END DO
C
      DO J = 1,NSOL
         DO I = 1,NSOL
            POW(I,J) = WP(I,J,IR_MATCH)
            QOW(I,J) = WQ(I,J,IR_MATCH)
         END DO
      END DO
C
      RETURN
C
99001 FORMAT (' P/C NOT CONV. IN <DIRAC> ',2I4,2X,F10.7,2X,2E12.4,3I2,
     &        '/2 ',A3)
99999 CONTINUE
      END
