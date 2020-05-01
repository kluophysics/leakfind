C*==dirbsbi.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSBI(GETIRRSOL,C,E,L,MJ,KAP1,KAP2,PIS,CG1,CG2,CG4,
     &                   CG5,CG8,V,B,Z,R,DRDI,DOVR,IRTOP,IRCUT,NPAN,PR,
     &                   QR,PI,QI,DP,DQ,AP,AQ,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and the inward integration is started analytically             *
C   *   the integration itself is done by the BURLISCH-STOER method    *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point NTOP           *
C   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
C   *   and    R/I standing for regular/irregular solution             *
C   *                                                                  *
C   *   NOTE: the routine expects the potential terms V and B          *
C   *         to be COMPLEX                                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,JLAG1
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--DIRBSBI27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLAG,MPSMAX,NPEMAX,NCFMAX
      PARAMETER (NLAG=3,MPSMAX=40,NPEMAX=4,NCFMAX=8)
C
C Dummy arguments
C
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ
      COMPLEX*16 E,PIS
      LOGICAL GETIRRSOL
      INTEGER IRTOP,KAP1,KAP2,L,NPAN,NRMAX,Z
      REAL*8 AP(2,2,NRMAX),AQ(2,2,NRMAX),DOVR(NRMAX),DRDI(NRMAX),
     &       R(NRMAX)
      COMPLEX*16 B(NRMAX),DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX),
     &           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX),V(NRMAX)
      INTEGER IRCUT(0:NPAN)
C
C Local variables
C
      REAL*8 A11,A12,A21,A22,CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),GAM(2)
     &       ,GPM,HBS,RPWGPM,RR,SK(2),SK1,SK2,TZ,X
      COMPLEX*16 AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BC(0:NPEMAX),BETA,
     &           BQQ,CFAC,DETD,DY(NCFMAX),DYSAV(:),EFAC,EMVPP,EMVQQ,
     &           PC(2,2,-NPEMAX:MPSMAX),QC(2,2,-NPEMAX:MPSMAX),
     &           VC(0:NPEMAX),W1,W2,W3,W4,W5,W6,W7,Y(NCFMAX),YERR(:),
     &           YM(:),YN(:),YSAV(:),YSEQ(:),ZZ
      COMPLEX*16 CJLZ
      INTEGER I,IP,IPAN,IR,IR_LOWER_BOUND,ISK1,ISK2,IV,J,K,LB(2),LB1,
     &        LB2,M,MPS,NPE,NSOL,NY
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE YM,YN,YERR,YSAV,YSEQ,DYSAV
      ALLOCATE (YERR(NCFMAX),YSAV(NCFMAX),YSEQ(NCFMAX),DYSAV(NCFMAX))
      ALLOCATE (YM(NCFMAX),YN(NCFMAX))
C
      CALL CINIT(NCFMAX,YM)
      CALL CINIT(NCFMAX,YN)
      CALL CINIT(NCFMAX,YERR)
      CALL CINIT(NCFMAX,YSAV)
      CALL CINIT(NCFMAX,YSEQ)
      CALL CINIT(NCFMAX,DYSAV)
C
      CSQR = C*C
      CFAC = PIS*C/(E+CSQR)
C
C ----- find   NPE  expansion coefficients for the potential and B-field
C
      NPE = 4
C
      TZ = DBLE(2*Z)
C
      DO IV = 1,NPE
         DO IR = 1,NPE
            CM(IR,IV) = R(IR)**(IV-1)
         END DO
      END DO
C
      CALL RINVGJ(CMI,CM,NPEMAX,NPE)
C
      DO IV = 1,NPE
         VC(IV-1) = 0.0D0
         DO IR = 1,NPE
            IF ( (.NOT.FINITE_NUCLEUS) .OR. (Z.EQ.0) ) THEN
               VC(IV-1) = VC(IV-1) + CMI(IV,IR)*(V(IR)+TZ/R(IR))
            ELSE
               VC(IV-1) = VC(IV-1) + CMI(IV,IR)*V(IR)
            END IF
         END DO
      END DO
C
      DO IV = 1,NPE
         BC(IV-1) = 0.0D0
         DO IR = 1,NPE
            BC(IV-1) = BC(IV-1) + CMI(IV,IR)*B(IR)
         END DO
      END DO
C
C    calculate g-coefficients of b-field
C
      ISK1 = ISIGN(1,KAP1)
      ISK2 = ISIGN(1,KAP2)
      SK1 = DBLE(ISK1)
      SK2 = DBLE(ISK2)
      LB1 = L - ISK1
      LB2 = L - ISK2
C
      CG1 = -MJ/(KAP1+0.5D0)
      CG5 = -MJ/(-KAP1+0.5D0)
      CGD(1) = CG1
      CGMD(1) = CG5
      KAP(1) = DBLE(KAP1)
C
      IF ( (.NOT.FINITE_NUCLEUS) .OR. (Z.EQ.0) ) THEN
         GAM(1) = DSQRT(KAP(1)**2-(TZ/C)**2)
      ELSE
         GAM(1) = DABS(KAP(1))
      END IF
C
      LB(1) = LB1
      SK(1) = SK1
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
         LB(2) = 0
         SK(2) = 0.0D0
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
         IF ( (.NOT.FINITE_NUCLEUS) .OR. (Z.EQ.0) ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/C)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
C
         LB(2) = LB2
         SK(2) = SK2
      END IF
C
      NSOLBS = NSOL
      EBS = E
C
      DO I = 1,2
         DO J = 1,2
            DO IP = -NPEMAX,MPSMAX
               PC(I,J,IP) = C0
               QC(I,J,IP) = C0
            END DO
         END DO
      END DO
C
C ======================================================================
      IF ( Z.NE.0 ) THEN
C
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KAP(J))-GAM(J))
            QC(J,J,0) = (KAP(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = C0
            QC(I,J,0) = C0
         END DO
C
C  determine higher expansion coefficients for the wave functions
C
         MPS = 40
C
         AA12 = -TZ/CSQR
         AA21 = TZ
         EMVQQ = (E-VC(0)+CSQR)/CSQR
         EMVPP = -E + VC(0)
         BQQ = BC(0)/CSQR
C
C-----------------------------------------------------------------------
C
         IF ( .NOT.FINITE_NUCLEUS ) THEN
C
            DO J = 1,NSOL
C
               DO M = 1,MPS
                  DO I = 1,NSOL
                     K = 3 - I
                     BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                     BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                     *CGO*PC(K,J,M-1)
                     DO IP = 1,NPE - 1
                        BB1 = BB1 + (-VC(IP)+BC(IP)*CGMD(I))
     &                        *QC(I,J,M-1-IP)/CSQR
                        BB2 = BB2 + (+VC(IP)+BC(IP)*CGD(I))
     &                        *PC(I,J,M-1-IP) + (+BC(IP)*CGO)
     &                        *PC(K,J,M-1-IP)
                     END DO
C
                     AA11 = GAM(J) + M + KAP(I)
                     AA22 = GAM(J) + M - KAP(I)
                     DETD = AA11*AA22 - AA12*AA21
                     PC(I,J,M) = (BB1*AA22-AA12*BB2)/DETD
                     QC(I,J,M) = (AA11*BB2-BB1*AA21)/DETD
                  END DO
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C expansion adapted for potentials with finite nucleus
C expansion of potential up to second order: V_O+V_1*R+V_2*R*R
C
         ELSE
C
            DO J = 1,NSOL
               I = 3 - J
C -------------------------------------------- arbitrary starting values
               IF ( KAP(J).GT.0 ) THEN
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
     &                 = (-W3*PC(I,J,0)+W4*PC(3-I,J,0))/A12
               END DO
            END DO
            DO J = 1,NSOL
               DO I = 1,NSOL
                  W1 = EMVQQ + BQQ*CGMD(I)
                  W3 = -EMVPP + BC(0)*CGD(I)
                  A11 = GAM(J) + KAP(I) + 2D0
                  A12 = GAM(J) - KAP(I) + 2D0
                  IF ( RNON0(A11) ) PC(I,J,2)
     &                 = (W1*QC(I,J,1)-W2*QC(I,J,0))/A11
                  IF ( RNON0(A12) ) QC(I,J,2)
     &                 = (-W3*PC(I,J,1)+W4*PC(3-I,J,1)+W5*PC(I,J,0))/A12
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
     &                    = (W1*QC(I,J,M-1)-W2*QC(I,J,M-2)
     &                    -W6*QC(I,J,M-3))/A21
                     IF ( RNON0(A22) ) QC(I,J,M)
     &                    = (-W3*PC(I,J,M-1)+W4*PC(3-I,J,M-1)
     &                    +W5*PC(I,J,M-2)+W7*PC(I,J,M-3))/A22
                  END DO
               END DO
            END DO
         END IF
C
C -------- perform summation over wave function - expansion coefficients
C --------------------------------------- for the first r - mesh - point
C
         IR = 1
         RR = R(IR)
C
         DO J = 1,NSOL
            RPWGPM = RR**GAM(J)
C
            DO I = 1,NSOL
               PR(I,J,IR) = PC(I,J,0)*RPWGPM
               QR(I,J,IR) = QC(I,J,0)*RPWGPM
               DP(I,J,IR) = PC(I,J,0)*RPWGPM*GAM(J)*DOVR(IR)
               DQ(I,J,IR) = QC(I,J,0)*RPWGPM*GAM(J)*DOVR(IR)
            END DO
C
            DO M = 1,MPS
               RPWGPM = RPWGPM*RR
               GPM = GAM(J) + M
C
               DO I = 1,NSOL
                  PR(I,J,IR) = PR(I,J,IR) + PC(I,J,M)*RPWGPM
                  QR(I,J,IR) = QR(I,J,IR) + QC(I,J,M)*RPWGPM
                  DP(I,J,IR) = DP(I,J,IR) + PC(I,J,M)
     &                         *RPWGPM*GPM*DOVR(IR)
                  DQ(I,J,IR) = DQ(I,J,IR) + QC(I,J,M)
     &                         *RPWGPM*GPM*DOVR(IR)
               END DO
C
            END DO
         END DO
C
C ======================================================================
C                                                     == EMPTY SPHERE ==
      ELSE
C                     assume constant pot: V=V(1)   ignore coupling: B=0
C
         IR = 1
         ZZ = CDSQRT(E-V(1))*R(IR)
         EFAC = (ZZ/R(IR))*C/(E+CSQR)
C
         DO J = 1,NSOL
            I = 3 - J
            PR(J,J,IR) = CJLZ(L,ZZ)*R(IR)
            QR(J,J,IR) = EFAC*SK(J)*CJLZ(LB(J),ZZ)*R(IR)*C
            DP(J,J,IR) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))*DRDI(IR)
            M = LB(J)
            DQ(J,J,IR) = EFAC*SK(J)
     &                   *(DBLE(M+1)*CJLZ(M,ZZ)-ZZ*CJLZ(M+1,ZZ))
     &                   *DRDI(IR)*C
C
            PR(I,J,IR) = C0
            QR(I,J,IR) = C0
            DP(I,J,IR) = C0
            DQ(I,J,IR) = C0
         END DO
C
      END IF
C
C ======================================================================
C
      NY = 0
      DO J = 1,NSOL
         DO I = 1,NSOL
            Y(NY+1) = PR(I,J,1)
            Y(NY+2) = QR(I,J,1)
            NY = NY + 2
         END DO
      END DO
C
      X = 1.0D0
      HBS = 1.0D0
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = 1,NPAN
C
         CALL DIRBSRADBI(X,Y,DY,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C
         IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN-1) + 2,MIN(IRCUT(IPAN),IRTOP)
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR-(NLAG+1)/2+1)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL DIRBSSTPBI(Y,DY,NY,X,HBS,B,V,R,DRDI,YM,YN,YERR,YSAV,
     &                      YSEQ,DYSAV,NRMAX,NCFMAX,AP,AQ)
C
            NY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PR(I,J,IR) = Y(NY+1)
                  QR(I,J,IR) = Y(NY+2)
                  NY = NY + 2
               END DO
            END DO
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------------------------------ to first mesh point of next panel
C -------------------------- derivatives DRDI have to be calculated anew
C ----------------------- because DRDI jumps when crossing panel borders
C ---------------------- update X for start of integration in next panel
C
         IF ( IPAN.NE.NPAN ) THEN
            IR = IRCUT(IPAN) + 1
            NY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PR(I,J,IR) = Y(NY+1)
                  QR(I,J,IR) = Y(NY+2)
                  NY = NY + 2
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C ######################################################################
C                          irregular solution
C ######################################################################
C
C     calculate initial values of the wavefunction at sphere boundary
C
      IR = IRTOP
C
      ZZ = PIS*R(IR)
C
      DO J = 1,NSOL
         PI(J,J,IR) = CJLZ(L,ZZ)*R(IR)
         QI(J,J,IR) = CFAC*SK(J)*CJLZ(LB(J),ZZ)*R(IR)*C
         DP(J,J,IR) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))*DRDI(IR)
C
         M = LB(J)
         DQ(J,J,IR) = CFAC*SK(J)*(DBLE(M+1)*CJLZ(M,ZZ)-ZZ*CJLZ(M+1,ZZ))
     &                *DRDI(IR)*C
C
         I = 3 - J
         PI(I,J,IR) = C0
         QI(I,J,IR) = C0
         DP(I,J,IR) = C0
         DQ(I,J,IR) = C0
      END DO
C
C =============================================================== IR ===
C
      NY = 0
      DO J = 1,NSOL
         DO I = 1,NSOL
            Y(NY+1) = PI(I,J,IRTOP)
            Y(NY+2) = QI(I,J,IRTOP)
            NY = NY + 2
         END DO
      END DO
C
      X = DBLE(IRTOP)
      HBS = -1.0D0
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = NPAN,1, - 1
C
         IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
C
         CALL DIRBSRADBI(X,Y,DY,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN) - 1,IRCUT(IPAN-1) + 1, - 1
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL DIRBSSTPBI(Y,DY,NY,X,HBS,B,V,R,DRDI,YM,YN,YERR,YSAV,
     &                      YSEQ,DYSAV,NRMAX,NCFMAX,AP,AQ)
C
            NY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PI(I,J,IR) = Y(NY+1)
                  QI(I,J,IR) = Y(NY+2)
                  NY = NY + 2
               END DO
            END DO
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------ to first mesh point of next panel (except for last panel)
C ---------------------- update X for start of integration in next panel
C
         IF ( IPAN.NE.1 ) THEN
            IR = IRCUT(IPAN-1)
            NY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PI(I,J,IR) = Y(NY+1)
                  QI(I,J,IR) = Y(NY+2)
                  NY = NY + 2
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      END
C*==dirbsstpbi.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSSTPBI(Y,DY,NY,X,HBS,B,V,R,DRDI,YM,YN,YERR,YSAV,
     &                      YSEQ,DYSAV,NRMAX,NCFMAX,AP,AQ)
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DY    for last mesh-point                        *
C   *   on exit:  X,Y,DY    updated for X = X(last) + HBS              *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *         no step size adjusted in case of no convergency > STOP   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:TOL_DIRBS
      IMPLICIT NONE
C*--DIRBSSTPBI540
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ISEQMAX,NUSE
      PARAMETER (ISEQMAX=36,NUSE=7)
      REAL*8 TINYBIT
      PARAMETER (TINYBIT=1.0D-20)
C
C Dummy arguments
C
      REAL*8 HBS,X
      INTEGER NCFMAX,NRMAX,NY
      REAL*8 AP(2,2,NRMAX),AQ(2,2,NRMAX),DRDI(NRMAX),R(NRMAX)
      COMPLEX*16 B(NRMAX),DY(NY),DYSAV(NCFMAX),V(NRMAX),Y(NY),
     &           YERR(NCFMAX),YM(NCFMAX),YN(NCFMAX),YSAV(NCFMAX),
     &           YSEQ(NCFMAX)
C
C Local variables
C
      REAL*8 ABSBB,ERRMAX,FXX(NUSE),H2MID,HMID,XEST,XMID,XX(ISEQMAX)
      COMPLEX*16 BB,BB1,CC,CRAT,DD(NCFMAX,NUSE),DYY,SWAP,VV,YY
      INTEGER I,ISTEP,IY,K,M1,NSEQ(ISEQMAX)
C
C*** End of declarations rewritten by SPAG
C
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,
     &     1024,1536,2048,3072,4096,6144,8192,12288,16384,24576,32768,
     &     49152,65536,98304,131072,196608,262144,393216,524288/
C
C
      DO IY = 1,NY
         YSAV(IY) = Y(IY)
         DYSAV(IY) = DY(IY)
      END DO
C
      DO I = 1,ISEQMAX
C
C------------------------------------------------------------------- MID
C
         HMID = HBS/NSEQ(I)
         DO IY = 1,NY
            YM(IY) = YSAV(IY)
            YN(IY) = YSAV(IY) + HMID*DYSAV(IY)
         END DO
         XMID = X + HMID
C
         CALL DIRBSRADBI(XMID,YN,YSEQ,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C
         H2MID = 2.D0*HMID
         DO ISTEP = 2,NSEQ(I)
            DO IY = 1,NY
               SWAP = YM(IY) + H2MID*YSEQ(IY)
               YM(IY) = YN(IY)
               YN(IY) = SWAP
            END DO
            XMID = XMID + HMID
            CALL DIRBSRADBI(XMID,YN,YSEQ,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
         END DO
         DO IY = 1,NY
            YSEQ(IY) = 0.5D0*(YM(IY)+YN(IY)+HMID*YSEQ(IY))
         END DO
C
C-----------------------------------------------------------------------
         XEST = HBS/NSEQ(I)
         XEST = XEST*XEST
C------------------------------------------------------------------- RZE
C
         XX(I) = XEST
         IF ( I.EQ.1 ) THEN
            DO IY = 1,NY
               Y(IY) = YSEQ(IY)
               DD(IY,1) = YSEQ(IY)
               YERR(IY) = YSEQ(IY)
            END DO
         ELSE
            M1 = MIN(I,NUSE)
            DO K = 1,M1 - 1
               FXX(K+1) = XX(I-K)/XEST
            END DO
            DO IY = 1,NY
               YY = YSEQ(IY)
               VV = DD(IY,1)
               CC = YY
               DD(IY,1) = YY
               DO K = 2,M1
                  BB1 = FXX(K)*VV
                  BB = BB1 - CC
                  ABSBB = ABS(DREAL(BB)) + ABS(DIMAG(BB))
                  IF ( ABSBB.GT.1D-40 ) THEN
                     BB = (CC-VV)/BB
                     DYY = CC*BB
                     CC = BB1*BB
                  ELSE
                     DYY = VV
                  END IF
                  VV = DD(IY,K)
                  DD(IY,K) = DYY
                  YY = YY + DYY
               END DO
               YERR(IY) = DYY
               Y(IY) = YY
            END DO
         END IF
C
C-----------------------------------------------------------------------
C
         DO IY = 1,NY
            CRAT = YERR(IY)/(Y(IY)+TINYBIT)
            ERRMAX = ABS(DREAL(CRAT)) + ABS(DIMAG(CRAT))
            IF ( ERRMAX.GT.TOL_DIRBS ) GOTO 100
         END DO
C
C------------------------------------------------------------- converged
C
         X = X + HBS
C
         CALL DIRBSRADBI(X,Y,DY,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<DIRBSSTP>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error : ',ERRMAX
      WRITE (6,*) 'tolerance TOL_DIRBS : ',TOL_DIRBS
      WRITE (6,*) 'grid position     X : ',X
      X = X + HBS
C
      CALL DIRBSRADBI(X,Y,DY,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C
      END
C*==dirbsradbi.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSRADBI(XBS,Y,DYDX,DRDI,B,V,R,NRMAX,NCFMAX,AP,AQ)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,JLAG1
      IMPLICIT NONE
C*--DIRBSRADBI703
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLAG
      PARAMETER (NLAG=3)
C
C Dummy arguments
C
      INTEGER NCFMAX,NRMAX
      REAL*8 XBS
      REAL*8 AP(2,2,NRMAX),AQ(2,2,NRMAX),DRDI(NRMAX),R(NRMAX)
      COMPLEX*16 B(NRMAX),DYDX(NCFMAX),V(NRMAX),Y(NCFMAX)
C
C Local variables
C
      COMPLEX*16 APBS(2,2),AQBS(2,2),BBS,BPP,BPP_CGO,BQQ,EMVPP,
     &           EMVPP_BPP_CGD(2),EMVQQ,EMVQQ_BQQ_CGMD(2),VBS
      REAL*8 DRDIBS,DROVR,DX,DXSQ,KAP_DROVR(2),RBS,W(NLAG)
      INTEGER I,IS,J,JS,K,KI,M
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-12 ) THEN
         I = NINT(XBS)
         VBS = V(I)
         BBS = B(I)
         RBS = R(I)
         DRDIBS = DRDI(I)
         DO IS = 1,2
            DO JS = 1,2
               APBS(IS,JS) = AP(IS,JS,I)
               AQBS(IS,JS) = AQ(IS,JS,I)
            END DO
         END DO
      ELSE
         DX = XBS - JLAG1
         DXSQ = DX*DX
         W(3) = 0.5D0*(-DX+DXSQ)
         W(1) = 1D0 - DX + W(3)
         W(2) = DX + DX - DXSQ
         VBS = 0D0
         BBS = 0D0
         RBS = 0D0
         APBS(1:2,1:2) = 0.0D0
         AQBS(1:2,1:2) = 0.0D0
         DRDIBS = 0D0
         J = JLAG1 - 1
         DO I = 1,NLAG
            J = J + 1
            VBS = VBS + W(I)*V(J)
            BBS = BBS + W(I)*B(J)
            RBS = RBS + W(I)*R(J)
            DRDIBS = DRDIBS + W(I)*DRDI(J)
            DO IS = 1,2
               DO JS = 1,2
                  APBS(IS,JS) = APBS(IS,JS) + W(I)*AP(IS,JS,J)
                  AQBS(IS,JS) = AQBS(IS,JS) + W(I)*AQ(IS,JS,J)
               END DO
            END DO
         END DO
      END IF
C-----------------------------------------------------------------------
C
      EMVPP = -(EBS-VBS)*DRDIBS
CCCCC EMVQQ = (EBS-VBS+CSQR)*DRDIBS/CSQRD
      EMVQQ = -EMVPP/CSQR + DRDIBS
C
      BPP = BBS*DRDIBS
CCCCC BQQ = BBS*DRDIBS/CSQR
      BQQ = BPP/CSQR
C
      DROVR = DRDIBS/RBS
C
      DO I = 1,NSOLBS
         KAP_DROVR(I) = KAP(I)*DROVR
         EMVQQ_BQQ_CGMD(I) = EMVQQ + BQQ*CGMD(I)
         EMVPP_BPP_CGD(I) = EMVPP + BPP*CGD(I)
      END DO
      BPP_CGO = BPP*CGO
      DO IS = 1,2
         DO JS = 1,2
            APBS(IS,JS) = APBS(IS,JS)*DRDIBS
            AQBS(IS,JS) = AQBS(IS,JS)*DRDIBS
         END DO
      END DO
C
      M = 0
      DO J = 1,NSOLBS
         DO I = 1,NSOLBS
            K = 3 - 4*(I-1)
            KI = 3 - I
C
CCCCC       DYDX(M+1) = -KAP(I)*Y(M+1)/RBS*DRDIBS + (EMVQQ+BQQ*CGMD(I))
CCCCC       DYDX(M+1) = -KAP(I)*Y(M+1)*DROVR + (EMVQQ+BQQ*CGMD(I))
CCCCC&                  *Y(M+2)
C
            DYDX(M+1) = -KAP_DROVR(I)*Y(M+1) + EMVQQ_BQQ_CGMD(I)*Y(M+2)
     &                  + APBS(I,I)*Y(M+1) + APBS(I,KI)*Y(M+K)
C
CCCCC       DYDX(M+2) = KAP(I)*Y(M+2)/RBS*DRDIBS + (EMVPP+BPP*CGD(I))
CCCCC       DYDX(M+2) = KAP(I)*Y(M+2)*DROVR + (EMVPP+BPP*CGD(I))*Y(M+1)
CCCCC&                  + BPP*CGO*Y(M+K)
C
            DYDX(M+2) = KAP_DROVR(I)*Y(M+2) + EMVPP_BPP_CGD(I)*Y(M+1)
     &                  + BPP_CGO*Y(M+K) + AQBS(I,I)*Y(M+2) + AQBS(I,KI)
     &                  *Y(M+K+1)
C
            M = M + 2
         END DO
      END DO
      END
