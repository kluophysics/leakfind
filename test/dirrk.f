C*==dirrk.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DIRRK(GETIRRSOL,C,E,L,MJ,KAP1,KAP2,PIS,CG1,CG2,CG4,CG5,
     &                 CG8,V,B,Z,R,DRDI,DOVR,IRTOP,IRCUT,NPAN,PR,QR,PI,
     &                 QI,DP,DQ,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and the inward integration is started analytically             *
C   *   the integration itself is done by the RUNGE-KUTTA method       *
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
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--DIRRK25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,NCFMAX
      PARAMETER (MPSMAX=40,NPEMAX=4,NCFMAX=8)
C
C Dummy arguments
C
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ
      COMPLEX*16 E,PIS
      LOGICAL GETIRRSOL
      INTEGER IRTOP,KAP1,KAP2,L,NPAN,NRMAX,Z
      COMPLEX*16 B(NRMAX),DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX),
     &           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX),V(NRMAX)
      REAL*8 DOVR(NRMAX),DRDI(NRMAX),R(NRMAX)
      INTEGER IRCUT(0:NPAN)
C
C Local variables
C
      REAL*8 A11,A12,A21,A22,CGD(2),CGMD(2),CGO,CM(NPEMAX,NPEMAX),
     &       CMI(NPEMAX,NPEMAX),CSQR,DRDIRK,DROVR,FDR(:),FDRMID(:),
     &       FOV(:),FOVMID(:),GAM(2),GPM,KAP(2),RPWGPM,RR,RRK,RY1,RY2,
     &       RY3,RY4,SK(2),SK1,SK2,TZ,XPP(:),XPPMID(:),XQQ(:),XQQMID(:),
     &       YPP(:),YPPMID(:),YQQ(:),YQQMID(:)
      COMPLEX*16 AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BC(0:NPEMAX),BETA,
     &           BPP,BQQ,CFAC,CPL(:),CPLMID(:),CY1,CY2,CY3,CY4,DETD,
     &           DYSAV(:),EFAC,EMVPP,EMVQQ,FBB(:),FBBMID(:),FVV(:),
     &           FVVMID(:),PC(2,2,-NPEMAX:MPSMAX),QC(2,2,-NPEMAX:MPSMAX)
     &           ,VC(0:NPEMAX),W1,W2,W3,W4,W5,W6,W7,XPQ(:),XPQMID(:),
     &           XQP(:),XQPMID(:),YERR(:),YM(:),YN(:),YPQ(:),YPQMID(:),
     &           YQP(:),YQPMID(:),YSAV(:),YSEQ(:),ZZ
      COMPLEX*16 CINTPOL1,CINTPOL2,CINTPOL3
      COMPLEX*16 CJLZ
      INTEGER I,IP,IPAN,IR,IRLIM1,IRLIM2,ISK1,ISK2,ISOL,IV,J,K,LB(2),
     &        LB1,LB2,M,MPS,NPE,NSOL
      REAL*8 RINTPOL1,RINTPOL2,RINTPOL3
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C             in line functions for interpolation
C-----------------------------------------------------------------------
      RINTPOL1(RY1,RY2,RY3) = RY2 + (-RY3+RY1+RY1+RY1-RY2-RY2)/8D0
      RINTPOL3(RY1,RY2,RY3) = RY2 + (RY3+RY3+RY3-RY1-RY2-RY2)/8D0
      RINTPOL2(RY1,RY2,RY3,RY4) = (-RY1+9D0*(RY2+RY3)-RY4)/16D0
C
      CINTPOL1(CY1,CY2,CY3) = CY2 + (-CY3+CY1+CY1+CY1-CY2-CY2)/8D0
      CINTPOL3(CY1,CY2,CY3) = CY2 + (CY3+CY3+CY3-CY1-CY2-CY2)/8D0
      CINTPOL2(CY1,CY2,CY3,CY4) = (-CY1+9D0*(CY2+CY3)-CY4)/16D0
C-----------------------------------------------------------------------
C
      ALLOCATABLE XPP,XPPMID,XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID
      ALLOCATABLE YPP,YPPMID,YQQ,YQQMID,YPQ,YPQMID,YQP,YQPMID
      ALLOCATABLE FVV,FBB,FVVMID,FBBMID,FDR,FOV,FDRMID,FOVMID
      ALLOCATABLE CPL,CPLMID
      ALLOCATABLE YM,YN,YERR,YSAV,YSEQ,DYSAV
C
      ALLOCATE (YERR(NCFMAX),YSAV(NCFMAX),YSEQ(NCFMAX),DYSAV(NCFMAX))
      ALLOCATE (YM(NCFMAX),YN(NCFMAX))
C
      ALLOCATE (FVV(NRMAX),FBB(NRMAX),FVVMID(NRMAX),FBBMID(NRMAX))
      ALLOCATE (FDR(NRMAX),FOV(NRMAX),FDRMID(NRMAX),FOVMID(NRMAX))
C
      CALL CINIT(NCFMAX,YM)
      CALL CINIT(NCFMAX,YN)
      CALL CINIT(NCFMAX,YERR)
      CALL CINIT(NCFMAX,YSAV)
      CALL CINIT(NCFMAX,YSEQ)
      CALL CINIT(NCFMAX,DYSAV)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      DO IR = 1,IRTOP
C
         RRK = R(IR)
         DRDIRK = DRDI(IR)
C
         FDR(IR) = DRDIRK
         FOV(IR) = DRDIRK/RRK
         FVV(IR) = -(E-V(IR))*DRDIRK
         FBB(IR) = B(IR)*DRDIRK
C
      END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C-----------------------------------------------------------------------
C                  interpolate for the mid point positions
C-----------------------------------------------------------------------
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = 1,NPAN
C
         IRLIM1 = IRCUT(IPAN-1) + 1
         IRLIM2 = MIN(IRCUT(IPAN),IRTOP) - 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN-1) + 1,MIN(IRCUT(IPAN),IRTOP) - 1
C
            IF ( IR.EQ.IRLIM1 ) THEN
C
               FDRMID(IR) = RINTPOL1(FDR(IR),FDR(IR+1),FDR(IR+2))
               FOVMID(IR) = RINTPOL1(FOV(IR),FOV(IR+1),FOV(IR+2))
               FVVMID(IR) = CINTPOL1(FVV(IR),FVV(IR+1),FVV(IR+2))
               FBBMID(IR) = CINTPOL1(FBB(IR),FBB(IR+1),FBB(IR+2))
C
            ELSE IF ( IR.LT.IRLIM2 ) THEN
C
               FDRMID(IR) = RINTPOL2(FDR(IR-1),FDR(IR),FDR(IR+1),
     &                      FDR(IR+2))
               FOVMID(IR) = RINTPOL2(FOV(IR-1),FOV(IR),FOV(IR+1),
     &                      FOV(IR+2))
               FVVMID(IR) = CINTPOL2(FVV(IR-1),FVV(IR),FVV(IR+1),
     &                      FVV(IR+2))
               FBBMID(IR) = CINTPOL2(FBB(IR-1),FBB(IR),FBB(IR+1),
     &                      FBB(IR+2))
C
            ELSE
C
               FDRMID(IR) = RINTPOL3(FDR(IR-1),FDR(IR),FDR(IR+1))
               FOVMID(IR) = RINTPOL3(FOV(IR-1),FOV(IR),FOV(IR+1))
               FVVMID(IR) = CINTPOL3(FVV(IR-1),FVV(IR),FVV(IR+1))
               FBBMID(IR) = CINTPOL3(FBB(IR-1),FBB(IR),FBB(IR+1))
C
            END IF
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
         IR = MIN(IRCUT(IPAN),IRTOP)
         FDRMID(IR) = 0D0
         FOVMID(IR) = 0D0
         FVVMID(IR) = 0D0
         FBBMID(IR) = 0D0
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
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
C
C ######################################################################
C                       irregular solution
C ######################################################################
C
C         calculate the initial values of the wavefunction
C                     at the sphere boundary
C
      IF ( GETIRRSOL ) THEN
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
            DQ(J,J,IR) = CFAC*SK(J)
     &                   *(DBLE(M+1)*CJLZ(M,ZZ)-ZZ*CJLZ(M+1,ZZ))
     &                   *DRDI(IR)*C
C
            I = 3 - J
            PI(I,J,IR) = C0
            QI(I,J,IR) = C0
            DP(I,J,IR) = C0
            DQ(I,J,IR) = C0
         END DO
C
      END IF
C
C ======================================================================
C
      ALLOCATE (XPP(NRMAX),XPPMID(NRMAX),XQQ(NRMAX),XQQMID(NRMAX))
      ALLOCATE (XPQ(NRMAX),XPQMID(NRMAX),XQP(NRMAX),XQPMID(NRMAX))
C
C ======================================================================
C                   NSOL = 1   NO coupling of wave functions
C ======================================================================
      IF ( NSOL.EQ.1 ) THEN
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = 1,IRTOP
C
C--------------------------------------------------- regular mesh points
            EMVPP = FVV(IR)
            EMVQQ = -EMVPP/CSQR + FDR(IR)
C
            BPP = FBB(IR)
            BQQ = BPP/CSQR
C
            DROVR = FOV(IR)
C
            XPP(IR) = -KAP(1)*DROVR
            XPQ(IR) = EMVQQ + BQQ*CGMD(1)
            XQP(IR) = EMVPP + BPP*CGD(1)
            XQQ(IR) = -XPP(IR)
C
C------------------------------------------------------------ mid points
            EMVPP = FVVMID(IR)
            EMVQQ = -EMVPP/CSQR + FDRMID(IR)
C
            BPP = FBBMID(IR)
            BQQ = BPP/CSQR
C
            DROVR = FOVMID(IR)
C
            XPPMID(IR) = -KAP(1)*DROVR
            XPQMID(IR) = EMVQQ + BQQ*CGMD(1)
            XQPMID(IR) = EMVPP + BPP*CGD(1)
            XQQMID(IR) = -XPPMID(IR)
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
         CALL SOLVE_RK2(GETIRRSOL,IRTOP,IRCUT,NPAN,PR,QR,PI,QI,XPP,
     &                  XPPMID,XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID,2,NRMAX)
C
         RETURN
C
      ELSE
C ======================================================================
C                   NSOL = 2    coupling of wave functions
C ======================================================================
C
         ALLOCATE (YPP(NRMAX),YPPMID(NRMAX),YQQ(NRMAX),YQQMID(NRMAX))
         ALLOCATE (YPQ(NRMAX),YPQMID(NRMAX),YQP(NRMAX),YQPMID(NRMAX))
         ALLOCATE (CPL(NRMAX),CPLMID(NRMAX))
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = 1,IRTOP
C
C--------------------------------------------------- regular mesh points
            EMVPP = FVV(IR)
            EMVQQ = -EMVPP/CSQR + FDR(IR)
C
            BPP = FBB(IR)
            BQQ = BPP/CSQR
C
            DROVR = FOV(IR)
C
            XPP(IR) = -KAP(1)*DROVR
            XPQ(IR) = EMVQQ + BQQ*CGMD(1)
            XQP(IR) = EMVPP + BPP*CGD(1)
            XQQ(IR) = -XPP(IR)
C
            YPP(IR) = -KAP(2)*DROVR
            YPQ(IR) = EMVQQ + BQQ*CGMD(2)
            YQP(IR) = EMVPP + BPP*CGD(2)
            YQQ(IR) = -YPP(IR)
C
            CPL(IR) = BPP*CGO
C
C------------------------------------------------------------ mid points
            EMVPP = FVVMID(IR)
            EMVQQ = -EMVPP/CSQR + FDRMID(IR)
C
            BPP = FBBMID(IR)
            BQQ = BPP/CSQR
C
            DROVR = FOVMID(IR)
C
            XPPMID(IR) = -KAP(1)*DROVR
            XPQMID(IR) = EMVQQ + BQQ*CGMD(1)
            XQPMID(IR) = EMVPP + BPP*CGD(1)
            XQQMID(IR) = -XPPMID(IR)
C
            YPPMID(IR) = -KAP(2)*DROVR
            YPQMID(IR) = EMVQQ + BQQ*CGMD(2)
            YQPMID(IR) = EMVPP + BPP*CGD(2)
            YQQMID(IR) = -YPPMID(IR)
C
            CPLMID(IR) = BPP*CGO
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
         DO ISOL = 1,NSOL
C
            CALL SOLVE_RK4(ISOL,GETIRRSOL,IRTOP,IRCUT,NPAN,PR,QR,PI,QI,
     &                     XPP,XPPMID,XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID,
     &                     YPP,YPPMID,YQQ,YQQMID,YPQ,YPQMID,YQP,YQPMID,
     &                     CPL,CPLMID,NRMAX)
C
         END DO
C
      END IF
C
C ======================================================================
C
      END
C*==solve_rk4.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE SOLVE_RK4(ISOL,GETIRRSOL,IRTOP,IRCUT,NPAN,PR,QR,PI,QI,
     &                     XPP,XPPMID,XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID,
     &                     YPP,YPPMID,YQQ,YQQMID,YPQ,YPQMID,YQP,YQPMID,
     &                     CPL,CPLMID,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine to solve radial differential equations        *
C   *   using the  RUNGE-KUTTA method                                  *
C   *                                                                  *
C   *   solve 4 coupled equations for P and Q                          *
C   *                                                                  *
C   *   NOTE: the wave functions at the first mesh points have to be   *
C   *         supplied,  i.e.  PR(1), QR(1)  and  PR(IRTOP), QR(IRTOP) *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SOLVE_RK4628
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL GETIRRSOL
      INTEGER IRTOP,ISOL,NPAN,NRMAX
      COMPLEX*16 CPL(NRMAX),CPLMID(NRMAX),PI(2,2,NRMAX),PR(2,2,NRMAX),
     &           QI(2,2,NRMAX),QR(2,2,NRMAX),XPQ(NRMAX),XPQMID(NRMAX),
     &           XQP(NRMAX),XQPMID(NRMAX),YPQ(NRMAX),YPQMID(NRMAX),
     &           YQP(NRMAX),YQPMID(NRMAX)
      INTEGER IRCUT(0:NPAN)
      REAL*8 XPP(NRMAX),XPPMID(NRMAX),XQQ(NRMAX),XQQMID(NRMAX),
     &       YPP(NRMAX),YPPMID(NRMAX),YQQ(NRMAX),YQQMID(NRMAX)
C
C Local variables
C
      INTEGER IPAN,IR
      COMPLEX*16 K11P,K11Q,K12P,K12Q,K13P,K13Q,K14P,K14Q,K21P,K21Q,K22P,
     &           K22Q,K23P,K23Q,K24P,K24Q,P1AUX,P2AUX,Q1AUX,Q2AUX
C
C*** End of declarations rewritten by SPAG
C
C ######################################################################
C                           regular solution
C ######################################################################
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = 1,NPAN
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN-1) + 1,MIN(IRCUT(IPAN),IRTOP) - 1
C
            K11P = XPP(IR)*PR(1,ISOL,IR) + XPQ(IR)*QR(1,ISOL,IR)
            K11Q = XQP(IR)*PR(1,ISOL,IR) + XQQ(IR)*QR(1,ISOL,IR)
     &             + CPL(IR)*PR(2,ISOL,IR)
C
            K21P = YPP(IR)*PR(2,ISOL,IR) + YPQ(IR)*QR(2,ISOL,IR)
            K21Q = YQP(IR)*PR(2,ISOL,IR) + YQQ(IR)*QR(2,ISOL,IR)
     &             + CPL(IR)*PR(1,ISOL,IR)
C
            P1AUX = PR(1,ISOL,IR) + 0.5D0*K11P
            Q1AUX = QR(1,ISOL,IR) + 0.5D0*K11Q
C
            P2AUX = PR(2,ISOL,IR) + 0.5D0*K21P
            Q2AUX = QR(2,ISOL,IR) + 0.5D0*K21Q
C
C-----------------------------------------------------------------------
C
            K12P = XPPMID(IR)*P1AUX + XPQMID(IR)*Q1AUX
            K12Q = XQPMID(IR)*P1AUX + XQQMID(IR)*Q1AUX + CPLMID(IR)
     &             *P2AUX
C
            K22P = YPPMID(IR)*P2AUX + YPQMID(IR)*Q2AUX
            K22Q = YQPMID(IR)*P2AUX + YQQMID(IR)*Q2AUX + CPLMID(IR)
     &             *P1AUX
C
            P1AUX = PR(1,ISOL,IR) + 0.5D0*K12P
            Q1AUX = QR(1,ISOL,IR) + 0.5D0*K12Q
C
            P2AUX = PR(2,ISOL,IR) + 0.5D0*K22P
            Q2AUX = QR(2,ISOL,IR) + 0.5D0*K22Q
C
C-----------------------------------------------------------------------
C
            K13P = XPPMID(IR)*P1AUX + XPQMID(IR)*Q1AUX
            K13Q = XQPMID(IR)*P1AUX + XQQMID(IR)*Q1AUX + CPLMID(IR)
     &             *P2AUX
C
            K23P = YPPMID(IR)*P2AUX + YPQMID(IR)*Q2AUX
            K23Q = YQPMID(IR)*P2AUX + YQQMID(IR)*Q2AUX + CPLMID(IR)
     &             *P1AUX
C
            P1AUX = PR(1,ISOL,IR) + K13P
            Q1AUX = QR(1,ISOL,IR) + K13Q
C
            P2AUX = PR(2,ISOL,IR) + K23P
            Q2AUX = QR(2,ISOL,IR) + K23Q
C
C-----------------------------------------------------------------------
C
            K14P = XPP(IR+1)*P1AUX + XPQ(IR+1)*Q1AUX
            K14Q = XQP(IR+1)*P1AUX + XQQ(IR+1)*Q1AUX + CPL(IR+1)*P2AUX
C
            K24P = YPP(IR+1)*P2AUX + YPQ(IR+1)*Q2AUX
            K24Q = YQP(IR+1)*P2AUX + YQQ(IR+1)*Q2AUX + CPL(IR+1)*P1AUX
C
C-----------------------------------------------------------------------
C
            PR(1,ISOL,IR+1) = PR(1,ISOL,IR)
     &                        + (K11P+K12P+K12P+K13P+K13P+K14P)/6D0
            QR(1,ISOL,IR+1) = QR(1,ISOL,IR)
     &                        + (K11Q+K12Q+K12Q+K13Q+K13Q+K14Q)/6D0
C
            PR(2,ISOL,IR+1) = PR(2,ISOL,IR)
     &                        + (K21P+K22P+K22P+K23P+K23P+K24P)/6D0
            QR(2,ISOL,IR+1) = QR(2,ISOL,IR)
     &                        + (K21Q+K22Q+K22Q+K23Q+K23Q+K24Q)/6D0
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------------------------------ to first mesh point of next panel
C -------------------------- derivatives DRDI have to be calculated anew
C ----------------------- because DRDI jumps when crossing panel borders
C
         IF ( IPAN.NE.NPAN ) THEN
            IR = IRCUT(IPAN) + 1
            PR(1,ISOL,IR) = PR(1,ISOL,IR-1)
            QR(1,ISOL,IR) = QR(1,ISOL,IR-1)
            PR(2,ISOL,IR) = PR(2,ISOL,IR-1)
            QR(2,ISOL,IR) = QR(2,ISOL,IR-1)
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
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = NPAN,1, - 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN),IRCUT(IPAN-1) + 2, - 1
C
            K11P = XPP(IR)*PI(1,ISOL,IR) + XPQ(IR)*QI(1,ISOL,IR)
            K11Q = XQP(IR)*PI(1,ISOL,IR) + XQQ(IR)*QI(1,ISOL,IR)
     &             + CPL(IR)*PI(2,ISOL,IR)
C
            K21P = YPP(IR)*PI(2,ISOL,IR) + YPQ(IR)*QI(2,ISOL,IR)
            K21Q = YQP(IR)*PI(2,ISOL,IR) + YQQ(IR)*QI(2,ISOL,IR)
     &             + CPL(IR)*PI(1,ISOL,IR)
C
            P1AUX = PI(1,ISOL,IR) - 0.5D0*K11P
            Q1AUX = QI(1,ISOL,IR) - 0.5D0*K11Q
C
            P2AUX = PI(2,ISOL,IR) - 0.5D0*K21P
            Q2AUX = QI(2,ISOL,IR) - 0.5D0*K21Q
C
C-----------------------------------------------------------------------
C
            K12P = XPPMID(IR-1)*P1AUX + XPQMID(IR-1)*Q1AUX
            K12Q = XQPMID(IR-1)*P1AUX + XQQMID(IR-1)
     &             *Q1AUX + CPLMID(IR-1)*P2AUX
C
            K22P = YPPMID(IR-1)*P2AUX + YPQMID(IR-1)*Q2AUX
            K22Q = YQPMID(IR-1)*P2AUX + YQQMID(IR-1)
     &             *Q2AUX + CPLMID(IR-1)*P1AUX
C
            P1AUX = PI(1,ISOL,IR) - 0.5D0*K12P
            Q1AUX = QI(1,ISOL,IR) - 0.5D0*K12Q
C
            P2AUX = PI(2,ISOL,IR) - 0.5D0*K22P
            Q2AUX = QI(2,ISOL,IR) - 0.5D0*K22Q
C
C-----------------------------------------------------------------------
C
            K13P = XPPMID(IR-1)*P1AUX + XPQMID(IR-1)*Q1AUX
            K13Q = XQPMID(IR-1)*P1AUX + XQQMID(IR-1)
     &             *Q1AUX + CPLMID(IR-1)*P2AUX
C
            K23P = YPPMID(IR-1)*P2AUX + YPQMID(IR-1)*Q2AUX
            K23Q = YQPMID(IR-1)*P2AUX + YQQMID(IR-1)
     &             *Q2AUX + CPLMID(IR-1)*P1AUX
C
            P1AUX = PI(1,ISOL,IR) - K13P
            Q1AUX = QI(1,ISOL,IR) - K13Q
C
            P2AUX = PI(2,ISOL,IR) - K23P
            Q2AUX = QI(2,ISOL,IR) - K23Q
C
C-----------------------------------------------------------------------
C
            K14P = XPP(IR-1)*P1AUX + XPQ(IR-1)*Q1AUX
            K14Q = XQP(IR-1)*P1AUX + XQQ(IR-1)*Q1AUX + CPL(IR-1)*P2AUX
C
            K24P = YPP(IR-1)*P2AUX + YPQ(IR-1)*Q2AUX
            K24Q = YQP(IR-1)*P2AUX + YQQ(IR-1)*Q2AUX + CPL(IR-1)*P1AUX
C
C-----------------------------------------------------------------------
C
            PI(1,ISOL,IR-1) = PI(1,ISOL,IR)
     &                        - (K11P+K12P+K12P+K13P+K13P+K14P)/6D0
            QI(1,ISOL,IR-1) = QI(1,ISOL,IR)
     &                        - (K11Q+K12Q+K12Q+K13Q+K13Q+K14Q)/6D0
C
            PI(2,ISOL,IR-1) = PI(2,ISOL,IR)
     &                        - (K21P+K22P+K22P+K23P+K23P+K24P)/6D0
            QI(2,ISOL,IR-1) = QI(2,ISOL,IR)
     &                        - (K21Q+K22Q+K22Q+K23Q+K23Q+K24Q)/6D0
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------ to first mesh point of next panel (except for last panel)
C
         IF ( IPAN.NE.1 ) THEN
            IR = IRCUT(IPAN-1)
            PI(1,ISOL,IR) = PI(1,ISOL,IR+1)
            PI(1,ISOL,IR) = PI(1,ISOL,IR+1)
            QI(2,ISOL,IR) = QI(2,ISOL,IR+1)
            QI(2,ISOL,IR) = QI(2,ISOL,IR+1)
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      END
