C*==dirbs_soc.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBS_SOC(GETIRRSOL,C,SOCSCL,E,L,MJ,KAP1,KAP2,P,CG1,
     &                     CG2,CG4,CG5,CG8,V,B,Z,R,DRDI,DOVR,IRTOP,
     &                     IRCUT,NPAN,DXP,PR,QR,PI,QI,DP,DQ,NRMAX,
     &                     GF_CONV_RH)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *   with the SOC-operator manipulated                              *
C   *                                                                  *
C   *   SOCSCL >= 0: scaling the SPIN-ORBIT-COUPLING                   *
C   *                                                                  *
C   *   SOCSCL <  0:  H_soc = xi [ (sig_x*l_x+sig*l_y) + sig_z*l_z ]   *
C   *                                                                  *
C   *   SOCSCL = -1  ==  IXY = 0 : neglect 1st term of SOC-operator    *
C   *   SOCSCL = -2: ==  IXY = 1 : neglect 2nd term of SOC-operator    *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   the inward integration is started analytically                 *
C   *                                                                  *
C   *   BURLISCH-STOER version                                         *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point IRTOP          *
C   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
C   *   and    R/I standing for regular/irregular solution             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,CGD,CGMD,CGO,CGOZ,CGZ,LMK,KPX,
     &    NSOLBS,JLAG1,K_SOC
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--DIRBS_SOC35
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
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ,SOCSCL
      COMPLEX*16 E,P
      LOGICAL GETIRRSOL,GF_CONV_RH
      INTEGER IRTOP,KAP1,KAP2,L,NPAN,NRMAX,Z
      COMPLEX*16 B(NRMAX),DP(2,2,NRMAX),DQ(2,2,NRMAX),DXP(2,2),
     &           PI(2,2,NRMAX),PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX)
     &           ,V(NRMAX)
      REAL*8 DOVR(NRMAX),DRDI(NRMAX),R(NRMAX)
      INTEGER IRCUT(0:NPAN)
C
C Local variables
C
      COMPLEX*16 AA11,AA12,AA21,AA22,BB1,BB2,BC(0:NPEMAX),BQQ,CFAC,
     &           DBDR(:),DETD,DVDR(:),DY(:),DYSAV(:),EMVPP,EMVQQ,
     &           PC(2,2,-NPEMAX:MPSMAX),QC(2,2,-NPEMAX:MPSMAX),S0,
     &           SPHFUNL,SPHFUNLB,SPHFUNLBP1,SPHFUNLP1,T0,VC(0:NPEMAX),
     &           Y(:),YERR(:),YM(:),YN(:),YSAV(:),YSEQ(:),ZTOP,ZZ
      REAL*8 AR1(:),AR2(:),CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),GAM(2),
     &       GPM,HBS,KAP(2),KPY(2),RPWGPM,RR,SK(2),SK1,SK2,TZ,X
      COMPLEX*16 CJLZ,CNLZ
      INTEGER I,IA_ERR,IP,IPAN,IR,IR_LOWER_BOUND,ISK1,ISK2,IV,IXY,J,K,
     &        LB(2),LB1,LB2,LBJ,M,MPS,N,NPE,NSOL,NY
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DBDR,DVDR,AR1,AR2
      ALLOCATABLE Y,YM,YN,YERR,YSAV,YSEQ,DYSAV,DY
C
      ALLOCATE (YERR(NCFMAX),YSAV(NCFMAX),YSEQ(NCFMAX),DYSAV(NCFMAX))
      ALLOCATE (Y(NCFMAX),DY(NCFMAX),YM(NCFMAX),YN(NCFMAX))
C
      ALLOCATE (AR1(NRMAX),AR2(NRMAX))
      ALLOCATE (DBDR(NRMAX),DVDR(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:dirac_soc2 -> DVDR'
C
      YM(:) = 0D0
      YN(:) = 0D0
      YERR(:) = 0D0
      YSAV(:) = 0D0
      YSEQ(:) = 0D0
      DYSAV(:) = 0D0
C
      CSQR = C*C
      CFAC = P*C/(E+CSQR)
C
C find   NPE  expansion coefficients for the potential and b-field
      NPE = 4
C
      TZ = DBLE(2*Z)
C
      DO IV = 1,NPE
         DO N = 1,NPE
            CM(N,IV) = R(N)**(IV-1)
         END DO
      END DO
C
      CALL RINVGJ(CMI,CM,NPEMAX,NPE)
C
      DO IV = 1,NPE
         VC(IV-1) = 0.0D0
         DO N = 1,NPE
            VC(IV-1) = VC(IV-1) + CMI(IV,N)*(V(N)+TZ/R(N))
         END DO
      END DO
C
      DO IV = 1,NPE
         BC(IV-1) = 0.0D0
         DO N = 1,NPE
            BC(IV-1) = BC(IV-1) + CMI(IV,N)*B(N)
         END DO
      END DO
C
C------------------------ IXY=0 ==> only ZZ-part, IXY=1 ==> only XY-part
      IF ( NINT(SOCSCL).EQ.-1 ) THEN
         K_SOC = 2
         IXY = 0
      ELSE IF ( NINT(SOCSCL).EQ.-2 ) THEN
         K_SOC = 2
         IXY = 1
      ELSE
         K_SOC = 1
         IXY = 999999
      END IF
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
      GAM(1) = DSQRT(KAP(1)**2-(TZ/C)**2)
      CGZ(:) = 999999D0
      IF ( IXY.EQ.0 ) THEN
         CGZ(1) = KAP(1) + 1 + MJ*CGD(1) - 0.5D0
      ELSE IF ( IXY.EQ.1 ) THEN
         CGZ(1) = -(MJ*CGD(1)-0.5D0)
      END IF
      LB(1) = LB1
      SK(1) = SK1
      IF ( DABS(MJ).GT.L ) THEN
         CG2 = 0.0D0
         CG4 = 0.0D0
         CG8 = 0.0D0
         NSOL = 1
         CGD(2) = 0.0D0
         CGO = 0.0D0
         CGOZ = 0.0D0
         CGMD(2) = 0.0D0
         GAM(2) = 0.0D0
         KAP(2) = 0.0D0
         CGZ(2) = 0.0D0
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
         GAM(2) = DSQRT(KAP(2)**2-(TZ/C)**2)
         IF ( IXY.EQ.0 ) THEN
            CGZ(2) = KAP(2) + 1 + MJ*CGD(2) - 0.5D0
            CGOZ = MJ*CGO
         ELSE IF ( IXY.EQ.1 ) THEN
            CGZ(2) = -(MJ*CGD(2)-0.5D0)
            CGOZ = -MJ*CGO
         END IF
         LB(2) = LB2
         SK(2) = SK2
      END IF
C
      IF ( K_SOC.EQ.1 ) THEN
C
         DO I = 1,NSOL
            KPX(I) = -1.0D0 + SOCSCL*DBLE(1+KAP(I))
            LMK(I) = DBLE(L*(L+1)) - KPX(I)*(KPX(I)+1.0D0)
C
            CGMD(I) = 0.0D0
C-------------------------------------- causes numerical inconsistencies
C        GAM(I) = DSQRT( KPX(I)**2 - (TZ/C)**2 )
C        KPY(I) = KPX(I)
C
            GAM(I) = DSQRT(KAP(I)**2-(TZ/C)**2)
            KPY(I) = KAP(I)
C
         END DO
C
      ELSE
C
         DO I = 1,NSOL
            KPX(I) = KAP(I)
            KPY(I) = KAP(I)
         END DO
C
         AR1(1:IRTOP) = DREAL(V(1:IRTOP))
         CALL DVDRSPLINE(AR1,R,AR2,IRTOP)
         DVDR(1:IRTOP) = AR2(1:IRTOP)
C
         AR1(1:IRTOP) = DREAL(B(1:IRTOP))
         CALL DERSPL(IRTOP,R,AR1,AR2)
         DBDR(1:IRTOP) = AR2(1:IRTOP)
C
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
      IF ( (Z.NE.0) .AND. (.NOT.FINITE_NUCLEUS) ) THEN
C
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KPY(J))-GAM(J))
            QC(J,J,0) = (KPY(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = C0
            QC(I,J,0) = C0
         END DO
C
C  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
C
         MPS = 40
C
         AA12 = -TZ/CSQR
         AA21 = TZ
         EMVQQ = (E-VC(0)+CSQR)/CSQR
         EMVPP = -E + VC(0)
         BQQ = BC(0)/CSQR
C
         DO J = 1,NSOL
C
            DO M = 1,MPS
               DO I = 1,NSOL
                  K = 3 - I
                  BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                  BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                  *CGO*PC(K,J,M-1)
                  DO IP = 1,NPE - 1
                     BB1 = BB1 + (-VC(IP)+BC(IP)*CGMD(I))*QC(I,J,M-1-IP)
     &                     /CSQR
                     BB2 = BB2 + (+VC(IP)+BC(IP)*CGD(I))*PC(I,J,M-1-IP)
     &                     + BC(IP)*CGO*PC(K,J,M-1-IP)
                  END DO
C
                  AA11 = GAM(J) + M + KPY(I)
                  AA22 = GAM(J) + M - KPY(I)
                  DETD = AA11*AA22 - AA12*AA21
                  PC(I,J,M) = (BB1*AA22-AA12*BB2)/DETD
                  QC(I,J,M) = (AA11*BB2-BB1*AA21)/DETD
               END DO
            END DO
C
         END DO
C
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
C                                  == EMPTY SPHERE  or FINITE NUCLEUS ==
      ELSE
C
C        assume constant pot: V=V(1)   ignore coupling: B=0
C
         T0 = E - V(1)
         S0 = (E-V(1))/CSQR + 1D0
C
         IR = 1
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PR(I,J,IR) = C0
               QR(I,J,IR) = C0
               DP(I,J,IR) = C0
               DQ(I,J,IR) = C0
            END DO
         END DO
C
         ZZ = CDSQRT(S0*T0)*R(IR)
C
         DO J = 1,NSOL
            PR(J,J,IR) = CJLZ(L,ZZ)*R(IR)
            DP(J,J,IR) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))*DRDI(IR)
C
            QR(J,J,IR) = (DP(J,J,IR)/DRDI(IR)+PR(J,J,IR)*(KPX(J)/R(IR)))
     &                   /S0
            DQ(J,J,IR) = QR(J,J,IR)*(KPX(J)/R(IR)) - PR(J,J,IR)*T0
         END DO
C
      END IF
C ===================================================================
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
         CALL DIRBS_SOCRAD(X,Y,DY,DRDI,B,V,DBDR,DVDR,R,NRMAX,NCFMAX)
C
         IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN-1) + 2,MIN(IRCUT(IPAN),IRTOP)
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR-(NLAG+1)/2+1)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL DIRBS_SOCSTP(Y,DY,NY,X,HBS,B,V,DBDR,DVDR,R,DRDI,YM,YN,
     &                        YERR,YSAV,YSEQ,DYSAV,NRMAX,NCFMAX)
C
            NY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PR(I,J,IR) = Y(NY+1)
                  QR(I,J,IR) = Y(NY+2)
                  DP(I,J,IR) = DY(NY+1)
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
                  DP(I,J,IR) = DY(NY+1)
                  NY = NY + 2
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO I = 1,2
         DO J = 1,2
            DXP(I,J) = DP(I,J,IRTOP)
         END DO
      END DO
C
C     the minor component for the soc-manipulated wf is meaningless
C     =>  set it to zero
C
      QR(:,:,:) = C0
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
      ZTOP = P*R(IR)
C
      DO J = 1,NSOL
C
         LBJ = LB(J)
C
C----------------------------------------- convention for Green function
         IF ( GF_CONV_RH ) THEN
C                                                basis functions R and H
            SPHFUNL = CI*P*(CJLZ(L,ZTOP)+CI*CNLZ(L,ZTOP))
            SPHFUNLP1 = CI*P*(CJLZ(L+1,ZTOP)+CI*CNLZ(L+1,ZTOP))
            SPHFUNLB = CI*P*(CJLZ(LBJ,ZTOP)+CI*CNLZ(LBJ,ZTOP))
            SPHFUNLBP1 = CI*P*(CJLZ(LBJ+1,ZTOP)+CI*CNLZ(LBJ+1,ZTOP))
C
         ELSE
C                                                basis functions Z and J
            SPHFUNL = CJLZ(L,ZTOP)
            SPHFUNLP1 = CJLZ(L+1,ZTOP)
            SPHFUNLB = CJLZ(LBJ,ZTOP)
            SPHFUNLBP1 = CJLZ(LBJ+1,ZTOP)
C
         END IF
C
         PI(J,J,IR) = SPHFUNL*R(IR)
         QI(J,J,IR) = CFAC*SK(J)*SPHFUNLB*R(IR)*C
         DP(J,J,IR) = (DBLE(L+1)*SPHFUNL-ZTOP*SPHFUNLP1)*DRDI(IR)
         DQ(J,J,IR) = CFAC*SK(J)*(DBLE(LBJ+1)*SPHFUNLB-ZTOP*SPHFUNLBP1)
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
         CALL DIRBS_SOCRAD(X,Y,DY,DRDI,B,V,DBDR,DVDR,R,NRMAX,NCFMAX)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN) - 1,IRCUT(IPAN-1) + 1, - 1
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL DIRBS_SOCSTP(Y,DY,NY,X,HBS,B,V,DBDR,DVDR,R,DRDI,YM,YN,
     &                        YERR,YSAV,YSEQ,DYSAV,NRMAX,NCFMAX)
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
C     the minor component for the soc-manipulated wf is meaningless
C     =>  set it to zero
C
      QI(:,:,:) = C0
C
      END
C*==dirbs_socstp.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBS_SOCSTP(Y,DY,NY,X,HBS,B,V,DBDR,DVDR,R,DRDI,YM,YN,
     &                        YERR,YSAV,YSEQ,DYSAV,NRMAX,NCFMAX)
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
C*--DIRBS_SOCSTP566
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
      COMPLEX*16 B(NRMAX),DBDR(NRMAX),DVDR(NRMAX),DY(NY),DYSAV(NCFMAX),
     &           V(NRMAX),Y(NY),YERR(NCFMAX),YM(NCFMAX),YN(NCFMAX),
     &           YSAV(NCFMAX),YSEQ(NCFMAX)
      REAL*8 DRDI(NRMAX),R(NRMAX)
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
         CALL DIRBS_SOCRAD(XMID,YN,YSEQ,DRDI,B,V,DBDR,DVDR,R,NRMAX,
     &                     NCFMAX)
C
         H2MID = 2.D0*HMID
         DO ISTEP = 2,NSEQ(I)
            DO IY = 1,NY
               SWAP = YM(IY) + H2MID*YSEQ(IY)
               YM(IY) = YN(IY)
               YN(IY) = SWAP
            END DO
            XMID = XMID + HMID
            CALL DIRBS_SOCRAD(XMID,YN,YSEQ,DRDI,B,V,DBDR,DVDR,R,NRMAX,
     &                        NCFMAX)
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
         CALL DIRBS_SOCRAD(X,Y,DY,DRDI,B,V,DBDR,DVDR,R,NRMAX,NCFMAX)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<DIRBS_SOCSTP>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error : ',ERRMAX
      WRITE (6,*) 'tolerance TOL_DIRBS : ',TOL_DIRBS
      WRITE (6,*) 'grid position     X : ',X
      X = X + HBS
C
      CALL DIRBS_SOCRAD(X,Y,DY,DRDI,B,V,DBDR,DVDR,R,NRMAX,NCFMAX)
C
      END
C*==dirbs_socrad.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBS_SOCRAD(XBS,Y,DYDX,DRDI,B,V,DBDR,DVDR,R,NRMAX,
     &                        NCFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,CGD,CGMD,CGO,CGOZ,CGZ,LMK,KPX,
     &    NSOLBS,JLAG1,K_SOC
      IMPLICIT NONE
C*--DIRBS_SOCRAD732
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
      COMPLEX*16 B(NRMAX),DBDR(NRMAX),DVDR(NRMAX),DYDX(NCFMAX),V(NRMAX),
     &           Y(NCFMAX)
      REAL*8 DRDI(NRMAX),R(NRMAX)
C
C Local variables
C
      COMPLEX*16 BBS,BPP,BPP_CGO,BQQ,DBDRBS,DVDRBS,EMVPP,
     &           EMVPP_BPP_CGD(2),EMVQQ,EMVQQ_BQQ_CGMD(2),SOCPP(2),VBS
      REAL*8 DRDIBS,DROVR,DX,DXSQ,KPX_DROVR(2),RBS,W(NLAG)
      INTEGER I,J,K,M
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-12 ) THEN
         I = NINT(XBS)
         VBS = V(I)
         BBS = B(I)
         RBS = R(I)
         DRDIBS = DRDI(I)
         DVDRBS = DVDR(I)
         DBDRBS = DBDR(I)
      ELSE
         DX = XBS - JLAG1
         DXSQ = DX*DX
         W(3) = 0.5D0*(-DX+DXSQ)
         W(1) = 1D0 - DX + W(3)
         W(2) = DX + DX - DXSQ
         VBS = 0D0
         BBS = 0D0
         RBS = 0D0
         DVDRBS = 0D0
         DBDRBS = 0D0
         DRDIBS = 0D0
         J = JLAG1 - 1
         DO I = 1,NLAG
            J = J + 1
            VBS = VBS + W(I)*V(J)
            BBS = BBS + W(I)*B(J)
            DVDRBS = DVDRBS + W(I)*DVDR(J)
            DBDRBS = DBDRBS + W(I)*DBDR(J)
            RBS = RBS + W(I)*R(J)
            DRDIBS = DRDIBS + W(I)*DRDI(J)
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
         KPX_DROVR(I) = KPX(I)*DROVR
         EMVQQ_BQQ_CGMD(I) = EMVQQ + BQQ*CGMD(I)
         EMVPP_BPP_CGD(I) = EMVPP + BPP*CGD(I)
      END DO
      BPP_CGO = BPP*CGO
C
C-----------------------------------------------------------------------
      IF ( K_SOC.EQ.1 ) THEN
C
         DO I = 1,NSOLBS
            SOCPP(I) = LMK(I)*DROVR**2/(EMVQQ+BQQ*CGMD(I))
         END DO
C
         M = 0
         DO J = 1,NSOLBS
            DO I = 1,NSOLBS
               K = 3 - 4*(I-1)
C
CCCCC       DYDX(M+1) = -KPX(I)*Y(M+1)/RBS*DRDIBS + (EMVQQ+BQQ*CGMD(I))
CCCCC       DYDX(M+1) = -KPX(I)*Y(M+1)*DROVR + (EMVQQ+BQQ*CGMD(I))
CCCCC&                  *Y(M+2)
C
               DYDX(M+1) = -KPX_DROVR(I)*Y(M+1) + EMVQQ_BQQ_CGMD(I)
     &                     *Y(M+2)
C
CCCCC       DYDX(M+2) = KPX(I)*Y(M+2)/RBS*DRDIBS + (EMVPP+BPP*CGD(I))
CCCCC       DYDX(M+2) = KPX(I)*Y(M+2)*DROVR + (EMVPP+BPP*CGD(I))*Y(M+1)
CCCCC&                  + BPP*CGO*Y(M+K)
C
               DYDX(M+2) = KPX_DROVR(I)*Y(M+2) + EMVPP_BPP_CGD(I)*Y(M+1)
     &                     + BPP_CGO*Y(M+K) + SOCPP(I)*Y(M+1)
               M = M + 2
C
            END DO
         END DO
C
C-----------------------------------------------------------------------
      ELSE
C
         DO I = 1,NSOLBS
            SOCPP(I) = DROVR*CSQR*(DVDRBS-DBDRBS*CGMD(I))
     &                 /(EBS-VBS+CSQR+BBS*CGMD(I))**2
         END DO
C
         M = 0
         DO J = 1,NSOLBS
            DO I = 1,NSOLBS
               K = 3 - 4*(I-1)
C
CCCCC       DYDX(M+1) = -KPX(I)*Y(M+1)/RBS*DRDIBS + (EMVQQ+BQQ*CGMD(I))
CCCCC       DYDX(M+1) = -KPX(I)*Y(M+1)*DROVR + (EMVQQ+BQQ*CGMD(I))
CCCCC&                  *Y(M+2)
C
               DYDX(M+1) = -KPX_DROVR(I)*Y(M+1) + EMVQQ_BQQ_CGMD(I)
     &                     *Y(M+2)
C
CCCCC       DYDX(M+2) = KPX(I)*Y(M+2)/RBS*DRDIBS + (EMVPP+BPP*CGD(I))
CCCCC       DYDX(M+2) = KPX(I)*Y(M+2)*DROVR + (EMVPP+BPP*CGD(I))*Y(M+1)
CCCCC&                  + BPP*CGO*Y(M+K)
C
               DYDX(M+2) = KPX_DROVR(I)*Y(M+2) + EMVPP_BPP_CGD(I)*Y(M+1)
     &                     + (BPP_CGO+SOCPP(3-I)*CGOZ)*Y(M+K) + SOCPP(I)
     &                     *CGZ(I)*Y(M+1)
               M = M + 2
            END DO
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      END
