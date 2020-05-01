C*==rwfrk.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RWFRK(GETIRRSOL,C,E,P,L,V,Z,R,DRDI,IRTOP,IRCUT,NPAN,PR,
     &                 QR,DPRIRTOP,PI,QI,NRMAX,GAMMA,IGAMMA,IREL,
     &                 GF_CONV_RH)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE RADIAL DIRAC EQUATIONS WITHIN THE         *
C   *   SCALAR RELATIVISTIC APPROXIMATION                              *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and the inward integration is started analytically             *
C   *   the integration itself is done by the RUNGE-KUTTA method       *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point NTOP  PR,QR    *
C   *   and PI,QI with P=r*g and Q=r*c*f and R/I standing for the      *
C   *   regular/irregular solution                                     *
C   *                                                                  *
C   *   NOTE: the routine expects the potential term V to be COMPLEX   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--RWFRK25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NPEMAX,NCOEFF
      PARAMETER (NPEMAX=8,NCOEFF=4)
C
C Dummy arguments
C
      REAL*8 C,GAMMA,IGAMMA
      COMPLEX*16 DPRIRTOP,E,P
      LOGICAL GETIRRSOL,GF_CONV_RH
      INTEGER IREL,IRTOP,L,NPAN,NRMAX,Z
      REAL*8 DRDI(NRMAX),R(NRMAX)
      INTEGER IRCUT(0:NPAN)
      COMPLEX*16 PI(NRMAX),PR(NRMAX),QI(NRMAX),QR(NRMAX),V(NRMAX)
C
C Local variables
C
      COMPLEX*16 CINTPOL1,CINTPOL2,CINTPOL3
      COMPLEX*16 CJLZ,CNLZ
      REAL*8 CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),DRDIRK,DROVR,LIM,
     &       QQOVPP,RRK,RY1,RY2,RY3,RY4,TZ,XPP(:),XPPMID(:),XQQ(:),
     &       XQQMID(:)
      COMPLEX*16 CSQR,CY1,CY2,CY3,CY4,DPD0,DPN,EMV0,EMVPP,EMVQQ,
     &           PC(-2:NCOEFF),PIP0,PJM2,PJM3,S,SPHFUNL,SPHFUNLP1,
     &           VC(0:NPEMAX),VRK,XPQ(:),XPQMID(:),XQP(:),XQPMID(:),ZTOP
      INTEGER IPAN,IR,IRLIM1,IRLIM2,J,NPE
      REAL*8 RINTPOL1,RINTPOL2,RINTPOL3
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
      ALLOCATE (XPP(NRMAX),XPPMID(NRMAX),XQQ(NRMAX),XQQMID(NRMAX))
      ALLOCATE (XPQ(NRMAX),XPQMID(NRMAX),XQP(NRMAX),XQPMID(NRMAX))
C
      CSQR = C*C
      TZ = 2.0D0*Z
C
      IF ( IREL.EQ.0 ) THEN
         LIM = 0.0D0
      ELSE
         LIM = 1.0D0
      END IF
C
      NPE = 3
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      DO IR = 1,IRTOP
C
         VRK = V(IR)
         RRK = R(IR)
         DRDIRK = DRDI(IR)
C
         IF ( IREL.EQ.0 ) THEN
            S = 1D0
         ELSE
            S = CSQR/((E-VRK)+CSQR)
         END IF
         EMVPP = -(E-VRK-S*DBLE(L*(L+1))/(RRK*RRK))*DRDIRK
         EMVQQ = DRDIRK/S
         DROVR = DRDIRK/RRK
C
         XPP(IR) = -(GAMMA-1.0D0)*DROVR
         XPQ(IR) = EMVQQ
         XQP(IR) = EMVPP
         XQQ(IR) = -(GAMMA+1.0D0)*DROVR
C
      END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
      QQOVPP = (GAMMA+1.0D0)/(GAMMA-1.0D0)
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
               XPPMID(IR) = RINTPOL1(XPP(IR),XPP(IR+1),XPP(IR+2))
               XPQMID(IR) = CINTPOL1(XPQ(IR),XPQ(IR+1),XPQ(IR+2))
               XQPMID(IR) = CINTPOL1(XQP(IR),XQP(IR+1),XQP(IR+2))
C              XQQMID(IR) = RINTPOL1(XQQ(IR),XQQ(IR+1),XQQ(IR+2))
C
            ELSE IF ( IR.LT.IRLIM2 ) THEN
C
               XPPMID(IR) = RINTPOL2(XPP(IR-1),XPP(IR),XPP(IR+1),
     &                      XPP(IR+2))
               XPQMID(IR) = CINTPOL2(XPQ(IR-1),XPQ(IR),XPQ(IR+1),
     &                      XPQ(IR+2))
               XQPMID(IR) = CINTPOL2(XQP(IR-1),XQP(IR),XQP(IR+1),
     &                      XQP(IR+2))
C              XQQMID(IR) = RINTPOL2(XQQ(IR-1),XQQ(IR),XQQ(IR+1),
C    &                      XQQ(IR+2))
C
            ELSE
C
               XPPMID(IR) = RINTPOL3(XPP(IR-1),XPP(IR),XPP(IR+1))
               XPQMID(IR) = CINTPOL3(XPQ(IR-1),XPQ(IR),XPQ(IR+1))
               XQPMID(IR) = CINTPOL3(XQP(IR-1),XQP(IR),XQP(IR+1))
C              XQQMID(IR) = RINTPOL3(XQQ(IR-1),XQQ(IR),XQQ(IR+1))
C
            END IF
C
            XQQMID(IR) = QQOVPP*XPPMID(IR)
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
         IR = MIN(IRCUT(IPAN),IRTOP)
         XPPMID(IR) = 0D0
         XPQMID(IR) = 0D0
         XQPMID(IR) = 0D0
         XQQMID(IR) = 0D0
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C ---------------------------- calculation of the potential coefficients
C
      DO J = 1,NPE
         DO IR = 1,NPE
            CM(IR,J) = DBLE(R(IR)**(J-1))
         END DO
      END DO
C
      CALL RINVGJ(CMI,CM,NPEMAX,NPE)
C
      DO J = 1,NPE
         VC(J-1) = 0.0D0
         DO IR = 1,NPE
            IF ( (.NOT.FINITE_NUCLEUS) .OR. (Z.EQ.0) ) THEN
               VC(J-1) = VC(J-1) + CMI(J,IR)*(V(IR)+TZ/R(IR))
            ELSE
               VC(J-1) = VC(J-1) + CMI(J,IR)*V(IR)
            END IF
         END DO
      END DO
C
      EMV0 = E - VC(0) + CSQR
      PJM3 = (E-VC(0))*EMV0**2/CSQR
      PJM2 = TZ*EMV0*(3.0D0*(E-VC(0))+CSQR)/CSQR
C
      DO J = -2,NCOEFF
         PC(J) = C0
      END DO
C
      PC(0) = 1.0D0
C
C -------------------------------------- non-relativistic initialisation
C
      IF ( IREL.EQ.0 ) THEN
C
         IF ( Z.EQ.0 ) THEN
            PC(2) = -(E-VC(0))*PC(0)/(2.0D0*(2.0D0*GAMMA+1.0D0))
            PC(3) = VC(1)*PC(0)/(6.0D0*(GAMMA+1.0D0))
            PC(4) = (-(E-VC(0))*PC(2)+VC(2)*PC(0))
     &              /(4.0D0*(2.0D0*GAMMA+3.0D0))
C
         ELSE
C
            PC(1) = -TZ/(2.0D0*GAMMA)*PC(0)
            PC(2) = (-TZ*PC(1)-(E-VC(0))*PC(0))
     &              /(2.0D0*(2.0D0*GAMMA+1.0D0))
            PC(3) = (-TZ*PC(2)-(E-VC(0))*PC(1)+VC(1)*PC(0))
     &              /(6.0D0*(GAMMA+1.0D0))
            PC(4) = (-TZ*PC(3)-(E-VC(0))*PC(2)+VC(1)*PC(1)+VC(2)*PC(0))
     &              /(4.0D0*(2.0D0*GAMMA+3.0D0))
         END IF
C
C ----------------------------------- scalar relativistic initialisation
C
      ELSE IF ( Z.EQ.0 ) THEN
         PC(1) = ((-VC(1)*L)/((E-VC(0)+CSQR)*2.0D0*(L+1.0D0)))*PC(0)
         PC(2) = (-PC(0)*(((E-VC(0))*(E-VC(0)+CSQR)**2)/CSQR+VC(2)*L)
     &           +VC(1)*(L+1.0D0)*PC(1))
     &           /(2.0D0*(E-VC(0)+CSQR)*(2.0D0*L+3.0D0))
      ELSE
C
         DO J = 1,NCOEFF
            PC(J) = ((EMV0*((J+GAMMA-1.0D0)*(J+GAMMA-2.0D0)-L*(L+1.0D0)+
     &              3.0D0*TZ**2/CSQR)-TZ**2)*PC(J-1)+PJM3*PC(J-3)
     &              +PJM2*PC(J-2))/(-1.0D0*TZ*J*(2.0D0*GAMMA+J))
         END DO
C
      END IF
C
C ------------------------------- power series to initialise Runge-Kutta
C
      PIP0 = PC(0)
      DPD0 = 0.0D0
C
      DO J = 1,NCOEFF
         PIP0 = PIP0 + PC(J)*R(1)**J
         DPD0 = DPD0 + DBLE(J)*PC(J)*R(1)**(J-1)
      END DO
C
      PR(1) = PIP0
      QR(1) = (DPD0+(GAMMA-1.0D0)*PIP0/R(1))*CSQR/((E-V(1))*LIM+CSQR)
C
C ======================================================================
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
         ZTOP = P*R(IR)
C
C        S = CSQR/(E*LIM+CSQR)
         S = CSQR/((E-V(IR))*LIM+CSQR)
C
C----------------------------------------- convention for Green function
         IF ( GF_CONV_RH ) THEN
C                                                basis functions R and H
            SPHFUNL = CI*P*(CJLZ(L,ZTOP)+CI*CNLZ(L,ZTOP))
            SPHFUNLP1 = CI*P*(CJLZ(L+1,ZTOP)+CI*CNLZ(L+1,ZTOP))
C
         ELSE
C                                                basis functions Z and J
            SPHFUNL = CJLZ(L,ZTOP)
            SPHFUNLP1 = CJLZ(L+1,ZTOP)
         END IF
C
         PI(IR) = SPHFUNL*R(IR)**(-IGAMMA+1.0D0)
         DPN = (L+1)*SPHFUNL - ZTOP*SPHFUNLP1
         DPN = DPN*R(IR)**(-IGAMMA) - IGAMMA/R(IR)*PI(IR)
         QI(IR) = (DPN*R(IR)-(-IGAMMA+1.0D0)*PI(IR))*S/R(IR)
C
      END IF
C
C ======================================================================
C
      CALL SOLVE_RK2(GETIRRSOL,IRTOP,IRCUT,NPAN,PR,QR,PI,QI,XPP,XPPMID,
     &               XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID,1,NRMAX)
C
C ======================================================================
C
C
C-------------------------- "radial" derivative dP/di at last mesh point
C
      IR = IRTOP
      DPRIRTOP = XPP(IR)*PR(IR) + XPQ(IR)*QR(IR)
C
      END
C*==solve_rk2.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SOLVE_RK2(GETIRRSOL,IRTOP,IRCUT,NPAN,PR,QR,PI,QI,XPP,
     &                     XPPMID,XQQ,XQQMID,XPQ,XPQMID,XQP,XQPMID,NDIM,
     &                     NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine to solve radial differential equations        *
C   *   using the  RUNGE-KUTTA method                                  *
C   *                                                                  *
C   *   solve 2 coupled equations for P and Q                          *
C   *                                                                  *
C   *   NOTE: the wave functions at the first mesh points have to be   *
C   *         supplied,  i.e.  PR(1), QR(1)  and  PR(IRTOP), QR(IRTOP) *
C   *                                                                  *
C   *         the wave functions have additional indices to allow the  *
C   *         routine to be used for the NON- as well as FULL          *
C   *         relativistic case  NDIM= 1 or 2, respectively            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SOLVE_RK2334
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL GETIRRSOL
      INTEGER IRTOP,NDIM,NPAN,NRMAX
      INTEGER IRCUT(0:NPAN)
      COMPLEX*16 PI(NDIM,NDIM,NRMAX),PR(NDIM,NDIM,NRMAX),
     &           QI(NDIM,NDIM,NRMAX),QR(NDIM,NDIM,NRMAX),XPQ(NRMAX),
     &           XPQMID(NRMAX),XQP(NRMAX),XQPMID(NRMAX)
      REAL*8 XPP(NRMAX),XPPMID(NRMAX),XQQ(NRMAX),XQQMID(NRMAX)
C
C Local variables
C
      INTEGER IPAN,IR
      COMPLEX*16 K1P,K1Q,K2P,K2Q,K3P,K3Q,K4P,K4Q,PAUX,QAUX
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
            K1P = XPP(IR)*PR(1,1,IR) + XPQ(IR)*QR(1,1,IR)
            K1Q = XQP(IR)*PR(1,1,IR) + XQQ(IR)*QR(1,1,IR)
C
            PAUX = PR(1,1,IR) + 0.5D0*K1P
            QAUX = QR(1,1,IR) + 0.5D0*K1Q
C
            K2P = XPPMID(IR)*PAUX + XPQMID(IR)*QAUX
            K2Q = XQPMID(IR)*PAUX + XQQMID(IR)*QAUX
C
            PAUX = PR(1,1,IR) + 0.5D0*K2P
            QAUX = QR(1,1,IR) + 0.5D0*K2Q
C
            K3P = XPPMID(IR)*PAUX + XPQMID(IR)*QAUX
            K3Q = XQPMID(IR)*PAUX + XQQMID(IR)*QAUX
C
            PAUX = PR(1,1,IR) + K3P
            QAUX = QR(1,1,IR) + K3Q
C
            K4P = XPP(IR+1)*PAUX + XPQ(IR+1)*QAUX
            K4Q = XQP(IR+1)*PAUX + XQQ(IR+1)*QAUX
C
            PR(1,1,IR+1) = PR(1,1,IR) + (K1P+K2P+K2P+K3P+K3P+K4P)/6D0
            QR(1,1,IR+1) = QR(1,1,IR) + (K1Q+K2Q+K2Q+K3Q+K3Q+K4Q)/6D0
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
            PR(1,1,IR) = PR(1,1,IR-1)
            QR(1,1,IR) = QR(1,1,IR-1)
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
            K1P = XPP(IR)*PI(1,1,IR) + XPQ(IR)*QI(1,1,IR)
            K1Q = XQP(IR)*PI(1,1,IR) + XQQ(IR)*QI(1,1,IR)
C
            PAUX = PI(1,1,IR) - 0.5D0*K1P
            QAUX = QI(1,1,IR) - 0.5D0*K1Q
C
            K2P = XPPMID(IR-1)*PAUX + XPQMID(IR-1)*QAUX
            K2Q = XQPMID(IR-1)*PAUX + XQQMID(IR-1)*QAUX
C
            PAUX = PI(1,1,IR) - 0.5D0*K2P
            QAUX = QI(1,1,IR) - 0.5D0*K2Q
C
            K3P = XPPMID(IR-1)*PAUX + XPQMID(IR-1)*QAUX
            K3Q = XQPMID(IR-1)*PAUX + XQQMID(IR-1)*QAUX
C
            PAUX = PI(1,1,IR) - K3P
            QAUX = QI(1,1,IR) - K3Q
C
            K4P = XPP(IR-1)*PAUX + XPQ(IR-1)*QAUX
            K4Q = XQP(IR-1)*PAUX + XQQ(IR-1)*QAUX
C
            PI(1,1,IR-1) = PI(1,1,IR) - (K1P+K2P+K2P+K3P+K3P+K4P)/6D0
            QI(1,1,IR-1) = QI(1,1,IR) - (K1Q+K2Q+K2Q+K3Q+K3Q+K4Q)/6D0
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------ to first mesh point of next panel (except for last panel)
C
         IF ( IPAN.NE.1 ) THEN
            IR = IRCUT(IPAN-1)
            PI(1,1,IR) = PI(1,1,IR+1)
            QI(1,1,IR) = QI(1,1,IR+1)
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      END
