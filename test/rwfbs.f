C*==rwfbs.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RWFBS(GETIRRSOL,C,E,P,L,V,Z,R,DRDI,IRTOP,IRCUT,NPAN,PR,
     &                 QR,DPRIRTOP,PI,QI,NRMAX,GAMMA,IGAMMA,IREL,
     &                 GF_CONV_RH)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE RADIAL DIRAC EQUATIONS WITHIN THE         *
C   *   SCALAR RELATIVISTIC APPROXIMATION                              *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and the inward integration is started analytically             *
C   *   the integration itself is done by the BURLISCH-STOER method    *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point NTOP  PR,QR    *
C   *   and PI,QI with P=r*g and Q=r*c*f and R/I standing for the      *
C   *   regular/irregular solution                                     *
C   *                                                                  *
C   *   NOTE: the routine expects the potential term V to be COMPLEX   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,LIM,LBS,JLAG1
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--RWFBS27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLAG,NPEMAX,NCOEFF,NCFMAX
      PARAMETER (NLAG=3,NPEMAX=8,NCOEFF=4,NCFMAX=2)
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
      COMPLEX*16 CJLZ,CNLZ
      REAL*8 CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),HBS,TZ,X
      COMPLEX*16 CSCL,DPD0,DPN,DY(NCFMAX),EMV0,PC(-2:NCOEFF),PIP0,PJM2,
     &           PJM3,S,SPHFUNL,SPHFUNLP1,VC(0:NPEMAX),Y(NCFMAX),ZTOP
      INTEGER IPAN,IR,IR_LOWER_BOUND,J,NPE,NY
C
C*** End of declarations rewritten by SPAG
C
      CSQR = C*C
      CSCL = 1.0D5
      EBS = E
      LBS = L
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
      EMV0 = EBS - VC(0) + CSQR
      PJM3 = (EBS-VC(0))*EMV0**2/CSQR
      PJM2 = TZ*EMV0*(3.0D0*(EBS-VC(0))+CSQR)/CSQR
C
      DO J = -2,NCOEFF
         PC(J) = C0
      END DO
C
      PC(0) = 1.0D0/CSCL
C
C --------------------------------------- non-relativistic initilisation
C
      IF ( IREL.EQ.0 ) THEN
C
         IF ( Z.EQ.0 ) THEN
            PC(2) = -(EBS-VC(0))*PC(0)/(2.0D0*(2.0D0*GAMMA+1.0D0))
            PC(3) = VC(1)*PC(0)/(6.0D0*(GAMMA+1.0D0))
            PC(4) = (-(EBS-VC(0))*PC(2)+VC(2)*PC(0))
     &              /(4.0D0*(2.0D0*GAMMA+3.0D0))
C
         ELSE
C
            PC(1) = -TZ/(2.0D0*GAMMA)*PC(0)
            PC(2) = (-TZ*PC(1)-(EBS-VC(0))*PC(0))
     &              /(2.0D0*(2.0D0*GAMMA+1.0D0))
            PC(3) = (-TZ*PC(2)-(EBS-VC(0))*PC(1)+VC(1)*PC(0))
     &              /(6.0D0*(GAMMA+1.0D0))
            PC(4) = (-TZ*PC(3)-(EBS-VC(0))*PC(2)+VC(1)*PC(1)+VC(2)*PC(0)
     &              )/(4.0D0*(2.0D0*GAMMA+3.0D0))
         END IF
C
C ------------------------------------ scalar relativistic initilisation
C
      ELSE IF ( Z.EQ.0 ) THEN
         PC(1) = ((-VC(1)*L)/((EBS-VC(0)+CSQR)*2.0D0*(L+1.0D0)))*PC(0)
         PC(2) = (-PC(0)*(((EBS-VC(0))*(EBS-VC(0)+CSQR)**2)/CSQR+VC(2)*L
     &           )+VC(1)*(L+1.0D0)*PC(1))
     &           /(2.0D0*(EBS-VC(0)+CSQR)*(2.0D0*L+3.0D0))
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
C ---------------------------- power series to initialise Bulirsch Stoer
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
      Y(1) = PR(1)
      Y(2) = QR(1)
C
      X = 1.0D0
      HBS = 1.0D0
      NY = 2
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = 1,NPAN
C
         CALL RWFBSRAD(X,Y,DY,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C
         IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN-1) + 2,MIN(IRCUT(IPAN),IRTOP)
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR-(NLAG+1)/2+1)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL RWFBSSTP(Y,DY,NY,X,HBS,V,R,DRDI,NRMAX,NCFMAX,GAMMA)
C
            PR(IR) = Y(1)
            QR(IR) = Y(2)
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
            PR(IR) = Y(1)
            QR(IR) = Y(2)
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IR = 1,IRTOP
         PR(IR) = PR(IR)*CSCL
         QR(IR) = QR(IR)*CSCL
      END DO
C
      DPRIRTOP = DY(1)*CSCL
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C ######################################################################
C                     irregular solution
C ######################################################################
C
C         calculate the initial values of the wavefunction
C                     at the sphere boundary
C
      IR = IRTOP
C
      ZTOP = P*R(IR)
CXXXXXS = CSQR/(EBS*LIM+CSQR)
      S = CSQR/((EBS-V(IR))*LIM+CSQR)
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
C ======================================================================
C
      Y(1) = PI(IRTOP)/CSCL
      Y(2) = QI(IRTOP)/CSCL
C
      X = DBLE(IRTOP)
      HBS = -1.0D0
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = NPAN,1, - 1
C
         IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
C
         CALL RWFBSRAD(X,Y,DY,DRDI,V,R,NRMAX,NCFMAX,IGAMMA)
C
C         if( ipan .eq. 1 ) then
C             write(*,*) '********************* L=',l
C             write(*,*) '*************** DY(1)=',DY(1)
C             dy(1) = DPN*drdi(irtop)/cscl
C             write(*,*) '*************** DY(1)=',DY(1)
C         end if
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN) - 1,IRCUT(IPAN-1) + 1, - 1
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL RWFBSSTP(Y,DY,NY,X,HBS,V,R,DRDI,NRMAX,NCFMAX,IGAMMA)
C
            PI(IR) = Y(1)
            QI(IR) = Y(2)
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
            PI(IR) = Y(1)
            QI(IR) = Y(2)
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IR = 1,IRTOP - 1
         PI(IR) = PI(IR)*CSCL
         QI(IR) = QI(IR)*CSCL
      END DO
C
      END
C*==rwfbsstp.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RWFBSSTP(Y,DYDX,NY,X,HBS,V,R,DRDI,NRMAX,NCFMAX,GAMMA)
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DXDY  for last mesh-point                        *
C   *   on exit:  X,Y,DXDY  updated for X = X(last) + HBS              *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *         no step size adjusted in case of no convergency > STOP   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:TOL_RWFBS
      IMPLICIT NONE
C*--RWFBSSTP323
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ISEQMAX,NUSE
      PARAMETER (ISEQMAX=35,NUSE=7)
      REAL*8 TINYBIT
      PARAMETER (TINYBIT=1.0D-20)
C
C Dummy arguments
C
      REAL*8 GAMMA,HBS,X
      INTEGER NCFMAX,NRMAX,NY
      REAL*8 DRDI(NRMAX),R(NRMAX)
      COMPLEX*16 DYDX(NY),V(NRMAX),Y(NY)
C
C Local variables
C
      REAL*8 ABSBB,ERRMAX,FF(NUSE),H2MID,HMID,XEST,XMID,XX(ISEQMAX)
      COMPLEX*16 BB,BB1,CC,CRAT,DD(NY,NUSE),DYSAV(NCFMAX),DYY,SWAP,VV,
     &           YERR(NCFMAX),YM(NY),YN(NY),YSAV(NCFMAX),YSEQ(NCFMAX),YY
      INTEGER I,ISTEP,IV,J,JJ,KK,MM,NSEQ(ISEQMAX)
C
C*** End of declarations rewritten by SPAG
C
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,
     &     1024,1536,2048,3072,4096,6144,8192,12288,16384,24576,32768,
     &     49152,65536,98304,131072,262144,524288,1048576/
C
      DO I = 1,NY
         YSAV(I) = Y(I)
         DYSAV(I) = DYDX(I)
      END DO
      DO I = 1,ISEQMAX
C
C------------------------------------------------------------------- MID
C
         HMID = HBS/NSEQ(I)
         DO IV = 1,NY
            YM(IV) = YSAV(IV)
            YN(IV) = YSAV(IV) + HMID*DYSAV(IV)
         END DO
         XMID = X + HMID
C
         CALL RWFBSRAD(XMID,YN,YSEQ,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C
         H2MID = HMID + HMID
C
         DO ISTEP = 2,NSEQ(I)
            DO IV = 1,NY
               SWAP = YM(IV) + H2MID*YSEQ(IV)
               YM(IV) = YN(IV)
               YN(IV) = SWAP
            END DO
            XMID = XMID + HMID
C
            CALL RWFBSRAD(XMID,YN,YSEQ,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C
         END DO
C
         DO IV = 1,NY
            YSEQ(IV) = 0.5D0*(YM(IV)+YN(IV)+HMID*YSEQ(IV))
         END DO
C
C-----------------------------------------------------------------------
C
         XEST = HBS/NSEQ(I)
         XEST = XEST*XEST
C
C------------------------------------------------------------------- RZE
C
         XX(I) = XEST
C
         IF ( I.EQ.1 ) THEN
            DO JJ = 1,NY
               Y(JJ) = YSEQ(JJ)
               DD(JJ,1) = YSEQ(JJ)
               YERR(JJ) = YSEQ(JJ)
            END DO
         ELSE
            MM = MIN(I,NUSE)
            DO KK = 1,MM - 1
               FF(KK+1) = XX(I-KK)/XEST
            END DO
            DO JJ = 1,NY
               YY = YSEQ(JJ)
               VV = DD(JJ,1)
               CC = YY
               DD(JJ,1) = YY
               DO KK = 2,MM
                  BB1 = FF(KK)*VV
                  BB = BB1 - CC
                  ABSBB = ABS(DREAL(BB)) + ABS(DIMAG(BB))
                  IF ( ABSBB.GT.1D-40 ) THEN
                     BB = (CC-VV)/BB
                     DYY = CC*BB
                     CC = BB1*BB
C
                  ELSE
C
                     DYY = VV
                  END IF
C
                  VV = DD(JJ,KK)
C
                  DD(JJ,KK) = DYY
                  YY = YY + DYY
               END DO
               YERR(JJ) = DYY
               Y(JJ) = YY
            END DO
C
         END IF
C
C-----------------------------------------------------------------------
C
         DO J = 1,NY
            CRAT = YERR(J)/(Y(J)+TINYBIT)
            ERRMAX = ABS(DREAL(CRAT)) + ABS(DIMAG(CRAT))
            IF ( ERRMAX.GT.TOL_RWFBS ) GOTO 100
         END DO
C
C------------------------------------------------------------- converged
C
         X = X + HBS
C
         CALL RWFBSRAD(X,Y,DYDX,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<RWFBSSTP>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error : ',ERRMAX
      WRITE (6,*) 'tolerance TOL_RWFBS : ',TOL_RWFBS
      WRITE (6,*) 'grid position     X : ',X
C
      X = X + HBS
C
      CALL RWFBSRAD(X,Y,DYDX,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C
C
      END
C*==rwfbsrad.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RWFBSRAD(XBS,Y,DYDX,DRDI,V,R,NRMAX,NCFMAX,GAMMA)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,LIM,LBS,JLAG1
      IMPLICIT NONE
C*--RWFBSRAD497
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
      REAL*8 GAMMA,XBS
      INTEGER NCFMAX,NRMAX
      REAL*8 DRDI(NRMAX),R(NRMAX)
      COMPLEX*16 DYDX(NCFMAX),V(NRMAX),Y(NCFMAX)
C
C Local variables
C
      REAL*8 DRDIBS,DROVR,DX,DXSQ,RBS,W(NLAG)
      COMPLEX*16 EMVPP,EMVQQ,S,VBS
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-14 ) THEN
         I = NINT(XBS)
         VBS = V(I)
         RBS = R(I)
         DRDIBS = DRDI(I)
      ELSE
         DX = XBS - JLAG1
         DXSQ = DX*DX
         W(3) = 0.5D0*(-DX+DXSQ)
         W(1) = 1D0 - DX + W(3)
         W(2) = DX + DX - DXSQ
         VBS = 0D0
         RBS = 0D0
         DRDIBS = 0D0
         J = JLAG1 - 1
         DO I = 1,NLAG
            J = J + 1
            VBS = VBS + W(I)*V(J)
            RBS = RBS + W(I)*R(J)
            DRDIBS = DRDIBS + W(I)*DRDI(J)
         END DO
      END IF
C
C-----------------------------------------------------------------------
C
      S = CSQR/((EBS-VBS)*LIM+CSQR)
      EMVPP = -(EBS-VBS-S*DBLE(LBS*(LBS+1.0D0))/RBS**2)*DRDIBS
      EMVQQ = DRDIBS/S
      DROVR = DRDIBS/RBS
C
      DYDX(1) = -(GAMMA-1.0D0)*Y(1)*DROVR + EMVQQ*Y(2)
      DYDX(2) = -(GAMMA+1.0D0)*Y(2)*DROVR + EMVPP*Y(1)
C
      END
