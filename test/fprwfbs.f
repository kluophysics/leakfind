C*==fprwfbs.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPRWFBS(IREL,SPNWGT,GETIRRSOL,C,E,P,CV0,VNS,BNS,LDAU,W,
     &                   NFPT,JRNS1,R,DRDI,VAMEG,IBLK,NSOLBLK,IKMSOLBLK,
     &                   IRCUT,NPAN,PRM,QRM,PIM,QIM,PHL,QHL,GAMMA,
     &                   IGAMMA,RGAM,RIGAM,DPRIRTOP,GF_CONV_RH)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE RADIAL DIRAC EQUATIONS WITHIN THE         *
C   *   SCALAR RELATIVISTIC APPROXIMATION                              *
C   *                                                                  *
C   *   NOTE: the routine expects the potential terms to be COMPLEX    *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,LIM,JLAG1,KNS
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_RMESH,ONLY:NPANMAX,NRMAX,JRNSMIN
      USE MOD_ANGMOM,ONLY:L_LM,NKMMAX,NL,NLMAX
      USE MOD_TYPES,ONLY:NCPLWFMAX,NLMFPMAX
      IMPLICIT NONE
C*--FPRWFBS23
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
      REAL*8 C
      COMPLEX*16 E,P,SPNWGT
      LOGICAL GETIRRSOL,GF_CONV_RH,LDAU
      INTEGER IBLK,IREL,JRNS1,NFPT,NPAN
      COMPLEX*16 BNS(JRNSMIN:NRMAX,NLMFPMAX),CV0(NRMAX),
     &           DPRIRTOP(NCPLWFMAX,NCPLWFMAX),PHL(NRMAX,NL),
     &           PIM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           PRM(NCPLWFMAX,NCPLWFMAX,NRMAX),QHL(NRMAX,NL),
     &           QIM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           QRM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VNS(JRNSMIN:NRMAX,NLMFPMAX),W(NCPLWFMAX,NCPLWFMAX)
      REAL*8 DRDI(NRMAX),GAMMA(0:NLMAX),IGAMMA(0:NLMAX),R(NRMAX),
     &       RGAM(1:NRMAX,0:NLMAX),RIGAM(1:NRMAX,0:NLMAX)
      INTEGER IKMSOLBLK(NKMMAX,NKMMAX),IRCUT(0:NPANMAX),NSOLBLK(NKMMAX)
C
C Local variables
C
      COMPLEX*16 CJLZ,CNLZ
      COMPLEX*16 CSUM,DPN,DY(:),S,SPHFUNL,SPHFUNLP1,VLL(:,:,:),Y(:),ZTOP
      REAL*8 GAMSOL(:),HBS,IGAMSOL(:),X
      INTEGER I,IA_ERR,IL,IPAN,IPOT,IR,IRBOT_PAN,IRTOP,IR_LOWER_BOUND,
     &        IY,J,L,LJ,LM,L_SOL(:),NSOL,NY,NYH
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE VLL,IGAMSOL,GAMSOL,L_SOL,DY,Y
C
      EBS = E
      CSQR = C*C
      IF ( IREL.EQ.0 ) THEN
         LIM = 0.0D0
      ELSE
         LIM = 1.0D0
      END IF
C
      NSOL = NSOLBLK(IBLK)
      NYH = NSOL*NSOL
      NY = NYH + NYH
      IRTOP = IRCUT(NPAN)
C
      ALLOCATE (IGAMSOL(NSOL),GAMSOL(NSOL),L_SOL(NSOL))
      ALLOCATE (VLL(NSOL,NSOL,JRNSMIN:NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPRWFBS -> VLL'
      ALLOCATE (DY(NKMMAX*NKMMAX),Y(NKMMAX*NKMMAX))
C
C=======================================================================
C                   set up potential
C=======================================================================
C     spherical part CV0 is already done in  <FPNRSSITE>
C
      DO IR = JRNS1,IRTOP
         DO I = 1,NSOL
            DO J = 1,NSOL
               CSUM = 0D0
               DO IPOT = 2,NFPT
                  CSUM = CSUM + VAMEG(I,J,IPOT,1)
     &                   *(VNS(IR,IPOT)+SPNWGT*BNS(IR,IPOT))
               END DO
               VLL(I,J,IR) = CSUM
            END DO
         END DO
      END DO
C
      IF ( LDAU ) THEN
         DO IR = JRNS1,IRTOP
            DO J = 1,NSOL
               DO I = 1,NSOL
                  VLL(I,J,IR) = VLL(I,J,IR) + W(J,I)
               END DO
            END DO
         END DO
      END IF
C
C=======================================================================
C      calculate regular radial wave functions for spherical regime
C=======================================================================
C
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*JRNS1,PRM)
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*JRNS1,QRM)
C
      DO I = 1,NSOLBLK(IBLK)
         LM = IKMSOLBLK(I,IBLK)
C
         L = L_LM(LM)
         GAMSOL(I) = GAMMA(L)
         IGAMSOL(I) = IGAMMA(L)
         L_SOL(I) = L
         IL = L + 1
C
         DO IR = 1,JRNS1
            PRM(I,I,IR) = PHL(IR,IL)*RGAM(IR,L)
            QRM(I,I,IR) = QHL(IR,IL)*RGAM(IR,L)
         END DO
C
      END DO
C
C=======================================================================
C       calculate radial wave functions for NON-spherical regime
C=======================================================================
C
      IY = 1
      DO J = 1,NSOL
         LJ = L_SOL(J)
         DO I = 1,NSOL
            Y(IY) = PRM(I,J,JRNS1)/RGAM(JRNS1,LJ)
            Y(IY+NYH) = QRM(I,J,JRNS1)/RGAM(JRNS1,LJ)
            IY = IY + 1
         END DO
      END DO
C
      X = DBLE(JRNS1)
      HBS = 1.0D0
      KNS = .TRUE.
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = 1,NPAN
C
         IF ( IPAN.GT.1 ) THEN
            IRBOT_PAN = IRCUT(IPAN-1) + 2
         ELSE
            IRBOT_PAN = JRNS1 + 1
         END IF
C
         CALL FPRWFBSRAD(GAMSOL,L_SOL,NSOL,X,Y,DY,NY,NYH,CV0,VLL,R,DRDI,
     &                   JRNSMIN,NRMAX)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRBOT_PAN,IRCUT(IPAN)
C
            JLAG1 = MAX(1,IR-(NLAG+1)/2+1)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL FPRWFBSSTP(GAMSOL,L_SOL,NSOL,X,HBS,Y,DY,NY,NYH,CV0,VLL,
     &                      R,DRDI,JRNSMIN,NRMAX)
C
            IY = 1
            DO J = 1,NSOL
               LJ = L_SOL(J)
               DO I = 1,NSOL
                  PRM(I,J,IR) = RGAM(IR,LJ)*Y(IY)
                  QRM(I,J,IR) = RGAM(IR,LJ)*Y(IY+NYH)
                  IY = IY + 1
               END DO
            END DO
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------------------------------ to first mesh point of next panel
C ---------------------- however, derivatives have to be calculated anew
C ----------------------- because DRDI jumps when crossing panel borders
C
         IF ( IPAN.NE.NPAN ) THEN
            IR = IRCUT(IPAN) + 1
            IY = 1
            DO J = 1,NSOL
               LJ = L_SOL(J)
               DO I = 1,NSOL
                  PRM(I,J,IR) = RGAM(IR,LJ)*Y(IY)
                  QRM(I,J,IR) = RGAM(IR,LJ)*Y(IY+NYH)
                  IY = IY + 1
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C-----------------------------------------------------------------------
C  store derivative of regular wave function at critical radius dP/di
C  NOTE the derivative DY refers to  \bar P = r^gamma + P
C-----------------------------------------------------------------------
C
      IY = 1
      DO J = 1,NSOL
         LJ = L_SOL(J)
         DO I = 1,NSOL
            DPRIRTOP(I,J) = PRM(I,J,IRTOP)*DRDI(IRTOP)*GAMMA(LJ)
     &                      /R(IRTOP) + DY(IY)*RGAM(IRTOP,LJ)
            IY = IY + 1
         END DO
      END DO
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C ######################################################################
C                     IRREGULAR SOLUTION
C ######################################################################
C
C         calculate the initial values of the wavefunction
C                     at the sphere boundary
C
      IR = IRTOP
C
      DO J = 1,NSOL
         DO I = 1,NSOL
            PIM(I,J,IR) = C0
            QIM(I,J,IR) = C0
         END DO
      END DO
C
      DO J = 1,NSOL
C
         L = L_SOL(J)
         ZTOP = P*R(IR)
         S = CSQR/(EBS*LIM+CSQR)
         S = CSQR/((EBS-CV0(IR))*LIM+CSQR)
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
         PIM(J,J,IR) = SPHFUNL*R(IR)**(-IGAMSOL(J)+1.0D0)
         DPN = (L+1)*SPHFUNL - ZTOP*SPHFUNLP1
         DPN = DPN*R(IR)**(-IGAMSOL(J)) - IGAMSOL(J)/R(IR)*PIM(J,J,IR)
         QIM(J,J,IR) = (DPN*R(IR)-(-IGAMSOL(J)+1.0D0)*PIM(J,J,IR))
     &                 *S/R(IR)
C
      END DO
C
C ======================================================================
C
      IY = 1
      DO J = 1,NSOL
         LJ = L_SOL(J)
         DO I = 1,NSOL
            Y(IY) = PIM(I,J,IR)
            Y(IY+NYH) = QIM(I,J,IR)
C
            PIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY)
            QIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY+NYH)
            IY = IY + 1
         END DO
      END DO
C
      X = DBLE(IRTOP)
      HBS = -1.0D0
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPAN = NPAN,1, - 1
C
         IF ( IPAN.GT.1 ) THEN
            IR_LOWER_BOUND = IRCUT(IPAN-1) + 1
         ELSE
            IR_LOWER_BOUND = JRNS1
         END IF
         JLAG1 = MAX(IR_LOWER_BOUND,IRCUT(IPAN)+1-(NLAG+1)/2)
C
         CALL FPRWFBSRAD(IGAMSOL,L_SOL,NSOL,X,Y,DY,NY,NYH,CV0,VLL,R,
     &                   DRDI,JRNSMIN,NRMAX)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRCUT(IPAN) - 1,IRCUT(IPAN-1) + 1, - 1
C
            JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            IF ( IR.EQ.(JRNS1-1) ) THEN
C
               KNS = .FALSE.
               IR_LOWER_BOUND = 1
C
               CALL FPRWFBSRAD(IGAMSOL,L_SOL,NSOL,X,Y,DY,NY,NYH,CV0,VLL,
     &                         R,DRDI,JRNSMIN,NRMAX)
C
               JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
               JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            END IF
C
            CALL FPRWFBSSTP(IGAMSOL,L_SOL,NSOL,X,HBS,Y,DY,NY,NYH,CV0,
     &                      VLL,R,DRDI,JRNSMIN,NRMAX)
C
            IY = 1
            DO J = 1,NSOL
               LJ = L_SOL(J)
               DO I = 1,NSOL
                  PIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY)
                  QIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY+NYH)
                  IY = IY + 1
               END DO
            END DO
C
         END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C --------------- if current panel is worked through, copy wavefunctions
C ------------ to first mesh point of next panel (except for last panel)
C
         IF ( IPAN.NE.1 ) THEN
            IR = IRCUT(IPAN-1)
            IY = 1
            DO J = 1,NSOL
               LJ = L_SOL(J)
               DO I = 1,NSOL
                  PIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY)
                  QIM(I,J,IR) = RIGAM(IR,LJ)*Y(IY+NYH)
                  IY = IY + 1
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DEALLOCATE (VLL,IGAMSOL,GAMSOL,L_SOL,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:FPRWFBS -> VLL'
C
      END
C*==fprwfbsstp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPRWFBSSTP(GAMSOL,L_SOL,NSOL,X,HBS,Y,DYDX,NY,NYH,CV0,
     &                      VLL,R,DRDI,JRNSMIN,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DXDY  for last mesh-point                        *
C   *   on exit:  X,Y,DXDY  updated for X = X(last) + HBS              *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:TOL_FPRWFBS
      IMPLICIT NONE
C*--FPRWFBSSTP384
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
      REAL*8 HBS,X
      INTEGER JRNSMIN,NRMAX,NSOL,NY,NYH
      COMPLEX*16 CV0(NRMAX),DYDX(NY),VLL(NSOL,NSOL,JRNSMIN:NRMAX),Y(NY)
      REAL*8 DRDI(NRMAX),GAMSOL(NSOL),R(NRMAX)
      INTEGER L_SOL(NSOL)
C
C Local variables
C
      REAL*8 ABSBB,ERRMAX,FF(NUSE),H2MID,HMID,XEST,XMID,XX(ISEQMAX)
      COMPLEX*16 BB,BB1,CC,CRAT,DD(:,:),DYSAV(:),DYY,SWAP,VV,YERR(:),
     &           YM(:),YN(:),YSAV(:),YSEQ(:),YY
      INTEGER I,ISTEP,IY,J,JJ,KK,MM,NSEQ(ISEQMAX)
C
C*** End of declarations rewritten by SPAG
C
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,
     &     1024,1536,2048,3072,4096,6144,8192,12288,16384,24576,32768,
     &     49152,65536,98304,131072,262144,524288,1048576/
C
      ALLOCATABLE DD,YM,YN,YERR,YSAV,YSEQ,DYSAV
      ALLOCATE (DD(NY,NUSE),YM(NY),YN(NY),YERR(NY),YSAV(NY))
      ALLOCATE (YSEQ(NY),DYSAV(NY))
C
      CALL CINIT(NY*NUSE,DD)
C
      DO I = 1,NY
         YSAV(I) = Y(I)
         DYSAV(I) = DYDX(I)
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
         CALL FPRWFBSRAD(GAMSOL,L_SOL,NSOL,XMID,YN,YSEQ,NY,NYH,CV0,VLL,
     &                   R,DRDI,JRNSMIN,NRMAX)
C
         H2MID = HMID + HMID
C
         DO ISTEP = 2,NSEQ(I)
            DO IY = 1,NY
               SWAP = YM(IY) + H2MID*YSEQ(IY)
               YM(IY) = YN(IY)
               YN(IY) = SWAP
            END DO
            XMID = XMID + HMID
C
            CALL FPRWFBSRAD(GAMSOL,L_SOL,NSOL,XMID,YN,YSEQ,NY,NYH,CV0,
     &                      VLL,R,DRDI,JRNSMIN,NRMAX)
C
         END DO
C
         DO IY = 1,NY
            YSEQ(IY) = 0.5D0*(YM(IY)+YN(IY)+HMID*YSEQ(IY))
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
            IF ( ERRMAX.GT.TOL_FPRWFBS ) GOTO 100
         END DO
C
C------------------------------------------------------------- converged
C
         X = X + HBS
C
         CALL FPRWFBSRAD(GAMSOL,L_SOL,NSOL,X,Y,DYDX,NY,NYH,CV0,VLL,R,
     &                   DRDI,JRNSMIN,NRMAX)
C
         DEALLOCATE (DD,YM,YN,YERR,YSAV,YSEQ,DYSAV)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<FPRWFBSSTP>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error   : ',ERRMAX
      WRITE (6,*) 'tolerance TOL_FPRWFBS : ',TOL_FPRWFBS
      WRITE (6,*) 'grid position       X : ',X
C
      END
C*==fprwfbsrad.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPRWFBSRAD(GAMSOL,L_SOL,NSOL,XBS,Y,DYDX,NY,NYH,CV0,VLL,
     &                      R,DRDI,JRNSMIN,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RADIAL_SOLVER,ONLY:EBS,CSQR,LIM,JLAG1,KNS
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--FPRWFBSRAD568
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
      INTEGER JRNSMIN,NRMAX,NSOL,NY,NYH
      REAL*8 XBS
      COMPLEX*16 CV0(NRMAX),DYDX(NY),VLL(NSOL,NSOL,JRNSMIN:NRMAX),Y(NY)
      REAL*8 DRDI(NRMAX),GAMSOL(NSOL),R(NRMAX)
      INTEGER L_SOL(NSOL)
C
C Local variables
C
      REAL*8 DRDIBS,DROVR,DX,DXSQ,GM1,GP1,RBS,WLAG(:),WVP,WVQ
      COMPLEX*16 EMVPP,EMVPP0,EMVQQ,S,SOVR2,SVP,SVQ,VBS,VLLI(:,:)
      INTEGER I,IP,IQ,J,JHIGH,JLOW,JX,K,KP,KP0,KQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WLAG,VLLI
      ALLOCATE (WLAG(NRMAX),VLLI(NSOL,NSOL))
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-12 ) THEN
         JX = NINT(XBS)
         VBS = CV0(JX)
         RBS = R(JX)
         DRDIBS = DRDI(JX)
         JLOW = 0
C
      ELSE
C
         JLOW = JLAG1
         JHIGH = JLAG1 + NLAG - 1
C
         DX = XBS - JLAG1
         DXSQ = DX*DX
         WLAG(JLOW+0) = 0.5D0*(2-3*DX+DXSQ)
         WLAG(JLOW+1) = 2*DX - DXSQ
         WLAG(JLOW+2) = 0.5D0*(-DX+DXSQ)
         VBS = 0D0
         RBS = 0D0
         DRDIBS = 0D0
         DO J = JLOW,JHIGH
            VBS = VBS + WLAG(J)*CV0(J)
            RBS = RBS + WLAG(J)*R(J)
            DRDIBS = DRDIBS + WLAG(J)*DRDI(J)
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C
      S = CSQR/((EBS-VBS)*LIM+CSQR)
      SOVR2 = (S/(RBS*RBS))*DRDIBS
      EMVPP0 = -(EBS-VBS)*DRDIBS
      EMVQQ = DRDIBS/S
      DROVR = DRDIBS/RBS
C
C---------------------------------------------------- NON-SPHERICAL CASE
C
      IF ( KNS ) THEN
C
         IF ( JLOW.EQ.0 ) THEN
C
            CALL ZCOPY(NYH,VLL(1,1,JX),1,VLLI,1)
C
         ELSE
C
            CALL FPRWFBSLAG(VLL(1,1,JLAG1),VLLI,WLAG(JLAG1),NYH,NLAG)
C
         END IF
C
C-----------------------------------------------------------------------
C
         WVP = DRDIBS
         WVQ = -DRDIBS/CSQR
C
         IP = 0
         KP0 = 0
         DO J = 1,NSOL
C
            GM1 = -(GAMSOL(J)-1.0D0)*DROVR
            GP1 = -(GAMSOL(J)+1.0D0)*DROVR
C
            DO I = 1,NSOL
               EMVPP = EMVPP0 + SOVR2*DBLE(L_SOL(I)*(L_SOL(I)+1))
               KP = KP0
               SVP = C0
               SVQ = C0
               DO K = 1,NSOL
                  KP = KP + 1
                  KQ = KP + NYH
                  SVP = SVP + VLLI(I,K)*Y(KP)
                  SVQ = SVQ + VLLI(I,K)*Y(KQ)
               END DO
C
               IP = IP + 1
               IQ = IP + NYH
               DYDX(IP) = GM1*Y(IP) + EMVQQ*Y(IQ) + WVQ*SVQ
               DYDX(IQ) = GP1*Y(IQ) + EMVPP*Y(IP) + WVP*SVP
            END DO
            KP0 = IP
         END DO
C
      ELSE
C
C-------------------------------------------------------- SPHERICAL CASE
C
         IP = 0
         DO J = 1,NSOL
C
            GM1 = -(GAMSOL(J)-1.0D0)*DROVR
            GP1 = -(GAMSOL(J)+1.0D0)*DROVR
C
            DO I = 1,NSOL
               EMVPP = EMVPP0 + SOVR2*DBLE(L_SOL(I)*(L_SOL(I)+1))
C
               IP = IP + 1
               IQ = IP + NYH
               DYDX(IP) = GM1*Y(IP) + EMVQQ*Y(IQ)
               DYDX(IQ) = GP1*Y(IQ) + EMVPP*Y(IP)
            END DO
         END DO
C
      END IF
C
      DEALLOCATE (WLAG,VLLI)
C
      END
C*==fprwfbslag.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPRWFBSLAG(VLL,VLLI,WLAG,NYH,NLAG)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--FPRWFBSLAG724
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLAG,NYH
      COMPLEX*16 VLL(NYH,NLAG),VLLI(NYH)
      REAL*8 WLAG(NLAG)
C
C Local variables
C
      INTEGER IY,J
C
C*** End of declarations rewritten by SPAG
C
      DO IY = 1,NYH
         VLLI(IY) = C0
      END DO
C
      DO J = 1,NLAG
         DO IY = 1,NYH
            VLLI(IY) = VLLI(IY) + WLAG(J)*VLL(IY,J)
         END DO
      END DO
C
      END
