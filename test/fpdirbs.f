C*==fpdirbs.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPDIRBS(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,GETIRRSOL,C,E,P,
     &                   V,B,VNS,BNS,LDAU,W,Z,NFPT,JRNS1,R,DRDI,DOVR,
     &                   VAMEG,VAMEF,IBLK,ISOLIKM,NSOLBLK,IKMSOLBLK,
     &                   IRCUT,NPAN,PRM,QRM,PIM,QIM,DPM,DQM,PR,QR,DP,DQ,
     &                   NCPLWFMAX,NPANMAX,NKMMAX,NLMFPMAX,JRNSMIN,
     &                   NRMAX,GF_CONV_RH)
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
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:KAPPA_IKM,MUEM05_IKM
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--FPDIRBS28
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
      COMPLEX*16 E,P
      LOGICAL GETIRRSOL,GF_CONV_RH,LDAU,LHS_SOL_EQ_RHS_SOL
      INTEGER IBLK,I_HAND_SIDE,JRNS1,JRNSMIN,NCPLWFMAX,NFPT,NKMMAX,
     &        NLMFPMAX,NPAN,NPANMAX,NRMAX,Z
      COMPLEX*16 B(NRMAX),BNS(JRNSMIN:NRMAX,NLMFPMAX),DP(2,2,NRMAX),
     &           DPM(NCPLWFMAX,NCPLWFMAX,NRMAX),DQ(2,2,NRMAX),
     &           DQM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           PIM(NCPLWFMAX,NCPLWFMAX,NRMAX),PR(2,2,NRMAX),
     &           PRM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           QIM(NCPLWFMAX,NCPLWFMAX,NRMAX),QR(2,2,NRMAX),
     &           QRM(NCPLWFMAX,NCPLWFMAX,NRMAX),V(NRMAX),
     &           VAMEF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VNS(JRNSMIN:NRMAX,NLMFPMAX),W(NCPLWFMAX,NCPLWFMAX)
      REAL*8 DOVR(NRMAX),DRDI(NRMAX),R(NRMAX)
      INTEGER IKMSOLBLK(NKMMAX,NKMMAX),IRCUT(0:NPANMAX),ISOLIKM(NKMMAX),
     &        NSOLBLK(NKMMAX)
C
C Local variables
C
      COMPLEX*16 ABF(:,:,:),ABG(:,:,:),AVF(:,:,:),AVG(:,:,:),CFAC,DY(:),
     &           EBS,FV,FVTR(:,:,:),GV,GVTR(:,:,:),PI(:,:,:),QI(:,:,:),
     &           SPHFUNL,SPHFUNLB,SPHFUNLBP1,SPHFUNLP1,Y(:),ZTOP
      REAL*8 CG1,CG2,CG4,CG5,CG8,CSQR,HBS,KAPFP(:),MJ,SK(2),SK1,SK2,
     &       SKFP(:),X
      COMPLEX*16 CJLZ,CNLZ
      INTEGER I,I1,IA_ERR,IFP(2),IKM1,IKM2,IPAN,IPOT,IR,IRBOT_PAN,IRTOP,
     &        IR_LOWER_BOUND,ISK1,ISK2,IY,J,JLAG1,KAP1,KAP2,L,LB(2),LB1,
     &        LB2,LBFP(:),LBJ,LFP(:),NSOL,NSOLSPH,NSPHER,NY,NYH
      INTEGER IKAPMUE
      LOGICAL KNS
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DY,LFP,LBFP,SKFP,KAPFP,GVTR,FVTR
      ALLOCATABLE PI,QI,Y
      ALLOCATABLE AVF,AVG,ABF,ABG
C
      ALLOCATE (AVF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (AVG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (ABF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (ABG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (DY(2*NKMMAX*NKMMAX),LFP(NKMMAX),LBFP(NCPLWFMAX))
      ALLOCATE (SKFP(NCPLWFMAX),KAPFP(NCPLWFMAX))
      ALLOCATE (PI(2,2,NRMAX),QI(2,2,NRMAX),Y(2*NKMMAX*NKMMAX))
C
      CSQR = C*C
      CFAC = P*C/(E+CSQR)
C
      NSOL = NSOLBLK(IBLK)
      NYH = NSOL*NSOL
      NY = NYH + NYH
      IRTOP = IRCUT(NPAN)
C
      ALLOCATE (GVTR(NSOL,NSOL,JRNSMIN:NRMAX))
      ALLOCATE (FVTR(NSOL,NSOL,JRNSMIN:NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPDIRBS -> FVTR'
C
C-----------------------------------------------------------------------
C      r-dependent potential matrix elements   < LAM | V(r) | LAM' >
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C                        LEFT HAND SIDE solution
C-----------------------------------------------------------------------
      IF ( .NOT.LHS_SOL_EQ_RHS_SOL .AND. I_HAND_SIDE.EQ.1 ) THEN
C
         DO I = 1,NSOL
            DO I1 = 1,NSOL
               DO IPOT = 2,NFPT
                  AVF(I,I1,IPOT) = DCONJG(VAMEF(I,I1,IPOT,1))
                  ABF(I,I1,IPOT) = DCONJG(VAMEF(I,I1,IPOT,2))
                  AVG(I,I1,IPOT) = DCONJG(VAMEG(I,I1,IPOT,1))
                  ABG(I,I1,IPOT) = DCONJG(VAMEG(I,I1,IPOT,2))
               END DO
            END DO
         END DO
C
      ELSE
C-----------------------------------------------------------------------
C                        RIGHT HAND SIDE solution
C-----------------------------------------------------------------------
C
         DO I = 1,NSOL
            DO I1 = 1,NSOL
               DO IPOT = 2,NFPT
                  AVF(I,I1,IPOT) = VAMEF(I,I1,IPOT,1)
                  ABF(I,I1,IPOT) = VAMEF(I,I1,IPOT,2)
                  AVG(I,I1,IPOT) = VAMEG(I,I1,IPOT,1)
                  ABG(I,I1,IPOT) = VAMEG(I,I1,IPOT,2)
               END DO
            END DO
         END DO
C
      END IF
C
      DO IR = JRNS1,IRTOP
C
         DO I = 1,NSOL
C
            DO I1 = 1,NSOL
               GV = 0D0
               FV = B(IR)*VAMEG(I,I1,1,2)
C
               DO IPOT = 2,NFPT
                  GV = GV - (AVF(I,I1,IPOT)*VNS(IR,IPOT)-ABF(I,I1,IPOT)
     &                 *BNS(IR,IPOT))
                  FV = FV + (AVG(I,I1,IPOT)*VNS(IR,IPOT)+ABG(I,I1,IPOT)
     &                 *BNS(IR,IPOT))
               END DO
               GVTR(I1,I,IR) = GV/CSQR
               FVTR(I1,I,IR) = FV
            END DO
C
            GVTR(I,I,IR) = GVTR(I,I,IR) + B(IR)*VAMEF(I,I,1,2)/CSQR
C
            IF ( LDAU ) THEN
               DO I1 = 1,NSOL
                  FVTR(I1,I,IR) = FVTR(I1,I,IR) + W(I1,I)
               END DO
            END IF
C
         END DO
      END DO
C
C=======================================================================
C      calculate regular radial wave functions for spherical regime
C=======================================================================
C
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*JRNS1,PRM)
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*JRNS1,QRM)
C
      NSPHER = JRNS1 + (NLAG+1)/2
C
      DO I1 = 1,NSOLBLK(IBLK)
C
         IKM1 = IKMSOLBLK(I1,IBLK)
         KAP1 = KAPPA_IKM(IKM1)
         MJ = MUEM05_IKM(IKM1) + 0.5D0
         KAPFP(I1) = DBLE(KAP1)
         ISK1 = ISIGN(1,KAP1)
C
         IF ( KAP1.LT.0 ) THEN
            CALL CINIT(2*2*NRMAX,PR)
            CALL CINIT(2*2*NRMAX,QR)
            CALL CINIT(2*2*NRMAX,DP)
            CALL CINIT(2*2*NRMAX,DQ)
            L = -KAP1 - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
            IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
            IFP(1) = ISOLIKM(IKM1)
            IFP(2) = ISOLIKM(IKM2)
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
C
            LB(1) = LB1
            SK(1) = SK1
            IF ( DABS(MJ).GT.L ) THEN
               NSOLSPH = 1
               LB(2) = 0
               SK(2) = 0.0D0
C
               CG2 = 0.0D0
               CG4 = 0.0D0
               CG8 = 0.0D0
C
            ELSE
C
               CG2 = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CG4 = -MJ/(KAP2+0.5D0)
               CG8 = -MJ/(-KAP2+0.5D0)
C
               NSOLSPH = 2
               LB(2) = LB2
               SK(2) = SK2
            END IF
C
            EBS = E
C
            DO I = 1,NSOLSPH
               SKFP(IFP(I)) = SK(I)
               LBFP(IFP(I)) = LB(I)
               LFP(IFP(I)) = L
            END DO
C
            CALL DIRBS(.FALSE.,C,E,L,MJ,KAP1,KAP2,P,CG1,CG2,CG4,CG5,CG8,
     &                 V,B,Z,R,DRDI,DOVR,NSPHER,IRCUT,NPAN,PR,QR,PI,QI,
     &                 DP,DQ,NRMAX,GF_CONV_RH)
C
            DO IR = 1,JRNS1
               DO J = 1,NSOLSPH
                  DO I = 1,NSOLSPH
                     PRM(IFP(I),IFP(J),IR) = PR(I,J,IR)
                     QRM(IFP(I),IFP(J),IR) = QR(I,J,IR)
                  END DO
               END DO
            END DO
C
         END IF
      END DO
C
C=======================================================================
C       calculate radial wave functions for NON-spherical regime
C=======================================================================
C
      IY = 1
      DO J = 1,NSOL
         DO I = 1,NSOL
            Y(IY) = PRM(I,J,JRNS1)
            Y(IY+NYH) = QRM(I,J,JRNS1)
            IY = IY + 1
         END DO
      END DO
C
      JLAG1 = 999999
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
         CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,X,Y,DY,DRDI,VAMEG,VAMEF,
     &                   GVTR,FVTR,NY,NYH,KNS,KAPFP,B,V,R,NLMFPMAX,
     &                   NCPLWFMAX,JRNSMIN,NRMAX)
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         DO IR = IRBOT_PAN,IRCUT(IPAN)
C
            JLAG1 = MAX(1,IR-(NLAG+1)/2+1)
            JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            CALL FPDIRBSSTP(EBS,CSQR,NSOL,JLAG1,Y,DY,NY,NYH,X,HBS,B,V,
     &                      GVTR,FVTR,KNS,KAPFP,R,DRDI,VAMEG,VAMEF,
     &                      NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C
            IY = 1
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PRM(I,J,IR) = Y(IY)
                  QRM(I,J,IR) = Y(IY+NYH)
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
               DO I = 1,NSOL
                  PRM(I,J,IR) = Y(IY)
                  QRM(I,J,IR) = Y(IY+NYH)
                  IY = IY + 1
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
      ZTOP = P*R(IR)
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
         L = LFP(J)
         LBJ = LBFP(J)
C
C----------------------------------------- convention for Green function
         IF ( GF_CONV_RH ) THEN
C                                               basis functions R and H
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
         END IF
C
         PIM(J,J,IR) = SPHFUNL*R(IR)
         QIM(J,J,IR) = CFAC*SKFP(J)*SPHFUNLB*R(IR)*C
C
         DPM(J,J,IR) = (DBLE(L+1)*SPHFUNL-ZTOP*SPHFUNLP1)*DRDI(IR)
         DQM(J,J,IR) = CFAC*SKFP(J)
     &                 *(DBLE(LBJ+1)*SPHFUNLB-ZTOP*SPHFUNLBP1)*DRDI(IR)
     &                 *C
C
      END DO
C
C ======================================================================
C
      IY = 1
      DO J = 1,NSOL
         DO I = 1,NSOL
            Y(IY) = PIM(I,J,IRTOP)
            Y(IY+NYH) = QIM(I,J,IRTOP)
            DY(IY) = DPM(I,J,IRTOP)
            DY(IY+NYH) = DQM(I,J,IRTOP)
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
         CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,X,Y,DY,DRDI,VAMEG,VAMEF,
     &                   GVTR,FVTR,NY,NYH,KNS,KAPFP,B,V,R,NLMFPMAX,
     &                   NCPLWFMAX,JRNSMIN,NRMAX)
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
               CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,X,Y,DY,DRDI,VAMEG,
     &                         VAMEF,GVTR,FVTR,NY,NYH,KNS,KAPFP,B,V,R,
     &                         NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C
               JLAG1 = MAX(IR_LOWER_BOUND,IR+2-(NLAG+1)/2)
               JLAG1 = MIN(JLAG1,IRCUT(IPAN)-NLAG+1)
C
            END IF
C
            CALL FPDIRBSSTP(EBS,CSQR,NSOL,JLAG1,Y,DY,NY,NYH,X,HBS,B,V,
     &                      GVTR,FVTR,KNS,KAPFP,R,DRDI,VAMEG,VAMEF,
     &                      NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C
            IY = 1
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PIM(I,J,IR) = Y(IY)
                  QIM(I,J,IR) = Y(IY+NYH)
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
               DO I = 1,NSOL
                  PIM(I,J,IR) = Y(IY)
                  QIM(I,J,IR) = Y(IY+NYH)
                  IY = IY + 1
               END DO
            END DO
            X = X + HBS
         END IF
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DEALLOCATE (DY,LFP,LBFP,SKFP,KAPFP,GVTR,FVTR)
C
      END
C*==fpdirbsstp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPDIRBSSTP(EBS,CSQR,NSOL,JLAG1,Y,DYDX,NY,NYH,X,HTRY,B,
     &                      V,GVTR,FVTR,KNS,KAP,R,DRDI,VAMEG,VAMEF,
     &                      NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DXDY  for last mesh-point                        *
C   *   on exit:  X,Y,DXDY  updated for X = X(last) + HTRY             *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:TOL_FPDIRBS
      IMPLICIT NONE
C*--FPDIRBSSTP492
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
      REAL*8 CSQR,HTRY,X
      COMPLEX*16 EBS
      INTEGER JLAG1,JRNSMIN,NCPLWFMAX,NLMFPMAX,NRMAX,NSOL,NY,NYH
      LOGICAL KNS
      COMPLEX*16 B(NRMAX),DYDX(NY),FVTR(NSOL,NSOL,JRNSMIN:NRMAX),
     &           GVTR(NSOL,NSOL,JRNSMIN:NRMAX),V(NRMAX),
     &           VAMEF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),Y(NY)
      REAL*8 DRDI(NRMAX),KAP(NCPLWFMAX),R(NRMAX)
C
C Local variables
C
      REAL*8 ABSBB,ERRMAX,FF(NUSE),H,H2MID,HMID,XEST,XMID,XX(ISEQMAX)
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
      H = HTRY
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
         HMID = H/NSEQ(I)
         DO IY = 1,NY
            YM(IY) = YSAV(IY)
            YN(IY) = YSAV(IY) + HMID*DYSAV(IY)
         END DO
         XMID = X + HMID
C
         CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,XMID,YN,YSEQ,DRDI,VAMEG,
     &                   VAMEF,GVTR,FVTR,NY,NYH,KNS,KAP,B,V,R,NLMFPMAX,
     &                   NCPLWFMAX,JRNSMIN,NRMAX)
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
            CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,XMID,YN,YSEQ,DRDI,VAMEG,
     &                      VAMEF,GVTR,FVTR,NY,NYH,KNS,KAP,B,V,R,
     &                      NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C
         END DO
C
         DO IY = 1,NY
            YSEQ(IY) = 0.5D0*(YM(IY)+YN(IY)+HMID*YSEQ(IY))
         END DO
C
C-----------------------------------------------------------------------
C
         XEST = H/NSEQ(I)
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
            IF ( ERRMAX.GT.TOL_FPDIRBS ) GOTO 100
         END DO
C
C------------------------------------------------------------- converged
C
         X = X + H
C
         CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,X,Y,DYDX,DRDI,VAMEG,VAMEF,
     &                   GVTR,FVTR,NY,NYH,KNS,KAP,B,V,R,NLMFPMAX,
     &                   NCPLWFMAX,JRNSMIN,NRMAX)
C
         DEALLOCATE (DD,YM,YN,YERR,YSAV,YSEQ,DYSAV)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<FPDIRBSSTP>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error   : ',ERRMAX
      WRITE (6,*) 'tolerance TOL_FPDIRBS : ',TOL_FPDIRBS
      WRITE (6,*) 'grid position       X : ',X
C
      X = X + H
C
      CALL FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,X,Y,DYDX,DRDI,VAMEG,VAMEF,
     &                GVTR,FVTR,NY,NYH,KNS,KAP,B,V,R,NLMFPMAX,NCPLWFMAX,
     &                JRNSMIN,NRMAX)
C
      END
C*==fpdirbsrad.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPDIRBSRAD(EBS,CSQR,NSOL,JLAG1,XBS,Y,DYDX,DRDI,VAMEG,
     &                      VAMEF,GVTR,FVTR,NY,NYH,KNS,KAP,B,V,R,
     &                      NLMFPMAX,NCPLWFMAX,JRNSMIN,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--FPDIRBSRAD689
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
      REAL*8 CSQR,XBS
      COMPLEX*16 EBS
      INTEGER JLAG1,JRNSMIN,NCPLWFMAX,NLMFPMAX,NRMAX,NSOL,NY,NYH
      LOGICAL KNS
      COMPLEX*16 B(NRMAX),DYDX(NY),FVTR(NSOL,NSOL,JRNSMIN:NRMAX),
     &           GVTR(NSOL,NSOL,JRNSMIN:NRMAX),V(NRMAX),
     &           VAMEF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),Y(NY)
      REAL*8 DRDI(NRMAX),KAP(NCPLWFMAX),R(NRMAX)
C
C Local variables
C
      COMPLEX*16 BBS,BPP,BQQ,EMVPP,EMVQQ,FF,FG(:),FVTRI(:,:),GF,GG,
     &           GVTRI(:,:),VBS
      REAL*8 DRDIBS,DX,DXSQ,RBS,WLAG(:)
      INTEGER I,I1,J,JHIGH,JLOW,JX,MF,MF1,MG,MG1,MGI,MGJ
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE FG,WLAG,FVTRI,GVTRI
      ALLOCATE (FG(NCPLWFMAX),WLAG(NRMAX),FVTRI(NSOL,NSOL))
      ALLOCATE (GVTRI(NSOL,NSOL))
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-12 ) THEN
         JX = NINT(XBS)
         VBS = V(JX)
         BBS = B(JX)
         RBS = R(JX)
         DRDIBS = DRDI(JX)
         JLOW = 0
C
      ELSE
         JLOW = JLAG1
         JHIGH = JLAG1 + NLAG - 1
C
         DX = XBS - JLAG1
         DXSQ = DX*DX
         WLAG(JLOW+0) = 0.5D0*(2-3*DX+DXSQ)
         WLAG(JLOW+1) = 2*DX - DXSQ
         WLAG(JLOW+2) = 0.5D0*(-DX+DXSQ)
         VBS = 0D0
         BBS = 0D0
         RBS = 0D0
         DRDIBS = 0D0
         DO J = JLOW,JHIGH
            VBS = VBS + WLAG(J)*V(J)
            BBS = BBS + WLAG(J)*B(J)
            RBS = RBS + WLAG(J)*R(J)
            DRDIBS = DRDIBS + WLAG(J)*DRDI(J)
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      EMVQQ = (EBS-VBS+CSQR)/CSQR
      EMVPP = -(EBS-VBS)
      BQQ = BBS/CSQR
      BPP = BBS
C
C---------------------------------- NOTE the order of the I- and J-loops
C
C---------------------------------------------------- NON-SPHERICAL CASE
C
      IF ( KNS ) THEN
C
         IF ( JLOW.EQ.0 ) THEN
C
            CALL ZCOPY(NYH,GVTR(1,1,JX),1,GVTRI,1)
            CALL ZCOPY(NYH,FVTR(1,1,JX),1,FVTRI,1)
C
         ELSE
C
            CALL FPDIRBSLAG(GVTR(1,1,JLAG1),FVTR(1,1,JLAG1),GVTRI,FVTRI,
     &                      WLAG(JLAG1),NYH,NLAG)
C
         END IF
C
         MGI = 0
C
         DO I = 1,NSOL
C
            GVTRI(I,I) = GVTRI(I,I) + EMVQQ
            FVTRI(I,I) = FVTRI(I,I) + EMVPP
C
            FF = KAP(I)/RBS
            GG = -FF
C
            MGI = MGI + 1
            MGJ = -NSOL
C
            DO J = 1,NSOL
               MGJ = MGJ + NSOL
               MG = MGI + MGJ
               MF = MG + NYH
C
               DYDX(MG) = GG*Y(MG)
               DYDX(MF) = FF*Y(MF)
C
               MG1 = MGJ
               DO I1 = 1,NSOL
                  MG1 = MG1 + 1
                  MF1 = MG1 + NYH
                  DYDX(MG) = DYDX(MG) + GVTRI(I1,I)*Y(MF1)
                  DYDX(MF) = DYDX(MF) + FVTRI(I1,I)*Y(MG1)
               END DO
C
            END DO
         END DO
C
      ELSE
C
C-------------------------------------------------------- SPHERICAL CASE
C
         MGI = 0
         DO I = 1,NSOL
            MGI = MGI + 1
            MGJ = -NSOL
C
            FF = KAP(I)/RBS
            GG = -FF
            GF = EMVQQ + BQQ*VAMEF(I,I,1,2)
C
            DO I1 = 1,NSOL
               FG(I1) = BPP*VAMEG(I,I1,1,2)
            END DO
            FG(I) = FG(I) + EMVPP
C
            DO J = 1,NSOL
               MGJ = MGJ + NSOL
               MG = MGI + MGJ
               MF = MG + NYH
C
               DYDX(MG) = GG*Y(MG) + GF*Y(MF)
               DYDX(MF) = FF*Y(MF)
C
               MG1 = MGJ
               DO I1 = 1,NSOL
                  MG1 = MG1 + 1
                  DYDX(MF) = DYDX(MF) + FG(I1)*Y(MG1)
               END DO
C
            END DO
         END DO
C
      END IF
C
      DO I = 1,NY
         DYDX(I) = DYDX(I)*DRDIBS
      END DO
C
      DEALLOCATE (FG,WLAG,FVTRI,GVTRI)
C
      END
C*==fpdirbslag.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPDIRBSLAG(GVTR,FVTR,GVTRI,FVTRI,WLAG,NYH,NLAG)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--FPDIRBSLAG879
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLAG,NYH
      COMPLEX*16 FVTR(NYH,NLAG),FVTRI(NYH),GVTR(NYH,NLAG),GVTRI(NYH)
      REAL*8 WLAG(NLAG)
C
C Local variables
C
      INTEGER IY,J
C
C*** End of declarations rewritten by SPAG
C
      DO IY = 1,NYH
         GVTRI(IY) = C0
         FVTRI(IY) = C0
      END DO
C
      DO J = 1,NLAG
         DO IY = 1,NYH
            GVTRI(IY) = GVTRI(IY) + WLAG(J)*GVTR(IY,J)
            FVTRI(IY) = FVTRI(IY) + WLAG(J)*FVTR(IY,J)
         END DO
      END DO
C
      END
