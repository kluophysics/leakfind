C*==fpdirborn.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPDIRBORN(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,IT,GETIRRSOL,
     &                     C,E,P,V,B,VNS,BNS,LDAU,W,NFPT,JRNS1,DRDI,
     &                     VAMEG,VAMEF,IBLK,ISOLIKM,NSOLBLK,IKMSOLBLK,
     &                     IRCUT,NPAN,PRM,QRM,PIM,QIM,PR,QR,DP,DQ,TBLK,
     &                     TBLK0,TSST0,SBLK,SBLK0,SSST0,CBLK,DBLK,GBLK,
     &                     NCWF_SOL_SPH,IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,
     &                     JG_SPH,JF_SPH,NCPLWFMAX,NPANMAX,NKMMAX,
     &                     NLMFPMAX,JRNSMIN,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *   using the Born series scheme and RH convention for the GF      *
C   *                                                                  *
C   *   *BLK  denotes matrices restricted to present block  IBLK       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:WKM2,IPIVKM,KAPPA_IKM,MUEM05_IKM
      USE MOD_CALCMODE,ONLY:TOL_FPDIRBORN,LIM_FPDIRBORN
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--FPDIRBORN25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EPS
      PARAMETER (EPS=1D-16)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 E,P
      LOGICAL GETIRRSOL,LDAU,LHS_SOL_EQ_RHS_SOL
      INTEGER IBLK,IT,I_HAND_SIDE,JRNS1,JRNSMIN,NCPLWFMAX,NFPT,NKMMAX,
     &        NLMFPMAX,NPAN,NPANMAX,NRMAX
      COMPLEX*16 B(NRMAX),BNS(JRNSMIN:NRMAX,NLMFPMAX),
     &           CBLK(NCPLWFMAX,NCPLWFMAX),DBLK(NCPLWFMAX,NCPLWFMAX),
     &           DP(2,2,NRMAX),DQ(2,2,NRMAX),GBLK(NCPLWFMAX,NCPLWFMAX),
     &           JF_SPH(NRMAX,2,NKMMAX),JG_SPH(NRMAX,2,NKMMAX),
     &           PIM(NCPLWFMAX,NCPLWFMAX,NRMAX),PR(2,2,NRMAX),
     &           PRM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           QIM(NCPLWFMAX,NCPLWFMAX,NRMAX),QR(2,2,NRMAX),
     &           QRM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           SBLK(NCPLWFMAX,NCPLWFMAX),SBLK0(NCPLWFMAX,NCPLWFMAX),
     &           SSST0(NKMMAX,NKMMAX),TBLK(NCPLWFMAX,NCPLWFMAX),
     &           TBLK0(NCPLWFMAX,NCPLWFMAX),TSST0(NKMMAX,NKMMAX),
     &           V(NRMAX),VAMEF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VNS(JRNSMIN:NRMAX,NLMFPMAX),W(NCPLWFMAX,NCPLWFMAX),
     &           ZF_SPH(NRMAX,2,NKMMAX),ZG_SPH(NRMAX,2,NKMMAX)
      REAL*8 DRDI(NRMAX)
      INTEGER IKMSOLBLK(NKMMAX,NKMMAX),IKM_CWFSOL_SPH(2,NKMMAX),
     &        IRCUT(0:NPANMAX),ISOLIKM(NKMMAX),NCWF_SOL_SPH(NKMMAX),
     &        NSOLBLK(NKMMAX)
C
C Local variables
C
      COMPLEX*16 ABF(:,:,:),ABG(:,:,:),AVF(:,:,:),AVG(:,:,:),A_KJ,B_KJ,
     &           CINTWGT,CPRE,C_KJ,DV_P(:,:,:),DV_Q(:,:,:),D_KJ,
     &           FUNA(:,:,:),FUNB(:,:,:),FV,GBLKINV(:,:),GKJ,GV,
     &           HBLK(:,:),INTA(:,:,:),INTA_LAST(:),INTB(:,:,:),
     &           INTB_LAST(:),PH(:,:,:),PH0(:,:,:),PR0(:,:,:),P_KJ,
     &           QH(:,:,:),QH0(:,:,:),QR0(:,:,:),Q_KJ,WINP(:,:),
     &           XIRR(:,:,:),XL,XP,XREG(:,:,:),YIRR(:,:,:),YREG(:,:,:)
      REAL*8 ADP,AXP,A_ERROR,B_ERROR,CSQR,C_ERROR,D_ERROR,MJ
      INTEGER I,I1,IA_ERR,IBORN,IFP(2),IKM,IKM1,IKM2,IKMCB(2),IPOT,IR,
     &        IRNSBOT,IRTOP,J,JKM,K,KAP1,KAP2,KBOT(:),KTOP(:),L,NSOL,
     &        NSOLSPH,NSOLSQ
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DV_Q,DV_P,FUNA,FUNB,KBOT,KTOP
      ALLOCATABLE PH,QH,INTA,INTB,XREG,XIRR,YREG,YIRR
      ALLOCATABLE PR0,QR0,PH0,QH0,GBLKINV,HBLK
      ALLOCATABLE INTA_LAST,INTB_LAST
      ALLOCATABLE AVF,AVG,ABF,ABG,WINP
C
      ALLOCATE (AVF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (AVG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (ABF(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (ABG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX))
      ALLOCATE (PH(2,2,NRMAX),QH(2,2,NRMAX))
      ALLOCATE (WINP(NCPLWFMAX,NCPLWFMAX))
      WINP = C0
C
C***********************************************************************
C
      CSQR = C*C
C
      NSOL = NSOLBLK(IBLK)
      NSOLSQ = NSOL*NSOL
C
      ALLOCATE (KBOT(NSOL),KTOP(NSOL))
      ALLOCATE (HBLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (GBLKINV(NCPLWFMAX,NCPLWFMAX))
C
      KBOT(1:NSOL) = 1
      KTOP(1:NSOL) = NSOL
C     KBOT(1:NSOL) = NSOL
C     KTOP(1:NSOL) = 1
C
      IRTOP = IRCUT(NPAN)
      IF ( MOD(IRCUT(1)-JRNS1,2).EQ.0 ) THEN
         IRNSBOT = JRNS1
      ELSE
         IRNSBOT = JRNS1 + 1
      END IF
C
      ALLOCATE (INTA_LAST(NSOL),INTB_LAST(NSOL))
      ALLOCATE (FUNA(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (FUNB(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (INTA(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (INTB(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (DV_Q(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (DV_P(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (XREG(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (XIRR(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (YREG(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (YIRR(NSOL,NSOL,IRNSBOT:NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPDIRBORN -> DV_P'
C
C-----------------------------------------------------------------------
C      r-dependent potential matrix elements   < LAM | V(r) | LAM' >
C-----------------------------------------------------------------------
C include ALL prefactors {..} and weights for simpson integration [..]:
C               { -ip     }  [ (1/12) dr/di ]   for  DV_P
C               { -ip/c^2 }  [ (1/12) dr/di ]   for  DV_Q
C
C NOTE:  the wave functions PHO contain the prefactor  ip
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C                        LEFT HAND SIDE solution
C-----------------------------------------------------------------------
      IF ( .NOT.LHS_SOL_EQ_RHS_SOL .AND. I_HAND_SIDE.EQ.1 ) THEN
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               DO IPOT = 2,NFPT
                  AVF(I,J,IPOT) = DCONJG(VAMEF(I,J,IPOT,1))
                  ABF(I,J,IPOT) = DCONJG(VAMEF(I,J,IPOT,2))
                  AVG(I,J,IPOT) = DCONJG(VAMEG(I,J,IPOT,1))
                  ABG(I,J,IPOT) = DCONJG(VAMEG(I,J,IPOT,2))
               END DO
               WINP(I,J) = W(J,I)
            END DO
         END DO
C
      ELSE
C-----------------------------------------------------------------------
C                        RIGHT HAND SIDE solution
C-----------------------------------------------------------------------
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               DO IPOT = 2,NFPT
                  AVF(I,J,IPOT) = VAMEF(I,J,IPOT,1)
                  ABF(I,J,IPOT) = VAMEF(I,J,IPOT,2)
                  AVG(I,J,IPOT) = VAMEG(I,J,IPOT,1)
                  ABG(I,J,IPOT) = VAMEG(I,J,IPOT,2)
               END DO
               WINP(I,J) = W(I,J)
            END DO
         END DO
C
      END IF
C
C
      CPRE = -1D0/12D0
C
      DO IR = IRNSBOT,IRTOP
C
         DO J = 1,NSOL
C
            DO I = 1,NSOL
               GV = 0D0
               FV = 0D0
C
               DO IPOT = 2,NFPT
                  GV = GV + (AVF(I,J,IPOT)*VNS(IR,IPOT)-ABF(I,J,IPOT)
     &                 *BNS(IR,IPOT))
                  FV = FV + (AVG(I,J,IPOT)*VNS(IR,IPOT)+ABG(I,J,IPOT)
     &                 *BNS(IR,IPOT))
               END DO
               DV_Q(I,J,IR) = GV/CSQR
               DV_P(I,J,IR) = FV
            END DO
C
            IF ( LDAU ) THEN
               DO I = 1,NSOL
                  DV_P(I,J,IR) = DV_P(I,J,IR) + WINP(I,J)
               END DO
            END IF
C
         END DO
C
         CINTWGT = CPRE*DRDI(IR)
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               DV_P(I,J,IR) = DV_P(I,J,IR)*CINTWGT
               DV_Q(I,J,IR) = DV_Q(I,J,IR)*CINTWGT
            END DO
         END DO
C
      END DO
C
C=======================================================================
C           calculate radial wave functions for spherical regime
C=======================================================================
C
C-----------------------------------------------------------------------
C          normalized   regular   (PR,QR) -> J - ip SUM H * t
C
      IF ( IBLK.EQ.1 ) CALL DIRBS_NORM(IT,C,E,P,V,B,IRTOP,PR,QR,PH,QH,
     &                                 DP,DQ,NCWF_SOL_SPH,
     &                                 IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,
     &                                 JG_SPH,JF_SPH,TSST0,SSST0,
     &                                 NSOLSPH)
C
      ALLOCATE (PR0(NSOL,NSOL,NRMAX))
      ALLOCATE (QR0(NSOL,NSOL,NRMAX))
      ALLOCATE (PH0(NSOL,NSOL,NRMAX))
      ALLOCATE (QH0(NSOL,NSOL,NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPDIRBORN -> QH0'
C
      PR0(1:NSOL,1:NSOL,1:NRMAX) = C0
      QR0(1:NSOL,1:NSOL,1:NRMAX) = C0
      PH0(1:NSOL,1:NSOL,1:NRMAX) = C0
      QH0(1:NSOL,1:NSOL,1:NRMAX) = C0
C
      DO I1 = 1,NSOLBLK(IBLK)
C
         IKM1 = IKMSOLBLK(I1,IBLK)
C
         KAP1 = KAPPA_IKM(IKM1)
         MJ = MUEM05_IKM(IKM1) + 0.5D0
C
         IF ( KAP1.LT.0 ) THEN
C
            L = -KAP1 - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
            IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
            IFP(1) = ISOLIKM(IKM1)
            IFP(2) = ISOLIKM(IKM2)
            IKMCB(1) = IKM1
            IKMCB(2) = IKM2
C
            IF ( ABS(MJ).GT.L ) THEN
               NSOLSPH = 1
            ELSE
               NSOLSPH = 2
            END IF
C
C-----------------------------------------------------------------------
C
            DO IR = 1,IRTOP
               DO J = 1,NSOLSPH
                  IKM = IKMCB(J)
                  DO I = 1,NSOLSPH
                     PR0(IFP(I),IFP(J),IR) = ZG_SPH(IR,I,IKM)
                     QR0(IFP(I),IFP(J),IR) = ZF_SPH(IR,I,IKM)
C
                     PH0(IFP(I),IFP(J),IR) = JG_SPH(IR,I,IKM)
                     QH0(IFP(I),IFP(J),IR) = JF_SPH(IR,I,IKM)
                  END DO
               END DO
            END DO
C
C-----------------------------------------------------------------------
            DO J = 1,NSOLSPH
               DO I = 1,NSOLSPH
                  KBOT(IFP(I)) = MIN(KBOT(IFP(I)),IFP(J))
                  KTOP(IFP(I)) = MAX(KTOP(IFP(I)),IFP(J))
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END IF
      END DO
C
C----------------------------------- note indexing of PR0, PH0, QR0, QH0
C
      XREG(1:NSOL,1:NSOL,IRNSBOT:NRMAX) = C0
      XIRR(1:NSOL,1:NSOL,IRNSBOT:NRMAX) = C0
      YREG(1:NSOL,1:NSOL,IRNSBOT:NRMAX) = C0
      YIRR(1:NSOL,1:NSOL,IRNSBOT:NRMAX) = C0
C
      DO IR = IRNSBOT,IRTOP
         DO J = 1,NSOL
            DO I = 1,NSOL
               DO K = KBOT(I),KTOP(I)
                  XREG(I,J,IR) = XREG(I,J,IR) + PR0(K,I,IR)*DV_P(K,J,IR)
                  XIRR(I,J,IR) = XIRR(I,J,IR) + PH0(K,I,IR)*DV_P(K,J,IR)
C
                  YREG(I,J,IR) = YREG(I,J,IR) + QR0(K,I,IR)*DV_Q(K,J,IR)
                  YIRR(I,J,IR) = YIRR(I,J,IR) + QH0(K,I,IR)*DV_Q(K,J,IR)
               END DO
            END DO
         END DO
      END DO
C
C=======================================================================
C               calculate     REGULAR    wave functions
C=======================================================================
C            use arrays   PIM, QIM   to calclate  \bar{P} and \bar{Q}
C
      PIM(1:NSOL,1:NSOL,1:IRTOP) = PR0(1:NSOL,1:NSOL,1:IRTOP)
      QIM(1:NSOL,1:NSOL,1:IRTOP) = QR0(1:NSOL,1:NSOL,1:IRTOP)
C
      INTA_LAST(1:NSOL) = 0D0
      INTB_LAST(1:NSOL) = 0D0
C
      DO IBORN = 1,LIM_FPDIRBORN
C
C----------------------------------------- set up matrices ABAR and BBAR
C                               use INWARD integration for ABAR !
C
         FUNA(1:NSOL,1:NSOL,IRNSBOT:IRTOP) = C0
         FUNB(1:NSOL,1:NSOL,IRNSBOT:IRTOP) = C0
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO K = 1,NSOL
                  P_KJ = PIM(K,J,IR)
                  Q_KJ = QIM(K,J,IR)
C
                  DO I = 1,NSOL
C
                     FUNA(I,J,IR) = FUNA(I,J,IR) + XIRR(I,K,IR)
     &                              *P_KJ + YIRR(I,K,IR)*Q_KJ
C
                     FUNB(I,J,IR) = FUNB(I,J,IR) + XREG(I,K,IR)
     &                              *P_KJ + YREG(I,K,IR)*Q_KJ
C
                  END DO
               END DO
            END DO
         END DO
C
         CALL FPINTVEC_INW(FUNA,INTA,NSOLSQ,IRNSBOT,IRTOP,NPAN,IRCUT)
C
         CALL FPINTVEC_OUT(FUNB,INTB,NSOLSQ,IRNSBOT,IRTOP,NPAN,IRCUT)
C
C------------------------------------------------- update wave functions
C
         PIM(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
     &      = PR0(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
         QIM(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
     &      = QR0(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO K = KBOT(J),KTOP(J)
C
                  A_KJ = INTA(K,J,IRNSBOT) - INTA(K,J,IR)
                  B_KJ = INTB(K,J,IR)
C
                  DO I = 1,NSOL
                     PIM(I,J,IR) = PIM(I,J,IR) - PR0(I,K,IR)
     &                             *A_KJ + PH0(I,K,IR)*B_KJ
C
                     QIM(I,J,IR) = QIM(I,J,IR) - QR0(I,K,IR)
     &                             *A_KJ + QH0(I,K,IR)*B_KJ
                  END DO
C
               END DO
            END DO
         END DO
C
C----------------------------------------------------- check convergency
C
         A_ERROR = 0D0
         B_ERROR = 0D0
         DO I = 1,NSOL
            XP = INTA(I,I,IRNSBOT)
            XL = INTA_LAST(I)
            INTA_LAST(I) = XP
            AXP = ABS(DREAL(XP)) + ABS(DIMAG(XP))
            IF ( AXP.GT.1D-16 ) THEN
               ADP = ABS(DREAL(XP-XL)) + ABS(DIMAG(XP-XL))
               A_ERROR = MAX(A_ERROR,ADP)
            END IF
C
            XP = INTB(I,I,IRTOP)
            XL = INTB_LAST(I)
            INTB_LAST(I) = XP
            AXP = ABS(DREAL(XP)) + ABS(DIMAG(XP))
            IF ( AXP.GT.1D-16 ) THEN
               ADP = ABS(DREAL(XP-XL)) + ABS(DIMAG(XP-XL))
               B_ERROR = MAX(B_ERROR,ADP)
            END IF
         END DO
         IF ( IPRINT.GT.2 ) WRITE (6,99002) IBORN,A_ERROR,B_ERROR
         IF ( A_ERROR.LE.TOL_FPDIRBORN .AND. B_ERROR.LE.TOL_FPDIRBORN )
     &        EXIT
C
      END DO
C
      IF ( IPRINT.GT.2 ) WRITE (6,99001) IT,'regular',LIM_FPDIRBORN,
     &                          A_ERROR,B_ERROR
C
C=======================================================================
C
C
C-----------------------------------------------------------------------
C                      force: \bar A = - D^T
C-----------------------------------------------------------------------
C
      DBLK(1:NSOL,1:NSOL) = -TRANSPOSE(INTA(1:NSOL,1:NSOL,IRNSBOT))
C
C-----------------------------------------------------------------------
C                  transformation matrix gamma  =  GBLK
C-----------------------------------------------------------------------
C
      GBLKINV(1:NSOL,1:NSOL) = -INTA(1:NSOL,1:NSOL,IRNSBOT)
C
      DO I = 1,NSOL
         GBLKINV(I,I) = 1D0 + GBLKINV(I,I)
      END DO
C
      CALL CMATINV3(NSOL,NCPLWFMAX,IPIVKM,GBLKINV,WKM2,GBLK)
C
C------------------------------ alternative inversion routine by H. AKAI
C
C     CALL CMATINV(NSOL,NCPLWFMAX,GBLKINV,GBLK)
C
C-----------------------------------------------------------------------
C          transform wave function   P = \bar{P} x GBLK
C                                    Q = \bar{Q} x GBLK
C-----------------------------------------------------------------------
C
C     CALL  CMATSTR('G-MATRIX',GBLK,NSOL,NCPLWFMAX,0,0,1,1.0D-9,6)
C
      PRM(:,:,:) = C0
      QRM(:,:,:) = C0
C
      DO J = 1,NSOL
         DO K = 1,NSOL
            GKJ = GBLK(K,J)
            IF ( ABS(DREAL(GKJ))+ABS(DIMAG(GKJ)).GT.EPS ) THEN
               DO IR = 1,IRTOP
                  DO I = 1,NSOL
                     PRM(I,J,IR) = PRM(I,J,IR) + PIM(I,K,IR)*GKJ
                     QRM(I,J,IR) = QRM(I,J,IR) + QIM(I,K,IR)*GKJ
                  END DO
               END DO
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C         single site t-matrix   t = t0 + delta t
C-----------------------------------------------------------------------
C
      DO J = 1,NSOL
         JKM = IKMSOLBLK(J,IBLK)
         DO I = 1,NSOL
            IKM = IKMSOLBLK(I,IBLK)
            TBLK0(I,J) = TSST0(IKM,JKM)
         END DO
      END DO
C
      HBLK(1:NSOL,1:NSOL) = -MATMUL(INTB(1:NSOL,1:NSOL,IRTOP),
     &                      GBLK(1:NSOL,1:NSOL))
C
      TBLK(1:NSOL,1:NSOL) = TBLK0(1:NSOL,1:NSOL) + HBLK(1:NSOL,1:NSOL)
C
C-----------------------------------------------------------------------
C                  alfa-matrix   a = s0 * delta a
C-----------------------------------------------------------------------
C
      DO J = 1,NSOL
         JKM = IKMSOLBLK(J,IBLK)
         DO I = 1,NSOL
            IKM = IKMSOLBLK(I,IBLK)
            SBLK0(I,J) = SSST0(IKM,JKM)
         END DO
      END DO
C
      SBLK(1:NSOL,1:NSOL) = MATMUL(SBLK0(1:NSOL,1:NSOL),GBLK(1:NSOL,1:
     &                      NSOL))
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C=======================================================================
C           calculate     IRRREGULAR    wave functions
C=======================================================================
C
C
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*IRTOP,PIM)
      CALL CINIT(NCPLWFMAX*NCPLWFMAX*IRTOP,QIM)
C
      PIM(1:NSOL,1:NSOL,1:IRTOP) = PH0(1:NSOL,1:NSOL,1:IRTOP)
      QIM(1:NSOL,1:NSOL,1:IRTOP) = QH0(1:NSOL,1:NSOL,1:IRTOP)
C
      INTA_LAST(1:NSOL) = 0D0
      INTB_LAST(1:NSOL) = 0D0
C
      DO IBORN = 1,LIM_FPDIRBORN
C
C----------------------------------------- set up matrices CBAR and DBAR
C                                                 use INWARD integration
C
         FUNA(1:NSOL,1:NSOL,IRNSBOT:IRTOP) = C0
         FUNB(1:NSOL,1:NSOL,IRNSBOT:IRTOP) = C0
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO K = 1,NSOL
C
                  P_KJ = PIM(K,J,IR)
                  Q_KJ = QIM(K,J,IR)
C
                  DO I = 1,NSOL
C
                     FUNA(I,J,IR) = FUNA(I,J,IR) + XIRR(I,K,IR)
     &                              *P_KJ + YIRR(I,K,IR)*Q_KJ
C
                     FUNB(I,J,IR) = FUNB(I,J,IR) + XREG(I,K,IR)
     &                              *P_KJ + YREG(I,K,IR)*Q_KJ
C
                  END DO
               END DO
            END DO
         END DO
C
         CALL FPINTVEC_INW(FUNA,INTA,NSOLSQ,IRNSBOT,IRTOP,NPAN,IRCUT)
C
         CALL FPINTVEC_INW(FUNB,INTB,NSOLSQ,IRNSBOT,IRTOP,NPAN,IRCUT)
C
C------------------------------------------------- update wave functions
C
         PIM(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
     &      = PH0(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
         QIM(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
     &      = QH0(1:NSOL,1:NSOL,IRNSBOT:IRTOP)
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO K = KBOT(J),KTOP(J)
                  C_KJ = INTA(K,J,IR)
                  D_KJ = -INTB(K,J,IR)
C
                  DO I = 1,NSOL
                     PIM(I,J,IR) = PIM(I,J,IR) + PR0(I,K,IR)
     &                             *C_KJ + PH0(I,K,IR)*D_KJ
                     QIM(I,J,IR) = QIM(I,J,IR) + QR0(I,K,IR)
     &                             *C_KJ + QH0(I,K,IR)*D_KJ
                  END DO
C
               END DO
            END DO
         END DO
C
C----------------------------------------------------- check convergency
C
         C_ERROR = 0D0
         D_ERROR = 0D0
         DO I = 1,NSOL
            XP = INTA(I,I,IRNSBOT)
            XL = INTA_LAST(I)
            INTA_LAST(I) = XP
            AXP = ABS(DREAL(XP)) + ABS(DIMAG(XP))
            IF ( AXP.GT.1D-16 ) THEN
               ADP = ABS(DREAL(XP-XL)) + ABS(DIMAG(XP-XL))
               C_ERROR = MAX(C_ERROR,ADP)
            END IF
C
            XP = INTB(I,I,IRNSBOT)
            XL = INTB_LAST(I)
            INTB_LAST(I) = XP
            AXP = ABS(DREAL(XP)) + ABS(DIMAG(XP))
            IF ( AXP.GT.1D-16 ) THEN
               ADP = ABS(DREAL(XP-XL)) + ABS(DIMAG(XP-XL))
               D_ERROR = MAX(D_ERROR,ADP)
            END IF
         END DO
         IF ( IPRINT.GT.2 ) WRITE (6,99002) IBORN,C_ERROR,D_ERROR
         IF ( C_ERROR.LE.TOL_FPDIRBORN .AND. D_ERROR.LE.TOL_FPDIRBORN )
     &        EXIT
C
      END DO
C
      IF ( IPRINT.GT.2 ) WRITE (6,99001) IT,'irregular',LIM_FPDIRBORN,
     &                          C_ERROR,D_ERROR
C
C------------------------------------------------- force: \bar A = - D^T
C     DBLK(1:NSOL,1:NSOL)  ==  - INTB(1:NSOL,1:NSOL,IRNSBOT)
C----------------------------------------------- complete wave functions
C
C     CALL  CMATSTR('D-MATRIX',DBLK,NSOL,NCPLWFMAX,0,0,1,1.0D-9,6)
C     CALL  CMATSTR('B-MATRIX',INTB(1,1,IRNSBOT),nsol,NCPLWFMAX,0,0,1,
C    &                   1.0D-9,6)
Cc      DBLK(1:NSOL,1:NSOL) = DBLK(1:NSOL,1:NSOL)
Cc     &                      + INTB(1:NSOL,1:NSOL,IRNSBOT)
Cc      CALL  CMATSTR('DIFFER  ',DBLK,NSOL,NCPLWFMAX,0,0,1,1.0D-9,6)
C
      DO IR = 1,IRNSBOT - 1
         DO J = 1,NSOL
            DO K = KBOT(J),KTOP(J)
               C_KJ = INTA(K,J,IRNSBOT)
               D_KJ = DBLK(K,J)
C
               DO I = 1,NSOL
                  PIM(I,J,IR) = PIM(I,J,IR) + PR0(I,K,IR)
     &                          *C_KJ + PH0(I,K,IR)*D_KJ
                  QIM(I,J,IR) = QIM(I,J,IR) + QR0(I,K,IR)
     &                          *C_KJ + QH0(I,K,IR)*D_KJ
               END DO
C
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C  weight for regular functions  R^0 in charge density CHI = GAM x c^T
C  store matrix c
C-----------------------------------------------------------------------
C
      CBLK(1:NSOL,1:NSOL) = INTA(1:NSOL,1:NSOL,IRNSBOT)
C
C      HBLK(1:NSOL,1:NSOL) = TRANSPOSE(INTA(1:NSOL,1:NSOL,IRNSBOT))
C
C      CBLK(1:NSOL,1:NSOL) = MATMUL(GBLK(1:NSOL,1:NSOL),HBLK(1:NSOL,1:
C     &                      NSOL))
C
C-----------------------------------------------------------------------
99001 FORMAT (/,5X,'<FPDIR_BORN>: IT =',I3,'  calculating the ',A,
     &        ' wave function',/,5X,' max. number of Born iterations (',
     &        I2,') exceeded',/,5X,'last errors:   A =',E14.5,'   B =',
     &        E14.5,/)
99002 FORMAT ('iter ',i3,20(2E25.15,2x))
      END
C*==fpintvec_out.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPINTVEC_OUT(F,C,N,IRNSBOT,IRTOP,NPAN,IRCUT)
C   ********************************************************************
C   *                                                                  *
C   *   perform OUTWARD integration for a N-dim vector function  F(i)  *
C   *   NOTE: the nxn matrix function in the calling routine           *
C   *         is treated as a vector function with N=n*n               *
C   *                                                                  *
C   *            C(i) = int[IRNSBOT..i] d i'  F(i')                    *
C   *                                                                  *
C   *   via Simpson rule for a mesh with NPAN kinks at IRCUT           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--FPINTVEC_OUT672
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRNSBOT,IRTOP,N,NPAN
      COMPLEX*16 C(N,IRNSBOT:IRTOP),F(N,IRNSBOT:IRTOP)
      INTEGER IRCUT(0:NPAN)
C
C Local variables
C
      COMPLEX*16 F0_X2,F0_X4,F0_X5,F1_X16,F1_X2,F1_X4,F1_X8,F2_X2,F2_X4
      INTEGER I,IPAN,IR,IRBOT_PAN
C
C*** End of declarations rewritten by SPAG
C
      C(1:N,IRNSBOT) = C0
C
C-----------------------------------------------------------------------
C                         loop over panels
C-----------------------------------------------------------------------
      DO IPAN = 1,NPAN
C
         IF ( IPAN.EQ.1 ) THEN
            IRBOT_PAN = IRNSBOT
         ELSE
            IRBOT_PAN = IRCUT(IPAN-1) + 1
            C(1:N,IRBOT_PAN) = C(1:N,IRBOT_PAN-1)
         END IF
C
C-------------------------------- deal with ODD number of steps in panel
         IF ( MOD(IRCUT(IPAN)-IRBOT_PAN,2).EQ.1 ) THEN
C
            IR = IRBOT_PAN
C
            DO I = 1,N
C
               F0_X2 = F(I,IR) + F(I,IR)
               F0_X4 = F0_X2 + F0_X2
               F0_X5 = F0_X4 + F(I,IR)
C
               F1_X2 = F(I,IR+1) + F(I,IR+1)
               F1_X4 = F1_X2 + F1_X2
               F1_X8 = F1_X4 + F1_X4
C
               C(I,IR+1) = C(I,IR) + F0_X5 + F1_X8 - F(I,IR+2)
C
            END DO
C
            IRBOT_PAN = IRBOT_PAN + 1
C
         END IF
C
         DO IR = IRBOT_PAN,IRCUT(IPAN) - 2,2
C
            DO I = 1,N
C
               F0_X2 = F(I,IR) + F(I,IR)
               F0_X4 = F0_X2 + F0_X2
               F0_X5 = F0_X4 + F(I,IR)
C
               F1_X2 = F(I,IR+1) + F(I,IR+1)
               F1_X4 = F1_X2 + F1_X2
               F1_X8 = F1_X4 + F1_X4
               F1_X16 = F1_X8 + F1_X8
C
               F2_X2 = F(I,IR+2) + F(I,IR+2)
               F2_X4 = F2_X2 + F2_X2
C
               C(I,IR+1) = C(I,IR) + F0_X5 + F1_X8 - F(I,IR+2)
               C(I,IR+2) = C(I,IR) + F0_X4 + F1_X16 + F2_X4
C
            END DO
C
         END DO
C
      END DO
C
      END
C*==fpintvec_inw.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPINTVEC_INW(F,C,N,IRNSBOT,IRTOP,NPAN,IRCUT)
C   ********************************************************************
C   *                                                                  *
C   *   perform INWARD integration for a N-dim vector function  F(i)   *
C   *   NOTE: the nxn matrix function in the calling routine           *
C   *         is treated as a vector function with N=n*n               *
C   *                                                                  *
C   *            C(i) = int[i..IRTOP] d i'  F(i')                      *
C   *                                                                  *
C   *   via Simpson rule for a mesh with NPAN kinks at IRCUT           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--FPINTVEC_INW778
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRNSBOT,IRTOP,N,NPAN
      COMPLEX*16 C(N,IRNSBOT:IRTOP),F(N,IRNSBOT:IRTOP)
      INTEGER IRCUT(0:NPAN)
C
C Local variables
C
      COMPLEX*16 F0_X2,F0_X4,F0_X5,F1_X16,F1_X2,F1_X4,F1_X8,F2_X2,F2_X4
      INTEGER I,IPAN,IR,IRBOT_PAN,IRTOP_PAN
C
C*** End of declarations rewritten by SPAG
C
      C(1:N,IRCUT(NPAN)) = C0
C
C-----------------------------------------------------------------------
C                         loop over panels
C-----------------------------------------------------------------------
      DO IPAN = NPAN,1, - 1
C
         IRTOP_PAN = IRCUT(IPAN)
         IF ( IPAN.EQ.1 ) THEN
            IRBOT_PAN = IRNSBOT
         ELSE
            IRBOT_PAN = IRCUT(IPAN-1) + 1
         END IF
         IF ( IPAN.NE.NPAN ) C(1:N,IRTOP_PAN) = C(1:N,IRTOP_PAN+1)
C
C-------------------------------- deal with ODD number of steps in panel
         IF ( MOD(IRTOP_PAN-IRBOT_PAN,2).EQ.1 ) THEN
C
            IR = IRTOP_PAN
C
            DO I = 1,N
C
               F0_X2 = F(I,IR) + F(I,IR)
               F0_X4 = F0_X2 + F0_X2
               F0_X5 = F0_X4 + F(I,IR)
C
               F1_X2 = F(I,IR-1) + F(I,IR-1)
               F1_X4 = F1_X2 + F1_X2
               F1_X8 = F1_X4 + F1_X4
C
               C(I,IR-1) = C(I,IR) + F0_X5 + F1_X8 - F(I,IR-2)
C
            END DO
C
            IRTOP_PAN = IRTOP_PAN - 1
C
         END IF
C
         DO IR = IRTOP_PAN,IRBOT_PAN + 2, - 2
C
            DO I = 1,N
C
               F0_X2 = F(I,IR) + F(I,IR)
               F0_X4 = F0_X2 + F0_X2
               F0_X5 = F0_X4 + F(I,IR)
C
               F1_X2 = F(I,IR-1) + F(I,IR-1)
               F1_X4 = F1_X2 + F1_X2
               F1_X8 = F1_X4 + F1_X4
               F1_X16 = F1_X8 + F1_X8
C
               F2_X2 = F(I,IR-2) + F(I,IR-2)
               F2_X4 = F2_X2 + F2_X2
C
               C(I,IR-1) = C(I,IR) + F0_X5 + F1_X8 - F(I,IR-2)
               C(I,IR-2) = C(I,IR) + F0_X4 + F1_X16 + F2_X4
C
            END DO
C
         END DO
C
      END DO
C
      END
