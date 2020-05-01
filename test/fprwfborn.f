C*==fprwfborn.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPRWFBORN(SPNWGT,GETIRRSOL,IT,VNS,BNS,LDAU,W,VAMEG,
     &                     IBLK,PRM,PIM,PRL,PHL,SBLK,SMT0L,TBLK,TMT0L,
     &                     CBLK,GBLK,L_LM,DPRIRTOP)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE RADIAL DIRAC EQUATIONS WITHIN THE         *
C   *   SCALAR RELATIVISTIC APPROXIMATION                              *
C   *                                                                  *
C   *   using the Born series scheme and RH convention for the GF      *
C   *                                                                  *
C   *   NOTE: the routine expects the potential terms to be COMPLEX    *
C   *                                                                  *
C   *   *BLK  denotes matrices restricted to present block  IBLK       *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CALCMODE,ONLY:TOL_FPRWFBORN,LIM_FPRWFBORN
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_RMESH,ONLY:NRMAX,JRNSMIN,JRNS1,JRCUT,NPAN,DRDI
      USE MOD_ANGMOM,ONLY:NKMMAX,NL,NLMAX,WKM2,IPIVKM
      USE MOD_TYPES,ONLY:NCPLWFMAX,NLMFPMAX,IMT,NFPT,IKMSOLBLK,NSOLBLK
      IMPLICIT NONE
C*--FPRWFBORN26
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
      LOGICAL GETIRRSOL,LDAU
      INTEGER IBLK,IT
      COMPLEX*16 SPNWGT
      COMPLEX*16 BNS(JRNSMIN:NRMAX,NLMFPMAX),CBLK(NCPLWFMAX,NCPLWFMAX),
     &           DPRIRTOP(NCPLWFMAX,NCPLWFMAX),GBLK(NCPLWFMAX,NCPLWFMAX)
     &           ,PHL(NRMAX,NL),PIM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           PRL(NRMAX,NL),PRM(NCPLWFMAX,NCPLWFMAX,NRMAX),
     &           SBLK(NCPLWFMAX,NCPLWFMAX),SMT0L(NLMAX),
     &           TBLK(NCPLWFMAX,NCPLWFMAX),TMT0L(NLMAX),
     &           VAMEG(NCPLWFMAX,NCPLWFMAX,NLMFPMAX,2),
     &           VNS(JRNSMIN:NRMAX,NLMFPMAX),W(NCPLWFMAX,NCPLWFMAX)
      INTEGER L_LM(NKMMAX)
C
C Local variables
C
      REAL*8 ADP,AXP,A_ERROR,B_ERROR,C_ERROR,D_ERROR
      COMPLEX*16 A_IJ,B_IJ,CINTWGT,CPRE,CSUM,C_IJ,DBLK(:,:),DV_P(:,:,:),
     &           D_IJ,FUNA(:,:,:),FUNB(:,:,:),GBLKINV(:,:),GKJ,HBLK(:,:)
     &           ,INTA(:,:,:),INTA_LAST(:),INTB(:,:,:),INTB_LAST(:),
     &           PH0(:,:),PR0(:,:),P_KJ,VLL(:,:,:),WGT_VLL,XIRR(:,:,:),
     &           XL,XP,XREG(:,:,:)
      INTEGER I,IA_ERR,IBORN,IL,IM,IPOT,IR,IRNSBOT,IRTOP,J,K,L,LM,NSOL,
     &        NSOLSQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DBLK,HBLK,GBLKINV,VLL,DV_P,FUNA,FUNB
      ALLOCATABLE INTA,INTB,XREG,XIRR,PR0,PH0,INTA_LAST,INTB_LAST
C
      IM = IMT(IT)
      NSOL = NSOLBLK(IBLK,IT)
      NSOLSQ = NSOL*NSOL
      IRTOP = JRCUT(NPAN(IM),IM)
C
      ALLOCATE (HBLK(NCPLWFMAX,NCPLWFMAX),DBLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (GBLKINV(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (VLL(NSOL,NSOL,JRNSMIN:NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPRWFBS -> VLL'
C
      IRTOP = JRCUT(NPAN(IM),IM)
      IF ( MOD(JRCUT(1,IM)-JRNS1(IM),2).EQ.0 ) THEN
         IRNSBOT = JRNS1(IM)
      ELSE
         IRNSBOT = JRNS1(IM) + 1
      END IF
C
      ALLOCATE (INTA_LAST(NSOL),INTB_LAST(NSOL))
      ALLOCATE (FUNA(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (FUNB(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (INTA(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (INTB(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (DV_P(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (XREG(NSOL,NSOL,IRNSBOT:NRMAX))
      ALLOCATE (XIRR(NSOL,NSOL,IRNSBOT:NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPRWFBORN -> DV_P'
C
C-----------------------------------------------------------------------
C      r-dependent potential matrix elements   < LAM | V(r) | LAM' >
C-----------------------------------------------------------------------
C include ALL prefactors {..} and weights for simpson integration [..]:
C               { -ip     }  [ (1/12) dr/di ]   for  DV_P
C
C NOTE:  the wave functions PHO contain the prefactor  ip
C
C-----------------------------------------------------------------------
C
      CPRE = -1D0/12D0
C
      DO IR = IRNSBOT,IRTOP
         DO I = 1,NSOL
            DO J = 1,NSOL
               CSUM = 0D0
               DO IPOT = 2,NFPT(IT)
                  CSUM = CSUM + VAMEG(I,J,IPOT,1)
     &                   *(VNS(IR,IPOT)+SPNWGT*BNS(IR,IPOT))
               END DO
               VLL(I,J,IR) = CSUM
C              VLL(I,J,IR) = 0d0
            END DO
         END DO
      END DO
C
      IF ( LDAU ) THEN
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO I = 1,NSOL
                  VLL(I,J,IR) = VLL(I,J,IR) + W(J,I)
               END DO
            END DO
         END DO
      END IF
C
      DO IR = IRNSBOT,IRTOP
         CINTWGT = CPRE*DRDI(IR,IM)
         DO I = 1,NSOL
            DO J = 1,NSOL
               WGT_VLL = CINTWGT*VLL(I,J,IR)
               DV_P(I,J,IR) = WGT_VLL
            END DO
         END DO
      END DO
C
C=======================================================================
C           radial wave functions for spherical potential
C=======================================================================
C
      ALLOCATE (PR0(NSOL,NRMAX))
      ALLOCATE (PH0(NSOL,NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FPRWFBORN -> QH0'
C
      PR0(1:NSOL,1:NRMAX) = C0
      PH0(1:NSOL,1:NRMAX) = C0
C
      DO I = 1,NSOL
         LM = IKMSOLBLK(I,IBLK,IT)
         L = L_LM(LM)
         IL = L + 1
C
         DO IR = 1,IRTOP
            PR0(I,IR) = PRL(IR,IL)
            PH0(I,IR) = PHL(IR,IL)
         END DO
C
      END DO
C=======================================================================
C                                            NOTE: PR0, PH0 are diagonal
C
      XREG(:,:,:) = C0
      XIRR(:,:,:) = C0
C
      DO IR = IRNSBOT,IRTOP
         DO J = 1,NSOL
            DO I = 1,NSOL
               XREG(I,J,IR) = XREG(I,J,IR) + PR0(I,IR)*DV_P(I,J,IR)
               XIRR(I,J,IR) = XIRR(I,J,IR) + PH0(I,IR)*DV_P(I,J,IR)
            END DO
         END DO
      END DO
C
C=======================================================================
C               calculate     REGULAR    wave functions
C=======================================================================
C                use array   PIM   to calclate  \bar{P}
C
      PIM(:,:,:) = C0
C
      DO IR = 1,IRTOP
         DO I = 1,NSOL
            PIM(I,I,IR) = PR0(I,IR)
         END DO
      END DO
C
      INTA_LAST(1:NSOL) = 0D0
      INTB_LAST(1:NSOL) = 0D0
C
      DO IBORN = 1,LIM_FPRWFBORN
C
C----------------------------------------- set up matrices ABAR and BBAR
C                                      use INWARD integration for ABAR !
C
         FUNA(:,:,:) = C0
         FUNB(:,:,:) = C0
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               DO K = 1,NSOL
C
                  P_KJ = PIM(K,J,IR)
C
                  DO I = 1,NSOL
C
                     FUNA(I,J,IR) = FUNA(I,J,IR) + XIRR(I,K,IR)*P_KJ
C
                     FUNB(I,J,IR) = FUNB(I,J,IR) + XREG(I,K,IR)*P_KJ
C
                  END DO
               END DO
            END DO
         END DO
C
         CALL FPINTVEC_INW(FUNA,INTA,NSOLSQ,IRNSBOT,IRTOP,NPAN(IM),
     &                     JRCUT(0,IM))
C
         CALL FPINTVEC_OUT(FUNB,INTB,NSOLSQ,IRNSBOT,IRTOP,NPAN(IM),
     &                     JRCUT(0,IM))
C
C------------------------------------------------- update wave functions
C                                                     r = r_ns ... r_top
C
         PIM(1:NSOL,1:NSOL,IRNSBOT:IRTOP) = C0
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               I = J
               PIM(I,I,IR) = PR0(I,IR)
C
               DO I = 1,NSOL
C
                  A_IJ = INTA(I,J,IRNSBOT) - INTA(I,J,IR)
                  B_IJ = INTB(I,J,IR)
C
                  PIM(I,J,IR) = PIM(I,J,IR) - PR0(I,IR)*A_IJ + PH0(I,IR)
     &                          *B_IJ
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
         IF ( A_ERROR.LE.TOL_FPRWFBORN .AND. B_ERROR.LE.TOL_FPRWFBORN )
     &        EXIT
C
      END DO
C
      IF ( IPRINT.GT.2 ) WRITE (6,99001) IT,'regular',LIM_FPRWFBORN,
     &                          A_ERROR,B_ERROR
C
C------------------------------------------------- update wave functions
C                                                         r = r ... r_ns
C
      PIM(1:NSOL,1:NSOL,1:(IRNSBOT-1)) = C0
C
      DO IR = 1,(IRNSBOT-1)
         DO I = 1,NSOL
            PIM(I,I,IR) = PR0(I,IR)
         END DO
      END DO
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
C-----------------------------------------------------------------------
C          transform wave function   P = \bar{P} x GBLK
C-----------------------------------------------------------------------
C
      PRM(:,:,:) = C0
C
      DO IR = 1,IRTOP
         DO J = 1,NSOL
            DO K = 1,NSOL
               GKJ = GBLK(K,J)
               IF ( ABS(DREAL(GKJ))+ABS(DIMAG(GKJ)).GT.EPS ) THEN
                  DO I = 1,NSOL
                     PRM(I,J,IR) = PRM(I,J,IR) + PIM(I,K,IR)*GKJ
                  END DO
               END IF
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C         single site t-matrix   t = t0 + delta t
C                  alfa-matrix   a = s0 * delta a
C-----------------------------------------------------------------------
C
      TBLK(1:NSOL,1:NSOL) = -MATMUL(INTB(1:NSOL,1:NSOL,IRTOP),
     &                      GBLK(1:NSOL,1:NSOL))
C
      DO I = 1,NSOL
C
         LM = IKMSOLBLK(I,IBLK,IT)
         IL = L_LM(LM) + 1
C
         TBLK(I,I) = TMT0L(IL) + TBLK(I,I)
C
         DO J = 1,NSOL
            SBLK(I,J) = SMT0L(IL)*GBLK(I,J)
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C         radial derivative at shere boundary from radial equations
C         not used any more - activate derivation for testing
C-----------------------------------------------------------------------
C
      DPRIRTOP(:,:) = C0
C
C     CALL CHECK_DFDI(PRM,NSOL,IRTOP,DPRIRTOP,NCPLWFMAX)
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C=======================================================================
C           calculate     IRRREGULAR    wave functions
C=======================================================================
C
      PIM(:,:,:) = C0
C
      DO IR = 1,IRTOP
         DO I = 1,NSOL
            PIM(I,I,IR) = PH0(I,IR)
         END DO
      END DO
C
      INTA_LAST(1:NSOL) = 0D0
      INTB_LAST(1:NSOL) = 0D0
C
      DO IBORN = 1,LIM_FPRWFBORN
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
C
                  DO I = 1,NSOL
C
                     FUNA(I,J,IR) = FUNA(I,J,IR) + XIRR(I,K,IR)*P_KJ
C
                     FUNB(I,J,IR) = FUNB(I,J,IR) + XREG(I,K,IR)*P_KJ
C
                  END DO
               END DO
            END DO
         END DO
C
         CALL FPINTVEC_INW(FUNA,INTA,NSOLSQ,IRNSBOT,IRTOP,NPAN(IM),
     &                     JRCUT(0,IM))
C
         CALL FPINTVEC_INW(FUNB,INTB,NSOLSQ,IRNSBOT,IRTOP,NPAN(IM),
     &                     JRCUT(0,IM))
C
C------------------------------------------------- update wave functions
C
         CALL CINIT(NCPLWFMAX*NCPLWFMAX*IRTOP,PIM)
C
         DO IR = IRNSBOT,IRTOP
            DO J = 1,NSOL
               PIM(J,J,IR) = PH0(J,IR)
               DO I = 1,NSOL
                  C_IJ = INTA(I,J,IR)
                  D_IJ = -INTB(I,J,IR)
C
                  PIM(I,J,IR) = PIM(I,J,IR) + PR0(I,IR)*C_IJ + PH0(I,IR)
     &                          *D_IJ
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
         IF ( C_ERROR.LE.TOL_FPRWFBORN .AND. D_ERROR.LE.TOL_FPRWFBORN )
     &        EXIT
C
      END DO
C
      IF ( IPRINT.GT.2 ) WRITE (6,99001) IT,'irregular',LIM_FPRWFBORN,
     &                          C_ERROR,D_ERROR
C
C------------------------------------------------- force: \bar A = - D^T
C     DBLK(1:NSOL,1:NSOL)  ==  - INTB(1:NSOL,1:NSOL,IRNSBOT)
C----------------------------------------------- complete wave functions
C
      DO IR = 1,IRNSBOT - 1
         DO J = 1,NSOL
            PIM(J,J,IR) = PH0(J,IR)
            DO I = 1,NSOL
               C_IJ = INTA(I,J,IRNSBOT)
               D_IJ = DBLK(I,J)
C
               PIM(I,J,IR) = PIM(I,J,IR) + PR0(I,IR)*C_IJ + PH0(I,IR)
     &                       *D_IJ
C
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C  weight for regular functions  R^0 in charge density CHI = GAM x c^T
C-----------------------------------------------------------------------
C
      HBLK(1:NSOL,1:NSOL) = TRANSPOSE(INTA(1:NSOL,1:NSOL,IRNSBOT))
C
      CBLK(1:NSOL,1:NSOL) = MATMUL(GBLK(1:NSOL,1:NSOL),HBLK(1:NSOL,1:
     &                      NSOL))
C
C-----------------------------------------------------------------------
99001 FORMAT (/,5X,'<FPRWS_BORN>: IT =',I3,'  calculating the ',A,
     &        ' wave function',/,5X,' max. number of Born iterations (',
     &        I2,') exceeded',/,5X,'last errors:   A =',E14.5,'   B =',
     &        E14.5,/)
99002 FORMAT ('iter ',i3,20(2E25.15,2x))
      END
C*==check_dfdi.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE CHECK_DFDI(F,NSOL,N,D,MSOL)
C   ********************************************************************
C   *                                                                  *
C   *                  get d F / d i   at point N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CHECK_DFDI509
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER MSOL,N,NSOL
      COMPLEX*16 D(MSOL,MSOL),F(NSOL,NSOL,N)
C
C Local variables
C
      COMPLEX*16 DFDI
      INTEGER I,J
      REAL*8 RAT
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NSOL
         DO J = 1,NSOL
            DFDI = ((11.D0*F(I,J,N)+9.D0*F(I,J,N-2))
     &             -(18.D0*F(I,J,N-1)+2.D0*F(I,J,N-3)))/6.0D0
            RAT = ABS(1D0-D(I,J)/DFDI)
            RAT = ABS(D(I,J)-DFDI)
            IF ( RAT.GT.1D-7 ) THEN
               WRITE (6,99001) I,J,D(I,J),DFDI,RAT,NSOL
               WRITE (6,99002)
               IF ( NSOL.NE.1 ) STOP 'XXXXXXXXXXX'
            END IF
         END DO
      END DO
C
99001 FORMAT (/,2I4,2E20.12,/,8x,2E20.12,F20.16,I7)
99002 FORMAT (4(1x,79('*'),/))
      END
