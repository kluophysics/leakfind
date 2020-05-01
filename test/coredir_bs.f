C*==coredir_bs.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COREDIR_BS(C,E,L,MJ,WAY,VV,BB,RC,DRDIC,DOVRC,NMATCH,
     &                      NZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                      CGMD,CGO,NRC,Z)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE RADIAL SPIN-POLARISED DIRAC EQUATIONS         *
C   *   FOR THE CORE WAVE FUNCTIONS                                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      IMPLICIT NONE
C*--COREDIR_BS14
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,INVMAX
      PARAMETER (MPSMAX=20,NPEMAX=20,INVMAX=3)
      REAL*8 TOLBS
      PARAMETER (TOLBS=2.0D-9)
      INTEGER NCFMAX,NLAG
      PARAMETER (NCFMAX=8,NLAG=3)
C
C COMMON variables
C
      REAL*8 CSQR,EBS,KAP(2),XCGD(2),XCGMD(2),XCGO
      INTEGER NRADBS,NSOLBS
      COMMON /COMCOR/ EBS,CSQR,XCGD,XCGMD,XCGO,KAP,NSOLBS,NRADBS
C
C Dummy arguments
C
      REAL*8 C,CGO,E,MJ
      INTEGER L,NMATCH,NRC,NZERO,Z
      CHARACTER*3 WAY
      REAL*8 BB(NRC),CGD(2),CGMD(2),DOVRC(NRC),DP(2,2,NRC),DQ(2,2,NRC),
     &       DRDIC(NRC),FC(2,2,NRC),GC(2,2,NRC),PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),RC(NRC),VV(NRC),WP(2,2,NRC),WQ(2,2,NRC)
C
C Local variables
C
      REAL*8 A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,
     &       BC(0:NPEMAX),BETA,BOVA,BQQ,CG1,CG2,CG4,CG5,CG8,DET,DIX,
     &       DMUE,DVC,DY(NCFMAX),EMVPP,EMVQQ,FY(NCFMAX),GAM(2),GPM,HBS,
     &       PC(2,2,0:MPSMAX),QC(2,2,0:MPSMAX),RPWGPM,RR,TZ,VC(0:NPEMAX)
     &       ,W1,W2,W3,W4,W5,W6,W7,X
      INTEGER I,IR,IV,J,K,KAP1,KAP2,M,MPS,NFY,NN,NSOL
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C      PARAMETER (TOLBS=2.0D-12)
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      DVC = C
      CSQR = DVC*DVC
C
C     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
C
      TZ = 2.0D0*DBLE(Z)
C
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         VC(0) = VV(1) - (-TZ)/RC(1)
      ELSE
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
C .
      KAP1 = -L - 1
      KAP2 = +L
C
      CG1 = -MJ/(KAP1+0.5D0)
      CG5 = -MJ/(-KAP1+0.5D0)
      CGD(1) = CG1
      CGMD(1) = CG5
      KAP(1) = DBLE(KAP1)
C MB
      IF ( .NOT.FINITE_NUCLEUS ) THEN
         GAM(1) = DSQRT(KAP(1)**2-(TZ/DVC)**2)
      ELSE
         GAM(1) = DABS(KAP(1))
      END IF
C MB
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
CMBA
         IF ( .NOT.FINITE_NUCLEUS ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/DVC)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
CMBE
      END IF
C
      DO I = 1,2
         XCGD(I) = CGD(I)
         XCGMD(I) = CGMD(I)
      END DO
      XCGO = CGO
      EBS = E
      NSOLBS = NSOL
C
      IF ( WAY.EQ.'INW' ) THEN
C
C #####################################################################
C #####################################################################
C #####################################################################
C
C             INWARD INTEGRATION
C
         DMUE = SQRT(-E-E*E/CSQR)
         BOVA = -DMUE/(1.0D0+E/CSQR)
C
         IR = NZERO
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
C
C =============================================================== n ====
         HBS = -1.0D0
C
         NFY = 0
         DO J = 1,NSOL
            DO I = 1,NSOL
               FY(NFY+1) = WP(I,J,NZERO)
               FY(NFY+2) = WQ(I,J,NZERO)
               DY(NFY+1) = DP(I,J,NZERO)
               DY(NFY+2) = DQ(I,J,NZERO)
               NFY = NFY + 2
            END DO
         END DO
         X = DBLE(NZERO)
C
         NRADBS = NZERO - NLAG + 1
C
         CALL DIRBSRAD_COR(X,FY,DY,DRDIC,BB,VV,RC,NRC,NCFMAX)
C
         DO IR = NZERO - 1,NMATCH, - 1
            NRADBS = MAX(1,IR-(NLAG+1)/2+1)
            IF ( (NRADBS+NLAG-1).GT.NZERO ) NRADBS = NZERO - NLAG + 1
C
            CALL DIRBSSTP_COR(FY,DY,NFY,X,HBS,TOLBS,BB,VV,RC,DRDIC,NRC,
     &                        NCFMAX)
C
            NFY = 0
            DO J = 1,NSOL
               DO I = 1,NSOL
                  WP(I,J,IR) = FY(NFY+1)
                  WQ(I,J,IR) = FY(NFY+2)
                  NFY = NFY + 2
               END DO
            END DO
C
         END DO
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
         DO IR = NMATCH,NZERO
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
               PIW(I,J) = WP(I,J,NMATCH)
               QIW(I,J) = WQ(I,J,NMATCH)
            END DO
         END DO
C
         GOTO 99999
C
      END IF
C
C #####################################################################
C #####################################################################
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
C MBA
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
CMBE
C
C  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
C  FOR THE FIRST NN R - MESH - POINTS
C
      NN = 1
      DO IR = 1,NN
         RR = RC(IR)
C
         DO J = 1,NSOL
            RPWGPM = RR**GAM(J)
C
            DO I = 1,NSOL
               WP(I,J,IR) = PC(I,J,0)*RPWGPM
               WQ(I,J,IR) = QC(I,J,0)*RPWGPM
C
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
C
                  DP(I,J,IR) = DP(I,J,IR) + PC(I,J,M)
     &                         *RPWGPM*GPM*DOVRC(IR)
                  DQ(I,J,IR) = DQ(I,J,IR) + QC(I,J,M)
     &                         *RPWGPM*GPM*DOVRC(IR)
               END DO
C
            END DO
C
         END DO
      END DO
C
C ======================================================================
C
      NFY = 0
      DO J = 1,NSOL
         DO I = 1,NSOL
C
            FY(NFY+1) = WP(I,J,1)
            FY(NFY+2) = WQ(I,J,1)
            DY(NFY+1) = DP(I,J,1)
            DY(NFY+2) = DQ(I,J,1)
            NFY = NFY + 2
         END DO
      END DO
      X = 1.0D0
      DIX = 1.0D0
C
      DO IR = 2,NMATCH
C
         NRADBS = MAX(1,IR-(NLAG+1)/2)
         IF ( (NRADBS+NLAG-1).GT.NMATCH ) NRADBS = NMATCH - NLAG + 1
C
         CALL DIRBSSTP_COR(FY,DY,NFY,X,DIX,TOLBS,BB,VV,RC,DRDIC,NRC,
     &                     NCFMAX)
C
         NFY = 0
         DO J = 1,NSOL
            DO I = 1,NSOL
               WP(I,J,IR) = FY(NFY+1)
               WQ(I,J,IR) = FY(NFY+2)
               NFY = NFY + 2
            END DO
         END DO
C
      END DO
C
C ======================================================================
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
C
      DO IR = 1,NMATCH
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
            POW(I,J) = WP(I,J,NMATCH)
            QOW(I,J) = WQ(I,J,NMATCH)
         END DO
      END DO
C
      RETURN
C
99999 CONTINUE
      END
C*==dirbsstp_cor.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSSTP_COR(Y,DYDX,NV,X,HTRY,TOL,B,V,R,DRDI,NRMAX,
     &                        NCFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DXDY  for last mesh-point                        *
C   *   on exit:  X,Y,DXDY  updated for X = X(last) + HTRY             *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *   note: don't set NUSE    > NUSEMAX in <DIRBSRZE_COR>            *
C   *         don't set ISEQMAX > ISEQMAX in <DIRBSRZE_COR>            *
C   *         no step size adjusted in case of no convergency > STOP   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--DIRBSSTP_COR510
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ISEQMAX,NUSE
      PARAMETER (ISEQMAX=30,NUSE=7)
      REAL*8 TINYBIT
      PARAMETER (TINYBIT=(1.0D-20))
C
C Dummy arguments
C
      REAL*8 HTRY,TOL,X
      INTEGER NCFMAX,NRMAX,NV
      REAL*8 B(NRMAX),DRDI(NRMAX),DYDX(NV),R(NRMAX),V(NRMAX),Y(NV)
C
C Local variables
C
      REAL*8 CRAT,DYSAV(NCFMAX),ERRMAX,H,XEST,XSAV,YERR(NCFMAX),
     &       YSAV(NCFMAX),YSEQ(NCFMAX)
      INTEGER I,J,NSEQ(ISEQMAX)
C
C*** End of declarations rewritten by SPAG
C
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,
     &     1024,1536,2048,3072,4096,6144,8192,12288,16384,24576,32768,
     &     49152,65536/
C
      H = HTRY
      XSAV = X
C
      DO I = 1,NV
         YSAV(I) = Y(I)
         DYSAV(I) = DYDX(I)
      END DO
      DO I = 1,ISEQMAX
C
         CALL DIRBSMID_COR(YSAV,DYSAV,NV,XSAV,H,NSEQ(I),YSEQ,B,V,R,DRDI,
     &                     NRMAX,NCFMAX)
         XEST = H/NSEQ(I)
         XEST = XEST*XEST
C
         CALL DIRBSRZE_COR(I,XEST,YSEQ,Y,YERR,NV,NUSE)
C
         DO J = 1,NV
            CRAT = YERR(J)/(Y(J)+TINYBIT)
C            ERRMAX = ABS(DREAL(CRAT))+ABS(DIMAG(CRAT))
            ERRMAX = ABS(CRAT)
            IF ( ERRMAX.GT.TOL ) GOTO 100
         END DO
C
C------------------------------------------------------------- converged
C
         X = X + H
C
         CALL DIRBSRAD_COR(X,Y,DYDX,DRDI,B,V,R,NRMAX,NCFMAX)
C
         RETURN
C
C--------------------------------------------------------- NOT converged
C
 100  END DO
C
      WRITE (6,*) '<DIRBSSTP_COR>  not converged after ',ISEQMAX,
     &            ' refinements'
      WRITE (6,*) 'step size will not be adjusted !!!!!!'
      WRITE (6,*) 'max. relative error : ',ERRMAX
      WRITE (6,*) 'tolerance             ',TOL
      WRITE (6,*) 'grid position  X      ',X
      END
C*==dirbsrad_cor.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSRAD_COR(XBS,Y,DYDX,DRDI,B,V,R,NRMAX,NCFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--DIRBSRAD_COR604
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLAG
      PARAMETER (NLAG=3)
C
C COMMON variables
C
      REAL*8 CGD(2),CGMD(2),CGO,CSQR,EBS,KAP(2)
      INTEGER NRADBS,NSOLBS
      COMMON /COMCOR/ EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,NRADBS
C
C Dummy arguments
C
      INTEGER NCFMAX,NRMAX
      REAL*8 XBS
      REAL*8 B(NRMAX),DRDI(NRMAX),DYDX(NCFMAX),R(NRMAX),V(NRMAX),
     &       Y(NCFMAX)
C
C Local variables
C
      REAL*8 BBS,BPP,BQQ,DRDIBS,DROVR,DX,DXSQ,EMVPP,EMVQQ,RBS,VBS,
     &       W(NLAG)
      INTEGER I,J,K,M
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(XBS-NINT(XBS)).LT.1.0D-12 ) THEN
         I = NINT(XBS)
         VBS = V(I)
         BBS = B(I)
         RBS = R(I)
         DRDIBS = DRDI(I)
      ELSE
         DX = XBS - NRADBS
         DXSQ = DX*DX
         W(1) = 0.5D0*(2-3*DX+DXSQ)
         W(2) = 2*DX - DXSQ
         W(3) = 0.5D0*(-DX+DXSQ)
         VBS = 0D0
         BBS = 0D0
         RBS = 0D0
         DRDIBS = 0D0
         J = NRADBS - 1
         DO I = 1,NLAG
            J = J + 1
            VBS = VBS + W(I)*V(J)
            BBS = BBS + W(I)*B(J)
            RBS = RBS + W(I)*R(J)
            DRDIBS = DRDIBS + W(I)*DRDI(J)
         END DO
      END IF
C-----------------------------------------------------------------------
C
      EMVPP = -(EBS-VBS)*DRDIBS
      EMVQQ = -EMVPP/CSQR + DRDIBS
C
      BPP = BBS*DRDIBS
      BQQ = BPP/CSQR
C
      DROVR = DRDIBS/RBS
C
      M = 0
      DO J = 1,NSOLBS
         DO I = 1,NSOLBS
            K = 3 - 4*(I-1)
C
            DYDX(M+1) = -KAP(I)*Y(M+1)*DROVR + (EMVQQ+BQQ*CGMD(I))
     &                  *Y(M+2)
C
            DYDX(M+2) = KAP(I)*Y(M+2)*DROVR + (EMVPP+BPP*CGD(I))*Y(M+1)
     &                  + BPP*CGO*Y(M+K)
            M = M + 2
         END DO
      END DO
      END
C*==dirbsmid_cor.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSMID_COR(Y,DYDX,NV,XS,HTOT,NSTEP,YOUT,B,V,R,DRDI,
     &                        NRMAX,NCFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   modified midpoint step to support the  Burlisch-Stoer method   *
C   *   on exit:  the incremented variable is in   YOUT                *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.3                            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--DIRBSMID_COR749
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 HTOT,XS
      INTEGER NCFMAX,NRMAX,NSTEP,NV
      REAL*8 B(NRMAX),DRDI(NRMAX),DYDX(NV),R(NRMAX),V(NRMAX),Y(NV),
     &       YOUT(NV)
C
C Local variables
C
      REAL*8 H,H2,SWAP,X,YM(NCFMAX),YN(NCFMAX)
      INTEGER I,IR
C
C*** End of declarations rewritten by SPAG
C
      H = HTOT/NSTEP
      DO I = 1,NV
         YM(I) = Y(I)
         YN(I) = Y(I) + H*DYDX(I)
      END DO
      X = XS + H
C
      CALL DIRBSRAD_COR(X,YN,YOUT,DRDI,B,V,R,NRMAX,NCFMAX)
C
      H2 = 2.D0*H
      DO IR = 2,NSTEP
         DO I = 1,NV
            SWAP = YM(I) + H2*YOUT(I)
            YM(I) = YN(I)
            YN(I) = SWAP
         END DO
         X = X + H
         CALL DIRBSRAD_COR(X,YN,YOUT,DRDI,B,V,R,NRMAX,NCFMAX)
      END DO
      DO I = 1,NV
         YOUT(I) = 0.5D0*(YM(I)+YN(I)+H*YOUT(I))
      END DO
C
      END
C*==dirbsrze_cor.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBSRZE_COR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
C   ********************************************************************
C   *                                                                  *
C   *   diagonal rational function extrapolation to support the        *
C   *   Burlisch-Stoer method                                          *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--DIRBSRZE_COR814
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NCFMAX,ISEQMAX,NUSEMAX
      PARAMETER (NCFMAX=8,ISEQMAX=30,NUSEMAX=7)
C
C Dummy arguments
C
      INTEGER IEST,NUSE,NV
      REAL*8 XEST
      REAL*8 DY(NV),YEST(NV),YZ(NV)
C
C Local variables
C
      REAL*8 B,B1,C,D(NCFMAX,NUSEMAX),DDY,FX(NUSEMAX),V,X(ISEQMAX),YY
      INTEGER J,K,M1
      SAVE B,B1,C,D,DDY,FX,J,K,M1,V,X,YY
C
C*** End of declarations rewritten by SPAG
C
C
C
C*** End of declarations rewritten by SPAG
C
      X(IEST) = XEST
      IF ( IEST.EQ.1 ) THEN
         DO J = 1,NV
            YZ(J) = YEST(J)
            D(J,1) = YEST(J)
            DY(J) = YEST(J)
         END DO
      ELSE
         M1 = MIN(IEST,NUSE)
         DO K = 1,M1 - 1
            FX(K+1) = X(IEST-K)/XEST
         END DO
         DO J = 1,NV
            YY = YEST(J)
            V = D(J,1)
            C = YY
            D(J,1) = YY
            DO K = 2,M1
               B1 = FX(K)*V
               B = B1 - C
               IF ( B.NE.0D0 ) THEN
                  B = (C-V)/B
                  DDY = C*B
                  C = B1*B
               ELSE
                  DDY = V
               END IF
               V = D(J,K)
               D(J,K) = DDY
               YY = YY + DDY
            END DO
            DY(J) = DDY
            YZ(J) = YY
         END DO
      END IF
      END
