C*==splinel.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SPLINEL(YVAL,XVAL,NVAL,YVALPA,YVALPB,CALC_DERIVS,L,
     &                   ALPH,BET,GAM,DELT,NL)
C    *******************************************************************
C    *                                                                 *
C    * CALCULATE GENERALIZED SPLINE INTERPOLANT YVALNEW ON             *
C    * ARBITRARY MESH XVALNEW OF FUNCTION YVAL GIVEN ON                *
C    * GENERAL, NONUNIFORM MESH XVAL                                   *
C    *                                                                 *
C    * YVALPA/YVALPB: 1ST DERIVATIVES OF FUNCTION YVAL                 *
C    *                                                                 *
C    * AT 1ST/LAST MESH POINT (SUPPLIED BY SPLINE IF REQUIRED)         *
C    *                                                                 *
C    * SPLINE BASIS FUNCTIONS TO BE SUPPLIED VIA STATEMENT FUNCTIONS   *
C    * IN INCLUDE FILE                                                 *
C    *                                                                 *
C    * SEE BULIRSCH/STOER, NUMERISCHE MATHEMATIK 1, SPRINGER, P. 101 FF*
C    *                                                                 *
C    * TH                                                     19/08/98 *
C    *                                                                 *
C    * modif. by db in order to calculate interpolated values of 2     *
C    * arguments function                                              *
C    *                                                                 *
C    * parameter LAMBDA removed as it was not used any more            *
C    *                                                                 *
C    *******************************************************************
C
      IMPLICIT NONE
C*--SPLINEL29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NVALMAX,KLMAX,KUMAX,KL,KU
      PARAMETER (NVALMAX=300,KLMAX=3,KUMAX=3,KL=1,KU=1)
C
C Dummy arguments
C
      LOGICAL CALC_DERIVS
      INTEGER L,NL,NVAL
      REAL*8 ALPH(NVAL-1,0:NL),BET(NVAL-1,0:NL),DELT(NVAL-1,0:NL),
     &       GAM(NVAL-1,0:NL),XVAL(NVAL),YVAL(NVAL),YVALPA(0:NL),
     &       YVALPB(0:NL)
C
C Local variables
C
      REAL*8 ALGS(2*KLMAX+KUMAX+1,NVAL),BLGS(NVAL),CDELT1(NVALMAX-1),
     &       CDELT2(NVALMAX-1),CGAM1(NVALMAX-1),CGAM2(NVALMAX-1),DETM,
     &       DH(NVALMAX-1),DPHI(NVALMAX-1),DPSI(NVALMAX-1),DY(NVALMAX-1)
     &       ,PHIP0(NVALMAX-1),PHIPH(NVALMAX-1),PHIPP0(NVALMAX-1),
     &       PHIPPH(NVALMAX-1),PSIP0(NVALMAX-1),PSIPH(NVALMAX-1),PSIPP,
     &       PSIPP0(NVALMAX-1),PSIPPH(NVALMAX-1),X
      INTEGER I,INFO,IPIV(NVAL)
      REAL*8 PHI,PHIP,PHIPP,PSI,PSIP
      EXTERNAL DGBTRF,DGBTRS
C
C*** End of declarations rewritten by SPAG
C
      PHI(X) = X*X*X
      PHIP(X) = 3.0D0*X*X
      PHIPP(X) = 6.0D0*X
C
      PSI(X) = X*X
      PSIP(X) = 2.0D0*X
      PSIPP = 2.0D0
C
      IF ( NVAL.GT.NVALMAX ) STOP 'in <SPLINEL>: NVAL too large'
C
C     TH CALCULATE FUNCTION DERIVATIVES AT BOTH ENDS OF THE TABLE IF
C     REQUIRED
C
C
      IF ( CALC_DERIVS ) THEN
         CALL RDERSPL_N(YVAL,YVALPA,YVALPB,XVAL,NVAL,L,NL)
      ELSE
         YVALPA(L) = 0.D0
         YVALPB(L) = 0.D0
      END IF
C
C TH SET UP AUXILIARY FUNCTIONS (SEE NOTES BY TH)
C
      DO I = 1,NVAL - 1
C
         DH(I) = XVAL(I+1) - XVAL(I)
         DY(I) = (YVAL(I+1)-YVAL(I))/DH(I)
         PHIPH(I) = PHIP(DH(I))
         PHIP0(I) = PHIP(0.0D0)
         PHIPPH(I) = PHIPP(DH(I))
         PHIPP0(I) = PHIPP(0.0D0)
C
         DPHI(I) = (PHI(DH(I))-PHI(0.0D0))/DH(I)
C
         PSIPH(I) = PSIP(DH(I))
         PSIP0(I) = PSIP(0.0D0)
         PSIPPH(I) = PSIPP
         PSIPP0(I) = PSIPP
C
         DPSI(I) = (PSI(DH(I))-PSI(0.0D0))/DH(I)
C
         DETM = PSIPP0(I)*PHIPPH(I) - PSIPPH(I)*PHIPP0(I)
C
         CGAM1(I) = PHIPPH(I)/DETM
         CGAM2(I) = -PHIPP0(I)/DETM
C
         CDELT1(I) = -PSIPPH(I)/DETM
         CDELT2(I) = PSIPP0(I)/DETM
C
      END DO
C
C     TH SET UP SYSTEM OF LINEAR EQUATIONS
C
      CALL RINIT((2*KLMAX+KUMAX+1)*NVAL,ALGS)
      CALL RINIT(NVAL,BLGS)
C
      ALGS(KL+1+KU+1-1,1) = (CGAM1(1)*(PSIP0(1)-DPSI(1))+CDELT1(1)
     &                      *(PHIP0(1)-DPHI(1)))
C
      ALGS(KL+1+KU+1-2,2) = (CGAM2(1)*(PSIP0(1)-DPSI(1))+CDELT2(1)
     &                      *(PHIP0(1)-DPHI(1)))
C
      BLGS(1) = YVALPA(L) - DY(1)
C
C
      DO I = 2,NVAL - 1
C
         ALGS(KL+1+KU+I-(I-1),(I-1)) = (CGAM1(I-1)*(PSIPH(I-1)-DPSI(I-1)
     &                                 )+CDELT1(I-1)
     &                                 *(PHIPH(I-1)-DPHI(I-1)))
C
         ALGS(KL+1+KU+I-I,I) = (CGAM2(I-1)*(PSIPH(I-1)-DPSI(I-1))+CDELT2
     &                         (I-1)*(PHIPH(I-1)-DPHI(I-1))-CGAM1(I)
     &                         *(PSIP0(I)-DPSI(I))-CDELT1(I)
     &                         *(PHIP0(I)-DPHI(I)))
C
         ALGS(KL+1+KU+I-(I+1),(I+1)) = (-CGAM2(I)*(PSIP0(I)-DPSI(I))-
     &                                 CDELT2(I)*(PHIP0(I)-DPHI(I)))
         BLGS(I) = DY(I) - DY(I-1)
C
      END DO
C
C
      ALGS(KL+1+KU+NVAL-(NVAL-1),(NVAL-1))
     &   = (CGAM1(NVAL-1)*(PSIPH(NVAL-1)-DPSI(NVAL-1))+CDELT1(NVAL-1)
     &   *(PHIPH(NVAL-1)-DPHI(NVAL-1)))
      ALGS(KL+1+KU+NVAL-NVAL,NVAL) = (CGAM2(NVAL-1)*(PSIPH(NVAL-1)-DPSI(
     &                               NVAL-1))+CDELT2(NVAL-1)
     &                               *(PHIPH(NVAL-1)-DPHI(NVAL-1)))
      BLGS(NVAL) = YVALPB(L) - DY(NVAL-1)
C
C     TH LU-FACTORIZE BAND MATRIX
C
      CALL DGBTRF(NVAL,NVAL,KL,KU,ALGS,2*KLMAX+KUMAX+1,IPIV,INFO)
C
C     TH SOLVE SYSTEM OF LINEAR EQUATIONS
C
      CALL DGBTRS('N',NVAL,KL,KU,1,ALGS,2*KLMAX+KUMAX+1,IPIV,BLGS,
     &            NVALMAX,INFO)
C
      DO I = 1,NVAL - 1
C
         GAM(I,L) = CGAM1(I)*BLGS(I) + CGAM2(I)*BLGS(I+1)
         DELT(I,L) = CDELT1(I)*BLGS(I) + CDELT2(I)*BLGS(I+1)
C
         ALPH(I,L) = YVAL(I) - GAM(I,L)*PSI(0.0D0) - DELT(I,L)
     &               *PHI(0.0D0)
         BET(I,L) = DY(I) - GAM(I,L)*DPSI(I) - DELT(I,L)*DPHI(I)
C
      END DO
      END
C*==rderspl_n.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE RDERSPL_N(F,RDIFF1,RDIFFN,XVAL,N,L,NL)
C    *******************************************************************
C    *                                                                 *
C    * CALCULATE DERIVATIVES RDIFF1 AND RDIFFN IN FIRST AND LAST POINT *
C    * OF FUNCTION F GIVEN ON GENERAL, NONUNIFORM MESH XVAL            *
C    * AS TYPICALLY REQUIRED BY SPLINE INTERPOLATION ROUTINES          *
C    *                                                                 *
C    * 4-POINT FORMULA, ERROR OF ~3RD ORDER IN                         *
C    * RECIPROCAL MESH DENSITY (SEE NOTES BY TH)                       *
C    *                                                                 *
C    *                                                                 *
C    * modified by db in order to be used to calculate rdiff1(0:nl)    *
C    * and rdiff(0:nl)                                                 *
C    *                                                                 *
C    *******************************************************************
C
      IMPLICIT NONE
C*--RDERSPL_N203
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L,N,NL
      REAL*8 F(N),RDIFF1(0:NL),RDIFFN(0:NL),XVAL(N)
C
C Local variables
C
      REAL*8 DIFF_3PT
      REAL*8 F0,F1,F2,F3,H1,H2,H3
      INTEGER I,I1,I2,I3
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( N.LT.4 ) STOP 
     &              '<RDERSPL_N>: mesh has to contain at least 4 points'
C
C     PURE FORWARD DIFFERENCE FOR 1ST POINT
C
      I = 1
      I1 = 2
      I2 = 3
      I3 = 4
C
      H1 = XVAL(I1) - XVAL(I)
      H2 = XVAL(I2) - XVAL(I)
      H3 = XVAL(I3) - XVAL(I)
      F0 = F(I)
      F1 = F(I1)
      F2 = F(I2)
      F3 = F(I3)
C
      RDIFF1(L) = DIFF_3PT(H1,H2,H3,F0,F1,F2,F3)
      IF ( ABS(RDIFF1(L)).LE.1D-5 ) RDIFF1(L) = 0.D0
C
C     PURE BACKWARD DIFFERENCE FOR LAST POINT
C
      I = N
C
      I1 = N - 3
      I2 = N - 2
      I3 = N - 1
C
      H1 = XVAL(I1) - XVAL(I)
      H2 = XVAL(I2) - XVAL(I)
      H3 = XVAL(I3) - XVAL(I)
      F0 = F(I)
      F1 = F(I1)
      F2 = F(I2)
      F3 = F(I3)
C
      RDIFFN(L) = DIFF_3PT(H1,H2,H3,F0,F1,F2,F3)
C
      IF ( ABS(RDIFFN(L)).LE.1D-5 ) RDIFF1(L) = 0.D0
      END
C*==diff_3pt.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      REAL*8 FUNCTION DIFF_3PT(H1,H2,H3,F0,F1,F2,F3)
      IMPLICIT NONE
C*--DIFF_3PT276
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 F0,F1,F2,F3,H1,H2,H3
C
C*** End of declarations rewritten by SPAG
C
      DIFF_3PT = -(H1*H2+H1*H3+H3*H2)/(H1*H2*H3)
     &           *F0 + H3*H2/(H3*H2-H1*H2+H1*H1-H1*H3)
     &           /H1*F1 - H3*H1/(-H1*H3+H1*H2-H2*H2+H3*H2)
     &           /H2*F2 + H1*H2/(-H1*H3+H1*H2+H3*H3-H3*H2)/H3*F3
C
      END
