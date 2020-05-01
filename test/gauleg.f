C*==gauleg.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C   ********************************************************************
C   *                                                                  *
C   *   FIND MESH AND WEIGHT FOR GAUSS-LEGENDRE QUADRATURE             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--GAULEG11
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EPS
      PARAMETER (EPS=3.D-14)
C
C Dummy arguments
C
      INTEGER N
      REAL*8 X1,X2
      REAL*8 W(N),X(N)
C
C Local variables
C
      INTEGER I,J,M
      REAL*8 P1,P2,P3,PP,XL,XM,Z,Z1
C
C*** End of declarations rewritten by SPAG
C
      M = (N+1)/2
      XM = 0.5D0*(X2+X1)
      XL = 0.5D0*(X2-X1)
      DO I = 1,M
         Z = COS(PI*(I-0.25D0)/(N+0.5D0))
 50      CONTINUE
         P1 = 1.D0
         P2 = 0.D0
         DO J = 1,N
            P3 = P2
            P2 = P1
            P1 = ((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
         END DO
         PP = N*(Z*P1-P2)/(Z*Z-1.D0)
         Z1 = Z
         Z = Z1 - P1/PP
         IF ( ABS(Z-Z1).GT.EPS ) GOTO 50
         X(I) = XM - XL*Z
         X(N+1-I) = XM + XL*Z
         W(I) = 2.D0*XL/((1.D0-Z*Z)*PP*PP)
         W(N+1-I) = W(I)
      END DO
      END
