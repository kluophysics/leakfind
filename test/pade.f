C*==padeapprox.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE PADEAPPROX(YE,E,Z,P,NPADMAX,N)
C   ********************************************************************
C   *                                                                  *
C   *    Calculation of a function value at E from                     *
C   *    given pade-coefficients  P(I,J)                               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--PADEAPPROX11
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 E,YE
      INTEGER N,NPADMAX
      COMPLEX*16 P(NPADMAX,NPADMAX),Z(NPADMAX)
C
C Local variables
C
      COMPLEX*16 A(0:NPADMAX),B(0:NPADMAX)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      A(0) = 0.0D0
      A(1) = P(1,1)
      B(0) = 1.D0
      B(1) = 1.D0
      DO I = 1,N - 1
         A(I+1) = A(I) + (E-Z(I))*P(I+1,I+1)*A(I-1)
         B(I+1) = B(I) + (E-Z(I))*P(I+1,I+1)*B(I-1)
      END DO
      YE = A(N)/B(N)
      END
C*==padecoeff.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE PADECOEFF(YPAD,Z,P,NPADMAX,N)
C   ********************************************************************
C   *                                                                  *
C   *   Recursion for Pade coefficients (J.Serene)                     *
C   *   see: Vidberg and Serene                                        *
C   *        J. of Low Temp. Physics, 29, p179 (1977)                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--PADECOEFF60
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,NPADMAX
      COMPLEX*16 P(NPADMAX,NPADMAX),YPAD(NPADMAX),Z(NPADMAX)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         P(1,J) = YPAD(J)
      END DO
C
      DO J = 2,N
         DO I = 2,J
            P(I,J) = (P(I-1,I-1)-P(I-1,J))/(Z(J)-Z(I-1))/P(I-1,J)
         END DO
      END DO
      END
C
