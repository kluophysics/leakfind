C*==cgc1.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION CGC1(J1,J,M1,M2)
C   ********************************************************************
C   *                                                                  *
C   *   clebsch gordon coefficients for  J2 = 1                        *
C   *                                                                  *
C   *   See  E.M.ROSE  RELAT. ELECTRON THEORY  TABLE 5.2               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CGC112
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER J,J1,M1,M2
      REAL*8 CGC1
C
C Local variables
C
      INTEGER M
      REAL*8 X
C
C*** End of declarations rewritten by SPAG
C
      CGC1 = 0.0D0
C
C------------------------------------------- account for selection rules
C
      IF ( J.LT.0 ) RETURN
      IF ( J1.LT.0 ) RETURN
      IF ( ABS(M1).GT.J1 ) RETURN
      IF ( ABS(M2).GT.1 ) RETURN
      IF ( ABS(J1-J).GT.1 ) RETURN
C
      M = M1 + M2
C
      IF ( J.EQ.(J1+1) ) THEN
         X = (2.0D0*J1+1.0D0)*(2.0D0*J1+2.0D0)
C
         IF ( M2.EQ.1 ) CGC1 = DSQRT((J1+M)*(J1+M+1)/X)
         IF ( M2.EQ.0 ) CGC1 = DSQRT((J1-M+1)*(J1+M+1)/(X*0.5D0))
         IF ( M2.EQ.-1 ) CGC1 = DSQRT((J1-M)*(J1-M+1)/X)
C
      ELSE IF ( J.EQ.J1 ) THEN
         X = 2.0D0*J1*(J1+1.0D0)
C
         IF ( M2.EQ.1 ) CGC1 = -DSQRT((J1+M)*(J1-M+1)/X)
         IF ( M2.EQ.0 ) CGC1 = M*DSQRT(1.0D0/(X*0.5D0))
         IF ( M2.EQ.-1 ) CGC1 = DSQRT((J1-M)*(J1+M+1)/X)
C
      ELSE IF ( J.EQ.(J1-1) ) THEN
         X = 2.0D0*J1*(2.0D0*J1+1.0D0)
C
         IF ( M2.EQ.1 ) CGC1 = DSQRT((J1-M)*(J1-M+1)/X)
         IF ( M2.EQ.0 ) CGC1 = -DSQRT((J1-M)*(J1+M)/(X*0.5D0))
         IF ( M2.EQ.-1 ) CGC1 = DSQRT((J1+M+1)*(J1+M)/X)
      END IF
C
      END
