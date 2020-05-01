C*==spline.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SPLINE(X,Y,IRTOP,YP1,YPN,Y2)
C
C Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
C y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
C the first derivative of the interpolating function at points 1 and n,
C respectively, this routine returns an array y2(1:n) of length n which
C contains the second derivatives of the interpolating function at the
C tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
C the routine is signaled to set the corresponding boundary condition for a
C natural spline with zero second derivative on that boundary.
C Parameter: nmax is the largest anticipiated value of n
C
      IMPLICIT NONE
C*--SPLINE15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NMAX
      PARAMETER (NMAX=500)
C
C Dummy arguments
C
      INTEGER IRTOP
      REAL*8 YP1,YPN
      REAL*8 X(IRTOP),Y(IRTOP),Y2(IRTOP)
C
C Local variables
C
      INTEGER IR
      REAL*8 P,QN,SIG,U(NMAX),UN
C
C*** End of declarations rewritten by SPAG
C
      IF ( YP1.GT..99D30 ) THEN
         Y2(1) = 0.0D0
         U(1) = 0.0D0
      ELSE
         Y2(1) = -0.5D0
         U(1) = (3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      END IF
C
      DO IR = 2,IRTOP - 1
         SIG = (X(IR)-X(IR-1))/(X(IR+1)-X(IR-1))
         P = SIG*Y2(IR-1) + 2.0D0
         Y2(IR) = (SIG-1.0D0)/P
         U(IR) = (6.0D0*((Y(IR+1)-Y(IR))/(X(IR+1)-X(IR))-(Y(IR)-Y(IR-1))
     &           /(X(IR)-X(IR-1)))/(X(IR+1)-X(IR-1))-SIG*U(IR-1))/P
      END DO
C
      IF ( YPN.GT..99D30 ) THEN
         QN = 0.0D0
         UN = 0.0D0
      ELSE
         QN = 0.5D0
         UN = (3.0D0/(X(IRTOP)-X(IRTOP-1)))
     &        *(YPN-(Y(IRTOP)-Y(IRTOP-1))/(X(IRTOP)-X(IRTOP-1)))
      END IF
C
      Y2(IRTOP) = (UN-QN*U(IRTOP-1))/(QN*Y2(IRTOP-1)+1.0D0)
C
      DO IR = IRTOP - 1,1, - 1
         Y2(IR) = Y2(IR)*Y2(IR+1) + U(IR)
      END DO
C
      END
