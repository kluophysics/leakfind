C*==strconfra.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION STRCONFRA(AA,X)
C   ********************************************************************
C   *                                                                  *
C   *         CONTINUED FRACTION   /D (18)/                            *
C   *                                                                  *
C   *         perform evaluation of continued fraction                 *
C   *         f(a,x) = [1/x+ (1-a)/1+ 1/x+ (2-a)/1+ ... ]              *
C   *         breaking sequence at (IMAX-a)/1+                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--STRCONFRA14
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER IMAX0
      PARAMETER (IMAX0=100)
C
C Dummy arguments
C
      REAL*8 AA,X
      REAL*8 STRCONFRA
C
C Local variables
C
      INTEGER I,IMAX,J
      REAL*8 WERT,WERT0
C
C*** End of declarations rewritten by SPAG
C
C
      IMAX = IMAX0 - 20
C
 100  CONTINUE
      IMAX = IMAX + 20
C
      WERT = IMAX/X
C
      I = IMAX
      DO J = 2,IMAX
         WERT = 1.0D0 + WERT
         WERT = X + (I-AA)/WERT
         WERT = (I-1)/WERT
         I = I - 1
      END DO
C
      WERT = 1.0D0 + WERT
      WERT = X + (1.0D0-AA)/WERT
      WERT = 1.0D0/WERT
C
      IF ( IMAX.EQ.IMAX0 ) THEN
         WERT0 = WERT
         GOTO 100
      END IF
      IF ( ABS((WERT-WERT0)/WERT).GT.1.0D-10 ) THEN
         WERT0 = WERT
         GOTO 100
      END IF
C
      STRCONFRA = WERT
      END
