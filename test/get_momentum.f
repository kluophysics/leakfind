C*==get_momentum.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GET_MOMENTUM(IREL,C,E,P)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the momentum P associated with energy E               *
C   *  take care of case Im(E) < 0                                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--GET_MOMENTUM11
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 E,P
      INTEGER IREL
C
C Local variables
C
      COMPLEX*16 CSQR
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------------- SCALAR- or FULLY relativistic
      IF ( IREL.NE.0 ) THEN
C
         CSQR = C*C
C
         IF ( DIMAG(E).GE.-1.0D-15 ) THEN
C
            P = SQRT(E*(1.0D0+E/CSQR))
C
         ELSE
C
            P = -DCONJG(SQRT(DCONJG(E)*(1.0D0+DCONJG(E)/CSQR)))
C
         END IF
C
C------------------------------------------------------ NON-relativistic
      ELSE IF ( DIMAG(E).GE.-1.0D-15 ) THEN
C
         P = SQRT(E)
C
      ELSE
C
         P = -DCONJG(SQRT(DCONJG(E)))
C
      END IF
C-----------------------------------------------------------------------
C
      END
