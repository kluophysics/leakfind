C*==cdjlzdz.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION CDJLZDZ(L,Z,MODE)
C   ********************************************************************
C   *                                                                  *
C   *     d j(L,Z) / dz    analytically                                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CDJLZDZ10
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L,MODE
      COMPLEX*16 Z
      COMPLEX*16 CDJLZDZ
C
C Local variables
C
      COMPLEX*16 CJLZ
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
      IF ( MODE.EQ.1 ) THEN
C
         IF ( L.EQ.0 ) THEN
C
            CDJLZDZ = L*CJLZ(L,Z)/Z - CJLZ(L+1,Z)
         ELSE
            CDJLZDZ = (L*CJLZ(L-1,Z)-(L+1)*CJLZ(L+1,Z))/DBLE(2*L+1)
            RETURN
         END IF
      ELSE IF ( MODE.EQ.2 ) THEN
C
         IF ( L.EQ.0 ) THEN
            CDJLZDZ = L*CJLZ(L,Z)/Z - CJLZ(L+1,Z)
         ELSE
            CDJLZDZ = CJLZ(L-1,Z) - (L+1)*CJLZ(L,Z)/Z
            RETURN
         END IF
      ELSE
         CDJLZDZ = L*CJLZ(L,Z)/Z - CJLZ(L+1,Z)
      END IF
      END
