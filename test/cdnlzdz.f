C*==cdnlzdz.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      FUNCTION CDNLZDZ(L,Z,MODE)
C   ********************************************************************
C   *                                                                  *
C   *     d n(L,Z) / dz    analytically                                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CDNLZDZ10
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L,MODE
      COMPLEX*16 Z
      COMPLEX*16 CDNLZDZ
C
C Local variables
C
      COMPLEX*16 CNLZ
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
            CDNLZDZ = L*CNLZ(L,Z)/Z - CNLZ(L+1,Z)
         ELSE
            CDNLZDZ = (L*CNLZ(L-1,Z)-(L+1)*CNLZ(L+1,Z))/DBLE(2*L+1)
            RETURN
         END IF
      ELSE IF ( MODE.EQ.2 ) THEN
C
         IF ( L.EQ.0 ) THEN
            CDNLZDZ = L*CNLZ(L,Z)/Z - CNLZ(L+1,Z)
         ELSE
            CDNLZDZ = CNLZ(L-1,Z) - (L+1)*CNLZ(L,Z)/Z
            RETURN
         END IF
      ELSE
         CDNLZDZ = L*CNLZ(L,Z)/Z - CNLZ(L+1,Z)
      END IF
      END
