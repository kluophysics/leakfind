C*==readstr.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE READSTR(IFIL,STR80,LNG)
C   ********************************************************************
C   *                                                                  *
C   *   read string STR80 of max length 80 and skip leading blanks     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--READSTR9
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,LNG
      CHARACTER*80 STR80
C
C Local variables
C
      INTEGER I,LNGMAX
C
C*** End of declarations rewritten by SPAG
C
C
C
      LNGMAX = 80
C
      DO I = 1,LNGMAX
         STR80(I:I) = ' '
      END DO
C
      READ (IFIL,FMT='(A)') STR80
C
      DO I = 1,LNGMAX
         IF ( STR80(1:1).NE.' ' ) EXIT
         STR80(1:LNGMAX-1) = STR80(2:LNGMAX)
         STR80(LNGMAX:LNGMAX) = ' '
      END DO
C
      LNG = INDEX(STR80,' ') - 1
C
      DO I = LNG + 1,LNGMAX
         STR80(I:I) = ' '
      END DO
C
      END
