C*==posfil.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSFIL(IFIL,STR,IPOS,ISTOP)
C   ********************************************************************
C   *                                                                  *
C   *   position file  IFIL  to line starting with string  STR         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--POSFIL9
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL,IPOS,ISTOP
      CHARACTER*10 STR
C
C Local variables
C
      CHARACTER*10 STRIN
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      REWIND IFIL
C
 100  CONTINUE
      READ (IFIL,FMT='(A)',END=200) STRIN
      IF ( STR.NE.STRIN ) GOTO 100
      IPOS = 1
      RETURN
C
 200  CONTINUE
      IPOS = 0
      IF ( ISTOP.NE.0 ) THEN
         WRITE (6,*) ' STOP IN <POSFIL>'
         WRITE (6,*) ' STRING ',STR,' NOT FOUND IN FILE ',IFIL
         STOP
      ELSE
         REWIND IFIL
         RETURN
      END IF
      END
