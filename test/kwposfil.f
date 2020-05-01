C*==kwposfil.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE KWPOSFIL(IFIL,KW,LKW,ISTOP,ICASE,IPOS)
C   ********************************************************************
C   *                                                                  *
C   *   position file  IFIL  to line starting with string  KW10        *
C   *   IN:   IFIL                                                     *
C   *         KW        keyword                                        *
C   *         LKW       length of keyword   L <= 10                    *
C   *         ISTOP     0=    CONTINUE if string is not found          *
C   *                   else= STOP                                     *
C   *         ICASE     0=    ACCEPT lower and upper case characters   *
C   *                   else= DON'T ACCEPT                             *
C   *   OUT:  IPOS      0= string NOT found    1=found                 *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--KWPOSFIL17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ICASE,IFIL,IPOS,ISTOP,LKW
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*10 KWIN,KWP
      INTEGER LKWP
C
C*** End of declarations rewritten by SPAG
C
C
C
      REWIND (UNIT=IFIL,ERR=300)
C
      IF ( LKW.GT.10 ) THEN
         WRITE (6,*) ' WARNING from <KWPOSFIL> '
         WRITE (6,*) ' 10 < LKW=',LKW
         LKWP = 10
      ELSE
         LKWP = LKW
      END IF
C
      KWP = KW
      IF ( ICASE.EQ.0 ) CALL CNVTOUC(KWP,1,LKWP)
C
 100  CONTINUE
      READ (IFIL,FMT='(A)',END=200) KWIN
C
      IF ( ICASE.EQ.0 ) CALL CNVTOUC(KWIN,1,LKWP)
C
      IF ( KWP(1:LKWP).NE.KWIN(1:LKWP) ) GOTO 100
      IPOS = 1
      RETURN
C
 200  CONTINUE
      IPOS = 0
C
      IF ( ISTOP.NE.0 ) THEN
         WRITE (6,*) ' STOP in <KWPOSFIL>'
         WRITE (6,*) ' KW ',KW,' not found in file ',IFIL
         STOP
      ELSE
         REWIND IFIL
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,*) ' WARNING from <KWPOSFIL> ********************'
      WRITE (6,*) ' unit ',IFIL,' not connected -- REWIND failed'
      IPOS = 0
      END
C
