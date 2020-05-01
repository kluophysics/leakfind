C*==getlrecreal.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GETLRECREAL(IOTMP,IPRINT,LRECREAL8)
C   ********************************************************************
C   *                                                                  *
C   *   determine the MACHINE DEPENDENT record length used for a       *
C   *   REAL double precission variable needed to open a               *
C   *   direct access file                                             *
C   *                                                                  *
C   *   trial and error: writing and reading of 2 reals should work    *
C   *                    writing and reading of 3 reals should fail    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--GETLRECREAL14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IOTMP,IPRINT,LRECREAL8
C
C Local variables
C
      INTEGER IOL
      CHARACTER*60 TEXT
      REAL*8 TOL,X1,X2,X3,XX1,XX2,XX3
C
C*** End of declarations rewritten by SPAG
C
      X1 = SQRT(2D0)
      X2 = 2.0D0*X1
      X3 = 3.0D0*X1
      TOL = 1D-15
C
C--------------------------------------------------------------- trial A
C
      LRECREAL8 = 8
      TEXT = 'a double precission real*8 is stored in 8 Bytes'
C
      IOL = 2*LRECREAL8
C
      OPEN (IOTMP,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='DIRECT',
     &      RECL=IOL)
C
      WRITE (IOTMP,REC=1,ERR=100) X1,X2
      WRITE (IOTMP,REC=2,ERR=100) X2,X3
      READ (IOTMP,REC=1,ERR=100) XX1,XX2
C
      IF ( ABS(X1-XX1).GT.TOL .OR. ABS(X2-XX2).GT.TOL ) THEN
         WRITE (6,99002) LRECREAL8
         WRITE (6,*) 'trial A: X1  X2   written  ',X1,X2
         WRITE (6,*) 'trial A: XX1 XX2  read     ',XX1,XX2
      END IF
C
      WRITE (IOTMP,REC=1,ERR=400) X1,X2,X3
      WRITE (IOTMP,REC=2,ERR=400) X2,X3,X2
      READ (IOTMP,REC=1,ERR=400) XX1,XX2,XX3
      CLOSE (IOTMP)
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,*) 'trial A: X1  X2  X3   written  ',X1,X2,X3
         WRITE (6,*) 'trial A: XX1 XX2 XX3  read     ',XX1,XX2,XX3
      END IF
C
C--------------------------------------------------------------- trial B
C
 100  CONTINUE
      LRECREAL8 = 2
      TEXT = 'a double precission real*8 is stored in 2 words a 32 bits'
      IF ( IPRINT.GT.0 ) WRITE (6,99002) LRECREAL8
C
      IOL = 2*LRECREAL8
C
      OPEN (IOTMP,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='DIRECT',
     &      RECL=IOL)
C
      WRITE (IOTMP,REC=1,ERR=200) X1,X2
      WRITE (IOTMP,REC=2,ERR=200) X2,X3
      READ (IOTMP,REC=1,ERR=200) XX1,XX2
      IF ( ABS(X1-XX1).GT.TOL .OR. ABS(X2-XX2).GT.TOL ) THEN
         WRITE (6,99002) LRECREAL8
         WRITE (6,*) 'trial B: X1  X2   written  ',X1,X2
         WRITE (6,*) 'trial B: XX1 XX2  read     ',XX1,XX2
      END IF
C
      WRITE (IOTMP,REC=1,ERR=400) X1,X2,X3
      WRITE (IOTMP,REC=2,ERR=400) X2,X3,X2
      READ (IOTMP,REC=1,ERR=400) XX1,XX2,XX3
      CLOSE (IOTMP)
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,*) 'trial B: X1  X2  X3   written  ',X1,X2,X3
         WRITE (6,*) 'trial B: XX1 XX2 XX3  read     ',XX1,XX2,XX3
      END IF
C
C--------------------------------------------------------------- trial C
C
 200  CONTINUE
      LRECREAL8 = 1
      TEXT = 'a double precission real*8 is stored in word a 64 bits'
      IF ( IPRINT.GT.0 ) WRITE (6,99002) LRECREAL8
C
      IOL = 2*LRECREAL8
C
      OPEN (IOTMP,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='DIRECT',
     &      RECL=IOL)
C
      WRITE (IOTMP,REC=1,ERR=300) X1,X2
      WRITE (IOTMP,REC=2,ERR=300) X2,X3
      READ (IOTMP,REC=1,ERR=300) XX1,XX2
      IF ( ABS(X1-XX1).GT.TOL .OR. ABS(X2-XX2).GT.TOL ) THEN
         WRITE (6,99002) LRECREAL8
         WRITE (6,*) 'trial C: X1  X2   written  ',X1,X2
         WRITE (6,*) 'trial C: XX1 XX2  read     ',XX1,XX2
      END IF
C
      WRITE (IOTMP,REC=1,ERR=400) X1,X2,X3
      WRITE (IOTMP,REC=2,ERR=400) X2,X3,X2
      READ (IOTMP,REC=1,ERR=400) XX1,XX2,XX3
      CLOSE (IOTMP)
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,*) 'trial C: X1  X2  X3   written  ',X1,X2,X3
         WRITE (6,*) 'trial C: XX1 XX2 XX3  read     ',XX1,XX2,XX3
      END IF
C-----------------------------------------------------------------------
 300  CONTINUE
      STOP 'in <GETLRECREAL>: record length parameter not found'
C
C-----------------------------------------------------------------------
 400  CONTINUE
      WRITE (6,99001) LRECREAL8,TEXT
      CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
99001 FORMAT (//,1X,79('*'),/,33X,'<GETLRECREAL>',/,1X,79('*'),//,10X,
     &        'record length parameter for REAL''s   LRECREAL8 =',I3,/,
     &        10X,A,/)
99002 FORMAT (//,1X,79('!'),/,10X,
     &        'warning from <GETLRECREAL8>:  trying  LRECREAL8 = ',I5)
      END
