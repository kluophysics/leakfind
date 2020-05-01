C*==aa0001.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      IMPLICIT NONE
C*--AA00013
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER I,I1,I2,J
      CHARACTER*80 LINE
      CHARACTER*20 MODUL,VAR
C
C*** End of declarations rewritten by SPAG
C
C
C
      OPEN (5,FILE='kkrmain.f')
      OPEN (7,FILE='dump_modules_insert')
C
 100  CONTINUE
      READ (5,'(A)') LINE
      IF ( INDEX(LINE,'      USE ').EQ.0 ) GOTO 100
C
C==============================================================
C
 200  CONTINUE
      IF ( INDEX(LINE,'      USE ').NE.0 ) THEN
         I1 = 11
         I2 = INDEX(LINE,'ONLY') - 2
         MODUL = LINE(I1:I2)
         WRITE (6,'(2A)') 'C      MODUL: ',MODUL
         WRITE (7,'(2A)') 'C      MODUL: ',MODUL
C
         I1 = INDEX(LINE,':') + 1
         LINE = LINE(I1:80)
      ELSE
         LINE = LINE(7:80)
      END IF
C
      J = 0
      DO I = 1,80
         IF ( LINE(I:I).NE.' ' .AND. LINE(I:I).NE.',' ) THEN
            J = J + 1
            VAR(J:J) = LINE(I:I)
         ELSE IF ( J.GT.0 ) THEN
            WRITE (6,99001) MODUL,VAR(1:J),VAR(1:J)
            WRITE (7,99001) MODUL,VAR(1:J),VAR(1:J)
            J = 0
         END IF
      END DO
C
C
      READ (5,'(A)') LINE
      IF ( INDEX(LINE,'IMPLICIT NONE').EQ.0 ) GOTO 200
      WRITE (6,*)
      WRITE (6,*) 'DUMP INSERT in dump_modules_insert'
      STOP 'finished'
C
99001 FORMAT ('      write(iotmp,*) ''',A,A,''', ',A)
      END
C
