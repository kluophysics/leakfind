C*==copfil.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COPFIL(N1,N2,LINE,STR5,M)
      IMPLICIT NONE
C*--COPFIL4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*80 LINE
      INTEGER M,N1,N2
      CHARACTER*5 STR5
C
C Local variables
C
      INTEGER I,NLINES
C
C*** End of declarations rewritten by SPAG
C
C   ********************************************************************
C   *                                                                  *
C   *  copy file  N1  to  N2  from actual position on                  *
C   *  if M > 0 :             copy just  M  lines                      *
C   *  else                                                            *
C   *  if STR5 =  '<EOF>':    copy until  EOF  is reached              *
C   *                         and rewind files                         *
C   *  if STR5 <> '<EOF>':    copy until a line starting with  STR5    *
C   *                         is found                                 *
C   *                                                                  *
C   ********************************************************************
C
      IF ( M.GT.0 ) THEN
         NLINES = M
      ELSE
         NLINES = 10000
      END IF
C
      IF ( STR5.EQ.'<EOF>' .OR. M.GT.0 ) THEN
         DO I = 1,NLINES
            READ (N1,'(A)',END=100) LINE
            WRITE (N2,'(A)') LINE
         END DO
C
      ELSE
C
         DO I = 1,10000
            READ (N1,'(A)',END=100) LINE
            WRITE (N2,'(A)') LINE
            IF ( LINE(1:5).EQ.STR5 ) RETURN
         END DO
      END IF
C
 100  CONTINUE
      IF ( STR5.EQ.'<EOF>' .AND. M.LE.0 ) THEN
         REWIND N1
         REWIND N2
      END IF
      END
