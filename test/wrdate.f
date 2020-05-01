C*==wrdate.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WRDATE(TEXT,IFIL)
C   ********************************************************************
C   *                                                                  *
C   *   write TEXT + date and time to chanel IFIL                      *
C   *                                                                  *
C   *   DATE_AND_TIME argument VALUES                                  *
C   *    (1) 4-digit year.                                             *
C   *    (2) month of the year.                                        *
C   *    (3) day of the month.                                         *
C   *    (4) time difference w.r.t. Coordinated Universal Time in min. *
C   *    (5) hour of the day (range 0 to 23)                           *
C   *    (6) minutes of the hour (range 0 to 59)                       *
C   *    (7) seconds of the minute (range 0 to 59)                     *
C   *    (8) milliseconds of the second (range 0 to 999)               *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,GIT_HASH,GIT_COMPILE_DATE,
     &    GIT_BRANCH
      USE MOD_CALCMODE,ONLY:TASK
      IMPLICIT NONE
C*--WRDATE22
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL
      CHARACTER*(*) TEXT
C
C Local variables
C
      CHARACTER*12 DATE,TIME,ZONE
      INTEGER LTXT,VALUES(8)
C
C*** End of declarations rewritten by SPAG
C
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
C
      LTXT = LEN_TRIM(TEXT)
C
      IF ( IFIL.NE.IFILBUILDBOT ) THEN
C
         IF ( LTXT.GT.0 ) THEN
            WRITE (IFIL,99001) TEXT(1:LTXT),VALUES(3),VALUES(2),
     &                         VALUES(1),VALUES(5),VALUES(6),VALUES(7)
         ELSE
            WRITE (IFIL,99002) VALUES(3),VALUES(2),VALUES(1),VALUES(5),
     &                         VALUES(6),VALUES(7)
         END IF
C
         IF ( IPRINT.GE.10 ) WRITE (6,*) 'DATE_AND_TIME',DATE,TIME,ZONE
C
      ELSE
C
         WRITE (IFILBUILDBOT,99003) TEXT(1:LTXT),VALUES(3),VALUES(2),
     &                              VALUES(1),VALUES(5),VALUES(6),
     &                              VALUES(7),TASK,
     &                              GIT_BRANCH(1:LEN_TRIM(GIT_BRANCH)),
     &                              GIT_HASH(1:LEN_TRIM(GIT_HASH)),
     &                              GIT_COMPILE_DATE
     &                              (1:LEN_TRIM(GIT_COMPILE_DATE))
C
      END IF
C
99001 FORMAT (10X,A,'  on  ',I4.2,'/',I2.2,'/',I4.4,'    at   ',I2.2,
     &        ':',I2.2,':',I2.2,/)
99002 FORMAT (10X,'on  ',I4.2,'/',I2.2,'/',I2.2,'    at   ',I2.2,':',
     &        I2.2,':',I2.2,/)
99003 FORMAT ('#',/,'#  run of ',A,'  on  ',I4.2,'/',I2.2,'/',I4.4,
     &        '    at   ',I2.2,':',I2.2,':',I2.2,/,'#',/,
     &        '#  TASK             ',A,/,'#  GIT_BRANCH       ',A,/,
     &        '#  GIT_HASH         ',A,/,'#  GIT_COMPILE_DATE ',A,/,'#')
C
      END
