C*==open_iotmp_scratch.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C   ********************************************************************
C   *                                                                  *
C   *  open the temporary UNFORMATTED scratch file  IOTMP              *
C   *                                                                  *
C   *  check first whether the file is already open - if yes: STOP     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--OPEN_IOTMP_SCRATCH12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IOTMP
      CHARACTER*(*) ROUTINE
C
C Local variables
C
      LOGICAL FILE_IS_OPEN
C
C*** End of declarations rewritten by SPAG
C
C-------------------------------------------- check whether file is open
C
      INQUIRE (IOTMP,OPENED=FILE_IS_OPEN)
C
      IF ( FILE_IS_OPEN ) THEN
C
         WRITE (6,99001) IOTMP
         CALL STOP_MESSAGE(ROUTINE,'scratch file already open')
C
      END IF
C
C----------------------------------------- open UNFORMATTED scratch file
C
      OPEN (UNIT=IOTMP,STATUS='SCRATCH',FORM='UNFORMATTED')
C
99001 FORMAT (//,10X,'<OPEN_IOTMP_SCRATCH> called for IOTMP =',I5,//)
      END
C*==open_iotmp_file.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM)
C   ********************************************************************
C   *                                                                  *
C   *  open the temporary FORMATTED file  IOTMP  with name  FILNAM     *
C   *                                                                  *
C   *  check first whether the file is already open - if yes: STOP     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--OPEN_IOTMP_FILE66
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) FILNAM,ROUTINE
      INTEGER IOTMP
C
C Local variables
C
      LOGICAL FILE_IS_OPEN
C
C*** End of declarations rewritten by SPAG
C
C-------------------------------------------- check whether file is open
C
      INQUIRE (IOTMP,OPENED=FILE_IS_OPEN)
C
      IF ( FILE_IS_OPEN ) THEN
C
         WRITE (6,99001) IOTMP
         CALL STOP_MESSAGE(ROUTINE,'scratch file already open')
C
      END IF
C
C--------------------------------------------- open named FORMATTED file
C
      OPEN (UNIT=IOTMP,FILE=FILNAM)
C
99001 FORMAT (//,10X,'<OPEN_IOTMP_SCRATCH> called for IOTMP =',I5,//)
      END
