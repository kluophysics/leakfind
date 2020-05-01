C*==dmft_close_iotmp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_CLOSE_IOTMP(IOTMP)
      IMPLICIT NONE
C*--DMFT_CLOSE_IOTMP4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IOTMP
C
C Local variables
C
      LOGICAL ISOPENED
C
C*** End of declarations rewritten by SPAG
C
C
C     IF UNIT IOTMP IS OPENED CLOSE IT
C
C
      INQUIRE (UNIT=IOTMP,EXIST=ISOPENED)
      IF ( ISOPENED ) CLOSE (IOTMP,ERR=99999)
C
C
99999 CONTINUE
      END
