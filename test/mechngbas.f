C*==mechngbas.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MECHNGBAS(MODE,ME,NKM,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *      change the basis for the matrix elements  ME                *
C   *                                                                  *
C   *      MODE = 'SPH>CAR' :     -1, 0, +1  =>  x, y, z               *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MECHNGBAS11
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI
      PARAMETER (CI=(0.0D0,1.0D0))
C
C Dummy arguments
C
      CHARACTER*7 MODE
      INTEGER NKM,NKMMAX
      COMPLEX*16 ME(NKMMAX,NKMMAX,3)
C
C Local variables
C
      COMPLEX*16 CIOVRT2,ME0,MEMIN,MEPLS
      INTEGER I,J
      REAL*8 ONEOVRT2
C
C*** End of declarations rewritten by SPAG
C
      ONEOVRT2 = 1D0/DSQRT(2.0D0)
      CIOVRT2 = CI*ONEOVRT2
C
C------------------------------ change from spherical to cartesian basis
C------------------------------             -1, 0, +1 =>  x, y, z
C
      IF ( MODE.EQ.'SPH>CAR' ) THEN
C
         DO I = 1,NKM
            DO J = 1,NKM
C
               MEMIN = ME(I,J,1)
               ME0 = ME(I,J,2)
               MEPLS = ME(I,J,3)
C
               ME(I,J,1) = (MEMIN-MEPLS)*ONEOVRT2
               ME(I,J,2) = (MEMIN+MEPLS)*CIOVRT2
               ME(I,J,3) = ME0
C
            END DO
         END DO
C
      ELSE
         STOP 'in <MECHNGBAS>   MODE not allowed'
      END IF
C
      END
