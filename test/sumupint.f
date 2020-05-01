C*==sumupint.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SUMUPINT(CSUM,VG,G,WG,VF,F,WF,N)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SUMUPINT8
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CSUM
      INTEGER N
      REAL*8 VF,VG
      COMPLEX*16 F(2,2),G(2,2)
      REAL*8 WF(2,2),WG(2,2)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
      CSUM = 0.0D0
      DO J = 1,N
         DO I = 1,N
            CSUM = CSUM + VG*G(I,J)*WG(I,J) + VF*F(I,J)*WF(I,J)
         END DO
      END DO
C
      END
