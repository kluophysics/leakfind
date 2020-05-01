C*==csumsslt.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      FUNCTION CSUMSSLT(X,IL,IT,I,J,NLMAX,NTMAX,NS1,N,NS2,M)
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CSUMSSLT5
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IL,IT,J,M,N,NLMAX,NS1,NS2,NTMAX
      COMPLEX*16 CSUMSSLT
      COMPLEX*16 X(NLMAX,NTMAX,2,N,2,M)
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER K,L
C
C*** End of declarations rewritten by SPAG
C
      CSUM = C0
      DO K = 1,NS1
         DO L = 1,NS2
            CSUM = CSUM + X(IL,IT,K,I,L,J)
         END DO
      END DO
      CSUMSSLT = CSUM
      END
