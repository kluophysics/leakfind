C*==viomix.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE VIOMIX(V1,V2,A,N,MESHR)
C----------------------------------------------------------------------
C
      IMPLICIT NONE
C*--VIOMIX6
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER MESHR,N
      REAL*8 A(2),V1(N,2),V2(N,2)
C
C Local variables
C
      REAL*8 B(2),VM,VP
      INTEGER K
C
C*** End of declarations rewritten by SPAG
C
C
C
      B(1) = 1D0 - A(1)
      B(2) = 1D0 - A(2)
      DO K = 1,MESHR
         VP = B(1)*(V1(K,1)+V1(K,2)) + A(1)*(V2(K,1)+V2(K,2))
         VM = B(2)*(V1(K,1)-V1(K,2)) + A(2)*(V2(K,1)-V2(K,2))
         V2(K,1) = 5D-1*(VP+VM)
         V2(K,2) = 5D-1*(VP-VM)
      END DO
      END
