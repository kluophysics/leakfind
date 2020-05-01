C*==rmserr.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RMSERR(V1,V2,RMS,DR,XR,MESHR)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*--RMSERR5
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER MESHR
      REAL*8 RMS
      REAL*8 DR(MESHR),V1(MESHR),V2(MESHR),XR(MESHR)
C
C Local variables
C
      REAL*8 ATVOL,VOLINS,VOLINT
      INTEGER INTS,K
C
C*** End of declarations rewritten by SPAG
C
C
C
      ATVOL = XR(MESHR)**3/3D0
      VOLINS = XR(MESHR-1)**3/3D0
      VOLINT = ATVOL - VOLINS
      RMS = 0D0
      INTS = MESHR - 2
      DO K = 2,INTS,2
         RMS = RMS + ((V1(K)-V2(K))*XR(K))**2*DR(K)
      END DO
      RMS = RMS*2D0
      DO K = 3,INTS,2
         RMS = RMS + ((V1(K)-V2(K))*XR(K))**2*DR(K)
      END DO
      RMS = RMS*2D0 + ((V1(1)-V2(1))*XR(1))**2*DR(1)
     &      + ((V1(MESHR-1)-V2(MESHR-1))*XR(MESHR-1))**2*DR(MESHR-1)
      RMS = RMS/3D0 + VOLINT*(V1(MESHR)-V2(MESHR))**2
      RMS = SQRT(RMS*ATVOL)
      END
