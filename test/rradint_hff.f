C*==rradint_hff.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RRADINT_HFF(IM,IRTOP,RAUX1,RME)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,R,RAUX2
      IMPLICIT NONE
C*--RRADINT_HFF12
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='RRADINT_HFF')
C
C Dummy arguments
C
      INTEGER IM,IRTOP
      REAL*8 RME
      REAL*8 RAUX1(NRMAX)
C
C Local variables
C
      REAL*8 DELTA,RLIM1,RLIM2,SX,SXX,SXY,SY,Y
      INTEGER IR,IRMAX,IRMIN,N
C
C*** End of declarations rewritten by SPAG
C
      CALL RRADINT_R(IM,RAUX1,RAUX2)
C
      RLIM1 = 0.5D-5
      RLIM2 = 0.5D-4
C     RLIM1 = 1D-4
C     RLIM2 = 5D-4
      IF ( R(1,IM).GT.RLIM1 ) THEN
         RLIM1 = R(1,IM)
         RLIM2 = R(20,IM)
      END IF
C
      IRMIN = 0
      IRMAX = 0
      DO IR = 1,IRTOP
         IF ( R(IR,IM).LE.RLIM1 ) IRMIN = IR
         IF ( R(IR,IM).GE.RLIM2 ) THEN
            IRMAX = IR
            EXIT
         END IF
      END DO
      IF ( IRMIN.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'IRMIN = 0')
      IF ( IRMIN.GT.IRTOP ) CALL STOP_MESSAGE(ROUTINE,'IRMIN > IRTOP')
      IF ( IRMAX.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'IRMAX = 0')
      IF ( IRMAX.GT.IRTOP ) CALL STOP_MESSAGE(ROUTINE,'IRMAX > IRTOP')
C
      N = 0
      SX = 0.0D0
      SXX = 0.0D0
      SY = 0.0D0
      SXY = 0.0D0
C
      DO IR = IRMIN,IRMAX
         Y = RAUX2(IRTOP) - RAUX2(IR)
         N = N + 1
         SX = SX + R(IR,IM)
         SXX = SXX + R(IR,IM)**2
         SY = SY + Y
         SXY = SXY + Y*R(IR,IM)
      END DO
C
      DELTA = N*SXX - SX*SX
C
      RME = (SXX*SY-SX*SXY)/DELTA
C
      END
