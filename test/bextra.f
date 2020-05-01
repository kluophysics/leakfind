C*==bextra.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION BEXTRA(R,INTR0R,INTR0S)
C   ********************************************************************
C   *                                                                  *
C   * extrapolate integral  INT(R..S) = INT(R0..S) - INT(R0..R) to R=0 *
C   * linear and quadratic extrapolation give more or less the same    *
C   * results for potentials of akai's program                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--BEXTRA13
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 INTR0S
      REAL*8 BEXTRA
      COMPLEX*16 INTR0R(NRMAX)
      REAL*8 R(NRMAX)
C
C Local variables
C
      REAL*8 A2,R1,R2,R3,Y1,Y10,Y2,Y3
      INTEGER I1,I2,I3
C
C*** End of declarations rewritten by SPAG
C
      I1 = 1
      I2 = 7
      I3 = 14
C
      Y10 = DBLE(INTR0S-INTR0R(I1))
      Y1 = 1.0D0
      Y2 = DBLE(INTR0S-INTR0R(I2))/Y10
      Y3 = DBLE(INTR0S-INTR0R(I3))/Y10
      R1 = 1.0D0
      R2 = R(I2)/R(I1)
      R3 = R(I3)/R(I1)
C
      A2 = (Y1*R2*R3**2+R1*R2**2*Y3+R1**2*Y2*R3-Y1*R2**2*R3-R1*Y2*R3**2-
     &     R1**2*R2*Y3)
     &     /(R2*R3**2+R1*R2**2+R1**2*R3-R2**2*R3-R1*R3**2-R1**2*R2)
C
      BEXTRA = A2*Y10
C
      END
C*==cbextra.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION CBEXTRA(R,INTR0R,INTR0S)
C   ********************************************************************
C   *                                                                  *
C   * extrapolate integral  INT(R..S) = INT(R0..S) - INT(R0..R) to R=0 *
C   * linear and quadratic extrapolation give more or less the same    *
C   * results for potentials of akai's program                         *
C   *                                                                  *
C   * COMPLEX VERSION - to calculate hyperfine fields                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CBEXTRA75
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 INTR0S
      COMPLEX*16 CBEXTRA
      COMPLEX*16 INTR0R(NRMAX)
      REAL*8 R(NRMAX)
C
C Local variables
C
      COMPLEX*16 A2,Y1,Y10,Y2,Y3
      INTEGER I1,I2,I3
      REAL*8 R1,R2,R3
C
C*** End of declarations rewritten by SPAG
C
      I1 = 1
      I2 = 7
      I3 = 14
C
      Y10 = INTR0S - INTR0R(I1)
      Y1 = 1.0D0
      Y2 = (INTR0S-INTR0R(I2))/Y10
      Y3 = (INTR0S-INTR0R(I3))/Y10
      R1 = 1.0D0
      R2 = R(I2)/R(I1)
      R3 = R(I3)/R(I1)
C
      A2 = (Y1*R2*R3**2+R1*R2**2*Y3+R1**2*Y2*R3-Y1*R2**2*R3-R1*Y2*R3**2-
     &     R1**2*R2*Y3)
     &     /(R2*R3**2+R1*R2**2+R1**2*R3-R2**2*R3-R1*R3**2-R1**2*R2)
C
      CBEXTRA = A2*Y10
C
      END
