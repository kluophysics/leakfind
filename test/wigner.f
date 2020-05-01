C*==wig_3j.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      FUNCTION WIG_3J(J1,J2,J3,M1,M2,M3)
C   ********************************************************************
C   *                                                                  *
C   *     Wigner 3j-Symbol for arbitrary quantum numbers j1,j2,...     *
C   *                                                                  *
C   *     SEE: EDMONDS EQ. (3.7.3)  PAGE 46                            *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--WIG_3J11
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J1,J2,J3,M1,M2,M3
      REAL*8 WIG_3J
C
C Local variables
C
      REAL*8 CGC_RACAH
C
C*** End of declarations rewritten by SPAG
C
C
      WIG_3J = (-1)**NINT(J1-J2-M3)*(1D0/SQRT(2D0*J3+1D0))
     &         *CGC_RACAH(J1,J2,J3,M1,M2,-M3)
C
      END
C*==wig_6j.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      FUNCTION WIG_6J(J1,J2,J12,J3,J,J23)
C   ********************************************************************
C   *                                                                  *
C   *     Wigner 6j-Symbol for arbitrary quantum numbers j1,j2,...     *
C   *                                                                  *
C   *     SEE: EDMONDS EQ. (6.1.5)  PAGE 92                            *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--WIG_6J52
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J,J1,J12,J2,J23,J3
      REAL*8 WIG_6J
C
C Local variables
C
      REAL*8 CGC_RACAH
      REAL*8 M,M1,M12,M2,M23,M3,RSUM
      INTEGER M1M05,M2M05,M3M05
C
C*** End of declarations rewritten by SPAG
C
      RSUM = 0.0D0
C
      DO M1M05 = NINT(-J1-0.5D0),NINT(J1-0.5D0)
         M1 = M1M05 + 0.5D0
C
         DO M2M05 = NINT(-J2-0.5D0),NINT(J2-0.5D0)
            M2 = M2M05 + 0.5D0
C
            DO M3M05 = NINT(-J3-0.5D0),NINT(J3-0.5D0)
               M3 = M3M05 + 0.5D0
C
               M = M1 + M2 + M3
               M12 = M1 + M2
               M23 = M - M1
C
               RSUM = RSUM + CGC_RACAH(J1,J2,J12,M1,M2,M12)
     &                *CGC_RACAH(J12,J3,J,M12,M3,M)
     &                *CGC_RACAH(J2,J3,J23,M2,M3,M23)
     &                *CGC_RACAH(J1,J23,J,M1,M23,M)
            END DO
         END DO
      END DO
C
      WIG_6J = RSUM*(-1)**NINT(J1+J2+J3+J)
     &         /SQRT((2D0*J12+1D0)*(2D0*J23+1D0))/(2D0*J+1)
C
      END
C*==wig_6j_racah.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      FUNCTION WIG_6J_RACAH(J1,J2,J3,L1,L2,L3)
C   ********************************************************************
C   *                                                                  *
C   *     Wigner 6j-Symbol for arbitrary quantum numbers j1,j2,...     *
C   *     according to the general expression due to Racah             *
C   *                                                                  *
C   *     SEE: EDMONDS EQ. (6.3.7)  PAGE 99                            *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:FACT
      IMPLICIT NONE
C*--WIG_6J_RACAH119
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J1,J2,J3,L1,L2,L3
      REAL*8 WIG_6J_RACAH
C
C Local variables
C
      REAL*8 A,B,C,RSUM,S,X
      REAL*8 DELTA,RFACT
      INTEGER N,N1,N2,N3,N4,N5,N6,N7,NBOT,NTOP
C
C*** End of declarations rewritten by SPAG
C
C INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
C
      RFACT(X) = FACT(NINT(X))
      DELTA(A,B,C) = SQRT(RFACT(A+B-C)*RFACT(A-B+C)*RFACT(-A+B+C)
     &               /RFACT(A+B+C+1))
C
      N1 = NINT(J1+J2+J3)
      N2 = NINT(J1+L2+L3)
      N3 = NINT(L1+J2+L3)
      N4 = NINT(L1+L2+J3)
      N5 = NINT(J1+J2+L1+L2)
      N6 = NINT(J2+J3+L2+L3)
      N7 = NINT(J3+J1+L3+L1)
      NBOT = MAX(N1,N2,N3,N4)
      NTOP = MIN(N5,N6,N7)
C
      N = NBOT + 1
      IF ( N.EQ.(2*(N/2)) ) THEN
         S = +1.0D0
      ELSE
         S = -1.0D0
      END IF
      RSUM = 0.0D0
C
      DO N = NBOT,NTOP
         S = -S
         RSUM = RSUM + S*FACT(N+1)
     &          /(FACT(N-N1)*FACT(N-N2)*FACT(N-N3)*FACT(N-N4)*FACT(N5-N)
     &          *FACT(N6-N)*FACT(N7-N))
      END DO
C
      WIG_6J_RACAH = RSUM*DELTA(J1,J2,J3)*DELTA(J1,L2,L3)
     &               *DELTA(L1,J2,L3)*DELTA(L1,L2,J3)
C
      END
