C*==cgc_racah.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION CGC_RACAH(J1,J2,J3,M1,M2,M3)
C   ********************************************************************
C   *                                                                  *
C   *     CLEBSCH GORDAN coefficients for integer AND half-integer     *
C   *     quantum numbers  J1,J2,J3,M1,M2,M3                           *
C   *     according to the formula of   RACAH                          *
C   *                                                                  *
C   *     see: M.E.Rose Elementary Theory of Angular Momentum          *
C   *          Eq. (3.19)                                              *
C   *          or Edmonds Eq. (3.6.11) page 45                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:FACT
      IMPLICIT NONE
C*--CGC_RACAH17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J1,J2,J3,M1,M2,M3
      REAL*8 CGC_RACAH
C
C Local variables
C
      INTEGER J,N,N1,N2,N3,N4,N5,NBOT,NTOP
      REAL*8 RFACT
      REAL*8 RSUM,S,VF,X,Y
C
C*** End of declarations rewritten by SPAG
C
C inline function    factorial for REAL argument
C
      RFACT(X) = FACT(NINT(X))
C
      CGC_RACAH = 0.0D0
C
      IF ( ABS(M3-(M1+M2)).GT.1.0D-6 ) RETURN
      IF ( ABS(J1-J2).GT.J3 ) RETURN
      IF ( (J1+J2).LT.J3 ) RETURN
      IF ( ABS(M1).GT.(J1+1.0D-6) ) RETURN
      IF ( ABS(M2).GT.(J2+1.0D-6) ) RETURN
      IF ( ABS(M3).GT.(J3+1.0D-6) ) RETURN
C
      DO J = ABS(NINT(2*(J1-J2))),NINT(2*(J1+J2)),2
         IF ( J.EQ.NINT(2*J3) ) GOTO 100
      END DO
      RETURN
C
C
 100  CONTINUE
      X = (2.0D0*J3+1.0D0)*RFACT(J1+J2-J3)*RFACT(J1-J2+J3)
     &    *RFACT(-J1+J2+J3)*RFACT(J1+M1)*RFACT(J1-M1)*RFACT(J2+M2)
     &    *RFACT(J2-M2)*RFACT(J3+M3)*RFACT(J3-M3)
C
      Y = RFACT(J1+J2+J3+1)
C
      VF = DSQRT(X/Y)
C
      N1 = NINT(J1+J2-J3)
      N2 = NINT(J1-M1)
      N3 = NINT(J2+M2)
      N4 = NINT(J3-J2+M1)
      N5 = NINT(J3-J1-M2)
      NTOP = MIN(N1,N2,N3)
      NBOT = MAX(0,-N4,-N5)
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
         Y = FACT(N)*FACT(N1-N)*FACT(N2-N)*FACT(N3-N)*FACT(N4+N)
     &       *FACT(N5+N)
         RSUM = RSUM + (S/Y)
      END DO
C
      CGC_RACAH = VF*RSUM
C
      END
C*==cgc_racah_int.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      FUNCTION CGC_RACAH_INT(L1,L2,L3,ML1,ML2,ML3)
C   ********************************************************************
C   *                                                                  *
C   *     CLEBSCH GORDAN coefficients for integer                      *
C   *     quantum numbers  L1,L2,L3,ML1,ML2,ML3J1,J2,J3,M1,M2,M3       *
C   *     according to the formula of   RACAH                          *
C   *                                                                  *
C   *     see: M.E.Rose Elementary Theory of Angular Momentum          *
C   *          Eq. (3.19)                                              *
C   *          or Edmonds Eq. (3.6.11) page 45                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CGC_RACAH_INT114
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L1,L2,L3,ML1,ML2,ML3
      REAL*8 CGC_RACAH_INT
C
C Local variables
C
      REAL*8 CGC_RACAH
      REAL*8 J1,J2,J3,M1,M2,M3
C
C*** End of declarations rewritten by SPAG
C
      J1 = DBLE(L1)
      J2 = DBLE(L2)
      J3 = DBLE(L3)
      M1 = DBLE(ML1)
      M2 = DBLE(ML2)
      M3 = DBLE(ML3)
C
      CGC_RACAH_INT = CGC_RACAH(J1,J2,J3,M1,M2,M3)
C
      END
