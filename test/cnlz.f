C*==cnlz.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      FUNCTION CNLZ(L,Z)
C   ********************************************************************
C   *                                                                  *
C   *     von NEUMANN  - FUNCTION  N(L,Z)  FOR COMPLEX ARGUMENT  Z     *
C   *                  see:  e.g. MERZBACHER EQ. (10.34)               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CNLZ11
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
      INTEGER LP2MAX
      PARAMETER (LP2MAX=25)
C
C Dummy arguments
C
      INTEGER L
      COMPLEX*16 Z
      COMPLEX*16 CNLZ
C
C Local variables
C
      COMPLEX*16 CJLZ
      REAL*8 DFAC
      COMPLEX*16 DT,S(LP2MAX),T,ZSQ
      INTEGER I,K,LLP1
C
C*** End of declarations rewritten by SPAG
C
C
C
      IF ( L.LT.0 ) THEN
         CNLZ = CJLZ(L+1,Z)
         RETURN
      END IF
C
      ZSQ = Z*Z
      LLP1 = L + L + 1
C
      IF ( ABS(ZSQ/DBLE(LLP1)).LE.10.D0 ) THEN
C
         DFAC = 1.0D0
         DO K = 3,LLP1,2
            DFAC = DFAC*DBLE(K)
         END DO
C
         DT = C1
         T = C1
C
         DO I = 2,400,2
            DT = -DT*ZSQ/DBLE(I*(I-LLP1))
            T = T + DT
            IF ( ABS(DT).LT.1.0D-10 ) EXIT
         END DO
C
         CNLZ = -T*Z**(-L-1)*DFAC/DBLE(LLP1)
C
      ELSE
         IF ( L.GT.23 ) STOP '<cnlz>: l too large'
         S(2) = COS(Z)
         IF ( L.LE.0 ) THEN
            CNLZ = -S(2)*Z**(-L-1)
            RETURN
         END IF
C
         S(1) = -SIN(Z)/Z
         DO I = 3,L + 2
            S(I) = S(I-1)*(2*I-5) - ZSQ*S(I-2)
         END DO
         CNLZ = -S(L+2)*Z**(-L-1)
C
      END IF
C
      END
