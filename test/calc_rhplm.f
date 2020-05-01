C*==calc_rhplm.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_RHPLM(X,Y,Z,HP,LHPMAX,NLMHPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   REAL  harmonic polynomials  HP(->r) = r**l * Ylm(^r)           *
C   *                                                                  *
C   *   for  ->r  = (x,y,z)  up to l_max = LHPMAX                      *
C   *                                                                  *
C   *   NOTE:  NLMHPMAX  = (LHPMAX+1)**2   !!!                         *
C   *                                                                  *
C   *   see: Williams, Janak and Moruzzi                               *
C   *        Phys.Rev. B6 p.4509 (1972)                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_ANGMOM,ONLY:FACT
      IMPLICIT NONE
C*--CALC_RHPLM19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LHPMAX,NLMHPMAX
      REAL*8 X,Y,Z
      REAL*8 HP(NLMHPMAX)
C
C Local variables
C
      REAL*8 AUX,C(:),F,F1,F10,F2,F20,F3,Q(:,:),RSQ,S(:),T(:,:),XY,ZSQ
      INTEGER I,J,JP1,L,LMAX_INIT,M0L
      SAVE C,Q,S,T
C
C*** End of declarations rewritten by SPAG
C
      DATA LMAX_INIT/ - 1/
C
      ALLOCATABLE C,S,T,Q
C
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
      IF ( LHPMAX.GT.LMAX_INIT ) THEN
C
         IF ( LMAX_INIT.GT.-1 ) DEALLOCATE (C,S,T,Q)
C
         LMAX_INIT = LHPMAX
C
         ALLOCATE (C(0:LMAX_INIT),T(0:LMAX_INIT,0:LMAX_INIT))
         ALLOCATE (S(0:LMAX_INIT),Q(0:LMAX_INIT,0:LMAX_INIT))
C
C ------------------------ pre factors for the real harmonic polynomials
C
         DO L = 0,LMAX_INIT
            DO J = 0,L
C
               IF ( J.EQ.0 ) THEN
                  AUX = 0.5D0
               ELSE
                  AUX = FACT(L-J)/FACT(L+J)
               END IF
C
               Q(J,L) = SQRT(AUX*(2.0D0*L+1)/(2.0D0*PI))
C
            END DO
         END DO
C
C --------------- argument-independent terms of the harmonic polynomials
C
C        T(L,J)
C        T(L,L) = (2*L-1)!!
         T(0,0) = 1.0D0
         T(1,1) = 1.0D0
C
         DO L = 1,LMAX_INIT
            T(L,L) = T(L-1,L-1)*(2.0D0*L-1.0D0)
         END DO
C
         C(0) = 1.0D0
         S(0) = 0.0D0
C
      END IF
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
C
      XY = X*X + Y*Y
      ZSQ = Z*Z
      RSQ = XY + ZSQ
C
C     /WJM (B16) (B17)/
C
      DO J = 0,(LHPMAX-1)
         JP1 = J + 1
         C(JP1) = X*C(J) - Y*S(J)
         S(JP1) = X*S(J) + Y*C(J)
      END DO
C
C     T(L,J)
      T(1,0) = Z
C
C     /WJM (B19)/
C
      F1 = Z
      F2 = 0
      F3 = 1.0D0
C
      DO L = 1,(LHPMAX-1)
         F1 = F1 + Z + Z
         F2 = F2 + RSQ
         F3 = F3 + 1.0D0
         T(L+1,0) = (F1*T(L,0)-F2*T(L-1,0))/F3
      END DO
C
C
C     /WJM (B20)/
C
      IF ( XY.GT.ZSQ ) THEN
C
         F20 = Z/XY
         F10 = 1.0D0 + Z*F20
C
         DO J = 0,(LHPMAX-1)
            JP1 = J + 1
            F1 = F10*(J+J+1)
            F2 = F20
C
            DO L = (J+2),LHPMAX
               F1 = F1 + F10
               F2 = F2 + F20
               T(L,JP1) = F1*T(L-1,J) - F2*T(L,J)
            END DO
         END DO
C
      ELSE
C
         F1 = -XY/Z
         F20 = RSQ/Z
C
         DO L = 2,LHPMAX
            J = L
            F2 = F20*(L+J)
            F3 = L - J
C
            DO I = 1,(L-1)
               J = J - 1
               F2 = F2 - F20
               F3 = F3 + 1.0D0
               T(L,J) = (F1*T(L,J+1)+F2*T(L-1,J))/F3
            END DO
         END DO
C
      END IF
C
C     /WJM (B10) - B(12)/
C
C     HP(1) = Q(0,0) * T(0,0)
C
      HP(1) = Q(0,0)
C
      M0L = 1
      DO L = 1,LHPMAX
         M0L = M0L + L + L
         HP(M0L) = Q(0,L)*T(L,0)
         DO J = 1,L
C
            F = Q(J,L)*T(L,J)
            HP(M0L+J) = F*C(J)
            HP(M0L-J) = F*S(J)
C
         END DO
      END DO
C
      IF ( M0L+LHPMAX.NE.NLMHPMAX ) THEN
         WRITE (*,*) 'M0L         ',M0L
         WRITE (*,*) 'LHPMAX      ',LHPMAX
         WRITE (*,*) 'M0L+LHPMAX  ',M0L + LHPMAX
         WRITE (*,*) 'NLMHPMAX    ',NLMHPMAX
         STOP 'in <CALC_RHPLM> NLMAX**2 != NLMMAX'
      END IF
C
      END
