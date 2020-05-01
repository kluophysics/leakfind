C*==calc_chplm.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALC_CHPLM(X,Y,Z,HP,LHPMAX,NLMHPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   COMPLEX  harmonic polynomials  HP(->r) = r**l * Ylm(^r)        *
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
C*--CALC_CHPLM19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LHPMAX,NLMHPMAX
      REAL*8 X,Y,Z
      COMPLEX*16 HP(NLMHPMAX)
C
C Local variables
C
      REAL*8 AUX,C(:),F1,F10,F2,F20,F3,Q(:,:),RSQ,S(:),T(:,:),WZ05,XY,
     &       ZSQ
      INTEGER I,J,JP1,L,LMAX_INIT,M0L
      SAVE C,Q,S,T
C
C*** End of declarations rewritten by SPAG
C
      DATA LMAX_INIT/ - 1/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C,S,T,Q
C
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
      IF ( LHPMAX.GT.LMAX_INIT ) THEN
         IF ( LMAX_INIT.GT.-1 ) DEALLOCATE (C,S,T,Q)
         LMAX_INIT = LHPMAX
         ALLOCATE (C(0:LMAX_INIT),T(0:LMAX_INIT,0:LMAX_INIT))
         ALLOCATE (S(0:LMAX_INIT),Q(-LMAX_INIT:LMAX_INIT,0:LMAX_INIT))
C
C ---------------------- pre factors for the complex spherical harmonics
C
         WZ05 = SQRT(1.0D0/2.0D0)
         DO L = 0,LMAX_INIT
            Q(0,L) = SQRT(0.5D0*(2.0D0*L+1)/(2.0D0*PI))
            DO J = 1,L
               AUX = FACT(L-J)/FACT(L+J)
               Q(-J,L) = SQRT(AUX*(2.0D0*L+1)/(2.0D0*PI))*WZ05
               Q(+J,L) = Q(-J,L)*(-1.0D0)**J
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
C
         DO J = 1,L
            HP(M0L+J) = Q(J,L)*T(L,J)*DCMPLX(C(J),S(J))
            HP(M0L-J) = Q(-J,L)*T(L,J)*DCMPLX(C(J),-S(J))
         END DO
      END DO
C
      IF ( M0L+LHPMAX.NE.NLMHPMAX )
     &      STOP 'in <CALC_CHPLM> NLMAX**2 != NLMMAX'
C
      END
