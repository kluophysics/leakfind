C*==beshan.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE BESHAN(HL,JL,NL,Z,LMAX)
C-----------------------------------------------------------------------
C  calculates spherical bessel, hankel and neumann functions
C  for the orders l .le. lmax.
C  For |z| .lt. 1 the taylor expansions of jl and nl are used.
C  For |z| .ge. 1 the explicit expressions for hl(+), hl(-) are used.
C
C                            R. Zeller   Jan. 1990
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*--BESHAN12
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE COMPLEX CI
      PARAMETER (CI=(0.0D0,1.0D0))
C
C Dummy arguments
C
      INTEGER LMAX
      DOUBLE COMPLEX Z
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C
C Local variables
C
      INTEGER L,M,N
      DOUBLE PRECISION RL,RN,RNM
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
C
C*** End of declarations rewritten by SPAG
C
      ZJ = 1.D0
      ZN = 1.D0
      Z2 = Z*Z
C-----------------------------------------------------------------------
      IF ( ABS(Z).LT.LMAX+1.D0 ) THEN
         DO L = 0,LMAX
            RL = L + L
            TERMJ = -0.5D0/(RL+3.D0)*Z2
            TERMN = 0.5D0/(RL-1.D0)*Z2
            JL(L) = 1.D0
            NL(L) = 1.D0
            DO N = 2,25
               JL(L) = JL(L) + TERMJ
               NL(L) = NL(L) + TERMN
               RN = N + N
               TERMJ = -TERMJ/(RL+RN+1.D0)/RN*Z2
               TERMN = TERMN/(RL-RN+1.D0)/RN*Z2
            END DO
            JL(L) = JL(L)*ZJ
            NL(L) = -NL(L)*ZN/Z
            HL(L) = JL(L) + NL(L)*CI
C
            ZJ = ZJ*Z/(RL+3.D0)
            ZN = ZN/Z*(RL+1.D0)
         END DO
C
C-----------------------------------------------------------------------
      ELSE
C
         DO L = 0,LMAX
            HL(L) = 0.D0
            NL(L) = 0.D0
            RNM = 1.D0
            DO M = 0,L
               HL(L) = HL(L) + RNM/(-CI*(Z+Z))**M
               NL(L) = NL(L) + RNM/(CI*(Z+Z))**M
               RNM = RNM*(L*L+L-M*M-M)/(M+1.D0)
            END DO
            HL(L) = HL(L)*(-CI)**L*EXP(CI*Z)/(CI*Z)
            NL(L) = NL(L)*CI**L*EXP(-CI*Z)/(-CI*Z)
            JL(L) = (HL(L)+NL(L))*0.5D0
            NL(L) = (HL(L)-JL(L))/CI
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      END
