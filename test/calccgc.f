C*==calccgc.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCCGC(CGC,NKMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
C   *                                                                  *
C   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
C   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
C   *   IS= 1/2  SPIN DOWN/UP                                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CALCCGC14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMAX,NKMPMAX
      REAL*8 CGC(NKMPMAX,2)
C
C Local variables
C
      INTEGER IKM,K,KAPPA,L,M,NL,NMUE
      REAL*8 J,MUE,TWOLP1
C
C*** End of declarations rewritten by SPAG
C
      NL = (NKMAX+1)/2
      IF ( NKMPMAX.NE.2*NL**2+2*NL ) THEN
         WRITE (6,99001) NKMAX,NKMPMAX,2*NL**2 + 2*NL
         STOP
      END IF
C
      IKM = 0
      DO K = 1,(NKMAX+1)
         L = K/2
         IF ( K.EQ.2*L ) THEN
            KAPPA = L
         ELSE
            KAPPA = -L - 1
         END IF
         J = ABS(KAPPA) - 0.5D0
         NMUE = NINT(2*J+1)
         MUE = -J - 1.0D0
         TWOLP1 = 2.0D0*L + 1.0D0
C
         IF ( KAPPA.LT.0 ) THEN
C
C     J = L + 1/2
            DO M = 1,NMUE
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L-MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = DSQRT((L+MUE+0.5D0)/TWOLP1)
            END DO
         ELSE
C     J = L - 1/2
            DO M = 1,NMUE
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L+MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = -DSQRT((L-MUE+0.5D0)/TWOLP1)
C
            END DO
         END IF
C
      END DO
99001 FORMAT (5X,'calling <CALCCGC> using inconsistent arguments',/,5X,
     &        'NKMAX    =',I5,/,5X,'NKMPMAX  =',I5,'  instead of: ',I5)
C
      END
