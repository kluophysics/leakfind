C*==gaunt_cylm.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      FUNCTION GAUNT_CYLM(L1,M1,L2,M2,L3,M3)
C   ********************************************************************
C   *                                                                  *
C   *     GAUNT COEFFICIENTS for complex spherical harmonics  Y[l,m]   *
C   *                                                                  *
C   *            G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]              *
C   *                                                                  *
C   * see: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM  Eq. (4.34)  *
C   *                                                                  *
C   * 26/01/95  HE                                                     *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--GAUNT_CYLM16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L1,L2,L3,M1,M2,M3
      REAL*8 GAUNT_CYLM
C
C Local variables
C
      REAL*8 CGC_RACAH
      REAL*8 G
C
C*** End of declarations rewritten by SPAG
C
      IF ( (L1.LT.0) .OR. (L2.LT.0) .OR. (L3.LT.0) ) THEN
         G = 0.0D0
      ELSE
         G = (DBLE(2*L2+1)*DBLE(2*L3+1)/(4.0D0*PI*DBLE(2*L1+1)))
     &       **0.5D0*CGC_RACAH(DBLE(L3),DBLE(L2),DBLE(L1),DBLE(M3),
     &       DBLE(M2),DBLE(M1))*CGC_RACAH(DBLE(L3),DBLE(L2),DBLE(L1),
     &       0.0D0,0.0D0,0.0D0)
      END IF
C
      IF ( ABS(G).LT.1D-12 ) THEN
         GAUNT_CYLM = 0D0
      ELSE
         GAUNT_CYLM = G
      END IF
C
      END
