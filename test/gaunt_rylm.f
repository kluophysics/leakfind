C*==gaunt_rylm.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      FUNCTION GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C   ********************************************************************
C   *                                                                  *
C   *    GAUNT COEFFICIENTS for   REAL   spherical harmonics  Y[l,m]   *
C   *                                                                  *
C   *            G = INT dr^  Y[l1,m1]  Y[l2,m2]  Y[l3,m3]             *
C   *                                                                  *
C   *     the selection rules are exploited                            *
C   *                                                                  *
C   *                   |l1-l2| <= l3 <= l1+l2                         *
C   *            |m1+m2| = |m3|  or  |m1-m2| = |m3|                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI,CI
      IMPLICIT NONE
C*--GAUNT_RYLM18
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L1,L2,L3,M1,M2,M3
      REAL*8 GAUNT_RYLM
C
C Local variables
C
      COMPLEX*16 CC1,CC2,CC3
      REAL*8 CGAUNT,G,WZ05,XL1,XL2,XL3,XM1,XM2,XM3
      REAL*8 CGC_RACAH
      INTEGER I1,I2,I3
C
C*** End of declarations rewritten by SPAG
C
      WZ05 = SQRT(0.5D0)
      G = 0.0D0
C
      IF ( (ABS(L1-L2).LE.L3) .AND. (L3.LE.(L1+L2)) ) THEN
         IF ( (ABS(M1+M2).EQ.ABS(M3)) .OR. (ABS(M1-M2).EQ.ABS(M3)) )
     &        THEN
C-----------------------------------------------------------------------
            XL1 = L1
            XL2 = L2
            XL3 = L3
C
            G = 0.0D0
            DO I1 = 1,2
               IF ( I1.EQ.1 ) THEN
                  XM1 = -ABS(M1)
                  IF ( M1.LT.0 ) CC1 = CI*WZ05
                  IF ( M1.EQ.0 ) CC1 = 1.0D0
                  IF ( M1.GT.0 ) CC1 = WZ05
               ELSE
                  XM1 = ABS(M1)
                  IF ( M1.LT.0 ) CC1 = -CI*WZ05*(-1)**ABS(M1)
                  IF ( M1.EQ.0 ) CC1 = 0.0D0
                  IF ( M1.GT.0 ) CC1 = WZ05*(-1)**ABS(M1)
               END IF
C
               DO I2 = 1,2
                  IF ( I2.EQ.1 ) THEN
                     XM2 = -ABS(M2)
                     IF ( M2.LT.0 ) CC2 = CI*WZ05
                     IF ( M2.EQ.0 ) CC2 = 1.0D0
                     IF ( M2.GT.0 ) CC2 = WZ05
                  ELSE
                     XM2 = ABS(M2)
                     IF ( M2.LT.0 ) CC2 = -CI*WZ05*(-1)**ABS(M2)
                     IF ( M2.EQ.0 ) CC2 = 0.0D0
                     IF ( M2.GT.0 ) CC2 = WZ05*(-1)**ABS(M2)
                  END IF
C
                  DO I3 = 1,2
                     IF ( I3.EQ.1 ) THEN
                        XM3 = -ABS(M3)
                        IF ( M3.LT.0 ) CC3 = CI*WZ05
                        IF ( M3.EQ.0 ) CC3 = 1.0D0
                        IF ( M3.GT.0 ) CC3 = WZ05
                     ELSE
                        XM3 = ABS(M3)
                        IF ( M3.LT.0 ) CC3 = -CI*WZ05*(-1)**ABS(M3)
                        IF ( M3.EQ.0 ) CC3 = 0.0D0
                        IF ( M3.GT.0 ) CC3 = WZ05*(-1)**ABS(M3)
                     END IF
C
                     IF ( (NINT(XM1)+NINT(XM2)+NINT(XM3)).EQ.0 ) THEN
C
                        CGAUNT = SQRT((2*XL1+1)*(2*XL2+1)
     &                           /(4*PI*(2*XL3+1)))*(-1)**(-NINT(XM3))
     &                           *CGC_RACAH(XL1,XL2,XL3,0.0D0,0.0D0,
     &                           0.0D0)
     &                           *CGC_RACAH(XL1,XL2,XL3,XM1,XM2,-XM3)
C
                        G = G + DREAL(CC1*CC2*CC3)*CGAUNT
                     END IF
                  END DO
               END DO
            END DO
C-----------------------------------------------------------------------
         END IF
      END IF
C
      IF ( ABS(G).LT.1D-12 ) THEN
         GAUNT_RYLM = 0D0
      ELSE
         GAUNT_RYLM = G
      END IF
      END
