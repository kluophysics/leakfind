C*==strgaunt.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE STRGAUNT(LMAX,ALAT,RGNT,NRGNT,IRGNT,CIPWL,NLMAX,IG123,
     &                    LRGNT12,LRGNT123)
C   ********************************************************************
C   *                                                                  *
C   *    CALCULATION OF GAUNT COEFFICIENTS                             *
C   *                                                                  *
C   *    G(L1,M1;L2,M2;L3,M3) == INT Y(L1,M1) * Y(L2,M2) * Y(L3,M3)    *
C   *                                                                  *
C   *                 for  REAL SPHERICAL HARMONICS                    *
C   *                                                                  *
C   *    SEE EDMONDS EQS. (4.6.3), (3.7.3) AND (3.6.11)                *
C   *                                                                  *
C   *    TO SET UP THE GLL' - MATRIX FROM THE DLM'S THE GAUNTS ARE     *
C   *    MULTIPLIED BY THE FACTOR 4*PI                                 *
C   *    THE ADDITIONAL FACTOR 2*PI/A CONVERTS THE GLL' - MATRIX       *
C   *    FROM D.U.'S TO A.U.'S                                         *
C   *                                                                  *
C   *    OUTPUT:  RGNT = (2*pi/a) * 4*pi * C(l1,m1,l2,m2,l3,m3)        *
C   *                                                                  *
C   *                                                                  *
C   *    THE GAUNTS ARE CALCULATED ONLY FOR  FOR LM2 <= LM1            *
C   *    FOR THE LOOP LM1=1,LM1MAX / LM2=1,LM1 / LM3=1,LM3MAX          *
C   *    NRGNT     SPECIFIES THE NUMBER OF NON-ZERO RGNT'S             *
C   *              FOR THE CURRENT LINEAR INDEX I == (LM1,LM2)         *
C   *    IRGNT(I)  GIVES THE  LM3  FOR THE I-TH NON-ZERO RGNT          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI,CI
      IMPLICIT NONE
C*--STRGAUNT32
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER IG123,LMAX,LRGNT12,LRGNT123,NLMAX
      COMPLEX*16 CIPWL((2*NLMAX)**2)
      INTEGER IRGNT(LRGNT123),NRGNT(LRGNT12)
      REAL*8 RGNT(LRGNT123)
C
C Local variables
C
      COMPLEX*16 CC1,CC2,CC3
      REAL*8 CGAUNT,PRE,RGAUNT,WZ05,XL1,XL2,XL3,XM1,XM2,XM3
      REAL*8 CGC_RACAH
      INTEGER I1,I2,I3,L,L1,L2,L3,LM,LM1,LM1LM2,LM2,M,M1,M2,M3
C
C*** End of declarations rewritten by SPAG
C
      LM = 0
      DO L = 0,2*LMAX
         DO M = -L,L
            LM = LM + 1
            CIPWL(LM) = CI**L
         END DO
      END DO
C
      WZ05 = SQRT(0.5D0)
      PRE = 4*PI*2*PI/ALAT
      LM1LM2 = 0
      IG123 = 0
C
C=======================================================================
C
      DO L1 = 0,LMAX
         DO M1 = -L1,L1
            XL1 = L1
            LM1 = L1*(L1+1) + M1 + 1
C
            DO L2 = 0,LMAX
               DO M2 = -L2,L2
                  XL2 = L2
                  LM2 = L2*(L2+1) + M2 + 1
C
C---> STORE ONLY GAUNT COEFFIENTS FOR  LM2 <= LM1
C
                  IF ( LM2.LE.LM1 ) THEN
                     LM1LM2 = LM1LM2 + 1
                     NRGNT(LM1LM2) = 0
C
                     DO L3 = ABS(L1-L2),(L1+L2)
                        DO M3 = -L3,L3
                           XL3 = L3
C
                           RGAUNT = 0.0D0
                           DO I1 = 1,2
                              IF ( I1.EQ.1 ) THEN
                                 XM1 = -ABS(M1)
                                 IF ( M1.LT.0 ) CC1 = CI*WZ05
                                 IF ( M1.EQ.0 ) CC1 = 1.0D0
                                 IF ( M1.GT.0 ) CC1 = WZ05
                              ELSE
                                 XM1 = ABS(M1)
                                 IF ( M1.LT.0 ) CC1 = -CI*WZ05*(-1)
     &                                **ABS(M1)
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
                                    IF ( M2.LT.0 ) CC2 = -CI*WZ05*(-1)
     &                                 **ABS(M2)
                                    IF ( M2.EQ.0 ) CC2 = 0.0D0
                                    IF ( M2.GT.0 ) CC2 = WZ05*(-1)
     &                                 **ABS(M2)
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
                                       IF ( M3.LT.0 )
     &                                    CC3 = -CI*WZ05*(-1)**ABS(M3)
                                       IF ( M3.EQ.0 ) CC3 = 0.0D0
                                       IF ( M3.GT.0 ) CC3 = WZ05*(-1)
     &                                    **ABS(M3)
                                    END IF
C
                                    IF ( (NINT(XM1)+NINT(XM2)+NINT(XM3))
     &                                 .EQ.0 ) THEN
C
                                       CGAUNT = SQRT((2*XL1+1)*(2*XL2+1)
     &                                    /(4*PI*(2*XL3+1)))*(-1)
     &                                    **(-NINT(XM3))
     &                                    *CGC_RACAH(XL1,XL2,XL3,0.0D0,
     &                                    0.0D0,0.0D0)
     &                                    *CGC_RACAH(XL1,XL2,XL3,XM1,
     &                                    XM2,-XM3)
C
                                       RGAUNT = RGAUNT + 
     &                                    DREAL(CC1*CC2*CC3)*CGAUNT
                                    END IF
                                 END DO
                              END DO
                           END DO
C
                           IF ( ABS(RGAUNT).GT.1.0D-8 ) THEN
                              IG123 = IG123 + 1
                              IF ( IG123.GT.LRGNT123 ) THEN
                                 WRITE (6,99002) LRGNT123
                                 STOP 'in <STRGAUNT>'
                              END IF
                              NRGNT(LM1LM2) = NRGNT(LM1LM2) + 1
                              RGNT(IG123) = PRE*RGAUNT
                              IRGNT(IG123) = L3*(L3+1) + M3 + 1
C
                           END IF
C
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C=======================================================================
C
      WRITE (6,99001) LMAX,IG123
C
99001 FORMAT (/,10X,'for LMAX =',I2,I6,' non-zero GAUNTS tabulated',/)
99002 FORMAT (/,10X,'LRGNT123 =',I6,'  too small ',
     &        'increase array - size  ',/,10X,
     &        '    100  for  l_max = 2     ',/,10X,
     &        '    400  for  l_max = 3     ',/,10X,
     &        '   1200  for  l_max = 4     ')
      END
