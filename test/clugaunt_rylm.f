C*==clugaunt_rylm.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUGAUNT_RYLM(LMAX,RGNT,NRGNT,IL3RGNT,LM3RGNT,NRGNT12,
     &                         NRGNT123,NRGNT12MAX,NRGNT123MAX)
C   ********************************************************************
C   *                                                                  *
C   *    calculation of modified   Gaunt coefficients                  *
C   *    for cluster calculations and REAL spherical harmonics         *
C   *                                                                  *
C   *    G(l1,m1;l2,m2;l3,m3) ==  - 4 PI * i^(l1-l2-l3)                *
C   *                            * INT Y(l1,m1)  Y(l2,m2)  Y(l3,m3)    *
C   *                                                                  *
C   *    the Gaunts are calculated only for  for LM2 <= LM1            *
C   *                                                                  *
C   *    for the loop LM1=1,LM1MAX / LM2=1,LM1 / LM3=1,LM3MAX          *
C   *    NRGNT       specifies the number of non-zero RGNT's           *
C   *                for the current linear index I == (LM1,LM2)       *
C   *    LM3RGNT(I)  gives the  LM3  for the I-th non-zero RGNT        *
C   *                the arraylength  NRGNT123MAX  has to be chosen    *
C   *                accordingly   --  see below                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--CLUGAUNT_RYLM25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LMAX,NRGNT12,NRGNT123,NRGNT123MAX,NRGNT12MAX
      INTEGER IL3RGNT(NRGNT123MAX),LM3RGNT(NRGNT123MAX),
     &        NRGNT(NRGNT12MAX)
      REAL*8 RGNT(NRGNT123MAX)
C
C Local variables
C
      REAL*8 G,PRE
      REAL*8 GAUNT_RYLM
      INTEGER IG123,L1,L2,L3,LEXP,LM1,LM1LM2,LM2,M1,M2,M3
C
C*** End of declarations rewritten by SPAG
C
      PRE = -4*PI
      LM1LM2 = 0
      IG123 = 0
C
C=======================================================================
C
      DO L1 = 0,LMAX
         DO M1 = -L1,L1
C
            LM1 = L1*(L1+1) + M1 + 1
C
C=======================================================================
            DO L2 = 0,LMAX
               DO M2 = -L2,L2
C
                  LM2 = L2*(L2+1) + M2 + 1
C
C***********************************************************************
                  IF ( LM2.LE.LM1 ) THEN
C
                     LM1LM2 = LM1LM2 + 1
                     NRGNT(LM1LM2) = 0
C
C=======================================================================
                     DO L3 = ABS(L1-L2),(L1+L2),2
                        IF ( MOD((L1-L2-L3),2).NE.0 ) STOP 
     &                       'in <CLUGAUNT_RYLM>:   L1-L2-L3 <> n*2'
C
                        LEXP = (L1-L2-L3)/2
C
                        DO M3 = -L3,L3
C
                           G = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C
                           IF ( ABS(G).GT.1.0D-8 ) THEN
                              IF ( LM1LM2.LE.NRGNT12MAX ) NRGNT(LM1LM2)
     &                             = NRGNT(LM1LM2) + 1
                              IG123 = IG123 + 1
C
                              IF ( IG123.LE.NRGNT123MAX ) THEN
                                 RGNT(IG123) = PRE*(-1)**LEXP*G
                                 LM3RGNT(IG123) = L3*(L3+1) + M3 + 1
                                 IL3RGNT(IG123) = L3 + 1
                              END IF
C
                           END IF
C
                        END DO
                     END DO
C=======================================================================
C
                  END IF
C***********************************************************************
C
               END DO
            END DO
C=======================================================================
C
         END DO
      END DO
C=======================================================================
C
      WRITE (6,99001) LMAX,IG123
C
      NRGNT12 = LM1LM2
      NRGNT123 = IG123
C
      IF ( NRGNT12.GT.NRGNT12MAX ) THEN
         WRITE (6,*) '   NRGNT12MAX =',NRGNT12MAX,'  too small'
         STOP 'in <CLUGAUNT_RYLM> '
      END IF
      IF ( NRGNT123.GT.NRGNT123MAX ) THEN
         WRITE (6,*) '   NRGNT123MAX =',NRGNT123MAX,'  too small'
         STOP 'in <CLUGAUNT_RYLM> '
      END IF
C
99001 FORMAT (/,10X,'for lmax =',I2,I6,' non-zero GAUNTS tabulated ',/)
      END
