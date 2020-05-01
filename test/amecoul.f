C*==amecoul.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMECOUL(AMECIG,AMECIF,IPRINT,NLMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements for the                  *
C   *                                                                  *
C   *                    COULOMB INTERACTION                           *
C   *                                                                  *
C   *  AMECIG(IKM1,IKM2,ILR) = < KAP1,MUE1 | Y[LR,MR] |  KAP2,MUE2>    *
C   *  AMECIF(IKM1,IKM2,ILR) = <-KAP1,MUE1 | Y[LR,MR] | -KAP2,MUE2>    *
C   *                                                                  *
C   *  22/01/98  HE                                                    *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_ANGMOM,ONLY:FACT,CGC
      IMPLICIT NONE
C*--AMECOUL17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NKMMAX,NLMAX
      REAL*8 AMECIF(NKMMAX,NKMMAX,2*NLMAX),AMECIG(NKMMAX,NKMMAX,2*NLMAX)
C
C Local variables
C
      REAL*8 A1,A2,A3,AME,B,BME,C,CJLJ,D1,D2,D3,DJ1,DJ2,DLR,DM1,DM2,F,H,
     &       W3J,X4PI
      REAL*8 CGC_RACAH,GAUNT_CYLM,WIG_3J
      INTEGER IKM1,IKM2,ILR,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,
     &        L2,LB1,LB2,LR,M1,M2,MR,MSM05,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      DATA DJ1/0D0/,DM1/0D0/,DJ2/0D0/,DM2/0D0/
C
      NK = 2*NLMAX - 1
C
      CALL RINIT(NKMMAX**2*2*NLMAX,AMECIG)
      CALL RINIT(NKMMAX**2*2*NLMAX,AMECIF)
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
            LB1 = L1 - 1
         ELSE
            KAP1 = -L1 - 1
            LB1 = L1 + 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
            IMKM1 = LB1*2*J1P05 + J1P05 + MUE1M05 + 1
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
                  LB2 = L2 - 1
               ELSE
                  KAP2 = -L2 - 1
                  LB2 = L2 + 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
                  IMKM2 = LB2*2*J2P05 + J2P05 + MUE2M05 + 1
C ----------------------------------------------------------------------
                  ILR = 0
                  DO LR = ABS(L1-L2),(L1+L2),2
                     ILR = LR + 1
                     DO MSM05 = -1,0
                        M1 = MUE1M05 - MSM05
                        M2 = MUE2M05 - MSM05
                        MR = M1 - M2
C
C                G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]
C
                        AMECIG(IKM1,IKM2,ILR) = AMECIG(IKM1,IKM2,ILR)
     &                     + CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                     *GAUNT_CYLM(L1,M1,LR,MR,L2,M2)
C
                        AMECIF(IKM1,IKM2,ILR) = AMECIF(IKM1,IKM2,ILR)
     &                     + CGC(IMKM1,MSM05+2)*CGC(IMKM2,MSM05+2)
     &                     *GAUNT_CYLM(LB1,M1,LR,MR,LB2,M2)
C
                     END DO
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
C=======================================================================
C               compare with expression of Shadwick et al.
C=======================================================================
C
      IF ( IPRINT.LE.2 ) RETURN
C
      WRITE (6,'("amcoul::")')
      X4PI = 1D0/SQRT(4*PI)
      H = 0.5D0
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
         ELSE
            KAP1 = -L1 - 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
               ELSE
                  KAP2 = -L2 - 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                  ILR = 0
                  DO LR = ABS(L1-L2),(L1+L2),2
                     ILR = LR + 1
C
                     AME = AMECIG(IKM1,IKM2,ILR)
C
C
                     MR = MUE1M05 - MUE2M05
C
                     F = 2*SQRT(2D0*LR+1D0)*(-1)**(J2P05+J1P05-1+LR)
C
C                     write(*,*) '###',J2P05,J1P05,LR,J2P05+J1P05-1-LR
                     IF ( ABS(J2P05-J1P05).LE.LR .AND. 
     &                    LR.LE.(J2P05+J1P05-1) ) THEN
C
                        A1 = FACT(J2P05+J1P05-1-LR)
                        A2 = FACT(J2P05-J1P05+LR)
                        A3 = FACT(J1P05-J2P05+LR)
                        B = FACT(J2P05+J1P05+LR)
C
                        C = FACT(INT((J2P05+J1P05+LR)/2D0))
                        D1 = FACT(INT((J2P05+J1P05-1-LR)/2D0))
                        D2 = FACT(INT((J1P05-J2P05+LR)/2D0))
                        D3 = FACT(INT((J2P05-J1P05+LR)/2D0))
C
                        CJLJ = F*SQRT(((A1*A2*A3)/B))*(C/(D1*D2*D3))
C
                        W3J = WIG_3J((J1P05-H),DBLE(LR),(J2P05-H),
     &                        -(MUE1M05+H),DBLE(MR),(MUE2M05+H))
C
                        BME = X4PI*CJLJ*(-1)**(MUE1M05+1)*W3J
C
                     ELSE
                        BME = 0D0
                     END IF
C
                     IF ( ABS(AME-BME).GT.1D-10 ) WRITE (6,99001) IKM1,
     &                    IKM2,LR,AME,BME,AME/BME
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C=======================================================================
C               compare with expression of EE
C=======================================================================
C
      WRITE (6,'("amcoul:: EE")')
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
         ELSE
            KAP1 = -L1 - 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
               ELSE
                  KAP2 = -L2 - 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                  ILR = 0
                  DO LR = ABS(L1-L2),(L1+L2),2
                     ILR = LR + 1
C
                     AME = AMECIG(IKM1,IKM2,ILR)
C
                     IF ( ABS(J2P05-J1P05).LE.LR .AND. 
     &                    LR.LE.(J2P05+J1P05-1) ) THEN
C
                        DJ1 = DBLE(J1P05) - 0.5D0
                        DJ2 = DBLE(J2P05) - 0.5D0
                        DLR = DBLE(LR)
                        DM1 = DBLE(MUE1M05) + 0.5D0
                        DM2 = DBLE(MUE2M05) + 0.5D0
C
                        A1 = (-1D0)**(MUE2M05+1)
                        A2 = (1D0+(-1)**(L1+L2+LR))/2D0
                        A3 = SQRT((2D0*DJ1+1D0)*(2D0*DJ2+1D0)/(4D0*PI)
     &                       /(2*LR+1D0))
C
                        D1 = CGC_RACAH(DJ1,DJ2,DLR,0.5D0,-0.5D0,0D0)
                        D2 = CGC_RACAH(DJ1,DJ2,DLR,-DM1,DM2,DM2-DM1)
C
                        BME = A1*A2*A3*D1*D2
C
                     ELSE
                        BME = 0D0
                     END IF
C
                     IF ( ABS(AME).GT.1D-10 .OR. ABS(BME).GT.1D-10 )
     &                    WRITE (6,99002) L1,DJ1,DM1,L2,DJ2,DM2,
     &                           DM2 - DM1,LR,AME,BME,AME/BME
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
C      STOP 'amecoul::'
C
99001 FORMAT (/,3I4,f20.12,/,12X,2F20.12)
99002 FORMAT (/,2(I4,2F5.1),f5.2,i4,3F12.8)
C
      END
