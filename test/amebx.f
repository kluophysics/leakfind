C*==amebx.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMEBX(AMEBXG,AMEBXF,NLMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements for the                  *
C   *                                                                  *
C   *            B_x  +  COULOMB INTERACTION                           *
C   *                                                                  *
C   *  AMEBXG(IKM1,IKM2,ILR) =                                         *
C   *                  < KAP1,MUE1 | SIGMA_Z (Y[LR,MR])* |  KAP2,MUE2> *
C   *                                                                  *
C   *  AMEBXF(IKM1,IKM2,ILR) =                                         *
C   *                  <-KAP1,MUE1 | SIGMA_Z (Y[LR,MR])* | -KAP2,MUE2> *
C   *                                                                  *
C   *  22/01/98  HE                                                    *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:CGC
      IMPLICIT NONE
C*--AMEBX19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX,NLMAX
      REAL*8 AMEBXF(NKMMAX,NKMMAX,2*NLMAX),AMEBXG(NKMMAX,NKMMAX,2*NLMAX)
C
C Local variables
C
      REAL*8 F
      REAL*8 GAUNT_CYLM
      INTEGER IKM1,IKM2,ILR,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,
     &        L2,LB1,LB2,LR,M1,M2,MR,MSM05,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      NK = 2*NLMAX - 1
C
      CALL RINIT(NKMMAX**2*2*NLMAX,AMEBXG)
      CALL RINIT(NKMMAX**2*2*NLMAX,AMEBXF)
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
                        MR = M2 - M1
                        F = (-1)**MR*(2*MSM05+1)
C
C                G = INT dr^  Y[l1,m1]* Y[l2,-m2] Y[l3,m3]
C
                        AMEBXG(IKM1,IKM2,ILR) = AMEBXG(IKM1,IKM2,ILR)
     &                     + F*CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                     *GAUNT_CYLM(L1,M1,LR,-MR,L2,M2)
C
                        AMEBXF(IKM1,IKM2,ILR) = AMEBXF(IKM1,IKM2,ILR)
     &                     + F*CGC(IMKM1,MSM05+2)*CGC(IMKM2,MSM05+2)
     &                     *GAUNT_CYLM(LB1,M1,LR,-MR,LB2,M2)
C
                     END DO
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
      END
