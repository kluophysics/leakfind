C*==meada_ameada_old_version.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MEADA_AMEADA_OLD_VERSION
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements                          *
C   *                                                                  *
C   *                               ->  ->                             *
C   *   A(LAM,LAM';POL) = <chi(LAM)|SIG*e(POL)|chi(LAM')>              *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(+) = [e(x) + i e(y)] / sqrt(2)  left  circ. positive helicity *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(-) = [e(x) - i e(y)] / sqrt(2)  right circ. negative helicity *
C   *                                                                  *
C   *   A1 = A( KAP,MUE;-KAP',MUE')   electric dipolar matrix elements *
C   *   A2 = A(-KAP,MUE; KAP',MUE')                                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:A1_ADA,A2_ADA,NKMMAX,AME_G,IMKM_IKM,NKM,ISMT
      IMPLICIT NONE
C*--MEADA_AMEADA_OLD_VERSION25
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IKM1,IKM2,IMKM1,IMKM2,IPOL
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.ALLOCATED(A1_ADA) )
     &     ALLOCATE (A1_ADA(NKMMAX,NKMMAX,3),A2_ADA(NKMMAX,NKMMAX,3))
C
      A1_ADA(1:NKMMAX,1:NKMMAX,1:3) = 0D0
      A2_ADA(1:NKMMAX,1:NKMMAX,1:3) = 0D0
C
      DO IKM1 = 1,NKM
         IMKM1 = IMKM_IKM(IKM1)
         DO IKM2 = 1,NKM
            IMKM2 = IMKM_IKM(IKM2)
            DO IPOL = 1,3
               IF ( IPOL.EQ.1 ) THEN
C--------------------------------------------------------------------(+)
                  A1_ADA(IKM1,IKM2,IPOL) = -AME_G(IKM1,IMKM2,3,ISMT)
                  A2_ADA(IKM1,IKM2,IPOL) = -AME_G(IMKM1,IKM2,3,ISMT)
               ELSE IF ( IPOL.EQ.2 ) THEN
C--------------------------------------------------------------------(-)
                  A1_ADA(IKM1,IKM2,IPOL) = +AME_G(IKM1,IMKM2,1,ISMT)
                  A2_ADA(IKM1,IKM2,IPOL) = +AME_G(IMKM1,IKM2,1,ISMT)
               ELSE
C--------------------------------------------------------------------(z)
                  A1_ADA(IKM1,IKM2,IPOL) = +AME_G(IKM1,IMKM2,2,ISMT)
                  A2_ADA(IKM1,IKM2,IPOL) = +AME_G(IMKM1,IKM2,2,ISMT)
               END IF
C-----------------------------------------------------------------------
            END DO
         END DO
      END DO
C
      END
C*==fpameada.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE FPAMEADA
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements                          *
C   *                                                                  *
C   *                                   ->  ->                         *
C   *   A(LAM,LAM';POL) = <chi(LAM)|Y_L SIG*e(POL)|chi(LAM')>          *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(+) = [e(x) + i e(y)] / sqrt(2)  left  circ. positive helicity *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(-) = [e(x) - i e(y)] / sqrt(2)  right circ. negative helicity *
C   *                                                                  *
C   *   A1 = A( KAP,MUE;-KAP',MUE')   electric dipolar matrix elements *
C   *   A2 = A(-KAP,MUE; KAP',MUE')                                    *
C   *   B1 = A( KAP,MUE;-KAP',MUE')   next higher terms in the         *
C   *   B2 = A(-KAP,MUE; KAP',MUE')   expansion of  exp(i*q*r)         *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   *  NOTE: the non-dipolar terms still have to be checked !!!!!!!!!  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:CGC,NLMAX,NKMMAX,NLAMEFPMAX,FPA1_ADA,FPA2_ADA,
     &    FPB1_ADA,FPB2_ADA
      USE MOD_CONSTANTS,ONLY:C0,PI
      IMPLICIT NONE
C*--FPAMEADA104
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 D,F,WZ2,WZ3,WZ4PI
      COMPLEX*16 FPA1_ADAC,FPA2_ADAC,FPB1_ADAC,FPB2_ADAC,YR(2)
      REAL*8 GAUNT_CYLM
      INTEGER IKM1,IKM2,ILA,IMKM1,IMKM2,IPOL(-1:1),IPOL1,J1P05,J2P05,K1,
     &        K2,KAP1,KAP2,L1,L2,LA,LAP,LB1,LB2,M,M1,M2,MA,MIPOL(3),
     &        MIPOL1,MSM05,MUE1M05,MUE2M05,N,NK
C
C*** End of declarations rewritten by SPAG
C
      DATA IPOL/2,3,1/
      DATA MIPOL/2,1,3/
C
      NK = 2*NLMAX - 1
      WZ2 = DSQRT(2.0D0)
      WZ4PI = DSQRT(4.0D0*PI)
      WZ3 = DSQRT(3.0D0)
C
      NLAMEFPMAX = 1 + 4*(NLMAX-1)
      IF ( .NOT.ALLOCATED(FPA1_ADA) ) THEN
         ALLOCATE (FPA1_ADA(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPA2_ADA(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPB1_ADA(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPB2_ADA(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
      END IF
C
      FPA1_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPA2_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPB1_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPB2_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
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
C
                  DO LA = 0,(2*NLMAX-1)
C
                     ILA = LA + 1
C
C------------- FPA1_ADA,FPA2_ADA:  for polarisations -,+  == IPOL1 = 2,1
C
                     DO M = -1, + 1,2
                        IPOL1 = IPOL(M)
C
                        MA = -M + MUE1M05 - MUE2M05
C
                        IF ( ABS(MA).LE.LA ) THEN
C
                           M1 = (2*MUE1M05+1-M)/2
                           M2 = (2*MUE2M05+1+M)/2
C
                           F = WZ2
                           FPA1_ADAC = F*CGC(IKM1,((3+M)/2))
     &                                 *CGC(IMKM2,((3-M)/2))
     &                                 *GAUNT_CYLM(L1,M1,LA,MA,LB2,M2)
C
                           FPA2_ADAC = F*CGC(IMKM1,((3+M)/2))
     &                                 *CGC(IKM2,((3-M)/2))
     &                                 *GAUNT_CYLM(LB1,M1,LA,MA,L2,M2)
C
                           CALL FPAMEAUX(FPA1_ADAC,YR,MA)
C
                           FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1) = YR(1)
                           FPA1_ADA(IKM1,IKM2,ILA,IPOL1,2) = YR(2)
C
                           CALL FPAMEAUX(FPA2_ADAC,YR,MA)
C
                           FPA2_ADA(IKM1,IKM2,ILA,IPOL1,1) = YR(1)
                           FPA2_ADA(IKM1,IKM2,ILA,IPOL1,2) = YR(2)
                        END IF
                     END DO
C
C----------------- FPA1_ADA,FPA2_ADA:  for polarisations z  == IPOL1 = 3
C
                     MA = MUE1M05 - MUE2M05
                     IPOL1 = IPOL(0)
C
                     IF ( ABS(MA).LE.ABS(LA) ) THEN
                        DO MSM05 = -1,0
                           M1 = MUE1M05 - MSM05
                           M2 = MUE2M05 - MSM05
                           F = 2.0D0*(MSM05+0.5D0)
                           FPA1_ADAC = F*CGC(IKM1,MSM05+2)
     &                                 *CGC(IMKM2,MSM05+2)
     &                                 *GAUNT_CYLM(L1,M1,LA,MA,LB2,M2)
C
                           FPA2_ADAC = F*CGC(IMKM1,MSM05+2)
     &                                 *CGC(IKM2,MSM05+2)
     &                                 *GAUNT_CYLM(LB1,M1,LA,MA,L2,M2)
C
                           CALL FPAMEAUX(FPA1_ADAC,YR,MA)
C
                           FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1)
     &                        = FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1) + YR(1)
                           FPA1_ADA(IKM1,IKM2,ILA,IPOL1,2)
     &                        = FPA1_ADA(IKM1,IKM2,ILA,IPOL1,2) + YR(2)
C
                           CALL FPAMEAUX(FPA2_ADAC,YR,MA)
C
                           FPA2_ADA(IKM1,IKM2,ILA,IPOL1,1)
     &                        = FPA2_ADA(IKM1,IKM2,ILA,IPOL1,1) + YR(1)
                           FPA2_ADA(IKM1,IKM2,ILA,IPOL1,2)
     &                        = FPA2_ADA(IKM1,IKM2,ILA,IPOL1,2) + YR(2)
C
                        END DO
                     END IF
C
C------------- FPB1_ADA,FPB2_ADA:  for polarisations -,+  == IPOL1 = 2,1
C
                     F = WZ2*WZ4PI/WZ3
C
                     DO M = -1, + 1,2
C
                        IPOL1 = IPOL(M)
                        MA = -M + MUE1M05 - MUE2M05
C
                        IF ( ABS(MA).LE.ABS(LA) ) THEN
C
                           FPB1_ADAC = C0
                           FPB2_ADAC = C0
C
                           DO LAP = ABS(LA-1),LA + 1
C
                              M1 = (2*MUE1M05+1-M)/2
                              M2 = (2*MUE2M05+1+M)/2
C
                              FPB1_ADAC = FPB1_ADAC + 
     &                           F*CGC(IKM1,((3+M)/2))
     &                           *CGC(IMKM2,((3-M)/2))
     &                           *GAUNT_CYLM(L1,M1,LAP,MA,LB2,M2)
     &                           *GAUNT_CYLM(LAP,MA,LA,MA,1,0)
C
                              FPB2_ADAC = FPB2_ADAC + 
     &                           F*CGC(IMKM1,((3+M)/2))
     &                           *CGC(IKM2,((3-M)/2))
     &                           *GAUNT_CYLM(LB1,M1,LAP,MA,L2,M2)
     &                           *GAUNT_CYLM(LAP,MA,LA,MA,1,0)
C
                           END DO
C
                           CALL FPAMEAUX(FPB1_ADAC,YR,MA)
C
                           FPB1_ADA(IKM1,IKM2,ILA,IPOL1,1) = YR(1)
                           FPB1_ADA(IKM1,IKM2,ILA,IPOL1,2) = YR(2)
C
                           CALL FPAMEAUX(FPB2_ADAC,YR,MA)
C
                           FPB2_ADA(IKM1,IKM2,ILA,IPOL1,1) = YR(1)
                           FPB2_ADA(IKM1,IKM2,ILA,IPOL1,2) = YR(2)
C
                        END IF
                     END DO
C
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
C symmetry check for FPA1_ADA/2 from <SCFBIAME> --- only for MUE1 = MUE2
C
      N = 0
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
C
                  IF ( (MUE1M05-MUE2M05).EQ.0 ) THEN
                     DO ILA = 1,2*NLMAX
                        DO IPOL1 = 1,3
                           MIPOL1 = MIPOL(IPOL1)
                           IF ( ABS(FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1))
     &                          .GT.1.0D-10 ) N = N + 1
                           IF ( IPOL1.NE.3 ) THEN
                              F = +1.0D0
                           ELSE
                              F = -1.0D0
                           END IF
                           D = ABS(FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1)
     &                         +F*FPA2_ADA(IKM2,IKM1,ILA,MIPOL1,1))
                           IF ( D.GT.1.D-10 ) WRITE (6,'(A,4I3,3E15.5)')
     &                           ' trouble in <FPAMEADA>',IKM1,IKM2,ILA,
     &                          IPOL1,FPA1_ADA(IKM1,IKM2,ILA,IPOL1,1),D
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
C
      END
C*==fpameaux.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE FPAMEAUX(YC,YR,M)
C   ********************************************************************
C   *                                                                  *
C   *    perform transformation from complex spherical harmonic        *
C   *    representation to real spherical harmonic representation      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,C1,C0
      IMPLICIT NONE
C*--FPAMEAUX376
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M
      COMPLEX*16 YC
      COMPLEX*16 YR(2)
C
C Local variables
C
      COMPLEX*16 FAC1,FACI,SQRT2
      INTEGER MABS
C
C*** End of declarations rewritten by SPAG
C
      SQRT2 = DCMPLX(DSQRT(2.0D0))
      FACI = CI/SQRT2
      FAC1 = C1/SQRT2
C
      MABS = ABS(M)
      IF ( MABS.EQ.0 ) THEN
         YR(1) = YC
         YR(2) = C0
      ELSE IF ( MABS.EQ.M ) THEN
         YR(1) = FACI*DCMPLX((-1)**(MABS+1))*YC
         YR(2) = FAC1*DCMPLX((-1)**MABS)*YC
      ELSE IF ( MABS.EQ.-M ) THEN
         YR(1) = FACI*YC
         YR(2) = FAC1*YC
      END IF
C
      END
C*==fpamenab.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE FPAMENAB
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed                   *
C   *                                                                  *
C   *   for the          NABLA * A -  form                             *
C   *                                                                  *
C   *                       ->   ->                                    *
C   *          A1 = < LAM| NAB . e(POL) | LAM'>     l = l'+1           *
C   *          A2 = < LAM| NAB . e(POL) | LAM'>     l = l'-1           *
C   *          B1 = < LAM| SIG . e(POL) |-LAM'>     <AMEADA>           *
C   *          B2 = <-LAM| SIG . e(POL) | LAM'>     <AMEADA>           *
C   *          C1 = < LAM| SIG x e(POL) |-LAM'>/i                      *
C   *          C2 = <-LAM| SIG x e(POL) | LAM'>/i                      *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(+) = [e(x) + i e(y)] / sqrt(2)  left  circ. positive helicity *
C   *   ps                                                             *
C   *  ->      ->       ->                                             *
C   *  e(-) = [e(x) - i e(y)] / sqrt(2)  right circ. negative helicity *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   *     (LA,MA)  quantum numbers for real sperical harmonics         *
C   *              in shape functions                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:CGC,NLMAX,NKMMAX,NLAMEFPMAX,FPA1_ADA,FPA2_ADA,
     &    FPAN1_NAB,FPAN2_NAB,FPAV1_NAB,FPAV2_NAB,FPAB1_NAB,FPAB2_NAB
      IMPLICIT NONE
C*--FPAMENAB455
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPAMENAB')
C
C Local variables
C
      REAL*8 CGC1,GAUNT_CYLM
      COMPLEX*16 FPAN1_NABC,FPAN2_NABC,YR(2)
      INTEGER IDXPOL(-1:+1),IKM1,IKM2,ILA,IPOL,J1P05,J2P05,K1,K2,KAP1,
     &        KAP2,L1,L2,LA,M,M1,M2,MA,MSM05,MUE1M05,MUE2M05,NK
      REAL*8 SGN(-1:+1)
C
C*** End of declarations rewritten by SPAG
C
      DATA IDXPOL/2,3,1/
      DATA SGN/ + 1.0D0, + 1.0D0, - 1.0D0/
C
      CALL TRACK_INFO(ROUTINE)
C
      NK = 2*NLMAX - 1
C
      CALL FPAMEADA
C
      NLAMEFPMAX = 1 + 4*(NLMAX-1)
      IF ( .NOT.ALLOCATED(FPAN1_NAB) ) THEN
         ALLOCATE (FPAN1_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPAN2_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPAV1_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPAV2_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPAB1_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
         ALLOCATE (FPAB2_NAB(NKMMAX,NKMMAX,NLAMEFPMAX,3,2))
      END IF
C
      FPAN1_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPAN2_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPAV1_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2)
     &   = FPA1_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2)
      FPAV2_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2)
     &   = FPA2_ADA(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2)
      FPAB1_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
      FPAB2_NAB(1:NKMMAX,1:NKMMAX,1:NLAMEFPMAX,1:3,1:2) = C0
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
C
C ----------------------------------------------------------------------
C
                  DO LA = 0,(2*NLMAX-1)
C
                     ILA = LA + 1
                     DO M = -1, + 1
                        IPOL = IDXPOL(M)
C
                        MA = -M + MUE1M05 - MUE2M05
C
                        IF ( ABS(MA).LE.LA ) THEN
C-----------------------------------------------------------------------
                           DO MSM05 = -1,0
                              M1 = MUE1M05 - MSM05
                              M2 = MUE2M05 - MSM05 + M
C
                              FPAN1_NABC = SGN(M)
     &                           *SQRT(DBLE(L2+1)/DBLE(2*L2+3))
     &                           *CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                           *CGC1(L2,L2+1,(MUE2M05-MSM05),M)
     &                           *GAUNT_CYLM(L1,M1,LA,MA,L2+1,M2)
C
                              FPAN2_NABC = SGN(M)
     &                           *SQRT(DBLE(L2)/DBLE(2*L2-1))
     &                           *CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                           *CGC1(L2,L2-1,(MUE2M05-MSM05),M)
     &                           *GAUNT_CYLM(L1,M1,LA,MA,L2-1,M2)
C
                              CALL FPAMEAUX(FPAN1_NABC,YR,MA)
C
                              FPAN1_NAB(IKM1,IKM2,ILA,IPOL,1)
     &                           = FPAN1_NAB(IKM1,IKM2,ILA,IPOL,1)
     &                           + YR(1)
                              FPAN1_NAB(IKM1,IKM2,ILA,IPOL,2)
     &                           = FPAN1_NAB(IKM1,IKM2,ILA,IPOL,2)
     &                           + YR(2)
C
                              CALL FPAMEAUX(FPAN2_NABC,YR,MA)
C
                              FPAN2_NAB(IKM1,IKM2,ILA,IPOL,1)
     &                           = FPAN2_NAB(IKM1,IKM2,ILA,IPOL,1)
     &                           + YR(1)
                              FPAN2_NAB(IKM1,IKM2,ILA,IPOL,2)
     &                           = FPAN2_NAB(IKM1,IKM2,ILA,IPOL,2)
     &                           + YR(2)
C
                           END DO
C-----------------------------------------------------------------------
C   AMEAC    ............. to be done
C-----------------------------------------------------------------------
                        END IF
                     END DO
                  END DO
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
      END
