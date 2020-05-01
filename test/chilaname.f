C*==chilaname.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHILANAME(AMEA1,AMEA2,AMEB1,AMEB2,AMEB1C,AMEB2C)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed                   *
C   *                                                                  *
C   *   for the          NABLA * A -  form                             *
C   *                                                                  *
C   *                         ->   ->                                  *
C   *          A1 = < LAM|   NAB . e(POL) | LAM'>     l = l'+1         *
C   *          A2 = < LAM|   NAB . e(POL) | LAM'>     l = l'-1         *
C   *          B1 = < LAM| y NAB . e(POL) | LAM'>     l = l'+1         *
C   *          B2 = < LAM| y NAB . e(POL) | LAM'>     l = l'-1         *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(+) = [e(x) + i e(y)] / sqrt(2)  left  circ. positive helicity *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(-) = [e(x) - i e(y)] / sqrt(2)  right circ. negative helicity *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:CGC,NKMMAX,NLMAX
      IMPLICIT NONE
C*--CHILANAME26
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 AMEA1(NKMMAX,NKMMAX,3),AMEA2(NKMMAX,NKMMAX,3),
     &       AMEB1(NKMMAX,NKMMAX,3,3),AMEB1C(NKMMAX,NKMMAX,3,3),
     &       AMEB2(NKMMAX,NKMMAX,3,3),AMEB2C(NKMMAX,NKMMAX,3,3)
C
C Local variables
C
      REAL*8 CCDN,CCUP,SGN(-1:+1),SRL1,SRL2
      REAL*8 CGC1,GAUNT_CYLM
      INTEGER IC,IDXPOL(-1:+1),IKM1,IKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,
     &        L2,MC,MM,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      DATA IDXPOL/2,3,1/,SGN/ + 1.0D0, + 1.0D0, - 1.0D0/
C
      NK = 2*NLMAX - 1
C
      CALL RINIT(NKMMAX**2*3,AMEA1)
      CALL RINIT(NKMMAX**2*3,AMEA2)
      CALL RINIT(NKMMAX**2*3*3,AMEB1)
      CALL RINIT(NKMMAX**2*3*3,AMEB2)
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
                  CCDN = CGC(IKM1,1)*CGC(IKM2,1)
                  CCUP = CGC(IKM1,2)*CGC(IKM2,2)
C
C-----------------------------------------------------------------------
                  MM = MUE1M05 - MUE2M05
                  IF ( ABS(MM).LE.1 ) THEN
C
                     SRL1 = SGN(MM)*SQRT(DBLE(L2+1)/DBLE(2*L2+3))
                     SRL2 = SGN(MM)*SQRT(DBLE(L2)/DBLE(2*L2-1))
C
                     IF ( L1.EQ.L2+1 ) AMEA1(IKM1,IKM2,IDXPOL(MM))
     &                    = SRL1*(CCDN*CGC1(L2,L2+1,(MUE2M05+1),MM)
     &                    +CCUP*CGC1(L2,L2+1,(MUE2M05+0),MM))
C
                     IF ( L1.EQ.L2-1 ) AMEA2(IKM1,IKM2,IDXPOL(MM))
     &                    = SRL2*(CCDN*CGC1(L2,L2-1,(MUE2M05+1),MM)
     &                    +CCUP*CGC1(L2,L2-1,(MUE2M05+0),MM))
                  END IF
C
C-----------------------------------------------------------------------
C
                  DO MM = -1, + 1
                     SRL1 = SGN(MM)*SQRT(DBLE(L2+1)/DBLE(2*L2+3))
                     SRL2 = SGN(MM)*SQRT(DBLE(L2)/DBLE(2*L2-1))
C
                     DO MC = -1, + 1
                        IC = MC + 2
C
                        AMEB1(IKM1,IKM2,IC,IDXPOL(MM))
     &                     = SRL1*(CCDN*CGC1(L2,L2+1,(MUE2M05+1),MM)
     &                     *GAUNT_CYLM(L1,(MUE1M05+1),1,MC,L2+1,
     &                     (MUE2M05+1+MM))
     &                     +CCUP*CGC1(L2,L2+1,(MUE2M05+0),MM)
     &                     *GAUNT_CYLM(L1,(MUE1M05+0),1,MC,L2+1,
     &                     (MUE2M05+0+MM)))
C
                        AMEB2(IKM1,IKM2,IC,IDXPOL(MM))
     &                     = SRL2*(CCDN*CGC1(L2,L2-1,(MUE2M05+1),MM)
     &                     *GAUNT_CYLM(L1,(MUE1M05+1),1,MC,L2-1,
     &                     (MUE2M05+1+MM))
     &                     +CCUP*CGC1(L2,L2-1,(MUE2M05+0),MM)
     &                     *GAUNT_CYLM(L1,(MUE1M05+0),1,MC,L2-1,
     &                     (MUE2M05+0+MM)))
C
                        AMEB1C(IKM1,IKM2,IC,IDXPOL(MM))
     &                     = SRL1*(CCDN*CGC1(L1,L1+1,(MUE1M05+1),-MM)
     &                     *GAUNT_CYLM(L1+1,(MUE1M05+1-MM),1,MC,L2,
     &                     (MUE2M05+1))
     &                     +CCUP*CGC1(L1,L1+1,(MUE1M05+0),-MM)
     &                     *GAUNT_CYLM(L1+1,(MUE1M05+0-MM),1,MC,L2,
     &                     (MUE2M05+0)))
C
                        AMEB2C(IKM1,IKM2,IC,IDXPOL(MM))
     &                     = SRL2*(CCDN*CGC1(L1,L1-1,(MUE1M05+1),-MM)
     &                     *GAUNT_CYLM(L1-1,(MUE1M05+1+MM),1,MC,L2,
     &                     (MUE2M05+1))
     &                     +CCUP*CGC1(L1,L1-1,(MUE1M05+0),-MM)
     &                     *GAUNT_CYLM(L1-1,(MUE1M05+0+MM),1,MC,L2,
     &                     (MUE2M05+0)))
C
                     END DO
                  END DO
C
C-----------------------------------------------------------------------
C
               END DO
            END DO
         END DO
      END DO
C
      END
