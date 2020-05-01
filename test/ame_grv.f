C*==ame_grv.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AME_GRV
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed                   *
C   *                                                                  *
C   *   for the         grad V * A -  form                             *
C   *                                                                  *
C   *                       ^    ^                                     *
C   *          GV = < LAM|  r  . e(POL) | LAM'>     AV_GRV             *
C   *          GB = < LAM| S_3 r.e(POL) | LAM'>     AB_GRV             *
C   *          B1 = < LAM| SIG . e(POL) |-LAM'>     B1_GRV             *
C   *          B2 = <-LAM| SIG . e(POL) | LAM'>     B2_GRV             *
C   *          C1 = < LAM| SIG x e(POL) |-LAM'>/i   C1_GRV             *
C   *          C2 = <-LAM| SIG x e(POL) | LAM'>/i   C2_GRV             *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(+) = [e(x) + i e(y)] / sqrt(2)  left  circ. positive helicity *
C   *                                                                  *
C   *  ->      ->       ->                                             *
C   *  e(-) = [e(x) - i e(y)] / sqrt(2)  right circ. negative helicity *
C   *                                                                  *
C   * 28/10/94  HE                                                     *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:AV_GRV,AB_GRV,AB1_GRV,AB2_GRV,AC1_GRV,AC2_GRV,
     &    A1_ADA,A2_ADA,AB1_NAB,AB2_NAB,NKMMAX,NLMAX,CGC
      IMPLICIT NONE
C*--AME_GRV30
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 CGC1
      REAL*8 F,F1,F2,SGN(-1:+1),TOL
      INTEGER I,IDXPOL(-1:+1),IFLAG,IKM1,IKM2,J,J1P05,J2P05,K1,K2,KAP1,
     &        KAP2,L1,L2,MM,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      DATA IDXPOL/2,3,1/,SGN/ + 1.0D0, + 1.0D0, - 1.0D0/
C
      NK = 2*NLMAX - 1
C
      IF ( .NOT.ALLOCATED(AV_GRV) ) THEN
         ALLOCATE (AV_GRV(NKMMAX,NKMMAX,3),AB_GRV(NKMMAX,NKMMAX,3))
         ALLOCATE (AB1_GRV(NKMMAX,NKMMAX,3),AB2_GRV(NKMMAX,NKMMAX,3))
         ALLOCATE (AC1_GRV(NKMMAX,NKMMAX,3),AC2_GRV(NKMMAX,NKMMAX,3))
      END IF
C
      CALL AME_NAB
C
      AB1_GRV(1:NKMMAX,1:NKMMAX,1:3) = A1_ADA(1:NKMMAX,1:NKMMAX,1:3)
      AB2_GRV(1:NKMMAX,1:NKMMAX,1:3) = A2_ADA(1:NKMMAX,1:NKMMAX,1:3)
      AC1_GRV(1:NKMMAX,1:NKMMAX,1:3) = AB1_NAB(1:NKMMAX,1:NKMMAX,1:3)
      AC2_GRV(1:NKMMAX,1:NKMMAX,1:3) = AB2_NAB(1:NKMMAX,1:NKMMAX,1:3)
C
      AV_GRV(1:NKMMAX,1:NKMMAX,1:3) = 0D0
      AB_GRV(1:NKMMAX,1:NKMMAX,1:3) = 0D0
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
                  MM = MUE1M05 - MUE2M05
C-----------------------------------------------------------------------
                  IF ( ABS(L1-L2).EQ.1 ) THEN
                     IF ( ABS(MM).LE.1 ) THEN
                        F = SGN(MM)*SQRT(DBLE(2*L2+1)/DBLE(2*L1+1))
     &                      *CGC1(L2,L1,0,0)
                        F1 = CGC(IKM1,1)*CGC(IKM2,1)
     &                       *CGC1(L2,L1,(MUE2M05+1),MM)
                        F2 = CGC(IKM1,2)*CGC(IKM2,2)
     &                       *CGC1(L2,L1,(MUE2M05+0),MM)
                        AV_GRV(IKM1,IKM2,IDXPOL(MM)) = F*(F2+F1)
                        AB_GRV(IKM1,IKM2,IDXPOL(MM)) = F*(F2-F1)
                     END IF
                  END IF
C-----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
C=======================================================================
C                  check the symmetry of the matrices
C=======================================================================
      TOL = 1D-10
      IFLAG = 0
C
      DO I = 1,NKMMAX
         DO J = 1,NKMMAX
            IF ( ABS(AV_GRV(I,J,1)-AV_GRV(J,I,2)).GT.TOL ) THEN
               WRITE (6,99001) 'V',1,I,J
               IFLAG = 1
            END IF
            IF ( ABS(AV_GRV(I,J,3)-AV_GRV(J,I,3)).GT.TOL ) THEN
               WRITE (6,99001) 'V',3,I,J
               IFLAG = 1
            END IF
            IF ( ABS(AB_GRV(I,J,1)-AB_GRV(J,I,2)).GT.TOL ) THEN
               WRITE (6,99001) 'V',1,I,J
               IFLAG = 1
            END IF
            IF ( ABS(AB_GRV(I,J,3)-AB_GRV(J,I,3)).GT.TOL ) THEN
               WRITE (6,99001) 'V',3,I,J
               IFLAG = 1
            END IF
         END DO
      END DO
C
      IF ( IFLAG.NE.0 ) STOP 'in <AMEGRV>'
C
C=======================================================================
99001 FORMAT ('TROUBLE for TYPE=',A,' IPOL=',I3,' I=',I3,' J=',I3)
C
      END
