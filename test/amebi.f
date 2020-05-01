C*==amebi.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMEBI
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements occuring with the        *
C   *                      BREIT INTERACTION                           *
C   *                                                                  *
C   *  A1 = < KAP1,MUE1|Y[LA,-M]*SIG[M]|-KAP2,MUE2>                    *
C   *  A2 = <-KAP1,MUE1|Y[LA,-M]*SIG[M]| KAP2,MUE2>                    *
C   *                                                                  *
C   *  G1 = G[1,-m;L1,M1;L1-1,M1M] * G[L2,M2;1,-m';L1-1,M1M]           *
C   *  G2 = G[1,-m;L2,M2;L2+1,M2P] * G[L1,M1;1,-m';L2+1,M2P]           *
C   *                                                                  *
C   *  26/01/95  HE                                                    *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:AMEBI1,AMEBI2,GBIL,GBIG,CGC,NLMAX,NLABIMAX,
     &    NKMMAX
      IMPLICIT NONE
C*--AMEBI19
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 D,F,WZ2
      REAL*8 GAUNT_CYLM
      INTEGER IA_ERR,IKM1,IKM2,ILA,ILY1,IMKM1,IMKM2,IOFF,J1P05,J2P05,K1,
     &        K2,KAP1,KAP2,L1,L2,LA,LB1,LB2,LY1,LY1M,LY2,LY2P,M,M1,M2,
     &        MA,MP,MSM05,MUE1M05,MUE2M05,N,NK
C
C*** End of declarations rewritten by SPAG
C
      IF ( ALLOCATED(GBIG) ) DEALLOCATE (GBIL,GBIG)
      ALLOCATE (GBIG(2,NLABIMAX,0:1,-1:+1,-1:+1))
      ALLOCATE (GBIL(2,NLABIMAX,0:1,-1:+1,-1:+1))
C
      IF ( ALLOCATED(AMEBI1) ) DEALLOCATE (AMEBI1,AMEBI2)
      ALLOCATE (AMEBI1(NKMMAX,NKMMAX,NLABIMAX,-1:+1))
      ALLOCATE (AMEBI2(NKMMAX,NKMMAX,NLABIMAX,-1:+1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: AMEBI -> AMEBI1'
C
      CALL RINIT(NKMMAX**2*NLABIMAX*3,AMEBI1)
      CALL RINIT(NKMMAX**2*NLABIMAX*3,AMEBI2)
C
      NK = 2*NLMAX - 1
      WZ2 = DSQRT(2.0D0)
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
         J1P05 = ABS(KAP1)
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
               J2P05 = ABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
                  IMKM2 = LB2*2*J2P05 + J2P05 + MUE2M05 + 1
C ----------------------------------------------------------------------
                  IF ( (MUE1M05-MUE2M05).EQ.0 .AND. L1.EQ.L2 ) THEN
                     DO LA = 1,(2*NLMAX-1),2
                        ILA = (LA+1)/2
C
                        DO M = -1, + 1,2
                           M1 = (2*MUE1M05+1-M)/2
                           M2 = (2*MUE2M05+1+M)/2
                           MA = -M
                           F = (-1)**INT((1+M)/2)*WZ2
                           AMEBI1(IKM1,IKM2,ILA,M)
     &                        = F*CGC(IKM1,((3+M)/2))
     &                        *CGC(IMKM2,((3-M)/2))
     &                        *GAUNT_CYLM(L1,M1,LA,MA,LB2,M2)
C
                           AMEBI2(IKM1,IKM2,ILA,M)
     &                        = F*CGC(IMKM1,((3+M)/2))
     &                        *CGC(IKM2,((3-M)/2))
     &                        *GAUNT_CYLM(LB1,M1,LA,MA,L2,M2)
                        END DO
C
                        DO MSM05 = -1,0
                           M1 = MUE1M05 + 1 - MSM05
                           M2 = MUE2M05 + 1 - MSM05
                           F = 2.0D0*(MSM05+0.5D0)
                           AMEBI1(IKM1,IKM2,ILA,0)
     &                        = AMEBI1(IKM1,IKM2,ILA,0)
     &                        + F*CGC(IKM1,MSM05+2)*CGC(IMKM2,MSM05+2)
     &                        *GAUNT_CYLM(L1,M1,LA,0,LB2,M2)
C
                           AMEBI2(IKM1,IKM2,ILA,0)
     &                        = AMEBI2(IKM1,IKM2,ILA,0)
     &                        + F*CGC(IMKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                        *GAUNT_CYLM(LB1,M1,LA,0,L2,M2)
                        END DO
                     END DO
                  END IF
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C ======================================================================
C           product of gaunt numbers for retardation term
C ======================================================================
      CALL RINIT(2*NLABIMAX*2*3*3,GBIL)
      CALL RINIT(2*NLABIMAX*2*3*3,GBIG)
C
      DO LY1 = 1,(2*NLMAX-1),2
         ILY1 = (LY1+1)/2
         LY1M = LY1 - 1
C
         DO IOFF = 0,1
            LY2 = LY1 - IOFF*2
            LY2P = LY2 + 1
C
            DO M = -1, + 1
C
               DO MP = -1, + 1
C
                  GBIL(1,ILY1,IOFF,M,MP)
     &               = GAUNT_CYLM(1,-M,LY1,-M,LY1M,0)
     &               *GAUNT_CYLM(LY2,-MP,1,-MP,LY1M,0)
C
                  GBIL(2,ILY1,IOFF,M,MP)
     &               = GAUNT_CYLM(1,-M,LY2,-MP,LY2P,(MP-M))
     &               *GAUNT_CYLM(LY1,-M,1,-MP,LY2P,(MP-M))
C
                  GBIG(1,ILY1,IOFF,M,MP)
     &               = GAUNT_CYLM(1,-M,LY1,-MP,LY1M,(MP-M))
     &               *GAUNT_CYLM(LY2,-M,1,-MP,LY1M,(MP-M))
C
                  GBIG(2,ILY1,IOFF,M,MP)
     &               = GAUNT_CYLM(1,-M,LY2,-M,LY2P,0)
     &               *GAUNT_CYLM(LY1,-MP,1,-MP,LY2P,0)
C
               END DO
            END DO
         END DO
      END DO
C
      WRITE (6,*) ' <SCFBIAME> checking symmetry of AMEBI1/2 '
      N = 0
      DO IKM1 = 1,NKMMAX
         DO IKM2 = 1,NKMMAX
            DO ILA = 1,NLABIMAX
               DO M = -1, + 1
                  IF ( ABS(AMEBI1(IKM1,IKM2,ILA,M)).GT.1.0D-10 )
     &                 N = N + 1
                  D = AMEBI1(IKM1,IKM2,ILA,M) - AMEBI2(IKM2,IKM1,ILA,-M)
                  IF ( ABS(D).GT.1.D-10 ) WRITE (6,'(A,4I3,2E15.5)')
     &                  ' trouble in <SCFBIAME>',IKM1,IKM2,ILA,M,
     &                 AMEBI1(IKM1,IKM2,ILA,M),D
               END DO
            END DO
         END DO
      END DO
      WRITE (6,*) ' <SCFBIAME> non-0  AMEBI1''s :',N
C
      END
