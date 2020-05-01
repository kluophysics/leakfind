C*==ame_rsig.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AME_RSIG
C   ********************************************************************
C   *                                                                  *
C   *   calculate angular matrix elements occuring for the             *
C   *           full conductivity tensor                               *
C   *                                                                  *
C   *  A_RSIG = < KAP1,MUE1| r_lam * SIG_lam' |KAP2,MUE2>              *
C   *                                                                  *
C   *  14/11/08  HE                                                    *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI,SQRT_2
      USE MOD_ANGMOM,ONLY:A_RSIG,A_SIGMA,CGC,NLMAX,NKMPMAX
      IMPLICIT NONE
C*--AME_RSIG15
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 AMP,APM,CGNT,F,NORM_Y1M,RSUM,SIG
      REAL*8 GAUNT_CYLM
      INTEGER IA_ERR,IKM1,IKM2,IMS1,IMS2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,
     &        L2,LAM1,LAM2,M1,M2,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (A_RSIG(NKMPMAX,NKMPMAX,-1:+1,-1:+1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: AME_RSIG -> AME_RSIG'
C
      A_RSIG(1:NKMPMAX,1:NKMPMAX,-1:+1,-1:+1) = 0D0
C
      NK = 2*NLMAX
      NORM_Y1M = DSQRT(4D0*PI/3D0)
C
C=======================================================================
      DO LAM1 = -1, + 1
         DO LAM2 = -1, + 1
C
C ----------------------------------------------------------------------
            IKM1 = 0
            DO K1 = 1,NK
               L1 = K1/2
               IF ( MOD(K1,2).EQ.0 ) THEN
                  KAP1 = L1
               ELSE
                  KAP1 = -L1 - 1
               END IF
               J1P05 = ABS(KAP1)
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
                     J2P05 = ABS(KAP2)
C
                     DO MUE2M05 = -J2P05,J2P05 - 1
                        IKM2 = IKM2 + 1
C
C ----------------------------------------------------------------------
                        RSUM = 0D0
C
                        DO IMS1 = 1,2
                           M1 = MUE1M05 - IMS1 + 2
C
                           DO IMS2 = 1,2
                              M2 = MUE2M05 - IMS2 + 2
C
                              CGNT = GAUNT_CYLM(L1,M1,1,LAM1,L2,M2)
                              SIG = A_SIGMA(IMS1,LAM2,IMS2)
C
                              RSUM = RSUM + CGC(IKM1,IMS1)
     &                               *CGC(IKM2,IMS2)*CGNT*SIG
C
                           END DO
                        END DO
C
                        A_RSIG(IKM1,IKM2,LAM1,LAM2) = NORM_Y1M*RSUM
C ----------------------------------------------------------------------
C
                     END DO
                  END DO
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END DO
      END DO
C
C ======================================================================
C
      F = SQRT_2*NORM_Y1M
C
C ======================================================================
C                       check for (-,+) and (+,-)
C=======================================================================
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
         ELSE
            KAP1 = -L1 - 1
         END IF
         J1P05 = ABS(KAP1)
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
               J2P05 = ABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
C
C ----------------------------------------------------------------------
                  IMS1 = 2
                  IMS2 = 1
                  M1 = MUE1M05 - IMS1 + 2
                  M2 = MUE2M05 - IMS2 + 2
                  CGNT = GAUNT_CYLM(L1,M1,1,-1,L2,M2)
C
                  AMP = -F*CGC(IKM1,IMS1)*CGC(IKM2,IMS2)*CGNT
C
C
                  IMS1 = 1
                  IMS2 = 2
                  M1 = MUE1M05 - IMS1 + 2
                  M2 = MUE2M05 - IMS2 + 2
                  CGNT = GAUNT_CYLM(L1,M1,1,+1,L2,M2)
C
                  APM = F*CGC(IKM1,IMS1)*CGC(IKM2,IMS2)*CGNT
C
                  IF ( ABS(A_RSIG(IKM1,IKM2,-1,+1)-AMP).GT.1D-8 .OR. 
     &                 ABS(A_RSIG(IKM1,IKM2,+1,-1)-APM).GT.1D-8 ) THEN
                     WRITE (6,*) ' TROUBLE in <AME_RSIG>'
                     WRITE (6,*) ' IKM1 = ',IKM1,'    IKM2 = ',IKM2
                     WRITE (6,*) ' A(-,+) ',A_RSIG(IKM1,IKM2,-1,+1),AMP
                     WRITE (6,*) ' A(+,-) ',A_RSIG(IKM1,IKM2,+1,-1),APM
                  END IF
C ----------------------------------------------------------------------
C
               END DO
            END DO
         END DO
      END DO
C
C ======================================================================
C
      END
