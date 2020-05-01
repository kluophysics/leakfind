C*==amey1m.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMEY1M
C   ********************************************************************
C   *                                                                  *
C   *   angular matrix element needed for electric dipole moment       *
C   *                                                                  *
C   *   RELATIVISTIC                                                   *
C   *                                                                  *
C   *   A_Y1M  (IKM1,IKM2,MR) = < KAP1,MUE1 | Y[1,MR] |  KAP2,MUE2>    *
C   *   A_Y1M_F(IKM1,IKM2,MR) = <-KAP1,MUE1 | Y[1,MR] | -KAP2,MUE2>    *
C   *                         = A_Y1M(IKM1,IKM2,MR)                    *
C   *   only A_Y1M is stored                                           *
C   *                                                                  *
C   *                                                                  *
C   *   NON-RELATIVISTIC                                               *
C   *                                                                  *
C   *   A_Y1M(LM1,LM2,MR) = < L1,M1 | Y[1,MR] | L2,M2>                 *
C   *                     = Gaunt coefficient                          *
C   *                                                                  *
C   *   Y_1m : REAL sperical harmonics for l=1, m=-1,0,+1              *
C   *                                                                  *
C   *   NOTE:  A_Y1M and A_Y1M_F are declared to be COMPLEX*16         *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:A_Y1M,CGC,NLMAX,NKMMAX,NLMMAX,NKMPMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:CI,C0
      IMPLICIT NONE
C*--AMEY1M30
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      COMPLEX*16 AMIN,APLS,A_Y1M_F(:,:,:)
      REAL*8 GAUNT_CYLM,GAUNT_RYLM
      INTEGER IKM1,IKM2,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,L2,
     &        L3,LB1,LB2,LM1,LM2,LR,M1,M2,M3,MR,MSM05,MUE1M05,MUE2M05,NK
      REAL*8 SQRT_2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE A_Y1M_F
C
C=======================================================================
C                            RELATIVISTIC
C=======================================================================
      IF ( IREL.EQ.3 ) THEN
C
         NK = 2*NLMAX
C
         IF ( ALLOCATED(A_Y1M) ) DEALLOCATE (A_Y1M)
C
         ALLOCATE (A_Y1M(NKMPMAX,NKMPMAX,-1:+1))
         ALLOCATE (A_Y1M_F(NKMPMAX,NKMPMAX,-1:+1))
C
         A_Y1M(1:NKMPMAX,1:NKMPMAX,-1:+1) = C0
         A_Y1M_F(1:NKMPMAX,1:NKMPMAX,-1:+1) = C0
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
                     LR = 1
                     DO MSM05 = -1,0
                        M1 = MUE1M05 - MSM05
                        M2 = MUE2M05 - MSM05
                        MR = M1 - M2
C
C                G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]
C
                        IF ( MR.GE.-1 .AND. MR.LE.1 ) THEN
C
                           A_Y1M(IKM1,IKM2,MR) = A_Y1M(IKM1,IKM2,MR)
     &                        + CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                        *GAUNT_CYLM(L1,M1,LR,MR,L2,M2)
C
                           A_Y1M_F(IKM1,IKM2,MR) = A_Y1M_F(IKM1,IKM2,MR)
     &                        + CGC(IMKM1,MSM05+2)*CGC(IMKM2,MSM05+2)
     &                        *GAUNT_CYLM(LB1,M1,LR,MR,LB2,M2)
C
                        END IF
C
                     END DO
C ----------------------------------------------------------------------
                  END DO
               END DO
            END DO
         END DO
C
C ----------------------------------------------------------------------
C               change to REAL spherical harmonics
C ----------------------------------------------------------------------
C
         SQRT_2 = SQRT(2D0)
C
         DO IKM1 = 1,NKMMAX
            DO IKM2 = 1,NKMMAX
C
               DO MR = -1, + 1
                  IF ( ABS(A_Y1M(IKM1,IKM2,MR)-A_Y1M_F(IKM1,IKM2,MR))
     &                 .GT.1D-10 ) THEN
                     WRITE (6,*) '<AMEY1M>: ',IKM1,IKM2,MR,
     &                           A_Y1M(IKM1,IKM2,MR),
     &                           A_Y1M_F(IKM1,IKM2,MR)
                     STOP
                  END IF
               END DO
C
               AMIN = A_Y1M(IKM1,IKM2,-1)
               APLS = A_Y1M(IKM1,IKM2,+1)
C
               A_Y1M(IKM1,IKM2,-1) = (AMIN+APLS)*CI/SQRT_2
               A_Y1M(IKM1,IKM2,+1) = (AMIN-APLS)/SQRT_2
C
            END DO
         END DO
C
         DEALLOCATE (A_Y1M_F)
C
C=======================================================================
C                          NON-RELATIVISTIC
C=======================================================================
C
      ELSE
C
         IF ( ALLOCATED(A_Y1M) ) DEALLOCATE (A_Y1M)
         ALLOCATE (A_Y1M(NLMMAX,NLMMAX,-1:+1))
C
         A_Y1M(1:NLMMAX,1:NLMMAX,-1:+1) = C0
C
         LM1 = 0
         DO L1 = 0,(NLMAX-1)
            DO M1 = -L1,L1
               LM1 = LM1 + 1
C
               LM2 = 0
               DO L2 = 0,(NLMAX-1)
                  DO M2 = -L2,L2
                     LM2 = LM2 + 1
C
                     L3 = 1
                     DO M3 = -L3, + L3
                        A_Y1M(LM1,LM2,M3)
     &                     = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
                     END DO
C
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
C=======================================================================
C
      END
