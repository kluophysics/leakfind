C*==dmft_spttma.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C     Within KKR we don't bother ourself with making SE full
C     (hermitian+antihermitian) we will get it by making
C     Kramers-Kronig transform in the driving routine.
C
      SUBROUTINE DMFT_SPTTMA(NREP,NRR,X,GF0,NZEROG,SE_AC,V,IPRT,TOL)
C
      USE MOD_CONSTANTS,ONLY:C0,CI,PI
      IMPLICIT NONE
C*--DMFT_SPTTMA10
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL IPRT
      INTEGER NREP,NRR
      REAL*8 TOL
      COMPLEX*16 GF0(NREP,NREP,NRR),SE_AC(NRR,NREP,NREP)
      INTEGER NZEROG(NREP,NREP)
      REAL*8 V(NREP**2,NREP**2),X(NRR)
C
C Local variables
C
      COMPLEX*16 CTMP,CWORK(:,:),CWORK1(:),PHI_AC(:,:,:),PHI_AC3(:,:,:),
     &           PHI_FC(:,:,:),PHI_FC1(:,:,:),PHI_FC2(:,:,:),
     &           PHI_HC(:,:,:)
      REAL*8 DX,TIME1,TIME2,TMP
      INTEGER I,IE,IND(:,:),INDINV(:,:),INFO,IPVT(:),IR,IR1,IR2,IR3,IR4,
     &        IRP,IX,IZ,K,LAMBDA,LAMBDA1,LAMBDA2,NREP2,NZERO2(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CWORK,CWORK1
      ALLOCATABLE PHI_AC,PHI_HC,PHI_FC,PHI_FC1,PHI_FC2
      ALLOCATABLE PHI_AC3
      ALLOCATABLE IPVT
      ALLOCATABLE IND,INDINV,NZERO2
C
      WRITE (6,99001)
C
      NREP2 = NREP**2
C
      ALLOCATE (IND(NREP,NREP))
      ALLOCATE (INDINV(NREP2,2))
C
C=======================================================
C  Find the index of the Fermi-level:
C=======================================================
C
      DO IE = 1,NRR
         IF ( X(IE).LT.0D0 ) IZ = IE
      END DO
      IF ( X(IZ+1).LT.(DABS(X(IZ))) ) IZ = IZ + 1
      DX = (X(NRR)-X(1))/(NRR-1D0)
C
C
      IF ( IPRT ) CALL CMATSTRUCT('TMA-SOLVER: input GF-matrix',
     &                            GF0(1,1,IZ+10),NREP,NREP,0,0,0,TOL,6)
      IF ( IPRT ) CALL CPU_TIME(TIME1)
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Create double index:
C=======================================================
      K = 0
      IND = 0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            K = K + 1
            IND(IR,IRP) = K
            INDINV(K,1) = IR
            INDINV(K,2) = IRP
         END DO
      END DO
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Calculate antihermitian complex phi:
C=======================================================
      CALL CPU_TIME(TIME1)
C
      ALLOCATE (PHI_AC(NRR,NREP2,NREP2),PHI_HC(NRR,NREP2,NREP2))
      ALLOCATE (CWORK1(NRR))
      PHI_AC(1:NRR,1:NREP2,1:NREP2) = C0
      PHI_HC(1:NRR,1:NREP2,1:NREP2) = C0
      CWORK1(1:NRR) = C0
C
      DO LAMBDA1 = 1,NREP2
         DO LAMBDA2 = LAMBDA1,NREP2
            IR1 = INDINV(LAMBDA1,1)
            IR2 = INDINV(LAMBDA1,2)
            IR3 = INDINV(LAMBDA2,1)
            IR4 = INDINV(LAMBDA2,2)
            IF ( NZEROG(IR1,IR3).EQ.1 .AND. NZEROG(IR2,IR4).EQ.1 ) THEN
               DO IE = 1,NRR
                  DO I = 1,NRR
                     IX = IE - I + IZ
                     IF ( IX.GE.1 .AND. IX.LE.NRR ) THEN
                        CWORK1(I) = GF0(IR1,IR3,I)*GF0(IR2,IR4,IX)
                     ELSE
                        CWORK1(I) = C0
                     END IF
                  END DO
                  CALL DMFT_CINTEGRAL(PHI_AC(IE,LAMBDA1,LAMBDA2),DX,
     &                                CWORK1,NRR,IZ,IE)
               END DO
            END IF
         END DO
      END DO
C
      PHI_AC(IZ,1:NREP2,1:NREP2) = C0
      DO LAMBDA1 = 1,NREP2
         DO LAMBDA2 = LAMBDA1,NREP2
            IR1 = INDINV(LAMBDA1,1)
            IR2 = INDINV(LAMBDA1,2)
            IR3 = INDINV(LAMBDA2,1)
            IR4 = INDINV(LAMBDA2,2)
            IF ( NZEROG(IR1,IR3).EQ.1 .AND. NZEROG(IR2,IR4).EQ.1 )
     &           CALL DMFT_HILBERTT(NRR,X,PHI_AC(1,LAMBDA1,LAMBDA2),
     &                              CWORK1,PHI_HC(1,LAMBDA1,LAMBDA2),DX)
         END DO
      END DO
      DEALLOCATE (CWORK1)
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Printout for checking:
C=======================================================
      IF ( IPRT ) THEN
         OPEN (1,FILE='phi_ah_diag_re.dat')
         DO IE = 1,NRR
            WRITE (1,'(1000E13.3)') X(IE),
     &                              ((DREAL(PHI_AC(IE,IR,IRP)),IRP=IR,
     &                              IR),IR=1,NREP2)
         END DO
         CLOSE (1)
         OPEN (1,FILE='phi_he_diag_im.dat')
         DO IE = 1,NRR
            WRITE (1,'(1000E13.3)') X(IE),
     &                              ((DIMAG(PHI_HC(IE,IR,IRP)),IRP=IR,
     &                              IR),IR=1,NREP2)
         END DO
         CLOSE (1)
         CALL CPU_TIME(TIME2)
         WRITE (6,'(/,10X,A,1X,F7.3,1X,A)') 'time for PHI-matrix:',
     &          TIME2 - TIME1,'[sec]'
         ALLOCATE (CWORK(NREP2,NREP2))
         CWORK(1:NREP2,1:NREP2) = PHI_AC(IZ+10,1:NREP2,1:NREP2)
         CALL CMATSTRUCT('TMA-SOLVER: PHI_AC',CWORK,NREP2/4,NREP2,0,0,0,
     &                   TOL,6)
         DEALLOCATE (CWORK)
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Calculate full complex phi:
C=======================================================
      ALLOCATE (PHI_FC(NRR,NREP2,NREP2))
C
      DO IR1 = 1,NREP
         DO IR2 = 1,NREP
            DO IR3 = 1,NREP
               DO IR4 = 1,NREP
                  LAMBDA1 = IND(IR1,IR2)
                  LAMBDA2 = IND(IR3,IR4)
                  IF ( LAMBDA2.GE.LAMBDA1 ) THEN
                     DO IE = 1,NRR
                        PHI_FC(IE,IND(IR1,IR2),IND(IR3,IR4))
     &                     = (PHI_AC(IE,LAMBDA1,LAMBDA2)
     &                     +PHI_HC(IE,LAMBDA1,LAMBDA2))/(PI*CI)
                     END DO
                  ELSE
                     DO IE = 1,NRR
                        PHI_FC(IE,IND(IR1,IR2),IND(IR3,IR4))
     &                     = (PHI_AC(IE,LAMBDA2,LAMBDA1)
     &                     +PHI_HC(IE,LAMBDA2,LAMBDA1))/(PI*CI)
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
C
      DEALLOCATE (PHI_AC)
      DEALLOCATE (PHI_HC)
C
C printout:
      IF ( IPRT ) THEN
         ALLOCATE (CWORK(NREP2,NREP2))
         CWORK(1:NREP2,1:NREP2) = PHI_FC(IZ+10,1:NREP2,1:NREP2)
         CALL CMATSTRUCT('TMA-SOLVER: PHI_FULL',CWORK,NREP2/4,NREP2,0,0,
     &                   0,TOL,6)
         DEALLOCATE (CWORK)
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Calculate full complex phi, ready to invert:
C=======================================================
      IF ( IPRT ) CALL CPU_TIME(TIME1)
C
      ALLOCATE (PHI_FC1(NRR,NREP2,NREP2))
      PHI_FC1 = C0
      DO LAMBDA1 = 1,NREP2
         DO LAMBDA2 = 1,NREP2
            DO LAMBDA = 1,NREP2
               IF ( DABS(V(LAMBDA1,LAMBDA)).GT.TOL ) THEN
                  DO IE = 1,NRR
                     PHI_FC1(IE,LAMBDA1,LAMBDA2)
     &                  = PHI_FC1(IE,LAMBDA1,LAMBDA2)
     &                  - V(LAMBDA1,LAMBDA)*PHI_FC(IE,LAMBDA,LAMBDA2)
                  END DO
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (PHI_FC)
C
      DO LAMBDA = 1,NREP2
         PHI_FC1(1:NRR,LAMBDA,LAMBDA) = PHI_FC1(1:NRR,LAMBDA,LAMBDA)
     &                                  + 1D0
      END DO
C
      IF ( IPRT ) THEN
         CALL CPU_TIME(TIME2)
         IF ( IPRT ) WRITE (6,'(/,10X,A,1X,F7.3,1X,A)')
     &                       'time for 1-VxPHI-matrix:',TIME2 - TIME1,
     &                      '[sec]'
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Invert full complex phi:
C=======================================================
      IF ( IPRT ) CALL CPU_TIME(TIME1)
C
      ALLOCATE (CWORK(NREP2,NREP2),CWORK1(NREP2))
      ALLOCATE (IPVT(NREP2))
      DO IE = 1,NRR
         DO I = 1,NREP2
            DO K = 1,NREP2
               CWORK(K,I) = PHI_FC1(IE,K,I)
            END DO
         END DO
         CALL ZGETRF(NREP2,NREP2,CWORK,NREP2,IPVT,INFO)
         CALL ZGETRI(NREP2,CWORK,NREP2,IPVT,CWORK1,NREP2,INFO)
         DO I = 1,NREP2
            DO K = 1,NREP2
               PHI_FC1(IE,K,I) = CWORK(K,I)
            END DO
         END DO
      END DO
      DEALLOCATE (IPVT,CWORK,CWORK1)
C
C printout:
      IF ( IPRT ) THEN
         CALL CPU_TIME(TIME2)
         WRITE (6,'(/,10X,A,1X,F7.3,1X,A)') 'time for inversion:',
     &          TIME2 - TIME1,'[sec]'
C
C
         ALLOCATE (CWORK(NREP2,NREP2))
         CWORK(1:NREP2,1:NREP2) = PHI_FC1(IZ+10,1:NREP2,1:NREP2)
         CALL CMATSTRUCT('TMA-SOLVER: (1-V*PHI_FC1)^{-1}',CWORK,NREP2/4,
     &                   NREP2,0,0,0,TOL,6)
         DEALLOCATE (CWORK)
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Calculate full complex T:
C=======================================================
      IF ( IPRT ) CALL CPU_TIME(TIME1)
C
      ALLOCATE (PHI_FC2(NRR,NREP2,NREP2))
      PHI_FC2(1:NRR,1:NREP2,1:NREP2) = C0
      DO LAMBDA1 = 1,NREP2
         DO LAMBDA = 1,NREP2
            IF ( DABS(V(LAMBDA,LAMBDA1)).GT.TOL ) THEN
               DO K = LAMBDA1,NREP2
                  DO IE = 1,NRR
                     PHI_FC2(IE,K,LAMBDA1) = PHI_FC2(IE,K,LAMBDA1)
     &                  + PHI_FC1(IE,K,LAMBDA)*V(LAMBDA,LAMBDA1)
                  END DO
                  PHI_FC2(1:NRR,LAMBDA1,K) = PHI_FC2(1:NRR,K,LAMBDA1)
               END DO
            END IF
         END DO
      END DO
      DEALLOCATE (PHI_FC1)
C
      PHI_FC2(IZ,1:NREP2,1:NREP2) = V(1:NREP2,1:NREP2) ! stress static limit
C
C
C trying to get rid of statics (1st order contributions) (testing!!!)
C$$$      DO IE=1,NRR
C$$$         PHI_FC2(IE,:,:)=PHI_FC2(IE,:,:)-V(:,:)
C$$$      ENDDO
C
C
C antihermitise T:
      PHI_FC2(1:NRR,1:NREP2,1:NREP2) = DCMPLX(0D0,DIMAG(PHI_FC2(1:NRR,1:
     &                                 NREP2,1:NREP2)))
C printout:
      IF ( IPRT ) THEN
         CALL CPU_TIME(TIME2)
         WRITE (6,'(/,10X,A,1X,F7.3,1X,A)') 'time for AH T-matrix:',
     &          TIME2 - TIME1,'[sec]'
         ALLOCATE (CWORK(NREP2,NREP2))
         CWORK(1:NREP2,1:NREP2) = PHI_FC2(IZ+10,1:NREP2,1:NREP2)
         CALL CMATSTRUCT('TMA-SOLVER: 1/2 (1-V*PHI_FC1)^{-1}*V',CWORK,
     &                   (NREP2)/2,NREP2,0,0,0,TOL,6)
         DEALLOCATE (CWORK)
      END IF
C
C make Tijkl=Tijkl-Tijlk and mark non-zero elements:
C
      ALLOCATE (NZERO2(NREP2,NREP2))
      ALLOCATE (PHI_AC3(NRR,NREP2,NREP2))
      PHI_AC3(1:NRR,1:NREP2,1:NREP2) = C0
      NZERO2(1:NREP2,1:NREP2) = 0
      DO IR1 = 1,NREP
         DO IR2 = 1,NREP
            DO IR3 = 1,NREP
               DO IR4 = 1,NREP
                  TMP = 0D0
                  DO IE = 1,NRR
                     PHI_AC3(IE,IND(IR1,IR2),IND(IR3,IR4))
     &                  = PHI_FC2(IE,IND(IR1,IR2),IND(IR3,IR4))
     &                  - PHI_FC2(IE,IND(IR1,IR2),IND(IR4,IR3))
                     TMP = TMP + 
     &                     CDABS(PHI_AC3(IE,IND(IR1,IR2),IND(IR3,IR4)))
                  END DO
                  IF ( TMP.GE.TOL ) NZERO2(IND(IR1,IR2),IND(IR3,IR4))
     &                 = 1
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (PHI_FC2)
C
C
C printout:
      IF ( IPRT ) THEN
         ALLOCATE (CWORK(NREP2,NREP2))
         CWORK(1:NREP2,1:NREP2) = PHI_AC3(IZ+10,1:NREP2,1:NREP2)
         CALL CMATSTRUCT('TMA-SOLVER: 1/2 T-MATRIX:',CWORK,(NREP2)/2,
     &                   NREP2,0,0,0,TOL,6)
         DEALLOCATE (CWORK)
         OPEN (1,FILE='TM_AH_DIAG_IM.DAT')
         DO IE = 1,NRR
            WRITE (1,'(1000E13.3)') X(IE),
     &                              ((DIMAG(PHI_AC3(IE,IR,IRP)),IRP=IR,
     &                              IR),IR=1,NREP2)
         END DO
         CLOSE (1)
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C  Calculate the antihermitian SELF-ENERGY:
C=======================================================
      IF ( IPRT ) CALL CPU_TIME(TIME1)
C
      ALLOCATE (CWORK1(NRR))
C
      SE_AC(1:NRR,1:NREP,1:NREP) = C0
      DO IR1 = 1,NREP
         DO IR2 = IR1,NREP
            DO IR3 = 1,NREP
               DO IR4 = 1,NREP
C
                  IF ( NZEROG(IR4,IR3).EQ.1 .AND. 
     &                 NZERO2(IND(IR1,IR3),IND(IR2,IR4)).EQ.1 ) THEN
                     DO IE = 1,NRR
                        DO I = 1,NRR
                           IX = IE + I - IZ
                           IF ( IX.GE.1 .AND. IX.LE.NRR ) THEN
                              CWORK1(I)
     &                           = PHI_AC3(IX,IND(IR1,IR3),IND(IR2,IR4))
     &                           *GF0(IR4,IR3,I)
                           ELSE
                              CWORK1(I) = C0
                           END IF
                        END DO  !I
C
                        IX = 2*IZ - IE
                        IF ( IX.LT.1 ) IX = 1
                        IF ( IX.GT.NRR ) IX = NRR
C
                        CALL DMFT_CINTEGRAL(CTMP,DX,CWORK1,NRR,IX,IZ)
                        IF ( IX.EQ.IZ ) CTMP = C0
                        SE_AC(IE,IR1,IR2) = SE_AC(IE,IR1,IR2) + CTMP
C
                     END DO     ! IE
                  END IF        ! NZEROG
C
               END DO      ! IR4
            END DO         ! IR3
            SE_AC(1:NRR,IR2,IR1) = SE_AC(1:NRR,IR1,IR2)
         END DO                 ! IR2
      END DO  ! IR1
C
C
      SE_AC(1:NRR,1:NREP,1:NREP) = SE_AC(1:NRR,1:NREP,1:NREP)/(PI*CI)
C
      SE_AC(1:NRR,1:NREP,1:NREP) = DCMPLX(0D0,DIMAG(SE_AC(1:NRR,1:NREP,1
     &                             :NREP)))
C
C
      DEALLOCATE (CWORK1)
      DEALLOCATE (NZERO2)
      DEALLOCATE (PHI_AC3)
      DEALLOCATE (IND,INDINV)
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Check printout:
C=======================================================
      IF ( IPRT ) THEN
         ALLOCATE (CWORK(NREP,NREP))
         DO I = 1,NREP
            DO K = 1,NREP
               CWORK(K,I) = SE_AC(IZ+10,K,I)
            END DO
         END DO
         CALL CMATSTRUCT('antihermitian SE-matrix',CWORK,NREP,NREP,0,0,
     &                   0,TOL,6)
         DEALLOCATE (CWORK)
         CALL CPU_TIME(TIME2)
         WRITE (6,'(/,10X,A,1X,F7.3,1X,A)')
     &           'time for AH self-energy-matrix:',TIME2 - TIME1,'[sec]'
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Set SELF-ENERGY to ZERO on the Fermi-level:
C=======================================================
      DO IR = 1,NREP
         DO IRP = 1,NREP
            CTMP = SE_AC(IZ,IR,IRP)
            DO IE = 1,NRR
               SE_AC(IE,IR,IRP) = SE_AC(IE,IR,IRP) - CTMP
            END DO
         END DO
      END DO
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C
      IF ( IPRT ) WRITE (6,'(10X,A)') 'TMA_SOLVER - FINISHED'
C
99001 FORMAT (/,1X,79('*'),/,35X,'<TMA_SOLVER>',/,1X,79('*'),/)
      END
