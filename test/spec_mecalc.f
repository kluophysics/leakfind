C*==spec_mecalc.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_MECALC(RMATO,RMATS,ZMAT0,ZMAT0INC,TSSTS,MSSTS,P,
     &                       GAMMA1,TAUMT,UMAT_VT_PES,CPAPROJ,CPAATOM,
     &                       VLM,AA,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,
     &                       AMAT2Y,AMAT2V,AMAT2,AMAT1T,AMAT2T,PFI,
     &                       CLIGHT,OMHAR,IP,IBLOCH,EHIGH,ELOW,USEEULER,
     &                       JE,ATA,LAYS,NATL,INITIAL,FINAL,IRSTATE,
     &                       NOUT1,RADMESHM)
C
C   ********************************************************************
C   *                                                                  *
C   *  driver to calculations of matrix elemnts for spec               *
C   *                                                                  *
C   ********************************************************************
C
C
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,NVFTPHOMAX,NTPHOMAX,NFULLPOT,
     &    MQK,XMAXE,EPS12
      USE MOD_CONSTANTS,ONLY:C0,C1,CI
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_THERMAL,ONLY:X_VFT,NVFO_Q,IVFT_VFOQ
      USE MOD_CPA,ONLY:NCPA
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_SPEC_POTLM,ONLY:ALPHA,BB
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SPEC_WAVE,ONLY:WFFM,WFFMZ,WFFXM,WFFYM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATA,FINAL,IBLOCH,INITIAL,IP,IRSTATE,JE,LAYS,NOUT1,USEEULER
      REAL*8 CLIGHT,OMHAR
      COMPLEX*16 EHIGH,ELOW,PFI
      COMPLEX*16 AA(3),AMAT1(XMAXE,4,NFULLPOT,2),
     &           AMAT1T(XMAXE,4,NFULLPOT,2),AMAT2(XMAXE,4,NFULLPOT),
     &           AMAT2T(XMAXE,4,NFULLPOT),
     &           CPAPROJ(LAYSM,NATLM,MQD,MQD,2,NVFTPHOMAX),
     &           GAMMA1(LAYSM,NATLM,3,MQD,MQD),
     &           MSSTS(NKMMAX,NKMMAX,NTMAX,2),P(2),
     &           RMATO(MQD,MQD,LAYSM,NATLM),RMATS(MQD,MQD,LAYSM,NATLM),
     &           TAUMT(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX+1),
     &           TSSTS(NKMMAX,NKMMAX,NTMAX,2),
     &           UMAT_VT_PES(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX),
     &           ZMAT0(MQD,MQD,LAYSM,NATLM),
     &           ZMAT0INC(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX+1)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        CPAATOM(LAYSM,NATLM),NATL(LAYSM)
      REAL*8 RADMESHM(NRMAX,NATLM,LAYSM)
C
C Local variables
C
      INTEGER ATOM,G1(:),G2(:),I,I1,I2,ICV,IO,IQ,ITMP,ITS,J,LAY,MAXG,
     &        PHSTATE,PMS,PMS1(:),PMS2(:)
      COMPLEX*16 C2(:,:),C3(:,:),C4(:,:),C5(:,:),CM1(:,:),EE,
     &           RDIP1M(:,:,:,:,:),RDIP2M(:,:,:,:),RMATOC(:,:,:,:,:),
     &           RMATOTP(:,:,:,:,:),RMATSC(:,:,:,:,:),RMATSTP(:,:,:,:,:)
     &           ,ZMAT0TP(:,:,:,:,:)
      REAL*8 CONZ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RDIP1M,RDIP2M,G1,G2,CM1,C2,C3,C4,C5
      ALLOCATABLE PMS1,PMS2
      ALLOCATABLE RMATOC,RMATOTP,RMATSC,RMATSTP,ZMAT0TP
C
C=======================================================================
C                  Init matrix elements
C=======================================================================
C
      ALLOCATE (RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX))
      ALLOCATE (RDIP2M(NFULLPOT,NRMAX,2,NTPHOMAX))
C
      RDIP1M = C0
      RDIP2M = C0
C
      ALLOCATE (G1(MQK),G2(MQK))
      G1 = 0
      G2 = 0
      MAXG = 0
      ALLOCATE (CM1(MQD,MQD),C2(MQD,MQD))
      ALLOCATE (C3(MQD,MQD),C4(MQD,MQD),C5(MQD,MQD))
      CM1 = C0
      C2 = C0
      C3 = C0
      C4 = C0
      C5 = C0
      ALLOCATE (PMS1(MQK),PMS2(MQK))
      ALLOCATE (RMATOTP(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (RMATSTP(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (RMATOC(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (RMATSC(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (ZMAT0TP(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      RMATOTP = C0
      RMATSTP = C0
      RMATOC = C0
      RMATSC = C0
      ZMAT0TP = C0
      RMATO = C0
      RMATS = C0
      ZMAT0 = C0
      ZMAT0INC = C0
C
C=======================================================================
C        calculation of angular matrix elements
C=======================================================================
C
C
      CALL DETMAT(AA,BB,VLM,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,
     &            AMAT2V,AMAT2,AMAT1T,AMAT2T,NATL,LAYS)
C=======================================================================
C
      LOOP_LAY:DO LAY = 1,LAYS
         LOOP_ATOM:DO ATOM = 1,NATL(LAY)
            IQ = IQ_SPR_LAYAT(ATOM,LAY)
C
            LOOP_IO:DO IO = 1,NVFO_Q(IQ)
               ICV = IO
               ITMP = IO
               IF ( THERMAL_VIBRA_FLUCT ) ITMP = 1
C
               IF ( NCPA.GT.0 ) THEN
                  ITS = IVFT_VFOQ(IO,IQ)
                  CONZ = X_VFT(ITS)
               ELSE IF ( NCPA.EQ.0 ) THEN
                  CONZ = 1.0D0
               END IF
C
               IF ( NCPA.GT.0 .AND. IP.GT.0 ) WRITE (6,*) 'CONZ',CONZ
C
C=======================================================================
C        determine radial parts of dipole operator
C=======================================================================
C
               CALL RDIPOLEM(LAY,ATOM,VLM,RDIP1M,RDIP2M,CLIGHT,OMHAR,
     &                       ALPHA(LAY,ATOM),ITMP)
C
C
C=======================================================================
C                   prepare single site wave functions
C=======================================================================
C
               DO PHSTATE = 1,2
                  IF ( PHSTATE.EQ.1 ) EE = ELOW
                  IF ( PHSTATE.EQ.2 ) EE = EHIGH
                  CALL FULLF(LAY,ATOM,ITMP,PHSTATE,EE,CLIGHT,USEEULER,
     &                       IBLOCH,JE,VLM,MAXG,G1,G2,
     &                       TSSTS(1,1,1,PHSTATE),MSSTS(1,1,1,PHSTATE),
     &                       P(PHSTATE))
               END DO
C=======================================================================
               CM1(1:MQD,1:MQD) = C0
               FORALL(I=1:MQD)CM1(I,I) = C1
               IF ( CPAATOM(LAY,ATOM).EQ.1 .AND. ATA.EQ.0 )
     &              CM1(1:MQD,1:MQD) = CPAPROJ(LAY,ATOM,1:MQD,1:MQD,1,IO
     &                                 )
               C2(1:MQD,1:MQD) = GAMMA1(LAY,ATOM,1,1:MQD,1:MQD)
C
               PMS1 = 0
               PMS2 = 0
               PMS = 0
               DO I = 1,MQD
                  DO J = 1,MQD
                     IF ( CDABS(CM1(I,J)).GT.1.D-16 ) THEN
                        PMS = PMS + 1
                        PMS1(PMS) = I
                        PMS2(PMS) = J
                     END IF
                  END DO
               END DO
C
               IF ( THERMAL_VIBRA_FLUCT ) THEN
                  C3(1:MQD,1:MQD) = CM1(1:MQD,1:MQD)
                  C4(1:MQD,1:MQD) = UMAT_VT_PES(LAY,ATOM,1:MQD,1:MQD,IO)
                  C5 = TRANSPOSE(C4)
                  CM1 = MATMUL(C5,C3)
               END IF
C
C
C=======================================================================
C                   calculate regular MELE
C=======================================================================
C
               IF ( IBLOCH.EQ.4 ) THEN
                  IF ( NCPA.EQ.0 ) THEN
                     CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,
     &                           AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,WFFYM,
     &                           LAY,ATOM,RMATOTP,RMATSTP,
     &                           ALPHA(LAY,ATOM),OMHAR,FINAL,IO,INITIAL,
     &                           AMAT1T,AMAT2T,WFFM,PFI,NCPA,ATA,ICV)
                  ELSE IF ( NCPA.GT.0 ) THEN
                     IF ( ATA.EQ.1 ) THEN
                        IF ( THERMAL_VIBRA_FLUCT ) THEN
                           CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                                 AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,
     &                                 WFFXM,WFFYM,LAY,ATOM,RMATOTP,
     &                                 RMATSTP,ALPHA(LAY,ATOM),OMHAR,
     &                                 FINAL,1,INITIAL,AMAT1T,AMAT2T,
     &                                 WFFM,PFI,NCPA,ATA,ICV)
                        ELSE
                           CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                                 AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,
     &                                 WFFXM,WFFYM,LAY,ATOM,RMATOTP,
     &                                 RMATSTP,ALPHA(LAY,ATOM),OMHAR,
     &                                 FINAL,IO,INITIAL,AMAT1T,AMAT2T,
     &                                 WFFM,PFI,NCPA,ATA,ICV)
                        END IF
                     ELSE IF ( ATA.EQ.0 ) THEN
                        IF ( CPAATOM(LAY,ATOM).EQ.1 ) THEN
                           IF ( THERMAL_VIBRA_FLUCT ) THEN
                              CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                           AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,
     &                           WFFYM,LAY,ATOM,RMATOTP,RMATSTP,
     &                           ALPHA(LAY,ATOM),OMHAR,FINAL,1,INITIAL,
     &                           AMAT1T,AMAT2T,WFFMZ,PFI,NCPA,ATA,ICV)
                           ELSE
                              CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                           AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,
     &                           WFFYM,LAY,ATOM,RMATOTP,RMATSTP,
     &                           ALPHA(LAY,ATOM),OMHAR,FINAL,IO,INITIAL,
     &                           AMAT1T,AMAT2T,WFFMZ,PFI,NCPA,ATA,ICV)
                           END IF
                        ELSE IF ( CPAATOM(LAY,ATOM).EQ.0 ) THEN
                           CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                                 AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,
     &                                 WFFXM,WFFYM,LAY,ATOM,RMATOTP,
     &                                 RMATSTP,ALPHA(LAY,ATOM),OMHAR,
     &                                 FINAL,IO,INITIAL,AMAT1T,AMAT2T,
     &                                 WFFM,PFI,NCPA,ATA,ICV)
                        END IF
                     END IF
                  END IF
               END IF
C JMJM MODIFIED BY JM: TO BE TESTED
C               IF ( IBLOCH.EQ.4 ) THEN
C                  IF ( NCPA.EQ.0 .OR. CPAATOM(LAY,ATOM).EQ.0 ) THEN
C                     CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,
C     &                           AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,WFFYM,
C     &                           LAY,ATOM,RMATOTP,RMATSTP,
C     &                           ALPHA(LAY,ATOM),OMHAR,FINAL,ITMP,
C     &                           INITIAL,AMAT1T,AMAT2T,WFFM,PFI,NCPA,
C     &                           ATA,ICV)
C                  ELSE IF ( NCPA.GT.0 ) THEN
C                     IF ( ATA.EQ.1 ) THEN
C                        CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,
C     &                              AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,
C     &                              WFFYM,LAY,ATOM,RMATOTP,RMATSTP,
C     &                              ALPHA(LAY,ATOM),OMHAR,FINAL,ITMP,
C     &                              INITIAL,AMAT1T,AMAT2T,WFFM,PFI,NCPA,
C     &                              ATA,ICV)
C                     ELSE IF ( ATA.EQ.0 ) THEN
C                        IF ( CPAATOM(LAY,ATOM).EQ.1 )
C     &                       CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
C     &                       AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,
C     &                       WFFYM,LAY,ATOM,RMATOTP,RMATSTP,
C     &                       ALPHA(LAY,ATOM),OMHAR,FINAL,ITMP,INITIAL,
C     &                       AMAT1T,AMAT2T,WFFMZ,PFI,NCPA,ATA,ICV)
C                     END IF
C                  END IF
C               END IF
C
C-------------------------------------------------------------------
               IF ( IBLOCH.EQ.2 )
     &              CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,
     &              AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,WFFYM,LAY,ATOM,
     &              RMATOTP,RMATSTP,ALPHA(LAY,ATOM),OMHAR,FINAL,ITMP,
     &              INITIAL,AMAT1T,AMAT2T,WFFMZ,PFI,NCPA,ATA,ICV)
C
               IF ( NCPA.GT.0 ) THEN
                  IF ( ATA.EQ.1 ) THEN
                     CALL ROMATMCPA(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,
     &                              AMAT1V,AMAT2V,RDIP1M,WFFM,WFFXM,
     &                              WFFYM,LAY,ATOM,RMATOC,RMATSC,
     &                              ALPHA(LAY,ATOM),OMHAR,FINAL,ITMP,
     &                              INITIAL,AMAT1T,AMAT2T,CM1,MAXG,WFFM,
     &                              G2,PMS,PMS1,PMS2,ICV)
                  ELSE IF ( ATA.EQ.0 ) THEN
                     IF ( CPAATOM(LAY,ATOM).EQ.1 ) THEN
                        CALL ROMATMCPA(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                                 AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,
     &                                 WFFXM,WFFYM,LAY,ATOM,RMATOC,
     &                                 RMATSC,ALPHA(LAY,ATOM),OMHAR,
     &                                 FINAL,ITMP,INITIAL,AMAT1T,AMAT2T,
     &                                 CM1,MAXG,WFFMZ,G2,PMS,PMS1,PMS2,
     &                                 ICV)
                     ELSE IF ( CPAATOM(LAY,ATOM).EQ.0 ) THEN
                        CALL ROMATMCPA(AMAT1,AMAT1X,AMAT2X,AMAT1Y,
     &                                 AMAT2Y,AMAT1V,AMAT2V,RDIP1M,WFFM,
     &                                 WFFXM,WFFYM,LAY,ATOM,RMATOC,
     &                                 RMATSC,ALPHA(LAY,ATOM),OMHAR,
     &                                 FINAL,IO,INITIAL,AMAT1T,AMAT2T,
     &                                 CM1,MAXG,WFFM,G2,PMS,PMS1,PMS2,
     &                                 ICV)
C
                     END IF
                  END IF
C
C
C=======================================================================
C                   sumup over concs and U_transformation
C                              regular MELE
C=======================================================================
C
                  IF ( THERMAL_VIBRA_FLUCT ) THEN
                     CM1(1:MQD,1:MQD) = RMATOC(1:MQD,1:MQD,LAY,ATOM,IO)
                     C4(1:MQD,1:MQD) = UMAT_VT_PES(LAY,ATOM,1:MQD,1:MQD,
     &                                 IO)
                     C5 = TRANSPOSE(C4)
                     C3 = MATMUL(CM1,C5)
                     RMATOC(1:MQD,1:MQD,LAY,ATOM,IO) = C3(1:MQD,1:MQD)
C
                     CM1(1:MQD,1:MQD) = UMAT_VT_PES(LAY,ATOM,1:MQD,1:MQD
     &                                  ,IO)
                     C4(1:MQD,1:MQD) = RMATSC(1:MQD,1:MQD,LAY,ATOM,IO)
                     C3 = MATMUL(CM1,C4)
                     RMATSC(1:MQD,1:MQD,LAY,ATOM,IO) = C3(1:MQD,1:MQD)
                  END IF
C
C=======================================================================
C                   PRINT OUT OF MELES
C=======================================================================
                  IF ( IP.GE.2 ) THEN
                     WRITE (NOUT1,99008) ATOM,LAY,IO
                     DO I1 = 1,MQD
                        DO I2 = 1,MQD
                           IF ( CDABS(RMATOC(I1,I2,LAY,ATOM,IO))
     &                          .GT.EPS12 ) WRITE (NOUT1,99007) I1,I2,
     &                          RMATOC(I1,I2,LAY,ATOM,IO)
                        END DO
                     END DO
                     WRITE (NOUT1,99009) ATOM,LAY,IO
                     DO I1 = 1,MQD
                        DO I2 = 1,MQD
                           IF ( CDABS(RMATSC(I1,I2,LAY,ATOM,IO))
     &                          .GT.EPS12 ) WRITE (NOUT1,99007) I1,I2,
     &                          RMATSC(I1,I2,LAY,ATOM,IO)
                        END DO
                     END DO
                     WRITE (NOUT1,99010) ATOM,LAY,IO
                     DO I1 = 1,MQD
                        DO I2 = 1,MQD
                           IF ( CDABS(RMATOTP(I1,I2,LAY,ATOM,IO))
     &                          .GT.EPS12 ) WRITE (NOUT1,99007) I1,I2,
     &                          RMATOTP(I1,I2,LAY,ATOM,IO)
                        END DO
                     END DO
                     WRITE (NOUT1,99011) ATOM,LAY,IO
                     DO I1 = 1,MQD
                        DO I2 = 1,MQD
                           IF ( CDABS(RMATSTP(I1,I2,LAY,ATOM,IO))
     &                          .GT.EPS12 ) WRITE (NOUT1,99007) I1,I2,
     &                          RMATSTP(I1,I2,LAY,ATOM,IO)
                        END DO
                     END DO
                  END IF
C=======================================================================
C
                  IF ( ATA.EQ.1 ) THEN
                     RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                  = RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                  + CONZ*RMATOC(1:MQD,1:MQD,LAY,ATOM,IO)
                     RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                  = RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                  + CONZ*RMATSC(1:MQD,1:MQD,LAY,ATOM,IO)
                  ELSE IF ( ATA.EQ.0 ) THEN
                     IF ( CPAATOM(LAY,ATOM).EQ.1 ) THEN
                        C3(1:MQD,1:MQD)
     &                     = CONZ*RMATOC(1:MQD,1:MQD,LAY,ATOM,IO)*CI/PFI
                        CM1 = MATMUL(C2,C3)
                        RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                     = RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                     + CM1(1:MQD,1:MQD)
C
                        C3(1:MQD,1:MQD)
     &                     = CONZ*RMATSC(1:MQD,1:MQD,LAY,ATOM,IO)*CI/PFI
                        CM1 = MATMUL(C3,C2)
                        RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                     = RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                     + CM1(1:MQD,1:MQD)
C
                     ELSE IF ( CPAATOM(LAY,ATOM).EQ.0 ) THEN
                        RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                     = RMATO(1:MQD,1:MQD,LAY,ATOM)
     &                     + CONZ*RMATOC(1:MQD,1:MQD,LAY,ATOM,IO)
                        RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                     = RMATS(1:MQD,1:MQD,LAY,ATOM)
     &                     + CONZ*RMATSC(1:MQD,1:MQD,LAY,ATOM,IO)
                     END IF
                  END IF
C
               ELSE IF ( NCPA.EQ.0 ) THEN
C
                  RMATO(1:MQD,1:MQD,LAY,ATOM)
     &               = RMATOTP(1:MQD,1:MQD,LAY,ATOM,1)
                  RMATS(1:MQD,1:MQD,LAY,ATOM)
     &               = RMATSTP(1:MQD,1:MQD,LAY,ATOM,1)
C
               END IF
C
C
C=======================================================================
C                   calculate irregular MELE
C=======================================================================
C
               CALL DOUBLERADM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,
     &                         AMAT2V,RDIP1M,WFFM,WFFXM,WFFYM,LAY,ATOM,
     &                         ZMAT0TP,ALPHA(LAY,ATOM),OMHAR,FINAL,
     &                         INITIAL,IRSTATE,RADMESHM,AMAT1T,AMAT2T,
     &                         ITMP,ICV)
C
C=======================================================================
C                   sumup over concs and U_transformation
C                              irregular MELE
C=======================================================================
C
               IF ( NCPA.GT.0 ) THEN
                  IF ( THERMAL_VIBRA_FLUCT ) THEN
C
                     CM1(1:MQD,1:MQD) = UMAT_VT_PES(LAY,ATOM,1:MQD,1:MQD
     &                                  ,IO)
                     C4(1:MQD,1:MQD) = ZMAT0TP(1:MQD,1:MQD,LAY,ATOM,IO)
                     C3 = MATMUL(CM1,C4)
C
                     CM1(1:MQD,1:MQD) = TRANSPOSE(UMAT_VT_PES(LAY,ATOM,1
     &                                  :MQD,1:MQD,IO))
                     C4 = MATMUL(C3,CM1)
                     ZMAT0(1:MQD,1:MQD,LAY,ATOM)
     &                  = ZMAT0(1:MQD,1:MQD,LAY,ATOM)
     &                  + CONZ*C4(1:MQD,1:MQD)
C
                  ELSE
                     ZMAT0(1:MQD,1:MQD,LAY,ATOM)
     &                  = ZMAT0(1:MQD,1:MQD,LAY,ATOM)
     &                  + CONZ*ZMAT0TP(1:MQD,1:MQD,LAY,ATOM,IO)
                  END IF
C
               ELSE IF ( NCPA.EQ.0 ) THEN
                  ZMAT0(1:MQD,1:MQD,LAY,ATOM)
     &               = ZMAT0TP(1:MQD,1:MQD,LAY,ATOM,1)
               END IF
C
            END DO LOOP_IO
C**************************************************************  LOOP_IO
C
C=======================================================================
C                   PRINT OUT OF MELES
C=======================================================================
            IF ( IP.GE.2 ) THEN
               WRITE (NOUT1,99005) ATOM,LAY
               DO I = 1,MQD
                  DO J = 1,MQD
                     IF ( CDABS(RMATO(I,J,LAY,ATOM)).GT.EPS12 )
     &                    WRITE (NOUT1,99001) I,J,RMATO(I,J,LAY,ATOM)
                  END DO
               END DO
               WRITE (NOUT1,99006) ATOM,LAY
               DO I = 1,MQD
                  DO J = 1,MQD
                     IF ( CDABS(RMATS(I,J,LAY,ATOM)).GT.EPS12 )
     &                    WRITE (NOUT1,99001) I,J,RMATS(I,J,LAY,ATOM)
                  END DO
               END DO
            END IF
C
            IF ( IP.GE.2 ) THEN
               WRITE (NOUT1,99003) ATOM,LAY
               DO I = 1,MQD
                  DO J = 1,MQD
                     IF ( CDABS(ZMAT0(I,J,LAY,ATOM)).GT.EPS12 )
     &                    WRITE (NOUT1,99001) I,J,ZMAT0(I,J,LAY,ATOM)
                  END DO
               END DO
            END IF
C=======================================================================
C
C
C
C=======================================================================
C                   calculate incoherent MELEs
C=======================================================================
C
            IF ( IBLOCH.EQ.4 ) THEN
C
               IF ( CPAATOM(LAY,ATOM).EQ.1 )
     &              CALL ZMATINCO(RMATOTP,RMATSTP,RMATO,RMATS,ZMAT0INC,
     &              TAUMT,LAY,ATOM,UMAT_VT_PES)
            ELSE IF ( IBLOCH.EQ.2 ) THEN
               CALL ZMATINCO(RMATOTP,RMATSTP,RMATO,RMATS,ZMAT0INC,TAUMT,
     &                       LAY,ATOM,UMAT_VT_PES)
            END IF
            IF ( IP.GE.2 ) THEN
               WRITE (6,*) 'CPAATM',LAY,ATOM
               IQ = IQ_SPR_LAYAT(ATOM,LAY)
               DO IO = 1,NVFO_Q(IQ) + 1
                  WRITE (NOUT1,99004) ATOM,LAY,IO
                  DO I = 1,MQD
                     DO J = 1,MQD
                        IF ( CDABS(ZMAT0INC(I,J,LAY,ATOM,IO)).GT.EPS12 )
     &                       WRITE (NOUT1,99002) I,J,IO,
     &                              ZMAT0INC(I,J,LAY,ATOM,IO)
                     END DO
                  END DO
               END DO
            END IF
C
C=======================================================================
C
         END DO LOOP_ATOM
C**************************************************************LOOP_ATOM
      END DO LOOP_LAY
C**************************************************************LOOP_LAY
C
C
C
99001 FORMAT (2I3,2(1x,e15.7))
99002 FORMAT (3I3,2(1x,e15.7))
99003 FORMAT ('coherent double matrix elements for atom',2x,i3,2x,
     &        'layer',2x,i3)
99004 FORMAT ('incoherent zmatrix for atom',2x,i3,2x,'layer',2x,i3,2x,
     &        'type',2x,i3)
99005 FORMAT ('coherent matrix elements rmato for atom',2x,i3,2x,
     &        'layer',2x,i3)
99006 FORMAT ('coherent matrix elements rmats for atom',2x,i3,2x,
     &        'layer',2x,i3)
99007 FORMAT (2I3,2(1x,e15.7))
99008 FORMAT ('radial matrix elements urmatoc for atom',2x,i3,2x,
     &        'layer',2x,i3,2x,'type',2x,i3)
99009 FORMAT ('radial matrix elements urmatsc for atom',2x,i3,2x,
     &        'layer',2x,i3,2x,'type',2x,i3)
99010 FORMAT ('radial matrix elements urmatotp for atom',2x,i3,2x,
     &        'layer',2x,i3,2x,'type',2x,i3)
99011 FORMAT ('radial matrix elements urmatstp for atom',2x,i3,2x,
     &        'layer',2x,i3,2x,'type',2x,i3)
      END
