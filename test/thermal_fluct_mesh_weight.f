C*==thermal_fluct_mesh_weight.f    processed by SPAG 6.70Rc at 11:32 on  7 Mar 2017
      SUBROUTINE THERMAL_FLUCT_MESH_WEIGHT(I_TEMP_LAT,TEMP_LAT,MM_MC,
     &   NT_MAGNETIC,IT_MAG,W_WEISS_T)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:KB_SI,EV_J,RY_EV,PI,DEG_ARC
      USE MOD_THERMAL,ONLY:FTET_FT,FPHI_FT,TET_FLUCT,PHI_FLUCT,
     &    NPHI_FLUCT,NTET_FLUCT,NVIBRA,X_VFT,X_FT,X_VT,NFLUCT,W0_FLUCT,
     &    W0_TET,NTET_FLUCT_POT
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC,X_CHEM_T
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SCF,ONLY:SCF_THETA_DEPENDENT_POT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_FLUCT_MESH_WEIGHT')
      LOGICAL USE_WEISS_AVRG
      PARAMETER (USE_WEISS_AVRG=.TRUE.)
C
C Dummy arguments
C
      INTEGER I_TEMP_LAT,NT_MAGNETIC
      REAL*8 TEMP_LAT
      INTEGER IT_MAG(NT)
      REAL*8 MM_MC(NT_MAGNETIC),W_WEISS_T(NTMAX)
C
C Local variables
C
      REAL*8 BETA_WEISS,DPHI,DTET,SUMW,SUMW0,TET_ARC,W0,W0_PHI,
     &       WEISS_AVRG,X_DUM
      INTEGER IFLUCT,IFT,IPHI,IT,ITET,ITT,IT_CHEM,IT_MC,IVFT,IVIBRA,IVT,
     &        N,NT_CHEM
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( SCF_THETA_DEPENDENT_POT ) THEN
C
C=======================================================================
C           potential depends on  IT_CHEM and THETA
C
C     all types IT are pseudo-types representing the chemical types
C     together with  NTET_FLUCT_POT  spin angles theta
C     IFLUCT runs only over the angle PHI, i.e. IPHI==IFLUCT
C=======================================================================
C
         IF ( NFLUCT.EQ.2 ) THEN
C
            CALL STOP_REGULAR(ROUTINE,'DLM with 2 spin projections')
C
         ELSE IF ( NTET_FLUCT.EQ.1 ) THEN
C
C            CALL STOP_REGULAR(ROUTINE,'Conical distribution of M')
            WRITE (6,'(A)') 'THETA-dependent potential'
C
         END IF
C
         NT_CHEM = NT/NTET_FLUCT_POT
C
         IFT = 0
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         DO IT_CHEM = 1,NT_CHEM
C
C ==================================== different models for fluctuations
C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
C
            DTET = PI/(NTET_FLUCT_POT-1)
            DPHI = 2.D0*PI/(NPHI_FLUCT)
            W0_PHI = DPHI
            SUMW = 0.D0
C
C --------------------------------------------- M directions over sphere
C ----------------------------------------------------------------------
            IF ( SUM(DABS(W_WEISS_T(1:NT))).LE.1.D-8 ) THEN
C
               DO ITET = 1,NTET_FLUCT_POT
                  SUMW0 = 0.D0
                  DO IPHI = 1,NFLUCT
C
                     IFT = IFT + 1
C                     FTET_FT(IFT) = 0.0D0
C SM:  FTET_FT(IFT) is used in CHRDNS  !
                     FTET_FT(IFT) = (ITET-1)*DTET/DEG_ARC
                     FPHI_FT(IFT) = (IPHI-1)*DPHI/DEG_ARC
C
                     TET_ARC = (ITET-1)*DTET
                     IF ( DABS(SIN(TET_ARC)).LE.1.D-8 ) THEN
                        W0_TET(ITET) = SIN(DTET/2.0D0)*0.5D0
                     ELSE
                        W0_TET(ITET) = (SIN(TET_ARC+DTET/2.0D0)+SIN(
     &                                 TET_ARC-DTET/2.0D0))*0.5D0
                     END IF
                     W0_TET(ITET) = W0_TET(ITET)*DTET
                     W0_FLUCT(IPHI) = W0_PHI/(2.D0*PI)
C
                     X_FT(IFT) = W0_TET(ITET)*W0_PHI
C
                     SUMW = SUMW + X_FT(IFT)
                     SUMW0 = SUMW0 + W0_FLUCT(IPHI)
C
                  END DO
C
                  W0_FLUCT(1:NFLUCT) = W0_FLUCT(1:NFLUCT)/SUMW0
C
                  IF ( ABS(1D0-SUM(W0_FLUCT(1:NFLUCT))).GT.1D-10 )
     &                 CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
               END DO
C
               N = NFLUCT*NTET_FLUCT_POT
C
               X_FT(IFT-N+1:IFT) = X_FT(IFT-N+1:IFT)/SUMW
C
               IF ( ABS(1D0-SUM(X_FT(IFT-N+1:IFT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
            ELSE
C
C--------------------------------------------------- AVERAGE WEISS FIELD
C
               IF ( USE_WEISS_AVRG ) THEN
C
                  WEISS_AVRG = 0.D0
                  SUMW0 = 0.D0
                  DO ITET = 1,NTET_FLUCT_POT
C
                     ITT = (IT_CHEM-1)*NT_CHEM + ITET
C
                     TET_ARC = (ITET-1)*DTET
                     IF ( DABS(SIN(TET_ARC)).LE.1.D-8 ) THEN
                        W0_TET(ITET) = SIN(DTET/2.0D0)*0.5D0
                     ELSE
                        W0_TET(ITET) = (SIN(TET_ARC+DTET/2.0D0)+SIN(
     &                                 TET_ARC-DTET/2.0D0))*0.5D0
                     END IF
                     W0_TET(ITET) = W0_TET(ITET)*DTET
C
                     SUMW0 = SUMW0 + W0_TET(ITET)
C
                     WEISS_AVRG = WEISS_AVRG + W_WEISS_T(ITT)
     &                            *W0_TET(ITET)
                  END DO
C
                  WEISS_AVRG = WEISS_AVRG/SUMW0
C
               END IF
C-----------------------------------------------------------------------
C --------------------------------------------- M directions over sphere
C
C ---------------------------------------------------- weighting factors
C
               DO ITET = 1,NTET_FLUCT_POT
C
                  IT = (IT_CHEM-1)*NT_CHEM + ITET
C
C
                  IF ( USE_WEISS_AVRG ) THEN
C------------------------------- for the case of integrated Weiss field
                     BETA_WEISS = (EV_J*RY_EV/KB_SI)*DABS(WEISS_AVRG)
     &                            /TEMP_LAT
                  ELSE
C--------------------------- for the case of type dependent Weiss field
                     BETA_WEISS = (EV_J*RY_EV/KB_SI)*W_WEISS_T(IT)
     &                            /TEMP_LAT
C
                  END IF
C
                  SUMW0 = 0.D0
                  DO IPHI = 1,NFLUCT
C
                     IFT = IFT + 1
C                     FTET_FT(IFT) = 0.0D0
C SM:  FTET_FT(IFT) is used in CHRDNS  !
                     FTET_FT(IFT) = (ITET-1)*DTET/DEG_ARC
                     FPHI_FT(IFT) = (IPHI-1)*DPHI/DEG_ARC
C
                     TET_ARC = (ITET-1)*DTET
                     IF ( DABS(SIN(TET_ARC)).LE.1.D-8 ) THEN
                        W0_TET(ITET) = SIN(DTET/2.0D0)*0.5D0
                     ELSE
                        W0_TET(ITET) = (SIN(TET_ARC+DTET/2.0D0)+SIN(
     &                                 TET_ARC-DTET/2.0D0))*0.5D0
                     END IF
                     W0_TET(ITET) = W0_TET(ITET)*DTET
C
                     W0_FLUCT(IPHI) = W0_PHI/(2.D0*PI)
C
                     IF ( USE_WEISS_AVRG ) THEN
C------------------------------- for the case of integrated Weiss field
                        X_DUM = BETA_WEISS*(COS(TET_ARC)-1.D0)
                     ELSE
C--------------------------- for the case of type dependent Weiss field
                        X_DUM = BETA_WEISS
                     END IF
C
                     X_FT(IFT) = W0_TET(ITET)*W0_PHI*EXP(X_DUM)
C
                     SUMW0 = SUMW0 + W0_FLUCT(IPHI)
                     SUMW = SUMW + X_FT(IFT)
C
                  END DO
C
                  W0_FLUCT(1:NFLUCT) = W0_FLUCT(1:NFLUCT)/SUMW0
C
                  IF ( ABS(1D0-SUM(W0_FLUCT(1:NFLUCT))).GT.1D-10 )
     &                 CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
               END DO
C
               N = NFLUCT*NTET_FLUCT_POT
C
               X_FT(IFT-N+1:IFT) = X_FT(IFT-N+1:IFT)/SUMW
C
               IF ( ABS(1D0-SUM(X_FT(IFT-N+1:IFT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
C ------------------------------------------- TEST OUTPUT: Normalization
C
               IF ( IPRINT.GT.1 ) THEN
                  WRITE (1112,*) I_TEMP_LAT,SUMW
                  SUMW = 0.D0
                  DO IFLUCT = 1,NFLUCT
                     TET_ARC = TET_FLUCT(IFLUCT)
                     SUMW = SUMW + X_FT(IFLUCT)
                     WRITE (9500+I_TEMP_LAT,'(2F30.10)') TET_ARC,
     &                      X_FT(IFLUCT)
                     WRITE (9300+I_TEMP_LAT,*) IFLUCT,X_FT(IFLUCT)
                  END DO
                  WRITE (1113,*) IT,I_TEMP_LAT,SUMW
               END IF
C
            END IF
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C           complete the weights of vibrations and fluctuations
C ----------------------------------------------------------------------
C
            DO ITET = 1,NTET_FLUCT_POT
               IT = (IT_CHEM-1)*NT_CHEM + ITET
               IFT = NFLUCT*(IT-1)
C
               CONC(IT) = 0
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
                  IVT = (IT-1)*NVIBRA
                  DO IVIBRA = 1,NVIBRA
                     IVT = IVT + 1
C
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
C
C     SM:  X_CHEM_T(IT) is not properly specified and has to be updated:
C          X_CHEM_T(IT) is the same as initial CONC(IT)
C           therefore  I use here  X_CHEM_T(IT)*NTET_FLUCT_POT
C
C                     X_VFT(IVFT) = X_VT(IVT)*X_FT(IFT)*X_CHEM_T(IT)
                     X_VFT(IVFT) = X_VT(IVT)*X_FT(IFT)
     &                             *(X_CHEM_T(IT)*NTET_FLUCT_POT)
                     CONC(IT) = CONC(IT) + X_FT(IFT)
                     IF ( IPRINT.GT.1 ) WRITE (9600,'(I5,5F15.8)') IVFT,
     &                    IFT,IT,X_VFT(IVFT),X_VT(IVT),X_FT(IFT),
     &                    X_CHEM_T(IT)*NTET_FLUCT_POT,CONC(IT)
C
                  END DO
               END DO
            END DO
C ----------------------------------------------------------------------
C
         END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      ELSE
C
C=======================================================================
C                potential depends only on  IT == IT_CHEM
C
C     all types IT represent only chemical types
C     IFLUCT runs over the angles THETA and PHI
C=======================================================================
C
C     DO IT = 1, NT_MAGNETIC             ! to be used for big cell
         IT_MC = 0
         IFT = 0
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         DO IT = 1,NT
C
            IT_MC = IT_MC + IT_MAG(IT)
C
C ==================================== different models for fluctuations
C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
            IF ( NFLUCT.EQ.2 ) THEN
C ------------------------------------------ DLM with 2 spin projections
C
               X_DUM = (1-MM_MC(IT_MC))/2.D0
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
                  FTET_FT(IFT) = TET_FLUCT(IFLUCT)
                  FPHI_FT(IFT) = PHI_FLUCT(IFLUCT)
                  X_FT(IFT) = DABS((IFLUCT-1)-X_DUM)
                  WRITE (5599,'(2I5, 3F10.5)') IFLUCT,IFT,X_FT(IFT),
     &                   FTET_FT(IFT),FPHI_FT(IFT)
               END DO
C
            ELSE IF ( NTET_FLUCT.EQ.1 ) THEN
C --------------------------------- Conical distribution of M directions
C
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
                  TET_FLUCT(IFLUCT) = ACOS(MM_MC(IT_MC))
                  FTET_FT(IFT) = TET_FLUCT(IFLUCT)/DEG_ARC
                  FPHI_FT(IFT) = PHI_FLUCT(IFLUCT)
                  X_FT(IFT) = 1.0/NPHI_FLUCT
               END DO
C
C --------------------------------------------- M directions over sphere
C ----------------------------------------------------------------------
            ELSE IF ( DABS(W_WEISS_T(IT)).LE.1.D-10 ) THEN
C
               DTET = PI/(NTET_FLUCT-1)
               DPHI = 2.D0*PI/(NPHI_FLUCT)
               SUMW0 = 0.D0
C
               DO IFLUCT = 1,NFLUCT
C
                  IFT = IFT + 1
                  FTET_FT(IFT) = TET_FLUCT(IFLUCT)
                  FPHI_FT(IFT) = PHI_FLUCT(IFLUCT)
                  TET_ARC = TET_FLUCT(IFLUCT)*DEG_ARC
                  W0 = SIN(TET_ARC)
                  W0 = W0*DTET*DPHI
C
                  W0 = W0/(4.D0*PI)
C
                  IF ( (ABS(TET_FLUCT(IFLUCT)).LT.1D-10) .OR. 
     &                 (ABS(TET_FLUCT(IFLUCT)-PI).LT.1D-10) )
     &                 W0 = W0/2.D0
C
                  SUMW0 = SUMW0 + W0
                  IF ( IPRINT.GT.1 ) WRITE (9500+I_TEMP_LAT,'(2F30.10)')
     &                 TET_ARC,W0
C
                  W0_FLUCT(IFLUCT) = W0
                  X_FT(IFT) = W0
C
               END DO
C
               IF ( IPRINT.GT.1 ) WRITE (1113,*) I_TEMP_LAT,SUMW0
C
               W0_FLUCT(1:NFLUCT) = W0_FLUCT(1:NFLUCT)/SUMW0
               X_FT(IFT-NFLUCT+1:IFT) = X_FT(IFT-NFLUCT+1:IFT)/SUMW0
C
               IF ( ABS(1D0-SUM(W0_FLUCT(1:NFLUCT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
               IF ( ABS(1D0-SUM(X_FT(IFT-NFLUCT+1:IFT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
            ELSE
C --------------------------------------------- M directions over sphere
C --------------------------------------------------- Partition function
C
               BETA_WEISS = (EV_J*RY_EV/KB_SI)*DABS(W_WEISS_T(IT))
     &                      /TEMP_LAT
C
C ---------------------------------------------------- weighting factors
C
               DTET = PI/(NTET_FLUCT-1)
               DPHI = 2.D0*PI/(NPHI_FLUCT)
               SUMW = 0.D0
               SUMW0 = 0.D0
C
               DO IFLUCT = 1,NFLUCT
C
                  IFT = IFT + 1
                  FTET_FT(IFT) = TET_FLUCT(IFLUCT)
                  FPHI_FT(IFT) = PHI_FLUCT(IFLUCT)
                  TET_ARC = TET_FLUCT(IFLUCT)*DEG_ARC
                  IF ( DABS(SIN(TET_ARC)).LE.1.D-8 ) THEN
                     W0 = SIN(DTET/2.0D0)*0.5D0
                  ELSE
                     W0 = (SIN(TET_ARC+DTET/2.0D0)
     &                    +SIN(TET_ARC-DTET/2.0D0))*0.5D0
                  END IF
                  W0 = W0*DTET*DPHI
C
                  X_DUM = BETA_WEISS*(COS(TET_ARC)-1.D0)
C
                  W0_FLUCT(IFLUCT) = W0/(4.D0*PI)
                  X_FT(IFT) = W0*EXP(X_DUM)
C
                  SUMW = SUMW + X_FT(IFT)
                  SUMW0 = SUMW0 + W0_FLUCT(IFLUCT)
C
               END DO
C
               IF ( IPRINT.GT.1 ) WRITE (1112,*) I_TEMP_LAT,SUMW
C
               X_FT(IFT-NFLUCT+1:IFT) = X_FT(IFT-NFLUCT+1:IFT)/SUMW
               W0_FLUCT(1:NFLUCT) = W0_FLUCT(1:NFLUCT)/SUMW0
C
               IF ( ABS(1D0-SUM(W0_FLUCT(1:NFLUCT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
               IF ( ABS(1D0-SUM(X_FT(IFT-NFLUCT+1:IFT))).GT.1D-10 )
     &              CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
C ------------------------------------------- TEST OUTPUT: Normalization
               IF ( IPRINT.GT.1 ) THEN
                  SUMW = 0.D0
                  DO IFLUCT = 1,NFLUCT
                     TET_ARC = TET_FLUCT(IFLUCT)
                     SUMW = SUMW + X_FT(IFLUCT)
                     WRITE (9500+I_TEMP_LAT,'(2F30.10)') TET_ARC,
     &                      X_FT(IFLUCT)
                     WRITE (9300+I_TEMP_LAT,*) IFLUCT,X_FT(IFLUCT)
                  END DO
                  WRITE (1113,*) IT,I_TEMP_LAT,SUMW
               END IF
C
            END IF
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C           complete the weights of vibrations and fluctuations
C ----------------------------------------------------------------------
C
            IFT = NFLUCT*(IT-1)
            DO IFLUCT = 1,NFLUCT
               IFT = IFT + 1
C
               IVT = (IT-1)*NVIBRA
               DO IVIBRA = 1,NVIBRA
                  IVT = IVT + 1
C
                  IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                   *NFLUCT + IFLUCT
C
                  X_VFT(IVFT) = X_VT(IVT)*X_FT(IFT)*CONC(IT)
                  IF ( IPRINT.GT.1 ) WRITE (9600,'(I5,3F15.8)') IVFT,
     &                 X_VFT(IVFT),X_VT(IVT),X_FT(IFT)
C
               END DO
            END DO
C ----------------------------------------------------------------------
C
         END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END IF
      END
