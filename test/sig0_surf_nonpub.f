C*==sig0_surf_nonpub.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG0_SURF_NONPUB(MIRR_VFTAB,MIRR_VFTBA,MREG_VFTAB,
     &                            MREG_VFTBA,REBS_VFTB,IMTAU_VFTA,
     &                            IMTSS_VFTA,RETAU_VFTB,RETSS_VFTB,
     &                            SIG0IR_BS,SIG0IR_SS,SIG0IR_BS_VFT,
     &                            SIG0IR_SS_VFT,SIG0IRSS_SS,SIG0IR_IRR,
     &                            SIG0IRSS_SS_VFT,SIG0IR_VFT,
     &                            SIG0IR_IRR_VFT,SIG0IR,N,M,MUE,NUE,
     &                            IVFT,IXX,NSPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *  non-public part of  SIG0_SURF                                   *
C   *  providing off-diagonal elements of the response tensors         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,WKM2
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_THERMAL,ONLY:NVFTMAX,X_VFT
      IMPLICIT NONE
C*--SIG0_SURF_NONPUB21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IVFT,IXX,M,MUE,N,NSPINPROJ,NUE
      COMPLEX*16 IMTAU_VFTA(NKMMAX,NKMMAX),IMTSS_VFTA(NKMMAX,NKMMAX),
     &           MIRR_VFTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           MIRR_VFTBA(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           MREG_VFTAB(NKMMAX,NKMMAX,3,NSPINPROJ),
     &           MREG_VFTBA(NKMMAX,NKMMAX,3,NSPINPROJ),
     &           REBS_VFTB(NKMMAX,NKMMAX),RETAU_VFTB(NKMMAX,NKMMAX),
     &           RETSS_VFTB(NKMMAX,NKMMAX),SIG0IR(3,3),SIG0IRSS_SS(3,3),
     &           SIG0IRSS_SS_VFT(3,3,NVFTMAX,2),SIG0IR_BS(3,3),
     &           SIG0IR_BS_VFT(3,3,NVFTMAX),SIG0IR_IRR(3,3),
     &           SIG0IR_IRR_VFT(3,3,NVFTMAX),SIG0IR_SS(3,3),
     &           SIG0IR_SS_VFT(3,3,NVFTMAX,2),SIG0IR_VFT(3,3,NVFTMAX)
C
C Local variables
C
      INTEGER I
      COMPLEX*16 TRACE
C
C*** End of declarations rewritten by SPAG
C
C-------------------------------- i ( J_m Im G j_n - j_n Im G J_m ) Re G
C-----------------------------------------------------------regular part
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),RETAU_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),WKM2,WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),RETAU_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),WKM2,WKM1)
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IR(MUE,NUE) = SIG0IR(MUE,NUE) + X_VFT(IVFT)*CI*TRACE
C
      SIG0IR_VFT(MUE,NUE,IVFT) = SIG0IR_VFT(MUE,NUE,IVFT) + X_VFT(IVFT)
     &                           *CI*TRACE
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C---------------------------------------------------------irregular part
C
      CALL CMATMUL(N,M,IMTAU_VFTA,MIRR_VFTBA(1,1,NUE,MUE,IXX),WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,IMTAU_VFTA,MIRR_VFTAB(1,1,MUE,NUE,IXX),WKM1)
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IR_IRR(MUE,NUE) = SIG0IR_IRR(MUE,NUE) + X_VFT(IVFT)*TRACE*CI
C
      SIG0IR_IRR_VFT(MUE,NUE,IVFT) = SIG0IR_IRR_VFT(MUE,NUE,IVFT)
     &                               + X_VFT(IVFT)*TRACE*CI
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-------------------------------- here some tools for analysis follow
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C-----------------------------------------------------------------------
C-------------------------------------------- single site contribution 1
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-----------------------------------------------------------------------
C----------------------------------------------------Z tau Z * ZJ   term
C
      SIG0IR_SS(MUE,NUE) = SIG0IR_SS(MUE,NUE) - X_VFT(IVFT)*TRACE*CI
      SIG0IR_SS_VFT(MUE,NUE,IVFT,1) = SIG0IR_SS_VFT(MUE,NUE,IVFT,1)
     &                                - X_VFT(IVFT)*TRACE*CI
C----------------------------------------------------Z tau Z * ZtZ  term
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),RETSS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),WKM2,WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),RETSS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),WKM2,WKM1)
C
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IR_SS(MUE,NUE) = SIG0IR_SS(MUE,NUE) + X_VFT(IVFT)*CI*TRACE
      SIG0IR_SS_VFT(MUE,NUE,IVFT,2) = SIG0IR_SS_VFT(MUE,NUE,IVFT,2)
     &                                + X_VFT(IVFT)*TRACE*CI
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C------------------------------------------ back scattering contribution
C-------------------------------------------- Z tau Z * Z(tau-t)Z   term
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),REBS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),WKM2,WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),REBS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTAU_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),WKM2,WKM1)
C
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IR_BS(MUE,NUE) = SIG0IR_BS(MUE,NUE) + X_VFT(IVFT)*CI*TRACE
      SIG0IR_BS_VFT(MUE,NUE,IVFT) = SIG0IR_BS_VFT(MUE,NUE,IVFT)
     &                              + X_VFT(IVFT)*CI*TRACE
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C-------------------------------------------- single site contribution 2
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C----------------------------------------------------Z t   Z * ZtZ  term
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,NUE,1),RETSS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTSS_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),WKM2,WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,MREG_VFTAB(1,1,MUE,IXX),RETSS_VFTB,WKM1)
C
      CALL CMATMUL(N,M,IMTSS_VFTA,WKM1,WKM2)
C
      CALL CMATMUL(N,M,MREG_VFTBA(1,1,NUE,1),WKM2,WKM1)
C
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IRSS_SS(MUE,NUE) = SIG0IRSS_SS(MUE,NUE) + X_VFT(IVFT)*CI*TRACE
      SIG0IRSS_SS_VFT(MUE,NUE,IVFT,2) = SIG0IRSS_SS_VFT(MUE,NUE,IVFT,2)
     &                                  + X_VFT(IVFT)*TRACE*CI
C----------------------------------------------------Z t   Z * ZJ   term
C
      CALL CMATMUL(N,M,IMTSS_VFTA,MIRR_VFTBA(1,1,NUE,MUE,IXX),WKM1)
C
      TRACE = C0
      DO I = 1,NKM
         TRACE = TRACE + WKM1(I,I)
      END DO
C
      CALL CMATMUL(N,M,IMTSS_VFTA,MIRR_VFTAB(1,1,MUE,NUE,IXX),WKM1)
C
      DO I = 1,NKM
         TRACE = TRACE - WKM1(I,I)
      END DO
C
      SIG0IRSS_SS(MUE,NUE) = SIG0IRSS_SS(MUE,NUE) - X_VFT(IVFT)*TRACE*CI
      SIG0IRSS_SS_VFT(MUE,NUE,IVFT,1) = SIG0IRSS_SS_VFT(MUE,NUE,IVFT,1)
     &                                  - X_VFT(IVFT)*TRACE*CI
C
      END
