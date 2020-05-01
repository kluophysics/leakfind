C*==mod_thermal.f    processed by SPAG 6.70Rc at 08:24 on 30 Jan 2013
      MODULE MOD_THERMAL
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables used for   THERMAL  calculations  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NVIBRA_TAB,NFLUCT_TAB,NTET_FLUCT,NPHI_FLUCT,NTET_FLUCT_POT
      PARAMETER (NVIBRA_TAB=6,NFLUCT_TAB=6)
      REAL*8 R0,R1,RM1
      PARAMETER (R0=0D0,R1=1D0,RM1=-1D0)
C
C Local variables
C
      COMPLEX*16 FMAT_FT(:,:,:),UMAT_VT(:,:,:)
      REAL*8 FPHI_FT(:),FTET_FT(:),SVEC_VIBRA0(3,NVIBRA_TAB),
     &       SVEC_VT(:,:),X_VFT(:),
     &       MFMAT_FT(:,:,:),X_FT(:),X_VT(:),TEMP_LAT,
     &       TEMP_LAT_MIN,TEMP_LAT_MAX,MNT_TEMP_MCS,
     &       NHAT_MAG(3),EHAT_FLUCT(:,:,:),
     &       TET_FLUCT(:),PHI_FLUCT(:),W0_FLUCT(:),W0_TET(:)
      INTEGER IQ_AVFT(:,:),IVFT_VFOQ(:,:),NAVFT(:),NFLUCT,NFT,NFTMAX,
     &        NVFO_Q(:),NVFT,NVFTMAX,NVIBFLU,NVIBRA,NVT,NVTMAX,
     &        I_TEMP_LAT,N_TEMP_LAT,
     &        NFLUCT_T(:),NVIBFLU_T(:)
      LOGICAL THERMAL_CHECK_FLUCT
      CHARACTER*10 FLUCT_DIR_SETTING
      SAVE FMAT_FT,MFMAT_FT,
     &     IQ_AVFT,IVFT_VFOQ,NAVFT,NVFO_Q,SVEC_VT,UMAT_VT,
     &     X_VFT,X_FT,X_VT,
     &     THERMAL_CHECK_FLUCT,FTET_FT,FPHI_FT, 
     &     TEMP_LAT_MIN, TEMP_LAT_MAX,TEMP_LAT,I_TEMP_LAT,N_TEMP_LAT, 
     &     FLUCT_DIR_SETTING,MNT_TEMP_MCS,NFLUCT_T,NVIBFLU_T,NTET_FLUCT,
     &     NPHI_FLUCT,EHAT_FLUCT,
     &     TET_FLUCT,PHI_FLUCT,W0_FLUCT,W0_TET,NTET_FLUCT_POT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE X_VFT,X_FT,X_VT
      ALLOCATABLE SVEC_VT,FTET_FT,FPHI_FT,NVFO_Q,IVFT_VFOQ,IQ_AVFT
      ALLOCATABLE NAVFT,FMAT_FT,MFMAT_FT,UMAT_VT
      ALLOCATABLE NFLUCT_T,NVIBFLU_T
      ALLOCATABLE EHAT_FLUCT
      ALLOCATABLE PHI_FLUCT,TET_FLUCT,W0_FLUCT,W0_TET
C     
      DATA NVFT/0/,NVFTMAX/0/,NVT/0/,NVTMAX/0/,NFT/0/,NFTMAX/0/,
     &     NVIBRA/1/,NFLUCT/1/,NVIBFLU/1/
      DATA SVEC_VIBRA0/R1,R0,R0,R0,R1,R0,R0,R0,R1,RM1,R0,R0,R0,RM1,R0,
     &     R0,R0,RM1/!,R1,R1,R1,R1,R1,RM1,R1,RM1,R1,R1,RM1,RM1,RM1,R1,R1,
!     &     RM1,R1,RM1,RM1,RM1,R1,RM1,RM1,RM1/
      DATA THERMAL_CHECK_FLUCT/.FALSE./
      DATA N_TEMP_LAT/1/,TEMP_LAT/0D0/
      DATA FLUCT_DIR_SETTING/'not set   '/
      DATA MNT_TEMP_MCS/0D0/,TEMP_LAT_MIN/9999D0/,TEMP_LAT_MAX/9999D0/
      DATA NHAT_MAG/0.0D0,0.0D0,1.0D0/
C
      END
