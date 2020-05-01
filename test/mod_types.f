C*==mod_types.f    processed by SPAG 6.70Rc at 18:51 on 23 Feb 2014
      MODULE MOD_TYPES
C   ********************************************************************
C   *                                                                  *
C   *  module to store all type information                            *
C   *                                                                  *
C   *  ALL data initialized in  <<MAINPROGRAM>>                        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLSHELLMAX,NCORMAX
      PARAMETER (NLSHELLMAX=15,NCORMAX=110)
C
C Local variables
C
      COMPLEX*16 DROT_T(:,:,:)
      REAL*8 ABIT(:,:,:,:),AOPT(:,:,:),BEXT,BNST(:,:,:),BT(:,:),
     &       CMNTT(:,:),CONC(:),X_CHEM_T(:),CTL(:,:),ECORTAB(:,:),ETOT,
     &       JDNST(:,:,:,:),MUEORB,MUESPN,NVALTOT,MROT_T(:,:,:),
     &       QEL(:),RHO2NS(:,:,:,:),RHOCHR(:,:),RHOCHRC(:,:),
     &       RHOORB(:,:),RHOSPN(:,:),RHOSPNC(:,:),SOCTL(:,:),
     &       VMTZ,VNST(:,:,:),VT(:,:),X_TCLU(:),
     &       JALF_LMCT(:,:,:,:),JORB_LMCT(:,:,:,:),
     &       OBS_T (:,:,:), OBS_LT (:,:,:,:), OBS  (:,:),
     &       OBS_TX(:,:,:), OBS_LTX(:,:,:,:), OBS_X(:,:),TOTDOS_BS_EF,
     &       OBS_T_GLO (:,:,:), OBS_GLO  (:,:), 
     &       OBS_TX_GLO(:,:,:), OBS_X_GLO(:,:), 
     &       W_WEISS_T(:),ECOR_LT(:,:),MTET_T(:),MPHI_T(:),MGAM_T(:),
     &       DOBS_LTEX(:,:,:,:,:),DOBS_TEX_GLO(:,:,:,:)
      INTEGER IKMCPLWF(:,:),IKMCPLWF_LA(:,:),IKMCPLWF_LB(:,:),
     &        IKMCPLWF_RA(:,:),IKMCPLWF_RB(:,:),IKMSOLBLK(:,:,:),IMT(:),
     &        ISOLIKM(:,:),ITBOT,ITTOP,KLMFP(:,:),LCXRAY(:),LMIFP(:,:),
     &        LOPT(:),LTXT_T(:),LTXT_TCLU(:),NAT(:),NA_TCLU(:),NBLK(:),
     &        NCORT(:),NCOR_TCLU(:),NCPLWFMAX,NCXRAY(:),NFPT(:),NKM_T(:)
     &        ,NLAFPMAX,NLFP,NLFPMAX,NLIN_T(:),NLMFP,NLMFPMAX,NLMFPT(:),
     &        NLT(:),NPOTMAX,NSEMCORSHLT(:),NSEMCORSHL_TCLU(:),
     &        NSOLBLK(:,:),NT,NTCLU,NTCLUMAX,NTHOST,NTMAX,NT_L,NVALT(:),
     &        NVAL_TCLU(:),Z(:),Z_TCLU(:),
     &        BISNST(:,:,:),BIST(:,:),VISNST(:,:,:),VIST(:,:)
      CHARACTER*2 SEMCORSHLT(:,:)
      CHARACTER*8 TXT_T(:),TXT_TCLU(:)
      COMPLEX*16 VAMEF(:,:,:,:,:,:),VAMEG(:,:,:,:,:,:)
      COMPLEX*16 JFLA(:,:,:),
     &           JGLA(:,:,:),
     &           ZFLA(:,:,:),
     &           ZGLA(:,:,:)
      COMPLEX*16 JFLB(:,:,:),JFRA(:,:,:),
     &           JGLB(:,:,:),JGRA(:,:,:),
     &           ZFLB(:,:,:),ZFRA(:,:,:),
     &           ZGLB(:,:,:),ZGRA(:,:,:),
     &       DOBS_T (:,:,:), DOBS_LT  (:,:,:,:), DOBS  (:,:),
     &       DOBS_TX(:,:,:), DOBS_LTX (:,:,:,:), DOBS_X(:,:),
     &       DOBS_T_GLO (:,:,:), DOBS_GLO  (:,:),
     &       DOBS_TX_GLO(:,:,:), DOBS_X_GLO(:,:),
     &           DOBS_BS_EF_TX(:,:,:),
     &           DOBS_BS_EF_LTX(:,:,:,:)
      SAVE ABIT,AOPT,BNST,BT,CMNTT,CONC,X_CHEM_T,CTL,ECORTAB,ECOR_LT,
     &     IKMCPLWF,IKMCPLWF_LA,IKMCPLWF_LB,IKMCPLWF_RA,IKMCPLWF_RB,
     &     IKMSOLBLK,IMT,ISOLIKM,ITBOT,ITTOP,JDNST,KLMFP,LCXRAY,LMIFP,
     &     LOPT,LTXT_T,LTXT_TCLU,NAT,NA_TCLU,NBLK,NCORT,NCOR_TCLU,
     &     NCXRAY,NFPT,NKM_T,NLIN_T,NLMFPT,NLT,NSEMCORSHLT,
     &     NSEMCORSHL_TCLU,NSOLBLK,NT,NTCLUMAX,NTMAX,NT_L,NVALT,NVALTOT,
     &     NVAL_TCLU,QEL,RHO2NS,RHOCHR,RHOCHRC,RHOORB,RHOSPN,
     &     RHOSPNC,SEMCORSHLT,SOCTL,TXT_T,TXT_TCLU,VAMEF,VAMEG,
     &     VMTZ,VNST,VT,X_TCLU,Z,Z_TCLU,JALF_LMCT,JORB_LMCT,
     &     W_WEISS_T
      SAVE JGLA,JFLA,ZFLA,ZGLA
      SAVE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
      SAVE BISNST,BIST,VISNST,VIST
      SAVE OBS_T  , OBS_LT   , OBS   ,
     &     OBS_TX , OBS_LTX  , OBS_X ,
     &     DOBS_T , DOBS_LT  , DOBS  ,
     &     DOBS_TX, DOBS_LTX , DOBS_X
      SAVE OBS_T_GLO  , OBS_GLO   ,
     &     OBS_TX_GLO , OBS_X_GLO ,
     &     DOBS_T_GLO , DOBS_GLO ,
     &     DOBS_TX_GLO, DOBS_X_GLO
      SAVE DOBS_LTEX,DOBS_TEX_GLO
      SAVE DOBS_BS_EF_TX,DOBS_BS_EF_LTX,TOTDOS_BS_EF
      SAVE DROT_T,MROT_T,MTET_T,MPHI_T,MGAM_T
C     
C*** End of declarations rewritten by SPAG

C---------------------------------------- variables depending only on NT
C
      ALLOCATABLE CONC,X_CHEM_T,ECORTAB,ECOR_LT
      ALLOCATABLE LOPT,NLT,NFPT,NLMFPT,NBLK
      ALLOCATABLE NKM_T,NLIN_T
      ALLOCATABLE IMT,LCXRAY,NCXRAY,NCORT,NVALT,TXT_T,LTXT_T,NAT,Z,QEL
      ALLOCATABLE X_TCLU,NVAL_TCLU,LTXT_TCLU,TXT_TCLU,Z_TCLU,NCOR_TCLU
      ALLOCATABLE NA_TCLU
      ALLOCATABLE W_WEISS_T
C     
C------------------------------- variables depending on semi core states
C
      ALLOCATABLE NSEMCORSHLT,SEMCORSHLT,NSEMCORSHL_TCLU
C
C-------------------------------------- variables depending on NT and NL
C
      ALLOCATABLE SOCTL,CTL,IKMSOLBLK,NSOLBLK,ISOLIKM,CMNTT
      ALLOCATABLE JALF_LMCT,JORB_LMCT
      ALLOCATABLE OBS_T,OBS_LT,OBS,OBS_TX,OBS_LTX,OBS_X
      ALLOCATABLE DOBS_T,DOBS_LT,DOBS
      ALLOCATABLE DOBS_TX,DOBS_LTX,DOBS_X
      ALLOCATABLE DOBS_BS_EF_TX,DOBS_BS_EF_LTX
      ALLOCATABLE OBS_T_GLO, OBS_GLO
      ALLOCATABLE OBS_TX_GLO,OBS_X_GLO
      ALLOCATABLE DOBS_T_GLO, DOBS_GLO
      ALLOCATABLE DOBS_TX_GLO,DOBS_X_GLO
      ALLOCATABLE DOBS_LTEX,DOBS_TEX_GLO
C     
C--------------------------------- variables depending on NT and NL(POT)
C
      ALLOCATABLE KLMFP,LMIFP,IKMCPLWF,VAMEF,VAMEG
C
C-------------------------------------- variables depending on NT and NR
C
      ALLOCATABLE RHOCHR,RHOCHRC,RHOSPN,RHOSPNC,RHOORB,VT,BT
C
C----------------------------- variables depending on NT, NR and NL(POT)
C
      ALLOCATABLE VNST,BNST,RHO2NS
      ALLOCATABLE BISNST,BIST,VISNST,VIST
C
C------------------------------ variables depending on NT, NR and NL(OP)
      ALLOCATABLE AOPT,ABIT,JDNST
C
C-------------------------------- wave functions and related quantities
C
      ALLOCATABLE JGLA,JFLA,ZFLA,ZGLA
      ALLOCATABLE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
      ALLOCATABLE IKMCPLWF_RA,IKMCPLWF_LA,IKMCPLWF_RB,IKMCPLWF_LB
      ALLOCATABLE DROT_T,MROT_T,MTET_T,MPHI_T,MGAM_T
C
      DATA NLAFPMAX/1/,NPOTMAX/1/,NCPLWFMAX/1/BEXT/0D0/,NTCLU/0/,
     &     NTHOST/999999/
      DATA ETOT/0D0/,MUEORB/0D0/,MUESPN/0D0/
      DATA NLFP/0/,NLMFP/0/,NLFPMAX/0/,NLMFPMAX/0/
C
      END
