C*==mod_angmom.f    processed by SPAG 6.70Rc at 10:07 on 24 Jan 2015
      MODULE MOD_ANGMOM
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tables dependent ONLY on                    *
C   *                                                                  *
C   *      NLMAX  and  NQMAX   and derived array dimensions            *
C   *                                                                  *
C   *  ALL data initialized by  <INIT_MOD_ANGMOM>                      *
C   *                                                                  *
C   ********************************************************************
C   *       relativistic MATRIX ELEMENTS and OBSERVABLES               *
C   ********************************************************************
C   *                                                                  *
C   *   1: IDOS  < LAM |      1      | LAM' >    (l,s)-resolved DOS    *
C   *   2: ISMT  < LAM | sigma(ipol) | LAM' >    spin moment           *
C   *   3: IOMT  < LAM |     l(ipol) | LAM' >    orbital moment        *
C   *   4: IHFF  < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
C   *   5: ISDM  < LAM |     T(ipol) | LAM' >    spin dipole moment    *
C   *   6: IKDS              (dummy)             kappa-resolved DOS    *
C   *   7: IBND  < LAM |      1      | LAM' >    band energy           *
C   *   8: IODN  < LAM |P_dn l(ipol) | LAM' >    orb. mnt. spin down   *
C   *   9: IOUP  < LAM |P_up l(ipol) | LAM' >    orb. mnt. spin up     *
C   *  10: ITRQ  < LAM |     t(ipol) | LAM' >    torque operator       *
C   *  11: IWFD                                  Weiss field           *
C   *                                                                  *
C   *                    SCF        ELSE         array size for:       *
C   *   NMEMAX   =        4           6          MEZZ, MEZJ            *
C   *   NOBSMAX  =       12          12          OBS_T, OBS_LT, etc.   *
C   *                                                                  *
C   *   IPOL= 1,2,3  ==  -1,0,+1  ==  (-),(z),(+)    NEW               *
C   *                                                                  *
C   *   IPOL= 1,2,3  ==  +1,-1,z                     OLD               *
C   *   !!!!!!!!! STILL IN USE FOR SPECTROSCOPY !!!!!!!!!!!!!!!!!!!!!! *
C   *                                                                  *
C   *   IKM = 2 * l * (j+1/2) + j + mj + 1                             *
C   *                                                                  *
C   *   list revised with version 7.4.0  Sept 2015                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NMEMAX,NMEHFMAX,NMVECMAX,NPOL,NPOLMAX,NOBSMAX,NLCORE
      PARAMETER (NMEMAX=9,NMEHFMAX=5,NMVECMAX=4,
     &           NOBSMAX=12,NPOL=3,NPOLMAX=NPOL)
      CHARACTER*4 GLO_TO_LOC,LOC_TO_GLO
      PARAMETER (GLO_TO_LOC='SAS+',LOC_TO_GLO='S+AS')
      INTEGER IDOS,ISMT,IOMT,IHFF,IBND,ISDM,IODN,IOUP,IKDS,ITRQ,IWFD
      PARAMETER (IDOS=1,ISMT=2,IOMT=3,IHFF=4,ISDM=5,IKDS=6,
     &           IBND=7,IODN=8,IOUP=9,ITRQ=10,IWFD=11)
      PARAMETER (NLCORE=4)
C
C Local variables
C
      REAL*8 A1_ADA(:,:,:),A2_ADA(:,:,:),AB1_GRV(:,:,:),AB1_NAB(:,:,:),
     &       AB2_GRV(:,:,:),AB2_NAB(:,:,:),AB_GRV(:,:,:),AC1_GRV(:,:,:),
     &       AC2_GRV(:,:,:),AMEBI1(:,:,:,:),AMEBI2(:,:,:,:),
     &       AMEOPC(:,:,:),AMEOPO(:,:,:),AME_F(:,:,:,:),AME_G(:,:,:,:),
     &       AN1_NAB(:,:,:),AN2_NAB(:,:,:),AV1_NAB(:,:,:),AV2_NAB(:,:,:)
     &       ,AV_GRV(:,:,:),A_RSIG(:,:,:,:),A_SIGMA(:,:,:),
     &       B1_ADA(:,:,:,:),B2_ADA(:,:,:,:),
     &       CGC(:,:),CGCFP(:,:),DBLFACT_L(:),FACT(:),
     &       GBIG(:,:,:,:,:),
     &       GBIL(:,:,:,:,:),QDIA(:),
     &       QMDIA(:),QMOFF(:),QOFF(:),SDIA(:),SMDIA(:),SMOFF(:),SOFF(:)
     &       ,AF1_SPIN_CURR_NAB(:,:,:,:),AF2_SPIN_CURR_NAB(:,:,:,:),
     &       AG1_SPIN_CURR_NAB(:,:,:,:),AG2_SPIN_CURR_NAB(:,:,:,:)
      COMPLEX*16 A1_SPIN_CURR_ALF(:,:,:,:),A2_SPIN_CURR_ALF(:,:,:,:),
     &           AF_RGNT(:,:,:),AG_RGNT(:,:,:),AME_RLM(:,:,:,:,:),
     &           A_RGNT(:,:,:),A_SIGY(:,:,:,:),A_SIG_RLM(:,:,:,:),
     &           A_Y1M(:,:,:),CREL(:,:),FPAB1_NAB(:,:,:,:,:),
     &           FPAB2_NAB(:,:,:,:,:),FPAN1_NAB(:,:,:,:,:),
     &           FPAN2_NAB(:,:,:,:,:),FPAV1_NAB(:,:,:,:,:),
     &           FPAV2_NAB(:,:,:,:,:),
     &           LMAT2(:,:),LMAT3(:,:),MEZJ(:,:,:,:)
     &           ,MEZZ(:,:,:,:),MIRR_2(:,:,:,:),MIRR_3(:,:,:,:),
     &           MIRR_4(:,:,:,:),MSSQ(:,:,:),MSST(:,:,:),MZAZB(:,:,:),
     &           MZBZA(:,:,:),RC(:,:),RREL(:,:),SSST(:,:,:),TAUQ(:,:,:),
     &           TAUT(:,:,:),TSSQ(:,:,:),TSST(:,:,:),U_CS(3,3),U_SC(3,3)
     &           ,WKM1(:,:),WKM2(:,:),WKM3(:,:),WKM4(:,:),WKM_LPK(:,:),
     &           WLM1(:,:),WLM2(:,:),WLM3(:,:),WXM1(:,:),WXM2(:,:),
     &           WXM3(:,:),WXM4(:,:)
      COMPLEX*16 FPA1_ADA(:,:,:,:,:),FPA2_ADA(:,:,:,:,:),
     &           FPB1_ADA(:,:,:,:,:),FPB2_ADA(:,:,:,:,:)
      INTEGER DELTA_L_EXT,IKM1LIN(:),IKM2LIN(:),IKMCPL_KM(:),IKMLLIM1(:)
     &        ,IKMLLIM2(:),IKM_KAP_MUEM05(:,:),ILMBOT_LM(:),ILMTOP_LM(:)
     &        ,IMKM_IKM(:),IND0Q(:),IPIVKM(:),IXM0_QCLU(:),IXM0_QTB(:),
     &        JP05_IKM(:),KAPPA_IKM(:),KAPTAB(:),LBTAB(:),LB_IKM(:),
     &        LINMAX,LTAB(:),L_IKM(:),L_LM(:),L_LMS(:),MUEM05_IKM(:),
     &        M_LM(:),M_LMS(:),NCGNT123TAB(18),NCPLWF(:),NCPLWF_LA(:),
     &        NCPLWF_LB(:),NCPLWF_RA(:),NCPLWF_RB(:),NK,NKKR,NKKR_CLU,
     &        NKM,NKMAX,NKMCPL_KM(:),NKMMAX,NKMPMAX,NKMP_EXT,NKMQ(:),
     &        NKM_EXT,NKQ(:),NK_EXT,NL,NLABIMAX,NLAMEFPMAX,NLIN,NLINQ(:)
     &        ,NLM,NLMAX,NLMMADMAX,NLMMAX,NLMQ(:),NLM_AME_RLM_EXT,
     &        NLM_EXT,NLQ(:),NL_AME_RLM_EXT,NL_EXT,NMUEMAX,NMUETAB(:),
     &        NRGNT123TAB(20),NSOLLM(:,:),NSPIN,NXM,NXMMAX,NXM_Q(:),
     &        NXM_QCLU(:),NXM_QTB(:)
      LOGICAL K_AME(:,:),K_AME_RLM(:,:,:),K_AME_RLM_Z(:,:,:),
     &        K_AME_Z(:,:)
      CHARACTER*3 TXT_OBS(NOBSMAX)
      CHARACTER*4 TXT_J(:)
      CHARACTER*1 TXT_L(:)
      SAVE A1_ADA,A1_SPIN_CURR_ALF,A2_ADA,A2_SPIN_CURR_ALF,AB1_GRV,
     &     AB1_NAB,AB2_GRV,AB2_NAB,AB_GRV,AC1_GRV,AC2_GRV,
     &     AF1_SPIN_CURR_NAB,AF2_SPIN_CURR_NAB,
     &     AF_RGNT,AG1_SPIN_CURR_NAB,
     &     AG2_SPIN_CURR_NAB,
     &     AG_RGNT,AMEBI1,AMEBI2,AMEOPC,AMEOPO,AME_F,AME_G,AME_RLM,
     &     AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,AV_GRV,A_RGNT,A_RSIG,A_SIGMA,
     &     A_SIGY,A_SIG_RLM,A_Y1M,B1_ADA,B2_ADA,CGC,CGCFP,CREL,
     &     DBLFACT_L,FACT,
     &     FPAB1_NAB,FPAB2_NAB,FPAN1_NAB,FPAN2_NAB,FPAV1_NAB,FPAV2_NAB,
     &     GBIG,GBIL,IKM1LIN,
     &     IKM2LIN,IKMCPL_KM,IKMLLIM1,IKMLLIM2,IKM_KAP_MUEM05,ILMBOT_LM,
     &     ILMTOP_LM,IMKM_IKM,IND0Q,IPIVKM,IXM0_QCLU,IXM0_QTB,JP05_IKM,
     &     KAPPA_IKM,KAPTAB,K_AME,K_AME_RLM,K_AME_RLM_Z,K_AME_Z,LBTAB,
     &     LB_IKM,LINMAX,LMAT2,LMAT3,LTAB,L_IKM,L_LM,L_LMS,MEZJ,MEZZ,
     &     MIRR_2,MIRR_3,MIRR_4,MSSQ,MSST,MUEM05_IKM,MZAZB,MZBZA,M_LM,
     &     M_LMS,NCPLWF,NCPLWF_LA,NCPLWF_LB,NCPLWF_RA,NCPLWF_RB,NK,NKKR,
     &     NKKR_CLU,NKM,NKMAX,NKMCPL_KM,NKMMAX,NKMPMAX,NKMP_EXT,NKMQ,
     &     NKM_EXT,NKQ,NK_EXT,NL,NLABIMAX,NLIN,NLINQ,NLM,NLMAX,
     &     NLMMADMAX,NLMMAX,NLMQ,NLM_AME_RLM_EXT,NLM_EXT,NLQ,
     &     NL_AME_RLM_EXT,NL_EXT,NMUEMAX,NMUETAB,NSPIN,NXM,NXMMAX,NXM_Q,
     &     NXM_QCLU,NXM_QTB,QDIA,QMDIA,QMOFF,QOFF,RC,RREL,SDIA,SMDIA,
     &     SMOFF,SOFF,SSST,TAUQ,TAUT,TSSQ,TSST,TXT_J,TXT_L,U_CS,U_SC,
     &     WKM1,WKM2,WKM3,WKM4,WKM_LPK,WLM1,WLM2,WLM3,WXM1,WXM2,WXM3,
     &     WXM4

      SAVE FPA1_ADA,FPA2_ADA,FPB1_ADA,FPB2_ADA
C
C*** End of declarations rewritten by SPAG
C
      DATA NCGNT123TAB/1,19,126,492,1453,3503,7450,14382,25655,43109,
     &     69034,106066,157607,227393,319920,440346,594733,789521/
C
      DATA NRGNT123TAB/1,15,96,388,1181,2917,6342,12452,22525,38289,
     &     61912,95914,143531,208371,294744,407644,552931,736829,966544,
     &     1250346/
      DATA NLAMEFPMAX/1/
      DATA DELTA_L_EXT/1/
      DATA TXT_OBS/'CHR','SMT','OMT','HFF','SDM','KDS',
     &             'BND','ODN','OUP','TRQ','WFD','***'/
C
      ALLOCATABLE FPA1_ADA,FPA2_ADA,FPB1_ADA,FPB2_ADA

      ALLOCATABLE TAUT,TSST,MSST,SSST
      ALLOCATABLE TAUQ,TSSQ,MSSQ
      ALLOCATABLE MEZJ,MEZZ
      ALLOCATABLE K_AME,K_AME_RLM,K_AME_Z,K_AME_RLM_Z
C
C---------------------------------------- variables depending only on NQ
C
      ALLOCATABLE NLQ,NKQ,NKMQ,NLMQ,IND0Q,NLINQ
      ALLOCATABLE IXM0_QCLU,NXM_QCLU
      ALLOCATABLE IXM0_QTB,NXM_QTB,NXM_Q
C
C---------------------------------------- variables depending only on NL
C
      ALLOCATABLE SDIA,SOFF,SMDIA,SMOFF,QDIA,QOFF,QMDIA,QMOFF
      ALLOCATABLE TXT_L,TXT_J,IKM1LIN,IKM2LIN,IKMLLIM1,IKMLLIM2,IMKM_IKM
      ALLOCATABLE IKM_KAP_MUEM05,L_LM,M_LM,L_LMS,M_LMS
      ALLOCATABLE KAPPA_IKM,MUEM05_IKM,L_IKM,JP05_IKM,LB_IKM
      ALLOCATABLE IKMCPL_KM,NKMCPL_KM
      ALLOCATABLE CGC,CGCFP,ILMBOT_LM,ILMTOP_LM
      ALLOCATABLE KAPTAB,LTAB,LBTAB,NMUETAB,NSOLLM,NCPLWF,FACT,DBLFACT_L
      ALLOCATABLE CREL,RC,RREL,A_RGNT,AG_RGNT,AF_RGNT
      ALLOCATABLE A_SIGMA,A_RSIG,A_SIGY,A_SIG_RLM
      ALLOCATABLE NCPLWF_RA,NCPLWF_LA,NCPLWF_RB,NCPLWF_LB
C
C----------------------------- angular matrix elements for NME operators
C
      ALLOCATABLE AME_G,AME_F,A_Y1M,AME_RLM
C
C-------------------------------- matrix elements for NPOL polarisations
C
      ALLOCATABLE MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4
C
C
C --------------- angular matrix elements for Brooks and Breit formalism
C
      ALLOCATABLE AMEOPC,AMEOPO,AMEBI1,AMEBI2,GBIG,GBIL
C
C ----------------------------- angular matrix elements for spectroscopy
C
      ALLOCATABLE A1_ADA,A2_ADA,B1_ADA,B2_ADA,AN1_NAB,AN2_NAB
      ALLOCATABLE AV1_NAB,AV2_NAB,AB1_NAB,AB2_NAB,AV_GRV,AB_GRV
      ALLOCATABLE AB1_GRV,AB2_GRV,AC1_GRV,AC2_GRV
      ALLOCATABLE FPAN1_NAB
      ALLOCATABLE FPAN2_NAB,FPAV1_NAB,FPAV2_NAB,FPAB1_NAB,FPAB2_NAB
      ALLOCATABLE A1_SPIN_CURR_ALF,A2_SPIN_CURR_ALF
      ALLOCATABLE AG1_SPIN_CURR_NAB,AG2_SPIN_CURR_NAB
      ALLOCATABLE AF1_SPIN_CURR_NAB,AF2_SPIN_CURR_NAB
C
C ----------------------------------- workspace for matrix inversion etc
C    to avoid conflicts use WKM_LPK ONLY as argument for LAPACK routines
C
      ALLOCATABLE WKM1,WKM2,WKM3,WKM4,WLM1,WLM2,WLM3,WKM_LPK,IPIVKM
      ALLOCATABLE WXM1,WXM2,WXM3,WXM4,LMAT2,LMAT3
C
      END
