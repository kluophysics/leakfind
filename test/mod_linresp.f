C*==mod_kspace.f    processed by SPAG 6.70Rc at 10:38 on 26 Feb 2013
      MODULE MOD_LINRESP
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables connected with                    *
C   *                                                                  *
C   *                                                                  *
C   *    H_PERT_LMTP   angular expansion      of perturbation   IPERT  *
C   *                                                                  *
C   *    K_PERT_LMTP   flag indicating H_PERT_LMTP(r,LM,IRERT,IT) != 0 *
C   *                                                                  *
C   *    AMEG_LMOP     angular matrix element for operator      IOPER  *
C   *    AMEF_LMOP     for large (G) and small (F) component           *
C   *                  for IREL < 3 only AG is used                    *
C   *                                                                  *
C   *  linear response calculations                                    *
C   *                                                                  *
C   *    HA_ZZ(r) = Int_{0}^{r}       Z(r',E_a)  Delta H  Z^x(r',E_b)  *
C   *    HB_JZ(r) = Int_{r}^{r_crit}  J(r',E_a)  Delta H  Z^x(r',E_b)  *
C   *    HC_ZJ(r) = Int_{r}^{r_crit}  Z(r',E_a)  Delta H  J^x(r',E_b)  *
C   *    HD_JJ(r) = Int_{r}^{r_crit}  J(r',E_a)  Delta H  J^x(r',E_b)  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
      LOGICAL LINRESP_CHECK_SUM_RULES,LINRESP_CHECK_WITH_MEZZ
c modified by XJQ; not_store_chiz and use_jtauk_square, for big system, don't store huge matrix chiz
c                  prt_sig_realz, print sig(E) for kappa program using Kubo-Greenwood
      logical not_store_chiz, use_jtauk_square, prt_sig_realz
      integer ie0_sig, ie1_sig
c end-mod-xjq
      PARAMETER (LINRESP_CHECK_SUM_RULES=.FALSE.,
     &           LINRESP_CHECK_WITH_MEZZ=.FALSE.)
C
      INTEGER IDOS,ISPN,IORB,IHFI,ICDIA,IKDIA
      PARAMETER (IDOS=1,ISPN=2,IORB=3,IHFI=4,ICDIA=6,IKDIA=7)

      INTEGER NL_CHI,NLM_CHI
      PARAMETER (NL_CHI=2,NLM_CHI=9)
C
      INTEGER NPERT,NOPER,NOBSE,IOPER_PERT(:),IOPER_OBSE(:),
     &        NZ12,NZ12MAX

      REAL*8 GNTTAB(:,:,:)

      REAL*8 HZ_PERT_LMTP(:,:,:,:),HX_PERT_LMTP(:,:,:,:)
      COMPLEX*16 AMEG_LMOP(:,:,:,:),AMEF_LMOP(:,:,:,:)
      LOGICAL K_PERT_LMTP(:,:,:)

      COMPLEX*16 SHAZ_ZZ_T(:,:,:,:),SHBZ_JZ_T(:,:,:,:),
     &           SHDZ_JJ_T(:,:,:,:),SHCZ_ZJ_T(:,:,:,:)
      COMPLEX*16 HAZ_ZZ_T(:,:,:,:,:),HBZ_JZ_T(:,:,:,:,:),
     &           HDZ_JJ_T(:,:,:,:,:),HCZ_ZJ_T(:,:,:,:,:)
      COMPLEX*16 HAX_ZZ_T(:,:,:,:,:),HBX_JZ_T(:,:,:,:,:),
     &           HDX_JJ_T(:,:,:,:,:),HCX_ZJ_T(:,:,:,:,:)

      COMPLEX*16 FAZ_ZZ_T(:,:,:,:,:,:),FBZ_JZ_T(:,:,:,:,:,:),
     &           FDZ_JJ_T(:,:,:,:,:,:),FCZ_ZJ_T(:,:,:,:,:,:)
      COMPLEX*16 FAX_ZZ_T(:,:,:,:,:,:),FBX_JZ_T(:,:,:,:,:,:),
     &           FDX_JJ_T(:,:,:,:,:,:),FCX_ZJ_T(:,:,:,:,:,:)
      COMPLEX*16 FZ_ZBZA_T(:,:,:,:,:,:)

      COMPLEX*16 MZBZA_O(:,:,:),
     &           MIRR2_OP(:,:,:,:),MIRR3_OP(:,:,:,:),
     &           MIRR4_OP(:,:,:,:)

      COMPLEX*16 QQ1(:,:,:,:,:,:,:),
     &           QQ2(:,:,:,:,:,:,:),
     &           QQ3(:,:,:,:,:,:,:),
     &           QQ4(:,:,:,:,:,:,:)

      REAL*8 RHO2NS_GG(:,:,:,:,:)
      REAL*8 CHI_TO(:,:)

      COMPLEX*16 DOBS_TO(:,:)

      COMPLEX*16 THZT_DIA_P(:,:,:),THZT_OFF_TP(:,:,:,:),
     &           THXT_DIA_P(:,:,:),THXT_OFF_TP(:,:,:,:)      

      COMPLEX*16 D0ZL(:,:,:,:),   D0ZL1(:,:,:,:),
     &           D1ZL(:,:,:,:),   D1ZL1(:,:,:,:),
     &            DZL(:,:,:,:),    DZL1(:,:,:,:),
     &          DIJZL(:,:,:,:,:),DIJZL1(:,:,:,:,:)
      COMPLEX*16 D0Z (:,:,:),     D0Z1 (:,:,:),
     &           D1Z (:,:,:),     D1Z1 (:,:,:),
     &            DZ (:,:,:),      DZ1 (:,:,:),
     &          DIJZ (:,:,:,:),  DIJZ1 (:,:,:,:)

      COMPLEX*16 D0XL(:,:,:,:),  
     &          DIJXL(:,:,:,:,:)
      COMPLEX*16 D0X (:,:,:),  D1X (:,:,:),  DX (:,:,:),   
     &          DIJX (:,:,:,:)

      COMPLEX*16 T0ZL(:,:,:,:),  
     &           T1ZL(:,:,:,:),  
     &            TZL(:,:,:,:),  
     &          TIJZL(:,:,:,:,:)
      COMPLEX*16 T0Z (:,:,:),    
     &           T1Z (:,:,:),    
     &            TZ (:,:,:),    
     &          TIJZ (:,:,:,:)
      COMPLEX*16 TZ_STD (:,:,:)
      COMPLEX*16 TZ_DK (:,:,:,:)
      COMPLEX*16 T1Z_DK (:,:,:,:)
      COMPLEX*16 T0Z_DK (:,:,:,:)

      COMPLEX*16 T0XL(:,:,:,:),  
     &          TIJXL(:,:,:,:,:)
      COMPLEX*16 T0X (:,:,:),    T1X (:,:,:),     TX (:,:,:),    
     &          TIJX (:,:,:,:)     

      COMPLEX*16 CHIZ(:,:,:),CHI_DK(:,:,:)
      COMPLEX*16 DDTAUTAUT(:,:),TKTKTT(:,:,:)

      INTEGER ITTA(:),ITTB(:),ITTC(:),ITTD(:),ITTQ1(:),ITTQ2(:),ITTMAX,
     &        JTT1(:),JTT2(:),JTTX(:),JTTMAX,NTTJ(:),NTKTKLIN,
     &        NTKTKMAX
      COMPLEX*16 WTTJ(:)

      INTEGER NLINCHIMAX,NLIN23_CHI,NLIN41_CHI,
     &        IKM1_CHI_LIN(:),IKM2_CHI_LIN(:),
     &        IKM3_CHI_LIN(:),IKM4_CHI_LIN(:),
     &        NAB_CHI_QQ(:,:),NCOLROW(:,:),LAMCOLROW(:,:,:) 

      REAL*8 QVEC_PERT(3) 

      LOGICAL ERYDA_EQ_ERYDB, QVEC_PERT_EQ_0VEC
      LOGICAL K_TAU(:,:,:),K_PHI(:,:,:,:,:),K_QQ(:,:,:,:,:)

C---------------------------------------------------- perturbation terms
      ALLOCATABLE HZ_PERT_LMTP,HX_PERT_LMTP,K_PERT_LMTP

      ALLOCATABLE SHAZ_ZZ_T,SHBZ_JZ_T,SHCZ_ZJ_T,SHDZ_JJ_T
      ALLOCATABLE HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T
      ALLOCATABLE HAX_ZZ_T,HBX_JZ_T,HCX_ZJ_T,HDX_JJ_T

      ALLOCATABLE FAZ_ZZ_T,FBZ_JZ_T,FCZ_ZJ_T,FDZ_JJ_T
      ALLOCATABLE FAX_ZZ_T,FBX_JZ_T,FCX_ZJ_T,FDX_JJ_T
      ALLOCATABLE FZ_ZBZA_T
      ALLOCATABLE QQ1,QQ2,QQ3,QQ4,K_QQ

C---------------------------------------------- angular matrix elkements
      ALLOCATABLE AMEG_LMOP,AMEF_LMOP

C-------------- pointers between operators and perturbation / observable
      ALLOCATABLE IOPER_PERT,IOPER_OBSE

C--------------------------------------------------- response quantities
      ALLOCATABLE RHO2NS_GG

      ALLOCATABLE MZBZA_O,MIRR2_OP,MIRR3_OP,MIRR4_OP

C-----------------------------------------------------------------------
C   variables used for BZ-integration Int d^3k Tau(k,E_a)*Tau(k+q,E_b)
C-----------------------------------------------------------------------

C------------------------------------------- variables depending on NZ12
      ALLOCATABLE CHIZ
C--------------------------------------------- variables depending on DK
      ALLOCATABLE CHI_DK
C------------------------------------------ variables depending on NTKTK
      ALLOCATABLE DDTAUTAUT,TKTKTT
C----------------------------------------- variables depending on ITTMAX
      ALLOCATABLE ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ
C----------------------------------------- variables depending on JTTMAX
      ALLOCATABLE JTT1,JTT2,JTTX,WTTJ
C---------------------------------------- variables depending on NLINCHI
      ALLOCATABLE IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN
C--------------------------------------------- variables depending on NQ
      ALLOCATABLE NAB_CHI_QQ
C-------------------------------------- variables depending on NQ and NL
      ALLOCATABLE LAMCOLROW,NCOLROW

C------------------------------------------- E-densities for observables
      ALLOCATABLE DOBS_TO
      ALLOCATABLE CHI_TO

C---------------------------------------------------- Gaunt coefficients
      ALLOCATABLE GNTTAB

C----------------------------------- combination of TAU and perturbation
      ALLOCATABLE THZT_OFF_TP,THXT_OFF_TP,THZT_DIA_P,THXT_DIA_P

      ALLOCATABLE D0ZL,   D0ZL1,
     &           D1ZL,   D1ZL1,
     &            DZL,    DZL1,
     &          DIJZL,DIJZL1
      ALLOCATABLE D0Z ,     D0Z1 ,
     &           D1Z ,     D1Z1 ,
     &            DZ ,      DZ1 ,
     &          DIJZ ,  DIJZ1 

      ALLOCATABLE D0XL,  
     &          DIJXL
      ALLOCATABLE D0X , D1X ,    DX ,    
     &          DIJX 

      ALLOCATABLE T0ZL,  
     &           T1ZL,  
     &            TZL,  
     &          TIJZL
      ALLOCATABLE T0Z ,    
     &           T1Z ,    
     &            TZ ,    
     &          TIJZ ,
     &          TZ_STD,
     &          TZ_DK,
     &          T0Z_DK,
     &          T1Z_DK

      ALLOCATABLE T0XL,  
     &          TIJXL
      ALLOCATABLE T0X ,    T1X ,   TX ,    
     &          TIJX      

      ALLOCATABLE K_TAU,K_PHI

      DATA NZ12/999999/,NZ12MAX/999999/
      DATA NPERT/1/,NOPER/1/,NOBSE/1/
      DATA QVEC_PERT_EQ_0VEC/.TRUE./

      SAVE NPERT,NOPER,NOBSE,IOPER_PERT,IOPER_OBSE,
     &   SHAZ_ZZ_T,SHBZ_JZ_T,SHCZ_ZJ_T,SHDZ_JJ_T,
     &   HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T,
     &   HAX_ZZ_T,HBX_JZ_T,HCX_ZJ_T,HDX_JJ_T,
     &   HZ_PERT_LMTP,HX_PERT_LMTP,K_PERT_LMTP,AMEG_LMOP,AMEF_LMOP,
     &   CHIZ,CHI_DK,DDTAUTAUT,TKTKTT,NZ12,NZ12MAX,ERYDA_EQ_ERYDB,
     &   ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,
     &   JTT1,JTT2,JTTX,WTTJ,ITTMAX,JTTMAX,NTKTKLIN,NTKTKMAX,
     &   NLINCHIMAX,NLIN23_CHI,NLIN41_CHI,
     &   IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN,
     &   NAB_CHI_QQ,LAMCOLROW,NCOLROW,RHO2NS_GG, 
     &   MZBZA_O,MIRR2_OP,MIRR3_OP,MIRR4_OP,
     &   DOBS_TO,CHI_TO
      SAVE FAZ_ZZ_T,FBZ_JZ_T,FCZ_ZJ_T,FDZ_JJ_T,
     &     FAX_ZZ_T,FBX_JZ_T,FCX_ZJ_T,FDX_JJ_T,FZ_ZBZA_T
      SAVE QQ1,QQ2,QQ3,QQ4,K_QQ

      SAVE GNTTAB

      SAVE QVEC_PERT_EQ_0VEC,QVEC_PERT

      SAVE THZT_OFF_TP,THXT_OFF_TP,THZT_DIA_P,THXT_DIA_P

      SAVE D0ZL,   D0ZL1,
     &           D1ZL,   D1ZL1,
     &            DZL,    DZL1,
     &          DIJZL,DIJZL1
      SAVE D0Z ,     D0Z1 ,
     &           D1Z ,     D1Z1 ,
     &            DZ ,      DZ1 ,
     &          DIJZ ,  DIJZ1 

      SAVE D0XL,  
     &          DIJXL
      SAVE D0X ,  D1X ,     DX ,    
     &          DIJX 

      SAVE T0ZL,  
     &           T1ZL,  
     &            TZL,  
     &          TIJZL
      SAVE T0Z ,    
     &           T1Z ,    
     &            TZ ,   
     &          TIJZ ,
     &          TZ_STD,
     &          TZ_DK,
     &          T1Z_DK,
     &          T0Z_DK

      SAVE T0XL,  
     &          TIJXL
      SAVE T0X ,    T1X ,   TX ,     
     &          TIJX      

      SAVE K_TAU,K_PHI

      END
