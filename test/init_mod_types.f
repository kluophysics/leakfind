C*==init_mod_types.f    processed by SPAG 6.70Rc at 17:48 on 10 Jun 2016
      SUBROUTINE INIT_MOD_TYPES(NLSHELLMAX,FULLPOT,FULLPOT5,SPHERCELL)
C   ********************************************************************
C   *                                                                  *
C   *  initialize all tables dependent ONLY on                         *
C   *                                                                  *
C   *      NQ,  NT  and  NL   and derived array dimensions             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRNSMIN,NRMAX
      USE MOD_SITES,ONLY:VLMMAD_HOST,NQMAX,KFP_LMQ,CMNTQ
      USE MOD_TYPES,ONLY:CONC,ECORTAB,ECOR_LT,LOPT,NLT,NFPT,NLMFPT,NBLK,
     &    IMT,LCXRAY,NCXRAY,NCORT,NVALT,TXT_T,LTXT_T,NAT,Z,NSEMCORSHLT,
     &    SEMCORSHLT,SOCTL,CTL,IKMSOLBLK,NSOLBLK,ISOLIKM,KLMFP,LMIFP,
     &    RHOCHR,RHOCHRC,RHOSPN,RHOSPNC,RHOORB,VT,BT,VNST,BNST,RHO2NS,
     &    VMTZ,QEL,ABIT,AOPT,NTMAX,NLIN_T,NKM_T,NLMFPMAX,CMNTT,OBS_T,
     &    OBS_LT,OBS,OBS_TX,OBS_LTX,OBS_X,DOBS_T,DOBS_LT,DOBS,DOBS_TX,
     &    DOBS_LTX,DOBS_X,DOBS_BS_EF_TX,DOBS_BS_EF_LTX,OBS_T_GLO,
     &    OBS_GLO,OBS_TX_GLO,OBS_X_GLO,DOBS_T_GLO,DOBS_GLO,DOBS_TX_GLO,
     &    DOBS_X_GLO,MTET_T,MPHI_T,MGAM_T,X_CHEM_T
      USE MOD_ANGMOM,ONLY:NOBSMAX,NLMAX,NKMMAX,NLABIMAX,NLCORE
      USE MOD_CALCMODE,ONLY:BLCOUPL,ORBPOL,BREITINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FULLPOT,FULLPOT5,SPHERCELL
      INTEGER NLSHELLMAX
C
C Local variables
C
      INTEGER IT
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------- variables depending only on NT
C
      ALLOCATE (NAT(NTMAX),Z(NTMAX),LCXRAY(NTMAX),NCXRAY(NTMAX))
      ALLOCATE (ECORTAB(120,NTMAX),ECOR_LT(NLCORE,NTMAX))
      ALLOCATE (CONC(NTMAX),X_CHEM_T(NTMAX),IMT(NTMAX),LOPT(NTMAX))
      ALLOCATE (NCORT(NTMAX),NVALT(NTMAX),LTXT_T(NTMAX),TXT_T(NTMAX))
      ALLOCATE (NLT(NTMAX),NFPT(NTMAX),NLMFPT(NTMAX),NBLK(NTMAX))
      ALLOCATE (QEL(NTMAX))
      ALLOCATE (NLIN_T(NTMAX),NKM_T(NTMAX))
C
      ALLOCATE (OBS_T(0:3,NOBSMAX,NTMAX))
      ALLOCATE (OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX))
      ALLOCATE (OBS_TX(0:3,NOBSMAX,NTMAX))
      ALLOCATE (OBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX))
      ALLOCATE (OBS(0:3,NOBSMAX),OBS_X(0:3,NOBSMAX))
      ALLOCATE (DOBS_T(0:3,NOBSMAX,NTMAX))
      ALLOCATE (DOBS_LT(0:3,NOBSMAX,NLMAX,NTMAX))
      ALLOCATE (DOBS_TX(0:3,NOBSMAX,NTMAX))
      ALLOCATE (DOBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX))
      ALLOCATE (DOBS(0:3,NOBSMAX),DOBS_X(0:3,NOBSMAX))
      ALLOCATE (DOBS_BS_EF_TX(0:3,NOBSMAX,NTMAX))
      ALLOCATE (DOBS_BS_EF_LTX(0:3,NOBSMAX,NLMAX,NTMAX))
C
      ALLOCATE (OBS_T_GLO(0:3,NOBSMAX,NTMAX))
      ALLOCATE (OBS_TX_GLO(0:3,NOBSMAX,NTMAX))
      ALLOCATE (OBS_GLO(0:3,NOBSMAX),OBS_X_GLO(0:3,NOBSMAX))
      ALLOCATE (DOBS_T_GLO(0:3,NOBSMAX,NTMAX))
      ALLOCATE (DOBS_TX_GLO(0:3,NOBSMAX,NTMAX))
      ALLOCATE (DOBS_GLO(0:3,NOBSMAX),DOBS_X_GLO(0:3,NOBSMAX))
C
      NFPT = 9999999
      NLMFPT = 9999999
      CONC = 999999D0
      X_CHEM_T = 999999D0
      ECORTAB = 999999D0
      NBLK = 999999
      NLT = 999999
      NCORT = 999999
      LTXT_T = 999999
      VMTZ = 0.0D0
      QEL(:) = 0.0D0
      OBS_LT(:,:,:,:) = 0.0D0
      OBS_T(:,:,:) = 0.0D0
C
      DO IT = 1,NTMAX
         NCXRAY(IT) = 0
         LCXRAY(IT) = 0
         NVALT(IT) = 0
         LOPT(IT) = -1
      END DO
C
C-------------------------- variables specify type dependent local frame
C
      ALLOCATE (MTET_T(NTMAX),MPHI_T(NTMAX),MGAM_T(NTMAX))
      MTET_T(:) = 0.0D0
      MPHI_T(:) = 0.0D0
      MGAM_T(:) = 0.0D0
C
C------------------------------- variables depending on semi core states
C
      ALLOCATE (NSEMCORSHLT(NTMAX),SEMCORSHLT(NLSHELLMAX,NTMAX))
      SEMCORSHLT = 'XX'
      DO IT = 1,NTMAX
         NSEMCORSHLT(IT) = 0
      END DO
C
C-------------------------------------- variables depending on NT and NL
      ALLOCATE (CTL(NTMAX,NLMAX),SOCTL(NTMAX,NLMAX))
C
C--------------------------------- variables depending on NT and NL(POT)
C
      ALLOCATE (KLMFP(NLMFPMAX,NTMAX),LMIFP(NLMFPMAX,NTMAX))
      ALLOCATE (ISOLIKM(NKMMAX,NTMAX),NSOLBLK(NKMMAX,NTMAX))
      ALLOCATE (IKMSOLBLK(NKMMAX,NKMMAX,NTMAX))
      KLMFP = 999999
      LMIFP = 999999
      NSOLBLK = 999999
      ISOLIKM = 999999
      IKMSOLBLK = 999999
C
C-------------------------------------- variables depending on NT and NR
C
      ALLOCATE (RHOCHR(NRMAX,NTMAX),RHOCHRC(NRMAX,NTMAX))
      ALLOCATE (RHOSPN(NRMAX,NTMAX),RHOSPNC(NRMAX,NTMAX))
      ALLOCATE (RHOORB(NRMAX,NTMAX),VT(NRMAX,NTMAX),BT(NRMAX,NTMAX))
      RHOCHR = 999999D0
      RHOSPN = 999999D0
      RHOORB = 999999D0
C
C----------------------------- variables depending on NT, NR and NL(POT)
C
      IF ( FULLPOT .OR. FULLPOT5 .OR. SPHERCELL ) THEN
         ALLOCATE (VNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (BNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (RHO2NS(NRMAX,NLMFPMAX,NTMAX,3))
         VNST(:,:,:) = 0D0
         BNST(:,:,:) = 0D0
         RHO2NS(:,:,:,:) = 0D0
      ELSE
         ALLOCATE (VNST(1,1,1),BNST(1,1,1))
         ALLOCATE (RHO2NS(1,1,1,1))
      END IF
C
C------------------------------ variables depending on NT, NR and NL(OP)
C
      IF ( BREITINT ) THEN
         ALLOCATE (ABIT(NRMAX,NLABIMAX,-1:1,NTMAX))
         ABIT(:,:,:,:) = 0D0
      ELSE
         ALLOCATE (ABIT(1,1,1,1))
      END IF
C-----------------------------------------------------------------------
C
      IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL ) THEN
         ALLOCATE (AOPT(NRMAX,2,NTMAX))
      ELSE
         ALLOCATE (AOPT(1,1,1))
      END IF
C
C--------------------------------- variables depending on NL(POT) and NT
C
      ALLOCATE (CMNTT(NLMFPMAX,NTMAX))
      CMNTT(:,:) = 0D0
C
C--------------------------------- variables depending on NL(POT) and NQ
C
      ALLOCATE (VLMMAD_HOST(NLMFPMAX,NQMAX),KFP_LMQ(NLMFPMAX,NQMAX))
      ALLOCATE (CMNTQ(NLMFPMAX,NQMAX))
C
      END
