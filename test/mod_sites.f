C*==mod_sites.f    processed by SPAG 6.70Rc at 17:25 on 17 Sep 2012
      MODULE MOD_SITES
C   ********************************************************************
C   *                                                                  *
C   *  module to store all lattice information                         *
C   *                                                                  *
C   *  ALL data initialized in  <<MAINPROGRAM>>                        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      COMPLEX*16 ABIMAD(:,:,:,:),DROTQ(:,:,:)
      REAL*8 ALFDEG,AVMAD(:,:,:,:),AVMAD_LI(:,:,:),AVMAD_LL(:,:,:),
     &       AVMAD_LR(:,:,:),AVMAD_RI(:,:,:),AVMAD_RL(:,:,:),
     &       AVMAD_RR(:,:,:),BETDEG,CMNTQ(:,:),DQBAS_QCLU(:,:),
     &       FORCE_QCLU(:,:),GAMDEG,MDIRQ(:,:),MROTQ(:,:,:),QBAS(:,:),
     &       QBAS0_QCLU(:,:),QBAS_QCLU(:,:),QMGAM(:),QMGAM_POTFIL(:),
     &       QMPHI(:),QMPHI_POTFIL(:),QMTET(:),QMTET_POTFIL(:),QMVEC(3),
     &       Q_LMQ(:,:),RNNQ(:),SHIFT_CLU,VLMMAD_BACK(:,:),
     &       VLMMAD_HOST(:,:),VMAD2D_A,VMAD2D_B
      INTEGER ICPA(:),IMQ(:),IM_QCLU(:),IQAT(:,:),IQBOT,IQBOT_CHI,
     &        IQBOT_CLU,IQBOT_HOST,IQBOT_L,IQBOT_R,IQBOT_TB,
     &        IQCLU_ATCLU(:,:),IQTOP,IQTOP_CHI,IQTOP_CLU,IQTOP_HOST,
     &        IQTOP_TB,IQ_LREF,IQ_QCLU(:),IQ_RREF,ITCLU_OQCLU(:,:),
     &        ITOQ(:,:),IT_OQCLU(:,:),KFP_LMQ(:,:),N5VEC_QCLU(:,:),
     &        NLMAD,NLMMAD,NLMQMAD,NLMVMAD,NLQMAD,NLVMAD,NNNQ(:),NOQ(:),
     &        NO_QCLU(:),NQ,NQCLU,NQCLUMAX,NQHOST,NQMAX,NQTB,NQ_CHI,
     &        NQ_I,NQ_L,NQ_R,NOMAX
      LOGICAL KMAD1M_Q(:),MAGROT_Q(:),OPTIMIZE_BASIS_GEOMETRY
      SAVE ABIMAD,AVMAD,AVMAD_LI,AVMAD_LL,AVMAD_LR,AVMAD_RI,AVMAD_RL,
     &     AVMAD_RR,CMNTQ,DQBAS_QCLU,DROTQ,FORCE_QCLU,ICPA,IMQ,IM_QCLU,
     &     IQAT,IQBOT_HOST,IQCLU_ATCLU,IQTOP,IQTOP_CLU,IQTOP_HOST,
     &     IQ_LREF,IQ_QCLU,IQ_RREF,ITCLU_OQCLU,ITOQ,IT_OQCLU,KFP_LMQ,
     &     KMAD1M_Q,MAGROT_Q,MDIRQ,MROTQ,N5VEC_QCLU,NNNQ,NOQ,NO_QCLU,NQ,
     &     NQCLUMAX,NQHOST,NQMAX,NQTB,NQ_I,NQ_L,NQ_R,QBAS,QBAS0_QCLU,
     &     QBAS_QCLU,QMGAM,QMGAM_POTFIL,QMPHI,QMPHI_POTFIL,QMTET,
     &     QMTET_POTFIL,Q_LMQ,RNNQ,VLMMAD_BACK,VLMMAD_HOST,IQBOT_CLU,
     &     NOMAX 
C
C*** End of declarations rewritten by SPAG
C
      DATA QMVEC/0.0D0,0.0D0,0.0D0/,SHIFT_CLU/0.0D0/
      DATA ALFDEG/0.0D0/,BETDEG/0.0D0/,GAMDEG/0.0D0/
      DATA IQBOT/1/,IQBOT_L/1/,IQBOT_R/1/,IQBOT_TB/0/,IQTOP_TB/0/
      DATA NQCLU/0/
      DATA NLMAD/0/,NLMMAD/0/,NLVMAD/0/,NLMVMAD/0/,NLQMAD/0/,NLMQMAD/0/
      DATA VMAD2D_A/0.0D0/,VMAD2D_B/0.0D0/
      DATA IQBOT_CHI/0/,IQTOP_CHI/0/,NQ_CHI/0/
      DATA OPTIMIZE_BASIS_GEOMETRY/.FALSE./
C
      ALLOCATABLE ABIMAD,DROTQ,AVMAD,DQBAS_QCLU,MDIRQ,MROTQ,QBAS
      ALLOCATABLE QMGAM,QMGAM_POTFIL,QMPHI,QMPHI_POTFIL
      ALLOCATABLE QMTET,QMTET_POTFIL,RNNQ,AVMAD_LI,AVMAD_LL
      ALLOCATABLE AVMAD_LR,AVMAD_RI,AVMAD_RL,AVMAD_RR
      ALLOCATABLE ICPA,IMQ,IM_QCLU,IQAT,IQCLU_ATCLU,IQ_QCLU
      ALLOCATABLE ITCLU_OQCLU,ITOQ,IT_OQCLU,KFP_LMQ,N5VEC_QCLU,NNNQ
      ALLOCATABLE MAGROT_Q,NOQ,NO_QCLU,QBAS_QCLU,FORCE_QCLU,QBAS0_QCLU
      ALLOCATABLE CMNTQ,Q_LMQ,KMAD1M_Q,VLMMAD_HOST,VLMMAD_BACK
C
      END
