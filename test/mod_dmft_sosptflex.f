C*==mod_dmft_sosptflex.f    processed by SPAG 6.70Rc at 15:43 on  2 Feb 2012
      MODULE MOD_DMFT_SOSPTFLEX
C**********************************************************
C
C modules for SPTF-FLEX+SO
C
C***********************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI
      PARAMETER (CI=(0.D0,1.D0))
      REAL*8 PI
      PARAMETER (PI=3.1415926535898D0)
C
C Local variables
C
      REAL*8 DOMEGA,DTIME,G_MU(:,:),MU_1,OMEGA(:),TIME(:)
      COMPLEX*16 ENDSIG,ENDSIGMA(:,:),PPH(:,:,:,:,:),PPP(:,:,:,:,:),
     &           TLAST(:,:,:,:),UC(:,:,:,:),UW0(:,:,:,:)
      INTEGER G_F,G_FIL(:,:),IPRT,MINOM(:),NLM,NLMS,NN,NNNN,NOM,NS
      CHARACTER*1 SIGNATURE
      SAVE DOMEGA,DTIME,ENDSIG,ENDSIGMA,G_F,G_FIL,G_MU,IPRT,MINOM,MU_1,
     &     NLM,NLMS,NN,NNNN,NOM,NS,OMEGA,PPH,PPP,SIGNATURE,TIME,TLAST,
     &     UC,UW0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ENDSIGMA,G_FIL,G_MU,MINOM,OMEGA,PPH,PPP,TIME,TLAST
      ALLOCATABLE UC,UW0
C
C*** End of declarations rewritten by SPAG
C
      END
C
C
