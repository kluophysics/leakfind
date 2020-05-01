C*==mod_dmft_ldau.f    processed by SPAG 6.70Rc at 14:35 on  2 Feb 2012
      MODULE MOD_DMFT_LDAU
C
C   ********************************************************************
C   *                                                                  *
C   *  module to store all informations concerning LDA+U or DMFT       *
C   *                                                                  *
C   *  ALL data initialized in  << DMFT_INIT  >>                       *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER*4 DMFTDBLC,UMODE
      LOGICAL DMFTPOT,DMFTSCF,SYMMETRISE_OCC,DMFT_FIX_DYN_SE
      LOGICAL DMFT_FIX_BAS, ZEROSIG_CONTN
      COMPLEX*16 DMFTSIG(:,:,:),DMFTSIGMA(:,:,:,:),EREFDMFT(:),
     &           EREFLDAU(:),SEBNST(:,:,:),SEBT(:,:),SEVNST(:,:,:),
     &           SEVT(:,:)
      CHARACTER*10 DMFTSOLVER
      REAL*8 DMFTTEMP,ELDAU(:),JEFF(:),UEFF(:),SIGTOL,WLDAU(:,:,:)
      INTEGER KSELF(:),IBASIS,IEREF
      SAVE DMFTSCF,DMFTSIG,DMFTSIGMA,DMFTSOLVER,DMFTTEMP,EREFDMFT,JEFF,
     &     KSELF,SEBNST,SEBT,SEVNST,SEVT,UEFF,IBASIS,UMODE,
     &     DMFTDBLC,SIGTOL,IEREF
      
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE DMFTSIG,EREFDMFT,EREFLDAU,WLDAU,KSELF,ELDAU
      ALLOCATABLE DMFTSIGMA,SEBT,SEVT,SEBNST,SEVNST,UEFF,JEFF
C
      DATA SYMMETRISE_OCC/.FALSE./
      DATA DMFT_FIX_DYN_SE/.TRUE./
      DATA DMFT_FIX_BAS/.FALSE./
C
      END
