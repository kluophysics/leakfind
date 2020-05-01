C*==mod_str.f    processed by SPAG 6.70Rc at 08:38 on 13 Sep 2012
      MODULE MOD_STR
C   ********************************************************************
C   *                                                                  *
C   *  module to store all data needed to set up the                   *
C   *  structure constants                                             *
C   *                                                                  *
C   *  NOTE: for USE_NEW_BBDD_VERSION=.TRUE.                           *
C   *        restrict site pairs to JQ <= IQ                           *
C   *        for QQMLRS, GGJLRS and IILERS dealt with in part CC       *
C   *                                                                  *
C   *  NOTE: for OPTIMIZE_LIMITS_JXXMAX =.TRUE.                        *
C   *        optimize the series limits  J13MAX  and  J22MAX           *
C   *        otherwise the defaults set in <STRINIT> are used          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL USE_NEW_BBDD_VERSION,OPTIMIZE_LIMITS_JXXMAX
      PARAMETER (USE_NEW_BBDD_VERSION=.TRUE.,
     &           OPTIMIZE_LIMITS_JXXMAX=.TRUE.)
      INTEGER J13MAX_LOWER_LIMIT,J13MAX_UPPER_LIMIT
      PARAMETER (J13MAX_LOWER_LIMIT=113,J13MAX_UPPER_LIMIT=300)
      INTEGER J22MAX_LOWER_LIMIT,J22MAX_UPPER_LIMIT
      PARAMETER (J22MAX_LOWER_LIMIT=30,J22MAX_UPPER_LIMIT=150)
C
C Local variables
C
      REAL*8 ALPHA0,CHP(:),ETA,GGJLRS(:,:,:),GMAXSQ,HP(:),M1PWL(:),
     &       PWEX2K1(:),PWEX2K2(:),PWEX2K3(:),QQPX(:),QQPY(:),QQPZ(:),
     &       RGNT(:),SHP(:)
      COMPLEX*16 CILMAT(:,:),D1TERM3(:),D300,DLLMMKE(:,:),EDU,
     &           EXPGNQ(:,:),IILERS(:,:,:),PWEXIK1(:),PWEXIK2(:),
     &           PWEXIK3(:),PWP(:),QQMLRS(:,:,:),SRREL(:,:,:),WK(:,:)
      INTEGER G1(:),G123MAX,G2(:),G3(:),IJQ(:,:),INDR(:,:),IRGNT(:),
     &        IRREL(:,:,:),J13MAX,J22MAX,LLARR,LLMAX,LRGNT12,LRGNT123,
     &        MMLLMAX,NGRL,NIJQ(:),NLLMMMAX,NQQP_STR,NQQP_STR_CC,
     &        NQQP_STR_RED,NRDL,NRGNT(:),NRREL(:,:),NSDL,R1(:),R123MAX,
     &        R2(:),R3(:),SMAX(:),STRMODE
      SAVE ALPHA0,CHP,CILMAT,D1TERM3,D300,DLLMMKE,EDU,ETA,EXPGNQ,G1,
     &     G123MAX,G2,G3,GGJLRS,GMAXSQ,HP,IILERS,IJQ,INDR,IRGNT,IRREL,
     &     J13MAX,J22MAX,LLARR,LLMAX,LRGNT12,LRGNT123,M1PWL,MMLLMAX,
     &     NGRL,NIJQ,NLLMMMAX,NQQP_STR,NQQP_STR_CC,NQQP_STR_RED,NRDL,
     &     NRGNT,NRREL,NSDL,PWEX2K1,PWEX2K2,PWEX2K3,PWEXIK1,PWEXIK2,
     &     PWEXIK3,PWP,QQMLRS,QQPX,QQPY,QQPZ,R1,R123MAX,R2,R3,RGNT,SHP,
     &     SMAX,SRREL,WK
C
C*** End of declarations rewritten by SPAG
C
      DATA STRMODE/1/
C
      ALLOCATABLE RGNT,IRGNT,NRGNT,IJQ,NIJQ,GGJLRS,PWP
      ALLOCATABLE SRREL,NRREL,IRREL
      ALLOCATABLE G1,G2,G3,R1,R2,R3,QQPZ,QQPY,QQPX,HP
      ALLOCATABLE EXPGNQ,INDR,SMAX,QQMLRS,CILMAT,M1PWL
      ALLOCATABLE D1TERM3,IILERS,DLLMMKE,CHP,SHP,WK
      ALLOCATABLE PWEX2K1,PWEX2K2,PWEX2K3,PWEXIK1,PWEXIK2,PWEXIK3
C
      END
