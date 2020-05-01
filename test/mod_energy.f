C*==mod_energy.f    processed by SPAG 6.55Rc at 15:03 on 18 Jun 2010
      MODULE MOD_ENERGY
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tvariables used for   CPA   calculations    *
C   *                                                                  *
C   *  ALL data initialized in  <INIT_MOD_ENERGY>                      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      logical lepath_contn,lactive_contn
      integer nepol_contn
      real*8 emin_contn,emax_contn,ime_contn
      REAL*8 EFERMI,EILOW,EIMAG,EMAX,EMIN,IMEMAX,IMEMIN,EWORK
      real*8 around_fermi, temp_around_fermi
      COMPLEX*16 ETAB(:,:),WETAB(:,:)
c      complex*16 we_dmft_aux(:)
      COMPLEX*16 ECHEBY1(:),ECHEBY2(:)
      COMPLEX*16 ECHEBY1_BOT,ECHEBY1_TOP
      COMPLEX*16 ECHEBY2_BOT,ECHEBY2_TOP
      REAL*8 TCHEBY1(:,:),TCHEBY2(:,:)
      INTEGER NCHEBY1,NCHEBY2
      COMPLEX*16 PHASA(:,:,:),PHASK(:),PHAST(:,:,:)
      COMPLEX*16 ZGCHEBY1(:,:,:,:),ZFCHEBY1(:,:,:,:)
      COMPLEX*16 JGCHEBY1(:,:,:,:),JFCHEBY1(:,:,:,:)
      COMPLEX*16 ZGCHEBY2(:,:,:,:),ZFCHEBY2(:,:,:,:)
      COMPLEX*16 JGCHEBY2(:,:,:,:),JFCHEBY2(:,:,:,:)
      INTEGER IGRID(2),NE,NEMAX,NEPANEL,NEPATH,NETAB(2),
     &        NEPOL,nepoldummy,NEFD0,NEFD1,NEFD2,NEFD3,necontn
      LOGICAL SEARCHEF,SPLITSS

      SAVE EILOW,EMIN,ETAB,IGRID,IMEMAX,IMEMIN,NEMAX,NEPATH,SEARCHEF,
     &     WETAB,NEPOL,NEFD1,NEFD2,NEFD3,emin_contn,emax_contn,necontn,
     &     ime_contn,lepath_contn
      SAVE PHASA,PHASK,PHAST
      SAVE ECHEBY1,ECHEBY2,TCHEBY1,TCHEBY2,NCHEBY1,NCHEBY2
      SAVE ZGCHEBY1,ZFCHEBY1,JGCHEBY1,JFCHEBY1
      SAVE ZGCHEBY2,ZFCHEBY2,JGCHEBY2,JFCHEBY2
      SAVE ECHEBY1_BOT,ECHEBY1_TOP
      SAVE ECHEBY2_BOT,ECHEBY2_TOP
C     
C*** End of declarations rewritten by SPAG
C
      DATA EFERMI/ - 9999.D0/,EMAX/ - 9999.D0/,EIMAG/0.01D0/,NETAB/0,0/,
     &     SPLITSS/.FALSE./,NE/0/,NEPANEL/1/,EWORK/ - 9999.D0/,
     &     NEPOL/5/,NEFD1/3/,NEFD2/20/,NEFD3/10/,EMIN_CONTN/-0.2/,
     &     lepath_contn/.false./,lactive_contn/.false./,
     &     ime_contn/0.0d0/,around_fermi/0.0d0/,temp_around_fermi/0.0d0/
C
C---------------------------------------- variables depending only on NE
      ALLOCATABLE ETAB,WETAB

C-------------------------------- variables connected with Lloyd formula
      ALLOCATABLE PHASA,PHASK,PHAST

C---- variables connected with Chebychev interpolation of wave functions
      ALLOCATABLE ECHEBY1,ECHEBY2,TCHEBY1,TCHEBY2
      ALLOCATABLE ZGCHEBY1,ZFCHEBY1,JGCHEBY1,JFCHEBY1
      ALLOCATABLE ZGCHEBY2,ZFCHEBY2,JGCHEBY2,JFCHEBY2
C
      END
