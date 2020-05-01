C*==mod_kspace.f    processed by SPAG 6.70Rc at 10:38 on 26 Feb 2013
      MODULE MOD_KSPACE
C   ********************************************************************
C   *                                                                  *
C   *  module to store all kspace information                          *
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
      INTEGER NGFEPMAX
      PARAMETER (NGFEPMAX=200)
C
C Local variables
C
      LOGICAL BNDSPLIT_BZ(:),CMPLX_KVEC,USE_TZ(:)
      REAL*8 BRA_KZ(:),GBAD(:,:),GFEP(:,:),IM_KTAB(:,:),KTAB(:,:),
     &       KTET(:,:),KVEC_KZ(:,:),WKSUM,WKTAB(:),W_KZ(:),W_TZ(:)
      INTEGER IBND_CTZ(:,:),IBZINT,IKCTET(:,:),IKS_BZ(:),IK_B1Z(:),
     &        IK_B2Z(:),IK_CTZ(:,:),ITBZ(:),JTBZ(:),NELMT,NGBAD,
     &        NGBADMAX,NGFEP,NKPTS0,NKTAB,NKTABINP,NKTABMAX,NKTET,
     &        NPTMAX,NPTMIN,NTAUUV(:),NTETS,NZOOM,QTBZ(:),UTAUUV(:),
     &        VTAUUV(:)
      COMPLEX*16 WTAUUV(:)
      SAVE BNDSPLIT_BZ,BRA_KZ,GBAD,GFEP,IBND_CTZ,IKCTET,IKS_BZ,IK_B1Z,
     &     IK_B2Z,IK_CTZ,IM_KTAB,ITBZ,JTBZ,KTAB,KTET,KVEC_KZ,NELMT,
     &     NGBAD,NGFEP,NKPTS0,NKTABINP,NTAUUV,QTBZ,USE_TZ,UTAUUV,VTAUUV,
     &     WKTAB,WTAUUV,W_KZ,W_TZ
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE KTAB,WKTAB,KTET,IKCTET,IM_KTAB
C
C-------------------------------------- variables for tetrahedron method
C
      ALLOCATABLE ITBZ,JTBZ,QTBZ
      ALLOCATABLE NTAUUV,WTAUUV,UTAUUV,VTAUUV
C
C---------------------------- variables depending on free electron poles
C
      ALLOCATABLE GBAD,GFEP
C
C------------------------------ variables used for zoom-in calclulations
C
      ALLOCATABLE IK_CTZ,KVEC_KZ,W_KZ,USE_TZ,W_TZ,IBND_CTZ,BRA_KZ
      ALLOCATABLE BNDSPLIT_BZ,IK_B1Z,IK_B2Z,IKS_BZ
C
      DATA IBZINT/2/,NKTAB/300/,NTETS/1/,NKTET/1/
      DATA NGBADMAX/1/,NPTMIN/300/,NPTMAX/500/,NKTABMAX/10000/
      DATA CMPLX_KVEC/.FALSE./,NZOOM/3/,WKSUM/ - 9999999D0/
C
      END
