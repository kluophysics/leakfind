C*==mod_scf.f    processed by SPAG 6.70Rc at 17:08 on 11 May 2014
      MODULE MOD_SCF
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables that control the                  *
C   *  SCF iterations                                                  *
C   *                                                                  *
C   ********************************************************************
C:
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ALAT_TAB(:),EFERMI_TAB(:),ETOT_TAB(:),MUEORB_TAB(:),
     &       MUESPN_TAB(:),RMSAVA,RMSAVB,RMSAVV,SCFMIX,SCFMIXOP,SCFSIM,
     &       SCFTOL,VOL_TAB(:),W_WEISS_ERROR
      INTEGER ITRSCF,LSCFALG,LSCFVXC,NSCFITER
      LOGICAL NFCOR,SCF_CHECK_SPLITSS,SCF_THETA_DEPENDENT_POT
      CHARACTER*10 POSVC,SCFALG,SCFSTATUS,SCFSTATUS_CLU,SCFSTATUS_HOST,
     &             SCFVXC
      SAVE ALAT_TAB,EFERMI_TAB,ETOT_TAB,LSCFALG,LSCFVXC,MUEORB_TAB,
     &     MUESPN_TAB,NSCFITER,RMSAVA,SCFMIX,SCFMIXOP,SCFTOL,VOL_TAB,
     &     W_WEISS_ERROR,SCF_THETA_DEPENDENT_POT
C
C*** End of declarations rewritten by SPAG
C
      DATA W_WEISS_ERROR/99.99D0/
      DATA RMSAVB/99.99D0/,RMSAVV/99.99D0/,ITRSCF/0/,SCFSIM/0D0/
      DATA SCFVXC/'VWN       '/,SCFALG/'BROYDEN2  '/
      DATA POSVC/'BM        '/,SCF_CHECK_SPLITSS/.TRUE./
      DATA NFCOR/.FALSE./
      DATA SCFSTATUS/'NOT SET   '/,SCFSTATUS_CLU/'NOT SET   '/
      DATA SCFSTATUS_HOST/'NOT SET   '/,SCF_THETA_DEPENDENT_POT/.FALSE./
C
      ALLOCATABLE ALAT_TAB,VOL_TAB,EFERMI_TAB,ETOT_TAB
      ALLOCATABLE MUEORB_TAB,MUESPN_TAB
C
      END
