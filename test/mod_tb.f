C*==mod_tb.f    processed by SPAG 6.55Rc at 16:27 on 10 Jul 2007
      MODULE MOD_TB
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables connected with                    *
C   *  TB calculations                                                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      LOGICAL CMPLX_KVEC,TBOPT_DECI_OUT,TBOPT_ONEBULK,TBOPT_SYMGK,
     &        VACFLAG(2)
      COMPLEX*16 FACTL(:,:)
      INTEGER ICHECK(:,:),IDECI,INVMOD,IREFQ(:),NKKRNR_TB,NKKR_RS,
     &        NKKR_TB,NKMSLAY,NPLAY,NREF,NSLAY_PER_PLAY,
     &     ISYM_IQ_R(:),IQ_L_IQ_R(:)
      REAL*8 RMTREF(:),VREF(:)
      CHARACTER*10 L_R_RELATION
      SAVE FACTL,ICHECK,IDECI,INVMOD,IREFQ,NKKRNR_TB,NKKR_RS,NKKR_TB,
     &     NKMSLAY,NPLAY,NREF,NSLAY_PER_PLAY,RMTREF,VREF,L_R_RELATION,
     &     ISYM_IQ_R,IQ_L_IQ_R
C
C*** End of declarations rewritten by SPAG
C
      DATA TBOPT_DECI_OUT/.FALSE./,TBOPT_SYMGK/.FALSE./,
     &     TBOPT_ONEBULK/.FALSE./,VACFLAG/.FALSE.,.FALSE./
     &     ,L_R_RELATION /'UNKNOWN   '/
C
      ALLOCATABLE RMTREF,VREF,IREFQ,FACTL,ICHECK,ISYM_IQ_R,IQ_L_IQ_R
C
      END
