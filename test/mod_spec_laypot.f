C*==mod_transpho_laypot.f    processed by SPAG 6.70Rc at 09:25 on  3 Feb 2012
      MODULE MOD_TRANSPHO_LAYPOT
C   ********************************************************************
C   *                                                                  *
C   *  module to store conversion tables between SPR and PHOTO         *
C   *  for potentials                                                  *
C   *                                                                  *
C   *  ALL data initialized in  <<transpho_potfull_spkkr>>             *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IBLOCHTMP,IQ_SPR_LAYAT(:,:)
      SAVE IBLOCHTMP,IQ_SPR_LAYAT
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATABLE IQ_SPR_LAYAT
C
      END
