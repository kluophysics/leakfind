C*==init_mod_tb.f    processed by SPAG 6.55Rc at 13:17 on 19 Nov 2008
      SUBROUTINE INIT_MOD_TB(NM,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TB,ONLY:NREF,RMTREF,VREF,IREFQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NM,NQMAX
C
C*** End of declarations rewritten by SPAG
C
      NREF = NM
C
      ALLOCATE (RMTREF(NREF),VREF(NREF),IREFQ(NQMAX))
C
      RMTREF(1:NREF) = 0D0
      VREF(1:NREF) = 0D0
C
      END
