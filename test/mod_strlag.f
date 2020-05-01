C*==mod_strlag.f    processed by SPAG 6.55Rc at 08:58 on  9 Jul 2007
      MODULE MOD_STRLAG
C   ********************************************************************
C   *                                                                  *
C   *  module to store all needed to interpolate the                   *
C   *  structure constants                                             *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      COMPLEX*16 DLMLAG(:,:,:,:),ESTRLAG(:),EXPGNQ_RED(:,:),WSTRLAG(:)
      INTEGER NGFEP,NSTRLAG
      SAVE DLMLAG,ESTRLAG,EXPGNQ_RED,NGFEP,NSTRLAG,WSTRLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE EXPGNQ_RED,ESTRLAG,WSTRLAG,DLMLAG
C
      END
