C*==init_mrot_frame.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE INIT_MROT_FRAME(CHANGE_FRAME,MROT_FRAME)
C   ********************************************************************
C   *                                                                  *
C   *   specify a special frame of reference system                    *
C   *   by the angles THETA and PHI giving the orientation of the      *
C   *   special z-axis w.r.t. the common GLOBAL frame of reference     *
C   *   and additional angle GAMMA may be given to orient the frame    *
C   *                                                                  *
C   *   NOTE:      ->z_global = MROT_frame * ->z_special               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FOUND_SECTION,FOUND_REAL
      IMPLICIT NONE
C*--INIT_MROT_FRAME16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CHANGE_FRAME
      REAL*8 MROT_FRAME(3,3)
C
C Local variables
C
      REAL*8 FRAMEGAM,FRAMEPHI,FRAMETET
C
C*** End of declarations rewritten by SPAG
C
      DATA FRAMETET/0D0/,FRAMEPHI/0D0/,FRAMEGAM/0D0/
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( .NOT.FOUND_SECTION ) RETURN
C
C=======================================================================
C
      CALL SECTION_SET_REAL('FRAMETET',FRAMETET,9999D0,0)
      CHANGE_FRAME = FOUND_REAL
      CALL SECTION_SET_REAL('FRAMEPHI',FRAMEPHI,9999D0,0)
      CHANGE_FRAME = CHANGE_FRAME .OR. FOUND_REAL
      CALL SECTION_SET_REAL('FRAMEGAM',FRAMEGAM,9999D0,0)
      CHANGE_FRAME = CHANGE_FRAME .OR. FOUND_REAL
C
      IF ( .NOT.CHANGE_FRAME ) RETURN
C
C=======================================================================
C
      CALL GETMROT(-FRAMEGAM,-FRAMETET,-FRAMEPHI,MROT_FRAME)
C
      WRITE (6,99001) FRAMEGAM,FRAMETET,FRAMEPHI,MROT_FRAME
C
      RETURN
C
99001 FORMAT (//,1X,79('*'),/,34X,'<INIT_MROT_FRAME>',/,1X,79('*'),//,
     &        10X,'special frame of reference specified by',//,10X,
     &        'GAM =',F8.3,3X,'TET =',F8.3,3X,'PHI =',F8.3,3X,//,10X,
     &        'rotation matrix:',//,
     &        3(30X,'(',F8.3,',',F8.3,',',F8.3,'  )',/),/)
      END
