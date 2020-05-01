C*==clu_chi_drive.f    processed by SPAG 6.70Rc at 08:57 on  8 Mar 2017
      SUBROUTINE CLU_CHI_DRIVE
C   ********************************************************************
C   *                                                                  *
C   *  driving routine for linear response calculations                *
C   *  according to TASK it calls                                      *
C   *                                                                  *
C   *  - CHI       magnetic susceptibility                             *
C   *  - CHIXAS    field induced XMXD                                  *
C   *  - GILBERT   Gilbert damping parameter                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:CHIZ
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_SITES,ONLY:NQ,IQBOT_CHI,IQTOP_CHI,NQ_CHI
      USE MOD_TYPES,ONLY:NTMAX,DOBS_LTEX,DOBS_TEX_GLO
      USE MOD_CALCMODE,ONLY:KKRMODE,TASK
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_ANGMOM,ONLY:NKMMAX,NOBSMAX,NLMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLU_CHI_DRIVE')
C
C Local variables
C
      INTEGER IA_ERR,M
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
C              regime of sites treated for response functions
C ======================================================================
C
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
         IQBOT_CHI = 1
         IQTOP_CHI = NQ
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
      END IF
C
      NQ_CHI = IQTOP_CHI - IQBOT_CHI + 1
C
C=======================================================================
C
      ALLOCATE (DOBS_LTEX(0:3,NOBSMAX,NLMAX,NTMAX,NEMAX))
      ALLOCATE (DOBS_TEX_GLO(0:3,NOBSMAX,NTMAX,NEMAX))
C
C-----------------------------------------------------------------------
C
      IF ( FULLPOT .AND. TASK.NE.'CHI       ' .AND. TASK(1:5)
     &     .NE.'SIGMA' .AND. TASK(1:6).NE.'MAGNET' )
     &     CALL STOP_MESSAGE(ROUTINE,
     &     'the chosen TASK does not YET work for FULL POTENTIAL mode')
C
C ======================================================================
C                       response functions
C ======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      WRITE (6,99001) TASK
C
      M = NKMMAX
C
      ALLOCATE (CHIZ(M*M*NQ_CHI,M*M*NQ_CHI,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CJI')
C
C ======================================================================
C                         perform requested TASK
C ======================================================================
C
C ======================================================================
C                        GILBERT DAMPING
C ======================================================================
C
      IF ( TASK.EQ.'GILBERT   ' ) THEN
C
         CALL GIL
C
C ======================================================================
      ELSE
C
         WRITE (6,*) '########################## TASK = ',TASK
         WRITE (6,*) '##  allowed: TASK = GILBERT'
         CALL STOP_MESSAGE(ROUTINE,'TASK NOT allowed')
      END IF
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
99001 FORMAT (/,10X,'Linear response calculation for  TASK = ',A,//)
C
      END
