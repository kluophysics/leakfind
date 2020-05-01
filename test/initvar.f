C*==initvar.f    processed by SPAG 6.70Rc at 08:48 on  8 Mar 2017
      SUBROUTINE INITVAR(NQ,TASK,POTFMTINP,POTFMTOUT,ITEST,TXTTEST,IREL,
     &                   ICORE,IHYPER,ISMQHFI,IXRAY,RDTAU,RDTAUMQ,RDDOS,
     &                   WRLOG,WRMAT,WRTAU,WRTAUMQ,USETAUNN,USELPLS1,
     &                   RUNELOOP,INITELOOP,SPLITSS,NOWRDOS,IBZINT,
     &                   KMROT,KMROT_COMMON,QMVEC,MDIRQ,ALFDEG,BETDEG,
     &                   GAMDEG,LLOYD,NONMAG,PLOTPRS,PLOT2DPR,BREITINT5,
     &                   NQMAX,ITESTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  read standard input file  IFILINP   and read all variables      *
C   *  that might influence the symmetry and the array sizes           *
C   *                                                                  *
C   *  this applies in particular to   TASK                            *
C   *  accordingly the variables associated with TASK are set as well  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY,NO_SYMMETRY_CLU,SYMMETRIZE_MSS,
     &    SYMMETRIZE_RHO,NO_SYMMETRY_LINRESP
      USE MOD_CALCMODE,ONLY:PROGNAME,SYMCHECK,IREL_HOST,SIGMA_PROJECT,
     &    SIGMA_LAYER,THERMAL_VIBRA_FLUCT,TUTORIAL,USE_FPCHR_ZJOVZZ,
     &    POTFIX,USE_CONST_POTENTIAL,ROTATE_LATTICE,BREAKPOINT,CALC_EFG,
     &    CALC_FORCE,SETPOT
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      USE MOD_FILES,ONLY:IPRINT,DEBUG,WRNOS,RDUMMY,WR_TRACK_INFO,
     &    FOUND_SECTION,FOUND_REAL,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_SITES,ONLY:NQCLU,QMPHI,QMTET,QMGAM
      USE MOD_FILES,ONLY:PLOT_JALF,CALC_JALF,CALC_JORB,PLOT_JORB
      USE MOD_CPA,ONLY:NCPA,USENLCPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INITVAR')
C
C Dummy arguments
C
      REAL*8 ALFDEG,BETDEG,GAMDEG
      LOGICAL BREITINT5,INITELOOP,LLOYD,NONMAG,NOWRDOS,RDDOS,RDTAU,
     &        RDTAUMQ,RUNELOOP,SPLITSS,USELPLS1,USETAUNN,WRLOG,WRMAT,
     &        WRTAU,WRTAUMQ
      INTEGER IBZINT,ICORE,IHYPER,IREL,ISMQHFI,ITEST,ITESTMAX,IXRAY,
     &        KMROT,KMROT_COMMON,NQ,NQMAX,POTFMTINP,POTFMTOUT
      CHARACTER*10 TASK
      REAL*8 MDIRQ(3,NQMAX),QMVEC(3)
      LOGICAL PLOT2DPR(0:3),PLOTPRS(0:3)
      CHARACTER*40 TXTTEST(0:ITESTMAX)
C
C Local variables
C
      LOGICAL ANGLE_INTEGRATED,FOUND,FOUND3,FOUND4
      REAL*8 EPS_MACH,MDIR(3),PI_CALC
      INTEGER IQ
      CHARACTER*3 MODE
      CHARACTER*10 SELFENERGY,STR10
C
C*** End of declarations rewritten by SPAG
C
      DATA MDIR/0.0D0,0.0D0,1.0D0/,MODE/'   '/
C
      CALL TRACK_INFO(ROUTINE)
C
C-------------------------------------------------------------------- PI
C                             check PI
C-------------------------------------------------------------------- PI
C
      EPS_MACH = EPSILON(PI)
      PI_CALC = ATAN(1D0)*4D0
C
      IF ( ABS(PI-PI_CALC).GT.2*EPS_MACH ) THEN
         WRITE (6,99007) PI,PI_CALC,EPS_MACH
         CALL STOP_MESSAGE(ROUTINE,'PI inconsistent')
      END IF
C
C
C-----------------------------------------------------------------------
C
C=======================================================================
C                     read  CONTROL  parameters
C=======================================================================
C
      WRLOG = .FALSE.
      DEBUG = .FALSE.
      WRMAT = .FALSE.
      WRTAU = .FALSE.
      RDTAU = .FALSE.
      USELPLS1 = .FALSE.
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('WRNOS',WRNOS)
         IF ( WRNOS ) LLOYD = .TRUE.
         CALL SECTION_FIND_KEYWORD('WRLOG',WRLOG)
         CALL SECTION_FIND_KEYWORD('DEBUG',DEBUG)
         CALL SECTION_FIND_KEYWORD('TRACK',WR_TRACK_INFO)
         IF ( DEBUG ) WR_TRACK_INFO = .TRUE.
         CALL SECTION_FIND_KEYWORD('WRMAT',WRMAT)
         CALL SECTION_FIND_KEYWORD('WRTAU',WRTAU)
         CALL SECTION_FIND_KEYWORD('RDTAU',RDTAU)
         CALL SECTION_FIND_KEYWORD('WRTAUMQ',WRTAUMQ)
         CALL SECTION_FIND_KEYWORD('RDTAUMQ',RDTAUMQ)
         CALL SECTION_FIND_KEYWORD('PLOTPOT',PLOTPRS(1))
         CALL SECTION_FIND_KEYWORD('PLOTRHO',PLOTPRS(2))
         CALL SECTION_FIND_KEYWORD('PLOTSFN',PLOTPRS(3))
         PLOTPRS(0) = PLOTPRS(1) .OR. PLOTPRS(2) .OR. PLOTPRS(3)
         CALL SECTION_FIND_KEYWORD('PLOTJALF',PLOT_JALF)
         CALL SECTION_FIND_KEYWORD('PLOTJORB',PLOT_JORB)
         IF ( PROGNAME(4:6).EQ.'SCF' ) THEN
            CALL SECTION_FIND_KEYWORD('EFG',CALC_EFG)
            CALL SECTION_FIND_KEYWORD('FORCE',CALC_FORCE)
         END IF
      END IF
C
      IF ( PROGNAME(4:6).EQ.'SCF' .AND. IPRINT.GE.1 ) WRLOG = .TRUE.
      WRLOG = .TRUE.
C
C=======================================================================
C                     read  TASK  parameters
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION .OR. (PLOTPRS(0) .OR. PLOT2DPR(0)) ) THEN
C
         CALL SECTION_FIND_KEYWORD('SCF',FOUND)
         IF ( FOUND ) TASK = 'SCF       '
C
         CALL SECTION_FIND_KEYWORD('WFPLOT',FOUND)
         IF ( FOUND ) TASK = 'WFPLOT    '
C
         CALL SECTION_FIND_KEYWORD('MEPLOT',FOUND)
         IF ( FOUND ) TASK = 'MEPLOT    '
C
         CALL SECTION_FIND_KEYWORD('MECHECK',FOUND)
         IF ( FOUND ) TASK = 'MECHECK   '
C
         CALL SECTION_FIND_KEYWORD('SOCPAR',FOUND)
         IF ( FOUND ) TASK = 'SOCPAR    '
C
         CALL SECTION_FIND_KEYWORD('PSHIFT',FOUND)
         IF ( FOUND ) TASK = 'PSHIFT    '
C
         CALL SECTION_FIND_KEYWORD('FSCAT',FOUND)
         IF ( FOUND ) TASK = 'FSCAT     '
C
         CALL SECTION_FIND_KEYWORD('PLOTPRS',FOUND)
         IF ( FOUND ) THEN
            TASK = 'PLOTPRS   '
            PLOTPRS(0:3) = .TRUE.
         ELSE
            CALL SECTION_FIND_KEYWORD('PLOTPOT',PLOTPRS(1))
            CALL SECTION_FIND_KEYWORD('PLOTRHO',PLOTPRS(2))
            CALL SECTION_FIND_KEYWORD('PLOTSFN',PLOTPRS(3))
            PLOTPRS(0) = PLOTPRS(1) .OR. PLOTPRS(2) .OR. PLOTPRS(3)
            IF ( PLOTPRS(0) ) TASK = 'PLOTPRS   '
         END IF
C
         CALL SECTION_FIND_KEYWORD('PLOT2DPOT',PLOT2DPR(1))
         CALL SECTION_FIND_KEYWORD('PLOT2DRHO',PLOT2DPR(2))
         PLOT2DPR(0) = PLOT2DPR(1) .OR. PLOT2DPR(2)
         CALL SECTION_FIND_KEYWORD('PLOT2DPR',FOUND)
         IF ( FOUND ) PLOT2DPR(0:2) = .TRUE.
         IF ( PLOT_JALF ) TASK = 'PLOT2DPR  '
         IF ( PLOT_JORB ) TASK = 'PLOT2DPR  '
C
         IF ( TASK.EQ.'WFPLOT    ' .OR. TASK.EQ.'MECHECK   ' .OR. 
     &        TASK.EQ.'MEPLOT    ' .OR. TASK.EQ.'SOCPAR    ' .OR. 
     &        TASK.EQ.'PSHIFT    ' .OR. TASK.EQ.'FSCAT     ' .OR. 
     &        TASK.EQ.'PLOTPRS   ' .OR. TASK.EQ.'PLOT2DPR  ' ) THEN
            RDTAU = .FALSE.
            WRTAU = .FALSE.
            WRTAUMQ = .FALSE.
            WRTAU = .FALSE.
            NOWRDOS = .TRUE.
            IHYPER = 0
            ICORE = 0
            RUNELOOP = .FALSE.
         END IF
C
         ANGLE_INTEGRATED = .TRUE.
         CALL SECTION_FIND_KEYWORD('VBXPS',FOUND)
         IF ( .NOT.FOUND ) THEN
            CALL SECTION_FIND_KEYWORD('VBPES',FOUND)
            IF ( .NOT.FOUND ) THEN
               CALL SECTION_FIND_KEYWORD('ARPES',FOUND)
               ANGLE_INTEGRATED = .NOT.FOUND
               IF ( .NOT.FOUND )
     &              CALL SECTION_FIND_KEYWORD('VBUPS',FOUND)
            END IF
         ELSE
            CALL SECTION_FIND_KEYWORD('ARPES',FOUND)
            ANGLE_INTEGRATED = .NOT.FOUND
         END IF
         IF ( FOUND ) THEN
            IF ( ANGLE_INTEGRATED ) THEN
               TASK = 'VBPES     '
               RDTAU = .TRUE.
            ELSE
               TASK = 'ARPES     '
               RDTAU = .FALSE.
            END IF
            WRTAU = .FALSE.
            IHYPER = 0
            ICORE = 0
            IF ( PROGNAME.EQ.'KKRGEN    ' ) THEN
               USELPLS1 = .TRUE.
               CALL SECTION_FIND_KEYWORD('USETAUNN',USETAUNN)
               IF ( .NOT.USETAUNN )
     &              CALL SECTION_FIND_KEYWORD('ARPES',USETAUNN)
            END IF
         END IF
C
         CALL SECTION_FIND_KEYWORD('BIS',FOUND)
         IF ( FOUND ) THEN
            TASK = 'BIS       '
            RDTAU = .TRUE.
            IHYPER = 0
            ICORE = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('CLXPS',FOUND)
         IF ( FOUND ) THEN
            NOWRDOS = .TRUE.
            TASK = 'CLXPS     '
            RDTAU = .TRUE.
            USELPLS1 = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('APS',FOUND)
         IF ( FOUND ) THEN
            TASK = 'APS       '
            RUNELOOP = .FALSE.
            USELPLS1 = .TRUE.
            RDDOS = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('AES',FOUND)
         IF ( FOUND ) THEN
            TASK = 'AES       '
            RUNELOOP = .FALSE.
            RDTAU = .TRUE.
            WRTAU = .FALSE.
            USELPLS1 = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('NRAES',FOUND)
         IF ( FOUND ) THEN
            TASK = 'NRAES     '
            RUNELOOP = .FALSE.
            USELPLS1 = .TRUE.
            RDDOS = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('RELAX',FOUND)
         IF ( FOUND ) THEN
            TASK = 'RELAX     '
            RUNELOOP = .FALSE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('XAS',FOUND)
         IF ( FOUND ) TASK = 'XAS       '
C
         CALL SECTION_FIND_KEYWORD('XES',FOUND)
         IF ( FOUND ) TASK = 'XES       '
C
         CALL SECTION_FIND_KEYWORD('XMO',FOUND)
         IF ( FOUND ) TASK = 'XMO       '
C
         CALL SECTION_FIND_KEYWORD('XRS',FOUND)
         IF ( FOUND ) TASK = 'XRS       '
C
         IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &        TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
            NOWRDOS = .TRUE.
            RDTAU = .TRUE.
            IXRAY = 1
         END IF
C
         CALL SECTION_FIND_KEYWORD('T1',FOUND)
         IF ( FOUND ) THEN
            TASK = 'T1        '
            RUNELOOP = .FALSE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('RHO',FOUND)
         IF ( FOUND ) THEN
            TASK = 'RHO       '
            ICORE = 0
            NOWRDOS = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('OPTICS',FOUND)
         IF ( .NOT.FOUND ) CALL SECTION_FIND_KEYWORD('OPT',FOUND)
         IF ( FOUND ) THEN
            IF ( PROGNAME.NE.'KKRCHI    ' ) CALL STOP_MESSAGE(ROUTINE,
     &           'TASK OPTICS has to be run by KKRCHI')
            TASK = 'OPTICS    '
            ICORE = 0
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('SIGMA',FOUND)
         IF ( FOUND .OR. TASK.EQ.'OPTICS    ' ) THEN
            IF ( PROGNAME.NE.'KKRCHI    ' ) CALL STOP_MESSAGE(ROUTINE,
     &           'TASK SIGMA has to be run by KKRCHI')
            TASK = 'SIGMA     '
            ICORE = 0
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
            CALL SECTION_FIND_KEYWORD('PROJECT',SIGMA_PROJECT)
            CALL SECTION_FIND_KEYWORD('LAYER',SIGMA_LAYER)
         END IF
C
         CALL SECTION_FIND_KEYWORD('NEGF',FOUND)
         IF ( FOUND ) THEN
            TASK = 'NEGF'
            ICORE = 0
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
            IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) NO_SYMMETRY = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('GILBERT',FOUND)
         IF ( FOUND ) THEN
            TASK = 'GILBERT  '
            ICORE = 0
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('STONER',FOUND)
         IF ( FOUND ) THEN
            TASK = 'STONER   '
            ICORE = 0
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('PHONONS',FOUND)
         IF ( FOUND ) THEN
            IF ( PROGNAME.NE.'KKRCHI    ' ) CALL STOP_MESSAGE(ROUTINE,
     &           'TASK PHONONS has to be run by KKRCHI')
            TASK = 'PHONONS  '
            ICORE = 0
            NOWRDOS = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('PUMP',FOUND)
         IF ( FOUND ) THEN
            IF ( PROGNAME.NE.'KKRCHI    ' ) CALL STOP_MESSAGE(ROUTINE,
     &           'TASK PUMP    has to be run by KKRCHI')
            TASK = 'PUMP     '
            ICORE = 0
            NOWRDOS = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('OPM',FOUND)
         IF ( FOUND ) TASK = 'OPM       '
C
         CALL SECTION_FIND_KEYWORD('CHI',FOUND)
         IF ( FOUND ) TASK = 'CHI       '
C
         CALL SECTION_FIND_KEYWORD('CHIDYN',FOUND)
         IF ( FOUND ) TASK = 'CHIDYN    '
C
C=======================================================================
C       deal with  SYMMETRY  in case of linear response via  KKRCHI
C=======================================================================
         IF ( PROGNAME.EQ.'KKRCHI    ' ) THEN
C
            IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) NO_SYMMETRY = .TRUE.
C
            CALL SECTION_FIND_KEYWORD('NOSYMSIG',NO_SYMMETRY_LINRESP)
C
            IF ( .NOT.NO_SYMMETRY_LINRESP )
     &            CALL SECTION_FIND_KEYWORD('NOSYMCHI',
     &           NO_SYMMETRY_LINRESP)
C
CDK------------- suppress globally symmetry if lin resp doesn't use it
CDK------------- TODO: use two different kmesh arrays later to allow for
CDK-------------      for more efficiency and freedom
            IF ( NO_SYMMETRY_LINRESP ) NO_SYMMETRY = .TRUE.
C
         END IF
C=======================================================================
C
         CALL SECTION_FIND_KEYWORD('MAGNET',FOUND)
         IF ( FOUND ) TASK = 'MAGNET    '
C
         CALL SECTION_FIND_KEYWORD('CHILANDAU',FOUND)
         IF ( FOUND ) TASK = 'CHILANDAU '
C
         CALL SECTION_FIND_KEYWORD('CHIXAS',FOUND)
         IF ( FOUND ) TASK = 'CHIXAS    '
C
         CALL SECTION_FIND_KEYWORD('EKREL',FOUND)
         IF ( FOUND ) THEN
            TASK = 'EKREL     '
            NOWRDOS = .TRUE.
            RUNELOOP = .FALSE.
            ICORE = 0
            IHYPER = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('BSF',FOUND)
         IF ( FOUND ) THEN
            TASK = 'BSF       '
            RDTAUMQ = .TRUE.
            NOWRDOS = .TRUE.
            RUNELOOP = .FALSE.
            ICORE = 0
            IHYPER = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('COMPTON',FOUND)
         IF ( FOUND ) THEN
            TASK = 'COMPTON   '
            RDTAUMQ = .TRUE.
            NOWRDOS = .TRUE.
            RUNELOOP = .FALSE.
            ICORE = 0
            IHYPER = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('NUCSSCPL',FOUND)
         IF ( FOUND ) THEN
            TASK = 'NUCSSCPL  '
            RUNELOOP = .TRUE.
            ICORE = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('FORCETEN',FOUND)
         IF ( FOUND ) THEN
            TASK = 'FORCETEN  '
            RUNELOOP = .TRUE.
            ICORE = 0
         END IF
C
         CALL SECTION_FIND_KEYWORD('BONDING',FOUND)
         IF ( FOUND ) THEN
            TASK = 'BONDING   '
            RDTAU = .FALSE.
            WRTAU = .FALSE.
            IHYPER = 1
            ICORE = 1
         END IF
C
         CALL SECTION_FIND_KEYWORD('TORQUE',FOUND)
         IF ( FOUND ) THEN
            TASK = 'TORQUE    '
            NOWRDOS = .TRUE.
            RDTAUMQ = .TRUE.
            WRTAUMQ = .FALSE.
            INITELOOP = .TRUE.
         END IF
C
         CALL SECTION_FIND_KEYWORD('TUTORIAL',FOUND)
         IF ( FOUND ) THEN
            TASK = 'TUTORIAL  '
            RUNELOOP = .TRUE.
            WRTAU = .FALSE.
            WRTAUMQ = .FALSE.
            TUTORIAL = .TRUE.
         END IF
C
C-------------------- allow tasks that need ELOOP to be done in parallel
C
         CALL SECTION_FIND_KEYWORD('FMAG',FOUND3)
         CALL SECTION_FIND_KEYWORD('JXC',FOUND4)
         CALL SECTION_FIND_KEYWORD('JALF',CALC_JALF)
         CALL SECTION_FIND_KEYWORD('JORB',CALC_JORB)
C
         IF ( FOUND3 .OR. FOUND4 .OR. CALC_JALF .OR. CALC_JORB ) THEN
            TASK = 'TL00000000'
            IF ( FOUND3 ) TASK = 'TL10000000'
            IF ( FOUND4 ) TASK = TASK(1:3)//'1'//TASK(5:10)
            NOWRDOS = .TRUE.
            RDTAUMQ = .FALSE.
            RDTAU = .FALSE.
            WRTAUMQ = .FALSE.
            WRTAU = .FALSE.
            ICORE = 0
            IHYPER = 0
         END IF
C
         CALL SECTION_SET_STRING('SELFENERGY',SELFENERGY,'9999',0)
C
C        IF( TASK(1:4).EQ.'NONE' .OR. TASK(1:2).EQ.'T1'    .OR.
C    &       TASK(1:3).EQ.'RHO'  .OR. TASK(1:5).EQ.'EKREL' .OR.
C    &       TASK(1:3).EQ.'BSF'  ) THEN
C           SELFENERGY = 'NONE'
C        END IF
C
      END IF
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('NOHFF',FOUND)
         IF ( FOUND ) IHYPER = 0
         CALL SECTION_FIND_KEYWORD('NOCORE',FOUND)
         IF ( FOUND ) ICORE = 0
         CALL SECTION_FIND_KEYWORD('FSOHFF',FOUND)
         IF ( FOUND ) ISMQHFI = 1
C        CALL SECTION_FIND_KEYWORD('EFG',FOUND)
C        IF ( FOUND ) ISMQHFI = 1
      END IF
C
      IF ( ISMQHFI.NE.0 ) SPLITSS = .TRUE.
C
      CALL INPUT_FIND_SECTION('ENERGY',0)
      IF ( FOUND_SECTION ) THEN
         IF ( .NOT.SPLITSS ) CALL SECTION_FIND_KEYWORD('SPLITSS',
     &        SPLITSS)
      END IF
C
C ---------------------------- reset  WRTAU   RDTAU   ICORE if necessary
C
      IF ( TASK(1:7).EQ.'MECHECK' .OR. TASK(1:5).EQ.'VBPES' .OR. 
     &     TASK(1:5).EQ.'ARPES' ) THEN
         IHYPER = 0
         ISMQHFI = 0
         ICORE = 0
         IXRAY = 0
      END IF
C
      IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &     TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' .OR. 
     &     IXRAY.NE.0 ) THEN
         ICORE = 0
         RDTAU = .TRUE.
         WRTAU = .FALSE.
      END IF
C
      IF ( (MODE.EQ.'IMP') .OR. (MODE.EQ.'imp') ) THEN
         RDTAU = .TRUE.
         WRTAU = .FALSE.
      END IF
C
C---------------------------------------------------------------- IHYPER
C
C   DEFAULT: no calculation of HFF with SCF calculation
C            force calculation by setting  HFF
C
      IF ( PROGNAME(4:6).EQ.'SCF' ) IHYPER = 0
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('HFF',FOUND)
         IHYPER = 1
      ELSE
         CALL INPUT_FIND_SECTION('SCF',0)
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('HFF',FOUND)
            IHYPER = 1
         END IF
      END IF
C
C----------------------------------------------------------------- ITEST
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_FIND_KEYWORD('NOSYM',FOUND)
         NO_SYMMETRY = FOUND .OR. NO_SYMMETRY_LINRESP
C
         CALL SECTION_FIND_KEYWORD('NOSYMCLU',NO_SYMMETRY_CLU)
C
         CALL SECTION_FIND_KEYWORD('NOSYMMSS',FOUND)
         IF ( FOUND ) SYMMETRIZE_MSS = .FALSE.
C
         CALL SECTION_FIND_KEYWORD('NOSYMRHO',FOUND)
         IF ( FOUND ) SYMMETRIZE_RHO = .FALSE.
C_
         CALL SECTION_SET_INTEGER('TEST',ITEST,9999,0)
         IF ( ITEST.LT.0 .OR. ITEST.GT.ITESTMAX ) THEN
            CALL STOP_MESSAGE(ROUTINE,
     &                        'ITEST out of range  1 ... ITESTMAX')
         ELSE IF ( (ITEST.EQ.6 .OR. ITEST.EQ.7) .AND. 
     &             PROGNAME.NE.'KKRGEN    ' ) THEN
            WRITE (6,99001) TXTTEST(ITEST)
            CALL STOP_MESSAGE(ROUTINE,
     &                        'test has to be performed using KKRGEN ')
         END IF
         IF ( ITEST.EQ.8 .AND. 
     &        .NOT.(PROGNAME.NE.'KKRGEN    ' .OR. PROGNAME.NE.
     &        'KKRSCF    ') ) THEN
            WRITE (6,99001) TXTTEST(ITEST)
            CALL STOP_MESSAGE(ROUTINE,
     &                'test has to be performed using KKRSCF or KKRGEN '
     &                )
         END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT X
C       setting a BREAKPOINT by setting  BREAK = X  in section CONTROL
C       forces additional output and stop of the program
C       see CASE list in <STOP_BREAKPOINT> for the various BREAKPONTs
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT X
C
         CALL SECTION_SET_INTEGER('BREAK',BREAKPOINT,9999,0)
C
      END IF
C
      IF ( (ITEST.EQ.2) .OR. (ITEST.EQ.3) ) THEN
         RUNELOOP = .TRUE.
         IBZINT = -1
         SYMCHECK = .TRUE.
      END IF
      IF ( ITEST.GT.ITESTMAX )
     &      CALL STOP_MESSAGE(ROUTINE,'ITEST should be <=5')
C
      IF ( RDTAU ) RUNELOOP = .FALSE.
C --------------------------------------------- allow RDTAU AND RUNELOOP
      CALL INPUT_FIND_SECTION('CONTROL',0)
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('RUNELOOP',FOUND)
         IF ( FOUND ) RUNELOOP = .TRUE.
      END IF
      IF ( .NOT.RUNELOOP ) INITELOOP = .FALSE.
C ----------- allow TAU for the final states to be calculated in <VBPES>
      IF ( TASK(1:5).EQ.'VBPES' .OR. TASK(1:5).EQ.'ARPES' )
     &     INITELOOP = .TRUE.
C
      IF ( TASK(1:3).EQ.'BSF' ) INITELOOP = .TRUE.
C
      IF ( TASK(1:6).EQ.'TORQUE' ) INITELOOP = .TRUE.
C
      IF ( TASK(1:7).EQ.'COMPTON' ) INITELOOP = .TRUE.
C
      IF ( TASK(1:5).EQ.'EKREL' ) INITELOOP = .TRUE.
C
      IF ( TASK(1:3).EQ.'AES' ) INITELOOP = .TRUE.
C
      IF ( TASK(1:5).EQ.'CLXPS' ) INITELOOP = .TRUE.
C
      IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &     TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' .OR. 
     &     IXRAY.NE.0 ) INITELOOP = .TRUE.
C
C=======================================================================
C                     read  MODE  parameters
C=======================================================================
C
      IREL = 3
C--------------- use NQCLU as indicator for embedded cluster calculation
      IF ( NQCLU.GT.0 ) IREL = IREL_HOST
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
C
C-----------------------------------------------------------------------
C                        Lloyd formula
C-----------------------------------------------------------------------
C     NOTE: LLOYD activated for  WRNOS
C
         IF ( .NOT.LLOYD ) CALL SECTION_FIND_KEYWORD('LLOYD',LLOYD)
C
C-----------------------------------------------------------------------
C              use ZJOVZZ to avoid irregular solution
C-----------------------------------------------------------------------
C
         CALL SECTION_FIND_KEYWORD('ZJOVZZ',USE_FPCHR_ZJOVZZ)
C
C-----------------------------------------------------------------------
C                        Breit interaction
C-----------------------------------------------------------------------
C
         CALL SECTION_FIND_KEYWORD('BREITINT',BREITINT5)
C
C-----------------------------------------------------------------------
C          thermal lattice vibrations or spin fluctuations
C-----------------------------------------------------------------------
C
         CALL SECTION_FIND_KEYWORD('THERMAL',THERMAL_VIBRA_FLUCT)
C
C                                          allow for OLD key word LATVIB
         IF ( .NOT.THERMAL_VIBRA_FLUCT )
     &         CALL SECTION_FIND_KEYWORD('LATVIB',THERMAL_VIBRA_FLUCT)
C
C-----------------------------------------------------------------------
C               calculational mode  --- fix IREL
C-----------------------------------------------------------------------
C
         CALL SECTION_FIND_KEYWORD('NREL',FOUND)
         IF ( FOUND ) IREL = 0
         CALL SECTION_FIND_KEYWORD('SREL',FOUND)
         IF ( FOUND ) IREL = 1
         CALL SECTION_FIND_KEYWORD('SP-SREL',FOUND)
         IF ( FOUND ) IREL = 2
C
      END IF
C
      IF ( TASK(1:5).EQ.'APS  ' ) IREL = MIN(1,IREL)
      IF ( TASK(1:5).EQ.'AES  ' ) IREL = MIN(3,IREL)
      IF ( TASK(1:5).EQ.'NRAES' ) IREL = MIN(1,IREL)
      IF ( TASK(1:5).EQ.'RELAX' ) IREL = MIN(1,IREL)
      IF ( TASK(1:6).EQ.'WFPLOT' ) IREL = 3
C
C-----------------------------------------------------------------------
C         fix common Euler angles for magnetic moment
C-----------------------------------------------------------------------
C
      IF ( .NOT.ROTATE_LATTICE ) THEN
C
         CALL INPUT_FIND_SECTION('MODE',0)
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_REAL('MALF',ALFDEG,9999D0,0)
            IF ( FOUND_REAL ) KMROT_COMMON = 1
            CALL SECTION_SET_REAL('MBET',BETDEG,9999D0,0)
            IF ( FOUND_REAL ) KMROT_COMMON = 1
            CALL SECTION_SET_REAL('MGAM',GAMDEG,9999D0,0)
            IF ( FOUND_REAL ) KMROT_COMMON = 1
C
            CALL SECTION_SET_REAL_ARRAY('MDIR',MDIR,N_FOUND,3,0,9999D0,
     &                                  0)
            IF ( FOUND_REAL_ARRAY ) THEN
               CALL CONVERT_CART_TO_SPHER(MDIR(1),MDIR(2),MDIR(3),
     &            RDUMMY,BETDEG,ALFDEG,.TRUE.)
               GAMDEG = 0.0D0
               KMROT_COMMON = 1
            END IF
         END IF
C
      ELSE
C
         KMROT_COMMON = 0
         KMROT = 0
C
      END IF
C
C-----------------------------------------------------------------------
C                               KMROT
C-----------------------------------------------------------------------
C  0:   no rotation of the magnetisation
C  1:   individual rotation of the magnetisation for every site IQ
C  2:   global COMMON rotation of the magnetisation
C  3:   spin spiral    Theta =  90 }    QMVEC <> 0-vector
C  4:   spin spiral    Theta <> 90 }    IREL = 2
C-----------------------------------------------------------------------
C
      IF ( (KMROT.NE.0) .AND. (KMROT_COMMON.NE.0) ) WRITE (6,99003)
C
      IF ( (KMROT.EQ.0) .AND. (KMROT_COMMON.NE.0) ) THEN
         DO IQ = 1,NQ
            QMPHI(IQ) = ALFDEG
            QMTET(IQ) = BETDEG
            QMGAM(IQ) = GAMDEG
         END DO
         KMROT = 2
      END IF
C
      DO IQ = 1,NQMAX
         MDIRQ(1,IQ) = 0.0D0
         MDIRQ(2,IQ) = 0.0D0
         MDIRQ(3,IQ) = 1.0D0
      END DO
C
      CALL INPUT_FIND_SECTION('MODE',0)
      IF ( FOUND_SECTION ) THEN
         DO IQ = 1,NQ
            STR10 = 'MDIR'
            CALL STRING_ADD_N(STR10,IQ)
            CALL SECTION_SET_REAL_ARRAY(STR10,MDIRQ(1,IQ),N_FOUND,3,1,
     &                                  9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
C
               KMROT = 2
C
               CALL CONVERT_CART_TO_SPHER(MDIRQ(1,IQ),MDIRQ(2,IQ),
     &            MDIRQ(3,IQ),RDUMMY,QMTET(IQ),QMPHI(IQ),.TRUE.)
C
            END IF
         END DO
C
C---------------------------------------------------------- spin spirals
         CALL SECTION_SET_REAL_ARRAY('QMVEC',QMVEC,N_FOUND,3,0,9999D0,0)
C
         IF ( FOUND_REAL_ARRAY ) THEN
            KMROT = 3
            DO IQ = 1,NQ
               IF ( ABS(90D0-QMTET(IQ)).GT.0.001D0 ) KMROT = 4
            END DO
         END IF
      END IF
C
      IF ( KMROT.NE.0 .AND. KMROT.LE.2 ) THEN
         KMROT = 2
         DO IQ = 2,NQ
            IF ( ABS(QMTET(IQ)-QMTET(1)).GT.1D-6 ) KMROT = 1
            IF ( ABS(QMPHI(IQ)-QMPHI(1)).GT.1D-6 ) KMROT = 1
         END DO
      END IF
C
C--------------------------- do a final test and adjust KMROT eventually
C
      CALL MROT_CHECK_SETTINGS(NQ,IREL,QMPHI,QMTET,QMGAM,KMROT,QMVEC,
     &                         NQMAX)
C
      IF ( (ITEST.EQ.6 .OR. ITEST.EQ.7) .AND. IREL.NE.3 ) THEN
         WRITE (6,99002) TXTTEST(ITEST)
         IREL = 3
      END IF
C
      IF ( ITEST.EQ.8 ) IREL = 0
C
      IF ( BREAKPOINT.NE.0 ) WRITE (6,99006) BREAKPOINT
C
C=======================================================================
C
C------------------------------------------- NONMAG / POTFMTOUT / SETPOT
C
      POTFMTOUT = MAX(5,POTFMTINP,POTFMTOUT)
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('NONMAG',NONMAG)
         CALL SECTION_SET_INTEGER('POTFMT',POTFMTOUT,9999,0)
         CALL SECTION_SET_STRING('SETPOT',SETPOT,'9999',0)
      END IF
C
      IF ( IREL.LE.1 ) NONMAG = .TRUE.
C
C=======================================================================
C         use a constant potential for testing purposes
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL('POTFIX',POTFIX,9999D0,0)
         USE_CONST_POTENTIAL = FOUND_REAL
         IF ( USE_CONST_POTENTIAL ) THEN
            WRITE (6,99004) POTFIX
            WRNOS = .TRUE.
            LLOYD = .TRUE.
            IF ( PROGNAME.NE.'KKRGEN    ' ) THEN
               WRITE (6,99005) PROGNAME
               CALL STOP_MESSAGE(ROUTINE,
     &                           'TASK POTFIX has to be run by KKRGEN')
            END IF
         END IF
      END IF
C
C=======================================================================
C        set USENLCPA - needed in case of SIGMA + NLCPA
C=======================================================================
      IF ( NCPA.NE.1 ) THEN
         CALL INPUT_FIND_SECTION('CPA',0)
         IF ( FOUND_SECTION )
     &        CALL SECTION_FIND_KEYWORD('NLCPA',USENLCPA)
      END IF
C=======================================================================
99001 FORMAT (2(/,1X,79('*')),/,10X,'test option selected:',/,10X,A40,
     &        2(/,1X,79('*')),/)
99002 FORMAT (2(/,1X,79('*')),/,10X,'test option selected:',/,10X,A40,
     &        //,10X,'IREL reset to 3 !!!!!!!!!!!',/,2(/,1X,79('*')),/)
99003 FORMAT (//,70('*'),/,10X,'WARNING from  <<INITALL>> ',/,10X,
     &        'rotation angles read from potential file',/,10X,
     &        'global rotation angles in input file ignored',/,70('*'),
     &        /)
99004 FORMAT (//,2(1X,79('T')),/,10X,'potential will be set to a ',
     &        'constant value for testing purposes',/,10X,'POTFIX = ',
     &        F10.5,/,2(1X,79('T')),//)
99005 FORMAT (//,10X,'PROGNAME = ',A,'  test runs only for  KKRGEN !',/)
99006 FORMAT (2(/,1X,79('B')),/,10X,'BREAKPOINT SET TO:',I5,//,10X,
     &  'setting a BREAKPOINT by setting  BREAK = X  in section CONTROL'
     &  ,/,10X,'forces additional output and stop of the program',/,10X,
     &  'see CASE list in <STOP_BREAKPOINT> for the various BREAKPONTs',
     &  //,2(/,1X,79('B')),/)
99007 FORMAT (10X,'tabulated      PI = ',E40.30,/,10X,
     &        'calculated     PI = ',E40.30,/,10X,
     &        'machine accuracy  = ',E40.30,/)
C
C
      END
