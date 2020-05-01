C*==init_mod_calcmode1.f    processed by SPAG 6.70Rc at 08:53 on  8 Mar 2017
      SUBROUTINE INIT_MOD_CALCMODE1(BEXT)
C   ********************************************************************
C   *                                                                  *
C   *  initialize variables that control the calculation mode          *
C   *                                                                  *
C   *  STEP 1:  ORBPOL, BLCOUPL, EXTFIELD                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FOUND_SECTION,FOUND_REAL
      USE MOD_CALCMODE,ONLY:ORBPOL,BLCOUPL,EXTFIELD,U_POT_SHIFT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BEXT
C
C Local variables
C
      LOGICAL FOUND
C
C*** End of declarations rewritten by SPAG
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_SET_REAL('BEXT',BEXT,9999D0,0)
         EXTFIELD = FOUND_REAL
         IF ( EXTFIELD ) THEN
            CALL SECTION_FIND_KEYWORD('NOBLC',FOUND)
            IF ( FOUND ) BLCOUPL = .FALSE.
         ELSE
            BLCOUPL = .FALSE.
         END IF
C
         CALL SECTION_SET_STRING('OP',ORBPOL,'9999',0)
C
C-----------------------------------------------------------------------
         CALL SECTION_SET_REAL('UVSHIFT',U_POT_SHIFT,9999D0,0)
C
      END IF
      END
C*==init_mod_calcmode2.f    processed by SPAG 6.70Rc at 08:53 on  8 Mar 2017
      SUBROUTINE INIT_MOD_CALCMODE2(ITEST,KMROT_COMMON,ALFDEG,BETDEG,
     &                              GAMDEG,BEXT,FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *  initialize variables that control the calculation mode          *
C   *                                                                  *
C   *  STEP 2                                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT
      USE MOD_FILES,ONLY:INFO,LINFO,FOUND_SECTION,FOUND_REAL,
     &    FOUND_STRING,FOUND_REAL_ARRAY,FOUND_STRING_ARRAY
      USE MOD_ANGMOM,ONLY:NL,NLMAX
      USE MOD_CONSTANTS,ONLY:C_AU
      USE MOD_TYPES,ONLY:LOPT,SOCTL,CTL,NT,Z,NTMAX
      USE MOD_CALCMODE,ONLY:SOLVER,ORBPOL,BLCOUPL,BREITINT,EXTFIELD,
     &    DMFT,KKRMODE,IREL,SOLVER_FP,TOL_DIRBS,TOL_FPDIRBS,
     &    TOL_FPDIRBORN,TOL_FPRWFBORN,TOL_RWFBS,TOL_FPRWFBS,
     &    LIM_FPDIRBORN,LIM_FPRWFBORN,LDAU,GF_CONV_RH,
     &    CALC_KINETIC_ENERGY_DENSITY
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_CALCMODE2')
C
C Dummy arguments
C
      REAL*8 ALFDEG,BETDEG,BEXT,GAMDEG
      LOGICAL FULLPOT5
      INTEGER ITEST,KMROT_COMMON
C
C Local variables
C
      REAL*8 AUXL(:),BT(NRMAX,NTMAX),SCLBXC(:)
      LOGICAL FOUND,MANPULBXC,MANPULC,MANPULSOC,UPDATE
      INTEGER I,IL,IT,NIN
      CHARACTER*2 SOCII(-2:-1)
      CHARACTER*10 STR10,WFACCURACY
      CHARACTER*5 STR5
      CHARACTER*1 STRLOPT(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA SOCII/'xy','zz'/
C
      ALLOCATABLE AUXL,STRLOPT,SCLBXC
C
      ALLOCATE (SCLBXC(NTMAX),STRLOPT(NTMAX),AUXL(NLMAX))
C
C-----------------------------------------------------------------------
      MANPULC = .FALSE.
      MANPULSOC = .FALSE.
      MANPULBXC = .FALSE.
C
      DO IT = 1,NTMAX
         SCLBXC(IT) = 1.0D0
         DO IL = 1,NLMAX
            SOCTL(IT,IL) = 1.0D0
            CTL(IT,IL) = C_AU
         END DO
      END DO
C
      SOLVER = 'BS        '
      IF ( BREITINT ) SOLVER = 'BS-BI     '
      IF ( ORBPOL.NE.'NONE' ) SOLVER = 'ABM-OP    '
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_STRING('SOLVER',SOLVER,
     &     '9999',0)
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_SET_REAL('BEXT',BEXT,9999D0,0)
         EXTFIELD = FOUND_REAL
         IF ( EXTFIELD ) THEN
            CALL SECTION_FIND_KEYWORD('NOBLC',FOUND)
            IF ( FOUND ) BLCOUPL = .FALSE.
         ELSE
            BLCOUPL = .FALSE.
         END IF
C
         CALL SECTION_SET_STRING('OP',ORBPOL,'9999',0)
C
C=======================================================================
C                         DMFT-calculational mode
C=======================================================================
C
         IF ( ORBPOL(1:4).EQ.'DMFT' ) THEN
            DMFT = .TRUE.
C
            SOLVER = 'BS        '
         ELSE IF ( ORBPOL(1:5).EQ.'LDA+U' ) THEN
            LDAU = .TRUE.
C
C           SOLVER = 'BORN      '
            SOLVER = 'BS        '
         END IF
C=======================================================================
C
         IF ( ORBPOL(1:4).NE.'NONE' .OR. BLCOUPL ) THEN
            IF ( BLCOUPL ) SOLVER = 'ABM-OP    '
            IF ( ORBPOL(1:6).EQ.'BROOKS' ) SOLVER = 'ABM-OP    '
C
            IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. ORBPOL(1:5)
     &           .EQ.'SIGMA' .OR. BLCOUPL ) THEN
C
               CALL SECTION_SET_STRING_ARRAY('LOPT',STRLOPT,NIN,NTMAX,0,
     &            '9999',0)
C
               IF ( .NOT.FOUND_STRING_ARRAY .OR. NIN.NE.NT ) THEN
                  DO IT = 1,NT
                     IF ( 22.LE.Z(IT) .AND. Z(IT).LE.28 ) LOPT(IT) = 2
                     IF ( 40.LE.Z(IT) .AND. Z(IT).LE.46 ) LOPT(IT) = 2
                     IF ( 72.LE.Z(IT) .AND. Z(IT).LE.78 ) LOPT(IT) = 2
                     IF ( 58.LE.Z(IT) .AND. Z(IT).LE.71 ) LOPT(IT) = 3
                     IF ( 89.LE.Z(IT) ) LOPT(IT) = 3
                  END DO
               ELSE
                  DO IT = 1,NT
                     IF ( STRLOPT(IT).EQ.'d' .OR. STRLOPT(IT).EQ.'2' )
     &                    THEN
                        LOPT(IT) = 2
                     ELSE IF ( STRLOPT(IT).EQ.'f' .OR. STRLOPT(IT)
     &                         .EQ.'3' ) THEN
                        LOPT(IT) = 3
                     ELSE
                        LOPT(IT) = -1
                     END IF
                  END DO
               END IF
            END IF
C
         END IF
C
         CALL SECTION_SET_REAL_ARRAY('BXC',SCLBXC,NIN,NTMAX,0,9999D0,0)
C
C=======================================================================
C  specify whether the SOC should be manipulated l- and type-dependent
C    MANPULC:      via CTL     scale the speed of light
C    MANPULSOC:    via SOCTL   scale strength of SOC
C                              or select between xy- and zz-SOC term
C=======================================================================
         IF ( IREL.EQ.3 ) THEN
C
C------------------------------ apply the same manipulation to ALL types
C
            STR10 = 'C'
            CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,NIN,NLMAX,0,9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
               DO IT = 1,NT
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        CTL(IT,IL) = CTL(IT,IL)/AUXL(IL)**0.5D0
                        MANPULC = .TRUE.
                     END IF
                  END DO
               END DO
            END IF
C
            STR10 = 'SOC'
            CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,NIN,NLMAX,0,9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
               DO IT = 1,NT
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-I  '
                     END IF
                     IF ( AUXL(IL).LT.0D0 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-II '
                     END IF
                  END DO
               END DO
            END IF
C
C------------------ apply the manipulation to individually for the types
C
            DO IT = 1,NT
               STR10 = 'C'
               CALL STRING_ADD_N(STR10,IT)
               CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,NIN,NLMAX,0,
     &            9999D0,0)
C
               IF ( FOUND_REAL_ARRAY ) THEN
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        CTL(IT,IL) = CTL(IT,IL)/AUXL(IL)**0.5D0
                        MANPULC = .TRUE.
                     END IF
                  END DO
               END IF
C
               STR10 = 'SOC'
               CALL STRING_ADD_N(STR10,IT)
               CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,NIN,NLMAX,0,
     &            9999D0,0)
C
               IF ( FOUND_REAL_ARRAY ) THEN
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-I  '
                     END IF
                     IF ( AUXL(IL).LT.0D0 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-II '
                     END IF
                  END DO
               END IF
            END DO
         END IF
C
      END IF
C
C--------------- set SOC = 0 for spin-polarized scalar relativistic case
C------------------------------------ including  spin sprirals KMROT=3,4
C
      IF ( IREL.EQ.2 ) THEN
         DO IT = 1,NT
            DO IL = 1,NLMAX
               SOCTL(IT,IL) = 0.0D0
            END DO
         END DO
      END IF
C
      WRITE (6,FMT='(/,1X,79(''m''),/)')
      IF ( IREL.EQ.0 ) THEN
         WRITE (6,99001) 'NON-relativistic calculation '
      ELSE IF ( IREL.EQ.1 ) THEN
         WRITE (6,99001) 'SCALAR-relativistic calculation '
      ELSE IF ( IREL.EQ.2 ) THEN
         WRITE (6,99001) 
     &                 'spin polarized SCALAR-relativistic calculation '
      ELSE
         WRITE (6,99001) 'FULLY relativistic calculation '
      END IF
C
      IF ( KMROT_COMMON.NE.0 ) WRITE (6,99002) ALFDEG,BETDEG,GAMDEG
C
      IF ( EXTFIELD ) THEN
         WRITE (6,'(/,10X,''an external magnetic field is present'')')
         WRITE (6,99003) 'BEXT(Ry):    ',BEXT
         IF ( BLCOUPL ) THEN
            WRITE (6,99012) 'B - l_z coupling'
         ELSE
            WRITE (6,99011) 'B - l_z coupling'
         END IF
         INFO = INFO(1:LINFO)//' B-ext'
         LINFO = LINFO + 6
      END IF
C
      IF ( ITEST.EQ.6 ) THEN
         SOLVER = 'BS        '
         MANPULC = .TRUE.
         DO IT = 1,NTMAX
            DO IL = 1,NLMAX
               CTL(IT,IL) = C_AU*1.0D4
            END DO
         END DO
      END IF
C
      IF ( MANPULC ) THEN
         WRITE (6,'(/,10X,''the speed of light is scaled with:'')')
         DO IT = 1,NT
            WRITE (6,99004) IT,(NL-1),((C_AU/CTL(IT,IL))**2,IL=1,NL)
         END DO
         INFO = INFO(1:LINFO)//' C var.'
         LINFO = LINFO + 7
      END IF
      IF ( MANPULSOC ) THEN
         IF ( SOLVER.EQ.'BS-SOC-I  ' ) THEN
            WRITE (6,'(/,10X,''the SOC is scaled with:'')')
            DO IT = 1,NT
               WRITE (6,99004) IT,(NL-1),(SOCTL(IT,IL),IL=1,NL)
               DO IL = 1,NL
                  IF ( SOCTL(IT,IL).LT.0D0 ) CALL STOP_MESSAGE(ROUTINE,
     &                 'all SOCTL should be >= 0 >>>> check input file')
               END DO
            END DO
            INFO = INFO(1:LINFO)//' SOC var.'
            LINFO = LINFO + 9
            SOLVER = 'BS-SOC    '
         ELSE
            WRITE (6,99005)
            DO IT = 1,NT
               WRITE (6,99006) IT,NL - 1,
     &                         (SOCII(NINT(SOCTL(IT,IL))),IL=1,NL)
               DO IL = 1,NL
                  IF ( SOCTL(IT,IL).GT.-0.1D0 )
     &                 CALL STOP_MESSAGE(ROUTINE,
     &                 'all SOCTL should be < 0 >>>> check input file')
               END DO
            END DO
            INFO = INFO(1:LINFO)//' SOC-II  '
            LINFO = LINFO + 9
            SOLVER = 'BS-SOC    '
         END IF
      END IF
      IF ( MANPULBXC ) THEN
         WRITE (6,'(/,10X,''the B_xc - field is scaled with:'')')
         WRITE (6,99007) (IT,SCLBXC(IT),IT=1,NT)
         INFO = INFO(1:LINFO)//' BXC var.'
         LINFO = LINFO + 9
         DO IT = 1,NT
            IF ( ABS(SCLBXC(IT)-1D0).GT.1D-5 ) THEN
               DO I = 1,NRMAX
                  BT(I,IT) = SCLBXC(IT)*BT(I,IT)
               END DO
            END IF
         END DO
      END IF
C
      IF ( BREITINT ) THEN
         WRITE (6,99012) 'BREIT interaction'
         INFO = INFO(1:LINFO)//' BREIT'
         LINFO = LINFO + 6
      ELSE
         WRITE (6,99011) 'BREIT interaction'
      END IF
C
      IF ( ORBPOL(1:4).EQ.'NONE' ) THEN
         WRITE (6,99011) 'ORBITAL POLARISATION'
      ELSE
         WRITE (6,99012) 'ORBITAL POLARISATION',ORBPOL
         IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
            WRITE (6,99008) (LOPT(IT),IT=1,NT)
         ELSE IF ( ORBPOL(1:6).EQ.'SIGMA' ) THEN
            WRITE (6,99008) (LOPT(IT),IT=1,NT)
C
         ELSE IF ( .NOT.DMFT .OR. .NOT.LDAU ) THEN
            WRITE (6,99010) ORBPOL
C
         ELSE
            WRITE (6,99010) ORBPOL
            STOP 'in <INIT_MOD_CALCMODE>'
         END IF
C
         IF ( BREITINT ) THEN
            WRITE (6,99009)
            BREITINT = .FALSE.
         END IF
         INFO = INFO(1:LINFO)//' OP:'//ORBPOL
         LINFO = LINFO + 14
      END IF
C
C=======================================================================
C                                                              IREL <= 2
      IF ( IREL.LE.2 ) THEN
C
         SOLVER_FP = 'BORN'
C
         IF ( IREL.EQ.0 ) THEN
            STR5 = 'NREL '
         ELSE
            STR5 = 'SREL '
         END IF
C
C-----------------------------------------------------------------------
         CALL INPUT_FIND_SECTION('CONTROL',0)
C
         IF ( FOUND_SECTION ) THEN
C
            CALL SECTION_SET_STRING('WFACCURACY',WFACCURACY,'9999',0)
C
            IF ( FOUND_STRING ) THEN
               SELECT CASE (WFACCURACY)
               CASE ('LOW')
                  TOL_RWFBS = 2.0D-6
                  TOL_FPRWFBS = 2.0D-6
                  TOL_FPRWFBORN = 1.0D-6
                  LIM_FPRWFBORN = 10
               CASE ('MEDIUM')
                  TOL_RWFBS = 2.0D-8
                  TOL_FPRWFBS = 2.0D-8
                  TOL_FPRWFBORN = 1.0D-10
                  LIM_FPRWFBORN = 10
               CASE ('HIGH')
                  TOL_RWFBS = 2.0D-12
                  TOL_FPRWFBS = 2.0D-12
                  TOL_FPRWFBORN = 1.0D-12
                  LIM_FPRWFBORN = 20
               CASE ('HIGHEST')
                  TOL_RWFBS = 2.0D-13
                  TOL_FPRWFBS = 2.0D-13
                  TOL_FPRWFBORN = 1.0D-14
                  LIM_FPRWFBORN = 40
               CASE DEFAULT
                  TOL_RWFBS = 2.0D-8
                  TOL_FPRWFBS = 2.0D-8
                  TOL_FPRWFBORN = 1.0D-10
                  LIM_FPRWFBORN = 10
               END SELECT
            END IF
C
         END IF
C
         CALL SECTION_SET_REAL('TOLBS',TOL_RWFBS,9999D0,0)
C
         IF ( FULLPOT .OR. FULLPOT5 ) THEN
C
            CALL SECTION_SET_STRING('SOLVER_FP',SOLVER_FP,'9999',0)
C
            IF ( SOLVER_FP.EQ.'BS' ) THEN
               CALL SECTION_SET_REAL('TOLFPBS',TOL_FPRWFBS,9999D0,0)
            ELSE
               CALL SECTION_SET_REAL('TOLFPBORN',TOL_FPRWFBORN,9999D0,0)
               CALL SECTION_SET_INTEGER('LIMFPBORN',LIM_FPRWFBORN,9999,
     &                                  0)
            END IF
         END IF
C
C-----------------------------------------------------------------------
C
         IF ( SOLVER.EQ.'RK' ) THEN
            WRITE (6,99013) STR5,SOLVER
         ELSE
            WRITE (6,99013) STR5,SOLVER,TOL_RWFBS
         END IF
C
         IF ( .NOT.(FULLPOT) .OR. .NOT.(FULLPOT5) ) THEN
            IF ( SOLVER_FP.EQ.'BS' ) THEN
               WRITE (6,99014) STR5,SOLVER_FP,TOL_FPRWFBS
            ELSE
               WRITE (6,99015) STR5,SOLVER_FP,TOL_FPRWFBORN,
     &                         LIM_FPRWFBORN
            END IF
         END IF
C
C ----------------- if no special solver hes been selected automatically
C --------------------------------------- allow to select a specific one
      ELSE IF ( (.NOT.MANPULSOC) .AND. (.NOT.BREITINT) .AND. 
     &          (ORBPOL(1:4).EQ.'NONE') ) THEN
C
         CALL INPUT_FIND_SECTION('CONTROL',0)
C
         IF ( FOUND_SECTION ) THEN
C
            CALL SECTION_SET_STRING('WFACCURACY',WFACCURACY,'9999',0)
C
            IF ( FOUND_STRING ) THEN
               SELECT CASE (WFACCURACY)
               CASE ('LOW')
                  TOL_DIRBS = 2.0D-6
                  TOL_FPDIRBS = 2.0D-6
                  TOL_FPDIRBORN = 1.0D-6
                  LIM_FPDIRBORN = 10
               CASE ('MEDIUM')
                  TOL_DIRBS = 2.0D-8
                  TOL_FPDIRBS = 2.0D-8
                  TOL_FPDIRBORN = 1.0D-10
                  LIM_FPDIRBORN = 10
               CASE ('HIGH')
                  TOL_DIRBS = 2.0D-12
                  TOL_FPDIRBS = 2.0D-12
                  TOL_FPDIRBORN = 1.0D-12
                  LIM_FPDIRBORN = 20
               CASE ('HIGHEST')
                  TOL_DIRBS = 2.0D-13
                  TOL_FPDIRBS = 2.0D-13
                  TOL_FPDIRBORN = 1.0D-14
                  LIM_FPDIRBORN = 40
               CASE DEFAULT
                  TOL_DIRBS = 2.0D-8
                  TOL_FPDIRBS = 2.0D-8
                  TOL_FPDIRBORN = 1.0D-10
                  LIM_FPDIRBORN = 10
               END SELECT
            END IF
C
            CALL SECTION_SET_STRING('SOLVER',SOLVER,'9999',0)
C
            IF ( SOLVER(1:2).EQ.'BS' ) THEN
               IF ( MANPULSOC ) THEN
                  SOLVER = 'BS-SOC    '
               ELSE
                  SOLVER = 'BS        '
               END IF
            END IF
C
            CALL SECTION_SET_REAL('TOLBS',TOL_DIRBS,9999D0,0)
C
            IF ( SOLVER.EQ.'RK' ) THEN
               WRITE (6,99013) 'Dirac',SOLVER
            ELSE
               WRITE (6,99013) 'Dirac',SOLVER,TOL_DIRBS
            END IF
C
C-----------------------------------------------------------------------
C                          FP-SPR solver
C-----------------------------------------------------------------------
C
            IF ( FULLPOT .OR. FULLPOT5 ) THEN
C
               SOLVER_FP = 'BORN'
C
               CALL SECTION_SET_STRING('SOLVER_FP',SOLVER_FP,'9999',0)
C
               IF ( SOLVER_FP.EQ.'BS' ) THEN
                  CALL SECTION_SET_REAL('TOLFPBS',TOL_FPDIRBS,9999D0,0)
                  WRITE (6,99014) 'Dirac',SOLVER_FP,TOL_FPDIRBS
               ELSE
                  CALL SECTION_SET_REAL('TOLFPBORN',TOL_FPDIRBORN,
     &                                  9999D0,0)
                  CALL SECTION_SET_INTEGER('LIMFPBORN',LIM_FPDIRBORN,
     &               9999,0)
C
                  WRITE (6,99015) 'Dirac',SOLVER_FP,TOL_FPDIRBORN,
     &                            LIM_FPDIRBORN
               END IF
C
            END IF
C
         END IF
C=======================================================================
      END IF
C=======================================================================
C
      WRITE (6,99016) SOLVER
      IF ( FULLPOT ) WRITE (6,99017) SOLVER_FP
C
C-----------------------------------------------------------------------
C                        convention for Green function
C-----------------------------------------------------------------------
C     DEFAULT: SOLVER_FP='BORN'  implies:   GF_CONV_RH = .TRUE. !!!!!
C              else                         GF_CONV_RH = .FALSE.
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('GF_CONV_RH',UPDATE)
         IF ( UPDATE ) GF_CONV_RH = .TRUE.
      END IF
C
      IF ( FULLPOT .OR. FULLPOT5 ) THEN
         IF ( SOLVER_FP.EQ.'BORN' ) GF_CONV_RH = .TRUE.
      END IF
C
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         KKRMODE = 'STANDARD-KKR'
      ELSE
         KKRMODE = 'TB-KKR    '
      END IF
C
      CALL INPUT_FIND_SECTION('TAU',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_STRING('KKRMODE',KKRMODE,
     &     '9999',0)
C
C------------------------------------------------ deal with old settings
C
      IF ( KKRMODE(1:10).EQ.'STANDARD  ' ) THEN
         KKRMODE = 'STANDARD-KKR'
      ELSE IF ( KKRMODE(1:10).EQ.'TB        ' ) THEN
         KKRMODE = 'TB-KKR    '
      END IF
C
C-----------------------------------------------------------------------
C                 calculation of kinetic energy density
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      CALL SECTION_FIND_KEYWORD('EKINDNS',CALC_KINETIC_ENERGY_DENSITY)
C
C=======================================================================
99001 FORMAT (10X,A,5I10)
99002 FORMAT (/,10X,
     &   'the magnetization is rotated globally using the Euler angles '
     &   ,//,10X,'ALF =',F8.3,3X,'BET =',F8.3,3X,'GAM =',F8.3)
99003 FORMAT (10X,A,5F10.6)
99004 FORMAT (10X,'IT=',I2,' L=0,..,',I1,':',20F7.3)
99005 FORMAT (/,10X,'the SOC is manipulated -- part of the SOC kept: ')
99006 FORMAT (10X,'IT=',I2,' L=0,..,',I1,':',20(2X,A))
99007 FORMAT (10X,'IT=',I2,':',F5.2,20(I4,':',F5.2))
99008 FORMAT (/,10X,'LOPT =',30I2,(16X,30I2))
99009 FORMAT (/,10X,'WARNING:  the BREIT interaction will be suppressed'
     &        )
99010 FORMAT (/,10X,'orbital polarisation scheme  ',A,'  not known')
99011 FORMAT (/,10X,'the  ',A,'  will be ignored')
99012 FORMAT (/,10X,'the  ',A,'  will be accounted for':,/,10X,
     &        'using the scheme  ',A,/)
99013 FORMAT (/,10X,'solver for     spherical ',A,' equation: ',A,:,3X,
     &        'TOL =',1P,E8.1)
99014 FORMAT (10X,'solver for NON-spherical ',A,' equation: ',A,:,3X,
     &        'TOL =',1P,E8.1)
99015 FORMAT (10X,'solver for NON-spherical ',A,' equation: ',A,:,3X,
     &        'TOL =',1P,E8.1,/,64X,'LIM = ',I7)
99016 FORMAT (/,10X,'solver for spherical equation:     ',A)
99017 FORMAT (10X,'solver for non-spherical equation: ',A)
      END
