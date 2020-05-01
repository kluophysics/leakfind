C*==thermal_init.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C   ********************************************************************
C   *                                                                  *
C   *  set up  NVIBRA  vectors     SVEC_VT         for displacements   *
C   *          NFLUCT  directions  FTET_FT,PHI_FT  for fluctuations    *
C   *          and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ       *
C   *                                                                  *
C   *  to do THERMAL calculations specify TEMPERATURE (range) in input *
C   *                                                                  *
C   *  NVIBRA NFLUCT set in input                                      *
C   *                                                                  *
C   *    N      N     no temperatrure effects                          *
C   *    Y      N     displacements                                    *
C   *    N      Y     fluctuations                                     *
C   *    Y      Y     displacements  AND  fluctuations                 *
C   *                                                                  *
C   *  at the moment NVIBRA=NVIBRA_TAB is used if NVIBRA>1 is set      *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   * IT_MAG(IT) = 1 when magnetic moment M(IT) > 0                    *
C   * IT_MAG(IT) = 0 when magnetic moment M(IT) = 0                    *
C   * NFLUCT_T(IT) - number of fluctuations applied to atom IT         *
C   * NVIBFLU_T(IT)- number of (fluctuations and vibrations)           *
C   *                applied to atom IT                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:IWFD
      USE MOD_THERMAL,ONLY:NVFT,NVFTMAX,NVT,NVTMAX,NFT,NFTMAX,SVEC_VT,
     &    SVEC_VIBRA0,NVIBFLU,NVIBRA,NVIBRA_TAB,FTET_FT,FPHI_FT,NFLUCT,
     &    THERMAL_CHECK_FLUCT,FLUCT_DIR_SETTING,NFLUCT_T,NVIBFLU_T,
     &    NPHI_FLUCT,NTET_FLUCT,PHI_FLUCT,TET_FLUCT,X_FT,EHAT_FLUCT,
     &    W0_FLUCT,W0_TET,NTET_FLUCT_POT
      USE MOD_TYPES,ONLY:NT,NTMAX,W_WEISS_T,OBS_T,Z
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_FILES,ONLY:IPRINT,IFILFLUCT,FLUCTFIL,POTFIL,LPOTFIL,
     &    FOUND_SECTION,FOUND_REAL,FOUND_INTEGER,FOUND_STRING,N_FOUND
      USE MOD_SCF,ONLY:SCFMIX,SCF_THETA_DEPENDENT_POT
c modified by XJQ: scf of vibrations
      use mod_scfvb_cpa_sigma,only:lscfvb
c end-mod-xjq
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT')
C
C Dummy arguments
C
      INTEGER I_TEMP_LAT,N_TEMP_LAT
      REAL*8 TEMP_LAT
C
C Local variables
C
c modified by XJQ: scf of vibrations
      logical found
c end-mod-xjq
      REAL*8 DMAGDT,MASS,MM_MC(:),MM_T0(:),RMSU,RMSU0_IT(:),RMSU_IT(:),
     &       RNORM,TCURIE,TEMP_LAT_MAX,TEMP_LAT_MIN,T_DEBYE,WFLUCT,
     &       W_WEISS_ERROR
      REAL*8 DNRM2,TAB_EST_W_WEISS
      LOGICAL FLAG_FLUCT,FLAG_FLUCT_FILE,FLAG_VIBRA,INITIALIZE,UDT
      INTEGER I,IFLUCT,IFT,IT,IT_MAG(NT),IVIBRA,IVT,NTET,NT_MAGNETIC
      SAVE DMAGDT,MASS,MM_MC,MM_T0,NT_MAGNETIC,RMSU0_IT,RMSU_IT,TCURIE,
     &     TEMP_LAT_MAX
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./,T_DEBYE/0D0/,FLAG_FLUCT_FILE/.FALSE./
      DATA TEMP_LAT_MIN/0D0/
C
      ALLOCATABLE MM_T0,MM_MC,RMSU0_IT,RMSU_IT
C
      CALL TRACK_INFO(ROUTINE)
C
      NT_MAGNETIC = NT
      DO IT = 1,NT
         IT_MAG(IT) = 1
      END DO
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (NFLUCT_T(NTMAX),NVIBFLU_T(NT))
         ALLOCATE (MM_MC(NTMAX),W_WEISS_T(NTMAX))
C
         NTET_FLUCT = 1
         NPHI_FLUCT = 1
         THERMAL_VIBRA_FLUCT = .FALSE.
C
         CALL INPUT_FIND_SECTION('TASK',0)
C
         IF ( FOUND_SECTION ) THEN
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C                         Thermal spin fluctuations
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
            CALL SECTION_SET_STRING('SETFLUCT',FLUCT_DIR_SETTING,'9999',
     &                              0)
            FLAG_FLUCT = FOUND_STRING
C
            IF ( .NOT.SCF_THETA_DEPENDENT_POT )
     &           CALL SECTION_FIND_KEYWORD('TETDEPPOT',
     &           SCF_THETA_DEPENDENT_POT)
C
            IF ( FLAG_FLUCT ) THEN
C
               WRITE (6,99001) FLUCT_DIR_SETTING
C
               CALL SECTION_SET_INTEGER('NFTET',NTET_FLUCT,9999,0)
               CALL SECTION_SET_INTEGER('NFPHI',NPHI_FLUCT,9999,0)
C
C-----------------------------------------------------------------------
C         reset  NTET_FLUCT  in case of theta-dependent potential
C-----------------------------------------------------------------------
C
               IF ( SCF_THETA_DEPENDENT_POT ) THEN
C
                  IF ( NT_MAGNETIC.LT.NTET_FLUCT ) THEN
C
                     CALL THERMAL_INIT_POT_THETA
C
                     WRITE (6,99012) NT*NTET_FLUCT,POTFIL(1:LPOTFIL-4)
     &                               //'_THETA.pot'
C
                     CALL STOP_REGULAR(ROUTINE,
     &                         'restart program with new potential file'
     &                         )
C
                  END IF
C
                  NTET_FLUCT_POT = NTET_FLUCT
                  NTET_FLUCT = 1
C
               ELSE
C
                  NTET_FLUCT_POT = 1
C
               END IF
C
               NFLUCT = NTET_FLUCT*NPHI_FLUCT
C
               SELECT CASE (FLUCT_DIR_SETTING)
C
C-----------------------------------------------------------------------
C                use temperature dependent magnetisation M(T)
C-----------------------------------------------------------------------
C
               CASE ('M_T       ')
C
                  CALL SECTION_SET_STRING('FLUCTFIL',FLUCTFIL,'9999',0)
C
                  FLAG_FLUCT = FOUND_STRING
                  FLAG_FLUCT_FILE = .TRUE.
C
                  ALLOCATE (MM_T0(NTMAX))
C
                  CALL THERMAL_INIT_FLUCT_M_T(TCURIE,N_TEMP_LAT,
     &               NT_MAGNETIC,MM_T0)
C
                  TEMP_LAT_MAX = TCURIE
C
C-----------------------------------------------------------------------
C                use fluctuation file from Monte Carlo simulations
C-----------------------------------------------------------------------
C
               CASE ('MCS       ')
C
                  CALL SECTION_SET_STRING('FLUCTFIL',FLUCTFIL,'9999',0)
C
                  FLAG_FLUCT = FOUND_STRING
                  FLAG_FLUCT_FILE = .TRUE.
C
                  CALL THERMAL_INIT_FLUCT_MCS(TEMP_LAT_MIN,TEMP_LAT_MAX,
     &               N_TEMP_LAT,NFLUCT)
C
C-----------------------------------------------------------------------
C                       use the RDLM method
C-----------------------------------------------------------------------
C
               CASE ('RDLM      ')
C
                  CALL SECTION_SET_REAL('TMPMIN',TEMP_LAT_MIN,9999D0,1)
                  CALL SECTION_SET_REAL('TMPMAX',TEMP_LAT_MAX,9999D0,1)
                  CALL SECTION_SET_INTEGER('NTMP',N_TEMP_LAT,9999,1)
C
                  CALL SECTION_SET_REAL('TCURIE',TCURIE,0D0,0)
C
                  WRITE (6,99010)
                  DO IT = 1,NT
C
                     W_WEISS_T(IT) = TAB_EST_W_WEISS(Z(IT))
                     OBS_T(0,IWFD,IT) = W_WEISS_T(IT)
                     WRITE (6,99011) IT,W_WEISS_T(IT)
C
                  END DO
C
                  CALL SCF_BROYDEN_W_WEISS(IPRINT,0,SCFMIX,OBS_T,
     &               W_WEISS_ERROR)
C
C-----------------------------------------------------------------------
C          use regular temperature mesh with linear M(T) = M - DMAGDT*T
C-----------------------------------------------------------------------
C
               CASE ('MLIN      ')
C
                  CALL SECTION_SET_REAL('TMPMIN',TEMP_LAT_MIN,9999D0,0)
C
                  CALL SECTION_SET_REAL('TMPMAX',TEMP_LAT_MAX,9999D0,0)
C
                  IF ( .NOT.FOUND_REAL ) THEN
                     TEMP_LAT_MAX = TEMP_LAT_MIN
                     N_TEMP_LAT = 1
                  ELSE
                     CALL SECTION_SET_INTEGER('NTMP',N_TEMP_LAT,9999,0)
                     IF ( .NOT.FOUND_INTEGER ) N_TEMP_LAT = 1
                  END IF
C
                  IF ( NPHI_FLUCT.GT.1 ) THEN
C
                     CALL SECTION_SET_REAL('DMAGDT',DMAGDT,9999D0,1)
C
                  ELSE
C
                     NFLUCT = 1
C
                  END IF
C
C-----------------------------------------------------------------------
C                       no proper value set
C-----------------------------------------------------------------------
C
               CASE DEFAULT
C
                  WRITE (6,99013) FLUCT_DIR_SETTING
C
                  CALL STOP_MESSAGE(ROUTINE,
     &                     'no proper input value for FLUCT_DIR_SETTING'
     &                     )
C
               END SELECT
C=======================================================================
C
               IF ( N_TEMP_LAT.LT.1 ) CALL STOP_MESSAGE(ROUTINE,
     &              'number of temperature steps NTMP < 1 !!!')
C
            END IF
C
C-----------------------------------------------------------------------
C               complete settings related to spin fluctuations
C                    also needed in case of  NFLUCT = 1
C-----------------------------------------------------------------------
C
C
C -------------------- in the case of a system with 'non-magnetic' atoms
C --------------------- the number of 'magnetic' atoms: NT_MAGNETIC < NT
C ------------ and fluctuations are applied only to the 'magnetic' atoms
C
            NFT = NFLUCT*NT_MAGNETIC + (NT-NT_MAGNETIC)
C
            ALLOCATE (PHI_FLUCT(NFLUCT),TET_FLUCT(NFLUCT))
            NTET = MAX(NTET_FLUCT,NTET_FLUCT_POT)
            ALLOCATE (W0_FLUCT(NFLUCT),W0_TET(NTET))
            W0_FLUCT(:) = 0D0
            W0_TET(:) = 0D0
C
            NFTMAX = NFT
C
            ALLOCATE (X_FT(NFTMAX))
C
            WFLUCT = 1D0/DBLE(NFLUCT)
C
C ---------- default setting for the weights of the fluctioan directions
C ------------------------------------------------- X_FT(IFT) = 1/NFLUCT
C
            DO IT = 1,NT
C
               IFT = (IT-1)*NFLUCT
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
                  X_FT(IFT) = WFLUCT
               END DO
C
            END DO
C
            W0_FLUCT(:) = WFLUCT
C
C ----- generate regular fluctuation grid for THETA and PHI as variables
C
            ALLOCATE (EHAT_FLUCT(3,NFLUCT,NTMAX))
C
            CALL THERMAL_FLUCT_GEN_MESH
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
C------------------------------------------- use regular temparture grid
C---------------------------- in case of RDLM or ONLY lattice vibrations
C
            IF ( .NOT.FLAG_FLUCT ) THEN
C
               CALL SECTION_SET_REAL('TMPMIN',TEMP_LAT_MIN,9999D0,0)
               UDT = FOUND_REAL
               THERMAL_VIBRA_FLUCT = FOUND_REAL
               CALL SECTION_SET_REAL('TMPMAX',TEMP_LAT_MAX,9999D0,0)
               UDT = UDT .AND. FOUND_REAL
               THERMAL_VIBRA_FLUCT = THERMAL_VIBRA_FLUCT .OR. FOUND_REAL
               CALL SECTION_SET_INTEGER('NTMP',N_TEMP_LAT,9999,0)
               UDT = UDT .AND. FOUND_INTEGER
               THERMAL_VIBRA_FLUCT = THERMAL_VIBRA_FLUCT .OR. FOUND_REAL
C
               IF ( THERMAL_VIBRA_FLUCT .AND. .NOT.UDT )
     &              CALL STOP_MESSAGE(ROUTINE,
     &                                'specify TMPMIN, TMPMAX, NTMP')
            END IF
C
            THERMAL_VIBRA_FLUCT = THERMAL_VIBRA_FLUCT .OR. FLAG_FLUCT
C
C=======================================================================
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                       Thermal lattice vibrations
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
            CALL SECTION_SET_REAL('TDEBYE',T_DEBYE,9999D0,0)
C
            ALLOCATE (RMSU0_IT(NTMAX),RMSU_IT(NTMAX))
            RMSU0_IT(:) = 0D0
            RMSU_IT(:) = 0D0
C
            IF ( T_DEBYE.LE.0 ) THEN
               CALL SECTION_SET_REAL_ARRAY('RMSU0',RMSU0_IT,N_FOUND,
     &            NTMAX,0,9999D0,0)
               T_DEBYE = DABS(T_DEBYE)
            END IF
C
            CALL THERMAL_RMSDISP(0D0,T_DEBYE,MASS,RMSU)
C
            CALL SECTION_SET_INTEGER('NVIBRA',NVIBRA,9999,0)
C
            WRITE (6,99004) NVIBRA,T_DEBYE
C
            IF ( FOUND_INTEGER ) THEN
               IF ( NVIBRA.LT.1 ) NVIBRA = 1
               IF ( NVIBRA.GT.1 ) NVIBRA = NVIBRA_TAB
            ELSE
               NVIBRA = 1
            END IF
C
            NVT = NVIBRA*NT
C
            FLAG_VIBRA = NVIBRA.NE.1
C
C ------------------------------------------- Thermal lattice vibrations
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C=======================================================================
C
            IF ( .NOT.FLAG_VIBRA .AND. .NOT.FLAG_FLUCT ) THEN
               WRITE (6,99013) FLUCT_DIR_SETTING
               WRITE (6,99014) N_TEMP_LAT,TEMP_LAT_MIN,TEMP_LAT_MAX
C
               IF ( TEMP_LAT_MIN.GT.1.0D-6 ) CALL STOP_ERROR(ROUTINE,
     &              'parameters for THERMAL calcs. incomplete')
C
            END IF
C
            CALL THERMAL_INIT_FLUCT_VIBRA_START
C
            ALLOCATE (SVEC_VT(3,NVTMAX),FTET_FT(NFTMAX),FPHI_FT(NFTMAX))
            WRITE (6,99003) N_TEMP_LAT,TEMP_LAT_MIN,TEMP_LAT_MAX
C
            IF ( NVIBRA.EQ.1 .AND. NFLUCT.EQ.1 )
     &           THERMAL_CHECK_FLUCT = .TRUE.
C
            SVEC_VT(:,:) = 0D0
            FTET_FT(:) = 0D0
            FPHI_FT(:) = 0D0
c modified by XJQ: scf of vibrations
            CALL SECTION_FIND_KEYWORD('SCFVB',FOUND)
            lscfvb = found
c end-mod-xjq
C                                                    THERMAL_VIBRA_FLUCT
C=======================================================================
C
         END IF
C                                                             TASK FOUND
C=======================================================================
C
         IF ( .NOT.THERMAL_VIBRA_FLUCT ) THEN
            WRITE (6,99002) 'NO thermal lattice '//
     &                   'vibrations NOR spin fluctuations to deal with'
C
            NVIBRA = 1
            NFLUCT = 1
            NVT = NT
            NFT = NT
            NVFT = NT
            NVTMAX = NT
            NFTMAX = NT
            NVFTMAX = NT
         END IF
C
         WRITE (6,99005) NVIBFLU,NT,NVT,NVTMAX,NFT,NFTMAX,NVFT,NVFTMAX
C
         INITIALIZE = .FALSE.
C
         RETURN
C
      END IF
C=======================================================================
C                             INITIALIZE                            END
C=======================================================================
C
C=======================================================================
C=======================================================================
C
      IF ( .NOT.THERMAL_VIBRA_FLUCT ) RETURN
C
C=======================================================================
C                 set temperature TEMP_LAT
C=======================================================================
C
      SELECT CASE (FLUCT_DIR_SETTING)
C
C-----------------------------------------------------------------------
      CASE ('M_T       ')
C
         IF ( FLAG_FLUCT_FILE ) THEN
            READ (IFILFLUCT,*) TEMP_LAT,(MM_MC(I),I=1,NT_MAGNETIC)
         ELSE
            CALL STOP_MESSAGE(ROUTINE,
     &                        'specify the M(T) file name in the input '
     &                        )
         END IF
C
         IF ( IPRINT.GT.0 ) WRITE (*,'(9F8.2)') TEMP_LAT,
     &                             (MM_MC(I),I=1,NT_MAGNETIC)
C
         DO I = 1,NT_MAGNETIC
            MM_MC(I) = MM_MC(I)/MM_T0(I)
         END DO
C
C-----------------------------------------------------------------------
      CASE ('MCS       ')
C
         IF ( FLAG_FLUCT_FILE ) THEN
            READ (IFILFLUCT,*) TEMP_LAT
         ELSE
            CALL STOP_MESSAGE(ROUTINE,
     &                     'specify the M-fluct file name in the input '
     &                     )
         END IF
C
C-----------------------------------------------------------------------
      CASE ('RDLM      ')
C
         IF ( N_TEMP_LAT.GT.1 ) THEN
            TEMP_LAT = TEMP_LAT_MIN + (I_TEMP_LAT-1)
     &                 *(TEMP_LAT_MAX-TEMP_LAT_MIN)/DBLE(N_TEMP_LAT-1)
         ELSE
            TEMP_LAT = TEMP_LAT_MIN
         END IF
C
C-----------------------------------------------------------------------
      CASE ('MLIN       ')
C
         IF ( N_TEMP_LAT.GT.1 ) THEN
            TEMP_LAT = TEMP_LAT_MIN + (I_TEMP_LAT-1)
     &                 *(TEMP_LAT_MAX-TEMP_LAT_MIN)/DBLE(N_TEMP_LAT-1)
         ELSE
            TEMP_LAT = TEMP_LAT_MIN
         END IF
C
         DO I = 1,NT_MAGNETIC
            MM_MC(I) = 1D0 - DMAGDT*TEMP_LAT
         END DO
C
      END SELECT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C --------------------------------- set vibration parameters for given T
C
      IF ( NVIBRA.GT.1 ) THEN
C
         IF ( TEMP_LAT.GT.1D0 ) THEN
            CALL THERMAL_RMSDISP(TEMP_LAT,T_DEBYE,MASS,RMSU)
         ELSE
            RMSU = 0D0
         END IF
C
         WRITE (6,99006) I_TEMP_LAT,TEMP_LAT,NVIBRA
         IF ( IPRINT.GE.0 .AND. I_TEMP_LAT.EQ.1 ) WRITE (6,99007)
C
         IVT = 0
         DO IT = 1,NT
            RMSU_IT(IT) = RMSU + RMSU0_IT(IT)
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               SVEC_VT(1:3,IVT) = SVEC_VIBRA0(1:3,IVIBRA)
               RNORM = RMSU_IT(IT)/DNRM2(3,SVEC_VT(1,IVT),1)
               CALL DSCAL(3,RNORM,SVEC_VT(1,IVT),1)
C
               IF ( IPRINT.GE.0 .AND. I_TEMP_LAT.EQ.1 ) WRITE (6,99009)
     &              IT,IVIBRA,SVEC_VT(1:3,IVT)
C
            END DO
         END DO
C
      ELSE
C
         IVT = 0
         DO IT = 1,NT
            IVT = IVT + 1
C
            SVEC_VT(1:3,IVT) = 0D0
C
         END DO
C
      END IF
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
      WRITE (6,99008) I_TEMP_LAT,TEMP_LAT,NFLUCT
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C ------------------------------- set fluctuation parameters for given T
C
      IF ( NFLUCT.GT.1 ) THEN
C
         SELECT CASE (FLUCT_DIR_SETTING)
C
C-----------------------------------------------------------------------
C                use temperature dependent magnetisation M(T)
C-----------------------------------------------------------------------
C
         CASE ('M_T       ')
C
            CALL THERMAL_WEISS_FIELD(TCURIE,TEMP_LAT,MM_MC,NT_MAGNETIC,
     &                               IT_MAG,W_WEISS_T)
C
            DO IT = 1,NT
               WRITE (6,99011) IT,W_WEISS_T(IT)
            END DO
C
C-----------------------------------------------------------------------
C                use fluctuation file from Monte Carlo simulations
C-----------------------------------------------------------------------
C
         CASE ('MCS       ')
C
            WRITE (6,*) ' !!! Under construction !!!'
C
C-----------------------------------------------------------------------
C                       use the RDLM method
C         take Weiss field from last temperature to start iteration
C-----------------------------------------------------------------------
C
         CASE ('RDLM      ')
C
            W_WEISS_T(:) = OBS_T(0,IWFD,:)
C
C-----------------------------------------------------------------------
C          use regular temperature mesh with linear M(T) = M - DMAGDT*T
C-----------------------------------------------------------------------
C
C        CASE ('MLIN      ')
C
         END SELECT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
         CALL THERMAL_FLUCT_MESH_WEIGHT(I_TEMP_LAT,TEMP_LAT,MM_MC,
     &                                  NT_MAGNETIC,IT_MAG,W_WEISS_T)
C
      ELSE
C
         IFT = 0
         DO IT = 1,NT
            IFT = IFT + 1
C
            FTET_FT(IFT) = 0D0
            FPHI_FT(IFT) = 0D0
C
         END DO
C
      END IF
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
C=======================================================================
C
99001 FORMAT (//,1X,79('*'),/,33X,'<THERMAL_INIT>',/,1X,79('*'),//,5X,
     &        'dealing with thermal spin fluctuations using',/,5X,
     &        'FLUCT_DIR_SETTING = ',A,/)
99002 FORMAT (//,1X,79('*'),/,33X,'<THERMAL_INIT>',/,1X,79('*'),//,5X,A,
     &        /)
99003 FORMAT (I8,' temperatures in the range  T = ',F8.2,' -- ',F8.2,
     &        ' [K]',/)
99004 FORMAT (5X,'No. of displacements:',6X,I6,/,5X,
     &        'Debye temperature:     ',F10.2,' [K]',/)
99005 FORMAT (5X,'setting of parameters and array dimensions',//,5X,
     &        'NVIBFLU  =',I6,6X,'NT      =',I6,/,5X,'NVT      =',I6,6X,
     &        'NVTMAX  =',I6,/,5X,'NFT      =',I6,6X,'NFTMAX  =',I6,/,
     &        5X,'NVFT     =',I6,6X,'NVFTMAX =',I6,/)
99006 FORMAT (//,1X,79('*'),/,33X,'<THERMAL_INIT>',/,1X,79('*'),//,5X,
     &        'dealing with thermal  vibrations  for:',//,5X,'ITMP =',
     &        I4,3X,'T = ',F8.3,' [K]',3X,'NVIBRA =',I3,/)
99007 FORMAT (10X,'IT  vib. dir.',8X,'RMS-distortion in  [a.u.] ')
99008 FORMAT (//,1X,79('*'),/,33X,'<THERMAL_INIT>',/,1X,79('*'),//,5X,
     &        'dealing with thermal fluctuations for:',//,5X,'ITMP =',
     &        I4,3X,'T = ',F8.3,' [K]',3X,15X,'          NFLUCT =',I5,/)
99009 FORMAT (I12,2X,I5,4X,3F12.5)
99010 FORMAT (/,5X,'initial values for the Weiss field:',/)
99011 FORMAT (5X,'IT =',I4,5X,'W_Weiss =',F10.4,' a.u.')
99012 FORMAT (//,1X,79('*'),/,33X,'<THERMAL_INIT>',/,1X,79('*'),//,10X,
     &        'this is the first program run with  TETDEPPOT  set',/,
     &        10X,'the start potentials for NT =',I5,'    pseudo types',
     &        /,10X,'have been written to   ',A,/,10X,
     &        'restart program using new potential file',/)
99013 FORMAT (10X,'FLUCT_DIR_SETTING: ',A,/,10X,
     &        'allowed input: SETFLUCT = M_T, MCS, RDLM, MLIN')
99014 FORMAT (10X,'N_TEMP_LAT:        ',I10,/,10X,'TEMP_LAT_MIN:      ',
     &        F10.3,/,10X,'TEMP_LAT_MAX:      ',F10.3,/)
      END
C*==thermal_indices.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INDICES(IVFT,IT,IFLUCT,IVIBRA,IVT,IFT)
C   ********************************************************************
C   *                                                                  *
C   *    reconstract from the combined index   IVFT                    *
C   *                                                                  *
C   *    IVFT = (NVIBFLU)*(IT-1) + NFLUCT*(IFLUCT+1) + IVIBRA          *
C   *                                                                  *
C   *    all other indices     IT,IFLUCT,IVIBRA,IVT,IFT                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_THERMAL,ONLY:NVIBFLU,NVIBRA,NFLUCT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFLUCT,IFT,IT,IVFT,IVIBRA,IVT
C
C*** End of declarations rewritten by SPAG
C
      IT = (IVFT/NVIBFLU) + 1
      IFLUCT = ((IVFT-NVIBFLU*(IT-1))/NFLUCT) + 1
      IVIBRA = IVFT - NVIBFLU*(IT-1) - NFLUCT*(NFLUCT-1)
C
      IVT = NVIBRA*(IT-1) + IVIBRA
      IFT = NFLUCT*(IT-1) + IFLUCT
C
      IF ( IVFT.NE.NVIBFLU*(IT-1)+NFLUCT*(NFLUCT-1)+IVIBRA .OR. 
     &     IVIBRA.GT.NVIBRA .OR. IFLUCT.GT.NFLUCT ) THEN
         WRITE (6,99001) IVFT,IT,IFLUCT,IVIBRA
         STOP
      END IF
C
99001 FORMAT (//,1X,79('#'),/,31X,'<THERMAL_INDICES>',/,1X,79('#'),//,
     &        10X,'indices connected with thermal effects inconsistent',
     &        /,10X,'input   IVFT    ',I10,/,10X,'        IT      ',I10,
     &        /,10X,'        IFLUCT  ',I10,/,10X,'        IVIBRA  ',I10,
     &        /,10X,'the program stops  ',/,1X,79('#'))
      END
C*==thermal_fluct_gen_mesh.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_FLUCT_GEN_MESH
C   ********************************************************************
C   *                                                                  *
C   *  NTET_FLUCT  -  number of THETA angles                           *
C   *  NPHI_FLUCT  -  number of PHI   angles                           *
C   *  NFLUCT -  number of all M directions                            *
C   *  TET_FLUCT(IFLUCT) - THETA angle for M direction IFLUCT in DEGs  *
C   *  PHI_FLUCT(IFLUCT) - PHI   angle for M direction IFLUCT in DEGs  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI,DEG_ARC
      USE MOD_TYPES,ONLY:NT
      USE MOD_THERMAL,ONLY:NFLUCT,TET_FLUCT,PHI_FLUCT,NPHI_FLUCT,
     &    NTET_FLUCT,EHAT_FLUCT,NTET_FLUCT_POT
      USE MOD_SCF,ONLY:SCF_THETA_DEPENDENT_POT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 COS_TET,DPHI,DTET,EHAT_X,EHAT_Y,EHAT_Z,SIN_TET,TET_ARC
      INTEGER IFLUCT,IPHI,IT,ITET,IT_CHEM,NT_CHEM
C
C*** End of declarations rewritten by SPAG
C
      IF ( NTET_FLUCT.GT.2 ) THEN
C
         DTET = 180.D0/(NTET_FLUCT-1)
         DPHI = 360.D0/(NPHI_FLUCT)
C
         IFLUCT = 0
         DO ITET = 1,NTET_FLUCT
            DO IPHI = 1,NPHI_FLUCT
C
               IFLUCT = IFLUCT + 1
               TET_FLUCT(IFLUCT) = (ITET-1)*DTET
               PHI_FLUCT(IFLUCT) = (IPHI-1)*DPHI
C
            END DO
         END DO
C
C----------------------------------------------------------------------
C          for testing purposes set ALL directions the same
C----------------------------------------------------------------------
      ELSE IF ( NTET_FLUCT.EQ.2 ) THEN
C
C---------------------------------------- z-direction
C        TET_FLUCT(1:NFLUCT) = 0.0D0
C        PHI_FLUCT(1:NFLUCT) = 0.01D0
C
C---------------------------------------- x-direction
         TET_FLUCT(1:NFLUCT) = 90.0D0
         PHI_FLUCT(1:NFLUCT) = 0.0D0
C
C---------------------------------------- y-direction
C        TET_FLUCT(1:NFLUCT) = 90.0D0
C        PHI_FLUCT(1:NFLUCT) = 90.0D0
C
      ELSE IF ( NTET_FLUCT.EQ.1 ) THEN
C
C-----------------------------------------------------------------------
C        theta-dependent potential with TET_FLUCT(IFLUCT) = 0.0D0
C   OR   conical distribution of moments
C        the common value for THETA will be changed according to M(T)
C-----------------------------------------------------------------------
C
         DPHI = 360.D0/NPHI_FLUCT
         IFLUCT = 0
C
         DO IPHI = 1,NPHI_FLUCT
C
            IFLUCT = IFLUCT + 1
            TET_FLUCT(IFLUCT) = 0.0D0
            PHI_FLUCT(IFLUCT) = (IPHI-1)*DPHI
C
         END DO
C
      END IF
C
      IF ( SCF_THETA_DEPENDENT_POT ) THEN
C
         NT_CHEM = NT/NTET_FLUCT_POT
         DTET = PI/(NTET_FLUCT_POT-1)
C
         DO IT_CHEM = 1,NT_CHEM
            DO ITET = 1,NTET_FLUCT_POT
               IT = (IT_CHEM-1)*NT_CHEM + ITET
               TET_ARC = (ITET-1)*DTET
               SIN_TET = SIN(TET_ARC)
               COS_TET = COS(TET_ARC)
C
               DO IFLUCT = 1,NFLUCT
                  EHAT_FLUCT(1,IFLUCT,IT)
     &               = COS(DEG_ARC*PHI_FLUCT(IFLUCT))*SIN_TET
                  EHAT_FLUCT(2,IFLUCT,IT)
     &               = SIN(DEG_ARC*PHI_FLUCT(IFLUCT))*SIN_TET
                  EHAT_FLUCT(3,IFLUCT,IT) = COS_TET
               END DO
            END DO
         END DO
C
      ELSE
C
         DO IFLUCT = 1,NFLUCT
C
            SIN_TET = SIN(DEG_ARC*TET_FLUCT(IFLUCT))
C
            EHAT_X = COS(DEG_ARC*PHI_FLUCT(IFLUCT))*SIN_TET
            EHAT_Y = SIN(DEG_ARC*PHI_FLUCT(IFLUCT))*SIN_TET
            EHAT_Z = COS(DEG_ARC*TET_FLUCT(IFLUCT))
C
            DO IT = 1,NT
               EHAT_FLUCT(1,IFLUCT,IT) = EHAT_X
               EHAT_FLUCT(2,IFLUCT,IT) = EHAT_Y
               EHAT_FLUCT(3,IFLUCT,IT) = EHAT_Z
            END DO
C
         END DO
C
      END IF
C
      END
C*==thermal_init_fluct_vibra_start.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INIT_FLUCT_VIBRA_START
C   ********************************************************************
C   *                                                                  *
C   *  set up  NFLUCT  directions  FTET_FT,PHI_FT  for fluctuations    *
C   *          and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ       *
C   *                                                                  *
C   *  NVIBRA NFLUCT set in input                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_THERMAL,ONLY:NVFT,NVFTMAX,NVT,NVTMAX,NFT,NFTMAX,X_VFT,
     &    X_FT,X_VT,NVFO_Q,IQ_AVFT,IVFT_VFOQ,NVIBFLU,NVIBRA,NAVFT,
     &    NFLUCT,NFLUCT_T
      USE MOD_CPA,ONLY:NCPA
      USE MOD_TYPES,ONLY:NT,CONC
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_FLUCT_VIBRA_START')
C
C Local variables
C
      INTEGER IFLUCT,IFT,IO,IQ,IT,IVFO,IVFT,IVIBFLU,IVIBRA,IVT
      REAL*8 WVIBRA
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      NVIBFLU = NVIBRA*NFLUCT
C
      NFLUCT_T(:) = 1
      NVFT = NVIBFLU*NT
      NFT = NFLUCT*NT
C
      NVFTMAX = NVFT
      NVTMAX = NVT
      NFTMAX = NFT
C
      ALLOCATE (X_VT(NVTMAX))
      ALLOCATE (X_VFT(NVFTMAX),NVFO_Q(NQMAX))
      ALLOCATE (IVFT_VFOQ(NVFTMAX,NQMAX),IQ_AVFT(NQMAX,NVFTMAX))
      ALLOCATE (NAVFT(NVFTMAX))
C
C-----------------------------------------------------------------------
      DO IQ = 1,NQ
         NVFO_Q(IQ) = NVIBFLU*NOQ(IQ)
         IVFO = 0
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            IVFT = (IT-1)*NVIBFLU
            DO IVIBFLU = 1,NVIBFLU
               IVFT = IVFT + 1
               IVFO = IVFO + 1
               IVFT_VFOQ(IVFO,IQ) = IVFT
            END DO
         END DO
      END DO
C
      NAVFT(1:NVFTMAX) = 0
      NCPA = 0
      DO IQ = 1,NQ
         DO IVFO = 1,NVFO_Q(IQ)
            IVFT = IVFT_VFOQ(IVFO,IQ)
            NAVFT(IVFT) = NAVFT(IVFT) + 1
            IQ_AVFT(NAVFT(IVFT),IVFT) = IQ
         END DO
         IF ( NVFO_Q(IQ).GT.1 ) THEN
            ICPA(IQ) = 1
            NCPA = NCPA + ICPA(IQ)
         END IF
      END DO
C
      WVIBRA = 1D0/DBLE(NVIBRA)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IVT = (IT-1)*NVIBRA
         DO IVIBRA = 1,NVIBRA
            IVT = IVT + 1
            X_VT(IVT) = WVIBRA
         END DO
C
         IFT = (IT-1)*NFLUCT
         DO IFLUCT = 1,NFLUCT
            IFT = IFT + 1
C
            IVT = (IT-1)*NVIBRA
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)*NFLUCT + IFLUCT
C
               X_VFT(IVFT) = X_VT(IVT)*X_FT(IFT)*CONC(IT)
C
            END DO
         END DO
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==thermal_init_fluct_m_t.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INIT_FLUCT_M_T(TCURIE,N_TEMP_LAT,NT_MAGNETIC,
     &                                  MM_T0)
C   ********************************************************************
C   *                                                                  *
C   * IT_MAG(IT) = 1 when magnetic moment M(IT) > 0                    *
C   * IT_MAG(IT) = 0 when magnetic moment M(IT) = 0                    *
C   * NFLUCT_T(IT) - number of fluctuations applied to atom IT         *
C   * NVIBFLU_T(IT)- number of (fluctuations and vibrations)           *
C   *                applied to atom IT                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NT,NTMAX
      USE MOD_FILES,ONLY:IFILFLUCT,FLUCTFIL,RDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_FLUCT_M_T')
C
C Dummy arguments
C
      INTEGER NT_MAGNETIC,N_TEMP_LAT
      REAL*8 TCURIE
      REAL*8 MM_T0(NTMAX)
C
C Local variables
C
      CHARACTER*40 CH_DUMMY
      INTEGER I,IT,IT_MAG(NT),I_DUMMY
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C ---------------------------------- Read M(T) information from data file
C
      OPEN (UNIT=IFILFLUCT,FILE=FLUCTFIL)
      READ (IFILFLUCT,*) NT_MAGNETIC
      READ (IFILFLUCT,*)
      READ (IFILFLUCT,*)
      DO IT = 1,NT
         READ (IFILFLUCT,*) I_DUMMY,CH_DUMMY,IT_MAG(IT)
         WRITE (6,*) I_DUMMY,CH_DUMMY,IT_MAG(IT)
      END DO
      READ (IFILFLUCT,*)
      READ (IFILFLUCT,*) N_TEMP_LAT
      READ (IFILFLUCT,*) TCURIE
C
      READ (IFILFLUCT,*) RDUMMY,(MM_T0(I),I=1,NT_MAGNETIC)
C
      WRITE (6,99001) FLUCTFIL,NT_MAGNETIC,N_TEMP_LAT,TCURIE,
     &                (MM_T0(I),I=1,NT_MAGNETIC)
C
C=======================================================================
99001 FORMAT (5X,'FLUCT_DIR_SETTING:  M_T    ',/,5X,
     &        'data file  MGNTFIL: ',A,/,5X,
     &        'No. of magnetic atoms:     ',I6,/,5X,
     &        'No. of temperature values: ',I6,/,5X,
     &        'Curie temperature:     ',F10.2,' [K]',/,5X,
     &        'Magnetic moments:      ',10F6.2,/)
C
      END
C*==thermal_init_fluct_mcs.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INIT_FLUCT_MCS(TEMP_LAT_MIN,TEMP_LAT_MAX,
     &                                  N_TEMP_LAT,NFLUCT)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILFLUCT,FLUCTFIL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_FLUCT_MCS')
C
C Dummy arguments
C
      INTEGER NFLUCT,N_TEMP_LAT
      REAL*8 TEMP_LAT_MAX,TEMP_LAT_MIN
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C-----------------------------------------------------------------------
C   open the file  FLUCTFIL with NFLUCT spin configurations generated
C   within the MC calculations for certain temperature
C
      OPEN (IFILFLUCT,FILE=FLUCTFIL)
      READ (IFILFLUCT,*) N_TEMP_LAT
      READ (IFILFLUCT,*) TEMP_LAT_MIN
      READ (IFILFLUCT,*) TEMP_LAT_MAX
      READ (IFILFLUCT,*) NFLUCT
C
      END
C*==thermal_init_fluct_rdlm.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_INIT_FLUCT_RDLM(TEMP_LAT_MIN,TEMP_LAT_MAX,
     &   N_TEMP_LAT,NFLUCT)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILFLUCT,FLUCTFIL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_FLUCT_MCS')
C
C Dummy arguments
C
      INTEGER NFLUCT,N_TEMP_LAT
      REAL*8 TEMP_LAT_MAX,TEMP_LAT_MIN
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C-----------------------------------------------------------------------
C   open the file  FLUCTFIL with NFLUCT spin configurations generated
C   within the MC calculations for certain temperature
C
      OPEN (IFILFLUCT,FILE=FLUCTFIL)
      READ (IFILFLUCT,*) N_TEMP_LAT
      READ (IFILFLUCT,*) TEMP_LAT_MIN
      READ (IFILFLUCT,*) TEMP_LAT_MAX
      READ (IFILFLUCT,*) NFLUCT
C
      END
C*==thermal_weiss_field.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_WEISS_FIELD(TCURIE,TEMP_LAT,MM_MC,NT_MAGNETIC,
     &                               IT_MAG,W_WEISS_T)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:KB_SI,EV_J,RY_EV
      USE MOD_TYPES,ONLY:NT,NTMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_WEISS_FIELD')
C
C Dummy arguments
C
      INTEGER NT_MAGNETIC
      REAL*8 TCURIE,TEMP_LAT
      INTEGER IT_MAG(NT)
      REAL*8 MM_MC(NT_MAGNETIC),W_WEISS_T(NTMAX)
C
C Local variables
C
      REAL*8 DELTA,DELTA0,DWEISS,LANGEVIN,MAGN_AV,SGN,X_DUM
      INTEGER IFIT,IT,IT_MC,IWEISS
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IT_MC = 0
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IT_MC = IT_MC + IT_MAG(IT)
C --------------------------------------------  Magnetization curve from
C------------------------------------------------  solution of Eq: x=L(x)
C
         MAGN_AV = 1.D0
         W_WEISS_T(IT) = 3.D0*TCURIE
         DWEISS = W_WEISS_T(IT)/1000
         DO IFIT = 1,1000
            X_DUM = MAGN_AV*W_WEISS_T(IT)/TEMP_LAT
            LANGEVIN = 1.D0/TANH(X_DUM) - 1.D0/X_DUM
            IF ( DABS(MAGN_AV-LANGEVIN).LE.1.D-6 ) EXIT
            MAGN_AV = LANGEVIN
C
         END DO
C
         IF ( IPRINT.GT.1 ) WRITE (5500,*) TEMP_LAT,MAGN_AV
C
C----------------------------- Magnetization fitting to MC or experiment
C -------------------------------------------------- to find Weiss field
C
         SGN = 1.D0
         W_WEISS_T(IT) = 3.D0*TCURIE
         DELTA0 = 10.D0
         DO IWEISS = 1,100000
C
            MAGN_AV = 1.D0
            DO IFIT = 1,10000
               X_DUM = MAGN_AV*W_WEISS_T(IT)/TEMP_LAT
               LANGEVIN = 1.D0/TANH(X_DUM) - 1.D0/X_DUM
               IF ( DABS(MAGN_AV-LANGEVIN).LE.1.D-11 ) EXIT
               MAGN_AV = LANGEVIN
               IF ( DABS(MAGN_AV).LE.1.D-8 ) MAGN_AV = 1.D-8
            END DO
C
            DELTA = DABS(MAGN_AV-MM_MC(IT_MC))
            IF ( DABS(DELTA).LE.1.D-11 ) EXIT
C
            IF ( (DELTA0-DELTA).LE.0 ) THEN
               DWEISS = DWEISS/2
               SGN = -SGN
            END IF
C
            DELTA0 = DELTA
            W_WEISS_T(IT) = W_WEISS_T(IT) + SGN*DWEISS
C
         END DO
C
         W_WEISS_T(IT) = MAGN_AV*W_WEISS_T(IT)*KB_SI/EV_J/RY_EV
C
C ----------------------------------------------------------------------
         IF ( IPRINT.GT.1 ) THEN
            WRITE (5501+IT_MC*10,*) TEMP_LAT,W_WEISS_T(IT),3.D0*TCURIE
            WRITE (5502+IT_MC*10,*) TEMP_LAT,MAGN_AV,MM_MC(IT_MC)
         END IF
C ----------------------------------------------------------------------
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==thermal_m_rdlm.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_M_RDLM(TCURIE,TEMP_LAT,MM_MC,NT_MAGNETIC,
     &                          IT_MAG,W_WEISS_T)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:KB_SI,EV_J,RY_EV
      USE MOD_TYPES,ONLY:NT,NTMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_M_RDLM')
C
C Dummy arguments
C
      INTEGER NT_MAGNETIC
      REAL*8 TCURIE,TEMP_LAT
      INTEGER IT_MAG(NT)
      REAL*8 MM_MC(NT_MAGNETIC),W_WEISS_T(NTMAX)
C
C Local variables
C
      INTEGER IFIT,IT,IT_MC
      REAL*8 LANGEVIN,MAGN_AV,X_DUM
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IT_MC = 0
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IT_MC = IT_MC + IT_MAG(IT)
C
C --------------------------------------------  Magnetization curve from
C------------------------------------------------  solution of Eq: x=L(x)
C
         MAGN_AV = 1.D0
         W_WEISS_T(IT) = 3.D0*TCURIE
         DO IFIT = 1,1000
            X_DUM = MAGN_AV*W_WEISS_T(IT)/TEMP_LAT
            LANGEVIN = 1.D0/TANH(X_DUM) - 1.D0/X_DUM
            IF ( DABS(MAGN_AV-LANGEVIN).LE.1.D-6 ) EXIT
            MAGN_AV = LANGEVIN
C
         END DO
         MM_MC(IT) = MAGN_AV
C
      END DO
      W_WEISS_T(IT) = W_WEISS_T(IT)*KB_SI/EV_J/RY_EV
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==thermal_cont.f    processed by SPAG 6.70Rc at 17:01 on  7 Mar 2017
      SUBROUTINE THERMAL_CONT
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NT,NTMAX,W_WEISS_T
      USE MOD_THERMAL,ONLY:TEMP_LAT,NFLUCT_T,NVIBFLU_T
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT')
C
C Local variables
C
      INTEGER IT,IT_MAG(NT),I_TEMP_LAT,NT_MAGNETIC
      REAL*8 MM_MC(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MM_MC
C
      CALL TRACK_INFO(ROUTINE)
C
C
      NT_MAGNETIC = NT
      DO IT = 1,NT
         IT_MAG(IT) = 1
      END DO
C
      IF ( .NOT.ALLOCATED(NFLUCT_T) ) ALLOCATE (NFLUCT_T(NTMAX))
      IF ( .NOT.ALLOCATED(NVIBFLU_T) ) ALLOCATE (NVIBFLU_T(NT))
      IF ( .NOT.ALLOCATED(MM_MC) ) ALLOCATE (MM_MC(NTMAX))
      IF ( .NOT.ALLOCATED(W_WEISS_T) ) ALLOCATE (W_WEISS_T(NTMAX))
C
      MM_MC(:) = 1.D0
      I_TEMP_LAT = 1
C
      WRITE (6,'(A,2F20.10)') '     TEMPERATURE = ',TEMP_LAT
      DO IT = 1,NT
         WRITE (6,'(A5,I5,A,2F20.10)') 'IT = ',IT,'   WEISS FIELD = ',
     &                                 W_WEISS_T(IT)
      END DO
C
C-----------------------------------------------------------------------
C
      CALL THERMAL_FLUCT_MESH_WEIGHT(I_TEMP_LAT,TEMP_LAT,MM_MC,
     &                               NT_MAGNETIC,IT_MAG,W_WEISS_T)
C
      END
