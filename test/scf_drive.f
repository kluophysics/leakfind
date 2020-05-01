C*==scf_drive.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCF_DRIVE(KBZI,KMOL,FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *  driving routine for SCF calculations                            *
C   *  that allows to vary specific systems parameters                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE,ALAT,ABAS
      USE MOD_SCF,ONLY:SCFSTATUS,ALAT_TAB,EFERMI_TAB,ETOT_TAB,
     &    MUEORB_TAB,MUESPN_TAB,VOL_TAB
      USE MOD_ENERGY,ONLY:EMIN,EMAX,EFERMI
      USE MOD_FILES,ONLY:IPRINT,POTFIL,LPOTFIL,IFILSCLALAT,DATSET,
     &    LDATSET,SYSTEM,LSYSTEM,WRPOT,FOUND_SECTION,FOUND_REAL
      USE MOD_TYPES,ONLY:ETOT,MUEORB,MUESPN
      USE MOD_CALCMODE,ONLY:KKRMODE,SCL_ALAT,IREL
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_MPI,ONLY:MPI_ID
      USE MOD_THERMAL,ONLY:X_VFT,UMAT_VT,NVIBRA,NFLUCT,FMAT_FT
c modified by XJQ: scf of vibrations
      use mod_scfvb_cpa_sigma,only:lscfvb,scfvb_cpa_sigma
c end-mod-xjq
      IMPLICIT NONE
C*--SCF_DRIVE21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF_DRIVE')
C
C Dummy arguments
C
      LOGICAL FULLPOT5
      INTEGER KBZI,KMOL
C
C Local variables
C
      LOGICAL DELTA_SCHEME,QAVAILABLE,UDT1
      COMPLEX*16 ERYD
      REAL*8 ETOP,PSCLMAX_VOL,SCL,SCL1,SCL2,SCL_VOL
      INTEGER ISCL,ISCL0,LPOTFIL0,NSCL
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      CALL INPUT_FIND_SECTION('SCF',0)
C
      IF ( FOUND_SECTION ) THEN
C
C=======================================================================
C scale lattice parameter ALAT and all dependent radial mesh parameters
C=======================================================================
C
         CALL SECTION_FIND_KEYWORD('SCL_ALAT',SCL_ALAT)
C
         IF ( SCL_ALAT ) THEN
C
            IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' .OR. 
     &           SYSTEM_DIMENSION(1:2).NE.'3D' ) THEN
               WRITE (6,99005)
               CALL STOP_MESSAGE(ROUTINE,
     &                         'SYSTEM_TYPE or SYSTEM_DIMENSION  not OK'
     &                         )
            END IF
C
            CALL SECTION_FIND_KEYWORD('DELTA',DELTA_SCHEME)
C
C------------------------------------------------------------------------
C                scale the volume in NSCL equidistant steps
C------------------------------------------------------------------------
            IF ( DELTA_SCHEME ) THEN
C
               NSCL = 7
               PSCLMAX_VOL = 6D0
               IF ( MOD(NSCL,2).EQ.0 ) NSCL = NSCL + 1
               IF ( NSCL.LE.3 ) CALL STOP_MESSAGE(ROUTINE,'NSCL>=3')
               ISCL0 = (NSCL+1)/2
C
            ELSE
C
C------------------------------------------------------------------------
C           scale the lattice parameter in NSCL equidistant steps
C------------------------------------------------------------------------
               CALL SECTION_SET_REAL('SCL1',SCL1,9999D0,0)
               UDT1 = FOUND_REAL
               CALL SECTION_SET_REAL('SCL2',SCL2,9999D0,0)
               CALL SECTION_SET_INTEGER('NSCL',NSCL,9999,1)
C
               IF ( .NOT.(UDT1 .AND. FOUND_REAL) )
     &              CALL STOP_MESSAGE(ROUTINE,
     &              'specify SCL1, SCL2, NSCL !!')
C
            END IF
C
            ALLOCATE (ALAT_TAB(NSCL),EFERMI_TAB(NSCL),ETOT_TAB(NSCL))
            ALLOCATE (VOL_TAB(NSCL),MUEORB_TAB(NSCL),MUESPN_TAB(NSCL))
C
            CALL SECTION_FIND_KEYWORD('WRPOT',WRPOT)
            LPOTFIL0 = LPOTFIL - 4
C
            OPEN (IFILSCLALAT,FILE=DATSET(1:LDATSET)//'_ALAT.dat')
C
            WRITE (IFILSCLALAT,99001) DATSET(1:LDATSET),
     &                                SYSTEM(1:LSYSTEM),IREL,ABAS,ALAT,
     &                                SCL1,SCL2,NSCL
            CALL FLUSH(IFILSCLALAT)
C
C=======================================================================
            DO ISCL = 1,NSCL
C
               POTFIL = POTFIL(1:LPOTFIL0)//'_ALAT_'
               CALL STRING_ADD_N(POTFIL,ISCL)
               LPOTFIL = LEN_TRIM(POTFIL)
               POTFIL = POTFIL(1:LPOTFIL)//'.pot'
               LPOTFIL = LPOTFIL + 4
C
C-----------------------------------------------------------------------
C                scale the volume in NSCL equidistant steps
C------------------------------------------------------------------------
C
               IF ( DELTA_SCHEME ) THEN
C
                  SCL_VOL = 1D0 + PSCLMAX_VOL*(ISCL-ISCL0)
     &                      /DBLE(NSCL-ISCL0)
C
                  SCL = SCL_VOL**(1.0/3.0)
C
C-----------------------------------------------------------------------
C           scale the lattice parameter in NSCL equidistant steps
C------------------------------------------------------------------------
C
               ELSE IF ( NSCL.GT.1 ) THEN
                  SCL = SCL1 + (SCL2-SCL1)*(ISCL-1)/DBLE(NSCL-1)
               ELSE
                  SCL = SCL1
               END IF
C-----------------------------------------------------------------------
C
C-------------------------------------------- enforce start from scratch
C
               IF ( FULLPOT .OR. FULLPOT5 ) THEN
                  FULLPOT5 = .TRUE.
                  FULLPOT = .FALSE.
               ELSE
                  FULLPOT5 = .FALSE.
                  FULLPOT = .FALSE.
               END IF
C
               SCFSTATUS = 'START'
C-----------------------------------------------------------------------
C
               CALL SCFSCALE_ALAT(SCL,QAVAILABLE)
C
               WRITE (6,99002) ISCL,SCL,ALAT
               IF ( WRPOT ) THEN
                  WRITE (6,99003) POTFIL(1:LPOTFIL)
               ELSE
                  WRITE (6,99004)
               END IF
C
               IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
                  ETOP = MAX(ABS(EMIN),ABS(EMAX))
C
                  CALL STRINIT(ETOP)
C
                  CALL STRCC(ERYD,.TRUE.)
C
               ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C
                  CALL INIT_MOD_TBCLU(IPRINT)
C
                  CALL INIT_TBCALC(IPRINT)
C
               ELSE
C
                  CALL STOP_MESSAGE(ROUTINE,'KKRMODE not found')
C
               END IF
C
               CALL SCF(KBZI,KMOL,FULLPOT5)
C
               ALAT_TAB(ISCL) = ALAT
               ETOT_TAB(ISCL) = ETOT
               EFERMI_TAB(ISCL) = EFERMI
               MUEORB_TAB(ISCL) = MUEORB
               MUESPN_TAB(ISCL) = MUESPN
C
C
            END DO
C=======================================================================
C
            IF ( MPI_ID.EQ.0 ) CALL SCFSCALE_ALAT_EOS(NSCL)
C
            WRITE (6,*) 
     &               'loop over scaled lattice parameter ALAT completed'
            STOP
C
         END IF
C=======================================================================
C
      END IF
C
C=======================================================================
C=======================================================================
C      standard run for a single setting of the system parameters
C=======================================================================
C=======================================================================
C
      POTFIL = POTFIL(1:LPOTFIL)//'_new'
      LPOTFIL = LPOTFIL + 4
C
c modified by XJQ: scfvb
      if(lscfvb) then
c
        call scfvb_cpa_sigma(kbzi,kmol,fullpot5)
c
      else
c
        CALL SCF(KBZI,KMOL,FULLPOT5)
c
      endif
c end-mod-xjq
C
C=======================================================================
99001 FORMAT ('#',/,'#  SCF-runs for ',/,'#  DATSET = ',A,/,
     &        '#  SYSTEM = ',A,/,'#  IREL   =',I2,/,'#  ABAS   = ',
     &        3F14.8,/,2('#',11X,3F14.8,/),
     &        '#  scaling the lattice parameter ALAT = ',F14.8,/,
     &        '#  using   SCL1 =',F6.3,'   SCL2 =',F6.3,'   NSCL =',I3,
     &        /,'#',/,
     &        '#     ALAT    EFERMI    TOTDOS    MUESPN    MUEORB',
     &        '                ETOT    SCLNOS      RMSV      RMSB')
99002 FORMAT (//,3(1X,79('*'),/),/,34X,'<SCF_DRIVE>',//,3(1X,79('*'),/),
     &        /,10X,'scaling the lattice parameter ALAT ',//,10X,'run',
     &        I3,':',6X,'SCL =',F6.3,6X,'ALAT = ',F10.6,/)
99003 FORMAT (10X,'potential will be written to ',A,/)
99004 FORMAT (10X,'NO potential will be written ',/)
99005 FORMAT (//,10X,'scaling of lattice parameter ALAT',/10X,
     &        'only for 3D bulk systems',/)
      END
