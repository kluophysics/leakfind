C*==init_mod_dmft_ldau.f    processed by SPAG 6.70Rc at 13:38 on 14 Dec 2016
      SUBROUTINE INIT_MOD_DMFT_LDAU(NE,ETAB,DMFTMIX,NEMAX)
C
C=======================================================================
C                Initialisation of DMFT (common for SCF, GEN, SPEC)
C                         DMFT-calculational mode
C
C     Readin input file
C     Allocate DMFTSIGMA and READ IT IN FROM FILE
C=======================================================================
C
      USE MOD_FILES,ONLY:IPRINT,FOUND_REAL,FOUND_INTEGER_ARRAY,
     &    FOUND_STRING_ARRAY
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,SUB_SYSTEM
      USE MOD_TYPES,ONLY:LOPT,NTMAX,ITBOT,ITTOP,TXT_T,Z,NT,NLMFPMAX
      USE MOD_CALCMODE,ONLY:DMFT,LDAU,ORBPOL,KKRMODE
      USE MOD_DMFT_LDAU,ONLY:DMFTTEMP,DMFTSCF,DMFTSIGMA,EREFDMFT,KSELF,
     &    UEFF,JEFF,DMFTSIG,SEBT,SEVT,SEVNST,SEBNST,EREFLDAU,DMFTSOLVER,
     &    ELDAU,DMFTDBLC,SYMMETRISE_OCC,DMFT_FIX_DYN_SE,IBASIS,UMODE,
     &    SIGTOL,IEREF,zerosig_contn
      USE MOD_SCF,ONLY:SCFMIX
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_RMESH,ONLY:NRMAX,JRNSMIN
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DMFTMIX
      INTEGER NE,NEMAX
      COMPLEX*16 ETAB(NE)
C
C Local variables
C
      INTEGER IT,ITBOTSAV,ITTOPSAV,NIN
      REAL*8 RAUX
      CHARACTER*10 STR10
      CHARACTER*1 STRLOPT(:)
      LOGICAL UDT
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE STRLOPT
C
      WRITE (6,99006)
C
C-----------------------------------------------------------------------
C       Allocate and initialise variables
C-----------------------------------------------------------------------
C
      ALLOCATE (STRLOPT(NTMAX))
      IF ( .NOT.ALLOCATED(DMFTSIGMA) )
     &     ALLOCATE (DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX))
      CALL CINIT(NKMMAX*NKMMAX*NTMAX*NEMAX,DMFTSIGMA)
      IF ( .NOT.ALLOCATED(DMFTSIG) )
     &     ALLOCATE (DMFTSIG(NKMMAX,NKMMAX,NTMAX))
      CALL CINIT(NKMMAX*NKMMAX*NTMAX,DMFTSIG)
      IF ( .NOT.ALLOCATED(KSELF) ) ALLOCATE (KSELF(NTMAX))
      KSELF(1:NTMAX) = 0
      IF ( .NOT.ALLOCATED(UEFF) ) ALLOCATE (UEFF(NTMAX),JEFF(NTMAX))
      UEFF(1:NTMAX) = 0.0D0
      JEFF(1:NTMAX) = 0.0D0
      IF ( .NOT.ALLOCATED(EREFDMFT) )
     &     ALLOCATE (EREFDMFT(NTMAX),EREFLDAU(NTMAX))
      IF ( LDAU ) THEN
         IF ( .NOT.ALLOCATED(ELDAU) ) ALLOCATE (ELDAU(NTMAX))
         ELDAU(1:NTMAX) = 0.0D0
      END IF
C
      IF ( .NOT.ALLOCATED(SEVT) )
     &     ALLOCATE (SEVT(NRMAX,NTMAX),SEBT(NRMAX,NTMAX))
      IF ( .NOT.ALLOCATED(SEVNST) ) THEN
         ALLOCATE (SEVNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (SEBNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
      END IF
      SEVT = C0
      SEBT = C0
      SEVNST = C0
      SEBNST = C0
C
      DMFTTEMP = 999999D0
      DMFTSOLVER(1:10) = '          '
C
      IF ( LDAU ) DMFTSOLVER(1:4) = 'LDAU'
C
      IF ( DMFT ) THEN
         DMFTSOLVER(1:4) = 'FLEX'
         IF ( ORBPOL(1:10).EQ.'DMFT-FLEX ' ) THEN
            DMFTSOLVER(1:4) = 'FLEX'
         ELSE IF ( ORBPOL(1:8).EQ.'DMFT-TMA' ) THEN
            DMFTSOLVER(1:4) = 'TMA '
         ELSE IF ( ORBPOL(1:10).EQ.'DMFT-FLEXN' ) THEN
            DMFTSOLVER(1:8) = 'FLEX-IDM'
         ELSE IF ( ORBPOL(1:8).EQ.'DMFT-3BS' ) THEN
            DMFTSOLVER(1:3) = '3BS'
         ELSE IF ( ORBPOL(1:8).EQ.'DMFT-ED ' ) THEN
            DMFTSOLVER(1:6) = 'ED    '
         ELSE IF ( ORBPOL(1:8).EQ.'DMFT-EDN' ) THEN
            DMFTSOLVER(1:6) = 'ED-IDM'
         END IF
      END IF
C
C-----------------------------------------------------------------------
C       Read in input file
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( DMFT ) THEN
C
         DMFTSCF = .FALSE.
         CALL SECTION_FIND_KEYWORD('DMFTSCF',UDT)
         IF ( UDT ) DMFTSCF = .TRUE.
C
c         CALL SECTION_SET_REAL('DMFTMIX',DMFTMIX,SCFMIX,0)
         DMFTMIX = 1.D0
C
         IF ( DMFTSOLVER(1:4).EQ.'FLEX' .OR. DMFTSOLVER(1:3)
     &        .EQ.'QMC' .OR. DMFTSOLVER(1:6).EQ.'ED-IDM' ) THEN
            CALL SECTION_SET_REAL('TEMP',DMFTTEMP,9999D0,1)
         ELSE IF ( DMFTSOLVER(1:3).EQ.'TMA' ) THEN
            DMFTTEMP = 0.0D0
         END IF
      END IF
c
      zerosig_contn=.false.
      call section_find_keyword('ZEROSIG_CONTN',zerosig_contn)
C
      DMFT_FIX_DYN_SE = .FALSE.
C
      CALL SECTION_FIND_KEYWORD('DMFT_FIX_DYN_SE',UDT)
      IF ( UDT ) DMFT_FIX_DYN_SE = .TRUE.
C
      SIGTOL = 1D-18
C
      CALL SECTION_SET_REAL('SIGTOL',SIGTOL,9999D0,0)
      IF ( FOUND_REAL ) SIGTOL = ABS(SIGTOL)
      IF ( .NOT.FOUND_REAL .AND. DMFTSCF ) SIGTOL = 1D-5
C
      IBASIS = 1
      CALL SECTION_SET_INTEGER('BASIS',IBASIS,9999,0)
      IF ( IBASIS.GT.4 ) IBASIS = 0
C
      SYMMETRISE_OCC = .FALSE.
      CALL SECTION_FIND_KEYWORD('SYMMETRISE_OCC',UDT)
      IF ( UDT ) SYMMETRISE_OCC = .TRUE.
      IF ( KKRMODE(1:16).EQ.'EMBEDDED-CLUSTER' )
     &     SYMMETRISE_OCC = .FALSE.
C     NOT YET BECAUSE OF GLOBAL COORD SYSTEM
C      SYMMETRISE_OCC = .FALSE.
C
      CALL SECTION_SET_STRING('DBLC',DMFTDBLC,'AMF ',0)
C
      IF ( ORBPOL(1:8).EQ.'LDA+U-AL' ) THEN
         DMFTDBLC = 'AAL '
      ELSE IF ( ORBPOL(1:8).EQ.'LDA+U-MF' ) THEN
         DMFTDBLC = 'AMF '
      END IF
C
      CALL SECTION_SET_STRING('UMODE',UMODE,'ROTI',0)
C
      CALL SECTION_SET_REAL('EREF',RAUX,9999D0,0)
      IF ( FOUND_REAL ) THEN
         EREFDMFT(1:NTMAX) = DCMPLX(RAUX,0.0D0)
         EREFLDAU(1:NTMAX) = DCMPLX(RAUX,0.0D0)
      ELSE
         EREFDMFT(1:NTMAX) = DCMPLX(0.7D0,0.0D0)
         EREFLDAU(1:NTMAX) = DCMPLX(0.7D0,0.0D0)
      END IF
C
      IEREF = 0
      CALL SECTION_SET_INTEGER('IEREF',IEREF,9999,0)
C
      IF ( IEREF.GT.0 ) IBASIS = 1
C
      IF ( IEREF.LT.0 ) WRITE (6,99010) 
     &                              'EREF taken from input: not updated'
      IF ( IEREF.EQ.0 ) WRITE (6,99010) 
     &                          'EREF from center of mass of SCF states'
      IF ( IEREF.EQ.1 ) WRITE (6,99010) 
     &                          'EREF from resonance of the phase shift'
C
C
      IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         ITBOT = 1
         ITTOP = NT
      END IF
      UEFF(ITBOT:ITTOP) = 0.0D0
      JEFF(ITBOT:ITTOP) = 0.0D0
C
      CALL SECTION_SET_INTEGER_ARRAY('KSELF',KSELF,NT,NTMAX,1,9999,0)
      IF ( .NOT.FOUND_INTEGER_ARRAY ) KSELF(ITBOT:ITTOP) = 0
C
      CALL SECTION_SET_REAL('UEFF',UEFF(1),9999D0,0)
      IF ( FOUND_REAL ) THEN
         KSELF(ITBOT:ITTOP) = 1
         UEFF(ITBOT:ITTOP) = UEFF(1)
      END IF
      CALL SECTION_SET_REAL('JEFF',JEFF(1),9999D0,0)
      IF ( FOUND_REAL ) JEFF(ITBOT:ITTOP) = JEFF(1)
C
!New JMinar 17.03.20 UEFF_Z26 sets U for Z=26 atoms same for J below
      DO IT = ITBOT,ITTOP
         STR10 = 'UEFF_Z'
         CALL STRING_ADD_N(STR10,Z(IT))
         CALL SECTION_SET_REAL(STR10,UEFF(IT),9999D0,0)
         CALL SECTION_SET_REAL(STR10,RAUX,9999D0,0)
C
         IF ( FOUND_REAL ) KSELF(IT) = 1
C
         STR10 = 'JEFF_Z'
         CALL STRING_ADD_N(STR10,Z(IT))
         CALL SECTION_SET_REAL(STR10,JEFF(IT),9999D0,0)
C
      END DO


C
      DO IT = ITBOT,ITTOP
         STR10 = 'UEFF'
         CALL STRING_ADD_N(STR10,IT)
         CALL SECTION_SET_REAL(STR10,UEFF(IT),9999D0,0)
         CALL SECTION_SET_REAL(STR10,RAUX,9999D0,0)
C
         IF ( FOUND_REAL ) KSELF(IT) = 1
C
         STR10 = 'JEFF'
         CALL STRING_ADD_N(STR10,IT)
         CALL SECTION_SET_REAL(STR10,JEFF(IT),9999D0,0)
C
      END DO
!End New JMinar 17.03.20

C
      CALL SECTION_SET_STRING_ARRAY('LOPT',STRLOPT,NIN,NTMAX,0,'9999',0)
C
      IF ( .NOT.FOUND_STRING_ARRAY .OR. NIN.NE.NT ) THEN
         DO IT = ITBOT,ITTOP
            IF ( 6.LE.Z(IT) .AND. Z(IT).LE.21 ) LOPT(IT) = 2
c            IF ( 6.LE.Z(IT) .AND. Z(IT).LE.21 ) LOPT(IT) = -1
            IF ( 22.LE.Z(IT) .AND. Z(IT).LE.71 ) LOPT(IT) = 2
            IF ( 58.LE.Z(IT) .AND. Z(IT).LE.71 ) LOPT(IT) = 3
            IF ( 89.LE.Z(IT) ) LOPT(IT) = 3
         END DO
      ELSE
         DO IT = ITBOT,ITTOP
            IF ( STRLOPT(IT).EQ.'S' .OR. STRLOPT(IT).EQ.'0' ) THEN
               LOPT(IT) = 0
            ELSE IF ( STRLOPT(IT).EQ.'P' .OR. STRLOPT(IT).EQ.'1' ) THEN
               LOPT(IT) = 1
            ELSE IF ( STRLOPT(IT).EQ.'D' .OR. STRLOPT(IT).EQ.'2' ) THEN
               LOPT(IT) = 2
            ELSE IF ( STRLOPT(IT).EQ.'F' .OR. STRLOPT(IT).EQ.'3' ) THEN
               LOPT(IT) = 3
            ELSE
               LOPT(IT) = -1
               KSELF(IT) = 0
            END IF
         END DO
      END IF
      DO IT = ITBOT,ITTOP
         IF ( LOPT(IT).LT.0 ) KSELF(IT) = 0
      END DO
C
      IF ( DMFT ) THEN
         WRITE (6,99004)
      ELSE IF ( LDAU ) THEN
         WRITE (6,99005)
      END IF
C
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).EQ.1 .AND. LOPT(IT).GT.0 ) WRITE (6,99003) IT,
     &        TXT_T(IT),LOPT(IT),UEFF(IT),JEFF(IT)
      END DO
C
      CALL DMFT_READSIG(NE,ETAB,IPRINT)
C
      IF ( DMFT ) THEN
         WRITE (6,99002) 'Temperature (K)): ',DMFTTEMP
         WRITE (6,99002) 'DMFTMIX         : ',DMFTMIX
         DO IT = ITBOT,ITTOP
            IF ( KSELF(IT).EQ.1 ) WRITE (6,99008) IT,DREAL(EREFDMFT(IT))
         END DO
         WRITE (6,99001) 'DBL-COUNTING    : ',DMFTDBLC
      ELSE IF ( LDAU ) THEN
         DO IT = ITBOT,ITTOP
            IF ( KSELF(IT).EQ.1 ) WRITE (6,99008) IT,DREAL(EREFDMFT(IT))
         END DO
         WRITE (6,99001) 'DBL-COUNTING    : ',DMFTDBLC
      END IF
      IF ( UMODE.EQ.'ROTI' ) THEN
         WRITE (6,99011) 
     &                 'Fully rotational invariant(Lichtenstein et al.)'
      ELSE IF ( UMODE.EQ.'DUDA' ) THEN
         WRITE (6,99011) 'Spherically avaraged (Dudarev et al.)'
      END IF
C
      IF ( SYMMETRISE_OCC ) WRITE (6,99007)
C
      IF ( IBASIS.EQ.0 ) WRITE (6,99009)
     &                           'BT, VNST=0 and shape funct. used'
      IF ( IBASIS.EQ.1 ) WRITE (6,99009)
     &                           'BT, VNST=0 and basis within RMT'
      IF ( IBASIS.EQ.2 ) WRITE (6,99009)
     &                           'BT, VNST=0 and basis within RWS'
      IF ( IBASIS.EQ.3 ) WRITE (6,99009) 
     &                     'BT, VNST=0 and basis within RMT, full GFMAT'
      IF ( IBASIS.EQ.4 ) WRITE (6,99009) 
     &                     'BT, VNST=0 and basis within RWS, full GFMAT'
C
      IF ( SYMMETRISE_OCC ) WRITE (6,99007)
C
      IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
      END IF
C
      WRITE (6,'(1X,79(''*''),/)')
C
99001 FORMAT (10X,A,A)
99002 FORMAT (10X,A,50F20.10)
99003 FORMAT (10X,'IT =',I3,2X,A,'  LOP =',I2,2X,'U_eff =',F5.2,' eV ',
     &        2X,'J_eff =',F5.2,' eV ')
C
99004 FORMAT (//,1X,79('*'),/,28X,'parameters for LDA+DMFT',/,1X,79('*')
     &        )
99005 FORMAT (//,1X,79('*'),/,28X,'parameters for LDA+U',/,1X,79('*'))
99006 FORMAT (//,1X,79('*'),/,32X,'<INIT_MOD_DMFT_LDAU>',/,1X,79('*'),/)
99007 FORMAT (10X,'Occupation matrix will be symmetrised ')
99008 FORMAT (10X,'LDAU reference energy for IT',I3,' = ',F16.12)
99009 FORMAT (10X,'Choose of localised basis functions:  ',A)
99010 FORMAT (10X,'Calculation of reference energy    :  ',A)
99011 FORMAT (10X,'LDA+U scheme:  ',A)
      END
