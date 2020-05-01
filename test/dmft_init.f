C*==dmft_init.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_INIT(NE,ETAB,DMFTMIX,DMFTDBLC,NEMAX,IPRINT)
C=======================================================================
C                Initialisation of DMFT (common for SCF, GEN, SPEC)
C                         DMFT-calculational mode
C
C     Readin input file
C     Allocate DMFTSIGMA and READ IT IN FROM FILE
C=======================================================================
C
C
      USE MOD_FILES,ONLY:FOUND_REAL,FOUND_INTEGER_ARRAY,
     &    FOUND_STRING_ARRAY
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,SUB_SYSTEM
      USE MOD_TYPES,ONLY:LOPT,NTMAX,ITBOT,ITTOP,TXT_T,Z,NT
      USE MOD_SCF,ONLY:SCFMIX
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_DMFT_LDAU,ONLY:JEFF,UEFF,KSELF,DMFTSIGMA,EREFDMFT,DMFTTEMP
      IMPLICIT NONE
C*--DMFT_INIT20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*4 DMFTDBLC
      REAL*8 DMFTMIX
      INTEGER IPRINT,NE,NEMAX
      COMPLEX*16 ETAB(NE)
C
C Local variables
C
      INTEGER IT,ITBOTSAV,ITTOPSAV,NIN
      REAL*8 RAUX
      CHARACTER*10 STR10
      CHARACTER*1 STRLOPT(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE STRLOPT
C
      ALLOCATE (STRLOPT(NTMAX))
C
      IF ( .NOT.ALLOCATED(DMFTSIGMA) ) THEN
         ALLOCATE (DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX))
         CALL CINIT(NKMMAX*NKMMAX*NTMAX*NEMAX,DMFTSIGMA)
      END IF
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      CALL SECTION_SET_REAL('DMFTMIX',DMFTMIX,SCFMIX,0)
C
      CALL SECTION_SET_STRING('DBLC',DMFTDBLC,'META',0)
C
      CALL SECTION_SET_REAL('EREF',RAUX,0.7D0,0)
      EREFDMFT = DCMPLX(RAUX,0.0D0)
C
      IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         ITBOT = 1
         ITTOP = NT
      END IF
C
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
      DO IT = ITBOT,ITTOP
         STR10 = 'UEFF'
         CALL STRING_ADD_N(STR10,IT)
         CALL SECTION_SET_REAL(STR10,UEFF(IT),9999D0,0)
C
         IF ( FOUND_REAL ) KSELF(IT) = 1
C
         STR10 = 'JEFF'
         CALL STRING_ADD_N(STR10,IT)
         CALL SECTION_SET_REAL(STR10,JEFF(IT),9999D0,0)
      END DO
C
      CALL SECTION_SET_STRING_ARRAY('LOPT',STRLOPT,NIN,NTMAX,0,'9999',0)
C
      IF ( .NOT.FOUND_STRING_ARRAY .OR. NIN.NE.NT ) THEN
         DO IT = ITBOT,ITTOP
            IF ( 22.LE.Z(IT) .AND. Z(IT).LE.28 ) LOPT(IT) = 2
            IF ( 40.LE.Z(IT) .AND. Z(IT).LE.46 ) LOPT(IT) = 2
            IF ( 72.LE.Z(IT) .AND. Z(IT).LE.78 ) LOPT(IT) = 2
            IF ( 58.LE.Z(IT) .AND. Z(IT).LE.71 ) LOPT(IT) = 3
            IF ( 89.LE.Z(IT) ) LOPT(IT) = 3
         END DO
      ELSE
         DO IT = ITBOT,ITTOP
            IF ( STRLOPT(IT).EQ.'s' .OR. STRLOPT(IT).EQ.'0' ) THEN
               LOPT(IT) = 0
            ELSE IF ( STRLOPT(IT).EQ.'p' .OR. STRLOPT(IT).EQ.'1' ) THEN
               LOPT(IT) = 1
            ELSE IF ( STRLOPT(IT).EQ.'d' .OR. STRLOPT(IT).EQ.'2' ) THEN
               LOPT(IT) = 2
            ELSE IF ( STRLOPT(IT).EQ.'f' .OR. STRLOPT(IT).EQ.'3' ) THEN
               LOPT(IT) = 3
            ELSE
               LOPT(IT) = -1
            END IF
         END DO
      END IF
      DO IT = ITBOT,ITTOP
         IF ( LOPT(IT).LT.0 ) KSELF(IT) = 0
      END DO
C
      WRITE (6,99004)
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).EQ.1 .AND. LOPT(IT).GT.0 ) WRITE (6,99003) IT,
     &        TXT_T(IT),LOPT(IT),UEFF(IT),JEFF(IT)
      END DO
C
      CALL DMFT_READSIG(NE,ETAB,IPRINT)
C
      IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
      END IF
C
      WRITE (6,'(/)')
      WRITE (6,99002) 'Temperature: ',DMFTTEMP
      WRITE (6,99002) 'EREF       : ',DREAL(EREFDMFT)
      WRITE (6,99002) 'DMFTMIX    : ',DMFTMIX
      WRITE (6,99001) 'DBL-COUNTING:       ',DMFTDBLC
      WRITE (6,'(1X,79(''*''),/)')
C
99001 FORMAT (10X,A,A)
99002 FORMAT (10X,A,5F20.10)
99003 FORMAT (10X,'IT =',I3,2X,A,'  LOP =',I2,2X,'U_eff =',F5.2,' eV ',
     &        2X,'J_eff =',F5.2,' eV ')
C
99004 FORMAT (//,1X,79('*'),/,28X,'parameters for LDA+DMFT',/,1X,79('*')
     &        )
      END
