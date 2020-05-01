C*==clugetlimits.f    processed by SPAG 6.70Rc at 08:54 on  8 Mar 2017
      SUBROUTINE CLUGETLIMITS(NQCLU,NMCLU)
C   ********************************************************************
C   *                                                                  *
C   *  get the number of additional atom types and meshes in a cluster *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:SHIFT_CLU
      USE MOD_FILES,ONLY:IPRINT,STRINP,POTFIL_CLU,LPOTFIL_CLU,
     &    IFILPOT_CLU,FOUND_SECTION,FOUND_STRING
      USE MOD_CALCMODE,ONLY:CLUTYPE,RELAX_CLU
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUGETLIMITS')
C
C Dummy arguments
C
      INTEGER NMCLU,NQCLU
C
C Local variables
C
      INTEGER IWRI
C
C*** End of declarations rewritten by SPAG
C
      NQCLU = 0
      NMCLU = 0
C
      CALL INPUT_FIND_SECTION('CONTROL',0)
C
      IF ( .NOT.FOUND_SECTION ) THEN
         IF ( IPRINT.GT.0 ) WRITE (6,*) 'section CONTROL not found'
         RETURN
      END IF
C
      CALL SECTION_SET_STRING('CLUTYPE',STRINP,'9999',0)
C
      IF ( .NOT.FOUND_STRING ) THEN
         IF ( IPRINT.GT.0 ) WRITE (6,*) ' NO cluster calculation'
         RETURN
      END IF
C
C--------------------------------------------------- Cluster calculation
C
      IF ( STRINP(1:8).NE.'embedded' ) THEN
         CALL STOP_MESSAGE(ROUTINE,
     &                     'CLUTYPE = embedded  expected in input file')
      ELSE
         CLUTYPE = 'embedded'
      END IF
C
      CALL SECTION_SET_STRING('POTFIL_CLU',POTFIL_CLU,'9999',1)
C
      LPOTFIL_CLU = LEN_TRIM(POTFIL_CLU)
C
      IF ( IPRINT.GE.1 ) THEN
         IWRI = 6
      ELSE
         IWRI = 0
      END IF
C
      OPEN (IFILPOT_CLU,STATUS='OLD',FILE=POTFIL_CLU(1:LPOTFIL_CLU),
     &      ERR=100)
C
      CALL READKWINT(IFILPOT_CLU,'NQCLU     ',NQCLU,0,IWRI,1)
      CALL READKWINT(IFILPOT_CLU,'NMCLU     ',NMCLU,0,IWRI,1)
C
      CLOSE (IFILPOT_CLU)
C
C-----------------------------------------------------------------------
C  * if no relaxation of the cluster site positions is done
C    NMCLU = 0; i.e. use meshes as for host atoms
C  * in case of relaxation:
C     - relaxation is allowed for the first time:
C       assume worst case, i.e.  NMCLU = NQCLU
C     - the previous run was done for relaxation: NMCLU should be fixed
C
      CALL INPUT_FIND_SECTION('SITES',0)
C
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_SET_REAL('SHIFT_CLU',SHIFT_CLU,9999D0,0)
C
         IF ( RELAX_CLU .AND. NMCLU.EQ.0 ) NMCLU = NQCLU
      END IF
C
      RETURN
C
C=======================================================================
 100  CONTINUE
      WRITE (6,99001) POTFIL_CLU(1:LPOTFIL_CLU)
      CALL STOP_MESSAGE(ROUTINE,
     &       'the requested potential file could not be opened  - check'
     &       )
C
99001 FORMAT (//,10X,'POTFIL_CLU : ',A,/)
      END
C*==clusymmetry.f    processed by SPAG 6.70Rc at 08:54 on  8 Mar 2017
      SUBROUTINE CLUSYMMETRY(NTCLU,NMCLU)
C   ********************************************************************
C   *                                                                  *
C   *  find the POINT SYMMETRY of a given cluster and fix the array    *
C   *  sizes  NTCLU  and  NMCLU that way                               *
C   *                                                                  *
C   *  for an embedded cluster the allowed symmtry operations are      *
C   *  a subset of the host symmetry operations                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NQCLU_EQCLU,IQCLU_EQCLU,NSYM,SYMACCEPTED,
     &    SYMSYMBL,SYMEULANG,ISYMGENQ,IQORGQP,IQREPQ,IQPSYMQ,NSFTSYMQ,
     &    IQREPMSYM,NO_SYMMETRY_CLU,NSYMACCEPTED_CLU,SYMACCEPTED_CLU
      USE MOD_SCF,ONLY:SCFSTATUS_CLU
      USE MOD_TYPES,ONLY:NTCLUMAX,X_TCLU,Z_TCLU,NA_TCLU
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CALCMODE,ONLY:IREL,KMROT,NONMAG,RELAX_CLU
      USE MOD_SITES,ONLY:QMPHI,QMTET,NOQ,NQCLU,NQHOST,NO_QCLU,
     &    QBAS0_QCLU,NQCLUMAX,ITCLU_OQCLU,IQCLU_ATCLU
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CHECK
      PARAMETER (CHECK=.TRUE.)
C
C Dummy arguments
C
      INTEGER NMCLU,NTCLU
C
C Local variables
C
      INTEGER IA,IFLAG,IO,IQ,IQ1,IQCLU,IQCLU_ATCLU0(:,:),ISYM,ITCLU,
     &        ITCLU0,ITCLU0_TCLU(:),ITCLU_OQCLU0(:,:),
     &        ITCLU_OQCLU_POTFIL(:,:),IWEDGEROT_TMP(NSYMMAX),NA_TCLU0(:)
     &        ,NMSYM_TMP,NO_QCLU_POTFIL(:),NSYMCRYSYS_TMP,NTCLU0_TMP,
     &        NTCLU_MIN,NTCLU_POTFIL,NTCLU_TMP,NWEDGE_TMP,
     &        SYMDET_TMP(NSYMMAX),Z_OQCLU(:,:)
      REAL*8 MROTK_TMP(3,3,NSYMMAX),MROTR_TMP(3,3,NSYMMAX),QMVEC_TMP(3),
     &       SYMTVEC_TMP(3,NSYMMAX),X_OQCLU(:,:)
      LOGICAL SYMCRYSYS_TMP(NSYMMAX),SYMUNITARY_TMP(NSYMMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE Z_OQCLU,X_OQCLU,NO_QCLU_POTFIL,ITCLU_OQCLU_POTFIL
      ALLOCATABLE IQCLU_ATCLU0,NA_TCLU0,ITCLU0_TCLU,ITCLU_OQCLU0
C
      WRITE (6,99003)
C
C ======================================================================
C             read  information specifying cluster occupation
C ======================================================================
C
      CALL CLUPOTRD(1)
C
      ALLOCATE (NO_QCLU_POTFIL(NQCLUMAX))
      ALLOCATE (ITCLU_OQCLU_POTFIL(NTCLUMAX,NQCLUMAX))
C
      NTCLU_POTFIL = NTCLU
      NO_QCLU_POTFIL(1:NQCLU) = NO_QCLU(1:NQCLU)
      ITCLU_OQCLU_POTFIL(1:NTCLU,1:NQCLU) = ITCLU_OQCLU(1:NTCLU,1:NQCLU)
C
C ======================================================================
C
      ALLOCATE (NQCLU_EQCLU(NQCLUMAX),IQCLU_EQCLU(NQCLUMAX,NQCLUMAX))
      ALLOCATE (IQCLU_ATCLU0(NQCLUMAX,NTCLUMAX))
      ALLOCATE (NA_TCLU0(NTCLUMAX),ITCLU0_TCLU(NTCLUMAX))
      ALLOCATE (ITCLU_OQCLU0(NTCLUMAX,NQCLUMAX))
      ALLOCATE (Z_OQCLU(NTCLUMAX,NQCLUMAX))
      ALLOCATE (X_OQCLU(NTCLUMAX,NQCLUMAX))
      ALLOCATE (IQCLU_ATCLU(NQCLUMAX,NTCLUMAX))
      ALLOCATE (NA_TCLU(NTCLUMAX))
C
C-----------------------------------------------------------------------
C              store MINIMUM original occupation information
C-----------------------------------------------------------------------
C
      DO IQCLU = 1,NQCLU
         IQ = NQHOST + IQCLU
C
         NOQ(IQ) = NO_QCLU(IQCLU)
C
         DO IO = 1,NO_QCLU(IQCLU)
C
            ITCLU = ITCLU_OQCLU(IO,IQCLU)
C
            Z_OQCLU(IO,IQCLU) = Z_TCLU(ITCLU)
            X_OQCLU(IO,IQCLU) = X_TCLU(ITCLU)
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C          transfer the occupation info to the standard variables
C-----------------------------------------------------------------------
C
      NTCLU_TMP = 0
      NA_TCLU0(1) = 0
C
      DO IQCLU = 1,NQCLU
         DO IO = 1,NO_QCLU(IQCLU)
C
            DO ITCLU = 1,NTCLU_TMP
               IF ( Z_OQCLU(IO,IQCLU).EQ.Z_TCLU(ITCLU) .AND. 
     &              ABS(X_OQCLU(IO,IQCLU)-X_TCLU(ITCLU)).LT.1D-8 ) THEN
C
                  NA_TCLU0(ITCLU) = NA_TCLU0(ITCLU) + 1
                  IA = NA_TCLU0(ITCLU)
                  IQCLU_ATCLU0(IA,ITCLU) = IQCLU
                  ITCLU_OQCLU0(IO,IQCLU) = ITCLU
                  GOTO 50
               END IF
            END DO
C
            NTCLU_TMP = NTCLU_TMP + 1
            NA_TCLU0(NTCLU_TMP) = 1
            IQCLU_ATCLU0(1,NTCLU_TMP) = IQCLU
            ITCLU_OQCLU0(IO,IQCLU) = NTCLU_TMP
            Z_TCLU(NTCLU_TMP) = Z_OQCLU(IO,IQCLU)
            X_TCLU(NTCLU_TMP) = X_OQCLU(IO,IQCLU)
C
 50      END DO
C
      END DO
C
      NTCLU_MIN = NTCLU_TMP
C
C***********************************************************************
      IF ( CHECK ) THEN
C
         WRITE (6,*)
         DO IQCLU = 1,NQCLU
            DO IO = 1,NO_QCLU(IQCLU)
               WRITE (6,99006) IQCLU,IO,ITCLU_OQCLU0(IO,IQCLU),
     &                         QBAS0_QCLU(1:3,IQCLU)
            END DO
         END DO
C
         WRITE (6,*)
C
         DO ITCLU = 1,NTCLU_TMP
            DO IA = 1,NA_TCLU0(ITCLU)
               WRITE (6,99005) ITCLU,IA,IQCLU_ATCLU0(IA,ITCLU),
     &                         Z_TCLU(ITCLU),X_TCLU(ITCLU)
            END DO
         END DO
C
      END IF
C***********************************************************************
C
C ======================================================================
C
      IF ( NO_SYMMETRY_CLU ) THEN
         NSYMCRYSYS_TMP = 1
         SYMCRYSYS_TMP(1) = .TRUE.
         SYMCRYSYS_TMP(2:NSYMMAX) = .FALSE.
      ELSE
         NSYMCRYSYS_TMP = NSYM
         SYMCRYSYS_TMP(1:NSYMMAX) = SYMACCEPTED(1:NSYMMAX)
      END IF
C
      QMVEC_TMP(1:3) = 0D0
C
      IQ1 = NQHOST + 1
C
      CALL SYMLATTICE(0,.TRUE.,IPRINT,NWEDGE_TMP,IWEDGEROT_TMP,NQCLU,
     &                QBAS0_QCLU,ABAS,NONMAG,IREL,KMROT,QMVEC_TMP,
     &                MROTR_TMP,MROTK_TMP,NSYM,NSYMCRYSYS_TMP,
     &                SYMACCEPTED_CLU,SYMUNITARY_TMP,SYMCRYSYS_TMP,
     &                SYMTVEC_TMP,SYMSYMBL,SYMEULANG,NOQ(IQ1),
     &                ITCLU_OQCLU0,ITCLU0_TCLU,NTCLU0_TMP,NTCLU_TMP,
     &                NA_TCLU0,IQCLU_ATCLU0,QMPHI(IQ1),QMTET(IQ1),
     &                IQORGQP(1,IQ1),IQPSYMQ(1,IQ1),NSFTSYMQ(1,1,IQ1),
     &                NSYMACCEPTED_CLU,SYMDET_TMP,ISYMGENQ(IQ1),
     &                IQREPQ(IQ1),NQCLU_EQCLU,IQCLU_EQCLU,NMSYM_TMP,
     &                IQREPMSYM(IQ1),NTCLUMAX,NQCLUMAX)
C
      DO IQ = NQHOST + 1,NQHOST + NQCLU
C
         IQREPQ(IQ) = NQHOST + IQREPQ(IQ)
         IQREPMSYM(IQ) = NQHOST + IQREPMSYM(IQ)
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED_CLU(ISYM) ) THEN
               IQORGQP(ISYM,IQ) = NQHOST + IQORGQP(ISYM,IQ)
               IQPSYMQ(ISYM,IQ) = NQHOST + IQPSYMQ(ISYM,IQ)
               NSFTSYMQ(1:3,ISYM,IQ) = 0
            ELSE
               IQORGQP(ISYM,IQ) = 999999
               IQPSYMQ(ISYM,IQ) = 999999
               NSFTSYMQ(1:3,ISYM,IQ) = 999999
            END IF
         END DO
C
      END DO
C
C***********************************************************************
      IF ( CHECK ) THEN
C
         WRITE (6,*)
C
C-------------------------------------------------- switch output ON/OFF
         IF ( NTCLU.LT.0 ) THEN
            DO ITCLU = 1,NTCLU
               DO IA = 1,NA_TCLU0(ITCLU)
                  WRITE (6,99002) ITCLU,IA,IQCLU_ATCLU0(IA,ITCLU),
     &                            ITCLU0_TCLU(ITCLU),Z_TCLU(ITCLU),
     &                            X_TCLU(ITCLU)
               END DO
            END DO
         END IF
C-------------------------------------------------- switch output ON/OFF
C
         DO IQCLU = 1,NQCLU
            DO IO = 1,NO_QCLU(IQCLU)
               WRITE (6,99001) IQCLU,IO,ITCLU_OQCLU0(IO,IQCLU)
            END DO
         END DO
C
      END IF
C***********************************************************************
C
      NTCLU = NTCLU_TMP
C
C-----------------------------------------------------------------------
C    CHECK consistency of read (potfil) and derived symmetry info
C-----------------------------------------------------------------------
      IF ( SCFSTATUS_CLU(1:5).NE.'START' ) THEN
C
         IFLAG = 0
         IF ( NTCLU_POTFIL.NE.NTCLU ) IFLAG = 1
C
         DO IQCLU = 1,NQCLU
            IF ( NO_QCLU_POTFIL(IQCLU).NE.NO_QCLU(IQCLU) ) IFLAG = 1
C
            DO IO = 1,NO_QCLU(IQCLU)
               IF ( ITCLU_OQCLU_POTFIL(IO,IQCLU)
     &              .NE.ITCLU_OQCLU(IO,IQCLU) ) IFLAG = 1
            END DO
C
         END DO
C
         IF ( IFLAG.NE.0 ) THEN
            WRITE (6,99007) NTCLU,NTCLU_POTFIL,NTCLU_MIN
            WRITE (6,99008) 'NO_QCLU   new versus potfil'
            WRITE (6,99009) (NO_QCLU(IQCLU),IQCLU=1,NQCLU)
            WRITE (6,99009) (NO_QCLU_POTFIL(IQCLU),IQCLU=1,NQCLU)
            WRITE (6,99008) 'ITCLU_OQCLU   new versus potfil'
            WRITE (6,99009) ((ITCLU_OQCLU(IO,IQCLU),IO=1,NO_QCLU(IQCLU))
     &                      ,IQ=1,NQCLU)
            WRITE (6,99009) ((ITCLU_OQCLU_POTFIL(IO,IQCLU),IO=1,
     &                      NO_QCLU_POTFIL(IQCLU)),IQ=1,NQCLU)
            WRITE (6,99010)
            STOP
         END IF
C
      END IF
C
C-----------------------------------------------------------------------
C            transfer occupation information to standard tables
C-----------------------------------------------------------------------
C
      DO ITCLU = 1,NTCLU_TMP
         ITCLU0 = ITCLU0_TCLU(ITCLU)
         Z_TCLU(ITCLU) = Z_TCLU(ITCLU0)
         X_TCLU(ITCLU) = X_TCLU(ITCLU0)
      END DO
C
      NA_TCLU(1:NTCLU) = NA_TCLU0(1:NTCLU)
      ITCLU_OQCLU(1:NTCLU,1:NQCLU) = ITCLU_OQCLU0(1:NTCLU,1:NQCLU)
      IQCLU_ATCLU(1:NQCLU,1:NTCLU) = IQCLU_ATCLU0(1:NQCLU,1:NTCLU)
C
C-----------------------------------------------------------------------
      IF ( .NOT.RELAX_CLU ) THEN
         NMCLU = 0
      ELSE
         NMCLU = 0
         DO IQ = NQHOST + 1,NQHOST + NQCLU
            IF ( IQREPQ(IQ).EQ.IQ ) NMCLU = NMCLU + 1
         END DO
      END IF
C-----------------------------------------------------------------------
      WRITE (6,99004) NTCLU,NMCLU
C
99001 FORMAT (10X,'IQCLU=',I4,2X,'IO=',I3,2X,'ITCLU=',I4)
99002 FORMAT (10X,'ITCLU=',I4,2X,'IA=',I3,2X,'IQCLU=',I4,2X,'IT0CLU=',
     &        I3,2X,'Z=',I3,2X,'X=',F6.3)
99003 FORMAT (2(/,1X,79('*')),/,33X,'<CLUSYMMETRY>',2(/,1X,79('*')),//,
     &        10X,'fixing the cluster symmetry and array sizes',/)
99004 FORMAT (/,10X,'NTCLU =',I6,/,10X,'NMCLU =',I6,/,2(/,1X,79('*')),/,
     &        17X,'<CLUSYMMETRY> done  --  switching back to host',
     &        2(/,1X,79('*')))
99005 FORMAT (10X,'ITCLU=',I4,2X,'IA=',I3,2X,'IQCLU=',I4,2X,'Z=',I3,2X,
     &        'X=',F6.3)
99006 FORMAT (10X,'IQCLU=',I4,2X,'IO=',I3,2X,'ITCLU=',I4,3X,'QBAS=',
     &        3F8.3)
99007 FORMAT (/,' ##### TROUBLE in <CLUSYMMETRY> ',46('#'),:,/,10X,
     &        'NTCLU        =',I3,/,10X,'NTCLU_POTFIL =',I3,/,10X,
     &        'NTCLU_MIN    =',I3)
99008 FORMAT (10X,A)
99009 FORMAT (10X,15I4)
99010 FORMAT (1X,71('#'),/)
      END
