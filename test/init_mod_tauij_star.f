C*==init_mod_tauij_star.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE INIT_MOD_TAUIJ_STAR
C   ********************************************************************
C   *                                                                  *
C   *   determine the number of TAU_ij to be calculated for a          *
C   *   given radius  CLURAD   around a lattice site IQ                *
C   *                                                                  *
C   *                                                                  *
C   *            TO BE REVISED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         *
C   *                                                                  *
C   *                                                                  *
C   *   I    numbers the unit cell (or bravais lattice vector) n       *
C   *        AND the lattice site IQ    I  ==  (n,IQ)                  *
C   *                                                                  *
C   *   if   CLURAD = 0  only the lattice sites within the             *
C   *   unit cell are considered                                       *
C   *                                                                  *
C   *   NQTAB1_TAUIJ       number of sites IQ in the central cell (n=0)*
C   *                      to be treated == inequivalent sites         *
C   *   IQ_QTAB1_TAUIJ     site index IQ (n=0) for i=1,NQTAB1_TAUIJ    *
C   *   NTAUIJ_QTAB2_TAUIJ(IQ,JQ) specifies the number of pairs        *
C   *                      (IQ,JQ) in list with JQ in any other cell n *
C   *                                                                  *
C   *   NTAUIJ             total number of TAU_ij                      *
C   *                      index IQ for site in the central cell (n=0) *
C   *                      index JQ for site in the cell at R_n(i)     *
C   *   N5VEC_TAUIJ(i)     R_n(i) = SUM(p) N5VEC_TAUIJ(i,p) * ABAS(p)  *
C   *   N123TAUIJMAX       largest  |N5VEC_TAUIJ|  that occurs         *
C   *   DQCLU_TAUIJ_CL(i)  distance between (IQ,0) and (JQ,n)          *
C   *   ITAUIJ(i)          i: TAU(nq;n'q') <-> ITAUJI(i): TAU(n'q';nq) *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,ADAINV_R,
     &    ADAINV_I,SYSTEM_DIMENSION
      USE MOD_CALCMODE,ONLY:MOL
      USE MOD_SYMMETRY,ONLY:NCL,IQ_MBCL
      USE MOD_TYPES,ONLY:Z,CONC,NT
      USE MOD_SITES,ONLY:NQMAX,QBAS,NQ_L,NQ_R,NQ,ITOQ,NOQ,IQBOT,IQTOP
      USE MOD_TAUIJ,ONLY:NSHLCLU_TAUIJ_CL,CLURAD_TAUIJ_CL,
     &    NQCLU_TAUIJ_CL,RQCLU_TAUIJ_CL,DQCLU_TAUIJ_CL,
     &    NSHLCLUMAX_TAUIJ_CL,IQ_TAUIJ_CL,N5VEC_TAUIJ_CL,
     &    NQCLUMAX_TAUIJ_CL,NQSHLCLU_TAUIJ_CL,ITAUIJ_LOOP,ITAUJI_LOOP,
     &    NLOOP_TAUIJ,N5VEC_TAUIJ,SELECTED_TAUIJ_CL
      USE MOD_FILES,ONLY:IPRINT,IOTMP,FOUND_SECTION,FOUND_REAL,
     &    FOUND_INTEGER
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NSELEXCMAX
      PARAMETER (NSELEXCMAX=100)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_TAUIJ_STAR')
C
C Local variables
C
      INTEGER ATOMNUMB
      REAL*8 CLURAD,CLURAD0,DVECIJ(3),DVECJI(3),D_TAUIJ
      REAL*8 DNRM2
      LOGICAL FOUND,IQ_SELECTED,TROUBLE
      INTEGER I,ICL,ILOOP,IPRINTLOC,IQ,IQCNTR,IQ_TAUIJ_TMP(:),ITAUIJ,
     &        ITAUIJ_ORIG,ITAUIJ_SORT,ITAUIJ_SORT_ORIG(:),ITAUJI,
     &        ITAUJI_ORIG,ITAUJI_SORT,J,JQ,JQCLU,JQ_TAUIJ_TMP(:),
     &        N5VEC_IQCLU(5),N5VEC_JQCLU(5),N5VEC_TAUIJ_TMP(:,:),
     &        NEXCLUDE,NQCLU,NQCLU_I,NQCLU_L,NQCLU_R,NSELECT,NSHLCLU,
     &        NSHLCLU0,NTAUIJMAX_TMP,NTAUIJ_TMP,
     &        Z_EXCLUDE_LIST(NSELEXCMAX),Z_SELECT_LIST(NSELEXCMAX)
      CHARACTER*2 STR2_LIST(NSELEXCMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE N5VEC_TAUIJ_TMP,IQ_TAUIJ_TMP,JQ_TAUIJ_TMP
      ALLOCATABLE ITAUIJ_SORT_ORIG
C
C-----------------------------------------------------------------------
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (NSHLCLU_TAUIJ_CL(NCL),CLURAD_TAUIJ_CL(NCL))
      ALLOCATE (NQCLU_TAUIJ_CL(NCL))
C
      IPRINTLOC = -2
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
C=======================================================================
C   select or exclude lattice sites by selecting or exluding elements
C=======================================================================
C
      NSELECT = 0
      NEXCLUDE = 0
C
      IF ( FOUND_SECTION ) THEN
C
         CALL SECTION_SET_STRING_ARRAY('SELECT',STR2_LIST,NSELECT,
     &                                 NSELEXCMAX,0,'9999',0)
C
         CALL SECTION_SET_STRING_ARRAY('EXCLUDE',STR2_LIST,NEXCLUDE,
     &                                 NSELEXCMAX,0,'9999',0)
C
         IF ( NSELECT.NE.0 .AND. NEXCLUDE.NE.0 )
     &        CALL STOP_MESSAGE(ROUTINE,
     &        'SELECT  AND EXCLUDE not allowed')
C
         DO I = 1,NSELECT
            Z_SELECT_LIST(I) = ATOMNUMB(STR2_LIST(I))
         END DO
         DO I = 1,NEXCLUDE
            Z_EXCLUDE_LIST(I) = ATOMNUMB(STR2_LIST(I))
         END DO
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'section TASK not found')
      END IF
C
C=======================================================================
C                   fix the cluster radius
C=======================================================================
C
C------------------------------------------------- set NSHLCLU or CLURAD
C
      NSHLCLU = 0
      CALL SECTION_SET_INTEGER('NSHLCLU',NSHLCLU,9999,0)
C
      IF ( .NOT.FOUND_INTEGER ) THEN
         CALL SECTION_SET_REAL('CLURAD',CLURAD,9999D0,0)
         IF ( .NOT.FOUND_REAL ) CALL STOP_MESSAGE(ROUTINE,
     &        'NSHLCLU nor CLURAD set in section TASK')
      END IF
C
C-----------------------------------------------------------------------
C     find the site-specific cluster radii  CLURAD and array sizes
C-----------------------------------------------------------------------
C
      CLURAD0 = CLURAD
      NSHLCLU0 = NSHLCLU
      NQCLUMAX_TAUIJ_CL = 0
      NSHLCLUMAX_TAUIJ_CL = 0
C
      DO ICL = 1,NCL
C
         IQCNTR = IQ_MBCL(1,ICL)
         CLURAD = CLURAD0
         NSHLCLU = NSHLCLU0
C
         CALL CLUSSITES(IOTMP,IPRINTLOC,MOL,SYSTEM_DIMENSION,ABAS,
     &                  ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,ADAINV_R,
     &                  QBAS,CLURAD,IQCNTR,NQCLU,NQCLU_L,NQCLU_I,
     &                  NQCLU_R,NSHLCLU,NQ,NQ_L,NQ_R,NQMAX)
         CLOSE (IOTMP)
C
         CLURAD_TAUIJ_CL(ICL) = CLURAD
         NSHLCLU_TAUIJ_CL(ICL) = NSHLCLU
         NQCLU_TAUIJ_CL(ICL) = NQCLU
         NQCLUMAX_TAUIJ_CL = MAX(NQCLUMAX_TAUIJ_CL,NQCLU)
         NSHLCLUMAX_TAUIJ_CL = MAX(NSHLCLUMAX_TAUIJ_CL,NSHLCLU)
C
         IF ( IPRINT.GT.0 ) WRITE (6,99002) ICL,CLURAD_TAUIJ_CL(ICL),
     &                             NSHLCLU_TAUIJ_CL(ICL),
     &                             NQCLU_TAUIJ_CL(ICL)
C
      END DO
C
C-----------------------------------------------------------------------
C                  allocate arrays in module
C-----------------------------------------------------------------------
C
      ALLOCATE (SELECTED_TAUIJ_CL(NQCLUMAX_TAUIJ_CL,NCL))
      ALLOCATE (RQCLU_TAUIJ_CL(3,NQCLUMAX_TAUIJ_CL,NCL))
      ALLOCATE (DQCLU_TAUIJ_CL(NQCLUMAX_TAUIJ_CL,NCL))
      ALLOCATE (IQ_TAUIJ_CL(NQCLUMAX_TAUIJ_CL,NCL))
      ALLOCATE (N5VEC_TAUIJ_CL(5,NQCLUMAX_TAUIJ_CL,NCL))
      ALLOCATE (NQSHLCLU_TAUIJ_CL(NSHLCLUMAX_TAUIJ_CL,NCL))
C
      SELECTED_TAUIJ_CL(1:NQCLUMAX_TAUIJ_CL,1:NCL) = .FALSE.
C
C-----------------------------------------------------------------------
C                 allocate temporary work storage
C-----------------------------------------------------------------------
      NTAUIJMAX_TMP = 2*NQCLUMAX_TAUIJ_CL*NCL
C
      ALLOCATE (N5VEC_TAUIJ_TMP(5,NTAUIJMAX_TMP))
      ALLOCATE (IQ_TAUIJ_TMP(NTAUIJMAX_TMP),JQ_TAUIJ_TMP(NTAUIJMAX_TMP))
C
      NLOOP_TAUIJ = 2*SUM(NQCLU_TAUIJ_CL(1:NCL))
      ALLOCATE (ITAUIJ_LOOP(NLOOP_TAUIJ),ITAUJI_LOOP(NLOOP_TAUIJ))
      ITAUIJ_LOOP(1:NLOOP_TAUIJ) = 0
      ITAUJI_LOOP(1:NLOOP_TAUIJ) = 0
C
      N5VEC_TAUIJ_TMP(1:5,1:NTAUIJMAX_TMP) = 999999
      IQ_TAUIJ_TMP(1:NTAUIJMAX_TMP) = 999999
      JQ_TAUIJ_TMP(1:NTAUIJMAX_TMP) = 999999
C
C-----------------------------------------------------------------------
C       run over the loops for which the TAUIJ's are needed lateron
C               and set up the list of non-equivalent TAUIJ's
C       NOTE: use the same loop structure everywhere because of the
C             linear pointer  ITAUIJ_LOOP -----> TAUIJ(ITAUIJ)
C-----------------------------------------------------------------------
C
C------------------------------------------------------- deal with TAUIJ
      NTAUIJ_TMP = 0
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQCNTR = IQ_MBCL(1,ICL)
C
         IQ_SELECTED = .FALSE.
C
C
         IF ( IQCNTR.GE.IQBOT .AND. IQCNTR.LE.IQTOP )
     &        CALL SITES_SELECT(NSELECT,Z_SELECT_LIST,NEXCLUDE,
     &        Z_EXCLUDE_LIST,NSELEXCMAX,Z,CONC,NT,NOQ(IQCNTR),
     &        ITOQ(1,IQCNTR),IQ_SELECTED)
C
         IF ( IQ_SELECTED ) THEN
C
            CLURAD = CLURAD_TAUIJ_CL(ICL)
            NSHLCLU = 0
C
            CALL CLUSSITES(IOTMP,IPRINTLOC,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,CLURAD,IQCNTR,NQCLU,NQCLU_L,
     &                     NQCLU_I,NQCLU_R,NSHLCLU,NQ,NQ_L,NQ_R,NQMAX)
C
            READ (IOTMP) ((RQCLU_TAUIJ_CL(J,I,ICL),J=1,3),DQCLU_TAUIJ_CL
     &                   (I,ICL),IQ_TAUIJ_CL(I,ICL),I=1,NQCLU),
     &                   (NQSHLCLU_TAUIJ_CL(I,ICL),I=1,NSHLCLU)
            READ (IOTMP) ((N5VEC_TAUIJ_CL(J,I,ICL),J=1,5),I=1,NQCLU)
C
            IF ( NQCLU.NE.NQCLU_TAUIJ_CL(ICL) )
     &           CALL STOP_MESSAGE(ROUTINE,' --> NQCLU')
C
            CLOSE (IOTMP)
C
            IQ = IQ_MBCL(1,ICL)
            N5VEC_IQCLU(1:5) = 0
C
            IF ( IPRINT.GT.0 ) WRITE (6,99003)
C
         END IF
C
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
            IF ( IQ_SELECTED ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               IF ( JQ.GE.IQBOT .AND. JQ.LE.IQTOP ) THEN
C
                  CALL SITES_SELECT(NSELECT,Z_SELECT_LIST,NEXCLUDE,
     &                              Z_EXCLUDE_LIST,NSELEXCMAX,Z,CONC,NT,
     &                              NOQ(JQ),ITOQ(1,JQ),
     &                              SELECTED_TAUIJ_CL(JQCLU,ICL))
               ELSE
C
                  SELECTED_TAUIJ_CL(JQCLU,ICL) = .FALSE.
C
               END IF
C
            ELSE
C
               SELECTED_TAUIJ_CL(JQCLU,ICL) = .FALSE.
C
            END IF
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               N5VEC_JQCLU(1:5) = N5VEC_TAUIJ_CL(1:5,JQCLU,ICL)
C
               ILOOP = ILOOP + 1
C
               CALL INIT_MOD_TAUIJ_LOOP(ITAUIJ_LOOP,ILOOP,NLOOP_TAUIJ,
     &                                  IQ,N5VEC_IQCLU,JQ,N5VEC_JQCLU,
     &                                  NTAUIJ_TMP,IQ_TAUIJ_TMP,
     &                                  JQ_TAUIJ_TMP,N5VEC_TAUIJ_TMP,
     &                                  NTAUIJMAX_TMP)
C
               IF ( IPRINT.GT.0 ) WRITE (6,99004) ILOOP,ICL,IQ,JQCLU,JQ,
     &              N5VEC_JQCLU(1:5),RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &              ITAUIJ_LOOP(ILOOP),NTAUIJ_TMP
C
            END IF
C-----------------------------------------------------------------------
         END DO
C
      END DO
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,*) ' '
         DO ITAUIJ = 1,NTAUIJ_TMP
C
            WRITE (6,99001) ITAUIJ,IQ_TAUIJ_TMP(ITAUIJ),
     &                      JQ_TAUIJ_TMP(ITAUIJ),
     &                      N5VEC_TAUIJ_TMP(1:5,ITAUIJ)
         END DO
         WRITE (6,*) ' '
      END IF
C
C------------------------------------------------------- deal with TAUJI
C   reverse order of I and J related arguments in INIT_MOD_TAUIJ_LOOP
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
         N5VEC_IQCLU(1:5) = 0
C
         IF ( IPRINT.GT.0 ) WRITE (6,99003)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
               N5VEC_JQCLU(1:5) = N5VEC_TAUIJ_CL(1:5,JQCLU,ICL)
C
               ILOOP = ILOOP + 1
C
               CALL INIT_MOD_TAUIJ_LOOP(ITAUJI_LOOP,ILOOP,NLOOP_TAUIJ,
     &                                  JQ,N5VEC_JQCLU,IQ,N5VEC_IQCLU,
     &                                  NTAUIJ_TMP,IQ_TAUIJ_TMP,
     &                                  JQ_TAUIJ_TMP,N5VEC_TAUIJ_TMP,
     &                                  NTAUIJMAX_TMP)
C
               IF ( IPRINT.GT.0 ) WRITE (6,99004) ILOOP,ICL,IQ,JQCLU,JQ,
     &              N5VEC_JQCLU(1:5),RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &              ITAUJI_LOOP(ILOOP),NTAUIJ_TMP
C
            END IF
C-----------------------------------------------------------------------
C
         END DO
C
      END DO
C
      WRITE (6,99009) NTAUIJ_TMP
C
C-----------------------------------------------------------------------
C               sort the list of non-equivalent TAUIJ's
C              to be optimized w.r.t. the BZ integration
C-----------------------------------------------------------------------
C
      ALLOCATE (ITAUIJ_SORT_ORIG(NTAUIJ_TMP))
      ITAUIJ_SORT_ORIG(1:NTAUIJ_TMP) = 999999
C
      CALL INIT_MOD_TAUIJ_TABLE(NTAUIJ_TMP,ITAUIJ_SORT_ORIG,
     &                          IQ_TAUIJ_TMP,JQ_TAUIJ_TMP,
     &                          N5VEC_TAUIJ_TMP,NTAUIJMAX_TMP)
C
C-----------------------------------------------------------------------
C         account for sorting the list of non-equivalent TAUIJ's
C-----------------------------------------------------------------------
C
      ILOOP = 0
C
      IF ( IPRINT.GT.0 ) WRITE (6,99005)
      DO ICL = 1,NCL
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               ILOOP = ILOOP + 1
C
               ITAUIJ_ORIG = ITAUIJ_LOOP(ILOOP)
               ITAUIJ_SORT = ITAUIJ_SORT_ORIG(ITAUIJ_ORIG)
               ITAUIJ_LOOP(ILOOP) = ITAUIJ_SORT
C
               ITAUJI_ORIG = ITAUJI_LOOP(ILOOP)
               ITAUJI_SORT = ITAUIJ_SORT_ORIG(ITAUJI_ORIG)
               ITAUJI_LOOP(ILOOP) = ITAUJI_SORT
C
               IF ( IPRINT.GT.0 ) THEN
                  WRITE (6,99006) 'IJ',ILOOP,ICL,JQCLU,ITAUIJ_ORIG,
     &                            ITAUIJ_SORT
                  WRITE (6,99006) 'JI',ILOOP,ICL,JQCLU,ITAUJI_ORIG,
     &                            ITAUJI_SORT
               END IF
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
      DEALLOCATE (ITAUIJ_SORT_ORIG)
      DEALLOCATE (N5VEC_TAUIJ_TMP,IQ_TAUIJ_TMP,JQ_TAUIJ_TMP)
C
C=======================================================================
C                check consistency of final tables
C=======================================================================
C
      TROUBLE = .FALSE.
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               ITAUJI = ITAUJI_LOOP(ILOOP)
C
C-----------------------------------------------------------------------
               FOUND = .FALSE.
               DO I = 1,5
                  FOUND = FOUND .OR. N5VEC_TAUIJ_CL(I,JQCLU,ICL)
     &                    .NE.N5VEC_TAUIJ(I,ITAUIJ)
               END DO
               IF ( FOUND ) THEN
                  WRITE (6,99007) 'TROUBLE: CHECK A: ICL IQ JQCLU ',ICL,
     &                            IQ,JQCLU,N5VEC_TAUIJ_CL(I,JQCLU,ICL),
     &                            N5VEC_TAUIJ(I,ITAUIJ)
                  TROUBLE = .TRUE.
               END IF
C
C-----------------------------------------------------------------------
C
               CALL RVECLCIVB(3,N5VEC_TAUIJ_CL(1,JQCLU,ICL),ABAS,DVECIJ)
C
               DVECIJ(1:3) = DVECIJ(1:3) + QBAS(1:3,JQ) - QBAS(1:3,IQ)
C
               FOUND = .FALSE.
               DO I = 1,3
                  FOUND = FOUND .OR. 
     &                    ABS(DVECIJ(I)-RQCLU_TAUIJ_CL(I,JQCLU,ICL))
     &                    .GT.1D-6
               END DO
               IF ( FOUND ) THEN
                  WRITE (6,99008) 'TROUBLE: CHECK B: ICL IQ JQCLU ',ICL,
     &                            IQ,JQCLU,DVECIJ(1:3),
     &                            RQCLU_TAUIJ_CL(1:3,JQCLU,ICL)
                  TROUBLE = .TRUE.
               END IF
C
C-----------------------------------------------------------------------
C
               CALL RVECLCIVB(3,N5VEC_TAUIJ(1,ITAUIJ),ABAS,DVECIJ)
C
               CALL RVECLCIVB(3,N5VEC_TAUIJ(1,ITAUJI),ABAS,DVECJI)
C
               DO I = 1,3
                  IF ( ABS(DVECIJ(I)+DVECJI(I)).GT.1D-6 ) THEN
                     WRITE (6,*) 'TROUBLE: RIJ-RJI: ITAUIJ ---> inverse'
     &                           ,ITAUIJ,ITAUJI
                     TROUBLE = .TRUE.
                  END IF
               END DO
C
C-----------------------------------------------------------------------
C
               DVECIJ(1:3) = DVECIJ(1:3) + QBAS(1:3,JQ) - QBAS(1:3,IQ)
C
               FOUND = .FALSE.
               DO I = 1,3
                  FOUND = FOUND .OR. 
     &                    ABS(DVECIJ(I)-RQCLU_TAUIJ_CL(I,JQCLU,ICL))
     &                    .GT.1D-6
               END DO
               IF ( FOUND ) THEN
                  WRITE (6,99008) 'TROUBLE: RIJ-RJI: ICL IQ JQCLU ',ICL,
     &                            IQ,JQCLU,DVECIJ(1:3),
     &                            RQCLU_TAUIJ_CL(1:3,JQCLU,ICL)
                  TROUBLE = .TRUE.
               END IF
C
               D_TAUIJ = DNRM2(3,DVECIJ,1)
C
               IF ( ABS(D_TAUIJ-DQCLU_TAUIJ_CL(JQCLU,ICL)).GT.1D-6 )
     &              THEN
                  TROUBLE = .TRUE.
                  WRITE (6,99008) 'TROUBLE: RIJ-RJI: ICL IQ JQCLU ',ICL,
     &                            IQ,JQCLU,D_TAUIJ,
     &                            DQCLU_TAUIJ_CL(JQCLU,ICL)
               END IF
C
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
      IF ( TROUBLE ) CALL STOP_MESSAGE(ROUTINE,'TROUBLE')
C
99001 FORMAT (3I5,5I3)
99002 FORMAT (/,5X,'ICL =',I3,'  CLURAD =',F10.4,'  NSHLCLU =',I6,
     &        '  NQCLU =',I6,/)
99003 FORMAT (/,' ILOOP  ICL   IQ   JQCLU   JQ        N5VEC',
     &        '             RQCLU        ITAUIJ')
99004 FORMAT (1X,3I5,3X,2I5,2X,5I3,2X,3F6.2,I8,I4)
99005 FORMAT (/,5X,'ILOOP  ICL   JQCLU   ITAUIJ (orig)   ITAUIJ (sort)')
99006 FORMAT (2X,A,1X,2I5,i8,I9,I16)
C
99007 FORMAT (5X,A,3I4,2X,5I3,/,50X,5I3,/)
99008 FORMAT (5X,A,3I4,2X,3F8.4,/,50X,3F8.4,/)
99009 FORMAT (/,1X,79('*'),/,29X,'<INIT_MOD_TAUIJ_STAR>',/,1X,79('*'),
     &        //,5X,'number of TAU(I,J)''s to be calculated:',I5,/)
C
      END
C*==sites_select.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SITES_SELECT(NSELECT,Z_SELECT_LIST,NEXCLUDE,
     &                        Z_EXCLUDE_LIST,NSELEXCMAX,Z,CONC,NT,NOQ,
     &                        ITOQ,IQ_SELECTED)
C   ********************************************************************
C   *                                                                  *
C   *   select or exclude a site according to occupation and           *
C   *   select or exclude lists                                        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-6)
C
C Dummy arguments
C
      LOGICAL IQ_SELECTED
      INTEGER NEXCLUDE,NOQ,NSELECT,NSELEXCMAX,NT
      REAL*8 CONC(NT)
      INTEGER ITOQ(NOQ),Z(NT),Z_EXCLUDE_LIST(NSELEXCMAX),
     &        Z_SELECT_LIST(NSELEXCMAX)
C
C Local variables
C
      INTEGER I,IO,IT
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C SELECT: select a site if it is occupied to an appreciable extent with
C         one of the selected elements; i.e.  with concentration > TOL
C
      IF ( NSELECT.GT.0 ) THEN
C
         IQ_SELECTED = .FALSE.
C
         DO IO = 1,NOQ
            IT = ITOQ(IO)
            DO I = 1,NSELECT
               IF ( Z(IT).EQ.Z_SELECT_LIST(I) ) THEN
                  IF ( CONC(IT).GT.TOL ) THEN
                     IQ_SELECTED = .TRUE.
                     RETURN
                  END IF
               END IF
            END DO
         END DO
C
C-----------------------------------------------------------------------
C EXCLUDE: exclude a site if it is occupied exclusivly with one of the
C          elements to be excluded; i.e with concentration > 1-TOL
C
      ELSE IF ( NEXCLUDE.GT.0 ) THEN
C
         IQ_SELECTED = .TRUE.
C
         DO IO = 1,NOQ
            IT = ITOQ(IO)
            DO I = 1,NEXCLUDE
               IF ( Z(IT).EQ.Z_EXCLUDE_LIST(I) ) THEN
                  IF ( CONC(IT).GT.(1D0-TOL) ) THEN
                     IQ_SELECTED = .FALSE.
                     RETURN
                  END IF
               END IF
            END DO
         END DO
C
C-----------------------------------------------------------------------
C DEFAULT: select site
C
      ELSE
C
         IQ_SELECTED = .TRUE.
C
      END IF
C-----------------------------------------------------------------------
      END
