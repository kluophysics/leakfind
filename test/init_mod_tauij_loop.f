C*==init_mod_tauij_loop.f    processed by SPAG 6.55Rc at 18:26 on 18 Mar 2008
      SUBROUTINE INIT_MOD_TAUIJ_LOOP(ITAUIJ_LOOP,ILOOP,NLOOP_TAUIJ,IQ,
     &                               N5VEC_IQ,JQ,N5VEC_JQ,NTAUIJ,
     &                               IQ_TAUIJ,JQ_TAUIJ,N5VEC_TAUIJ,
     &                               NTAUIJMAX)
C   ********************************************************************
C   *                                                                  *
C   *   look up whether present site pair (IQ,N5VEC_IQ)-(JQ,N5VEC_JQ)  *
C   *   is in list of non-equivalent  TAUIJ's                          *
C   *   if not: extent list                                            *
C   *   return index  ITAUIJ  for present loop position  ILOOP         *
C   *                                                                  *
C   *   NOTE: the list of TAUIJ depends only on  N5VEC_JQ - N5VEC_IQ   *
C   *         according to the definition of TAUIJ                     *
C   *         TAUIJ = INT d3k  TAU[IQ,JQ](k) * exp(ik*(R_JQ - R_IQ))   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ILOOP,IQ,JQ,NLOOP_TAUIJ,NTAUIJ,NTAUIJMAX
      INTEGER IQ_TAUIJ(NTAUIJMAX),ITAUIJ_LOOP(NLOOP_TAUIJ),
     &        JQ_TAUIJ(NTAUIJMAX),N5VEC_IQ(5),N5VEC_JQ(5),
     &        N5VEC_TAUIJ(5,NTAUIJMAX)
C
C Local variables
C
      INTEGER DEL_IJ,IBV,ITAUIJ
C
C*** End of declarations rewritten by SPAG
C
      DO ITAUIJ = 1,NTAUIJ
         IF ( IQ.EQ.IQ_TAUIJ(ITAUIJ) ) THEN
            IF ( JQ.EQ.JQ_TAUIJ(ITAUIJ) ) THEN
               DO IBV = 1,5
                  DEL_IJ = N5VEC_JQ(IBV) - N5VEC_IQ(IBV)
                  IF ( DEL_IJ.NE.N5VEC_TAUIJ(IBV,ITAUIJ) ) GOTO 100
               END DO
C--------------------------------------------- entry found in the table
               ITAUIJ_LOOP(ILOOP) = ITAUIJ
               RETURN
C-----------------------------------------------------------------------
            END IF
         END IF
C
 100  END DO
C
C----------------------------------------------------------------------
C                     new entry to the table
C----------------------------------------------------------------------
C
      NTAUIJ = NTAUIJ + 1
C
      IF ( NTAUIJ.GT.NTAUIJMAX ) STOP 
     &                   'in <INIT_MOD_TAUIJ_LOOP>:  NTAUIJMAX exceeded'
C
      ITAUIJ = NTAUIJ
      IQ_TAUIJ(ITAUIJ) = IQ
      JQ_TAUIJ(ITAUIJ) = JQ
      DO IBV = 1,5
         DEL_IJ = N5VEC_JQ(IBV) - N5VEC_IQ(IBV)
         N5VEC_TAUIJ(IBV,ITAUIJ) = DEL_IJ
      END DO
C
      ITAUIJ_LOOP(ILOOP) = ITAUIJ
C
      END
C*==init_mod_tauij_table.f    processed by SPAG 6.55Rc at 18:26 on 18 Mar 2008
      SUBROUTINE INIT_MOD_TAUIJ_TABLE(NTAUIJ_TMP,ITAUIJ_SORT_ORIG,
     &                                IQ_TAUIJ_TMP,JQ_TAUIJ_TMP,
     &                                N5VEC_TAUIJ_TMP,NTAUIJMAX_TMP)
C   ********************************************************************
C   *                                                                  *
C   *   sort the list of non-equivalent TAUIJ's                        *
C   *   to be optimized w.r.t. the BZ integration                      *
C   *   and set up the additional tables for that                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_TAUIJ,ONLY:NTAUIJ,IQ_QTAB1_TAUIJ,JQ_QTAB2_TAUIJ,
     &    NQTAB2_TAUIJ,JQTAB_TAUIJMAX,IQ_TAUIJ,JQ_TAUIJ,N5VEC_TAUIJ,
     &    N123TAUIJMAX,NTAUIJ_QTAB2_TAUIJ,NQTAB1_TAUIJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NTAUIJMAX_TMP,NTAUIJ_TMP
      INTEGER IQ_TAUIJ_TMP(NTAUIJMAX_TMP),ITAUIJ_SORT_ORIG(NTAUIJ_TMP),
     &        JQ_TAUIJ_TMP(NTAUIJMAX_TMP),
     &        N5VEC_TAUIJ_TMP(5,NTAUIJMAX_TMP)
C
C Local variables
C
      INTEGER CHKSUMI1,CHKSUMI2,CHKSUMJ1,CHKSUMJ2,CHKSUMN1(5),
     &        CHKSUMN2(5),I,IQ,IQTAB,IQ_LAST,IQ_MIN,ITAUIJ,ITAUIJ1,
     &        ITAUIJ2,ITAUIJ_TAB,ITAUIJ_TMP,JQ,JQTAB,JQ_LAST,JQ_MIN,N,
     &        N5_MIN,NQTAB2_TAUIJ_TMP(:),NTAUIJ_DONE
      LOGICAL KDONE(:),SELECT(:),TROUBLE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KDONE,SELECT,NQTAB2_TAUIJ_TMP
C
      ALLOCATE (KDONE(NTAUIJ_TMP),SELECT(NTAUIJ_TMP))
C
      NTAUIJ = NTAUIJ_TMP
C
      TROUBLE = .FALSE.
C
C-----------------------------------------------------------------------
C                    find array sizes
C-----------------------------------------------------------------------
C
      ALLOCATE (NQTAB2_TAUIJ_TMP(NTAUIJ))
C
      KDONE(1:NTAUIJ) = .FALSE.
      NQTAB1_TAUIJ = 0
      NQTAB2_TAUIJ_TMP(1:NTAUIJ) = 0
      IQ_LAST = 0
      JQ_LAST = 0
      JQTAB_TAUIJMAX = 0
C
      DO ITAUIJ_TAB = 1,NTAUIJ
C
         SELECT(1:NTAUIJ) = .NOT.KDONE(1:NTAUIJ)
C
         IQ_MIN = 10000
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) IQ_MIN = MIN(IQ_MIN,IQ_TAUIJ_TMP(I))
         END DO
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               IF ( IQ_MIN.NE.IQ_TAUIJ_TMP(I) ) SELECT(I) = .FALSE.
            END IF
         END DO
C
         JQ_MIN = 10000
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) JQ_MIN = MIN(JQ_MIN,JQ_TAUIJ_TMP(I))
         END DO
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               IF ( JQ_MIN.NE.JQ_TAUIJ_TMP(I) ) SELECT(I) = .FALSE.
            END IF
         END DO
C
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               ITAUIJ_TMP = I
               KDONE(ITAUIJ_TMP) = .TRUE.
C
               IQ = IQ_TAUIJ_TMP(ITAUIJ_TMP)
               JQ = JQ_TAUIJ_TMP(ITAUIJ_TMP)
C
               IF ( IQ.GT.IQ_LAST ) THEN
                  NQTAB1_TAUIJ = NQTAB1_TAUIJ + 1
                  JQ_LAST = 0
               END IF
               IQTAB = NQTAB1_TAUIJ
C
               IF ( JQ.GT.JQ_LAST ) NQTAB2_TAUIJ_TMP(IQTAB)
     &              = NQTAB2_TAUIJ_TMP(IQTAB) + 1
               JQTAB = NQTAB2_TAUIJ_TMP(IQTAB)
C
               JQTAB_TAUIJMAX = MAX(JQTAB_TAUIJMAX,
     &                          NQTAB2_TAUIJ_TMP(IQTAB))
C
               IQ_LAST = IQ
               JQ_LAST = JQ
               GOTO 100
            END IF
         END DO
         STOP '<INIT_MOD_TAUIJ_TABLE> - LOOP 1: NO ITAUIJ found'
 100  END DO
C
      DEALLOCATE (NQTAB2_TAUIJ_TMP)
C
C-----------------------------------------------------------------------
C                   allocate storage and swap information
C-----------------------------------------------------------------------
C
      ALLOCATE (NQTAB2_TAUIJ(NQTAB1_TAUIJ))
      ALLOCATE (IQ_QTAB1_TAUIJ(NQTAB1_TAUIJ))
      ALLOCATE (JQ_QTAB2_TAUIJ(NQTAB1_TAUIJ,JQTAB_TAUIJMAX))
      ALLOCATE (NTAUIJ_QTAB2_TAUIJ(NQTAB1_TAUIJ,JQTAB_TAUIJMAX))
      ALLOCATE (IQ_TAUIJ(NTAUIJ),JQ_TAUIJ(NTAUIJ))
      ALLOCATE (N5VEC_TAUIJ(5,NTAUIJ))
C
      KDONE(1:NTAUIJ) = .FALSE.
      NQTAB2_TAUIJ(1:NQTAB1_TAUIJ) = 0
      IQ_QTAB1_TAUIJ(1:NQTAB1_TAUIJ) = 999999
      JQ_QTAB2_TAUIJ(1:NQTAB1_TAUIJ,1:JQTAB_TAUIJMAX) = 999999
      NTAUIJ_QTAB2_TAUIJ(1:NQTAB1_TAUIJ,1:JQTAB_TAUIJMAX) = 0
      IQ_TAUIJ(1:NTAUIJ) = 999999
      JQ_TAUIJ(1:NTAUIJ) = 999999
      N5VEC_TAUIJ(1:5,1:NTAUIJ) = 999999
C
      NQTAB1_TAUIJ = 0
      NTAUIJ_DONE = 0
      IQ_LAST = 0
      JQ_LAST = 0
C
      DO ITAUIJ_TAB = 1,NTAUIJ
C
         SELECT(1:NTAUIJ) = .NOT.KDONE(1:NTAUIJ)
C
         IQ_MIN = 10000
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) IQ_MIN = MIN(IQ_MIN,IQ_TAUIJ_TMP(I))
         END DO
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               IF ( IQ_MIN.NE.IQ_TAUIJ_TMP(I) ) SELECT(I) = .FALSE.
            END IF
         END DO
C
         JQ_MIN = 10000
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) JQ_MIN = MIN(JQ_MIN,JQ_TAUIJ_TMP(I))
         END DO
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               IF ( JQ_MIN.NE.JQ_TAUIJ_TMP(I) ) SELECT(I) = .FALSE.
            END IF
         END DO
C
         DO N = 1,5
            N5_MIN = 10000
            DO I = 1,NTAUIJ
               IF ( SELECT(I) ) N5_MIN = MIN(N5_MIN,N5VEC_TAUIJ_TMP(N,I)
     &              )
            END DO
            DO I = 1,NTAUIJ
               IF ( SELECT(I) ) THEN
                  IF ( N5_MIN.NE.N5VEC_TAUIJ_TMP(N,I) ) SELECT(I)
     &                 = .FALSE.
               END IF
            END DO
         END DO
C
         DO I = 1,NTAUIJ
            IF ( SELECT(I) ) THEN
               ITAUIJ_TMP = I
               KDONE(ITAUIJ_TMP) = .TRUE.
               NTAUIJ_DONE = NTAUIJ_DONE + 1
               ITAUIJ = NTAUIJ_DONE
C
               IQ = IQ_TAUIJ_TMP(ITAUIJ_TMP)
               JQ = JQ_TAUIJ_TMP(ITAUIJ_TMP)
C
               IF ( IQ.GT.IQ_LAST ) THEN
                  NQTAB1_TAUIJ = NQTAB1_TAUIJ + 1
                  JQ_LAST = 0
               END IF
               IQTAB = NQTAB1_TAUIJ
C
               IF ( JQ.GT.JQ_LAST ) NQTAB2_TAUIJ(IQTAB)
     &              = NQTAB2_TAUIJ(IQTAB) + 1
               JQTAB = NQTAB2_TAUIJ(IQTAB)
C
               IQ_QTAB1_TAUIJ(IQTAB) = IQ
               JQ_QTAB2_TAUIJ(IQTAB,JQTAB) = JQ
C
               NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
     &            = NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB) + 1
C
               ITAUIJ_SORT_ORIG(ITAUIJ_TMP) = ITAUIJ
C
               IQ_TAUIJ(ITAUIJ) = IQ
               JQ_TAUIJ(ITAUIJ) = JQ
               N5VEC_TAUIJ(1:5,ITAUIJ) = N5VEC_TAUIJ_TMP(1:5,ITAUIJ_TMP)
C
               IQ_LAST = IQ
               JQ_LAST = JQ
C
               IF ( IPRINT.GT.0 ) WRITE (6,99001) ITAUIJ_TAB,ITAUIJ_TMP,
     &              ITAUIJ,IQ_TAUIJ(ITAUIJ),JQ_TAUIJ(ITAUIJ),
     &              N5VEC_TAUIJ(1:5,ITAUIJ)
C
               GOTO 200
            END IF
         END DO
         STOP '<INIT_MOD_TAUIJ_TABLE> - LOOP 2: NO ITAUIJ found'
 200  END DO
C
      IF ( NTAUIJ_DONE.NE.NTAUIJ ) THEN
         WRITE (6,*) '<INIT_MOD_TAUIJ_TABLE>: NTAUIJ_DONE .NE. NTAUIJ'
         TROUBLE = .TRUE.
      END IF
C
C------------------ check with the loop structure for the BZ integration
C
      CHKSUMI1 = 0
      CHKSUMJ1 = 0
      CHKSUMN1(1:5) = 0
      ITAUIJ2 = 0
      N123TAUIJMAX = 0
C
      DO IQTAB = 1,NQTAB1_TAUIJ
         IQ = IQ_QTAB1_TAUIJ(IQTAB)
         DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
            JQ = JQ_QTAB2_TAUIJ(IQTAB,JQTAB)
            ITAUIJ1 = ITAUIJ2 + 1
            ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
            DO ITAUIJ = ITAUIJ1,ITAUIJ2
               CHKSUMI1 = CHKSUMI1 + IQ
               CHKSUMJ1 = CHKSUMJ1 + JQ
               CHKSUMN1(1:5) = CHKSUMN1(1:5) + N5VEC_TAUIJ(1:5,ITAUIJ)
               DO I = 1,5
                  N123TAUIJMAX = MAX(ABS(N5VEC_TAUIJ(I,ITAUIJ)),
     &                           N123TAUIJMAX)
               END DO
            END DO
         END DO
      END DO
C
      CHKSUMI2 = SUM(IQ_TAUIJ_TMP(1:NTAUIJ))
      CHKSUMJ2 = SUM(JQ_TAUIJ_TMP(1:NTAUIJ))
C
      DO I = 1,5
         CHKSUMN2(I) = SUM(N5VEC_TAUIJ(I,1:NTAUIJ))
         IF ( CHKSUMN1(I).NE.CHKSUMN2(I) ) TROUBLE = .TRUE.
      END DO
C
      IF ( TROUBLE ) THEN
         WRITE (6,*) 'CHECKSUM IQ ',CHKSUMI1,CHKSUMI2
         WRITE (6,*) 'CHECKSUM JQ ',CHKSUMJ1,CHKSUMJ2
         DO I = 1,5
            WRITE (6,*) 'CHECKSUM N5 ',CHKSUMN1(I),CHKSUMN2(I)
         END DO
      END IF
C
      IF ( ITAUIJ2.NE.NTAUIJ ) THEN
         WRITE (6,*) '<INIT_MOD_TAUIJ_TABLE>: ITAUIJ2 .NE. NTAUIJ'
         TROUBLE = .TRUE.
      END IF
C
      IF ( TROUBLE ) STOP 'trouble in <INIT_MOD_TAUIJ_TABLE>'
C
99001 FORMAT ('selected TAUIJ:',3I4,2x,2I3,2x,5I3,/)
      END
