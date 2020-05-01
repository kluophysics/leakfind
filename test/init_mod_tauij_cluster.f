C*==init_mod_tauij_cluster.f    processed by SPAG 6.70Rc at 09:27 on 13 Oct 2014
      SUBROUTINE INIT_MOD_TAUIJ_CLUSTER
C   ********************************************************************
C   *                                                                  *
C   *   determine the number of TAU_ij to be calculated for a          *
C   *   given embedded cluster                                         *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQCLU,IQ_QCLU,N5VEC_QCLU
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_TAUIJ,ONLY:ITAUIJ_LOOP,NLOOP_TAUIJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_TAUIJ_CLUSTER')
C
C Local variables
C
      INTEGER ILOOP,IQ,IQCLU,IQ_TAUIJ_TMP(:),ITAUIJ_ORIG,ITAUIJ_SORT,
     &        ITAUIJ_SORT_ORIG(:),JQ,JQCLU,JQ_TAUIJ_TMP(:),
     &        N5VEC_IQCLU(5),N5VEC_JQCLU(5),N5VEC_TAUIJ_TMP(:,:),
     &        NTAUIJMAX_TMP,NTAUIJ_TMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE N5VEC_TAUIJ_TMP,IQ_TAUIJ_TMP,JQ_TAUIJ_TMP
      ALLOCATABLE ITAUIJ_SORT_ORIG
C
      CALL TRACK_INFO(ROUTINE)
C
C-----------------------------------------------------------------------
C                 allocate temporary work storage
C-----------------------------------------------------------------------
      NTAUIJMAX_TMP = NQCLU*NQCLU
C
      ALLOCATE (N5VEC_TAUIJ_TMP(5,NTAUIJMAX_TMP))
      ALLOCATE (IQ_TAUIJ_TMP(NTAUIJMAX_TMP),JQ_TAUIJ_TMP(NTAUIJMAX_TMP))
C
      NLOOP_TAUIJ = NQCLU*NQCLU
      ALLOCATE (ITAUIJ_LOOP(NLOOP_TAUIJ))
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
      NTAUIJ_TMP = 0
      ILOOP = 0
C
      DO JQCLU = 1,NQCLU
         JQ = IQ_QCLU(JQCLU)
         N5VEC_JQCLU(1:5) = N5VEC_QCLU(1:5,JQCLU)
C
         DO IQCLU = 1,NQCLU
            IQ = IQ_QCLU(IQCLU)
            N5VEC_IQCLU(1:5) = N5VEC_QCLU(1:5,IQCLU)
C
            ILOOP = ILOOP + 1
C
            CALL INIT_MOD_TAUIJ_LOOP(ITAUIJ_LOOP,ILOOP,NLOOP_TAUIJ,IQ,
     &                               N5VEC_IQCLU,JQ,N5VEC_JQCLU,
     &                               NTAUIJ_TMP,IQ_TAUIJ_TMP,
     &                               JQ_TAUIJ_TMP,N5VEC_TAUIJ_TMP,
     &                               NTAUIJMAX_TMP)
C
            IF ( IPRINT.GT.0 ) WRITE (6,'(A,2(2I3,2x),I4)')
     &                                 ' cluster sites ',IQCLU,IQ,JQCLU,
     &                                JQ,ITAUIJ_LOOP(ILOOP)
         END DO
      END DO
C
      WRITE (6,99001) NTAUIJ_TMP,NQCLU
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
      DO JQCLU = 1,NQCLU
         DO IQCLU = 1,NQCLU
            ILOOP = ILOOP + 1
C
            ITAUIJ_ORIG = ITAUIJ_LOOP(ILOOP)
            ITAUIJ_SORT = ITAUIJ_SORT_ORIG(ITAUIJ_ORIG)
            ITAUIJ_LOOP(ILOOP) = ITAUIJ_SORT
C
         END DO
      END DO
C
      DEALLOCATE (ITAUIJ_SORT_ORIG)
      DEALLOCATE (N5VEC_TAUIJ_TMP,IQ_TAUIJ_TMP,JQ_TAUIJ_TMP)
C
99001 FORMAT (/,1X,79('*'),/,28X,'<INIT_MOD_TAUIJ_CLUSTER>',/,1X,79('*')
     &        ,//,5X,'number of TAU(I,J)''s to be calculated:',I5,
     &        '   for  NQCLU =',I5,/)
      END
