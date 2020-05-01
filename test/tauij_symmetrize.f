C*==tauij_symmetrize.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_SYMMETRIZE(MCLU)
C   ********************************************************************
C   *                                                                  *
C   *   determine set of equivalent cluster site pairs within TAU_ij   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQCLU
      USE MOD_SYMMETRY,ONLY:DROT,ICL_IJQCLU,NCL_IJQCLU,ISYM_GEN_IJQCLU,
     &    SYMUNITARY,IQCLU_REP_IJQCLU,JQCLU_REP_IJQCLU
      USE MOD_ANGMOM,ONLY:NKKR_CLU,NXM,NKMMAX,WKM1,WKM2
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--TAUIJ_SYMMETRIZE16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_SYMMETRIZE')
      LOGICAL CHECK_SYMMETRISATION
      PARAMETER (CHECK_SYMMETRISATION=.FALSE.)
      REAL*8 TOL_SYMSUMRTR,THRESH_SYMSUMRTR
      PARAMETER (TOL_SYMSUMRTR=1D-8,THRESH_SYMSUMRTR=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 MCLU(NKKR_CLU,NKKR_CLU)
C
C Local variables
C
      CHARACTER*1 CNT
      INTEGER ICL,IQCLU,IQCLU_REP,IRELEFF,ISYM,JQCLU,JQCLU_REP,M,N,
     &        N_EQUIV_PAIRS,N_ISYM_EQ_1_DEVIATE,N_ISYM_EQ_1_SAME,N_SAME
      COMPLEX*16 MCLU_IJ(:,:),WXREP(:,:),WXREP_ORG(:,:)
      LOGICAL SAME
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MCLU_IJ,WXREP,WXREP_ORG
C
      CALL TRACK_INFO(ROUTINE)
C
      N = NXM
      M = NKMMAX
C
      IF ( IREL.LE.2 ) THEN
         IRELEFF = 1
      ELSE
         IRELEFF = 3
      END IF
C
      N_SAME = 0
      N_ISYM_EQ_1_SAME = 0
      N_ISYM_EQ_1_DEVIATE = 0
C
      ALLOCATE (WXREP(M,M),WXREP_ORG(M,M),MCLU_IJ(M,M))
C
C=======================================================================
C                  scan all classes of cluster site pairs
C=======================================================================
C
      LOOP_ICL:DO ICL = 1,NCL_IJQCLU
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK_SYMMETRISATION ) WRITE (6,99003) ICL
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         WXREP(1:N,1:N) = C0
C
C-----------------------------------------------------------------------
C    average inverse TAU-matrices for the equivalent pairs of a class
C-----------------------------------------------------------------------
         N_EQUIV_PAIRS = 0
C
         LOOP_IQCLU_1:DO IQCLU = 1,NQCLU
            LOOP_JQCLU_1:DO JQCLU = 1,NQCLU
C
               IF ( ICL.NE.ICL_IJQCLU(IQCLU,JQCLU) ) CYCLE LOOP_JQCLU_1
C
               N_EQUIV_PAIRS = N_EQUIV_PAIRS + 1
C
               CALL TAUIJ_GET_MCLUIJ(MCLU,IQCLU,JQCLU,MCLU_IJ)
C
               IF ( N_EQUIV_PAIRS.EQ.1 ) THEN
                  IQCLU_REP = IQCLU_REP_IJQCLU(IQCLU,JQCLU)
                  JQCLU_REP = JQCLU_REP_IJQCLU(IQCLU,JQCLU)
                  IF ( CHECK_SYMMETRISATION ) WXREP_ORG(1:N,1:N)
     &                 = MCLU_IJ(1:N,1:N)
               END IF
C
C-------- all pairs of a class have to have the same representative pair
               IF ( IQCLU_REP.NE.IQCLU_REP_IJQCLU(IQCLU,JQCLU) .OR. 
     &              JQCLU_REP.NE.JQCLU_REP_IJQCLU(IQCLU,JQCLU) )
     &              CALL STOP_MESSAGE(ROUTINE,
     &                                'IQCLU_REP .ne. IQCLU_REP_IJQCLU')
C
C.......................................................................
C                rotate MCLU_IJ   TO   representative site pair
C.......................................................................
C
               ISYM = ISYM_GEN_IJQCLU(IQCLU,JQCLU)
C
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C
               CALL ZGEMM('C',CNT,N,N,N,C1,DROT(1,1,ISYM),M,MCLU_IJ,M,
     &                    C0,WKM1,M)
C
               CALL ZGEMM('N','N',N,N,N,C1,WKM1,M,DROT(1,1,ISYM),M,C0,
     &                    WKM2,M)
C.......................................................................
C
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               IF ( CHECK_SYMMETRISATION ) THEN
                  WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),ICL,ISYM,
     &                            IQCLU_REP,JQCLU_REP,IQCLU,JQCLU
C
                  CALL CMATSTRUCT('original matrix B = MCLU_IJ before '
     &                            //'transfer to representative pair',
     &                            MCLU_IJ,N,M,IRELEFF,IRELEFF,1,
     &                            THRESH_SYMSUMRTR,6)
C
                  CALL CMATCMP(N,M,IRELEFF,
     &                         'WXREP - representative pair',WXREP_ORG,
     &                         'WKM2  - transferred pair   ',WKM2,
     &                         THRESH_SYMSUMRTR,TOL_SYMSUMRTR,SAME)
C
                  IF ( SAME ) THEN
                     N_SAME = N_SAME + 1
                     IF ( ISYM.EQ.1 )
     &                    N_ISYM_EQ_1_SAME = N_ISYM_EQ_1_SAME + 1
                  ELSE
                     IF ( ISYM.EQ.1 )
     &                    N_ISYM_EQ_1_DEVIATE = N_ISYM_EQ_1_DEVIATE + 1
                  END IF
               END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
               WXREP(1:N,1:N) = WXREP(1:N,1:N) + WKM2(1:N,1:N)
C
            END DO LOOP_JQCLU_1
         END DO LOOP_IQCLU_1
C
         WXREP(1:N,1:N) = WXREP(1:N,1:N)/DBLE(N_EQUIV_PAIRS)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         IF ( CHECK_SYMMETRISATION ) THEN
            CALL CMATSTRUCT('symmetrized representative matrix A',WXREP,
     &                      N,M,IRELEFF,IRELEFF,1,THRESH_SYMSUMRTR,6)
C
            CALL CMATCMP(N,M,IRELEFF,
     &                   'WXREP_ORG - original    representative pair',
     &                   WXREP_ORG,
     &                   'WXREP     - symmetrized representative pair',
     &                   WXREP,THRESH_SYMSUMRTR,TOL_SYMSUMRTR,SAME)
C
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C-----------------------------------------------------------------------
C  transfer representative matrix to all equivalent pairs of a class
C-----------------------------------------------------------------------
         LOOP_IQCLU_2:DO IQCLU = 1,NQCLU
            LOOP_JQCLU_2:DO JQCLU = 1,NQCLU
C
               IF ( ICL.NE.ICL_IJQCLU(IQCLU,JQCLU) ) CYCLE LOOP_JQCLU_2
C
C.......................................................................
C                rotate MCLU_IJ   FROM   representative site pair
C.......................................................................
C
               ISYM = ISYM_GEN_IJQCLU(IQCLU,JQCLU)
C
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C
               CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),M,WXREP,M,C0,
     &                    WKM1,M)
C
               CALL ZGEMM('N','C',N,N,N,C1,WKM1,M,DROT(1,1,ISYM),M,C0,
     &                    MCLU_IJ,M)
C.......................................................................
C
               CALL TAUIJ_PUT_MCLUIJ(MCLU,IQCLU,JQCLU,MCLU_IJ)
C
            END DO LOOP_JQCLU_2
         END DO LOOP_IQCLU_2
C
      END DO LOOP_ICL
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHECK_SYMMETRISATION ) WRITE (6,99002) NQCLU**2,NCL_IJQCLU,
     &     N_SAME,N_ISYM_EQ_1_SAME,N_ISYM_EQ_1_DEVIATE
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
99001 FORMAT (//,10X,'<',A,'>  calling <CMATCMP> for class ICL =',I4,
     &        '   ISYM =',I4,/,10X,'IQCLU_REP =',I4,'  JQCLU_REP =',I4,
     &        6X,'IQCLU =',I4,'  JQCLU =',I4,/,10X,
     &        'comparing original and symmetrized matrices')
99002 FORMAT (//,10X,'CHECK_SYMMETRISATION:',/,10X,
     &        'total number of site pairs           NQCLU*NQCLU: ',I6,/,
     &        10X,'inequivalent site pairs (classes)    NCL_IJQCLU:  ',
     &        I6,/,10X,
     &        'coinciding matrices                               ',I6,/,
     &        10X,'coincidences for ISYM = 1                         ',
     &        I6,/,10X,
     &        'deviations   for ISYM = 1                         ',I6,/)
99003 FORMAT (//,1X,79('C'),/,10X,
     &        'checking symmetry of matrices for class ICL =',I4,/,1X,
     &        79('C'),/)
      END
C*==tauij_get_mcluij.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_GET_MCLUIJ(MCLU,IQCLU,JQCLU,MCLU_IJ)
C   ********************************************************************
C   *                                                                  *
C   *   get sub block TAU_ij from cluster matrix MCLU                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:IXM0_QCLU,NKKR_CLU,NXM_QCLU,NKMMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--TAUIJ_GET_MCLUIJ250
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_GET_MCLUIJ')
C
C Dummy arguments
C
      INTEGER IQCLU,JQCLU
      COMPLEX*16 MCLU(NKKR_CLU,NKKR_CLU),MCLU_IJ(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,II0,J,JJ,JJ0,NI,NJ
C
C*** End of declarations rewritten by SPAG
C
      IF ( IPRINT.GE.5 ) CALL TRACK_INFO(ROUTINE)
C
      II0 = IXM0_QCLU(IQCLU)
      NI = NXM_QCLU(IQCLU)
C
      JJ0 = IXM0_QCLU(JQCLU)
      NJ = NXM_QCLU(JQCLU)
C
C-----------------------------------------------------------------------
      DO J = 1,NJ
         JJ = JJ0 + J
         DO I = 1,NI
            MCLU_IJ(I,J) = MCLU(II0+I,JJ)
         END DO
      END DO
C-----------------------------------------------------------------------
C
      END
C*==tauij_put_mcluij.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_PUT_MCLUIJ(MCLU,IQCLU,JQCLU,MCLU_IJ)
C   ********************************************************************
C   *                                                                  *
C   *   put sub block TAU_ij into cluster matrix MCLU                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:IXM0_QCLU,NKKR_CLU,NXM_QCLU,NKMMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--TAUIJ_PUT_MCLUIJ314
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_PUT_MCLUIJ')
C
C Dummy arguments
C
      INTEGER IQCLU,JQCLU
      COMPLEX*16 MCLU(NKKR_CLU,NKKR_CLU),MCLU_IJ(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,II0,J,JJ,JJ0,NI,NJ
C
C*** End of declarations rewritten by SPAG
C
      IF ( IPRINT.GE.5 ) CALL TRACK_INFO(ROUTINE)
C
      II0 = IXM0_QCLU(IQCLU)
      NI = NXM_QCLU(IQCLU)
C
      JJ0 = IXM0_QCLU(JQCLU)
      NJ = NXM_QCLU(JQCLU)
C
C-----------------------------------------------------------------------
      DO J = 1,NJ
         JJ = JJ0 + J
         DO I = 1,NI
            MCLU(II0+I,JJ) = MCLU_IJ(I,J)
         END DO
      END DO
C-----------------------------------------------------------------------
C
      END
