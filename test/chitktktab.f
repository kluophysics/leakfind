C*==chitktktab.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHITKTKTAB(NAB_CHI_QQ,NCOLROW,LAMCOLROW,KTKTK,
     &                      ERYDA_EQ_ERYDB,NTKTKLIN,NLIN41_CHI,
     &                      NLIN23_CHI)
C   ********************************************************************
C   *                                                                  *
C   *   set up the table of indices for the itegrals of the type       *
C   *                                                                  *
C   *     Int d^3k  TAU(k)[L1,L2] * TAU(k)[L3,L4]                      *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   the subroutine sets the array sizes NLINCHIMAX,NTKTKMAX        *
C   *   to start these are fixed according to the worst case           *
C   *   at the end the parameters and the corresponding arrays are     *
C   *   written to the transfer file  IOTMP                            *
C   *                                                                  *
C   ********************************************************************
C
C     NXM = NKM  IREL  = 3
C
      USE MOD_ANGMOM,ONLY:NL,NKM,NLIN,IKM2LIN,IKM1LIN,NKMMAX,NXM
      USE MOD_FILES,ONLY:IPRINT,IOTMP
      USE MOD_SITES,ONLY:NQMAX,ICPA,IQBOT_CHI,IQTOP_CHI
      USE MOD_CPA,ONLY:NCPA
      USE MOD_KSPACE,ONLY:NELMT,ITBZ,JTBZ,QTBZ
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--CHITKTKTAB29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHITKTKTAB')
      LOGICAL CHECK
      PARAMETER (CHECK=.TRUE.)
C
C Dummy arguments
C
      LOGICAL ERYDA_EQ_ERYDB
      INTEGER KTKTK,NLIN23_CHI,NLIN41_CHI,NTKTKLIN
      INTEGER LAMCOLROW(NKMMAX,NKMMAX,NQMAX),NAB_CHI_QQ(NQMAX,NQMAX),
     &        NCOLROW(NKMMAX,NQMAX)
C
C Local variables
C
      REAL*8 DMAT(NXM,NXM),MJ,MLAMBDA(NKMMAX),RMJ1,RMJ4
      LOGICAL FLAG(0:(1+(1+(1+NKM)*NKM)*NKM)*NKM)
      INTEGER I,IABCD,IA_ERR,ICDAB,ICOLROW,IDA,IDB,IDC,IDD,II,IKM,IKM1,
     &        IKM1_CHI_LIN(:),IKM2_CHI_LIN(:),IKM3_CHI_LIN(:),IKM4,
     &        IKM4_CHI_LIN(:),IL,ILAMA,ILAMB,ILAMC,ILAMD,IQ,ITT,J,JM05,
     &        JQ,JQTOP,K1,K2,K3,K4,L,L1,L4,LAMA,LAMB,LAMC,LAMD,LIN,
     &        LIN23,LIN41,LLAMBDA(NKMMAX),MJM05,N,NFLAG,NKM1,NKM2,NKM3,
     &        NLINCHIMAX,NTKTKLIM,NTKTKMAX,NTKTK_QQ(:,:)
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------- variables depending on NLINCHI
      ALLOCATABLE NTKTK_QQ
      ALLOCATABLE IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN
C
C***********************************************************************
      IF ( IREL.NE.3 ) THEN
         WRITE (6,99001) IREL
         STOP
      END IF
C***********************************************************************
C
      ALLOCATE (NTKTK_QQ(NQMAX,NQMAX))
C
      NLINCHIMAX = NKM*NKM
      ALLOCATE (IKM1_CHI_LIN(NLINCHIMAX))
      ALLOCATE (IKM2_CHI_LIN(NLINCHIMAX))
      ALLOCATE (IKM3_CHI_LIN(NLINCHIMAX))
      ALLOCATE (IKM4_CHI_LIN(NLINCHIMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 )
     &      CALL STOP_MESSAGE(ROUTINE,'allocate IKM3_CHI_LIN')
C
      DO I = 1,NLINCHIMAX
         IKM1_CHI_LIN(I) = 0
         IKM2_CHI_LIN(I) = 0
         IKM3_CHI_LIN(I) = 0
         IKM4_CHI_LIN(I) = 0
      END DO
C
      NTKTKMAX = 0
C
      WRITE (6,99002) NCPA,IQBOT_CHI,IQTOP_CHI,
     &                (ICPA(IQ),IQ=IQBOT_CHI,IQTOP_CHI)
      NKM1 = NKM
      NKM2 = NKM1*NKM
      NKM3 = NKM2*NKM
      NFLAG = (1+(1+(1+NKM)*NKM)*NKM)*NKM
C
C ------------------------ fix structure of projection matrices D and D~
C ------------------------------- assuming the same structure as for TAU
C
C NCOLROW(I)    :  number of non-0 columns (rows) for row (column) I
C LAMCOLROW(J,I):  index LAM of Jth non-0 column (row) in row (column) I
C
      DO IQ = 1,NQMAX
         DO IKM = 1,NKMMAX
            NCOLROW(IKM,IQ) = 0
         END DO
      END DO
C
      DO N = 1,NELMT
         IQ = QTBZ(N)
         I = ITBZ(N)
         J = JTBZ(N)
         NCOLROW(I,IQ) = NCOLROW(I,IQ) + 1
         LAMCOLROW(NCOLROW(I,IQ),I,IQ) = J
      END DO
C ----------------------------------------------------------------------
C                  no CPA and no special case: use structure of 1-matrix
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IF ( ICPA(IQ).EQ.0 ) THEN
            DO I = 1,NKM
               NCOLROW(I,IQ) = 1
               LAMCOLROW(NCOLROW(I,IQ),I,IQ) = I
            END DO
         END IF
      END DO
C
      IF ( CHECK .OR. IPRINT.GE.2 ) THEN
         WRITE (6,99003)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            WRITE (6,*) ' '
            DO I = 1,NKM
               WRITE (6,99004) IQ,I,NCOLROW(I,IQ),
     &                         (LAMCOLROW(N,I,IQ),N=1,NCOLROW(I,IQ))
            END DO
C
            DMAT(:,:) = 0D0
C
            DO I = 1,NXM
               DO ICOLROW = 1,NCOLROW(I,IQ)
                  DMAT(I,LAMCOLROW(ICOLROW,I,IQ)) = 1D0
               END DO
            END DO
C
            IF ( IREL.EQ.3 ) THEN
               CALL RMATSTRUCT('structure of D-matrix',DMAT,NXM,NXM,3,3,
     &                         0,1D-8,6)
            ELSE
               CALL RMATSTRUCT('structure of D-matrix',DMAT,NXM,NXM,1,1,
     &                         0,1D-8,6)
            END IF
            WRITE (6,*) ' '
C
         END DO
C
      END IF
C
C
C ======================================================================
C
      IKM = 0
      DO IL = 1,NL
         L = IL - 1
C
         DO JM05 = MAX(0,L-1),L
C
            DO MJM05 = -JM05 - 1,JM05
               MJ = DBLE(MJM05) + 0.5D0
C
               IKM = IKM + 1
               LLAMBDA(IKM) = IL - 1
               MLAMBDA(IKM) = MJ
C
            END DO
C
         END DO
C
      END DO
C
C ----------------------------------------------------------------------
C  selection rules scheme for perturbation  B
C                    |l-l'| = 0   |mj-mj'| = 0
C ----------------------------------------------------------------------
      NLIN23_CHI = NLIN
      IF ( NLIN.GT.NLINCHIMAX )
     &      CALL STOP_MESSAGE(ROUTINE,' NLIN > NLINCHIMAX (A1)')
      DO LIN23 = 1,NLIN23_CHI
         IKM2_CHI_LIN(LIN23) = IKM1LIN(LIN23)
         IKM3_CHI_LIN(LIN23) = IKM2LIN(LIN23)
      END DO
C
      IF ( KTKTK.EQ.0 ) THEN
C ----------------------------------------------------------------------
C  selection rules scheme for observable    A
C                    |l-l'| = 0   |mj-mj'| = 0
C ----------------------------------------------------------------------
         NLIN41_CHI = NLIN
         IF ( NLIN.GT.NLINCHIMAX )
     &         CALL STOP_MESSAGE(ROUTINE,' NLIN > NLINCHIMAX (A2)')
         DO LIN41 = 1,NLIN41_CHI
            IKM4_CHI_LIN(LIN41) = IKM1LIN(LIN41)
            IKM1_CHI_LIN(LIN41) = IKM2LIN(LIN41)
         END DO
C ----------------------------------------------------------------------
      ELSE IF ( KTKTK.EQ.1 ) THEN
C ----------------------------------------------------------------------
C  coupling scheme for non-spherical representation of observable   A
C  cubic system with         |l-l'| = 0,2,4,..  |mj-mj'| = 0,4,8,....
C ----------------------------------------------------------------------
         NLIN41_CHI = 0
         DO IKM4 = 1,NKM
C
            L4 = LLAMBDA(IKM4)
            RMJ4 = MLAMBDA(IKM4)
C
            DO IKM1 = 1,NKM
C
               L1 = LLAMBDA(IKM1)
               RMJ1 = MLAMBDA(IKM1)
C
               IF ( MOD(L4-L1,2).EQ.0 ) THEN
                  IF ( ABS(MOD(NINT(RMJ4-RMJ1),4)).EQ.0 ) THEN
                     NLIN41_CHI = NLIN41_CHI + 1
                     IF ( NLIN41_CHI.GT.NLINCHIMAX )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                    ' NLIN41_CHI > NLINCHIMAX (B)')
                     IKM4_CHI_LIN(NLIN41_CHI) = IKM4
                     IKM1_CHI_LIN(NLIN41_CHI) = IKM1
                  END IF
               END IF
            END DO
         END DO
C ----------------------------------------------------------------------
      ELSE IF ( KTKTK.EQ.2 ) THEN
C ----------------------------------------------------------------------
C  coupling scheme for the full susceptibility tensor
C                    |l-l'| = 0   |mj-mj'| = 0, 1
C  to be applied to the operator A term 41
C  and              perturbation B term 23
C ----------------------------------------------------------------------
         NLIN41_CHI = 0
         DO IKM4 = 1,NKM
C
            L4 = LLAMBDA(IKM4)
            RMJ4 = MLAMBDA(IKM4)
C
            DO IKM1 = 1,NKM
C
               L1 = LLAMBDA(IKM1)
               RMJ1 = MLAMBDA(IKM1)
C
               IF ( (L4-L1).EQ.0 ) THEN
                  IF ( ABS(RMJ4-RMJ1).LE.1 ) THEN
                     NLIN41_CHI = NLIN41_CHI + 1
                     IF ( NLIN41_CHI.GT.NLINCHIMAX )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                    ' NLIN41_CHI > NLINCHIMAX (C)')
                     IKM4_CHI_LIN(NLIN41_CHI) = IKM4
                     IKM1_CHI_LIN(NLIN41_CHI) = IKM1
                     IKM2_CHI_LIN(NLIN41_CHI) = IKM4
                     IKM3_CHI_LIN(NLIN41_CHI) = IKM1
                  END IF
               END IF
            END DO
         END DO
C
         NLIN23_CHI = NLIN41_CHI
C ----------------------------------------------------------------------
      ELSE IF ( KTKTK.EQ.3 ) THEN
C ----------------------------------------------------------------------
C  coupling scheme for the conductivity tensor
C                    |l-l'| = 1   |mj-mj'| = 0, 1
C  to be applied to the operator A term 41
C  and              perturbation B term 23
C ----------------------------------------------------------------------
         NLIN41_CHI = 0
         DO IKM4 = 1,NKM
C
            L4 = LLAMBDA(IKM4)
            RMJ4 = MLAMBDA(IKM4)
C
            DO IKM1 = 1,NKM
C
               L1 = LLAMBDA(IKM1)
               RMJ1 = MLAMBDA(IKM1)
C
               IF ( ABS(L4-L1).EQ.1 ) THEN
                  IF ( ABS(RMJ4-RMJ1).LE.1 ) THEN
                     IF ( NLIN41_CHI.GT.NLINCHIMAX )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                    'NLIN41_CHI > NLINCHIMAX (D)')
                     NLIN41_CHI = NLIN41_CHI + 1
                     IKM4_CHI_LIN(NLIN41_CHI) = IKM4
                     IKM1_CHI_LIN(NLIN41_CHI) = IKM1
                     IKM2_CHI_LIN(NLIN41_CHI) = IKM4
                     IKM3_CHI_LIN(NLIN41_CHI) = IKM1
                  END IF
               END IF
            END DO
         END DO
C
         NLIN23_CHI = NLIN41_CHI
C ----------------------------------------------------------------------
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'check setting of KTKTK (1)')
      END IF
C
      IF ( KTKTK.EQ.0 ) THEN
         WRITE (6,99005)
      ELSE IF ( KTKTK.EQ.1 ) THEN
         WRITE (6,99006)
      ELSE IF ( KTKTK.EQ.2 ) THEN
         WRITE (6,99007)
      ELSE IF ( KTKTK.EQ.3 ) THEN
         WRITE (6,99008)
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'check setting of KTKTK (2)')
      END IF
C
      IF ( CHECK .OR. IPRINT.GT.0 ) THEN
C
         DMAT(:,:) = 0D0
         DO LIN = 1,NLIN41_CHI
            DMAT(IKM4_CHI_LIN(LIN),IKM1_CHI_LIN(LIN)) = 1
         END DO
C
         CALL RMATSTRUCT('coupling scheme for term 41',DMAT,NKM,NKM,3,3,
     &                   0,1D-8,6)
         WRITE (6,*) ' '
C
         DMAT(:,:) = 0D0
         DO LIN = 1,NLIN23_CHI
            DMAT(IKM2_CHI_LIN(LIN),IKM3_CHI_LIN(LIN)) = 1
         END DO
C
         CALL RMATSTRUCT('coupling scheme for term 23',DMAT,NKM,NKM,3,3,
     &                   0,1D-8,6)
         WRITE (6,*) ' '
C
      END IF
C
      ITT = 0
C
C=======================================================================
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IF ( ERYDA_EQ_ERYDB ) THEN
            JQTOP = IQ
         ELSE
            JQTOP = IQTOP_CHI
         END IF
C
         DO JQ = IQBOT_CHI,JQTOP
C
            NTKTK_QQ(IQ,JQ) = 0
C
            DO I = 0,NFLAG
               FLAG(I) = .FALSE.
            END DO
C
            NAB_CHI_QQ(IQ,JQ) = 0
            NTKTKLIM = 0
C
C---------------------------------------------------------------------41
            DO LIN41 = 1,NLIN41_CHI
C
               K4 = IKM4_CHI_LIN(LIN41)
               K1 = IKM1_CHI_LIN(LIN41)
C
C---------------------------------------------------------------------23
               DO LIN23 = 1,NLIN23_CHI
C
                  K2 = IKM2_CHI_LIN(LIN23)
                  K3 = IKM3_CHI_LIN(LIN23)
C
                  NAB_CHI_QQ(IQ,JQ) = NAB_CHI_QQ(IQ,JQ) + 1
C
C-----------------------------------------------------------------------
                  DO IDA = 1,NCOLROW(K1,IQ)
                     LAMA = LAMCOLROW(IDA,K1,IQ)
                     ILAMA = (LAMA-1)*NKM3
C
                     DO IDB = 1,NCOLROW(K2,JQ)
                        LAMB = LAMCOLROW(IDB,K2,JQ)
                        ILAMB = ILAMA + (LAMB-1)*NKM2
C
                        DO IDC = 1,NCOLROW(K3,JQ)
                           LAMC = LAMCOLROW(IDC,K3,JQ)
                           ILAMC = ILAMB + (LAMC-1)*NKM
C
                           DO IDD = 1,NCOLROW(K4,IQ)
                              LAMD = LAMCOLROW(IDD,K4,IQ)
                              IABCD = ILAMC + (LAMD-1)
C
                              NTKTKLIM = NTKTKLIM + 1
C
                              FLAG(IABCD) = .TRUE.
C
                           END DO
                        END DO
                     END DO
                  END DO
C
               END DO
C---------------------------------------------------------------------23
            END DO
C---------------------------------------------------------------------41
C
            DO I = 0,NFLAG
               IF ( FLAG(I) ) THEN
                  ILAMA = I/NKM3
                  II = I - ILAMA*NKM3
                  ILAMB = II/NKM2
                  II = II - ILAMB*NKM2
                  ILAMC = II/NKM
                  ILAMD = II - ILAMC*NKM
C
                  IF ( ERYDA_EQ_ERYDB .AND. IQ.EQ.JQ ) THEN
                     ICDAB = ILAMC*NKM3 + ILAMD*NKM2 + ILAMA*NKM + ILAMB
                     FLAG(ICDAB) = .FALSE.
                     CYCLE
                  END IF
C
                  NTKTK_QQ(IQ,JQ) = NTKTK_QQ(IQ,JQ) + 1
C
                  ITT = ITT + 1
C
               END IF
            END DO
C
            WRITE (6,99009) IQ,JQ,NAB_CHI_QQ(IQ,JQ),NLIN**2,NTKTKLIM,
     &                      NTKTK_QQ(IQ,JQ)
            NTKTKMAX = MAX(NTKTKMAX,NTKTK_QQ(IQ,JQ))
         END DO
      END DO
C=======================================================================
C
      NTKTKLIN = ITT
C
      WRITE (6,99010) NTKTKMAX
C
C=======================================================================
C       write tables to temporary transfer file - read in <CHIDRIVE>
C ======================================================================
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      NLINCHIMAX = MAX(NLIN41_CHI,NLIN23_CHI)
      WRITE (IOTMP) NLINCHIMAX,NTKTKMAX
      WRITE (IOTMP) (IKM1_CHI_LIN(I),IKM2_CHI_LIN(I),IKM3_CHI_LIN(I),
     &              IKM4_CHI_LIN(I),I=1,NLINCHIMAX)
C
C=======================================================================
C
      DEALLOCATE (IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN)
C
99001 FORMAT (//,10X,'<CHITKTKTAB> called for IREL =',I3,//,10X,
     &        'the routine should be used exclusivly for',
     &        ' the fully relativisitc mode (IREL=3)',//)
99002 FORMAT (//,1X,79('*'),/,29X,'<CHITKTKTAB>',/,1X,79('*'),//,5X,
     &        'setting index tables for',/,5X,
     &        'BZ-integral   Int d^3k  TAU(k)[L1,L2] * TAU(k)[L3,L4]',
     &        //,5X,'NCPA                ',I12,/,5X,//,5X,
     &        'for sites IQ =      ',I12,' ...',I5,/,5X,
     &        'ICPA(IQ)            ',9X,20I3,/)
99003 FORMAT (/,5X,'structure of projection matrices D and D~ ')
99004 FORMAT (5X,'site',I3,'  row',I3,'  with',I3,'   columns:',3X,30I3)
99005 FORMAT (/,5X,'standard settings for CHI-calculation ',/,5X,
     &        'selection rules for observable   A (term 41)',
     &        '  |l-l''| = 0  |mj-mj''| = 0',/,5X,
     &        'selection rules for perturbation B (term 23)',
     &        '  |l-l''| = 0  |mj-mj''| = 0',/)
99006 FORMAT (/,5X,'calculate non-spherical CHI(r)  ',/,5X,
     &        'selection rules for observable   A (term 41)',
     &        '  |l-l''| = 0,2,4,... |mj-mj''| = 0,4,8,...',' (cubic)',
     &        5X,'selection rules for perturbation B (term 23)',
     &        '  |l-l''| = 0  |mj-mj''| = 0',/)
99007 FORMAT (/,5X,'calculate full CHI-tensor  ',/,5X,
     &        'selection rules for observable   A (term 41)',
     &        '  |l-l''| = 0  |mj-mj''| = 0, 1',/,5X,
     &        'selection rules for perturbation B (term 23)',
     &        '  |l-l''| = 0  |mj-mj''| = 0, 1',/)
99008 FORMAT (/,5X,'calculate full conductivity tensor  ',/,5X,
     &        'selection rules for observable   A (term 41)',
     &        '  |l-l''| = 1  |mj-mj''| = 0, 1',/,5X,
     &        'selection rules for perturbation B (term 23)',
     &        '  |l-l''| = 1  |mj-mj''| = 0, 1',/)
99009 FORMAT (5X,'sites',31X,'IQ =',I3,'    JQ =',I3,/,5X,
     &        'number of LAMA-LAMB  combinations   NAB_QQ   =',I12,/,5X,
     &        '                                    NLIN**2  =',I12,/,5X,
     &        'number of D*D*D*D*TT products       NTKTKLIM =',I12,/,5X,
     &        'number of different TT-terms        NTKTK    =',I12,/,5X,
     &        /)
99010 FORMAT (5X,'array size                          NTKTKMAX =',I12,/)
      END
