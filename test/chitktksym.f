C*==chitktksym.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHITKTKSYM(SYMACCEPTED,NSYMACCEPTED)
C   ********************************************************************
C   *                                                                  *
C   *   set up the table of indices for the itegrals of the type       *
C   *                                                                  *
C   *     Int d^3k  TAU(k,z1)[L1,L2;q1,q2] * TAU(k,z2)[L3,L4;q2,q1]    *
C   *                                                                  *
C   *  assuming calculation of the symmetry coefficients for           *
C   *  ALL NON-vanishing matrix elements                               *
C   *                                                                  *
C   *  ERYDA_EQ_ERYDB  = TRUE make use of the fact that TAU(k,z1) and  *
C   *                    TAU(k,z2) may be interchanged                 *
C   *            FALSE   this applies for SIGMA - calculation          *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  the tables ITT* and JTT* are passed via module MOD_SYMMETRY     *
C   *  the array sizes  ITTMAX and JTTMAX are set to the actual        *
C   *  length of the tables                                            *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  NOTE: the linear indices JTT* refer to the elements of the      *
C   *        KKR matrix. For TB calculations this is set up only       *
C   *        for the I-zone. therefore there is an offset in the       *
C   *        indices for the Q-dependent blocks by (IQBOT_TB-1)        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,ITTA,ITTB,ITTC,ITTD,ITTQ1,
     &    ITTQ2,NTTJ,JTTMAX,ITTMAX,NTKTKLIN,NTKTKMAX,JTT1,JTT2,JTTX,WTTJ
      USE MOD_FILES,ONLY:IOTMP,IOTMP1,IOTMP2
      USE MOD_SYMMETRY,ONLY:NSYM,SYMUNITARY,IQORGQP,DROT,NSYMMAX
      USE MOD_CALCMODE,ONLY:IREL,KMROT,KKRMODE
      USE MOD_ANGMOM,ONLY:NL,NLM,NKM,NXM,NKMQ,NLMQ,NKKR,IND0Q
      USE MOD_SITES,ONLY:NQ,IQBOT_TB,IQBOT_CHI,IQTOP_CHI
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      IMPLICIT NONE
C*--CHITKTKSYM39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHITKTKSYM')
      REAL*8 TOLZ
      PARAMETER (TOLZ=1D-8)
      LOGICAL NO_IMAG_WTT
      PARAMETER (NO_IMAG_WTT=.FALSE.)
C
C Dummy arguments
C
      INTEGER NSYMACCEPTED
      LOGICAL SYMACCEPTED(NSYMMAX)
C
C Local variables
C
      COMPLEX*16 DD(:,:,:),F1,F2,WTT,WTTJ_JTT,WW(:,:,:,:,:,:),ZARG
      CHARACTER*80 FNAMTTTAB
      INTEGER HLIM(:),I,I1,IA,IAB,IA_ERR,IB,IC,ICD,ICTOP,ID,ID_SEND,
     &        ILOOP,IND0OFFSET,IND0QX(:),IPROC,IQ,IQ01,IQ02,IQ1,IQ2,
     &        ISYM,ITST1,ITST2,ITST3,ITT,ITT1_PROC(:),IWORK(:),J,J1,J12,
     &        J1H,J1L,J2,J2H,J2L,J3,J34,J3H,J3L,J4,J4H,J4L,JADD(:,:,:),
     &        JBDD(:,:,:),JQ,JTT,JTT1_JTT,JTT1_PROC(:),JTT2_JTT,
     &        JTTX_JTT,K,L,LFN,LLIM(:),M,MJM05,N,NDD(:,:),NI,NJ,NN1,NN2,
     &        NQ2,NSCAN,NTERMSUM,NTERMSUMTOP,NTKTKQQ(:,:),N_ITT_PROC(:),
     &        N_JTT_PROC(:)
      LOGICAL LARGEZ
      LOGICAL LTST,MAGNETIC,SYM_CONNECTED(:,:)
      REAL*8 RJ,TOLTT,WA,WI,WNORM,WR
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
      ALLOCATABLE LLIM,HLIM
      ALLOCATABLE NDD,DD,WW,JADD,JBDD,IND0QX
      ALLOCATABLE SYM_CONNECTED,NTKTKQQ
      ALLOCATABLE ITT1_PROC,JTT1_PROC,N_ITT_PROC,N_JTT_PROC,IWORK
C
C-----------------------------------------------------------------------
C
      LARGEZ(ZARG) = ABS(DREAL(ZARG)) + ABS(DIMAG(ZARG)).GT.TOLZ
C
C-----------------------------------------------------------------------
C
      ALLOCATE (NTKTKQQ(NQ,NQ),IND0QX(NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NTKTKQQ')
C
C-----------------------------------------------------------------------
C
      IF ( IREL.LE.2 .AND. NXM.NE.NLM )
     &      CALL STOP_MESSAGE(ROUTINE,'IREL <= 2  and  NXM <> NLM')
      IF ( IREL.EQ.3 .AND. NXM.NE.NKM )
     &      CALL STOP_MESSAGE(ROUTINE,'IREL = 3   and  NXM <> NKM')
C
      ALLOCATE (NDD(NXM*NXM,NSYMMAX),LLIM(NXM),HLIM(NXM))
C
C-----------------------------------------------------------------------
C
      ALLOCATE (SYM_CONNECTED(NQ,NQ))
      SYM_CONNECTED(:,:) = .FALSE.
      DO IQ1 = IQBOT_CHI,IQTOP_CHI
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               IQ01 = IQORGQP(ISYM,IQ1)
               SYM_CONNECTED(IQ1,IQ01) = .TRUE.
               SYM_CONNECTED(IQ01,IQ1) = .TRUE.
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      ALLOCATE (DD(2*NL*2*NL,NXM*NXM,NSYMMAX))
      ALLOCATE (JADD(2*NL*2*NL,NXM*NXM,NSYMMAX))
      ALLOCATE (JBDD(2*NL*2*NL,NXM*NXM,NSYMMAX))
      IQ1 = IQBOT_CHI
      IQ2 = IQTOP_CHI
      ALLOCATE (WW(NXM,NXM,NXM,NXM,IQ1:IQ2,IQ1:IQ2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WW')
C
      WRITE (6,99003) ERYDA_EQ_ERYDB
C
C***********************************************************************
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
         IND0QX(1:NQ) = IND0Q(1:NQ)
C
C***********************************************************************
      ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
         IND0OFFSET = (IQBOT_TB-1)*NXM
C
         IND0QX(1:NQ) = 999999
         IND0QX(IQBOT_CHI:IQTOP_CHI) = IND0Q(IQBOT_CHI:IQTOP_CHI)
     &                                 - IND0OFFSET
C
C***********************************************************************
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
      END IF
C***********************************************************************
C
C ======================================================================
C           try to read symmetry tables from file    TKTK_tables
C                create table if parameters are inconsistent
C ======================================================================
C
      FNAMTTTAB = 'TKTK_tables'
      LFN = 11
      OPEN (IOTMP,FILE=FNAMTTTAB(1:LFN),FORM='UNFORMATTED')
C
      READ (IOTMP,END=100,ERR=100) LTST,ITST1,ITST2,ITST3
      IF ( .NOT.(LTST .NEQV. ERYDA_EQ_ERYDB) ) THEN
         IF ( ITST1.EQ.IQBOT_CHI .AND. ITST2.EQ.IQTOP_CHI ) THEN
            IF ( ITST3.EQ.NL ) THEN
               READ (IOTMP,END=100,ERR=100) ITST1,ITST2
               IF ( ITST1.EQ.NXM ) THEN
                  IF ( ITST2.EQ.NSYMACCEPTED ) THEN
                     READ (IOTMP,END=100,ERR=100) ITST1,ITST2
                     IF ( ITST1.EQ.NKKR ) THEN
                        IF ( ITST2.EQ.NSYM ) THEN
C
                           DO ISYM = 1,NSYM
                              READ (IOTMP,END=100,ERR=100) LTST
                              IF ( LTST .NEQV. SYMACCEPTED(ISYM) )
     &                             GOTO 100
                              READ (IOTMP,END=100,ERR=100) LTST
                              IF ( LTST .NEQV. SYMUNITARY(ISYM) )
     &                             GOTO 100
                           END DO
C
                           READ (IOTMP,END=100,ERR=100) ITTMAX,JTTMAX,
     &                           NSCAN,NTERMSUMTOP,
     &                           ((NTKTKQQ(IQ,JQ),JQ=IQBOT_CHI,
     &                           IQTOP_CHI),IQ=IQBOT_CHI,IQTOP_CHI)
C
                           NTKTKLIN = ITTMAX
                           ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX))
                           ALLOCATE (ITTC(ITTMAX),ITTD(ITTMAX))
                           ALLOCATE (ITTQ1(ITTMAX),ITTQ2(ITTMAX))
                           ALLOCATE (NTTJ(ITTMAX),STAT=IA_ERR)
                           IF ( IA_ERR.NE.0 )
     &                           CALL STOP_MESSAGE(ROUTINE,'ALLOC: ITTD'
     &                          )
C
                           DO I = 1,NTKTKLIN
                              READ (IOTMP) ITTA(I),ITTB(I),ITTC(I),
     &                              ITTD(I),ITTQ1(I),ITTQ2(I),NTTJ(I)
                           END DO
C
                           ALLOCATE (WTTJ(JTTMAX),JTT1(JTTMAX))
                           ALLOCATE (JTT2(JTTMAX),JTTX(JTTMAX))
                           DO J = 1,JTTMAX
                              READ (IOTMP,END=100,ERR=100) JTT1(J),
     &                              JTT2(J),JTTX(J),WTTJ(J)
                           END DO
C
C-------------------------------------- TKTK tables successfully read in
C
                           CLOSE (IOTMP)
C
                           WRITE (6,99001) FNAMTTTAB(1:LFN)
C
                           WRITE (6,99005) (NQ**2)*(NXM**4),NSCAN,
     &                            NTKTKLIN,ITTMAX
                           DO IQ = IQBOT_CHI,IQTOP_CHI
                              WRITE (6,99006) IQ,IQTOP_CHI,
     &                               (NTKTKQQ(IQ,JQ),JQ=IQBOT_CHI,
     &                               IQTOP_CHI)
                              DO JQ = IQBOT_CHI,IQTOP_CHI
                                 NTKTKMAX = MAX(NTKTKMAX,NTKTKQQ(IQ,JQ))
                              END DO
                           END DO
                           WRITE (6,99007) NTKTKMAX,JTTMAX,NTERMSUMTOP
C
                           RETURN
C
C-------------------------------------- TKTK tables successfully read in
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C ======================================================================
C
C
C ======================================================================
C              TKTK tables not available or not consistent
C                             create tables
C ======================================================================
C
C ----------------------------------------------------------------------
C        rewind temporary file to store  JTT1,JTT2,JTTX, WTTJ
C ----------------------------------------------------------------------
C
C
 100  CONTINUE
      REWIND IOTMP
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP1)
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP2)
C
C ----------------------------------------------------------------------
C
      MAGNETIC = .FALSE.
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
            IF ( .NOT.SYMUNITARY(ISYM) ) MAGNETIC = .TRUE.
         END IF
      END DO
C
      WRITE (6,99002) FNAMTTTAB(1:LFN)
C
C ----------------------------------------------------------------------
C        find the bounding indices  LLIM  and  HLIM  for IKM-loops
C ----------------------------------------------------------------------
      IF ( IREL.LE.2 ) THEN
         I = 0
         DO L = 0,(NL-1)
            DO M = -L,L
               I = I + 1
               LLIM(I) = L**2 + 1
               HLIM(I) = L**2 + 2*L + 1
            END DO
         END DO
      ELSE
         I = 0
         DO K = 1,(2*NL-1)
            L = K/2
            IF ( MOD(K,2).EQ.1 ) THEN
               RJ = L + 0.5D0
            ELSE
               RJ = L - 0.5D0
            END IF
            DO MJM05 = NINT(-RJ-0.5D0),NINT(RJ-0.5D0)
               I = I + 1
               LLIM(I) = NINT(L*2*(RJ+0.5D0)+1)
               HLIM(I) = NINT(L*2*(RJ+0.5D0)+2*RJ+1)
            END DO
         END DO
      END IF
      WRITE (6,99004) 'low   ',(LLIM(I),I=1,NXM)
      WRITE (6,99004) 'high  ',(HLIM(I),I=1,NXM)
C
C ======================================================================
C   store the products    D * D+  of 2 rotation matrices for later use
C
C ================================================================ IA IB
      DO IA = 1,NXM
         IAB = (IA-1)*NXM
         DO IB = 1,NXM
            IAB = IAB + 1
C ======================================================================
            DO ISYM = 1,NSYM
               IF ( SYMACCEPTED(ISYM) ) THEN
C
                  N = 0
C ================================================================ J1 J2
                  DO J1 = LLIM(IA),HLIM(IA)
                     F1 = DROT(IA,J1,ISYM)
                     IF ( LARGEZ(F1) ) THEN
                        DO J2 = LLIM(IB),HLIM(IB)
                           F2 = DROT(IB,J2,ISYM)
                           IF ( LARGEZ(F2) ) THEN
C ----------------------------------------------------------------------
                              N = N + 1
                              DD(N,IAB,ISYM) = F1*DCONJG(F2)
                              JADD(N,IAB,ISYM) = J1
                              JBDD(N,IAB,ISYM) = J2
C ----------------------------------------------------------------------
                           END IF
                        END DO
                     END IF
                  END DO
C ================================================================ J1 J2
C
               END IF
               NDD(IAB,ISYM) = N
            END DO
C ================================================================= ISYM
         END DO
      END DO
C ================================================================ IA IB
C
      WNORM = 1D0/DBLE(NSYMACCEPTED)
      TOLTT = 1D-6
C
      NSCAN = 0
      NTERMSUMTOP = 0
C
      ITT = 0
      JTT = 0
C
      ILOOP = -1
C
C ======================================================================
C  determine symmetry coefficients for ALL NON-VANISHING matrix elements
C ======================================================================
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ1 = IQBOT_CHI,IQTOP_CHI
         IF ( IREL.EQ.2 ) THEN
            NN1 = NLMQ(IQ1)
         ELSE
            NN1 = NKMQ(IQ1)
         END IF
C
         IF ( ERYDA_EQ_ERYDB ) THEN
            NQ2 = IQ1
         ELSE
            NQ2 = IQTOP_CHI
         END IF
C
         DO IQ2 = IQBOT_CHI,NQ2
            IF ( IREL.EQ.2 ) THEN
               NN2 = NLMQ(IQ2)
            ELSE
               NN2 = NKMQ(IQ2)
            END IF
            NTKTKQQ(IQ1,IQ2) = ITT
C ================================================================ IA IB
            DO IA = 1,NN1
               IAB = (IA-1)*NXM
               DO IB = 1,NN2
                  IAB = IAB + 1
C ================================================================ IC ID
                  IF ( (IQ1.EQ.IQ2) .AND. ERYDA_EQ_ERYDB ) THEN
                     ICTOP = IA
                  ELSE
                     ICTOP = NN2
                  END IF
                  DO IC = 1,ICTOP
                     ICD = (IC-1)*NXM
                     DO ID = 1,NN1
                        ICD = ICD + 1
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
                        IF ( MPI ) THEN
                           ILOOP = ILOOP + 1
                           IF ( MOD(ILOOP,NPROCS).NE.MPI_ID ) CYCLE
                        END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                        IF ( MAGNETIC ) THEN
                           J4L = MIN(LLIM(ID),LLIM(IC))
                           J3L = MIN(LLIM(ID),LLIM(IC))
                           J2L = MIN(LLIM(IA),LLIM(IB))
                           J1L = MIN(LLIM(IA),LLIM(IB))
                           J4H = MAX(HLIM(ID),HLIM(IC))
                           J3H = MAX(HLIM(ID),HLIM(IC))
                           J2H = MAX(HLIM(IA),HLIM(IB))
                           J1H = MAX(HLIM(IA),HLIM(IB))
                        ELSE
                           J4L = LLIM(ID)
                           J3L = LLIM(IC)
                           J2L = LLIM(IB)
                           J1L = LLIM(IA)
                           J4H = HLIM(ID)
                           J3H = HLIM(IC)
                           J2H = HLIM(IB)
                           J1H = HLIM(IA)
                        END IF
C ======================================================================
C                      scan all symmetry operations
C                and sum the products    D * D+ * D * D+
C ======================================================================
                        DO IQ02 = IQBOT_CHI,IQTOP_CHI
                           DO IQ01 = IQBOT_CHI,IQTOP_CHI
C
                              IF ( .NOT.SYM_CONNECTED(IQ1,IQ01) .AND. 
     &                             .NOT.SYM_CONNECTED(IQ2,IQ01) ) CYCLE
                              IF ( .NOT.SYM_CONNECTED(IQ1,IQ02) .AND. 
     &                             .NOT.SYM_CONNECTED(IQ2,IQ02) ) CYCLE
C
                              DO J4 = J4L,J4H
                                 DO J3 = J3L,J3H
                                    DO J2 = J2L,J2H
                                       DO J1 = J1L,J1H
                                         WW(J1,J2,J3,J4,IQ01,IQ02) = C0
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
C
C ================================================================= ISYM
                        DO ISYM = 1,NSYM
                           IF ( SYMACCEPTED(ISYM) ) THEN
                              IQ01 = IQORGQP(ISYM,IQ1)
                              IQ02 = IQORGQP(ISYM,IQ2)
C
                              I = 0
C ================================================================ J1 J2
                              DO J12 = 1,NDD(IAB,ISYM)
                                 J1 = JADD(J12,IAB,ISYM)
                                 J2 = JBDD(J12,IAB,ISYM)
C ================================================================ J3 J4
                                 DO J34 = 1,NDD(ICD,ISYM)
                                    J3 = JADD(J34,ICD,ISYM)
                                    J4 = JBDD(J34,ICD,ISYM)
C
                                    IF ( SYMUNITARY(ISYM) ) THEN
C
                                       WW(J1,J2,J3,J4,IQ01,IQ02)
     &                                    = WW(J1,J2,J3,J4,IQ01,IQ02)
     &                                    + DD(J12,IAB,ISYM)
     &                                    *DD(J34,ICD,ISYM)
C
                                    ELSE
C
                                       WW(J2,J1,J4,J3,IQ02,IQ01)
     &                                    = WW(J2,J1,J4,J3,IQ02,IQ01)
     &                                    + DD(J12,IAB,ISYM)
     &                                    *DD(J34,ICD,ISYM)
C
                                    END IF
C
                                 END DO
C ================================================================ J3 J4
                              END DO
C ================================================================ J1 J2
C
                           END IF
                        END DO
C ================================================================= ISYM
C
                        NTERMSUM = 0
                        NSCAN = NSCAN + 1
C
C ======================================================================
C                now look for non-0 weigths and store them
C ============================================================ IQ01 IQ02
                        DO IQ02 = IQBOT_CHI,IQTOP_CHI
                           DO IQ01 = IQBOT_CHI,IQTOP_CHI
C
                              IF ( .NOT.SYM_CONNECTED(IQ1,IQ01) .AND. 
     &                             .NOT.SYM_CONNECTED(IQ2,IQ01) ) CYCLE
                              IF ( .NOT.SYM_CONNECTED(IQ1,IQ02) .AND. 
     &                             .NOT.SYM_CONNECTED(IQ2,IQ02) ) CYCLE
C
C ================================================================ J1 J2
                              DO J1 = J1L,J1H
                                 DO J2 = J2L,J2H
C ================================================================ J3 J4
                                    DO J3 = J3L,J3H
                                       DO J4 = J4L,J4H
C
                                         WTT = WW(J1,J2,J3,J4,IQ01,IQ02)
                                         WR = DREAL(WTT)
                                         WI = DIMAG(WTT)
                                         WA = ABS(WR) + ABS(WI)
C-----------------------------------------------------------------------
                                         IF ( WA.GT.TOLTT ) THEN
C
                                         IF ( KMROT.EQ.0 ) THEN
                                         IF ( ABS(WI).GT.TOLTT .AND. 
     &                                      NO_IMAG_WTT )
     &                                      WRITE (6,99008) IA,IB,IC,ID,
     &                                      J1,J2,J3,J4,WTT
                                         END IF
                                         NTERMSUM = NTERMSUM + 1
C
                                         JTT = JTT + 1
C
                                         JTT1_JTT = (IND0QX(IQ02)+J2-1)
     &                                      *NKKR + IND0QX(IQ01) + J1
                                         JTT2_JTT = (IND0QX(IQ01)+J4-1)
     &                                      *NKKR + IND0QX(IQ02) + J3
                                         JTTX_JTT = (IND0QX(IQ02)+J3-1)
     &                                      *NKKR + IND0QX(IQ01) + J4
                                         WTTJ_JTT = WTT*WNORM
C
                                         WRITE (IOTMP1) JTT1_JTT,
     &                                      JTT2_JTT,JTTX_JTT,WTTJ_JTT
C
                                         END IF
C-----------------------------------------------------------------------
C
                                       END DO
                                    END DO
C ================================================================ J3 J4
                                 END DO
                              END DO
C ================================================================ J1 J2
                           END DO
                        END DO
C ============================================================ IQ01 IQ02
C
                        NTERMSUMTOP = MAX(NTERMSUM,NTERMSUMTOP)
C
                        IF ( NTERMSUM.GT.0 ) THEN
C
                           ITT = ITT + 1
C
                           WRITE (IOTMP2) IA,IB,IC,ID,IQ1,IQ2,NTERMSUM
C
                        END IF
C
                     END DO
                  END DO
C ================================================================ IC ID
               END DO
            END DO
C ================================================================ IA IB
            NTKTKQQ(IQ1,IQ2) = ITT - NTKTKQQ(IQ1,IQ2)
         END DO
      END DO
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      NTKTKLIN = ITT
      ITTMAX = ITT
C
      JTTMAX = JTT
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         ALLOCATE (IWORK(NQ*NQ))
         CALL DRV_MPI_REDUCE_I(NTKTKQQ(1,1),IWORK(1),NQ*NQ)
         DEALLOCATE (IWORK)
C
         ALLOCATE (IWORK(0:NPROCS))
         ALLOCATE (N_ITT_PROC(0:NPROCS-1),ITT1_PROC(0:NPROCS-1))
         ALLOCATE (N_JTT_PROC(0:NPROCS-1),JTT1_PROC(0:NPROCS-1))
         N_ITT_PROC(:) = 0
         N_JTT_PROC(:) = 0
         N_ITT_PROC(MPI_ID) = ITTMAX
         N_JTT_PROC(MPI_ID) = JTTMAX
C
         CALL DRV_MPI_REDUCE_I(N_ITT_PROC(0),IWORK(0),NPROCS)
         CALL DRV_MPI_REDUCE_I(N_JTT_PROC(0),IWORK(0),NPROCS)
         DEALLOCATE (IWORK)
C
         CALL DRV_MPI_BARRIER
C
         IF ( MPI_ID.EQ.0 ) THEN
C
            ITT1_PROC(0) = 1
            JTT1_PROC(0) = 1
C
            DO IPROC = 1,(NPROCS-1)
C
               ITT1_PROC(IPROC) = ITT1_PROC(IPROC-1)
     &                            + N_ITT_PROC(IPROC-1)
               JTT1_PROC(IPROC) = JTT1_PROC(IPROC-1)
     &                            + N_JTT_PROC(IPROC-1)
               ITTMAX = ITTMAX + N_ITT_PROC(IPROC)
               JTTMAX = JTTMAX + N_JTT_PROC(IPROC)
C
            END DO
C
            NTKTKLIN = ITTMAX
C
         END IF
C
         CALL DRV_MPI_BARRIER
C
         CALL DRV_MPI_BCAST_I(0,N_ITT_PROC(0),NPROCS)
         CALL DRV_MPI_BCAST_I(0,N_JTT_PROC(0),NPROCS)
         CALL DRV_MPI_BCAST_I(0,ITT1_PROC(0),NPROCS)
         CALL DRV_MPI_BCAST_I(0,JTT1_PROC(0),NPROCS)
         CALL DRV_MPI_BCAST_I(0,ITTMAX,1)
         CALL DRV_MPI_BCAST_I(0,JTTMAX,1)
         CALL DRV_MPI_BCAST_I(0,NTKTKLIN,1)
C
C
C ----------------------------------------------------------------------
C             get the variables JTT* and ITT*  from file
C ----------------------------------------------------------------------
C
         ALLOCATE (WTTJ(JTTMAX),JTT1(JTTMAX))
         ALLOCATE (JTT2(JTTMAX),JTTX(JTTMAX))
C
         WTTJ(1:JTTMAX) = C0
         JTT1(1:JTTMAX) = 0
         JTT2(1:JTTMAX) = 0
         JTTX(1:JTTMAX) = 0
C
         REWIND IOTMP1
         DO J = JTT1_PROC(MPI_ID),JTT1_PROC(MPI_ID) + N_JTT_PROC(MPI_ID)
     &      - 1
            READ (IOTMP1) JTT1(J),JTT2(J),JTTX(J),WTTJ(J)
         END DO
C
         ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX),ITTD(ITTMAX))
         ALLOCATE (ITTQ1(ITTMAX),ITTQ2(ITTMAX),NTTJ(ITTMAX),STAT=IA_ERR)
C
         NTTJ(1:ITTMAX) = 0
         ITTA(1:ITTMAX) = 0
         ITTB(1:ITTMAX) = 0
         ITTC(1:ITTMAX) = 0
         ITTD(1:ITTMAX) = 0
         ITTQ1(1:ITTMAX) = 0
         ITTQ2(1:ITTMAX) = 0
C
         REWIND IOTMP2
         DO I = ITT1_PROC(MPI_ID),ITT1_PROC(MPI_ID) + N_ITT_PROC(MPI_ID)
     &      - 1
            READ (IOTMP2) ITTA(I),ITTB(I),ITTC(I),ITTD(I),ITTQ1(I),
     &                    ITTQ2(I),NTTJ(I)
         END DO
C
         CALL DRV_MPI_BARRIER
C
         DO ID_SEND = 0,(NPROCS-1)
C
            J1 = JTT1_PROC(ID_SEND)
            NJ = N_JTT_PROC(ID_SEND)
            CALL DRV_MPI_BCAST_I(ID_SEND,JTT1(J1),NJ)
            CALL DRV_MPI_BCAST_I(ID_SEND,JTT2(J1),NJ)
            CALL DRV_MPI_BCAST_I(ID_SEND,JTTX(J1),NJ)
            CALL DRV_MPI_BCAST_C(ID_SEND,WTTJ(J1),NJ)
C
            I1 = ITT1_PROC(ID_SEND)
            NI = N_ITT_PROC(ID_SEND)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTA(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTB(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTC(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTD(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTQ1(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,ITTQ2(I1),NI)
            CALL DRV_MPI_BCAST_I(ID_SEND,NTTJ(I1),NI)
C
         END DO
C
         CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ELSE
C
C ----------------------------------------------------------------------
C             get the variables JTT* and ITT*  from file
C ----------------------------------------------------------------------
C
         ALLOCATE (WTTJ(JTTMAX),JTT1(JTTMAX))
         ALLOCATE (JTT2(JTTMAX),JTTX(JTTMAX))
C
         REWIND IOTMP1
         DO J = 1,JTTMAX
            READ (IOTMP1) JTT1(J),JTT2(J),JTTX(J),WTTJ(J)
         END DO
C
         ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX),ITTD(ITTMAX))
         ALLOCATE (ITTQ1(ITTMAX),ITTQ2(ITTMAX),NTTJ(ITTMAX),STAT=IA_ERR)
C
         REWIND IOTMP2
         DO I = 1,NTKTKLIN
            READ (IOTMP2) ITTA(I),ITTB(I),ITTC(I),ITTD(I),ITTQ1(I),
     &                    ITTQ2(I),NTTJ(I)
         END DO
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      CLOSE (IOTMP1)
      CLOSE (IOTMP2)
C
      IF ( MPI_ID.EQ.0 ) THEN
C
C ----------------------------------------------------------------------
C                      print out information
C ----------------------------------------------------------------------
C
         WRITE (6,99005) (NQ**2)*(NXM**4),NSCAN,NTKTKLIN,ITTMAX
         DO IQ = IQBOT_CHI,IQTOP_CHI
            WRITE (6,99006) IQ,IQTOP_CHI,
     &                      (NTKTKQQ(IQ,JQ),JQ=IQBOT_CHI,IQTOP_CHI)
            DO JQ = IQBOT_CHI,IQTOP_CHI
               NTKTKMAX = MAX(NTKTKMAX,NTKTKQQ(IQ,JQ))
            END DO
         END DO
         WRITE (6,99007) NTKTKMAX,JTTMAX,NTERMSUMTOP
C
C ----------------------------------------------------------------------
C   now write the COMPLETE tables to file  TKTK_tables  for later use
C ----------------------------------------------------------------------
C
         REWIND IOTMP
C
         WRITE (IOTMP) ERYDA_EQ_ERYDB,IQBOT_CHI,IQTOP_CHI,NL
         WRITE (IOTMP) NXM,NSYMACCEPTED
         WRITE (IOTMP) NKKR,NSYM
C
         DO ISYM = 1,NSYM
            WRITE (IOTMP) SYMACCEPTED(ISYM)
            WRITE (IOTMP) SYMUNITARY(ISYM)
         END DO
C
         WRITE (IOTMP) ITTMAX,JTTMAX,NSCAN,NTERMSUMTOP,
     &                 ((NTKTKQQ(IQ,JQ),JQ=IQBOT_CHI,IQTOP_CHI),
     &                 IQ=IQBOT_CHI,IQTOP_CHI)
C
         DO I = 1,NTKTKLIN
            WRITE (IOTMP) ITTA(I),ITTB(I),ITTC(I),ITTD(I),ITTQ1(I),
     &                    ITTQ2(I),NTTJ(I)
         END DO
C
         DO J = 1,JTTMAX
            WRITE (IOTMP) JTT1(J),JTT2(J),JTTX(J),WTTJ(J)
         END DO
C
      END IF
C
      CLOSE (IOTMP)
C
C=======================================================================
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      DEALLOCATE (DD,WW,JADD,JBDD,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
      RETURN
99001 FORMAT (5X,'symmetry-derived coefficients for BZ-integral',
     &        ' available',/,5X,'data read from file ',A,/)
99002 FORMAT (5X,'symmetry-derived coefficients for BZ-integral',
     &        ' NOT available',/,5X,
     &        'data will be created and stored on file ',A,/)
99003 FORMAT (//,1X,79('*'),/,29X,'<CHITKTKSYM>',/,1X,79('*'),//,5X,
     &        'BZ-integral   Int d^3k  TAU(k)[L1,L2] * TAU(k)[L3,L4]',
     &        //,5X,'for   (E_a=E_b) = ',L1,/)
99004 FORMAT (5X,'index-range ',A,/,(5X,24I3))
99005 FORMAT (/,5X,'max. number of TAU(k)*TAU(k) elements     =',I12,/,
     &        5X,'number of elements scanned       NSCAN    =',I12,/,5X,
     &        'total number of non-0 elements   NTKTKLIN =',I12,/,5X,
     &        'array size                       ITTMAX   =',I12,/,5X,
     &        'number of non-0 elements         NTKTKQQ')
99006 FORMAT (5X,'for IQ=',I2,'  JQ=1,..,',I2,5X,30I7)
99007 FORMAT (5X,'size for site dependent array    NTKTKMAX =',I12,//,
     &        5X,'number of terms                  JTTMAX   =',I12,/,5X,
     &        'maximum number of terms to sum   NTERMSTOP=',I12,//)
99008 FORMAT ('<CHITKTKSYM> Im(W) > 0',8I3,2E12.5)
      END
