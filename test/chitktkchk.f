C*==chitktkchk.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHITKTKCHK(ALLTKTK,ERYD,P,IPRINT,ICPAFLAG,CPACHNG,
     &                      ITCPA,ICPACONV,TSST,MSST,TSSQ,MSSQ,TAUQ,
     &                      TAUQZ)
C   ********************************************************************
C   *                                                                  *
C   *  the BZ-integration for the expression  TAU(k)*TAU(k) is done    *
C   *  twice for one representative arbitrary k-vector  KVEC           *
C   *                                                                  *
C   *  run A: KVEC is rotated according to the symmetry group          *
C   *         of the crystal system (i.e. 48 times for cubic lattices) *
C   *         and MTTA is summed over                                  *
C   *                                                                  *
C   *  run B: KVEC is rotated according to the rotations that create   *
C   *         the irreducible wedge. only the non-0 elements of MTTB   *
C   *         are summed using the index and weight table produced     *
C   *         by <CHITKTKSYM>                                          *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  in part 2 the subroutine   <CHIKLOOPS>  is checked              *
C   *                                                                  *
C   *  run B: set up the symmetry reduced k-mesh via <KMESHS>          *
C   *         run <CHIKLOOPS> making use of the symmetry weights WTTJ  *
C   *  run A: set up the full k-mesh via <KMESHS> by setting NSYM=1    *
C   *         run <CHIKLOOPS> setting the symmetry weights WTTJ=1      *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array for run B !        *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   * 04/07/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NZ12,NZ12MAX,ERYDA_EQ_ERYDB,CHIZ,NTKTKLIN,
     &    NTTJ,JTT1,JTT2,JTTX,WTTJ,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,
     &    ITTMAX
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,DROTQ
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM
      USE MOD_SYMMETRY,ONLY:NWEDGE,NSYMACCEPTED,IQORGQP,DROT,SYMUNITARY,
     &    SYMDET,SYMACCEPTED,MROTK,MROTR,IWEDGEROT,SYMCRYSYS,NSYM
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CHITKTKCHK49
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHITKTKCHK')
      REAL*8 TOLA,TOLAB
      PARAMETER (TOLA=1D-5,TOLAB=1D-7)
C
C Dummy arguments
C
      LOGICAL ALLTKTK
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUM1,CSUMX,MAUX(:,:),
     &           MTTA(NKMMAX*NKMMAX*NQMAX,NKMMAX*NKMMAX*NQMAX,2),
     &           MTTB(NKMMAX*NKMMAX*NQMAX,NKMMAX*NKMMAX*NQMAX,2),
     &           MTTBLIN(ITTMAX,2),TAUK(:,:),TAUKLIN(:),
     &           TAUQA(NKMMAX,NKMMAX,NQMAX),TTA,TTB,WT
      INTEGER I,I0,IA,IA_ERR,IB,IC,ID,IK,INFO,IPIV(NKKR),IQ,IROT,ISYM,
     &        IWEDGE,IX,IZ,J,J0,JQ,JTT,KADI,KBCJ,KCBJ,KDAI,NERRAB(2),
     &        NERRTAU,NKA,NKMSQ,NKPTS0,NKTAB,NKTABMAX,NON0A,NSYMDUM
      REAL*8 KTAB(:,:),KVEC(3),KVECP(3),MVEC(3),MVECP(3),RVEC(3),
     &       RVECP(3),SCLA,SCLB,WKTAB(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUK,TAUKLIN,KTAB,WKTAB
C
      ALLOCATE (TAUKLIN(NKKR*NKKR))
      ALLOCATE (MAUX(NKKR,NKKR),TAUK(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUKLIN')
C
      WRITE (6,99009)
C
      NKMSQ = NKM*NKM
C=======================================================================
      RVEC(1) = 0.412D0
      RVEC(2) = 0.231D0
      RVEC(3) = 0.123D0
      KVEC(1) = 0.412D0
      KVEC(2) = 0.231D0
      KVEC(3) = 0.123D0
      MVEC(1) = 0.0D0
      MVEC(2) = 0.0D0
      MVEC(3) = 1.0D0
C=======================================================================
      WRITE (6,99005)
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,RVEC,1,0D0,RVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,MVEC,1,0D0,MVECP,1)
            CALL DSCAL(3,DBLE(SYMDET(ISYM)),MVECP,1)
C
            WRITE (6,99004) ISYM,SYMUNITARY(ISYM),'M',MVECP,'R',RVECP,
     &                      'k',KVECP
C
         END IF
C
      END DO
C=======================================================================
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
C=======================================================================
C                straightforward integration   run A
C=======================================================================
C
      WRITE (6,99007)
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      CALL CINIT(NKMMAX**4*NQMAX**2*2,MTTA)
C
      IK = 0
      DO ISYM = 1,NSYM
         IF ( SYMCRYSYS(ISYM) ) THEN
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
C
            IK = IK + 1
            WRITE (6,99006) 'A',IK,KVECP
C
            CALL STRSET(IK,KVECP,TAUK,MAUX,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,NKMMAX)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
            DO JQ = 1,NQ
               J0 = IND0Q(JQ)
               DO IQ = 1,NQ
                  I0 = IND0Q(IQ)
                  DO ID = 1,NKMQ(IQ)
                     DO IC = 1,NKMQ(JQ)
                        DO IB = 1,NKMQ(JQ)
                           DO IA = 1,NKMQ(IQ)
                              KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
                              KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
                              MTTA(KADI,KBCJ,1) = MTTA(KADI,KBCJ,1)
     &                           + TAUK(I0+IA,J0+IB)*TAUK(J0+IC,I0+ID)
                              MTTA(KADI,KBCJ,2) = MTTA(KADI,KBCJ,2)
     &                           + TAUK(I0+IA,J0+IB)
     &                           *DCONJG(TAUK(I0+ID,J0+IC))
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
         END IF
C
      END DO
      NKA = IK
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      SCLA = 1D0/DBLE(NKA)
      NON0A = 0
      DO KADI = 1,NKMSQ*NQ
         DO KBCJ = 1,NKMSQ*NQ
            MTTA(KADI,KBCJ,1) = MTTA(KADI,KBCJ,1)*SCLA
            MTTA(KADI,KBCJ,2) = MTTA(KADI,KBCJ,2)*SCLA
            IF ( ABS(MTTA(KADI,KBCJ,1)).GT.TOLA ) NON0A = NON0A + 1
         END DO
      END DO
      WRITE (6,99001) NON0A,TOLA,(NKMSQ*NQ)**2
C
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      CALL CINIT(NKMMAX**4*NQMAX**2*2,MTTB)
      CALL CINIT(ITTMAX*2,MTTBLIN)
C
      IK = 0
      DO IWEDGE = 1,NWEDGE
C
         IROT = IWEDGEROT(IWEDGE)
C
         CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KVEC,1,0D0,KVECP,1)
C
         IK = IK + 1
         WRITE (6,99006) 'B',IK,KVECP
C
         CALL STRSET(IK,KVECP,TAUKLIN,MAUX,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLIN,MSSQ,NQMAX,NKKR,NKMMAX)
C
         CALL ZGETRF(NKKR,NKKR,TAUKLIN,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUKLIN,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C--------------------------------------------------------- store TAU*TAU
         JTT = 0
         DO I = 1,NTKTKLIN
            CSUM1 = C0
            CSUMX = C0
            DO J = 1,NTTJ(I)
               JTT = JTT + 1
               WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
               CSUM1 = CSUM1 + WT*TAUKLIN(JTT2(JTT))
               CSUMX = CSUMX + WT*DCONJG(TAUKLIN(JTTX(JTT)))
            END DO
            MTTBLIN(I,1) = MTTBLIN(I,1) + CSUM1
            MTTBLIN(I,2) = MTTBLIN(I,2) + CSUMX
         END DO
C-----------------------------------------------------------------------
      END DO
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      SCLB = 1D0/DBLE(NWEDGE)
C
      DO I = 1,NTKTKLIN
C
         IA = ITTA(I)
         IB = ITTB(I)
         IC = ITTC(I)
         ID = ITTD(I)
         IQ = ITTQ1(I)
         JQ = ITTQ2(I)
         KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
         KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
         DO IZ = 1,2
            MTTB(KADI,KBCJ,IZ) = MTTBLIN(I,IZ)*SCLB
            IF ( ERYDA_EQ_ERYDB ) THEN
               KCBJ = (IC-1)*NKM + IB + (JQ-1)*NKMSQ
               KDAI = (IA-1)*NKM + IA + (IQ-1)*NKMSQ
               MTTB(KCBJ,KDAI,IZ) = MTTB(KADI,KBCJ,IZ)
            END IF
         END DO
C
      END DO
C
C=======================================================================
C
      NERRAB(1) = 0
      NERRAB(2) = 0
      IF ( ALLTKTK ) THEN
C=======================================================================
         WRITE (6,99010)
         DO JQ = 1,NQ
            J0 = IND0Q(JQ)
            DO IQ = 1,NQ
               I0 = IND0Q(IQ)
               DO ID = 1,NKMQ(IQ)
                  DO IC = 1,NKMQ(JQ)
                     DO IB = 1,NKMQ(JQ)
                        DO IA = 1,NKMQ(IQ)
                           DO IZ = 1,2
                              KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
                              KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
                              TTA = MTTA(KADI,KBCJ,IZ)
                              TTB = MTTB(KADI,KBCJ,IZ)
                              IF ( ABS(TTA)+ABS(TTB).GT.TOLAB ) THEN
C                                 WRITE (6,99003) IA,IB,IC,ID,IQ,JQ,TTA,
C     &                                  TTB,IZ
                                 IF ( ABS(TTA-TTB).GT.TOLAB ) NERRAB(IZ)
     &                                = NERRAB(IZ) + 1
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C=======================================================================
      ELSE
C=======================================================================
         WRITE (6,99011)
         DO I = 1,NTKTKLIN
C
            IA = ITTA(I)
            IB = ITTB(I)
            IC = ITTC(I)
            ID = ITTD(I)
            IQ = ITTQ1(I)
            JQ = ITTQ2(I)
C
            DO IZ = 1,2
               KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
               KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
               TTA = MTTA(KADI,KBCJ,IZ)
               TTB = MTTB(KADI,KBCJ,IZ)
               IF ( ABS(TTA)+ABS(TTB).GT.TOLAB ) THEN
                  IF ( ABS(TTA-TTB).GT.TOLAB ) THEN
                     NERRAB(IZ) = NERRAB(IZ) + 1
                     WRITE (6,99002) IA,IB,IC,ID,IQ,JQ,TTA,TTB,IZ
                  END IF
               END IF
            END DO
         END DO
C=======================================================================
      END IF
C
      WRITE (6,99008) 1,NON0A,NTKTKLIN,TOLAB,NERRAB
C
C
C=======================================================================
C                               PART  2
C
C             check the k-integration routine  <CHIKLOOPS>
C
C=======================================================================
C
      NZ12 = 2
      NCPA = 0
C
      NERRAB(1) = 0
      NERRAB(2) = 0
C
      NKPTS0 = 100
C
C ------------------------- run B:    set up the symmetry reduced k-mesh
C ------------------- run <CHIKLOOPS> making use of the symmetry weights
C
      CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB')
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
      CALL CHIKLOOPS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,NQ,NKM,NKKR,NKMQ,
     &               IND0Q,KMROT,DROTQ,CPATOL,NCPA,ICPA,ICPAALG,ITCPA,
     &               ITCPAMAX,ICPACONV,DROT,IQORGQP,SYMUNITARY,
     &               SYMACCEPTED,NOQ,ITOQ,CONC,WKTAB,KTAB,NKTAB,NSYM,
     &               NSYMACCEPTED,NTMAX,NQMAX,NKMMAX,NKTABMAX,TSST,MSST,
     &               TSSQ,MSSQ,TAUQ,TAUQZ)
C
      MTTB(:,:,:) = CHIZ(:,:,:)
C
C ------------------------------------- run A:    set up the full k-mesh
C -------------------- run <CHIKLOOPS> setting the symmetry weights to 1
C
      NSYMDUM = 1
C
      CALL KMESHS(IOTMP,NSYMDUM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      DEALLOCATE (KTAB,WKTAB)
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB')
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
      DO I = 1,NTKTKLIN
C
         IA = ITTA(I)
         IB = ITTB(I)
         IC = ITTC(I)
         ID = ITTD(I)
         IQ = ITTQ1(I)
         JQ = ITTQ2(I)
C
         NTTJ(I) = 1
         JTT1(I) = (IND0Q(JQ)+IB-1)*NKKR + IND0Q(IQ) + IA
         JTT2(I) = (IND0Q(IQ)+ID-1)*NKKR + IND0Q(JQ) + IC
         JTTX(I) = (IND0Q(JQ)+IC-1)*NKKR + IND0Q(IQ) + ID
         WTTJ(I) = 1D0
C
      END DO
C
      CALL CHIKLOOPS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,NQ,NKM,NKKR,NKMQ,
     &               IND0Q,KMROT,DROTQ,CPATOL,NCPA,ICPA,ICPAALG,ITCPA,
     &               ITCPAMAX,ICPACONV,DROT,IQORGQP,SYMUNITARY,
     &               SYMACCEPTED,NOQ,ITOQ,CONC,WKTAB,KTAB,NKTAB,NSYM,
     &               NSYMACCEPTED,NTMAX,NQMAX,NKMMAX,NKTABMAX,TSST,MSST,
     &               TSSQ,MSSQ,TAUQA,TAUQZ)
C
      MTTA(:,:,:) = CHIZ(:,:,:)
C
C ---------------------------------- compare the calculated TAU-matrices
C
      NERRTAU = 0
      DO IQ = 1,NQ
         DO IB = 1,NKMQ(IQ)
            DO IA = 1,NKMQ(IQ)
               TTA = TAUQA(IA,IB,IQ)
               TTB = TAUQ(IA,IB,IQ)
               IF ( ABS(TTA)+ABS(TTB).GT.TOLAB ) THEN
                  IF ( ABS(TTA-TTB).GT.TOLAB ) THEN
                     NERRTAU = NERRTAU + 1
                     WRITE (6,99003) IA,IB,IQ,TTA,TTB
                  END IF
               END IF
            END DO
         END DO
      END DO
C
      IF ( ALLTKTK ) THEN
C=======================================================================
         WRITE (6,99010)
         DO JQ = 1,NQ
            J0 = IND0Q(JQ)
            DO IQ = 1,NQ
               I0 = IND0Q(IQ)
               DO ID = 1,NKMQ(IQ)
                  DO IC = 1,NKMQ(JQ)
                     DO IB = 1,NKMQ(JQ)
                        DO IA = 1,NKMQ(IQ)
                           DO IZ = 1,2
                              KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
                              KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
                              TTA = MTTA(KADI,KBCJ,IZ)
                              TTB = MTTB(KADI,KBCJ,IZ)
                              IF ( ABS(TTA)+ABS(TTB).GT.TOLAB ) THEN
C                                 WRITE (6,99003) IA,IB,IC,ID,IQ,JQ,TTA,
C     &                                  TTB
                                 IF ( ABS(TTA-TTB).GT.TOLAB ) NERRAB(IZ)
     &                                = NERRAB(IZ) + 1
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C=======================================================================
      ELSE
C=======================================================================
         WRITE (6,99011)
         DO I = 1,NTKTKLIN
C
            IA = ITTA(I)
            IB = ITTB(I)
            IC = ITTC(I)
            ID = ITTD(I)
            IQ = ITTQ1(I)
            JQ = ITTQ2(I)
C
            DO IZ = 1,2
               KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
               KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
               TTA = MTTA(KADI,KBCJ,IZ)
               TTB = MTTB(KADI,KBCJ,IZ)
               IF ( ABS(TTA)+ABS(TTB).GT.TOLAB ) THEN
                  IF ( ABS(TTA-TTB).GT.TOLAB ) THEN
                     NERRAB(IZ) = NERRAB(IZ) + 1
                     WRITE (6,99002) IA,IB,IC,ID,IQ,JQ,TTA,TTB,IZ
                  END IF
               END IF
            END DO
         END DO
C=======================================================================
      END IF
C
      WRITE (6,99008) 2,NON0A,NTKTKLIN,TOLAB,NERRAB,NERRTAU
C
C=======================================================================
      DEALLOCATE (MAUX,TAUK,TAUKLIN,KTAB,WKTAB)
C
      STOP 'TEST 3 completed'
99001 FORMAT (/,5X,'run A: ',I8,'  elements >',F12.8,'  out of ',I8,/)
99002 FORMAT (4I3,2x,2I3,2x,2E17.8,'  A',/,22x,2E17.8,'  B  ERROR',I4,/)
99003 FORMAT (2I3,2x,1I3,2x,2E17.8,' T A',/,13x,2E17.8,' T B  ERROR',/)
99004 FORMAT (' ISYM ',I3,' unitary ',L1,3X,A,1X,3F5.1,:,2X,A,1X,3F5.1,
     &        :,2X,A,1X,3F5.1)
99005 FORMAT (5X,'action of the symmetry operations on vector ->R',
     &        ' in real space',/,5X,'vector ->k in reciprocal space',
     &        ' and magnetisation ->M',/,5X,
     &        'antiunitary rotations of ->k include inversion I',/)
99006 FORMAT (5x,'  k-loop    run   ',A,I4,' ->k ',3F8.4)
99007 FORMAT (//,5x,'running k-loops:   A (straightforward)  and ',
     &        ' B (symmetry based)',/)
99008 FORMAT (/,5X,'<CHITKTKCHK>  part ',I2,' finished',/,5X,
     &        'run A:   number of non-0 elements found    ',I8,/,5X,
     &        'run B:   number of non-0 elements expected ',I8,/,5X,
     &        'for tolerance                          ',F12.8,//,5X,
     &        'number of deviations run A and B           ',2I8,//,:,5X,
     &        'number of deviations run A and B for TAU   ',2I8,//)
99009 FORMAT (//,1X,79('*'),/,29X,'<CHITKTKCHK>',/,1X,79('*'),//,5X,
     &        'check BZ-integral   ',
     &        'Int d^3k  TAU(k)[L1,L2] * TAU(k)[L3,L4]',/)
99010 FORMAT (/,5X,'checking ALL NON-VANISHING matrix elements')
99011 FORMAT (/,5X,'checking SELECTED matrix elements')
      END
