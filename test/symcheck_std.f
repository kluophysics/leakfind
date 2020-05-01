C*==symcheck_std.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMCHECK_STD(ERYD,P,ICPAFLAG,CPACHNG,IPRINT,ITCPA,
     &                        ICPACONV,TSST,MSST,TSSQ,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *   check the symmetry properties exploited for BZ integration     *
C   *                                                                  *
C   *           called for   ITEST = 2                                 *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  PART 1:                                                         *
C   *                                                                  *
C   *  check   for each symmetry operation R:                          *
C   *          TAU[q'q'](Rk) = R TAU[qq](k) R+                         *
C   *            G[q'q'](Rk) = R   G[qq](k) R+       for q' = R q      *
C   *            m[q'q']     = R   m[qq]    R+                         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  PART 2:                                                         *
C   *                                                                  *
C   *  the BZ-integration for the site-diagonal TAUQ is done twice     *
C   *  for one representative arbitrary k-vector  KVEC                 *
C   *                                                                  *
C   *  run A: KVEC is rotated according to the symmetry group          *
C   *         of the crystal system (i.e. 48 times for cubic lattices) *
C   *         and TAUQA is summed over                                 *
C   *                                                                  *
C   *  run B: KVEC is rotated according to the rotations that create   *
C   *         the irreducible wedge. only the non-0 elements of TAUQB  *
C   *         are summed using the index and weight table produced     *
C   *         by <SYMLATTICE>                                          *
C   *                                                                  *
C   *  for each IQ the matrix TAUQA is printed as well as the          *
C   *  difference   TAUQB - TAUQA, i.e. the 2nd matrix should be blank *
C   *------------------------------------------------------------------*
C   *  PART 3:                                                         *
C   *                                                                  *
C   *  run C: KVEC is rotated according to the rotations that create   *
C   *         the irreducible wedge. the remaining wedges are          *
C   *         accounted for by multiplication with the rotation        *
C   *         matriced DROT                                            *
C   *                                                                  *
C   *  for each IQ the printed matrix  TAUQC - TAUQA  should be blank  *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  PART 4:                                                         *
C   *                                                                  *
C   *  <KMESHS>  and  <TAU_STD_KPOINTS>                                *
C   *  are run WITH and WITHOUT  making use of symmetry                *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  PART 5:                                                         *
C   *                                                                  *
C   *  check the relation   TAU(q',q') = R TAU(q,q) R+   for q' = R q  *
C   *                                                                  *
C   * 04/07/2000 HE - update 02/07/2011                                *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_TYPES,ONLY:NTMAX,CONC
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_KSPACE,ONLY:NELMT,WTAUUV,VTAUUV,UTAUUV,NTAUUV,QTBZ,JTBZ,
     &    ITBZ,NKTAB,KTAB,WKTAB
      USE MOD_SITES,ONLY:NQ,NQMAX,NOQ,ITOQ,ICPA
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM
      USE MOD_SYMMETRY,ONLY:SYMUNITARY,SYMDET,IQORGQP,DROT,NSYMACCEPTED,
     &    SYMACCEPTED,MROTK,MROTR,IWEDGEROT,SYMCRYSYS,NSYM,NWEDGE,
     &    SYMTAUTET_AVAILABLE,IQREPQ,IQPSYMQ
      USE MOD_CALCMODE,ONLY:KMROT,IREL,ITEST
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--SYMCHECK_STD74
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMCHECK_STD')
      REAL*8 TOLR,TOLM,TOL3
      PARAMETER (TOLR=1.0D-8,TOLM=1.0D-8,TOL3=1.0D-8)
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      COMPLEX*16 CSUM,GQ(:,:,:),GQP(:,:,:),MAUX(:,:),MQ(:,:,:),
     &           MQP(:,:,:),MSSQ0(:,:,:),MSSQ_4A(:,:,:),MSSQ_4B(:,:,:),
     &           PHASK(1),SCLC,TA,TAUK(:,:),TAUQA(:,:,:),TAUQB(:,:,:),
     &           TAUQC(:,:,:),TC,TSSQ0(:,:,:),W1(:,:),W2(:,:),WG(:,:),
     &           WM(:,:)
      REAL*8 ERRMAX,KVEC(3),KVECP(3),MVEC(3),MVECP(3),RVEC(3),RVECP(3),
     &       SCLA,SCLB
      INTEGER I,I1,I3,IA_ERR,IE,IK,INFO,IO,IPIV(:),IQ,IQP,IROT,ISYM,IT,
     &        IWEDGE,IX,J,J1,M,NELMTMAX,NERR2,NERR3,NERR4,NERR5,NERRG,
     &        NERRGS,NERRGSS,NERRM,NERRMS,NERRMSS,NERRT,NERRTS,NERRTSS,
     &        NKA,NKPTS0,NKTABMAX,NSYMACCEPTED_SAV,NSYMDUM,NSYMLOC,
     &        NSYM_SAV,SNTAUUVMAX
      CHARACTER*5 STRERR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUK,IPIV,TAUQA,TAUQB,TAUQC,W1,W2,GQ,MQ,WG
      ALLOCATABLE WM,GQP,MQP,MSSQ0,TSSQ0,MSSQ_4A,MSSQ_4B
C
C=======================================================================
C           supply tables for tetrahedron integration in PART 2
C=======================================================================
C
      DEALLOCATE (ITBZ,JTBZ,QTBZ,NTAUUV,UTAUUV,VTAUUV,WTAUUV)
C
      CALL SYMTAU(ITEST)
C
C-----------------------------------------------------------------------
      IF ( SYMTAUTET_AVAILABLE ) THEN
C
         REWIND (IOTMP)
         READ (IOTMP) NELMTMAX,SNTAUUVMAX
         WRITE (6,99023) ' TAU-elements     NELMT: ',NELMT,NELMTMAX
C
         NELMT = NELMTMAX
C
         ALLOCATE (ITBZ(NELMTMAX),JTBZ(NELMTMAX))
         ALLOCATE (QTBZ(NELMTMAX),NTAUUV(NELMTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: QTBZ')
C
         ALLOCATE (UTAUUV(SNTAUUVMAX))
         ALLOCATE (VTAUUV(SNTAUUVMAX),WTAUUV(SNTAUUVMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WTAUUV')
C
         READ (IOTMP) (ITBZ(I),JTBZ(I),NTAUUV(I),QTBZ(I),I=1,NELMTMAX)
         READ (IOTMP) (UTAUUV(I),VTAUUV(I),WTAUUV(I),I=1,SNTAUUVMAX)
         CLOSE (IOTMP)
C-----------------------------------------------------------------------
      ELSE
         CLOSE (IOTMP)
      END IF
C
C=======================================================================
C
      WRITE (6,99021)
C
C=======================================================================
C            in case of Full Potential: symmetrize m-matrix
C=======================================================================
C
      IF ( FULLPOT .AND. NQ.LT.0 ) THEN
C
         ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX))
         DO IQ = 1,NQ
            IF ( IQ.EQ.IQREPQ(IQ) ) THEN
               DO IO = 1,NOQ(IQ)
C
                  IT = ITOQ(IO,IQ)
C
                  NSYMLOC = 0
                  W1(:,:) = C0
C
                  DO ISYM = 1,NSYM
                     IF ( SYMACCEPTED(ISYM) ) THEN
                        IF ( IQ.EQ.IQPSYMQ(ISYM,IQ) ) THEN
                           NSYMLOC = NSYMLOC + 1
C
                           CALL ROTATEGEN(TSST(1,1,IT),W1,NKM,
     &                        DROT(1,1,ISYM),SYMUNITARY(ISYM),NKMMAX)
C
                        END IF
                     END IF
                  END DO
C
                  W1(1:NKM,1:NKM) = W1(1:NKM,1:NKM)/DBLE(NSYMLOC)
C
                  W2(1:NKM,1:NKM) = TSST(1:NKM,1:NKM,IT)
     &                              - W1(1:NKM,1:NKM)
C
                  TSST(1:NKM,1:NKM,IT) = W1(1:NKM,1:NKM)
C
                  ERRMAX = 0D0
                  DO I = 1,NKM
                     DO J = 1,NKM
                        ERRMAX = MAX(ERRMAX,ABS(W2(I,J)))
                     END DO
                  END DO
C
                  WRITE (6,99024) IQ,NSYMLOC,ERRMAX
C
                  CALL CMATSTRUCT('t - <t>          ',W2,NKM,NKMMAX,
     &                            IREL,IREL,0,1D-12,6)
                  CALL CMATSTRUCT('    <t>          ',W1,NKM,NKMMAX,
     &                            IREL,IREL,0,1D-12,6)
C
               END DO
            END IF
         END DO
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
         DEALLOCATE (W1,W2)
      END IF
C=======================================================================
C
      M = NKMMAX
C
      ALLOCATE (MQ(M,M,NQMAX),WG(M,M),WM(M,M),GQP(M,M,NQMAX))
      ALLOCATE (TAUQC(M,M,NQMAX),W2(M,M),GQ(M,M,NQMAX))
      ALLOCATE (TAUQA(M,M,NQMAX),TAUQB(M,M,NQMAX),W1(M,M))
      ALLOCATE (MAUX(NKKR,NKKR),MQP(M,M,NQMAX))
      ALLOCATE (TAUK(NKKR,NKKR),IPIV(NKKR))
      ALLOCATE (TSSQ0(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSSQ_4A(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSSQ_4B(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSSQ0(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      TSSQ0(:,:,:) = TSSQ(:,:,:)
      MSSQ0(:,:,:) = MSSQ(:,:,:)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      NSYM_SAV = NSYM
C
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
C                               PART  0
C=======================================================================
C
      WRITE (6,99001) 0,' check results of symmetry operations R'
      WRITE (6,99010)
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,RVEC,1,0D0,RVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,MVEC,1,0D0,MVECP,1)
            CALL DSCAL(3,DBLE(SYMDET(ISYM)),MVECP,1)
C
            WRITE (6,99003) ISYM,SYMUNITARY(ISYM),'M',MVECP,'R',RVECP,
     &                      'k',KVECP
C
         END IF
C
      END DO
C
C=======================================================================
C                               PART  1
C
C      check   for each symmetry operation R:
C              TAU[q'q'](Rk) = R TAU[qq](k) R+
C                G[q'q'](Rk) = R   G[qq](k) R+       for q' = R q
C                m[q'q']     = R   m[qq]    R+
C=======================================================================
C
      WRITE (6,99001) 1,' check results of symmetry operations R'
      WRITE (6,99011) KVEC
C
C=======================================================================
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      CALL CINIT(M*M*NQ,TAUQ)
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      NERRTSS = 0
      NERRMSS = 0
      NERRGSS = 0
      WRITE (6,*) ' '
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
C-NOTE:  MROTK includes the INVERSION for ANTI - unitary transformations
C
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
C
C-----------------------------------------------------------------------
C                  calculate   TAU[q'q'](Rk)
C-----------------------------------------------------------------------
C
            CALL STRSET(1,KVECP,MAUX,TAUK,P)
C
            DO IQP = 1,NQ
               I1 = IND0Q(IQP) + 1
               DO J = 1,NKMQ(IQP)
                  J1 = IND0Q(IQP) + J
                  CALL ZCOPY(NKMQ(IQP),TAUK(I1,J1),1,GQP(1,J,IQP),1)
                  CALL ZCOPY(NKMQ(IQP),MSSQ(1,J,IQP),1,MQP(1,J,IQP),1)
               END DO
            END DO
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
            DO IQP = 1,NQ
               I1 = IND0Q(IQP) + 1
               DO J = 1,NKMQ(IQP)
                  J1 = IND0Q(IQP) + J
                  CALL ZCOPY(NKMQ(IQP),TAUK(I1,J1),1,TAUQ(1,J,IQP),1)
               END DO
            END DO
C-----------------------------------------------------------------------
C               get  TAU[qq](k)  for all sites  q
C-----------------------------------------------------------------------
C
            CALL STRSET(1,KVEC,MAUX,TAUK,P)
C
            DO IQP = 1,NQ
               I1 = IND0Q(IQP) + 1
               DO J = 1,NKMQ(IQP)
                  J1 = IND0Q(IQP) + J
                  CALL ZCOPY(NKMQ(IQP),TAUK(I1,J1),1,GQ(1,J,IQP),1)
                  CALL ZCOPY(NKMQ(IQP),MSSQ(1,J,IQP),1,MQ(1,J,IQP),1)
               END DO
            END DO
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C-----------------------------------------------------------------------
C            perform rotations   R TAU[qq](k) R+
C            for site q for which   q' = R q
C-----------------------------------------------------------------------
C
            NERRTS = 0
            NERRMS = 0
            NERRGS = 0
            DO IQP = 1,NQ
C
               IQ = IQORGQP(ISYM,IQP)
C
               I1 = IND0Q(IQ) + 1
               DO J = 1,NKMQ(IQ)
                  J1 = IND0Q(IQ) + J
                  CALL ZCOPY(NKMQ(IQ),TAUK(I1,J1),1,W1(1,J),1)
               END DO
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C------------------------------------------------------------------- TAU
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,W1,M,
     &                    C0,W2,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),M,
     &                    C0,W1,M)
C................................................................... Mss
C
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                    MQ(1,1,IQ),M,C0,W2,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),M,
     &                    C0,WM,M)
C..................................................................... G
C
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                    GQ(1,1,IQ),M,C0,W2,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),M,
     &                    C0,WG,M)
C-----------------------------------------------------------------------
C
               NERRT = 0
               NERRM = 0
               NERRG = 0
               DO J = 1,NKM
                  DO I = 1,NKM
                     IF ( ABS(W1(I,J)-TAUQ(I,J,IQP)).GT.TOLR ) THEN
                        NERRT = NERRT + 1
                        WRITE (6,99002) ISYM,IQP,I,J,'TAU(R k)   ',
     &                                  TAUQ(I,J,IQP),'R TAU(k) R+',
     &                                  W1(I,J),W1(I,J) - TAUQ(I,J,IQP)
                     END IF
                     IF ( ABS(WM(I,J)-MQP(I,J,IQP)).GT.TOLR ) THEN
                        NERRM = NERRM + 1
                        WRITE (6,99002) ISYM,IQP,I,J,'  m(R k)   ',
     &                                  MQP(I,J,IQP),'R   m(k) R+',
     &                                  WM(I,J),WM(I,J) - MQP(I,J,IQP)
                     END IF
                     IF ( ABS(WG(I,J)-GQP(I,J,IQP)).GT.TOLR ) THEN
                        NERRG = NERRG + 1
                        WRITE (6,99002) ISYM,IQP,I,J,'  G(R k)   ',
     &                                  GQP(I,J,IQP),'R   G(k) R+',
     &                                  WG(I,J),WG(I,J) - GQP(I,J,IQP)
                     END IF
                  END DO
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C
            NERRTS = NERRTS + NERRT
            NERRMS = NERRMS + NERRM
            NERRGS = NERRGS + NERRG
            NERRTSS = NERRTSS + NERRTS
            NERRMSS = NERRMSS + NERRMS
            NERRGSS = NERRGSS + NERRGS
C
            IF ( NERRTS+NERRMS+NERRGS.GT.0 ) THEN
               STRERR = 'ERROR'
            ELSE
               STRERR = 'OK   '
            END IF
            WRITE (6,99009) ISYM,SYMUNITARY(ISYM),'K',KVECP,NERRTS,
     &                      NERRMS,NERRGS,STRERR
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      WRITE (6,99014) 1,NERRTSS + NERRMSS + NERRGSS
C
C=======================================================================
C                               PART  2
C=======================================================================
C
      TAUQA(1:NKM,1:NKM,1:NQ) = C0
      NERR2 = 0
      WRITE (6,99001) 2,'   check integration of TAU - matrix'
C
      WRITE (6,99004)
      IF ( KMROT.EQ.0 ) THEN
C
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         CALL CINIT(M*M*NQ,TAUQA)
C
         IK = 0
         DO ISYM = 1,NSYM
            IF ( SYMCRYSYS(ISYM) ) THEN
C
               CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,KVEC,1,0D0,
     &                    KVECP,1)
C
               IK = IK + 1
               WRITE (6,99012) 'A',IK,KVECP
C
               CALL STRSET(IK,KVECP,MAUX,TAUK,P)
C
               CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
               CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
               CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
               DO IQ = 1,NQ
                  I1 = IND0Q(IQ) + 1
                  DO J = 1,NKMQ(IQ)
                     J1 = IND0Q(IQ) + J
                     CALL ZAXPY(NKMQ(IQ),C1,TAUK(I1,J1),1,TAUQA(1,J,IQ),
     &                          1)
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
C
         TAUQA(1:NKM,1:NKM,1:NQ) = TAUQA(1:NKM,1:NKM,1:NQ)*SCLA
C
         WRITE (6,*) ' '
C
         IF ( SYMTAUTET_AVAILABLE ) THEN
C
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            CALL CINIT(M*M*NQ,TAUQB)
C
            IK = 0
            DO IWEDGE = 1,NWEDGE
C
               IROT = IWEDGEROT(IWEDGE)
C
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KVEC,1,0D0,
     &                    KVECP,1)
C
               IK = IK + 1
               WRITE (6,99012) 'B',IK,KVECP
C
               CALL STRSET(IK,KVECP,MAUX,TAUK,P)
C
               CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
               CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
               CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C------------------------------------------------------------ store TAUQ
               I = 0
               DO I1 = 1,NELMT
                  CSUM = C0
                  DO I3 = 1,NTAUUV(I1)
                     I = I + 1
                     CSUM = CSUM + WTAUUV(I)*TAUK(UTAUUV(I),VTAUUV(I))
                  END DO
C
                  TAUQB(ITBZ(I1),JTBZ(I1),QTBZ(I1))
     &               = TAUQB(ITBZ(I1),JTBZ(I1),QTBZ(I1)) + CSUM
               END DO
C
C-----------------------------------------------------------------------
            END DO
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
            SCLB = 1D0/DBLE(NWEDGE)
C
            DO IQ = 1,NQ
C
               WRITE (6,99013) 2,IQ
C
               DO J = 1,NKMQ(IQ)
                  DO I = 1,NKMQ(IQ)
                     TAUQB(I,J,IQ) = TAUQB(I,J,IQ)*SCLB - TAUQA(I,J,IQ)
                     IF ( ABS(TAUQB(I,J,IQ)).GT.TOLM ) NERR2 = NERR2 + 1
                  END DO
               END DO
C
               CALL CMATSTRUCT('TAU(IQ)  run A   ',TAUQA(1,1,IQ),
     &                         NKMQ(IQ),M,IREL,IREL,0,TOLM,6)
C
               IF ( KMROT.EQ.0 ) THEN
                  CALL CMATSTRUCT('TAU(IQ) run B - A  SHOULD BE BLANK',
     &                            TAUQB(1,1,IQ),NKMQ(IQ),M,IREL,IREL,0,
     &                            TOLM,6)
C
               ELSE
                  NERR2 = 9999
               END IF
            END DO
C
            WRITE (6,99014) 2,NERR2
C
         END IF
C
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'KMROT not allowed')
C
      END IF
C
C
      IF ( .NOT.SYMTAUTET_AVAILABLE ) WRITE (6,99005)
C
C=======================================================================
C                               PART  3
C=======================================================================
C
      WRITE (6,99001) 3,'   check integration of TAU - matrix'
      WRITE (6,99006)
C
      NERR3 = 0
      CALL CINIT(M*M*NQ,TAUQB)
      CALL CINIT(M*M*NQ,TAUQC)
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IK = 0
      DO IWEDGE = 1,NWEDGE
C
         IROT = IWEDGEROT(IWEDGE)
C
         CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,KVEC,1,0D0,KVECP,1)
C
         IK = IK + 1
         WRITE (6,99012) 'C',IK,KVECP
C
         CALL STRSET(IK,KVECP,MAUX,TAUK,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
         CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C------------------------------------------------------------ store TAUQ
         DO IQ = 1,NQ
            I1 = IND0Q(IQ) + 1
            DO J = 1,NKMQ(IQ)
               J1 = IND0Q(IQ) + J
               CALL ZAXPY(NKMQ(IQ),C1,TAUK(I1,J1),1,TAUQB(1,J,IQ),1)
            END DO
         END DO
C-----------------------------------------------------------------------
      END DO
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
            IF ( SYMUNITARY(ISYM) ) THEN
               CNT = 'N'
            ELSE
               CNT = 'T'
            END IF
C
            DO IQ = 1,NQ
C
               IQP = IQORGQP(ISYM,IQ)
C-----------------------------------------------------------------------
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                    TAUQB(1,1,IQP),M,C0,W1,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W1,M,DROT(1,1,ISYM),M,
     &                    C1,TAUQC(1,1,IQ),M)
C.......................................................................
            END DO
         END IF
      END DO
C
      SCLC = 1D0/DBLE(NWEDGE*NSYMACCEPTED)
C
      DO IQ = 1,NQ
C
         WRITE (6,99013) 3,IQ
C
         DO J = 1,NKMQ(IQ)
            DO I = 1,NKMQ(IQ)
               TAUQC(I,J,IQ) = TAUQC(I,J,IQ)*SCLC - TAUQA(I,J,IQ)
               TA = TAUQA(I,J,IQ)
               TC = TAUQC(I,J,IQ)
               IF ( ABS(TA).GT.TOL3 ) THEN
                  IF ( ABS(TC/TA).GT.TOL3 ) NERR3 = NERR3 + 1
               ELSE IF ( ABS(TC).GT.TOL3 ) THEN
                  NERR3 = NERR3 + 1
               END IF
C
            END DO
         END DO
C
         CALL CMATSTRUCT('TAU(IQ) run C - A  SHOULD BE BLANK',
     &                   TAUQC(1,1,IQ),NKMQ(IQ),M,IREL,IREL,0,TOL3,6)
C
      END DO
C
C
      WRITE (6,99014) 3,NERR3
C
C=======================================================================
C                               PART  4
C=======================================================================
C
      WRITE (6,99001) 4,'   check integration of TAU - matrix'
      WRITE (6,99007)
C
      IE = 1
      NERR4 = 0
      NKPTS0 = 100
C
      CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      IF ( ALLOCATED(KTAB) ) DEALLOCATE (KTAB,WKTAB)
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB')
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
C----------------------------------- start with same MSSQ in case of CPA
      MSSQ(:,:,:) = MSSQ0(:,:,:)
      TSSQ(:,:,:) = TSSQ0(:,:,:)
C
      CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                     ICPACONV,CONC,NOQ,ITOQ,PHASK,IE,NTMAX,TSST,
     &                     MSST,TSSQ,MSSQ,TAUQA)
C
      MSSQ_4A(:,:,:) = MSSQ(:,:,:)
C
      WRITE (6,99008)
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
C----------------------------------------------------- suppress symmetry
      NSYM = NSYMDUM
      NSYMACCEPTED_SAV = NSYMACCEPTED
      NSYMACCEPTED = NSYMDUM
C
C----------------------------------- start with same MSSQ in case of CPA
      MSSQ(:,:,:) = MSSQ0(:,:,:)
      TSSQ(:,:,:) = TSSQ0(:,:,:)
C
      CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                     ICPACONV,CONC,NOQ,ITOQ,PHASK,IE,NTMAX,TSST,
     &                     MSST,TSSQ,MSSQ,TAUQB)
C
      MSSQ_4B(:,:,:) = MSSQ(:,:,:)
C
      NSYM = NSYM_SAV
      NSYMACCEPTED = NSYMACCEPTED_SAV
C
      DO IQ = 1,NQ
C
         WRITE (6,99013) 4,IQ
C
         DO J = 1,NKMQ(IQ)
            DO I = 1,NKMQ(IQ)
               TAUQB(I,J,IQ) = TAUQB(I,J,IQ) - TAUQA(I,J,IQ)
               IF ( ABS(TAUQB(I,J,IQ)).GT.TOLM ) NERR4 = NERR4 + 1
C               IF ( ABS(TAUQB(I,J,IQ)).GT.TOLM )
C     & TAUQB(I,J,IQ) = TAUQB(I,J,IQ) / TAUQA(I,J,IQ)
            END DO
         END DO
C
         CALL CMATSTRUCT('TAU(IQ)  run A   ',TAUQA(1,1,IQ),NKMQ(IQ),M,
     &                   IREL,IREL,0,TOLM,6)
C
         CALL CMATSTRUCT('TAU(IQ) run B - A  SHOULD BE BLANK',
     &                   TAUQB(1,1,IQ),NKMQ(IQ),M,IREL,IREL,0,TOLM,6)
C
C
         IF ( ICPA(IQ).NE.0 ) THEN
C
            DO J = 1,NKMQ(IQ)
               DO I = 1,NKMQ(IQ)
                  MSSQ_4B(I,J,IQ) = MSSQ_4B(I,J,IQ) - MSSQ_4A(I,J,IQ)
                  IF ( ABS(MSSQ_4B(I,J,IQ)).GT.TOLM ) NERR4 = NERR4 + 1
C               IF ( ABS(MSSQ_4B(I,J,IQ)).GT.TOLM )
C     & MSSQ_4B(I,J,IQ) = MSSQ_4B(I,J,IQ) / MSSQ_4A(I,J,IQ)
               END DO
            END DO
C
            CALL CMATSTRUCT('mss(IQ)  run A   ',MSSQ_4A(1,1,IQ),NKMQ(IQ)
     &                      ,M,IREL,IREL,0,TOLM,6)
C
            CALL CMATSTRUCT('mss(IQ) run B - A  SHOULD BE BLANK',
     &                      MSSQ_4B(1,1,IQ),NKMQ(IQ),M,IREL,IREL,0,TOLM,
     &                      6)
C
         END IF
C
      END DO
C
      WRITE (6,99014) 4,NERR4
C
C=======================================================================
C                               PART  5
C      check the relation   TAU(q',q') = R TAU(q,q) R+   for q' = R q
C=======================================================================
C
C----------------------------------- start with same MSSQ in case of CPA
      MSSQ(:,:,:) = MSSQ0(:,:,:)
      TSSQ(:,:,:) = TSSQ0(:,:,:)
C
      WRITE (6,99001) 5,'check symmetry of site-diagonal TAU(q,q)'
      WRITE (6,99022)
C
      NERR5 = 0
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            DO IQP = 1,NQ
C
               IQ = IQORGQP(ISYM,IQP)
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C-----------------------------------------------------------------------
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                    TAUQA(1,1,IQ),M,C0,W2,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),M,
     &                    C0,W1,M)
C-----------------------------------------------------------------------
C
               DO J = 1,NKM
                  DO I = 1,NKM
                     IF ( ABS(W1(I,J)-TAUQA(I,J,IQP)).GT.TOLR ) THEN
                        NERR5 = NERR5 + 1
                        WRITE (6,99002) ISYM,IQP,I,J,'TAU(q'',q'')   ',
     &                                  TAUQA(I,J,IQP),'R TAU(q,q) R+',
     &                                  W1(I,J),W1(I,J) - TAUQA(I,J,IQP)
                     END IF
                  END DO
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      WRITE (6,99014) 5,NERR5
C
C=======================================================================
C
      WRITE (6,99015) IREL,NQ,KMROT,NSYM,NSYMACCEPTED,NWEDGE
      IF ( KMROT.EQ.0 ) WRITE (6,99016) NELMT
      WRITE (6,99017) NERRTSS,NERRMSS,NERRGSS
      IF ( KMROT.EQ.0 ) WRITE (6,99018) NERR2
      WRITE (6,99020) NERR3,NERR4,NERR5
      IF ( .NOT.SYMTAUTET_AVAILABLE ) WRITE (6,99019)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) THEN
         WRITE (IFILBUILDBOT,99025) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'PART 1',REAL(NERRTSS),REAL(NERRMSS)
     &                              ,REAL(NERRGSS)
         IF ( KMROT.EQ.0 ) WRITE (IFILBUILDBOT,99026)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),'PART 2',
     &                            REAL(NERR2)
         WRITE (IFILBUILDBOT,99026) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'PART 3',REAL(NERR3)
         WRITE (IFILBUILDBOT,99026) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'PART 4',REAL(NERR4)
         WRITE (IFILBUILDBOT,99026) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'PART 5',REAL(NERR5)
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (MAUX,TAUK,IPIV,TAUQA,TAUQB,TAUQC,W1,W2,GQ,MQ,WG)
      DEALLOCATE (WM,GQP,MQP,KTAB,WKTAB)
C
      STOP 'TEST 2 completed'
C
C=======================================================================
99001 FORMAT (//,1X,79('='),/,28X,'<SYMCHECK_STD>   PART',I3,/,20X,A,/,
     &        1X,79('='),/)
99002 FORMAT (I3,' Q',I3,2X,2I3,3X,A,2E13.5,/,19X,A,2E13.5,2F20.12)
99003 FORMAT (3X,'ISYM ',I3,' unitary ',L1,3X,A,1X,3F5.1,:,2X,A,1X,
     &        3F5.1,:,2X,A,1X,3F5.1)
99004 FORMAT (/,10X,'run A:  run through ALL equivalent ->k - vectors ',
     &        /,18X,'that are compatible with the crystal system',//,
     &        10X,
     &     'run B:  run through all NWEDGE NON-equivalent ->k - vectors'
     &     ,/,18X,'and use   TAU = SUM W(->k) * TAUK(->k)',/,18X,
     &     'as it is done in the tetrahedron integration method '/)
99005 FORMAT (/,10X,'symmetry information for ',/,
     &        'tetrahedron integration method not available',/,10X,
     &        'this part has been skipped ',/)
99006 FORMAT (/,10X,
     &     'run C:  run through all NWEDGE NON-equivalent ->k - vectors'
     &     ,/,18X,'and use   TAU = SUM(SYM) R  [ SUM(k) TAUK(->k)] R+',
     &     /,18X,
     &     'as it is done in the special points integration method '/)
99007 FORMAT (/,10X,'run A:  run <TAU_STD_KPOINTS> ',
     &        'making use of symmetry')
99008 FORMAT (/,10X,'run B:  run <TAU_STD_KPOINTS> ',
     &        'WITHOUT making use of symmetry  ---  NSYM = 0')
99009 FORMAT (3X,'ISYM ',I3,' unitary ',L1,3X,A,1X,3F5.1,:,4X,'T:',I4,
     &        '   m:',I4,'   G:',I4,3X,A)
99010 FORMAT (/,5X,'action of the symmetry operations  R  on:',/,7X,
     &        '- magnetisation vector ->M',/,7X,
     &        '- vector ->R in real space',/,7X,
     &        '- vector ->k in reciprocal space',/,5X,
     &        'antiunitary rotations of ->k include inversion I',/)
99011 FORMAT (/,5X,'action of the symmetry operations  R  on:',/,7X,
     &        '(T): scattering path operator TAU(->k)   ',5X,
     &        'TAU(Rk) = R TAU(k) R+',/,7X,
     &        '(m): single site m-matrix m_q            ',5X,
     &        '      m = R m R+',/,7X,
     &        '(G): KKR structure constant matrix G(->k)',5X,
     &        '  G(Rk) = R G(k) R+',/,5X,'for arbitrary vector ->k = (',
     &        2(F5.2,','),F5.2,' )',/)
99012 FORMAT (10x,'k-loop    run   ',A,I4,' ->k ',3F8.4)
99013 FORMAT (//,14x,53('='),/,14x,I2,
     &        ' check integration of TAU - matrix    for  IQ =',I3,/,
     &        14x,53('='),/)
99014 FORMAT (/,10X,'<SYMCHECK_STD>  finished PART',I2,:,5X,
     &        'ERROR code ',3I5)
99015 FORMAT (//,14X,53('='),/,34X,'SUMMARY',/,14X,53('='),//,10X,
     &        'run parameters:',//,10X,'IREL             ',I5,/,10X,
     &        'NQ               ',I5,/,10X,'KMROT            ',I5,/,10X,
     &        'NSYM             ',I5,/,10X,'NSYMACCEPTED     ',I5,/,10X,
     &        'NWEDGE           ',I5)
99016 FORMAT (10X,'NELMT            ',I5)
99017 FORMAT (/,10X,'ERROR codes:',//,10X,
     &        'PART 1  TAU(Rk) = R TAU(k) R+',24X,3I5)
99018 FORMAT (10X,'PART 2  TAU = SUM W(->k) * TAUK(->k)',17X,I5)
99019 FORMAT (10X,'PART 2  has been skipped -- SYM info not available')
99020 FORMAT (10X,'PART 3  TAU = SUM(SYM) R [ SUM(k) TAUK(->k)] R+',6X,
     &        I5,/,10X,'PART 4  run <TAU_STD_KPOINTS> making use of',
     &        ' symmetry ',I5,/,10X,'PART 5  TAU(q'',q'') =',
     &        ' R TAU(q,q) R+','   for q'' = R q',4X,I5,/)
99021 FORMAT (//,1X,79('*'),/,33X,'<SYMCHECK_STD>',/,1X,79('*'),//,21X,
     &        'BZ-integral   Int d^3k  TAU(k)[L1,L2] ')
99022 FORMAT (/,10X,'check the relation     TAU(q'',q'') = ',
     &        'R TAU(q,q) R+   for q'' = R q',/)
99023 FORMAT (10X,A,10I5)
99024 FORMAT (//,14x,53('='),/,14x,
     &        '  compare symmetrized    m - matrix    for  IQ =',I3,/,
     &        14x,53('='),//,10X,
     &        'number of local symmetry operations   ',I5,/,10X,
     &        'maximum difference  | m - <m> |  found',F20.12)
99025 FORMAT ('# BUILDBOT: ',A,': ',A,/,(1PE22.14),/,(1PE22.14),/,
     &        (1PE22.14))
99026 FORMAT ('# BUILDBOT: ',A,': ',A,/,(1PE22.14))
      END
