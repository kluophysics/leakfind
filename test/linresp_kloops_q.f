C*==linresp_kloops_q.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_KLOOPS_Q(ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQBZ,
     &                            MSSQA,MSSQB)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array                    *
C   *------------------------------------------------------------------*
C   *  WKSUM is determined in <INIT_MOD_KSPACE>                        *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,CHIZ,NZ12,NZ12MAX,ITTA,ITTB,
     &    ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,JTTX,WTTJ,NTKTKLIN,
     &    QVEC_PERT
      USE MOD_KSPACE,ONLY:WKTAB,NKTAB,KTAB,WKSUM
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM,WKM1
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--LINRESP_KLOOPS_Q35
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_KLOOPS_Q')
C
C Dummy arguments
C
      COMPLEX*16 ERYDA,ERYDB,PA,PB
      COMPLEX*16 MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           TAUQA(NKMMAX,NKMMAX,NQMAX),TAUQB(NKMMAX,NKMMAX,NQMAX),
     &           TAUQBZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX)
C
C Local variables
C
      CHARACTER*1 CNT
      COMPLEX*16 CSCL,CSUM1,CSUMX,CWK,CWORK(:,:),MAUX(:,:),SUMQA(:,:,:),
     &           SUMQB(:,:,:),TAUKLINA(:),TAUKLINB(:),TKTKQQ(:,:),WT
      INTEGER I,I1,IA_ERR,IK,INFO,IPIV(:),IPROCK(:),IQ,IQP,ISYM,J,JQ,
     &        JTT,K14I,K23J,L1,L2,L3,L4,M,N,NKMSQ,Z2
      REAL*8 WK
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUKLINA,TAUKLINB,SUMQA,SUMQB,TKTKQQ,IPIV
      ALLOCATABLE IPROCK,CWORK
C
      ALLOCATE (SUMQA(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TKTKQQ(NTKTKLIN,NZ12),IPIV(NKKR))
C
      M = NKKR
      IF ( .NOT.ERYDA_EQ_ERYDB )
     &     ALLOCATE (TAUKLINB(M*M),SUMQB(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MAUX(M,M),TAUKLINA(M*M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUKLINA')
C
      IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL = 2')
      IF ( NZ12.LT.1 .OR. NZ12.GT.2 )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12 not set properly')
      IF ( NZ12.NE.NZ12MAX )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12MAX not set properly')
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., NKTAB
C      over processors;   IK=NKTAB  dealt with by IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCK(NKTAB))
C
      CALL MPI_DISTRIBUTE(IPROCK,NKTAB,MPI_KLOOP,'K')
C
C ------------ calculate energy - dependent terms of structure constants
C
      CALL STRCC(ERYDA,.FALSE.)
C
      TKTKQQ(1:NTKTKLIN,1:NZ12) = C0
C
C------------------------- assume the same L-expansion for every site IQ
      NKM = NKMQ(1)
      DO IQ = 2,NQ
         IF ( NKMQ(IQ).NE.NKM )
     &         CALL STOP_MESSAGE(ROUTINE,'NKMQ(IQ)<>NKM')
      END DO
C
      SUMQA(:,:,:) = C0
      IF ( .NOT.ERYDA_EQ_ERYDB ) SUMQB(:,:,:) = C0
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCK(IK) .OR. MPI_ELOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            WK = WKTAB(IK)
            CWK = DCMPLX(WK,0D0)
C
C-------------------------------------------------------- calculate TAUA
C
            IF ( .NOT.ERYDA_EQ_ERYDB ) CALL STRCC(ERYDA,.FALSE.)
C
            CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLINA,PA)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLINA,MSSQA,NQMAX,NKKR,NKMMAX)
C
            CALL ZGETRF(NKKR,NKKR,TAUKLINA,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUKLINA,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C----------------------------------------------------------- store TAUQA
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               I1 = IND0Q(IQ)*NKKR + IND0Q(IQ) + 1 - NKKR
               N = NKMQ(IQ)
               DO J = 1,N
                  I1 = I1 + NKKR
                  CALL ZAXPY(N,CWK,TAUKLINA(I1),1,SUMQA(1,J,IQ),1)
               END DO
            END DO
C
            IF ( .NOT.(ERYDA_EQ_ERYDB) ) THEN
C
C-------------------------------------------------------- calculate TAUB
C
               CALL STRCC(ERYDB,.FALSE.)
C
               CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLINB,PB)
C
               CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLINB,MSSQB,NQMAX,NKKR,
     &                     NKMMAX)
C
               CALL ZGETRF(NKKR,NKKR,TAUKLINB,NKKR,IPIV,INFO)
               CALL ZGETRI(NKKR,TAUKLINB,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C----------------------------------------------------------- store TAUQB
C
               DO IQ = IQBOT_CHI,IQTOP_CHI
                  I1 = IND0Q(IQ)*NKKR + IND0Q(IQ) + 1 - NKKR
                  N = NKMQ(IQ)
                  DO J = 1,N
                     I1 = I1 + NKKR
                     CALL ZAXPY(N,CWK,TAUKLINB(I1),1,SUMQB(1,J,IQ),1)
                  END DO
               END DO
C
C------------------------------------------------------- store TAUA*TAUB
C
               IF ( NZ12.EQ.1 ) THEN
                  JTT = 0
                  DO I = 1,NTKTKLIN
                     CSUM1 = C0
                     DO J = 1,NTTJ(I)
                        JTT = JTT + 1
                        WT = WTTJ(JTT)*TAUKLINA(JTT1(JTT))
                        CSUM1 = CSUM1 + WT*TAUKLINB(JTT2(JTT))
                     END DO
                     TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                  END DO
               ELSE
                  JTT = 0
                  DO I = 1,NTKTKLIN
                     CSUM1 = C0
                     CSUMX = C0
                     DO J = 1,NTTJ(I)
                        JTT = JTT + 1
                        WT = WTTJ(JTT)*TAUKLINA(JTT1(JTT))
                        CSUM1 = CSUM1 + WT*TAUKLINB(JTT2(JTT))
                        CSUMX = CSUMX + WT*DCONJG(TAUKLINB(JTTX(JTT)))
                     END DO
                     TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                     TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
                  END DO
               END IF
C
C-----------------------------------------------------------------------
C                       ERYDA = ERYDB
C------------------------------------------------------- store TAUA*TAUA
C
            ELSE IF ( NZ12.EQ.1 ) THEN
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     WT = WTTJ(JTT)*TAUKLINA(JTT1(JTT))
                     CSUM1 = CSUM1 + WT*TAUKLINA(JTT2(JTT))
                  END DO
                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
               END DO
            ELSE
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  CSUMX = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     WT = WTTJ(JTT)*TAUKLINA(JTT1(JTT))
                     CSUM1 = CSUM1 + WT*TAUKLINA(JTT2(JTT))
                     CSUMX = CSUMX + WT*DCONJG(TAUKLINA(JTTX(JTT)))
                  END DO
                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                  TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
               END DO
C-----------------------------------------------------------------------
C
            END IF
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_KLOOP ) THEN
C
         CALL DRV_MPI_BARRIER
C
C        use TAUQA and TAUQB as work space for transfer
         M = NKMMAX*NKMMAX*NQMAX
C
         CALL DRV_MPI_REDUCE_C(SUMQA(1,1,1),TAUQA(1,1,1),M)
         IF ( .NOT.ERYDA_EQ_ERYDB ) CALL DRV_MPI_REDUCE_C(SUMQB(1,1,1),
     &        TAUQB(1,1,1),M)
C
         M = NTKTKLIN*NZ12
         ALLOCATE (CWORK(NTKTKLIN,NZ12))
         CALL DRV_MPI_REDUCE_C(TKTKQQ(1,1),CWORK(1,1),M)
         DEALLOCATE (CWORK)
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C
         TAUQA(:,:,:) = C0
         IF ( .NOT.ERYDA_EQ_ERYDB ) TAUQB(:,:,:) = C0
C
C-----------------------------------------------------------------------
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
               DO IQ = IQBOT_CHI,IQTOP_CHI
C
                  N = NKMQ(IQ)
                  IQP = IQORGQP(ISYM,IQ)
C-----------------------------------------------------------------------
                  CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),NKMMAX,
     &                       SUMQA(1,1,IQP),NKMMAX,C0,WKM1,NKMMAX)
C
                  CALL ZGEMM('N','C',N,N,N,C1,WKM1,NKMMAX,DROT(1,1,ISYM)
     &                       ,NKMMAX,C1,TAUQA(1,1,IQ),NKMMAX)
C
C.......................................................................
C
                  IF ( .NOT.ERYDA_EQ_ERYDB ) THEN
                     CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),NKMMAX,
     &                          SUMQB(1,1,IQP),NKMMAX,C0,WKM1,NKMMAX)
C
                     CALL ZGEMM('N','C',N,N,N,C1,WKM1,NKMMAX,
     &                          DROT(1,1,ISYM),NKMMAX,C1,TAUQB(1,1,IQ),
     &                          NKMMAX)
                  END IF
C-----------------------------------------------------------------------
               END DO
            END IF
         END DO
C
         CSCL = 1D0/DBLE(NSYMACCEPTED*WKSUM)
C
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO J = 1,NKMQ(IQ)
               CALL ZSCAL(NKMQ(IQ),CSCL,TAUQA(1,J,IQ),1)
            END DO
         END DO
C
         IF ( ERYDA_EQ_ERYDB ) THEN
C
            TAUQB(:,:,IQBOT_CHI:IQTOP_CHI)
     &         = TAUQA(:,:,IQBOT_CHI:IQTOP_CHI)
C
         ELSE
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               DO J = 1,NKMQ(IQ)
                  CALL ZSCAL(NKMQ(IQ),CSCL,TAUQB(1,J,IQ),1)
               END DO
            END DO
C
         END IF
C
         CSCL = 1D0/WKSUM
         DO Z2 = 1,NZ12
            CALL ZSCAL(NTKTKLIN,CSCL,TKTKQQ(1,Z2),1)
         END DO
C
C
C-----------------------------------------------------------------------
C                   set up the matrix   TAUQBZ
C-----------------------------------------------------------------------
         IF ( NZ12.EQ.1 ) THEN
C
            IF ( ERYDA_EQ_ERYDB ) THEN
C
               DO IQ = IQBOT_CHI,IQTOP_CHI
                  DO J = 1,NKMQ(IQ)
                     DO I = 1,NKMQ(IQ)
                        TAUQBZ(I,J,IQ,1) = TAUQA(I,J,IQ)
                     END DO
                  END DO
               END DO
C
            ELSE
C
               DO IQ = IQBOT_CHI,IQTOP_CHI
                  DO J = 1,NKMQ(IQ)
                     DO I = 1,NKMQ(IQ)
                        TAUQBZ(I,J,IQ,1) = TAUQB(I,J,IQ)
                     END DO
                  END DO
               END DO
C
            END IF
C
         ELSE IF ( ERYDA_EQ_ERYDB ) THEN
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               DO J = 1,NKMQ(IQ)
                  DO I = 1,NKMQ(IQ)
                     TAUQBZ(I,J,IQ,1) = TAUQA(I,J,IQ)
                     TAUQBZ(I,J,IQ,2) = DCONJG(TAUQA(J,I,IQ))
                  END DO
               END DO
            END DO
C
         ELSE
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               DO J = 1,NKMQ(IQ)
                  DO I = 1,NKMQ(IQ)
                     TAUQBZ(I,J,IQ,1) = TAUQB(I,J,IQ)
                     TAUQBZ(I,J,IQ,2) = DCONJG(TAUQB(J,I,IQ))
                  END DO
               END DO
            END DO
C
C
         END IF
C
C-----------------------------------------------------------------------
C               set up the super-matrix  CHI(K14I,k2)
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
         NKMSQ = NKM*NKM
C
         DO Z2 = 1,NZ12
            CALL CINIT(NKMMAX**4*NQMAX**2,CHIZ(1,1,Z2))
         END DO
C
         DO I = 1,NTKTKLIN
C
            L1 = ITTA(I)
            L2 = ITTB(I)
            L3 = ITTC(I)
            L4 = ITTD(I)
            IQ = ITTQ1(I)
            JQ = ITTQ2(I)
C
            K14I = (L1-1)*NKM + L4 + (IQ-1)*NKMSQ
            K23J = (L2-1)*NKM + L3 + (JQ-1)*NKMSQ
C
C-----------------------------------------------------------------------
C
            DO Z2 = 1,NZ12
C
               IF ( IQ.EQ.JQ ) THEN
                  CHIZ(K14I,K23J,Z2) = TKTKQQ(I,Z2) - TAUQA(L1,L2,IQ)
     &                                 *TAUQBZ(L3,L4,IQ,Z2)
               ELSE
                  CHIZ(K14I,K23J,Z2) = TKTKQQ(I,Z2)
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
C
         END DO
C
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_KLOOP ) CALL DRV_MPI_BARRIER
C
      DEALLOCATE (MAUX,TAUKLINA,SUMQA,TKTKQQ)
      IF ( .NOT.ERYDA_EQ_ERYDB ) DEALLOCATE (TAUKLINB,SUMQB)
C
      END
