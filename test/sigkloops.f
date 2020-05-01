C*==sigkloops.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGKLOOPS(ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQBZ,MSSQA,
     &                     MSSQB)
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
     &    ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,JTTX,WTTJ,NTKTKLIN
      USE MOD_KSPACE,ONLY:WKTAB,NKTAB,KTAB,WKSUM
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM,WKM1
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_CALCMODE,ONLY:IREL
c modified by XJQ: parallel on k-point
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq
      IMPLICIT NONE
C*--SIGKLOOPS34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGKLOOPS')
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
     &        JTT,K14I,K23J,K32J,K41I,L1,L2,L3,L4,M,N,NKMSQ,Z2
      REAL*8 WK
c modified by XJQ: parallel on ik
      logical lparalleled, lopen
      integer level, ik0, ik1, nkp, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      real*8 t0_sigk, t1_sigk, wksum_tmp
      complex*16, dimension(:,:), allocatable :: tktkqq_tmp
      complex*16, dimension(:,:,:), allocatable :: sumqa_tmp, sumqb_tmp
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUKLINA,TAUKLINB,SUMQA,SUMQB,TKTKQQ,IPIV
      ALLOCATABLE IPROCK,CWORK
C
      CALL TRACK_INFO(ROUTINE)
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
c modified by XJQ: parallel on k-point
      level=toplevel+1
c
      call get_comm_level(level,'parent',lparalleled,
     &                    parent_comm,nprocs_parent,parent_rank)
      if(parent_rank==0) call cpu_time(t0_sigk)
c
c      CALL MPI_DISTRIBUTE(IPROCK,NKTAB,MPI_KLOOP,'K')
      call mpi_multilevel_distribute(routine,level,nktab,
     &                               lparalleled,ik0,ik1)
c
      allocate(sumqa_tmp(nkmmax,nkmmax,nqmax))
      allocate(tktkqq_tmp(ntktklin,nz12))
      if(.not. eryda_eq_erydb) allocate(sumqb_tmp(nkmmax,nkmmax,nqmax))
c end-mod-xjq
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C ------------ calculate energy - dependent terms of structure constants
C
      CALL STRCC(ERYDA,.FALSE.)
C
c modified by XJQ: parallel on k-point
c      TKTKQQ(1:NTKTKLIN,1:NZ12) = C0
      tktkqq_tmp(1:ntktklin,1:nz12) = c0
c end-mod-xjq
C
C------------------------- assume the same L-expansion for every site IQ
      NKM = NKMQ(1)
      DO IQ = 2,NQ
         IF ( NKMQ(IQ).NE.NKM )
     &         CALL STOP_MESSAGE(ROUTINE,'NKMQ(IQ)<>NKM')
      END DO
C
c modified by XJQ: parallel on k-point
c      SUMQA(:,:,:) = C0
c      IF ( .NOT.ERYDA_EQ_ERYDB ) SUMQB(:,:,:) = C0
      sumqa_tmp(:,:,:) = c0
      if(.not. eryda_eq_erydb) sumqb_tmp(:,:,:) = c0
c end-mod-xjq
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ: parallel on k-point
c         IF ( MPI_ID.NE.IPROCK(IK) .AND. .NOT.MPI_ELOOP ) CYCLE
         if((ik < ik0 .or. ik > ik1) .and. .not. mpi_eloop) cycle
         write(*,*) 'ik=',ik,',ik0=',ik0,',ik1=',ik1
c end-mod-xjq
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
c modified by XJQ: parallel on k-point
c               CALL ZAXPY(N,CWK,TAUKLINA(I1),1,SUMQA(1,J,IQ),1)
               CALL ZAXPY(N,CWK,TAUKLINA(I1),1,SUMQA_TMP(1,J,IQ),1)
c end-mod-xjq
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
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLINB,MSSQB,NQMAX,NKKR,NKMMAX)
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
c modified by XJQ: parallel on k-point
c                  CALL ZAXPY(N,CWK,TAUKLINB(I1),1,SUMQB(1,J,IQ),1)
                  CALL ZAXPY(N,CWK,TAUKLINB(I1),1,SUMQB_TMP(1,J,IQ),1)
c end-mod-xjq
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
c modified by XJQ: parallel on k-point
c                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                  TKTKQQ_TMP(I,1) = TKTKQQ_TMP(I,1) + CSUM1*WK
c end-mod-xjq  
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
c modified by XJQ: parallel on k-point
c                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
c                  TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
                  TKTKQQ_TMP(I,1) = TKTKQQ_TMP(I,1) + CSUM1*WK
                  TKTKQQ_TMP(I,NZ12) = TKTKQQ_TMP(I,NZ12) + CSUMX*WK
c end-mod-xjq
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
c modified by XJQ: parallel on k-point
c               TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
               TKTKQQ_TMP(I,1) = TKTKQQ_TMP(I,1) + CSUM1*WK
c end-mod-xjq
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
c modified by XJQ: parallel on k-point
c               TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
c               TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
               TKTKQQ_TMP(I,1) = TKTKQQ_TMP(I,1) + CSUM1*WK
               TKTKQQ_TMP(I,NZ12) = TKTKQQ_TMP(I,NZ12) + CSUMX*WK
c end-mod-xjq
            END DO
C-----------------------------------------------------------------------
C
         END IF
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ: parallel on k-point
c      IF ( MPI_KLOOP ) THEN
cC
c         CALL DRV_MPI_BARRIER
cC
cC        use TAUQA and TAUQB as work space for transfer
c         M = NKMMAX*NKMMAX*NQMAX
cC
c         CALL DRV_MPI_REDUCE_C(SUMQA(1,1,1),TAUQA(1,1,1),M)
c         IF ( .NOT.ERYDA_EQ_ERYDB ) CALL DRV_MPI_REDUCE_C(SUMQB(1,1,1),
c     &        TAUQB(1,1,1),M)
cC
c         M = NTKTKLIN*NZ12
c         ALLOCATE (CWORK(NTKTKLIN,NZ12))
c         CALL DRV_MPI_REDUCE_C(TKTKQQ(1,1),CWORK(1,1),M)
c         DEALLOCATE (CWORK)
cC
c      END IF
      if(.not. lparalleled) then
        sumqa(:,:,:) = sumqa_tmp(:,:,:)
        tktkqq(:,:) = tktkqq_tmp(:,:)
        if(.not. eryda_eq_erydb) sumqb(:,:,:) = sumqb_tmp(:,:,:)
      else
        call mpi_barrier(parent_comm,mpierr)
        call get_comm_level(level,'inter ',lparalleled,
     &                      inter_comm,nprocs_inter,inter_rank)
c
        call mpi_reduce(sumqa_tmp,sumqa,nkmmax*nkmmax*nqmax,
     &                  mpi_double_complex,mpi_sum,0,inter_comm,mpierr)
        call mpi_reduce(tktkqq_tmp,tktkqq,ntktklin*nz12,
     &                  mpi_double_complex,mpi_sum,0,inter_comm,mpierr)
        if(.not. eryda_eq_erydb)
     &    call mpi_reduce(sumqb_tmp,sumqb,nkmmax*nkmmax*nqmax,
     &                    mpi_double_complex,mpi_sum,
     &                    0,inter_comm,mpierr)
        call mpi_multilevel_free(level)
      endif
      if(parent_rank==0) then
        call cpu_time(t1_sigk)
        write(*,*)
        write(*,*) 'part1 of sigkloops takes ',
     &              t1_sigk-t0_sigk,' seconds'
        write(*,*)
      endif
      deallocate(sumqa_tmp)
      deallocate(tktkqq_tmp)
      if(.not. eryda_eq_erydb) deallocate(sumqb_tmp)
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c XJQ: since chiz is allocated only for mpi_id 0, we cannot parallel the part below
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
            IF ( NZ12.EQ.2 .OR. .NOT.ERYDA_EQ_ERYDB ) THEN
C
               DO Z2 = 1,NZ12
C
C
                  IF ( IQ.EQ.JQ ) THEN
                     CHIZ(K14I,K23J,Z2) = TKTKQQ(I,Z2) - TAUQA(L1,L2,IQ)
     &                  *TAUQBZ(L3,L4,IQ,Z2)
                  ELSE
                     CHIZ(K14I,K23J,Z2) = TKTKQQ(I,Z2)
                  END IF
               END DO
C
C-----------------------------------------------------------------------
            ELSE
C
               CHIZ(K14I,K23J,1) = TKTKQQ(I,1)
C
               K32J = (L3-1)*NKM + L2 + (JQ-1)*NKMSQ
               K41I = (L4-1)*NKM + L1 + (IQ-1)*NKMSQ
C
               CHIZ(K32J,K41I,1) = TKTKQQ(I,1)
C
               IF ( IQ.EQ.JQ ) THEN
C
                  CHIZ(K14I,K23J,1) = CHIZ(K14I,K23J,1)
     &                                - TAUQA(L1,L2,IQ)*TAUQA(L3,L4,IQ)
C
                  IF ( L1.NE.L3 .OR. L2.NE.L4 ) CHIZ(K32J,K41I,1)
     &                 = CHIZ(K32J,K41I,1) - TAUQA(L1,L2,IQ)
     &                 *TAUQA(L3,L4,IQ)
C
               END IF
C
            END IF
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
