C*==tau_std_kpoints.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                           ICPACONV,CONC,NOQ,ITOQ,PHASK,IE,NTMAX,
     &                           TSST,MSST,TSSQ,MSSQ,TAUQ)
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
C   *    - using LU-decomposition for matrix inversion                 *
C   *    - LLOYD = .TRUE.: the k-dependent terms in the Lloyd formula  *
C   *                      are calculated                              *
C   *                                                                  *
C   *   UPDATE 2012/2013                                               *
C   *    - paralellisation of the k-loop for   MPI_KLOOP = .T.         *
C   *    - import of type dependent variables via argument list        *
C   *      to allow calculations for vibrations etc                    *
C   *      via the CPA (alloy analogy model) using  "pseudo types"     *
C   *                                                                  *
C   *   UPDATE 2014                                                    *
C   *      old settings restored                                       *
C   *                                                                  *
C   *      real      CONC   ==  CONC        alloys                     *
C   *      alloys    NOQ    ==  NOQ         + vibrations               *
C   *                ITOQ   ==  ITOQ         + spin fluctuations ...   *
C   *                NTMAX  ==  NTMAX                                  *
C   *                MSST   ==  MSST                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_KSPACE,ONLY:WKTAB,NGFEP,GFEP,NKTAB,KTAB
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,ICPAALG,NCPA
      USE MOD_SITES,ONLY:NQ,NQMAX,ICPA,DROTQ
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR
      USE MOD_CALCMODE,ONLY:ITEST,LLOYD,IREL,KMROT
      USE MOD_CONSTANTS,ONLY:C0,PI
c modified by XJQ: parallel on k-point
c      USE MOD_MPI,ONLY:MPI_ID,MPI,MPI_KLOOP,MPI_ELOOP
      USE MOD_MPI,ONLY:NPROCS,MPI_ID,MPI_KLOOP,MPI_ELOOP
      use mod_files,only:iotmp
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq
      IMPLICIT NONE
C*--TAU_STD_KPOINTS48
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_STD_KPOINTS')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IE,IPRINT,ITCPA,NTMAX
      REAL*8 CONC(NTMAX)
      INTEGER ITOQ(NTMAX,NQMAX),NOQ(NQMAX)
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(IE),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 CPACORR,CPAERR,CPAERRL,GX,GY,GZ,KN2,WK,WKSUM
      COMPLEX*16 CWORK(:,:,:),DETK,DETKFE,DETKGT,EDU,MAUX(:,:),
     &           MSSQLU(:,:,:),SUMQ(:,:,:),TAUK(:,:),UDIA(:),W1(:,:)
      INTEGER I,I1,IA_ERR,IK,IPROC_K(:),IQ,J,M,N
c modified by XJQ: parallel on ik
      logical lparalleled, lopen, lexist
      integer level, ik0, ik1, nkp, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row, ikm
      real*8 t0_tau, t1_tau, wksum_tmp
      complex*16 :: detkfe_tmp, detkgt_tmp
      complex*16, dimension(:,:,:), allocatable :: sumq_tmp
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,MAUX,SUMQ,W1,MSSQLU,UDIA,CWORK,IPROC_K
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (MSSQLU(NKMMAX,NKMMAX,NQMAX),UDIA(NKKR))
      ALLOCATE (SUMQ(NKMMAX,NKMMAX,NQMAX),W1(NKMMAX,NKMMAX))
      ALLOCATE (TAUK(NKKR,NKKR),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUK')
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., NKTAB
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C------------------------------------ ensure proper setting of MPI_ELOOP
c      IF ( MPI .AND. .NOT.MPI_KLOOP ) MPI_ELOOP = .TRUE.
      IF ( NPROCS>1 .AND. .NOT.MPI_KLOOP ) MPI_ELOOP = .TRUE.
C
      ALLOCATE (IPROC_K(NKTAB))
C
c modified by XJQ: parallel on k-point
c if energies have been also paralleled using mpi_multilevel_distribute,
c setting lothermpi is not needed, but currently it is not yet done
      level=toplevel+1
      if(mpi_eloop) lothermpi=.true.
      call mpi_multilevel_distribute(routine,level,nktab,
     &                               lparalleled,ik0,ik1)
c
      if(lparalleled) then
        call get_comm_level(level,'parent',lparalleled,
     &                      parent_comm,nprocs_parent,parent_rank)
        call mpi_barrier(parent_comm,mpierr)
      endif
      call cpu_time(t0_tau)
      allocate(sumq_tmp(nkmmax,nkmmax,nqmax))
c end-mod-xjq
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      EDU = ERYD/(2*PI/ALAT)**2
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ICPACONV = 0
      CPAERRL = 1.0D+6
      ITCPA = 0
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      TAUQ(1:NKMMAX,1:NKMMAX,1:NQ) = C0
c modified by XJQ: parallel on k-point
c      SUMQ(1:NKMMAX,1:NKMMAX,1:NQ) = C0
      sumq_tmp(1:nkmmax,1:nkmmax,1:nq) = c0
c end-mod-xjq
C
      IF ( LLOYD ) THEN
C
         DO IQ = 1,NQ
            N = NKMQ(IQ)
            DO J = 1,N
               CALL ZCOPY(N,MSSQ(1,J,IQ),1,MSSQLU(1,J,IQ),1)
            END DO
            CALL CINVLU(MSSQLU(1,1,IQ),W1,N,NKMMAX)
         END DO
C
      END IF
C
c modified by XJQ: parallel on k-point
c      DETKGT = C0
c      DETKFE = C0
      detkgt_tmp = c0
      detkfe_tmp = c0
      PHASK(IE) = C0
c      WKSUM = 0D0
      wksum_tmp = 0d0
c end-mod-xjq
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ: parallel on k-point
c         IF ( MPI_ID.EQ.IPROC_K(IK) ) THEN
         if(ik >= ik0 .and. ik <=ik1) then
            if(lparalleled) write(*,*) 'mpi_id=',mpi_id,
     &        ',ik=',ik,',ik0=',ik0,',ik1=',ik1
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            CALL STRSET(IK,KTAB(1,IK),TAUK,MAUX,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,MAUX,MSSQ,NQMAX,NKKR,NKMMAX)
C
            WK = WKTAB(IK)
c modified by XJQ: parallel on k-point
c            WKSUM = WKSUM + WK
            wksum_tmp = wksum_tmp + wk
c end-mod-xjq
C
c modified by XJQ: parallel on k-point
c            CALL CINVDIABLK(NKKR,MAUX,UDIA,SUMQ,NKMMAX,NQ,IND0Q,NKMQ,
c     &                      .TRUE.,WK)
            CALL CINVDIABLK(NKKR,MAUX,UDIA,SUMQ_TMP,NKMMAX,NQ,IND0Q,
     &                      NKMQ,.TRUE.,WK)
c end-mod-xjq
C
C--------------------------------------------------------- LLOYD FORMULA
            IF ( LLOYD ) THEN
               DETK = C0
               DO IQ = 1,NQ
                  I1 = IND0Q(IQ)
                  N = NKMQ(IQ)
                  DO J = 1,N
                     I = I1 + J
                     DETK = DETK + LOG(UDIA(I)/MSSQLU(J,J,IQ))
                  END DO
               END DO
c modified by XJQ: parallel on k-point
c               DETKGT = DETKGT - WK*DETK
               DETKGT_TMP = DETKGT_TMP - WK*DETK
c end-mod-xjq
C
               DETK = C0
               DO I = 1,NGFEP
                  GX = GFEP(1,I)
                  GY = GFEP(2,I)
                  GZ = GFEP(3,I)
                  KN2 = (KTAB(1,IK)+GX)**2 + (KTAB(2,IK)+GY)
     &                  **2 + (KTAB(3,IK)+GZ)**2
                  DETK = DETK + LOG(KN2-EDU)
               END DO
C
c modified by XJQ: parallel on k-point
c               DETKFE = DETKFE - WK*DETK
               DETKFE_TMP = DETKFE_TMP - WK*DETK
c end-mod-xjq
            END IF
C-----------------------------------------------------------------------
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c      IF ( MPI_KLOOP ) THEN
c modified by XJQ: parallel on k-point
cC
c         CALL DRV_MPI_BARRIER
cC
c         M = NKMMAX*NKMMAX*NQMAX
c         ALLOCATE (CWORK(NKMMAX,NKMMAX,NQMAX))
c         CALL DRV_MPI_REDUCE_C(SUMQ(1,1,1),CWORK(1,1,1),M)
cC
c         CALL DRV_MPI_REDUCE_C(DETKFE,CWORK(1,1,1),1)
c         CALL DRV_MPI_REDUCE_C(DETKGT,CWORK(1,1,1),1)
c         CALL DRV_MPI_REDUCE_R(WKSUM,WK,1)
cC
c         DEALLOCATE (CWORK)
         
c      END IF
      if(.not. lparalleled) then
        sumq(:,:,:) = sumq_tmp(:,:,:)
        wksum = wksum_tmp
        detkfe = detkfe_tmp
        detkgt = detkgt_tmp
        lothermpi=.false.
      else
c collecting is safe, since mpi_reduce + mpi_bcase is used
        call mpi_barrier(parent_comm,mpierr)
        call get_comm_level(level,'inter ',lparalleled,
     &                      inter_comm,nprocs_inter,inter_rank)
c
        call mpi_reduce(sumq_tmp,sumq,nkmmax*nkmmax*nqmax,
     &                  mpi_double_complex,mpi_sum,0,inter_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_reduce(wksum_tmp,wksum,1,mpi_double_precision,
     &                  mpi_sum,0,inter_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_reduce(detkfe_tmp,detkfe,1,mpi_double_complex,
     &                  mpi_sum,0,inter_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_reduce(detkgt_tmp,detkgt,1,mpi_double_complex,
     &                  mpi_sum,0,inter_comm,mpierr)
      endif
c debug
      if(.false.) then
      inquire(iotmp,opened=lopen)
      if(lopen) stop 'iotmp is opened, try another for debugging'
      if(mpi_id==0) then
        open(iotmp,file='debug_sumq.dat',status='replace')
        do ikm=1,nkmmax
          write(iotmp,'(i4,6e15.8)') ikm,sumq(ikm,ikm,1),
     &      sumq(ikm,ikm,nq/2+1),sumq(ikm,ikm,nq)
        enddo
        close(iotmp)
      endif
      endif ! .false.
c debug
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C
         DETKFE = DETKFE/DBLE(WKSUM)
         DETKGT = DETKGT/DBLE(WKSUM)
C
         IF ( IREL.GE.3 ) THEN
            PHASK(IE) = DETKGT + 2*DETKFE
         ELSE
            PHASK(IE) = DETKGT + DETKFE
         END IF
C
         CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,WKSUM,SUMQ,TAUQ,W1,NQ,
     &                  NKMQ,DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                  NSYMACCEPTED,NQMAX,NKMMAX)
C
         IF ( NCPA.GT.0 ) THEN
C
            IF ( ICPAALG.EQ.1 ) THEN
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                       NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,
     &                       NTMAX,NQMAX,NKMMAX)
C
            ELSE IF ( ICPAALG.EQ.2 ) THEN
               CALL CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,
     &                        NKMQ,NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,
     &                        TAUQ,KMROT,DROTQ,NTMAX,NQMAX,NKMMAX)
            ELSE
               CALL CPABROYDEN(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,
     &                         NQ,NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUQ,
     &                         KMROT,DROTQ,NTMAX,NQMAX,NKMMAX)
            END IF
C
            CALL SYMSUMRTR(.FALSE.,.TRUE.,.FALSE.,1D0,SUMQ,MSSQ,W1,NQ,
     &                     NKMQ,DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,
     &                     NSYM,NSYMACCEPTED,NQMAX,NKMMAX)
C
            IF ( IPRINT.GE.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
            IF ( CPAERR.LE.CPATOL ) THEN
               ICPACONV = 1
               IF ( IPRINT.GT.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &              CPACHNG
            ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
               WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 1
            ELSE IF ( CPAERR.GT.20000*CPAERRL ) THEN
               WRITE (6,99003) ITCPA
               WRITE (6,99004) CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 2
            ELSE
               CPAERRL = CPAERR
            END IF
C
         END IF
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
c modified by XJQ: parallel on k-point
c      IF ( MPI_KLOOP ) THEN
c         CALL DRV_MPI_BARRIER
c         CALL DRV_MPI_BCAST_C(0,MSSQ(1,1,1),NKMMAX*NKMMAX*NQ)
c         CALL DRV_MPI_BCAST_C(0,TAUQ(1,1,1),NKMMAX*NKMMAX*NQ)
c         CALL DRV_MPI_BCAST_I(0,ICPAFLAG,1)
c         CALL DRV_MPI_BCAST_I(0,ICPACONV,1)        
c      END IF
c because the part above was run only on mpi_id 0, we can change in future
      call cpu_time(t1_tau)
      write(*,*) 'tau takes ',t1_tau-t0_tau,' s in mpi_id ',mpi_id
      if(lparalleled) then
        call mpi_barrier(parent_comm,mpierr)
        call mpi_bcast(mssq,nkmmax*nkmmax*nqmax,mpi_double_complex,
     &                 0,parent_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_bcast(tauq,nkmmax*nkmmax*nqmax,mpi_double_complex,
     &                 0,parent_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_bcast(icpaflag,1,mpi_int,0,parent_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_bcast(icpaconv,1,mpi_int,0,parent_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
      endif
c end-mod-xjq
C
      IF ( ICPAFLAG.EQ.0 .AND. ICPACONV.EQ.0 .AND. NCPA.GT.0 ) GOTO 100
c
c modified by XJQ: parallel on k-point
      deallocate(sumq_tmp)
      if(lparalleled) call mpi_multilevel_free(level)
c end-mod-xjq
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 ) THEN
C
         IF ( ITEST.EQ.5 ) THEN
            WRITE (77,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(PHASK(IE))
            WRITE (78,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKFE)
            WRITE (79,'(2f14.5,3x,2f14.5)') DREAL(ERYD),DIMAG(DETKGT)
         END IF
C
         IF ( IPRINT.GE.2 .AND. NCPA.GT.0 ) WRITE (12,99005) ERYD,
     &        ICPAFLAG,CPAERR,CPACORR,CPACHNG
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      DEALLOCATE (TAUK,MAUX,SUMQ,W1,MSSQLU,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'deallocate')
C
      IF ( MPI_KLOOP ) CALL DRV_MPI_BARRIER
C
C ----------------------------------------------------------------------
99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
99005 FORMAT ('E ',2F10.5,' CPA ',I5,3E12.5)
      END
