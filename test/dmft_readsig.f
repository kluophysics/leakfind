C*==dmft_readsig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_READSIG(NEGF,EGF,IPRINT)
C
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_CALCMODE,ONLY:IREL,LDAU
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,LOPT,LTXT_T,ITBOT,ITTOP,NT
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0,RY_EV
      USE MOD_MPI,ONLY:MPI_ID
      USE MOD_DMFT_LDAU,ONLY:KSELF,DMFTSIGMA,EREFLDAU,EREFDMFT,DMFTTEMP,
     &    DMFTSOLVER
c modified by XJQ: parallel on types
      USE MOD_MPI,ONLY:NPROCS,MPI_ID
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq
      IMPLICIT NONE
C*--DMFT_READSIG15
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NEGF
      COMPLEX*16 EGF(NEGF)
C
C Local variables
C
      COMPLEX*16 CWORK1(:,:),CWORK2(:,:)
      CHARACTER*80 FILNAM,STRINP
      INTEGER IE,IT,ITBOTSAV,ITTOPSAV,LFILNAM,LL,LMP
      LOGICAL READSIG
c modified by XJQ: parallel on types using my own subroutines
      character*40 routine
      parameter ( routine='dmft_readsig' )
      logical lparalleled, lopen
      integer it0_distr, it1_distr, nt_distr,
     &        it0_collct, it1_collct, nt_collct, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      real*8 :: t0_readsig, t1_readsig
      integer :: nt_dmft, itdmft_it(itbot:ittop), it_itdmft(1:ntmax)
      integer, dimension(:), allocatable :: disp_proc, nt_rev_proc
      complex*16 :: eref_tmp(1:ntmax)
      complex*16, dimension(:,:,:), allocatable :: dmftsig_tmp
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CWORK1,CWORK2
C
      READSIG = .TRUE.
C
      IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         ITBOT = 1
         ITTOP = NT
      END IF
      WRITE (6,'(/)')
c modified by XJQ: parallel on types
c parallel on types where there are electron correlations
      itdmft_it(itbot:ittop)=0
      it_itdmft(:)=0
      nt_dmft=0
      do it=itbot,ittop
        if(kself(it)==1) then
          nt_dmft = nt_dmft + 1
          itdmft_it(it) = nt_dmft
          it_itdmft(nt_dmft) = it
        endif
      enddo
c
      call get_comm_level(1,'parent',lparalleled,
     &                    parent_comm,nprocs_parent,parent_rank)
      if(nprocs_parent>1) call mpi_barrier(parent_comm,mpierr)
      call cpu_time(t0_readsig)
c
      call mpi_multilevel_distribute(routine,1,nt_dmft,
     &  lparalleled,it0_distr,it1_distr)
c end-mod-xjq
      DO IT = ITBOT,ITTOP
C
c modified by XJQ: parallel on types
         IF ( (KSELF(IT).EQ.1 .OR. LDAU .AND. LOPT(IT).GE.0) .and.
     &        itdmft_it(it)>=it0_distr .and. 
     &        itdmft_it(it)<=it1_distr ) THEN
c end-mod-xjq
            FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'_SIGUNFORM'//'.sig'
            LFILNAM = LTXT_T(IT) + 10 + 4
            INQUIRE (FILE=FILNAM(1:LFILNAM),EXIST=READSIG)
            DMFTSIGMA(1:NKM,1:NKM,IT,1:NEGF) = C0
C
            DMFTSIGMA(1:NKM,1:NKM,IT,1:NEGF) = C0
C
            IF ( READSIG ) THEN
C
               CALL DMFT_CLOSE_IOTMP(IOTMP)
C
               OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),STATUS='UNKNOWN',
     &               FORM='unformatted')
               REWIND IOTMP
C=======================================================
C                   *    *    *
C=======================================================
C     Read DMFT-FLEX self-energy from unformatted file:
               IF ( DMFTSOLVER(1:4).EQ.'FLEX' ) THEN
                  CALL DMFT_FLEXREADSIG(IOTMP,DMFTSIGMA,EGF,NKM,NEGF,IT,
     &                                  IREL,DMFTTEMP,LOPT(IT),NKMMAX,
     &                                  NEGF,NTMAX,EREFDMFT(IT))
C     Read DMFT-TMA self-energy from unformatted file:
               ELSE IF ( DMFTSOLVER(1:3).EQ.'TMA' ) THEN
                  CALL DMFT_SPTTMAREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,
     &               LOPT(IT),NKM,NKMMAX,NEGF,NTMAX,EREFDMFT(IT))
               ELSE IF ( LDAU ) THEN
                  CALL DMFT_LDAUREADSIG(IOTMP,DMFTSIGMA,IT,NEGF,LOPT(IT)
     &                                  ,NKM,NKMMAX,NTMAX,EREFLDAU(IT))
               ELSE IF ( DMFTSOLVER(1:3).EQ.'3BS' ) THEN
                  CALL DMFT_3BSREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,
     &                                 LOPT(IT),NKM,NKMMAX,NEGF,NTMAX,
     &                                 EFERMI)
               ELSE IF ( DMFTSOLVER(1:2).EQ.'ED' ) THEN
                  CALL DMFT_EDREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,
     &                                LOPT(IT),NKM,NKMMAX,NEGF,NTMAX,
     &                                EFERMI)
               ELSE
                  STOP '<DMFT_READSIG> Not yet implemented'
               END IF           ! FLEX,TMA,etc..
C
               WRITE (6,*) '         SELF ENERGY READ IN FROM FILE: ',
     &                     FILNAM(1:LFILNAM)
C
C
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
               CALL DMFT_CLOSE_IOTMP(IOTMP)
            ELSE
               WRITE (6,'(1X,79(''*''),/)')
               WRITE (6,*) 'SELF ENERGY READ IN FROM FILE: ',
     &                     FILNAM(1:LFILNAM)
               WRITE (6,*)
               WRITE (6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE (6,*) 'Warning: DMFT self-energy file not found'
               WRITE (6,*) '         self-energy set to zero        '
               WRITE (6,*) '                                        '
               WRITE (6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE (6,'(1X,79(''*''),/)')
C
               DMFTSIGMA(1:NKM,1:NKM,IT,1:NEGF) = C0
            END IF
C
         END IF                 ! KSELF
      END DO                    ! IT
c modified by XJQ: parallel on types
      if(lparalleled) then
        call mpi_barrier(parent_comm,mpierr)
        call get_comm_level(1,'inter ',lparalleled,
     &                      inter_comm,nprocs_inter,inter_rank)
c colloct it0_collct, it1_collct and nt_collct for further collecting
        allocate(disp_proc(0:nprocs_inter-1))
        allocate(nt_rev_proc(0:nprocs_inter-1))
        it0_collct = it_itdmft(it0_distr)
        it1_collct = it_itdmft(it1_distr)
        nt_collct = it1_collct - it0_collct + 1
        call collect_varstartnvar(1,inter_comm,nprocs_inter,
     &                            it0_collct,nt_collct,
     &                            disp_proc(0:),nt_rev_proc(0:))
c
        call mpi_barrier(parent_comm,mpierr)
        allocate(dmftsig_tmp(nkmmax,nkmmax,itbot:ittop))
        dmftsig_tmp(:,:,itbot:ittop)=c0
        call mpi_type_contiguous(nkmmax*nkmmax,mpi_double_complex,
     &                           mpi_row,mpierr)
        call mpi_type_commit(mpi_row,mpierr)
        do ie=1,negf
          call mpi_gatherv(dmftsigma(1,1,it0_collct,ie),nt_collct,
     &                     mpi_row,dmftsig_tmp(1,1,1),
     &                     nt_rev_proc(0:),disp_proc(0:),
     &                     mpi_row,0,inter_comm,mpierr)
          call mpi_barrier(parent_comm,mpierr)
          if(parent_rank==0)
     &      dmftsigma(:,:,itbot:ittop,ie) = dmftsig_tmp(:,:,itbot:ittop)
        enddo
c
        eref_tmp(:)=c0
        call mpi_gatherv(erefdmft(it0_collct),nt_collct,
     &                   mpi_double_complex,eref_tmp(1),
     &                   nt_rev_proc(0:),disp_proc(0:),
     &                   mpi_double_complex,0,inter_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        if(parent_rank==0) erefdmft(:) = eref_tmp(:)
c
        call mpi_barrier(parent_comm,mpierr)
        deallocate(dmftsig_tmp)
        deallocate(disp_proc)
        deallocate(nt_rev_proc)
        call mpi_bcast(erefdmft(1),ntmax,mpi_double_complex,
     &                 0,parent_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
        call mpi_bcast(dmftsigma(1,1,1,1),nkmmax*nkmmax*ntmax*negf
     &    ,mpi_double_complex,0,parent_comm,mpierr)
c
        call mpi_barrier(parent_comm,mpierr)
        call mpi_multilevel_free(1)
      endif
c end-mod-xjq
C
C=======================================================
C                   *    *    *
C=======================================================
C Print out the self-energy read:
C=======================================================
      IF ( IPRINT.GE.1 .AND. MPI_ID.EQ.0 ) THEN
c
         call system("mkdir diag-sig_read")
C
         CALL DMFT_CLOSE_IOTMP(IOTMP)
C
         DO IT = ITBOT,ITTOP
            STRINP = 'ImSig_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &               (1:LTXT_T(IT))//'_1.dat'
C
            LL = 6 + LDATSET + 1 + LTXT_T(IT) + 6
            OPEN (UNIT=IOTMP,FILE='diag-sig_read/'//STRINP(1:LL))
            DO IE = 1,NEGF
               WRITE (IOTMP,'(40f20.10)') DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                (DIMAG(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),LMP=5,9)
            END DO
            CLOSE (IOTMP)
         END DO
C
         CALL DMFT_CLOSE_IOTMP(IOTMP)
C
         DO IT = ITBOT,ITTOP
            STRINP = 'ImSig_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &               (1:LTXT_T(IT))//'_2.dat'
C
            LL = 6 + LDATSET + 1 + LTXT_T(IT) + 6
            OPEN (UNIT=IOTMP,FILE='diag-sig_read/'//STRINP(1:LL))
            DO IE = 1,NEGF
               WRITE (IOTMP,'(40f20.10)') DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                (DIMAG(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),LMP=14,18)
            END DO
            CLOSE (IOTMP)
C
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         DO IT = ITBOT,ITTOP
            STRINP = 'ReSig_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &               (1:LTXT_T(IT))//'_1.dat'
C
            LL = 6 + LDATSET + 1 + LTXT_T(IT) + 6
            OPEN (UNIT=IOTMP,FILE='diag-sig_read/'//STRINP(1:LL))
            DO IE = 1,NEGF
               WRITE (IOTMP,'(40f20.10)') DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                (DREAL(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),LMP=5,9)
            END DO
            CLOSE (IOTMP)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         DO IT = ITBOT,ITTOP
            STRINP = 'ReSig_'//DATSET(1:LDATSET)//'_'//TXT_T(IT)
     &               (1:LTXT_T(IT))//'_2.dat'
C
            LL = 6 + LDATSET + 1 + LTXT_T(IT) + 6
            OPEN (UNIT=IOTMP,FILE='diag-sig_read/'//STRINP(1:LL))
            DO IE = 1,NEGF
               WRITE (IOTMP,'(40f20.10)') DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                (DREAL(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),LMP=14,18)
            END DO
            CLOSE (IOTMP)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C=======================================================
C                   *    *    *
C=======================================================
C     TRANSFORM SELF-ENERGY into proper representation
C=======================================================
C
C     2012: all matrices assumed to have dimension NKMMAX
C
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).EQ.1 .OR. LDAU ) THEN
            IF ( IREL.EQ.3 ) THEN
               IF ( .NOT.ALLOCATED(CWORK1) )
     &              ALLOCATE (CWORK1(NKMMAX,NKMMAX),
     &              CWORK2(NKMMAX,NKMMAX))
               CWORK1 = C0
               CWORK2 = C0
               DO IE = 1,NEGF
C     Self-energy from the solver in the l,ml,ms representation
C
                  CWORK1(1:NKMMAX,1:NKMMAX)
     &               = DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
                  IF ( IPRINT.GT.0 .AND. MPI_ID.EQ.0 ) THEN
                     WRITE (6,*) '********************** IE=',IE
                     CALL CMATSTRUCT('DMFTSIGMA_CLM',CWORK1(1,1),NKM,
     &                               NKMMAX,IREL,IREL,0,1D-8,6)
                  END IF
C
                  IF ( DMFTSOLVER(1:3).EQ.'3BS' ) THEN
C                     CALL CHANGEREP(NKM,NKMMAX,CWORK1,'RLM>REL',CWORK2)
                     CALL CHANGEREP(NKM,NKMMAX,CWORK1,'CLM>REL',CWORK2)
                  ELSE IF ( DMFTSOLVER(1:3).EQ.'ED ' ) THEN
                     CALL CHANGEREP(NKM,NKMMAX,CWORK1,'RLM>REL',CWORK2)
                  ELSE
                     CALL CHANGEREP(NKM,NKMMAX,CWORK1,'CLM>REL',CWORK2)
                  END IF
C
                  DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
     &               = CWORK2(1:NKMMAX,1:NKMMAX)
                  IF ( IPRINT.GT.0 .AND. MPI_ID.EQ.0 ) THEN
                     WRITE (6,*) '********************** IE=',IE
                     CALL CMATSTRUCT('DMFTSIGMA_REL',
     &                               DMFTSIGMA(1,1,IT,IE),NKM,NKMMAX,
     &                               IREL,IREL,0,1D-8,6)
                  END IF
               END DO
               IF ( IPRINT.GT.1 .AND. MPI_ID.EQ.0 ) THEN
                  CALL DMFT_CLOSE_IOTMP(IOTMP)
                  OPEN (UNIT=IOTMP,FILE='sigma_dos1_imlam.dat')
                  CLOSE (IOTMP,STATUS='DELETE')
                  OPEN (UNIT=IOTMP,FILE='sigma_dos1_imlam.dat')
                  DO IE = 1,NEGF
                     WRITE (IOTMP,'(40f20.10)')
     &                      DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                      (DIMAG(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),
     &                      LMP=9,18)
                  END DO
                  WRITE (IOTMP,*) 
     &                     '###########################################'
                  WRITE (IOTMP,*) 
     &                     '###########################################'
                  CALL DMFT_CLOSE_IOTMP(IOTMP)
                  OPEN (UNIT=IOTMP,FILE='sigma_dos1_relam.dat')
                  CLOSE (IOTMP,STATUS='DELETE')
                  OPEN (UNIT=IOTMP,FILE='sigma_dos1_relam.dat')
                  DO IE = 1,NEGF
                     WRITE (IOTMP,'(40f20.10)')
     &                      DREAL((EGF(IE)-EFERMI)*RY_EV),
     &                      (DREAL(DMFTSIGMA(LMP,LMP,IT,IE)*RY_EV),
     &                      LMP=9,18)
                  END DO
                  WRITE (IOTMP,*) 
     &                     '###########################################'
                  WRITE (IOTMP,*) 
     &                     '###########################################'
                  CALL DMFT_CLOSE_IOTMP(IOTMP)
               END IF
            ELSE
               IF ( .NOT.ALLOCATED(CWORK1) )
     &              ALLOCATE (CWORK1(NKMMAX,NKMMAX),
     &              CWORK2(NKMMAX,NKMMAX))
               CWORK1 = C0
               CWORK2 = C0
               DO IE = 1,NEGF
C     Self-energy from the solver in the CLM representation
C
                  CWORK1(1:NKMMAX,1:NKMMAX)
     &               = DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
C
                  IF ( DMFTSOLVER(1:3).EQ.'ED ' ) THEN
                     CWORK2 = CWORK1
                  ELSE
                     CALL CHANGEREP(NKM,NKMMAX,CWORK1,'CLM>RLM',CWORK2)
                  END IF
C
                  DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
     &               = CWORK2(1:NKMMAX,1:NKMMAX)
C
                  IF ( IPRINT.GT.0 .AND. MPI_ID.EQ.0 ) THEN
                     WRITE (6,*) '********************** IE=',IE
                     CALL CMATSTRUCT('DMFTSIGMA_CLM',
     &                               DMFTSIGMA(1,1,IT,IE),NKM,NKMMAX,
     &                               IREL,IREL,0,1D-8,6)
                  END IF
               END DO
            END IF
C
         END IF                 !KSELF
      END DO                    !IT
C=======================================================
C                   *    *    *
C=======================================================
      IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
      END IF
      END
C
