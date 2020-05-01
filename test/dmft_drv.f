C*==dmft_drv.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_DRV(ITRSCF,GFMAT,DMFTSIGMA,WEGF,EGFOLD,EGFNEW,
     &                    DMFTMIX,IPRINT,EREFDMFT,OBS_LT,SCLNOS)
C
C   ********************************************************************
C   *                                                                  *
C   *   driver to run the various DMFT-solvers                         *
C   *   according to the calculation mode                              *
C   *                                                                  *
C   *   DMFT solvers:                                                  *
C   *                DMFT_SOSPTFLEX: L. Pourovskii SOC+FLEX            *
C   *                DMFT-TMA: S. Chadov  TMA on real axis (all modes) *
C   *                                                                  *
C   *   INPUT: GFMAT: Greens function matrix in l,ml,ms represantation *
C   *                 in:  REL -complex spherical harmonics            *
C   *                 in:  SREL-real spherical harmonics               *
C   *   OUTPUT: DMFTSIGMA: Self-energy on complex contour in the       *
C   *           representation of IREL (IREL=2 RLM, IREL=3 REL)        *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NOBSMAX,NLMAX,NKMMAX
      USE MOD_ENERGY,ONLY:EFERMI,NEMAX,NETAB,ETAB
      USE MOD_TYPES,ONLY:NTMAX,LTXT_T,TXT_T,NT,LOPT,ITBOT,ITTOP
      USE MOD_FILES,ONLY:IFILGFWF,FOUND_SECTION
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_DMFT_LDAU,ONLY:DMFTTEMP,DMFTSCF,KSELF,DMFTSOLVER,JEFF,
     &    UEFF,DMFTDBLC,SIGTOL
      USE MOD_CONSTANTS,ONLY:C0,PI
c modified by XJQ: parallel on types using my own subroutines
c                  confict with mpi_types=.true.
      use mod_energy,only:igrid,nefd1,nefd2,nefd3,nepol,necontn,
     &                    lactive_contn
      USE MOD_MPI,ONLY:NPROCS,MPI_ID
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq 
      IMPLICIT NONE
C*--DMFT_DRV34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DMFT_DRV')
C
C Dummy arguments
C
      REAL*8 DMFTMIX,SCLNOS
      INTEGER IPRINT,ITRSCF
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),EGFNEW(NEMAX),
     &           EGFOLD(NEMAX),EREFDMFT(NTMAX),
     &           GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX),WEGF(NEMAX)
      REAL*8 OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CEMIG(:),CSIGSAV(:,:,:,:),CWKT(:),CWORK(:,:,:),
     &           CWORK1(:,:),CWORK2(:,:),CWORK3(:,:),
     &           DMFTSIGMAOM(:,:,:,:),EOM(:),EREFLDAU(:),ERYD,
     &           GFMAT0(:,:,:,:),GFMATOM(:,:,:,:),LDAUSIGMA(:,:,:,:),
     &           SIGLDAU(:,:,:),WETABOM(:)
      CHARACTER*4 DBLCSAV
      LOGICAL DMFT_GALMIG,LDAMP,LSSITE,MPI_TYPES,NOPADE,TMA_NOBATH,
     &        TMA_NOCCSCL,TMA_SPINFLIPOFF,TMA_STATIC,UDT
c original code define rwkt as integer, it is probably wrong
      REAL*8 ED_EDC,ED_EMAX,ED_EMIN,EILOW,ELDAU(:),EMAX,EMIN,FDAMP,
     &       KDAMP,KPENALTY,NBATH_PORB,RTMP,SIGERR,SIGERRT(:),TIME,
     &       TIME0,TIME_ALL,TIME_ALL0,TMA_REFACTOR,UDMFT(:,:,:,:,:),
     &       RWKT(:)
      INTEGER I,I1,I2,I3,I4,IE,IFLEX,INC,IODMFT,IPRINTDMFT,IPROC,
     &        IPROCE(:),IPROCOM(:),IPROCT(:),IT,ITR,ITRDMFT,IVEC,J,JE,
     &        LWKKMTE,NBATH,NELEC,NENVEXTRA,NITERDMFT,NLM,NOM,NREP,
     &        NSWEEP,NTAU
      CHARACTER*80 STR
      CHARACTER*3 TMA_DBLC
      SAVE CSIGSAV
c modified by XJQ: parallel on types using my own subroutines
      logical lparalleled, lopen
      integer it0_distr, it1_distr, nt_distr, 
     &        it0_collct, it1_collct, nt_collct, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      real*8 :: t0_dmft, t1_dmft
      integer :: nt_dmft, itdmft_it(itbot:ittop), it_itdmft(1:ntmax)
      integer, dimension(:), allocatable :: disp_proc, nt_rev_proc
      complex*16, dimension(:,:,:), allocatable :: dmftsig_tmp
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CWORK,CWORK1,UDMFT,CEMIG,CWORK2,CWORK3,CSIGSAV
      ALLOCATABLE IPROCE,CWKT,IPROCT,LDAUSIGMA,EREFLDAU,ELDAU
      ALLOCATABLE SIGLDAU,GFMAT0,GFMATOM,EOM,WETABOM,DMFTSIGMAOM
      ALLOCATABLE IPROCOM,RWKT,SIGERRT
C
C
      CALL CPU_TIME(TIME_ALL0)
C
      ALLOCATE (UDMFT(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX,NTMAX))
      ALLOCATE (CWORK(NKMMAX,NKMMAX,NEMAX))
      ALLOCATE (CEMIG(NT),IPROCE(NEMAX),IPROCT(NTMAX))
      ALLOCATE (GFMAT0(NKMMAX,NKMMAX,NEMAX,NTMAX))
C
      ALLOCATE (SIGERRT(NTMAX))
      ALLOCATE (SIGLDAU(NKMMAX,NKMMAX,NTMAX))
      SIGLDAU = C0
      CALL CINIT(NKMMAX*NKMMAX*NEMAX,CWORK)
      CALL CINIT(NT,CEMIG)
      WRITE (6,99001) TRIM(DMFTSOLVER)
C
      SIGERRT = 0.0D0
      SIGERR = 0.0D0
C     Prepared for lloyd
C      IF ( LLOYD ) GFMAT = GFMAT*SCLNOS
      IF ( .FALSE. ) GFMAT = GFMAT*SCLNOS
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( DMFTSCF ) THEN
C     MPI  over energy
         IPROCE(:) = 0
         IPROCT(:) = 0
C
         IF ( NPROCS.GT.1 ) THEN
            IPROC = 0
            INC = 1
            DO JE = 1,NETAB(1) - 1
               IE = JE
               IPROC = IPROC + INC
               IF ( IPROC.EQ.NPROCS ) THEN
                  IPROC = NPROCS - 1
                  INC = -1
               ELSE IF ( IPROC.EQ.-1 ) THEN
                  IPROC = 0
                  INC = 1
               END IF
               IPROCE(IE) = IPROC
            END DO
         END IF
      ELSE
C     MPI over types
         IPROCT(1:NTMAX) = 0
C
         IF ( NPROCS.GT.1 ) THEN
            IF ( NPROCS.GT.(ITTOP-ITBOT+1) ) IPROCT(1:NTMAX) = -1
            IPROCT(ITTOP) = 0
C
            IPROC = 0
            INC = 1
            DO IT = ITBOT,ITTOP - 1
               IPROC = IPROC + INC
               IF ( IPROC.EQ.NPROCS ) THEN
                  IPROC = NPROCS - 1
                  INC = -1
               ELSE IF ( IPROC.EQ.-1 ) THEN
                  IPROC = 0
                  INC = 1
               END IF
               IPROCT(IT) = IPROC
            END DO
            MPI_TYPES = .TRUE.
            IF ( MPI_TYPES .AND. IPRINT.GT.0 ) THEN
               WRITE (6,*) 'DMFT: MPI Parallelisation over TYPES'
               DO IT = ITBOT,ITTOP
                  WRITE (6,*) 'IT=',IT,'KSELF=',KSELF(IT),'IPROC=',
     &                        IPROCT(IT)
               END DO
            END IF
C
C           Paralelisation over types is prepared but not yet tested.
            IPROCT = 0
            MPI_TYPES = .FALSE.
C FLEX-IDM Solver has its own parallelisation
c pade_averaging has mpi, change it temporarily
c            IF ( DMFTSOLVER(1:8).EQ.'FLEX-IDM' ) THEN
            IF ( DMFTSOLVER(1:4).EQ.'FLEX' ) THEN
               IPROCT = MPI_ID
               MPI_TYPES = .FALSE.
            END IF
         END IF
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C
C================================================================
C  Consistency checks:
C               FLEX     SOFLEX      TMA
C   ASA IREL=2   X         --        X
C       IREL=3   --        --        X
C
C   FP  IREL=2   X         --        X
C       IREL=3   --        X         X
C
C     FLEX: Is using its own calculation of U(1,2,3,4)
C     SO-FLEX: Is using DMFT_U_MATRIX for U(1,2,3,4)
C     TMA: In FP case using DMFT_U_MATRIX and in ASA own
C
      IF ( IREL.LE.1 ) STOP 'DMFT only for IREL GE 2'
C=======================================================
      IFLEX = 1
      TMA_REFACTOR = 1D0
      IPRINTDMFT = 0
      TMA_STATIC = .FALSE.
      TMA_SPINFLIPOFF = .FALSE.
      TMA_NOBATH = .FALSE.
      DMFT_GALMIG = .FALSE.
      TMA_NOCCSCL = .FALSE.
      LDAMP = .FALSE.
      KDAMP = 0.25D0
      LSSITE = .FALSE.
      NITERDMFT = 1
      IODMFT = 387
      NBATH_PORB = 2.0D0
      NBATH = 2
      KPENALTY = 0.0D0
      NELEC = 5
      NENVEXTRA = 0
      NOPADE = .FALSE.
      MPI_TYPES = .FALSE.
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
C
         IF ( ORBPOL(1:9).EQ.'DMFT-FLEX' ) THEN
            CALL SECTION_SET_INTEGER('NMATSUB',NOM,2048,0)
            UDT = .FALSE.
            CALL SECTION_FIND_KEYWORD('DOTMA',UDT)
            IF ( UDT ) IFLEX = 0
         END IF
C
         IF ( DMFTSCF ) THEN
            WRITE (6,99006)
            CALL SECTION_SET_INTEGER('NITER',NITERDMFT,10,0)
         END IF
C
C
         CALL SECTION_FIND_KEYWORD('DMFT_GALMIG',UDT)
         IF ( UDT ) DMFT_GALMIG = .TRUE.
C
         CALL SECTION_SET_INTEGER('IPRINTDMFT',IPRINTDMFT,IPRINT,0)
C
         CALL SECTION_SET_REAL('TMA_REFACTOR',TMA_REFACTOR,1D0,0)
C
         CALL SECTION_FIND_KEYWORD('TMA_STATIC',UDT)
         IF ( UDT ) TMA_STATIC = .TRUE.
C
         CALL SECTION_FIND_KEYWORD('TMA_SPINFLIPOFF',UDT)
         IF ( UDT ) TMA_SPINFLIPOFF = .TRUE.
C
C
         CALL SECTION_SET_STRING('TMA_DBLC',TMA_DBLC,'AMF',0)
C
         CALL SECTION_FIND_KEYWORD('TMA_NOBATH',UDT)
         IF ( UDT ) TMA_NOBATH = .TRUE.
C
         CALL SECTION_FIND_KEYWORD('TMA_NOCCSCL',UDT)
         IF ( UDT ) TMA_NOCCSCL = .TRUE.
         CALL SECTION_FIND_KEYWORD('DAMP',UDT)
         IF ( UDT ) LDAMP = .TRUE.
         CALL SECTION_SET_REAL('KDAMP',KDAMP,0.25D0,0)
C
         IF ( ORBPOL(1:8).EQ.'DMFT-QMC' ) THEN
            CALL SECTION_SET_INTEGER('NMATSUB',NOM,2048,0)
            WRITE (6,*) 'Number of matsubara NOM=',NOM
            CALL SECTION_SET_INTEGER('NSWEEP',NSWEEP,5000,0)
            WRITE (6,*) 'Number of qmc sweeps NSWEEP=',NSWEEP
            CALL SECTION_SET_INTEGER('NTAU',NTAU,60,0)
            WRITE (6,*) 'Number of qmc sweeps NTAU=',NTAU
         END IF
C
C
         IF ( ORBPOL(1:8).EQ.'DMFT-ED ' ) THEN
            CALL SECTION_SET_REAL('ED_EDC',ED_EDC,15.0D0,0)
            CALL SECTION_SET_REAL('ED_EMIN',ED_EMIN,DREAL(EGFOLD(1)),0)
            CALL SECTION_SET_REAL('ED_EMAX',ED_EMAX,
     &                            DREAL(EGFOLD(NETAB(1))),0)
         END IF
         IF ( ORBPOL(1:8).EQ.'DMFT-EDN' ) THEN
            CALL SECTION_SET_INTEGER('NMATSUB',NOM,2048,0)
            CALL SECTION_SET_REAL('NBATH_PORB',NBATH_PORB,2.0D0,0)
            CALL SECTION_SET_INTEGER('NBATH',NBATH,2,0)
            CALL SECTION_SET_REAL('KPENALTY',KPENALTY,0.0D0,0)
            CALL SECTION_SET_INTEGER('NELEC',NELEC,7,0)
            CALL SECTION_SET_INTEGER('NENVEXTRA',NENVEXTRA,0,0)
            CALL SECTION_FIND_KEYWORD('NOPADE',UDT)
            IF ( UDT ) NOPADE = .TRUE.
C
         END IF
C
      END IF
C
      IF ( NPROCS>1 ) CALL DRV_MPI_BARRIER
C
      DO ITRDMFT = 1,NITERDMFT
C
         IF ( DMFTSCF ) THEN
            WRITE (6,*) 'CALCULATING Green''s function matrix'
            GFMAT = C0
            DO IE = 1,NETAB(1)
C
C     MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
               IF ( MPI_ID.EQ.IPROCE(IE) ) THEN
C     MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
                  ERYD = ETAB(IE,1)
                  IF ( DIMAG(ERYD).GT.4.0D0 ) LSSITE = .TRUE.
                  WRITE (6,*) 'E=',IE,ERYD,LSSITE
C
                  CALL DMFT_TAU_DRIVE(ERYD,IE,LSSITE,GFMAT,DMFTSIGMA,
     &                                NEMAX)
               END IF
            END DO
            IF ( NPROCS>1 ) THEN
               IF ( ALLOCATED(CWKT) ) DEALLOCATE (CWKT)
               ALLOCATE (CWKT(NKMMAX))
               LWKKMTE = NKMMAX
               DO IVEC = 1,NKMMAX
                  DO IE = 1,NEMAX
                     DO IT = ITBOT,ITTOP
                        CALL DRV_MPI_REDUCE_C(GFMAT(1,IVEC,IE,IT),
     &                     CWKT(1),LWKKMTE)
                     END DO
                  END DO
               END DO
               DEALLOCATE (CWKT)
C
               CALL DRV_MPI_BARRIER
C
               CALL DRV_MPI_BCAST_C(0,GFMAT(1,1,1,1),
     &                              NKMMAX*NKMMAX*NTMAX*NETAB(1))
C
            END IF
C
C
         END IF
         IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
C
            ALLOCATE (LDAUSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX))
            LDAUSIGMA(:,:,:,:) = C0
            ALLOCATE (EREFLDAU(NTMAX),ELDAU(NTMAX))
            ELDAU = 0.0D0
            EREFLDAU = C0
            DBLCSAV = DMFTDBLC
            DMFTDBLC(1:4) = '    '
            DMFTDBLC(1:3) = TMA_DBLC(1:3)
C
            GFMAT0 = C0
C
            DO IT = ITBOT,ITTOP
               IF ( KSELF(IT).EQ.1 )
     &              CALL DMFT_CALCG0(GFMAT,DMFTSIGMA,NKM,NEMAX,LOPT(IT),
     &              NKMMAX,GFMAT0,NTMAX,IT)
            END DO
C
C
C
            GFMAT0 = GFMAT
            CALL LDAU_DRV(GFMAT0,LDAUSIGMA,WEGF,EREFLDAU,ELDAU,OBS_LT,
     &                    1.0D0,.FALSE.)
C
            DMFTDBLC = DBLCSAV
            IF ( .NOT.ALLOCATED(CWORK1) )
     &           ALLOCATE (CWORK1(NKMMAX,NKMMAX))
            IF ( .NOT.ALLOCATED(CWORK2) )
     &           ALLOCATE (CWORK2(NKMMAX,NKMMAX))
            CWORK1 = C0
            CWORK2 = C0
C
            IF ( IREL.EQ.3 ) THEN
               STR(1:7) = 'REL>CLM'
            ELSE
               STR(1:7) = 'RLM>CLM'
            END IF
C
            DO IT = 1,NTMAX
               CWORK1(1:NKMMAX,1:NKMMAX)
     &            = LDAUSIGMA(1:NKMMAX,1:NKMMAX,IT,1)
C
               CALL CHANGEREP(NKM,NKMMAX,CWORK1(1,1),STR(1:7),CWORK2)
C
               SIGLDAU(1:NKMMAX,1:NKMMAX,IT) = CWORK2(1:NKMMAX,1:NKMMAX)
            END DO
            DEALLOCATE (CWORK1,CWORK2)
         END IF
C
         IF ( NOPADE ) THEN
            WRITE (6,*) 'CALCULATING Greens function matrix: MATSUBARA'
            ALLOCATE (GFMATOM(NKMMAX,NKMMAX,NOM,NTMAX))
            ALLOCATE (EOM(NOM),WETABOM(NOM))
            ALLOCATE (DMFTSIGMAOM(NKMMAX,NKMMAX,NTMAX,NOM))
            ALLOCATE (IPROCOM(NOM))
            DO IE = 1,NOM
               IPROCOM(IE) = 0
            END DO
            IF ( NPROCS.GT.1 ) THEN
               IPROC = 0
               INC = 1
               DO JE = 1,NOM - 1
                  IE = JE
                  IPROC = IPROC + INC
                  IF ( IPROC.EQ.NPROCS ) THEN
                     IPROC = NPROCS - 1
                     INC = -1
                  ELSE IF ( IPROC.EQ.-1 ) THEN
                     IPROC = 0
                     INC = 1
                  END IF
                  IPROCOM(IE) = IPROC
               END DO
            END IF
C
            DMFTSIGMAOM = C0
            EMIN = DREAL(EGFOLD(NETAB(1)))
            EMAX = DREAL(EGFOLD(NETAB(1)))
            CALL EPATH(9,EMIN,EMAX,DMFTTEMP,NOM,EOM,WETABOM,EILOW,
     &                 IPRINT,NOM)
            GFMATOM = C0
            DO IE = 1,NOM
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
               IF ( MPI_ID.EQ.IPROCOM(IE) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                  ERYD = EOM(IE)
                  IF ( DIMAG(ERYD).GT.-2.0D0 ) LSSITE = .TRUE.
                  WRITE (6,*) 'E=',IE,ERYD,LSSITE
                  CALL DMFT_TAU_DRIVE(ERYD,IE,LSSITE,GFMATOM,
     &                                DMFTSIGMAOM,NOM)
               END IF
            END DO
            IF ( NPROCS>1 ) THEN
               IF ( ALLOCATED(CWKT) ) DEALLOCATE (CWKT)
               ALLOCATE (CWKT(NKMMAX))
               LWKKMTE = NKMMAX
               DO IVEC = 1,NKMMAX
                  DO IE = 1,NOM
                     DO IT = ITBOT,ITTOP
                        CALL DRV_MPI_REDUCE_C(GFMATOM(1,IVEC,IE,IT),
     &                     CWKT(1),LWKKMTE)
                     END DO
                  END DO
               END DO
               DEALLOCATE (CWKT)
C
               CALL DRV_MPI_BARRIER
C
               CALL DRV_MPI_BCAST_C(0,GFMATOM(1,1,1,1),
     &                              NKMMAX*NKMMAX*NTMAX*NOM)
C
            END IF
         END IF
C
         IF ( DMFTSCF ) THEN
            ITR = ITRDMFT
         ELSE
            ITR = ITRSCF
         END IF
C
C=======================================================
         UDMFT(1:2*NLMAX,1:2*NLMAX,1:2*NLMAX,1:2*NLMAX,1:NTMAX) = 0.0D0
C
         IF ( IREL.GE.2 ) THEN
            DO IT = ITBOT,ITTOP
               IF ( KSELF(IT).GT.0 ) THEN
                  I1 = ITBOT
                  I2 = ITTOP
                  ITBOT = IT
                  ITTOP = IT
                  CALL DMFT_U_MATRIX(EREFDMFT,IPRINTDMFT,JEFF,IFILGFWF,
     &                               UEFF,UDMFT)
                  ITBOT = I1
                  ITTOP = I2
C
               END IF
            END DO
         END IF
C=======================================================
C                   *    *    *
C=======================================================
C
C=======================================================
C Call the DMFT-solver:
C=======================================================
C
         IF ( .NOT.ALLOCATED(CSIGSAV) ) THEN
            ALLOCATE (CSIGSAV(NKMMAX,NKMMAX,NTMAX,NETAB(1)))
            IF ( .NOT.ALLOCATED(CWORK1) )
     &           ALLOCATE (CWORK1(NKMMAX,NKMMAX))
            IF ( IREL.EQ.3 ) THEN
               STR(1:7) = 'REL>CLM'
            ELSE
               STR(1:7) = 'RLM>CLM'
            END IF
            DO IT = 1,NTMAX
               DO IE = 1,NETAB(1)
                  CWORK(1:NKMMAX,1:NKMMAX,IE)
     &               = DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
C
                  CALL CHANGEREP(NKM,NKMMAX,CWORK(1,1,IE),STR(1:7),
     &                           CWORK1)
C
                  CSIGSAV(1:NKMMAX,1:NKMMAX,IT,IE)
     &               = CWORK1(1:NKMMAX,1:NKMMAX)
               END DO
            END DO
         END IF
c modified by XJQ: parallel on types using my own subroutines
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
     &                       parent_comm,nprocs_parent,parent_rank)
         if(nprocs_parent>1) call mpi_barrier(parent_comm,mpierr)
         if(mpi_id==0) call system("mkdir diag-sig")
         if(mpi_id==0) call system("mkdir diag-gf")
         call cpu_time(t0_dmft)
c
         if(mpi_types) lothermpi=.true. ! mpi_types is used by Jan but not XJQ
         call mpi_multilevel_distribute(routine,1,nt_dmft,
     &     lparalleled,it0_distr,it1_distr)
c end-mod-xjq
         DMFTSIGMA(:,:,:,:) = C0
         DO IT = ITBOT,ITTOP
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
c modified by XJQ: parallel on types using my own subroutines
c if mpi_types, use Jan's parallel, otherwise, use XJQ's
            IF ( KSELF(IT).EQ.1 ) THEN
               IF ( MPI_ID.EQ.IPROCT(IT) .and.
     &              itdmft_it(it)>=it0_distr .and. 
     &              itdmft_it(it)<=it1_distr) THEN
                  CALL CPU_TIME(TIME0)
                  write(*,*) 'it=',it,' itdmft_it=',itdmft_it(it),
     &                       ' it0_distr=',it0_distr,
     &                       ' it1_distr=',it1_distr
c end-mod-xjq
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
                  DO IE = 1,NETAB(1)
c                     DO I = 1,NKMMAX
c                        DO J = 1,NKMMAX
c                           CWORK(J,I,IE) = CSIGSAV(J,I,IT,IE)
                     cwork(:,:,ie)=csigsav(:,:,it,ie)
c                        END DO
c                     END DO
                  END DO
C
                  WRITE (6,99002) IT,LOPT(IT),UEFF(IT),JEFF(IT)
C
C Call the solver:
                  IF ( DMFTSOLVER(1:4).EQ.'FLEX' .AND. DMFTSOLVER(1:8)
     &                 .NE.'FLEX-IDM' ) THEN
C
                     CALL DMFT_DRV_SOSPTFLEX(GFMAT(1,1,1,IT),CWORK,NKM,
     &                  NETAB(1),LOPT(IT),UEFF(IT),JEFF(IT),
     &                  UDMFT(1,1,1,1,IT),DMFTTEMP,DMFTDBLC,NOM,IODMFT,
     &                  IPRINTDMFT,IFLEX,TXT_T(IT),LTXT_T(IT),ITR,
     &                  DMFTMIX,NKMMAX,NLMAX,WEGF,EGFOLD,EGFNEW,
     &                  EREFDMFT(IT),IT,NTMAX,TMA_STATIC,TMA_NOBATH,
     &                  TMA_DBLC,OBS_LT(0,1,1,IT),TMA_NOCCSCL)
C
                  ELSE IF ( DMFTSOLVER(1:3).EQ.'TMA' ) THEN
                     CALL DMFT_DRV_SPTTMA(ITR,NKMMAX,NKM,NETAB(1),
     &                  DMFTMIX,LOPT(IT),EGFOLD,EGFNEW,EREFDMFT(IT),
     &                  WEGF,GFMAT(1,1,1,IT),CWORK,TXT_T(IT),LTXT_T(IT),
     &                  EFERMI,IPRINTDMFT,UDMFT(1,1,1,1,IT),NLMAX,
     &                  DMFTDBLC,TMA_SPINFLIPOFF,TMA_REFACTOR,
     &                  TMA_STATIC,UEFF(IT),JEFF(IT),TMA_DBLC,
     &                  TMA_NOBATH,IT,NTMAX,OBS_LT(0,1,1,IT),
     &                  TMA_NOCCSCL,SIGLDAU(1,1,IT))
C
                  ELSE IF ( DMFTSOLVER(1:8).EQ.'FLEX-IDM' ) THEN
                     CALL DMFT_DRV_SOSPTFLEX_IDM(GFMAT(1,1,1,IT),CWORK,
     &                  NKM,NETAB(1),LOPT(IT),UEFF(IT),JEFF(IT),
     &                  UDMFT(1,1,1,1,IT),DMFTTEMP,DMFTDBLC,NOM,IODMFT,
     &                  IPRINTDMFT,IFLEX,TXT_T(IT),LTXT_T(IT),ITR,
     &                  DMFTMIX,NKMMAX,NLMAX,WEGF,EGFOLD,EGFNEW,
     &                  EREFDMFT(IT),IT,NTMAX,TMA_STATIC,TMA_NOBATH,
     &                  TMA_DBLC,OBS_LT(0,1,1,IT),TMA_NOCCSCL,
     &                  SIGLDAU(1,1,IT))
C
                  ELSE IF ( DMFTSOLVER(1:6).EQ.'ED-IDM' ) THEN
                     CALL STOP_MESSAGE(ROUTINE,'ED SOLVER NOT YET')
C                     CALL DMFT_DRV_ED_IDM(GFMAT(1,1,1,IT),CWORK,
C     &                    NKM,NETAB(1),LOPT(IT),UEFF(IT),JEFF(IT),
C     &                   UDMFT(1,1,1,1,IT),DMFTTEMP,DMFTDBLC,NOM,IODMFT,
C     &                    IPRINTDMFT,IFLEX,TXT_T(IT),LTXT_T(IT),ITR,
C     &                    DMFTMIX,NKMMAX,NLMAX,WEGF,EGFOLD,EGFNEW,
C     &                    EREFDMFT(IT),IT,NTMAX,TMA_STATIC,TMA_NOBATH,
C     &             TMA_DBLC,OBS_LT(0,1,1,IT),TMA_NOCCSCL,SIGLDAU(1,1,
C     &                    IT),NBATH_PORB,NBATH,KPENALTY,NELEC,
C     &                    NENVEXTRA,NOPADE,GFMATOM(1,1,1,IT),EOM)
C
                  ELSE IF ( DMFTSOLVER(1:2).EQ.'ED' ) THEN
                     CALL STOP_MESSAGE(ROUTINE,'ED SOLVER NOT YET')
                     IF ( .NOT.ALLOCATED(CWORK1) )
     &                    ALLOCATE (CWORK1(NKMMAX,NKMMAX))
                     IF ( .NOT.ALLOCATED(CWORK2) )
     &                    ALLOCATE (CWORK2(NKMMAX,NKMMAX))
                     DO IE = 1,NETAB(1)
                        CWORK1(:,:) = GFMAT(:,:,IE,IT)
                        CALL CHANGEREP(NKM,NKMMAX,CWORK1,'CLM>RLM',
     &                                 CWORK2)
                        GFMAT(:,:,IE,IT) = CWORK2(:,:)
                     END DO
C                     CALL DMFT_DRV_ED(GFMAT(1,1,1,IT),CWORK,EGFOLD,
C     &                    EGFNEW,WEGF,
C     &                    NETAB(1),NKMMAX,ED_EDC,LOPT(IT),NKM,DMFTMIX,
C     &                    ITR,IPRINTDMFT,tma_noccscl,TMA_NOBATH,
C     &                    nos_lt(1,IT),nlmax,IODMFT,TXT_T(IT),
C     &                    LTXT_T(IT))
                     DO IE = 1,NETAB(1)
                        CWORK1(:,:) = CWORK(:,:,IE)
                        CALL CHANGEREP(NKM,NKMMAX,CWORK1,'RLM>CLM',
     &                                 CWORK2)
                        CWORK(:,:,IE) = CWORK2(:,:)
                     END DO
C
                  ELSE IF ( ORBPOL(1:8).NE.'DMFT-QMC' ) THEN
                     WRITE (6,*) 'STOP:',ORBPOL(1:10)
                     CALL STOP_MESSAGE(ROUTINE,'SOLVER: NOT YET')
                  END IF
C
C     Calculate the Galitsky-Migdal contribution on the semicircle:
C
                  IF ( DMFT_GALMIG ) THEN
C
                     I1 = LOPT(IT)**2 + 1
                     I2 = I1 + (LOPT(IT)*2+1) - 1
                     I3 = NKM/2 + I1
                     I4 = I3 + (LOPT(IT)*2+1) - 1
C
                     NLM = LOPT(IT)*2 + 1
                     NREP = 2*NLM
C
                     IF ( .NOT.ALLOCATED(CWORK1) )
     &                    ALLOCATE (CWORK1(NKMMAX,NKMMAX))
C
                     IF ( .NOT.ALLOCATED(CWORK2) )
     &                    ALLOCATE (CWORK2(NREP,NREP))
                     IF ( .NOT.ALLOCATED(CWORK3) )
     &                    ALLOCATE (CWORK3(NREP,NREP))
C
C
                     DO IE = 1,NETAB(1)
C
                        CWORK1(1:NREP/2,1:NREP/2)
     &                     = GFMAT(I1:I2,I1:I2,IE,IT)
C
                        CWORK1(NREP/2+1:NREP,NREP/2+1:NREP)
     &                     = GFMAT(I3:I4,I3:I4,IE,IT)
C
                        CWORK1(1:NREP/2,NREP/2+1:NREP)
     &                     = GFMAT(I1:I2,I3:I4,IE,IT)
C
                        CWORK1(NREP/2+1:NREP,1:NREP/2)
     &                     = GFMAT(I3:I4,I3:I4,IE,IT)
C
                        CWORK2(1:NREP/2,1:NREP/2)
     &                     = CWORK(I1:I2,I1:I2,IE)
C
                        CWORK2(NREP/2+1:NREP,NREP/2+1:NREP)
     &                     = CWORK(I3:I4,I3:I4,IE)
C
                        CWORK2(1:NREP/2,NREP/2+1:NREP)
     &                     = CWORK(I1:I2,I3:I4,IE)
C
                        CWORK2(NREP/2+1:NREP,1:NREP/2)
     &                     = CWORK(I3:I4,I3:I4,IE)
C
                        CWORK3(1:NREP,1:NREP)
     &                     = MATMUL(CWORK1(1:NREP,1:NREP),
     &                     CWORK2(1:NREP,1:NREP))
C
                        DO I1 = 1,NREP
                           CEMIG(IT) = CEMIG(IT) + CWORK3(I1,I1)
     &                                 *WEGF(IE)
                        END DO
                     END DO
                     CEMIG(IT) = -(CEMIG(IT)/PI)
                  END IF
C
C=======================================================
C Change the representation from lms to kappa-mue
C=======================================================
C
                  IF ( .NOT.ALLOCATED(CWORK1) )
     &                 ALLOCATE (CWORK1(NKMMAX,NKMMAX))
C
                  DO IE = 1,NETAB(1)
                     if(lactive_contn .and. ie > nefd1+nefd2+nefd3
     &                 .and. ie <= netab(1)-nepol) cycle
                     CWORK1(1:NKMMAX,1:NKMMAX)
     &                  = CSIGSAV(1:NKMMAX,1:NKMMAX,IT,IE)
     &                  - CWORK(1:NKMMAX,1:NKMMAX,IE)
                     SIGERRT(IT) = MAX(SIGERRT(IT),MAXVAL(ABS(CWORK1)))
                  END DO
                  CSIGSAV(:,:,IT,1:NETAB(1)) = CWORK(:,:,1:NETAB(1))
                  IF ( IREL.EQ.3 ) THEN
C
                     IF ( .NOT.ALLOCATED(CWORK1) )
     &                    ALLOCATE (CWORK1(NKMMAX,NKMMAX))
C
                     CWORK1 = C0
                     DO IE = 1,NETAB(1)
C
                        CALL CHANGEREP(NKM,NKMMAX,CWORK(1,1,IE),
     &                                 'CLM>REL',CWORK1)
C
                        DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
     &                     = CWORK1(1:NKMMAX,1:NKMMAX)
                        IF ( IPRINTDMFT.GT.0 .AND. IE.EQ.0 ) THEN
                           WRITE (6,*) 'IT=',IT,'IE=',IE
                           CALL CMATSTRUCT('DMFTSIGMA_REL',
     &                        DMFTSIGMA(1,1,IT,IE),NKMMAX,NKMMAX,IREL,
     &                        IREL,0,1D-8,6)
                        END IF
                     END DO
                  ELSE
                     IF ( .NOT.ALLOCATED(CWORK1) )
     &                    ALLOCATE (CWORK1(NKMMAX,NKMMAX))
C
                     CWORK1 = C0
                     DO IE = 1,NETAB(1)
C
                        CALL CHANGEREP(NKM,NKMMAX,CWORK(1,1,IE),
     &                                 'CLM>RLM',CWORK1)
C
                        DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
     &                     = CWORK1(1:NKMMAX,1:NKMMAX)
                        IF ( IPRINTDMFT.GT.0 .AND. IE.EQ.0 ) THEN
                           WRITE (6,*) 'IT=',IT,'IE=',IE
                           CALL CMATSTRUCT('DMFTSIGMA_RLM',
     &                        DMFTSIGMA(1,1,IT,IE),NKMMAX,NKMMAX,IREL,
     &                        IREL,0,1D-8,6)
                        END IF
                     END DO
C
                  END IF
                  CALL CPU_TIME(TIME)
                  WRITE (6,99010) IT,TIME - TIME0
C=======================================================
C     *    *    *
C=======================================================
               END IF           !MPI
            END IF              !KSELF
C
         END DO                 !IT
c modified by XJQ: cpu time
         call cpu_time(t1_dmft)
         write(*,*) 'dmft takes ',t1_dmft-t0_dmft,' s in mpi_id ',mpi_id
c end-mod-xjq
         IF ( NPROCS>1 ) THEN
            CALL DRV_MPI_BARRIER
            IF ( MPI_TYPES ) THEN
               IF ( ALLOCATED(CWKT) ) DEALLOCATE (CWKT)
               ALLOCATE (CWKT(NKMMAX))
               LWKKMTE = NKMMAX
               DO IVEC = 1,NKMMAX
                  DO IE = 1,NEMAX
                     DO IT = ITBOT,ITTOP
                        CALL DRV_MPI_REDUCE_C(DMFTSIGMA(1,IVEC,IT,IE),
     &                     CWKT(1),LWKKMTE)
                     END DO
                  END DO
               END DO
               DEALLOCATE (CWKT)
               IF ( ALLOCATED(RWKT) ) DEALLOCATE (RWKT)
               ALLOCATE (RWKT(NTMAX))
               CALL DRV_MPI_REDUCE_R(SIGERRT,RWKT,NTMAX)
c modified by XJQ: parallel on types using my own subroutines
            elseif(lparalleled) then
              call mpi_barrier(parent_comm,mpierr)
              call get_comm_level(1,'inter ',lparalleled,
     &                            inter_comm,nprocs_inter,inter_rank)
c colloct it0_collct, it1_collct and nt_collct for further collecting
              allocate(disp_proc(0:nprocs_inter-1))
              allocate(nt_rev_proc(0:nprocs_inter-1))
              it0_collct = it_itdmft(it0_distr)
              it1_collct = it_itdmft(it1_distr)
              nt_collct = it1_collct - it0_collct + 1
              call collect_varstartnvar(1,inter_comm,nprocs_inter,
     &                                  it0_collct,nt_collct,
     &                                  disp_proc(0:),nt_rev_proc(0:))
c
              call mpi_barrier(parent_comm,mpierr)
              allocate(dmftsig_tmp(nkmmax,nkmmax,itbot:ittop))
              dmftsig_tmp(:,:,itbot:ittop) = c0
              call mpi_type_contiguous(nkmmax*nkmmax,mpi_double_complex,
     &                                 mpi_row,mpierr)
              call mpi_type_commit(mpi_row,mpierr)
              do ie=1,nemax
                call mpi_gatherv(dmftsigma(1,1,it0_collct,ie),nt_collct,
     &                           mpi_row,dmftsig_tmp(1,1,1),
     &                           nt_rev_proc(0:),disp_proc(0:),
     &                           mpi_row,0,inter_comm,mpierr)
                call mpi_barrier(parent_comm,mpierr)
                if(parent_rank==0) 
     &            dmftsigma(:,:,itbot:ittop,ie) = 
     &              dmftsig_tmp(:,:,itbot:ittop)
              enddo
c
              call mpi_barrier(parent_comm,mpierr)
              IF ( ALLOCATED(RWKT) ) DEALLOCATE (RWKT)
              ALLOCATE (RWKT(NTMAX))
              rwkt(:)=c0
              call mpi_gatherv(sigerrt(it0_collct),nt_collct,
     &                         mpi_double_precision,rwkt(1),
     &                         nt_rev_proc(0:),disp_proc(0:),
     &                         mpi_double_precision,
     &                         0,inter_comm,mpierr)
              call mpi_barrier(parent_comm,mpierr)
              if(parent_rank==0) sigerrt(:) = rwkt(:)
c
              call mpi_barrier(parent_comm,mpierr)
              deallocate(dmftsig_tmp)
              deallocate(disp_proc)
              deallocate(nt_rev_proc)
              lothermpi=.false.
              call mpi_multilevel_free(1)
c end-mod-xjq
            END IF
C
            CALL DRV_MPI_BARRIER
C
            CALL DRV_MPI_BCAST_C(0,DMFTSIGMA(1,1,1,1),
     &                           NKMMAX*NKMMAX*NTMAX*NETAB(1))
         ENDIF
CC=======================================================
C-----------------------------------------------------------------------
C               apply damping to the self energy
C-----------------------------------------------------------------------
         IF ( LDAMP ) THEN
            FDAMP = 1D0 - EXP(-KDAMP*ITRSCF)
            IF ( FDAMP.LT.0.97D0 ) THEN
               WRITE (6,99005) FDAMP
               DO IT = 1,NTMAX
                  DO IE = 1,NETAB(1)
                     DO I = 1,NKMMAX
                        DO J = 1,NKMMAX
                           DMFTSIGMA(J,I,IT,IE)
     &                        = FDAMP*DMFTSIGMA(J,I,IT,IE)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END IF
C-----------------------------------------------------------------------
C
C ITRDMFT
C
         IF ( DMFTSCF ) THEN
            WRITE (6,99003)
            SIGERR = 0.0D0
            RTMP = 0.0D0
            DO IT = ITBOT,ITTOP
               IF ( KSELF(IT).GT.0 ) THEN
                  SIGERR = SIGERR + SIGERRT(IT)
                  RTMP = RTMP + 1.0D0
                  WRITE (6,99009) IT,SIGERRT(IT)
               END IF
            END DO
            SIGERR = SIGERR/RTMP
            IF ( SIGERR.LE.SIGTOL ) THEN
               WRITE (6,99007) ITRDMFT,SIGERR,
     &                         'DMFT - cycle converged !!!!!!!!! '
            ELSE IF ( ITRDMFT.LE.NITERDMFT ) THEN
               WRITE (6,99007) ITRDMFT,SIGERR,
     &                         'DMFT - cycle not converged  -  continue'
               CYCLE
            ELSE
               WRITE (6,99007) ITRDMFT,SIGERR,
     &               'DMFT - cycle not converged - NITERDMFT exhausted '
            END IF
            EXIT
         END IF
      END DO
C=======================================================
C Print Galitsky-Migdal energy, obtained on semicircle:
C=======================================================
      IF ( DMFT_GALMIG ) THEN
         WRITE (*,99004)
         DO IT = ITBOT,ITTOP
            WRITE (*,'(10X,2I5,2X,2F9.4)') ITRSCF,IT,DREAL(CEMIG(IT)),
     &             DIMAG(CEMIG(IT))
         END DO
         DO IT = ITBOT + 1,ITTOP
            CEMIG(ITBOT) = CEMIG(ITBOT) + CEMIG(IT)
         END DO
         WRITE (*,'(10X,30A)') (('-'),I=1,30)
         WRITE (*,'(10X,I5,A7,2F9.4)') ITRSCF,'  E_GM:',DREAL(CEMIG(1)),
     &                                 DIMAG(CEMIG(1))
         WRITE (*,'(10X,69A)') (('-'),I=1,69)
      END IF
C
      DEALLOCATE (CEMIG)
C=======================================================
      WRITE (6,99003)
C
C=======================================================
C Calculate new reference energy
C=======================================================
C
C
      CALL DMFT_LDAU_EREF(OBS_LT,EREFDMFT)
C
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).EQ.1 ) WRITE (6,99008) IT,DREAL(EREFDMFT(IT))
      END DO
      WRITE (6,'(/)')
C
      CALL CPU_TIME(TIME_ALL)
      WRITE (6,99011) TIME_ALL - TIME_ALL0
C
C      IF ( SIGERR.LE.SIGTOL .AND. ITR.GT.1 ) DMFT_FIX_DYN_SE = .TRUE.
      IF ( DMFTSCF ) STOP
C
      RETURN
C
99001 FORMAT (/,1X,79('*'),/,25X,'DMFT SOLVER:    ',A,/,1X,79('*'))
99002 FORMAT (/,1X,5X,'DMFT calculation for IT= ',I3,'  LOP =',I2,2X,
     &        'U_eff =',F5.2,' eV ',2X,'J_eff =',F5.2,' eV ')
99003 FORMAT (/,1X,79('*'),/)
99004 FORMAT (/,10X,69('-'),/,10X,'Galitsky-Migdal contribution (Ry):',
     &        /,10X,30('-'),/,10X,'ITRSCF  IT    Re E_GM  Im E_GM   ')
99005 FORMAT (10X,'damping factor for SCF - cycle    FDAMP =',F6.3,/)
99006 FORMAT (/,1X,79('*'),/,25X,'ONLY DMFT SELF CONSISTENCY',/,79('*'))
99007 FORMAT (/,I4,' ERR',1E10.3,/,12x,A,/,1x,79('*'),/)
99008 FORMAT (10X,'DMFT reference energy for IT',I3,' =',F11.8)
99009 FORMAT (5X,'error for type',i3,':  Sigma = ',1P,D11.4)
99010 FORMAT (/,1X,5X,'DMFT execution time for type',i3,':',F14.3,
     &        ' secs',/)
99011 FORMAT (/,1X,5X,'DMFT total execution time :',F14.3,' secs',/)
C
      END
C
