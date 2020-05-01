C*==straa.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE STRAA(IPRINT,RMAX,NRDL,SMAX,NGRL,NQQP_STR,LLMAX,ETA,
     &                 ALPHA0,ABAS,BBAS,QQPX,QQPY,QQPZ,EXPGNQ,QQMLRS,G1,
     &                 G2,G3,G123MAX,R1,R2,R3,R123MAX,NSDL,INDR,HP,
     &                 GGJLRS,NGRLMAX,NRDLMAX,J13MAX,J22MAX,NQQP_STRMAX,
     &                 LLARR,NLLMMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  CALCULATE ALL QUANTITIES INDEPENDENT OF  ENERGY  AND  K         *
C   *                THIS ROUTINE IS CALLED ONLY ONCE                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_STR,ONLY:OPTIMIZE_LIMITS_JXXMAX,J13MAX_LOWER_LIMIT,
     &    J13MAX_UPPER_LIMIT,J22MAX_LOWER_LIMIT,J22MAX_UPPER_LIMIT
      USE MOD_CONSTANTS,ONLY:C0,CI,PI
      USE MOD_ENERGY,ONLY:EMIN,EMAX
      USE MOD_LATTICE,ONLY:ALAT
c modified by XJQ: parallel on NQQP_STR
      use mod_files, only : iotmp
      use mpi
      use mod_mpi_multilevels
c end-mod-xjq
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRAA')
      REAL*8 TOL
      PARAMETER (TOL=1D-14)
C
C Dummy arguments
C
      REAL*8 ALPHA0,ETA,RMAX
      INTEGER G123MAX,IPRINT,J13MAX,J22MAX,LLARR,LLMAX,NGRL,NGRLMAX,
     &        NLLMMMAX,NQQP_STR,NQQP_STRMAX,NRDL,NRDLMAX,NSDL,R123MAX
      REAL*8 ABAS(3,3),BBAS(3,3),
     &       GGJLRS(-J22MAX_UPPER_LIMIT:LLARR,NSDL,NQQP_STR),
     &       HP(NLLMMMAX),QQPX(NQQP_STRMAX),QQPY(NQQP_STRMAX),
     &       QQPZ(NQQP_STRMAX)
      COMPLEX*16 EXPGNQ(NGRL,NQQP_STR),QQMLRS(NLLMMMAX,NSDL,NQQP_STR)
      INTEGER G1(NGRLMAX),G2(NGRLMAX),G3(NGRLMAX),INDR(NSDL,NQQP_STR),
     &        R1(NRDLMAX),R2(NRDLMAX),R3(NRDLMAX),SMAX(NQQP_STR)
C
C Local variables
C
      COMPLEX*16 ALPHA,CADD,CSUM,D300,EDU,EDULIM(2),EHOCHJ
      REAL*8 ATVOL,D1TERM1,F,FAKTOR,GX,GY,GZ,Q1,Q2,Q3,RELDIFF,RS,RSQUAD,
     &       RX,RX0,RY,RY0,RZ,RZ0
      LOGICAL C_LT_EPS
      INTEGER I,I1,IE,IQQP,J,J13,J22,LL,LL_MIN_J22,MM,MMLL,N,S
      REAL*8 STRCONFRA
c modified by XJQ: parallel on NQQP_STR
      logical lparalleled, lopen
      integer iqqp0, iqqp1, nqqp, nqqp_max, mpierr, 
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      real*8 t0_straa, t1_straa
      integer :: r123max_tmp
      integer :: nrev, iqq0_tmp
      integer, dimension(:), allocatable :: disp_proc, nqqp_proc
      integer, dimension(:,:), allocatable :: indr_tmp!(nsdl,nqqp_str)
      real*8, dimension(:,:,:), allocatable :: ggjlrs_tmp!(-j22max_upper_limit:llarr,nsdl,nqqp_str)
      complex*16, dimension(:,:,:), allocatable :: qqmlrs_tmp!(nllmmmax,nsdl,nqqp_str)
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
c     DO J = 1,NQQP_STR
c        DO I = 1,NSDL
c           INDR(I,J) = 0
c        END DO
c     END DO
C
      INDR(1:NSDL,1:NQQP_STR) = 0

      EDULIM(1) = (EMIN-0.2D0)/(2*PI/ALAT)**2
      EDULIM(2) = (EMAX+0.2D0)/(2*PI/ALAT)**2
C
C  =====================================================================
C                      ********
C                      * DLM1 *
C                      ********
C
      RX = ABAS(2,2)*ABAS(3,3) - ABAS(3,2)*ABAS(2,3)
      RY = ABAS(3,2)*ABAS(1,3) - ABAS(1,2)*ABAS(3,3)
      RZ = ABAS(1,2)*ABAS(2,3) - ABAS(2,2)*ABAS(1,3)
C
      ATVOL = ABS(ABAS(1,1)*RX+ABAS(2,1)*RY+ABAS(3,1)*RZ)*(2*PI)**3
C
      D1TERM1 = -4.0D0*PI/ATVOL
C
      G123MAX = 0
      DO N = 1,NGRL
         G123MAX = MAX(G123MAX,ABS(G1(N)),ABS(G2(N)),ABS(G3(N)))
C
         GX = G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)*BBAS(1,3)
         GY = G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)*BBAS(2,3)
         GZ = G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)*BBAS(3,3)
C
         F = D1TERM1*EXP(-(GX**2+GY**2+GZ**2)/ETA)
C
         DO IQQP = 1,NQQP_STR
            EXPGNQ(N,IQQP) = F*CDEXP(CI*2*PI*(GX*QQPX(IQQP)+GY*QQPY(IQQP
     &                       )+GZ*QQPZ(IQQP)))
         END DO
      END DO
C
C  =====================================================================
C                      ********
C                      * DLM2 *
C                      ********
C
      Q1 = -SQRT(ETA/PI)*0.5D0
C
c modified by XJQ: parallel on NQQP_STR
      R123MAX_TMP = 0
c
      call get_comm_level(1,'parent',lparalleled,
     &                    parent_comm,nprocs_parent,parent_rank)
      if(parent_rank==0) call cpu_time(t0_straa)
c when nqqp_str<nprocs, which usually not huge, not necessary to parallel
      if(nprocs_parent < nqqp_str) then
      call mpi_multilevel_distribute(routine,1,nqqp_str,
     &                               lparalleled,iqqp0,iqqp1)
      else
      lparalleled=.false.
      iqqp0=1
      iqqp1=nqqp_str
      endif
c
      allocate( indr_tmp( 1:nsdl, iqqp0:iqqp1 ) )
c allocate one more iqqp for later collecting
      allocate( ggjlrs_tmp( -j22max_upper_limit:llarr, 1:nsdl,
     &                      iqqp0:iqqp1+1 ) )
      allocate( qqmlrs_tmp( 1:nllmmmax, 1:nsdl, iqqp0:iqqp1 ) )
c end-mod-xjq
C
c     DO IQQP = 1,NQQP_STR
      DO IQQP = iqqp0, iqqp1
c modified by XJQ: parallel on NQQP_STR
         if(iqqp >= iqqp0 .and. iqqp <= iqqp1) then
c end-mod-xjq
C
C ---> suppress contribution for ->R=0 and IQQP=1 (IQ=IQP)
C      using the fact that the ->R's have been ordered
         IF ( IQQP.EQ.1 ) THEN
            SMAX(1) = SMAX(1) - 1
            I1 = 2
         ELSE
            I1 = 1
         END IF
         S = 0
C
         DO I = I1,NRDL
C
            RX0 = R1(I)*ABAS(1,1) + R2(I)*ABAS(1,2) + R3(I)*ABAS(1,3)
     &            - QQPX(IQQP)
            RY0 = R1(I)*ABAS(2,1) + R2(I)*ABAS(2,2) + R3(I)*ABAS(2,3)
     &            - QQPY(IQQP)
            RZ0 = R1(I)*ABAS(3,1) + R2(I)*ABAS(3,2) + R3(I)*ABAS(3,3)
     &            - QQPZ(IQQP)
            RX = 2*PI*RX0
            RY = 2*PI*RY0
            RZ = 2*PI*RZ0
C
            RSQUAD = RX**2 + RY**2 + RZ**2
            RS = DSQRT(RX0**2+RY0**2+RZ0**2)
C
            IF ( RS.LE.RMAX ) THEN
C  ---------------------------------------------------------------------
C                         ->R    ACCEPTED
C  ---------------------------------------------------------------------
c modified by XJQ: parallel on NQQP_STR
c               R123MAX = MAX(R123MAX,ABS(R1(I)),ABS(R2(I)),ABS(R3(I)))
               R123MAX_TMP = MAX(R123MAX_TMP,ABS(R1(I)),ABS(R2(I)),
     &                           ABS(R3(I)))
C
               S = S + 1
c modified by XJQ: parallel on NQQP_STR
c               INDR(S,IQQP) = I
               INDR_TMP(S,IQQP) = I
c end-mod-xjq
C
               CALL CALC_RHPLM(RX,RY,RZ,HP,LLMAX,NLLMMMAX)
C
C              Q2 = (-ETA/2.0D0)**LL
               Q2 = 1.0D0/(-ETA/2.0D0)
               Q3 = EXP(-ETA*RSQUAD/4.0D0)
C
               MMLL = 0
C
               DO LL = 0,LLMAX
                  Q2 = Q2*(-ETA/2.0D0)
                  F = Q1*Q2*Q3
C
                  DO MM = -LL,LL
                     MMLL = MMLL + 1
c modified by XJQ: parallel on NQQP_STR
c                     QQMLRS(MMLL,S,IQQP) = F*HP(MMLL)
                     QQMLRS_TMP(MMLL,S,IQQP) = F*HP(MMLL)
c end-mod-xjq
                  END DO
               END DO
C
C
C  /D (23)/
C
               DO LL_MIN_J22 = -J22MAX_UPPER_LIMIT,LLMAX
c modified by XJQ: parallel on NQQP_STR
c                  GGJLRS(LL_MIN_J22,S,IQQP)
                  GGJLRS_TMP(LL_MIN_J22,S,IQQP)
c end-mod-xjq
     &               = STRCONFRA((LL_MIN_J22+0.5D0),(RSQUAD*ETA/4.0D0))
               END DO
C
            END IF
         END DO
         IF ( S.NE.SMAX(IQQP) ) THEN
            WRITE (6,*) ' WARNING FROM <STRAA> '
            WRITE (6,*) ' IQQP = ',IQQP
            WRITE (6,*) ' S    = ',S
            WRITE (6,*) ' SMAX = ',SMAX(IQQP)
         END IF
c modified by XJQ: parallel on NQQP_STR
         endif ! iqqp >= iqqp0 .and. iqqp <= iqqp1
c end-mod-xjq
      END DO
c modified by XJQ: parallel on NQQP_STR, collecting processes
      if(.not. lparalleled) then
        r123max = r123max_tmp
        indr(:,:) = indr_tmp(:,:)
        qqmlrs(:,:,:) = qqmlrs_tmp(:,:,:)
        ggjlrs(:,:,:) = ggjlrs_tmp(:,:,:)
        deallocate(indr_tmp)
        deallocate(ggjlrs_tmp)
        deallocate(qqmlrs_tmp)
      else
c        call get_comm_level(1,'parent',lparalleled,
c     &                      parent_comm,nprocs_parent,parent_rank)
c
        call mpi_barrier(parent_comm,mpierr)
        call get_comm_level(1,'inter ',lparalleled,
     &                      inter_comm,nprocs_inter,inter_rank)
c colloct iqqp0, iqqp1 and nqqp for further collecting
        allocate(disp_proc(0:nprocs_inter-1))
        allocate(nqqp_proc(0:nprocs_inter-1))
        nqqp = iqqp1-iqqp0+1
        call collect_varstartnvar(1,inter_comm,nprocs_inter,
     &                            iqqp0,nqqp,
     &                            disp_proc(0:),nqqp_proc(0:))
c r123max
        call mpi_allreduce(r123max_tmp,r123max,1,mpi_int,
     &                     mpi_max,inter_comm,mpierr)
c collect smax in inter_comm
        if(iqqp0/=1) smax(1) = smax(1) - 1
c collect indr in inter_comm
        call mpi_barrier(parent_comm,mpierr)
        call mpi_type_contiguous(nsdl,mpi_int,mpi_row,mpierr)
        call mpi_type_commit(mpi_row,mpierr)
c        call mpi_gatherv(indr_tmp(1,iqqp0),nqqp,mpi_row,
c     &                   indr(1,1),nqqp_proc(0:),disp_proc(0:),
c     &                   mpi_row,0,inter_comm,mpierr)
        call mpi_allgatherv(indr_tmp(1,iqqp0),nqqp,mpi_row,
     &                      indr(1,1),nqqp_proc(0:),disp_proc(0:),
     &                      mpi_row,inter_comm,mpierr)
        call mpi_type_free(mpi_row,mpierr)
        deallocate(indr_tmp)
c collect qqmlrs in inter_comm
        call mpi_barrier(parent_comm,mpierr)
        call mpi_type_contiguous(nllmmmax*nsdl,mpi_double_complex,
     &                           mpi_row,mpierr)
        call mpi_type_commit(mpi_row,mpierr)
c        call mpi_gatherv(qqmlrs_tmp(1,1,iqqp0),nqqp,mpi_row,
c     &                   qqmlrs(1,1,1),nqqp_proc(0:),disp_proc(0:),
c     &                   mpi_row,0,inter_comm,mpierr)
        call mpi_allgatherv(qqmlrs_tmp(1,1,iqqp0),nqqp,mpi_row,
     &                      qqmlrs(1,1,1),nqqp_proc(0:),disp_proc(0:),
     &                      mpi_row,inter_comm,mpierr)
        call mpi_type_free(mpi_row,mpierr)
        deallocate(qqmlrs_tmp)
c collect ggjlrs in inter_comm
        call mpi_barrier(parent_comm,mpierr)
        call mpi_type_contiguous((llarr+j22max_upper_limit+1)*nsdl,
     &                           mpi_double_precision,mpi_row,mpierr)
        call mpi_type_commit(mpi_row,mpierr)
c        call mpi_gatherv(ggjlrs_tmp(-j22max_upper_limit,1,iqqp0),nqqp,
c     &                   mpi_row,ggjlrs(-j22max_upper_limit,1,1),
c     &                   nqqp_proc(0:),disp_proc(0:),
c     &                   mpi_row,0,inter_comm,mpierr)
c allgatherv for only 1 qqp to reduce the needed memory
        call mpi_allreduce(nqqp,nqqp_max,1,mpi_int,mpi_max,
     &                     inter_comm,mpierr)
        do iqqp=iqqp0,iqqp0+nqqp_max-1
        if(iqqp<=iqqp1) then
          nrev=1
        else
          nrev=0
        endif
        call collect_varstartnvar(1,inter_comm,nprocs_inter,
     &                            iqqp,nrev,
     &                            disp_proc(0:),nqqp_proc(0:))
c        call mpi_allgatherv(ggjlrs_tmp(-j22max_upper_limit,1,iqqp0),
c     &                      nqqp,
c     &                      mpi_row,ggjlrs(-j22max_upper_limit,1,1),
c     &                      nqqp_proc(0:),disp_proc(0:),
c     &                      mpi_row,inter_comm,mpierr)
        call mpi_allgatherv(ggjlrs_tmp(-j22max_upper_limit,1,iqqp),
     &                      1,
     &                      mpi_row,ggjlrs(-j22max_upper_limit,1,1),
     &                      nqqp_proc(0:),disp_proc(0:),
     &                      mpi_row,inter_comm,mpierr)
c        call mpi_barrier(parent_comm,mpierr)
        enddo
        call mpi_type_free(mpi_row,mpierr)
        deallocate(ggjlrs_tmp)
c bcast to all processes in parent_comm
c        call mpi_barrier(parent_comm,mpierr)
c        call mpi_bcast(r123max,1,mpi_int,0,parent_comm,mpierr)
c        call mpi_barrier(parent_comm,mpierr)
c        call mpi_bcast(indr,nsdl*nqqp_str,mpi_int,0,parent_comm,mpierr)
c matrices are too large to bcast once
c        call mpi_barrier(parent_comm,mpierr)
c        do iqqp = 1, nqqp_str
c          call mpi_bcast(qqmlrs(:,:,iqqp),nllmmmax*nsdl,
c     &                   mpi_double_complex,0,parent_comm,mpierr)
c        enddo
c        call mpi_barrier(parent_comm,mpierr)
c        do iqqp = 1, nqqp_str
c          call mpi_bcast(ggjlrs(:,:,iqqp),
c     &                   (llarr+j22max_upper_limit+1)*nsdl,
c     &                   mpi_double_precision,0,parent_comm,mpierr)
c        enddo
        call mpi_barrier(parent_comm,mpierr)
c free mpi_distribution on iqqp
        call mpi_multilevel_free(1)
      endif ! lparallel
c debug
      if(.false.) then
      if(parent_rank==0) write(*,*) 'nsdl = ',nsdl
      inquire(iotmp,opened=lopen)
      if(lopen) stop 'iotmp is opened, try another for debugging'
      if(nprocs_parent==1) then
        open(iotmp,file='debug_indr_nompi.dat',status='replace')
        do iqqp=1,nqqp_str
          write(iotmp,'(2i6)') indr(nsdl/2,iqqp),indr(nsdl,iqqp)
        enddo
        close(iotmp)
      endif
      if(nprocs_parent>2) then
        if(parent_rank==0) then
          open(iotmp,file='debug_indr_rank0.dat',status='replace')
          do iqqp=1,nqqp_str
            write(iotmp,'(2i6)') indr(nsdl/2,iqqp),indr(nsdl,iqqp)
          enddo
          close(iotmp)
        endif
        if(parent_rank==(nprocs_parent/2)) then
          open(iotmp,file='debug_indr_rankmid.dat',status='replace')
          do iqqp=1,nqqp_str
            write(iotmp,'(2i6)') indr(nsdl/2,iqqp),indr(nsdl,iqqp)
          enddo
          close(iotmp)
        endif
        if(parent_rank==nprocs_parent-1) then
          open(iotmp,file='debug_indr_rankend.dat',status='replace')
          do iqqp=1,nqqp_str
            write(iotmp,'(2i6)') indr(nsdl/2,iqqp),indr(nsdl,iqqp)
          enddo
          close(iotmp)
        endif
      endif
      endif ! .false.
c debug
      call mpi_barrier(parent_comm,mpierr)
      if(parent_rank==0) then
        call cpu_time(t1_straa)
        write(*,*)
        write(*,*) 'straa takes ',t1_straa-t0_straa,' seconds'
        write(*,*)
      endif
c end-mod-xjq
C
      J22MAX = J22MAX_LOWER_LIMIT
C
      IF ( IPRINT.GE.3 ) WRITE (6,99004) J22MAX_UPPER_LIMIT,J22MAX
C
C--------------------------------------------- search optimum for J22MAX
      IF ( OPTIMIZE_LIMITS_JXXMAX ) THEN
C
         RELDIFF = 0D0
         DO IQQP = 1,NQQP_STR
            DO S = 1,SMAX(IQQP)
               DO IE = 1,2
C
                  EDU = EDULIM(IE)
C
                  DO LL = 0,LLMAX
                     FAKTOR = 1.0D0
                     EHOCHJ = 1.0D0
                     CSUM = C0
                     DO J22 = 0,J22MAX_UPPER_LIMIT
                        LL_MIN_J22 = LL - J22
C
                        CADD = EHOCHJ*GGJLRS(LL_MIN_J22,S,IQQP)*FAKTOR
                        CSUM = CSUM + CADD
C
                        IF ( C_LT_EPS((CADD/CSUM),TOL) ) THEN
                           J22MAX = MAX(J22MAX,J22)
                           EXIT
                        END IF
C
                        FAKTOR = FAKTOR/(ETA*(J22+1.0D0))
                        EHOCHJ = EHOCHJ*EDU
                     END DO
                     RELDIFF = MAX(RELDIFF,ABS(CADD/CSUM))
                  END DO
C
               END DO
            END DO
         END DO
         IF ( RELDIFF.GT.TOL ) THEN
            WRITE (6,99005) EMIN,EMAX,ETA,TOL,RELDIFF,
     &                      'J22MAX   J22MAX ',J22MAX
C
            CALL STOP_MESSAGE(ROUTINE,'no upper limit J22MAX found')
         END IF
C
      END IF
C--------------------------------------- search optimum for J22MAX - end
C
C  =====================================================================
C                      ********
C                      * DLM3 *
C                      ********
C  /D (8),(13)/
C
      ALPHA0 = SQRT(ETA)/(2.0D0*PI)
C
C-----------------------------------------------------------------------
C      set upper limit J13MAX for series expansion - see <STRCC>
C-----------------------------------------------------------------------
C
      J13MAX = J13MAX_LOWER_LIMIT
C
C--------------------------------------------- search optimum for J13MAX
      IF ( OPTIMIZE_LIMITS_JXXMAX ) THEN
C
         DO IE = 1,2
C
            EDU = EDULIM(IE)
C
            D300 = 0.0D0
            EHOCHJ = 1.0D000
            ALPHA = ALPHA0
            J13 = -1
 20         CONTINUE
            J13 = J13 + 1
            D300 = ALPHA*EHOCHJ + D300
            ALPHA = ALPHA*(2.0D0*J13-1.0D0)
     &              /(ETA*(J13+1.0D0)*(2.0D0*J13+1.0D0))
            EHOCHJ = EHOCHJ*EDU
C
C     prevent floating point underflow
            IF ( C_LT_EPS(EHOCHJ,1D-50) ) EHOCHJ = C0
C
            IF ( J13.LT.J13MAX_UPPER_LIMIT ) THEN
C
               IF ( .NOT.C_LT_EPS((ALPHA*EHOCHJ/D300),TOL) ) GOTO 20
               J13MAX = MAX(J13MAX_LOWER_LIMIT,J13)
C
            ELSE
C
               WRITE (6,99005) EMIN,EMAX,ETA,TOL,ABS(ALPHA*EHOCHJ/D300),
     &                         'J13MAX    ',J13MAX
               CALL STOP_MESSAGE(ROUTINE,'no upper limit J13MAX found')
C
            END IF
C
         END DO
C
      END IF
C--------------------------------------- search optimum for J13MAX - end
C
C  =====================================================================
      IF ( IPRINT.GE.3 ) THEN
         WRITE (6,'(''  LATTICE VECTORS OF THE RECIPROCAL LATTICE'',/)')
         DO N = 1,NGRL
            WRITE (6,99001) N,G1(N),G2(N),G3(N)
         END DO
C
         WRITE (6,'(//,''  LATTICE VECTORS OF THE DIRECT LATTICE '',/)')
         DO S = 1,NRDL
            WRITE (6,99002) S,R1(S),R2(S),R3(S)
         END DO
C
         WRITE (6,'(//,''  R-INDEX TABLE FOR OFF-DIAG. ELEMENTS  '',/)')
         DO S = 1,NSDL
            WRITE (6,99003) S,(INDR(S,IQQP),IQQP=1,NQQP_STR)
         END DO
      END IF
C
C  =====================================================================
99001 FORMAT ('  ->GN(',I3,')=  (',I3,',',I3,',',I3,')')
99002 FORMAT ('  ->RS(',I3,') =  (',I3,',',I3,',',I3,')')
99003 FORMAT ('  INDR(',I3,'):',30I4)
99004 FORMAT (/,10X,'fixing upper boundary J22MAX for series expansion',
     &        /,10X,'LIMIT   ',I8,10X,'START   ',I8)
99005 FORMAT (/,10X,'<STRAA> - setting boundaries for series expansion',
     &        /,10X,'EMIN      ',F10.5,/,10X,'EMAX      ',F10.5,/,10X,
     &        'ETA       ',F10.5,/,10X,'TOL       ',E13.4,/,10X,
     &        'rel. diff.',E13.4,/,10X,A,I10)
      END
