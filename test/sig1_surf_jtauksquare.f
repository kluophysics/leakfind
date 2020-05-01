C*==sig1_surf.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_SURF_JTAUKSQUARE(
     &             ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQBZ,MSSQA,MSSQB,
     &             MAQAB,MBQAB,MCQAB,MDQAB,SIG1Q,KEY,MAQBA,
     &             MCQBA,MBQBA,MDQBA,NSPINPROJ)
C***********************************************************************
C   Not Use CHIZ (so works only when No CPA)
C   The multiply operators for sigCB = K*[ 2*L*(L*Q)^2 + (L*Q)^2 ]
C   should be much smaller than that for sigkloop+sig1 and sig1_notstorechiz
C
C   sig1 = sig1_cpa + sig1_k
C   sig1_k = J_{L4,L1,IQ} TAUK_{L1,IQ,L2,JQ} TAUK_{L3,JQ,L4,IQ} J_{L2,L3,JQ}
C
C   reorder the summation, we will get
C   sig1_k = 0.25 * sum_k sum_{L4,IQ,L2,JQ} [ sigCB - conjg(sigDD) - sigAA ]
C   sigCB = JCTAUK_{L4,IQ,L2,JQ} JBTAUK_{L2,JQ,L4,IQ}
C   sigDD = JD1TAUK_{L4,IQ,L2,JQ} JD2TAUK_{L2,JQ,L4,IQ}
C   JCTAUK_{L4,IQ,L2,JQ} = { sum_{L1} [ JC_{L4,L1,IQ} TAUK_{L1,IQ,L2,JQ} ] }
C***********************************************************************
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,NZ12,NZ12MAX
      USE MOD_KSPACE,ONLY:WKTAB,NKTAB,KTAB,WKSUM
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,LMAT3,IND0Q,NKMQ,NKKR,WKM1
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SOTSI,EESI,
     &    IRESPONSE_SOT,IRESPONSE_EDELSTEIN
      USE MOD_CONSTANTS,ONLY:C0,CI,C1
      USE MOD_CALCMODE,ONLY:PUBLIC_VERSION,IREL
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      use mod_files,only:iotmp
      use mpi
      use mod_mpi_multilevels
      IMPLICIT NONE
C*--SIG1_SURF27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG1_SURF_JTAUKSQUARE')
C
C Dummy arguments
C
      COMPLEX*16 ERYDA,ERYDB,PA,PB
      COMPLEX*16 MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           TAUQA(NKMMAX,NKMMAX,NQMAX),TAUQB(NKMMAX,NKMMAX,NQMAX),
     &           TAUQBZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX)
      CHARACTER*1 KEY
      INTEGER NSPINPROJ
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           SIG1Q(3,3,NSPINPROJ,NQMAX,NQMAX)
C
C Local variables
C
      INTEGER I,IA_ERR,IQ,IQCHI,ISPINPROJ,ISPR,J,JQ,JQCHI,K1,K2,L1,L2,
     &        L3,L4,MUE,NKMSQ,NUE,Z1,Z2,IPIV(:),INFO
      INTEGER JTT10,JTT20,JTT1,JTT2,JTTX,K14I,K23J
      COMPLEX*16 CWK,MAUX(:,:),
     &           TAUK1,TAUK2,TAUKLINA(:),TAUKLINB(:),TAUKX,CHIZ(:,:,:)
      COMPLEX*16 MBQAB_L(:,:,:,:,:),MBQBA_L(:,:,:,:,:),
     &           MCQAB_L(:,:,:,:,:),MCQBA_L(:,:,:,:,:),MDQABX(:,:,:,:,:)
     &           ,MDQAB_L(:,:,:,:,:),MDQBAX(:,:,:,:,:),
     &           MDQBA_L(:,:,:,:,:),S11_1,S11_2,S12_1,S12_2,S21_1,S21_2,
     &           S22_1,S22_2,SIG1IIQ(:,:,:,:),SIG1IRQ(:,:,:,:),
     &           SUMX_1(3,3,2,2),SUMX_2(3,3,2,2),SUM_A(3,3),SUM_ALL(3,3)
     &           ,SUM_CB(3,3),SUM_D(3,3),S_A,S_CB,S_D,TMP_1(3,3),
     &           TMP_2(3,3)
      complex*16 ctmp1, ctmp2, ctmp3
c parallel
      logical lparalleled, lopen
      integer level, ik0, ik1, nkp, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      integer ik
      real*8 t0_sig1, t1_sig1, wk
      complex*16, dimension(:,:,:,:,:), allocatable :: sig1q_tmp
c jtauk
      complex*16 zdotu
      complex*16,dimension(3,3,nspinproj,nqmax) :: sig1_cpa
      complex*16,dimension(nkmmax,nkmmax) :: tauka_ij, taukx_ij,
     &  taukb_ji, j1tauk, j2tauk_t
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIG1IIQ,SIG1IRQ,MDQBAX,MDQABX
      ALLOCATABLE MBQAB_L,MBQBA_L,MCQAB_L,MCQBA_L,MDQAB_L,MDQBA_L
      ALLOCATABLE MAUX,TAUKLINA,TAUKLINB,CHIZ,IPIV
C
      CALL TRACK_INFO(ROUTINE)
C
      NKMSQ = NKM*NKM
C
      ALLOCATE (MDQBAX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQABX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
      ALLOCATE (MBQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MBQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MCQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MCQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
C-----------------------------------------------------------------------
C  multiply the averaged MEs with proper LMAT to take into account
C  that CHIZ does not incude the LMATs
C-----------------------------------------------------------------------
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               MBQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(MBQAB(:,:,MUE,ISPINPROJ,IQ),LMAT3)
               MBQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MBQBA(:,:,MUE,ISPINPROJ,IQ))
C
               MCQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(MCQAB(:,:,MUE,ISPINPROJ,IQ),LMAT3)
               MCQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MCQBA(:,:,MUE,ISPINPROJ,IQ))
C
               MDQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQAB(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
               MDQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQBA(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
            END DO
         END DO
      END DO
C
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               DO I = 1,NKM
                  DO J = 1,NKM
                     MDQBAX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQBA_L(J,I,MUE,ISPINPROJ,IQ))
C
                     MDQABX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQAB_L(J,I,MUE,ISPINPROJ,IQ))
                  END DO
               END DO
C
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
c sig1_cpa, see explanations of the computing process in sig1_k
C-----------------------------------------------------------------------
      DO IQ = IQBOT_CHI,IQTOP_CHI
c
        DO ISPR = 1,NSPR
          ISPINPROJ = LIST_ISPR(ISPR)
c
          DO NUE = 1,3
          DO MUE = 1,3
            s_a=c0
            s_d=c0
            s_cb=c0
c
            call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                 c1,MCQBA_L(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                 tauqa(1,1,iq),nkmmax,c0,j1tauk(1,1),nkmmax)
            call zgemm('t','t',nkmmax,nkmmax,nkmmax,
     &                 c1,tauqbz(1,1,iq,2),nkmmax,
     &                 MBQAB_L(1,1,NUE,1,IQ),nkmmax,
     &                 c0,j2tauk_t(1,1),nkmmax)
            do L2=1,NKM
              s_cb=s_cb-zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
            enddo
c
            call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                 c1,MDQBAX(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                 tauqa(1,1,iq),nkmmax,c0,j1tauk(1,1),nkmmax)
            call zgemm('t','t',nkmmax,nkmmax,nkmmax,
     &                 c1,tauqbz(1,1,iq,1),nkmmax,
     &                 MDQABX(1,1,NUE,1,IQ),nkmmax,
     &                 c0,j2tauk_t(1,1),nkmmax)
            do L2=1,NKM
              S_D=S_D-zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
            enddo
c
            call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                 c1,MAQBA(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                 tauqa(1,1,iq),nkmmax,c0,j1tauk(1,1),nkmmax)
            call zgemm('t','t',nkmmax,nkmmax,nkmmax,
     &                 c1,tauqbz(1,1,iq,1),nkmmax,
     &                 MAQAB(1,1,NUE,1,IQ),nkmmax,
     &                 c0,j2tauk_t(1,1),nkmmax)
            do L2=1,NKM
              S_A=S_A-zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
            enddo
c
            SUM_A(MUE,NUE) = 0.25D0*S_A
            SUM_D(MUE,NUE) = 0.25D0*S_D
            SUM_CB(MUE,NUE) = 0.25D0*S_CB
          enddo
          enddo
          sig1_cpa(:,:,ispinproj,iq) = 2*sum_cb(:,:) 
     &      - dconjg(sum_d(:,:)) - sum_a(:,:)
        enddo
      enddo
C-----------------------------------------------------------------------
c allocate, check and mpi for sig1_k
C-----------------------------------------------------------------------
      ALLOCATE (IPIV(NKKR))
c      IF ( .NOT.ERYDA_EQ_ERYDB )
c     &     ALLOCATE (TAUKLINB(NKKR*NKKR))
      ALLOCATE (TAUKLINB(NKKR*NKKR))
      ALLOCATE (MAUX(NKKR,NKKR),TAUKLINA(NKKR*NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUKLINA')
c
      IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL = 2')
      IF ( NZ12.LT.1 .OR. NZ12.GT.2 )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12 not set properly')
      IF ( NZ12.NE.NZ12MAX )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12MAX not set properly')
      if(nz12.ne.2) stop 'temporarily force nz12 to be 2'
c
      level=toplevel+1
c
      call get_comm_level(level,'parent',lparalleled,
     &                    parent_comm,nprocs_parent,parent_rank)
      call mpi_barrier(parent_comm,mpierr)
      if(parent_rank==0) call cpu_time(t0_sig1)
c
      call mpi_multilevel_distribute(routine,level,nktab,
     &                               lparalleled,ik0,ik1)
c
      ALLOCATE(sig1q_tmp(3,3,NSPINPROJ,NQMAX,NQMAX))
      sig1q_tmp(:,:,:,:,:) = c0
c
      CALL STRCC(ERYDA,.FALSE.)
c
      if(parent_rank==0) write(*,*) 'wksum = ',wksum
      call mpi_barrier(parent_comm,mpierr)
C-----------------------------------------------------------------------
C tauk
C-----------------------------------------------------------------------
ckkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      DO IK = 1,NKTAB
         if((ik >= ik0 .and. ik <= ik1) .or. mpi_eloop) then
c tauka
         WK = WKTAB(IK)
         IF(WK .ne. 1.) then
           write(*,*) 'ik = ',ik,', wk = ',wk
           stop
         endif
         write(*,*) 'ik=',ik,',ik0=',ik0,',ik1=',ik1
         CWK = DCMPLX(WK,0D0)
         IF ( .NOT.ERYDA_EQ_ERYDB ) CALL STRCC(ERYDA,.FALSE.)
         CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLINA,PA)
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLINA,MSSQA,NQMAX,NKKR,NKMMAX)
         CALL ZGETRF(NKKR,NKKR,TAUKLINA,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUKLINA,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
c taukb
         IF ( .NOT.(ERYDA_EQ_ERYDB) ) THEN
            CALL STRCC(ERYDB,.FALSE.)
            CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLINB,PB)
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLINB,MSSQB,NQMAX,NKKR,
     &                  NKMMAX)
            CALL ZGETRF(NKKR,NKKR,TAUKLINB,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUKLINB,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
         ELSE
            TAUKLINB(:) = TAUKLINA(:)
         ENDIF
cqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
         DO JQ = IQBOT_CHI,IQTOP_CHI
           JQCHI = JQ - IQBOT_CHI + 1
         DO IQ = IQBOT_CHI,IQTOP_CHI
           IQCHI = IQ - IQBOT_CHI + 1
c
c store tauk as (L1,L2,IQ,JQ)
           jtt10 = (IND0Q(JQCHI)-1)*NKKR + IND0Q(IQCHI)
           jtt20 = (IND0Q(IQCHI)-1)*NKKR + IND0Q(JQCHI)
           DO L2 = 1,NKM
             JTT1 = jtt10 + L2*NKKR
             jtt2 = jtt20 + L2*NKKR
           DO L1 = 1,NKM
             JTT1 = JTT1 + 1
             jtt2 = jtt2 + 1
             tauka_ij(L1,L2)=tauklina(jtt1)
             taukx_ij(L1,L2) = dconjg(tauklinb(jtt1))
             taukb_ji(L1,L2) = tauklinb(jtt2)
           enddo ! l1
           enddo ! l2
C-----------------------------------------------------------------------
C sig1_k
C-----------------------------------------------------------------------
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
           DO ISPR = 1,NSPR
             ISPINPROJ = LIST_ISPR(ISPR)
c
           DO NUE = 1,3
           DO MUE = 1,3
c
             S_A = C0
             S_D = C0
             S_CB = C0
C
c sigCB = MCQBA_L(L4,L1,MUE,ISPINPROJ,IQ) * tauka(L1,L2,IQ,JQ) *
c         conjg( taukb(L4,L3,IQ,JQ) ) * MBQAB_L(L2,L3,NUE,1,JQ)
             call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                  c1,MCQBA_L(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                  tauka_ij(1,1),nkmmax,c0,j1tauk(1,1),nkmmax)
             call zgemm('n','t',nkmmax,nkmmax,nkmmax,
     &                  c1,taukx_ij(1,1),nkmmax,
     &                  MBQAB_L(1,1,NUE,1,JQ),nkmmax,
     &                  c0,j2tauk_t(1,1),nkmmax)
             do L2=1,NKM
               S_CB=S_CB + zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
             enddo
C
c sigDD = MDQBAX(L4,L1,MUE,ISPINPROJ,IQ) * tauka(L1,L2,IQ,JQ) *
c         taukb(L3,L4,JQ,IQ) * MDQABX(L2,L3,NUE,1,JQ)
             call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                  c1,MDQBAX(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                  tauka_ij(1,1),nkmmax,c0,j1tauk(1,1),nkmmax)
             call zgemm('t','t',nkmmax,nkmmax,nkmmax,
     &                  c1,taukb_ji(1,1),nkmmax,
     &                  MDQABX(1,1,NUE,1,JQ),nkmmax,
     &                  c0,j2tauk_t(1,1),nkmmax)
             do L2=1,NKM
               S_D=S_D + zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
             enddo
C
c sigAA = MAQBA(L4,L1,MUE,ISPINPROJ,IQ) * tauka(L1,L2,IQ,JQ) *
c         taukb(L3,L4,JQ,IQ) * MAQAB(L2,L3,NUE,1,JQ)
             call zgemm('n','n',nkmmax,nkmmax,nkmmax,
     &                  c1,MAQBA(1,1,MUE,ISPINPROJ,IQ),nkmmax,
     &                  tauka_ij(1,1),nkmmax,c0,j1tauk(1,1),nkmmax)
             call zgemm('t','t',nkmmax,nkmmax,nkmmax,
     &                  c1,taukb_ji(1,1),nkmmax,
     &                  MAQAB(1,1,NUE,1,JQ),nkmmax,
     &                  c0,j2tauk_t(1,1),nkmmax)
             do L2=1,NKM
               S_A=S_A + zdotu(nkmmax,j1tauk(1,L2),1,j2tauk_t(1,L2),1)
             enddo
C
             SUM_A(MUE,NUE) = 0.25D0*S_A
             SUM_D(MUE,NUE) = 0.25D0*S_D
             SUM_CB(MUE,NUE) = 0.25D0*S_CB
           END DO
           END DO ! mue, nue
C
           SUM_ALL = 2*SUM_CB - DCONJG(SUM_D) - SUM_A
C
           SIG1Q_TMP(1:3,1:3,ISPINPROJ,IQ,JQ) = 
     &       SIG1Q_TMP(1:3,1:3,ISPINPROJ,IQ,JQ) 
     &       + SUM_ALL(1:3,1:3) / wksum
C
           END DO ! ISPR
         END DO ! JQ
         END DO ! IQ
      endif ! ik
      enddo ! ik
c
c collect sig1q
      if(.not. lparalleled) then
        sig1q(:,:,:,:,:) = sig1q_tmp(:,:,:,:,:)
      else
        call mpi_barrier(parent_comm,mpierr)
        call get_comm_level(level,'inter ',lparalleled,
     &                      inter_comm,nprocs_inter,inter_rank)
c
        call mpi_reduce(sig1q_tmp,sig1q,3*3*NSPINPROJ*NQMAX*NQMAX,
     &                  mpi_double_complex,mpi_sum,0,inter_comm,mpierr)
        call mpi_barrier(parent_comm,mpierr)
c free mpi_distribution on ik
        call mpi_multilevel_free(level)
      endif
      if(parent_rank==0) then
        call cpu_time(t1_sig1)
        write(*,*)
        write(*,*) 'sig1_surf_jtauksquare takes ',
     &              t1_sig1-t0_sig1,' seconds'
        write(*,*)
      endif
      deallocate(sig1q_tmp)
C-----------------------------------------------------------------------
c sig1 = sig1_cpa + sum_k sig1_k
C-----------------------------------------------------------------------
      DO IQ = IQBOT_CHI,IQTOP_CHI
        sig1q(:,:,:,iq,iq) = sig1q(:,:,:,iq,iq) + sig1_cpa(:,:,:,iq)
      enddo
ckkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ, JQ = ',2I3,/,10X,21('='),/)
99002 FORMAT (/,10X,'sigma 1 ',5X,A,/)
99003 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2)
99004 FORMAT (3(F14.6,F12.6))
C99005 FORMAT (//,1X,79('*'),/,34X,'<SIG1_SURF>',/,1X,79('*'),//,10X,A,/)
99005 FORMAT (//,1X,79('*'),/,34X,'<',A9,'>',/,1X,79('*'),//,10X,A,/)
99006 FORMAT (/,10X,'(z1,z2) = ',2I2,/)
99007 FORMAT (10X,3F14.6)
99008 FORMAT (/,'    SI units [',A,']')
99009 FORMAT (/,40('*'),/,40('*'))
99010 FORMAT (3('(',F14.6,',',F12.6,')'))
      END
