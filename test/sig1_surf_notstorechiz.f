C*==sig1_surf.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_SURF_NOTSTORECHIZ(
     &             ERYDA,PA,ERYDB,PB,TAUQA,TAUQB,TAUQBZ,MSSQA,MSSQB,
     &             MAQAB,MBQAB,MCQAB,MDQAB,SIG1Q,KEY,MAQBA,
     &             MCQBA,MBQBA,MDQBA,NSPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,NZ12,NZ12MAX
      USE MOD_KSPACE,ONLY:WKTAB,NKTAB,KTAB,WKSUM
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,LMAT3,IND0Q,NKMQ,NKKR,WKM1
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SOTSI,EESI,
     &    IRESPONSE_SOT,IRESPONSE_EDELSTEIN
      USE MOD_CONSTANTS,ONLY:C0,CI
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
      PARAMETER (ROUTINE='SIG1_SURF_NOTSTORECHIZ')
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
      REAL*8 PRESI
      CHARACTER*40 SIUNITS
      complex*16 ctmp1, ctmp2, ctmp3
c
      logical lparalleled, lopen
      integer level, ik0, ik1, nkp, mpierr,
     &        parent_comm, nprocs_parent, parent_rank,
     &        inter_comm,nprocs_inter,inter_rank,
     &        mpi_row
      integer ik
      real*8 t0_sig1, t1_sig1, wk
      complex*16, dimension(:,:,:,:,:), allocatable :: sig1q_tmp
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
      IF ( KEY.EQ.'N' ) WRITE (6,99005) ROUTINE,
     &                                  'without vertex-corrections'
      IF ( KEY.EQ.'V' ) WRITE (6,99005) ROUTINE,
     &                                  'including vertex-corrections'
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
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
ckkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
      ALLOCATE (IPIV(NKKR))
      ALLOCATE (CHIZ(NKM**2,NKM**2,NZ12))
c      IF ( .NOT.ERYDA_EQ_ERYDB )
c     &     ALLOCATE (TAUKLINB(NKKR*NKKR))
      ALLOCATE (TAUKLINB(NKKR*NKKR))
      ALLOCATE (MAUX(NKKR,NKKR),TAUKLINA(NKKR*NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUKLINA')
C
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
c
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
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
         IQCHI = IQ - IQBOT_CHI + 1
c
c         if(parent_rank==0) write(*,*) 'iq = ',iq
C
         DO JQ = IQBOT_CHI,IQTOP_CHI
         JQCHI = JQ - IQBOT_CHI + 1
c debug
            if(mpi_id==0) write(*,*) 'iq= ',iq,'jq= ',jq
c debug
C
            jtt10 = (IND0Q(JQCHI)-1)*NKKR + IND0Q(IQCHI)
            jtt20 = (IND0Q(IQCHI)-1)*NKKR + IND0Q(JQCHI)
            DO L4 = 1,NKM
               JTT2 = jtt20 + L4*NKKR
               JTTX = jtt10 + L4
               DO L3 = 1,NKM
                  JTT2 = JTT2 + 1
                  JTTX = JTTX + NKKR
                  K23J = L3 - NKM
C
                  TAUK2 = TAUKLINB(JTT2)
                  TAUKX = DCONJG(TAUKLINB(JTTX))
C
                  DO L2 = 1,NKM
                     JTT1 = jtt10 + L2*NKKR
                     K23J = K23J + NKM
                     K14I = L4 - NKM
                     DO L1 = 1,NKM
                        JTT1 = JTT1 + 1
                        K14I = K14I + NKM
C
                        TAUK1 = TAUKLINA(JTT1)
c
                        if(iq==jq) then
                          CHIZ(K14I,K23J,1) = TAUK1*TAUK2 
     &                      - TAUQA(L1,L2,IQ)*TAUQBZ(L3,L4,IQ,1)
                          CHIZ(K14I,K23J,2) = TAUK1*TAUKX 
     &                      - TAUQA(L1,L2,IQ)*TAUQBZ(L3,L4,IQ,2)
                        else
                          CHIZ(K14I,K23J,1) = TAUK1*TAUK2
                          CHIZ(K14I,K23J,2) = TAUK1*TAUKX
                        endif
C
                     END DO
                  END DO
               END DO
            END DO
c
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
            DO ISPR = 1,NSPR
               ISPINPROJ = LIST_ISPR(ISPR)
C
               IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
C---------- torkance prefactor for SI output (C * m)
                  PRESI = SOTSI
                  SIUNITS = '10**-30 C*m'
               ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
C---------- Edelstein prefactor for SI output (m/V)
                  PRESI = EESI
                  SIUNITS = 'm/V'
               ELSE
C---------- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
                  PRESI = 1D-8*CONSI
                  SIUNITS = '1/(muOhm*cm)'
               END IF
C
c               WRITE (6,99009)
c               WRITE (6,'(/,A80)') STR_ISP_PROJ(ISPINPROJ)
c               WRITE (6,99009)
C
               DO NUE = 1,3
C
                  DO MUE = 1,3
C
                     S_A = C0
                     S_D = C0
                     S_CB = C0
C
                     k2 = 0
c                     K2 = (JQCHI-1)*NKMSQ
                     DO L2 = 1,NKM
                        DO L3 = 1,NKM
                           K2 = K2 + 1
c
                           ctmp1 = MBQAB_L(L2,L3,NUE,1,JQ)
                           ctmp2 = MDQABX(L2,L3,NUE,1,JQ)
                           ctmp3 = MAQAB(L2,L3,NUE,1,JQ)
c
c                           K1 = (IQCHI-1)*NKMSQ
                           k1 = 0
                           DO L1 = 1,NKM
                              DO L4 = 1,NKM
                                 K1 = K1 + 1
C
                                 S_CB = S_CB + 
     &                                  MCQBA_L(L4,L1,MUE,ISPINPROJ,IQ)
     &                                  *CHIZ(K1,K2,2)*ctmp1
C
                                 S_D = S_D + MDQBAX(L4,L1,MUE,ISPINPROJ,
     &                                 IQ)*CHIZ(K1,K2,1)*ctmp2
C
                                 S_A = S_A + MAQBA(L4,L1,MUE,ISPINPROJ,
     &                                 IQ)*CHIZ(K1,K2,1)*ctmp3
C
                              END DO
C
                           END DO
                        END DO
                     END DO
C
                     SUM_A(MUE,NUE) = 0.25D0*S_A
                     SUM_D(MUE,NUE) = 0.25D0*S_D
                     SUM_CB(MUE,NUE) = 0.25D0*S_CB
C
                  END DO
               END DO
C
               SUM_ALL = 2*SUM_CB - DCONJG(SUM_D) - SUM_A
c debug
               if(mpi_id==0) then
                 do nue=1,3
                   write(*,'(a,6e16.8)') 'sum_cb = ',sum_cb(1:3,nue)
                 enddo
                 do nue=1,3
                   write(*,'(a,6e16.8)') 'sum_d  = ',sum_d (1:3,nue)
                 enddo
                 do nue=1,3
                   write(*,'(a,6e16.8)') 'sum_a  = ',sum_a (1:3,nue)
                 enddo
               endif
c debug
C
C-----------------------------------------------------------------------
               SIG1Q_TMP(1:3,1:3,ISPINPROJ,IQ,JQ) = 
     &           SIG1Q_TMP(1:3,1:3,ISPINPROJ,IQ,JQ) 
     &           + SUM_ALL(1:3,1:3) / wksum
C
C-----------------------------------------------------------------------
C
            END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         END DO ! JQ
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      endif ! ik
      enddo ! ik
c
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
        write(*,*) 'sig1_surf_notstorechiz takes ',
     &              t1_sig1-t0_sig1,' seconds'
        write(*,*)
      endif
      deallocate(sig1q_tmp)
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
