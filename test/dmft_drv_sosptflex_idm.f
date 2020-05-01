C*==dmft_drv_sosptflex_idm.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
      SUBROUTINE DMFT_DRV_SOSPTFLEX_IDM(GFMAT,DMFTSIGMA,NKM,NEGF,LOP,
     &                                  UEFF,JEFF,UDMFT,DMFTTEMP,
     &                                  DMFTDBLC,NOM,IOTMP,IPRINTDMFT,
     &                                  IFLEX,TXT_T,LTXT_T,ITRSCF,
     &                                  DMFTMIX,NKMMAX,NLMAX,WEGFOLD,
     &                                  EGFOLD,EGFNEW,EREFDMFT,IT,NTMAX,
     &                                  TMA_STATIC,TMA_NOBATH,TMA_DBLC,
     &                                  DOSINT,TMA_NOCCSCL,LDAUSIGMA)
C
C
C
      USE GREEN_SPTF_INTERFACE,ONLY:SPTF
      USE MOD_MPI,ONLY:MPI_ID
      USE MOD_DMFT_LDAU,ONLY:DMFT_FIX_DYN_SE, zerosig_contn
      use mod_energy,only:igrid,nefd1,nefd2,nefd3,nepol,necontn,
     &                    lactive_contn
      use mod_mpi_multilevels
      IMPLICIT NONE
C*--DMFT_DRV_SOSPTFLEX_IDM18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EVOVRY,PI,TOL
      PARAMETER (EVOVRY=13.605826D0,PI=3.141592653589793238462643D0,
     &           TOL=1.0D-10)
      COMPLEX*16 C0,C1,CI
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
C
C Dummy arguments
C
      CHARACTER*4 DMFTDBLC
      REAL*8 DMFTMIX,DMFTTEMP,JEFF,UEFF
      COMPLEX*16 EREFDMFT
      INTEGER IFLEX,IOTMP,IPRINTDMFT,IT,ITRSCF,LOP,LTXT_T,NEGF,NKM,
     &        NKMMAX,NLMAX,NOM,NTMAX
      CHARACTER*3 TMA_DBLC
      LOGICAL TMA_NOBATH,TMA_NOCCSCL,TMA_STATIC
      CHARACTER*8 TXT_T
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NEGF),EGFNEW(NEGF),EGFOLD(NEGF)
     &           ,GFMAT(NKM,NKM,NEGF),LDAUSIGMA(NKMMAX,NKMMAX),
     &           WEGFOLD(NEGF)
      REAL*8 DOSINT(NLMAX),UDMFT(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX)
C
C Local variables
C
c modified by XJQ: replace first continuation by an integral
c      real*8 reegf, imegf, aegf
c      complex*16 kernalomz
c end-mod-xjq
      REAL*8 BETA,DE,DOMEGA,DTIME,EFERMINEW,EFERMIOLD,OCC(:,:),OMEGA(:),
     &       RWORK,TEMP,TIME(:),TMP,UJJ,UR(:,:,:,:),UUU
      COMPLEX*16 CSUM,CWORK(:),CWORK1(:,:),CWORK2(:,:),CWORKO(:),
     &           DCTMP(:,:),EMIGD,EREAL(:),GOM(:,:,:),GOMIN(:,:,:),
     &           GOMWORK(:,:,:),GZ(:,:,:),HAMIL(:,:),OVLP(:,:),
     &           SEREAL(:,:,:),SEZ(:,:,:,:),SEZNEW(:,:,:),SIGSTAT(:,:),
     &           SIG_HF(:,:),
     &           SOM(:,:,:),SOMSAV(:,:,:,:),SOMWORK(:,:,:),
     &           UCREP(:,:,:,:)
      CHARACTER*80 FILNAM
      INTEGER I1,I2,I3,I4,ICALL,IE,IOM,IR,IRP,IS1,IS2,IZERO(:,:),
     &        LFILNAM,LM,LMP,LMPP,LMPPP,MINOM(:),MIOM,NEREAL,NLM,NOMIN,
     &        NOMLOC,NREP
      LOGICAL VERBOSE
c modified by XJQ: mpi_multilevels
      integer ie0_contn, ie1_contn
      logical lpositive
      logical lparalleled
      integer level, parent_comm, nprocs_parent, parent_rank, mpierr
      complex*16, allocatable :: epath_tmp(:)
c end-mod-xjq
      SAVE IZERO,SEZ,SOMSAV
C
C*** End of declarations rewritten by SPAG
C
C
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE GZ,CWORK1,CWORK2,SEZ,SOM,SEZNEW,SOMSAV
      ALLOCATABLE OMEGA,TIME,MINOM,IZERO,GOM,CWORK,CWORKO,GOMIN
c      ALLOCATABLE cwork_egf_reord, cwork_gz_reord
      ALLOCATABLE SEREAL,EREAL,UR,GOMWORK,SOMWORK,UCREP
      ALLOCATABLE OCC,SIGSTAT
C
      ALLOCATABLE SIG_HF,DCTMP,HAMIL,OVLP
C
      NEREAL = 301
C
C
      NLM = LOP*2 + 1
      NREP = 2*NLM
C
      TEMP = DMFTTEMP*0.6333659D-5
      UUU = UEFF/EVOVRY
      UJJ = JEFF/EVOVRY
C
      EFERMIOLD = DREAL(EGFOLD(NEGF))
      EFERMINEW = DREAL(EGFNEW(NEGF))
C
      IF ( IPRINTDMFT.EQ.1 ) WRITE (6,*) TMA_DBLC
      IF ( IFLEX.EQ.0 ) WRITE (6,99002) 
     &                'Only TMA-approximation used (PH-channel ignored)'
C
C
C
      IF ( ICALL.EQ.0 ) then
        if(.not.allocated(sez)) ALLOCATE (SEZ(NREP,NREP,NEGF,NTMAX))
      endif
C
      ALLOCATE (GZ(NREP,NREP,NEGF))
      I1 = LOP**2 + 1
      I2 = I1 + NLM - 1
      I3 = NKM/2 + I1
      I4 = I3 + NLM - 1
C
      SEZ(1:NLM,1:NLM,1:NEGF,IT) = DMFTSIGMA(I1:I2,I1:I2,1:NEGF)
      SEZ(NLM+1:NREP,NLM+1:NREP,1:NEGF,IT)
     &   = DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
      SEZ(1:NLM,NLM+1:NREP,1:NEGF,IT) = DMFTSIGMA(I1:I2,I3:I4,1:NEGF)
      SEZ(NLM+1:NREP,1:NLM,1:NEGF,IT) = DMFTSIGMA(I3:I4,I1:I2,1:NEGF)
C
      GZ(1:NLM,1:NLM,1:NEGF) = GFMAT(I1:I2,I1:I2,1:NEGF)         !1_1
      GZ(NLM+1:NREP,NLM+1:NREP,1:NEGF) = GFMAT(I3:I4,I3:I4,1:NEGF)
                                                                 !2_2
      GZ(1:NLM,NLM+1:NREP,1:NEGF) = GFMAT(I1:I2,I3:I4,1:NEGF)    !1_2
      GZ(NLM+1:NREP,1:NLM,1:NEGF) = GFMAT(I3:I4,I1:I2,1:NEGF)    !2_1
C
C
C============================================================
C============================================================
C     Check which elements of GZ are non_zero
C============================================================
      IF ( .NOT.ALLOCATED(IZERO) ) THEN
         ALLOCATE (IZERO(NREP,NREP))
         IZERO(1:NREP,1:NREP) = 0
         DO LM = 1,NREP
            DO LMP = 1,NREP
               CSUM = C0
               DO IE = 1,NEGF
                  CSUM = CSUM + GZ(LMP,LM,IE)
               END DO
               IF ( ABS(CSUM).GE.TOL ) IZERO(LMP,LM) = 1
            END DO
         END DO
      END IF
C
      DO IR = 1,NREP
         DO IRP = 1,NREP
            IF ( IZERO(IR,IRP).EQ.0 ) GZ(IR,IRP,1:NEGF) = C0
         END DO
      END DO
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C Calculate occupation numbers:
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      ALLOCATE (OCC(NREP,NREP))
      TMP = 0D0
      OCC(1:NREP,1:NREP) = 0D0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            IF ( IZERO(IR,IRP).EQ.1 ) THEN
               DO IE = 1,NEGF
                  OCC(IR,IRP) = OCC(IR,IRP)
     &                          - DIMAG(GZ(IR,IRP,IE)*WEGFOLD(IE))
               END DO
            END IF
         END DO
         TMP = TMP + OCC(IR,IR)
      END DO
      OCC(1:NREP,:NREP) = OCC(1:NREP,1:NREP)/PI
      TMP = TMP/PI
C
      IF ( IPRINTDMFT.EQ.1 ) THEN
         CALL RMATSTRUCT('occupation numbers:',OCC,NREP,NREP,0,0,0,TOL,
     &                   6)
         IF ( MPI_ID.EQ.0 ) THEN
            WRITE (*,'(/,10X,A,F8.4)') 'Trace:',TMP
         ENDIF
      END IF
C
C============================================================
C     G_BATH=G^-1+SIGMA
C============================================================
C
C
      ALLOCATE (CWORK1(NREP,NREP),CWORK2(NREP,NREP))
      allocate (cwork(nrep))
c
      if(lactive_contn) then
        ie0_contn = nefd1+nefd2+nefd3+1
        ie1_contn = negf-nepol
      else
        ie0_contn = 1
        ie1_contn = negf
      endif
C
      IF ( .NOT.TMA_NOBATH ) THEN
         DO IE = 1,NEGF
c            CWORK1(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,IE)
c            CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
c            CWORK1(1:NREP,1:NREP) = CWORK2(1:NREP,1:NREP)
c     &                              + SEZ(1:NREP,1:NREP,IE,IT)
c            CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
c            GZ(1:NREP,1:NREP,IE) = CWORK2(1:NREP,1:NREP)
            if( zerosig_contn .and. lactive_contn .and.
     &          ie>=ie0_contn .and. ie<=ie1_contn ) cycle
            CWORK2(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,IE)
            call cmatinv2(nrep,nrep,cwork,cwork1,cwork2)
            cwork1(1:NREP,1:NREP) = CWORK2(1:NREP,1:NREP)
     &                              + SEZ(1:NREP,1:NREP,IE,IT)
            call cmatinv2(nrep,nrep,cwork,cwork2,cwork1)
            GZ(1:NREP,1:NREP,IE) = CWORK1(1:NREP,1:NREP)
         END DO
      END IF
      deallocate (cwork)
      DEALLOCATE (CWORK1,CWORK2)
C
      IF ( IPRINTDMFT.GT.0 ) THEN
         DO LM = 1,NREP
            WRITE (6,99001) (IZERO(LM,LMP),LMP=1,NREP)
         END DO
      END IF
      IF ( TMA_NOCCSCL ) THEN
         TMP = DOSINT(LOP+1)/TMP
         DO IR = 1,NREP
            OCC(IR,IR) = OCC(IR,IR)*TMP
            GZ(IR,IR,:) = GZ(IR,IR,:)*TMP
         END DO
         TMP = 0.0D0
         DO IR = 1,NREP
            TMP = TMP + OCC(IR,IR)
         END DO
         WRITE (*,'(/,20X,A,F8.4)') 'Trace: RESCALED',TMP
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C
C============================================================
C     Create matsubara mesh
C============================================================
C    Be carefull, SOC-SPTFLEX (Leonid Pourovski) is using
C     different numbering of  matsubara poles (using only positive
C     matsubara) as SPTFLEX  (Lulu Chioncel).  That's why
C     NOMLOC=NOM*2 in here to have in SOC-SPTFLEX NOM as
C     given in input
C
      NOMLOC = 2*NOM
C
      ALLOCATE (OMEGA(NOMLOC),TIME(NOMLOC),MINOM(NOMLOC))
      BETA = 1.0D0/TEMP
      DOMEGA = TEMP*PI
      DTIME = 1.D0/NOMLOC/TEMP*2.D0
      CALL DMFT_MATSUB(DTIME,DOMEGA,TIME,OMEGA,MINOM,NOMLOC)
C
C============================================================
C     Create Green's function on matsubara (pade approx)
C============================================================
      ALLOCATE (GOM(NREP,NREP,NOMLOC))
      ALLOCATE (GOMIN(NREP,NREP,NOMLOC))
      ALLOCATE (CWORK(NEGF))
      ALLOCATE (CWORKO(NOMLOC))
      GOM(1:NREP,1:NREP,1:NOMLOC) = C0
      GOMIN(1:NREP,1:NREP,1:NOMLOC) = C0
      CWORK(1:NEGF) = C0
      CWORKO(1:NOMLOC) = C0
C
      DO LM = 1,NREP
         DO LMP = 1,NREP
            IF ( IZERO(LMP,LM).EQ.1 ) THEN
c
c modified by XJQ: more stable pade using code of Johan Schott and Elin Lundin,
c                  their subroutine pade is paralleled so that speed is not slow,
c                  subroutine pade_averaging_interface is an interface for sprkkr
c
c instead of Pade continuation, use
c g(iwn) = int dw a(w) / ( iwn - w ) - 2pi*ci*res. 
c this way failed in a trial. does it work when w is complex ? should a(w) be normalized ?
c
c                 cworko(1:nomloc) = c0
c                 do iom=1,nomloc/2
c                   do ie=1,negf-nepol
c                     reegf = real(egfold(ie)) - EFERMIOLD
c                     imegf = aimag(egfold(ie))
c                     kernalomz = - ( ci*( omega(2*iom)-imegf ) + reegf )
c     &                 / ( ( omega(2*iom)-imegf )**2d0 + reegf**2d0 )
c                     aegf = -aimag(gz(lmp,lm,ie) * we_dmft_aux(ie)) /pi
c                     cworko(2*iom) = cworko(2*iom) + aegf * kernalomz
c                   enddo
c                   if(iom .le. nepol) then
c                     aegf = - aimag(gz(lmp,lm,negf-iom+1)) / pi
c                     cworko(2*iom) = cworko(2*iom) - 2.0d0*pi*ci*aegf
c                   endif
c                 enddo
c
               CWORK(1:NEGF) = GZ(LMP,LM,1:NEGF)
c               CALL DMFT_PADEZO(EFERMIOLD,NEGF,EGFOLD,CWORK,
c     &                          OMEGA,MINOM,CWORKO,NOMLOC)
               call pade_averaging_interface(egfold,ie0_contn,ie1_contn,
     &           cwork,efermiold,temp,nomloc,cworko,1,it,lm,lmp,0,iotmp)
c               if(.not.lactive_contn) then
c               call pade_averaging_interface(egfold,1,negf,
c     &           cwork,efermiold,temp,nomloc,cworko,1,it,lm,lmp,0,iotmp)
c               else
c               call pade_averaging_interface(egfold,nefd1+nefd2+nefd3+1,
c     &           negf-nepol,cwork,efermiold,temp,nomloc,cworko,
c     &           1,it,lm,lmp,0,iotmp)
c               endif
C
               do iom=1,nomloc/2
                 GOM(LMP,LM,2*IOM) = CWORKO(2*IOM)
               enddo
c the following copy part may be not needed, it can be closed by "if(.false.) then"
               if(.false.) then
               if(igrid(1)==11) then
                 do iom=1,nepol
                   if(mpi_id==0 .and. lm==lmp) then
                     write(6,'(a3,i3,a4,i3,a4,2e16.8)')
     &                 'lm=',lm,'iom=',iom,
     &                 'err=',(CWORKO(2*IOM)-GZ(LM,LM,negf-iom+1))
                   endif
                   CWORKO(2*IOM) = GZ(LMP,LM,negf-iom+1)
                   GOM(LMP,LM,2*IOM) = CWORKO(2*IOM)
                 enddo
               else
               DO IOM = 1,NOMLOC/2
                  DO IE = 1,NEGF
C    IF SOME KKR ENERGY POINTS DONE ON  MATSUBARA JUST TAKE THEM OVER
                     IF ( ABS(EFERMIOLD-DREAL(EGFOLD(IE))).LE.1D-8 .AND. 
     &                    ABS(OMEGA(2*IOM)-DIMAG(EGFOLD(IE))).LT.1D-8 ) 
     &               THEN
                       if(mpi_id==0) then
                         write(6,'(a5,i3,a6,i3,a6,2e16.8)') 
     &                     'lm = ',lm,'iom = ',iom,
     &                     'err = ',(CWORKO(2*IOM)-GZ(LMP,LM,IE))
                       endif 
                       CWORKO(2*IOM) = GZ(LMP,LM,IE)
                       GOM(LMP,LM,2*IOM) = CWORKO(2*IOM)
                     ENDIF
                  END DO
               END DO
               endif
               endif
c
c end-mod-xjq
c
            END IF
         END DO
      END DO
      GOMIN(1:NREP,1:NREP,1:NOMLOC) = GOM(1:NREP,1:NREP,1:NOMLOC)
      DEALLOCATE (CWORK)
      DEALLOCATE (CWORKO)
c modified by XJQ: mpi_multilevels
       level=toplevel+1
       call get_comm_level(level,'parent',lparalleled,
     &                     parent_comm,nprocs_parent,parent_rank)
c end-mod-xjq
      IF ( parent_rank.EQ.0 .AND. IPRINTDMFT.GT.0 ) THEN
c
        CALL DMFT_CLOSE_IOTMP(IOTMP)
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//'_gz_im.dat')
        CLOSE (IOTMP,STATUS='DELETE')
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//'_gz_im.dat')
        DO IE=1,NEGF
          WRITE (IOTMP,'(40f20.10)') EGFOLD(IE),
     &          (DIMAG(gz(LM,LM,IE)),LM=1,NREP)
        END DO
        CALL DMFT_CLOSE_IOTMP(IOTMP)
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//'_gz_re.dat')
        CLOSE (IOTMP,STATUS='DELETE')
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//'_gz_re.dat')
        DO IE=1,NEGF
          WRITE (IOTMP,'(40f20.10)') EGFOLD(IE),
     &          (DREAL(gz(LM,LM,IE)),LM=1,NREP)
        END DO
c
        CALL DMFT_CLOSE_IOTMP(IOTMP)
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//
     &                        '_gomin_im.dat')
        CLOSE (IOTMP,STATUS='DELETE')
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//
     &                        '_gomin_im.dat')
        DO IOM = 2,NOMLOC/2,2
          WRITE (IOTMP,'(40f20.10)') OMEGA(IOM),
     &          (DIMAG(GOMIN(LM,LM,IOM)),LM=1,NREP)
        END DO
        CALL DMFT_CLOSE_IOTMP(IOTMP)
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//
     &                        '_gomin_re.dat')
        CLOSE (IOTMP,STATUS='DELETE')
        OPEN (UNIT=IOTMP,FILE='diag-gf/'//TXT_T(1:LTXT_T)//
     &                         '_gomin_re.dat')
        DO IOM = 2,NOMLOC/2,2
          WRITE (IOTMP,'(40f20.10)') OMEGA(IOM),
     &          (DREAL(GOMIN(LM,LM,IOM)),LM=1,NREP)
        END DO
        CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C
C
C=================================================================
C     Calculate self-energy  on matsubara frequencies
C=================================================================
C
      ALLOCATE (SOM(NREP,NREP,NOMLOC))
      SOM(1:NREP,1:NREP,1:NOMLOC) = C0
      IF ( ICALL.EQ.0 ) THEN
         ALLOCATE (SOMSAV(NREP,NREP,NOMLOC,NTMAX))
         SOMSAV = C0
      END IF
C
C============================================================
C     Setup 4-index U matrix
C============================================================
      ALLOCATE (UR(NLM,NLM,NLM,NLM))
      UR(1:NLM,1:NLM,1:NLM,1:NLM) = 0.0D0
      ALLOCATE (UCREP(NREP,NREP,NREP,NREP))
      UCREP(1:NREP,1:NREP,1:NREP,1:NREP) = C0
      IF ( IPRINTDMFT.GT.0 ) WRITE (6,*) 
     &                            'Calculation of the 4-index U matrix:'
C
      DO LM = 1,NLM
         DO LMP = 1,NLM
            DO LMPP = 1,NLM
               DO LMPPP = 1,NLM
                  UR(LMPPP,LMPP,LMP,LM) = UDMFT(LMPPP,LMP,LMPP,LM)
               END DO
            END DO
         END DO
      END DO
C
      DO IS1 = 1,2
         DO IS2 = 1,2
            UCREP((IS1-1)*NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM,(IS1-1)
     &         *NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM)
     &         = UR(1:NLM,1:NLM,1:NLM,1:NLM)
         END DO
      END DO
C
C
C
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C Calculate static shifts /LSDA+U: Around Mean-Field part:
C  Czyzik and Zawatzky, PRB49(1994)14211: formula (4)/:
C TMP is taken as orbital average for each spin:
C
      ALLOCATE (SIGSTAT(NREP,NREP))
      SIGSTAT(1:NREP,1:NREP) = C0
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
         I1 = LOP**2 + 1
         I2 = I1 + NLM - 1
         I3 = NKM/2 + I1
         I4 = I3 + NLM - 1
         SIGSTAT(1:NLM,1:NLM) = (LDAUSIGMA(I1:I2,I1:I2))
         SIGSTAT(NLM+1:NREP,NLM+1:NREP) = (LDAUSIGMA(I3:I4,I3:I4))
         SIGSTAT(1:NLM,NLM+1:NREP) = (LDAUSIGMA(I1:I2,I3:I4))
         SIGSTAT(NLM+1:NREP,1:NLM) = (LDAUSIGMA(I3:I4,I1:I2))
         IF ( IPRINTDMFT.GT.0 ) CALL CMATSTRUCT('NEWtic self-energy:',
     &        SIGSTAT,NREP,NREP,0,0,0,TOL,6)
      END IF
C DMFTDBLC='LDAU'
C
C
C
      IF ( DMFT_FIX_DYN_SE ) TMA_STATIC = .TRUE.
      IF ( .NOT.TMA_STATIC ) THEN
C
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
         IF ( IPRINTDMFT.GT.0 ) THEN
            WRITE (6,*) 'Ueff (Ry)=',UUU,'Jeff (Ry)=',UJJ
            WRITE (6,*) 'Uij ='
            WRITE (6,'(5f9.4)') ((UR(LM,LMP,LM,LMP),LM=1,NLM),LMP=1,NLM)
         END IF
         NOMIN = 0
         DO IOM = 2,NOMLOC/2,2
            NOMIN = NOMIN + 1
         END DO
C
         ALLOCATE (GOMWORK(NREP,NREP,1:NOMIN))
         ALLOCATE (SOMWORK(NREP,NREP,1:NOMIN))
C
         GOMWORK(1:NREP,1:NREP,1:NOMIN) = C0
         SOMWORK(1:NREP,1:NREP,1:NOMIN) = C0
C
         NOMIN = 0
         DO IOM = 2,NOMLOC/2,2
            NOMIN = NOMIN + 1
            GOMWORK(1:NREP,1:NREP,NOMIN) = GOM(1:NREP,1:NREP,IOM)
         END DO
C
C         OPEN (UNIT=47,FILE='gf.dat')
C         WRITE (47,*) 'U=(Ry)',UUU,'J=(Ry)',UJJ
C         WRITE (47,*) 'TEMP=',TEMP
C         WRITE (47,*) 'NOM=',NOMIN,'NREP=',NREP
C         DO IOM = 1,NOMIN
C            WRITE (47,'(I4,200f16.10)') IOM,
C     &                                  ((GOMWORK(LM,LMP,IOM),LM=1,NREP)
C     &                                  ,LMP=1,NREP)
C         END DO
C         CLOSE (47)
C
         ALLOCATE (SIG_HF(NREP,NREP),DCTMP(NREP,NREP),HAMIL(NREP,NREP),
     &             OVLP(NREP,NREP))
C
         SIG_HF(:,:) = C0
         DCTMP(:,:) = C0
         HAMIL(:,:) = C0
         OVLP(:,:) = C0
         DO LM = 1,NREP
            OVLP(LM,LM) = C1
         END DO
         VERBOSE = .FALSE.
         IF ( IPRINTDMFT.EQ.1 ) VERBOSE = .TRUE.
C
         CALL SPTF(GOMWORK,SOMWORK,SIG_HF,UCREP,BETA,DCTMP,HAMIL,OVLP,
     &             VERBOSE)
C
C         OPEN (UNIT=14,FILE='som.dat')
C         DO IOM = 1,NOMIN
C            WRITE (14,'(I4,200f16.10)') IOM,
C     &           ((SOMWORK(LM,LMP,IOM),LM=1,NREP)
C     &           ,LMP=1,NREP)
C         END DO
C
C         CALL FLUSH(14)
C         CLOSE (14)
         NOMIN = 0
         DO IOM = 2,NOMLOC/2,2
            NOMIN = NOMIN + 1
            SOM(1:NREP,1:NREP,IOM) = SOMWORK(1:NREP,1:NREP,NOMIN)
         END DO
C
         DO IOM = 1,NOMLOC/2
            MIOM = NOMLOC/2 + IOM
            SOM(1:NREP,1:NREP,MIOM) = DCONJG(SOM(1:NREP,1:NREP,IOM))
         END DO
C
         IF ( IPRINTDMFT.GT.0 .AND. MPI_ID.EQ.0 )
     &         CALL CMATSTRUCT('Sigma at first omega',SOM(1,1,2),NREP,
     &        NREP,0,0,0,1D-10,6)
C
         DO LM = 1,NREP
            DO LMP = 1,NREP
               IF ( IZERO(LM,LMP).EQ.0 ) SOM(LM,LMP,:) = C0
            END DO
         END DO
         SOMSAV(:,:,:,IT) = SOM(:,:,:)
C  TMA_STATIC
      ELSE
         SOM(:,:,:) = C0
         IF ( DMFT_FIX_DYN_SE ) SOM(:,:,:) = SOMSAV(:,:,:,IT)
      END IF
C
C double-counting::
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
         DO LM = 1,NREP
            IF ( IZERO(LM,LM).NE.0 ) TMP = DREAL(SOM(LM,LM,2))
            SOM(LM,LM,1:NOMLOC) = SOM(LM,LM,1:NOMLOC) - TMP
         END DO
         DO IOM = 1,NOMLOC/2
            SOM(1:NREP,1:NREP,2*IOM) = SOM(1:NREP,1:NREP,2*IOM)
     &                                 + SIGSTAT(1:NREP,1:NREP)
         END DO
      END IF
C
C         IF (.NOT. FIXSE) THEN
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Set up the total energy correction, Galitskii-Migdal formula
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EMIGD = C0
      DO LM = 1,NREP
         DO LMP = 1,NREP
C            IF ( LM.EQ.LMP ) THEN
            DO IOM = 1,NOMLOC/2
               EMIGD = EMIGD + TEMP*SOM(LM,LMP,2*IOM)
     &                 *GOMIN(LMP,LMP,2*IOM)/2.0D0
            END DO
C            END IF
         END DO
      END DO
C
      IF ( MPI_ID.EQ.0 ) THEN
         IF ( IPRINTDMFT.GT.0 ) THEN
            WRITE (6,99003)
            WRITE (6,99005) DREAL(EMIGD),DIMAG(EMIGD)
            WRITE (6,99004)
         END IF
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     print-out of the self energy
C
      IF ( parent_rank.EQ.0 ) THEN
         IF ( IPRINTDMFT.GT.0 ) THEN
            CALL DMFT_CLOSE_IOTMP(IOTMP)
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &                            '_sig_matsub_im.dat')
            CLOSE (IOTMP,STATUS='DELETE')
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &                            '_sig_matsub_im.dat')
            DO IOM = 2,NOMLOC/2,2
               WRITE (IOTMP,'(40f20.10)') OMEGA(IOM),
     &               (DIMAG(SOM(LM,LM,IOM)),LM=1,NREP)
            END DO
            WRITE (IOTMP,*) '##'
            CALL DMFT_CLOSE_IOTMP(IOTMP)
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &                            '_sig_matsub_re.dat')
            CLOSE (IOTMP,STATUS='DELETE')
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &                            '_sig_matsub_re.dat')
            DO IOM = 2,NOMLOC/2,2
               WRITE (IOTMP,'(40f20.10)') OMEGA(IOM),
     &               (DREAL(SOM(LM,LM,IOM)),LM=1,NREP)
            END DO
            WRITE (IOTMP,*) '##'
            CALL DMFT_CLOSE_IOTMP(IOTMP)
         END IF
      ENDIF
C============================================================
C     Create self energy on complex contour (pade approx)
c     Creat self energy on real energy axis for check
C============================================================
      ALLOCATE (SEZNEW(NREP,NREP,NEGF))
      ALLOCATE (SEREAL(NREP,NREP,NEREAL))
      ALLOCATE (EREAL(NEREAL))
      ALLOCATE (CWORK(NEGF+NEREAL))
      ALLOCATE (CWORKO(NOMLOC/4))
      allocate (epath_tmp(NEGF+NEREAL))
      SEZNEW(1:NREP,1:NREP,1:NEGF) = C0
      SEREAL(1:NREP,1:NREP,1:NEREAL) = C0
      CWORK(1:NEGF+NEREAL) = C0
      CWORKO(1:NOMLOC/4) = C0

      epath_tmp(1:negf)=egfnew(1:negf)
      de=3d0/(nereal-1)
      ereal(1)=dcmplx(eferminew-1.5d0,0d0)
      do ie=2,nereal
        ereal(ie) = ereal(ie-1) + de
      enddo
      epath_tmp(negf+1:negf+nereal) = ereal(1:nereal)

      DO LM = 1,NREP
         DO LMP = 1,NREP
            IF ( IZERO(LMP,LM).EQ.1 ) THEN
               DO IOM = 1,NOMLOC/4
                  CWORKO(IOM) = SOM(LMP,LM,2*IOM)
               END DO
c
c modified by XJQ: more stable pade using code of Johan Schott and Elin Lundin
c
c               CALL DMFT_PADEOZ(OMEGA,CWORKO,EFERMINEW,NEGF,EGFNEW,
c     &                          CWORK,NOMLOC)
               call pade_averaging_interface(epath_tmp,1,negf+nereal,
     &           cwork,efermiold,temp,nomloc/4,cworko,
     &           2,it,lm,lmp,0,iotmp)
c
c end-mod-xjq
c
               SEZNEW(LMP,LM,1:NEGF) = CWORK(1:NEGF)
c we don't need self-energy of those energies having zero weight of charge density
c with zere self-energy, g_btah maybe more stable
               if(zerosig_contn .and. lactive_contn) 
     &           seznew(lmp,lm,ie0_contn:ie1_contn) = c0
               sereal(lmp,lm,1:nereal) = cwork(negf+1:negf+nereal)
C
            END IF
         END DO
      END DO
      DEALLOCATE (CWORK)
C         END IF
C============================================================
C      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
C         DO IE=1,NEGF
C      SEZNEW(:,:,IE)=SEZNEWS(:,:,IE) +SIGSTAT(1:NREP,1:NREP)
C         END DO
C      END IF
C
C============================================================
C     Mixing of the self energy
C============================================================
C
      IF ( ITRSCF.EQ.1 ) THEN
C
         DMFTSIGMA(I1:I2,I1:I2,1:NEGF) = SEZNEW(1:NLM,1:NLM,1:NEGF)
                                            !1_1
         DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
     &      = SEZNEW(NLM+1:NREP,NLM+1:NREP,1:NEGF)     !2_2
         DMFTSIGMA(I1:I2,I3:I4,1:NEGF) = SEZNEW(1:NLM,NLM+1:NREP,1:NEGF)
                                                  !1_2
         DMFTSIGMA(I3:I4,I1:I2,1:NEGF) = SEZNEW(NLM+1:NREP,1:NLM,1:NEGF)
                                                  !2_1
C
         SEZ(1:NREP,1:NREP,1:NEGF,IT) = SEZNEW(1:NREP,1:NREP,1:NEGF)
      ELSE
         DMFTSIGMA(I1:I2,I1:I2,1:NEGF) = SEZ(1:NLM,1:NLM,1:NEGF,IT)
     &      *(1D0-DMFTMIX) + SEZNEW(1:NLM,1:NLM,1:NEGF)*DMFTMIX
                                                    !1_1
         DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
     &      = SEZ(NLM+1:NREP,NLM+1:NREP,1:NEGF,IT)*(1D0-DMFTMIX)
     &      + SEZNEW(NLM+1:NREP,NLM+1:NREP,1:NEGF)*DMFTMIX    !2_2
         DMFTSIGMA(I1:I2,I3:I4,1:NEGF) = SEZ(1:NLM,NLM+1:NREP,1:NEGF,IT)
     &      *(1D0-DMFTMIX) + SEZNEW(1:NLM,NLM+1:NREP,1:NEGF)*DMFTMIX   !1_2
         DMFTSIGMA(I3:I4,I1:I2,1:NEGF) = SEZ(NLM+1:NREP,1:NLM,1:NEGF,IT)
     &      *(1D0-DMFTMIX) + SEZNEW(NLM+1:NREP,1:NLM,1:NEGF)*DMFTMIX
                                                                   !2_1
C
         SEZ(1:NREP,1:NREP,1:NEGF,IT) = SEZ(1:NREP,1:NREP,1:NEGF,IT)
     &                                  *(1D0-DMFTMIX)
     &                                  + SEZNEW(1:NREP,1:NREP,1:NEGF)
     &                                  *DMFTMIX
      END IF
C============================================================
C     print-out of the self energy
C============================================================
       if(parent_rank==0) then
         FILNAM = TXT_T(1:LTXT_T)//'_SIGUNFORM'//'.sig'
         LFILNAM = LTXT_T + 10 + 4
C
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),STATUS='UNKNOWN',
     &         FORM='unformatted')
C
         REWIND IOTMP
         WRITE (IOTMP) NOMLOC/2,DMFTTEMP
         WRITE (IOTMP) (((SOM(LM,LMP,2*IOM),LM=1,NREP),
     &                 LMP=1,NREP),IOM=1,NOMLOC/2)
         WRITE (IOTMP) EFERMIOLD
         WRITE (IOTMP) EREFDMFT
         CLOSE (IOTMP)
c
c         IF ( IPRINTDMFT.GT.0 ) THEN
            CALL DMFT_CLOSE_IOTMP(IOTMP)
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &           '_sigma_realz_im.dat')
            CLOSE (IOTMP,STATUS='DELETE')
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &           '_sigma_realz_im.dat')
            DO IE = 1,NEREAL
               WRITE (IOTMP,'(40f20.10)') EREAL(IE),
     &               (DIMAG(SEREAL(LMP,LMP,IE)),LMP=1,
     &                                    NREP)
            END DO
            CALL DMFT_CLOSE_IOTMP(IOTMP)
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &           '_sigma_realz_re.dat')
            CLOSE (IOTMP,STATUS='DELETE')
            OPEN (UNIT=IOTMP,FILE='diag-sig/'//TXT_T(1:LTXT_T)//
     &           '_sigma_realz_re.dat')
            DO IE = 1,NEREAL
               WRITE (IOTMP,'(40f20.10)') EREAL(IE),
     &               (DREAL(SEREAL(LMP,LMP,IE)),LMP=1,
     &                                    NREP)
            END DO
            CALL DMFT_CLOSE_IOTMP(IOTMP)
c         END IF
C
      endif ! parent_rank==0
C
      ICALL = ICALL + 1
      RETURN
C
99001 FORMAT (1x,20(I2))
99002 FORMAT (1x,A)
99003 FORMAT (/,1X,'===== Calculated on matsubara       ======',/,5X,
     &        '        EMIGDAL')
99004 FORMAT (/,1X,'===== Calculated in time domain     ======',/,5X,
     &        'TOTAL         FIRST       SECOND  CONTRIBUTION ')
99005 FORMAT (5X,5F14.6)
      END
