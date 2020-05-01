C*==dmft_drv_sosptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_DRV_SOSPTFLEX(GFMAT,DMFTSIGMA,NKM,NEGF,LOP,UEFF,
     &                              JEFF,UDMFT,DMFTTEMP,DMFTDBLC,NOM,
     &                              IOTMP,IPRINTDMFT,IFLEX,TXT_T,LTXT_T,
     &                              ITRSCF,DMFTMIX,NKMMAX,NLMAX,WEGFOLD,
     &                              EGFOLD,EGFNEW,EREFDMFT,IT,NTMAX,
     &                              TMA_STATIC,TMA_NOBATH,TMA_DBLC,
     &                              DOSINT,TMA_NOCCSCL)
C
C
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      use mod_energy,only:igrid,nefd1,nefd2,nefd3,nepol,necontn,
     &                    lactive_contn
      use mod_mpi,only:mpi_id
      IMPLICIT NONE
C*--DMFT_DRV_SOSPTFLEX13
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1.0D-10)
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
     &           ,GFMAT(NKM,NKM,NEGF),WEGFOLD(NEGF)
      REAL*8 DOSINT(NLMAX),UDMFT(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,CWORK(:),CWORK1(:,:),CWORK2(:,:),CWORKO(:),EMIGD,
     &           EREAL(:),GOM(:,:,:),GOMIN(:,:,:),GOMWORK(:,:,:),
     &           GZ(:,:,:),SEREAL(:,:,:),SEZ(:,:,:,:),SEZNEW(:,:,:),
     &           SOM(:,:,:),SOMWORK(:,:,:)
      INTEGER DC,I1,I2,I3,I4,ICALL,IE,IOM,IR,IR1,IR2,IRP,IS1,IS2,
     &        IZERO(:,:),IZEROSIG(:,:),LFILNAM,LM,LMP,LMPP,LMPPP,
     &        MINOM(:),MIOM,NEREAL,NLM,NOMIN,NOMLOC,NREP
      REAL*8 DE,DOMEGA,DTIME,EFERMINEW,EFERMIOLD,EMIGDR(3),OCC(:,:),
     &       OMEGA(:),RWORK,SIGSTAT(:,:),TEMP,TIME(:),TMP,TMP1,TMP2,UJJ,
     &       UR(:,:,:,:),URREP(:,:,:,:),UUU
      CHARACTER*80 FILNAM
      SAVE IZERO,IZEROSIG,SEZ
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE GZ,CWORK1,CWORK2,SEZ,SOM,SEZNEW
      ALLOCATABLE OMEGA,TIME,MINOM,IZERO,GOM,CWORK,CWORKO,GOMIN
      ALLOCATABLE SEREAL,EREAL,UR,GOMWORK,SOMWORK,URREP
      ALLOCATABLE IZEROSIG
      ALLOCATABLE OCC,SIGSTAT
C
      NEREAL = 500
C
C
      NLM = LOP*2 + 1
      NREP = 2*NLM
C
      TEMP = DMFTTEMP*0.6333659D-5
      UUU = UEFF/RY_EV
      UJJ = JEFF/RY_EV
C
      EFERMIOLD = DREAL(EGFOLD(NEGF))
      EFERMINEW = DREAL(EGFNEW(NEGF))
C
      IF ( IFLEX.EQ.0 ) WRITE (6,99002) 
     &                'Only TMA-approximation used (PH-channel ignored)'
C
C
C
      IF ( ICALL.EQ.0 ) ALLOCATE (SEZ(NREP,NREP,NEGF,NTMAX))
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
         ALLOCATE (IZEROSIG(NREP,NREP))
         IZEROSIG(1:NREP,1:NREP) = 0
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
C============================================================
C     G_BATH=G^-1+SIGMA
C============================================================
C
C
      ALLOCATE (CWORK1(NREP,NREP),CWORK2(NREP,NREP))
C
      IF ( .NOT.TMA_NOBATH ) THEN
         DO IE = 1,NEGF
            CWORK1(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,IE)
            CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
            CWORK1(1:NREP,1:NREP) = CWORK2(1:NREP,1:NREP)
     &                              + SEZ(1:NREP,1:NREP,IE,IT)
            CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
            GZ(1:NREP,1:NREP,IE) = CWORK2(1:NREP,1:NREP)
         END DO
      END IF
      DEALLOCATE (CWORK1,CWORK2)
C
C
C
      IF ( IPRINTDMFT.GT.0 ) THEN
         DO LM = 1,NREP
            WRITE (6,99001) (IZERO(LM,LMP),LMP=1,NREP)
         END DO
      END IF
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
         WRITE (*,'(/,10X,A,F8.4)') 'Trace:',TMP
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
               CWORK(1:NEGF) = GZ(LMP,LM,1:NEGF)
c
c modified by XJQ: more stable pade using code of Johan Schott and Elin Lundin,
c                  their subroutine pade is paralleled so that speed is not slow,
c                  subroutine pade_averaging_interface is an interface for sprkkr
c
               if(.not. lactive_contn) then
c                 CALL DMFT_PADEZO(EFERMIOLD,NEGF,EGFOLD,CWORK,
c     &                            OMEGA,MINOM,CWORKO,NOMLOC)
                 call pade_averaging_interface(egfold,1,negf,cwork,
     &             efermiold,temp,nomloc,cworko,1,it,lm,lmp,0,iotmp)
               else
                 call pade_averaging_interface(egfold,
     &             nefd1+nefd2+nefd3+1,negf,cwork,
     &             efermiold,temp,nomloc,cworko,1,it,lm,lmp,0,iotmp)
               endif
C
c the following copy part may be not needed, it can be closed by "if(.false.) then"
               if(.true.) then
               if(igrid(1)==11) then
                 do iom=1,nepol
                   CWORKO(2*IOM) = GZ(LMP,LM,negf-iom+1)
                   GOM(LMP,LM,2*IOM) = CWORKO(2*IOM)
                 enddo
                 do iom=nepol+1,NOMLOC/2
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
                     ENDIF
                  END DO
                  GOM(LMP,LM,2*IOM) = CWORKO(2*IOM)
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
C
C
C=================================================================
C     Calculate self-energy  on matsubara frequencies
C=================================================================
C
      ALLOCATE (SOM(NREP,NREP,NOMLOC))
      SOM(1:NREP,1:NREP,1:NOMLOC) = C0
C
      DC = 1
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) DC = 1
      IF ( DMFTDBLC(1:4).EQ.'META' ) DC = 1
      IF ( DMFTDBLC(1:4).EQ.'ASPN' ) DC = 4
      IF ( DMFTDBLC(1:4).EQ.'METC' ) DC = 5
C============================================================
C     Setup 4-index U matrix
C============================================================
      ALLOCATE (UR(NLM,NLM,NLM,NLM))
      UR(1:NLM,1:NLM,1:NLM,1:NLM) = 0.0D0
      ALLOCATE (URREP(NREP,NREP,NREP,NREP))
      URREP(1:NREP,1:NREP,1:NREP,1:NREP) = 0.0D0
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
            URREP((IS1-1)*NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM,(IS1-1)
     &         *NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM)
     &         = UR(1:NLM,1:NLM,1:NLM,1:NLM)
         END DO
      END DO
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
      SIGSTAT(1:NREP,1:NREP) = 0D0
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
C
C
         TMP1 = 0D0
         TMP2 = 0D0
         DO IR = 1,NREP/2
            TMP1 = TMP1 + OCC(IR,IR)
            TMP2 = TMP2 + OCC(IR+NREP/2,IR+NREP/2)
         END DO
         TMP1 = 2D0*TMP1/NREP
         TMP2 = 2D0*TMP2/NREP
C
         DO IR = 1,NREP
            DO IRP = 1,NREP
C full static electron interaction:
C     SIGMA_{12}=\sum_{34}(V_{1324}-V_{1342})*n_{34}
               DO IR1 = 1,NREP
                  DO IR2 = 1,NREP
                     SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP)
     &                                 + (URREP(IR,IR1,IRP,IR2)
     &                                 -URREP(IR,IR1,IR2,IRP))
     &                                 *OCC(IR1,IR2)
                  END DO
               END DO
C subtract double-counting /Tr<SIGMA>/:
               IF ( IR.EQ.IRP ) THEN
                  IF ( TMA_DBLC(1:3).EQ.'AMF' ) THEN
C        AMF case:
C        SIGMA_{11}=SIGMA_{11}-\sum_{2}(V_{1212}-V_{1221})*<n_{2}>
                     DO IR1 = 1,NREP/2
                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                     - (URREP(IR,IR1,IR,IR1)-URREP(IR,IR1,IR1,IR))
     &                     *TMP1
                     END DO
                     DO IR1 = NREP/2 + 1,NREP
                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                     - (URREP(IR,IR1,IR,IR1)-URREP(IR,IR1,IR1,IR))
     &                     *TMP2
                     END DO
                  ELSE IF ( TMA_DBLC(1:3).EQ.'AAL' ) THEN
C        AAL case:
C        SIGMA_{11}=SIGMA_{11}-U(N-1/2)+J(N_{1}-1/2)
                     SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                                - UUU*((TMP1+TMP2)*NREP/2-0.5D0)
                     IF ( IR.LE.NREP/2 ) THEN
                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                     + UJJ*(TMP1*NREP/2-0.5D0)
                     ELSE
                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                     + UJJ*(TMP2*NREP/2-0.5D0)
                     END IF
                  ELSE
                     STOP 'DMFT_DRV_SOSPTFLEX: CHECK TMA_DBLC !'
                  END IF
               END IF
                  ! IR=IRP
            END DO
         END DO
C
C         DO IR = 1,NREP/2
C            DO IR1 = 1,NREP/2
C               SIGSTAT(IR,IR1) = SIGSTAT(IR,IR1)
C     &                          + (OCC(IR1+NREP/2,IR1+NREP/2)-TMP2)
C     &                          *UR(IR,IR1,IR,IR1)
C               IF ( IR1.NE.IR ) SIGSTAT(IR,IR1) = SIGSTAT(IR,IR1)
C     &              + (UR(IR,IR1,IR,IR1)-UR(IR,IR,IR1,IR1))
C     &              *(OCC(IR1,IR1)-TMP1)
C               SIGSTAT(IR+NREP/2,IR1+NREP/2)
C     &            = SIGSTAT(IR+NREP/2,IR1+NREP/2) + (OCC(IR1,IR1)-TMP1)
C     &            *UR(IR,IR1,IR,IR1)
C               IF ( IR1.NE.IR ) SIGSTAT(IR+NREP/2,IR1+NREP/2)
C     &              = SIGSTAT(IR+NREP/2,IR1+NREP/2)
C     &              + (UR(IR,IR1,IR,IR1)-UR(IR,IR,IR1,IR1))
C     &              *(OCC(IR1+NREP/2,IR1+NREP/2)-TMP2)
C            END DO
C         END DO
         IF ( IPRINTDMFT.GT.0 ) CALL RMATSTRUCT('Static self-energy:',
     &        SIGSTAT,NREP,NREP,0,0,0,TOL,6)
C
      END IF
            ! DMFTDBLC='LDAU'
C
C
C
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
         ALLOCATE (GOMWORK(NREP,NREP,0:NOMLOC/2+1))
         ALLOCATE (SOMWORK(NREP,NREP,0:NOMLOC/2+1))
         GOMWORK(1:NREP,1:NREP,0:NOMLOC/2+1) = C0
         SOMWORK(1:NREP,1:NREP,0:NOMLOC/2+1) = C0
         NOMIN = 0
         DO IOM = 2,NOMLOC/2,2
            GOMWORK(1:NREP,1:NREP,NOMIN) = GOM(1:NREP,1:NREP,IOM)
            NOMIN = NOMIN + 1
         END DO
C         OPEN (UNIT=47,FILE='gf.dat')
C         WRITE (47,*) 'U=(Ry)',UUU,'J=(Ry)',UJJ
C         WRITE (47,*) 'TEMP=',TEMP
C         WRITE (47,*) 'NOM=',NOMIN,'NREP=',NREP
C         DO IOM = 0,NOMIN
C            WRITE (47,'(I4,200f16.10)') IOM,
C     &                                  ((GOMWORK(LM,LMP,IOM),LM=1,NREP)
C     &                                  ,LMP=1,NREP)
C         END DO
C         CLOSE (47)
C
         CALL DMFT_SOSPTFLEX(GOMWORK,SOMWORK,UR,NREP,NREP/2,NOMIN,
     &                       EMIGDR,TEMP,DC,IPRINTDMFT,SIGSTAT)
C
         NOMIN = 0
C
C         OPEN (UNIT=14,FILE='som.dat')
C         DO IOM = 0,NOMIN
C            WRITE (14,'(I4,200f16.10)') IOM,
C     &                                  ((SOMWORK(LM,LMP,IOM),LM=1,NREP)
C     &                                  ,LMP=1,NREP)
C         END DO
C         CALL FLUSH(14)
C         CLOSE (14)
         DO IOM = 2,NOMLOC/2,2
            DO LM = 1,NREP
               DO LMP = 1,NREP
                  SOM(LMP,LM,IOM) = SOMWORK(LMP,LM,NOMIN)
               END DO
            END DO
            NOMIN = NOMIN + 1
         END DO
C
C
         DO IOM = 1,NOMLOC/2
            MIOM = NOMLOC/2 + IOM
            SOM(1:NREP,1:NREP,MIOM) = DCONJG(SOM(1:NREP,1:NREP,IOM))
         END DO
C
         DO LM = 1,NREP
            DO LMP = 1,NREP
               CSUM = C0
               DO IE = 1,NOMLOC/4
                  CSUM = CSUM + SOM(LMP,LM,2*IE)
               END DO
               IF ( ABS(CSUM).GE.TOL ) IZEROSIG(LMP,LM) = 1
            END DO
         END DO
         IF ( IPRINTDMFT.GT.0 ) THEN
            DO LM = 1,NREP
               WRITE (6,99001) (IZEROSIG(LM,LMP),LMP=1,NREP)
            END DO
         END IF
         DO LM = 1,NREP
            DO LMP = 1,NREP
               IF ( IZERO(LM,LMP).EQ.0 ) SOM(LM,LMP,:) = C0
            END DO
         END DO
C  TMA_STATIC
      ELSE
         WRITE (6,*) 'TMA_STATIC=TRUE, e.g. Only Static part used'
         SOM(:,:,:) = C0
      END IF
C
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
      IF ( IPRINTDMFT.GT.0 ) THEN
         WRITE (6,99003)
         WRITE (6,99005) DREAL(EMIGD),DIMAG(EMIGD)
         WRITE (6,99004)
         WRITE (6,99005) EMIGDR(3),EMIGDR(1),EMIGDR(2)
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     print-out of the self energy
C
      IF ( MPI_ID==0 .and. IPRINTDMFT.GT.0 ) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_matsub_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_matsub_im.dat')
         DO IOM = 2,NOMLOC,2
            WRITE (IOTMP,'(11f20.10)') OMEGA(IOM),
     &                                 (DIMAG(SOM(LM,LM,IOM)),LM=1,NREP)
         END DO
         WRITE (IOTMP,*) '##'
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_matsub_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_matsub_re.dat')
         DO IOM = 2,NOMLOC,2
            WRITE (IOTMP,'(11f20.10)') OMEGA(IOM),
     &                                 (DREAL(SOM(LM,LM,IOM)),LM=1,NREP)
         END DO
         WRITE (IOTMP,*) '##'
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C============================================================
C     Create self energy on complex contour (pade approx)
C============================================================
      ALLOCATE (CWORK(NEGF))
      ALLOCATE (CWORKO(NOMLOC/4))
      ALLOCATE (SEZNEW(NREP,NREP,NEGF))
      SEZNEW(1:NREP,1:NREP,1:NEGF) = C0
      CWORK(1:NEGF) = C0
      CWORKO(1:NOMLOC/4) = C0
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
               call pade_averaging_interface(egfnew,1,negf,cwork,
     &           eferminew,temp,nomloc/4,cworko,2,it,lm,lmp,0,iotmp)
c
c end-mod-xjq
c
               SEZNEW(LMP,LM,1:NEGF) = CWORK(1:NEGF)
C
C    IF SOME KKR ENERGY POINTS DONE ON  MATSUBARA JUST TAKE THEM OVER
               DO IOM = 1,NOMLOC/2
                  DO IE = 1,NEGF
                     IF ( ABS(EFERMINEW-DREAL(EGFNEW(IE))).LT.1D-8 .AND. 
     &                    ABS(OMEGA(2*IOM)-DIMAG(EGFNEW(IE))).LT.1D-8 )
     &                    SEZNEW(LMP,LM,IE) = SOM(LMP,LM,2*IOM)
                  END DO
               END DO
            END IF
         END DO
      END DO
      DEALLOCATE (CWORK)
C============================================================
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
C============================================================
C     Create self energy on real axis (just for check)
C============================================================
      IF ( MPI_ID==0 .and. IPRINTDMFT.GT.1 ) THEN
         ALLOCATE (SEREAL(NREP,NREP,NEREAL))
         ALLOCATE (EREAL(NEREAL))
         ALLOCATE (CWORK(NEREAL))
         SEREAL(1:NREP,1:NREP,1:NEREAL) = C0
         CWORK(1:NEREAL) = C0
         CWORKO(1:NEREAL) = C0
         DE = (2.0D0+1.0D0)/DBLE(NEREAL-1)
         DO IE = 1,NEREAL
            RWORK = -1.0D0 + DE*DBLE(IE-1)
            EREAL(IE) = DCMPLX(RWORK,0.001D0)
         END DO
         DO LM = 1,NREP
            DO LMP = 1,NREP
               IF ( IZERO(LMP,LM).EQ.1 ) THEN
                  DO IOM = 1,NOMLOC/2
                     CWORKO(IOM) = SOM(LMP,LM,IOM*2)
                  END DO
c                   CALL DMFT_PADEOZ(OMEGA,CWORKO,EFERMINEW,NEREAL,
c     &                              EREAL,CWORK,NOMLOC)
                   call pade_averaging_interface(ereal,1,nereal,cwork,
     &               eferminew,temp,nomloc/2,cworko,2,it,lm,lmp,0,iotmp)
                  SEREAL(LMP,LM,1:NEREAL) = CWORK(1:NEREAL)
               END IF
            END DO
         END DO
      END IF
C============================================================
      ICALL = ICALL + 1
C
C============================================================
C     print-out of the self energy
C============================================================
      if( MPI_ID==0 ) then
      FILNAM = TXT_T(1:LTXT_T)//'_SIGUNFORM'//'.sig'
      LFILNAM = LTXT_T + 10 + 4
C
      OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),STATUS='UNKNOWN',
     &      FORM='unformatted')
C
      REWIND IOTMP
      WRITE (IOTMP) NOMLOC/2,DMFTTEMP
      WRITE (IOTMP) (((SOM(LM,LMP,2*IOM),LM=1,NREP),LMP=1,NREP),IOM=1,
     &              NOMLOC/2)
      WRITE (IOTMP) EFERMIOLD
      WRITE (IOTMP) EREFDMFT
      CLOSE (IOTMP)
      endif
C
C
      IF ( MPI_ID==0 .and. IPRINTDMFT.GT.0 .and. .false.) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_egfnew_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_egfnew_im.dat')
         DO IE = 1,NEGF
            WRITE (IOTMP,'(30f20.10)') EGFNEW(IE),
     &                                 (DIMAG(DMFTSIGMA(LMP,LMP,IE)),
     &                                 LMP=1,NKM)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_egfnew_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sig_egfnew_re.dat')
         DO IE = 1,NEGF
            WRITE (IOTMP,'(30f20.10)') EGFNEW(IE),
     &                                 (DREAL(DMFTSIGMA(LMP,LMP,IE)),
     &                                 LMP=1,NKM)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C
      IF ( MPI_ID==0 .and. IPRINTDMFT.GT.1 ) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sigma_realz_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sigma_realz_im.dat')
         DO IE = 1,NEREAL
            WRITE (IOTMP,'(40f20.10)') EREAL(IE),
     &                                 (DIMAG(SEREAL(LMP,LMP,IE)),LMP=1,
     &                                 NREP)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sigma_realz_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE=TXT_T(1:LTXT_T)//'_sigma_realz_re.dat')
         DO IE = 1,NEREAL
            WRITE (IOTMP,'(40f20.10)') EREAL(IE),
     &                                 (DREAL(SEREAL(LMP,LMP,IE)),LMP=1,
     &                                 NREP)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C
99001 FORMAT (1x,20(I2))
99002 FORMAT (1x,A)
99003 FORMAT (/,1X,'===== Calculated on matsubara       ======',/,5X,
     &        '        EMIGDAL')
99004 FORMAT (/,1X,'===== Calculated in time domain     ======',/,5X,
     &        'TOTAL         FIRST       SECOND  CONTRIBUTION ')
99005 FORMAT (5X,5F14.6)
      END
