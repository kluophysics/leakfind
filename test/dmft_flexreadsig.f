C*==dmft_flexreadsig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_FLEXREADSIG(IOTMP,DMFTSIGMA,EGF,NKM,NEGF,IT,IREL,
     &                            DMFTTEMP,LOP,NKMMAX,NEMAX,NTMAX,
     &                            EREFDMFT)
      USE MOD_CONSTANTS,ONLY:C0,PI
      use mod_files,only:iotmp2
      USE MOD_TYPES,ONLY:TXT_T,LTXT_T
      use mod_dmft_ldau,only:zerosig_contn
      use mod_energy,only:igrid,nefd1,nefd2,nefd3,nepol,necontn,
     &                    lactive_contn
      use mod_mpi_multilevels
      IMPLICIT NONE
C*--DMFT_FLEXREADSIG7
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EVOVRY,TOL
      PARAMETER (EVOVRY=13.605826D0,TOL=1.0D-10)
C
C Dummy arguments
C
      REAL*8 DMFTTEMP
      COMPLEX*16 EREFDMFT
      INTEGER IOTMP,IREL,IT,LOP,NEGF,NEMAX,NKM,NKMMAX,NTMAX
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),EGF(NEMAX)
C
C Local variables
C
      COMPLEX*16 CTMP,CWORK(:),CWORKO(:),SEZ(:,:,:,:),SOM(:,:,:,:)
      LOGICAL DMFTMINB
      REAL*8 DOMEGA,DTIME,EFERMIOLD,OMEGA(:),RTMP,TEMP,TIME(:)
      INTEGER I1,I2,I3,I4,IE,IKM,IKM1,IKMP,IOM,IS,IZERO(:,:,:),LM,LMP,
     &        MINOM(:),NLM,NLMOP,NOM,NOMHALF,NSPIN
      character*80 filnam
      integer lfilnam, ns_tmp, ne_tmp
      integer level, parent_comm, nprocs_parent, parent_rank
      logical lexist, lparalleled
      complex*16, dimension(:), allocatable :: egf_tmp
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      ALLOCATABLE OMEGA,TIME,MINOM,SOM,SEZ,IZERO
      ALLOCATABLE CWORK,CWORKO
      DMFTMINB = .FALSE.
      REWIND (IOTMP)
      READ (IOTMP) NOMHALF,RTMP
      IF ( ABS(DMFTTEMP-RTMP).GT.TOL ) THEN
         WRITE (6,99001)
         WRITE (6,99002) DMFTTEMP,RTMP
         WRITE (6,99003)
         DMFTTEMP = RTMP
      END IF
      IF ( DMFTMINB ) WRITE (6,*) 'TEST OPTION: SIGMA INTERCHANGE SPINS'
C
      IF ( IREL.LT.2 ) STOP '<DMFT> only for irel >= 2'
C
      IF ( IREL.EQ.2 ) THEN
         NSPIN = 2
         NLM = LOP*2 + 1
         I1 = LOP**2 + 1
         I2 = I1 + NLM - 1
         I3 = NKM/2 + I1
         I4 = I3 + NLM - 1
      ELSE
         NSPIN = 1
         NLM = 2*(LOP*2+1)
         NLMOP = LOP*2 + 1
         I1 = 2*LOP**2 + 1
      END IF
C============================================================
C     Create matsubara mesh
C============================================================
      NOM = NOMHALF*2
c      TEMP = DMFTTEMP/11605.0D0/EVOVRY
      TEMP = DMFTTEMP*0.6333659D-5
      DTIME = 1.D0/NOM/TEMP*2.D0
      DOMEGA = TEMP*PI
C
      ALLOCATE (TIME(NOM),OMEGA(NOM),MINOM(NOM))
C
      CALL DMFT_MATSUB(DTIME,DOMEGA,TIME,OMEGA,MINOM,NOM)
C============================================================
C     Read in self-energy from scf-cycle
C============================================================
      ALLOCATE (SOM(NLM,NLM,NOM,NSPIN))
C
      READ (IOTMP) ((((SOM(LM,LMP,2*IOM,IS),LM=1,NLM),LMP=1,NLM),IOM=1,
     &             NOM/2),IS=1,NSPIN)
      READ (IOTMP) EFERMIOLD
      READ (IOTMP,ERR=100) EREFDMFT
      GOTO 200
 100  CONTINUE
      EREFDMFT = DCMPLX(0.7D0,0.0D0)
C
C============================================================
C     Check which elements of SOM are non_zero
C============================================================
 200  CONTINUE
      ALLOCATE (IZERO(NLM,NLM,NSPIN))
      IZERO = 0
      DO IS = 1,NSPIN
         DO LM = 1,NLM
            DO LMP = 1,NLM
               CTMP = C0
               DO IOM = 1,NOM/2
                  CTMP = CTMP + SOM(LM,LMP,2*IOM,IS)
               END DO
               IF ( ABS(CTMP).GT.TOL ) IZERO(LM,LMP,IS) = 1
            END DO
         END DO
      END DO
C============================================================
C     Create self energy on complex contour (pade approx)
C============================================================
      ALLOCATE (SEZ(NLM,NLM,NEGF,NSPIN))
      FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'_SIG_EPATH_UNFORM'//'.sig'
      LFILNAM = LTXT_T(IT) + 17 + 4
      inquire(file=filnam(1:lfilnam),exist=lexist)
      if(lexist) then
        OPEN (IOTMP2,FILE=FILNAM(1:LFILNAM),STATUS='OLD',
     &        FORM='unformatted')
        read(iotmp2) ns_tmp, ne_tmp
        if(ns_tmp==nspin .and. ne_tmp==negf) then
          allocate(egf_tmp(1:negf))
          do is=1,nspin
            do ie=1,negf
              read(iotmp2) egf_tmp(ie),((sez(lm,lmp,ie,is),lm=1,nlm),
     &                                 lmp=1,nlm)
              if( abs( egf_tmp(ie) - egf(ie) ) > 1e-8 ) then
                lexist=.false.
                exit
              endif
            enddo
          enddo
          deallocate(egf_tmp)
        else
          lexist=.false.
        endif
        close(iotmp2)
      endif
c
      if(.not.lexist) then 
c
        ALLOCATE (CWORK(NEGF))
        ALLOCATE (CWORKO(NOM/4))
        SEZ = C0
        CWORK = C0
        CWORKO = C0
        DO IS = 1,NSPIN
           DO LM = 1,NLM
              DO LMP = 1,NLM
                 IF ( IZERO(LMP,LM,IS).EQ.1 ) THEN
                    DO IOM = 1,NOM/4
                       CWORKO(IOM) = SOM(LMP,LM,2*IOM,IS)
                    END DO
c                    CALL DMFT_PADEOZ(OMEGA,CWORKO,EFERMIOLD,NEGF,EGF,
c     &                               CWORK,NOM)
                    call pade_averaging_interface(egf,1,negf,cwork,
     &                efermiold,temp,nom/4,cworko,2,it,lm,lmp,0,iotmp2)
                    do ie=1,negf
                      SEZ(LMP,LM,ie,IS) = CWORK(ie)
                    enddo
c we don't need self-energy of those energies having zero weight of charge density
c with zere self-energy, g_btah maybe more stable
                    if(zerosig_contn .and. lactive_contn) sez(lmp,lm,
     &                nefd1+nefd2+nefd3+1:negf-nepol,is) = c0
                 END IF
              END DO
           END DO
        END DO      
        DEALLOCATE (CWORK,CWORKO)
c
c write sigma on real-energy axis to unformatted files
        level = toplevel + 1
        call get_comm_level(level,'parent',lparalleled,
     &    parent_comm,nprocs_parent,parent_rank)
        if(parent_rank==0) then
          FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'_SIG_EPATH_UNFORM'//'.sig'
          LFILNAM = LTXT_T(IT) + 17 + 4
          OPEN (IOTMP2,FILE=FILNAM(1:LFILNAM),STATUS='replace',
     &          FORM='unformatted')
          write(iotmp2) nspin,negf
          do is=1,nspin
            do ie=1,negf
              write(iotmp2) egf(ie),((sez(lm,lmp,ie,is),lm=1,nlm),
     &                              lmp=1,nlm)
            enddo
          enddo
          close(iotmp2)
        endif
      endif
C============================================================
C
C
      IF ( IREL.LE.2 ) THEN
         DO IS = 1,NSPIN
            IF ( IS.EQ.1 ) THEN
               IKM1 = I1
            ELSE
               IKM1 = I3
            END IF
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  IKM = IKM1 + LM - 1
                  IKMP = IKM1 + LMP - 1
                  DO IE = 1,NEGF
                     DMFTSIGMA(IKM,IKMP,IT,IE) = SEZ(LM,LMP,IE,IS)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         I1 = LOP**2 + 1
         I2 = I1 + NLMOP - 1
         I3 = NKM/2 + I1
         I4 = I3 + NLMOP - 1
         DO IE = 1,NEGF
            IF ( DMFTMINB ) THEN
               DMFTSIGMA(I1:I2,I1:I2,IT,IE) = SEZ(1:NLMOP,1:NLMOP,IE,1)
               DMFTSIGMA(I3:I4,I3:I4,IT,IE) = SEZ(1:NLMOP,1:NLMOP,IE,1)
               DMFTSIGMA(I1:I2,I3:I4,IT,IE) = SEZ(1:NLMOP,1:NLMOP,IE,1)
               DMFTSIGMA(I3:I4,I1:I2,IT,IE) = SEZ(1:NLMOP,1:NLMOP,IE,1)
            ELSE
               DMFTSIGMA(I1:I2,I1:I2,IT,IE) = SEZ(1:NLMOP,1:NLMOP,IE,1)
               DMFTSIGMA(I3:I4,I3:I4,IT,IE)
     &            = SEZ(NLMOP+1:NLM,NLMOP+1:NLM,IE,1)
               DMFTSIGMA(I1:I2,I3:I4,IT,IE)
     &            = SEZ(1:NLMOP,NLMOP+1:NLM,IE,1)
               DMFTSIGMA(I3:I4,I1:I2,IT,IE)
     &            = SEZ(NLMOP+1:NLM,1:NLMOP,IE,1)
            END IF
         END DO
      END IF
C============================================================
C
99001 FORMAT (/,1X,79('*'),/,25X,'WARNING: DMFT_FLEXREADSIG',A,/)
99002 FORMAT ('from input: DMFTTEMP=',F10.4,' from file:',F10.4,/,
     &        'Temprature from file will be used')
99003 FORMAT (/,1X,79('*'),/)
C
      END
