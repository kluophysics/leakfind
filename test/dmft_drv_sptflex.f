C*==dmft_drv_sptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_DRV_SPTFLEX(GFMAT,DMFTSIGMA,NKM,NEGF,EGFOLD,
     &                            EGFNEW,LOP,UEFF,JEFF,DMFTTEMP,
     &                            DMFTDBLC,NOM,NSPIN,IOTMP,IPRINT,IFLEX,
     &                            TXT_T,LTXT_T,ITRSCF,DMFTMIX,NKMMAX)
C
C
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      IMPLICIT NONE
C*--DMFT_DRV_SPTFLEX10
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
      INTEGER IFLEX,IOTMP,IPRINT,ITRSCF,LOP,LTXT_T,NEGF,NKM,NKMMAX,NOM,
     &        NSPIN
      CHARACTER*8 TXT_T
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NEGF),EGFNEW(NEGF),EGFOLD(NEGF)
     &           ,GFMAT(NKM,NKM,NEGF)
C
C Local variables
C
      COMPLEX*16 CSUM,CWORK(:),CWORK1(:,:),CWORK2(:,:),CWORKO(:),EMIGD,
     &           EREAL(:),GOM(:,:,:,:),GOMIN(:,:,:,:),GZ(:,:,:,:),
     &           SEREAL(:,:,:,:),SEZ(:,:,:,:),SEZNEW(:,:,:,:),
     &           SOM(:,:,:,:)
      REAL*8 DE,DOMEGA,DTIME,EFERMINEW,EFERMIOLD,FAC1,FAC2,OMEGA(:),
     &       RWORK,TEMP,TIME(:),UJJ,UUU
      CHARACTER*80 FILNAM
      INTEGER I1,I2,I3,I4,IE,IKM,IKM1,IKMP,IOM,IPRINTDMFT,IS,
     &        IZERO(:,:,:),LFILNAM,LM,LMP,MINOM(:),NEREAL,NLM,NLMQU,
     &        NLMSQ
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE GZ,CWORK1,CWORK2,SEZ,SOM,SEZNEW
      ALLOCATABLE OMEGA,TIME,MINOM,IZERO,GOM,CWORK,CWORKO,GOMIN
      ALLOCATABLE SEREAL,EREAL
C
      NEREAL = 500
C
      IPRINTDMFT = IPRINT
C
      NLM = LOP*2 + 1
      NLMSQ = NLM**2
      NLMQU = NLMSQ**2
C
      TEMP = DMFTTEMP/11605.0D0/RY_EV
      UUU = UEFF/RY_EV
      UJJ = JEFF/RY_EV
C
      EFERMIOLD = DREAL(EGFOLD(NEGF))
      EFERMINEW = DREAL(EGFNEW(NEGF))
C
      IF ( IFLEX.EQ.0 ) WRITE (6,99003) 
     &                      'Only TMA-approximation used (FLEX ignored)'
C
C
      ALLOCATE (GZ(NLM,NLM,NEGF,NSPIN))
      ALLOCATE (SEZ(NLM,NLM,NEGF,NSPIN))
C
      SEZ(1:NLM,1:NLM,1:NEGF,1:NSPIN) = C0
      I1 = LOP**2 + 1
      I2 = I1 + NLM - 1
      I3 = NKM/2 + I1
      I4 = I3 + NLM - 1
      GZ(1:NLM,1:NLM,1:NEGF,1) = GFMAT(I1:I2,I1:I2,1:NEGF)
      IF ( NSPIN.EQ.2 ) GZ(1:NLM,1:NLM,1:NEGF,2)
     &     = GFMAT(I3:I4,I3:I4,1:NEGF)
      SEZ(1:NLM,1:NLM,1:NEGF,1) = DMFTSIGMA(I1:I2,I1:I2,1:NEGF)
      IF ( NSPIN.EQ.2 ) SEZ(1:NLM,1:NLM,1:NEGF,2)
     &     = DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
C============================================================
C     G_BATH=G^-1+SIGMA
C============================================================
      ALLOCATE (CWORK1(NLM,NLM),CWORK2(NLM,NLM))
C
      DO IS = 1,NSPIN
         DO IE = 1,NEGF
C
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  CWORK1(LMP,LM) = GZ(LMP,LM,IE,IS)
               END DO
            END DO
            CALL CMATINV(NLM,NLM,CWORK1,CWORK2)
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  CWORK1(LMP,LM) = CWORK2(LMP,LM) + SEZ(LMP,LM,IE,IS)
               END DO
            END DO
            CALL CMATINV(NLM,NLM,CWORK1,CWORK2)
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  GZ(LMP,LM,IE,IS) = CWORK2(LMP,LM)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (CWORK1,CWORK2)
C============================================================
C============================================================
C     Check which elements of GZ are non_zero
C============================================================
      ALLOCATE (IZERO(NLM,NLM,NSPIN))
      IZERO(1:NLM,1:NLM,2) = 0
      DO IS = 1,NSPIN
         DO LM = 1,NLM
            DO LMP = 1,NLM
               CSUM = C0
               DO IE = 1,NEGF
                  CSUM = CSUM + GZ(LMP,LM,IE,IS)
               END DO
               IF ( ABS(CSUM).GE.TOL ) IZERO(LMP,LM,IS) = 1
            END DO
         END DO
      END DO
      IF ( IPRINTDMFT.GT.0 ) THEN
         DO IS = 1,NSPIN
            WRITE (6,99001) IS
            DO LM = 1,NLM
               WRITE (6,99002) (IZERO(LM,LMP,IS),LMP=1,NLM)
            END DO
         END DO
      END IF
C
C
C============================================================
C     Create matsubara mesh
C============================================================
      ALLOCATE (OMEGA(NOM),TIME(NOM),MINOM(NOM))
      DOMEGA = TEMP*PI
      DTIME = 1.D0/NOM/TEMP*2.D0
      CALL DMFT_MATSUB(DTIME,DOMEGA,TIME,OMEGA,MINOM,NOM)
C
C============================================================
C     Create Green's function on matsubara (pade approx)
C============================================================
      ALLOCATE (GOM(NLM,NLM,NOM,NSPIN))
      ALLOCATE (GOMIN(NLM,NLM,NOM,NSPIN))
      ALLOCATE (CWORK(NEGF))
      ALLOCATE (CWORKO(NOM))
      GOM(1:NLM,1:NLM,1:NOM,1:2) = C0
      GOMIN(1:NLM,1:NLM,1:NOM,1:2) = C0
      CWORK(1:NEGF) = C0
      CWORKO(1:NOM) = C0
C
      DO IS = 1,NSPIN
         DO LM = 1,NLM
            DO LMP = 1,NLM
               IF ( IZERO(LMP,LM,IS).EQ.1 ) THEN
                  CWORK(1:NEGF) = GZ(LMP,LM,1:NEGF,IS)
                  CALL DMFT_PADEZO(EFERMIOLD,NEGF,EGFOLD,CWORK,OMEGA,
     &                             MINOM,CWORKO,NOM)
                  DO IOM = 1,NOM/2
                     GOM(LMP,LM,2*IOM,IS) = CWORKO(2*IOM)
                  END DO
               END IF
            END DO
         END DO
      END DO
      GOMIN(1:NLM,1:NLM,1:NOM,1:2) = GOM(1:NLM,1:NLM,1:NOM,1:2)
      DEALLOCATE (CWORK)
      DEALLOCATE (CWORKO)
C
C=================================================================
C     Calculate self-energy  on matsubara frequencies
C=================================================================
C
      ALLOCATE (SOM(NLM,NLM,NOM,NSPIN))
      SOM(1:NLM,1:NLM,1:NOM,1:2) = C0
C
      CALL DMFT_SPTFLEX(UUU,UJJ,TEMP,GOM,SOM,MINOM,OMEGA,DTIME,DOMEGA,
     &                  NSPIN,NLM,LOP,NLMSQ,NLMQU,NOM,DMFTDBLC,
     &                  IPRINTDMFT,IFLEX)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Set up the total energy correction, Galitskii-Migdal formula
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EMIGD = C0
      DO IS = 1,NSPIN
         DO LM = 1,NLM
            DO LMP = 1,NLM
               IF ( LM.EQ.LMP ) THEN
                  DO IOM = 1,NOM/2
                     EMIGD = EMIGD + TEMP*SOM(LM,LMP,2*IOM,IS)
     &                       *GOMIN(LM,LMP,2*IOM,IS)/2.0D0
                  END DO
               END IF
            END DO
         END DO
      END DO
C
      WRITE (6,99004)
      WRITE (6,99005) DREAL(EMIGD)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     print-out of the self energy
C
      IF ( IPRINTDMFT.GT.0 ) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE='sig_matsub_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sig_matsub_im.dat')
         DO IS = 1,NSPIN
            DO IOM = 1,NOM/2
               WRITE (IOTMP,'(11f20.10)') OMEGA(IOM*2),
     &                (DIMAG(SOM(LM,LM,IOM*2,IS)),LM=1,5)
            END DO
            WRITE (IOTMP,*) '##'
         END DO                 ! lz
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE='sig_matsub_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sig_matsub_re.dat')
         DO IS = 1,NSPIN
            DO IOM = 1,NOM/2
               WRITE (IOTMP,'(11f20.10)') OMEGA(IOM*2),
     &                (DREAL(SOM(LM,LM,IOM*2,IS)),LM=1,5)
            END DO
            WRITE (IOTMP,*) '##'
         END DO                 ! lz
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C============================================================
C     Create self energy on complex contour (pade approx)
C============================================================
      ALLOCATE (CWORK(NEGF))
      ALLOCATE (CWORKO(NOM))
      ALLOCATE (SEZNEW(NLM,NLM,NEGF,NSPIN))
      SEZNEW(1:NLM,1:NLM,1:NEGF,1:NSPIN) = C0
      CWORK(1:NEGF) = C0
      CWORKO(1:NOM) = C0
      DO IS = 1,NSPIN
         DO LM = 1,NLM
            DO LMP = 1,NLM
               IF ( IZERO(LMP,LM,IS).EQ.1 ) THEN
                  DO IOM = 1,NOM/2
                     CWORKO(IOM) = SOM(LMP,LM,2*IOM,IS)
                  END DO
                  CALL DMFT_PADEOZ(OMEGA,CWORKO,EFERMINEW,NEGF,EGFNEW,
     &                             CWORK,NOM)
                  SEZNEW(LMP,LM,1:NEGF,IS) = CWORK(1:NEGF)
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (CWORK)
C============================================================
C
C============================================================
C     Mixing of the self energy
C============================================================
      FAC1 = 1.D0 - DMFTMIX
      FAC2 = DMFTMIX
      IF ( ITRSCF.EQ.1 ) THEN
         FAC1 = 0.0D0
         FAC2 = 1.0D0
      END IF
C
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
                  DMFTSIGMA(IKM,IKMP,IE) = FAC1*SEZ(LM,LMP,IE,IS)
     &               + FAC2*SEZNEW(LM,LMP,IE,IS)
               END DO
            END DO
         END DO
      END DO
C============================================================
C============================================================
C     Create self energy on real axix (just for check)
C============================================================
      IF ( IPRINTDMFT.GT.0 ) THEN
         ALLOCATE (SEREAL(NLM,NLM,NEREAL,NSPIN))
         ALLOCATE (EREAL(NEREAL))
         ALLOCATE (CWORK(NEREAL))
         SEREAL(1:NLM,1:NLM,1:NEREAL,1:NSPIN) = C0
         CWORK(1:NEREAL) = C0
         CWORKO(1:NOM) = C0
         DE = (2.0D0+1.0D0)/DBLE(NEREAL-1)
         DO IE = 1,NEREAL
            RWORK = -1.0D0 + DE*DBLE(IE-1)
            EREAL(IE) = DCMPLX(RWORK,0.001D0)
         END DO
         DO IS = 1,NSPIN
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  IF ( IZERO(LMP,LM,IS).EQ.1 ) THEN
                     DO IOM = 1,NOM/2
                        CWORKO(IOM) = SOM(LMP,LM,IOM*2,IS)
                     END DO
                     CALL DMFT_PADEOZ(OMEGA,CWORKO,EFERMINEW,NEREAL,
     &                                EREAL,CWORK,NOM)
                     SEREAL(LMP,LM,1:NEREAL,IS) = CWORK(1:NEREAL)
                  END IF
               END DO
            END DO
         END DO
      END IF
C============================================================
C
C============================================================
C     print-out of the self energy
C============================================================
      FILNAM = TXT_T(1:LTXT_T)//'_SIGUNFORM'//'.sig'
      LFILNAM = LTXT_T + 10 + 4
C
      OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),STATUS='UNKNOWN',
     &      FORM='unformatted')
C
      REWIND IOTMP
      WRITE (IOTMP) NOM/2,DMFTTEMP
      WRITE (IOTMP) ((((SOM(LM,LMP,2*IOM,IS),LM=1,NLM),LMP=1,NLM),IOM=1,
     &              NOM/2),IS=1,NSPIN)
      WRITE (IOTMP) EFERMIOLD
      CALL DMFT_CLOSE_IOTMP(IOTMP)
C
C
      IF ( IPRINTDMFT.GT.0 ) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE='sig_egfnew_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sig_egfnew_im.dat')
         DO IE = 1,NEGF
            WRITE (IOTMP,'(30f20.10)') (EGFNEW(IE)-EFERMINEW)*RY_EV,
     &                                 (DIMAG(DMFTSIGMA(LMP,LMP,IE))
     &                                 *RY_EV,LMP=1,NKM)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE='sig_egfnew_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sig_egfnew_re.dat')
         DO IE = 1,NEGF
            WRITE (IOTMP,'(30f20.10)') (EGFNEW(IE)-EFERMINEW)*RY_EV,
     &                                 (DREAL(DMFTSIGMA(LMP,LMP,IE))
     &                                 *RY_EV,LMP=1,NKM)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C
      IF ( IPRINTDMFT.GT.0 ) THEN
         CALL DMFT_CLOSE_IOTMP(IOTMP)
C
         OPEN (UNIT=IOTMP,FILE='sigma_realz_im.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sigma_realz_im.dat')
         DO IE = 1,NEREAL
            WRITE (IOTMP,'(40f20.10)') (EREAL(IE)-EFERMINEW)*RY_EV,
     &                                 ((DIMAG(SEREAL(LMP,LMP,IE,IS))
     &                                 *RY_EV,LMP=1,NLM),IS=1,NSPIN)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
         OPEN (UNIT=IOTMP,FILE='sigma_realz_re.dat')
         CLOSE (IOTMP,STATUS='DELETE')
         OPEN (UNIT=IOTMP,FILE='sigma_realz_re.dat')
         DO IE = 1,NEREAL
            WRITE (IOTMP,'(40f20.10)') (EREAL(IE)-EFERMINEW)*RY_EV,
     &                                 ((DREAL(SEREAL(LMP,LMP,IE,IS))
     &                                 *RY_EV,LMP=1,NLM),IS=1,NSPIN)
         END DO
         CALL DMFT_CLOSE_IOTMP(IOTMP)
      END IF
C
99001 FORMAT (/,1x,'Non-zero elements of Gz for IS=',I2)
99002 FORMAT (1x,20(I2))
99003 FORMAT (1x,A)
99004 FORMAT (/,1X,'Migdal contributions to the total energy',/,5X,
     &        '        EMIGDAL')
99005 FORMAT (5X,5F14.6)
      END
