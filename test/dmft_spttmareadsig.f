C*==dmft_spttmareadsig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_SPTTMAREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,LOP,NKM,
     &                              NKMMAX,NEMAX,NTMAX,EREFDMFT)
C
C*****************************************************************
C     CIRCLEMESH:  semicircle mesh, read from unformatted file
C     SIGMA:   self-energy on the semicircle, read from
C              unformatted file
C     CPADE:       Pade coefficients
C     TMA_STATIC:  logical variable, read from unformatted file:
C                  if .TRUE., Im(SIGMA)=0, Re(Sigma)=const(energy)
C*****************************************************************
C
C
      USE MOD_ENERGY,ONLY:EFERMI
      IMPLICIT NONE
C*--DMFT_SPTTMAREADSIG17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 EREFDMFT
      INTEGER IOTMP,IT,LOP,NEGF,NEMAX,NKM,NKMMAX,NTMAX
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),EGF(NEMAX)
C
C Local variables
C
      COMPLEX*16 CIRCLEMESH(:),CPADE(:,:),CTMP,CWORK(:),CWORK1(:),
     &           SE_AR(:,:,:),SIGMA(:,:,:)
      INTEGER I1,I2,IE,IR,IRP,IZERO(:,:),NREP,NRR
      REAL*8 SIGSTAT(:,:),XX(:)
      LOGICAL TMA_STATIC,WRSIGKK
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATABLE CPADE,SIGMA,IZERO,CIRCLEMESH,SIGSTAT,CWORK,XX,SE_AR
      ALLOCATABLE CWORK1
C
      WRSIGKK = .FALSE.
C
      READ (IOTMP) NRR,NREP,TMA_STATIC
C
C
      ALLOCATE (CIRCLEMESH(NRR),SIGMA(NREP,NREP,NRR))
      ALLOCATE (IZERO(NREP,NREP),SIGSTAT(NREP,NREP))
      ALLOCATE (CPADE(NRR,NRR),CWORK(NRR))
C
      READ (IOTMP) CIRCLEMESH,SIGMA,IZERO
      WRITE (6,*) SIGMA
      READ (IOTMP,ERR=100) EREFDMFT
      READ (IOTMP,ERR=200) WRSIGKK
C
      GOTO 200
 100  CONTINUE
      EREFDMFT = DCMPLX(0.7D0,0.0D0)
C
C
C      WRSIGKK=.FALSE.
 200  CONTINUE
      IF ( .NOT.WRSIGKK ) THEN
         DO IR = 1,NREP
            SIGSTAT(IR,IR) = DREAL(SIGMA(IR,IR,NRR))
         END DO
C
C
         DO IR = 1,NREP
            DO IRP = IR,NREP
               IF ( IR.LE.NREP/2 ) I1 = LOP**2 + IR
               IF ( IR.GT.NREP/2 ) I1 = NKM/2 + LOP**2 + IR - NREP/2
               IF ( IRP.LE.NREP/2 ) I2 = LOP**2 + IRP
               IF ( IRP.GT.NREP/2 ) I2 = NKM/2 + LOP**2 + IRP - NREP/2
               IF ( TMA_STATIC ) THEN
                  DMFTSIGMA(I1,I2,IT,1:NEGF) = SIGMA(IRP,IR,NRR)
                  DMFTSIGMA(I2,I1,IT,1:NEGF) = SIGMA(IR,IRP,NRR)
               ELSE
                  CTMP = (0D0,0D0)
                  IF ( IZERO(IRP,IR).EQ.1 ) THEN
C make Pade from semicircle
C                  IF ( IR.EQ.IRP ) SIGMA(IR,IR,1:NRR)
C     &                 = SIGMA(IR,IR,1:NRR) - SIGSTAT(IR,IR)
                     CWORK(1:NRR) = SIGMA(IRP,IR,1:NRR)
                     CALL PADECOEFF(CWORK,CIRCLEMESH,CPADE,NRR,NRR-1)
C                 NRR-1: avoid SIGMA(Ef)=0
                     DO IE = 1,NEGF
                        CALL PADEAPPROX(CTMP,EGF(IE),CIRCLEMESH,CPADE,
     &                                  NRR,NRR)
C                     IF ( IR.EQ.IRP ) CTMP = CTMP + SIGSTAT(IR,IR)
                        DMFTSIGMA(I1,I2,IT,IE) = CTMP
                        DMFTSIGMA(I2,I1,IT,IE) = CTMP
                     END DO
                  END IF
C IZERO
               END IF
C TMA_STATIC
            END DO
         END DO
      ELSE
         READ (IOTMP) SIGSTAT
         READ (IOTMP) NRR
         ALLOCATE (XX(NRR))
         READ (IOTMP) XX
         ALLOCATE (SE_AR(NRR,NREP,NREP))
         READ (IOTMP) SE_AR
C
         ALLOCATE (CWORK1(NEGF))
         CWORK1(1:NEGF) = (0D0,0D0)
         WRITE (6,*) 'KK will be used to read in SIG'
         XX(1:NRR) = XX(1:NRR) + EFERMI
         DO IR = 1,NREP
            DO IRP = IR,NREP
               IF ( IR.LE.NREP/2 ) I1 = LOP**2 + IR
               IF ( IR.GT.NREP/2 ) I1 = NKM/2 + LOP**2 + IR - NREP/2
               IF ( IRP.LE.NREP/2 ) I2 = LOP**2 + IRP
               IF ( IRP.GT.NREP/2 ) I2 = NKM/2 + LOP**2 + IRP - NREP/2
C
               IF ( IZERO(IRP,IR).EQ.1 ) THEN
C
                  CALL DMFT_CAUCHYT(NRR,XX,SE_AR(1,IRP,IR),NEGF,EGF,
     &                              CWORK1)
                  CWORK1(1:NEGF) = CWORK1(1:NEGF) - CWORK1(NEGF)
C
                  DMFTSIGMA(I1,I2,IT,1:NEGF) = CWORK1(1:NEGF)
     &               + SIGSTAT(IR,IRP)
                  DMFTSIGMA(I2,I1,IT,1:NEGF) = CWORK1(1:NEGF)
     &               + SIGSTAT(IRP,IR)
               END IF
            END DO
         END DO
      END IF
C
C
      DEALLOCATE (SIGMA,IZERO,SIGSTAT,CIRCLEMESH,CPADE,CWORK)
C
C
      END
