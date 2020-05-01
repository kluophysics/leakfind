C*==dmft_3bsreadsig.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DMFT_3BSREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,LOP,NKM,
     &                           NKMMAX,NEMAX,NTMAX,EFERMI)
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
      IMPLICIT NONE
C*--DMFT_3BSREADSIG15
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EFERMI
      INTEGER IOTMP,IT,LOP,NEGF,NEMAX,NKM,NKMMAX,NTMAX
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),EGF(NEMAX)
C
C Local variables
C
      REAL*8 CTMP(:),EIN(:),RTMP,RTMP2
      COMPLEX*16 CWORK(:,:,:)
      INTEGER I1,I2,IE,IR,IRP,NEIN,NREP
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CWORK,EIN,CTMP
C
      READ (IOTMP) NEIN,NREP
C
      ALLOCATE (EIN(NEIN),CWORK(NREP,NREP,NEIN),CTMP(NEIN))
      READ (IOTMP) EIN,CWORK
      EIN(:) = EIN(:) + EFERMI
C
C      ALLOCATE (ROT(NLM3BS,NLM3BS))
C      ALLOCATE (CROT(5,5))
C      CALL ROTMAT_RLM(NL3BS,45.0D0,90.0D0,45.0D0,ROT,NLM3BS)
CC
C      DO IR = 1,NLM3BS
C         DO IRP = 1,NLM3BS
C            W1(IRP,IR) = DCMPLX(ROT(IRP,IR),0.0D0)
C         END DO
C      END DO
C      CROT(1:5,1:5) = W1(5:9,5:9)
C      ALLOCATE (W2(5,5))
C      ALLOCATE (W3(5,5))
C      DO IE = 1,NEIN
C         W2(:,:) = CWORK(1:5,1:5,IE)
C         CALL ZGEMM('N','N',5,5,5,C1,CROT,5,W2,5,C0,W3,5)
CC
C         CALL ZGEMM('N','C',5,5,5,C1,W3,5,CROT,5,C0,W2,5)
C         CWORK(1:5,1:5,IE) = W2
C         W2(:,:) = CWORK(6:10,6:10,IE)
C         CALL ZGEMM('N','N',5,5,5,C1,CROT,5,W2,5,C0,W3,5)
CC
C         CALL ZGEMM('N','C',5,5,5,C1,W3,5,CROT,5,C0,W2,5)
C         CWORK(6:10,6:10,IE) = W2
C      END DO
      DO IE = 1,NEGF
         DO IR = 1,NREP
            DO IRP = 1,NREP
               IF ( IR.LE.NREP/2 ) I1 = LOP**2 + IR
               IF ( IR.GT.NREP/2 ) I1 = NKM/2 + LOP**2 + IR - NREP/2
               IF ( IRP.LE.NREP/2 ) I2 = LOP**2 + IRP
               IF ( IRP.GT.NREP/2 ) I2 = NKM/2 + LOP**2 + IRP - NREP/2
               IF ( DREAL(EGF(IE)).LT.EIN(1) .OR. DREAL(EGF(IE))
     &              .GT.EIN(NEIN) ) THEN
                  DMFTSIGMA(I1,I2,IT,IE) = DCMPLX(0.0D0,0.0D0)
               ELSE
                  CTMP(1:NEIN) = DREAL(CWORK(IR,IRP,1:NEIN))
                  RTMP = YLAG(DREAL(EGF(IE)),EIN,CTMP,0,3,NEIN)
                  CTMP(1:NEIN) = DIMAG(CWORK(IR,IRP,1:NEIN))
                  RTMP2 = YLAG(DREAL(EGF(IE)),EIN,CTMP,0,3,NEIN)
                  DMFTSIGMA(I1,I2,IT,IE) = DCMPLX(RTMP,RTMP2)
               END IF
            END DO
         END DO
      END DO
C
      END
