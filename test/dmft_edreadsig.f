C*==dmft_edreadsig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_EDREADSIG(IOTMP,DMFTSIGMA,EGF,IT,NEGF,LOP,NKM,
     &                          NKMMAX,NEMAX,NTMAX,EFERMI)
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
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C*--DMFT_EDREADSIG17
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
C
C
      ALLOCATABLE CWORK,EIN,CTMP
C
C
      READ (IOTMP) NEIN,NREP
      ALLOCATE (EIN(NEIN),CWORK(NREP,NREP,NEIN),CTMP(NEIN))
      READ (IOTMP) EIN
      READ (IOTMP) CWORK
      EIN(:) = EIN(:)/RY_EV + EFERMI
C
      DO IE = 1,NEGF
         DO IR = 1,NREP
            DO IRP = 1,NREP
               IF ( IR.LE.NREP/2 ) I1 = LOP**2 + IR
               IF ( IR.GT.NREP/2 ) I1 = NKM/2 + LOP**2 + IR - NREP/2
               IF ( IRP.LE.NREP/2 ) I2 = LOP**2 + IRP
               IF ( IRP.GT.NREP/2 ) I2 = NKM/2 + LOP**2 + IRP - NREP/2
               IF ( IR.EQ.IRP ) THEN
                  IF ( DREAL(EGF(IE)).LT.EIN(1) .OR. DREAL(EGF(IE))
     &                 .GT.EIN(NEIN) ) THEN
                     DMFTSIGMA(I1,I2,IT,IE) = DCMPLX(0.0D0,0.0D0)
                  ELSE
                     CTMP(1:NEIN) = DREAL(CWORK(IR,IRP,1:NEIN))
                     RTMP = YLAG(DREAL(EGF(IE)),EIN,CTMP,0,3,NEIN)
                     CTMP(1:NEIN) = DIMAG(CWORK(IR,IRP,1:NEIN))
                     RTMP2 = YLAG(DREAL(EGF(IE)),EIN,CTMP,0,3,NEIN)
                     DMFTSIGMA(I1,I2,IT,IE) = DCMPLX(RTMP,RTMP2)
                  END IF
               END IF
            END DO
         END DO
      END DO
C
      END
