C*==dmft_ldaureadsig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_LDAUREADSIG(IOTMP,DMFTSIGMA,IT,NEGF,LOPT,NKM,
     &                            NKMMAX,NTMAX,EREFLDAU)
      IMPLICIT NONE
C*--DMFT_LDAUREADSIG5
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 EREFLDAU
      INTEGER IOTMP,IT,LOPT,NEGF,NKM,NKMMAX,NTMAX
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEGF)
C
C Local variables
C
      INTEGER I1,I2,I3,I4,IE,LM,LMP,NLM,NREP
      COMPLEX*16 RWORK(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE RWORK
C
      REWIND (IOTMP)
      READ (IOTMP) NREP
      ALLOCATE (RWORK(NREP,NREP))
      READ (IOTMP) ((RWORK(LM,LMP),LM=1,NREP),LMP=1,NREP)
      READ (IOTMP) EREFLDAU
C
      NLM = LOPT*2 + 1
      I1 = LOPT**2 + 1
      I2 = I1 + NLM - 1
      I3 = NKM/2 + I1
      I4 = I3 + NLM - 1
C
      DO IE = 1,NEGF
         DMFTSIGMA(I1:I2,I1:I2,IT,IE) = RWORK(1:NLM,1:NLM)
         DMFTSIGMA(I3:I4,I3:I4,IT,IE) = RWORK(NLM+1:NREP,NLM+1:NREP)
         DMFTSIGMA(I1:I2,I3:I4,IT,IE) = RWORK(1:NLM,NLM+1:NREP)
         DMFTSIGMA(I3:I4,I1:I2,IT,IE) = RWORK(NLM+1:NREP,1:NLM)
      END DO
C
      END