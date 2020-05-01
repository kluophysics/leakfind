C*==dmft_calcg0.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_CALCG0(GFMAT,DMFTSIGMA,NKM,NEGF,LOP,NKMMAX,GFMAT0,
     &                       NTMAX,IT)
C
C
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--DMFT_CALCG08
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
      INTEGER IT,LOP,NEGF,NKM,NKMMAX,NTMAX
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEGF),
     &           GFMAT(NKM,NKM,NEGF,NTMAX),GFMAT0(NKM,NKM,NEGF,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,CWORK1(:,:),CWORK2(:,:),GZ(:,:,:),SEZ(:,:,:,:)
      INTEGER I1,I2,I3,I4,IE,IR,IRP,IZERO(:,:),LM,LMP,NLM,NREP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CWORK1,CWORK2,GZ,SEZ,IZERO
C
      NLM = LOP*2 + 1
      NREP = 2*NLM
      ALLOCATE (SEZ(NREP,NREP,NEGF,NTMAX))
      ALLOCATE (GZ(NREP,NREP,NEGF))
      I1 = LOP**2 + 1
      I2 = I1 + NLM - 1
      I3 = NKM/2 + I1
      I4 = I3 + NLM - 1
C
      SEZ(1:NLM,1:NLM,1:NEGF,IT) = DMFTSIGMA(I1:I2,I1:I2,IT,1:NEGF)
      SEZ(NLM+1:NREP,NLM+1:NREP,1:NEGF,IT)
     &   = DMFTSIGMA(I3:I4,I3:I4,IT,1:NEGF)
      SEZ(1:NLM,NLM+1:NREP,1:NEGF,IT) = DMFTSIGMA(I1:I2,I3:I4,IT,1:NEGF)
      SEZ(NLM+1:NREP,1:NLM,1:NEGF,IT) = DMFTSIGMA(I3:I4,I1:I2,IT,1:NEGF)
C
      GZ(1:NLM,1:NLM,1:NEGF) = GFMAT(I1:I2,I1:I2,1:NEGF,IT)         !1_1
      GZ(NLM+1:NREP,NLM+1:NREP,1:NEGF) = GFMAT(I3:I4,I3:I4,1:NEGF,IT)
                                                                 !2_2
      GZ(1:NLM,NLM+1:NREP,1:NEGF) = GFMAT(I1:I2,I3:I4,1:NEGF,IT)    !1_2
      GZ(NLM+1:NREP,1:NLM,1:NEGF) = GFMAT(I3:I4,I1:I2,1:NEGF,IT)    !2_1
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
C============================================================
C     G_BATH=G^-1+SIGMA
C============================================================
C
C
      ALLOCATE (CWORK1(NREP,NREP),CWORK2(NREP,NREP))
C
      DO IE = 1,NEGF
         CWORK1(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,IE)
         CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
         CWORK1(1:NREP,1:NREP) = CWORK2(1:NREP,1:NREP)
     &                           + SEZ(1:NREP,1:NREP,IE,IT)
         CALL CMATINV(NREP,NREP,CWORK1,CWORK2)
         GZ(1:NREP,1:NREP,IE) = CWORK2(1:NREP,1:NREP)
      END DO
      DEALLOCATE (CWORK1,CWORK2)
C
      GFMAT0(I1:I2,I1:I2,1:NEGF,IT) = GZ(1:NLM,1:NLM,1:NEGF)
      GFMAT0(I3:I4,I3:I4,1:NEGF,IT) = GZ(NLM+1:NREP,NLM+1:NREP,1:NEGF)
      GFMAT0(I1:I2,I3:I4,1:NEGF,IT) = GZ(1:NLM,NLM+1:NREP,1:NEGF)
      GFMAT0(I3:I4,I1:I2,1:NEGF,IT) = GZ(NLM+1:NREP,1:NLM,1:NEGF)
      END
