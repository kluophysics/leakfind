C*==scfmad2d_ab.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD2D_AB
C   ********************************************************************
C   *                                                                  *
C   *    prepare caucltaion of fix the Madelung parameters             *
C   *    VMAD2D_A and VMAD2D_B  for 2D systems                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ_L,NQ_I,QBAS,AVMAD_LI,AVMAD_RI,CMNTQ,NLMQMAD,
     &    VMAD2D_A,VMAD2D_B,IQ_LREF,IQ_RREF
      USE MOD_SITES,ONLY:AVMAD_LR,AVMAD_LL,AVMAD_RR,AVMAD_RL,AVMAD,NQ_R
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      IMPLICIT NONE
C*--SCFMAD2D_AB17
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      LOGICAL INITIALIZE
      INTEGER IQ,IQ_R1,IQ_R2,JQ,LM
      REAL*8 V00MAD_LB,V00MAD_LL,V00MAD_LR,V00MAD_RB,V00MAD_RL,V00MAD_RR
      SAVE V00MAD_LB,V00MAD_LL,V00MAD_LR,V00MAD_RB,V00MAD_RL,V00MAD_RR
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      IQ_R1 = NQ_L + NQ_I + 1
      IQ_R2 = NQ_L + NQ_I + NQ_R
C
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
      IF ( INITIALIZE ) THEN
C
C----------------------------- Madelung contribution to potential:  V^LB
C
         IQ = IQ_LREF
         V00MAD_LB = 0.0D0
         V00MAD_LL = 0.0D0
         DO JQ = 1,NQ_L
            DO LM = 1,NLMQMAD
               V00MAD_LB = V00MAD_LB + AVMAD(IQ,JQ,1,LM)*CMNTQ(LM,JQ)
               V00MAD_LL = V00MAD_LL + AVMAD_LL(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
         V00MAD_LR = 0.0D0
         DO JQ = IQ_R1,IQ_R2
            DO LM = 1,NLMQMAD
               V00MAD_LR = V00MAD_LR + AVMAD_LR(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
C----------------------------- Madelung contribution to potential:  V^RB
C
         IQ = IQ_R1
         V00MAD_RB = 0.0D0
         DO JQ = IQ_R1,IQ_R2
            DO LM = 1,NLMQMAD
               WRITE (*,*) AVMAD(IQ,JQ,1,LM),CMNTQ(LM,JQ)
               V00MAD_RB = V00MAD_RB + AVMAD(IQ,JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
         IQ = IQ_RREF
         V00MAD_RR = 0.0D0
         DO JQ = IQ_R1,IQ_R2
            DO LM = 1,NLMQMAD
               V00MAD_RR = V00MAD_RR + AVMAD_RR(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
         V00MAD_RL = 0.0D0
         DO JQ = 1,NQ_L
            DO LM = 1,NLMQMAD
               V00MAD_RL = V00MAD_RL + AVMAD_RL(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
C--------------- L- and R-bulk related Madelung constants no more needed
         DEALLOCATE (AVMAD_LL,AVMAD_LR,AVMAD_RL,AVMAD_RR)
C
         INITIALIZE = .FALSE.
C
      END IF
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
C
C------------------------------------------------- Calculation of B term
C
      VMAD2D_B = V00MAD_LB - V00MAD_LL - V00MAD_LR
C
      DO JQ = NQ_L + 1,NQ_L + NQ_I
         DO LM = 1,NLMQMAD
            VMAD2D_B = VMAD2D_B - AVMAD_LI(JQ,1,LM)*CMNTQ(LM,JQ)
         END DO
      END DO
C
      VMAD2D_B = VMAD2D_B/SQRT_4PI
C
C------------------------------------------------- Calculation of A term
      IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
C------------- at the boundary to vacuum d V / dr  should be zero
C
         VMAD2D_A = 0.0D0
C
      ELSE
C
         VMAD2D_A = V00MAD_RB - (V00MAD_RL+V00MAD_RR)
C
         DO JQ = NQ_L + 1,NQ_L + NQ_I
            DO LM = 1,NLMQMAD
               VMAD2D_A = VMAD2D_A - AVMAD_RI(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
         VMAD2D_A = VMAD2D_A - V00MAD_LB + (V00MAD_LL+V00MAD_LR)
C
         DO JQ = NQ_L + 1,NQ_L + NQ_I
            DO LM = 1,NLMQMAD
               VMAD2D_A = VMAD2D_A + AVMAD_LI(JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
C
         VMAD2D_A = VMAD2D_A/(QBAS(3,NQ_L+NQ_I+1)-QBAS(3,NQ_L))/SQRT_4PI
      END IF
C
      IF ( IPRINT.GE.0 ) THEN
         WRITE (6,99001)
         WRITE (6,99002) V00MAD_LB,V00MAD_RB,VMAD2D_A,VMAD2D_B
      END IF
C
99001 FORMAT (/,1X,'VMAD contributions for LIV or VIV              ',/,
     &        5X,'V00MAD_LB V00MAD_RB VMAD2D_A   VMAD2D_B')
99002 FORMAT (2X,4F10.5)
C
      END
