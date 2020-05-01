C*==types_swap.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TYPES_SWAP(NT0,IT0_TAUX)
C   ********************************************************************
C   *                                                                  *
C   *   swap type dependent information if symmetry has changed        *
C   *                                                                  *
C   *   NOTE: the symmetry programm has split group of equivalent      *
C   *         atom types leading to NT instead of NT0 types            *
C   *         the new types have been added at the end of the          *
C   *         tables -- to avoid conflicts in case of TB or cluster    *
C   *         calculations the new types are inserted into the table   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT,NQ
      USE MOD_RMESH,ONLY:JRWS,NRMAX
      USE MOD_TYPES,ONLY:NT,CONC,IMT,NVALT,Z,RHOCHR,RHOCHRC,RHOSPN,
     &    RHOSPNC,VT,BT,AOPT,NAT
      USE MOD_CALCMODE,ONLY:BLCOUPL,ORBPOL
      IMPLICIT NONE
C*--TYPES_SWAP21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NT0
      INTEGER IT0_TAUX(NT)
C
C Local variables
C
      REAL*8 AOP_TAUX(:,:,:),B_TAUX(:,:),RHOCHRC_TAUX(:,:),
     &       RHOCHR_TAUX(:,:),RHOSPNC_TAUX(:,:),RHOSPN_TAUX(:,:),
     &       V_TAUX(:,:),X_TAUX(:)
      INTEGER IA,IM_TAUX(:),IQ_ATAUX(:,:),IT,IT0,JT,NA_TAUX(:),
     &        NVAL_TAUX(:),NX,Z_TAUX(:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE IQ_ATAUX,NA_TAUX
      ALLOCATABLE X_TAUX,Z_TAUX,IM_TAUX,NVAL_TAUX
      ALLOCATABLE V_TAUX,RHOCHR_TAUX,RHOCHRC_TAUX
      ALLOCATABLE B_TAUX,RHOSPN_TAUX,RHOSPNC_TAUX,AOP_TAUX
C
      WRITE (6,99001)
C
      ALLOCATE (X_TAUX(NT),Z_TAUX(NT),IM_TAUX(NT),NVAL_TAUX(NT))
      X_TAUX(1:NT) = CONC(1:NT)
      Z_TAUX(1:NT) = Z(1:NT)
      IM_TAUX(1:NT) = IMT(1:NT)
      NVAL_TAUX(1:NT) = NVALT(1:NT)
C
      ALLOCATE (IQ_ATAUX(NQ,NT),NA_TAUX(NT))
      NA_TAUX(1:NT) = NAT(1:NT)
      IQ_ATAUX(1:NQ,1:NT) = IQAT(1:NQ,1:NT)
C
      ALLOCATE (V_TAUX(NRMAX,NT),B_TAUX(NRMAX,NT))
      ALLOCATE (RHOCHR_TAUX(NRMAX,NT),RHOCHRC_TAUX(NRMAX,NT))
      ALLOCATE (RHOSPN_TAUX(NRMAX,NT),RHOSPNC_TAUX(NRMAX,NT))
C
      V_TAUX(1:NRMAX,1:NT) = VT(1:NRMAX,1:NT)
      B_TAUX(1:NRMAX,1:NT) = BT(1:NRMAX,1:NT)
      RHOCHR_TAUX(1:NRMAX,1:NT) = RHOCHR(1:NRMAX,1:NT)
      RHOSPN_TAUX(1:NRMAX,1:NT) = RHOSPN(1:NRMAX,1:NT)
      RHOCHRC_TAUX(1:NRMAX,1:NT) = RHOCHRC(1:NRMAX,1:NT)
      RHOSPNC_TAUX(1:NRMAX,1:NT) = RHOSPNC(1:NRMAX,1:NT)
C
      IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL ) THEN
         ALLOCATE (AOP_TAUX(NRMAX,1:2,NT))
         AOP_TAUX(1:NRMAX,1:2,1:NT) = AOPT(1:NRMAX,1:2,1:NT)
      END IF
C
      JT = 0
C=======================================================================
      DO IT0 = 1,NT0
C
         DO IT = 1,NT
            IF ( IT0_TAUX(IT).EQ.IT0 ) THEN
               JT = JT + 1
C
C-----------------------------------------------------------------------
C  the occupation parameters NAT and IQ_ATAUX were already updated by
C  symmetry program -- adjust to new numbering
               NAT(JT) = NA_TAUX(IT)
               DO IA = 1,NAT(JT)
                  IQAT(IA,JT) = IQ_ATAUX(IA,IT)
               END DO
C-----------------------------------------------------------------------
C
               Z(JT) = Z_TAUX(IT0)
               CONC(JT) = X_TAUX(IT0)
               IMT(JT) = IM_TAUX(IT0)
               NVALT(JT) = NVAL_TAUX(IT0)
C
               NX = JRWS(IMT(JT))
C
               VT(1:NX,JT) = V_TAUX(1:NX,IT0)
               BT(1:NX,JT) = B_TAUX(1:NX,IT0)
               RHOCHR(1:NX,JT) = RHOCHR_TAUX(1:NX,IT0)
               RHOSPN(1:NX,JT) = RHOSPN_TAUX(1:NX,IT0)
               RHOCHRC(1:NX,JT) = RHOCHRC_TAUX(1:NX,IT0)
               RHOSPNC(1:NX,JT) = RHOSPNC_TAUX(1:NX,IT0)
               IF ( ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL )
     &              AOPT(1:NX,1:2,JT) = AOP_TAUX(1:NX,1:2,IT0)
            END IF
C
         END DO
      END DO
C=======================================================================
C
99001 FORMAT (/,1X,79('*'),/,34X,'<TYPES_SWAP>',/,1X,79('*'),//,10X,
     &        'the symmetry is reduced compared to potential file',/,
     &        10X,'IT-dependent tables will be updated   ',/)
      END
