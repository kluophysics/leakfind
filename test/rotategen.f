C*==rotategen.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE ROTATEGEN(T1,T2,N,ROT,UNITARY,M)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the rotated matrix                                    *
C   *                                                                  *
C   *                  T2 = R T1 R+                                    *
C   *                                                                  *
C   *  allow for anti-unitary symmetry operations R                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--ROTATEGEN15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ROTATEGEN')
C
C Dummy arguments
C
      INTEGER M,N
      LOGICAL UNITARY
      COMPLEX*16 ROT(M,M),T1(M,M),T2(M,M)
C
C Local variables
C
      CHARACTER*1 CNT
      INTEGER IA_ERR
      COMPLEX*16 W1(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1
C
      ALLOCATE (W1(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: W1')
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
      IF ( UNITARY ) THEN
         CNT = 'N'
      ELSE
         CNT = 'T'
      END IF
C
      CALL ZGEMM('N',CNT,N,N,N,C1,ROT,M,T1,M,C0,W1,M)
C
      CALL ZGEMM('N','C',N,N,N,C1,W1,M,ROT,M,C1,T2,M)
C
      DEALLOCATE (W1,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
      END
