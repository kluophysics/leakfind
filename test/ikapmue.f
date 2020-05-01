C*==ikapmue.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      FUNCTION IKAPMUE(KAPPA,MUEM05)
C   ********************************************************************
C   *                                                                  *
C   *  INDEXING OF MATRIX-ELEMENTS:                                    *
C   *                                                                  *
C   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--IKAPMUE12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAPPA,MUEM05
      INTEGER IKAPMUE
C
C Local variables
C
      INTEGER JP05,L
C
C*** End of declarations rewritten by SPAG
C
      JP05 = IABS(KAPPA)
C
      IF ( KAPPA.LT.0 ) THEN
         L = -KAPPA - 1
      ELSE
         L = KAPPA
      END IF
C
      IKAPMUE = 2*L*JP05 + JP05 + MUEM05 + 1
C
      END
C*==fun_l_lm.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      FUNCTION FUN_L_LM(LM)
C   ********************************************************************
C   *                                                                  *
C   *  return  l  for given combined index  LM                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--FUN_L_LM58
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LM
      INTEGER FUN_L_LM
C
C Local variables
C
      INTEGER I,L,M
C
C*** End of declarations rewritten by SPAG
C
      I = 0
      DO L = 0,100
         DO M = -L,L
            I = I + 1
            IF ( I.EQ.LM ) THEN
               FUN_L_LM = L
               RETURN
            END IF
         END DO
      END DO
C
      WRITE (6,*) 'ERROR in <FUN_L_LM>:  l not found for LM:',LM
C
      STOP
C
      END
C*==fun_m_lm.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      FUNCTION FUN_M_LM(LM)
C   ********************************************************************
C   *                                                                  *
C   *  return  l  for given combined index  LM                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--FUN_M_LM109
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LM
      INTEGER FUN_M_LM
C
C Local variables
C
      INTEGER I,L,M
C
C*** End of declarations rewritten by SPAG
C
      I = 0
      DO L = 0,100
         DO M = -L,L
            I = I + 1
            IF ( I.EQ.LM ) THEN
               FUN_M_LM = M
               RETURN
            END IF
         END DO
      END DO
C
      WRITE (6,*) 'ERROR in <FUN_M_LM>:  l not found for LM:',LM
C
      STOP
C
      END
