C*==strsmat.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE STRSMAT(SRREL,NRREL,IRREL)
C   ********************************************************************
C   *                                                                  *
C   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
C   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
C   *                                                                  *
C   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:RREL,NKMMAX,NLM,NKM
      IMPLICIT NONE
C*--STRSMAT14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
      COMPLEX*16 SRREL(2,2,NKMMAX)
C
C Local variables
C
      INTEGER LAM,LR,NS1,NS2
C
C*** End of declarations rewritten by SPAG
C
C ----------------------------------------------------------------------
C RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
C ----------------------------------------------------------------------
C
C     ---------------------------------------------------
C     store the elements of  RREL
C     ---------------------------------------------------
C
      DO LAM = 1,NKM
         NS1 = 0
         NS2 = 0
C
         DO LR = 1,2*NLM
            IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
               IF ( LR.LE.NLM ) THEN
                  NS1 = NS1 + 1
                  IF ( NS1.GT.2 ) STOP ' IN <STRSMAT>   NS1 > 2'
                  SRREL(NS1,1,LAM) = RREL(LR,LAM)
                  IRREL(NS1,1,LAM) = LR
               ELSE
                  NS2 = NS2 + 1
                  IF ( NS2.GT.2 ) STOP ' IN <STRSMAT>   NS2 > 2'
                  SRREL(NS2,2,LAM) = RREL(LR,LAM)
                  IRREL(NS2,2,LAM) = LR - NLM
               END IF
            END IF
         END DO
C
         NRREL(1,LAM) = NS1
         NRREL(2,LAM) = NS2
      END DO
C
      END
