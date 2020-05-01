C*==sigmanipul.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGMANIPUL(MREL,KEY)
C   ********************************************************************
C   *                                                                  *
C   *     manipulate the matrix  MREL  in (kappa,mue)-representation   *
C   *                                                                  *
C   *     KEY = 1:      set spin-flip blocks to 0                      *
C   *     KEY = 2:  1 + set spin down block  to 0                      *
C   *     KEY = 3:  1 + set spin up   block  to 0                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,WKM1
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--SIGMANIPUL16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KEY
      COMPLEX*16 MREL(NKMMAX,NKMMAX)
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C
      CALL CHANGEREP(NKM,NKMMAX,MREL,'REL>RLM',WKM1)
C
C---------------------------------------------- suppress spin-flip terms
      WKM1(1:NLM,NLM+1:NKM) = C0
      WKM1(NLM+1:NKM,1:NLM) = C0
C
      IF ( KEY.EQ.2 ) THEN
C---------------------------------------------- suppress spin down terms
C
         WKM1(1:NLM,1:NLM) = C0
C
      ELSE IF ( KEY.EQ.3 ) THEN
C---------------------------------------------- suppress spin up   terms
C
         WKM1(NLM+1:NKM,NLM+1:NKM) = C0
C
      END IF
C
      CALL CHANGEREP(NKM,NKMMAX,WKM1,'RLM>REL',MREL)
C
      END
