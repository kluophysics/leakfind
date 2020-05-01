C*==specmain.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SPECMAIN
C
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine program for       KKRSPEC                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:MEZJ,MEZZ,MSSQ,SSST,TAUQ,TAUT
      USE MOD_ENERGY,ONLY:PHASK
      IMPLICIT NONE
C*--SPECMAIN13
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,99001)
C
      CALL TRANSPHO_RSLAB(MEZJ,MEZZ,MSSQ,SSST,TAUQ,TAUT,PHASK)
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*              ****   *****   *****   ******                 *'
     &  ,/,10X,
     &  '*             *    *  *    *  *       *                      *'
     &  ,/,10X,
     &  '*             *       *    *  *       *                      *'
     &  ,/,10X,
     &  '*              ****   *****   *****   *                      *'
     &  ,/,10X,
     &  '*                  *  *       *       *                      *'
     &  ,/,10X,
     &  '*             *    *  *       *       *                      *'
     &  ,/,10X,
     &  '*              ****   *       ******  ******                 *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
      END
