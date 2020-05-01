C*==dmft_matsub.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_MATSUB(DTIME,DOMEGA,TIME,OMEGA,MINOM,NOM)
C
      IMPLICIT NONE
C*--DMFT_MATSUB5
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DOMEGA,DTIME
      INTEGER NOM
      INTEGER MINOM(NOM)
      REAL*8 OMEGA(NOM),TIME(NOM)
C
C Local variables
C
      INTEGER IO,IOM
C
C*** End of declarations rewritten by SPAG
C
C
C===========================================================
C      domega=temp*pi
C      dtime=1.d0/nom/temp*2.d0
C      Ordering of the "periodic" Matsubara frequencies:
C      0,piT,2piT,..... Nom/2piT,-(Nom/2-1),....,-piT
C===========================================================
      DO IOM = 1,NOM
         IO = IOM - 1
         IF ( IOM.GT.NOM/2 ) IO = IOM - NOM - 1
         OMEGA(IOM) = DOMEGA*IO
         TIME(IOM) = DTIME*IO
      END DO
C===========================================================
C------- Define "minus" array for Tau
C------- if i=0,...,N-1 then  -i=N-i
C------- if i=1,...,N   then  -i=N+2-i
C===========================================================
      MINOM(1) = 1
      DO IOM = 2,NOM
         MINOM(IOM) = NOM + 2 - IOM
      END DO
C
      END
