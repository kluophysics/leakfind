C*==chiinthvv3.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTHVV3(IM,CIAB,AG,BG,WGTG,AMEG,AF,BF,WGTF,AMEF,
     &                      WGTR,NSOL,IRTOP)
C   ********************************************************************
C   *                                                                  *
C   *      CUMULATIVE INTEGRAL FROM 0 TO R_I OF PRODUCT OF             *
C   *        REGULAR/IRREGULAR SOLUTIONS OF DIRAC EQUATION             *
C   *        (E.G. Z(A)Z(B), Z(A)J(B), J(A)Z(B) UND J(A)J(B))          *
C   *        BY CALCULATING THE INTEGRALS G*G AND F*F                  *
C   *        AND PROPERLY SUMMING OVER INTERNAL COUPLINGS.             *
C   *                                                                  *
C   *      HF/MD: 12/96                                                *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CHIINTHVV317
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,NSOL
      REAL*8 WGTF,WGTG
      COMPLEX*16 AF(NRMAX,2,2),AG(NRMAX,2,2),BF(NRMAX,2,2),BG(NRMAX,2,2)
     &           ,CIAB(2,2,NRMAX),WGTR(NRMAX)
      REAL*8 AMEF(2,2),AMEG(2,2)
C
C Local variables
C
      COMPLEX*16 AFBF(2,2),AGBG(2,2),CINT(NRMAX),CSUM
      INTEGER K1,K2
C
C*** End of declarations rewritten by SPAG
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
C
            CALL CHIINTHVV2(IM,AG(1,1,K1),BG(1,1,K2),AGBG,NSOL,NSOL,
     &                      IRTOP,CINT,WGTR)
            CALL CHIINTHVV2(IM,AF(1,1,K1),BF(1,1,K2),AFBF,NSOL,NSOL,
     &                      IRTOP,CINT,WGTR)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
            CALL SUMUPINT(CSUM,WGTG,AGBG,AMEG,WGTF,AFBF,AMEF,NSOL)
C
            CIAB(K1,K2,IRTOP) = CSUM
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
         END DO
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      END
