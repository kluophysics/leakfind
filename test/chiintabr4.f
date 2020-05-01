C*==chiintabr4.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTABR4(CSUM,AG,BG,WGTG,AGBG,AMEG,AF,BF,WGTF,AFBF,
     &                      AMEF,WGTR,NSOL,JTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *      INTEGRAL FROM 0 TO R_JTOP OF PRODUCT OF                     *
C   *        REGULAR/IRREGULAR SOLUTIONS OF DIRAC EQUATION             *
C   *        (E.G. Z(A)Z(B), Z(A)J(B), J(A)Z(B) UND J(A)J(B))          *
C   *        BY CALCULATING THE INTEGRALS G*G AND F*F                  *
C   *        AND PROPERLY SUMMING OVER INTERNAL COUPLINGS.             *
C   *      USED FOR CALCULATING DOUBLE INTEGRAL IN CHIRADINT WITH      *
C   *      INNER R-DEPENDENT INTEGRAL AS WEIGHT FUNCTION.              *
C   *                                                                  *
C   *      HF/MD: 12/96                                                *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHIINTABR418
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CSUM
      INTEGER JTOP,NRMAX,NSOL
      REAL*8 WGTF,WGTG
      COMPLEX*16 AF(NRMAX,2),AFBF(2,2,NRMAX),AG(NRMAX,2),AGBG(2,2,NRMAX)
     &           ,BF(NRMAX,2),BG(NRMAX,2),WGTR(NRMAX)
      REAL*8 AMEF(2,2),AMEG(2,2)
C
C*** End of declarations rewritten by SPAG
C
      CALL CHIINTABR2(AG,BG,AGBG,AF,BF,AFBF,WGTR,NSOL,NSOL,JTOP,NRMAX)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
      CALL SUMUPINT(CSUM,WGTG,AGBG(1,1,JTOP),AMEG,WGTF,AFBF(1,1,JTOP),
     &              AMEF,NSOL)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C
      END
