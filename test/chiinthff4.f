C*==chiinthff4.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTHFF4(IM,CSUM,AG,BG,AF,BF,WGTHF,RMEHF,AMEHF,WGTR,
     &                      NSOL,IRTOP)
C   ********************************************************************
C   *                                                                  *
C   *      INTEGRAL FROM 0 TO R_JTOP OF PRODUCT OF                     *
C   *        REGULAR/IRREGULAR SOLUTIONS OF DIRAC EQUATION             *
C   *        (E.G. Z(A)Z(B), Z(A)J(B), J(A)Z(B) UND J(A)J(B))          *
C   *        BY CALCULATING THE INTEGRALS G*G AND F*F                  *
C   *        AND PROPERLY SUMMING OVER INTERNAL COUPLINGS.             *
C   *      USED FOR CALCULATING DOUBLE INTEGRAL IN AMEHFIRADINT WITH   *
C   *      INNER R-DEPENDENT INTEGRAL AS WEIGHT FUNCTION.              *
C   *                                                                  *
C   *      HF/MD: 12/96                                                *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CHIINTHFF420
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CSUM
      INTEGER IM,IRTOP,NSOL
      REAL*8 WGTHF
      COMPLEX*16 AF(NRMAX,2),AG(NRMAX,2),BF(NRMAX,2),BG(NRMAX,2),
     &           RMEHF(NCPLWFMAX,NCPLWFMAX),WGTR(NRMAX)
      REAL*8 AMEHF(2,2)
C
C Local variables
C
      COMPLEX*16 CINT(NRMAX)
C
C*** End of declarations rewritten by SPAG
C
      CALL CHIINTHFF2(IM,AG,AF,BG,BF,RMEHF,NSOL,NSOL,IRTOP,CINT,WGTR)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
      CALL SUMUPINT(CSUM,WGTHF,RMEHF,AMEHF,0.0D0,RMEHF,AMEHF,NSOL)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C
      END
