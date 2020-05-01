C*==chiintabr3in.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTABR3IN(CIAB,AG,BG,WGTG,AGBG,AMEG,AF,BF,WGTF,AFBF,
     &                        AMEF,WGTR,NSOL,JTOP,NRMAX)
C
C
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
      IMPLICIT NONE
C*--CHIINTABR3IN18
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NRMAX,NSOL
      REAL*8 WGTF,WGTG
      COMPLEX*16 AF(NRMAX,2,2),AFBF(2,2,NRMAX),AG(NRMAX,2,2),
     &           AGBG(2,2,NRMAX),BF(NRMAX,2,2),BG(NRMAX,2,2),
     &           CIAB(2,2,NRMAX),WGTR(NRMAX)
      REAL*8 AMEF(2,2),AMEG(2,2)
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER I,K1,K2
C
C*** End of declarations rewritten by SPAG
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
C
            CALL CHIINTABR2IN(AG(1,1,K1),BG(1,1,K2),AGBG,AF(1,1,K1),
     &                        BF(1,1,K2),AFBF,WGTR,NSOL,NSOL,JTOP,NRMAX)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
            DO I = 1,JTOP
C
               CALL SUMUPINT(CSUM,WGTG,AGBG(1,1,I),AMEG,WGTF,AFBF(1,1,I)
     &                       ,AMEF,NSOL)
               CIAB(K1,K2,I) = CSUM
C
            END DO
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
         END DO
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      END
