C*==chiintghff3.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTGHFF3(CIHF,AG,BG,AF,BF,WGTHF,AMEHF,WGTR,NSOL,
     &                       JTOP,NRMAX)
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
C*--CHIINTGHFF316
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NRMAX,NSOL
      REAL*8 WGTHF
      COMPLEX*16 AF(NRMAX,2,2),AG(NRMAX,2,2),BF(NRMAX,2,2),BG(NRMAX,2,2)
     &           ,CIHF(2,2,NRMAX),WGTR(NRMAX)
      REAL*8 AMEHF(2,2)
C
C Local variables
C
      COMPLEX*16 CINTHF(2,2,NRMAX),CSUM
      INTEGER I,K1,K2
C
C*** End of declarations rewritten by SPAG
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
C
            CALL CHIINTGHFF2(AG(1,1,K1),AF(1,1,K1),BG(1,1,K2),BF(1,1,K2)
     &                       ,CINTHF,NSOL,NSOL,JTOP,WGTR,NRMAX)
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
            DO I = 1,JTOP
               CALL SUMUPINT(CSUM,WGTHF,CINTHF(1,1,I),AMEHF,0.0D0,
     &                       CINTHF(1,1,I),AMEHF,NSOL)
               CIHF(K1,K2,I) = CSUM
            END DO
C
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
         END DO
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      END
