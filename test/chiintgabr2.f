C*==chiintgabr2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTGABR2(AG,BG,AGBG,AF,BF,AFBF,WGTR,NKA,NKB,JTOP,
     &                       NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  SET UP INTEGRAND    FROM 1 TO JTOP FOR PRODUCT OF LARGE/SMALL   *
C   *    COMPONENTS OF REG/IRR. SOLUTIONS OF DIRAC EQUATION            *
C   *                                                                  *
C   *          N: ODD (FOR LATER SIMPSON INTEGRATION)                  *
C   *          FX = AG*BG*WGTR   and   FX = AF*BF*WGTR                 *
C   *                                                                  *
C   *  HF:  2/99                                                       *
C   ********************************************************************
      IMPLICIT NONE
C*--CHIINTGABR215
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,NKA),AFBF(2,2,NRMAX),AG(NRMAX,NKA),
     &           AGBG(2,2,NRMAX),BF(NRMAX,NKB),BG(NRMAX,NKB),WGTR(NRMAX)
C
C Local variables
C
      INTEGER I,KA,KB
C
C*** End of declarations rewritten by SPAG
C
      IF ( MOD(JTOP,2).EQ.0 ) STOP '<CINTABR>  JTOP is even !!!'
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            DO I = 1,JTOP
               AGBG(KA,KB,I) = AG(I,KA)*BG(I,KB)*WGTR(I)
               AFBF(KA,KB,I) = AF(I,KA)*BF(I,KB)*WGTR(I)
            END DO
         END DO
      END DO
      END
