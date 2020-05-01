C*==chiintabr2in.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTABR2IN(AG,BG,AGBG,AF,BF,AFBF,WGTR,NKA,NKB,JTOP,
     &                        NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  CUMULATIVE INTEGRAL FROM 1 TO JTOP FOR PRODUCT OF LARGE/SMALL   *
C   *    COMPONENTS OF REG/IRR. SOLUTIONS OF DIRAC EQUATION            *
C   *                                                                  *
C   *  SIMPSON - INTEGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP  *
C   *  AND EQUIDISTANT MESH    I                                       *
C   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
C   *                                                                  *
C   *  SIMPSON STEPS FOR ODD INDICES, BACK STEP WITH TRAPEZOID RULE    *
C   *    FOR EVEN INDICES                                              *
C   *                                                                  *
C   *          FX = AG*BG*WGTR   and   FX = AF*BF*WGTR                 *
C   *                                                                  *
C   *  HF/MD: 12/96                                                    *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CHIINTABR2IN22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 WSIMP,WTRPZ
      PARAMETER (WSIMP=1D0/3D0,WTRPZ=1D0/2D0)
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
      COMPLEX*16 INTFF(2,2,NRMAX),INTGG(2,2,NRMAX)
C
C*** End of declarations rewritten by SPAG
C
      IF ( MOD(JTOP,2).EQ.0 ) STOP '<CINTABR>  JTOP is even !!!'
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB,JTOP) = C0
            AFBF(KA,KB,JTOP) = C0
            DO I = 1,JTOP
               INTGG(KA,KB,I) = AG(I,KA)*BG(I,KB)*WGTR(I)
               INTFF(KA,KB,I) = AF(I,KA)*BF(I,KB)*WGTR(I)
            END DO
         END DO
      END DO
C
      DO I = JTOP,3, - 2
         DO KB = 1,NKB
            DO KA = 1,NKA
               AGBG(KA,KB,I-2) = AGBG(KA,KB,I)
     &                           + WSIMP*(INTGG(KA,KB,I-2)+
     &                           4D0*INTGG(KA,KB,I-1)+INTGG(KA,KB,I))
               AGBG(KA,KB,I-1) = AGBG(KA,KB,I)
     &                           + WTRPZ*(INTGG(KA,KB,I-1)+
     &                           INTGG(KA,KB,I))
               AFBF(KA,KB,I-2) = AFBF(KA,KB,I)
     &                           + WSIMP*(INTFF(KA,KB,I-2)+
     &                           4D0*INTFF(KA,KB,I-1)+INTFF(KA,KB,I))
               AFBF(KA,KB,I-1) = AFBF(KA,KB,I)
     &                           + WTRPZ*(INTFF(KA,KB,I-1)+
     &                           INTFF(KA,KB,I))
            END DO
         END DO
      END DO
      END
