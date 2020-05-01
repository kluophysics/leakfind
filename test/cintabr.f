C*==cintabwr.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINTABWR(AG,BG,AGBG,AF,BF,AFBF,WR,NKA,NKB,JTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  integration routine for complex integrands from x = 1 ... JTOP  *
C   *                                                                  *
C   *            FX = AG*BG*RPW   and   FX = AF*BF*RPW                 *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CINTABWR12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,2),AFBF(2,2),AG(NRMAX,2),AGBG(2,2),BF(NRMAX,2)
     &           ,BG(NRMAX,2)
      REAL*8 WR(NRMAX)
C
C Local variables
C
      COMPLEX*16 BFW,BGW
      INTEGER I,KA,KB
C
C*** End of declarations rewritten by SPAG
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB) = C0
            AFBF(KA,KB) = C0
         END DO
      END DO
C
      DO I = 1,JTOP
         DO KB = 1,NKB
            BGW = BG(I,KB)*WR(I)
            BFW = BF(I,KB)*WR(I)
            DO KA = 1,NKA
               AGBG(KA,KB) = AGBG(KA,KB) + AG(I,KA)*BGW
               AFBF(KA,KB) = AFBF(KA,KB) + AF(I,KA)*BFW
            END DO
         END DO
      END DO
C
      END
C*==cintabr.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINTABR(AG,BG,AGBG,AF,BF,AFBF,RPW,NKA,NKB,JTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  SIMPSON - INTERGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP *
C   *  AND EQUIDISTANT MESH    I                                       *
C   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
C   *                                                                  *
C   *            FX = AG*AG*RPW   and   FX = AF*AF*RPW                 *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CINTABR72
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,2),AFBF(2,2),AG(NRMAX,2),AGBG(2,2),BF(NRMAX,2)
     &           ,BG(NRMAX,2)
      REAL*8 RPW(NRMAX)
C
C Local variables
C
      COMPLEX*16 BFW,BGW
      REAL*8 F,SIMP
      INTEGER I,IST,KA,KB
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( MOD(JTOP,2).EQ.0 ) THEN
         IST = 2
      ELSE
         IST = 1
      END IF
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB) = AG(IST,KA)*BG(IST,KB)*RPW(IST)
            AFBF(KA,KB) = AF(IST,KA)*BF(IST,KB)*RPW(IST)
         END DO
      END DO
C
      SIMP = -1.0D0
      DO I = IST + 1,JTOP - 1
         SIMP = -SIMP
         F = (3.0D0+SIMP)*RPW(I)
         DO KB = 1,NKB
            BGW = BG(I,KB)*F
            BFW = BF(I,KB)*F
            DO KA = 1,NKA
               AGBG(KA,KB) = AGBG(KA,KB) + AG(I,KA)*BGW
               AFBF(KA,KB) = AFBF(KA,KB) + AF(I,KA)*BFW
            END DO
         END DO
      END DO
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB) = (AGBG(KA,KB)+AG(JTOP,KA)*BG(JTOP,KB)*RPW(JTOP)
     &                    )/3.0D0
            AFBF(KA,KB) = (AFBF(KA,KB)+AF(JTOP,KA)*BF(JTOP,KB)*RPW(JTOP)
     &                    )/3.0D0
         END DO
      END DO
C
      END
