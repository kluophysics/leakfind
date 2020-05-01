C*==cintabr1.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINTABR1(AG,BG,AGBG,AF,BF,AFBF,WINTR,NKA,NKB,JBOT,JTOP,
     &                    NRMAX,NCPLWFMAX)
C   ********************************************************************
C   *                                                                  *
C   *  perform a radial integration with the r-dependent weights       *
C   *  set before                                                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CINTABR111
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JBOT,JTOP,NCPLWFMAX,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,NCPLWFMAX),AFBF(NCPLWFMAX,NCPLWFMAX),
     &           AG(NRMAX,NCPLWFMAX),AGBG(NCPLWFMAX,NCPLWFMAX),
     &           BF(NRMAX,NCPLWFMAX),BG(NRMAX,NCPLWFMAX)
      REAL*8 WINTR(NRMAX)
C
C Local variables
C
      COMPLEX*16 BFW(NRMAX),BGW(NRMAX),SUMF,SUMG
      INTEGER I,KA,KB
C
C*** End of declarations rewritten by SPAG
C
      DO KB = 1,NKB
C
         DO I = JBOT,JTOP
            BGW(I) = BG(I,KB)*WINTR(I)
            BFW(I) = BF(I,KB)*WINTR(I)
         END DO
C
         DO KA = 1,NKA
            SUMG = 0D0
            SUMF = 0D0
            DO I = JBOT,JTOP
               SUMG = SUMG + AG(I,KA)*BGW(I)
               SUMF = SUMF + AF(I,KA)*BFW(I)
            END DO
            AGBG(KA,KB) = SUMG
            AFBF(KA,KB) = SUMF
         END DO
C
      END DO
C
      END
