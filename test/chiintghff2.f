C*==chiintghff2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTGHFF2(AG,AF,BG,BF,CINTHF,NKA,NKB,JTOP,WGTR,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--CHIINTGHFF210
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,NKA),AG(NRMAX,NKA),BF(NRMAX,NKB),BG(NRMAX,NKB)
     &           ,CINTHF(2,2,NRMAX),WGTR(NRMAX)
C
C Local variables
C
      INTEGER I,KA,KB
C
C*** End of declarations rewritten by SPAG
C
      DO KB = 1,NKB
         DO KA = 1,NKA
C
            DO I = 1,JTOP
               CINTHF(KA,KB,I) = (AG(I,KA)*BF(I,KB)+AF(I,KA)*BG(I,KB))
     &                           *WGTR(I)
            END DO
C
         END DO
      END DO
C
      END
