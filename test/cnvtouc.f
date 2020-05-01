C*==cnvtouc.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CNVTOUC(STR,I1,I2)
      IMPLICIT NONE
C*--CNVTOUC4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I1,I2
      CHARACTER*(*) STR
C
C Local variables
C
      INTEGER I,IA,INC,IZ,J
C
C*** End of declarations rewritten by SPAG
C
C   ********************************************************************
C   *                                                                  *
C   *   convert string STR(L1:L2)  to upper case                       *
C   *                                                                  *
C   ********************************************************************
C
      IA = ICHAR('a')
      IZ = ICHAR('z')
      INC = ICHAR('A') - IA
C
      DO I = I1,I2
         J = ICHAR(STR(I:I))
         IF ( J.GE.IA .AND. J.LE.IZ ) STR(I:I) = CHAR(J+INC)
      END DO
      END
