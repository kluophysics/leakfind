C*==dmft_integral.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_INTEGRAL(RINT,DX,F,N,II1,II2)
      IMPLICIT NONE
C*--DMFT_INTEGRAL4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DX,RINT
      INTEGER II1,II2,N
      REAL*8 F(N)
C
C Local variables
C
      REAL*8 H
      INTEGER I,I1,I2,N1,NSTEP
C
C*** End of declarations rewritten by SPAG
C
C Take the integral using Gregorie method (5-th precision order).
C If number of points is less than NSTEP, then routine uses boxes.
C
C
C Dummy arguments
C
C
C
C Local variables
C
C
      NSTEP = 1
C
      H = DX
      IF ( II2.LT.II1 ) THEN
         I1 = II2
         I2 = II1
      ELSE
         I1 = II1
         I2 = II2
      END IF
C
      IF ( I2-I1+1.GT.NSTEP ) THEN
         RINT = 0.375D0*(F(I1)+F(I2)) + (7.D0/6.D0)*(F(I1+1)+F(I2-1))
     &          + (23.D0/24.D0)*(F(I1+2)+F(I2-2))
         N1 = I2 - 3
         DO I = I1 + 3,N1
            RINT = RINT + F(I)
         END DO
      ELSE IF ( I2-I1.LE.NSTEP .AND. I2.GT.I1 ) THEN
         RINT = 0D0
         DO I = I1 + 1,I2 - 1
            RINT = RINT + F(I)
         END DO
         RINT = RINT + (F(I1)+F(I2))/2D0
      ELSE
         RINT = 0D0
      END IF
      RINT = RINT*H
      END
C*==dmft_cintegral.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C$$$      SUBROUTINE CINTEGRALXX(CINT,DX,F,N,II1,II2)
C$$$      IMPLICIT NONE
C$$$C Dummy arguments
C$$$      REAL*8 DX
C$$$      COMPLEX*16 CINT,F(N)
C$$$      INTEGER II1,II2,N
C$$$C Local variables
C$$$      COMPLEX*16 C0
C$$$      REAL*8 H
C$$$      INTEGER I,I1,I2,N1
C$$$C
C$$$      C0=dcmplx(0d0,0d0)
C$$$      CINT=C0
C$$$      H = DX
C$$$      IF ( II2.LT.II1 ) THEN
C$$$         I1 = II2
C$$$         I2 = II1
C$$$      ELSE
C$$$         I1 = II1
C$$$         I2 = II2
C$$$      END IF
C$$$      DO I=I1,I2-1
C$$$         CINT=CINT+F(I)
C$$$      ENDDO
C$$$      CINT=CINT*DX
C$$$      END
C
C
      SUBROUTINE DMFT_CINTEGRAL(CINT,DX,F,N,II1,II2)
      IMPLICIT NONE
C*--DMFT_CINTEGRAL106
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CINT
      REAL*8 DX
      INTEGER II1,II2,N
      COMPLEX*16 F(N)
C
C Local variables
C
      COMPLEX*16 C0
      REAL*8 H
      INTEGER I,I1,I2,N1,NSTEP
C
C*** End of declarations rewritten by SPAG
C
C Takes the integral using Gregorie method (5-th precision order).
C If number of points is less than NSTEP, then routine uses boxes.
C
C Dummy arguments
C
C
C
C
C Local variables
C
C
C
      C0 = DCMPLX(0D0,0D0)
      NSTEP = 5
C
      H = DX
      IF ( II2.LT.II1 ) THEN
         I1 = II2
         I2 = II1
      ELSE
         I1 = II1
         I2 = II2
      END IF
C
      IF ( I2-I1+1.GT.NSTEP ) THEN
         CINT = 0.375D0*(F(I1)+F(I2)) + (7.D0/6.D0)*(F(I1+1)+F(I2-1))
     &          + (23.D0/24.D0)*(F(I1+2)+F(I2-2))
         N1 = I2 - 3
         DO I = I1 + 3,N1
            CINT = CINT + F(I)
         END DO
      ELSE IF ( I2-I1.LE.NSTEP ) THEN
         CINT = C0
         DO I = I1 + 1,I2 - 1
            CINT = CINT + F(I)
         END DO
         CINT = CINT + (F(I1)+F(I2))/2D0
      ELSE
         CINT = C0
      END IF
      CINT = CINT*H
      END
C*==dmft_cintegral1.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
      SUBROUTINE DMFT_CINTEGRAL1(CINT,DX,F,N,II1,II2)
      IMPLICIT NONE
C*--DMFT_CINTEGRAL1185
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CINT
      REAL*8 DX
      INTEGER II1,II2,N
      COMPLEX*16 F(N)
C
C Local variables
C
      COMPLEX*16 C0
      REAL*8 H
      INTEGER I1,I2
C
C*** End of declarations rewritten by SPAG
C
C
C
      C0 = DCMPLX(0D0,0D0)
C
      IF ( II1.GT.N ) II1 = N
      IF ( II1.LT.N ) II1 = 1
      IF ( II2.GT.N ) II2 = N
      IF ( II2.LT.N ) II2 = 1
C
      IF ( ABS(II2-II1).GE.2 ) THEN
         H = DX
         IF ( II2.LT.II1 ) THEN
            I1 = II2
            I2 = II1
         ELSE
            I1 = II1
            I2 = II2
         END IF
         CALL DMFT_QTS1G(CINT,H,F(I1),I2-I1+1)
      ELSE
         CINT = C0
      END IF
C
C
      END
C*==dmft_qts1g.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_QTS1G(RINT,H,F,N)
      IMPLICIT NONE
C*--DMFT_QTS1G245
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H
      INTEGER N
      COMPLEX*16 RINT
      COMPLEX*16 F(N)
C
C Local variables
C
      COMPLEX*16 F12
      INTEGER I,N2,N3
C
C*** End of declarations rewritten by SPAG
C
      RINT = DCMPLX(0.D0,0.D0)
      N2 = 1
      N3 = N/2
      N3 = 2*N3
      IF ( N3.EQ.N ) THEN
         F12 = (5.D0*F(1)+15.D0*F(2)-5.D0*F(3)+F(4))/16.D0
         RINT = (F(1)+4.D0*F12+F(2))/2.D0
         N2 = 2
      END IF
      RINT = RINT + F(N2) - F(N)
      N2 = N2 + 1
      DO I = N2,N,2
         RINT = RINT + F(I)*4.D0 + F(I+1)*2.D0
      END DO
      RINT = RINT*H/3.D0
      END
C
C
C
C
