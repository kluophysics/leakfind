C*==dmft_interpolatec.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_INTERPOLATEC(N0,X0,Y0,N,X,Y)
      IMPLICIT NONE
C*--DMFT_INTERPOLATEC4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,N0
      REAL*8 X(N),X0(N0)
      COMPLEX*16 Y(N),Y0(N0)
C
C Local variables
C
      REAL*8 AUX(:),AUX0(:),WORK(:),Y20(:)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE AUX0,AUX,WORK,Y20
C
      ALLOCATE (AUX0(N0),AUX(N),WORK(N0),Y20(N0))
C
      AUX0(1:N0) = DREAL(Y0(1:N0))
      CALL DMFT_INTERPOLATE(N0,X0,AUX0,Y20,WORK,N,X,AUX)
      DO I = 1,N
         Y(I) = DCMPLX(AUX(I),0D0)
      END DO
      AUX0(1:N0) = DIMAG(Y0(1:N0))
      CALL DMFT_INTERPOLATE(N0,X0,AUX0,Y20,WORK,N,X,AUX)
      DO I = 1,N
         Y(I) = Y(I) + DCMPLX(0D0,AUX(I))
      END DO
      DEALLOCATE (AUX0,AUX,WORK,Y20)
C
      END
C*==dmft_interpolate.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_INTERPOLATE(N0,X0,Y0,Y20,WORK,N,X,Y)
C PARABOLA SPLINE INTERPOLATION
      IMPLICIT NONE
C*--DMFT_INTERPOLATE56
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,N0
      DOUBLE PRECISION WORK(N0),X(N),X0(N0),Y(N),Y0(N0),Y20(N0)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C     driver for routine splint, which calls spline
      CALL DMFT_SPLINE(X0,Y0,N0,0D0,0D0,Y20,WORK)
      DO I = 1,N
         IF ( X(I).LT.X0(1) .OR. X(I).GT.X0(N0) ) THEN
            Y(I) = 0D0
         ELSE
            CALL DMFT_SPLINT(X0,Y0,Y20,N0,X(I),Y(I))
         END IF
      END DO
      END
C*==dmft_linnterp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
      SUBROUTINE DMFT_LINNTERP(N0,X0,Y0,X,Y)
C LINEAR INTERPOLATION
      IMPLICIT NONE
C*--DMFT_LINNTERP100
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N0
      DOUBLE PRECISION X,Y
      DOUBLE PRECISION X0(N0),Y0(N0)
C
C Local variables
C
      INTEGER I,IMIN
      DOUBLE PRECISION XMAX,XMIN,YMAX,YMIN
C
C*** End of declarations rewritten by SPAG
C
      IF ( X.LT.X0(1) .OR. X.GT.X0(N0) ) THEN
         Y = 0D0
      ELSE
         DO I = 1,N0 - 1
            IF ( X0(I).LE.X ) IMIN = I
         END DO
         XMAX = X0(IMIN+1)
         YMAX = Y0(IMIN+1)
         XMIN = X0(IMIN)
         YMIN = Y0(IMIN)
         Y = (YMAX*(X-XMIN)+YMIN*(XMAX-X))/(XMAX-XMIN)
      END IF
C
      END
C*==dmft_splint.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
      SUBROUTINE DMFT_SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
C*--DMFT_SPLINT146
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      DOUBLE PRECISION X,Y
      DOUBLE PRECISION XA(N),Y2A(N),YA(N)
C
C Local variables
C
      DOUBLE PRECISION A,B,H
      INTEGER K,KHI,KLO
C
C*** End of declarations rewritten by SPAG
C
      KLO = 1
      KHI = N
 100  CONTINUE
      IF ( KHI-KLO.GT.1 ) THEN
         K = (KHI+KLO)/2
         IF ( XA(K).GT.X ) THEN
            KHI = K
         ELSE
            KLO = K
         END IF
         GOTO 100
      END IF
      H = XA(KHI) - XA(KLO)
      IF ( ABS(H).LE.1.D-12 ) STOP 'bad xa input in splint'
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))
     &    *(H**2)/6D0
      END
C*==dmft_spline.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_SPLINE(X,Y,N,YP1,YPN,Y2,U)
      IMPLICIT NONE
C*--DMFT_SPLINE198
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      DOUBLE PRECISION YP1,YPN
      DOUBLE PRECISION U(N),X(N),Y(N),Y2(N)
C
C Local variables
C
      INTEGER I,K
      DOUBLE PRECISION P,QN,SIG,UN
C
C*** End of declarations rewritten by SPAG
C
      IF ( YP1.GT..99D30 ) THEN
         Y2(1) = 0D0
         U(1) = 0D0
      ELSE
         Y2(1) = -0.5D0
         U(1) = (3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      END IF
      DO I = 2,N - 1
         SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
         P = SIG*Y2(I-1) + 2.D0
         Y2(I) = (SIG-1.D0)/P
         U(I) = (6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X
     &          (I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO
      IF ( YPN.GT..99D30 ) THEN
         QN = 0D0
         UN = 0D0
      ELSE
         QN = 0.5D0
         UN = (3D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      END IF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1D0)
      DO K = N - 1,1, - 1
         Y2(K) = Y2(K)*Y2(K+1) + U(K)
      END DO
      END
