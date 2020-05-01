C*==dmft_kkt_imre.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_KKT_IMRE(N,X,YIM,WORK,YRE,DX)
      IMPLICIT NONE
C*--DMFT_KKT_IMRE4
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DX
      INTEGER N
      REAL*8 WORK(N),X(N),YIM(N),YRE(N)
C
C Local variables
C
      REAL*8 DELTA2,FFF,TMP
      INTEGER I,IE
C
C*** End of declarations rewritten by SPAG
C
C Dummy arguments
C Local variables
C
      TMP = 1D0/(4D0*DATAN(1D0))
      DELTA2 = DX*DX
      DO IE = 1,N
         DO I = 1,N
            FFF = X(IE) - X(I)
            WORK(I) = FFF/(FFF*FFF+DELTA2)*YIM(I)
         END DO
         CALL DMFT_INTEGRAL(YRE(IE),DX,WORK,N,1,N)
         YRE(IE) = -YRE(IE)*TMP
      END DO
      END
C*==dmft_kkt_extrapol.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_KKT_EXTRAPOL(N,X,YIM,NC,Z,YRES)
      IMPLICIT NONE
C*--DMFT_KKT_EXTRAPOL51
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,NC
      REAL*8 X(N),YIM(N)
      COMPLEX*16 YRES(NC),Z(NC)
C
C Local variables
C
      COMPLEX*16 C0
      REAL*8 DELTA,SHIFT,TMP
      INTEGER IC,IE
C
C*** End of declarations rewritten by SPAG
C
C
C     Dummy arguments
C     Local variables
C
      C0 = DCMPLX(0D0,0D0)
      TMP = 0.5D0/(4*DATAN(1D0))
      DO IC = 1,NC
         YRES(IC) = C0
         DO IE = 1,N - 1
            DELTA = 0D0
            IF ( CDABS(X(IE)-Z(IC)).GE.1D-20 .AND. CDABS(X(IE+1)-Z(IC))
     &           .GE.1D-20 ) DELTA = 1D-20
            YRES(IC) = YRES(IC) + (YIM(IE)+YIM(IE+1))
     &                 *LOG((X(IE+1)-Z(IC))/(X(IE)-Z(IC)+DELTA))
         END DO
         YRES(IC) = YRES(IC)*TMP
      END DO
C
      SHIFT = DREAL(YRES(NC))
      DO IC = 1,NC
         YRES(IC) = DCMPLX(DREAL(YRES(IC))-SHIFT,DIMAG(YRES(IC)))
      END DO
      YRES(NC) = C0
      END
C*==dmft_hilbertt.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_HILBERTT(N,X,YINP,WORK,YOUT,DX)
      IMPLICIT NONE
C*--DMFT_HILBERTT109
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DX
      INTEGER N
      COMPLEX*16 WORK(N),YINP(N),YOUT(N)
      REAL*8 X(N)
C
C Local variables
C
      COMPLEX*16 C1,CTMP
      REAL*8 DELTA2,FFF
      INTEGER I,IE
C
C*** End of declarations rewritten by SPAG
C
C Dummy arguments
C Local variables
      C1 = DCMPLX(0D0,1D0)
      CTMP = C1/(4D0*DATAN(1D0))
      DO IE = 1,N
         DO I = 1,N
            DELTA2 = 0D0
            IF ( IE.EQ.I ) DELTA2 = 1D-10
            FFF = X(IE) - X(I)
            WORK(I) = FFF/(FFF*FFF+DELTA2)*YINP(I)
         END DO
         CALL DMFT_CINTEGRAL(YOUT(IE),DX,WORK,N,1,N)
      END DO
      YOUT(1:N) = YOUT(1:N)*CTMP
      END
C*==dmft_cauchyt.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
      SUBROUTINE DMFT_CAUCHYT(N,X,Y,NC,Z,YRES)
      IMPLICIT NONE
C*--DMFT_CAUCHYT161
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,NC
      REAL*8 X(N)
      COMPLEX*16 Y(N),YRES(NC),Z(NC)
C
C Local variables
C
      COMPLEX*16 C0,C1,CTMP
      INTEGER IC,IE
C
C*** End of declarations rewritten by SPAG
C
C
C Dummy arguments
C Local variables
C$$$      COMPLEX*16 CWORK(:)
C$$$      REAL*8 DX
C
C$$$      ALLOCATABLE CWORK
C
      C0 = DCMPLX(0D0,0D0)
      C1 = DCMPLX(0D0,1D0)
      CTMP = 0.5D0/3.141592653589793238462643D0/C1
      YRES(1:NC) = C0
C$$$      DX = (X(N)-X(1))/(N-1D0)
C
C$$$      ALLOCATE (CWORK(N))
C
      DO IC = 1,NC
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C method 1:
         DO IE = 1,N - 1
            IF ( CDABS(X(IE)-Z(IC)).GE.1D-20 .AND. CDABS(X(IE+1)-Z(IC))
     &           .GE.1D-20 ) YRES(IC) = YRES(IC) + (Y(IE)+Y(IE+1))
     &                                  *LOG((X(IE+1)-Z(IC))
     &                                  /(X(IE)-Z(IC)))
         END DO
C method 2:
C$$$         DO IE = 1,N
C$$$            CWORK(IE) = 2d0*Y(IE)/(X(IE)-Z(IC))
C$$$         END DO
C$$$         CALL DMFT_CINTEGRAL(YRES(IC),DX,CWORK,N,1,N)
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END DO
C
      YRES(1:NC) = YRES(1:NC)*CTMP
C
C$$$      DEALLOCATE (CWORK)
C
C
C
C
      END
