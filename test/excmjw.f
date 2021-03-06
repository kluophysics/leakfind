C*==excmjw.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCMJW(RHO4PI,V,EXC,WEXC,IRTOP,NRMAX,IXC)
C   ********************************************************************
C   *                                                                  *
C   *     sdf parametrization by von Moruzzi, Janak, Williams  (MJW)   *
C   *                                                                  *
C   *     Calculated Electronic Properties of Metals                   *
C   *     Pergamon, New York, 1978                                     *
C   *                                                                  *
C   *     RHO4PI(IR,1) = 4 pi (rho(up)+rho(dn)                         *
C   *     RHO4PI(IR,2) = 4 pi (rho(up)-rho(dn)                         *
C   *                                                                  *
C   *     IXC = 0  exchange only                                       *
C   *     IXC = 1  exchange + correlation                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--EXCMJW18
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,IXC,NRMAX
      REAL*8 EXC(NRMAX),RHO4PI(NRMAX,2),V(NRMAX,2),WEXC(NRMAX,2)
C
C Local variables
C
      REAL*8 A,C,C1,C2,CF,CFL,CNU,CP,CPL,ECF,ECPM,EXPM,FCF,FCPM,FF,FXPM,
     &       GAMMA,RF,RP,RS,RT,S,T,TCF,U,X,XD,XU,Y
      REAL*8 F,G
      INTEGER IR
C
C*** End of declarations rewritten by SPAG
C
      DATA T/.333333333333D0/,A/7.93700526D-1/,GAMMA/5.129762803D0/,CP,
     &     RP,CF,RF/4.5D-2,21D0,2.25D-2,52.916684096D0/
C    &    ,cp,rp,cf,rf/5.04d-2, 30d0,  2.54d-2,     75d0/
C
      G(X) = LOG(1D0+1D0/X)
      F(C,X) = ((C*X-1D0)*X+5D-1)*X + C - T
C
C          a=2d0**(-t)
C          gamma=(4d0/3d0)*a/(1d0-a)
C      (  for mjw parameter  cf=cp/2 and rf=2**(4/3)*a/(1-a) holds )
C
C-----------------------------------------------------------------------
C                            exchange only
C-----------------------------------------------------------------------
      IF ( IXC.EQ.0 ) THEN
C
         DO IR = 1,IRTOP
            IF ( RHO4PI(IR,1).GT.1D-20 ) THEN
C
               RS = (3D0/RHO4PI(IR,1))**T
               RT = RHO4PI(IR,2)/RHO4PI(IR,1)
               RT = MIN(1D0,MAX(-1D0,RT))
C
               S = 5D-1*(1D0+RT)
               FF = (S**(4D0/3D0)+(1D0-S)**(4D0/3D0)-A)/(1D0-A)
               X = RS/RP
               Y = RS/RF
               FXPM = -1.221774115422D0/RS
               EXPM = FXPM*.75D0
               EXC(IR) = EXPM + FXPM*FF/GAMMA
               C1 = FXPM
               XU = C1*(1D0+RT)**T
               XD = C1*(1D0-RT)**T
               V(IR,1) = V(IR,1) + XU
               V(IR,2) = V(IR,2) + XD
C
               U = WEXC(IR,1)
               WEXC(IR,1) = RHO4PI(IR,1)*(U+EXC(IR))
               WEXC(IR,2) = RHO4PI(IR,1)*U - 3D0*(EXC(IR)-XU)
     &                      *(RHO4PI(IR,1)+RHO4PI(IR,2))
     &                      *5D-1 - 3D0*(EXC(IR)-XD)
     &                      *(RHO4PI(IR,1)-RHO4PI(IR,2))*5D-1
C
            ELSE
               WEXC(IR,1:2) = 0D0
            END IF
         END DO
C
C-----------------------------------------------------------------------
C                      exchange + correlation
C-----------------------------------------------------------------------
      ELSE
C
         DO IR = 1,IRTOP
            IF ( RHO4PI(IR,1).GT.1D-20 ) THEN
C
               RS = (3D0/RHO4PI(IR,1))**T
               RT = RHO4PI(IR,2)/RHO4PI(IR,1)
               RT = MIN(1D0,MAX(-1D0,RT))
C
               S = 5D-1*(1D0+RT)
               FF = (S**(4D0/3D0)+(1D0-S)**(4D0/3D0)-A)/(1D0-A)
               X = RS/RP
               Y = RS/RF
               FXPM = -1.221774115422D0/RS
               EXPM = FXPM*.75D0
               CPL = G(X)
               CFL = G(Y)
               FCPM = -CP*CPL
               FCF = -CF*CFL
               ECPM = -CP*F(CPL,X)
               ECF = -CF*F(CFL,Y)
               CNU = GAMMA*(ECF-ECPM)
               EXC(IR) = EXPM + ECPM + (FXPM+CNU)*FF/GAMMA
               C1 = FXPM + CNU
               C2 = FCPM - CNU
               TCF = (FCF-FCPM-4D0*(ECF-ECPM)/3D0)*FF
               XU = C1*(1D0+RT)**T + C2 + TCF
               XD = C1*(1D0-RT)**T + C2 + TCF
               V(IR,1) = V(IR,1) + XU
               V(IR,2) = V(IR,2) + XD
C
               U = WEXC(IR,1)
               WEXC(IR,1) = RHO4PI(IR,1)*(U+EXC(IR))
               WEXC(IR,2) = RHO4PI(IR,1)*U - 3D0*(EXC(IR)-XU)
     &                      *(RHO4PI(IR,1)+RHO4PI(IR,2))
     &                      *5D-1 - 3D0*(EXC(IR)-XD)
     &                      *(RHO4PI(IR,1)-RHO4PI(IR,2))*5D-1
C
            ELSE
               WEXC(IR,1:2) = 0D0
            END IF
         END DO
C
      END IF
C
      END
