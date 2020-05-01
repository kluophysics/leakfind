C*==excvwn.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCVWN(RHO4PI,V,EXC,WEXC,IRTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *     sdf parametrization by Vosko, Wilk and Nusair  (VWN)         *
C   *                                                                  *
C   *     RHO4PI(IR,1) = 4 pi (rho(up)+rho(dn)                         *
C   *     RHO4PI(IR,2) = 4 pi (rho(up)-rho(dn)                         *
C   *                                                                  *
C   *     calling  UXCOR  coded by M. Mannien                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--EXCVWN15
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,NRMAX
      REAL*8 EXC(NRMAX),RHO4PI(NRMAX,2),V(NRMAX,2),WEXC(NRMAX,2)
C
C Local variables
C
      INTEGER IR
      REAL*8 RS,RT,U,VXCDN,VXCUP
C
C*** End of declarations rewritten by SPAG
C
      DO IR = 1,IRTOP
         IF ( RHO4PI(IR,1).GT.1D-20 ) THEN
C
            RS = (3D0/RHO4PI(IR,1))**.333333333333D0
            RT = RHO4PI(IR,2)/RHO4PI(IR,1)
            RT = MIN(1D0,MAX(-1D0,RT))
C
            CALL UXCOR(RS,RT,VXCUP,VXCDN,EXC(IR))
C
            V(IR,1) = V(IR,1) + VXCUP
            V(IR,2) = V(IR,2) + VXCDN
C
            U = WEXC(IR,1)
            WEXC(IR,1) = RHO4PI(IR,1)*(U+EXC(IR))
            WEXC(IR,2) = RHO4PI(IR,1)*U - 3D0*(EXC(IR)-VXCUP)
     &                   *(RHO4PI(IR,1)+RHO4PI(IR,2))
     &                   *5D-1 - 3D0*(EXC(IR)-VXCDN)
     &                   *(RHO4PI(IR,1)-RHO4PI(IR,2))*5D-1
         ELSE
            WEXC(IR,1:2) = 0D0
         END IF
      END DO
      END
