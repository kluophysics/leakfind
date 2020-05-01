C*==posani_gfactor.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANI_GFACTOR(RHO_EL,GFACTOR,IRTOP)
C   ********************************************************************
C   *                                                                  *
C   *     electronic density                                           *
C   *     RHO_EL(IR) = rho(up)+rho(dn)                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      IMPLICIT NONE
C*--POSANI_GFACTOR12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP
      REAL*8 GFACTOR(IRTOP),RHO_EL(IRTOP)
C
C Local variables
C
      INTEGER IR
      REAL*8 RS,RSPW2,RSPW3,RSPW3OV2,RSPW5OV2,SQRT_RS,THIRD
C
C*** End of declarations rewritten by SPAG
C
      THIRD = 1D0/3D0
C
      DO IR = 1,IRTOP
         IF ( RHO_EL(IR).GT.1D-20 ) THEN
C
            RS = (3D0/(CONST_4PI*RHO_EL(IR)))**THIRD
C
            SQRT_RS = SQRT(RS)
            RSPW3OV2 = RS*SQRT_RS
            RSPW5OV2 = RS*RSPW3OV2
            RSPW2 = RS*RS
            RSPW3 = RSPW2*RS
C
            GFACTOR(IR) = 1D0 + 1.23D0*RS + 0.8295D0*RSPW3OV2 - 
     &                    1.26*RSPW2 + 0.3286D0*RSPW5OV2 + RSPW3/6D0
C
         ELSE
C
            GFACTOR(IR) = 0D0
C
         END IF
      END DO
C
      END
