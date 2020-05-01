C*==vcepbm.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE VCEPBM(RHO4PI,V,IRTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *     LDA parametrization for the                                  *
C   *     electron-positron correlation potential                      *
C   *     by Boronski and Nieminen  PRB 34 p3820 (1986) (Appendix)     *
C   *                                                                  *
C   *     electronic densities                                         *
C   *     RHO4PI(IR,1) = 4 pi (rho(up)+rho(dn)                         *
C   *     RHO4PI(IR,2) = 4 pi (rho(up)-rho(dn)                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      IMPLICIT NONE
C*--VCEPBM17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,NRMAX
      REAL*8 RHO4PI(NRMAX,2),V(NRMAX,2)
C
C Local variables
C
      INTEGER IR
      REAL*8 LNRS,RHO,RS,RS25,THIRD,VC
C
C*** End of declarations rewritten by SPAG
C
      THIRD = 1D0/3D0
C
      DO IR = 1,IRTOP
         IF ( RHO4PI(IR,1).GT.1D-20 ) THEN
C
            RS = (3D0/RHO4PI(IR,1))**THIRD
C
            IF ( RS.LE.0.302D0 ) THEN
C
               LNRS = LOG(RS)
C
               VC = -1.56D0/SQRT(RS) + (0.051D0*LNRS-0.081D0)
     &              *LNRS + 1.14D0
C
            ELSE IF ( RS.LE.0.56D0 ) THEN
C
               VC = -0.92305D0 - 0.05459D0/(RS*RS)
C
            ELSE IF ( RS.LE.8D0 ) THEN
C
               RS25 = RS + 2.5D0
C
               VC = -0.6298D0 - 13.15111D0/(RS25*RS25) + 2.8655D0/RS25
C
            ELSE
C
               RHO = RHO4PI(IR,1)/CONST_4PI
C
               VC = -0.524D0 - 179856.2768D0*RHO*RHO + 186.4207D0*RHO
C
            END IF
C
            V(IR,1) = V(IR,1) + VC
            V(IR,2) = V(IR,1)
C
         END IF
      END DO
      END
