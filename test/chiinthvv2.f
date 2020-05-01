C*==chiinthvv2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTHVV2(IM,AG,BG,RMEVV,NKA,NKB,IRTOP,FX,WGTR)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CHIINTHVV212
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,NKA,NKB
      COMPLEX*16 AG(NRMAX,NKA),BG(NRMAX,NKA),FX(IRTOP),RMEVV(2,2),
     &           WGTR(NRMAX)
C
C Local variables
C
      INTEGER IR,KA,KB
C
C*** End of declarations rewritten by SPAG
C
      DO KB = 1,NKB
         DO KA = 1,NKA
C
            DO IR = 1,IRTOP
               FX(IR) = AG(IR,KA)*BG(IR,KB)*WGTR(IR)
            END DO
C
            CALL CRADINT_HFF(IM,IRTOP,FX,RMEVV(KA,KB))
C
         END DO
      END DO
C
      END
