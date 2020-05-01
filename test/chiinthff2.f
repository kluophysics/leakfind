C*==chiinthff2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTHFF2(IM,AG,AF,BG,BF,RMEHF,NKA,NKB,IRTOP,FX,WGTR)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   *  TODO: use cradint_hff as auxilary routine                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CHIINTHFF215
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,NKA,NKB
      COMPLEX*16 AF(NRMAX,NKA),AG(NRMAX,NKA),BF(NRMAX,NKB),BG(NRMAX,NKB)
     &           ,FX(IRTOP),RMEHF(NCPLWFMAX,NCPLWFMAX),WGTR(NRMAX)
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
               FX(IR) = (AG(IR,KA)*BF(IR,KB)+AF(IR,KA)*BG(IR,KB))
     &                  *WGTR(IR)
            END DO
C
            CALL CRADINT_HFF(IM,IRTOP,FX,RMEHF(KA,KB))
C
         END DO
      END DO
C
      END
