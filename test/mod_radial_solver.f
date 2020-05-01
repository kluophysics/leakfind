C*==mod_radial_solver.f    processed by SPAG 6.70Rc at 17:54 on 22 Jun 2015
C
      MODULE MOD_RADIAL_SOLVER
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tvariables used for   RADIAL SOLVERS        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
      REAL*8 CGD(2),CGMD(2),CGO,CSQR,KAP(2)
      REAL*8 LIM
      REAL*8 CGOZ,CGZ(2),KPX(2),LMK(2)
C
      COMPLEX*16 EBS
C
      INTEGER JLAG1,NSOLBS
      INTEGER K_SOC
      INTEGER LBS

      LOGICAL KNS
C
      SAVE EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,JLAG1
      SAVE LIM
      SAVE KNS
      SAVE CGOZ,CGZ,KPX,LMK
      SAVE K_SOC
      SAVE LBS
C                                        
      END

