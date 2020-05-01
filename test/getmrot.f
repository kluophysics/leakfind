C*==getmrot.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GETMROT(ALF,BET,GAM,MROT)
C   ********************************************************************
C   *                                                                  *
C   *  get the 3x3 real space rotation matrix  MROT(alf,bet,gam)[at]   *
C   *  corresponding to the Euler angles  ALF, BET, GAM (in degree)    *
C   *  according to the ACTIVE-TEMPORARY AXIS convention  (Rose)       *
C   *                                                                  *
C   *  make use of the relation:                                       *
C   *                                                                  *
C   *      MROT(a,b,g)[at] = M(a,z)[af] * M(b,y)[af] * M(g,z)[af]      *
C   *                                                                  *
C   *  using the step by step decomposition in the ACTIVE-FIXED conv.  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--GETMROT19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALF,BET,GAM
      REAL*8 MROT(3,3)
C
C Local variables
C
      REAL*8 CA,CB,CG,MA(3,3),MB(3,3),MBG(3,3),MG(3,3),SA,SB,SG
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------- rotation by ALF around fixed z-axis
      SA = SIN(ALF*PI/180D0)
      CA = COS(ALF*PI/180D0)
C
      MA(1:3,1:3) = 0D0
      MA(1,1) = +CA
      MA(1,2) = -SA
      MA(2,1) = +SA
      MA(2,2) = +CA
      MA(3,3) = 1D0
C
C----------------------------------- rotation by BET around fixed y-axis
      SB = SIN(BET*PI/180D0)
      CB = COS(BET*PI/180D0)
C
      MB(1:3,1:3) = 0D0
      MB(1,1) = +CB
      MB(1,3) = +SB
      MB(3,1) = -SB
      MB(3,3) = +CB
      MB(2,2) = 1D0
C
C----------------------------------- rotation by GAM around fixed z-axis
      SG = SIN(GAM*PI/180D0)
      CG = COS(GAM*PI/180D0)
C
      MG(1:3,1:3) = 0D0
      MG(1,1) = +CG
      MG(1,2) = -SG
      MG(2,1) = +SG
      MG(2,2) = +CG
      MG(3,3) = 1D0
C
C---------------------------------------------- combining the 3 rotations
C   MROT(alf,bet,gam)[at]  = M(alf,z)[af] * M(bet,y)[af] * M(gam,z)[af]
C
      MBG(1:3,1:3) = MATMUL(MB(1:3,1:3),MG(1:3,1:3))
C
      MROT(1:3,1:3) = MATMUL(MA(1:3,1:3),MBG(1:3,1:3))
C
      END
