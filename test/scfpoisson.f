C*==scfpoisson.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFPOISSON(R,V,Z,DR,XR,IRTOP)
C   ********************************************************************
C   *                                                                  *
C   *  this subroutine solves poisson equation by use of the           *
C   *  simpson's integration method. it is rather a straightfoward     *
C   *  way, though not so elegant. an advantage is the method is       *
C   *  quite tough and has never failed even in the cases of           *
C   *  singular charge distributions of muonic atoms.                  *
C   *                                                                  *
C   *  based on H. AKAI's routine <POISNA>                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS
      IMPLICIT NONE
C*--SCFPOISSON17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,Z
      REAL*8 DR(IRTOP),R(IRTOP),V(IRTOP),XR(IRTOP)
C
C Local variables
C
      REAL*8 F,FP,RNUC,RZ,W1(2000),W2(2000)
      INTEGER IR
      REAL*8 POT_NUC,RNUCTAB
C
C*** End of declarations rewritten by SPAG
C
      DO IR = 1,IRTOP
         W1(IR) = 2D0*XR(IR)*DR(IR)*R(IR)
         W2(IR) = W1(IR)*XR(IR)
      END DO
      V(1) = 0D0
      DO IR = 1,IRTOP - 2,2
         V(1) = V(1) + W1(IR) + 4D0*W1(IR+1) + W1(IR+2)
      END DO
      V(1) = V(1)/3D0
      F = 0D0
      DO IR = 2,IRTOP - 1,2
         FP = F + (5D0*W2(IR-1)-W2(IR+1)+8D0*W2(IR))/12D0
         V(IR) = V(IR-1) - F/XR(IR-1) + FP/XR(IR)
     &           - (5D0*W1(IR-1)-W1(IR+1)+8D0*W1(IR))/12D0
         FP = F + (W2(IR-1)+4D0*W2(IR)+W2(IR+1))/3D0
         V(IR+1) = V(IR-1) - F/XR(IR-1) + FP/XR(IR+1)
     &             - (W1(IR-1)+4D0*W1(IR)+W1(IR+1))/3D0
         F = FP
      END DO
C
      IF ( .NOT.FINITE_NUCLEUS .OR. Z.EQ.0 ) THEN
C
         DO IR = 1,IRTOP
            V(IR) = V(IR) - 2D0*Z/XR(IR)
         END DO
C
      ELSE
C
         RZ = DBLE(Z)
         RNUC = RNUCTAB(Z)
         DO IR = 1,IRTOP
            IF ( RNUC.GE.XR(IR) ) THEN
               V(IR) = V(IR) + POT_NUC(XR(IR),RZ)
            ELSE
               V(IR) = V(IR) - 2D0*Z/XR(IR)
            END IF
         END DO
C
      END IF
C
      END
