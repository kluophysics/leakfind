C*==metorque.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE METORQUE(DVSPIN,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the torque matrix elements              *
C   *  used for Gilbert damping                                        *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,AME_G,IOMT,ISMT,WKM3,WKM4
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--METORQUE14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DVSPIN(NKMMAX,NKMMAX,3),MSST(NKMMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 J(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE J
C
C-----------------------------------------------------------------------
      ALLOCATE (J(NKMMAX,NKMMAX,3))
C-----------------------------------------------------------------------
C
C-- calculate J(IPOL)  angular matrix elements of TOTAL angular momentum
C--------------------------------- spherical coordinates   (-), (0), (+)
C
      J(1:NKM,1:NKM,:) = AME_G(1:NKM,1:NKM,:,IOMT)
     &                   + 0.5D0*AME_G(1:NKM,1:NKM,:,ISMT)
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      CALL CMAT_CONVERT_POLAR(J,'S>C')
C
C----------------------------------- prepare derivatives uJ*msst-msst*uJ
C
      CALL CMATMUL(NKM,NKMMAX,J(1,1,1),MSST(1,1),WKM3)
      CALL CMATMUL(NKM,NKMMAX,MSST(1,1),J(1,1,1),WKM4)
C
      DVSPIN(1:NKM,1:NKM,1) = (WKM3(1:NKM,1:NKM)-WKM4(1:NKM,1:NKM))
C
      CALL CMATMUL(NKM,NKMMAX,J(1,1,2),MSST(1,1),WKM3)
      CALL CMATMUL(NKM,NKMMAX,MSST(1,1),J(1,1,2),WKM4)
C
      DVSPIN(1:NKM,1:NKM,2) = (WKM3(1:NKM,1:NKM)-WKM4(1:NKM,1:NKM))
C
      DVSPIN(1:NKM,1:NKM,3) = C0
C
      END
