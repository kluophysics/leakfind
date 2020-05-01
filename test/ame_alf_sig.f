C*==ame_alf_sig.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AME_ALF_SIG
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements  needed  for the         *
C   *                                                                  *
C   *   Bargmann-Wigner spin-polarization operator:                    *
C   *                                                                  *
C   *   \beta\Sigma_ipolspin \alpha_ipol                               *
C   *                                                                  *
C   *   NB: They are prepared such that, after CMAT_CONVERT_POLAR      *
C   *       from spherical to cartesian coordinates is performed,      *
C   *       the proper MEs are obtained. This is done to stay          *
C   *       consistently within the spherical coordinates.             *
C   *                                                                  *
C   *   ipol     = 1,2,3  ==  (-),(0),(+):                             *
C   *                                                                  *
C   *   NB!: The following convention is used for transforming from    *
C   *        cartesian to (special) spherical coordinates              *
C   *                                                                  *
C   *  ->        ->       ->                                           *
C   *  e(+) = - [e(x) + i e(y)] / sqrt(2) left circ. positive helicity *
C   *                                                                  *
C   *  ->        ->       ->                                           *
C   *  e(-) = + [e(x) - i e(y)] / sqrt(2) right circ. negative helicity*
C   *                                                                  *
C   *  DK, 08/2013, 01/2015                                            *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKMP_EXT,A1_SPIN_CURR_ALF,AME_G,ISMT,IMKM_IKM,
     &    A2_SPIN_CURR_ALF,NPOL
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_2
      IMPLICIT NONE
C*--AME_ALF_SIG33
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 AMESIGA1(:,:,:),AMESIGA2(:,:,:),AME_NU_NU(:,:),WZ2
      INTEGER M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AME_NU_NU,AMESIGA1,AMESIGA2
C
C=======================================================================
C angular matrix elements for ALFA-related term in spin current density
C=======================================================================
C
      IF ( .NOT.ALLOCATED(A1_SPIN_CURR_ALF) ) THEN
         ALLOCATE (A1_SPIN_CURR_ALF(NKMP_EXT,NKMP_EXT,NPOL,NPOL))
         ALLOCATE (A2_SPIN_CURR_ALF(NKMP_EXT,NKMP_EXT,NPOL,NPOL))
      END IF
C
      A1_SPIN_CURR_ALF(:,:,:,:) = C0
      A2_SPIN_CURR_ALF(:,:,:,:) = C0
C
C=======================================================================
C                          work space
C=======================================================================
C
      ALLOCATE (AMESIGA1(NKMP_EXT,NKMP_EXT,NPOL))
      ALLOCATE (AMESIGA2(NKMP_EXT,NKMP_EXT,NPOL))
C
      ALLOCATE (AME_NU_NU(NKMP_EXT,NKMP_EXT))
C
      WZ2 = SQRT_2
C
      AMESIGA1(:,:,:) = 0D0
      AMESIGA2(:,:,:) = 0D0
      AME_NU_NU(:,:) = 0D0
C
C=======================================================================
C
      DO N = 1,NKMP_EXT
         DO M = 1,NKMP_EXT
            AMESIGA1(N,M,1) = AME_G(N,IMKM_IKM(M),3,ISMT)
            AMESIGA1(N,M,2) = AME_G(N,IMKM_IKM(M),1,ISMT)
            AMESIGA1(N,M,3) = AME_G(N,IMKM_IKM(M),2,ISMT)
            AMESIGA2(N,M,1) = AME_G(IMKM_IKM(N),M,3,ISMT)
            AMESIGA2(N,M,2) = AME_G(IMKM_IKM(N),M,1,ISMT)
            AMESIGA2(N,M,3) = AME_G(IMKM_IKM(N),M,2,ISMT)
            IF ( N.EQ.IMKM_IKM(M) ) AME_NU_NU(N,M) = 1D0
         END DO
      END DO
C
      AME_NU_NU(:,:) = AME_NU_NU(:,:)/WZ2
C
C ----------------------------------------------------------------------
C ---- refined/redefined spin-polarization operator --------------------
C ---- set xx,yy,zz component to zero ----------------------------------
C ---- this takes the definition of the Bargmann-Wigner / Vernes
C ---- inspired form and modifies it
C ---- overall: real parts of tensors are unaffected, imaginary parts
C ------------- now vanish
      AME_NU_NU(:,:) = 0D0
C ----------------------------------------------------------------------
C
      N = NKMP_EXT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ------------- spherical components contributing to x-spin-polarization
C
      A1_SPIN_CURR_ALF(1:N,1:N,3,1) = AMESIGA1(1:N,1:N,3)
     &                                /WZ2 - AME_NU_NU(1:N,1:N)
      A1_SPIN_CURR_ALF(1:N,1:N,1,1) = AMESIGA1(1:N,1:N,3)
     &                                /WZ2 + AME_NU_NU(1:N,1:N)
C
      A2_SPIN_CURR_ALF(1:N,1:N,1,1) = -AMESIGA2(1:N,1:N,3)
     &                                /WZ2 - AME_NU_NU(1:N,1:N)
C
      A2_SPIN_CURR_ALF(1:N,1:N,3,1) = -AMESIGA2(1:N,1:N,3)
     &                                /WZ2 + AME_NU_NU(1:N,1:N)
C
      A1_SPIN_CURR_ALF(1:N,1:N,2,1) = (AMESIGA1(1:N,1:N,2)+AMESIGA1(1:N,
     &                                1:N,1))/WZ2
C
      A2_SPIN_CURR_ALF(1:N,1:N,2,1) = (-AMESIGA2(1:N,1:N,1)-AMESIGA2(1:N
     &                                ,1:N,2))/WZ2
C
C ----------------------------------------------------------------------
C ------------- spherical components contributing to y-spin-polarization
C
      A1_SPIN_CURR_ALF(1:N,1:N,3,2) = -AMESIGA1(1:N,1:N,3)
     &                                /WZ2/CI - AME_NU_NU(1:N,1:N)*CI
      A1_SPIN_CURR_ALF(1:N,1:N,1,2) = AMESIGA1(1:N,1:N,3)
     &                                /WZ2/CI - AME_NU_NU(1:N,1:N)*CI
C
      A2_SPIN_CURR_ALF(1:N,1:N,3,2) = AMESIGA2(1:N,1:N,3)
     &                                /WZ2/CI + AME_NU_NU(1:N,1:N)*CI
C
      A2_SPIN_CURR_ALF(1:N,1:N,1,2) = -AMESIGA2(1:N,1:N,3)
     &                                /WZ2/CI + AME_NU_NU(1:N,1:N)*CI
C
C
      A1_SPIN_CURR_ALF(1:N,1:N,2,2) = (-AMESIGA1(1:N,1:N,2)+AMESIGA1(1:N
     &                                ,1:N,1))/(WZ2*CI)
C
      A2_SPIN_CURR_ALF(1:N,1:N,2,2) = (AMESIGA2(1:N,1:N,2)-AMESIGA2(1:N,
     &                                1:N,1))/(WZ2*CI)
C ----------------------------------------------------------------------
C ------------- spherical components contributing to z-spin-polarization
C
      A1_SPIN_CURR_ALF(1:N,1:N,3,3) = AMESIGA1(1:N,1:N,1)
      A1_SPIN_CURR_ALF(1:N,1:N,1,3) = -AMESIGA1(1:N,1:N,2)
C
      A2_SPIN_CURR_ALF(1:N,1:N,3,3) = -AMESIGA2(1:N,1:N,1)
      A2_SPIN_CURR_ALF(1:N,1:N,1,3) = AMESIGA2(1:N,1:N,2)
C
      A1_SPIN_CURR_ALF(1:N,1:N,2,3) = AME_NU_NU(1:N,1:N)*WZ2
      A2_SPIN_CURR_ALF(1:N,1:N,2,3) = -AME_NU_NU(1:N,1:N)*WZ2
C
      END
