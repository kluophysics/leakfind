C*==init_mod_files.f    processed by SPAG 6.70Rc at 03:10 on 16 Jun 2015
C     *==init_mod_files.f    processed by SPAG 6.70Rc at 08:45 on 27 Mar 2014
      SUBROUTINE INIT_MOD_FILES
C   ********************************************************************
C   *                                                                  *
C   *  initialize some variables needed for string processing          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:ICHAR_UCA,ICHAR_UCZ,ICHAR_LCA,ICHAR_LCZ,
     &    ICHAR_0,ICHAR_9,ICHAR_BLANK,ICHAR_UNDERSCORE,IFILCBWF0,
     &    IFILCBWF
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
      ICHAR_UCA = ICHAR('A')
      ICHAR_UCZ = ICHAR('Z')
      ICHAR_LCA = ICHAR('a')
      ICHAR_LCZ = ICHAR('z')
      ICHAR_0 = ICHAR('0')
      ICHAR_9 = ICHAR('9')
      ICHAR_BLANK = ICHAR(' ')
      ICHAR_UNDERSCORE = ICHAR('_')
C
C----------------------------------------------------------------------- 
C                 see MOD_FILES for the setting scheme
C----------------------------------------------------------------------- 
C
      IFILCBWF = IFILCBWF0 + 1
C
      END
C*==set_ifil_lhs.f    processed by SPAG 6.70Rc at 03:10 on 16 Jun 2015
      SUBROUTINE SET_IFIL_LHS(IFIL_RHS,IFIL_LHS)
C   ********************************************************************
C   *                                                                  *
C   *         set the chanel IFIL_LHS for the LHS wave functions       *
C   *  related to the chanel IFIL_RHS for the RHS wave functions       *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      USE MOD_FILES,ONLY:IFILCORWF,IFILGFWF,IFILLDAU
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL_LHS,IFIL_RHS
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------------------------
C     core wave function or reference wave function (LDA+U, DMFT)
C------------------------------------------------------------------------
C
      IF ( (IFIL_RHS.EQ.IFILCORWF) .OR. (IFIL_RHS.EQ.IFILGFWF) .OR. 
     &     (IFIL_RHS.EQ.IFILLDAU) ) THEN
C
         IFIL_LHS = IFIL_RHS
C
         RETURN
C
      END IF
C
C------------------------------------------------------------------------
C              ordinary conduction band wave function
C------------------------------------------------------------------------
C
      IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
         IFIL_LHS = IFIL_RHS
C
      ELSE
C
         IFIL_LHS = IFIL_RHS + 1
C
      END IF
C
      END
C*==set_ifil_sph.f    processed by SPAG 6.70Rc at 03:10 on 16 Jun 2015
      SUBROUTINE SET_IFIL_SPH(IFIL_RHS,IFIL_SPH)
C   ********************************************************************
C   *                                                                  *
C   *         set the chanel IFIL_SPH for the SPH wave functions       *
C   *  related to the chanel IFIL_RHS for the RHS wave functions       *
C   *                                                                  *
C   *    this chanel is used for the spherical wave functions          *
C   *    if the BORN solver is used. For the spherical solutions       *
C   *    one has always  LHS = RHS.                                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL_RHS,IFIL_SPH
C
C*** End of declarations rewritten by SPAG
C
      IFIL_SPH = IFIL_RHS + 2
C
      END
