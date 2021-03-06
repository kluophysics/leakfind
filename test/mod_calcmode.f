C*==mod_calcmode.f    processed by SPAG 6.70Rc at 15:56 on  5 Feb 2013
      MODULE MOD_CALCMODE
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables that control the                  *
C   *  calculation mode                                                *
C   *                                                                  *
C   *  ALL data initialized in  <INIT_MOD_CALCMODE>                    *
C   *                       or  <INITVAR>                              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C PARAMETER definitions
C
      LOGICAL BORN_USE_SPH
      PARAMETER (BORN_USE_SPH=.TRUE.)
      LOGICAL FP_USE_DIRECTIONS
      PARAMETER (FP_USE_DIRECTIONS=.FALSE.)
      LOGICAL MIGRATE
      PARAMETER (MIGRATE=.FALSE.)
      LOGICAL PUBLIC_VERSION
      PARAMETER (PUBLIC_VERSION=.FALSE.)
C      PARAMETER (PUBLIC_VERSION=.TRUE.)
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      LOGICAL BLCOUPL,BREITINT,DMFT,EXTFIELD,FULLCHARGE,GF_CONV_RH,LDAU,
     &        LHS_SOL_EQ_RHS_SOL,LLOYD,ME_CC_BRA_RWF,MOL,NONMAG,
     &        RELAX_CLU,ROTATE_LATTICE,SCL_ALAT,SEMICORE,SIGMA_LAYER,
     &        SIGMA_PROJECT,STORE_SFN,SYMCHECK,THERMAL_VIBRA_FLUCT,
     &        TUTORIAL,UPDATE_EFERMI,USE_CONST_POTENTIAL,
     &        USE_FPCHR_ZJOVZZ,U_CALCULATION,L_PROJECTED_ME,
     &        CALC_KINETIC_ENERGY_DENSITY,MOMENTS_ROTATED,
     &        CALC_EFG,CALC_FORCE,WAVE_FUNCTIONS_AVAILABLE,
     &        VMTZ_AVERAGE_SURFACE
      INTEGER BREAKPOINT,ICORE,IHYPER,IREL,IREL_HOST,ISMQHFI,ITEST,
     &        IT_CENTER_U_CALCULATION,IXRAY,KMROT,
     &        LIM_FPDIRBORN,LIM_FPRWFBORN
      CHARACTER*10 CLUTYPE,ORBPOL,PROGNAME,SELFENERGY,SOLVER,
     &             SOLVER_FP,TASK
      CHARACTER*30 KKRMODE,SETPOT
      CHARACTER*20 TXTKMROT(0:4)
      REAL*8 POTFIX,TOL_DIRBS,TOL_FPDIRBORN,TOL_FPDIRBS,TOL_FPRWFBORN,
     &       TOL_FPRWFBS,TOL_RWFBS,U_POT_SHIFT
      SAVE BLCOUPL,EXTFIELD,IREL_HOST,PROGNAME,SEMICORE,SOLVER,
     &     L_PROJECTED_ME,CALC_KINETIC_ENERGY_DENSITY
      SAVE CALC_EFG,CALC_FORCE,WAVE_FUNCTIONS_AVAILABLE
      SAVE VMTZ_AVERAGE_SURFACE,SETPOT
C     
C*** End of declarations rewritten by SPAG
C
      DATA KMROT/0/,MOMENTS_ROTATED/.FALSE./
      DATA ME_CC_BRA_RWF/.FALSE./
      DATA THERMAL_VIBRA_FLUCT/.FALSE./
      DATA TUTORIAL/.FALSE./,
     &     LHS_SOL_EQ_RHS_SOL/.TRUE./,L_PROJECTED_ME/.TRUE./
      DATA ORBPOL/'NONE'/,DMFT/.FALSE./
      DATA IREL/3/,SELFENERGY/'NONE'/,ICORE/1/,IHYPER/1/,ISMQHFI/0/
      DATA NONMAG/.TRUE./,IXRAY/0/,ITEST/0/,MOL/.FALSE./
      DATA LLOYD/.FALSE./,KKRMODE/'STANDARD-KKR'/
      DATA LDAU/.FALSE./,BREITINT/.FALSE./
      DATA CLUTYPE/'NONE'/,RELAX_CLU/.FALSE./,SOLVER_FP/'BORN'/
      DATA UPDATE_EFERMI/.TRUE./,SYMCHECK/.FALSE./
      DATA TOL_RWFBS/2.0D-8/,TOL_FPRWFBS/2.0D-8/
      DATA TOL_DIRBS/2.0D-8/,TOL_FPDIRBS/2.0D-8/
      DATA TOL_FPDIRBORN/1.0D-10/,TOL_FPRWFBORN/1.0D-10/
      DATA LIM_FPDIRBORN/10/,LIM_FPRWFBORN/10/
      DATA USE_FPCHR_ZJOVZZ/.FALSE./,SCL_ALAT/.FALSE./
      DATA SIGMA_PROJECT/.FALSE./
      DATA STORE_SFN/.FALSE./,GF_CONV_RH/.FALSE./
      DATA TASK/'NONE'/,POTFIX/0D0/,USE_CONST_POTENTIAL/.FALSE./
      DATA FULLCHARGE/.FALSE./,BREAKPOINT/0/
      DATA SIGMA_LAYER/.FALSE./,ROTATE_LATTICE/.FALSE./
      DATA U_CALCULATION/.FALSE./,U_POT_SHIFT/0D0/
      DATA IT_CENTER_U_CALCULATION/0/
      DATA CALC_KINETIC_ENERGY_DENSITY/.FALSE./
      DATA TXTKMROT/'no M-rotation       ','M-rotation IQ-dept. ',
     &     'global M-rotation   ','spin spiral  TET=90 ',
     &     'spin spiral  TET<>90'/
      DATA CALC_EFG/.FALSE./,CALC_FORCE/.FALSE./
      DATA WAVE_FUNCTIONS_AVAILABLE/.FALSE./
      DATA VMTZ_AVERAGE_SURFACE/.FALSE./
      DATA SETPOT/'NO'/
      END
