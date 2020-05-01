C*==mod_constants.f    processed by SPAG 6.55Rc at 22:16 on 14 May 2007
      MODULE MOD_CONSTANTS
C   ********************************************************************
C   *                                                                  *
C   *             DON'T CHANGE THE FORMAT OF THIS FILE                 *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *  module to store all   CONSTANTS  in cgs units                   *
C   *                                                                  *
C   *         Fundamental Physical Constants --- Complete Listing      *
C   *                                                                  *
C   *           From:  http://physics.nist.gov/constants  21/12/16     *
C   *                                                                  *
C   *   Avogadro constant                                              *
C   *   NA_SI                       6.022 140 857 D+23     1 / mol     *
C   *   NA_CGS                                                         *
C   *                                                                  *
C   *   Bohr magneton                                                  *
C   *   MB_SI                       927.400 999 4 D-26     J / T       *
C   *   MB_CGS                      9.274 009 994 D-21     erg / G     *
C   *                                                                  *
C   *   Bohr radius                                                    *
C   *   A0_SI                       0.529 177 210 67 D-10  m           *
C   *   A0_CGS                      0.529 177 210 67 D-08  cm          *
C   *                                                                  *
C   *   nuclear magneton                                               *
C   *   MN_SI                       5.050 783 699 D-27     J / T       *
C   *                                                                  *
C   *   Boltzmann constant                                             *
C   *   KB_SI                       1.380 648 52 D-23      J / K       *
C   *   KB_CGS                      1.380 648 52 D-16      erg / K     *
C   *                                                                  *
C   *   electron mass                                                  *
C   *   M0_SI                       9.109 383 56 D-31      kg          *
C   *   M0_CGS                      0.910 938 356 D-27     g           *
C   *                                                                  *
C   *   proton mass                                                    *
C   *   MP_SI                       1.672 621 898 D-27     kg          *
C   *                                                                  *
C   *   electron volt                                                  *
C   *   EV_J                        1.602 176 6208 D-19    J           *
C   *   EV_ERG                      1.602 176 6208 D-12    erg         *
C   *                                                                  *
C   *   elementary charge                                              *
C   *   E0_SI                       1.602 176 6208 D-19    C           *
C   *                                                                  *
C   *   fine-structure constant                                        *
C   *   ALPHA_FS                    7.297 352 5664 D-3                 *
C   *                                                                  *
C   *                                                                  *
C   *   Planck constant over 2 pi                                      *
C   *   HBAR_SI                     1.054 571 800 D-34     J s         *
C   *   HBAR_CGS                    1.054 571 800 D-27     erg s       *
C   *                                                                  *
C   *   Rydberg constant times hc                                      *
C   *   RY_EV                       13.605 693 009         eV          *
C   *                                                                  *
C   *   speed of light in vacuum                                       *
C   *   C_SI                        299 792 458 (exact)    m / s       *
C   *   C_CGS                                                          *
C   *   C_AU_NIST                   2 / ALPHA_FS                       *
C   *                                                                  *
C   *   classical electron radius                                      *
C   *   RE_CGS = E0_CGS^2 / M0_CGS * C_CGS^2                           *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   CONVERSION FACTORS                                             *
C   *                                                                  *
C   *   convert J to erg                                               *
C   *   J_ERG                       1 D+7                  erg / J     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C                  frequently used numerical constants
C
C     Y00     spherical harmonics for L=(0,0)
C     F00SF   shape function for L=(0,0) in the muffin tin regime
C
C-----------------------------------------------------------------------
C
      COMPLEX*16 CI,C0,C1
      PARAMETER (CI=(0.0D0,1.0D0),C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))
C
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
      COMPLEX*16 CI2PI
      PARAMETER (CI2PI=CI*2*PI)
      REAL*8 CONST_2PI,CONST_4PI
      PARAMETER (CONST_2PI = 2.0D0*PI,CONST_4PI = 4.0D0*PI)
      REAL*8 SQRT_PI,SQRT_4PI
      PARAMETER (SQRT_PI = SQRT(PI),SQRT_4PI = SQRT(CONST_4PI))
      REAL*8 CONST_4PIOV3,SQRT_4PIOV3
      PARAMETER (CONST_4PIOV3 = CONST_4PI/3.0D0)
      PARAMETER (SQRT_4PIOV3 = SQRT(CONST_4PIOV3))
C      
      REAL*8 SQRT_2,Y00,F00SF
      PARAMETER (SQRT_2 = SQRT(2.0D0))
      PARAMETER (Y00=SQRT(1.0D0/(4.0D0*PI)),F00SF = SQRT_4PI)
C      
      REAL*8 DEG_ARC,ARC_DEG
      PARAMETER (DEG_ARC=PI/180.0D0,ARC_DEG=180.0D0/PI)
C
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C     fundamental constants in SI units from NIST table (see above)
C
      REAL*8 ALPHA_FS,A0_SI,
     &       C_SI,E0_SI,HBAR_SI,KB_SI,MB_SI,MN_SI,M0_SI,NA_SI,MP_SI
C
      PARAMETER ( ALPHA_FS  = 7.2973525664D-3,
     &            A0_SI     = 0.52917721067D-10,
     &            C_SI      = 299792458D0,
     &            E0_SI     = 1.6021766208D-19,
     &            HBAR_SI   = 1.054571800D-34,
     &            KB_SI     = 1.38064852D-23,
     &            MB_SI     = 927.4009994D-26,
     &            MN_SI     = 5.050783699D-27, 
     &            M0_SI     = 9.10938356D-31,
     &            MP_SI     = 1.672621898D-27,
     &            NA_SI     = 6.022140857D+23)
C
C     factors to change units
C
      REAL*8 J_ERG,EV_J,EV_ERG,RY_ERG,RY_EV,M_CM,T_GAUSS,M_ANG
C
      PARAMETER (   J_ERG    = 1D+7,
     &              EV_J     = E0_SI,
     &              RY_EV    = 13.605693009D0,
     &              EV_ERG   = EV_J*J_ERG,
     &              RY_ERG   = RY_EV*EV_ERG,
     &              M_CM     = 1D+2,
     &              M_ANG    = 1D+10,
     &              T_GAUSS  = 1D+4)
C
C     changing units for fundamental constants to be used in program
C
      REAL*8 A0_CGS,A0_ANG,C_AU,C_HU,C_CGS,E0_CGS,HBAR_CGS,KB_CGS,
     &       MB_CGS,MN_CGS,M0_CGS,NA_CGS,RE_CGS,HU_RY,HU_EV
C
      PARAMETER ( C_AU      = 2.0D0/ALPHA_FS,
     &            C_HU      = 1.0D0/ALPHA_FS,
     &            C_CGS     = C_SI*M_CM,
     &            HU_RY     = 2.0D0,
     &            HU_EV     = HU_RY*RY_EV, 
     &            A0_CGS    = A0_SI*M_CM,
     &            A0_ANG    = A0_SI*M_ANG,
     &            E0_CGS    = E0_SI*C_SI*10,
     &            HBAR_CGS  = HBAR_SI*J_ERG,
     &            KB_CGS    = KB_SI*J_ERG,
     &            MB_CGS    = MB_SI*J_ERG/T_GAUSS,
     &            MN_CGS    = MN_SI*J_ERG/T_GAUSS,
     &            M0_CGS    = M0_SI*1D+03,
     &            NA_CGS    = NA_SI,
     &            RE_CGS    = E0_CGS**2 / (M0_CGS * C_CGS**2)    
     &           )
C
C     factor to convert B-field from a.u. to cgs units i.e. to Gauss
C
      REAL*8 B_AU2CGS
      PARAMETER ( B_AU2CGS  = E0_CGS/(A0_CGS*A0_CGS) )
C
C     factor to convert susceptibility from a.u. to cm^3/mol
C
      REAL*8 CHI_AU2CGS
      PARAMETER (CHI_AU2CGS = NA_CGS*MB_CGS**2/RY_ERG )
C
      END
