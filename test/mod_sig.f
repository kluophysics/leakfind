C*==mod_sig.f    processed by SPAG 6.70Rc at 08:36 on 20 Jul 2015
      MODULE MOD_SIG
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NSPINPROJ
      PARAMETER (NSPINPROJ=11)
      INTEGER IRESPONSE_ADA_ADA
      INTEGER IRESPONSE_SOT,IRESPONSE_EDELSTEIN
      PARAMETER (IRESPONSE_ADA_ADA=1)
      PARAMETER (IRESPONSE_SOT=8,IRESPONSE_EDELSTEIN=10)
C
C Local variables
C
      REAL*8 CONSI,CONSI2,EESI,EESI2,EE_PREFAC_AU,SIG_PREFAC_AU,
     &       SOTSI,SOTSI2,SOT_PREFAC_AU
      LOGICAL FLAG_ISPINPROJ(:),LCHECKME_SPIN,LJTOZ,LSP_NO_ALPHA,
     &        LSP_NO_NABLA
      INTEGER LIST_ISPR(:),NKTABSIG,NSPR
      CHARACTER*10 SIG_MODE
      CHARACTER*80 STR_ISP_PROJ(11)
      SAVE CONSI,CONSI2,EESI,EESI2,EE_PREFAC_AU,FLAG_ISPINPROJ,LIST_ISPR
     &     ,NKTABSIG,NSPR,SIG_MODE,SIG_PREFAC_AU,SOTSI,SOTSI2,
     &     SOT_PREFAC_AU
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE LIST_ISPR,FLAG_ISPINPROJ
C
C |-------+------------------------------------------------------|
C | IRESP | type                                                 |
C |-------+------------------------------------------------------|
C |     1 | alpha-alpha                                          |
C |     2 | Bargmann-Wigner Polarization  x                      |
C |     3 | Bargmann-Wigner Polarization  y                      |
C |     4 | Bargmann-Wigner Polarization  z                      |
C |     5 | Bargmann-Wigner Polarization  x, Nabla in spinpol    |
C |     6 | Bargmann-Wigner Polarization  y, Nabla in spinpol    |
C |     7 | Bargmann-Wigner Polarization  z, Nabla in spinpol    |
C |     8 | SOT spin-orbit torque                                |
C |     9 | Nabla-Nabla                                          |
C |    10 | sigma-alpha (Edelstein effect)                       |
C |-------+------------------------------------------------------|
C
      DATA STR_ISP_PROJ/
     &     ' alpha   alpha ',
     &     ' Polarization Bargmann-Wigner x',
     &     ' Polarization Bargmann-Wigner y',
     &     ' Polarization Bargmann-Wigner z',
     &     ' Polarization Bargmann-Wigner x,  Nabla',
     &     ' Polarization Bargmann-Wigner y,  Nabla ',
     &     ' Polarization Bargmann-Wigner z,  Nabla',
     &     ' SOT spin-orbit torque',' Nabla   Nabla',
     &     ' sigma   alpha   (Edelstein effect)',
     &     ' '/
C
      DATA LCHECKME_SPIN/.FALSE./
      DATA LSP_NO_ALPHA/.FALSE./
      DATA LSP_NO_NABLA/.FALSE./
      DATA LJTOZ/.FALSE./
      END
