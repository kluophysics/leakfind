C*==mod_cpa.f    processed by SPAG 6.55Rc at 21:28 on  3 Apr 2007
      MODULE MOD_CPA
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tvariables used for   CPA   calculations    *
C   *                                                                  *
C   *  ALL data initialized in  <INIT_MOD_CPA>                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--MOD_CPA12
      REAL*8 CPAMIX,CPATOL,ALPHASRO
      INTEGER NCPA,CPALVL,CPASTART,ICPAALG,ITCPAMAX
      LOGICAL NLCPAWRCFG,USENLCPA,WRCPA
      CHARACTER*80 CPAFIL
C
      DATA CPALVL/1/,USENLCPA/.FALSE./,NLCPAWRCFG/.FALSE./
      DATA CPAMIX/1D0/,ALPHASRO/0D0/,WRCPA/.FALSE./
C
      SAVE CPAFIL,CPAMIX,CPATOL,CPALVL,CPASTART,ICPAALG,ITCPAMAX,
     &   NLCPAWRCFG,USENLCPA,WRCPA,ALPHASRO,NCPA
C
      END
