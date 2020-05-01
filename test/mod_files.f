C*==mod_files.f    processed by SPAG 6.70Rc at 15:46 on 24 Sep 2012
      MODULE MOD_FILES
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables connected with files              *
C   *                                                                  *
C   *  the variables IFIL... specify the reserved I/O chanel           *
C   *                                                                  *
C   ********************************************************************
C
C   IFILPOT_CLU   3    potential file for embedded cluster
C   IFILPOT       4    potential file used for input and output
C   IFILINP       5    input file with instructions set by user
C   IFILTAU       9    file to write/read site-diagonal tau
C   IFILTAUIJ     8    file to write/read site-off-diagonal tau
C   IFILDOS      10    DOS file
C   IFILNOS      11    NOS file
C   IFILDAT      12    DAT file to tabulated results
C   IFILUCALC    13    DAT file to calculate U parameter
C   IFILFLUCT    16    DAT file to import T-dependent spin fluctuations
C                      from Monte Carlo simulations (FLUCT) or
C                      other sources       
C   IFILMEZZL    52    file to pass MEZZL,MEZJL 
C                      <FPSSITE> -> <FPCALCDOS>, <FPCHRDNS>
C   IFILLOG      53    SCF-log file
C   IFILSCLALAT  54    log file for SCL_ALAT
C   IFILBREAK    57    scratch file used for BREAK
C   IFILBUILDBOT 58    file used for BUILDBOT tests
C
C   --------------------------------------------------------------------
C   Files connected with SPEC routines
C   reserved range  60 - 69
C   --------------------------------------------------------------------
C
C   IFILSPECOU1  60 output for results
C   IFILSPECOU2  61 output for results
C   IFILSPECINP  62 input file (old format)
C   IFILSPECSTR  63 structure  (old format)
C   IFILSPECBAR  64 barier     (old format)
C   IFILSPECPOT  65 potential  (old format)
C   IFILSPECIO1  66 temporary file
C   IFILSPECIO2  67 temporary file
C   IFILSPECIO3  68 temporary file
C   IFILSPECOU3  69 output file (old format)
C
C   --------------------------------------------------------------------
C   Files connected with wave functions
C   reserved range  1000 - 1999
C   --------------------------------------------------------------------
C
C   IFILWF0      1000
C   IFILCORWF    IFILWF0 + 1    core wave function
C   IFILGFWF     IFILWF0 + 5    reference wave function for GF-matrix
C   IFILLDAU     IFILWF0 + 6    reference wave function for LDA+U
C
C   IFILCBWF0    1100
C
C   standard case: 1 energy z
C
C   IFILCBWF     IFILCBWF0 + 1  RHS valence band wave function
C   IFILCBWF_LHS IFILCBWF0 + 2  LHS valence band wave function
C                               left hand side solution
C   IFILCBWF_SPH IFILCBWF0 + 3  spherical solution in full potential case
C
C   special case: several energies IE simultanously only E
C
C   IFIL = IFILCBWF0 + 10*IE + 1 RHS for E 
C                            + 2 LHS for E  
C                            + 3 SPH for E  
C
C        >> IFILCBWF_INCREMENT_ENERGY=10
C
C   special case: several energies IE simultanously for  E  and  E*
C
C   IFIL = IFILCBWF0 + 10*IE + 1 RHS for E
C                            + 2 LHS for E 
C                            + 3 SPH for E 
C                            + 4 RHS for E*
C                            + 5 LHS for E*  
C                            + 6 SPH for E*  
C
C        >> IFILCBWF_INCREMENT_CC=3
C
C   NOTE     IFILCBWF is used as variable in OPM routines
C            ----->  init in DATA statement below
C
C   --------------------------------------------------------------------
C   temporary files
C   reserved range  80 - 89
C   --------------------------------------------------------------------
C
C   IOTMP        80    temporary I/O file
C   IOTMP1       81    temporary I/O file
C   IOTMP2       82    temporary I/O file
C   IOTMP3       83    temporary I/O file
C   IOTMP4       84    temporary I/O file
C
C   --------------------------------------------------------------------
C   Files used for application of the Broyden method
C   --------------------------------------------------------------------
C
C   IFILBROY_RHOPOT  86 iteration of potential or charge
C   IFILBROY_WEISS   87 iteration of Weiss field
C   IFILBROY_ANG     88 iteration of Weiss field
C   IFILBROY_CPA     89 iteration of CPA
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LCMDMAX,LSTRMAX,IFILPOT_CLU,IFILPOT,IFILINP,IFILTAU,
     &        IFILDOS,IFILNOS,IFILDAT,IFILUCALC,IFILLOG,IFILSCLALAT,
     &        IFILSPECOU1,IFILSPECOU2,IFILSPECINP,IFILSPECSTR,
     &        IFILSPECBAR,IFILSPECPOT,IFILSPECIO1,IFILSPECIO2,
     &        IFILSPECIO3,IFILSPECOU3,IFILCORWF,IFILCBWF_LHS,
     &        IFILCBWF_SPH,IFILGFWF,IFILLDAU,IOTMP,IOTMP1,IOTMP2,IOTMP3,
     &        IOTMP4,IFILMEZZL,IFILBREAK,IFILFLUCT,IFILCBWF0,
     &        IFILCBWF_INCREMENT_ENERGY,IFILCBWF_INCREMENT_CC,
     &        IFILBUILDBOT,LCMD
      INTEGER IFILBROY_RHOPOT,IFILBROY_WEISS,IFILBROY_ANG,IFILBROY_CPA
      CHARACTER*10 PACKAGE
!REC increased LCMDMAX from 500
      PARAMETER (LCMDMAX=9999,LSTRMAX=80,IFILPOT_CLU=3,IFILPOT=4,
     &           Ifilinp=5,IFILTAU=9,IFILDOS=10,IFILNOS=11,IFILDAT=12,
     &           IFILUCALC=13,IFILFLUCT=16,
     &           IFILMEZZL=52,IFILLOG=53,IFILSCLALAT=54,
     &           IFILBREAK=57,IFILBUILDBOT=58,
     &           IFILSPECOU1=60,
     &           IFILSPECOU2=61,IFILSPECINP=62,IFILSPECSTR=63,
     &           IFILSPECBAR=64,IFILSPECPOT=65,IFILSPECIO1=66,
     &           IFILSPECIO2=67,IFILSPECIO3=68,IFILSPECOU3=69,
     &           IFILCORWF=71,
     &           IFILGFWF=75,IFILLDAU=76,IOTMP=80,IOTMP1=81,IOTMP2=82,
     &           IOTMP3=83,IOTMP4=84,
     &           IFILBROY_RHOPOT=86,IFILBROY_WEISS=87,
     &           IFILBROY_ANG=88,IFILBROY_CPA=89,
     &           IFILCBWF0=1100,
     &           IFILCBWF_INCREMENT_ENERGY=10,
     &     IFILCBWF_INCREMENT_CC=3)
C     
C Local variables
C
      INTEGER ICHAR_UNDERSCORE,IDUMMY,IFILDUMP,IFILTAUIJ,
     &        IFILTAUIJ_TMP,IPRINT,LDATFIL,LDATSET,LDATSET0,
     &        LINFO,LPOTFIL,LPOTFIL_CLU,LRECREAL8,LSFNFIL,LSYSTEM,
     &        LTAUIJFIL,LTITLE,RECLNGWF,RECLNG_TAUIJ,LFLUCTFIL,
     &        IFILCBWF,LFILNAM,POTFMTOUT,ICHAR_0,NDUMMY,RECLNGWF_SPH,
     &        ICHAR_9,ICHAR_BLANK,ICHAR_LCA,ICHAR_LCZ,ICHAR_UCA,
     &        ICHAR_UCZ,WALL_TIME_PROGRAM_START,WALL_TIME_LAST_CALL,
     &        IPOS_KEYWORD,N_FOUND
      COMPLEX*16 CDUMMY
      CHARACTER*(LCMDMAX) CMD00,CMDUC
      CHARACTER*80 DATFIL,DATSET,DATSET0,DOSFIL,INFO,NOSFIL,POTFIL,
     &             POTFIL_CLU,SFNFIL,SYSTEM,TAUFIL,TAUIJFIL,TITLE,
     &             FILNAM,FLUCTFIL
      LOGICAL DEBUG,LDUMMY,NOWRDOS,PLOT2DPR(0:3),PLOTPRS(0:3)
     &        ,RDDOS,RDTAU,RDTAUMQ,
     &        WRKAPDOS,WRLOG,WRMAT,WRNOS,WRPOLAR,WRPOT,WRTAU,WRTAUIJ,
     &        WRTAUMQ,WR_TRACK_INFO,CALC_JALF,CALC_JORB,
     &        PLOT_JALF,PLOT_JORB,WRBUILDBOT
      LOGICAL FOUND_SECTION,FOUND_KEYWORD,FOUND_REAL,FOUND_INTEGER,
     &        FOUND_STRING,FOUND_REAL_ARRAY,FOUND_INTEGER_ARRAY,
     &        FOUND_STRING_ARRAY,FOUND_QUANTUM_NUMBER
      REAL*8 CPU_TIME_PROGRAM_START,CPU_TIME_LAST_CALL,RDUMMY
      CHARACTER*10 S10DUMMY
      CHARACTER*(LSTRMAX) STRINP
      CHARACTER*1 CPOL_CART(3)
      CHARACTER*1 CPOL_SPHR(3)
      CHARACTER*255 GIT_HASH,GIT_COMPILE_DATE,GIT_BRANCH

      SAVE CDUMMY,CMD00,CMDUC,DATFIL,DATSET,DATSET0,DOSFIL,ICHAR_0,
     &     ICHAR_9,ICHAR_BLANK,ICHAR_LCA,ICHAR_LCZ,ICHAR_UCA,ICHAR_UCZ,
     &     ICHAR_UNDERSCORE,IDUMMY,INFO,LDATFIL,LDATSET,LDATSET0,LDUMMY,
     &     LFILNAM,LINFO,LPOTFIL,LPOTFIL_CLU,LRECREAL8,LSFNFIL,LSYSTEM,
     &     LTAUIJFIL,LTITLE,NDUMMY,NOSFIL,POTFIL,POTFIL_CLU,RDTAU,
     &     RDTAUMQ,RDUMMY,RECLNGWF,RECLNG_TAUIJ,S10DUMMY,SFNFIL,STRINP,
     &     SYSTEM,TAUFIL,TAUIJFIL,TITLE,WRLOG,WRMAT,WRTAU,WRTAUMQ,
     &     WR_TRACK_INFO,CPOL_CART,CPOL_SPHR,
     &     PLOT_JALF,CALC_JALF,CALC_JORB,PLOT_JORB,FLUCTFIL,LFLUCTFIL,
     &     GIT_HASH,GIT_COMPILE_DATE,GIT_BRANCH,
     &     IFILCBWF_SPH,IFILCBWF_LHS,WRBUILDBOT,PACKAGE,LCMD
      SAVE CPU_TIME_PROGRAM_START,CPU_TIME_LAST_CALL
      SAVE WALL_TIME_PROGRAM_START,WALL_TIME_LAST_CALL
      SAVE FOUND_SECTION,FOUND_KEYWORD,FOUND_REAL,FOUND_INTEGER,
     &     FOUND_STRING,FOUND_REAL_ARRAY,FOUND_INTEGER_ARRAY,
     &     FOUND_STRING_ARRAY,FOUND_QUANTUM_NUMBER,
     &        IPOS_KEYWORD,N_FOUND
C
C*** End of declarations rewritten by SPAG
C
      DATA IPRINT/0/,NOWRDOS/.FALSE./,WRKAPDOS/.FALSE./,WRPOLAR/.FALSE./
      DATA RDDOS/.FALSE./,WRPOT/.TRUE./,WRNOS/.FALSE./
      DATA WR_TRACK_INFO/.FALSE./,WRBUILDBOT/.TRUE./
      DATA PLOTPRS/.FALSE.,.FALSE.,.FALSE.,.FALSE./,WRTAUIJ/.FALSE./
      DATA PLOT2DPR/.FALSE.,.FALSE.,.FALSE.,.FALSE./,POTFMTOUT/8/
      DATA CALC_JALF/.FALSE./,CALC_JORB/.FALSE./
      DATA PLOT_JALF/.FALSE./,PLOT_JORB/.FALSE./
      DATA IFILTAUIJ_TMP/18/,IFILTAUIJ/8/,IFILDUMP/900/,DEBUG/.FALSE./
      DATA RECLNGWF_SPH/0/
      DATA IFILCBWF/1101/,LFLUCTFIL/0/
      DATA CPOL_CART/'x','y','z'/
      DATA CPOL_SPHR/'-','0','+'/
      DATA GIT_HASH/' not initialized'/
      DATA GIT_BRANCH/' not initialized'/
      DATA GIT_COMPILE_DATE/' not initialized'/
      DATA PACKAGE/'SPR-KKR   '/
C     
      END
