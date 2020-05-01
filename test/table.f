C*==tabmuqin.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TABMUQIN(Z,A,MUNUC,QNUC,INUC,IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *   get nuclear magnetic moment in units of the nuclear magneton   *
C   *                  MUNUC = mu/mu_N                                 *
C   *   get nuclear quadrupole moment QNUC in barn                     *
C   *   get nuclear spin INUC                                          *
C   *   taken from Lechner, Physikalisch-chemische Daten D'Ans/Lax-    *
C   *   Taschenbuch fuer Chemiker und Physiker Bd.1, 4.Aufl.),         *
C   *   Springer, Berlin, 1992.                                        *
C   *   For input atomic number Z, the isotope listed has mass         *
C   *   number A and is the one with the highest natural abundance     *
C   *   (with the longest lifetime, resp.)                             *
C   *   that has a non-vanishing magnetic moment.                      *
C   *                                                  HF 05/10/95     *
C   ********************************************************************
      USE MOD_TABLES,ONLY:TAB_CHSYM
      IMPLICIT NONE
C*--TABMUQIN20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER A,IPRINT,Z
      REAL*8 INUC,MUNUC,QNUC
C
C Local variables
C
      INTEGER ATAB(0:110),I
      CHARACTER*2 CHSYM
      REAL*8 ITAB(0:110),MUTAB(0:110),QTAB(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA ATAB/0,1,3,7,9,11,13,14,17,19,21,23,25,27,29,31,33,35,39,39,
     &     43,45,47,51,53,55,57,59,61,63,67,69,73,75,77,81,83,85,87,89,
     &     91,93,97,99,101,103,105,109,111,115,119,121,125,127,129,133,
     &     137,139,0,141,143,147,147,153,157,159,163,165,167,169,173,
     &     175,179,181,183,187,189,193,195,197,199,205,207,209,0,0,0,0,
     &     0,227,229,231,235,0,239,241,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      DATA ITAB/0.0D0,0.5D0,0.5D0,1.5D0,1.5D0,1.5D0,0.5D0,1.0D0,2.5D0,
     &     0.5D0,1.5D0,1.5D0,2.5D0,2.5D0,0.5D0,0.5D0,1.5D0,1.5D0,3.5D0,
     &     1.5D0,3.5D0,3.5D0,2.5D0,3.5D0,1.5D0,2.5D0,0.5D0,3.5D0,1.5D0,
     &     1.5D0,2.5D0,1.5D0,4.5D0,1.5D0,0.5D0,1.5D0,4.5D0,2.5D0,0.5D0,
     &     0.5D0,2.5D0,4.5D0,2.5D0,4.5D0,2.5D0,0.5D0,2.5D0,0.5D0,0.5D0,
     &     4.5D0,0.5D0,2.5D0,0.5D0,2.5D0,0.5D0,3.5D0,1.5D0,3.5D0,0.0D0,
     &     2.5D0,3.5D0,3.5D0,3.5D0,2.5D0,1.5D0,1.5D0,2.5D0,3.5D0,3.5D0,
     &     0.5D0,2.5D0,3.5D0,4.5D0,3.5D0,0.5D0,2.5D0,1.5D0,1.5D0,0.5D0,
     &     1.5D0,0.5D0,0.5D0,0.5D0,4.5D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     1.5D0,2.5D0,1.5D0,3.5D0,0.0D0,0.5D0,2.5D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0/
      DATA (MUTAB(I),I=1,40)/2.792846D0, - 2.127624D0,3.256424D0,
     &      - 1.177800D0,2.688637D0,0.702411D0,0.403761D0, - 1.893790D0,
     &      2.628866D0, - 0.661796D0, - 2.217520D0, - 0.855450D0,
     &      3.641504D0, - 0.555290D0,1.131600D0,0.643821D0,0.821874D0,
     &      - 1.300000D0,0.391466D0, - 1.317260D0,4.756483D0,
     &      - 0.788480D0,5.151400D0, - 0.474540D0,3.453200D0,0.090623D0,
     &      4.627000D0, - 0.750020D0,2.223300D0,0.875478D0,2.016590D0,
     &      - 0.879467D0,1.439470D0,0.535060D0,2.270560D0, - 0.970669D0,
     &      1.353030D0, - 1.092820D0, - 0.137415D0, - 1.303620D0/
      DATA (MUTAB(I),I=41,80)/6.170500D0, - 0.933500D0,5.684700D0,
     &      - 0.718800D0, - 0.088400D0, - 0.642000D0, - 0.130691D0,
     &      - 0.594886D0,5.540800D0, - 1.047280D0,3.363400D0,
     &      - 0.888280D0,2.813270D0, - 0.777976D0,2.582023D0,0.937365D0,
     &      2.783200D0,0.0D0,4.136000D0, - 1.065000D0,2.580000D0,
     &      - 0.814800D0,1.533000D0, - 0.339800D0,2.014000D0,0.672600D0,
     &      4.173000D0, - 0.566500D0, - 0.231600D0, - 0.679890D0,
     &      2.232700D0, - 0.640900D0,2.370000D0,0.117785D0,3.219700D0,
     &      0.659933D0,0.159100D0,0.609490D0,0.148158D0,0.505885D0/
      DATA (MUTAB(I),I=81,110)/1.638213D0,0.582190D0,4.110600D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,1.100000D0,0.460000D0,2.310000D0,
     &      - 0.350000D0,0.0D0,0.203000D0,1.610000D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0/
      DATA (QTAB(I),I=1,40)/0.0D0,0.0D0, - 0.040000D0,0.053000D0,
     &      0.040000D0,0.0D0,0.019300D0, - 0.026000D0,0.0D0,0.102900D0,
     &      0.108000D0,0.220000D0,0.150000D0,0.0D0,0.0D0, - 0.064000D0,
     &      - 0.082490D0,0.0D0,0.054000D0,0.230000D0, - 0.220000D0,
     &      0.290000D0, - 0.051500D0, - 0.028500D0,0.330000D0,0.0D0,
     &      0.420000D0,0.162000D0, - 0.222000D0,0.150000D0,0.168000D0,
     &      - 0.190000D0,0.290000D0,0.0D0,0.270000D0,0.260000D0,
     &      0.273000D0,0.150000D0,0.0D0,0.0D0/
      DATA (QTAB(I),I=41,80)/ - 0.280000D0,0.200000D0,0.340000D0,
     &      0.440000D0,0.0D0,0.660000D0,0.0D0,0.0D0,0.861000D0,0.0D0,
     &      - 0.330000D0,0.0D0, - 0.789000D0,0.0D0, - 0.003000D0,
     &      0.340000D0,0.200000D0,0.0D0, - 0.041000D0, - 0.560000D0,
     &      0.660000D0, - 0.180000D0,3.920000D0,1.340000D0,1.340000D0,
     &      2.510000D0,2.730000D0,2.827000D0,0.0D0,2.800000D0,
     &      5.680000D0,5.100000D0,3.440000D0,0.0D0,2.220000D0,
     &      0.800000D0,0.700000D0,0.0D0,0.594000D0,0.0D0/
      DATA (QTAB(I),I=81,110)/0.0D0,0.0D0, - 0.460000D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,1.700000D0,4.400000D0,0.0D0,4.300000D0,
     &      0.0D0,0.0D0,4.900000D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
C
      IF ( Z.GE.0 .AND. Z.LE.110 ) THEN
         A = ATAB(Z)
         MUNUC = MUTAB(Z)
         INUC = ITAB(Z)
         QNUC = QTAB(Z)
         IF ( IPRINT.GT.0 ) THEN
            CHSYM = TAB_CHSYM(Z)
            WRITE (6,99001) MUNUC,A,CHSYM
            WRITE (6,99002) INUC,QNUC
         END IF
         IF ( A.NE.0 ) RETURN
         IF ( IPRINT.GT.0 ) WRITE (6,99003)
      END IF
      WRITE (6,*) ' TABMUQIN: warning for atomic number = ',Z
99001 FORMAT (' nuclear magnetic moment mu/mu_N =',F8.5,'  for ',I3,A2)
99002 FORMAT (' nuclear spin I =',F5.1,
     &        ' nuclear quadrupole moment Q = ',F9.6,' [barn]')
99003 FORMAT (' no stable isotope with non-vanishing mu')
      END
C*==tabgamman.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TABGAMMAN(Z,GAMN)
C   ********************************************************************
C   *                                                                  *
C   *     get nuclear gyro-magnetic ratio  GAMN=gamma_n in  [1/s*G]    *
C   *                                                                  *
C   *     GN = gamma_n/2*pi [kHz/G] taken from Carter et al.           *
C   *          Prog. Mater. Sci. Sci. I, 1977, Tab. 9.2a               *
C   *          the isotope listed has mass number  MN  and is          *
C   *          normally that with the highest natural abundance        *
C   *                                                                  *
C   *                                                         03/08/93 *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_TABLES,ONLY:TAB_CHSYM
      IMPLICIT NONE
C*--TABGAMMAN145
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 GAMN
      INTEGER Z
C
C Local variables
C
      CHARACTER*2 CHSYM
      REAL*8 GN(0:110)
      INTEGER I,MN(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA MN/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &     63,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,103,105,109,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,195,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
C
      DATA (GN(I),I=1,40)/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,1.12850D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0/
      DATA (GN(I),I=41,80)/0.0D0,0.0D0,0.0D0,0.0D0,0.13400D0,0.19500D0,
     &      0.19808D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.90940D0,0.0D0,0.0D0/
      DATA (GN(I),I=81,110)/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0/
C
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         GAMN = GN(Z)*2*PI*1000
         CHSYM = TAB_CHSYM(Z)
         WRITE (6,99001) GN(Z),MN(Z),CHSYM
         IF ( ABS(GAMN).GT.1D-4 ) RETURN
         WRITE (6,*) ' table uncomplete '
      END IF
      WRITE (6,*) ' atomic number = ',Z
      STOP ' in <TABGAMMAN>'
99001 FORMAT (' nuclear gyro-magnetic ratio  g_n/2*pi=',F8.5,
     &        ' [kHz/G]    for ',I3,A2)
      END
C*==tabrwshls.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABRWSHLS(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return the Wigner-Seitz-radius RWS for atomic number  Z      *
C   *     experimental values taken from H.L.Skrivers's book           *
C   *                                                                  *
C   *                                                         03/08/93 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABRWSHLS215
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      REAL*8 TABRWSHLS
C
C Local variables
C
      INTEGER I
      REAL*8 TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA (TABLE(I),I=0,40)/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,4.862D0,4.122D0,3.427D0,3.052D0,2.818D0,2.684D0,
     &      2.699D0,2.662D0,2.621D0,2.602D0,2.669D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,0.0D0,0.0D0,0.0D0,5.197D0,4.494D0,3.761D0,3.347D0/
      DATA (TABLE(I),I=41,80)/3.071D0,2.922D0,2.840D0,2.791D0,2.809D0,
     &      2.873D0,3.005D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      5.656D0,4.652D0,3.920D0,3.800D0,3.818D0,3.804D0,3.783D0,
     &      3.768D0,4.263D0,3.764D0,3.720D0,3.704D0,3.687D0,3.668D0,
     &      3.649D0,4.052D0,3.624D0,3.301D0,3.069D0,2.945D0,2.872D0,
     &      2.825D0,2.835D0,2.897D0,3.002D0,0.0D0/
      DATA (TABLE(I),I=81,110)/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      5.900D0,4.790D0,3.900D0,3.756D0,3.430D0,3.221D0,3.140D0,
     &      3.181D0,3.614D0,3.641D0,3.550D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &      0.0D0,3.500D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
C
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABRWSHLS = TABLE(Z)
         IF ( ABS(TABRWSHLS).GT.1D-5 ) RETURN
         WRITE (6,*) ' table uncomplete '
      END IF
      WRITE (6,*) ' atomic number = ',Z
      STOP ' in <TABRWSHLS>'
      END
C*==tabestmom.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABESTMOM(Z)
C   ********************************************************************
C   *                                                                  *
C   *  return the estimated spin magnetic moment for atomic number  Z  *
C   *                                                                  *
C   *                                                         26/08/93 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABESTMOM275
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      REAL*8 TABESTMOM
C
C Local variables
C
      REAL*8 TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA TABLE/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.6D0,4.0D0,2.1D0,1.6D0,0.6D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,
     &     2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,7.0D0,6.0D0,5.0D0,4.0D0,3.0D0,
     &     2.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0/
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABESTMOM = TABLE(Z)
         RETURN
      END IF
      WRITE (6,*) ' atomic number = ',Z
      STOP ' in <TABESTMOM>'
      END
C*==tab_curie_temp.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TAB_CURIE_TEMP(Z)
C   ********************************************************************
C   *                                                                  *
C   *  return the experimental Curie temperature in [K] for given Z    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--TAB_CURIE_TEMP328
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      REAL*8 TAB_CURIE_TEMP
C
C Local variables
C
      REAL*8 TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA TABLE/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.00D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.000D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.043D3,
     &     1.4D3,6.27D2,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,2.92D2,0.0D0,8.8D1,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &     0.0D0,0.0D0,0.0D0,0.0D0/
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TAB_CURIE_TEMP = TABLE(Z)
         RETURN
      END IF
      WRITE (6,*) ' atomic number = ',Z
      STOP ' in <TABESTMOM>'
      END
C*==tab_est_w_weiss.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TAB_EST_W_WEISS(Z)
C   ********************************************************************
C   *                                                                  *
C   *  return the estimated Weiss field for atomic number  Z           *
C   *  for all magnetic elements  -  otherwise:  w = 0.0               *
C   *  based on the experimental Curie temperature T_C                 *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:KB_SI,EV_J,RY_EV
      IMPLICIT NONE
C*--TAB_EST_W_WEISS387
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      REAL*8 TAB_EST_W_WEISS
C
C Local variables
C
      REAL*8 TAB_CURIE_TEMP
      REAL*8 TC
C
C*** End of declarations rewritten by SPAG
C
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
C
         TC = TAB_CURIE_TEMP(Z)
C ---------------------------------------------- Weiss field in Ry units
         TAB_EST_W_WEISS = 3.D0*TC*KB_SI/EV_J/RY_EV
C
         RETURN
      END IF
      WRITE (6,*) ' atomic number = ',Z
      STOP ' in <TABESTMOM>'
      END
C*==tabncore.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABNCORE(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return number of core electrons for atomic number  Z         *
C   *                                                                  *
C   *                                                         26/03/04 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABNCORE434
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      INTEGER TABNCORE
C
C Local variables
C
      INTEGER I,TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA (TABLE(I),I=0,10)/0,0,0,2,2,2,2,2,2,4,4/
      DATA (TABLE(I),I=11,18)/10,10,10,10,10,10,12,12/ 
      DATA (TABLE(I),I=19,36)/18,18,18,18,18,18,18,18,18,18,18,18,18,28,
     &      28,28,30,30/ 
      DATA (TABLE(I),I=37,54)/36,36,36,36,36,36,36,36,36,36,36,36,46,46,
     &      46,46,48,48/
      DATA (TABLE(I),I=55,86)/54,54,54,54,54,54,54,54,54,54,54,54,54,54,
     &      54,54,54,68,68,68,68,68,68,68,68,68,78,78,78,78,80,80/
      DATA (TABLE(I),I=87,110)/86,86,86,86,86,86,86,86,86,86,86,86,86,
     &      86,86,86,86,86,86,86,86,86,86,86/
C
      If ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABNCORE = TABLE(Z)
      ELSE
         WRITE (6,*) ' atomic number = ',Z
         STOP ' in <TABNCORE>'
      END IF
      END
C*==tabnval.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABNVAL(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return number of valence electrons for atomic number  Z      *
C   *                                                                  *
C   *                                                         03/08/93 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABNVAL487
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      INTEGER TABNVAL
C
C Local variables
C
      INTEGER TABNCORE
C
C*** End of declarations rewritten by SPAG
C
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABNVAL = Z - TABNCORE(Z)
      ELSE
         WRITE (6,*) ' atomic number = ',Z
         STOP ' in <TABNVAL>'
      END IF
      END
C*==tabnlval.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABNLVAL(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return the minimum number of (l_max+1) needed to represent   *
C   *     the valence states for atomic number  Z                      *
C   *                                                                  *
C   *     f-electrons are assumed to be core electrons                 *
C   *                                                         26/08/93 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABNLVAL531
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      INTEGER TABNLVAL
C
C Local variables
C
      INTEGER TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA TABLE/1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
     &     3,3,3,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,
     &     3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,
     &     2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABNLVAL = TABLE(Z)
      ELSE
         WRITE (6,*) ' atomic number = ',Z
         STOP ' in <TABNLVAL>'
      END IF
      END
C*==tabnsemcor.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION TABNSEMCOR(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return number of semi-core electrons for atomic number  Z    *
C   *                                                                  *
C   *                                                         03/08/93 *
C   ********************************************************************
      IMPLICIT NONE
C*--TABNSEMCOR577
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      REAL*8 TABNSEMCOR
C
C Local variables
C
      INTEGER TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA TABLE/0,0,0,0,0,0,0,0,0,0,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     &     6,6,6,10,10,10,10,10,10,10,6,6,6,6,6,6,6,6,6,6,6,10,10,10,10,
     &     10,10,10,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     &     10,10,10,10,10,10,10,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     &     6,6,6,6/
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABNSEMCOR = TABLE(Z)
      ELSE
         WRITE (6,*) ' atomic number = ',Z
         STOP ' in <TABNSEMCOR>'
      END IF
      END
C*==atomnumb.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      FUNCTION ATOMNUMB(CHSYM)
C   ********************************************************************
C   *                                                                  *
C   *     return atomic number Z for chemical symbol CHSYM             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TABLES,ONLY:TAB_CHSYM
      IMPLICIT NONE
C*--ATOMNUMB624
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*2 CHSYM
      INTEGER ATOMNUMB
C
C Local variables
C
      CHARACTER*2 CHSYMP
      INTEGER Z
C
C*** End of declarations rewritten by SPAG
C
C
      CHSYMP = CHSYM
C
      CALL STRING_CONVERT_TO_LC(CHSYMP(2:2))
C
      IF ( CHSYM(2:2).EQ.' ' .OR. CHSYM(2:2).EQ.'_' .OR. CHSYM(2:2)
     &     .EQ.'.' ) CHSYMP = CHSYM(1:1)//' '
C
      ATOMNUMB = 0
      DO Z = 0,104
         IF ( CHSYMP.EQ.TAB_CHSYM(Z) ) THEN
            ATOMNUMB = Z
            RETURN
         END IF
      END DO
      WRITE (6,*) ' atomic number not found for  CHSYM=',CHSYMP
      STOP
      END
C*==pot_nuc.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      REAL*8 FUNCTION POT_NUC(R,Z)
C   ********************************************************************
C   *                                                                  *
C   *  this function returns the potential inside the nucleus          *
C   *  assuming a homogenously charged sphere                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--POT_NUC678
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 R,Z
C
C Local variables
C
      REAL*8 RNUCTAB
      REAL*8 R_NUC
C
C*** End of declarations rewritten by SPAG
C
      R_NUC = RNUCTAB(NINT(Z))
C
      POT_NUC = (Z/R_NUC**3)*(R**2-3*R_NUC**2)
C
      END
C*==rnuctab.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
C
      REAL*8 FUNCTION RNUCTAB(Z)
C   ********************************************************************
C   *                                                                  *
C   *  this function contains the nuclear radii for many    Z          *
C   *  calculated  from  1.128*a**(1/3)fm   (here saved in a.u.)       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--RNUCTAB719
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
C
C Local variables
C
      REAL*8 RNUC(0:109)
C
C*** End of declarations rewritten by SPAG
C
      DATA RNUC/0.0D0,2.688280D-05,3.394993D-05,3.893984D-05,
     &     4.293464D-05,4.632524D-05,4.930297D-05,5.197748D-05,
     &     5.441818D-05,5.667210D-05,5.877274D-05,6.074494D-05,
     &     6.260766D-05,6.437575D-05,6.606108D-05,6.767330D-05,
     &     6.922041D-05,7.070905D-05,7.214488D-05,7.353271D-05,
     &     7.487670D-05,7.618046D-05,7.744715D-05,7.867954D-05,
     &     7.988009D-05,8.105101D-05,8.219423D-05,8.331154D-05,
     &     8.440452D-05,8.547461D-05,8.652311D-05,8.755123D-05,
     &     8.856005D-05,8.955057D-05,9.052372D-05,9.148035D-05,
     &     9.242124D-05,9.334712D-05,9.425867D-05,9.515651D-05,
     &     9.604123D-05,9.691338D-05,9.777347D-05,9.862198D-05,
     &     9.945937D-05,1.002860D-04,1.011024D-04,1.019088D-04,
     &     1.027057D-04,1.034933D-04,1.042720D-04,1.050420D-04,
     &     1.058038D-04,1.065574D-04,1.073032D-04,1.080414D-04,
     &     1.087723D-04,1.094961D-04,1.102130D-04,1.109231D-04,
     &     1.116267D-04,1.123240D-04,1.130152D-04,1.137003D-04,
     &     1.143796D-04,1.150532D-04,1.157213D-04,1.163840D-04,
     &     1.170414D-04,1.176938D-04,1.183411D-04,1.189835D-04,
     &     1.196212D-04,1.202542D-04,1.208827D-04,1.215068D-04,
     &     1.221265D-04,1.227420D-04,1.233534D-04,1.239607D-04,
     &     1.245640D-04,1.251634D-04,1.257591D-04,1.263510D-04,
     &     1.269393D-04,1.275240D-04,1.281052D-04,1.286830D-04,
     &     1.292574D-04,1.298285D-04,1.303964D-04,1.309611D-04,
     &     1.315227D-04,1.320813D-04,1.326369D-04,1.331895D-04,
     &     1.337392D-04,1.342861D-04,1.348302D-04,1.353715D-04,
     &     1.359102D-04,1.364462D-04,1.369796D-04,1.375105D-04,
     &     1.380389D-04,1.385648D-04,1.390882D-04,1.396093D-04,
     &     1.401280D-04,1.406444D-04/
      RNUCTAB = RNUC(Z)
      END
