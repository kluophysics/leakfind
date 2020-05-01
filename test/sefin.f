C*==sefin.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFIN(IECURR,ERYD,EFERMI,NT,RHOCHR,RHOSPN,KSEFIN,SEBT,
     &                 SEVT,JRWS,IMT,NRMAX,NMMAX,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   set up the spherical complex self energy   SEVT,  SEBT         *
C   *   for  final states accoding to the parameter KSEFIN             *
C   *   using the FEFF-routines of J.J. Rehr  et al.                   *
C   *                                                                  *
C   ********************************************************************
C   * modifications:  - the shift of the potential by EREF is          *
C   *                   suppressed                                     *
C   *                 - all code related to EREF, RSO, GAMACH removed  *
C   *                 - a spin dependent self energy is created        *
C   *                   that is NOT added to the potential             *
C   *                                                                  *
C   *   This and all called subroutines use atomic (HARTREE) units     *
C   *   All input  to and output  from   SEFIN      is in Ry.          *
C   *                                                                  *
C   ********************************************************************
C   *   INPUT                                                          *
C   *   IT, IE      used only for debug and labels.                    *
C   *   NRMAX       number of points in current r-grid                 *
C   *                                                                  *
C   *   KSEFIN       1  Hedin-Lunqvist + const real & imag part        *
C   *               2  Dirac-Hara + const real & imag part             *
C   *               3  Dirac-Hara + HL imag part                       *
C   *                             + const real & imag part             *
C   *   IFIRST      first entry flag, set to zero before first call    *
C   *               for each unique potential                          *
C   *               see VXCEFR and VXCEFI below                        *
C   *   JTOP        index of first interstitial point in current       *
C   *               r grid                                             *
C   *   EM          current energy grid point                          *
C   *               for E<= EFERMI   SEVT = SEBT = 0                   *
C   *   XMU         fermi level                                        *
C   *   VI0         const imag part to subtract from potential         *
C   *   RHO(nrmax)  electron density                                   *
C   *                                                                  *
C   *   OUTPUT                                                         *
C   *   SEVT, SEBT  spin averaged and spin dependent part of           *
C   *               energy dependent self energy                       *
C   *                                                                  *
C   *   WORKSPACE                                                      *
C   *   VXCEFR and VXCEFI are calculated only on first entry for a     *
C   *   particular unique potential, re-used on subsequent entries.    *
C   *   VXCEFR(nrmax)  real part of xc at fermi level                  *
C   *   VXCEFI(nrmax)  imag part of xc at fermi level                  *
C   *                                                                  *
C   *   kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37      *
C   *                                                                  *
C   *   fine structure alpha                                           *
C   *   speed of light in louck's units (rydbergs?)                    *
C   *                                                                  *
C   *   First calculate vxc to correct the local momentum dispersion   *
C   *   relation, delta = vxc(e,k) - vxc(mu,k), and                    *
C   *             p^2 = k^2 -mu + kf^2 - delta.                        *
C   *                                                                  *
C   *   In jr theory, V(E,r) = Vcoul(r) + Vxc(e,r)                     *
C   *                        = Vcoul(r) + Vxcgs(r) + delta(e,r).       *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI,CONST_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION THIRD,FA
      PARAMETER (THIRD=1D0/3D0,FA=1.919158292677512811D0)
C
C Dummy arguments
C
      REAL*8 EFERMI
      COMPLEX*16 ERYD
      INTEGER IECURR,KSEFIN,NMMAX,NRMAX,NT,NTMAX
      INTEGER IMT(NTMAX),JRWS(NMMAX)
      REAL*8 RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX)
      COMPLEX*16 SEBT(NRMAX,NTMAX),SEVT(NRMAX,NTMAX)
C
C Local variables
C
      DOUBLE COMPLEX DELTA,SIGMA(NRMAX,2)
      REAL*8 EM,RHO(:),VXCEFI(:,:),VXCEFR(:,:),VZ
      INTEGER I,IA_ERR,ICUSP,IFIRST,IM,IR,IS,IT,JTOP
      DOUBLE PRECISION RS,VI0,VXCI,VXCR,XF,XK,XK2,XMU
      SAVE VXCEFI,VXCEFR
      EXTERNAL SEFINEDP,SEFINIMHL,SEFINRHL
C
C*** End of declarations rewritten by SPAG
C
      DATA IFIRST/0/
C
      ALLOCATABLE VXCEFI,VXCEFR,RHO
C
      ALLOCATE (RHO(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SEFIN -> NRMAX'
C
      IFIRST = IFIRST + 1
      IF ( IFIRST.EQ.1 ) THEN
         ALLOCATE (VXCEFI(NRMAX,NT),VXCEFR(NRMAX,NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:SEFIN -> VXCEFR'
      END IF
C
      IF ( KSEFIN.LE.0 ) RETURN
      IF ( KSEFIN.GT.3 ) STOP 'in <SEFIN> ......  KSEFIN > 3'
C
C***********************************************************************
C                     interface to Rehr's subroutines
C***********************************************************************
C
      EM = DREAL(ERYD)
      XMU = EFERMI
      VI0 = 0D0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IM = IMT(IT)
         JTOP = JRWS(IM)
C
C------------------------------------------ spin loop: IS=1: down, 2: up
         DO IS = 1,2
            VZ = 2D0*DBLE(IS) - 3.D0
C
            DO IR = 1,JTOP
               RHO(IR) = RHOCHR(IR,IT) - VZ*RHOSPN(IR,IT)
               RHO(IR) = RHO(IR)/CONST_4PI
            END DO
C
C***********************************************************************
C                                                                      *
C        Add the self energy correction
C
            DO I = 1,JTOP
C
               RS = (3/(4*PI*RHO(I)))**THIRD
C           xf = 1.9191.../rs
               XF = FA/RS
C
C           vxc_mu indep of energy, calc only once
C           Calculate vxc at fermi level e = mu, j.m. 1/12/89
C
               IF ( IFIRST.EQ.1 ) THEN
C
                  XK = XF*1.00001D0
                  IF ( KSEFIN.EQ.1 )
     &                 CALL SEFINRHL(RS,XK,VXCEFR(I,IT),VXCEFI(I,IT))
                  IF ( KSEFIN.EQ.2 )
     &                 CALL SEFINEDP(RS,XK,VI0,VXCEFR(I,IT),VXCEFI(I,IT)
     &                 )
                  IF ( KSEFIN.EQ.3 ) THEN
                     CALL SEFINEDP(RS,XK,VI0,VXCEFR(I,IT),VXCEFI(I,IT))
                     CALL SEFINIMHL(RS,XK,VXCEFI(I,IT),ICUSP)
                  END IF
C
               END IF
C
C           xk2 is the local momentum squared, p^2 = k^2 - mu + kf^2,
C           k^2 represents energy measured from vacuum.
C           See formula 2.15 in Lee and Beni's paper with the last 2
C           terms neglected.  (complete reference?)
C
               XK2 = EM + XF**2 - XMU
C
               IF ( XK2.LT.0 ) THEN
                  DELTA = 0D0
                  GOTO 10
C                  PRINT *,'i, jtop',I,JTOP
C                  PRINT *,'rs, densty(i)',RS,RHO(I)
C                  PRINT *,'xf, fa',XF,FA
C                  PRINT *,'em, xmu, xk2',EM,XMU,XK2
C                  STOP 'SEFINXCPOT-1'
               END IF
C
               XK = SQRT(XK2)
               IF ( KSEFIN.EQ.1 ) CALL SEFINRHL(RS,XK,VXCR,VXCI)
               IF ( KSEFIN.EQ.2 ) CALL SEFINEDP(RS,XK,VI0,VXCR,VXCI)
               IF ( KSEFIN.EQ.3 ) THEN
                  CALL SEFINEDP(RS,XK,VI0,VXCR,VXCI)
                  CALL SEFINIMHL(RS,XK,VXCI,ICUSP)
               END IF
C
               DELTA = DCMPLX(VXCR-VXCEFR(I,IT),VXCI-VXCEFI(I,IT))
C
C           Correct local momentum according to the formula
C           p^2 = k^2 - mu + kf^2 - delta.  Note that imag part
C           of delta is ignored, since xk2 is a real quantity.
C
Ccccccccccccccccc               XK2 = EM + XF**2 - XMU - DREAL(DELTA)
               IF ( XK2.LT.0 ) THEN
                  PRINT *,' XK2, I, IECURR, IT, NT',XK2,I,IECURR,IT,NT
                  PRINT *,'EM, XF**2, XMU, DELTA'
                  PRINT *,EM,XF**2,XMU,DELTA
                  STOP 'SEFINXCPOT-2'
C
               END IF
C
               XK = SQRT(XK2)
C
C           recalculate vxc(e,k) and vxc(mu,k) with the corrected
C           local momentum
C
               IF ( KSEFIN.EQ.1 ) CALL SEFINRHL(RS,XK,VXCR,VXCI)
               IF ( KSEFIN.EQ.2 ) CALL SEFINEDP(RS,XK,VI0,VXCR,VXCI)
               IF ( KSEFIN.EQ.3 ) THEN
                  CALL SEFINEDP(RS,XK,VI0,VXCR,VXCI)
                  CALL SEFINIMHL(RS,XK,VXCI,ICUSP)
               END IF
C
C           delta corrected calculated with new local momentum
C
               DELTA = DCMPLX(VXCR-VXCEFR(I,IT),VXCI-VXCEFI(I,IT))
C
C           Note multiplication by 2 in the exchange correlation part to
C           to convert it to RYDBERG units.
C
 10            CONTINUE
               SIGMA(I,IS) = 2*DELTA
C
            END DO
C                                                                      *
C***********************************************************************
C                     interface to Rehr's subroutines
C***********************************************************************
         END DO
C
C------------------------------------------ spin loop: IS=1: down, 2: up
         DO IR = 1,JTOP
            SEVT(IR,IT) = (SIGMA(IR,2)+SIGMA(IR,1))/2D0
            SEBT(IR,IT) = (SIGMA(IR,2)-SIGMA(IR,1))/2D0
         END DO
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
C*==sefinrhl.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFINRHL(RS,XK,ERL,EIM)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *   input:  rs, xk                                                 *
C   *   output: erl, eim                                               *
C   *                                                                  *
C   *   This is a new hl subroutine, using interpolation for the       *
C   *   real part while the imaginary part is calculated analytically. *
C   *   It uses hl to calculate values at the mesh points for the      *
C   *   inter-polation of the real part. The imaginary part is         *
C   *   calculated using subroutine imhl.                              *
C   *                                                                  *
C   *   written by jose mustre                                         *
C   *   polynomial in rs has a 3/2 power term. j.m.                    *
C   *                                                                  *
C   *                                                                  *
C   *   for the right branch the interpolation has the form:           *
C   *   hl(rs,x) = e/x + f/x**2 + g/x**3                               *
C   *   where e is known and                                           *
C   *      f = sum (i=1,3) ff(i) rs**(i+1)/2                           *
C   *      g = sum (i=1,3) gg(i) rs**(i+1)/2                           *
C   *                                                                  *
C   *                                                                  *
C   *   lrs=number of rs panels, in this case one has 4 panels         *
C   *   nrs=number of standard rs values, also order of rs expansion   *
C   *   if you change nrs you need to change the expansion of hl       *
C   *   in powers of rs that only has 3 terms!                         *
C   *   nleft=number of coefficients for x<x0                          *
C   *   nright=number of coefficients for x>x0                         *
C   *                                                                  *
C   *   kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37      *
C   *                                                                  *
C   *   fine structure alpha                                           *
C   *   speed of light in louck's units (rydbergs?)                    *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LRS,NRS,NLEFT,NRIGHT
      PARAMETER (LRS=4,NRS=3,NLEFT=4,NRIGHT=2)
      DOUBLE PRECISION FA
      PARAMETER (FA=1.919158292677512811D0)
C
C Dummy arguments
C
      DOUBLE PRECISION EIM,ERL,RS,XK
C
C Local variables
C
      DOUBLE PRECISION CLEFT(NLEFT),CRIGHT(NRIGHT),EEE,EF,
     &                 RCFL(LRS,NRS,NLEFT),RCFR(LRS,NRS,NRIGHT),RKF,WP,
     &                 XX
      INTEGER ICUSP,J,MRS
      EXTERNAL SEFINIMHL
C
C*** End of declarations rewritten by SPAG
C
      DATA RCFR/ - 0.173963D+00, - 0.173678D+00, - 0.142040D+00,
     &     - 0.101030D+00, - 0.838843D-01, - 0.807046D-01,
     &     - 0.135577D+00, - 0.177556D+00, - 0.645803D-01,
     &     - 0.731172D-01, - 0.498823D-01, - 0.393108D-01,
     &     - 0.116431D+00, - 0.909300D-01, - 0.886979D-01,
     &     - 0.702319D-01,0.791051D-01, - 0.359401D-01, - 0.379584D-01,
     &     - 0.419807D-01, - 0.628162D-01,0.669257D-01,0.667119D-01,
     &     0.648175D-01/
      DATA RCFL/0.590195D+02,0.478860D+01,0.812813D+00,0.191145D+00,
     &     - 0.291180D+03, - 0.926539D+01, - 0.858348D+00,
     &     - 0.246947D+00,0.363830D+03,0.460433D+01,0.173067D+00,
     &     0.239738D-01, - 0.181726D+03, - 0.169709D+02, - 0.409425D+01,
     &     - 0.173077D+01,0.886023D+03,0.301808D+02,0.305836D+01,
     &     0.743167D+00, - 0.110486D+04, - 0.149086D+02, - 0.662794D+00,
     &     - 0.100106D+00,0.184417D+03,0.180204D+02,0.450425D+01,
     &     0.184349D+01, - 0.895807D+03, - 0.318696D+02, - 0.345827D+01,
     &     - 0.855367D+00,0.111549D+04,0.156448D+02,0.749582D+00,
     &     0.117680D+00, - 0.620411D+02, - 0.616427D+01, - 0.153874D+01,
     &     - 0.609114D+00,0.300946D+03,0.109158D+02,0.120028D+01,
     &     0.290985D+00, - 0.374494D+03, - 0.535127D+01, - 0.261260D+00,
     &     - 0.405337D-01/
C
C     calculate hl using interpolation coefficients
C
      RKF = FA/RS
      EF = RKF**2/2
      WP = SQRT(3/RS**3)
      CALL SEFINIMHL(RS,XK,EIM,ICUSP)
C
C     eim already has a factor of ef in it j.m.
C     eim also gives the position of the cusp
C
      XX = XK/RKF
C     set to fermi level if below fermi level
      IF ( XX.LT.1.00001D0 ) XX = 1.00001D0
C     calculate right hand side coefficients
      IF ( RS.LT.0.2D0 ) THEN
         MRS = 1
C
      ELSE IF ( RS.LT.1.0D0 ) THEN
         MRS = 2
C
      ELSE IF ( RS.LT.5.0D0 ) THEN
         MRS = 3
C
      ELSE
         MRS = 4
      END IF
C
      DO J = 1,NRIGHT
         CRIGHT(J) = RCFR(MRS,1,J)*RS + RCFR(MRS,2,J)*RS*SQRT(RS)
     &               + RCFR(MRS,3,J)*RS**2
      END DO
      EEE = -PI*WP/(4*RKF*EF)
C
      IF ( ICUSP.NE.1 ) THEN
         DO J = 1,NLEFT
            CLEFT(J) = RCFL(MRS,1,J)*RS + RCFL(MRS,2,J)*RS**1.5D0 + 
     &                 RCFL(MRS,3,J)*RS**2
         END DO
         ERL = CLEFT(1)
         DO J = 2,NLEFT
            ERL = ERL + CLEFT(J)*XX**(J-1)
         END DO
C
      ELSE
C        right branch
         ERL = EEE/XX
         DO J = 1,NRIGHT
            ERL = ERL + CRIGHT(J)/XX**(J+1)
         END DO
      END IF
C
      ERL = ERL*EF
C
      END
C*==sefinimhl.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFINIMHL(RS,XK,EIM,ICUSP)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *      what is xk?  k**2 - mu + kf**2?                             *
C   *                                                                  *
C   *  written by j. mustre (march 1988)                               *
C   *  code is based on analytical expression derived by john rehr.    *
C   *  it leaves the real part, calculated in rhl unchanged.           *
C   *                                                                  *
C   *  modified by j. rehr  (oct 1991) - adds quinn approximation      *
C   *  for losses due to electron-hole pairs below the plasmon turn    *
C   *  on see new subroutine quinn.f, which incorporates r. albers     *
C   *  coding of j.j. quinn's approximations for details.              *
C   *                                                                  *
C   *  kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37       *
C   *                                                                  *
C   *  fine structure alpha                                            *
C   *  speed of light in louck's units (rydbergs?)                     *
C   *  alph is Hedin-Lundquist parameter                               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION FA,ALPH
      PARAMETER (FA=1.919158292677512811D0,ALPH=4.0D0/3.0D0)
C
C Dummy arguments
C
      DOUBLE PRECISION EIM,RS,XK
      INTEGER ICUSP
C
C Local variables
C
      DOUBLE PRECISION D1,D2,D3,EF,EI,Q2,QMINUS,QPLUS,QU,RAD,WP,XF,XK0,
     &                 XS
      INTEGER ICOUNT
      DOUBLE PRECISION SEFINFFQ
      EXTERNAL SEFINCUBIC,SEFINFFQ,SEFINQUINN
C
C*** End of declarations rewritten by SPAG
C
      DATA ICOUNT/0/
C
      ICUSP = 0
      XF = FA/RS
      EF = XF**2/2
C
C     xk0 is xk normalized by k fermi.
      XK0 = XK/XF
C     set to fermi level if below fermi level
      IF ( XK0.LT.1.00001D0 ) XK0 = 1.00001D0
C
C     wp is given in units of the fermi energy in the formula below.
      WP = SQRT(3/RS**3)/EF
      XS = WP**2 - (XK0**2-1)**2
C
      EIM = 0
      IF ( XS.LT.0.D0 ) THEN
         Q2 = SQRT((SQRT(ALPH**2-4*XS)-ALPH)/2)
         QU = MIN(Q2,(1+XK0))
         D1 = QU - (XK0-1)
         IF ( D1.GT.0 ) EIM = SEFINFFQ(QU,EF,XK,WP,ALPH)
     &                        - SEFINFFQ(XK0-1,EF,XK,WP,ALPH)
C
      END IF
C
      CALL SEFINCUBIC(XK0,WP,ALPH,RAD,QPLUS,QMINUS)
C
      IF ( RAD.LE.0 ) THEN
         D2 = QPLUS - (XK0+1)
         IF ( D2.GT.0 ) EIM = EIM + SEFINFFQ(QPLUS,EF,XK,WP,ALPH)
     &                        - SEFINFFQ(XK0+1,EF,XK,WP,ALPH)
C
         D3 = (XK0-1) - QMINUS
         IF ( D3.GT.0 ) THEN
            EIM = EIM + SEFINFFQ(XK0-1,EF,XK,WP,ALPH)
     &            - SEFINFFQ(QMINUS,EF,XK,WP,ALPH)
C           beginning of the imaginary part and position of the cusp x0
            ICUSP = 1
         END IF
C
      END IF
C
      CALL SEFINQUINN(XK0,RS,WP,EF,EI)
      IF ( EIM.GE.EI ) EIM = EI
C
      ICOUNT = ICOUNT + 1
C
      END
C*==sefinedp.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFINEDP(RS,XK,VI0,VR,VI)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION FA
      PARAMETER (FA=1.919158292677512811D0)
C
C Dummy arguments
C
      DOUBLE PRECISION RS,VI,VI0,VR,XK
C
C Local variables
C
      DOUBLE PRECISION C,X,XF
C
C*** End of declarations rewritten by SPAG
C
C     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
C
C     fine structure alpha
C     speed of light in louck's units (rydbergs?)
C
      XF = FA/RS
C
C     p = sqrt (k^2 + kf^2) is the local momentum, and x = p / kf
C     Reference formula 23 in Role of Inelastic effects in EXAFS
C     by Rehr and Chou. EXAFS1 conference editted by Bianconi.
C     x is local momentum in units of fermi momentum
C
      X = XK/XF
      X = X + 1.0D-5
C     set to fermi level if below fermi level
      IF ( X.LT.1.00001D0 ) X = 1.00001D0
      C = ABS((1+X)/(1-X))
      C = LOG(C)
      VR = -(XF/PI)*(1+C*(1-X**2)/(2*X))
C
C     Note vi=vi0/2 to have both real and imaginary part in hartrees
C     to be consistent with  other subroutines.
      VI = VI0/2
C
      END
C*==sefinffq.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      DOUBLE PRECISION FUNCTION SEFINFFQ(Q,EF,XK,WP,ALPH)
C   ********************************************************************
C   *                                                                  *
C   *    input:  q, wp, alph, ef, xk                                   *
C   *             q is dimensionless, normalized to fermi momentum     *
C   *            xk is momentum in invBohrs                            *
C   *    output: ffq only                                              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION ALPH,EF,Q,WP,XK
C
C Local variables
C
      DOUBLE PRECISION WQ
C
C*** End of declarations rewritten by SPAG
C
      WQ = SQRT(WP**2+ALPH*Q**2+Q**4)
      SEFINFFQ = (WP+WQ)/(Q**2) + ALPH/(2*WP)
      SEFINFFQ = ((EF*WP)/(4*XK))*LOG(SEFINFFQ)
C
      END
C*==sefincubic.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFINCUBIC(XK0,WP,ALPH,RAD,QPLUS,QMINUS)
C   ********************************************************************
C   *                                                                  *
C   *     input:  xk0, wp, alph                                        *
C   *     output: rad, qplus, qminus                                   *
C   *     this subroutine finds the roots of the equation              *
C   *     4xk0 * q^3  +  (alph-4xk0^2) * q^2  +  wp^2 = 0              *
C   *     see abramowitz and stegun pg 17 for formulae.                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION THREE,THIRD
      PARAMETER (THREE=3,THIRD=1D0/3D0)
C
C Dummy arguments
C
      DOUBLE PRECISION ALPH,QMINUS,QPLUS,RAD,WP,XK0
C
C Local variables
C
      DOUBLE PRECISION A0,A1,A2,Q,R
      DOUBLE COMPLEX S1,S13
C
C*** End of declarations rewritten by SPAG
C
      A2 = (ALPH/(4*XK0**2)-1)*XK0
      A0 = WP**2/(4*XK0)
      A1 = 0
      Q = A1/3 - A2**2/9
      R = (A1*A2-3*A0)/6 - A2**3/27
      RAD = Q**3 + R**2
      IF ( RAD.GT.0 ) THEN
         QPLUS = 0
         QMINUS = 0
         RETURN
C
      END IF
C
      S13 = DCMPLX(R,SQRT(-RAD))
      S1 = S13**THIRD
      QPLUS = DREAL(2*S1-A2/3)
      QMINUS = DREAL(-(S1-SQRT(THREE)*DIMAG(S1)+A2/3))
C
      END
C*==sefinquinn.f    processed by SPAG 6.70Rc at 21:16 on 19 Dec 2016
      SUBROUTINE SEFINQUINN(X,RS,WP,EF,EI)
C   ********************************************************************
C   *                                                                  *
C   *     input  x, rs, wp, ef                                         *
C   *     output ei                                                    *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *   quinn: calculates low energy gamma (approx. prop. to e**2)     *
C   *   formula taken from john j. quinn, phys. rev. 126,              *
C   *   1453 (1962); equation (7).                                     *
C   *   a cut-off is set up at quinn's cutoff + ef = ekc               *
C   *   it is arounded inverted step function (a fermi function)       *
C   *   theta = 1/( 1 + exp((e-ekc)/gam)) )                            *
C   *   where the rounding factor gam is set to be about 0.3 ekc.      *
C   *   modified by j. rehr (oct 1991) based on coding of r. albers    *
C   *   subroutines quinn.f and quinnc.f                               *
C   *                                                                  *
C   *   variables:                                                     *
C   *      x  = p/pf                                                   *
C   *      rs = ws density parameter                                   *
C   *      ei = imaginary self energy                                  *
C   *      pfqryd = quinn's prefactor in atomic-rydberg units          *
C   *      wkc = quinn's plasmon threshold                             *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *   kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37      *
C   *                                                                  *
C   *   fine structure alpha                                           *
C   *   speed of light in louck's units (rydbergs?)                    *
C   *                                                                  *
C   *   calculate quinn prefactor in atomin Hartree units              *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI,SQRT_PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION FA,ALPHAQ
      PARAMETER (FA=1.919158292677512811D0,ALPHAQ=1/FA)
C
C Dummy arguments
C
      DOUBLE PRECISION EF,EI,RS,WP,X
C
C Local variables
C
      DOUBLE PRECISION ARG,EABS,EKC,F,GAM,PFQ,TEMP1,TEMP2,WKC
C
C*** End of declarations rewritten by SPAG
C
      PFQ = SQRT_PI/(32*(ALPHAQ*RS)**1.5D0)
      TEMP1 = ATAN(SQRT(PI/(ALPHAQ*RS)))
      TEMP2 = SQRT(ALPHAQ*RS/PI)/(1+ALPHAQ*RS/PI)
      PFQ = PFQ*(TEMP1+TEMP2)
C
C     calculate quinn cutoff
C     wkc = quinn's plasmon threshold
C     wkc is cut-off of quinn, pr126, 1453, 1962, eq. (11)
C     in formulae below wp=omegap/ef
      WKC = (SQRT(1+WP)-1)**2
      WKC = (1+(6.D0/5.D0)*WKC/WP**2)*WP*EF
C
C     we add fermi energy to get correct energy for
C     plasma excitations to turn on
      EKC = WKC + EF
C
C     calculate gamma
C     gamryd = 2 * (pfqryd/x) * (x**2-1)**2
      GAM = (PFQ/X)*(X**2-1)**2
C
C     put in fermi function cutoff
      EABS = EF*X**2
      ARG = (EABS-EKC)/(0.3D0*EKC)
      F = 0
      IF ( ARG.LT.80 ) F = 1/(1+EXP(ARG))
C
      EI = -GAM*F/2
C
      END
