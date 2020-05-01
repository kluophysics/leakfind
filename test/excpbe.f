C*==excpbe.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCPBE(V,E,W,IRTOP,NRMAX,RHO,AGRDRHO,LAPRHO,GDGAG,RHOU,
     &                  AGRDRHOU,LAPRHOU,GDGAGU,RHOD,AGRDRHOD,LAPRHOD,
     &                  GDGAGD)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine for   PBE  and  EV-GGA  subroutines              *
C   *                                                                  *
C   *   RHO                 RHO without any additional factors         *
C   *   AGRDRHO             | grad RHO |                               *
C   *   LAPRHO              grad^2 RHO                                 *
C   *   GDGAG               grad RHO . grad | grad RHO |               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SCF,ONLY:SCFVXC
      IMPLICIT NONE
C*--EXCPBE19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EXC_PBE_EV_GGA')
      REAL*8 THRD,THRD2
      PARAMETER (THRD=0.333333333333D0,THRD2=0.666666666667D0)
C
C Dummy arguments
C
      INTEGER IRTOP,NRMAX
      REAL*8 AGRDRHO(NRMAX),AGRDRHOD(NRMAX),AGRDRHOU(NRMAX),E(NRMAX),
     &       GDGAG(NRMAX),GDGAGD(NRMAX),GDGAGU(NRMAX),LAPRHO(NRMAX),
     &       LAPRHOD(NRMAX),LAPRHOU(NRMAX),RHO(NRMAX),RHOD(NRMAX),
     &       RHOU(NRMAX),V(NRMAX,2),W(NRMAX,2)
C
C Local variables
C
      REAL*8 AGRDSP,CONF,CONRS,D,EC,ECRS,ECZET,EX,EXC,FK,G,GDGAGDSP,
     &       LAPSP,RS,SK,SS,TT,U,UU,VCDN,VCUP,VV,VX,VXCDN,VXCUP,WW,XD,
     &       XU,ZET
      INTEGER IR,JSP,LLDA
C
C*** End of declarations rewritten by SPAG
C
      CONF = (3.D0*PI**2)**THRD
      CONRS = (3.D0/(4.D0*PI))**THRD
      LLDA = 0
C
      DO IR = IRTOP + 1,NRMAX
         AGRDRHO(IR) = 0.0D0
         LAPRHO(IR) = 0.0D0
         V(IR,1) = 0.0D0
         V(IR,2) = 0.0D0
         W(IR,1) = 0.0D0
         W(IR,2) = 0.0D0
         E(IR) = 0.0D0
      END DO
C
      DO IR = 1,IRTOP
         EXC = 0.0D0
         VXCUP = 0.0D0
         VXCDN = 0.0D0
         IF ( RHOU(IR).GT.1D-12 .AND. RHOD(IR).GT.1D-12 ) THEN
C
C ------------------------------------ begin the spin loop for exchange
            IF ( RHOU(IR).LE.1D-6 .OR. RHOD(IR).LE.1D-6 ) LLDA = 1
            DO JSP = 1,2
               IF ( JSP.EQ.1 ) THEN
                  D = 2.0D0*RHOU(IR)
                  AGRDSP = 2.0D0*AGRDRHOU(IR)
                  GDGAGDSP = 2.0D0*2.0D0*GDGAGU(IR)
                  LAPSP = 2.0D0*LAPRHOU(IR)
               ELSE
                  D = 2.0D0*RHOD(IR)
                  AGRDSP = 2.0D0*AGRDRHOD(IR)
                  GDGAGDSP = 2.0D0*2.0D0*GDGAGD(IR)
                  LAPSP = 2.0D0*LAPRHOD(IR)
               END IF
C
               FK = CONF*D**THRD
               SS = AGRDSP/(D*2.0D0*FK)
               UU = GDGAGDSP/(D**2*(2.0D0*FK)**3)
               VV = LAPSP/(D*(2.0D0*FK)**2)
C
               IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
                  CALL EXCPBEX(D,SS,UU,VV,EX,VX,LLDA)
               ELSE IF ( SCFVXC(1:6).EQ.'EV-GGA' ) THEN
                  CALL EV92EX(D,SS,UU,VV,EX,VX)
               END IF
C
               EXC = EXC + EX*(D/2.0D0)/RHO(IR)
               IF ( JSP.EQ.1 ) VXCUP = VX
               IF ( JSP.EQ.2 ) VXCDN = VX
            END DO
C --------------------------------------------------------- correlation
            D = RHO(IR)
            ZET = (RHOU(IR)-RHOD(IR))/RHO(IR)
            RS = CONRS/D**THRD
C
            FK = 1.91915829D0/RS
            SK = DSQRT(4.D0*FK/PI)
            G = ((1.D0+ZET)**THRD2+(1.D0-ZET)**THRD2)/2.D0
            TT = AGRDRHO(IR)/(D*2.0D0*SK*G)
            UU = GDGAG(IR)/(D**2*(2.0D0*SK*G)**3)
            VV = LAPRHO(IR)/(D*(2.0D0*SK*G)**2)
            WW = AGRDRHO(IR)*(AGRDRHOU(IR)-AGRDRHOD(IR)-ZET*AGRDRHO(IR))
     &           /(D**2*(2.0D0*SK*G)**2)
C
            IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
               CALL EXCPBEC(RS,ZET,TT,UU,VV,WW,EC,VCUP,VCDN,LLDA)
            ELSE IF ( SCFVXC(1:6).EQ.'EV-GGA' ) THEN
               CALL CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET)
            ELSE
               CALL STOP_MESSAGE(ROUTINE,'SCFVXC not allowed')
            END IF
C
            EXC = EXC + EC
            VXCUP = VXCUP + VCUP
            VXCDN = VXCDN + VCDN
C
C ------------------------------------------------ convert from H to Ry
            EXC = 2.0D0*EXC
            XU = 2.0D0*VXCUP
            XD = 2.0D0*VXCDN
C
            V(IR,1) = V(IR,1) + XU
            V(IR,2) = V(IR,2) + XD
            U = W(IR,1)
            W(IR,1) = RHO(IR)*(U+EXC)
            W(IR,2) = RHO(IR)*U - 3D0*(EXC-XU)*RHOU(IR) - 3D0*(EXC-XD)
     &                *RHOD(IR)
            E(IR) = EXC
         END IF
      END DO
      END
C*==ev92ex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EV92EX(RHO,S,U,V,EX,VX)
C  This subroutine evaluates the exchange-only energy per electron
C  and exchange-only potential for the GGA presented in the paper
C  E. Engel, S. H. Vosko, Phys. Rev. B47, 13164 (1993), Eq.(21).
C  Its structure and use are completely equivalent to the subroutine
C  EXCH of the PW91 subroutine package available from J. P. Perdew.
C  Thus this subroutine may replace EXCH in order to be used together
C  with the PW91 form for the correlation energy functional in the
C  local density approximation, i.e. eliminating Perdew's subroutine
C  CORGGA.
C
C  INPUT:  D = DENSITY
C          S = ABS(GRAD D)/(2*KF*D)
C          U = (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C          V = (LAPLACIAN D)/(D*(2*KF)**2)
C  OUTPUT:  EX = EXCHANGE ENERGY PER ELECTRON (in atomic units)
C           VX = EXCHANGE POTENTIAL (in atomic units)
      IMPLICIT NONE
C*--EV92EX171
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 CVXLDA,CEXLDA,A1,A2,A3,B1,B2,B3
      PARAMETER (CVXLDA=-0.9847450218426966D0,CEXLDA=-0.738558766D0,
     &           A1=1.64712732976569D0,A2=.9801182732515644D0,
     &           A3=1.739938440729993D-02,B1=1.52367053964223D0,
     &           B2=.3672289494300232D0,B3=1.128174995137635D-02)
C
C Dummy arguments
C
      REAL*8 EX,RHO,S,U,V,VX
C
C Local variables
C
      REAL*8 A,AP,APP,B,BP,BPP,DROOT,F,FP,FPP,RA,RAP,RB,RBP,XI
C
C*** End of declarations rewritten by SPAG
C
C Dummy arguments
C
C Local variables
C
C
      IF ( RHO.GT.0.D0 ) THEN
C
         DROOT = RHO**0.3333333333333333D0
         XI = S**2
C
         A = 1.D0 + A1*XI + A2*XI**2 + A3*XI**3
         AP = A1 + 2.D0*A2*XI + 3.D0*A3*XI**2
         APP = 2.D0*A2 + 6.D0*A3*XI
         B = 1.D0 + B1*XI + B2*XI**2 + B3*XI**3
         BP = B1 + 2.D0*B2*XI + 3.D0*B3*XI**2
         BPP = 2.D0*B2 + 6.D0*B3*XI
C
         RA = AP/A
         RAP = APP/A - RA*RA
         RB = BP/B
         RBP = BPP/B - RB*RB
         F = A/B
         FP = F*(RA-RB)
         FPP = FP*(RA-RB) + F*(RAP-RBP)
C
         VX = CVXLDA*DROOT*(F-1.5D0*V*FP-(3.D0*U-4.D0*S**3)*S*FPP)
         EX = CEXLDA*DROOT*F
C
      ELSE
C
         VX = 0.D0
         EX = 0.D0
C
      END IF
C
      END
C*==excpbex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCPBEX(RHO,S,U,V,EX,VX,LLDA)
C----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  K Burke's modification of PW91 codes, May 14, 1996
C  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
C----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
C   (for U,V, see PW86(24))
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
C----------------------------------------------------------------------
C References:
C [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
C     Phys. Rev. Lett. 77, 3865 (1996).
C [b] J.P. Perdew and Y. Wang,
C     Phys. Rev. B33, 8800 (1986); B40, 3399 (1989) (E).
C----------------------------------------------------------------------
C Formulas:
C   	e_x[unif]=ax*rho^(4/3)  [LDA]
C ax = -0.75*(3/pi)^(1/3)
C	e_x[PBE]=e_x[unif]*FxPBE(s)
C	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
C uk, ul defined after [a](13)
C----------------------------------------------------------------------
C
C  All input and output is in atomic units!
C
C  Modifications by: E. Engel
C  Last revision:    May 9, 2001
Cengel
      IMPLICIT NONE
C*--EXCPBEX277
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 THRD,THRD4,AX,UM,UK,UL
      PARAMETER (THRD=1.D0/3.D0,THRD4=4.D0/3.D0,
     &           AX=-0.738558766382022405884230032680836D0,
     &           UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
C
C Dummy arguments
C
      REAL*8 EX,RHO,S,U,V,VX
      INTEGER LLDA
C
C Local variables
C
      REAL*8 EXUNIF,FS,FSS,FXPBE,P0,S2
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C construct LDA exchange energy density
      EXUNIF = AX*RHO**THRD
      IF ( LLDA.EQ.1 ) THEN
         EX = EXUNIF
         VX = EX*THRD4
         RETURN
      END IF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C construct PBE enhancement factor
      S2 = S*S
      P0 = 1.D0 + UL*S2
      FXPBE = 1D0 + UK - UK/P0
      EX = EXUNIF*FXPBE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
C  find first and second derivatives of Fx w.r.t s.
C  Fs=(1/s)*d FxPBE/ ds
C  Fss=d Fs/ds
      FS = 2.D0*UK*UL/(P0*P0)
      FSS = -4.D0*UL*S*FS/P0
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C calculate potential from [b](24)
      VX = EXUNIF*(THRD4*FXPBE-(U-THRD4*S2*S)*FSS-V*FS)
      END
C*==corlsd.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE CORLSD(RTRS,ZET,EC,VCUP,VCDN,ECRS,ECZET)
C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
C  POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
      IMPLICIT NONE
C*--CORLSD349
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 GAM,FZZ,THRD,THRD4
      PARAMETER (GAM=0.5198421D0,FZZ=1.709921D0,THRD=1.D0/3.D0,
     &           THRD4=4.D0/3.D0)
C
C Dummy arguments
C
      REAL*8 EC,ECRS,ECZET,RTRS,VCDN,VCUP,ZET
C
C Local variables
C
      REAL*8 ALFM,ALFRSM,COMM,EP,EPRS,EU,EURS,F,FZ,P1,P2,PP1,PP2,Z4
C
C*** End of declarations rewritten by SPAG
C
C      DATA GAM,FZZ/0.5198421D0,1.709921D0/
C      DATA THRD,THRD4/0.333333333333D0,1.333333333333D0/
C Dummy variables
C Local variables
C
      P1 = (1.D0+ZET)
      P2 = (1.D0-ZET)
      IF ( ABS(P1).GT.1D-12 ) THEN
         PP1 = P1**THRD
         P1 = PP1*PP1*PP1*PP1
      ELSE
         PP1 = 0
         P1 = 0
      END IF
      IF ( ABS(P2).GT.1D-12 ) THEN
         PP2 = P2**THRD
         P2 = PP2*PP2*PP2*PP2
      ELSE
         PP2 = 0
         P2 = 0
      END IF
      F = (P1+P2-2.D0)/GAM
C
      CALL EXCGCOR2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &              0.49294D0,RTRS,EU,EURS)
      CALL EXCGCOR2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &              0.62517D0,RTRS,EP,EPRS)
      CALL EXCGCOR2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     &              0.49671D0,RTRS,ALFM,ALFRSM)
C
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4) + EP*F*Z4 - ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.D0-Z4)/FZZ
C      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      FZ = THRD4*(PP1-PP2)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)
     &        + FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
      COMM = EC - RTRS*ECRS/3.D0 - ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      END
C*==excpbec.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCPBEC(RS,ZETA,T,UU,VV,WW,EC,VCUP,VCDN,LLDA)
C  This subroutine evaluates the correlation energy per particle and
C  spin-up and spin-dn correlation potentials within the Perdew-Burke-
C  Ernzerhof GGA. It is a slightly modified version of K. Burke's
C  official PBE subroutine.
C
C  Input:  RS   = WIGNER-SEITZ RADIUS = ( 3 / (4*PI*(DUP+DDN)) )**(1/3)
C          ZETA = RELATIVE SPIN POLARIZATION = (DUP-DDN)/(DUP+DDN)
C          T    = ABS(GRAD D) / ( (2*SK*G) * D )
C          UU   = (GRAD D)*GRAD(ABS(GRAD D)) / ( (2*SK*G)**3 * D**2 )
C          VV   = (LAPLACIAN D) / ( (2*SK*G)**2 * D )
C          WW   = (GRAD D)*(GRAD ZETA) / ( (2*SK*G)**2 * D )
C  where:  FK   = LOCAL FERMI MOMENTUM = (3*PI**2*(DUP+DDN))**(1/3)
C          SK   = LOCAL SCREENING MOMENTUM = (4*FK/PI)**(1/2)
C
C  Output: EC   = correlation energy per particle
C          VCUP = spin-up correlation potential
C          VCDN = spin-dn correlation potential
C
C  All input and output is in atomic units!
C
C References:
C [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
C     Phys. Rev. Lett. 77, 3865 (1996).
C [b] J. P. Perdew, K. Burke, and Y. Wang,
C     Phys. Rev. B54, 16533 (1996).
C [c] J. P. Perdew and Y. Wang,
C     Phys. Rev. B45, 13244 (1992).
C
C
C  Last revision:    May 9, 2001
C  Written by:       K. Burke, May 14, 1996.
C  Modifications by: E. Engel
Cengel
      IMPLICIT NONE
C*--EXCPBEC461
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 THRD,THRDM,THRD2,SIXTHM,THRD4,GAM,FZZ,GAMMA,BET,DELT,ETA
      PARAMETER (THRD=1.D0/3.D0,THRDM=-THRD,THRD2=2.D0*THRD,
     &           SIXTHM=THRDM/2.D0,THRD4=4.D0*THRD,
     &           GAM=0.5198420997897463295344212145565D0,
     &           FZZ=8.D0/(9.D0*GAM),
     &           GAMMA=0.03109069086965489503494086371273D0,
     &           BET=0.06672455060314922D0,DELT=BET/GAMMA,ETA=1.D-12)
C
C Dummy arguments
C
      REAL*8 EC,RS,T,UU,VCDN,VCUP,VV,WW,ZETA
      INTEGER LLDA
C
C Local variables
C
      REAL*8 ALFM,ALFRSM,B,B2,BEC,BG,COMM,ECRS,ECZETA,EP,EPRS,EU,EURS,F,
     &       FAC,FACT0,FACT1,FACT2,FACT3,FACT5,FZ,G,G3,G4,GZ,H,HB,HBT,
     &       HRS,HRST,HT,HTT,HZ,HZT,PON,PREF,Q4,Q5,Q8,Q9,RSTHRD,RTRS,T2,
     &       T4,T6,Z4
      EXTERNAL EXCGCOR2
C
C*** End of declarations rewritten by SPAG
C
C thrd*=various multiples of 1/3
C numbers for use in LSD energy spin-interpolation formula, [c](9).
C      GAM= 2^(4/3)-2
C      FZZ=f''(0)= 8/(9*GAM)
C numbers for construction of PBE
C      gamma=(1-log(2))/pi^2
C      bet=coefficient in gradient expansion for correlation, [a](4).
C      eta=small number to stop d phi/ dzeta from blowing up at
C          |zeta|=1.
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C find LSD energy contributions, using [c](10) and Table I[c].
C EU=unpolarized LSD correlation energy
C EURS=dEU/drs
C EP=fully polarized LSD correlation energy
C EPRS=dEP/drs
C ALFM=-spin stiffness, [c](3).
C ALFRSM=-dalpha/drs
C F=spin-scaling factor from [c](9).
C construct ec, using [c](8)
      IF ( RS.LT.3.D5 ) THEN
         RTRS = SQRT(RS)
         CALL EXCGCOR2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &                 0.49294D0,RTRS,EU,EURS)
         CALL EXCGCOR2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,
     &                 3.3662D0,0.62517D0,RTRS,EP,EPRS)
         CALL EXCGCOR2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,
     &                 0.88026D0,0.49671D0,RTRS,ALFM,ALFRSM)
         Z4 = ZETA**4
         F = ((1.D0+ZETA)**THRD4+(1.D0-ZETA)**THRD4-2.D0)/GAM
         EC = EU*(1.D0-F*Z4) + EP*F*Z4 - ALFM*F*(1.D0-Z4)/FZZ
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C LSD potential from [c](A1)
C ECRS = dEc/drs [c](A2)
C ECZETA=dEc/dzeta [c](A3)
C FZ = dF/dzeta [c](A4)
         ECRS = EURS*(1.D0-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.D0-Z4)/FZZ
         FZ = THRD4*((1.D0+ZETA)**THRD-(1.D0-ZETA)**THRD)/GAM
         ECZETA = 4.D0*(ZETA**3)*F*(EP-EU+ALFM/FZZ)
     &            + FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
         COMM = EC - RS*ECRS/3.D0 - ZETA*ECZETA
         VCUP = COMM + ECZETA
         VCDN = COMM - ECZETA
         IF ( LLDA.EQ.1 ) RETURN
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C PBE correlation energy
C G=phi(zeta), given after [a](3)
C DELT=bet/gamma
C B=A of [a](8)
         G = ((1.D0+ZETA)**THRD2+(1.D0-ZETA)**THRD2)/2.D0
         G3 = G**3
         PON = -EC/(G3*GAMMA)
         B = DELT/(EXP(PON)-1.D0)
         B2 = B*B
         T2 = T*T
         T4 = T2*T2
         Q4 = 1.D0 + B*T2
         Q5 = 1.D0 + B*T2 + B2*T4
         H = G3*(BET/DELT)*LOG(1.D0+DELT*Q4*T2/Q5)
         EC = EC + H
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
         G4 = G3*G
         T6 = T4*T2
         RSTHRD = RS/3.D0
         GZ = (((1.D0+ZETA)**2+ETA)**SIXTHM-((1.D0-ZETA)**2+ETA)
     &        **SIXTHM)/3.D0
         FAC = DELT/B + 1.D0
         BG = -3.D0*B2*EC*FAC/(BET*G4)
         BEC = B2*FAC/(BET*G3)
         Q8 = Q5*Q5 + DELT*Q4*Q5*T2
         Q9 = 1.D0 + 2.D0*B*T2
         HB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
         HRS = -RSTHRD*HB*BEC*ECRS
         FACT0 = 2.D0*DELT - 6.D0*B
         FACT1 = Q5*Q9 + Q4*Q9*Q9
         HBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
         HRST = RSTHRD*T2*HBT*BEC*ECRS
         HZ = 3.D0*GZ*H/G + HB*(BG*GZ+BEC*ECZETA)
         HT = 2.D0*BET*G3*Q9/Q8
         HZT = 3.D0*GZ*HT/G + HBT*(BG*GZ+BEC*ECZETA)
         FACT2 = Q4*Q5 + B*T2*(Q4*Q9+Q5)
         FACT3 = 2.D0*B*Q5*Q9 + DELT*FACT2
         HTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
         COMM = H + HRS + HRST + T2*HT/6.D0 + 7.D0*T2*T*HTT/6.D0
         PREF = HZ - GZ*T2*HT/G
         FACT5 = GZ*(2.D0*HT+T*HTT)/G
         COMM = COMM - PREF*ZETA - UU*HTT - VV*HT - WW*(HZT-FACT5)
         VCUP = VCUP + COMM + PREF
         VCDN = VCDN + COMM - PREF
      ELSE
         VCUP = 0.D0
         VCDN = 0.D0
      END IF
      END
C*==excgcor2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EXCGCOR2(A,A1,B1,B2,B3,B4,RTRS,GG,GGRS)
C----------------------------------------------------------------------
C######################################################################
C----------------------------------------------------------------------
C slimmed down version of GCOR used in PW91 routines, to interpolate
C LSD correlation energy, as given by (10) of
C J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
C K. Burke, May 11, 1996.
      IMPLICIT NONE
C*--EXCGCOR2612
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,A1,B1,B2,B3,B4,GG,GGRS,RTRS
C
C Local variables
C
      REAL*8 Q0,Q1,Q2,Q3
C
C*** End of declarations rewritten by SPAG
C
      Q0 = -2.D0*A*(1.D0+A1*RTRS*RTRS)
      Q1 = 2.D0*A*RTRS*(B1+RTRS*(B2+RTRS*(B3+B4*RTRS)))
      Q2 = LOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RTRS+2.D0*B2+RTRS*(3.D0*B3+4.D0*B4*RTRS))
      GGRS = -2.D0*A*A1*Q2 - Q0*Q3/(Q1*(1.D0+Q1))
      END
