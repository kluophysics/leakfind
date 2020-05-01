C*==uxcor.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE UXCOR(RS,S,UXC1,UXC2,EXC)
C-----------------------------------------------------------------------
C    ---- this subroutine was coded by m. mannien ----
C    subroutine ( m. manninen )  to calculate energy and potential
C    from ceperley-alder ( parametrization of vosko, wilk and nusair )
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*--UXCOR9
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EXC,RS,S,UXC1,UXC2
C
C Local variables
C
      REAL*8 AF,AP,ATNF,ATNP,BETA,BF,BP,CBRT1,CBRT2,CF,CF1,CF2,CF3,CP,
     &       CP1,CP2,CP3,DBETA,DFS,DUC,DUC1,DUC2,EC,ECF,ECP,FS,QF,QP,S4,
     &       TF1,TP1,UC0,UC1,UC10,UC2,UC20,UCF,UCP,X,XF0,XFX,XP0,XPX
C
C*** End of declarations rewritten by SPAG
C
      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814D0, - 0.10498D0,
     &     3.72744D0,12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0,
     &     - 0.031167608D0/
      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907D0, - 0.32500D0,
     &     7.06042D0,18.0578D0,4.7309269D0,2.9847935D0,2.7100059D0,
     &     - 0.1446006D0/
      X = SQRT(RS)
      XPX = X*X + BP*X + CP
      XFX = X*X + BF*X + CF
      S4 = S**4 - 1D0
      CBRT1 = (1D0+S)**(1D0/3D0)
      CBRT2 = (1D0-S)**(1D0/3D0)
      FS = ((1D0+S)**(4D0/3D0)+(1D0-S)**(4D0/3D0)-2D0)
     &     /(2D0**(4D0/3D0)-2D0)
      BETA = 1D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
      DFS = 4D0/3D0*(CBRT1-CBRT2)/(2D0**(4D0/3D0)-2D0)
      DBETA = -(0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
      ATNP = ATAN(QP/(2D0*X+BP))
      ATNF = ATAN(QF/(2D0*X+BF))
      ECP = AP*(LOG(X*X/XPX)+CP1*ATNP-CP3*(LOG((X-XP0)**2/XPX)+CP2*ATNP)
     &      )
      ECF = AF*(LOG(X*X/XFX)+CF1*ATNF-CF3*(LOG((X-XF0)**2/XFX)+CF2*ATNF)
     &      )
      EC = ECP + FS*(ECF-ECP)*(1D0+S4*BETA)
      TP1 = (X*X+BP*X)/XPX
      TF1 = (X*X+BF*X)/XFX
      UCP = ECP - AP/3D0*(1D0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))
      UCF = ECF - AF/3D0*(1D0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))
      UC0 = UCP + (UCF-UCP)*FS
      UC10 = UC0 - (ECF-ECP)*(S-1D0)*DFS
      UC20 = UC0 - (ECF-ECP)*(S+1D0)*DFS
      DUC = (UCF-UCP)*BETA*S4*FS + (ECF-ECP)*(-RS/3D0)*DBETA*S4*FS
      DUC1 = DUC - (ECF-ECP)*BETA*(S-1D0)*(4D0*S**3*FS+S4*DFS)
      DUC2 = DUC - (ECF-ECP)*BETA*(S+1D0)*(4D0*S**3*FS+S4*DFS)
      UC1 = UC10 + DUC1
      UC2 = UC20 + DUC2
      UXC1 = UC1 - 1.221774D0/RS*CBRT1
      UXC2 = UC2 - 1.221774D0/RS*CBRT2
      EXC = EC - 0.9163306D0/RS - 0.2381735D0/RS*FS
      END
