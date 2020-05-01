C*==chivwn.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIVWN(RHOCHR,RHOSPN,EMM,ENM,ENN,IT,JTOP)
C   ********************************************************************
C   *                                                                  *
C   *  SUBROUTINE TO CALCULATE THE CURVATURE OF EXC WITH RESPECT       *
C   *  TO THE SPIN AND CHARGE DENSITIES AT S >= 0                      *
C   *  NEEDED FOR THE CALCULATION OF THE ENHANCEMENT OF THE            *
C   *  SPIN SUSCEPTIBILITY BY EXCHANGE-CORRELATION EFFECTS             *
C   *  ( BASED ON THE ROUTINE UXCOR ORIGINALLLY CODED BY M. MANNIEN    *
C   *  USING THE PARAMETRIZATION OF VOSKO, WILK AND NUSAIR FOR THE     *
C   *  EXCHANGE-CORRELATION ENERGY, COMMENTS REFER TO NOTATION AND     *
C   *  EQUATIONS IN THEIR ORIGINAL PAPER CAN.J.PHYS. 58, 1200 (1980) ) *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,LTXT_T
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET,IPRINT
      IMPLICIT NONE
C*--CHIVWN21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT,JTOP
      REAL*8 EMM(NRMAX),ENM(NRMAX),ENN(NRMAX),RHOCHR(NRMAX),
     &       RHOSPN(NRMAX)
C
C Local variables
C
      REAL*8 AF,ALPHAB,ALPHACV1,ALPHACV2,ALPHAV,ALPHCB1,ALPHCB2,ALPHXB,
     &       ALPHXV,AP,ATNF,ATNP,BETA,BF,BP,CBRT1,CBRT2,CF,CF1,CF2,CF3,
     &       CP,CP1,CP2,CP3,D2BETA,D2EFRHO,D2EPRHO,D2FS,D2FSX,DBETA,
     &       DBETAX,DCBRT1,DCBRT2,DDUC1X,DDUC2X,DDUCX,DEFRHO,DEFX,
     &       DEPRHO,DEPX,DFS,DFSS,DFSX,DS4X,DSX,DTF1X,DTP1X,DUC0X,
     &       DUC10X,DUC1RHO,DUC1X,DUC20X,DUC2RHO,DUC2X,DUCFX,DUCPX,
     &       DUX1RHO,DUX1X,DUX2RHO,DUX2X,DUXC1RHO,DUXC2RHO,ECF,ECP,FS,
     &       QF,QP,RHO,RHOS,RS,S,S4,TF1,TP1,UCF,UCP,X,XF0,XFX
      CHARACTER*80 FILNAM
      INTEGER I,LFN
      LOGICAL NEGVALS
      REAL*8 XP0,XPX
C
C*** End of declarations rewritten by SPAG
C
      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814D0, - 0.10498D0,
     &     3.72744D0,12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0,
     &     - 0.031167608D0/
      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907D0, - 0.32500D0,
     &     7.06042D0,18.0578D0,4.7309269D0,2.9847935D0,2.7100059D0,
     &     - 0.1446006D0/
C
C
      NEGVALS = .FALSE.
C
      IF ( IPRINT.GT.0 ) THEN
         CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_VWN_A-xc.dat',
     &                 13,FILNAM,LFN,1,IOTMP,'A-xc file ',10,NTMAX)
         WRITE (IOTMP,99001)
      END IF
C
      DO I = 1,JTOP
C
         RHO = RHOCHR(I)/(4.0D0*PI)
         RHOS = RHOSPN(I)/(4.0D0*PI)
         RS = (3D0/(4D0*PI*RHO))**(1D0/3D0)
         S = RHOS/RHO
C
C -------------------------------------------------------- VWN variables
C
         X = SQRT(RS)
         XPX = X*X + BP*X + CP
         XFX = X*X + BF*X + CF
         S4 = S**4 - 1D0
         CBRT1 = (1D0+S)**(1D0/3D0)
         CBRT2 = (1D0-S)**(1D0/3D0)
C
C ---------------------------------------------------------------- ds/dx
C
         DSX = -(S/RHO)*(-6D0*RHO/X)
C
C --------------------------------------------------------------- ds4/dx
C
         DS4X = 4D0*S**3*DSX
C
C -------------------------------------------- d(cbrt1)/ds , d(cbrt2)/ds
C
         DCBRT1 = (1D0/3D0)*(1D0+S)**(-2D0/3D0)
         DCBRT2 = -(1D0/3D0)*(1D0-S)**(-2D0/3D0)
C
C ----------------------------------------------------------------- f(s)
C
         FS = ((1D0+S)**(4D0/3D0)+(1D0-S)**(4D0/3D0)-2D0)
     &        /(2D0**(4D0/3D0)-2D0)
C
C ----------------------------------------------------------- d(f(s))/ds
C
         DFS = 4D0/3D0*(CBRT1-CBRT2)/(2D0**(4D0/3D0)-2D0)
         DFSS = DFS
C
C ----------------------------------------------------------- d(f(s))/dx
C
         DFSX = (4D0/3D0)*DSX*((1D0+S)**(1D0/3D0)-(1D0-S)**(1D0/3D0))
     &          /(2D0**(4D0/3D0)-2D0)
C
C ------------------------------------------------------- d^2(f(s))/ds^2
C
         D2FS = 4D0/9D0*((1.D0+S)**(-2D0/3D0)+(1.D0-S)**(-2D0/3D0))
     &          /(2D0**(4D0/3D0)-2D0)
C
C ------------------------------------------------------- d^2(f(s))/dsdx
C
         D2FSX = 4D0/9D0*DSX*((1.D0+S)**(-2D0/3D0)+(1.D0-S)**(-2D0/3D0))
     &           /(2D0**(4D0/3D0)-2D0)
C
C ----------------------------------------------------------------- beta
C
         BETA = 1D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
C
C ------------------------------------------------------------- dbeta/dx
C
         DBETAX = -2.D0*X*(0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
C
C ----------------------------------------------------------- d(beta)/ds
C
         DBETA = -(0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
         D2BETA = -(0.27402D0-1.591D0/X*X)*BETA**2 - 
     &            (0.27402D0*X+0.09873D0+1.591D0/X)*2*BETA*DBETAX
C
         ATNP = ATAN(QP/(2D0*X+BP))
         ATNF = ATAN(QF/(2D0*X+BF))
         ECP = AP*(LOG(X*X/XPX)
     &         +CP1*ATNP-CP3*(LOG((X-XP0)**2/XPX)+CP2*ATNP))
         ECF = AF*(LOG(X*X/XFX)
     &         +CF1*ATNF-CF3*(LOG((X-XF0)**2/XFX)+CF2*ATNF))
         TP1 = (X*X+BP*X)/XPX
         TF1 = (X*X+BF*X)/XFX
         UCP = ECP - AP/3D0*(1D0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))
         UCF = ECF - AF/3D0*(1D0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))
C
C ------------------------------------------------------------- d2e/dndn
C
         DEPRHO = -AP/3D0*(1D0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))/RHO
         DEPX = (-6D0*RHO/X)*DEPRHO
         DEFRHO = -AF/3D0*(1D0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))/RHO
         DEFX = (-6D0*RHO/X)*DEFRHO
C
C
         DTP1X = (2D0*X+BP)*CP/XPX**2
         DTF1X = (2D0*X+BF)*CF/XFX**2
C
         D2EPRHO = -
     &             AP/3D0*(-DTP1X-CP3*(-XP0/(X-XP0)-DTP1X-XP0*((CP-X*X)/
     &             XPX**2)))
         D2EFRHO = -
     &             AF/3D0*(-DTF1X-CF3*(-XF0/(X-XF0)-DTF1X-XF0*((CF-X*X)/
     &             XFX**2)))
C
         DUCPX = DEPX + D2EPRHO
         DUCFX = DEFX + D2EFRHO
C
         DUC0X = DUCPX + (DUCFX-DUCPX)*FS + (UCF-UCP)*DFSX
         DUC10X = DUC0X - (DEFX-DEPX)*(S-1D0)*DFSX - (ECF-ECP)
     &            *DSX*DFSX - (ECF-ECP)*(S-1D0)*D2FSX
         DUC20X = DUC0X - (DEFX-DEPX)*(S+1D0)*DFSX - (ECF-ECP)
     &            *DSX*DFSX - (ECF-ECP)*(S+1D0)*D2FSX
C
         DDUCX = (DUCFX-DUCPX)*BETA*S4*FS + (UCF-UCP)*DBETAX*S4*FS + 
     &           (UCF-UCP)*BETA*DS4X*FS + (UCF-UCP)*BETA*S4*DFSX + 
     &           (DEFX-DEPX)*(-RS/3D0)*DBETA*S4*FS + (ECF-ECP)
     &           *(-X*2D0/3D0)*DBETA*S4*FS + (ECF-ECP)*(-RS/3D0)
     &           *D2BETA*S4*FS + (ECF-ECP)*(-RS/3D0)*DBETA*DS4X*FS + 
     &           (ECF-ECP)*(-RS/3D0)*DBETA*S4*DFSX
C
         DDUC1X = DDUCX - (DEFX-DEPX)*BETA*(S-1D0)*(4D0*S**3*FS+S4*DFS)
     &            - (ECF-ECP)*DBETAX*(S-1D0)*(4D0*S**3*FS+S4*DFS)
     &            - (ECF-ECP)*BETA*DSX*(4D0*S**3*FS+S4*DFS) - (ECF-ECP)
     &            *BETA*(S-1D0)
     &            *(12D0*S**2*DSX*FS+4D0*S**3*DFSX+DS4X*DFS+S4*D2FSX)
C
         DDUC2X = DDUCX - (DEFX-DEPX)*BETA*(S+1D0)*(4D0*S**3*FS+S4*DFS)
     &            - (ECF-ECP)*DBETAX*(S+1D0)*(4D0*S**3*FS+S4*DFS)
     &            - (ECF-ECP)*BETA*DSX*(4D0*S**3*FS+S4*DFS) - (ECF-ECP)
     &            *BETA*(S+1D0)
     &            *(12D0*S**2*DSX*FS+4D0*S**3*DFSX+DS4X*DFS+S4*D2FSX)
C
         DUC1X = DUC10X + DDUC1X
         DUC2X = DUC20X + DDUC2X
         DUC1RHO = (-X/(6D0*RHO))*DUC1X
         DUC2RHO = (-X/(6D0*RHO))*DUC2X
C
         DUX1X = (1.221774D0/3D0)*DSX*(1+S)**(-2D0/3D0)
     &           /RS - 2.D0*1.221774D0*(1+S)**(1D0/3D0)/X**3
C
         DUX2X = -(1.221774D0/3D0)*DSX*(1-S)**(-2D0/3D0)
     &           /RS - 2.D0*1.221774D0*(1-S)**(1D0/3D0)/X**3
C
         DUX1RHO = (-X/(6D0*RHO))*DUX1X
         DUX2RHO = (-X/(6D0*RHO))*DUX2X
C
         DUXC1RHO = DUC1RHO - DUX1RHO
         DUXC2RHO = DUC2RHO - DUX2RHO
C
         ENN(I) = (DUXC1RHO+DUXC2RHO)/2.D0
C
C ------------------------------------------------------------- d2e/dndm
C ----------------------------------------------- exchange contributions
C
         ALPHXV = -1.221774D0/RS*(DCBRT1+DCBRT2)/2.D0
C
C ----------------------------- correlation contrribution to dVxc: 1 & 2
C
         ALPHACV1 = (UCF-UCP)*DFSS - (ECF-ECP)*DFSS - (ECF-ECP)*S*D2FS
         ALPHACV2 = (UCF-UCP)*BETA*4.D0*S**3*FS + (UCF-UCP)
     &              *BETA*S4*DFSS - (ECF-ECP)*(RS/3.D0)
     &              *DBETA*4.D0*S**3.D0*FS - (ECF-ECP)*(RS/3.D0)
     &              *DBETA*S4*DFSS - (ECF-ECP)
     &              *BETA*(16.D0*S**3*FS+4.D0*S**4*DFSS+
     &              (5D0*S**4.D0-1.D0)*DFSS+(S**5.D0-S)*D2FS)
         ALPHAV = ALPHACV1 + ALPHACV2 + ALPHXV
         ENM(I) = -ALPHAV/RHO
C
C ------------------------------------------------------------- d2e/dmdm
C ----------------------------------------------- exchange contributions
C
         ALPHXB = -1.221774D0/RS*(DCBRT1-DCBRT2)/2.D0
C
C ------------------------------- correlation contribution to Bxc: 1 & 2
C
         ALPHCB1 = (ECF-ECP)*D2FS
         ALPHCB2 = (ECF-ECP)
     &             *BETA*(12.D0*S*S*FS+8.D0*S**3.D0*DFSS+(S**4.D0-1)
     &             *D2FS)
         ALPHAB = ALPHXB + ALPHCB1 + ALPHCB2
         EMM(I) = -ALPHAB/RHO
C
C *** check for negative value
         IF ( EMM(I).LT.0.0D0 ) THEN
            EMM(I) = 0.0D0
            NEGVALS = .TRUE.
         END IF
C
         IF ( IPRINT.GT.0 ) WRITE (IOTMP,'(1I4,1X,6E12.4)') I,RHOCHR(I),
     &                             RS,EMM(I),ENM(I),ENN(I)
      END DO
C
      IF ( NEGVALS ) WRITE (6,*) 
     &        'WARNING !!! NEGATIVE VALUES FOR STONER KERNEL IN CHIVWN '
     &        ,'FOR IT=',IT
C
      CLOSE (IOTMP)
99001 FORMAT ('# FILE fort.(520+IT)',/,
     &        '# CHIVWN: I,RHOCHR,RS,EMM,ENM,ENN')
C
      END
