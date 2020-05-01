C*==scfscale_alat_eos.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE SCFSCALE_ALAT_EOS(NSCL)
C   ********************************************************************
C   *                                                                  *
C   *   plot the Equation of State  in case of   <SCFSCALE_ALAT>       *
C   *                                                                  *
C   *   the fit to the Murnaghan EOS is based on some driver routines  *
C   *   by   A. V. Postnikov  and  the MINPACK routine LMDIF           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_SCF,ONLY:ALAT_TAB,EFERMI_TAB,ETOT_TAB,MUEORB_TAB,
     &    MUESPN_TAB,VOL_TAB
      USE MOD_FILES,ONLY:DATSET,LDATSET,SYSTEM,LSYSTEM
      USE MOD_CONSTANTS,ONLY:RY_EV,EV_J,A0_SI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NLEGMAX,NPLOT
      PARAMETER (NLEGMAX=3,NPLOT=50)
C
C Dummy arguments
C
      INTEGER NSCL
C
C Local variables
C
      REAL*8 AA,ALAT_MIN,BB,B_MBAR,CC,DIAG(:),DX,EE,ENTHA,EPSFCN,FACTOR,
     &       FJAC(:,:),FTOL,PRESS,P_AU2MBAR,QTF(4),VOL_ALAT,WA1(4),
     &       WA2(4),WA3(4),WA4(:),X(4),X12,X23,X31,XMAX,XMIN,XPLOT(:),
     &       XX,Y(:),YMAX(2),YMIN(2),YPLOT(:)
      CHARACTER*80 FILNAM
      INTEGER IFIL,IMID,INFO,IPLOT,IPVT(4),ISCL,LDFJAC,LFILNAM,MAXFEV,
     &        MODE,N,NFEV,NGRAPH,NPRINT
      LOGICAL K_MURNAGHAN_FIT
      CHARACTER*20 LEG(NLEGMAX)
      EXTERNAL EOS_FIT_MURNAGHAN,EOS_MURNAGHAN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DIAG,FJAC,WA4,Y
      ALLOCATABLE XPLOT,YPLOT
C
      IF ( NSCL.LE.1 ) RETURN
C
      WRITE (6,99006)
C
      IFIL = 7
      ALLOCATE (XPLOT(NPLOT),YPLOT(NPLOT))
      ALLOCATE (DIAG(NSCL),FJAC(NSCL,4),WA4(NSCL),Y(NSCL))
C
C--- convert pressure from a.u. (Ry/a_0^3) to Pascal and finally to Mbar
C
      P_AU2MBAR = RY_EV*EV_J/(A0_SI)**3*1D-11
C
C=======================================================================
C                  magnetic moments and Fermi energy
C=======================================================================
C
      YMIN(1:2) = +1D+30
      YMAX(1:2) = -1D+30
C
      DO ISCL = 1,NSCL
         YMIN(1) = MIN(YMIN(1),EFERMI_TAB(ISCL))
         YMAX(1) = MAX(YMAX(1),EFERMI_TAB(ISCL))
      END DO
C
      IF ( IREL.GE.2 ) THEN
         DO ISCL = 1,NSCL
            YMIN(2) = MIN(YMIN(2),MUESPN_TAB(ISCL))
            YMIN(2) = MIN(YMIN(2),MUEORB_TAB(ISCL))
            YMAX(2) = MAX(YMAX(2),MUESPN_TAB(ISCL))
            YMAX(2) = MAX(YMAX(2),MUEORB_TAB(ISCL))
         END DO
         NGRAPH = 2
      ELSE
         NGRAPH = 1
      END IF
C
C ----------------------------------------------------------------------
C
      CALL XMGRHEAD(DATSET,LDATSET,'m-vs-a',6,'XXX',0,FILNAM,80,LFILNAM,
     &              IFIL,NGRAPH,ALAT_TAB(1),1,ALAT_TAB(NSCL),1,YMIN(1),
     &              1,YMAX(1),1,YMIN(2),1,YMAX(2),1,
     &              'lattice parameter a (a.u.)',26,'E!sFermi!N (Ry)',
     &              15,'!xm!F (!xm!F!sB!N)',18,
     &              'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &              25+LSYSTEM,
     &              'magnetic moments and Fermi energy as function of a'
     &              ,50,.FALSE.)
C
C----------------------------------------- upper panel: magnetic moments
C
      IF ( IREL.GE.2 ) THEN
C
         LEG(1) = '!xm!F!sspin!N '
         LEG(2) = '!xm!F!sorb!N '
C
         CALL XMGRLEG1(IFIL,1,2,LEG,0.18D0,0.83D0)
C
         CALL XMGRTABLE(1,0,ALAT_TAB,MUESPN_TAB,1.0D0,NSCL,IFIL)
         CALL XMGRTABLE(1,1,ALAT_TAB,MUEORB_TAB,1.0D0,NSCL,IFIL)
C
      END IF
C
C-------------------------------------------------- lower panel: E_Fermi
C
      CALL XMGRTABLE(0,0,ALAT_TAB,EFERMI_TAB,1.0D0,NSCL,IFIL)
C
      WRITE (6,99007) 'moments and E_Fermi',FILNAM(1:LFILNAM)
C
      CLOSE (IFIL)
C
C=======================================================================
C             fit of total energy to  Murnaghan  EOS
C=======================================================================
C
      CALL RVECSPAT(ABAS(1,1),ABAS(1,2),ABAS(1,3),VOL_ALAT,1)
C
      DO ISCL = 1,NSCL
         VOL_TAB(ISCL) = VOL_ALAT*ALAT_TAB(ISCL)**3
      END DO
C
      IF ( NSCL.LT.5 ) THEN
         WRITE (6,*) ' Number of data points =',NSCL,
     &               ' is too small for a reasonable fit'
         K_MURNAGHAN_FIT = .FALSE.
      ELSE
         K_MURNAGHAN_FIT = .TRUE.
      END IF
C
C-----------------------------------------------------------------------
      IF ( K_MURNAGHAN_FIT ) THEN
C
C --- initialize parameters to call subr. lmdif:
C     three-point parabolic fit: E = AA*vol**2 + BB*vol + CC
C
         IMID = NSCL/2
         X12 = VOL_TAB(1) - VOL_TAB(IMID)
         X23 = VOL_TAB(IMID) - VOL_TAB(NSCL)
         X31 = VOL_TAB(NSCL) - VOL_TAB(1)
         AA = -ETOT_TAB(1)/X12/X31 - ETOT_TAB(IMID)/X12/X23 - 
     &        ETOT_TAB(NSCL)/X31/X23
         BB = ETOT_TAB(1)*(VOL_TAB(IMID)+VOL_TAB(NSCL))/X12/X31 + 
     &        ETOT_TAB(IMID)*(VOL_TAB(1)+VOL_TAB(NSCL))/X12/X23 + 
     &        ETOT_TAB(NSCL)*(VOL_TAB(1)+VOL_TAB(IMID))/X31/X23
         CC = -ETOT_TAB(1)*VOL_TAB(IMID)*VOL_TAB(NSCL)/X12/X31 - 
     &        ETOT_TAB(IMID)*VOL_TAB(1)*VOL_TAB(NSCL)/X12/X23 - 
     &        ETOT_TAB(NSCL)*VOL_TAB(1)*VOL_TAB(IMID)/X31/X23
C
C     initial guess of Murnaghan parameters from 2d order E(V) fit:
C
         X(1) = CC - 0.25D0*BB*BB/AA
                                   !  min. energy
         X(2) = -0.5D0*BB/AA       !  equilibrium volume
         X(3) = -BB                !  bulk modulus
         X(4) = -10.0D0            !  dB/dP
         B_MBAR = X(3)*P_AU2MBAR
         ALAT_MIN = (X(2)/VOL_ALAT)**(1D0/3D0)
C
         WRITE (6,99001)
         WRITE (6,99004) ' ',X(2),' ',ALAT_MIN,' ',X(1),' ',B_MBAR,' ',
     &                   X(4)
C
         FTOL = 1.D-9
         MAXFEV = 10000
         EPSFCN = 0.001D0
         MODE = 1
         DIAG(1:NSCL) = 1D0
         NPRINT = 1
         FACTOR = 100.D0
         LDFJAC = 100
         N = 4
C
         CALL LMDIF(EOS_FIT_MURNAGHAN,NSCL,N,X,Y,FTOL,FTOL,FTOL,MAXFEV,
     &              EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,
     &              LDFJAC,IPVT,QTF,WA1,WA2,WA3,WA4)
C
         WRITE (6,99002) INFO
         WRITE (6,99003)
C
         DO IPLOT = 1,NSCL
            WRITE (6,99005) VOL_TAB(IPLOT),ETOT_TAB(IPLOT),Y(IPLOT)
         END DO
C
C  After fitting: x(1): min. energy in initial units
C                 x(2): equilibrium volume in initial units
C                 x(4): d(Bulk modulus)/d(Pressure) - dimensionless
C                 x(3): bulk modulus and needs conversion
C
         ALAT_MIN = (X(2)/VOL_ALAT)**(1D0/3D0)
C
         B_MBAR = X(3)*P_AU2MBAR
         WRITE (6,99004) ' ',X(2),' ',ALAT_MIN,' ',X(1),' ',B_MBAR,' ',
     &                   X(4)
C         WRITE (IFILSCLALAT,99004) '#',X(2),'#',ALAT_MIN,'#',X(1),'#',
C     &                             B_MBAR,'#',X(4)
C      IFILSCLALAT in #FILES#
C-----------------------------------------------------------------------
CC.....create plotfile for Murnaghan fit:
C      XMIN = VOL_TAB(1)
C      XMAX = VOL_TAB(NSCL)
C      DX = (XMAX-XMIN)/(NSCL-1)
C      WRITE (6,99005)
C      DO IPLOT = 1,NSCL
C         XX = XMIN + (IPLOT-1)*DX
C         CALL EOS_MURNAGHAN(X,XX,EE,PRESS,ENTHA)
C         WRITE (6,99009) XX,EE,PRESS,ENTHA
C      END DO
C
C99005 FORMAT (/,10X,'total energy fitted to Murnaghan EOS:',/10X,
C     &        '  Volume(a.u.^3)',7x,'Energy(Ry)',5x,'Pressure (GPa)',6x,
C     &        'Enthalpy(Ry)')
      END IF
C-----------------------------------------------------------------------
C
      YMIN(1:2) = +1D+30
      YMAX(1:2) = -1D+30
C
      DO ISCL = 1,NSCL
         YMIN(1) = MIN(YMIN(1),ETOT_TAB(ISCL))
         YMAX(1) = MAX(YMAX(1),ETOT_TAB(ISCL))
      END DO
C
C ----------------------------------------------------------------------
C
      CALL XMGRHEAD(DATSET,LDATSET,'EOS',3,'XXX',0,FILNAM,80,LFILNAM,
     &              IFIL,1,ALAT_TAB(1),1,ALAT_TAB(NSCL),1,YMIN(1),1,
     &              YMAX(1),1,YMIN(2),1,YMAX(2),1,
     &              'lattice parameter a (a.u.)',26,'E!stot!N (Ry)',13,
     &              '!xm!F (!xm!F!sB!N)',18,
     &              'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &              25+LSYSTEM,'total energy E!stot!N fitted to '//
     &              'Murnaghan equation of state',59,.FALSE.)
C
      WRITE (IFIL,99008) 's0 symbol 1'
      WRITE (IFIL,99008) 's0 symbol fill pattern 1'
      WRITE (IFIL,99008) 's0 line type 0'
      WRITE (IFIL,99008) 'yaxis  label place 0.000000, 0.200000'
      WRITE (IFIL,99008) 'yaxis  ticklabel prec 9'
      WRITE (IFIL,99008) 'view 0.250000, 0.150000, 0.950000, 0.850000'
C
      LEG(1) = 'calculated'
      LEG(2) = 'Murnaghan fit'
C
      CALL XMGRLEG1(IFIL,0,2,LEG,0.4D0,0.8D0)
C
C
      CALL XMGRTABLE(0,0,ALAT_TAB,ETOT_TAB,1.0D0,NSCL,IFIL)
C
C-----------------------------------------------------------------------
      IF ( K_MURNAGHAN_FIT ) THEN
C
         XMIN = VOL_TAB(1)*0.95D0
         XMAX = VOL_TAB(NSCL)*1.05D0
         DX = (XMAX-XMIN)/(NPLOT-1)
         DO IPLOT = 1,NPLOT
            XX = XMIN + (IPLOT-1)*DX
C
            CALL EOS_MURNAGHAN(X,XX,EE,PRESS,ENTHA)
C
            XPLOT(IPLOT) = (XX/VOL_ALAT)**(1D0/3D0)
            YPLOT(IPLOT) = EE
         END DO
C
         CALL XMGRTABLE(0,1,XPLOT,YPLOT,1.0D0,NPLOT,IFIL)
C
      END IF
C-----------------------------------------------------------------------
C
      WRITE (6,99007) 'E_tot and Murnaghan fit',FILNAM(1:LFILNAM)
C
      CLOSE (IFIL)
C
99001 FORMAT (/,10X,'fit to Murnaghan equation of state',//,10X,
     &        'rough estimate from the 2d order polynomial fit:')
99002 FORMAT (/,10X,'info code from fitting routine:',i4)
99003 FORMAT (/,10X,'input V (a.u.)^3',6x,'input E (Ry)',4x,
     &        'Absolute error')
99004 FORMAT (/,A,9X,'equilibrium volume           ',F16.6,' (a.u.)^3'/,
     &        A,9X,'equilibrium lattice parameter',F16.6,' a.u.',/,A,9X,
     &        'minimum energy               ',F16.6,' Ry',/,A,9X,
     &        'bulk modulus                 ',F16.6,' Mbar',/,A,9X,
     &        'BP                           ',F16.6)
99005 FORMAT (8X,4F18.8)
99006 FORMAT (//,1X,79('*'),/,30X,'<SCFSCALE_ALAT_EOS>',/,1X,79('*'),//)
99007 FORMAT (/,10X,'xmgrace-file for ',A,': ',A,/)
99008 FORMAT ('@   ',A)
C
      END
C*==eos_fit_murnaghan.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE EOS_FIT_MURNAGHAN(M,N,X,FVEC,IFLAG)
C   ********************************************************************
C     fit to the Murnaghan equation of state.
C
C     On fitting, the parameters x(1)...x(4) are:
C        x(1) : minimum energy E
C        x(2) : volume V at equilibrium
C        x(3) : bulk modulus B = V*d^2(E)/d(V)^2 in corresponding units
C        x(4) : d(B)/d(P), dimensionless
C   ********************************************************************
C
      USE MOD_SCF,ONLY:ETOT_TAB,VOL_TAB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFLAG,M,N
      REAL*8 FVEC(*),X(*)
C
C Local variables
C
      INTEGER I
      REAL*8 XH
C
C*** End of declarations rewritten by SPAG
C
      IFLAG = 0
C
      IF ( N.EQ.3 ) THEN
         XH = X(4)
         X(4) = 5.0D0
      END IF
      DO I = 1,M
         FVEC(I) = X(1) + X(3)*VOL_TAB(I)/(X(4)*(X(4)-1.D0))
     &             *(X(4)*(1.D0-X(2)/VOL_TAB(I))+(X(2)/VOL_TAB(I))**X(4)
     &             -1.D0) - ETOT_TAB(I)
      END DO
      IF ( N.EQ.3 ) X(4) = XH
      END
C*==eos_murnaghan.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE EOS_MURNAGHAN(X,VOL,ENE,PRESS,ENTHA)
C   ********************************************************************
C     fit to the Murnaghan equation of state.
C
C     On fitting, the parameters x(1)...x(4) are:
C        x(1) : minimum energy E in Ry
C        x(2) : volume V at equilibrium in a.u.^3
C        x(3) : bulk modulus B = V*d^2(E)/d(V)^2 in corresponding units
C        x(4) : d(B)/d(P), dimensionless
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:RY_EV,EV_J,A0_SI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ENE,ENTHA,PRESS,VOL
      REAL*8 X(4)
C
C Local variables
C
      REAL*8 P_AU2MBAR
C
C*** End of declarations rewritten by SPAG
C
C --- fit to energy:
      ENE = X(1) + X(3)*VOL/(X(4)*(X(4)-1.D0))
     &      *(X(4)*(1.D0-X(2)/VOL)+(X(2)/VOL)**X(4)-1.D0)
C
C --- fit to pressure:
      PRESS = X(3)/X(4)*((X(2)/VOL)**X(4)-1.D0)
C
C --- enthalpy:
      ENTHA = ENE + PRESS*VOL
C
C--- convert pressure from a.u. (Ry/a_0^3) to Pascal and finally to Mbar
C
      P_AU2MBAR = RY_EV*EV_J/(A0_SI)**3*1D-11
C
      PRESS = PRESS*P_AU2MBAR
C
      END
C*==lmdif.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE LMDIF(FCN,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,DIAG,
     &                 MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,IPVT,
     &                 QTF,WA1,WA2,WA3,WA4)
C   ********************************************************************
C
C     the purpose of lmdif is to minimize the sum of the squares of
C     m nonlinear functions in n variables by a modification of
C     the levenberg-marquardt algorithm. the user must provide a
C     subroutine which calculates the functions. the jacobian is
C     then calculated by a forward-difference approximation.
C
C     the subroutine statement is
C
C       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
C                        diag,mode,factor,nprint,info,nfev,fjac,
C                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
C
C     where
C
C       fcn is the name of the user-supplied subroutine which
C         calculates the functions. fcn must be declared
C         in an external statement in the user calling
C         program, and should be written as follows.
C
C         subroutine fcn(m,n,x,fvec,iflag)
C         integer m,n,iflag
C         double precision x(n),fvec(m)
C         ----------
C         calculate the functions at x and
C         return this vector in fvec.
C         ----------
C         return
C         end
C
C         the value of iflag should not be changed by fcn unless
C         the user wants to terminate execution of lmdif.
C         in this case set iflag to a negative integer.
C
C       m is a positive integer input variable set to the number
C         of functions.
C
C       n is a positive integer input variable set to the number
C         of variables. n must not exceed m.
C
C       x is an array of length n. on input x must contain
C         an initial estimate of the solution vector. on output x
C         contains the final estimate of the solution vector.
C
C       fvec is an output array of length m which contains
C         the functions evaluated at the output x.
C
C       ftol is a nonnegative input variable. termination
C         occurs when both the actual and predicted relative
C         reductions in the sum of squares are at most ftol.
C         therefore, ftol measures the relative error desired
C         in the sum of squares.
C
C       xtol is a nonnegative input variable. termination
C         occurs when the relative error between two consecutive
C         iterates is at most xtol. therefore, xtol measures the
C         relative error desired in the approximate solution.
C
C       gtol is a nonnegative input variable. termination
C         occurs when the cosine of the angle between fvec and
C         any column of the jacobian is at most gtol in absolute
C         value. therefore, gtol measures the orthogonality
C         desired between the function vector and the columns
C         of the jacobian.
C
C       maxfev is a positive integer input variable. termination
C         occurs when the number of calls to fcn is at least
C         maxfev by the end of an iteration.
C
C       epsfcn is an input variable used in determining a suitable
C         step length for the forward-difference approximation. this
C         approximation assumes that the relative errors in the
C         functions are of the order of epsfcn. if epsfcn is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
C
C       diag is an array of length n. if mode = 1 (see
C         below), diag is internally set. if mode = 2, diag
C         must contain positive entries that serve as
C         multiplicative scale factors for the variables.
C
C       mode is an integer input variable. if mode = 1, the
C         variables will be scaled internally. if mode = 2,
C         the scaling is specified by the input diag. other
C         values of mode are equivalent to mode = 1.
C
C       factor is a positive input variable used in determining the
C         initial step bound. this bound is set to the product of
C         factor and the euclidean norm of diag*x if nonzero, or else
C         to factor itself. in most cases factor should lie in the
C         interval (.1,100.). 100. is a generally recommended value.
C
C       nprint is an integer input variable that enables controlled
C         printing of iterates if it is positive. in this case,
C         fcn is called with iflag = 0 at the beginning of the first
C         iteration and every nprint iterations thereafter and
C         immediately prior to return, with x and fvec available
C         for printing. if nprint is not positive, no special calls
C         of fcn with iflag = 0 are made.
C
C       info is an integer output variable. if the user has
C         terminated execution, info is set to the (negative)
C         value of iflag. see description of fcn. otherwise,
C         info is set as follows.
C
C         info = 0  improper input parameters.
C
C         info = 1  both actual and predicted relative reductions
C                   in the sum of squares are at most ftol.
C
C         info = 2  relative error between two consecutive iterates
C                   is at most xtol.
C
C         info = 3  conditions for info = 1 and info = 2 both hold.
C
C         info = 4  the cosine of the angle between fvec and any
C                   column of the jacobian is at most gtol in
C                   absolute value.
C
C         info = 5  number of calls to fcn has reached or
C                   exceeded maxfev.
C
C         info = 6  ftol is too small. no further reduction in
C                   the sum of squares is possible.
C
C         info = 7  xtol is too small. no further improvement in
C                   the approximate solution x is possible.
C
C         info = 8  gtol is too small. fvec is orthogonal to the
C                   columns of the jacobian to machine precision.
C
C       nfev is an integer output variable set to the number of
C         calls to fcn.
C
C       fjac is an output m by n array. the upper n by n submatrix
C         of fjac contains an upper triangular matrix r with
C         diagonal elements of nonincreasing magnitude such that
C
C                t     t           t
C               p *(jac *jac)*p = r *r,
C
C         where p is a permutation matrix and jac is the final
C         calculated jacobian. column j of p is column ipvt(j)
C         (see below) of the identity matrix. the lower trapezoidal
C         part of fjac contains information generated during
C         the computation of r.
C
C       ldfjac is a positive integer input variable not less than m
C         which specifies the leading dimension of the array fjac.
C
C       ipvt is an integer output array of length n. ipvt
C         defines a permutation matrix p such that jac*p = q*r,
C         where jac is the final calculated jacobian, q is
C         orthogonal (not stored), and r is upper triangular
C         with diagonal elements of nonincreasing magnitude.
C         column j of p is column ipvt(j) of the identity matrix.
C
C       qtf is an output array of length n which contains
C         the first n elements of the vector (q transpose)*fvec.
C
C       wa1, wa2, and wa3 are work arrays of length n.
C
C       wa4 is a work array of length m.
C
C     subprograms called
C
C       user-supplied ...... fcn
C
C       minpack-supplied ... dpmpar,fdjac2,lmpar,qrfac
C
C       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
C
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
C
C     NIST Guide to Available Math Software.
C     Fullsource for module LMDIF from package MINPACK.
C     Retrieved from NETLIB on Tue Dec  1 11:06:31 1998.
C
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPSFCN,FACTOR,FTOL,GTOL,XTOL
      INTEGER INFO,LDFJAC,M,MAXFEV,MODE,N,NFEV,NPRINT
      REAL*8 DIAG(N),FJAC(LDFJAC,N),FVEC(M),QTF(N),WA1(N),WA2(N),WA3(N),
     &       WA4(M),X(N)
      INTEGER IPVT(N)
C
C Local variables
C
      REAL*8 ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM,ONE,P0001,P1,
     &       P25,P5,P75,PAR,PNORM,PRERED,RATIO,RSUM,TEMP,TEMP1,TEMP2,
     &       XNORM,ZERO
      REAL*8 DNRM2,DPMPAR
      INTEGER I,IFLAG,ITER,J,L
      LOGICAL REQU0,RNON0
      EXTERNAL FCN
C
C*** End of declarations rewritten by SPAG
C
      DATA ONE,P1,P5,P25,P75,P0001,ZERO/1.0D0,1.0D-1,5.0D-1,2.5D-1,
     &     7.5D-1,1.0D-4,0.0D0/
C
C     epsmch is the machine precision.
C
      EPSMCH = DPMPAR(1)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
C
C     check the input parameters for errors.
C
      IF ( N.GT.0 .AND. M.GE.N .AND. LDFJAC.GE.M .AND. 
     &     FTOL.GE.ZERO .AND. XTOL.GE.ZERO .AND. GTOL.GE.ZERO .AND. 
     &     MAXFEV.GT.0 .AND. FACTOR.GT.ZERO ) THEN
         IF ( MODE.EQ.2 ) THEN
            DO J = 1,N
               IF ( DIAG(J).LE.ZERO ) GOTO 100
            END DO
         END IF
C
C     evaluate the function at the starting point
C     and calculate its norm.
C
         IFLAG = 1
         CALL FCN(M,N,X,FVEC,IFLAG)
         NFEV = 1
         IF ( IFLAG.GE.0 ) THEN
            FNORM = DNRM2(M,FVEC,1)
C
C     initialize levenberg-marquardt parameter and iteration counter.
C
            PAR = ZERO
            ITER = 1
C
C     beginning of the outer loop.
C
C
C        calculate the jacobian matrix.
C
 20         CONTINUE
            IFLAG = 2
            CALL FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
            NFEV = NFEV + N
            IF ( IFLAG.GE.0 ) THEN
C
C        if requested, call fcn to enable printing of iterates.
C
               IF ( NPRINT.GT.0 ) THEN
                  IFLAG = 0
                  IF ( MOD(ITER-1,NPRINT).EQ.0 )
     &                 CALL FCN(M,N,X,FVEC,IFLAG)
                  IF ( IFLAG.LT.0 ) GOTO 100
               END IF
C
C        compute the qr factorization of the jacobian.
C
               CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
C
C        on the first iteration and if mode is 1, scale according
C        to the norms of the columns of the initial jacobian.
C
               IF ( ITER.EQ.1 ) THEN
                  IF ( MODE.NE.2 ) THEN
                     DO J = 1,N
                        DIAG(J) = WA2(J)
C                       IF ( WA2(J).EQ.ZERO ) DIAG(J) = ONE
                        IF ( REQU0(WA2(J)) ) DIAG(J) = ONE
                     END DO
                  END IF
C
C        on the first iteration, calculate the norm of the scaled x
C        and initialize the step bound delta.
C
                  DO J = 1,N
                     WA3(J) = DIAG(J)*X(J)
                  END DO
                  XNORM = DNRM2(N,WA3,1)
                  DELTA = FACTOR*XNORM
C                 IF ( DELTA.EQ.ZERO ) DELTA = FACTOR
                  IF ( REQU0(DELTA) ) DELTA = FACTOR
               END IF
C
C        form (q transpose)*fvec and store the first n components in
C        qtf.
C
               DO I = 1,M
                  WA4(I) = FVEC(I)
               END DO
               DO J = 1,N
C                 IF ( FJAC(J,J).NE.ZERO ) THEN
                  IF ( RNON0(FJAC(J,J)) ) THEN
                     RSUM = ZERO
                     DO I = J,M
                        RSUM = RSUM + FJAC(I,J)*WA4(I)
                     END DO
                     TEMP = -RSUM/FJAC(J,J)
                     DO I = J,M
                        WA4(I) = WA4(I) + FJAC(I,J)*TEMP
                     END DO
                  END IF
                  FJAC(J,J) = WA1(J)
                  QTF(J) = WA4(J)
               END DO
C
C        compute the norm of the scaled gradient.
C
               GNORM = ZERO
C              IF ( FNORM.NE.ZERO ) THEN
               IF ( RNON0(FNORM) ) THEN
                  DO J = 1,N
                     L = IPVT(J)
C                    IF ( WA2(L).NE.ZERO ) THEN
                     IF ( RNON0(WA2(L)) ) THEN
                        RSUM = ZERO
                        DO I = 1,J
                           RSUM = RSUM + FJAC(I,J)*(QTF(I)/FNORM)
                        END DO
                        GNORM = DMAX1(GNORM,DABS(RSUM/WA2(L)))
                     END IF
                  END DO
               END IF
C
C        test for convergence of the gradient norm.
C
               IF ( GNORM.LE.GTOL ) INFO = 4
               IF ( INFO.EQ.0 ) THEN
C
C        rescale if necessary.
C
                  IF ( MODE.NE.2 ) THEN
                     DO J = 1,N
                        DIAG(J) = DMAX1(DIAG(J),WA2(J))
                     END DO
                  END IF
C
C        beginning of the inner loop.
C
C           determine the levenberg-marquardt parameter.
C
 25               CONTINUE
                  CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,
     &                       WA2,WA3,WA4)
C
C           store the direction p and x + p. calculate the norm of p.
C
                  DO J = 1,N
                     WA1(J) = -WA1(J)
                     WA2(J) = X(J) + WA1(J)
                     WA3(J) = DIAG(J)*WA1(J)
                  END DO
                  PNORM = DNRM2(N,WA3,1)
C
C           on the first iteration, adjust the initial step bound.
C
                  IF ( ITER.EQ.1 ) DELTA = DMIN1(DELTA,PNORM)
C
C           evaluate the function at x + p and calculate its norm.
C
                  IFLAG = 1
                  CALL FCN(M,N,WA2,WA4,IFLAG)
                  NFEV = NFEV + 1
                  IF ( IFLAG.GE.0 ) THEN
                     FNORM1 = DNRM2(M,WA4,1)
C
C           compute the scaled actual reduction.
C
                     ACTRED = -ONE
                     IF ( P1*FNORM1.LT.FNORM ) ACTRED = ONE - 
     &                    (FNORM1/FNORM)**2
C
C           compute the scaled predicted reduction and
C           the scaled directional derivative.
C
                     DO J = 1,N
                        WA3(J) = ZERO
                        L = IPVT(J)
                        TEMP = WA1(L)
                        DO I = 1,J
                           WA3(I) = WA3(I) + FJAC(I,J)*TEMP
                        END DO
                     END DO
                     TEMP1 = DNRM2(N,WA3,1)/FNORM
                     TEMP2 = (DSQRT(PAR)*PNORM)/FNORM
                     PRERED = TEMP1**2 + TEMP2**2/P5
                     DIRDER = -(TEMP1**2+TEMP2**2)
C
C           compute the ratio of the actual to the predicted
C           reduction.
C
                     RATIO = ZERO
C                    IF ( PRERED.NE.ZERO ) RATIO = ACTRED/PRERED
                     IF ( RNON0(PRERED) ) RATIO = ACTRED/PRERED
C
C           update the step bound.
C
                     IF ( RATIO.LE.P25 ) THEN
                        IF ( ACTRED.GE.ZERO ) TEMP = P5
                        IF ( ACTRED.LT.ZERO )
     &                       TEMP = P5*DIRDER/(DIRDER+P5*ACTRED)
                        IF ( P1*FNORM1.GE.FNORM .OR. TEMP.LT.P1 )
     &                       TEMP = P1
                        DELTA = TEMP*DMIN1(DELTA,PNORM/P1)
                        PAR = PAR/TEMP
C                    ELSE IF ( PAR.EQ.ZERO .OR. RATIO.GE.P75 ) THEN
                     ELSE IF ( REQU0(PAR) .OR. RATIO.GE.P75 ) THEN
                        DELTA = PNORM/P5
                        PAR = P5*PAR
                     END IF
C
C           test for successful iteration.
C
                     IF ( RATIO.GE.P0001 ) THEN
C
C           successful iteration. update x, fvec, and their norms.
C
                        DO J = 1,N
                           X(J) = WA2(J)
                           WA2(J) = DIAG(J)*X(J)
                        END DO
                        DO I = 1,M
                           FVEC(I) = WA4(I)
                        END DO
                        XNORM = DNRM2(N,WA2,1)
                        FNORM = FNORM1
                        ITER = ITER + 1
                     END IF
C
C           tests for convergence.
C
                     IF ( DABS(ACTRED).LE.FTOL .AND. 
     &                    PRERED.LE.FTOL .AND. P5*RATIO.LE.ONE )
     &                    INFO = 1
                     IF ( DELTA.LE.XTOL*XNORM ) INFO = 2
                     IF ( DABS(ACTRED).LE.FTOL .AND. 
     &                    PRERED.LE.FTOL .AND. P5*RATIO.LE.ONE .AND. 
     &                    INFO.EQ.2 ) INFO = 3
                     IF ( INFO.EQ.0 ) THEN
C
C           tests for termination and stringent tolerances.
C
                        IF ( NFEV.GE.MAXFEV ) INFO = 5
                        IF ( DABS(ACTRED).LE.EPSMCH .AND. 
     &                       PRERED.LE.EPSMCH .AND. P5*RATIO.LE.ONE )
     &                       INFO = 6
                        IF ( DELTA.LE.EPSMCH*XNORM ) INFO = 7
                        IF ( GNORM.LE.EPSMCH ) INFO = 8
                        IF ( INFO.EQ.0 ) THEN
C
C           end of the inner loop. repeat if iteration unsuccessful.
C
C
C        end of the outer loop.
C
                           IF ( RATIO.LT.P0001 ) GOTO 25
                           GOTO 20
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     termination, either normal or user imposed.
C
 100  CONTINUE
      IF ( IFLAG.LT.0 ) INFO = IFLAG
      IFLAG = 0
      IF ( NPRINT.GT.0 ) CALL FCN(M,N,X,FVEC,IFLAG)
C
      END
C*==dpmpar.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      REAL*8 FUNCTION DPMPAR(I)
C   ********************************************************************
C
C     This function provides double precision machine parameters
C     when the appropriate set of data statements is activated (by
C     removing the c from column 1) and all other data statements are
C     rendered inactive. Most of the parameter values were obtained
C     from the corresponding Bell Laboratories Port Library function.
C
C     The function statement is
C
C       double precision function dpmpar(i)
C
C     where
C
C       i is an integer input variable set to 1, 2, or 3 which
C         selects the desired machine parameter. If the machine has
C         t base b digits and its smallest and largest exponents are
C         emin and emax, respectively, then these parameters are
C
C         dpmpar(1) = b**(1 - t), the machine precision,
C
C         dpmpar(2) = b**(emin - 1), the smallest magnitude,
C
C         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
C
C     Argonne National Laboratory. MINPACK Project. November 1996.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
C
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I
C
C Local variables
C
      REAL*8 DMACH(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA DMACH(1)/2.22044604926D-16/
      DATA DMACH(2)/2.22507385852D-308/
      DATA DMACH(3)/1.79769313485D+308/
C
      DPMPAR = DMACH(I)
C
      END
C*==fdjac2.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
C   ********************************************************************
C
C     this subroutine computes a forward-difference approximation
C     to the m by n jacobian matrix associated with a specified
C     problem of m functions in n variables.
C
C     the subroutine statement is
C
C       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
C
C     where
C
C       fcn is the name of the user-supplied subroutine which
C         calculates the functions. fcn must be declared
C         in an external statement in the user calling
C         program, and should be written as follows.
C
C         subroutine fcn(m,n,x,fvec,iflag)
C         integer m,n,iflag
C         double precision x(n),fvec(m)
C         ----------
C         calculate the functions at x and
C         return this vector in fvec.
C         ----------
C         return
C         end
C
C         the value of iflag should not be changed by fcn unless
C         the user wants to terminate execution of fdjac2.
C         in this case set iflag to a negative integer.
C
C       m is a positive integer input variable set to the number
C         of functions.
C
C       n is a positive integer input variable set to the number
C         of variables. n must not exceed m.
C
C       x is an input array of length n.
C
C       fvec is an input array of length m which must contain the
C         functions evaluated at x.
C
C       fjac is an output m by n array which contains the
C         approximation to the jacobian matrix evaluated at x.
C
C       ldfjac is a positive integer input variable not less than m
C         which specifies the leading dimension of the array fjac.
C
C       iflag is an integer variable which can be used to terminate
C         the execution of fdjac2. see description of fcn.
C
C       epsfcn is an input variable used in determining a suitable
C         step length for the forward-difference approximation. this
C         approximation assumes that the relative errors in the
C         functions are of the order of epsfcn. if epsfcn is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
C
C       wa is a work array of length m.
C
C     subprograms called
C
C       user-supplied ...... fcn
C
C       minpack-supplied ... dpmpar
C
C       fortran-supplied ... dabs,dmax1,dsqrt
C
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
C
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPSFCN
      INTEGER IFLAG,LDFJAC,M,N
      REAL*8 FJAC(LDFJAC,N),FVEC(M),WA(M),X(N)
C
C Local variables
C
      REAL*8 DPMPAR
      REAL*8 EPS,EPSMCH,H,TEMP
      INTEGER I,J
      LOGICAL REQU0
C
C*** End of declarations rewritten by SPAG
C
C      DATA ZERO/0.0D0/
C
C     epsmch is the machine precision.
C
      EPSMCH = DPMPAR(1)
C
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
      DO J = 1,N
         TEMP = X(J)
         H = EPS*DABS(TEMP)
C        IF ( H.EQ.ZERO ) H = EPS
         IF ( REQU0(H) ) H = EPS
         X(J) = TEMP + H
         CALL FCN(M,N,X,WA,IFLAG)
         IF ( IFLAG.LT.0 ) EXIT
         X(J) = TEMP
         DO I = 1,M
            FJAC(I,J) = (WA(I)-FVEC(I))/H
         END DO
      END DO
C
      END
C*==lmpar.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,WA1,WA2)
C   ********************************************************************
C
C     given an m by n matrix a, an n by n nonsingular diagonal
C     matrix d, an m-vector b, and a positive number delta,
C     the problem is to determine a value for the parameter
C     par such that if x solves the system
C
C           a*x = b ,     sqrt(par)*d*x = 0 ,
C
C     in the least squares sense, and dxnorm is the euclidean
C     norm of d*x, then either par is zero and
C
C           (dxnorm-delta) .le. 0.1*delta ,
C
C     or par is positive and
C
C           abs(dxnorm-delta) .le. 0.1*delta .
C
C     this subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     qr factorization, with column pivoting, of a. that is, if
C     a*p = q*r, where p is a permutation matrix, q has orthogonal
C     columns, and r is an upper triangular matrix with diagonal
C     elements of nonincreasing magnitude, then lmpar expects
C     the full upper triangle of r, the permutation matrix p,
C     and the first n components of (q transpose)*b. on output
C     lmpar also provides an upper triangular matrix s such that
C
C            t   t                   t
C           p *(a *a + par*d*d)*p = s *s .
C
C     s is employed within lmpar and may be of separate interest.
C
C     only a few iterations are generally needed for convergence
C     of the algorithm. if, however, the limit of 10 iterations
C     is reached, then the output par will contain the best
C     value obtained so far.
C
C     the subroutine statement is
C
C       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
C                        wa1,wa2)
C
C     where
C
C       n is a positive integer input variable set to the order of r.
C
C       r is an n by n array. on input the full upper triangle
C         must contain the full upper triangle of the matrix r.
C         on output the full upper triangle is unaltered, and the
C         strict lower triangle contains the strict upper triangle
C         (transposed) of the upper triangular matrix s.
C
C       ldr is a positive integer input variable not less than n
C         which specifies the leading dimension of the array r.
C
C       ipvt is an integer input array of length n which defines the
C         permutation matrix p such that a*p = q*r. column j of p
C         is column ipvt(j) of the identity matrix.
C
C       diag is an input array of length n which must contain the
C         diagonal elements of the matrix d.
C
C       qtb is an input array of length n which must contain the first
C         n elements of the vector (q transpose)*b.
C
C       delta is a positive input variable which specifies an upper
C         bound on the euclidean norm of d*x.
C
C       par is a nonnegative variable. on input par contains an
C         initial estimate of the levenberg-marquardt parameter.
C         on output par contains the final estimate.
C
C       x is an output array of length n which contains the least
C         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
C         for the output par.
C
C       sdiag is an output array of length n which contains the
C         diagonal elements of the upper triangular matrix s.
C
C       wa1 and wa2 are work arrays of length n.
C
C     subprograms called
C
C       minpack-supplied ... dpmpar,qrsolv
C
C       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
C
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
C
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DELTA,PAR
      INTEGER LDR,N
      REAL*8 DIAG(N),QTB(N),R(LDR,N),SDIAG(N),WA1(N),WA2(N),X(N)
      INTEGER IPVT(N)
C
C Local variables
C
      REAL*8 DNRM2,DPMPAR
      REAL*8 DWARF,DXNORM,FP,GNORM,P001,P1,PARC,PARL,PARU,RSUM,TEMP,ZERO
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      LOGICAL REQU0
C
C*** End of declarations rewritten by SPAG
C
      DATA P1,P001,ZERO/1.0D-1,1.0D-3,0.0D0/
C
C     dwarf is the smallest positive magnitude.
C
      DWARF = DPMPAR(2)
C
C     compute and store in x the gauss-newton direction. if the
C     jacobian is rank-deficient, obtain a least squares solution.
C
      NSING = N
      DO J = 1,N
         WA1(J) = QTB(J)
C        IF ( R(J,J).EQ.ZERO .AND. NSING.EQ.N ) NSING = J - 1
         IF ( REQU0(R(J,J)) .AND. NSING.EQ.N ) NSING = J - 1
         IF ( NSING.LT.N ) WA1(J) = ZERO
      END DO
      IF ( NSING.GE.1 ) THEN
         DO K = 1,NSING
            J = NSING - K + 1
            WA1(J) = WA1(J)/R(J,J)
            TEMP = WA1(J)
            JM1 = J - 1
            IF ( JM1.GE.1 ) THEN
               DO I = 1,JM1
                  WA1(I) = WA1(I) - R(I,J)*TEMP
               END DO
            END IF
         END DO
      END IF
      DO J = 1,N
         L = IPVT(J)
         X(L) = WA1(J)
      END DO
C
C     initialize the iteration counter.
C     evaluate the function at the origin, and test
C     for acceptance of the gauss-newton direction.
C
      ITER = 0
      DO J = 1,N
         WA2(J) = DIAG(J)*X(J)
      END DO
      DXNORM = DNRM2(N,WA2,1)
      FP = DXNORM - DELTA
      IF ( FP.LE.P1*DELTA ) THEN
C
C     termination.
C
         IF ( ITER.EQ.0 ) PAR = ZERO
      ELSE
C
C     if the jacobian is not rank deficient, the newton
C     step provides a lower bound, parl, for the zero of
C     the function. otherwise set this bound to zero.
C
         PARL = ZERO
         IF ( NSING.GE.N ) THEN
            DO J = 1,N
               L = IPVT(J)
               WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
            END DO
            DO J = 1,N
               RSUM = ZERO
               JM1 = J - 1
               IF ( JM1.GE.1 ) THEN
                  DO I = 1,JM1
                     RSUM = RSUM + R(I,J)*WA1(I)
                  END DO
               END IF
               WA1(J) = (WA1(J)-RSUM)/R(J,J)
            END DO
            TEMP = DNRM2(N,WA1,1)
            PARL = ((FP/DELTA)/TEMP)/TEMP
         END IF
C
C     calculate an upper bound, paru, for the zero of the function.
C
         DO J = 1,N
            RSUM = ZERO
            DO I = 1,J
               RSUM = RSUM + R(I,J)*QTB(I)
            END DO
            L = IPVT(J)
            WA1(J) = RSUM/DIAG(L)
         END DO
         GNORM = DNRM2(N,WA1,1)
         PARU = GNORM/DELTA
C        IF ( PARU.EQ.ZERO ) PARU = DWARF/DMIN1(DELTA,P1)
         IF ( REQU0(PARU) ) PARU = DWARF/DMIN1(DELTA,P1)
C
C     if the input par lies outside of the interval (parl,paru),
C     set par to the closer endpoint.
C
         PAR = DMAX1(PAR,PARL)
         PAR = DMIN1(PAR,PARU)
C        IF ( PAR.EQ.ZERO ) PAR = GNORM/DXNORM
         IF ( REQU0(PAR) ) PAR = GNORM/DXNORM
C
C     beginning of an iteration.
C
 50      CONTINUE
         ITER = ITER + 1
C
C        evaluate the function at the current value of par.
C
C        IF ( PAR.EQ.ZERO ) PAR = DMAX1(DWARF,P001*PARU)
         IF ( REQU0(PAR) ) PAR = DMAX1(DWARF,P001*PARU)
         TEMP = DSQRT(PAR)
         DO J = 1,N
            WA1(J) = TEMP*DIAG(J)
         END DO
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SDIAG,WA2)
         DO J = 1,N
            WA2(J) = DIAG(J)*X(J)
         END DO
         DXNORM = DNRM2(N,WA2,1)
         TEMP = FP
         FP = DXNORM - DELTA
C
C        if the function is small enough, accept the current value
C        of par. also test for the exceptional cases where parl
C        is zero or the number of iterations has reached 10.
C
C        IF ( DABS(FP).LE.P1*DELTA .OR. PARL.EQ.ZERO .AND.
         IF ( DABS(FP).LE.P1*DELTA .OR. REQU0(PARL) .AND. 
     &        FP.LE.TEMP .AND. TEMP.LT.ZERO .OR. ITER.EQ.10 ) THEN
            IF ( ITER.EQ.0 ) PAR = ZERO
         ELSE
C
C        compute the newton correction.
C
            DO J = 1,N
               L = IPVT(J)
               WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
            END DO
            DO J = 1,N
               WA1(J) = WA1(J)/SDIAG(J)
               TEMP = WA1(J)
               JP1 = J + 1
               IF ( N.GE.JP1 ) THEN
                  DO I = JP1,N
                     WA1(I) = WA1(I) - R(I,J)*TEMP
                  END DO
               END IF
            END DO
            TEMP = DNRM2(N,WA1,1)
            PARC = ((FP/DELTA)/TEMP)/TEMP
C
C        depending on the sign of the function, update parl or paru.
C
            IF ( FP.GT.ZERO ) PARL = DMAX1(PARL,PAR)
            IF ( FP.LT.ZERO ) PARU = DMIN1(PARU,PAR)
C
C        compute an improved estimate for par.
C
            PAR = DMAX1(PARL,PAR+PARC)
C
C        end of an iteration.
C
            GOTO 50
         END IF
      END IF
C
C     last card of subroutine lmpar.
C
      END
C*==qrfac.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
C   ********************************************************************
C
C     this subroutine uses householder transformations with column
C     pivoting (optional) to compute a qr factorization of the
C     m by n matrix a. that is, qrfac determines an orthogonal
C     matrix q, a permutation matrix p, and an upper trapezoidal
C     matrix r with diagonal elements of nonincreasing magnitude,
C     such that a*p = q*r. the householder transformation for
C     column k, k = 1,2,...,min(m,n), is of the form
C
C                           t
C           i - (1/u(k))*u*u
C
C     where u has zeros in the first k-1 positions. the form of
C     this transformation and the method of pivoting first
C     appeared in the corresponding linpack subroutine.
C
C     the subroutine statement is
C
C       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
C
C     where
C
C       m is a positive integer input variable set to the number
C         of rows of a.
C
C       n is a positive integer input variable set to the number
C         of columns of a.
C
C       a is an m by n array. on input a contains the matrix for
C         which the qr factorization is to be computed. on output
C         the strict upper trapezoidal part of a contains the strict
C         upper trapezoidal part of r, and the lower trapezoidal
C         part of a contains a factored form of q (the non-trivial
C         elements of the u vectors described above).
C
C       lda is a positive integer input variable not less than m
C         which specifies the leading dimension of the array a.
C
C       pivot is a logical input variable. if pivot is set true,
C         then column pivoting is enforced. if pivot is set false,
C         then no column pivoting is done.
C
C       ipvt is an integer output array of length lipvt. ipvt
C         defines the permutation matrix p such that a*p = q*r.
C         column j of p is column ipvt(j) of the identity matrix.
C         if pivot is false, ipvt is not referenced.
C
C       lipvt is a positive integer input variable. if pivot is false,
C         then lipvt may be as small as 1. if pivot is true, then
C         lipvt must be at least n.
C
C       rdiag is an output array of length n which contains the
C         diagonal elements of r.
C
C       acnorm is an output array of length n which contains the
C         norms of the corresponding columns of the input matrix a.
C         if this information is not needed, then acnorm can coincide
C         with rdiag.
C
C       wa is a work array of length n. if pivot is false, then wa
C         can coincide with rdiag.
C
C     subprograms called
C
C       minpack-supplied ... dpmpar
C
C       fortran-supplied ... dmax1,dsqrt,min0
C
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
C
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LDA,LIPVT,M,N
      LOGICAL PIVOT
      REAL*8 A(LDA,N),ACNORM(N),RDIAG(N),WA(N)
      INTEGER IPVT(LIPVT)
C
C Local variables
C
      REAL*8 AJNORM,EPSMCH,ONE,P05,RSUM,TEMP,ZERO
      REAL*8 DNRM2,DPMPAR
      INTEGER I,J,JP1,K,KMAX,MINMN
      LOGICAL REQU0,RNON0
C
C*** End of declarations rewritten by SPAG
C
      DATA ONE,P05,ZERO/1.0D0,5.0D-2,0.0D0/
C
C     epsmch is the machine precision.
C
      EPSMCH = DPMPAR(1)
C
C     compute the initial column norms and initialize several arrays.
C
      DO J = 1,N
         ACNORM(J) = DNRM2(M,A(1,J),1)
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF ( PIVOT ) IPVT(J) = J
      END DO
C
C     reduce a to r with householder transformations.
C
      MINMN = MIN0(M,N)
      DO J = 1,MINMN
         IF ( PIVOT ) THEN
C
C        bring the column of largest norm into the pivot position.
C
            KMAX = J
            DO K = J,N
               IF ( RDIAG(K).GT.RDIAG(KMAX) ) KMAX = K
            END DO
            IF ( KMAX.NE.J ) THEN
               DO I = 1,M
                  TEMP = A(I,J)
                  A(I,J) = A(I,KMAX)
                  A(I,KMAX) = TEMP
               END DO
               RDIAG(KMAX) = RDIAG(J)
               WA(KMAX) = WA(J)
               K = IPVT(J)
               IPVT(J) = IPVT(KMAX)
               IPVT(KMAX) = K
            END IF
         END IF
C
C        compute the householder transformation to reduce the
C        j-th column of a to a multiple of the j-th unit vector.
C
         AJNORM = DNRM2(M-J+1,A(J,J),1)
C        IF ( AJNORM.NE.ZERO ) THEN
         IF ( RNON0(AJNORM) ) THEN
            IF ( A(J,J).LT.ZERO ) AJNORM = -AJNORM
            DO I = J,M
               A(I,J) = A(I,J)/AJNORM
            END DO
            A(J,J) = A(J,J) + ONE
C
C        apply the transformation to the remaining columns
C        and update the norms.
C
            JP1 = J + 1
            IF ( N.GE.JP1 ) THEN
               DO K = JP1,N
                  RSUM = ZERO
                  DO I = J,M
                     RSUM = RSUM + A(I,J)*A(I,K)
                  END DO
                  TEMP = RSUM/A(J,J)
                  DO I = J,M
                     A(I,K) = A(I,K) - TEMP*A(I,J)
                  END DO
C                 IF ( .NOT.(.NOT.PIVOT .OR. RDIAG(K).EQ.ZERO) ) THEN
                  IF ( .NOT.(.NOT.PIVOT .OR. REQU0(RDIAG(K))) ) THEN
                     TEMP = A(J,K)/RDIAG(K)
                     RDIAG(K) = RDIAG(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
                     IF ( P05*(RDIAG(K)/WA(K))**2.LE.EPSMCH ) THEN
                        RDIAG(K) = DNRM2(M-J,A(JP1,K),1)
                        WA(K) = RDIAG(K)
                     END IF
                  END IF
               END DO
            END IF
         END IF
         RDIAG(J) = -AJNORM
      END DO
C
C     last card of subroutine qrfac.
C
      END
C*==qrsolv.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)
C   ********************************************************************
C
C     given an m by n matrix a, an n by n diagonal matrix d,
C     and an m-vector b, the problem is to determine an x which
C     solves the system
C
C           a*x = b ,     d*x = 0 ,
C
C     in the least squares sense.
C
C     this subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     qr factorization, with column pivoting, of a. that is, if
C     a*p = q*r, where p is a permutation matrix, q has orthogonal
C     columns, and r is an upper triangular matrix with diagonal
C     elements of nonincreasing magnitude, then qrsolv expects
C     the full upper triangle of r, the permutation matrix p,
C     and the first n components of (q transpose)*b. the system
C     a*x = b, d*x = 0, is then equivalent to
C
C                  t       t
C           r*z = q *b ,  p *d*p*z = 0 ,
C
C     where x = p*z. if this system does not have full rank,
C     then a least squares solution is obtained. on output qrsolv
C     also provides an upper triangular matrix s such that
C
C            t   t               t
C           p *(a *a + d*d)*p = s *s .
C
C     s is computed within qrsolv and may be of separate interest.
C
C     the subroutine statement is
C
C       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
C
C     where
C
C       n is a positive integer input variable set to the order of r.
C
C       r is an n by n array. on input the full upper triangle
C         must contain the full upper triangle of the matrix r.
C         on output the full upper triangle is unaltered, and the
C         strict lower triangle contains the strict upper triangle
C         (transposed) of the upper triangular matrix s.
C
C       ldr is a positive integer input variable not less than n
C         which specifies the leading dimension of the array r.
C
C       ipvt is an integer input array of length n which defines the
C         permutation matrix p such that a*p = q*r. column j of p
C         is column ipvt(j) of the identity matrix.
C
C       diag is an input array of length n which must contain the
C         diagonal elements of the matrix d.
C
C       qtb is an input array of length n which must contain the first
C         n elements of the vector (q transpose)*b.
C
C       x is an output array of length n which contains the least
C         squares solution of the system a*x = b, d*x = 0.
C
C       sdiag is an output array of length n which contains the
C         diagonal elements of the upper triangular matrix s.
C
C       wa is a work array of length n.
C
C     subprograms called
C
C       fortran-supplied ... dabs,dsqrt
C
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
C
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LDR,N
      REAL*8 DIAG(N),QTB(N),R(LDR,N),SDIAG(N),WA(N),X(N)
      INTEGER IPVT(N)
C
C Local variables
C
      REAL*8 COTAN,P25,P5,QTBPJ,RCOS,RSIN,RSUM,RTAN,TEMP,ZERO
      INTEGER I,J,JP1,K,KP1,L,NSING
      LOGICAL REQU0,RNON0
C
C*** End of declarations rewritten by SPAG
C
      DATA P5,P25,ZERO/5.0D-1,2.5D-1,0.0D0/
C
C     copy r and (q transpose)*b to preserve input and initialize s.
C     in particular, save the diagonal elements of r in x.
C
      DO J = 1,N
         DO I = J,N
            R(I,J) = R(J,I)
         END DO
         X(J) = R(J,J)
         WA(J) = QTB(J)
      END DO
C
C     eliminate the diagonal matrix d using a givens rotation.
C
      DO J = 1,N
C
C        prepare the row of d to be eliminated, locating the
C        diagonal element using p from the qr factorization.
C
         L = IPVT(J)
C        IF ( DIAG(L).NE.ZERO ) THEN
         IF ( RNON0(DIAG(L)) ) THEN
            DO K = J,N
               SDIAG(K) = ZERO
            END DO
            SDIAG(J) = DIAG(L)
C
C        the transformations to eliminate the row of d
C        modify only a single element of (q transpose)*b
C        beyond the first n, which is initially zero.
C
            QTBPJ = ZERO
            DO K = J,N
C
C           determine a givens rotation which eliminates the
C           appropriate element in the current row of d.
C
C              IF ( SDIAG(K).NE.ZERO ) THEN
               IF ( RNON0(SDIAG(K)) ) THEN
                  IF ( DABS(R(K,K)).GE.DABS(SDIAG(K)) ) THEN
                     RTAN = SDIAG(K)/R(K,K)
                     RCOS = P5/DSQRT(P25+P25*RTAN**2)
                     RSIN = RCOS*RTAN
                  ELSE
                     COTAN = R(K,K)/SDIAG(K)
                     RSIN = P5/DSQRT(P25+P25*COTAN**2)
                     RCOS = RSIN*COTAN
                  END IF
C
C           compute the modified diagonal element of r and
C           the modified element of ((q transpose)*b,0).
C
                  R(K,K) = RCOS*R(K,K) + RSIN*SDIAG(K)
                  TEMP = RCOS*WA(K) + RSIN*QTBPJ
                  QTBPJ = -RSIN*WA(K) + RCOS*QTBPJ
                  WA(K) = TEMP
C
C           accumulate the tranformation in the row of s.
C
                  KP1 = K + 1
                  IF ( N.GE.KP1 ) THEN
                     DO I = KP1,N
                        TEMP = RCOS*R(I,K) + RSIN*SDIAG(I)
                        SDIAG(I) = -RSIN*R(I,K) + RCOS*SDIAG(I)
                        R(I,K) = TEMP
                     END DO
                  END IF
               END IF
            END DO
         END IF
C
C        store the diagonal element of s and restore
C        the corresponding diagonal element of r.
C
         SDIAG(J) = R(J,J)
         R(J,J) = X(J)
      END DO
C
C     solve the triangular system for z. if the system is
C     singular, then obtain a least squares solution.
C
      NSING = N
      DO J = 1,N
C        IF ( SDIAG(J).EQ.ZERO .AND. NSING.EQ.N ) NSING = J - 1
         IF ( REQU0(SDIAG(J)) .AND. NSING.EQ.N ) NSING = J - 1
         IF ( NSING.LT.N ) WA(J) = ZERO
      END DO
      IF ( NSING.GE.1 ) THEN
         DO K = 1,NSING
            J = NSING - K + 1
            RSUM = ZERO
            JP1 = J + 1
            IF ( NSING.GE.JP1 ) THEN
               DO I = JP1,NSING
                  RSUM = RSUM + R(I,J)*WA(I)
               END DO
            END IF
            WA(J) = (WA(J)-RSUM)/SDIAG(J)
         END DO
      END IF
C
C     permute the components of z back to components of x.
C
      DO J = 1,N
         L = IPVT(J)
         X(L) = WA(J)
      END DO
C
      END
