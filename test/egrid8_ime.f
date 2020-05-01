C*==egrid8_ime.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE EGRID8_IME(EPS,EMIN,EMAX,NE,SETEPS)
C   ********************************************************************
C        Calculated EPS for energy grid 8 to fit requested EIMAG
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FOUND_SECTION,FOUND_REAL
      IMPLICIT NONE
C*--EGRID8_IME9
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EGRID8_IME')
C
C Dummy arguments
C
      REAL*8 EMAX,EMIN,EPS
      INTEGER NE
      LOGICAL SETEPS
C
C Local variables
C
      REAL*8 EIMAG,EPSMAX,EPSMIN,TOL
      INTEGER IFAIL
      REAL*8 IMEEPSDIF,ZBRENT
      EXTERNAL IMEEPSDIF
C
C*** End of declarations rewritten by SPAG
C
      SETEPS = .FALSE.
C
      CALL INPUT_FIND_SECTION('ENERGY',0)
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_REAL('IME',EIMAG,9999D0,0)
C
      IF ( FOUND_REAL .AND. EIMAG.GT.1.0D-9 ) THEN
         EPSMIN = 1.0D-12
         EPSMAX = 1.0D+12
         TOL = MAX(ABS(EIMAG)*0.01D0,1.0D-11)
         EPS = ZBRENT(IMEEPSDIF,EPSMIN,EPSMAX,TOL,EMIN,EMAX,EIMAG,NE,
     &         IFAIL)
C
         IF ( IFAIL.NE.0 ) THEN
            WRITE (6,'(A,/,3X,A,A,F15.8)')
     &              'Setting EPS for GRID=8 from ImE failed !',
     &             'Re-try it with a more reasonable (smaller)',
     &             ' value of ImE than',EIMAG
            CALL STOP_MESSAGE(ROUTINE,'EPS from ImE failed')
         END IF
C
         SETEPS = .TRUE.
      END IF
C
      END
C*==imeepsdif.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      FUNCTION IMEEPSDIF(EPS,EMIN,EMAX,EIMAG,NE)
C **********************************************************
C Different between imaginary part of ETAB for given EPS
C   and requested EIMAG
C **********************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,PI
      IMPLICIT NONE
C*--IMEEPSDIF81
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='IMEEPSDIF')
C
C Dummy arguments
C
      REAL*8 EIMAG,EMAX,EMIN,EPS
      INTEGER NE
      REAL*8 IMEEPSDIF
C
C Local variables
C
      COMPLEX*16 ETABNE
      INTEGER IA_ERR
      REAL*8 PHI1,PHI2,PHINE,R,WDUM(:),Y1,Y2,YYNE,Z(:),Z0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WDUM,Z
C
      ALLOCATE (WDUM(NE),Z(NE),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WDUM, Z')
C
      CALL GAULEG(-1.0D0,1.0D0,Z,WDUM,NE)
C
      PHI1 = PI
      PHI2 = 0.0D0
C
      R = (EMAX-EMIN)/2.0D0
      Z0 = EMIN + R
C
      Y1 = -LOG((EPS+PHI1)/EPS)
      Y2 = -LOG((EPS+PHI2)/EPS)
C
      YYNE = 0.5D0*(Y2-Y1)*Z(NE) + 0.5D0*(Y2+Y1)
      PHINE = EPS*(EXP(-YYNE)-1.0D0)
      ETABNE = Z0 + R*EXP(CI*PHINE)
C
      IMEEPSDIF = DIMAG(ETABNE) - EIMAG
C
      DEALLOCATE (WDUM,Z,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      END
C*==zbrent.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      FUNCTION ZBRENT(FUNC,X1,X2,TOL,EMIN,EMAX,EIMAG,NE,IFAIL)
C ********************************************************************
C Searches roots of function FUNC in interval X1, X2 with accuracy TOL
C    using the Brent's method.
C    Parameters are maximum number of iteration ITMAX
C    and machine precision  EPS.
C ********************************************************************
C
      IMPLICIT NONE
C*--ZBRENT154
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ITMAX
      PARAMETER (ITMAX=100)
      REAL*8 EPS,ZERO,ZERO_THRESHOLD
      PARAMETER (EPS=3.D-15,ZERO=0.0D0,ZERO_THRESHOLD=1D-20)
C
C Dummy arguments
C
      REAL*8 EIMAG,EMAX,EMIN,TOL,X1,X2
      INTEGER IFAIL,NE
      REAL*8 FUNC
      REAL*8 ZBRENT
C
C Local variables
C
      REAL*8 A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      INTEGER ITER
      EXTERNAL FUNC
C
C*** End of declarations rewritten by SPAG
C
      IFAIL = 0
C
      A = X1
      B = X2
      FA = FUNC(A,EMIN,EMAX,EIMAG,NE)
      FB = FUNC(B,EMIN,EMAX,EIMAG,NE)
      IF ( FB*FA.GT.ZERO ) THEN
         WRITE (6,*)
         WRITE (6,*) 'ZBPRENT error:  No sign change inside interval'
         WRITE (6,*) '   X1 =',X1,'     FA =',FA
         WRITE (6,*) '   X2 =',X2,'     FB =',FB
         WRITE (6,*)
         IFAIL = 1
         RETURN
      END IF
C
      FC = FB
      DO ITER = 1,ITMAX
         IF ( FB*FC.GT.ZERO ) THEN
            C = A                ! Rename A, B, C and adjust
            FC = FA              !   bounding interval D
            D = B - A
            E = D
         END IF
         IF ( ABS(FC).LT.ABS(FB) ) THEN
            A = B
            B = C
            C = A
            FA = FB
            FB = FC
            FC = FA
         END IF
C
         TOL1 = 2.0D0*EPS*ABS(B) + 0.5D0*TOL      ! Convergence check
         XM = 0.5D0*(C-B)
         IF ( ABS(XM).LE.TOL1 .OR. ABS(FB).LT.ZERO_THRESHOLD ) THEN
            ZBRENT = B
            RETURN
         END IF
C
C Attempt inverse quadratic interpolation
C
         IF ( ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB) ) THEN
            S = FB/FA
            IF ( ABS(A-C).LT.ZERO_THRESHOLD ) THEN
               P = 2.0D0*XM*S
               Q = 1.0D0 - S
            ELSE
               Q = FA/FC
               R = FB/FC
               P = S*(2.0D0*XM*Q*(Q-R)-(B-A)*(R-1.0D0))
               Q = (Q-1.0D0)*(R-1.0D0)*(S-1.0D0)
            END IF
            IF ( P.GT.ZERO ) Q = -Q          ! Check whether in bounds
            P = ABS(P)
            IF ( 2.0D0*P.LT.MIN(3.0D0*XM*Q-ABS(TOL1*Q),ABS(E*Q)) ) THEN
               E = D        ! Accept interpolation
               D = P/Q
            ELSE
               D = XM       ! Interpolation failed, use bisection
               E = D
            END IF
C
C
C Bounds decreasing too slowly, use bisection
C
         ELSE
            D = XM
            E = D
         END IF
         A = B              ! Move last best guess to A.
         FA = FB
         IF ( ABS(D).GT.TOL1 ) THEN  ! Evaluate new trial root
            B = B + D
         ELSE
            B = B + SIGN(TOL1,XM)
         END IF
         FB = FUNC(B,EMIN,EMAX,EIMAG,NE)
      END DO
C
      WRITE (6,*)
      WRITE (6,*) 'ZBRENT error:  Exceeding maximum iterations !'
      WRITE (6,*) '  Possible action:  Increase parameter ITMAX'
      WRITE (6,*)
      IFAIL = 2
C
      END
C
