C*==scfmad2d_smat.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC1,VEC2,RM2,NRMAX,NSHLR,
     &                         NSR,GN2,NGMAX,NSHLG,NSG,SMAT,VOLUC_2D,
     &                         LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C **********************************************************************
C *                                                                    *
C *   calculation of lattice sums for l .le. 2*nlfp :                  *
C *                                                                    *
C *                    ylm( q(i) - q(j) + rm )                         *
C *         sum      ===========================                       *
C *                  | q(i) - q(j) + rm |**(l+1)                       *
C *                                                                    *
C *          - summed over all 2D lattice vectors rm  -                *
C *                                                                    *
C *   ylm       : real spherical harmic to given l,m                   *
C *                                                                    *
C *   The sum is done different in the plane (qi-qj)z = 0              *
C *   and out of the plane. In plane an Ewald procedure similar        *
C *   to the 3d is used and we perform 2 sums (real and reciprocal)    *
C *   the l= 2,4 m=0 terms are calculated with a different method      *
C *                                                                    *
C *   The l=0 term is calculated with a extra factror sqrt(4*pi) this  *
C *   is for transparency reasons (so that the correction terms        *
C *   *r=0,g=0* can be followed in the program)                        *
C *   Literature : lm = (0,0), (1,0) terms :PRB 40, 12164 (1989)       *
C *                                         PRB 49, 2721 (1994)        *
C *                                         PRB 47, 16525 (1993)       *
C *                                         Zimman , p.39-40           *
C *       l=2,4 (m=0) terms are done with recursive diferentiation     *
C *                                                  v. 16.8.99        *
C *       The l multipoles are treated using the expansion             *
C *       for a complex plane wave.                                    *
C *       eq.(33) , M. Weinert, J. Math Phys. 22, 2439 (1981)          *
C *                                                                    *
C *                                                                    *
C **********************************************************************
      USE MOD_CONSTANTS,ONLY:PI,CONST_4PI,CONST_2PI,CI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD2D_SMAT')
      REAL*8 BOUND
      PARAMETER (BOUND=1.0D-9)
C
C Dummy arguments
C
      REAL*8 ALAT,VOLUC_2D
      INTEGER LMAX_RGNTMAD,NGMAX,NLMAD_LOC,NLMMAX_RGNTMAD,NRMAX,NSHLG,
     &        NSHLR
      REAL*8 GN2(2,*),RM2(2,*),SMAT(NLMMAX_RGNTMAD),VEC1(3),VEC2(3)
      INTEGER NSG(*),NSR(*)
C
C Local variables
C
      REAL*8 ALPHA,BETA,CON,CON1,DFAC(0:2*LMAX_RGNTMAD+1),DOT1,DQ1,DQ2,
     &       DQ3,DQDOTG,EXPBSQ,G(0:LMAX_RGNTMAD),G1,G2,G3,GA,
     &       GAL(0:LMAX_RGNTMAD),GI(0:4),GR(0:4),LAMBDA,
     &       PREF0(0:LMAX_RGNTMAD),R,R1,R2,R3,RFAC,S,SIGNRZ,
     &       SIGNRZL(0:LMAX_RGNTMAD),STEST0,YLM(NLMMAX_RGNTMAD),
     &       YLMPREF(0:LMAX_RGNTMAD),
     &       YLMPREF1(0:LMAX_RGNTMAD,0:LMAX_RGNTMAD)
      COMPLEX*16 APREFMM,APREFPP,BFAC,CIM(0:LMAX_RGNTMAD),
     &           EXPONL(0:LMAX_RGNTMAD),FACTEXP,PREF2(LMAX_RGNTMAD),
     &           S0(NLMMAX_RGNTMAD),SIMAG,STEST(NLMMAX_RGNTMAD),
     &           STESTNEW(NLMMAX_RGNTMAD)
      REAL*8 DERFC
      INTEGER I,IM,IR,L,LM,LMAX,LMMAX,M
      EXTERNAL FPLANEG,FPLANER
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C Factorial
C
      DFAC(0) = 1
      DO L = 1,2*LMAX_RGNTMAD + 1
         DFAC(L) = DFAC(L-1)*L
      END DO
      DO L = 0,LMAX_RGNTMAD
         PREF0(L) = 0D0
      END DO
C
      LMAX = 2*(NLMAD_LOC-1)
      LMMAX = (LMAX+1)**2
C
      IF ( NLMAD_LOC.GT.1 ) THEN
         PREF0(2) = SQRT(5D0/PI)/2D0/2D0
         PREF0(4) = 3D0*SQRT(9D0/PI)/16D0/9D0
      END IF
C
      SMAT(1:NLMMAX_RGNTMAD) = 0D0
C
C --> choose proper splitting parameter
C
      LAMBDA = SQRT(PI)/ALAT
C
      DQ1 = (VEC2(1)-VEC1(1))*ALAT   ! SCALE WITH ALAT
      DQ2 = (VEC2(2)-VEC1(2))*ALAT
      DQ3 = (VEC2(3)-VEC1(3))*ALAT
C
C --> Initialise
C
      DO LM = 1,LMMAX
         STEST(LM) = 0D0
         STESTNEW(LM) = 0D0
      END DO
C
C---> Add correction if rz = 0
C
      IF ( DABS(DQ3).LT.1D-6 ) THEN
         STEST(1) = STEST(1) - 2D0*LAMBDA/SQRT(PI) - 2D0*SQRT(PI)
     &              /LAMBDA/VOLUC_2D
         STESTNEW(1) = STESTNEW(1) - 2D0*LAMBDA/SQRT(PI) - 2D0*SQRT(PI)
     &                 /LAMBDA/VOLUC_2D
C
         IF ( (DQ1*DQ1+DQ2*DQ2).GT.1D-6 ) THEN
            STEST(1) = STEST(1) + 2D0*LAMBDA/SQRT(PI)
            STESTNEW(1) = STESTNEW(1) + 2D0*LAMBDA/SQRT(PI)
         END IF
      ELSE
C
C---> Add correction if rz<> 0
C
         STEST(1) = STEST(1) - DABS(DQ3)*CONST_4PI/2D0/VOLUC_2D
         STEST(3) = STEST(3) - DABS(DQ3)/DQ3*SQRT(3D0*CONST_4PI)
     &              /2D0/VOLUC_2D
C  -d/dz .... the correction for higher l vanishes...
      END IF
C
C **********************************************************************
C ******************    I N-P L A N E      M = 0  **********************
C **********************************************************************
      IF ( DABS(DQ3).LT.1D-6 ) THEN
C
C --> Real space sum
C
C ======================================================================
         DO IR = 1,NRMAX
            R1 = DQ1 - RM2(1,IR)
            R2 = DQ2 - RM2(2,IR)
            R3 = 0D0
            R = SQRT(R1*R1+R2*R2)
C-----------------------------------------------------------------------
            IF ( R.GT.1D-8 ) THEN
               ALPHA = LAMBDA*R
               CALL FPLANER(ALPHA,GR,R)
               DO L = 0,4
                  LM = L*(L+1) + 1 ! m =0
                  STEST(LM) = STEST(LM) + GR(L)
               END DO
               CALL CALC_RHPLM(R1/R,R2/R,R3/R,YLM,LMAX_RGNTMAD,
     &                         NLMMAX_RGNTMAD)
               CALL GAMFC(ALPHA,G,LMAX_RGNTMAD,R)
               YLM(1) = 1D0     ! just definition matter
C
               DO L = 0,LMAX
                  RFAC = G(L)/SQRT(PI)
                  DO M = -L,L
                     LM = L*(L+1) + M + 1
                     STESTNEW(LM) = STESTNEW(LM) + YLM(LM)*RFAC
                  END DO
               END DO
C
               IF ( IR.EQ.(NRMAX-NSR(NSHLR)) ) THEN
C     keep the value before the last shell to test convergence
                  DO L = 0,LMAX
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        S0(LM) = STESTNEW(LM)
                     END DO
                  END DO
               END IF
            END IF  ! r <> 0
C-----------------------------------------------------------------------
         END DO                 ! ir loop
C ======================================================================
C
C --> Check convergence
C
         S = 0D0
C ======================================================================
         DO L = 0,LMAX
            DO M = -L,L
               LM = L*(L+1) + M + 1
               STEST0 = ABS(S0(LM)-STESTNEW(LM))
               IF ( S.LT.STEST0 ) S = STEST0
            END DO
         END DO
C ======================================================================
         IF ( S.GT.BOUND ) WRITE (6,FMT=99001) ABS(S)
C
C --> Sum in reciprocal lattice
C
         CON = CONST_4PI/2D0/VOLUC_2D
C ======================================================================
         DO IM = 1,NGMAX
            G1 = GN2(1,IM)
            G2 = GN2(2,IM)
            G3 = 0D0
            GA = SQRT(G1*G1+G2*G2)
C ----------------------------------------------------------------------
            DOT1 = DQ1*G1 + DQ2*G2
            CALL FPLANEG(LAMBDA,GI,PREF0,LMAX_RGNTMAD,GA,VOLUC_2D)
            SIMAG = EXP(CI*DOT1)
            DO L = 0,4
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + GI(L)*SIMAG
            END DO
C ----------------------------------------------------------------------
            IF ( GA.GT.1D-6 ) THEN
               CALL CALC_RHPLM(G1/GA,G2/GA,G3/GA,YLM,LMAX_RGNTMAD,
     &                         NLMMAX_RGNTMAD)
               BETA = GA/LAMBDA
               EXPBSQ = DERFC(BETA/2D0)
C
               BFAC = CON*SIMAG*EXPBSQ
               STESTNEW(1) = STESTNEW(1) + BFAC/GA
C
               DO L = 0,LMAX
                  IF ( L.NE.0 ) THEN
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        STESTNEW(LM) = STESTNEW(LM) + YLM(LM)
     &                                 *BFAC*GA**(L-1)
                     END DO
                  END IF
                  BFAC = BFAC/CI/DBLE(2*L+1)
               END DO
            END IF
C ----------------------------------------------------------------------
            IF ( IM.EQ.(NGMAX-NSG(NSHLG)) ) THEN
C     keep the value before the last shell to test convergence
               DO LM = 1,LMMAX
                  S0(LM) = STESTNEW(LM)
               END DO
            END IF
         END DO
C ======================================================================
C
C --> Check convergence
C
         DO LM = 1,LMMAX
            STEST0 = ABS(S0(LM)-STESTNEW(LM))
            IF ( S.LT.STEST0 ) S = STEST0
         END DO
C
C --> Correction due to r=0 term only for DRn = 0
C
C ----------------------------------------------------------------------
         IF ( (DQ1*DQ1+DQ2*DQ2).LT.1D-6 ) THEN
            DO L = 2,4,2
               PREF0(L) = PREF0(L)*DFAC(L)/DFAC(L/2)/(L+1)
            END DO
C
            I = 1
            DO L = 2,4,2
               I = I + 1
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + (-1)**I*2D0/SQRT(PI)*PREF0(L)
     &                     *LAMBDA**(L+1)
            END DO
         END IF
C ----------------------------------------------------------------------
         DO L = 2,4,2
            LM = L*(L+1) + 1
            STESTNEW(LM) = STEST(LM)
         END DO
C --> end of correction
C ----------------------------------------------------------------------
C
         IF ( S.GT.BOUND ) WRITE (6,FMT=99001) ABS(S)
C
         DO LM = 1,LMMAX
            STEST(LM) = STESTNEW(LM)
         END DO
C **********************************************************************
      ELSE
C **********************************************************************
C ************************* OUT OF THE PLANE ***************************
C **********************************************************************
C
C --> Prepare arrays for speed up
C
         SIGNRZ = DQ3/DABS(DQ3)
         CON1 = CONST_2PI/VOLUC_2D
         DO L = 0,LMAX
            YLMPREF(L) = SQRT((2D0*L+1D0)/CONST_4PI)/DFAC(L)*CON1
            SIGNRZL(L) = (-SIGNRZ)**L
            CIM(L) = (-CI)**L
            DO M = 1,L
               YLMPREF1(L,M) = CON1*SQRT(DBLE(2*L+1)/2D0/CONST_4PI/DFAC(
     &                         L+M)/DFAC(L-M))
            END DO
         END DO
         YLMPREF(0) = 1D0*CON1
C
C --> Sum in reciprocal space
C
C ======================================================================
         DO I = 2,NGMAX
C
C   Exclude the origin all terms vanish for g=0
C   except the l = 0 components which are treated
C   sepparately. (look at the begining of the sub)
C
            G1 = GN2(1,I)
            G2 = GN2(2,I)
            GA = SQRT(G1*G1+G2*G2)
            DO L = 0,LMAX
               GAL(L) = GA**L
            END DO
C
            EXPBSQ = EXP(-GA*DABS(DQ3))
            DQDOTG = (DQ1*G1+DQ2*G2)
            FACTEXP = EXP(CI*DQDOTG)*EXPBSQ/GA
            EXPONL(0) = 1D0
            DO L = 1,LMAX
               EXPONL(L) = (G1+CI*G2)/GA*EXPONL(L-1)    ! exp(i*m*fi)
            END DO
C
C     In case rz < 0 then multiply by (-1)**(l-m)
C     (M. Weinert J. Math Phys. 22, 2439 (1981) formula 33
C      compare also formula A9)
C
            DO L = 0,LMAX
C     m = 0
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + YLMPREF(L)*GAL(L)*SIGNRZL(L)
     &                     *FACTEXP
C     m <> 0
               DO M = 1,L
                  PREF2(M) = CIM(M)*YLMPREF1(L,M)*GAL(L)
     &                       *SIGNRZ**M*SIGNRZL(L)
C
C Go from the <usual> Jackson Ylm to the ones we use
C
                  APREFPP = (1D0/EXPONL(M)+EXPONL(M))
                  APREFMM = (1D0/EXPONL(M)-EXPONL(M))*CI
C     m > 0
                  LM = L*(L+1) + M + 1
                  STEST(LM) = STEST(LM) + PREF2(M)*APREFPP*FACTEXP
C     m < 0
                  LM = L*(L+1) - M + 1
                  STEST(LM) = STEST(LM) + PREF2(M)*APREFMM*FACTEXP
               END DO
            END DO
C ----------------------------------------------------------------------
            IF ( I.EQ.(NGMAX-NSG(NSHLG)) ) THEN
C     keep the value before the last shell to test convergence
               DO LM = 1,LMMAX
                  S0(LM) = STEST(LM)
               END DO
            END IF
C ----------------------------------------------------------------------
         END DO
C ======================================================================
C
C --> Check convergence
C
         S = 0D0
         DO LM = 2,LMMAX
            STEST0 = ABS(S0(LM)-STEST(LM))
            IF ( S.LT.STEST0 ) S = STEST0
         END DO
         IF ( S.GT.BOUND ) WRITE (6,FMT=99001) ABS(S)
      END IF
C **********************************************************************
C
      DO LM = 1,LMMAX
         IF ( ABS(DIMAG(STEST(LM))).GT.BOUND ) THEN
            WRITE (6,*) ' ERROR: Imaginary contribution',
     &                  ' to REAL lattice sum'
            STOP
         END IF
         SMAT(LM) = DBLE(STEST(LM))
         STEST(LM) = 0.0D0
      END DO
C
      SMAT(1) = SMAT(1)/SQRT(CONST_4PI)
C
99001 FORMAT (5X,'WARNING : Convergence of 2D-sum is ',1P,D9.2,
     &        ' > 1D-8',/,15X,
     &        'You should use more lattice vectors (RMAX/GMAX)')
      END
C*==scfmad2d_coef.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE SCFMAD2D_COEF(LINTERFACE,NLMAD_LOC,A,B,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL LINTERFACE
      INTEGER NLMAD_LOC,NLMMAD_LOC,NLMMAX_RGNTMAD,NRGNTMAD,NRGNTMADMAX
      REAL*8 A(NLMMAD_LOC,NLMMAD_LOC),B(NLMMAD_LOC),RGNTMAD(NRGNTMADMAX)
     &       ,SMAT(NLMMAX_RGNTMAD)
      INTEGER IRGNTMAD(NRGNTMADMAX,3)
C
C Local variables
C
      REAL*8 DFAC(0:(NLMAD_LOC-1),0:(NLMAD_LOC-1))
      INTEGER I,L,L1,L2,LM1,LM2,LM3,LOFLM(NLMMAX_RGNTMAD),M
C
C*** End of declarations rewritten by SPAG
C
C --> determine the l-value for given lm
C
      I = 1
      DO L = 0,2*(NLMAD_LOC-1)
         DO M = -L,L
            LOFLM(I) = L
            I = I + 1
         END DO
      END DO
C
C --> calculate:                             (2*(l+l')-1)!!
C                 dfac(l,l') = 4pi**2 *  ----------------------
C                                        (2*l+1)!! * (2*l'+1)!!
C
      DFAC(0,0) = CONST_4PI*CONST_4PI
      DO L1 = 1,(NLMAD_LOC-1)
         DFAC(L1,0) = DFAC(L1-1,0)*DBLE(2*L1-1)/DBLE(2*L1+1)
         DFAC(0,L1) = DFAC(L1,0)
         DO L2 = 1,L1
            DFAC(L1,L2) = DFAC(L1,L2-1)*DBLE(2*(L1+L2)-1)/DBLE(2*L2+1)
            DFAC(L2,L1) = DFAC(L1,L2)
         END DO
      END DO
C
C --> initialize
C
      DO LM1 = 1,NLMMAD_LOC
         DO LM2 = 1,NLMMAD_LOC
            A(LM1,LM2) = 0.0D0
         END DO
      END DO
C
C --> calculate a(lm1,lm2)
C
      DO I = 1,NRGNTMAD
C
         LM1 = IRGNTMAD(I,1)
         LM2 = IRGNTMAD(I,2)
         LM3 = IRGNTMAD(I,3)
         L1 = LOFLM(LM1)
         L2 = LOFLM(LM2)
C
C --> this loop has to be calculated only for l1+l2=l3
C
         A(LM1,LM2) = A(LM1,LM2) + 2.0D0*DFAC(L1,L2)*SMAT(LM3)
     &                *RGNTMAD(I)
      END DO
C
      IF ( LINTERFACE ) RETURN
C
C --> initialize
C
      DO LM1 = 1,NLMMAD_LOC
         B(LM1) = 0.0D0
      END DO
C
C --> calculate b(lm1)
C
      DO LM1 = 1,NLMMAD_LOC
         L1 = LOFLM(LM1)
         B(LM1) = B(LM1) - 2.0D0*CONST_4PI/DBLE(2*L1+1)*SMAT(LM1)
      END DO
C
      END
C*==fplaner.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE FPLANER(ALPHA,G,R)
C ************************************************
C This sub calculates the derivatives of the real
C space contribution to the ewald sum .
C
C              l
C             d     erfc(lambda*sqrt(d*d+z*z))
C      lim    --   ------------------------
C      z->0     l        sqrt(d*d+z*z)
C             dz
C
C  Up to l = 4 (l=1,3,5,7 etc vanish)
C
C ************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI,SQRT_PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,R
      REAL*8 G(0:4)
C
C Local variables
C
      REAL*8 DERFC
      REAL*8 ER,EX,LAMBDA,PREF
      INTEGER L
C
C*** End of declarations rewritten by SPAG
C
      DO L = 0,4
         G(L) = 0.D0
      END DO
C
      ER = DERFC(ALPHA)
      EX = EXP(-ALPHA*ALPHA)
      LAMBDA = ALPHA/R
C
      G(0) = ER/R
C
      PREF = SQRT(5.D0/PI)/4.D0
      G(2) = -PREF*(ER/R/R/R+EX*2.D0*LAMBDA/R/R/SQRT_PI)
C
      PREF = 3.D0*SQRT(9.D0/PI)/16.D0/9.D0
      G(4) = PREF*(9.D0*ER+EX*(12.D0*ALPHA**3+18.D0*ALPHA)/SQRT_PI)
     &       /R/R/R/R/R
C
      END
C*==fplaneg.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE FPLANEG(LAMBDA,G,PREF,LMAX,GA,VOLUC_2D)
C **********************************************************************
C  This subroutine calculates the derivatives of the inverse
C  space contribution to the ewald sum
C
C
C     l
C    d   pi*( exp( gz)*erfc(g/lam/lam + 2*z)*lam/2
C           + exp(-gz)*erfc(g/lam/lam - 2*z)*lam/2 ) / (g Voluc_2d)
C    ---
C      l
C    dz
C
C   And the limit z -> 0 is taken (lam is the lambda parameter)
C
C *********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI,SQRT_PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 GA,LAMBDA,VOLUC_2D
      INTEGER LMAX
      REAL*8 G(0:4),PREF(0:LMAX)
C
C Local variables
C
      REAL*8 ALPHA,ER,EX
      REAL*8 DERFC
      INTEGER L
C
C*** End of declarations rewritten by SPAG
C
      DO L = 0,4
         G(L) = 0.D0
      END DO
      DO L = 0,LMAX
         PREF(L) = 0.D0
      END DO
C
      ALPHA = GA/2.D0/LAMBDA
      ER = DERFC(ALPHA)
      EX = EXP(-ALPHA*ALPHA)
C
      IF ( ABS(GA).GT.1.D-6 ) THEN
         G(0) = 2.D0*PI/VOLUC_2D*ER/GA
      ELSE
         G(0) = 0.D0
      END IF
C
      PREF(2) = SQRT(5.D0/PI)/2.D0/2.D0
      G(2) = PREF(2)/VOLUC_2D*(2.D0*PI*GA*ER-EX*4.D0*SQRT_PI*LAMBDA)
C
      PREF(4) = 3.D0*SQRT(9.D0/PI)/16.D0/9.D0
      G(4) = PREF(4)
     &       /VOLUC_2D*(2.D0*PI*GA*GA*GA*ER-EX*4.D0*SQRT_PI*(LAMBDA*GA*
     &       GA-2.D0*LAMBDA**3))
C
      END
C*==lattice2d.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE LATTICE2D(RMAX,GMAX,ALAT,ABAS_2D,BBAS_2D,NGMAX,NRMAX,
     &                     NSHLG,NSHLR,NSG,NSR,GN,RM,IPRINT,NMAXD,ISHLD)
C **********************************************************************
C *                                                                    *
C *  generate lattice vectors of direct and reciprocal space from      *
C *  basic translation vectors  ABAS_2D                                *
C *                                                                    *
C *  ALAT            : lattice constant                                *
C *  ABAS_2D(i,j)    : i=x,y,z j= 1,2,3 Bravais lattice basis vectors  *
C *                    *** in a.u. ****                                *
C *  rmax            : maximum radius in real space        (input)     *
C *  gmax            : maximum radius in reciprocal space  (input)     *
C *  ngmax           : Number of reciprocal lattice vectors            *
C *  gn(2,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
C *  nrmax           : Number of real lattice vectors                  *
C *  rm(2,nmaxd)     : x,y,z  of real space vectors                    *
C *  nshlg           : shells in reciprocal space                      *
C *  nshlr           : shells in real space                            *
C *  nsg,nsr         : integer arrays, number of atoms in each shell   *
C *                                                                    *
C *  The routine has been brought to a form which is very similar to   *
C *  LATTICE2D -- from which it has been originally derived            *
C *  Dimension of arrays GN,RM changed from (4,*) to (2,*), the 4th    *
C *  one it is used only locally (GNR/RMR) -- only GN/RM(2,*) are      *
C *  actually needed in SCFMAD2D_SMAT                                  *
C *                                                                    *
C **********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALAT,GMAX,RMAX
      INTEGER IPRINT,ISHLD,NGMAX,NMAXD,NRMAX,NSHLG,NSHLR
      REAL*8 ABAS_2D(3,3),BBAS_2D(3,3),GN(2,NMAXD),RM(2,NMAXD)
      INTEGER NSG(ISHLD),NSR(ISHLD)
C
C Local variables
C
      REAL*8 A,ABSG(3),ABSGM,ABSR(3),ABSRM,AG,AR,B,BG(3,3),BR(3,3),
     &       CJ(4,NMAXD),DA,DB,GMAX_AU,GNR(NMAXD),GX,GY,RMAX_AU,
     &       RMR(NMAXD),RX,RY,VMIN
      INTEGER I,K,L,M,N,N1,NG,NR,NSH,NSHL,NUMG,NUMGH,NUMR,NUMRH
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,'(5X,2A,/)') '<LATTICE2D> : ',
     &                    'generating direct/reciprocal lattice vectors'
C
      RMAX_AU = RMAX*ALAT
      GMAX_AU = GMAX/ALAT
C
      WRITE (6,FMT=99001) RMAX_AU,GMAX_AU
C
C ======================================================================
C
C --> basic trans. vectors and basis vectors
C
      DO I = 1,3
         BR(1,I) = ABAS_2D(1,I)*ALAT
         BR(2,I) = ABAS_2D(2,I)*ALAT
      END DO
C ======================================================================
C
C --> generate primitive vectors BG of reciprocal space
C
      DO I = 1,3
         BG(1,I) = BBAS_2D(1,I)*2D0*PI/ALAT
         BG(2,I) = BBAS_2D(2,I)*2D0*PI/ALAT
      END DO
C ======================================================================
C
C --> estimate no. of lattice vectors
C
      DO I = 1,3
         ABSR(I) = SQRT(BR(1,I)**2+BR(2,I)**2)
         ABSG(I) = SQRT(BG(1,I)**2+BG(2,I)**2)
      END DO
C
      ABSRM = MAX(ABSR(1),ABSR(2))
      ABSGM = MAX(ABSG(1),ABSG(2))
      ABSRM = 2.0D0*PI/ABSRM
      ABSGM = 2.0D0*PI/ABSGM
      NUMR = 2*(IDINT(RMAX_AU/ABSGM)+1) + 1
      NUMG = 2*(IDINT(GMAX_AU/ABSRM)+1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
C
C **********************************************************************
C                 generate lattice vectors of real space
C **********************************************************************
C
      NR = 0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO L = 1,NUMR
         A = DBLE(L-NUMRH)
         DO M = 1,NUMR
            B = DBLE(M-NUMRH)
C ----------------------------------------------------------------------
            RX = A*BR(1,1) + B*BR(1,2)
            RY = A*BR(2,1) + B*BR(2,2)
            AR = SQRT(RX*RX+RY*RY)
C ----------------------------------------------------------------------
            IF ( AR.LE.RMAX_AU ) THEN
               NR = NR + 1
               IF ( NR.GT.NMAXD ) THEN
                  WRITE (6,*) 
     &                      ' ERROR: Dimension NMAXD in inc.p too small'
     &                      ,NR,NMAXD
                  STOP
               END IF
               CJ(1,NR) = RX
               CJ(2,NR) = RY
               CJ(3,NR) = 0D0
               CJ(4,NR) = AR
            END IF
         END DO
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NRMAX = NR
C ======================================================================
C
C --> sort vectors in order of increasing absolute value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO K = 1,NR
         VMIN = RMAX_AU + 1.0D0
         DO N = 1,NR
            IF ( CJ(4,N)-VMIN.LT.0D0 ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         RM(1,K) = CJ(1,N1)
         RM(2,K) = CJ(2,N1)
         RMR(K) = CJ(4,N1)
         DB = VMIN
C ----------------------------------------------------------------------
         IF ( DB.GT.DA+1.D-06 ) THEN
            NSH = NSH + 1
            IF ( NSH.GT.ISHLD ) THEN
               WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',
     &                     NSH,ISHLD
               STOP
            END IF
C
            NSR(NSH) = NSHL
            NSHL = 0
            DA = DB
         END IF
C ----------------------------------------------------------------------
         CJ(4,N1) = RMAX_AU + 1.0D0
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.ISHLD ) THEN
         WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',NSH,
     &               ISHLD
         STOP
      END IF
C
      NSR(NSH) = NSHL
      NSHLR = NSH
      IF ( NSHLR.LE.1 ) STOP ' ERROR: cut-off radius RMAX_AU too small '
C
C **********************************************************************
C
C **********************************************************************
C                 generate lattice vectors of real space
C **********************************************************************
C
      NG = 0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO L = 1,NUMG
         A = DBLE(L-NUMGH)
         DO M = 1,NUMG
            B = DBLE(M-NUMGH)
C ----------------------------------------------------------------------
            GX = A*BG(1,1) + B*BG(1,2)
            GY = A*BG(2,1) + B*BG(2,2)
            AG = SQRT(GX*GX+GY*GY)
C ----------------------------------------------------------------------
            IF ( AG.LE.GMAX_AU ) THEN
               NG = NG + 1
               IF ( NG.GT.NMAXD ) THEN
                  WRITE (6,*) 
     &                      ' ERROR: Dimension NMAXD in inc.p too small'
     &                      ,NG,NMAXD
                  STOP
               END IF
               CJ(1,NG) = GX
               CJ(2,NG) = GY
               CJ(3,NG) = 0D0
               CJ(4,NG) = AG
            END IF
         END DO
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NGMAX = NG
C ======================================================================
C
C --> sort vectors in order of increasing abs. value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO K = 1,NG
         VMIN = GMAX_AU + 1.0D0
         DO N = 1,NG
            IF ( CJ(4,N).LT.VMIN ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         GN(1,K) = CJ(1,N1)
         GN(2,K) = CJ(2,N1)
         GNR(K) = CJ(4,N1)
         DB = VMIN
C ----------------------------------------------------------------------
         IF ( DB.GT.DA+1.D-07 ) THEN
            NSH = NSH + 1
            IF ( NSH.GT.ISHLD ) THEN
               WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',
     &                     NSH,ISHLD
               STOP
            END IF
C
            NSG(NSH) = NSHL
            NSHL = 0
            DA = DB
         END IF
C ----------------------------------------------------------------------
         CJ(4,N1) = GMAX_AU + 1.0D0
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.ISHLD ) THEN
         WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',NSH,
     &               ISHLD
         STOP
      END IF
C
      NSG(NSH) = NSHL
      NSHLG = NSH
      IF ( NSHLG.LE.1 ) STOP ' ERROR: cut-off radius GMAX_AU too small '
C
C **********************************************************************
C
C ----------------------------------------------------------------------
      WRITE (6,FMT=99002)
      WRITE (6,FMT=99003) 'Direct  lattice',NRMAX,NSHLR,RMR(NRMAX)
      WRITE (6,FMT=99003) 'Recipr. lattice',NGMAX,NSHLG,GNR(NGMAX)
      WRITE (6,FMT=99004)
C
      IF ( IPRINT.LT.3 ) RETURN
C ----------------------------------------------------------------------
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      K = 0
      WRITE (6,FMT=99005) 'real-space'
      DO L = 1,NSHLR
         WRITE (6,99006) L,NSR(L),RMR(K+1),(RM(M,K+1),M=1,2)
         DO N = 2,NSR(L)
            WRITE (6,FMT=99007) (RM(M,K+N),M=1,2)
         END DO
         IF ( L.NE.NSHLR ) WRITE (6,99008)
         K = K + NSR(L)
      END DO
      WRITE (6,99009)
      K = 0
      WRITE (6,FMT=99005) 'reciprocal'
      DO L = 1,NSHLG
         WRITE (6,99006) L,NSG(L),GNR(K+1),(GN(M,K+1),M=1,2)
         DO N = 2,NSG(L)
            WRITE (6,FMT=99007) (GN(M,K+N),M=1,2)
         END DO
         IF ( L.NE.NSHLG ) WRITE (6,99008)
         K = K + NSG(L)
      END DO
      WRITE (6,99009)
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (10X,'R max =',F10.5,' (a.u.)',/,10X,'G max =',F10.5,
     &        ' (1/a.u.)',/)
99002 FORMAT (10X,'               vectors  shells  max. R ',/,10X,
     &        '               ------------------------------')
99003 FORMAT (10X,A,I7,2X,I6,2X,F9.5)
99004 FORMAT (10X,'               ------------------------------',/)
99005 FORMAT (10X,45('+'),/,13X,'generated ',A,' lattice vectors',/,10X,
     &        45('+'),/,10X,'shell Nvec    radius          x         y',
     &        /,10X,45('-'))
99006 FORMAT (10X,I5,I5,F12.6,2X,2F10.5)
99007 FORMAT (34X,2F10.5)
99008 FORMAT (13X,42('-'))
99009 FORMAT (10X,45('+'),/)
      END
C*==madelgaunt.f    processed by SPAG 6.70Rc at 21:31 on 19 Dec 2016
      SUBROUTINE MADELGAUNT(NLMAD_LOC,RGNTMAD,IRGNTMAD,NRGNTMAD,
     &                      NRGNTMADMAX)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLMAD_LOC,NRGNTMAD,NRGNTMADMAX
      INTEGER IRGNTMAD(NRGNTMADMAX,3)
      REAL*8 RGNTMAD(NRGNTMADMAX)
C
C Local variables
C
      REAL*8 CLECG,FACTOR,S
      REAL*8 GAUNT_RYLM
      INTEGER I,L1,L2,L3,M1,M1A,M1S,M2,M2A,M2S,M3,M3A,M3S
C
C*** End of declarations rewritten by SPAG
C
C --> set up of the gaunt coefficients with an index field
C     recognize that they are needed here only for l3=l1+l2
C
      I = 1
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L1 = 0,(NLMAD_LOC-1)
         DO L2 = 0,(NLMAD_LOC-1)
            L3 = L1 + L2
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO M1 = -L1,L1
               DO M2 = -L2,L2
                  DO M3 = -L3,L3
                     M1S = SIGN(1,M1)
                     M2S = SIGN(1,M2)
                     M3S = SIGN(1,M3)
C **********************************************************************
                     IF ( M1S*M2S*M3S.GE.0 ) THEN
                        M1A = ABS(M1)
                        M2A = ABS(M2)
                        M3A = ABS(M3)
C
                        FACTOR = 0.0D0
                        IF ( M1A+M2A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(3*M3S+SIGN(1,-M3))/8.0D0
                        IF ( M1A-M2A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(M1S)/4.0D0
                        IF ( M2A-M1A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(M2S)/4.0D0
C ======================================================================
                        IF ( ABS(FACTOR).GT.1D-12 ) THEN
C
                           IF ( M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR. 
     &                          M1S*M3S.NE.1 ) FACTOR = -FACTOR
C
                           S = GAUNT_RYLM(L1,M1A,L2,M2A,L3,M3A)
C
                           CLECG = S*FACTOR
C ----------------------------------------------------------------------
                           IF ( ABS(CLECG).GT.1.D-10 ) THEN
                              IF ( I.GT.NRGNTMADMAX ) THEN
                                 WRITE (6,FMT='(2I10)') I,NRGNTMADMAX
                                 STOP ' Dim stop in MADELGAUNT '
                              END IF
                              RGNTMAD(I) = CLECG
                              IRGNTMAD(I,1) = L1*(L1+1) + M1 + 1
                              IRGNTMAD(I,2) = L2*(L2+1) + M2 + 1
                              IRGNTMAD(I,3) = L3*(L3+1) + M3 + 1
                              I = I + 1
                           END IF
C ----------------------------------------------------------------------
                        END IF
C ======================================================================
                     END IF
C **********************************************************************
                  END DO
               END DO
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      NRGNTMAD = I - 1
      END
