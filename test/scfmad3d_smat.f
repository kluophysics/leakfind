C*==scfmad3d_smat.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD3D_SMAT(ALAT,ABAS,QBAS,NQ_LOC,IPRINT,NLMAD,SMAT,
     &                         L3MAX,NLM3MAX)
C   ********************************************************************
C   *                                                                  *
C   * generate lattice vectors of direct and reciprocal space          *
C   * from basic translation vectors br                                *
C   *                                                                  *
C   * -----------------------------------------------------------------*
C   *                                                                  *
C   * calculation of lattice sums for l <= 2*lpot :                    *
C   *                                                                  *
C   *                  ylm( q(i) - q(j) - rm )                         *
C   *       sum     ----------------------------                       *
C   *                | q(i) - q(j) - rm |**(l+1)                       *
C   *                                                                  *
C   *        - summed over all lattice vectors rm  -                   *
C   *                                                                  *
C   * ylm       : real spherical harmic to given l,m                   *
C   * q(i),q(j) : basis vectors of the unit cell                       *
C   *                                                                  *
C   * in the case of i = j rm = 0 is omitted .                         *
C   *                                                                  *
C   * the ewald method is used to perform the lattice                  *
C   * summations the splitting parameter lamda is set                  *
C   * equal sqrt(pi)/alat (alat is the lattice constant) .             *
C   *                                                                  *
C   * if the contribution of the last shell of the direct              *
C   * and the reciprocal lattice is greater than 1.0e-8 a              *
C   * message is written                                               *
C   *                                                                  *
C   * based on   B. Drittler's  routines                               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CONST_4PI,PI,CI,SQRT_PI
      IMPLICIT NONE
C*--SCFMAD3D_SMAT38
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD3D_SMAT')
      REAL*8 BOUND,TOL,NTRY
      PARAMETER (BOUND=1.0D-10,TOL=1D-12,NTRY=150)
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER IPRINT,L3MAX,NLM3MAX,NLMAD,NQ_LOC
      REAL*8 ABAS(3,3),QBAS(3,NQ_LOC),SMAT(NQ_LOC,NQ_LOC,NLM3MAX)
C
C Local variables
C
      REAL*8 A,ABSG(3),ABSGM,ABSR(3),ABSRM,AG,ALPHA,AR,B,BETA,BG(3,3),
     &       BR(3,3),C,CJ(:,:),DA,DB,DQDOTG,DQIJ(3),EXPBSQ,FINC,G(:),GA,
     &       GDIR(3),GMAX,GMAX_RED,GN(:,:),GX,GY,GZ,LAMDA,QI(:,:),R,
     &       RDIR(3),RFAC,RM(:,:),RMAX,RMAX_RED,RX,RY,RZ,S,SMATPRV,VMIN,
     &       VOL,YLM(:)
      COMPLEX*16 BFAC,STEST(:)
      REAL*8 DDOT,DNRM2
      INTEGER I,I1,I2,IA_ERR,II,IQ,ITRY,J,JQ,K,L,LM,LMAX,LMMAX,M,N,N1,
     &        NG,NGE,NGS,NMAXD,NR,NRE,NRS,NSG(:),NSH,NSHL,NSHLD,NSHLG,
     &        NSHLR,NSR(:),NSTART,NUMG,NUMGH,NUMR,NUMRH
      EXTERNAL GAMFC
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE GN,RM,CJ,NSG,NSR,G,YLM,STEST,QI
C
      CALL TRACK_INFO(ROUTINE)
C
      SMAT(1:NQ_LOC,1:NQ_LOC,1:NLM3MAX) = 0D0
C
      WRITE (6,FMT=99004) TOL
C
C-----------------------------------------------------------------------
C               initialize real spherical harmonics
C
      ALLOCATE (YLM(NLM3MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD3D_SMAT->G'
C
C-----------------------------------------------------------------------
C
      ALLOCATE (G(0:L3MAX),STEST(NLM3MAX),QI(3,NQ_LOC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 .OR. NLM3MAX.LT.0 ) STOP 'alloc:SCFMAD3D_SMAT->G'
C
      NMAXD = 20000
      NSHLD = 500
C
      ITRY = 0
C     FINC = 1.05D0 original setting -- problems with Si - to be checked
      FINC = 2.05D0
      RMAX_RED = 4.0D0/FINC
      GMAX_RED = 28.5D0/FINC
C
      WRITE (6,FMT=99008)
      DO J = 1,3
         WRITE (6,FMT=99010) (ABAS(I,J),I=1,3)
      END DO
C
      BR(1:3,1:3) = ABAS(1:3,1:3)*ALAT
C
      WRITE (6,FMT=99009)
      DO J = 1,NQ_LOC
         WRITE (6,FMT=99010) (QBAS(I,J),I=1,3)
      END DO
C
C---> generate primitive vectors BG of reciprocal space
C
      DO I = 1,3
         I1 = 1 + MOD(I,3)
         I2 = 1 + MOD(I1,3)
C
         BG(1,I) = BR(2,I1)*BR(3,I2) - BR(2,I2)*BR(3,I1)
         BG(2,I) = BR(3,I1)*BR(1,I2) - BR(3,I2)*BR(1,I1)
         BG(3,I) = BR(1,I1)*BR(2,I2) - BR(1,I2)*BR(2,I1)
      END DO
C
      VOL = ABS(BR(1,1)*BG(1,1)+BR(2,1)*BG(2,1)+BR(3,1)*BG(3,1))
C
      WRITE (6,FMT=99011)
      DO I = 1,3
         BG(1:3,I) = BG(1:3,I)/VOL*2.0D0*PI
         WRITE (6,FMT=99010) BG(1,I),BG(2,I),BG(3,I)
      END DO
C
C-----------------------------------------------------------------------
 100  CONTINUE
      ITRY = ITRY + 1
      IF ( ITRY.GT.NTRY ) THEN
         WRITE (6,99001) NTRY
         STOP
      END IF
C
      RMAX_RED = RMAX_RED*FINC
      GMAX_RED = GMAX_RED*FINC
C
 200  CONTINUE
      IF ( ALLOCATED(GN) ) DEALLOCATE (NSG,NSR,GN,RM,CJ)
C
      ALLOCATE (NSG(NSHLD),NSR(NSHLD),STAT=IA_ERR)
      ALLOCATE (GN(4,NMAXD),RM(4,NMAXD),CJ(4,NMAXD),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD3D_SMAT->GN'
C
      RMAX = RMAX_RED*ALAT
      GMAX = GMAX_RED/ALAT
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99006) ALAT,RMAX,GMAX
C
      QI(1:3,1:NQ_LOC) = QBAS(1:3,1:NQ_LOC)*ALAT
C
C---> estimate no. of lattice vectors
C
      DO I = 1,3
         ABSR(I) = SQRT(BR(1,I)**2+BR(2,I)**2+BR(3,I)**2)
         ABSG(I) = SQRT(BG(1,I)**2+BG(2,I)**2+BG(3,I)**2)
      END DO
      ABSRM = MAX(ABSR(1),ABSR(2),ABSR(3))
      ABSGM = MAX(ABSG(1),ABSG(2),ABSG(3))
      ABSRM = 2.0D0*PI/ABSRM
      ABSGM = 2.0D0*PI/ABSGM
      NUMR = 2*(INT(RMAX/ABSGM)+1) + 1
      NUMG = 2*(INT(GMAX/ABSRM)+1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99017) NUMR,NUMG
C
C-----------------------------------------------------------------------
C              generate lattice vectors of real space
C-----------------------------------------------------------------------
C
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99012)
      NR = 0
      DO L = 1,NUMR
         A = DBLE(L-NUMRH)
         DO M = 1,NUMR
            B = DBLE(M-NUMRH)
            DO N = 1,NUMR
               C = DBLE(N-NUMRH)
               RX = A*BR(1,1) + B*BR(1,2) + C*BR(1,3)
               RY = A*BR(2,1) + B*BR(2,2) + C*BR(2,3)
               RZ = A*BR(3,1) + B*BR(3,2) + C*BR(3,3)
               AR = SQRT(RX*RX+RY*RY+RZ*RZ)
               IF ( AR.LE.RMAX ) THEN
                  NR = NR + 1
                  IF ( NR.GT.NMAXD ) THEN
                     NMAXD = INT(1.2D0*NMAXD)
                     WRITE (6,99018) 'NMAXD',NMAXD
                     GOTO 200
                  END IF
                  CJ(1,NR) = RX
                  CJ(2,NR) = RY
                  CJ(3,NR) = RZ
                  CJ(4,NR) = AR
               END IF
            END DO
         END DO
      END DO
C
      IF ( NR.GT.NMAXD ) STOP 'in <SCFMAD3D_SMAT>:  1 - latvec '
C
C---> sort vectors in order of increasing absolute value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
      DO K = 1,NR
         VMIN = RMAX + 1.0D0
         DO N = 1,NR
            IF ( CJ(4,N).LT.VMIN ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         RM(1:4,K) = CJ(1:4,N1)
         DB = VMIN
         IF ( DB.GT.DA+1.D-06 ) THEN
C
            NSH = NSH + 1
            IF ( NSH.GT.NSHLD ) THEN
               NSHLD = INT(1.2D0*NSHLD)
               WRITE (6,99018) 'NSHLD',NSHLD
               GOTO 200
            END IF
            IF ( IPRINT.GT.0 ) WRITE (6,FMT=99014) NSH,NSHL
            NSR(NSH) = NSHL
            IF ( IPRINT.GT.0 ) WRITE (6,FMT=99013) K,RM(1,K),RM(2,K),
     &                                RM(3,K),DB
            NSHL = 0
            DA = DB
C
         ELSE IF ( IPRINT.GT.0 ) THEN
            WRITE (6,FMT=99013) K,RM(1,K),RM(2,K),RM(3,K),DB
         END IF
C
         CJ(4,N1) = RMAX + 1.0D0
      END DO
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.NSHLD ) THEN
         NSHLD = INT(1.2D0*NSHLD)
         WRITE (6,99018) 'NSHLD',NSHLD
         GOTO 200
      END IF
      NSR(NSH) = NSHL
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99014) NSH,NSHL
C
      IF ( NSH.GT.NSHLD ) STOP 'in <SCFMAD3D_SMAT>:  2 - latvec '
C
      NSHLR = NSH
C
C-----------------------------------------------------------------------
C              generate lattice vectors of reciprocal space
C-----------------------------------------------------------------------
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99015)
      NG = 0
      DO L = 1,NUMG
C
         A = DBLE(L-NUMGH)
         DO M = 1,NUMG
            B = DBLE(M-NUMGH)
C
            DO N = 1,NUMG
               C = DBLE(N-NUMGH)
               GX = A*BG(1,1) + B*BG(1,2) + C*BG(1,3)
               GY = A*BG(2,1) + B*BG(2,2) + C*BG(2,3)
               GZ = A*BG(3,1) + B*BG(3,2) + C*BG(3,3)
               AG = SQRT(GX*GX+GY*GY+GZ*GZ)
               IF ( AG.LE.GMAX ) THEN
                  NG = NG + 1
                  IF ( NG.GT.NMAXD ) THEN
                     NMAXD = INT(1.2D0*NMAXD)
                     WRITE (6,99018) 'NMAXD',NMAXD
                     GOTO 200
                  END IF
                  CJ(1,NG) = GX
                  CJ(2,NG) = GY
                  CJ(3,NG) = GZ
                  CJ(4,NG) = AG
               END IF
            END DO
         END DO
      END DO
C
      IF ( NG.GT.NMAXD ) STOP 'in <SCFMAD3D_SMAT>:  3 - latvec '
C
C---> sort vectors in order of increasing abs. value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
      DO K = 1,NG
         VMIN = GMAX + 1.0D0
         DO N = 1,NG
            IF ( CJ(4,N).LT.VMIN ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         GN(1:4,K) = CJ(1:4,N1)
         DB = VMIN
         IF ( DB.GT.DA+1.D-07 ) THEN
C
            NSH = NSH + 1
            IF ( NSH.GT.NSHLD ) THEN
               NSHLD = INT(1.2D0*NSHLD)
               WRITE (6,99018) 'NSHLD',NSHLD
               GOTO 200
            END IF
            IF ( IPRINT.GT.0 ) WRITE (6,FMT=99014) NSH,NSHL
            IF ( IPRINT.GT.0 ) WRITE (6,FMT=99013) K,GN(1,K),GN(2,K),
     &                                GN(3,K),DB
            NSG(NSH) = NSHL
            NSHL = 0
            DA = DB
C
         ELSE IF ( IPRINT.GT.0 ) THEN
            WRITE (6,FMT=99013) K,GN(1,K),GN(2,K),GN(3,K),DB
         END IF
C
         CJ(4,N1) = GMAX + 1.0D0
      END DO
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.NSHLD ) THEN
         NSHLD = INT(1.2D0*NSHLD)
         WRITE (6,99018) 'NSHLD',NSHLD
         GOTO 200
      END IF
      NSG(NSH) = NSHL
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99014) NSH,NSHL
      IF ( NSH.GT.NSHLD ) STOP 'in <SCFMAD3D_SMAT>:  4 - latvec '
      NSHLG = NSH
      IF ( IPRINT.GT.0 ) WRITE (6,FMT=99016) NR,NG,VOL
C
C=======================================================================
C                      perform summations
C=======================================================================
C
      LMAX = 2*(NLMAD-1)
      LMMAX = (LMAX+1)*(LMAX+1)
C
C---> choose proper splitting parameter
C
      LAMDA = SQRT_PI/ALAT
C
C---> loop over atoms per unit cell
C
      DO IQ = 1,NQ_LOC
         DO JQ = 1,NQ_LOC
C
            DQIJ(1:3) = QI(1:3,IQ) - QI(1:3,JQ)
C
            STEST(1) = -SQRT(CONST_4PI)/VOL/(4.0D0*LAMDA*LAMDA)
            DO LM = 2,LMMAX
               STEST(LM) = 0.0D0
            END DO
C
C---> exclude the origin and add correction if IQ.eq.JQ
C
            IF ( IQ.EQ.JQ ) THEN
               STEST(1) = STEST(1) - LAMDA/PI
               NSTART = 2
            ELSE
               NSTART = 1
            END IF
C
C---> loop first over n-1 shells of real and recipro. lattice - then
C      add the contribution of the last shells to see convergence
C
            DO II = 1,2
               IF ( II.EQ.1 ) THEN
                  NRS = NSTART
                  NGS = 2
                  NRE = NR - NSR(NSHLR)
                  NGE = NG - NSG(NSHLG)
C
               ELSE
                  NRS = NRE + 1
                  NGS = NGE + 1
                  NRE = NR
                  NGE = NG
               END IF
C
C---> sum over real lattice
C
               DO I = NRS,NRE
                  RDIR(1:3) = DQIJ(1:3) - RM(1:3,I)
                  R = DNRM2(3,RDIR,1)
                  RDIR(1:3) = RDIR(1:3)/R
C
                  CALL CALC_RHPLM(RDIR(1),RDIR(2),RDIR(3),YLM,L3MAX,
     &                            NLM3MAX)
C
                  ALPHA = LAMDA*R
C
                  CALL GAMFC(ALPHA,G,LMAX,R)
C
                  DO L = 0,LMAX
C
                     RFAC = G(L)/SQRT_PI
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        STEST(LM) = STEST(LM) + YLM(LM)*RFAC
                     END DO
                  END DO
C
               END DO
C
C---> sum over reciprocal lattice
C
               DO I = NGS,NGE
                  GDIR(1:3) = GN(1:3,I)
                  GA = DNRM2(3,GDIR,1)
                  GDIR(1:3) = GDIR(1:3)/GA
C
                  CALL CALC_RHPLM(GDIR(1),GDIR(2),GDIR(3),YLM,L3MAX,
     &                            NLM3MAX)
C
                  BETA = GA/LAMDA
                  EXPBSQ = EXP(BETA*BETA/4.0D0)
                  DQDOTG = DDOT(3,DQIJ,1,GN(1,I),1)
                  BFAC = CONST_4PI*EXP(CI*DQDOTG)/(GA*GA*EXPBSQ*VOL)
C
                  DO L = 0,LMAX
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        STEST(LM) = STEST(LM) + YLM(LM)*BFAC
                     END DO
                     BFAC = BFAC*GA/CI/DBLE(2*L+1)
                  END DO
C
               END DO
C
               IF ( II.EQ.1 ) THEN
C
                  DO LM = 1,LMMAX
                     IF ( ABS(DIMAG(STEST(LM))).GT.BOUND ) GOTO 300
                     SMAT(IQ,JQ,LM) = DBLE(STEST(LM))
                     STEST(LM) = 0.0D0
                  END DO
C
               ELSE
C
C---> test convergence
C if the relative change due to last contribution is larger than TOL
C the radii RMAX and GMAX are BOTH increased and the summation is redone
C
                  DO LM = 1,LMMAX
                     SMATPRV = SMAT(IQ,JQ,LM)
                     S = DBLE(STEST(LM))
                     SMAT(IQ,JQ,LM) = SMAT(IQ,JQ,LM) + S
                     IF ( ABS(S).GT.BOUND ) WRITE (6,FMT=99005) IQ,JQ,
     &                    LM,ABS(S)
                     IF ( ABS(SMATPRV).GT.1D-12 ) THEN
                        IF ( ABS(S/SMATPRV).GT.TOL ) GOTO 100
                     END IF
                  END DO
               END IF
C
            END DO
C
            IF ( IPRINT.GT.0 ) THEN
               S = 1D-8
               WRITE (6,99002) S
               DO LM = 1,LMMAX
                  IF ( ABS(SMAT(IQ,JQ,LM)).GT.S ) WRITE (6,99003) IQ,JQ,
     &                 LM,SMAT(IQ,JQ,LM)
               END DO
            END IF
C
         END DO
      END DO
C
      WRITE (6,99007) ITRY,ALAT,RMAX,GMAX,NUMR,NUMG
C
      DEALLOCATE (G,YLM,STEST,GN,RM,CJ,NSG,NSR,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD3D_SMAT->NSR'
C
      RETURN
C
C-----------------------------------------------------------------------
 300  CONTINUE
      STOP 
     &  'in <SCFMAD3D_SMAT>: imaginary contribution to real lattice sum'
C
99001 FORMAT (/,10X,'<SCFMAD3D_SMAT>: number of attempts ',
     &        'exceeded maximum',I3,/)
99002 FORMAT (/,10X,'elements of SMAT >',1P,E9.1,//,10X,
     &        '  IQ  JQ  LM          SMAT ')
99003 FORMAT (10X,3I4,F25.15)
99004 FORMAT (//,1X,79('*'),/,32X,'<SCFMAD3D_SMAT>',/,1X,79('*'),//,10X,
     &        'calculate Madelung potential matrices ',//,10X,
     &        'relative tolerance imposed ',1P,1E9.1,/)
99005 FORMAT (1x,' convergence of SMAT(',i2,i2,i4,') : ',D12.5,
     &        ' is less than 1.0e-8 - use more lattice vectors ')
99006 FORMAT (/,10X,'lattice constant    ALAT:',F20.15,/,10X,'RMAX:',
     &        F10.5,5X,'GMAX:',F10.5)
99007 FORMAT (/,10X,'summations converged after',I4,' updates for',//,
     &        10X,'lattice constant    ALAT:',F20.15,/,10X,'RMAX:',
     &        F10.5,5X,'GMAX:',F10.5,/,10X,'NUMR:',I10,5X,'NUMG:',I10,/)
99008 FORMAT (/,10X,'primitive vectors for Bravais lattice',/)
99009 FORMAT (/,10X,'basis vectors in units of a'/)
99010 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
99011 FORMAT (/,10X,'primitive vectors of reciprocal space',/,10X,
     &        '       in units of 2*pi/a',/)
99012 FORMAT (/,10X,' real space lattice vectors ',/)
99013 FORMAT (10X,I5,4F10.5)
99014 FORMAT (10X,' shell no.',I5,' contains',I5,' points ',/)
99015 FORMAT (/,10X,'reciprocal space lattice vectors ',/)
99016 FORMAT (/,10X,'no. of lattice vectors   : ',I5,/,10X,
     &        'no. of rec. lat. vectors : ',I5,/,10X,
     &        'volume of the unit cell  : ',F10.5,/)
99017 FORMAT (/,10X,'NUMR: ',I3,5X,'NUMG: ',I3,/)
99018 FORMAT (/,10X,'INFO from <SCFMAD3D_SMAT>:  ',A,' set to ',I8)
      END
C*==gamfc.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE GAMFC(ALPHA,GLH,LMAX,R)
C   ********************************************************************
C   *                                                                  *
C   *  calculation of convergence function                             *
C   *                                                                  *
C   *   glh = i(alpha,l)/r**(l+1)*sqrt(pi)                             *
C   *                                                                  *
C   *  with                                                            *
C   *        alpha = r times the splitting paramter lamda              *
C   *  and                                                             *
C   *        i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *                   *
C   *                                                                  *
C   *                     sum ( 2**i * x**(2i-1) / (2i-1)!! )          *
C   *                   1..i..l                                        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--GAMFC556
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,R
      INTEGER LMAX
      REAL*8 GLH(0:LMAX)
C
C Local variables
C
      REAL*8 ARG,ERFCEX,EX,FACL,FEX
      INTEGER L
C
C*** End of declarations rewritten by SPAG
C
      ARG = ALPHA*ALPHA
C
      CALL COMERRFUN(ALPHA,EX,ERFCEX)
C
      GLH(0) = ERFCEX
C
      FACL = 2.0D0*ALPHA
C
C---> recursion
C
      DO L = 1,LMAX
         GLH(L) = GLH(L-1) + FACL
         FACL = FACL*ARG/(DBLE(L)+0.5D0)
      END DO
      FEX = 1.0D0/EXP(ARG)
      DO L = 0,LMAX
         FEX = FEX/R
         GLH(L) = GLH(L)*FEX
      END DO
C
      END
C*==comerrfun.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE COMERRFUN(Z,EXMZZ,ERFCEX)
C **********************************************************************
C *                                                                    *
C *      calculates complementary error function                       *
C *                 multiplied by sqrt(pi) * exp(z*z)                  *
C *          z .......: argument of error function complement          *
C *          exmzz ...: exp(-z*z)                                      *
C *                                                                    *
C * In order to avoid numerical problems with big z exmzz=exp(-z*z)    *
C * is returned                                                        *
C *                                                                    *
C **********************************************************************
      USE MOD_CONSTANTS,ONLY:SQRT_PI
      IMPLICIT NONE
C*--COMERRFUN620
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 BOUND
      PARAMETER (BOUND=3.D-11)
C
C Dummy arguments
C
      REAL*8 ERFCEX,EXMZZ,Z
C
C Local variables
C
      REAL*8 ERF1,EXZZ,F,FA,Q,RATIO,TERM,U,UA,V,X,XA,Y,Z2,ZZ
C
C*** End of declarations rewritten by SPAG
C
      ZZ = Z*Z
      IF ( ZZ.GT.150D0 ) THEN
         EXMZZ = 0.0D0
      ELSE
         EXMZZ = EXP(-ZZ)
      END IF
C
C  chose algorithm
      IF ( Z.LT.1.5D0 ) THEN
C
         EXZZ = EXP(ZZ)
         Z2 = 2.D0*ZZ
C  series expansion of error function: abramowitz p.297 eq. (7.1.6)
         ERF1 = Z
         RATIO = 1.D0
         TERM = Z
 50      CONTINUE
         RATIO = RATIO + 2.D0
         TERM = TERM*Z2/RATIO
         ERF1 = ERF1 + TERM
         IF ( TERM.GT.BOUND ) GOTO 50
         ERFCEX = SQRT_PI*EXZZ - 2.D0*ERF1
C
      ELSE
C  continued fraction expansion : abramowitz p.298, eqn. (7.1.14)
         U = 1.D0
         V = 0.D0
         X = Z
         Y = 1.D0
         Q = .5D0
         F = (U+V*Q)/(X+Y*Q)
 100     CONTINUE
         UA = U
         U = U*Z + V*Q
         V = UA
         XA = X
         X = X*Z + Y*Q
         Y = XA
         Q = Q + 0.5D0
         FA = F
         F = (U+V*Q)/(X+Y*Q)
         IF ( ABS(FA-F).GT.BOUND*F ) GOTO 100
         ERFCEX = F
      END IF
C
      END
