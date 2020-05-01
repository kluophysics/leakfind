C*==xmatr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XMATR(POS,GEVEN,GODD,E,ELM,VPI,EMACH,AK1,AK2,APQ,BPQ,
     &                 NATIN,NALM,NREL,CLIGHT,ITT,GODEV,GEVOD,NEV,NOD)
C**************************************************************
C       xmatkc calculates the matrix describing multiple
C       scattering within a layer, returning it as: xalm,
C       corresponding to odd l+m, with lm=(10),(2-1),(21),
C       .... followed by even l+m, with
C       lm=(00),(1-1),(11),(2-2),....
C       the program assumes that the layer is a bravais lattice.
C       the diagonal blocks in atomic subscript are calculated
C       via the original xmatk from peover ,xmatkc provides off-diagonal
C       elements and (rumpled layer)'s contribution
C       the summation over the lattice follows the ewald split
C       method suggested by kambe. the code is based on xmatkc
C       from f.maca & m.scheffler cpc 1989. the dimension are fixed
C       for  lmax=4. emach is the machine accuracy: 1.0d-6 .
C
C                                                             *
C subroutines and functions called from this routine          *
C                                                             *
C zcerf      sphrkc    trans                                  *
C                                                             *
C**************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQ,MLQNAT,MLZP,PI,CZERO,CONE,
     &    CIMAG
      USE MOD_SPEC_GEOM,ONLY:AR1,AR2,TV
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLR,MLT,MLA,MLC,MX2,LMAXM,L2MAXM,NDLMM
      PARAMETER (MLR=(ML+1)*ML/2,MLT=ML*(ML-1)/2,MLA=2*ML-1,MLC=4*ML-3,
     &           MX2=NATLM*ML,LMAXM=ML-1,L2MAXM=2*LMAXM,NDLMM=(L2MAXM+1)
     &           *(L2MAXM+1))
C
C Dummy arguments
C
      REAL*8 AK1,AK2,APQ,BPQ,CLIGHT,E,EMACH,VPI
      INTEGER ITT,NALM,NATIN,NEV,NOD,NREL
      REAL*8 ELM(MLZP),POS(3,NATLM,LAYSM)
      COMPLEX*16 GEVEN(NEV,NEV),GEVOD(NEV,NOD),GODD(NOD,NOD),
     &           GODEV(NOD,NEV)
C
C Local variables
C
      COMPLEX*16 A,AA,AB,ACC,AF(:),AGK(:),ALM,ALPHA,CEN,CF,CFF,CI,CP,
     &           CTAL,CTT,CUNIT,CX,CZ,CZ2,D(:),DLM(:),FAC(:),GK,GKCIJ(:)
     &           ,GKK,GKN(:),GP,GPSQ,KANT,KAPPA,KAPSQ,KNSQ,PREF(:),RTA,
     &           SD,STT,U,U1,U2,UODD1,UODD2,W,W1,W2,WW,XALM(:,:),XPA,
     &           XPK,XPM(:),YLMA(:),Z,Z2,ZZ
      REAL*8 AB1,AB2,AC,ACSQ,AD,AK(2),AKPT(2),AL,AN,AN1,AN2,AP,AP1,AP2,
     &       AR,B,B1(2),B2(2),C,DENOM(200),DMEL,DNORM,DVEL,FOUR,PI4,
     &       R(2),RTPI,RTV,S(3),TAL,TEST,TEST1,TEST2
      INTEGER I,I1,I2,I3,I4,IA,ICGL,ICOUNT,IE,IETWO,IF2,II,II1,III,IL,
     &        IL2,IL3,IOTWO,ISS,J,J1,J2,JCOUNT,JE,JETWO,JJ,JOTWO,JS4,
     &        JS5,K,KK,L,L1,L11,L2,L2MAX,L3,L30,LA1,LA11,LAF1,LB1,LB11,
     &        LL1,LL2,LMAX,LOOP,LOTWO,LST,M,M1,M2,M3,MM,N,N1,NA,NCOUNT,
     &        NDIM,NM,NN,NN1,NNDLM,NSUDE
      COMPLEX*16 ZCERF
      EXTERNAL ZCERF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE D,AF,FAC,AGK,DLM,GKN,XPM,PREF,XALM,YLMA,GKCIJ
      ALLOCATE (D(MLA),AF(MX2),FAC(MLC),AGK(MLA),DLM(NDLMM))
      ALLOCATE (GKN(MLA),XPM(MLA),PREF(NDLMM),XALM(MLQNAT,MLQNAT))
      ALLOCATE (YLMA(NDLMM),GKCIJ(MLA))
C
      CUNIT = CONE
      CI = CIMAG
      RTPI = SQRT(PI)
      FOUR = 4.D0
      PI4 = FOUR*PI
      DO I = 1,MX2
         AF(I) = DCMPLX(0.D0,-1.D0)
      END DO
      CEN = DCMPLX(E,VPI+0.5D0*EMACH)
      KAPSQ = 2.D0*CEN
      IF ( NREL.EQ.1 ) KAPSQ = CEN*(2.D0+CEN/(CLIGHT*CLIGHT))
      KAPPA = CDSQRT(KAPSQ)
      AK(1) = AK1 + APQ
      AK(2) = AK2 + BPQ
      LMAX = ML - 1
      LAF1 = ML
      L2MAX = 2*ML - 2
      LL2 = L2MAX + 1
      NDIM = MLQ
      FAC(1) = 1.D0
      II = L2MAX + L2MAX
      DO I = 1,II
         FAC(I+1) = DBLE(I)*FAC(I)
      END DO
      NNDLM = LL2*LL2
      ALPHA = TV*KAPSQ/PI4
      AL = CDABS(ALPHA)
      IF ( EXP(AL)*EMACH-5.D-5.GT.0. ) THEN
         AL = DLOG(5.D-5/EMACH)
         ALPHA = DCMPLX(SIGN(1.D0,E)*AL,0.D0)
      END IF
      RTA = CDSQRT(ALPHA)
      DO I = 1,NALM
         DO J = 1,NALM
            XALM(I,J) = CZERO
         END DO
      END DO
      DO I = 1,MLR
         DO J = 1,MLR
            XALM(I,J) = GEVEN(I,J)
         END DO
      END DO
      DO I = 1,MLT
         DO J = 1,MLT
            IE = I + MLR
            JE = J + MLR
            XALM(IE,JE) = GODD(I,J)
         END DO
      END DO
      DO IA = 2,NATIN
         DO I = 1,MLQ
            DO J = 1,MLQ
               JE = J + (IA-1)*MLQ
               IE = I + JE - J
               XALM(IE,JE) = XALM(I,J)
            END DO
         END DO
      END DO
      IF ( NATIN.NE.1 ) THEN
C
C     offdiagonal blocs of xalm(atomi,atomj) are further calculated.
C
         LOOP = NATIN*NATIN
         DO NCOUNT = 1,LOOP
            JCOUNT = 1 + INT(1.D0*(NCOUNT-1)/NATIN)
            ICOUNT = NCOUNT - (JCOUNT-1)*NATIN
            IF ( ICOUNT.EQ.JCOUNT ) CYCLE
            S(1) = POS(1,ICOUNT,ITT) - POS(1,JCOUNT,ITT)
            S(2) = POS(2,ICOUNT,ITT) - POS(2,JCOUNT,ITT)
            S(3) = POS(3,ICOUNT,ITT) - POS(3,JCOUNT,ITT)
            DO I = 1,NNDLM
               DLM(I) = CZERO
            END DO
C
C         dlm1, the sum over reciprocal lattice vectors, is
C         calculated first. the prefactor 'p1' is tabulated
C         for even values of l+|m|, thus lm=(00),(11),(20),
C         (22),....(l2max,l2max) followed by odd values of l+m
C         (10),(21),(30),(32),.... the factorial factor 'f1' is simulta-
C         neously tabulated in denom for all values of n=0,...,l-|m|
C         and s=n,...,min(2*n,l-|m|)
C
            II = 0
            K = 1
            KK = 1
 20         CONTINUE
            AP1 = -2.D0/TV
            AP2 = -1.D0
            CF = CI/KAPPA
            DO L = 1,LL2
               AP1 = AP1/2.D0
               AP2 = AP2 + 2.D0
               CP = CF
               MM = 1 + II
               IF ( MOD(L,2).EQ.0 ) MM = 2 - II
               IF ( MM.EQ.2 ) CP = CI*CP
               NN = L - MM + 3
               DO M = MM,L,2
                  J1 = L + M - 1
                  J2 = L - M + 1
C                  AP = AP1*CDSQRT(AP2*FAC(J1)*FAC(J2))
                  AP = AP1*SQRT(AP2*REAL(FAC(J1))*REAL(FAC(J2)))
                  PREF(KK) = AP*CP
                  CP = -CP
                  KK = KK + 1
                  NN = NN - 2
                  DO I = 1,NN
                     JS5 = NN
                     IF ( NN.GT.I+I-1 ) JS5 = I + I - 1 - II
                     IF ( JS5.NE.0 ) THEN
                        JS4 = JS5
 22                     CONTINUE
                        IF ( JS4-2.LT.I ) THEN
                           DO ISS = JS4,JS5,2
                              I1 = I + I - ISS
                              I2 = ISS - I + 1
                              I3 = (NN-ISS)/2 + 1
                              I4 = (NN+M+M-ISS)/2
                              DENOM(K)
     &                           = 1.D0/(REAL(FAC(I1))*REAL(FAC(I2))
     &                           *REAL(FAC(I3))*REAL(FAC(I4)))
                              K = K + 1
                           END DO
                        ELSE
                           JS4 = JS4 - 2
                           GOTO 22
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            IF ( II.EQ.1 ) THEN
C
C         the reciprocal lattice is defined by b1,b2. the summation
C         begins with the origin point of the lattice, and continues
C         in steps of 8*n1 points, each step involving the perimeter
C         of a parallelogram of lattice points about the origin, of
C         side 2*n1+1. each step begins at label 9.
C         akpt = the current lattice vector in the sum
C
               RTV = 2.D0*PI/TV
               B1(1) = -AR1(2)*RTV
               B1(2) = AR1(1)*RTV
               B2(1) = -AR2(2)*RTV
               B2(2) = AR2(1)*RTV
               TEST1 = 1.D06
               III = 1
               N1 = -1
 30            CONTINUE
               N1 = N1 + 1
               NA = N1 + N1 + III
               AN1 = DBLE(N1)
               AN2 = (-AN1) - 1.D0
               DO I1 = 1,NA
                  AN2 = AN2 + 1.D0
                  DO I2 = 1,4
                     AN = AN1
                     AN1 = -AN2
                     AN2 = AN
                     AB1 = AN1*B1(1) + AN2*B2(1)
                     AB2 = AN1*B1(2) + AN2*B2(2)
                     AKPT(1) = AK(1) + AB1
                     AKPT(2) = AK(2) + AB2
C
C          gkcij(j) contains values of (kappa*cijz)**(2n-s)
C          gkn(n) contains values of (gp/kappa)**(2*n-1)*del<n,p>
C          where l=0,l2max; m=-l,l; n=0,l-|m|; i=l-s; s=n,min(2*n,l-|m|)
C          del are the integrals evaluated by recurrence from the values
C          for n=0 and n=1 which can be expressed in terms of
C          complex*16 error function cerf
C
                     ACSQ = AKPT(1)*AKPT(1) + AKPT(2)*AKPT(2)
                     GPSQ = KAPSQ - ACSQ
                     AC = SQRT(ACSQ)
                     GP = SQRT(GPSQ)
                     XPK = CZERO
                     GK = CZERO
                     GKK = CUNIT
                     IF ( AC-EMACH.GT.0.D0 ) THEN
                        XPK = DCMPLX(AKPT(1)/AC,AKPT(2)/AC)
                        GK = AC/KAPPA
                        GKK = GPSQ/KAPSQ
                     END IF
                     XPM(1) = CUNIT
                     AGK(1) = CUNIT
                     GKCIJ(1) = CUNIT
                     DO I = 2,LL2
                        XPM(I) = XPM(I-1)*XPK
                        GKCIJ(I) = GKCIJ(I-1)*KAPPA*S(1)
                        AGK(I) = AGK(I-1)*GK
                     END DO
                     CF = KAPPA/GP
                     ZZ = -ALPHA*GKK
                     CZ = DCMPLX(0.D0,-1.D0)*SQRT((-ZZ))
                     Z2 = GPSQ*S(1)*S(1)
                     CZ2 = SQRT(Z2)
                     CX = EXP((-ZZ)+DCMPLX(0.25D0,0.D0)*Z2/ZZ)
                     U1 = DCMPLX(-0.5D0,0.D0)*CZ2/CZ + CI*CZ
                     U2 = DCMPLX(0.5D0,0.D0)*CZ2/CZ + CI*CZ
                     IF ( DBLE(CZ).LT.0.D0 .AND. DBLE(ZZ).GT.150.D0 )
     &                    THEN
                        D(1) = DCMPLX(2.D0*RTPI,0.D0)
                        IF ( DIMAG(ZZ).GT.0.D0 ) D(1) = D(1)
     &                       *EXP(DCMPLX(0.D0,(-2.D0*DIMAG(ZZ))))
                        CX = CZERO
                     ELSE
                        W1 = ZCERF(U1,EMACH)
                        W2 = ZCERF(U2,EMACH)
                        D(1) = DCMPLX(0.5D0,0.D0)*RTPI*CX*(W1+W2)
C                        IF ( S(1).NE.0.D0 ) THEN
                        IF ( ABS(S(1)).GT.1.D-16 ) THEN
                           D(2) = CI*RTPI/CZ2*CX*(W1-W2)
                           DO I = 3,LL2
                              CX = CX/ZZ
                              B = 2.5D0 - DBLE(I)
                              D(I) = DCMPLX(B,0.D0)*D(I-1) - D(I-2)
     &                               + CX*CZ
                              D(I) = 4.D0*D(I)/Z2
                           END DO
                           GOTO 32
                        END IF
                     END IF
                     DO I = 2,LL2
                        CX = CX/ZZ
                        B = 1.D0/(1.5D0-DBLE(I))
                        D(I) = B*(D(I-1)-CX*CZ)
                     END DO
 32                  CONTINUE
                     GKN(1) = CF*D(1)
                     DO I = 2,LL2
                        CF = CF*GKK
                        GKN(I) = CF*D(I)
                     END DO
C
C                 the contribution to the sum dlm1 for a particular
C                 reciprocal lattice vector is now accumulated into
C                 the elements of dlm. note special action if ac=0
C
                     II = 0
                     K = 1
                     KK = 1
                     NSUDE = LL2*(LL2+1)/2
 34                  CONTINUE
                     LST = 1 + II
                     DO L = LST,LL2
                        MM = 1 + II
                        IF ( MOD(L,2).EQ.0 ) MM = 2 - II
                        N = ((L-II)*(L-II)+MM)/2 + NSUDE*II
                        NN = L - MM + 3
                        DO M = MM,L,2
                           ACC = CZERO
                           NN = NN - 2
                           IF ( AC-EMACH.LE.0.D0 ) THEN
                              IF ( M.EQ.1 ) THEN
                                 JS4 = L/2 + 1
                                 IF ( MOD(L,2).EQ.1 ) JS4 = (L+1)/2
                                 DO I = JS4,L
                                    IF2 = I + I - L
                                    I3 = L - I + 1
                                    ACC = ACC + GKN(I)*GKCIJ(IF2)
     &                                 /FAC(I3)/FAC(IF2)
                                 END DO
                              END IF
                              GOTO 38
                           END IF
                           DO I = 1,NN
                              JS5 = NN
                              IF ( NN.GT.I+I-1 ) JS5 = I + I - 1 - II
                              IF ( JS5.NE.0 ) THEN
                                 JS4 = JS5
 36                              CONTINUE
                                 IF ( JS4-2.LT.I ) THEN
                                    DO ISS = JS4,JS5,2
                                       JJ = I + I - ISS
                                       IL = L - ISS + 1
                                       ACC = ACC + DENOM(K)*AGK(IL)
     &                                    *GKN(I)*GKCIJ(JJ)
                                       K = K + 1
                                    END DO
                                 ELSE
                                    JS4 = JS4 - 2
                                    GOTO 36
                                 END IF
                              END IF
                           END DO
 38                        CONTINUE
                           ACC = PREF(KK)*ACC
                           TAL = AKPT(1)*S(2) + AKPT(2)*S(3)
                           CTAL = DCMPLX(COS(TAL),SIN(TAL))
                           ACC = ACC*CTAL
                           IF ( AC-EMACH.GT.0.D0 ) THEN
                              DLM(N) = DLM(N) + ACC/XPM(M)
                              IF ( M.EQ.1 ) GOTO 40
                           END IF
                           NM = N - M + 1
                           DLM(NM) = DLM(NM) + ACC*XPM(M)
 40                        CONTINUE
                           KK = KK + 1
                           N = N + 1
                        END DO
                     END DO
                     IF ( II.NE.1 ) THEN
                        II = 1
                        GOTO 34
C
C                 end of k//-loop
C
                     ELSE IF ( III.GT.0 ) THEN
                        EXIT
                     END IF
                  END DO
                  III = 0
               END DO
               TEST2 = 0.D0
               DO I = 1,NNDLM
                  DNORM = CDABS(DLM(I))
                  TEST2 = TEST2 + DNORM*DNORM
               END DO
               TEST = ABS((TEST2-TEST1)/TEST1)
               TEST1 = TEST2
               IF ( TEST-EMACH.GT.0.D0 ) THEN
                  IF ( N1.LT.20 ) GOTO 30
                  WRITE (NOUT1,99004) N1
               END IF
C
C         finally the elements of dlm are multiplied by the
C         factor (-1.0)**((m+|m|)/2)
C
               DO L = 2,LL2,2
                  N = L*L/2 + 1
                  DO M = 2,L,2
                     DLM(N) = -DLM(N)
                     N = N + 1
                  END DO
               END DO
               DO L = 2,LL2,2
                  N = L*L/2 + 1 + NSUDE
                  DO M = 1,L,2
                     DLM(N) = -DLM(N)
                     N = N + 1
                  END DO
               END DO
C
C         r = the current lattice vector in the sum
C         ar= mod(r)
C         icgl is a number for in- or excluding the cell at the origin
C
               PREF(1) = -0.5D0*KAPPA/RTPI
               ICGL = 0
            ELSE
               II = 1
               GOTO 20
            END IF
 60         CONTINUE
            ICGL = ICGL + 1
            N1 = 0
 80         CONTINUE
            N1 = N1 + 1
            NA = N1 + N1
            AN1 = DBLE(N1)
            AN2 = (-AN1) - 1.D0
            DO I1 = 1,NA
               AN2 = AN2 + 1.D0
               DO I2 = 1,4
                  AN = AN1
                  AN1 = -AN2
                  AN2 = AN
                  R(1) = AN1*AR1(1) + AN2*AR2(1)
                  R(2) = AN1*AR1(2) + AN2*AR2(2)
                  IF ( ICGL.LE.1 ) THEN
                     R(1) = 0.D0
                     R(2) = 0.D0
                  END IF
                  AD = AK(1)*R(1) + AK(2)*R(2)
                  SD = EXP((-AD*CI))
                  R(1) = R(1) + S(2)
                  R(2) = R(2) + S(3)
C
C        the arguments theta and fi of the spherical harmonics are
C        calculated : cff= exp(i*fi), ctt = cos(theta), stt = sin(theta)
C
                  C = R(1)*R(1) + R(2)*R(2)
                  AR = SQRT(C+S(1)*S(1))
                  B = 0.D0
                  CFF = CUNIT
                  IF ( C.GT.EMACH ) THEN
                     B = SQRT(C)
                     CFF = DCMPLX(R(1)/B,R(2)/B)
                  END IF
                  CTT = DCMPLX(S(1)/AR,0.D0)
                  STT = DCMPLX(B/AR,0.D0)
                  CALL SPHRKC(YLMA,CTT,STT,CFF,L2MAX)
C
C                 for each lattice vector the integral 'u' is obtained
C                 from the recurrence relation in l suggested by
C                 kambe. u1 and u2 are the initial terms of this
C                 recurrence, for l=-1 and l=0, and they are evaluated
C                 in terms of the complex error function zcerf
C
                  KANT = 0.5D0*AR*KAPPA
                  KNSQ = KANT*KANT
                  Z = CI*KANT/RTA
                  ZZ = RTA - Z
                  Z = RTA + Z
                  WW = ZCERF((-ZZ),EMACH)
                  W = ZCERF(Z,EMACH)
                  AA = 0.5D0*RTPI*(W-WW)/CI
                  AB = 0.5D0*RTPI*(W+WW)
                  A = ALPHA - KNSQ/ALPHA
                  XPA = EXP(A)
                  U1 = AA*XPA
                  U2 = AB*XPA/KANT
                  UODD1 = U2
C
C                 the contribution to dlm2 from a particular
C                 lattice vector is accumulated into the elements of dlm.
C                 this procedure includes the term (kant**l), and
C                 the recurrence for the integral 'u'
C                 dvel = 2**lq      lq is orbital impulsmoment
C                 dmel = (-1)**mq   mq is magnetic quantum number
C
                  II = 0
 85               CONTINUE
                  LST = 1 + II
                  AL = (-0.5D0) + DBLE(II)
                  DVEL = 0.5D0*(1.D0+DBLE(II))
                  CP = RTA
                  IF ( II.EQ.1 ) CP = CP/ALPHA
                  CF = CUNIT
                  IF ( II.EQ.1 ) CF = -2.D0*KANT*CF
                  DO L = LST,LL2
                     DVEL = DVEL*2.D0
                     MM = 1 + II
                     IF ( MOD(L,2).EQ.0 ) MM = 2 - II
                     N = ((L-II)*(L-II)+MM)/2 + NSUDE*II
                     DMEL = 3.D0 - 2.D0*DBLE(MM)
                     DO M = MM,L,2
                        NN = (L-1)*(L-1) + L - 1 + M
                        NN1 = NN - M - M + 2
                        ACC = PREF(1)*U2*CF*SD/DVEL*DMEL
                        DLM(N) = DLM(N) + ACC*YLMA(NN1)
                        IF ( M.NE.1 ) THEN
                           NM = N - M + 1
                           DLM(NM) = DLM(NM) + ACC*YLMA(NN)
                        END IF
                        N = N + 1
                     END DO
                     AL = AL + 1.D0
                     CP = CP/ALPHA
                     U = (AL*U2-U1+CP*XPA)/KNSQ
                     U1 = U2
                     IF ( L.EQ.1 ) UODD2 = U
                     U2 = U
                     CF = -2.D0*KANT*CF
                  END DO
                  IF ( II.NE.1 ) THEN
                     II = 1
                     U1 = UODD1
                     U2 = UODD2
                     GOTO 85
                  ELSE IF ( ICGL.EQ.1 ) THEN
                     GOTO 60
                  END IF
               END DO
            END DO
C
C         after each step of the summation a test on the
C         convergence of the elements of dlm is made
C
            TEST2 = 0.D0
            DO I = 1,NNDLM
               DNORM = CDABS(DLM(I))
               TEST2 = TEST2 + DNORM*DNORM
            END DO
            TEST = ABS((TEST2-TEST1)/TEST1)
            TEST1 = TEST2
            IF ( TEST-EMACH.GT.0.D0 ) THEN
               IF ( N1.LT.20 ) GOTO 80
               WRITE (NOUT1,99005) N1
            END IF
            IF ( IP.GE.5 ) THEN
               WRITE (NOUT1,99001) ITT,CEN
               WRITE (NOUT1,99002) NCOUNT
               WRITE (NOUT1,99003) (DLM(J),J=1,NNDLM)
C
C         summation over the clebsch-gordon type coefficients
C         elm proceeds, first for even, and then for odd l+m.
C         this gives kambe's elements  xalm(l2,m2;l3,m3)  which,
C         combined with the phase shift terms in af, give the
C         elements  xalm(l3,m3;l2,m2)
C
            END IF
            K = 1
            II = 0
 100        CONTINUE
            L30 = LMAX + II
            I = 1
            DO IL2 = 1,L30
               L2 = IL2 - II
               M2 = (-L2) + 1 - II
               DO I2 = 1,IL2
                  J = 1
                  DO IL3 = 1,L30
                     L3 = IL3 - II
                     M3 = (-L3) + 1 - II
                     DO I3 = 1,IL3
                        ALM = CZERO
                        LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                        LB1 = L2 + L3
                        N = (LB1*(LB1+2)+M2-M3+2)/2
                        NN = 2*LB1
                        LB11 = LB1 + 1
                        LA11 = LA1 + 1
                        DO L1 = LA11,LB11,2
                           ALM = ALM + ELM(K)*DLM(N)
                           N = N - NN
                           NN = NN - 4
                           K = K + 1
                        END DO
                        ALM = ALM/KAPPA
                        IF ( II.LE.0 ) THEN
                           JOTWO = J + (JCOUNT-1)*NDIM + MLR
                           IOTWO = I + (ICOUNT-1)*NDIM + MLR
                           LOTWO = L3 + 1 + (JCOUNT-1)*LAF1
                           XALM(JOTWO,IOTWO) = -AF(LOTWO)*ALM
                        ELSE
                           JETWO = J + (JCOUNT-1)*NDIM
                           IETWO = I + (ICOUNT-1)*NDIM
                           LOTWO = L3 + 1 + (JCOUNT-1)*LAF1
                           XALM(JETWO,IETWO) = -AF(LOTWO)*ALM
                        END IF
                        M3 = M3 + 2
                        J = J + 1
                     END DO
                  END DO
                  M2 = M2 + 2
                  I = I + 1
               END DO
            END DO
            IF ( II.LE.0 ) THEN
               II = 1
               GOTO 100
            ELSE IF ( ABS(S(1)).GE.EMACH ) THEN
               II = 0
               II1 = 1
 110           CONTINUE
               L30 = LMAX + II
               LL1 = LMAX + II1
               I = 1
               DO IL2 = 1,L30
                  L2 = IL2 - II
                  M2 = (-L2) + 1 - II
                  DO I2 = 1,IL2
                     J = 1
                     DO IL3 = 1,LL1
                        L3 = IL3 - II1
                        M3 = (-L3) + 1 - II1
                        DO I3 = 1,IL3
                           ALM = CZERO
                           LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                           LB1 = L2 + L3
                           LB11 = LB1 + 1
                           LA11 = LA1 + 1
                           NN = LA11 + LB11 - 1
                           M1 = M2 - M3
                           DO L1 = LA11,LB11,2
                              L11 = NN - L1
                              N = NSUDE + (L11*L11+M1+1)/2
                              ALM = ALM + ELM(K)*DLM(N)
                              K = K + 1
                           END DO
                           ALM = ALM/KAPPA
                           IF ( II.LE.0 ) THEN
                              JOTWO = J + (JCOUNT-1)*NDIM
                              IOTWO = I + (ICOUNT-1)*NDIM + MLR
                              LOTWO = L3 + 1 + (JCOUNT-1)*LAF1
                              XALM(JOTWO,IOTWO) = -AF(LOTWO)*ALM
                           ELSE
                              JETWO = J + (JCOUNT-1)*NDIM + MLR
                              IETWO = I + (ICOUNT-1)*NDIM
                              LOTWO = L3 + 1 + (JCOUNT-1)*LAF1
                              XALM(JETWO,IETWO) = -AF(LOTWO)*ALM
                           END IF
                           M3 = M3 + 2
                           J = J + 1
                        END DO
                     END DO
                     M2 = M2 + 2
                     I = I + 1
                  END DO
               END DO
               IF ( II.LE.0 ) THEN
                  II = 1
                  II1 = 0
                  GOTO 110
               END IF
            END IF
         END DO
C
C     transformation from even1,odd1,even2,odd2,... order
C     to odd1,odd2...even1,even2.... order
C
      END IF
      CALL TRANS(XALM,NALM,MLR,MLT,MLQ,NATIN)
      DO I = 1,NOD
         DO J = 1,NOD
            GODD(I,J) = XALM(I,J)
         END DO
      END DO
      DO I = 1,NEV
         DO J = 1,NEV
            GEVEN(I,J) = XALM(I+NOD,J+NOD)
         END DO
      END DO
      DO I = 1,NOD
         DO J = 1,NEV
            GODEV(I,J) = XALM(I,J+NOD)
         END DO
      END DO
      DO I = 1,NEV
         DO J = 1,NOD
            GEVOD(I,J) = XALM(I+NOD,J)
         END DO
      END DO
      RETURN
99001 FORMAT (' dmatr : layer type it=',i2,' at energy e= ',2F12.4)
99002 FORMAT (' columm : ',i2)
99003 FORMAT (6E12.4)
99004 FORMAT (' ** dlm1 not converged in xmatr by n1 =',i3)
99005 FORMAT (' ** dlm2 not converged in xmatr by n1 =',i3)
      END
C*==akmj0.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE AKMJ0(MAXL,E,GANZ,CLIGHT,IREL,UGSP,UGSM,IPOL,IVOR,IAN,
     &                 YP1,YP2,YS1,YS2,EPOSM,EPOSP,IT,IL,KGT,DIMA,DIMB,
     &                 DIMC)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C sphrm4                                                      *
C**************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQ,LL,LG,LH,CZERO,CONE,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLY
      PARAMETER (MLY=2*ML+1)
C
C Dummy arguments
C
      REAL*8 CLIGHT
      INTEGER DIMA,DIMB,DIMC,GANZ,IAN,IL,IPOL,IREL,IT,IVOR,MAXL
      COMPLEX*16 E
      COMPLEX*16 EPOSM(NATLM,LH,LAYSM),EPOSP(NATLM,LH,LAYSM),
     &           UGSM(DIMA,DIMB),UGSP(DIMA,DIMB),YP1(ML,MLY,DIMC,NATLM),
     &           YP2(ML,MLY,DIMC,NATLM),YS1(ML,MLY,DIMC,NATLM),
     &           YS2(ML,MLY,DIMC,NATLM)
      REAL*8 KGT(2,LG)
C
C Local variables
C
      COMPLEX*16 CI,IC(:),K0BTRS,KZ,KZZ(2),YLM(:)
      INTEGER I1G,I2G,IA,IG,IKZ,J,L1,M,M1
      REAL*8 KXX,KYY
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IC,YLM
      ALLOCATE (IC(ML),YLM(MLQ))
C
C*** End of declarations rewritten by SPAG
C
      K0BTRS = 2.D0*E
      IF ( IREL.EQ.2 ) K0BTRS = 2.D0*E + E*E/(CLIGHT*CLIGHT)
      CI = CIMAG
      IC(1) = CONE
      DO IA = 1,MAXL
         IC(IA+1) = IC(IA)*CI
      END DO
      DO IG = 1,GANZ
         IF ( IREL.EQ.1 ) THEN
            I1G = IG
            I2G = IG
         ELSE
            I1G = 2*IG
            I2G = 2*IG - 1
         END IF
         KXX = KGT(1,IG)
         KYY = KGT(2,IG)
         KZZ(1) = (-1.D0)*CDSQRT(K0BTRS-DCMPLX(KXX*KXX+KYY*KYY,0.D0))
         KZZ(2) = -KZZ(1)
         DO IKZ = 1,2
            KZ = KZZ(IKZ)
            CALL SPHRM4(KXX,KYY,KZ,YLM)
            J = 1
            DO L1 = 1,ML
               M1 = 2*L1 - 1
               DO M = 1,M1
                  IF ( IKZ.NE.1 ) THEN
                     IF ( IPOL.EQ.1 ) THEN
                        IF ( IG.EQ.1 ) YP1(L1,M,IL,IAN) = CZERO
                     ELSE
                        IF ( IG.EQ.1 ) YP2(L1,M,IL,IAN) = CZERO
                     END IF
                     IF ( IVOR.EQ.1 ) THEN
                        IF ( IPOL.EQ.1 ) THEN
                           YP1(L1,M,IL,IAN) = YP1(L1,M,IL,IAN)
     &                        + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSP(IG,IL)
     &                        *(-1)**(IREL-1)*EPOSP(IAN,I2G,IT)
                        ELSE
                           YP2(L1,M,IL,IAN) = YP2(L1,M,IL,IAN)
     &                        + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSP(IG,IL)
     &                        *(-1)**(IREL-1)*EPOSP(IAN,I2G,IT)
                        END IF
                     ELSE IF ( IPOL.EQ.1 ) THEN
                        YP1(L1,M,IL,IAN) = YP1(L1,M,IL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSP(IG,IL)
     &                     *EPOSP(IAN,I2G,IT)
                     ELSE
                        YP2(L1,M,IL,IAN) = YP2(L1,M,IL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSP(IG,IL)
     &                     *EPOSP(IAN,I2G,IT)
                     END IF
                     IF ( IREL.NE.1 ) THEN
                        IF ( IPOL.EQ.1 ) THEN
                           IF ( IG.EQ.1 ) YP1(L1,M,IL+LL,IAN) = CZERO
                           YP1(L1,M,IL+LL,IAN) = YP1(L1,M,IL+LL,IAN)
     &                        + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                        *UGSP(IG+GANZ,IL)*EPOSP(IAN,I1G,IT)
                        ELSE
                           IF ( IG.EQ.1 ) YP2(L1,M,IL+LL,IAN) = CZERO
                           YP2(L1,M,IL+LL,IAN) = YP2(L1,M,IL+LL,IAN)
     &                        + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                        *UGSP(IG+GANZ,IL)*EPOSP(IAN,I1G,IT)
                        END IF
                     END IF
                     GOTO 5
                  END IF
                  IF ( IPOL.EQ.1 ) THEN
                     IF ( IG.EQ.1 ) YS1(L1,M,IL,IAN) = CZERO
                  ELSE
                     IF ( IG.EQ.1 ) YS2(L1,M,IL,IAN) = CZERO
                  END IF
                  IF ( IVOR.EQ.1 ) THEN
                     IF ( IPOL.EQ.1 ) THEN
                        YS1(L1,M,IL,IAN) = YS1(L1,M,IL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSM(IG,IL)
     &                     *(-1)**(IREL-1)*EPOSM(IAN,I2G,IT)
                     ELSE
                        YS2(L1,M,IL,IAN) = YS2(L1,M,IL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)*UGSM(IG,IL)
     &                     *(-1)**(IREL-1)*EPOSM(IAN,I2G,IT)
                     END IF
                  ELSE IF ( IPOL.EQ.1 ) THEN
                     YS1(L1,M,IL,IAN) = YS1(L1,M,IL,IAN)
     &                                  + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                                  *UGSM(IG,IL)*EPOSM(IAN,I2G,IT)
                  ELSE
                     YS2(L1,M,IL,IAN) = YS2(L1,M,IL,IAN)
     &                                  + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                                  *UGSM(IG,IL)*EPOSM(IAN,I2G,IT)
                  END IF
                  IF ( IREL.NE.1 ) THEN
                     IF ( IPOL.EQ.1 ) THEN
                        IF ( IG.EQ.1 ) YS1(L1,M,IL+LL,IAN) = CZERO
                        YS1(L1,M,IL+LL,IAN) = YS1(L1,M,IL+LL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                     *UGSM(IG+GANZ,IL)*EPOSM(IAN,I1G,IT)
                     ELSE
                        IF ( IG.EQ.1 ) YS2(L1,M,IL+LL,IAN) = CZERO
                        YS2(L1,M,IL+LL,IAN) = YS2(L1,M,IL+LL,IAN)
     &                     + YLM(J+M1-2*M+1)*(-1)**(J+1)
     &                     *UGSM(IG+GANZ,IL)*EPOSM(IAN,I1G,IT)
                     END IF
                  END IF
 5                CONTINUE
                  J = J + 1
               END DO
            END DO
         END DO
      END DO
C
      END
C*==akmji.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE AKMJI(K,MAXL,IRO,IREL,XED,XOD,IAN,NAT,IED,IOD,IEOD,
     &                 IOED,JLA)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C clegor                                                      *
C**************************************************************
C
      USE MOD_SPEC_COE1,ONLY:LSPOO,IIOO,POOO,POSOO,IGZOO,IROOO,IAOO
      USE MOD_SPEC_COE2,ONLY:IAOOS,LSPEE,IIEE,PEEE,PESEE,IGZEE,IROEE,
     &    IAEE
      USE MOD_SPEC_COE3,ONLY:IAEES,LSPOE,IIOE,POOE,PESOE,IGZOE,IROOE,
     &    IAOE
      USE MOD_SPEC_COE4,ONLY:IAOES,LSPEO,IIEO,PEEO,POSEO,IGZEO,IROEO,
     &    IAEO
      USE MOD_SPEC_COE5,ONLY:IAEOS,COO,CEE,COE,CEO
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,KME,KMO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IAN,IED,IEOD,IOD,IOED,IREL,IRO,JLA,K,MAXL,NAT,XED,XOD
C
C Local variables
C
      REAL*8 C,KREAL,MU,MUS
      REAL*8 CLEGOR
      INTEGER IANS,IGZ,II,II1,II2,INDE1,INDE2,INDES1,INDES2,INDO1,INDO2,
     &        INDOS1,INDOS2,ISIS,IU,IZ,JU,KAP1,KAP2,KS,LD,LKS,LS,MAXL1,
     &        MM,MMS,NSIS,PE,PES,PO,POS,XED1,XOD1
C
C*** End of declarations rewritten by SPAG
C
      MAXL1 = MAXL + 1
C
C     IF ( K.LT.0 ) THEN
C        LD = ((-K)) - 1
C        MUMIN = K + 0.5D0
C     ELSE
C        !  LD = ?
C        MUMIN = ((-K)) + 0.5D0
C     END IF
C
      IF ( K.LT.0 ) THEN
         LD = ((-K)) - 1
      ELSE
         LD = K
      END IF
C
      IF ( IREL.EQ.1 ) THEN
         KAP1 = 2*K + 1
      ELSE IF ( IREL.EQ.2 ) THEN
         KAP1 = 2*ABS(K)
      END IF
C
C     MUMAX = -MUMIN
C     MUEND = MUMAX
C     IF ( IREL.EQ.1 ) MUEND = MUMAX - 1.D0
C     DO MU = MUMIN,MUEND
C
      IGZ = 1
      DO II1 = 1,KAP1
         IF ( IREL.EQ.1 ) THEN
            MU = (-ABS(K)+II1-1)*1.0D0
         ELSE IF ( IREL.EQ.2 ) THEN
            MU = -ABS(K) + 0.5D0 + (II1*1.D0-1.0D0)
         END IF
C
         MM = IDINT(MU+0.5D0)
         PE = 0
         PO = 0
         XED1 = XED/NAT
         INDE1 = 1 + (IAN-1)*XED1
         INDE2 = INDE1 + XED1 - 1
         XOD1 = XOD/NAT
         INDO1 = 1 + (IAN-1)*XOD1
         INDO2 = INDO1 + XOD1 - 1
         IF ( IREL.EQ.2 ) THEN
            DO IU = INDE1,INDE2
               IF ( K-IDINT(KME(IU,1)).EQ.0 .AND. IDINT(MU-KME(IU,2))
     &              .EQ.0 ) PE = IU
               IF ( K-IDINT(KMO(IU,1)).EQ.0 .AND. IDINT(MU-KMO(IU,2))
     &              .EQ.0 ) PO = IU
            END DO
         ELSE
            DO IU = INDE1,INDE2
               IF ( LD-LME(IU,1).EQ.0 .AND. MM.EQ.LME(IU,2) ) PE = IU
            END DO
            DO IU = INDO1,INDO2
               IF ( LD-LMO(IU,1).EQ.0 .AND. MM.EQ.LMO(IU,2) ) PO = IU
            END DO
         END IF
         DO IANS = 1,NAT
            INDES1 = 1 + (IANS-1)*XED1
            INDES2 = INDES1 + XED1 - 1
            INDOS1 = 1 + (IANS-1)*XOD1
            INDOS2 = INDOS1 + XOD1 - 1
            DO LKS = 1,MAXL1
               DO NSIS = 1,IREL
                  ISIS = (-1)*(-1)**NSIS
                  KS = ISIS*(LKS-1) + (ISIS-1)/2
                  KREAL = DBLE(KS)
                  IF ( IREL.EQ.1 ) KS = -LKS
                  IF ( KS.NE.0 ) THEN
C
C                    IF ( KS.LT.0 ) THEN
C                       MUMINS = KS + 0.5D0
C                       LS = ((-KS)) - 1
C                    ELSE
C                       MUMINS = ((-KS)) + 0.5D0
C                       LS = KS
C                    END IF
C                    MUMAXS = -MUMINS
C                    MUSAN = MUMINS + 2.D0 - DBLE(IREL)
C
                     IF ( KS.LT.0 ) THEN
                        LS = ((-KS)) - 1
                     ELSE
                        LS = KS
                     END IF
C
                     IF ( IREL.EQ.1 ) THEN
                        KAP2 = 2*KS + 1
                     ELSE IF ( IREL.EQ.2 ) THEN
                        KAP2 = 2*ABS(KS)
                     END IF
C
C                    DO MUS = MUSAN,MUMAXS
C
                     IZ = 3 - IREL
                     DO II2 = 1,KAP2
                        IF ( IREL.EQ.1 ) THEN
                           MUS = (-ABS(KS)+II2-1)*1.0D0
                        ELSE IF ( IREL.EQ.2 ) THEN
                           MUS = -ABS(KS) + 0.5D0 + (II2*1.D0-1.0D0)
                        END IF
C
                        C = 1.D0
                        IF ( IRO.EQ.1 .AND. IREL.EQ.2 )
     &                       C = CLEGOR(KREAL,MUS,0.5D0)
                        IF ( IRO.EQ.2 .AND. IREL.EQ.2 )
     &                       C = CLEGOR(KREAL,MUS,-0.5D0)
                        MMS = INT(MUS-0.5D0)
                        PES = 0
                        POS = 0
                        IF ( IREL.EQ.2 ) THEN
                           DO JU = INDES1,INDES2
                              IF ( KS-IDINT(KME(JU,1)).EQ.0 .AND. 
     &                             INT(MUS-KME(JU,2)).EQ.0 ) PES = JU
                              IF ( KS-IDINT(KMO(JU,1)).EQ.0 .AND. 
     &                             INT(MUS-KMO(JU,2)).EQ.0 ) POS = JU
                           END DO
                        ELSE
                           DO JU = INDES1,INDES2
                              IF ( LS-LME(JU,1).EQ.0 .AND. 
     &                             MMS.EQ.LME(JU,2) ) PES = JU
                           END DO
                           DO JU = INDOS1,INDOS2
                              IF ( LS-LMO(JU,1).EQ.0 .AND. 
     &                             MMS.EQ.LMO(JU,2) ) POS = JU
                           END DO
                        END IF
                        II = IZ + IRO - 1
                        IF ( KS.LT.0 ) II = IZ + IRO - 2
                        IF ( KS.GE.0 .OR. II.NE.2*LS+2 .OR. IRO.NE.2 )
     &                       THEN
                           IF ( II.NE.0 ) THEN
                              IF ( LS.EQ.0 ) II = 1
                              IF ( PE.NE.0 .AND. PES.NE.0 ) THEN
                                 LSPEE(IED) = LS + 1
                                 IIEE(IED) = II
                                 PEEE(IED) = PE
                                 PESEE(IED) = PES
                                 CEE(IED) = C
                                 IGZEE(IED) = IGZ + JLA
                                 IROEE(IED) = IRO
                                 IAEE(IED) = IAN
                                 IAEES(IED) = IANS
                                 IED = IED + 1
                              ELSE IF ( PO.NE.0 .AND. POS.NE.0 ) THEN
                                 LSPOO(IOD) = LS + 1
                                 IIOO(IOD) = II
                                 POOO(IOD) = PO
                                 POSOO(IOD) = POS
                                 COO(IOD) = C
                                 IGZOO(IOD) = IGZ + JLA
                                 IROOO(IOD) = IRO
                                 IAOO(IOD) = IAN
                                 IAOOS(IOD) = IANS
                                 IOD = IOD + 1
                              ELSE IF ( PE.NE.0 .AND. POS.NE.0 ) THEN
                                 LSPEO(IEOD) = LS + 1
                                 IIEO(IEOD) = II
                                 PEEO(IEOD) = PE
                                 POSEO(IEOD) = POS
                                 CEO(IEOD) = C
                                 IGZEO(IEOD) = IGZ + JLA
                                 IROEO(IEOD) = IRO
                                 IAEO(IEOD) = IAN
                                 IAEOS(IEOD) = IANS
                                 IEOD = IEOD + 1
                              ELSE IF ( PO.NE.0 .AND. PES.NE.0 ) THEN
                                 LSPOE(IOED) = LS + 1
                                 IIOE(IOED) = II
                                 POOE(IOED) = PO
                                 PESOE(IOED) = PES
                                 COE(IOED) = C
                                 IGZOE(IOED) = IGZ + JLA
                                 IROOE(IOED) = IRO
                                 IAOE(IOED) = IAN
                                 IAOES(IOED) = IANS
                                 IOED = IOED + 1
                              END IF
                           END IF
                           IZ = IZ + 1
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END DO
         IGZ = IGZ + 1
      END DO
      END
C*==atphas.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE ATPHAS(POS,EPOSLM,EPOSLP,EPOSHM,EPOSHP,E,GANZ,LHL,IREL,
     &                  CLIGHT,IT,NATL,KGT)
      USE MOD_SPEC,ONLY:LAYSM,NATLM,LG,LH,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT
      COMPLEX*16 E
      INTEGER GANZ,IREL,IT,LHL
      COMPLEX*16 EPOSHM(NATLM,LH,LAYSM),EPOSHP(NATLM,LH,LAYSM),
     &           EPOSLM(NATLM,LH,LAYSM),EPOSLP(NATLM,LH,LAYSM)
      REAL*8 KGT(2,LG),POS(3,NATLM,LAYSM)
      INTEGER NATL(LAYSM)
C
C Local variables
C
      COMPLEX*16 CI,K0BTRS,KGZ
      INTEGER IAN,IG,IM,L,NAT
C
C*** End of declarations rewritten by SPAG
C
      CI = CIMAG
      K0BTRS = 2.D0*E
      IF ( IREL.EQ.2 ) K0BTRS = 2.D0*E + E*E/(CLIGHT*CLIGHT)
      NAT = NATL(IT)
      DO IAN = 1,NAT
         IG = 1
         DO L = 1,GANZ
            KGZ = CDSQRT(K0BTRS-(KGT(1,L)**2+KGT(2,L)**2))
            DO IM = 1,IREL
               IF ( LHL.EQ.1 ) THEN
                  EPOSHM(IAN,IG,IT) = CDEXP(CI*(KGT(1,L)*POS(2,IAN,IT)+
     &                                KGT(2,L)*POS(3,IAN,IT)
     &                                -KGZ*POS(1,IAN,IT)))
                  EPOSHP(IAN,IG,IT) = CDEXP(CI*(KGT(1,L)*POS(2,IAN,IT)+
     &                                KGT(2,L)*POS(3,IAN,IT)
     &                                +KGZ*POS(1,IAN,IT)))
               ELSE IF ( LHL.EQ.2 ) THEN
                  EPOSLM(IAN,IG,IT) = CDEXP(CI*(KGT(1,L)*POS(2,IAN,IT)+
     &                                KGT(2,L)*POS(3,IAN,IT)
     &                                -KGZ*POS(1,IAN,IT)))
                  EPOSLP(IAN,IG,IT) = CDEXP(CI*(KGT(1,L)*POS(2,IAN,IT)+
     &                                KGT(2,L)*POS(3,IAN,IT)
     &                                +KGZ*POS(1,IAN,IT)))
               END IF
               IG = IG + 1
            END DO
         END DO
      END DO
      END
C*==blm.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      REAL*8 FUNCTION BLM(IL1,IM1,IL2,IM2,IL3,IM3)
C     /****************************************************************/
C     # purpose       : calculate gaunt coefficients                   *
C                                                                      *
C     # function called from this routine                              *
C       fac                                                            *
C     /****************************************************************/
      USE MOD_SPEC,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IL1,IL2,IL3,IM1,IM2,IM3
C
C Local variables
C
      REAL*8 A,AA
      REAL*8 FAC
      INTEGER I,I1,I2,I3,I4,I5,I6,IS,IT,L,L1,L2,L3,M,M1,M2,M3
      EXTERNAL FAC
C
C*** End of declarations rewritten by SPAG
C
      BLM = 0.D0
      L1 = IL1
      L2 = IL2
      L3 = IL3
      M1 = IABS(IM1)
      M2 = IABS(IM2)
      M3 = IABS(IM3)
      IF ( L1.LT.M1 .OR. L2.LT.M2 .OR. L3.LT.M3 ) RETURN
C
      M = IABS(IM1+IM2+IM3) + MOD(L1+L2+L3,2)
      IF ( M.NE.0 .OR. L1.GT.L2+L3 .OR. L1.LT.IABS(L2-L3) ) RETURN
      L = L1
      M = M1
      M1 = MAX0(M1,M2,M3)
      IF ( M1.NE.M ) THEN
         IF ( M1.NE.M2 ) THEN
            IF ( M1.EQ.M3 ) THEN
               L1 = L3
               L3 = L
               M3 = M
               GOTO 100
            END IF
         END IF
         L1 = L2
         L2 = L
         M2 = M
      END IF
 100  CONTINUE
      IF ( L2.LT.L3 ) THEN
         L = L2
         M = M2
         L2 = L3
         M2 = M3
         L3 = L
         M3 = M
      END IF
C
      IS = (L1+L2+L3)/2
C
      BLM = (-1.D0)**(IS-L2-M3+M1)
     &      *SQRT(DBLE((2*L1+1)*(2*L2+1)*(2*L3+1))/PI)/DBLE(2*(2*IS+1))
     &      *SQRT(FAC(L1+M1,M1+M1))/FAC(IS-L1,IS-L1)
     &      *SQRT(FAC(L2+M2,M2+M2))/FAC(IS-L2,IS-L2)
     &      *SQRT(FAC(L3+M3,M3+M3))
     &      *FAC(MAX0(L2+L3-M1,L2-L3+M1),IABS(2*(L3-M1)))
     &      **ISIGN(1,L3-M1)
      IF ( L3.EQ.0 ) RETURN
      DO I = 1,L3
         BLM = BLM/DBLE(2*(2*(IS-I)+1))
      END DO
      A = 1.D0
      AA = 1.D0
      I1 = L1 + M1
      I2 = MAX0(1,L2+L3-M1)
      I3 = 0
      I4 = L1 - M1
      I5 = L2 - L3 + M1
      I6 = L3 - M3
      IT = MIN0(I4,I6)
      IF ( IT.EQ.0 ) RETURN
      DO I = 1,IT
         I1 = I1 + 1
         I3 = I3 + 1
         I5 = I5 + 1
         AA = -AA*DBLE(I1*I4*I6)/DBLE(I2*I3*I5)
         A = A + AA
         I2 = I2 - 1
         I4 = I4 - 1
         I6 = I6 - 1
      END DO
      BLM = A*BLM
C
      END
C*==block.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE BLOCK(X,NX,MX,IS,JS,IFLM,IOE,DLM1,CELM,NCLM)
C
      USE MOD_SPEC,ONLY:NATLM,ML,MLZP,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LMAXM,L2MAXM,NFM,NDLMM
      PARAMETER (LMAXM=ML-1,L2MAXM=2*LMAXM,NFM=NATLM*NATLM-NATLM+1,
     &           NDLMM=(L2MAXM+1)*(L2MAXM+1))
C
C Dummy arguments
C
      INTEGER IFLM,IOE,IS,JS,MX,NCLM,NX
      REAL*8 CELM(MLZP)
      COMPLEX*16 DLM1(NDLMM,NFM),X(NX,MX)
C
C Local variables
C
      COMPLEX*16 ACC
      INTEGER I,I1,I2,IL1,IL2,IL3,IX,J,JX,K,L1,L2,L32,LA3,LB3,LKX,LM,
     &        LMAX,M1,M2,M3
C
C*** End of declarations rewritten by SPAG
C
      K = 0
      LMAX = ML - 1
      IF ( IOE.EQ.1 ) K = NCLM
      LKX = LMAX + IOE
      I = 0
      DO IL2 = 1,LKX
         L2 = IL2 - IOE
         M2 = ((-L2)) - 1 - IOE
         DO I2 = 1,IL2
            M2 = M2 + 2
            I = I + 1
            J = 0
            DO IL1 = 1,LKX
               L1 = IL1 - IOE
               M1 = ((-L1)) - 1 - IOE
               DO I1 = 1,IL1
                  M1 = M1 + 2
                  J = J + 1
                  M3 = M1 - M2
                  LA3 = MAX0(IABS(L2-L1),IABS(M3))
                  LB3 = L2 + L1
                  L32 = 2*LB3
                  LM = (LB3*(LB3+2)+2-M3)/2
                  ACC = CZERO
                  DO IL3 = LA3,LB3,2
                     K = K + 1
                     ACC = ACC + CELM(K)*DLM1(LM,IFLM)
                     LM = LM - L32
                     L32 = L32 - 4
                  END DO
                  JX = J + IS
                  IX = I + JS
                  X(JX,IX) = ACC
               END DO
            END DO
         END DO
      END DO
      END
C*==celmg.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CELMG(CELM,YLM,FAC2,FAC1,ICE,NCLM)
C     *********************************************************
C     subroutines and functions called from this routine      *
C     blm                                                     *
C     *********************************************************
C
      USE MOD_SPEC,ONLY:ML,MLZP,PI
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLA,MLB
      PARAMETER (MLA=2*ML-1,MLB=MLA*MLA)
C
C Dummy arguments
C
      INTEGER ICE,NCLM
      REAL*8 CELM(MLZP),FAC1(MLA),FAC2(MLB),YLM(MLB)
C
C Local variables
C
      REAL*8 A,ALM,ASG,B,CL,CM
      REAL*8 BLM
      INTEGER I2,I3,II,IL2,IL3,K,L,L1,L11,L2,L2MAX,L3,L5,LA1,LA11,LB1,
     &        LB11,LM,LM2,LM3,LMAX,LN,LO,LP,LQ,M,M1,M2,M3,NFAC
C
C*** End of declarations rewritten by SPAG
C
      LMAX = ML - 1
      NFAC = 4*LMAX + 2
      IF ( NFAC.GE.50 ) THEN
         WRITE (6,99001) LMAX
         STOP
      END IF
      L2MAX = LMAX + LMAX
C
C       the array ylm is first loaded with spherical
C       harmonics, arguments theta=pi/2.0, fi=0.d0
C
      LM = 0
      CL = 0.D0
      A = 1.D0
      B = 1.D0
      ASG = 1.D0
      L1 = L2MAX + 1
C
C     multiplicative factors required
C
      DO L = 1,L1
         FAC1(L) = ASG*SQRT((2.D0*CL+1.D0)*A/(4.D0*PI*B*B))
         CM = -CL
         LN = L + L - 1
         DO M = 1,LN
            LO = LM + M
            FAC2(LO) = SQRT((CL+1.D0+CM)*(CL+1.D0-CM)/((2.D0*CL+3.D0)*(
     &                 2.D0*CL+1.D0)))
            CM = CM + 1.D0
         END DO
         CL = CL + 1.D0
         A = A*2.D0*CL*(2.D0*CL-1.D0)/4.D0
         B = B*CL
         ASG = -ASG
         LM = LM + LN
      END DO
C
C       first all the ylm for m=+-l and m=+-(l-1) are
C       are calculated by explicit formulae
C
      LM = 1
      CL = 1.D0
      ASG = -1.D0
      YLM(1) = FAC1(1)
      DO L = 1,L2MAX
         LN = LM + L + L + 1
         YLM(LN) = FAC1(L+1)
         YLM(LM+1) = ASG*FAC1(L+1)
         YLM(LN-1) = 0.D0
         YLM(LM+2) = 0.D0
         CL = CL + 1.D0
         ASG = -ASG
         LM = LN
      END DO
C
C       using ylm and yl(m-1) in a recurrence relation
C       yl(m+1) is calculated
C
      LM = 1
      L1 = L2MAX - 1
      DO L = 1,L1
         LN = L + L - 1
         LM2 = LM + LN + 4
         LM3 = LM - LN
         DO M = 1,LN
            LO = LM2 + M
            LP = LM3 + M
            LQ = LM + M + 1
            YLM(LO) = -(FAC2(LP)*YLM(LP))/FAC2(LQ)
         END DO
         LM = LM + L + L + 1
      END DO
      IF ( IP.GT.4 .AND. ICE.EQ.1 ) WRITE (NOUT1,99002) LMAX
      IF ( IP.GT.4 .AND. ICE.EQ.2 ) WRITE (NOUT1,99003) LMAX
      IF ( IP.GT.4 ) WRITE (NOUT1,99004)
      K = 1
      II = 0
 100  CONTINUE
      L1 = LMAX + II
      DO IL2 = 1,L1
         L2 = IL2 - II
         M2 = (-L2) + 1 - II
         DO I2 = 1,IL2
            DO IL3 = 1,L1
               L3 = IL3 - II
               M3 = (-L3) + 1 - II
               DO I3 = 1,IL3
                  LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                  LB1 = L2 + L3
                  LA11 = LA1 + 1
                  LB11 = LB1 + 1
                  M1 = M3 - M2
                  DO L11 = LA11,LB11,2
                     L5 = LA11 + LB11 - L11 - 1
                     L = (L2-L5-L3)/2 + M1 + M3
                     M = L5*(L5+1) - M1 + 1
                     ALM = (-1.D0)**L*4.D0*PI*BLM(L5,M1,L2,M2,L3,(-M3))
                     CELM(K) = ALM
C
                     IF ( ICE.EQ.1 ) CELM(K) = YLM(M)*ALM
                     IF ( IP.GT.4 ) WRITE (NOUT1,99005) K,L5,M1,L2,M2,
     &                    L3,M3,CELM(K)
C
                     K = K + 1
                  END DO
                  M3 = M3 + 2
               END DO
            END DO
            M2 = M2 + 2
         END DO
      END DO
      IF ( II.GT.0 ) THEN
         RETURN
      ELSE
         II = 1
         NCLM = K - 1
         GOTO 100
      END IF
99001 FORMAT ('0*** lmax =',i3,' too big to be handled by blm',
     &        ' routine')
99002 FORMAT ('1clm''s for xmat',5x,'lmax =',i3)
99003 FORMAT ('1elm''s for xmatk',5x,'lmax =',i3)
99004 FORMAT ('0',2x,'k',6x,'l1 m1 l2 m2 l3 m3',/)
99005 FORMAT (' ',i5,5x,6I3,5x,e15.5)
      END
C*==celmr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CELMR(CELM,YLM,FAC2,FAC1,ICE)
C**************************************************************
C                                                             *
C subroutines and functions called from this routine          *
C                                                             *
C blm                                                         *
C                                                             *
C**************************************************************
C
      USE MOD_SPEC,ONLY:ML,MLZP,PI
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLA,MLB
      PARAMETER (MLA=2*ML-1,MLB=MLA*MLA)
C
C Dummy arguments
C
      INTEGER ICE
      REAL*8 CELM(MLZP),FAC1(MLA),FAC2(MLB),YLM(MLB)
C
C Local variables
C
      REAL*8 A,ALM,ASG,B,CL,CM
      REAL*8 BLM
      INTEGER I2,I3,II,II1,IL2,IL3,K,L,L1,L11,L2,L2MAX,L3,LA1,LA11,LB1,
     &        LB11,LL1,LLTVAR,LM,LM2,LM3,LMAX,LN,LO,LP,LQ,M,M1,M2,M3,
     &        NFAC
C
C*** End of declarations rewritten by SPAG
C
      LMAX = ML - 1
      NFAC = 4*LMAX + 2
      IF ( NFAC.GE.50 ) THEN
         WRITE (NOUT1,99001) LMAX
         STOP
      END IF
      L2MAX = LMAX + LMAX
      LM = 0
      CL = 0.D0
      A = 1.D0
      B = 1.D0
      ASG = 1.D0
      LLTVAR = L2MAX + 1
      DO L = 1,LLTVAR
         FAC1(L) = ASG*SQRT((2.D0*CL+1.D0)*A/(4.D0*PI*B*B))
         CM = -CL
         LN = L + L - 1
         DO M = 1,LN
            LO = LM + M
            FAC2(LO) = SQRT((CL+1.D0+CM)*(CL+1.D0-CM)/((2.D0*CL+3.D0)*(
     &                 2.D0*CL+1.D0)))
            CM = CM + 1.D0
         END DO
         CL = CL + 1.D0
         A = A*2.D0*CL*(2.D0*CL-1.D0)/4.D0
         B = B*CL
         ASG = -ASG
         LM = LM + LN
      END DO
C
C      first all the ylm for m=+-l and m=+-(l-1) are
C      are calculated by explicit formulae
C
      LM = 1
      CL = 1.D0
      ASG = -1.D0
      YLM(1) = FAC1(1)
      DO L = 1,L2MAX
         LN = LM + L + L + 1
         YLM(LN) = FAC1(L+1)
         YLM(LM+1) = ASG*FAC1(L+1)
         YLM(LN-1) = 0.D0
         YLM(LM+2) = 0.D0
         CL = CL + 1.D0
         ASG = -ASG
         LM = LN
      END DO
C
C       using ylm and yl(m-1) in a recurrence relation
C       yl(m+1) is calculated
C
      LM = 1
      LLTVAR = L2MAX - 1
      DO L = 1,LLTVAR
         LN = L + L - 1
         LM2 = LM + LN + 4
         LM3 = LM - LN
         DO M = 1,LN
            LO = LM2 + M
            LP = LM3 + M
            LQ = LM + M + 1
            YLM(LO) = -(FAC2(LP)*YLM(LP))/FAC2(LQ)
         END DO
         LM = LM + L + L + 1
      END DO
      IF ( IP.GT.4 .AND. ICE.EQ.1 ) WRITE (NOUT1,99002) LMAX
      IF ( IP.GT.4 .AND. ICE.EQ.2 ) WRITE (NOUT1,99003) LMAX
      IF ( IP.GT.4 ) WRITE (NOUT1,99004)
      K = 1
      II = 0
 100  CONTINUE
      LLTVAR = LMAX + II
      DO IL2 = 1,LLTVAR
         L2 = IL2 - II
         M2 = (-L2) + 1 - II
         DO I2 = 1,IL2
            DO IL3 = 1,LLTVAR
               L3 = IL3 - II
               M3 = (-L3) + 1 - II
               DO I3 = 1,IL3
                  LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                  LB1 = L2 + L3
                  LA11 = LA1 + 1
                  LB11 = LB1 + 1
                  M1 = M3 - M2
                  DO L11 = LA11,LB11,2
                     L1 = LA11 + LB11 - L11 - 1
                     L = (L2-L1-L3)/2 + M1 + M3
                     M = L1*(L1+1) - M1 + 1
                     ALM = (-1.D0)**L*4.D0*PI*BLM(L1,M1,L2,M2,L3,(-M3))
                     CELM(K) = ALM
                     IF ( ICE.EQ.1 ) CELM(K) = YLM(M)*ALM
                     IF ( IP.GT.4 ) WRITE (NOUT1,99005) K,L1,M1,L2,M2,
     &                    L3,M3,CELM(K)
                     K = K + 1
                  END DO
                  M3 = M3 + 2
               END DO
            END DO
            M2 = M2 + 2
         END DO
      END DO
      IF ( II.GT.0 ) THEN
         II = 0
         II1 = 1
 150     CONTINUE
         LLTVAR = LMAX + II
         LL1 = LMAX + II1
         DO IL2 = 1,LLTVAR
            L2 = IL2 - II
            M2 = (-L2) + 1 - II
            DO I2 = 1,IL2
               DO IL3 = 1,LL1
                  L3 = IL3 - II1
                  M3 = (-L3) + 1 - II1
                  DO I3 = 1,IL3
                     LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                     LB1 = L2 + L3
                     LA11 = LA1 + 1
                     LB11 = LB1 + 1
                     M1 = M3 - M2
                     DO L11 = LA11,LB11,2
                        L1 = LA11 + LB11 - L11 - 1
                        L = (L2-L1-L3)/2 + M1 + M3
                        M = L1*(L1+1) - M1 + 1
                        ALM = (-1.D0)
     &                        **L*4.D0*PI*BLM(L1,M1,L2,M2,L3,(-M3))
                        CELM(K) = ALM
                        IF ( ICE.EQ.1 ) CELM(K) = YLM(M)*ALM
                        IF ( IP.GT.4 ) WRITE (NOUT1,99005) K,L1,M1,L2,
     &                       M2,L3,M3,CELM(K)
                        K = K + 1
                     END DO
                     M3 = M3 + 2
                  END DO
               END DO
               M2 = M2 + 2
            END DO
         END DO
         IF ( II.GT.0 ) THEN
            RETURN
         ELSE
            II = 1
            II1 = 0
            GOTO 150
         END IF
      ELSE
         II = 1
         GOTO 100
      END IF
99001 FORMAT ('0*** lmax =',i3,' too big to be handled by blm',
     &        ' routine')
99002 FORMAT ('1clm''s for xmat',5x,'lmax =',i3)
99003 FORMAT ('1elm''s for xmatk',5x,'lmax =',i3)
99004 FORMAT ('0',2x,'k',6x,'l1 m1 l2 m2 l3 m3',/)
99005 FORMAT (' ',i5,5x,6I4,5x,e15.5)
      END
C*==check1.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CHECK1(NSPIN,NREL,GANZ,MAXG)
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,MAXG,NREL,NSPIN
C
C*** End of declarations rewritten by SPAG
C
      IF ( NREL.GT.1 ) THEN
         WRITE (NOUT1,99001) NREL
         STOP
      END IF
      IF ( NSPIN.NE.1 .AND. NSPIN.NE.2 ) THEN
         WRITE (NOUT1,99002) NSPIN
         STOP
      END IF
      IF ( GANZ.GT.MAXG ) THEN
         WRITE (NOUT1,99003) GANZ,MAXG
         STOP
      END IF
      RETURN
99001 FORMAT (1x,'nrel= ',1I3,'not allowed ')
99002 FORMAT (1x,'nspin= ',1I3,'not allowed ')
99003 FORMAT (1x,'ganz= ',1I4,' is to large.')
      END
C*==clegor.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      REAL*8 FUNCTION CLEGOR(KAPPA,MU,S)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KAPPA,MU,S
C
C*** End of declarations rewritten by SPAG
C
      CLEGOR = 0.D0
C
      IF ( KAPPA.GT.0.0D0 ) THEN
C         IF ( S.EQ.0.5D0 ) THEN
         IF ( ABS(S-0.5D0).LE.1.0D-16 ) THEN
            CLEGOR = -SQRT((KAPPA-MU+0.5D0)/(2.D0*KAPPA+1.D0))
C         ELSE IF ( S.EQ.(-0.5D0) ) THEN
         ELSE IF ( ABS(S+0.5D0).LE.1.0D-16 ) THEN
            CLEGOR = SQRT((KAPPA+MU+0.5D0)/(2.D0*KAPPA+1.D0))
         END IF
      END IF
      IF ( KAPPA.LT.0.0D0 ) THEN
C         IF ( S.EQ.0.5D0 ) THEN
         IF ( ABS(S-0.5D0).LE.1.0D-16 ) THEN
            CLEGOR = SQRT((KAPPA-MU+0.5D0)/(2.D0*KAPPA+1.D0))
C         ELSE IF ( S.EQ.(-0.5D0) ) THEN
         ELSE IF ( ABS(S+0.5D0).LE.1.0D-16 ) THEN
            CLEGOR = SQRT((KAPPA+MU+0.5D0)/(2.D0*KAPPA+1.D0))
         END IF
      END IF
C
      END
C*==clmgo.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CLMGO(CL)
      USE MOD_SPEC,ONLY:ML,MLSP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CL(MLSP)
C
C Local variables
C
      REAL*8 FL,FLB,FMU,T1,T2
      INTEGER IN,IS,ISI,IT,JM,K,KN,KS,L,LK,MAXL1,MN,MU,NSI
C
C*** End of declarations rewritten by SPAG
C
C      integer mlr, mlt, mla, mlb
C      parameter(mlr=(ml+1)*ml/2,mlt=ml*(ml-1)/2)
C      parameter (mla = 2*ml - 1, mlb = mla*mla)
C
      KN = 0
      MAXL1 = ML
      DO LK = 1,MAXL1
         DO NSI = 1,2
            ISI = (-1)**(NSI+1)
            K = ISI*(LK-1) + (ISI-1)/2
            IF ( K.NE.0 ) THEN
               KS = -K/IABS(K)
               L = K
               IF ( K.LT.0 ) L = ((-K)) - 1
               FL = DBLE(L)
               FLB = DBLE(2*L+1)
               JM = 2*IABS(K)
               MU = JM + 1
               DO MN = 1,JM
                  MU = MU - 2
                  FMU = DBLE(MU)*0.5D0
                  T1 = FL + FMU + 0.5D0
                  T2 = FL - FMU + 0.5D0
                  DO IN = 1,2
                     IS = (-1)**(IN+1)
                     KN = KN + 1
                     IT = K*IS
                     IF ( IT.LT.0 ) THEN
                        CL(KN) = SQRT(T1/FLB)
                     ELSE
                        CL(KN) = DBLE(KS)*SQRT(T2/FLB)
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
C
      END
C*==cmatpr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CMATPR(T,M,IDIMVAR,NOUT1)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IDIMVAR,NOUT1
      CHARACTER*8 T
      COMPLEX*16 M(IDIMVAR,IDIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      WRITE (NOUT1,99001) T
      DO I = 1,IDIMVAR
         WRITE (NOUT1,99002) (M(I,J),J=1,IDIMVAR)
      END DO
      RETURN
99001 FORMAT (1x,'output matrix ',1A8)
99002 FORMAT (5(e12.5,e12.5))
      END
C*==copy.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE COPY(X,Y,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),Y(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Y(I,J) = X(I,J)
         END DO
      END DO
      END
C*==xcopr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XCOPR(X,XEE,XOO,XEO,XOE,DIM1,DIM2,DIM3,IR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIM1,DIM2,DIM3,IR
      COMPLEX*16 X(DIM1,DIM1),XEE(DIM2,DIM2),XEO(DIM2,DIM3),
     &           XOE(DIM3,DIM2),XOO(DIM3,DIM3)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      IF ( IR.EQ.1 ) THEN
         DO I = 1,DIM3
            DO J = 1,DIM3
               X(I,J) = XOO(I,J)
            END DO
         END DO
         DO I = 1,DIM2
            DO J = 1,DIM2
               X(I+DIM3,J+DIM3) = XEE(I,J)
            END DO
         END DO
         DO I = 1,DIM2
            DO J = 1,DIM3
               X(I+DIM3,J) = XEO(I,J)
            END DO
         END DO
         DO I = 1,DIM3
            DO J = 1,DIM2
               X(I,J+DIM3) = XOE(I,J)
            END DO
         END DO
      ELSE IF ( IR.EQ.2 ) THEN
         DO I = 1,DIM3
            DO J = 1,DIM3
               XOO(I,J) = X(I,J)
            END DO
         END DO
         DO I = 1,DIM2
            DO J = 1,DIM2
               XEE(I,J) = X(I+DIM3,J+DIM3)
            END DO
         END DO
         DO I = 1,DIM2
            DO J = 1,DIM3
               XEO(I,J) = X(I+DIM3,J)
            END DO
         END DO
         DO I = 1,DIM3
            DO J = 1,DIM2
               XOE(I,J) = X(I,J+DIM3)
            END DO
         END DO
      END IF
      END
C*==fac.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
C     ***************************************************************
C
      REAL*8 FUNCTION FAC(N,M)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
C
C Local variables
C
      INTEGER I,M1
C
C*** End of declarations rewritten by SPAG
C
      FAC = 1.D0
      IF ( M.EQ.0 ) RETURN
      M1 = N - M + 1
      DO I = M1,N
         FAC = FAC*DBLE(I)
      END DO
C
      END
C*==zcerf.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
C     ***************************************************************
C
      COMPLEX*16 FUNCTION ZCERF(Z,EMACH)
      USE MOD_SPEC,ONLY:PI,CZERO,CONE,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EMACH
      COMPLEX*16 Z
C
C Local variables
C
      REAL*8 ABSZ,ABTERM,API,AX,AY,EPS,FACT,FACTD,FACTN,Q,RTPI,TEST,X,Y
      COMPLEX*16 CI,H1,H2,H3,SUMVAR,TERM1,TERM2,U1,U2,U3,XZZS,ZZ,ZZS
      INTEGER N,NN
C
C*** End of declarations rewritten by SPAG
C
C     cerf, given complex argument z, provides the complex
C     error function:     w(z)=exp(-z**2)*(1.0-erf(-i*z))
C     the evaluation always takes place in the first quadrant.
C     one of three methods is employed depending on the size
C     of the argument: a power series, a recurrence based on
C     continued fractions theory, or an asymptotic series.
C
C
      EPS = 5.D0*EMACH
      API = 1.D0/PI
      CI = CIMAG
      ABSZ = CDABS(Z)
C      IF ( ABSZ.EQ.0. ) THEN
      IF ( ABS(ABSZ).LT.1.0D-16 ) THEN
         ZCERF = CONE
         RETURN
C
C       the argument is translated to the first quadrant from
C       the nn'th quadrant, before the method for the function
C       evaluation is chosen
C
      END IF
      X = DBLE(Z)
      Y = DIMAG(Z)
      AX = ABS(X)
      AY = ABS(Y)
      ZZ = DCMPLX(AX,AY)
      ZZS = ZZ*ZZ
      NN = 1
C      IF ( X.NE.AX ) NN = 2
      IF ( ABS(X-AX).GT.1.0D-16 ) NN = 2
C      IF ( Y.NE.AY ) NN = 5 - NN
      IF ( ABS(Y-AY).GT.1.0D-16 ) NN = 5 - NN
      IF ( ABSZ.LE.10.D0 ) THEN
         IF ( AY.GE.1.D0 .OR. ABSZ.GE.4.D0 ) THEN
C
C       continued fractions theory: w(z) is related to the limiting
C       value of u<n,z>/h<n,z>, where u and h obey the same
C       recurrence relation in n. see faddeeva and terent'ev:
C       tables of values of w(z) for complex arguments,pergamon,
C       n.y. 1961
C
            TERM2 = DCMPLX(1.D6,0.D0)
            Q = 1.D0
            H1 = CONE
            H2 = 2.D0*ZZ
            U1 = CZERO
            RTPI = 2.D0*SQRT(PI)
            U2 = DCMPLX(RTPI,0.D0)
 20         CONTINUE
            TERM1 = TERM2
            DO N = 1,5
               H3 = H2*ZZ - Q*H1
               U3 = U2*ZZ - Q*U1
               H1 = H2
               H2 = 2.D0*H3
               U1 = U2
               U2 = 2.D0*U3
               Q = Q + 1.D0
            END DO
            TERM2 = U3/H3
            TEST = CDABS((TERM2-TERM1)/TERM1)
            IF ( TEST-EPS.LT.0.D0 ) THEN
               ZCERF = API*CI*TERM2
               GOTO 100
            ELSE IF ( Q-60.D0.LE.0.D0 ) THEN
               GOTO 20
            END IF
         END IF
C
C       power series: see abramowitz and stegun's handbook of
C       mathematical functions, p297
C
         Q = 1.D0
         FACTN = -1.D0
         FACTD = 1.D0
         TERM1 = ZZ
         SUMVAR = ZZ
 50      CONTINUE
         DO N = 1,5
            FACTN = FACTN + 2.D0
            FACTD = FACTD + 2.D0
            FACT = FACTN/(Q*FACTD)
            TERM1 = FACT*ZZS*TERM1
            SUMVAR = SUMVAR + TERM1
            Q = Q + 1.D0
         END DO
         ABTERM = CDABS(TERM1)
         IF ( ABTERM-EPS.GE.0.D0 ) THEN
            IF ( Q-100.D0.LT.0.D0 ) GOTO 50
         END IF
         FACT = 2.D0*SQRT(API)
         SUMVAR = FACT*CI*SUMVAR
         XZZS = CDEXP((-ZZS))
         ZCERF = XZZS + XZZS*SUMVAR
         GOTO 100
C
C       asymptotic series: see abramowitz and stegun, p328
C
      END IF
      ZCERF = 0.5124242D0/(ZZS-0.2752551D0)
     &        + 0.05176536D0/(ZZS-2.724745D0)
      ZCERF = CI*ZZ*ZCERF
C
C       symmetry relations are now used to transform the function
C       back to quadrant nn
C
 100  CONTINUE
      IF ( NN.EQ.1 ) GOTO 99999
      IF ( NN.EQ.2 ) THEN
      ELSE IF ( NN.EQ.3 ) THEN
         ZCERF = 2.D0*CDEXP((-ZZS)) - ZCERF
         GOTO 99999
      ELSE
         ZCERF = 2.D0*CDEXP((-ZZS)) - ZCERF
      END IF
      ZCERF = DCONJG(ZCERF)
      RETURN
99999 CONTINUE
      END
C*==gikam.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE GIKAM(GODD,GEVEN,APQ,BPQ,E,VPI,AK1,AK2,CELM,EMACH,
     &                 CLIGHT,NREL,NEV,NOD)
C     /****************************************************************/
C     # subroutines and functions called from this routine             *
C       zcerf                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:ML,MLQ,MLZP,PI,EPS12,CZERO,CONE,CIMAG
      USE MOD_SPEC_GEOM,ONLY:AR1,AR2,TV
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLA,MLD,MLC,MLE
      PARAMETER (MLA=2*ML-1,MLD=ML*(2*ML-1),MLC=4*ML-3,MLE=ML*(ML+1)
     &           *(2*ML+1)/6)
C
C Dummy arguments
C
      REAL*8 AK1,AK2,APQ,BPQ,CLIGHT,E,EMACH,VPI
      INTEGER NEV,NOD,NREL
      REAL*8 CELM(MLZP)
      COMPLEX*16 GEVEN(NEV,NEV),GODD(NOD,NOD)
C
C Local variables
C
      REAL*8 AB1,AB2,AC,ACSQ,AD,AK(2),AKPT(2),AL,AN,AN1,AN2,AP,AP1,AP2,
     &       AR,B,B1(2),B2(2),DENOM(:),DIZZ,DNORM,FACT(:),R(2),RTPI,RTV,
     &       TEST,TEST1,TEST2
      COMPLEX*16 CA,CAA,CAB,CAC,CAF(:),CAGK(:),CALM,CALPHA,CANT,CAPSQ,
     &           CEN,CF,CFF,CGAM,CGK,CGKK,CGKN(:),CKAPPA,CNSQ,CP,CPA,
     &           CPK,CPSQ,CREF(:),CSD,CT,CTA,CTAI,CU,CU1,CU2,CW,CWW,CX,
     &           CXPM(:),CZ,DLM(:),Z,ZZ
      INTEGER I,I1,I2,I3,II,IL,IL2,IL3,IN,IT,J,J1,J2,K,KK,L,L1,L2,L2MAX,
     &        L3,L6,LA1,LA11,LB1,LB11,LL2,LLL,M,M2,M3,MM,N,N1,N1MAX,NA,
     &        NI,NM,NN,NNDLM
      COMPLEX*16 ZCERF
      EXTERNAL ZCERF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CAF,DLM,CAGK,FACT,CREF,CGKN,CXPM,DENOM
      ALLOCATE (CAF(ML),DLM(MLD),CAGK(MLA),FACT(MLC),CREF(MLQ))
      ALLOCATE (CGKN(ML),CXPM(MLA),DENOM(MLE))
C
C
      N1MAX = 20
      RTPI = SQRT(PI)
C
      DO I = 1,ML
         CAF(I) = DCMPLX(0.D0,-1.D0)
      END DO
C
      CEN = DCMPLX(E,VPI+0.5D0*EMACH)
      CAPSQ = 2.D0*CEN
      IF ( NREL.EQ.1 ) CAPSQ = CEN*DCMPLX(2.D0,0.D0)
     &                         + CEN/(CLIGHT*CLIGHT)
      CKAPPA = CDSQRT(CAPSQ)
      AK(1) = AK1 + APQ
      AK(2) = AK2 + BPQ
      L2MAX = (ML-1) + (ML-1)
      LL2 = L2MAX + 1
C
C     prepare factorials
      FACT(1) = 1.D0
      NI = L2MAX + L2MAX
      DO I = 1,NI
         FACT(I+1) = DBLE(I)*FACT(I)
      END DO
C
C     reset dlm
      NNDLM = L2MAX*(L2MAX+3)/2 + 1
      DO I = 1,NNDLM
         DLM(I) = CZERO
      END DO
C
      CALPHA = CAPSQ*DCMPLX(TV/(4.D0*PI),0.D0)
      AL = CDABS(CALPHA)
      IF ( EXP(AL)*EMACH.GT.5.D-5 ) THEN
         AL = DLOG(5.0D-5/EMACH)
         CALPHA = DCMPLX(SIGN(1.D0,E)*AL,0.D0)
      END IF
C
      CTA = CDSQRT(CALPHA)
      K = 1
      KK = 1
      AP1 = -2.D0/TV
      AP2 = -1.D0
      CF = CIMAG/CKAPPA
      DO L = 1,LL2
         AP1 = AP1/2.D0
         AP2 = AP2 + 2.D0
         CFF = CF
         MM = 1
         IF ( MOD(L,2).EQ.0 ) THEN
            MM = 2
            CFF = CIMAG*CFF
         END IF
         NN = (L-MM)/2 + 2
         DO M = MM,L,2
            J1 = L + M - 1
            J2 = L - M + 1
            AP = AP1*SQRT(AP2*FACT(J1)*FACT(J2))
            CREF(KK) = AP*CFF
            CFF = -CFF
            KK = KK + 1
            NN = NN - 1
            DO I = 1,NN
               I1 = I
               I2 = NN - I + 1
               I3 = NN + M - I
               DENOM(K) = 1.D0/(FACT(I1)*FACT(I2)*FACT(I3))
               K = K + 1
            END DO
         END DO
      END DO
      RTV = 2.D0*PI/TV
      B1(1) = -AR1(2)*RTV
      B1(2) = AR1(1)*RTV
      B2(1) = -AR2(2)*RTV
      B2(2) = AR2(1)*RTV
      TEST1 = 1.D6
      II = 1
C
      DO N1 = 0,N1MAX
         NA = N1 + N1 + II
         AN1 = DBLE(N1)
         AN2 = (-AN1) - 1.D0
         DO I1 = 1,NA
            AN2 = AN2 + 1.D0
            DO I2 = 1,4
               AN = AN1
               AN1 = -AN2
               AN2 = AN
               AB1 = AN1*B1(1) + AN2*B2(1)
               AB2 = AN1*B1(2) + AN2*B2(2)
               AKPT(1) = AK(1) + AB1
               AKPT(2) = AK(2) + AB2
               ACSQ = AKPT(1)*AKPT(1) + AKPT(2)*AKPT(2)
               CPSQ = CAPSQ - DCMPLX(ACSQ,0.D0)
               AC = SQRT(ACSQ)
               CP = CDSQRT(CPSQ)
               CPK = CZERO
               CGK = CZERO
               CGKK = CONE
               IF ( AC.GT.EMACH ) THEN
                  CPK = DCMPLX(AKPT(1)/AC,AKPT(2)/AC)
                  CGK = AC/CKAPPA
                  CGKK = CPSQ/CAPSQ
               END IF
               CXPM(1) = CONE
               CAGK(1) = CONE
               DO I = 2,LL2
                  CXPM(I) = CXPM(I-1)*CPK
                  CAGK(I) = CAGK(I-1)*CGK
               END DO
               CF = CKAPPA/CP
               ZZ = -CALPHA*CGKK
               CZ = CDSQRT((-ZZ))
               Z = -CIMAG*CZ
               DIZZ = DIMAG(CZ)
               IF ( DIZZ.LT.0.D0 .AND. DBLE(ZZ).GT.150.D0 ) THEN
                  CGAM = DCMPLX(2.D0*RTPI,0.0D0)
                  IF ( DIZZ.GT.0.0D0 )
     &                 CGAM = CGAM*CDEXP(DCMPLX(0.D0,(-2.D0*DIZZ)))
                  CX = CZERO
               ELSE
                  CX = CDEXP((-ZZ))
                  CGAM = RTPI*CX*ZCERF(CZ,EMACH)
               END IF
               CGKN(1) = CF*CGAM
               CT = Z
               B = 0.5D0
               LLL = L2MAX/2 + 1
               DO I = 2,LLL
                  CT = CT/ZZ
                  B = B - 1.D0
                  CGAM = (CGAM-CX*CT)/B
                  CF = CF*CGKK
                  CGKN(I) = CF*CGAM
               END DO
               K = 1
               KK = 1
               DO L = 1,LL2
                  MM = 1
                  IF ( MOD(L,2).EQ.0 ) MM = 2
                  N = (L*L+MM)/2
                  NN = (L-MM)/2 + 2
                  DO M = MM,L,2
                     CAC = CZERO
                     NN = NN - 1
                     IL = L
                     DO I = 1,NN
                        CAC = CAC + DENOM(K)*CAGK(IL)*CGKN(I)
                        IL = IL - 2
                        K = K + 1
                     END DO
                     CAC = CREF(KK)*CAC
                     IF ( AC.GT.1.D-6 ) THEN
                        DLM(N) = DLM(N) + CAC/CXPM(M)
                        IF ( M.NE.1 ) THEN
                           NM = N - M + 1
                           DLM(NM) = DLM(NM) + CAC*CXPM(M)
                        END IF
                     ELSE
                        NM = N - M + 1
                        DLM(NM) = DLM(NM) + CAC*CXPM(M)
                     END IF
                     KK = KK + 1
                     N = N + 1
                  END DO
               END DO
               IF ( II.GT.0 ) EXIT
            END DO
            II = 0
         END DO
C
         TEST2 = 0.D0
         DO I = 1,NNDLM
            DNORM = CDABS(DLM(I))
            TEST2 = TEST2 + DNORM*DNORM
         END DO
         TEST = ABS((TEST2-TEST1)/TEST1)
         TEST1 = TEST2
         IF ( TEST.LE.EMACH ) EXIT
         IF ( N1.EQ.N1MAX ) WRITE (6,99001) N1
      END DO
C
      KK = 1
      AP1 = TV/(4.D0*PI)
      CF = CAPSQ/CIMAG
      DO L = 1,LL2
         CFF = CF
         MM = 1
         IF ( MOD(L,2).EQ.0 ) THEN
            MM = 2
            CFF = -CIMAG*CFF
         END IF
         J1 = (L-MM)/2 + 1
         J2 = J1 + MM - 1
         IN = J1 + L - 2
         AP2 = AP1*((-1.D0)**IN)
         DO M = MM,L,2
            AP = AP2/(FACT(J1)*FACT(J2))
            CREF(KK) = AP*CFF*CREF(KK)
            J1 = J1 - 1
            J2 = J2 + 1
            AP2 = -AP2
            CFF = -CFF
            KK = KK + 1
         END DO
      END DO
C
      DO N1 = 1,N1MAX
         NA = N1 + N1
         AN1 = DBLE(N1)
         AN2 = (-AN1) - 1.D0
         DO I1 = 1,NA
            AN2 = AN2 + 1.D0
            DO I2 = 1,4
               AN = AN1
               AN1 = -AN2
               AN2 = AN
               R(1) = AN1*AR1(1) + AN2*AR2(1)
               R(2) = AN1*AR1(2) + AN2*AR2(2)
               AR = SQRT(R(1)*R(1)+R(2)*R(2))
               CPK = DCMPLX(R(1)/AR,R(2)/AR)
               CXPM(1) = CONE
               DO I = 2,LL2
                  CXPM(I) = CXPM(I-1)*CPK
               END DO
               AD = AK(1)*R(1) + AK(2)*R(2)
               CSD = CDEXP((-AD*CIMAG))
               CANT = 0.5D0*AR*CKAPPA
               CNSQ = CANT*CANT
               Z = CIMAG*CANT/CTA
               ZZ = CTA - Z
               Z = CTA + Z
               CWW = ZCERF((-ZZ),EMACH)
               CW = ZCERF(Z,EMACH)
               CAA = 0.5D0*RTPI*(CW-CWW)/CIMAG
               CAB = 0.5D0*RTPI*(CW+CWW)
               CA = CALPHA - CNSQ/CALPHA
               CPA = CDEXP(CA)
               CU1 = CAA*CPA
               CU2 = CAB*CPA/CANT
               KK = 1
               AL = -0.5D0
               CFF = CTA
               CF = CONE
               DO L = 1,LL2
                  MM = 1
                  IF ( MOD(L,2).EQ.0 ) MM = 2
                  N = (L*L+MM)/2
                  DO M = MM,L,2
                     CAC = CREF(KK)*CU2*CF*CSD
                     DLM(N) = DLM(N) + CAC/CXPM(M)
                     IF ( M.NE.1 ) THEN
                        NM = N - M + 1
                        DLM(NM) = DLM(NM) + CAC*CXPM(M)
                     END IF
                     KK = KK + 1
                     N = N + 1
                  END DO
                  AL = AL + 1.D0
                  CFF = CFF/CALPHA
                  CU = (AL*CU2-CU1+CFF*CPA)/CNSQ
                  CU1 = CU2
                  CU2 = CU
                  CF = CANT*CF
               END DO
            END DO
         END DO
C
C         after each step of the summation a test on the
C         convergence of the elements of dlm is made
C
         TEST2 = 0.D0
         DO I = 1,NNDLM
            DNORM = CDABS(DLM(I))
            TEST2 = TEST2 + DNORM*DNORM
         END DO
C
         TEST = ABS((TEST2-TEST1)/TEST1)
         TEST1 = TEST2
         IF ( TEST.LE.EMACH ) EXIT
         IF ( N1.EQ.N1MAX ) WRITE (6,99002) N1
      END DO
C
      CPA = CDEXP((-CALPHA))
      CTAI = CONE/(RTPI*CTA)
      CAC = CKAPPA*(CIMAG*(CPA-ZCERF(CTA,EMACH))-CTAI)/CPA
      AP = -0.5D0/RTPI
      DLM(1) = DLM(1) + AP*CAC
C
C       finally the elements of dlm are multiplied by the
C       factor (-1.0)**((m+!m!)/2)
C
      DO L = 2,LL2,2
         N = L*L/2 + 1
         DO M = 2,L,2
            DLM(N) = -DLM(N)
            N = N + 1
         END DO
      END DO
C
      K = 1
      DO II = 0,1
         L1 = ML - 1 + II
         I = 1
         DO IL2 = 1,L1
            L2 = IL2 - II
            M2 = ((-L2)) + 1 - II
            DO I2 = 1,IL2
               J = 1
               DO IL3 = 1,L1
                  L3 = IL3 - II
                  M3 = ((-L3)) + 1 - II
                  DO I3 = 1,IL3
                     CALM = CZERO
                     LA1 = MAX0(IABS(L2-L3),IABS(M2-M3))
                     LB1 = L2 + L3
                     N = (LB1*(LB1+2)+M2-M3+2)/2
                     NN = 2*LB1
                     LB11 = LB1 + 1
                     LA11 = LA1 + 1
                     DO L6 = LA11,LB11,2
                        CALM = CALM + CELM(K)*DLM(N)
                        N = N - NN
                        NN = NN - 4
                        K = K + 1
                     END DO
                     CALM = CALM/CKAPPA
                     IF ( I.EQ.J ) CALM = CALM + CIMAG
                     IF ( II.LE.0 ) THEN
                        GODD(J,I) = -CAF(L3+1)*CALM
                     ELSE
                        GEVEN(J,I) = -CAF(L3+1)*CALM
                     END IF
                     M3 = M3 + 2
                     J = J + 1
                  END DO
               END DO
               M2 = M2 + 2
               I = I + 1
            END DO
         END DO
      END DO
C
      IT = 1
C
      DO J = 1,MLD
         IF ( CDABS(DLM(J)).LT.EPS12 ) DLM(J) = CZERO
      END DO
C
      IF ( IP.GT.4 ) WRITE (NOUT1,99003) IT,CEN
      IF ( IP.GT.4 ) WRITE (NOUT1,99004) (DLM(J),J=1,MLD)
      RETURN
C
99001 FORMAT (2x,'** dlm1s not converged in gikam by n1 =',i3)
99002 FORMAT (2x,'** dlm2s not converged in gikam by n1 =',i3)
99003 FORMAT (1x,'dmatk1 : layer type it=',i2,' at energy e= ',2F12.4)
99004 FORMAT (6E12.4)
      END
C*==gikamm.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE GIKAMM(IT,E,VPI,CLIGHT,NATL,POS,NNSK,NNSR,NSK,NSR,RLVS,
     &                  RLVX,RLVY,VS,VX,VY,AKX,AKY,NREL,DLM1)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C zcerf                                                        *
C**************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,NPM,PI,CZERO,CONE,CIMAG
      USE MOD_SPEC_GEOM,ONLY:TV
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NFM,LMXM3,LMAXM,L2MAXM,L4MAXM,NDLMM
      PARAMETER (NFM=NATLM*NATLM-NATLM+1,LMXM3=ML*ML*ML,LMAXM=ML-1,
     &           L2MAXM=2*LMAXM,L4MAXM=L2MAXM*2,NDLMM=(L2MAXM+1)
     &           *(L2MAXM+1))
C
C Dummy arguments
C
      REAL*8 AKX,AKY,CLIGHT,E,VPI
      INTEGER IT,NNSK,NNSR,NREL
      COMPLEX*16 DLM1(NDLMM,NFM)
      INTEGER NATL(LAYSM),NSK(NPM),NSR(NPM)
      REAL*8 POS(3,NATLM,LAYSM),RLVS(NPM),RLVX(NPM),RLVY(NPM),VS(NPM),
     &       VX(NPM),VY(NPM)
C
C Local variables
C
      COMPLEX*16 ACC,ALPHA,C0,C1,CEN,CI,DLM2(:,:),DLM3(:),F1,F2,F3,F4,
     &           F5,F55,F6,F66,F7,GAM,GAMN(:),GK,GKN(:),GKZ,GKZ2,I0,IM,
     &           IX,KA,KASQ,PHI,PHIM(:),PREF1(:),PREF2(:),RTAL
      REAL*8 AL,B,DEN(:),DISTX,DISTY,DNORM,EMACH,FAC(:),FOUR,G1,G2,G3,
     &       G4,GKP,GKP2,GX,GY,H,HALF,ONE,PI4,PREF,R,RSQ,RTPI,RX,RXX,RY,
     &       RYY,TEST,TESTA,TESTB,TWO,ZERO
      INTEGER I,I1,I2,I3,IA,IA1,IB,IDEN,IFL,IL,ILM,IPREF,IRINGA,IRINGB,
     &        J,JA,JA1,KA1,L,L1,L2MAX,L4MAX,LC1,LC2,LM,LMAX,LMBOT,LMM,
     &        LMM2,LMTOP,LPM,LPM2,M,MPHI,N,N1,NAT,NDLM,NF,NMAX
      COMPLEX*16 ZCERF
      EXTERNAL ZCERF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FAC,DEN,GKN,DLM2,DLM3,GAMN,PHIM,PREF1,PREF2
      ALLOCATE (FAC(0:L4MAXM),DEN(LMXM3),GKN(0:L2MAXM))
      ALLOCATE (DLM2(NDLMM,NFM),DLM3(NDLMM),GAMN(0:LMAXM))
      ALLOCATE (PHIM(0:L4MAXM),PREF1(NDLMM),PREF2(NDLMM))
C
      EMACH = 1.D-9
      ZERO = 0.D0
      HALF = 0.5D0
      ONE = 1.D0
      TWO = 2.D0
      FOUR = 4.D0
      C0 = CZERO
      C1 = CONE
      CI = CIMAG
      RTPI = SQRT(PI)
      PI4 = FOUR*PI
      CEN = DCMPLX(E,VPI+0.5D0*EMACH)
      KASQ = 2.D0*CEN
      IF ( NREL.EQ.1 ) KASQ = CEN*(2.D0+CEN/(CLIGHT*CLIGHT))
      KA = CDSQRT(KASQ)
      LMAX = ML - 1
      L2MAX = 2*(ML-1)
      L4MAX = L2MAX + L2MAX
      NDLM = (L2MAX+1)*(L2MAX+2)/2
      NAT = NATL(IT)
      NF = NAT*NAT - NAT + 1
      DO I = 1,NF
         DO J = 1,NDLM
            DLM1(J,I) = C0
            DLM2(J,I) = C0
         END DO
      END DO
      FAC(0) = ONE
      DO I = 1,L4MAX
         FAC(I) = FAC(I-1)*DBLE(I)
      END DO
      ALPHA = TV*KASQ/PI4
      AL = CDABS(ALPHA)
      IF ( EXP(AL)*EMACH-5.D-5.GT.0. ) THEN
         AL = DLOG(5.D-5/EMACH)
         ALPHA = DCMPLX(SIGN(1.D0,E)*AL,0.D0)
      END IF
      RTAL = CDSQRT(ALPHA)
      F1 = -CI/(TV*KASQ)
      F2 = -CI
      G1 = TWO
      G2 = -ONE
      G3 = -ONE
      IPREF = 0
      IDEN = 0
      DO L = 0,L2MAX
         G1 = G1*HALF
         G2 = G2 + TWO
         G3 = -G3
         G4 = G3
         F2 = F2*CI
         F3 = F2
         DO M = -L,L,2
            LPM = L + M
            LMM = L - M
            IPREF = IPREF + 1
            PREF = G1*SQRT(G2*FAC(LPM)*FAC(LMM))
            PREF1(IPREF) = PREF*F1*F3
            LPM2 = LPM/2
            LMM2 = LMM/2
            PREF = PREF/(FAC(LMM2)*FAC(LPM2))
            PREF2(IPREF) = -PREF*G4/PI4
            F3 = -F3
            G4 = -G4
            IF ( M.LE.0 ) THEN
               DO N1 = 0,LPM2
                  I2 = LMM2 - N1
                  I3 = LPM2 - N1
                  IDEN = IDEN + 1
                  DEN(IDEN) = ONE/(FAC(N1)*FAC(I2)*FAC(I3))
               END DO
            END IF
         END DO
      END DO
      TESTB = ZERO
      IRINGB = 1
      IB = 0
      DO I1 = 1,NNSK
         DO I2 = 1,NSK(I1)
            IB = IB + 1
C
C     prepare local arrays for dlm1 contributions
C     gkn holds (k''+g)**l for l=0 --> l2max
C     phim holds exp(-i*m(phi(k''+g))) for m=-l2max --> 0
C
            GX = AKX + RLVX(IB)
            GY = AKY + RLVY(IB)
            GKP2 = GX*GX + GY*GY
            GKZ2 = KASQ - GKP2
            GKP = SQRT(GKP2)
            GKZ = CDSQRT(GKZ2)
            IF ( GKP.LT.100.E0*EMACH ) THEN
               PHI = C1
               GK = C0
            ELSE
               PHI = DCMPLX(GX/GKP,GY/GKP)
               GK = GKP/KA
            END IF
            GKN(0) = C1
            PHIM(L2MAX) = C1
            DO I = 1,L2MAX
               GKN(I) = GKN(I-1)*GK
               PHIM(L2MAX-I) = PHIM(L2MAX-I+1)*PHI
            END DO
C
C     gamm holds (kgz/kappa)**(2n-1) for n=0 --> lmax
C
            F1 = GKZ2/KASQ
            F2 = KA/GKZ
            F3 = -ALPHA*F1
            F4 = CDSQRT((-F3))
            F5 = -CI*F4
            F6 = CDEXP((-F3))
            GAM = RTPI*F6*ZCERF(F4,EMACH)
            GAMN(0) = F2*GAM
            F7 = F5
            B = HALF
            DO L = 1,LMAX
               F7 = F7/F3
               B = B - ONE
               GAM = (GAM-F6*F7)/B
               F2 = F2*F1
               GAMN(L) = F2*GAM
            END DO
C
C     update contributions to dlm1's for this k''+g
C
            IDEN = 0
            LMTOP = 1
            DO L = 0,L2MAX
               L1 = L + 1
               LMBOT = LMTOP - 1
               LMTOP = L1*(L1+1)/2 + 1
               LMM = LMBOT
               LPM = LMTOP
               DO M = -L,0,2
                  NMAX = (L+M)/2
                  IL = L + 2
                  ACC = C0
                  DO N = 0,NMAX
                     IL = IL - 2
                     IDEN = IDEN + 1
                     ACC = ACC + DEN(IDEN)*GKN(IL)*GAMN(N)
                  END DO
                  MPHI = L2MAX + M
                  PHI = PHIM(MPHI)
                  LMM = LMM + 1
                  DLM3(LMM) = ACC*PHI
                  DLM1(LMM,1) = DLM1(LMM,1) + DLM3(LMM)
                  IF ( M.NE.0 ) THEN
                     LPM = LPM - 1
                     DLM3(LPM) = ACC/PHI
                     DLM1(LPM,1) = DLM1(LPM,1) + DLM3(LPM)
                  END IF
               END DO
            END DO
            IF ( NAT.GT.1 ) THEN
               LC1 = 0
               DO I = 1,NAT - 1
                  DO J = I + 1,NAT
                     LC1 = LC1 + 2
                     LC2 = LC1 + 1
                     DISTX = POS(2,J,IT) - POS(2,I,IT)
                     DISTY = POS(3,J,IT) - POS(3,I,IT)
                     F1 = CDEXP(CI*(GX*DISTX+GY*DISTY))
                     F2 = DCONJG(F1)
                     LM = 0
                     DO L = 0,L2MAX
                        DO M = -L,L,2
                           LM = LM + 1
                           DLM1(LM,LC1) = DLM3(LM)*F1 + DLM1(LM,LC1)
                           DLM1(LM,LC2) = DLM3(LM)*F2 + DLM1(LM,LC2)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END DO
C         IRINGA = RLVS(IB)/RLVS(2)
         IRINGA = INT(RLVS(IB)/RLVS(2))
         IF ( IRINGA.GE.IRINGB+1 ) THEN
            TESTA = ZERO
            DO J = 1,NF
               DO I = 1,NDLM
                  DNORM = CDABS(DLM1(I,J))
                  TESTA = TESTA + DNORM*DNORM
               END DO
            END DO
            TEST = ABS(TESTB-TESTA)/(TESTA*DBLE(NF*NDLM))
            IF ( TEST.LT.EMACH ) GOTO 100
            IRINGB = IRINGA
            TESTB = TESTA
         END IF
      END DO
      WRITE (NOUT1,99001) I1
      STOP
 100  CONTINUE
      IF ( IP.GT.1 ) WRITE (NOUT1,'(1x,i3,"  recips used in dlm1 sum")')
     &                      IB
      IRINGB = 1
      TESTB = ZERO
      IB = 0
      DO I1 = 1,NNSR
         DO I2 = 1,NSR(I1)
            IB = IB + 1
            RX = VX(IB)
            RY = VY(IB)
            H = AKX*RX + AKY*RY
            F1 = CDEXP((-CI*H))
            LC1 = 0
            DO IA = 1,NAT
               DO JA = IA,NAT
                  IF ( IA.NE.JA .OR. JA.LE.1 ) THEN
                     IA1 = IA
                     JA1 = JA
                     DO IFL = 1,2
                        IF ( IA.NE.JA .OR. IFL.NE.2 ) THEN
                           RXX = RX + POS(2,JA1,IT) - POS(2,IA1,IT)
                           RYY = RY + POS(3,JA1,IT) - POS(3,IA1,IT)
                           RSQ = RXX*RXX + RYY*RYY
                           R = SQRT(RSQ)
                           LC1 = LC1 + 1
                           IF ( R.GE.100.E0*EMACH ) THEN
C
C     set up initial values of recurrence relation
C     integrals im,i0 for calculation of ip
C
                              F3 = KA*R/TWO
                              F2 = F3*F3
                              F4 = CDEXP(ALPHA-F2/ALPHA)
                              F5 = CI*F3/RTAL
                              F55 = RTAL + F5
                              F6 = ZCERF(F55,EMACH)
                              F66 = ((-RTAL)) + F5
                              F7 = ZCERF(F66,EMACH)
                              F5 = (F6+F7)*HALF
                              F6 = (F6-F7)*HALF/CI
                              I0 = RTPI*F4*F5/F3
                              IM = RTPI*F4*F6
                              GKN(0) = I0*F1
                              PHI = DCMPLX(RXX/R,RYY/R)
                              PHIM(L2MAX) = C1
                              G1 = HALF
                              F5 = PHI
                              F6 = F3
                              F7 = F4/RTAL
                              DO L = 1,L2MAX
                                 PHIM(L2MAX-L) = F5
                                 PHIM(L2MAX+L) = DCONJG(F5)
                                 IX = (G1*I0-IM+F7)/F2
                                 GKN(L) = IX*F6*F1
                                 F7 = F7/ALPHA
                                 IM = I0
                                 I0 = IX
                                 F6 = F6*F3
                                 F5 = F5*PHI
                                 G1 = G1 + ONE
                              END DO
                              ILM = 0
                              DO L = 0,L2MAX
                                 F2 = GKN(L)
                                 DO M = -L,L,2
                                    ILM = ILM + 1
                                    MPHI = L2MAX + M
                                    DLM2(ILM,LC1) = DLM2(ILM,LC1)
     &                                 + F2*PHIM(MPHI)
                                 END DO
                              END DO
                              KA1 = IA1
                              IA1 = JA1
                              JA1 = KA1
                           END IF
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END DO
C         IRINGA = VS(IB)/VS(2)
         IRINGA = INT(VS(IB)/VS(2))
         IF ( IRINGA.GE.IRINGB+1 ) THEN
            TESTA = ZERO
            DO J = 1,NF
               DO I = 1,NDLM
                  DNORM = CDABS(DLM2(I,J))
                  TESTA = TESTA + DNORM*DNORM
               END DO
            END DO
            TEST = ABS(TESTB-TESTA)/(TESTA*DBLE(NF*NDLM))
            IF ( TEST.LT.EMACH ) GOTO 200
            IRINGB = IRINGA
            TESTB = TESTA
         END IF
      END DO
      WRITE (NOUT1,99002) I1
      STOP
 200  CONTINUE
      IF ( IP.GT.1 ) WRITE (NOUT1,'(1x,i3," vectors used in dlm2 sum/")'
     &                      ) I
C
C     combine dlm1's ad dlm2's and store result in dlm1
C
      DO I = 1,NF
         DO J = 1,NDLM
            DLM1(J,I) = PREF1(J)*DLM1(J,I) + PREF2(J)*DLM2(J,I)
         END DO
      END DO
C
C     add in correction to l=0,m=0 term in a-a lattice sums
C
      F1 = CDEXP((-ALPHA))
      F2 = 1.D0/(RTPI*RTAL)
      F3 = (CI*(F1-ZCERF(RTAL,EMACH))-F2)/F1
      DLM1(1,1) = DLM1(1,1) - HALF*F3/RTPI
      IF ( IP.GT.3 ) THEN
         WRITE (NOUT1,99003) IT,CEN
         DO I = 1,NF
            WRITE (NOUT1,99004) I
            WRITE (NOUT1,99005) (DLM1(J,I),J=1,NDLM)
         END DO
      END IF
      RETURN
99001 FORMAT (' ** dlm1 not converged in gikamm by i1 =',i3)
99002 FORMAT (' ** dlm2 not converged in gikamm by i1 =',i3)
99003 FORMAT (' dmatkm : layer type it=',i2,' at energy e= ',2F12.4)
99004 FORMAT (' columm : ',i2)
99005 FORMAT (6E12.4)
      END
C*==hank1.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE HANK1(DKR,XJ,XN,XH1)
C     /****************************************************************/
C     # purpose      :                                                 *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MLP,CZERO,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DKR
      COMPLEX*16 XH1(MLP),XJ(MLP),XN(MLP)
C
C Local variables
C
      COMPLEX*16 ALPHA,XJ0,XJ1,XJ2
      INTEGER L,L1,L2,LMAX
C
C*** End of declarations rewritten by SPAG
C
      LMAX = INT(1.05*CDABS(DKR)+25.D0)
      LMAX = MIN(248,LMAX)
C
      XJ2 = CZERO
      XJ1 = DCMPLX(1.D-29,0.D0)
      DO L = 1,LMAX
         L1 = LMAX + 1 - L
         XJ0 = DBLE(2*L1+1)*XJ1/DKR - XJ2
         IF ( L1.LE.MLP ) XJ(L1) = XJ0
         XJ2 = XJ1
         XJ1 = XJ0
      END DO
C
      ALPHA = CDSIN(DKR)/(DKR*XJ(1))
      DO L = 1,MLP
         XJ(L) = ALPHA*XJ(L)
      END DO
      XN(1) = -CDCOS(DKR)/DKR
      XN(2) = (-CDSIN(DKR)/DKR) - CDCOS(DKR)/DKR**2
      L2 = MLP - 2
      DO L = 1,L2
         XN(L+2) = DBLE(2*L+1)*XN(L+1)/DKR - XN(L)
      END DO
      DO L = 1,MLP
         XH1(L) = XJ(L) + CIMAG*XN(L)
      END DO
C
      END
C*==hank2.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE HANK2(Z,JX)
C     /****************************************************************/
C     # purpose      :                                                 *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MLP,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 Z
      COMPLEX*16 JX(MLP)
C
C Local variables
C
      INTEGER I,J
      REAL*8 I1,JFAC,L1,SIGNVAR
      COMPLEX*16 Z2,ZS,ZSUM
C
C*** End of declarations rewritten by SPAG
C
      L1 = 1.0D0
      IF ( CDABS(Z).LT.1.D-6 ) THEN
         JX(1) = CONE
         DO I = 2,MLP
            JX(I) = CZERO
         END DO
      ELSE
         DO I = 1,MLP
            ZSUM = CONE
            Z2 = 0.5D0*Z*Z
            SIGNVAR = -1.0D0
            I1 = 2.0D0*DBLE(I-1) + 3.0D0
            J = 1
            JFAC = 1.0D0
 20         CONTINUE
            ZS = Z2**J/(I1*JFAC)
            IF ( CDABS(ZS).LT.1.0D-6 ) THEN
               JX(I) = (Z**(I-1)/L1)*ZSUM
               L1 = L1*(2.0D0*DBLE(I)+1.0D0)
            ELSE
               ZSUM = ZSUM + SIGNVAR*ZS
               SIGNVAR = -SIGNVAR
               I1 = I1*(2.0D0*DBLE(I-1)+2.0D0*DBLE(J)+3.0D0)
               J = J + 1
               JFAC = JFAC*DBLE(J)
               GOTO 20
            END IF
         END DO
      END IF
C
      END
C*==ilmkms.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE ILMKMS()
      USE MOD_SPEC,ONLY:NATLM,ML,MLQNAT
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,KME,KMO,SW,TAUW,LKAE,LKAO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER EVEN,EVPOS,IAN,KAI,L,M,MAXL,MUI,ODPOS
      REAL*8 KA,MU
C
C*** End of declarations rewritten by SPAG
C
      KME = 0.0D0
      KMO = 0.0D0
      EVPOS = 1
      ODPOS = 1
      MAXL = ML - 1
      DO IAN = 1,NATLM
         DO L = 0,MAXL
            DO M = -L,L
               IF ( MOD(L+M,2).EQ.1 ) THEN
                  LMO(ODPOS,1) = L
                  LMO(ODPOS,2) = M
                  ODPOS = ODPOS + 1
               ELSE
                  LME(EVPOS,1) = L
                  LME(EVPOS,2) = M
                  EVPOS = EVPOS + 1
               END IF
            END DO
         END DO
      END DO
      EVPOS = 1
      ODPOS = 1
      DO IAN = 1,NATLM
         DO KAI = (-MAXL) - 1,MAXL
            KA = DBLE(KAI)
            DO MUI = 1,2*ABS(KAI)
               MU = ((-ABS(KA))) - 0.5D0 + DBLE(MUI)
               IF ( MOD(IDINT(ABS(KA+MU+0.5D0)),2).EQ.1 ) THEN
                  EVEN = 0
               ELSE
                  EVEN = 1
               END IF
               IF ( EVEN.EQ.1 .AND. KA.GT.0.D0 .OR. EVEN.EQ.0 .AND. 
     &              KA.LT.0.D0 ) THEN
                  KMO(ODPOS,1) = KA
                  KMO(ODPOS,2) = MU
                  ODPOS = ODPOS + 1
               ELSE
                  KME(EVPOS,1) = KA
                  KME(EVPOS,2) = MU
                  EVPOS = EVPOS + 1
               END IF
            END DO
         END DO
      END DO
      SW(1) = 0.5D0
      SW(2) = -0.5D0
      TAUW(1) = -1.D0
      TAUW(2) = 1.D0
      DO KAI = 1,MLQNAT
         KA = KME(KAI,1)
         IF ( KA.LT.0.D0 ) THEN
            LKAE(KAI) = ((-IDINT(KA))) - 1
         ELSE
            LKAE(KAI) = IDINT(KA)
         END IF
      END DO
      DO KAI = 1,MLQNAT
         KA = KMO(KAI,1)
         IF ( KA.LT.0.D0 ) THEN
            LKAO(KAI) = ((-IDINT(KA))) - 1
         ELSE
            LKAO(KAI) = IDINT(KA)
         END IF
      END DO
      END
C*==inpcnt.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE INPCNT(ISEQ,LAYER1,GANZ,LAYS,LAYB,LAYP)
C
      USE MOD_SPEC,ONLY:LAYSM,LL
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,LAYB,LAYER1,LAYP,LAYS
      INTEGER ISEQ(0:LAYSM)
C
C Local variables
C
      INTEGER I,IHE,IS,LAYER,LAYSEQ(0:LL+1),NOB
C
C*** End of declarations rewritten by SPAG
C
C     nob:maximum number of beams used in expansion of plane waves
C
      NOB = GANZ
      IF ( IP.GT.0 ) WRITE (NOUT1,99002) NOB
C
C     layseq : scattering types of the photoemitting layers
C
      LAYER = -1
      DO IS = 0,LAYB - 1
         LAYER = LAYER + 1
         LAYSEQ(LAYER) = ISEQ(IS)
      END DO
      LAYER = LAYB - 1
 100  CONTINUE
      DO IHE = 1,LAYER1
         DO IS = LAYB,LAYS
            LAYER = LAYER + 1
            IF ( LAYER.GT.LAYP ) GOTO 200
            LAYSEQ(LAYER) = ISEQ(IS)
         END DO
      END DO
      GOTO 100
 200  CONTINUE
      IF ( IP.GT.0 ) WRITE (NOUT1,99001) (I,LAYSEQ(I),I=0,LAYP)
      RETURN
C
99001 FORMAT (1x,'photo emitting layer number=',i3,' is of type=',i3)
99002 FORMAT (1x,'beams ',i3,' to be used',
     &        ' in plane wave expansion in interstitial region')
      END
C*==intpol.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE INTPOL(INTG,POLG,SG,INTFAK,POL0G,INTGPM,GANZ)
      USE MOD_SPEC,ONLY:LG,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ
      COMPLEX*16 INTFAK(LG),INTG(LG),INTGPM(LG),POLG(LG,3),SG(LG,2,2)
      REAL*8 POL0G(LG,3)
C
C Local variables
C
      INTEGER DIM2,G,I,J,W
      REAL*8 POL0H(3)
      COMPLEX*16 RHO0(2,2),RHOG(2,2),RHOH(2,2),S(2,2),SKR(2,2),SPRHOG
C
C*** End of declarations rewritten by SPAG
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C rhoset    mult                                              *
C**************************************************************
C
C
      DO I = 1,2
         DO J = 1,2
            S(I,J) = CZERO
            SKR(I,J) = CZERO
         END DO
      END DO
C
      DO G = 1,LG
         INTG(G) = CZERO
         INTGPM(G) = CZERO
         DO I = 1,3
            POLG(G,I) = CZERO
         END DO
      END DO
C
      DO I = 1,3
         POL0H(I) = 0.0D0
      END DO
C
      DIM2 = 2
      DO G = 1,LG
         IF ( G.LE.GANZ ) THEN
            S(1,1) = SG(G,1,1)
            S(1,2) = SG(G,1,2)
            S(2,1) = SG(G,2,1)
            S(2,2) = SG(G,2,2)
            SKR(1,1) = DCONJG(S(1,1))
            SKR(1,2) = DCONJG(S(2,1))
            SKR(2,1) = DCONJG(S(1,2))
            SKR(2,2) = DCONJG(S(2,2))
            DO W = 1,2
               IF ( W.EQ.1 ) THEN
                  POL0H(1) = POL0G(G,1)
                  POL0H(2) = POL0G(G,2)
                  POL0H(3) = POL0G(G,3)
               ELSE IF ( W.EQ.2 ) THEN
                  POL0H(1) = -POL0G(G,1)
                  POL0H(2) = -POL0G(G,2)
                  POL0H(3) = -POL0G(G,3)
               END IF
               CALL RHOSET(RHO0,POL0H)
               CALL MULT(RHO0,SKR,RHOH,DIM2)
               CALL MULT(S,RHOH,RHOG,DIM2)
               SPRHOG = RHOG(1,1) + RHOG(2,2)
               IF ( W.EQ.1 ) THEN
                  INTG(G) = INTFAK(G)*SPRHOG
               ELSE IF ( W.EQ.2 ) THEN
                  INTGPM(G) = INTFAK(G)*SPRHOG
               END IF
               IF ( W.EQ.1 ) THEN
                  IF ( CDABS(SPRHOG).GT.1.0D-12 ) THEN
                     POLG(G,1) = (RHOG(2,1)+RHOG(1,2))/SPRHOG
                     POLG(G,2) = (DCMPLX(0.0D0,1.0D0)*(RHOG(1,2)-RHOG(2,
     &                           1)))/SPRHOG
                     POLG(G,3) = (RHOG(1,1)-RHOG(2,2))/SPRHOG
                  END IF
               END IF
            END DO
         END IF
      END DO
      END
C*==iposi.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE IPOSI
      USE MOD_SPEC,ONLY:NATLM,ML,MLQNAT
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,KME,KMO,SW,LKAE,LKAO
      USE MOD_SPEC_POSI,ONLY:PYLME,PYKMSE,PHLME,PXE,PTKE1,PTKE2,PYLMO,
     &    PYKMSO,PHLMO,PXO,PTKO1,PTKO2,PYLMME,PYLMMO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLR,MLT,MLRNAT,MLTNAT
      PARAMETER (MLR=(ML+1)*ML/2,MLT=ML*(ML-1)/2,MLRNAT=MLR*NATLM,
     &           MLTNAT=MLT*NATLM)
C
C Local variables
C
      INTEGER KA,KM,L,LM,M,ME,MO,S
      INTEGER SUCHE
C
C*** End of declarations rewritten by SPAG
C
C Local variables
C
C
C PARAMETER definitions
C
C
C*** End of declarations rewritten by SPAG
C
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C suche                                                       *
C**************************************************************
C
      IF ( .NOT.ALLOCATED(PHLME) ) THEN
         ALLOCATE (PHLME(MLRNAT),PHLMO(MLTNAT),PTKE1(MLQNAT),
     &             PTKE2(MLQNAT))
         ALLOCATE (PTKO1(MLQNAT),PTKO2(MLQNAT),PXE(MLQNAT,2),
     &             PXO(MLQNAT,2))
         ALLOCATE (PYKMSE(MLQNAT,2),PYKMSO(MLQNAT,2),PYLME(MLRNAT))
         ALLOCATE (PYLMME(MLRNAT),PYLMMO(MLTNAT),PYLMO(MLTNAT))
      END IF
C
      DO LM = 1,MLRNAT
         PYLME(LM) = LME(LM,2) + 1 + LME(LM,1)*(LME(LM,1)+1)
         PYLMME(LM) = ((-LME(LM,2))) + 1 + LME(LM,1)*(LME(LM,1)+1)
         PHLME(LM) = LME(LM,1) + 1
      END DO
      DO LM = 1,MLTNAT
         PYLMO(LM) = LMO(LM,2) + 1 + LMO(LM,1)*(LMO(LM,1)+1)
         PYLMMO(LM) = ((-LMO(LM,2))) + 1 + LMO(LM,1)*(LMO(LM,1)+1)
         PHLMO(LM) = LMO(LM,1) + 1
      END DO
      DO KM = 1,MLQNAT
         DO S = 1,2
            L = LKAE(KM)
            M = IDINT(KME(KM,2)-SW(S))
            IF ( ABS(M).LE.L ) THEN
               PYKMSE(KM,S) = M + 1 + L*(L+1)
            ELSE
               PYKMSE(KM,S) = 0
            END IF
         END DO
      END DO
      DO KM = 1,MLQNAT
         DO S = 1,2
            L = LKAO(KM)
            M = IDINT(KMO(KM,2)-SW(S))
            IF ( ABS(M).LE.L ) THEN
               PYKMSO(KM,S) = M + 1 + L*(L+1)
            ELSE
               PYKMSO(KM,S) = 0
            END IF
         END DO
      END DO
      DO KM = 1,MLQNAT
         L = LKAE(KM)
         ME = IDINT(KME(KM,2)-SW(1))
         MO = IDINT(KME(KM,2)-SW(2))
         IF ( ABS(ME).LE.L ) THEN
            PXE(KM,2) = SUCHE(L,ME,LME,MLRNAT)
         ELSE
            PXE(KM,2) = 0
         END IF
         IF ( ABS(MO).LE.L ) THEN
            PXO(KM,1) = SUCHE(L,MO,LMO,MLTNAT)
         ELSE
            PXO(KM,1) = 0
         END IF
      END DO
      DO KM = 1,MLQNAT
         L = LKAO(KM)
         ME = IDINT(KMO(KM,2)-SW(2))
         MO = IDINT(KMO(KM,2)-SW(1))
         IF ( ABS(ME).LE.L ) THEN
            PXE(KM,1) = SUCHE(L,ME,LME,MLRNAT)
         ELSE
            PXE(KM,1) = 0
         END IF
         IF ( ABS(MO).LE.L ) THEN
            PXO(KM,2) = SUCHE(L,MO,LMO,MLTNAT)
         ELSE
            PXO(KM,2) = 0
         END IF
      END DO
      DO KM = 1,MLQNAT
         KA = IDINT(KME(KM,1))
         IF ( KA.GT.0 ) THEN
            PTKE1(KM) = KA + 1
            PTKE2(KM) = 2
         ELSE
            PTKE1(KM) = -KA
            PTKE2(KM) = 1
         END IF
      END DO
      DO KM = 1,MLQNAT
         KA = IDINT(KMO(KM,1))
         IF ( KA.GT.0 ) THEN
            PTKO1(KM) = KA + 1
            PTKO2(KM) = 2
         ELSE
            PTKO1(KM) = -KA
            PTKO2(KM) = 1
         END IF
      END DO
      END
C*==ixcon.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE IXCON(CL,CLX,KMW,SW)
C**************************************************************
C                                                             *
C subroutines and functions called from this routine          *
C                                                             *
C clegor                                                      *
C                                                             *
C**************************************************************
C
      USE MOD_SPEC,ONLY:MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CL(MLQNAT,2),CLX(MLQNAT,MLQNAT,2)
      REAL*8 KMW(MLQNAT,2),SW(2)
C
C Local variables
C
      REAL*8 CLEGOR
      INTEGER KM,KMS,S
C
C*** End of declarations rewritten by SPAG
C
      DO KM = 1,MLQNAT
         DO S = 1,2
            CL(KM,S) = DCMPLX(CLEGOR(KMW(KM,1),KMW(KM,2),SW(S)),0.D0)
         END DO
      END DO
      DO S = 1,2
         DO KM = 1,MLQNAT
            DO KMS = 1,MLQNAT
               CLX(KM,KMS,S) = CL(KM,S)*CL(KMS,S)
            END DO
         END DO
      END DO
      END
C*==ixcona.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE IXCONA(CLX1,CLX2,CLX,KMW1,SW,KMW2)
C**************************************************************
C                                                             *
C subroutines and functions called from this routine          *
C                                                             *
C clegor                                                      *
C                                                             *
C**************************************************************
C
      USE MOD_SPEC,ONLY:MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CLX(MLQNAT,MLQNAT,2),CLX1(MLQNAT,2),CLX2(MLQNAT,2)
      REAL*8 KMW1(MLQNAT,2),KMW2(MLQNAT,2),SW(2)
C
C Local variables
C
      REAL*8 CLEGOR
      INTEGER KM,KMS,S
C
C*** End of declarations rewritten by SPAG
C
      DO KM = 1,MLQNAT
         DO S = 1,2
            CLX1(KM,S) = DCMPLX(CLEGOR(KMW1(KM,1),KMW1(KM,2),SW(S)),
     &                   0.D0)
         END DO
      END DO
      DO KM = 1,MLQNAT
         DO S = 1,2
            CLX2(KM,S) = DCMPLX(CLEGOR(KMW2(KM,1),KMW2(KM,2),SW(S)),
     &                   0.D0)
         END DO
      END DO
      DO S = 1,2
         DO KM = 1,MLQNAT
            DO KMS = 1,MLQNAT
               CLX(KM,KMS,S) = CLX1(KM,S)*CLX2(KMS,S)
            END DO
         END DO
      END DO
      END
C*==k0set.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE K0SET(K0,K03,K0BTRS,ME,MO,MSE,MSO,MCON,KGFAK1,GN,GINT,
     &                 E,THETA,PHI,CLE,CLO,GANZ,NREL,LHL,TEST,XED,XOD,
     &                 INTFAK,KGM,CLIGHT,IBLOCH,VX,VY,VS,NNSR,NSR,RLVX,
     &                 RLVY,RLVS,NNSK,NSK,KGFAK2,BRUSEP,NEV,NOD,IT,
     &                 KGZBLK,KGT)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C realvec   rezvec    sphrm4                                  *
C**************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,MLQ,NPM,LG,MLQNAT,LH,PI,CZERO,CIMAG
      USE MOD_SPEC_GEOM,ONLY:SEP,TV
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,TAUW,LKAE,LKAO
      USE MOD_SPEC_POSI,ONLY:PYLME,PYKMSE,PYLMO,PYKMSO,PYLMME,PYLMMO
      USE MOD_SPEC_QXQY,ONLY:QX,QY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT,PHI,TEST,THETA
      COMPLEX*16 E,K03,K0BTRS
      INTEGER GANZ,IBLOCH,IT,LHL,NEV,NNSK,NNSR,NOD,NREL,XED,XOD
      REAL*8 BRUSEP(3),GN(LG,2),K0(2),KGT(2,LG),RLVS(NPM),RLVX(NPM),
     &       RLVY(NPM),VS(NPM),VX(NPM),VY(NPM)
      COMPLEX*16 CLE(MLQNAT,2),CLO(MLQNAT,2),INTFAK(LG),
     &           KGFAK1(LG,2,LAYSM),KGFAK2(LG,2,LAYSM),KGM(LG),
     &           KGZBLK(LH),MCON(LG),ME(LG,2,MLQNAT,2),MO(LG,2,MLQNAT,2)
     &           ,MSE(LG,2,MLQNAT,2),MSO(LG,2,MLQNAT,2)
      INTEGER GINT(NPM,2),NSK(NPM),NSR(NPM)
C
C Local variables
C
      COMPLEX*16 CI,KGT3,KGTAS1,KGTAS2,KGZ,PI8SQR,PICON,VORZ,YLM(:),
     &           YLMM(:)
      INTEGER G,I,KM,L,LM,S,TAU
      REAL*8 TEST1,XXX,YYY
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE YLM,YLMM
      ALLOCATE (YLM(MLQ),YLMM(MLQ))
C
      CI = CIMAG
      PI8SQR = DCMPLX(8.D0*PI*PI,0.D0)
      K0BTRS = 2.D0*E
      IF ( NREL.EQ.1 ) K0BTRS = 2.D0*E + (E*E)/CLIGHT**2
      IF ( NREL.EQ.0 ) THEN
         TEST1 = 2.D0*TEST
      ELSE
         TEST1 = 2.D0*TEST + TEST*TEST/(CLIGHT*CLIGHT)
      END IF
      IF ( LHL.EQ.1 ) THEN
         K0(1) = -SQRT(TEST1)*COS(PHI)*SIN(THETA)
         K0(2) = -SQRT(TEST1)*SIN(PHI)*SIN(THETA)
      ELSE IF ( LHL.EQ.2 ) THEN
C        K0(1) = SQRT(TEST1)*COS(PHI)*SIN(THETA)
C        K0(2) = SQRT(TEST1)*SIN(PHI)*SIN(THETA)
         K0(1) = SQRT(TEST1)*COS(PHI)*SIN(THETA) - QX
         K0(2) = SQRT(TEST1)*SIN(PHI)*SIN(THETA) - QY
C         CALL KFOLD(K0,K1,BBAS)
C         K0 = K1
      END IF
      IF ( IBLOCH.EQ.0 ) THEN
         K0(1) = SQRT(TEST1)*COS(PHI)*SIN(THETA)
         K0(2) = SQRT(TEST1)*SIN(PHI)*SIN(THETA)
      END IF
      K03 = CDSQRT(K0BTRS-(K0(1)**2+K0(2)**2))
      PICON = PI8SQR/(CDSQRT(K0BTRS)*TV)
C
C     selecting real*8 and reciprocal lattice vectors
C
      CALL REALVEC(VX,VY,VS,NNSR,NSR)
      CALL REZVEC(RLVX,RLVY,RLVS,NNSK,NSK,GINT)
C
      DO I = 1,LG
         GN(I,1) = RLVX(I)
         GN(I,2) = RLVY(I)
      END DO
C
      DO G = 1,GANZ
         KGT(1,G) = K0(1) + GN(G,1)
         KGT(2,G) = K0(2) + GN(G,2)
         KGZ = CDSQRT(K0BTRS-(KGT(1,G)**2+KGT(2,G)**2))
         KGZBLK(G) = KGZ
         KGZBLK(G+GANZ) = KGZBLK(G)
         KGM(G) = -KGZ
         MCON(G) = PICON/KGZ
         INTFAK(G) = DCMPLX(DBLE(KGZ)/DBLE(K03),0.D0)
         DO TAU = 1,2
            VORZ = DCMPLX(TAUW(TAU),0.D0)
            KGT3 = VORZ*KGZ
            KGTAS1 = KGT(1,G)*DCMPLX(SEP(2,IT),0.D0) + KGT(2,G)
     &               *DCMPLX(SEP(3,IT),0.D0)
     &               + KGT3*DCMPLX(SEP(1,IT),0.D0)
            KGTAS2 = KGT(1,G)*DCMPLX(BRUSEP(2),0.D0) + KGT(2,G)
     &               *DCMPLX(BRUSEP(3),0.D0)
     &               + KGT3*DCMPLX(BRUSEP(1),0.D0)
            KGFAK1(G,TAU,IT) = CDEXP(VORZ*CI*KGTAS1)
            KGFAK2(G,TAU,IT) = CDEXP(VORZ*CI*KGTAS2)
            XXX = KGT(1,G)
            YYY = KGT(2,G)
            CALL SPHRM4(XXX,YYY,KGT3,YLM)
            DO LM = 1,NEV
               YLMM(PYLME(LM)) = DCMPLX((-1.D0)**LME(LM,2),0.D0)
     &                           *YLM(PYLMME(LM))
            END DO
            DO LM = 1,NOD
               YLMM(PYLMO(LM)) = DCMPLX((-1.D0)**LMO(LM,2),0.D0)
     &                           *YLM(PYLMMO(LM))
            END DO
            IF ( NREL.EQ.0 ) THEN
               DO LM = 1,NEV
                  L = LME(LM,1)
                  ME(G,TAU,LM,1) = CI**L*YLMM(PYLME(LM))
                  MSE(G,TAU,LM,1) = CI**(-L)*YLM(PYLME(LM))
               END DO
               DO LM = 1,NOD
                  L = LMO(LM,1)
                  MO(G,TAU,LM,1) = CI**L*YLMM(PYLMO(LM))
                  MSO(G,TAU,LM,1) = CI**(-L)*YLM(PYLMO(LM))
               END DO
            ELSE IF ( NREL.EQ.1 ) THEN
               DO S = 1,2
                  DO KM = 1,XED
                     IF ( PYKMSE(KM,S).NE.0 ) THEN
                        ME(G,TAU,KM,S) = CI**LKAE(KM)*CLE(KM,S)
     &                     *YLMM(PYKMSE(KM,S))
                        MSE(G,TAU,KM,S) = CI**(-LKAE(KM))*CLE(KM,S)
     &                     *YLM(PYKMSE(KM,S))
                     ELSE
                        ME(G,TAU,KM,S) = CZERO
                        MSE(G,TAU,KM,S) = CZERO
                     END IF
                  END DO
                  DO KM = 1,XOD
                     IF ( PYKMSO(KM,S).NE.0 ) THEN
                        MO(G,TAU,KM,S) = CI**LKAO(KM)*CLO(KM,S)
     &                     *YLMM(PYKMSO(KM,S))
                        MSO(G,TAU,KM,S) = CI**(-LKAO(KM))*CLO(KM,S)
     &                     *YLM(PYKMSO(KM,S))
                     ELSE
                        MO(G,TAU,KM,S) = CZERO
                        MSO(G,TAU,KM,S) = CZERO
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
      END
C*==kgtset.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE KGTSET(K0,GN,GANZ,KGT,EVACH,TAUW,THETAS,PHIS,ESCAPE)
      USE MOD_SPEC,ONLY:LG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EVACH
      INTEGER GANZ
      INTEGER ESCAPE(LG)
      REAL*8 GN(LG,2),K0(2),KGT(LG,2,3),PHIS(LG),TAUW(2),THETAS(LG)
C
C Local variables
C
      REAL*8 ARCTAN
      REAL*8 ARG
      INTEGER G,TAU
C
C*** End of declarations rewritten by SPAG
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C arctan                                                      *
C**************************************************************
C
      DO G = 1,LG
         DO TAU = 1,2
            KGT(G,TAU,1) = 0.0
            KGT(G,TAU,2) = 0.0
            KGT(G,TAU,3) = 0.0
         END DO
      END DO
C
      DO G = 1,LG
         IF ( G.LE.GANZ ) THEN
            DO TAU = 1,2
               KGT(G,TAU,1) = K0(1) + GN(G,1)
               KGT(G,TAU,2) = K0(2) + GN(G,2)
               ARG = 2*EVACH - KGT(G,TAU,1)**2 - KGT(G,TAU,2)**2
               IF ( ARG.GE.0.D0 ) THEN
                  KGT(G,TAU,3) = TAUW(TAU)*SQRT(ARG)
                  ESCAPE(G) = 1
               ELSE
                  KGT(G,TAU,3) = TAUW(TAU)*SQRT((-ARG))
                  ESCAPE(G) = 0
               END IF
            END DO
            IF ( ARG.GE.0.D0 ) THEN
               THETAS(G) = ACOS(SQRT(ARG)/SQRT(2*EVACH))
C               IF ( KGT(G,1,1).NE.0.D0 ) THEN
               IF ( ABS(KGT(G,1,1)).GT.1.0D-16 ) THEN
                  PHIS(G) = ARCTAN(KGT(G,1,1),KGT(G,1,2))
               ELSE IF ( KGT(G,1,2).GT.0.D0 ) THEN
                  PHIS(G) = 2.D0*ATAN(1.D0)
               ELSE IF ( KGT(G,1,2).LT.0.D0 ) THEN
                  PHIS(G) = -2.D0*ATAN(1.D0)
               ELSE
                  PHIS(G) = 0.D0
               END IF
            ELSE
               THETAS(G) = 0.D0
               PHIS(G) = 0.D0
            END IF
         ELSE
            THETAS(G) = 0.D0
            PHIS(G) = 0.D0
         END IF
      END DO
      END
C*==arctan.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
C     ***************************************************************
C
      REAL*8 FUNCTION ARCTAN(X,Y)
      USE MOD_SPEC,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 X,Y
C
C*** End of declarations rewritten by SPAG
C
      IF ( X.GE.0.D0 .AND. Y.GE.0.D0 ) ARCTAN = DATAN(Y/X)
      IF ( X.LT.0.D0 .AND. Y.GE.0.D0 ) ARCTAN = DATAN(Y/X) + PI
      IF ( X.LT.0.D0 .AND. Y.LT.0.D0 ) ARCTAN = DATAN(Y/X) - PI
      IF ( X.GE.0.D0 .AND. Y.LT.0.D0 ) ARCTAN = DATAN(Y/X)
      END
C*==mmat.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE MMAT(M,TAU,TAUS,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,
     &                XOD,QD,SD,GANZ,EPOSM,EPOSP,IT,NATL,IR,SUMEO,SUMOE,
     &                IBXY,XEEH,XEOH,XOEH,XOOH)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,LG,MLQNAT,LH,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,IBXY,IR,IT,QD,SD,TAU,TAUS,XED,XOD
      COMPLEX*16 EPOSM(NATLM,LH,LAYSM),EPOSP(NATLM,LH,LAYSM),M(QD,QD),
     &           MCON(LG),ME(LG,2,MLQNAT,2),MO(LG,2,MLQNAT,2),
     &           MSE(LG,2,MLQNAT,2),MSO(LG,2,MLQNAT,2),SUME(MLQNAT),
     &           SUMEO(MLQNAT),SUMO(MLQNAT),SUMOE(MLQNAT),XEEH(XED,XED),
     &           XEOH(XED,XOD),XOEH(XOD,XED),XOOH(XOD,XOD)
      INTEGER NATL(LAYSM),PQ(LG,2)
C
C Local variables
C
      INTEGER G,GS,I,IANF,IEND,IG,IGS,INAT,IS,NAT,S,SS,XEDH,XODH
      COMPLEX*16 KRON,SH,SUMSE,SUMSEO,SUMSO,SUMSOE
C
C*** End of declarations rewritten by SPAG
C
      NAT = NATL(IT)
      XEDH = XED/NAT
      XODH = XOD/NAT
      IGS = 1
      DO GS = 1,GANZ
         DO SS = 1,SD
            DO I = 1,XED
               SUME(I) = CZERO
               DO INAT = 1,NAT
                  IANF = 1 + XEDH*(INAT-1)
                  IEND = XEDH*INAT
                  SH = CZERO
                  DO IS = IANF,IEND
                     SH = SH + MSE(GS,TAUS,IS,SS)*XEEH(I,IS)
                  END DO
                  IF ( TAUS.EQ.1 ) THEN
                     SUME(I) = SUME(I) + SH/EPOSM(INAT,IGS,IT)
                  ELSE IF ( TAUS.EQ.2 ) THEN
                     SUME(I) = SUME(I) + SH/EPOSP(INAT,IGS,IT)
                  END IF
               END DO
            END DO
            IF ( (IR.GT.0) .OR. (IBXY.NE.0) ) THEN
               DO I = 1,XOD
                  SUMOE(I) = CZERO
                  DO INAT = 1,NAT
                     IANF = 1 + XEDH*(INAT-1)
                     IEND = XEDH*INAT
                     SH = CZERO
                     DO IS = IANF,IEND
                        SH = SH + MSE(GS,TAUS,IS,SS)*XOEH(I,IS)
                     END DO
                     IF ( TAUS.EQ.1 ) THEN
                        SUMOE(I) = SUMOE(I) + SH/EPOSM(INAT,IGS,IT)
                     ELSE IF ( TAUS.EQ.2 ) THEN
                        SUMOE(I) = SUMOE(I) + SH/EPOSP(INAT,IGS,IT)
                     END IF
                  END DO
               END DO
            END IF
            DO I = 1,XOD
               SUMO(I) = CZERO
               DO INAT = 1,NAT
                  IANF = 1 + XODH*(INAT-1)
                  IEND = XODH*INAT
                  SH = CZERO
                  DO IS = IANF,IEND
                     SH = SH + MSO(GS,TAUS,IS,SS)*XOOH(I,IS)
                  END DO
                  IF ( TAUS.EQ.1 ) THEN
                     SUMO(I) = SUMO(I) + SH/EPOSM(INAT,IGS,IT)
                  ELSE IF ( TAUS.EQ.2 ) THEN
                     SUMO(I) = SUMO(I) + SH/EPOSP(INAT,IGS,IT)
                  END IF
               END DO
            END DO
            IF ( (IR.GT.0) .OR. (IBXY.NE.0) ) THEN
               DO I = 1,XED
                  SUMEO(I) = CZERO
                  DO INAT = 1,NAT
                     IANF = 1 + XODH*(INAT-1)
                     IEND = XODH*INAT
                     SH = CZERO
                     DO IS = IANF,IEND
                        SH = SH + MSO(GS,TAUS,IS,SS)*XEOH(I,IS)
                     END DO
                     IF ( TAUS.EQ.1 ) THEN
                        SUMEO(I) = SUMEO(I) + SH/EPOSM(INAT,IGS,IT)
                     ELSE IF ( TAUS.EQ.2 ) THEN
                        SUMEO(I) = SUMEO(I) + SH/EPOSP(INAT,IGS,IT)
                     END IF
                  END DO
               END DO
            END IF
            IG = 1
            DO G = 1,GANZ
               DO S = 1,SD
                  IF ( G.EQ.GS .AND. S.EQ.SS .AND. TAU.EQ.TAUS ) THEN
                     KRON = CONE
                  ELSE
                     KRON = CZERO
                  END IF
                  SUMSE = CZERO
                  DO INAT = 1,NAT
                     IANF = 1 + XEDH*(INAT-1)
                     IEND = XEDH*INAT
                     SH = CZERO
                     DO I = IANF,IEND
                        SH = SH + ME(G,TAU,I,S)*SUME(I)
                     END DO
                     IF ( TAU.EQ.1 ) THEN
                        SUMSE = SUMSE + SH*EPOSM(INAT,IG,IT)
                     ELSE IF ( TAU.EQ.2 ) THEN
                        SUMSE = SUMSE + SH*EPOSP(INAT,IG,IT)
                     END IF
                  END DO
                  IF ( (IR.GT.0) .OR. (IBXY.NE.0) ) THEN
                     SUMSEO = CZERO
                     DO INAT = 1,NAT
                        IANF = 1 + XEDH*(INAT-1)
                        IEND = XEDH*INAT
                        SH = CZERO
                        DO I = IANF,IEND
                           SH = SH + ME(G,TAU,I,S)*SUMEO(I)
                        END DO
                        IF ( TAU.EQ.1 ) THEN
                           SUMSEO = SUMSEO + SH*EPOSM(INAT,IG,IT)
                        ELSE IF ( TAU.EQ.2 ) THEN
                           SUMSEO = SUMSEO + SH*EPOSP(INAT,IG,IT)
                        END IF
                     END DO
                  END IF
                  SUMSO = CZERO
                  DO INAT = 1,NAT
                     IANF = 1 + XODH*(INAT-1)
                     IEND = XODH*INAT
                     SH = CZERO
                     DO I = IANF,IEND
                        SH = SH + MO(G,TAU,I,S)*SUMO(I)
                     END DO
                     IF ( TAU.EQ.1 ) THEN
                        SUMSO = SUMSO + SH*EPOSM(INAT,IG,IT)
                     ELSE IF ( TAU.EQ.2 ) THEN
                        SUMSO = SUMSO + SH*EPOSP(INAT,IG,IT)
                     END IF
                  END DO
                  IF ( (IR.GT.0) .OR. (IBXY.NE.0) ) THEN
                     SUMSOE = CZERO
                     DO INAT = 1,NAT
                        IANF = 1 + XODH*(INAT-1)
                        IEND = XODH*INAT
                        SH = CZERO
                        DO I = IANF,IEND
                           SH = SH + MO(G,TAU,I,S)*SUMOE(I)
                        END DO
                        IF ( TAU.EQ.1 ) THEN
                           SUMSOE = SUMSOE + SH*EPOSM(INAT,IG,IT)
                        ELSE IF ( TAU.EQ.2 ) THEN
                           SUMSOE = SUMSOE + SH*EPOSP(INAT,IG,IT)
                        END IF
                     END DO
                  END IF
                  M(PQ(GS,SS),PQ(G,S)) = KRON + MCON(GS)*(SUMSE+SUMSO)
                  IF ( (IR.GT.0) .OR. (IBXY.NE.0) ) M(PQ(GS,SS),PQ(G,S))
     &                 = M(PQ(GS,SS),PQ(G,S)) + MCON(GS)*(SUMSOE+SUMSEO)
                  IG = IG + 1
               END DO
            END DO
            IGS = IGS + 1
         END DO
      END DO
      END
C*==mmath.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE MMATH(C1,C4,GANZ,QD)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,QD
      COMPLEX*16 C1(QD,QD),C4(QD,QD)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,GANZ
         DO J = 1,GANZ
            C4(I,J) = C1(I,J)
            C4(I+GANZ,J+GANZ) = C1(I+GANZ,J+GANZ)
         END DO
      END DO
      DO I = 1,GANZ
         DO J = 1,GANZ
            C4(I,J+GANZ) = -1.D0*C1(I,J+GANZ)
            C4(I+GANZ,J) = -1.D0*C1(I+GANZ,J)
         END DO
      END DO
      END
C*==polkmp.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE POLKMP(POLG,INTG,SG,INTFAK,POL0,POL0G,POL0L,GANZ,TYP,
     &                  THETA,PHI,THETAS,PHIS,INTGPM,PK,PKGM,PN,PKGMN,
     &                  PKN)
      USE MOD_SPEC,ONLY:LG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,TYP
      REAL*8 PHI,THETA
      COMPLEX*16 INTFAK(LG),INTG(LG),INTGPM(LG),POLG(LG,3),SG(LG,2,2)
      REAL*8 PHIS(LG),PK(LG),PKGM(LG),PKGMN(LG),PKN(LG),PN(LG),POL0(3),
     &       POL0G(LG,3),POL0L(3),THETAS(LG)
C
C Local variables
C
      REAL*8 EK(3),EKGM(:,:),EKGMN(:,:),EKGMNB,EKHELP(:,:),EKN(:,:),
     &       EKNB,EN(:,:),ENB,POL(3)
      INTEGER G,I
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE EN,EKN,EKGM,EKGMN,EKHELP
      ALLOCATE (EN(LG,3),EKN(LG,3),EKGM(LG,3),EKGMN(LG,3))
      ALLOCATE (EKHELP(LG,3))
C
C*** End of declarations rewritten by SPAG
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C kreuz     intpol    dot                                     *
C**************************************************************
C
C     ek:    direction of the incomming k-vector
C     ekgm:  direction of the corresponding outgoing k-vector
C     en:    defines normal to the scattering plane (ek x ekgm)
C     ekgmn: ekgm x en
C     ekn:   ek x en
C     (ek,en,ekn) and (ekgm,en,ekgmn) define right-handed ons
C
      DO G = 1,LG
         EKHELP(G,1) = 0.0D0
         EKHELP(G,2) = 0.0D0
         EKHELP(G,3) = 0.0D0
         EKGM(G,1) = 0.0D0
         EKGM(G,2) = 0.0D0
         EKGM(G,3) = 0.0D0
         POL0G(G,1) = 0.0D0
         POL0G(G,2) = 0.0D0
         POL0G(G,3) = 0.0D0
         POL(1) = 0.0D0
         POL(2) = 0.0D0
         POL(3) = 0.0D0
      END DO
C
      EK(1) = -COS(PHI)*SIN(THETA)
      EK(2) = -SIN(PHI)*SIN(THETA)
      EK(3) = COS(THETA)
      DO G = 1,LG
         IF ( G.LE.GANZ ) THEN
            EKHELP(G,1) = EK(1)
            EKHELP(G,2) = EK(2)
            EKHELP(G,3) = EK(3)
            EKGM(G,1) = COS(PHIS(G))*SIN(THETAS(G))
            EKGM(G,2) = SIN(PHIS(G))*SIN(THETAS(G))
            EKGM(G,3) = -COS(THETAS(G))
            CALL KREUZ(EKHELP,EKGM,EN,G)
            CALL KREUZ(EKGM,EN,EKGMN,G)
            CALL KREUZ(EKHELP,EN,EKN,G)
            ENB = SQRT(EN(G,1)**2+EN(G,2)**2+EN(G,3)**2)
            EKGMNB = SQRT(EKGMN(G,1)**2+EKGMN(G,2)**2+EKGMN(G,3)**2)
            EKNB = SQRT(EKN(G,1)**2+EKN(G,2)**2+EKN(G,3)**2)
            DO I = 1,3
C               IF ( ENB.NE.0.D0 ) EN(G,I) = EN(G,I)/ENB
C               IF ( EKGMNB.NE.0.D0 ) EKGMN(G,I) = EKGMN(G,I)/EKGMNB
C               IF ( EKNB.NE.0.D0 ) EKN(G,I) = EKN(G,I)/EKNB
               IF ( ABS(ENB).GT.1.0D-16 ) EN(G,I) = EN(G,I)/ENB
               IF ( ABS(EKGMNB).GT.1.0D-16 ) EKGMN(G,I) = EKGMN(G,I)
     &              /EKGMNB
               IF ( ABS(EKNB).GT.1.0D-16 ) EKN(G,I) = EKN(G,I)/EKNB
            END DO
         END IF
      END DO
C
      DO G = 1,LG
C
         POL0G(G,1) = POL0(1)
         POL0G(G,2) = POL0(2)
         POL0G(G,3) = POL0(3)
C
C     transformation of the polarisation vector from laboratory
C     into crystal fixed coordinate system
C
C         IF (ABS(POL0(1))+ABS(POL0(2))+ABS(POL0(3)).EQ.0.0) THEN
         IF ( ABS(POL0(1))+ABS(POL0(2))+ABS(POL0(3)).LE.1.0D-16 ) THEN
            IF ( TYP.EQ.1 .OR. TYP.EQ.2 ) THEN
               IF ( G.LE.GANZ ) THEN
                  POL0G(G,1) = POL0L(1)*EN(G,1) + POL0L(2)*EKN(G,1)
     &                         + POL0L(3)*EK(1)
                  POL0G(G,2) = POL0L(1)*EN(G,2) + POL0L(2)*EKN(G,2)
     &                         + POL0L(3)*EK(2)
                  POL0G(G,3) = POL0L(1)*EN(G,3) + POL0L(2)*EKN(G,3)
     &                         + POL0L(3)*EK(3)
               END IF
            END IF
         END IF
      END DO
C
C     calculation of intensity and polarisation
C
      CALL INTPOL(INTG,POLG,SG,INTFAK,POL0G,INTGPM,GANZ)
      DO G = 1,LG
         IF ( G.LE.GANZ ) THEN
            POL(1) = DBLE(POLG(G,1))
            POL(2) = DBLE(POLG(G,2))
            POL(3) = DBLE(POLG(G,3))
            CALL DOT(POL,EKHELP,PK,G)
            CALL DOT(POL,EKGM,PKGM,G)
            CALL DOT(POL,EN,PN,G)
            CALL DOT(POL,EKGMN,PKGMN,G)
            CALL DOT(POL,EKN,PKN,G)
         END IF
      END DO
      END
C*==pqset.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE PQSET(PQ,GANZ)
      USE MOD_SPEC,ONLY:LG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ
      INTEGER PQ(LG,2)
C
C Local variables
C
      INTEGER G,S
C
C*** End of declarations rewritten by SPAG
C
      DO G = 1,LG
         DO S = 1,2
            PQ(G,S) = G + (S-1)*GANZ
         END DO
      END DO
      END
C*==prop.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE PROP(PR,TAU,KGFAK,PQ,QD,SD,GANZ,LAYS,ISEQ)
      USE MOD_SPEC,ONLY:LG,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,ISEQ,LAYS,QD,SD,TAU
      COMPLEX*16 KGFAK(LG,2,LAYS),PR(QD,QD)
      INTEGER PQ(LG,2)
C
C Local variables
C
      INTEGER G,GS,S,SS
C
C*** End of declarations rewritten by SPAG
C
C
      DO G = 1,GANZ
         DO GS = 1,GANZ
            DO S = 1,SD
               DO SS = 1,SD
                  IF ( G.EQ.GS .AND. S.EQ.SS ) THEN
                     PR(PQ(G,S),PQ(GS,SS)) = KGFAK(G,TAU,ISEQ)
                  ELSE
                     PR(PQ(G,S),PQ(GS,SS)) = CZERO
                  END IF
               END DO
            END DO
         END DO
      END DO
      END
C*==realvec.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE REALVEC(VX,VY,VS,NRR,NSR)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C relord                                                      *
C**************************************************************
      USE MOD_SPEC,ONLY:NPM
      USE MOD_SPEC_GEOM,ONLY:AR1,AR2
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NRR
      INTEGER NSR(NPM)
      REAL*8 VS(NPM),VX(NPM),VY(NPM)
C
C Local variables
C
      REAL*8 AN,AN1,AN2,TOL
      INTEGER I,I1,I2,IJ,IRING,J,N1,NA
C
C*** End of declarations rewritten by SPAG
C
      TOL = 1.D-3
      VX(1) = 0.D0
      VY(1) = 0.D0
      VS(1) = 0.D0
      IJ = 1
      N1 = 0
 100  CONTINUE
      N1 = N1 + 1
      NA = N1 + N1
      AN1 = DBLE(N1)
      AN2 = ((-AN1)) - 1.D0
      DO I1 = 1,NA
         AN2 = AN2 + 1.D0
         DO I2 = 1,4
            AN = AN1
            AN1 = -AN2
            AN2 = AN
            IJ = IJ + 1
            VX(IJ) = AN1*AR1(1) + AN2*AR2(1)
            VY(IJ) = AN1*AR1(2) + AN2*AR2(2)
            VS(IJ) = SQRT(VX(IJ)*VX(IJ)+VY(IJ)*VY(IJ))
            NRR = IJ
         END DO
      END DO
C
C     n1 = nnr = number of rings in the real*8 space lattice
C
      IF ( N1.LT.20 ) GOTO 100
      CALL RELORD(VX,VY,VS,NRR)
      J = 0
      IRING = 0
      DO I = 1,NRR - 1
         IRING = IRING + 1
         IF ( ABS(VS(I+1)-VS(I)).GT.TOL ) THEN
            J = J + 1
            NSR(J) = IRING
            IRING = 0
         END IF
      END DO
      DO I = 1,NRR
         IF ( IP.GT.3 ) WRITE (NOUT1,99001) I,VX(I),VY(I),VS(I)
      END DO
C
      RETURN
99001 FORMAT (1x,'real*8 gv',i4,2x,3F12.5)
      END
C*==reflex.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE REFLEX(EVACH,K0,RFXANZ,GRFX)
      USE MOD_SPEC_GEOM,ONLY:RAR
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MAXRFX
      PARAMETER (MAXRFX=100)
C
C Dummy arguments
C
      REAL*8 EVACH
      INTEGER RFXANZ
      INTEGER GRFX(MAXRFX,2)
      REAL*8 K0(2)
C
C Local variables
C
      REAL*8 G1,G2,TEST
      INTEGER I,IS,IZ,J,JS,JZ,LAENG,RFX
C
C*** End of declarations rewritten by SPAG
C
      LAENG = 2
 100  CONTINUE
      RFX = 0
      IS = 0
      JS = 0
      DO IZ = 0,LAENG
         I = IZ*(-1)**IZ
         IS = IS + I
         JS = 0
         DO JZ = 0,LAENG
            J = JZ*(-1)**JZ
            JS = JS + J
            G1 = IS*RAR(1) + JS*RAR(3)
            G2 = IS*RAR(2) + JS*RAR(4)
            TEST = 2*EVACH - (K0(1)+G1)**2 - (K0(2)+G2)**2
            IF ( TEST.GE.0.D0 ) THEN
               RFX = RFX + 1
               GRFX(RFX,1) = IS
               GRFX(RFX,2) = JS
            END IF
         END DO
      END DO
      IF ( RFX.EQ.RFXANZ ) THEN
         IF ( IP.GT.2 ) WRITE (NOUT1,99001)
         IF ( IP.GT.2 ) WRITE (NOUT1,99002)
     &                         (GRFX(I,1),GRFX(I,2),I=1,RFXANZ)
         RETURN
      ELSE
         RFXANZ = RFX
         LAENG = LAENG + 2
         GOTO 100
      END IF
99001 FORMAT (1x,'values of possible leed-reflexes')
99002 FORMAT (20I4)
      END
C*==relord.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RELORD(VX,VY,VS,N)
      USE MOD_SPEC,ONLY:NPM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 VS(NPM),VX(NPM),VY(NPM)
C
C Local variables
C
      INTEGER I,IM,J,K,M
      REAL*8 T
C
C*** End of declarations rewritten by SPAG
C
C     shell ordering routine for real*8 space lattice vectors
C
C
      M = N
 100  CONTINUE
      M = M/2
      K = N - M
      DO J = 1,K
         I = J
 150     CONTINUE
         IM = I + M
         IF ( VS(IM).LT.VS(I) ) THEN
            T = VS(I)
            VS(I) = VS(IM)
            VS(IM) = T
            T = VX(I)
            VX(I) = VX(IM)
            VX(IM) = T
            T = VY(I)
            VY(I) = VY(IM)
            VY(IM) = T
            I = I - M
            IF ( I.GE.1 ) GOTO 150
         END IF
      END DO
      IF ( M.GT.0 ) GOTO 100
      END
C*==rezord.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE REZORD(VX,VY,VS,N,GINT)
      USE MOD_SPEC,ONLY:NPM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      INTEGER GINT(NPM,2)
      REAL*8 VS(NPM),VX(NPM),VY(NPM)
C
C Local variables
C
      INTEGER I,IM,J,K,M,S1,S2
      REAL*8 T
C
C*** End of declarations rewritten by SPAG
C
C     shell ordering routine for the reciprocal lattice vectors
C
      M = N
 100  CONTINUE
      M = M/2
      K = N - M
      DO J = 1,K
         I = J
 150     CONTINUE
         IM = I + M
         IF ( VS(IM).LT.VS(I) ) THEN
            T = VS(I)
            VS(I) = VS(IM)
            VS(IM) = T
            T = VX(I)
            VX(I) = VX(IM)
            VX(IM) = T
            T = VY(I)
            VY(I) = VY(IM)
            VY(IM) = T
            S1 = GINT(I,1)
            S2 = GINT(I,2)
            GINT(I,1) = GINT(IM,1)
            GINT(I,2) = GINT(IM,2)
            GINT(IM,1) = S1
            GINT(IM,2) = S2
            I = I - M
            IF ( I.GE.1 ) GOTO 150
         END IF
      END DO
      IF ( M.GT.0 ) GOTO 100
      END
C*==rezvec.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE REZVEC(RLVX,RLVY,RLVS,NQQ,NSK,GINT)
C**************************************************************
C                                                             *
C subroutines and functions called from this routine          *
C                                                             *
C rezord                                                      *
C                                                             *
C**************************************************************
C
      USE MOD_SPEC,ONLY:NPM
      USE MOD_SPEC_GEOM,ONLY:RAR
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NQQ
      INTEGER GINT(NPM,2),NSK(NPM)
      REAL*8 RLVS(NPM),RLVX(NPM),RLVY(NPM)
C
C Local variables
C
      REAL*8 AB1,AB2,AN,AN1,AN2,TOL
      INTEGER I,I1,I2,II,IJ,IRING,J,N1,NA
C
C*** End of declarations rewritten by SPAG
C
      TOL = 1.D-3
      II = 1
      IJ = 0
      N1 = -1
 100  CONTINUE
      N1 = N1 + 1
      NA = N1 + N1 + II
      AN1 = DBLE(N1)
      AN2 = ((-AN1)) - 1.D0
      DO I1 = 1,NA
         AN2 = AN2 + 1.D0
         DO I2 = 1,4
            AN = AN1
            AN1 = -AN2
            AN2 = AN
            AB1 = AN1*RAR(3) - AN2*RAR(1)
            AB2 = AN1*RAR(4) - AN2*RAR(2)
            IJ = IJ + 1
            RLVX(IJ) = AB1
            RLVY(IJ) = AB2
            RLVS(IJ) = SQRT(RLVX(IJ)*RLVX(IJ)+RLVY(IJ)*RLVY(IJ))
            GINT(IJ,1) = -INT(AN2)
            GINT(IJ,2) = INT(AN1)
            NQQ = IJ
            IF ( II.GT.0 ) EXIT
         END DO
         II = 0
      END DO
      IF ( N1.LT.20 ) GOTO 100
      CALL REZORD(RLVX,RLVY,RLVS,NQQ,GINT)
      J = 0
      IRING = 0
      DO I = 1,NQQ - 1
         IRING = IRING + 1
         IF ( ABS(RLVS(I+1)-RLVS(I)).GT.TOL ) THEN
            J = J + 1
            NSK(J) = IRING
            IRING = 0
         END IF
      END DO
      DO I = 1,NQQ
         IF ( IP.GT.3 ) WRITE (NOUT1,99001) I,GINT(I,1),GINT(I,2),
     &                         RLVX(I),RLVY(I),RLVS(I)
      END DO
C
      RETURN
99001 FORMAT (1x,'rezp gv',i4,2x,i3,2x,i3,2x,3F12.5)
      END
C*==rhoset.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RHOSET(RHO,POL)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 POL(3)
      COMPLEX*16 RHO(2,2)
C
C*** End of declarations rewritten by SPAG
C
      RHO(1,1) = 0.5D0*DCMPLX(1.D0+POL(3),0.D0)
      RHO(1,2) = 0.5D0*DCMPLX(POL(1),(-POL(2)))
      RHO(2,1) = DCONJG(RHO(1,2))
      RHO(2,2) = 0.5D0*DCMPLX(1.D0-POL(3),0.D0)
      END
C*==rmesh.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RMESH(RAD,RN,H,NMESH,IAN,IT)
C     /****************************************************************/
C     # purpose      : calculate radial mesh on mesh                   *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H
      INTEGER IAN,IT,NMESH
      REAL*8 RAD(NMESH),RN(NATLM,LAYSM)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NMESH
         RAD(I) = RN(IAN,IT)*EXP(H*DBLE(I-NMESH))
      END DO
C
      END
C*==rmeshda.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RMESHDA(RAD,RN,H,NMESH,ADD,IAN,IT)
C     /****************************************************************/
C     # purpose      : calculate radial mesh on mesh + add             *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ADD,IAN,IT,NMESH
      REAL*8 H
      REAL*8 RAD(NMESH+ADD),RN(NATLM,LAYSM)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NMESH + ADD
         RAD(I) = RN(IAN,IT)*EXP(H*(I-NMESH))
      END DO
C
      END
C*==rmesh1.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RMESH1(RAD1,RN,H,NMESH1,ATOM,LAY)
C     /****************************************************************/
C     # purpose      : calculate radial mesh on mesh+5                 *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,RSTEPP5
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,LAY,NMESH1
      REAL*8 H
      REAL*8 RAD1(RSTEPP5),RN(NATLM,LAYSM)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NMESH1
         RAD1(I) = RN(ATOM,LAY)*EXP(H*DBLE(6-I))
      END DO
C
      END
C*==simps.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SIMPS(H,F,QMT)
C     /****************************************************************/
C     # purpose      : integration over mesh                           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:RSTEP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H
      COMPLEX*16 QMT
      COMPLEX*16 F(RSTEP)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      QMT = F(1)
      DO I = 1,(RSTEP-2),2
         QMT = QMT + 4.D0*F(I+1) + 2.D0*F(I+2)
      END DO
      QMT = (H/3.D0)*(QMT-F(RSTEP))
C
      END
C*==sdouble.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SDOUBLE(QPP1,QPM1,QMP1,QMM1,QPP2,QPM2,QMP2,QMM2,C1ONE,
     &                   C1,C2,C3,C4,DIMVAR)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 C1(DIMVAR,DIMVAR),C1ONE(DIMVAR,DIMVAR),
     &           C2(DIMVAR,DIMVAR),C3(DIMVAR,DIMVAR),C4(DIMVAR,DIMVAR),
     &           QMM1(DIMVAR,DIMVAR),QMM2(DIMVAR,DIMVAR),
     &           QMP1(DIMVAR,DIMVAR),QMP2(DIMVAR,DIMVAR),
     &           QPM1(DIMVAR,DIMVAR),QPM2(DIMVAR,DIMVAR),
     &           QPP1(DIMVAR,DIMVAR),QPP2(DIMVAR,DIMVAR)
C
C*** End of declarations rewritten by SPAG
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C mult      sub       inv       add       copy                *
C**************************************************************
C
C     calculation of qpp2*(inverse(1-qmp1*qpm2)) in c1
C
      CALL MULT(QMP1,QPM2,C3,DIMVAR)
      CALL SUB(C1ONE,C3,C4,DIMVAR)
      CALL INV(C4,C3,DIMVAR)
      CALL MULT(QPP2,C3,C1,DIMVAR)
C
C     calculation of qmm1*(inverse(1-qpm2*qmp1)) in c2
C
      CALL MULT(QPM2,QMP1,C3,DIMVAR)
      CALL SUB(C1ONE,C3,C4,DIMVAR)
      CALL INV(C4,C3,DIMVAR)
      CALL MULT(QMM1,C3,C2,DIMVAR)
C
C     calculation of qpm
C
      CALL MULT(QPM2,QPP1,C3,DIMVAR)
      CALL MULT(C2,C3,C4,DIMVAR)
      CALL ADD(QPM1,C4,C3,DIMVAR)
      CALL COPY(C3,QPM1,DIMVAR)
C
C     calculation of qmp
C
      CALL MULT(QMP1,QMM2,C3,DIMVAR)
      CALL MULT(C1,C3,C4,DIMVAR)
      CALL ADD(QMP2,C4,C3,DIMVAR)
      CALL COPY(C3,QMP1,DIMVAR)
C
C     calculation of qpp
C
      CALL MULT(C1,QPP1,C3,DIMVAR)
      CALL COPY(C3,QPP1,DIMVAR)
C
C     calculation of qmm
C
      CALL MULT(C2,QMM2,C3,DIMVAR)
      CALL COPY(C3,QMM1,DIMVAR)
      END
C*==sgset.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SGSET(SG,QPM,NREL,PQ,GANZ,QD)
      USE MOD_SPEC,ONLY:LG,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,NREL,QD
      INTEGER PQ(LG,2)
      COMPLEX*16 QPM(QD,QD),SG(LG,2,2)
C
C Local variables
C
      INTEGER G,I,J
C
C*** End of declarations rewritten by SPAG
C
      DO G = 1,LG
         DO I = 1,2
            DO J = 1,2
               SG(G,I,J) = CZERO
            END DO
         END DO
      END DO
C
      IF ( NREL.EQ.1 ) THEN
         DO G = 1,LG
            IF ( G.LE.GANZ ) THEN
               SG(G,1,1) = QPM(PQ(G,1),PQ(1,1))
               SG(G,1,2) = QPM(PQ(G,1),PQ(1,2))
               SG(G,2,1) = QPM(PQ(G,2),PQ(1,1))
               SG(G,2,2) = QPM(PQ(G,2),PQ(1,2))
            END IF
         END DO
      END IF
      IF ( NREL.EQ.0 ) THEN
         DO G = 1,LG
            IF ( G.LE.GANZ ) THEN
               SG(G,1,1) = QPM(PQ(G,1),PQ(1,1))
               SG(G,1,2) = CZERO
               SG(G,2,1) = CZERO
               SG(G,2,2) = QPM(PQ(G,1),PQ(1,1))
            END IF
         END DO
      END IF
      END
C*==splout.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SPLOUT(INTG,INTGPM,POLG,TYP,GANZ,ISTR,GINT,N1,N2,EVAC,
     &                  NREL,NTEST,THETA)
      USE MOD_SPEC,ONLY:NPM,LG,HARTRE,PI
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT,IFILSPECOU1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SPLEED')
C
C Dummy arguments
C
      REAL*8 EVAC,THETA
      INTEGER GANZ,N1,N2,NREL,NTEST,TYP
      INTEGER GINT(NPM,2),ISTR(2)
      COMPLEX*16 INTG(LG),INTGPM(LG),POLG(LG,3)
C
C Local variables
C
      REAL*8 ASYMM(:),ENERGY,INTPM(:),P(:,:),RINTG(:)
      INTEGER G,I
C
C*** End of declarations rewritten by SPAG
C
C Dummy arguments
C
C
C Local variables
C
C
      ALLOCATABLE P,RINTG,ASYMM,INTPM
      ALLOCATE (P(3,LG),RINTG(LG),ASYMM(LG),INTPM(LG))
C
C*** End of declarations rewritten by SPAG
C
      DO G = 1,LG
         RINTG(G) = 0.0D0
         INTPM(G) = 0.0D0
         ASYMM(G) = 0.0D0
         DO I = 1,3
            P(I,G) = 0.0D0
         END DO
      END DO
C
      ENERGY = EVAC*HARTRE
C
      IF ( N1.EQ.1 .AND. N2.EQ.1 .AND. NTEST.EQ.1 ) THEN
         IF ( TYP.EQ.0 ) THEN
            WRITE (NOUT1,99001)
         ELSE IF ( TYP.EQ.1 ) THEN
            WRITE (NOUT1,99002)
         ELSE IF ( TYP.EQ.2 ) THEN
            WRITE (NOUT1,99003)
         END IF
      END IF
C
      IF ( TYP.EQ.0 ) THEN
         DO G = 1,LG
            IF ( G.LE.GANZ ) THEN
               RINTG(G) = DBLE(INTG(G))
               IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &              WRITE (*,99005) ENERGY,RINTG(G)
               IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &              WRITE (NOUT1,99005) ENERGY,RINTG(G)
            END IF
         END DO
      END IF
C
      IF ( TYP.EQ.1 .OR. TYP.EQ.2 ) THEN
         DO G = 1,LG
            IF ( G.LE.GANZ ) THEN
               RINTG(G) = DBLE(INTG(G))
               INTPM(G) = DBLE(INTGPM(G))
               P(1,G) = DBLE(POLG(G,1))
               P(2,G) = DBLE(POLG(G,2))
               P(3,G) = DBLE(POLG(G,3))
C               IF ( P(1,G)**2+P(2,G)**2+P(3,G)**2.NE.0.D0 ) THEN
               IF ( P(1,G)**2+P(2,G)**2+P(3,G)**2.GT.1.0D-16 ) THEN
                  ASYMM(G) = (RINTG(G)-INTPM(G))/(RINTG(G)+INTPM(G))
                  ASYMM(G) = ASYMM(G)
     &                       /SQRT(P(1,G)**2+P(2,G)**2+P(3,G)**2)
               ELSE
                  ASYMM(G) = 0.D0
               END IF
            END IF
         END DO
         IF ( NREL.EQ.1 ) THEN
            DO G = 1,LG
               IF ( G.LE.GANZ ) THEN
                  IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &                 WRITE (*,99004) ENERGY,RINTG(G),ASYMM(G),P(1,G),
     &                                 P(2,G),P(3,G)
                  IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &                 WRITE (IFILSPECOU1,99004) ENERGY,THETA*180.0/PI,
     &                        RINTG(G),ASYMM(G),P(1,G),P(2,G),P(3,G)
C     &         WRITE (333,99004) ENERGY,THETA*180.0/PI,RINTG(G),
C     &                        ASYMM(G),P(1,G),P(2,G),P(3,G)
                  IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &                 WRITE (NOUT1,99004) ENERGY,RINTG(G),ASYMM(G),
     &                        P(1,G),P(2,G),P(3,G)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                  IF ( WRBUILDBOT .AND. G.EQ.1 )
     &                 WRITE (IFILBUILDBOT,99006)
     &                 ROUTINE(1:LEN_TRIM(ROUTINE)),10,ENERGY,THETA,
     &                 RINTG(1:10)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               END IF
            END DO
         ELSE IF ( NREL.EQ.0 ) THEN
            DO G = 1,LG
               IF ( G.LE.GANZ ) THEN
                  IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &                 WRITE (*,99005) ENERGY,RINTG(G)
                  IF ( GINT(G,1).EQ.ISTR(1) .AND. GINT(G,2).EQ.ISTR(2) )
     &                 WRITE (NOUT1,99005) ENERGY,RINTG(G)
               END IF
            END DO
         END IF
      END IF
C
      RETURN
99001 FORMAT (' i(e)-diagram'/)
99002 FORMAT (' rotation diagram'/)
99003 FORMAT (' scattering-angle diagram,'/)
99004 FORMAT (e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,
     &        e12.5)
99005 FORMAT (d12.5,2x,d12.5)
99006 FORMAT ('# BUILDBOT: ',A,': (RINTG(G))',',G=1,LG)  for LG=',I5,4x,
     &        ' ENERGY = ',F10.6,4x,'THETA =',F10.6,/,(1PE22.14))
C
      END
C*==streu.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE STREU(TE,TO,DP,IREL,MAXLP1,NATL,IT,XED,XOD)
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQNAT,CZERO,CONE
      USE MOD_SPEC_POSI,ONLY:PHLME,PTKE1,PTKE2,PHLMO,PTKO1,PTKO2
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IREL,IT,MAXLP1,XED,XOD
      COMPLEX*16 DP(ML,2,NATLM,LAYSM),TE(MLQNAT),TO(MLQNAT)
      INTEGER NATL(LAYSM)
C
C Local variables
C
      COMPLEX*16 CTWO,CTWOI,DPH(ML,2,NATLM,LAYSM)
      INTEGER I,IAN,KM,LM,LM1,LM2,NAT,S,XEDH,XODH
C
C*** End of declarations rewritten by SPAG
C
      CTWO = DCMPLX(2.D0,0.D0)
      CTWOI = DCMPLX(0.D0,2.D0)
      DO I = 1,MLQNAT
         TO(I) = CZERO
         TE(I) = CZERO
      END DO
      IF ( IREL.EQ.1 ) THEN
         NAT = NATL(IT)
         DO IAN = 1,NAT
            DO LM = 1,MAXLP1
               DPH(LM,1,IAN,IT) = (CDEXP(CTWOI*DP(LM,1,IAN,IT))-CONE)
     &                            /CTWO
            END DO
         END DO
         NAT = NATL(IT)
         LM1 = 1
         LM2 = 1
         DO IAN = 1,NAT
            XEDH = XED/NAT
            XODH = XOD/NAT
            DO LM = 1,XEDH
               TE(LM1) = DPH(PHLME(LM1),1,IAN,IT)
               LM1 = LM1 + 1
            END DO
            DO LM = 1,XODH
               TO(LM2) = DPH(PHLMO(LM2),1,IAN,IT)
               LM2 = LM2 + 1
            END DO
         END DO
      ELSE IF ( IREL.EQ.2 ) THEN
         NAT = NATL(IT)
         DO IAN = 1,NAT
            DO KM = 1,MAXLP1
               DO S = 1,2
                  IF ( KM.EQ.1 .AND. S.EQ.2 ) THEN
                     DPH(KM,S,IAN,IT) = CZERO
                  ELSE
                     DPH(KM,S,IAN,IT) = (CDEXP(CTWOI*DP(KM,S,IAN,IT))-
     &                                  CONE)/CTWO
                  END IF
               END DO
            END DO
         END DO
         NAT = NATL(IT)
         LM1 = 1
         LM2 = 1
         DO IAN = 1,NAT
            XEDH = XED/NAT
            XODH = XOD/NAT
            DO KM = 1,XEDH
               TE(LM1) = DPH(PTKE1(KM),PTKE2(KM),IAN,IT)
               LM1 = LM1 + 1
            END DO
            DO KM = 1,XODH
               TO(LM2) = DPH(PTKO1(KM),PTKO2(KM),IAN,IT)
               LM2 = LM2 + 1
            END DO
         END DO
      END IF
      END
C*==suche.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      INTEGER FUNCTION SUCHE(L,M,LMW,ILM)
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ILM,L,M
      INTEGER LMW(ILM,2)
C
C Local variables
C
      INTEGER LM,OK
C
C*** End of declarations rewritten by SPAG
C
      OK = 0
      LM = 0
 100  CONTINUE
      LM = LM + 1
      IF ( LM.LE.ILM ) THEN
         IF ( L.EQ.LMW(LM,1) .AND. M.EQ.LMW(LM,2) ) THEN
            SUCHE = LM
            OK = 1
         END IF
         IF ( OK.EQ.0 ) GOTO 100
      END IF
      IF ( OK.EQ.0 ) THEN
         WRITE (NOUT1,99001) L,M
         STOP
      END IF
      RETURN
99001 FORMAT (1x,'position of ',2I5,' not found')
      END
C*==xmat.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XMAT(X,NX,MX,IT,IOE,NOE,NATL,DLM1,CELM,NCLM)
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLZP,CONE,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LMAXM,L2MAXM,NFM,NDLMM
      PARAMETER (LMAXM=ML-1,L2MAXM=2*LMAXM,NFM=NATLM*NATLM-NATLM+1,
     &           NDLMM=(L2MAXM+1)*(L2MAXM+1))
C
C Dummy arguments
C
      INTEGER IOE,IT,MX,NCLM,NOE,NX
      REAL*8 CELM(MLZP)
      COMPLEX*16 DLM1(NDLMM,NFM),X(NX,MX)
      INTEGER NATL(LAYSM)
C
C Local variables
C
      INTEGER I,I1,IFLM,IROW,IS,J,J1,JS,K,L,L2,LMAX,M,N,NAT
C
C*** End of declarations rewritten by SPAG
C
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C block                                                       *
C**************************************************************
C
C
      LMAX = ML - 1
      IS = 0
      JS = 0
      IFLM = 1
      CALL BLOCK(X,NX,MX,IS,JS,IFLM,IOE,DLM1,CELM,NCLM)
      NAT = NATL(IT)
      IF ( NAT.GT.1 ) THEN
         DO K = 2,NAT
            IS = IS + NOE
            JS = JS + NOE
            DO J = 1,NOE
               J1 = JS + J
               DO I = 1,NOE
                  I1 = IS + I
                  X(I1,J1) = X(I,J)
               END DO
            END DO
         END DO
         DO I = 1,NAT - 1
            IS = (I-1)*NOE
            DO J = I + 1,NAT
               JS = (J-1)*NOE
               IFLM = IFLM + 1
               CALL BLOCK(X,NX,MX,IS,JS,IFLM,IOE,DLM1,CELM,NCLM)
               IFLM = IFLM + 1
               CALL BLOCK(X,NX,MX,JS,IS,IFLM,IOE,DLM1,CELM,NCLM)
            END DO
         END DO
      END IF
      IROW = 0
      DO K = 1,NAT
         DO L = 0,LMAX
            L2 = L + IOE
            IF ( L2.NE.0 ) THEN
               DO M = 1,L2
                  IROW = IROW + 1
                  DO N = 1,MX
                     X(IROW,N) = X(IROW,N)*CIMAG
                  END DO
                  X(IROW,IROW) = X(IROW,IROW) - CONE
               END DO
            END IF
         END DO
      END DO
      END
C*==xpen.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XPEN(XP,T,GI,XD)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER XD
      COMPLEX*16 GI(XD,XD),T(XD),XP(XD,XD)
C
C Local variables
C
      INTEGER LM,LMS
C
C*** End of declarations rewritten by SPAG
C
      DO LM = 1,XD
         DO LMS = 1,XD
            XP(LM,LMS) = T(LM)*GI(LM,LMS)
         END DO
      END DO
      END
C*==xpenr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XPENR(XP,T,GI,XD1,XD2)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER XD1,XD2
      COMPLEX*16 GI(XD1,XD2),T(XD1),XP(XD1,XD2)
C
C Local variables
C
      INTEGER LM,LMS
C
C*** End of declarations rewritten by SPAG
C
      DO LM = 1,XD1
         DO LMS = 1,XD2
            XP(LM,LMS) = T(LM)*GI(LM,LMS)
         END DO
      END DO
      END
C*==xrelp.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XRELP(X,TK,CLX,IKM,GI1,PX1,ILM1,GI2,PX2,ILM2,I,J,NAT)
C
      USE MOD_SPEC,ONLY:MLQNAT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IKM,ILM1,ILM2,J,NAT
      COMPLEX*16 CLX(MLQNAT,MLQNAT,2),GI1(ILM1,ILM1),GI2(ILM2,ILM2),
     &           TK(IKM),X(IKM,IKM)
      INTEGER PX1(MLQNAT,2),PX2(MLQNAT,2)
C
C Local variables
C
      INTEGER I1,I2,IE1,IE2,IND1,IND2,IND3,IND4,INDEXVAR,KM,KMS,NM,NMS
      COMPLEX*16 SUM1,SUM2
C
C*** End of declarations rewritten by SPAG
C
C
      INDEXVAR = IKM/NAT
      IE1 = (I-1)*INDEXVAR
      IE2 = (J-1)*INDEXVAR
      IND1 = (I-1)*(ILM1/NAT)
      IND2 = (J-1)*(ILM1/NAT)
      IND3 = (I-1)*(ILM2/NAT)
      IND4 = (J-1)*(ILM2/NAT)
      DO KM = 1,INDEXVAR
         DO KMS = 1,INDEXVAR
            NM = KM + IE1
            NMS = KMS + IE2
            IF ( PX1(NM,1).NE.0 .AND. PX1(NMS,1).NE.0 ) THEN
               I1 = PX1(NM,1) + IND1
               I2 = PX1(NMS,1) + IND2
               SUM1 = CLX(NM,NMS,2)*GI1(I1,I2)
            ELSE
               SUM1 = CZERO
            END IF
            IF ( PX2(NM,2).NE.0 .AND. PX2(NMS,2).NE.0 ) THEN
               I1 = PX2(NM,2) + IND3
               I2 = PX2(NMS,2) + IND4
               SUM2 = CLX(NM,NMS,1)*GI2(I1,I2)
            ELSE
               SUM2 = CZERO
            END IF
            X(NM,NMS) = TK(NM)*(SUM1+SUM2)
         END DO
      END DO
C
      END
C*==xrelpr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XRELPR(X,TK,CLX,IKM,GI1,PX1,ILM1,GI2,PX2,ILM2,I,J,NAT)
C
      USE MOD_SPEC,ONLY:MLQNAT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IKM,ILM1,ILM2,J,NAT
      COMPLEX*16 CLX(MLQNAT,MLQNAT,2),GI1(ILM2,ILM1),GI2(ILM1,ILM2),
     &           TK(IKM),X(IKM,IKM)
      INTEGER PX1(MLQNAT,2),PX2(MLQNAT,2)
C
C Local variables
C
      INTEGER I1,I2,IE1,IE2,IND1,IND2,IND3,IND4,INDEXVAR,KM,KMS,NM,NMS
      COMPLEX*16 SUM1,SUM2
C
C*** End of declarations rewritten by SPAG
C
      INDEXVAR = IKM/NAT
      IE1 = (I-1)*INDEXVAR
      IE2 = (J-1)*INDEXVAR
      IND1 = (I-1)*(ILM2/NAT)
      IND2 = (J-1)*(ILM1/NAT)
      IND3 = (I-1)*(ILM1/NAT)
      IND4 = (J-1)*(ILM2/NAT)
      DO KM = 1,INDEXVAR
         DO KMS = 1,INDEXVAR
            NM = KM + IE1
            NMS = KMS + IE2
            IF ( PX2(NM,2).NE.0 .AND. PX1(NMS,2).NE.0 ) THEN
               I1 = PX2(NM,2) + IND1
               I2 = PX1(NMS,2) + IND2
               SUM1 = CLX(NM,NMS,1)*GI1(I1,I2)
            ELSE
               SUM1 = CZERO
            END IF
            IF ( PX1(NM,1).NE.0 .AND. PX2(NMS,1).NE.0 ) THEN
               I1 = PX1(NM,1) + IND3
               I2 = PX2(NMS,1) + IND4
               SUM2 = CLX(NM,NMS,2)*GI2(I1,I2)
            ELSE
               SUM2 = CZERO
            END IF
            X(NM,NMS) = TK(NM)*(SUM1+SUM2)
         END DO
      END DO
      END
C*==xrelfe.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XRELFE(X,CLX,IKM,GI1,PX1,ILM1,GI2,PX2,ILM2,I,J,NAT)
C
      USE MOD_SPEC,ONLY:MLQNAT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IKM,ILM1,ILM2,J,NAT
      COMPLEX*16 CLX(MLQNAT,MLQNAT,2),GI1(ILM1,ILM1),GI2(ILM2,ILM2),
     &           X(IKM,IKM)
      INTEGER PX1(MLQNAT,2),PX2(MLQNAT,2)
C
C Local variables
C
      INTEGER I1,I2,IE1,IE2,IND1,IND2,IND3,IND4,INDEXVAR,KM,KMS,NM,NMS
      COMPLEX*16 SUM1,SUM2
C
C*** End of declarations rewritten by SPAG
C
C
      INDEXVAR = IKM/NAT
      IE1 = (I-1)*INDEXVAR
      IE2 = (J-1)*INDEXVAR
      IND1 = (I-1)*(ILM1/NAT)
      IND2 = (J-1)*(ILM1/NAT)
      IND3 = (I-1)*(ILM2/NAT)
      IND4 = (J-1)*(ILM2/NAT)
      DO KM = 1,INDEXVAR
         DO KMS = 1,INDEXVAR
            NM = KM + IE1
            NMS = KMS + IE2
            IF ( PX1(NM,1).NE.0 .AND. PX1(NMS,1).NE.0 ) THEN
               I1 = PX1(NM,1) + IND1
               I2 = PX1(NMS,1) + IND2
               SUM1 = CLX(NM,NMS,2)*GI1(I1,I2)
            ELSE
               SUM1 = CZERO
            END IF
            IF ( PX2(NM,2).NE.0 .AND. PX2(NMS,2).NE.0 ) THEN
               I1 = PX2(NM,2) + IND3
               I2 = PX2(NMS,2) + IND4
               SUM2 = CLX(NM,NMS,1)*GI2(I1,I2)
            ELSE
               SUM2 = CZERO
            END IF
            X(NM,NMS) = (SUM1+SUM2)
         END DO
      END DO
C
      END
C*==xrelfer.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE XRELFER(X,CLX,IKM,GI1,PX1,ILM1,GI2,PX2,ILM2,I,J,NAT)
C
      USE MOD_SPEC,ONLY:MLQNAT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IKM,ILM1,ILM2,J,NAT
      COMPLEX*16 CLX(MLQNAT,MLQNAT,2),GI1(ILM2,ILM1),GI2(ILM1,ILM2),
     &           X(IKM,IKM)
      INTEGER PX1(MLQNAT,2),PX2(MLQNAT,2)
C
C Local variables
C
      INTEGER I1,I2,IE1,IE2,IND1,IND2,IND3,IND4,INDEXVAR,KM,KMS,NM,NMS
      COMPLEX*16 SUM1,SUM2
C
C*** End of declarations rewritten by SPAG
C
      INDEXVAR = IKM/NAT
      IE1 = (I-1)*INDEXVAR
      IE2 = (J-1)*INDEXVAR
      IND1 = (I-1)*(ILM2/NAT)
      IND2 = (J-1)*(ILM1/NAT)
      IND3 = (I-1)*(ILM1/NAT)
      IND4 = (J-1)*(ILM2/NAT)
      DO KM = 1,INDEXVAR
         DO KMS = 1,INDEXVAR
            NM = KM + IE1
            NMS = KMS + IE2
            IF ( PX2(NM,2).NE.0 .AND. PX1(NMS,2).NE.0 ) THEN
               I1 = PX2(NM,2) + IND1
               I2 = PX1(NMS,2) + IND2
               SUM1 = CLX(NM,NMS,1)*GI1(I1,I2)
            ELSE
               SUM1 = CZERO
            END IF
            IF ( PX1(NM,1).NE.0 .AND. PX2(NMS,1).NE.0 ) THEN
               I1 = PX1(NM,1) + IND3
               I2 = PX2(NMS,1) + IND4
               SUM2 = CLX(NM,NMS,2)*GI2(I1,I2)
            ELSE
               SUM2 = CZERO
            END IF
            X(NM,NMS) = (SUM1+SUM2)
         END DO
      END DO
      END
C*==trans.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE TRANS(XALM,NALM,NEVEN,NODD,NDIM,NATIN)
      USE MOD_SPEC,ONLY:MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NALM,NATIN,NDIM,NEVEN,NODD
      COMPLEX*16 XALM(MLQNAT,MLQNAT)
C
C Local variables
C
      INTEGER I,II,ISH,J,JJ,JSH,K,KD,KE,KO,NODDNAT
      COMPLEX*16 XALM1(MLQNAT,MLQNAT)
C
C*** End of declarations rewritten by SPAG
C
C
      NODDNAT = NODD*NATIN
      DO JJ = 1,NATIN
         JSH = (JJ-1)*NDIM
         DO II = 1,NATIN
            ISH = (II-1)*NDIM
            DO I = 1,NEVEN
               DO J = 1,NEVEN
                  XALM1(J+JSH+NODD,I+ISH+NODD) = XALM(J+JSH,I+ISH)
               END DO
            END DO
            DO I = 1,NEVEN
               DO J = NEVEN + 1,NDIM
                  XALM1(J+JSH-NEVEN,I+ISH+NODD) = XALM(J+JSH,I+ISH)
               END DO
            END DO
            DO I = NEVEN + 1,NDIM
               DO J = 1,NEVEN
                  XALM1(J+JSH+NODD,I+ISH-NEVEN) = XALM(J+JSH,I+ISH)
               END DO
            END DO
            DO I = NEVEN + 1,NDIM
               DO J = NEVEN + 1,NDIM
                  XALM1(J+JSH-NEVEN,I+ISH-NEVEN) = XALM(J+JSH,I+ISH)
               END DO
            END DO
         END DO
      END DO
      DO J = 1,NALM
         DO I = 1,NODD
            XALM(I,J) = XALM1(I,J)
         END DO
         DO I = 1,NEVEN
            II = I - NEVEN + NALM
            XALM(II,J) = XALM1(II,J)
         END DO
      END DO
      DO K = 2,NATIN
         KO = (K-1)*NODD
         KE = (K-2)*NEVEN + NODDNAT
         KD = (K-2)*NDIM
         DO J = 1,NALM
            DO I = 1,NODD
               II = I + KD
               XALM(I+KO,J) = XALM1(II+NDIM,J)
            END DO
            DO I = 1,NEVEN
               II = I + KD
               XALM(I+KE,J) = XALM1(II+NODD,J)
            END DO
         END DO
      END DO
      DO J = 1,NODD
         DO I = 1,NALM
            XALM1(I,J) = XALM(I,J)
         END DO
      END DO
      DO J = 1,NEVEN
         JJ = J - NEVEN + NALM
         DO I = 1,NALM
            XALM1(I,JJ) = XALM(I,JJ)
         END DO
      END DO
      DO K = 2,NATIN
         KO = (K-1)*NODD
         KE = (K-2)*NEVEN + NODDNAT
         KD = (K-2)*NDIM
         DO J = 1,NODD
            JJ = J + KD
            DO I = 1,NALM
               XALM1(I,J+KO) = XALM(I,JJ+NDIM)
            END DO
         END DO
         DO J = 1,NEVEN
            JJ = J + KD
            DO I = 1,NALM
               XALM1(I,J+KE) = XALM(I,JJ+NODD)
            END DO
         END DO
      END DO
      DO I = 1,NALM
         DO J = 1,NALM
            XALM(J,I) = XALM1(J,I)
         END DO
      END DO
      END
C*==sphrkc.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SPHRKC(YLM,CT,ST,CF,LMAX)
      USE MOD_SPEC,ONLY:ML,PI,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLA,MLB
      PARAMETER (MLA=2*ML-1,MLB=MLA*MLA)
C
C Dummy arguments
C
      COMPLEX*16 CF,CT,ST
      INTEGER LMAX
      COMPLEX*16 YLM(MLB)
C
C Local variables
C
      REAL*8 A,ASG,B,CL,CM,FAC1(:),FAC2(:),FAC3(:)
      INTEGER L,LLL,LM,LM2,LM3,LN,LO,LP,LQ,M
      COMPLEX*16 SA,SF
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE FAC1,FAC2,FAC3
      ALLOCATE (FAC1(MLA),FAC2(MLB),FAC3(MLA))
C
C*** End of declarations rewritten by SPAG
C
C     sphrkc,given ct=cos(theta),st=sin(theta),and
C     cf=exp(i*fi),calculates all the ylm(theta,fi)
C     up to l=lmax.subscripts are ordered thus:lm=(00),
C     (1-1),(10),(11),(2-2),....
C     dimensioned for lmax=4
C
      LM = 0
      CL = 0.D0
      A = 1.D0
      B = 1.D0
      ASG = 1.D0
      LLL = LMAX + 1
C
C     multiplicative factors required
C
      DO L = 1,LLL
         FAC1(L) = ASG*SQRT((2.D0*CL+1.D0)*A/(4.D0*PI*B*B))
         FAC3(L) = SQRT(2.D0*CL)
         CM = -CL
         LN = L + L - 1
         DO M = 1,LN
            LO = LM + M
            FAC2(LO) = SQRT((CL+1.D0+CM)*(CL+1.D0-CM)/((2.D0*CL+3.D0)*(
     &                 2.D0*CL+1.D0)))
            CM = CM + 1.D0
         END DO
         CL = CL + 1.D0
         A = A*2.D0*CL*(2.D0*CL-1.D0)/4.D0
         B = B*CL
         ASG = -ASG
         LM = LM + LN
      END DO
C
C      first all the ylm for m=+-l and m=+-(l-1) are
C      are calculated by explicit formulae
C
      LM = 1
      CL = 1.D0
      ASG = -1.D0
      SF = CF
      SA = CONE
      YLM(1) = DCMPLX(FAC1(1),0.D0)
      DO L = 1,LMAX
         LN = LM + L + L + 1
         YLM(LN) = FAC1(L+1)*SA*SF*ST
         YLM(LM+1) = ASG*FAC1(L+1)*SA*ST/SF
         YLM(LN-1) = -FAC3(L+1)*FAC1(L+1)*SA*SF*CT/CF
         YLM(LM+2) = ASG*FAC3(L+1)*FAC1(L+1)*SA*CT*CF/SF
         SA = ST*SA
         SF = SF*CF
         CL = CL + 1.D0
         ASG = -ASG
         LM = LN
      END DO
C
C       using ylm and yl(m-1) in a recurrence relation
C       yl(m+1) is calculated
C
      LM = 1
      LLL = LMAX - 1
      DO L = 1,LLL
         LN = L + L - 1
         LM2 = LM + LN + 4
         LM3 = LM - LN
         DO M = 1,LN
            LO = LM2 + M
            LP = LM3 + M
            LQ = LM + M + 1
            YLM(LO) = -(FAC2(LP)*YLM(LP)-CT*YLM(LQ))/FAC2(LQ)
         END DO
         LM = LM + L + L + 1
      END DO
      END
C*==ckmulm.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE CKMULM(IKM,TSKEL,TSKOL,TSKEOL,TSKOEL,TSKEH,TSKOH,
     &                  TSKEOH,TSKOEH,ISTATE,IT,IAN,GAMMA1)
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MLQ,MQD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IAN,IKM,ISTATE,IT
      COMPLEX*16 GAMMA1(LAYSM,NATLM,3,MQD,MQD),TSKEH(IKM,IKM,LAYSM),
     &           TSKEL(IKM,IKM,LAYSM),TSKEOH(IKM,IKM,LAYSM),
     &           TSKEOL(IKM,IKM,LAYSM),TSKOEH(IKM,IKM,LAYSM),
     &           TSKOEL(IKM,IKM,LAYSM),TSKOH(IKM,IKM,LAYSM),
     &           TSKOL(IKM,IKM,LAYSM)
C
C Local variables
C
      COMPLEX*16 GAMMAE(:,:),GAMMAEO(:,:),GAMMAO(:,:),GAMMAOE(:,:)
      INTEGER I,IANF,IEND,J,K,L
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE GAMMAE,GAMMAO,GAMMAEO,GAMMAOE
      ALLOCATE (GAMMAE(MLQ,MLQ),GAMMAO(MLQ,MLQ),GAMMAEO(MLQ,MLQ))
      ALLOCATE (GAMMAOE(MLQ,MLQ))
C
C*** End of declarations rewritten by SPAG
C
C*****changes parameters of t-matrix from mu,kappa to l,m ***************
C
      CALL GAMEO(GAMMAE,GAMMAO,GAMMAEO,GAMMAOE,ISTATE,IT,IAN,GAMMA1)
      IANF = 1 + MLQ*(IAN-1)
      IEND = MLQ*IAN
      K = 0
      DO I = IANF,IEND
         K = K + 1
         L = 0
         DO J = IANF,IEND
            L = L + 1
            IF ( ISTATE.EQ.1 ) THEN
               TSKEL(I,J,IT) = GAMMAE(K,L)
               TSKOL(I,J,IT) = GAMMAO(K,L)
               TSKEOL(I,J,IT) = GAMMAEO(K,L)
               TSKOEL(I,J,IT) = GAMMAOE(K,L)
            ELSE
               TSKEH(I,J,IT) = GAMMAE(K,L)
               TSKOH(I,J,IT) = GAMMAO(K,L)
               TSKEOH(I,J,IT) = GAMMAEO(K,L)
               TSKOEH(I,J,IT) = GAMMAOE(K,L)
            END IF
         END DO
      END DO
      END
C*==gameo.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE GAMEO(GAMMAE,GAMMAO,GAMMAEO,GAMMAOE,ISTATE,IT,IAN,
     &                 GAMMA1)
C
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQ,MQD,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IAN,ISTATE,IT
      COMPLEX*16 GAMMA1(LAYSM,NATLM,3,MQD,MQD),GAMMAE(MLQ,MLQ),
     &           GAMMAEO(MLQ,MLQ),GAMMAO(MLQ,MLQ),GAMMAOE(MLQ,MLQ)
C
C Local variables
C
      INTEGER I,IP,J,JP,MIT
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MLQ
         DO J = 1,MLQ
            GAMMAE(I,J) = CZERO
            GAMMAO(I,J) = CZERO
            GAMMAEO(I,J) = CZERO
            GAMMAOE(I,J) = CZERO
         END DO
      END DO
C
      MIT = ML*ML + ML
      DO I = 1,MQD
         DO J = 1,MQD
            IP = INT(DBLE(I+1)/2.D0)
            JP = INT(DBLE(J+1)/2.D0)
            IF ( I.LE.(MIT) ) THEN
               IF ( MOD(I,2).NE.0 ) THEN
                  IF ( (MOD(J,2).NE.0 .AND. J.LE.MIT) .OR. 
     &                 (MOD(J,2).EQ.0 .AND. J.GT.MIT) ) THEN
                     GAMMAO(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
                  ELSE
                     GAMMAOE(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
                  END IF
               ELSE IF ( (MOD(J,2).EQ.0 .AND. J.LE.MIT) .OR. 
     &                   (MOD(J,2).NE.0 .AND. J.GT.MIT) ) THEN
                  GAMMAE(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
               ELSE
                  GAMMAEO(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
               END IF
            ELSE IF ( MOD(I,2).NE.0 ) THEN
               IF ( (MOD(J,2).NE.0 .AND. J.GT.MIT) .OR. 
     &              (MOD(J,2).EQ.0 .AND. J.LE.MIT) ) THEN
                  GAMMAE(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
               ELSE
                  GAMMAEO(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
               END IF
            ELSE IF ( (MOD(J,2).EQ.0 .AND. J.GT.MIT) .OR. 
     &                (MOD(J,2).NE.0 .AND. J.LE.MIT) ) THEN
               GAMMAO(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
            ELSE
               GAMMAOE(IP,JP) = GAMMA1(IT,IAN,ISTATE,I,J)
            END IF
         END DO
      END DO
C
      END
C*==cleg.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      REAL*8 FUNCTION CLEG(KAPPA,MU,S)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KAPPA,MU,S
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(MU).GE.ABS(KAPPA) ) THEN
         CLEG = 0.D0
      ELSE IF ( KAPPA.GT.0.D0 ) THEN
         IF ( S.GT.0.D0 ) THEN
C           if (mu.gt.0.d0) go to 190
            CLEG = -SQRT((KAPPA-MU+0.5D0)/(2.D0*KAPPA+1.D0))
         ELSE IF ( S.LT.0.D0 ) THEN
C              if (mu.lt.0.d0) go to 190
            CLEG = SQRT((KAPPA+MU+0.5D0)/(2.D0*KAPPA+1.D0))
         END IF
      ELSE IF ( KAPPA.LT.0.D0 ) THEN
         IF ( S.GT.0.D0 ) THEN
C              if (mu.lt.0.d0) go to 190
            CLEG = SQRT((-KAPPA+MU-0.5D0)/(-2.D0*KAPPA-1.D0))
         ELSE IF ( S.LT.0.D0 ) THEN
C                 if (mu.gt.0.d0) go to 190
            CLEG = SQRT((-KAPPA-MU-0.5D0)/(-2.D0*KAPPA-1.D0))
         END IF
      END IF
C
      END
C*==rind.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE RIND(ML,MQD,KAP,MUD,NRL,NRM,NRS,NR,REL)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ML,MQD
      INTEGER KAP(MQD),MUD(MQD),NR(MQD),NRL(MQD),NRM(MQD),REL(MQD,2)
      REAL*8 NRS(MQD)
C
C Local variables
C
      INTEGER II,KAPI,LMAX,MUI,NRZ,ZSPZ
C
C*** End of declarations rewritten by SPAG
C
C
      LMAX = ML - 1
      II = 0
      DO KAPI = -LMAX - 1,LMAX
         IF ( KAPI.NE.0 ) THEN
            DO MUI = -ABS(2*KAPI) + 1,ABS(2*KAPI) - 1,2
               II = II + 1
               KAP(II) = KAPI
               MUD(II) = MUI
               NRS(II) = -SIGN(1,KAPI)*SIGN(1,MUI)/2.D0
               NRM(II) = INT(MUI/2.D0-NRS(II))
               NRL(II) = INT(KAPI)
               IF ( KAPI.LT.0 ) NRL(II) = -KAPI - 1
               NR(II) = NRL(II)**2 + NRL(II) + NRM(II) + 1
               ZSPZ = INT(-NRS(II)+1.5D0)
               NRZ = NR(II)
               REL(NRZ,ZSPZ) = II
            END DO
         END IF
      END DO
C
C     do i=1,mqd
C        write (*,*) i,kap(i),mud(i),nrl(i),nrm(i),nrs(i)
C     enddo
C
      END
C*==clgordon.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      REAL*8 FUNCTION CLGORDON(J,M,S)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 J,M,S
C
C Local variables
C
      REAL*8 J1,RES
C
C*** End of declarations rewritten by SPAG
C
C      IF ( J.EQ.1.5D0 ) THEN
      IF ( ABS(J-1.5D0).LE.1.0D-16 ) THEN
         J1 = J - 0.5D0
C         IF ( S.EQ.0.5D0 ) THEN
         IF ( ABS(S-0.5D0).LE.1.0D-16 ) THEN
            RES = SQRT((J1+M+0.5D0)/(2*J1+1))
C         ELSE IF ( S.EQ.-0.5D0 ) THEN
         ELSE IF ( ABS(S+0.5D0).LE.1.0D-16 ) THEN
            RES = SQRT((J1-M+0.5D0)/(2*J1+1))
         ELSE
            WRITE (*,*) 'wrong input in clgordon'
            STOP
         END IF
C      ELSE IF ( J.EQ.0.5D0 ) THEN
      ELSE IF ( ABS(J-0.5D0).LE.1.0D-16 ) THEN
         J1 = J + 0.5D0
C         IF ( S.EQ.0.5D0 ) THEN
         IF ( ABS(S-0.5D0).LE.1.0D-16 ) THEN
            RES = -SQRT((J1-M+0.5D0)/(2*J1+1))
C         ELSE IF ( S.EQ.-0.5D0 ) THEN
         ELSE IF ( ABS(S+0.5D0).LE.1.0D-16 ) THEN
            RES = SQRT((J1+M+0.5D0)/(2*J1+1))
         ELSE
            WRITE (*,*) 'wrong input in clgordon'
            STOP
         END IF
      ELSE
         WRITE (*,*) 'wrong input in clgordon'
         STOP
      END IF
C
      CLGORDON = RES
C
      END
C*==bcopy.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE BCOPY(X,Y,DIM1,DIM2)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIM1,DIM2
      COMPLEX*16 X(DIM1,DIM1),Y(DIM2,DIM2)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      IF ( DIM2.GT.DIM1 ) THEN
         WRITE (*,*) 'error in subroutine bcopy. aborting ...'
         STOP
      END IF
      DO I = 1,DIM2
         DO J = 1,DIM2
            Y(I,J) = X(I,J)
         END DO
      END DO
      END
C*==copyts.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE COPYTS(TS,Y,DIMVAR,ISEQ)
C     /****************************************************************/
C     # purpose:       copy t-matrix from tse/tso to work arrays       *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR,ISEQ
      COMPLEX*16 TS(MLQNAT,MLQNAT,LAYSM),Y(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Y(I,J) = TS(I,J,ISEQ)
         END DO
      END DO
      END
C*==copytsr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE COPYTSR(TS,Y,DIM1,DIM2,ISEQ)
C     /****************************************************************/
C     # purpose:       copy t-matrix from tse/tso to work arrays       *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIM1,DIM2,ISEQ
      COMPLEX*16 TS(MLQNAT,MLQNAT,LAYSM),Y(DIM1,DIM2)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIM1
         DO J = 1,DIM2
            Y(I,J) = TS(I,J,ISEQ)
         END DO
      END DO
      END
C*==copyq.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE COPYQ(Q,Y,DIMVAR,ISEQ)
C     /****************************************************************/
C     # purpose:       copy q-matrix to work arrays                    *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,LH
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR,ISEQ
      COMPLEX*16 Q(LH,LH,LAYSM),Y(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Y(I,J) = Q(I,J,ISEQ)
         END DO
      END DO
      END
C*==copybq.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE COPYBQ(X,Q,DIMVAR,ISEQ)
C     /****************************************************************/
C     # purpose:       copy work arrays to q-matrix                    *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,LH
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR,ISEQ
      COMPLEX*16 Q(LH,LH,LAYSM),X(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Q(I,J,ISEQ) = X(I,J)
         END DO
      END DO
      END
C*==tcopy.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE TCOPY(X,XF,DIMVAR,ISEQ)
C     /****************************************************************/
C     # purpose:       copy work arrays to x-matrix                    *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR,ISEQ
      COMPLEX*16 X(DIMVAR,DIMVAR),XF(MLQNAT,MLQNAT,LAYSM)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            XF(I,J,ISEQ) = X(I,J)
         END DO
      END DO
      END
C*==tcopyr.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE TCOPYR(X,XF,DIM1,DIM2,ISEQ)
C     /****************************************************************/
C     # purpose:       copy work arrays to x-matrix                    *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIM1,DIM2,ISEQ
      COMPLEX*16 X(DIM1,DIM2),XF(MLQNAT,MLQNAT,LAYSM)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIM1
         DO J = 1,DIM2
            XF(I,J,ISEQ) = X(I,J)
         END DO
      END DO
      END
C*==kfold.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE KFOLD(K1,K2,BBAS)
C   ********************************************************************
C   *  DK, January 2006, Modified by JM for SPEC proposes              *
C   *                                                                  *
C   *  Checks whether k1 can be translated into k2 by a reciprocal     *
C   *  lattice vector                                                  *
C   *                                                                  *
C   *                T                                                 *
C   * +  uses  :    A  *  B = I                                        *
C
C   *                          3                                       *
C   *                                                                  *
C   *          where A is matrix of column vectors of real space unit  *
C   *                     cell in units of alat                        *
C   *                B is matrix of column vectors of reciprocal space *
C   *                     unit cell in units  2Pi/alat                 *
C   *                                                                  *
C   * + checks :   whether                                             *
C   *                                                                  *
C   *                                                                  *
C   *   [ m       -1                T                                  *
C   *     n    = B   * (k1 - k2) = A * (k1 - k2)  is an integer vector *
C   *     p ]                                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BBAS(2,2),K1(2),K2(2)
C
C Local variables
C
      REAL*8 AR(2,2),DET,KD(2),MNP(2)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DET = 1.0D0/(BBAS(1,1)*BBAS(2,2)-BBAS(1,2)*BBAS(2,1))
      AR(1,1) = BBAS(2,2)
      AR(1,2) = -BBAS(1,2)
      AR(2,1) = -BBAS(2,1)
      AR(2,2) = BBAS(1,1)
      AR = AR*DET
C
C
      KD(1:2) = K1(1:2)
C
      CALL DGEMV('N',2,2,1D0,AR,2,KD,1,0D0,MNP,1)
C
      K2 = 0.0D0
      DO I = 1,2
         K2(I) = K1(I) - NINT(MNP(I))*BBAS(I,1) - NINT(MNP(I))*BBAS(I,2)
      END DO
C
      END
