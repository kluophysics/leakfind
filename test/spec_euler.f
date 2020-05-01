C*==selecteuler.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     euler.f
C
C     subroutines for rotation of wave functions etc:
C         selecteuler
C         rotbfield
C         euler
C         rotmatjb
C     functions
C         factr
C         factorial
C
      SUBROUTINE SELECTEULER(IBLOCH,EB,BB,NPARA,IPARA,MCD,USEEULER,ROT,
     &                       ROTM1,IBXY)
C     /****************************************************************/
C     # purpose       : prepare rotation for core-level calculation    *
C                       and set selection rules                        *
C                                                                      *
C     # parameter     :                                                *
C         useeuler = 0 for valence band calculations                   *
C                     (spleed, band structure, or ups (vbrun))         *
C         useeuler = 1 for core level calculations with in plane       *
C                      magnetisation (xps, xas, xes, or aes)           *
C                      may also be used in valence band calculations   *
C         npara    = 0 switch off selection rules                      *
C                  = 1 switch on  selection rules                      *
C         ipara    = 0 switch off selection rules                      *
C                  = 1 paramagnetic case                               *
C                  = 2 fullpot and ferromagnetic case                  *
C                  = 3 ferromagnetic corelevel case (useeuler=1)       *
C     # note          : useeuler is preset in inputs                   *
C                                                                      *
C     # calls the following subroutines:                               *
C       euler     rotmatjb                                             *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,EPS12
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IBLOCH,IBXY,IPARA,MCD,NPARA,USEEULER
      REAL*8 BB(3),EB(3)
      COMPLEX*16 ROT(MQD,MQD),ROTM1(MQD,MQD)
C
C Local variables
C
      REAL*8 BFIELD,BXY,BZ,EUA,EUB,EUC
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C     initialize matrices
      EUA = 0.D0
      EUB = 0.D0
      EUC = 0.D0
      CALL ROTMATJB(EUA,EUB,EUC,ROT,ROTM1)
C
C     set b-field components
      DO I = 1,3
         BB(I) = EB(I)
      END DO
      BFIELD = SQRT(EB(1)**2+EB(2)**2+EB(3)**2)
      BXY = SQRT(EB(1)**2+EB(2)**2)
      BZ = EB(3)
C
C     pre-set selection rules
      IPARA = NPARA
      IF ( BFIELD.GT.EPS12 ) THEN
C         ferromagnetic calculation
         IF ( IPARA.NE.0 ) IPARA = 2
      ELSE
C         paramagnetic calculation
         USEEULER = 0
         IF ( IPARA.NE.0 ) IPARA = 1
         IF ( MCD.NE.0 ) THEN
            MCD = 0
            WRITE (*,99001)
            WRITE (NOUT1,99001)
         END IF
      END IF
C
      IF ( BXY.GT.EPS12 ) THEN
         IBXY = 1
      ELSE IF ( MCD.EQ.2 ) THEN
         MCD = 0
         WRITE (*,99002)
         WRITE (NOUT1,99002)
      END IF
C
C     final settig of calculation type
      SELECT CASE (IBLOCH)
      CASE (0,4,7,8,9)
C         valence band calculations
         IF ( USEEULER.NE.0 ) USEEULER = 1
      CASE DEFAULT
C         corelevel-calculations
         IF ( BXY.GT.0. ) USEEULER = 1
C         IF ( USEEULER.NE.0 .AND. BZ.NE.0. ) USEEULER = 1
         IF ( USEEULER.NE.0 .AND. ABS(BZ).GT.1.0D-16 ) USEEULER = 1
      END SELECT
C
      IF ( USEEULER.EQ.1 ) THEN
C         prepare b field for rotation
         EB(1) = 0.D0
         EB(2) = 0.D0
         EB(3) = BFIELD
         IF ( IPARA.NE.0 ) IPARA = 3
         CALL EULER(BB,EUA,EUB,EUC)
         CALL ROTMATJB(EUA,EUB,EUC,ROT,ROTM1)
      END IF
C
      RETURN
99001 FORMAT (1x,'switched mcd off, b-field is zero')
99002 FORMAT (1x,'switched mcd=2 off, in plane b-field is zero')
      END
C*==rotbfield.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ROTBFIELD(USEEULER,MCD,EB,BB,ROT,ROTM1)
C     /****************************************************************/
C     purpose       : change direction of b-field for mcd calculations *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER MCD,USEEULER
      REAL*8 BB(3),EB(3)
      COMPLEX*16 ROT(MQD,MQD),ROTM1(MQD,MQD)
C
C Local variables
C
      REAL*8 BFIELD,BINPL,EUA,EUB,EUC,MU
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      WRITE (NOUT1,99001)
      WRITE (NOUT1,99002) BB
C
      IF ( MCD.EQ.2 ) THEN
C         special mcd -> rotate in plane b by pi/2
         MU = ATAN2(BB(2),BB(1)) + 2.D0*ATAN(1.D0)
         BINPL = SQRT(BB(1)**2+BB(2)**2)
         BB(1) = BINPL*COS(MU)
         BB(2) = BINPL*SIN(MU)
         BB(3) = BB(3)
      ELSE IF ( MCD.EQ.1 ) THEN
C         regular mcd -> invert b-field b -> -b
         BB(1) = -BB(1)
         BB(2) = -BB(2)
         BB(3) = -BB(3)
      END IF
      DO I = 1,3
         EB(I) = BB(I)
      END DO
C
      WRITE (NOUT1,99003) BB
C
      IF ( USEEULER.NE.0 ) THEN
C         new rotational matrix
         BFIELD = SQRT(BB(1)**2+BB(2)**2+BB(3)**2)
         EB(1) = 0.D0
         EB(2) = 0.D0
         EB(3) = BFIELD
         CALL EULER(BB,EUA,EUB,EUC)
         CALL ROTMATJB(EUA,EUB,EUC,ROT,ROTM1)
      END IF
C
      RETURN
99001 FORMAT (2x,'rotate b-field:')
99002 FORMAT (2x,'old b-field:',3(1x,e11.4))
99003 FORMAT (2x,'new b-field:',3(1x,e11.4))
      END
C*==euler.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE EULER(BB,EUA,EUB,EUC)
C     /****************************************************************/
C     # purpose:        calculate the euler angles for rotation        *
C                       by a given b-field (bx,by,bz)                  *
C     # note:           separate cases for some zero components are    *
C                       used because atan2 is not as accurate          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:EPS12
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EUA,EUB,EUC
      REAL*8 BB(3)
C
C Local variables
C
      REAL*8 ALPHA,BETA,BX,BXY,BY,BZ,DPI,GAMMA,MU,RD,THETA
      INTEGER IBX,IBY,IBZ
C
C*** End of declarations rewritten by SPAG
C
      DPI = 4.D0*ATAN(1.D0)
      RD = 180.D0/DPI
      BX = BB(1)
      BY = BB(2)
      BZ = BB(3)
      IBX = 0
      IBY = 0
      IBZ = 0
      IF ( ABS(BX).GT.EPS12 ) IBX = 1
      IF ( ABS(BY).GT.EPS12 ) IBY = 1
      IF ( ABS(BZ).GT.EPS12 ) IBZ = 1
C
      IF ( IBX.EQ.1 .AND. IBY.EQ.0 .AND. IBZ.EQ.0 ) THEN
C     b = (bx,0,0)
         GAMMA = DPI/2.D0
         BETA = DPI/2.D0
         ALPHA = DPI - (1.D0-SIGN(1.D0,BX))*DPI/2.D0
      ELSE IF ( IBX.EQ.0 .AND. IBY.EQ.1 .AND. IBZ.EQ.0 ) THEN
C     b = (0,by,0)
         GAMMA = DPI/2.D0
         BETA = DPI/2.D0
         ALPHA = DPI - SIGN(1.D0,BY)*DPI/2.D0
      ELSE IF ( IBX.EQ.0 .AND. IBY.EQ.0 .AND. IBZ.EQ.1 ) THEN
C     b = (0,0,bz)
         GAMMA = 0.D0
         BETA = (1.D0-SIGN(1.D0,BZ))*DPI/2.D0
         ALPHA = 0.D0
      ELSE IF ( IBX.EQ.1 .AND. IBY.EQ.1 .AND. IBZ.EQ.0 ) THEN
C     b = (bx,by,0)
         MU = DATAN(BY/BX)
         GAMMA = DPI/2.D0
         BETA = DPI/2.D0
         ALPHA = DPI - MU
      ELSE IF ( (IBX.EQ.1 .OR. IBY.EQ.1) .AND. IBZ.EQ.1 ) THEN
C     b = (bx,by,bz) or (bx,0,bz) or (0,by,bz)
         MU = DATAN2(BY,BX)
         BXY = SQRT(BX**2+BY**2)
         THETA = DATAN2(BZ,BXY)
         GAMMA = DPI/2.D0
         BETA = DPI/2.D0 - THETA
         ALPHA = DPI - MU
      ELSE
C     b = 0
         GAMMA = 0.D0
         BETA = 0.D0
         ALPHA = 0.D0
         WRITE (*,99001)
      END IF
C
C     euler angles for rotation
      EUA = ALPHA*RD
      EUB = BETA*RD
      EUC = GAMMA*RD
C
      IF ( IP.GE.1 ) WRITE (NOUT1,99002) EUA,EUB,EUC
C
      RETURN
C
99001 FORMAT ('b=0: euler set to zero')
99002 FORMAT ('  euler angles set to:',3(1x,f7.2))
      END
C*==rotmatjb.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ROTMATJB(EUA,EUB,EUC,ROT,ROTM1)
C     /****************************************************************/
C     # purpose      : determine general rotation matrix for integer   *
C                      and half integer j (rose 1967: elementary       *
C                      theory of angular momentum, p. 48 ff)           *
C                                                                      *
C     # parameter    :  eua, eub, euc = euler angles                   *
C                                                                      *
C     # uses the subroutines and functions:                            *
C       kapmue2k        factr      factorial                           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:ML,MQD,EPS12,CZERO,CIMAG
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EUA,EUB,EUC
      COMPLEX*16 ROT(MQD,MQD),ROTM1(MQD,MQD)
C
C Local variables
C
      REAL*8 A,B,C,J,M,MS,RD,SK,SQ,SR
      REAL*8 FACTORIAL,FACTR
      INTEGER I,ICOUNT,IJ,IM,IMS,IN1,IN2,IXC,IXS,K,K1,K2,K3,KAPPA
      INTEGER KAPMUE2K
      EXTERNAL FACTORIAL,FACTR,KAPMUE2K
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MQD
         DO K = 1,MQD
            ROT(I,K) = CZERO
            ROTM1(I,K) = CZERO
         END DO
      END DO
C
      RD = 4.D0*ATAN(1.D0)/180.D0
      A = EUA*RD
      B = EUB*RD
      C = EUC*RD
C
      DO KAPPA = -ML,ML - 1
         IF ( KAPPA.NE.0 ) THEN
            IJ = 2*IABS(KAPPA) - 1
            J = IJ/2.D0
            DO IM = -IJ,IJ,2
               M = IM/2.D0
               IN1 = KAPMUE2K(KAPPA,M,ML)
               DO IMS = -IJ,IJ,2
                  MS = IMS/2.D0
                  IN2 = KAPMUE2K(KAPPA,MS,ML)
C
                  SQ = DSQRT(FACTR(J+M)*FACTR(J-M)*FACTR(J+MS)
     &                 *FACTR(J-MS))
                  SK = 0.D0
                  DO K = 0,10
C                 sum over kappa:
C                 ===============
                     K1 = IDINT(J+M) - K
                     K2 = IDINT(J-MS) - K
                     K3 = IDINT(MS-M) + K
                     IF ( K1.GE.0 .AND. K2.GE.0 .AND. K3.GE.0 ) THEN
                        SR = (-1.D0)**K*SQ/(FACTORIAL(K)*FACTORIAL(K1)*
     &                       FACTORIAL(K2)*FACTORIAL(K3))
                        IXS = 2*K + IDINT(MS-M)
C                        IF ( B.EQ.0. .AND. IXS.EQ.0 ) THEN
                        IF ( ABS(B).LE.1.0D-16 .AND. IXS.EQ.0 ) THEN
                           SK = SK + SR
                        ELSE
                           IXC = IDINT(M-MS) + IJ - 2*K
                           SK = SK + SR*COS(B/2.D0)**IXC*(-SIN(B/2.D0))
     &                          **IXS
                        END IF
                     END IF
                  END DO
                  ROT(IN1,IN2) = SK*CDEXP(-CIMAG*(MS*A+M*C))
                  IF ( CDABS(ROT(IN1,IN2)).LT.EPS12 ) ROT(IN1,IN2)
     &                 = CZERO
               END DO
            END DO
         END IF
      END DO
C
      DO IN1 = 1,MQD
         DO IN2 = 1,MQD
            ROTM1(IN1,IN2) = DCONJG(ROT(IN2,IN1))
         END DO
      END DO
C
      IF ( IP.GT.4 ) THEN
         WRITE (NOUT1,99001)
         ICOUNT = 0
         DO IN1 = 1,MQD
            DO IN2 = 1,MQD
C               IF ( CDABS(ROT(IN1,IN2)).NE.0.0 ) THEN
               IF ( CDABS(ROT(IN1,IN2)).GT.1.0D-16 ) THEN
                  ICOUNT = ICOUNT + 1
                  WRITE (NOUT1,99002) ICOUNT,IN1,IN2,ROT(IN1,IN2),
     &                                ROT(IN2,IN1),ROTM1(IN1,IN2),
     &                                ROTM1(IN2,IN1)
               END IF
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (1x,'rotation matrix:   rot',63x,'rotm1')
99002 FORMAT (1x,i4,3x,i3,2x,i3,4x,2E15.7,2x,2E15.7,4x,2E15.7,2x,2E15.7)
      END
C*==factr.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION FACTR(X)
C     /****************************************************************/
C     # purpose      : calculate factorial n! for real arguments       *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 X
C
C Local variables
C
      INTEGER I
      REAL*8 RES
C
C*** End of declarations rewritten by SPAG
C
      IF ( X.LT.0. ) THEN
         RES = -1.D0
      ELSE
         RES = 1.D0
         DO I = 1,IDINT(X)
            RES = RES*DBLE(I)
         END DO
      END IF
      FACTR = RES
C
      END
C*==factorial.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION FACTORIAL(J)
C     /****************************************************************/
C     # purpose      : calculate factorial n! for integer argument n   *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER J
C
C Local variables
C
      INTEGER I
      REAL*8 RES
C
C*** End of declarations rewritten by SPAG
C
      IF ( J.LT.0 ) THEN
C         j!=infinity:
C         ============
         RES = -1.D0
      ELSE
         RES = 1.D0
         DO I = 1,J
            RES = RES*DBLE(I)
         END DO
      END IF
      FACTORIAL = RES
C
      END
