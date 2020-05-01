C*==scfchkcore.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFCHKCORE(IT,IM,IRTOP,IPRINT,RC,DRDIC,DOVRC,VV,BB,FC,
     &                      GC,DP,DQ,WP,WQ,Z,NRC,NSOL,ECOR,KAPCOR,NQN,
     &                      LQN,NRCMAX)
C   ********************************************************************
C   *                                                                  *
C   *   find the energy eigenvalues for core states                    *
C   *   with  quantum numbers  NQN  LQN  for KAPPA = L +/- 1/2         *
C   *   for   a spin-independent potential                             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:R2DRDI
      USE MOD_ANGMOM,ONLY:TXT_L,TXT_J
      IMPLICIT NONE
C*--SCFCHKCORE15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 C,UNEND,TOLVAR
      PARAMETER (C=2.0D0*1.370373D+02,UNEND=600.0D0,TOLVAR=1.0D-6)
      INTEGER ITERMAX
      PARAMETER (ITERMAX=200)
C
C Dummy arguments
C
      INTEGER IM,IPRINT,IRTOP,IT,LQN,NQN,NRC,NRCMAX,NSOL,Z
      REAL*8 BB(NRCMAX),DOVRC(NRCMAX),DP(2,2,NRCMAX),DQ(2,2,NRCMAX),
     &       DRDIC(NRCMAX),ECOR(2),FC(2,2,NRCMAX),GC(2,2,NRCMAX),
     &       RC(NRCMAX),VV(NRCMAX),WP(2,2,NRCMAX),WQ(2,2,NRCMAX)
      INTEGER KAPCOR(2)
C
C Local variables
C
      REAL*8 CGD(2),CGMD(2),CGO,DEC,EC,ELIM,MJ,NORM,PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),VAL
      INTEGER IR,IR_MATCH,IR_ZERO,IR_ZERO_TRIED,ISTART,ITER,KAP(2),KAP1,
     &        KAP2,L,LLL,MUEM05,NN,NODE,NPOSEWARN,S
C
C*** End of declarations rewritten by SPAG
C
      NPOSEWARN = 0
C
C                       --------------------------
C                       INITIALIZE QUANTUM NUMBERS
C                       --------------------------
C
      L = LQN
      MUEM05 = -1
      MJ = MUEM05 + 0.5D0
C
      KAP1 = -L - 1
      KAP2 = L
      IF ( L.EQ.0 ) KAP2 = KAP1
      KAP(1) = KAP1
      KAP(2) = KAP2
C
      LLL = L*(L+1)
      IF ( ABS(MJ).GT.L ) THEN
         NSOL = 1
      ELSE
         NSOL = 2
      END IF
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO S = 1,NSOL
         IR_ZERO_TRIED = 0
C                                      --------------------
C                                         FIND  E-LIMIT
C                                      --------------------
         IF ( LLL.EQ.0 ) THEN
            ELIM = -2*DBLE(Z*Z)/(1.5D0*NQN*NQN)
         ELSE
            ELIM = VV(1) + LLL/RC(1)**2
            DO IR = 2,NRC
               VAL = VV(IR) + LLL/RC(IR)**2
               IF ( VAL.LE.ELIM ) ELIM = VAL
            END DO
         END IF
C
         EC = -DBLE(Z*Z)/(2.0D0*NQN*NQN)
C
         ISTART = 1
 50      CONTINUE
         IF ( EC.LE.ELIM ) EC = ELIM*0.7D0
C
C                                      --------------------
C                                         FIND    IR_ZERO
C                                      --------------------
         DO IR = 1,(NRC-1)
            IF ( (VV(IR)-EC)*RC(IR)**2.GT.UNEND ) THEN
               IF ( MOD(IR,2).EQ.0 ) THEN
                  IR_ZERO = IR + 1
               ELSE
                  IR_ZERO = IR
               END IF
               GOTO 100
            END IF
         END DO
         IR_ZERO = NRC - 1
         IR_ZERO_TRIED = IR_ZERO_TRIED + 1
         IF ( IR_ZERO_TRIED.GT.20 ) THEN
            NSOL = -1
            RETURN
         END IF
C                                      --------------------
C                                         FIND    NMATCH
C                                      --------------------
 100     CONTINUE
         IR = IR_ZERO + 1
         DO NN = 1,IR_ZERO
            IR = IR - 1
            IF ( (VV(IR)+LLL/RC(IR)**2-EC).LT.0.0D0 ) THEN
               IR_MATCH = IR
               GOTO 150
            END IF
         END DO
         IF ( IPRINT.GT.0 ) WRITE (6,99004) IT,NQN,L,EC
C
 150     CONTINUE
         CALL CORE_DIR_ABM(IT,C,EC,L,MJ,'OUT',VV,BB,RC,DRDIC,DOVRC,
     &                     IR_MATCH,IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                     PIW,QIW,CGD,CGMD,CGO,NRC)
C
         NODE = 0
         DO IR = 2,IR_MATCH
            IF ( GC(S,S,IR)*GC(S,S,IR-1).LT.0.0D0 ) NODE = NODE + 1
         END DO
C
         IF ( IPRINT.GE.2 ) WRITE (6,99003) IT,NQN,L,KAP(S),(2*MUEM05+1)
     &                             ,0,EC,IR_MATCH,RC(IR_MATCH),IR_ZERO,
     &                             RC(IR_ZERO),NODE,
     &                             (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1))
C
         IF ( NODE.NE.(NQN-L-1) ) THEN
            IF ( NODE.GT.(NQN-L-1) ) THEN
               EC = 1.2D0*EC
            ELSE
               EC = 0.8D0*EC
            END IF
            GOTO 50
         ELSE IF ( (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1).LE.0.0D0) .OR. 
     &             (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-1).GE.1.0D0) ) THEN
            EC = 0.9D0*EC
            GOTO 50
         END IF
C
         CALL CORE_DIR_ABM(IT,C,EC,L,MJ,'INW',VV,BB,RC,DRDIC,DOVRC,
     &                     IR_MATCH,IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                     PIW,QIW,CGD,CGMD,CGO,NRC)
C
         ITER = 0
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 200     CONTINUE
         ITER = ITER + 1
C                         ----------------------------------
C                         CHECK WHETHER NUMBER OF NODES O.K.
C                         ----------------------------------
         IF ( ITER.GT.1 ) THEN
            NODE = 0
            DO IR = 2,(IR_MATCH-1)
               IF ( GC(S,S,IR)*GC(S,S,IR-1).LT.0.0D0 ) NODE = NODE + 1
            END DO
            IF ( IPRINT.GE.3 ) WRITE (6,99003) IT,NQN,L,KAP(S),
     &                                (2*MUEM05+1),ITER,EC,IR_MATCH,
     &                                RC(IR_MATCH),IR_ZERO,RC(IR_ZERO),
     &                                NODE,
     &                                (GC(S,S,IR_MATCH)/GC(S,S,IR_MATCH-
     &                                1))
C
            IF ( NODE.NE.(NQN-L-1) ) THEN
               IF ( NODE.GT.(NQN-L-1) ) THEN
                  EC = 1.2D0*EC
               ELSE
                  EC = 0.8D0*EC
               END IF
               ISTART = ISTART + 1
               IF ( ISTART.LT.20 ) GOTO 50
            END IF
         END IF
C
C                          --------------------------------
C                             CALCULATE THE EIGENVALUE
C                          USING THE CONVENTIONAL ALGORITHM
C                          --------------------------------
C
         CALL CORE_DIR_ABM(IT,C,EC,L,MJ,'OUT',VV,BB,RC,DRDIC,DOVRC,
     &                     IR_MATCH,IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                     PIW,QIW,CGD,CGMD,CGO,NRC)
C
         CALL CORE_DIR_ABM(IT,C,EC,L,MJ,'INW',VV,BB,RC,DRDIC,DOVRC,
     &                     IR_MATCH,IR_ZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,
     &                     PIW,QIW,CGD,CGMD,CGO,NRC)
C
         NORM = POW(S,S)/PIW(S,S)
         DO IR = IR_MATCH,IR_ZERO
            GC(S,S,IR) = GC(S,S,IR)*NORM
            FC(S,S,IR) = FC(S,S,IR)*NORM
         END DO
C
         NORM = 0.0D0
         DO IR = 3,MIN(IR_ZERO,IRTOP),2
            NORM = NORM + R2DRDI(IR,IM)*(GC(S,S,IR)**2+FC(S,S,IR)**2)
     &             + 4.D0*R2DRDI(IR-1,IM)
     &             *(GC(S,S,IR-1)**2+FC(S,S,IR-1)**2) + R2DRDI(IR-2,IM)
     &             *(GC(S,S,IR-2)**2+FC(S,S,IR-2)**2)
         END DO
         NORM = NORM/3.0D0
C
         DEC = POW(S,S)*(QOW(S,S)-RC(IR_MATCH)*C*FC(S,S,IR_MATCH))/NORM
C
         EC = EC + DEC
C
         IF ( EC.GT.0.0D0 ) THEN
            IF ( IPRINT.GE.1 ) WRITE (6,*) ' warning from <CORE> E=',EC,
     &                                IT,NQN,L
            EC = -0.2D0
            NPOSEWARN = NPOSEWARN + 1
            IF ( NPOSEWARN.GT.5 ) THEN
               NSOL = -1
               RETURN
            END IF
         END IF
C
         IF ( IPRINT.GE.2 ) WRITE (6,99002) EC,DEC
C
C--------------------------------------  check relative change in energy
C
         IF ( ITER.LT.ITERMAX ) THEN
            IF ( ABS(DEC/EC).GT.TOLVAR ) GOTO 200
C
            IF ( IPRINT.GT.0 ) WRITE (6,99001) IT,NQN,TXT_L(L),
     &                                TXT_J(IABS(KAP(S))),(2*MUEM05+1),
     &                                KAP(S),ITER,EC
C
            KAPCOR(S) = KAP(S)
            ECOR(S) = EC
         ELSE
            NSOL = -1
            RETURN
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
99001 FORMAT (2I4,A1,A4,I3,'/2',2I4,2X,F15.8,F17.3,:,F17.3,/)
99002 FORMAT (' E=',F14.7,' CORR ',E11.4)
99003 FORMAT (/,' IT=',I2,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,
     &        '/2',/,' E(',I2,')   =',F15.5,/,' NMATCH  =',I5,'    R=',
     &        F10.5,/,' IR_ZERO =',I5,'    R=',F10.5,/,' NODES   =',I5,
     &        '  RAT=',E11.4)
99004 FORMAT (//,'  STOP IN <<SCFCHKCORE>>',/,'  IT=',I2,' NQN=',I2,
     &        ' L=',I2,/,'  no matching-radius found for  EC=',F10.3)
      END
