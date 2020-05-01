C*==spec_init_transpho.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_INIT_TRANSPHO(NKMMAX,NLMAX)
C
C
      USE MOD_SPEC_TRANSFO,ONLY:NB,REL,SB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKMMAX,NLMAX
C
C Local variables
C
      INTEGER KAP(:),MUD(:),NR(:),NRL(:),NRM(:)
      REAL*8 NRS(:)
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATABLE KAP,MUD,NR,NRL,NRM,NRS
      ALLOCATE (NB(NKMMAX),SB(NKMMAX),NR(NKMMAX),KAP(NKMMAX))
      ALLOCATE (REL(NKMMAX,2),MUD(NKMMAX),NRL(NKMMAX),NRM(NKMMAX))
      ALLOCATE (NRS(NKMMAX))
C
      CALL TRANS_HEJB1(NLMAX,NKMMAX,KAP,MUD,NRL,NRM,NRS,NR,REL)
      CALL TRANS_HEJB2(NB,SB,NKMMAX,NLMAX)
C
      END
C*==spec_transpho.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
      SUBROUTINE SPEC_TRANSPHO(CM,CB,NKMMAX,MODE,P)
C
      USE MOD_SPEC_TRANSFO,ONLY:NB,REL,SB
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER MODE,NKMMAX
      COMPLEX*16 P
      COMPLEX*16 CB(NKMMAX,NKMMAX),CM(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,I11,J,J11
C
C*** End of declarations rewritten by SPAG
C
C
      DO I = 1,NKMMAX
         I11 = REL(NB(I),SB(I))
         DO J = 1,NKMMAX
            J11 = REL(NB(J),SB(J))
            CB(I11,J11) = CM(I,J)
         END DO
      END DO
C
      IF ( MODE.EQ.1 ) CB(:,:) = -CI*CB(:,:)*P
C
      END
C*==trans_hejb1.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE TRANS_HEJB1(ML,MQD,KAP,MUD,NRL,NRM,NRS,NR,REL)
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
C*==trans_hejb2.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
      SUBROUTINE TRANS_HEJB2(NB,SB,MQD,ML)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ML,MQD
      INTEGER NB(MQD),SB(MQD)
C
C Local variables
C
      INTEGER A,C,D,E,F,I,J,MLS,S
C
C*** End of declarations rewritten by SPAG
C
C     /****************************************************************/
C     # purpose      :                                                 *
C       transforms kappa = -4, -3,...., 3, 4 arrangement sb(mqd) in
C                  kappa = -1, 1, -2, 2, ... arrangement nb(mqd)
C
C     /****************************************************************/
C     /* output */
C
C     /* local */
C
      E = 1
      C = 0
      A = 1
      S = 2
      D = 0
      MLS = 2*ML*ML
      F = 2*ML
      DO I = 1,F
         S = 3 - S
         DO J = 1,A
            D = D + 1
            IF ( D.GT.MLS ) GOTO 99999
            IF ( S.EQ.2 ) THEN
               C = C + 1
               NB(D) = C
            ELSE IF ( S.EQ.1 ) THEN
               IF ( J.LE.INT(A/2) ) NB(D) = E + INT(A/2+1)
               IF ( J.GT.INT(A/2) ) NB(D) = E - INT(A/2)
               E = E + 1
            END IF
            SB(D) = 3 - S
C           sb(d) = s
         END DO
         A = A + 1
      END DO
C
99999 CONTINUE
      END
C
