C*==strcc.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE STRCC(ERYD,RESTART)
C   ********************************************************************
C   *                                                                  *
C   *    CALCULATE ALL QUANTITIES WHICH DEPEND ON THE ENERGY           *
C   *                                                                  *
C   *    the  terms  QQMLRS( MMLL, S, IQQP )  set up in  <STRAA>       *
C   *    are multiplied by  IILERS(LL,S,IQQP)  to save storage         *
C   *    this term is stored and removed at the next call of           *
C   *    <STRCC>                                                       *
C   *                                                                  *
C   *    RESTART = .T. : ICALL is set to 0 and the routine is left     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_STR,ONLY:ALPHA0,D1TERM3,QQMLRS,GGJLRS,IILERS,NSDL,SMAX,
     &    LLMAX,LLARR,J13MAX,J22MAX,D300,ETA,EDU,PWP,STRMODE,NQQP_STR_CC
      USE MOD_CONSTANTS,ONLY:C0,CI,PI
      IMPLICIT NONE
C*--STRCC21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      LOGICAL RESTART
C
C Local variables
C
      REAL*8 ALPHA,FAKTOR
      COMPLEX*16 CAUX,EHOCHJ,EPWMLLH(0:LLARR),PDU
      INTEGER ICALL,IQQP,J13,J22,LL,LL_MIN_J22,MM,MMLL,S
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
C  ===============================================================
C                        interpolation
C  ===============================================================
      IF ( STRMODE.EQ.2 ) THEN
         CALL STRCCLAG(ERYD)
         RETURN
      END IF
C  ===============================================================
C
      IF ( RESTART ) THEN
         ICALL = 0
         RETURN
      END IF
C
      IF ( .NOT.RESTART ) THEN
C
         EDU = ERYD/(2*PI/ALAT)**2
         PDU = CDSQRT(EDU)
C
         CAUX = PDU
         MMLL = 0
         DO LL = 0,LLMAX
            CAUX = CAUX/PDU
            DO MM = -LL, + LL
               MMLL = MMLL + 1
               PWP(MMLL) = CAUX
            END DO
         END DO
C
C  ===============================================================
C                      ********
C                      * DLM1 *
C                      ********
C  /D (6)/
C
         EPWMLLH(0) = 1.0D0
         D1TERM3(0) = CDEXP(EDU/ETA)
C
         DO LL = 1,LLMAX
            EPWMLLH(LL) = EPWMLLH(LL-1)*CI/PDU
            D1TERM3(LL) = D1TERM3(LL-1)/PDU
         END DO
C
C  ===============================================================
C                      ********
C                      * DLM2 *
C                      ********
C
C  ---------------------------------------------------------------
         IF ( ICALL.EQ.0 ) THEN
            CALL CINIT((1+LLARR)*NSDL*NQQP_STR_CC,IILERS)
         ELSE
C     remove the energy-dependent factor IILERS from the last run
C
            DO IQQP = 1,NQQP_STR_CC
               DO S = 1,SMAX(IQQP)
C
                  MMLL = 0
                  DO LL = 0,LLMAX
C
                     DO MM = -LL, + LL
                        MMLL = MMLL + 1
                        QQMLRS(MMLL,S,IQQP) = QQMLRS(MMLL,S,IQQP)
     &                     /IILERS(LL,S,IQQP)
                     END DO
                  END DO
               END DO
            END DO
         END IF
C  ---------------------------------------------------------------
C
         ICALL = ICALL + 1
C
         DO IQQP = 1,NQQP_STR_CC
            DO S = 1,SMAX(IQQP)
               DO LL = 0,LLMAX
                  FAKTOR = 1.0D0
                  EHOCHJ = 1.0D0
                  IILERS(LL,S,IQQP) = 0.0D0
                  DO J22 = 0,J22MAX
                     LL_MIN_J22 = LL - J22
                     IILERS(LL,S,IQQP) = IILERS(LL,S,IQQP)
     &                  + EHOCHJ*(GGJLRS(LL_MIN_J22,S,IQQP)*FAKTOR)
                     FAKTOR = FAKTOR/(ETA*(J22+1.0D0))
                     EHOCHJ = EHOCHJ*EDU
                  END DO
C
C     /D (20) AND (22)/
                  IILERS(LL,S,IQQP) = IILERS(LL,S,IQQP)*EPWMLLH(LL)
               END DO
            END DO
         END DO
C
C  ---------------------------------------------------------------
C     multiply the energy-dependent factor IILERS for the
C     current energy  ERYD
C
         DO IQQP = 1,NQQP_STR_CC
            DO S = 1,SMAX(IQQP)
C
               MMLL = 0
               DO LL = 0,LLMAX
C
                  DO MM = -LL, + LL
                     MMLL = MMLL + 1
                     QQMLRS(MMLL,S,IQQP) = QQMLRS(MMLL,S,IQQP)
     &                  *IILERS(LL,S,IQQP)
                  END DO
               END DO
            END DO
         END DO
C
C  ===============================================================
C                      ********
C                      * DLM3 *
C                      ********
C    /D (13)/
C
         D300 = 0.0D0
         EHOCHJ = 1.0D000
         ALPHA = ALPHA0
         J13 = -1
 50      CONTINUE
         J13 = J13 + 1
         D300 = ALPHA*EHOCHJ + D300
         ALPHA = ALPHA*(2.0D0*J13-1.0D0)
     &           /(ETA*(J13+1.0D0)*(2.0D0*J13+1.0D0))
         EHOCHJ = EHOCHJ*EDU
C     prevent floating point underflow
         IF ( ABS(EHOCHJ).LT.1D-50 ) EHOCHJ = C0
C
         IF ( CDABS(ALPHA*EHOCHJ/D300).GT.1.0D-10 ) GOTO 50
         IF ( J13.LT.J13MAX ) GOTO 50
C
      END IF
C
      END
