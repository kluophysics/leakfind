C*==sig0_optics.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG0_OPTICS(S10AQAB,S10BQAB,S10CQAB,S10DQAB,S2AQAB,
     &                       S2BQAB,S2CQAB,S2DQAB,S3AQAB,S3BQAB,S3CQAB,
     &                       S3DQAB,S4AQAB,S4BQAB,S4CQAB,S4DQAB,S10AQBA,
     &                       S10BQBA,S10CQBA,S10DQBA,S2AQBA,S2BQBA,
     &                       S2CQBA,S2DQBA,S3AQBA,S3BQBA,S3CQBA,S3DQBA,
     &                       S4AQBA,S4BQBA,S4CQBA,S4DQBA,SIG0Q_OPT,
     &                       NSPINPROJ,NOMEGA,OMEGATAB,ERYDA,ERYDB,
     &                       EPS_OPT,GAM_OPT,IEA,IEB)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *                J_m Im TAU00(t) J_n Im TAU00(t)                   *
C   *                                                                  *
C   * the prefactor 1/4 is included here in the S-terms                *
C   * constants are added in <SIG_SUM> when converting to SI units     *
C   *                                                                  *
C   ********************************************************************
C   *                         IZA       IZB                            *
C   *     METYPE = A:        2  (+)    2  (+)                          *
C   *              B:        2  (+)    1  (-)                          *
C   *              C:        1  (-)    2  (+)                          *
C   *              D:        1  (-)    1  (-)                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:WETAB,NETAB
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SITES,ONLY:NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,SIG_MODE
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C*--SIG0_OPTICS34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG0_OPTICS')
C
C Dummy arguments
C
      REAL*8 EPS_OPT,GAM_OPT
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IEA,IEB,NOMEGA,NSPINPROJ
      REAL*8 OMEGATAB(NOMEGA)
      COMPLEX*16 S10AQAB(3,3,NSPINPROJ,NQMAX),
     &           S10AQBA(3,3,NSPINPROJ,NQMAX),
     &           S10BQAB(3,3,NSPINPROJ,NQMAX),
     &           S10BQBA(3,3,NSPINPROJ,NQMAX),
     &           S10CQAB(3,3,NSPINPROJ,NQMAX),
     &           S10CQBA(3,3,NSPINPROJ,NQMAX),
     &           S10DQAB(3,3,NSPINPROJ,NQMAX),
     &           S10DQBA(3,3,NSPINPROJ,NQMAX),
     &           S2AQAB(3,3,NSPINPROJ,NQMAX),S2AQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S2BQAB(3,3,NSPINPROJ,NQMAX),
     &           S2BQBA(3,3,NSPINPROJ,NQMAX),S2CQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S2CQBA(3,3,NSPINPROJ,NQMAX),
     &           S2DQAB(3,3,NSPINPROJ,NQMAX),S2DQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S3AQAB(3,3,NSPINPROJ,NQMAX),
     &           S3AQBA(3,3,NSPINPROJ,NQMAX),S3BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S3BQBA(3,3,NSPINPROJ,NQMAX),
     &           S3CQAB(3,3,NSPINPROJ,NQMAX),S3CQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S3DQAB(3,3,NSPINPROJ,NQMAX),
     &           S3DQBA(3,3,NSPINPROJ,NQMAX),S4AQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4AQBA(3,3,NSPINPROJ,NQMAX),
     &           S4BQAB(3,3,NSPINPROJ,NQMAX),S4BQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S4CQAB(3,3,NSPINPROJ,NQMAX),
     &           S4CQBA(3,3,NSPINPROJ,NQMAX),S4DQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4DQBA(3,3,NSPINPROJ,NQMAX),
     &           SIG0Q_OPT(3,3,NSPINPROJ,NQMAX,NOMEGA)
C
C Local variables
C
      COMPLEX*16 ALF,AMIN,APLS,BET,BMIN,BPLS,SUM1_A(3,3),SUM1_B(3,3),
     &           SUM1_C(3,3),SUM1_D(3,3),SUM2_A(3,3),SUM2_B(3,3),
     &           SUM2_C(3,3),SUM2_D(3,3),W_AMIN,W_APLS,W_BMIN,W_BPLS,
     &           X1_A,X1_B,X1_C,X1_D,X2_A,X2_B,X2_C,X2_D
      INTEGER ILOOP_OBSE,IOBSE,IOMEGA,IQ
      REAL*8 X
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      WRITE (6,99001)
C
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      DO ILOOP_OBSE = 1,NSPR
         IOBSE = LIST_ISPR(ILOOP_OBSE)
         WRITE (6,99003)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(IOBSE)
         WRITE (6,99003)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
C
            WRITE (6,99002) IQ
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C
            SUM1_A(1:3,1:3) = S10AQAB(1:3,1:3,IOBSE,IQ)
     &                        + S2AQAB(1:3,1:3,IOBSE,IQ)
     &                        + S3AQAB(1:3,1:3,IOBSE,IQ)
     &                        + S4AQAB(1:3,1:3,IOBSE,IQ)
C
            SUM1_B(1:3,1:3) = S10BQAB(1:3,1:3,IOBSE,IQ)
     &                        + S2BQAB(1:3,1:3,IOBSE,IQ)
     &                        + S3BQAB(1:3,1:3,IOBSE,IQ)
     &                        + S4BQAB(1:3,1:3,IOBSE,IQ)
C
            SUM1_C(1:3,1:3) = S10CQAB(1:3,1:3,IOBSE,IQ)
     &                        + S2CQAB(1:3,1:3,IOBSE,IQ)
     &                        + S3CQAB(1:3,1:3,IOBSE,IQ)
     &                        + S4CQAB(1:3,1:3,IOBSE,IQ)
C
            SUM1_D(1:3,1:3) = S10DQAB(1:3,1:3,IOBSE,IQ)
     &                        + S2DQAB(1:3,1:3,IOBSE,IQ)
     &                        + S3DQAB(1:3,1:3,IOBSE,IQ)
     &                        + S4DQAB(1:3,1:3,IOBSE,IQ)
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C
            SUM2_A(1:3,1:3) = S10AQBA(1:3,1:3,IOBSE,IQ)
     &                        + S2AQBA(1:3,1:3,IOBSE,IQ)
     &                        + S3AQBA(1:3,1:3,IOBSE,IQ)
     &                        + S4AQBA(1:3,1:3,IOBSE,IQ)
C
            SUM2_B(1:3,1:3) = S10BQBA(1:3,1:3,IOBSE,IQ)
     &                        + S2BQBA(1:3,1:3,IOBSE,IQ)
     &                        + S3BQBA(1:3,1:3,IOBSE,IQ)
     &                        + S4BQBA(1:3,1:3,IOBSE,IQ)
C
            SUM2_C(1:3,1:3) = S10CQBA(1:3,1:3,IOBSE,IQ)
     &                        + S2CQBA(1:3,1:3,IOBSE,IQ)
     &                        + S3CQBA(1:3,1:3,IOBSE,IQ)
     &                        + S4CQBA(1:3,1:3,IOBSE,IQ)
C
            SUM2_D(1:3,1:3) = S10DQBA(1:3,1:3,IOBSE,IQ)
     &                        + S2DQBA(1:3,1:3,IOBSE,IQ)
     &                        + S3DQBA(1:3,1:3,IOBSE,IQ)
     &                        + S4DQBA(1:3,1:3,IOBSE,IQ)
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            IF ( IPRINT.GT.3 ) THEN
C
               X = 0.25D0
               WRITE (6,99004) 'S10AQAB',IQ
               CALL PR_COND_TENSOR(S10AQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S2AQAB',IQ
               CALL PR_COND_TENSOR(S2AQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S3AQAB',IQ
               CALL PR_COND_TENSOR(S3AQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S4AQAB',IQ
               CALL PR_COND_TENSOR(S4AQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S10BQAB',IQ
               CALL PR_COND_TENSOR(S10BQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S2BQAB',IQ
               CALL PR_COND_TENSOR(S2BQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S3BQAB',IQ
               CALL PR_COND_TENSOR(S3BQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S4BQAB',IQ
               CALL PR_COND_TENSOR(S4BQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S10DQAB',IQ
               CALL PR_COND_TENSOR(S10DQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S2DQAB',IQ
               CALL PR_COND_TENSOR(S2DQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S3DQAB',IQ
               CALL PR_COND_TENSOR(S3DQAB(1:3,1:3,IOBSE,IQ),X)
               WRITE (6,99004) 'S4DQAB',IQ
               CALL PR_COND_TENSOR(S4DQAB(1:3,1:3,IOBSE,IQ),X)
C
               WRITE (6,*) ' '
               WRITE (6,99004) 'SUM1_A',IQ
               CALL PR_COND_TENSOR(SUM1_A,X)
               WRITE (6,99004) 'SUM1_B',IQ
               CALL PR_COND_TENSOR(SUM1_B,X)
               WRITE (6,99004) 'SUM1_C',IQ
               CALL PR_COND_TENSOR(SUM1_C,X)
               WRITE (6,99004) 'SUM1_D',IQ
               CALL PR_COND_TENSOR(SUM1_D,X)
               WRITE (6,99004) 'SUM2_A',IQ
               CALL PR_COND_TENSOR(SUM2_A,X)
               WRITE (6,99004) 'SUM2_B',IQ
               CALL PR_COND_TENSOR(SUM2_B,X)
               WRITE (6,99004) 'SUM2_C',IQ
               CALL PR_COND_TENSOR(SUM2_C,X)
               WRITE (6,99004) 'SUM2_D',IQ
               CALL PR_COND_TENSOR(SUM2_D,X)
C
            END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
            APLS = ERYDA
            AMIN = DCONJG(APLS)
            W_APLS = WETAB(IEA,1)
            W_AMIN = DCONJG(W_APLS)
            BPLS = ERYDB
            BMIN = DCONJG(BPLS)
            W_BPLS = WETAB(IEB,2)
            W_BMIN = DCONJG(W_BPLS)
C
C-----------------------------------------------------------------------
            LOOP_IOMEGA:DO IOMEGA = 1,NOMEGA
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
               IF ( SIG_MODE.EQ.'OPTICS-ABS' ) THEN
C
                  IF ( (IEA-IEB+NETAB(2)).NE.IOMEGA ) CYCLE LOOP_IOMEGA
C
                  ALF = APLS - BPLS - CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA)
C
                  X1_A = 0.25D0*W_BPLS/OMEGATAB(IOMEGA)
                  X2_A = X1_A
                  X1_B = X1_A
                  X2_B = X1_A
C
                  X1_C = 0.0D0
                  X2_C = 0.0D0
                  X1_D = 0.0D0
                  X2_D = 0.0D0
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
               ELSE
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
                  ALF = APLS - BPLS - CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + BPLS - APLS + CI*GAM_OPT
                  X1_A = 0.25D0*W_APLS*W_BPLS/(ALF*BET)
C
                  ALF = APLS - BPLS + CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + APLS - BPLS + CI*GAM_OPT
                  X2_A = 0.25D0*W_APLS*W_BPLS/(ALF*BET)
C
                  ALF = APLS - BMIN - CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + BMIN - APLS + CI*GAM_OPT
                  X1_B = 0.25D0*W_APLS*W_BMIN/(ALF*BET)
C
                  ALF = APLS - BMIN + CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + APLS - BMIN + CI*GAM_OPT
                  X2_B = 0.25D0*W_APLS*W_BMIN/(ALF*BET)
C
                  ALF = AMIN - BPLS - CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + BPLS - AMIN + CI*GAM_OPT
                  X1_C = 0.25D0*W_AMIN*W_BPLS/(ALF*BET)
C
                  ALF = AMIN - BPLS + CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + AMIN - BPLS + CI*GAM_OPT
                  X2_C = 0.25D0*W_AMIN*W_BPLS/(ALF*BET)
C
                  ALF = AMIN - BMIN - CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + BMIN - AMIN + CI*GAM_OPT
                  X1_D = 0.25D0*W_AMIN*W_BMIN/(ALF*BET)
C
                  ALF = AMIN - BMIN + CI*EPS_OPT
                  BET = OMEGATAB(IOMEGA) + AMIN - BMIN + CI*GAM_OPT
                  X2_D = 0.25D0*W_AMIN*W_BMIN/(ALF*BET)
C
               END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
               SIG0Q_OPT(1:3,1:3,IOBSE,IQ,IOMEGA)
     &            = SIG0Q_OPT(1:3,1:3,IOBSE,IQ,IOMEGA)
     &            + X1_A*SUM1_A(1:3,1:3) - X1_B*SUM1_B(1:3,1:3)
     &            - X1_C*SUM1_C(1:3,1:3) + X1_D*SUM1_D(1:3,1:3)
     &            + X2_A*SUM2_A(1:3,1:3) - X2_B*SUM2_B(1:3,1:3)
     &            - X2_C*SUM2_C(1:3,1:3) + X2_D*SUM2_D(1:3,1:3)
C
            END DO LOOP_IOMEGA
C-----------------------------------------------------------------------
C
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      END DO ! ILOOP_OBSE
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
99001 FORMAT (//,1X,79('*'),/,34X,'<SIG0_OPTICS>',/,1X,79('*'),/)
99002 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99003 FORMAT (/,37('*'),' SIG0_OPTICS ',37('*'))
99004 FORMAT (/,'#### <SIG0_OPTICS> ',A,' IQ:',I3)
      END
