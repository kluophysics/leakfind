C*==sig1_optics.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_OPTICS(MAQAB,MBQAB,MCQAB,MDQAB,MAQBA,MCQBA,MBQBA,
     &                       MDQBA,NSPINPROJ,NOMEGA,OMEGATAB,ERYDA,
     &                       ERYDB,EPS_OPT,GAM_OPT,IEA,IEB,SIG1Q_OPT,
     &                       KEY)
C   ********************************************************************
C   *                                                                  *
C   *   calculate  TRACE jbar(mue,z2,z1)*mat*jbar(nue,z1,z2)           *
C   *   (including Vertex-corrections)                                 *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   *        IPERT = 1 i.e. the perturbation is represented by the     *
C   *        electrical currect density operator                       *
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
      USE MOD_LINRESP,ONLY:CHIZ
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,LMAT3
      USE MOD_SITES,ONLY:NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,SIG_MODE
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--SIG1_OPTICS39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG1_OPTICS')
      INTEGER IPERT
      PARAMETER (IPERT=1)
C
C Dummy arguments
C
      REAL*8 EPS_OPT,GAM_OPT
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IEA,IEB,NOMEGA,NSPINPROJ
      CHARACTER*1 KEY
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           SIG1Q_OPT(3,3,NSPINPROJ,NQMAX,NQMAX,NOMEGA)
      REAL*8 OMEGATAB(NOMEGA)
C
C Local variables
C
      COMPLEX*16 ALF,AMIN,APLS,BET,BMIN,BPLS,L_MCQAB(:,:,:,:,:),
     &           L_MCQBA(:,:,:,:,:),MBQAB_L(:,:,:,:,:),
     &           MBQBA_L(:,:,:,:,:),MDQABX(:,:,:,:,:),MDQAB_L(:,:,:,:,:)
     &           ,MDQBAX(:,:,:,:,:),MDQBA_L(:,:,:,:,:),S1_A,S1_B,S1_C,
     &           S1_D,S2_A,S2_B,S2_C,S2_D,SUM1_A(3,3),SUM1_B(3,3),
     &           SUM1_C(3,3),SUM1_D(3,3),SUM2_A(3,3),SUM2_B(3,3),
     &           SUM2_C(3,3),SUM2_D(3,3),W_AMIN,W_APLS,W_BMIN,W_BPLS,
     &           X1_A,X1_B,X1_C,X1_D,X2_A,X2_B,X2_C,X2_D
      INTEGER I,ILOOP_OBSE,IOBSE,IOMEGA,IQ,IQCHI,J,JQ,JQCHI,K1,K2,L1,L2,
     &        L3,L4,MUE,NKMSQ,NUE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MDQBAX,MDQABX
      ALLOCATABLE MBQAB_L,MBQBA_L,MDQAB_L,MDQBA_L
      ALLOCATABLE L_MCQBA,L_MCQAB
C
      CALL TRACK_INFO(ROUTINE)
C
      NKMSQ = NKM*NKM
C
      ALLOCATE (MDQBAX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQABX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
      ALLOCATE (MBQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MBQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (L_MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (L_MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
      IF ( KEY.EQ.'N' ) WRITE (6,99002) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                  'without vertex-corrections'
      IF ( KEY.EQ.'V' ) WRITE (6,99002) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                  'including vertex-corrections'
C
C-----------------------------------------------------------------------
C  multiply the averaged MEs with proper LMAT to take into account
C  that CHIZ does not incude the LMATs
C-----------------------------------------------------------------------
C
      DO ILOOP_OBSE = 1,NSPR
         IOBSE = LIST_ISPR(ILOOP_OBSE)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               MBQAB_L(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(MBQAB(:,:,MUE,IOBSE,IQ),LMAT3)
               MBQBA_L(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(MBQBA(:,:,MUE,IOBSE,IQ),LMAT3)
C
               L_MCQAB(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(LMAT3,MCQAB(:,:,MUE,IOBSE,IQ))
               L_MCQBA(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(LMAT3,MCQBA(:,:,MUE,IOBSE,IQ))
C
               MDQAB_L(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQAB(:,:,MUE,IOBSE,IQ),LMAT3))
               MDQBA_L(:,:,MUE,IOBSE,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQBA(:,:,MUE,IOBSE,IQ),LMAT3))
C
            END DO
         END DO
      END DO
C
      DO ILOOP_OBSE = 1,NSPR
         IOBSE = LIST_ISPR(ILOOP_OBSE)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               DO I = 1,NKM
                  DO J = 1,NKM
                     MDQBAX(I,J,MUE,IOBSE,IQ)
     &                  = DCONJG(MDQBA_L(J,I,MUE,IOBSE,IQ))
C
                     MDQABX(I,J,MUE,IOBSE,IQ)
     &                  = DCONJG(MDQAB_L(J,I,MUE,IOBSE,IQ))
                  END DO
               END DO
C
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
C
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      DO ILOOP_OBSE = 1,NSPR
         IOBSE = LIST_ISPR(ILOOP_OBSE)
C
         WRITE (6,99003)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(IOBSE)
         WRITE (6,99003)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
            IQCHI = IQ - IQBOT_CHI + 1
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            DO JQ = IQBOT_CHI,IQTOP_CHI
               JQCHI = JQ - IQBOT_CHI + 1
C
               WRITE (6,99001) IQ,JQ
C
               DO NUE = 1,3
C
                  DO MUE = 1,3
C
                     S1_A = C0
                     S1_B = C0
                     S1_C = C0
                     S1_D = C0
                     S2_A = C0
                     S2_B = C0
                     S2_C = C0
                     S2_D = C0
C
                     K1 = (IQCHI-1)*NKMSQ
                     DO L1 = 1,NKM
                        DO L4 = 1,NKM
                           K1 = K1 + 1
C
                           K2 = (JQCHI-1)*NKMSQ
                           DO L2 = 1,NKM
                              DO L3 = 1,NKM
                                 K2 = K2 + 1
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C
                                 S1_A = S1_A + MAQBA(L4,L1,MUE,IOBSE,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MAQAB(L2,L3,NUE,IPERT,JQ)
C
                                 S1_B = S1_B + 
     &                                  L_MCQBA(L4,L1,MUE,IOBSE,IQ)
     &                                  *CHIZ(K1,K2,2)
     &                                  *MBQAB_L(L2,L3,NUE,IPERT,JQ)
C
                                 S1_C = S1_C + 
     &                                  L_MCQAB(L4,L1,MUE,IPERT,IQ)
     &                                  *CHIZ(K1,K2,2)
     &                                  *MBQBA_L(L2,L3,NUE,IOBSE,JQ)
C
                                 S1_D = S1_D + 
     &                                  MDQBAX(L4,L1,MUE,IOBSE,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MDQABX(L2,L3,NUE,IPERT,JQ)
C
C11111111111111111111111111111111111111111111111111111111111111111111111
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C
                                 S2_A = S2_A + MAQBA(L4,L1,MUE,IPERT,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MAQAB(L2,L3,NUE,IOBSE,JQ)
C
                                 S2_B = S2_B + 
     &                                  L_MCQBA(L4,L1,MUE,IPERT,IQ)
     &                                  *CHIZ(K1,K2,2)
     &                                  *MBQAB_L(L2,L3,NUE,IOBSE,JQ)
C
                                 S2_C = S2_C + 
     &                                  L_MCQAB(L4,L1,MUE,IOBSE,IQ)
     &                                  *CHIZ(K1,K2,2)
     &                                  *MBQBA_L(L2,L3,NUE,IPERT,JQ)
C
                                 S2_D = S2_D + 
     &                                  MDQBAX(L4,L1,MUE,IPERT,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MDQABX(L2,L3,NUE,IOBSE,JQ)
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C
                              END DO
C
                           END DO
                        END DO
                     END DO
C
                     SUM1_A(MUE,NUE) = S1_A
                     SUM1_B(MUE,NUE) = S1_B
                     SUM1_C(MUE,NUE) = S1_C
                     SUM1_D(MUE,NUE) = DCONJG(S1_D)
C
                     SUM2_A(MUE,NUE) = S2_A
                     SUM2_B(MUE,NUE) = S2_B
                     SUM2_C(MUE,NUE) = S2_C
                     SUM2_D(MUE,NUE) = DCONJG(S2_D)
C
                  END DO
               END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               IF ( IPRINT.GE.3 ) THEN
                  WRITE (6,*) 'SUMS 1'
                  WRITE (6,*) 'SUM1_A'
                  WRITE (6,99004) ((SUM1_A(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM1_B'
                  WRITE (6,99004) ((SUM1_B(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM1_C'
                  WRITE (6,99004) ((SUM1_C(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM1_D'
                  WRITE (6,99004) ((SUM1_D(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUMS 2'
                  WRITE (6,*) 'SUM2_A'
                  WRITE (6,99004) ((SUM2_A(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM2_B'
                  WRITE (6,99004) ((SUM2_B(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM2_C'
                  WRITE (6,99004) ((SUM2_C(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM2_D'
                  WRITE (6,99004) ((SUM2_D(MUE,NUE),NUE=1,3),MUE=1,3)
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
                     IF ( (IEA-IEB+NETAB(2)).NE.IOMEGA )
     &                    CYCLE LOOP_IOMEGA
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
                  SIG1Q_OPT(:,:,IOBSE,IQ,JQ,IOMEGA)
     &               = SIG1Q_OPT(:,:,IOBSE,IQ,JQ,IOMEGA)
     &               + X1_A*SUM1_A(1:3,1:3) - X1_B*SUM1_B(1:3,1:3)
     &               - X1_C*SUM1_C(1:3,1:3) + X1_D*SUM1_D(1:3,1:3)
     &               + X2_A*SUM2_A(1:3,1:3) - X2_B*SUM2_B(1:3,1:3)
     &               - X2_C*SUM2_C(1:3,1:3) + X2_D*SUM2_D(1:3,1:3)
C
               END DO LOOP_IOMEGA
C-----------------------------------------------------------------------
C
            END DO ! JQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      END DO ! ILOOP_OBSE
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ, JQ = ',2I3,/,10X,21('='),/)
99002 FORMAT (//,1X,79('*'),/,34X,'<',A,'>',/,1X,79('*'),//,10X,A,/)
99003 FORMAT (/,40('*'),/,40('*'))
99004 FORMAT (3('(',F14.6,',',F12.6,')'))
      END
