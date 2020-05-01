C*==sack0.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK0(IQ,INT0,SLN,ALN,FLN,R13LN,R23SQ,M2SUM,R23PW,
     &                 KFP_L1,WINT,AFRLN,USUM,VR23,WR23PW,IR_START_MIN,
     &                 JR_START,JRCRIT,JM,LMAX_FP,LMAX_RHOTIL,NLAM)
C   ********************************************************************
C   *                                                                  *
C   *                evaluate integral   INT0                          *
C   *                                                                  *
C   *  see: Ch. Zecha, Thesis                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:KFP_LMQ
      USE MOD_TYPES,ONLY:NLMFPMAX
      USE MOD_RMESH,ONLY:NRMAX,NPAN,JRCUT
      IMPLICIT NONE
C*--SACK017
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,IR_START_MIN,JM,JRCRIT,JR_START,LMAX_FP,LMAX_RHOTIL,
     &        NLAM
      REAL*8 R13LN
      REAL*8 AFRLN(0:NLAM),ALN(-NLAM:NLAM),FLN(0:NLAM),INT0(NLMFPMAX),
     &       M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,NLMFPMAX)
     &       ,R23PW(IR_START_MIN:NRMAX,-(LMAX_RHOTIL+1):+(LMAX_RHOTIL+1)
     &       ),R23SQ(0:NLAM,NRMAX),USUM(0:NLAM),VR23(NRMAX),WINT(NRMAX),
     &       WR23PW(NRMAX)
      INTEGER KFP_L1(0:LMAX_FP),SLN(-NLAM:NLAM)
C
C Local variables
C
      REAL*8 DDOT
      INTEGER IPAN,JR,JR1,JR2,L1,L2,L3,L3O,L3U,LAMBD3P1,LAMBDAP1,
     &        LAMBDAP2,LM1,LM10,M1,M1LAM3P1,U,UO
      REAL*8 VSUM
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C           weight for radial integration via trapez rule
C-----------------------------------------------------------------------
C
      WINT(JR_START) = 0.0D0
      DO IPAN = NPAN(JM),1, - 1
C
         JR1 = JRCUT(IPAN-1,JM) + 1
         JR2 = JRCUT(IPAN,JM)
C
         IF ( JR2.EQ.JR_START ) EXIT
C
         IF ( JR1.LT.JR_START ) JR1 = JR_START
C
         DO JR = JR1 + 1,JR2 - 1
            WINT(JR) = 1.0D0
         END DO
C
         WINT(JR1) = 0.5D0
         WINT(JR2) = 0.5D0
C
         IF ( JR1.EQ.JR_START ) EXIT
C
      END DO
C
C-----------------------------------------------------------------------
C                    evaluate integral   INT0
C-----------------------------------------------------------------------
C
      INT0(:) = 0D0
C
      DO L1 = 0,LMAX_FP
         IF ( KFP_L1(L1).NE.0 ) THEN
C
            LM10 = L1*L1 + L1 + 1
C
            UO = L1 + LMAX_FP + 1
            DO U = 0,UO
               AFRLN(U) = ALN(L1-U) + FLN(U) + DBLE(U)*R13LN
            END DO
C
            DO L2 = 0,LMAX_RHOTIL
               L3U = ABS(L1-L2)
               L3O = L1 + L2
C
               LAMBDAP1 = MAX(L1,L2)
               M1LAM3P1 = (-1)**MIN(L1,L2)
C
               DO JR = JR_START,JRCRIT
                  WR23PW(JR) = WINT(JR)*R23PW(JR,-L2-1)
               END DO
C
               DO L3 = L3U,L3O,2
C
C  LAMBDAP1  =  1 + (l1+l2+l3)/2
                  LAMBDAP1 = LAMBDAP1 + 1
C
C  LAMBD3P1  =  1 + (l1+l2-l3)/2
                  LAMBD3P1 = LAMBDAP1 - L3
C
C  M1LAM3P1  =  (-1)**(1+(l1+l2-l3)/2)
                  M1LAM3P1 = -M1LAM3P1
C
                  CALL SACK_USUM(USUM,SLN,ALN,FLN,AFRLN,LAMBDAP1,
     &                           LAMBD3P1,M1LAM3P1,L1,L2,NLAM)
C
                  LAMBDAP2 = LAMBDAP1 + 1
C
                  DO JR = JR_START,JRCRIT
                     VSUM = DDOT(LAMBDAP2,USUM(0),1,R23SQ(0,JR),1)
                     VR23(JR) = VSUM*WR23PW(JR)
                  END DO
C
                  DO M1 = -L1,L1
                     LM1 = LM10 + M1
C
                     IF ( KFP_LMQ(LM1,IQ).NE.0 ) INT0(LM1) = INT0(LM1)
     &                    + DDOT(JRCRIT-JR_START+1,VR23(JR_START),1,
     &                    M2SUM(JR_START,L3,L2,LM1),1)
                  END DO
C
               END DO
C
            END DO
C
         END IF
      END DO
C
      END
C*==sack1.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK1(IQ,INT1R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                 JR_START,JRCRIT,JM,LMAX_FP,LMAX_RHOTIL,NL_RHOTIL,
     &                 NLAM)
C   ********************************************************************
C   *                                                                  *
C   *                evaluate integral   INT1                          *
C   *                                                                  *
C   *  see: Ch. Zecha, Thesis                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:KFP_LMQ
      USE MOD_TYPES,ONLY:NLFPMAX,NLMFPMAX
      USE MOD_RMESH,ONLY:NRMAX,NPAN,JRCUT
      IMPLICIT NONE
C*--SACK1158
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,IR_START_MIN,JM,JRCRIT,JR_START,LMAX_FP,LMAX_RHOTIL,
     &        NLAM,NL_RHOTIL
      REAL*8 DFAC(0:NL_RHOTIL,0:NL_RHOTIL),
     &       INT1R(IR_START_MIN:NRMAX,NLMFPMAX),
     &       M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,NLMFPMAX)
     &       ,R23PW(IR_START_MIN:NRMAX,-(LMAX_RHOTIL+1):+(LMAX_RHOTIL+1)
     &       )
      INTEGER KFP_L1(0:LMAX_FP)
C
C Local variables
C
      REAL*8 DF,FAKTOR,L2SUM(:,:)
      INTEGER JR,L1,L2,L3,LM1,LM10,M1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE L2SUM
C
      ALLOCATE (L2SUM(IR_START_MIN:NRMAX,-(NLFPMAX-1):(NLFPMAX-1)))
C
C*** End of declarations rewritten by SPAG
C
      DO L1 = 0,LMAX_FP
         IF ( KFP_L1(L1).NE.0 ) THEN
C
            LM10 = L1*L1 + L1 + 1
C
            L2SUM(:,:) = 0D0
C
            DO L2 = 0,L1
               L3 = L1 - L2
C
               DF = 4.0D0*DFAC(L2,L3)*(-1D0)**L2
C
               DO JR = JR_START,JRCRIT
C
C              4 * DFAC(L2,L3) * (-1)**L2 * (r'/r'')**L2
C
                  FAKTOR = DF*R23PW(JR,L2)
C
                  DO M1 = -L1,L1
                     LM1 = LM10 + M1
C
                     IF ( KFP_LMQ(LM1,IQ).NE.0 ) L2SUM(JR,M1)
     &                    = L2SUM(JR,M1) + FAKTOR*M2SUM(JR,L3,L2,LM1)
                  END DO
C
               END DO
C
            END DO
C
            DO M1 = -L1,L1
               LM1 = LM10 + M1
C
               IF ( KFP_LMQ(LM1,IQ).NE.0 )
     &              CALL SACK_RINT(L2SUM(IR_START_MIN,M1),
     &              INT1R(IR_START_MIN,LM1),NPAN(JM),JRCUT(0,JM),
     &              IR_START_MIN,JR_START)
            END DO
C
         END IF
      END DO
C
      END
C*==sack2.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK2(IQ,INT2R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                 JR_START,JRCRIT,JM,LMAX_FP,LMAX_RHOTIL,NL_RHOTIL,
     &                 NLAM)
C   ********************************************************************
C   *                                                                  *
C   *                evaluate integral   INT2                          *
C   *                                                                  *
C   *  see: Ch. Zecha, Thesis                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:KFP_LMQ
      USE MOD_TYPES,ONLY:NLFPMAX,NLMFPMAX
      USE MOD_RMESH,ONLY:NRMAX,NPAN,JRCUT
      IMPLICIT NONE
C*--SACK2255
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,IR_START_MIN,JM,JRCRIT,JR_START,LMAX_FP,LMAX_RHOTIL,
     &        NLAM,NL_RHOTIL
      REAL*8 DFAC(0:NL_RHOTIL,0:NL_RHOTIL),
     &       INT2R(IR_START_MIN:NRMAX,NLMFPMAX),
     &       M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,NLMFPMAX)
     &       ,R23PW(IR_START_MIN:NRMAX,-(LMAX_RHOTIL+1):+(LMAX_RHOTIL+1)
     &       )
      INTEGER KFP_L1(0:LMAX_FP)
C
C Local variables
C
      REAL*8 DF,FAKTOR,L2SUM(:,:)
      INTEGER JR,L1,L2,L3,LM1,LM10,M1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE L2SUM
C
      ALLOCATE (L2SUM(IR_START_MIN:NRMAX,-(NLFPMAX-1):(NLFPMAX-1)))
C
      DO L1 = 0,LMAX_FP
         IF ( KFP_L1(L1).NE.0 ) THEN
C
            LM10 = L1*L1 + L1 + 1
C
            L2SUM(:,:) = 0D0
C
            DO L2 = L1,LMAX_RHOTIL
               L3 = L2 - L1
               DF = 4.0D0*DFAC(L1,L3)
C
               DO JR = JR_START,JRCRIT
C
                  FAKTOR = DF*R23PW(JR,-L2-1)
C
                  DO M1 = -L1,L1
                     LM1 = LM10 + M1
C
                     IF ( KFP_LMQ(LM1,IQ).NE.0 ) L2SUM(JR,M1)
     &                    = L2SUM(JR,M1) + FAKTOR*M2SUM(JR,L3,L2,LM1)
                  END DO
C
               END DO
C
            END DO
C
            DO M1 = -L1,L1
               LM1 = LM10 + M1
C
               IF ( KFP_LMQ(LM1,IQ).NE.0 )
     &              CALL SACK_RINT(L2SUM(IR_START_MIN,M1),
     &              INT2R(IR_START_MIN,LM1),NPAN(JM),JRCUT(0,JM),
     &              IR_START_MIN,JR_START)
            END DO
C
         END IF
      END DO
C
      END
C*==sack3.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK3(IQ,INT3R,DFAC,M2SUM,R23PW,KFP_L1,IR_START_MIN,
     &                 JR_START,JRCRIT,JM,LMAX_FP,LMAX_RHOTIL,NL_RHOTIL,
     &                 NLAM)
C   ********************************************************************
C   *                                                                  *
C   *                evaluate integral   INT3                          *
C   *                                                                  *
C   *  see: Ch. Zecha, Thesis                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:KFP_LMQ
      USE MOD_TYPES,ONLY:NLFPMAX,NLMFPMAX
      USE MOD_RMESH,ONLY:NRMAX,NPAN,JRCUT
      IMPLICIT NONE
C*--SACK3347
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,IR_START_MIN,JM,JRCRIT,JR_START,LMAX_FP,LMAX_RHOTIL,
     &        NLAM,NL_RHOTIL
      REAL*8 DFAC(0:NL_RHOTIL,0:NL_RHOTIL),
     &       INT3R(IR_START_MIN:NRMAX,NLMFPMAX),
     &       M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,NLMFPMAX)
     &       ,R23PW(IR_START_MIN:NRMAX,-(LMAX_RHOTIL+1):+(LMAX_RHOTIL+1)
     &       )
      INTEGER KFP_L1(0:LMAX_FP)
C
C Local variables
C
      REAL*8 DF,FAKTOR,L2SUM(:,:)
      INTEGER JR,L1,L2,L3,LM1,LM10,M1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE L2SUM
C
      ALLOCATE (L2SUM(IR_START_MIN:NRMAX,-(NLFPMAX-1):(NLFPMAX-1)))
C
      DO L1 = 0,LMAX_FP
         IF ( KFP_L1(L1).NE.0 ) THEN
C
            LM10 = L1*L1 + L1 + 1
C
            L2SUM(:,:) = 0D0
C
            DO L2 = 0,LMAX_RHOTIL
               L3 = L1 + L2
               DF = 4.0D0*DFAC(L1,L2)
C
               DO JR = JR_START,JRCRIT
C
                  FAKTOR = DF*R23PW(JR,L2)
C
                  DO M1 = -L1,L1
                     LM1 = LM10 + M1
C
                     IF ( KFP_LMQ(LM1,IQ).NE.0 ) L2SUM(JR,M1)
     &                    = L2SUM(JR,M1) + FAKTOR*M2SUM(JR,L3,L2,LM1)
                  END DO
C
               END DO
C
            END DO
C
            DO M1 = -L1,L1
               LM1 = LM10 + M1
C
               IF ( KFP_LMQ(LM1,IQ).NE.0 )
     &              CALL SACK_RINT(L2SUM(IR_START_MIN,M1),
     &              INT3R(IR_START_MIN,LM1),NPAN(JM),JRCUT(0,JM),
     &              IR_START_MIN,JR_START)
C
            END DO
C
         END IF
      END DO
C
      END
C*==sack_usum.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_USUM(USUM,SLN,ALN,FLN,AFRLN,LAMBDAP1,LAMBD3P1,
     &                     M1LAM3P1,L1,L2,NLAM)
C   ********************************************************************
C   *                                                                  *
C   *  Ch. Zecha   Thesis                                              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SACK_USUM434
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L1,L2,LAMBD3P1,LAMBDAP1,M1LAM3P1,NLAM
      REAL*8 AFRLN(0:NLAM),ALN(-NLAM:NLAM),FLN(0:NLAM),USUM(0:NLAM)
      INTEGER SLN(-NLAM:NLAM)
C
C Local variables
C
      REAL*8 AFLNV,EAF,RSUM,SLNV,VZ
      INTEGER LAMBP1MV,U,UPV,V
C
C*** End of declarations rewritten by SPAG
C
      DO V = 0,LAMBDAP1
C
         LAMBP1MV = LAMBDAP1 - V
         SLNV = SLN(L2-V)
         AFLNV = ALN(L2-V) + FLN(V)
         RSUM = 0.0D0
C
         DO U = 0,LAMBP1MV
            UPV = U + V
C
            VZ = DBLE(M1LAM3P1*SLN(L1-U)*SLNV*SLN(UPV-LAMBD3P1))
C
            EAF = EXP(AFRLN(U)+AFLNV+ALN(UPV-LAMBD3P1)+FLN(LAMBDAP1-UPV)
     &            )
C
            RSUM = RSUM + VZ*EAF
         END DO
C
         USUM(V) = RSUM
      END DO
C
      END
C*==sack_cysum.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_CYSUM(M3SUM,RHAT3,LMRGNT_FPTILLAM,
     &                      NRGNT_FPTILLAM_LM1,RGNT_FPTILLAM,
     &                      NRGNT_FPTILLAM,NLMFP,NLM_RHOTIL,NLAM)
C   ********************************************************************
C   *                                                                  *
C   * calculates SUM Gaunt(L,L',L'') * RealSpherHarm(L'', hat(r)'')    *
C   *            m''                                                   *
C   *                                                                  *
C   *  Ch. Zecha   Thesis                                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:L_LM
      USE MOD_TYPES,ONLY:NLMFPMAX
      IMPLICIT NONE
C*--SACK_CYSUM500
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLAM,NLMFP,NLM_RHOTIL,NRGNT_FPTILLAM
      INTEGER LMRGNT_FPTILLAM(NRGNT_FPTILLAM,3),
     &        NRGNT_FPTILLAM_LM1(0:NLMFPMAX)
      REAL*8 M3SUM(0:(NLAM-1),NLM_RHOTIL,NLMFPMAX),
     &       RGNT_FPTILLAM(NRGNT_FPTILLAM),RHAT3(3)
C
C Local variables
C
      INTEGER J,L3,LM1,LM2,LM3,NLAMSQ
      REAL*8 RYLM(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RYLM
C
      NLAMSQ = NLAM*NLAM
C
      ALLOCATE (RYLM(NLAMSQ))
C
      M3SUM(:,:,:) = 0D0
C
      CALL CALC_RHPLM(RHAT3(1),RHAT3(2),RHAT3(3),RYLM,(NLAM-1),NLAMSQ)
C
      DO LM1 = 1,NLMFP
         DO J = NRGNT_FPTILLAM_LM1(LM1-1) + 1,NRGNT_FPTILLAM_LM1(LM1)
            LM2 = LMRGNT_FPTILLAM(J,2)
            LM3 = LMRGNT_FPTILLAM(J,3)
            L3 = L_LM(LM3)
C
            M3SUM(L3,LM2,LM1) = M3SUM(L3,LM2,LM1) + RGNT_FPTILLAM(J)
     &                          *RYLM(LM3)
         END DO
      END DO
C
      END
C*==sack_mmsum.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_MMSUM(M2SUM,M3SUM,R2RHOTIL,KRHOTIL_LMQ,
     &                      IR_START_MIN,LMAX_FP,LMAX_RHOTIL,NLM_RHOTIL,
     &                      NLAM)
C   ********************************************************************
C   *                                                                  *
C   * calculates SUM n^j_L' (r') * M3SUM                               *
C   *             m'                                                   *
C   *    uses M3SUM from <SACK_CYSUM>                                  *
C   *                                                                  *
C   *  Ch. Zecha   Thesis                                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NLMFPMAX,NLMFP
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--SACK_MMSUM569
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SACK_MMSUM')
C
C Dummy arguments
C
      INTEGER IR_START_MIN,LMAX_FP,LMAX_RHOTIL,NLAM,NLM_RHOTIL
      INTEGER KRHOTIL_LMQ(NLM_RHOTIL)
      REAL*8 M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,NLMFPMAX)
     &       ,M3SUM(0:(NLAM-1),NLM_RHOTIL,NLMFPMAX),
     &       R2RHOTIL(NRMAX,NLM_RHOTIL)
C
C Local variables
C
      REAL*8 FAKTOR
      INTEGER L2,L3,LM1,LM2,M2
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      IF ( NLMFP.NE.(LMAX_FP+1)**2 )
     &      CALL STOP_MESSAGE(ROUTINE,'NLMFP .NE. (LMAX_FP+1)**2')
C
      M2SUM(IR_START_MIN:NRMAX,0:(NLAM-1),0:LMAX_RHOTIL,1:NLMFPMAX)
     &   = 0D0
C
      DO LM1 = 1,NLMFP
C
         DO L3 = 0,(NLAM-1)
            LM2 = 0
C
            DO L2 = 0,LMAX_RHOTIL
               DO M2 = -L2,L2
                  LM2 = LM2 + 1
C
                  IF ( KRHOTIL_LMQ(LM2).NE.0 ) THEN
C
                     FAKTOR = M3SUM(L3,LM2,LM1)
C
                     CALL DAXPY(NRMAX-IR_START_MIN+1,FAKTOR,
     &                          R2RHOTIL(IR_START_MIN,LM2),1,
     &                          M2SUM(IR_START_MIN,L3,L2,LM1),1)
C
                  END IF
C
               END DO
            END DO
         END DO
      END DO
C
      END
C*==sack_fln.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_FLN(SLN,ALN,FLN,NLAM)
C   ********************************************************************
C   *                                                                  *
C   *  Ch. Zecha   Thesis                                              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SACK_FLN657
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLAM
      REAL*8 ALN(-NLAM:NLAM),FLN(0:NLAM)
      INTEGER SLN(-NLAM:NLAM)
C
C Local variables
C
      INTEGER JR
C
C*** End of declarations rewritten by SPAG
C
      SLN(0) = 1
      ALN(0) = 0.0D0
      FLN(0) = 0.0D0
C
      DO JR = 1,NLAM
         SLN(JR) = 1
         ALN(JR) = ALN(JR-1) + LOG(DBLE(-0.5D0+JR))
         SLN(-JR) = -SLN(-JR+1)
         ALN(-JR) = -ALN(JR)
         FLN(JR) = FLN(JR-1) - LOG(DBLE(JR))
      END DO
C
      END
C*==sack_facul.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_FACUL(DFAC,NL_RHOTIL)
C   ********************************************************************
C   *                                                                  *
C   *                            ( 2l + 2l' - 1 )!!                    *
C   *         DFAC(l,l') =   ---------------------------               *
C   *                        ( 2l + 1 )!!  ( 2l' + 1 )!!               *
C   *                                                                  *
C   *  Ch. Zecha   Thesis Eq. (6.18)                                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SACK_FACUL710
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NL_RHOTIL
      REAL*8 DFAC(0:NL_RHOTIL,0:NL_RHOTIL)
C
C Local variables
C
      INTEGER L0,L1
C
C*** End of declarations rewritten by SPAG
C
      DFAC(0:NL_RHOTIL,0:NL_RHOTIL) = 0.0D0
C
      DO L1 = 0,NL_RHOTIL
         DFAC(0,L1) = 1.0D0/DBLE(2*L1+1)
         DO L0 = 1,L1
            DFAC(L0,L1) = DFAC(L0-1,L1)*DBLE(2*(L0+L1)-1)/DBLE(2*L0+1)
         END DO
      END DO
C
      DO L1 = 0,NL_RHOTIL
         DO L0 = L1 + 1,NL_RHOTIL
            DFAC(L0,L1) = DFAC(L1,L0)
         END DO
      END DO
C
      END
C*==sack_rint.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SACK_RINT(F,FINT,NPAN,IRCUT,IR_START_MIN,JR_START)
C   ********************************************************************
C   *                                                                  *
C   *  this subroutine does an inwards integration of a function       *
C   *  with kinks                                                      *
C   *                                                                  *
C   *                           r_crit                                 *
C   *                 FINT(r) = S         f(r') dr'                    *
C   *                           r_start                                *
C   *                                                                  *
C   *  at each kink the integration is restarted.                      *
C   *  Restarting with 3 or 4 point lagrangian integration             *
C   *  yields bad results (peaks located at the kinks).                *
C   *  Therefore an extended trapezoidal rule is used throughout.      *
C   *  Integration is performed inwards from the critcal radius        *
C   *  only down to  JR_START                                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,NPANMAX
      IMPLICIT NONE
C*--SACK_RINT774
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IR_START_MIN,JR_START,NPAN
      REAL*8 F(IR_START_MIN:NRMAX),FINT(IR_START_MIN:NRMAX)
      INTEGER IRCUT(0:NPANMAX)
C
C Local variables
C
      INTEGER IPAN,JR,JR1,JR2
C
C*** End of declarations rewritten by SPAG
C
C---> loop over kinks
C
      DO IPAN = NPAN,1, - 1
C
         JR1 = IRCUT(IPAN)
         JR2 = IRCUT(IPAN-1) + 1
C
         IF ( IPAN.EQ.NPAN ) THEN
C
            FINT(JR1) = 0.0D0
C
         ELSE
C
            FINT(JR1) = FINT(JR1+1)
C
         END IF
C
         IF ( JR1.EQ.JR_START ) RETURN
C
         IF ( JR2.LT.JR_START ) JR2 = JR_START
C
         DO JR = JR1 - 1,JR2, - 1
            FINT(JR) = FINT(JR+1) + 0.5D0*(F(JR+1)+F(JR))
         END DO
C
         IF ( JR2.EQ.JR_START ) RETURN
C
      END DO
C
      END
