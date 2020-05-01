C*==wavfun_wgt_rad_flm_sra.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_WGT_RAD_FLM_SRA(IT,ZGA,JGA,ZGB,JGB)
C   ********************************************************************
C   *                                                                  *
C   *  weight the wave function ZA with the shape function f_LM        *
C   *                                                                  *
C   *           ZB_L1  = Sum{L2,Lf} f_Lf ZA_L2 C(L1,Lf,L2)  R2DRDI     *
C   *                                                                  *
C   *  and the radial factor R2DRDI for integration                    *
C   *                                                                  *
C   *  NOTE: weighting of the wave functions with the shape functions  *
C   *        must not change the coupling scheme, i.e. the settings    *
C   *        for NCPLWF and IKMCPLWF are unchanged                     *
C   *                                                                  *
C   *  scalar relativistic version                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NCPLWF,NLM,NSPIN,A_RGNT
      USE MOD_RMESH,ONLY:NRMAX,NSF,LMISF,FLMSF,JRCRI,JRMT,KLMSF,R2DRDI,
     &    FULLPOT,JRWS
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF,NLMFP
      IMPLICIT NONE
C*--WAVFUN_WGT_RAD_FLM_SRA25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
      COMPLEX*16 JGA(NRMAX,NCPLWFMAX,NKMMAX),JGB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 FLM_RGNT,RGNT
      INTEGER IA,IB,ILMA,ILMB,IM,IR,IRCRIT,IRSF,IRTOP_SPHERE,IS,ISF,
     &        JLMS,LMSF,LMSOFF
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRCRIT = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRCRIT = JRWS(IM)
         NSF(IM) = 1
         LMISF(1,IM) = 1
         KLMSF(1,IM) = 1
      END IF
C
      ZGB(:,:,:) = C0
      JGB(:,:,:) = C0
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO IS = 1,NSPIN
C
         LMSOFF = NLM*(IS-1)
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         DO JLMS = LMSOFF + 1,LMSOFF + NLM
C
C-----------------------------------------------------------------------
            DO IB = 1,NCPLWF(JLMS)
               ILMB = IKMCPLWF(IB,JLMS)
C
               DO IA = 1,NCPLWF(JLMS)
                  ILMA = IKMCPLWF(IA,JLMS)
C
                  DO ISF = 1,NSF(IM)
                     LMSF = LMISF(ISF,IM)
                     IF ( LMSF.LE.NLMFP .AND. KLMSF(LMSF,IM).NE.0 ) THEN
C
                        RGNT = DREAL(A_RGNT(ILMB,ILMA,LMSF))
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                        IF ( ABS(RGNT).GT.1D-6 ) THEN
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
                           IF ( LMSF.EQ.1 ) THEN
C
C                             FLM_RGNT = F00SF*RGNT = 1
C
                              DO IR = 1,IRTOP_SPHERE
C
                                 ZGB(IR,IB,JLMS) = ZGB(IR,IB,JLMS)
     &                              + ZGA(IR,IA,JLMS)*R2DRDI(IR,IM)
C
                                 JGB(IR,IB,JLMS) = JGB(IR,IB,JLMS)
     &                              + JGA(IR,IA,JLMS)*R2DRDI(IR,IM)
C
                              END DO
                           END IF
C
C-----------------------------------------------------------------------
C                         interstitial regime
C-----------------------------------------------------------------------
C
                           DO IR = IRTOP_SPHERE + 1,IRCRIT
                              IRSF = IR - IRTOP_SPHERE
                              FLM_RGNT = FLMSF(IRSF,ISF,IM)
     &                           *RGNT*R2DRDI(IR,IM)
C
                              ZGB(IR,IB,JLMS) = ZGB(IR,IB,JLMS)
     &                           + ZGA(IR,IA,JLMS)*FLM_RGNT
C
                              JGB(IR,IB,JLMS) = JGB(IR,IB,JLMS)
     &                           + JGA(IR,IA,JLMS)*FLM_RGNT
C
                           END DO
C
                        END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                     END IF
                  END DO
C                                                                     SF
C-----------------------------------------------------------------------
C
               END DO
            END DO
C                                                              ILMA ILMB
C-----------------------------------------------------------------------
C
         END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      END
C*==wavfun_wgt_rad_flm_rel.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_WGT_RAD_FLM_REL(IT,ZGA,ZFA,JGA,JFA,ZGB,ZFB,JGB,
     &                                  JFB,NCPLWF,IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *  weight the wave function ZA with the shape function f_LM        *
C   *                                                                  *
C   *           ZB_L1  = Sum{L2,Lf} f_Lf ZA_L2 C(L1,Lf,L2)  R2DRDI     *
C   *                                                                  *
C   *  and the radial factor R2DRDI for integration                    *
C   *                                                                  *
C   *  NOTE: weighting of the wave functions with the shape functions  *
C   *        must not change the coupling scheme, i.e. the settings    *
C   *        for NCPLWF and IKMCPLWF are unchanged                     *
C   *                                                                  *
C   *  relativistic version                                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,AG_RGNT,AF_RGNT
      USE MOD_RMESH,ONLY:NRMAX,JRCRI,JRMT,NSF,LMISF,R2DRDI,FLMSF,KLMSF,
     &    FULLPOT,JRWS
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,NLMFP
      IMPLICIT NONE
C*--WAVFUN_WGT_RAD_FLM_REL175
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
      COMPLEX*16 JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),ZFA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 AME_F,AME_G,FLM_AF_RGNT,FLM_AG_RGNT
      REAL*8 FLM_RWGT
      INTEGER IA,IB,IKMA,IKMB,IM,IR,IRCRIT,IRSF,IRTOP_SPHERE,ISF,LAM,
     &        LMSF
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRCRIT = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRCRIT = JRWS(IM)
         NSF(IM) = 1
         LMISF(1,IM) = 1
         KLMSF(1,IM) = 1
      END IF
C
      JFB(:,:,:) = C0
      JGB(:,:,:) = C0
      ZFB(:,:,:) = C0
      ZGB(:,:,:) = C0
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                                    LAM
      LOOP_LAM:DO LAM = 1,NKM
C
C-----------------------------------------------------------------------
         LOOP_IB:DO IB = 1,NCPLWF(LAM)
            IKMB = IKMCPLWF(IB,LAM)
C
            LOOP_IA:DO IA = 1,NCPLWF(LAM)
               IKMA = IKMCPLWF(IA,LAM)
C
               LOOP_ISF:DO ISF = 1,NSF(IM)
C
                  LMSF = LMISF(ISF,IM)
C
                  IF ( LMSF.GT.NLMFP ) CYCLE LOOP_ISF
                  IF ( KLMSF(LMSF,IM).NE.1 ) CYCLE LOOP_ISF
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                  AME_G = AG_RGNT(IKMB,IKMA,LMSF)
                  AME_F = AF_RGNT(IKMB,IKMA,LMSF)
C                 AME_F = AF_RGNT(IKMB,IKMA,LMSF) = AME_G
C
                  IF ( ABS(AME_G).LT.1D-6 ) CYCLE LOOP_ISF
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
                  IF ( LMSF.EQ.1 ) THEN
C
C                    FLM_A_RGNT = F00SF * A_RGNT = 1
C
                     DO IR = 1,IRTOP_SPHERE
                        ZGB(IR,IB,LAM) = ZGB(IR,IB,LAM) + ZGA(IR,IA,LAM)
     &                     *R2DRDI(IR,IM)
C
                        ZFB(IR,IB,LAM) = ZFB(IR,IB,LAM) + ZFA(IR,IA,LAM)
     &                     *R2DRDI(IR,IM)
C
                        JGB(IR,IB,LAM) = JGB(IR,IB,LAM) + JGA(IR,IA,LAM)
     &                     *R2DRDI(IR,IM)
C
                        JFB(IR,IB,LAM) = JFB(IR,IB,LAM) + JFA(IR,IA,LAM)
     &                     *R2DRDI(IR,IM)
                     END DO
C
                  END IF
C
C-----------------------------------------------------------------------
C                         interstitial regime
C-----------------------------------------------------------------------
C
                  DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                     IRSF = IR - IRTOP_SPHERE
C
                     FLM_RWGT = FLMSF(IRSF,ISF,IM)*R2DRDI(IR,IM)
C
                     FLM_AG_RGNT = FLM_RWGT*AME_G
                     FLM_AF_RGNT = FLM_RWGT*AME_F
C
                     ZGB(IR,IB,LAM) = ZGB(IR,IB,LAM) + ZGA(IR,IA,LAM)
     &                                *FLM_AG_RGNT
C
                     ZFB(IR,IB,LAM) = ZFB(IR,IB,LAM) + ZFA(IR,IA,LAM)
     &                                *FLM_AF_RGNT
C
                     JGB(IR,IB,LAM) = JGB(IR,IB,LAM) + JGA(IR,IA,LAM)
     &                                *FLM_AG_RGNT
C
                     JFB(IR,IB,LAM) = JFB(IR,IB,LAM) + JFA(IR,IA,LAM)
     &                                *FLM_AF_RGNT
C
                  END DO
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
               END DO LOOP_ISF
C-----------------------------------------------------------------------
C
            END DO LOOP_IA
         END DO LOOP_IB
C-----------------------------------------------------------------------
C
      END DO LOOP_LAM
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END
C*==wavfun_rad_wgt_rel.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WAVFUN_RAD_WGT_REL(IT,ZGA,ZFA,JGA,JFA,NCPLWF,KEY)
C   ********************************************************************
C   *                                                                  *
C   *  weight the wave functions Z and J with  W  for integration      *
C   *                                                                  *
C   *  KEY='R':        W =  R2DRDI                                     *
C   *                                                                  *
C   *  KEY='I':        W =  R2DRDI_W_RADINT = R2DRDI * W_RADINT        *
C   *                                                                  *
C   *  the INTEGRATION scheme is fixed when setting W_RADINT           *
C   *  so far the Simpson scheme is used throughout                    *
C   *                                                                  *
C   *  relativistic version                                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM
      USE MOD_RMESH,ONLY:NRMAX,JRCRI,NSF,LMISF,R2DRDI,KLMSF,FULLPOT,
     &    JRWS,R2DRDI_W_RADINT
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX
      IMPLICIT NONE
C*--WAVFUN_RAD_WGT_REL338
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_RAD_WGT_REL')
C
C Dummy arguments
C
      INTEGER IT
      CHARACTER*1 KEY
      COMPLEX*16 JFA(NRMAX,NCPLWFMAX,NKMMAX),JGA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX)
      INTEGER NCPLWF(NKMMAX)
C
C Local variables
C
      INTEGER IA,IM,IR,IRCRIT,LAM
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRCRIT = JRCRI(IM)
      ELSE
         IRCRIT = JRWS(IM)
         NSF(IM) = 1
         LMISF(1,IM) = 1
         KLMSF(1,IM) = 1
      END IF
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                                    LAM
      LOOP_LAM:DO LAM = 1,NKM
C
C-----------------------------------------------------------------------
         LOOP_IA:DO IA = 1,NCPLWF(LAM)
C
            IF ( KEY.EQ.'R' ) THEN
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
               DO IR = 1,IRCRIT
C
                  ZGA(IR,IA,LAM) = ZGA(IR,IA,LAM)*R2DRDI(IR,IM)
                  ZFA(IR,IA,LAM) = ZFA(IR,IA,LAM)*R2DRDI(IR,IM)
                  JGA(IR,IA,LAM) = JGA(IR,IA,LAM)*R2DRDI(IR,IM)
                  JFA(IR,IA,LAM) = JFA(IR,IA,LAM)*R2DRDI(IR,IM)
C
               END DO
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
            ELSE IF ( KEY.EQ.'I' ) THEN
CIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
               DO IR = 1,IRCRIT
C
                  ZGA(IR,IA,LAM) = ZGA(IR,IA,LAM)*R2DRDI_W_RADINT(IR,IM)
                  ZFA(IR,IA,LAM) = ZFA(IR,IA,LAM)*R2DRDI_W_RADINT(IR,IM)
                  JGA(IR,IA,LAM) = JGA(IR,IA,LAM)*R2DRDI_W_RADINT(IR,IM)
                  JFA(IR,IA,LAM) = JFA(IR,IA,LAM)*R2DRDI_W_RADINT(IR,IM)
C
               END DO
C
CIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
            ELSE
               CALL STOP_MESSAGE(ROUTINE,'KEY NOT FOUND')
            END IF
C
         END DO LOOP_IA
C-----------------------------------------------------------------------
C
      END DO LOOP_LAM
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END
