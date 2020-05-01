C*==linresp_density_sra.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_DENSITY_SRA(IT,IRTOP,WE,TMATA,TMATB,
     &                               THZT_OFF_TP,THZT_DIA_P,THXT_OFF_TP,
     &                               THXT_DIA_P,ZGA,JGA,ZGB,JGB)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <LINRESP_DENSITY>                       *
C   *                                                                  *
C   *   dealing with SCALAR RELATIVISTIC case                          *
C   *                                                                  *
C   *  NOTE: for IREL <= 2 (non- and scalar relativistic) the          *
C   *        minor component has no meaning and is not written/read    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T,HAX_ZZ_T,
     &    HBX_JZ_T,HCX_ZJ_T,HDX_JJ_T,AMEG_LMOP,NPERT,NOBSE,RHO2NS_GG,
     &    IOPER_OBSE,CHI_TO
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NCPLWFMAX,KLMFP,IKMCPLWF,NT,NTMAX
      USE MOD_ANGMOM,ONLY:NLM,L_LM,NSPIN,NKMMAX,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--LINRESP_DENSITY_SRA25
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOLTAU
      PARAMETER (TOLTAU=1D-12)
C
C Dummy arguments
C
      INTEGER IRTOP,IT
      COMPLEX*16 WE
      COMPLEX*16 JGA(NRMAX,NCPLWFMAX,NKMMAX),JGB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,THXT_DIA_P(NKMMAX,NKMMAX,NPERT),
     &           THXT_OFF_TP(NKMMAX,NKMMAX,NTMAX,NPERT),
     &           THZT_DIA_P(NKMMAX,NKMMAX,NPERT),
     &           THZT_OFF_TP(NKMMAX,NKMMAX,NTMAX,NPERT),
     &           TMATA(NKMMAX,NKMMAX),TMATB(NKMMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 AME_G,RHO_GG,SPNWGT,WGG(:)
      COMPLEX*16 HX,HZ,JGBJGA,JGBZGA,SUM_HX,SUM_HZ,THTJJ(:),THTJZ(:),
     &           THTZJ(:),THTZZ(:),TMAT12,TMAT34,WGHTGG,ZGBJGA,ZGBZGA
      INTEGER I1,I4,ILM1,ILM4,IOBSE,IOPER,IPERT,IR,IS,JLMS1,JLMS4,JT,K1,
     &        K2,K3,K4,L1,L3,L4,LM3,LMSOFF,M3
      LOGICAL JLMS1_EQ_JLMS4,TAU_NE_0
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*40 ROUTINE
C      PARAMETER (ROUTINE='LINRESP_DENSITY_SRA')
C
      ALLOCATABLE WGG,THTZJ,THTZZ,THTJZ,THTJJ
C
      ALLOCATE (WGG(NRMAX))
      ALLOCATE (THTZJ(NRMAX),THTZZ(NRMAX),THTJZ(NRMAX),THTJJ(NRMAX))
C
      WGHTGG = WE/PI
C
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                   loop over the various perturbation terms
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IPERT = 1,NPERT
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         DO IS = 1,NSPIN
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
            LMSOFF = NLM*(IS-1)
C
C=======================================================================
C
            DO JLMS4 = LMSOFF + 1,LMSOFF + NLM
               K4 = JLMS4
               DO JLMS1 = LMSOFF + 1,LMSOFF + NLM
                  K1 = JLMS1
C
                  JLMS1_EQ_JLMS4 = JLMS1.EQ.JLMS4
C
                  TAU_NE_0 = ABS(DREAL(TMATB(JLMS1,JLMS4)))
     &                       .GT.TOLTAU .OR. 
     &                       ABS(DIMAG(TMATB(JLMS1,JLMS4))).GT.TOLTAU
C
                  IF ( JLMS1_EQ_JLMS4 .OR. TAU_NE_0 ) THEN
C
C-----------------------------------------------------------------------
C      set up the radial weights for the radial wave functions
C-----------------------------------------------------------------------
                     SUM_HZ = THZT_DIA_P(JLMS1,JLMS4,IPERT)
                     SUM_HX = THXT_DIA_P(JLMS1,JLMS4,IPERT)
                     DO JT = 1,NT
                        SUM_HZ = SUM_HZ + THZT_OFF_TP(JLMS1,JLMS4,JT,
     &                           IPERT)
                        SUM_HX = SUM_HX + THXT_OFF_TP(JLMS1,JLMS4,JT,
     &                           IPERT)
                     END DO
C
                     DO IR = 1,IRTOP
C
                        HZ = SUM_HZ + HDZ_JJ_T(IR,JLMS1,JLMS4,IT,IPERT)
C
                        HX = SUM_HX + HDX_JJ_T(IR,JLMS1,JLMS4,IT,IPERT)
C
                        THTZZ(IR) = HZ + HX*CHI_TO(IT,IPERT)
C
                        THTZJ(IR) = 0D0
C
                        THTJZ(IR) = 0D0
C
                        THTJJ(IR) = HAZ_ZZ_T(IR,JLMS1,JLMS4,IT,IPERT)
     &                              + HAX_ZZ_T(IR,JLMS1,JLMS4,IT,IPERT)
     &                              *CHI_TO(IT,IPERT)
C
                     END DO
C
                     DO K3 = LMSOFF + 1,LMSOFF + NLM
C
                        TMAT34 = TMATB(K3,K4)
C
                        IF ( ABS(DREAL(TMAT34)).GT.TOLTAU .OR. 
     &                       ABS(DIMAG(TMAT34)).GT.TOLTAU ) THEN
C
                           DO IR = 1,IRTOP
C
                              HZ = HAZ_ZZ_T(IR,K1,K3,IT,IPERT)
                              HX = HAX_ZZ_T(IR,K1,K3,IT,IPERT)
C
                              THTJZ(IR) = THTJZ(IR)
     &                           - (HZ+HX*CHI_TO(IT,IPERT))*TMAT34
C
                              HZ = HBZ_JZ_T(IR,K1,K3,IT,IPERT)
                              HX = HBX_JZ_T(IR,K1,K3,IT,IPERT)
C
                              THTZZ(IR) = THTZZ(IR)
     &                           - (HZ+HX*CHI_TO(IT,IPERT))*TMAT34
C
                           END DO
C
                        END IF
                     END DO
C
                     DO K2 = LMSOFF + 1,LMSOFF + NLM
C
                        TMAT12 = TMATA(K1,K2)
C
                        IF ( ABS(DREAL(TMAT12)).GT.TOLTAU .OR. 
     &                       ABS(DIMAG(TMAT12)).GT.TOLTAU ) THEN
C
                           DO IR = 1,IRTOP
C
                              HZ = HAZ_ZZ_T(IR,K2,K4,IT,IPERT)
                              HX = HAX_ZZ_T(IR,K2,K4,IT,IPERT)
C
                              THTZJ(IR) = THTZJ(IR)
     &                           - TMAT12*(HZ+HX*CHI_TO(IT,IPERT))
C
                              HZ = HCZ_ZJ_T(IR,K2,K4,IT,IPERT)
                              HX = HCX_ZJ_T(IR,K2,K4,IT,IPERT)
C
                              THTZZ(IR) = THTZZ(IR)
     &                           - TMAT12*(HZ+HX*CHI_TO(IT,IPERT))
C
                           END DO
C
                        END IF
                     END DO
C
C-----------------------------------------------------------------------
C              include the factor WE/PI for energy integration
C-----------------------------------------------------------------------
C
                     DO IR = 1,IRTOP
                        THTZZ(IR) = WGHTGG*THTZZ(IR)
                        THTZJ(IR) = WGHTGG*THTZJ(IR)
                        THTJZ(IR) = WGHTGG*THTJZ(IR)
                        THTJJ(IR) = WGHTGG*THTJJ(IR)
                     END DO
C
C-----------------------------------------------------------------------
C
                     DO I4 = 1,NCPLWF(JLMS4)
                        ILM4 = IKMCPLWF(I4,JLMS4)
                        L4 = L_LM(ILM4)
C
                        DO I1 = 1,NCPLWF(JLMS1)
                           ILM1 = IKMCPLWF(I1,JLMS1)
                           L1 = L_LM(ILM1)
C
C------------------------------------------------------------------------
C    product of radial wave functions and weights ZB(4) * ZA(1) * WAB(r)
C    all terms may contribute to the induced densities
C    NOTE: the angular matrix element to be added is real
C          for the scalar relativistic case
C------------------------------------------------------------------------
C
                           DO IR = 1,IRTOP
                              ZGBZGA = ZGB(IR,I4,JLMS4)*ZGA(IR,I1,JLMS1)
                              ZGBJGA = ZGB(IR,I4,JLMS4)*JGA(IR,I1,JLMS1)
                              JGBZGA = JGB(IR,I4,JLMS4)*ZGA(IR,I1,JLMS1)
                              JGBJGA = JGB(IR,I4,JLMS4)*JGA(IR,I1,JLMS1)
C
                              WGG(IR) = DIMAG(ZGBZGA*THTZZ(IR))
     &                                  + DIMAG(ZGBJGA*THTJZ(IR))
     &                                  + DIMAG(JGBZGA*THTZJ(IR))
     &                                  + DIMAG(JGBJGA*THTJJ(IR))
                           END DO
C
                           DO L3 = ABS(L4-L1),L1 + L4,2
                              LM3 = L3*L3
                              DO M3 = -L3,L3
                                 LM3 = LM3 + 1
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
                                 DO IOBSE = 1,NOBSE
C
                                    IOPER = IOPER_OBSE(IOBSE)
C
C-----------------------------------------------------------------------
C      check angular matrix element  AME_OP  and
C      KPERT flag for contribution LM to observable   IOBSE
C-----------------------------------------------------------------------
C
                                    AME_G = DREAL
     &                                 (AMEG_LMOP(ILM4,ILM1,LM3,IOPER))
                                    IF ( IOBSE.EQ.2 )
     &                                 AME_G = SPNWGT*AME_G
C
                                    IF ( ABS(AME_G).GT.1D-6 .AND. 
     &                                 KLMFP(LM3,IT).NE.0 ) THEN
C
C-----------------------------------------------------------------------
C                          induced densities
C-----------------------------------------------------------------------
C
                                       DO IR = 1,IRTOP
C
                                         RHO_GG = AME_G*WGG(IR)
C
                                         RHO2NS_GG(IR,LM3,IT,IOBSE,
     &                                      IPERT)
     &                                      = RHO2NS_GG(IR,LM3,IT,IOBSE,
     &                                      IPERT) + RHO_GG
C
                                       END DO
C
                                    END IF
C
                                 END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
                              END DO
                           END DO
C
                        END DO
                     END DO
C-----------------------------------------------------------------------
C
                  END IF
               END DO
            END DO
C ======================================================================
C
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      END DO
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      END
C*==linresp_density_rel.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_DENSITY_REL(IT,IRTOP,WE,TMATA,TMATB,
     &                               THZT_OFF_TP,THZT_DIA_P,THXT_OFF_TP,
     &                               THXT_DIA_P,ZGA,ZFA,ZGB,ZFB,JGA,JFA,
     &                               JGB,JFB)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <LINRESP_DENSITY>                       *
C   *                                                                  *
C   *   dealing with FULLY RELATIVISTIC case                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T,HAX_ZZ_T,
     &    HBX_JZ_T,HCX_ZJ_T,HDX_JJ_T,AMEG_LMOP,AMEF_LMOP,NPERT,NOBSE,
     &    RHO2NS_GG,IOPER_OBSE,CHI_TO
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NCPLWFMAX,KLMFP,IKMCPLWF,NT,NTMAX
      USE MOD_ANGMOM,ONLY:NKM,L_IKM,NKMMAX,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--LINRESP_DENSITY_REL315
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOLTAU
      PARAMETER (TOLTAU=1D-12)
C
C Dummy arguments
C
      INTEGER IRTOP,IT
      COMPLEX*16 WE
      COMPLEX*16 JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),
     &           THXT_DIA_P(NKMMAX,NKMMAX,NPERT),
     &           THXT_OFF_TP(NKMMAX,NKMMAX,NTMAX,NPERT),
     &           THZT_DIA_P(NKMMAX,NKMMAX,NPERT),
     &           THZT_OFF_TP(NKMMAX,NKMMAX,NTMAX,NPERT),
     &           TMATA(NKMMAX,NKMMAX),TMATB(NKMMAX,NKMMAX),
     &           ZFA(NRMAX,NCPLWFMAX,NKMMAX),ZFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 AME_F,AME_G,HX,HZ,JFBJFA,JFBZFA,JGBJGA,JGBZGA,SUM_HX,
     &           SUM_HZ,THTJJ(:),THTJZ(:),THTZJ(:),THTZZ(:),TMAT12,
     &           TMAT34,WFF(:),WGG(:),WGHTGG,ZFBJFA,ZFBZFA,ZGBJGA,ZGBZGA
      INTEGER I1,I4,IKM1,IKM4,IOBSE,IOPER,IPERT,IR,JT,K1,K2,K3,K4,L1,L3,
     &        L4,LAMA1,LAMB4,LM3,M3
      LOGICAL LAMA1_EQ_LAMB4,TAU_NE_0
      REAL*8 RHO_FF,RHO_GG
C
C*** End of declarations rewritten by SPAG
C
C      CHARACTER*40 ROUTINE
C      PARAMETER (ROUTINE='LINRESP_DENSITY_REL')
C
      ALLOCATABLE WGG,WFF,THTZJ,THTZZ,THTJZ,THTJJ
C
      ALLOCATE (WGG(NRMAX),WFF(NRMAX))
      ALLOCATE (THTZJ(NRMAX),THTZZ(NRMAX),THTJZ(NRMAX),THTJJ(NRMAX))
C
      WGHTGG = WE/PI
C
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                   loop over the various perturbation terms
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IPERT = 1,NPERT
C
C=======================================================================
C
         DO LAMB4 = 1,NKM
            K4 = LAMB4
            DO LAMA1 = 1,NKM
               K1 = LAMA1
C
               LAMA1_EQ_LAMB4 = LAMA1.EQ.LAMB4
C
               TAU_NE_0 = ABS(DREAL(TMATB(LAMA1,LAMB4))).GT.TOLTAU .OR. 
     &                    ABS(DIMAG(TMATB(LAMA1,LAMB4))).GT.TOLTAU
C
               IF ( LAMA1_EQ_LAMB4 .OR. TAU_NE_0 ) THEN
C
C-----------------------------------------------------------------------
C      set up the radial weights for the radial wave functions
C-----------------------------------------------------------------------
                  SUM_HZ = THZT_DIA_P(LAMA1,LAMB4,IPERT)
                  SUM_HX = THXT_DIA_P(LAMA1,LAMB4,IPERT)
                  DO JT = 1,NT
                     SUM_HZ = SUM_HZ + THZT_OFF_TP(LAMA1,LAMB4,JT,IPERT)
                     SUM_HX = SUM_HX + THXT_OFF_TP(LAMA1,LAMB4,JT,IPERT)
                  END DO
C
                  DO IR = 1,IRTOP
C
                     HZ = SUM_HZ + HDZ_JJ_T(IR,LAMA1,LAMB4,IT,IPERT)
C
                     HX = SUM_HX + HDX_JJ_T(IR,LAMA1,LAMB4,IT,IPERT)
C
                     THTZZ(IR) = HZ + HX*CHI_TO(IT,IPERT)
C
                     THTZJ(IR) = 0D0
C
                     THTJZ(IR) = 0D0
C
                     THTJJ(IR) = HAZ_ZZ_T(IR,LAMA1,LAMB4,IT,IPERT)
     &                           + HAX_ZZ_T(IR,LAMA1,LAMB4,IT,IPERT)
     &                           *CHI_TO(IT,IPERT)
C
                  END DO
C
                  DO K3 = 1,NKM
C
                     TMAT34 = TMATB(K3,K4)
C
                     IF ( ABS(DREAL(TMAT34)).GT.TOLTAU .OR. 
     &                    ABS(DIMAG(TMAT34)).GT.TOLTAU ) THEN
C
                        DO IR = 1,IRTOP
C
                           HZ = HAZ_ZZ_T(IR,K1,K3,IT,IPERT)
                           HX = HAX_ZZ_T(IR,K1,K3,IT,IPERT)
C
                           THTJZ(IR) = THTJZ(IR)
     &                                 - (HZ+HX*CHI_TO(IT,IPERT))*TMAT34
C
                           HZ = HBZ_JZ_T(IR,K1,K3,IT,IPERT)
                           HX = HBX_JZ_T(IR,K1,K3,IT,IPERT)
C
                           THTZZ(IR) = THTZZ(IR)
     &                                 - (HZ+HX*CHI_TO(IT,IPERT))*TMAT34
C
                        END DO
C
                     END IF
                  END DO
C
                  DO K2 = 1,NKM
C
                     TMAT12 = TMATA(K1,K2)
C
                     IF ( ABS(DREAL(TMAT12)).GT.TOLTAU .OR. 
     &                    ABS(DIMAG(TMAT12)).GT.TOLTAU ) THEN
C
                        DO IR = 1,IRTOP
C
                           HZ = HAZ_ZZ_T(IR,K2,K4,IT,IPERT)
                           HX = HAX_ZZ_T(IR,K2,K4,IT,IPERT)
C
                           THTZJ(IR) = THTZJ(IR)
     &                                 - TMAT12*(HZ+HX*CHI_TO(IT,IPERT))
C
                           HZ = HCZ_ZJ_T(IR,K2,K4,IT,IPERT)
                           HX = HCX_ZJ_T(IR,K2,K4,IT,IPERT)
C
                           THTZZ(IR) = THTZZ(IR)
     &                                 - TMAT12*(HZ+HX*CHI_TO(IT,IPERT))
C
                        END DO
C
                     END IF
                  END DO
C
C-----------------------------------------------------------------------
C              include the factor WE/PI for energy integration
C-----------------------------------------------------------------------
C
                  DO IR = 1,IRTOP
                     THTZZ(IR) = WGHTGG*THTZZ(IR)
                     THTZJ(IR) = WGHTGG*THTZJ(IR)
                     THTJZ(IR) = WGHTGG*THTJZ(IR)
                     THTJJ(IR) = WGHTGG*THTJJ(IR)
                  END DO
C
C-----------------------------------------------------------------------
C
                  DO I4 = 1,NCPLWF(LAMB4)
                     IKM4 = IKMCPLWF(I4,LAMB4)
                     L4 = L_IKM(IKM4)
C
                     DO I1 = 1,NCPLWF(LAMA1)
                        IKM1 = IKMCPLWF(I1,LAMA1)
                        L1 = L_IKM(IKM1)
C
C------------------------------------------------------------------------
C    product of radial wave functions and weights ZB(4) * ZA(1) * WAB(r)
C    all terms may contribute to the induced densities
C    NOTE: the angular matrix element to be added is complex
C          for the fully relativistic case
C------------------------------------------------------------------------
C
                        DO IR = 1,IRTOP
                           ZGBZGA = ZGB(IR,I4,LAMB4)*ZGA(IR,I1,LAMA1)
                           ZGBJGA = ZGB(IR,I4,LAMB4)*JGA(IR,I1,LAMA1)
                           JGBZGA = JGB(IR,I4,LAMB4)*ZGA(IR,I1,LAMA1)
                           JGBJGA = JGB(IR,I4,LAMB4)*JGA(IR,I1,LAMA1)
C
                           WGG(IR) = ZGBZGA*THTZZ(IR) + ZGBJGA*THTJZ(IR)
     &                               + JGBZGA*THTZJ(IR)
     &                               + JGBJGA*THTJJ(IR)
C
                           ZFBZFA = ZFB(IR,I4,LAMB4)*ZFA(IR,I1,LAMA1)
                           ZFBJFA = ZFB(IR,I4,LAMB4)*JFA(IR,I1,LAMA1)
                           JFBZFA = JFB(IR,I4,LAMB4)*ZFA(IR,I1,LAMA1)
                           JFBJFA = JFB(IR,I4,LAMB4)*JFA(IR,I1,LAMA1)
C
                           WFF(IR) = ZFBZFA*THTZZ(IR) + ZFBJFA*THTJZ(IR)
     &                               + JFBZFA*THTZJ(IR)
     &                               + JFBJFA*THTJJ(IR)
                        END DO
C
                        DO L3 = ABS(L4-L1),L1 + L4,2
                           LM3 = L3*L3
                           DO M3 = -L3,L3
                              LM3 = LM3 + 1
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
                              DO IOBSE = 1,NOBSE
C
                                 IOPER = IOPER_OBSE(IOBSE)
C
C-----------------------------------------------------------------------
C      check angular matrix element  AME_OP  and
C      KPERT flag for contribution LM to observable   IOBSE
C-----------------------------------------------------------------------
C
                                 AME_G = AMEG_LMOP(IKM4,IKM1,LM3,IOPER)
                                 AME_F = AMEF_LMOP(IKM4,IKM1,LM3,IOPER)
C
                                 IF ( ABS(AME_G).GT.1D-6 .AND. 
     &                                KLMFP(LM3,IT).NE.0 ) THEN
C
C-----------------------------------------------------------------------
C                          induced densities
C-----------------------------------------------------------------------
C
                                    DO IR = 1,IRTOP
C
                                       RHO_GG = DIMAG(AME_G*WGG(IR))
                                       RHO_FF = DIMAG(AME_F*WFF(IR))
C
                                       RHO2NS_GG(IR,LM3,IT,IOBSE,IPERT)
     &                                    = RHO2NS_GG(IR,LM3,IT,IOBSE,
     &                                    IPERT) + RHO_GG + RHO_FF
C
                                    END DO
C
                                 END IF
C
                              END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
                           END DO
                        END DO
C
                     END DO
                  END DO
C-----------------------------------------------------------------------
C
               END IF
            END DO
         END DO
C ======================================================================
C
      END DO
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      END
