C*==comptonme0.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COMPTONME0(C,ZG,ZF,JG,JF,I_PW_ZR,I_ZL_PW,I_JL_PW,
     &                      F_PW_ZR,F_ZL_PW,F_JL_PW,ZGR,ZFR,JGR,JFR,ZGL,
     &                      ZFL,JGL,JFL,NCPLWF_T,IKMCPLWF_T,TASK,NP,
     &                      PFEMAX,ZGPOS,RME_PW_ZR_P,RME_ZL_PW_P,
     &                      RME_IRR_P,ZECH,ALPH,BET,GAM,DELT,NPMAX,
     &                      PPOINTS)
C     ******************************************************************
C     *                                                                *
C     *                         PART 0                                 *
C     *                                                                *
C     *   speed up the calculation of the matrix elements by           *
C     *   interpolating the p-dependent part                           *
C     *                                                                *
C     *----------------------------------------------------------------*
C     *                                                                *
C     *   calculate the RADIAL matrix elements                         *
C     *   for the momentum representation as a function of p           *
C     *                                                                *
C     *            RME_PW_ZR_P = < p m_s | Z_Lam^R >                   *
C     *                                                                *
C     *            RME_ZL_PW_P = < Z_Lam^L | p m_s >                   *
C     *                                                                *
C     ******************************************************************
C
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      USE MOD_TYPES,ONLY:NT,NTMAX,NLT,IMT,NCPLWFMAX
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R,JRCRI,FULLPOT
      USE MOD_ANGMOM,ONLY:NL,NKMMAX,NKM,L_IKM,LB_IKM,KAPPA_IKM
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--COMPTONME034
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTONME0')
C
C Dummy arguments
C
      REAL*8 C,PFEMAX
      INTEGER NP,NPMAX,PPOINTS
      CHARACTER*10 TASK
      REAL*8 ALPH(PPOINTS,0:NL),BET(PPOINTS,0:NL),DELT(PPOINTS,0:NL),
     &       GAM(PPOINTS,0:NL),ZECH(PPOINTS+1)
      COMPLEX*16 F_JL_PW(NRMAX),F_PW_ZR(NRMAX),F_ZL_PW(NRMAX),
     &           I_JL_PW(NRMAX),I_PW_ZR(NRMAX),I_ZL_PW(NRMAX),
     &           JF(NRMAX,NCPLWFMAX,NKM),JFL(NRMAX,NCPLWFMAX,NKM),
     &           JFR(NRMAX,NCPLWFMAX,NKM),JG(NRMAX,NCPLWFMAX,NKM),
     &           JGL(NRMAX,NCPLWFMAX,NKM),JGR(NRMAX,NCPLWFMAX,NKM),
     &           RME_IRR_P(NPMAX,NKMMAX,NKMMAX,NTMAX),
     &           RME_PW_ZR_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX),
     &           RME_ZL_PW_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX),
     &           ZF(NRMAX,NCPLWFMAX,NKM),ZFL(NRMAX,NCPLWFMAX,NKM),
     &           ZFR(NRMAX,NCPLWFMAX,NKM),ZG(NRMAX,NCPLWFMAX,NKM),
     &           ZGL(NRMAX,NCPLWFMAX,NKM),ZGPOS(NRMAX,NTMAX),
     &           ZGR(NRMAX,NCPLWFMAX,NKM)
      INTEGER IKMCPLWF_T(NCPLWFMAX,NKMMAX,NTMAX),NCPLWF_T(NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 AG(:,:),AG0(:,:),AJ_LF,AJ_LG,F1,F2,FIRR(:),IIRR
      REAL*8 CSQR,EPFE,PFE,WBET,WDELT,WGAM,WKP0,WKP1,WKP2,XIN,YY
      INTEGER IFIL_LHS,IFIL_RHS,IKM1,IKM2,IM,IP,IR,IRTOP,IT,J,JMESH,K1,
     &        K2,L,LAM,LF,LG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE  FIRR,AG,AG0
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (FIRR(NRMAX))
      ALLOCATE (AG(NRMAX,0:NL),AG0(NRMAX,0:NL))
C
      CSQR = C**2
C
      RME_IRR_P(:,:,:,:) = C0
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
C!!!!!!!!!!!!!!!!!!!!! TO BE CHECKED
C
Cc      CALL SET_IFIL_LHS(IFIL_RHS,IFIL_LHS)
C
      IFIL_RHS = IFILCBWF
      IFIL_LHS = IFILCBWF
C
C=======================================================================
C
C  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IM = IMT(IT)
C
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
         CALL WAVFUN_READ_REL(IFIL_RHS,IT,1,ZG,ZF,JG,JF,IRTOP,
     &                        NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
         CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZG,ZF,JG,JF,ZGR,ZFR,JGR,JFR,
     &                               NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
         IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
            ZGL(:,:,:) = ZGR(:,:,:)
            ZFL(:,:,:) = ZFR(:,:,:)
            JGL(:,:,:) = JGR(:,:,:)
            JFL(:,:,:) = JFR(:,:,:)
C
         ELSE
C
            CALL WAVFUN_READ_REL(IFIL_LHS,IT,1,ZG,ZF,JG,JF,IRTOP,
     &                           NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
            CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZG,ZF,JG,JF,ZGL,ZFL,JGL,JFL,
     &                                  NCPLWF_T(1,IT),
     &                                  IKMCPLWF_T(1,1,IT))
C
         END IF
C
C  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
         DO IP = 1,NP
C
            PFE = PFEMAX*DBLE(IP-1)/DBLE(NP-1)
C
            EPFE = (CSQR/2D0)*(SQRT(1D0+4D0*PFE**2/CSQR)-1D0)
C
            WKP0 = C*PFE/(EPFE+CSQR)
C
C=======================================================================
C    get AG0 = j_l(pr) by spline interpolation
C=======================================================================
C
            DO IR = 1,IRTOP
               IF ( PFE*R(IRTOP,IM)-ZECH(PPOINTS+1).GT.1D-14 )
     &              CALL STOP_MESSAGE(ROUTINE,'Error in meshes')
               XIN = PFE*R(IR,IM)
               JMESH = 1
               DO J = 1,PPOINTS
                  IF ( XIN.GE.ZECH(J+1) ) JMESH = JMESH + 1
               END DO
               IF ( XIN.GE.ZECH(PPOINTS+1) ) JMESH = PPOINTS
               IF ( XIN.LT.ZECH(1) ) JMESH = 1
C
C---------- define the weights for calculating the intrepolated function
C
               WBET = XIN - ZECH(JMESH)
               WGAM = WBET*WBET
               WDELT = WGAM*WBET
C
               DO L = 0,NL
                  YY = ALPH(JMESH,L) + BET(JMESH,L)*WBET + GAM(JMESH,L)
     &                 *WGAM + DELT(JMESH,L)*WDELT
C
                  AG0(IR,L) = YY
                  AG(IR,L) = AG0(IR,L)
               END DO
            END DO
C
C=======================================================================
C          include positronic wave function for  TASK = 'POSANI'
C=======================================================================
            IF ( TASK(1:6).EQ.'POSANI' ) THEN
               DO L = 0,(NLT(IT)-1)
                  DO IR = 1,IRTOP
                     AG(IR,L) = AG0(IR,L)*ZGPOS(IR,IT)
                  END DO
               END DO
            END IF
C=======================================================================
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO LAM = 1,NKM
C
C***********************************************************************
C           REGULAR PART   < PW | Z > and  < Z | PW >
C***********************************************************************
C
C K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2
               DO K2 = 1,NCPLWF_T(LAM,IT)
C
                  IKM2 = IKMCPLWF_T(K2,LAM,IT)
C
                  LG = L_IKM(IKM2)
                  LF = LB_IKM(IKM2)
                  WKP2 = WKP0*SIGN(1,KAPPA_IKM(IKM2))
C
                  DO IR = 1,IRTOP
C
                     AJ_LG = AG(IR,LG)
                     AJ_LF = AG(IR,LF)*WKP2
C
                     F_PW_ZR(IR) = AJ_LG*ZGR(IR,K2,LAM)
     &                             + AJ_LF*ZFR(IR,K2,LAM)
                     F_ZL_PW(IR) = ZGL(IR,K2,LAM)*AJ_LG + ZFL(IR,K2,LAM)
     &                             *AJ_LF
                     F_JL_PW(IR) = JGL(IR,K2,LAM)*AJ_LG + JFL(IR,K2,LAM)
     &                             *AJ_LF
                  END DO
C
                  CALL CRADINT_R(IM,F_PW_ZR,I_PW_ZR)
                  CALL CRADINT_R(IM,F_ZL_PW,I_ZL_PW)
C
                  RME_PW_ZR_P(IP,K2,LAM,IT) = I_PW_ZR(IRTOP)
                  RME_ZL_PW_P(IP,K2,LAM,IT) = I_ZL_PW(IRTOP)
C
C***********************************************************************
C           IRRREGULAR PART   < PW | ZJ + JZ | PW >
C***********************************************************************
C           NOTE: a summation is done w.r.t to the index LAM
C
C K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1
                  IF ( .NOT.SPLITSS ) THEN
C
                     CALL CRADINT_INW_R(IM,F_JL_PW,I_JL_PW)
C
                     DO K1 = 1,NCPLWF_T(LAM,IT)
C
                        IKM1 = IKMCPLWF_T(K1,LAM,IT)
C
                        LG = L_IKM(IKM1)
                        LF = LB_IKM(IKM1)
                        WKP1 = WKP0*SIGN(1,KAPPA_IKM(IKM1))
C
                        DO IR = 1,IRTOP
C
                           AJ_LG = AG(IR,LG)
                           AJ_LF = AG(IR,LF)*WKP1
C
                           F1 = (AJ_LG*ZGR(IR,K1,LAM)
     &                          +AJ_LF*ZFR(IR,K1,LAM))*I_JL_PW(IR)
C
                           F2 = (AJ_LG*JGR(IR,K1,LAM)
     &                          +AJ_LF*JFR(IR,K1,LAM))*I_ZL_PW(IR)
C
                           FIRR(IR) = F1 + F2
C
                        END DO
C
                        CALL CRADINT(IM,FIRR,IIRR)
C
                        RME_IRR_P(IP,IKM1,IKM2,IT)
     &                     = RME_IRR_P(IP,IKM1,IKM2,IT) + IIRR
C
                     END DO
                  END IF
C K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1K1
C
               END DO
C K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2
C
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         END DO
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
