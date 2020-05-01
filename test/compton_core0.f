C*==compton_core0.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COMPTON_CORE0(C,NCPLWF_CORT,IKMCPLWF_CORT,TASK,NP,
     &                         PFEMAX,ZGPOS,RME_COR_P,ZECH,ALPH,BET,GAM,
     &                         DELT,NPMAX,PPOINTS)
C     ******************************************************************
C     *                                                                *
C     *                    CORE CONTRIBUTION                           *
C     *                                                                *
C     *     to Compton or positron or positron annihilation spectrum   *
C     *                                                                *
C     *----------------------------------------------------------------*
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
C     *            RME_COR_P = < p m_s | Z_Lam >                       *
C     *                                                                *
C     ******************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRWS,R,JRCRI,FLMSF,JRMT,R2DRDI
      USE MOD_ANGMOM,ONLY:NL,LB_IKM,L_IKM,KAPPA_IKM
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,Z,IMT,ITBOT,ITTOP,NLT,NCORMAX,
     &    LCXRAY,NCXRAY,NCORT,TXT_T
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:Y00,RY_EV
      IMPLICIT NONE
C*--COMPTON_CORE034
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTON_CORE0')
      INTEGER NLSHELLMAX,NCSTMAX
      PARAMETER (NLSHELLMAX=15,NCSTMAX=14)
C
C Dummy arguments
C
      REAL*8 C,PFEMAX
      INTEGER NP,NPMAX,PPOINTS
      CHARACTER*10 TASK
      REAL*8 ALPH(PPOINTS,0:NL),BET(PPOINTS,0:NL),DELT(PPOINTS,0:NL),
     &       GAM(PPOINTS,0:NL),ZECH(PPOINTS+1)
      INTEGER IKMCPLWF_CORT(NCPLWFMAX,NCORMAX,NTMAX),
     &        NCPLWF_CORT(NCORMAX,NTMAX)
      COMPLEX*16 RME_COR_P(NPMAX,NCPLWFMAX,NCORMAX,NTMAX),
     &           ZGPOS(NRMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 AG(:,:),AG0(:,:),AJ_LF,AJ_LG,CAUX,F_AUX(:)
      REAL*8 BCOR(:),BCORS(:),CSQR,ECOR(:),EPFE,FCOR(:,:,:),FLM_A_RGNT,
     &       GCOR(:,:,:),PFE,RINT(:),SZCOR(:),WBET,WDELT,WGAM,WKP0,WKP1,
     &       XIN,XNORM(2),YY
      INTEGER IA_ERR,ICOR,ICOR0,ICST,IKM,IKMCOR(:,:),ILSHELL,IM,IP,IR,
     &        IRCRIT,IRSF,IRTOP,IRTOP_SPHERE,IT,IT_CORE,IZERO(:),J,
     &        JMESH,KAPCOR(:),KK,L,LF,LG,LQNTAB(NLSHELLMAX),MM05COR(:),
     &        NCST,NKPCOR(:),NLSHELL,NQNTAB(NLSHELLMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AG,AG0,F_AUX,RINT
      ALLOCATABLE BCOR,BCORS,NKPCOR
      ALLOCATABLE FCOR,GCOR,ECOR,SZCOR,IZERO,KAPCOR,MM05COR,IKMCOR
C
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (BCOR(NTMAX),BCORS(NTMAX))
      ALLOCATE (IZERO(NCSTMAX),SZCOR(NCSTMAX),ECOR(NCSTMAX))
      ALLOCATE (KAPCOR(NCSTMAX),MM05COR(NCSTMAX))
      ALLOCATE (IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX))
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GCOR')
      ALLOCATE (AG(NRMAX,0:NL),AG0(NRMAX,0:NL))
      ALLOCATE (F_AUX(NRMAX),RINT(NRMAX))
C
      CSQR = C**2
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         WRITE (6,99001) IT,TXT_T(IT)
C
         IM = IMT(IT)
C
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
            IRTOP_SPHERE = JRMT(IM)
            IRCRIT = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
            IRTOP_SPHERE = JRWS(IM)
            IRCRIT = JRWS(IM)
         END IF
C
         NLSHELL = 0
         IF ( Z(IT).GT.2 ) NLSHELL = 1
         IF ( Z(IT).GT.10 ) NLSHELL = 3
         IF ( Z(IT).GT.18 ) NLSHELL = 5
         IF ( Z(IT).GT.30 ) NLSHELL = 6
         IF ( Z(IT).GT.36 ) NLSHELL = 8
         IF ( Z(IT).GT.48 ) NLSHELL = 9
         IF ( Z(IT).GT.54 ) NLSHELL = 11
         IF ( Z(IT).GT.70 ) NLSHELL = 12
         IF ( Z(IT).GT.80 ) NLSHELL = 13
         IF ( Z(IT).GT.86 ) NLSHELL = 15
C
         IT_CORE = IT
C
         ICOR0 = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         LOOP_ILSHELL:DO ILSHELL = 1,NLSHELL
C
            NCXRAY(IT) = NQNTAB(ILSHELL)
            LCXRAY(IT) = LQNTAB(ILSHELL)
C
            NCST = 4*LCXRAY(IT) + 2
C
            IF ( NCST.GT.NCSTMAX )
     &            CALL STOP_MESSAGE(ROUTINE,' NCST > NCSTMAX')
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,IT_CORE,BCOR,BCORS,NCSTMAX)
C
            ICOR = ICOR0
C
            DO ICST = 1,NCST
C
               ICOR = ICOR + 1
C
               NCPLWF_CORT(ICOR,IT) = NKPCOR(ICST)
C
               DO KK = 1,NKPCOR(ICST)
C
                  DO IR = 1,IRTOP
                     RINT(IR) = R2DRDI(IR,IM)
     &                          *(GCOR(IR,KK,ICST)**2+FCOR(IR,KK,ICST)
     &                          **2)
                  END DO
C
                  CALL RRADINT(IM,RINT,XNORM(KK))
C
                  IKMCPLWF_CORT(KK,ICOR,IT) = IKMCOR(ICST,KK)
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
C                          FLM_A_RGNT = F00SF*A_RGNT = 1
C
                  DO IR = 1,IRTOP_SPHERE
                     FLM_A_RGNT = R2DRDI(IR,IM)
                     GCOR(IR,KK,ICST) = GCOR(IR,KK,ICST)*FLM_A_RGNT
                     FCOR(IR,KK,ICST) = FCOR(IR,KK,ICST)*FLM_A_RGNT
                  END DO
C
C--------------------------------------------------- interstitial regime
C
                  DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                     IRSF = IR - IRTOP_SPHERE
                     FLM_A_RGNT = FLMSF(IRSF,1,IM)*Y00*R2DRDI(IR,IM)
C
                     GCOR(IR,KK,ICST) = GCOR(IR,KK,ICST)*FLM_A_RGNT
                     FCOR(IR,KK,ICST) = FCOR(IR,KK,ICST)*FLM_A_RGNT
C
                  END DO
C
               END DO
C
               WRITE (6,99002) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                         (2*MM05COR(ICST)+1),IKMCOR(ICST,1),
     &                         XNORM(1),ECOR(ICST),ECOR(ICST)*RY_EV,
     &                         SZCOR(ICST),IZERO(ICST)
               IF ( NKPCOR(ICST).EQ.2 ) WRITE (6,99003) IKMCOR(ICST,2),
     &              XNORM(2)
C
            END DO
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
     &                 CALL STOP_MESSAGE(ROUTINE,'Error in meshes')
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
                     YY = ALPH(JMESH,L) + BET(JMESH,L)
     &                    *WBET + GAM(JMESH,L)*WGAM + DELT(JMESH,L)
     &                    *WDELT
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
               ICOR = ICOR0
C
               DO ICST = 1,NCST
C
                  ICOR = ICOR + 1
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
                  DO KK = 1,NKPCOR(ICST)
C
                     IKM = IKMCOR(ICST,KK)
C
                     LG = L_IKM(IKM)
                     LF = LB_IKM(IKM)
                     WKP1 = WKP0*SIGN(1,KAPPA_IKM(IKM))
C
                     DO IR = 1,IRTOP
C
                        AJ_LG = AG(IR,LG)
                        AJ_LF = AG(IR,LF)*WKP1
C
                        F_AUX(IR) = AJ_LG*GCOR(IR,KK,ICST)
     &                              + AJ_LF*FCOR(IR,KK,ICST)
                     END DO
C
                     CALL CRADINT(IM,F_AUX,CAUX)
C
                     RME_COR_P(IP,KK,ICOR,IT) = CAUX
C
                  END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
               END DO
C
            END DO
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
            ICOR0 = ICOR0 + NCST
C
         END DO LOOP_ILSHELL
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         IF ( ICOR0.NE.NCORT(IT) )
     &         CALL STOP_MESSAGE(ROUTINE,'ICOR0 .NE. NCORT(IT)')
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      WRITE (6,'(//)')
C
99001 FORMAT (//,10X,'core states for IT =',I4,4X,A,//,
     &        ' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99002 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99003 FORMAT (22X,I4,F12.6)
      END
