C*==linresp_init.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_INIT(LINRESP_TASK)
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine for phonon spectra                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,C0
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,NSF,LMISF,NPAN,JRCUT,R2DRDI,
     &    FLMSF,NLMSFMAX,NLSF
      USE MOD_FILES,ONLY:WRTAU,RDTAU,IFILCBWF,LRECREAL8,WR_TRACK_INFO
      USE MOD_ANGMOM,ONLY:NL,NLM,NKM,NKMMAX,L_LM,WKM1,WKM2,IMKM_IKM,
     &    AME_G,AME_F
      USE MOD_TYPES,ONLY:NTMAX,NLFP,NLMFP,NLMFPMAX,NT,IMT,NLMFPT,
     &    NCPLWFMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_LINRESP,ONLY:NPERT,NOPER,IOPER_PERT,NOBSE,IOPER_OBSE,
     &    SHAZ_ZZ_T,SHBZ_JZ_T,SHCZ_ZJ_T,SHDZ_JJ_T,HAZ_ZZ_T,HBZ_JZ_T,
     &    HCZ_ZJ_T,HDZ_JJ_T,HAX_ZZ_T,HBX_JZ_T,HCX_ZJ_T,HDX_JJ_T,
     &    HZ_PERT_LMTP,HX_PERT_LMTP,AMEG_LMOP,AMEF_LMOP,K_PERT_LMTP,
     &    RHO2NS_GG,MZBZA_O,MIRR2_OP,MIRR3_OP,MIRR4_OP,D0Z,D1Z,DZ,DIJZ,
     &    D0Z1,D1Z1,DZ1,DIJZ1,D0X,DIJX,T0Z,T1Z,TZ,TIJZ,T0X,TIJX,D1X,T1X,
     &    DX,TX,IDOS,ISPN,IORB,DOBS_TO,CHI_TO,GNTTAB,TZ_STD,TZ_DK,
     &    T1Z_DK,T0Z_DK
      IMPLICIT NONE
C*--LINRESP_INIT27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_INIT')
C
C Dummy arguments
C
      CHARACTER*10 LINRESP_TASK
C
C Local variables
C
      REAL*8 G123,HX_LM(:,:),HZ_LM(:,:)
      REAL*8 GAUNT_RYLM
      INTEGER IA_ERR,IKM1,IKM2,IM,IMKM1,IMKM2,IOBSE,IOPER,IPERT,
     &        IPERT_SIG_POL,IPOL,IR,IRCRIT,IRMTIN,IRSF,ISF,IT,L1,L2,L3,
     &        LDV,LH,LM1,LM2,LM3,LMDV,LMFP,LMH,LMSF,LSF,M1,M2,M3,MH,
     &        RECLNGWF
      LOGICAL INITIALIZE
      COMPLEX*16 MAUX(:,:)
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE HZ_LM,HX_LM,MAUX
C
      CALL TRACK_INFO(ROUTINE)
C
      WR_TRACK_INFO = .TRUE.
C
      ALLOCATE (HZ_LM(NRMAX,NLMFPMAX),HX_LM(NRMAX,NLMFPMAX))
      ALLOCATE (MAUX(NKMMAX,NKMMAX))
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SUSC
      IF ( LINRESP_TASK.EQ.'SUSC_STAT ' ) THEN
C
         NOPER = 3
         NOBSE = 3
         NPERT = 3
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM MAGNET
      ELSE IF ( LINRESP_TASK.EQ.'MAGNET    ' ) THEN
C
C-----------------------------------------------------------------------
C        IPERT
C
C        1        DOS
C        2 -  4   SPN
C        5 -  7   ORB
C        8 - 16   ( (Y_LM1 x sigma_IPOL), LM1=2,4), IPOL=1,3)
C       17 - 19             (sigma_IPOL),           IPOL=1,3)
C       20 - 22             (nabla_IPOL),           IPOL=1,3)
C
C                  LM1 =  2,  3,  4
C                   m  = -1,  0, +1
C                   c  =  y,  z,  x
C-----------------------------------------------------------------------
C
         NOPER = 19
         NOBSE = 1
         NPERT = 22
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP PHONON
      ELSE IF ( LINRESP_TASK.EQ.'PHONONS   ' ) THEN
C
         NOPER = 3
         NOBSE = 1
         NPERT = 1
C
      ELSE
C
         WRITE (6,*) '     LINRESP_TASK = ',LINRESP_TASK
         CALL STOP_MESSAGE(ROUTINE,'LINRESP_TASK unknown')
C
      END IF
C
C=======================================================================
C                           INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C-----------------------------------------------------------------------
C             Gaunt coefficients for REAL spherical harmonics
C-----------------------------------------------------------------------
C
         ALLOCATE (GNTTAB(NLMFPMAX,NLMFPMAX,NLMSFMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate GNTTAB')
         GNTTAB(:,:,:) = 0D0
C
         LM1 = 0
         DO L1 = 0,(NLFP-1)
            DO M1 = -L1, + L1
               LM1 = LM1 + 1
C
               LM2 = 0
               DO L2 = 0,(NLFP-1)
                  DO M2 = -L2, + L2
                     LM2 = LM2 + 1
C
                     LM3 = 0
                     DO L3 = 0,4*(NL-1)
                        DO M3 = -L3, + L3
                           LM3 = LM3 + 1
C
                           GNTTAB(LM1,LM2,LM3)
     &                        = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C
                        END DO
                     END DO
C
                  END DO
               END DO
C
            END DO
         END DO
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
C=======================================================================
C             angular matrix element for operator  IOPER
C             AMEG_LMOP(NKMMAX,NKMMAX,NLMFPMAX,NOPER)
C ======================================================================
C
      ALLOCATE (AMEG_LMOP(NKMMAX,NKMMAX,NLMFPMAX,NOPER))
C
      IOPER = 1
C
      AMEG_LMOP(:,:,:,:) = C0
C
      IF ( IREL.LE.2 ) THEN
C
         DO LMFP = 1,NLMFP
            DO LM1 = 1,NLM
               DO LM2 = 1,NLM
                  AMEG_LMOP(LM1,LM2,LMFP,IOPER)
     &               = DCMPLX(GNTTAB(LM1,LM2,LMFP),0D0)
               END DO
            END DO
         END DO
C
      ELSE
C
         ALLOCATE (AMEF_LMOP(NKMMAX,NKMMAX,NLMFPMAX,NOPER))
         AMEF_LMOP(:,:,:,:) = C0
C
         DO LMFP = 1,NLMFP
C
            WKM1(:,:) = C0
            DO LM1 = 1,NLM
               DO LM2 = 1,NLM
                  WKM1(LM1,LM2) = DCMPLX(GNTTAB(LM1,LM2,LMFP),0D0)
                  WKM1(NLM+LM1,NLM+LM2) = WKM1(LM1,LM2)
               END DO
            END DO
C
            CALL CHANGEREP(NKM,NKMMAX,WKM1,'RLM>REL',WKM2)
C
            AMEG_LMOP(1:NKM,1:NKM,LMFP,IOPER) = WKM2(1:NKM,1:NKM)
C
C-----------------------------------------------------------------------
C     set up the angular matrix element for the small component
C-----------------------------------------------------------------------
C
            DO IKM1 = 1,NKM
               IMKM1 = IMKM_IKM(IKM1)
               DO IKM2 = 1,NKM
                  IMKM2 = IMKM_IKM(IKM2)
C
                  IF ( IMKM1.LE.NKM .AND. IMKM2.LE.NKM )
     &                 AMEF_LMOP(IKM1,IKM2,LMFP,IOPER)
     &                 = AMEG_LMOP(IMKM1,IMKM2,LMFP,IOPER)
C
               END DO
            END DO
         END DO
C-----------------------------------------------------------------------
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SUSC
         IF ( LINRESP_TASK.EQ.'SUSC_STAT ' ) THEN
C
            DO LMFP = 1,NLMFP
C
               AMEG_LMOP(1:NKM,1:NKM,LMFP,ISPN)
     &            = MATMUL(AMEG_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &            AME_G(1:NKM,1:NKM,2,ISPN))
C
               AMEF_LMOP(1:NKM,1:NKM,LMFP,ISPN)
     &            = MATMUL(AMEF_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &            AME_F(1:NKM,1:NKM,2,ISPN))
C
               AMEG_LMOP(1:NKM,1:NKM,LMFP,IORB)
     &            = MATMUL(AMEG_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &            AME_G(1:NKM,1:NKM,2,IORB))
C
               AMEF_LMOP(1:NKM,1:NKM,LMFP,IORB)
     &            = MATMUL(AMEF_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &            AME_F(1:NKM,1:NKM,2,IORB))
C
            END DO
C
            LMFP = 1
            CALL CMATSTRUCT('AMEG (DOS)',AMEG_LMOP(1,1,LMFP,IDOS),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (DOS)',AMEF_LMOP(1,1,LMFP,IDOS),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
            CALL CMATSTRUCT('AMEG (SPN)',AMEG_LMOP(1,1,LMFP,ISPN),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (SPN)',AMEF_LMOP(1,1,LMFP,ISPN),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
            CALL CMATSTRUCT('AMEG (ORB)',AMEG_LMOP(1,1,LMFP,IORB),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (ORB)',AMEF_LMOP(1,1,LMFP,IORB),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
         END IF
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM MAGNET
         IF ( LINRESP_TASK.EQ.'MAGNET    ' ) THEN
C
C-----------------------------------------------------------------------
C        IPERT
C
C        1        DOS
C        2 -  4   SPN
C        5 -  7   ORB
C        8 - 16   ( (Y_LM1 x sigma_IPOL), LM1=2,4), IPOL=1,3)
C       17 - 19             (sigma_IPOL),           IPOL=1,3)
C       20 - 22             (nabla_IPOL),           IPOL=1,3)
C
C                  LM1 =  2,  3,  4
C                   m  = -1,  0, +1
C                   c  =  y,  z,  x
C-----------------------------------------------------------------------
C
            DO LMFP = 1,NLMFP
C
               IPERT = 1
C
C------------------------------------------------------------------- SPN
C
C                                        sigma_lamba  lambda = -1, 0, +1
               DO IPOL = 1,3
C
                  IPERT = IPERT + 1
C
                  AMEG_LMOP(1:NKM,1:NKM,LMFP,IPERT)
     &               = MATMUL(AMEG_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &               AME_G(1:NKM,1:NKM,IPOL,ISPN))
C
                  AMEF_LMOP(1:NKM,1:NKM,LMFP,IPERT)
     &               = MATMUL(AMEF_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &               AME_F(1:NKM,1:NKM,IPOL,ISPN))
C
               END DO
C
C------------------------------------------------------------------- ORB
C
C                                            l_lamba  lambda = -1, 0, +1
               DO IPOL = 1,3
C
                  IPERT = IPERT + 1
C
                  AMEG_LMOP(1:NKM,1:NKM,LMFP,IPERT)
     &               = MATMUL(AMEG_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &               AME_G(1:NKM,1:NKM,IPOL,IORB))
C
                  AMEF_LMOP(1:NKM,1:NKM,LMFP,IPERT)
     &               = MATMUL(AMEF_LMOP(1:NKM,1:NKM,LMFP,IDOS),
     &               AME_F(1:NKM,1:NKM,IPOL,IORB))
               END DO
C
C--------------------------------------------------------------- r x ALF
C
C                                        sigma_lamba  lambda = -1, 0, +1
               DO IPOL = 1,3
C
                  IPERT_SIG_POL = 1 + IPOL
C
C                             real spherical harmonic Y_1m m = -1, 0, +1
                  DO LM1 = 2,4
C
C    <LAM| Y_1m sig_lambda |LAM'> =
C                    SUM[LAM'']  <LAM|Y_1m|LAM''><LAM''|sig_lambda|LAM'>
C
                     MAUX(1:NKM,1:NKM)
     &                  = MATMUL(AMEG_LMOP(1:NKM,1:NKM,LM1,IDOS),
     &                  AMEG_LMOP(1:NKM,1:NKM,LMFP,IPERT_SIG_POL))
C
C     AMEG_{LAM,LAM'} =  < LAM | Y_1m sig_lambda |-LAM'>
C     AMEF_{LAM,LAM'} =  <-LAM | Y_1m sig_lambda | LAM'>
C
                     IPERT = IPERT + 1
C
                     DO IKM1 = 1,NKM
C
                        IMKM1 = IMKM_IKM(IKM1)
C
                        DO IKM2 = 1,NKM
C
                           IMKM2 = IMKM_IKM(IKM2)
C
                           IF ( IMKM2.LE.NKM )
     &                          AMEG_LMOP(IKM1,IKM2,LMFP,IPERT)
     &                          = MAUX(IKM1,IMKM2)
C
                           IF ( IMKM1.LE.NKM )
     &                          AMEF_LMOP(IKM1,IKM2,LMFP,IPERT)
     &                          = MAUX(IMKM1,IKM2)
C
                        END DO
                     END DO
C
                  END DO
C
               END DO
C
C------------------------------------------------------------------- ALF
C
C                                        sigma_lamba  lambda = -1, 0, +1
               DO IPOL = 1,3
C
                  IPERT_SIG_POL = 1 + IPOL
C
                  IPERT = IPERT + 1
C
                  DO IKM1 = 1,NKM
C
                     IMKM1 = IMKM_IKM(IKM1)
C
                     DO IKM2 = 1,NKM
C
                        IMKM2 = IMKM_IKM(IKM2)
C
                        IF ( IMKM2.LE.NKM )
     &                       AMEG_LMOP(IKM1,IKM2,LMFP,IPERT)
     &                       = AMEG_LMOP(IKM1,IMKM2,LMFP,IPERT_SIG_POL)
C
                        IF ( IMKM1.LE.NKM )
     &                       AMEF_LMOP(IKM1,IKM2,LMFP,IPERT)
     &                       = AMEG_LMOP(IMKM1,IKM2,LMFP,IPERT_SIG_POL)
C
                     END DO
                  END DO
C
               END DO
C
            END DO
C
C---------------------------------------------------- print z-components
            LMFP = 1
            CALL CMATSTRUCT('AMEG (DOS)',AMEG_LMOP(1,1,LMFP,IDOS),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (DOS)',AMEF_LMOP(1,1,LMFP,IDOS),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
            CALL CMATSTRUCT('AMEG (SPN)',AMEG_LMOP(1,1,LMFP,3),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (SPN)',AMEF_LMOP(1,1,LMFP,3),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
            CALL CMATSTRUCT('AMEG (ORB)',AMEG_LMOP(1,1,LMFP,6),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
            CALL CMATSTRUCT('AMEF (ORB)',AMEF_LMOP(1,1,LMFP,6),NKM,
     &                      NKMMAX,3,3,1,1D-8,6)
C
         END IF
C
C-----------------------------------------------------------------------
C              regular matrix elements for orbial current
C-----------------------------------------------------------------------
C
CCCCCCCCCCCCCCCCCCCCC         CALL AMENAB
         CALL AME_NAB
         CALL FPAMENAB
C
         WRITE (6,*) '# regular matrix elements for orbial current DONE'
C
      END IF
C=======================================================================
C
C
C=======================================================================
C                       perturbation terms
C=======================================================================
C
      ALLOCATE (HZ_PERT_LMTP(NRMAX,NLMFPMAX,NTMAX,NPERT))
      ALLOCATE (HX_PERT_LMTP(NRMAX,NLMFPMAX,NTMAX,NPERT))
      ALLOCATE (K_PERT_LMTP(NLMFPMAX,NTMAX,NPERT))
C
      IF ( LINRESP_TASK.EQ.'MAGNET    ' ) THEN
         ALLOCATE (SHAZ_ZZ_T(NKMMAX,NKMMAX,NTMAX,NPERT))
         ALLOCATE (SHBZ_JZ_T(NKMMAX,NKMMAX,NTMAX,NPERT))
         ALLOCATE (SHCZ_ZJ_T(NKMMAX,NKMMAX,NTMAX,NPERT))
         ALLOCATE (SHDZ_JJ_T(NKMMAX,NKMMAX,NTMAX,NPERT))
      END IF
C
      ALLOCATE (HAZ_ZZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HBZ_JZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HCZ_ZJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HDZ_JJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),STAT=IA_ERR)
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate HDZ_JJ')
C
      ALLOCATE (HAX_ZZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HBX_JZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HCX_ZJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (HDX_JJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),STAT=IA_ERR)
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate HDX_JJ')
C
      HZ_PERT_LMTP(:,1:1,1:NT,1:NPERT) = SQRT_4PI
      HZ_PERT_LMTP(:,2:NLMFPMAX,1:NT,1:NPERT) = 0D0
C
      HX_PERT_LMTP(:,1:1,1:NT,1:NPERT) = SQRT_4PI
      HX_PERT_LMTP(:,2:NLMFPMAX,1:NT,1:NPERT) = 0D0
C
C
      K_PERT_LMTP(:,:,:) = .FALSE.
      K_PERT_LMTP(1,1:NT,1:NPERT) = .TRUE.
C
      DO IT = 1,NT
C
         IM = IMT(IT)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCUT(NPAN(IM),IM)
C
         DO IPERT = 1,NPERT
C
            HZ_LM(:,:) = 0D0
            HX_LM(:,:) = 0D0
C
C---------------------------------- convolute HPERT with shape functions
C
C------------------------------------- perturbation HLM in the MT regime
C
            DO LMDV = 1,NLMFPT(IT)
               IF ( K_PERT_LMTP(LMDV,IT,IPERT) ) THEN
                  DO IR = 1,IRMTIN
                     HZ_LM(IR,LMDV) = HZ_PERT_LMTP(IR,LMDV,IT,IPERT)
     &                                *R2DRDI(IR,IM)
                     HX_LM(IR,LMDV) = HX_PERT_LMTP(IR,LMDV,IT,IPERT)
     &                                *R2DRDI(IR,IM)
                  END DO
               END IF
            END DO
C
C------------------------------------- perturbation HLM in the IS regime
C
            DO LMDV = 1,NLMFPT(IT)
               LDV = L_LM(LMDV)
               IF ( K_PERT_LMTP(LMDV,IT,IPERT) ) THEN
C
                  DO ISF = 1,NSF(IM)
                     LMSF = LMISF(ISF,IM)
                     LSF = L_LM(LMSF)
                     IF ( LSF.LE.(NLSF-1) ) THEN
C
                        DO LH = ABS(LSF-LDV),MIN((LSF+LDV),(NLFP-1)),2
                           LMH = LH*LH
                           DO MH = -LH,LH
                              LMH = LMH + 1
C
                              IF ( LMH.NE.(LH*LH+LH+MH+1) )
     &                             CALL STOP_MESSAGE(ROUTINE,'# LMH #')
C
                              G123 = GNTTAB(LMDV,LMH,LMSF)
C
                              IF ( ABS(G123).GT.1D-10 ) THEN
C
                                 K_PERT_LMTP(LMH,IT,IPERT) = .TRUE.
C
                                 DO IR = IRMTIN + 1,IRCRIT
                                    IRSF = IR - IRMTIN
C
                                    HZ_LM(IR,LMH) = HZ_LM(IR,LMH)
     &                                 + HZ_PERT_LMTP(IR,LMDV,IT,IPERT)
     &                                 *R2DRDI(IR,IM)*FLMSF(IRSF,ISF,IM)
     &                                 *G123
C
                                    HX_LM(IR,LMH) = HX_LM(IR,LMH)
     &                                 + HX_PERT_LMTP(IR,LMDV,IT,IPERT)
     &                                 *R2DRDI(IR,IM)*FLMSF(IRSF,ISF,IM)
     &                                 *G123
C
                                 END DO
C
                              END IF
C
                           END DO
                        END DO
C
                     END IF
                  END DO
C
               END IF
C
            END DO
C
            HZ_PERT_LMTP(1:IRCRIT,1:NLMFPT(IT),IT,IPERT)
     &         = HZ_LM(1:IRCRIT,1:NLMFPT(IT))
            HX_PERT_LMTP(1:IRCRIT,1:NLMFPT(IT),IT,IPERT)
     &         = HX_LM(1:IRCRIT,1:NLMFPT(IT))
C
         END DO
      END DO
C
C ======================================================================
C                           observables
C ======================================================================
C
      ALLOCATE (RHO2NS_GG(NRMAX,NLMFPMAX,NTMAX,NOBSE,NPERT))
C
      ALLOCATE (MIRR2_OP(NKMMAX,NKMMAX,NOBSE,NPERT))
      ALLOCATE (MIRR3_OP(NKMMAX,NKMMAX,NOBSE,NPERT))
      ALLOCATE (MIRR4_OP(NKMMAX,NKMMAX,NOBSE,NPERT))
      ALLOCATE (MZBZA_O(NKMMAX,NKMMAX,NOBSE))
C
C ======================================================================
C     pointer associating   operator  with   perturbation term  IPERT
C ======================================================================
      ALLOCATE (IOPER_PERT(NPERT))
C
      DO IPERT = 1,NPERT
         IOPER_PERT(IPERT) = IPERT
      END DO
C
C ======================================================================
C     pointer associating   operator  with   observable  IOBSE
C ======================================================================
C
      ALLOCATE (IOPER_OBSE(NOBSE))
C
      DO IOBSE = 1,NOBSE
         IOPER_OBSE(IOBSE) = IOBSE
      END DO
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate THXT_DIA')
C
C ======================================================================
C     partial response terms used for integration
C ======================================================================
C      ALLOCATE (D0ZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (D1ZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (DZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (DIJZL(NLMAX,NTMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (D0ZL1(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (D1ZL1(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (DZL1(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (DIJZL1(NLMAX,NTMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (D0XL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (DIJXL(NLMAX,NTMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (T0ZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (T1ZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (TZL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (TIJZL(NLMAX,NTMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (T0XL(NLMAX,NTMAX,NOBSE,NPERT))
C      ALLOCATE (TIJXL(NLMAX,NTMAX,NTMAX,NOBSE,NPERT))
C
      ALLOCATE (D0Z(NTMAX,NOBSE,NPERT))
      ALLOCATE (D1Z(NTMAX,NOBSE,NPERT))
      ALLOCATE (DZ(NTMAX,NOBSE,NPERT))
      ALLOCATE (DIJZ(NTMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (D0Z1(NTMAX,NOBSE,NPERT))
      ALLOCATE (D1Z1(NTMAX,NOBSE,NPERT))
      ALLOCATE (DZ1(NTMAX,NOBSE,NPERT))
      ALLOCATE (DIJZ1(NTMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (D0X(NTMAX,NOBSE,NPERT))
      ALLOCATE (D1X(NTMAX,NOBSE,NPERT))
      ALLOCATE (DX(NTMAX,NOBSE,NPERT))
      ALLOCATE (DIJX(NTMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (T0Z(NTMAX,NOBSE,NPERT))
      ALLOCATE (T1Z(NTMAX,NOBSE,NPERT))
      ALLOCATE (TZ(NTMAX,NOBSE,NPERT))
      ALLOCATE (TIJZ(NTMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (T0X(NTMAX,NOBSE,NPERT))
      ALLOCATE (T1X(NTMAX,NOBSE,NPERT))
      ALLOCATE (TX(NTMAX,NOBSE,NPERT))
      ALLOCATE (TIJX(NTMAX,NTMAX,NOBSE,NPERT))
      ALLOCATE (TZ_STD(NTMAX,NOBSE,NPERT))
      ALLOCATE (TZ_DK(NTMAX,NOBSE,NPERT,3))
      ALLOCATE (T0Z_DK(NTMAX,NOBSE,NPERT,3))
      ALLOCATE (T1Z_DK(NTMAX,NOBSE,NPERT,3))
C
      ALLOCATE (DOBS_TO(NTMAX,NOBSE))
      ALLOCATE (CHI_TO(NTMAX,NOBSE))
C ======================================================================
C
      WRTAU = .FALSE.
      RDTAU = .FALSE.
C
      CLOSE (IFILCBWF)
C
      IF ( FULLPOT ) THEN
         RECLNGWF = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         RECLNGWF = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
C-----------------------------------------------------------------------
C
      END
