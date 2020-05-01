C*==chrdns_rel_fp.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_REL_FP(NDENS,IT,IRTOP,WE,IFILCBWF,TMAT_LOC,
     &                         RHO2NSX,ZGRA,ZFRA,JGRA,JFRA,ZGLA,ZFLA,
     &                         JGLA,JFLA)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <CHRDNS>                                *
C   *                                                                  *
C   *   dealing with FULLY RELATIVISTIC case                           *
C   *                                                                  *
C   *  - if  LHS_SOL_EQ_RHS_SOL = .FALSE.                              *
C   *    ZGLA,ZFLA,JGLA,JFLA correspond to the LHS solution            *
C   *    ZGRA,ZFRA,JGRA,JFRA correspond to the RHS solution            *
C   *                                                                  *
C   *    otherwise LHS = RHS solution with ZG1 = ZG2 etc               *
C   *                                                                  *
C   *  - symmetrization (ZJ+JZ)/2 not needed                           *
C   *                                                                  *
C   *  - if LHS_SOL_EQ_RHS_SOL = .TRUE.                                *
C   *    the matrix ZETA is equal to the UNIT MATRIX                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILCBWF_SPH,IFILCBWF_LHS
      USE MOD_CALCMODE,ONLY:SOLVER_FP,LHS_SOL_EQ_RHS_SOL,BORN_USE_SPH
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX,NCPLWFMAX,KLMFP,IMT,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NKM,NL,NKMMAX,CGCFP,L_IKM,LB_IKM,MUEM05_IKM,
     &    IMKM_IKM,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI,SQRT_2,C0
      USE MOD_RMESH,ONLY:NRMAX,JRNS1
      IMPLICIT NONE
C*--CHRDNS_REL_FP32
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHRDNS_REL_FP')
      LOGICAL SYMMETRIZE_ZJ_JZ
      PARAMETER (SYMMETRIZE_ZJ_JZ=.FALSE.)
      REAL*8 TOLTAU
      PARAMETER (TOLTAU=1D-12)
C
C Dummy arguments
C
      INTEGER IFILCBWF,IRTOP,IT,NDENS
      COMPLEX*16 WE
      COMPLEX*16 JFLA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGLA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGRA(NRMAX,NCPLWFMAX,NKMMAX),TMAT_LOC(NKMMAX,NKMMAX),
     &           ZFLA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZFRA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGLA(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
      REAL*8 RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      COMPLEX*16 CHI(:,:),CHI12,GAMMA(:,:),GAMMA_LHS(:,:),GTIL(:,:),
     &           JF_SPH(:,:,:),JG_SPH(:,:,:),RHOC(:,:),TAU21,WEIGHTF(3),
     &           WEIGHTG(3),WETAU,WMAT(:,:),WZFJF(3),WZFJF_SPH(3),
     &           WZFZF(3),WZFZF_SPH(3),WZGJG(3),WZGJG_SPH(3),WZGZG(3),
     &           WZGZG_SPH(3),ZETA(:,:),ZETA12,ZFJF(:),ZFJF_SPH(:),
     &           ZFZF(:),ZFZF_SPH(:),ZF_SPH(:,:,:),ZGJG(:),ZGJG_SPH(:),
     &           ZGZG(:),ZGZG_SPH(:),ZG_SPH(:,:,:)
      REAL*8 FAC1,FACF(2),FACG(2),FCF,FCG,GAUNTF,GAUNTG,MJ3,MJ4,MS
      REAL*8 GAUNT_CYLM
      INTEGER I,I1,ICWF,ICWF1,ICWF2,IFLAG,IKM1,IKM2,IKM3,IKM4,
     &        IKM_CWFSOL_SPH(:,:),IMKM3,IMKM4,IR,IRNSBOT,IS,J,L0,L3,L4,
     &        LB3,LB4,LM01,LM02,ML0,N,NCWF_SOL_SPH(:)
      LOGICAL IKM1_EQ_IKM2,IKM3_EQ_IKM_CWFSOL_SPH1,
     &        SPHER_SOL_CONTRIBUTES,TAU_NE_0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JF_SPH,JG_SPH,ZF_SPH,ZG_SPH
      ALLOCATABLE CHI,GAMMA,IKM_CWFSOL_SPH,NCWF_SOL_SPH
      ALLOCATABLE ZFJF,ZFZF,ZGJG,ZGZG
      ALLOCATABLE ZFJF_SPH,ZFZF_SPH,ZGJG_SPH,ZGZG_SPH
      ALLOCATABLE RHOC,GTIL,WMAT,GAMMA_LHS,ZETA
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (ZFJF(NRMAX),ZFZF(NRMAX),ZGJG(NRMAX),ZGZG(NRMAX))
      ALLOCATE (RHOC(NRMAX,3))
C
C-----------------------------------------------------------------------
C                            Born solver
C-----------------------------------------------------------------------
C
      IRNSBOT = JRNS1(IMT(IT))
C
      IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
C
         ALLOCATE (GTIL(NKMMAX,NKMMAX),WMAT(NKMMAX,NKMMAX))
         ALLOCATE (GAMMA(NKMMAX,NKMMAX),CHI(NKMMAX,NKMMAX))
         ALLOCATE (GAMMA_LHS(NKMMAX,NKMMAX),ZETA(NKMMAX,NKMMAX))
         ALLOCATE (NCWF_SOL_SPH(NKMMAX),IKM_CWFSOL_SPH(2,NKMMAX))
         ALLOCATE (ZG_SPH(NRMAX,2,NKMMAX),JG_SPH(NRMAX,2,NKMMAX))
         ALLOCATE (ZF_SPH(NRMAX,2,NKMMAX),JF_SPH(NRMAX,2,NKMMAX))
         ALLOCATE (ZGJG_SPH(NRMAX),ZGZG_SPH(NRMAX))
         ALLOCATE (ZFJF_SPH(NRMAX),ZFZF_SPH(NRMAX))
C
         READ (IFILCBWF_SPH,REC=IT) ((CHI(I,J),I=1,NKM),J=1,NKM),
     &                              ((ZETA(I,J),I=1,NKM),J=1,NKM),
     &                              ((GAMMA(I,J),I=1,NKM),J=1,NKM),
     &                              ((GAMMA_LHS(I,J),I=1,NKM),J=1,NKM),
     &                              (NCWF_SOL_SPH(I),I=1,NKM),
     &                              ((IKM_CWFSOL_SPH(I,J),I=1,2),J=1,
     &                              NKM),
     &                              (((ZG_SPH(IR,I,J),ZF_SPH(IR,I,J),
     &                              JG_SPH(IR,I,J),JF_SPH(IR,I,J),IR=1,
     &                              IRTOP),I=1,2),J=1,NKM),IFLAG
C
         IF ( IFLAG.NE.IRTOP )
     &         CALL STOP_MESSAGE(ROUTINE,'IFLAG <> IRTOP')
C
C--------------------------------- GTIL = GAMMA * TMAT_LOC * (GAMMA^x)^T
C
         N = NKM
         GTIL(1:N,1:N) = TRANSPOSE(GAMMA_LHS(1:N,1:N))
         WMAT(1:N,1:N) = MATMUL(TMAT_LOC(1:N,1:N),GTIL(1:N,1:N))
         GTIL(1:N,1:N) = MATMUL(GAMMA(1:N,1:N),WMAT(1:N,1:N))
C
C------------------------------------ combine backscattering term (GTIL)
C                                   with regular single site term  (CHI)
C
         CHI(1:N,1:N) = CHI(1:N,1:N) - GTIL(1:N,1:N)
C
      END IF
C----------------------------------------------------------- Born solver
C
C ------------------------------------ read in wavefunctions for type IT
C
CHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
C
      CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZGLA,ZFLA,JGLA,JFLA,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
C     NOTE: for  LHS_SOL_EQ_RHS_SOL = .TRUE.
C
C     arrays ZG1 and ZG2 are identical in the argument list
C     of the calling routine  <CHRDNS>  -->>  no copying necessary !
C
C         ZGRA(:,:,:) = ZGLA(:,:,:)
C         ZFRA(:,:,:) = ZFLA(:,:,:)
C         JGRA(:,:,:) = ZGLA(:,:,:)
C         JFRA(:,:,:) = ZFLA(:,:,:)
C
      IF ( .NOT.(LHS_SOL_EQ_RHS_SOL) )
     &     CALL WAVFUN_READ_REL(IFILCBWF_LHS,IT,1,ZGRA,ZFRA,JGRA,JFRA,
     &     IRTOP,NCPLWF,IKMCPLWF)
C
CHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
C
C ----------------------------------------------------------------------
C
      DO IKM1 = 1,NKM
         DO IKM2 = 1,NKM
C
            IKM1_EQ_IKM2 = IKM1.EQ.IKM2
C
            TAU_NE_0 = ABS(DREAL(TMAT_LOC(IKM1,IKM2))).GT.TOLTAU .OR. 
     &                 ABS(DIMAG(TMAT_LOC(IKM1,IKM2))).GT.TOLTAU
C
            IF ( IKM1_EQ_IKM2 .OR. TAU_NE_0 ) THEN
C
               TAU21 = TMAT_LOC(IKM2,IKM1)
C
               WETAU = WE*TAU21
C
               DO I = 1,NCPLWF(IKM1)
                  IKM3 = IKMCPLWF(I,IKM1)
                  IMKM3 = IMKM_IKM(IKM3)
                  L3 = L_IKM(IKM3)
                  LB3 = LB_IKM(IKM3)
                  MJ3 = MUEM05_IKM(IKM3) + 0.5D0
C
C----------------------------------------------------------- Born solver
C
                  IKM3_EQ_IKM_CWFSOL_SPH1 = .FALSE.
C
                  IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
C
                     CHI12 = CHI(IKM1,IKM2)
                     ZETA12 = ZETA(IKM1,IKM2)
C
                     DO ICWF = 1,NCWF_SOL_SPH(IKM1)
                        IF ( IKM3.EQ.IKM_CWFSOL_SPH(ICWF,IKM1) ) THEN
                           IKM3_EQ_IKM_CWFSOL_SPH1 = .TRUE.
                           ICWF1 = ICWF
                           EXIT
                        END IF
                     END DO
                  END IF
C----------------------------------------------------------- Born solver
C
                  DO J = 1,NCPLWF(IKM2)
                     IKM4 = IKMCPLWF(J,IKM2)
                     IMKM4 = IMKM_IKM(IKM4)
                     L4 = L_IKM(IKM4)
                     LB4 = LB_IKM(IKM4)
                     MJ4 = MUEM05_IKM(IKM4) + 0.5D0
C
                     DO IS = 1,2
                        FACG(IS) = CGCFP(IKM3,IS)*CGCFP(IKM4,IS)
                        FACF(IS) = CGCFP(IMKM3,IS)*CGCFP(IMKM4,IS)
                     END DO
C
                     IF ( SYMMETRIZE_ZJ_JZ ) THEN
                        DO IR = 1,IRTOP
                           ZGZG(IR) = ZGRA(IR,I,IKM1)*ZGLA(IR,J,IKM2)
                           ZFZF(IR) = ZFRA(IR,I,IKM1)*ZFLA(IR,J,IKM2)
C
                           ZGJG(IR) = (ZGRA(IR,I,IKM1)*JGLA(IR,J,IKM2)+
     &                                JGRA(IR,I,IKM1)*ZGLA(IR,J,IKM2))
     &                                *0.5D0
                           ZFJF(IR) = (ZFRA(IR,I,IKM1)*JFLA(IR,J,IKM2)+
     &                                JFRA(IR,I,IKM1)*ZFLA(IR,J,IKM2))
     &                                *0.5D0
C
                        END DO
                     ELSE
                        DO IR = 1,IRTOP
                           ZGZG(IR) = ZGRA(IR,I,IKM1)*ZGLA(IR,J,IKM2)
                           ZFZF(IR) = ZFRA(IR,I,IKM1)*ZFLA(IR,J,IKM2)
C
                           ZGJG(IR) = ZGRA(IR,I,IKM1)*JGLA(IR,J,IKM2)
                           ZFJF(IR) = ZFRA(IR,I,IKM1)*JFLA(IR,J,IKM2)
                        END DO
                     END IF
C
C----------------------------------------------------------- Born solver
C------ use_spher_sol_A implies: SOLVER_FP.EQ.'BORN'
C                          .AND. BORN_USE_SPH
C
                     SPHER_SOL_CONTRIBUTES = .FALSE.
C
                     IF ( IKM3_EQ_IKM_CWFSOL_SPH1 ) THEN
C
                        DO ICWF = 1,NCWF_SOL_SPH(IKM2)
                           IF ( IKM4.EQ.IKM_CWFSOL_SPH(ICWF,IKM2) ) THEN
                              SPHER_SOL_CONTRIBUTES = .TRUE.
                              ICWF2 = ICWF
                              GOTO 2
                           END IF
                        END DO
                        GOTO 4
C
C
 2                      CONTINUE
                        DO IR = 1,IRNSBOT
                           ZGZG_SPH(IR) = ZG_SPH(IR,ICWF1,IKM1)
     &                        *ZG_SPH(IR,ICWF2,IKM2)
                           ZFZF_SPH(IR) = ZF_SPH(IR,ICWF1,IKM1)
     &                        *ZF_SPH(IR,ICWF2,IKM2)
                           ZGJG_SPH(IR) = ZG_SPH(IR,ICWF1,IKM1)
     &                        *JG_SPH(IR,ICWF2,IKM2)
                           ZFJF_SPH(IR) = ZF_SPH(IR,ICWF1,IKM1)
     &                        *JF_SPH(IR,ICWF2,IKM2)
                        END DO
C
                     END IF
C----------------------------------------------------------- Born solver
C
C
 4                   CONTINUE
                     ML0 = NINT(MJ3-MJ4)
C
                     DO L0 = MIN(ABS(L3-L4),ABS(LB3-LB4)),
     &                  MIN(MAX((L3+L4),(LB3+LB4)),2*NL-2)
C
                        IF ( ABS(ML0).LE.L0 ) THEN
C
C ------------------------------ index used for REAL spherical harmonics
C
                           LM01 = L0**2 + L0 + 1 - ABS(ML0)
                           LM02 = L0**2 + L0 + 1 + ABS(ML0)
C
C-----------------------------------------------------------------------
                           IF ( KLMFP(LM01,IT).NE.0 .OR. KLMFP(LM02,IT)
     &                          .NE.0 ) THEN
C
                              CALL CINIT(3,WEIGHTG)
                              CALL CINIT(3,WEIGHTF)
C
                              DO IS = 1,2
C
                                 MS = DBLE(IS) - 1.5D0
C
                                 GAUNTG = GAUNT_CYLM(L0,ML0,L3,
     &                              NINT(MJ3-MS),L4,NINT(-MJ4+MS))
C
                                 GAUNTF = GAUNT_CYLM(L0,ML0,LB3,
     &                              NINT(MJ3-MS),LB4,NINT(-MJ4+MS))
C
                                 FCG = FACG(IS)*(-1.0D0)**(NINT(MJ4-MS))
     &                                 *GAUNTG
C
                                 FCF = FACF(IS)*(-1.0D0)**(NINT(MJ4-MS))
     &                                 *GAUNTF
C
                                 WEIGHTG(1) = WEIGHTG(1) + FCG*(1.0D0)
                                 WEIGHTG(2) = WEIGHTG(2) + FCG*(2*MS)
                                 WEIGHTG(3) = WEIGHTG(3) + FCG*(MJ3-MS)
C
                                 WEIGHTF(1) = WEIGHTF(1) + FCF*(1.0D0)
                                 WEIGHTF(2) = WEIGHTF(2) + FCF*(2*MS)
                                 WEIGHTF(3) = WEIGHTF(3) + FCF*(MJ3-MS)
C
                              END DO
C
                              FAC1 = -1D0/PI
                              IF ( ML0.NE.0 ) FAC1 = FAC1/SQRT_2
                              IF ( ML0.GT.0 ) FAC1 = FAC1*(-1.0D0)**ML0
C
                              DO I1 = 1,NDENS
                                 WZGZG(I1) = FAC1*WEIGHTG(I1)*WETAU
                                 WZFZF(I1) = FAC1*WEIGHTF(I1)*WETAU
                                 WZGJG(I1) = FAC1*WEIGHTG(I1)*WE
                                 WZFJF(I1) = FAC1*WEIGHTF(I1)*WE
                              END DO
C----------------------------------------------- account for beta-matrix
                              DO I1 = 2,NDENS
                                 WZFZF(I1) = -WZFZF(I1)
                                 WZFJF(I1) = -WZFJF(I1)
                              END DO
C
C=======================================================================
C                      Born solver AND BORN_USE_SPH
C=======================================================================
C
                              IF ( SOLVER_FP.EQ.'BORN' .AND. 
     &                             BORN_USE_SPH ) THEN
C
                                 DO I1 = 1,NDENS
                                    WZGZG_SPH(I1) = WZGJG(I1)*CHI12
                                    WZFZF_SPH(I1) = WZGJG(I1)*CHI12
                                    WZGJG_SPH(I1) = WZGJG(I1)*ZETA12
                                    WZFJF_SPH(I1) = WZGJG(I1)*ZETA12
                                 END DO
C
C -------------------- set up radial functions including  irregular part
C
                                 IF ( IKM1_EQ_IKM2 ) THEN
C
C-------------------- back scattering, regular and irregular single site
C
                                    IF ( SPHER_SOL_CONTRIBUTES ) THEN
                                       DO I1 = 1,NDENS
                                         DO IR = 1,IRNSBOT - 1
                                         RHOC(IR,I1) = -WZGZG_SPH(I1)
     &                                      *ZGZG_SPH(IR)
     &                                      - WZFZF_SPH(I1)*ZFZF_SPH(IR)
     &                                      - WZGJG_SPH(I1)*ZGJG_SPH(IR)
     &                                      - WZFJF_SPH(I1)*ZFJF_SPH(IR)
                                         END DO
                                       END DO
                                    ELSE
                                       RHOC(1:(IRNSBOT-1),1:NDENS) = C0
                                    END IF
C
C------------------------------------------------- standard for r > r_ns
C
                                    DO I1 = 1,NDENS
                                       DO IR = IRNSBOT,IRTOP
C
                                         RHOC(IR,I1) = WZGZG(I1)
     &                                      *ZGZG(IR) + WZFZF(I1)
     &                                      *ZFZF(IR) - WZGJG(I1)
     &                                      *ZGJG(IR) - WZFJF(I1)
     &                                      *ZFJF(IR)
C
                                       END DO
                                    END DO
C
                                 ELSE
C
C -------------------- set up radial functions excluding  irregular part
C
C------------------------------- regular single site and back scattering
C
                                    IF ( SPHER_SOL_CONTRIBUTES ) THEN
                                       DO I1 = 1,NDENS
                                         DO IR = 1,IRNSBOT - 1
                                         RHOC(IR,I1) = -WZGZG_SPH(I1)
     &                                      *ZGZG_SPH(IR)
     &                                      - WZFZF_SPH(I1)*ZFZF_SPH(IR)
     &                                      - WZGJG_SPH(I1)*ZGJG_SPH(IR)
     &                                      - WZFJF_SPH(I1)*ZFJF_SPH(IR)
                                         END DO
                                       END DO
                                    ELSE
                                       RHOC(1:(IRNSBOT-1),1:NDENS) = C0
                                    END IF
C
C------------------------------------------------- standard for r > r_ns
C
                                    DO I1 = 1,NDENS
                                       DO IR = IRNSBOT,IRTOP
C
                                         RHOC(IR,I1) = WZGZG(I1)
     &                                      *ZGZG(IR) + WZFZF(I1)
     &                                      *ZFZF(IR)
C
                                       END DO
                                    END DO
C
                                 END IF
C
C=======================================================================
C                  standard set up of radial functions
C=======================================================================
C
C -------------------- set up radial functions including  irregular part
C
                              ELSE IF ( IKM1_EQ_IKM2 ) THEN
C
C--------------------------------------- single site and back scattering
C
                                 DO I1 = 1,NDENS
                                    DO IR = 1,IRTOP
                                       RHOC(IR,I1) = WZGZG(I1)*ZGZG(IR)
     &                                    + WZFZF(I1)*ZFZF(IR)
     &                                    - WZGJG(I1)*ZGJG(IR)
     &                                    - WZFJF(I1)*ZFJF(IR)
                                    END DO
                                 END DO
C
                              ELSE
C
C -------------------- set up radial functions excluding  irregular part
C
C------------------------------------------------------- back scattering
C
                                 DO I1 = 1,NDENS
                                    DO IR = 1,IRTOP
                                       RHOC(IR,I1) = WZGZG(I1)*ZGZG(IR)
     &                                    + WZFZF(I1)*ZFZF(IR)
                                    END DO
                                 END DO
C
                              END IF
C=============================================== radial functions set up
C=======================================================================
C
C collect contributions to RHO2NS - expanded in REAL spherical harmonics
C
C -------------------------------------------------------------- ML0 < 0
C
                              IF ( ML0.LT.0 ) THEN
C
                                 IF ( KLMFP(LM01,IT).NE.0 ) THEN
                                    DO I1 = 1,NDENS
                                       DO IR = 1,IRTOP
                                         RHO2NSX(IR,LM01,IT,I1)
     &                                      = RHO2NSX(IR,LM01,IT,I1)
     &                                      + DREAL(RHOC(IR,I1))
                                       END DO
                                    END DO
                                 END IF
C
                                 IF ( KLMFP(LM02,IT).NE.0 ) THEN
                                    DO I1 = 1,NDENS
                                       DO IR = 1,IRTOP
                                         RHO2NSX(IR,LM02,IT,I1)
     &                                      = RHO2NSX(IR,LM02,IT,I1)
     &                                      + DIMAG(RHOC(IR,I1))
                                       END DO
                                    END DO
                                 END IF
C
C -------------------------------------------------------------- ML0 > 0
C
                              ELSE IF ( ML0.GT.0 ) THEN
C
                                 IF ( KLMFP(LM01,IT).NE.0 ) THEN
                                    DO IR = 1,IRTOP
                                       RHO2NSX(IR,LM01,IT,1)
     &                                    = RHO2NSX(IR,LM01,IT,1)
     &                                    - DREAL(RHOC(IR,1))
                                       RHO2NSX(IR,LM01,IT,2)
     &                                    = RHO2NSX(IR,LM01,IT,2)
     &                                    - DREAL(RHOC(IR,2))
                                    END DO
                                    IF ( NDENS.EQ.3 ) THEN
                                       DO IR = 1,IRTOP
                                         RHO2NSX(IR,LM01,IT,3)
     &                                      = RHO2NSX(IR,LM01,IT,3)
     &                                      - DREAL(RHOC(IR,3))
                                       END DO
                                    END IF
                                 END IF
C
                                 IF ( KLMFP(LM02,IT).NE.0 ) THEN
                                    DO IR = 1,IRTOP
                                       RHO2NSX(IR,LM02,IT,1)
     &                                    = RHO2NSX(IR,LM02,IT,1)
     &                                    + DIMAG(RHOC(IR,1))
                                       RHO2NSX(IR,LM02,IT,2)
     &                                    = RHO2NSX(IR,LM02,IT,2)
     &                                    + DIMAG(RHOC(IR,2))
                                    END DO
                                    IF ( NDENS.EQ.3 ) THEN
                                       DO IR = 1,IRTOP
                                         RHO2NSX(IR,LM02,IT,3)
     &                                      = RHO2NSX(IR,LM02,IT,3)
     &                                      + DIMAG(RHOC(IR,3))
                                       END DO
                                    END IF
                                 END IF
C
C -------------------------------------------------------------- ML0 = 0
C
                              ELSE
                                 DO I1 = 1,NDENS
                                    DO IR = 1,IRTOP
                                       RHO2NSX(IR,LM01,IT,I1)
     &                                    = RHO2NSX(IR,LM01,IT,I1)
     &                                    + DIMAG(RHOC(IR,I1))
                                    END DO
                                 END DO
                              END IF
C
                           END IF
C-----------------------------------------------------------------------
                        END IF
                     END DO
                  END DO
               END DO
C
            END IF
         END DO
      END DO
C ======================================================================
C
      END
C*==chrdns_sra_fp.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHRDNS_SRA_FP(IT,IRTOP,WE,IFILCBWF,TMAT_LOC,RHO2NSX,
     &                         ZGRA,JGRA)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for   <CHRDNS>                                *
C   *                                                                  *
C   *   dealing with SCALAR RELATIVISTIC case                          *
C   *                                                                  *
C   *  NOTE: for IREL <= 2 (non- and scalar relativistic) the          *
C   *        minor component has no meaning and is not written/read    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL,SOLVER_FP,BORN_USE_SPH
      USE MOD_FILES,ONLY:IFILCBWF_SPH
      USE MOD_RMESH,ONLY:NRMAX,JRNS1
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX,NCPLWFMAX,KLMFP,NLT,IMT,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NLM,L_LM,M_LM,NSPIN,NKMMAX,NL,NCPLWF
      USE MOD_CONSTANTS,ONLY:PI,SQRT_4PI
      IMPLICIT NONE
C*--CHRDNS_SRA_FP575
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHRDNS_SRA_FP')
      REAL*8 TOLTAU
      PARAMETER (TOLTAU=1D-12)
C
C Dummy arguments
C
      INTEGER IFILCBWF,IRTOP,IT
      COMPLEX*16 WE
      COMPLEX*16 JGRA(NRMAX,NCPLWFMAX,NKMMAX),TMAT_LOC(NKMMAX,NKMMAX),
     &           ZGRA(NRMAX,NCPLWFMAX,NKMMAX)
      REAL*8 RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3)
C
C Local variables
C
      COMPLEX*16 CHIS(:,:,:),GAMMAS(:,:,:),GTIL(:,:),H0LS(:,:,:),
     &           R0LS(:,:,:),R1R2,TAU,WGHT,WGHTZJ,WGHTZZ,WMAT(:,:)
      REAL*8 GAUNT_RYLM
      REAL*8 GNT,RHORH,RHORR,RHOZJ,RHOZZ,SPNWGT,WRH(:),WRR(:),WZJ(:),
     &       WZZ(:)
      INTEGER I,I1,I2,IL,IL1,IL2,ILM1,ILM2,IR,IRNSBOT,IS,J,JLM1,JLMS1,
     &        JLMS2,L1,L2,L3,LM3,LMSOFF,M1,M2,M3,N
      LOGICAL JLMS1_EQ_JLMS2,TAU_NE_0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE R0LS,H0LS,CHIS,GAMMAS,WRR,WRH,WZJ,WZZ,GTIL,WMAT
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (WZZ(NRMAX),WZJ(NRMAX))
C
      WGHTZJ = -WE/PI
C
C-----------------------------------------------------------------------
C                            Born solver
C-----------------------------------------------------------------------
C
      IRNSBOT = JRNS1(IMT(IT))
C
      IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
C
         ALLOCATE (WRR(NRMAX),WRH(NRMAX))
         ALLOCATE (GTIL(NLM,NLM),WMAT(NLM,NLM))
         ALLOCATE (R0LS(NRMAX,NL,NSPIN),GAMMAS(NLM,NLM,NSPIN))
         ALLOCATE (H0LS(NRMAX,NL,NSPIN),CHIS(NLM,NLM,NSPIN))
C
         READ (IFILCBWF_SPH,REC=IT) (((CHIS(I,J,IS),I=1,NLM),J=1,NLM),((
     &                              GAMMAS(I,J,IS),I=1,NLM),J=1,NLM),
     &                              ((R0LS(IR,IL,IS),H0LS(IR,IL,IS),
     &                              IR=1,IRTOP),IL=1,NL),IS=1,NSPIN)
C
C----------------- spherical part of single site contribution  r <= r_ns
C
         DO IS = 1,NSPIN
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
C
            WRH(:) = 0D0
C
            DO IL = 1,NLT(IT)
C
               WGHT = -(WE/PI)*(2*(IL-1)+1)/SQRT_4PI
C
               DO IR = 1,IRNSBOT - 1
                  WRH(IR) = WRH(IR)
     &                      + DIMAG(WGHT*R0LS(IR,IL,IS)*H0LS(IR,IL,IS))
               END DO
C
            END DO
C
            DO IR = 1,IRNSBOT - 1
               RHORH = WRH(IR)
C
               RHO2NSX(IR,1,IT,1) = RHO2NSX(IR,1,IT,1) - RHORH
               RHO2NSX(IR,1,IT,2) = RHO2NSX(IR,1,IT,2) - SPNWGT*RHORH
            END DO
C
         END DO
C
      END IF
C----------------------------------------------------------- Born solver
C
C ----------------------------------------- read in wavefunctions for IT
C
      CALL WAVFUN_READ_SRA(IFILCBWF,IT,1,ZGRA,JGRA,IRTOP,NCPLWF,
     &                     IKMCPLWF)
C
C-----------------------------------------------------------------------
C
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
C----------------------------------------------------------- Born solver
C                                      GTIL = GAMMA * TMAT_LOC * GAMMA^T
         IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
            I1 = LMSOFF + 1
            I2 = LMSOFF + NLM
            N = NLM
            GTIL(1:N,1:N) = TRANSPOSE(GAMMAS(1:N,1:N,IS))
            WMAT(1:N,1:N) = MATMUL(TMAT_LOC(I1:I2,I1:I2),GTIL(1:N,1:N))
            GTIL(1:N,1:N) = MATMUL(GAMMAS(1:N,1:N,IS),WMAT(1:N,1:N))
C
C------------------------------------ combine backscattering term (GTIL)
C                                   with regular single site term  (CHI)
C
CcccccccccccCHIS(1:N,1:N,IS) = CHIS(1:N,1:N,IS) - GTIL(1:N,1:N)
C
         END IF
C-----------------------------------------------------------------------
C=======================================================================
C
         DO JLMS2 = LMSOFF + 1,LMSOFF + NLM
            DO JLMS1 = LMSOFF + 1,LMSOFF + NLM
               JLM1 = JLMS1 - LMSOFF
C
               JLMS1_EQ_JLMS2 = JLMS1.EQ.JLMS2
C
               TAU_NE_0 = ABS(DREAL(TMAT_LOC(JLMS1,JLMS2)))
     &                    .GT.TOLTAU .OR. 
     &                    ABS(DIMAG(TMAT_LOC(JLMS1,JLMS2))).GT.TOLTAU
C
               IF ( JLMS1_EQ_JLMS2 .OR. TAU_NE_0 ) THEN
C
                  TAU = TMAT_LOC(JLMS1,JLMS2)
C
                  WGHTZZ = -WE*TAU/PI
C
C-----------------------------------------------------------------------
C
                  DO I2 = 1,NCPLWF(JLMS2)
                     ILM2 = IKMCPLWF(I2,JLMS2)
                     L2 = L_LM(ILM2)
                     M2 = M_LM(ILM2)
                     IL2 = L2 + 1
C
                     DO I1 = 1,NCPLWF(JLMS1)
                        ILM1 = IKMCPLWF(I1,JLMS1)
                        L1 = L_LM(ILM1)
                        M1 = M_LM(ILM1)
                        IL1 = L1 + 1
C
                        DO IR = 1,IRTOP
                           WZZ(IR) = DIMAG(ZGRA(IR,I1,JLMS1)*ZGRA(IR,I2,
     &                               JLMS2)*WGHTZZ)
                           WZJ(IR) = DIMAG(ZGRA(IR,I1,JLMS1)*JGRA(IR,I2,
     &                               JLMS2)*WGHTZJ)
                        END DO
C
C----------------------------------------------------------- Born solver
                        IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH )
     &                       THEN
                           IF ( ILM1.EQ.JLM1 ) THEN
C
                              WGHT = -(WE/PI)*CHIS(ILM1,ILM2,IS)
C
                              DO IR = 1,IRNSBOT - 1
                                 R1R2 = R0LS(IR,IL1,IS)*R0LS(IR,IL2,IS)
                                 WRR(IR) = DIMAG(WGHT*R1R2)
                              END DO
C
                           END IF
                        END IF
C----------------------------------------------------------- Born solver
C
                        DO L3 = ABS(L2-L1),L1 + L2,2
                           LM3 = L3*L3
                           DO M3 = -L3,L3
                              LM3 = LM3 + 1
                              IF ( KLMFP(LM3,IT).NE.0 ) THEN
C
                                 GNT = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C
                                 IF ( ABS(GNT).GT.1D-6 ) THEN
C
C ======================================================================
C                radial functions excluding  irregular part
C ======================================================================
C
                                    IF ( .NOT.JLMS1_EQ_JLMS2 ) THEN
C
                                       DO IR = 1,IRTOP
                                         RHOZZ = GNT*WZZ(IR)
                                         RHO2NSX(IR,LM3,IT,1)
     &                                      = RHO2NSX(IR,LM3,IT,1)
     &                                      + RHOZZ
                                         RHO2NSX(IR,LM3,IT,2)
     &                                      = RHO2NSX(IR,LM3,IT,2)
     &                                      + SPNWGT*RHOZZ
                                       END DO
C
C
C ======================================================================
C                radial functions including  irregular part
C ======================================================================
C
C-----------------------------------------------------------------------
C                            Born solver
C-----------------------------------------------------------------------
                                    ELSE IF ( SOLVER_FP.EQ.'BORN' .AND. 
     &                                 BORN_USE_SPH ) THEN
C
C---------------- backscattering and single site contribution  r >= r_ns
C
                                       DO IR = IRNSBOT,IRTOP
                                         RHOZZ = GNT*WZZ(IR)
                                         RHOZJ = GNT*WZJ(IR)
                                         RHO2NSX(IR,LM3,IT,1)
     &                                      = RHO2NSX(IR,LM3,IT,1)
     &                                      + RHOZZ - RHOZJ
                                         RHO2NSX(IR,LM3,IT,2)
     &                                      = RHO2NSX(IR,LM3,IT,2)
     &                                      + SPNWGT*(RHOZZ-RHOZJ)
                                       END DO
C
C-------------------------------- backscattering contribution  r <= r_ns
C
                                       DO IR = 1,IRNSBOT - 1
                                         RHOZZ = GNT*WZZ(IR)
                                         RHO2NSX(IR,LM3,IT,1)
     &                                      = RHO2NSX(IR,LM3,IT,1)
     &                                      + RHOZZ
                                         RHO2NSX(IR,LM3,IT,2)
     &                                      = RHO2NSX(IR,LM3,IT,2)
     &                                      + SPNWGT*RHOZZ
                                       END DO
C
C------------- NON-spherical part of single site contribution  r <= r_ns
C
                                       IF ( ILM1.EQ.JLM1 ) THEN
C
                                         DO IR = 1,IRNSBOT - 1
                                         RHORR = GNT*WRR(IR)
C
                                         RHO2NSX(IR,LM3,IT,1)
     &                                      = RHO2NSX(IR,LM3,IT,1)
     &                                      - RHORR
                                         RHO2NSX(IR,LM3,IT,2)
     &                                      = RHO2NSX(IR,LM3,IT,2)
     &                                      - SPNWGT*RHORR
                                         END DO
C
                                       END IF
C
                                    ELSE
C
C-----------------------------------------------------------------------
C                         standard charge construction
C-----------------------------------------------------------------------
C
                                       DO IR = 1,IRTOP
                                         RHOZZ = GNT*WZZ(IR)
                                         RHOZJ = GNT*WZJ(IR)
                                         RHO2NSX(IR,LM3,IT,1)
     &                                      = RHO2NSX(IR,LM3,IT,1)
     &                                      + RHOZZ - RHOZJ
                                         RHO2NSX(IR,LM3,IT,2)
     &                                      = RHO2NSX(IR,LM3,IT,2)
     &                                      + SPNWGT*(RHOZZ-RHOZJ)
                                       END DO
C
                                    END IF
C
                                 END IF
C
                              END IF
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
      END
