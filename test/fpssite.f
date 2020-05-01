C*==fpssite.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE FPSSITE(IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,ERYD,P,
     &                   IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C   ********************************************************************
C   *                                                                  *
C   * ASSIGN QUANTUM NUMBERS AND CALL ROUTINE TO SOLVE                 *
C   * COUPLED PARTIAL DIFFERENTIAL RADIAL DIRAC EQUATIONS. THE         *
C   * RESULTING WAVEFUNCTIONS ARE USED TO CALCULATE T-MATRICES IN      *
C   * THE KAPPA-MU REPRESENTATION                                      *
C   *                                                                  *
C   * + CALCULATION OF THE RADIAL INTEGRALS                            *
C   *   [ G1*G2 + F1*F2 ] R**2 DR                                      *
C   *                                                                  *
C   * FOR IHYPER <> 0 :                                                *
C   * CALCULATION OF THE HYPERFINE MATRIX ELEMENTS                     *
C   *                                                                  *
C   * 08/07/10 HE write sequence for wave functions changed            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FILNAM,LFILNAM,IOTMP,IOTMP1,IFILCBWF_LHS,
     &    IFILCBWF_SPH,IFILGFWF,IFILMEZZL
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NLMAX,NCPLWF,L_IKM,KAPPA_IKM,
     &    NKM,WKM1,WKM2,IPIVKM,NLM
      USE MOD_RMESH,ONLY:JRNSMIN,NPANMAX,NRMAX,NPAN,JRNS1,JRCUT,JRCRI,
     &    DRDIOVR,DRDI,R
      USE MOD_TYPES,ONLY:CTL,NTMAX,NCPLWFMAX,NLMFPMAX,VAMEF,VAMEG,
     &    IKMSOLBLK,LMIFP,NBLK,NFPT,NSOLBLK,ISOLIKM,IKMCPLWF,NLMFPT,
     &    KLMFP,LOPT,BNST,VNST,Z,BT,VT,IMT,NT,ITBOT,ITTOP,NLT,NKM_T
      USE MOD_CALCMODE,ONLY:SOLVER_FP,GF_CONV_RH,BREAKPOINT,
     &    LHS_SOL_EQ_RHS_SOL,PROGNAME,IHYPER,L_PROJECTED_ME,KMROT
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,SEVT,SEBT,SEBNST,SEVNST
      USE MOD_CONSTANTS,ONLY:CI,C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPSSITE')
      LOGICAL WRONSKI
      PARAMETER (WRONSKI=.FALSE.)
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      LOGICAL GETIRRSOL
      INTEGER IFILSS,IPRINT,IWRIRRWF,IWRREGWF
      CHARACTER*10 ORBPOLSS
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ABLK(:,:),ARG,BBLK(:,:),BBLKINV(:,:),CBLK(:,:),
     &           CBNST(:,:),CBT(:),CHI(:,:),CHL,CHLM(:),CHLP(:),CJL,
     &           CJL_J,CMAT(:,:),CMAT_LHS(:,:),CNL,CREF,CRSQ,CSUM,
     &           CVNST(:,:),CVT(:),DBLK(:,:),DMAT(:,:),DMAT_LHS(:,:),
     &           DP(:,:,:),DPM(:,:,:),DQ(:,:,:),DQM(:,:,:),D_WRON(:,:,:)
     &           ,GAMMA(:,:),GAMMA_LHS(:,:),GBLK(:,:),HBLK(:,:),
     &           JF(:,:,:),JFL(:,:,:),JF_SPH(:,:,:),JG(:,:,:),JGL(:,:,:)
     &           ,JG_SPH(:,:,:),MEZJL(:,:,:,:),MEZZL(:,:,:,:),PFAC,
     &           PIM(:,:,:),PR(:,:,:),PRM(:,:,:),QIM(:,:,:),QR(:,:,:),
     &           QRM(:,:,:),RTOP,SBLK(:,:),SBLK0(:,:),SSST0(:,:),
     &           TBLK(:,:),TBLK0(:,:),TSST0(:,:),TSST_BORN(:,:),
     &           TSST_LHS(:,:),W(:,:),WFAC,ZETA(:,:),ZF(:,:,:),
     &           ZFL(:,:,:),ZF_SPH(:,:,:),ZG(:,:,:),ZGL(:,:,:),
     &           ZG_SPH(:,:,:)
      REAL*8 C,CSQR,SK1
      LOGICAL CALC_GFWF,DUMP_SSITE,LDAU,MODIFYV
      COMPLEX*16 CJLZ,CNLZ
      INTEGER I,IA_ERR,IBLK,IBLK1,IBLK2,ICPL,IFILSS_XHS,IFMT,IKM,IKM1,
     &        IKMBOT,IKMTOP,IKM_CWFSOL_SPH(:,:),IM,INFO,IPIV(:),IPOT,IR,
     &        IRTOP,ISK1,ISOL,IT,I_HAND_SIDE,J,J1,J2,JKM,JKM1,JKM2,K,
     &        KAP1,L,L1,LB1,LJ,LM,NCWF_SOL_SPH(:),NSOL,N_HAND_SIDE
C
C*** End of declarations rewritten by SPAG
C
      DATA DUMP_SSITE/.FALSE./
C
      ALLOCATABLE W,DP,CBT,ABLK,BBLK,CBLK,DBLK,HBLK,GBLK
      ALLOCATABLE CVT,BBLKINV,IPIV,CBNST,CVNST
      ALLOCATABLE PIM,PRM,QIM,QRM,DPM,DQM,CHLM,CHLP
      ALLOCATABLE DQ,QR,PR,TBLK,SBLK,TBLK0,SBLK0
      ALLOCATABLE ZF,ZG,JF,JG,MEZJL,MEZZL
      ALLOCATABLE ZGL,ZFL,JGL,JFL
      ALLOCATABLE TSST_LHS,CMAT,CMAT_LHS,DMAT,DMAT_LHS,GAMMA_LHS
      ALLOCATABLE TSST0,SSST0,JF_SPH,JG_SPH,ZF_SPH,ZG_SPH
      ALLOCATABLE IKM_CWFSOL_SPH,NCWF_SOL_SPH
      ALLOCATABLE TSST_BORN,CHI,GAMMA,ZETA,D_WRON
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (W(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (DP(2,2,NRMAX),CBT(NRMAX),ABLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (BBLK(NCPLWFMAX,NCPLWFMAX),CVT(NRMAX))
      ALLOCATE (CBLK(NCPLWFMAX,NCPLWFMAX),HBLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (DBLK(NCPLWFMAX,NCPLWFMAX),GBLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (BBLKINV(NCPLWFMAX,NCPLWFMAX),IPIV(NCPLWFMAX))
      ALLOCATE (CBNST(JRNSMIN:NRMAX,NLMFPMAX))
      ALLOCATE (CVNST(JRNSMIN:NRMAX,NLMFPMAX))
      ALLOCATE (CHLM(0:NLMAX),CHLP(0:NLMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (DQ(2,2,NRMAX),QR(2,2,NRMAX),PR(2,2,NRMAX))
      ALLOCATE (TBLK(NCPLWFMAX,NCPLWFMAX),TBLK0(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (SBLK(NCPLWFMAX,NCPLWFMAX),SBLK0(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX))
C
      ALLOCATE (PIM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      ALLOCATE (PRM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      ALLOCATE (QIM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      ALLOCATE (QRM(NCPLWFMAX,NCPLWFMAX,NRMAX),STAT=IA_ERR)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) THEN
         IPRINT = 5
         DUMP_SSITE = .TRUE.
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( IPRINT.GT.3 ) CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP1,
     &     'zzzzzz_check-matrices')
C=======================================================================
C
C-----------------------------------------------------------------------
C              set the calculation mode and file chanels
C-----------------------------------------------------------------------
      IF ( IFILSS.EQ.IFILGFWF ) THEN
C
         CALC_GFWF = .TRUE.
C
      ELSE
C
         CALC_GFWF = .FALSE.
C
         CALL SET_IFIL_LHS(IFILSS,IFILCBWF_LHS)
         CALL SET_IFIL_SPH(IFILSS,IFILCBWF_SPH)
C
      END IF
C-----------------------------------------------------------------------
C
CHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      ALLOCATE (JFL(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGL(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFL(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGL(NRMAX,NCPLWFMAX,NKMMAX))
C
      IF ( LHS_SOL_EQ_RHS_SOL .OR. CALC_GFWF ) THEN
         N_HAND_SIDE = 1
      ELSE
         N_HAND_SIDE = 2
         ALLOCATE (TSST_LHS(NKMMAX,NKMMAX))
      END IF
      ALLOCATE (DMAT_LHS(NKMMAX,NKMMAX))
      ALLOCATE (GAMMA_LHS(NKMMAX,NKMMAX))
      ALLOCATE (CMAT_LHS(NKMMAX,NKMMAX))
      IF ( WRONSKI ) THEN
         IF ( N_HAND_SIDE.EQ.1 ) THEN
            ALLOCATE (JFL(NRMAX,NCPLWFMAX,NKMMAX))
            ALLOCATE (JGL(NRMAX,NCPLWFMAX,NKMMAX))
            ALLOCATE (ZFL(NRMAX,NCPLWFMAX,NKMMAX))
            ALLOCATE (ZGL(NRMAX,NCPLWFMAX,NKMMAX))
            ALLOCATE (TSST_LHS(NKMMAX,NKMMAX))
         END IF
         ALLOCATE (D_WRON(NKMMAX,NKMMAX,NTMAX))
      END IF
CHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      IF ( SOLVER_FP.NE.'BS' .AND. SOLVER_FP.NE.'BORN' )
     &     CALL STOP_MESSAGE(ROUTINE,'solver  '//SOLVER_FP//
     &                       '  not found')
C-----------------------------------------------------------------------
C
      IF ( SOLVER_FP.EQ.'BS' ) THEN
         ALLOCATE (DPM(NCPLWFMAX,NCPLWFMAX,NRMAX))
         ALLOCATE (DQM(NCPLWFMAX,NCPLWFMAX,NRMAX),STAT=IA_ERR)
         ALLOCATE (NCWF_SOL_SPH(1),IKM_CWFSOL_SPH(1,1))
         ALLOCATE (ZG_SPH(1,1,1),JG_SPH(1,1,1))
         ALLOCATE (ZF_SPH(1,1,1),JF_SPH(1,1,1))
         ALLOCATE (CHI(NKMMAX,NKMMAX))
      ELSE
         ALLOCATE (DMAT(NKMMAX,NKMMAX),TSST_BORN(NKMMAX,NKMMAX))
         ALLOCATE (CMAT(NKMMAX,NKMMAX),ZETA(NKMMAX,NKMMAX))
         ALLOCATE (GAMMA(NKMMAX,NKMMAX),CHI(NKMMAX,NKMMAX))
         ALLOCATE (TSST0(NKMMAX,NKMMAX),SSST0(NKMMAX,NKMMAX))
         ALLOCATE (ZG_SPH(NRMAX,2,NKMMAX),JG_SPH(NRMAX,2,NKMMAX))
         ALLOCATE (ZF_SPH(NRMAX,2,NKMMAX),JF_SPH(NRMAX,2,NKMMAX))
         ALLOCATE (NCWF_SOL_SPH(NKMMAX))
         ALLOCATE (IKM_CWFSOL_SPH(2,NKMMAX),STAT=IA_ERR)
      END IF
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc: PIM')
      NCWF_SOL_SPH(:) = 999999
      IKM_CWFSOL_SPH(:,:) = 999999
      ZG_SPH(1,1,1) = 999999D0
      ZF_SPH(1,1,1) = 999999D0
      JG_SPH(1,1,1) = 999999D0
      JF_SPH(1,1,1) = 999999D0
C
C=======================================================================
C           prepare calculation of l-projected matrix elements
C=======================================================================
C
      L_PROJECTED_ME = .TRUE.
C
      IF ( PROGNAME(4:6).EQ.'SCF' ) THEN
         IF ( IHYPER.NE.0 ) THEN
            L_PROJECTED_ME = .TRUE.
         ELSE
            L_PROJECTED_ME = .FALSE.
         END IF
      END IF
C
      IF ( L_PROJECTED_ME ) THEN
         ALLOCATE (MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
         ALLOCATE (MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
      ELSE
         ALLOCATE (MEZJL(1,1,1,1),MEZZL(1,1,1,1))
      END IF
C
      IF ( L_PROJECTED_ME .AND. .NOT.CALC_GFWF ) REWIND IFILMEZZL
C
C=======================================================================
C
      C = CTL(1,1)
      CSQR = C*C
C
C     calculate relativistic momentum
C
      CALL GET_MOMENTUM(3,C,ERYD,P)
C
      PFAC = P/(1.0D0+ERYD/CSQR)
      WFAC = C*(1.0D0+ERYD/CSQR)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRTOP = JRCRI(IM)
C
         DO I = 1,IRTOP
            DRDIOVR(I,IM) = DRDI(I,IM)/R(I,IM)
         END DO
C
C HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
         DO I_HAND_SIDE = 1,N_HAND_SIDE
C
            MEZZ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
            MEZJ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
C
            ZG(:,:,:) = C0
            ZF(:,:,:) = C0
            JG(:,:,:) = C0
            JF(:,:,:) = C0
C
            TSST(:,:,IT) = C0
            MSST(:,:,IT) = C0
            SSST(:,:,IT) = C0
            DO J = 1,NKMMAX
C----------- fill diagonal with 1 to allow matrix inversion for NLT < NL
C                 NOTE: this dummy setting will be overwritten up to NLT
               TSST(J,J,IT) = 1D0
               MSST(J,J,IT) = 1D0
               SSST(J,J,IT) = 1D0
            END DO
C
C------------------------------------------------------------------ BORN
            IF ( SOLVER_FP.EQ.'BORN' ) THEN
               TSST_BORN(:,:) = C0
               CHI(:,:) = C0
               GAMMA(:,:) = C0
               CMAT(:,:) = C0
               DMAT(:,:) = C0
            END IF
C------------------------------------------------------------------ BORN
C
            RTOP = R(IRTOP,IM)
            ARG = P*RTOP
C
            DO L1 = 0,NLT(IT)
               CJL = CJLZ(L1,ARG)
               CNL = CNLZ(L1,ARG)
               CHLP(L1) = CJL + CI*CNL
               CHLM(L1) = CJL - CI*CNL
            END DO
C
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            DO IBLK = 1,NBLK(IT)
C
               DO J = 1,NSOLBLK(IBLK,IT)
                  DO I = 1,NSOLBLK(IBLK,IT)
                     W(I,J) = C0
                  END DO
               END DO
C
               MODIFYV = .FALSE.
C
               IF ( ORBPOLSS.NE.'NONE      ' ) THEN
                  NSOL = NSOLBLK(IBLK,IT)
                  DO I = 1,NSOL
                     IKM = IKMSOLBLK(I,IBLK,IT)
                     L = L_IKM(IKM)
                     IF ( L.EQ.LOPT(IT) ) THEN
                        IF ( ORBPOLSS(1:5).EQ.'SIGMA' .OR. ORBPOLSS(1:4)
     &                       .EQ.'DMFT' .OR. ORBPOLSS(1:5).EQ.'LDA+U' )
     &                       MODIFYV = .TRUE.
                        EXIT
                     END IF
                  END DO
               END IF
C
               IF ( MODIFYV ) THEN
C
                  IF ( ORBPOLSS(1:4).EQ.'DMFT' .OR. ORBPOLSS(1:5)
     &                 .EQ.'LDA+U' ) THEN
C
                     LDAU = .FALSE.
                     IKMBOT = 2*LOPT(IT)**2 + 1
                     IKMTOP = IKMBOT - 1 + 2*(2*LOPT(IT)+1)
C
                     CALL CINIT(NCPLWFMAX*NCPLWFMAX,W)
C
                     DO J = 1,NSOL
                        JKM = IKMSOLBLK(J,IBLK,IT)
                        IF ( JKM.GE.IKMBOT .AND. JKM.LE.IKMTOP ) THEN
                           DO I = 1,NSOL
                              IKM = IKMSOLBLK(I,IBLK,IT)
                              IF ( IKM.GE.IKMBOT .AND. IKM.LE.IKMTOP )
     &                             THEN
                                 LDAU = .TRUE.
                                 W(I,J) = DMFTSIG(IKM,JKM,IT)
                              END IF
                           END DO
                        END IF
                     END DO
C
                     DO IR = 1,JRCRI(IM)
                        CVT(IR) = VT(IR,IT)
                        CBT(IR) = BT(IR,IT)
                     END DO
C
                     DO IPOT = 2,NFPT(IT)
                        LM = LMIFP(IPOT,IT)
                        DO IR = JRNS1(IM),JRCRI(IM)
                           CVNST(IR,IPOT) = VNST(IR,LM,IT)
                           CBNST(IR,IPOT) = BNST(IR,LM,IT)
                        END DO
                     END DO
                  ELSE IF ( ORBPOLSS(1:5).EQ.'SIGMA' ) THEN
                     DO IR = 1,JRCRI(IM)
                        CVT(IR) = VT(IR,IT) + SEVT(IR,IT)
                        CBT(IR) = BT(IR,IT) + SEBT(IR,IT)
                     END DO
                     DO LM = 2,NLMFPT(IT)
                        IF ( KLMFP(LM,IT).NE.0 ) THEN
                           DO IR = JRNS1(IM),JRCRI(IM)
                              CVNST(IR,LM) = VNST(IR,LM,IT)
     &                           + SEVNST(IR,LM,IT)
                              CBNST(IR,LM) = BNST(IR,LM,IT)
     &                           + SEBNST(IR,LM,IT)
                           END DO
                        END IF
                     END DO
                     LDAU = .FALSE.
                  ELSE
                     CALL STOP_MESSAGE(ROUTINE,
     &                                 'MODIFYV = T AND ORBPOL != SIGMA'
     &                                 )
C
                  END IF
C
                  IF ( SOLVER_FP.EQ.'BS' ) THEN
C
                     CALL FPDIRBS(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,
     &                            GETIRRSOL,C,ERYD,P,CVT,CBT,
     &                            CVNST(JRNSMIN,1),CBNST(JRNSMIN,1),
     &                            LDAU,W,Z(IT),NFPT(IT),JRNS1(IM),
     &                            R(1,IM),DRDI(1,IM),DRDIOVR(1,IM),
     &                            VAMEG(1,1,1,1,IBLK,IT),
     &                            VAMEF(1,1,1,1,IBLK,IT),IBLK,
     &                            ISOLIKM(1,IT),NSOLBLK(1,IT),
     &                            IKMSOLBLK(1,1,IT),JRCUT(0,IM),NPAN(IM)
     &                            ,PRM,QRM,PIM,QIM,DPM,DQM,PR,QR,DP,DQ,
     &                            NCPLWFMAX,NPANMAX,NKMMAX,NLMFPMAX,
     &                            JRNSMIN,NRMAX,GF_CONV_RH)
C
                  ELSE
C
                     CALL FPDIRBORN(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,IT,
     &                              GETIRRSOL,C,ERYD,P,CVT,CBT,
     &                              CVNST(JRNSMIN,1),CBNST(JRNSMIN,1),
     &                              LDAU,W,NFPT(IT),JRNS1(IM),DRDI(1,IM)
     &                              ,VAMEG(1,1,1,1,IBLK,IT),
     &                              VAMEF(1,1,1,1,IBLK,IT),IBLK,
     &                              ISOLIKM(1,IT),NSOLBLK(1,IT),
     &                              IKMSOLBLK(1,1,IT),JRCUT(0,IM),
     &                              NPAN(IM),PRM,QRM,PIM,QIM,PR,QR,DP,
     &                              DQ,TBLK,TBLK0,TSST0,SBLK,SBLK0,
     &                              SSST0,CBLK,DBLK,GBLK,NCWF_SOL_SPH,
     &                              IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,JG_SPH,
     &                              JF_SPH,NCPLWFMAX,NPANMAX,NKMMAX,
     &                              NLMFPMAX,JRNSMIN,NRMAX)
                  END IF
C
               ELSE
C
                  LDAU = .FALSE.
C
                  DO IR = 1,JRCRI(IM)
                     CVT(IR) = VT(IR,IT)
                     CBT(IR) = BT(IR,IT)
                  END DO
C
                  DO IPOT = 2,NFPT(IT)
                     LM = LMIFP(IPOT,IT)
                     DO IR = JRNS1(IM),JRCRI(IM)
                        CVNST(IR,IPOT) = VNST(IR,LM,IT)
                        CBNST(IR,IPOT) = BNST(IR,LM,IT)
                     END DO
                  END DO
C
                  IF ( SOLVER_FP.EQ.'BS' ) THEN
C
                     CALL FPDIRBS(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,
     &                            GETIRRSOL,C,ERYD,P,CVT,CBT,
     &                            CVNST(JRNSMIN,1),CBNST(JRNSMIN,1),
     &                            LDAU,W,Z(IT),NFPT(IT),JRNS1(IM),
     &                            R(1,IM),DRDI(1,IM),DRDIOVR(1,IM),
     &                            VAMEG(1,1,1,1,IBLK,IT),
     &                            VAMEF(1,1,1,1,IBLK,IT),IBLK,
     &                            ISOLIKM(1,IT),NSOLBLK(1,IT),
     &                            IKMSOLBLK(1,1,IT),JRCUT(0,IM),NPAN(IM)
     &                            ,PRM,QRM,PIM,QIM,DPM,DQM,PR,QR,DP,DQ,
     &                            NCPLWFMAX,NPANMAX,NKMMAX,NLMFPMAX,
     &                            JRNSMIN,NRMAX,GF_CONV_RH)
C
                  ELSE
C
                     CALL FPDIRBORN(LHS_SOL_EQ_RHS_SOL,I_HAND_SIDE,IT,
     &                              GETIRRSOL,C,ERYD,P,CVT,CBT,
     &                              CVNST(JRNSMIN,1),CBNST(JRNSMIN,1),
     &                              LDAU,W,NFPT(IT),JRNS1(IM),DRDI(1,IM)
     &                              ,VAMEG(1,1,1,1,IBLK,IT),
     &                              VAMEF(1,1,1,1,IBLK,IT),IBLK,
     &                              ISOLIKM(1,IT),NSOLBLK(1,IT),
     &                              IKMSOLBLK(1,1,IT),JRCUT(0,IM),
     &                              NPAN(IM),PRM,QRM,PIM,QIM,PR,QR,DP,
     &                              DQ,TBLK,TBLK0,TSST0,SBLK,SBLK0,
     &                              SSST0,CBLK,DBLK,GBLK,NCWF_SOL_SPH,
     &                              IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,JG_SPH,
     &                              JF_SPH,NCPLWFMAX,NPANMAX,NKMMAX,
     &                              NLMFPMAX,JRNSMIN,NRMAX)
                  END IF
C
               END IF
C
C-----------------------------------------------------------------------
C
               NSOL = NSOLBLK(IBLK,IT)
C
C------------------------------------------------------------------ BORN
C                                                store specific matrices
               IF ( SOLVER_FP.EQ.'BORN' ) THEN
                  DO J = 1,NSOL
                     JKM = IKMSOLBLK(J,IBLK,IT)
                     DO I = 1,NSOL
                        IKM = IKMSOLBLK(I,IBLK,IT)
                        TSST_BORN(IKM,JKM) = TBLK(I,J)
                        GAMMA(IKM,JKM) = GBLK(I,J)
                        CMAT(IKM,JKM) = CBLK(I,J)
                        DMAT(IKM,JKM) = DBLK(I,J)
                     END DO
                  END DO
Cc                dblk(:,:) = tblk(:,:)
               END IF
C------------------------------------------------------------------ BORN
C
C
C-----------------------------------------------------------------------
C        GENERATE T-MATRIX FROM CALCULATED WAVE FUNCTIONS
C-----------------------------------------------------------------------
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     IKM1 = IKMSOLBLK(I,IBLK,IT)
                     L1 = L_IKM(IKM1)
                     KAP1 = KAPPA_IKM(IKM1)
                     ISK1 = ISIGN(1,KAP1)
                     SK1 = DBLE(ISK1)
                     LB1 = L1 - ISK1
C
                     ABLK(I,J) = -(CHLM(L1)*QRM(I,J,IRTOP)-PFAC*SK1*CHLM
     &                           (LB1)*PRM(I,J,IRTOP))
C
                     TBLK(I,J) = (CHLP(L1)*QRM(I,J,IRTOP)-PFAC*SK1*CHLP(
     &                           LB1)*PRM(I,J,IRTOP))
C
                  END DO
               END DO
C
               CALL CMATINV3(NSOL,NCPLWFMAX,IPIVKM,TBLK,WKM2,BBLKINV)
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     BBLK(I,J) = ABLK(I,J) - TBLK(I,J)
                  END DO
               END DO
C
               CALL ZGEMM('N','N',NSOL,NSOL,NSOL,CI/(2*P),BBLK,
     &                    NCPLWFMAX,BBLKINV,NCPLWFMAX,C0,TBLK,NCPLWFMAX)
C
               DO J = 1,NSOL
                  IKM = IKMSOLBLK(J,IBLK,IT)
                  DO I = 1,NSOL
                     IKM1 = IKMSOLBLK(I,IBLK,IT)
                     TSST(IKM1,IKM,IT) = TBLK(I,J)
                  END DO
               END DO
C
Cc               if ( solver_fp.eq.'BORN' ) then
Cc                     tsst(:,:,it) = tsst_born(:,:)
Cc                     tblk(:,:) = dblk(:,:)
Cc               end if
C
C---------------------------- get MSST by matrix inversion MSST = 1/TSST
C
               CALL CMATINV3(NSOL,NCPLWFMAX,IPIVKM,TBLK,WKM2,BBLKINV)
C
               DO J = 1,NSOL
                  IKM = IKMSOLBLK(J,IBLK,IT)
                  DO I = 1,NSOL
                     IKM1 = IKMSOLBLK(I,IBLK,IT)
                     MSST(IKM1,IKM,IT) = BBLKINV(I,J)
                  END DO
               END DO
C
C------------------------------------------------------------------ BORN
               IF ( SOLVER_FP.EQ.'BORN      ' ) THEN
                  DO J = 1,NSOL
                     IKM = IKMSOLBLK(J,IBLK,IT)
                     NCPLWF(IKM) = NSOL
                     DO I = 1,NSOL
                        ICPL = I
                        IKM1 = IKMSOLBLK(I,IBLK,IT)
                        IKMCPLWF(ICPL,IKM) = IKM1
                        DO IR = 1,IRTOP
                           ZG(IR,ICPL,IKM) = PRM(I,J,IR)/R(IR,IM)
                           ZF(IR,ICPL,IKM) = QRM(I,J,IR)/(R(IR,IM)*C)
                           JG(IR,ICPL,IKM) = PIM(I,J,IR)/R(IR,IM)
                           JF(IR,ICPL,IKM) = QIM(I,J,IR)/(R(IR,IM)*C)
                        END DO
                     END DO
                  END DO
C------------------------------------------------------------------ BORN
C
               ELSE
C
C-------------------------------------------------------------- START BS
C
C-----------------------------------------------------------------------
C               transformation of wave functions
C-----------------------------------------------------------------------
C set up matrix of coefficients
C wave functions at r_crit:  BBLK - auxilary temporary solutions
C                            ABLK - scattering normalisation
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
                  DO I = 1,NSOL
                     IKM1 = IKMSOLBLK(I,IBLK,IT)
                     L1 = L_IKM(IKM1)
                     CJL = 0.5D0*(CHLP(L1)+CHLM(L1))
                     CHL = CHLP(L1)
                     DO J = 1,NSOL
                        BBLK(I,J) = PRM(I,J,IRTOP)/RTOP
                        ABLK(I,J) = CJL*BBLKINV(I,J)
                        IF ( J.EQ.I ) ABLK(I,J) = ABLK(I,J) - CI*P*CHL
                     END DO
                  END DO
C
C ----------------------------------------------------------------------
C                                                basis functions R and H
                  IF ( GF_CONV_RH ) THEN
                     HBLK(1:NSOL,1:NSOL)
     &                  = MATMUL(ABLK(1:NSOL,1:NSOL),TBLK(1:NSOL,1:NSOL)
     &                  )
                     ABLK(1:NSOL,1:NSOL) = HBLK(1:NSOL,1:NSOL)
                  END IF
C
C ----------------- get the transformation matrix X with BBLK * X = ABLK
C     ZGETRF: LU-decomposition of BBLK   ZGETRS: gives X in ABLK on exit
C
                  CALL ZGETRF(NSOL,NSOL,BBLK,NCPLWFMAX,IPIV,INFO)
                  CALL ZGETRS('N',NSOL,NSOL,BBLK,NCPLWFMAX,IPIV,ABLK,
     &                        NCPLWFMAX,INFO)
C
                  DO J = 1,NSOL
                     IKM = IKMSOLBLK(J,IBLK,IT)
                     NCPLWF(IKM) = NSOL
                     DO I = 1,NSOL
                        ICPL = I
                        IKM1 = IKMSOLBLK(I,IBLK,IT)
                        IKMCPLWF(ICPL,IKM) = IKM1
                        DO IR = 1,IRTOP
                           DO J1 = 1,NSOL
                              ZG(IR,ICPL,IKM) = ZG(IR,ICPL,IKM)
     &                           + PRM(I,J1,IR)*ABLK(J1,J)/R(IR,IM)
                              ZF(IR,ICPL,IKM) = ZF(IR,ICPL,IKM)
     &                           + QRM(I,J1,IR)*ABLK(J1,J)/(R(IR,IM)*C)
                           END DO
                           JG(IR,ICPL,IKM) = PIM(I,J,IR)/R(IR,IM)
                           JF(IR,ICPL,IKM) = QIM(I,J,IR)/(R(IR,IM)*C)
                        END DO
                     END DO
                  END DO
               END IF
C---------------------------------------------------------------- END BS
C
C --------------------- set alpha matrix SSST required for Lloyd formula
C
               DO J = 1,NSOL
                  JKM = IKMSOLBLK(J,IBLK,IT)
                  LJ = L_IKM(JKM)
                  ARG = P*R(1,IM)
                  CJL_J = CJLZ(LJ,ARG)
C
                  DO I = 1,NSOL
                     IKM = IKMSOLBLK(I,IBLK,IT)
C
                     SSST(IKM,JKM,IT) = ZG(1,I,JKM)/CJL_J
C
                  END DO
               END DO
C
            END DO
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C                  test wave functions:  WRONSKI
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
            IF ( WRONSKI ) THEN
               ZGL(:,:,:) = ZG(:,:,:)
               ZFL(:,:,:) = ZF(:,:,:)
               JGL(:,:,:) = JG(:,:,:)
               JFL(:,:,:) = JF(:,:,:)
C
               TSST_LHS(:,:) = TSST(:,:,IT)
            END IF
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
C
C-----------------------------------------------------------------------
C                        LEFT HAND SIDE solution
C-----------------------------------------------------------------------
            IF ( .NOT.LHS_SOL_EQ_RHS_SOL .AND. I_HAND_SIDE.EQ.1 .AND. 
     &           .NOT.CALC_GFWF ) THEN
C
               ZGL(:,:,:) = ZG(:,:,:)
               ZFL(:,:,:) = ZF(:,:,:)
               JGL(:,:,:) = JG(:,:,:)
               JFL(:,:,:) = JF(:,:,:)
C
               TSST_LHS(:,:) = TSST(:,:,IT)
C
C------------------------------------------------------------------ BORN
               IF ( SOLVER_FP.EQ.'BORN' ) THEN
                  GAMMA_LHS(:,:) = GAMMA(:,:)
                  CMAT_LHS(:,:) = CMAT(:,:)
                  DMAT_LHS(:,:) = DMAT(:,:)
               END IF
C------------------------------------------------------------------ BORN
C
               IFILSS_XHS = IFILCBWF_LHS
C
            ELSE
C-----------------------------------------------------------------------
C                        RIGHT HAND SIDE solution
C-----------------------------------------------------------------------
C
               IFILSS_XHS = IFILSS
C
            END IF
C
Cc      if ( solver_fp.eq.'BS' ) then
Cc         write(200+i_hand_side,*) zg,zf,jg,jf,
Cc     &     tsst(:,:,it),msst(:,:,it),ssst(:,:,it)
Cc      else
Cc         read(200+i_hand_side,*) zg,zf,jg,jf,
Cc     &     tsst(:,:,it),msst(:,:,it),ssst(:,:,it)
Cc
Cc            if (  i_hand_side.eq.1  ) then
Cc               zg_lhs(:,:,:) = zg(:,:,:)
Cc               zf_lhs(:,:,:) = zf(:,:,:)
Cc               jg_lhs(:,:,:) = jg(:,:,:)
Cc               jf_lhs(:,:,:) = jf(:,:,:)
Cc
Cc               tsst_lhs(:,:) = tsst(:,:,it)
Cc            end if
Cc
Cc      end if
C
C=======================================================================
C                       write matrix elements
C=======================================================================
C
            IF ( I_HAND_SIDE.EQ.N_HAND_SIDE ) THEN
C
C-----------------------------------------------------------------------
C    SOLVER_FP = BORN:  write information for spherical solutions
C-----------------------------------------------------------------------
C------------------------------------------------------------ START BORN
               IF ( SOLVER_FP.EQ.'BORN' .AND. IWRREGWF.NE.0 .AND. 
     &              .NOT.CALC_GFWF ) THEN
C
                  IF ( LHS_SOL_EQ_RHS_SOL ) THEN
                     GAMMA_LHS(:,:) = GAMMA(:,:)
                     CMAT_LHS(:,:) = CMAT(:,:)
                     DMAT_LHS(:,:) = DMAT(:,:)
                  END IF
C
                  WKM1(1:NKM,1:NKM) = TRANSPOSE(CMAT_LHS(1:NKM,1:NKM))
C
                  CHI(1:NKM,1:NKM) = MATMUL(GAMMA(1:NKM,1:NKM),WKM1(1:
     &                               NKM,1:NKM))
C
                  IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
                     ZETA(1:NKM,1:NKM) = C0
                     DO I = 1,NKM
                        ZETA(I,I) = 1D0
                     END DO
C
                  ELSE
C
                     WKM1(1:NKM,1:NKM)
     &                  = TRANSPOSE(DMAT_LHS(1:NKM,1:NKM))
                     DO I = 1,NKM
                        WKM1(I,I) = 1D0 + WKM1(I,I)
                     END DO
C
                     ZETA(1:NKM,1:NKM)
     &                  = MATMUL(GAMMA(1:NKM,1:NKM),WKM1(1:NKM,1:NKM))
C
                  END IF
C
                  IF ( IPRINT.GT.3 )
     &                  CALL CMATSTRUCT('GAMMA * (1+D^x)^T = 1  ?',ZETA,
     &                 NKM,NKMMAX,3,3,1,TOL,IOTMP1)
C
C- - - - - - - - - - - - - - force the zeta-matrix to be the unit matrix
                  ZETA(1:NKM,1:NKM) = C0
                  DO I = 1,NKM
                     ZETA(I,I) = 1D0
                  END DO
C
C- - - - - - - - - - - - provide properly normalized spherical solutions
                  DO J = 1,NKM
                     IF ( NCWF_SOL_SPH(J).GT.2 )
     &                    CALL STOP_MESSAGE(ROUTINE,'NCWF_SOL_SPH>2')
                     DO I = 1,NCWF_SOL_SPH(J)
                        DO IR = 1,IRTOP
                           ZG_SPH(IR,I,J) = ZG_SPH(IR,I,J)/R(IR,IM)
                           JG_SPH(IR,I,J) = JG_SPH(IR,I,J)/R(IR,IM)
                           ZF_SPH(IR,I,J) = ZF_SPH(IR,I,J)/(R(IR,IM)*C)
                           JF_SPH(IR,I,J) = JF_SPH(IR,I,J)/(R(IR,IM)*C)
                        END DO
                     END DO
                  END DO
C
                  WRITE (IFILCBWF_SPH,REC=IT)
     &                   ((CHI(I,J),I=1,NKM),J=1,NKM),
     &                   ((ZETA(I,J),I=1,NKM),J=1,NKM),
     &                   ((GAMMA(I,J),I=1,NKM),J=1,NKM),
     &                   ((GAMMA_LHS(I,J),I=1,NKM),J=1,NKM),
     &                   (NCWF_SOL_SPH(I),I=1,NKM),
     &                   ((IKM_CWFSOL_SPH(I,J),I=1,2),J=1,NKM),
     &                   (((ZG_SPH(IR,I,J),ZF_SPH(IR,I,J),JG_SPH(IR,I,J)
     &                   ,JF_SPH(IR,I,J),IR=1,IRTOP),I=1,2),J=1,NKM),
     &                   IRTOP
C
               END IF
C--------------------------------------------------- spherical solutions
C-------------------------------------------------------------- END BORN
C
               IF ( LHS_SOL_EQ_RHS_SOL .OR. CALC_GFWF ) THEN
C
                  ZGL(:,:,:) = ZG(:,:,:)
                  ZFL(:,:,:) = ZF(:,:,:)
                  JGL(:,:,:) = JG(:,:,:)
                  JFL(:,:,:) = JF(:,:,:)
C
               END IF
C
               CALL FPMATELM(IT,IRTOP,ZGL,ZFL,JGL,JFL,ZG,ZF,JG,JF,CHI,
     &                       NCWF_SOL_SPH,IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,
     &                       JG_SPH,JF_SPH,MEZZ,MEZJ,MEZZL,MEZJL)
C
            END IF
C
C ------------------------------------------- write wavefunctions for IT
C
            DO ISOL = 1,NKM_T(IT)
C
               IF ( IWRREGWF.NE.0 )
     &              WRITE (IFILSS_XHS,REC=ISOL+(IT-1)*NKM) IT,'REG',
     &              ISOL,IRTOP,NCPLWF(ISOL),
     &              (IKMCPLWF(K,ISOL),K=1,NCPLWF(ISOL)),
     &              ((ZG(I,K,ISOL),I=1,IRTOP),K=1,NCPLWF(ISOL)),
     &              ((ZF(I,K,ISOL),I=1,IRTOP),K=1,NCPLWF(ISOL))
C
               IF ( IWRIRRWF.NE.0 )
     &              WRITE (IFILSS_XHS,REC=ISOL+(IT-1+NT)*NKM) IT,'IRR',
     &              ISOL,IRTOP,
     &              ((JG(I,K,ISOL),I=1,IRTOP),K=1,NCPLWF(ISOL)),
     &              ((JF(I,K,ISOL),I=1,IRTOP),K=1,NCPLWF(ISOL))
C
            END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                  CHECK
            IF ( IPRINT.GT.3 ) THEN
C
               WRITE (6,99003) 'atom type   ',IT,ERYD,CTL(IT,1)
               CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,3,3,1,
     &                         1.0D-9,6)
C------------------------------------------------------------------ BORN
               IF ( SOLVER_FP.EQ.'BORN' )
     &              CALL CMATSTRUCT('T-BORN  ',TSST_BORN(1,1),NKM,
     &              NKMMAX,3,3,1,1.0D-9,6)
C------------------------------------------------------------------ BORN
               CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,
     &                         3,3,1,1.0D-9,6)
               CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,
     &                         3,3,1,1.0D-9,6)
            END IF
C
            IF ( DUMP_SSITE ) THEN
               IF ( IT.EQ.1 .AND. I_HAND_SIDE.EQ.1 ) THEN
                  LFILNAM = LEN_TRIM(SOLVER_FP)
                  FILNAM = 'zzzzzz_ssite_REL_'//SOLVER_FP(1:LFILNAM)
                  LFILNAM = LFILNAM + 17
                  IF ( GF_CONV_RH ) THEN
                     FILNAM = FILNAM(1:LFILNAM)//'_RH.dat'
                  ELSE
                     FILNAM = FILNAM(1:LFILNAM)//'_ZJ.dat'
                  END IF
                  CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,
     &                                 FILNAM(1:(LFILNAM+7)))
               END IF
               WRITE (IOTMP,99003) 'atom type   ',IT,ERYD,CTL(IT,1),
     &                             SOLVER_FP,GF_CONV_RH
               IFMT = 1
               CALL CMATSTRUCT('T-MATRIX      (kappa,mu)',TSST(1,1,IT),
     &                         NKM,NKMMAX,3,3,1,TOL,IOTMP)
               CALL CHANGEREP(NKM,NKMMAX,TSST(1,1,IT),'REL>RLM',WKM2)
               CALL CMATSTRUCT('T-MATRIX       (l,ml,ms)',WKM2,IFMT*NLM,
     &                         NKMMAX,IFMT,IFMT,1,TOL,IOTMP)
               DO I = 1,IFMT*NLM
                  DO J = 1,I
                     WRITE (IOTMP,'(2I3,2F20.12,2X,A)') I,J,WKM2(I,J),
     &                      SOLVER_FP
                     WRITE (IOTMP,'(2I3,2F20.12,2X,A)') J,I,WKM2(J,I),
     &                      SOLVER_FP
                  END DO
               END DO
               CALL CMATSTRUCT('MEZZ-MATRIX   (kappa,mu)',MEZZ(1,1,IT,1)
     &                         ,NKM,NKMMAX,3,3,1,TOL,IOTMP)
               CALL CHANGEREP(NKM,NKMMAX,MEZZ(1,1,IT,1),'REL>RLM',WKM2)
               CALL CMATSTRUCT('MEZZ-MATRIX    (l,ml,ms)',WKM2,IFMT*NLM,
     &                         NKMMAX,IFMT,IFMT,1,TOL,IOTMP)
               CALL CMATSTRUCT('MEZJ-MATRIX    (kappa,mu)',
     &                         MEZJ(1,1,IT,1),NKM,NKMMAX,3,3,1,TOL,
     &                         IOTMP)
               CALL CHANGEREP(NKM,NKMMAX,MEZJ(1,1,IT,1),'REL>RLM',WKM2)
               CALL CMATSTRUCT('MEZJ-MATRIX    (l,ml,ms)',WKM2,IFMT*NLM,
     &                         NKMMAX,IFMT,IFMT,1,TOL,IOTMP)
               DO I = 1,NKM
                  DO J = 1,NKM
                     WRITE (IOTMP,99001) I,J,TSST(I,J,IT),MEZZ(I,J,IT,1)
     &                      ,MEZJ(I,J,IT,1)
                  END DO
               END DO
C
               IF ( .NOT.LHS_SOL_EQ_RHS_SOL .AND. I_HAND_SIDE.EQ.2 )
     &              THEN
C
C------------------------------------------------------------------ BORN
                  IF ( SOLVER_FP.EQ.'BORN' ) THEN
C
                     WKM1(1:NKM,1:NKM)
     &                  = TRANSPOSE(DMAT_LHS(1:NKM,1:NKM))
                     DO I = 1,NKM
                        WKM1(I,I) = 1D0 + WKM1(I,I)
                     END DO
                     CALL CMATMUL(NKM,NKMMAX,GAMMA,WKM1,WKM2)
C
                     CALL CMATSTRUCT('GAMMA * (1+D^x)^T = 1  ?',WKM2,
     &                               NKM,NKMMAX,3,3,1,TOL,IOTMP)
                  END IF
C------------------------------------------------------------------ BORN
C
                  WKM1(1:NKM,1:NKM) = TRANSPOSE(TSST_LHS(1:NKM,1:NKM))
                  WKM2(1:NKM,1:NKM) = WKM1(1:NKM,1:NKM)
     &                                - TSST(1:NKM,1:NKM,IT)
C
                  CALL CMATSTRUCT('t-matrix  t^xT - t = 0 ?',WKM2,NKM,
     &                            NKMMAX,3,3,1,TOL,IOTMP)
C
               END IF
C
               IF ( IT.EQ.ITTOP .AND. I_HAND_SIDE.EQ.N_HAND_SIDE ) THEN
                  CLOSE (IOTMP)
                  WRITE (6,99004) FILNAM(1:(LFILNAM+7))
               END IF
            END IF
C                                                                  CHECK
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
         END DO
C HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
C
C-----------------------------------------------------------------------
C               store l-projected matrix elements
C-----------------------------------------------------------------------
         IF ( L_PROJECTED_ME .AND. .NOT.CALC_GFWF ) WRITE (IFILMEZZL)
     &        IT,MEZZL,MEZJL
C-----------------------------------------------------------------------
C
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C                  test wave functions:  WRONSKI
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
         IF ( WRONSKI ) THEN
C
            WRITE (6,99006) 'atom type   ',IT,ERYD,CTL(IT,1),SOLVER_FP,
     &                      GF_CONV_RH
C
            D_WRON(:,:,IT) = 0D0
C
            DO IBLK2 = 1,NBLK(IT)
               DO J2 = 1,NSOL
                  JKM2 = IKMSOLBLK(J2,IBLK2,IT)
C
                  DO IBLK1 = 1,NBLK(IT)
                     IF ( IBLK1.NE.IBLK2 ) CYCLE
                     DO J1 = 1,NSOL
                        JKM1 = IKMSOLBLK(J1,IBLK1,IT)
C
                        IF ( JKM1.EQ.JKM2 ) THEN
                           CREF = 1D0
                        ELSE
                           CREF = 0D0
                        END IF
                        D_WRON(JKM1,JKM2,IT) = CREF
C
                        DO IR = 1,IRTOP
                           CRSQ = WFAC*R(IR,IM)**2
                           CSUM = C0
                           DO I = 1,NSOL
                              CSUM = CSUM + 
     &                               (ZF(IR,I,JKM1)*JGL(IR,I,JKM2)
     &                               -ZG(IR,I,JKM1)*JFL(IR,I,JKM2))*CRSQ
                           END DO
                           IF ( IPRINT.GE.1 ) WRITE (6,99002) IT,IBLK1,
     &                          IBLK2,J1,JKM1,J2,JKM2,IR,R(IR,IM),CSUM
C
                           IF ( ABS(CREF-D_WRON(JKM1,JKM2,IT))
     &                          .LT.ABS(CREF-CSUM) )
     &                          D_WRON(JKM1,JKM2,IT) = CSUM
C
                        END DO
                        WRITE (6,*) ' '
                     END DO
                  END DO
               END DO
            END DO
         END IF
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C         rotate ALL matrices from the LOCAL to the GLOBAL frame
C         starting with version 7.4.0 MEZZ and MEZJ are excluded
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      IF ( KMROT.NE.0 ) CALL SSITE_ROTATE(TSST,MSST,SSST)
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C                 symetrize ALL single site matrices
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
CC    CALL SSITE_SYMMETRIZE('TSST',TSST)
CC    CALL SSITE_SYMMETRIZE('MSST',MSST)
CC    CALL SSITE_SYMMETRIZE('SSST',SSST)
C
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C                  test wave functions:  WRONSKI
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      IF ( WRONSKI ) THEN
         DO IT = ITBOT,ITTOP
            DO JKM2 = 1,NKM_T(IT)
               DO JKM1 = 1,NKM_T(IT)
                  IF ( JKM1.EQ.JKM2 ) THEN
                     CREF = 1D0
                  ELSE
                     CREF = 0D0
                  END IF
                  IF ( ABS(CREF-D_WRON(JKM1,JKM2,IT)).GT.2D-12 )
     &                 WRITE (6,99007) IT,JKM1,JKM2,D_WRON(JKM1,JKM2,IT)
               END DO
            END DO
         END DO
         IF ( WRONSKI ) CALL STOP_REGULAR(ROUTINE,'WRONSKIAN tabulated')
      END IF
Cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
      IF ( IPRINT.GT.3 ) THEN
         WRITE (6,99005)
         CLOSE (IOTMP1)
      END IF
      IF ( IPRINT.GT.3 ) CLOSE (IOTMP1)
C
      DEALLOCATE (W,DP,CBT,ABLK,BBLK,CVT,BBLKINV,IPIV)
      DEALLOCATE (PIM,PRM,QIM,QRM,CBNST,CVNST)
      IF ( SOLVER_FP.EQ.'BS' ) DEALLOCATE (DPM,DQM)
C
C??????????????????????????????????????????????????????????????????????
      I = -100
      IF ( I.GT.0 ) THEN
         WRITE (I,'(''MEZZ'')')
         WRITE (I,'(2E25.14)') MEZZ(:,:,1,1:3)
         WRITE (I,'(''MEZJ'')')
         WRITE (I,'(2E25.14)') MEZJ(:,:,1,1:3)
         STOP
      END IF
C??????????????????????????????????????????????????????????????????????
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4,5,6
      IF ( BREAKPOINT.EQ.4 ) CALL STOP_BREAKPOINT(ROUTINE)
      IF ( BREAKPOINT.EQ.5 .OR. BREAKPOINT.EQ.6 )
     &     CALL BREAKPOINT_5_6_SSITE(ERYD,TSST,MSST,SSST,MEZZ,MEZJ)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4,5,6
C
C=======================================================================
99001 FORMAT (2I3,2X,'TSST ',2E25.15,/,8X,'MEZZ ',2E25.15,/,8X,'MEZJ ',
     &        2E25.15,/)
99002 FORMAT (' IT=',I2,'  IBLK=',I2,I4,'  J1=',I2,I4,'  J2=',I2,I4,I7,
     &        F12.8,(2X,2F18.12))
99003 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),//,10X,'energy',12X,
     &        2F12.6,/,10X,'speed of light',f16.6,/,:,10X,
     &        'solver            ',A,/,10X,'GF_CONV_RH        ',L10,/)
99004 FORMAT (//,5X,'>>>>>  DUMP written to file  ',A,//)
99005 FORMAT (//,5X,'test results for GAMMA * (1+D^x)^T = 1  ?',/,5X,
     &        'written ti file:  zzzzz_check-matrices',/)
99006 FORMAT (//,1X,79('W'),/,10X,A,I4,/,1X,79('W'),//,10X,'energy',12X,
     &        2F12.6,/,10X,'speed of light',f16.6,/,:,10X,
     &        'solver            ',A,/,10X,'GF_CONV_RH        ',L10,/,
     &        10X,'check WRONSKIAN',//)
99007 FORMAT (I10,2I3,2F18.12)
      END
C*==fpmatelm.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE FPMATELM(IT,IRTOP,ZGL,ZFL,JGL,JFL,ZGR,ZFR,JGR,JFR,CHI,
     &                    NCWF_SOL_SPH,IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,
     &                    JG_SPH,JF_SPH,MEZZ,MEZJ,MEZZL,MEZJL)
C   ********************************************************************
C   *                                                                  *
C   *  calculate all necessary matrix elements                         *
C   *               MEZZ =  < Z_LAM | A | Z_LAM' >                     *
C   *               MEZJ =  < Z_LAM | A | J_LAM' >                     *
C   *  with A = 1, sig_z, l_z, H_hf  (index IME)                       *
C   *                                                                  *
C   *  L_PROJECTED_ME = T:  calculate the l-diagonal part MEZZL, MEZJL *
C   *                     to be used in the FP-DOS calculation         *
C   *                                                                  *
C   *  NOTE: if  LHS_SOL_EQ_RHS_SOL = .FALSE.                          *
C   *        ZGL,ZFL,JGL,JFL correspond to the LHS solution            *
C   *        ZGR,ZFR,JGR,JFR correspond to the RHS solution            *
C   *                                                                  *
C   *        otherwise LHS = RHS solution with ZGL = ZGR etc           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IHYPER,SOLVER_FP,NONMAG,L_PROJECTED_ME,
     &    BORN_USE_SPH
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IMT,IKMCPLWF,NLT,NKM_T,Z
      USE MOD_RMESH,ONLY:R,DRDI,JRCUT,FLMSF,LMISF,JRNS1,NRMAX,NSF,
     &    NSFMAX,W_RADINT
      USE MOD_ANGMOM,ONLY:NMEMAX,NLMAX,NKMMAX,IMKM_IKM,AME_G,NCPLWF,NL,
     &    AG_RGNT,NLM_AME_RLM_EXT,IDOS
      USE MOD_CONSTANTS,ONLY:C1,C0,SQRT_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPMATELM')
      INTEGER NOBS
      PARAMETER (NOBS=3)
      LOGICAL CHECK_OBS_ME
      PARAMETER (CHECK_OBS_ME=.FALSE.)
C
C Dummy arguments
C
      INTEGER IRTOP,IT
      COMPLEX*16 CHI(NKMMAX,NKMMAX),JFL(NRMAX,NCPLWFMAX,NKMMAX),
     &           JFR(NRMAX,NCPLWFMAX,NKMMAX),JF_SPH(NRMAX,2,NKMMAX),
     &           JGL(NRMAX,NCPLWFMAX,NKMMAX),JGR(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JG_SPH(NRMAX,2,NKMMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           ZFL(NRMAX,NCPLWFMAX,NKMMAX),ZFR(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZF_SPH(NRMAX,2,NKMMAX),ZGL(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGR(NRMAX,NCPLWFMAX,NKMMAX),ZG_SPH(NRMAX,2,NKMMAX)
      INTEGER IKM_CWFSOL_SPH(2,NKMMAX),NCWF_SOL_SPH(NKMMAX)
C
C Local variables
C
      REAL*8 ABS_DELTA,RSCAL(:,:),RSCALHFF(:,:),RSCALHFFNUC(:,:),TOL,
     &       WINTR(:,:)
      LOGICAL ADD_TO_MEZZL,CPLMAT(:,:,:),CPLMES(:,:),NON0_MATELM
      COMPLEX*16 CSCL,CSUM,MEZJ_SPHER,MZBJA(:,:,:,:),MZBZA(:,:,:,:),
     &           RHOC(:,:),TMPINT(:,:),ZFJF_SPH(:),ZFZF_SPH(:),
     &           ZGJG_SPH(:),ZGZG_SPH(:)
      INTEGER IA_ERR,ICPLWF1,ICPLWF2,IKM,IKM1,IKM2,IKM3,IKM4,
     &        IKM_MATELM(:,:),IL,IM,IMATELM,IME,IR,IRBOT,IRNSBOT,ISF,J,
     &        JKM1,JKM2,LM,MIKM3,MIKM4,N,NKAME,NKMAME,NMATELM(:),NMETOP,
     &        NOBSE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RSCALHFF,RSCALHFFNUC,RSCAL,WINTR
      ALLOCATABLE CPLMAT,CPLMES
      ALLOCATABLE TMPINT,IKM_MATELM,NMATELM
      ALLOCATABLE ZFJF_SPH,ZFZF_SPH,ZGJG_SPH,ZGZG_SPH,RHOC
      ALLOCATABLE MZBJA,MZBZA
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (IKM_MATELM(NKMMAX,NKMMAX),NMATELM(NKMMAX))
      ALLOCATE (RSCALHFF(NRMAX,NSFMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RSCALHFF')
      ALLOCATE (RSCALHFFNUC(NRMAX,NSFMAX),RSCAL(NRMAX,NSFMAX))
      ALLOCATE (WINTR(NRMAX,NSFMAX))
      ALLOCATE (TMPINT(NCPLWFMAX,NCPLWFMAX))
C
      IM = IMT(IT)
C
C------------------------------------------------------------------BORN
      IF ( SOLVER_FP.EQ.'BORN' ) THEN
         ALLOCATE (ZGJG_SPH(NRMAX),ZGZG_SPH(NRMAX))
         ALLOCATE (ZFJF_SPH(NRMAX),ZFZF_SPH(NRMAX))
         ALLOCATE (RHOC(NRMAX,3))
         IRNSBOT = JRNS1(IM)
C------------------------------------------------------------------ BORN
      ELSE
C-------------------------------------------------------------------- BS
         ALLOCATE (ZGJG_SPH(1),ZGZG_SPH(1))
         ALLOCATE (ZFJF_SPH(1),ZFZF_SPH(1))
         ALLOCATE (RHOC(1,1))
         IRNSBOT = -999999
      END IF
C-------------------------------------------------------------------- BS
C
C------------------------------------------- NKAME should be > NK=2*NL-1
      NKAME = 2*NLT(IT) + 1
      NKMAME = 2*((NKAME+1)/2)**2
C
      IF ( NKM_T(IT).GT.NKMAME )
     &      CALL STOP_MESSAGE(ROUTINE,'NKM_T(IT) > NKAME')
C
      ALLOCATE (CPLMES(NKMAME,NKMAME))
      ALLOCATE (CPLMAT(NKMAME,NKMAME,NSFMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CPLMES')
C
C-----------------------------------------------------------------------
      IF ( L_PROJECTED_ME ) THEN
         MEZZL(:,:,:,:) = C0
         MEZJL(:,:,:,:) = C0
      END IF
C-----------------------------------------------------------------------
C
      CALL RINIT(NRMAX*NSFMAX,RSCAL)
      CALL RINIT(NRMAX*NSFMAX,RSCALHFFNUC)
      CALL RINIT(NRMAX*NSFMAX,RSCALHFF)
C
      RHOC(:,:) = C0
C
      IF ( NONMAG ) THEN
         NMETOP = 1
         NOBSE = 1
      ELSE
         NMETOP = 3
         NOBSE = 4
      END IF
C
C---- consider all possible couplings due to the angular matrix elements
C----------------- of the various operators ignoring the shape functions
C
      CPLMES(:,:) = .FALSE.
      DO IME = 1,NOBSE
         DO IKM2 = 1,NKMAME
            DO IKM1 = 1,NKMAME
               IF ( ABS(AME_G(IKM1,IKM2,2,IME)).GT.1D-8 )
     &              CPLMES(IKM1,IKM2) = .TRUE.
            END DO
         END DO
      END DO
C
      CPLMAT(:,:,:) = .FALSE.
C
C-------- fold coupling due to operator angular matrix elements (CPLMAT)
C---------------------------- with that due to shape functions (AG_RGNT)
C
      DO ISF = 1,NSF(IM)
C
         LM = LMISF(ISF,IM)
C
C-----------------------------------------------------------------------
         IF ( LM.GT.NLM_AME_RLM_EXT ) CYCLE
C-----------------------------------------------------------------------
C
         DO IKM1 = 1,NKMAME
            DO IKM3 = 1,NKMAME
               IF ( ABS(AG_RGNT(IKM1,IKM3,LM)).GT.1D-12 ) THEN
                  DO IKM2 = 1,NKMAME
C
                     CPLMAT(IKM1,IKM2,ISF) = CPLMAT(IKM1,IKM2,ISF) .OR. 
     &                  CPLMES(IKM3,IKM2)
C
                  END DO
               END IF
            END DO
         END DO
C
      END DO
C
C---------------------------- find possible nonvanishing matrix elements
C
      DO IKM1 = 1,NKM_T(IT)
         NMATELM(IKM1) = 0
         DO IKM2 = 1,NKM_T(IT)
C
            NON0_MATELM = .FALSE.
C
            DO ICPLWF1 = 1,NCPLWF(IKM1)
               IKM3 = IKMCPLWF(ICPLWF1,IKM1)
               MIKM3 = IMKM_IKM(IKM3)
               DO ICPLWF2 = 1,NCPLWF(IKM2)
                  IKM4 = IKMCPLWF(ICPLWF2,IKM2)
                  MIKM4 = IMKM_IKM(IKM4)
                  DO ISF = 1,NSF(IM)
C
                     IF ( CPLMAT(IKM3,IKM4,ISF) .OR. 
     &                    CPLMAT(MIKM3,MIKM4,ISF) ) NON0_MATELM = .TRUE.
                  END DO
               END DO
            END DO
C
            IF ( NON0_MATELM ) THEN
               NMATELM(IKM1) = NMATELM(IKM1) + 1
               IKM_MATELM(NMATELM(IKM1),IKM1) = IKM2
            END IF
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C          set up scaling factor for r-integration
C-----------------------------------------------------------------------
C
C for 0th and 1st panel, only shape function for lm = 1 (l = 0,ml = 0)
C with FLMSF(IR,ISF=1) = SQRT(4 PI)
C
      DO ISF = 1,NSF(IM)
C
         LM = LMISF(ISF,IM)
C
         IF ( LM.EQ.1 ) THEN
            DO J = 1,JRCUT(1,IM)
C
               RSCAL(J,ISF) = R(J,IM)*R(J,IM)*DRDI(J,IM)*SQRT_4PI
               RSCALHFF(J,ISF) = DRDI(J,IM)*SQRT_4PI
C
            END DO
         ELSE
            DO J = 1,JRCUT(1,IM)
C
               RSCAL(J,ISF) = 0D0
               RSCALHFF(J,ISF) = 0D0
C
            END DO
         END IF
C
C for remaining panels, all shape functions come into play
C
         DO J = JRCUT(1,IM) + 1,IRTOP
C
            RSCAL(J,ISF) = R(J,IM)*R(J,IM)*DRDI(J,IM)
     &                     *FLMSF(J-JRCUT(1,IM),ISF,IM)
C
            RSCALHFF(J,ISF) = DRDI(J,IM)*FLMSF(J-JRCUT(1,IM),ISF,IM)
C
C              RSCALHFFNUC(J,ISF) = (R(J,im)/RNUC)**3*DRDI(J,im)
C     &                            *FLMSF(J-JRCUT(1,im),ISF,im)
C
         END DO
C
         CALL RVECWGT(IRTOP,RSCAL(1,ISF),W_RADINT(1,IM),WINTR(1,ISF))
C
      END DO
C
C=======================================================================
C                                                              JKM1 JKM2
      DO JKM1 = 1,NKM_T(IT)
         DO IMATELM = 1,NMATELM(JKM1)
            JKM2 = IKM_MATELM(IMATELM,JKM1)
C
C------------------------------------------------------------- ISF START
            DO ISF = 1,NSF(IM)
C
               LM = LMISF(ISF,IM)
C
C------------- skip loop if there is no coupling JKM1-JKM2 for given ISF
C
               DO ICPLWF1 = 1,NCPLWF(JKM1)
                  IKM3 = IKMCPLWF(ICPLWF1,JKM1)
                  MIKM3 = IMKM_IKM(IKM3)
C
                  DO ICPLWF2 = 1,NCPLWF(JKM2)
                     IKM4 = IKMCPLWF(ICPLWF2,JKM2)
                     MIKM4 = IMKM_IKM(IKM4)
C
                     IF ( CPLMAT(IKM3,IKM4,ISF) .OR. 
     &                    CPLMAT(MIKM3,MIKM4,ISF) ) GOTO 10
C
                  END DO
               END DO
C
               CYCLE
C
C-----------------------------------------------------------------------
 10            CONTINUE
               IF ( ISF.EQ.1 ) THEN
                  IRBOT = 1
               ELSE
                  IRBOT = JRCUT(1,IM) + 1
               END IF
C
               ADD_TO_MEZZL = L_PROJECTED_ME .AND. ISF.EQ.1
C
               CALL CINTABRSUM(ADD_TO_MEZZL,IT,JKM1,JKM2,LM,
     &                         ZGL(1,1,JKM1),ZGR(1,1,JKM2),JGL(1,1,JKM1)
     &                         ,JGR(1,1,JKM2),ZFL(1,1,JKM1),
     &                         ZFR(1,1,JKM2),JFL(1,1,JKM1),JFR(1,1,JKM2)
     &                         ,RHOC,WINTR(1,ISF),MEZZ,MEZJ,MEZZL,MEZJL,
     &                         NCPLWF(JKM1),NCPLWF(JKM2),IRBOT,IRTOP,
     &                         IRNSBOT)
C
C------------------------------------------------------------------ BORN
               IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH .AND. 
     &              ISF.EQ.1 ) CALL CINTABRSUM_SPHER(JKM1,JKM2,LM,CHI,
     &              NCWF_SOL_SPH,IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,JG_SPH,
     &              JF_SPH,ZFJF_SPH,ZFZF_SPH,ZGJG_SPH,ZGZG_SPH,RHOC,
     &              WINTR,IRNSBOT,NMETOP)
C------------------------------------------------------------------ BORN
C
C------------------------------------------------------- HYPERFINE FIELD
C
               IF ( IHYPER.NE.0 .AND. Z(IT).NE.0 )
     &              CALL CINTHFFFP(ADD_TO_MEZZL,IT,JKM1,JKM2,LM,
     &              ZGL(1,1,JKM1),ZGR(1,1,JKM2),JGL(1,1,JKM1),
     &              ZFL(1,1,JKM1),ZFR(1,1,JKM2),JFL(1,1,JKM1),TMPINT,
     &              RSCALHFF(1,ISF),MEZZ,MEZJ,MEZZL,MEZJL)
C-----------------------------------------------------------------------
C
            END DO
C--------------------------------------------------------------- ISF END
C
         END DO
      END DO
C                                                              JKM1 JKM2
C=======================================================================
C
C=======================================================================
C           Born solver: renormalize irregular matrix elements
C                        to account for spherical solutions
C=======================================================================
C------------------------------------------------------------ START BORN
      IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
         DO IME = 1,NMETOP
C
            MEZJ_SPHER = C0
            DO IR = 1,IRTOP
               MEZJ_SPHER = MEZJ_SPHER + RHOC(IR,IME)
            END DO
C
            CSUM = C0
            DO IKM = 1,NKM_T(IT)
               CSUM = CSUM + MEZJ(IKM,IKM,IT,IME)
            END DO
            CSCL = C1
            IF ( ABS(CSUM).GT.1D-12 ) CSCL = MEZJ_SPHER/CSUM
            ABS_DELTA = ABS(1D0-CSCL)
            IF ( IME.EQ.IDOS ) THEN
               TOL = 1D-6
            ELSE
               TOL = 5D-3
            END IF
            IF ( ABS_DELTA.GT.TOL ) THEN
               IF ( ABS_DELTA.GT.0.01D0 ) THEN
                  WRITE (6,99001) 'WARNING',IT,IME,ABS_DELTA
               ELSE
                  WRITE (6,99001) 'Info   ',IT,IME,ABS_DELTA
               END IF
            END IF
C
            N = NKM_T(IT)
            MEZJ(1:N,1:N,IT,IME) = CSCL*MEZJ(1:N,1:N,IT,IME)
C
            IF ( L_PROJECTED_ME ) THEN
               DO IL = 1,NL
                  MEZJL(1:N,1:N,IL,IME) = CSCL*MEZJL(1:N,1:N,IL,IME)
               END DO
            END IF
C
         END DO
      END IF
C-------------------------------------------------------------- END BORN
C=======================================================================
C
C???????????????????????????????????????????????????????????????????????
C                     CHECK MATRIX ELEMENTS
C???????????????????????????????????????????????????????????????????????
      IF ( .NOT.CHECK_OBS_ME ) RETURN
C
      ALLOCATE (MZBJA(NKMMAX,NKMMAX,3,NOBS))
      ALLOCATE (MZBZA(NKMMAX,NKMMAX,3,NOBS))
C
      CALL CALC_OBS_ME(IT,ZGL,ZFL,JGL,JFL,NCPLWF,IKMCPLWF,ZGR,ZFR,JGR,
     &                 JFR,NCPLWF,IKMCPLWF,MEZZ,MEZJ,MZBZA,MZBJA,NOBS,
     &                 'XXX',CHECK_OBS_ME)
C
C???????????????????????????????????????????????????????????????????????
C
      DEALLOCATE (RSCALHFF,RSCALHFFNUC,RSCAL,WINTR,CPLMAT,CPLMES)
C
99001 FORMAT (5X,A,' from <FPMATELM> for IT=',I3,' IME=',I2,
     &        '    | 1-CSCL | =',2F14.10)
      END
C*==cinthfffp.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE CINTHFFFP(ADD_TO_MEZZL,IT,JKM1,JKM2,LM,ZGL,ZGR,JGL,ZFL,
     &                     ZFR,JFL,RMEHF,RSCAL,MEZZ,MEZJ,MEZZL,MEZJL)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:B_AU2CGS
      USE MOD_ANGMOM,ONLY:NCPLWF,NMEMAX,NKMMAX,NLMAX,L_IKM,IHFF,AME_RLM
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IKMCPLWF,IMT
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL ADD_TO_MEZZL
      INTEGER IT,JKM1,JKM2,LM
      COMPLEX*16 JFL(NRMAX,NCPLWFMAX),JGL(NRMAX,NCPLWFMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           RMEHF(NCPLWFMAX,NCPLWFMAX),ZFL(NRMAX,NCPLWFMAX),
     &           ZFR(NRMAX,NCPLWFMAX),ZGL(NRMAX,NCPLWFMAX),
     &           ZGR(NRMAX,NCPLWFMAX)
      REAL*8 RSCAL(NRMAX)
C
C Local variables
C
      COMPLEX*16 AMEHFF,ZJHFF,ZZHFF
      INTEGER I,IKM1,IKM2,IL,IM,IRTOP,J,L1,L2
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IRTOP = JRCRI(IM)
C
      CALL CINTHFF(IM,IRTOP,ZGL,ZFL,ZGR,ZFR,RMEHF,NCPLWF(JKM1),
     &             NCPLWF(JKM2),RSCAL)
C
      DO I = 1,NCPLWF(JKM1)
         IKM1 = IKMCPLWF(I,JKM1)
         L1 = L_IKM(IKM1)
         IL = L1 + 1
         DO J = 1,NCPLWF(JKM2)
            IKM2 = IKMCPLWF(J,JKM2)
            L2 = L_IKM(IKM2)
C
            AMEHFF = AME_RLM(IKM1,IKM2,2,LM,IHFF)
C
            ZZHFF = B_AU2CGS*RMEHF(I,J)*AMEHFF
C
            MEZZ(JKM1,JKM2,IT,IHFF) = MEZZ(JKM1,JKM2,IT,IHFF) + ZZHFF
C
            IF ( ADD_TO_MEZZL .AND. L1.EQ.L2 ) MEZZL(JKM1,JKM2,IL,IHFF)
     &           = MEZZL(JKM1,JKM2,IL,IHFF) + ZZHFF
C
         END DO
      END DO
C
      IF ( JKM1.EQ.JKM2 ) THEN
C
         CALL CINTHFF(IM,IRTOP,ZGL,ZFL,JGL,JFL,RMEHF,NCPLWF(JKM1),
     &                NCPLWF(JKM2),RSCAL)
C
         DO I = 1,NCPLWF(JKM1)
            IKM1 = IKMCPLWF(I,JKM1)
            L1 = L_IKM(IKM1)
            IL = L1 + 1
            DO J = 1,NCPLWF(JKM2)
               IKM2 = IKMCPLWF(J,JKM2)
               L2 = L_IKM(IKM2)
C
               AMEHFF = AME_RLM(IKM1,IKM2,2,LM,IHFF)
C
               ZJHFF = B_AU2CGS*RMEHF(I,J)*AMEHFF
C
               MEZJ(JKM1,JKM2,IT,IHFF) = MEZJ(JKM1,JKM2,IT,IHFF) + ZJHFF
C
               IF ( ADD_TO_MEZZL .AND. L1.EQ.L2 )
     &              MEZJL(JKM1,JKM2,IL,IHFF) = MEZJL(JKM1,JKM2,IL,IHFF)
     &              + ZJHFF
            END DO
         END DO
      END IF
C
      END
C*==sumupintfp.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE SUMUPINTFP(CSUM,VG,G,WG,VF,F,WF,IKMCPL1,IKMCPL2,NCPL1,
     &                      NCPL2,NSF,NKMMAX,NKMAME,NSFMAX,NCPLWFMAX)
C   ********************************************************************
C   *                                                                  *
C   * routine to sum up spin-angular components of single site         *
C   * matrix elements                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:IMKM_IKM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CSUM,VF,VG
      INTEGER NCPL1,NCPL2,NCPLWFMAX,NKMAME,NKMMAX,NSF,NSFMAX
      COMPLEX*16 F(NKMMAX,NKMMAX,NSFMAX),G(NKMMAX,NKMMAX,NSFMAX),
     &           WF(NKMAME,NKMAME,NSFMAX),WG(NKMAME,NKMAME,NSFMAX)
      INTEGER IKMCPL1(NCPLWFMAX),IKMCPL2(NCPLWFMAX)
C
C Local variables
C
      INTEGER I,IKM1,IKM2,IMKM1,IMKM2,ISF,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,NCPL1
         IKM1 = IKMCPL1(J)
         IMKM1 = IMKM_IKM(IKM1)
C
         DO I = 1,NCPL2
            IKM2 = IKMCPL2(I)
            IMKM2 = IMKM_IKM(IKM2)
C
            DO ISF = 1,NSF
C
               CSUM = CSUM + VG*G(IKM1,IKM2,ISF)*WG(IKM1,IKM2,ISF)
     &                + VF*F(IKM1,IKM2,ISF)*WF(IMKM1,IMKM2,ISF)
            END DO
         END DO
      END DO
C
      END
C*==cintabrsum.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE CINTABRSUM(ADD_TO_MEZZL,IT,JKM1,JKM2,LM,ZGL,ZGR,JGL,
     &                      JGR,ZFL,ZFR,JFL,JFR,RHOC,WINTR,MEZZ,MEZJ,
     &                      MEZZL,MEZJL,NK1,NK2,IRBOT,IRTOP,IRNSBOT)
C   ********************************************************************
C   *                                                                  *
C   *  perform a radial integration with the r-dependent weights       *
C   *  set already                                                     *
C   *                                                                  *
C   *  NOTE: if  LHS_SOL_EQ_RHS_SOL = .FALSE.                          *
C   *        ZGL,ZFL,JGL,JFL correspond to the LHS solution            *
C   *        ZGR,ZFR,JGR,JFR correspond to the RHS solution            *
C   *                                                                  *
C   *        otherwise LHS = RHS solution with ZGL = ZGR etc           *
C   *                                                                  *
C   *  NOTE: symmetrization (ZJ+JZ)/2 not needed                       *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IKMCPLWF
      USE MOD_CALCMODE,ONLY:SOLVER_FP,BORN_USE_SPH
      USE MOD_ANGMOM,ONLY:NMEMAX,NLMAX,NKMMAX,IMKM_IKM,L_IKM,IDOS,ISMT,
     &    IOMT,AME_RLM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL SYMMETRIZE_ZJ_JZ
      PARAMETER (SYMMETRIZE_ZJ_JZ=.FALSE.)
C
C Dummy arguments
C
      LOGICAL ADD_TO_MEZZL
      INTEGER IRBOT,IRNSBOT,IRTOP,IT,JKM1,JKM2,LM,NK1,NK2
      COMPLEX*16 JFL(NRMAX,NCPLWFMAX),JFR(NRMAX,NCPLWFMAX),
     &           JGL(NRMAX,NCPLWFMAX),JGR(NRMAX,NCPLWFMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),RHOC(NRMAX,3),
     &           ZFL(NRMAX,NCPLWFMAX),ZFR(NRMAX,NCPLWFMAX),
     &           ZGL(NRMAX,NCPLWFMAX),ZGR(NRMAX,NCPLWFMAX)
      REAL*8 WINTR(NRMAX)
C
C Local variables
C
      COMPLEX*16 ALZFF,ALZGG,AR1FF,AR1GG,ASZFF,ASZGG,JFLW(:),JGLW(:),
     &           SUMZFJF,SUMZFZF,SUMZGJG,SUMZGZG,ZFLW(:),ZGLW(:),ZJDOS,
     &           ZJOMT,ZJSMT,ZZDOS,ZZOMT,ZZSMT
      LOGICAL CNON0
      INTEGER IKM1,IKM2,IL,IMKM1,IMKM2,IR,K1,K2,L1,L2,NR
      COMPLEX*16 ZDOTU
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JFLW,JGLW,ZFLW,ZGLW
C
      IF ( SYMMETRIZE_ZJ_JZ ) THEN
         ALLOCATE (JFLW(NRMAX),JGLW(NRMAX),ZFLW(NRMAX),ZGLW(NRMAX))
      ELSE
         ALLOCATE (ZFLW(NRMAX),ZGLW(NRMAX))
      END IF
C
C-----------------------------------------------------------------------
C
      NR = IRTOP - IRBOT + 1
C
      DO K1 = 1,NK1
         IKM1 = IKMCPLWF(K1,JKM1)
         IMKM1 = IMKM_IKM(IKM1)
         L1 = L_IKM(IKM1)
         IL = L1 + 1
C
         IF ( SYMMETRIZE_ZJ_JZ ) THEN
            DO IR = IRBOT,IRTOP
               ZGLW(IR) = ZGL(IR,K1)*WINTR(IR)
               ZFLW(IR) = ZFL(IR,K1)*WINTR(IR)
               JGLW(IR) = JGL(IR,K1)*WINTR(IR)
               JFLW(IR) = JFL(IR,K1)*WINTR(IR)
            END DO
         ELSE
            DO IR = IRBOT,IRTOP
               ZGLW(IR) = ZGL(IR,K1)*WINTR(IR)
               ZFLW(IR) = ZFL(IR,K1)*WINTR(IR)
            END DO
         END IF
C
         DO K2 = 1,NK2
            IKM2 = IKMCPLWF(K2,JKM2)
            IMKM2 = IMKM_IKM(IKM2)
            L2 = L_IKM(IKM2)
C
            AR1GG = AME_RLM(IKM1,IKM2,2,LM,IDOS)
            AR1FF = AME_RLM(IMKM1,IMKM2,2,LM,IDOS)
            ASZGG = AME_RLM(IKM1,IKM2,2,LM,ISMT)
            ASZFF = AME_RLM(IMKM1,IMKM2,2,LM,ISMT)
            ALZGG = AME_RLM(IKM1,IKM2,2,LM,IOMT)
            ALZFF = AME_RLM(IMKM1,IMKM2,2,LM,IOMT)
C
            IF ( CNON0(AR1GG) .OR. CNON0(ASZGG) .OR. CNON0(AR1FF) .OR. 
     &           CNON0(ASZFF) ) THEN
C
               SUMZGZG = ZDOTU(NR,ZGLW(IRBOT),1,ZGR(IRBOT,K2),1)
               SUMZFZF = ZDOTU(NR,ZFLW(IRBOT),1,ZFR(IRBOT,K2),1)
C
               ZZDOS = AR1GG*SUMZGZG + AR1FF*SUMZFZF
C
               ZZSMT = ASZGG*SUMZGZG - ASZFF*SUMZFZF
C
               ZZOMT = ALZGG*SUMZGZG - ALZFF*SUMZFZF
C
               MEZZ(JKM1,JKM2,IT,IDOS) = MEZZ(JKM1,JKM2,IT,IDOS) + ZZDOS
               MEZZ(JKM1,JKM2,IT,ISMT) = MEZZ(JKM1,JKM2,IT,ISMT) + ZZSMT
               MEZZ(JKM1,JKM2,IT,IOMT) = MEZZ(JKM1,JKM2,IT,IOMT) + ZZOMT
C
               IF ( ADD_TO_MEZZL .AND. L1.EQ.L2 ) THEN
                  MEZZL(JKM1,JKM2,IL,IDOS) = MEZZL(JKM1,JKM2,IL,IDOS)
     &               + ZZDOS
                  MEZZL(JKM1,JKM2,IL,ISMT) = MEZZL(JKM1,JKM2,IL,ISMT)
     &               + ZZSMT
                  MEZZL(JKM1,JKM2,IL,IOMT) = MEZZL(JKM1,JKM2,IL,IOMT)
     &               + ZZOMT
               END IF
C
               IF ( JKM1.EQ.JKM2 ) THEN
C
                  IF ( SYMMETRIZE_ZJ_JZ ) THEN
                     SUMZGJG = (ZDOTU(NR,ZGLW(IRBOT),1,JGR(IRBOT,K2),1)
     &                         +ZDOTU(NR,JGLW(IRBOT),1,ZGR(IRBOT,K2),1))
     &                         *0.5D0
                     SUMZFJF = (ZDOTU(NR,ZFLW(IRBOT),1,JFR(IRBOT,K2),1)
     &                         +ZDOTU(NR,JFLW(IRBOT),1,ZFR(IRBOT,K2),1))
     &                         *0.5D0
                  ELSE
                     SUMZGJG = ZDOTU(NR,ZGLW(IRBOT),1,JGR(IRBOT,K2),1)
                     SUMZFJF = ZDOTU(NR,ZFLW(IRBOT),1,JFR(IRBOT,K2),1)
                  END IF
C
                  ZJDOS = AR1GG*SUMZGJG + AR1FF*SUMZFJF
C
                  ZJSMT = ASZGG*SUMZGJG - ASZFF*SUMZFJF
C
                  ZJOMT = ALZGG*SUMZGJG - ALZFF*SUMZFJF
C
                  MEZJ(JKM1,JKM2,IT,IDOS) = MEZJ(JKM1,JKM2,IT,IDOS)
     &               + ZJDOS
                  MEZJ(JKM1,JKM2,IT,ISMT) = MEZJ(JKM1,JKM2,IT,ISMT)
     &               + ZJSMT
                  MEZJ(JKM1,JKM2,IT,IOMT) = MEZJ(JKM1,JKM2,IT,IOMT)
     &               + ZJOMT
C
                  IF ( ADD_TO_MEZZL .AND. L1.EQ.L2 ) THEN
                     MEZJL(JKM1,JKM2,IL,IDOS) = MEZJL(JKM1,JKM2,IL,IDOS)
     &                  + ZJDOS
                     MEZJL(JKM1,JKM2,IL,ISMT) = MEZJL(JKM1,JKM2,IL,ISMT)
     &                  + ZJSMT
                     MEZJL(JKM1,JKM2,IL,IOMT) = MEZJL(JKM1,JKM2,IL,IOMT)
     &                  + ZJOMT
                  END IF
C
C-----------------------------------------------------------------------
C                            Born solver
C-----------------------------------------------------------------------
C
C------------------------------------------------------------ START BORN
                  IF ( SOLVER_FP.EQ.'BORN' .AND. BORN_USE_SPH ) THEN
C
                     DO IR = MAX(IRNSBOT,IRBOT),IRTOP
C                    DO IR = IRBOT,IRTOP
                        IF ( SYMMETRIZE_ZJ_JZ ) THEN
                           SUMZGJG = (ZGLW(IR)*JGR(IR,K2)+JGLW(IR)
     &                               *ZGR(IR,K2))*0.5D0
                           SUMZFJF = (ZFLW(IR)*JFR(IR,K2)+JFLW(IR)
     &                               *ZFR(IR,K2))*0.5D0
                        ELSE
                           SUMZGJG = ZGLW(IR)*JGR(IR,K2)
                           SUMZFJF = ZFLW(IR)*JFR(IR,K2)
                        END IF
C
                        RHOC(IR,1) = RHOC(IR,1) + AR1GG*SUMZGJG + 
     &                               AR1FF*SUMZFJF
C
                        RHOC(IR,2) = RHOC(IR,2) + ASZGG*SUMZGJG - 
     &                               ASZFF*SUMZFJF
C
                        RHOC(IR,3) = RHOC(IR,3) + ALZGG*SUMZGJG - 
     &                               ALZFF*SUMZFJF
                     END DO
C-----------------------------------------------------------------------
C
                  END IF
C-------------------------------------------------------------- END BORN
C
               END IF
C
            END IF
C
         END DO
C
      END DO
C
      END
C*==cintabrsum_spher.f    processed by SPAG 6.70Rc at 11:01 on 29 Jan 2017
      SUBROUTINE CINTABRSUM_SPHER(JKM1,JKM2,LM,CHI,NCWF_SOL_SPH,
     &                            IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,JG_SPH,
     &                            JF_SPH,ZFJF_SPH,ZFZF_SPH,ZGJG_SPH,
     &                            ZGZG_SPH,RHOC,WINTR,IRNSBOT,NMETOP)
C   ********************************************************************
C   *                                                                  *
C   *    deal with the contribution of the spherical solutions         *
C   *    for r < r_ns  and  SOLVER_FP = BORN                           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKMMAX,IMKM_IKM,IDOS,ISMT,IOMT,AME_RLM
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRNSBOT,JKM1,JKM2,LM,NMETOP
      COMPLEX*16 CHI(NKMMAX,NKMMAX),JF_SPH(NRMAX,2,NKMMAX),
     &           JG_SPH(NRMAX,2,NKMMAX),RHOC(NRMAX,3),ZFJF_SPH(NRMAX),
     &           ZFZF_SPH(NRMAX),ZF_SPH(NRMAX,2,NKMMAX),ZGJG_SPH(NRMAX),
     &           ZGZG_SPH(NRMAX),ZG_SPH(NRMAX,2,NKMMAX)
      INTEGER IKM_CWFSOL_SPH(2,NKMMAX),NCWF_SOL_SPH(NKMMAX)
      REAL*8 WINTR(NRMAX)
C
C Local variables
C
      COMPLEX*16 CHI12,WZFJF_SPH(3),WZFZF_SPH(3),WZGJG_SPH(3),
     &           WZGZG_SPH(3)
      LOGICAL CNON0
      INTEGER ICWF1,ICWF2,IKM3,IKM4,IME,IMKM3,IMKM4,IR
      LOGICAL JKM1_EQ_JKM2
C
C*** End of declarations rewritten by SPAG
C
      JKM1_EQ_JKM2 = JKM1.EQ.JKM2
      CHI12 = CHI(JKM1,JKM2)
C
C=======================================================================
      DO ICWF1 = 1,NCWF_SOL_SPH(JKM1)
         IKM3 = IKM_CWFSOL_SPH(ICWF1,JKM1)
         IMKM3 = IMKM_IKM(IKM3)
C
C ----------------------------------------------------------------------
         DO ICWF2 = 1,NCWF_SOL_SPH(JKM2)
            IKM4 = IKM_CWFSOL_SPH(ICWF2,JKM2)
            IMKM4 = IMKM_IKM(IKM4)
C
C ----------------------------------------------------------------------
C
            WZGJG_SPH(IDOS) = AME_RLM(IKM3,IKM4,2,LM,IDOS)
            WZFJF_SPH(IDOS) = AME_RLM(IMKM3,IMKM4,2,LM,IDOS)
            WZGJG_SPH(ISMT) = AME_RLM(IKM3,IKM4,2,LM,ISMT)
            WZFJF_SPH(ISMT) = AME_RLM(IMKM3,IMKM4,2,LM,ISMT)
            WZGJG_SPH(IOMT) = AME_RLM(IKM3,IKM4,2,LM,IOMT)
            WZFJF_SPH(IOMT) = AME_RLM(IMKM3,IMKM4,2,LM,IOMT)
C
            DO IME = 1,NMETOP
               WZGZG_SPH(IME) = WZGJG_SPH(IME)*CHI12
               WZFZF_SPH(IME) = WZFJF_SPH(IME)*CHI12
            END DO
C
            IF ( CNON0(WZGJG_SPH(IDOS)) .OR. CNON0(WZGJG_SPH(ISMT)) .OR. 
     &           CNON0(WZFJF_SPH(IDOS)) .OR. CNON0(WZFJF_SPH(ISMT)) )
     &           THEN
C
C -------------------- set up radial functions including  irregular part
C
               DO IR = 1,IRNSBOT - 1
                  ZGZG_SPH(IR) = ZG_SPH(IR,ICWF1,JKM1)
     &                           *ZG_SPH(IR,ICWF2,JKM2)*WINTR(IR)
                  ZFZF_SPH(IR) = ZF_SPH(IR,ICWF1,JKM1)
     &                           *ZF_SPH(IR,ICWF2,JKM2)*WINTR(IR)
                  ZGJG_SPH(IR) = ZG_SPH(IR,ICWF1,JKM1)
     &                           *JG_SPH(IR,ICWF2,JKM2)*WINTR(IR)
                  ZFJF_SPH(IR) = ZF_SPH(IR,ICWF1,JKM1)
     &                           *JF_SPH(IR,ICWF2,JKM2)*WINTR(IR)
               END DO
C
               IF ( JKM1_EQ_JKM2 ) THEN
C
                  DO IME = 1,NMETOP
                     DO IR = 1,IRNSBOT - 1
                        RHOC(IR,IME) = RHOC(IR,IME) + WZGZG_SPH(IME)
     &                                 *ZGZG_SPH(IR) + WZFZF_SPH(IME)
     &                                 *ZFZF_SPH(IR) + WZGJG_SPH(IME)
     &                                 *ZGJG_SPH(IR) + WZFJF_SPH(IME)
     &                                 *ZFJF_SPH(IR)
                     END DO
                  END DO
C
               ELSE
C
                  DO IME = 1,NMETOP
                     DO IR = 1,IRNSBOT - 1
                        RHOC(IR,IME) = RHOC(IR,IME) + WZGZG_SPH(IME)
     &                                 *ZGZG_SPH(IR) + WZFZF_SPH(IME)
     &                                 *ZFZF_SPH(IR)
                     END DO
                  END DO
C
               END IF
C
            END IF
C
         END DO
      END DO
C ======================================================================
C
      END
