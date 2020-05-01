C*==ssite.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SSITE(IWRREGWF,IWRIRRWF,IFILSS,CALCINT,GETIRRSOL,ERYD,
     &                 P,IPRINT,NKMSS,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C   ********************************************************************
C   *                                                                  *
C   * ASSIGN QUANTUM NUMBERS AND CALL ROUTINE TO SOLVE                 *
C   * 8 COUPLED PARTIAL DIFFERENTIAL RADIAL DIRAC EQUATIONS. THE       *
C   * RESULTING WAVEFUNCTIONS ARE USED TO CALCULATE T-MATRICES IN      *
C   * THE KAPPA-MU REPRESENTATION                                      *
C   *                                                                  *
C   * + CALCULATION OF THE RADIAL INTEGRALS                            *
C   *   [ G1*G2 + F1*F2 ] R**2 DR                                      *
C   *                                                                  *
C   * FOR IHYPER <> 0 :                                                *
C   * CALCULATION OF THE HYPERFINE MATRIX ELEMENTS                     *
C   *                                                                  *
C   * RYD-UNITS USED THROUGHOUT                                        *
C   *                                                                  *
C   * NOTE: to save storage force  JG/JF  and  PR/QR  to share the     *
C   *       same storage by corresponding argument list in CALL ....   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FILNAM,LFILNAM,IOTMP
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R2DRDI_W_RADINT,DRDI,DRDIOVR,R,
     &    FINITE_NUCLEUS
      USE MOD_ANGMOM,ONLY:IKM1LIN,NMEMAX,NKMMAX,LINMAX,NLMAX,AME_G,
     &    IKM2LIN,AMEBI1,AMEBI2,IDOS,ISMT,IOMT,IHFF,ISDM,IKDS
      USE MOD_CALCMODE,ONLY:IHYPER,SOLVER,BREITINT,GF_CONV_RH,
     &    BREAKPOINT,KMROT
      USE MOD_TYPES,ONLY:SOCTL,NTMAX,LOPT,IMT,Z,BT,VT,CTL,ABIT,AOPT,NLT,
     &    NLIN_T,NKM_T,ITBOT,ITTOP,NT,NCPLWFMAX
      USE MOD_CONSTANTS,ONLY:B_AU2CGS,CI,C0
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,SEVT,SEBT
      IMPLICIT NONE
C*--SSITE36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SSITE')
      LOGICAL WRONSKI
      PARAMETER (WRONSKI=.FALSE.)
      REAL*8 F1
      PARAMETER (F1=1.0D0)
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      COMPLEX*16 ERYD,P
      INTEGER IFILSS,IPRINT,IWRIRRWF,IWRREGWF,NKMSS
      CHARACTER*10 ORBPOLSS
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 A(2,2),ARG,B1,B2,BZJ(:),BZZ(:),CB(:),CHL,CHLB1,CHLB2,
     &           CJL,CJL1,CJLB1,CJLB2,CJLP,CNL,CNLB1,CNLB2,CNLP,CRSQ,
     &           CSUM,CV(:),DET,DXP(2,2),DZJ(:),DZZ(:),F11,F12,F21,F22,
     &           G11,G11P,G12,G12P,G21,G21P,G22,G22P,GAM(2,2),
     &           GAMINV(2,2),JF(:,:,:),JG(:,:,:),KZJ(:),KZZ(:),MAUX(:,:)
     &           ,MSST2(2,2),NORM,OZJ(:),OZZ(:),PFAC,PI(:,:,:),PR(:,:,:)
     &           ,QI(:,:,:),QR(:,:,:),RMEHF(NCPLWFMAX,NCPLWFMAX),
     &           RMEHF1(2,2),RMEHF2(2,2),SIG(2,2),SSSLIN(:),SZJ(:),
     &           SZZ(:),TSSLIN(:),TSST2(2,2),TZJ(:),TZZ(:),W(2,2),
     &           WRON(2,2),XSST2(2,2),ZF(:,:,:),ZFJF(2,2),ZFZF(2,2),
     &           ZG(:,:,:),ZGJG(2,2),ZGZG(2,2)
      REAL*8 AP(:,:,:),AQ(:,:,:),C,CFF(2,2),CFG(2,2),CG1,CG2,CG4,CG5,
     &       CG8,CGF(2,2),CGG(2,2),CHF(2,2),CSQR,CTF(2,2),CTG(2,2),
     &       DROVRN(:),MJ,R1M(2,2),RKD(2,2),RNUC,SCL,SK1,SK2,TDIA1,
     &       TDIA2,TOFF
      COMPLEX*16 CDJLZDZ,CDNLZDZ,CJLZ,CNLZ
      LOGICAL DUMP_SSITE,MODIFYV
      INTEGER I,I1,I2,I3,IKM(2),IKM1,IKM2,IKMBOT,IKMI,IKMJ,IKMTOP,IL,
     &        ILA,IM,IMKM(2),IMKMI,INFO,IPIV(:),IR,IRTOP,ISK1,ISK2,IT,J,
     &        JLIM,JRCUT(0:1),K,K1,K2,KAP1,KAP2,KC,L,L1,LB1,LB2,LIN,M,
     &        MUM05,N,NKM,NPAN,NSOL,NSTEP,NVIEW
      INTEGER IKAPMUE
      REAL*8 RNUCTAB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,IPIV,CV,CB,JF,JG,PI,PR,QI,QR,ZF,ZG
      ALLOCATABLE AP,AQ,DROVRN
      ALLOCATABLE BZJ,BZZ,DZJ,DZZ,KZJ,KZZ,OZJ,OZZ
      ALLOCATABLE SSSLIN,SZJ,SZZ,TSSLIN,TZJ,TZZ
C
      DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
      DATA RKD/1.0D0,0.0D0,0.0D0, - 1.0D0/
      DATA CG1/0D0/,CG2/0D0/,CG4/0D0/,CG5/0D0/,CG8/0D0/
      DATA DUMP_SSITE/.FALSE./
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (MAUX(NKMMAX,NKMMAX),IPIV(NKMMAX))
      ALLOCATE (CB(NRMAX),CV(NRMAX))
      ALLOCATE (JF(NRMAX,2,2),JG(NRMAX,2,2),PI(2,2,NRMAX),PR(2,2,NRMAX))
      ALLOCATE (QI(2,2,NRMAX),QR(2,2,NRMAX),ZF(NRMAX,2,2),ZG(NRMAX,2,2))
      ALLOCATE (AP(2,2,NRMAX),AQ(2,2,NRMAX),DROVRN(NRMAX))
      ALLOCATE (BZJ(LINMAX),BZZ(LINMAX),DZJ(LINMAX),DZZ(LINMAX))
      ALLOCATE (KZJ(LINMAX),KZZ(LINMAX),OZJ(LINMAX),OZZ(LINMAX))
      ALLOCATE (SSSLIN(LINMAX),SZJ(LINMAX),SZZ(LINMAX),TSSLIN(LINMAX))
      ALLOCATE (TZJ(LINMAX),TZZ(LINMAX))
C
      NKM = NKMSS
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) THEN
         IPRINT = 5
         DUMP_SSITE = .TRUE.
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C-----------------------------------------------------------------------
      IF ( GF_CONV_RH .AND. 
     &     .NOT.(SOLVER.EQ.'BS       ' .OR. SOLVER.EQ.'BS-SOC    ') )
     &     THEN
         WRITE (6,99008) SOLVER
         CALL STOP_MESSAGE(ROUTINE,'GF_CONV_RH=.T. is not consistent')
      END IF
C-----------------------------------------------------------------------
C
      C = CTL(1,1)
      CSQR = C*C
C
C     calculate relativistic momentum
C
      CALL GET_MOMENTUM(3,C,ERYD,P)
C
      PFAC = P/(1.0D0+ERYD/CSQR)
C
      IF ( IPRINT.GT.0 ) WRITE (6,99001) ERYD
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         MEZZ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
         MEZJ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
C
         DO J = 1,NKMMAX
            DO I = 1,NKMMAX
               TSST(I,J,IT) = C0
               MSST(I,J,IT) = C0
               SSST(I,J,IT) = C0
            END DO
C----------- fill diagonal with 1 to allow matrix inversion for NLT < NL
C                 NOTE: this dummy setting will be overwritten up to NLT
            TSST(J,J,IT) = 1D0
            MSST(J,J,IT) = 1D0
            SSST(J,J,IT) = 1D0
         END DO
C
         DO I = 1,LINMAX
            DZZ(I) = C0
            DZJ(I) = C0
            SZZ(I) = C0
            SZJ(I) = C0
            OZZ(I) = C0
            OZJ(I) = C0
            BZZ(I) = C0
            BZJ(I) = C0
            TZZ(I) = C0
            TZJ(I) = C0
            KZZ(I) = C0
            KZJ(I) = C0
         END DO
C
         IF ( IPRINT.GT.0 ) WRITE (6,99002) 'SOLVER',IT,SOLVER,
     &                             (SOCTL(IT,IL),IL=1,NLMAX)
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         NPAN = 1
         JRCUT(0) = 0
         JRCUT(1) = IRTOP
C
         IF ( FINITE_NUCLEUS ) THEN
            RNUC = RNUCTAB(Z(IT))
            IR = 1
            DO WHILE ( R(IR,IM).LT.RNUC )
               IR = IR + 1
            END DO
C
            JLIM = IR
            IF ( MOD(JLIM,2).EQ.0 ) JLIM = JLIM - 1
            RNUC = R(JLIM,IM)
C
            DO IR = 1,IRTOP
               DROVRN(IR) = (R(IR,IM)/RNUC)**3*DRDI(IR,IM)
            END DO
         END IF
C
         DO IR = 1,IRTOP
            CV(IR) = VT(IR,IT)
            CB(IR) = BT(IR,IT)
         END DO
C
         LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,(NLT(IT)-1)
            IL = L + 1
            C = CTL(IT,IL)
            CSQR = C*C
C
            MODIFYV = .FALSE.
            IF ( L.EQ.LOPT(IT) ) THEN
               IF ( ORBPOLSS(1:5).EQ.'SIGMA' .OR. ORBPOLSS(1:4)
     &              .EQ.'DMFT' .OR. ORBPOLSS(1:5).EQ.'LDA+U' )
     &              MODIFYV = .TRUE.
            END IF
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
            ISK1 = ISIGN(1,KAP1)
            ISK2 = ISIGN(1,KAP2)
            SK1 = DBLE(ISK1)
            SK2 = DBLE(ISK2)
            L1 = L
            LB1 = L - ISK1
            LB2 = L - ISK2
C
            ARG = P*R(IRTOP,IM)
            CJL = CJLZ(L1,ARG)
            CJLB1 = CJLZ(LB1,ARG)
            CJLB2 = CJLZ(LB2,ARG)
            CNL = CNLZ(L1,ARG)
            CNLB1 = CNLZ(LB1,ARG)
            CNLB2 = CNLZ(LB2,ARG)
            CHL = CJL + CI*CNL
            CHLB1 = CJLB1 + CI*CNLB1
            CHLB2 = CJLB2 + CI*CNLB2
C
            CJL1 = CJLZ(L,P*R(1,IM))
C
            IF ( SOLVER(1:6).EQ.'BS-SOC' ) THEN
               CJLP = CDJLZDZ(L1,ARG,1)*P
               CNLP = CDNLZDZ(L1,ARG,1)*P
            END IF
            IF ( MODIFYV ) THEN
               DO I = 1,IRTOP
                  CV(I) = VT(I,IT) + SEVT(I,IT)
                  CB(I) = BT(I,IT) + SEBT(I,IT)
               END DO
            END IF
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MUM05 = -L - 1, + L
               MJ = DBLE(MUM05) + 0.5D0
C
C-----------------------------------------------------------------------
C NO COUPLING FOR:  ABS(MUE)= J   +  J=L+1/2 == KAP=-L-1
               IF ( ABS(MJ).GE.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C-----------------------------------------------------------------------
               IKM1 = IKAPMUE(KAP1,MUM05)
               IKM2 = IKAPMUE(KAP2,MUM05)
               IKM(1) = IKM1
               IKM(2) = IKM2
               IMKM(1) = IKAPMUE(-KAP1,MUM05)
               IMKM(2) = IKAPMUE(-KAP2,NINT(MJ-0.5D0))
C
C-----------------------------------------------------------------------
C
               CALL RINIT(4,CGF)
               CALL RINIT(4,CFF)
               CALL RINIT(4,CTF)
C
               DO I = 1,NSOL
                  IKMI = IKM(I)
                  DO J = 1,NSOL
                     IKMJ = IKM(J)
                     CGG(I,J) = AME_G(IKMI,IKMJ,2,ISMT)
                     CFG(I,J) = AME_G(IKMI,IKMJ,2,IOMT)
                     CHF(I,J) = AME_G(IKMI,IKMJ,2,IHFF)
                     CTG(I,J) = AME_G(IKMI,IKMJ,2,ISDM)
                  END DO
C
                  IMKMI = IMKM(I)
                  CGF(I,I) = AME_G(IMKMI,IMKMI,2,ISMT)
                  CFF(I,I) = AME_G(IMKMI,IMKMI,2,IOMT)
                  CTF(I,I) = AME_G(IMKMI,IMKMI,2,ISDM)
               END DO
C
C
C-----------------------------------------------------------------------
C
               CALL RINIT(2*2*NRMAX,AP)
               CALL RINIT(2*2*NRMAX,AQ)
C
               IF ( BREITINT ) THEN
                  DO J = 1,NSOL
                     IKMJ = IKM(J)
                     DO I = 1,NSOL
                        IKMI = IKM(I)
                        DO ILA = 1,(2*L+1+1)/2
                           DO M = -1, + 1,2
                              DO N = 1,IRTOP
                                 AP(I,J,N) = AP(I,J,N)
     &                              + ABIT(N,ILA,M,IT)
     &                              *AMEBI2(IKMI,IKMJ,ILA,M)
                                 AQ(I,J,N) = AQ(I,J,N)
     &                              + ABIT(N,ILA,M,IT)
     &                              *AMEBI1(IKMI,IKMJ,ILA,M)
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
C------------------------------------------------- scale by: e/c in a.u.
                  SCL = SQRT(2.0D0)/C
                  AP(1:2,1:2,1:NRMAX) = SCL*AP(1:2,1:2,1:NRMAX)
                  AQ(1:2,1:2,1:NRMAX) = SCL*AQ(1:2,1:2,1:NRMAX)
               END IF
C
C-----------------------------------------------------------------------
C
               W(1:2,1:2) = C0
               IF ( L.EQ.LOPT(IT) .AND. 
     &              (ORBPOLSS(1:4).EQ.'DMFT' .OR. ORBPOLSS(1:5)
     &              .EQ.'LDA+U') ) THEN
                  IKMBOT = 2*LOPT(IT)**2 + 1
                  IKMTOP = IKMBOT - 1 + 2*(2*LOPT(IT)+1)
                  DO I = 1,NSOL
                     IKMI = IKM(I)
                     IF ( IKMI.GE.IKMBOT .AND. IKMI.LE.IKMTOP ) THEN
                        DO J = 1,NSOL
                           IKMJ = IKM(J)
                           IF ( IKMJ.GE.IKMBOT .AND. IKMJ.LE.IKMTOP )
     &                          W(I,J) = DMFTSIG(IKMI,IKMJ,IT)
                        END DO
                     END IF
                  END DO
               END IF
C
C-----------------------------------------------------------------------
C
               ZG(:,:,:) = C0
               ZF(:,:,:) = C0
               JG(:,:,:) = C0
               JF(:,:,:) = C0
               PR(:,:,:) = C0
               QR(:,:,:) = C0
               PI(:,:,:) = C0
               QI(:,:,:) = C0
C
               IF ( SOLVER.EQ.'BS' .AND. 
     &              .NOT.(ORBPOLSS(1:4).EQ.'DMFT' .OR. ORBPOLSS(1:5)
     &              .EQ.'LDA+U') ) THEN
C
                  CALL DIRBS(GETIRRSOL,CTL(IT,IL),ERYD,L,MJ,KAP1,KAP2,P,
     &                       CG1,CG2,CG4,CG5,CG8,CV,CB,Z(IT),R(1,IM),
     &                       DRDI(1,IM),DRDIOVR(1,IM),IRTOP,JRCUT,NPAN,
     &                       PR,QR,PI,QI,ZG,ZF,NRMAX,GF_CONV_RH)
C
               ELSE IF ( SOLVER.EQ.'BS' .AND. 
     &                   (ORBPOLSS(1:4).EQ.'DMFT' .OR. ORBPOLSS(1:5)
     &                   .EQ.'LDA+U') ) THEN
C
                  CALL DIRBSSIG(GETIRRSOL,CTL(IT,IL),ERYD,L,MJ,KAP1,
     &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,CV,CB,Z(IT),
     &                          R(1,IM),DRDI(1,IM),DRDIOVR(1,IM),IRTOP,
     &                          JRCUT,NPAN,PR,QR,PI,QI,ZG,ZF,NRMAX,
     &                          GF_CONV_RH,W)
C
               ELSE IF ( SOLVER.EQ.'RK        ' ) THEN
C
                  CALL DIRRK(GETIRRSOL,CTL(IT,IL),ERYD,L,MJ,KAP1,KAP2,P,
     &                       CG1,CG2,CG4,CG5,CG8,CV,CB,Z(IT),R(1,IM),
     &                       DRDI(1,IM),DRDIOVR(1,IM),IRTOP,JRCUT,NPAN,
     &                       PR,QR,PI,QI,ZG,ZF,NRMAX)
C
               ELSE IF ( SOLVER.EQ.'BS-BI     ' ) THEN
C
                  CALL DIRBSBI(GETIRRSOL,CTL(IT,IL),ERYD,L,MJ,KAP1,KAP2,
     &                         P,CG1,CG2,CG4,CG5,CG8,CV,CB,Z(IT),R(1,IM)
     &                         ,DRDI(1,IM),DRDIOVR(1,IM),IRTOP,JRCUT,
     &                         NPAN,PR,QR,PI,QI,ZG,ZF,AP,AQ,NRMAX)
C
               ELSE IF ( SOLVER.EQ.'ABM-OP    ' ) THEN
C
                  CALL DIRABMOP(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
     &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,VT(1,IT),
     &                          BT(1,IT),AOPT,Z(IT),R(1,IM),DRDI(1,IM),
     &                          DRDIOVR(1,IM),IRTOP,PR,QR,PI,QI,ZG,ZF,
     &                          AP,AQ,LOPT(IT),NTMAX,NRMAX)
C
               ELSE IF ( SOLVER(1:6).EQ.'BS-SOC' ) THEN
C
                  CALL DIRBS_SOC(GETIRRSOL,CTL(IT,IL),SOCTL(IT,IL),ERYD,
     &                           L,MJ,KAP1,KAP2,P,CG1,CG2,CG4,CG5,CG8,
     &                           CV,CB,Z(IT),R(1,IM),DRDI(1,IM),
     &                           DRDIOVR(1,IM),IRTOP,JRCUT,NPAN,DXP,PR,
     &                           QR,PI,QI,ZG,ZF,NRMAX,GF_CONV_RH)
C
               ELSE
                  WRITE (6,99008) SOLVER
                  CALL STOP_MESSAGE(ROUTINE,'solver not found')
               END IF
C
C  wavefunctions at the muffin-tin-radius
C
               G11 = PR(1,1,IRTOP)/R(IRTOP,IM)
               G12 = PR(1,2,IRTOP)/R(IRTOP,IM)
               G21 = PR(2,1,IRTOP)/R(IRTOP,IM)
               G22 = PR(2,2,IRTOP)/R(IRTOP,IM)
               F11 = QR(1,1,IRTOP)/(R(IRTOP,IM)*C)
               F12 = QR(1,2,IRTOP)/(R(IRTOP,IM)*C)
               F21 = QR(2,1,IRTOP)/(R(IRTOP,IM)*C)
               F22 = QR(2,2,IRTOP)/(R(IRTOP,IM)*C)
C
               IF ( SOLVER(1:6).EQ.'BS-SOC' ) THEN
                  G11P = (DXP(1,1)/DRDI(IRTOP,IM)-G11)/R(IRTOP,IM)
                  G21P = (DXP(2,1)/DRDI(IRTOP,IM)-G21)/R(IRTOP,IM)
                  G12P = (DXP(1,2)/DRDI(IRTOP,IM)-G12)/R(IRTOP,IM)
                  G22P = (DXP(2,2)/DRDI(IRTOP,IM)-G22)/R(IRTOP,IM)
               END IF
C
C ------- the minor component for the soc-manipulated wf is meaningless
C
               IF ( SOLVER(1:6).EQ.'BS-SOC' ) THEN
                  CALL CINIT(2*2*NRMAX,QR)
                  CALL CINIT(2*2*NRMAX,QI)
               END IF
C
C      COSD= NL * C * F11 - PFAC * SK1 * NLB1 * G11
C      SIND= JL * C * F11 - PFAC * SK1 * JLB1 * G11
C
C -------------------------------------------------------------------
C       T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
C -------------------------------------------------------------------
C
               CNL = (CHL-CJL)/CI
               CNLB1 = (CHLB1-CJLB1)/CI
               CNLB2 = (CHLB2-CJLB2)/CI
C
               IF ( SOLVER(1:6).EQ.'BS-SOC' ) THEN
                  GAM(1,1) = CJL*G11P - CJLP*G11
                  GAM(2,1) = CJL*G21P - CJLP*G21
                  GAM(1,2) = CJL*G12P - CJLP*G12
                  GAM(2,2) = CJL*G22P - CJLP*G22
C
                  SIG(1,1) = CNL*G11P - CNLP*G11
                  SIG(2,1) = CNL*G21P - CNLP*G21
                  SIG(1,2) = CNL*G12P - CNLP*G12
                  SIG(2,2) = CNL*G22P - CNLP*G22
               ELSE
                  GAM(1,1) = CJL*C*F11 - PFAC*SK1*CJLB1*G11
                  GAM(2,1) = CJL*C*F21 - PFAC*SK2*CJLB2*G21
                  GAM(1,2) = CJL*C*F12 - PFAC*SK1*CJLB1*G12
                  GAM(2,2) = CJL*C*F22 - PFAC*SK2*CJLB2*G22
C
                  SIG(1,1) = CNL*C*F11 - PFAC*SK1*CNLB1*G11
                  SIG(2,1) = CNL*C*F21 - PFAC*SK2*CNLB2*G21
                  SIG(1,2) = CNL*C*F12 - PFAC*SK1*CNLB1*G12
                  SIG(2,2) = CNL*C*F22 - PFAC*SK2*CNLB2*G22
               END IF
C
               GAMINV(1:NSOL,1:NSOL) = GAM(1:NSOL,1:NSOL)
               CALL ZGETRF(NSOL,NSOL,GAMINV,2,IPIV,INFO)
               CALL ZGETRI(NSOL,GAMINV,2,IPIV,MAUX,2*2,INFO)
C
               DO I2 = 1,NSOL
                  DO I1 = 1,NSOL
                     CSUM = 0.0D0
                     DO I3 = 1,NSOL
                        CSUM = CSUM + SIG(I1,I3)*GAMINV(I3,I2)
                     END DO
                     XSST2(I1,I2) = P*CSUM
                  END DO
               END DO
C
               DO I1 = 1,NSOL
                  DO I2 = 1,NSOL
                     MSST2(I1,I2) = -XSST2(I1,I2)
                  END DO
                  MSST2(I1,I1) = MSST2(I1,I1) + CI*P
               END DO
C
               TSST2(1:NSOL,1:NSOL) = MSST2(1:NSOL,1:NSOL)
               CALL ZGETRF(NSOL,NSOL,TSST2,2,IPIV,INFO)
               CALL ZGETRI(NSOL,TSST2,2,IPIV,MAUX,2*2,INFO)
C
C-----------------------------------------------------------------------
C
C   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
C
               CGG(1,1) = CG1
               CGG(1,2) = CG2
               CGG(2,1) = CG2
               CGG(2,2) = CG4
               CALL RINIT(4,CGF)
               CGF(1,1) = CG5
               CGF(2,2) = CG8
C
C   COEFFICIENTS TO CALCULATE THE SPIN  DIPOLAR MOMENT TZ
C
               TDIA1 = 2*MJ/DBLE((2*L1+1)*(2*LB1+1))
               TDIA2 = 2*MJ/DBLE((2*L1+1)*(2*LB2+1))
               TOFF = -SQRT((L1+0.5D0)**2-MJ**2)/DBLE(2*L1+1)
C
               CTG(1,1) = 0.5D0*(CG1-3.0D0*TDIA1)
               CTG(1,2) = 0.5D0*(CG2-3.0D0*TOFF)
               CTG(2,1) = 0.5D0*(CG2-3.0D0*TOFF)
               CTG(2,2) = 0.5D0*(CG4-3.0D0*TDIA2)
               CALL RINIT(4,CTF)
C
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C-----------------------------------------------------------------------
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
               CHF(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
               CHF(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
               IF ( NSOL.EQ.2 ) THEN
                  CHF(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
                  CHF(2,1) = CHF(1,2)
               END IF
C-----------------------------------------------------------------------
CALCULATE RADIAL INTEGRALS  UP TO   OR RWS(IRTOP=JRWS)
C
C
               IF ( NSOL.EQ.1 ) THEN
C====================================================================
C NO COUPLING TO OTHER SCATTERING CHANNELS
C REGULAR PART    Z*Z    Z = (GRA,FRA)
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
                  NORM = (P*CNL-CJL*XSST2(1,1))/G11
C                                                basis functions R and H
                  IF ( GF_CONV_RH ) NORM = NORM*TSST2(1,1)
C
                  DO I = 1,IRTOP
                     ZG(I,1,1) = (PR(1,1,I)/R(I,IM))*NORM
                     ZF(I,1,1) = (QR(1,1,I)/R(I,IM)/C)*NORM
                     JG(I,1,1) = PI(1,1,I)/R(I,IM)
                     JF(I,1,1) = QI(1,1,I)/R(I,IM)/C
                  END DO
C
                  IF ( IWRREGWF.NE.0 )
     &                 WRITE (IFILSS,REC=IKM1+(IT-1)*NKM) IT,'REG',IKM1,
     &                 IRTOP,NSOL,IKM1,((ZG(I,K,1),I=1,IRTOP),K=1,NSOL),
     &                 ((ZF(I,K,1),I=1,IRTOP),K=1,NSOL)
C
                  IF ( IWRIRRWF.NE.0 )
     &                 WRITE (IFILSS,REC=IKM1+(IT-1+NT)*NKM) IT,'IRR',
     &                 IKM1,IRTOP,((JG(I,K,1),I=1,IRTOP),K=1,NSOL),
     &                 ((JF(I,K,1),I=1,IRTOP),K=1,NSOL)
C
C============================================== NO COUPLING = END ===
               ELSE
C====================================================================
C COUPLING OF TWO SCATTERING CHANNELS
C   Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
C              INDEX 2: BOUNDARY CONDITION
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
                  DET = G11*G22 - G12*G21
C
COEFFICIENTS TO GET:   Z(K1,K1)  Z(K2,K1)
                  B1 = P*CNL - XSST2(1,1)*CJL
                  B2 = -XSST2(2,1)*CJL
                  A(1,1) = (G22*B1-G12*B2)/DET
                  A(2,1) = (G11*B2-G21*B1)/DET
C
COEFFICIENTS TO GET:   Z(K1,K2)  Z(K2,K2)
                  B1 = -XSST2(1,2)*CJL
                  B2 = P*CNL - XSST2(2,2)*CJL
                  A(1,2) = (G22*B1-G12*B2)/DET
                  A(2,2) = (G11*B2-G21*B1)/DET
C
C                                                basis functions R and H
                  IF ( GF_CONV_RH ) THEN
                     GAM(1:2,1:2) = MATMUL(A(1:2,1:2),TSST2(1:2,1:2))
                     A(1:2,1:2) = GAM(1:2,1:2)
                  END IF
C
CALCULATE FUNCTIONS: Z(K1,K1), Z(K2,K1), Z(K1,K2), Z(K2,K2)
                  DO K = 1,NSOL
                     DO I = 1,IRTOP
                        ZG(I,1,K) = (PR(1,1,I)*A(1,K)+PR(1,2,I)*A(2,K))
     &                              /R(I,IM)
                        ZF(I,1,K) = (QR(1,1,I)*A(1,K)+QR(1,2,I)*A(2,K))
     &                              /R(I,IM)/C
C
                        ZG(I,2,K) = (PR(2,1,I)*A(1,K)+PR(2,2,I)*A(2,K))
     &                              /R(I,IM)
                        ZF(I,2,K) = (QR(2,1,I)*A(1,K)+QR(2,2,I)*A(2,K))
     &                              /R(I,IM)/C
                     END DO
                  END DO
                  DO K = 1,NSOL
                     KC = 3 - K
                     DO I = 1,IRTOP
                        JG(I,K,K) = PI(K,K,I)/R(I,IM)
                        JF(I,K,K) = QI(K,K,I)/R(I,IM)/C
                        JG(I,KC,K) = PI(KC,K,I)/R(I,IM)
                        JF(I,KC,K) = QI(KC,K,I)/R(I,IM)/C
                     END DO
                  END DO
C
C-----------------------------------------------------------------------
                  IF ( IWRREGWF.NE.0 ) THEN
C solution 1
                     WRITE (IFILSS,REC=IKM1+(IT-1)*NKM) IT,'REG',IKM1,
     &                      IRTOP,NSOL,IKM1,IKM2,
     &                      ((ZG(I,K,1),I=1,IRTOP),K=1,NSOL),
     &                      ((ZF(I,K,1),I=1,IRTOP),K=1,NSOL)
C
C solution 2
                     WRITE (IFILSS,REC=IKM2+(IT-1)*NKM) IT,'REG',IKM2,
     &                      IRTOP,NSOL,IKM1,IKM2,
     &                      ((ZG(I,K,2),I=1,IRTOP),K=1,NSOL),
     &                      ((ZF(I,K,2),I=1,IRTOP),K=1,NSOL)
                  END IF
C
                  IF ( IWRIRRWF.NE.0 ) THEN
C solution 1
                     WRITE (IFILSS,REC=IKM1+(IT-1+NT)*NKM) IT,'IRR',
     &                      IKM1,IRTOP,((JG(I,K,1),I=1,IRTOP),K=1,NSOL),
     &                      ((JF(I,K,1),I=1,IRTOP),K=1,NSOL)
C
C solution 2
                     WRITE (IFILSS,REC=IKM2+(IT-1+NT)*NKM) IT,'IRR',
     &                      IKM2,IRTOP,((JG(I,K,2),I=1,IRTOP),K=1,NSOL),
     &                      ((JF(I,K,2),I=1,IRTOP),K=1,NSOL)
C
                  END IF
C================================================= COUPLING = END ===
               END IF
C
C
C
CALCULATE SUM OF INTEGRALS TO BE MULTIPLIED TO   TAU(K1,K2)
               DO K1 = 1,NSOL
                  DO K2 = 1,NSOL
C
                     LIN = LIN + 1
                     TSSLIN(LIN) = TSST2(K1,K2)
                     SSSLIN(LIN) = ZG(1,K1,K2)/CJL1
                     IF ( CALCINT ) THEN
C REGULAR PART    Z*Z
C
                        CALL CINTABWR(ZG(1,1,K1),ZG(1,1,K2),ZGZG,
     &                                ZF(1,1,K1),ZF(1,1,K2),ZFZF,
     &                                R2DRDI_W_RADINT(1,IM),NSOL,NSOL,
     &                                IRTOP,NRMAX)
C
                        CALL SUMUPINT(DZZ(LIN),F1,ZGZG,R1M,F1,ZFZF,R1M,
     &                                NSOL)
                        CALL SUMUPINT(SZZ(LIN),F1,ZGZG,CGG,-F1,ZFZF,CGF,
     &                                NSOL)
                        CALL SUMUPINT(OZZ(LIN),F1,ZGZG,CFG,-F1,ZFZF,CFF,
     &                                NSOL)
C
                        CALL SUMUPINT(TZZ(LIN),F1,ZGZG,CTG,-F1,ZFZF,CTF,
     &                                NSOL)
                        CALL SUMUPINT(KZZ(LIN),F1,ZGZG,RKD,F1,ZFZF,RKD,
     &                                NSOL)
C
C-----------------------------------------------------------------------
                        IF ( IHYPER.EQ.1 .AND. Z(IT).NE.0 ) THEN
C
                           CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),ZF(1,1,K1),
     &                                  ZG(1,1,K2),ZF(1,1,K2),RMEHF,
     &                                  NSOL,NSOL,DRDI(1,IM))
C
                           IF ( FINITE_NUCLEUS ) THEN
C calculates integrals inside nucleus but up to now only
C approximately because jlim is not the nuclear radius
C the same arguments are valid for the irregular parts below
                              CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),
     &                           ZF(1,1,K1),ZG(1,1,K2),ZF(1,1,K2),
     &                           RMEHF1,NSOL,NSOL,DRDI(1,IM))
                              CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),
     &                           ZF(1,1,K1),ZG(1,1,K2),ZF(1,1,K2),
     &                           RMEHF2,NSOL,NSOL,DRDI(1,IM))
                              DO I = 1,NSOL
                                 DO J = 1,NSOL
                                    RMEHF(I,J) = RMEHF(I,J)
     &                                 - RMEHF1(I,J) + RMEHF2(I,J)
                                 END DO
                              END DO
                           END IF
C
                           CALL SUMUPINT(BZZ(LIN),B_AU2CGS,RMEHF,CHF,
     &                        0.0D0,RMEHF,CHF,NSOL)
C
                        ELSE
C
                           BZZ(LIN) = 0D0
C
                        END IF
C-----------------------------------------------------------------------
C
C IRREGULAR PART    Z*J
C THE  ENDING  A (B)  STANDS FOR THE DOMINATING (DIMINATED)
C SET OF SPIN-ANGULAR-CHAR:  I.E.  J==J(A,A)  FOR R>RMT
C
                        IF ( K1.EQ.K2 ) THEN
C
                           CALL CINTABWR(ZG(1,1,K1),JG(1,1,K1),ZGJG,
     &                        ZF(1,1,K1),JF(1,1,K1),ZFJF,
     &                        R2DRDI_W_RADINT(1,IM),NSOL,NSOL,IRTOP,
     &                        NRMAX)
C
                           CALL SUMUPINT(DZJ(LIN),F1,ZGJG,R1M,F1,ZFJF,
     &                        R1M,NSOL)
                           CALL SUMUPINT(SZJ(LIN),F1,ZGJG,CGG,-F1,ZFJF,
     &                        CGF,NSOL)
                           CALL SUMUPINT(OZJ(LIN),F1,ZGJG,CFG,-F1,ZFJF,
     &                        CFF,NSOL)
C
                           CALL SUMUPINT(TZJ(LIN),F1,ZGJG,CTG,-F1,ZFJF,
     &                        CTF,NSOL)
                           CALL SUMUPINT(KZJ(LIN),F1,ZGJG,RKD,F1,ZFJF,
     &                        RKD,NSOL)
C
C-----------------------------------------------------------------------
                           IF ( IHYPER.EQ.1 .AND. Z(IT).NE.0 ) THEN
C
                              CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),
     &                           ZF(1,1,K1),JG(1,1,K1),JF(1,1,K1),RMEHF,
     &                           NSOL,NSOL,DRDI(1,IM))
C
                              IF ( FINITE_NUCLEUS ) THEN
C calculates integrals inside nucleus but up to now only
C approximately because jlim is not the nuclear radius
C the same arguments are valid for the irregular parts below
                                 CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),
     &                              ZF(1,1,K1),JG(1,1,K1),JF(1,1,K1),
     &                              RMEHF1,NSOL,NSOL,DRDI(1,IM))
                                 CALL CINTHFF(IM,IRTOP,ZG(1,1,K1),
     &                              ZF(1,1,K1),JG(1,1,K1),JF(1,1,K1),
     &                              RMEHF2,NSOL,NSOL,DRDI(1,IM))
                                 DO I = 1,NSOL
                                    DO J = 1,NSOL
                                       RMEHF(I,J) = RMEHF(I,J)
     &                                    - RMEHF1(I,J) + RMEHF2(I,J)
                                    END DO
                                 END DO
                              END IF
C
                              CALL SUMUPINT(BZJ(LIN),B_AU2CGS,RMEHF,CHF,
     &                           0.0D0,RMEHF,CHF,NSOL)
                           END IF
                        END IF
C-----------------------------------------------------------------------
                     END IF
C           ! OF IF (.CALCINT.)
                  END DO
               END DO
C
Check WRONSKI-relationship
C
               IF ( WRONSKI ) THEN
                  WRITE (6,99003) IT,L,NINT(2*MJ)
                  NSTEP = 20
                  NVIEW = 10
                  I = 0
 5                CONTINUE
                  IF ( I.LT.NVIEW .OR. I.GE.(IRTOP-NVIEW) ) THEN
                     I = I + 1
                  ELSE IF ( I.LT.(IRTOP-NVIEW-NSTEP) ) THEN
                     I = I + NSTEP
                  ELSE
                     I = IRTOP - NVIEW
                  END IF
                  IF ( I.LE.IRTOP ) THEN
                     CRSQ = (1.0D0+ERYD/CSQR)*C*R(I,IM)**2
C
                     WRON(1,1) = ZF(I,1,1)*JG(I,1,1) - ZG(I,1,1)
     &                           *JF(I,1,1) + ZF(I,2,1)*JG(I,2,1)
     &                           - ZG(I,2,1)*JF(I,2,1)
                     IF ( NSOL.EQ.2 ) THEN
                        WRON(2,2) = ZF(I,1,2)*JG(I,1,2) - ZG(I,1,2)
     &                              *JF(I,1,2) + ZF(I,2,2)*JG(I,2,2)
     &                              - ZG(I,2,2)*JF(I,2,2)
                        WRON(2,1) = ZF(I,1,2)*JG(I,1,1) - ZG(I,1,2)
     &                              *JF(I,1,1) + ZF(I,2,2)*JG(I,2,1)
     &                              - ZG(I,2,2)*JF(I,2,1)
                        WRON(2,1) = ZF(I,1,1)*JG(I,1,2) - ZG(I,1,1)
     &                              *JF(I,1,2) + ZF(I,2,1)*JG(I,2,2)
     &                              - ZG(I,2,1)*JF(I,2,2)
                     END IF
                     WRITE (6,99004) R(I,IM),I,
     &                               ((WRON(K1,K2)*CRSQ,K1=1,NSOL),K2=1,
     &                               NSOL)
                     GOTO 5
                  END IF
               END IF
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         DO LIN = 1,NLIN_T(IT)
            I1 = IKM1LIN(LIN)
            I2 = IKM2LIN(LIN)
            TSST(I1,I2,IT) = TSSLIN(LIN)
            SSST(I1,I2,IT) = SSSLIN(LIN)
         END DO
C
         MSST(1:NKM_T(IT),1:NKM_T(IT),IT)
     &      = TSST(1:NKM_T(IT),1:NKM_T(IT),IT)
C
         CALL ZGETRF(NKM_T(IT),NKM_T(IT),MSST(1,1,IT),NKMMAX,IPIV,INFO)
         CALL ZGETRI(NKM_T(IT),MSST(1,1,IT),NKMMAX,IPIV,MAUX,
     &               NKMMAX*NKMMAX,INFO)
C
         DO LIN = 1,NLIN_T(IT)
            I1 = IKM1LIN(LIN)
            I2 = IKM2LIN(LIN)
            MEZZ(I1,I2,IT,IDOS) = DZZ(LIN)
            MEZJ(I1,I2,IT,IDOS) = DZJ(LIN)
            MEZZ(I1,I2,IT,ISMT) = SZZ(LIN)
            MEZJ(I1,I2,IT,ISMT) = SZJ(LIN)
            MEZZ(I1,I2,IT,IOMT) = OZZ(LIN)
            MEZJ(I1,I2,IT,IOMT) = OZJ(LIN)
            MEZZ(I1,I2,IT,IHFF) = BZZ(LIN)
            MEZJ(I1,I2,IT,IHFF) = BZJ(LIN)
            MEZZ(I1,I2,IT,ISDM) = TZZ(LIN)
            MEZJ(I1,I2,IT,ISDM) = TZJ(LIN)
            MEZZ(I1,I2,IT,IKDS) = KZZ(LIN)
            MEZJ(I1,I2,IT,IKDS) = KZJ(LIN)
         END DO
C
         IF ( IPRINT.GE.3 ) THEN
            WRITE (6,99006) 'atom type   ',IT,ERYD,CTL(IT,1)
C
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,3,3,1,
     &                      1.0D-9,6)
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,6)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,6)
         END IF
         IF ( DUMP_SSITE ) THEN
            IF ( IT.EQ.1 ) THEN
               LFILNAM = LEN_TRIM(SOLVER)
               FILNAM = 'zzzzzz_ssite_REL_'//SOLVER(1:LFILNAM)
               LFILNAM = LFILNAM + 17
               IF ( GF_CONV_RH ) THEN
                  FILNAM = FILNAM(1:LFILNAM)//'_RH.dat'
               ELSE
                  FILNAM = FILNAM(1:LFILNAM)//'_ZJ.dat'
               END IF
               CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:(LFILNAM+7)))
            END IF
            WRITE (IOTMP,99006) 'atom type   ',IT,ERYD,CTL(IT,1),SOLVER,
     &                          GF_CONV_RH
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,3,3,1,
     &                      1.0D-9,IOTMP)
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,IOTMP)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,IOTMP)
            DO I = 1,NKM
               DO J = 1,NKM
                  WRITE (IOTMP,99005) I,J,TSST(I,J,IT),MEZZ(I,J,IT,1),
     &                                MEZJ(I,J,IT,1)
               END DO
            END DO
C
            IF ( IT.EQ.ITTOP ) THEN
               CLOSE (IOTMP)
               WRITE (6,99007) FILNAM(1:(LFILNAM+7))
            END IF
         END IF
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
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      IF ( WRONSKI ) CALL STOP_REGULAR(ROUTINE,'WRONSKI - test done ')
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) CALL STOP_BREAKPOINT(ROUTINE)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C
99001 FORMAT (/,2X,'SSITE:   ERYD = ',2F10.6)
99002 FORMAT (2X,A,I3,4X,A,10F10.4)
99003 FORMAT (/,' IT=',I2,' L=',I2,' MJ=',I2,'/2')
99004 FORMAT (F11.8,I4,1X,4(2F15.12,:,2X))
99005 FORMAT (2I3,2X,'TSST ',2E25.15,/,8X,'MEZZ ',2E25.15,/,8X,'MEZJ ',
     &        2E25.15,/)
99006 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),//,10X,'energy',12X,
     &        2F12.6,/,10X,'speed of light',f16.6,/,10X,
     &        'solver            ',A,/,10X,'GF_CONV_RH        ',L10,/)
99007 FORMAT (//,5X,'>>>>>  DUMP written to file  ',A,//)
99008 FORMAT (//,10X,'problem with SOLVER = ',A,/)
      END
