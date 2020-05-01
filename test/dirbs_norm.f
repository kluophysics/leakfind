C*==dirbs_norm.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DIRBS_NORM(IT,C,E,P,V,B,IRTOP,PR,QR,PH,QH,DP,DQ,
     &                      NCWF_SOL_SPH,IKM_CWFSOL_SPH,ZG_SPH,ZF_SPH,
     &                      JG_SPH,JF_SPH,TSST0,SSST0,NSOLSPH)
C   ********************************************************************
C   *                                                                  *
C   *      OUTPUT: normalized wave functions                           *
C   *                                                                  *
C   *       normalized   regular   (PR,QR) -> J - ip SUM H * t         *
C   *       normalized irregular   (PH,QH) -> H                        *
C   *                                                                  *
C   *      solutions to the single site problem for                    *
C   *      spherical potential functions V(r) and B(r)                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX,DRDI,DRDIOVR,R,NPAN,JRCUT
      USE MOD_TYPES,ONLY:Z,NLT,IMT
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_ANGMOM,ONLY:IKM1LIN,IKM2LIN,NKMMAX
      IMPLICIT NONE
C*--DIRBS_NORM22
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 E,P
      INTEGER IRTOP,IT,NSOLSPH
      COMPLEX*16 B(NRMAX),DP(2,2,NRMAX),DQ(2,2,NRMAX),
     &           JF_SPH(NRMAX,2,NKMMAX),JG_SPH(NRMAX,2,NKMMAX),
     &           PH(2,2,NRMAX),PR(2,2,NRMAX),QH(2,2,NRMAX),QR(2,2,NRMAX)
     &           ,SSST0(NKMMAX,NKMMAX),TSST0(NKMMAX,NKMMAX),V(NRMAX),
     &           ZF_SPH(NRMAX,2,NKMMAX),ZG_SPH(NRMAX,2,NKMMAX)
      INTEGER IKM_CWFSOL_SPH(2,NKMMAX),NCWF_SOL_SPH(NKMMAX)
C
C Local variables
C
      COMPLEX*16 A(2,2),ARG,B1,B2,CHL,CHLB1,CHLB2,CJL,CJL1,CJLB1,CJLB2,
     &           CNL,CNLB1,CNLB2,CRSQ,CSUM,DET,F11,F12,F21,F22,G11,G12,
     &           G21,G22,GAM(2,2),GAMINV(2,2),MAUX(2,2),MSST2(2,2),PFAC,
     &           R_NORM,SIG(2,2),TSST2(2,2),WRON(2,2),XSST2(2,2)
      REAL*8 CG1,CG2,CG4,CG5,CG8,CSQR,MJ,RTOP,SK1,SK2
      COMPLEX*16 CJLZ,CNLZ
      INTEGER I,I1,I2,I3,IKM,IKM1,IKM2,IKMCB(2),IM,INFO,IPIV(2),IR,ISK1,
     &        ISK2,ISOL,J,K,K1,K2,KAP1,KAP2,L,L1,LB1,LB2,LIN,MUM05,
     &        NSTEP,NVIEW
      INTEGER IKAPMUE
      LOGICAL WRONSKI
C
C*** End of declarations rewritten by SPAG
C
      WRONSKI = .FALSE.
C
      IM = IMT(IT)
      RTOP = R(IRTOP,IM)
C
      CSQR = C*C
      PFAC = P/(1.0D0+E/CSQR)
C
      TSST0(:,:) = C0
      SSST0(:,:) = C0
C
      LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L = 0,(NLT(IT)-1)
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         DO MUM05 = -L - 1, + L
            MJ = DBLE(MUM05) + 0.5D0
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
            IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
            IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
            IKMCB(1) = IKM1
            IKMCB(2) = IKM2
C
            CG1 = -MJ/(KAP1+0.5D0)
            CG5 = -MJ/(-KAP1+0.5D0)
C
            IF ( ABS(MJ).GT.L ) THEN
               NSOLSPH = 1
C
               CG2 = 0.0D0
               CG4 = 0.0D0
               CG8 = 0.0D0
C
            ELSE
C
               NSOLSPH = 2
C
               CG2 = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CG4 = -MJ/(KAP2+0.5D0)
               CG8 = -MJ/(-KAP2+0.5D0)
            END IF
C
            DO ISOL = 1,NSOLSPH
               IKM = IKMCB(ISOL)
               NCWF_SOL_SPH(IKM) = NSOLSPH
               IKM_CWFSOL_SPH(1:NSOLSPH,IKM) = IKMCB(1:NSOLSPH)
            END DO
C
C-----------------------------------------------------------------------
C        unnormalized   regular   (PR,QR)    (temporary)
C          normalized irregular   (PH,QH) -> H
C
            CALL CINIT(2*2*NRMAX,PR)
            CALL CINIT(2*2*NRMAX,QR)
            CALL CINIT(2*2*NRMAX,PH)
            CALL CINIT(2*2*NRMAX,QH)
C
            CALL DIRBS(.TRUE.,C,E,L,MJ,KAP1,KAP2,P,CG1,CG2,CG4,CG5,CG8,
     &                 V,B,Z(IT),R(1,IM),DRDI(1,IM),DRDIOVR(1,IM),IRTOP,
     &                 JRCUT(0,IM),NPAN(IM),PR,QR,PH,QH,DP,DQ,NRMAX,
     &                 .TRUE.)
C
            ISK1 = ISIGN(1,KAP1)
            ISK2 = ISIGN(1,KAP2)
            L1 = L
            LB1 = L - ISK1
            LB2 = L - ISK2
C
            ARG = P*RTOP
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
C-----------------------------------------------------------------------
C
C  wavefunctions at the muffin-tin-radius
C
            G11 = PR(1,1,IRTOP)/RTOP
            G12 = PR(1,2,IRTOP)/RTOP
            G21 = PR(2,1,IRTOP)/RTOP
            G22 = PR(2,2,IRTOP)/RTOP
            F11 = QR(1,1,IRTOP)/(RTOP*C)
            F12 = QR(1,2,IRTOP)/(RTOP*C)
            F21 = QR(2,1,IRTOP)/(RTOP*C)
            F22 = QR(2,2,IRTOP)/(RTOP*C)
C
C -------------------------------------------------------------------
C       T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
C -------------------------------------------------------------------
C
            SK1 = DBLE(ISK1)
            SK2 = DBLE(ISK2)
C
            CNL = (CHL-CJL)/CI
            CNLB1 = (CHLB1-CJLB1)/CI
            CNLB2 = (CHLB2-CJLB2)/CI
C
            GAM(1,1) = CJL*C*F11 - PFAC*SK1*CJLB1*G11
            GAM(2,1) = CJL*C*F21 - PFAC*SK2*CJLB2*G21
            GAM(1,2) = CJL*C*F12 - PFAC*SK1*CJLB1*G12
            GAM(2,2) = CJL*C*F22 - PFAC*SK2*CJLB2*G22
C
            SIG(1,1) = CNL*C*F11 - PFAC*SK1*CNLB1*G11
            SIG(2,1) = CNL*C*F21 - PFAC*SK2*CNLB2*G21
            SIG(1,2) = CNL*C*F12 - PFAC*SK1*CNLB1*G12
            SIG(2,2) = CNL*C*F22 - PFAC*SK2*CNLB2*G22
C
            CALL CMATCOP(NSOLSPH,2,GAM,GAMINV)
            CALL ZGETRF(NSOLSPH,NSOLSPH,GAMINV,2,IPIV,INFO)
            CALL ZGETRI(NSOLSPH,GAMINV,2,IPIV,MAUX,2*2,INFO)
C
            DO I2 = 1,NSOLSPH
               DO I1 = 1,NSOLSPH
                  CSUM = 0.0D0
                  DO I3 = 1,NSOLSPH
                     CSUM = CSUM + SIG(I1,I3)*GAMINV(I3,I2)
                  END DO
                  XSST2(I1,I2) = P*CSUM
               END DO
            END DO
C
            DO I1 = 1,NSOLSPH
               DO I2 = 1,NSOLSPH
                  MSST2(I1,I2) = -XSST2(I1,I2)
               END DO
               MSST2(I1,I1) = MSST2(I1,I1) + CI*P
            END DO
C
            CALL CMATCOP(NSOLSPH,2,MSST2,TSST2)
            CALL ZGETRF(NSOLSPH,NSOLSPH,TSST2,2,IPIV,INFO)
            CALL ZGETRI(NSOLSPH,TSST2,2,IPIV,MAUX,2*2,INFO)
C
            IF ( NSOLSPH.EQ.1 ) THEN
C====================================================================
C NO COUPLING TO OTHER SCATTERING CHANNELS
C
               R_NORM = (CJL-CI*P*CHL*TSST2(1,1))/G11
C
               DO IR = 1,IRTOP
                  ZG_SPH(IR,1,IKM1) = PR(1,1,IR)*R_NORM
                  ZF_SPH(IR,1,IKM1) = QR(1,1,IR)*R_NORM
               END DO
C
C============================================== NO COUPLING = END ===
            ELSE
C====================================================================
C COUPLING OF TWO SCATTERING CHANNELS
C   Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
C              INDEX 2: BOUNDARY CONDITION
C
               DET = G11*G22 - G12*G21
C
COEFFICIENTS TO GET:   R(K1,K1)  R(K2,K1)
               B1 = CJL - CI*P*CHL*TSST2(1,1)
               B2 = -CI*P*CHL*TSST2(2,1)
               A(1,1) = (G22*B1-G12*B2)/DET
               A(2,1) = (G11*B2-G21*B1)/DET
C
COEFFICIENTS TO GET:   R(K1,K2)  R(K2,K2)
               B1 = -CI*P*CHL*TSST2(1,2)
               B2 = CJL - CI*P*CHL*TSST2(2,2)
               A(1,2) = (G22*B1-G12*B2)/DET
               A(2,2) = (G11*B2-G21*B1)/DET
C
CALCULATE FUNCTIONS: R(K1,K1), R(K2,K1), R(K1,K2), R(K2,K2)
               DO IR = 1,IRTOP
                  DO K = 1,NSOLSPH
C
                     IKM = IKMCB(K)
C
                     ZG_SPH(IR,1,IKM) = PR(1,1,IR)*A(1,K) + PR(1,2,IR)
     &                                  *A(2,K)
                     ZG_SPH(IR,2,IKM) = PR(2,1,IR)*A(1,K) + PR(2,2,IR)
     &                                  *A(2,K)
C
                     ZF_SPH(IR,1,IKM) = QR(1,1,IR)*A(1,K) + QR(1,2,IR)
     &                                  *A(2,K)
                     ZF_SPH(IR,2,IKM) = QR(2,1,IR)*A(1,K) + QR(2,2,IR)
     &                                  *A(2,K)
C
                  END DO
               END DO
C
C================================================= COUPLING = END ===
            END IF
C
            DO I = 1,NSOLSPH
               DO J = 1,NSOLSPH
                  IKM = IKMCB(J)
                  LIN = LIN + 1
C
                  I1 = IKM1LIN(LIN)
                  I2 = IKM2LIN(LIN)
                  TSST0(I1,I2) = TSST2(I,J)
                  SSST0(I1,I2) = ZG_SPH(1,I,IKM)/CJL1
C
                  DO IR = 1,IRTOP
                     JG_SPH(IR,I,IKM) = PH(I,J,IR)
                     JF_SPH(IR,I,IKM) = QH(I,J,IR)
                  END DO
C
               END DO
            END DO
C
Check WRONSKI-relationship
C
            IF ( WRONSKI ) THEN
C
               CRSQ = (1.0D0+E/CSQR)*P*CI
C
               WRITE (6,99001) IT,L,NINT(2*MJ)
               NSTEP = 20
               NVIEW = 10
               NVIEW = 235
               IR = 0
 10            CONTINUE
               IF ( IR.LT.NVIEW .OR. IR.GE.(IRTOP-NVIEW) ) THEN
                  IR = IR + 1
               ELSE IF ( IR.LT.(IRTOP-NVIEW-NSTEP) ) THEN
                  IR = IR + NSTEP
               ELSE
                  IR = IRTOP - NVIEW
               END IF
               IF ( IR.LE.IRTOP ) THEN
C
                  WRON(1,1) = QR(1,1,IR)*PR(1,1,IR) - PR(1,1,IR)
     &                        *QR(1,1,IR) + QR(2,1,IR)*PR(2,1,IR)
     &                        - PR(2,1,IR)*QR(2,1,IR)
                  WRON(2,2) = QR(1,2,IR)*PR(1,2,IR) - PR(1,2,IR)
     &                        *QR(1,2,IR) + QR(2,2,IR)*PR(2,2,IR)
     &                        - PR(2,2,IR)*QR(2,2,IR)
                  WRON(2,1) = QR(1,2,IR)*PR(1,1,IR) - PR(1,2,IR)
     &                        *QR(1,1,IR) + QR(2,2,IR)*PR(2,1,IR)
     &                        - PR(2,2,IR)*QR(2,1,IR)
                  WRON(1,2) = QR(1,1,IR)*PR(1,2,IR) - PR(1,1,IR)
     &                        *QR(1,2,IR) + QR(2,1,IR)*PR(2,2,IR)
     &                        - PR(2,1,IR)*QR(2,2,IR)
C
                  WRITE (6,99002) R(IR,IM),IR,
     &                            ((WRON(K1,K2)*CRSQ,K1=1,NSOLSPH),K2=1,
     &                            NSOLSPH)
                  IF ( ABS(1D0-R(IR+1,IM)/R(IR,IM)).LT.1D-7 )
     &                 WRITE (6,*)
                  GOTO 10
               END IF
            END IF
C
         END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      IF ( WRONSKI ) STOP '<DIRBS_NORM>     WRONSKI - test done '
C
99001 FORMAT (/,' IT=',I2,' L=',I2,' MJ=',I2,'/2')
99002 FORMAT (F11.8,I4,1X,4(2F15.12,:,2X))
      END
