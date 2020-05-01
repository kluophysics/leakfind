C*==negf_tmat_sa.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NEGF_TMAT_SA(IFIL_ZPLS,IFIL_ZMIN,IFILBAR_ZPLS,
     &                        IFILBAR_ZMIN,DELT1_LBAR,DELT2_LBAR,
     &                        IQBOT_LBAR,IQTOP_LBAR,DELT1_RBAR,
     &                        DELT2_RBAR,IQBOT_RBAR,IQTOP_RBAR,
     &                        BAR_FLAG_Q,LBAR_FLAG_Q,RBAR_FLAG_Q,VBAR,
     &                        TSS_LBAR,TSS_RBAR,TSST,MEZZ,MEZZ_BAR,
     &                        MEZZMP)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the t-matrix                                          *
C   *                                                                  *
C   *   Delta t(z,z') = < R(z) | Delta V | R_bar(z') >                 *
C   *                                                                  *
C   *  for all combinations of energy arguments z and z'               *
C   *                                                                  *
C   *                                                                  *
C   *   MEZZMP = < R(-) | A | R(+) >  for A = 1, Sigma_z               *
C   *                                                                  *
C   *                                                                  *
C   *  NOTE: the routine assumes that the RH-convention is used        *
C   *  for the Green's function                                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:ITOQ,NOQ,NQ
      USE MOD_ANGMOM,ONLY:NLQ,NKMMAX,NKM,NMEMAX,IKM1LIN,IKM2LIN,WKM1,
     &    WKM2,WKM3,AME_G,NCPLWF
      USE MOD_RMESH,ONLY:R2DRDI,JRWS,NRMAX
      USE MOD_TYPES,ONLY:IMT,NTMAX,VT,BT,NCPLWFMAX,IKMCPLWF
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_CALCMODE,ONLY:GF_CONV_RH
      IMPLICIT NONE
C*--NEGF_TMAT_SA34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NEGF_TMAT_SA')
      INTEGER ISMT
      PARAMETER (ISMT=2)
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
      INTEGER CHECK_INTEGRATION
      PARAMETER (CHECK_INTEGRATION=0)
      LOGICAL CHECK_DT1,CHECK_T1_VS_T2
      PARAMETER (CHECK_DT1=.TRUE.,CHECK_T1_VS_T2=.TRUE.)
C
C Dummy arguments
C
      INTEGER IFILBAR_ZMIN,IFILBAR_ZPLS,IFIL_ZMIN,IFIL_ZPLS,IQBOT_LBAR,
     &        IQBOT_RBAR,IQTOP_LBAR,IQTOP_RBAR
      REAL*8 VBAR
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 DELT1_LBAR(NKMMAX,NKMMAX,2,2,IQBOT_LBAR:IQTOP_LBAR),
     &           DELT1_RBAR(NKMMAX,NKMMAX,2,2,IQBOT_RBAR:IQTOP_RBAR),
     &           DELT2_LBAR(NKMMAX,NKMMAX,2,2,IQBOT_LBAR:IQTOP_LBAR),
     &           DELT2_RBAR(NKMMAX,NKMMAX,2,2,IQBOT_RBAR:IQTOP_RBAR),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ_BAR(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX),
     &           TSS_LBAR(NKMMAX,NKMMAX,IQBOT_LBAR:IQTOP_LBAR),
     &           TSS_RBAR(NKMMAX,NKMMAX,IQBOT_RBAR:IQTOP_RBAR)
C
C Local variables
C
      COMPLEX*16 BZBARZ(:,:),BZZBAR(:,:),DDT1(:,:),DZZMP(:,:),JF(:,:,:),
     &           JFBAR(:,:,:),JFPLS(:,:,:),JG(:,:,:),JGBAR(:,:,:),
     &           JGPLS(:,:,:),RDDT1(:,:),SZZMP(:,:),VZBARZ(:,:),
     &           VZZBAR(:,:),ZF(:,:,:),ZFBAR(:,:,:),ZFBARZF(2,2),
     &           ZFPLS(:,:,:),ZFZFBAR(2,2),ZFZFPLS(2,2),ZG(:,:,:),
     &           ZGBAR(:,:,:),ZGBARZG(2,2),ZGPLS(:,:,:),ZGZGBAR(2,2),
     &           ZGZGPLS(2,2)
      REAL*8 CGF(2,2),CGG(2,2),DB_R2DRDI(:),DV_R2DRDI(:),MJ
      INTEGER I,I1,I2,IA_ERR,IFIL,IFILBAR,IKM1,IKM2,IKMI,IKMJ,IKM_SOL(2)
     &        ,IM,IMKMI,IMKM_SOL(2),IO,IQ,IRTOP,IT,IZ,IZBAR,J,K1,K2,
     &        KAP1,KAP2,KK,KKA,KKB,L,LAM1,LAM2,LIN,MJM05,NSOL
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DZZMP,SZZMP,VZZBAR,VZBARZ,DV_R2DRDI,DB_R2DRDI
      ALLOCATABLE BZBARZ,BZZBAR,DDT1,RDDT1
      ALLOCATABLE ZF,ZG,ZGBAR,ZFBAR,ZGPLS,ZFPLS
      ALLOCATABLE JF,JG,JGBAR,JFBAR,JGPLS,JFPLS
C
      ALLOCATE (DV_R2DRDI(NRMAX))
      ALLOCATE (DZZMP(NKMMAX,NKMMAX),SZZMP(NKMMAX,NKMMAX))
      ALLOCATE (VZZBAR(NKMMAX,NKMMAX),VZBARZ(NKMMAX,NKMMAX))
      ALLOCATE (BZZBAR(NKMMAX,NKMMAX),BZBARZ(NKMMAX,NKMMAX))
      ALLOCATE (DB_R2DRDI(NRMAX),DDT1(NKMMAX,NKMMAX))
      ALLOCATE (RDDT1(NKMMAX,NKMMAX))
      ALLOCATE (JFBAR(1,1,1))
      ALLOCATE (JGBAR(1,1,1))
      ALLOCATE (ZFBAR(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGBAR(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JFPLS(1,1,1))
      ALLOCATE (JGPLS(1,1,1))
      ALLOCATE (ZFPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JF(1,1,1))
      ALLOCATE (JG(1,1,1))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IF ( CHECK_INTEGRATION.EQ.1 ) THEN
         IFIL_ZPLS = IFILBAR_ZPLS
         IFIL_ZMIN = IFILBAR_ZMIN
      ELSE IF ( CHECK_INTEGRATION.EQ.2 ) THEN
         IFILBAR_ZPLS = IFIL_ZPLS
         IFILBAR_ZMIN = IFIL_ZMIN
      END IF
C
      WRITE (6,*) 'CHECK_INTEGRATION.EQ.',CHECK_INTEGRATION
C
      IF ( .NOT.GF_CONV_RH ) STOP 
     &                         '<NEGF_TMAT>: GF convention has to be RH'
C
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      LOOP_IQ:DO IQ = IQBOT_LBAR,IQTOP_RBAR
         LOOP_IO:DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
Csw            IRTOP = JRMT(IM)
C
            IF ( CHECK_INTEGRATION.EQ.0 ) THEN
C
Csw               DV_R2DRDI(1:IRTOP) = (VBAR-VT(1:IRTOP,IT))
               DV_R2DRDI(1:IRTOP) = (VT(1:IRTOP,IT)-VBAR)
     &                              *R2DRDI(1:IRTOP,IM)
               DB_R2DRDI(1:IRTOP) = -BT(1:IRTOP,IT)*R2DRDI(1:IRTOP,IM)
C
            ELSE
               DV_R2DRDI(1:IRTOP) = R2DRDI(1:IRTOP,IM)
               DB_R2DRDI(1:IRTOP) = R2DRDI(1:IRTOP,IM)
            END IF
C
            LOOP_IZ:DO IZ = 1,2
C
               IF ( IZ.EQ.1 ) THEN
                  IFIL = IFIL_ZPLS
               ELSE
                  IFIL = IFIL_ZMIN
               END IF
C
               CALL WAVFUN_READ_REL(IFIL,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                              IKMCPLWF)
C
               LOOP_IZBAR:DO IZBAR = 1,2
                  IF ( IZBAR.EQ.1 ) THEN
                     IFILBAR = IFILBAR_ZPLS
                  ELSE
                     IFILBAR = IFILBAR_ZMIN
                  END IF
C
                  CALL WAVFUN_READ_REL(IFILBAR,IT,0,ZGBAR,ZFBAR,JGBAR,
     &                                 JFBAR,IRTOP,NCPLWF,IKMCPLWF)
C
                  IF ( IZ.EQ.2 .AND. IZBAR.EQ.2 )
     &                 CALL WAVFUN_READ_REL(IFIL_ZPLS,IT,0,ZGPLS,ZFPLS,
     &                 JGPLS,JFPLS,IRTOP,NCPLWF,IKMCPLWF)
C
                  IF ( CHECK_INTEGRATION.EQ.3 ) THEN
                     ZG(:,:,:) = ZGPLS(:,:,:)
                     ZF(:,:,:) = ZFPLS(:,:,:)
                  END IF
C
                  WKM1(:,:) = C0
                  WKM2(:,:) = C0
                  WKM3(:,:) = C0
                  DDT1(:,:) = C0
                  RDDT1(:,:) = C0
                  VZZBAR(:,:) = C0
                  VZBARZ(:,:) = C0
                  BZZBAR(:,:) = C0
                  BZBARZ(:,:) = C0
                  DZZMP(:,:) = C0
                  SZZMP(:,:) = C0
C
                  LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
                  DO L = 0,(NLQ(IQ)-1)
C
                     KAP1 = -L - 1
                     KAP2 = L
                     IF ( L.EQ.0 ) KAP2 = KAP1
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                     DO MJM05 = -L - 1, + L
                        MJ = DBLE(MJM05) + 0.5D0
C
                        IF ( ABS(MJ).GE.DBLE(L) ) THEN
                           NSOL = 1
                        ELSE
                           NSOL = 2
                        END IF
C-----------------------------------------------------------------------
                        IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
                        IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C-----------------------------------------------------------------------
                        IKM_SOL(1) = IKM1
                        IKM_SOL(2) = IKM2
                        IMKM_SOL(1) = IKAPMUE(-KAP1,NINT(MJ-0.5D0))
                        IMKM_SOL(2) = IKAPMUE(-KAP2,NINT(MJ-0.5D0))
C
C-----------------------------------------------------------------------
C
                        CALL RINIT(4,CGF)
C
                        DO I = 1,NSOL
                           IKMI = IKM_SOL(I)
                           DO J = 1,NSOL
                              IKMJ = IKM_SOL(J)
                              CGG(I,J) = AME_G(IKMI,IKMJ,2,ISMT)
                           END DO
C
                           IMKMI = IMKM_SOL(I)
                           CGF(I,I) = AME_G(IMKMI,IMKMI,2,ISMT)
                        END DO
C
C-----------------------------------------------------------------------
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
                        DO K1 = 1,NSOL
                           LAM1 = IKM_SOL(K1)
C
                           DO K2 = 1,NSOL
                              LAM2 = IKM_SOL(K2)
                              LIN = LIN + 1
C
                              I1 = IKM1LIN(LIN)
                              I2 = IKM2LIN(LIN)
C
                              IF ( BAR_FLAG_Q(IQ) ) THEN
C
C***********************************************************************
C                            <Z(z)|Delta V|Z_bar(z')>
C***********************************************************************
C
                                 CALL CINTABR(ZG(1,1,LAM1),
     &                              ZGBAR(1,1,LAM2),ZGZGBAR,ZF(1,1,LAM1)
     &                              ,ZFBAR(1,1,LAM2),ZFZFBAR,DV_R2DRDI,
     &                              NSOL,NSOL,IRTOP,NRMAX)
C
                                 DO KK = 1,NSOL
                                    VZZBAR(I1,I2) = VZZBAR(I1,I2)
     &                                 + ZGZGBAR(KK,KK) + ZFZFBAR(KK,KK)
                                 END DO
C
                                 CALL CINTABR(ZG(1,1,LAM1),
     &                              ZGBAR(1,1,LAM2),ZGZGBAR,ZF(1,1,LAM1)
     &                              ,ZFBAR(1,1,LAM2),ZFZFBAR,DB_R2DRDI,
     &                              NSOL,NSOL,IRTOP,NRMAX)
C
                                 DO KKA = 1,NSOL
                                    DO KKB = 1,NSOL
                                       BZZBAR(I1,I2) = BZZBAR(I1,I2)
     &                                    + CGG(KKA,KKB)
     &                                    *ZGZGBAR(KKA,KKB)
     &                                    - CGF(KKA,KKB)
     &                                    *ZFZFBAR(KKA,KKB)
                                    END DO
                                 END DO
C
C***********************************************************************
C                            <Z_bar(z)|Delta V|Z(z')>
C***********************************************************************
C
                                 CALL CINTABR(ZGBAR(1,1,LAM1),
     &                              ZG(1,1,LAM2),ZGBARZG,ZFBAR(1,1,LAM1)
     &                              ,ZF(1,1,LAM2),ZFBARZF,DV_R2DRDI,
     &                              NSOL,NSOL,IRTOP,NRMAX)
C
                                 DO KK = 1,NSOL
                                    VZBARZ(I1,I2) = VZBARZ(I1,I2)
     &                                 + ZGBARZG(KK,KK) + ZFBARZF(KK,KK)
                                 END DO
C
                                 CALL CINTABR(ZGBAR(1,1,LAM1),
     &                              ZG(1,1,LAM2),ZGBARZG,ZFBAR(1,1,LAM1)
     &                              ,ZF(1,1,LAM2),ZFBARZF,DB_R2DRDI,
     &                              NSOL,NSOL,IRTOP,NRMAX)
C
                                 DO KKA = 1,NSOL
                                    DO KKB = 1,NSOL
                                       BZBARZ(I1,I2) = BZBARZ(I1,I2)
     &                                    + CGG(KKA,KKB)
     &                                    *ZGBARZG(KKA,KKB)
     &                                    - CGF(KKA,KKB)
     &                                    *ZFBARZF(KKA,KKB)
                                    END DO
                                 END DO
C***********************************************************************
C                            <Z(-)|Z(+)>
C***********************************************************************
C
                              ELSE IF ( IZ.EQ.2 .AND. IZBAR.EQ.2 ) THEN
C
                                 CALL CINTABR(ZG(1,1,LAM1),
     &                              ZGPLS(1,1,LAM2),ZGZGPLS,ZF(1,1,LAM1)
     &                              ,ZFPLS(1,1,LAM2),ZFZFPLS,
     &                              R2DRDI(1,IM),NSOL,NSOL,IRTOP,NRMAX)
C
                                 DO KK = 1,NSOL
                                    DZZMP(I1,I2) = DZZMP(I1,I2)
     &                                 + ZGZGPLS(KK,KK) + ZFZFPLS(KK,KK)
                                 END DO
C
                                 DO KKA = 1,NSOL
                                    DO KKB = 1,NSOL
                                       SZZMP(I1,I2) = SZZMP(I1,I2)
     &                                    + CGG(KKA,KKB)
     &                                    *ZGZGPLS(KKA,KKB)
     &                                    - CGF(KKA,KKB)
     &                                    *ZFZFPLS(KKA,KKB)
                                    END DO
                                 END DO
C
                              END IF
C
                           END DO
                        END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
                     END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                  END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
                  IF ( BAR_FLAG_Q(IQ) ) THEN
C
                     WKM1(:,:) = VZZBAR(:,:) + BZZBAR(:,:)
                     WKM2(:,:) = VZBARZ(:,:) + BZBARZ(:,:)
C
                     IF ( LBAR_FLAG_Q(IQ) ) THEN
C
                        DELT1_LBAR(:,:,IZ,IZBAR,IQ) = WKM1(:,:)
                        DELT2_LBAR(:,:,IZBAR,IZ,IQ) = WKM2(:,:)
C
Csw                        WKM3(:,:) = TSS_LBAR(:,:,IQ) - TSST(:,:,IT)
                        WKM3(:,:) = TSST(:,:,IT) - TSS_LBAR(:,:,IQ)
C
                     ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
C
                        DELT1_RBAR(:,:,IZ,IZBAR,IQ) = WKM1(:,:)
                        DELT2_RBAR(:,:,IZBAR,IZ,IQ) = WKM2(:,:)
C
Csw                        WKM3(:,:) = TSS_RBAR(:,:,IQ) - TSST(:,:,IT)
                        WKM3(:,:) = TSST(:,:,IT) - TSS_RBAR(:,:,IQ)
C
                     ELSE
                        STOP '<NEGF_DELT1>:  L/RBAR_FLAG_Q ??? '
                     END IF
C
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                     check consistency of data sets
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
                     IF ( CHECK_DT1 .OR. CHECK_T1_VS_T2 ) THEN
                        WRITE (100,*) '########### ',IQ,IZ,IZBAR
                        CALL CMATSTRUCT('DT(1) numerical ',WKM1,NKM,
     &                                  NKMMAX,3,3,0,TOL,100)
C
                     END IF
C
                     IF ( CHECK_T1_VS_T2 ) THEN
                        WRITE (101,*) '########### ',IQ,IZ,IZBAR
                        CALL CMATSTRUCT('DT(2) numerical ',WKM2,NKM,
     &                                  NKMMAX,3,3,0,TOL,101)
                     END IF
C
                     IF ( CHECK_DT1 ) THEN
C
                        IF ( IZ.EQ.1 .AND. IZBAR.EQ.1 ) THEN
C
                           DDT1 = WKM1(:,:) - WKM3(:,:)
                           DO I = 1,NKM
                              DO J = 1,NKM
                                 IF ( ABS(DREAL(WKM3(I,J))).GT.0D0 )
     &                                THEN
                                    RDDT1(I,J)
     &                                 = DCMPLX((DREAL(DDT1(I,J))
     &                                 /DREAL(WKM3(I,J))),
     &                                 (DIMAG(DDT1(I,J))
     &                                 /DIMAG(WKM3(I,J))))
                                 ELSE IF ( ABS(DREAL(DDT1(I,J)))
     &                              .GT.1D-12 ) THEN
                                    WRITE (6,*) '<NEGF_TMAT>: SYMMETRY?'
                                    RDDT1(I,J)
     &                                 = DCMPLX(DREAL(DDT1(I,J)),
     &                                 DIMAG(DDT1(I,J)))
                                 END IF
                              END DO
                           END DO
C
                           CALL CMATSTRUCT('DT(1) numerical ',WKM1,NKM,
     &                        NKMMAX,3,3,0,TOL,6)
                           CALL CMATSTRUCT('TSS - TBAR      ',WKM3,NKM,
     &                        NKMMAX,3,3,0,TOL,6)
                           CALL CMATSTRUCT('DDT(1)',DDT1,NKM,NKMMAX,3,3,
     &                        0,TOL,6)
                           CALL CMATSTRUCT('relative error',RDDT1,NKM,
     &                        NKMMAX,3,3,0,TOL,6)
C
                           IF ( LBAR_FLAG_Q(IQ) ) THEN
                              CALL CMATSTRUCT('TBAR  LEFT      ',
     &                           TSS_LBAR(1,1,IQ),NKM,NKMMAX,3,3,0,TOL,
     &                           6)
                           ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
                              CALL CMATSTRUCT('TBAR  RIGHT     ',
     &                           TSS_RBAR(1,1,IQ),NKM,NKMMAX,3,3,0,TOL,
     &                           6)
                           END IF
C
                           CALL CMATSTRUCT('TSS             ',
     &                        TSST(1,1,IT),NKM,NKMMAX,3,3,0,TOL,6)
C
                        END IF
C
                        IF ( CHECK_INTEGRATION.EQ.1 ) THEN
                           WKM1(:,:) = MEZZ_BAR(:,:,IT,1)
                           WKM2(:,:) = VZZBAR(:,:) - WKM1(:,:)
                           CALL CMATSTRUCT('VZZBAR',VZZBAR,NKM,NKMMAX,3,
     &                        3,0,1D-8,6)
                           CALL CMATSTRUCT('MEZZ_BAR_1',WKM1,NKM,NKMMAX,
     &                        3,3,0,1D-8,6)
                           CALL CMATSTRUCT('VZZBAR - MEZZ_BAR_1',WKM2,
     &                        NKM,NKMMAX,3,3,0,1D-12,6)
C
                           WKM1(:,:) = MEZZ_BAR(:,:,IT,2)
                           WKM2(:,:) = BZZBAR(:,:) - WKM1(:,:)
                           CALL CMATSTRUCT('BZZBAR',BZZBAR,NKM,NKMMAX,3,
     &                        3,0,1D-8,6)
                           CALL CMATSTRUCT('MEZZ_BAR_2',WKM1,NKM,NKMMAX,
     &                        3,3,0,1D-8,6)
                           CALL CMATSTRUCT('BZZBAR - MEZZ_BAR_2',WKM2,
     &                        NKM,NKMMAX,3,3,0,1D-12,6)
C
                        ELSE IF ( CHECK_INTEGRATION.EQ.2 ) THEN
C
                           WKM1(:,:) = MEZZ(:,:,IT,1)
                           WKM2(:,:) = VZZBAR(:,:) - WKM1(:,:)
                           CALL CMATSTRUCT('VZZBAR',VZZBAR,NKM,NKMMAX,3,
     &                        3,0,1D-8,6)
                           CALL CMATSTRUCT('MEZZ_1',WKM1,NKM,NKMMAX,3,3,
     &                        0,1D-8,6)
                           CALL CMATSTRUCT('VZZBAR - MEZZ_1',WKM2,NKM,
     &                        NKMMAX,3,3,0,1D-12,6)
C
                           WKM1(:,:) = MEZZ(:,:,IT,2)
                           WKM2(:,:) = BZZBAR(:,:) - WKM1(:,:)
                           CALL CMATSTRUCT('BZZBAR',BZZBAR,NKM,NKMMAX,3,
     &                        3,0,1D-8,6)
                           CALL CMATSTRUCT('MEZZ_2',WKM1,NKM,NKMMAX,3,3,
     &                        0,1D-8,6)
                           CALL CMATSTRUCT('BZZBAR - MEZZ_2',WKM2,NKM,
     &                        NKMMAX,3,3,0,1D-12,6)
C
C
                        END IF
                     END IF
C
                  ELSE IF ( IZ.EQ.2 .AND. IZBAR.EQ.2 ) THEN
C
                     MEZZMP(:,:,IT,1) = DZZMP(:,:)
                     MEZZMP(:,:,IT,2) = SZZMP(:,:)
C
                     IF ( CHECK_INTEGRATION.EQ.3 ) THEN
                        CALL CMATSTRUCT('MEZZ_1',MEZZ(:,:,IT,1),NKM,
     &                                  NKMMAX,3,3,0,1D-8,16)
                        CALL CMATSTRUCT('MEZZMP_1',MEZZMP(:,:,IT,1),NKM,
     &                                  NKMMAX,3,3,0,1D-8,17)
                        CALL CMATSTRUCT('MEZZ_2',MEZZ(:,:,IT,2),NKM,
     &                                  NKMMAX,3,3,0,1D-8,16)
                        CALL CMATSTRUCT('MEZZMP_2',MEZZMP(:,:,IT,2),NKM,
     &                                  NKMMAX,3,3,0,1D-8,17)
C                       STOP 'TESTING MEZZPP(1/2) vs. MEZZ(1/2)'
                     END IF
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
                  END IF
C
               END DO LOOP_IZBAR
C
            END DO LOOP_IZ
C
         END DO LOOP_IO
      END DO LOOP_IQ
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      IF ( CHECK_T1_VS_T2 ) THEN
         DO IQ = IQBOT_LBAR,IQTOP_RBAR
            IF ( LBAR_FLAG_Q(IQ) ) THEN
C
               WRITE (6,*) ' '
               WRITE (6,*) 'LEFT barrier IQ=',IQ
               CALL NEGF_TCHK('1-11','1-22',DELT1_LBAR(1,1,1,1,IQ),
     &                        DELT1_LBAR(1,1,2,2,IQ),' *L')
               CALL NEGF_TCHK('2-11','2-22',DELT2_LBAR(1,1,1,1,IQ),
     &                        DELT2_LBAR(1,1,2,2,IQ),' *L')
C
               CALL NEGF_TCHK('1-11','2-11',DELT1_LBAR(1,1,1,1,IQ),
     &                        DELT2_LBAR(1,1,1,1,IQ),'T  ')
               CALL NEGF_TCHK('1-22','2-22',DELT1_LBAR(1,1,2,2,IQ),
     &                        DELT2_LBAR(1,1,2,2,IQ),'T  ')
C
               CALL NEGF_TCHK('1-12','1-21',DELT1_LBAR(1,1,1,2,IQ),
     &                        DELT1_LBAR(1,1,2,1,IQ),' *L')
               CALL NEGF_TCHK('2-12','2-21',DELT2_LBAR(1,1,1,2,IQ),
     &                        DELT2_LBAR(1,1,2,1,IQ),' *L')
C
               CALL NEGF_TCHK('1-12','2-12',DELT1_LBAR(1,1,1,2,IQ),
     &                        DELT2_LBAR(1,1,1,2,IQ),'T*L')
               CALL NEGF_TCHK('1-12','2-21',DELT1_LBAR(1,1,1,2,IQ),
     &                        DELT2_LBAR(1,1,2,1,IQ),'T  ')
C
            ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
C
               WRITE (6,*) ' '
               WRITE (6,*) 'RIGHT barrier IQ=',IQ
               CALL NEGF_TCHK('1-11','1-22',DELT1_RBAR(1,1,1,1,IQ),
     &                        DELT1_RBAR(1,1,2,2,IQ),' *L')
               CALL NEGF_TCHK('2-11','2-22',DELT2_RBAR(1,1,1,1,IQ),
     &                        DELT2_RBAR(1,1,2,2,IQ),' *L')
C
               CALL NEGF_TCHK('1-11','2-11',DELT1_RBAR(1,1,1,1,IQ),
     &                        DELT2_RBAR(1,1,1,1,IQ),'T  ')
               CALL NEGF_TCHK('1-22','2-22',DELT1_RBAR(1,1,2,2,IQ),
     &                        DELT2_RBAR(1,1,2,2,IQ),'T  ')
C
               CALL NEGF_TCHK('1-12','1-21',DELT1_RBAR(1,1,1,2,IQ),
     &                        DELT1_RBAR(1,1,2,1,IQ),' *L')
               CALL NEGF_TCHK('2-12','2-21',DELT2_RBAR(1,1,1,2,IQ),
     &                        DELT2_RBAR(1,1,2,1,IQ),' *L')
C
               CALL NEGF_TCHK('1-12','2-12',DELT1_RBAR(1,1,1,2,IQ),
     &                        DELT2_RBAR(1,1,1,2,IQ),'T*L')
               CALL NEGF_TCHK('1-12','2-21',DELT1_RBAR(1,1,1,2,IQ),
     &                        DELT2_RBAR(1,1,2,1,IQ),'T  ')
C
            END IF
         END DO
C        STOP 'TTTTTT     <NEGF_TMAT>:  CHECK_T1_VS_T2    TTTTTT'
      END IF
C
      DEALLOCATE (ZF,ZG,ZGBAR,ZFBAR,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc <NEGF_DELT1>  -> ZFB'
      IF ( CHECK_INTEGRATION.NE.0 )
     &      STOP '<NEGF_DELT1> integration checked'
C
      END
C*==negf_tchk.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NEGF_TCHK(IZ1,IZ2,A,B,C)
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,WKM2,WKM3
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--NEGF_TCHK581
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*3 C
      CHARACTER*4 IZ1,IZ2
      COMPLEX*16 A(NKMMAX,NKMMAX),B(NKMMAX,NKMMAX)
C
C Local variables
C
      REAL*8 A1,A2,ADIFF,DEL,PWL,RDEL
      INTEGER I,IKM,J,L
C
C*** End of declarations rewritten by SPAG
C
      WKM1 = B
C
      IF ( C(1:1).EQ.'T' ) WKM1 = TRANSPOSE(WKM1)
      IF ( C(2:2).EQ.'*' ) WKM1 = DCONJG(WKM1)
C
      IF ( C(3:3).EQ.'L' ) THEN
C
         WKM2 = C0
         WKM3 = C0
         IKM = 0
         DO L = 0,12
            PWL = (-1D0)**L
            DO I = 1,2*(2*L+1)
               IKM = IKM + 1
               WKM2(IKM,IKM) = PWL
            END DO
            IF ( IKM.GE.NKM ) EXIT
         END DO
C
         IF ( IPRINT.GE.1 ) THEN
            CALL CMATSTRUCT('A',A(:,:),NKM,NKMMAX,3,3,0,1D-8,7777)
            CALL CMATSTRUCT('B',B(:,:),NKM,NKMMAX,3,3,0,1D-8,7777)
            CALL CMATSTRUCT('WKM1',WKM2(:,:),NKM,NKMMAX,3,3,0,1D-8,7777)
            CALL CMATSTRUCT('L REL',WKM2(:,:),NKM,NKMMAX,3,3,0,1D-8,
     &                      7777)
Csw            CALL CHANGEREP(NKM,NKMMAX,WKM2,'REL>RLM',L2)
Csw            CALL  CMATSTR('L RLM',L2(:,:),NKM,NKMMAX,2,2,0,1D-8,7777)
         END IF
C
         WKM3(1:NKM,1:NKM) = MATMUL(WKM2(1:NKM,1:NKM),WKM1(1:NKM,1:NKM))
         WKM1(1:NKM,1:NKM) = MATMUL(WKM3(1:NKM,1:NKM),WKM2(1:NKM,1:NKM))
C
      END IF
C
      RDEL = 0D0
      DEL = 0D0
      DO I = 1,NKM
         DO J = 1,NKM
            ADIFF = ABS(WKM1(I,J)-A(I,J))
            A1 = ABS(WKM1(I,J))
            A2 = ABS(A(I,J))
            IF ( ADIFF.GT.1D-12 ) THEN
               DEL = MAX(DEL,ADIFF)
               RDEL = MAX(RDEL,ADIFF/MAX(A1,A2))
            END IF
         END DO
      END DO
C
      WRITE (6,99001) IZ1,IZ2,C,DEL,RDEL
C
99001 FORMAT ('CHECK  Del T',A,'    Del T',A,'   ',A,3X,2F12.8)
C
      END
