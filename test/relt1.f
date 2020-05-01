C*==relt1.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RELT1(TAUT,CALCINT,GETIRRSOL,IPRINT,TSST,MSST,MEZZ,
     &                 MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *    calculates    relativistic T1-time                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRWS
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_SITES,ONLY:IQAT
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_ENERGY,ONLY:EFERMI,EIMAG
      USE MOD_TYPES,ONLY:NTMAX,IMT,Z,BT,NT
      USE MOD_ANGMOM,ONLY:LTAB,NMEMAX,NKMAX,NKMMAX,IKM2LIN,IKM1LIN,NKMQ,
     &    NLQ,NKQ,NKM,NK,NL,KAPTAB,LBTAB,NPOL
      USE MOD_CONSTANTS,ONLY:B_AU2CGS,HBAR_CGS,KB_CGS,PI,RY_ERG
      IMPLICIT NONE
C*--RELT120
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER IPRINT
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 AME1(NKMMAX,NKMMAX,NPOL),AME2(NKMMAX,NKMMAX,NPOL),BSUM,
     &       BT1(NKMMAX,NKMMAX),DOSK(NKMAX),DOSL(NKMAX),GAMN,I21,I30,J1,
     &       J2,L1,LB2,MJ,MJ1,MJ2,MJIREP,PREFAC,RAT,RATDOS(NKMAX),
     &       RATKK(NKMAX,NKMAX),RATL(NKMAX),RMEDOS,RTOT,T01,T23
      COMPLEX*16 CEFERMI,ME(NKMMAX,NKMMAX,NPOL),P,
     &           SSST(NKMMAX,NKMMAX,NTMAX),W2(NKMMAX,NKMMAX)
      REAL*8 CGC_RACAH
      CHARACTER*4 CHK(10),CHL(10)
      INTEGER I,I1,I2,IBTF(NKMAX),IBTI(NKMAX),IK,IK2,IKM1,IKM2,IKMMAX,
     &        IL,IL0,IL1,IL2,IL3,IQ,IT,J,K,K0,K1,K2,K3,KAP1,KAP2,KAPS,L,
     &        LAM0,LAM1,LAM2,LAM3,LIN,M0,M1,M2,M3,MJ1M05,MJ2M05,MJM05,
     &        MM,NN
      LOGICAL NONMAG
C
C*** End of declarations rewritten by SPAG
C
C ......................................................................
C
      DATA CHL/' s  ',' p  ',' d  ',' f  ',' g  ','    ','    ','    ',
     &     '    ','    '/
      DATA CHK/'s1/2','p1/2','p3/2','d3/2','d5/2','f5/2','f7/2','g7/2',
     &     'g9/2','h9/2'/
C
C ......................................................................
C
      WRITE (6,99001)
C
C --------------- index table for representative hyperfine field mj=+1/2
      MJIREP = 0.5D0
C
C-----------------------------------------------------------------------
      BSUM = 0.0D0
      DO IT = 1,NT
         DO J = 1,JRWS(IMT(IT))
            BSUM = BSUM + ABS(BT(J,IT))
         END DO
      END DO
      IF ( ABS(BSUM/DBLE(NT)).GT.0.00001D0 ) THEN
         NONMAG = .FALSE.
      ELSE
         NONMAG = .TRUE.
      END IF
C
C ------------------------------------------------ set up angular matrix
      AME1(:,:,:) = 0D0
      AME2(:,:,:) = 0D0
C
      I1 = 0
      DO K1 = 1,NK
         KAP1 = KAPTAB(K1)
         L1 = LTAB(K1)
         J1 = ABS(KAP1) - 0.5D0
C
         DO MJ1M05 = NINT(-J1-0.5D0),NINT(J1-0.5D0)
            MJ1 = DBLE(MJ1M05) + 0.5D0
            I1 = I1 + 1
            I2 = 0
C
            DO K2 = 1,NK
               KAP2 = KAPTAB(K2)
               LB2 = LBTAB(K2)
               KAPS = KAP1 + KAP2
               J2 = ABS(KAP2) - 0.5D0
C
               DO MJ2M05 = NINT(-J2-0.5D0),NINT(J2-0.5D0)
                  MJ2 = DBLE(MJ2M05) + 0.5D0
                  I2 = I2 + 1
                  AME1(I1,I2,2) = 0.0D0
                  IF ( ABS(MJ2-MJ1-1.0D0).LT..00001D0 ) THEN
                     IF ( KAP1.EQ.KAP2 ) THEN
                        AME1(I1,I2,2) = 4.0D0*KAP1/(4.0D0*KAP1**2-1.0D0)
     &                                  *SQRT(KAP1**2-(MJ1+0.5D0)**2)
                     ELSE IF ( ABS(KAPS).EQ.1 ) THEN
                        AME1(I1,I2,2) = -1.0D0/(2.0D0*KAP1-KAPS)
     &                                  *SQRT((KAP1-KAPS*(MJ1+0.5D0))
     &                                  *(KAP1-KAPS*(MJ1+1.5D0)))
                     END IF
                  END IF
C
                  W2(I1,I2) = SQRT(2.0D0*(2D0*LB2+1D0)/(2D0*L1+1D0))
     &                        *CGC_RACAH(LB2,1.0D0,L1,0.0D0,0.0D0,0.0D0)
     &                        *(SQRT(2.0D0)
     &                        *CGC_RACAH(L1,0.5D0,J1,(MJ1+0.5D0),-0.5D0,
     &                        MJ1)*CGC_RACAH(LB2,0.5D0,J2,(MJ2-0.5D0),
     &                        +0.5D0,MJ2)
     &                        *CGC_RACAH(LB2,1.0D0,L1,(MJ1+0.5D0),0.0D0,
     &                        (MJ1+0.5D0))
     &                        -CGC_RACAH(L1,0.5D0,J1,(MJ1-0.5D0),+0.5D0,
     &                        MJ1)*CGC_RACAH(LB2,0.5D0,J2,(MJ2-0.5D0),
     &                        +0.5D0,MJ2)
     &                        *CGC_RACAH(LB2,1.0D0,L1,(MJ1+0.5D0),
     &                        -1.0D0,(MJ1-0.5D0))
     &                        +CGC_RACAH(L1,0.5D0,J1,(MJ1+0.5D0),-0.5D0,
     &                        MJ1)*CGC_RACAH(LB2,0.5D0,J2,(MJ2+0.5D0),
     &                        -0.5D0,MJ2)
     &                        *CGC_RACAH(LB2,1.0D0,L1,(MJ1+1.5D0),
     &                        -1.0D0,(MJ1+0.5D0)))
                  IF ( ABS(AME1(I1,I2,2)).GT.1D-8 ) THEN
                     IF ( ABS(AME1(I1,I2,2)-W2(I1,I2)).GT.1D-6 )
     &                    WRITE (6,'(A,2I4,4F10.5)')
     &                     '*** PROBLEM with A for (I1,I2)',I1,I2,
     &                    AME1(I1,I2,2),DREAL(W2(I1,I2))
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C
      DO J = 1,NKMMAX
         DO I = 1,NKMMAX
            AME2(I,J,2) = AME1(I,J,2)
         END DO
      END DO
C
      IPRINT = 1
      IF ( IPRINT.GE.1 ) THEN
         NN = 2*NL**2
         MM = NKMMAX
         CALL RMATSTRUCT('A1(-) ',AME1(1,1,2),NN,MM,3,3,0,1D-8,6)
         CALL RMATSTRUCT('A2(-) ',AME2(1,1,2),NN,MM,3,3,0,1D-8,6)
      END IF
      IPRINT = 0
C
      CEFERMI = DCMPLX(EFERMI,EIMAG)
C
      WRITE (6,'(A10,2F10.6)') '   EFERMI ',CEFERMI
      IF ( NONMAG ) THEN
         WRITE (6,'(10X,A)') 'paramagnetic system'
      ELSE
         WRITE (6,'(10X,A)') 'ferromagnetic system'
      END IF
C
C
C ------------------------------ calculate WF and DOS - radial integrals
      CALL SSITE(1,1,IFILCBWF,CALCINT,GETIRRSOL,CEFERMI,P,IPRINT,NKM,
     &           TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         IKMMAX = NKMQ(IQ)
C
         WRITE (6,99003) IT
C
C ======================================================================
C                            matrix elements
C ======================================================================
C
         CALL RELT1ME(1,IKMMAX,IFILCBWF,'REG',1,IKMMAX,IFILCBWF,'REG',
     &                IT,AME1,AME2,ME,MJIREP,IBTF,IBTI,BT1)
C
         CALL CMATSTRUCT('MEHF  ',ME(1,1,2),NN,MM,3,3,0,1D-8,6)
C
         WRITE (6,99002)
C
C ----------------------------------------------------------- initialize
         RTOT = 0.0D0
         DO IK = 1,NK
            DOSK(IK) = 0.0D0
            DOSL(IK) = 0.0D0
            RATL(IK) = 0.0D0
            RATDOS(IK) = 0.0D0
C
            DO IK2 = 1,NK
               RATKK(IK,IK2) = 0.0D0
            END DO
         END DO
C
C -------------------------- store component projected tau-matrix in  W2
C
         W2(:,:) = TAUT(:,:,IT)
C
         CALL TABGAMMAN(Z(IT),GAMN)
C
         PREFAC = 4D0*PI*KB_CGS*HBAR_CGS*(GAMN*B_AU2CGS/(2D0*PI*RY_ERG))
     &            **2
C
         LAM0 = 0
         DO K0 = 1,NKQ(IQ)
            IL0 = LTAB(K0) + 1
            DO M0 = 1,2*ABS(KAPTAB(K0))
               LAM0 = LAM0 + 1
C
               LAM1 = 0
               DO K1 = 1,NKQ(IQ)
                  IL1 = LTAB(K1) + 1
                  DO M1 = 1,2*ABS(KAPTAB(K1))
                     LAM1 = LAM1 + 1
C
                     LAM2 = 0
                     DO K2 = 1,NKQ(IQ)
                        IL2 = LTAB(K2) + 1
                        DO M2 = 1,2*ABS(KAPTAB(K2))
                           LAM2 = LAM2 + 1
C
                           LAM3 = 0
                           DO K3 = 1,NKQ(IQ)
                              IL3 = LTAB(K3) + 1
                              DO M3 = 1,2*ABS(KAPTAB(K3))
                                 LAM3 = LAM3 + 1
C
                                 T01 = DIMAG(W2(LAM0,LAM1))
                                 T23 = DIMAG(W2(LAM2,LAM3))
                                 I30 = DREAL(ME(LAM3,LAM0,2))
                                 I21 = DREAL(ME(LAM2,LAM1,2))
C
                                 RAT = PREFAC*T01*T23*I30*I21
C
                                 IF ( IL0.EQ.IL1 .AND. IL1.EQ.IL2 .AND. 
     &                                IL2.EQ.IL3 ) RATL(IL0) = RATL(IL0)
     &                                + RAT
C
                                 IF ( K0.EQ.K1 .AND. K2.EQ.K3 )
     &                                RATKK(K0,K2) = RATKK(K0,K2) + RAT
C
                                 RTOT = RTOT + RAT
C
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         LIN = 0
         DO L = 0,NLQ(IQ) - 1
            IL = L + 1
C
            DO MJM05 = -L - 1, + L
               MJ = DBLE(MJM05) + 0.5D0
C
               DO K1 = 1,2
                  DO K2 = 1,2
                     LIN = LIN + 1
                     IKM1 = IKM1LIN(LIN)
                     IKM2 = IKM2LIN(LIN)
                     IF ( K1.EQ.K2 .AND. K1.EQ.1 ) THEN
                        IK = 2*L + 1
                        DOSK(IK) = DOSK(IK)
     &                             - DIMAG(W2(IKM1,IKM2)*MEZZ(IKM1,IKM2,
     &                             IT,1)/PI)
     &                             + DIMAG(MEZJ(IKM1,IKM2,IT,1)/PI)
                     END IF
                     IF ( K1.EQ.K2 .AND. K1.EQ.2 ) THEN
                        IK = 2*L
                        DOSK(IK) = DOSK(IK)
     &                             - DIMAG(W2(IKM1,IKM2)*MEZZ(IKM1,IKM2,
     &                             IT,1)/PI)
     &                             + DIMAG(MEZJ(IKM1,IKM2,IT,1)/PI)
                     END IF
C
                     DOSL(IL) = DOSL(IL)
     &                          - DIMAG(W2(IKM1,IKM2)*MEZZ(IKM1,IKM2,IT,
     &                          1)/PI)
C
                     IF ( K1.EQ.K2 ) DOSL(IL) = DOSL(IL)
     &                    + DIMAG(MEZJ(IKM1,IKM2,IT,1)/PI)
                     IF ( ABS(MJ).GT.DBLE(L) ) GOTO 20
                  END DO
               END DO
C
C
 20         END DO
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
         DO K = 1,NKQ(IQ)
            IL = K/2 + 1
            J1 = ABS(KAPTAB(K)) - 0.5D0
C
            IF ( IL.EQ.1 ) THEN
               RMEDOS = BT1(IBTF(K),IBTI(K))/(-2.D0*B_AU2CGS/3.D0)
            ELSE
               RMEDOS = BT1(IBTF(K),IBTI(K))
     &                  /(2.D0*B_AU2CGS/(2D0*KAPTAB(K)+2D0))
            END IF
C
            RATDOS(IL) = RATDOS(IL) + PREFAC*PI**2*(2*J1+1)
     &                   /(6*J1*(J1+1))*(DOSK(K)*RMEDOS)**2
C
            IF ( 2*(IL-1).EQ.K ) THEN
               RMEDOS = BT1(IBTF(K),IBTI(K+1))
     &                  /(2.D0*B_AU2CGS/(KAPTAB(K)+KAPTAB(K+1)+2))
               RATDOS(IL) = RATDOS(IL)
     &                      + PREFAC*PI**2*1.0D0/(3.0D0*(J1+1))*DOSK(K)
     &                      *DOSK(K+1)*RMEDOS**2
C
            END IF
         END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         WRITE (6,99004) (CHK(K2),K2=1,NKQ(IQ))
C
         DO K1 = 1,NKQ(IQ)
            WRITE (6,99005) CHK(K1),
     &                      (BT1(IBTF(K1),IBTI(K2))*1D-06,K2=1,NKQ(IQ))
         END DO
C
         WRITE (6,99006) (CHK(K2),K2=1,NKQ(IQ))
         DO K1 = 1,NKQ(IQ)
            WRITE (6,99005) CHK(K1),(RATKK(K1,K2),K2=1,NKQ(IQ))
         END DO
C
         WRITE (6,99007) (DOSK(IK),IK=1,NKQ(IQ))
C
         WRITE (6,99008) (CHL(IL),IL=1,NLQ(IQ))
         WRITE (6,99009) (DOSL(IL),IL=1,NLQ(IQ))
         WRITE (6,99010) (RATL(IL),IL=1,NLQ(IQ))
         WRITE (6,99011) (RATDOS(IL),IL=1,NLQ(IQ))
C
         WRITE (6,99012) RTOT
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*        *****   ******  *              *******    **        *'
     &  ,/,10X,
     &  '*        *    *  *       *                 *      * *        *'
     &  ,/,10X,
     &  '*        *    *  *       *                 *     *  *        *'
     &  ,/,10X,
     &  '*        *****   *****   *        ***      *        *        *'
     &  ,/,10X,
     &  '*        *  *    *       *                 *        *        *'
     &  ,/,10X,
     &  '*        *   *   *       *                 *        *        *'
     &  ,/,10X,
     &  '*        *    *  ******  ******            *        *        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,10X,'spin-flip hyperfine - matrix elements evaluated',/)
99003 FORMAT (/,' nuclear spin lattice relaxation rate',
     &        '  R = 1/T_1*T  for component',I3,/,1X,67('='),/)
99004 FORMAT (/,' hyperfine field  B   [MG]:',3X,7(8X,A4))
99005 FORMAT (27X,A4,7F12.4)
99006 FORMAT (/,' K-resolved rate   [1/s*K]:',3X,7(8X,A4))
99007 FORMAT (' K-resolved DOS    [1/Ry]: ',4X,7F12.4)
99008 FORMAT (/,' l-resolved quantities     ',11X,A4,7(14X,A4,6X))
99009 FORMAT (' DOS               [1/Ry]: ',4X,F12.4,7(6X,F12.4,6X))
99010 FORMAT (' rate              [1/s*K]:',4X,F12.4,7(6X,F12.4,6X))
99011 FORMAT (' rate (DOS-form.)  [1/s*K]:',4X,F12.4,7(6X,F12.4,6X))
C99012 FORMAT (' hyperfine mat.element in [au]:',7F12.4)
99012 FORMAT (' ------------------------------------------',/,
     &        ' R  total          [1/s*K]:    ',F12.4,/)
      END
C*==relt1me.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RELT1ME(IKMI1,IKMI2,IFILI,WFTI,IKMF1,IKMF2,IFILF,WFTF,
     &                   IT,AME1,AME2,ME,MJIREP,IBTF,IBTI,BT1)
C   ********************************************************************
C   *                                                                  *
C   * calculate radial matrix element for relat. spin lattice rate     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX
      USE MOD_ANGMOM,ONLY:NKM,NKMAX,NKMMAX,NPOL
      USE MOD_RMESH,ONLY:R2DRDI,DRDI,JRWS,NRMAX
      USE MOD_CONSTANTS,ONLY:B_AU2CGS
      IMPLICIT NONE
C*--RELT1ME412
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFILF,IFILI,IKMF1,IKMF2,IKMI1,IKMI2,IT
      REAL*8 MJIREP
      CHARACTER*3 WFTF,WFTI
      REAL*8 AME1(NKMMAX,NKMMAX,NPOL),AME2(NKMMAX,NKMMAX,NPOL),
     &       BT1(NKMMAX,NKMMAX)
      INTEGER IBTF(NKMAX),IBTI(NKMAX)
      COMPLEX*16 ME(NKMMAX,NKMMAX,NPOL)
C
C Local variables
C
      COMPLEX*16 CINTF(NRMAX),CINTG(NRMAX),FRFR,GRGR,NORMF,NORMI,
     &           RMEHF(NCPLWFMAX,NCPLWFMAX),ZFF(NRMAX,2),ZFI(NRMAX,2),
     &           ZGF(NRMAX,2),ZGI(NRMAX,2)
      INTEGER IKMF(2),IKMI(2),IKMTF,IKMTI,IM,IPOL,IPOL1,IPOL2,IR,IRTOP,
     &        ITF,ITI,JF,JI,K,KAPF(2),KAPI(2),KF,KI,KKF,KKI,LF,LI,NSOLF,
     &        NSOLI
      REAL*8 MJF,MJI,RME(NKMMAX,NKMMAX)
      CHARACTER*3 STR3
C
C*** End of declarations rewritten by SPAG
C
C conversion factor for hyperfine fields from A.U. to GAUSS
C                                      electron charge     in esu
C                                      bohr-radius         in cm
C-----------------------------------------------------------------------
      IPOL1 = 2
      IPOL2 = 2
C
      CALL CINIT(NKMMAX*NKMMAX*NPOL,ME)
      CALL RINIT(NKMMAX*NKMMAX,BT1)
      CALL RINIT(NKMMAX*NKMMAX,RME)
C-----------------------------------------------------------------------
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
C
      DO IKMTI = IKMI1,IKMI2
         READ (IFILI,REC=IKMTI+(IT-1)*NKM) ITI,LI,MJI,NSOLI,STR3,
     &         (KAPI(K),IKMI(K),(ZGI(IR,K),ZFI(IR,K),IR=1,IRTOP),K=1,
     &         NSOLI)
         IF ( IT.NE.ITI ) STOP ' IT(INI) <> IT in <RELT1>'
         IF ( STR3.NE.WFTI ) STOP 'WFT(INI) not consistent in <RELT1>'
C
         IF ( IKMTI.NE.IKMI(1) ) WRITE (6,*)
     &                                   'IIII #######################',
     &                                  IKMTI,IKMI(1)
         DO IR = 1,IRTOP
            CINTG(IR) = ZGI(IR,1)*ZGI(IR,1)*R2DRDI(IR,IM)
            CINTF(IR) = ZFI(IR,1)*ZFI(IR,1)*R2DRDI(IR,IM)
         END DO
C
         CALL CRADINT(IM,CINTG,GRGR)
         CALL CRADINT(IM,CINTF,FRFR)
C
         NORMI = SQRT(GRGR+FRFR)
         KKI = LI + ABS(KAPI(1))
C
         DO IKMTF = IKMF1,IKMF2
            READ (IFILF,REC=IKMTF+(IT-1)*NKM) ITF,LF,MJF,NSOLF,STR3,
     &            (KAPF(K),IKMF(K),(ZGF(IR,K),ZFF(IR,K),IR=1,IRTOP),K=1,
     &            NSOLF)
            IF ( IT.NE.ITF ) STOP 'IT(FIN) <> IT in <RELT1>'
            IF ( STR3.NE.WFTF ) STOP 
     &                              'WFT(FIN) not consistent in <RELT1>'
C
            IF ( IKMTF.NE.IKMF(1) ) WRITE (6,*)
     &            'FFFF #######################',IKMTF,IKMF(1)
            DO IR = 1,IRTOP
               CINTG(IR) = ZGF(IR,1)*ZGF(IR,1)*R2DRDI(IR,IM)
               CINTF(IR) = ZFF(IR,1)*ZFF(IR,1)*R2DRDI(IR,IM)
            END DO
C
            CALL CRADINT(IM,CINTG,GRGR)
            CALL CRADINT(IM,CINTF,FRFR)
C
            NORMF = SQRT(GRGR+FRFR)
            KKF = LF + ABS(KAPF(1))
C
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  DO IPOL = IPOL1,IPOL2
                     IF ( ABS(AME1(JF,JI,IPOL)).GT.1D-8 ) GOTO 20
                     IF ( ABS(AME2(JF,JI,IPOL)).GT.1D-8 ) GOTO 20
                  END DO
               END DO
            END DO
C -------------------------------------- all angular matrix elements = 0
            CYCLE
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
C
 20         CONTINUE
            CALL CINTHFF(IM,IRTOP,ZGF,ZFF,ZGI,ZFI,RMEHF,NSOLF,NSOLI,
     &                   DRDI(1,IM))
C
C -------------------------------------- calculate total matrix elements
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  DO IPOL = IPOL1,IPOL2
C
                     ME(IKMTF,IKMTI,IPOL) = ME(IKMTF,IKMTI,IPOL)
     &                  + RMEHF(KF,KI)*AME1(JF,JI,IPOL)
C
                  END DO
               END DO
            END DO
C
            RME(IKMTF,IKMTI) = DREAL(RMEHF(1,1)/(NORMF*NORMI))
C
            IF ( KAPF(1).EQ.-1 .AND. KAPI(1).EQ.-1 ) THEN
               BT1(IKMTF,IKMTI) = RME(IKMTF,IKMTI)*(-2D0/3D0)*B_AU2CGS
            ELSE
               BT1(IKMTF,IKMTI) = RME(IKMTF,IKMTI)
     &                            *(2D0/DBLE(KAPF(1)+KAPI(1)+2))
     &                            *B_AU2CGS
            END IF
C
            IF ( NINT(MJI-MJIREP).EQ.0 ) THEN
               IBTF(KKF) = IKMTF
               IBTI(KKI) = IKMTI
            END IF
C
            IF ( NINT(MJI-MJF).NE.+1 ) THEN
               WRITE (6,*) ' WARNING ************************'
               WRITE (6,*) ' MJ-selection rules ??????????????'
               WRITE (6,*) ' MJF=',MJF,'  MJI=',MJI
            END IF
C
         END DO
C
      END DO
C
      END
