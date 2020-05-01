C*==gfunmat.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE GFUNMAT(IFILB,IFILA,IFILA_LHS,IPRINT,ERYD,MEZZ,TAUT,
     &                   ZGA1,ZGA2,ZFA1,ZFA2,JGA1,JGA2,JFA1,JFA2,MC,IE,
     &                   GFMAT)
C   ********************************************************************
C   *                                                                  *
C   *   subroutine to calculate the Green's function matrix            *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKMQ,NMEMAX,NKMMAX,NKM,NLM
      USE MOD_SITES,ONLY:IQAT,NQMAX,NQ,DROTQ
      USE MOD_CALCMODE,ONLY:IREL,KMROT,LHS_SOL_EQ_RHS_SOL
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS,R2DRDI,NSF,NSFMAX,JRCUT,
     &    JRCRI,NMMAX
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,SYSTEM,LDATSET,IDUMMY
      USE MOD_TYPES,ONLY:TXT_T,NCPLWFMAX,NTMAX,IMT,NT,LTXT_T,ITBOT,ITTOP
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_DMFT_LDAU,ONLY:KSELF,SYMMETRISE_OCC
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GFUNMAT')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IE,IFILA,IFILA_LHS,IFILB,IPRINT,MC
      COMPLEX*16 GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX),JFA1(NRMAX,MC),
     &           JFA2(NRMAX,MC),JGA1(NRMAX,MC),JGA2(NRMAX,MC),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),ZFA1(NRMAX,MC),ZFA2(NRMAX,MC)
     &           ,ZGA1(NRMAX,MC),ZGA2(NRMAX,MC)
C
C Local variables
C
      LOGICAL AME_NE_0(:,:,:)
      COMPLEX*16 AUX,CINT(:),FB1,FB2,GB1,GB2,IZ1J(:),IZ1Z(:),JZ2,
     &           MEAB2(:,:),MEB1A(:,:),MEIRR(:,:),NORM,R1MCLEB(:,:,:),
     &           T1(:,:,:),T2(:,:,:),W1(:,:),W2(:,:),W3(:,:),X1,X2,X3,
     &           Z1F(:,:,:),Z1G(:,:,:),Z1J(:),Z1Z(:),Z2F(:,:,:),
     &           Z2G(:,:,:),Z2Z(:),ZZ2
      REAL*8 D12,D13,MAXDEL,RELGNT,RSUM,RWGT,WINTR(:,:)
      CHARACTER*80 FILNAM
      INTEGER I,I1,I2,I3,I4,IA_ERR,IFLAG1,IFLAG2,IFLAG_CHECK,IKM,IKMA(:)
     &        ,IKMCPLWFB(:,:),IKMEND,IKMME(:,:),IKMSTART,IKMTA,IKMTB1,
     &        IKMTB2,IKM_KA1,IKM_KA2,IKM_KB1,IKM_KB2,IM,IQ,IRBOT,IRTOP,
     &        IRTOPTMP,ISF,ISPIN,IT,ITA,ITB,IWI,IWR,J,K,KA1,KA2,KB,KB1,
     &        KB2,LFN,LTXTGF,M,N,NKAME,NKMAME,NME(NKMMAX),NSFTMP(:),
     &        NSOLA,NSOLIKMB(:),NSPIN
      CHARACTER*3 STR3
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NSOLIKMB,IKMCPLWFB,IKMA
C
      ALLOCATABLE Z1Z,Z1J,Z2Z,IZ1Z,IZ1J
      ALLOCATABLE Z1F,Z1G,Z2G,Z2F,CINT
      ALLOCATABLE MEIRR,MEAB2,MEB1A,T1,T2
      ALLOCATABLE R1MCLEB,WINTR,IKMME,AME_NE_0,NSFTMP
      ALLOCATABLE W1,W2,W3
C
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX),W3(NKMMAX,NKMMAX))
      ALLOCATE (CINT(NRMAX))
      ALLOCATE (Z1Z(NRMAX),Z1J(NRMAX),Z2Z(NRMAX))
      ALLOCATE (IZ1Z(NRMAX),IZ1J(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: Z1J')
      ALLOCATE (T1(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (T2(NKMMAX,NKMMAX,NQMAX))
C
      N = NKM
      ALLOCATE (Z1F(NRMAX,MC,N))
      ALLOCATE (Z1G(NRMAX,MC,N))
      ALLOCATE (Z2G(NRMAX,MC,N))
      ALLOCATE (Z2F(NRMAX,MC,N),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: Z2FB')
C
      ALLOCATE (MEIRR(NKMMAX,NKMMAX))
      ALLOCATE (MEAB2(NKMMAX,NKMMAX),MEB1A(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEIRR')
C
      ALLOCATE (NSOLIKMB(NKM),IKMCPLWFB(MC,NKM),IKMA(MC))
C
      IF ( IREL.LT.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL < 2')
      IF ( IREL.EQ.2 .AND. .NOT.FULLPOT )
     &      CALL STOP_MESSAGE(ROUTINE,'ONLY FULLPOT mode for IREL = 2 ')
C
      ALLOCATE (WINTR(NRMAX,NSFMAX))
      ALLOCATE (NSFTMP(NMMAX))
      NSFTMP(1:NMMAX) = 0
C
C------------------------------------------- NKAME should be > NK=2*NL-1
      IF ( IREL.GE.3 ) THEN
C
         NKAME = 2*NINT(SQRT(DBLE(NKM/2))) + 1
         NKMAME = 2*((NKAME+1)/2)**2
C
      ELSE
C
         NKMAME = NKM
C
      END IF
C
      IF ( NKM.GT.NKMAME ) CALL STOP_MESSAGE(ROUTINE,'NKM > NKMAME')
      ALLOCATE (IKMME(NKMMAX,NKMMAX),AME_NE_0(NKMAME,NKMAME,NSFMAX))
      ALLOCATE (R1MCLEB(NKMAME,NKMAME,NSFMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CPLMAT')
C
      IF ( CHECK .AND. .NOT.LHS_SOL_EQ_RHS_SOL ) WRITE (6,*)
     &      '<gfunmat>CHECK NOT YET FOR LHS'
C
C-----------------------------------------------------------------------
C
      IF ( IPRINT.GT.0 ) THEN
         LTXTGF = 12
      ELSE
         LTXTGF = 0
      END IF
      I1 = 5
      I2 = 9
      I3 = NKM/2 + I1
      I4 = NKM/2 + I2
C
      IFLAG_CHECK = 0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).GT.0 ) THEN
C
            IQ = IQAT(1,IT)
C
C=======================================================================
C                      calculate matrix elements
C=======================================================================
C
            IM = IMT(IT)
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
C
            CALL CINIT(NKMMAX*NKMMAX,MEIRR)
            CALL CINIT(NKMMAX*NKMMAX,MEAB2)
            CALL CINIT(NKMMAX*NKMMAX,MEB1A)
C
            IWR = 200 + IT
            IWI = 300 + IT
C
            IF ( IE.EQ.1 .AND. .NOT.CHECK ) THEN
               IF ( IPRINT.GT.0 ) THEN
                  WRITE (6,99003) IT,TXT_T(IT)
                  CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,
     &                          '_GR.dat',7,FILNAM,LFN,2,IWR,
     &                          'Re G  file',10,NTMAX)
                  CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,
     &                          '_GI.dat',7,FILNAM,LFN,2,IWI,
     &                          'Im G  file',10,NTMAX)
C
                  WRITE (IWR,99004) 'REAL',SYSTEM(1:LSYSTEM),IT,
     &                              TXT_T(IT)
                  WRITE (IWI,99004) 'IMAGINARY',SYSTEM(1:LSYSTEM),IT,
     &                              TXT_T(IT)
               END IF
            END IF
C
            M = NKMMAX
            N = NKMQ(IQ)
C
C-----------------------------------------------------------------------
C                    read basis set wave functions
C
C    NOTE: these have been calculated for    B = 0
C          thus the overlap integral matrix MEZZ has to be diagonal
C-----------------------------------------------------------------------
C
            Z2G = C0
            Z1G = C0
C
            Z2F = C0
            Z1F = C0
C
            DO IKM = 1,NKMQ(IQ)
C
C
               IF ( IREL.EQ.3 ) THEN
                  READ (IFILB,REC=IKM+(IT-1)*NKM) ITB,STR3,IDUMMY,
     &                  IRTOPTMP,NSOLIKMB(IKM),
     &                  (IKMCPLWFB(J,IKM),J=1,NSOLIKMB(IKM)),
     &                  ((Z2G(I,K,IKM),I=1,IRTOPTMP),K=1,NSOLIKMB(IKM)),
     &                  ((Z2F(I,K,IKM),I=1,IRTOPTMP),K=1,NSOLIKMB(IKM))
               ELSE
                  READ (IFILB,REC=IKM+(IT-1)*NKM) ITB,STR3,IDUMMY,
     &                  NSOLIKMB(IKM),
     &                  (IKMCPLWFB(J,IKM),J=1,NSOLIKMB(IKM)),
     &                  ((Z2G(I,K,IKM),I=1,IRTOP),K=1,NCPLWFMAX)
                  Z2F = C0
C
               END IF
C
               IF ( IT.NE.ITB .OR. STR3.NE.'REG' )
     &              CALL STOP_MESSAGE(ROUTINE,'WF(FIN) inconsistent')
C
               IF ( .NOT.CHECK ) THEN
C
                  NORM = 1D0/SQRT(MEZZ(IKM,IKM,IT,1))
C
                  DO KB = 1,NSOLIKMB(IKM)
                     DO I = 1,IRTOP
                        Z2G(I,KB,IKM) = Z2G(I,KB,IKM)*NORM
                        Z1G(I,KB,IKM) = DCONJG(Z2G(I,KB,IKM))
C
                        Z2F(I,KB,IKM) = Z2F(I,KB,IKM)*NORM
                        Z1F(I,KB,IKM) = DCONJG(Z2F(I,KB,IKM))
                     END DO
                  END DO
C
               ELSE
                  KB = NSOLIKMB(IKM)
                  Z1G(1:IRTOP,1:KB,IKM) = Z2G(1:IRTOP,1:KB,IKM)
                  Z1F(1:IRTOP,1:KB,IKM) = Z2F(1:IRTOP,1:KB,IKM)
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
C             initialize parameters for integration
C-----------------------------------------------------------------------
C
            IF ( .NOT.(FULLPOT) ) THEN
C
C --------------- only L- (LAMBDA-) diagonal terms occur for IREL= 2 (3)
C
               DO J = 1,IRTOP
                  WINTR(J,1) = R2DRDI(J,IM)
               END DO
C
               CALL CINIT(NKMAME*NKMAME,R1MCLEB)
C
               DO I = 1,NKMAME
                  R1MCLEB(I,I,1) = 1D0
               END DO
C
               NME(1:NKM) = 1
               DO IKM = 1,NKM
                  IKMME(1,IKM) = IKM
               END DO
C
C--------------------------------------------------------------- FULLPOT
            ELSE IF ( IREL.LE.2 ) THEN
C
               CALL GFUNMAT_INIT_SRA(IM,NKMAME,NSOLIKMB,IKMCPLWFB,NME,
     &                               IKMME,R1MCLEB,WINTR)
C
            ELSE
C
               CALL GFUNMAT_INIT_REL(IM,NKMAME,NSOLIKMB,IKMCPLWFB,NME,
     &                               IKMME,R1MCLEB,WINTR)
C
            END IF
C
            IF ( FULLPOT ) THEN
               NSFTMP(1:NMMAX) = NSF(1:NMMAX)
               DO ISF = 1,NSFTMP(IM)
                  DO J = 1,NKMAME
                     DO I = 1,NKMAME
                        AME_NE_0(I,J,ISF) = ABS(R1MCLEB(I,J,ISF))
     &                     .GT.1D-6
                     END DO
                  END DO
               END DO
            ELSE
               NSFTMP(1:NMMAX) = 1
               DO ISF = 1,NSFTMP(IM)
                  DO J = 1,NKMAME
                     DO I = 1,NKMAME
                        AME_NE_0(I,J,ISF) = ABS(R1MCLEB(I,J,ISF))
     &                     .GT.1D-6
                     END DO
                  END DO
               END DO
C
            END IF
C
C-----------------------------------------------------------------------
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C---------- loop over Green's function states (Lambda')
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
            IF ( IREL.EQ.2 ) THEN
               NSPIN = 2
            ELSE
               NSPIN = 1
            END IF
C
            DO ISPIN = 1,NSPIN
               IKMSTART = NLM*(ISPIN-1) + 1
               IF ( IREL.EQ.2 ) THEN
                  IKMEND = IKMSTART + NLM - 1
               ELSE
                  IKMEND = NKMQ(IQ)
               END IF
C
               DO IKMTA = IKMSTART,IKMEND
C
                  IF ( IREL.EQ.3 ) THEN
                     READ (IFILA,REC=IKMTA+(IT-1)*NKM) ITA,STR3,IDUMMY,
     &                     IRTOPTMP,NSOLA,(IKMCPLWFB(J,IKMTA),J=1,NSOLA)
     &                     ,((ZGA2(I,K),I=1,IRTOPTMP),K=1,NSOLA),
     &                     ((ZFA2(I,K),I=1,IRTOPTMP),K=1,NSOLA)
                     IF ( FULLPOT ) THEN
                        IF ( .NOT.LHS_SOL_EQ_RHS_SOL )
     &                       READ (IFILA_LHS,REC=IKMTA+(IT-1)*NKM) ITA,
     &                       STR3,IDUMMY,IRTOPTMP,NSOLA,
     &                       (IKMCPLWFB(J,IKMTA),J=1,NSOLA),
     &                       ((ZGA1(I,K),I=1,IRTOPTMP),K=1,NSOLA),
     &                       ((ZFA1(I,K),I=1,IRTOPTMP),K=1,NSOLA)
                     END IF
                  ELSE
                     READ (IFILA,REC=IKMTA+(IT-1)*NKM) ITA,STR3,IDUMMY,
     &                     NSOLA,
     &                     (IKMCPLWFB(J,IKMTA),J=1,NSOLA),
     &                     ((ZGA2(I,K),I=1,IRTOP),K=1,NCPLWFMAX)
                     ZFA2 = C0
                  END IF
C
                  DO KA1 = 1,NCPLWFMAX
                     IKMA(KA1) = IKMCPLWFB(KA1,IKMTA)
                  END DO
C
                  IF ( NSOLIKMB(IKMTA).NE.NSOLA ) THEN
                     WRITE (6,*) 'NSOLIKMB(IKMTA)=',NSOLIKMB(IKMTA),
     &                           '  <>  NSOLA =',NSOLA
                     CALL STOP_MESSAGE(ROUTINE,
     &                                 'WF coupling inconsistent')
                  END IF
C
                  IF ( IT.NE.ITA .OR. STR3.NE.'REG' )
     &                 CALL STOP_MESSAGE(ROUTINE,'WF(INI) inconsistent')
C
                  IF ( IREL.EQ.3 ) THEN
                     READ (IFILA,REC=IKMTA+(IT-1+NT)*NKM) ITA,STR3,
     &                     IDUMMY,IRTOPTMP,
     &                     ((JGA2(I,K),I=1,IRTOPTMP),K=1,NSOLA),
     &                     ((JFA2(I,K),I=1,IRTOPTMP),K=1,NSOLA)
                     IF ( .NOT.LHS_SOL_EQ_RHS_SOL )
     &                    READ (IFILA_LHS,REC=IKMTA+(IT-1+NT)*NKM) ITA,
     &                    STR3,IDUMMY,IRTOPTMP,
     &                    ((JGA1(I,K),I=1,IRTOPTMP),K=1,NSOLA),
     &                    ((JFA1(I,K),I=1,IRTOPTMP),K=1,NSOLA)
                  ELSE
                     READ (IFILA,REC=IKMTA+(IT-1+NT)*NKM) ITA,STR3,
     &                     IDUMMY,((JGA2(I,K),I=1,IRTOP),K=1,NSOLA)
C
                     JFA2 = C0
                  END IF
                  IF ( IT.NE.ITA )
     &                  CALL STOP_MESSAGE(ROUTINE,'IT(INI) <> IT')
                  IF ( STR3.NE.'IRR' .AND. .NOT.CHECK )
     &                 CALL STOP_MESSAGE(ROUTINE,'WFT(INI) <> IRR')
C
                  IF ( CHECK ) THEN
                     JGA2(1:IRTOP,1:NSOLA) = ZGA2(1:IRTOP,1:NSOLA)
                     JFA2(1:IRTOP,1:NSOLA) = ZFA2(1:IRTOP,1:NSOLA)
                  END IF
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB IKMTB1
C------------------------------- loop over FIRST basis function (Lambda)
C
                  DO IKMTB1 = IKMSTART,IKMEND
C
                     DO I = 1,IRTOP
                        Z1Z(I) = C0
                        Z1J(I) = C0
                        Z2Z(I) = C0
                     END DO
C
                     IFLAG1 = 0
C
                     DO ISF = 1,NSFTMP(IM)
C
                        IF ( ISF.EQ.1 ) THEN
                           IRBOT = 1
                        ELSE
                           IRBOT = JRCUT(1,IM) + 1
                        END IF
C
                        DO KA1 = 1,NSOLA
                           IKM_KA1 = IKMA(KA1)
                           DO KB1 = 1,NSOLIKMB(IKMTB1)
                              IKM_KB1 = IKMCPLWFB(KB1,IKMTB1)
C
                              IF ( AME_NE_0(IKM_KA1,IKM_KB1,ISF) ) THEN
C
                                 IFLAG1 = 1
C
                                 RELGNT = DREAL
     &                              (R1MCLEB(IKM_KA1,IKM_KB1,ISF))
                                 DO I = IRBOT,IRTOP
                                    RWGT = RELGNT*WINTR(I,ISF)
C
                                    GB1 = RWGT*Z1G(I,KB1,IKMTB1)
                                    FB1 = RWGT*Z1F(I,KB1,IKMTB1)
C
                                    Z1Z(I) = Z1Z(I) + GB1*ZGA2(I,KA1)
     &                                 + FB1*ZFA2(I,KA1)
                                    Z1J(I) = Z1J(I) + GB1*JGA2(I,KA1)
     &                                 + FB1*JFA2(I,KA1)
C
                                    GB2 = RWGT*Z2G(I,KB1,IKMTB1)
                                    FB2 = RWGT*Z2F(I,KB1,IKMTB1)
C
                                    Z2Z(I) = Z2Z(I) + ZGA1(I,KA1)
     &                                 *GB2 + ZFA1(I,KA1)*FB2
                                 END DO
C
                              END IF
                           END DO
                        END DO
C
                     END DO
C
C################################################################ IFLAG1
C
                     IF ( IFLAG1.NE.0 ) THEN
C
                        CALL CRADINT_R(IM,Z1Z,IZ1Z)
                        CALL CRADINT_INW_R(IM,Z1J,IZ1J)
                        CALL CRADINT(IM,Z2Z,AUX)
C
                        MEB1A(IKMTB1,IKMTA) = MEB1A(IKMTB1,IKMTA)
     &                     + IZ1Z(IRTOP)
C
                        MEAB2(IKMTA,IKMTB1) = MEAB2(IKMTA,IKMTB1) + AUX
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB IKMTB2
C---------------------------- loop over SECOND basis function (Lambda'')
C
                        DO IKMTB2 = IKMSTART,IKMEND
C
                           DO I = 1,IRTOP
                              CINT(I) = C0
                           END DO
C
                           IFLAG2 = 0
C
                           DO ISF = 1,NSFTMP(IM)
C
                              IF ( ISF.EQ.1 ) THEN
                                 IRBOT = 1
                              ELSE
                                 IRBOT = JRCUT(1,IM) + 1
                              END IF
C
                              DO KA2 = 1,NSOLA
                                 IKM_KA2 = IKMA(KA2)
                                 DO KB2 = 1,NSOLIKMB(IKMTB2)
                                    IKM_KB2 = IKMCPLWFB(KB2,IKMTB2)
C
                                    IF ( AME_NE_0(IKM_KA2,IKM_KB2,ISF) )
     &                                 THEN
C
                                       IFLAG2 = 1
C
                                       RELGNT = DREAL
     &                                    (R1MCLEB(IKM_KA2,IKM_KB2,ISF))
                                       DO I = IRBOT,IRTOP
                                         RWGT = RELGNT*WINTR(I,ISF)
C
                                         GB2 = RWGT*Z2G(I,KB2,IKMTB2)
                                         FB2 = RWGT*Z2F(I,KB2,IKMTB2)
C
                                         JZ2 = JGA1(I,KA2)
     &                                      *GB2 + JFA1(I,KA2)*FB2
                                         ZZ2 = ZGA1(I,KA2)
     &                                      *GB2 + ZFA1(I,KA2)*FB2
C
                                         CINT(I) = CINT(I) + IZ1Z(I)
     &                                      *JZ2 + IZ1J(I)*ZZ2
                                       END DO
C
                                    END IF
                                 END DO
                              END DO
C
                           END DO
C
                           IF ( IFLAG2.NE.0 ) THEN
                              CALL CRADINT(IM,CINT,AUX)
C
                              MEIRR(IKMTB1,IKMTB2)
     &                           = MEIRR(IKMTB1,IKMTB2) + AUX
C
                           END IF
C
                        END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB IKMTB2
C
                     END IF
C################################################################ IFLAG1
C
                  END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB IKMTB1
C
               END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
            END DO
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C-----------------------------------------------------------------------
C                   test consitency of matrix elements
C-----------------------------------------------------------------------
C
            IF ( CHECK ) THEN
C
               WRITE (6,99001) IT,NSFTMP(IM)
C
               CALL CMATSTRUCT('MEZZ ',MEZZ(1,1,IT,1),N,M,3,3,1,1D-8,6)
C
               MAXDEL = 0.0D0
               DO I = 1,N
                  DO J = 1,N
                     X1 = MEZZ(I,J,IT,1)
                     X2 = MEB1A(I,J)
                     X3 = MEAB2(I,J)
                     IF ( ABS(X1)+ABS(X2)+ABS(X3).GT.1D-8 ) THEN
                        D12 = ABS(X1-X2)/(ABS(X1)+ABS(X2))
                        D13 = ABS(X1-X3)/(ABS(X1)+ABS(X3))
                        WRITE (6,99006) ' '
                        WRITE (6,99006) 'MEZZ  ',I,J,X1
                        WRITE (6,99006) 'MEB1A ',I,J,X2,D12
                        WRITE (6,99006) 'MEAB2 ',I,J,X3,D13
                        IF ( D12.GT.1D-12 .OR. D13.GT.1D-12 ) THEN
                           WRITE (6,99002)
                           IFLAG_CHECK = 1
                           MAXDEL = MAX(MAXDEL,D12)
                           MAXDEL = MAX(MAXDEL,D13)
                        END IF
                     END IF
                  END DO
               END DO
               IF ( IFLAG_CHECK.EQ.1 ) WRITE (6,*) 'MAXIMAL DIFFERENCE',
     &              MAXDEL
C
               CALL CMATSTRUCT('MEIRR',MEIRR,N,M,3,3,1,1D-8,6)
               MAXDEL = 0.0D0
C
               RSUM = 0D0
               DO I = 1,N
                  DO J = 1,N
                     X2 = 0D0
                     DO K = 1,N
                        X2 = X2 + MEB1A(I,K)*MEAB2(K,J)
                     END DO
                     X1 = MEIRR(I,J)
                     RSUM = RSUM + ABS(X1-X2)
                     IF ( ABS(X1).GT.1D-8 .OR. ABS(X2).GT.1D-8 ) THEN
                        WRITE (6,*)
                        WRITE (6,99006) 'IRR',I,J,X1
                        D12 = ABS(X1-X2)/(ABS(X1)+ABS(X2))
                        WRITE (6,99006) 'REG',I,J,X2,D12
                        IF ( D12.GT.1D-12 ) THEN
                           WRITE (6,99002)
                           IFLAG_CHECK = 1
                           MAXDEL = MAX(MAXDEL,D12)
                        END IF
                     END IF
                  END DO
               END DO
               IF ( IFLAG_CHECK.EQ.1 ) WRITE (6,*)
     &               'MAXIMAL DIFFERENCE IRR',MAXDEL
C
            END IF
C
C
C=======================================================================
C                 set up    Green's function matrix
C=======================================================================
C
            CALL ZGEMM('N','N',N,N,N,C1,MEB1A,M,TAUT(1,1,IT),M,C0,W1,M)
C
            CALL ZGEMM('N','N',N,N,N,C1,W1,M,MEAB2,M,C0,W2,M)
C
            DO J = 1,N
               CALL ZAXPY(N,-C1,MEIRR(1,J),1,W2(1,J),1)
            END DO
C
C
            IF ( SYMMETRISE_OCC ) THEN
C
               T1 = C0
               T2 = C0
C
               IF ( KMROT.NE.0 ) THEN
C
                  M = NKMMAX
                  IF ( IREL.NE.2 ) THEN
                     N = NKMQ(IQ)
                  ELSE
                     N = NKM
                  END IF
C
                  CALL ROTATE(W2,'L->G',W3,N,DROTQ(1,1,IQ),M)
                  T1(:,:,IQ) = W3
                  CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1.0D0,T1,T2,W1,
     &                           NQ,NKMQ,DROT,IQORGQP,SYMUNITARY,
     &                           SYMACCEPTED,NSYM,NSYMACCEPTED,NQMAX,
     &                           NKMMAX)
                  W3 = T2(:,:,IQ)
                  CALL ROTATE(W3,'G->L',W2,N,DROTQ(1,1,IQ),M)
C
               ELSE
                  T1(:,:,IQ) = W2
                  CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1.0D0,T1,T2,W1,
     &                           NQ,NKMQ,DROT,IQORGQP,SYMUNITARY,
     &                           SYMACCEPTED,NSYM,NSYMACCEPTED,NQMAX,
     &                           NKMMAX)
                  W2 = T2(:,:,IQ)
               END IF
            END IF
C
C-----------------------------------------------------------------------
C   change representation to non-relativistic / REAL spher. harmonics
C-----------------------------------------------------------------------
C
            IF ( IREL.EQ.3 ) THEN
C
               CALL CHANGEREP(NKM,NKMMAX,W2,'REL>CLM',GFMAT(1,1,IE,IT))
C
               IF ( LTXTGF.GT.0 ) CALL CMATSTRUCT('GF (l,ml,ms)',W2,N,M,
     &              IREL,IREL,0,1D-8,6)
            ELSE
C
               CALL CHANGEREP(NKM,NKMMAX,W2,'RLM>CLM',GFMAT(1,1,IE,IT))
C
               IF ( LTXTGF.GT.0 ) CALL CMATSTRUCT('GF (l,ml,ms)',W2,N,M,
     &              IREL,IREL,0,1D-8,6)
C
            END IF
C
C-----------------------------------------------------------------------
C                      print matrix if requested
C-----------------------------------------------------------------------
C
            IF ( IPRINT.GT.0 .AND. .NOT.CHECK ) THEN
               WRITE (IWR,99005) ERYD,(DREAL(GFMAT(I,I,IE,IT)),I=I1,I2),
     &                           (DREAL(GFMAT(I,I,IE,IT)),I=I3,I4)
               WRITE (IWI,99005) ERYD,(DIMAG(GFMAT(I,I,IE,IT)),I=I1,I2),
     &                           (DIMAG(GFMAT(I,I,IE,IT)),I=I3,I4)
            END IF
C
C
         END IF
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( CHECK ) THEN
         WRITE (6,*)
         IF ( IFLAG_CHECK.NE.0 ) THEN
            CALL STOP_MESSAGE(ROUTINE,
     &                        'check completed  WITHOUT success !!')
         ELSE
            CALL STOP_REGULAR(ROUTINE,'check completed successfully')
         END IF
      END IF
C
      DEALLOCATE (CINT,Z1Z,Z1J,IZ1Z,IZ1J,Z1F,Z1G,Z2G,Z2F)
C
      RETURN
99001 FORMAT (//,5X,'test consitency of matrix elements for IT  =',I3,/,
     &        5X,'with number of shape functions         NSF =',I3,/)
99002 FORMAT (80('*'))
99003 FORMAT (/,10X,'structure of Green''s function matrix for IT = ',
     &        I3,3X,A)
99004 FORMAT ('#       Green''s function matrix  ',A,' part',/,
     &        '#       system  ',A,/,'#       type    IT = ',I3,3X,A)
99005 FORMAT (66F20.10)
99006 FORMAT (2X,A,2I3,6F22.12)
      END
C*==gfunmat_init_rel.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE GFUNMAT_INIT_REL(IM,NKMAME,NCPLWF,IKMCPLWF,NME,IKMME,
     &                            R1MCLEB,WINTR)
C   ********************************************************************
C   *                                                                  *
C   *  initialize the parameters to calculate the overlap integrals    *
C   *                                                                  *
C   *                 -- RELATIVISTIC CASE --                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NSF,NRMAX,NSFMAX,DRDI,FLMSF,R,LMISF,JRCRI,JRCUT
      USE MOD_ANGMOM,ONLY:IMKM_IKM,NKM,NKMMAX,AG_RGNT,NLM_AME_RLM_EXT
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GFUNMAT_INIT_REL')
C
C Dummy arguments
C
      INTEGER IM,NKMAME
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),IKMME(NKMMAX,NKMMAX),
     &        NCPLWF(NKMMAX),NME(NKMMAX)
      COMPLEX*16 R1MCLEB(NKMAME,NKMAME,NSFMAX)
      REAL*8 WINTR(NRMAX,NSFMAX)
C
C Local variables
C
      LOGICAL CPLMAT(:,:,:)
      LOGICAL C_GT_EPS
      INTEGER IA_ERR,ICPLWF1,ICPLWF2,IKM1,IKM2,IKM3,IKM4,IRTOP,ISF,J,LM,
     &        MIKM3,MIKM4
      REAL*8 RSCAL(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CPLMAT,RSCAL
C
C------------------------------------------- NKAME should be > NK=2*NL-1
C
      ALLOCATE (CPLMAT(NKMAME,NKMAME,NSFMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CPLMAT')
C
      ALLOCATE (RSCAL(NRMAX))
C
C-----------------------------------------------------------------------
C
      NME(:) = 0
      CALL RINIT(NRMAX,RSCAL)
      CALL CINIT(NKMAME*NKMAME*NSFMAX,R1MCLEB)
C
      CPLMAT(:,:,:) = .FALSE.
C
      IRTOP = JRCRI(IM)
C
C calculate spin-angular matrix elements for shape functions
C
      DO ISF = 1,NSF(IM)
C
         LM = LMISF(ISF,IM)
C
         IF ( LM.GT.NLM_AME_RLM_EXT ) CYCLE
C
         DO IKM1 = 1,NKM
            DO IKM2 = 1,NKM
               IF ( C_GT_EPS(AG_RGNT(IKM1,IKM2,LM),1D-12) ) THEN
                  R1MCLEB(IKM1,IKM2,ISF) = AG_RGNT(IKM1,IKM2,LM)
                  CPLMAT(IKM1,IKM2,ISF) = .TRUE.
               END IF
            END DO
         END DO
C
      END DO
C
C determine matrix elements to be calculated, i.e. the nonvanishing ones
C
      DO IKM1 = 1,NKM
         DO IKM2 = 1,NKM
            DO ICPLWF1 = 1,NCPLWF(IKM1)
               IKM3 = IKMCPLWF(ICPLWF1,IKM1)
               MIKM3 = IMKM_IKM(IKM3)
C
               DO ICPLWF2 = 1,NCPLWF(IKM2)
                  IKM4 = IKMCPLWF(ICPLWF2,IKM2)
                  MIKM4 = IMKM_IKM(IKM4)
C
                  DO ISF = 1,NSF(IM)
C
                     IF ( CPLMAT(IKM3,IKM4,ISF) .OR. 
     &                    CPLMAT(MIKM3,MIKM4,ISF) ) THEN
                        NME(IKM1) = NME(IKM1) + 1
                        IKMME(NME(IKM1),IKM1) = IKM2
                        GOTO 50
                     END IF
                  END DO
               END DO
            END DO
 50      END DO
      END DO
C
C------------------------------- set up scaling factor for r-integration
C
C for 0th and 1st panel, only shape function for lm = 1 (l = 0,ml = 0)
C with flmsf(j,1,im) = SQRT_4PI
C
      DO ISF = 1,NSF(IM)
C
         LM = LMISF(ISF,IM)
C
         IF ( LM.EQ.1 ) THEN
            DO J = 1,JRCUT(1,IM)
               RSCAL(J) = R(J,IM)*R(J,IM)*DRDI(J,IM)*SQRT_4PI
            END DO
         ELSE
            DO J = 1,JRCUT(1,IM)
               RSCAL(J) = 0D0
            END DO
         END IF
C
C-------------- for remaining panels, all shape functions come into play
C
         DO J = JRCUT(1,IM) + 1,IRTOP
            RSCAL(J) = R(J,IM)*R(J,IM)*DRDI(J,IM)
     &                 *FLMSF(J-JRCUT(1,IM),ISF,IM)
         END DO
C
         WINTR(1:IRTOP,ISF) = RSCAL(1:IRTOP)
C
      END DO
C
      END
C*==gfunmat_init_sra.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE GFUNMAT_INIT_SRA(IM,NKMAME,NCPLWF,IKMCPLWF,NME,IKMME,
     &                            R1MCLEB,WINTR)
C   ********************************************************************
C   *                                                                  *
C   *  initialize the parameters to calculate the overlap integrals    *
C   *                                                                  *
C   *                -- NON-RELATIVISTIC CASE --                       *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NSF,NRMAX,NSFMAX,DRDI,FLMSF,R,LMISF,JRCRI,JRCUT
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,NLM,NKMMAX
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 RZERO
      PARAMETER (RZERO=0.0D0)
C
C Dummy arguments
C
      INTEGER IM,NKMAME
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),IKMME(NKMMAX,NKMMAX),
     &        NCPLWF(NKMMAX),NME(NKMMAX)
      COMPLEX*16 R1MCLEB(NKMAME,NKMAME,NSFMAX)
      REAL*8 WINTR(NRMAX,NSFMAX)
C
C Local variables
C
      REAL*8 GAUNT_RYLM
      REAL*8 GNT,RSCAL(NRMAX)
      INTEGER ICPLWF1,ICPLWF2,IRTOP,ISF,J,L3,L4,LM,LM1,LM2,LM3,LM4,LMSF,
     &        LSF,M3,M4,MSF
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C
      NME(:) = 0
      CALL RINIT(NRMAX,RSCAL)
      CALL CINIT(NKMAME*NKMAME*NSFMAX,R1MCLEB)
C
      IRTOP = JRCRI(IM)
C
C determine matrix elements to be calculated, i.e. the nonvanishing ones
C
      DO LM1 = 1,NLM
         NME(LM1) = 0
         DO LM2 = 1,NLM
            DO ICPLWF1 = 1,NCPLWF(LM1)
               LM3 = IKMCPLWF(ICPLWF1,LM1)
               L3 = L_LM(LM3)
               M3 = M_LM(LM3)
C
               DO ICPLWF2 = 1,NCPLWF(LM2)
                  LM4 = IKMCPLWF(ICPLWF2,LM2)
                  L4 = L_LM(LM4)
                  M4 = M_LM(LM4)
C
                  DO ISF = 1,NSF(IM)
                     LMSF = LMISF(ISF,IM)
                     LSF = L_LM(LMSF)
                     MSF = M_LM(LMSF)
C
                     GNT = GAUNT_RYLM(L3,M3,L4,M4,LSF,MSF)
C
                     IF ( ABS(GNT).GT.1D-6 ) THEN
C
                        R1MCLEB(LM3,LM4,ISF) = GNT
C
                        NME(LM1) = NME(LM1) + 1
                        IKMME(NME(LM1),LM1) = LM2
C
                        IF ( IREL.GE.1 ) THEN
C
                           R1MCLEB(NLM+LM3,NLM+LM4,ISF) = GNT
C
                           NME(NLM+LM1) = NME(NLM+LM1) + 1
                           IKMME(NME(NLM+LM1),NLM+LM1) = NLM + LM2
C
                        END IF
C
                        GOTO 50
                     END IF
                  END DO
               END DO
            END DO
 50      END DO
      END DO
C
C------------------------------- set up scaling factor for r-integration
C
C for 0th and 1st panel, only shape function for lm = 1 (l = 0,ml = 0)
C with FLMSF(j,1,im) = SQRT_4PI
C
      DO ISF = 1,NSF(IM)
C
         LM = LMISF(ISF,IM)
C
         IF ( LM.EQ.1 ) THEN
            DO J = 1,JRCUT(1,IM)
               RSCAL(J) = R(J,IM)*R(J,IM)*DRDI(J,IM)*SQRT_4PI
            END DO
         ELSE
            DO J = 1,JRCUT(1,IM)
               RSCAL(J) = RZERO
            END DO
         END IF
C
C-------------- for remaining panels, all shape functions come into play
C
         DO J = JRCUT(1,IM) + 1,IRTOP
            RSCAL(J) = R(J,IM)*R(J,IM)*DRDI(J,IM)
     &                 *FLMSF(J-JRCUT(1,IM),ISF,IM)
         END DO
C
         WINTR(1:IRTOP,ISF) = RSCAL(1:IRTOP)
C
      END DO
C
      END
