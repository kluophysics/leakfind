C*==nucsscpl.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE NUCSSCPL(IECURR,MEZZ,TAUQ,MSSQ,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the nuclear spin-spin coupling tensor  A_ij           *
C   *                                                                  *
C   *   KB        1.3807D-16        erg/K                              *
C   *   HBAR      1.05459D-27       erg*s                              *
C   *   RY        13.6058D0         erg                                *
C   *   MN        5.05082D-24       erg/G                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NCL,IQ_MBCL,NSYM,SYMACCEPTED,MROTR
      USE MOD_TAUIJ,ONLY:TAUIJ,DQCLU_TAUIJ_CL,NTAUIJ,ITAUIJ_LOOP,
     &    ITAUJI_LOOP,IQ_TAUIJ_CL,NQCLU_TAUIJ_CL,N5VEC_TAUIJ_CL,
     &    RQCLU_TAUIJ_CL,SELECTED_TAUIJ_CL
      USE MOD_SITES,ONLY:NQMAX,IQAT,ITOQ,NOQ,ICPA,NQ,QBAS
      USE MOD_ANGMOM,ONLY:AME_G,NMEMAX,NKMMAX,NKMQ,NKM,NL,IHFF
      USE MOD_RMESH,ONLY:NRMAX,JRWS,FINITE_NUCLEUS,DRDI
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_TYPES,ONLY:Z,NTMAX,CONC,IMT,NT,LTXT_T,TXT_T,NCPLWFMAX
      USE MOD_CALCMODE,ONLY:IREL,NONMAG,KMROT
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET,IFILCBWF,
     &    IPRINT
      USE MOD_CONSTANTS,ONLY:B_AU2CGS,HBAR_CGS,KB_CGS,PI,C0,RY_ERG,
     &    MN_CGS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NUCSSCPL')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      REAL*8 A11,AIJ,AITJT(:,:,:,:,:),CNUC(NTMAX),DITJT(:,:,:),
     &       DVECITJT(:,:,:,:),FEV,FHZ,FTP,INUC(NTMAX),MAT1(3,3),
     &       MAT2(3,3),MNUC(NTMAX),QNUC(NTMAX),R11,RAT,RIJ,
     &       RITJT(:,:,:,:,:),RVECP(3),XA,XB,XMAX,XMIN,XVEC(:),YMAX,
     &       YMIN,YVEC(:)
      CHARACTER*10 AIJ1(3,3),AIJ2(3,3),FMTH,FMTK,RIJ1(3,3),RIJ2(3,3),
     &             ZEROH,ZEROK
      COMPLEX*16 AIJINT(:,:,:,:,:),CSUM,DMAMC(:,:),DMATT(:,:,:),
     &           DTILT(:,:,:),ME(:,:,:,:),MTIJ(:,:,:),MTJI(:,:,:),
     &           RIJINT(:,:,:,:,:),RMEHF(NCPLWFMAX,NCPLWFMAX),W1(:,:),
     &           W2(:,:),WEKP,ZF(:,:,:),ZG(:,:,:)
      INTEGER ANUC(NTMAX),I,I2,IA_ERR,ICL,IJM,IKM1,IKM2,IKMCB(2,NKMMAX),
     &        IKMT1,IKMT2,ILJ,ILOOP,IM,IO,IPOL,IQ,IRTOP,ISYM,IT,ITAUIJ,
     &        ITAUJI,J,J1,J2,JJM,JLJ,JO,JPOL,JQ,JQCLU,JT,K1,K2,LFILNAM,
     &        LFN,M,N,NIJ,NITJT(:,:),NPOL,NSAMER,NSOLCB(NKMMAX),
     &        NTROUBLE,NX
      LOGICAL CHECK,SAMEA,SAMER
      CHARACTER*80 FILNAM
      CHARACTER*20 LEG(NT)
      LOGICAL RVEC_SAME
      SAVE AIJINT,RIJINT
C
C*** End of declarations rewritten by SPAG
C
      DATA ZEROH/'    0     '/,FMTH/'(F10.4) '/,ZEROK/'    0     '/,
     &     FMTK/'(F9.3,1X) '/
C
      ALLOCATABLE ZG,ZF,ME,MTIJ,MTJI,W1,W2,AIJINT,RIJINT
      ALLOCATABLE DMAMC,DTILT,DMATT
      ALLOCATABLE AITJT,RITJT,DITJT,NITJT,DVECITJT,XVEC,YVEC
C
      IF ( .NOT.NONMAG ) CALL STOP_MESSAGE(ROUTINE,'ferromagn. system')
      IF ( IREL.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'IREL != 3 ')
      IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KMROT != 0 ')
      IF ( NEPATH.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NEPATH != 1 ')
      IF ( FINITE_NUCLEUS )
     &      CALL STOP_MESSAGE(ROUTINE,'FINITE_NUCLEUS = .TRUE.')
      IF ( IGRID(1).NE.5 ) CALL STOP_MESSAGE(ROUTINE,'IGRID != 5 ')
C
      M = NKMMAX
C
      ALLOCATE (ME(M,M,3,NT))
      ALLOCATE (W1(M,M),W2(M,M),MTIJ(M,M,3),MTJI(M,M,3))
      ALLOCATE (DMAMC(M,M),DMATT(M,M,NT),DTILT(M,M,NT))
      ALLOCATE (ZG(NRMAX,2,M),ZF(NRMAX,2,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZF')
C
      WEKP = -WETAB(IECURR,1)/PI
C
C=======================================================================
C                     initialize all if  IECURR = 1
C=======================================================================
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (AIJINT(3,3,NTAUIJ,NT,NT))
         ALLOCATE (RIJINT(3,3,NTAUIJ,NT,NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: AIJINT')
C
         CALL CINIT(3*3*NTAUIJ*NT*NT,AIJINT)
         CALL CINIT(3*3*NTAUIJ*NT*NT,RIJINT)
C
      END IF
C=======================================================================
C
      CHECK = .TRUE.
C     CHECK = .FALSE.
      NPOL = 3
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         M = NKMMAX
         N = NKMQ(IQAT(1,IT))
         IM = IMT(IT)
         IM = 1
         IRTOP = JRWS(IM)
C
C--------------------------- calculate projection matrices DMAT and DTIL
C
         IQ = IQAT(1,IT)
C
         IF ( ICPA(IQ).NE.0 ) CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),
     &        DTILT(1,1,IT),DMAMC,N,MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
C=======================================================================
C                      calculate matrix elements
C=======================================================================
C
         CALL CINIT(M*M*3*NT,ME)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,0,ZG,ZF,ZG,ZF,IRTOP,NSOLCB,
     &                        IKMCB)
C
         DO IKMT2 = 1,N
C
            DO IKMT1 = 1,N
C
               DO K2 = 1,NSOLCB(IKMT2)
                  J2 = IKMCB(K2,IKMT2)
                  DO K1 = 1,NSOLCB(IKMT1)
                     J1 = IKMCB(K1,IKMT1)
                     DO IPOL = 1,NPOL
                        IF ( ABS(AME_G(J1,J2,IPOL,IHFF)).GT.1D-8 )
     &                       GOTO 10
                     END DO
                  END DO
               END DO
C -------------------------------------- all angular matrix elements = 0
               CYCLE
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
C
C
 10            CONTINUE
               CALL CINTHFF(IM,IRTOP,ZG(1,1,IKMT1),ZF(1,1,IKMT1),
     &                      ZG(1,1,IKMT2),ZF(1,1,IKMT2),RMEHF,
     &                      NSOLCB(IKMT1),NSOLCB(IKMT2),DRDI(1,IM))
C
C -------------------------------------- calculate total matrix elements
C
               DO K2 = 1,NSOLCB(IKMT2)
                  IKM2 = IKMCB(K2,IKMT2)
C
                  DO K1 = 1,NSOLCB(IKMT1)
                     IKM1 = IKMCB(K1,IKMT1)
C
                     DO IPOL = 1,NPOL
                        ME(IKMT1,IKMT2,IPOL,IT)
     &                     = ME(IKMT1,IKMT2,IPOL,IT)
     &                     + B_AU2CGS*AME_G(IKM1,IKM2,IPOL,IHFF)
     &                     *RMEHF(K1,K2)
                     END DO
C
                  END DO
C
               END DO
C
            END DO
C
         END DO
C
         IF ( CHECK ) THEN
            J = 0
            DO JLJ = 1,NL
               DO JJM = 1,4*JLJ - 2
                  J = J + 1
                  I = 0
                  DO ILJ = 1,NL
                     DO IJM = 1,4*ILJ - 2
                        I = I + 1
                        XA = ABS(MEZZ(I,J,IT,IHFF))
                        XB = ABS(ME(I,J,2,IT))
                        IF ( (XA+XB).GT.1D-8 ) THEN
                           RAT = ABS(1D0-MEZZ(I,J,IT,IHFF)/ME(I,J,2,IT))
                           IF ( RAT.GT.1D-10 .AND. 
     &                          (ABS(JLJ-ILJ).NE.2 .OR. XA.GT.1D-8) )
     &                          WRITE (6,99002) I,J,MEZZ(I,J,IT,IHFF),
     &                                 ME(I,J,2,IT),RAT
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
C
C----------------------------- from now on:  CARTESIAN basis used for ME
C
         CALL MECHNGBAS('SPH>CAR',ME(1,1,1,IT),NKM,NKMMAX)
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
C
      N = NKM
      M = NKMMAX
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               ITAUJI = ITAUJI_LOOP(ILOOP)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
C
                  DO JO = 1,NOQ(JQ)
                     JT = ITOQ(JO,JQ)
C
                     IF ( ICPA(JQ).NE.0 ) THEN
                        CALL CMATMUL(N,M,TAUIJ(1,1,1,ITAUIJ),
     &                               DTILT(1,1,JT),W1)
                     ELSE
                        CALL CMATCOP(N,M,TAUIJ(1,1,1,ITAUIJ),W1)
                     END IF
C
                     IF ( ICPA(IQ).NE.0 ) THEN
                        CALL CMATMUL(N,M,DMATT(1,1,IT),W1,W2)
                        DO IPOL = 1,3
                           CALL CMATMUL(N,M,ME(1,1,IPOL,IT),W2,
     &                                  MTIJ(1,1,IPOL))
                        END DO
                     ELSE
                        DO IPOL = 1,3
                           CALL CMATMUL(N,M,ME(1,1,IPOL,IT),W1,
     &                                  MTIJ(1,1,IPOL))
                        END DO
                     END IF
C
C
                     IF ( ICPA(IQ).NE.0 ) THEN
                        CALL CMATMUL(N,M,TAUIJ(1,1,1,ITAUJI),
     &                               DTILT(1,1,IT),W1)
                     ELSE
                        CALL CMATCOP(N,M,TAUIJ(1,1,1,ITAUJI),W1)
                     END IF
C
                     IF ( ICPA(JQ).NE.0 ) THEN
                        CALL CMATMUL(N,M,DMATT(1,1,JT),W1,W2)
                        DO IPOL = 1,3
                           CALL CMATMUL(N,M,ME(1,1,IPOL,JT),W2,
     &                                  MTJI(1,1,IPOL))
                        END DO
                     ELSE
                        DO IPOL = 1,3
                           CALL CMATMUL(N,M,ME(1,1,IPOL,JT),W1,
     &                                  MTJI(1,1,IPOL))
                        END DO
                     END IF
C
                     DO IPOL = 1,3
                        DO JPOL = 1,3
C
C--------------------------------------- FULL nuclear spin-spin coupling
C
                           CALL CMATMUL(N,M,MTIJ(1,1,IPOL),
     &                                  MTJI(1,1,JPOL),W1)
C
                           CSUM = C0
                           DO I = 1,NKM
                              CSUM = CSUM + W1(I,I)
                           END DO
C
                           AIJINT(IPOL,JPOL,ITAUIJ,IT,JT)
     &                        = AIJINT(IPOL,JPOL,ITAUIJ,IT,JT)
     &                        + WEKP*CSUM
C
C
C------------- Ruderman-Kittel nuclear spin-spin coupling; i.e. s-s-part
C
                           CALL CMATMUL(2,M,MTIJ(1,1,IPOL),
     &                                  MTJI(1,1,JPOL),W1)
C
                           CSUM = C0
                           DO I = 1,2
                              CSUM = CSUM + W1(I,I)
                           END DO
C
                           RIJINT(IPOL,JPOL,ITAUIJ,IT,JT)
     &                        = RIJINT(IPOL,JPOL,ITAUIJ,IT,JT)
     &                        + WEKP*CSUM
C
                        END DO
                     END DO
C
                  END DO
               END DO
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
C=======================================================================
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'RK.dat')
 100  CONTINUE
      READ (IOTMP,*,END=200)
      GOTO 100
 200  CONTINUE
      WRITE (IOTMP,'(2e20.7)') DIMAG(RIJINT(1,1,1,1,1))
C
      IF ( IECURR.NE.NETAB(1) ) THEN
         DEALLOCATE (ZG,ZF,ME,MTIJ,MTJI,W1,W2)
         DEALLOCATE (DMAMC,DTILT,DMATT,TAUIJ)
         DEALLOCATE (AITJT,RITJT,DITJT,NITJT,DVECITJT,XVEC,YVEC)
         RETURN
      END IF
C
C=======================================================================
C                              write results
C=======================================================================
C
      FILNAM = DATSET(1:LDATSET)//'_A_ij.dat'
      LFN = LDATSET + 9
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP+1,'tab_aij.dat')
C
      WRITE (6,99006)
      WRITE (IOTMP,99006)
C
      WRITE (IOTMP,99008)
      DO IQ = 1,NQ
         WRITE (IOTMP,99009) IQ,NOQ(IQ),
     &                       (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,NOQ(IQ)
     &                       )
      END DO
C
      NIJ = SUM(NQCLU_TAUIJ_CL(1:NCL))
C
      ALLOCATE (XVEC(3*NIJ),YVEC(3*NIJ),DVECITJT(3,NIJ,NT,NT))
      ALLOCATE (NITJT(NT,NT),DITJT(NIJ,NT,NT))
      ALLOCATE (AITJT(3,3,NIJ,NT,NT),RITJT(3,3,NIJ,NT,NT),STAT=IA_ERR)
C
      NITJT(1:NT,1:NT) = 0
C
      DO IT = 1,NT
         CALL TABMUQIN(Z(IT),ANUC(IT),MNUC(IT),QNUC(IT),INUC(IT),IPRINT)
         CNUC(IT) = MNUC(IT)*MN_CGS
         I2 = NINT(INUC(IT)*2)
         IF ( MOD(I2,2).EQ.0 ) THEN
            WRITE (6,99003) IT,TXT_T(IT),MNUC(IT),NINT(INUC(IT)),' '
         ELSE
            WRITE (6,99003) IT,TXT_T(IT),MNUC(IT),I2,'/2'
         END IF
      END DO
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               ITAUJI = ITAUJI_LOOP(ILOOP)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
C
                  DO JO = 1,NOQ(JQ)
                     JT = ITOQ(JO,JQ)
C
C-------------------------------------------------- exclude on-site term
C
                     IF ( DQCLU_TAUIJ_CL(JQCLU,ICL).GE.1D-6 ) THEN
C
                        FEV = -CNUC(IT)*CNUC(JT)/(2*RY_ERG)
                        FHZ = FEV/(HBAR_CGS*2*PI)
                        FTP = 1.0D9*FEV/KB_CGS
C
                        WRITE (6,99004) IQ,IT,JQ,JT,(QBAS(I,IQ),I=1,3),
     &                                  (QBAS(I,JQ),I=1,3)
C
                        WRITE (IOTMP,99004) IQ,IT,JQ,JT,
     &                         (QBAS(I,IQ),I=1,3),(QBAS(I,JQ),I=1,3)
C
                        NITJT(IT,JT) = NITJT(IT,JT) + 1
                        I = NITJT(IT,JT)
C
                        DO JPOL = 1,3
                           DO IPOL = 1,3
                              AITJT(IPOL,JPOL,I,IT,JT)
     &                           = DIMAG(AIJINT(IPOL,JPOL,ITAUIJ,IT,JT))
C
                              RITJT(IPOL,JPOL,I,IT,JT)
     &                           = DIMAG(RIJINT(IPOL,JPOL,ITAUIJ,IT,JT))
                           END DO
                        END DO
C
                        A11 = AITJT(1,1,I,IT,JT)
                        R11 = RITJT(1,1,I,IT,JT)
                        DO JPOL = 1,3
                           DO IPOL = 1,3
C
                              AIJ = AITJT(IPOL,JPOL,I,IT,JT)
                              IF ( ABS(AIJ/A11).LT.1D-10 ) THEN
                                 AITJT(IPOL,JPOL,I,IT,JT) = 0D0
                                 AIJ1(IPOL,JPOL) = ZEROH
                                 AIJ2(IPOL,JPOL) = ZEROK
                              ELSE
                                 WRITE (AIJ1(IPOL,JPOL),FMT=FMTH)
     &                                  FHZ*AIJ
                                 WRITE (AIJ2(IPOL,JPOL),FMT=FMTK)
     &                                  FTP*AIJ
                              END IF
C
                              RIJ = RITJT(IPOL,JPOL,I,IT,JT)
                              IF ( ABS(RIJ/R11).LT.1D-10 ) THEN
                                 RITJT(IPOL,JPOL,I,IT,JT) = 0D0
                                 RIJ1(IPOL,JPOL) = ZEROH
                                 RIJ2(IPOL,JPOL) = ZEROK
                              ELSE
                                 WRITE (RIJ1(IPOL,JPOL),FMT=FMTH)
     &                                  FHZ*RIJ
                                 WRITE (RIJ2(IPOL,JPOL),FMT=FMTK)
     &                                  FTP*RIJ
                              END IF
                           END DO
                        END DO
C
                        DVECITJT(1:3,I,IT,JT)
     &                     = RQCLU_TAUIJ_CL(1:3,JQCLU,ICL)
                        DITJT(I,IT,JT) = DQCLU_TAUIJ_CL(JQCLU,ICL)
C
                        WRITE (IOTMP+1,'(f10.5,2x,2f15.5)')
     &                         2*DITJT(I,IT,JT),FTP*A11,FTP*R11
C
C-----------------------------------------------------------------------
C
                        WRITE (6,99005) ITAUIJ,ITAUJI,
     &                                  N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                                  RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                                  DQCLU_TAUIJ_CL(JQCLU,ICL),
     &                                  ((RIJ1(IPOL,JPOL),JPOL=1,3),
     &                                  (RIJ2(IPOL,JPOL),JPOL=1,3),
     &                                  IPOL=1,3),
     &                                  ((AIJ1(IPOL,JPOL),JPOL=1,3),
     &                                  (AIJ2(IPOL,JPOL),JPOL=1,3),
     &                                  IPOL=1,3)
C
                        WRITE (IOTMP,99005) ITAUIJ,ITAUJI,
     &                         N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                         RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                         DQCLU_TAUIJ_CL(JQCLU,ICL),
     &                         ((RIJ1(IPOL,JPOL),JPOL=1,3),
     &                         (RIJ2(IPOL,JPOL),JPOL=1,3),IPOL=1,3),
     &                         ((AIJ1(IPOL,JPOL),JPOL=1,3),
     &                         (AIJ2(IPOL,JPOL),JPOL=1,3),IPOL=1,3)
                     END IF
C
C-----------------------------------------------------------------------
                  END DO
               END DO
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
      WRITE (6,99007) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C                      find suitable plot windows
C-----------------------------------------------------------------------
C
      XMIN = 1D10
      XMAX = 0D0
      DO IT = 1,NT
         LEG(IT) = TXT_T(IT)
         DO JT = 1,NT
            DO I = 1,NITJT(IT,JT)
               XMIN = MIN(XMIN,DITJT(I,IT,JT))
               XMAX = MAX(XMAX,DITJT(I,IT,JT))
            END DO
         END DO
      END DO
C
C-------------------------------------------- start plot at X = R_ij = 0
      XMIN = 0D0
C
      DO IT = 1,NT
C
         YMIN = 0D0
         YMAX = 0D0
C
         DO JT = 1,NT
            DO I = 1,NITJT(IT,JT)
               DO JPOL = 1,3
                  DO IPOL = 1,3
                     YMIN = MIN(YMIN,AITJT(IPOL,JPOL,I,IT,JT))
                     YMAX = MAX(YMAX,AITJT(IPOL,JPOL,I,IT,JT))
                     YMIN = MIN(YMIN,RITJT(IPOL,JPOL,I,IT,JT))
                     YMAX = MAX(YMAX,RITJT(IPOL,JPOL,I,IT,JT))
                  END DO
               END DO
            END DO
         END DO
C
         YMIN = YMAX*FTP
         YMAX = YMIN*FTP
C
         CALL XMGRHEAD(DATSET,LDATSET,'Aij',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,'distance R!sij!N (a.u.)',
     &                 23,'A!sij!N (nK) ',13,' ',0,
     &                 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &                 25+LSYSTEM,
     &                 'nuclear spin-spin coupling parameters'//
     &                 ' A!sij!N for '//TXT_T(IT)(1:LTXT_T(IT))
     &                 //' at center',(50+LTXT_T(IT)+10),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,1,NT,0,LEG,LEG)
C
         CALL XMGRPOINTS(IOTMP,1,NT,0,2,1,1)
C
         DO JT = 1,NT
            NX = 0
            DO I = 1,NITJT(IT,JT)
               DO IPOL = 1,3
                  NX = NX + 1
                  XVEC(NX) = DITJT(I,IT,JT)
                  YVEC(NX) = AITJT(IPOL,IPOL,I,IT,JT)*FTP
               END DO
            END DO
C
            CALL XMGRTABLE(0,(IT-1),XVEC,YVEC,1.0D0,NX,IOTMP)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '    results written to xmgrace file ',
     &               FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C
         CALL XMGRHEAD(DATSET,LDATSET,'Rij',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,'distance R!sij!N (a.u.)',
     &                 23,'R!sij!N (nK) ',13,' ',0,
     &                 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &                 25+LSYSTEM,
     &                 'Ruderman-Kittel coupling parameters'//
     &                 ' R!sij!N for '//TXT_T(IT)(1:LTXT_T(IT))
     &                 //' at center',(48+LTXT_T(IT)+10),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,1,NT,0,LEG,LEG)
C
         CALL XMGRPOINTS(IOTMP,1,NT,0,2,1,1)
C
         DO JT = 1,NT
            NX = 0
            DO I = 1,NITJT(IT,JT)
               DO IPOL = 1,3
                  NX = NX + 1
                  XVEC(NX) = DITJT(I,IT,JT)
                  YVEC(NX) = RITJT(IPOL,IPOL,I,IT,JT)*FTP
               END DO
            END DO
C
            CALL XMGRTABLE(0,(IT-1),XVEC,YVEC,1.0D0,NX,IOTMP)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '    results written to xmgrace file ',
     &               FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IOTMP)
C
      END DO
C
C=======================================================================
C                          perform symmetry check
C=======================================================================
      WRITE (6,*) ' '
      WRITE (6,*) '    performing symmetry check '
      WRITE (6,*) ' '
C
      NTROUBLE = 0
      NSAMER = 0
C
      DO IT = 1,NT
         DO JT = 1,NT
            DO I = 1,NITJT(IT,JT)
               DO J = I + 1,NITJT(IT,JT)
                  IF ( ABS(DITJT(I,IT,JT)-DITJT(J,IT,JT)).LT.1D-8 ) THEN
                     DO ISYM = 1,NSYM
                        IF ( SYMACCEPTED(ISYM) ) THEN
C
                           CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,
     &                                DVECITJT(1,I,IT,JT),1,0D0,RVECP,1)
C
                           SAMER = RVEC_SAME(3,RVECP,DVECITJT(1,J,IT,JT)
     &                             ,1D-6)
C
                           IF ( SAMER ) THEN
                              NSAMER = NSAMER + 1
C
                              CALL DGEMM('N','N',3,3,3,1D0,
     &                           AITJT(1,1,J,IT,JT),3,MROTR(1,1,ISYM),3,
     &                           0D0,MAT1,3)
C
                              CALL DGEMM('T','N',3,3,3,1D0,
     &                           MROTR(1,1,ISYM),3,MAT1,3,0D0,MAT2,3)
C
                              SAMEA = RVEC_SAME(9,AITJT(1,1,I,IT,JT),
     &                                MAT2,ABS(MAT2(1,1))*1D-5)
C
                              IF ( .NOT.SAMEA ) THEN
                                 WRITE (6,99001) IT,JT,I,J,ISYM
                                 NTROUBLE = NTROUBLE + 1
                              END IF
C
                           END IF
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
C
      WRITE (6,*) ' '
      WRITE (6,*) '    number of coupled Rij''s  NSAMER   = ',NSAMER
      WRITE (6,*) '    number of deviations     NTROUBLE = ',NTROUBLE
      WRITE (6,*) ' '
C
C-----------------------------------------------------------------------
99001 FORMAT ('#############  WARNING for IT,JT:',2I3,'  I,J:',2I3,
     &        '   ISYM:',I3)
99002 FORMAT (5X,'checking hyperfine matrix element (z)-pol',/,2I3,
     &        2E17.8,2X,/,6X,2E17.8,2X,'relat. deviat.',E12.3,/)
99003 FORMAT (/,10X,'nuclear data for atom type',I3,3X,A,/,10X,
     &        'magnetic moment   m =',F8.5,' m_N',/,10X,
     &        'spin              I =',I3,A,/)
99004 FORMAT (/,5X,'IQ =',I3,' IT =',I3,20X,' JQ =',I3,' JT =',I3,/,5X,
     &        2(3X,'->Q = (',F7.3,',',F7.3,',',F7.3,')',5X))
99005 FORMAT (//,5X,'IJ =',I3,2X,'JI =',I3,2X,'->N =',3I3,2X,'->DR =',
     &        3F7.3,2X,'D =',F7.3,/,/19X,'A_ij [Hz] ',25X,'A_ij [nK]',
     &        //,5X,3X,3A,5X,3A,/,5X,'RK ',3A,5X,3A,/,5X,'   ',3A,5X,3A,
     &        :,//,5X,'   ',3A,5X,3A,/,5X,'tot',3A,5X,3A,/,5X,'   ',3A,
     &        5X,3A)
99006 FORMAT (/,1X,79('*'),/,35X,'<NUCSSCPL>',/,28X,
     &        'XC-coupling constants A_ij',/,1X,79('*'),/)
99007 FORMAT (/,5X,'results written to file: ',A,/)
99008 FORMAT (/,10X,'number of sites   NQ = ',I3,/,10X,
     &        'number of types   NT = ',I3,/,10X,'site occupation:')
99009 FORMAT (10X,2I4,10(I3,F6.3))
      END
