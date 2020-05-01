C*==clupotrd.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUPOTRD(STEP)
C   ********************************************************************
C   *                                                                  *
C   *  read the cluster potential information                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:A5BAS
      USE MOD_SCF,ONLY:SCFSTATUS_CLU
      USE MOD_CALCMODE,ONLY:ORBPOL,KMROT
      USE MOD_SITES,ONLY:NQCLU,QMPHI,QMTET,N5VEC_QCLU,IQ_QCLU,NQCLUMAX,
     &    DQBAS_QCLU,QBAS_QCLU,FORCE_QCLU,IM_QCLU,NO_QCLU,NQHOST,NQ_L,
     &    NQ_R,ITCLU_OQCLU,QBAS,QBAS0_QCLU
      USE MOD_RMESH,ONLY:JRCUT,RMESHTYPE,DX,NMCLU,JRNS1,JRCRI,NPAN,JRMT,
     &    JRWS,R,RMT,RWS,FULLPOT,NMHOST
      USE MOD_FILES,ONLY:IPRINT,INFO,LINFO,SYSTEM,LSYSTEM,TITLE,LTITLE,
     &    IFILPOT_CLU,POTFIL_CLU,LPOTFIL_CLU,IDUMMY
      USE MOD_TYPES,ONLY:NTCLU,NTHOST,AOPT,VT,BT,VNST,BNST,RHOCHR,
     &    RHOSPN,QEL,OBS_T,IMT,LOPT,NFPT,NLMFPT,LMIFP,KLMFP,NLMFPMAX,
     &    NTCLUMAX,X_TCLU,NVAL_TCLU,TXT_TCLU,Z_TCLU,NCOR_TCLU,
     &    NSEMCORSHL_TCLU,NCORT,NVALT,TXT_T,Z,NSEMCORSHLT,CONC,RHO2NS
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,IHFF
      IMPLICIT NONE
C*--CLUPOTRD25
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER STEP
C
C Local variables
C
      LOGICAL CHECK,TROUBLE
      CHARACTER*40 FMT03,FMT05,FMT06,FMT07,FMT10,FMT11
      CHARACTER*80 HEADER,LINE
      INTEGER I,IFP,IM,IO,IOPOT,IPAN,IPOS,IQ,IQCLU,IT,ITCLU,IWRI,JTOP,
     &        LM,MS
      REAL*8 RVEC(3)
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C >> CONC
C        FMT02 = '(I10,3I10,10(I6,F8.5))'
C >> QBAS,VREF,RMTREF,VLMMAD_HOST,QMTET,QMPHI
      FMT03 = '(I10,1P,3E22.14)'
C >> VLMMAD_HOST
C        FMT04 = '(6X,2I4,1P,3E22.14)'
C >> R1,DX,RMT,RWS
      FMT05 = '(1X,I4,1P,2E22.14,2(I5,1P,E22.14))'
C >> NVALT
      FMT06 = '(1X,I4,5X,A,I10,I10,I10,9X,I4)'
C >> V,B,...
      FMT07 = '(1P,5E22.14)'
C >> QMVEC
C         FMT08 = '(A10,3(1P,E22.14))'
C >> ECORE
C         FMT09 = '(3X,I2,1P,2E22.14)'
C >> N5VEC_QCLU,IQ_QCLU
      FMT10 = '(I10,I9,4I4,I10)'
C >> IM_QCLU,NOCC,ITCLU_OCC,X_OCC
      FMT11 = '(3I10,(I10,F10.5))'
C-----------------------------------------------------------------------
C
      TROUBLE = .FALSE.
      IF ( IPRINT.GE.3 ) THEN
         CHECK = .TRUE.
         IWRI = 6
      ELSE
         CHECK = .FALSE.
         IWRI = 0
      END IF
      LINE = '          '
C-----------------------------------------------------------------------
      IOPOT = IFILPOT_CLU
C
      OPEN (IOPOT,FILE=POTFIL_CLU(1:LPOTFIL_CLU),ERR=100)
C
      REWIND IOPOT
      CALL READKWSTR(IOPOT,'HEADER    ',HEADER,LINE,80,IWRI,1)
      WRITE (6,'(10X,A)') HEADER
      LINE = '          '
      CALL READKWSTR(IOPOT,'TITLE     ',TITLE,LINE,80,IWRI,1)
      LTITLE = LEN_TRIM(TITLE)
      CALL READKWSTR(IOPOT,'SYSTEM    ',SYSTEM,LINE,80,IWRI,1)
      LSYSTEM = LEN_TRIM(SYSTEM)
C
C ======================================================================
C            STEP 1:    information specifying cluster occupation
C ======================================================================
C
      IF ( STEP.EQ.1 ) THEN
C
         CALL READKWINT(IOPOT,'NQCLU     ',NQCLU,0,IWRI,1)
         CALL READKWINT(IOPOT,'NTCLU     ',NTCLU,0,IWRI,1)
         CALL READKWINT(IOPOT,'NMCLU     ',NMCLU,0,IWRI,1)
C
C-----------------------------------------------------------------------
         NQCLUMAX = NQCLU
         NTCLUMAX = NTCLU
C-----------------------------------------------------------------------
C
         ALLOCATE (N5VEC_QCLU(5,NQCLUMAX),IQ_QCLU(NQCLUMAX))
         ALLOCATE (FORCE_QCLU(3,NQCLUMAX),IM_QCLU(NQCLUMAX))
         ALLOCATE (DQBAS_QCLU(3,NQCLUMAX),QBAS_QCLU(3,NQCLUMAX))
         ALLOCATE (NO_QCLU(NQCLUMAX),QBAS0_QCLU(3,NQCLUMAX))
         ALLOCATE (ITCLU_OQCLU(NTCLUMAX,NQCLUMAX),X_TCLU(NTCLUMAX))
C
         ALLOCATE (TXT_TCLU(NTCLUMAX),Z_TCLU(NTCLUMAX))
         ALLOCATE (NCOR_TCLU(NTCLUMAX),NVAL_TCLU(NTCLUMAX))
         ALLOCATE (NSEMCORSHL_TCLU(NTCLUMAX))
C
C ======================================================================
C
         CALL READKWSTR(IOPOT,'INFO      ',INFO,LINE,80,IWRI,0)
         LINFO = LEN_TRIM(INFO)
         CALL READKWSTR(IOPOT,'SCFSTATUS ',SCFSTATUS_CLU,LINE,10,IWRI,0)
C
C ======================================================================
         CALL POSFIL(IOPOT,'SITES     ',IPOS,1)
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQCLU = 1,NQCLU
            READ (IOPOT,*) IDUMMY,N5VEC_QCLU(1:5,IQCLU),IQ_QCLU(IQCLU)
            IF ( CHECK ) WRITE (6,FMT=FMT10) IDUMMY,
     &                          N5VEC_QCLU(1:5,IQCLU),IQ_QCLU(IQCLU)
            IF ( IQ_QCLU(IQCLU).LE.NQ_L ) THEN
               TROUBLE = .TRUE.
               WRITE (6,*) 
     &                  'No cluster sites within LEFT BULK allowed !!!!'
            END IF
            IF ( IQ_QCLU(IQCLU).GE.(NQHOST-NQ_R+1) ) THEN
               TROUBLE = .TRUE.
               WRITE (6,*) 
     &                 'No cluster sites within RIGHT BULK allowed !!!!'
            END IF
         END DO
C
C ======================================================================
         CALL POSFIL(IOPOT,'SHIFTS    ',IPOS,1)
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQCLU = 1,NQCLU
            READ (IOPOT,*) IDUMMY,DQBAS_QCLU(1:3,IQCLU)
            IF ( CHECK ) WRITE (6,FMT03) IDUMMY,DQBAS_QCLU(1:3,IQCLU)
C
            CALL RVECLCIVB(5,N5VEC_QCLU(1,IQCLU),A5BAS,RVEC)
C
            QBAS0_QCLU(1:3,IQCLU) = RVEC(1:3) + QBAS(1:3,IQ_QCLU(IQCLU))
C
            QBAS_QCLU(1:3,IQCLU) = QBAS0_QCLU(1:3,IQCLU)
     &                             + DQBAS_QCLU(1:3,IQCLU)
C
         END DO
C
C ======================================================================
         CALL POSFIL(IOPOT,'FORCES    ',IPOS,0)
         IF ( IPOS.NE.0 ) THEN
            READ (IOPOT,'(A80)') LINE
            IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
            READ (IOPOT,'(A80)') LINE
            IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
            DO IQCLU = 1,NQCLU
               READ (IOPOT,*) IDUMMY,FORCE_QCLU(1:3,IQCLU)
               IF ( CHECK ) WRITE (6,FMT03) IDUMMY,FORCE_QCLU(1:3,IQCLU)
            END DO
         ELSE
            FORCE_QCLU(1:3,1:NQCLU) = 0D0
         END IF
C
C ======================================================================
C
         CALL POSFIL(IOPOT,'OCCUPATION',IPOS,1)
C
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO IQCLU = 1,NQCLU
            READ (IOPOT,*) IDUMMY,IM_QCLU(IQCLU),NO_QCLU(IQCLU),
     &                     (ITCLU_OQCLU(IO,IQCLU),
     &                     X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                     NO_QCLU(IQCLU))
C
            IF ( CHECK ) WRITE (6,FMT=FMT11) IDUMMY,IM_QCLU(IQCLU),
     &                          NO_QCLU(IQCLU),
     &                          (ITCLU_OQCLU(IO,IQCLU),X_TCLU
     &                          (ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                          NO_QCLU(IQCLU))
         END DO
C
C ======================================================================
C
         CALL POSFIL(IOPOT,'MAGNETISAT',IPOS,1)
C
         READ (IOPOT,'(A80)') LINE
         IF ( SCFSTATUS_CLU(1:5).EQ.'START' ) READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         KMROT = 0
         DO IQ = NQHOST + 1,NQHOST + NQCLU
            READ (IOPOT,*) IDUMMY,QMTET(IQ),QMPHI(IQ)
            IF ( CHECK ) WRITE (6,FMT=FMT03) IDUMMY,QMTET(IQ),QMPHI(IQ)
            IF ( ABS(QMTET(IQ)).GT.1D-6 ) KMROT = 1
            IF ( ABS(QMPHI(IQ)).GT.1D-6 ) KMROT = 1
         END DO
C
C ======================================================================
C
         IF ( SCFSTATUS_CLU(1:5).EQ.'START' ) THEN
C
C ======================================================================
            DO IQCLU = 1,NQCLU
C
C               IQEXT = NQHOST + IQCLU
C               IQHOST = IQ_QCLU(IQCLU)
C
C ----------------------------------------------------------------------
C
Cc               IMHOST = IMQ(IQHOST)
Cc               IMQ(IQEXT) = IMHOST
Cc               IMEXT = IMQ(IQEXT)
C
Cc               R(1,IMEXT) = R(1,IMHOST)
Cc               DX(IMEXT) = DX(IMHOST)
Cc               JRMT(IMEXT) = JRMT(IMHOST)
Cc               RMT(IMEXT) = RMT(IMHOST)
Cc               JRWS(IMEXT) = JRWS(IMHOST)
Cc               RWS(IMEXT) = RWS(IMHOST)
C ----------------------------------------------------------------------
Cc               IF ( FULLPOT ) THEN
Cc                  JRNS1(IMEXT) = JRNS1(IMHOST)
Cc                  JRCRI(IMEXT) = JRCRI(IMHOST)
Cc                  NPAN(IMEXT) = NPAN(IMHOST)
Cc                  DO IPAN = 1,NPAN(IMEXT)
Cc                     JRCUT(IPAN,IMEXT) = JRCUT(IPAN,IMHOST)
Cc                  END DO
Cc               END IF
C
            END DO
C
C ======================================================================
C                                                  cluster specific mesh
         ELSE IF ( NMCLU.NE.0 ) THEN
C
            CALL POSFIL(IOPOT,'MESH INFOR',IPOS,1)
            CALL READKWSTR(IOPOT,'MESH-TYPE ',RMESHTYPE,'            ',
     &                     12,IWRI,1)
C ----------------------------------------------------------------------
            IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
               READ (IOPOT,'(A80)') LINE
               IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
               DO IM = NMHOST + 1,NMHOST + NMCLU
                  READ (IOPOT,*) IDUMMY,R(1,IM),DX(IM),JRMT(IM),RMT(IM),
     &                           JRWS(IM),RWS(IM)
                  IF ( CHECK ) WRITE (6,FMT=FMT05) IDUMMY,R(1,IM),DX(IM)
     &                                ,JRMT(IM),RMT(IM),JRWS(IM),RWS(IM)
               END DO
            ELSE
               READ (6,*) STR10,RMESHTYPE
               STOP 'in <POTRD>  mesh-type  not allowed'
            END IF
C ----------------------------------------------------------------------
            IF ( FULLPOT ) THEN
               READ (IOPOT,*) STR10
               DO IM = NMHOST + 1,NMHOST + NMCLU
                  READ (IOPOT,*) IDUMMY,JRNS1(IM),JRCRI(IM),NPAN(IM),
     &                           (JRCUT(IPAN,IM),IPAN=1,NPAN(IM))
               END DO
            END IF
C ======================================================================
C
         END IF
C
C ======================================================================
         CALL POSFIL(IOPOT,'TYPES     ',IPOS,1)
C
C                                                  cluster specific mesh
C ======================================================================
C
         READ (IOPOT,'(A80)') LINE
         IF ( CHECK ) WRITE (6,'(A)') LINE(1:LEN_TRIM(LINE))
         DO ITCLU = 1,NTCLU
            READ (IOPOT,*) IDUMMY,TXT_TCLU(ITCLU),Z_TCLU(ITCLU),
     &                     NCOR_TCLU(ITCLU),NVAL_TCLU(ITCLU),
     &                     NSEMCORSHL_TCLU(ITCLU)
            IF ( CHECK ) WRITE (6,FMT=FMT06) ITCLU,TXT_TCLU(ITCLU),
     &                          Z_TCLU(ITCLU),NCOR_TCLU(ITCLU),
     &                          NVAL_TCLU(ITCLU),NSEMCORSHL_TCLU(ITCLU)
         END DO
C======================================================================-
C
      END IF
C
      IF ( TROUBLE ) STOP 'in <CLUPOTRD>'
      IF ( STEP.EQ.1 .OR. SCFSTATUS_CLU(1:5).EQ.'START' ) RETURN
C
C======================================================================-
C                       pass information to global tables
C======================================================================-
      DO ITCLU = 1,NTCLU
         IT = NTHOST + ITCLU
         TXT_T(IT) = TXT_TCLU(ITCLU)
         Z(IT) = Z_TCLU(ITCLU)
         CONC(IT) = X_TCLU(ITCLU)
         NCORT(IT) = NCOR_TCLU(ITCLU)
         NVALT(IT) = NVAL_TCLU(ITCLU)
         NSEMCORSHLT(IT) = NSEMCORSHL_TCLU(ITCLU)
      END DO
C
C ======================================================================
C            STEP 2:    potential, charge, etc
C ======================================================================
C
      CALL POSFIL(IOPOT,'POTENTIAL ',IPOS,1)
      DO IT = NTHOST + 1,NTHOST + NTCLU
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            JTOP = JRCRI(IM)
         ELSE
            JTOP = JRWS(IM)
         END IF
         READ (IOPOT,*) STR10,IDUMMY
         READ (IOPOT,FMT=FMT07) (VT(I,IT),I=1,JTOP)
         READ (IOPOT,FMT=FMT07) (BT(I,IT),I=1,JTOP)
         IF ( CHECK ) WRITE (6,*) 'IT,IM,JTOP',IT,IM,JTOP
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            DO LM = 1,NLMFPMAX
               KLMFP(LM,IT) = 0
            END DO
            KLMFP(1,IT) = 1
            LMIFP(1,IT) = 1
            READ (IOPOT,*) STR10,NFPT(IT)
            READ (IOPOT,*) STR10,NLMFPT(IT)
            DO IFP = 2,NFPT(IT)
               READ (IOPOT,*) STR10,LM
               KLMFP(LM,IT) = 1
               LMIFP(IFP,IT) = LM
               READ (IOPOT,FMT=FMT07) (VNST(I,LM,IT),I=JRNS1(IM),JTOP)
               READ (IOPOT,FMT=FMT07) (BNST(I,LM,IT),I=JRNS1(IM),JTOP)
            END DO
         END IF
C ----------------------------------------------------------------------
         READ (IOPOT,*) STR10
      END DO
C ======================================================================
      CALL POSFIL(IOPOT,'CHARGE    ',IPOS,1)
      DO IT = NTHOST + 1,NTHOST + NTCLU
         READ (IOPOT,*) STR10,IDUMMY
         IM = IMT(IT)
         READ (IOPOT,FMT=FMT07) (RHOCHR(I,IT),I=1,JRWS(IM))
         READ (IOPOT,FMT=FMT07) (RHOSPN(I,IT),I=1,JRWS(IM))
         READ (IOPOT,*) STR10
C ----------------------------------------------------------------------
C       RHO2NS added to potfile Nov. 2009 - allow reading older potfiles
         IF ( FULLPOT .AND. STR10.NE.'==========' ) THEN
            READ (IOPOT,*)
            READ (IOPOT,*)
            DO LM = 1,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  READ (IOPOT,*)
                  READ (IOPOT,FMT=FMT07) (RHO2NS(I,LM,IT,1),I=1,JTOP)
                  READ (IOPOT,FMT=FMT07) (RHO2NS(I,LM,IT,2),I=1,JTOP)
               END IF
            END DO
            READ (IOPOT,*) STR10
         END IF
C ----------------------------------------------------------------------
      END DO
C ======================================================================
      CALL POSFIL(IOPOT,'MOMENTS   ',IPOS,1)
      DO IT = NTHOST + 1,NTHOST + NTCLU
         READ (IOPOT,*) STR10,IDUMMY
         READ (IOPOT,FMT=FMT07) QEL(IT),OBS_T(0,IDOS,IT),
     &                          OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),
     &                          OBS_T(0,IHFF,IT)
C
         READ (IOPOT,*) STR10
      END DO
C ======================================================================
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
C        CALL POSFIL(IOPOT,'VECTOR-POTENTIAL',IPOS,1)
         CALL POSFIL(IOPOT,'VECTOR-POT',IPOS,1)
         DO IT = NTHOST + 1,NTHOST + NTCLU
            READ (IOPOT,*) STR10,IDUMMY
            READ (IOPOT,*) STR10,LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
               IM = IMT(IT)
               DO MS = 1,2
                  READ (IOPOT,FMT=FMT07) (AOPT(I,MS,IT),I=1,JRWS(IM))
               END DO
            END IF
            READ (IOPOT,*) STR10
         END DO
CCC ======================================================================
CC      ELSE IF ( ORBPOL(1:4).EQ.'DMFT' ) THEN
CC         CALL POSFIL(IOPOT,'SELF-ENERG',IPOS,1)
CC         DO IT = NTHOST + 1,NTHOST + NTCLU
CC            READ (IOPOT,*) STR10,IDUMMY
CC            READ (IOPOT,*) STR10,LOPT(IT)
CCC ----------------------------------------------------------------------
CC            IF ( LOPT(IT).GT.0 ) THEN
CC            END IF
CC            READ (IOPOT,*) STR10
CC         END DO
CCC ======================================================================
      END IF
C
      CLOSE (IOPOT)
C
      RETURN
C ======================================================================
 100  CONTINUE
      WRITE (6,*) 'ERROR opening  POTFIL_CLU=',POTFIL_CLU(1:LPOTFIL_CLU)
      STOP 'in <CLURDPOT>'
C
      END
C*==clupotwr.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUPOTWR
C   ********************************************************************
C   *                                                                  *
C   *  write the cluster potential information                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NETAB,EFERMI
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,SYSTEM_DIMENSION
      USE MOD_SCF,ONLY:SCFSTATUS_CLU,SCFVXC,SCFTOL,SCFMIX,RMSAVB,RMSAVV,
     &    SCFALG,ITRSCF
      USE MOD_CALCMODE,ONLY:ORBPOL,BLCOUPL,BREITINT,EXTFIELD,PROGNAME,
     &    NONMAG,LLOYD,SEMICORE
      USE MOD_SITES,ONLY:NQCLU,QMPHI,QMTET,N5VEC_QCLU,IQ_QCLU,QBAS_QCLU,
     &    DQBAS_QCLU,IM_QCLU,NO_QCLU,NQHOST,ITCLU_OQCLU,FORCE_QCLU
      USE MOD_RMESH,ONLY:JRCUT,RMESHTYPE,DX,NMCLU,JRNS1,JRCRI,NPAN,JRMT,
     &    JRWS,R,RMT,RWS,FULLPOT,NMHOST,SPHERCELL
      USE MOD_FILES,ONLY:INFO,LINFO,SYSTEM,LSYSTEM,TITLE,LTITLE,
     &    IFILPOT_CLU,POTFIL_CLU,LPOTFIL_CLU
      USE MOD_TYPES,ONLY:NTCLU,NTHOST,AOPT,VT,BT,VNST,BNST,RHOCHR,
     &    RHOSPN,QEL,OBS_T,IMT,LOPT,NFPT,NLMFPT,KLMFP,VMTZ,BEXT,X_TCLU,
     &    NVAL_TCLU,TXT_TCLU,Z_TCLU,NCOR_TCLU,NSEMCORSHL_TCLU,RHO2NS
      USE MOD_KSPACE,ONLY:NKTAB,IBZINT
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,IHFF
      IMPLICIT NONE
C*--CLUPOTWR461
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER*40 FMT03,FMT05,FMT06,FMT07,FMT08,FMT11,FMT12
      INTEGER I,IM,IO,IOPOT,IPAN,IQ,IQCLU,IT,ITCLU,JTOP,LM,MS
C
C*** End of declarations rewritten by SPAG
C
      SCFSTATUS_CLU = 'ITERATING'
C
C-----------------------------------------------------------------------
C >> CONC
C     FMT02 = '(I10,3I10,10(I6,F8.5))'
C >> QBAS,VREF,RMTREF,VLMMAD_HOST,QMTET,QMPHI
      FMT03 = '(I10,1P,3E22.14)'
C >> VLMMAD_HOST
C       FMT04 = '(6X,2I4,1P,3E22.14)'
C >> R1,DX,RMT,RWS
      FMT05 = '(1X,I4,1P,2E22.14,2(I5,1P,E22.14))'
C >> NVALT
      FMT06 = '(1X,I4,5X,A,I10,I10,I10,9X,I4)'
C >> V,B,...
      FMT07 = '(1P,5E22.14)'
C >> QMVEC
      FMT08 = '(A10,3(1P,E22.14))'
C >> ECORE
C        FMT09 = '(3X,I2,1P,2E22.14)'
C >> N5VEC_QCLU,IQ_QCLU
C        FMT10 = '(I10,I9,4I4,I10)'
C >> IM_QCLU,NOCC,ITCLU_OCC,X_OCC
      FMT11 = '(3I10,(I10,F10.5))'
C >> IM,JRNS1,JRCRI,NPAN,JRCUT
      FMT12 = '(I5,3I6,20I5)'
C-----------------------------------------------------------------------
C
      IOPOT = IFILPOT_CLU
C
      OPEN (UNIT=IOPOT,STATUS='UNKNOWN',FILE=POTFIL_CLU(1:LPOTFIL_CLU))
C
C ======================================================================
      WRITE (IOPOT,99015) 'HEADER    ',
     &                    'SPR-KKR cluster dataset created by ',PROGNAME
      WRITE (IOPOT,99014)
      WRITE (IOPOT,99011) 'TITLE     ',TITLE(1:LTITLE)
      WRITE (IOPOT,99010) 'SYSTEM    ',SYSTEM(1:LSYSTEM)
C ======================================================================
      WRITE (IOPOT,99003) 'GLOBAL SYSTEM PARAMETER'
      WRITE (IOPOT,99013) 'NQCLU     ',NQCLU
      WRITE (IOPOT,99013) 'NTCLU     ',NTCLU
      WRITE (IOPOT,99013) 'NMCLU     ',NMCLU
C ======================================================================
      WRITE (IOPOT,99003) 'SCF-INFO  '
      WRITE (IOPOT,99011) 'INFO      ',INFO(1:LINFO)
      WRITE (IOPOT,99011) 'SCFSTATUS ',SCFSTATUS_CLU
      WRITE (IOPOT,99012) 'FULLPOT   ',FULLPOT
      IF ( FULLPOT ) WRITE (IOPOT,99012) 'SPHERCELL ',SPHERCELL
C
      WRITE (IOPOT,99012) 'BREITINT  ',BREITINT
      WRITE (IOPOT,99012) 'NONMAG    ',NONMAG
      WRITE (IOPOT,99010) 'ORBPOL    ',ORBPOL
      WRITE (IOPOT,99012) 'EXTFIELD  ',EXTFIELD
      WRITE (IOPOT,99012) 'BLCOUPL   ',BLCOUPL
      WRITE (IOPOT,FMT=FMT08) 'BEXT      ',BEXT
      WRITE (IOPOT,99012) 'SEMICORE  ',SEMICORE
      WRITE (IOPOT,99012) 'LLOYD     ',LLOYD
      WRITE (IOPOT,99013) 'NE        ',NETAB
      WRITE (IOPOT,99013) 'IBZINT    ',IBZINT
      WRITE (IOPOT,99013) 'NKTAB     ',NKTAB
      WRITE (IOPOT,99010) 'XC-POT    ',SCFVXC
      WRITE (IOPOT,99010) 'SCF-ALG   ',SCFALG
      WRITE (IOPOT,99013) 'SCF-ITER  ',ITRSCF
      WRITE (IOPOT,FMT=FMT08) 'SCF-MIX   ',SCFMIX
      WRITE (IOPOT,FMT=FMT08) 'SCF-TOL   ',SCFTOL
      WRITE (IOPOT,FMT=FMT08) 'RMSAVV    ',RMSAVV
      WRITE (IOPOT,FMT=FMT08) 'RMSAVB    ',RMSAVB
      WRITE (IOPOT,FMT=FMT08) 'EF        ',EFERMI
      WRITE (IOPOT,FMT=FMT08) 'VMTZ      ',VMTZ
      WRITE (IOPOT,99003) 'LATTICE   '
      WRITE (IOPOT,99011) 'SYSDIM    ',SYSTEM_DIMENSION
      WRITE (IOPOT,99011) 'SYSTYPE   ',SYSTEM_TYPE
C ======================================================================
      WRITE (IOPOT,99003) 'SITES     '
      WRITE (IOPOT,99002) 'positions of cluster sites'//
     &            '       (QX,QY,QZ) given only for info include SHIFTS'
      WRITE (IOPOT,99002) 
     &                   '     IQCLU            N5VEC_IQCLU          IQ'
     &                   //'       QX        QY        QZ'
      DO IQCLU = 1,NQCLU
         WRITE (IOPOT,99001) IQCLU,N5VEC_QCLU(1:5,IQCLU),IQ_QCLU(IQCLU),
     &                       QBAS_QCLU(1:3,IQCLU)
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'SHIFTS    '
      WRITE (IOPOT,99002) 
     &'shift from original lattice site positions in case of relaxation'
      WRITE (IOPOT,99002) 
     &           '     IQCLU        DQX             DQY             DQZ'
      DO IQCLU = 1,NQCLU
         WRITE (IOPOT,FMT=FMT03) IQCLU,DQBAS_QCLU(1:3,IQCLU)
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'FORCES    '
      WRITE (IOPOT,99002) 'forces acting on the cluster atoms '
      WRITE (IOPOT,99002) 
     &           '     IQCLU        FQX             FQY             FQZ'
      DO IQCLU = 1,NQCLU
         WRITE (IOPOT,FMT=FMT03) IQCLU,FORCE_QCLU(1:3,IQCLU)
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'OCCUPATION'
      WRITE (IOPOT,99002) 
     &              '     IQCLU    IMQCLU    NOQCLU   ITOQCLU      CONC'
      DO IQCLU = 1,NQCLU
         WRITE (IOPOT,FMT=FMT11) IQCLU,IM_QCLU(IQCLU),NO_QCLU(IQCLU),
     &                           (ITCLU_OQCLU(IO,IQCLU),
     &                           X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                           NO_QCLU(IQCLU))
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'MAGNETISATION DIRECTION'
      WRITE (IOPOT,99004)
      DO IQ = NQHOST + 1,NQHOST + NQCLU
         WRITE (IOPOT,FMT=FMT03) IQ,QMTET(IQ),QMPHI(IQ)
      END DO
C ======================================================================
C                                                  cluster specific mesh
      IF ( NMCLU.NE.0 ) THEN
C
         WRITE (IOPOT,99003) 'MESH INFORMATION'
         WRITE (IOPOT,99010) 'MESH-TYPE ',RMESHTYPE
         IF ( RMESHTYPE.EQ.'EXPONENTIAL ' ) THEN
            WRITE (IOPOT,99005)
            DO IM = NMHOST + 1,NMHOST + NMCLU
               WRITE (IOPOT,FMT=FMT05) IM,R(1,IM),DX(IM),JRMT(IM),
     &                                 RMT(IM),JRWS(IM),RWS(IM)
            END DO
         ELSE
            WRITE (6,99011) 'MESH-TYPE ',RMESHTYPE
            STOP 'in <CLUPOTWR>  mesh-type  not allowed'
         END IF
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            WRITE (IOPOT,99006)
            DO IM = NMHOST + 1,NMHOST + NMCLU
               WRITE (IOPOT,FMT=FMT12) IM,JRNS1(IM),JRCRI(IM),NPAN(IM),
     &                                 (JRCUT(IPAN,IM),IPAN=1,NPAN(IM))
            END DO
         END IF
C ----------------------------------------------------------------------
      END IF
C                                                  cluster specific mesh
C ======================================================================
C
C ======================================================================
      WRITE (IOPOT,99003) 'TYPES'
      WRITE (IOPOT,99007)
      DO ITCLU = 1,NTCLU
         WRITE (IOPOT,FMT=FMT06) ITCLU,TXT_TCLU(ITCLU),Z_TCLU(ITCLU),
     &                           NCOR_TCLU(ITCLU),NVAL_TCLU(ITCLU),
     &                           NSEMCORSHL_TCLU(ITCLU)
      END DO
C======================================================================-
      WRITE (IOPOT,99003) 'POTENTIAL'
      DO IT = NTHOST + 1,NTHOST + NTCLU
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            JTOP = JRCRI(IM)
         ELSE
            JTOP = JRWS(IM)
         END IF
         WRITE (IOPOT,99013) 'TYPE      ',IT
         WRITE (IOPOT,FMT=FMT07) (VT(I,IT),I=1,JTOP)
         WRITE (IOPOT,FMT=FMT07) (BT(I,IT),I=1,JTOP)
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            WRITE (IOPOT,99013) 'NFP       ',NFPT(IT)
            WRITE (IOPOT,99013) 'NLMFPT    ',NLMFPT(IT)
            DO LM = 2,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  WRITE (IOPOT,99013) 'LM        ',LM
                  WRITE (IOPOT,FMT=FMT07)
     &                   (VNST(I,LM,IT),I=JRNS1(IM),JTOP)
                  WRITE (IOPOT,FMT=FMT07)
     &                   (BNST(I,LM,IT),I=JRNS1(IM),JTOP)
               END IF
            END DO
         END IF
C ----------------------------------------------------------------------
         WRITE (IOPOT,99008)
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'CHARGE'
      DO IT = NTHOST + 1,NTHOST + NTCLU
         WRITE (IOPOT,99013) 'TYPE      ',IT
         IM = IMT(IT)
         WRITE (IOPOT,FMT=FMT07) (RHOCHR(I,IT),I=1,JRWS(IM))
         WRITE (IOPOT,FMT=FMT07) (RHOSPN(I,IT),I=1,JRWS(IM))
C ----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            WRITE (IOPOT,99009)
            WRITE (IOPOT,99013) 'NFP       ',NFPT(IT)
            WRITE (IOPOT,99013) 'NLMFPT    ',NLMFPT(IT)
            DO LM = 1,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  WRITE (IOPOT,99013) 'LM        ',LM
                  WRITE (IOPOT,FMT=FMT07) (RHO2NS(I,LM,IT,1),I=1,JTOP)
                  WRITE (IOPOT,FMT=FMT07) (RHO2NS(I,LM,IT,2),I=1,JTOP)
               END IF
            END DO
         END IF
C ----------------------------------------------------------------------
         WRITE (IOPOT,99008)
      END DO
C ======================================================================
      WRITE (IOPOT,99003) 'MOMENTS'
      DO IT = NTHOST + 1,NTHOST + NTCLU
         WRITE (IOPOT,99013) 'TYPE      ',IT
         WRITE (IOPOT,FMT=FMT07) QEL(IT),OBS_T(0,IDOS,IT),
     &                           OBS_T(0,ISMT,IT),OBS_T(0,IOMT,IT),
     &                           OBS_T(0,IHFF,IT)
C
         WRITE (IOPOT,99008)
      END DO
C ======================================================================
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
         WRITE (IOPOT,99003) 'VECTOR-POTENTIAL'
         DO IT = NTHOST + 1,NTHOST + NTCLU
            WRITE (IOPOT,99013) 'TYPE      ',IT
            WRITE (IOPOT,99013) 'LOPT      ',LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
               IM = IMT(IT)
               DO MS = 1,2
                  WRITE (IOPOT,FMT=FMT07) (AOPT(I,MS,IT),I=1,JRWS(IM))
               END DO
            END IF
            WRITE (IOPOT,99008)
         END DO
C ======================================================================
      ELSE IF ( ORBPOL(1:4).EQ.'DMFT' ) THEN
         WRITE (IOPOT,99003) 'SELF-ENERGY'
         DO IT = NTHOST + 1,NTHOST + NTCLU
            WRITE (IOPOT,99013) 'TYPE      ',IT
            WRITE (IOPOT,99013) 'LOPT      ',LOPT(IT)
C ----------------------------------------------------------------------
            IF ( LOPT(IT).GT.0 ) THEN
            END IF
            WRITE (IOPOT,99008)
         END DO
C ======================================================================
      END IF
C ======================================================================
      WRITE (IOPOT,99003) 'CORE STATES'
C      DO IT = NTHOST + 1,NTHOST + NTCLU
C         WRITE (IOPOT,99015) 'TYPE      ',IT
C         WRITE (IOPOT,99015) 'NCORE     ',NCORE(IT)
C         IF ( NCORE(IT).GT.0 ) THEN
C            WRITE (IOPOT,99011)
C            DO IS = 1,NS
C               IQ = (IT-1)*NS + IS
C               DO IM = 1,NCORE(IT)
C                  WRITE (IOPOT,99025) LCORE(IM,IT),ECORE(IM,IQ)
C               END DO
C            END DO
C         END IF
C         WRITE (IOPOT,99010)
C      END DO
C ======================================================================
C
      CLOSE (IOPOT)
C
C ======================================================================
99001 FORMAT (I10,I9,4I4,I10,3X,3F10.5)
99002 FORMAT (A)
99003 FORMAT (79('*'),/,A)
99004 FORMAT (6X,'  IQ','      QMTET     ','      QMPHI ')
99005 FORMAT (1X,'  IM ','     R(1)     ','       DX        ',' JRMT',
     &        '      RMT       ',' JRWS','      RWS')
99006 FORMAT (3X,'IM',' JRNS1',' JRCRI','  NPAN','  JRCUT')
99007 FORMAT ('ITCLU',' TXT_TCLU','    Z_TCLU',' NCOR_TCLU',
     &        ' NVAL_TCLU',' NSEMCORSHL_TCLU')
99008 FORMAT (79('='))
99009 FORMAT (79('-'))
99010 FORMAT (A10,A)
99011 FORMAT (A10,'''',A,'''')
99012 FORMAT (A10,L1)
99013 FORMAT (A10,6I10)
99014 FORMAT (79('*'))
99015 FORMAT (79('*'),/,A,'''',A,A,'''')
      END
