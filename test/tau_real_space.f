C*==tau_real_space.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE TAU_REAL_SPACE(KTAUIJ,IPRINT,ERYD,TSSQ,MSSQ,TAUQ,TSST,
     &                          MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the TAU - matrix for a finite cluster of atoms        *
C   *  by inverting the real space KKR matrix  M = [1/t - G0]          *
C   *                                                                  *
C   *  the routine uses the  REAL spherical harmonics representation   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  revision 2014:                                                  *
C   *                                                                  *
C   *  NSCT == NQ  ISCT == IQ  =>  ISCT, NSCT removed                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_THERMAL,ONLY:UMAT_VT,FMAT_FT,X_VFT,NVIBRA,NFLUCT,NVIBFLU
      USE MOD_CONSTANTS,ONLY:C1,C0,A0_ANG
      USE MOD_SYMMETRY,ONLY:NSYM,SYMUNITARY,SYMACCEPTED,IQORGQP,DROT,
     &    NSYMACCEPTED,NSYMMAX
      USE MOD_CPA,ONLY:ITCPAMAX,NCPA,CPATOL
      USE MOD_FILES,ONLY:IOTMP,FOUND_SECTION,FOUND_INTEGER
      USE MOD_LATTICE,ONLY:ABAS,ALAT,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,
     &    ADAINV_R,ADAINV_I,SYSTEM_DIMENSION
      USE MOD_SITES,ONLY:NQ,NQMAX,ICPA,DROTQ,QBAS,IQAT,ITOQ,NOQ,NQ_L,
     &    NQ_R,NQHOST
      USE MOD_ANGMOM,ONLY:NL,NKMMAX,NLMAX,NKMQ,NLQ,NKM,NLM,NRGNT123TAB,
     &    WKM1,WKM2,IPIVKM
      USE MOD_TYPES,ONLY:NTMAX,CONC,NT
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_REAL_SPACE')
      REAL*8 TOL,C
      PARAMETER (TOL=1.0D-6,C=2.0D0*1.370373D+02)
      LOGICAL CHECKMG0,CHECKTAU
      PARAMETER (CHECKMG0=.FALSE.,CHECKTAU=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IPRINT,KTAUIJ
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CLUIPH(:,:),DMAT_VFT(:,:),DROTQCLU(:,:,:),DTIL_VFT(:,:)
     &           ,HLTAB(2*NLMAX-1),JLTAB(2*NLMAX-1),MSSQCLU(:,:,:),
     &           MSSTCLU(:,:,:),MSS_VFT(:,:),NLTAB(2*NLMAX-1),P,
     &           TAUQCLU(:,:,:),TAUT(:,:,:),TSSQCLU(:,:,:),
     &           TSSTCLU(:,:,:),TSS_FT(:,:),TSS_VFT(:,:),
     &           W1(NKMMAX,NKMMAX),WA(:,:),WB(:,:),WC(:,:),WD(:,:)
      REAL*8 CLURAD,CLURYLM_NNP(:,:),CPACHNG,DIR(3),DQCLU(:),DSSP(:),
     &       FL1L2(:,:),NORM,RGNT_CLU(:),RQCLU(:,:),RSSP(:,:),S,TEMP_LAT
      INTEGER I,I0,I1,I2,I2B,IA,IA_ERR,ICALL,ICPACLU(:),ICPAFLAG,ID,
     &        IDIR,IFLUCT,IFT,IG,IGIJ(:,:),IL3RGNT(:),IND0QCLU(:),
     &        IND0QCLU0(:),IO,IPIV(:),IQ,IQ1,IQ2,IQA,IQB,IQCLU,IQCLU0,
     &        IQCLU1,IQCLU2,IQCLUB,IQCLUORG(:,:),IQCLUTAB(:,:),
     &        IQCLU_ISTQ(:,:),IQCNTR,IQ_IQCLU(:),ISDA4(:,:),ISITE,ISSP,
     &        ISSP2AB(:),ISSP2BB(:),ISSP4(:,:),ISSP5(:,:),ISSPDIR(:),
     &        ISYM,IT,IT0,ITCNTR,ITIMP,ITOQCLU(:,:),IVFT,IVIBRA,IVT,IX,
     &        I_TEMP_LAT,J,J0,JQ,JQCLU,JQCLU1,JQCLU2,JQCLUTAB(:,:),
     &        JSDA4(:,:),L1,L2,L3MAX,LM1,LM2,LM3RGNT(:),LMAX,M,M1,M2,N,
     &        NGIJ,NI,NIJSTAB(:),NIJSTABMAX,NJ,NKKR
      LOGICAL IMPURITY,MOL,SAMENLQ,SWAPTAUQ,SYMMETRIZE
      INTEGER NKMOUT,NKMPROJ,NKMQCLU(:),NLM3MAX,NLMQCLU(:),NLOUT,
     &        NLQCLU(:),NOQCLU(:),NQCLU,NQCLU_I,NQCLU_IQ(:),NQCLU_L,
     &        NQCLU_R,NRGNT(:),NRGNT12,NRGNT123,NRGNT123MAX,NRGNT12MAX,
     &        NSHLCLU,NSSP1,NSSP2A,NSSP2B,NSSP4(:),NSSP4MAX,NSSPABS,
     &        NSSPDIR,N_TEMP_LAT
      SAVE CLUIPH,CLURAD,CLURYLM_NNP,DQCLU,DROTQCLU,DSSP,FL1L2,ICPACLU,
     &     IGIJ,IL3RGNT,IMPURITY,IND0QCLU,IND0QCLU0,IQ1,IQ2,IQCLUORG,
     &     IQCLUTAB,IQCLU_ISTQ,IQCNTR,IQ_IQCLU,ISDA4,ISSP2AB,ISSP2BB,
     &     ISSP4,ISSP5,ISSPDIR,ITCNTR,ITIMP,ITOQCLU,JQCLUTAB,JSDA4,
     &     L3MAX,LM3RGNT,LMAX,MOL,NGIJ,NIJSTAB,NIJSTABMAX,NKKR,NKMOUT,
     &     NKMQCLU,NLM3MAX,NLMQCLU,NLQCLU,NOQCLU,NQCLU,NQCLU_IQ,NRGNT,
     &     NRGNT12,NRGNT123,NRGNT123MAX,NRGNT12MAX,NSHLCLU,NSSP1,NSSP2A,
     &     NSSP2B,NSSP4,NSSP4MAX,NSSPABS,NSSPDIR,RGNT_CLU,RQCLU,RSSP,
     &     SAMENLQ,SWAPTAUQ,TAUT
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/,IA_ERR/0/
C
      ALLOCATABLE RQCLU,DQCLU,IQ_IQCLU,IQCLUORG
      ALLOCATABLE RGNT_CLU,IL3RGNT,LM3RGNT,NRGNT,CLURYLM_NNP
      ALLOCATABLE IQCLU_ISTQ,NQCLU_IQ
      ALLOCATABLE WA,WB,WC,WD,IPIV,NLQCLU,NLMQCLU,NKMQCLU
      ALLOCATABLE RSSP,DSSP,NIJSTAB,IQCLUTAB,JQCLUTAB,ISSP2AB,ISSP2BB
      ALLOCATABLE NSSP4,ISSP4,ISDA4,JSDA4,ISSPDIR,ISSP5
      ALLOCATABLE CLUIPH,IGIJ,IND0QCLU,IND0QCLU0,FL1L2
      ALLOCATABLE MSSQCLU,MSSTCLU,TSSQCLU,TSSTCLU,TAUQCLU
      ALLOCATABLE NOQCLU,ITOQCLU,ICPACLU,DROTQCLU
      ALLOCATABLE TAUT,TSS_FT,TSS_VFT,MSS_VFT,DMAT_VFT,DTIL_VFT
C
      ALLOCATE (TSS_FT(NKMMAX,NKMMAX),TSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (DMAT_VFT(NKMMAX,NKMMAX),DTIL_VFT(NKMMAX,NKMMAX))
C
C=======================================================================
      ICALL = ICALL + 1
Cc      SYMMETRIZE = .TRUE.
      SYMMETRIZE = .FALSE.
C
C   ********************************************************************
C                   INITIALISATION - START
C   ********************************************************************
C
      IF ( ICALL.EQ.1 ) THEN
C
C---------------------------------- DEAL WITH TEMPERATURE - IF REQUESTED
C
         TEMP_LAT = 0D0
         N_TEMP_LAT = 1
         I_TEMP_LAT = 1
C
         CALL THERMAL_INIT(0,N_TEMP_LAT,TEMP_LAT)
C
C ======================================================================
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C            and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ
C
         CALL THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
         M = NKMMAX
C
         ALLOCATE (TAUT(M,M,NT),STAT=IA_ERR)
C
         WRITE (6,99006)
C
C ----------------------------------------------------------------------
C
         L3MAX = 2*(NLMAX-1)
         NLM3MAX = (L3MAX+1)*(L3MAX+1)
C
         WRITE (6,99004)
C
C--------- initialize TAUQ -- necessary if NLQ differs from site to site
C
         CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C
         CALL INPUT_FIND_SECTION('TAU',0)
C ------------------ allow to use the old section CONTROL instead of TAU
         IF ( .NOT.FOUND_SECTION ) CALL INPUT_FIND_SECTION('CONTROL',0)
C
         NLOUT = 4
C
         CALL SECTION_SET_INTEGER('NLOUT',NLOUT,9999,0)
C
         IF ( NL.LT.4 ) NLOUT = 3
         NLOUT = MIN(NL,NLOUT)
C
         NKMOUT = 2*NLOUT**2
C
         CALL SECTION_FIND_KEYWORD('IMPURITY',IMPURITY)
         CALL SECTION_FIND_KEYWORD('MOL',MOL)
C
         IF ( IMPURITY ) CALL SECTION_SET_INTEGER('ITIMP',ITIMP,9999,1)
C
         IF ( MOL ) THEN
            IQCNTR = 1
            ITCNTR = ITOQ(1,IQCNTR)
            IQ1 = 1
            IQ2 = NQ
         ELSE
            CALL SECTION_SET_INTEGER('IQCNTR',IQCNTR,9999,1)
            ITCNTR = ITOQ(1,IQCNTR)
            IQ1 = IQCNTR
            IQ2 = IQCNTR
         END IF
C
C------------------------------------------------- set NSHLCLU or CLURAD
C
         NSHLCLU = 0
         CALL SECTION_SET_INTEGER('NSHLCLU',NSHLCLU,9999,0)
C
         IF ( .NOT.FOUND_INTEGER )
     &        CALL SECTION_SET_REAL('CLURAD',CLURAD,9999D0,1)
C
C ----------------------------------------------------------------------
C
         WRITE (6,99003)
C
C------------------------------------------- set up factor (-1)**(l1-l2)
C
         ALLOCATE (FL1L2((2*NLMAX-1)**2,(2*NLMAX-1)**2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FL1L2')
C
         LM1 = 0
         DO L1 = 0,2*(NLMAX-1)
            DO M1 = -L1, + L1
               LM1 = LM1 + 1
               LM2 = 0
               DO L2 = 0,2*(NLMAX-1)
                  DO M2 = -L2, + L2
                     LM2 = LM2 + 1
                     FL1L2(LM1,LM2) = (-1D0)**(L1-L2)
                  END DO
               END DO
            END DO
         END DO
C
C-----------------------------------------------------------------------
C                 - allocate storage for real GAUNT coefficients
C                 - set up modified GAUNT's
C-----------------------------------------------------------------------
C
         N = NRGNT123TAB(NL)
C
         NRGNT123MAX = N
         ALLOCATE (RGNT_CLU(N),IL3RGNT(N),LM3RGNT(N),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ALLOC: NRGNT123MAX')
C
         NRGNT12MAX = NL**4
         ALLOCATE (NRGNT(NRGNT12MAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ALLOC: NRGNT12MAX')
C
         LMAX = NL - 1
C
         CALL CLUGAUNT_RYLM(LMAX,RGNT_CLU,NRGNT,IL3RGNT,LM3RGNT,NRGNT12,
     &                      NRGNT123,NRGNT12MAX,NRGNT123MAX)
C
C-----------------------------------------------------------------------
C           - determine cluster size  NQCLU and number of
C             atomic shells  NSSHCLUS according to the geometry
C           - allocate related storage and read data
C-----------------------------------------------------------------------
C
         CALL CLUSSITES(IOTMP,IPRINT,MOL,SYSTEM_DIMENSION,ABAS,ABAS_L,
     &                  ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,ADAINV_R,QBAS,
     &                  CLURAD,IQCNTR,NQCLU,NQCLU_L,NQCLU_I,NQCLU_R,
     &                  NSHLCLU,NQHOST,NQ_L,NQ_R,NQMAX)
C
C -------------------------------------------------------------- GENERAL
C
         WRITE (6,99001)
         WRITE (6,99005) 'cluster sites    ','NQCLU   ',NQCLU
         WRITE (6,99005) 'cluster shells   ','NSHLCLU ',NSHLCLU
C
         ALLOCATE (IQCLUORG(NSYMMAX,NQCLU))
         ALLOCATE (RQCLU(3,NQCLU),DQCLU(NQCLU))
         ALLOCATE (IQ_IQCLU(NQCLU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQCLUS')
C
C--------- Additional precautions regarding input of cluster coordinates
C
         READ (IOTMP,IOSTAT=IA_ERR,ERR=50,END=50)
     &         ((RQCLU(J,I),J=1,3),DQCLU(I),IQ_IQCLU(I),I=1,NQCLU)
 50      CONTINUE
         IF ( IA_ERR.NE.0 ) THEN
            WRITE (6,'(/,a,2i5)') 'ia_err =',IA_ERR
            CALL FLUSH(6)
            WRITE (6,'(a,i5)') 'NQCLU =',NQCLU
            DO I = 1,NQCLU
               WRITE (6,'(i4,4x,3f10.5,f12.3,i7)') I,(RQCLU(J,I),J=1,3),
     &                DQCLU(I),IQ_IQCLU(I)
            END DO
            CALL STOP_MESSAGE(ROUTINE,'READ: RQCLU etc')
         END IF
         CLOSE (IOTMP)
C
C----------------- Print finite cluster for which TAU matrix is evaluated
C
         WRITE (6,'(//,a,i5,/)') 'Real-space cluster employed: NQCLU =',
     &                           NQCLU
         WRITE (6,'(a,4x,9(a))') ' IQCLU','  RQCLU(1)','  RQCLU(2)',
     &                           '  RQCLU(3)','       DQCLU',
     &                           '   IQ_IQCLU','  ITOQ'
         DO IQCLU = 1,NQCLU
            WRITE (6,'(i6,4x,3f10.5,f12.3,i9,3x,11i4)') IQCLU,
     &             (RQCLU(J,IQCLU),J=1,3),DQCLU(IQCLU),IQ_IQCLU(IQCLU),
     &             (ITOQ(IO,IQ_IQCLU(IQCLU)),IO=1,
     &             MIN(11,NOQ(IQ_IQCLU(IQCLU))))
            IF ( NOQ(IQ_IQCLU(IQCLU)).GT.11 ) WRITE (6,'(64x,11i4)')
     &           (ITOQ(IO,IQ_IQCLU(IQCLU)),IO=12,NOQ(IQ_IQCLU(IQCLU)))
         END DO
         WRITE (6,'(//)')
C
C--------------------------------------------------------- supply NLQCLU
C
         ALLOCATE (NLQCLU(NQCLU))
C
         SWAPTAUQ = .FALSE.
C
         DO IQCLU = 1,NQCLU
C
            IF ( IQ_IQCLU(IQCLU).NE.IQCLU ) SWAPTAUQ = .TRUE.
C
            NLQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))
C
            IQA = IQ_IQCLU(IQCLU)
            DO ISYM = 1,NSYM
               IF ( SYMACCEPTED(ISYM) ) THEN
                  IQB = IQORGQP(ISYM,IQA)
                  DO JQCLU = 1,NQCLU
                     IF ( IQB.EQ.IQ_IQCLU(JQCLU) ) IQCLUB = JQCLU
                  END DO
                  IQCLUORG(ISYM,IQCLUB) = IQCLU
Cc                  IQCLUORG(ISYM,IQCLU) = IQCLUB
               END IF
            END DO
         END DO
C
C--------------------------------------- set up tables specifying blocks
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
         CALL CLUSIT2(IOTMP,NQCLU,RQCLU,IPRINT,NQCLU)
C
         REWIND IOTMP
         READ (IOTMP) NSSP1,NIJSTABMAX,NSSP2A,NSSP2B,NSSPABS,NSSP4MAX,
     &                NSSPDIR
C
         ALLOCATE (RSSP(3,0:NSSP1),IQCLUTAB(NIJSTABMAX,NSSP1))
         ALLOCATE (DSSP(0:NSSP1),JQCLUTAB(NIJSTABMAX,NSSP1))
         ALLOCATE (NIJSTAB(NSSP1),ISSP2AB(NSSP2A),ISSP2BB(NSSP2A))
         ALLOCATE (NSSP4(NSSPABS),ISSP4(NSSP4MAX,NSSPABS))
         ALLOCATE (ISDA4(NSSP4MAX,NSSPABS),ISSPDIR(NSSPDIR))
         ALLOCATE (JSDA4(NSSP4MAX,NSSPABS))
         ALLOCATE (ISSP5(NSSP4MAX,NSSPABS),STAT=IA_ERR)
C
C***********************************************************************
         RSSP = 999999
         DSSP = 999999
         IQCLUTAB = 999999
         JQCLUTAB = 999999
         NIJSTAB = 999999
         ISSP2AB = 999999
         ISSP2BB = 999999
         NSSP4 = 999999
         ISSP4 = 999999
         ISDA4 = 999999
         JSDA4 = 999999
         ISSPDIR = 999999
C
C        NKKR=999999
C        NSSP1=999999
C        NSSP2A=999999
C        NSSP2B=999999
C        NSSPDIR=999999
C        NSSPABS=999999
C
         ISSP5 = 999999
         NLQCLU = 999999
C***********************************************************************
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSP5')
C
         READ (IOTMP) ((RSSP(J,I),J=1,3),DSSP(I),NIJSTAB(I),I=1,NSSP1),
     &                ((IQCLUTAB(J,I),JQCLUTAB(J,I),J=1,NIJSTAB(I)),I=1,
     &                NSSP1),(ISSP2AB(I),ISSP2BB(I),I=1,NSSP2A),
     &                (NSSP4(I),I=1,NSSPABS),
     &                ((ISSP4(J,I),ISDA4(J,I),JSDA4(J,I),J=1,NSSP4(I)),
     &                I=1,NSSPABS),(ISSPDIR(I),I=1,NSSPDIR),
     &                ((ISSP5(J,I),J=1,NSSP4(I)),I=1,NSSPABS)
         CLOSE (IOTMP)
C
         WRITE (6,99005) 'S-S'' combinations','NSSP1   ',NSSP1
         WRITE (6,99005) 'S-S'' lengths     ','NSSPABS ',NSSPABS
         WRITE (6,99005) 'S-S'' directions  ','NSSPDIR ',NSSPDIR
         WRITE (6,99005) 'l-expansion      ','NL      ',NL,NLMAX
         WRITE (6,99005) 'Gaunt table I    ','NRGNT12 ',NRGNT12,
     &                   NRGNT12MAX
         WRITE (6,99005) 'Gaunt table II   ','NRGNT123',NRGNT123,
     &                   NRGNT123MAX
C
         ALLOCATE (IND0QCLU0(NQCLU),IGIJ(NQCLU,NQCLU),IND0QCLU(NQCLU))
         ALLOCATE (IQCLU_ISTQ(NQCLU,NQ),NQCLU_IQ(NQ))
         ALLOCATE (NOQCLU(NQCLU),NLMQCLU(NQCLU),NKMQCLU(NQCLU))
         ALLOCATE (ITOQCLU(NTMAX,NQCLU),ICPACLU(NQCLU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ITOQCLU')
         NKMQCLU = 999999
         IND0QCLU0 = 999999
C
         IF ( KMROT.GT.0 .AND. IREL.GE.3 ) THEN
            ALLOCATE (DROTQCLU(NKMMAX,NKMMAX,NQCLU),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROTQCLU')
         ELSE
            IF ( KMROT.NE.0 .AND. IREL.LE.2 ) THEN
               WRITE (6,*) 'inconsistent parameter settings'
               WRITE (6,*) 'IREL = ',IREL,'   KMROT = ',KMROT
               CALL STOP_MESSAGE(ROUTINE,
     &                'rotation of magn. (KMROT <> 0) only for IREL = 3'
     &                )
               STOP ' in <TAU_REAL_SPACE>'
            END IF
            ALLOCATE (DROTQCLU(1,1,1),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROTQCLU')
         END IF
C
C--------------------- set up tables specifying which cluster site IQCLU
C----------------------------------------- is occupied by scatterer IQ
C
         NQCLU_IQ(:) = 0
C
         DO IQCLU = 1,NQCLU
            IQ = IQ_IQCLU(IQCLU)
            NQCLU_IQ(IQ) = NQCLU_IQ(IQ) + 1
            IQCLU_ISTQ(NQCLU_IQ(IQ),IQ) = IQCLU
         END DO
C
C ----------------------------------------------------------------------
C
         DO IQ = 1,NQ
C
            DO ISITE = 1,NQCLU_IQ(IQ)
               IQCLU = IQCLU_ISTQ(ISITE,IQ)
               ICPACLU(IQCLU) = ICPA(IQ)
               NOQCLU(IQCLU) = NOQ(IQ)
               DO IO = 1,NOQ(IQ)
                  ITOQCLU(IO,IQCLU) = ITOQ(IO,IQ)
               END DO
C
               IF ( KMROT.NE.0 .AND. IREL.GE.3 )
     &              CALL CHANGEREP(NKM,NKMMAX,DROTQ(1,1,IQ),'REL>RLM',
     &              DROTQCLU(1,1,IQCLU))
C
            END DO
C
         END DO
C
C------------------------------------------ transform rotation matrices
C
         IF ( SYMMETRIZE .AND. IREL.GE.3 ) THEN
C
            DO ISYM = 1,NSYM
C
               CALL CHANGEREP(NKM,NKMMAX,DROT(1,1,ISYM),'REL>RLM',WKM2)
C
               CALL CMATCOP(NKM,NKMMAX,WKM2,DROT(1,1,ISYM))
C
            END DO
C
         END IF
C
C---------------------------------------------------- set up table  IGIJ
         IG = 0
         DO IA = 2,NSSPABS
            DO ID = 1,NSSP4(IA)
               IG = IG + 1
               IQCLU = ISDA4(ID,IA)
               JQCLU = JSDA4(ID,IA)
               IGIJ(IQCLU,JQCLU) = IG
            END DO
         END DO
C
         DO I2B = 1,NSSP2B
            IG = IG + 1
            IQCLU = IQCLUTAB(1,ISSP2BB(I2B))
            JQCLU = JQCLUTAB(1,ISSP2BB(I2B))
            IGIJ(IQCLU,JQCLU) = IG
         END DO
         NGIJ = IG
C
         DO I1 = 2,NSSP1
            IQCLU1 = IQCLUTAB(1,I1)
            JQCLU1 = JQCLUTAB(1,I1)
C
            DO I2 = 2,NIJSTAB(I1)
               IQCLU2 = IQCLUTAB(I2,I1)
               JQCLU2 = JQCLUTAB(I2,I1)
               IGIJ(IQCLU2,JQCLU2) = IGIJ(IQCLU1,JQCLU1)
            END DO
         END DO
         WRITE (6,99005) 'number of G(I,J) ','NGIJ    ',NGIJ
         WRITE (6,99002)
C
         IF ( IPRINT.GT.1 ) THEN
            DO IQCLU = 1,NQCLU
               DO JQCLU = 1,IQCLU - 1
                  WRITE (6,'(10X,A,2I4,2(A,I4))') 'IGIJ',IQCLU,JQCLU,
     &                   '  i,j:',IGIJ(IQCLU,JQCLU),' j,i:',
     &                   IGIJ(JQCLU,IQCLU)
               END DO
            END DO
         END IF
C
C --------------------------------------------------------- NLM-blocking
C
         NKKR = 0
         NLM = NL**2
         SAMENLQ = .TRUE.
C
         DO IQCLU = 1,NQCLU
            NLQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))
            NLMQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))**2
            IF ( NLQCLU(IQCLU).NE.NL ) SAMENLQ = .FALSE.
C
            NKMQCLU(IQCLU) = 2*NLMQCLU(IQCLU)
            IF ( IQCLU.EQ.1 ) THEN
               IND0QCLU(1) = 0
               IND0QCLU0(1) = 0
            ELSE
               IND0QCLU(IQCLU) = IND0QCLU(IQCLU-1) + NLMQCLU(IQCLU-1)
               IND0QCLU0(IQCLU) = IND0QCLU0(IQCLU-1) + NLM
            END IF
            NKKR = NKKR + NLMQCLU(IQCLU)
         END DO
C
C-------------------------------------- tabulate the spherical harmonics
C
         ALLOCATE (CLURYLM_NNP(NLM3MAX,NSSPDIR),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ALLOC: CLURYLM_NNP')
         CLURYLM_NNP = 999999
C
         ALLOCATE (CLUIPH(2*NLMAX-1,NSSPABS),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CLUIPH')
         CLUIPH = 999999
C
         DO IDIR = 1,NSSPDIR
            ISSP = ISSPDIR(IDIR)
            IF ( DABS(DSSP(ISSP)).GE.1D-6 ) THEN
               DO I = 1,3
                  DIR(I) = RSSP(I,ISSP)/DSSP(ISSP)
               END DO
            ELSE
               DIR(1:3) = 0D0
            END IF
C
            CLURYLM_NNP(1:NLM3MAX,IDIR) = CLURYLM_NNP(1:NLM3MAX,1)
C
            IF ( IDIR.NE.1 ) CALL CALC_RHPLM(DIR(1),DIR(2),DIR(3),
     &           CLURYLM_NNP(1,IDIR),L3MAX,NLM3MAX)
C
         END DO
C
C----------------------------------------------------- write RASMOL data
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'rasmol_cluster.pdb')
         WRITE (IOTMP,FMT=99009)
C
         S = 3*A0_ANG*ALAT
         I = 0
         DO IQCLU = 1,NQCLU
            IQ = IQ_IQCLU(IQCLU)
            I = I + 1
            WRITE (IOTMP,FMT=99010) I,I,(S*RQCLU(IX,IQCLU),IX=1,3),
     &                              DBLE(IQ)*3D0/DBLE(NQ)
         END DO
C
         CLOSE (IOTMP)
C
C--------------------------------------------------- write RASMOL script
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'rasmol_cluster.ras')
         WRITE (IOTMP,*) 'load ''rasmol_cluster.pdb'' '
         WRITE (IOTMP,*) 'color temperature'
         WRITE (IOTMP,*) 'set fontsize 20'
         WRITE (IOTMP,*) 'cpk   450'
         WRITE (IOTMP,*) 'select all '
         WRITE (IOTMP,*) 'set axes on '
C
         CLOSE (IOTMP)
C
      END IF
C **********************************************************************
C                   INITALISATION - END
C **********************************************************************
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
      CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
      CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
      ALLOCATE (TSSQCLU(NKMMAX,NKMMAX,NQCLU))
      ALLOCATE (MSSQCLU(NKMMAX,NKMMAX,NQCLU))
      ALLOCATE (WA(NKKR,NKKR),WB(NKKR,NKKR))
      ALLOCATE (WC(NKKR,NKKR),WD(NKKR,NKKR),IPIV(NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WC')
C
      IF ( NCPA.EQ.0 ) THEN
         ALLOCATE (MSSTCLU(1,1,1),TSSTCLU(1,1,1))
      ELSE
         ALLOCATE (TSSTCLU(NKMMAX,NKMMAX,NT))
         ALLOCATE (MSSTCLU(NKMMAX,NKMMAX,NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MSSTCLU')
C
         N = NKM
         M = NKMMAX
         DO IT = 1,NT
C
            CALL CHANGEREP(N,M,MSST(1,1,IT),'REL>RLM',MSSTCLU(1,1,IT))
C
            CALL CHANGEREP(N,M,TSST(1,1,IT),'REL>RLM',TSSTCLU(1,1,IT))
C
            IVT = (IT-1)*NVIBRA
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               WKM1(1:M,1:M) = UMAT_VT(1:M,1:M,IVT)
C
               CALL CHANGEREP(N,M,WKM1,'REL>RLM',UMAT_VT(1,1,IVT))
            END DO
C
         END DO
C
      END IF
C
C ----------------------------------------------------------------------
C----------------------------------------------------- scan scatterer IQ
C ------------------------------------------- MSSQ(IQ) -> MSSQCLU(IQCLU)
C-------------- transform MSS from (KAPPA,MUE) to (L,ML,MS) in necessary
C--------------------------------- copy MSS for equivalent cluster sites
C
      DO IQ = 1,NQ
C
         IF ( NQCLU_IQ(IQ).LE.0 ) CYCLE
C
         IF ( IREL.GE.2 ) THEN
C
            IQCLU0 = IQCLU_ISTQ(1,IQ)
C
            CALL CHANGEREP(NKM,NKMMAX,MSSQ(1,1,IQ),'REL>RLM',
     &                     MSSQCLU(1,1,IQCLU0))
C
            DO ISITE = 2,NQCLU_IQ(IQ)
               IQCLU = IQCLU_ISTQ(ISITE,IQ)
               MSSQCLU(1:NKM,1:NKM,IQCLU) = MSSQCLU(1:NKM,1:NKM,IQCLU0)
            END DO
C
         ELSE
C
            DO ISITE = 1,NQCLU_IQ(IQ)
               IQCLU = IQCLU_ISTQ(ISITE,IQ)
               MSSQCLU(1:NKM,1:NKM,IQCLU) = MSSQ(1:NKM,1:NKM,IQ)
            END DO
C
         END IF
C
      END DO
C
C------------------------ calculate free electron Greens function matrix
C                                                    stored as  WA = -G0
C
      IF ( CHECKMG0 ) CALL CLUCHKG0(WA,NKKR,NQCLU,NLMQCLU,NLQCLU,DSSP,
     &                              ISDA4,NSSP2A,ISSP2AB,ISSP2BB,ISSP4,
     &                              ISSP5,IQCLUTAB,JSDA4,JQCLUTAB,
     &                              NIJSTAB,NIJSTABMAX,NSSP4,NSSP4MAX,
     &                              NSSP1,NSSP2B,NSSPABS,NSSPDIR,NRGNT,
     &                              IL3RGNT,LM3RGNT,RGNT_CLU,LMAX,ERYD,
     &                              P,C,ALAT,NLM,IND0QCLU,IND0QCLU0,
     &                              SAMENLQ,CLUIPH,CLURYLM_NNP,JLTAB,
     &                              NLTAB,HLTAB,FL1L2,NLMAX,NLM3MAX,
     &                              NRGNT12MAX,NRGNT123MAX)
C
      WA = 0D0
C
      CALL CLUG0MAT(WA,NKKR,NQCLU,NLMQCLU,DSSP,ISDA4,NSSP2A,ISSP2AB,
     &              ISSP2BB,ISSP4,ISSP5,IQCLUTAB,JSDA4,JQCLUTAB,NIJSTAB,
     &              NIJSTABMAX,NSSP4,NSSP4MAX,NSSP1,NSSP2B,NSSPABS,
     &              NSSPDIR,NRGNT,IL3RGNT,LM3RGNT,RGNT_CLU,LMAX,ERYD,P,
     &              C,ALAT,NLM,IND0QCLU,IND0QCLU0,SAMENLQ,CLUIPH,
     &              CLURYLM_NNP,JLTAB,NLTAB,HLTAB,FL1L2,NLMAX,NLM3MAX,
     &              NRGNT12MAX,NRGNT123MAX)
C
C------------------------------------------ write Greens function matrix
C
      IF ( KTAUIJ.GE.2 ) THEN
C
         OPEN (UNIT=18,STATUS='SCRATCH',FORM='UNFORMATTED')
C
         DO JQ = 1,NQ
            IF ( NQCLU_IQ(JQ).NE.1 )
     &            CALL STOP_MESSAGE(ROUTINE,'NQCLU_IQ(JQ).NE.1')
            JQCLU = IQCLU_ISTQ(1,JQ)
            J0 = IND0QCLU(JQCLU)
            NJ = NLMQCLU(JQCLU)
            DO IQ = 1,NQ
               IF ( NQCLU_IQ(IQ).NE.1 )
     &               CALL STOP_MESSAGE(ROUTINE,'NQCLU_IQ(IQ).NE.1')
               IQCLU = IQCLU_ISTQ(1,IQ)
               I0 = IND0QCLU(IQCLU)
               NI = NLMQCLU(IQCLU)
C
               WRITE (18) ((-WA(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
C
            END DO
         END DO
C
      END IF
C
C=======================================================================
      IF ( SWAPTAUQ ) THEN
C
         ALLOCATE (TAUQCLU(NKMMAX,NKMMAX,NQCLU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUQCLU')
C
         IQ1 = 1
         IQ2 = NQCLU
C
         IF ( CHECKTAU ) CALL CLUTAUCHK(IPRINT,IREL,NKKR,NLM,IQ1,IQ2,
     &                                  NLMQCLU,NKMQCLU,CONC,WA,WB,WC,
     &                                  WD,IPIV,ITCPAMAX,NOQCLU,ITOQCLU,
     &                                  CPATOL,NCPA,ICPACLU,CPACHNG,
     &                                  ICPAFLAG,NSYM,DROT,IQCLUORG,
     &                                  SYMACCEPTED,SYMUNITARY,
     &                                  SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                                  TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,
     &                                  TAUQCLU)
C
         CALL CLUTAUMAT(KTAUIJ,IPRINT,IREL,NKKR,NLM,IQ1,IQ2,NLMQCLU,
     &                  NKMQCLU,CONC,WA,WB,WC,WD,IPIV,ITCPAMAX,NOQCLU,
     &                  ITOQCLU,CPATOL,NCPA,ICPACLU,CPACHNG,ICPAFLAG,
     &                  NSYM,NSYMACCEPTED,DROT,IQCLUORG,SYMACCEPTED,
     &                  SYMUNITARY,SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                  TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,TAUQCLU)
C
         DO IQCLU = 1,NQCLU
            IQ = IQ_IQCLU(IQCLU)
            N = NKMQ(IQ)
            IF ( IQCNTR.EQ.IQ ) THEN
               TAUQ(1:N,1:N,IQ) = TAUQCLU(1:N,1:N,IQCLU)
            ELSE
               TAUQ(1:N,1:N,IQ) = C0
            END IF
         END DO
C
         DEALLOCATE (TAUQCLU)
C
      ELSE
C
         IF ( CHECKTAU ) CALL CLUTAUCHK(IPRINT,IREL,NKKR,NLM,IQ1,IQ2,
     &                                  NLMQCLU,NKMQCLU,CONC,WA,WB,WC,
     &                                  WD,IPIV,ITCPAMAX,NOQCLU,ITOQCLU,
     &                                  CPATOL,NCPA,ICPACLU,CPACHNG,
     &                                  ICPAFLAG,NSYM,DROT,IQCLUORG,
     &                                  SYMACCEPTED,SYMUNITARY,
     &                                  SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                                  TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,
     &                                  TAUQ)
C
C
         WRITE (*,*) 'BBBBBBBBBBB CLUTAUMAT',I_TEMP_LAT,N_TEMP_LAT,
     &               TEMP_LAT,SWAPTAUQ
         CALL CLUTAUMAT(KTAUIJ,IPRINT,IREL,NKKR,NLM,IQ1,IQ2,NLMQCLU,
     &                  NKMQCLU,CONC,WA,WB,WC,WD,IPIV,ITCPAMAX,NOQCLU,
     &                  ITOQCLU,CPATOL,NCPA,ICPACLU,CPACHNG,ICPAFLAG,
     &                  NSYM,NSYMACCEPTED,DROT,IQCLUORG,SYMACCEPTED,
     &                  SYMUNITARY,SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                  TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,TAUQ)
C
      END IF
C=======================================================================
C
C---------------------------------------- write off-diagonal  TAU-matrix
C
      IF ( KTAUIJ.NE.0 ) THEN
C
         OPEN (UNIT=19,STATUS='SCRATCH',FORM='UNFORMATTED')
         WRITE (19) TAUQ
C
         DO JQ = 1,NQ
            IF ( NQCLU_IQ(JQ).NE.1 )
     &            CALL STOP_MESSAGE(ROUTINE,'NQCLU_IQ(JQ).NE.1')
            JQCLU = IQCLU_ISTQ(1,JQ)
            J0 = IND0QCLU(JQCLU)
            NJ = NLMQCLU(JQCLU)
            DO IQ = 1,NQ
               IF ( NQCLU_IQ(IQ).NE.1 )
     &               CALL STOP_MESSAGE(ROUTINE,'NQCLU_IQ(IQ).NE.1')
               IQCLU = IQCLU_ISTQ(1,IQ)
               I0 = IND0QCLU(IQCLU)
               NI = NLMQCLU(IQCLU)
C
               IF ( IREL.LE.2 ) THEN
                  WRITE (19) ((WA(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
                  WRITE (19) ((WB(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
               ELSE
                  WRITE (19) ((WA(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
                  WRITE (19) ((WB(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
                  WRITE (19) ((WC(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
                  WRITE (19) ((WD(I,J),I=I0+1,I0+NI),J=J0+1,J0+NJ)
               END IF
C
            END DO
         END DO
C
      END IF
C
C------------------------ transform  TAUQ  from (L,ML,MS) to (KAPPA,MUE)
C
      IF ( IREL.GE.2 ) THEN
         DO IQ = MAX(1,IQ1),MIN(IQ2,NQ)
C
            CALL CHANGEREP(NKM,NKMMAX,TAUQ(1,1,IQ),'RLM>REL',WKM2)
C
            CALL CMATCOP(NKM,NKMMAX,WKM2,TAUQ(1,1,IQ))
C
         END DO
      END IF
C
C=======================================================================
C                                 CPA
C=======================================================================
C
      IF ( NCPA.GT.0 ) THEN
C
C--------------------- transform UMAT back to (kappa,mue)-representation
         N = NKM
         M = NKMMAX
C
         DO IT = 1,NT
C
            IVT = (IT-1)*NVIBRA
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               WKM1(1:M,1:M) = UMAT_VT(1:M,1:M,IVT)
C
               CALL CHANGEREP(N,M,WKM1,'RLM>REL',UMAT_VT(1,1,IVT))
            END DO
C
         END DO
C
C----------------------------------------------------- scan scatterer IQ
C ------------------------------------------- MSSQCLU(IQCLU) -> MSSQ(IQ)
C-------------- transform MSS from (KAPPA,MUE) to (L,ML,MS) in necessary
C--------------------------------- copy MSS for equivalent cluster sites
C
         DO IQ = 1,NQ
C
            IQCLU = IQCLU_ISTQ(1,IQ)
C
            IF ( IREL.GE.2 ) THEN
C
               CALL CHANGEREP(NKM,NKMMAX,MSSQCLU(1,1,IQCLU),'RLM>REL',
     &                        MSSQ(1,1,IQ))
C
            ELSE
C
               MSSQ(1:N,1:N,IQ) = MSSQCLU(1:N,1:N,IQCLU)
C
            END IF
C
         END DO
C
C ------------------------------------------------ project on components
C ------------------------- THERMAL_VIBRA_FLUCT: transfer to local frame
C
         NKMPROJ = NKMOUT
C
         IQ = IQCNTR
C
         LOOP_IO:DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
            TAUT(:,:,IT) = 0D0
            NORM = 0D0
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
            IFT = (IT-1)*NFLUCT
            LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
               IFT = IFT + 1
C
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
               IF ( NFLUCT.EQ.1 ) THEN
C
                  TSS_FT(:,:) = TSST(:,:,IT)
C
               ELSE
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           TSST(1,1,IT),'SAS+',TSS_FT)
               END IF
C
C============================================================= IVIBRA ==
C                                             perform local displacement
               IVT = (IT-1)*NVIBRA
               LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                  IVT = IVT + 1
C
                  IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                   *NFLUCT + IFLUCT
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                  IF ( NVIBRA.EQ.1 ) THEN
C
                     TSS_VFT(:,:) = TSS_FT(:,:)
C
                  ELSE
C
                     CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,'UAUT',
     &                                 TSS_VFT)
C
                  END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                  IF ( NVIBFLU.EQ.1 ) THEN
C
                     MSS_VFT(:,:) = MSST(:,:,IT)
C
                  ELSE
C
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WKM1,MSS_VFT)
C
                  END IF
C
C---------------------------- get projection matrices DMAT_vft, DTIL_vft
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMAT_VFT,DTIL_VFT,WKM1,NKM,
     &                         MSSQ(1,1,IQ),MSS_VFT,NKMMAX)
C
C-------------------------------------------- TAUT_vft = TAUQ * DTIL_vft
C
                  CALL ZGEMM('N','N',NKM,NKM,NKM,C1,TAUQ(1,1,IQ),NKMMAX,
     &                       DTIL_VFT,NKMMAX,C0,WKM1,NKMMAX)
C
                  TAUT(:,:,IT) = TAUT(:,:,IT) + X_VFT(IVFT)*WKM1(:,:)
C
                  NORM = NORM + X_VFT(IVFT)
C
C-----------------------------------------------------------------------
C
               END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
            END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
            TAUT(:,:,IT) = TAUT(:,:,IT)/NORM
C
         END DO LOOP_IO
C================================================================= IO ==
C
C
C             IF (KMROT.NE.0)
C    &             CALL STOP_MESSAGE(ROUTINE,'ROT not implemented yet')
C
C ------------------------- write projected TAU to TAU file for types IT
C
         DO IO = 1,NOQ(IQCNTR)
            IT = ITOQ(IO,IQCNTR)
            W1(1:NKMPROJ,1:NKMPROJ) = TAUT(1:NKMPROJ,1:NKMPROJ,IT)
            IF ( .NOT.MOL ) THEN
               WRITE (9,99007) ERYD,IT,IQ,ICPAFLAG,CPACHNG
               DO J = 1,NKMOUT
                  DO I = 1,NKMOUT
                     IF ( I.EQ.J ) THEN
                        WRITE (9,99008) I,J,W1(I,J)
                     ELSE IF ( CDABS(W1(I,J)).GT.TOL ) THEN
                        WRITE (9,99008) I,J,W1(I,J)
                     END IF
                  END DO
               END DO
            END IF
         END DO
C
C---------------------------------------------------------------- NO CPA
C
      ELSE IF ( .NOT.MOL ) THEN
C
         IF ( KMROT.NE.0 ) THEN
            DO IQ = MAX(1,IQ1),MIN(IQ2,NQ)
C
               CALL ROTATE(TAUQ(1,1,IQ),'G->L',WKM2,NKM,DROTQ(1,1,IQ),
     &                     NKMMAX)
C
               DO J = 1,NKM
                  DO I = 1,NKM
                     TAUQ(I,J,IQ) = WKM2(I,J)
                  END DO
               END DO
            END DO
         END IF
C
         IF ( .NOT.IMPURITY ) THEN
            IT0 = ITCNTR
         ELSE
            IT0 = ITIMP
         END IF
C
         WRITE (9,99007) ERYD,IT0,IQCNTR
         DO J = 1,NKMOUT
            DO I = 1,NKMOUT
               IF ( CDABS(TAUQ(I,J,IQCNTR)).GT.TOL ) WRITE (9,99008) I,
     &              J,TAUQ(I,J,IQCNTR)
            END DO
         END DO
      ELSE
         DO IT = 1,NT
            IQ = IQAT(1,IT)
            WRITE (9,99007) ERYD,IT,IQ
            DO J = 1,NKMOUT
               DO I = 1,NKMOUT
                  IF ( CDABS(TAUQ(I,J,IQ)).GT.TOL ) WRITE (9,99008) I,J,
     &                 TAUQ(I,J,IQ)
               END DO
            END DO
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      DEALLOCATE (MSSQCLU,WA,WB,WC,WD,IPIV,MSSTCLU)
C
C=======================================================================
99001 FORMAT (/,1X,79('*'))
99002 FORMAT (1X,79('*'),/)
99003 FORMAT (/,10X,'TAU will be evaluated via matrix inversion')
99004 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*  *****   *       *     *   ****   *******  *******  *****  *'
     &  ,/,10X,
     &  '* *     *  *       *     *  *    *     *     *        *    * *'
     &  ,/,10X,
     &  '* *        *       *     *  *          *     *        *    * *'
     &  ,/,10X,
     &  '* *        *       *     *   ****      *     ****     *****  *'
     &  ,/,10X,
     &  '* *        *       *     *       *     *     *        * *    *'
     &  ,/,10X,
     &  '* *     *  *       *     *  *    *     *     *        *   *  *'
     &  ,/,10X,
     &  '*  *****   ******   *****    ****      *     *******  *    * *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99005 FORMAT (10X,A,9X,A,I10,:,7X,'(',I7,')')
99006 FORMAT (/,10X,'cluster configuration stored in data-file ',
     &        ' rasmol_cluster.pdb',/,10X,
     &        'view via:   rasmol  -script rasmol_cluster.ras',/)
99007 FORMAT (/,80('*')/,2F21.15,' RYD   TAU FOR IT=',I5,'  IQ=',I5,:,
     &        '  CPA:',I2,F15.8)
99008 FORMAT (2I5,1P,4E22.14)
99009 FORMAT (
     &       'HEADER    cluster confguration set up by <TAU_REAL_SPACE>'
     &       ,/,'SOURCE    SPRKKR - program       ',/,
     &       'AUTHOR    H. Ebert               ',/,
     &       'REMARK    None                   ')
99010 FORMAT ('ATOM   ',I4,'           ',I4,'    ',3F8.3,'  0.00',F8.3)
      END
