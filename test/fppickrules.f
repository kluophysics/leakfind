C*==fppickrules.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE FPPICKRULES(MODE,NMMAX,TXTT,LTXTT,IMQ,IQAT,NAT,KFP_LMQ,
     &                       KLMFP,NLMFPT,NFPT,LMIFP,ISFLM,KLMSF,NSF,
     &                       LMISF,NLMFPMAX,NLMSFMAX)
C   ********************************************************************
C   *                                                                  *
C   *   the subroutine finds all non-zero terms in an expansion        *
C   *                                                                  *
C   *               A(->r) = SUM(L) A_L(r) Y_L(^r)                     *
C   *                                                                  *
C   * - the symmetry of the system has to be fixed !!!!!               *
C   *   in particular SYMACCEPTED, NTMAX, IMQ and IQAT are set         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   * - the rotation matrices for the 'empty' lattice (needed to  set  *
C   *   up the shape functions) are determined first                   *
C   * - the 'picking rules' for the shape functions are fixed for      *
C   *   all sites IQ with info stored in  KSF_LMQ                      *
C   * - the 'picking rules' for the potential expansion are fixed for  *
C   *   all sites IQ with info stored in  KFP_LMQ                      *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  the second step depends on the requested operation MODE         *
C   *                                                                  *
C   *  MODE = DIMEN    find the array sizes  NSFMAX and  NFPMAX        *
C   *         CHECK    check whether the input is compatible           *
C   *         SETUP    set up the picking rules                        *
C   *         SHAPE    set up the picking rules                        *
C   *                  print out tables ONLY for shape functions       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_CALCMODE,ONLY:IREL,KMROT,NONMAG,MOL,RELAX_CLU
      USE MOD_SYMMETRY,ONLY:NSYM,SYMACCEPTED,MROTR,NSYMMAX
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX,NTHOST,NT
      USE MOD_SITES,ONLY:NQ,QBAS,IQBOT,IQTOP,NQMAX,QMVEC,QMPHI,QMTET,
     &    NQCLU,QBAS_QCLU,QBAS0_QCLU,DQBAS_QCLU,NQHOST,IQ_QCLU,NQ_L,NQ_R
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,BRAVAIS,SYSTEM_DIMENSION,ABAS,
     &    ABAS_L,ABAS_R,ABAS_I,ADAINV_L,ADAINV_I,ADAINV_R
      USE MOD_ANGMOM,ONLY:NL,L_LM,M_LM
      USE MOD_RMESH,ONLY:SPHERCELL
      USE MOD_FILES,ONLY:IPRINT,IOTMP,IDUMMY,RDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPPICKRULES')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      CHARACTER*5 MODE
      INTEGER NLMFPMAX,NLMSFMAX,NMMAX
      INTEGER IMQ(NQ),IQAT(NQMAX,NTMAX),ISFLM(NLMSFMAX,NMMAX),
     &        KFP_LMQ(NLMFPMAX,NQMAX),KLMFP(NLMFPMAX,NTMAX),
     &        KLMSF(NLMSFMAX,NMMAX),LMIFP(NLMFPMAX,NTMAX),
     &        LMISF(NLMSFMAX,NMMAX),LTXTT(NTMAX),NAT(NTMAX),NFPT(NTMAX),
     &        NLMFPT(NTMAX),NSF(NMMAX)
      CHARACTER*8 TXTT(NTMAX)
C
C Local variables
C
      REAL*8 CLURAD_SYM,DCELL_SYM,DCLU_SYM,DINV(:,:),DROT_RLM(:,:),
     &       DROT_RLM_TAB(:,:,:),MROTK_TMP(3,3,NSYMMAX),
     &       MROTR_TMP(3,3,NSYMMAX),QMPHI_SYM(:),QMTET_SYM(:),
     &       QMVEC_TMP(3),RMAT(:,:),RQCLU_SYM(:,:),RVEC(3),
     &       SYMEULANG(3,NSYMMAX),SYMTVEC(3,NSYMMAX),W1(:,:)
      REAL*8 DNRM2
      LOGICAL DONEQ(:),EXTEND_TABLES,KDONE_SYM(:),SYMACCEPTED_Q(:,:),
     &        SYMACCEPTED_SF(NSYMMAX),SYMACCEPTED_SFQ(:,:),
     &        SYMCRYSYS(NSYMMAX),SYMCRYSYS_TMP(NSYMMAX),
     &        SYMUNITARY_Q(:,:),SYMUNITARY_SFQ(:,:),
     &        SYMUNITARY_TMP(NSYMMAX)
      INTEGER I,I1,I2,I3,IA,IANG,IA_ERR,IEQU,IEQUORG,IFLAG,IFP,IM,IMBOT,
     &        IMTOP,IND0EQU(:),INON0(:),IPRINT_TMP,IQ,IQ1,IQAT_TMP(:,:),
     &        IQCLU,IQCLU_ATCLU0(:,:),IQCLU_ATCLU0_SYM(:,:),
     &        IQCLU_EQCLU(:,:),IQCLU_SYM,IQCNTR_SYM,IQEQ_TMP(:,:),
     &        IQLOOP,IQORG,IQORGQP_TMP(:,:),IQPSYMQ_TMP(:,:),
     &        IQREPMSYM_TMP(:),IQREPQ_TMP(:),IQ_IQCLU_SYM(:),IREL_TMP,
     &        IRMAT,ISF,ISYM,ISYMGENQ_TMP(:),IT,IT0_TMP(:),ITCLU,
     &        ITCLU_OQCLU0(:,:),ITOQ_TMP(:,:),IWEDGEROT_TMP(NSYMMAX),J,
     &        JEQU,JQ,JQCLU,JRMAT,KAUX_LMQ(:),KFP_LMQ_LOC(:,:),
     &        KKFP_LMQ(:,:),KKSF_LMQ(:,:),KMROT_TMP,KSF_LMQ(:,:),L,LM,
     &        LMTAB1(:),LMTAB2(:),LMTAB3(:),M,N1,N2,N3,NAT_SYM(:),
     &        NAT_TMP(:),NCHECK,NFPLIM,NFP_Q(:),NFP_Q_LOC(:),NKAUX,
     &        NKORG,NLMFP_Q(:),NLMFP_TMP,NLMSF_Q(:),NLMSF_TMP,NLM_TMP,
     &        NL_TMP,NMSYM_TMP
      INTEGER NOQ_TMP(:),NQCLUMAX,NQCLUMAX_SYM,NQCLU_EQCLU(:),
     &        NQCLU_I_SYM,NQCLU_L_SYM,NQCLU_R_SYM,NQCLU_SYM,NQEQMAX,
     &        NQEQ_TMP(:),NRMAT,NRMATMAX,NSFLIM,NSFTSYMQ_TMP(:,:,:),
     &        NSF_Q(:),NSHLCLU_SYM,NSYMACCEPTED_Q(:),NSYMACCEPTED_SF,
     &        NSYMACCEPTED_SFQ(:),NSYMCRYSYS,NSYMCRYSYS_TMP,NSYMH,
     &        NT0_TMP,NTCLUMAX,NTCLUMAX_SYM,NT_TMP,NWEDGE_TMP,
     &        SYMDET_TMP(NSYMMAX)
      LOGICAL RVEC_SAME
      CHARACTER*4 SYMSYMBL(NSYMMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QMPHI_SYM,QMTET_SYM,IQCLU_ATCLU0_SYM
      ALLOCATABLE RQCLU_SYM,IQ_IQCLU_SYM,KDONE_SYM,NAT_SYM
      ALLOCATABLE SYMACCEPTED_Q,SYMUNITARY_Q,NSYMACCEPTED_Q
      ALLOCATABLE SYMACCEPTED_SFQ,SYMUNITARY_SFQ,NSYMACCEPTED_SFQ
      ALLOCATABLE W1 ,NQCLU_EQCLU,IQCLU_EQCLU,DROT_RLM_TAB
      ALLOCATABLE NSFTSYMQ_TMP,NQEQ_TMP,IQEQ_TMP,IQREPMSYM_TMP
      ALLOCATABLE IQREPQ_TMP,IQCLU_ATCLU0,ITCLU_OQCLU0,DROT_RLM
      ALLOCATABLE ISYMGENQ_TMP,DONEQ,RMAT,IND0EQU,INON0,DINV
      ALLOCATABLE KSF_LMQ,KKSF_LMQ,KKFP_LMQ,KAUX_LMQ
      ALLOCATABLE LMTAB1,LMTAB2,LMTAB3
      ALLOCATABLE IQORGQP_TMP,IQAT_TMP,IQPSYMQ_TMP,ITOQ_TMP,NFP_Q
      ALLOCATABLE NLMFP_Q,NLMSF_Q,NOQ_TMP,NSF_Q,NAT_TMP,IT0_TMP
      ALLOCATABLE KFP_LMQ_LOC,NFP_Q_LOC
C
C **********************************************************************
C **********************************************************************
C                         GENERAL PART
C **********************************************************************
C **********************************************************************
C
      IMBOT = NMMAX
      IMTOP = 0
      DO IQ = IQBOT,IQTOP
         IMBOT = MIN(IMBOT,IMQ(IQ))
         IMTOP = MAX(IMTOP,IMQ(IQ))
      END DO
C
C ----------------------------------------------------------------------
C     allow to deal with full potential arrays  AND  shape functions
C     with              NLFP = 2 * NL -1        and  NLSF = 4 * NL - 3
C
      NL_TMP = 4*NL - 3
      NLM_TMP = NL_TMP**2
C
C ----------------------------------------------------------------------
C        allocate local temporary arrays according to worst case
C ----------------------------------------------------------------------
C
      ALLOCATE (KKFP_LMQ(NLMFPMAX,NQMAX),KAUX_LMQ(NLMFPMAX))
      ALLOCATE (KSF_LMQ(NLMSFMAX,NQMAX),KKSF_LMQ(NLMSFMAX,NQMAX))
C
C ----------------------------------------------------------------------
      ALLOCATE (NSFTSYMQ_TMP(3,NSYMMAX,NQMAX))
      ALLOCATE (NQEQ_TMP(NQMAX),IQEQ_TMP(NQMAX,NQMAX))
      ALLOCATE (IQREPMSYM_TMP(NQMAX),IQREPQ_TMP(NQMAX))
      ALLOCATE (ISYMGENQ_TMP(NQMAX))
C
C
      ALLOCATE (LMTAB1(NLM_TMP),LMTAB2(NLM_TMP),LMTAB3(NLM_TMP))
      ALLOCATE (W1(NLM_TMP,NLM_TMP),DINV(NLM_TMP,NLM_TMP),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: LLIM')
C
      ALLOCATE (SYMACCEPTED_Q(NSYMMAX,NQ),SYMUNITARY_Q(NSYMMAX,NQ))
      ALLOCATE (NSYMACCEPTED_Q(NQ))
      ALLOCATE (SYMACCEPTED_SFQ(NSYMMAX,NQ),SYMUNITARY_SFQ(NSYMMAX,NQ))
      ALLOCATE (NSYMACCEPTED_SFQ(NQ))
      ALLOCATE (IQORGQP_TMP(NSYMMAX,NQMAX),IQAT_TMP(NQMAX,NTMAX))
      ALLOCATE (IQPSYMQ_TMP(NSYMMAX,NQMAX),ITOQ_TMP(NTMAX,NQMAX))
      ALLOCATE (NFP_Q(NQMAX),NLMFP_Q(NQMAX),NLMSF_Q(NQMAX))
      ALLOCATE (NOQ_TMP(NQMAX),NSF_Q(NQMAX),IT0_TMP(NTMAX))
      ALLOCATE (NAT_TMP(NTMAX))
      ALLOCATE (KFP_LMQ_LOC(NLMFPMAX,NQMAX),NFP_Q_LOC(NQMAX))
C
C ----------------------------------------------------------------------
C                         INITIALIZE
C ----------------------------------------------------------------------
C
      IQORGQP_TMP(1:NSYMMAX,1:NQMAX) = 0
      IQAT_TMP(1:NQMAX,1:NTMAX) = 0
      IQPSYMQ_TMP(1:NSYMMAX,1:NQMAX) = 0
      ITOQ_TMP(1:NTMAX,1:NQMAX) = 0
      NFP_Q(1:NQMAX) = 0
      NLMFP_Q(1:NQMAX) = 0
      NLMSF_Q(1:NQMAX) = 0
      NOQ_TMP(1:NQMAX) = 0
      NSF_Q(1:NQMAX) = 0
      IT0_TMP(1:NTMAX) = 0
      NAT_TMP(1:NTMAX) = 0
C
      NLMSF_TMP = NLM_TMP
C
      NLMFP_TMP = MIN((2*NL-1)**2,NLM_TMP)
C
      NSF_Q(1:NQMAX) = 0
      NLMSF_Q(1:NQMAX) = NLMSF_TMP
      KFP_LMQ(1:NLMFPMAX,1:NQMAX) = 1
      KKFP_LMQ(1:NLMFPMAX,1:NQMAX) = 1
C
      NFP_Q(1:NQMAX) = 0
      NLMFP_Q(1:NQMAX) = NLMFP_TMP
      KSF_LMQ(1:NLMSFMAX,1:NQMAX) = 1
      KKSF_LMQ(1:NLMSFMAX,1:NQMAX) = 1
C
C
      IF ( NLMSFMAX.LT.NLMSF_TMP .OR. NLMFPMAX.LT.NLMFP_TMP ) THEN
         WRITE (6,99001) NLMSFMAX,NLMSF_TMP,NLMFPMAX,NLMFP_TMP
         CALL STOP_MESSAGE(ROUTINE,'NLMSFMAX.LT.NLMSF_TMP')
      END IF
C
C ----------------------------------------------------------------------
C
      IF ( MODE.EQ.'DIMEN' ) THEN
         WRITE (6,99016) 'fix dimensions via'
      ELSE IF ( MODE.EQ.'CHECK' ) THEN
         WRITE (6,99016) 'check'
      ELSE
         WRITE (6,99016) 'find'
      END IF
C
      IPRINT_TMP = IPRINT - 1
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
      IF ( MPI_ID.EQ.0 ) THEN
C
         ALLOCATE (DROT_RLM(NLM_TMP,NLM_TMP))
         ALLOCATE (DROT_RLM_TAB(NLM_TMP,NLM_TMP,NSYMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROT_RLM')
C
C ----------------------------------------------------------------------
C                         supply the Euler angles
C ----------------------------------------------------------------------
C
         CALL SYMINIT(IPRINT_TMP,BRAVAIS,NSYM,NSYMCRYSYS,SYMCRYSYS,
     &                SYMSYMBL,SYMEULANG)
C
C ======================================================================
C         get ALL rotation matrices  DROT  up to  NL_SF = 4*NL - 3
C         assuming a BARE LATTICE occupied with the same atom type
C         all allowed   NSYMACCEPTED_SF   symmetry operations for the
C         shape functions are indicated by    SYMACCEPTED_SF
C ======================================================================
C
         IREL_TMP = 0
         KMROT_TMP = 0
C
         NOQ_TMP(1:NQMAX) = 1
         ITOQ_TMP(1,1:NQMAX) = 1
         NT_TMP = 1
C
C-----------------------------------------------------------------------
C                            PERIODIC SYSTEM
C-----------------------------------------------------------------------
         IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
            NAT_TMP(1) = NQ
            DO IA = 1,NAT_TMP(1)
               IQAT_TMP(IA,1) = IA
            END DO
C
            CALL SYMLATTICE(0,.FALSE.,IPRINT_TMP,NWEDGE_TMP,
     &                      IWEDGEROT_TMP,NQ,QBAS,ABAS,.TRUE.,IREL_TMP,
     &                      KMROT_TMP,QMVEC,MROTR_TMP,MROTK_TMP,NSYM,
     &                      NSYMCRYSYS,SYMACCEPTED_SF,SYMUNITARY_TMP,
     &                      SYMCRYSYS,SYMTVEC,SYMSYMBL,SYMEULANG,
     &                      NOQ_TMP,ITOQ_TMP,IT0_TMP,NT0_TMP,NT_TMP,
     &                      NAT_TMP,IQAT_TMP,QMPHI,QMTET,IQORGQP_TMP,
     &                      IQPSYMQ_TMP,NSFTSYMQ_TMP,NSYMACCEPTED_SF,
     &                      SYMDET_TMP,ISYMGENQ_TMP,IQREPQ_TMP,NQEQ_TMP,
     &                      IQEQ_TMP,NMSYM_TMP,IQREPMSYM_TMP,NTMAX,
     &                      NQMAX)
C
C22222222222222222222222222222222222222222222222222222222222222222222222
            DO IQ = IQBOT,IQTOP
C
               NAT_TMP(1) = NQ
               DO IA = 1,NAT_TMP(1)
                  IQAT_TMP(IA,1) = IA
               END DO
C
               CALL SYMSITES(IQ,MOL,IPRINT_TMP,NQ,QBAS,ABAS,.TRUE.,0,
     &                       KMROT,QMVEC,MROTR,NSYM,SYMSYMBL,1,1,
     &                       IQAT_TMP,NAT_TMP,QMPHI,QMTET,NTMAX,NQMAX,
     &                       NSYMACCEPTED_SFQ(IQ),SYMACCEPTED_SFQ(1,IQ),
     &                       SYMUNITARY_SFQ(1,IQ))
C
               CALL SYMSITES(IQ,MOL,IPRINT_TMP,NQ,QBAS,ABAS,NONMAG,IREL,
     &                       KMROT,QMVEC,MROTR,NSYM,SYMSYMBL,ITBOT,
     &                       ITTOP,IQAT,NAT,QMPHI,QMTET,NTMAX,NQMAX,
     &                       NSYMACCEPTED_Q(IQ),SYMACCEPTED_Q(1,IQ),
     &                       SYMUNITARY_Q(1,IQ))
C
            END DO
C22222222222222222222222222222222222222222222222222222222222222222222222
C
         ELSE
CTEST
C            iprint=1
CTEST
C
C-----------------------------------------------------------------------
C                            CLUSTER SYSTEM
C-----------------------------------------------------------------------
C
C
            NQCLUMAX = NQMAX - NQHOST
            NTCLUMAX = NTMAX - NTHOST
C
            ALLOCATE (NQCLU_EQCLU(NQCLUMAX),
     &                IQCLU_EQCLU(NQCLUMAX,NQCLUMAX))
            ALLOCATE (IQCLU_ATCLU0(NQCLUMAX,NTCLUMAX))
            ALLOCATE (ITCLU_OQCLU0(NTCLUMAX,NQCLUMAX))
C
            DO IQCLU = 1,NQCLUMAX
               NOQ_TMP(IQCLU) = 1
               ITCLU_OQCLU0(1,IQCLU) = 1
            END DO
C
            NAT_TMP(1) = NQCLU
            DO IA = 1,NAT_TMP(1)
               IQCLU_ATCLU0(IA,1) = IA
            END DO
C
            NSYMCRYSYS_TMP = NSYM
            SYMCRYSYS_TMP(1:NSYMMAX) = SYMACCEPTED(1:NSYMMAX)
            QMVEC_TMP(1:3) = 0D0
C
            IQ1 = NQHOST + 1
C
            CALL SYMLATTICE(0,.TRUE.,IPRINT_TMP,NWEDGE_TMP,
     &                      IWEDGEROT_TMP,NQCLU,QBAS_QCLU,ABAS,.TRUE.,
     &                      IREL_TMP,KMROT_TMP,QMVEC_TMP,MROTR_TMP,
     &                      MROTK_TMP,NSYM,NSYMCRYSYS_TMP,
     &                      SYMACCEPTED_SF,SYMUNITARY_TMP,SYMCRYSYS_TMP,
     &                      SYMTVEC,SYMSYMBL,SYMEULANG,NOQ_TMP(IQ1),
     &                      ITCLU_OQCLU0,IT0_TMP,NT0_TMP,NT_TMP,NAT_TMP,
     &                      IQCLU_ATCLU0,QMPHI(IQ1),QMTET(IQ1),
     &                      IQORGQP_TMP(1,IQ1),IQPSYMQ_TMP(1,IQ1),
     &                      NSFTSYMQ_TMP(1,1,IQ1),NSYMACCEPTED_SF,
     &                      SYMDET_TMP,ISYMGENQ_TMP(IQ1),IQREPQ_TMP(IQ1)
     &                      ,NQCLU_EQCLU,IQCLU_EQCLU,NMSYM_TMP,
     &                      IQREPMSYM_TMP(IQ1),NTCLUMAX,NQCLUMAX)
C
C
            DO IQCLU = 1,NQCLU
               IQ = NQHOST + IQCLU
               NQEQ_TMP(IQ) = NQCLU_EQCLU(IQCLU)
C
               DO IEQU = 1,NQCLU_EQCLU(IQCLU)
                  IQEQ_TMP(IEQU,IQ) = NQHOST + IQCLU_EQCLU(IEQU,IQCLU)
               END DO
            END DO
C
            DO IQCLU = 1,NQCLU
               IQ = NQHOST + IQCLU
               DO ISYM = 1,NSYM
                  IQORGQP_TMP(ISYM,IQ) = NQHOST + IQORGQP_TMP(ISYM,IQ)
               END DO
            END DO
C
C22222222222222222222222222222222222222222222222222222222222222222222222
            DCELL_SYM = 0D0
            DO I1 = -1, + 1
               DO I2 = -1, + 1
                  DO I3 = -1, + 1
                     CALL RVECLCIB(I1,I2,I3,ABAS,RVEC)
                     DCELL_SYM = MAX(DCELL_SYM,DNRM2(3,RVEC,1))
                  END DO
               END DO
            END DO
C
            DO IQCLU = 1,NQCLU
C
               IQ = NQHOST + IQCLU
C
               IF ( IPRINT.GT.0 ) WRITE (6,99002) IQCLU,IQ
C
               IQCNTR_SYM = IQ_QCLU(IQCLU)
C
               DCLU_SYM = 0D0
               DO JQCLU = 1,NQCLU
                  RVEC(1:3) = QBAS0_QCLU(1:3,JQCLU)
     &                        - QBAS0_QCLU(1:3,IQCLU)
                  DCLU_SYM = MAX(DCLU_SYM,DNRM2(3,RVEC,1))
               END DO
C
C------------------------------------------------------ find NN distance
C
               CLURAD_SYM = 0D0
               NSHLCLU_SYM = 2
C
               CALL CLUSSITES(IOTMP,IPRINT-2,MOL,SYSTEM_DIMENSION,ABAS,
     &                        ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                        ADAINV_R,QBAS,CLURAD_SYM,IQCNTR_SYM,
     &                        NQCLU_SYM,NQCLU_L_SYM,NQCLU_I_SYM,
     &                        NQCLU_R_SYM,NSHLCLU_SYM,NQHOST,NQ_L,NQ_R,
     &                        NQMAX)
C
               IF ( ALLOCATED(RQCLU_SYM) )
     &              DEALLOCATE (RQCLU_SYM,IQ_IQCLU_SYM,KDONE_SYM)
C
               ALLOCATE (IQ_IQCLU_SYM(NQCLU_SYM),KDONE_SYM(NQCLU_SYM))
               ALLOCATE (RQCLU_SYM(3,NQCLU_SYM),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: RQCLU_SYM')
C
               READ (IOTMP) ((RQCLU_SYM(J,I),J=1,3),RDUMMY,IQ_IQCLU_SYM(
     &                      I),I=1,NQCLU_SYM),(IDUMMY,I=1,NSHLCLU_SYM)
               CLOSE (IOTMP)
C
C----------------------- set up cluster within radius  4.0 x NN distance
C
C              CLURAD_SYM = MAX(DCLU_SYM,DCELL_SYM)*1.0001D0
C
               CLURAD_SYM = DNRM2(3,RQCLU_SYM(1,2),1)*4.0D0
               NSHLCLU_SYM = 0
C
               CALL CLUSSITES(IOTMP,IPRINT-2,MOL,SYSTEM_DIMENSION,ABAS,
     &                        ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                        ADAINV_R,QBAS,CLURAD_SYM,IQCNTR_SYM,
     &                        NQCLU_SYM,NQCLU_L_SYM,NQCLU_I_SYM,
     &                        NQCLU_R_SYM,NSHLCLU_SYM,NQHOST,NQ_L,NQ_R,
     &                        NQMAX)
C
               IF ( ALLOCATED(RQCLU_SYM) )
     &              DEALLOCATE (RQCLU_SYM,IQ_IQCLU_SYM,KDONE_SYM)
C
               ALLOCATE (IQ_IQCLU_SYM(NQCLU_SYM),KDONE_SYM(NQCLU_SYM))
               ALLOCATE (RQCLU_SYM(3,NQCLU_SYM),STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: RQCLU_SYM')
C
               READ (IOTMP) ((RQCLU_SYM(J,I),J=1,3),RDUMMY,IQ_IQCLU_SYM(
     &                      I),I=1,NQCLU_SYM),(IDUMMY,I=1,NSHLCLU_SYM)
               CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C
               NQCLUMAX_SYM = NQCLU_SYM
               NTCLUMAX_SYM = NTMAX*NQCLU_SYM
C
               IF ( IQCLU.GT.1 ) DEALLOCATE (NAT_SYM,IQCLU_ATCLU0_SYM,
     &              QMPHI_SYM,QMTET_SYM)
C
               ALLOCATE (IQCLU_ATCLU0_SYM(NQCLUMAX_SYM,NTCLUMAX_SYM))
               ALLOCATE (QMPHI_SYM(NQCLUMAX_SYM),QMTET_SYM(NQCLUMAX_SYM)
     &                   )
               ALLOCATE (NAT_SYM(NTMAX))
C
C-----------------------------------------------------------------------
C   identify the sites in the 'real' embedded cluster with a
C   corresponding site in the auxilary symmetry cluster and
C   shift the position for the later one to get the symmetry correct
C-----------------------------------------------------------------------
C
               NCHECK = 0
               KDONE_SYM(1:NQCLU_SYM) = .FALSE.
               DO JQCLU = 1,NQCLU
                  RVEC(1:3) = QBAS0_QCLU(1:3,JQCLU)
     &                        - QBAS0_QCLU(1:3,IQCLU)
C
                  DO IQCLU_SYM = 1,NQCLU_SYM
                     IF ( .NOT.KDONE_SYM(IQCLU_SYM) ) THEN
                        IF ( RVEC_SAME(3,RVEC,RQCLU_SYM(1,IQCLU_SYM),
     &                       1D-8) ) THEN
                           RQCLU_SYM(1:3,IQCLU_SYM)
     &                        = RQCLU_SYM(1:3,IQCLU_SYM)
     &                        + DQBAS_QCLU(1:3,JQCLU)
     &                        - DQBAS_QCLU(1:3,IQCLU)
C
                           IQ_IQCLU_SYM(IQCLU_SYM) = IQ_QCLU(JQCLU)
                           KDONE_SYM(IQCLU_SYM) = .TRUE.
                           NCHECK = NCHECK + 1
                        END IF
                     END IF
                  END DO
C
               END DO
C
               IF ( CHECK ) THEN
                  WRITE (*,*) '************* IQCLU         ',IQCLU
                  WRITE (*,*) '************* IQ            ',IQ
                  WRITE (*,*) '************* DCELL_SYM     ',DCELL_SYM
                  WRITE (*,*) '************* DCLU_SYM      ',DCLU_SYM
                  WRITE (*,*) '************* CLURAD_SYM    ',CLURAD_SYM
                  WRITE (*,*) '************* SYSDIM        ',
     &                        SYSTEM_DIMENSION
                  WRITE (*,*) '************* MOL           ',MOL
                  WRITE (*,*) '************* IQCNTR_SYM    ',IQCNTR_SYM
                  WRITE (*,*) '************* NQHOST        ',NQHOST
                  WRITE (*,*) '************* NSHLCLU_SYM   ',NSHLCLU_SYM
                  WRITE (*,*) '************* NQCLU_SYM     ',NQCLU_SYM
                  WRITE (*,*) '************* NCHECK        ',NCHECK
                  WRITE (*,*) '************* NTCLUMAX_SYM  ',
     &                        NTCLUMAX_SYM
                  WRITE (*,*) '************* NQCLUMAX_SYM  ',
     &                        NQCLUMAX_SYM
                  WRITE (*,*) 
     &                   '#############################################'
                  WRITE (*,*) ' '
               END IF
               IF ( NCHECK.NE.NQCLU )
     &               CALL STOP_MESSAGE(ROUTINE,'NCHECK .NE.  NQCLU')
C-----------------------------------------------------------------------
C
               NAT_SYM(1) = NQCLU_SYM
               DO IA = 1,NAT_TMP(1)
                  IQCLU_ATCLU0(IA,1) = IA
               END DO
               QMPHI_SYM(1:NQCLUMAX_SYM) = 0D0
               QMTET_SYM(1:NQCLUMAX_SYM) = 0D0
C
               IF ( CHECK ) WRITE (*,*) 
     &                            '################################ AAA'
               CALL SYMSITES(1,.TRUE.,IPRINT_TMP,NQCLU_SYM,RQCLU_SYM,
     &                       ABAS,.TRUE.,0,KMROT_TMP,QMVEC_TMP,
     &                       MROTR_TMP,NSYM,SYMSYMBL,1,1,
     &                       IQCLU_ATCLU0_SYM,NAT_SYM,QMPHI_SYM(IQ1),
     &                       QMTET_SYM(IQ1),NTCLUMAX_SYM,NQCLUMAX_SYM,
     &                       NSYMACCEPTED_SFQ(IQ),SYMACCEPTED_SFQ(1,IQ),
     &                       SYMUNITARY_SFQ(1,IQ))
               IF ( CHECK ) WRITE (*,*) 
     &                            '################################ BBB'
C
C            NAT_TMP(1) = NQCLU
C            DO IA = 1,NAT_TMP(1)
C               IQCLU_ATCLU0(IA,1) = IA
C            END DO
C
C            CALL SYMSITES(IQCLU,.TRUE.,IPRINT_TMP,NQCLU,QBAS_QCLU,ABAS,
C     &                    .TRUE.,0,KMROT_TMP,QMVEC_TMP,MROTR_TMP,NSYM,
C     &                    SYMSYMBL,1,1,IQCLU_ATCLU0,NAT_TMP,QMPHI(IQ1),
C     &                    QMTET(IQ1),NTCLUMAX,NQCLUMAX,
C     &                    NSYMACCEPTED_SFQ(IQ),SYMACCEPTED_SFQ(1,IQ),
C     &                    SYMUNITARY_SFQ(1,IQ))
C
               IF ( IPRINT.GT.0 ) WRITE (6,99003) IQCLU,IQ
C
               DO ITCLU = 1,NT - NTHOST
                  IT = NTHOST + ITCLU
                  NAT_TMP(ITCLU) = NAT(IT)
                  DO IA = 1,NAT_TMP(ITCLU)
                     IQCLU_ATCLU0(IA,ITCLU) = IQAT(IA,IT) - NQHOST
                  END DO
               END DO
C
               CALL SYMSITES(IQCLU,.TRUE.,IPRINT_TMP,NQCLU,QBAS_QCLU,
     &                       ABAS,NONMAG,IREL,KMROT_TMP,QMVEC_TMP,
     &                       MROTR_TMP,NSYM,SYMSYMBL,1,NT-NTHOST,
     &                       IQCLU_ATCLU0,NAT_TMP,QMPHI(IQ1),QMTET(IQ1),
     &                       NTCLUMAX,NQCLUMAX,NSYMACCEPTED_Q(IQ),
     &                       SYMACCEPTED_Q(1,IQ),SYMUNITARY_Q(1,IQ))
C
C
            END DO
C22222222222222222222222222222222222222222222222222222222222222222222222
C
C
            DEALLOCATE (NQCLU_EQCLU,IQCLU_EQCLU)
C
         END IF
C-----------------------------------------------------------------------
C           preparations for PERIODIC or CLUSTER SYSTEM done
C-----------------------------------------------------------------------
C
         DEALLOCATE (ISYMGENQ_TMP,NSFTSYMQ_TMP,IQREPMSYM_TMP,IQREPQ_TMP)
C
C-----------------------------------------------------------------------
C
         ALLOCATE (DONEQ(NQMAX),IND0EQU(NQMAX))
C
         NQEQMAX = MAXVAL(NQEQ_TMP(IQBOT:IQTOP))
C
         NRMATMAX = NLMSF_TMP*NQEQMAX
         ALLOCATE (RMAT(NRMATMAX,NRMATMAX),INON0(NRMATMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) WRITE (6,*) '<FPPICKRULES> NRMATMAX=',
     &                                  NRMATMAX
C
         IND0EQU(1) = 0
         DO IEQU = 2,NQEQMAX
            IND0EQU(IEQU) = IND0EQU(IEQU-1) + NLMSF_TMP
         END DO
C
         NSYMH = NSYM/2
C
C-----------------------------------------------------------------------
C                     create matrix for inversion
C-----------------------------------------------------------------------
         CALL RINIT(NLM_TMP*NLM_TMP,DINV)
C
         I = 0
         DO L = 0,(NL_TMP-1)
            DO M = -L,L
               I = I + 1
               DINV(I,I) = (-1.0D0)**L
            END DO
         END DO
C
C ******************************************************** ISYM -- START
C                    scan ALL symmetry operations
C ******************************************************** ISYM -- START
         LOOP_SCAN_ALL_SYMMETRY_OPERATIONS:DO ISYM = 1,NSYM
            DO IQ = IQBOT,IQTOP
               DONEQ(IQ) = .FALSE.
            END DO
C
C-------------------- get rotation matrices for REAL spherical harmonics
C
            IF ( ISYM.LE.NSYMH ) THEN
               IANG = ISYM
C
               CALL ROTMAT_RYLM(NL_TMP,NLM_TMP,SYMEULANG(1,IANG),
     &                          SYMEULANG(2,IANG),SYMEULANG(3,IANG),
     &                          DROT_RLM)
C
            ELSE
               IANG = ISYM - NSYMH
C
               CALL ROTMAT_RYLM(NL_TMP,NLM_TMP,SYMEULANG(1,IANG),
     &                          SYMEULANG(2,IANG),SYMEULANG(3,IANG),W1)
C
               CALL RMATMUL(NLM_TMP,NLM_TMP,W1,DINV,DROT_RLM)
C
            END IF
C
            DROT_RLM_TAB(1:NLM_TMP,1:NLM_TMP,ISYM)
     &         = DROT_RLM(1:NLM_TMP,1:NLM_TMP)
C
C22222222222222222222222222222222222222222222222222222222222222222222222
C               use scheme 2  working with local symmetry
C22222222222222222222222222222222222222222222222222222222222222222222222
C
            DO IQ = IQBOT,IQTOP
               IF ( SYMACCEPTED_Q(ISYM,IQ) .OR. SYMACCEPTED_SFQ(ISYM,IQ)
     &              ) THEN
C
C-----------------------------------------------------------------------
C         deal with homogenous equations   (R-1) * V = 0
C       transform the coefficient matrix (R-1) to (tri)diagonal form
C-----------------------------------------------------------------------
C
                  NRMAT = NLMSF_TMP
C
                  RMAT(1:NRMAT,1:NRMAT) = DROT_RLM(1:NRMAT,1:NRMAT)
                  DO I = 1,NRMAT
                     RMAT(I,I) = RMAT(I,I) - 1D0
                  END DO
C
                  CALL FPPICKRULES_HOM_EQU(RMAT,INON0,NRMAT,NRMATMAX)
C
C-----------------------------------------------------------------------
C                    deal with shape functions
C-----------------------------------------------------------------------
C
                  IF ( SYMACCEPTED_SFQ(ISYM,IQ) ) THEN
                     DO LM = 1,NLMSF_TMP
                        IF ( INON0(LM).EQ.1 ) KKSF_LMQ(LM,IQ) = 0
                     END DO
                  END IF
C
C-----------------------------------------------------------------------
C                    deal with potential and charge expansion
C-----------------------------------------------------------------------
C
                  IF ( SYMACCEPTED_Q(ISYM,IQ) ) THEN
                     DO LM = 1,NLMFP_TMP
                        IF ( INON0(LM).EQ.1 ) KKFP_LMQ(LM,IQ) = 0
                     END DO
                  END IF
C
               END IF
            END DO
C22222222222222222222222222222222222222222222222222222222222222222222222
C
            IF ( SYMACCEPTED_SF(ISYM) ) THEN
C
C ****************************************************** IQLOOP -- START
               DO IQLOOP = IQBOT,IQTOP
                  IF ( .NOT.DONEQ(IQLOOP) ) THEN
C
C-----------------------------------------------------------------------
C         deal with homogenous equations   (R-1) * V = 0
C       transform the coefficient matrix (R-1) to (tri)diagonal form
C-----------------------------------------------------------------------
C
                     CALL RINIT(NRMATMAX*NRMATMAX,RMAT)
                     NRMAT = NLMSF_TMP*NQEQ_TMP(IQLOOP)
C
                     DO IEQU = 1,NQEQ_TMP(IQLOOP)
C
                        IQ = IQEQ_TMP(IEQU,IQLOOP)
                        IQORG = IQORGQP_TMP(ISYM,IQ)
                        DONEQ(IQ) = .TRUE.
C
                        DO JEQU = 1,NQEQ_TMP(IQLOOP)
                           JQ = IQEQ_TMP(JEQU,IQLOOP)
                           IF ( JQ.EQ.IQORG ) THEN
                              IEQUORG = JEQU
                              GOTO 2
                           END IF
                        END DO
                        CALL STOP_MESSAGE(ROUTINE,'IEQUORG not found')
C
 2                      CONTINUE
                        JRMAT = IND0EQU(IEQUORG)
                        DO J = 1,NLM_TMP
                           JRMAT = JRMAT + 1
                           IRMAT = IND0EQU(IEQU)
                           DO I = 1,NLM_TMP
                              IRMAT = IRMAT + 1
                              RMAT(IRMAT,JRMAT) = DROT_RLM(I,J)
                           END DO
                        END DO
C
                     END DO
C
                     DO I = 1,NRMAT
                        RMAT(I,I) = RMAT(I,I) - 1D0
                     END DO
C
                     CALL FPPICKRULES_HOM_EQU(RMAT,INON0,NRMAT,NRMATMAX)
C
C-----------------------------------------------------------------------
C                    deal with shape functions
C-----------------------------------------------------------------------
C
                     DO IEQU = 1,NQEQ_TMP(IQLOOP)
                        IQ = IQEQ_TMP(IEQU,IQLOOP)
C
                        DO LM = 1,NLMSF_TMP
                           I = IND0EQU(IEQU) + LM
                           IF ( INON0(I).EQ.1 ) KSF_LMQ(LM,IQ) = 0
                        END DO
C
                     END DO
C
C-----------------------------------------------------------------------
C                    deal with potential and charge expansion
C-----------------------------------------------------------------------
C
                     IF ( SYMACCEPTED(ISYM) ) THEN
                        DO IEQU = 1,NQEQ_TMP(IQLOOP)
                           IQ = IQEQ_TMP(IEQU,IQLOOP)
C
                           DO LM = 1,NLMFP_TMP
                              I = IND0EQU(IEQU) + LM
                              IF ( INON0(I).EQ.1 ) KFP_LMQ(LM,IQ) = 0
                           END DO
C
                        END DO
                     END IF
C
                  END IF
               END DO
C ******************************************************** IQLOOP -- END
C
            END IF
         END DO LOOP_SCAN_ALL_SYMMETRY_OPERATIONS
C ********************************************************** ISYM -- END
C
C
C=======================================================================
C                  check equivalency of lattice sites
C=======================================================================
C
C ******************************************************** ISYM -- START
         DO ISYM = 1,NSYM
C
            DO IQ = IQBOT,IQTOP
               IF ( SYMACCEPTED(ISYM) ) THEN
                  IQORG = IQORGQP_TMP(ISYM,IQ)
C
                  CALL FPPICKRULES_ROT_KLM(NKAUX,KAUX_LMQ,NKORG,
     &               KFP_LMQ(1,IQORG),DROT_RLM_TAB(1,1,ISYM),NLMFP_TMP,
     &               NLM_TMP)
                  IF ( IQ.NE.IQORG ) THEN
C
                     N1 = 0
                     DO LM = 1,NLMFP_TMP
                        IF ( KFP_LMQ(LM,IQ).NE.0 ) N1 = N1 + 1
                     END DO
C
                     N1 = 0
                     N2 = 0
                     N3 = 0
                     IFLAG = 0
                     DO LM = 1,NLMFP_TMP
                        IF ( KFP_LMQ(LM,IQ).NE.KAUX_LMQ(LM) ) THEN
                           IF ( KFP_LMQ(LM,IQ).EQ.1 ) THEN
                              IFLAG = 2
                           ELSE
                              IFLAG = 1
                           END IF
                        END IF
                        IF ( KFP_LMQ(LM,IQ).NE.0 ) THEN
                           N1 = N1 + 1
                           LMTAB1(N1) = LM
                        END IF
                        IF ( KAUX_LMQ(LM).NE.0 ) THEN
                           N2 = N2 + 1
                           LMTAB2(N2) = LM
                        END IF
                        IF ( KFP_LMQ(LM,IQORG).NE.0 ) THEN
                           N3 = N3 + 1
                           LMTAB3(N3) = LM
                        END IF
                     END DO
C
                     IF ( IPRINT_TMP.GT.0 .OR. IFLAG.EQ.2 ) THEN
                        WRITE (6,99009) ISYM,IQORG,IQ
                        WRITE (6,99010) 'LM(ORG)  ',(LMTAB3(I),I=1,N3)
                        WRITE (6,99010) 'LM(ROT)  ',(LMTAB2(I),I=1,N2)
                        WRITE (6,99010) 'LM(IQ)   ',(LMTAB1(I),I=1,N1)
                        IF ( IFLAG.EQ.1 ) WRITE (6,99013)
                        IF ( IFLAG.EQ.2 ) WRITE (6,99014)
                     END IF
C
                  END IF
C
               END IF
            END DO
C
         END DO
C ********************************************************** ISYM -- END
C
C ======================================================================
C                        set  NSF_Q  and  NFP_Q
C ======================================================================
C
         DO IQ = IQBOT,IQTOP
C
            IF ( IPRINT.GT.0 ) WRITE (6,99004) IQ
C
C-----------------------------------------------------------------------
C                    deal with shape functions
C-----------------------------------------------------------------------
C
            DO LM = 1,NLMSF_TMP
               IF ( KSF_LMQ(LM,IQ).EQ.1 ) THEN
                  NSF_Q(IQ) = NSF_Q(IQ) + 1
                  IF ( IPRINT.GT.0 ) WRITE (6,'(10X,A,4I4)')
     &                  'shape functions',NSF_Q(IQ),LM,L_LM(LM),M_LM(LM)
               END IF
C
               IF ( KSF_LMQ(LM,IQ).NE.KKSF_LMQ(LM,IQ) ) THEN
C
C------------------------------------------------------ KSF=0 and KKSF=1
                  IF ( KSF_LMQ(LM,IQ).EQ.0 ) THEN
C
                     IF ( IPRINT.GT.0 ) WRITE (6,99007) IQ,LM,'KSF_LMQ',
     &                    KSF_LMQ(LM,IQ),KKSF_LMQ(LM,IQ)
C------------------------------------------------------ KSF=1 and KKSF=0
                  ELSE
C
                     IF ( IPRINT.GT.0 ) WRITE (6,99008) IQ,LM,'KSF_LMQ',
     &                    KSF_LMQ(LM,IQ),KKSF_LMQ(LM,IQ)
C
                     NSF_Q(IQ) = NSF_Q(IQ) - 1
                     KSF_LMQ(LM,IQ) = 0
                  END IF
C
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
C                    deal with potential and charge expansion
C-----------------------------------------------------------------------
C                   in doubt: accept contribution !
C
            IF ( IPRINT.GT.0 ) WRITE (6,*) '  '
C
            DO LM = 1,NLMFP_TMP
               IF ( KFP_LMQ(LM,IQ).EQ.1 ) THEN
                  NFP_Q(IQ) = NFP_Q(IQ) + 1
                  IF ( IPRINT.GT.0 ) WRITE (6,'(10X,A,4I4)')
     &                  'full potential ',NFP_Q(IQ),LM,L_LM(LM),M_LM(LM)
               END IF
C
               IF ( KFP_LMQ(LM,IQ).NE.KKFP_LMQ(LM,IQ) ) THEN
C
C------------------------------------------------------ KFP=0 and KKFP=1
                  IF ( KFP_LMQ(LM,IQ).EQ.0 ) THEN
C
Ccc                     WRITE (6,99007) IQ,LM,'KFP_LMQ',KFP_LMQ(LM,IQ),
Ccc     &                               KKFP_LMQ(LM,IQ)
C
C accept scheme 2
                     NFP_Q(IQ) = NFP_Q(IQ) + 1
                     KFP_LMQ(LM,IQ) = 1
                     WRITE (6,*) '##################   KFP=0 and KKFP=1'
Ccc                     STOP
C
C------------------------------------------------------ KFP=1 and KKFP=0
                  ELSE
C
                     IF ( IPRINT.GT.0 ) WRITE (6,99008) IQ,LM,'KFP_LMQ',
     &                    KFP_LMQ(LM,IQ),KKFP_LMQ(LM,IQ)
C
C accept scheme 1 - do not change setting for KFP_LMQ
C
                     WRITE (6,*) '##################   KFP=1 and KKFP=0'
                     KKFP_LMQ(LM,IQ) = 1
C
                  END IF
C
               END IF
C
            END DO
         END DO
C
C ======================================================================
C
         DEALLOCATE (RMAT,W1,DROT_RLM,DINV,INON0,IND0EQU,DONEQ)
C
      END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
C
         CALL DRV_MPI_BCAST_I(0,IQORGQP_TMP(1,1),NSYMMAX*NQMAX)
         CALL DRV_MPI_BCAST_I(0,IQAT_TMP(1,1),NQMAX*NTMAX)
         CALL DRV_MPI_BCAST_I(0,IQPSYMQ_TMP(1,1),NSYMMAX*NQMAX)
         CALL DRV_MPI_BCAST_I(0,ITOQ_TMP(1,1),NTMAX*NQMAX)
         CALL DRV_MPI_BCAST_I(0,NFP_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NLMFP_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NLMSF_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NOQ_TMP(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NSF_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,IT0_TMP(1),NTMAX)
         CALL DRV_MPI_BCAST_I(0,NAT_TMP(1),NTMAX)
C
         CALL DRV_MPI_BCAST_I(0,NSF_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NLMSF_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,KFP_LMQ(1,1),NLMFPMAX*NQMAX)
         CALL DRV_MPI_BCAST_I(0,KKFP_LMQ(1,1),NLMFPMAX*NQMAX)
C
         CALL DRV_MPI_BCAST_I(0,NFP_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,NLMFP_Q(1),NQMAX)
         CALL DRV_MPI_BCAST_I(0,KSF_LMQ(1,1),NLMSFMAX*NQMAX)
         CALL DRV_MPI_BCAST_I(0,KKSF_LMQ(1,1),NLMSFMAX*NQMAX)
C
      END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C ======================================================================
C       transfer the picking rules to the local frame if necessary
C
C       in this case the quantities *_LOC refer to the local frame
C       with the magnetic moment rotated and therefore correspond to
C       the type dependent situation
C ======================================================================
C
      IF ( KMROT.NE.0 .AND. IREL.EQ.3 ) THEN
C
         CALL FPPICKRULES_LOCAL(NFP_Q,KFP_LMQ,NFP_Q_LOC,KFP_LMQ_LOC,
     &                          NLMFPMAX,NSF_Q,KSF_LMQ,NLMSFMAX,NQMAX)
C
      ELSE
C
         NFP_Q_LOC(:) = NFP_Q(:)
         KFP_LMQ_LOC(:,:) = KFP_LMQ(:,:)
C
      END IF
C
C **********************************************************************
C **********************************************************************
C                       MODE SPECIFIC PART
C **********************************************************************
C **********************************************************************
C
      IF ( IPRINT.GT.0 ) WRITE (6,99005) MODE
C
C
C ======================================================================
C             find DIMENSION of arrays from the  PICKING RULES    DIMEN
C ======================================================================
      IF ( MODE.EQ.'DIMEN' ) THEN
C
         NSFLIM = 0
         NFPLIM = 0
         DO IQ = IQBOT,IQTOP
C
            ISF = 0
            DO LM = 1,NLMSF_Q(IQ)
               IF ( KSF_LMQ(LM,IQ).NE.0 ) ISF = ISF + 1
            END DO
            NSFLIM = MAX(NSFLIM,ISF)
C
            IFP = 0
            DO LM = 1,NLMFPMAX
               IF ( KFP_LMQ_LOC(LM,IQ).NE.0 ) IFP = IFP + 1
            END DO
            NFPLIM = MAX(NFPLIM,IFP)
C
         END DO
C
C------------------------------------------ export results via NSF, NFPT
         NSF(1) = NSFLIM
         NFPT(1) = NFPLIM
C
         DEALLOCATE (NQEQ_TMP,IQEQ_TMP)
C
         RETURN
C
C ======================================================================
C             CHECK present settings for   PICKING RULES          CHECK
C ======================================================================
      ELSE IF ( MODE.EQ.'CHECK' ) THEN
C
         DO IT = ITBOT,ITTOP
C
            IFLAG = 0
C------------------------------ deal with potential and charge expansion
            IQ = IQAT(1,IT)
C
            IF ( NFPT(IT).GT.NFP_Q_LOC(IQ) .OR. NLMFPT(IT)
     &           .GT.NLMFP_Q(IQ) ) THEN
               WRITE (6,99015) 
     &                     'input file supplied more data than expected'
               WRITE (6,99011) 'IT    ',IT,' IQ    ',IQ
               WRITE (6,99011) 'NFPT  ',NFPT(IT),' NFP_Q_LOC  ',
     &                         NFP_Q_LOC(IQ)
               WRITE (6,99011) 'NLMFPT',NLMFPT(IT),' NLMFP_Q',
     &                         NLMFP_Q(IQ)
               IFLAG = 1
            END IF
C
C ------------------------ overwrite input tables with symmetry settings
            EXTEND_TABLES = .FALSE.
            DO LM = 1,NLMFP_Q(IQ)
               IF ( KLMFP(LM,IT).NE.KFP_LMQ_LOC(LM,IQ) ) THEN
                  IF ( KLMFP(LM,IT).GT.KFP_LMQ_LOC(LM,IQ) ) THEN
                     IFLAG = 2
                     WRITE (6,99015) 'input data inconsistent'
                     WRITE (6,99012) 'IT      IQ     LM',IT,IQ,LM
                     WRITE (6,99012) 'KLMFP   KFP_LMQ_LOC  ',
     &                               KLMFP(LM,IT),KFP_LMQ_LOC(LM,IQ)
                  END IF
                  EXTEND_TABLES = .TRUE.
               END IF
            END DO
C
            IF ( EXTEND_TABLES .AND. IFLAG.EQ.0 ) THEN
               WRITE (6,99015) 'extending tables for KLMFP and LMIFP'//
     &                         ' -- original settings:'
               WRITE (6,99019) 'atom type',IT,TXTT(IT)(1:LTXTT(IT))
               WRITE (6,99019) 'NFP ',NFPT(IT),'mesh IM ',IMQ(IQ)
               WRITE (6,99020) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                         M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
C
               KLMFP(1:NLMFPMAX,IT) = KFP_LMQ_LOC(1:NLMFPMAX,IQ)
               NLMFPT(IT) = NLMFP_Q(IQ)
               NFPT(IT) = NFP_Q_LOC(IQ)
               IFP = 0
               DO LM = 1,NLMFPMAX
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     IFP = IFP + 1
                     LMIFP(IFP,IT) = LM
                  END IF
               END DO
            END IF
C
            IFP = 0
            DO LM = 1,NLMFPMAX
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  IFP = IFP + 1
                  IF ( LMIFP(IFP,IT).NE.LM ) THEN
                     WRITE (6,99015) 'input data inconsistent'
                     WRITE (6,99012) 'IT      IFP      ',IT,IFP
                     WRITE (6,99012) 'LMIFP   LM       ',LMIFP(IFP,IT),
     &                               LM
                     IFLAG = 1
                  END IF
               END IF
            END DO
C
            IF ( .NOT.SPHERCELL ) THEN
C--------------------------------------------- deal with shape functions
C
               IQ = IQAT(1,IT)
               IM = IMQ(IQ)
C
               IF ( NSF(IM).NE.NSF_Q(IQ) ) THEN
                  IF ( NSF(IM).GT.NSF_Q(IQ) .OR. IPRINT.GT.0 ) THEN
                     WRITE (6,99015) 'SHAPE FUNCTION data inconsistent:'
     &                               //'  NSF(IM) .NE. NSF_Q(IQ)'
                     WRITE (6,99012) 'IM    IQ      ',IM,IQ
                     WRITE (6,99012) 'NSF   NSF_Q   ',NSF(IM),NSF_Q(IQ)
                     IFLAG = 1
                  END IF
               END IF
               ISF = 0
               DO LM = 1,NLMSF_Q(IQ)
                  IF ( KLMSF(LM,IM).NE.KSF_LMQ(LM,IQ) ) THEN
                     IF ( (KLMSF(LM,IM).EQ.1 .AND. KSF_LMQ(LM,IQ).EQ.0)
     &                    .OR. IPRINT.GT.0 ) THEN
                        WRITE (6,99015) 
     &                                'SHAPE FUNCTION data inconsistent'
                        WRITE (6,99012) 'IM      IQ     LM',IM,IQ,LM
                        WRITE (6,99012) 'KLMSF   KSF_LMQ  ',KLMSF(LM,IM)
     &                                  ,KSF_LMQ(LM,IQ)
                        IFLAG = 1
                     END IF
                  END IF
                  IF ( KLMSF(LM,IM).NE.0 ) THEN
                     ISF = ISF + 1
                     IF ( LMISF(ISF,IM).NE.LM ) THEN
                        WRITE (6,99015) 
     &                          'SHAPE FUNCTION input data inconsistent'
                        WRITE (6,99012) 'IM      ISF      ',IM,ISF
                        WRITE (6,99012) 'LMISF   LM       ',
     &                                  LMISF(ISF,IM),LM
                        IFLAG = 1
                     END IF
                  END IF
               END DO
            END IF
C
         END DO
C
         DEALLOCATE (NQEQ_TMP,IQEQ_TMP)
C
         WRITE (6,99017)
         DO IM = IMBOT,IMTOP
            WRITE (6,*) ' '
            WRITE (6,99018) 'mesh',IM
            WRITE (6,99019) 'NSF ',NSF(IM)
            WRITE (6,99020) (LMISF(ISF,IM),L_LM(LMISF(ISF,IM)),
     &                      M_LM(LMISF(ISF,IM)),ISF=1,NSF(IM))
         END DO
C
         WRITE (6,99021)
         DO IT = ITBOT,ITTOP
            IQ = IQAT(1,IT)
            IM = IMQ(IQ)
            WRITE (6,*) ' '
            WRITE (6,99019) 'atom type',IT,TXTT(IT)(1:LTXTT(IT))
            WRITE (6,99019) 'NFP ',NFPT(IT),'mesh IM ',IM
            WRITE (6,99020) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                      M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
         END DO
C
         IFLAG = 0
         IF ( IFLAG.EQ.0 ) THEN
            WRITE (6,99023)
            WRITE (6,99006) MODE
            RETURN
         ELSE
            CALL STOP_MESSAGE(ROUTINE,'IFLAG.NE.0')
         END IF
C
C ======================================================================
C                       SETUP   of   PICKING RULES                SETUP
C ======================================================================
      ELSE IF ( MODE.EQ.'SETUP' .OR. MODE.EQ.'SHAPE' ) THEN
C
         DO LM = 1,NLMFPMAX
            DO IT = ITBOT,ITTOP
               KLMFP(LM,IT) = 0
            END DO
         END DO
         IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' .OR. RELAX_CLU )
     &        THEN
            DO LM = 1,NLMSFMAX
               DO IM = 1,NMMAX
                  KLMSF(LM,IM) = 0
               END DO
            END DO
         END IF
C
         DO IT = ITBOT,ITTOP
C------------------------------ deal with potential and charge expansion
            IQ = IQAT(1,IT)
            NFPT(IT) = NFP_Q_LOC(IQ)
            NLMFPT(IT) = NLMFP_Q(IQ)
            IFP = 0
            DO LM = 1,NLMFPMAX
               KLMFP(LM,IT) = KFP_LMQ_LOC(LM,IQ)
               IF ( KLMFP(LM,IT).NE.0 ) THEN
                  IFP = IFP + 1
                  LMIFP(IFP,IT) = LM
               END IF
            END DO
C
C--------------------------------------------- deal with shape functions
C       no cluster-relaxation: use shape functions from host calculation
C
            IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' .OR. 
     &           RELAX_CLU ) THEN
C
               IQ = IQAT(1,IT)
               IM = IMQ(IQ)
               IF ( .NOT.SPHERCELL ) THEN
                  NSF(IM) = NSF_Q(IQ)
                  ISF = 0
                  DO LM = 1,NLMSF_Q(IQ)
                     KLMSF(LM,IM) = KSF_LMQ(LM,IQ)
                     IF ( KLMSF(LM,IM).NE.0 ) THEN
                        ISF = ISF + 1
                        LMISF(ISF,IM) = LM
                        ISFLM(LM,IM) = ISF
                     END IF
                  END DO
               ELSE
                  NSF(IM) = 1
                  LMISF(1,IM) = 1
                  ISFLM(1,IM) = 1
                  KLMSF(1,IM) = 1
               END IF
            END IF
C
         END DO
C
         IF ( IPRINT.GE.0 ) THEN
            WRITE (6,99017)
            DO IM = IMBOT,IMTOP
               WRITE (6,*) ' '
               WRITE (6,99018) 'mesh',IM
               WRITE (6,99019) 'NSF ',NSF(IM)
               WRITE (6,99020) (LMISF(ISF,IM),L_LM(LMISF(ISF,IM)),
     &                         M_LM(LMISF(ISF,IM)),ISF=1,NSF(IM))
            END DO
         END IF
C
         IF ( MODE.EQ.'SETUP' ) THEN
C
            WRITE (6,99021)
            DO IT = ITBOT,ITTOP
               IQ = IQAT(1,IT)
               IM = IMQ(IQ)
               WRITE (6,*) ' '
               WRITE (6,99019) 'atom type',IT,TXTT(IT)(1:LTXTT(IT))
               WRITE (6,99019) 'NFP ',NFPT(IT),'mesh IM ',IM
               WRITE (6,99020) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                         M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
            END DO
C
         END IF
C
         WRITE (6,99006) MODE
C
         DEALLOCATE (NQEQ_TMP,IQEQ_TMP)
C
CTEST
C         IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) STOP
CTEST
C ======================================================================
      ELSE
         WRITE (6,99022) MODE
         STOP
      END IF
C ======================================================================
C
      IF ( MODE.EQ.'CHECK' ) STOP
C
C ======================================================================
99001 FORMAT (10X,'<FPPICKINGRULES> ',/,10X,'NLMSFMAX  = ',I5,10X,
     &        'NLMSF_TMP   = ',I5,/,10X,'NLMFPMAX  = ',I5,10X,
     &        'NLMFP_TMP   = ',I5,/)
99002 FORMAT (//,10X,'<FPPICKRULES>:  finding symmetry of ',
     &        'shape functions',/,10X,'for  IQCLU =',I5,'      IQ = ',
     &        I5,/)
99003 FORMAT (//,10X,'<FPPICKRULES>:  finding symmetry of',
     &        ' full potential terms',/,10X,'for  IQCLU =',I5,
     &        '      IQ = ',I5,/)
99004 FORMAT (//,10X,'set  NSF_Q  and  NFP_Q  for IQ = ',I4,//,28X,
     &        'I  LM   L   M ')
99005 FORMAT (/,10X,'<FPPICKRULES>:starting specific part for MODE = ',
     &        A,/)
99006 FORMAT (/,10X,'<FPPICKRULES> done for MODE = ',A,/)
C99007 FORMAT ('TROUBLE ',/,'TROUBLE IQ=',I3,'  LM=',I3,/,'TROUBLE ',A,
99007 FORMAT ('INFO ',/,'INFO IQ=',I3,'  LM=',I3,/,'INFO ',A,
     &        '(LM,IQ)  for scheme 1 vs. 2:',2I3)
99008 FORMAT ('INFO ',/,'INFO IQ=',I3,'  LM=',I3,/,'INFO ',A,
     &        '(LM,IQ)  for scheme 1 vs. 2:',2I3)
99009 FORMAT (/,5X,'ISYM =',I3,20X,' IQORG =',I3,'  -->  IQ =',I3)
99010 FORMAT (5X,A,(20I3))
99011 FORMAT (10X,A,I5,5X,A,10I5)
99012 FORMAT (10X,A,10I5)
99013 FORMAT (5X,55('-'))
99014 FORMAT (' ##### TROUBLE ',65('#'))
99015 FORMAT (/,' ##### TROUBLE in <FPPICKRULES> ',48('#'),:,/,10X,A)
99016 FORMAT (2(/,1X,79('*')),/,32X,'<FPPICKRULES>',/,18X,A,
     &        ' the full potential picking rules',2(/,1X,79('*')),/)
99017 FORMAT (/,10X,'shape functions:')
99018 FORMAT (10X,A,I4,:,4X,'JRCUT   ',I5,5X,'R =',F12.8,5X,F12.8)
99019 FORMAT (10X,A,I4,4X,A,I5)
99020 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99021 FORMAT (/,10X,'full potential parameters:')
99022 FORMAT (/,'STOP in <FPPICKRULES>   MODE=',A,'  not allowed')
99023 FORMAT (/,10X,'settings for picking rules from input    OK',/)
      END
C*==fpdimensions.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE FPDIMENSIONS(IQAT,NAT,NTMAX,NFPLIM,NSFLIM)
C   ********************************************************************
C   *                                                                  *
C   *   driving routine to call FPPICKRULES to get array sizes         *
C   *   for full potential calculations                                *
C   *                                                                  *
C   *   NOTE: IQAT and NTMAX are already fixed                         *
C   *         NFPLIM and NSFLIM are still to be determined             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NLMSFMAX
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NLMFPMAX
      USE MOD_SITES,ONLY:NQMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPDIMENSIONS')
C
C Dummy arguments
C
      INTEGER NFPLIM,NSFLIM,NTMAX
      INTEGER IQAT(NQMAX,NTMAX),NAT(NTMAX)
C
C Local variables
C
      INTEGER IA_ERR,IM_QAUX(:),ISF_LMMAUX(:,:),IT,KFP_LMQ(:,:),
     &        KFP_LMTAUX(:,:),KSF_LMTAUX(:,:),LM_FPTAUX(:,:),
     &        LM_SFMAUX(:,:),LTXT_TAUX(:),NFP_TAUX(:),NLMFP_TAUX(:),
     &        NMMAX,NSFMAX,NSF_MAUX(:)
      CHARACTER*8 TXT_TAUX(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TXT_TAUX,LTXT_TAUX,IM_QAUX,ISF_LMMAUX,KSF_LMTAUX
      ALLOCATABLE KFP_LMTAUX,LM_FPTAUX,LM_SFMAUX
      ALLOCATABLE KFP_LMQ,NLMFP_TAUX,NFP_TAUX,NSF_MAUX
C
      ALLOCATE (TXT_TAUX(NTMAX),LTXT_TAUX(NTMAX))
      ALLOCATE (IM_QAUX(NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IM_QAUX')
C
C ----------------------------------------------------------------------
C     allow to deal with full potential arrays  AND  shape functions
C     with              NLFP = 2 * NL -1        and  NLSF = 4 * NL - 3
C
C ----------------------------------------------------------------------
      NSFMAX = NLMSFMAX
      NMMAX = 1
C-------------------------------------------------------- dummy settings
      DO IT = ITBOT,ITTOP
         TXT_TAUX(IT) = 'XXXX'
         LTXT_TAUX(IT) = 4
      END DO
      IM_QAUX(1:NQMAX) = 1
C ----------------------------------------------------------------------
C
      ALLOCATE (NFP_TAUX(NTMAX),NSF_MAUX(NMMAX),KFP_LMQ(NLMFPMAX,NQMAX))
      ALLOCATE (LM_SFMAUX(NSFMAX,NMMAX),NLMFP_TAUX(NTMAX))
      ALLOCATE (KFP_LMTAUX(NLMFPMAX,NTMAX),LM_FPTAUX(NLMFPMAX,NTMAX))
      ALLOCATE (ISF_LMMAUX(NLMSFMAX,NMMAX))
      ALLOCATE (KSF_LMTAUX(NLMSFMAX,NMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISF_LMMAUX')
C
      NSF_MAUX(:) = 0
      NFP_TAUX(:) = 0
      NLMFP_TAUX(:) = 0
      KFP_LMTAUX(:,:) = 0
      LM_FPTAUX(:,:) = 0
      KSF_LMTAUX(:,:) = 0
      LM_SFMAUX(:,:) = 0
C
      CALL FPPICKRULES('DIMEN',NMMAX,TXT_TAUX,LTXT_TAUX,IM_QAUX,IQAT,
     &                 NAT,KFP_LMQ,KFP_LMTAUX,NLMFP_TAUX,NFP_TAUX,
     &                 LM_FPTAUX,ISF_LMMAUX,KSF_LMTAUX,NSF_MAUX,
     &                 LM_SFMAUX,NLMFPMAX,NLMSFMAX)
C
      NSFLIM = NSF_MAUX(1)
      NFPLIM = NFP_TAUX(1)
C
      WRITE (6,99001) NFPLIM,NSFLIM
C
      DEALLOCATE (TXT_TAUX,LTXT_TAUX,IM_QAUX,ISF_LMMAUX,KSF_LMTAUX)
      DEALLOCATE (KFP_LMTAUX,LM_FPTAUX,LM_SFMAUX)
      DEALLOCATE (NLMFP_TAUX,NFP_TAUX,NSF_MAUX)
C
99001 FORMAT (/,1X,79('*'),/,35X,'<FPDIMENSIONS>',/,18X,
     &        'array sizes for full potential calculations',/,1X,79('*')
     &        ,//,10X,'NFPLIM =',I4,5X,'NSFLIM =',I4,//)
C
      END
C*==fppickrules_rot_klm.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE FPPICKRULES_ROT_KLM(N_LMQ,K_LMQ,N_LMQORG,K_LMQORG,DROT,
     &                               NLMTOP,NLMTOPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   if a rotation DROT moves a site q_org to new position q        *
C   *   then the (l,m) expansion will change from K_LMQORG to K_LMQ    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLMTOP,NLMTOPMAX,N_LMQ,N_LMQORG
      REAL*8 DROT(NLMTOPMAX,NLMTOPMAX)
      INTEGER K_LMQ(NLMTOP),K_LMQORG(NLMTOP)
C
C Local variables
C
      INTEGER LM1,LM2
C
C*** End of declarations rewritten by SPAG
C
      DO LM1 = 1,NLMTOP
         K_LMQ(LM1) = 0
         DO LM2 = 1,NLMTOP
            IF ( ABS(DROT(LM1,LM2)).GT.1D-6 .AND. K_LMQORG(LM2).NE.0 )
     &           K_LMQ(LM1) = 1
         END DO
      END DO
C
      N_LMQ = 0
      DO LM1 = 1,NLMTOP
         IF ( K_LMQ(LM1).NE.0 ) N_LMQ = N_LMQ + 1
      END DO
C
      N_LMQORG = 0
      DO LM1 = 1,NLMTOP
         IF ( K_LMQORG(LM1).NE.0 ) N_LMQORG = N_LMQORG + 1
      END DO
C
      END
C*==fppickrules_hom_equ.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE FPPICKRULES_HOM_EQU(RMAT,INON0,NRMAT,NRMATMAX)
C   ********************************************************************
C   *                                                                  *
C   *   deal with homogenous equations   (R-1) * V = 0                 *
C   *   transform the coefficient matrix (R-1) to (tri)diagonal form   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-6)
C
C Dummy arguments
C
      INTEGER NRMAT,NRMATMAX
      INTEGER INON0(NRMATMAX)
      REAL*8 RMAT(NRMATMAX,NRMATMAX)
C
C Local variables
C
      REAL*8 ARG,RSUM,SCL
      INTEGER I,J,L
      LOGICAL NONZERO
C
C*** End of declarations rewritten by SPAG
C
      NONZERO(ARG) = ABS(ARG).GT.TOL
C
C----------------------------------------- whip out elements in column j
C
      DO J = 1,NRMAT
         IF ( NONZERO(RMAT(J,J)) ) THEN
            DO I = 1,J - 1
               IF ( NONZERO(RMAT(I,J)) ) THEN
                  SCL = RMAT(I,J)/RMAT(J,J)
                  DO L = J,NRMAT
                     RMAT(I,L) = RMAT(I,L) - SCL*RMAT(J,L)
                  END DO
               END IF
            END DO
            DO I = J + 1,NRMAT
               IF ( NONZERO(RMAT(I,J)) ) THEN
                  SCL = RMAT(I,J)/RMAT(J,J)
                  DO L = J,NRMAT
                     RMAT(I,L) = RMAT(I,L) - SCL*RMAT(J,L)
                  END DO
               END IF
            END DO
         ELSE
            RSUM = 0D0
            DO I = 1,NRMAT
               RSUM = RSUM + ABS(RMAT(I,J))
            END DO
C                        IF ( NONZERO(RSUM) ) WRITE (6,*)
C     &                        '################ J RSUM',J,RSUM
         END IF
      END DO
C
C-----------------------------------------------------------------------
C    if a row  i  of the (tri)diagonal matrix derived from (R-1)'
C    has a non-0 diagonal element and 0-s otherwise
C    the corresponding expansion term SF(i), V(i) ... has to be 0
C-----------------------------------------------------------------------
C
      DO I = 1,NRMAT
         INON0(I) = 0
         IF ( NONZERO(RMAT(I,I)) ) THEN
            INON0(I) = 1
            DO J = I + 1,NRMAT
               IF ( NONZERO(RMAT(I,J)) ) THEN
                  INON0(I) = 0
                  EXIT
               END IF
            END DO
         END IF
      END DO
C
      END
