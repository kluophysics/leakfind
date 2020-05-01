C*==chidrive.f    processed by SPAG 6.70Rc at 09:01 on  8 Mar 2017
      SUBROUTINE CHIDRIVE(KBZI,IPRINT,WRMAT,TASK)
C   ********************************************************************
C   *                                                                  *
C   *  driving routine for linear response calculations                *
C   *  according to TASK it calls                                      *
C   *                                                                  *
C   *  - CHI       magnetic susceptibility                             *
C   *  - CHIXAS    field induced XMXD                                  *
C   *  - GILBERT   Gilbert damping parameter                           *
C   *  - SIGMA     residual resistivity                                *
C   *  - OPT       optical conductivity  (TO BE DONE)                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_LINRESP,ONLY:CHIZ,DDTAUTAUT,TKTKTT,NZ12,NZ12MAX,
     &    ERYDA_EQ_ERYDB,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,
     &    JTTX,WTTJ,ITTMAX,JTTMAX,NTKTKLIN,NTKTKMAX,IKM1_CHI_LIN,
     &    IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN,NAB_CHI_QQ,LAMCOLROW,
     &    NCOLROW,NLINCHIMAX,NLIN23_CHI,NLIN41_CHI,
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
     &    not_store_chiz, use_jtauk_square,
     &    prt_sig_realz, ie0_sig, ie1_sig
c end-mod-xjq
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,BRAVAIS,ABAS,BBAS
      USE MOD_CPA,ONLY:ALPHASRO,CPALVL,USENLCPA,ICPAALG,NCPA
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_KSPACE,ONLY:IBZINT,NKPTS0
      USE MOD_SYMMETRY,ONLY:MROTK,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    NO_SYMMETRY_LINRESP,NSYMMAX
      USE MOD_SITES,ONLY:NQ,NQMAX,NOQ,ICPA,IQBOT_CHI,IQTOP_CHI,NQ_CHI,
     &    IQBOT_TB,IQBOT,IQTOP,NQTB
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC,ITBOT,ITTOP,DOBS_LTEX,
     &    DOBS_TEX_GLO
      USE MOD_CALCMODE,ONLY:UPDATE_EFERMI,SIGMA_LAYER,KKRMODE,
     &    ME_CC_BRA_RWF,PUBLIC_VERSION
      USE MOD_ENERGY,ONLY:NEMAX,NETAB
      USE MOD_FILES,ONLY:LDATSET,DATSET,IOTMP,FOUND_SECTION
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,NLQ,NKM,NOBSMAX,NLMAX
      USE MOD_CONSTANTS,ONLY:CI,PI
      USE MOD_MPI,ONLY:MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHIDRIVE')
C
C Dummy arguments
C
      INTEGER IPRINT,KBZI
      CHARACTER*10 TASK
      LOGICAL WRMAT
C
C Local variables
C
      LOGICAL ALLTKTK,FOUND,LCPAMILLS,LIRK,LIRREDKN,LNLCPAAVG,
     &        LNLCPASYMAVG,NEGF_MODE_SA,SAMENLQ,SYMACCEPTEDKN(:,:),
     &        SYMACCEPTEDMOD(NSYMMAX),SYM_IRR(:,:,:),USESTDKMESH
      REAL*8 CALPHA,FALPHA,KD(3),KD0(3),KNNLCPA(:,:),KNRIJ,KTABKN(:,:,:)
     &       ,PCFG(:),RIJ,RQNLCPA(:,:),WKTABKN(:,:),WSUM
      INTEGER CHIPRINT,I,IALPHA,IA_ERR,ICFG,IESORT(:),IK,IKN,IKNM,
     &        IND0QCLU(:),IOCC,IOCCCFGQ(:,:),IPRINTL,IPROCE(:),IQ,IQCLU,
     &        IQCPA,IQ_IQCLU(:),ISYM,JQ,JQCLU,KTKTK,L1,L2,L3,L4,M,
     &        MAXNKNIRMU,NCFG,NDIMCLU,NDIMCLUSQ,NKMQCLU(:),NKN,
     &        NKNIRMU(:),NKPTS00,NKRED,NKTABKN(:),NKTABKNMAX,NLMQCLU(:),
     &        NQNLCPA,NQNLCPATAB(14,3),NSYMACCEPTEDKN(:),NSYMACCEPTEDMOD
      COMPLEX*16 EIKNRIJ(:,:,:,:)
      INTEGER NLCPACONF
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
      real*8 mem_chiz, max_mem_chiz
      parameter (max_mem_chiz=50.0)!GB
c end-mod-xjq
C
C*** End of declarations rewritten by SPAG
C
      DATA USESTDKMESH/.TRUE./
      DATA NQNLCPATAB/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
     &     1,4,2,0,0,0,0,0,0,0,0,0,0,0,8,32,16/
      DATA CHIPRINT/0/
C
C ======================================================================
C
      ALLOCATABLE IESORT,IPROCE
C
C----------------------------------------------------------------- NLCPA
      ALLOCATABLE IQ_IQCLU,IND0QCLU,RQNLCPA,KNNLCPA,NKMQCLU,NLMQCLU
      ALLOCATABLE KTABKN,WKTABKN,NKTABKN,SYMACCEPTEDKN,NKNIRMU
      ALLOCATABLE NSYMACCEPTEDKN,SYM_IRR,IOCCCFGQ,PCFG
      ALLOCATABLE EIKNRIJ
C
C ======================================================================
C              regime of sites treated for response functions
C ======================================================================
C
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
         IQBOT_CHI = 1
         IQTOP_CHI = NQ
      ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
         IQBOT_CHI = IQBOT_TB
         IQTOP_CHI = IQBOT_TB + NQTB - 1
      ELSE
C
         CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
      END IF
C
      NQ_CHI = IQTOP_CHI - IQBOT_CHI + 1
C
C=======================================================================
C
      ALLOCATE (NAB_CHI_QQ(NQMAX,NQMAX))
      ALLOCATE (NCOLROW(NKMMAX,NQMAX))
      ALLOCATE (LAMCOLROW(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: LAMCOLROW')
C
      ALLOCATE (IESORT(NEMAX),IPROCE(NEMAX))
C
      ALLOCATE (DOBS_LTEX(0:3,NOBSMAX,NLMAX,NTMAX,NEMAX))
      ALLOCATE (DOBS_TEX_GLO(0:3,NOBSMAX,NTMAX,NEMAX))
C
C-----------------------------------------------------------------------
C
      IF ( IBZINT.EQ.0 .AND. KBZI.EQ.1 ) CALL STOP_MESSAGE(ROUTINE,
     &     'TAU  via  cluster  calculation BUT KBZI=1')
      IF ( IBZINT.NE.0 .AND. KBZI.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'TAU  via  BZ-integration  BUT KBZI=0')
      IF ( FULLPOT .AND. TASK.NE.'PHONONS   ' .AND. 
     &     TASK.NE.'CHI       ' .AND. TASK(1:5).NE.'SIGMA' .AND. 
     &     TASK(1:6).NE.'MAGNET' .AND. TASK(1:6).NE.'CHIDYN' )
     &     CALL STOP_MESSAGE(ROUTINE,
     &       'the chosen TASK does not YET work for FULL POTENTIAL mode'
     &       )
C
      UPDATE_EFERMI = .FALSE.
C
C ======================================================================
C                          CHI-LANDAU
C                   real space formulation
C ======================================================================
C
      IF ( TASK.EQ.'CHILANDAU ' ) THEN
C
         CALL CHILANDAU(IESORT,IPROCE)
C
         STOP
      END IF
C
C ======================================================================
C                       response functions
C ======================================================================
C
      ALLTKTK = .TRUE.
      ERYDA_EQ_ERYDB = .FALSE.
C
      IF ( TASK(1:5).EQ.'SIGMA' ) THEN
C
         ERYDA_EQ_ERYDB = .FALSE.
         KTKTK = 3
         NZ12 = 2
C
      ELSE IF ( TASK.EQ.'GILBERT   ' ) THEN
C
         ERYDA_EQ_ERYDB = .FALSE.
         KTKTK = 2
         NZ12 = 2
C
      ELSE IF ( TASK.EQ.'PUMP      ' ) THEN
C
         ERYDA_EQ_ERYDB = .FALSE.
         KTKTK = 0
         NZ12 = 1
C
      ELSE IF ( TASK.EQ.'PHONONS   ' ) THEN
C
         ERYDA_EQ_ERYDB = .TRUE.
         KTKTK = 0
         NZ12 = 1
C
      ELSE
C
         ERYDA_EQ_ERYDB = .TRUE.
         KTKTK = 0
         NZ12 = 1
C
      END IF
C
      NZ12MAX = NZ12
C
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('TASK',1)
c modified by XJQ; not_store_chiz, for big system, don't store huge matrix chiz
c                  prt_sig_realz, print sig(E) for kappa program using Kubo-Greenwood
      prt_sig_realz = .false.
      CALL SECTION_FIND_KEYWORD('PRT_SIG_REALZ',prt_sig_realz)
      call section_set_integer('IE0_SIG',ie0_sig,1,0)
      call section_set_integer('IE1_SIG',ie1_sig,netab(1),0)
      not_store_chiz = .false.
      use_jtauk_square = .false.
      CALL SECTION_FIND_KEYWORD('NOT_STORE_CHIZ',not_store_chiz)
      CALL SECTION_FIND_KEYWORD('USE_JTAUK_SQUARE',use_jtauk_square)
      if(not_store_chiz .and. use_jtauk_square)
     &  stop 'not_store_chiz .and. use_jtauk_square'
c end-mod-xjq
C
C---------------------------------- overwrite default for ERYDA_EQ_ERYDB
      CALL SECTION_FIND_KEYWORD('Z1EQZ2',FOUND)
      IF ( FOUND ) ERYDA_EQ_ERYDB = .TRUE.
C
C------------------------------------------- overwrite default for KTKTK
      CALL SECTION_SET_INTEGER('KTKTK',KTKTK,9999,0)
      IF ( TASK(1:3).EQ.'CHI' ) CALL SECTION_SET_INTEGER('CHIPRINT',
     &     CHIPRINT,IPRINT,0)
C
      WRITE (6,99006) TASK,ALLTKTK,ERYDA_EQ_ERYDB,KTKTK
C
C ======================================================================
C                         NLCPA - INITIALISATION
C ======================================================================
      IF ( USENLCPA ) THEN
C
C        --------------------------------------------- set default modes
         IPRINTL = IPRINT
C        --- use irreducible set of  k
         LIRK = .TRUE.
C        --- use irreducible set of  kn
         LIRREDKN = .TRUE.
C        --- use Mills for mixing
         LCPAMILLS = ICPAALG.EQ.1
C        --- average site-diagonal blocks of several matrices
         LNLCPAAVG = .TRUE.
C        --- symetrise the average of site-diagonal blocks
         LNLCPASYMAVG = .TRUE.
C        --- assign default value to alphasro
         ALPHASRO = 0D0
C        --------------------------------- update modes using input file
         CALL INPUT_FIND_SECTION('CPA',0)
C
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('NLCPA',USENLCPA)
            IF ( USENLCPA ) THEN
               CALL SECTION_FIND_KEYWORD('NSYMAVG',LNLCPASYMAVG)
               CALL SECTION_FIND_KEYWORD('NIRRBZ',LIRK)
               CALL SECTION_FIND_KEYWORD('NIRRKN',LIRREDKN)
               CALL SECTION_SET_INTEGER('PRINT',IPRINTL,9999,0)
               CALL SECTION_SET_REAL('alpha',ALPHASRO,9999D0,0)
               LNLCPASYMAVG = .NOT.LNLCPASYMAVG
               LIRREDKN = .NOT.LIRREDKN
               LIRK = .NOT.LIRK
            END IF
         END IF
C        ---------------------------------------------------------------
         IF ( NCPA.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NCPA <> 1')
         IQCPA = 0
         DO IQ = 1,NQ
            IF ( ICPA(IQ).EQ.1 ) IQCPA = IQ
         END DO
         IF ( IQCPA.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'IQCPA not set')
         IF ( NOQ(IQCPA).NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'NOQ(IQCPA) .NE. 2')
         IF ( IPRINTL.GT.4 ) WRITE (6,'("nlcpa:: IQCPA,ICPA",5i5)')
     &                              ICPA,IQCPA
C
         NQNLCPA = NQNLCPATAB(BRAVAIS,CPALVL)
         NKN = NQNLCPA
         NDIMCLU = NKM*NQNLCPA
C
         IF ( NQNLCPA.GT.30 ) THEN
            WRITE (6,99004)
            STOP
         END IF
C
         IF ( NQNLCPA.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'NQNLCPA = 0')
C
         ALLOCATE (KNNLCPA(3,NQNLCPA))
         ALLOCATE (RQNLCPA(3,NQNLCPA),IQ_IQCLU(NQNLCPA),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RQNLCPA')
         ALLOCATE (IND0QCLU(NQNLCPA),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IND0QCLU')
         ALLOCATE (NLMQCLU(NQNLCPA),NKMQCLU(NQNLCPA),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NLMQCLU')
C
C        ---------------------------------- get cluster R and Kn vectors
C
         CALL NLCPASET(BRAVAIS,CPALVL,NQNLCPA,RQNLCPA,KNNLCPA)
C
C        ----------------------------- determine matrix offset variables
C
         DO IQCLU = 1,NQNLCPA
            IQ_IQCLU(IQCLU) = IQCPA
         END DO
C
         SAMENLQ = .TRUE.
C
         DO IQCLU = 1,NQNLCPA
            IQ = IQ_IQCLU(IQCLU)
            IF ( IPRINTL.GT.4 ) WRITE (6,'(5x,"iqclu,iq,nqnlcpa",10i5)')
     &                                 IQCLU,IQ,NQNLCPA
            NLMQCLU(IQCLU) = NLQ(IQ)**2
            NKMQCLU(IQCLU) = NKMQ(IQ)
            IF ( NLMQCLU(IQCLU).NE.NLMQCLU(1) ) SAMENLQ = .FALSE.
C
            IF ( IQCLU.EQ.1 ) THEN
               IND0QCLU(1) = 0
            ELSE
               IND0QCLU(IQCLU) = IND0QCLU(IQCLU-1) + NKMQCLU(IQCLU-1)
            END IF
         END DO
         IF ( .NOT.SAMENLQ )
     &        CALL STOP_MESSAGE(ROUTINE,'SAMENLQ = FALSE')
C
C        ---------------------------------- create coarse grained k-mesh
C
         WRITE (*,*) '############################## NKPTS0 ',NKPTS0,NKN
         NKPTS00 = NKPTS0/NKN
C
         ALLOCATE (NKTABKN(NKN),SYMACCEPTEDKN(NSYMMAX,NKN))
         ALLOCATE (NSYMACCEPTEDKN(NKN))
C
         IF ( LIRK ) THEN
C
            DO IKN = 1,NKN
               SYMACCEPTEDKN(1:NSYMMAX,IKN) = SYMACCEPTED(1:NSYMMAX)
               NSYMACCEPTEDKN(IKN) = NSYMACCEPTED
            END DO
C
C-----------------------------------------------------------------------
C              use largest COMMON symmetry group of ALL tiles
C-----------------------------------------------------------------------
            IF ( NKN.GT.1 ) THEN
C
               CALL NLCPAKMESH(BRAVAIS,CPALVL,BBAS,NKTABKN,NKN,IPRINTL,
     &                         MROTK,LIRK,NSYM,SYMACCEPTEDKN,
     &                         NSYMACCEPTEDKN,KNNLCPA,IOTMP,NKPTS00,
     &                         NKRED)
C
               NSYMACCEPTEDKN(1) = 0
               DO ISYM = 1,NSYM
                  DO IKN = 2,NKN
                     SYMACCEPTEDKN(ISYM,1) = SYMACCEPTEDKN(ISYM,1) .AND. 
     &                  SYMACCEPTEDKN(ISYM,IKN)
                  END DO
                  IF ( SYMACCEPTEDKN(ISYM,1) ) NSYMACCEPTEDKN(1)
     &                 = NSYMACCEPTEDKN(1) + 1
               END DO
C
               DO IKN = 2,NKN
                  SYMACCEPTEDKN(1:NSYMMAX,IKN)
     &               = SYMACCEPTEDKN(1:NSYMMAX,1)
                  NSYMACCEPTEDKN(IKN) = NSYMACCEPTEDKN(1)
               END DO
C
            END IF
C
         ELSE
C
            SYMACCEPTEDKN(1,1:NKN) = .TRUE.
            SYMACCEPTEDKN(2:NSYMMAX,1:NKN) = .FALSE.
            NSYMACCEPTEDKN(1:NKN) = 1
         END IF
C
         SYMACCEPTEDMOD(1:NSYMMAX) = SYMACCEPTEDKN(1:NSYMMAX,1)
         NSYMACCEPTEDMOD = NSYMACCEPTEDKN(1)
C
         WRITE (6,99002) NKN
C
         CLOSE (IOTMP)
C
C        -------------------------------- determine reducible Kn vectors
         CALL REDUCEKN(NKN,KNNLCPA,SYMACCEPTEDMOD,MROTK,NSYM,IPRINTL,
     &                 IOTMP,LIRREDKN,ABAS)
         ALLOCATE (NKNIRMU(NKN))
C        -----------------------------------------  read infos from disk
         REWIND (IOTMP)
         READ (IOTMP) (NKNIRMU(IKN),IKN=1,NKN)
         MAXNKNIRMU = MAXVAL(NKNIRMU)
         ALLOCATE (SYM_IRR(NKN,MAXNKNIRMU,NSYM))
         DO IKN = 1,NKN
            IF ( NKNIRMU(IKN).GT.1 ) READ (IOTMP)
     &           (SYM_IRR(IKN,I,1:NSYM),I=1,NKNIRMU(IKN)-1)
         END DO
         CLOSE (IOTMP)
C
         IF ( LIRREDKN ) WRITE (6,99005) NKN
C        ---------------------------------------------------------------
         IF ( NKN.EQ.1 .AND. USESTDKMESH ) THEN
C
            CALL KMESHS(IOTMP,NSYM,SYMACCEPTEDMOD,MROTK,IPRINT,NKPTS0,
     &                  BBAS,DATSET,LDATSET)
C
         ELSE
C
            CALL NLCPAKMESH(BRAVAIS,CPALVL,BBAS,NKTABKN,NKN,IPRINTL,
     &                      MROTK,LIRK,NSYM,SYMACCEPTEDKN,
     &                      NSYMACCEPTEDKN,KNNLCPA,IOTMP,NKPTS00,NKRED)
C
         END IF
C
         REWIND (IOTMP)
         READ (IOTMP) (NKTABKN(IKN),IKN=1,NKN)
         NKTABKNMAX = MAXVAL(NKTABKN(1:NKN))
         ALLOCATE (KTABKN(3,NKTABKNMAX,NKN))
         ALLOCATE (WKTABKN(NKTABKNMAX,NKN),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WKTABKN')
C
         READ (IOTMP) (((KTABKN(I,IK,IKN),I=1,3),WKTABKN(IK,IKN),IK=1,
     &                NKTABKN(IKN)),IKN=1,NKN)
         CLOSE (IOTMP)
C
         WRITE (6,99003)
         WRITE (6,'(5x,65("-"))')
         DO IKN = 1,NKN
            WRITE (6,'(2x,i5," | ",3f10.6," | ",i5," | ",i5)') IKN,
     &             (KNNLCPA(I,IKN),I=1,3),NKTABKN(IKN),
     &             NSYMACCEPTEDKN(IKN)
         END DO
C
         NCFG = 2**NQNLCPA
C
C------------------------------------------------------ occupation table
C -----  constructing pcfg(x_A,x_B,alpha) harwired for NQNLCPA=2 or 4
C -----  prescription
C -----        site     pcfg
C -----        1  2
C ----- occup  A  A     x_A * x_B + f(alphasro)
C ----- occup  A  B     x_A * x_B + f(alphasro)
C -----       . . .     . . .
C ----- f(alphsro) is a hardwired function for the moment
C
         ALLOCATE (IOCCCFGQ(NCFG,NQNLCPA),PCFG(NCFG))
C
         WRITE (6,99001) ALPHASRO
         WSUM = 0D0
         CALPHA = 3D0*(CONC(1)**2D0)*(CONC(2)**2D0)
C
         DO ICFG = 1,NCFG
            IALPHA = 0
            FALPHA = 1D0/4D0
            PCFG(ICFG) = 1D0
C
            DO IQCLU = 1,NQNLCPA
               IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
               IOCCCFGQ(ICFG,IQCLU) = IOCC
               PCFG(ICFG) = PCFG(ICFG)*CONC(IOCC)
               IF ( IOCC.EQ.1 ) IALPHA = IALPHA + 1
            END DO
C ----- hardwired for 2 site cluster
            IF ( NQNLCPA.EQ.2 ) THEN
               IF ( IALPHA.EQ.1 ) FALPHA = -1D0/4D0
            END IF
C alternative hardwiring for 2 site cluster
            IF ( NQNLCPA.EQ.2 .AND. ALPHASRO.GT.0.01D0 ) THEN
C               DO ICFG = 1,NCFG
C                  PCFG(ICFG) = 0.0d0
C               END DO
               PCFG(1:NCFG) = 0.0D0
C
               PCFG(1) = CONC(1)
               PCFG(4) = CONC(2)
            END IF
            IF ( NQNLCPA.EQ.2 .AND. ALPHASRO.LT.-0.01D0 ) THEN
               IF ( CONC(1).GE.0.5D0 ) THEN
                  PCFG(1) = 2.0D0*CONC(1) - 1.0D0
                  PCFG(2) = 1.0D0 - CONC(1)
                  PCFG(3) = 1.0D0 - CONC(1)
                  PCFG(4) = 0.0D0
               END IF
C
               IF ( CONC(1).LT.0.5D0 ) THEN
                  PCFG(1) = 0.0D0
                  PCFG(2) = CONC(1)
                  PCFG(3) = CONC(1)
                  PCFG(4) = 1.0D0 - 2.0D0*CONC(1)
               END IF
            END IF
C
C
C
C            IF( NQNLCPA.EQ.2.AND.CONC(1).EQ.0.5d0)
            IF ( NQNLCPA.EQ.2 .AND. CONC(1).GE.0.5D0 .AND. CONC(2)
     &           .GE.0.5D0 ) THEN
               PCFG(1) = CONC(1)**2 + ALPHASRO*1D0/4D0
               PCFG(2) = CONC(1)*CONC(2) - ALPHASRO*1D0/4D0
               PCFG(3) = CONC(1)*CONC(2) - ALPHASRO*1D0/4D0
               PCFG(4) = CONC(2)**2 + ALPHASRO*1D0/4D0
               WRITE (*,*) 'PCFG'
               WRITE (*,*) PCFG(1:4)
            END IF
C
C
C ----- hardwired for 4 site cluster (see total E paper Derwyn)
            IF ( NQNLCPA.EQ.4 ) THEN
               IF ( IALPHA.EQ.0 .OR. IALPHA.EQ.4 ) FALPHA = CALPHA
               IF ( IALPHA.EQ.2 ) FALPHA = -CALPHA/3D0
               IF ( IALPHA.EQ.1 .OR. IALPHA.EQ.3 ) FALPHA = 0D0
            END IF
C
C            PCFG(ICFG) = PCFG(ICFG) + FALPHA*ALPHASRO
            WSUM = WSUM + PCFG(ICFG)
         END DO
C
         PCFG(1:NCFG) = PCFG(1:NCFG)/WSUM
C alternative hardwiring for 4 site cluster
         IF ( NQNLCPA.EQ.4 .AND. ALPHASRO.GT.0.01D0 ) THEN
            DO ICFG = 1,NCFG
               PCFG(ICFG) = 0.0D0
            END DO
            PCFG(1) = CONC(1)
            PCFG(16) = CONC(2)
         END IF
         IF ( NQNLCPA.EQ.4 .AND. ALPHASRO.LT.-0.01D0 ) THEN
            WSUM = 0D0
            DO ICFG = 1,NCFG
               IALPHA = 0
               DO IQCLU = 1,NQNLCPA
                  IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
                  IF ( IOCC.EQ.1 ) IALPHA = IALPHA + 1
               END DO
               IF ( CONC(1).GE.0.5D0 .AND. CONC(1).LT.0.75D0 ) THEN
                  IF ( IALPHA.EQ.0 .OR. IALPHA.EQ.4 .OR. IALPHA.EQ.1 )
     &                 PCFG(ICFG) = 0.0D0
                  IF ( IALPHA.EQ.2 ) PCFG(ICFG) = (3.0D0-4.0D0*CONC(1))
     &                 /6.0D0
                  IF ( IALPHA.EQ.3 ) PCFG(ICFG) = CONC(1) - 0.5D0
               END IF
               IF ( CONC(1).GT.0.25D0 .AND. CONC(1).LE.0.5D0 ) THEN
                  IF ( IALPHA.EQ.0 .OR. IALPHA.EQ.4 .OR. IALPHA.EQ.3 )
     &                 PCFG(ICFG) = 0.0D0
                  IF ( IALPHA.EQ.2 ) PCFG(ICFG) = (4.0D0*CONC(1)-1.0D0)
     &                 /6.0D0
                  IF ( IALPHA.EQ.1 ) PCFG(ICFG) = 0.5D0 - CONC(1)
               END IF
               IF ( CONC(1).GE.0.75D0 .AND. CONC(1).LE.1.0D0 ) THEN
                  IF ( IALPHA.EQ.0 .OR. IALPHA.EQ.1 .OR. IALPHA.EQ.2 )
     &                 PCFG(ICFG) = 0.0D0
                  IF ( IALPHA.EQ.3 ) PCFG(ICFG) = CONC(2)
                  IF ( IALPHA.EQ.4 ) PCFG(ICFG) = 4.0D0*CONC(1) - 3.0D0
               END IF
               IF ( CONC(1).LE.0.25D0 .AND. CONC(1).GE.0.0D0 ) THEN
                  IF ( IALPHA.EQ.4 .OR. IALPHA.EQ.3 .OR. IALPHA.EQ.2 )
     &                 PCFG(ICFG) = 0.0D0
                  IF ( IALPHA.EQ.1 ) PCFG(ICFG) = CONC(1)
                  IF ( IALPHA.EQ.0 ) PCFG(ICFG) = 1.0D0 - 4.0D0*CONC(1)
               END IF
               WSUM = WSUM + PCFG(ICFG)
            END DO
            PCFG(1:NCFG) = PCFG(1:NCFG)/WSUM
         END IF
C
         DO ICFG = 1,NCFG
            WRITE (*,*) ICFG,PCFG(ICFG)
         END DO
C
C
C-----------continuous variation of ALPHASRO from -1/3 to 1--------
C
C         IF( NQNLCPA.EQ.4.AND.CONC(1).EQ.0.75D0) THEN
C            write(*,*)'TEST A'
C               PCFG(2) = CONC(2)*(CONC(1))**3D0*(1D0-ALPHASRO)**3D0
C               PCFG(3) = CONC(2)*(CONC(1))**3D0*(1D0-ALPHASRO)**3D0
C               PCFG(4) = CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(5) = CONC(2)*(CONC(1))**3D0*(1D0-ALPHASRO)**3D0
C               PCFG(6) =  CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(7) =  CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(8) = CONC(2)*CONC(1)*(1D0-ALPHASRO)*(1-CONC(1)*(1D0-
C     &                   ALPHASRO))**2D0
C               PCFG(9) = CONC(2)*(CONC(1))**3D0*(1D0-ALPHASRO)**3D0
C               PCFG(10) =  CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(11) =  CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(12) = CONC(2)*CONC(1)*(1D0-ALPHASRO)*(1-CONC(1)*(1D0
C     &                    -ALPHASRO))**2D0
C               PCFG(13) =  CONC(2)*(CONC(1))**2D0*
C     &                   (1D0-ALPHASRO)**2D0*(1-CONC(1)*(1D0-ALPHASRO))
C               PCFG(14) = CONC(2)*CONC(1)*(1D0-ALPHASRO)*(1-CONC(1)*(1D0
C     &                    -ALPHASRO))**2D0
C               PCFG(15) = CONC(2)*CONC(1)*(1D0-ALPHASRO)*(1-CONC(1)*(1D0
C     &                    -ALPHASRO))**2D0
C               PCFG(16) = CONC(2)*(1-CONC(1)*(1D0-ALPHASRO))**3D0
C               PCFG(1) =  1D0-4D0*PCFG(2)-6D0*PCFG(4)-4D0*PCFG(8)-
C     &                    PCFG(16)
C            END IF
C
C------------------------------------------------------------
C
         DO ICFG = 1,NCFG
            IF ( PCFG(ICFG).LT.10D0**(-5D0) ) PCFG(ICFG) = 0D0
         END DO
C
         WSUM = 0
         DO ICFG = 1,NCFG
            WRITE (*,*) ICFG,PCFG(ICFG)
            WSUM = WSUM + PCFG(ICFG)
         END DO
         WRITE (*,*) 'WSUM',WSUM
C
C         IF ( IPRINTL.GT.3 ) THEN
C            WRITE (*,'("   cfg        pgamma")')
C            FMTSTR = '(2i3,f10.6)'
C            IF ( NQNLCPA.EQ.4 ) WRITE (FMTSTR(2:2),'(a1)') '4'
C            DO ICFG = 1,NCFG
C               WRITE (*,*) (IOCCCFGQ(ICFG,IQCLU),IQCLU=1,NQNLCPA),
C     &              PCFG(ICFG)
C            END DO
C         END IF
C
C
C        ------------------ set up the phase factors for coarse graining
         ALLOCATE (EIKNRIJ(NQNLCPA,MAXNKNIRMU,NQNLCPA,NQNLCPA),
     &             STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: EIKNRIJ')
C
         DO IKN = 1,NKN
C           --- take care of reducible kvecs obtained by rot of irr kvec
            DO IKNM = 1,NKNIRMU(IKN)
C
               KD = KNNLCPA(1:3,IKN)
C
               IF ( IKNM.GT.1 ) THEN
                  DO ISYM = 1,NSYM
                     IF ( SYM_IRR(IKN,IKNM-1,ISYM) ) THEN
C
                        KD0(1:3) = KD(1:3)
                        CALL DGEMV('T',3,3,1D0,MROTK(1,1,ISYM),3,KD0,1,
     &                             0D0,KD,1)
C
                        GOTO 10
                     END IF
                  END DO
                  STOP
               END IF
C
 10            CONTINUE
               DO JQCLU = 1,NQNLCPA
                  DO IQCLU = 1,NQNLCPA
                     KNRIJ = 0D0
                     DO I = 1,3
                        RIJ = RQNLCPA(I,IQCLU) - RQNLCPA(I,JQCLU)
                        KNRIJ = KNRIJ + KD(I)*RIJ
                     END DO
                     EIKNRIJ(IKN,IKNM,IQCLU,JQCLU) = EXP(CI*2*PI*KNRIJ)
                  END DO
               END DO
            END DO
         END DO
C
      ELSE
C------------------------------------------- NO NLCPA: use full symmetry
C
         SYMACCEPTEDMOD(1:NSYMMAX) = SYMACCEPTED(1:NSYMMAX)
         NSYMACCEPTEDMOD = NSYMACCEPTED
C
      END IF
C
C ======================================================================
C                       response functions
C ======================================================================
C
      IF ( TASK(1:3).EQ.'CHI' ) THEN
C---------------------------- restrict TAU(k)*TAU(k) via selection rules
C
         CALL CHITKTKTAB(NAB_CHI_QQ,NCOLROW,LAMCOLROW,KTKTK,
     &                   ERYDA_EQ_ERYDB,NTKTKLIN,NLIN41_CHI,NLIN23_CHI)
C
         REWIND (IOTMP)
         READ (IOTMP) NLINCHIMAX,NTKTKMAX
C
         ALLOCATE (IKM1_CHI_LIN(NLINCHIMAX),IKM2_CHI_LIN(NLINCHIMAX))
         ALLOCATE (IKM3_CHI_LIN(NLINCHIMAX),IKM4_CHI_LIN(NLINCHIMAX),
     &             STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ALLOC: IKM3_CHI_LIN')
C
         READ (IOTMP) (IKM1_CHI_LIN(I),IKM2_CHI_LIN(I),IKM3_CHI_LIN(I),
     &                IKM4_CHI_LIN(I),I=1,NLINCHIMAX)
         CLOSE (IOTMP)
C
      ELSE
C------------------------------- consider ALL TAU(k)*TAU(k) combinations
C
         NTKTKMAX = NKM**4
C
      END IF
C
C-----------------------------------------------------------------------
C                     2D layered systems via Landauer
C-----------------------------------------------------------------------
C
      IF ( SIGMA_LAYER .OR. SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         WRITE (6,*) '<CHIDRIVE>: IT-/IQBOT -- IT-/IQTOP = 1 -- NT/NQ'
C
         ITBOT = 1
         ITTOP = NT
         IQBOT = 1
         IQTOP = NQ
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         JTTMAX = 1
C
         ALLOCATE (WTTJ(JTTMAX))
         ALLOCATE (JTT1(JTTMAX),JTT2(JTTMAX),JTTX(JTTMAX))
C
         ITTMAX = NTKTKMAX*NQ*NQ
         ALLOCATE (NTTJ(ITTMAX))
         ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX))
         ALLOCATE (ITTD(ITTMAX),ITTQ1(ITTMAX),ITTQ2(ITTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ITTD')
C
         I = 0
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO JQ = IQBOT_CHI,IQTOP_CHI
               DO L1 = 1,NKM
                  DO L4 = 1,NKM
                     DO L2 = 1,NKM
                        DO L3 = 1,NKM
                           I = I + 1
                           ITTA(I) = L1
                           ITTB(I) = L2
                           ITTC(I) = L3
                           ITTD(I) = L4
                           ITTQ1(I) = IQ
                           ITTQ2(I) = JQ
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE IF ( .NOT.NO_SYMMETRY_LINRESP ) THEN
C
C-----------------------------------------------------------------------
C
         CALL CHITKTKSYM(SYMACCEPTEDMOD,NSYMACCEPTEDMOD)
C
      ELSE
C
C-----------------------------------------------------------------------
C
         WRITE (6,99007)
C
      END IF
C-----------------------------------------------------------------------
C
      M = NKMMAX
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C NOTE: + Save RAM space by allocating CHIZ by master process only
C         (but only in case of dealing with one energy point)
C
c modified by XJQ; not_store_chiz or use_jtauk_square, for big system, don't store huge matrix chiz
      mem_chiz = real(M*M*NQ_CHI)**2.0*NZ12MAX*16.0/1024.0/1024.0/1024.0
      if(mpi_id==0) then
        write(*,*)
        write(*,*) 'M=',M,',NQ=',NQ_CHI,',NZ12=',NZ12MAX
        write(*,*) 'complex chiz( M*M*NQ_CHI, M*M*NQ_CHI, NZ12MAX )'
        write(*,*) 'min mem for sigkloop is 2 * mem of chiz: ',
     &             2.0*mem_chiz,' GB'
        write(*,*)
      endif
      if( mem_chiz > max_mem_chiz .and. .not.not_store_chiz .and.
     &    .not.use_jtauk_square ) 
     &  stop 'too large memory of chiz, set NOT_STORE_CHIZ or 
     &        USE_JTAUK_SQUARE'
      if( not_store_chiz .or. use_jtauk_square ) then ! do nothing
        if(mpi_id==0) write(*,*) 'use sig1_surf_notstorechiz'
        if(mpi_id==0) write(*,*)
c      IF ( ((TASK.EQ.'SIGMA' .AND. (NETAB(1).EQ.1)) .OR. 
      ELSEIF ( ((TASK.EQ.'SIGMA' .AND. (NETAB(1).EQ.1)) .OR. 
c end-mod-xjq
     &     TASK.EQ.'PHONONS') .AND. MPI_ID.NE.0 .AND. KKRMODE(1:2)
     &     .NE.'TB' ) THEN
C
         ALLOCATE (CHIZ(1,1,1),STAT=IA_ERR)
C
      ELSE IF ( TASK.NE.'NEGF' ) THEN
C
         ALLOCATE (CHIZ(M*M*NQ_CHI,M*M*NQ_CHI,NZ12MAX),STAT=IA_ERR)
C
      END IF
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CHIZ')
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( TASK(1:3).EQ.'CHI' ) THEN
         ALLOCATE (DDTAUTAUT(NTKTKMAX,NTMAX))
         ALLOCATE (TKTKTT(NTKTKMAX,NTMAX,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TKTKTT')
      END IF
C
C ======================================================================
C                   mode for matrix element calculation
C
C   apply Complex Conjugation to radial wave function of <BRA| state
C   for matrix elements specific to spectroscopy
C ======================================================================
C
      IF ( TASK.EQ.'CHIXAS    ' ) THEN
C
         ME_CC_BRA_RWF = .TRUE.
C
      ELSE
C
         ME_CC_BRA_RWF = .FALSE.
C
      END IF
C
C ======================================================================
C                           INITALISATION - END
C ======================================================================
C
C
C ======================================================================
C                         perform requested TASK
C ======================================================================
C
C ======================================================================
C                        MAGNETIC SUSCEPTIBILITY
C ======================================================================
C
      IF ( TASK.EQ.'CHI       ' ) THEN
C
         IF ( FULLPOT ) THEN
C
C---------------------------------------- select calculation mode STATIC
            CALL SECTION_FIND_KEYWORD('STAT',FOUND)
            FOUND = .TRUE.
C
            IF ( FOUND ) THEN
C
               IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &           'PUBLIC VERSION:  TASK = CHI AND FULLPOT not supported'
     &           )
C
               CALL LINRESP_SUSC_STAT
C
            END IF
C
         ELSE
C
            CALL CHISUSC(IPROCE,CHIPRINT)
C
         END IF
C
C ======================================================================
C                        MAGNETISATION
C ======================================================================
      ELSE IF ( TASK.EQ.'MAGNET    ' ) THEN
C
         CALL LINRESP_MAGNET
C
C ======================================================================
      ELSE IF ( TASK.EQ.'CHIXAS    ' ) THEN
C
         CALL CHIXAS
C
C ======================================================================
C                        TRANSPORT CALCULATIONS
C ======================================================================
C
      ELSE IF ( TASK.EQ.'SIGMA     ' ) THEN
C
C-----------------------------------------------------------------------
C                     2D layered systems via Landauer
C-----------------------------------------------------------------------
         IF ( SIGMA_LAYER ) THEN
C
            IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &           'PUBLIC VERSION:  SIGMA_LAYER not supported')
C
            CALL SIGLAYER
C
C-----------------------------------------------------------------------
C                          Kubo formalism
C-----------------------------------------------------------------------
C
         ELSE IF ( USENLCPA ) THEN
C
            CALL SIGNLC(IPRINT,WRMAT,NDIMCLUSQ,IND0QCLU,RQNLCPA,KNNLCPA,
     &                  NSYMACCEPTEDKN,KTABKN,WKTABKN,NKTABKN,
     &                  SYMACCEPTEDKN,NKNIRMU,SYM_IRR,PCFG,MAXNKNIRMU,
     &                  EIKNRIJ,NCFG,NDIMCLU,NKN,NKTABKNMAX,NQNLCPA,
     &                  LCPAMILLS,LNLCPASYMAVG,LNLCPAAVG,IPRINTL,IQCPA)
C
C
         ELSE
C
            CALL SIG
C
         END IF
C
C ======================================================================
C                  Non-Equilibrium Green's Function Mode
C ======================================================================
C
      ELSE IF ( TASK.EQ.'NEGF' ) THEN
C
         IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &        'PUBLIC VERSION:  TASK = NEGF not supported')
C
         IF ( SYSTEM_DIMENSION(1:2).NE.'2D' )
     &         CALL STOP_MESSAGE(ROUTINE,'SYSTEM_DIMENSION <> 2D ')
C
         NEGF_MODE_SA = .FALSE.
         CALL SECTION_FIND_KEYWORD('MODE_SA',NEGF_MODE_SA)
C
         IF ( NEGF_MODE_SA ) THEN
            CALL NEGFDRIVE_SA
         ELSE
            CALL NEGFDRIVE
         END IF
C
C ======================================================================
C                        GILBERT DAMPING
C ======================================================================
C
      ELSE IF ( TASK.EQ.'GILBERT   ' ) THEN
C
         IF ( USENLCPA ) THEN
C
            CALL GILNLC(IPRINT,WRMAT,NDIMCLUSQ,IND0QCLU,RQNLCPA,KNNLCPA,
     &                  NSYMACCEPTEDKN,KTABKN,WKTABKN,NKTABKN,
     &                  SYMACCEPTEDKN,NKNIRMU,SYM_IRR,PCFG,MAXNKNIRMU,
     &                  EIKNRIJ,NCFG,NDIMCLU,NKN,NKTABKNMAX,NQNLCPA,
     &                  LCPAMILLS,LNLCPASYMAVG,LNLCPAAVG,IPRINTL,IQCPA)
C
         ELSE
C
            CALL GIL
C
         END IF
C
C ======================================================================
C                             PUMP
C ======================================================================
C
      ELSE IF ( TASK.EQ.'PUMP      ' ) THEN
C
         IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &        'PUBLIC VERSION:  TASK = PUMP not supported')
C
         CALL LINRESP_PUMP
C
C ======================================================================
C                             PHONONS
C ======================================================================
C
      ELSE IF ( TASK.EQ.'PHONONS   ' ) THEN
C
         IF ( PUBLIC_VERSION ) CALL STOP_REGULAR(ROUTINE,
     &        'PUBLIC VERSION:  TASK = PHONONS not supported')
C
C---------------------- select mode SINGLE particle energy approximation
C
         CALL SECTION_FIND_KEYWORD('SINGLE',FOUND)
C
         IF ( FOUND ) CALL LINRESP_PHON_SINGLE
C
C-------------------------------------------------- full SCF calculation
C
         CALL LINRESP_PHONONS
C
         STOP
C
C ======================================================================
      ELSE
C
         WRITE (6,*) '########################## TASK = ',TASK
         WRITE (6,*) '##  allowed: TASK = CHI, CHIXAS, CHILANDAU, SIGMA'
         WRITE (6,*) '##  allowed: TASK = GILBERT'
         CALL STOP_MESSAGE(ROUTINE,'TASK NOT allowed')
      END IF
C
      DEALLOCATE (IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN)
      DEALLOCATE (ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ)
      DEALLOCATE (JTT1,JTT2,JTTX,WTTJ,CHIZ,NAB_CHI_QQ)
      DEALLOCATE (LAMCOLROW,NCOLROW,IESORT,IPROCE)
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
99001 FORMAT (/,5X,'SRO - parameter ALPHA           ',F6.3)
99002 FORMAT (//,1X,79('*'),/,25X,'NLCPA - coarse graining of BZ',/,1X,
     &        79('*'),//,5X,'number of tiles ',i5,/)
99003 FORMAT (/,'     ikn       kx        ky        kz',
     &        '      # of ks   nsym')
99004 FORMAT ('nlcpa:: too many atoms in cluster, function nlcpaconf ',
     &        'will fail on 32 bit machine')
99005 FORMAT (/,5x,'Number of irreducible tiles  :',i5,/)
99006 FORMAT (/,10X,'Linear response calculation for  TASK = ',A,//,10X,
     &        'ALLTKTK   = ',L3,/,10X,'(E_a=E_b) = ',L3,/,10X,
     &        'KTKTK     = ',I3,/)
99007 FORMAT (/,10X,'no symmtry used for linear response calculation',/)
      END
