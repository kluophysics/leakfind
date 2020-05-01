C*==nlcpa.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPA(NLCPAOUT,NLCPAMAD,ICPAFLAG,CPAERR,IPRINT,WRMAT,
     &                 WRTAUMQ,IFILTAU,TAUQ,ITCPA,ICPACONV,MSSQ,MSST,
     &                 TAUT,IECURR,ERYD,P,MEZZ,MEZJ,NKPTS0)
C   ********************************************************************
C   *                                                                  *
C   *  DK, HE, January 2006                                            *
C   *                                                                  *
C   *  Non-local CPA                                                   *
C   *  -------------                                                   *
C   *  following NLCPA algorithm of D.A. Rowlands et al.               *
C   *  Phys. Rev. B 67, 115109 (2003)                                  *
C   *                                                                  *
C   *  with modifications                                              *
C   *    + formulation which exludes the calculation and use of        *
C   *      real space structure constants                              *
C   *                                                                  *
C   *  NLCPAOUT determines additional output created     PROGRAM       *
C   *    0      none                                                   *
C   *    1      correction to the Madelung potential       SCF         *
C   *    2      expectation values for each config.        GEN         *
C   *           written to file if IECURR = NETAB                      *
C   *    3      DOS for each configuration                 GEN         *
C   *           written to file if IECURR = NETAB                      *
C   *                                                                  *
C   *  ICPAALG = 1   MILLS algorithm                                   *
C   *            4   simple mixing                                     *
C   *                                                                  *
C   ********************************************************************
C
C    ------------------------------------------------------------ to do
C    + put in alphasro appropriately
C    + allow for hexagonal lattice
C    + check again convention for rotation
C    + adjust write statements/printlevel
C    + recover cpa for nlcpalvl = 1
C    ------------------------------------------------------------ to do
C
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      USE MOD_ENERGY,ONLY:EFERMI,WETAB,NETAB
      USE MOD_SYMMETRY,ONLY:DROT,MROTK,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,IQORGQP,NSYMMAX
      USE MOD_LATTICE,ONLY:BRAVAIS,ABAS,BBAS
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKMQ,NLQ,NKM
      USE MOD_TYPES,ONLY:NTMAX,CONC
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET,
     &    FOUND_SECTION
      USE MOD_CPA,ONLY:CPALVL,ALPHASRO,CPAMIX,CPATOL,ICPAALG,ITCPAMAX,
     &    NCPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NLCPA')
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C Dummy arguments
C
      REAL*8 CPAERR,NLCPAMAD
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IECURR,IFILTAU,IPRINT,ITCPA,NKPTS0,
     &        NLCPAOUT
      LOGICAL WRMAT,WRTAUMQ
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 CALPHA,FALPHA,KNNLCPA(:,:),KTABKN(:,:,:),PCFG(:),
     &       RQNLCPA(:,:),WGIT(:),WKTABKN(:,:),WSUM
      COMPLEX*16 CMATTRC
      COMPLEX*16 CSUM,EVD(:,:,:),EVDZJ(:,:),EVI(:,:,:),MMAT(:,:),
     &           MUEHAT(:,:),OMEGAHAT(:,:),TAUHAT(:,:),TAUII(:,:),
     &           TAUIMP(:,:),W1(:,:),W1HAT(:,:)
      INTEGER I,I0,I0_QU(:,:),I1_Q(:),I1_QU(:,:),I2_Q(:),I2_QU(:,:),
     &        IALPHA,IA_ERR,ICALL,ICFG,IK,IKN,IME,IND0QCLU(:),IOCC,
     &        IOCCCFGQ(:,:),IPRINTL,IQ,IQCLU,IQCPA,IQ_IQCLU(:),IT,IU,
     &        IWRI,M,MAXNKNIRMU,N,NAQ,NAQU,NA_Q(:),NCFG,NDIMCLU,
     &        NDIMCLUSQ,NKMQCLU(:),NKN,NKNIRMU(:),NKPTS00,NKRED,
     &        NKTABKN(:),NKTABKNMAX,NLMQCLU(:),NME,NQNLCPA,
     &        NQNLCPATAB(14,3),NSYMACCEPTEDKN(:),NU
      LOGICAL LCPAMILLS,LIRK,LIRREDKN,LNLCPAAVG,LNLCPASYMAVG,SAMENLQ,
     &        SYMACCEPTEDKN(:,:),SYM_IRR(:,:,:),USENLCPA,USESTDKMESH
      INTEGER NLCPACONF
      SAVE EVI,IND0QCLU,IOCCCFGQ,IPRINTL,IQCPA,IQ_IQCLU,KNNLCPA,KTABKN,
     &     LCPAMILLS,LIRK,LIRREDKN,LNLCPAAVG,LNLCPASYMAVG,MAXNKNIRMU,
     &     NCFG,NDIMCLU,NKMQCLU,NKN,NKNIRMU,NKTABKN,NKTABKNMAX,NLMQCLU,
     &     NME,NQNLCPA,NSYMACCEPTEDKN,PCFG,RQNLCPA,SYMACCEPTEDKN,
     &     SYM_IRR,WKTABKN
C
C*** End of declarations rewritten by SPAG
C
      DATA USESTDKMESH/.TRUE./
      DATA ICALL/0/
      DATA NQNLCPATAB/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
     &     1,4,2,0,0,0,0,0,0,0,0,0,0,0,8,32,16/
C
      ALLOCATABLE IQ_IQCLU,IND0QCLU,RQNLCPA,KNNLCPA,NKMQCLU,NLMQCLU
      ALLOCATABLE TAUHAT,OMEGAHAT,MUEHAT,NSYMACCEPTEDKN
      ALLOCATABLE KTABKN,WKTABKN,NKTABKN,SYMACCEPTEDKN,NKNIRMU
      ALLOCATABLE SYM_IRR,IOCCCFGQ,PCFG
      ALLOCATABLE MMAT,TAUIMP,W1HAT,WGIT,EVI,EVD,EVDZJ,TAUII,W1
      ALLOCATABLE NA_Q,I0_QU,I1_QU,I2_QU,I1_Q,I2_Q
C
C=======================================================================
      NU = NQNLCPA
C
      ALLOCATE (NA_Q(NQ),I0_QU(NQ,NU),I1_QU(NQ,NU),I2_QU(NQ,NU))
      ALLOCATE (I1_Q(NQ),I2_Q(NQ))
      NAQ = 0
      DO IQ = 1,NQ
         NA_Q(IQ) = NKMQ(IQ)
         NAQ = NAQ + NA_Q(IQ)
         IF ( IQ.EQ.1 ) THEN
            I1_Q(IQ) = 0
         ELSE
            I1_Q(IQ) = I2_Q(IQ-1) + 1
         END IF
         I2_Q(IQ) = I1_Q(IQ) + NA_Q(IQ)
      END DO
      NAQU = NAQ*NU
C
      DO IU = 1,NU
         DO IQ = 1,NQ
            IF ( IQ.EQ.1 ) THEN
               I0_QU(IQ,IU) = (IU-1)*NAQ
            ELSE
               I0_QU(IQ,IU) = I0_QU(IQ-1,IU) + NA_Q(IQ-1)
            END IF
            I1_QU(IQ,IU) = I0_QU(IQ,IU) + 1
            I2_QU(IQ,IU) = I0_QU(IQ,IU) + NA_Q(IQ)
         END DO
      END DO
C
C   ********************************************************************
C                   INITIALISATION - START
C   ********************************************************************
C
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
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
C
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
         IF ( NCPA.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NCPA.NE.1')
         DO IQ = 1,NQ
            IF ( ICPA(IQ).EQ.1 ) IQCPA = IQ
         END DO
         IF ( NOQ(IQCPA).NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'NOQ(IQCPA).NE.2')
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
         IF ( NQNLCPA.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'NQNLCPA.EQ.0')
C
         ALLOCATE (IND0QCLU(NQNLCPA),IQ_IQCLU(NQNLCPA))
         ALLOCATE (NLMQCLU(NQNLCPA),NKMQCLU(NQNLCPA))
         ALLOCATE (KNNLCPA(3,NQNLCPA),RQNLCPA(3,NQNLCPA),STAT=IA_ERR)
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
         IF ( .NOT.SAMENLQ ) CALL STOP_MESSAGE(ROUTINE,'.NOT.SAMENLQ')
C
C        ---------------------------------- create coarse grained k-mesh
C
         ALLOCATE (SYMACCEPTEDKN(NSYMMAX,NKN))
         ALLOCATE (NSYMACCEPTEDKN(NKN))
C
         IF ( LIRK ) THEN
            DO IKN = 1,NKN
               SYMACCEPTEDKN(1:NSYMMAX,IKN) = SYMACCEPTED(1:NSYMMAX)
               NSYMACCEPTEDKN(IKN) = NSYMACCEPTED
            END DO
         ELSE
            SYMACCEPTEDKN(1,1:NKN) = .TRUE.
            SYMACCEPTEDKN(2:NSYMMAX,1:NKN) = .FALSE.
            NSYMACCEPTEDKN(1:NKN) = 1
         END IF
C
         NKPTS00 = NKPTS0/NKN
C
         ALLOCATE (NKTABKN(NKN),STAT=IA_ERR)
C
         WRITE (6,99002) NKN
C
C        -------------------------------- determine reducible Kn vectors
         CALL REDUCEKN(NKN,KNNLCPA,SYMACCEPTED,MROTK,NSYM,IPRINTL,IOTMP,
     &                 LIRREDKN,ABAS)
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
            CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &                  DATSET,LDATSET)
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
         ALLOCATE (WKTABKN(NKTABKNMAX,NKN))
         ALLOCATE (KTABKN(3,NKTABKNMAX,NKN),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NKTABKN')
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
         WRITE (6,'(72("="))')
C
C---------------------------------------------------- expectation values
C
         IF ( NLCPAOUT.LE.1 ) THEN
            NME = 1
         ELSE IF ( NLCPAOUT.EQ.2 ) THEN
            NME = 2
         ELSE IF ( NLCPAOUT.EQ.3 ) THEN
            NME = 4
         END IF
         IF ( NME.GT.NMEMAX ) CALL STOP_MESSAGE(ROUTINE,'NME.GT.NMEMAX')
C
         NCFG = 2**NQNLCPA
C
         ALLOCATE (EVI(NCFG,NQNLCPA,NME),STAT=IA_ERR)
         CALL CINIT(NCFG*NQNLCPA*NME,EVI)
C
         NLCPAMAD = 0D0
C
C------------------------------------------------------ occupation table
C -----  constructing pcfg(x_A,x_B,alpha) hardwired for NQNLCPA=2 or 4
C -----  prescription
C -----        site     pcfg
C -----        1  2
C ----- occup  A  A     x_A * x_B + f(alphasro)
C ----- occup  A  B     x_A * x_B + f(alphasro)
C -----       . . .     . . .
C ----- f(alphsro) is a hardwired function for the moment
C
         ALLOCATE (IOCCCFGQ(NCFG,NQNLCPA),PCFG(NCFG),STAT=IA_ERR)
C
         WRITE (6,99001) ALPHASRO
         WSUM = 0D0
         CALPHA = 3D0*(CONC(1)**2D0)*(CONC(2)**2D0)
         DO ICFG = 1,NCFG
            IALPHA = 0
            FALPHA = 1D0/4D0
            PCFG(ICFG) = 1D0
            DO IQCLU = 1,NQNLCPA
               IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
               IOCCCFGQ(ICFG,IQCLU) = IOCC
               PCFG(ICFG) = PCFG(ICFG)*CONC(IOCC)
               IF ( IOCC.EQ.1 ) IALPHA = IALPHA + 1
            END DO
            IF ( NQNLCPA.EQ.1 ) FALPHA = 0D0
C ----- hardwired for 2 site cluster
            IF ( NQNLCPA.EQ.2 ) THEN
               IF ( IALPHA.EQ.1 ) FALPHA = -1D0/4D0
            END IF
C ----- hardwired for 4 site cluster (see total E paper Derwyn)
            IF ( NQNLCPA.EQ.4 ) THEN
               IF ( IALPHA.EQ.0 .OR. IALPHA.EQ.4 ) FALPHA = CALPHA
               IF ( IALPHA.EQ.2 ) FALPHA = -CALPHA/3D0
               IF ( IALPHA.EQ.1 .OR. IALPHA.EQ.3 ) FALPHA = 0D0
            END IF
            PCFG(ICFG) = PCFG(ICFG) + FALPHA*ALPHASRO
            WSUM = WSUM + PCFG(ICFG)
         END DO
C
         PCFG(1:NCFG) = PCFG(1:NCFG)/WSUM
C
         IF ( IPRINTL.GT.3 ) THEN
            WRITE (6,'("   cfg        pgamma")')
Cc            FMTSTR = '(2i3,f10.6)'
Cc            IF ( NQNLCPA.EQ.4 ) WRITE (FMTSTR(2:2),'(a1)') '4'
            DO ICFG = 1,NCFG
Cc               WRITE (6,FMTSTR(1:10)) (IOCCCFGQ(ICFG,IQCLU),IQCLU=1,
               WRITE (6,*) (IOCCCFGQ(ICFG,IQCLU),IQCLU=1,NQNLCPA),
     &                     PCFG(ICFG)
            END DO
         END IF
C
      END IF
C **********************************************************************
C                   INITALISATION - END
C **********************************************************************
C
C      WRITE (6,*) 'iprintl',IPRINTL,IPRINT
      IF ( IPRINTL.GT.3 ) CALL CMATSTRUCT('TAUQ Start',TAUQ(1,1,1),NKM,
     &     NKMMAX,IREL,IREL,0,1D-8,6)
C
      ALLOCATE (OMEGAHAT(NDIMCLU,NDIMCLU))
      ALLOCATE (TAUHAT(NDIMCLU,NDIMCLU))
      ALLOCATE (MUEHAT(NDIMCLU,NDIMCLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MUEHAT')
C
C     ------------------------------------------- do the NLCPA iteration
C
      CALL NLCPAITER(CPAERR,P,IPRINTL,NKM,NQ,NKMQ,TAUQ,DROT,IQORGQP,
     &               SYMUNITARY,MSSQ,MSST,ITOQ,IREL,NSYM,NTMAX,NQMAX,
     &               NKMMAX,NKTABKNMAX,MROTK,MUEHAT,OMEGAHAT,TAUHAT,
     &               NDIMCLU,IND0QCLU,NQNLCPA,NKTABKN,KTABKN,WKTABKN,
     &               NKN,NSYMACCEPTEDKN,SYMACCEPTEDKN,KNNLCPA,RQNLCPA,
     &               IQCPA,CPATOL,CPAMIX,NKNIRMU,SYM_IRR,MAXNKNIRMU,
     &               LCPAMILLS,LNLCPASYMAVG,LNLCPAAVG,ITCPA,ITCPAMAX,
     &               ICPAFLAG,ICPACONV,PCFG,NCFG,NAQ,NA_Q,NU,NAQU,I0_QU,
     &               I1_QU,I2_QU,I1_Q,I2_Q)
C
C   ********************************************************************
C       determine TAUT and calculate  NME  configuration specific
C           expectation values as controlled via  NLCPAOUT
C   ********************************************************************
C
      ALLOCATE (EVD(NCFG,NQNLCPA,NME),EVDZJ(2,NME),WGIT(NTMAX))
      ALLOCATE (W1(NKMMAX,NKMMAX),TAUII(NKMMAX,NKMMAX))
      M = NDIMCLU
      ALLOCATE (W1HAT(M,M),MMAT(M,M),TAUIMP(M,M),STAT=IA_ERR)
C
      CALL RINIT(NTMAX,WGIT)
      NDIMCLUSQ = NDIMCLU*NDIMCLU
      CALL CINIT(NDIMCLUSQ,MMAT)
C
      M = NDIMCLU
      N = NKMQ(IQCPA)
C
C------------------------------------------------ irregular contribution
      DO IME = 1,NME
         DO IOCC = 1,2
            IT = ITOQ(IOCC,IQCPA)
            EVDZJ(IOCC,IME) = CPRE*CMATTRC(N,NKMMAX,MEZJ(1,1,IT,IME))
         END DO
      END DO
C
C     ----------------------------  zero out TAUT blocks for nlcpa IT''s
C
      DO IOCC = 1,2
         IT = ITOQ(IOCC,IQCPA)
         TAUT(1:N,1:N,IT) = C0
      END DO
C
C=================================================================== CGF
      DO ICFG = 1,NCFG
C        --------------- fill up cavity with configuration of scatterers
C        ----------------------- (determine MMAT for configuration ICFG)
         DO IQCLU = 1,NQNLCPA
            IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
            IT = ITOQ(IOCC,IQCPA)
            I0 = IND0QCLU(IQCLU)
            MMAT(I0+1:I0+N,I0+1:I0+N) = MSST(1:N,1:N,IT)
            IF ( IPRINTL.GT.5 ) WRITE (6,'(5x,"icfg,it,iqclu",9i5)')
     &                                 ICFG,IT,IQCLU
         END DO
C
         W1HAT(1:M,1:M) = MMAT(1:M,1:M) - OMEGAHAT(1:M,1:M)
C
         CALL CMATINV(M,M,W1HAT,TAUIMP)
C
C        -------- get contribution to TAUT(IT) from cluster site IQCLU=1
C
         IOCC = NLCPACONF(ICFG,NQNLCPA,1,1,2)
         IT = ITOQ(IOCC,IQCPA)
         TAUT(1:N,1:N,IT) = TAUT(1:N,1:N,IT) + PCFG(ICFG)
     &                      *TAUIMP(1:N,1:N)
         WGIT(IT) = WGIT(IT) + PCFG(ICFG)
C
C-----------------------------------------------------------------------
         DO IQCLU = 1,NQNLCPA
            IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
            IT = ITOQ(IOCC,IQCPA)
            I0 = IND0QCLU(IQCLU)
            TAUII(1:N,1:N) = TAUIMP(I0+1:I0+N,I0+1:I0+N)
C
            DO IME = 1,NME
C
               CALL CMATMUL(N,NKMMAX,MEZZ(1,1,IT,IME),TAUII,W1)
C
               CSUM = CPRE*CMATTRC(N,NKMMAX,W1)
C
               EVD(ICFG,IQCLU,IME) = CSUM - EVDZJ(IOCC,IME)
               EVI(ICFG,IQCLU,IME) = EVI(ICFG,IQCLU,IME)
     &                               + WETAB(IECURR,1)
     &                               *EVD(ICFG,IQCLU,IME)
C
            END DO
C
         END DO
C-----------------------------------------------------------------------
C
      END DO
C=================================================================== CGF
C
C     --------------------------------------------------- normalise TAUT
      DO IOCC = 1,2
         IT = ITOQ(IOCC,IQCPA)
         TAUT(1:N,1:N,IT) = TAUT(1:N,1:N,IT)/WGIT(IT)
         IF ( IPRINTL.GT.3 ) THEN
            WRITE (6,'(5x,"nlcpa:: taut for it",i5)') IT
            CALL CMATSTRUCT('TAUT      (KAP,MUE)',TAUT(1,1,IT),N,NKMMAX,
     &                      IREL,IREL,0,1D-8,6)
         END IF
      END DO
C
C=T=T=T=T=T=T=T=T=T=T=T=T=T=T=T=T=T=T=T= END  conf average of TAUT: Loop
C
      IF ( NLCPAOUT.GE.2 ) CALL NLCPAWR(NLCPAOUT,IOTMP,DATSET,LDATSET,
     &                                  SYSTEM,LSYSTEM,IQCPA,ITOQ,CONC,
     &                                  EFERMI,IECURR,ERYD,NETAB,EVI,
     &                                  EVD,PCFG,IOCCCFGQ,NQNLCPA,NCFG,
     &                                  NME,NTMAX,NQMAX)
C
      IF ( WRTAUMQ ) CALL NLCTAUIO(IFILTAU,WRTAUMQ,IPRINT,IECURR,ERYD,
     &                             TAUQ,MSSQ,OMEGAHAT,MUEHAT,NQ,NKMQ,
     &                             CPAERR,ICPAFLAG,IQCPA,NDIMCLU,NQMAX,
     &                             NKMMAX)
C
      IF ( WRMAT ) THEN
         IWRI = 6
         IWRI = 100
         WRITE (IWRI,*) 'NLCPA'
         CALL NLCDUMPTAU(IECURR,ERYD,IWRI,MUEHAT,OMEGAHAT,TAUHAT,IREL,N,
     &                   NDIMCLU)
      END IF
C
      DEALLOCATE (OMEGAHAT,MMAT,TAUIMP,W1HAT,WGIT,EVD,EVDZJ,TAUII,W1)
      DEALLOCATE (MUEHAT,TAUHAT)
99001 FORMAT (/,5X,'SRO - parameter ALPHA           ',F6.3)
99002 FORMAT (//,1X,79('*'),/,25X,'NLCPA - coarse graining of BZ',/,1X,
     &        79('*'),//,5X,'number of tiles ',i5,/)
99003 FORMAT (/,'     ikn       kx        ky        kz',
     &        '      # of ks   nsym')
99004 FORMAT ('NLCPA:: too many atoms in cluster, function nlcpaconf ',
     &        'will fail on 32 bit machine')
99005 FORMAT (/,5x,'Number of irreducible tiles  :',i5,/)
      END
C*==nlcpaiter.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPAITER(CPAERR,P,IPRINT,NKM,NQ,NKMQ,TAUQ,DROT,
     &                     IQORGQP,SYMUNITARY,MSSQ,MSST,ITOQ,IREL,NSYM,
     &                     NTMAX,NQMAX,NKMMAX,NKTABMAX,MROTK,MUEHAT,
     &                     OMEGAHAT,TAUHAT,NDIMCLU,IND0QCLU,NQNLCPA,
     &                     NKTABKN,KTABKN,WKTABKN,NKN,NSYMACCEPTEDKN,
     &                     SYMACCEPTEDKN,KNNLCPA,RQNLCPA,IQCPA,CPATOL,
     &                     CPAMIX,NKNIRMU,SYM_IRR,MAXNKNIRMU,LCPAMILLS,
     &                     LNLCPASYMAVG,LNLCPAAVG,ITCPA,ITCPAMAX,
     &                     ICPAFLAG,ICPACONV,PCFG,NCFG,NAQ,NA_Q,NU,NAQU,
     &                     I0_QU,I1_QU,I2_QU,I1_Q,I2_Q)
C   ********************************************************************
C   *                                                                  *
C   *  DK, HE, January 2006                                            *
C   *                                                                  *
C   *   + Implementation of the NLCPA iterative algorithm to           *
C   *     determine effective NLCPA medium                             *
C   *                                                                  *
C   *                   after D. A. Rowlands, PRB 67, 115109 (2003)    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      USE MOD_CONSTANTS,ONLY:C0,CI,PI
      USE MOD_ANGMOM,ONLY:NKKR
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NLCPAITER')
C
C Dummy arguments
C
      REAL*8 CPAERR,CPAMIX,CPATOL
      INTEGER ICPACONV,ICPAFLAG,IPRINT,IQCPA,IREL,ITCPA,ITCPAMAX,
     &        MAXNKNIRMU,NAQ,NAQU,NCFG,NDIMCLU,NKM,NKMMAX,NKN,NKTABMAX,
     &        NQ,NQMAX,NQNLCPA,NSYM,NTMAX,NU
      LOGICAL LCPAMILLS,LNLCPAAVG,LNLCPASYMAVG
      COMPLEX*16 P
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),MUEHAT(NDIMCLU,NDIMCLU),
     &           OMEGAHAT(NDIMCLU,NDIMCLU),TAUHAT(NDIMCLU,NDIMCLU),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER I0_QU(NQ,NU),I1_Q(NQ),I1_QU(NQ,NU),I2_Q(NQ),I2_QU(NQ,NU),
     &        IND0QCLU(NQNLCPA),IQORGQP(NSYMMAX,NQMAX),ITOQ(NTMAX,NQMAX)
     &        ,NA_Q(NQ),NKMQ(NQMAX),NKNIRMU(NKN),NKTABKN(NKN),
     &        NSYMACCEPTEDKN(NKN)
      REAL*8 KNNLCPA(3,NQNLCPA),KTABKN(3,NKTABMAX,NKN),
     &       MROTK(3,3,NSYMMAX),PCFG(NCFG),RQNLCPA(3,NQNLCPA),
     &       WKTABKN(NKTABMAX,NKN)
      LOGICAL SYMACCEPTEDKN(NSYMMAX,NKN),SYMUNITARY(NSYMMAX),
     &        SYM_IRR(NKN,MAXNKNIRMU,NSYM)
C
C Local variables
C
      COMPLEX*16 CWGT,EIKNRIJ(:,:,:,:),ERRMILLS(:,:),MAUX(:,:),MMAT(:,:)
     &           ,MSSQKN(:,:,:),MUEHATKN(:,:),MUEHATMIX(:,:),SUMQ(:,:,:)
     &           ,TAUHATBZI(:,:),TAUIMP(:,:),TAUK(:,:),TAUKN(:,:),
     &           TAUKND(:,:),W1(:,:),W1HAT(:,:),W2HAT(:,:),XIMP(:,:)
      INTEGER I,I0,IA_ERR,ICALL,ICFG,IKN,IKND,IKNM,IOCC,IPRINTL,IQ,
     &        IQCLU,ISYM,IT,IU,J,J0,JQCLU,JU,M,NDIMCLUSQ
      REAL*8 KD(3),KD0(3),KNRIJ,RIJ
      INTEGER NLCPACONF
      SAVE EIKNRIJ
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE MMAT,MUEHATMIX,MSSQKN,MAUX,SUMQ,TAUKND
      ALLOCATABLE TAUIMP,MUEHATKN,TAUKN,W1HAT,W2HAT,W1,EIKNRIJ,TAUK
      ALLOCATABLE TAUHATBZI,XIMP,ERRMILLS
C
      ALLOCATE (TAUK(NKKR,NKKR),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUK MAUX')
C
      ALLOCATE (W1(NKMMAX,NKMMAX),SUMQ(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SUMQ W1')
C
      IPRINTL = IPRINT
C
C
      NU = NQNLCPA
C
C=======================================================================
      ICALL = ICALL + 1
C
C   ********************************************************************
C                   INITIALISATION - START
C   ********************************************************************
C
      IF ( ICALL.EQ.1 ) THEN
C
C        ------------------ set up the phase factors for coarse graining
         ALLOCATE (EIKNRIJ(NU,MAXNKNIRMU,NU,NU),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NLMQCLU')
C
         DO IKN = 1,NKN
C           --- take care of reducible kvecs obtained by rot of irr kvec
            DO IKNM = 1,NKNIRMU(IKN)
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
 10            CONTINUE
               DO JU = 1,NU
                  DO IU = 1,NU
                     KNRIJ = 0D0
                     DO I = 1,3
                        RIJ = RQNLCPA(I,IU) - RQNLCPA(I,JU)
                        KNRIJ = KNRIJ + KD(I)*RIJ
                     END DO
                     EIKNRIJ(IKN,IKNM,IU,JU) = EXP(CI*2*PI*KNRIJ)
                  END DO
               END DO
            END DO
         END DO
C-----------------------------------------------------------------------
      END IF
C **********************************************************************
C                   INITALISATION - END
C **********************************************************************
C
      ALLOCATE (MUEHATKN(NAQ,NAQ),TAUKN(NAQ,NAQ),TAUKND(NAQ,NAQ))
      TAUKN(1:NAQ,1:NAQ) = C0
C
      M = NKMMAX
      ALLOCATE (MSSQKN(M,M,NQMAX))
      CALL CINIT(NKMMAX*NKMMAX*NQMAX,MSSQKN)
C
      ALLOCATE (W1HAT(NAQU,NAQU),TAUHATBZI(NAQU,NAQU),XIMP(NAQU,NAQU))
      ALLOCATE (MUEHATMIX(NAQU,NAQU),MMAT(NAQU,NAQU),W2HAT(NAQU,NAQU))
      ALLOCATE (TAUIMP(NAQU,NAQU),ERRMILLS(NAQU,NAQU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: W1HAT')
C
      TAUIMP(1:NAQU,1:NAQU) = C0
      MMAT(1:NAQU,1:NAQU) = C0
C
      NDIMCLUSQ = NDIMCLU*NDIMCLU
C
C------------------------------------------------- initialize MUEHAT (1)
C
      MUEHAT(1:NAQU,1:NAQU) = C0
      MUEHATMIX(1:NAQU,1:NAQU) = C0
C
      DO IU = 1,NU
         DO IQ = 1,NQ
            I0 = I0_QU(IQ,IU)
C
            DO J = 1,NKM
               DO I = 1,NKM
                  MUEHAT(I0+I,I0+J) = MSSQ(I,J,IQ)
                  MUEHATMIX(I0+I,I0+J) = MSSQ(I,J,IQ)
               END DO
            END DO
         END DO
      END DO
C
C============================================================ NLCPA-LOOP
C
      ITCPA = 0
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      TAUHAT(1:NAQU,1:NAQU) = C0
C--------------------------------------------------------- KN-LOOP START
      IKND = 0
      DO IKN = 1,NKN
         IKND = IKND + 1
C
C---------------------------------------------------------- MUEHATKN (3)
C
         MUEHATKN(1:NAQ,1:NAQ) = C0
         DO JU = 1,NU
            J0 = I0_QU(1,JU)
            DO IU = 1,NU
               I0 = I0_QU(1,IU)
               CWGT = 1D0/EIKNRIJ(IKN,1,IU,JU)/DBLE(NU)
               DO J = 1,NAQ
                  CALL ZAXPY(NAQ,CWGT,MUEHAT(I0+1,J0+J),1,MUEHATKN(1,J),
     &                       1)
               END DO
            END DO
         END DO
C------------------------------------------------------------- TAUKN (4)
C----------------------------------------  do (BZ) integration over tile
         IF ( IPRINTL.GT.4 ) THEN
            WRITE (6,'("kn = ",4i5)') IKN
            CALL CMATSTRUCT('MUEHATKN',MUEHATKN,NKM,NAQ,0,0,0,1D-8,6)
         END IF
C
C-------------------------------------------------- add MUEHATKN to MSSQ
         DO IQ = 1,NQ
            M = NKMQ(IQ)
            MSSQKN(1:M,1:M,IQ) = MUEHATKN(I1_Q(IQ):I2_Q(IQ),I1_Q(IQ):
     &                           I2_Q(IQ))
         END DO
C
C##         write(6,'(3f10.6)') ktabkn(:,1:nktabkn(ikn),ikn)
         CALL NLCPABZINT(P,TAUQ,TAUK,MAUX,SUMQ,W1,DROT,IQORGQP,
     &                   SYMUNITARY,SYMACCEPTEDKN(1,IKN),MSSQKN,
     &                   WKTABKN(1,IKN),KTABKN(1,1,IKN),NKTABKN(IKN),
     &                   NSYM,NSYMACCEPTEDKN(IKN),NKTABMAX)
C
         IF ( IPRINTL.GT.3 ) THEN
            WRITE (6,'(70("-"))')
            WRITE (6,'("nlcpa:: TAU for ikn ",9i5)') IKN
            CALL CMATSTRUCT('TAUQ         ',TAUQ(1,1,1),NKM,NKMMAX,IREL,
     &                      IREL,0,1D-8,6)
         END IF
C
C-------------------------------------------------------------------- IQ
C        ------------------------ obtain tau's for symmetry related kn's
         DO IQ = 1,NQ
C
            TAUKND(1:NKM,1:NKM) = TAUQ(1:NKM,1:NKM,IQ)
            DO IKNM = 1,NKNIRMU(IKN)
               IF ( IKNM.GT.1 ) THEN
                  IKND = IKND + 1
C        ------------------------------------------- do rotation of tauq
                  DO ISYM = 1,NSYM
                     IF ( SYM_IRR(IKN,IKNM-1,ISYM) ) THEN
C                     CALL ROTATE(TAUKND(1,1),'G->L',TAUKN,NKM,
C     &                           DROT(1,1,ISYM),NKMMAX)
                        CALL ROTATETAUKN(TAUKND,TAUKN,W1,DROT,NKM,
     &                     NKMMAX,SYMUNITARY,ISYM)
                        IF ( IPRINTL.GT.3 ) THEN
                           WRITE (6,
     &                            '("nlcpa:: TAU for rotated ikn ",9i5)'
     &                            ) IKN
                           CALL CMATSTRUCT('TAUQ         ',TAUKN(1,1),
     &                        NKM,NKMMAX,IREL,IREL,0,1D-8,6)
                        END IF
                        EXIT
                     END IF
                  END DO
                  TAUQ(1:NKM,1:NKM,IQ) = TAUKN(1:NKM,1:NKM)
               END IF
C
               DO JQCLU = 1,NQNLCPA
                  J0 = IND0QCLU(JQCLU)
                  DO IQCLU = 1,NQNLCPA
                     I0 = IND0QCLU(IQCLU)
C
                     CWGT = EIKNRIJ(IKN,IKNM,IQCLU,JQCLU)/DBLE(NQNLCPA)
                     DO J = 1,NKM
                        CALL ZAXPY(NKM,CWGT,TAUQ(1,J,IQCPA),1,
     &                             TAUHAT(I0+1,J0+J),1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C-------------------------------------------------------------------- IQ
C
      END DO
C----------------------------------------------------------- KN-LOOP END
C
C----------------------- enforce symmetry, equalize site diagonal blocks
      IF ( LNLCPAAVG ) CALL NLCPAAVG(TAUHAT,NQNLCPA,IND0QCLU,DROT,
     &                               SYMUNITARY,NSYM,NKMMAX,NKM,IQCPA,
     &                               NDIMCLU,NQMAX,NKMQ,IQORGQP,
     &                               LNLCPASYMAVG)
C-----------------------------------------------------------------------
      IF ( IPRINTL.GT.3 ) THEN
         WRITE (6,'(70("-"))')
         WRITE (6,'("nlcpa:: tauhat (I,J)")')
         CALL CMATSTRUCT('TAUHAT',TAUHAT,NDIMCLU,NDIMCLU,0,0,0,1D-8,6)
      END IF
C
      M = NDIMCLU
      TAUHATBZI(1:M,1:M) = TAUHAT(1:M,1:M)
C
C---------------------------------------------------------- OMEGAHAT (5)
C
      CALL CMATINV(NDIMCLU,NDIMCLU,TAUHAT,W1HAT)
C
      M = NDIMCLU
      OMEGAHAT(1:M,1:M) = MUEHAT(1:M,1:M) - W1HAT(1:M,1:M)
C
      IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('OMEGAHAT',OMEGAHAT,NDIMCLU,
     &     NDIMCLU,0,0,0,1D-8,6)
C
C---------------------------------------------------------- imp-LOOP (6)
C
      CALL CINIT(NDIMCLUSQ,TAUHAT)
C
      M = NDIMCLU
C
      IF ( LCPAMILLS ) CALL CINIT(NDIMCLUSQ,ERRMILLS)
C
      DO ICFG = 1,NCFG
         IF ( IPRINTL.GT.2 ) WRITE (6,99001) ICFG
C        ------------------------- determine MMAT for configuration ICFG
         DO IQCLU = 1,NQNLCPA
            IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
            IT = ITOQ(IOCC,IQCPA)
            I0 = IND0QCLU(IQCLU)
            MMAT(I0+1:I0+NKM,I0+1:I0+NKM) = MSST(1:NKM,1:NKM,IT)
         END DO
C
         IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('mmat',MMAT,NDIMCLU,
     &        NDIMCLU,0,0,0,1D-8,6)
C        ---------------------------------------------------------------
         IF ( .NOT.LCPAMILLS ) THEN
            W1HAT(1:M,1:M) = MMAT(1:M,1:M) - OMEGAHAT(1:M,1:M)
C
            IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('W1HAT',W1HAT,NDIMCLU,
     &           NDIMCLU,0,0,0,1D-8,6)
C
            CALL CMATINV(NDIMCLU,NDIMCLU,W1HAT,TAUIMP)
C
            IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('TAUIMP',TAUIMP,NDIMCLU,
     &           NDIMCLU,0,0,0,1D-8,6)
C
            CWGT = PCFG(ICFG)
            CALL ZAXPY(NDIMCLUSQ,CWGT,TAUIMP,1,TAUHAT,1)
C        --------------------------------------------- prepare CPA-Mills
         ELSE
            W1HAT(1:M,1:M) = MMAT(1:M,1:M) - MUEHATMIX(1:M,1:M)
            CALL CMATINV(NDIMCLU,NDIMCLU,W1HAT,W2HAT)
            W1HAT(1:M,1:M) = W2HAT(1:M,1:M) + TAUHATBZI(1:M,1:M)
            CALL CMATINV(NDIMCLU,NDIMCLU,W1HAT,XIMP)
            ERRMILLS(1:M,1:M) = ERRMILLS(1:M,1:M) - PCFG(ICFG)
     &                          *XIMP(1:M,1:M)
         END IF
      END DO
C------------------------------------------------------------ MUEHAT (7)
C
      IF ( .NOT.LCPAMILLS ) THEN
         IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('TAUHAT',TAUHAT,NDIMCLU,
     &        NDIMCLU,0,0,0,1D-8,6)
C----------------------- enforce symmetry, equalize site diagonal blocks
C##       IF ( LNLCPAAVG ) then
C##          CALL NLCPAAVG(TAUHAT,NQNLCPA,IND0QCLU,DROT,SYMUNITARY,NSYM
C##   &           ,NKMMAX,NKM,IQCPA,NDIMCLU,NQMAX,NKMQ,IQORGQP
C##   &           ,LNLCPASYMAVG)
C##
C##          IF ( IPRINTL.GT.4 ) CALL  CMATSTR('TAUHAT_SYM',TAUHAT
C##   &           ,NDIMCLU,NDIMCLU,0,0,0,1D-8,6)
C##       end if
C
         W1HAT(:,:) = TAUHAT(:,:)
C
         CALL CMATINV(NDIMCLU,NDIMCLU,W1HAT,W2HAT)
         IF ( IPRINTL.GT.4 ) CALL CMATSTRUCT('w2hat',W2HAT,NDIMCLU,
     &        NDIMCLU,0,0,0,1D-8,6)
C
         M = NDIMCLU
         MUEHAT(1:M,1:M) = W2HAT(1:M,1:M) + OMEGAHAT(1:M,1:M)
C
         IF ( IPRINTL.GT.3 ) CALL CMATSTRUCT('m_NLCPA',MUEHAT,NDIMCLU,
     &        NDIMCLU,0,0,0,1D-8,6)
C
C        --------------------- enforce symmetry to stabilize convergence
C        --------------------- equalize site diagonal blocks
         IF ( LNLCPAAVG ) CALL NLCPAAVG(MUEHAT,NQNLCPA,IND0QCLU,DROT,
     &                                  SYMUNITARY,NSYM,NKMMAX,NKM,
     &                                  IQCPA,NDIMCLU,NQMAX,NKMQ,
     &                                  IQORGQP,LNLCPASYMAVG)
C
         IF ( IPRINTL.GT.3 ) CALL CMATSTRUCT('m_NLCPA_symm',MUEHAT,
     &        NDIMCLU,NDIMCLU,0,0,0,1D-8,6)
C
      END IF
C
C------------------------------------------------- check convergence (8)
C     implement mixing and error check
C----------------------------------------------------------------------
C
      IF ( LCPAMILLS ) THEN
C        ---------------------------------------------------------------
C        ------------------------------------------ CPA Mills alogorithm
C        ---------------------------------------------------------------
C
C        --------------- enforce symmetry, equalize site diagonal blocks
         IF ( LNLCPAAVG ) CALL NLCPAAVG(ERRMILLS,NQNLCPA,IND0QCLU,DROT,
     &                                  SYMUNITARY,NSYM,NKMMAX,NKM,
     &                                  IQCPA,NDIMCLU,NQMAX,NKMQ,
     &                                  IQORGQP,LNLCPASYMAVG)
C
         IF ( IPRINTL.GT.3 ) CALL CMATSTRUCT('errmills',ERRMILLS,
     &        NDIMCLU,NDIMCLU,0,0,0,1D-8,6)
C
         CALL CMATMUL(M,M,ERRMILLS,TAUHATBZI,W1HAT)
C
         DO I = 1,M
            W1HAT(I,I) = W1HAT(I,I) + 1D0
         END DO
C
         CALL CMATINV(NDIMCLU,NDIMCLU,W1HAT,W2HAT)
C
         CALL CMATMUL(M,M,W2HAT,ERRMILLS,W1HAT)
C
C        ---------------------------------------------------- MUEHAT (7)
         MUEHAT(1:M,1:M) = MUEHATMIX(1:M,1:M) - W1HAT(1:M,1:M)
C
         IF ( IPRINTL.GT.3 ) CALL CMATSTRUCT('w1hat',W1HAT,NDIMCLU,
     &        NDIMCLU,0,0,0,1D-8,6)
         CPAERR = 0D0
         DO J = 1,NDIMCLU
            DO I = 1,NDIMCLU
               CPAERR = MAX(CPAERR,ABS(MUEHAT(I,J)-MUEHATMIX(I,J)))
            END DO
         END DO
         MUEHATMIX(1:M,1:M) = MUEHAT(1:M,1:M)
C         WRITE (6,'("itcpa, nlcpaerr",i5,e20.6)') ITCPA,CPAERR
      ELSE
C        ---------------------------------------------------------------
C        ------------------------------------------------- simple mixing
C        ---------------------------------------------------------------
         CPAERR = 0D0
         DO J = 1,NDIMCLU
            DO I = 1,NDIMCLU
               CPAERR = MAX(CPAERR,ABS(MUEHAT(I,J)-MUEHATMIX(I,J)))
            END DO
         END DO
         WRITE (6,'("itcpa, nlcpaerr",i5,e20.6)') ITCPA,CPAERR
C        ---------------------- do mixing on inverse of cluster t-matrix
         M = NDIMCLU
         MUEHAT(1:M,1:M) = MUEHATMIX(1:M,1:M)*(1D0-CPAMIX)
     &                     + MUEHAT(1:M,1:M)*CPAMIX
C
         MUEHATMIX(1:M,1:M) = MUEHAT(1:M,1:M)
      END IF
C
      IF ( IPRINT.GE.1 ) WRITE (6,99004) CPAERR
C
      IF ( CPAERR.LE.CPATOL ) THEN
         ICPACONV = 1
         IF ( IPRINTL.GT.0 ) WRITE (6,99002) ITCPA,CPAERR
      ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
         WRITE (6,99003) ITCPA,CPAERR
         ICPAFLAG = 1
      ELSE
         GOTO 100
      END IF
C======================================================  END  NLCPA-LOOP
C
      M = NDIMCLU
      TAUHAT(1:M,1:M) = TAUHATBZI(1:M,1:M)
C
      DEALLOCATE (MMAT,MUEHATMIX,MSSQKN,MAUX,SUMQ,TAUKND)
      DEALLOCATE (TAUIMP,MUEHATKN,TAUKN,W1HAT,W2HAT,W1,TAUK)
      DEALLOCATE (TAUHATBZI,XIMP,ERRMILLS)
C
99001 FORMAT (5x,'nlcpaiter:: imp-loop : icfg',i5)
99002 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,:,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99003 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,:,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
C99004FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
C    &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,:,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
C99006 FORMAT ('E ',2F10.5,' CPA ',I5,3E12.5)
      END
C*==nlcpaavg.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPAAVG(MATINOUT,NQNLCPA,IND0QCLU,DROT,SYMUNITARY,
     &                    NSYM,NKMMAX,NKM,IQCPA,NDIMCLU,NQMAX,NKMQ,
     &                    IQORGQP,LNLCPASYMAVG)
C   ********************************************************************
C   *                                                                  *
C   *  DK, HE, February                                                *
C   *                                                                  *
C   *    + average the site diagonal blocks of the cluster matrices    *
C   *    + optionally symmetrise                                       *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQCPA,NDIMCLU,NKM,NKMMAX,NQMAX,NQNLCPA,NSYM
      LOGICAL LNLCPASYMAVG
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),MATINOUT(NDIMCLU,NDIMCLU)
      INTEGER IND0QCLU(NQNLCPA),IQORGQP(NSYMMAX,NQMAX),NKMQ(NQMAX)
      LOGICAL SYMUNITARY(NSYMMAX)
C
C Local variables
C
      INTEGER I0,IQCLU,IQORG(NSYMMAX),ISYM,N,NSYMACCEPTEDAVG
      LOGICAL SYMACCEPTEDAVG(NSYMMAX)
      COMPLEX*16 WR1(:,:),WR2(:,:),WR3(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WR1,WR2,WR3
C
      ALLOCATE (WR1(NKMMAX,NKMMAX),WR2(NKMMAX,NKMMAX))
      ALLOCATE (WR3(NKMMAX,NKMMAX))
C
C     - take into account possibilty of having several CPA sites in cell
      NSYMACCEPTEDAVG = 0
      DO ISYM = 1,NSYM
         IF ( IQCPA.EQ.IQORGQP(ISYM,IQCPA) ) THEN
            SYMACCEPTEDAVG(ISYM) = .TRUE.
            NSYMACCEPTEDAVG = NSYMACCEPTEDAVG + 1
         ELSE
            SYMACCEPTEDAVG(ISYM) = .FALSE.
         END IF
      END DO
C     -------------------------------------- sum up site diagonal blocks
      N = NKM
      WR1(1:N,1:N) = (0D0,0D0)
      DO IQCLU = 1,NQNLCPA
         I0 = IND0QCLU(IQCLU)
         WR1(1:N,1:N) = WR1(1:N,1:N) + MATINOUT(I0+1:I0+N,I0+1:I0+N)
      END DO
C
      IQORG(1:NSYMMAX) = IQORGQP(1:NSYMMAX,IQCPA)
C      CALL  CMATSTR('wr1',wr1,nkmmax,nkmmax,0,0,0,1D-8,6)
C     ------------------------------------------------------- symmetrise
      IF ( LNLCPASYMAVG ) CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,
     &     DBLE(NQNLCPA),WR1,WR2,WR3,1,NKMQ(IQCPA),DROT,IQORG,
     &     SYMUNITARY,SYMACCEPTEDAVG,NSYM,NSYMACCEPTEDAVG,1,NKMMAX)
C     -------------------------- put blocks back onto site-diagonal part
      N = NKM
      DO IQCLU = 1,NQNLCPA
         I0 = IND0QCLU(IQCLU)
         IF ( .NOT.LNLCPASYMAVG ) THEN
            MATINOUT(I0+1:I0+N,I0+1:I0+N) = WR1(1:N,1:N)/DBLE(NQNLCPA)
         ELSE
            MATINOUT(I0+1:I0+N,I0+1:I0+N) = WR2(1:N,1:N)
         END IF
      END DO
C
      DEALLOCATE (WR1,WR2,WR3)
C
      END
C*==nlcpabzint.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPABZINT(P,TAUQ,TAUK,MAUX,SUMQ,W1,DROT,IQORGQP,
     &                      SYMUNITARY,SYMACCEPTED,MSSQ,WKTAB,KTAB,
     &                      NKTAB,NSYM,NSYMACCEPTED,NKTABMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SITES,ONLY:NQ,NQMAX
      USE MOD_ANGMOM,ONLY:NKKR,NKMQ,IND0Q,NKMMAX
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKTAB,NKTABMAX,NSYM,NSYMACCEPTED
      COMPLEX*16 P
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),MAUX(NKKR,NKKR),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),SUMQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUK(NKKR,NKKR),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           W1(NKMMAX,NKMMAX)
      INTEGER IQORGQP(NSYMMAX,NQMAX)
      REAL*8 KTAB(3,NKTABMAX),WKTAB(NKTABMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      COMPLEX*16 CWK
      INTEGER I1,IK,IQ,J,J1,N
      REAL*8 WK,WKSUM
C
C*** End of declarations rewritten by SPAG
C
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
      CALL CINIT(NKMMAX*NKMMAX*NQ,SUMQ)
C
      WKSUM = 0D0
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      DO IK = 1,NKTAB
C
         CALL STRSET(IK,KTAB(1,IK),TAUK,MAUX,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,MAUX,MSSQ,NQMAX,NKKR,NKMMAX)
C
         CALL CINVLU(MAUX,TAUK,NKKR,NKKR)
C
C------------------------------------------------------------ store TAUQ
C
         WK = WKTAB(IK)
         WKSUM = WKSUM + WK
         CWK = DCMPLX(WK,0D0)
         DO IQ = 1,NQ
            I1 = IND0Q(IQ) + 1
            N = NKMQ(IQ)
            DO J = 1,N
               J1 = IND0Q(IQ) + J
               CALL ZAXPY(N,CWK,TAUK(I1,J1),1,SUMQ(1,J,IQ),1)
            END DO
         END DO
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,WKSUM,SUMQ,TAUQ,W1,NQ,NKMQ,
     &               DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM,
     &               NSYMACCEPTED,NQMAX,NKMMAX)
C
      END
C*==nlcpaset.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPASET(BRAVAIS,CPALVL,NQNLCPA,RQNLCPA,KNNLCPA)
C   ********************************************************************
C   *                                                                  *
C   *  get the NLCPA cluster data                                      *
C   *   RQNLCPA  site vectors in (a)                                   *
C   *   KNNLCPA  K_n vectors  in (2PI/a)                               *
C   *                                                                  *
C   *  for the cubic lattices see: Rowlands, PRB 67, 115109 (2003)     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BRAVAIS,CPALVL,NQNLCPA
      REAL*8 KNNLCPA(3,NQNLCPA),RQNLCPA(3,NQNLCPA)
C
C Local variables
C
      INTEGER I,ISC
      REAL*8 KSC(3,8),RSC(3,8),V000(3),V001(3),V010(3),V0HH(3),V100(3),
     &       VH0H(3),VHH0(3),VHHH(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA V000/0D0,0D0,0D0/,V100/1D0,0D0,0D0/,V010/0D0,1D0,0D0/
      DATA V001/0D0,0D0,1D0/
      DATA VHHH/0.5D0,0.5D0,0.5D0/,V0HH/0D0,0.5D0,0.5D0/
      DATA VH0H/0.5D0,0D0,0.5D0/,VHH0/0.5D0,0.5D0,0D0/
C
      DATA RSC/0D0,0D0,0D0,0D0,0D0,1D0,0D0,1D0,0D0,0D0,1D0,1D0,1D0,0D0,
     &     0D0,1D0,0D0,1D0,1D0,1D0,0D0,1D0,1D0,1D0/
      DATA KSC/0D0,0D0,0D0,0D0,0D0,0.5D0,0D0,0.5D0,0D0,0D0,0.5D0,0.5D0,
     &     0.5D0,0D0,0D0,0.5D0,0D0,0.5D0,0.5D0,0.5D0,0D0,0.5D0,0.5D0,
     &     0.5D0/
C
C------------------------------------- recovering the CPA ( CPALVL = 1 )
      IF ( CPALVL.EQ.1 ) THEN
         KNNLCPA = 0D0
         RQNLCPA = 0D0
         RETURN
      END IF
C-------------------------------------------------------------------- SC
      IF ( BRAVAIS.EQ.12 ) THEN
         IF ( CPALVL.LE.3 ) THEN
            CALL RVECCOP(3*NQNLCPA,RSC,RQNLCPA)
            CALL RVECCOP(3*NQNLCPA,KSC,KNNLCPA)
         END IF
C------------------------------------------------------------------- FCC
      ELSE IF ( BRAVAIS.EQ.13 ) THEN
         IF ( CPALVL.EQ.2 ) THEN
            CALL RVECCOP(3,V000,RQNLCPA(1,1))
            CALL RVECCOP(3,VH0H,RQNLCPA(1,2))
            CALL RVECCOP(3,VHH0,RQNLCPA(1,3))
            CALL RVECCOP(3,V0HH,RQNLCPA(1,4))
C
            CALL RVECCOP(3,V000,KNNLCPA(1,1))
            CALL RVECCOP(3,V100,KNNLCPA(1,2))
            CALL RVECCOP(3,V010,KNNLCPA(1,3))
            CALL RVECCOP(3,V001,KNNLCPA(1,4))
         ELSE IF ( CPALVL.EQ.3 ) THEN
            DO ISC = 1,8
               DO I = 1,3
                  RQNLCPA(I,ISC+0) = V000(I) + RSC(I,ISC)
                  RQNLCPA(I,ISC+8) = VH0H(I) + RSC(I,ISC)
                  RQNLCPA(I,ISC+16) = VHH0(I) + RSC(I,ISC)
                  RQNLCPA(I,ISC+24) = V0HH(I) + RSC(I,ISC)
C
                  KNNLCPA(I,ISC+0) = V000(I) + KSC(I,ISC)
                  KNNLCPA(I,ISC+8) = V100(I) + KSC(I,ISC)
                  KNNLCPA(I,ISC+16) = V010(I) + KSC(I,ISC)
                  KNNLCPA(I,ISC+24) = V001(I) + KSC(I,ISC)
               END DO
            END DO
         END IF
C------------------------------------------------------------------- BCC
      ELSE IF ( BRAVAIS.EQ.14 ) THEN
         IF ( CPALVL.EQ.2 ) THEN
            CALL RVECCOP(3,V000,RQNLCPA(1,1))
            CALL RVECCOP(3,VHHH,RQNLCPA(1,2))
C
            CALL RVECCOP(3,V000,KNNLCPA(1,1))
            CALL RVECCOP(3,V100,KNNLCPA(1,2))
         ELSE IF ( CPALVL.EQ.3 ) THEN
            DO ISC = 1,8
               DO I = 1,3
                  RQNLCPA(I,ISC+0) = V000(I) + RSC(I,ISC)
                  RQNLCPA(I,ISC+8) = VHHH(I) + RSC(I,ISC)
C
                  KNNLCPA(I,ISC+0) = V000(I) + KSC(I,ISC)
                  KNNLCPA(I,ISC+8) = V100(I) + KSC(I,ISC)
               END DO
            END DO
         END IF
      END IF
C
      END
C*==reducekn.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE REDUCEKN(NKN,KNNLCPA,SYMACCEPTED,MROTK,NSYM,IPRINT,
     &                    IOTMP,LIRREDKN,ABAS)
C   ********************************************************************
C   *  DK, HE, January 2006                                            *
C   *                                                                  *
C   *  Reduce the set of Kn vectors pointing to the tiles to a set     *
C   *  of irreducible Kn's                                             *
C   *                                                                  *
C   *  keep track of the rotaions for the BZ integration               *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='REDUCEKN')
C
C Dummy arguments
C
      INTEGER IOTMP,IPRINT,NKN,NSYM
      LOGICAL LIRREDKN
      REAL*8 ABAS(3,3),KNNLCPA(3,NKN),MROTK(3,3,NSYMMAX)
      LOGICAL SYMACCEPTED(NSYMMAX)
C
C Local variables
C
      LOGICAL EVER_USED,KNREDUCIBLE(:),SYMMETRY_REL(:,:,:),
     &        SYM_IRR(:,:,:)
      INTEGER I,IKN1,IKN2,IRREDKN,ISYM,NKNIRMU(:),STARVECS(:)
      REAL*8 KD(3),KD1(3),KNNLCPAIRR(:,:)
      LOGICAL KVECTEQ
C
C*** End of declarations rewritten by SPAG
C
C
C*** Start of declarations rewritten by SPAG
C PARAMETER definitions
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KNREDUCIBLE,SYMMETRY_REL,STARVECS
      ALLOCATABLE KNNLCPAIRR,SYM_IRR,NKNIRMU
C
      ALLOCATE (KNREDUCIBLE(NKN),STARVECS(NKN))
      ALLOCATE (SYMMETRY_REL(NKN,NKN,NSYM))
C
      IF ( .NOT.LIRREDKN ) THEN
         ALLOCATE (NKNIRMU(NKN))
         NKNIRMU(1:NKN) = 1
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
         WRITE (IOTMP) (NKNIRMU(IKN1),IKN1=1,NKN)
         RETURN
      END IF
C
      STARVECS(1:NKN) = 0
      SYMMETRY_REL = .FALSE.
      KNREDUCIBLE = .FALSE.
      IRREDKN = 1
C
      IF ( IPRINT.GT.4 ) THEN
         WRITE (6,'(/,5x,"k1 is transformed into k2 by sym - trafo")')
         WRITE (6,'(8x,"k1   k2   sym-trafo")')
         WRITE (6,'(8x,20("-"))')
      END IF
C
      DO IKN1 = 1,NKN - 1
         KD1 = KNNLCPA(1:3,IKN1)
         IRREDKN = IRREDKN + 1
         DO IKN2 = IKN1 + 1,NKN
            IF ( .NOT.KNREDUCIBLE(IKN2) ) THEN
               EVER_USED = .FALSE.
               DO ISYM = 1,NSYM
                  IF ( SYMACCEPTED(ISYM) ) THEN
                     CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,
     &                          KNNLCPA(1,IKN2),1,0D0,KD,1)
                     IF ( KVECTEQ(KD,KD1,ABAS) ) THEN
C                        WRITE (6,'(2(3f8.3,"|"),l1)') kd,kd1
C     &                       ,kvecteq(ABAS,kd,kd1)
                        EVER_USED = EVER_USED .OR. .TRUE.
                        IF ( IPRINT.GT.4 ) WRITE (6,'(5x,6i5)') IKN1,
     &                       IKN2,ISYM
                        SYMMETRY_REL(IKN1,IKN2,ISYM) = .TRUE.
                     END IF
                  END IF
               END DO
               IF ( EVER_USED ) THEN
                  IRREDKN = IRREDKN - 1
                  KNREDUCIBLE(IKN2) = .TRUE.
                  STARVECS(IKN1) = STARVECS(IKN1) + 1
               END IF
            END IF
         END DO
      END DO
C
      ALLOCATE (KNNLCPAIRR(3,IRREDKN))
      I = MAXVAL(STARVECS(1:NKN))
      ALLOCATE (SYM_IRR(IRREDKN,I,NSYM))
      ALLOCATE (NKNIRMU(IRREDKN))
C
C     ---------------------- compress information about irreducible kn's
      NKNIRMU(1:IRREDKN) = 1
C
      IRREDKN = 0
      DO IKN1 = 1,NKN
         IF ( .NOT.KNREDUCIBLE(IKN1) ) THEN
            IRREDKN = IRREDKN + 1
            KNNLCPAIRR(1:3,IRREDKN) = KNNLCPA(1:3,IKN1)
            IF ( STARVECS(IKN1).GT.0 ) THEN
               I = 1
               DO IKN2 = IKN1 + 1,NKN
                  IF ( ANY(SYMMETRY_REL(IKN1,IKN2,1:NSYM)) ) THEN
                     NKNIRMU(IRREDKN) = NKNIRMU(IRREDKN) + 1
                     SYM_IRR(IRREDKN,I,1:NSYM)
     &                  = SYMMETRY_REL(IKN1,IKN2,1:NSYM)
                     I = I + 1
                  END IF
               END DO
            END IF
         END IF
      END DO
C
      WRITE (6,*)
      WRITE (6,'(5x,"knir multpl. #   kvec / symm trafo  ")')
      WRITE (6,'(5x,65("-"))')
      DO IKN1 = 1,IRREDKN
         DO I = 1,NKNIRMU(IKN1)
            IF ( I.EQ.1 ) THEN
               WRITE (6,'(5x,3(i3," |"),3f8.4)') IKN1,NKNIRMU(IKN1),I,
     &                KNNLCPAIRR(1:3,IKN1)
            ELSE
               WRITE (6,'(5x,2(4x,"|"),i3," |  ",48l1)') I,
     &                SYM_IRR(IKN1,I-1,1:NSYM)
            END IF
         END DO
      END DO
C     ------------------------- copy knnlcpairr into knnlcpa, update nkn
C     ----------------- keep info about multiplicity and symmetry trafos
C     ------- communicate the latter info via hard drive to main program
      NKN = IRREDKN
      DO IKN1 = 1,IRREDKN
         KNNLCPA(1:3,IKN1) = KNNLCPAIRR(1:3,IKN1)
      END DO
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      WRITE (IOTMP) (NKNIRMU(IKN1),IKN1=1,IRREDKN)
      DO IKN1 = 1,IRREDKN
         IF ( NKNIRMU(IKN1).GT.1 ) WRITE (IOTMP)
     &        (SYM_IRR(IKN1,I,1:NSYM),I=1,NKNIRMU(IKN1)-1)
      END DO
C
      DEALLOCATE (KNREDUCIBLE,SYMMETRY_REL,STARVECS)
      DEALLOCATE (KNNLCPAIRR,SYM_IRR,NKNIRMU)
C
      END
C*==nlcpaconf.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      INTEGER FUNCTION NLCPACONF(ICONF,NC,ISITE,INDEX1,INDEX2)
C   ********************************************************************
C   *  DK, January 2006                                                *
C   *                                                                  *
C   *   nlcpaconf.f -                                                  *
C   *                                                                  *
C   *   description - returns the occupation (index1 or index2) of a   *
C   *                 site isite (isite = 1, ..nc) of a cluster        *
C   *                 consisting of nc atoms                           *
C   *               - uses intrinsic function btest for bit testing    *
C   *               - for cluster sizes larger than 31 this            *
C   *                 routine will fail on a 32 bit machine!           *
C   *                                                                  *
C   *   created - Tue Jan 10 13:21:46 CET 2006                         *
C   *   author  - Diemo Koedderitzsch <dk@eeb05>                       *
C   *                                                                  *
C   *  in                                                              *
C   *  --                                                              *
C   *  nc     . . .  number of atoms in cluster                        *
C   *  iconf  . . .  number of the configuration,iconf = 1,2, .. 2**nc *
C   *  isite  . . .  cluster site                                      *
C   *  index# . . .  number denoting constituent (A or B type atom)    *
C   *                                                                  *
C   *  out                                                             *
C   *  ---                                                             *
C   *  type (index1 or index2) which occupies isite in cluster         *
C   *  of size nc                                                      *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NLCPACONF')
C
C Dummy arguments
C
      INTEGER ICONF,INDEX1,INDEX2,ISITE,NC
C
C*** End of declarations rewritten by SPAG
C
      IF ( ICONF.GT.2**NC ) CALL STOP_MESSAGE(ROUTINE,
     &     'ICONF higher than # of allowed configurations!')
C
      NLCPACONF = INDEX1
      IF ( BTEST(ICONF-1,ISITE-1) ) NLCPACONF = INDEX2
C
C##C ----------------------------------------------------    test nclpaconf
C##C     j = # of sites in cluster, should be smaller than 31 for 32 bit
C##C     machine, because 2**31 already too large to be computed
C##      IF ( IPRINTL.GT.3 ) THEN
C##         J = 4
C##         IF ( J.GT.BIT_SIZE(J)-2 ) STOP
C##     &                        'nlcpa:: nqcluster higher than bit size !'
C##C
C##         DO K = 1,2**J
C##            WRITE (6,'("nlcpa:: nlcpaconf::", 99i5)') K,
C##     &             (NLCPACONF(K,J,I,1,2),I=1,J)
C##         END DO
C##      END IF
C##C ----------------------------------------------------    test nclpaconf
      END
C*==nlcpakmesh.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPAKMESH(BRAVAIS,CPALVL,BBAS,NKTABKN,NKN,IPRINT,
     &                      MROTK,LIRK,NSYM,SYMACCEPTEDKN,
     &                      NSYMACCEPTEDKN,KNNLCPA,IOTMP,NKPTS00,NKRED)
C   ********************************************************************
C   *                                                                  *
C   *  DK, January 2006                                                *
C   *                                                                  *
C   *  + determine kpoints for coarse grained tiles                    *
C   *  + reduce them to irreducible set if LIRK is true                *
C   *  + dump sets into file, write a Gnuplot driver file for          *
C   *    visualisation                                                 *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NLCPAKMESH')
C
C Dummy arguments
C
      INTEGER BRAVAIS,CPALVL,IOTMP,IPRINT,NKN,NKPTS00,NKRED,NSYM
      LOGICAL LIRK
      REAL*8 BBAS(3,3),KNNLCPA(3,NKN),MROTK(3,3,NSYMMAX)
      INTEGER NKTABKN(NKN),NSYMACCEPTEDKN(NKN)
      LOGICAL SYMACCEPTEDKN(NSYMMAX,NKN)
C
C Local variables
C
      LOGICAL ALREADYASS,SYM_EVER_ACC(:)
      REAL*8 CSM(3,3),CVEC(3,3),KD(3),KD1(3),KDL1,KTABD(:,:),
     &       KTABKN(:,:,:),OFFSET(3),TOL,TVEC(3),WK,WKTABKN(:,:),
     &       WTABD(:)
      CHARACTER*50 CSTR1
      CHARACTER*4000 CSTR2
      CHARACTER*100 CSTR3
      INTEGER DIVISIONS,I,I1,I2,I3,IK,IK1,IKN,ILS,IP,ISYM,J,NIK,
     &        NKTABKNMAX
      CHARACTER*15 KPTILEFILE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KTABD,WTABD,SYM_EVER_ACC,KTABKN,WKTABKN
C
      NKTABKNMAX = INT(DBLE(NKPTS00)*1.1D0)
      DIVISIONS = INT(DBLE(NKTABKNMAX)**(1D0/3D0))
      IF ( MOD(DIVISIONS,2).EQ.1 ) DIVISIONS = DIVISIONS + 1
      NKTABKNMAX = DIVISIONS**3
      WRITE (6,99004) NKTABKNMAX,DIVISIONS
      ALLOCATE (KTABKN(3,NKTABKNMAX,NKN),WKTABKN(NKTABKNMAX,NKN))
C
C-----------------------------------------------------------------------
C  DETERMINE KTABKN FOR SPECIAL CASE OF HAVING BCC or FCC
C-----------------------------------------------------------------------
C
      TVEC(1) = 0.5D0
      TVEC(2) = 0.5D0
      TVEC(3) = 0.5D0
C
      CVEC(1:3,1:3) = 0D0
      IF ( BRAVAIS.EQ.13 .OR. BRAVAIS.EQ.14 ) THEN
         DO I = 1,3
            CVEC(I,I) = 1D0
         END DO
         IF ( CPALVL.EQ.3 ) THEN
            CVEC(1:3,1:3) = CVEC(1:3,1:3)/2D0
            TVEC(1:3) = TVEC(1:3)/2D0
         END IF
      ELSE
         WRITE (6,*) 'nlcpacgkmesh:: lattice or cpalvl not impl. yet'
      END IF
C      -------------------------------------------------- has to be even
      TOL = 1D0/DBLE(100*DIVISIONS)
      CSM(1:3,1:3) = CVEC(1:3,1:3)/DBLE(DIVISIONS)
      OFFSET(1:3) = (CSM(1:3,1)+CSM(1:3,2)+CSM(1:3,3))/2D0
C
      KTABKN = 0D0
C
      NKRED = 0
C
      DO I1 = 0,DIVISIONS - 1
         DO I2 = 0,DIVISIONS - 1
            DO I3 = 0,DIVISIONS - 1
               NKRED = NKRED + 1
               KTABKN(1:3,NKRED,1) = I1*CSM(1:3,1) + I2*CSM(1:3,2)
     &                               + I3*CSM(1:3,3)
            END DO
         END DO
      END DO
C     -------------------------------------  Gamma points causes trouble
C##?? ---  because of shift this is not important !!!
C      KTABKN(1:3,1,1) = 1D-8
C
      DO J = 1,NKRED
         KTABKN(:,J,1) = KTABKN(:,J,1) + OFFSET(:) - TVEC(:)
      END DO
C## testing - not yet general
      IF ( CPALVL.EQ.1 ) THEN
         DO IK = 1,NKRED
            KTABKN(1:3,IK+NKRED,1) = KTABKN(1:3,IK,1) + KNNLCPA(1:3,2)
         END DO
         WKTABKN = 1D0
         NKTABKN = 2*NKRED
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'IBZ_bcc_1_cube.dat')
         DO IK = 1,2*NKRED
            WRITE (IOTMP,'(9F10.6)') (KTABKN(J,IK,1),J=1,3)
         END DO
         CLOSE (IOTMP)
         RETURN
      END IF
C
C     ----- produce kpoints for other tiles by adding transl vec kn(ikn)
      DO IKN = 2,NKN
         DO IK = 1,NKRED
            KTABKN(1:3,IK,IKN) = KTABKN(1:3,IK,1) + KNNLCPA(1:3,IKN)
         END DO
      END DO
C     -------------------------------- dump reducible kpoints into files
      IF ( IPRINT.GT.4 ) THEN
         KPTILEFILE = 'tiles_red.dat'
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,KPTILEFILE)
         DO IKN = 1,NKN
            DO IK = 1,NKRED
               WRITE (IOTMP,'(9f10.6)') (KTABKN(J,IK,IKN),J=1,3)
            END DO
            WRITE (IOTMP,'("#")')
            WRITE (IOTMP,'(" ")')
         END DO
         CLOSE (IOTMP)
      END IF
C
      WRITE (6,'(5x,"Number of reducible kvectors :",i5,/)') NKRED
      WKTABKN = 1D0
      NKTABKN = NKRED
C     ------------------------------------ reduce # of kpoints to irr kp
      IF ( LIRK ) THEN
C        ------------------ reduce the number of k-vecs for central tile
         ALLOCATE (KTABD(3,NKTABKNMAX))
         ALLOCATE (WTABD(NKTABKNMAX))
         ALLOCATE (SYM_EVER_ACC(NSYM))
         SYM_EVER_ACC(1:NSYM) = .TRUE.
C
         IF ( IPRINT.GT.4 ) THEN
            KPTILEFILE = 'tiles_irr.dat'
            CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,KPTILEFILE)
         END IF
C
         WRITE (6,'("     ikn  sum of wgts  symm trafos")')
         WRITE (6,'(5x,65("-"))')
C        BEGIN RTRTRTRTRTRTRTRT -------------- loop over reducible tiles
         DO IKN = 1,NKN
CC##            WRITE (6,'("kn  : ",3f10.6)') (KNNLCPA(I,IKN),I=1,3)
            NIK = 1
            KTABD(1:3,1) = KTABKN(1:3,1,IKN)
            WTABD(1) = 1D0
            SYM_EVER_ACC(1:NSYM) = .FALSE.
C           E is always accepted
            SYM_EVER_ACC(1) = .TRUE.
C           kkkkkkkkkkkkkkkkkkkkkkkkkk----- loop over reducible k-points
            DO IK = 2,NKRED
               ALREADYASS = .FALSE.
               KD1(1:3) = KTABKN(1:3,IK,IKN)
               KDL1 = DOT_PRODUCT(KD1,KD1)
C              ---------------------------------- loop over irr k-points
               DO IK1 = 1,NIK
                  KD(1:3) = KTABD(1:3,IK1)
                  IF ( ABS(DOT_PRODUCT(KD,KD)-KDL1).LT.TOL ) THEN
                     DO ISYM = 2,NSYM
                        IF ( SYMACCEPTEDKN(ISYM,IKN) ) THEN
                           CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,
     &                                KTABD(1,IK1),1,0D0,KD,1)
                           KD(1:3) = KD(1:3) - KD1(1:3)
                           IF ( DOT_PRODUCT(KD,KD).LT.1D-6 ) THEN
                              SYM_EVER_ACC(ISYM) = SYM_EVER_ACC(ISYM)
     &                           .OR. .TRUE.
                              WTABD(IK1) = WTABD(IK1) + 1D0
                              ALREADYASS = .TRUE.
                              GOTO 10
                           END IF
                        END IF
                     END DO
                  END IF
               END DO
 10            CONTINUE
               IF ( .NOT.ALREADYASS ) THEN
                  NIK = NIK + 1
                  KTABD(1:3,NIK) = KTABKN(1:3,IK,IKN)
                  WTABD(NIK) = 1D0
               END IF
            END DO
C           kkkkkkkkkkkkkkkkkkkkkkkkkk----- loop over reducible k-points
            KTABKN(1:3,1:NIK,IKN) = KTABD(1:3,1:NIK)
            WKTABKN(1:NIK,IKN) = WTABD(1:NIK)
            NKTABKN(IKN) = NIK
            SYMACCEPTEDKN(1:NSYM,IKN) = SYM_EVER_ACC(1:NSYM)
C     WRITE (6,*) wktabkn(1:nik,1 )
            NSYMACCEPTEDKN(IKN) = 0
            DO ISYM = 1,NSYM
               IF ( SYMACCEPTEDKN(ISYM,IKN) ) NSYMACCEPTEDKN(IKN)
     &              = NSYMACCEPTEDKN(IKN) + 1
            END DO
            WK = 0D0
            DO I = 1,NIK
               WK = WK + WKTABKN(I,IKN)
            END DO
            WRITE (6,'(2x,i5,f10.3,5x,48l1)') IKN,WK,SYM_EVER_ACC
C
            IF ( IPRINT.GT.4 ) THEN
               DO IK = 1,NIK
                  WRITE (IOTMP,'(9f10.6)') (KTABKN(J,IK,IKN),J=1,3)
               END DO
               WRITE (IOTMP,'("#")')
               WRITE (IOTMP,'(" ")')
            END IF
         END DO
         IF ( IPRINT.GT.4 ) CLOSE (IOTMP)
C        END RTRTRTRTRTRTRTRTRTR ------------- loop over reducible tiles
      END IF
C
C     ----------------------------------------------- write gnuplot file
      IF ( IPRINT.GT.4 ) THEN
         OPEN (UNIT=IOTMP,FILE='plot_kp.gp',STATUS='replace',
     &         FORM='formatted')
C
         WRITE (IOTMP,99001)
C     -------------------------- put reciprocal space  unit cell vectors
         CSTR3 = 'set arrow    lw 3 lt 3 from'
         WRITE (CSTR3(11:12),'(i2)') 1
         WRITE (CSTR3(30:100),99002) (BBAS(1,I),I=1,3)
         WRITE (IOTMP,'(a100)') CSTR3
         WRITE (CSTR3(11:12),'(i2)') 2
         WRITE (CSTR3(30:100),99002) (BBAS(2,I),I=1,3)
         WRITE (IOTMP,'(a100)') CSTR3
         WRITE (CSTR3(11:12),'(i2)') 3
         WRITE (CSTR3(30:100),99002) (BBAS(3,I),I=1,3)
         WRITE (IOTMP,'(a100)') CSTR3
C     --------------------------------------------- label the kn vectors
         CSTR3 = 'set label "     "  at  '
         WRITE (CSTR3(63:69),'(a6)') 'center'
         WRITE (CSTR3(70:90),'(a20)') 'font "Helvetica,24"'
         DO IKN = 1,NKN
            WRITE (CSTR3(13:15),'(i2)') IKN
            WRITE (CSTR3(30:62),'(2(f10.6,","),f10.6)') KNNLCPA(1:3,IKN)
            WRITE (IOTMP,'(a100)') CSTR3
         END DO
         WRITE (IOTMP,'("sp 0 notitle")')
         WRITE (IOTMP,99003)
C     ---------------------------------------- plotting the irr k points
         CSTR1 = ',''tiles_irr.dat'' every :::  ::   w p ps 3 t "  "'
         ILS = 50
         CSTR2 = 'splot'
         IP = 7
         DO IKN = 1,NKN
            IF ( IKN.GT.09 ) WRITE (CSTR1(27:28),'(i2)') IKN - 1
            IF ( IKN.LT.10 ) WRITE (CSTR1(27:27),'(i1)') IKN - 1
            IF ( IKN.GT.09 ) WRITE (CSTR1(31:32),'(i2)') IKN - 1
            IF ( IKN.LT.10 ) WRITE (CSTR1(31:31),'(i1)') IKN - 1
            IF ( IKN.GT.09 ) WRITE (CSTR1(46:47),'(i2)') IKN
            IF ( IKN.LT.10 ) WRITE (CSTR1(47:47),'(i1)') IKN
            CSTR2 = CSTR2(1:IP+ILS)//CSTR1
            IP = IP + ILS
         END DO
         WRITE (CSTR2(58:58),'(" ")')
         WRITE (IOTMP,'(a4000)') CSTR2
         WRITE (IOTMP,99003)
C     ---------------------------------------- plotting the red k points
         CSTR1 = ',''tiles_red.dat'' every :::  ::   w p ps 3 t "  "'
         DO IKN = 1,NKN
            IF ( IKN.GT.09 ) WRITE (CSTR1(27:28),'(i2)') IKN - 1
            IF ( IKN.LT.10 ) WRITE (CSTR1(27:27),'(i1)') IKN - 1
            IF ( IKN.GT.09 ) WRITE (CSTR1(31:32),'(i2)') IKN - 1
            IF ( IKN.LT.10 ) WRITE (CSTR1(31:31),'(i1)') IKN - 1
            IF ( IKN.GT.09 ) WRITE (CSTR1(46:47),'(i2)') IKN
            IF ( IKN.LT.10 ) WRITE (CSTR1(47:47),'(i1)') IKN
            CSTR2 = CSTR2(1:IP+ILS)//CSTR1
            IP = IP + ILS
         END DO
C
         WRITE (IOTMP,'(a4000)') CSTR2
         WRITE (IOTMP,99003)
C
         CLOSE (IOTMP)
C
      END IF
C -------------------- write k-points to file to be read in main program
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      WRITE (IOTMP) (NKTABKN(IKN),IKN=1,NKN)
      WRITE (IOTMP) (((KTABKN(I,IK,IKN),I=1,3),WKTABKN(IK,IKN),IK=1,
     &              NKTABKN(IKN)),IKN=1,NKN)
C
      DEALLOCATE (KTABD,WTABD,SYM_EVER_ACC,KTABKN,WKTABKN)
C
99001 FORMAT ('set xrange[-2:2]; set yrange[-2:2]; set zrange[-2:2]')
99002 FORMAT ('0,0,0',' to ',F10.6,',',F10.6,',',F10.6)
99003 FORMAT ('pause -1')
99004 FORMAT (5x,'nktabmaxkn, divisions',8x,':',2I10)
      END
C*==kvecteq.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      LOGICAL FUNCTION KVECTEQ(K1,K2,ABAS)
C   ********************************************************************
C   *  DK, January 2006                                                *
C   *                                                                  *
C   *  Checks whether k1 can be translated into k2 by a reciprocal     *
C   *  lattice vector                                                  *
C   *                                                                  *
C   *                T                                                 *
C   * +  uses  :    A  *  B = I                                        *
C   *                          3                                       *
C   *                                                                  *
C   *          where A is matrix of column vectors of real space unit  *
C   *                     cell in units of alat                        *
C   *                B is matrix of column vectors of reciprocal space *
C   *                     unit cell in units  2Pi/alat                 *
C   *                                                                  *
C   * + checks :   whether                                             *
C   *                                                                  *
C   *                                                                  *
C   *   [ m       -1                T                                  *
C   *     n    = B   * (k1 - k2) = A * (k1 - k2)  is an integer vector *
C   *     p ]                                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      REAL*8 ABAS(3,3),K1(3),K2(3)
C
C Local variables
C
      REAL*8 AR(3,3),KD(3),MNP(3),MNPDIFF
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      AR(1:3,1:3) = ABAS(1:3,1:3)
C
      KVECTEQ = .FALSE.
C
      KD(1:3) = K2(1:3) - K1(1:3)
C
      CALL DGEMV('N',3,3,1D0,AR,3,KD,1,0D0,MNP,1)
C
      MNPDIFF = 0D0
      DO I = 1,3
         MNPDIFF = MNPDIFF + ABS(MNP(I)-NINT(MNP(I)))
      END DO
C
      IF ( MNPDIFF.LT.TOL ) KVECTEQ = .TRUE.
C
      END
C*==nlcpawr.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPAWR(NLCPAOUT,IOTMP,DATSET,LDATSET,SYSTEM,LSYSTEM,
     &                   IQCPA,ITOQ,CONC,EFERMI,IECURR,ERYD,NETAB,EVI,
     &                   EVD,PCFG,IOCCCFGQ,NQNLCPA,NCFG,NME,NTMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *  NLCPAOUT determines additional output created     PROGRAM       *
C   *    2      expectation values for each config.        GEN         *
C   *           written to file if IECURR = NETAB                      *
C   *    3      DOS for each configuration                 GEN         *
C   *           written to file if IECURR = NETAB                      *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:TXT_T,LTXT_T
      USE MOD_ANGMOM,ONLY:NSPIN,IDOS,ISMT
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NLCPAWR')
C
C Dummy arguments
C
      CHARACTER*80 DATSET,SYSTEM
      REAL*8 EFERMI
      COMPLEX*16 ERYD
      INTEGER IECURR,IOTMP,IQCPA,LDATSET,LSYSTEM,NCFG,NLCPAOUT,NME,
     &        NQMAX,NQNLCPA,NTMAX
      REAL*8 CONC(NTMAX),PCFG(NCFG)
      COMPLEX*16 EVD(NCFG,NQNLCPA,NME),EVI(NCFG,NQNLCPA,NME)
      INTEGER IOCCCFGQ(NCFG,NQNLCPA),ITOQ(NTMAX,NQMAX),NETAB(2)
C
C Local variables
C
      REAL*8 DOS(:,:,:,:),MDNS,NDNS,WAVG,X(:),XAVG(1),XCFG,XMAX,XMIN,
     &       Y(:),YAVG(1),YMAX,YMIN,YSCL(4)
      CHARACTER*80 DTXT1,DTXT2,FILNAM,HTXT(4),YTXT
      CHARACTER*7 EXT(4)
      INTEGER IA_ERR,ICALL,ICFG,IE,IME,IOCC,IOCC1,IOCC2,IQCLU,IS,IT,IT1,
     &        ITA,IX,IY,IY1,IY2,LDTXT1,LDTXT2,LFILNAM,LHTXT(4),LS,LYTXT,
     &        LYTXT1(4),LYTXT2(4),NE,NX,NXMAX,NYOCC(2),NYOCCMAX
      CHARACTER*20 STR20
      CHARACTER*30 YTXT1(4),YTXT2(4)
      SAVE DOS,NE,NXMAX,NYOCC,NYOCCMAX,X,Y
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
      DATA EXT/'CFG_CHR','CFG_SMT','CFG_OMT','CFG_HFF'/
      DATA YSCL/1D0,1D0,1D0,1D-3/
      DATA YTXT1/'q!m{1}!s','!xm!0!m{1}!s','!xm!0!m{1}!s','B!m{1}!s'/
      DATA YTXT2/'  !N (e)!M{1}!Sval','  !N (!xm!0!sB!N)!M{1}!Sspin',
     &     '  !N (!xm!0!sB!N)!M{1}!Sorb','  !N (kG)!M{1}!Shf'/
      DATA LYTXT1/8,12,12,8/,LYTXT2/19,28,27,18/
      DATA HTXT/'valence band charge for       ',
     &     'spin magnetic moment for      ',
     &     'orbtal magnetic moment for    ',
     &     'magnetic hyperfine field for  '/,LHTXT/24,25,27,29/
C
      ALLOCATABLE DOS,X,Y
C
C   ********************************************************************
C                   INITIALISATION - START
C   ********************************************************************
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
         IF ( NLCPAOUT.EQ.2 ) THEN
            NE = NETAB(1)
C
            DO IOCC = 1,2
               NYOCC(IOCC) = 0
            END DO
            DO ICFG = 1,NCFG
               DO IQCLU = 1,NQNLCPA
                  IOCC = IOCCCFGQ(ICFG,IQCLU)
                  NYOCC(IOCC) = NYOCC(IOCC) + 1
               END DO
            END DO
            NYOCCMAX = MAX(NYOCC(1),NYOCC(2))
C
            ALLOCATE (DOS(NE,NYOCCMAX,2,2),X(NETAB(1)),Y(1),STAT=IA_ERR)
         ELSE
            NE = 1
            ALLOCATE (DOS(1,1,1,1))
C
            NXMAX = NCFG*NQNLCPA
            ALLOCATE (X(NXMAX),Y(NXMAX),STAT=IA_ERR)
         END IF
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DOS')
      END IF
C   ********************************************************************
C
C=======================================================================
C                            write  DOS
C=======================================================================
      IF ( NLCPAOUT.EQ.2 ) THEN
C
         IY1 = 0
         IY2 = 0
         DO IQCLU = 1,NQNLCPA
            DO ICFG = 1,NCFG
               IOCC = IOCCCFGQ(ICFG,IQCLU)
               NDNS = DIMAG(EVD(ICFG,IQCLU,IDOS))
               MDNS = DIMAG(EVD(ICFG,IQCLU,ISMT))
               IF ( IOCC.EQ.1 ) THEN
                  IY1 = IY1 + 1
                  DOS(IECURR,IY1,1,IOCC) = (NDNS+MDNS)*0.5D0
                  DOS(IECURR,IY1,2,IOCC) = (NDNS-MDNS)*0.5D0
               ELSE
                  IY2 = IY2 + 1
                  DOS(IECURR,IY2,1,IOCC) = (NDNS+MDNS)*0.5D0
                  DOS(IECURR,IY2,2,IOCC) = (NDNS-MDNS)*0.5D0
               END IF
            END DO
         END DO
         X(IECURR) = (DREAL(ERYD)-EFERMI)*RY_EV
C
         IF ( IECURR.NE.NETAB(1) ) RETURN
C
C ----------------------------------------------------------------------
C                              XMGRACE - output
C ----------------------------------------------------------------------
         XMIN = X(1)
         XMAX = X(NETAB(1))
C
         DO IOCC = 1,2
            IT = ITOQ(IOCC,IQCPA)
C
C
            YMAX = 0D0
            DO IS = 1,NSPIN
               DO IY = 1,NYOCC(IOCC)
                  DO IE = 1,NETAB(1)
                     YMAX = MAX(YMAX,DOS(IE,IY,IS,IOCC))
                  END DO
               END DO
            END DO
C
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(E) (sts./eV)'
            LS = LTXT_T(IT) + 15
            IF ( NSPIN.EQ.1 ) THEN
               DTXT1 = 'n!s'//STR20(1:LS)
               LDTXT1 = 3 + LS
               DTXT2 = ' '
               LDTXT2 = 1
            ELSE
               DTXT1 = 'n!m{1}!S!UP!M{1}!N!s'//STR20(1:LS)
               LDTXT1 = 20 + LS
               DTXT2 = 'n!m{1}!S!DN!M{1}!N!s'//STR20(1:LS)
               LDTXT2 = 20 + LS
            END IF
C
            CALL XMGRHEAD(DATSET,LDATSET,'CFG_DOS',7,TXT_T(IT),
     &                    LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,NSPIN,XMIN,
     &                    1,XMAX,1,0.0D0,0,YMAX,1,0.0D0,0,YMAX,1,
     &                    'energy (eV)',11,DTXT1,LDTXT1,DTXT2,LDTXT2,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'NLCPA-DOS of '//TXT_T(IT)
     &                    (1:LTXT_T(IT))//' in '//SYSTEM(1:LSYSTEM),
     &                    (13+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
C         CALL XMGRCURVES(IOTMP,NSPIN,NCURVES,NCURVES,2,1,0)
C
            DO IS = 1,NSPIN
               DO IY = 1,NYOCC(IOCC)
                  CALL XMGRTABLE((IS-1),IY-1,X,DOS(1,IY,IS,IOCC),1D0,NE,
     &                           IOTMP)
               END DO
            END DO
C
            WRITE (6,*) '  '
            WRITE (6,*) '   DOS written to file  ',FILNAM(1:LFILNAM)
            CLOSE (IOTMP)
C
         END DO
C=======================================================================
C                      write  expectation values
C=======================================================================
      ELSE
C
         IF ( IECURR.NE.NETAB(1) ) RETURN
C
         ITA = ITOQ(1,IQCPA)
C
         DO IME = 1,NME
            DO IOCC1 = 1,2
               IT1 = ITOQ(IOCC1,IQCPA)
               ITA = IT1
C
               YMIN = 0D0
               YMAX = 0D0
               YAVG(1) = 0D0
               WAVG = 0D0
               IX = 0
               DO ICFG = 1,NCFG
                  XCFG = 0D0
                  DO IQCLU = 1,NQNLCPA
                     IOCC2 = IOCCCFGQ(ICFG,IQCLU)
                     IF ( IOCC1.EQ.IOCC2 ) XCFG = XCFG + 1D0
                  END DO
C                 IF( IT1 .NE. ITA ) XCFG = NQNLCPA + 1 - XCFG
                  DO IQCLU = 1,NQNLCPA
                     IOCC2 = IOCCCFGQ(ICFG,IQCLU)
                     IF ( IOCC1.EQ.IOCC2 ) THEN
                        IX = IX + 1
                        Y(IX) = DIMAG(EVI(ICFG,IQCLU,IME))
                        X(IX) = XCFG
                        YAVG(1) = YAVG(1) + PCFG(ICFG)*Y(IX)
                        WAVG = WAVG + PCFG(ICFG)
                        YMIN = MIN(YMIN,Y(IX))
                        YMAX = MAX(YMAX,Y(IX))
                     END IF
                  END DO
               END DO
C
               YAVG(1) = YAVG(1)/WAVG
               XAVG(1) = CONC(IT1)*NQNLCPA
C              IF( IT1 .NE. ITA ) XAVG(1) = NQNLCPA + 1 - XAVG(1)
C
               YMIN = YMIN*YSCL(IME)
               YMAX = YMAX*YSCL(IME)
               NX = IX
               IF ( NX.GT.NXMAX )
     &              CALL STOP_MESSAGE(ROUTINE,'NX > NXMAX')
               XMIN = 0.0001D0
               XMAX = DBLE(NQNLCPA) + 0.9999D0
C
               YTXT = YTXT1(IME)(1:LYTXT1(IME))//TXT_T(IT1)
     &                (1:LTXT_T(IT1))//YTXT2(IME)(1:LYTXT2(IME))
               LYTXT = LYTXT1(IME) + LTXT_T(IT1) + LYTXT2(IME)
C
               CALL XMGRHEAD(DATSET,LDATSET,EXT(IME),7,TXT_T(IT1),
     &                       LTXT_T(IT1),FILNAM,80,LFILNAM,IOTMP,1,XMIN,
     &                       0,XMAX,0,YMIN,1,YMAX,1,YMIN,0,YMAX,0,
     &                       'cluster occupation by '//TXT_T(ITA)
     &                       (1:LTXT_T(IT1)),(22+LTXT_T(IT1)),YTXT,
     &                       LYTXT,' ',0,
     &                       'SPR-KKR NLCPA-calculations for '//
     &                       SYSTEM(1:LSYSTEM),31+LSYSTEM,HTXT(IME)
     &                       (1:LHTXT(IME))//TXT_T(IT1)(1:LTXT_T(IT1)),
     &                       LHTXT(IME)+LTXT_T(IT1),.FALSE.)
C
C              CALL XMGRLEGEND(IOTMP,1,NT,0,LEG,LEG)
C
               CALL XMGRPOINTS(IOTMP,1,2,0,0,1,1)
C
               CALL XMGRTICK(IOTMP,0,'x','major',1D0)
               CALL XMGRTICK(IOTMP,0,'x','minor',1D0)
C
               CALL XMGRTABLE(0,0,X,Y,YSCL(IME),NX,IOTMP)
C
               CALL XMGRTABLE(0,1,XAVG,YAVG,YSCL(IME),1,IOTMP)
C
               WRITE (6,*) ' '
               WRITE (6,*) '    results written to xmgrace file ',
     &                     FILNAM(1:LFILNAM)
               WRITE (6,*) ' '
               CLOSE (IOTMP)
C
            END DO
         END DO
C
C=======================================================================
      END IF
      END
C*==rotatetaukn.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ROTATETAUKN(TAUKND,TAUKN,W1,DROT,NKM,NKMMAX,SYMUNITARY,
     &                       ISYM)
C   ********************************************************************
C   *                                                                  *
C   *   Rotate tau-matrix and symmetrise                               *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ISYM,NKM,NKMMAX
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),TAUKN(NKMMAX,NKMMAX),
     &           TAUKND(NKMMAX,NKMMAX),W1(NKMMAX,NKMMAX)
      LOGICAL SYMUNITARY(NSYMMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      INTEGER N
C
C*** End of declarations rewritten by SPAG
C
      IF ( SYMUNITARY(ISYM) ) THEN
         CNT = 'N'
      ELSE
         CNT = 'T'
      END IF
C
      N = NKM
C
      CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),NKMMAX,TAUKND(1,1),
     &           NKMMAX,C0,W1,NKMMAX)
C
      CALL ZGEMM('N','C',N,N,N,C1,W1,NKMMAX,DROT(1,1,ISYM),NKMMAX,C0,
     &           TAUKN(1,1),NKMMAX)
C
      END
C*==nlcpaproj.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE NLCPAPROJ(ICFG,IQCPA,MSST,OMEGAHAT,TAUIMP,W1HAT,
     &                     NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,NTMAX,
     &                     NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  perform the   NLCPA  projection procedure                       *
C   *                                                                  *
C   *  N: dimension of the block matrices w.r.t. (kappa,mue)           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ICFG,IQCPA,N,NDIMCLU,NKMMAX,NQMAX,NQNLCPA,NTMAX
      INTEGER IND0QCLU(NQNLCPA),ITOQ(NTMAX,NQMAX)
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),OMEGAHAT(NDIMCLU,NDIMCLU),
     &           TAUIMP(NDIMCLU,NDIMCLU),W1HAT(NDIMCLU,NDIMCLU)
C
C Local variables
C
      INTEGER I0,IOCC,IQCLU,IT,M
      INTEGER NLCPACONF
C
C*** End of declarations rewritten by SPAG
C
      M = NDIMCLU
C
      W1HAT(1:M,1:M) = -OMEGAHAT(1:M,1:M)
C
      DO IQCLU = 1,NQNLCPA
C
         IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
C
         IT = ITOQ(IOCC,IQCPA)
         I0 = IND0QCLU(IQCLU)
C
         W1HAT(I0+1:I0+N,I0+1:I0+N) = MSST(1:N,1:N,IT)
     &                                + W1HAT(I0+1:I0+N,I0+1:I0+N)
      END DO
C
      CALL CMATINV(M,M,W1HAT,TAUIMP)
C
      END
