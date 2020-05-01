C*==vbpes.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE VBPES(TASK,TSST,MSST,NE,MEZZ,MEZJ,IPRINT,PHASK,TAUQ,
     &                 TAUT,CALCINT,GETIRRSOL,NPOLMAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of the      VB - PES       spectra             *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *   VBPES   angle integrated valence band photo emission spectrum  *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   ARPES   angle resolved valence band photo emission spectrum    *
C   *                                                                  *
C   *   NK = 1     scan E and fix k_fin = (E_fin - E_vac)^1/2          *
C   *  (default)                                                       *
C   *                                                                  *
C   *   NK > 1     scan E and k_ini                                    *
C   *              specify path for k_ini via                          *
C   *              - KPATH = KEY    see settings in <KDIRTAB>          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NETAB,NEMAX,ETAB,EFERMI,EMAX,EMIN,IGRID
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_CPA,ONLY:NCPA
      USE MOD_CALCMODE,ONLY:DMFT,ORBPOL
      USE MOD_TYPES,ONLY:NT,NCPLWFMAX,NTMAX,SOCTL,CTL,NLT
      USE MOD_SITES,ONLY:NQ,NQMAX,IQAT
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NLMAX,NLINQ,IND0Q,NKKR,NKMQ,NLQ,
     &    NKM,NLM,NL,WKM1,WKM2
      USE MOD_FILES,ONLY:SYSTEM,DATSET,LDATSET,LRECREAL8,LSYSTEM,WRTAU,
     &    IFILTAU,IFILDOS,FOUND_REAL,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_CONSTANTS,ONLY:C1,C0,CI,A0_ANG,PI,RY_EV,C_AU
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='VBPES')
      REAL*8 R05
      PARAMETER (R05=1D0/1.4142135623730951D0)
      COMPLEX*16 CR05,IR05,MIR05
      PARAMETER (CR05=C1*R05,IR05=CI*R05,MIR05=-IR05)
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER IPRINT,NE,NPOLMAX
      CHARACTER*10 TASK
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),PHASK(NEMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CCYY(:,:,:,:),ERYD,ETABFIN(:,:),ETABINI(:,:),
     &           MBAR(:,:,:),MCS(3,3),MEBAR(:,:,:,:),MEIRR(:,:,:,:,:),
     &           MEREG(:,:,:,:),MG0(3,3),MG1(3,3),MG2(3,3),MGEO(3,3),
     &           MIRR(:,:,:,:),MREG(:,:,:),MSSQ(:,:,:),MSSQFIN(:,:,:),
     &           MSSTFIN(:,:,:),SA(:,:,:),SL(:,:,:,:),SSST(:,:,:),
     &           ST(:,:,:,:),TIE(:,:),TM(:,:,:),TMBAR(:,:,:),TSSQ(:,:,:)
     &           ,TSSQFIN(:,:,:),TSSTFIN(:,:,:)
      CHARACTER*1 CHPOL(5)
      REAL*8 DEFEP0,DEFEP0EV,EPHOT,EPHOTEV,EVAC,EWORK,EWORKEV,IMEFIN,
     &       KUNIT,MJ,MJLAM(NKMMAX),MSPN(2),NSPIN(3),PHIPHOT,QHAT(3),
     &       QPHOT,QVECPHOT(3),RJ,RNORM,SK,TETPHOT,TIME,TIME0
      REAL*8 DNRM2
      CHARACTER*80 FILNAM,SPEC,TAUFIL,TAUFNNFIL
      INTEGER I,IAME,IA_ERR,IE,IFILFIN,IFILFINA,IFILINI,IFILINIA,
     &        IFILTAUFIN,IFILTAUINI,IFILTAUNN,IOL,IPRINTTAU,IQ,IT,IWME,
     &        K,KAP,KGEO,L,LAM,LFN,LLAM(NKMMAX),LSEPHOTEV,LSP,
     &        LTAUFNNFIL,M,MJM05,MODE,NKMFIN,NKMINI,NPOL,NQ9,NT9
      CHARACTER*3 MEFORM
      CHARACTER*10 ORBPOLFIN,SEPHOTEV,STR10
      CHARACTER*4 STR4
      LOGICAL TAUNNAVAILABLE,USETAUNN,USETFIN
C
C*** End of declarations rewritten by SPAG
C
      DATA IAME/0/
      DATA CHPOL/'+','-','z','x','y'/
      DATA MCS/CR05,MIR05,C0,CR05,IR05,C0,C0,C0,C1/
      DATA MG0/C1,C0,C0,C0,C1,C0,C0,C0,C1/
      DATA MG1/C0,C0,C1,CR05,CR05,C0,IR05,MIR05,C0/
      DATA MG2/IR05,MIR05,C0,C0,C0,C1,CR05,CR05,C0/
      DATA EWORKEV/0D0/
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE TMBAR
      ALLOCATABLE MBAR,MREG,MIRR,MEIRR
      ALLOCATABLE CCYY,ETABINI,ETABFIN,MEREG,MSSQ,SA,SL,SSST,ST
      ALLOCATABLE TIE,TM,MEBAR,TSSTFIN,MSSTFIN,TSSQ,TSSQFIN,MSSQFIN
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C ----------------------------- increase angular momentum expansion by 1
C --------------------- to deal with final states with l_fin = l_ini + 1
      NL = 0
C
C*** Start of declarations rewritten by SPAG
C
      NKKR = 0
      DO IQ = 1,NQ
         NLQ(IQ) = NLQ(IQ) + 1
         NL = MAX(NL,NLQ(IQ))
         IF ( NLQ(IQ).GT.NLMAX )
     &        CALL STOP_MESSAGE(ROUTINE,'NLQ > NLMAX')
         NKMQ(IQ) = 2*NLQ(IQ)**2
         NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
         NKKR = NKKR + 2*NLQ(IQ)**2
         IF ( IQ.EQ.1 ) THEN
            IND0Q(IQ) = 0
         ELSE
            IND0Q(IQ) = IND0Q(IQ-1) + 2*NLQ(IQ-1)**2
         END IF
      END DO
C
      DO IT = 1,NT
         NLT(IT) = NLQ(IQAT(1,IT))
         CTL(IT,NLT(IT)) = CTL(IT,NLT(IT)-1)
         SOCTL(IT,NLT(IT)) = SOCTL(IT,NLT(IT)-1)
      END DO
C
      LAM = 0
      MSPN(1) = -0.5D0
      MSPN(2) = 0.5D0
      DO K = 1,2*NLMAX - 1
         IF ( MOD(K,2).EQ.0 ) THEN
            L = K/2
            KAP = +L
         ELSE
            L = (K-1)/2
            KAP = -L - 1
         END IF
         SK = DBLE(SIGN(1,KAP))
         RJ = DBLE(L) - SK/2D0
         DO MJM05 = NINT(-RJ-0.5D0),NINT(RJ-0.5D0)
            MJ = DBLE(MJM05) + 0.5D0
            LAM = LAM + 1
            LLAM(LAM) = L
            MJLAM(LAM) = MJ
         END DO
      END DO
      NLM = NL**2
      NKM = 2*NLM
      NKMINI = 2*(NL-1)**2
      NKMFIN = NKM
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      M = NKMMAX
C
      ALLOCATE (MSSQ(M,M,NQMAX),MSSQFIN(M,M,NQMAX))
      ALLOCATE (TSSQ(M,M,NQMAX),TSSQFIN(M,M,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MSSQ')
      ALLOCATE (TSSTFIN(M,M,NTMAX),MSSTFIN(M,M,NTMAX))
      ALLOCATE (SSST(M,M,NTMAX),TM(M,M,NT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: SSST')
      ALLOCATE (TMBAR(M,M,NT),MBAR(M,M,3),MREG(M,M,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TMBAR')
      ALLOCATE (MEBAR(M,M,3,NT),MEREG(M,M,3,NT))
      ALLOCATE (SL(NLMAX,2,2,NPOLMAX),SA(2,2,NPOLMAX))
      ALLOCATE (ST(2,2,NPOLMAX,NTMAX),TIE(M,M))
      ALLOCATE (CCYY(M,2,M,2),ETABINI(NEMAX,2),ETABFIN(NEMAX,2))
      ALLOCATE (MEIRR(M,M,M,3,NT),MIRR(M,M,3,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEIRR')
C
      CALL CINIT(1,SL)
C
      IF ( TASK(3:5).NE.'PES' ) CALL STOP_MESSAGE(ROUTINE,'TASK ???')
C=======================================================================
C
      CALL CINIT(NKM*NKM*2*2,CCYY)
C
      CALL CPU_TIME(TIME0)
C
C=======================================================================
C                     Reading input file
C=======================================================================
C
      NE = NETAB(1)
      KUNIT = 2D0*PI/ALAT
      IPRINTTAU = 0
      GETIRRSOL = .TRUE.
      IWME = 0
C
      IF ( DMFT ) THEN
         ORBPOLFIN = 'NONE      '
      ELSE
         ORBPOLFIN = ORBPOL
      END IF
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      CALL SECTION_SET_REAL('EPHOT',EPHOTEV,1253.6D0,0)
C
      IF ( NINT(EPHOTEV).EQ.1 ) EPHOTEV = 1253.6D0
      IF ( NINT(EPHOTEV).EQ.2 ) EPHOTEV = 1486.6D0
C
      EPHOT = EPHOTEV/RY_EV
C
      CALL SECTION_SET_REAL('DEFEP0',DEFEP0EV,0.0D0,0)
C
      DEFEP0 = DEFEP0EV/RY_EV
C
      CALL SECTION_SET_REAL_ARRAY('QVEC',QHAT,N_FOUND,3,0,9999D0,0)
C
      IF ( FOUND_REAL_ARRAY ) THEN
         RNORM = 1D0/DNRM2(3,QHAT,1)
         CALL DSCAL(3,RNORM,QHAT,1)
      ELSE
         CALL RINIT(3,QHAT)
      END IF
      QPHOT = EPHOT/C_AU
C
      WRITE (SEPHOTEV,FMT='(I10)') INT(EPHOTEV)
      CALL STRING_TRIM_LEFT(SEPHOTEV)
      LSEPHOTEV = LEN_TRIM(SEPHOTEV)
C
      CALL SECTION_SET_STRING('ME',MEFORM,'ADA',0)
      CALL SECTION_SET_INTEGER('MODE',MODE,2,0)
C
      IF ( FULLPOT ) MODE = 2
C
      CALL SECTION_SET_REAL('EWORK',EWORKEV,9999D0,0)
C
      EWORK = EWORKEV/RY_EV
      EVAC = EFERMI + EWORK
C
      KGEO = 0
      CALL SECTION_SET_INTEGER('geometry',KGEO,9999,0)
C
      IF ( KGEO.EQ.0 ) THEN
         CALL ZCOPY(3*3,MG0,1,MGEO,1)
         QHAT(3) = 1D0
      ELSE IF ( KGEO.EQ.1 ) THEN
         CALL ZGEMM('N','N',3,3,3,C1,MG1,3,MCS,3,C0,MGEO,3)
         IF ( IPRINT.GT.0 ) CALL CMATSTRUCT('MG1  ',MG1,3,3,0,0,0,1D-8,
     &        6)
         IF ( IPRINT.GT.0 ) CALL CMATSTRUCT('MGEO ',MGEO,3,3,0,0,0,1D-8,
     &        6)
         QHAT(1) = 1D0
      ELSE IF ( KGEO.EQ.2 ) THEN
         CALL ZGEMM('N','N',3,3,3,C1,MG2,3,MCS,3,C0,MGEO,3)
         IF ( IPRINT.GT.0 ) CALL CMATSTRUCT('MG2  ',MG2,3,3,0,0,0,1D-8,
     &        6)
         IF ( IPRINT.GT.0 ) CALL CMATSTRUCT('MGEO ',MGEO,3,3,0,0,0,1D-8,
     &        6)
         QHAT(2) = 1D0
      ELSE
         CALL SECTION_SET_REAL('TETPH',TETPHOT,9999D0,0)
         IF ( FOUND_REAL ) THEN
            CALL SECTION_SET_REAL('PHIPH',PHIPHOT,9999D0,1)
            KGEO = 3
         ELSE
            CALL SECTION_SET_REAL_ARRAY('QVECPH',QVECPHOT,N_FOUND,3,0,
     &                                  9999D0,0)
            IF ( .NOT.FOUND_REAL_ARRAY )
     &            CALL STOP_MESSAGE(ROUTINE,'QPHOT not found')
            KGEO = 4
         END IF
C
         CALL POLCONV(TETPHOT,PHIPHOT,QVECPHOT,QHAT,KGEO,MGEO)
C
      END IF
      IF ( IPRINT.GE.0 ) CALL CMATSTRUCT('MGEO ',MGEO,3,3,0,0,0,1D-8,6)
C
      CALL SECTION_SET_REAL_ARRAY('NSPIN',NSPIN,N_FOUND,3,0,9999D0,0)
C
      IF ( FOUND_REAL_ARRAY ) THEN
         RNORM = 1D0/DNRM2(3,NSPIN,1)
         CALL DSCAL(3,RNORM,NSPIN,1)
      ELSE
         CALL DCOPY(3,QHAT,1,NSPIN,1)
         CALL DSCAL(3,-1D0,NSPIN,1)
      END IF
C
C=======================================================================
C                        TAU for initial states
C=======================================================================
C  - scan file 9 for TAU_ini and set  E-mesh
C  - reposition file 9
C
      IF ( TASK(1:5).NE.'ARPES' ) THEN
C
         WRITE (6,*) ' reading TAU-file for initial states'
C
         CALL READTAU(9,ERYD,0,NE,TAUT,0,NT,.TRUE.,WKM1,WKM2,0,NQ,0,
     &                NKMMAX,NKMMAX,IPRINTTAU)
C
         IF ( NE.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'NE=0 in initial state TAU-file')
C
         DO IE = 1,NE
            CALL READTAU(9,ERYD,IE,NE,TAUT(1,1,1),1,NT,.TRUE.,WKM1,WKM2,
     &                   0,NQ,1,NKMINI,NKMMAX,IPRINTTAU)
            ETAB(IE,1) = ERYD
            ETAB(IE,2) = ERYD
         END DO
         NETAB(1) = NE
         IGRID(1) = 3
C
         REWIND 9
         READ (9,'(////,10X,I3,5X,I3)') NT9,NQ9
         IF ( NQ.NE.NQ9 ) CALL STOP_MESSAGE(ROUTINE,'NQ <> NQ9')
         READ (9,*) (STR4,I=1,(NT9+NQ9))
         IF ( (NT9+NQ9).LT.1 ) WRITE (6,*) '<VBPES>: ',STR4
C
      ELSE IF ( NCPA.NE.0 ) THEN
C
         WRITE (6,*) ' reading TAU-file for initial states'
C
         CALL READTAU(9,ERYD,0,NE,WKM1,0,NT,.FALSE.,WKM1,WKM2,0,NQ,0,
     &                NKMMAX,NKMMAX,IPRINT)
C
         IF ( NE.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'NE=0 in initial state TAU-file')
C
         DO IE = 1,NE
            ETAB(IE,2) = ETAB(IE,1)
         END DO
         NETAB(1) = NE
         IGRID(1) = 3
C
      END IF
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
C=======================================================================
C                       TAU for final states
C=======================================================================
C  - set up E-mesh for initial and final states
C  - if switch  USETAUNN  is set
C      - check whether file with  TAU_fin  is available
C      - if not: calculate  TAU_fin  by calling   <ELOOP>
C
      CALL SECTION_SET_REAL('IMEFIN',IMEFIN,DIMAG(ETAB(NE,1)),0)
C
      DO IE = 1,NE
         ETABINI(IE,1) = ETAB(IE,1)
         ETABFIN(IE,1) = DCMPLX(DREAL(ETAB(IE,1)),IMEFIN) + EPHOT
      END DO
C
      CALL SECTION_FIND_KEYWORD('USETAUNN',USETAUNN)
C
      USETFIN = .FALSE.
      CALL SECTION_FIND_KEYWORD('USETFIN',USETFIN)
      IF ( USETFIN ) USETAUNN = .FALSE.
C
      IF ( USETAUNN ) THEN
C
         WRITE (6,99002)
C
         DO IE = 1,NE
            ETABFIN(IE,1) = DCMPLX(DREAL(ETAB(IE,1)),IMEFIN) + EPHOT
         END DO
C
         TAUFNNFIL = DATSET(1:LDATSET)//'.tau.'//'nnF_'//
     &               SEPHOTEV(1:LSEPHOTEV)
         FILNAM = DATSET(1:LDATSET)//'.dos.'//'nnF_'//
     &            SEPHOTEV(1:LSEPHOTEV)
         LTAUFNNFIL = LDATSET + 5 + LSEPHOTEV + 4
C
         IFILTAUNN = 20
         INQUIRE (FILE=TAUFNNFIL(1:LTAUFNNFIL),EXIST=TAUNNAVAILABLE)
C
         OPEN (IFILTAUNN,FILE=TAUFNNFIL(1:LTAUFNNFIL))
C
         IF ( TAUNNAVAILABLE ) THEN
C
            WRITE (6,99003) TAUFNNFIL(1:LTAUFNNFIL)
C
            READ (IFILTAUNN,'(////,10X,I3,5X,I3)') NT9,NQ9
            IF ( NQ.NE.NQ9 ) CALL STOP_MESSAGE(ROUTINE,'NQ <> NQ9')
            READ (IFILTAUNN,*) (STR4,I=1,(NT9+NQ9))
            IF ( (NT9+NQ9).LT.1 ) WRITE (6,*) '<VBPES>: ',STR4
C
         ELSE
C
C----------------------- use IFILTAU temporarily for final state TAU-file
            CLOSE (IFILTAUNN)
            CLOSE (IFILTAU)
C
            OPEN (IFILTAU,FILE=TAUFNNFIL(1:LTAUFNNFIL))
C
            WRITE (6,99004)
            WRITE (6,99006) 'TAU - file',IFILTAUNN,
     &                      TAUFNNFIL(1:LTAUFNNFIL)
C
            WRITE (6,'(10X,A,A)') 'DOS - file:  (10) ',
     &                            FILNAM(1:LTAUFNNFIL)
            OPEN (UNIT=10,FILE=FILNAM(1:LTAUFNNFIL))
C
            CALL WRHEAD(IFILDOS,FILNAM,'DOS       ',NETAB(1))
            WRITE (IFILDOS,99011) 'DOS-FMT:  ','OLD-SPRKKR'
C
C
C-------------------------- use ETAB temporarily for final states energy
            ETAB(1:NE,1) = ETABFIN(1:NE,1)
C
            WRTAU = .TRUE.
C
            CALL ELOOP(IPRINT,PHASK,TAUQ,TAUT,CALCINT,GETIRRSOL)
C
            ETAB(1:NE,1) = ETABINI(1:NE,1)
C
            TAUFIL = DATSET(1:LDATSET)//'.tau'
C
C
C--------------------------------- recover IFILTAU for initial TAU-file
            CLOSE (IFILTAU)
            OPEN (UNIT=IFILTAU,FILE=TAUFIL)
            REWIND IFILTAU
C
            CALL READTAU(IFILTAU,ERYD,0,NE,TAUT,0,NT,(NCPA.EQ.0),WKM1,
     &                   WKM2,0,NQ,0,NKMMAX,NKMMAX,IPRINTTAU)
C
            OPEN (IFILTAUNN,FILE=TAUFNNFIL(1:LTAUFNNFIL))
C
         END IF
C
         REWIND IFILTAUNN
         READ (IFILTAUNN,'(////,10X,I3,5X,I3)') NT9,NQ9
         READ (IFILTAUNN,*) (STR4,I=1,(NT9+NQ))
         IF ( (NT9+NQ9).LT.1 ) WRITE (6,*) '<VBPES>: ',STR4
         IF ( NQ.NE.NQ9 ) CALL STOP_MESSAGE(ROUTINE,'NQ <> NQ9')
C
      END IF
C
C
C=======================================================================
C                  Initialisation of files
C=======================================================================
C
      OPEN (UNIT=99,FILE='arpes.dat')
      OPEN (UNIT=98,FILE='arpes_spol.dat')
      OPEN (UNIT=40,FILE='spin.matrix.xps')
      CLOSE (87)
      CLOSE (88)
      CLOSE (85)
      CLOSE (86)
      CLOSE (82)
      CLOSE (81)
      IFILINI = 87
      IFILFIN = 88
      IFILINIA = 85
      IFILFINA = 86
      IFILTAUINI = 82
      IFILTAUFIN = 81
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
         OPEN (IFILINI,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
         OPEN (IFILFIN,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
      ELSE
         OPEN (UNIT=IFILINI,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=(4+2*(2+4*NRMAX))*2*LRECREAL8)
         OPEN (UNIT=IFILFIN,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=(4+2*(2+4*NRMAX))*2*LRECREAL8)
         OPEN (UNIT=IFILINIA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=(4+2*(2+4*NRMAX))*2*LRECREAL8)
         OPEN (UNIT=IFILFINA,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=(4+2*(2+4*NRMAX))*2*LRECREAL8)
      END IF
      OPEN (UNIT=IFILTAUFIN,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(NKMMAX*NKMMAX*NTMAX*2*LRECREAL8))
      OPEN (UNIT=IFILTAUINI,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(NKMMAX*NKMMAX*NQMAX*NTMAX*2*LRECREAL8)
     &      )
C
      IF ( LDATSET.NE.0 ) THEN
         FILNAM = DATSET(1:LDATSET)//'.xps'
         LFN = LDATSET + 4
      ELSE
         FILNAM = 'vb-xps'
         LFN = 6
      END IF
      IF ( TASK(1:5).EQ.'ARPES' ) THEN
         STR10(1:10) = 'VB-AR-PES '
         SPEC = ' VB-AR-PES of '//SYSTEM(1:MIN(69,LSYSTEM))
         LSP = 14 + MIN(69,LSYSTEM)
      ELSE
         STR10(1:10) = 'VB-PES    '
         SPEC = ' VB-PES of '//SYSTEM(1:MIN(69,LSYSTEM))
         LSP = 11 + MIN(69,LSYSTEM)
      END IF
C
      OPEN (UNIT=7,FILE=FILNAM(1:LFN))
      WRITE (6,'(10X,A,A,/)') 'SPEC-FILE :  ( 7) ',FILNAM(1:LFN)
C
      CALL WRHEAD(7,FILNAM,STR10,NE)
C
      NPOL = 3
C
      WRITE (6,99005)
C
      WRITE (6,99016) 'SPECTRUM  ',SPEC(1:LSP)
      WRITE (6,99016) 'MEFORM    ',MEFORM
      WRITE (6,99012) 'MODE      ',MODE
      WRITE (6,99007) 'USETAUNN  ',USETAUNN
      WRITE (6,99007) 'USETFIN   ',USETFIN
      WRITE (6,99016) ' '
      WRITE (6,99014) 'E_ph      ',EPHOTEV,' eV    ',EPHOT,' Ry'
      WRITE (6,99014) 'q_ph      ',QPHOT,' a.u.  ',QPHOT/KUNIT,' 2*pi/a'
      WRITE (6,99015) '->q       ',QHAT
      WRITE (6,99013) 'NPOL      ',NPOL,(CHPOL(I),I=1,NPOL)
      WRITE (6,99016) ' '
      WRITE (6,99012) 'NE        ',NE
      WRITE (6,99014) 'EMIN      ',EMIN*RY_EV,' eV    ',EMIN,' Ry'
      WRITE (6,99014) 'EMAX      ',EMAX*RY_EV,' eV    ',EMAX,' Ry'
      WRITE (6,99014) 'E_Fermi   ',EFERMI*RY_EV,' eV    ',EFERMI,' Ry'
      WRITE (6,99014) 'E_work    ',EWORK*RY_EV,' eV    ',EWORK,' Ry'
      WRITE (6,99014) 'E_vac     ',EVAC*RY_EV,' eV    ',EVAC,' Ry'
      WRITE (6,99014) 'DEL_E0    ',DEFEP0EV,' eV    ',DEFEP0,' Ry'
      WRITE (6,99014) 'Im(E)_fin ',IMEFIN*RY_EV,' eV    ',IMEFIN,' Ry'
      WRITE (6,99016) ' '
      WRITE (6,99014) 'kunit     ',KUNIT/A0_ANG,' 1/Ang ',KUNIT,' 1/a0'
      WRITE (6,99012) 'NLM       ',NLM
      WRITE (6,99012) 'NKMINI    ',NKMINI
      WRITE (6,99012) 'NKMFIN    ',NKMFIN
      WRITE (6,99012) 'GEOMETRY  ',KGEO
      WRITE (6,99015) '->n(spin) ',NSPIN
C
      WRITE (7,99011) 'SPEC-FILE ',FILNAM(1:LFN)
      WRITE (7,99011) 'SPECTRUM  ',SPEC(1:LSP)
      WRITE (7,99011) 'MEFORM    ',MEFORM
      WRITE (7,99008) 'MODE      ',MODE
      WRITE (7,99010) 'E_ph (eV) ',EPHOTEV
      WRITE (7,99009) 'NPOL      ',NPOL,(CHPOL(I),I=1,NPOL)
      WRITE (7,99008) 'NE        ',NE
C
C ======================================================================
C                calculate angular matrix elements
C ======================================================================
C
      IF ( IPRINT.GT.0 ) IWME = 1
C
      CALL AME_INIT(MEFORM,IWME)
C
C-----------------------------------------------------------------------
C                               omit  AME for KINI(IK)=0 .OR. KFIN(IK)=0
C
      IF ( IAME.GT.0 ) CALL STOP_MESSAGE(ROUTINE,'IAME.GT.0')
C
C-----------------------------------------------------------------------
C
      IF ( TASK(1:5).EQ.'ARPES' ) THEN
C
         CALL STOP_MESSAGE(ROUTINE,
     &              'ARPES no more supplied within KKRGEN - use KKRSPEC'
     &              )
C
      ELSE
C
         CALL VBPESINT(MGEO,NE,MEZZ,MEZJ,NKMINI,NKMFIN,LLAM,MJLAM,MSPN,
     &                 IFILINI,IFILFIN,IFILTAUNN,MEFORM,IPRINT,CALCINT,
     &                 GETIRRSOL,MODE,NPOL,ORBPOLFIN,USETAUNN,KGEO,
     &                 NSPIN,TMBAR,MBAR,MREG,MIRR,MEIRR,CCYY,ETABFIN,
     &                 MEREG,SA,ST,TIE,TM,MEBAR,NPOLMAX,TSST,MSST,SSST,
     &                 TSSQ,MSSQ,TAUT,TSSTFIN,MSSTFIN,TSSQFIN,MSSQFIN)
C
      END IF
C
      CALL CPU_TIME(TIME)
      WRITE (6,99001) TIME - TIME0
C
      DEALLOCATE (TMBAR,MBAR,MREG,MIRR,MEIRR)
      DEALLOCATE (CCYY,ETABFIN,MEREG,MSSQFIN,MSSQ,SA,SL,SSST,ST)
      DEALLOCATE (TIE,TM,MEBAR,TSSTFIN)
C
      STOP
C
C-----------------------------------------------------------------------
C-----------------------------------                    FORMAT statments
C
99001 FORMAT (/,10X,'execution time for <VBPES>: ',F14.3,' secs',/)
99002 FORMAT (/,5X,'USETAUNN set! ',/,5X,
     &        'single scatterer approximation will not be used')
99003 FORMAT (/,10X,'File  ',A,'   for TAUnn(final) available !',/)
99004 FORMAT (/,10X,'TAUnn(fin) NOT available !',
     &        ' >>>  will be calculated ')
99005 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      *     *  *****          *****   ******    ****        *'
     &  ,/,10X,
     &  '*      *     *  *    *         *    *  *        *    *       *'
     &  ,/,10X,
     &  '*       *   *   *    *         *    *  *        *            *'
     &  ,/,10X,
     &  '*       *   *   *****    ***   *****   *****     ****        *'
     &  ,/,10X,
     &  '*        * *    *    *         *       *             *       *'
     &  ,/,10X,
     &  '*        * *    *    *         *       *        *    *       *'
     &  ,/,10X,
     &  '*         *     *****          *       ******    ****        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99006 FORMAT (10X,A,':  (',I2,') ',A)
99007 FORMAT (10X,A10,L10)
99008 FORMAT (A10,I10)
99009 FORMAT (A10,I10,5(:,'  (',A,')'))
99010 FORMAT (A10,F10.5)
99011 FORMAT (A10,A,A)
99012 FORMAT (10X,A10,I10)
99013 FORMAT (10X,A10,I10,5(:,'  (',A,')'))
99014 FORMAT (10X,A10,4(F10.5,A))
99015 FORMAT (10X,A10,' ( ',F6.3,',',F6.3,',',F6.3,' )')
99016 FORMAT (10X,A10,A,A)
      END
