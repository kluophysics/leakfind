C*==aes.f    processed by SPAG 6.70Rc at 15:53 on 19 Dec 2016
      SUBROUTINE AES(IPRINT,TSST,MSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,
     &               KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,
     &               BCORS,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of     aes  -  spectra                         *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   * run SSITE always for the  NL+1  angular momentum expansion       *
C   *                                                                  *
C   * jm 10.05.1999                                                    *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EFERMI
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R2DRDI
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_ANGMOM,ONLY:NKM,NMEMAX,NKMMAX,NLMAX,CGC,NLINQ,NKMQ,NLQ
      USE MOD_SITES,ONLY:NQ,IQAT
      USE MOD_FILES,ONLY:DATSET,SYSTEM,LDATSET,LSYSTEM,IFILBUILDBOT,
     &    WRBUILDBOT
      USE MOD_TYPES,ONLY:SOCTL,NTMAX,LCXRAY,NCXRAY,IMT,NLT,NT,LTXT_T,
     &    TXT_T,CTL
      USE MOD_CONSTANTS,ONLY:RY_EV,C_AU
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='AES')
      INTEGER IFILCBWF0
      PARAMETER (IFILCBWF0=60)
C
C Dummy arguments
C
      INTEGER IPRINT,ITXRAY,NCSTMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 AMECIF(:,:,:),AMECIG(:,:,:),DE,DEME,EE,EI,EPSD,IMEFIN,JMC,
     &       MJ,P3,Q4,REME,RINT(:),RJ,SK,WA,WB,WC,WD,WE,WF,WGTE,XNORM(2)
      LOGICAL CALCINT,GETIRRSOL,PRINTAME,PRINTMELE,READMELE
      CHARACTER*2 CL
      COMPLEX*16 DIFME(:,:,:),EA,EB,EBOT,EC,ED,EME(:),ERYD,ETABFIN(:,:),
     &           IAES(:,:,:),ME(:,:,:,:,:),MMED(:,:,:),MMEE(:,:,:),
     &           MMETILDD(:,:,:),MMETILDE(:,:,:),P,RME,
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAU(:,:),TAUB(:,:),TAUC(:,:),
     &           TAUTLEED(:,:),WAES(2),WRK1(:,:),WRK2(:,:),WTMP(2)
      CHARACTER*80 FILNAM,SPEC
      INTEGER I,IA_ERR,IB,IC,ICST,ICSTME,ID,IE,IEB,IEC,IED,IEME,IEME30,
     &        IEME40,IFIL,IFILB,IFILC,IFILD,IFILFIN,IFILME,IFILVAL,II,
     &        IKMB,IKMC,IKMD,IL,ILAM,ILAMP,IM,IQ,IS,IT,ITX,J,K,KAP,L,
     &        LAM,LAM1,LAM2,LAM3,LAMMAGCPL(:),LAMP,LAMP1,LAMP2,LAMP3,
     &        LFN,LSP,LSUBSH(0:4),MJM05,MS,N,NCST,NE,NEBOT,NEFIN,NEINT,
     &        NEME,NEMEINP,NEMEMAX,NKMB,NKMC,NKMD,NL,NLM,NQ9,NT9,
     &        NTXRSGRP
      CHARACTER*1 SHELL(5)
      CHARACTER*4 STRFILE,STRME(12)
      CHARACTER*3 SUBSH(0:4),SUBSHP(0:4)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IAES,MMED,MMEE,RINT,TAUB,TAUC,MMETILDD,MMETILDE
      ALLOCATABLE TAUTLEED,DIFME,LAMMAGCPL,ME,AMECIF,AMECIG
      ALLOCATABLE ETABFIN,TAU,WRK1,WRK2,EME
C
      DATA STRME/'ME1','ME2','ME3','ME4','ME5','ME6','ME7','ME8','ME9',
     &     'ME10','ME11','ME12'/,IA_ERR/0/
      DATA SHELL/'K','L','M','N','O'/
      DATA LSUBSH/1,3,3,3,3/
      DATA SUBSH/'1  ','2,3','4,5','6,7','8,9'/
      DATA SUBSHP/'1  ','23 ','45 ','67 ','89 '/
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATE (IAES(NCSTMAX,2,(2*NEMAX-1)),RINT(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> IAES'
C
      ALLOCATE (MMED(NKMMAX,NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> MMED'
C
      ALLOCATE (MMEE(NKMMAX,NKMMAX,NKMMAX),LAMMAGCPL(NKMMAX),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> MMEE'
C
      ALLOCATE (TAUB(NKMMAX,NKMMAX),TAUC(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> TAUC'
C
      ALLOCATE (MMETILDD(NKMMAX,NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> MMETILDD'
C
      ALLOCATE (MMETILDE(NKMMAX,NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> MMETILDE'
C
      ALLOCATE (TAUTLEED(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> TAUTLEED'
C
      ALLOCATE (DIFME(NKMMAX,NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> DIFME'
C
      ALLOCATE (AMECIF(NKMMAX,NKMMAX,2*NLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> AMECIF'
C
      ALLOCATE (AMECIG(NKMMAX,NKMMAX,2*NLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> AMECIG'
C
      ALLOCATE (ETABFIN((2*NEMAX-1),2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> ECORTAB'
C
      ALLOCATE (TAU(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> NKMQVAL'
C
      ALLOCATE (WRK1(NKMMAX,NKMMAX),WRK2(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> WRK2'
C
C==================================================================start
      WRITE (6,99002)
C-----------------------------------------------------reading from input
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      CALL SECTION_SET_INTEGER('IT',ITXRAY,1,0)
      IT = ITXRAY
C
      CALL SECTION_SET_INTEGER('NEME',NEMEINP,10,0)
      NEMEMAX = NEMEINP
C
      ALLOCATE (ME(NKMMAX,NKMMAX,NKMMAX,NEMEMAX,NEMEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> ME'
C
      ALLOCATE (EME(NEMEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:aes -> ME'
C
      CALL SECTION_FIND_KEYWORD('PRINTME',PRINTMELE)
C
      CALL SECTION_FIND_KEYWORD('PRINTAME',PRINTAME)
C
      CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
      IQ = IQAT(1,IT)
      IM = IMT(IT)
C
      IF ( LDATSET.NE.0 ) THEN
         FILNAM = DATSET(1:LDATSET)//TXT_T(IT)(1:LTXT_T(IT))//'.'
      ELSE
         FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'.'
      END IF
      LFN = LDATSET + LTXT_T(IT) + 1
C
      FILNAM = FILNAM(1:LFN)//SHELL(NCXRAY(IT))
      LFN = LFN + 1
      SPEC = '  '//SHELL(NCXRAY(IT))
C
      IF ( NCXRAY(IT).NE.1 ) THEN
         SPEC = SPEC(1:3)//SUBSH(LCXRAY(IT))
         LSP = 3 + LSUBSH(LCXRAY(IT))
         FILNAM = FILNAM(1:LFN)//SUBSHP(LCXRAY(IT))
         LFN = LFN + MAX(1,LSUBSH(LCXRAY(IT))-1)
      ELSE
         LSP = 3
      END IF
C
      FILNAM = FILNAM(1:LFN)//'.aes'
      SPEC = SPEC(1:LSP)//' - aes spectrum of '//TXT_T(IT)(1:LTXT_T(IT))
     &       //' in  '//SYSTEM(1:LSYSTEM)
C
      OPEN (UNIT=7,FILE=FILNAM)
      WRITE (6,'(10x,a,a,/)') 'spec-file :  ( 7) ',FILNAM
      WRITE (6,'(a)') SPEC
      WRITE (6,99003) IT,NCXRAY(IT),LCXRAY(IT)
C
C ----------------------------- increase angular momentum expansion by 1
C --------------------- to deal with final states with l_fin = l_ini + 1
      NL = 0
      DO IQ = 1,NQ
         NLQ(IQ) = NLQ(IQ) + 1
         NL = MAX(NL,NLQ(IQ))
         IF ( NLQ(IQ).GT.NLMAX ) STOP 'in <aes> nlq(iq)>>nlmax'
         NKMQ(IQ) = 2*NLQ(IQ)**2
         NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
      END DO
      DO ITX = 1,NT
         NLT(ITX) = NLQ(IQAT(1,ITX))
         CTL(IT,NLT(ITX)) = CTL(IT,NLT(ITX)-1)
         SOCTL(IT,NLT(ITX)) = SOCTL(IT,NLT(ITX)-1)
      END DO
      LAM = 0
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
            JMC = DBLE(L) + SK/2D0
            IF ( ABS(MJ).LT.L ) THEN
               LAMMAGCPL(LAM) = NINT(2*L*(JMC+0.5D0)+JMC+MJ+1)
            ELSE
               LAMMAGCPL(LAM) = 0
            END IF
         END DO
      END DO
      NLM = NL**2
      NKM = 2*NLM
      NKMB = 2*(NL-1)**2
      NKMC = NKMB
      NKMD = NKM
C ======================================================================
C read tau to get vb-info
C ======================================================================
      IFILVAL = 99
      OPEN (UNIT=IFILVAL,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(NKMMAX*NKMMAX*2*8))
C
      REWIND 9
      NE = 0
      ITX = 0
C
      CALL READTAU(9,ERYD,0,NE,TAU,0,NT,.TRUE.,WRK1,WRK2,0,NQ,0,NKMMAX,
     &             NKMMAX,IPRINT)
C
      IF ( NE.EQ.0 ) THEN
         WRITE (6,*) ' something is wrong with the TAU-file: NE=0'
         STOP ' in  <AES>'
      END IF
C
      DO IE = 1,NE
         CALL READTAU(9,ERYD,IE,NE,TAU,IT,NT,.TRUE.,WRK1,WRK2,0,NQ,1,
     &                NKMB,NKMMAX,IPRINT)
         ETAB(IE,1) = ERYD
         ETAB(IE,2) = ERYD
         WRITE (IFILVAL,REC=IE) TAU
      END DO
C
      REWIND 9
      READ (9,'(////,10X,I3,5X,I3)') NT9,NQ9
      IF ( NQ.NE.NQ9 ) STOP ' in <AES>:   NQ <> NQ9'
      DO I = 1,(NT9+NQ9)
         READ (9,*)
      END DO
C
      REWIND 9
      EBOT = ETAB(1,1)
C==================================================================
      CALL WRHEAD(7,FILNAM,'SP-AES    ',NE)
C
      NTXRSGRP = 1
      WRITE (7,99012) 'NTXRSGRP  ',NTXRSGRP
      WRITE (7,99013) 'SPECTRUM  ',SPEC(1:LSP)
      WRITE (7,99012) 'IT        ',IT
      WRITE (7,99012) 'NCXRAY    ',NCXRAY(IT)
      WRITE (7,99012) 'LCXRAY    ',LCXRAY(IT)
C ======================================================================
C prepare core states
C corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecoreco
      NCST = 4*LCXRAY(IT) + 2
      IF ( NCST.GT.NCSTMAX ) THEN
         WRITE (6,99001) LCXRAY(IT),NCSTMAX
         STOP
      END IF
C
      CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &          IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
      DO IFIL = 6,7
         WRITE (IFIL,99006) NCST,(NKPCOR(ICST),ICST=1,NCST)
         WRITE (IFIL,99007)
      END DO
      DO ICST = 1,NCST
C
         DO K = 1,NKPCOR(ICST)
            DO N = 1,JRWS(IM)
               RINT(N) = R2DRDI(N,IM)
     &                   *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
            END DO
            CALL RRADINT(IM,RINT,XNORM(K))
         END DO
C
         DO IFIL = 6,7
            WRITE (IFIL,99008) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                         (2*MM05COR(ICST)+1),IKMCOR(ICST,1),
     &                         XNORM(1),ECOR(ICST),ECOR(ICST)*RY_EV,
     &                         SZCOR(ICST),IZERO(ICST)
            IF ( NKPCOR(ICST).EQ.2 ) WRITE (IFIL,99009) IKMCOR(ICST,2),
     &           XNORM(2)
         END DO
C
      END DO
C   corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecore
      WRITE (6,99004) NE
      WRITE (7,99005) NE,1
C  =====================================================================
C         set up energy table for matrix elements
C  =====================================================================
      NEME = MIN(NEMEINP,NE)
      DEME = (EFERMI-DREAL(EBOT))/DBLE(NEME-1)
      CALCINT = .FALSE.
      GETIRRSOL = .FALSE.
      DO IE = 1,NEME
         EME(IE) = EBOT + DEME*DBLE(IE-1)
      END DO
C  =====================================================================
C                    opening of files
C  =====================================================================
      IFILFIN = 98
      IFILD = IFILCBWF0
      IFILB = IFILCBWF0 + 1
      IFILC = IFILCBWF0 + 2
      OPEN (UNIT=IFILFIN,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(4+2*(2+4*NRMAX))*2*8)
      OPEN (UNIT=IFILD,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(4+2*(2+4*NRMAX))*2*8)
      OPEN (UNIT=IFILB,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(4+2*(2+4*NRMAX))*2*8)
      OPEN (UNIT=IFILC,STATUS='scratch',FORM='unformatted',
     &      ACCESS='direct',RECL=(4+2*(2+4*NRMAX))*2*8)
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C                   read TAU for initial states
C  tautautautautautautautautautautautautautautautautautautautautautautau
      CALL CINIT(NKMMAX**2,TAU)
C
      WRITE (6,'(/,5x,A)') 'Reading TAU for valence states from file 9'
C
C      DO IE = 1,NE
C         CALL READTAU(9,ERYD,IE,NE,TAU,1,NT,.TRUE.,WRK1,WRK2,0,
C     &        NQ,1,NKMB,NKMMAX,IPRINT)
C         WRITE (IFILVAL,REC=IE) TAU
C      END DO
C
      DE = DREAL(ETAB(2,1)) - DREAL(ETAB(1,1))
      WRITE (6,99014) NE,DREAL(ETAB(1,1)),DREAL(ETAB(1,1))*RY_EV,
     &                DREAL(ETAB(NE,1)),DREAL(ETAB(NE,1))*RY_EV,DE,
     &                DE*RY_EV
C  tautautautautautautautautautautautautautautautautautautautautautautau
C
C=======================================================================
C               Calculate angular matrix elements
C  angularangularangularangularangularangularangularangularangularangula
C
      WRITE (6,*) 'Angular matrix elements will be calculated'
C
      CALL AMECOUL(AMECIG,AMECIF,IPRINT,NLMAX,NKMMAX)
C
      IF ( PRINTAME ) THEN
         DO IL = 1,2*NLMAX
            WRITE (6,*) 'Angular matrix elements for LR=',IL
            CALL RMATSTRUCT('AMECIG',AMECIG(1,1,IL),NKMD,NKMMAX,0,0,1,
     &                      1D-8,6)
            CALL RMATSTRUCT('AMECIF',AMECIF(1,1,IL),NKMD,NKMMAX,0,0,1,
     &                      1D-8,6)
         END DO
      END IF
C  angularangularangularangularangularangularangularangularangularangula
C
C ======================================================================
C               calculate spr-aes-spectrum
C ======================================================================
      WRITE (6,*)
      WRITE (6,*) 'CALCULATING SPR-AES SPECTRUM'
      WRITE (6,*)
      CALL CINIT(NCSTMAX*2*(2*NEMAX-1),IAES)
Ccorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecor
C
      DO ICST = 1,NCST
         WRITE (6,*) '       core state ',ICST
         WRITE (6,*) '       -------------'
         CALL CINIT(NKMMAX*NKMMAX*NKMMAX*NEMEMAX*NEMEMAX,ME)
         CALL CINIT(NKMMAX*NKMMAX*NKMMAX,MMED)
         CALL CINIT(NKMMAX*NKMMAX*NKMMAX,MMEE)
C     reading me from file if...
         IFILME = 23
         INQUIRE (FILE='ME',EXIST=READMELE)
         IF ( READMELE .EQV. .TRUE. ) THEN
            STRFILE = 'ME'
            GOTO 50
         END IF
         STRFILE = STRME(ICST)
         INQUIRE (FILE=STRFILE,EXIST=READMELE)
 50      CONTINUE
         IF ( READMELE ) PRINTMELE = .NOT.READMELE
         IF ( READMELE ) OPEN (UNIT=IFILME,FILE=STRFILE)
         IF ( READMELE ) THEN
            WRITE (6,*)
            WRITE (6,*) '     reading meles from file !!!!!!!!!!! '
            DO I = 1,6
               READ (IFILME,*)
            END DO
            IEC = 0
            IEB = 0
 60         CONTINUE
            IEC = IEC + 1
            IF ( IEC.GT.NEME ) GOTO 120
 80         CONTINUE
            IEB = IEB + 1
            DO I = 1,2
               READ (IFILME,*)
            END DO
            DO J = 1,1000000
               READ (IFILME,99011,ERR=100) ICSTME,IB,IC,ID,RME
               IF ( ICSTME.EQ.ICST ) ME(IB,IC,ID,IEB,IEC) = RME
            END DO
 100        CONTINUE
            IF ( IEC.LE.NEME ) THEN
               IF ( IEB.LT.NEME ) GOTO 80
               IEB = 0
               GOTO 60
            END IF
 120        CONTINUE
            CLOSE (IFILME)
C
            WRITE (6,99015) NEME
C
            WRITE (6,99016) DREAL(EBOT),DREAL(EBOT)*RY_EV,EFERMI,
     &                      EFERMI*RY_EV,DEME,DEME*RY_EV
            CLOSE (IFILME)
         ELSE
C=======================================================================
C                Calculation of matrix elements
CmememmememeCmememmemememmememememememememmememememmememememmememememmem
            WRITE (6,*)
            WRITE (6,*) ' CALCULATING ME for (IEA,IEB,IEC,IED)'
C
            DO IEC = 1,NEME
               EC = EBOT + DEME*(IEC-1)
               CALL SSITE(1,1,IFILC,CALCINT,GETIRRSOL,EC,P,IPRINT,NKM,
     &                    TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
               DO IEB = 1,NEME
                  EB = EBOT + DEME*(IEB-1)
                  EA = ECOR(ICST)
                  EE = DREAL(EB+EC-EA)
                  EI = DIMAG(EB)
                  ED = DCMPLX(EE,EI)
                  IF ( IPRINT.GE.0 ) WRITE (6,99010) ICST,IEB,IEC,
     &                 IEC + IEB - 1
C
                  CALL SSITE(1,1,IFILB,CALCINT,GETIRRSOL,EB,P,IPRINT,
     &                       NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
                  CALL SSITE(1,1,IFILD,CALCINT,GETIRRSOL,ED,P,IPRINT,
     &                       NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
                  CALL MECOUL(ICST,GCOR,FCOR,MM05COR,NKPCOR,IKMCOR,
     &                        LCXRAY,1,NKMB,IFILB,1,NKMC,IFILC,1,NKMD,
     &                        IFILD,IT,ME,IEB,IEC,AMECIG,AMECIF,NEMEMAX,
     &                        NCSTMAX)
               END DO
            END DO
C     mememmemememmememememememememmememememmememememmememememmeme
C     writing of matrix elements
C     mememmemememmememememememememmememememmememememmememememmeme
            IF ( PRINTMELE ) THEN
               IFILME = 24
               STRFILE = STRME(ICST)
               OPEN (UNIT=IFILME,FILE=STRFILE)
               WRITE (6,*)
               WRITE (6,*) 
     &                   ' matrix elements m(c,v,v,l) output in file me'
               WRITE (IFILME,'(a)') 
     &                          '    coulomb matrix elements m(c,v,v,l)'
               WRITE (IFILME,'(10x,a,i2)') 'k_c   = 1,',NCST
               WRITE (IFILME,'(10x,a,i3)') 'k_v   = 1,',NKMB
               WRITE (IFILME,'(10x,a,i3)') 'k_l   = 1,',NKMD
               WRITE (IFILME,*)
               WRITE (IFILME,'(2a)') 
     &                       '*****************************************'
     &                       ,'*****************************'
               DO IEC = 1,NEME
                  DO IEB = 1,NEME
                     EC = EBOT + DEME*DBLE(IEC-1)
                     EB = EBOT + DEME*DBLE(IEB-1)
                     WRITE (IFILME,'(a6,f10.6,a6,f10.6)') ' eb = ',
     &                      DREAL(EB),' ec = ',DREAL(EC)
                     WRITE (IFILME,'(6x,a17)') 'k_c k_vb k_vc k_l'
                     DO ID = 1,NKMD
                        DO IC = 1,NKMC
                           DO IB = 1,NKMB
                              IF ( CDABS(ME(IB,IC,ID,IEB,IEC))
     &                             .GT.1.D-14 ) WRITE (IFILME,99011)
     &                             ICST,IB,IC,ID,ME(IB,IC,ID,IEB,IEC)
                           END DO
                        END DO
                     END DO
                     WRITE (IFILME,'(2a)')
     &                       '*****************************************'
     &                      ,'*****************************'
                  END DO
               END DO
               CLOSE (IFILME)
            END IF
         END IF
C     mememmemememmememememememememmememememmememememmememememmeme
C         IF ( (2*NE-1).GT.NEMAX ) THEN
C            WRITE (6,*) ' < aes > : increase nemax to ',2*NE - 1
C            STOP
C         END IF
         IMEFIN = DIMAG(ETAB(1,1))
         NEFIN = 2*NE - 1
         DO IE = NEFIN,1, - 1
            ETABFIN(IE,1) = DCMPLX(2*EFERMI-ECOR(ICST),IMEFIN)
     &                      - DCMPLX((NEFIN-IE)*DE,0.D0)
         END DO
C     iebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebie
         DO IED = 1,2*NE - 1
            WRITE (6,*) '  IED=',IED,'  ED=',DREAL(ETABFIN(IED,1))
            IF ( IED.LE.NE ) THEN
               NEBOT = 1
               NEINT = IED
            ELSE IF ( IED.GT.NE ) THEN
               NEBOT = IED + 1 - NE
               NEINT = NE
            END IF
            ED = ETABFIN(IED,1)
C     iecieciecieciecieciecieciecieciecieciecieciecieciecieciecieciecie
            DO IEB = NEBOT,NEINT
               IEC = IED + 1 - IEB
C
               IF ( (IEB.EQ.NEBOT) .OR. (IEB.EQ.NEINT) ) THEN
                  WGTE = DE
               ELSE
                  WGTE = 2D0*DE
               END IF
               EB = EBOT + DBLE(IEB-1)*DE
               EC = EBOT + DBLE(IEC-1)*DE
C-----------------interpolate the calculated MELEs to the vb-mesh
C     interpolationinterpolationinterpolationinterpolationinterpolation
               IEME30 = NEME
               IEME40 = NEME
               DO IEME = 1,NEME
                  REME = DREAL(EME(IEME))
                  IF ( (IEME30.EQ.NEME) .AND. (REME.GT.DREAL(EB)) )
     &                 IEME30 = MAX(2,IEME-1)
                  IF ( (IEME40.EQ.NEME) .AND. (REME.GT.DREAL(EC)) )
     &                 IEME40 = MAX(2,IEME-1)
               END DO
               IEME30 = MIN(IEME30,NEME-1)
               IEME40 = MIN(IEME40,NEME-1)
               P3 = (DREAL(EB)-DREAL(EME(IEME30)))/DEME
               Q4 = (DREAL(EC)-DREAL(EME(IEME40)))/DEME
C
               WA = Q4*(Q4-1D0)/2D0
               WB = P3*(P3-1D0)/2D0
               WC = 1D0 + P3*Q4 - P3**2 - Q4**2
               WD = P3*(P3-2D0*Q4+1D0)/2D0
               WE = Q4*(Q4-2D0*P3+1D0)/2D0
               WF = P3*Q4
               DO IKMB = 1,NKMB
                  DO IKMC = 1,NKMC
                     DO IKMD = 1,NKMD
                        MMED(IKMB,IKMC,IKMD)
     &                     = ME(IKMB,IKMC,IKMD,IEME30+0,IEME40-1)
     &                     *WA + ME(IKMB,IKMC,IKMD,IEME30-1,IEME40+0)
     &                     *WB + ME(IKMB,IKMC,IKMD,IEME30+0,IEME40+0)
     &                     *WC + ME(IKMB,IKMC,IKMD,IEME30+1,IEME40+0)
     &                     *WD + ME(IKMB,IKMC,IKMD,IEME30+0,IEME40+1)
     &                     *WE + ME(IKMB,IKMC,IKMD,IEME30+1,IEME40+1)*WF
                        MMEE(IKMB,IKMC,IKMD)
     &                     = ME(IKMB,IKMC,IKMD,IEME40+0,IEME30-1)
     &                     *WA + ME(IKMB,IKMC,IKMD,IEME40-1,IEME30+0)
     &                     *WB + ME(IKMB,IKMC,IKMD,IEME40+0,IEME30+0)
     &                     *WC + ME(IKMB,IKMC,IKMD,IEME40+1,IEME30+0)
     &                     *WD + ME(IKMB,IKMC,IKMD,IEME40+0,IEME30+1)
     &                     *WE + ME(IKMB,IKMC,IKMD,IEME40+1,IEME30+1)*WF
                     END DO
                  END DO
               END DO
C     interpolationinterpolationinterpolationinterpolationinterpolation
Cttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
C                   get TAU for final states
C                single-scatterer approximation
C     tautautautautautautautautautautautautautautautautautautautautautau
               ERYD = ETABFIN(IED,1)
               CALL CINIT(NKMMAX*NKMMAX*NTMAX,TSST)
               CALL SSITE(1,1,IFILFIN,CALCINT,GETIRRSOL,ERYD,P,IPRINT,
     &                    NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
               CALL CINIT(NKMMAX*NKMMAX,TAUTLEED)
               DO I = 1,NKMD
                  DO II = 1,NKMD
                     TAUTLEED(I,II) = TSST(I,II,IT)
                  END DO
               END DO
C     tautautautautautautautautautautautautautautautautautautautautautau
               READ (IFILVAL,REC=IEB) TAUB
               READ (IFILVAL,REC=IEC) TAUC
C
C-------tild
               CALL CINIT(NKMMAX*NKMMAX*NKMMAX,MMETILDD)
               CALL CINIT(NKMMAX*NKMMAX*NKMMAX,MMETILDE)
               DO IKMB = 1,NKMB
                  DO IKMC = 1,NKMC
                     DO ILAM = 1,NKMD
                        DO ILAMP = 1,NKMD
                           MMETILDD(IKMB,IKMC,ILAM)
     &                        = MMETILDD(IKMB,IKMC,ILAM)
     &                        + DCONJG(TAUTLEED(ILAMP,ILAM))
     &                        *MMED(IKMB,IKMC,ILAMP)
                           MMETILDE(IKMB,IKMC,ILAM)
     &                        = MMETILDE(IKMB,IKMC,ILAM)
     &                        + DCONJG(TAUTLEED(ILAMP,ILAM))
     &                        *MMEE(IKMB,IKMC,ILAMP)
                        END DO
                     END DO
                  END DO
               END DO
C-------tild
               CALL CINIT(NKMMAX*NKMMAX*NKMMAX,DIFME)
               DO LAM1 = 1,NKMB
                  DO LAM2 = 1,NKMC
                     DO LAM3 = 1,NKMD
                        DIFME(LAM1,LAM2,LAM3)
     &                     = (MMETILDD(LAM1,LAM2,LAM3)
     &                     -MMETILDE(LAM2,LAM1,LAM3))
                     END DO
                  END DO
               END DO
               CALL CINIT(2,WTMP)
               DO LAMP2 = 1,NKMB
                  DO LAMP3 = 1,NKMB
                     DO LAMP = 1,NKMC
                        DO LAMP1 = 1,NKMC
                           CALL CINIT(2,WAES)
C     sumation over lam(LEED)
                           DO LAM1 = 1,NKMD
                              DO LAM2 = 1,NKMD
C
                                 IF ( (LAM1.EQ.LAM2) .OR. 
     &                                (LAM2.EQ.LAMMAGCPL(LAM1)) ) THEN
                                    WAES(1) = WAES(1) + CGC(LAM1,1)
     &                                 *CGC(LAM2,1)
     &                                 *DIFME(LAMP2,LAMP,LAM1)
     &                                 *DCONJG(DIFME(LAMP3,LAMP1,LAM2))
                                    WAES(2) = WAES(2) + CGC(LAM1,2)
     &                                 *CGC(LAM2,2)
     &                                 *DIFME(LAMP2,LAMP,LAM1)
     &                                 *DCONJG(DIFME(LAMP3,LAMP1,LAM2))
C
C     sumation over lam(LEED)
                                 END IF
                              END DO
                           END DO
                           DO MS = 1,2
                              WTMP(MS) = WTMP(MS)
     &                           + DIMAG(TAUC(LAMP,LAMP1))
     &                           *DIMAG(TAUB(LAMP2,LAMP3))*WAES(MS)
                           END DO
C
C     sumation over lam(valence)
                        END DO
                     END DO
                  END DO
               END DO
               DO MS = 1,2
                  IAES(ICST,MS,IED) = IAES(ICST,MS,IED) + WGTE*WTMP(MS)
               END DO
C     iebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebie
C             end if
            END DO
C     iecieciecieciecieciecieciecieciecieciecieciecieciecieciecieciecie
            EPSD = 16.D0*(DREAL(ETABFIN(IED,1))+C_AU**2)
     &             /(2.D0*DREAL(ETABFIN(IED,1))+C_AU**2)
            DO MS = 1,2
               IAES(ICST,MS,IED) = IAES(ICST,MS,IED)*EPSD
            END DO
C
            IF ( IPRINT.GT.0 ) THEN
               WRITE (6,*) IAES(ICST,1,IED)
               WRITE (6,*) IAES(ICST,2,IED)
            END IF
Cwrite spectrum
            WRITE (7,'(1x,e12.5,2x,e12.5,a)') DREAL(ETABFIN(IED,1)),
     &             DREAL(ETABFIN(IED,1)) - DREAL(ETABFIN(NEFIN,1))
            WRITE (7,'(20x,6(1x,e18.10))') (IAES(ICST,IS,IED),IS=1,2),
     &             IAES(ICST,1,IED) + IAES(ICST,2,IED)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT .AND. IED.LE.3 ) WRITE (IFILBUILDBOT,99017)
     &           ROUTINE(1:LEN_TRIM(ROUTINE)),IT,ICST,IED,ETABFIN(IED,1)
     &           ,(IAES(ICST,IS,IED),IS=1,2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         END DO
C
Ccorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecor
      END DO
C
      DEALLOCATE (IAES,MMED,MMEE,RINT,TAUB,TAUC,MMETILDD,MMETILDE)
      DEALLOCATE (TAUTLEED,DIFME,LAMMAGCPL,ME,AMECIF,AMECIG)
      DEALLOCATE (ETABFIN,TAU,WRK1,WRK2,EME)
C
99001 FORMAT (//,60('*'),/,10x,' stop in <aes>',/,10x,'lcxray=',i3,
     &        ' too large for ncstmax=',i3)
99002 FORMAT (' ',//,10x,62('*'),/,10x,'*',60x,'*',/,10x,
     &  '*        ****   *****            ***    ******   ****        *'
     &  ,/,10x,
     &  '*       *    *  *    *          *   *   *       *    *       *'
     &  ,/,10x,
     &  '*       *       *    *         *     *  *       *            *'
     &  ,/,10x,
     &  '*        ****   *****    ***   *******  *****    ****        *'
     &  ,/,10x,
     &  '*            *  *              *     *  *            *       *'
     &  ,/,10x,
     &  '*       *    *  *              *     *  *       *    *       *'
     &  ,/,10x,
     &  '*        ****   *              *     *  ******   ****        *'
     &  ,/,10x,'*',60x,'*',/,10x,62('*'),//)
C----------------------------------------------------------------formats
99003 FORMAT (/,4X,' CORE QUANTUM-NUMBERS  FOR  IT=',I2,':   N=',I2,
     &        '  L=',I2,/)
99004 FORMAT (' ',//,10X,'NUMBER OF TABULATED ENERGIES',I5,/)
99005 FORMAT (/,' number of energies',i5,/,' output format ifmt',i5,/)
99006 FORMAT (//,' CORE STATES :',//,' NCST:  ',I4,/,' NKPCOR:',20I4)
99007 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(Ev)    <SIGMA_z>  I0')
99008 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99009 FORMAT (22X,I4,F12.6)
99010 FORMAT (28X,'(',3(I3,','),I3,')')
99011 FORMAT (I9,I5,I5,I4,2X,2(G21.14))
99012 FORMAT (A10,I10)
99013 FORMAT (A10,A)
99014 FORMAT (' TAU:  NUMBER OF ENERGIES READ IN  NE=',I3,/,
     &        '       FROM   E(1)  =',F8.3,' RY    (',F8.3,' EV)',/,
     &        '       TO     E(NE) =',F8.3,' RY    (',F8.3,' EV)',/,
     &        '       WITH   DE    =',F8.3,' RY    (',F8.3,' EV)')
99015 FORMAT (' ME: WAVE FUNCTIONS READ IN FOR NEME=',I3,' ENERGIES')
99016 FORMAT (' ME:   FROM   EBOT  =',F8.3,' RY    (',F8.3,' EV)',/,
     &        '       TO     EFERMI=',F8.3,' RY    (',F8.3,' EV)',/,
     &        '       WITH   DEME  =',F8.3,' RY    (',F8.3,' EV)')
C----------------------------------------------------------------formats
99017 FORMAT ('# BUILDBOT: ',A,': AES intensity for IT =',I5,'  ICST =',
     &        I3,/,'#',10X,'energy  IED =',I5,' ETABFIN = ',2F10.6,/,
     &        (1PE22.14))
      END
