C*==nraps.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NRAPS(TSST,MSST,SSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,
     &                 KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,
     &                 BCORS,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of  spin-resolved  APS - spectra               *
C   *                                                                  *
C   * 10/06/96  HE+VP                                                  *
C   * 28/02/97  VP    The {ls,l's'} cross section introduced (WLS)     *
C   *                 The matrix elements and the cross sections are   *
C   *                 written in GNUPLOT-format files if OUTGNU = TRUE *
C   *                 (l,s) - resolved spectra are calculated if       *
C   *                 CALCLS = TRUE                                    *
C   *                                                                  *
C   * 07/04/97  VP    the DOS above E_F is broadened with a linear     *
C   *                 energy-dependent Lorenzian BEFORE calculating    *
C   *                 spectra if BRDDOS = TRUE                         *
C   *                 the energy range is extended below E_F and above *
C   *                 the calculated limit                             *
C   *                                                                  *
C   * 29/04/97  VP    the OUTGNU, CALCLS and GAMMA (broadening         *
C   *                 parameter) are read in from the input file       *
C   *                 ( 'task aps' - line )                            *
C   * 12/05/2000 VP   MD and ME arrays restructured to 6 dimensions,   *
C   *                 i.e., IL3 and IL4 are combined in one index,     *
C   *                 mapped by the function INDIL(IL3,IL4)            *
C   *                 the fundamental constants PI, E**2 and so on     *
C   *                 were explicitly introduced                       *
C   *                                                                  *
C   * NOTE: NL, NLQ, VT and BT  in  MODULES  temporarily overwritten   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,IFILCBWF,IFILBUILDBOT,WRBUILDBOT
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EMAX,EFERMI
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SITES,ONLY:NQ
      USE MOD_ANGMOM,ONLY:TXT_L,NMEMAX,NKMMAX,NLMAX,NL,NSPIN
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,INFO,SYSTEM,LDATSET,RECLNGWF,
     &    FOUND_REAL
      USE MOD_TYPES,ONLY:NTMAX,LCXRAY,NCXRAY,BT,VT,IMT,NT,LTXT_T,TXT_T,
     &    NLT
      USE MOD_CONSTANTS,ONLY:RY_EV,PI
      IMPLICIT NONE
C*--NRAPS47
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NRAPS')
      INTEGER NEMEMAX
      PARAMETER (NEMEMAX=20)
      REAL*8 TWOPIOVHBAR,FOURPIE,FOURPIESQ
      PARAMETER (TWOPIOVHBAR=2D0*PI,FOURPIE=4D0*PI**2*2D0,
     &           FOURPIESQ=FOURPIE**2)
C
C Dummy arguments
C
      INTEGER ITXRAY,NCSTMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 ADDRANGE,APSINT,BS(:,:),BT0(:,:),DE,DEBOT1,
     &       DECIAPS(:,:,:,:,:),DEF1,DEL,DELE,DEME,DOS(:,:,:,:),E1,E2,
     &       E3,E4,EA,EB,EFIN,EONSET,GAMMAV,HEADER,IAPS(:,:,:),
     &       MD(:,:,:,:,:,:,:),MDXME,ME(:,:,:,:,:,:,:),P3,Q4,REME,
     &       RINT(:),RNORM,RPW(:,:),RSUM,SDD(:,:,:,:,:,:),
     &       SDE(:,:,:,:,:,:),SEE(:,:,:,:,:,:),SG(:),SIGNIS,SL(:),
     &       TDD(:,:,:,:),TDE(:,:,:,:),TEE(:,:,:,:),TG(:),TL(:),VECBR(:)
     &       ,VGL,VS(:,:),VT0(:,:),W1W3R2,W1W4R2,WA,WB,WC,WD,WE,WF,
     &       WF1(:),WF2(:),WF3(:),WF4(:),WGTE3,WLORTAB(:),WLS(:,:,:,:),
     &       XNORM(2)
      LOGICAL BRDDOS,CALCLS,FND,OUTGNU
      COMPLEX*16 CDUM1,CINT(:),CNORM,CTGDEL,EBOT,EBRD(:),EDOS,
     &           EME(NEMEMAX),ERYD,P,WRCA(:,:,:)
      CHARACTER*2 CL
      CHARACTER*80 FILNAM,FILNAMDEC,SPEC
      INTEGER I,IA_ERR,ICST,IE,IE2,IE3,IE4,IEME,IEME30,IEME40,IFIL,
     &        IFIL34WF,IFILLEED,IL,IL2,IL3,IL4,ILA,ILAM,ILAMP,ILB,ILM,
     &        ILP,IM,IQP,IR,IREC,IRELIN,IRTOP,IS,IS1,IS13,IS14,IS2,IS23,
     &        IS24,IS3,ISA,ISB,ISP,IT,ITP,J,JCST,JS,K,KGT,L1,L2,L3,L4,
     &        LAM,LAMP,LFN,LP,LSP,LSUBSH(0:4),N,NADD,NCST,NE,NE3,NEBRD,
     &        NEME,NL0,NLLEED,NLP,NLT0(:),NLTLEED(:),NSOLP,NTXRSGRP
      INTEGER INDIL
      CHARACTER*1 SHELL(5),TXTS(2)
      CHARACTER*3 STR3,SUBSH(0:4),SUBSHP(0:4)
      LOGICAL TRIANGLE
      REAL*8 WIG_3J,WIG_6J_RACAH
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE EBRD,IAPS,CINT,WRCA,RINT,VECBR,MD,ME,SG,TG,SL,TL
      ALLOCATABLE WF1,WF2,WF3,WF4,SDD,SDE,SEE,TDD,TDE,TEE,DOS,DECIAPS
      ALLOCATABLE WLS,RPW,WLORTAB,BS,VS,VT0,BT0
      ALLOCATABLE NLTLEED,NLT0
C
      DATA SHELL/'K','L','M','N','O'/
      DATA LSUBSH/1,3,3,3,3/
      DATA SUBSH/'1  ','2,3','4,5','6,7','8,9'/
      DATA SUBSHP/'1  ','23 ','45 ','67 ','89 '/
      DATA TXTS/'U','D'/
C
      TRIANGLE(L1,L2,L3) = (L1.GE.ABS(L3-L2)) .AND. (L1.LE.(L3+L2))
     &                     .AND. (MOD((L1+L2+L3),2).EQ.0)
C
      INDIL(IL3,IL4) = IL3 + NLMAX*(IL4-1)
C
      ALLOCATE (BT0(NRMAX,NTMAX),BS(NRMAX,NTMAX))
      ALLOCATE (VT0(NRMAX,NTMAX),VS(NRMAX,NTMAX))
      ALLOCATE (IAPS(2*NEMAX-1,0:NLMAX,0:2),RINT(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> IAPS'
C
      ALLOCATE (WRCA(NRMAX,2,2),CINT(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> WRCB'
C
      ALLOCATE (VECBR(NEMAX+101),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> VECBR'
C
      NEME = NEMEMAX
      ALLOCATE (MD(NLMAX,NLMAX**2,2*NLMAX,2,2,NEME,NEME))
      ALLOCATE (ME(NLMAX,NLMAX**2,2*NLMAX,2,2,NEME,NEME),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> ME'
C
      ALLOCATE (SG(NRMAX),TG(NRMAX),SL(NRMAX),TL(NRMAX))
      ALLOCATE (WF1(NRMAX),WF2(NRMAX),WF3(NRMAX),WF4(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> WF4'
C
      ALLOCATE (SDD(NLMAX,NLMAX,2,2,NEMEMAX,NEMEMAX))
      ALLOCATE (SDE(NLMAX,NLMAX,2,2,NEMEMAX,NEMEMAX))
      ALLOCATE (SEE(NLMAX,NLMAX,2,2,NEMEMAX,NEMEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> SEE'
C
      ALLOCATE (TDD(NLMAX,NLMAX,2,2),TDE(NLMAX,NLMAX,2,2))
      ALLOCATE (TEE(NLMAX,NLMAX,2,2),EBRD(NEMAX+101),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> TEE'
C
      ALLOCATE (DOS(NEMAX+101,NLMAX,2,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> DOS'
C
      ALLOCATE (DECIAPS(2*NEMAX-1,NLMAX,NLMAX,2,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraps -> DECIAPS'
C
      ALLOCATE (WLS(NLMAX,NLMAX,2,2),RPW(NRMAX,2*NLMAX),NLT0(NTMAX))
      ALLOCATE (WLORTAB(NEMAX+101),NLTLEED(NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:nraes -> NLTLEED'
C
      NEME = MIN(NEMEMAX,20)
C
      IFILLEED = 55
      OPEN (UNIT=IFILLEED,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
      IFIL34WF = 56
      OPEN (UNIT=IFIL34WF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
C ======================================================================
      WRITE (6,99002)
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      CALL SECTION_SET_INTEGER('IT',ITXRAY,1,0)
      IT = ITXRAY
C
      CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
C
      DO IR = 1,NRMAX
         RPW(IR,1) = R(IR,IM)
         DO IL = 2,2*NLMAX
            RPW(IR,IL) = RPW(IR,IL-1)*R(IR,IM)
         END DO
      END DO
C
      IF ( LDATSET.NE.0 ) THEN
         FILNAM = DATSET(1:LDATSET)//'.'//TXT_T(IT)(1:LTXT_T(IT))//'.'
         LFN = LDATSET + LTXT_T(IT) + 2
      ELSE
         FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'.'
         LFN = LTXT_T(IT) + 1
      END IF
C
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
      FILNAM = FILNAM(1:LFN)//'.aps'
      SPEC = SPEC(1:LSP)//' - APS spectrum of '//TXT_T(IT)//' in  '//
     &       SYSTEM(1:LSYSTEM)
C
      OPEN (UNIT=7,FILE=FILNAM)
      WRITE (6,'(10X,A,A,/)') 'SPEC-FILE :  ',FILNAM
C
      WRITE (6,'(A)') SPEC
      WRITE (6,99003) IT,NCXRAY(IT),LCXRAY(IT)
C
C-----------------------------------------------------------------------
C
      CALCLS = .FALSE.
      CALL SECTION_FIND_KEYWORD('CALCLS',FND)
      IF ( FND ) CALCLS = .TRUE.
C
      OUTGNU = .FALSE.
      CALL SECTION_FIND_KEYWORD('OUTGNU',FND)
      IF ( FND ) OUTGNU = .TRUE.
C
      BRDDOS = .FALSE.
      GAMMAV = -0.1D0
C
      CALL SECTION_SET_REAL('GAMMA',GAMMAV,9999D0,0)
      IF ( FOUND_REAL ) BRDDOS = .TRUE.
      IF ( GAMMAV.LE.0D0 ) BRDDOS = .FALSE.
C
      WRITE (6,*)
      WRITE (6,'(10X,A)') '********************************************'
      WRITE (6,'(10X,A,L2)') 'GNUPLOT files output                 = ',
     &                       OUTGNU
      WRITE (6,'(10X,A,L2)') '(l,s) - resolved spectra calculation = ',
     &                       CALCLS
      WRITE (6,'(10X,A,L2)') 'DOS E-dependent linear broadening    = ',
     &                       BRDDOS
      IF ( BRDDOS ) WRITE (6,'(10X,A,F6.3,A)')
     &                      '                          with gamma = ',
     &                     GAMMAV,' eV'
      WRITE (6,'(10X,A)') '********************************************'
      WRITE (6,*)
C
      ADDRANGE = 0.0D0
      NADD = 0
C
C-------------------------------------------    extend E-range below E_F
C
      IF ( BRDDOS ) THEN
         ADDRANGE = 5.0D0/RY_EV
         GAMMAV = GAMMAV/RY_EV
      END IF
      EBOT = EFERMI - ADDRANGE
      DEME = (EMAX-DREAL(EBOT))/DBLE(NEME-1)
C
      REWIND 10
      DO I = 1,6
         READ (10,99013)
      END DO
      READ (10,99014) NE
      PRINT *,'   read  NRAPS  '
C
      CALL WRHEAD(7,FILNAM,'NRAPS     ',NE)
C
      IF ( IREL.GT.1 ) STOP ' in <NRAPS>:  IREL > 1 !!!!!!!!'
C
      NTXRSGRP = 1
      WRITE (7,99018) 'NTXRSGRP  ',NTXRSGRP
      WRITE (7,99019) 'SPECTRUM  ',SPEC(1:LSP)
      WRITE (7,99018) 'IT        ',IT
      WRITE (7,99018) 'NCXRAY    ',NCXRAY(IT)
      WRITE (7,99018) 'LCXRAY    ',LCXRAY(IT)
C
C ======================================================================
C                                                     prepare core state
      NCST = 4*LCXRAY(IT) + 2
C
      IF ( NCST.GT.NCSTMAX ) THEN
         WRITE (6,99001) LCXRAY(IT),NCSTMAX
         STOP
      END IF
C
      PRINT *,'  core  '
C
      CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &          IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
      DO IFIL = 6,7
         WRITE (IFIL,99004) NCST,(NKPCOR(ICST),ICST=1,NCST)
         WRITE (IFIL,99005)
      END DO
C
      E1 = -1D+10
      DO JCST = 1,NCST
         IF ( ECOR(JCST).GT.E1 ) THEN
            E1 = ECOR(JCST)
            ICST = JCST
         END IF
      END DO
C
      DEF1 = EFERMI - E1
C
      DEBOT1 = DREAL(EBOT) - E1
C
      VGL = 0D0
      DO K = 1,NKPCOR(ICST)
         DO N = 1,JRWS(IM)
            RINT(N) = R2DRDI(N,IM)*(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
         END DO
         CALL RRADINT(IM,RINT,XNORM(K))
         IF ( XNORM(K).GT.VGL ) KGT = K
      END DO
C
      DO IR = 1,IRTOP
         RINT(IR) = R2DRDI(IR,IM)*GCOR(IR,KGT,ICST)**2
      END DO
C
      CALL RRADINT(IM,RINT,RSUM)
C
      RNORM = 1D0/SQRT(RSUM)
C
      DO IR = 1,IRTOP
         WF1(IR) = GCOR(IR,KGT,ICST)*RNORM
      END DO
C
      DO IFIL = 6,7
         WRITE (IFIL,99006) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                      (2*MM05COR(ICST)+1),IKMCOR(ICST,1),XNORM(1),
     &                      ECOR(ICST),ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                      IZERO(ICST)
         IF ( NKPCOR(ICST).EQ.2 ) WRITE (IFIL,99007) IKMCOR(ICST,2),
     &        XNORM(2)
      END DO
C
      WRITE (6,*) ' CORE :',ICST,KGT,RNORM,E1,DEF1,DEBOT1
C
C ======================================================================
C                                set up energy table for matrix elements
C                          calculate and store normalized wave functions
      NLLEED = NL + 1
      NL0 = NL
      IF ( NLLEED.GT.NLMAX ) STOP ' NRAPS:  NLLEED > NLMAX '
C
      NLTLEED(IT) = NLT(IT) + 1
      NLT0(IT) = NLT(IT)
C
      CALL RVECCOP(NRMAX*NT,VT,VT0)
      CALL RVECCOP(NRMAX*NT,BT,BT0)
C
      VS(:,:) = 0.0D0
      BS(:,:) = 0.0D0
C
C-------------------------------------------------------------------- IS
      DO IS = 1,2
C
         SIGNIS = 2.0D0*(DBLE(IS)-1.5D0)
         DO IR = 1,IRTOP
            VS(IR,IT) = VT(IR,IT) - SIGNIS*BT(IR,IT)
            BS(IR,IT) = 0D0
         END DO
C
C-------------------------------------------------------------------- IE
C----------------------------------------- deal with states close the EF
         DO IE = 1,NEME
C
            EME(IE) = EBOT + DEME*DBLE(IE-1)
            ERYD = EME(IE)
C
            P = SQRT(ERYD)
C
            NL = NLLEED
            NLT(IT) = NLTLEED(IT)
            CALL RVECCOP(NRMAX*NT,VS,VT)
            CALL RVECCOP(NRMAX*NT,BS,BT)
C
            IRELIN = IREL
            IREL = 0
            CALL NRSSITE(1,0,IFILCBWF,ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,
     &                   MEZJ)
            IREL = IRELIN
C
            NL = NL0
            NLT(IT) = NLT0(IT)
            CALL RVECCOP(NRMAX*NT,VT0,VT)
            CALL RVECCOP(NRMAX*NT,BT0,BT)
C
            DO IL = 1,NLLEED
               READ (IFILCBWF,REC=IL+(IT-1)*NLLEED) ITP,LP,NLP,NSOLP,
     &               STR3,((WRCA(IR,1,1),IR=1,IRTOP),JS=1,NSOLP)
C
               IF ( STR3.NE.'NRR' .OR. IT.NE.ITP ) THEN
                  WRITE (6,*) '##### TROUBLE reading WF:',ITP,LP,NLP,
     &                        NSOLP,STR3,WRCA(1,1,1)
                  STOP 'in <NRAPS>'
               END IF
C                  IF ( IWRREGWF.NE.0 ) WRITE (IFILSS,REC=IL+(IT-1)*NL)
C     &                 IT,L,NL,NSPIN,'NRR',
C     &                 ((ZGS(I,IL,JS),I=1,IRTOP),JS=1,NSPIN)
C
               DO N = 1,JRWS(IM)
                  CINT(N) = R2DRDI(N,IM)*WRCA(N,1,1)**2
               END DO
C
               CALL CRADINT(IM,CINT,CDUM1)
C
               CNORM = 1.0D0/CDUM1
C
               IREC = 2*((IE-1)*NLLEED+(IL-1)) + IS
               WRITE (IFIL34WF,REC=IREC) IT,IL,IS,
     &                (DREAL(WRCA(IR,1,1)*CNORM),IR=1,IRTOP)
            END DO
C
         END DO
C
C-------------------------------------------------------------------- IE
C----------------------------------------------- deal with LEED - states
         DO IE2 = 1,2*NEME - 1
C
            ERYD = EBOT + DEME*DBLE(IE2-1) + DEBOT1
            P = SQRT(ERYD)
C
            NL = NLLEED
            NLT(IT) = NLTLEED(IT)
            CALL RVECCOP(NRMAX*NT,VS,VT)
            CALL RVECCOP(NRMAX*NT,BS,BT)
C
            IRELIN = IREL
            IREL = 0
            CALL NRSSITE(1,0,IFILCBWF,ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,
     &                   MEZJ)
            IREL = IRELIN
C
            NL = NL0
            NLT(IT) = NLT0(IT)
            CALL RVECCOP(NRMAX*NT,VT0,VT)
            CALL RVECCOP(NRMAX*NT,BT0,BT)
C
            DO IL = 1,NLLEED
               READ (IFILCBWF,REC=IL+(IT-1)*NLLEED) ITP,LP,NLP,NSOLP,
     &               STR3,((WRCA(IR,1,1),IR=1,IRTOP),JS=1,NSOLP)
C
               ILM = IL**2
               CTGDEL = -1D0/(P*TSST(ILM,ILM,IT)) + DCMPLX(0D0,1D0)
               DEL = ATAN(1D0/DREAL(CTGDEL))
               CNORM = TSST(ILM,ILM,IT)*CDEXP(-DCMPLX(0D0,DEL))
C
               IREC = 2*((IE2-1)*NLLEED+(IL-1)) + IS
               WRITE (IFILLEED,REC=IREC) IT,IL,IS,
     &                (DREAL(WRCA(IR,1,1)*CNORM),IR=1,IRTOP)
            END DO
C
         END DO
C-------------------------------------------------------------------- IE
C
         WRITE (6,99008) IS,NEME
      END DO
C-------------------------------------------------------------------- IS
C
      WRITE (6,99009) DREAL(EBOT),DREAL(EBOT)*RY_EV,EMAX,EMAX*RY_EV,
     &                DEME,DEME*RY_EV,DEBOT1,DEBOT1*RY_EV,DEF1,
     &                DEF1*RY_EV
C
C ======================================================================
C                   calculate matrix elements
C
C       MD =  W3J(L1,LAM,L3) * W3J(L2,LAM,L4) * I(1,2,LAM,3,4)
C       ME =  W3J(L1,LAM,L4) * W3J(L2,LAM,L3) * I(1,2,LAM,4,3)
C
C                      sum up matrix elements
C
C       SDD =  SUM(L2,LAM)      1/(2L3+1) ....  MD**2
C       SEE =  SUM(L2,LAM)      1/(2L3+1) ....  ME**2
C
      CALL RINIT(NLMAX*NLMAX*NLMAX*2*NLMAX*2*2*NEMEMAX*NEMEMAX,MD)
      CALL RINIT(NLMAX*NLMAX*NLMAX*2*NLMAX*2*2*NEMEMAX*NEMEMAX,ME)
      CALL RINIT(NLMAX*NLMAX*2*2*NEMEMAX*NEMEMAX,SDD)
      CALL RINIT(NLMAX*NLMAX*2*2*NEMEMAX*NEMEMAX,SEE)
      CALL RINIT(NLMAX*NLMAX*2*2*NEMEMAX*NEMEMAX,SDE)
C
      L1 = LCXRAY(IT)
C
C-------------------------------------------------------------------- MD
C
C
      DO IS13 = 1,2
         DO IS24 = 1,2
C
C-------------------------------------------------------- IE3
            DO IE3 = 1,NEME
C---------------------------------------------------- IE4
               DO IE4 = 1,NEME
C
                  IE2 = IE3 + IE4 - 1
C
                  DO IL4 = 1,NL
                     L4 = IL4 - 1
                     IREC = 2*((IE4-1)*NLLEED+(IL4-1)) + IS24
                     READ (IFIL34WF,REC=IREC) ITP,ILP,ISP,
     &                     (WF4(IR),IR=1,IRTOP)
C
                     DO IL3 = 1,NL
                        L3 = IL3 - 1
                        IREC = 2*((IE3-1)*NLLEED+(IL3-1)) + IS13
                        READ (IFIL34WF,REC=IREC) ITP,ILP,ISP,
     &                        (WF3(IR),IR=1,IRTOP)
C
                        DO IL2 = 1,NLLEED
                           L2 = IL2 - 1
                           IREC = 2*((IE2-1)*NLLEED+(IL2-1)) + IS24
                           READ (IFILLEED,REC=IREC) ITP,ILP,ISP,
     &                           (WF2(IR),IR=1,IRTOP)
C
                           DO LAM = ABS(L3-L1),(L3+L1),2
                              IF ( TRIANGLE(L2,LAM,L4) ) THEN
                                 ILAM = LAM + 1
C
                                 DO IR = 1,IRTOP
                                    W1W3R2 = WF1(IR)*WF3(IR)
     &                                 *R2DRDI(IR,IM)
                                    TL(IR) = W1W3R2*RPW(IR,ILAM)
                                    TG(IR) = W1W3R2/RPW(IR,ILAM+1)
                                 END DO
C
                                 CALL RRADINT_R(IM,TL,SL)
                                 CALL RRADINT_R(IM,TG,SG)
C
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
                                 DO IR = 1,IRTOP
                                    TL(IR) = SL(IR)/RPW(IR,ILAM+1)
     &                                 + (SG(IRTOP)-SG(IR))*RPW(IR,ILAM)
                                 END DO
C
                                 DO IR = 1,IRTOP
                                    RINT(IR) = WF2(IR)*WF4(IR)*TL(IR)
     &                                 *R2DRDI(IR,IM)
                                 END DO
C
                                 CALL RRADINT(IM,RINT,APSINT)
C
                                 MD(IL2,INDIL(IL3,IL4),ILAM,IS13,IS24,
     &                              IE3,IE4)
     &                              = WIG_3J(DBLE(L1),DBLE(LAM),DBLE(L3)
     &                              ,0D0,0D0,0D0)
     &                              *WIG_3J(DBLE(L2),DBLE(LAM),DBLE(L4),
     &                              0D0,0D0,0D0)*APSINT
C
                                 SDD(IL3,IL4,IS13,IS24,IE3,IE4)
     &                              = SDD(IL3,IL4,IS13,IS24,IE3,IE4)
     &                              + (DBLE((2*L1+1)*(2*L2+1))
     &                              /DBLE(2*LAM+1))
     &                              *(MD(IL2,INDIL(IL3,IL4),ILAM,IS13,
     &                              IS24,IE3,IE4)**2)
C
                              END IF
C-------------------------------------- triangle l2,lam,l4
                           END DO
C
                        END DO
C
                        SDD(IL3,IL4,IS13,IS24,IE3,IE4)
     &                     = FOURPIESQ*SDD(IL3,IL4,IS13,IS24,IE3,IE4)
C
                     END DO
                  END DO
C
               END DO
C---------------------------------------------------- IE4
            END DO
C-------------------------------------------------------- IE3
C
         END DO
      END DO
C
C
C-------------------------------------------------------------------- MD
C
C
C-------------------------------------------------------------------- ME
C---------------------- compared with MD:  interchange indices  3  and 4
C
C
      DO IS14 = 1,2
         DO IS23 = 1,2
C
C-------------------------------------------------------- IE4
            DO IE4 = 1,NEME
C---------------------------------------------------- IE3
               DO IE3 = 1,NEME
C
                  IE2 = IE3 + IE4 - 1
C
                  DO IL3 = 1,NL
                     L3 = IL3 - 1
                     IREC = 2*((IE3-1)*NLLEED+(IL3-1)) + IS23
                     READ (IFIL34WF,REC=IREC) ITP,ILP,ISP,
     &                     (WF3(IR),IR=1,IRTOP)
C
                     DO IL4 = 1,NL
                        L4 = IL4 - 1
                        IREC = 2*((IE4-1)*NLLEED+(IL4-1)) + IS14
                        READ (IFIL34WF,REC=IREC) ITP,ILP,ISP,
     &                        (WF4(IR),IR=1,IRTOP)
C
                        DO IL2 = 1,NLLEED
                           L2 = IL2 - 1
                           IREC = 2*((IE2-1)*NLLEED+(IL2-1)) + IS23
                           READ (IFILLEED,REC=IREC) ITP,ILP,ISP,
     &                           (WF2(IR),IR=1,IRTOP)
C
                           DO LAM = ABS(L4-L1),(L4+L1),2
                              IF ( TRIANGLE(L2,LAM,L3) ) THEN
                                 ILAM = LAM + 1
C
                                 DO IR = 1,IRTOP
                                    W1W4R2 = WF1(IR)*WF4(IR)
     &                                 *R2DRDI(IR,IM)
                                    TL(IR) = W1W4R2*RPW(IR,ILAM)
                                    TG(IR) = W1W4R2/RPW(IR,ILAM+1)
                                 END DO
C
                                 CALL RRADINT_R(IM,TL,SL)
                                 CALL RRADINT_R(IM,TG,SG)
C
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
                                 DO IR = 1,IRTOP
                                    TL(IR) = SL(IR)/RPW(IR,ILAM+1)
     &                                 + (SG(IRTOP)-SG(IR))*RPW(IR,ILAM)
                                 END DO
C
                                 DO IR = 1,IRTOP
                                    RINT(IR) = WF2(IR)*WF3(IR)*TL(IR)
     &                                 *R2DRDI(IR,IM)
                                 END DO
C
                                 CALL RRADINT(IM,RINT,APSINT)
C
                                 ME(IL2,INDIL(IL4,IL3),ILAM,IS14,IS23,
     &                              IE4,IE3)
     &                              = WIG_3J(DBLE(L1),DBLE(LAM),DBLE(L4)
     &                              ,0D0,0D0,0D0)
     &                              *WIG_3J(DBLE(L2),DBLE(LAM),DBLE(L3),
     &                              0D0,0D0,0D0)*APSINT
C
                                 SEE(IL4,IL3,IS14,IS23,IE4,IE3)
     &                              = SEE(IL4,IL3,IS14,IS23,IE4,IE3)
     &                              + (DBLE((2*L1+1)*(2*L2+1))
     &                              /DBLE(2*LAM+1))
     &                              *(ME(IL2,INDIL(IL4,IL3),ILAM,IS14,
     &                              IS23,IE4,IE3)**2)
C
                              END IF
C-------------------------------------- triangle l2,lam,l3
                           END DO
C
                        END DO
C
                        SEE(IL4,IL3,IS14,IS23,IE4,IE3)
     &                     = FOURPIESQ*SEE(IL4,IL3,IS14,IS23,IE4,IE3)
C
                     END DO
                  END DO
C
               END DO
C---------------------------------------------------- IE3
            END DO
C-------------------------------------------------------- IE4
C
         END DO
      END DO
C
C
C-------------------------------------------------------------------- ME
C
C
C ======================================================================
C                      sum up matrix elements
C
C       SDE =  SUM(L2,LAM,LAM') 1/(2L3+1) ....  MD * ME
C
C
      L1 = LCXRAY(IT)
C
C------------------------------------------------------------------- SDE
      DO IS = 1,2
C
C-------------------------------------------------------- IE3
         DO IE3 = 1,NEME
C---------------------------------------------------- IE4
            DO IE4 = 1,NEME
               DO IL4 = 1,NL
                  L4 = IL4 - 1
C
                  DO IL3 = 1,NL
                     L3 = IL3 - 1
C
                     DO IL2 = 1,NLLEED
                        L2 = IL2 - 1
C
                        DO LAM = ABS(L3-L1),(L3+L1),2
                           IF ( TRIANGLE(L2,LAM,L4) ) THEN
                              ILAM = LAM + 1
C
                              DO LAMP = ABS(L4-L1),(L4+L1),2
                                 IF ( TRIANGLE(L2,LAMP,L3) ) THEN
                                    ILAMP = LAMP + 1
C
                                    MDXME = MD(IL2,INDIL(IL3,IL4),ILAM,
     &                                 IS,IS,IE3,IE4)
     &                                 *ME(IL2,INDIL(IL4,IL3),ILAMP,IS,
     &                                 IS,IE4,IE3)
C
C
                                    SDE(IL3,IL4,IS,IS,IE3,IE4)
     &                                 = SDE(IL3,IL4,IS,IS,IE3,IE4)
     &                                 + DBLE((2*L1+1)*(2*L2+1))
     &                                 *DBLE((-1)**(L3+LAMP+L1))
     &                                 *WIG_6J_RACAH(DBLE(L3),DBLE(LAMP)
     &                                 ,DBLE(L2),DBLE(L4),DBLE(LAM),
     &                                 DBLE(L1))*MDXME
C
                                 END IF
C-------------------------------------- triangle l2,lamp,l3
                              END DO
C
                           END IF
C-------------------------------------- triangle l2,lam,l4
                        END DO
                     END DO
C
                     SDE(IL3,IL4,IS,IS,IE3,IE4)
     &                  = FOURPIESQ*SDE(IL3,IL4,IS,IS,IE3,IE4)
C
                  END DO
               END DO
C
            END DO
C---------------------------------------------------- IE4
         END DO
C-------------------------------------------------------- IE3
C
      END DO
C
C
C
C ======================================================================
C  GNU OUTPUT D, E, DE and W
C
      IF ( OUTGNU ) THEN
         CALL RINIT(NLMAX*NLMAX*2*2,WLS)
         WRITE (6,*) '  GNUPLOT - input files will be created '
         DO IS14 = 1,2
            DO IS23 = 1,2
               DO IL4 = 1,NL
                  DO IL3 = 1,NL
C
                     FILNAMDEC = 'D'//TXT_L(IL3)//TXT_L(IL4)//TXTS(IS23)
     &                           //TXTS(IS14)
                     OPEN (51,FILE=FILNAMDEC)
                     FILNAMDEC = 'E'//TXT_L(IL3)//TXT_L(IL4)//TXTS(IS23)
     &                           //TXTS(IS14)
                     OPEN (52,FILE=FILNAMDEC)
                     FILNAMDEC = 'DE'//TXT_L(IL3)//TXT_L(IL4)
     &                           //TXTS(IS23)//TXTS(IS14)
                     OPEN (53,FILE=FILNAMDEC)
                     FILNAMDEC = 'W'//TXT_L(IL3)//TXT_L(IL4)//TXTS(IS23)
     &                           //TXTS(IS14)
                     OPEN (54,FILE=FILNAMDEC)
C
                     DO IE4 = 1,NEME
                        E4 = DREAL(EME(IE4))
C
                        DO IE3 = 1,NEME
                           E3 = DREAL(EME(IE3))
C
                           WRITE (51,99010) E3,E4,
     &                            SDD(IL3,IL4,IS23,IS14,IE3,IE4)
                           WRITE (52,99010) E3,E4,
     &                            SEE(IL3,IL4,IS23,IS14,IE3,IE4)
                           WRITE (53,99010) E3,E4,
     &                            SDE(IL3,IL4,IS23,IS14,IE3,IE4)
C
C ======================================================================
C                      Compute the cross sections W(L3,L4,IS3,IS4,E3,E4)
C
                           WLS(IL3,IL4,IS23,IS14)
     &                        = SDD(IL3,IL4,IS23,IS14,IE3,IE4)
     &                        + SEE(IL4,IL3,IS14,IS23,IE4,IE3)
C
                           IF ( IS23.EQ.IS14 ) WLS(IL3,IL4,IS23,IS23)
     &                          = WLS(IL3,IL4,IS23,IS23)
     &                          - 2D0*SDE(IL3,IL4,IS23,IS23,IE3,IE4)
                           IF ( WLS(IL3,IL4,IS23,IS14).LT.0D0 ) THEN
                              WRITE (6,'(A)')
     &                                ' NEGATIVE APS CROSS-SECTION ? '
                              WRITE (6,99011) E3,E4,IL3,IL4,IS23,IS14,
     &                               WLS(IL3,IL4,IS23,IS14)
                           END IF
                           WRITE (54,99010) E3,E4,WLS(IL3,IL4,IS23,IS14)
C
C
                        END DO
C
                        WRITE (51,99010)
                        WRITE (52,99010)
                        WRITE (53,99010)
                        WRITE (54,99010)
C
                     END DO
C
                  END DO
               END DO
            END DO
         END DO
C
         CLOSE (51)
         CLOSE (52)
         CLOSE (53)
         CLOSE (54)
C
      END IF
C
C  GNU OUTPUT D, E, DE and W
C ======================================================================
C
C
C ======================================================================
C                                                               read DOS
C
      REWIND 10
      WRITE (6,'(/,5X,A)') 'reading spin-polarized DOS from file 10'
C
      READ (10,99013) HEADER
      READ (10,99013) HEADER
      READ (10,99013) SYSTEM
      READ (10,99014) NQ
      READ (10,99014) NT
      READ (10,99014)
      READ (10,99014) NE
      READ (10,99014) IRELIN
      READ (10,99015) EFIN
      READ (10,99013) INFO
      READ (10,99013)
C
      WRITE (6,99012) HEADER,SYSTEM,NQ,NT,NSPIN,NE,IRELIN,EFIN,INFO
      READ (10,99013)
      DO IQP = 1,NQ
         READ (10,99013)
      END DO
C
      READ (10,99013)
      DO ITP = 1,NT
         READ (10,99013)
      END DO
C
      READ (10,99013)
C
      DO IE = 1,NE
C
         READ (10,99016) EDOS,
     &                   (((DOS(IE,IL,IS,ITP),IL=1,NL),IS=1,2),ITP=1,NT)
         ETAB(IE,1) = EDOS
      END DO
C
      DE = DREAL(ETAB(2,1)-ETAB(1,1))
C
      WRITE (6,99017) NE,DREAL(ETAB(1,1)),DREAL(ETAB(1,1))*RY_EV,
     &                DREAL(ETAB(NE,1)),DREAL(ETAB(NE,1))*RY_EV,DE,
     &                DE*RY_EV
C
C                                                           end read DOS
C ======================================================================
C
C
C ======================================================================
C                                   linear E-dependent broadening of DOS
      IF ( BRDDOS ) THEN
C
         WRITE (6,'(A)') ' The DOS will be broadened BEFORE convolution'
C
         DELE = DREAL(ETAB(2,1)-ETAB(1,1))
         NADD = NINT(ADDRANGE/DELE) + 1
C
         IF ( NADD+1.GT.101 ) STOP 'NADD+1 > 101 IN < NRAPS >'
C
C
C-----------------------------------------------------------------------
C              SET ENERGY-RANGE FOR BROADENING
C
C     1..NADD                5 eV BELOW E_F (DELE - DISTANCED)
C     NADD+1                 E_F - 0.001 RY (IT IS ASSUMED THAT
C                            ETAB(1,1) = EFERMI)
C     NADD+2..NE+NADD+1      CALCULATED DOS RANGE
C     NE+NADD+1..NE+2NADD+1  5 eV ABOVE CALCULATED LIMIT (DELE - DIST.)
C
C-----------------------------------------------------------------------
C
         DO IE = NADD,1, - 1
            EBRD(IE) = ETAB(1,1) - (NADD+1-IE)*DELE
         END DO
         EBRD(NADD+1) = ETAB(1,1) - 0.001D0
         DO IE = NADD + 2,NADD + NE + 1
            EBRD(IE) = ETAB(IE-NADD-1,1)
         END DO
         DO IE = 1,NADD
            EBRD(IE+NADD+NE+1) = ETAB(NE,1) + IE*DELE
         END DO
         ITP = IT
         DO IL = 1,NL
            DO IS = 1,2
C
C     VECBR contains the extended DOS entering in broadening routine
C
               DO J = 1,NADD + 1
                  VECBR(J) = 0D0
               END DO
               DO J = NADD + 2,NADD + NE + 1
                  VECBR(J) = DOS(J-NADD-1,IL,IS,ITP)
               END DO
               DO J = 1,NADD
                  VECBR(J+NADD+NE+1) = DOS(NE,IL,IS,ITP)
               END DO
C
C     perform linear E-dependent broadening of DOS
C
               NEBRD = 2*NADD + NE + 1
               DO I = 1,NEBRD
                  WLORTAB(I) = 0.0D0
                  IF ( DREAL(EBRD(I)).GE.EFIN ) WLORTAB(I)
     &                 = GAMMAV*DREAL(EBRD(I))
               END DO
               CALL VECLORBRD(EBRD,VECBR,NEBRD,NEMAX+101,WLORTAB)
C
C     overwrite broadened DOS from EBOT to EMAX
C
               DO IE = 1,NADD
                  DOS(IE,IL,IS,ITP) = VECBR(IE)
               END DO
               DO IE = NADD + 2,NE + NADD + 1
                  DOS(IE-1,IL,IS,ITP) = VECBR(IE)
               END DO
            END DO
         END DO
C
C     adjust E-range
C
         DO IE = 1,NADD
            ETAB(IE,1) = EBRD(IE)
         END DO
         DO IE = NADD + 2,NE + NADD + 1
            ETAB(IE-1,1) = EBRD(IE)
         END DO
         NE = NE + NADD
C
C
C     write broadened DOS
C
         FILNAMDEC = DATSET(1:LDATSET)//'.dos.brd'
         OPEN (12,FILE=FILNAMDEC)
         CALL WRHEAD(12,FILNAMDEC,'SRAPSDOS  ',NE)
         WRITE (12,'(A10,A,A)') 'DOS-FMT:  ','OLD-SPRKKR'
         DO IE = 1,NE
            WRITE (12,99016) ETAB(IE,1),
     &                       ((DOS(IE,IL,IS,IT),IL=1,NL),IS=1,2)
         END DO
         CLOSE (12)
C
      END IF
C
C                                   linear E-dependent broadening of DOS
C ======================================================================
C
C ======================================================================
C
C                      CALCULATE APS - SPECTRUM
C
C ======================================================================
C
C
      CALL RINIT((2*NEMAX-1)*(1+NLMAX)*(1+2),IAPS)
      IF ( CALCLS ) THEN
         CALL RINIT((2*NEMAX-1)*NLMAX*NLMAX*2*2,DECIAPS)
         CALL RINIT(NLMAX*NLMAX*2*2,WLS)
      END IF
C
      DO IE2 = 1,2*NE - 1
         E2 = 2.D0*DREAL(EBOT) - E1 + (IE2-1)*DE
         EA = DREAL(EBOT)
         EB = E2 + E1 - DREAL(EBOT)
C
         NE3 = NINT((EB-EA)/DE) + 1
         DO IE3 = 1,NE3
            IE4 = NE3 + 1 - IE3
            E3 = EA + (IE3-1)*DE
            E4 = EA + (IE4-1)*DE
C
            IF ( (IE3.EQ.1) .OR. (IE3.EQ.NE3) ) THEN
               WGTE3 = DE*TWOPIOVHBAR
            ELSE
               WGTE3 = 2D0*DE*TWOPIOVHBAR
            END IF
C
C----------------------- interpolate SDD, SEE, SDE from coarse ME E-mesh
C----------------------------------------------- to the finer DOS E-mesh
            IEME30 = NEME
            IEME40 = NEME
            DO IEME = 1,NEME
               REME = DREAL(EME(IEME))
               IF ( (IEME30.EQ.NEME) .AND. (REME.GT.E3) )
     &              IEME30 = MAX(2,IEME-1)
               IF ( (IEME40.EQ.NEME) .AND. (REME.GT.E4) )
     &              IEME40 = MAX(2,IEME-1)
            END DO
            IEME30 = MIN(IEME30,NEME-1)
            IEME40 = MIN(IEME40,NEME-1)
C
            P3 = (E3-DREAL(EME(IEME30)))/DEME
            Q4 = (E4-DREAL(EME(IEME40)))/DEME
C
            WA = Q4*(Q4-1D0)/2D0
            WB = P3*(P3-1D0)/2D0
            WC = 1D0 + P3*Q4 - P3**2 - Q4**2
            WD = P3*(P3-2D0*Q4+1D0)/2D0
            WE = Q4*(Q4-2D0*P3+1D0)/2D0
            WF = P3*Q4
C
            DO ISB = 1,2
               DO ISA = 1,2
                  DO ILB = 1,NL
                     DO ILA = 1,NL
C
                        TDD(ILA,ILB,ISA,ISB)
     &                     = SDD(ILA,ILB,ISA,ISB,IEME30+0,IEME40-1)
     &                     *WA + SDD(ILA,ILB,ISA,ISB,IEME30-1,IEME40+0)
     &                     *WB + SDD(ILA,ILB,ISA,ISB,IEME30+0,IEME40+0)
     &                     *WC + SDD(ILA,ILB,ISA,ISB,IEME30+1,IEME40+0)
     &                     *WD + SDD(ILA,ILB,ISA,ISB,IEME30+0,IEME40+1)
     &                     *WE + SDD(ILA,ILB,ISA,ISB,IEME30+1,IEME40+1)
     &                     *WF
C
                        TEE(ILA,ILB,ISA,ISB)
     &                     = SEE(ILA,ILB,ISA,ISB,IEME40-1,IEME30+0)
     &                     *WA + SEE(ILA,ILB,ISA,ISB,IEME40+0,IEME30-1)
     &                     *WB + SEE(ILA,ILB,ISA,ISB,IEME40+0,IEME30+0)
     &                     *WC + SEE(ILA,ILB,ISA,ISB,IEME40+0,IEME30+1)
     &                     *WD + SEE(ILA,ILB,ISA,ISB,IEME40+1,IEME30+0)
     &                     *WE + SEE(ILA,ILB,ISA,ISB,IEME40+1,IEME30+1)
     &                     *WF
C
                        TDE(ILA,ILB,ISA,ISB)
     &                     = SDE(ILA,ILB,ISA,ISB,IEME30+0,IEME40-1)
     &                     *WA + SDE(ILA,ILB,ISA,ISB,IEME30-1,IEME40+0)
     &                     *WB + SDE(ILA,ILB,ISA,ISB,IEME30+0,IEME40+0)
     &                     *WC + SDE(ILA,ILB,ISA,ISB,IEME30+1,IEME40+0)
     &                     *WD + SDE(ILA,ILB,ISA,ISB,IEME30+0,IEME40+1)
     &                     *WE + SDE(ILA,ILB,ISA,ISB,IEME30+1,IEME40+1)
     &                     *WF
C
                     END DO
                  END DO
               END DO
            END DO
C
C
            IF ( CALCLS ) THEN
               DO ILA = 1,NL
                  DO ILB = 1,NL
                     DO ISA = 1,2
                        DO ISB = 1,2
C
                           WLS(ILA,ILB,ISA,ISB) = TDD(ILA,ILB,ISA,ISB)
     &                        + TEE(ILB,ILA,ISB,ISA)
C
                           IF ( ISA.EQ.ISB ) WLS(ILA,ILB,ISA,ISA)
     &                          = WLS(ILA,ILB,ISA,ISA)
     &                          - 2D0*TDE(ILA,ILB,ISA,ISA)
                           IF ( ABS(WLS(ILA,ILB,ISA,ISB)).LT.1D-14 )
     &                          WLS(ILA,ILB,ISA,ISB) = 0D0
                        END DO
                     END DO
                  END DO
               END DO
            END IF
C
            DO IL = 1,NL
               DO IS2 = 1,2
                  IL4 = IL
                  DO IL3 = 1,NL
                     DO IS3 = 1,2
                        IAPS(IE2,IL,IS2) = IAPS(IE2,IL,IS2)
     &                     + WGTE3*DOS(IE4,IL4,IS2,IT)
     &                     *DOS(IE3,IL3,IS3,IT)*TDD(IL3,IL4,IS3,IS2)
C
                     END DO
                  END DO
C
                  DO IL3 = 1,NL
                     DO IS3 = 1,2
                        IAPS(IE2,IL,IS2) = IAPS(IE2,IL,IS2)
     &                     + WGTE3*DOS(IE4,IL4,IS3,IT)
     &                     *DOS(IE3,IL3,IS2,IT)*TEE(IL4,IL3,IS3,IS2)
C
                     END DO
                  END DO
C
                  DO IL3 = 1,NL
                     IAPS(IE2,IL,IS2) = IAPS(IE2,IL,IS2)
     &                                  - 2D0*WGTE3*DOS(IE4,IL4,IS2,IT)
     &                                  *DOS(IE3,IL3,IS2,IT)
     &                                  *TDE(IL3,IL4,IS2,IS2)
C
                  END DO
C
               END DO
            END DO
C
            IF ( CALCLS ) THEN
               DO IL = 1,NL
                  DO ILP = 1,NL
                     DO IS = 1,2
                        DO ISP = 1,2
                           DECIAPS(IE2,IL,ILP,IS,ISP)
     &                        = DECIAPS(IE2,IL,ILP,IS,ISP)
     &                        + WGTE3*DOS(IE3,IL,IS,IT)
     &                        *DOS(IE4,ILP,ISP,IT)*WLS(IL,ILP,IS,ISP)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
C
         END DO
      END DO
C
C ======================================================================
C                            fill up arrays
      DO IE2 = 1,2*NE - 1
         DO IL = 1,NL
            DO IS = 1,2
               IAPS(IE2,0,IS) = IAPS(IE2,0,IS) + IAPS(IE2,IL,IS)
            END DO
            IAPS(IE2,IL,0) = IAPS(IE2,IL,1) + IAPS(IE2,IL,2)
         END DO
         IAPS(IE2,0,0) = IAPS(IE2,0,1) + IAPS(IE2,0,2)
      END DO
C
C ======================================================================
C                            write spectrum
C the onset corresponds to E2-EFERMI = EFERMI-E1
C
C
      EONSET = 2.D0*EFERMI - E1
C
      DO IE = 1,2*NE - 1
         E2 = (2.D0*DREAL(EBOT)-E1+(IE-1)*DE) - EONSET
         WRITE (7,'(30E12.5)') E2*RY_EV,(IAPS(IE,0,IS),IS=0,2),
     &                         ((IAPS(IE,IL,IS),IL=1,NL),IS=0,2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT .AND. IE.LE.3 ) WRITE (IFILBUILDBOT,99020)
     &        ROUTINE(1:LEN_TRIM(ROUTINE)),IT,IE,E2*RY_EV,
     &        (IAPS(IE,0,IS),IS=0,2),((IAPS(IE,IL,IS),IL=1,NL),IS=0,2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      END DO
      CLOSE (7)
C
      IF ( CALCLS ) THEN
         WRITE (6,*) '  Output of (l,s)-decomposed spectra '
C
         DO IL3 = 1,NL
            DO IL4 = 1,NL
               FILNAMDEC = FILNAM(1:LFN+4)//'.'//TXT_L(IL3)//TXT_L(IL4)
               OPEN (UNIT=7,FILE=FILNAMDEC)
C
               CALL WRHEAD(7,FILNAMDEC,'NRAPSDEC  ',NE)
C
               NTXRSGRP = 1
               WRITE (7,99018) 'NTXRSGRP  ',NTXRSGRP
               WRITE (7,99019) 'SPECTRUM  ',SPEC(1:LSP)//'_'//TXT_L(IL3)
     &                         //TXT_L(IL4)
               WRITE (7,99018) 'IT        ',IT
               WRITE (7,99018) 'NCXRAY    ',NCXRAY(IT)
               WRITE (7,99018) 'LCXRAY    ',LCXRAY(IT)
               WRITE (7,99004) NCST,(NKPCOR(ICST),ICST=1,NCST)
               WRITE (7,99005)
               E1 = -1D+10
               DO JCST = 1,NCST
                  IF ( ECOR(JCST).GT.E1 ) THEN
                     E1 = ECOR(JCST)
                     ICST = JCST
                  END IF
               END DO
               WRITE (7,99006) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                         (2*MM05COR(ICST)+1),IKMCOR(ICST,1),
     &                         XNORM(1),ECOR(ICST),ECOR(ICST)*RY_EV,
     &                         SZCOR(ICST),IZERO(ICST)
               IF ( NKPCOR(ICST).EQ.2 ) WRITE (7,99007) IKMCOR(ICST,2),
     &              XNORM(2)
C
               DO IE = 1,2*NE - 1
                  E2 = (2.D0*DREAL(EBOT)-E1+(IE-1)*DE) - EONSET
                  WRITE (7,'(30(E12.5))') E2*RY_EV,
     &                   ((DECIAPS(IE,IL3,IL4,IS1,IS2),IS1=1,2),IS2=1,2)
               END DO
               CLOSE (7)
            END DO
         END DO
      END IF
C
      DEALLOCATE (EBRD,IAPS,CINT,WRCA,RINT,VECBR,MD,ME,SG,TG,SL,TL)
      DEALLOCATE (WF1,WF2,WF3,WF4,SDD,SDE,SEE,TDD,TDE,TEE,DOS,DECIAPS)
      DEALLOCATE (WLS,RPW,WLORTAB,BS,VS)
      DEALLOCATE (NLTLEED)
C
      STOP
C ======================================================================
99001 FORMAT (//,60('*'),/,10X,' STOP in <NRAPS>',/,10X,'LCXRAY=',I3,
     &        ' too large for NCSTMAX=',I3)
99002 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*       *    *  *****            ***    *****    ****        *'
     &  ,/,10X,
     &  '*       **   *  *    *          *   *   *    *  *    *       *'
     &  ,/,10X,
     &  '*       * *  *  *    *         *     *  *    *  *            *'
     &  ,/,10X,
     &  '*       *  * *  *****    ***   *******  *****    ****        *'
     &  ,/,10X,
     &  '*       *   **  *  *           *     *  *            *       *'
     &  ,/,10X,
     &  '*       *    *  *   *          *     *  *       *    *       *'
     &  ,/,10X,
     &  '*       *    *  *    *         *     *  *        ****        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
C
99003 FORMAT (/,4X,' CORE QUANTUM-NUMBERS  FOR  IT=',I2,':   N=',I2,
     &        '  L=',I2,/)
99004 FORMAT (//,' CORE STATES :',//,' NCST:  ',I4,/,' NKPCOR:',20I4)
99005 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99006 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99007 FORMAT (22X,I4,F12.6)
C
99008 FORMAT (' ME:   IS=',I2,' wave functions calculated for NEME=',I3,
     &        ' energies')
99009 FORMAT (' ME:   from   EMIN   =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '       to     EMAX   =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '       with   DEME   =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '              DEBOT1 =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '       and    DEF1   =',F8.3,' Ry    (',F8.3,' eV)')
C
99010 FORMAT (2F7.4,E12.5)
99011 FORMAT (' E3 = ',F8.3,' E4 = ',F8.3,' IL3 = ',I2,' IL4 = ',I2,
     &        ' IS3 = ',I2,' IS4 = ',I2,/,
     &        '              W(LS,L''S'') = ',F16.10,/,
     &        ' *************************************************')
99012 FORMAT (A,/,A,/,5X,'NQ=',I3,'  NT=',I3,'  NS=',I3,'  NE=',I3,/,5X,
     &        'IREL=',I3,'  EFERMI=',F7.4,/,A)
99013 FORMAT (10X,A80)
99014 FORMAT (10X,I10)
99015 FORMAT (10X,F10.5)
99016 FORMAT (8E10.4,:,/,(10X,7E10.4))
99017 FORMAT (' DOS:  number of energies read in  NE=',I3,/,
     &        '       from   E(1)  =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '       to     E(NE) =',F8.3,' Ry    (',F8.3,' eV)',/,
     &        '       with   DE    =',F8.3,' Ry    (',F8.3,' eV)')
99018 FORMAT (A10,I10)
99019 FORMAT (A10,A)
99020 FORMAT ('# BUILDBOT: ',A,': APS intensity for IT =',I5,/,'#',10X,
     &        'energy  IE =',I5,' E2 = ',2F10.6,/,(1PE22.14))
      END
C*==veclorbrd.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE VECLORBRD(E,F0,NE,NEMAX,WLORTAB)
C **********************************************************************
C *                                                                    *
C *  lorentzian-broadening of the function   F0   width: G0            *
C *                                                                    *
C **********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--VECLORBRD1301
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NE,NEMAX
      COMPLEX*16 E(NEMAX)
      REAL*8 F0(NEMAX),WLORTAB(NEMAX)
C
C Local variables
C
      REAL*8 FB(NEMAX),G,W,X0,X1,X2
      INTEGER IE,J
C
C*** End of declarations rewritten by SPAG
C
      DO IE = 1,NE
         FB(IE) = 0.0D0
      END DO
C
      DO IE = 1,NE
C
         X0 = DREAL(E(IE))
C
         DO J = 1,NE - 1
            X1 = DREAL(E(J)) + 0.000001D0
            X2 = DREAL(E(J+1)) - 0.000001D0
            G = (WLORTAB(J)+WLORTAB(J+1))/2.D0
            IF ( DABS(G).LT.1D-12 ) THEN
               W = 0.0D0
            ELSE
               W = ATAN((X2-X1)/(G+(X2-X0)*(X1-X0)/G))/2.D0
            END IF
            FB(IE) = FB(IE) + W*(F0(J)+F0(J+1))/PI
C
         END DO
C
      END DO
C
      DO IE = 1,NE
         F0(IE) = FB(IE)
      END DO
      END
