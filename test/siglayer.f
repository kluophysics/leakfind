C*==siglayer.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE SIGLAYER
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the conductance of layered systems      *
C   *                                                                  *
C   *  NOTE:                                                           *
C   *      *  the Green's function convention is   RH                  *
C   *      *  the conductance is calculated using the                  *
C   *         (L,ML,MS) representation - even for IREL = 3             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_THERMAL,ONLY:UMAT_VT
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,NCPA
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,NLMMAX,NSPIN,NXM,MEZJ,MEZZ,
     &    TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI,IQBOT_L,IQBOT_R,
     &    IQBOT_TB,NQTB,IQBOT,IQTOP
      USE MOD_TYPES,ONLY:NT,NTMAX,ITBOT,ITTOP
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,LCHECKME_SPIN,LJTOZ,
     &    NSPINPROJ
      USE MOD_CALCMODE,ONLY:ORBPOL,SIGMA_PROJECT,THERMAL_VIBRA_FLUCT,
     &    DMFT,LDAU,KKRMODE,IREL,GF_CONV_RH
      USE MOD_ENERGY,ONLY:ETAB,NEPATH,NETAB,EFERMI,NEMAX,EIMAG
      USE MOD_FILES,ONLY:IFILCBWF,RECLNGWF,IPRINT,WRTAUMQ,IFILTAU,
     &    TAUFIL,FOUND_SECTION
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,NPROCS,MPI_ID,MPI_ELOOP
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY
      USE MOD_CONSTANTS,ONLY:CI
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGLAYER')
C
C Local variables
C
      LOGICAL ALLIJ,CALCINT,CHECK,FOUND,GETIRRSOL,INTEG
      REAL*8 CPACHNG,CPACHNGMAX,DERYD,DMFTMIX,ERYDTOP,TIME,TIME0,TMP
      COMPLEX*16 DSSQX(:,:,:),ERYD,ERYD9,ERYDA,ERYDACC,ERYDB,ERYDBCC,
     &           GREF_I1(:,:,:),MAQAB(:,:,:,:,:),MAQBA(:,:,:,:,:),
     &           MBQAB(:,:,:,:,:),MBQBA(:,:,:,:,:),MCQAB(:,:,:,:,:),
     &           MCQBA(:,:,:,:,:),MDQAB(:,:,:,:,:),MDQBA(:,:,:,:,:),
     &           MEOFF_AHE_T(:,:,:,:,:),MEOFF_BARG_T(:,:,:,:,:),
     &           ME_SZ_T(:,:,:,:,:),MREG_TAB(:,:,:,:,:,:,:),
     &           MREG_TBA(:,:,:,:,:,:,:),MSSQA(:,:,:),MSSQB(:,:,:),
     &           MSSQX(:,:,:),MSSTA(:,:,:),MSSTB(:,:,:),P,PA,PB,
     &           S10AQAB(:,:,:,:),S10BQAB(:,:,:,:),S10CQAB(:,:,:,:),
     &           S10DQAB(:,:,:,:),S10TAB(:,:,:,:,:,:),S2AQAB(:,:,:,:),
     &           S2BQAB(:,:,:,:),S2CQAB(:,:,:,:),S2DQAB(:,:,:,:),
     &           S2TAB(:,:,:,:,:,:),S3AQAB(:,:,:,:),S3BQAB(:,:,:,:),
     &           S3CQAB(:,:,:,:),S3DQAB(:,:,:,:),S3TAB(:,:,:,:,:,:),
     &           S4AQAB(:,:,:,:),S4BQAB(:,:,:,:),S4CQAB(:,:,:,:),
     &           S4DQAB(:,:,:,:),S4TAB(:,:,:,:,:,:),TAUTA(:,:,:),
     &           TAUTB(:,:,:),TSSQA(:,:,:),TSSQB(:,:,:),TSSQX(:,:,:),
     &           TSSREF(:,:,:),TSSTA(:,:,:),TSSTB(:,:,:),WA(:,:),WX(:,:)
      INTEGER IA_ERR,ICPACONV,ICPAFLAG,IE,IEA,IEB,IEPATH,IESORT(:),
     &        IFILA,IFILACC,IFILB,IFILBCC,IJ,ILOOP_SIGMA_PROJECT,INC,
     &        IOFF,IPOL,IPROC,IPROCE(:),IQ,IQTB_IGIJ(:),IREF,ISPIN,
     &        ISPINPROJ,ISPR,ITCPA,ITMP,IWRI,IWRIRRWF,IWRREGWF,JE,JQ,
     &        JQTB_IGIJ(:),M,NCPAFAIL,NEA,NEB,NEINTEG,NEMAXSAV,NEMAXTMP,
     &        NETAU,NGIJ,NSP,NTMP,OUT
C
C*** End of declarations rewritten by SPAG
C
      DATA OUT/6/,CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/
      DATA ERYD9/(999999D0,999999D0)/
C
      ALLOCATABLE S10AQAB,S2AQAB,S3AQAB,S4AQAB
      ALLOCATABLE S10BQAB,S2BQAB,S3BQAB,S4BQAB
      ALLOCATABLE S10CQAB,S2CQAB,S3CQAB,S4CQAB
      ALLOCATABLE S10DQAB,S2DQAB,S3DQAB,S4DQAB
      ALLOCATABLE S10TAB,S2TAB,S3TAB,S4TAB
C
      ALLOCATABLE MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA
      ALLOCATABLE IPROCE,IESORT
      ALLOCATABLE MEOFF_AHE_T,MEOFF_BARG_T,ME_SZ_T
      ALLOCATABLE TSSTA,MSSTA,TSSQA,MSSQA,TAUTB
      ALLOCATABLE TSSTB,MSSTB,TSSQB,MSSQB,TAUTA
      ALLOCATABLE MREG_TAB
      ALLOCATABLE MREG_TBA
C
      ALLOCATABLE TSSREF,GREF_I1,TSSQX,DSSQX,MSSQX,WX,WA
      ALLOCATABLE IQTB_IGIJ,JQTB_IGIJ
C
      M = NKMMAX
C
      ALLOCATE (TSSQA(M,M,NQMAX),MSSQA(M,M,NQMAX))
      ALLOCATE (TSSQB(M,M,NQMAX),MSSQB(M,M,NQMAX))
C
      ALLOCATE (MAQAB(M,M,3,NSPINPROJ,NQMAX),
     &          MAQBA(M,M,3,NSPINPROJ,NQMAX))
      ALLOCATE (MBQAB(M,M,3,NSPINPROJ,NQMAX),
     &          MBQBA(M,M,3,NSPINPROJ,NQMAX))
      ALLOCATE (MCQAB(M,M,3,NSPINPROJ,NQMAX),
     &          MCQBA(M,M,3,NSPINPROJ,NQMAX))
      ALLOCATE (MDQAB(M,M,3,NSPINPROJ,NQMAX),
     &          MDQBA(M,M,3,NSPINPROJ,NQMAX))
C
      NSP = NSPINPROJ
C
      ALLOCATE (S10AQAB(3,3,NSP,NQMAX))
      ALLOCATE (S10BQAB(3,3,NSP,NQMAX))
      ALLOCATE (S10CQAB(3,3,NSP,NQMAX))
      ALLOCATE (S10DQAB(3,3,NSP,NQMAX))
      ALLOCATE (S2AQAB(3,3,NSP,NQMAX))
      ALLOCATE (S2BQAB(3,3,NSP,NQMAX))
      ALLOCATE (S2CQAB(3,3,NSP,NQMAX))
      ALLOCATE (S2DQAB(3,3,NSP,NQMAX))
      ALLOCATE (S3AQAB(3,3,NSP,NQMAX))
      ALLOCATE (S3BQAB(3,3,NSP,NQMAX))
      ALLOCATE (S3CQAB(3,3,NSP,NQMAX))
      ALLOCATE (S3DQAB(3,3,NSP,NQMAX))
      ALLOCATE (S4AQAB(3,3,NSP,NQMAX))
      ALLOCATE (S4BQAB(3,3,NSP,NQMAX))
      ALLOCATE (S4CQAB(3,3,NSP,NQMAX))
      ALLOCATE (S4DQAB(3,3,NSP,NQMAX),STAT=IA_ERR)
C
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc MAQAB')
C
C      CHECK = .TRUE.
      CHECK = .FALSE.
C
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' )
     &      CALL STOP_MESSAGE(ROUTINE,'KKRMOD = STANDARD')
      IF ( SYSTEM_DIMENSION(1:2).NE.'2D' )
     &      CALL STOP_MESSAGE(ROUTINE,'SYSTEM_DIMENSION <> 2D ')
      IF ( SYSTEM_TYPE(1:3).NE.'LIR' .AND. SYSTEM_TYPE(1:3).NE.'LIV' )
     &     CALL STOP_MESSAGE(ROUTINE,'SYSTEM_TYPE <> LIR or LIV')
      IF ( IREL.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'IREL <> 3 ')
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      IF ( NO_SYMMETRY ) WRITE (6,*) 'NO_SYMMETRY'
C
      GF_CONV_RH = .TRUE.
C      GF_CONV_RH = .false.
C
C ======================================================================
C          prepare the k-mesh without symmetry restrictions
C ======================================================================
C
      CALL INIT_MOD_TAUIJ_KMESH
C
C ----------------------------------------------------------------------
C
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TSSREF')
C
      ALLOCATE (WX(NXM,NXM))
      ALLOCATE (WA(NKMMAX,NKMMAX))
      ALLOCATE (MSSQX(NXM,NXM,NQ))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TSSREF')
C
      NKMSLAY = NSLAY_PER_PLAY*NXM
C
C=======================================================================
C                  set pairs of layers to be treated
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      CALL SECTION_FIND_KEYWORD('ALLIJ',ALLIJ)
C
      IF ( .NOT.ALLIJ ) THEN
C
         CALL SECTION_SET_INTEGER('NGIJ',NGIJ,9999,1)
         ALLOCATE (IQTB_IGIJ(NGIJ),JQTB_IGIJ(NGIJ))
C
         CALL SECTION_SET_INTEGER_ARRAY('IQ_IGIJ',IQTB_IGIJ,NGIJ,NGIJ,1,
     &                                  0,1)
         CALL SECTION_SET_INTEGER_ARRAY('JQ_IGIJ',JQTB_IGIJ,NGIJ,NGIJ,1,
     &                                  0,1)
C
      ELSE
C
         NGIJ = NQTB*NQTB
         ALLOCATE (IQTB_IGIJ(NGIJ),JQTB_IGIJ(NGIJ))
C
         IJ = 0
         DO IQ = IQBOT_TB,IQBOT_TB + NQTB - 1
            DO JQ = IQBOT_TB,IQBOT_TB + NQTB - 1
               IJ = IJ + 1
               IQTB_IGIJ(IJ) = IQ
               JQTB_IGIJ(IJ) = JQ
            END DO
         END DO
C
      END IF
C
      IQTB_IGIJ(1:NGIJ) = IQTB_IGIJ(1:NGIJ) - IQBOT_TB + 1
      JQTB_IGIJ(1:NGIJ) = JQTB_IGIJ(1:NGIJ) - IQBOT_TB + 1
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( DMFT .OR. LDAU ) CALL INIT_MOD_DMFT_LDAU(NETAB(1),ETAB(1,1),
     &     DMFTMIX,NEMAX)
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      TMP = 0D0
      NTMP = 1
C
      CALL THERMAL_INIT(0,NTMP,TMP)
C
C=======================================================================
C
      M = NKMMAX
C
      ALLOCATE (MSSTA(M,M,NTMAX),MSSTB(M,M,NTMAX))
      ALLOCATE (TAUTA(M,M,NTMAX),TAUTB(M,M,NTMAX))
      ALLOCATE (TSSTA(M,M,NTMAX),TSSTB(M,M,NTMAX))
      ALLOCATE (MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2))
      ALLOCATE (MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2))
      ALLOCATE (MEOFF_AHE_T(NKMMAX,NKMMAX,3,3,NTMAX))
      ALLOCATE (MEOFF_BARG_T(NKMMAX,NKMMAX,3,3,NTMAX))
      ALLOCATE (ME_SZ_T(NKMMAX,NKMMAX,3,3,NTMAX),STAT=IA_ERR)
C
C=======================================================================
      WRITE (6,99001)
C ======================================================================
C
      ALLOCATE (LIST_ISPR(NSPINPROJ))
C
      NSPR = NSPINPROJ
      LIST_ISPR(1:11) = (/1,2,3,4,5,6,7,8,9,10,11/)
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('SIGSPROJ',FOUND)
         IF ( FOUND ) CALL SECTION_SET_INTEGER_ARRAY('SIGSPROJ',
     &        LIST_ISPR,NSPINPROJ,NSPINPROJ,0,9999,0)
         CALL SECTION_FIND_KEYWORD('CHECKME',FOUND)
         IF ( FOUND ) THEN
            WRITE (6,'(10X,"<SIG_SPIN> checking MEs")')
            LCHECKME_SPIN = .TRUE.
            NSPR = 1
            LIST_ISPR(1:11) = (/1,0,0,0,0,0,0,0,0,0,0/)
            CALL SECTION_FIND_KEYWORD('JTOZ',FOUND)
            IF ( FOUND ) THEN
               LJTOZ = .TRUE.
               WRITE (6,'(10X,"<SIG_SPIN> will overwrite Js with Zs")')
            END IF
         END IF
      END IF
C
      WRITE (6,'(10x,"<SIG_SPIN> number of choosen operators: ",i5,/)')
     &       NSPR
C
      DO ISPR = 1,NSPR
         WRITE (6,'(10X,"<SIG_SPIN> operators: ",I3,1X,A80)') ISPR,
     &          STR_ISP_PROJ(LIST_ISPR(ISPR))
      END DO
C
C ======================================================================
C
      CALL CPU_TIME(TIME0)
C
      IF ( IBZINT.NE.2 .AND. IBZINT.NE.6 .AND. IBZINT.NE.4 ) THEN
         WRITE (6,'("<SIG_SPIN>:  IBZINT   <>   2 | 4 | 6")')
         STOP
      END IF
C
      IF ( IBZINT.EQ.4 ) THEN
C
         ERYDTOP = 1.4D0*EFERMI
C
         CALL FRELPOLE(ERYDTOP,IPRINT)
      END IF
C
      NCPAFAIL = 0
      ITCPAMAX = 30
C
      CALL READTAU(9,ERYD9,0,NETAU,TAUT,0,NT,.FALSE.,TAUQ(1,1,1),
     &             MSSQ(1,1,1),1,NQ,0,NKM,NKMMAX,IPRINT)
C
      DERYD = 0.001D0
C
      IF ( NETAB(1).GT.1 ) THEN
         NEINTEG = NETAB(1)
      ELSE
         NEINTEG = 0
      END IF
      NEB = NEINTEG + 1
      NEA = 0
C
      IF ( .NOT.(MPI) ) THEN
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .FALSE.
      ELSE IF ( NEB.EQ.1 ) THEN
         MPI_KLOOP = .TRUE.
         MPI_ELOOP = .FALSE.
      ELSE
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .TRUE.
      END IF
C
      NEMAXTMP = NEB*3 + 1
      NEMAXSAV = NEMAX
C
      IFILA = IFILCBWF
      IFILB = IFILCBWF + 1
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%    TAUQ       MSSQ                                   BLOCK BEGIN
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=======================================================================
C   calculate first TAUQ and MSSQ and store them for E-loops B and A
C=======================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      OPEN (IFILB,STATUS='SCRATCH',FORM='UNFORMATTED',ACCESS='DIRECT',
     &      RECL=RECLNGWF)
C
      IWRI = 100
C
      IF ( NETAU.EQ.0 .AND. NTMP.EQ.1 ) THEN
C
         CALL THERMAL_INIT(1,NTMP,TMP)
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop                               START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
         NEPATH = 1
         IEPATH = 1
         CPACHNGMAX = 0D0
C
         IE = 0
         DO IEB = 1,NEB
C
            IF ( IEB.LT.NEB ) THEN
               INTEG = .TRUE.
               ERYDB = ETAB(IEB,IEPATH)
               NEA = 2
C               if (LCHECKME_SPIN) erydb = dreal(erydb)
            ELSE
               INTEG = .FALSE.
               ERYDB = EFERMI + CI*EIMAG
               NEA = 0
            END IF
C
C
            DO IEA = 0,NEA
               IE = IE + 1
C
               IF ( INTEG ) THEN
                  ERYDA = ERYDB + (IEA-1.5D0)*DERYD
               ELSE
                  ERYDA = EFERMI + CI*EIMAG
               END IF
C
               IF ( IEA.EQ.0 ) THEN
                  ERYD = ERYDB
                  ERYDA = 0D0
               ELSE
                  ERYD = ERYDA
               END IF
C
               WRITE (6,99005) EFERMI,IEA,ERYDA,IEB,ERYDB,IE,ERYD
               ICPAFLAG = 0
               CPACHNG = 0.0D0
C
C ===================================== solve SS - differential equation
C
               IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &              = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
Csbsbsbsbsbsbsbsbsbsbsbsbsbbsbsbsbsbsbsb
               ITBOT = 1
               ITTOP = NT
               IQBOT = 1
               IQTOP = NQ
Csbsbsbsbsbsbsbsbsbsbsbsbsbbsbsbsbsbsbsb
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,
     &                       GETIRRSOL,ERYD,P,IPRINT,TSST,MSST,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
               CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
               NEMAX = NEMAXTMP
C
C ------------- for tau-matrix calculation,  here we make it temporarily
C ------------- parallel in k (even if we are parallel in E in sigma
C ------------  calculation)
               IF ( MPI_ELOOP ) MPI_KLOOP = .TRUE.
C
C***********************************************************************
C
C        no vibrations:
C        CONC   == X_VFT
C        NOQ    == NVO_Q
C        ITOQ   == IVFT_VFOQ
C        NTMAX  == NVFTMAX
C
               CALL TAU_TB(ERYD,P,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,
     &                     MSST,TSSQ,MSSQ,TAUQ)
C***********************************************************************
C
               IF ( MPI_ELOOP ) MPI_KLOOP = .FALSE.
C
               NEMAX = NEMAXSAV
C
               IF ( ICPAFLAG.NE.0 ) THEN
                  NCPAFAIL = NCPAFAIL + 1
                  CPACHNGMAX = MAX(CPACHNGMAX,ABS(CPACHNG))
               END IF
C
C------------------------ to avoid in asynchronisations using NFS sytems
C------------------------ we use this somewhat awkward scheme
               IF ( MPI ) THEN
                  CALL DRV_MPI_BARRIER
                  IF ( MPI_ID.EQ.0 ) THEN
                     CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                            TAUT)
                     CALL FLUSH(IFILTAU)
                  ELSE
                     CLOSE (IFILTAU)
                  END IF
                  CALL DRV_MPI_BARRIER
                  IF ( MPI_ID.NE.0 ) OPEN (UNIT=IFILTAU,FILE=TAUFIL)
                  CALL DRV_MPI_BARRIER
               ELSE
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
               END IF
C
C ======================================================================
C
               IF ( IPRINT.GE.3 .OR. CHECK )
     &              CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,TAUQ)
C
C-----------------------------------------------------------------------
               IF ( IE.EQ.NETAB(IEPATH) ) THEN
C
                  IF ( NCPAFAIL.NE.0 ) THEN
                     WRITE (OUT,99002) NCPAFAIL,CPATOL,CPACHNGMAX
                  ELSE IF ( NCPA.NE.0 ) THEN
                     WRITE (OUT,99003)
                  END IF
C
               END IF
C
C ======================================================================
C
               IF ( MPI ) CALL DRV_MPI_BARRIER
C
            END DO ! NEA
         END DO ! NEB
C
         CALL READTAU(9,ERYD9,0,NETAU,TAUT,0,NT,.FALSE.,TAUQ(1,1,1),
     &                MSSQ(1,1,1),1,NQ,0,NKM,NKMMAX,IPRINT)
      END IF
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop                                 END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%    TAUQ       MSSQ                                   BLOCK END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C              loop to get the SPIN-projected conductivities       START
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      ILOOP_SIGMA_PROJECT = 0
C
C
C
      IF ( SIGMA_PROJECT ) ILOOP_SIGMA_PROJECT = ILOOP_SIGMA_PROJECT + 1
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C***********************************************************************
C***********************************************************************
C****************   START OF SIGMA CALCULATION   ***********************
C***********************************************************************
C***********************************************************************
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      E-Punkte IEB = 1...,(NEB-1) fuer B-Schleife
C       auf Prozessoren verteilen wie in SCF
C      IEB=NEB auf IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCE(NEB),IESORT(NEB))
      DO IE = 1,NEB
         IESORT(IE) = IE
         IPROCE(IE) = 0
      END DO
C
      IF ( MPI_ELOOP ) THEN
         IPROC = 0
         INC = 1
         DO JE = 1,NEB - 1
            IE = IESORT(JE)
            IPROC = IPROC + INC
            IF ( IPROC.EQ.NPROCS ) THEN
               IPROC = NPROCS - 1
               INC = -1
            ELSE IF ( IPROC.EQ.-1 ) THEN
               IPROC = 0
               INC = 1
            END IF
            IPROCE(IE) = IPROC
         END DO
C
         CALL DRV_MPI_BARRIER
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C=======================================================================
C=======================================================================
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== START
C
CDK ToDo for Fermi sea calculation one needs to properly take care of
C        sequential run -- vibrations haveto be taken into account
C        --------------- sanity check -----------
      IF ( THERMAL_VIBRA_FLUCT .AND. MPI_ELOOP ) THEN
         WRITE (6,99006)
         STOP
      END IF
C
      DO ITMP = 1,NTMP
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop 2   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         WRTAUMQ = .FALSE.
C
         NEPATH = 1
         IEPATH = 1
C
         ICPAFLAG = 0
         CPACHNG = 0.0D0
C
C ======================================================================
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C            and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ
C
         CALL THERMAL_INIT(ITMP,NTMP,TMP)
C
C
         IF ( MPI ) CALL DRV_MPI_BARRIER
C
         IE = 0
C_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IE
         DO IEB = 1,NEB
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ID.EQ.IPROCE(IEB) .OR. MPI_KLOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C               ERYD = EFERMI + CI*EIMAG
C               WRITE (6,99001) EFERMI,IEB,ERYD
C
C
               IE = IE + 1
C
               IF ( IEB.LT.NEB ) THEN
                  INTEG = .TRUE.
                  ERYDB = ETAB(IEB,IEPATH)
                  NEA = 2
C                  if (LCHECKME_SPIN) erydb = dreal(erydb)
               ELSE
                  INTEG = .FALSE.
                  ERYDB = EFERMI + CI*EIMAG
                  NEA = 1
               END IF
C
C ===================================== solve SS - differential equation
C
               IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &              = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILB,GETIRRSOL,
     &                       ERYDB,PB,IPRINT,TSSTB,MSSTB,SSST,MEZZ,MEZJ,
     &                       ORBPOL)
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
               IF ( THERMAL_VIBRA_FLUCT ) CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
               CALL MSSINIT(TSSTB,MSSTB,TSSQB,MSSQB)
C
               IF ( MPI_ELOOP .OR. THERMAL_VIBRA_FLUCT ) THEN
C
                  NEMAX = NEMAXTMP
C
                  CALL TAU_TB(ERYDB,PB,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,
     &                        TSSTB,MSSTB,TSSQB,MSSQB,TAUQ)
C
                  NEMAX = NEMAXSAV
C
               ELSE
C
                  DO IQ = IQBOT_CHI,IQTOP_CHI
                     CALL READTAU(9,ERYD9,IE,NETAU,TAUTB,0,NT,.FALSE.,
     &                            TAUQ(1,1,IQ),MSSQB(1,1,IQ),IQ,NQ,1,
     &                            NKM,NKMMAX,IPRINT)
                     IF ( ABS(ERYDB-ERYD9).GT.1D-8 ) STOP 
     &                    'in <SIG_SPIN>: >>> FILE 9 ERYDB <> ERYD9'
                  END DO
C
               END IF
C
C ----------------------------------------------- project TAU on type IT
C
               CALL PROJTAU(ICPAFLAG,CPACHNG,ERYDB,MSSTB,MSSQB,TAUQ,
     &                      TAUTB)
C
C_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IE
               DO IEA = 1,NEA
C
                  IF ( INTEG ) THEN
                     IE = IE + 1
                     ERYDA = ERYDB + (IEA-1.5D0)*DERYD
                  ELSE
                     ERYDA = EFERMI + CI*EIMAG
                  END IF
C
                  WRITE (6,99005) EFERMI,IEA,ERYDA,IEB,ERYDB,IE
C
C ======================================================================
C ===================================== solve SS - differential equation
C
                  IF ( DMFT .OR. LDAU )
     &                 DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &                 = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
                  CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILA,
     &                          GETIRRSOL,ERYDA,PA,IPRINT,TSSTA,MSSTA,
     &                          SSST,MEZZ,MEZJ,ORBPOL)
C
                  CALL MSSINIT(TSSTA,MSSTA,TSSQA,MSSQA)
C
                  IF ( .NOT.(INTEG) ) THEN
C
                     DO IQ = IQBOT_CHI,IQTOP_CHI
C                       TAUQ(1,1,IQ)
                        MSSQA(:,:,IQ) = MSSQB(:,:,IQ)
                     END DO
C
                  ELSE IF ( MPI_ELOOP ) THEN
C
                     NEMAX = NEMAXTMP
C
                     CALL TAU_TB(ERYDA,PA,ICPAFLAG,CPACHNG,ITCPA,
     &                           ICPACONV,TSSTA,MSSTA,TSSQA,MSSQA,TAUQ)
C
                     NEMAX = NEMAXSAV
C
                  ELSE
C
                     DO IQ = IQBOT_CHI,IQTOP_CHI
                        CALL READTAU(9,ERYD9,IE,NETAU,TAUTA,0,NT,
     &                               .FALSE.,TAUQ(1,1,IQ),MSSQA(1,1,IQ),
     &                               IQ,NQ,1,NKM,NKMMAX,IPRINT)
                        IF ( ABS(ERYDA-ERYD9).GT.1D-8 ) STOP 
     &                       'in <SIG>: >>> FILE 9    ERYDA <> ERYD9'
                     END DO
C
                  END IF
C
C
C ----------------------------------------------- project TAU on type IT
C
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYDA,MSSTA,MSSQA,TAUQ,
     &                         TAUTA)
C
C***********************************************************************
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
                  IF ( .NOT.MPI_KLOOP .OR. MPI_ID.EQ.0 ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C ======================================================================
C
                     IF ( IPRINT.GE.3 .OR. CHECK )
     &                    CALL DUMPTAU(IE,ERYDA,IWRI,MSST,MSSQA,TAUTA,
     &                    TAUQ)
C
C ======================================================================
C
                     ERYDACC = DCONJG(ERYDA)
                     ERYDBCC = DCONJG(ERYDB)
                     IFILACC = IFILA
                     IFILBCC = IFILB
C
                     CALL SIGME(ILOOP_SIGMA_PROJECT,ERYDA,ERYDACC,ERYDB,
     &                          ERYDBCC,IFILA,IFILACC,IFILB,IFILBCC,
     &                          MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,
     &                          MDQAB,MDQBA,MREG_TAB,MREG_TBA,S10TAB,
     &                          S2TAB,S3TAB,S4TAB,S10AQAB,S10BQAB,
     &                          S10CQAB,S10DQAB,S2AQAB,S2BQAB,S2CQAB,
     &                          S2DQAB,S3AQAB,S3BQAB,S3CQAB,S3DQAB,
     &                          S4AQAB,S4BQAB,S4CQAB,S4DQAB,UMAT_VT,
     &                          MSST,TSST,MSSQA,TAUQ,UMAT_VT,MSST,TSST,
     &                          MSSQB,TAUQ,NSPINPROJ,MEOFF_BARG_T,
     &                          MEOFF_AHE_T,ME_SZ_T)
C
C **********************************************************************
C    supply ssite t-matrix  TSSREF of the  NREF  TB reference atoms
C **********************************************************************
C         supply the real space TB Green''s function matrices   GREF_I1
C            for the inequivalent   NCLU_REF  reference clusters
C **********************************************************************
C
                     CALL TBGREFRS(ERYD,P,TSSREF,GREF_I1)
C
C **********************************************************************
                     DO ISPIN = 1,NSPIN
C
C ----------------------------------------------------------------------
                        IF ( IREL.LE.2 ) THEN
C
                           IOFF = NLM*(ISPIN-1)
C
                           DO IQ = 1,NQ
                              IREF = IREFQ(IQ)
C
                              MSSQX(1:NLM,1:NLM,IQ)
     &                           = MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,
     &                           IQ)
C
                              WX(1:NLM,1:NLM) = MSSQX(1:NLM,1:NLM,IQ)
C
                              CALL CMATINV(NLM,NLM,WX,TSSQX(1,1,IQ))
C
                              WX(1:NLM,1:NLM) = TSSQX(1:NLM,1:NLM,IQ)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
C
                              CALL CMATINV(NLM,NLM,WX,DSSQX(1,1,IQ))
C
                           END DO
C
C ----------------------------------------------------------------------
                        ELSE IF ( IREL.EQ.3 ) THEN
C
Csbswsbswsbsw     complex conjugation + scaling of matrixelements
Csbswsbswsbsw     to be in line with Voicu''s TBKKR routines
C
                           MAQAB = DCONJG(MAQAB)
                           MAQAB = MAQAB*CI*SQRT(2.0D0)*2.0D0
C
                           DO IQ = 1,NQ
                              IREF = IREFQ(IQ)
C
                              CALL CHANGEREP(NKM,NKMMAX,MSSQ(1,1,IQ),
     &                           'REL>RLM',WA)
C
                              MSSQX(1:NKM,1:NKM,IQ) = WA(1:NKM,1:NKM)
C
                              CALL CMATINV(NKM,NKMMAX,WA,TSSQX(1,1,IQ))
C
                              WA(1:NKM,1:NKM) = TSSQX(1:NKM,1:NKM,IQ)
                              WA(1:NLM,1:NLM) = WA(1:NLM,1:NLM)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
                              WA(NLM+1:NKM,NLM+1:NKM)
     &                           = WA(NLM+1:NKM,NLM+1:NKM)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
C
                              CALL CMATINV(NKM,NKMMAX,WA,DSSQX(1,1,IQ))
C
                              DO IPOL = 1,3
                                 DO ISPINPROJ = 1,NSPINPROJ
C
                                    CALL CHANGEREP(NKM,NKMMAX,
     &                                 MAQAB(1,1,IPOL,ISPINPROJ,IQ),
     &                                 'REL>RLM',WA)
C
                                    MAQAB(1:NKM,1:NKM,IPOL,ISPINPROJ,IQ)
     &                                 = WA(1:NKM,1:NKM)
C
                                 END DO
                              END DO
C
                           END DO
C
C ----------------------------------------------------------------------
                        ELSE
                           CALL STOP_MESSAGE(ROUTINE,'IREL ??????')
                        END IF
C
C ----------------------------------------------------------------------
C
                        CALL SIGLAYKLOOP(DSSQX(1,1,IQBOT_TB),
     &                     DSSQX(1,1,IQBOT_L),DSSQX(1,1,IQBOT_R),
     &                     GREF_I1,ICLU_REF_IQ(IQBOT_TB),ISPIN,NSPIN,
     &                     MAQAB,NGIJ,IQTB_IGIJ,JQTB_IGIJ)
C
                     END DO
C **********************************************************************
C
C=======================================================================
                  END IF ! MPI
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
               END DO ! IEA
C_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IEA_IE
            END IF ! MPI
         END DO ! IEB
C_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IEB_IE
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.0 ) THEN
C
            CALL CPU_TIME(TIME)
            WRITE (6,99004) TIME - TIME0
         END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO ! ITMP
C
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== END =
C=======================================================================
C=======================================================================
C
      DEALLOCATE (MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA)
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
C   ====================================================================
C
99001 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*            ****   ***   ****   *     *    ***              *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  **   *    *   *             *'
     &  ,/,10X,
     &  '*           *        *   *       * * * *  *     *            *'
     &  ,/,10X,
     &  '*            ****    *   *  ***  *  *  *  *******            *'
     &  ,/,10X,
     &  '*                *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*            ****   ***   ****   *     *  *     *            *'
     &  ,/,10X,
     &  '*         *         ***    *     *  ******  ******           *'
     &  ,/,10X,
     &  '*         *        *   *    *   *   *       *     *          *'
     &  ,/,10X,
     &  '*         *       *     *    * *    *       *     *          *'
     &  ,/,10X,
     &  '*         *       *******     *     ****    ******           *'
     &  ,/,10X,
     &  '*         *       *     *     *     *       *   *            *'
     &  ,/,10X,
     &  '*         *       *     *     *     *       *    *           *'
     &  ,/,10X,
     &  '*         ******  *     *     *     ******  *     *          *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,1X,79('*'),/,10X,'CPA not converged for',I3,
     &        ' energies:',/,10X,'tolerance for CPA-cycle:',F15.7,/,10X,
     &        'maximum deviation:      ',F15.7,/,1X,79('*'),/)
99003 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99004 FORMAT (/,5X,'execution time for <SIG_SPIN>:',F14.3,' secs',/)
99005 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,25X,'EFERMI =',F10.5,/,10X,'IEA  =',I4,5X,
     &        'ERYDA  =',2F10.5,/,10X,'IEB  =',I4,5X,'ERYDB  =',2F10.5/,
     &        10X,'IE   =',I4,5X,:,'ERYD   =',2F10.5)
99006 FORMAT ('<SIG_SPIN> sequential mode and finite temp not yet ',
     &        'available')
      END
C*==siglaykloop.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE SIGLAYKLOOP(DSSQX,DSSQX_L,DSSQX_R,GREF_I1,
     &                       ICLU_REF_IQTB,ISPIN,NSPIN,MEQ,NGIJ,
     &                       IQTB_IGIJ,JQTB_IGIJ)
C   ********************************************************************
C   *                                                                  *
C   *  perform the loop over the k-mesh and calculate the current      *
C   *                                                                  *
C   *  NOTE:                                                           *
C   *      *  the Green's function convention is   RH                  *
C   *      *  the conductance is calculated using the                  *
C   *         (L,ML,MS) representation - even for IREL = 3             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_SIG,ONLY:NSPINPROJ
      USE MOD_SITES,ONLY:NQMAX,NQ_L,NQ_R,NQTB
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_ANGMOM,ONLY:NKMMAX,NLMMAX,NKM,NLM,NXM
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,KTABTAUIJ,WKTABTAUIJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGLAYKLOOP')
      LOGICAL USE_FULL_SYMMETRY
      PARAMETER (USE_FULL_SYMMETRY=.TRUE.)
C
C Dummy arguments
C
      INTEGER ISPIN,NGIJ,NSPIN
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_L(NXM,NXM,NQ_L),
     &           DSSQX_R(NXM,NXM,NQ_R),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           MEQ(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX)
      INTEGER ICLU_REF_IQTB(NQTB),IQTB_IGIJ(NGIJ),JQTB_IGIJ(NGIJ)
C
C Local variables
C
      COMPLEX*16 GLLIJ(:,:),GLLIJ_CONJG(:,:),GLLKE(:,:),GREFKE(:,:)
      INTEGER I,IA_ERR,IK,IKTOP,IPROCK(:),IQTB,IU,IU0,IV,IV0,IX,IX0,J,
     &        JQTB,JU,JU0,JV,JV0,JX
      REAL*8 IM_KVEC(3),KVEC(3),KVEC_K(:,:),WK,WKSUM,WK_K(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE IPROCK,GLLKE,GREFKE,GLLIJ,GLLIJ_CONJG,KVEC_K,WK_K
C
      ALLOCATE (GLLKE(NKKR_TB,NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate GLLKE')
      IF ( IREL.GT.2 ) THEN
         ALLOCATE (GREFKE(NKKRNR_TB,NKKRNR_TB),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate GREFKE')
      ELSE IF ( NKKR_TB.NE.NKKRNR_TB ) THEN
         CALL STOP_MESSAGE(ROUTINE,'NKKR_TB <> NKKRNR_TB')
      END IF
C
      ALLOCATE (GLLIJ_CONJG(NXM,NXM),GLLIJ(NXM,NXM),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate GLLIJ')
C
      IF ( USE_FULL_SYMMETRY ) THEN
         IKTOP = NKTAB
         WKSUM = SUM(WKTABTAUIJ(1:NKTABTAUIJ))
         ALLOCATE (KVEC_K(3,IKTOP),WK_K(IKTOP))
         KVEC_K(1:3,1:IKTOP) = KTAB(1:3,1:IKTOP)
         WK_K(1:IKTOP) = WKTABTAUIJ(1:IKTOP)
      ELSE
         IKTOP = NKTABTAUIJ
         WKSUM = SUM(WKTAB(1:NKTAB))
         ALLOCATE (KVEC_K(3,IKTOP),WK_K(IKTOP))
         KVEC_K(1:3,1:IKTOP) = KTABTAUIJ(1:3,1:IKTOP)
         WK_K(1:IKTOP) = WKTAB(1:IKTOP)
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., IKTOP
C      over processors;   IK=IKTOP  dealt with by IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCK(IKTOP))
      CALL MPI_DISTRIBUTE(IPROCK,IKTOP,MPI_KLOOP,'K')
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C                                                          K-points loop
      DO IK = 1,IKTOP
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCK(IK) .OR. MPI_ELOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IF ( USE_FULL_SYMMETRY ) THEN
               WK = WKTAB(IK)
               KVEC(1:3) = KTAB(1:3,IK)
            ELSE
               WK = WKTABTAUIJ(IK)
               KVEC(1:3) = KTABTAUIJ(1:3,IK)
            END IF
C
C-----------------------------------------------------------------------
C    Fourier transformation, set KKR matrix M = [-(t)^-1 + G^r]
C-----------------------------------------------------------------------
C
            IF ( IREL.LE.2 ) THEN
               CALL TBGREFKE(KVEC,IM_KVEC,NQTB,GREF_I1,GLLKE,
     &                       ICLU_REF_IQTB,NKKRNR_RSMAX)
C
            ELSE
C
               CALL TBGREFKE(KVEC,IM_KVEC,NQTB,GREF_I1,GREFKE,
     &                       ICLU_REF_IQTB,NKKRNR_RSMAX)
C
               CALL CINIT(NKKR_TB*NKKR_TB,GLLKE)
C
               JU0 = -NKM
               JV0 = -NLM
               DO JQTB = 1,NQTB
                  JU0 = JU0 + NKM
                  JV0 = JV0 + NLM
                  DO J = 1,NLM
                     JU = JU0 + J
                     JV = JV0 + J
C
                     IU0 = -NKM
                     IV0 = -NLM
                     DO IQTB = 1,NQTB
                        IU0 = IU0 + NKM
                        IV0 = IV0 + NLM
                        DO I = 1,NLM
                           IU = IU0 + I
                           IV = IV0 + I
                           GLLKE(IU,JU) = GREFKE(IV,JV)
                           GLLKE(NLM+IU,NLM+JU) = GREFKE(IV,JV)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
C
C-----------------------------------------------------------------------
C              call decimation routine if requested
C-----------------------------------------------------------------------
C
            IF ( IDECI.EQ.1 ) CALL TBDECIMATE(GLLKE,DSSQX_L,DSSQX_R,
     &           VACFLAG,FACTL,NQ_L,NQ_R,NXM)
C
C-----------------------------------------------------------------------
C    Construct the matrix M=[-(t-t_ref)^-1 + G_ref] and store it
C    in the same matrix GLLKE where  G_ref  was stored.
C-----------------------------------------------------------------------
C
            IX0 = -NXM
            DO IQTB = 1,NQTB
               IX0 = IX0 + NXM
               JX = IX0
               DO J = 1,NXM
                  JX = JX + 1
                  IX = IX0
                  DO I = 1,NXM
                     IX = IX + 1
                     GLLKE(IX,JX) = GLLKE(IX,JX) - DSSQX(I,J,IQTB)
                  END DO
               END DO
            END DO
C
C-----------------------------------------------------------------------
C     Perform the inversion of matrix M
C     the output is the scattering path operator -TAU_DELTA(k)
C-----------------------------------------------------------------------
C
C NOTE: DSSQX dummy argument - not used for INVMOD <> 4
C
            IF ( INVMOD.EQ.4 ) CALL STOP_MESSAGE(ROUTINE,'INVMOD = 4')
C
            CALL TBINVERSION(GLLKE,INVMOD,ICHECK,NXM,DSSQX,WK)
C
C-----------------------------------------------------------------------
C  calculate the conductance according to Landauer-Buettiker formalism
C-----------------------------------------------------------------------
C
            CALL SIGLAYLB(IK,IKTOP,KVEC_K,WK_K,WKSUM,GLLKE,GLLIJ,
     &                    GLLIJ_CONJG,MEQ,DSSQX,ISPIN,NSPIN,IQTB_IGIJ,
     &                    JQTB_IGIJ,NGIJ)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
      IF ( ALLOCATED(GLLKE) ) THEN
         DEALLOCATE (GLLKE,STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'dealloc GLLKE')
      END IF
      IF ( ALLOCATED(GLLIJ) ) THEN
         DEALLOCATE (GLLIJ,STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'dealloc GLLIJ')
      END IF
C
      END
C*==siglaylb.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE SIGLAYLB(IK,IKTOP,KVEC_K,WK_K,WKSUM,GLLKE,GLLIJ,
     &                    GLLIJ_CONJG,MEQ,DSSQX,ISPIN,NSPIN,IQTB_IGIJ,
     &                    JQTB_IGIJ,NGIJ)
C   ********************************************************************
C   *                                                                  *
C   * This subroutine calculates the current using the Landauer        *
C   * formula. It uses as input the current matrix elements and the    *
C   * k-resolved green's function.                                     *
C   *                                                                  *
C   *  NOTE:                                                           *
C   *      *  the Green's function convention is   RH                  *
C   *      *  the conductance is calculated using the                  *
C   *         (L,ML,MS) representation - even for IREL = 3             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_SIG,ONLY:NSPINPROJ
      USE MOD_SITES,ONLY:NQTB,QBAS,IQBOT_TB,NQMAX
      USE MOD_ANGMOM,ONLY:NKMMAX,NXM,WKM1,WKM2,WKM3
      USE MOD_CONSTANTS,ONLY:C0,CI,PI,C1
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TB,ONLY:NKKR_TB
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGLAYLB')
      REAL*8 TWOPI
      PARAMETER (TWOPI=2D0*PI)
C
C Dummy arguments
C
      INTEGER IK,IKTOP,ISPIN,NGIJ,NSPIN
      REAL*8 WKSUM
      COMPLEX*16 DSSQX(NKMMAX,NKMMAX,NQTB),GLLIJ(NXM,NXM),
     &           GLLIJ_CONJG(NXM,NXM),GLLKE(NKKR_TB,NKKR_TB),
     &           MEQ(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX)
      INTEGER IQTB_IGIJ(NGIJ),JQTB_IGIJ(NGIJ)
      REAL*8 KVEC_K(3,IKTOP),WK_K(IKTOP)
C
C Local variables
C
      COMPLEX*16 CSCL,CURR_IJ_K(3),CURR_SUM_IJ_K(3)
      REAL*8 CURR_K(:,:),DQJI(3),RWORK(:),SIG0Q(:,:,:),SIG1QQ(:,:,:,:),
     &       SUMCURR(3)
      REAL*8 DDOT
      INTEGER IDEBUG,IGIJ,IKP,ILM,IPOL,IQ,IQBOT_CHI,IQTB,IQTOP_CHI,JLM,
     &        JQ,JQTB,LM1,LM2,N
      LOGICAL INITIALIZE
      SAVE CURR_K,SIG0Q,SIG1QQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CURR_K,RWORK,SIG0Q,SIG1QQ
C
      DATA IDEBUG/0/,INITIALIZE/.TRUE./
C
      IF ( ISPIN.NE.NSPIN ) CALL STOP_MESSAGE(ROUTINE,'ISPIN <> NSPIN')
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (SIG0Q(3,3,NQMAX),SIG1QQ(3,3,NQMAX,NQMAX))
         SIG0Q(1:3,1:3,1:NQMAX) = 0D0
         SIG1QQ(1:3,1:3,1:NQMAX,1:NQMAX) = 0D0
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF ( IK.EQ.1 ) THEN
C
         ALLOCATE (CURR_K(IKTOP,3))
         CURR_K(1:IKTOP,1:3) = 0.D0
C
      END IF
C
C ======================================================================
C                      LOOP OVER PAIRS OF LAYERS
C ======================================================================
C
      CURR_SUM_IJ_K(1:3) = C0
C
      LOOP_IGIJ:DO IGIJ = 1,NGIJ
C
         IQTB = IQTB_IGIJ(IGIJ)
         JQTB = JQTB_IGIJ(IGIJ)
C
         ILM = NKMMAX*(IQTB-1) + 1
         DO LM1 = 1,NKMMAX
            JLM = NKMMAX*(JQTB-1) + LM1
            CALL ZCOPY(NXM,GLLKE(ILM,JLM),1,GLLIJ(1,LM1),1)
         END DO
C
         IQ = IQTB + IQBOT_TB - 1
         JQ = JQTB + IQBOT_TB - 1
C
         DQJI(1:3) = QBAS(1:3,JQ) - QBAS(1:3,IQ)
C
Csbsw      Took away factor of 2pi/a to be in line with Voicu's TBKKR
         CSCL = EXP(-CI*TWOPI*DDOT(3,KVEC_K(1,IK),1,DQJI,1))
                                                            !*TWOPI/ALAT
C
         GLLIJ(1:NXM,1:NXM) = CSCL*GLLIJ(1:NXM,1:NXM)
C
C-----------------------------------------------------------------------
C    convert the site-off-diagonal blocks from TAU_Delta to G_ij
C-----------------------------------------------------------------------
C
C------------- G_ij = 1/(delta t_i) * TAU_Delta *1/(delta t_j)
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,GLLIJ,NXM,DSSQX(1,1,JQTB),
     &              NXM,C0,WKM1,NXM)
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WKM1,NXM,
     &              C0,GLLIJ,NXM)
C
C-----------------------------------------------------------------------
         IF ( IDEBUG.NE.0 ) THEN
            WRITE (6,*) ' IGIJ = ',IGIJ,' IQTB = ',IQTB,' JQTB = ',JQTB
            CALL CMATSTRUCT(' GLL ',GLLIJ,NXM,NXM,IREL,IREL,1,1.D-8,6)
         END IF
C-----------------------------------------------------------------------
C
         DO LM2 = 1,NXM
            DO LM1 = 1,NXM
               GLLIJ_CONJG(LM2,LM1) = DCONJG(GLLIJ(LM1,LM2))
            END DO
         END DO
C
C-----------------------------------------------------------------------
C            calculate current and sum over pairs IQ - JQ
C-----------------------------------------------------------------------
C
         CALL SIGLAYLBCURR(IQ,JQ,GLLIJ,GLLIJ_CONJG,MEQ,CURR_IJ_K,WKM1,
     &                     WKM2,WKM3)
C
         CURR_SUM_IJ_K(1:3) = CURR_SUM_IJ_K(1:3) + CURR_IJ_K(1:3)
C
         DO IPOL = 1,3
            SIG1QQ(IPOL,IPOL,IQ,JQ) = SIG1QQ(IPOL,IPOL,IQ,JQ)
     &                                + DREAL(CURR_IJ_K(IPOL))*WK_K(IK)
     &                                /WKSUM/TWOPI
         END DO
C
      END DO LOOP_IGIJ
C                                                                 pairs
C ======================================================================
C
      CURR_K(IK,1:3) = MAX(0D0,DREAL(CURR_SUM_IJ_K(1:3)))
C
C ======================================================================
      IF ( IK.LT.IKTOP ) RETURN
C ======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_KLOOP ) THEN
C
         CALL DRV_MPI_BARRIER
C
         ALLOCATE (RWORK(IKTOP*3))
C
         CALL DRV_MPI_REDUCE_R(CURR_K,RWORK(1),IKTOP*3)
C
         DEALLOCATE (RWORK)
C
         N = 3*3*NQMAX*NQMAX
C
         ALLOCATE (RWORK(N))
C
         CALL DRV_MPI_REDUCE_R(SIG1QQ,RWORK(1),N)
C
         DEALLOCATE (RWORK)
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C                                                           last k-point
      IF ( MPI_ID.EQ.0 ) THEN
C
         SUMCURR(1:3) = 0D0
         DO IPOL = 1,3
            DO IKP = 1,IKTOP
               CURR_K(IKP,IPOL) = CURR_K(IKP,IPOL)*WK_K(IKP)/WKSUM/TWOPI
Csbssw      CURR_K(IKP) = CURR_K(IKP)
Csbssw      SUMCURR = SUMCURR + CURR_K(IKP)
               SUMCURR(IPOL) = SUMCURR(IPOL) + CURR_K(IKP,IPOL)
            END DO
         END DO
C
         WRITE (*,*) '# of kpts',IKTOP
C        write (6,*) "x-pol"
         WRITE (6,99001) WKSUM,SUMCURR(1),SUMCURR(2),SUMCURR(3)
C        write (6,*) "y-pol"
C        WRITE (6,7777) WKSUM,SUMCURR(2)
C        write (6,*) "z-pol"
C        WRITE (6,7777) WKSUM,SUMCURR(3)
C
C ======================================================================
C        write layer resolved conductivity in case of layered system
C ======================================================================
C
         IQBOT_CHI = IQBOT_TB
         IQTOP_CHI = IQBOT_TB + NQTB - 1
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'sigma_layer.dat')
C
         WRITE (IOTMP,99003) IQBOT_CHI,IQTOP_CHI
         WRITE (IOTMP,99002) 0,((SIG0Q(IPOL,IPOL,IQ),IQ=IQBOT_CHI,
     &                       IQTOP_CHI),IPOL=1,3)
C
         WRITE (IOTMP,99004) IQBOT_CHI,IQTOP_CHI
         DO JQ = IQBOT_CHI,IQTOP_CHI
C
            WRITE (IOTMP,99002) JQ,
     &                          ((SIG1QQ(IPOL,IPOL,IQ,JQ),IQ=IQBOT_CHI,
     &                          IQTOP_CHI),IPOL=1,3)
         END DO
C
         CLOSE (IOTMP)
C
      END IF
C                                                           last k-point
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
99001 FORMAT (/,10X,'VOLUME_sbz           :',F16.10,/,10X,
     &        'Gx/A_BZ/atom (e^2/h)  :',F16.10,/,10X,
     &        'Gy/A_BZ/atom (e^2/h)  :',F16.10,/,10X,
     &        'Gz/A_BZ/atom (e^2/h)  :',F16.10,/)
99002 FORMAT (I3,300E13.5)
99003 FORMAT ('#  sigma0_I   for layers I=',I3,'  to  ',I3,
     &        '     first columns xx, then yy and last zz',/,
     &        '#  dummy index 0 in first column indicates',/,
     &        '#  layer diagonal term')
99004 FORMAT ('#  sigma0_IJ  for layers I=',I3,'  to  ',I3,
     &        '     first columns xx, then yy and last zz',/,
     &        '#  1 line per neighboring layer J ',
     &        '(given in first column) ')
      END
C*==siglaylbcurr.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
C
      SUBROUTINE SIGLAYLBCURR(IQ,JQ,GLLIJ,GLLIJ_CONJG,MEQ,CURR_IJ_K,
     &                        JONE,JTWO,JTEMP)
C   ********************************************************************
C   *                                                                  *
C   *   The sub calculates the contribution to the current for a       *
C   *   single I,J pair and a single k-point                           *
C   *   Input                                                          *
C   *   IQTB, JQTB  : pair of atom layers                              *
C   *       MEQ    : current matrix elements                          *
C   *       GLLIJ   : Green's function matrix element                  *
C   *   Output                                                         *
C   *       CURR_IJ_K  : Current contribution                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SIG,ONLY:NSPINPROJ
      USE MOD_ANGMOM,ONLY:NKMMAX,NXM
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_SITES,ONLY:NQMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,JQ
      COMPLEX*16 CURR_IJ_K(3),GLLIJ(NXM,NXM),GLLIJ_CONJG(NXM,NXM),
     &           JONE(NXM,NXM),JTEMP(NXM),JTWO(NXM,NXM),
     &           MEQ(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX)
C
C Local variables
C
      INTEGER IPOL,ISPINPROJ,LM,LM1,LM2,LM3,LMA,LMB
C
C*** End of declarations rewritten by SPAG
C
      ISPINPROJ = 1
C
      CURR_IJ_K(1:3) = C0
C
      DO IPOL = 1,3
C
C-----------------------------------------------------------------------
C        calculate the current between atom layers IQ and JQ
C-----------------------------------------------------------------------
C
         DO LMA = 1,NXM
            JTEMP(LMA) = C0
            DO LMB = 1,NXM
               JONE(LMA,LMB) = C0
               JTWO(LMA,LMB) = C0
            END DO
         END DO
C
         DO LM2 = 1,NXM
            DO LM1 = 1,NXM
               DO LM = 1,NXM
                  JONE(LM2,LM1) = JONE(LM2,LM1)
     &                            + MEQ(LM,LM2,IPOL,ISPINPROJ,IQ)
     &                            *GLLIJ(LM,LM1)
               END DO
            END DO
         END DO
C
         DO LM1 = 1,NXM
            DO LM2 = 1,NXM
               DO LM3 = 1,NXM
                  JTWO(LM1,LM2) = JTWO(LM1,LM2)
     &                            + MEQ(LM1,LM3,IPOL,ISPINPROJ,JQ)
     &                            *GLLIJ_CONJG(LM3,LM2)
C                        DCONJG(GLLIJ(LM2,LM3))
               END DO
            END DO
         END DO
C
         DO LM2 = 1,NXM
            DO LM1 = 1,NXM
               JTEMP(LM2) = JTEMP(LM2) + JONE(LM2,LM1)*JTWO(LM1,LM2)
            END DO
         END DO
C
         DO LM2 = 1,NXM
            CURR_IJ_K(IPOL) = CURR_IJ_K(IPOL) + JTEMP(LM2)
         END DO
C
      END DO
      END
