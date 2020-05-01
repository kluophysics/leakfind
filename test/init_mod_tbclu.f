C*==init_tbcalc.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE INIT_TBCALC(IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *  initialize the parameters for a TB calculation                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NLMMAX,NL,NKM,NLM,NKMMAX,IXM0_QTB,NXM_QTB,
     &    NKMQ,NLMQ
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF,NQCLU_ICLU_REF,
     &    ICLU_GEO_ICLU_REF,RQCLU_GEO,IQTBORGQTBP,IQTB_IQCLU_ICLU_REF,
     &    IQ_IQCLU_ICLU_REF,NQCLU_GEOMAX
      USE MOD_SYMMETRY,ONLY:IQORGQP,NSYMMAX
      USE MOD_TB,ONLY:NREF,NPLAY,NSLAY_PER_PLAY,VACFLAG,FACTL,INVMOD,
     &    ICHECK,NKKR_TB,NKKRNR_TB,L_R_RELATION,ISYM_IQ_R,IQ_L_IQ_R,
     &    IDECI
      USE MOD_SITES,ONLY:IQBOT_L,IQBOT_R,IQBOT_TB,IQTOP_TB,NQ,NQTB,
     &    IQBOT,IQTOP,NQ_I,NQ_L,NQ_R,NOQ,ITOQ,QBAS
      USE MOD_TYPES,ONLY:CONC,NVALT,NVALTOT,ITBOT,ITTOP,NT
      USE MOD_FILES,ONLY:IOTMP,FOUND_SECTION,FOUND_INTEGER
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE,SUB_SYSTEM,
     &    ABAS_I,ABAS_L
      USE MOD_CALCMODE,ONLY:IREL,PROGNAME
      USE MOD_SCF,ONLY:SCFSTATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_TBCALC')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER IPRINT
C
C Local variables
C
      REAL*8 D,QBAS_R(3),SUML(0:NL)
      REAL*8 DNRM2
      COMPLEX*16 ERYD,GGX,GREF_I1(:,:,:),P,TSSREF(:,:,:)
      INTEGER I,I0,IA_ERR,ICLU_GEO,ICLU_REF,II,II0,IO,IPLAY,IQ,
     &        IQCLU_REF,IQTB,IQ_L,IQ_R,IT,J,JJ,L,NII
      CHARACTER*35 INVTAB(0:4)
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      DATA INVTAB/'FULL MATRIX                        ',
     &     'BANDED MATRIX (slab)               ',
     &     'BANDED + CORNERS MATRIX            ',
     &     'BANDED MATRIX (NOT USED)           ',
     &     'FULL MATRIX: ONLY SITE DIAG. ELEM. '/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TSSREF,GREF_I1
C
C=======================================================================
C            DEALLOCATE arrays in modules if necessary
C=======================================================================
      IF ( ALLOCATED(FACTL) ) DEALLOCATE (FACTL,ICHECK,IQTBORGQTBP,
     &     IQTB_IQCLU_ICLU_REF,ISYM_IQ_R,IQ_L_IQ_R)
C=======================================================================
C
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TSSREF')
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
C
C **********************************************************************
C    calculate the real space TB Green''s function matrices   GREF_I1
C **********************************************************************
C
      IF ( CHECK .OR. IPRINT.GE.1 ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'norm_G_n0.dat')
C
         ERYD = 0.8D0
         P = SQRT(ERYD)
C
         CALL TBGREFRS(ERYD,P,TSSREF,GREF_I1)
C
         DO ICLU_REF = 1,NCLU_REF
C
            ICLU_GEO = ICLU_GEO_ICLU_REF(ICLU_REF)
C
            DO IQCLU_REF = 1,NQCLU_ICLU_REF(ICLU_REF)
C
               SUML(0:NL) = 0D0
               DO L = 0,NL - 1
                  I0 = L**2
                  II0 = (IQCLU_REF-1)*NLM + I0
                  NII = 2*L + 1
                  DO JJ = I0 + 1,I0 + NII
                     DO II = II0 + 1,II0 + NII
                        GGX = GREF_I1(II,JJ,ICLU_REF)
     &                        *DCONJG(GREF_I1(II,JJ,ICLU_REF))
                        SUML(L) = SUML(L) + DREAL(GGX)
                     END DO
                  END DO
                  SUML(L) = SQRT(SUML(L))
               END DO
C
               IF ( CHECK .OR. IPRINT.GE.1 ) THEN
C
                  D = DNRM2(3,RQCLU_GEO(1,IQCLU_REF,ICLU_GEO),1)
C
                  WRITE (IOTMP,'(10E15.6)') D,(SUML(L),L=0,NL-1)
C
C ----------------------------------------------------------------------
C                              XMGRACE - output
C ----------------------------------------------------------------------
C         XMIN = X(1)
C         XMAX = X(NETAB(1))
CC
C         NS = 2
CC
C         DO IOCC = 1,2
C            IT = ITOQ(IOCC,IQCPA)
CC
C            YMAX = 0D0
C            DO IS = 1,NS
C               DO IY = 1,NYOCC(IOCC)
C                  DO IE = 1,NETAB(1)
C                     YMAX = MAX(YMAX,DOS(IE,IY,IS,IOCC))
Cc                  END DO
C               END DO
C            END DO
CC
C            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(E) (sts./eV)'
C            LS = LTXT_T(IT) + 15
C            IF ( NS.EQ.1 ) THEN
C               DTXT1 = 'n!s'//STR20(1:LS)
C               LDTXT1 = 3 + LS
C               DTXT2 = ' '
C               LDTXT2 = 1
C            ELSE
C               DTXT1 = 'n!m{1}!S!UP!M{1}!N!s'//STR20(1:LS)
C               LDTXT1 = 20 + LS
C               DTXT2 = 'n!m{1}!S!DN!M{1}!N!s'//STR20(1:LS)
C               LDTXT2 = 20 + LS
C            END IF
CC
C        CALL XMGRHEAD(DATSET,LDATSET,'CFG_DOS',7,TXT_T(IT),LTXT_T(IT),
C     &                    FILNAM,80,LFILNAM,IOTMP,NSPIN,XMIN,1,XMAX,1,
C     &                    0.0D0,0,YMAX,1,0.0D0,0,YMAX,1,'energy (eV)',
C     &                    11,DTXT1,LDTXT1,DTXT2,LDTXT2,
C     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
C     &                    ,25+LSYSTEM,'NLCPA-DOS of '//TXT_T(IT)
C     &                    (1:LTXT_T(IT))//' in '//SYSTEM(1:LSYSTEM),
C     &                    (13+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
CC
CC         CALL XMGRCURVES(IOTMP,NSPIN,NCURVES,NCURVES,2,1,0)
CC
C            DO IS = 1,NSPIN
C               DO IY = 1,NYOCC(IOCC)
C                  CALL XMGRTABLE((IS-1),IY-1,X,DOS(1,IY,IS,IOCC),1D0,NE,
C     &                           IOTMP)
C               END DO
C            END DO
CC
C            WRITE (6,*) '  '
C            WRITE (6,*) '   DOS written to file  ',FILNAM(1:LFILNAM)
C            CLOSE (IOTMP)
C
               END IF
C
            END DO
         END DO
C
      END IF
      IF ( CHECK .OR. IPRINT.GE.1 ) CLOSE (IOTMP)
C
      IF ( CHECK .OR. IPRINT.GE.1 ) WRITE (*,*)
     &      '********************** TESTING'
C
C=======================================================================
C                 initialize remaining TB - parameters
C=======================================================================
C
      IQBOT_L = 1
      IQBOT_R = 1
C
C=======================================================================
C                              3D-bulk
C=======================================================================
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C
         SUB_SYSTEM = 'BULK'
C
         IQBOT_TB = 1
         NQTB = NQ
         IQBOT = 1
         IQTOP = NQ
C
C------------------------------------------------------ dummy allocation
         ALLOCATE (ISYM_IQ_R(1),IQ_L_IQ_R(1))
C
C=======================================================================
C                            2D-layered system
C=======================================================================
      ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
C-------------------------- force to deal with BULK systems if necessary
         IF ( PROGNAME(4:6).EQ.'SCF' ) THEN
            IF ( SCFSTATUS(1:5).EQ.'START' .AND. SYSTEM_TYPE(1:3)
     &           .EQ.'VIV' ) THEN
               SUB_SYSTEM = 'I-ZONE'
            ELSE IF ( SCFSTATUS(1:5).EQ.'START' .OR. 
     &                SCFSTATUS.EQ.'ITR-L-BULK' ) THEN
               SUB_SYSTEM = 'L-BULK'
            ELSE IF ( SCFSTATUS.EQ.'ITR-R-BULK' ) THEN
               SUB_SYSTEM = 'R-BULK'
            END IF
         END IF
C
         IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
            NQTB = NQ_I
            IQBOT_TB = NQ_L + 1
            IQBOT = NQ_L + 1
            IQTOP = NQ_L + NQ_I
            VACFLAG(1) = SYSTEM_TYPE(1:1).EQ.'V'
            VACFLAG(2) = SYSTEM_TYPE(3:3).EQ.'V'
         ELSE IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
            NQTB = NQ_L
            IQBOT_TB = 1
            IQBOT = 1
            IQTOP = NQ_L
         ELSE IF ( SUB_SYSTEM(1:6).EQ.'R-BULK' ) THEN
            NQTB = NQ_R
            IQBOT_TB = NQ - NQ_R + 1
            IQBOT = NQ - NQ_R + 1
            IQTOP = NQ
         ELSE
            CALL STOP_MESSAGE(ROUTINE,'2D SUB_SYSTEM ???')
         END IF
C
         IQTOP_TB = IQBOT_TB + NQTB + 1
C
C--------------------------------------------------------------- NVALTOT
C
         ITBOT = NT
         ITTOP = 1
         NVALTOT = 0D0
         DO IQ = IQBOT,IQTOP
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               NVALTOT = NVALTOT + CONC(IT)*NVALT(IT)
               ITBOT = MIN(ITBOT,IT)
               ITTOP = MAX(ITTOP,IT)
            END DO
         END DO
C
C-----------------------------------------------------------------------
C                find out relation between L- and R-BULK
C-----------------------------------------------------------------------
         ALLOCATE (ISYM_IQ_R(NQ),IQ_L_IQ_R(NQ))
         ISYM_IQ_R = 0
         IQ_L_IQ_R = 0
C
         IF ( L_R_RELATION.EQ.'UNKNOWN   ' ) THEN
            L_R_RELATION = 'DIFFERENT '
            IF ( NQ_L.EQ.NQ_R ) THEN
               L_R_RELATION = 'IDENTICAL '
               DO IQ_L = 1,NQ_L
                  IQ_R = IQ_L
                  IQ = IQ_R + NQ_L + NQ_I
                  QBAS_R(1:3) = QBAS(1:3,IQ) - ABAS_L(1:3,3)
     &                          - ABAS_I(1:3,3)
                  IF ( .NOT.RVEC_SAME(3,QBAS_R,QBAS(1,IQ_L),1D-5) ) THEN
                     L_R_RELATION = 'DIFFERENT '
                     EXIT
                  END IF
                  IQ_L_IQ_R(IQ) = IQ_L
                  ISYM_IQ_R(IQ) = 1
               END DO
            END IF
         END IF
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'SYSTEM_DIMENSION <> 3, 2 ')
      END IF
C
C
      IF ( IREL.EQ.2 ) THEN
         NKKR_TB = NQTB*NLM
      ELSE
         NKKR_TB = NQTB*NKM
      END IF
      NKKRNR_TB = NQTB*NLM
C
      ALLOCATE (NXM_QTB(NQTB),IXM0_QTB(NQTB))
C
      DO IQTB = 1,NQTB
         IF ( IREL.EQ.2 ) THEN
            NXM_QTB(IQTB) = NLMQ(IQBOT_TB-1+IQTB)
         ELSE
            NXM_QTB(IQTB) = NKMQ(IQBOT_TB-1+IQTB)
         END IF
C
         IF ( IQTB.EQ.1 ) THEN
            IXM0_QTB(IQTB) = 0
         ELSE
            IXM0_QTB(IQTB) = IXM0_QTB(IQTB-1) + NXM_QTB(IQTB-1)
         END IF
      END DO
C
C-----------------------------------------------------------------------
C----------------------------------------------------------  TBINVERSION
C
C  INVMOD = 0  ----> total inversion scheme
C  INVMOD = 1  ----> band matrix inversion scheme
C  INVMOD = 2  ----> corner band matrix inversion scheme
C  INVMOD = 3  ----> sparse matrix inversion scheme
C  INVMOD = 4  ----> total inversion calculating only
C                    diagonal blocks of tau
C
C  NPLAY           number of principle layers
C  NSLAY_PER_PLAY  number of single layers per principle layers
C
C
      INVMOD = 0
      IDECI = 0
C
      CALL INPUT_FIND_SECTION('TAU',0)
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
         NPLAY = 1
         NSLAY_PER_PLAY = NQ
C
         IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER('NSLAY_PER_PLAY',
     &        NSLAY_PER_PLAY,9999,0)
C
         IF ( FOUND_INTEGER ) THEN
            IF ( MOD(NQ,NSLAY_PER_PLAY).NE.0 )
     &           CALL STOP_MESSAGE(ROUTINE,
     &           'NQ/NSLAY_PER_PLAY NOT INTEGER')
            NPLAY = NQ/NSLAY_PER_PLAY
         END IF
C
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
C
         INVMOD = 1
         NSLAY_PER_PLAY = NQ_L
         IDECI = 1
C
         IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER('NSLAY_PER_PLAY',
     &        NSLAY_PER_PLAY,9999,0)
C
         IF ( MOD(NQ_I,NSLAY_PER_PLAY).NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'NQ_I/NSLAY_PER_PLAY NOT INTEGER')
         IF ( MOD(NSLAY_PER_PLAY,NQ_L).NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'NSLAY_PER_PLAY/NQ_L NOT INTEGER')
         IF ( MOD(NSLAY_PER_PLAY,NQ_R).NE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'CHECK NSLAY_PER_PLAY/NQ_R NOT INTEGER')
         NPLAY = NQ_I/NSLAY_PER_PLAY
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
C
         NPLAY = 1
         NSLAY_PER_PLAY = NQ_L
C
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'R-BULK' ) THEN
         NPLAY = 1
         NSLAY_PER_PLAY = NQ_R
      END IF
C
      IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER('INVMOD',INVMOD,
     &     9999,0)
C
      ALLOCATE (FACTL(NKMMAX,NKMMAX),ICHECK(NPLAY,NPLAY),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ICHECK')
C
      ICHECK(1:NPLAY,1:NPLAY) = 0
C
      IF ( INVMOD.EQ.1 ) THEN
         DO IPLAY = 1,NPLAY
            ICHECK(IPLAY,IPLAY) = 1
         END DO
      ELSE IF ( INVMOD.EQ.2 ) THEN
         DO I = 1,NPLAY
            DO J = 1,NPLAY
               IF ( (I.EQ.J) .OR. (J.EQ.NPLAY) .OR. (I.EQ.NPLAY) ) THEN
                  ICHECK(I,J) = 1
               ELSE
                  ICHECK(I,J) = 0
               END IF
            END DO
         END DO
      END IF
C
C=======================================================================
C          tables for shifted site indices    IQ --> IQTB
C=======================================================================
      ALLOCATE (IQTBORGQTBP(NSYMMAX,NQTB))
      ALLOCATE (IQTB_IQCLU_ICLU_REF(NQCLU_GEOMAX,NCLU_REF))
C
      DO IQ = IQBOT,IQTOP
         IQTB = IQ - IQBOT + 1
         IQTBORGQTBP(1:NSYMMAX,IQTB) = IQORGQP(1:NSYMMAX,IQ) - IQBOT + 1
      END DO
C
      DO ICLU_REF = 1,NCLU_REF
         IQTB_IQCLU_ICLU_REF(1:NQCLU_GEOMAX,ICLU_REF)
     &      = IQ_IQCLU_ICLU_REF(1:NQCLU_GEOMAX,ICLU_REF) - IQBOT + 1
      END DO
C=======================================================================
C
      WRITE (6,99001) SYSTEM_DIMENSION,SUB_SYSTEM,ITBOT,ITTOP,IQBOT,
     &                IQTOP,NQTB,NVALTOT,INVTAB(INVMOD),NPLAY,
     &                NSLAY_PER_PLAY,VACFLAG
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) WRITE (6,99002) L_R_RELATION
C
99001 FORMAT (/,1X,79('*'),/,33X,'<INIT_TBCALC>',/,1X,79('*'),//,10X,
     &        'parameters for TB-calculations',/,10X,
     &        'global system dimension            ',A,/,10X,
     &        'TB subsystem to be treated         ',A,/,10X,
     &        'range of global atom types        ITBOT =',I3,
     &        '  --  ITTOP =',I3,/,10X,
     &        'range of global lattice sites     IQBOT =',I3,
     &        '  --  IQTOP =',I3,/,10X,
     &        'sites in TB subsystem             NQTB =',I3,/,10X,
     &        'valence electrons in TB subsystem NVALTOT = ',F10.3,/,
     &        10X,'inversion mode                    ',A,/,10X,
     &        'principle layers                  NPLAY   =',I3,/,10X,
     &        'single in principle layer         NSLAY_PER_PLAY =',I3,/,
     &        10X,'vacuum flags                      VACFLAG =',2L2,/)
99002 FORMAT (10X,'L- to R-BULK relation             ',A)
C
      END
C*==init_mod_tbclu.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE INIT_MOD_TBCLU(IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *  initialize ALL possible TB clusters and                         *
C   *  all structural data needed to set up the  TB G-matrix           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:RMT
      USE MOD_FILES,ONLY:IOTMP,FOUND_INTEGER
      USE MOD_LATTICE,ONLY:ALAT,ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,
     &    ADAINV_R,ADAINV_I,SYSTEM_DIMENSION
      USE MOD_SITES,ONLY:NQ,NQMAX,QBAS,NQ_L,NQ_I,NQ_R
      USE MOD_ANGMOM,ONLY:L_LM,NL,NLMAX,NLQ,NLM,NRGNT123TAB
      USE MOD_CALCMODE,ONLY:MOL
      USE MOD_TBCLU,ONLY:CLURAD_RS,IQREP_ICLU_GEO,NSHLCLU_GEO,
     &    NQCLU_GEOMAX,NQCLU_GEO,IQ_IQCLU_GEO_IQCNTR,ICLU_GEO_IQ,
     &    NCLU_GEO,RSSP_RS,RQCLU_GEO,DSSP_RS,NIJSTAB_RS,IL3RGNT,LM3RGNT,
     &    ISSP2AB_RS,ISSP2BB_RS,IQCLUTAB_RS,JQCLUTAB_RS,NSSP4_RS,
     &    ISSP4_RS,ISDA4_RS,JSDA4_RS,ISSPDIR_RS,ISSP5_RS,NSSP1_RS,
     &    NSSP2A_RS,NSSP2B_RS,NSSPDIR_RS,NSSPABS_RS,CLURYLM_RS,
     &    CLUIPH_RS,IND0QCLU_GEO,IND0QCLU0_RS,FL1L2,NIJSTABMAX,NSSP4MAX,
     &    NKKR_RS,NKKRNR_RSMAX,NLMQCLU_GEO,NCLU_REF_ICLU_GEO,
     &    ICLU_REF_ICLU_GEO,ICLU_GEO_ICLU_REF,IREF_IQCLU_ICLU_REF,
     &    NCLU_REF,ICLU_REF_IQ,NQCLU_ICLU_REF,IQ_IQCLU_ICLU_REF,
     &    R0ICLU_GEO,RI0CLU_GEO,RGNT_CLU,NRGNT,NRGNT12MAX,NRGNT123MAX
      USE MOD_TB,ONLY:NREF,RMTREF,VREF,IREFQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_TBCLU')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER IPRINT
C
C Local variables
C
      REAL*8 CLURAD,DQCLU(:),QBAS_R(:,:),RQCLU(:,:),RQTMP(:,:,:)
      REAL*8 DNRM2
      LOGICAL FOUND
      INTEGER I,I1,I2,I3,I4,I5,I6,I7,IA_ERR,ICLU,ICLU_GEO,ICLU_GEO_SAME,
     &        ICLU_REF,IFLAG,II,IL,IM,IQ,IQBOT_L,IQBOT_R,IQCLU,
     &        IQCLU_GEO,IQCNTR,IQREP_ICLU_REF(:),IQTABTMP(:,:),IQTOP_L,
     &        IQTOP_R,IQ_IQCLU(:),IQ_R,IREF,IREF_IQCLU_IQ(:,:),
     &        IREF_IQCLU_IQ_TMP(:,:),J,JQ,L1,L2,LM,LM1,LM2,LMAX,LRASFIL,
     &        M,M1,M2,N,NKMQCLU(:),NLQCLU(:),NN,NQCLU,NQCLU_I,
     &        NQCLU_IQ(:),NQCLU_L,NQCLU_R,NQTMP(:),NRGNT12,NRGNT123,
     &        NSHLCLU,NSSP1MAX,NSSP2AMAX,NSSP2BMAX,NSSPABSMAX,NSSPDIRMAX
      CHARACTER*80 RASFIL
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      DATA IA_ERR/0/
C
      ALLOCATABLE RQCLU,DQCLU,IQ_IQCLU,RQTMP,NQTMP,IQTABTMP,QBAS_R
      ALLOCATABLE NLQCLU,NKMQCLU,NQCLU_IQ,IREF_IQCLU_IQ
      ALLOCATABLE IQREP_ICLU_REF,IREF_IQCLU_IQ_TMP
C
C=======================================================================
C               set defaults for TB scatterers if necessary
C=======================================================================
      DO IREF = 1,NREF
         IM = IREF
         IF ( ABS(RMTREF(IREF)).LT.1D-8 ) RMTREF(IREF) = RMT(IM)
         IF ( ABS(VREF(IREF)).LT.1D-8 ) VREF(IREF) = 4D0
      END DO
C=======================================================================
C
      WRITE (6,99011)
C
      CALL INPUT_FIND_SECTION('TAU',0)
C
C------------------------------------------ set NSHLCLU_GEO or CLURAD_RS
C
      NSHLCLU_GEO = 0
      CALL SECTION_SET_INTEGER('NSHLCLU',NSHLCLU_GEO,9999,0)
C
      CLURAD_RS = 1.5D0
C         IF ( .NOT.FOUND_INTEGER  ) CALL STOP_MESSAGE(ROUTINE,
C     &       '<INIT_MOD_TBCLU>: cluster radius  CLURAD not specified ')
      IF ( .NOT.FOUND_INTEGER )
     &     CALL SECTION_SET_REAL('CLURAD',CLURAD_RS,9999D0,0)
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
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NRGNT123MAX')
C
      NRGNT12MAX = NL**4
      ALLOCATE (NRGNT(NRGNT12MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NRGNT12MAX')
C
      LMAX = NL - 1
C
      CALL CLUGAUNT_RYLM(LMAX,RGNT_CLU,NRGNT,IL3RGNT,LM3RGNT,NRGNT12,
     &                   NRGNT123,NRGNT12MAX,NRGNT123MAX)
C
      WRITE (6,99012) 'l-expansion      ','NL      ',NL,NLMAX
      WRITE (6,99012) 'Gaunt table I    ','NRGNT12 ',NRGNT12,NRGNT12MAX
      WRITE (6,99012) 'Gaunt table II   ','NRGNT123',NRGNT123,
     &                NRGNT123MAX
C
C------------------------------------------- check consistency of tables
      IFLAG = 0
      DO I = 1,NRGNT123
         LM = LM3RGNT(I)
         IL = IL3RGNT(I)
         IF ( IL.NE.L_LM(LM)+1 ) THEN
            WRITE (*,*) '#########  IL LM L_LM(LM):',IL,LM,L_LM(LM)
            IFLAG = 1
         END IF
      END DO
      IF ( IFLAG.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'IL <> L_LM(LM) + 1')
C
C--------------------------- set up factors i**(l1-l2) and (-1)**(l1-l2)
C
      M = (2*NLMAX-1)**2
C
      ALLOCATE (FL1L2(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FL1L2')
      CALL RINIT(M*M,FL1L2)
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
C=======================================================================
C        find the number  NCLU_GEO  of inequivalent TB-clusters
C              according to the geometry to set up G0
C=======================================================================
C      - create a cluster around each lattice site IQ
C      - compare the resulting clusters
C      - keep only NCLU_GEO inequivalent TB clusters
C=======================================================================
C
      ALLOCATE (NQCLU_IQ(NQMAX))
C
C---------------------------------- set range for sites in L- and R-BULK
C------------------------------------- provide shifted QBAS_R for R-BULK
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
         IQBOT_L = 1
         IQTOP_L = NQ_L
C
         IQBOT_R = NQ - NQ_R + 1
         IQTOP_R = NQ
         ALLOCATE (QBAS_R(3,NQ_R))
         DO IQ_R = 1,NQ_R
            IQ = IQ_R + NQ_L + NQ_I
            QBAS_R(1:3,IQ_R) = QBAS(1:3,IQ) - ABAS_L(1:3,3)
     &                         - ABAS_I(1:3,3)
         END DO
      ELSE
         IQBOT_L = 0
         IQTOP_L = 0
         IQBOT_R = 0
         IQTOP_R = 0
      END IF
C
      NCLU_GEO = 0
C
      DO IQCNTR = 1,NQ
C
C------------------------------------- set cluster parameters to default
C
         CLURAD = CLURAD_RS
         NSHLCLU = NSHLCLU_GEO
C
         IF ( IQCNTR.GE.IQBOT_L .AND. IQCNTR.LE.IQTOP_L ) THEN
C
C---------------------------------------------------- restrict to L-BULK
            CALL CLUSSITES(IOTMP,IPRINT,MOL,'3D        ',ABAS_L,ABAS_L,
     &                     ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,ADAINV_R,
     &                     QBAS,CLURAD,IQCNTR,NQCLU,NQCLU_L,NQCLU_I,
     &                     NQCLU_R,NSHLCLU,NQ_L,0,0,NQMAX)
C
         ELSE IF ( IQCNTR.GE.IQBOT_R .AND. IQCNTR.LE.IQTOP_R ) THEN
C
C---------------------------------------------------- restrict to R-BULK
C
            IQ_R = IQCNTR - NQ_L - NQ_I
C
            CALL CLUSSITES(IOTMP,IPRINT,MOL,'3D        ',ABAS_R,ABAS_L,
     &                     ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,ADAINV_R,
     &                     QBAS_R,CLURAD,IQ_R,NQCLU,NQCLU_L,NQCLU_I,
     &                     NQCLU_R,NSHLCLU,NQ_R,0,0,NQMAX)
C
         ELSE
C
C------------------------------------------------ deal with total system
C-------------------------------------- to have I-ZONE properly embedded
C
            CALL CLUSSITES(IOTMP,IPRINT,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,CLURAD,IQCNTR,NQCLU,NQCLU_L,
     &                     NQCLU_I,NQCLU_R,NSHLCLU,NQ,NQ_L,NQ_R,NQMAX)
C
         END IF
C
C -------------------------------------------- allocate storage and read
C
         IF ( ALLOCATED(RQCLU) ) DEALLOCATE (RQCLU,DQCLU,IQ_IQCLU)
C
         ALLOCATE (RQCLU(3,NQCLU),DQCLU(NQCLU),IQ_IQCLU(NQCLU))
C
         READ (IOTMP) ((RQCLU(J,I),J=1,3),DQCLU(I),IQ_IQCLU(I),I=1,
     &                NQCLU)
         CLOSE (IOTMP)
C
         IF ( IQCNTR.GE.IQBOT_R .AND. IQCNTR.LE.IQTOP_R )
     &        IQ_IQCLU(1:NQCLU) = IQ_IQCLU(1:NQCLU) + NQ_L + NQ_I
C
         IF ( IPRINT.GE.5 ) THEN
            WRITE (6,99007) IQCNTR
            DO I = 1,NQCLU
               WRITE (6,99008) I,(RQCLU(J,I),J=1,3),DQCLU(I),IQ_IQCLU(I)
            END DO
         END IF
C
C--------------------------------------------------- store first cluster
         IF ( NCLU_GEO.EQ.0 ) THEN
            NCLU_GEO = 1
            NQCLU_GEOMAX = 10*NQCLU
C
            ALLOCATE (RQTMP(3,NQCLU_GEOMAX,NQ))
            ALLOCATE (NQTMP(NQCLU_GEOMAX),IQREP_ICLU_GEO(NQCLU_GEOMAX))
            ALLOCATE (IQTABTMP(NQCLU_GEOMAX,NQ),ICLU_GEO_IQ(NQMAX))
            ALLOCATE (IREF_IQCLU_IQ_TMP(NQCLU_GEOMAX,NQ))
C
            RQTMP(1:3,1:NQCLU,NCLU_GEO) = RQCLU(1:3,1:NQCLU)
            NQTMP(NCLU_GEO) = NQCLU
            IQREP_ICLU_GEO(NCLU_GEO) = IQCNTR
            IQTABTMP(1:NQCLU,IQCNTR) = IQ_IQCLU(1:NQCLU)
            ICLU_GEO_IQ(IQCNTR) = NCLU_GEO
C---------------------------------------- compare with previous clusters
         ELSE
            ICLU_GEO_SAME = 0
            DO ICLU_GEO = 1,NCLU_GEO
               IF ( NQTMP(ICLU_GEO).EQ.NQCLU ) THEN
                  IF ( RVEC_SAME(3*NQCLU,RQTMP(1,1,ICLU_GEO),RQCLU,1D-6)
     &                 ) THEN
                     ICLU_GEO_SAME = ICLU_GEO
                     EXIT
                  END IF
               END IF
            END DO
C
            IF ( ICLU_GEO_SAME.NE.0 ) THEN
               ICLU_GEO_IQ(IQCNTR) = ICLU_GEO_SAME
               IQTABTMP(1:NQCLU,IQCNTR) = IQ_IQCLU(1:NQCLU)
               WRITE (6,99009) IQCNTR,IQREP_ICLU_GEO(ICLU_GEO),ICLU_GEO
C
            ELSE
C
               NCLU_GEO = NCLU_GEO + 1
               IF ( NCLU_GEO.GT.NQCLU_GEOMAX )
     &              CALL STOP_MESSAGE(ROUTINE,'NCLU_GEO > NQCLU_GEOMAX')
               RQTMP(1:3,1:NQCLU,NCLU_GEO) = RQCLU(1:3,1:NQCLU)
               NQTMP(NCLU_GEO) = NQCLU
               IQREP_ICLU_GEO(NCLU_GEO) = IQCNTR
               IQTABTMP(1:NQCLU,IQCNTR) = IQ_IQCLU(1:NQCLU)
               ICLU_GEO_IQ(IQCNTR) = NCLU_GEO
            END IF
         END IF
C-----------------------------------------------------------------------
C
         NQCLU_IQ(IQCNTR) = NQCLU
C
         DO IQCLU = 1,NQCLU
            IQ = IQ_IQCLU(IQCLU)
            IREF_IQCLU_IQ_TMP(IQCLU,IQCNTR) = IREFQ(IQ)
         END DO
C
      END DO
C
C
C=======================================================================
C              find actual array sizes needed, allocate  and
C            swap data from temporary variables to proper ones
C=======================================================================
C
      NQCLU_GEOMAX = 0
      DO ICLU_GEO = 1,NCLU_GEO
         NQCLU_GEOMAX = MAX(NQCLU_GEOMAX,NQTMP(ICLU_GEO))
      END DO
C
      ALLOCATE (RI0CLU_GEO(3,NQCLU_GEOMAX,NCLU_GEO))
      ALLOCATE (R0ICLU_GEO(3,NQCLU_GEOMAX,NCLU_GEO))
      ALLOCATE (RQCLU_GEO(3,NQCLU_GEOMAX,NCLU_GEO),NQCLU_GEO(NCLU_GEO))
      ALLOCATE (IQ_IQCLU_GEO_IQCNTR(NQCLU_GEOMAX,NQMAX))
C
      NQCLU_GEO(1:NCLU_GEO) = NQTMP(1:NCLU_GEO)
      DO ICLU_GEO = 1,NCLU_GEO
         NQCLU = NQCLU_GEO(ICLU_GEO)
         RQCLU_GEO(1:3,1:NQCLU,ICLU_GEO) = RQTMP(1:3,1:NQCLU,ICLU_GEO)
         R0ICLU_GEO(1:3,1:NQCLU,ICLU_GEO) = RQTMP(1:3,1:NQCLU,ICLU_GEO)
         RI0CLU_GEO(1:3,1:NQCLU,ICLU_GEO) = -RQTMP(1:3,1:NQCLU,ICLU_GEO)
      END DO
C
C
      IQ_IQCLU_GEO_IQCNTR(1:NQCLU_GEOMAX,1:NQ)
     &   = IQTABTMP(1:NQCLU_GEOMAX,1:NQ)
C
      WRITE (6,99010) NCLU_GEO,(NQCLU_GEO(ICLU_GEO),ICLU_GEO=1,NCLU_GEO)
C
      DEALLOCATE (RQTMP,NQTMP,IQTABTMP)
C
C
      ALLOCATE (IREF_IQCLU_IQ(NQCLU_GEOMAX,NQ))
C
      IREF_IQCLU_IQ(1:NQCLU_GEOMAX,1:NQ)
     &   = IREF_IQCLU_IQ_TMP(1:NQCLU_GEOMAX,1:NQ)
C
C=====================================================================IQ
C
      IF ( CHECK ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'SPRKKR-TB-cluster')
C
         WRITE (IOTMP,'(''NCLU_GEO ='',I4)') NCLU_GEO
C
         DO ICLU_GEO = 1,NCLU_GEO
            WRITE (IOTMP,'(''NQCLU_GEO ='',I4)') NQCLU_GEO(ICLU_GEO)
            DO IQCLU = 1,NQCLU_GEO(ICLU_GEO)
               WRITE (6,99008) IQCLU,(RQCLU(J,IQCLU),J=1,3),
     &                         DNRM2(3,RQCLU(1,IQCLU),1)
            END DO
         END DO
C
         CLOSE (IOTMP)
      END IF
C
C=======================================================================
C             set up tables specifying blocks of G-matrix
C=======================================================================
C
C-----------------------------------------------------------------------
C               find first the maximum array sizes
C-----------------------------------------------------------------------
      NSSP1MAX = 0
      NIJSTABMAX = 0
      NSSP2AMAX = 0
      NSSP2BMAX = 0
      NSSPABSMAX = 0
      NSSP4MAX = 0
      NSSPDIRMAX = 0
C
      DO ICLU_GEO = 1,NCLU_GEO
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
         CALL CLUSIT2(IOTMP,NQCLU_GEO(ICLU_GEO),RQCLU_GEO(1,1,ICLU_GEO),
     &                0,NQCLU_GEO(ICLU_GEO))
C
         REWIND IOTMP
         READ (IOTMP) I1,I2,I3,I4,I5,I6,I7
C
         NSSP1MAX = MAX(NSSP1MAX,I1)
         NIJSTABMAX = MAX(NIJSTABMAX,I2)
         NSSP2AMAX = MAX(NSSP2AMAX,I3)
         NSSP2BMAX = MAX(NSSP2BMAX,I4)
         NSSPABSMAX = MAX(NSSPABSMAX,I5)
         NSSP4MAX = MAX(NSSP4MAX,I6)
         NSSPDIRMAX = MAX(NSSPDIRMAX,I7)
C
         CLOSE (IOTMP)
C
      END DO
C
C-----------------------------------------------------------------------
C       allocate storage for geometry information on all TB clusters
C-----------------------------------------------------------------------
C
      ALLOCATE (RSSP_RS(3,0:NSSP1MAX,NCLU_GEO))
      ALLOCATE (DSSP_RS(0:NSSP1MAX,NCLU_GEO))
      ALLOCATE (IQCLUTAB_RS(NIJSTABMAX,NSSP1MAX,NCLU_GEO))
      ALLOCATE (JQCLUTAB_RS(NIJSTABMAX,NSSP1MAX,NCLU_GEO))
      ALLOCATE (NIJSTAB_RS(NSSP1MAX,NCLU_GEO))
      ALLOCATE (ISSP2AB_RS(NSSP2AMAX,NCLU_GEO))
      ALLOCATE (ISSP2BB_RS(NSSP2AMAX,NCLU_GEO))
      ALLOCATE (NSSP4_RS(NSSPABSMAX,NCLU_GEO))
      ALLOCATE (ISSP4_RS(NSSP4MAX,NSSPABSMAX,NCLU_GEO))
      ALLOCATE (ISDA4_RS(NSSP4MAX,NSSPABSMAX,NCLU_GEO))
      ALLOCATE (JSDA4_RS(NSSP4MAX,NSSPABSMAX,NCLU_GEO))
      ALLOCATE (ISSPDIR_RS(NSSPDIRMAX,NCLU_GEO),NKKR_RS(NCLU_GEO))
      ALLOCATE (NSSP1_RS(NCLU_GEO),NSSP2A_RS(NCLU_GEO))
      ALLOCATE (NSSP2B_RS(NCLU_GEO))
      ALLOCATE (NSSPDIR_RS(NCLU_GEO),NSSPABS_RS(NCLU_GEO))
      ALLOCATE (NLMQCLU_GEO(NQCLU_GEOMAX,NCLU_GEO))
      ALLOCATE (ISSP5_RS(NSSP4MAX,NSSPABSMAX,NCLU_GEO),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ISSP5')
      ALLOCATE (CLUIPH_RS(2*NLMAX-1,NSSPABSMAX,NCLU_GEO))
      ALLOCATE (CLURYLM_RS((2*NLMAX-1)**2,NSSPDIRMAX,NCLU_GEO),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CLUIPH')
C
C------------------------------------------ local temporary dummy arrays
C
      ALLOCATE (NLQCLU(NQCLU_GEOMAX),NKMQCLU(NQCLU_GEOMAX))
      ALLOCATE (IND0QCLU0_RS(NQCLU_GEOMAX,NCLU_GEO))
      ALLOCATE (IND0QCLU_GEO(NQCLU_GEOMAX,NCLU_GEO),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ITOQCLU')
C
C***********************************************************************
      RSSP_RS = 999999
      DSSP_RS = 999999
      IQCLUTAB_RS = 999999
      JQCLUTAB_RS = 999999
      NIJSTAB_RS = 999999
      ISSP2AB_RS = 999999
      ISSP2BB_RS = 999999
      NSSP4_RS = 999999
      ISSP4_RS = 999999
      ISDA4_RS = 999999
      JSDA4_RS = 999999
      ISSPDIR_RS = 999999
      NKKR_RS = 999999
      NSSP1_RS = 999999
      NSSP2A_RS = 999999
      NSSP2B_RS = 999999
      NSSPDIR_RS = 999999
      NSSPABS_RS = 999999
      NLMQCLU_GEO = 999999
      ISSP5_RS = 999999
      CLUIPH_RS = 999999
      CLURYLM_RS = 999999
      NLQCLU = 999999
      NKMQCLU = 999999
      IND0QCLU0_RS = 999999
      IND0QCLU_GEO = 999999
C***********************************************************************
C
C-----------------------------------------------------------------------
C               set geometry information on all TB clusters
C-----------------------------------------------------------------------
C
      NKKRNR_RSMAX = 0
C
      DO ICLU_GEO = 1,NCLU_GEO
C
         IQ = IQREP_ICLU_GEO(ICLU_GEO)
C
         RASFIL = 'rasmol_TB_cluster_'
         CALL STRING_ADD_N(RASFIL,ICLU_GEO)
         LRASFIL = LEN_TRIM(RASFIL)
C
         CALL CLUINITSTEP2(IPRINT,IOTMP,RASFIL,LRASFIL,ALAT,
     &                     RQCLU_GEO(1,1,ICLU_GEO),RSSP_RS(1,0,ICLU_GEO)
     &                     ,DSSP_RS(0,ICLU_GEO),
     &                     IQCLUTAB_RS(1,1,ICLU_GEO),
     &                     JQCLUTAB_RS(1,1,ICLU_GEO),
     &                     NIJSTAB_RS(1,ICLU_GEO),ISSP2AB_RS(1,ICLU_GEO)
     &                     ,ISSP2BB_RS(1,ICLU_GEO),NSSP4_RS(1,ICLU_GEO),
     &                     ISSP4_RS(1,1,ICLU_GEO),ISDA4_RS(1,1,ICLU_GEO)
     &                     ,JSDA4_RS(1,1,ICLU_GEO),
     &                     ISSPDIR_RS(1,ICLU_GEO),ISSP5_RS(1,1,ICLU_GEO)
     &                     ,IND0QCLU0_RS(1,ICLU_GEO),
     &                     IND0QCLU_GEO(1,ICLU_GEO),NQ,
     &                     IQ_IQCLU_GEO_IQCNTR(1,IQ),NL,NLQ,NLQCLU,
     &                     NLMQCLU_GEO(1,ICLU_GEO),NKMQCLU,
     &                     CLURYLM_RS(1,1,ICLU_GEO),NSSP1_RS(ICLU_GEO),
     &                     NSSP2A_RS(ICLU_GEO),NSSP2B_RS(ICLU_GEO),
     &                     NSSPDIR_RS(ICLU_GEO),NSSPABS_RS(ICLU_GEO),
     &                     NKKR_RS(ICLU_GEO),NSSP1MAX,NSSP2AMAX,
     &                     NSSP4MAX,NSSPABSMAX,NSSPDIRMAX,NIJSTABMAX,
     &                     NQCLU_GEO(ICLU_GEO),NLMAX,NQMAX)
C
C--------------------- use the same l_max for all sites and all clusters
C
         NKKR_RS(ICLU_GEO) = NLM*NQCLU_GEO(ICLU_GEO)
         NKKRNR_RSMAX = MAX(NKKRNR_RSMAX,NKKR_RS(ICLU_GEO))
         DO IQCLU_GEO = 1,NQCLU_GEO(ICLU_GEO)
            NLMQCLU_GEO(IQCLU_GEO,ICLU_GEO) = NLM
         END DO
C
      END DO
C
C=======================================================================
C        find the number  NCLU_REF  of inequivalent TB-clusters
C              according to the occupation to set up GREF_I1
C=======================================================================
      ALLOCATE (IQREP_ICLU_REF(NQMAX),ICLU_REF_IQ(NQMAX))
C
      IQ = 1
      ICLU_REF = 1
      NCLU_REF = ICLU_REF
      IQREP_ICLU_REF(NCLU_REF) = IQ
      ICLU_REF_IQ(IQ) = NCLU_REF
C
      DO IQ = 1,NQ
C
         DO ICLU_REF = 1,NCLU_REF
            JQ = IQREP_ICLU_REF(ICLU_REF)
C
            IF ( NQCLU_IQ(IQ).EQ.NQCLU_IQ(JQ) ) THEN
C------------------ compare cluster site occupation with reference atoms
               DO IQCLU = 1,NQCLU_IQ(IQ)
                  IF ( IREF_IQCLU_IQ(IQCLU,IQ)
     &                 .NE.IREF_IQCLU_IQ(IQCLU,JQ) ) GOTO 50
               END DO
C------------------- cluster around sites IQ and JQ have SAME occupation
               ICLU_REF_IQ(IQ) = ICLU_REF
               GOTO 100
C
C-------------- cluster around sites IQ and JQ have DIFFERENT occupation
            END IF
C
 50      END DO
C--------------------------------------------------- add new REF cluster
         NCLU_REF = NCLU_REF + 1
         IQREP_ICLU_REF(NCLU_REF) = IQ
         ICLU_REF_IQ(IQ) = NCLU_REF
C
 100  END DO
C
C=======================================================================
C        find the number of  NCLU_REF_ICLU_GEO  reference clusters
C             connected with a geometry cluster  ICLU_GEO
C=======================================================================
C
      ALLOCATE (NCLU_REF_ICLU_GEO(NCLU_GEO))
      ALLOCATE (ICLU_REF_ICLU_GEO(NCLU_REF,NCLU_GEO))
C
      NCLU_REF_ICLU_GEO(1:NCLU_GEO) = 0
      ICLU_REF_ICLU_GEO(1:NCLU_REF,1:NCLU_GEO) = 999999
C
      DO IQ = 1,NQ
         ICLU_GEO = ICLU_GEO_IQ(IQ)
         ICLU_REF = ICLU_REF_IQ(IQ)
C
         FOUND = .FALSE.
         DO ICLU = 1,NCLU_REF_ICLU_GEO(ICLU_GEO)
            IF ( ICLU_REF_ICLU_GEO(ICLU,ICLU_GEO).EQ.ICLU_REF )
     &           FOUND = .TRUE.
         END DO
         IF ( .NOT.FOUND ) THEN
            NN = NCLU_REF_ICLU_GEO(ICLU_GEO) + 1
            NCLU_REF_ICLU_GEO(ICLU_GEO) = NN
            ICLU_REF_ICLU_GEO(NN,ICLU_GEO) = ICLU_REF
         END IF
C
      END DO
C
      ALLOCATE (ICLU_GEO_ICLU_REF(NCLU_REF))
C
      DO ICLU_GEO = 1,NCLU_GEO
         DO II = 1,NCLU_REF_ICLU_GEO(ICLU_GEO)
            ICLU_REF = ICLU_REF_ICLU_GEO(II,ICLU_GEO)
            ICLU_GEO_ICLU_REF(ICLU_REF) = ICLU_GEO
         END DO
      END DO
C
      ALLOCATE (IREF_IQCLU_ICLU_REF(NQCLU_GEOMAX,NCLU_REF))
      ALLOCATE (NQCLU_ICLU_REF(NCLU_REF))
      ALLOCATE (IQ_IQCLU_ICLU_REF(NQCLU_GEOMAX,NCLU_REF))
C
      DO ICLU_REF = 1,NCLU_REF
         IQCNTR = IQREP_ICLU_REF(ICLU_REF)
C
         IREF_IQCLU_ICLU_REF(1:NQCLU_GEOMAX,ICLU_REF)
     &      = IREF_IQCLU_IQ(1:NQCLU_GEOMAX,IQCNTR)
C
         ICLU_GEO = ICLU_GEO_ICLU_REF(ICLU_REF)
C
         NQCLU_ICLU_REF(ICLU_REF) = NQCLU_IQ(IQCNTR)
C
         IQ_IQCLU_ICLU_REF(1:NQCLU_GEOMAX,ICLU_REF)
     &      = IQ_IQCLU_GEO_IQCNTR(1:NQCLU_GEOMAX,IQCNTR)
C
      END DO
C
C
C=======================================================================
C                         write overview
C=======================================================================
C
      WRITE (6,99003)
      IF ( CHECK .OR. IPRINT.GE.1 ) THEN
         DO IQ = 1,NQ
            WRITE (6,99002) IQ,NQCLU_IQ(IQ),ICLU_GEO_IQ(IQ),
     &                      ICLU_REF_IQ(IQ)
            IF ( IPRINT.GE.0 ) THEN
               WRITE (6,99001) IQ
               DO IQCLU = 1,NQCLU_IQ(IQ)
                  WRITE (6,99002) IQ,IQCLU,IREF_IQCLU_IQ(IQCLU,IQ)
               END DO
            END IF
         END DO
      END IF
C
      WRITE (6,99005)
C
      DO ICLU_GEO = 1,NCLU_GEO
         NN = NCLU_REF_ICLU_GEO(ICLU_GEO)
         WRITE (6,99004) ICLU_GEO,NN,
     &                   (ICLU_REF_ICLU_GEO(II,ICLU_GEO),II=1,NN)
      END DO
C
      WRITE (6,99006)
C
      DO ICLU_REF = 1,NCLU_REF
         WRITE (6,99004) ICLU_REF,ICLU_GEO_ICLU_REF(ICLU_REF)
      END DO
C
C=======================================================================
99001 FORMAT (/,10X,'occupation of cluster around site IQCNTR = ',I5,/,
     &        10X,'IQCNTR          IQCLU           IREF  ')
99002 FORMAT (10X,I5,3I15)
99003 FORMAT (/,10X,'TB cluster information ',//,10X,
     &        'IQCNTR          NQCLU         ICLU_GEO       ICLU_REF')
99004 FORMAT (10X,I5,I11,6X,8I5,:,/,(32X,8I5))
99005 FORMAT (/,10X,'geometry - reference cluster relationship',//,10X,
     &        'ICLU_GEO   NCLU_REF   ICLU_REF')
99006 FORMAT (/,10X,'ICLU_REF   ICLU_GEO')
99007 FORMAT (/,10X,'cluster around site IQ =',I3,/,
     &     '     IQCLU      RX        RY        RZ          D        IQ'
     &     )
99008 FORMAT (I10,2X,3F10.5,2X,F10.5,I5)
99009 FORMAT (/,10X,'cluster around site IQ =',I4,' same as for JQ =',
     &        I4,':   ICLU_GEO = ',I4)
99010 FORMAT (/,10X,'number of TB-clusters     NCLU_GEO  =',I5,/,10X,
     &        'with cluster sites        NQCLU_GEO =',5I5,/,:,(10X,14I5)
     &        )
99011 FORMAT (/,1X,79('*'),/,27X,'<INIT_MOD_TBCLU>',/,1X,79('*'),/)
99012 FORMAT (10X,A,9X,A,I10,:,7X,'(',I7,')')
      END
C*==cluinitstep2.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE CLUINITSTEP2(IPRINT,IOTMP,RASFIL,LRASFIL,ALAT,RQCLU,
     &                        RSSP,DSSP,IQCLUTAB,JQCLUTAB,NIJSTAB,
     &                        ISSP2AB,ISSP2BB,NSSP4,ISSP4,ISDA4,JSDA4,
     &                        ISSPDIR,ISSP5,IND0QCLU0,IND0QCLU,NQ,
     &                        IQ_IQCLU,NL,NLQ,NLQCLU,NLMQCLU,NKMQCLU,
     &                        CLURYLM_NNP,NSSP1,NSSP2A,NSSP2B,NSSPDIR,
     &                        NSSPABS,NKKR,NSSP1MAX,NSSP2AMAX,NSSP4MAX,
     &                        NSSPABSMAX,NSSPDIRMAX,NIJSTABMAX,NQCLU,
     &                        NLMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:A0_ANG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUINITSTEP2')
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER IOTMP,IPRINT,LRASFIL,NIJSTABMAX,NKKR,NL,NLMAX,NQ,NQCLU,
     &        NQMAX,NSSP1,NSSP1MAX,NSSP2A,NSSP2AMAX,NSSP2B,NSSP4MAX,
     &        NSSPABS,NSSPABSMAX,NSSPDIR,NSSPDIRMAX
      CHARACTER*80 RASFIL
      REAL*8 CLURYLM_NNP((2*NLMAX-1)**2,NSSPDIRMAX),DSSP(0:NSSP1MAX),
     &       RQCLU(3,NQCLU),RSSP(3,0:NSSP1MAX)
      INTEGER IND0QCLU(NQCLU),IND0QCLU0(NQCLU),
     &        IQCLUTAB(NIJSTABMAX,NSSP1MAX),IQ_IQCLU(NQCLU),
     &        ISDA4(NSSP4MAX,NSSPABSMAX),ISSP2AB(NSSP2AMAX),
     &        ISSP2BB(NSSP2AMAX),ISSP4(NSSP4MAX,NSSPABSMAX),
     &        ISSP5(NSSP4MAX,NSSPABSMAX),ISSPDIR(NSSPDIRMAX),
     &        JQCLUTAB(NIJSTABMAX,NSSP1MAX),JSDA4(NSSP4MAX,NSSPABSMAX),
     &        NIJSTAB(NSSP1MAX),NKMQCLU(NQCLU),NLMQCLU(NQCLU),NLQ(NQMAX)
     &        ,NLQCLU(NQCLU),NSSP4(NSSPABSMAX)
C
C Local variables
C
      REAL*8 DIR(3),S
      INTEGER I,I1,I2,I2B,IA,ID,IDIR,IG,IGIJ(:,:),IQ,IQCLU,IQCLU1,
     &        IQCLU2,ISSP,IX,J,JQCLU,JQCLU1,JQCLU2,L3MAX,NGIJ,
     &        NIJSTABMAX_DUM,NLM,NLM3MAX,NSSP4MAX_DUM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IGIJ
C
      ALLOCATE (IGIJ(NQCLU,NQCLU))
C
C--------------------------------------------------------- supply NLQCLU
C
      DO IQCLU = 1,NQCLU
         NLQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))
      END DO
C
C--------------------------------------- set up tables specifying blocks
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
      CALL CLUSIT2(IOTMP,NQCLU,RQCLU,IPRINT,NQCLU)
C
      REWIND IOTMP
      READ (IOTMP) NSSP1,NIJSTABMAX_DUM,NSSP2A,NSSP2B,NSSPABS,
     &             NSSP4MAX_DUM,NSSPDIR
C
      READ (IOTMP) ((RSSP(J,I),J=1,3),DSSP(I),NIJSTAB(I),I=1,NSSP1),
     &             ((IQCLUTAB(J,I),JQCLUTAB(J,I),J=1,NIJSTAB(I)),I=1,
     &             NSSP1),(ISSP2AB(I),ISSP2BB(I),I=1,NSSP2A),
     &             (NSSP4(I),I=1,NSSPABS),
     &             ((ISSP4(J,I),ISDA4(J,I),JSDA4(J,I),J=1,NSSP4(I)),I=1,
     &             NSSPABS),(ISSPDIR(I),I=1,NSSPDIR),
     &             ((ISSP5(J,I),J=1,NSSP4(I)),I=1,NSSPABS)
      CLOSE (IOTMP)
C
      WRITE (6,99002) 'S-S'' combinations','NSSP1   ',NSSP1
      WRITE (6,99002) 'S-S'' lengths     ','NSSPABS ',NSSPABS
      WRITE (6,99002) 'S-S'' directions  ','NSSPDIR ',NSSPDIR
      IF ( IPRINT.GE.3 ) THEN
         WRITE (6,99002) 'S-S''             ','NIJSTAB ',NIJSTABMAX_DUM
         WRITE (6,99002) 'S-S''             ','NSSP4   ',NSSP4MAX_DUM
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
      WRITE (6,99002) 'number of G(I,J) ','NGIJ    ',NGIJ
      WRITE (6,99001)
C
      IF ( IPRINT.GT.1 ) THEN
         DO IQCLU = 1,NQCLU
            DO JQCLU = 1,IQCLU - 1
               WRITE (6,'(10X,A,2I4,2(A,I4))') 'IGIJ',IQCLU,JQCLU,
     &                '  i,j:',IGIJ(IQCLU,JQCLU),' j,i:',
     &                IGIJ(JQCLU,IQCLU)
            END DO
         END DO
      END IF
C
C --------------------------------------------------------- NLM-blocking
C
      NKKR = 0
      NLM = NL**2
C
      DO IQCLU = 1,NQCLU
         NLQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))
         NLMQCLU(IQCLU) = NLQ(IQ_IQCLU(IQCLU))**2
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
      L3MAX = 2*(NLMAX-1)
      NLM3MAX = (L3MAX+1)*(L3MAX+1)
C
      DO IDIR = 1,NSSPDIR
         ISSP = ISSPDIR(IDIR)
         IF ( DABS(DSSP(ISSP)).GE.1D-6 ) THEN
            DO I = 1,3
               DIR(I) = RSSP(I,ISSP)/DSSP(ISSP)
            END DO
         ELSE
            DO I = 1,3
               DIR(I) = 0D0
            END DO
         END IF
C
         CLURYLM_NNP(1:NLM3MAX,IDIR) = CLURYLM_NNP(1:NLM3MAX,1)
C
         IF ( IDIR.NE.1 ) CALL CALC_RHPLM(DIR(1),DIR(2),DIR(3),
     &        CLURYLM_NNP(1,IDIR),L3MAX,NLM3MAX)
C
      END DO
      IF ( IPRINT.GE.1 ) THEN
C
C----------------------------------------------------- write RASMOL data
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,RASFIL(1:LRASFIL)//'.pdb')
         WRITE (IOTMP,FMT=99004)
C
         S = 3*A0_ANG*ALAT
         I = 0
         DO IQCLU = 1,NQCLU
            IQ = IQ_IQCLU(IQCLU)
            I = I + 1
            WRITE (IOTMP,FMT=99003) I,I,(S*RQCLU(IX,IQCLU),IX=1,3),
     &                              DBLE(IQ)*3D0/DBLE(NQ)
         END DO
C
         CLOSE (IOTMP)
C
C--------------------------------------------------- write RASMOL script
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,RASFIL(1:LRASFIL)//'.ras')
         WRITE (IOTMP,*) 'load '''//RASFIL(1:LRASFIL)//'.pdb'' '
         WRITE (IOTMP,*) 'color temperature'
         WRITE (IOTMP,*) 'set fontsize 20'
         WRITE (IOTMP,*) 'cpk   450'
         WRITE (IOTMP,*) 'select all '
         WRITE (IOTMP,*) 'set axes on '
C
         CLOSE (IOTMP)
C
         WRITE (6,99005) RASFIL(1:LRASFIL)//'.ras'
      END IF
C
C=======================================================================
99001 FORMAT (1X,79('*'),/)
99002 FORMAT (10X,A,9X,A,I10,:,7X,'(',I7,')')
99003 FORMAT ('ATOM   ',I4,'           ',I4,'    ',3F8.3,'  0.00',F8.3)
99004 FORMAT ('HEADER    cluster confguration set up by <CLUSTER>',/,
     &        'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99005 FORMAT (/,10X,'cluster configuration stored in data-file ',
     &        ' rasmol_cluster.pdb',/,10X,
     &        'view via:   rasmol  -script ',A,/)
      END
