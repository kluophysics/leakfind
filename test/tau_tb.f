C*==tau_tb.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_TB(ERYD,P,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,
     &                  MSST,TSSQ,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   *  - call <TBKLOOP> to perform the appropriate TB k-integration    *
C   *    and return TAU                                                *
C   *                                                                  *
C   *  - perform the CPA-cyle if necessary                             *
C   *                                                                  *
C   *  NOTE: this routines works for IREL=0,1,2 and 3                  *
C   *        in particular for the spin-polarised scalar relativistic  *
C   *        case many temporary arrays are needed "                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX,NCPA
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,ICPA,IQBOT_L,IQBOT_R,
     &    IQBOT_TB,NQTB,IQBOT,IQTOP
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY,NSYMMAX
      USE MOD_TYPES,ONLY:CONC,NTMAX,NT,ITBOT,ITTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,NLMMAX,NKMQ,NLMQ,NSPIN,NXM
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF,IQTBORGQTBP
      USE MOD_CALCMODE,ONLY:KMROT,IREL
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      USE MOD_SCF,ONLY:SCFSTATUS
      IMPLICIT NONE
C*--TAU_TB34
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 CPACORR,CPAERR,CPAERRL
      COMPLEX*16 DROTX(:,:,:),DSSQX(:,:,:),GREF_I1(:,:,:),MSSQX(:,:,:),
     &           MSSTX(:,:,:),TAUQAUX(:,:,:),TAUQX(:,:,:),TSSQX(:,:,:),
     &           TSSREF(:,:,:),TSSTX(:,:,:),WA(:,:),WB(:,:),WX(:,:)
      INTEGER IA_ERR,IOFF,IQ,IQ1,IQ2,IREF,ISPIN,IT
      LOGICAL INITIALIZE
      SAVE DROTX
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE TSSREF,GREF_I1,TSSQX,DSSQX,TAUQX,WX,WA,WB
      ALLOCATABLE DROTX,MSSTX,TSSTX,MSSQX,TAUQAUX
C
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
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
C
      NKMSLAY = NSLAY_PER_PLAY*NXM
C
C ----------------------------------------------------------------------
      IF ( INITIALIZE ) THEN
         IF ( IREL.EQ.2 ) THEN
            ALLOCATE (DROTX(NXM,NXM,NSYMMAX))
            DROTX(1:NLM,1:NLM,1:NSYMMAX) = DROT(1:NLM,1:NLM,1:NSYMMAX)
         ELSE
            ALLOCATE (DROTX(1,1,1))
         END IF
         INITIALIZE = .FALSE.
      END IF
C ----------------------------------------------------------------------
C
      ALLOCATE (TAUQX(NXM,NXM,NQ),WX(NXM,NXM))
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      IF ( IREL.EQ.2 ) THEN
         ALLOCATE (TSSTX(NLM,NLM,NT),MSSQX(NLM,NLM,NQ))
         ALLOCATE (MSSTX(NLM,NLM,NT))
      END IF
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C
      IQ1 = IQBOT
      IQ2 = IQTOP
C
C ----------------------------------------------------------------------
C     LIR and LIV: deal with left bulk first
C ----------------------------------------------------------------------
      IF ( SYSTEM_TYPE(1:3).EQ.'LIR' .AND. 
     &     SCFSTATUS.NE.'ITR-L-BULK' .AND. SCFSTATUS.NE.'START' ) THEN
         IQ1 = 1
         IQ2 = NQ
      ELSE IF ( SYSTEM_TYPE(1:3).EQ.'LIV' .AND. 
     &          SCFSTATUS.NE.'ITR-L-BULK' .AND. SCFSTATUS.NE.'START' )
     &          THEN
         IQ1 = 1
         IQ2 = IQTOP
      END IF
C ----------------------------------------------------------------------
C
C **********************************************************************
      DO ISPIN = 1,NSPIN
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         ICPACONV = 0
         CPAERRL = 1.0D+6
         ITCPA = 0
 50      CONTINUE
         ITCPA = ITCPA + 1
C
C ----------------------------------------------------------------------
         IF ( IREL.LE.2 ) THEN
C
            IOFF = NLM*(ISPIN-1)
C
            DO IQ = IQ1,IQ2
C
               IREF = IREFQ(IQ)
C
               WX(1:NLM,1:NLM) = MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,
     &                           IQ)
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
C
            IF ( IREL.EQ.2 ) THEN
               DO IQ = IQ1,IQ2
                  MSSQX(1:NLM,1:NLM,IQ)
     &               = MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,IQ)
               END DO
               DO IT = ITBOT,ITTOP
                  MSSTX(1:NLM,1:NLM,IT)
     &               = MSST(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,IT)
               END DO
            END IF
C
C ----------------------------------------------------------------------
         ELSE IF ( IREL.EQ.3 ) THEN
C
            DO IQ = IQ1,IQ2
               IREF = IREFQ(IQ)
C
               CALL CHANGEREP(NKM,NKMMAX,MSSQ(1,1,IQ),'REL>RLM',WA)
C
               CALL CMATINV(NKM,NKMMAX,WA,TSSQX(1,1,IQ))
C
               WA(1:NKM,1:NKM) = TSSQX(1:NKM,1:NKM,IQ)
               WA(1:NLM,1:NLM) = WA(1:NLM,1:NLM)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
               WA(NLM+1:NKM,NLM+1:NKM) = WA(NLM+1:NKM,NLM+1:NKM)
     &            - TSSREF(1:NLM,1:NLM,IREF)
C
               CALL CMATINV(NKM,NKMMAX,WA,DSSQX(1,1,IQ))
C
C
            END DO
C
C ----------------------------------------------------------------------
         ELSE
            STOP 'IREL'
         END IF
C
C ----------------------------------------------------------------------
C
         CALL TBKLOOP(TAUQX(1,1,IQBOT_TB),TSSQX(1,1,IQBOT_TB),
     &                DSSQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_L),
     &                DSSQX(1,1,IQBOT_R),GREF_I1,ICLU_REF_IQ(IQBOT_TB),
     &                WA,WB,NXM)
C
C ----------------------------------------------------------------------
C
         IF ( IREL.NE.2 ) THEN
C
            IF ( .NOT.ALLOCATED(TAUQAUX) )
     &           ALLOCATE (TAUQAUX(NKMMAX,NKMMAX,NQ))
C
            IF ( IREL.LT.3 ) THEN
C----------------- use auxilary array having the proper dimension NKMMAX
C
               DO IQ = IQBOT,IQTOP
                  TAUQAUX(1:NXM,1:NXM,IQ) = TAUQX(1:NXM,1:NXM,IQ)
               END DO
C
            ELSE
C--------------------- transform matrices to relativistic representation
C
               DO IQ = IQBOT,IQTOP
                  WA(1:NKM,1:NKM) = TAUQX(1:NKM,1:NKM,IQ)
C
                  CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
                  TAUQAUX(1:NKM,1:NKM,IQ) = WB(1:NKM,1:NKM)
               END DO
C
            END IF
C
            CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,
     &                     TAUQAUX(1,1,IQBOT_TB),TAUQ(1,1,IQBOT_TB),WA,
     &                     NQTB,NKMQ(IQBOT_TB),DROT,IQTBORGQTBP,
     &                     SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                     NQTB,NKMMAX)
C
         ELSE
C
            CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,
     &                     TAUQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_TB),WX,
     &                     NQTB,NLMQ(IQBOT_TB),DROTX,IQTBORGQTBP,
     &                     SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                     NQTB,NLM)
C
            TAUQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,1:NQ)
     &         = DSSQX(1:NLM,1:NLM,1:NQ)
C
         END IF
C
         IF ( NCPA.GT.0 ) THEN
C
            IF ( IREL.NE.2 ) THEN
C
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,
     &                       ICPA(IQBOT_TB),NQTB,NKMQ(IQBOT_TB),
     &                       NOQ(IQBOT_TB),ITOQ(1,IQBOT_TB),CONC,TSST,
     &                       MSST,TSSQ(1,1,IQBOT_TB),MSSQ(1,1,IQBOT_TB),
     &                       TAUQ(1,1,IQBOT_TB),NTMAX,NQTB,NKMMAX)
C
               CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,
     &                        TAUQX(1,1,IQBOT_TB),MSSQ(1,1,IQBOT_TB),WA,
     &                        NQTB,NKMQ(IQBOT_TB),DROT,IQTBORGQTBP,
     &                        SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                        NQTB,NKMMAX)
            ELSE
C
               TAUQX(1:NLM,1:NLM,1:NQ)
     &            = TAUQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,1:NQ)
C
               IF ( KMROT.NE.0 ) STOP '<TAU_TB>: KMROT <> 0'
C
               CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,
     &                       ICPA(IQBOT_TB),NQTB,NLMQ(IQBOT_TB),
     &                       NOQ(IQBOT_TB),ITOQ(1,IQBOT_TB),CONC,TSSTX,
     &                       MSSTX,TSSQX(1,1,IQBOT_TB),
     &                       MSSQX(1,1,IQBOT_TB),TAUQX(1,1,IQBOT_TB),
     &                       NTMAX,NQTB,NLMMAX)
C
               CALL SYMSUMRTR(.TRUE.,.TRUE.,.FALSE.,1D0,
     &                        DSSQX(1,1,IQBOT_TB),MSSQX(1,1,IQBOT_TB),
     &                        WX,NQTB,NLMQ(IQBOT_TB),DROTX,IQTBORGQTBP,
     &                        SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                        NQTB,NLMMAX)
C
               MSSQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,1:NQ)
     &            = MSSQX(1:NLM,1:NLM,1:NQ)
C
            END IF
C
            IF ( IPRINT.GE.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
            IF ( CPAERR.LE.CPATOL ) THEN
               ICPACONV = 1
               IF ( IPRINT.GT.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &              CPACHNG
            ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
               WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 1
            ELSE IF ( CPAERR.GT.20000*CPAERRL ) THEN
               WRITE (6,99003) ITCPA
               WRITE (6,99004) CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 2
            ELSE
               CPAERRL = CPAERR
               GOTO 50
            END IF
C
         END IF
C
      END DO
C
C **********************************************************************
C
      DEALLOCATE (TSSREF)
C
C ----------------------------------------------------------------------
99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
      END
C*==tbtssref.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TBTSSREF(ERYD,VREF,RMTREF,NL,TSSREF,NLMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculates analytically the single-site scattering matrix for a *
C   *  constant potential Vo at energy E                               *
C   *                                                                  *
C   *  input: potential radius RMTREF (R)                              *
C   *         energy ERYD (E)                                          *
C   *         potential constant value VREF (Vo)                       *
C   * output: single-site matrix TSSREF                                *
C   *                                                                  *
C   *   tmat(l) =   aj(l+1,aR)j(l,bR) - bj(l,aR)j(l+1,bR)              *
C   *             - ------------------------------------               *
C   *               j(l,bR)h(l+1,aR)a - bh(l,aR)j(l+1,bR)              *
C   *                                                                  *
C   *   a = sqrt(E),  b = sqrt(E-Vo)                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:CI
      USE MOD_TYPES,ONLY:CTL
      IMPLICIT NONE
C*--TBTSSREF356
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER NL,NLMMAX
      REAL*8 RMTREF,VREF
      COMPLEX*16 TSSREF(NLMMAX,NLMMAX)
C
C Local variables
C
      COMPLEX*16 A1,B1,BESSJW1(0:NL),BESSJW2(0:NL),BESSYW1(0:NL),
     &           BESSYW2(0:NL),CI_OV_P_FEG,HANKWS1(0:NL),HANKWS2(0:NL),
     &           P_FEG,P_TB,TLL
      REAL*8 C
      INTEGER IL,L,LM,M
      EXTERNAL BESHAN,CINIT
C
C*** End of declarations rewritten by SPAG
C
      CALL CINIT(NLMMAX*NLMMAX,TSSREF)
C
C---------------------------------------------------- calculate momentum
C
      C = CTL(1,1)
C
      CALL GET_MOMENTUM(IREL,C,ERYD,P_FEG)
      CALL GET_MOMENTUM(IREL,C,(ERYD-VREF),P_TB)
C
      IF ( DIMAG(ERYD).GE.-1.0D-15 ) THEN
         P_FEG = SQRT(ERYD)
         P_TB = SQRT(ERYD-VREF)
      ELSE
         P_FEG = -DCONJG(SQRT(DCONJG(ERYD)))
         P_TB = -DCONJG(SQRT(DCONJG(ERYD-VREF)))
      END IF
C
      A1 = P_FEG*RMTREF
      B1 = P_TB*RMTREF
      CI_OV_P_FEG = CI/P_FEG
C
      CALL BESHAN(HANKWS1,BESSJW1,BESSYW1,A1,NL)
      CALL BESHAN(HANKWS2,BESSJW2,BESSYW2,B1,NL)
C
      LM = 0
      DO IL = 1,NL
         L = IL - 1
         A1 = P_FEG*BESSJW1(L+1)*BESSJW2(L) - P_TB*BESSJW1(L)
     &        *BESSJW2(L+1)
C
         B1 = P_FEG*HANKWS1(L+1)*BESSJW2(L) - P_TB*HANKWS1(L)
     &        *BESSJW2(L+1)
C
         TLL = -CI_OV_P_FEG*A1/B1
         DO M = -L,L
            LM = LM + 1
            TSSREF(LM,LM) = TLL
         END DO
      END DO
C
      END
C*==tbgrefrs.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TBGREFRS(ERYD,P,TSSREF,GREF_I1)
C   ********************************************************************
C   *                                                                  *
C   *  supply ssite t-matrix  TSSREF of the  NREF  TB reference atoms  *
C   *  and the real space TB Green''s function matrices   GREF_I1      *
C   *  for the inequivalent   NCLU_REF  reference clusters             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1,C_AU
      USE MOD_TB,ONLY:NREF,RMTREF,VREF
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_ANGMOM,ONLY:NL,NLMMAX,NLMAX,NLM
      USE MOD_TBCLU,ONLY:NQCLU_GEO,NCLU_GEO,DSSP_RS,NIJSTAB_RS,
     &    ISSP2AB_RS,ISSP2BB_RS,IQCLUTAB_RS,JQCLUTAB_RS,NSSP4_RS,
     &    ISSP4_RS,ISDA4_RS,JSDA4_RS,ISSP5_RS,NSSP1_RS,NSSP2A_RS,
     &    NSSP2B_RS,NSSPDIR_RS,NSSPABS_RS,CLURYLM_RS,CLUIPH_RS,
     &    IND0QCLU_GEO,IND0QCLU0_RS,NKKR_RS,NKKRNR_RSMAX,NIJSTABMAX,
     &    NSSP4MAX,NLMQCLU_GEO,NCLU_REF_ICLU_GEO,ICLU_REF_ICLU_GEO,
     &    IREF_IQCLU_ICLU_REF,NCLU_REF,NRGNT,NRGNT12MAX,NRGNT123MAX,
     &    IL3RGNT,LM3RGNT,RGNT_CLU,FL1L2
      IMPLICIT NONE
C*--TBGREFRS453
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL SAMENLQ
      PARAMETER (SAMENLQ=.TRUE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      COMPLEX*16 GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           TSSREF(NLMMAX,NLMMAX,NREF)
C
C Local variables
C
      COMPLEX*16 BKJ,HLTAB(2*NLMAX-1),JLTAB(2*NLMAX-1),MG0(:,:),
     &           NLTAB(2*NLMAX-1),WG(:,:)
      INTEGER ICLU,ICLU_GEO,ICLU_REF,II,INFO,IPIV(:),IREF,ISUM,J,JJ,
     &        JQCLU,K,KK,KK0,L3MAX,LMAX,N,NLM3MAX,NLMQCLU,NQCLU
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
      ALLOCATABLE MG0,IPIV,WG
C
      ALLOCATE (WG(NKKRNR_RSMAX,NKKRNR_RSMAX))
      ALLOCATE (MG0(NKKRNR_RSMAX,NKKRNR_RSMAX),IPIV(NKKRNR_RSMAX))
C-----------------------------------------------------------------------
C
      LMAX = NL - 1
      L3MAX = 2*(NLMAX-1)
      NLM3MAX = (L3MAX+1)*(L3MAX+1)
C
C **********************************************************************
C    supply ssite t-matrix  TSSREF of the  NREF  TB reference atoms
C **********************************************************************
C
      DO IREF = 1,NREF
         CALL TBTSSREF(ERYD,VREF(IREF),RMTREF(IREF),NL,TSSREF(1,1,IREF),
     &                 NLMMAX)
      END DO
C
C **********************************************************************
C         supply the real space TB Green''s function matrices   GREF_I1
C            for the inequivalent   NCLU_REF  reference clusters
C **********************************************************************
C
      ISUM = 0
      DO ICLU_GEO = 1,NCLU_GEO
         ISUM = ISUM + NCLU_REF_ICLU_GEO(ICLU_GEO)
      END DO
      IF ( ISUM.NE.NCLU_REF ) WRITE (6,*) '<TAU_TB> ISUM <> NCLU_REF'
C
      DO ICLU_GEO = 1,NCLU_GEO
C
C-----------------------------------------------------------------------
C     get real space free electron Green''s function    MG0 = - G0
C         for inequivalent clusters concerning  ONLY  geometry
C-----------------------------------------------------------------------
C
         CALL CLUG0MAT(MG0,NKKRNR_RSMAX,NQCLU_GEO(ICLU_GEO),
     &                 NLMQCLU_GEO(1,ICLU_GEO),DSSP_RS(0,ICLU_GEO),
     &                 ISDA4_RS(1,1,ICLU_GEO),NSSP2A_RS(ICLU_GEO),
     &                 ISSP2AB_RS(1,ICLU_GEO),ISSP2BB_RS(1,ICLU_GEO),
     &                 ISSP4_RS(1,1,ICLU_GEO),ISSP5_RS(1,1,ICLU_GEO),
     &                 IQCLUTAB_RS(1,1,ICLU_GEO),JSDA4_RS(1,1,ICLU_GEO),
     &                 JQCLUTAB_RS(1,1,ICLU_GEO),NIJSTAB_RS(1,ICLU_GEO),
     &                 NIJSTABMAX,NSSP4_RS(1,ICLU_GEO),NSSP4MAX,
     &                 NSSP1_RS(ICLU_GEO),NSSP2B_RS(ICLU_GEO),
     &                 NSSPABS_RS(ICLU_GEO),NSSPDIR_RS(ICLU_GEO),NRGNT,
     &                 IL3RGNT,LM3RGNT,RGNT_CLU,LMAX,ERYD,P,C_AU,ALAT,
     &                 NLM,IND0QCLU_GEO(1,ICLU_GEO),
     &                 IND0QCLU0_RS(1,ICLU_GEO),SAMENLQ,
     &                 CLUIPH_RS(1,1,ICLU_GEO),CLURYLM_RS(1,1,ICLU_GEO),
     &                 JLTAB,NLTAB,HLTAB,FL1L2,NLMAX,NLM3MAX,NRGNT12MAX,
     &                 NRGNT123MAX)
C
C ======================================================================
         DO ICLU = 1,NCLU_REF_ICLU_GEO(ICLU_GEO)
C
            ICLU_REF = ICLU_REF_ICLU_GEO(ICLU,ICLU_GEO)
C
C-----------------------------------------------------------------------
C    multiply with TB t-matrix  WG := -(-G0) * t_REF  = G0 * t_REF
C-----------------------------------------------------------------------
C
            CALL CINIT(NKKRNR_RSMAX*NKKRNR_RSMAX,WG)
C
            NQCLU = NQCLU_GEO(ICLU_GEO)
            NLMQCLU = NLM*NQCLU
C
            JJ = 0
            DO JQCLU = 1,NQCLU
               IREF = IREF_IQCLU_ICLU_REF(JQCLU,ICLU_REF)
C
               KK0 = IND0QCLU0_RS(JQCLU,ICLU_GEO)
               DO J = 1,NLM
                  JJ = JJ + 1
C
                  KK = KK0
                  DO K = 1,NLM
                     KK = KK + 1
                     BKJ = TSSREF(K,J,IREF)
                     DO II = 1,NLMQCLU
                        WG(II,JJ) = WG(II,JJ) - MG0(II,KK)*BKJ
                     END DO
                  END DO
C
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                  WG := -1 + WK = -1 + G0 * t_REF
C-----------------------------------------------------------------------
C
            DO JJ = 1,NLMQCLU
               WG(JJ,JJ) = -C1 + WG(JJ,JJ)
            END DO
C
C-----------------------------------------------------------------------
C              solve ( 1 - t_REF * G0) * G_REF =   G0
C              or    (-1 + t_REF * G0) * G_REF = - G0
C
C             solve only for first QCLU-column with width NLM
C             i.e. get   GREF_I1  for  JQCLU=1  and  IQCLU=1,...,NQCLU
C-----------------------------------------------------------------------
C
            N = NKKR_RS(ICLU_GEO)
C
            CALL ZGETRF(N,N,WG,NKKRNR_RSMAX,IPIV,INFO)
C
            GREF_I1(1:NKKRNR_RSMAX,1:NLM,ICLU_REF)
     &         = MG0(1:NKKRNR_RSMAX,1:NLM)
C
            CALL ZGETRS('N',N,NLM,WG,NKKRNR_RSMAX,IPIV,
     &                  GREF_I1(1,1,ICLU_REF),NKKRNR_RSMAX,INFO)
C
         END DO
C ======================================================================
C
      END DO
C
      END
