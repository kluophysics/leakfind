C*==tau_tb_1k.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAU_TB_1K(ERYD,P,MSSQ,TAUQ,CALCGREF,KTABIN)
C
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   *  - call <TBKLOOP> to perform the appropriate TB                  *
C   *    and return TAU(k) for one k-point                             *
C   *                                                                  *
C   ********************************************************************
C
C
C
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_L,IQBOT_R,IQBOT_TB,IQBOT,IQTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,NLMMAX,NSPIN,NXM
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB
      IMPLICIT NONE
C*--TAU_TB_1K24
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CALCGREF
      COMPLEX*16 ERYD,P
      REAL*8 KTABIN(1:3)
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 DSSQX(:,:,:),GREF_I1(:,:,:),TAUQX(:,:,:),TSSQX(:,:,:),
     &           TSSREF(:,:,:),WA(:,:),WB(:,:),WX(:,:)
      INTEGER IA_ERR,IOFF,IQ,IQ1,IQ2,IREF,ISPIN,NKTABSAV
      REAL*8 KTABSAV(:,:),WKTABSAV(:)
      SAVE GREF_I1,TSSREF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TSSREF,GREF_I1,TSSQX,DSSQX,TAUQX,WX,WA,WB
      ALLOCATABLE KTABSAV,WKTABSAV
C
      IF ( .NOT.ALLOCATED(GREF_I1) ) THEN
         ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
         ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
      END IF
C
C
      ALLOCATE (KTABSAV(1:3,1:NKTAB),WKTABSAV(1:NKTAB))
C
      NKTABSAV = NKTAB
      KTABSAV(1:3,1:NKTAB) = KTAB(1:3,1:NKTAB)
      WKTABSAV(1:NKTAB) = WKTAB(1:NKTAB)
C
      NKTAB = 1
      KTAB(1:3,1) = KTABIN(1:3)
      WKTAB(1) = 1.0D0
C
C **********************************************************************
C    supply ssite t-matrix  TSSREF of the  NREF  TB reference atoms
C **********************************************************************
C         supply the real space TB Green''s function matrices   GREF_I1
C            for the inequivalent   NCLU_REF  reference clusters
C **********************************************************************
C
      IF ( CALCGREF ) THEN
C
         CALL TBGREFRS(ERYD,P,TSSREF,GREF_I1)
C
         RETURN
C
      END IF
C
C **********************************************************************
C
      NKMSLAY = NSLAY_PER_PLAY*NXM
C
C ----------------------------------------------------------------------
      ALLOCATE (TAUQX(NXM,NXM,NQ),WX(NXM,NXM))
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
C
      IQ1 = IQBOT
      IQ2 = IQTOP
      IF ( SYSTEM_TYPE(1:3).EQ.'LIR' ) THEN
         IQ1 = 1
         IQ2 = NQ
      ELSE IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
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
C ----------------------------------------------------------------------
         IF ( IREL.LE.2 ) THEN
C
            IOFF = NLM*(ISPIN-1)
C
            DO IQ = IQ1,IQ2
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
            IF ( IREL.EQ.3 ) THEN
               DO IQ = IQBOT,IQTOP
                  WA(1:NKM,1:NKM) = TAUQX(1:NKM,1:NKM,IQ)
C
                  CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
                  TAUQ(1:NKM,1:NKM,IQ) = WB(1:NKM,1:NKM)
               END DO
            END IF
            IF ( IREL.LE.1 ) TAUQ(1:NLM,1:NLM,1:NQ)
     &           = TAUQX(1:NLM,1:NLM,1:NQ)
         ELSE
C
            TAUQ(IOFF+1:IOFF+NLM,IOFF+1:IOFF+NLM,1:NQ)
     &         = TAUQX(1:NLM,1:NLM,1:NQ)
C
         END IF
C
      END DO
C
      NKTAB = NKTABSAV
      KTAB(1:3,1:NKTAB) = KTABSAV(1:3,1:NKTAB)
      WKTAB(1:NKTAB) = WKTABSAV(1:NKTAB)
C
C **********************************************************************
C
      END
