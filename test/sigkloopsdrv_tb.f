C*==sigkloopsdrv_tb.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGKLOOPSDRV_TB(ERYD,P,MSSQ,TAUQBZ)
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NZ12MAX
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_L,IQBOT_R,IQBOT_TB,IQBOT,IQTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,RREL,NLMMAX,NSPIN,NXM
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--SIGKLOOPSDRV_TB18
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQBZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX)
C
C Local variables
C
      COMPLEX*16 DSSQX(:,:,:),GREF_I1(:,:,:),TAUQX(:,:,:),TSSQX(:,:,:),
     &           TSSREF(:,:,:),USSQX(:,:,:),VSSQX(:,:,:),WA(:,:),WB(:,:)
     &           ,WX(:,:)
      INTEGER IA_ERR,IOFF,IQ,IREF,ISPIN,M
      LOGICAL INITIALIZE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TSSREF,GREF_I1,USSQX,VSSQX,TSSQX,DSSQX,TAUQX,WX,WA,WB
C
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
C
      DATA INITIALIZE/.TRUE./
C
      IF ( NXM.NE.NKMMAX ) STOP '<SIGKLOOPSDRV_TB>:  NXM <> NKMMAX '
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C ----------------------------------------------------------------------
C          prepare the k-mesh without symmetry restrictions
C ----------------------------------------------------------------------
C
         CALL INIT_MOD_TAUIJ_KMESH
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
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
C
      ALLOCATE (TAUQX(NXM,NXM,NQ),WX(NXM,NXM))
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      ALLOCATE (USSQX(NXM,NXM,NQ),VSSQX(NXM,NXM,NQ))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
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
C ----------------------------------------------------------------------
         ELSE IF ( IREL.EQ.3 ) THEN
C
            M = NKMMAX
C
            DO IQ = IQBOT,IQTOP
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
C - LEFT  transformation matrix U = RREL^+ * t * 1/(delta t)
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,RREL,M,TSSQX(1,1,IQ),M,
     &                    C0,WA,M)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,DSSQX(1,1,IQ),
     &                    NXM,C0,USSQX(1,1,IQ),NXM)
C
C - RIGHT transformation matrix V = 1/(delta t) * t * RREL
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQ),NXM,
     &                    TSSQX(1,1,IQ),NXM,C0,WA,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,M,RREL,M,C0,
     &                    VSSQX(1,1,IQ),M)
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
         CALL SIGKLOOPS_TB(TAUQX(1,1,IQBOT_TB),TSSQX(1,1,IQBOT_TB),
     &                     DSSQX(1,1,IQBOT_TB),USSQX(1,1,IQBOT_TB),
     &                     VSSQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_L),
     &                     DSSQX(1,1,IQBOT_R),GREF_I1,
     &                     ICLU_REF_IQ(IQBOT_TB),WA,WB,TAUQBZ)
C
      END DO
C **********************************************************************
C
      DEALLOCATE (TSSREF)
C
      END
