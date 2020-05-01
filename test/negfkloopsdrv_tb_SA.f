C*==negfkloopsdrv_tb_sa.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NEGFKLOOPSDRV_TB_SA(IFIL_ZPLS,IFIL_ZMIN,IFILBAR_ZPLS,
     &                               IFILBAR_ZMIN,ERYD,P,MSSQ,TAUQ,
     &                               TAUQ_BAR,GLESQ,IQBOT_LBAR,
     &                               IQTOP_LBAR,IQBOT_RBAR,IQTOP_RBAR,
     &                               BAR_FLAG_Q,LBAR_FLAG_Q,RBAR_FLAG_Q,
     &                               TSS_LBAR,TSS_RBAR,DELT1_LBAR,
     &                               DELT1_RBAR,DELT2_LBAR,DELT2_RBAR,
     &                               VBAR,TSST,MEZZ,MEZZ_BAR,MEZZMP,CCZ,
     &                               CHECK_ZPM,LMAT3)
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_L,IQBOT_R,IQBOT_TB,IQBOT,IQTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,RREL,NLMMAX,NSPIN,NXM,NMEMAX,
     &    NL,WKM1,WKM2
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_CONSTANTS,ONLY:C0,C1,CI
      IMPLICIT NONE
C*--NEGFKLOOPSDRV_TB_SA27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER IZPLS,IZMIN
      PARAMETER (IZPLS=1,IZMIN=2)
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM
      COMPLEX*16 ERYD,P
      INTEGER IFILBAR_ZMIN,IFILBAR_ZPLS,IFIL_ZMIN,IFIL_ZPLS,IQBOT_LBAR,
     &        IQBOT_RBAR,IQTOP_LBAR,IQTOP_RBAR
      REAL*8 VBAR
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 DELT1_LBAR(NKMMAX,NKMMAX,2,2,IQBOT_LBAR:IQTOP_LBAR),
     &           DELT1_RBAR(NKMMAX,NKMMAX,2,2,IQBOT_RBAR:IQTOP_RBAR),
     &           DELT2_LBAR(NKMMAX,NKMMAX,2,2,IQBOT_LBAR:IQTOP_LBAR),
     &           DELT2_RBAR(NKMMAX,NKMMAX,2,2,IQBOT_RBAR:IQTOP_RBAR),
     &           GLESQ(NKMMAX,NKMMAX,NQMAX),LMAT3(NKM,NKM),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ_BAR(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQ_BAR(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
     &           ,TSS_LBAR(NKMMAX,NKMMAX,IQBOT_LBAR:IQTOP_LBAR),
     &           TSS_RBAR(NKMMAX,NKMMAX,IQBOT_RBAR:IQTOP_RBAR)
C
C Local variables
C
      COMPLEX*16 DSSQX(:,:,:),DSSQX_BAR(:,:,:),GLESQX(:,:,:),
     &           GREF_I1(:,:,:),GSSQX_LBAR(:,:,:),GSSQX_RBAR(:,:,:),
     &           LMAT(:,:),LMAT2(NKM,NKM),TAUQX(:,:,:),TAUQX_BAR(:,:,:),
     &           TSSQX(:,:,:),TSSQX_BAR(:,:,:),TSSREF(:,:,:),
     &           USSQX(:,:,:),USSQX_BAR(:,:,:),USSQX_LBAR(:,:,:,:),
     &           USSQX_RBAR(:,:,:,:),VSSQX(:,:,:),VSSQX_BAR(:,:,:),
     &           VSSQX_LBAR(:,:,:,:),VSSQX_RBAR(:,:,:,:),WA(:,:),WB(:,:)
     &           ,WX(:,:)
      INTEGER I,IA_ERR,IL,IOFF,IQ,IREF,ISPIN,IZ,JZ,M,N
      LOGICAL INITIALIZE
      REAL*8 M1PWL
      SAVE LMAT
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE TSSREF,GREF_I1,USSQX,VSSQX,TSSQX,DSSQX,TAUQX,WX,WA,WB
      ALLOCATABLE USSQX_BAR,VSSQX_BAR,TSSQX_BAR,DSSQX_BAR,TAUQX_BAR
      ALLOCATABLE USSQX_LBAR,VSSQX_LBAR,USSQX_RBAR,VSSQX_RBAR,GLESQX
      ALLOCATABLE GSSQX_LBAR,GSSQX_RBAR,LMAT
C
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
C
      DATA INITIALIZE/.TRUE./
C
      IF ( NXM.NE.NKMMAX ) STOP '<NEGFKLOOPSDRV_TB>:  NXM <> NKMMAX '
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
         ALLOCATE (LMAT(NXM,NXM))
C
         LMAT(:,:) = C0
         LMAT2(:,:) = C0
C
         IF ( IREL.LE.2 ) THEN
            N = 1
         ELSE
            N = 2
         END IF
         I = 0
         M1PWL = -1D0
         DO IL = 1,NL
            M1PWL = -M1PWL
            DO M = 1,N*(2*IL-1)
               I = I + 1
               LMAT(I,I) = M1PWL
            END DO
         END DO
Csw
         LMAT3(:,:) = LMAT(:,:)
Csw
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
C=======================================================================
C             calculate   Delta t(1) and (2)
C=======================================================================
C
      CALL NEGF_TMAT_SA(IFIL_ZPLS,IFIL_ZMIN,IFILBAR_ZPLS,IFILBAR_ZMIN,
     &                  DELT1_LBAR,DELT2_LBAR,IQBOT_LBAR,IQTOP_LBAR,
     &                  DELT1_RBAR,DELT2_RBAR,IQBOT_RBAR,IQTOP_RBAR,
     &                  BAR_FLAG_Q,LBAR_FLAG_Q,RBAR_FLAG_Q,VBAR,
     &                  TSS_LBAR,TSS_RBAR,TSST,MEZZ,MEZZ_BAR,MEZZMP)
C
C **********************************************************************
C
      NKMSLAY = NSLAY_PER_PLAY*NXM
C
C ----------------------------------------------------------------------
C
      ALLOCATE (USSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR,2))
      ALLOCATE (USSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR,2))
      ALLOCATE (VSSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR,2))
      ALLOCATE (VSSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR,2))
      ALLOCATE (GSSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR))
      ALLOCATE (GSSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR))
C
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      ALLOCATE (GLESQX(NXM,NXM,NQ))
      ALLOCATE (TAUQX(NXM,NXM,NQ),WX(NXM,NXM))
      ALLOCATE (USSQX(NXM,NXM,NQ),VSSQX(NXM,NXM,NQ))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ))
      ALLOCATE (TAUQX_BAR(NXM,NXM,NQ))
      ALLOCATE (USSQX_BAR(NXM,NXM,NQ),VSSQX_BAR(NXM,NXM,NQ))
      ALLOCATE (TSSQX_BAR(NXM,NXM,NQ),DSSQX_BAR(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
C
C **********************************************************************
      DO ISPIN = 1,NSPIN
C
C=======================================================================
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
C ----------------------------------------------------------------------
C                          barrier sites
C ----------------------------------------------------------------------
C
               IF ( LBAR_FLAG_Q(IQ) ) THEN
                  TSSQX_BAR(1:NLM,1:NLM,IQ) = TSS_LBAR(1:NLM,1:NLM,IQ)
               ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
                  TSSQX_BAR(1:NLM,1:NLM,IQ) = TSS_RBAR(1:NLM,1:NLM,IQ)
               ELSE
                  TSSQX_BAR(1:NLM,1:NLM,IQ) = TSSQX(1:NLM,1:NLM,IQ)
               END IF
C
               WX(1:NLM,1:NLM) = TSSQX_BAR(1:NLM,1:NLM,IQ)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
C
               CALL CMATINV(NLM,NLM,WX,DSSQX_BAR(1,1,IQ))
C
C ----------------------------------------------------------------------
C
            END DO
C=======================================================================
         ELSE IF ( IREL.EQ.3 ) THEN
C
            M = NKMMAX
C
            DO IQ = IQBOT,IQTOP
               IREF = IREFQ(IQ)
C
Csw
               CALL CHANGEREP(NKM,NKMMAX,LMAT,'REL>RLM',LMAT2)
Csw
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
C ----------------------------------------------------------------------
C                          barrier sites
C ----------------------------------------------------------------------
C
               IF ( LBAR_FLAG_Q(IQ) ) THEN
                  WA(1:NKM,1:NKM) = TSS_LBAR(1:NKM,1:NKM,IQ)
                  CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                           TSSQX_BAR(1,1,IQ))
               ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
                  WA(1:NKM,1:NKM) = TSS_RBAR(1:NKM,1:NKM,IQ)
                  CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                           TSSQX_BAR(1,1,IQ))
               ELSE
                  TSSQX_BAR(1:NKM,1:NKM,IQ) = TSSQX(1:NKM,1:NKM,IQ)
               END IF
C
               WA(1:NKM,1:NKM) = TSSQX_BAR(1:NKM,1:NKM,IQ)
               WA(1:NLM,1:NLM) = WA(1:NLM,1:NLM)
     &                           - TSSREF(1:NLM,1:NLM,IREF)
               WA(NLM+1:NKM,NLM+1:NKM) = WA(NLM+1:NKM,NLM+1:NKM)
     &            - TSSREF(1:NLM,1:NLM,IREF)
C
               CALL CMATINV(NKM,NKMMAX,WA,DSSQX_BAR(1,1,IQ))
C
C - LEFT  transformation matrix U = RREL^+ * t * 1/(delta t)
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,RREL,M,
     &                    TSSQX_BAR(1,1,IQ),M,C0,WA,M)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,
     &                    DSSQX_BAR(1,1,IQ),NXM,C0,USSQX_BAR(1,1,IQ),
     &                    NXM)
C
C - RIGHT transformation matrix V = 1/(delta t) * t * RREL
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQ),NXM,
     &                    TSSQX_BAR(1,1,IQ),NXM,C0,WA,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,M,RREL,M,C0,
     &                    VSSQX_BAR(1,1,IQ),M)
C
C
C ----------------------------------------------------------------------
C     transformation matrices U(+), U(-), V(+), V(-) for barrier sites
C             and single site contribution to  gamma(+,-)
C ----------------------------------------------------------------------
C
               IF ( LBAR_FLAG_Q(IQ) ) THEN
C
                  DO IZ = 1,2
                     DO JZ = 1,2
C
                        WA(1:NKM,1:NKM)
     &                     = DELT1_LBAR(1:NKM,1:NKM,IZ,JZ,IQ)
                        CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                                 DELT1_LBAR(1,1,IZ,JZ,IQ))
C
                        WA(1:NKM,1:NKM)
     &                     = DELT2_LBAR(1:NKM,1:NKM,IZ,JZ,IQ)
                        CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                                 DELT2_LBAR(1,1,IZ,JZ,IQ))
C
                     END DO
                  END DO
C
C------------------------------------------------------------------ U(+)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT1_LBAR(1,1,IZPLS,IZPLS,IQ),NXM,
     &                       DSSQX_BAR(1,1,IQ),NXM,C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,CI,DSSQX(1,1,IQ),NXM,
     &                       WA,NXM,C0,USSQX_LBAR(1,1,IQ,IZPLS),NXM)
C
C------------------------------------------------------------------ U(-)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX_BAR(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT1_LBAR(1,1,IZPLS,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,CI,DSSQX(1,1,IQ),NXM,
     &                       WA,NXM,C0,USSQX_LBAR(1,1,IQ,IZMIN),NXM)
C
C------------------------------------------------------------------ V(+)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT2_LBAR(1,1,IZPLS,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQ),
     &                       NXM,WA,NXM,C0,VSSQX_LBAR(1,1,IQ,IZPLS),NXM)
C
C------------------------------------------------------------------ V(-)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT2_LBAR(1,1,IZMIN,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WA,NXM,C0,
     &                       WB,NXM)
C
                  CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQ),
     &                       NXM,WB,NXM,C0,VSSQX_LBAR(1,1,IQ,IZMIN),NXM)
C
C-------------------------------------------------------------- gss(+,-)
C
                  GSSQX_LBAR(1:NXM,1:NXM,IQ)
     &               = CI*(DELT1_LBAR(1:NXM,1:NXM,IZPLS,IZMIN,IQ)
     &               -DELT2_LBAR(1:NXM,1:NXM,IZPLS,IZMIN,IQ))
C
Csw                  IF ( IPRINT.GE.1 ) THEN
                  WRITE (41,*) 'IQ =',IQ
                  CALL CMATSTRUCT('GSSQX_LBAR',GSSQX_LBAR(1,1,IQ),18,18,
     &                            2,2,1,1.D-8,41)
Csw                  END IF
C
C--------------------- multiply from left x(+) coming from outer G^+ (+)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQ),NXM,
     &                       GSSQX_LBAR(1,1,IQ),NXM,C0,WA,NXM)
C
C-- and from right with x(-) = L x(+)^dagger L coming from outer G^+ (-)
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WKM1,NXM)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,
     &                       C0,WKM2,NXM)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,WKM2,NXM,C0,
     &                       GSSQX_LBAR(1,1,IQ),NXM)
C
               ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
C
                  DO IZ = 1,2
                     DO JZ = 1,2
C
                        WA(1:NKM,1:NKM)
     &                     = DELT1_RBAR(1:NKM,1:NKM,IZ,JZ,IQ)
                        CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                                 DELT1_RBAR(1,1,IZ,JZ,IQ))
C
                        WA(1:NKM,1:NKM)
     &                     = DELT2_RBAR(1:NKM,1:NKM,IZ,JZ,IQ)
                        CALL CHANGEREP(NKM,NKMMAX,WA,'REL>RLM',
     &                                 DELT2_RBAR(1,1,IZ,JZ,IQ))
C
                     END DO
                  END DO
C
C------------------------------------------------------------------ U(+)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT1_RBAR(1,1,IZPLS,IZPLS,IQ),NXM,
     &                       DSSQX_BAR(1,1,IQ),NXM,C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,CI,DSSQX(1,1,IQ),NXM,
     &                       WA,NXM,C0,USSQX_RBAR(1,1,IQ,IZPLS),NXM)
C
C------------------------------------------------------------------ U(-)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX_BAR(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT1_RBAR(1,1,IZPLS,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,CI,DSSQX(1,1,IQ),NXM,
     &                       WA,NXM,C0,USSQX_RBAR(1,1,IQ,IZMIN),NXM)
C
C------------------------------------------------------------------ V(+)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT2_RBAR(1,1,IZPLS,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQ),
     &                       NXM,WA,NXM,C0,VSSQX_RBAR(1,1,IQ,IZPLS),NXM)
C
C------------------------------------------------------------------ V(-)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WB,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,
     &                       DELT2_RBAR(1,1,IZMIN,IZMIN,IQ),NXM,WB,NXM,
     &                       C0,WA,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WA,NXM,C0,
     &                       WB,NXM)
C
                  CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQ),
     &                       NXM,WB,NXM,C0,VSSQX_RBAR(1,1,IQ,IZMIN),NXM)
C
C-------------------------------------------------------------- gss(+,-)
C
                  GSSQX_RBAR(1:NXM,1:NXM,IQ)
     &               = CI*(DELT1_RBAR(1:NXM,1:NXM,IZPLS,IZMIN,IQ)
     &               -DELT2_RBAR(1:NXM,1:NXM,IZPLS,IZMIN,IQ))
C
Csw                  IF ( IPRINT.GE.1 ) THEN
                  WRITE (42,*) 'IQ =',IQ
                  CALL CMATSTRUCT('GSSQX_RBAR',GSSQX_RBAR(1,1,IQ),18,18,
     &                            2,2,1,1.D-8,42)
Csw                  END IF
C
C--------------------- multiply from left x(+) coming from outer G^+ (+)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQ),NXM,
     &                       GSSQX_RBAR(1,1,IQ),NXM,C0,WA,NXM)
C
C-- and from right with x(-) = L x(+)^dagger L coming from outer G^+ (-)
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,IQ),NXM,C0,WKM1,NXM)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,
     &                       C0,WKM2,NXM)
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,WKM2,NXM,C0,
     &                       GSSQX_RBAR(1,1,IQ),NXM)
C
               END IF
C
C ----------------------------------------------------------------------
C
            END DO
C
C ----------------------------------------------------------------------
         ELSE
            STOP 'IREL'
         END IF
C
C ----------------------------------------------------------------------
         WRITE (6,99001) 'calling NEGFKLOOPS_SA'
C
         CALL NEGFKLOOPS_TB_SA(GLESQX(1,1,IQBOT_TB),TAUQX(1,1,IQBOT_TB),
     &                         TSSQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_TB),
     &                         DSSQX(1,1,IQBOT_L),DSSQX(1,1,IQBOT_R),
     &                         TAUQX_BAR(1,1,IQBOT_TB),
     &                         TSSQX_BAR(1,1,IQBOT_TB),
     &                         DSSQX_BAR(1,1,IQBOT_TB),USSQX_LBAR,
     &                         VSSQX_LBAR,GSSQX_LBAR,USSQX_RBAR,
     &                         VSSQX_RBAR,GSSQX_RBAR,GREF_I1,
     &                         ICLU_REF_IQ(IQBOT_TB),WA,WB,LMAT,GLESQ,
     &                         TAUQ,TAUQ_BAR,IQBOT_LBAR,IQTOP_LBAR,
     &                         IQBOT_RBAR,IQTOP_RBAR,LBAR_FLAG_Q,
     &                         RBAR_FLAG_Q,BAR_FLAG_Q,CCZ,CHECK_ZPM)
C
      END DO
C **********************************************************************
C
      DEALLOCATE (TSSREF)
C
99001 FORMAT (/,79('*'),/,10X,A,/,79('*'),/)
      END
