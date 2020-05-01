C*==negfkloopsdrv_tb.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE NEGFKLOOPSDRV_TB(IFIL_ZPLS,IFIL_ZMIN,ERYD,P,MSSQ,TAUQ,
     &                            GLESQ,IQBOT_LBAR,IQTOP_LBAR,
     &                            IQBOT_RBAR,IQTOP_RBAR,BAR_FLAG_Q,
     &                            LBAR_FLAG_Q,RBAR_FLAG_Q,MEZZ,MEZZMM,
     &                            MEZZMP,CCZ,CHECK_ZPM,LMAT3,STT,
     &                            TRANSMISSION)
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,RREL,NLMMAX,NSPIN,NXM,NMEMAX,NL
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_L,IQBOT_R,IQBOT_TB,IQBOT,IQTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C*--NEGFKLOOPSDRV_TB23
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM,STT,TRANSMISSION
      COMPLEX*16 ERYD,P
      INTEGER IFIL_ZMIN,IFIL_ZPLS,IQBOT_LBAR,IQBOT_RBAR,IQTOP_LBAR,
     &        IQTOP_RBAR
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 GLESQ(NKMMAX,NKMMAX,NQMAX),LMAT3(NKM,NKM),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 DELTA(:,:,:),DSSQX(:,:,:),GLESQX(:,:,:),GREF_I1(:,:,:),
     &           LMAT(:,:),LMAT2(NKM,NKM),MEZZPM(:,:,:,:),TAUQX(:,:,:),
     &           TSSQX(:,:,:),TSSREF(:,:,:),USSQX(:,:,:),VSSQX(:,:,:),
     &           WA(:,:),WB(:,:),WX(:,:)
      INTEGER I,IA_ERR,IL,IOFF,IQ,IREF,ISPIN,M,N
      LOGICAL INITIALIZE
      REAL*8 M1PWL
      SAVE LMAT
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE DELTA,DSSQX,GLESQX,GREF_I1,LMAT,MEZZPM
      ALLOCATABLE TAUQX,TSSQX,TSSREF,USSQX,VSSQX,WA,WB,WX
C
      ALLOCATE (DELTA(NLMMAX,NLMMAX,NTMAX))
      ALLOCATE (GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF))
      ALLOCATE (MEZZPM(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (TSSREF(NLMMAX,NLMMAX,NREF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
C
      DATA INITIALIZE/.TRUE./
C
      IF ( NXM.NE.NKMMAX ) STOP '<NEGFKLOOPSDRV_TB_MO>:  NXM <> NKMMAX '
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
C
         LMAT3(:,:) = LMAT(:,:)
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
      WRITE (6,99001) 'calling TBGREFRS'
C
      CALL TBGREFRS(ERYD,P,TSSREF,GREF_I1)
C
C=======================================================================
C             calculate   Delta t(1) and (2) for SA approach
C             get MEZZPM for MO approach
C=======================================================================
      WRITE (6,99001) 'calling NEGF_TMAT_MO'
C
      CALL NEGF_TMAT(IFIL_ZPLS,IFIL_ZMIN,IQBOT_LBAR,IQTOP_RBAR,MEZZ,
     &               MEZZMM,MEZZMP,MEZZPM,STT,DELTA)
C
C **********************************************************************
C
      NKMSLAY = NSLAY_PER_PLAY*NXM
C
C ----------------------------------------------------------------------
C
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      ALLOCATE (GLESQX(NXM,NXM,NQ))
      ALLOCATE (TAUQX(NXM,NXM,NQ),WX(NXM,NXM))
      ALLOCATE (USSQX(NXM,NXM,NQ),VSSQX(NXM,NXM,NQ))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ))
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
               CALL CHANGEREP(NKM,NKMMAX,LMAT,'REL>RLM',LMAT2)
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
         WRITE (6,99001) 'calling NEGFKLOOPS_MO'
C
         CALL NEGFKLOOPS_TB(GLESQX(1,1,IQBOT_TB),TAUQX(1,1,IQBOT_TB),
     &                      TSSQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_TB),
     &                      DSSQX(1,1,IQBOT_L),DSSQX(1,1,IQBOT_R),
     &                      GREF_I1,ICLU_REF_IQ(IQBOT_TB),WA,WB,LMAT,
     &                      GLESQ,TAUQ,IQBOT_LBAR,IQTOP_LBAR,IQBOT_RBAR,
     &                      IQTOP_RBAR,LBAR_FLAG_Q,RBAR_FLAG_Q,
     &                      BAR_FLAG_Q,CCZ,CHECK_ZPM,MEZZ,MEZZMM,MEZZMP,
     &                      MEZZPM,TRANSMISSION,STT,DELTA)
C
      END DO
C **********************************************************************
C
      DEALLOCATE (TSSREF)
C
99001 FORMAT (/,79('*'),/,10X,A,/,79('*'),/)
      END
