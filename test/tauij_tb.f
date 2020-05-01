C*==tauij_tb.f    processed by SPAG 6.70Rc at 22:18 on 20 Dec 2016
      SUBROUTINE TAUIJ_TB(ERYD,P,MSSQ)
C   ********************************************************************
C   *                                                                  *
C   *  - call <TBGREFRS> to supply ssite t-matrix  TSSREF              *
C   *    and the real space TB Green''s function matrices   GREF_I1    *
C   *                                                                  *
C   *  - call <TAUIJ_TB_KLOOP> to perform the TB k-integration         *
C   *    and return TAUIJ                                              *
C   *                                                                  *
C   *  NOTE: this routines works for IREL=0,1,2 and 3                  *
C   *        in particular for the spin-polarised scalar relativistic  *
C   *        case many temporary arrays are needed                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_L,IQBOT_R,IQBOT_TB,IQBOT,IQTOP
      USE MOD_TB,ONLY:NREF,IREFQ,NKMSLAY,NSLAY_PER_PLAY
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM,NKM,NLMMAX,NSPIN,NXM
      USE MOD_TBCLU,ONLY:ICLU_REF_IQ,NKKRNR_RSMAX,NCLU_REF
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TAUIJ,ONLY:TAUIJ,NTAUIJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 DSSQX(:,:,:),GREF_I1(:,:,:),TSSQX(:,:,:),TSSREF(:,:,:),
     &           WA(:,:),WB(:,:),WX(:,:)
      INTEGER IA_ERR,IOFF,IQ,IREF,ISPIN,ITAUIJ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TSSREF,GREF_I1,TSSQX,DSSQX,WX,WA,WB
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
      ALLOCATE (WX(NXM,NXM))
      ALLOCATE (WA(NKMMAX,NKMMAX),WB(NKMMAX,NKMMAX))
      ALLOCATE (TSSQX(NXM,NXM,NQ),DSSQX(NXM,NXM,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '        allocate TSSREF'
C
C-----------------------------------------------------------------------
C                           initialize
C-----------------------------------------------------------------------
C
      IF ( .NOT.ALLOCATED(TAUIJ) ) THEN
         ALLOCATE (TAUIJ(NXM,NXM,NSPIN,NTAUIJ),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc: TAUIJ_STD->TAUIJ'
      END IF
C
C **********************************************************************
      DO ISPIN = 1,NSPIN
C
C ----------------------------------------------------------------------
         IF ( IREL.LE.2 ) THEN
C
            IOFF = NLM*(ISPIN-1)
C
            DO IQ = IQBOT,IQTOP
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
            END DO
C
C ----------------------------------------------------------------------
         ELSE
            STOP 'IREL'
         END IF
C
C ----------------------------------------------------------------------
C
         CALL TAUIJ_TB_KLOOP(TSSQX(1,1,IQBOT_TB),DSSQX(1,1,IQBOT_TB),
     &                       DSSQX(1,1,IQBOT_L),DSSQX(1,1,IQBOT_R),
     &                       GREF_I1,ICLU_REF_IQ(IQBOT_TB),WA,WB,ISPIN)
C
C
C
C ----------------------------------------------------------------------
C
         IF ( IREL.NE.2 ) THEN
C
            IF ( IREL.EQ.3 ) THEN
               DO ITAUIJ = 1,NTAUIJ
                  WA(1:NXM,1:NXM) = TAUIJ(1:NXM,1:NXM,ISPIN,ITAUIJ)
C
                  CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
                  TAUIJ(1:NXM,1:NXM,ISPIN,ITAUIJ) = WB(1:NKM,1:NKM)
               END DO
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
      END
C*==tauij_tb_kloop.f    processed by SPAG 6.70Rc at 22:18 on 20 Dec 2016
      SUBROUTINE TAUIJ_TB_KLOOP(TSSQX,DSSQX,DSSQX_L,DSSQX_R,GREF_I1,
     &                          ICLU_REF_IQTB,WA,WB,ISPIN)
C   ********************************************************************
C   *                                                                  *
C   *  perform the loop over k-points                                  *
C   *                                                                  *
C   *  INPUT:                                                          *
C   *          GREF_I1  G-matrix of the real space TB refernce cluster *
C   *                                                                  *
C   *  the ending X indicates matrices with:                           *
C   *  - the angular momentum indexing depends on the calculation mode *
C   *    NXM = NLM for IREL=0,1,2                                      *
C   *        = NKM for IREL=3                                          *
C   *  - the site index runs over the present TB-sites                 *
C   *    IQTB = 1,..,NQTB with                                         *
C   *                                                                  *
C   *    NQTB = NQ   dealing with a 3D-system                          *
C   *           NQ_L dealing with the left host                        *
C   *           NQ_I dealing with the interaction zone or slab         *
C   *           NQ_R dealing with the right host                       *
C   *                                                                  *
C   *  NOTE: the routine calculates the site-diagonal  TAU_Delta       *
C   *        according to the TB-scheme. After completing the k-loop   *
C   *        the result is converted to the standard   TAUQIJ          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_CONSTANTS,ONLY:C0,C1,CI2PI
      USE MOD_TAUIJ,ONLY:NTAUIJ,NTAUIJ_QTAB2_TAUIJ,IQ_QTAB1_TAUIJ,
     &    NQTAB1_TAUIJ,N5VEC_TAUIJ,N123TAUIJMAX,TAUIJ,NKTABTAUIJ,
     &    KTABTAUIJ,JQ_QTAB2_TAUIJ,NQTAB2_TAUIJ,ITAUIJ_LOOP
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NLMMAX,NKM,NLM,NXM,IXM0_QTB,NXM_QTB
      USE MOD_SITES,ONLY:NQ_L,NQ_R,NQTB,IQBOT_TB,NQCLU,QBAS
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK,NPLAY
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ISPIN
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_L(NXM,NXM,NQ_L),
     &           DSSQX_R(NXM,NXM,NQ_R),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           TSSQX(NXM,NXM,NQTB),WA(NXM,NXM),WB(NXM,NXM)
      INTEGER ICLU_REF_IQTB(NQTB)
C
C Local variables
C
      REAL*8 DDOT
      REAL*8 DQJI(3),IM_KVEC(3),KVEC(3),NORM,WK
      COMPLEX*16 EXIK1,EXIK2,EXIK3,EXIKDQ,EXIKR,EXIKRS,GLLKE(:,:),
     &           GREFKE(:,:),PWEXIK1(:),PWEXIK2(:),PWEXIK3(:)
      INTEGER I,I1,IA_ERR,IK,ILOOP,IPW,IQ,IQ0_TB,IQCLU,IQTAB,IQTB,
     &        ITAUIJ,ITAUIJ1,ITAUIJ2,ITAUIJ_DIAG(:),IU,IU0,IV,IV0,IX,
     &        IX0,J,J0,J1,JQ,JQCLU,JQTAB,JQTB,JU,JU0,JV,JV0,JX,N1,N2,N3,
     &        NI,NJ,NSUM
C
C*** End of declarations rewritten by SPAG
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE GLLKE,GREFKE
      ALLOCATABLE PWEXIK1,PWEXIK2,PWEXIK3
      ALLOCATABLE ITAUIJ_DIAG
C
      IF ( IPRINT.GE.1 ) WRITE (6,*) 
     &              '>>> TAUIJ_TB_KLOOP: integrate over BZ and get GFUN'
C----------------------------------------------------------------------
C
C
      ALLOCATE (GLLKE(NKKR_TB,NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '<TAUIJ_TB_KLOOP> allocate GLLKE'
      IF ( IREL.GT.2 ) THEN
         ALLOCATE (GREFKE(NKKRNR_TB,NKKRNR_TB),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP '<TAUIJ_TB_KLOOP> allocate GREFKE'
      ELSE IF ( NKKR_TB.NE.NKKRNR_TB ) THEN
         STOP '<TAUIJ_TB_KLOOP>:  NKKR_TB <> NKKRNR_TB'
      END IF
      ALLOCATE (PWEXIK1(-N123TAUIJMAX:N123TAUIJMAX))
      ALLOCATE (PWEXIK2(-N123TAUIJMAX:N123TAUIJMAX))
      ALLOCATE (PWEXIK3(-N123TAUIJMAX:N123TAUIJMAX))
      ALLOCATE (ITAUIJ_DIAG(NTAUIJ))
      ITAUIJ_DIAG(1:NTAUIJ) = 0
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Modify ICHECK needed in the case of band inversion
C               (all off-diagonal blocks will be calculated)
C               (it can be optimised for blocs that are actually needed)
C
      IF ( INVMOD.GE.1 .AND. INVMOD.LE.2 ) THEN
         DO N1 = 1,NPLAY
            DO N2 = 1,NPLAY
               ICHECK(N2,N1) = 1
            END DO
         END DO
      END IF
C
      TAUIJ(1:NXM,1:NXM,ISPIN,1:NTAUIJ) = C0
C
      PWEXIK1(0) = 1.0D0
      PWEXIK2(0) = 1.0D0
      PWEXIK3(0) = 1.0D0
C
      IQ0_TB = IQBOT_TB - 1
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C                                                          K-points loop
      DO IK = 1,NKTABTAUIJ
C
         WK = 0D0
C
         KVEC(1:3) = KTABTAUIJ(1:3,IK)
C
C-----------------------------------------------------------------------
C -> Fourier transformation, set KKR matrix M = [-(t)^-1 + G^r]
C-----------------------------------------------------------------------
C
         IF ( IREL.LE.2 ) THEN
            CALL TBGREFKE(KVEC,IM_KVEC,NQTB,GREF_I1,GLLKE,ICLU_REF_IQTB,
     &                    NKKRNR_RSMAX)
C
         ELSE
C
            CALL TBGREFKE(KVEC,IM_KVEC,NQTB,GREF_I1,GREFKE,
     &                    ICLU_REF_IQTB,NKKRNR_RSMAX)
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
     &        VACFLAG,FACTL,NQ_L,NQ_R,NXM)
C
C-----------------------------------------------------------------------
C -> Construct the matrix M=[-(t-t_ref)^-1 + G_ref] and store it
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
C --> Perform the inversion of matrix M
C     the output is the scattering path operator -TAU_DELTA(k)
C-----------------------------------------------------------------------
C
C NOTE: DSSQX dummy argument - not used for INVMOD <> 4
C
         IF ( INVMOD.EQ.4 ) STOP 'INVMOD = 4 '
C
         CALL TBINVERSION(GLLKE,INVMOD,ICHECK,NXM,DSSQX,WK)
C
C -------------------------- table to calculate   exp( i ->k * ->R[ij] )
C
         EXIK1 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,1),1))
         EXIK2 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,2),1))
         EXIK3 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,3),1))
C
         DO IPW = 1,N123TAUIJMAX
            PWEXIK1(IPW) = PWEXIK1(IPW-1)*EXIK1
            PWEXIK1(-IPW) = 1.0D0/PWEXIK1(IPW)
C
            PWEXIK2(IPW) = PWEXIK2(IPW-1)*EXIK2
            PWEXIK2(-IPW) = 1.0D0/PWEXIK2(IPW)
C
            PWEXIK3(IPW) = PWEXIK3(IPW-1)*EXIK3
            PWEXIK3(-IPW) = 1.0D0/PWEXIK3(IPW)
         END DO
C
C-----------------------------------------------------------------------
C                     get the off site-diagonal blocks
C-----------------------------------------------------------------------
C
         ITAUIJ2 = 0
         DO IQTAB = 1,NQTAB1_TAUIJ
            IQTB = IQ_QTAB1_TAUIJ(IQTAB) - IQ0_TB
            IQ = IQTB
            NI = NXM_QTB(IQTB)
            I1 = IXM0_QTB(IQTB) + 1
C
            DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
               JQTB = JQ_QTAB2_TAUIJ(IQTAB,JQTAB) - IQ0_TB
               JQ = JQTB
C
               NJ = NXM_QTB(JQTB)
               J0 = IXM0_QTB(JQTB)
               ITAUIJ1 = ITAUIJ2 + 1
               ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
               DQJI(1:3) = QBAS(1:3,JQ) - QBAS(1:3,IQ)
               EXIKDQ = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,DQJI(1),1)
     &                  )
C
               DO ITAUIJ = ITAUIJ1,ITAUIJ2
                  N1 = N5VEC_TAUIJ(1,ITAUIJ)
                  N2 = N5VEC_TAUIJ(2,ITAUIJ)
                  N3 = N5VEC_TAUIJ(3,ITAUIJ)
C
                  EXIKRS = PWEXIK1(N1)*PWEXIK2(N2)*PWEXIK3(N3)
C
                  EXIKR = EXIKRS*EXIKDQ
C
C
                  DO J = 1,NJ
                     J1 = J0 + J
                     CALL ZAXPY(NI,EXIKR,GLLKE(I1,J1),1,
     &                          TAUIJ(1,J,ISPIN,ITAUIJ),1)
                  END DO
C
               END DO
            END DO
         END DO
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
      NORM = -1D0/DBLE(NKTABTAUIJ)
      TAUIJ(1:NXM,1:NXM,ISPIN,1:NTAUIJ)
     &   = NORM*TAUIJ(1:NXM,1:NXM,ISPIN,1:NTAUIJ)
C
      IF ( NQCLU.EQ.0 ) THEN
         ITAUIJ2 = 0
         DO IQTAB = 1,NQTAB1_TAUIJ
            IQ = IQ_QTAB1_TAUIJ(IQTAB)
            DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
               JQ = JQ_QTAB2_TAUIJ(IQTAB,JQTAB)
               ITAUIJ1 = ITAUIJ2 + 1
               ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
               DO ITAUIJ = ITAUIJ1,ITAUIJ2
                  N1 = N5VEC_TAUIJ(1,ITAUIJ)
                  N2 = N5VEC_TAUIJ(2,ITAUIJ)
                  N3 = N5VEC_TAUIJ(3,ITAUIJ)
                  NSUM = ABS(N1) + ABS(N2) + ABS(N3)
                  IF ( IQ.EQ.JQ .AND. NSUM.EQ.0 ) ITAUIJ_DIAG(ITAUIJ)
     &                 = 1
               END DO
            END DO
         END DO
      ELSE
         ILOOP = 0
         DO JQCLU = 1,NQCLU
            DO IQCLU = 1,NQCLU
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               IF ( IQCLU.EQ.JQCLU ) ITAUIJ_DIAG(ITAUIJ) = 1
            END DO
         END DO
      END IF
C
      ITAUIJ2 = 0
      DO IQTAB = 1,NQTAB1_TAUIJ
         IQTB = IQ_QTAB1_TAUIJ(IQTAB) - IQ0_TB
         IQ = IQTB
         NI = NXM_QTB(IQTB)
         I1 = IXM0_QTB(IQTB) + 1
         DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
            JQTB = JQ_QTAB2_TAUIJ(IQTAB,JQTAB) - IQ0_TB
            JQ = JQTB
            NJ = NXM_QTB(JQTB)
            J0 = IXM0_QTB(JQTB)
            ITAUIJ1 = ITAUIJ2 + 1
            ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
            DO ITAUIJ = ITAUIJ1,ITAUIJ2
C
C------------- G = 1/(delta t) * TAU_Delta *1/(delta t) - 1/(delta t)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TAUIJ(1,1,ISPIN,ITAUIJ)
     &                    ,NXM,DSSQX(1,1,JQTB),NXM,C0,WB,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WB,
     &                    NXM,C0,WA,NXM)
C
               IF ( ITAUIJ_DIAG(ITAUIJ).EQ.1 ) WA(1:NXM,1:NXM)
     &              = WA(1:NXM,1:NXM) - DSSQX(1:NXM,1:NXM,IQTB)
C
C--------------------------------------------------- TAU = t * G * t + t
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,TSSQX(1,1,JQTB),
     &                    NXM,C0,WB,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TSSQX(1,1,IQTB),NXM,WB,
     &                    NXM,C0,TAUIJ(1:NXM,1:NXM,ISPIN,ITAUIJ),NXM)
C
               IF ( ITAUIJ_DIAG(ITAUIJ).EQ.1 )
     &              TAUIJ(1:NXM,1:NXM,ISPIN,ITAUIJ)
     &              = TAUIJ(1:NXM,1:NXM,ISPIN,ITAUIJ)
     &              + TSSQX(1:NXM,1:NXM,IQTB)
C
            END DO
         END DO
      END DO
C
      END
