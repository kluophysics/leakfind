C*==clug0mat.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUG0MAT(MG0MAT,NKKR,NQCLU,NLMQCLU,DSSP,ISDA4,NSSP2A,
     &                    ISSP2AB,ISSP2BB,ISSP4,ISSP5,IQCLUTAB,JSDA4,
     &                    JQCLUTAB,NIJSTAB,NIJSTABMAX,NSSP4,NSSP4MAX,
     &                    NSSP1,NSSP2B,NSSPABS,NSSPDIR,NRGNT,IL3RGNT,
     &                    LM3RGNT,RGNT_CLU,LMAX,ERYD,P,C,ALAT,NLM,
     &                    IND0QCLU,IND0QCLU0,SAMENLQ,CLUIPH,CLURYLM_NNP,
     &                    JLTAB,NLTAB,HLTAB,FL1L2,NLMAX,NLM3MAX,
     &                    NRGNT12MAX,NRGNT123MAX)
C   ********************************************************************
C   *                                                                  *
C   *     set up the real space free electron G-matrix                 *
C   *     in the NON-relativistic (l,m_l)-representation               *
C   *     using  REAL  spherical harmonics                             *
C   *     NOTE:    - G(i,L;i',L')   is stored in   MG0MAT              *
C   *                                                                  *
C   *    SAMENLQ = TRUE   for all sites the same l-cut off is used     *
C   *                                                                  *
C   *    SAMENLQ = FALSE  different sites may have different l-cut off *
C   *                     the set of inequivalent blocks is stored     *
C   *                     with maximum l-cut off in  MG0TMP            *
C   *                     then the reduced blocks are copied           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C*--CLUG0MAT29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUG0MAT')
C
C Dummy arguments
C
      REAL*8 ALAT,C
      COMPLEX*16 ERYD,P
      INTEGER LMAX,NIJSTABMAX,NKKR,NLM,NLM3MAX,NLMAX,NQCLU,NRGNT123MAX,
     &        NRGNT12MAX,NSSP1,NSSP2A,NSSP2B,NSSP4MAX,NSSPABS,NSSPDIR
      LOGICAL SAMENLQ
      COMPLEX*16 CLUIPH(2*NLMAX-1,NSSPABS),HLTAB(2*NLMAX-1),
     &           JLTAB(2*NLMAX-1),MG0MAT(NKKR,NKKR),NLTAB(2*NLMAX-1)
      REAL*8 CLURYLM_NNP(NLM3MAX,NSSPDIR),DSSP(0:NSSP1),
     &       FL1L2((2*NLMAX-1)**2,(2*NLMAX-1)**2),RGNT_CLU(NRGNT123MAX)
      INTEGER IL3RGNT(NRGNT123MAX),IND0QCLU(NQCLU),IND0QCLU0(NQCLU),
     &        IQCLUTAB(NIJSTABMAX,NSSP1),ISDA4(NSSP4MAX,NSSPABS),
     &        ISSP2AB(NSSP2A),ISSP2BB(NSSP2A),ISSP4(NSSP4MAX,NSSPABS),
     &        ISSP5(NSSP4MAX,NSSPABS),JQCLUTAB(NIJSTABMAX,NSSP1),
     &        JSDA4(NSSP4MAX,NSSPABS),LM3RGNT(NRGNT123MAX),
     &        NIJSTAB(NSSP1),NLMQCLU(NQCLU),NRGNT(NRGNT12MAX),
     &        NSSP4(NSSPABS)
C
C Local variables
C
      COMPLEX*16 ARG,GH,MG0TMP(:,:)
      REAL*8 G
      INTEGER I,I01,I02,I0A,I0B,I1,I11,I12,I2,I2B,IA,IA_ERR,ID,IL,IL3,
     &        ILM12,ILM123,ILM3,IQ,IQ2,J,J01,J02,J0A,J0B,J2,JQ,JQ2,LM1,
     &        LM2,LM3,N,NLMI2,NLMIB,NLMJ2,NLMJB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MG0TMP
C
      CALL CINIT(NKKR*NKKR,MG0MAT)
C
C---------------------------------------------------- calculate momentum
C
      CALL GET_MOMENTUM(IREL,C,ERYD,P)
C
C----------------------------------------- and tabulate Hankel functions
C------------------------------------------ omit site-diagonal case IA=1
C
      DO IA = 2,NSSPABS
C
         ARG = P*DSSP(ISSP4(1,IA))*ALAT
C
         CALL BESHAN(HLTAB,JLTAB,NLTAB,ARG,2*LMAX)
C
         DO IL = 1,(2*LMAX+1)
            CLUIPH(IL,IA) = CI*P*HLTAB(IL)
         END DO
C
      END DO
C
      IF ( SAMENLQ ) THEN
C=======================================================================
C
C----------------------------- set up all inequivalent G(S,S'') - blocks
C------------------------------------------ omit site-diagonal case IA=1
C
         ILM12 = 0
         ILM123 = 0
C
         DO LM1 = 1,NLM
            DO LM2 = 1,LM1
C
               ILM12 = ILM12 + 1
               DO ILM3 = 1,NRGNT(ILM12)
C
                  ILM123 = ILM123 + 1
                  IL3 = IL3RGNT(ILM123)
                  LM3 = LM3RGNT(ILM123)
                  G = RGNT_CLU(ILM123)
                  DO IA = 2,NSSPABS
                     GH = G*CLUIPH(IL3,IA)
                     DO ID = 1,NSSP4(IA)
                        I = IND0QCLU(ISDA4(ID,IA)) + LM1
                        J = IND0QCLU(JSDA4(ID,IA)) + LM2
                        N = ISSP5(ID,IA)
                        MG0MAT(I,J) = MG0MAT(I,J)
     &                                - GH*CLURYLM_NNP(LM3,N)
                     END DO
                  END DO
C
               END DO
            END DO
         END DO
C
C----------------------- fill up matrices G(IS,LM1;JS,LM2) for LM2 > LM1
C------------------------ use inversion symmetry to get G(JS,LM1;IS,LM2)
C
         DO I2B = 1,NSSP2B
            I0A = IND0QCLU(IQCLUTAB(1,ISSP2AB(I2B)))
            J0A = IND0QCLU(JQCLUTAB(1,ISSP2AB(I2B)))
            I0B = IND0QCLU(IQCLUTAB(1,ISSP2BB(I2B)))
            J0B = IND0QCLU(JQCLUTAB(1,ISSP2BB(I2B)))
C
            DO J = 1,NLM
               DO I = 1,J
C
                  MG0MAT(I0A+I,J0A+J) = FL1L2(I,J)*MG0MAT(I0A+J,J0A+I)
C
                  MG0MAT(I0B+I,J0B+J) = MG0MAT(I0A+J,J0A+I)
                  MG0MAT(I0B+J,J0B+I) = MG0MAT(I0A+I,J0A+J)
C
               END DO
            END DO
C
         END DO
C
C------------------------------------------------- copy identical blocks
C
         DO I1 = 2,NSSP1
            I01 = IND0QCLU(IQCLUTAB(1,I1))
            J01 = IND0QCLU(JQCLUTAB(1,I1))
C
            DO I2 = 2,NIJSTAB(I1)
               I02 = IND0QCLU(IQCLUTAB(I2,I1))
               J02 = IND0QCLU(JQCLUTAB(I2,I1))
C
               MG0MAT((I02+1):(I02+NLM),(J02+1):(J02+NLM))
     &            = MG0MAT((I01+1):(I01+NLM),(J01+1):(J01+NLM))
C
            END DO
         END DO
C
C=======================================================================
C
      ELSE
C
C=======================================================================
C
         N = NLM*NQCLU
C
         ALLOCATE (MG0TMP(N,N),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MG0TMP')
C
         CALL CINIT(N*N,MG0TMP)
C
C----------------------------- set up all inequivalent G(S,S'') - blocks
C------------------------------------------ omit site-diagonal case IA=1
C
         ILM12 = 0
         ILM123 = 0
         DO LM1 = 1,NLM
            DO LM2 = 1,LM1
C
               ILM12 = ILM12 + 1
               DO ILM3 = 1,NRGNT(ILM12)
C
                  ILM123 = ILM123 + 1
                  IL3 = IL3RGNT(ILM123)
                  LM3 = LM3RGNT(ILM123)
                  G = RGNT_CLU(ILM123)
                  DO IA = 2,NSSPABS
                     GH = G*CLUIPH(IL3,IA)
                     DO ID = 1,NSSP4(IA)
                        I = IND0QCLU0(ISDA4(ID,IA)) + LM1
                        J = IND0QCLU0(JSDA4(ID,IA)) + LM2
                        N = ISSP5(ID,IA)
                        MG0TMP(I,J) = MG0TMP(I,J)
     &                                - GH*CLURYLM_NNP(LM3,N)
                     END DO
                  END DO
C
               END DO
            END DO
         END DO
C
C----------------------- fill up matrices G(IS,LM1;JS,LM2) for LM2 > LM1
C------------------------ use inversion symmetry to get G(JS,LM1;IS,LM2)
C
         DO I2B = 1,NSSP2B
            I0A = IND0QCLU0(IQCLUTAB(1,ISSP2AB(I2B)))
            J0A = IND0QCLU0(JQCLUTAB(1,ISSP2AB(I2B)))
            I0B = IND0QCLU0(IQCLUTAB(1,ISSP2BB(I2B)))
            J0B = IND0QCLU0(JQCLUTAB(1,ISSP2BB(I2B)))
C
            DO J = 1,NLM
               DO I = 1,NLM
C
                  MG0MAT(I0A+I,J0A+J) = FL1L2(I,J)*MG0MAT(I0A+J,J0A+I)
C
                  MG0MAT(I0B+I,J0B+J) = MG0MAT(I0A+J,J0A+I)
                  MG0MAT(I0B+J,J0B+I) = MG0MAT(I0A+I,J0A+J)
C
               END DO
            END DO
         END DO
C
C------------------------------------------------- copy identical blocks
C
         DO I2B = 1,NSSP2B
            IQ = IQCLUTAB(1,ISSP2AB(I2B))
            JQ = JQCLUTAB(1,ISSP2AB(I2B))
C
            I0A = IND0QCLU0(IQ)
            J0A = IND0QCLU0(JQ)
C
            I0B = IND0QCLU(IQ)
            J0B = IND0QCLU(JQ)
            NLMIB = NLMQCLU(IQ)
            NLMJB = NLMQCLU(JQ)
C
            DO J = 1,NLMJB
               DO I = 1,NLMIB
                  MG0MAT(I0B+I,J0B+J) = MG0TMP(I0A+I,J0A+J)
               END DO
            END DO
         END DO
C
C
         DO I2B = 1,NSSP2B
            IQ = IQCLUTAB(1,ISSP2BB(I2B))
            JQ = JQCLUTAB(1,ISSP2BB(I2B))
C
            I0A = IND0QCLU0(IQ)
            J0A = IND0QCLU0(JQ)
C
            I0B = IND0QCLU(IQ)
            J0B = IND0QCLU(JQ)
            NLMIB = NLMQCLU(IQ)
            NLMJB = NLMQCLU(JQ)
C
            DO J = 1,NLMJB
               DO I = 1,NLMIB
                  MG0MAT(I0B+I,J0B+J) = MG0TMP(I0A+I,J0A+J)
               END DO
            END DO
         END DO
C
         DO I1 = 2,NSSP1
            I11 = IND0QCLU0(IQCLUTAB(1,I1)) + 1
            J01 = IND0QCLU0(JQCLUTAB(1,I1))
C
            DO I2 = 2,NIJSTAB(I1)
               IQ2 = IQCLUTAB(I2,I1)
               JQ2 = JQCLUTAB(I2,I1)
               I12 = IND0QCLU(IQ2) + 1
               J2 = IND0QCLU(JQ2)
               NLMI2 = NLMQCLU(IQ2)
               NLMJ2 = NLMQCLU(JQ2)
               DO J = J01 + 1,J01 + NLMJ2
C
                  J2 = J2 + 1
                  CALL ZCOPY(NLMI2,MG0TMP(I11,J),1,MG0MAT(I12,J2),1)
C
               END DO
            END DO
         END DO
C
         DEALLOCATE (MG0TMP)
C=======================================================================
      END IF
C
      END
C*==cluchkg0.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUCHKG0(WA,NKKR,NQCLU,NLMQCLU,NLQCLU,DSSP,ISDA4,
     &                    NSSP2A,ISSP2AB,ISSP2BB,ISSP4,ISSP5,IQCLUTAB,
     &                    JSDA4,JQCLUTAB,NIJSTAB,NIJSTABMAX,NSSP4,
     &                    NSSP4MAX,NSSP1,NSSP2B,NSSPABS,NSSPDIR,NRGNT,
     &                    IL3RGNT,LM3RGNT,RGNT_CLU,LMAX,ERYD,P,C,ALAT,
     &                    NLM,IND0QCLU,IND0QCLU0,SAMENLQ,CLUIPH,
     &                    CLURYLM_NNP,JLTAB,NLTAB,HLTAB,FL1L2,NLMAX,
     &                    NLM3MAX,NRGNT12MAX,NRGNT123MAX)
C   ********************************************************************
C   *                                                                  *
C   *          check the set of the G-matrix                           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CLUCHKG0322
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='clug0mat.f')
      REAL*8 TOL
      PARAMETER (TOL=1D-12)
C
C Dummy arguments
C
      REAL*8 ALAT,C
      COMPLEX*16 ERYD,P
      INTEGER LMAX,NIJSTABMAX,NKKR,NLM,NLM3MAX,NLMAX,NQCLU,NRGNT123MAX,
     &        NRGNT12MAX,NSSP1,NSSP2A,NSSP2B,NSSP4MAX,NSSPABS,NSSPDIR
      LOGICAL SAMENLQ
      COMPLEX*16 CLUIPH(2*NLMAX-1,NSSPABS),HLTAB(2*NLMAX-1),
     &           JLTAB(2*NLMAX-1),NLTAB(2*NLMAX-1),WA(NKKR,NKKR)
      REAL*8 CLURYLM_NNP((2*NLMAX-1)**2,NSSPDIR),DSSP(0:NSSP1),
     &       FL1L2((2*NLMAX-1)**2,(2*NLMAX-1)**2),RGNT_CLU(NRGNT123MAX)
      INTEGER IL3RGNT(NRGNT123MAX),IND0QCLU(NQCLU),IND0QCLU0(NQCLU),
     &        IQCLUTAB(NIJSTABMAX,NSSP1),ISDA4(NSSP4MAX,NSSPABS),
     &        ISSP2AB(NSSP2A),ISSP2BB(NSSP2A),ISSP4(NSSP4MAX,NSSPABS),
     &        ISSP5(NSSP4MAX,NSSPABS),JQCLUTAB(NIJSTABMAX,NSSP1),
     &        JSDA4(NSSP4MAX,NSSPABS),LM3RGNT(NRGNT123MAX),
     &        NIJSTAB(NSSP1),NLMQCLU(NQCLU),NLQCLU(NQCLU),
     &        NRGNT(NRGNT12MAX),NSSP4(NSSPABS)
C
C Local variables
C
      INTEGER I,IA,IA_ERR,IB,II,ILM,IQCLU,J,JA,JB,JJ,JLM,JQCLU,LMB,N,
     &        NDEL,NDIM,NSIG
      REAL*8 TIME,TIME0
      COMPLEX*16 WB(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WB
C
      N = NLM*NQCLU
      ALLOCATE (WB(N,N),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WB')
C
      WRITE (6,99002) 'matrix dimension       N',N
C
C=======================================================================
      IF ( SAMENLQ ) THEN
C
         CALL CPU_TIME(TIME0)
C
         CALL CLUG0MAT(WA,NKKR,NQCLU,NLMQCLU,DSSP,ISDA4,NSSP2A,ISSP2AB,
     &                 ISSP2BB,ISSP4,ISSP5,IQCLUTAB,JSDA4,JQCLUTAB,
     &                 NIJSTAB,NIJSTABMAX,NSSP4,NSSP4MAX,NSSP1,NSSP2B,
     &                 NSSPABS,NSSPDIR,NRGNT,IL3RGNT,LM3RGNT,RGNT_CLU,
     &                 LMAX,ERYD,P,C,ALAT,NLM,IND0QCLU,IND0QCLU0,
     &                 SAMENLQ,CLUIPH,CLURYLM_NNP,JLTAB,NLTAB,HLTAB,
     &                 FL1L2,NLMAX,NLM3MAX,NRGNT12MAX,NRGNT123MAX)
C
         CALL CPU_TIME(TIME)
C
         WRITE (6,99001) TIME - TIME0,SAMENLQ
C
         SAMENLQ = .FALSE.
         CALL CPU_TIME(TIME0)
C
         CALL CLUG0MAT(WB,NKKR,NQCLU,NLMQCLU,DSSP,ISDA4,NSSP2A,ISSP2AB,
     &                 ISSP2BB,ISSP4,ISSP5,IQCLUTAB,JSDA4,JQCLUTAB,
     &                 NIJSTAB,NIJSTABMAX,NSSP4,NSSP4MAX,NSSP1,NSSP2B,
     &                 NSSPABS,NSSPDIR,NRGNT,IL3RGNT,LM3RGNT,RGNT_CLU,
     &                 LMAX,ERYD,P,C,ALAT,NLM,IND0QCLU,IND0QCLU0,
     &                 SAMENLQ,CLUIPH,CLURYLM_NNP,JLTAB,NLTAB,HLTAB,
     &                 FL1L2,NLMAX,NLM3MAX,NRGNT12MAX,NRGNT123MAX)
C
         CALL CPU_TIME(TIME)
C
         WRITE (6,99001) TIME - TIME0,SAMENLQ
C
         WRITE (6,*)
         WRITE (6,99002) 'comparing matrices '
         WRITE (6,99002) 'A calculated with old mode for SAMENLQ = T'
         WRITE (6,99002) 'B calculated with new mode for SAMENLQ = F'
         WRITE (6,*)
C
         NDEL = 0
         NSIG = 0
         DO J = 1,NKKR
            DO I = 1,NKKR
C
               IF ( ABS(WA(I,J)-WB(I,J)).GT.TOL ) NDEL = NDEL + 1
C
               IF ( ABS(WA(I,J)).GT.TOL ) THEN
                  IF ( ABS(WA(I,J)+WB(I,J)).LT.TOL ) NSIG = NSIG + 1
               END IF
C
            END DO
         END DO
C
         WRITE (6,99002) 'different elements  DEL ',NDEL
         WRITE (6,99002) 'different signs     SIG ',NSIG
         WRITE (6,99002) 'number of all elements  ',NKKR**2
C
C
      END IF
C=======================================================================
C
      WRITE (6,*)
      WRITE (6,99002) 'reducing dimension for every 3rd block'
      WRITE (6,*)
      NDIM = 0
      DO IQCLU = 1,NQCLU
         IF ( MOD(IQCLU,3).EQ.0 ) NLQCLU(IQCLU) = 2
C
         NLMQCLU(IQCLU) = NLQCLU(IQCLU)**2
C
         IF ( IQCLU.EQ.1 ) THEN
            IND0QCLU(1) = 0
            IND0QCLU0(1) = 0
         ELSE
            IND0QCLU(IQCLU) = IND0QCLU(IQCLU-1) + NLMQCLU(IQCLU-1)
            IND0QCLU0(IQCLU) = IND0QCLU0(IQCLU-1) + NLM
         END IF
         NDIM = NDIM + NLMQCLU(IQCLU)
      END DO
      WRITE (6,99002) 'old dimension   NKKR  ',NKKR
      WRITE (6,99002) 'new dimension     NDIM  ',NDIM
C
      SAMENLQ = .FALSE.
C
      CALL CPU_TIME(TIME0)
C
      CALL CLUG0MAT(WB,NKKR,NQCLU,NLMQCLU,DSSP,ISDA4,NSSP2A,ISSP2AB,
     &              ISSP2BB,ISSP4,ISSP5,IQCLUTAB,JSDA4,JQCLUTAB,NIJSTAB,
     &              NIJSTABMAX,NSSP4,NSSP4MAX,NSSP1,NSSP2B,NSSPABS,
     &              NSSPDIR,NRGNT,IL3RGNT,LM3RGNT,RGNT_CLU,LMAX,ERYD,P,
     &              C,ALAT,NLM,IND0QCLU,IND0QCLU0,SAMENLQ,CLUIPH,
     &              CLURYLM_NNP,JLTAB,NLTAB,HLTAB,FL1L2,NLMAX,NLM3MAX,
     &              NRGNT12MAX,NRGNT123MAX)
C
      CALL CPU_TIME(TIME)
C
      WRITE (6,99001) TIME - TIME0,SAMENLQ
C
      WRITE (6,*)
      WRITE (6,99002) 'writing matrices on file '
      WRITE (6,99002) '80:  old full matrix A'
      WRITE (6,99002) '81:  new reduced matrix B'
      WRITE (6,*)
C
      DO I = 1,NKKR
         DO J = 1,NKKR
            WRITE (80,'(2i4,2e22.10)') I,J,WA(I,J)
         END DO
      END DO
C
      I = 0
      DO IQCLU = 1,NQCLU
         II = IND0QCLU0(IQCLU)
         DO ILM = 1,NLMQCLU(IQCLU)
            I = I + 1
            II = II + 1
C
            J = 0
            DO JQCLU = 1,NQCLU
               JJ = IND0QCLU0(JQCLU)
               DO JLM = 1,NLMQCLU(JQCLU)
                  J = J + 1
                  JJ = JJ + 1
C
                  WRITE (81,'(2i4,2e22.10)') II,JJ,WB(I,J)
C
               END DO
            END DO
         END DO
      END DO
C
      WRITE (6,*)
      WRITE (6,99002) 'reducing old full matrix A'
      WRITE (6,99002) 'by deleting columns and rows'
      WRITE (6,*)
C
      JB = 0
      DO IQCLU = 1,NQCLU
         JA = IND0QCLU0(IQCLU)
         DO LMB = 1,NLMQCLU(IQCLU)
            JA = JA + 1
            JB = JB + 1
            DO I = 1,NKKR
               WA(I,JB) = WA(I,JA)
            END DO
         END DO
      END DO
C
      IB = 0
      DO IQCLU = 1,NQCLU
         IA = IND0QCLU0(IQCLU)
         DO LMB = 1,NLMQCLU(IQCLU)
            IA = IA + 1
            IB = IB + 1
            DO J = 1,NDIM
               WA(IB,J) = WA(IA,J)
            END DO
         END DO
      END DO
C
      WRITE (6,*)
      WRITE (6,99002) 'comparing reduced matrix A with new matrix B'
      WRITE (6,*)
C
      NDEL = 0
      NSIG = 0
      DO J = 1,NDIM
         DO I = 1,NDIM
C
            IF ( ABS(WA(I,J)-WB(I,J)).GT.TOL ) NDEL = NDEL + 1
C
            IF ( ABS(WA(I,J)).GT.TOL ) THEN
               IF ( ABS(WA(I,J)+WB(I,J)).LT.TOL ) NSIG = NSIG + 1
            END IF
C
         END DO
      END DO
C
      WRITE (6,99002) 'different elements  DEL ',NDEL
      WRITE (6,99002) 'different signs     SIG ',NSIG
      WRITE (6,99002) 'number of all elements  ',NDIM**2
C
      DEALLOCATE (WB)
C
      STOP
C
99001 FORMAT (/,10X,'<CLUCHKG0>:  time to set up G-set',F6.3,' secs',
     &        '  for SAMENLQ = ',L1,/)
99002 FORMAT (10X,'<CLUCHKG0>:    ',A,I10)
      END
