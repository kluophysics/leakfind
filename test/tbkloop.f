C*==tbkloop.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBKLOOP(TAUQX,TSSQX,DSSQX,DSSQX_L,DSSQX_R,GREF_I1,
     &                   ICLU_REF_IQTB,WA,WB,NXM)
C   ********************************************************************
C   *                                                                  *
C   *  perform the loop over k-points                                  *
C   *                                                                  *
C   *  INPUT:                                                          *
C   *          GREF_I1  G-matrix of the real space TB refernce cluster *
C   *                                                                  *
C   *  OUTPUT:                                                         *
C   *          TAUQX  the TAU-matrix obtained by BZ-integration        *
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
C   *        the result is converted to the site-diagonal   TAUQ       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB,IM_KTAB,CMPLX_KVEC
      USE MOD_CALCMODE,ONLY:IREL,SYMCHECK
      USE MOD_ANGMOM,ONLY:NLMMAX,NKM,NLM
      USE MOD_SITES,ONLY:NQ_L,NQ_R,NQTB
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NXM
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_L(NXM,NXM,NQ_L),
     &           DSSQX_R(NXM,NXM,NQ_R),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           TAUQX(NXM,NXM,NQTB),TSSQX(NXM,NXM,NQTB),WA(NXM,NXM),
     &           WB(NXM,NXM)
      INTEGER ICLU_REF_IQTB(NQTB)
C
C Local variables
C
      COMPLEX*16 GLLKE(:,:),GREFKE(:,:)
      INTEGER I,IA_ERR,IK,IK0,IKVEC,IL,IL0,IQTB,IX,IX0,J,JK,JK0,JL,JL0,
     &        JQTB,JX
      REAL*8 IM_KVEC(3),KVEC(3),WK,WKSUM
C
C*** End of declarations rewritten by SPAG
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE GLLKE,GREFKE
C
      IF ( IPRINT.GE.1 ) WRITE (6,*) 
     &                     '>>> TBKLOOP: integrate over BZ and get GFUN'
C----------------------------------------------------------------------
C
C
      ALLOCATE (GLLKE(NKKR_TB,NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '      <TBKLOOP> allocate GLLKE'
      IF ( IREL.GT.2 ) THEN
         ALLOCATE (GREFKE(NKKRNR_TB,NKKRNR_TB),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP '      <TBKLOOP> allocate GREFKE'
      ELSE IF ( NKKR_TB.NE.NKKRNR_TB ) THEN
         STOP '<TBKLOOP>:  NKKR_TB <> NKKRNR_TB'
      END IF
C
      WKSUM = SUM(WKTAB(1:NKTAB))
C
      CALL CINIT(NXM*NXM*NQTB,TAUQX)
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C                                                          K-points loop
      DO IKVEC = 1,NKTAB
C
         WK = WKTAB(IKVEC)/WKSUM
C
         KVEC(1:3) = KTAB(1:3,IKVEC)
         IF ( CMPLX_KVEC ) IM_KVEC(1:3) = IM_KTAB(1:3,IKVEC)
C
C-----------------------------------------------------------------------
C    Fourier transformation, set KKR matrix M = [-(t)^-1 + G^r]
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
C
            CALL CINIT(NKKR_TB*NKKR_TB,GLLKE)
C
            JK0 = -NKM
            JL0 = -NLM
            DO JQTB = 1,NQTB
               JK0 = JK0 + NKM
               JL0 = JL0 + NLM
               DO J = 1,NLM
                  JK = JK0 + J
                  JL = JL0 + J
C
                  IK0 = -NKM
                  IL0 = -NLM
                  DO IQTB = 1,NQTB
                     IK0 = IK0 + NKM
                     IL0 = IL0 + NLM
                     DO I = 1,NLM
                        IK = IK0 + I
                        IL = IL0 + I
                        GLLKE(IK,JK) = GREFKE(IL,JL)
                        GLLKE(NLM+IK,NLM+JK) = GREFKE(IL,JL)
                     END DO
                  END DO
               END DO
            END DO
C
         END IF
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C             write  G_ref[q,q](k)   for checking the symmetry
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( SYMCHECK ) THEN
            JK0 = -NKM
            JL0 = -NLM
            DO IQTB = 1,NQTB
               JK0 = JK0 + NKM
               JL0 = JL0 + NLM
               DO J = 1,NLM
                  JK = JK0 + J
                  JL = JL0 + J
C
                  IK0 = -NKM
                  IL0 = -NLM
C
                  IK0 = IK0 + NKM
                  IL0 = IL0 + NLM
                  DO I = 1,NLM
                     IK = IK0 + I
                     IL = IL0 + I
                     WRITE (IOTMP) IQTB,I,J,GLLKE(IK,JK)
                  END DO
C
               END DO
            END DO
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C-----------------------------------------------------------------------
C              call decimation routine if requested
C-----------------------------------------------------------------------
C
         IF ( IDECI.EQ.1 ) CALL TBDECIMATE(GLLKE,DSSQX_L,DSSQX_R,
     &        VACFLAG,FACTL,NQ_L,NQ_R,NXM)
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
C
         CALL TBINVERSION(GLLKE,INVMOD,ICHECK,NXM,TAUQX,WK)
C
C
C-----------------------------------------------------------------------
C           get the site-diagonal blocks and convert to TAU(k)
C-----------------------------------------------------------------------
C
         IF ( INVMOD.NE.4 ) THEN
            IX0 = -NXM
            DO IQTB = 1,NQTB
               IX0 = IX0 + NXM
               JX = IX0
               DO J = 1,NXM
                  JX = JX + 1
                  IX = IX0
                  DO I = 1,NXM
                     IX = IX + 1
                     TAUQX(I,J,IQTB) = TAUQX(I,J,IQTB) - WK*GLLKE(IX,JX)
                  END DO
               END DO
            END DO
         END IF
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C             write  TAU[q,q](k)   for checking the symmetry
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( SYMCHECK ) THEN
            IX0 = -NXM
            DO IQTB = 1,NQTB
               IX0 = IX0 + NXM
               JX = IX0
               DO J = 1,NXM
                  JX = JX + 1
                  IX = IX0
                  DO I = 1,NXM
                     IX = IX + 1
                     WRITE (IOTMP) IQTB,I,J, - GLLKE(IX,JX)
                  END DO
               END DO
            END DO
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      END DO
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C-----------------------------------------------------------------------
C    convert the site-diagonal blocks from TAU_Delta via G to TAU
C-----------------------------------------------------------------------
C
      DO IQTB = 1,NQTB
C
C------------- G = 1/(delta t) * TAU_Delta *1/(delta t) - 1/(delta t)
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TAUQX(1,1,IQTB),NXM,
     &              DSSQX(1,1,IQTB),NXM,C0,WB,NXM)
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WB,NXM,
     &              C0,WA,NXM)
C
         WB(1:NXM,1:NXM) = WA(1:NXM,1:NXM) - DSSQX(1:NXM,1:NXM,IQTB)
C
C--------------------------------------------------- TAU = t * G * t + t
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WB,NXM,TSSQX(1,1,IQTB),NXM,
     &              C0,WA,NXM)
C
         CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TSSQX(1,1,IQTB),NXM,WA,NXM,
     &              C0,WB,NXM)
C
         DO J = 1,NXM
            DO I = 1,NXM
               TAUQX(I,J,IQTB) = WB(I,J) + TSSQX(I,J,IQTB)
            END DO
         END DO
C
      END DO
C
      END
C*==tbgrefke.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBGREFKE(KVEC,IM_KVEC,NQTB,GREF_I1,GREFKE,
     &                    ICLU_REF_IQTB,NKKRNR_RSMAX)
C   ********************************************************************
C   *                                                                  *
C   *  Performs Fourier transformation for                             *
C   *  NON-relativistic reference Green function                       *
C   *                                                                  *
C   *  Modifications according to H. Hoehler                           *
C   *     define Fourier transformation as                             *
C   *                                                                  *
C   *              /        n   1                                      *
C   *   G mu mu'= | sum_n G mu mu' exp(-ikR ) +                        *
C   *     L  L'   \         L   L'         n                           *
C   *                                                                  *
C   *                       n   1                \    1                *
C   *               sum_n G mu mu' exp(-ik(-R )) | * ---               *
C   *                       L   L'           n  /     2                *
C   *                                                                  *
C   *   this operation has to be done to satisfy the point symmetry;   *
C   *   the application of the fourier transformation is just an       *
C   *   approximation for the tb system, since the transl. invariance  *
C   *   is not satisfied --> force it by R, -R                         *
C   *                                                                  *
C   *      Here we do   --                  nn'                        *
C   *                   \                   ii'          ii'           *
C   *                   /  exp(+ik(R  -R ))G   (E)  =   G   (k,E)      *
C   *                   --          n'  n   LL'          LL'           *
C   *                   n'                                             *
C   *                                                                  *
C   *     a minus sign must be included here because                   *
C   *     RI0CLU_GEO is not symmetric around each atom.                *
C   *     The minus comes from the fact that the repulsive             *
C   *     potential GF is calculated for 0n and not n0                 *
C   *     and that is why we need a minus sign extra!                  *
C   *                                                                  *
C   *                                                                  *
C   *   SB, HE 15/03/2013                                              *
C   *                                                                  *
C   *   The prefactor of ->k changed to   exp(+ik(Rn'-Rn))Gnn'         *
C   *   standard KKR. This is important only for site-off diagonal     *
C   *   properties.                                                    *
C   *                                                                  *
C   *   The consistency of this change with the TB-formulation         *
C   *   should still be checked.                                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:CMPLX_KVEC
      USE MOD_CONSTANTS,ONLY:C0,CI,CI2PI
      USE MOD_ANGMOM,ONLY:NLMMAX,NLM
      USE MOD_TB,ONLY:NKKRNR_TB
      USE MOD_TBCLU,ONLY:RI0CLU_GEO,ICLU_GEO_ICLU_REF,NQCLU_ICLU_REF,
     &    IQ_IQCLU_GEO_IQCNTR,NCLU_REF
      USE MOD_SITES,ONLY:IQBOT_TB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKKRNR_RSMAX,NQTB
      COMPLEX*16 GREFKE(NKKRNR_TB,NKKRNR_TB),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF)
      INTEGER ICLU_REF_IQTB(NQTB)
      REAL*8 IM_KVEC(3),KVEC(3)
C
C Local variables
C
      COMPLEX*16 CKI,EIKR,EIKRM,GIJ,TT
      INTEGER I,ICLU_GEO,ICLU_REF,II,II0,IICLU0,IQCLU_REF,IQTB,J,JJ,JJ0,
     &        JQ,JQTB
C
C*** End of declarations rewritten by SPAG
C
      CALL CINIT(NKKRNR_TB*NKKRNR_TB,GREFKE)
C
C  AENDERUNG              SB, HE 15/03/2013   "- CI2PI"  ==>>  "+ CI2PI"
C  RUECKGAENGIG GEMACHT UND SVEN'S SETZUNG k=-k VERWENDET
C
      KVEC(:) = -KVEC(:)
C
C ======================================================================
      DO JQTB = 1,NQTB
C
         ICLU_REF = ICLU_REF_IQTB(JQTB)
         ICLU_GEO = ICLU_GEO_ICLU_REF(ICLU_REF)
         JQ = JQTB + IQBOT_TB - 1
C
         JJ0 = (JQTB-1)*NLM
C
C ----------------------------------------------------------------------
C          loop over cluster sites around TB-site JQTB
C ----------------------------------------------------------------------
C
         DO IQCLU_REF = 1,NQCLU_ICLU_REF(ICLU_REF)
C
            IQTB = IQ_IQCLU_GEO_IQCNTR(IQCLU_REF,JQ) - IQBOT_TB + 1
C
            IF ( IQTB.GE.1 .AND. IQTB.LE.NQTB ) THEN
C
C                                              --> for COMPLEX k-vectors
C                         SB, HE 15/03/2013   "- CI2PI"  ==>>  "+ CI2PI"
               IF ( CMPLX_KVEC ) THEN
C
                  TT = C0
                  DO I = 1,3
                     CKI = KVEC(I) + CI*IM_KVEC(I)
                     TT = TT - CI2PI*CKI*RI0CLU_GEO(I,IQCLU_REF,
     &                    ICLU_GEO)
                  END DO
                  EIKR = EXP(TT)*0.5D0
C     EIKR = (EIKR+DCONJG(EIKR))*0.5D0
                  EIKRM = DCONJG(EIKR)
               ELSE
C                                                 --> for REAL k-vectors
C                         SB, HE 15/03/2013   "- CI2PI"  ==>>  "+ CI2PI"
                  TT = C0
                  DO I = 1,3
                     TT = TT - CI2PI*KVEC(I)
     &                    *RI0CLU_GEO(I,IQCLU_REF,ICLU_GEO)
                  END DO
C
C                  EIKR = (EXP(TT)+EXP(-TT))*0.5D0
                  EIKR = (EXP(TT))*0.5D0
                  EIKRM = (EXP(-TT))*0.5D0
C
C
               END IF
C
               IICLU0 = (IQCLU_REF-1)*NLM
C
               II0 = (IQTB-1)*NLM
C
               DO J = 1,NLM
                  JJ = JJ0 + J
                  DO I = 1,NLM
                     GIJ = GREF_I1(IICLU0+I,J,ICLU_REF)
                     II = II0 + I
                     GREFKE(II,JJ) = GREFKE(II,JJ) + EIKR*GIJ
                     GREFKE(JJ,II) = GREFKE(JJ,II) + EIKRM*GIJ
                  END DO
               END DO
            END IF
C
         END DO
C ----------------------------------------------------------------------
C
      END DO
C     ================================================================
      END
C*==tbdecimate.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBDECIMATE(GLLKE,DSSQX_L,DSSQX_R,VACFLAG,FACTL,NQ_L,
     &                      NQ_R,NXM)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TB,ONLY:NKKR_TB,NKMSLAY,NPLAY,NSLAY_PER_PLAY,TBOPT_ONEBULK
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NQ_L,NQ_R,NXM
      COMPLEX*16 DSSQX_L(NXM,NXM,NQ_L),DSSQX_R(NXM,NXM,NQ_R),
     &           FACTL(NXM,NXM),GLLKE(NKKR_TB,NKKR_TB)
      LOGICAL VACFLAG(2)
C
C Local variables
C
      COMPLEX*16 A1(:,:),AN(:,:),B1(:,:),BN(:,:),C1(:,:),CN(:,:),X1(:,:)
     &           ,XN(:,:)
      REAL*8 ERRMAX
      INTEGER ICHCK,IHOST,II1,II2,IL1,IL2,IP1,IP1T,IP2,IP2T,ITERMAX,
     &        LDI1,LDI1T,LDI2,LDI2T,LM1,LM2
      LOGICAL INITIALIZE
      SAVE A1,AN,B1,BN,C1,CN,ERRMAX,ICHCK,ITERMAX,X1,XN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE A1,B1,C1,X1,AN,BN,CN,XN
C
      DATA INITIALIZE/.TRUE./
C
C ----------------------------------------------------------------------
C     Parameters for the "decimation" technique.
      IF ( INITIALIZE ) THEN
         ITERMAX = 300
         ERRMAX = 1.0D-18
         ICHCK = 1
         ALLOCATE (A1(NKMSLAY,NKMSLAY),AN(NKMSLAY,NKMSLAY))
         ALLOCATE (B1(NKMSLAY,NKMSLAY),BN(NKMSLAY,NKMSLAY))
         ALLOCATE (C1(NKMSLAY,NKMSLAY),CN(NKMSLAY,NKMSLAY))
         ALLOCATE (X1(NKMSLAY,NKMSLAY),XN(NKMSLAY,NKMSLAY))
C
         INITIALIZE = .FALSE.
      END IF
C ----------------------------------------------------------------------
      IF ( .NOT.VACFLAG(1) ) THEN
C
C Get the matrix B1
C
         CALL TBBOFM(1,1,B1,NKMSLAY,GLLKE,NKKR_TB)
C
C Now Subtract t-mat of left host
         DO IP1 = 1,NSLAY_PER_PLAY
            IHOST = MOD(IP1-1,NQ_L) + 1
            DO LM1 = 1,NXM
               DO LM2 = 1,NXM
                  IL1 = NXM*(IP1-1) + LM1
                  IL2 = NXM*(IP1-1) + LM2
                  B1(IL1,IL2) = (B1(IL1,IL2)-DSSQX_L(LM1,LM2,IHOST))
               END DO
            END DO
         END DO
C
         CALL TBBOFM(1,2,C1,NKMSLAY,GLLKE,NKKR_TB)
         CALL TBBOFM(2,1,A1,NKMSLAY,GLLKE,NKKR_TB)
C
C     it performs the 'space decimation' iterative procedure.
         CALL TBSURFGF(A1,B1,C1,X1,ITERMAX,ERRMAX,ICHCK)
C     adds to the matrix GLLKE the elements that couples the
C     interface to the two half-spaces.
         DO IP1 = 1,NSLAY_PER_PLAY
            DO IP2 = 1,NSLAY_PER_PLAY
               II1 = IP1
               II2 = IP2
               DO LM1 = 1,NXM
                  DO LM2 = 1,NXM
                     LDI1 = NXM*(IP1-1) + LM1
                     IL1 = NXM*(II1-1) + LM1
                     LDI2 = NXM*(IP2-1) + LM2
                     IL2 = NXM*(II2-1) + LM2
                     GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - X1(LDI1,LDI2)
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      IF ( .NOT.VACFLAG(2) ) THEN
C
C  If 'ONEBULK' is activated then it calculates the xn decimated element
C  from the x1 element: this is just in the case of equal bulks on the
C
         IF ( .NOT.TBOPT_ONEBULK ) THEN
C
C     Get the matrix BN
C
            CALL TBBOFM(NPLAY,NPLAY,BN,NKMSLAY,GLLKE,NKKR_TB)
C
C Now Substract t-mat right host
C Notes : the indexing is easier like that
C
            DO IP1 = 1,NSLAY_PER_PLAY
               IHOST = MOD(IP1-1,NQ_R) + 1
               DO LM1 = 1,NXM
                  DO LM2 = 1,NXM
                     IL1 = NXM*(IP1-1) + LM1
                     IL2 = NXM*(IP1-1) + LM2
                     BN(IL1,IL2) = (BN(IL1,IL2)-DSSQX_R(LM1,LM2,IHOST))
                  END DO
               END DO
            END DO
C
            CALL TBBOFM(NPLAY,NPLAY-1,AN,NKMSLAY,GLLKE,NKKR_TB)
            CALL TBBOFM(NPLAY-1,NPLAY,CN,NKMSLAY,GLLKE,NKKR_TB)
C
C     it performs the 'space decimation' iterative procedure.
C
            CALL TBSURFGF(CN,BN,AN,XN,ITERMAX,ERRMAX,ICHCK)
C
         ELSE
C
            DO IP1 = 1,NSLAY_PER_PLAY
               DO IP2 = 1,NSLAY_PER_PLAY
                  IP1T = (NSLAY_PER_PLAY+1) - IP2
                  IP2T = (NSLAY_PER_PLAY+1) - IP1
                  DO LM1 = 1,NXM
                     DO LM2 = 1,NXM
                        LDI1 = NXM*(IP1-1) + LM1
                        LDI2 = NXM*(IP2-1) + LM2
                        LDI1T = NXM*(IP1T-1) + LM2
                        LDI2T = NXM*(IP2T-1) + LM1
                        XN(LDI1T,LDI2T) = FACTL(LM1,LM2)*X1(LDI1,LDI2)
                     END DO
                  END DO
               END DO
            END DO
C
         END IF
C=======================================================================
C     adds to the matrix GLLKE the elements that couples the
C     interface to the two half-spaces.
C=======================================================================
C
         DO IP1 = 1,NSLAY_PER_PLAY
            DO IP2 = 1,NSLAY_PER_PLAY
               II1 = (NPLAY-1)*NSLAY_PER_PLAY + IP1
               II2 = (NPLAY-1)*NSLAY_PER_PLAY + IP2
               DO LM1 = 1,NXM
                  DO LM2 = 1,NXM
                     LDI1 = NXM*(IP1-1) + LM1
                     IL1 = NXM*(II1-1) + LM1
                     LDI2 = NXM*(IP2-1) + LM2
                     IL2 = NXM*(II2-1) + LM2
                     GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - XN(LDI1,LDI2)
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      END
C*==tbbofm.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBBOFM(PL1,PL2,BLOCK,NSIZE,GIN,NKKR_TB)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKKR_TB,NSIZE,PL1,PL2
      COMPLEX*16 BLOCK(NSIZE,NSIZE),GIN(NKKR_TB,NKKR_TB)
C
C Local variables
C
      INTEGER I1,I1S,I2,I2S
C
C*** End of declarations rewritten by SPAG
C
      I1S = (PL1-1)*NSIZE
      I2S = (PL2-1)*NSIZE
      DO I2 = 1,NSIZE
         DO I1 = 1,NSIZE
            BLOCK(I1,I2) = GIN(I1S+I1,I2S+I2)
         END DO
      END DO
C
      END
C*==tbsurfgf.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBSURFGF(ML,M0,MR,X,ITERMAX,ERRMAX,ICHCK)
C   ********************************************************************
C   *                                                                  *
C   *  solve surface Green's function: F(X)=ML*(M0-X)**(-1)*MR         *
C   *  method: decimation technique                                    *
C   *                                                                  *
C   *  input:  ML,M0,MR - complex rectangular matrices                 *
C   *  output: X        - result, matrix of same type as before        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_TB,ONLY:NKMSLAY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ERRMAX
      INTEGER ICHCK,ITERMAX
      COMPLEX*16 M0(NKMSLAY,NKMSLAY),ML(NKMSLAY,NKMSLAY),
     &           MR(NKMSLAY,NKMSLAY),X(NKMSLAY,NKMSLAY)
C
C Local variables
C
      COMPLEX*16 AA(:,:),ALFA(:,:),BB(:,:),BETA(:,:),CC(:,:),CUNIT(:,:),
     &           EPS(:,:),TEMPIN(:,:),TEMPOUT(:,:),Y1(:,:),Y2(:,:)
      REAL*8 ERR,RSUM,XIM,XRE
      INTEGER I,INFO,IPROBLEM,IPVT(:),ITER,J,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AA,BB,CC,EPS,ALFA,BETA,CUNIT,TEMPIN,TEMPOUT,Y1,Y2,IPVT
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.ALLOCATED(AA) ) THEN
         ALLOCATE (AA(NKMSLAY,NKMSLAY),ALFA(NKMSLAY,NKMSLAY))
         ALLOCATE (BB(NKMSLAY,NKMSLAY),BETA(NKMSLAY,NKMSLAY))
         ALLOCATE (CC(NKMSLAY,NKMSLAY),CUNIT(NKMSLAY,NKMSLAY))
         ALLOCATE (EPS(NKMSLAY,NKMSLAY),TEMPIN(NKMSLAY,NKMSLAY))
         ALLOCATE (TEMPOUT(NKMSLAY,NKMSLAY),Y1(NKMSLAY,NKMSLAY))
         ALLOCATE (Y2(NKMSLAY,NKMSLAY),IPVT(NKMSLAY))
      END IF
C
C
      CALL CINIT(NKMSLAY*NKMSLAY,CUNIT)
      DO N = 1,NKMSLAY
         CUNIT(N,N) = C1
      END DO
C
      CALL ZCOPY(NKMSLAY*NKMSLAY,M0,1,EPS,1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,ML,1,ALFA,1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,MR,1,BETA,1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,M0,1,X,1)
C
      ITER = 1
C
 100  CONTINUE
      CALL ZCOPY(NKMSLAY*NKMSLAY,EPS,1,Y1,1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,Y1,1,TEMPIN,1)
      CALL ZGETRF(NKMSLAY,NKMSLAY,TEMPIN,NKMSLAY,IPVT,INFO)
C
C     aa = eps^-1 * alfa
      CALL ZCOPY(NKMSLAY*NKMSLAY,ALFA,1,TEMPOUT,1)
      CALL ZGETRS('N',NKMSLAY,NKMSLAY,TEMPIN,NKMSLAY,IPVT,TEMPOUT,
     &            NKMSLAY,INFO)
      CALL ZCOPY(NKMSLAY*NKMSLAY,TEMPOUT,1,AA,1)
C
C     bb = eps^-1 * beta
C
      CALL ZCOPY(NKMSLAY*NKMSLAY,BETA,1,TEMPOUT,1)
      CALL ZGETRS('N',NKMSLAY,NKMSLAY,TEMPIN,NKMSLAY,IPVT,TEMPOUT,
     &            NKMSLAY,INFO)
      CALL ZCOPY(NKMSLAY*NKMSLAY,TEMPOUT,1,BB,1)
C
C     alfa_new = alfa * aa
C
      CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,ALFA,NKMSLAY,AA,
     &           NKMSLAY,C0,Y1,NKMSLAY)
C
C     beta_new = beta * bb
C
      CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,BETA,NKMSLAY,BB,
     &           NKMSLAY,C0,Y2,NKMSLAY)
C
C     cc = - alfa * bb
C
      CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,ALFA,NKMSLAY,BB,
     &           NKMSLAY,C0,CC,NKMSLAY)
C
C     x_new = x + cc
C
      CALL ZAXPY(NKMSLAY*NKMSLAY,C1,CC,1,X,1)
C
C     cc = eps + cc
C
      CALL ZAXPY(NKMSLAY*NKMSLAY,C1,CC,1,EPS,1)
C
C     eps_new = cc - beta * aa
C
      CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,BETA,NKMSLAY,AA,
     &           NKMSLAY,C1,EPS,NKMSLAY)
C
      CALL ZCOPY(NKMSLAY*NKMSLAY,Y1,1,ALFA,1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,Y2,1,BETA,1)
C
      RSUM = 0.D0
      IPROBLEM = 0
      DO I = 1,NKMSLAY
         DO J = 1,NKMSLAY
            XRE = DBLE(ALFA(I,J))
            XIM = DIMAG(ALFA(I,J))
            IF ( ABS(XRE).LT.1D20 .OR. ABS(XIM).LT.1D20 ) THEN
               RSUM = RSUM + XRE*XRE + XIM*XIM
            ELSE
               IPROBLEM = IPROBLEM + 1
            END IF
         END DO
      END DO
C
      ERR = DSQRT(RSUM)
C
      IF ( ERR.LT.ERRMAX .OR. ITER.GT.ITERMAX .OR. IPROBLEM.GT.0 ) THEN
C
C
         IF ( IPROBLEM.GT.0 ) THEN
            WRITE (6,*) 
     &                '<TBSURFGF>: WARNING: Decimation seems to diverge'
            CALL ZCOPY(NKMSLAY*NKMSLAY,M0,1,X,1)
         END IF
C
         CALL ZCOPY(NKMSLAY*NKMSLAY,X,1,TEMPIN,1)
         CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT,1,TEMPOUT,1)
         CALL ZGETRF(NKMSLAY,NKMSLAY,TEMPIN,NKMSLAY,IPVT,INFO)
         CALL ZGETRS('N',NKMSLAY,NKMSLAY,TEMPIN,NKMSLAY,IPVT,TEMPOUT,
     &               NKMSLAY,INFO)
         CALL ZCOPY(NKMSLAY*NKMSLAY,TEMPOUT,1,X,1)
C
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,X,NKMSLAY,MR,
     &              NKMSLAY,C0,TEMPIN,NKMSLAY)
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,ML,NKMSLAY,
     &              TEMPIN,NKMSLAY,C0,X,NKMSLAY)
C
         IF ( ITER.GT.ITERMAX ) WRITE (6,FMT=
     &              '('' itermax too small.  iter='',i5,'' ERR='',F9.6)'
     &              ) ITER,ERR
         IF ( ICHCK.EQ.0 ) RETURN
      ELSE
         ITER = ITER + 1
         GOTO 100
      END IF
C
      END
C*==tbinversion.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBINVERSION(GLLKE,INVMOD,ICHECK,NXM,TAUQX,WK)
C   ********************************************************************
C   *                                                                  *
C   *  This subroutine calculates the inversion of a matrix            *
C   *  in 4 different ways depending on the form of the matrix         *
C   *                                                                  *
C   *      INVMOD = 0  ----> total inversion scheme                    *
C   *      INVMOD = 1  ----> band matrix inversion scheme              *
C   *      INVMOD = 2  ----> corner band matrix inversion scheme       *
C   *      INVMOD = 3  ----> sparse matrix inversion scheme            *
C   *      INVMOD = 4  ----> total inversion calculating just          *
C   *                        the diagonal blocks of a matrix           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQTB
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_TB,ONLY:NKKR_TB,NKMSLAY,NPLAY,NSLAY_PER_PLAY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INVMOD,NXM
      REAL*8 WK
      COMPLEX*16 GLLKE(NKKR_TB,NKKR_TB),TAUQX(NXM,NXM,NQTB)
      INTEGER ICHECK(NPLAY,NPLAY)
C
C Local variables
C
      COMPLEX*16 GDI(:,:,:),GDOW(:,:,:),GTEMP(:,:),GUP(:,:,:),UDIA(:)
      INTEGER I,I1,II1,II2,IL1,IL2,IND0QLOC(:),INFO,IP1,IP2,IPVT(:),
     &        IQTB,J,LDI1,LDI2,LM1,LM2,NKMQLOC(:)
      LOGICAL INITIALIZE
      SAVE IND0QLOC,NKMQLOC,UDIA
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
C
      ALLOCATABLE GDI,GDOW,GTEMP,GUP,IPVT,UDIA,IND0QLOC,NKMQLOC
C
      ALLOCATE (GDI(NKMSLAY,NKMSLAY,NPLAY),GDOW(NKMSLAY,NKMSLAY,NPLAY),
     &          GTEMP(NKKR_TB,NKKR_TB),GUP(NKMSLAY,NKMSLAY,NPLAY),
     &          IPVT(NKKR_TB))
C
C------------------------------------------------ total matrix inversion
C
      IF ( INVMOD.EQ.0 ) THEN
C
         DO I = 1,NKKR_TB
            DO J = 1,NKKR_TB
               GTEMP(I,J) = C0
               IF ( I.EQ.J ) GTEMP(I,J) = C1
            END DO
         END DO
C
         CALL ZGETRF(NKKR_TB,NKKR_TB,GLLKE,NKKR_TB,IPVT,INFO)
         CALL ZGETRS('N',NKKR_TB,NKKR_TB,GLLKE,NKKR_TB,IPVT,GTEMP,
     &               NKKR_TB,INFO)
C
         CALL ZCOPY(NKKR_TB*NKKR_TB,GTEMP,1,GLLKE,1)
C
      ELSE IF ( INVMOD.EQ.4 ) THEN
C
         IF ( INITIALIZE ) THEN
            ALLOCATE (UDIA(NKKR_TB),IND0QLOC(NQTB),NKMQLOC(NQTB))
C
            DO IQTB = 1,NQTB
               NKMQLOC(IQTB) = NXM
               IF ( IQTB.EQ.1 ) THEN
                  IND0QLOC(IQTB) = 0
               ELSE
                  IND0QLOC(IQTB) = IND0QLOC(IQTB-1) + NXM
               END IF
            END DO
            INITIALIZE = .FALSE.
         END IF
C
         CALL CINVDIABLK(NKKR_TB,GLLKE,UDIA,TAUQX,NXM,NQTB,IND0QLOC,
     &                   NKMQLOC,.FALSE.,-WK)
C
C
C------------------------------------------  slab or supercell inversion
      ELSE IF ( (INVMOD.GE.1) .AND. (INVMOD.LE.2) ) THEN
C
C
         DO I1 = 1,NPLAY
            DO IP1 = 1,NSLAY_PER_PLAY
               DO IP2 = 1,NSLAY_PER_PLAY
                  II1 = (I1-1)*NSLAY_PER_PLAY + IP1
                  II2 = (I1-1)*NSLAY_PER_PLAY + IP2
                  DO LM1 = 1,NXM
                     DO LM2 = 1,NXM
                        LDI1 = NXM*(IP1-1) + LM1
                        IL1 = NXM*(II1-1) + LM1
                        LDI2 = NXM*(IP2-1) + LM2
                        IL2 = NXM*(II2-1) + LM2
                        GDI(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C     this part now is correct also for supercell geometry
C
C---> upper linear part
         DO I1 = 1,NPLAY
            DO IP1 = 1,NSLAY_PER_PLAY
               DO IP2 = 1,NSLAY_PER_PLAY
                  DO LM1 = 1,NXM
                     DO LM2 = 1,NXM
                        LDI1 = NXM*(IP1-1) + LM1
                        LDI2 = NXM*(IP2-1) + LM2
                        IF ( I1.LE.(NPLAY-1) ) THEN
                           II1 = (I1-1)*NSLAY_PER_PLAY + IP1
                           II2 = I1*NSLAY_PER_PLAY + IP2
                           IL1 = NXM*(II1-1) + LM1
                           IL2 = NXM*(II2-1) + LM2
                           GUP(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
                        ELSE
                           II1 = IP1
                           II2 = (NPLAY-1)*NSLAY_PER_PLAY + IP2
                           IL1 = NXM*(II1-1) + LM1
                           IL2 = NXM*(II2-1) + LM2
                           GDOW(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C
C---> lower linear part
         DO I1 = 1,NPLAY
            DO IP1 = 1,NSLAY_PER_PLAY
               DO IP2 = 1,NSLAY_PER_PLAY
                  DO LM1 = 1,NXM
                     DO LM2 = 1,NXM
                        LDI1 = NXM*(IP1-1) + LM1
                        LDI2 = NXM*(IP2-1) + LM2
                        IF ( I1.LE.(NPLAY-1) ) THEN
                           II1 = I1*NSLAY_PER_PLAY + IP1
                           II2 = (I1-1)*NSLAY_PER_PLAY + IP2
                           IL1 = NXM*(II1-1) + LM1
                           IL2 = NXM*(II2-1) + LM2
                           GDOW(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
                        ELSE
                           II1 = (NPLAY-1)*NSLAY_PER_PLAY + IP1
                           II2 = IP2
                           IL1 = NXM*(II1-1) + LM1
                           IL2 = NXM*(II2-1) + LM2
                           GUP(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C
         IF ( INVMOD.EQ.1 ) THEN
C
            CALL TBINVSLAB(GDI,GUP,GDOW,GLLKE,ICHECK)
C
C          write (6,*) '-------slab calculation--------'
C
         ELSE IF ( INVMOD.EQ.2 ) THEN
                                  ! supercell geometry inversion
C
            CALL TBINVSUPERCELL(GDI,GUP,GDOW,GLLKE,ICHECK)
C
C          write (6,*) '-------supercell calculation--------'
C
         END IF
C
C     NOT YET IMPLEMENTED!!!!!!!!!
C
      END IF
C
      END
C*==tbinvslab.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBINVSLAB(GDI,GUP,GDOW,GIN,ICHECK)
C   ********************************************************************
C   *                                                                  *
C   *    ---> ALGORITM FOR SLAB GEOMETRY                               *
C   *                                                                  *
C   *    ---> factorization D ^-1 = (prod L) * M * (prod U)            *
C   *                                                                  *
C   *         see notes R. Zeller                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_TB,ONLY:NKMSLAY,NKKR_TB,NPLAY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 GDI(NKMSLAY,NKMSLAY,NPLAY),GDOW(NKMSLAY,NKMSLAY,NPLAY),
     &           GIN(NKKR_TB,NKKR_TB),GUP(NKMSLAY,NKMSLAY,NPLAY)
      INTEGER ICHECK(NPLAY,NPLAY)
C
C Local variables
C
      COMPLEX*16 CUNIT(:,:),DINVER(:,:,:),DMAT(:,:,:),E(:,:),F(:,:),
     &           G(:,:),GDIOLD(:,:,:)
      INTEGER I,INFO,IPVT(:),IROW,J,LM,N
      SAVE CUNIT,DINVER,DMAT,E,F,G,GDIOLD,IPVT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CUNIT,DINVER,DMAT,E,F,G,GDIOLD,IPVT
      IF ( .NOT.ALLOCATED(CUNIT) ) THEN
         ALLOCATE (CUNIT(NKMSLAY,NKMSLAY),DINVER(NKMSLAY,NKMSLAY,NPLAY))
         ALLOCATE (DMAT(NKMSLAY,NKMSLAY,NPLAY),E(NKMSLAY,NKMSLAY))
         ALLOCATE (F(NKMSLAY,NKMSLAY),G(NKMSLAY,NKMSLAY))
         ALLOCATE (GDIOLD(NKMSLAY,NKMSLAY,NPLAY))
         ALLOCATE (IPVT(NKMSLAY))
      END IF
C
      CALL CINIT(NKMSLAY*NKMSLAY,E)
      CALL CINIT(NKMSLAY*NKMSLAY,F)
      CALL CINIT(NKMSLAY*NKMSLAY,G)
      CALL CINIT(NKMSLAY*NKMSLAY,CUNIT)
      CALL CINIT(NKMSLAY*NKMSLAY*NPLAY,DMAT)
C
      DO N = 1,NKMSLAY
         CUNIT(N,N) = C1
      END DO
C
      DO N = 1,NPLAY
         CALL ZCOPY(NKMSLAY*NKMSLAY,GDI(1,1,N),1,GDIOLD(1,1,N),1)
      END DO
C
C
C---> calculate D_1 = (M_11)**(-1)
C
C
      CALL ZCOPY(NKMSLAY*NKMSLAY,GDI(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT,1,DMAT(1,1,1),1)
C
      CALL ZGETRF(NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,INFO)
      CALL ZGETRS('N',NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,DMAT(1,1,1),
     &            NKMSLAY,INFO)
C
C---> claculate D_N (2 <= N <= NPLAY)
C
      DO N = 2,NPLAY
C
C---> F = D(N-1) * M1(N-1)
C
C
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,DMAT(1,1,N-1),
     &              NKMSLAY,GUP(1,1,N-1),NKMSLAY,C0,F(1,1),NKMSLAY)
C
C
C---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
C
         CALL ZCOPY(NKMSLAY*NKMSLAY,GDI(1,1,N),1,E(1,1),1)
         CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,DMAT(1,1,N),1)
C
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,GDOW(1,1,N-1),
     &              NKMSLAY,F(1,1),NKMSLAY,C1,E(1,1),NKMSLAY)
C
         CALL ZCOPY(NKMSLAY*NKMSLAY,E(1,1),1,DINVER(1,1,N),1)
C
         CALL ZGETRF(NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,INFO)
         CALL ZGETRS('N',NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,DMAT(1,1,N)
     &               ,NKMSLAY,INFO)
C
      END DO
C
C     At this point the matrix DMAT(ndim,ndim,NPLAY) contains the
C     matrices [of dimension (ndim,ndim)]  D^n, n=1,..,NPLAY
C
C---> calculate Z_n for 1 =< n <= n-1
C
      DO N = NPLAY,1,(-1)
C
         IF ( N.EQ.NPLAY ) THEN
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,DMAT(1,1,NPLAY),1,E(1,1),1)
C
            CALL TBBTOM(NPLAY,NPLAY,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,DMAT(1,1,NPLAY),1,GDI(1,1,NPLAY),
     &                 1)
         ELSE
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,GDOW(1,1,N),
     &                 NKMSLAY,DMAT(1,1,N),NKMSLAY,C0,F(1,1),NKMSLAY)
C
            CALL TBBOFM(N+1,N+1,E,NKMSLAY,GIN,NKKR_TB)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,E(1,1),
     &                 NKMSLAY,F(1,1),NKMSLAY,C0,G(1,1),NKMSLAY)
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,GUP(1,1,N),
     &                 NKMSLAY,G(1,1),NKMSLAY,C0,F(1,1),NKMSLAY)
C
            DO LM = 1,NKMSLAY
               F(LM,LM) = C1 + F(LM,LM)
            END DO
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,DMAT(1,1,N),
     &                 NKMSLAY,F(1,1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
            CALL TBBTOM(N,N,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,E(1,1),1,GDI(1,1,N),1)
C
         END IF
C
C     here start the two loops on the row index,
C     in order to span all the matrix
C     and to calculate just the blocks that are needed
C     for the construction of the cluster of green's function
C
         IF ( ICHECK(N,N).NE.0 ) THEN
C
            IF ( N.NE.1 ) THEN
C
C     this is the loop for element G_ij with i<j
C
               DO IROW = (N-1),1,(-1)
C
                  IF ( ICHECK(IROW,N).EQ.1 ) THEN
C
                     CALL TBBOFM(IROW+1,N,E,NKMSLAY,GIN,NKKR_TB)
C
                     CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                          GUP(1,1,IROW),NKMSLAY,E(1,1),NKMSLAY,C0,
     &                          F(1,1),NKMSLAY)
C
                     CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                          DMAT(1,1,IROW),NKMSLAY,F(1,1),NKMSLAY,
     &                          C0,E(1,1),NKMSLAY)
C
                     CALL TBBTOM(IROW,N,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
                  END IF
C
C
               END DO
            END IF
C
C
            IF ( N.NE.NPLAY ) THEN
C
C     this is the loop for element G_ij with i>j
C
               DO IROW = N + 1,NPLAY,1
C
                  IF ( ICHECK(IROW,N).EQ.1 ) THEN
C
                     CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,E(1,1),1)
C
                     CALL TBBOFM(IROW,IROW,F,NKMSLAY,GIN,NKKR_TB)
C
                     CALL ZGETRF(NKMSLAY,NKMSLAY,F(1,1),NKMSLAY,IPVT,
     &                           INFO)
                     CALL ZGETRS('N',NKMSLAY,NKMSLAY,F(1,1),NKMSLAY,
     &                           IPVT,E(1,1),NKMSLAY,INFO)
C
                     DO I = 1,NKMSLAY
                        DO J = 1,NKMSLAY
                           F(I,J) = GDIOLD(I,J,IROW)
     &                              - (DINVER(I,J,IROW)-E(I,J))
                        END DO
                     END DO
C
                     CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,E(1,1),1)
C
                     CALL ZGETRF(NKMSLAY,NKMSLAY,F(1,1),NKMSLAY,IPVT,
     &                           INFO)
                     CALL ZGETRS('N',NKMSLAY,NKMSLAY,F(1,1),NKMSLAY,
     &                           IPVT,E(1,1),NKMSLAY,INFO)
C
                     CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                          E(1,1),NKMSLAY,GDOW(1,1,IROW-1),NKMSLAY,
     &                          C0,F(1,1),NKMSLAY)
C
                     CALL TBBOFM(IROW-1,N,E,NKMSLAY,GIN,NKKR_TB)
C
C
                     CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                          F(1,1),NKMSLAY,E(1,1),NKMSLAY,C0,G(1,1),
     &                          NKMSLAY)
C
                     CALL TBBTOM(IROW,N,G,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
C
                  END IF
C
C
               END DO
            END IF
         END IF
C
      END DO
C
      END
C*==tbinvsupercell.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBINVSUPERCELL(M2,M1,M3,GIN,ICHECK)
C   ********************************************************************
C   *                                                                  *
C   *  ---> ALGORITHM FOR SUPERCELL GEOMETRY                           *
C   *                                                                  *
C   *  ---> factorization D ^-1 = (prod L) * M * (prod U)              *
C   *                                                                  *
C   *       see notes R. Zeller                                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TB,ONLY:NKMSLAY,NKKR_TB,NPLAY
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 GIN(NKKR_TB,NKKR_TB),M1(NKMSLAY,NKMSLAY,NPLAY),
     &           M2(NKMSLAY,NKMSLAY,NPLAY),M3(NKMSLAY,NKMSLAY,NPLAY)
      INTEGER ICHECK(NPLAY,NPLAY)
C
C Local variables
C
      COMPLEX*16 A(:,:),B(:,:,:),C(:,:,:),CUNIT(:,:),D(:,:,:),E(:,:),
     &           F(:,:),G(:,:)
      INTEGER ICOL,INFO,IPVT(:),IROW,LM,N,NL
      SAVE A,B,C,CUNIT,D,E,F,G,IPVT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE A,B,C,D,E,F,G,CUNIT,IPVT
      IF ( .NOT.ALLOCATED(A) ) THEN
         ALLOCATE (A(NKMSLAY,NKMSLAY),B(NKMSLAY,NKMSLAY,NPLAY))
         ALLOCATE (C(NKMSLAY,NKMSLAY,NPLAY),CUNIT(NKMSLAY,NKMSLAY))
         ALLOCATE (D(NKMSLAY,NKMSLAY,NPLAY),E(NKMSLAY,NKMSLAY))
         ALLOCATE (F(NKMSLAY,NKMSLAY),G(NKMSLAY,NKMSLAY))
         ALLOCATE (IPVT(NKMSLAY))
      END IF
C ---> START OF THE FACTORIZATION L * M * U
C
      CALL CINIT(NKMSLAY*NKMSLAY,A)
      CALL CINIT(NKMSLAY*NKMSLAY*NPLAY,B)
      CALL CINIT(NKMSLAY*NKMSLAY*NPLAY,C)
      CALL CINIT(NKMSLAY*NKMSLAY*NPLAY,D)
      CALL CINIT(NKMSLAY*NKMSLAY,E)
      CALL CINIT(NKMSLAY*NKMSLAY,F)
      CALL CINIT(NKMSLAY*NKMSLAY,G)
      CALL CINIT(NKMSLAY*NKMSLAY,CUNIT)
C
C ----------------------------------------------------------------------
C
C ---> cunit = complex unity matrix of order NDIM
C
      DO N = 1,NKMSLAY
         CUNIT(N,N) = C1
      END DO
C
C
      CALL ZCOPY(NKMSLAY*NKMSLAY,M2(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT,1,D(1,1,1),1)
      CALL ZGETRF(NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,INFO)
      CALL ZGETRS('N',NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,D(1,1,1),
     &            NKMSLAY,INFO)
C
C
      NL = NPLAY
C
C
      IF ( NL.NE.1 ) THEN
C
         CALL ZCOPY(NKMSLAY*NKMSLAY,M2(1,1,NL),1,A(1,1),1)
         CALL ZCOPY(NKMSLAY*NKMSLAY,M1(1,1,NL),1,B(1,1,1),1)
         CALL ZCOPY(NKMSLAY*NKMSLAY,M3(1,1,NL),1,C(1,1,1),1)
C
C
C ----------------------------------------------------------------------
C
C ---> 2 <= N < NL-1
C
         IF ( NL.NE.2 ) THEN
C
            IF ( NL.NE.3 ) THEN
C
               DO N = 2,NL - 2
C
C
C ---> E = D(N-1) * C(N-1)
C
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       D(1,1,N-1),NKMSLAY,C(1,1,N-1),NKMSLAY,C0,
     &                       E(1,1),NKMSLAY)
C
C
C ---> F = D(N-1) * M1(N-1)
C
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       D(1,1,N-1),NKMSLAY,M1(1,1,N-1),NKMSLAY,C0,
     &                       F(1,1),NKMSLAY)
C
C
C ---> A = A - B(N-1)*D(N-1)*C(N-1)
C
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       B(1,1,N-1),NKMSLAY,E(1,1),NKMSLAY,C1,A(1,1)
     &                       ,NKMSLAY)
C
C
C ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       B(1,1,N-1),NKMSLAY,F(1,1),NKMSLAY,C0,
     &                       B(1,1,N),NKMSLAY)
C
C
C ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)
C
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       M3(1,1,N-1),NKMSLAY,E(1,1),NKMSLAY,C0,
     &                       C(1,1,N),NKMSLAY)
C
C
C ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
C
C
                  CALL ZCOPY(NKMSLAY*NKMSLAY,M2(1,1,N),1,E(1,1),1)
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       M3(1,1,N-1),NKMSLAY,F(1,1),NKMSLAY,C1,
     &                       E(1,1),NKMSLAY)
                  CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,D(1,1,N),1)
                  CALL ZGETRF(NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,INFO)
                  CALL ZGETRS('N',NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,
     &                        D(1,1,N),NKMSLAY,INFO)
C
C
               END DO
            END IF
C
C ----------------------------------------------------------------------
C
C ---> N = NL - 1
C
C
C
            N = NL - 1
C
C
C ---> E = D(N-1) * C(N-1)
C
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,N-1),
     &                 NKMSLAY,C(1,1,N-1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
C
C ---> F = D(N-1) * M1(N-1)
C
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,N-1),
     &                 NKMSLAY,M1(1,1,N-1),NKMSLAY,C0,F(1,1),NKMSLAY)
C
C
C ---> A = A - B(N-1)*D(N-1)*C(N-1)
C
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,B(1,1,N-1),
     &                 NKMSLAY,E(1,1),NKMSLAY,C1,A(1,1),NKMSLAY)
C
C
C ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)
C
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,M3(1,1,N),1,B(1,1,N),1)
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,B(1,1,N-1),
     &                 NKMSLAY,F(1,1),NKMSLAY,C1,B(1,1,N),NKMSLAY)
C
C
C ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)
C
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,M1(1,1,N),1,C(1,1,N),1)
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,M3(1,1,N-1),
     &                 NKMSLAY,E(1,1),NKMSLAY,C1,C(1,1,N),NKMSLAY)
C
C
C ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
C
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,M2(1,1,N),1,E(1,1),1)
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,M3(1,1,N-1),
     &                 NKMSLAY,F(1,1),NKMSLAY,C1,E(1,1),NKMSLAY)
            CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,D(1,1,N),1)
            CALL ZGETRF(NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,INFO)
            CALL ZGETRS('N',NKMSLAY,NKMSLAY,E(1,1),NKMSLAY,IPVT,D(1,1,N)
     &                  ,NKMSLAY,INFO)
         END IF
C
C
C ----------------------------------------------------------------------
C
C ---> N = NL
C
C
         N = NL
C
C
C ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1
C
C
C ---> E = D(NL-1) * C(NL-1)
C
C
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,NL-1),
     &              NKMSLAY,C(1,1,NL-1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
C
C ---> A = A - B(NL-1) * E
C
C
         CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,B(1,1,NL-1),
     &              NKMSLAY,E(1,1),NKMSLAY,C1,A(1,1),NKMSLAY)
C
C
C ---> D(NL) = (A)^-1
C
C
         CALL ZCOPY(NKMSLAY*NKMSLAY,CUNIT(1,1),1,D(1,1,NL),1)
         CALL ZGETRF(NKMSLAY,NKMSLAY,A(1,1),NKMSLAY,IPVT,INFO)
         CALL ZGETRS('N',NKMSLAY,NKMSLAY,A(1,1),NKMSLAY,IPVT,D(1,1,NL),
     &               NKMSLAY,INFO)
      END IF
C
C ----------------------------------------------------------------------
C
C --->  END OF FACTORIZATION
C
C ----------------------------------------------------------------------
C
C
C ---> HERE THE LOOP OVER THE DIAGONAL ELEMENTS STARTS.
C ---> THE PROGRAM CHECKS IF ICHECK(N,N) = 1 AND, IF SO,
C ---> IT CALCULATES THE DIAGONAL BLOCK (N,N)
C ---> THEN IT MAKES TWO LOOPS, ONE OVER A ROW INDEX `IROW`
C ---> AND THE OTHER OVER A COLUMN INDEX `ICOL`, AND CHECKS
C ---> WHICH ARE THE ELEMENTS THAT HAS TO BE CALCULATED.
C ---> (THE ONES FOR WHICH ICHECK = 1)
C
C
C ---> IT STARTS THE LOOP OVER N
C
C
      DO N = NL,1,(-1)          ! START OF THE LOOP OVER THE DIAGONAL
C
         IF ( N.EQ.NL ) THEN
C
C ---> GTOT(NL,NL) = D(NL)
C
            CALL ZCOPY(NKMSLAY*NKMSLAY,D(1,1,NL),1,E(1,1),1)
C
            CALL TBBTOM(NL,NL,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
C
         ELSE IF ( N.EQ.NL-1 ) THEN
C
            IF ( ICHECK(NL-1,NL-1).EQ.1 ) THEN
C
C ---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
C
C
               CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,B(1,1,NL-1)
     &                    ,NKMSLAY,D(1,1,NL-1),NKMSLAY,C0,E(1,1),
     &                    NKMSLAY)
C
               CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,NL),
     &                    NKMSLAY,E(1,1),NKMSLAY,C0,F(1,1),NKMSLAY)
C
               CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,C(1,1,NL-1)
     &                    ,NKMSLAY,F(1,1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
               DO LM = 1,NKMSLAY
                  E(LM,LM) = C1 + E(LM,LM)
               END DO
C
               CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,NL-1)
     &                    ,NKMSLAY,E(1,1),NKMSLAY,C0,F(1,1),NKMSLAY)
C
               CALL TBBTOM(NL-1,NL-1,F,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
            END IF
C
C
         ELSE IF ( ICHECK(N,N).EQ.1 ) THEN
C
C
C ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
C                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
C                  - Z(N,NL)*B(N)*D(N)
C
            CALL TBBOFM(NL,N+1,F,NKMSLAY,GIN,NKKR_TB)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,C(1,1,N),
     &                 NKMSLAY,F(1,1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
            CALL TBBOFM(N+1,N+1,F,NKMSLAY,GIN,NKKR_TB)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,M1(1,1,N),
     &                 NKMSLAY,F(1,1),NKMSLAY,C1,E(1,1),NKMSLAY)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,M3(1,1,N),
     &                 NKMSLAY,D(1,1,N),NKMSLAY,C0,F(1,1),NKMSLAY)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,E(1,1),
     &                 NKMSLAY,F(1,1),NKMSLAY,C0,G(1,1),NKMSLAY)
C
            DO LM = 1,NKMSLAY
               G(LM,LM) = C1 + G(LM,LM)
            END DO
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,D(1,1,N),
     &                 NKMSLAY,G(1,1),NKMSLAY,C0,E(1,1),NKMSLAY)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,B(1,1,N),
     &                 NKMSLAY,D(1,1,N),NKMSLAY,C0,F(1,1),NKMSLAY)
C
            CALL TBBOFM(N,NL,G,NKMSLAY,GIN,NKKR_TB)
C
            CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,G(1,1),
     &                 NKMSLAY,F(1,1),NKMSLAY,C1,E(1,1),NKMSLAY)
C
            CALL TBBTOM(N,N,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
         END IF
C
C
         DO IROW = (N-1),1,(-1) ! LOOP OVER THE ROW FOR THE COLUMN N
C
            IF ( ICHECK(IROW,N).EQ.1 ) THEN
C
               IF ( (N.EQ.NL) .AND. (IROW.EQ.(NL-1)) ) THEN
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       D(1,1,NL-1),NKMSLAY,C(1,1,NL-1),NKMSLAY,C0,
     &                       F(1,1),NKMSLAY)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,F(1,1),
     &                       NKMSLAY,D(1,1,NL),NKMSLAY,C0,G(1,1),
     &                       NKMSLAY)
C
                  CALL TBBTOM(NL-1,NL,G,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
               ELSE
C
C     M(I,I+1) * Z(I+1,J)
C
                  CALL TBBOFM(IROW+1,N,E,NKMSLAY,GIN,NKKR_TB)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       M1(1,1,IROW),NKMSLAY,E(1,1),NKMSLAY,C0,
     &                       G(1,1),NKMSLAY)
C
C     M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)
C
                  CALL TBBOFM(NL,N,E,NKMSLAY,GIN,NKKR_TB)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       C(1,1,IROW),NKMSLAY,E(1,1),NKMSLAY,C1,
     &                       G(1,1),NKMSLAY)
C
C     -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       D(1,1,IROW),NKMSLAY,G(1,1),NKMSLAY,C0,
     &                       E(1,1),NKMSLAY)
C
                  CALL TBBTOM(IROW,N,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
               END IF
C
            END IF
C
         END DO                 ! LOOP OVER THE ROW FOR THE COLUMN N
C
C
         DO ICOL = (N-1),1,(-1) ! LOOP OVER THE COLUMN FOR THE ROW N
C
            IF ( ICHECK(N,ICOL).EQ.1 ) THEN
C
               IF ( (N.EQ.NL) .AND. (ICOL.EQ.(NL-1)) ) THEN
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,
     &                       B(1,1,NL-1),NKMSLAY,D(1,1,NL-1),NKMSLAY,C0,
     &                       E(1,1),NKMSLAY)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,
     &                       D(1,1,NL),NKMSLAY,E(1,1),NKMSLAY,C0,G(1,1),
     &                       NKMSLAY)
C
                  CALL TBBTOM(NL,NL-1,G,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
               ELSE
C
C     Z(I,J+1) * M(J+1,J)
C
                  CALL TBBOFM(N,ICOL+1,E,NKMSLAY,GIN,NKKR_TB)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,E(1,1),
     &                       NKMSLAY,M3(1,1,ICOL),NKMSLAY,C0,G(1,1),
     &                       NKMSLAY)
C
C     Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)
C
                  CALL TBBOFM(N,NL,E,NKMSLAY,GIN,NKKR_TB)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,C1,E(1,1),
     &                       NKMSLAY,B(1,1,ICOL),NKMSLAY,C1,G(1,1),
     &                       NKMSLAY)
C
C     -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)
C
                  CALL ZGEMM('N','N',NKMSLAY,NKMSLAY,NKMSLAY,-C1,G(1,1),
     &                       NKMSLAY,D(1,1,ICOL),NKMSLAY,C0,E(1,1),
     &                       NKMSLAY)
C
                  CALL TBBTOM(N,ICOL,E,NKMSLAY,GIN,NKKR_TB,.FALSE.)
C
               END IF
C
            END IF
C
         END DO                 ! LOOP OVER THE COLUMN FOR THE ROW N
C
C
      END DO                    ! END OF THE LOOP OVER THE DIAGONAL
C
      END
C*==tbbtom.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBBTOM(PL1,PL2,BLOCK,NSIZE,GIN,NKKR_TB,LSUB)
C   ********************************************************************
C   *                                                                  *
C   *   This subroutine copies or subtracts a block to a matrix        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL LSUB
      INTEGER NKKR_TB,NSIZE,PL1,PL2
      COMPLEX*16 BLOCK(NSIZE,NSIZE),GIN(NKKR_TB,NKKR_TB)
C
C Local variables
C
      INTEGER I1,I1S,I2,I2S
C
C*** End of declarations rewritten by SPAG
C
C
      I1S = (PL1-1)*NSIZE
      I2S = (PL2-1)*NSIZE
      IF ( LSUB ) THEN
         DO I1 = 1,NSIZE
            DO I2 = 1,NSIZE
               GIN(I1S+I1,I2S+I2) = GIN(I1S+I1,I2S+I2) - BLOCK(I1,I2)
            END DO
         END DO
      ELSE
         DO I1 = 1,NSIZE
            DO I2 = 1,NSIZE
               GIN(I1S+I1,I2S+I2) = BLOCK(I1,I2)
            END DO
         END DO
      END IF
C
      END
