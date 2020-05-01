C*==sigkloops_tb.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE SIGKLOOPS_TB(TAUQX,TSSQX,DSSQX,USSQX,VSSQX,DSSQX_L,
     &                        DSSQX_R,GREF_I1,ICLU_REF_IQTB,WA,WB,
     &                        TAUQBZ)
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
C   *  NOTE: the routine calculates the TAU_Delta according to the     *
C   *        the TB-scheme. WITHIN the k-loop  TAU_Delta(k) is         *
C   *        converted to TAU(k)                                       *
C   *                                                                  *
C   *  NOTE: the k-mesh without symmetry restrictions is used          *
C   *        if  USE_FULL_SYMMETRY=.FALSE.                             *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:CHIZ,NZ12,NZ12MAX,ITTA,ITTB,ITTC,ITTD,ITTQ1,
     &    ITTQ2,NTTJ,JTT1,JTT2,JTTX,WTTJ,NTKTKLIN
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,KTABTAUIJ
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:C0,C1,CI2PI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NLMMAX,NKM,NLM,NXM,NKMMAX,NKMQ
      USE MOD_SITES,ONLY:NQ_L,NQ_R,NQTB,IQBOT_TB,NQMAX,IQBOT_CHI,
     &    IQTOP_CHI,NQ_CHI,QBAS
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF,IQTBORGQTBP
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK,NPLAY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGKLOOPS_TB')
      LOGICAL USE_FULL_SYMMETRY
      PARAMETER (USE_FULL_SYMMETRY=.TRUE.)
C
C Dummy arguments
C
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_L(NXM,NXM,NQ_L),
     &           DSSQX_R(NXM,NXM,NQ_R),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           TAUQBZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),TAUQX(NXM,NXM,NQTB)
     &           ,TSSQX(NXM,NXM,NQTB),USSQX(NXM,NXM,NQTB),
     &           VSSQX(NXM,NXM,NQTB),WA(NXM,NXM),WB(NXM,NXM)
      INTEGER ICLU_REF_IQTB(NQTB)
C
C Local variables
C
      COMPLEX*16 AUXQ1(:,:,:,:),AUXQ2(:,:,:),CSCL,CSUM1,CSUMX,CWORK(:,:)
     &           ,EXIKDQ_QTB(:,:),GLLKE(:,:),GREFKE(:,:),TAUKLIN(:),
     &           TKTKQQ(:,:),WT
      REAL*8 DDOT
      REAL*8 DQ_QTB(:,:,:),IM_KVEC(3),KVEC(3),WK,WKSUM
      INTEGER I,IA_ERR,IK,IKTOP,IPROCK(:),IQ,IQCHI,IQTB,IU,IU0,IV,IV0,
     &        IX,IX0,J,JQ,JQCHI,JQTB,JTT,JU,JU0,JV,JV0,JX,K1,K2,L1,L2,
     &        L3,L4,M,N1,N2,NKMSQ,Z2
C
C*** End of declarations rewritten by SPAG
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE DQ_QTB,EXIKDQ_QTB
      ALLOCATABLE AUXQ1,AUXQ2
      ALLOCATABLE GLLKE,GREFKE,TAUKLIN,TKTKQQ,IPROCK,CWORK
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (AUXQ2(NXM,NXM,NQTB),AUXQ1(NXM,NXM,NQTB,2))
      ALLOCATE (DQ_QTB(3,NQTB,NQTB),EXIKDQ_QTB(NQTB,NQTB))
C
C----------------------------------------------------------------------
C     list of inter site vectors      DQ[IQTB,JQTB] = q[JQ] - q[IQ]
C     used for phase factor   EXIKDQ_QTB(IQTB,JQTB)
C----------------------------------------------------------------------
      DO JQTB = 1,NQTB
         JQ = IQBOT_TB - 1 + JQTB
         DO IQTB = 1,NQTB
            IQ = IQBOT_TB - 1 + IQTB
            DQ_QTB(1:3,IQTB,JQTB) = QBAS(1:3,JQ) - QBAS(1:3,IQ)
         END DO
      END DO
C
C----------------------------------------------------------------------
C
      ALLOCATE (TKTKQQ(NTKTKLIN,NZ12))
      TKTKQQ(:,:) = C0
C
      CHIZ(:,:,:) = C0
C
      ALLOCATE (TAUKLIN(NKKR_TB*NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUKLIN')
      TAUKLIN = C0
C
      ALLOCATE (GLLKE(NKKR_TB,NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GLLKE')
      IF ( IREL.GT.2 ) THEN
         ALLOCATE (GREFKE(NKKRNR_TB,NKKRNR_TB),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: GREFKE')
      ELSE IF ( NKKR_TB.NE.NKKRNR_TB ) THEN
         CALL STOP_MESSAGE(ROUTINE,'NKKR_TB <> NKKRNR_TB')
      END IF
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
      CALL CINIT(NXM*NXM*NQTB,TAUQX)
C
      IF ( USE_FULL_SYMMETRY ) THEN
         IKTOP = NKTAB
         WKSUM = SUM(WKTAB(1:NKTAB))
      ELSE
         IKTOP = NKTABTAUIJ
         WKSUM = DBLE(IKTOP)
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., IKTOP
C      over processors;   IK=IKTOP  dealt with by IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCK(NKTAB))
      CALL MPI_DISTRIBUTE(IPROCK,NKTAB,MPI_KLOOP,'?')
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
               WK = 1D0/WKSUM
               KVEC(1:3) = KTABTAUIJ(1:3,IK)
            END IF
C
C----------------------------------------------------------------------
C   phase factor   EXIKDQ_QTB(IQTB,JQTB) = exp{-2*PI*k*(q[JQ] - q[IQ])}
C   see NOTE in <SIGKLOOPSTRANS_TB>
C----------------------------------------------------------------------
            DO JQTB = 1,NQTB
               DO IQTB = 1,JQTB - 1
                  EXIKDQ_QTB(IQTB,JQTB)
     &               = -CDEXP(-CI2PI*DDOT(3,KVEC,1,DQ_QTB(1,IQTB,JQTB),
     &               1))
                  EXIKDQ_QTB(JQTB,IQTB) = DCONJG(EXIKDQ_QTB(IQTB,JQTB))
               END DO
            END DO
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
C           get the site-diagonal blocks and convert to TAU(k)
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
                     TAUQX(I,J,IQTB) = TAUQX(I,J,IQTB) - WK*GLLKE(IX,JX)
     &                                 /WKSUM
                  END DO
               END DO
            END DO
C
C-------------------------------------- transform TAU_Delta(k) => TAU(k)
C
            CALL SIGKLOOPSTRANS_TB(GLLKE,TAUKLIN,USSQX,VSSQX,WA,WB,
     &                             NKKR_TB,EXIKDQ_QTB)
C
C------------------------------------------------------- store TAUA*TAUB
C
            IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
               CALL SIGKLOOPSCHI(WK,WKSUM,TAUKLIN,CHIZ,NKKR_TB,NZ12MAX)
C
            ELSE
C
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  CSUMX = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
                     CSUM1 = CSUM1 + WT*TAUKLIN(JTT2(JTT))
                     CSUMX = CSUMX + WT*DCONJG(TAUKLIN(JTTX(JTT)))
                  END DO
                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                  TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
               END DO
C-----------------------------------------------------------------------
C
            END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_KLOOP ) THEN
C
         CALL DRV_MPI_BARRIER
C
C        use TAUQBZ as work space for transfer
         M = NXM*NXM*NQTB
         CALL DRV_MPI_REDUCE_C(TAUQX(1,1,1),TAUQBZ(1,1,1,1),M)
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
            M = NKMMAX*NKMMAX*NQ_CHI*NKMMAX*NKMMAX*NQ_CHI
            ALLOCATE (CWORK(M,1))
            CALL DRV_MPI_REDUCE_C(CHIZ(1,1,1),CWORK(1,1),M)
            CALL DRV_MPI_REDUCE_C(CHIZ(1,1,2),CWORK(1,1),M)
            DEALLOCATE (CWORK)
         ELSE
            M = NTKTKLIN*NZ12
            ALLOCATE (CWORK(NTKTKLIN,NZ12))
            CALL DRV_MPI_REDUCE_C(TKTKQQ(1,1),CWORK(1,1),M)
            DEALLOCATE (CWORK)
         END IF
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C
C-----------------------------------------------------------------------
C    set up the auxilary matrices to correct the site-diagonal case
C-----------------------------------------------------------------------
C
         DO IQTB = 1,NQTB
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TAUQX(1,1,IQTB),NXM,
     &                 VSSQX(1,1,IQTB),NXM,C0,WB,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX(1,1,IQTB),NXM,WB,
     &                 NXM,C0,AUXQ1(1,1,IQTB,1),NXM)
C
         END DO
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
            AUXQ2(1:NXM,1:NXM,1:NQTB) = AUXQ1(1:NXM,1:NXM,1:NQTB,1)
C
         ELSE
C
            CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,AUXQ1(1,1,1,1),
     &                     AUXQ2(1,1,1),WA,NQTB,NKMQ(IQBOT_TB),DROT,
     &                     IQTBORGQTBP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                     NSYMACCEPTED,NQTB,NKMMAX)
C
         END IF
C
         DO IQTB = 1,NQTB
C
            DO J = 1,NXM
               DO I = 1,NXM
                  AUXQ1(I,J,IQTB,1) = AUXQ2(I,J,IQTB)
                  AUXQ1(I,J,IQTB,2) = DCONJG(AUXQ2(J,I,IQTB))
               END DO
            END DO
C
         END DO
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
     &                 DSSQX(1,1,IQTB),NXM,C0,WB,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WB,
     &                 NXM,C0,WA,NXM)
C
            WB(1:NXM,1:NXM) = WA(1:NXM,1:NXM) - DSSQX(1:NXM,1:NXM,IQTB)
C
C--------------------------------------------------- TAU = t * G * t + t
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WB,NXM,TSSQX(1,1,IQTB),
     &                 NXM,C0,WA,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TSSQX(1,1,IQTB),NXM,WA,
     &                 NXM,C0,WB,NXM)
C
            DO J = 1,NXM
               DO I = 1,NXM
                  TAUQX(I,J,IQTB) = WB(I,J) + TSSQX(I,J,IQTB)
               END DO
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C             change representation amd symmetrize
C-----------------------------------------------------------------------
         IF ( IREL.EQ.3 ) THEN
            DO IQTB = 1,NQTB
               WA(1:NKM,1:NKM) = TAUQX(1:NKM,1:NKM,IQTB)
C
               CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
               TAUQX(1:NKM,1:NKM,IQTB) = WB(1:NKM,1:NKM)
            END DO
         END IF
C
         CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,TAUQX(1,1,1),
     &                  AUXQ2(1,1,1),WA,NQTB,NKMQ(IQBOT_TB),DROT,
     &                  IQTBORGQTBP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                  NSYMACCEPTED,NQTB,NKMMAX)
C
C-----------------------------------------------------------------------
C                   set up the matrix   TAUQBZ
C-----------------------------------------------------------------------
C
         DO IQTB = 1,NQTB
C
            IQ = IQBOT_TB - 1 + IQTB
            DO J = 1,NXM
               DO I = 1,NXM
                  TAUQBZ(I,J,IQ,1) = AUXQ2(I,J,IQTB)
                  TAUQBZ(I,J,IQ,2) = DCONJG(AUXQ2(J,I,IQTB))
               END DO
            END DO
C
         END DO
C
C=======================================================================
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
            NKMSQ = NKM*NKM
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               IQCHI = IQ - IQBOT_CHI + 1
C
               DO JQ = IQBOT_CHI,IQTOP_CHI
                  JQCHI = JQ - IQBOT_CHI + 1
C
                  K1 = (IQCHI-1)*NKMSQ
                  DO L1 = 1,NKM
                     DO L4 = 1,NKM
                        K1 = K1 + 1
C
                        K2 = (JQCHI-1)*NKMSQ
                        DO L2 = 1,NKM
                           DO L3 = 1,NKM
                              K2 = K2 + 1
C
                              IF ( IQ.EQ.JQ ) THEN
C
                                 CHIZ(K1,K2,1) = CHIZ(K1,K2,1)
     &                              - AUXQ1(L1,L2,IQCHI,1)
     &                              *AUXQ1(L3,L4,IQCHI,1)
C
                                 CHIZ(K1,K2,2) = CHIZ(K1,K2,2)
     &                              - AUXQ1(L1,L2,IQCHI,1)
     &                              *AUXQ1(L3,L4,IQCHI,2)
C
                              END IF
C
                           END DO
                        END DO
                     END DO
                  END DO
C
               END DO
            END DO
C
         ELSE
C=======================================================================
C
            CSCL = 1D0/WKSUM
            DO Z2 = 1,NZ12
               CALL ZSCAL(NTKTKLIN,CSCL,TKTKQQ(1,Z2),1)
            END DO
C
C-----------------------------------------------------------------------
C               set up the super-matrix  CHI(K1,K2)
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
            NKMSQ = NKM*NKM
C
            DO I = 1,NTKTKLIN
C
               L1 = ITTA(I)
               L2 = ITTB(I)
               L3 = ITTC(I)
               L4 = ITTD(I)
               IQ = ITTQ1(I)
               JQ = ITTQ2(I)
               IQCHI = IQ - IQBOT_CHI + 1
               JQCHI = JQ - IQBOT_CHI + 1
               IQTB = IQ - IQBOT_TB + 1
C
               K1 = (L1-1)*NKM + L4 + (IQCHI-1)*NKMSQ
               K2 = (L2-1)*NKM + L3 + (JQCHI-1)*NKMSQ
C
               DO Z2 = 1,NZ12
C
                  IF ( IQ.EQ.JQ ) THEN
                     CHIZ(K1,K2,Z2) = TKTKQQ(I,Z2) - AUXQ1(L1,L2,IQTB,1)
     &                                *AUXQ1(L3,L4,IQTB,Z2)
                  ELSE
                     CHIZ(K1,K2,Z2) = TKTKQQ(I,Z2)
                  END IF
               END DO
            END DO
C
         END IF
C
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_KLOOP ) CALL DRV_MPI_BARRIER
C
      END
