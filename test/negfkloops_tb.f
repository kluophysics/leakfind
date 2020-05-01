C*==negfkloops_tb.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOPS_TB(GLESQX,TAUQX,TSSQX,DSSQX,DSSQX_L,DSSQX_R,
     &                         GREF_I1,ICLU_REF_IQTB,WA,WB,LMAT3,GLESQ,
     &                         TAUQ,IQBOT_LBAR,IQTOP_LBAR,IQBOT_RBAR,
     &                         IQTOP_RBAR,LBAR_FLAG_Q,RBAR_FLAG_Q,
     &                         BAR_FLAG_Q,CCZ,CHECK_ZPM,MEZZ,MEZZMM,
     &                         MEZZMP,MEZZPM,TRANSMISSION,STT,DELTA)
C   ********************************************************************
C   *                                                                  *
C   *  perform the loop over k-points                                  *
C   *                                                                  *
C   *  INPUT:                                                          *
C   *         GREF_I1  G-matrix of the real space TB reference cluster *
C   *                                                                  *
C   *  OUTPUT:                                                         *
C   *         TAUQX  the TAU-matrix obtained by BZ-integration         *
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
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NLMMAX,NKM,NLM,NXM,NKMMAX,NKMQ,WKM1,WKM2,
     &    NMEMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:C0,C1,CI2PI
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,WRBUILDBOT
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB
      USE MOD_LATTICE,ONLY:VOLUC_2D
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_SITES,ONLY:NQ_L,NQ_R,NQTB,IQBOT_TB,IQTOP_TB,QBAS,NQ
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,KTABTAUIJ
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK,NPLAY,NKMSLAY
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF,IQTBORGQTBP
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL USE_FULL_SYMMETRY
      PARAMETER (USE_FULL_SYMMETRY=.TRUE.)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NEGFKLOOPS_TB')
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM,STT,TRANSMISSION
      INTEGER IQBOT_LBAR,IQBOT_RBAR,IQTOP_LBAR,IQTOP_RBAR
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 DELTA(NLMMAX,NLMMAX,NTMAX),DSSQX(NXM,NXM,NQTB),
     &           DSSQX_L(NXM,NXM,NQ_L),DSSQX_R(NXM,NXM,NQ_R),
     &           GLESQ(NKMMAX,NKMMAX,NQ),GLESQX(NXM,NXM,NQTB),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),LMAT3(NXM,NXM),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZPM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQ),TAUQX(NXM,NXM,NQTB),
     &           TSSQX(NXM,NXM,NQTB),WA(NXM,NXM),WB(NXM,NXM)
      INTEGER ICLU_REF_IQTB(NQTB)
C
C Local variables
C
      COMPLEX*16 AUXQ(:,:,:),EXIKDQ_QTB(:,:),GLESQX_TR(:,:,:),GLLKE(:,:)
     &           ,GLLKE_BAR(:,:),GLLKE_BAR_TR(:,:),GREFKE(:,:),
     &           LMAT2(NXM,NXM),SIGMA_L(NKMSLAY,NKMSLAY),
     &           SIGMA_R(NKMSLAY,NKMSLAY),TINT(3),TRANS(3),WIJ(:,:),
     &           WJI(:,:)
      REAL*8 DDOT
      REAL*8 DQ_QTB(:,:,:),IM_KVEC(3),KVEC(3),WK,WKSUM
      INTEGER I,IA_ERR,IK,IKTOP,IPROCK(:),IQ,IQTB,IU,IU0,IV,IV0,IX,IX0,
     &        J,JQ,JQTB,JU,JU0,JV,JV0,JX,KM1_LASTPLAY,M,N1,N2
C
C*** End of declarations rewritten by SPAG
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE AUXQ,WIJ,WJI
      ALLOCATABLE DQ_QTB,EXIKDQ_QTB
      ALLOCATABLE GLLKE,GREFKE,GLLKE_BAR,GLLKE_BAR_TR,IPROCK
      ALLOCATABLE GLESQX_TR
C
      ALLOCATE (AUXQ(NXM,NXM,NQTB),WIJ(NXM,NXM),WJI(NXM,NXM))
      ALLOCATE (DQ_QTB(3,NQTB,NQTB),EXIKDQ_QTB(NQTB,NQTB))
      IF ( STT ) ALLOCATE (GLESQX_TR(NXM,NXM,NQTB))
C
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
      ALLOCATE (GLLKE(NKKR_TB,NKKR_TB))
      ALLOCATE (GLLKE_BAR(NKKR_TB,NKKR_TB))
      ALLOCATE (GLLKE_BAR_TR(NKKR_TB,NKKR_TB),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP '<TAUIJ_TB_KLOOP> allocate GLLKE'
      IF ( IREL.GT.2 ) THEN
         ALLOCATE (GREFKE(NKKRNR_TB,NKKRNR_TB),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP '<TAUIJ_TB_KLOOP> allocate GREFKE'
      ELSE IF ( NKKR_TB.NE.NKKRNR_TB ) THEN
         STOP '<TAUIJ_TB_KLOOP>:  NKKR_TB <> NKKRNR_TB'
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
      TAUQX(:,:,:) = C0
      GLESQX(:,:,:) = C0
      IF ( STT ) GLESQX_TR(:,:,:) = C0
      TINT(:) = C0
C
      IF ( USE_FULL_SYMMETRY ) THEN
         IKTOP = NKTAB
         WKSUM = SUM(WKTAB(1:NKTAB))
      ELSE
         IKTOP = NKTABTAUIJ
         WKSUM = DBLE(IKTOP)
      END IF
C
      CALL CHANGEREP(NKM,NKMMAX,LMAT3,'REL>RLM',LMAT2)
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
C     IKTOP = 1
C     WKSUM = 1D0
      WRITE (6,*) 'VOLUC_2D =',VOLUC_2D,'[ALAT^2]'
C
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
               WK = 1D0
               KVEC(1:3) = KTABTAUIJ(1:3,IK)
            END IF
C
            IF ( IKTOP.EQ.1 ) KVEC(1:3) = (/0.0D0,0.0D0,0D0/)
C
            IF ( IPRINT.GE.1 ) THEN
               WRITE (345,*) 'IK =',IK
               WRITE (345,*) 'KVEC =',KVEC
               WRITE (345,*) 'WK =',WK
            END IF
C
C----------- taking now care of the site-diagonal case & the sign of tau
            DO JQTB = 1,NQTB
               DO IQTB = 1,JQTB
                  EXIKDQ_QTB(IQTB,JQTB)
     &               = -CDEXP(-CI2PI*DDOT(3,KVEC,1,DQ_QTB(1,IQTB,JQTB),
     &               1))
                  EXIKDQ_QTB(JQTB,IQTB) = DCONJG(EXIKDQ_QTB(IQTB,JQTB))
               END DO
            END DO
C
C-----------------------------------------------------------------------
C    Fourier transformation, set KKR matrix M = [-(t)^-1 + G^r]
CSW                                              ???????
C-----------------------------------------------------------------------
C
            IF ( IREL.LE.2 ) THEN
C
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
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM_MO('GREFKE',111,GLLKE,EXIKDQ_QTB,
     &           LMAT2,CCZ)
C
C-----------------------------------------------------------------------
C              call decimation routine if requested
C              store self-energies of left & right surface SIGMA_L/R
C-----------------------------------------------------------------------
C
            IF ( IDECI.EQ.1 ) CALL TBDECIMATE_MO(GLLKE,DSSQX_L,DSSQX_R,
     &           VACFLAG,FACTL,NQ_L,NQ_R,NXM,SIGMA_L,SIGMA_R)
C
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM_MO('GREFKE_DEC',1111,GLLKE,
     &           EXIKDQ_QTB,LMAT2,CCZ)
C
CSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSW
C  copy SIGMA_Xs to GLLKE_BAR (used as work array)
C
            GLLKE_BAR(:,:) = C0
            GLLKE_BAR_TR(:,:) = C0
C
            KM1_LASTPLAY = (NKKR_TB-NKMSLAY) + 1
C
            GLLKE_BAR(1:NKMSLAY,1:NKMSLAY) = SIGMA_L(:,:)
C
            GLLKE_BAR(KM1_LASTPLAY:NKKR_TB,KM1_LASTPLAY:NKKR_TB)
     &         = SIGMA_R(:,:)
C
            IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM_MO('SIGMA_X',777,GLLKE_BAR,
     &           EXIKDQ_QTB,LMAT2,CCZ)
C
CSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSWSW
C-----------------------------------------------------------------------
C    Construct the matrix M=[-(t-t_ref)^-1 + G_ref] and store it
C    in the same matrix GLLKE where  G_ref  was stored.
C-----------------------------------------------------------------------
C
            IX0 = -NXM
            DO IQTB = 1,NQTB
               IQ = IQBOT_TB - 1 + IQTB
               IX0 = IX0 + NXM
               JX = IX0
C
               DO J = 1,NXM
                  JX = JX + 1
                  IX = IX0
                  DO I = 1,NXM
                     IX = IX + 1
                     GLLKE(IX,JX) = GLLKE(IX,JX) - DSSQX(I,J,IQTB)
                  END DO
               END DO
C
            END DO
C
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM_MO('GLLKE',112,GLLKE,EXIKDQ_QTB,
     &           LMAT2,CCZ)
C
C
C-----------------------------------------------------------------------
C     Perform the inversion of matrix M
C     the output is the scattering path operator -TAU_DELTA(k)
C-----------------------------------------------------------------------
C
C NOTE: DSSQX dummy argument - not used for INVMOD <> 4
C
            IF ( INVMOD.EQ.4 ) STOP 'INVMOD = 4 '
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
C
               IF ( CHECK_ZPM .AND. IK.EQ.IKTOP ) THEN
                  IF ( CCZ ) THEN
C
                     WKM1(:,:) = TAUQX(:,:,IQTB)
C
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM1,
     &                          NXM,C0,WKM2,NXM)
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,LMAT2,
     &                          NXM,C0,WKM1,NXM)
C
                     WRITE (91,*) 'IQTB =',IQTB
                     CALL CMATSTRUCT('TAUQX (-)',WKM1(:,:),18,18,2,2,1,
     &                               1.D-8,91)
C
                  ELSE
C
                     WRITE (91,*) 'IQTB =',IQTB
                     CALL CMATSTRUCT('TAUQX (+)',TAUQX(:,:,IQTB),18,18,
     &                               2,2,1,1.D-8,91)
                  END IF
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
C    convert self-energies for the barriers to renormalized G_surf
C-----------------------------------------------------------------------
C
C                                                           left surface
C
            CALL NEGFKLOOP_GBAR_MO(GLLKE_BAR,DSSQX,NKKR_TB,EXIKDQ_QTB,
     &                             WIJ,WJI,IQBOT_LBAR,IQTOP_LBAR,CCZ,
     &                             LMAT2,IK,IKTOP,CHECK_ZPM,MEZZ,MEZZMM,
     &                             MEZZMP,MEZZPM,TRANSMISSION,STT,
     &                             GLLKE_BAR_TR)
C
C                                                          right surface
C
            CALL NEGFKLOOP_GBAR_MO(GLLKE_BAR,DSSQX,NKKR_TB,EXIKDQ_QTB,
     &                             WIJ,WJI,IQBOT_RBAR,IQTOP_RBAR,CCZ,
     &                             LMAT2,IK,IKTOP,CHECK_ZPM,MEZZ,MEZZMM,
     &                             MEZZMP,MEZZPM,.FALSE.,STT,
     &                             GLLKE_BAR_TR)
C
C-----------------------------------------------------------------------
C                    calculate lesser Green's function   G<   for DOS
C                    or calculate Transmission for TRANMISSION and STT
C-----------------------------------------------------------------------
C
            CALL NEGFKLOOP_GLES_MO(GLESQX,GLLKE,GLLKE_BAR,NKKR_TB,
     &                             EXIKDQ_QTB,WIJ,WJI,BAR_FLAG_Q,
     &                             LBAR_FLAG_Q,RBAR_FLAG_Q,CCZ,LMAT2,IK,
     &                             IKTOP,CHECK_ZPM,WK,WKSUM,
     &                             TRANSMISSION,DSSQX,KVEC,TRANS,STT)
C
            IF ( STT ) THEN
C
C-----------------------------------------------------------------------
C                    calculate lesser Green's function   G<_tr   for STT
C-----------------------------------------------------------------------
C
               CALL NEGFKLOOP_GLES_MO(GLESQX_TR,GLLKE,GLLKE_BAR_TR,
     &                                NKKR_TB,EXIKDQ_QTB,WIJ,WJI,
     &                                BAR_FLAG_Q,LBAR_FLAG_Q,
     &                                RBAR_FLAG_Q,CCZ,LMAT2,IK,IKTOP,
     &                                CHECK_ZPM,WK,WKSUM,.FALSE.,DSSQX,
     &                                KVEC,TRANS,STT)
C
C-----------------------------------------------------------------------
C                    calculate the spin transfer torque from
C                    G<_tr(k_II,E_F) and DELTA_IQ
C-----------------------------------------------------------------------
C
               CALL NEGFKLOOP_STT(GLESQX_TR,NKKR_TB,DELTA,LMAT2,IK,
     &                            IKTOP,WK,WKSUM,KVEC,TRANS,STT)
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
C        use TAUQ as work space for transfer
         M = NXM*NXM*NQTB
         CALL DRV_MPI_REDUCE_C(TAUQX(1,1,1),TAUQ(1,1,1),M)
         CALL DRV_MPI_REDUCE_C(GLESQX(1,1,1),TAUQ(1,1,1),M)
         IF ( STT ) CALL DRV_MPI_REDUCE_C(GLESQX_TR(1,1,1),TAUQ(1,1,1),
     &        M)
         CALL DRV_MPI_REDUCE_C(TRANS(1),TINT(1),3)
C
         IF ( TRANSMISSION ) THEN
            WRITE (6,FMT='(/,18X,A,10X,A)') 'Re','Im'
            WRITE (6,FMT='(A,2F12.6)') 'T(E) total: ',TINT(1)
            WRITE (6,FMT='(A,2F12.6)') 'T(E) spin1: ',TINT(2)
            WRITE (6,FMT='(A,2F12.6)') 'T(E) spin2: ',TINT(3)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT .AND. MPI_ID.EQ.0 )
     &           WRITE (IFILBUILDBOT,99001) ROUTINE(1:LEN_TRIM(ROUTINE))
     &           ,DREAL(TINT(1:3)),DIMAG(TINT(1:3))
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
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
C
C ----------------------------------------------------------------------
C
               IF ( IPRINT.GE.1 ) THEN
                  WRITE (291,*) 'IQTB =',IQTB
                  CALL CMATSTRUCT('TAUQX (+)',TAUQX(:,:,IQTB),18,18,3,3,
     &                            1,1.D-8,291)
               END IF
C
            END DO
         END IF
C
         CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,TAUQX(1,1,1),
     &                  AUXQ(1,1,1),WA,NQTB,NKMQ(IQBOT_TB),DROT,
     &                  IQTBORGQTBP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                  NSYMACCEPTED,NQTB,NKMMAX)
C
C-----------------------------------------------------------------------
C                   set up the matrix   TAUQ
C-----------------------------------------------------------------------
C
         DO IQTB = 1,NQTB
C
            IQ = IQBOT_TB - 1 + IQTB
            DO J = 1,NXM
               DO I = 1,NXM
                  TAUQ(I,J,IQ) = AUXQ(I,J,IQTB)
               END DO
            END DO
C
            IF ( IPRINT.GE.1 ) THEN
               WRITE (391,*) 'IQTB =',IQTB
               CALL CMATSTRUCT('TAUQ (+)',TAUQ(:,:,IQ),18,18,3,3,1,
     &                         1.D-8,391)
            END IF
C
         END DO
C
C-----------------------------------------------------------------------
C                   set up the matrix   GLESQ
C-----------------------------------------------------------------------
C
         DO IQ = IQBOT_TB,IQTOP_TB
C
            IQTB = IQ - IQBOT_TB + 1
C
            CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,LMAT2,
     &                 NXM,C0,WA,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,GLESQX(1,1,IQTB),NXM,WA,
     &                 NXM,C0,WB,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WB,
     &                 NXM,C0,WA,NXM)
C
            IF ( STT ) WA(:,:) = GLESQX_TR(:,:,IQTB)
CSW-------- the DSSQXes have already been attached inside the k-loop
C
            CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
            GLESQ(1:NKM,1:NKM,IQ) = WB(1:NKM,1:NKM)
C
C ----------------------------------------------------------------------
CSW NO symmetrization of glesq!!?
C
         END DO
C
C=======================================================================
C
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_KLOOP ) CALL DRV_MPI_BARRIER
C
99001 FORMAT ('# BUILDBOT: ',A,':  T(E) total, up & dn (3xRe, 3xIm)',/,
     &        (1PE22.14))
      END
C*==negfkloop_gbar_mo.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOP_GBAR_MO(TAUDELTAK_BAR,DSSQX,NKKR_TB,
     &                             EXIKDQ_QTB,WIJ,WJI,IQBOT_XBAR,
     &                             IQTOP_XBAR,CCZ,LMAT2,IK,IKTOP,
     &                             CHECK_ZPM,MEZZ,MEZZMM,MEZZMP,MEZZPM,
     &                             TRANSMISSION,STT,TAUDELTAK_BAR_TR)
C   ********************************************************************
C   *                                                                  *
C   *  Calculate the broadening function gamma from the self-energies  *
C   *  of the left and right leads (surface Green's function approach  *
C   *  by M. Ogura)                                                    *
C   *                                                                  *
C   *       sigma(i,j;+) c(j;+,-) - [sigma(j,i;+) c(i;+,-)]^(T*)       *
C   *                                                                  *
C   *  or for transmission                                             *
C   *                                                                  *
C   *       sigma_L(i,j;-) d(j;-,+) - [sigma_L(j,i;-) d(i;-,+)]^(T*)   *
C   *                                                                  *
C   *  - TAUDELTAK_BAR in:     self-energies of left & right surface   *
C   *                 out:     renormalized gamma                      *
C   *    TAUDELTAK_BAR_TR:     used for STT only, when two different   *
C   *                          self energies/gammas are needed         *
C   *  - called for X = L (left) and R (right) barrier separately      *
C   *  - the phase factor exp(ik Delta R) is added                     *
C   *                                                                  *
C   *  index IZ = 1 --> E(+)  = 2 --> E(-)                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NXM,WKM1,WKM2,WKM3,IPIVKM,NMEMAX
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SITES,ONLY:NQTB,IQBOT_TB
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CONLY,PLSONLY,MINONLY
      PARAMETER (CONLY=.FALSE.,PLSONLY=.FALSE.,MINONLY=.FALSE.)
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM,STT,TRANSMISSION
      INTEGER IK,IKTOP,IQBOT_XBAR,IQTOP_XBAR,NKKR_TB
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),EXIKDQ_QTB(NQTB,NQTB),
     &           LMAT2(NXM,NXM),MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZPM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TAUDELTAK_BAR(NKKR_TB,NKKR_TB),
     &           TAUDELTAK_BAR_TR(NKKR_TB,NKKR_TB),WIJ(NXM,NXM),
     &           WJI(NXM,NXM)
C
C Local variables
C
      COMPLEX*16 DTSSQX(:,:,:),GAMMAK(:,:),GAMMAKT(:,:),MEZZI2(:,:),
     &           MEZZIINV(:,:),MEZZIMP2(:,:),MEZZIPM2(:,:),MEZZJ2(:,:),
     &           MEZZJINV(:,:),MEZZJMP2(:,:),MEZZJPM2(:,:),WKMT1(:,:),
     &           WKMT2(:,:),XGAMMAKX(:,:)
      INTEGER ICALL,IQ,IQTB,IX0,IX1,IX2,JQ,JQTB,JX0,JX1,JX2,L
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/1/
C
      ALLOCATABLE GAMMAK,GAMMAKT,XGAMMAKX
      ALLOCATABLE MEZZIINV,MEZZJINV,MEZZI2,MEZZJ2
      ALLOCATABLE MEZZIMP2,MEZZIPM2,MEZZJMP2,MEZZJPM2,DTSSQX,WKMT1,WKMT2
C
      ALLOCATE (GAMMAK(NXM,NXM),XGAMMAKX(NXM,NXM))
      ALLOCATE (MEZZIINV(NXM,NXM),MEZZJINV(NXM,NXM))
      ALLOCATE (MEZZI2(NXM,NXM),MEZZJ2(NXM,NXM))
      ALLOCATE (MEZZIMP2(NXM,NXM),MEZZIPM2(NXM,NXM))
      ALLOCATE (MEZZJMP2(NXM,NXM),MEZZJPM2(NXM,NXM))
      ALLOCATE (DTSSQX(NXM,NXM,NQTB))
      IF ( STT ) THEN
         ALLOCATE (GAMMAKT(NXM,NXM))
         ALLOCATE (WKMT1(NKMMAX,NKMMAX),WKMT2(NKMMAX,NKMMAX))
      END IF
C
C*** End of declarations rewritten by SPAG
C
      IF ( IK.EQ.IKTOP ) WRITE (6,*) 'ICALL = ',ICALL
      ICALL = ICALL + 1
      IF ( ICALL.EQ.3 ) ICALL = 1
C
      IF ( PLSONLY .AND. IK.EQ.1 ) WRITE (6,*)
     &      '+++++++++++++++++    PLSONLY    +++++++++++++++++'
      IF ( MINONLY .AND. IK.EQ.1 ) WRITE (6,*)
     &      '-----------------    MINONLY    -----------------'
      IF ( CONLY .AND. IK.EQ.1 ) WRITE (6,*) 
     &               'ccccccccccccccccc     CONLY     ccccccccccccccccc'
C
C----------------------------------------------------------------------
C
      DO JQ = IQBOT_XBAR,IQTOP_XBAR
         JQTB = JQ - IQBOT_TB + 1
         JX0 = (JQTB-1)*NXM
         JX1 = JX0 + 1
         JX2 = JX0 + NXM
C
CSW ToDo: Rewrite gammas using the calculated MEZZ(-,-) matrix elements
CSW ToDo: Calculate pure gammas first and then, where necessary attach
CSW ToDo: the Xes (outside of this routine)??
         CALL CHANGEREP(NKM,NKMMAX,MEZZMP(1,1,JQ,1),'REL>RLM',MEZZJMP2)
         CALL CHANGEREP(NKM,NKMMAX,MEZZPM(1,1,JQ,1),'REL>RLM',MEZZJPM2)
C------------------------------------------------------- get (1/MEZZ(j))
         CALL CHANGEREP(NKM,NKMMAX,MEZZ(1,1,JQ,1),'REL>RLM',MEZZJ2)
         CALL CMATINV3(NXM,NKMMAX,IPIVKM,MEZZJ2,WKM1,MEZZJINV)
         CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,JQTB),WKM1,
     &                 DTSSQX(1,1,JQTB))
C
         IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP ) THEN
            WRITE (*,*) 'JQ =',JQ
            CALL CMATSTRUCT('DTSSQX',DTSSQX(:,:,JQTB),18,18,2,2,1,1.D-8,
     &                      6)
            CALL CMATSTRUCT('MEZZ',MEZZJ2(:,:),18,18,2,2,1,1.D-8,6)
            CALL CMATSTRUCT('MEZZJINV',MEZZJINV(:,:),18,18,2,2,1,1.D-8,
     &                      6)
            CALL CMATSTRUCT('MEZZMP',MEZZJMP2(:,:),18,18,2,2,1,1.D-8,6)
            CALL CMATSTRUCT('MEZZPM',MEZZJPM2(:,:),18,18,2,2,1,1.D-8,6)
            CALL CMATSTRUCT('LMAT2',LMAT2(:,:),18,18,2,2,1,1.D-8,6)
            CALL CMATSTRUCT('DSSQX',DSSQX(:,:,JQTB),18,18,2,2,1,1.D-8,6)
         END IF
C
         DO IQ = IQBOT_XBAR,JQ
            IQTB = IQ - IQBOT_TB + 1
            IX0 = (IQTB-1)*NXM
            IX1 = IX0 + 1
            IX2 = IX0 + NXM
C
            CALL CHANGEREP(NKM,NKMMAX,MEZZMP(1,1,IQ,1),'REL>RLM',
     &                     MEZZIMP2)
            CALL CHANGEREP(NKM,NKMMAX,MEZZPM(1,1,IQ,1),'REL>RLM',
     &                     MEZZIPM2)
C------------------------------------------------------- get (1/MEZZ(i))
            CALL CHANGEREP(NKM,NKMMAX,MEZZ(1,1,IQ,1),'REL>RLM',MEZZI2)
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,MEZZI2,WKM1,MEZZIINV)
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,IQTB),WKM1,
     &                    DTSSQX(1,1,IQTB))
C
CSW in TAUDELTAK_BAR the self-energies of the left and right surface
CSW are stored in the 1st and last NKMSLAY^2 blocks
C
            WIJ(1:NXM,1:NXM) = TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
     &                         *EXIKDQ_QTB(IQTB,JQTB)
            WJI(1:NXM,1:NXM) = TAUDELTAK_BAR(JX1:JX2,IX1:IX2)
     &                         *EXIKDQ_QTB(JQTB,IQTB)
C
            IF ( CONLY ) THEN
               WIJ(:,:) = C0
               WJI(:,:) = C0
               DO L = 1,NXM
                  WIJ(L,L) = C1
                  WJI(L,L) = C1
               END DO
            END IF
C
            IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP ) THEN
               WRITE (*,*) 'JQ =',JQ,'JQTB =',JQTB,'IQ =',IQ,'IQTB =',
     &                     IQTB
               CALL CMATSTRUCT('SIGMA_X(I,J)',WIJ(:,:),18,18,2,2,1,
     &                         1.D-8,6)
               WRITE (*,*) 'PHASE FACTOR = ',EXIKDQ_QTB(IQTB,JQTB)
               CALL CMATSTRUCT('SIGMA_X(J,I)',WJI(:,:),18,18,2,2,1,
     &                         1.D-8,6)
               WRITE (*,*) 'PHASE FACTOR = ',EXIKDQ_QTB(JQTB,IQTB)
            END IF
C
            IF ( TRANSMISSION ) THEN
C
               IF ( .NOT.STT .AND. IK.EQ.IKTOP ) WRITE (6,*)
     &               'SIGMA_L for Transmission'
               IF ( STT .AND. IK.EQ.IKTOP ) WRITE (6,*)
     &               'SIGMA_L for Transmission in STT'
C
C----------------------------------------------------------------------
C IJ:   L x(i,+)^(T*) L MEZZMP(i) MEZZ(i)^-1 dt(i,+) sigma_x(i,j;+) L
C     - [L sigma_x(j,i;+)^(T*) dt(j,+)^(T*) (MEZZ(j)^-1)^(T*) MEZZMP(j)^(T*)
C       L x(j,+) L]
C----------------------------------------------------------------------
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                    DSSQX(1,1,IQTB),NXM,C0,WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,C0,
     &                    WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,MEZZIMP2,NXM,
     &                    C0,WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,MEZZIINV,NXM,
     &                    C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,
     &                    DTSSQX(1,1,IQTB),NXM,C0,WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,WIJ,NXM,C0,
     &                    WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,C0,
     &                    WKM1,NXM)
C
C-----------------------------------------------------------------------
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,WJI,NXM,C0,
     &                    WKM2,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,
     &                    DTSSQX(1,1,JQTB),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,MEZZJINV,NXM,
     &                    C0,WKM2,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,MEZZJMP2,NXM,
     &                    C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,C0,
     &                    WKM2,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,
     &                    DSSQX(1,1,JQTB),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,C0,
     &                    WKM2,NXM)
C
               IF ( PLSONLY ) WKM1(:,:) = C0
               IF ( MINONLY ) WKM2(:,:) = C0
C
               IF ( STT ) THEN
                  WKMT1(:,:) = WKM1(:,:)
                  WKMT2(:,:) = WKM2(:,:)
               END IF
C
            END IF
C
            IF ( .NOT.TRANSMISSION .OR. STT ) THEN
C
               IF ( .NOT.TRANSMISSION .AND. .NOT.STT .AND. IK.EQ.IKTOP )
     &              WRITE (6,*) 
     &                   'SIGMA L+R for DOS or SIGMA_R for Transmission'
               IF ( STT .AND. TRANSMISSION .AND. IK.EQ.IKTOP )
     &              WRITE (6,*) 'SIGMA_L for G^<_tr'
               IF ( STT .AND. .NOT.TRANSMISSION .AND. IK.EQ.IKTOP )
     &              WRITE (6,*) 'SIGMA_R for Transmission in STT      ',
     &                          'SIGMA_R for G^<_tr (= 0)'
C
C----------------------------------------------------------------------
C IJ:   sigma_x(i,j;+) dt(j,+) MEZZ(j)^-1 MEZZPM(j) L x(j,+)^(T*)
C     - [dt(i,+) MEZZ(i)^-1 MEZZPM(i) L x(i,+)^(T*)]^(T*) sigma_x(j,i;+)^(T*)
C
C    =  sigma_x(i,j;+) dt(j,+) MEZZ(j)^-1 MEZZPM(j) L x(j,+)^(T*)
C     - x(i,+) L^(T*) MEZZPM(i)^(T*) [MEZZ(i)^-1]^(T*) dt(i,+)^(T*) sigma_x(j,i;+)^(T*)
C----------------------------------------------------------------------
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WIJ,NXM,
     &                    DTSSQX(1,1,JQTB),NXM,C0,WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,MEZZJINV,NXM,
     &                    C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,MEZZJPM2,NXM,
     &                    C0,WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,C0,
     &                    WKM3,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,
     &                    DSSQX(1,1,JQTB),NXM,C0,WKM1,NXM)
C
C-----------------------------------------------------------------------
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,
     &                    LMAT2,NXM,C0,WKM2,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,MEZZIPM2,NXM,
     &                    C0,WKM3,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,MEZZIINV,NXM,
     &                    C0,WKM2,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,
     &                    DTSSQX(1,1,IQTB),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,WJI,NXM,C0,
     &                    WKM2,NXM)
C
               IF ( MINONLY ) WKM1(:,:) = C0
               IF ( PLSONLY ) WKM2(:,:) = C0
C
            END IF
C
            IF ( STT ) THEN
C
               IF ( .NOT.TRANSMISSION ) THEN
C------------ this surely can be done better ...
                  WKMT1(1:NXM,1:NXM) = WKM1(1:NXM,1:NXM)
                  WKMT2(1:NXM,1:NXM) = WKM2(1:NXM,1:NXM)
C------------ If STT and .NOT. TRANSMISSION we are on the right surface,
C------------ whose self energy is not needed for G^<_tr!
CSW                  WKM1(:,:) = C0
CSW                  WKM2(:,:) = C0
               END IF
C
               TAUDELTAK_BAR_TR(IX1:IX2,JX1:JX2) = WKM1(1:NXM,1:NXM)
     &            - WKM2(1:NXM,1:NXM)
               WKM1(1:NXM,1:NXM) = WKMT1(1:NXM,1:NXM)
               WKM2(1:NXM,1:NXM) = WKMT2(1:NXM,1:NXM)
C
            END IF
C
            TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = WKM1(1:NXM,1:NXM)
     &         - WKM2(1:NXM,1:NXM)
C
            IF ( TRANSMISSION ) THEN
C--------------------------------------------------- get pure gamma(-,+)
C
               WKM3(1:NXM,1:NXM) = TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,C0,
     &                    WKM1,NXM)
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DTSSQX(1,1,IQTB),NXM,
     &                    WKM1,NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,C0,
     &                    WKM1,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,
     &                    DTSSQX(1,1,JQTB),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,C0,
     &                    WKM1,NXM)
C
               TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = WKM1(1:NXM,1:NXM)
C
            END IF
C
            IF ( IQ.NE.JQ ) THEN
C
               IF ( TRANSMISSION ) THEN
C
                  IF ( .NOT.STT .AND. IK.EQ.IKTOP ) WRITE (6,*)
     &                  'SIGMA_L for Transmission'
                  IF ( STT .AND. IK.EQ.IKTOP ) WRITE (6,*)
     &                  'SIGMA_L for Transmission in STT'
C
C----------------------------------------------------------------------
C JI:   L x(j,+)^(T*) L MEZZMP(j) MEZZ(j)^-1 dt(j,+) sigma_x(j,i;+) L
C     - [L sigma_x(i,j;+)^(T*) dt(i,+)^(T*) (MEZZ(i)^-1)^(T*) MEZZMP(i)^(T*)
C       L x(i,+) L]
C----------------------------------------------------------------------
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,
     &                       DSSQX(1,1,JQTB),NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,
     &                       C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,MEZZJMP2,
     &                       NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,MEZZJINV,
     &                       NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,
     &                       DTSSQX(1,1,JQTB),NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,WJI,NXM,C0,
     &                       WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,
     &                       C0,WKM1,NXM)
C
C-----------------------------------------------------------------------
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,WIJ,NXM,
     &                       C0,WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,
     &                       DTSSQX(1,1,IQTB),NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,MEZZIINV,
     &                       NXM,C0,WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,MEZZIMP2,
     &                       NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,
     &                       C0,WKM2,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,
     &                       DSSQX(1,1,IQTB),NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,
     &                       C0,WKM2,NXM)
C
                  IF ( PLSONLY ) WKM1(:,:) = C0
                  IF ( MINONLY ) WKM2(:,:) = C0
C
                  IF ( STT ) THEN
                     WKMT1(:,:) = WKM1(:,:)
                     WKMT2(:,:) = WKM2(:,:)
                  END IF
C
               END IF
C
               IF ( .NOT.TRANSMISSION .OR. STT ) THEN
C
                  IF ( .NOT.TRANSMISSION .AND. .NOT.STT .AND. 
     &                 IK.EQ.IKTOP ) WRITE (6,*)
     &                  'SIGMA L+R for DOS or SIGMA_R for Transmission'
                  IF ( STT .AND. TRANSMISSION .AND. IK.EQ.IKTOP )
     &                 WRITE (6,*) 'SIGMA_L for G^<_tr'
                  IF ( STT .AND. .NOT.TRANSMISSION .AND. IK.EQ.IKTOP )
     &                 WRITE (6,*) 
     &                            'SIGMA_R for Transmission in STT     '
     &                            ,' SIGMA_R for G^<_tr (= 0)'
C
C----------------------------------------------------------------------
C JI:   sigma_x(j,i;+) dt(i,+) MEZZ(i)^-1 MEZZPM(i) L x(i,+)^(T*)
C     - [dt(j,+) MEZZ(j)^-1 MEZZPM(j) L x(j,+)^(T*)]^(T*) sigma_x(i,j;+)^(T*)
C
C    =  sigma_x(j,i;+) dt(i,+) MEZZ(i)^-1 MEZZPM(i) L x(i,+)^(T*)
C     - x(j,+) L^(T*) MEZZPM(j)^(T*) [MEZZ(j)^-1]^(T*) dt(j,+)^(T*) sigma_x(i,j;+)^(T*)
C----------------------------------------------------------------------
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WJI,NXM,
     &                       DTSSQX(1,1,IQTB),NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,MEZZIINV,
     &                       NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,MEZZIPM2,
     &                       NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,LMAT2,NXM,
     &                       C0,WKM3,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,
     &                       DSSQX(1,1,IQTB),NXM,C0,WKM1,NXM)
C
C-----------------------------------------------------------------------
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,DSSQX(1,1,JQTB),NXM,
     &                       LMAT2,NXM,C0,WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,MEZZJPM2,
     &                       NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,MEZZJINV,
     &                       NXM,C0,WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,
     &                       DTSSQX(1,1,JQTB),NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM3,NXM,WIJ,NXM,C0,
     &                       WKM2,NXM)
C
                  IF ( MINONLY ) WKM1(:,:) = C0
                  IF ( PLSONLY ) WKM2(:,:) = C0
C
               END IF
C
               IF ( STT ) THEN
C
                  IF ( .NOT.TRANSMISSION ) THEN
C------------ this surely can be done better ...
                     WKMT1(1:NXM,1:NXM) = WKM1(1:NXM,1:NXM)
                     WKMT2(1:NXM,1:NXM) = WKM2(1:NXM,1:NXM)
C------------ If STT and .NOT. TRANSMISSION we are on the right surface,
C------------ whose self energy is not needed for G^<_tr!
CSW                     WKM1(:,:) = C0
CSW                     WKM2(:,:) = C0
                  END IF
C
                  TAUDELTAK_BAR_TR(JX1:JX2,IX1:IX2) = WKM1(1:NXM,1:NXM)
     &               - WKM2(1:NXM,1:NXM)
                  WKM1(1:NXM,1:NXM) = WKMT1(1:NXM,1:NXM)
                  WKM2(1:NXM,1:NXM) = WKMT2(1:NXM,1:NXM)
C
               END IF
C
               TAUDELTAK_BAR(JX1:JX2,IX1:IX2) = WKM1(1:NXM,1:NXM)
     &            - WKM2(1:NXM,1:NXM)
C
               IF ( TRANSMISSION ) THEN
C--------------------------------------------------- get pure gamma(-,+)
C
                  WKM3(1:NXM,1:NXM) = TAUDELTAK_BAR(JX1:JX2,IX1:IX2)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,
     &                       C0,WKM1,NXM)
C
                  CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DTSSQX(1,1,JQTB),
     &                       NXM,WKM1,NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,
     &                       C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,
     &                       DTSSQX(1,1,IQTB),NXM,C0,WKM3,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,
     &                       C0,WKM1,NXM)
C
                  TAUDELTAK_BAR(JX1:JX2,IX1:IX2) = WKM1(1:NXM,1:NXM)
C
               END IF
C
            END IF
C
C===================================================  extract pure gamma
C
C----------------------------- copy block of interest from TAUDELTAK_BAR
            XGAMMAKX(1:NXM,1:NXM) = TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
C
C--------------------------------------------------- get (1/x) = delta t
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,IQTB),WKM1,WKM2)
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,JQTB),WKM1,WKM3)
C
            IF ( TRANSMISSION ) THEN
               GAMMAK(1:NXM,1:NXM) = XGAMMAKX(1:NXM,1:NXM)
               IF ( STT ) THEN
                  GAMMAKT(1:NXM,1:NXM) = GAMMAK(1:NXM,1:NXM)
                  XGAMMAKX(1:NXM,1:NXM)
     &               = TAUDELTAK_BAR_TR(IX1:IX2,JX1:JX2)
C
C--------- multiply from left & right with 1/x^nu & L (1/x^nu')^dagger L
Cleft
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,
     &                       XGAMMAKX(1,1),NXM,C0,WKM1,NXM)
Cright
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,
     &                       C0,WKM2,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,WKM2,NXM,
     &                       C0,GAMMAK,NXM)
C
               END IF
            ELSE
C
C--------- multiply from left & right with 1/x^nu & L (1/x^nu')^dagger L
Cleft
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,XGAMMAKX(1,1),
     &                    NXM,C0,WKM1,NXM)
Cright
               CALL ZGEMM('N','C',NXM,NXM,NXM,C1,LMAT2,NXM,WKM3,NXM,C0,
     &                    WKM2,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,WKM2,NXM,C0,
     &                    GAMMAK,NXM)
C
               IF ( STT ) THEN
                  GAMMAKT(1:NXM,1:NXM) = GAMMAK(1:NXM,1:NXM)
                  GAMMAK(1:NXM,1:NXM) = TAUDELTAK_BAR_TR(IX1:IX2,JX1:JX2
     &                                  )
C---------------------------------- is zero anyway, so no need to purify
               END IF
C
            END IF
C
            IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP ) THEN
               WRITE (43,*) 'IQ =',IQ,'JQ =',JQ
               CALL CMATSTRUCT('GAMMAK',GAMMAK,18,18,2,2,1,1.D-8,43)
               IF ( STT ) THEN
                  WRITE (53,*) 'IQ =',IQ,'JQ =',JQ
                  CALL CMATSTRUCT('GAMMAKT',GAMMAKT,18,18,2,2,1,1.D-8,
     &                            53)
               END IF
            END IF
C
C===================================================  extract pure gamma
C
         END DO
      END DO
C
      IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP )
     &      CALL WR_GORTAUK_ZPM_MO('UGV_BAR',114,TAUDELTAK_BAR,
     &     EXIKDQ_QTB,LMAT2,CCZ)
      IF ( STT .AND. IPRINT.GE.1 .AND. IK.EQ.IKTOP )
     &     CALL WR_GORTAUK_ZPM_MO('UGV_BAR_TR',115,TAUDELTAK_BAR_TR,
     &                            EXIKDQ_QTB,LMAT2,CCZ)
C
      END
C*==negfkloop_gles_mo.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOP_GLES_MO(GLESQX,TAUDELTAK,UGV_BAR,NKKR_TB,
     &                             EXIKDQ_QTB,WMPNP,WMN,BAR_FLAG_Q,
     &                             LBAR_FLAG_Q,RBAR_FLAG_Q,CCZ,LMAT2,IK,
     &                             IKTOP,CHECK_ZPM,WK,WKSUM,
     &                             TRANSMISSION,DSSQX,KVEC,TRANS,STT)
C   ********************************************************************
C   *                                                                  *
C   *  Construct the lesser Green's function from the gammas and taus  *
C   *                                                                  *
C   *   - For DOS G^< is calculated                                    *
C   *                                                                  *
C   *   - For TRANSMISSION and STT the expression                      *
C   *                                                                  *
C   * T(E,k) = -sum_ijkl gam_L,ij(+,k) Gjk(+,k) gam_R,kl(-,k) Gli(-,k) *
C   *                                                                  *
C   *     is evaluated                                                 *
C   *                                                                  *
C   *   - For STT additionally G^<_tr is calculated, with sigma_R = 0  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NLM,NLMMAX,NXM,WKM1,WKM2,WKM3
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,WRBUILDBOT
      USE MOD_LATTICE,ONLY:VOLUC_2D
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP
      USE MOD_SITES,ONLY:IQBOT_TB,IQTOP_TB,NQ,NQTB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL NO_GAMMA
      PARAMETER (NO_GAMMA=.FALSE.)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NEGFKLOOP_GLES_MO')
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM,STT,TRANSMISSION
      INTEGER IK,IKTOP,NKKR_TB
      REAL*8 WK,WKSUM
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),EXIKDQ_QTB(NQTB,NQTB),
     &           GLESQX(NXM,NXM,NQTB),LMAT2(NXM,NXM),
     &           TAUDELTAK(NKKR_TB,NKKR_TB),TRANS(3),
     &           UGV_BAR(NKKR_TB,NKKR_TB),WMN(NXM,NXM),WMPNP(NXM,NXM)
      REAL*8 KVEC(3)
C
C Local variables
C
      COMPLEX*16 CMATTRC
      INTEGER ICALL,IFIL,IIQ,IIQTB,IIX0,IIX1,IIX2,IQ,IQTB,IX0,IX1,IX2,
     &        IXX,JJQ,JJQTB,JJX0,JJX1,JJX2,JQ,JQTB,JX0,JX1,JX2
      COMPLEX*16 TRANSK,TRANSKDN,TRANSKUP,TRANSM(:,:),TRANSMDN(:,:),
     &           TRANSMUP(:,:)
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE TRANSM,TRANSMDN,TRANSMUP
C
      ALLOCATE (TRANSM(NXM,NXM),TRANSMDN(NLM,NLM),TRANSMUP(NLM,NLM))
C
C*** End of declarations rewritten by SPAG
C
      TRANSK = C0
      TRANSKUP = C0
      TRANSKDN = C0
      ICALL = ICALL + 1
C
      IF ( ICALL.EQ.1 ) TRANS(1:3) = C0
C
C======================================================================
C
C------------------------------------------------------------------- mu'
      LOOP_JJQ:DO JJQ = IQBOT_TB,IQTOP_TB
         JJQTB = JJQ - IQBOT_TB + 1
         JJX0 = (JJQTB-1)*NXM
         JJX1 = JJX0 + 1
         JJX2 = JJX0 + NXM
C
C-------------------------------------------------------------------- mu
         LOOP_IIQ:DO IIQ = IQBOT_TB,IQTOP_TB
            IF ( TRANSMISSION ) THEN
               IF ( .NOT.(LBAR_FLAG_Q(IIQ) .AND. LBAR_FLAG_Q(JJQ)) )
     &              CYCLE LOOP_IIQ
            ELSE IF ( IIQ.NE.JJQ ) THEN
               CYCLE LOOP_IIQ
            END IF
            IIQTB = IIQ - IQBOT_TB + 1
            IIX0 = (IIQTB-1)*NXM
            IIX1 = IIX0 + 1
            IIX2 = IIX0 + NXM
C
C----------------------------------------------------------------------
C                        BARRIERS
C----------------------------------------------------------------------
C
C------------------------------------------------------------------- nu'
            LOOP_JQ:DO JQ = IQBOT_TB,IQTOP_TB
               IF ( .NOT.BAR_FLAG_Q(JQ) ) CYCLE LOOP_JQ
               JQTB = JQ - IQBOT_TB + 1
               JX0 = (JQTB-1)*NXM
               JX1 = JX0 + 1
               JX2 = JX0 + NXM
C
               WMPNP(1:NXM,1:NXM) = TAUDELTAK(JJX1:JJX2,JX1:JX2)
     &                              *EXIKDQ_QTB(JJQTB,JQTB)
C
               IF ( TRANSMISSION .OR. STT ) THEN
CSW------------ this step is done after the k-loop for DOS calculations!
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,JJQTB),
     &                       NXM,WMPNP,NXM,C0,WKM1,NXM)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM1,NXM,
     &                       C0,WMPNP,NXM)
C
               END IF
C
C-------------------------------------------------------------------- nu
               LOOP_IQ:DO IQ = IQBOT_TB,IQTOP_TB
                  IF ( TRANSMISSION ) THEN
                     IF ( .NOT.((LBAR_FLAG_Q(IIQ) .AND. LBAR_FLAG_Q(JJQ)
     &                    ) .AND. 
     &                    (RBAR_FLAG_Q(IQ) .AND. RBAR_FLAG_Q(JQ))) )
     &                    CYCLE LOOP_IQ
                  ELSE IF ( .NOT.((LBAR_FLAG_Q(IQ) .AND. LBAR_FLAG_Q(JQ)
     &                      ) .OR. 
     &                      (RBAR_FLAG_Q(IQ) .AND. RBAR_FLAG_Q(JQ))) )
     &                      THEN
                     CYCLE LOOP_IQ
                  END IF
                  IQTB = IQ - IQBOT_TB + 1
                  IX0 = (IQTB-1)*NXM
                  IX1 = IX0 + 1
                  IX2 = IX0 + NXM
C
                  WMN(1:NXM,1:NXM) = TAUDELTAK(IIX1:IIX2,IX1:IX2)
     &                               *EXIKDQ_QTB(IIQTB,IQTB)
C
                  IF ( TRANSMISSION .OR. STT ) THEN
CSW------------ this step is done after the k-loop for DOS calculations!
C
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IIQTB),
     &                          NXM,WMN,NXM,C0,WKM1,NXM)
C
                     WMN(1:NXM,1:NXM) = WKM1(1:NXM,1:NXM)
C
                  END IF
C
                  IF ( NO_GAMMA ) THEN
                     UGV_BAR(IX1:IX2,JX1:JX2) = C0
                     DO IXX = IX1,IX2
                        UGV_BAR(IXX,IXX) = C1
                     END DO
                  END IF
C
                  WKM1(1:NXM,1:NXM) = UGV_BAR(IX1:IX2,JX1:JX2)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WMN,NXM,WKM1,NXM,C0,
     &                       WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,WMPNP,NXM,
     &                       C0,WKM1,NXM)
C
                  IF ( TRANSMISSION ) THEN
C
                     WKM2(1:NXM,1:NXM) = UGV_BAR(JJX1:JJX2,IIX1:IIX2)
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,WKM1,
     &                          NXM,C0,WKM3,NXM)
C
                     IF ( IPRINT.GE.1 .AND. IK.EQ.IKTOP ) THEN
                        WRITE (44,*) 'IQ =',IQ,'JQ =',JQ,'IIQ =',IIQ,
     &                               'JJQ =',JJQ
                        CALL CMATSTRUCT('TIJIIJJ',WKM3,18,18,2,2,1,
     &                                  1.D-8,44)
                     END IF
C
                     TRANSM(1:NXM,1:NXM) = WKM3(1:NXM,1:NXM)
                     TRANSMUP(1:NLM,1:NLM) = WKM3(1:NLM,1:NLM)
                     TRANSMDN(1:NLM,1:NLM) = WKM3(NLM+1:NXM,NLM+1:NXM)
C
CSW--------------- dividing by area of 2D UC to get T(E) per unit area!?
C
                     TRANSK = TRANSK - CMATTRC(NKM,NKMMAX,TRANSM)
                     TRANS(1) = TRANS(1) - CMATTRC(NKM,NKMMAX,TRANSM)
     &                          *WK/WKSUM/VOLUC_2D
C
                     TRANSKUP = TRANSKUP - CMATTRC(NLM,NLMMAX,TRANSMUP)
                     TRANS(2) = TRANS(2) - CMATTRC(NLM,NLMMAX,TRANSMUP)
     &                          *WK/WKSUM/VOLUC_2D
C
                     TRANSKDN = TRANSKDN - CMATTRC(NLM,NLMMAX,TRANSMDN)
                     TRANS(3) = TRANS(3) - CMATTRC(NLM,NLMMAX,TRANSMDN)
     &                          *WK/WKSUM/VOLUC_2D
C
                  ELSE
C
                     GLESQX(1:NXM,1:NXM,IIQTB)
     &                  = GLESQX(1:NXM,1:NXM,IIQTB) + WKM1(1:NXM,1:NXM)
     &                  *WK/WKSUM
CSW  &               + CI*WKM1(1:NXM,1:NXM)*WK/WKSUM
CSW--------------------------------- some i''s & signs cancel each other
C
                  END IF
C
               END DO LOOP_IQ
            END DO LOOP_JQ
C
C======================================================================
C
         END DO LOOP_IIQ
      END DO LOOP_JJQ
C======================================================================
C
      IF ( TRANSMISSION ) THEN
C
         IFIL = 234 + MPI_ID
         WRITE (IFIL,FMT='(I10,2F12.6,3F12.6)') IK,KVEC(1:2),
     &          DREAL(TRANSK),DREAL(TRANSKUP),DREAL(TRANSKDN)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT .AND. MPI_ID.EQ.0 ) WRITE (IFILBUILDBOT,99001)
     &        ROUTINE(1:LEN_TRIM(ROUTINE)),IK,KVEC(1:2),DREAL(TRANSK),
     &        DREAL(TRANSKUP),DREAL(TRANSKDN)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C
         IF ( IK.EQ.IKTOP ) THEN
            IF ( .NOT.MPI_KLOOP ) THEN
               WRITE (6,FMT='(/,18X,A,10X,A)') 'Re','Im'
               WRITE (6,FMT='(A,2F12.6)') 'T(E) total: ',TRANS(1)
               WRITE (6,FMT='(A,2F12.6)') 'T(E) spin1: ',TRANS(2)
               WRITE (6,FMT='(A,2F12.6)') 'T(E) spin2: ',TRANS(3)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99002)
     &                                  ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                  DREAL(TRANS(1:3)),
     &                                  DIMAG(TRANS(1:3))
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            END IF
         END IF
      END IF
C
      IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &      CALL WR_GORTAUK_ZPM_MO('TAUDELTAK',115,TAUDELTAK,EXIKDQ_QTB,
     &     LMAT2,CCZ)
C
99001 FORMAT ('# BUILDBOT: ',A,':  Re(T(k_II)) for IK =',I2,' (',2F10.6,
     &        '): tot, up & dn',/,(1PE22.14))
99002 FORMAT ('# BUILDBOT: ',A,':  T(E) total, up & dn (3xRe, 3xIm)',/,
     &        (1PE22.14))
      END
C*==wr_gortauk_zpm_mo.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE WR_GORTAUK_ZPM_MO(STR,IFIL,GORTAUK,EXIKDQ_QTB,LMAT2,
     &                             CCZ)
C   ********************************************************************
C   *                                                                  *
C   *   Write out GREFKE, GLLKE, or TAUDELTAK block-wise for z & z*    *
C   *   to check wether                                                *
C   *                                                                  *
C   *     X_LL'^nn' (z) = L_LL' * (X_LL'^n'n (z*))^dagger * L_LL'      *
C
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NXM
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SITES,ONLY:NQTB,IQBOT_TB
      USE MOD_TB,ONLY:NKKR_TB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CCZ
      INTEGER IFIL
      CHARACTER*(*) STR
      COMPLEX*16 EXIKDQ_QTB(NQTB,NQTB),GORTAUK(NKKR_TB,NKKR_TB),
     &           LMAT2(NXM,NXM)
C
C Local variables
C
      INTEGER IQ,IQTB,IQTOP_TB,IX0,IX1,IX2,JQ,JQTB,JX0,JX1,JX2
      COMPLEX*16 PHFCTR(NQTB,NQTB),WORK1(NXM,NXM),WORK2(NXM,NXM)
C
C*** End of declarations rewritten by SPAG
C
      IQTOP_TB = IQBOT_TB + NQTB - 1
C
      PHFCTR(:,:) = EXIKDQ_QTB(:,:)
      IF ( STR(1:5).EQ.'UGV_B' .OR. STR(1:5).EQ.'TAUDE' ) PHFCTR(:,:)
     &     = C1
C
      IF ( CCZ ) THEN
C======================================================================
C
C-------------------------------------------------------------------- mu
         DO IQ = IQBOT_TB,IQTOP_TB
            IQTB = IQ - IQBOT_TB + 1
            IX0 = (IQTB-1)*NXM
            IX1 = IX0 + 1
            IX2 = IX0 + NXM
C
C------------------------------------------------------------------- mu'
            DO JQ = IQBOT_TB,IQTOP_TB
               JQTB = JQ - IQBOT_TB + 1
               JX0 = (JQTB-1)*NXM
               JX1 = JX0 + 1
               JX2 = JX0 + NXM
C
               IF ( IQ.EQ.JQ ) THEN
                  WORK1(1:NXM,1:NXM) = GORTAUK(JX1:JX2,IX1:IX2)
               ELSE
                  WORK1(1:NXM,1:NXM) = GORTAUK(JX1:JX2,IX1:IX2)
     &                                 *PHFCTR(JQTB,IQTB)
               END IF
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WORK1,NXM,C0,
     &                    WORK2,NXM)
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WORK2,NXM,LMAT2,NXM,C0,
     &                    WORK1,NXM)
C
               WRITE (IFIL,*) STR,' for jq =',JQ,'iq =',IQ
               CALL CMATSTRUCT('GORTAUK_MPM (-)',WORK1(:,:),18,18,2,2,1,
     &                         1.D-8,IFIL)
               IF ( IQ.NE.JQ ) WRITE (IFIL,*) 'PHASE FACTOR =',
     &                                PHFCTR(JQTB,IQTB)
               WRITE (IFIL,*) ''
C
            END DO
         END DO
C======================================================================
C
      ELSE
C
C======================================================================
C
C------------------------------------------------------------------- mu'
         DO JQ = IQBOT_TB,IQTOP_TB
            JQTB = JQ - IQBOT_TB + 1
            JX0 = (JQTB-1)*NXM
            JX1 = JX0 + 1
            JX2 = JX0 + NXM
C
C-------------------------------------------------------------------- mu
            DO IQ = IQBOT_TB,IQTOP_TB
               IQTB = IQ - IQBOT_TB + 1
               IX0 = (IQTB-1)*NXM
               IX1 = IX0 + 1
               IX2 = IX0 + NXM
C
               IF ( IQ.EQ.JQ ) THEN
                  WORK1(1:NXM,1:NXM) = GORTAUK(JX1:JX2,IX1:IX2)
               ELSE
                  WORK1(1:NXM,1:NXM) = GORTAUK(JX1:JX2,IX1:IX2)
     &                                 *PHFCTR(JQTB,IQTB)
               END IF
C
               WRITE (IFIL,*) STR,' for jq =',JQ,'iq =',IQ
               CALL CMATSTRUCT('GORTAUK_MPM (+)',WORK1(:,:),18,18,2,2,1,
     &                         1.D-8,IFIL)
               IF ( IQ.NE.JQ ) WRITE (IFIL,*) 'PHASE FACTOR =',
     &                                PHFCTR(JQTB,IQTB)
               WRITE (IFIL,*) ''
C
            END DO
         END DO
C======================================================================
      END IF
C
      END
C*==tbdecimate_mo.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE TBDECIMATE_MO(GLLKE,DSSQX_L,DSSQX_R,VACFLAG,FACTL,NQ_L,
     &                         NQ_R,NXM,SIGMA_L,SIGMA_R)
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
     &           FACTL(NXM,NXM),GLLKE(NKKR_TB,NKKR_TB),
     &           SIGMA_L(NKMSLAY,NKMSLAY),SIGMA_R(NKMSLAY,NKMSLAY)
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
C
CSW copy self-energy of left surface
         SIGMA_L(:,:) = X1(:,:)
CSW copy self-energy of left surface
C
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
CSW copy self-energy of right surface
            SIGMA_R(:,:) = XN(:,:)
CSW copy self-energy of right surface
C
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
CSW copy self-energy of right surface
            SIGMA_R(:,:) = XN(:,:)
CSW copy self-energy of right surface
C
C
         END IF
C=======================================================================
C     adds to the matrix GLLKE the elements that couple the
C     interfaces to the two half-spaces.
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
