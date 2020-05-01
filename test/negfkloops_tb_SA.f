C*==negfkloops_tb_sa.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOPS_TB_SA(GLESQX,TAUQX,TSSQX,DSSQX,DSSQX_L,
     &                            DSSQX_R,TAUQX_BAR,TSSQX_BAR,DSSQX_BAR,
     &                            USSQX_LBAR,VSSQX_LBAR,GSSQX_LBAR,
     &                            USSQX_RBAR,VSSQX_RBAR,GSSQX_RBAR,
     &                            GREF_I1,ICLU_REF_IQTB,WA,WB,LMAT3,
     &                            GLESQ,TAUQ,TAUQ_BAR,IQBOT_LBAR,
     &                            IQTOP_LBAR,IQBOT_RBAR,IQTOP_RBAR,
     &                            LBAR_FLAG_Q,RBAR_FLAG_Q,BAR_FLAG_Q,
     &                            CCZ,CHECK_ZPM)
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
      USE MOD_KSPACE,ONLY:NKTAB,WKTAB,KTAB
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,KTABTAUIJ
      USE MOD_SYMMETRY,ONLY:DROT,NSYMACCEPTED,NSYM,SYMACCEPTED,
     &    SYMUNITARY
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:C0,C1,CI2PI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NLMMAX,NKM,NLM,NXM,NKMMAX,NKMQ,WKM1,WKM2
      USE MOD_SITES,ONLY:NQ_L,NQ_R,NQTB,IQBOT_TB,NQMAX,QBAS,NQ
      USE MOD_TBCLU,ONLY:NKKRNR_RSMAX,NCLU_REF,IQTBORGQTBP
      USE MOD_TB,ONLY:NKKRNR_TB,NKKR_TB,VACFLAG,IDECI,INVMOD,FACTL,
     &    ICHECK,NPLAY
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL USE_FULL_SYMMETRY
      PARAMETER (USE_FULL_SYMMETRY=.TRUE.)
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM
      INTEGER IQBOT_LBAR,IQBOT_RBAR,IQTOP_LBAR,IQTOP_RBAR
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_BAR(NXM,NXM,NQTB),
     &           DSSQX_L(NXM,NXM,NQ_L),DSSQX_R(NXM,NXM,NQ_R),
     &           GLESQ(NKMMAX,NKMMAX,NQ),GLESQX(NXM,NXM,NQTB),
     &           GREF_I1(NKKRNR_RSMAX,NLMMAX,NCLU_REF),
     &           GSSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR),
     &           GSSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR),
     &           LMAT3(NXM,NXM),TAUQ(NKMMAX,NKMMAX,NQ),
     &           TAUQX(NXM,NXM,NQTB),TAUQX_BAR(NXM,NXM,NQTB),
     &           TAUQ_BAR(NKMMAX,NKMMAX,NQMAX),TSSQX(NXM,NXM,NQTB),
     &           TSSQX_BAR(NXM,NXM,NQTB),
     &           USSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR,2),
     &           USSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR,2),
     &           VSSQX_LBAR(NXM,NXM,IQBOT_LBAR:IQTOP_LBAR,2),
     &           VSSQX_RBAR(NXM,NXM,IQBOT_RBAR:IQTOP_RBAR,2),WA(NXM,NXM)
     &           ,WB(NXM,NXM)
      INTEGER ICLU_REF_IQTB(NQTB)
C
C Local variables
C
      COMPLEX*16 AUXQ(:,:,:),EXIKDQ_QTB(:,:),GLLKE(:,:),GLLKE_BAR(:,:),
     &           GREFKE(:,:),LMAT2(NXM,NXM),WIJ(:,:),WJI(:,:)
      REAL*8 DDOT
      REAL*8 DQ_QTB(:,:,:),IM_KVEC(3),KVEC(3),WK,WKSUM
      INTEGER I,IA_ERR,IK,IKTOP,IPROCK(:),IQ,IQTB,IU,IU0,IV,IV0,IX,IX0,
     &        J,JQ,JQTB,JU,JU0,JV,JV0,JX,M,N1,N2
C
C*** End of declarations rewritten by SPAG
C
Csw      PARAMETER (USE_FULL_SYMMETRY=.FALSE.)
C
      DATA IM_KVEC/0D0,0D0,0D0/
C
      ALLOCATABLE DQ_QTB,EXIKDQ_QTB
      ALLOCATABLE AUXQ,WIJ,WJI
      ALLOCATABLE GLLKE,GREFKE,GLLKE_BAR,IPROCK
C
      ALLOCATE (AUXQ(NXM,NXM,NQTB),WIJ(NXM,NXM),WJI(NXM,NXM))
      ALLOCATE (DQ_QTB(3,NQTB,NQTB),EXIKDQ_QTB(NQTB,NQTB))
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
      ALLOCATE (GLLKE_BAR(NKKR_TB,NKKR_TB),STAT=IA_ERR)
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
      TAUQX_BAR(:,:,:) = C0
      GLESQX(:,:,:) = C0
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
Csw
C     IKTOP=1
C     WKSUM=1D0
      WRITE (6,*) 'IKTOP =',IKTOP
                                 !actual number of points, nowhere given?
      WRITE (6,*) 'WKSUM =',WKSUM
                                 !already written out somewhere
Csw
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
Csw
            IF ( IKTOP.EQ.1 ) KVEC(1:3) = (/0.5D0,0.5D0,0D0/)
C
            IF ( IPRINT.GE.1 ) THEN
               WRITE (345,*) 'IK =',IK
               WRITE (345,*) 'KVEC =',KVEC
               WRITE (345,*) 'WK =',WK
            END IF
Csw
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
Csw
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM('GREFKE',111,GLLKE,EXIKDQ_QTB,
     &           LMAT2,CCZ)
Csw
C-----------------------------------------------------------------------
C              call decimation routine if requested
C-----------------------------------------------------------------------
C
            IF ( IDECI.EQ.1 ) CALL TBDECIMATE(GLLKE,DSSQX_L,DSSQX_R,
     &           VACFLAG,FACTL,NQ_L,NQ_R,NXM)
Csw
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &            CALL WR_GORTAUK_ZPM('GREFKE_DEC',1111,GLLKE,
     &           EXIKDQ_QTB,LMAT2,CCZ)
Csw
C
            GLLKE_BAR(:,:) = GLLKE(:,:)
C
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
               IF ( LBAR_FLAG_Q(IQ) ) THEN
C
                  DO J = 1,NXM
                     JX = JX + 1
                     IX = IX0
                     DO I = 1,NXM
                        IX = IX + 1
                        GLLKE(IX,JX) = GLLKE(IX,JX) - DSSQX(I,J,IQTB)
                        GLLKE_BAR(IX,JX) = GLLKE_BAR(IX,JX)
     &                     - DSSQX_BAR(I,J,IQTB)
                     END DO
                  END DO
C
               ELSE IF ( RBAR_FLAG_Q(IQ) ) THEN
C
                  DO J = 1,NXM
                     JX = JX + 1
                     IX = IX0
                     DO I = 1,NXM
                        IX = IX + 1
                        GLLKE(IX,JX) = GLLKE(IX,JX) - DSSQX(I,J,IQTB)
                        GLLKE_BAR(IX,JX) = GLLKE_BAR(IX,JX)
     &                     - DSSQX_BAR(I,J,IQTB)
                     END DO
                  END DO
C
               ELSE
C
                  DO J = 1,NXM
                     JX = JX + 1
                     IX = IX0
                     DO I = 1,NXM
                        IX = IX + 1
                        GLLKE(IX,JX) = GLLKE(IX,JX) - DSSQX(I,J,IQTB)
                        GLLKE_BAR(IX,JX) = GLLKE_BAR(IX,JX)
     &                     - DSSQX_BAR(I,J,IQTB)
                     END DO
                  END DO
C
               END IF
C
Csw
               IF ( IK.EQ.IKTOP ) THEN
                  WRITE (654,*) 'IQTB =',IQTB
                  CALL CMATSTRUCT('DSSQX',DSSQX(:,:,IQTB),18,18,2,2,1,
     &                            1.D-8,654)
               END IF
Csw
            END DO
Csw
            IF ( CHECK_ZPM .AND. IK.EQ.IKTOP ) THEN
               CALL WR_GORTAUK_ZPM('GLLKE',112,GLLKE,EXIKDQ_QTB,LMAT2,
     &                             CCZ)
C
               CALL WR_GORTAUK_ZPM('GLLKE_BAR',113,GLLKE_BAR,EXIKDQ_QTB,
     &                             LMAT2,CCZ)
            END IF
Csw
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
            CALL TBINVERSION(GLLKE_BAR,INVMOD,ICHECK,NXM,DSSQX_BAR,WK)
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
                     TAUQX_BAR(I,J,IQTB) = TAUQX_BAR(I,J,IQTB)
     &                  - WK*GLLKE_BAR(IX,JX)/WKSUM
                  END DO
               END DO
Csw
CSW not reasonable(?) in MPI case!! at least it goes wrong...
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
                     WKM1(:,:) = TAUQX_BAR(:,:,IQTB)
C
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM1,
     &                          NXM,C0,WKM2,NXM)
                     CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,LMAT2,
     &                          NXM,C0,WKM1,NXM)
C
                     WRITE (92,*) 'IQTB =',IQTB
                     CALL CMATSTRUCT('TAUQX_BAR (-)',WKM1(:,:),18,18,2,
     &                               2,1,1.D-8,92)
C
                  ELSE
C
                     WRITE (91,*) 'IQTB =',IQTB
                     CALL CMATSTRUCT('TAUQX (+)',TAUQX(:,:,IQTB),18,18,
     &                               2,2,1,1.D-8,91)
C
                     WRITE (92,*) 'IQTB =',IQTB
                     CALL CMATSTRUCT('TAUQX_BAR (+)',TAUQX_BAR(:,:,IQTB)
     &                               ,18,18,2,2,1,1.D-8,92)
                  END IF
               END IF
Csw
            END DO
C
C-----------------------------------------------------------------------
C    convert TAU_Delta  for the barriers to renormalized  G_bar
C-----------------------------------------------------------------------
C
C                                                           left barrier
C
            CALL NEGFKLOOP_GBAR(GLLKE_BAR,USSQX_LBAR,VSSQX_LBAR,
     &                          GSSQX_LBAR,DSSQX,DSSQX_BAR,NKKR_TB,
     &                          EXIKDQ_QTB,WIJ,WJI,IQBOT_LBAR,
     &                          IQTOP_LBAR,CCZ,LMAT2,IK,IKTOP,CHECK_ZPM)
C
C                                                          right barrier
C
            CALL NEGFKLOOP_GBAR(GLLKE_BAR,USSQX_RBAR,VSSQX_RBAR,
     &                          GSSQX_RBAR,DSSQX,DSSQX_BAR,NKKR_TB,
     &                          EXIKDQ_QTB,WIJ,WJI,IQBOT_RBAR,
     &                          IQTOP_RBAR,CCZ,LMAT2,IK,IKTOP,CHECK_ZPM)
C
C-----------------------------------------------------------------------
C                    calculate lesser Green's function   G<
C-----------------------------------------------------------------------
C
            CALL NEGFKLOOP_GLES(GLESQX,GLLKE,GLLKE_BAR,NKKR_TB,
     &                          EXIKDQ_QTB,WIJ,WJI,BAR_FLAG_Q,
     &                          LBAR_FLAG_Q,RBAR_FLAG_Q,CCZ,LMAT2,IK,
     &                          IKTOP,CHECK_ZPM,WK,WKSUM)
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
         CALL DRV_MPI_REDUCE_C(TAUQX_BAR(1,1,1),TAUQ(1,1,1),M)
         CALL DRV_MPI_REDUCE_C(GLESQX(1,1,1),TAUQ(1,1,1),M)
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
C ----------------------------------------------------------------------
C                          barrier system
C ----------------------------------------------------------------------
C
C------------- G = 1/(delta t) * TAU_Delta *1/(delta t) - 1/(delta t)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TAUQX_BAR(1,1,IQTB),NXM,
     &                 DSSQX_BAR(1,1,IQTB),NXM,C0,WB,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX_BAR(1,1,IQTB),NXM,
     &                 WB,NXM,C0,WA,NXM)
C
            WB(1:NXM,1:NXM) = WA(1:NXM,1:NXM)
     &                        - DSSQX_BAR(1:NXM,1:NXM,IQTB)
C
C--------------------------------------------------- TAU = t * G * t + t
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WB,NXM,TSSQX_BAR(1,1,IQTB)
     &                 ,NXM,C0,WA,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,TSSQX_BAR(1,1,IQTB),NXM,
     &                 WA,NXM,C0,WB,NXM)
C
            DO J = 1,NXM
               DO I = 1,NXM
                  TAUQX_BAR(I,J,IQTB) = WB(I,J) + TSSQX_BAR(I,J,IQTB)
               END DO
            END DO
Csw
            WRITE (191,*) 'IQTB =',IQTB
            CALL CMATSTRUCT('TAUQX (+)',TAUQX(:,:,IQTB),18,18,2,2,1,
     &                      1.D-8,191)
C
            WRITE (192,*) 'IQTB =',IQTB
            CALL CMATSTRUCT('TAUQX_BAR (+)',TAUQX_BAR(:,:,IQTB),18,18,2,
     &                      2,1,1.D-8,192)
Csw
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
C                          barrier system
C
               WA(1:NKM,1:NKM) = TAUQX_BAR(1:NKM,1:NKM,IQTB)
C
               CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
               TAUQX_BAR(1:NKM,1:NKM,IQTB) = WB(1:NKM,1:NKM)
C ----------------------------------------------------------------------
Csw
               WRITE (291,*) 'IQTB =',IQTB
               CALL CMATSTRUCT('TAUQX (+)',TAUQX(:,:,IQTB),18,18,3,3,1,
     &                         1.D-8,291)
C
               WRITE (292,*) 'IQTB =',IQTB
               CALL CMATSTRUCT('TAUQX_BAR (+)',TAUQX_BAR(:,:,IQTB),18,
     &                         18,3,3,1,1.D-8,292)
Csw
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
Csw
            WRITE (391,*) 'IQTB =',IQTB
            CALL CMATSTRUCT('TAUQ (+)',TAUQ(:,:,IQ),18,18,3,3,1,1.D-8,
     &                      391)
Csw
         END DO
C
C ----------------------------------------------------------------------
C                          barrier system
C ----------------------------------------------------------------------
C
         CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,1D0,TAUQX_BAR(1,1,1),
     &                  AUXQ(1,1,1),WA,NQTB,NKMQ(IQBOT_TB),DROT,
     &                  IQTBORGQTBP,SYMUNITARY,SYMACCEPTED,NSYM,
     &                  NSYMACCEPTED,NQTB,NKMMAX)
C
C-----------------------------------------------------------------------
C                   set up the matrix   TAUQ_BAR
C-----------------------------------------------------------------------
C
         DO IQTB = 1,NQTB
C
            IQ = IQBOT_TB - 1 + IQTB
            DO J = 1,NXM
               DO I = 1,NXM
                  TAUQ_BAR(I,J,IQ) = AUXQ(I,J,IQTB)
               END DO
            END DO
C
Csw
            WRITE (392,*) 'IQTB =',IQTB
            CALL CMATSTRUCT('TAUQ_BAR (+)',TAUQ_BAR(:,:,IQ),18,18,3,3,1,
     &                      1.D-8,392)
Csw
         END DO
C
C-----------------------------------------------------------------------
C                   set up the matrix   GLESQ
C-----------------------------------------------------------------------
C
         DO IQ = 1,NQ
C
            IQTB = IQ - IQBOT_TB + 1
C
            IF ( IQ.LE.IQTOP_LBAR .OR. IQ.GE.IQBOT_RBAR ) THEN
C
               GLESQ(:,:,IQ) = C0
C
            ELSE
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,
     &                    LMAT2,NXM,C0,WA,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,GLESQX(1,1,IQTB),NXM,
     &                    WA,NXM,C0,WB,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,DSSQX(1,1,IQTB),NXM,WB,
     &                    NXM,C0,WA,NXM)
C
               CALL CHANGEREP(NKM,NKMMAX,WA,'RLM>REL',WB)
C
               GLESQ(1:NKM,1:NKM,IQ) = WB(1:NKM,1:NKM)
C
C ----------------------------------------------------------------------
Csw NO symmetrization of glesq!!?
C
            END IF
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
      END
C*==negfkloop_gbar.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOP_GBAR(TAUDELTAK_BAR,USSQX_XBAR,VSSQX_XBAR,
     &                          GSSQX_XBAR,DSSQX,DSSQX_BAR,NKKR_TB,
     &                          EXIKDQ_QTB,WIJ,WJI,IQBOT_XBAR,
     &                          IQTOP_XBAR,CCZ,LMAT2,IK,IKTOP,CHECK_ZPM)
C   ********************************************************************
C   *                                                                  *
C   *  transfer  TAU_Delta  to  auxilary matrix by applying the        *
C   *  transformation matrices U and V and sum the 2 contributions     *
C   *                                                                  *
C   *    u(i;+) tau(i,j;+) v(j,+)  -  u(i;-) tau(j,i;+)^(T*) v(j,-)    *
C   *                                                                  *
C   *  - TAUDELTAK_BAR in:     TAU_Delta(k,E)                          *
C   *                 out:     renormalized G_bar(k,E)                 *
C   *  - called for X= L (left and R (right) barrier separately        *
C   *  - the phase factor exp(ik Delta R) is added                     *
C   *  - the single site contribution GAM is added                     *
C   *                                                                  *
C   *  index IZ = 1 --> E(+)  = 2 --> E(-)                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_ANGMOM,ONLY:NKMMAX,NXM,WKM1,WKM2,WKM3,IPIVKM
      USE MOD_SITES,ONLY:NQTB,IQBOT_TB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL SSONLY,BSONLY,PLSONLY,NOBSSD,UVUNITY,ONLYBSSD,BSCORRSS
      PARAMETER (SSONLY=.FALSE.,BSONLY=.FALSE.,PLSONLY=.FALSE.,
     &           NOBSSD=.FALSE.,UVUNITY=.FALSE.,ONLYBSSD=.FALSE.,
     &           BSCORRSS=.FALSE.)
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM
      INTEGER IK,IKTOP,IQBOT_XBAR,IQTOP_XBAR,NKKR_TB
      COMPLEX*16 DSSQX(NXM,NXM,NQTB),DSSQX_BAR(NXM,NXM,NQTB),
     &           EXIKDQ_QTB(NQTB,NQTB),
     &           GSSQX_XBAR(NXM,NXM,IQBOT_XBAR:IQTOP_XBAR),
     &           LMAT2(NXM,NXM),TAUDELTAK_BAR(NKKR_TB,NKKR_TB),
     &           USSQX_XBAR(NXM,NXM,IQBOT_XBAR:IQTOP_XBAR,2),
     &           VSSQX_XBAR(NXM,NXM,IQBOT_XBAR:IQTOP_XBAR,2),
     &           WIJ(NXM,NXM),WJI(NXM,NXM)
C
C Local variables
C
      COMPLEX*16 GAMMAK(:,:),XGAMMAKX(:,:),XINV(:,:)
      INTEGER I,IQ,IQTB,IX0,IX1,IX2,JQ,JQTB,JX0,JX1,JX2
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE GAMMAK,XGAMMAKX,XINV
C
      ALLOCATE (GAMMAK(NXM,NXM),XGAMMAKX(NXM,NXM),XINV(NXM,NXM))
C
C*** End of declarations rewritten by SPAG
C
      IF ( SSONLY .AND. IK.EQ.1 ) WRITE (6,*)
     &      'SSSSSSSSSSSSSSSSS    SSONLY    SSSSSSSSSSSSSSSSS'
      IF ( BSONLY .AND. IK.EQ.1 ) WRITE (6,*)
     &      'BBBBBBBBBBBBBBBBB    BSONLY    BBBBBBBBBBBBBBBBB'
      IF ( SSONLY .AND. BSONLY ) STOP '<NEGFKLOOP_GBAR>: SS or BS only?'
      IF ( PLSONLY .AND. IK.EQ.1 ) WRITE (6,*)
     &      '+++++++++++++++++    PLSONLY    +++++++++++++++++'
      IF ( .NOT.BSONLY .AND. PLSONLY )
     &      STOP '<NEGFKLOOP_GBAR>: PLSONLY applies to BS-term!'
      IF ( NOBSSD .AND. IK.EQ.1 ) WRITE (6,*)
     &      'OOOOOOOOOOOOOOOOO    NOBSSD     OOOOOOOOOOOOOOOOO'
      IF ( .NOT.BSONLY .AND. NOBSSD )
     &      STOP '<NEGFKLOOP_GBAR>: NOBSSD applies to BS-term!'
      IF ( ONLYBSSD .AND. IK.EQ.1 ) WRITE (6,*)
     &      'YYYYYYYYYYYYYYYY    ONLYBSSD     YYYYYYYYYYYYYYYY'
      IF ( .NOT.BSONLY .AND. ONLYBSSD )
     &      STOP '<NEGFKLOOP_GBAR>: ONLYBSSD applies to BS-term!'
      IF ( BSCORRSS .AND. IK.EQ.1 ) WRITE (6,*)
     &      'CCCCCCCCCCCCCCCC    BSCORRSS     CCCCCCCCCCCCCCCC'
      IF ( .NOT.BSONLY .AND. BSCORRSS )
     &      STOP '<NEGFKLOOP_GBAR>: BSCORRSS applies to BS-term!'
      IF ( BSCORRSS .AND. .NOT.ONLYBSSD ) STOP 
     &     '<NEGFKLOOP_GBAR>: BSCORRSS applies to site-diag. g_BS!'
C
C----------------------------------------------------------------------
C
      DO JQ = IQBOT_XBAR,IQTOP_XBAR
         JQTB = JQ - IQBOT_TB + 1
         JX0 = (JQTB-1)*NXM
         JX1 = JX0 + 1
         JX2 = JX0 + NXM
C
C------------------------------------------------------------- get (1/x)
         CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX_BAR(1,1,JQTB),WKM1,XINV)
C
         DO IQ = IQBOT_XBAR,JQ
            IQTB = IQ - IQBOT_TB + 1
            IX0 = (IQTB-1)*NXM
            IX1 = IX0 + 1
            IX2 = IX0 + NXM
C
            WIJ(1:NXM,1:NXM) = TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
     &                         *EXIKDQ_QTB(IQTB,JQTB)
            WJI(1:NXM,1:NXM) = TAUDELTAK_BAR(JX1:JX2,IX1:IX2)
     &                         *EXIKDQ_QTB(JQTB,IQTB)
C
C----------------------------- subtract 1/x(+) from site diagonal blocks
C
            IF ( IQ.EQ.JQ ) THEN
               IF ( BSCORRSS ) THEN
                  WIJ(:,:) = C0
                  WJI(:,:) = C0
               END IF
C
               WIJ(:,:) = WIJ(:,:) - XINV(:,:)
               WJI(:,:) = WJI(:,:) - XINV(:,:)
            END IF
C
            IF ( UVUNITY ) THEN
               USSQX_XBAR(:,:,:,:) = C0
               VSSQX_XBAR(:,:,:,:) = C0
               DO I = 1,NXM
                  USSQX_XBAR(I,I,:,:) = C1
                  VSSQX_XBAR(I,I,:,:) = C1
               END DO
            END IF
C
C----------------------------------------------------------------------
C  IJ:   u(i;+) tau(i,j;+) v(j,+)  -  u(i;-) tau(j,i;+)^(T*) v(j,-)
C----------------------------------------------------------------------
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WIJ,NXM,
     &                 VSSQX_XBAR(1,1,JQ,1),NXM,C0,WKM3,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX_XBAR(1,1,IQ,1),NXM,
     &                 WKM3,NXM,C0,WKM1,NXM)
C
            CALL ZGEMM('C','N',NXM,NXM,NXM,C1,WJI,NXM,
     &                 VSSQX_XBAR(1,1,JQ,2),NXM,C0,WKM3,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX_XBAR(1,1,IQ,2),NXM,
     &                 WKM3,NXM,C0,WKM2,NXM)
C
Csw
            IF ( PLSONLY ) WKM2(:,:) = C0
Csw
            TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = WKM1(1:NXM,1:NXM)
     &         - WKM2(1:NXM,1:NXM)
C
Csw
            IF ( SSONLY ) TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = C0
            IF ( NOBSSD ) THEN
               IF ( IQ.EQ.JQ ) TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = C0
            ELSE IF ( ONLYBSSD ) THEN
               IF ( IQ.NE.JQ ) TAUDELTAK_BAR(IX1:IX2,JX1:JX2) = C0
            END IF
Csw
C----------------------------------------------------------------------
C  JI:   u(j;+) tau(j,i;+) v(i,+)  -  u(j;-) tau(i,j;+)^(T*) v(i,-)
C----------------------------------------------------------------------
            IF ( IQ.NE.JQ ) THEN
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WJI,NXM,
     &                    VSSQX_XBAR(1,1,IQ,1),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX_XBAR(1,1,JQ,1),
     &                    NXM,WKM3,NXM,C0,WKM1,NXM)
C
               CALL ZGEMM('C','N',NXM,NXM,NXM,C1,WIJ,NXM,
     &                    VSSQX_XBAR(1,1,IQ,2),NXM,C0,WKM3,NXM)
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX_XBAR(1,1,JQ,2),
     &                    NXM,WKM3,NXM,C0,WKM2,NXM)
C
Csw
               IF ( PLSONLY ) WKM2(:,:) = C0
Csw
               TAUDELTAK_BAR(JX1:JX2,IX1:IX2) = WKM1(1:NXM,1:NXM)
     &            - WKM2(1:NXM,1:NXM)
C
Csw
               IF ( SSONLY ) TAUDELTAK_BAR(JX1:JX2,IX1:IX2) = C0
               IF ( ONLYBSSD ) TAUDELTAK_BAR(JX1:JX2,IX1:IX2) = C0
Csw
C----------------------------------------------------------------------
C                   single site contribution
C----------------------------------------------------------------------
            ELSE IF ( .NOT.BSONLY ) THEN
C
               CALL ZGEMM('N','N',NXM,NXM,NXM,C1,GSSQX_XBAR(1,1,IQ),NXM,
     &                    LMAT2,NXM,C0,WKM1,NXM)
C
               TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
     &            = TAUDELTAK_BAR(IX1:IX2,JX1:JX2) + WKM1(1:NXM,1:NXM)
C
C----------------------------------------------------------------------
            END IF
C
C===================================================  extract pure gamma
C
C----------------------------- copy block of interest from TAUDELTAK_BAR
            XGAMMAKX(1:NXM,1:NXM) = TAUDELTAK_BAR(IX1:IX2,JX1:JX2)
C
C------------------------------------------------------------- get (1/x)
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,IQTB),WKM1,WKM2)
            CALL CMATINV3(NXM,NKMMAX,IPIVKM,DSSQX(1,1,JQTB),WKM1,WKM3)
C
C--------- multiply from left & right with 1/x^nu & L (1/x^nu')^dagger L
Cleft
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM2,NXM,XGAMMAKX(1,1),
     &                 NXM,C0,WKM1,NXM)
Cright
            CALL ZGEMM('C','N',NXM,NXM,NXM,C1,WKM3,NXM,LMAT2,NXM,C0,
     &                 WKM2,NXM)
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,LMAT2,NXM,WKM2,NXM,C0,
     &                 WKM3,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WKM1,NXM,WKM3,NXM,C0,
     &                 GAMMAK,NXM)
C
            IF ( IK.EQ.IKTOP ) THEN
               WRITE (43,*) 'IQ =',IQ,'JQ =',JQ
               CALL CMATSTRUCT('GAMMAK',GAMMAK,18,18,2,2,1,1.D-8,43)
            END IF
C
C===================================================  extract pure gamma
C
         END DO
      END DO
Csw
      IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &     CALL WR_GORTAUK_ZPM('UGV_BAR',114,TAUDELTAK_BAR,EXIKDQ_QTB,
     &     LMAT2,CCZ)
Csw
C
      END
C*==negfkloop_gles.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE NEGFKLOOP_GLES(GLESQX,TAUDELTAK,UGV_BAR,NKKR_TB,
     &                          EXIKDQ_QTB,WMPNP,WMN,BAR_FLAG_Q,
     &                          LBAR_FLAG_Q,RBAR_FLAG_Q,CCZ,LMAT2,IK,
     &                          IKTOP,CHECK_ZPM,WK,WKSUM)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_ANGMOM,ONLY:NXM,WKM1,WKM2
      USE MOD_SITES,ONLY:IQBOT_TB,IQTOP_TB,NQ,NQTB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CCZ,CHECK_ZPM
      INTEGER IK,IKTOP,NKKR_TB
      REAL*8 WK,WKSUM
      LOGICAL BAR_FLAG_Q(NQ),LBAR_FLAG_Q(NQ),RBAR_FLAG_Q(NQ)
      COMPLEX*16 EXIKDQ_QTB(NQTB,NQTB),GLESQX(NXM,NXM,NQTB),
     &           LMAT2(NXM,NXM),TAUDELTAK(NKKR_TB,NKKR_TB),
     &           UGV_BAR(NKKR_TB,NKKR_TB),WMN(NXM,NXM),WMPNP(NXM,NXM)
C
C Local variables
C
      INTEGER IIQ,IIQTB,IIX0,IIX1,IIX2,IQ,IQTB,IX0,IX1,IX2,JJQ,JJQTB,
     &        JJX0,JJX1,JJX2,JQ,JQTB,JX0,JX1,JX2
C
C*** End of declarations rewritten by SPAG
C
C======================================================================
C
C------------------------------------------------------------------- mu'
      LOOP_JJQ:DO JJQ = IQBOT_TB,IQTOP_TB
         IF ( BAR_FLAG_Q(JJQ) ) CYCLE LOOP_JJQ
         JJQTB = JJQ - IQBOT_TB + 1
         JJX0 = (JJQTB-1)*NXM
         JJX1 = JJX0 + 1
         JJX2 = JJX0 + NXM
C
C-------------------------------------------------------------------- mu
         LOOP_IIQ:DO IIQ = JJQ,JJQ
            IF ( BAR_FLAG_Q(IIQ) ) CYCLE LOOP_IIQ
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
C-------------------------------------------------------------------- nu
               LOOP_IQ:DO IQ = IQBOT_TB,IQTOP_TB
                  IF ( .NOT.((LBAR_FLAG_Q(IQ) .AND. LBAR_FLAG_Q(JQ))
     &                 .OR. (RBAR_FLAG_Q(IQ) .AND. RBAR_FLAG_Q(JQ))) )
     &                 CYCLE LOOP_IQ
                  IQTB = IQ - IQBOT_TB + 1
                  IX0 = (IQTB-1)*NXM
                  IX1 = IX0 + 1
                  IX2 = IX0 + NXM
C
                  WMN(1:NXM,1:NXM) = TAUDELTAK(IIX1:IIX2,IX1:IX2)
     &                               *EXIKDQ_QTB(IIQTB,IQTB)
C
                  WKM1(1:NXM,1:NXM) = UGV_BAR(IX1:IX2,JX1:JX2)
C
                  CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WMN,NXM,WKM1,NXM,C0,
     &                       WKM2,NXM)
C
                  CALL ZGEMM('N','C',NXM,NXM,NXM,C1,WKM2,NXM,WMPNP,NXM,
     &                       C0,WKM1,NXM)
C
                  GLESQX(1:NXM,1:NXM,IIQTB) = GLESQX(1:NXM,1:NXM,IIQTB)
     &               + WKM1(1:NXM,1:NXM)*WK/WKSUM
Csw     &               + CI*WKM1(1:NXM,1:NXM)*WK/WKSUM
CSW-- no factor CI: we have it pure real (on the L-diagonal) like MEZZMP
CSW----- either leave it like that or "apply CI" twice at correct places
CSW----------------------------------- + sign should have to be switched
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
Csw
      IF ( CHECK_ZPM .AND. IK.EQ.IKTOP )
     &     CALL WR_GORTAUK_ZPM('TAUDELTAK',115,TAUDELTAK,EXIKDQ_QTB,
     &     LMAT2,CCZ)
Csw
C
      END
C*==wr_gortauk_zpm.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE WR_GORTAUK_ZPM(STR,IFIL,GORTAUK,EXIKDQ_QTB,LMAT2,CCZ)
C   ********************************************************************
C   *                                                                  *
C   *   Write out GREFKE, GLLKE, or TAUDELTAK block-wise for z & z*    *
C   *   to check wether                                                *
C   *                                                                  *
C   *     X_LL'^nn' (z) = L_LL' * (X_LL'^n'n (z*))^dagger * L_LL'      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_ANGMOM,ONLY:NXM
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
