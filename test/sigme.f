C*==sigme.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGME(ILOOP_SIGMA_PROJECT,ERYDA,ERYDACC,ERYDB,ERYDBCC,
     &                 IFILA,IFILACC,IFILB,IFILBCC,MAQAB,MAQBA,MBQAB,
     &                 MBQBA,MCQAB,MCQBA,MDQAB,MDQBA,MREG_TAB,MREG_TBA,
     &                 S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,S10AQAB,
     &                 S10BQAB,S10CQAB,S10DQAB,S2AQAB,S2BQAB,S2CQAB,
     &                 S2DQAB,S3AQAB,S3BQAB,S3CQAB,S3DQAB,S4AQAB,S4BQAB,
     &                 S4CQAB,S4DQAB,UMAT_VTA,MSSTA,TSSTA,MSSQA,
     &                 TAUQA_TMP,UMAT_VTB,MSSTB,TSSTB,MSSQB,TAUQB_TMP,
     &                 NSPINPROJ,MEOFF_BARG_T,MEOFF_AHE_T,ME_SZ_T)
C   ********************************************************************
C   *                                                                  *
C   *   - calculate the SIGMA matrix elements                          *
C   *   - transform the matrix elements for different E-arguments      *
C   *                                                                  *
C   *   the l-expansion is set via   NKM   for all atomic types to NL  *
C   *                                                                  *
C   *   ALL matrix elements are                                        *
C   *   - calculated in the local frame of reference                   *
C   *   - the polarisation is transferred to cartesian coordinates     *
C   *   - stored as type IT dependent quantities                       *
C   *   - the transformation according to thermal fluctuations         *
C   *     and/or vibrations is done ON THE FLY                         *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *  KMROT                                                           *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site IQ   *
C   *  2: global COMMON rotation of the magnetisation                  *
C   *  3: spin spiral    Theta =  90    not allowed                    *
C   *  4: spin spiral    Theta <> 90    not allowed                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,
c modified by XJQ: scfvb
     &    lrecreal8
c end-mod-xjq
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,WKM2,LMAT3,IPIVKM
      USE MOD_SITES,ONLY:NQMAX,MROTQ,DROTQ,IQAT,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:CTL,NTMAX,ITBOT,ITTOP
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_SIG,ONLY:LCHECKME_SPIN,NSPR,FLAG_ISPINPROJ,LIST_ISPR
      USE MOD_CALCMODE,ONLY:MOMENTS_ROTATED
      USE MOD_THERMAL,ONLY:NVT,NFT,FMAT_FT,MFMAT_FT,NVIBFLU,IQ_AVFT,
     &    NVIBRA,NFLUCT,NVFTMAX,
c modified by XJQ: scfvb
     &    i_temp_lat, n_temp_lat
      use mod_scfvb_cpa_sigma,only:lscfvb,ifilscftv,ifilmreg,ifilmirr
c end-mod-xjq
      IMPLICIT NONE
C*--SIGME46
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGME')
      LOGICAL CHECK,CCF
      PARAMETER (CHECK=.FALSE.,CCF=.FALSE.)
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      COMPLEX*16 ERYDA,ERYDACC,ERYDB,ERYDBCC
      INTEGER IFILA,IFILACC,IFILB,IFILBCC,ILOOP_SIGMA_PROJECT,NSPINPROJ
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MEOFF_AHE_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           MEOFF_BARG_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           ME_SZ_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2),
     &           MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2),
     &           MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           S10AQAB(3,3,NSPINPROJ,NQMAX),
     &           S10BQAB(3,3,NSPINPROJ,NQMAX),
     &           S10CQAB(3,3,NSPINPROJ,NQMAX),
     &           S10DQAB(3,3,NSPINPROJ,NQMAX),
     &           S10_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S2AQAB(3,3,NSPINPROJ,NQMAX),S2BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S2CQAB(3,3,NSPINPROJ,NQMAX),
     &           S2DQAB(3,3,NSPINPROJ,NQMAX),
     &           S2_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S3AQAB(3,3,NSPINPROJ,NQMAX),S3BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S3CQAB(3,3,NSPINPROJ,NQMAX),
     &           S3DQAB(3,3,NSPINPROJ,NQMAX),
     &           S3_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S4AQAB(3,3,NSPINPROJ,NQMAX),S4BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4CQAB(3,3,NSPINPROJ,NQMAX),
     &           S4DQAB(3,3,NSPINPROJ,NQMAX),
     &           S4_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           TAUQA_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TAUQB_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TSSTA(NKMMAX,NKMMAX,NTMAX),TSSTB(NKMMAX,NKMMAX,NTMAX),
     &           UMAT_VTA(NKMMAX,NKMMAX,NVT),UMAT_VTB(NKMMAX,NKMMAX,NVT)
C
C Local variables
C
c modified by XJQ: scfvb
      integer reclng_tmat, reclng_mreg, reclng_mirr
c end-mod-xjq
      REAL*8 C
      COMPLEX*16 CMATTRC
      COMPLEX*16 CSUM,EATMP,EBTMP,MIRR_2(:,:,:,:,:),
     &           MIRR_2_FTAB(:,:,:,:,:),MIRR_2_TMP(:,:,:,:),
     &           MIRR_2_VFTAB(:,:,:,:,:),MIRR_3(:,:,:,:,:),
     &           MIRR_3_FTAB(:,:,:,:,:),MIRR_3_TMP(:,:,:,:),
     &           MIRR_3_VFTAB(:,:,:,:,:),MIRR_4(:,:,:,:,:),
     &           MIRR_4_FTAB(:,:,:,:,:),MIRR_4_TMP(:,:,:,:),
     &           MIRR_4_VFTAB(:,:,:,:,:),MREG_FTAB(:,:,:,:),
     &           MREG_FTBA(:,:,:,:),MREG_TAB_TMP(:,:,:),
     &           MREG_TBA_TMP(:,:,:),MREG_VFTAB(:,:,:,:),
     &           MREG_VFTBA(:,:,:,:),MSS_VFTA(:,:),MSS_VFTB(:,:),
     &           TAU_VFTA(:,:),TAU_VFTAZ(:,:,:),TAU_VFTB(:,:),
     &           TAU_VFTBZ(:,:,:),TSS_FTA(:,:),TSS_FTB(:,:),
     &           TSS_VFTA(:,:),TSS_VFTB(:,:),UMAT_VTAZ(:,:,:),
     &           UMAT_VTBZ(:,:,:)
      CHARACTER*50 FMT_CHECK
      INTEGER I,IA_ERR,IFILATMP,IFILBTMP,IFLUCT,IFT,IQ,ISP,IT,IVFT,
     &        IVIBRA,IVT,IXX,IZA,IZB,J,M,MUE,N,NUE
      LOGICAL IM_ERYDA_NE_0,IM_ERYDB_NE_0,INITIALIZE
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE MIRR_2,MIRR_3,MIRR_4
      ALLOCATABLE MREG_TAB_TMP,MREG_TBA_TMP
      ALLOCATABLE MIRR_2_TMP,MIRR_3_TMP,MIRR_4_TMP
      ALLOCATABLE TAU_VFTAZ,TAU_VFTBZ,UMAT_VTAZ,UMAT_VTBZ
      ALLOCATABLE TAU_VFTA,TAU_VFTB
      ALLOCATABLE TSS_FTA,TSS_VFTA,MSS_VFTA
      ALLOCATABLE TSS_FTB,TSS_VFTB,MSS_VFTB
      ALLOCATABLE MIRR_2_FTAB,MIRR_3_FTAB,MIRR_4_FTAB
      ALLOCATABLE MREG_FTAB,MREG_FTBA
      ALLOCATABLE MIRR_2_VFTAB,MIRR_3_VFTAB,MIRR_4_VFTAB
      ALLOCATABLE MREG_VFTAB,MREG_VFTBA
C
C
      CALL TRACK_INFO(ROUTINE)
C
c modified by XJQ: scfvb
c
      if(lscfvb .and. i_temp_lat==1) then
        reclng_tmat = 2*lrecreal8 * nkmmax * nkmmax
        reclng_mreg = 2*lrecreal8 * nkmmax * nkmmax * 3 * 2
        reclng_mirr = 2*lrecreal8 * nkmmax * nkmmax * 3 * 3 * 3
        open(unit=ifilscftv,file='scftv.dat',form='unformatted',
     &       access='direct',recl=reclng_tmat,status='old')
        open(unit=ifilmreg,file='mreg.dat',form='unformatted',
     &       access='direct',recl=reclng_mreg,status='old')
        open(unit=ifilmirr,file='mirr.dat',form='unformatted',
     &       access='direct',recl=reclng_mirr,status='old')
      endif
c end-mod-xjq
c
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FMT_CHECK = '(///,''CHECK '',50(''>''),4X,A,//,(2F26.12))'
      IF ( CHECK ) THEN
         WRITE (6,FMT=FMT_CHECK) 'UMAT_VTA   ',UMAT_VTA
         WRITE (6,FMT=FMT_CHECK) 'MSSTA      ',MSSTA
         WRITE (6,FMT=FMT_CHECK) 'TSSTA      ',TSSTA
         WRITE (6,FMT=FMT_CHECK) 'MSSQA      ',MSSQA
         WRITE (6,FMT=FMT_CHECK) 'TAUQA_TMP  ',TAUQA_TMP
C
         WRITE (6,FMT=FMT_CHECK) 'UMAT_VTB   ',UMAT_VTB
         WRITE (6,FMT=FMT_CHECK) 'MSSTB      ',MSSTB
         WRITE (6,FMT=FMT_CHECK) 'TSSTB      ',TSSTB
         WRITE (6,FMT=FMT_CHECK) 'MSSQB      ',MSSQB
         WRITE (6,FMT=FMT_CHECK) 'TAUQB_TMP  ',TAUQB_TMP
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C=======================================================================
C     Checks on MEs
C=======================================================================
C
      IF ( CCF ) CALL STOP_MESSAGE(ROUTINE,
     &                             'CCF = .TRUE. for TRANSPORT ?!?!')
C
      IF ( LCHECKME_SPIN ) WRITE (6,99003)
C
      M = NKMMAX
      N = NTMAX
C
      ALLOCATE (TSS_FTA(NKMMAX,NKMMAX),TSS_FTB(NKMMAX,NKMMAX))
      ALLOCATE (TSS_VFTA(NKMMAX,NKMMAX),TSS_VFTB(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFTA(NKMMAX,NKMMAX),MSS_VFTB(NKMMAX,NKMMAX))
      ALLOCATE (MIRR_2_FTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_3_FTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_4_FTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MREG_FTAB(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MREG_FTBA(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MIRR_2_VFTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_3_VFTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MIRR_4_VFTAB(NKMMAX,NKMMAX,3,3,NSPINPROJ))
      ALLOCATE (MREG_VFTAB(NKMMAX,NKMMAX,3,NSPINPROJ))
      ALLOCATE (MREG_VFTBA(NKMMAX,NKMMAX,3,NSPINPROJ))
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C  rotation from local to global frame for non-collinear magnetisation
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      IF ( MOMENTS_ROTATED ) THEN
C
         ALLOCATE (MREG_TAB_TMP(NKMMAX,NKMMAX,3))
         ALLOCATE (MREG_TBA_TMP(NKMMAX,NKMMAX,3))
         ALLOCATE (MIRR_2_TMP(NKMMAX,NKMMAX,3,3))
         ALLOCATE (MIRR_3_TMP(NKMMAX,NKMMAX,3,3))
         ALLOCATE (MIRR_4_TMP(NKMMAX,NKMMAX,3,3))
C
      END IF
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C ------------------------------------ calculate angular matrix elements
C
         CALL AME_ALF_SIG
C
         CALL AME_NAB_SIG
C
         CALL AME_NAB
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C                               INITIALIZE                           END
C=======================================================================
C
      MAQAB(:,:,:,:,:) = C0
      MAQBA(:,:,:,:,:) = C0
      MBQAB(:,:,:,:,:) = C0
      MBQBA(:,:,:,:,:) = C0
      MCQAB(:,:,:,:,:) = C0
      MCQBA(:,:,:,:,:) = C0
      MDQAB(:,:,:,:,:) = C0
      MDQBA(:,:,:,:,:) = C0
C
      S10AQAB(:,:,:,:) = C0
      S10BQAB(:,:,:,:) = C0
      S10CQAB(:,:,:,:) = C0
      S10DQAB(:,:,:,:) = C0
C
      S2AQAB(:,:,:,:) = C0
      S2BQAB(:,:,:,:) = C0
      S2CQAB(:,:,:,:) = C0
      S2DQAB(:,:,:,:) = C0
C
      S3AQAB(:,:,:,:) = C0
      S3BQAB(:,:,:,:) = C0
      S3CQAB(:,:,:,:) = C0
      S3DQAB(:,:,:,:) = C0
C
      S4AQAB(:,:,:,:) = C0
      S4BQAB(:,:,:,:) = C0
      S4CQAB(:,:,:,:) = C0
      S4DQAB(:,:,:,:) = C0
C
      CALL CINIT(M*M*3*3*NTMAX,MEOFF_BARG_T)
      CALL CINIT(M*M*3*3*NTMAX,MEOFF_AHE_T)
      CALL CINIT(M*M*3*3*NTMAX,ME_SZ_T)
C
C=======================================================================
C
      M = NKMMAX
C
      ALLOCATE (TAU_VFTA(NKMMAX,NKMMAX),TAU_VFTB(NKMMAX,NKMMAX))
      ALLOCATE (MIRR_2(M,M,3,3,NSPINPROJ))
      ALLOCATE (MIRR_3(M,M,3,3,NSPINPROJ))
      ALLOCATE (MIRR_4(M,M,3,3,NSPINPROJ))
      ALLOCATE (UMAT_VTAZ(NKMMAX,NKMMAX,2))
      ALLOCATE (UMAT_VTBZ(NKMMAX,NKMMAX,2))
      ALLOCATE (TAU_VFTAZ(NKMMAX,NKMMAX,2))
      ALLOCATE (TAU_VFTBZ(NKMMAX,NKMMAX,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate MREG_T')
C
      S10_VFTAB(:,:,:,:,:,:) = C0
      S2_VFTAB(:,:,:,:,:,:) = C0
      S3_VFTAB(:,:,:,:,:,:) = C0
      S4_VFTAB(:,:,:,:,:,:) = C0
C
C--------------------------------------------------- work on real E-axis
      IM_ERYDA_NE_0 = .FALSE.
      IM_ERYDB_NE_0 = .FALSE.
C----------------------------------------------- work in complex E-plane
      IF ( ABS(DIMAG(ERYDA-ERYDACC)).GT.TOL ) IM_ERYDA_NE_0 = .TRUE.
      IF ( ABS(DIMAG(ERYDB-ERYDBCC)).GT.TOL ) IM_ERYDB_NE_0 = .TRUE.
C
      IF ( IM_ERYDA_NE_0 .OR. IM_ERYDB_NE_0 ) WRITE (6,99002) ROUTINE
C
      MREG_TAB(:,:,:,:,:,:,:) = C0
      MREG_TBA(:,:,:,:,:,:,:) = C0
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
         IF ( IQ.LT.IQBOT_CHI .OR. IQ.GT.IQTOP_CHI ) CYCLE LOOP_IT
C
         M = NKMMAX
         N = NKM
C
         IF ( IPRINT.GE.3 ) WRITE (6,99001) IT
C
         C = CTL(IT,1)
C
CIZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA
         LOOP_IZA:DO IZA = 1,2
C
            IF ( IZA.EQ.1 ) THEN
               EATMP = ERYDACC
               IFILATMP = IFILACC
            ELSE
               EATMP = ERYDA
               IFILATMP = IFILA
            END IF
C
CIZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB
            LOOP_IZB:DO IZB = 1,2
C
               IF ( IZB.EQ.1 ) THEN
                  EBTMP = ERYDBCC
                  IFILBTMP = IFILBCC
               ELSE
                  EBTMP = ERYDB
                  IFILBTMP = IFILB
               END IF
C
C-----------------------------------------------------------------------
C     MEs for current density:                        alpha-alpha form
C-----------------------------------------------------------------------
C
               IXX = 1
C
               CALL ME_ALF_ALF(1,NKM,IFILATMP,EATMP,1,NKM,IFILBTMP,
     &                         EBTMP,CCF,IT,
     &                         MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &                         MREG_TBA(1,1,1,IXX,IT,IZB,IZA),
     &                         MIRR_2(1,1,1,1,IXX),MIRR_3(1,1,1,1,IXX),
     &                         MIRR_4(1,1,1,1,IXX),C,K_CALC_ME)
C
C-----------------------------------------------------------------------
C     MEs for Bargmann-Wigner Polarization operator:  alpha-alpha form
C-----------------------------------------------------------------------
C
               IXX = 2
C
               IF ( FLAG_ISPINPROJ(IXX) )
     &              CALL ME_SPIN_ALF_ALF(IFILATMP,EATMP,IFILBTMP,EBTMP,
     &              IT,MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &              MREG_TBA(1,1,1,IXX,IT,IZB,IZA),MIRR_2(1,1,1,1,IXX),
     &              MIRR_3(1,1,1,1,IXX),MIRR_4(1,1,1,1,IXX),C)
C
C-----------------------------------------------------------------------
C     MEs for Bargmann-Wigner Polarization operator:  nabla-alpha form
C-----------------------------------------------------------------------
C
               IXX = 5
C
               IF ( FLAG_ISPINPROJ(IXX) )
     &              CALL ME_SPIN_NAB_ALF(IFILATMP,EATMP,IFILBTMP,EBTMP,
     &              IT,MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &              MREG_TBA(1,1,1,IXX,IT,IZB,IZA),MIRR_2(1,1,1,1,IXX),
     &              MIRR_3(1,1,1,1,IXX),MIRR_4(1,1,1,1,IXX),C)
C
C-----------------------------------------------------------------------
C     MEs for SOT spin-orbit torque
C-----------------------------------------------------------------------
C
               IXX = 8
C
               IF ( FLAG_ISPINPROJ(IXX) )
     &              CALL ME_SOT_ALF(IFILATMP,EATMP,IFILBTMP,EBTMP,IT,
     &              MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &              MREG_TBA(1,1,1,IXX,IT,IZB,IZA),MIRR_2(1,1,1,1,IXX),
     &              MIRR_3(1,1,1,1,IXX),MIRR_4(1,1,1,1,IXX),C)
C
C-----------------------------------------------------------------------
C     MEs for menab
C-----------------------------------------------------------------------
C
               IXX = 9
C
               IF ( FLAG_ISPINPROJ(IXX) )
     &              CALL ME_NAB_NAB(1,NKM,IFILATMP,EATMP,1,NKM,IFILBTMP,
     &              EBTMP,CCF,IT,MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &              MREG_TBA(1,1,1,IXX,IT,IZB,IZA),MIRR_2(1,1,1,1,IXX),
     &              MIRR_3(1,1,1,1,IXX),MIRR_4(1,1,1,1,IXX),C,K_CALC_ME)
C
C-----------------------------------------------------------------------
C     MEs for SIGMA (Edelstein effect)
C-----------------------------------------------------------------------
C
               IXX = 10
C
               IF ( FLAG_ISPINPROJ(IXX) )
     &              CALL ME_SIG_ALF(IFILATMP,EATMP,IFILBTMP,EBTMP,IT,
     &              MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &              MREG_TBA(1,1,1,IXX,IT,IZB,IZA),MIRR_2(1,1,1,1,IXX),
     &              MIRR_3(1,1,1,1,IXX),MIRR_4(1,1,1,1,IXX),C)
C
C-----------------------------------------------------------------------
C
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
               IF ( IZA.EQ.2 .AND. IZB.EQ.2 .AND. 1.EQ.0 )
     &              CALL SIGME_DEBUG(1,IT,MREG_TAB(1,1,1,1,IT,IZA,IZB),
     &                               MREG_TBA(1,1,1,1,IT,IZB,IZA),
     &                               MIRR_2,MIRR_3,MIRR_4,NSPINPROJ)
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
C
C
C=======================================================================
C             transform and renormalize matrix elements
C=======================================================================
C
               LOOP_IXX_A:DO ISP = 1,NSPR
C
                  IXX = LIST_ISPR(ISP)
C
C-----------------------------------------------------------------------
C           transformation from spherical to cartesian coordinates
C              from now on        ipol = 1,2,3  ==  (x),(y),(z)
C-----------------------------------------------------------------------
C
                  CALL CMAT_CONVERT_POLAR(MREG_TAB(1,1,1,IXX,IT,IZA,IZB)
     &               ,'S>C')
C
                  CALL CMAT_CONVERT_POLAR(MREG_TBA(1,1,1,IXX,IT,IZB,IZA)
     &               ,'S>C')
C
                  CALL CMAT_CONVERT_POLAR2(MIRR_2(1,1,1,1,IXX),'S>C')
                  CALL CMAT_CONVERT_POLAR2(MIRR_3(1,1,1,1,IXX),'S>C')
                  CALL CMAT_CONVERT_POLAR2(MIRR_4(1,1,1,1,IXX),'S>C')
C
C???????????????????????????????????????????????????????????????????????
                  IF ( ISP.LT.0 ) THEN
C
                     WRITE (*,*) '####################### ISP, IXX',ISP,
     &                           IXX
                     IQ = 100
                     WRITE (IQ,*) '####################### ISP IXX',ISP,
     &                            IXX
                     WRITE (IQ,*) 'MZAZB'
                     WRITE (IQ,'(2e22.14)')
     &                      MREG_TAB(:,:,:,IXX,IT,IZA,IZB)
                     WRITE (IQ,*) 'MZBZA'
                     WRITE (IQ,'(2e22.14)')
     &                      MREG_TBA(:,:,:,IXX,IT,IZB,IZA)
                     WRITE (IQ,*) 'MIRR_2'
                     WRITE (IQ,'(2e22.14)') MIRR_2(:,:,:,:,IXX)
                     WRITE (IQ,*) 'MIRR_3'
                     WRITE (IQ,'(2e22.14)') MIRR_3(:,:,:,:,IXX)
                     WRITE (IQ,*) 'MIRR_4'
                     WRITE (IQ,'(2e22.14)') MIRR_4(:,:,:,:,IXX)
                     IF ( ISP.EQ.8 ) STOP
                     CYCLE LOOP_IXX_A
C
                  END IF
C???????????????????????????????????????????????????????????????????????
C
C-----------------------------------------------------------------------
C  rotation from local to global frame for non-collinear magnetisation
C-----------------------------------------------------------------------
C
                  IF ( MOMENTS_ROTATED ) THEN
C
                     IQ = IQAT(1,IT)
C
                     MREG_TAB_TMP(:,:,:)
     &                  = MREG_TAB(:,:,:,IXX,IT,IZA,IZB)
                     MREG_TBA_TMP(:,:,:)
     &                  = MREG_TBA(:,:,:,IXX,IT,IZB,IZA)
C
                     MIRR_2_TMP(:,:,:,:) = MIRR_2(:,:,:,:,IXX)
                     MIRR_3_TMP(:,:,:,:) = MIRR_3(:,:,:,:,IXX)
                     MIRR_4_TMP(:,:,:,:) = MIRR_4(:,:,:,:,IXX)
C
                     CALL ME_ROTATE_REG(IQ,NQMAX,DROTQ,MROTQ,
     &                                  MREG_TAB_TMP,MREG_TBA_TMP,
     &                                  MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &                                  MREG_TBA(1,1,1,IXX,IT,IZB,IZA))
C
                     CALL ME_ROTATE_IRR(IQ,NQMAX,DROTQ,MROTQ,MIRR_2_TMP,
     &                                  MIRR_3_TMP,MIRR_4_TMP,
     &                                  MIRR_2(1,1,1,1,IXX),
     &                                  MIRR_3(1,1,1,1,IXX),
     &                                  MIRR_4(1,1,1,1,IXX))
C
                  END IF
C
C----------------------------------------------------------------------
C   in case of energy on real axis multiply by appropriate
C   phase factors (LMAT) to obtain MEs
C   M(E+,E-), M(E-,E+) and M(E-,E-)
C----------------------------------------------------------------------
C
                  DO MUE = 1,3
C
                     IF ( (.NOT.IM_ERYDA_NE_0) .AND. IZA.EQ.1 ) THEN
C
                        CALL CMATMUL(N,M,LMAT3,
     &                               MREG_TAB(1,1,MUE,IXX,IT,IZA,IZB),
     &                               WKM1)
                        MREG_TAB(:,:,MUE,IXX,IT,IZA,IZB) = WKM1(:,:)
C
                        CALL CMATMUL(N,M,
     &                               MREG_TBA(1,1,MUE,IXX,IT,IZB,IZA),
     &                               LMAT3,WKM1)
                        MREG_TBA(:,:,MUE,IXX,IT,IZB,IZA) = WKM1(:,:)
C
                        DO NUE = 1,3
C
                           CALL CMATMUL(N,M,LMAT3,
     &                                  MIRR_3(1,1,MUE,NUE,IXX),WKM1)
C
                           CALL CMATMUL(N,M,WKM1,LMAT3,
     &                                  MIRR_3(1,1,MUE,NUE,IXX))
C
                        END DO
C
                     END IF
C
                     IF ( (.NOT.IM_ERYDB_NE_0) .AND. IZB.EQ.1 ) THEN
C
                        CALL CMATMUL(N,M,
     &                               MREG_TAB(1,1,MUE,IXX,IT,IZA,IZB),
     &                               LMAT3,WKM1)
                        MREG_TAB(:,:,MUE,IXX,IT,IZA,IZB) = WKM1(:,:)
C
                        CALL CMATMUL(N,M,LMAT3,
     &                               MREG_TBA(1,1,MUE,IXX,IT,IZB,IZA),
     &                               WKM1)
                        MREG_TBA(:,:,MUE,IXX,IT,IZB,IZA) = WKM1(:,:)
C
                        DO NUE = 1,3
C
                           CALL CMATMUL(N,M,LMAT3,
     &                                  MIRR_2(1,1,MUE,NUE,IXX),WKM1)
C
                           CALL CMATMUL(N,M,WKM1,LMAT3,
     &                                  MIRR_2(1,1,MUE,NUE,IXX))
C
                        END DO
C
                     END IF
C
                  END DO
C----------------------------------------------------------------------
C
               END DO LOOP_IXX_A
c modified by XJQ: scfvb
c
c               if(lscfvb)
c     &           call checkmregmirr_noscfvb(iza,izb,it,
c     &             mreg_tab(:,:,:,1,it,iza,izb),
c     &             mreg_tba(:,:,:,1,it,izb,iza),
c     &             mirr_2(:,:,:,:,1),mirr_3(:,:,:,:,1),
c     &             mirr_4(:,:,:,:,1))
c
c end-mod-xjq
C=======================================================================
C
C
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
               IF ( IZA.EQ.1 .AND. IZB.EQ.1 .AND. 1.EQ.0 )
     &              CALL SIGME_DEBUG(1,IT,MREG_TAB(1,1,1,1,IT,IZA,IZB),
     &                               MREG_TBA(1,1,1,1,IT,IZB,IZA),
     &                               MIRR_2,MIRR_3,MIRR_4,NSPINPROJ)
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
C
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C  manipulate the matrix elements to get SPIN-projected conductivities
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
               IF ( ILOOP_SIGMA_PROJECT.NE.0 ) THEN
                  CALL STOP_MESSAGE(ROUTINE,
     &       'SPIN-projected conductivities not supported at the moment'
     &       )
C
                  DO MUE = 1,3
C
                     CALL SIGMANIPUL(MREG_TAB(1,1,MUE,IXX,IT,IZA,IZB),
     &                               ILOOP_SIGMA_PROJECT)
                     CALL SIGMANIPUL(MREG_TBA(1,1,MUE,IXX,IT,IZB,IZA),
     &                               ILOOP_SIGMA_PROJECT)
C
C                  DO NUE = 1,3
C                  CALL SIGMANIPUL(MIRRTAB(1,1,1,MUE,NUE,IXX,IT),
C     &                    ILOOP_SIGMA_PROJECT)
C                  END DO
C
                  END DO
               END IF
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
               IFT = (IT-1)*NFLUCT
               LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
                  IF ( NFLUCT.EQ.1 ) THEN
C
                     TSS_FTA(:,:) = TSSTA(:,:,IT)
                     TSS_FTB(:,:) = TSSTB(:,:,IT)
C
                     DO ISP = 1,NSPR
                        IXX = LIST_ISPR(ISP)
                        MREG_FTAB(:,:,:,IXX)
     &                     = MREG_TAB(:,:,:,IXX,IT,IZA,IZB)
                        MREG_FTBA(:,:,:,IXX)
     &                     = MREG_TBA(:,:,:,IXX,IT,IZB,IZA)
                        MIRR_2_FTAB(:,:,:,:,IXX) = MIRR_2(:,:,:,:,IXX)
                        MIRR_3_FTAB(:,:,:,:,IXX) = MIRR_3(:,:,:,:,IXX)
                        MIRR_4_FTAB(:,:,:,:,IXX) = MIRR_4(:,:,:,:,IXX)
                     END DO
C
                  ELSE
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSSTA(1,1,IT),'SAS+',TSS_FTA)
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSSTB(1,1,IT),'SAS+',TSS_FTB)
C
                     DO ISP = 1,NSPR
                        IXX = LIST_ISPR(ISP)
C
                        CALL ME_ROTATE_REG(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                     MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &                     MREG_TBA(1,1,1,IXX,IT,IZB,IZA),
     &                     MREG_FTAB(1,1,1,IXX),MREG_FTBA(1,1,1,IXX))
C
                        CALL ME_ROTATE_IRR(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                     MIRR_2(1,1,1,1,IXX),MIRR_3(1,1,1,1,IXX),
     &                     MIRR_4(1,1,1,1,IXX),MIRR_2_FTAB(1,1,1,1,IXX),
     &                     MIRR_3_FTAB(1,1,1,1,IXX),
     &                     MIRR_4_FTAB(1,1,1,1,IXX))
C
                     END DO
C
                  END IF
C
C============================================================= IVIBRA ==
C                                             perform local displacement
                  IVT = (IT-1)*NVIBRA
                  LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                     IVT = IVT + 1
C
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
C
                     IF ( NVIBRA.EQ.1 ) THEN
C
                        TSS_VFTA(:,:) = TSS_FTA(:,:)
                        TSS_VFTB(:,:) = TSS_FTB(:,:)
C
                        DO ISP = 1,NSPR
                           IXX = LIST_ISPR(ISP)
                           MREG_VFTAB(:,:,:,IXX) = MREG_FTAB(:,:,:,IXX)
                           MREG_VFTBA(:,:,:,IXX) = MREG_FTBA(:,:,:,IXX)
                           MIRR_2_VFTAB(:,:,:,:,IXX)
     &                        = MIRR_2_FTAB(:,:,:,:,IXX)
                           MIRR_3_VFTAB(:,:,:,:,IXX)
     &                        = MIRR_3_FTAB(:,:,:,:,IXX)
                           MIRR_4_VFTAB(:,:,:,:,IXX)
     &                        = MIRR_4_FTAB(:,:,:,:,IXX)
                        END DO
C
                     ELSE
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C---------------------------- NB: umat and t at energies above real axis
C
c modified by XJQ: scfvb
                        call read_scftv(ivt,TSS_FTA)
c end-mod-xjq
                        CALL CMAT_U_TRANS(UMAT_VTA(1,1,IVT),TSS_FTA,
     &                     'UAUT',TSS_VFTA)
                        CALL CMAT_U_TRANS(UMAT_VTB(1,1,IVT),TSS_FTB,
     &                     'UAUT',TSS_VFTB)
C
C------------------------------- MREG_VFT: M_vft = U_v * M_ft * U_v^(-1)
C
                        CALL GET_MATZ_UMAT(UMAT_VTA(1,1,IVT),UMAT_VTAZ)
                        CALL GET_MATZ_UMAT(UMAT_VTB(1,1,IVT),UMAT_VTBZ)
C
                        DO ISP = 1,NSPR
                           IXX = LIST_ISPR(ISP)
c modified by XJQ: scfvb
                           call read_mregmirr_scfvb(iza,izb,ivt,
     &                       mreg_ftab(:,:,:,ixx),mreg_ftba(:,:,:,ixx),
     &                       mirr_2_ftab(:,:,:,:,ixx),
     &                       mirr_3_ftab(:,:,:,:,ixx),
     &                       mirr_4_ftab(:,:,:,:,ixx))
c end-mod-xjq
                           CALL ME_UTRANS_REG(UMAT_VTAZ(1,1,IZA),
     &                        UMAT_VTBZ(1,1,IZB),MREG_FTAB(1,1,1,IXX),
     &                        MREG_FTBA(1,1,1,IXX),MREG_VFTAB(1,1,1,IXX)
     &                        ,MREG_VFTBA(1,1,1,IXX))
C
                           CALL ME_UTRANS_IRR(UMAT_VTAZ(1,1,IZA),
     &                        UMAT_VTBZ(1,1,IZB),
     &                        MIRR_2_FTAB(1,1,1,1,IXX),
     &                        MIRR_3_FTAB(1,1,1,1,IXX),
     &                        MIRR_4_FTAB(1,1,1,1,IXX),
     &                        MIRR_2_VFTAB(1,1,1,1,IXX),
     &                        MIRR_3_VFTAB(1,1,1,1,IXX),
     &                        MIRR_4_VFTAB(1,1,1,1,IXX))
                        END DO
C
                     END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C------ NB: umat and t at energies above real axis --> m above real axis
C
                     IF ( NVIBFLU.EQ.1 ) THEN
                        MSS_VFTA(:,:) = MSSTA(:,:,IT)
                        MSS_VFTB(:,:) = MSSTB(:,:,IT)
                     ELSE
                        CALL CMATINV3(N,M,IPIVKM,TSS_VFTA,WKM1,MSS_VFTA)
                        CALL CMATINV3(N,M,IPIVKM,TSS_VFTB,WKM1,MSS_VFTB)
                     END IF
C
C-------------------------------------- get projected TAU matrix TAU_vft
C-------------------------------------------------- TAU(t) = TAU * D~(t)
C------- NB: m, tau at energy above real axis --> tau(t) above real axis
C
                     IQ = IQ_AVFT(1,IVFT)
C
                     CALL GET_TAUT(MSS_VFTA,MSSQA(1,1,IQ),
     &                             TAUQA_TMP(1,1,IQ),TAU_VFTA)
C
                     CALL GET_TAUT(MSS_VFTB,MSSQB(1,1,IQ),
     &                             TAUQB_TMP(1,1,IQ),TAU_VFTB)
C
C-----------------------------------------------------------------------
C perform proper transformation for z -> z*    TAU(z*) = L TAU(z)^+ L
C-----------------------------------------------------------------------
C
                     CALL GET_MATZ(TAU_VFTA,TAU_VFTAZ)
                     CALL GET_MATZ(TAU_VFTB,TAU_VFTBZ)
C
C-----------------------------------------------------------------------
C
                     LOOP_IXX:DO ISP = 1,NSPR
                        IXX = LIST_ISPR(ISP)
C
                        DO NUE = 1,3
                           DO MUE = 1,3
C
C------------------------------------------------------------------- S10
C
                              CALL CMATMUL(N,M,MREG_VFTBA(1,1,MUE,IXX),
     &                           TAU_VFTAZ(1,1,IZA),WKM1)
C
                              CALL CMATMUL(N,M,WKM1,
     &                           MREG_VFTAB(1,1,NUE,IXX),WKM2)
C
                              CALL CMATMUL(N,M,WKM2,TAU_VFTBZ(1,1,IZB),
     &                           WKM1)
C
                              S10_VFTAB(MUE,NUE,IXX,IVFT,IZA,IZB)
     &                           = CMATTRC(N,M,WKM1)
C
C-------------------------------------------------------------------- S2
C
                              CALL CMATMUL(N,M,
     &                           MIRR_2_VFTAB(1,1,MUE,NUE,IXX),
     &                           TAU_VFTBZ(1,1,IZB),WKM1)
C
                              S2_VFTAB(MUE,NUE,IXX,IVFT,IZA,IZB)
     &                           = -CMATTRC(N,M,WKM1)
C
C-------------------------------------------------------------------- S3
C
                              CALL CMATMUL(N,M,
     &                           TAU_VFTAZ(1:NKM,1:NKM,IZA),
     &                           MIRR_3_VFTAB(1,1,MUE,NUE,IXX),WKM1)
C
                              S3_VFTAB(MUE,NUE,IXX,IVFT,IZA,IZB)
     &                           = -CMATTRC(N,M,WKM1)
C
C-------------------------------------------------------------------- S4
C
                              CSUM = C0
                              DO J = 1,NKM
                                 DO I = 1,NKM
                                    CSUM = CSUM + 
     &                                 MIRR_4_VFTAB(I,J,MUE,NUE,IXX)
                                 END DO
                              END DO
C
                              S4_VFTAB(MUE,NUE,IXX,IVFT,IZA,IZB) = CSUM
C
                           END DO
                        END DO
C
                     END DO LOOP_IXX
C
                  END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
               END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
            END DO LOOP_IZB
CIZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB_IZB
C
         END DO LOOP_IZA
CIZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA_IZA
C
      END DO LOOP_IT
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( LCHECKME_SPIN ) THEN
         CALL SIG_SPIN_CHECKME(ERYDA,ERYDB,CCF,IFILA,IFILB)
         CALL STOP_REGULAR(ROUTINE,'checked MEs ---> STOP')
      END IF
C
C=======================================================================
C       determine the site dependent average for the
C       the auxilary response terms  S10, S2, S3 and S4
C       and matrix elements MAQAB, MBQAB, MCQAB, and MDQAB
C=======================================================================
C
c modified by XJQ: scfvb
      if(lscfvb .and. i_temp_lat==n_temp_lat) close(ifilmirr)
c end-mod-xjq
      CALL SIGME_AUX(MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA,
     &               MREG_TAB,MREG_TBA,S10_VFTAB,S2_VFTAB,S3_VFTAB,
     &               S4_VFTAB,S10AQAB,S10BQAB,S10CQAB,S10DQAB,S2AQAB,
     &               S2BQAB,S2CQAB,S2DQAB,S3AQAB,S3BQAB,S3CQAB,S3DQAB,
     &               S4AQAB,S4BQAB,S4CQAB,S4DQAB,MSSTA,TSSTA,MSSQA,
     &               TAUQA_TMP,MSSTB,TSSTB,MSSQB,TAUQB_TMP,TSS_FTA,
     &               TSS_FTB,TSS_VFTA,TSS_VFTB,MSS_VFTA,MSS_VFTB,
     &               UMAT_VTA,UMAT_VTB)
C
      IF ( CHECK ) CALL STOP_REGULAR(ROUTINE,'CHECK completed')
c modified by XJQ: scfvb
      if(lscfvb .and. i_temp_lat==n_temp_lat) then
        close(ifilscftv)
        close(ifilmreg)
      endif
c end-mod-xjq
C
99001 FORMAT (//,' J - matrix elements for component IT=',I3)
99002 FORMAT (/,5X,'<',A8,'>',' Info: ',
     &        'MEs calculated for complex energies')
99003 FORMAT (80('*'),/,'<SIGME> check MEs',/,80('*'))
      END
C*==sigme_debug.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGME_DEBUG(KDEBUG,IT,MREG_TAB,MREG_TBA,MIRR_2,MIRR_3,
     &                       MIRR_4,NSPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *  print matrix elements for debugging purposes ONLY               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      IMPLICIT NONE
C*--SIGME_DEBUG848
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TWR
      PARAMETER (TWR=1D-8)
C
C Dummy arguments
C
      INTEGER IT,KDEBUG,NSPINPROJ
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           MIRR_3(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           MIRR_4(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ),
     &           MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ)
C
C Local variables
C
      INTEGER FILN,I,ISP,J
C
C*** End of declarations rewritten by SPAG
C
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
C
      IF ( KDEBUG.EQ.1 ) THEN
         OPEN (UNIT=666,FILE='MEs_sigma_sp_alpha.txt')
         OPEN (UNIT=667,FILE='MEs_sigma_sp_nabla.txt')
         DO ISP = 1,7
            FILN = 666
            IF ( ISP.GT.4 ) FILN = 667
            DO J = 1,3
C
               WRITE (FILN,*) 'mreg_tba it jpol isp',IT,J,ISP
               CALL CMATSTRUCT('MREG_TBA',MREG_TBA(1,1,J,ISP),NKM,NKM,3,
     &                         3,0,TWR,FILN)
C
C                       WRITE (FILN,*)  'mreg_tab it jpol isp',IT,J,ISP
C                       CALL  CMATSTR('MREG_TAB',MREG_TAB(1,1,J,ISP),
C    &                               NKM,NKM,3,3,0,TWR,FILN)
C
               DO I = 1,3
                  WRITE (FILN,*) 'mirr_2 it ipol jpol isp',IT,I,J,ISP
                  CALL CMATSTRUCT('MIRR_2',MIRR_2(1,1,I,J,ISP),NKM,NKM,
     &                            3,3,0,TWR,FILN)
C
                  WRITE (FILN,*) 'mirr_3 it ipol jpol isp',IT,I,J,ISP
                  CALL CMATSTRUCT('MIRR_3',MIRR_3(1,1,I,J,ISP),NKM,NKM,
     &                            3,3,0,TWR,FILN)
C
                  WRITE (FILN,*) 'mirr_4 it ipol jpol isp',IT,I,J,ISP
                  CALL CMATSTRUCT('MIRR_4',MIRR_4(1,1,I,J,ISP),NKM,NKM,
     &                            3,3,0,TWR,FILN)
               END DO
            END DO
         END DO
         CLOSE (666)
         CLOSE (667)
C
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
C
      ELSE IF ( KDEBUG.EQ.2 ) THEN
C
         DO J = 1,3
            WRITE (888,*) 'mreg_tba it jpol',IT,J
            CALL CMATSTRUCT('MREG_TBA',MREG_TBA(1,1,J,4),NKM,NKM,3,3,0,
     &                      TWR,888)
         END DO
         DO J = 1,3
            WRITE (888,*) 'mreg_tab it jpol',IT,J
            CALL CMATSTRUCT('MREG_TAB',MREG_TAB(1,1,J,4),NKM,NKM,3,3,0,
     &                      TWR,888)
         END DO
C
      END IF
C
CDEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__DEBUG__D
C
      END
