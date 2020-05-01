C*==fp_pot_local.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FP_POT_LOCAL
C   ********************************************************************
C   *                                                                  *
C   *   the non-spherical potentials refer to the LOCAL frames         *
C   *   aligned to the magnetic moments. Rotate the potentials         *
C   *   if the orientation set in the input differs from that          *
C   *   specified in the potential file                                *
C   *                                                                  *
C   *   PRESENT STATUS: rotation only accepted if angles are 0         *
C   *                   in the potential file - otherwise  STOP        *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:DMFT,ORBPOL
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SITES,ONLY:QMPHI,QMTET,QMGAM,IQAT,QMTET_POTFIL,
     &    QMPHI_POTFIL,QMGAM_POTFIL
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_RMESH,ONLY:NRMAX,JRNSMIN
      USE MOD_TYPES,ONLY:NFPT,ITBOT,ITTOP,KLMFP,LMIFP,VNST,BNST,RHO2NS,
     &    NLFP,NLMFP,NLMFPMAX
      USE MOD_DMFT_LDAU,ONLY:SEVNST,SEBNST
      IMPLICIT NONE
C*--FP_POT_LOCAL24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Local variables
C
      REAL*8 BNST_OLD(:,:,:),DROT_FP(:,:),D_21,RHO2NS_OLD(:,:,:,:),
     &       VNST_OLD(:,:,:)
      INTEGER IA_ERR,IFP,IQ,IT,KFP_LMT_OLD(:,:),L2,LM1,LM2,
     &        LM_FPT_OLD(:,:),M2,N,NFPTMAX_NEW,NFPTMAX_OLD,NFPT_OLD(:)
      LOGICAL KLMFP_NE_0,MROTATION_DIFFER,MROTATION_POTFIL
      COMPLEX*16 SEBNST_OLD(:,:,:),SEVNST_OLD(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE VNST_OLD,BNST_OLD,RHO2NS_OLD,SEVNST_OLD,SEBNST_OLD
      ALLOCATABLE KFP_LMT_OLD,NFPT_OLD,LM_FPT_OLD,DROT_FP
C
      MROTATION_POTFIL = .FALSE.
      MROTATION_DIFFER = .FALSE.
C
C=======================================================================
C                         check angle settings
C=======================================================================
C
      DO IT = ITBOT,ITTOP
         IQ = IQAT(1,IT)
C
         IF ( ABS(QMTET_POTFIL(IQ)).GT.TOL .OR. ABS(QMPHI_POTFIL(IQ))
     &        .GT.TOL .OR. ABS(QMGAM_POTFIL(IQ)).GT.TOL )
     &        MROTATION_POTFIL = .TRUE.
C
         IF ( ABS(QMTET_POTFIL(IQ)-QMTET(IQ)).GT.TOL .OR. 
     &        ABS(QMPHI_POTFIL(IQ)-QMPHI(IQ)).GT.TOL .OR. 
     &        ABS(QMGAM_POTFIL(IQ)-QMGAM(IQ)).GT.TOL )
     &        MROTATION_DIFFER = .TRUE.
C
      END DO
C
      IF ( .NOT.MROTATION_DIFFER ) RETURN
C
      WRITE (6,99002)
      DO IT = ITBOT,ITTOP
         IQ = IQAT(1,IT)
         WRITE (6,99003) IQ,QMPHI_POTFIL(IQ),QMTET_POTFIL(IQ),
     &                   QMGAM_POTFIL(IQ),QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
      END DO
C
      IF ( MROTATION_POTFIL ) STOP 
     &          'in <FP_POT_LOCAL>: non-0 angles in potfile not allowed'
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      IF ( NLMFPMAX.NE.NLMFP ) STOP '<FP_POT_LOCAL>: NLMFPMAX <> NLMFP'
C
C ----------------------------------------------------------------------
C
      ALLOCATE (DROT_FP(NLMFP,NLMFP),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: FP_POT_LOCAL -> DROT_FP'
C
C-----------------------------------------------------------------------
C    swap the full potential information on temporary arrays  *_OLD
C-----------------------------------------------------------------------
C
      ALLOCATE (VNST_OLD(JRNSMIN:NRMAX,NLMFPMAX,ITBOT:ITTOP))
      ALLOCATE (BNST_OLD(JRNSMIN:NRMAX,NLMFPMAX,ITBOT:ITTOP))
      ALLOCATE (RHO2NS_OLD(NRMAX,NLMFPMAX,ITBOT:ITTOP,3))
      ALLOCATE (NFPT_OLD(ITBOT:ITTOP))
      ALLOCATE (KFP_LMT_OLD(NLMFPMAX,ITBOT:ITTOP))
      ALLOCATE (LM_FPT_OLD(NLMFPMAX,ITBOT:ITTOP))
C
      VNST_OLD(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
     &   = VNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
      BNST_OLD(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
     &   = BNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
      RHO2NS_OLD(1:NRMAX,1:NLMFPMAX,ITBOT:ITTOP,1:3)
     &   = RHO2NS(1:NRMAX,1:NLMFPMAX,ITBOT:ITTOP,1:3)
C
      NFPT_OLD(ITBOT:ITTOP) = NFPT(ITBOT:ITTOP)
      KFP_LMT_OLD(1:NLMFPMAX,ITBOT:ITTOP)
     &   = KLMFP(1:NLMFPMAX,ITBOT:ITTOP)
      LM_FPT_OLD(1:NLMFPMAX,ITBOT:ITTOP) = LMIFP(1:NLMFPMAX,ITBOT:ITTOP)
C
      VNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP) = 0D0
      BNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP) = 0D0
      RHO2NS(1:NRMAX,1:NLMFPMAX,ITBOT:ITTOP,1:3) = 0D0
C
      NFPT(ITBOT:ITTOP) = 0
      KLMFP(1:NLMFPMAX,ITBOT:ITTOP) = 0
      LMIFP(1:NLMFPMAX,ITBOT:ITTOP) = 0
C
      IF ( ORBPOL(1:6).EQ.'SIGMA ' .OR. DMFT ) THEN
         ALLOCATE (SEVNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP))
         ALLOCATE (SEBNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP))
         ALLOCATE (SEVNST_OLD(JRNSMIN:NRMAX,NLMFPMAX,ITBOT:ITTOP))
         ALLOCATE (SEBNST_OLD(JRNSMIN:NRMAX,NLMFPMAX,ITBOT:ITTOP))
         SEVNST_OLD(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
     &      = SEVNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
         SEBNST_OLD(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
     &      = SEBNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP)
         SEVNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP) = 0D0
         SEBNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITBOT:ITTOP) = 0D0
      END IF
C
C-----------------------------------------------------------------------
      NFPTMAX_OLD = 0
      DO IT = ITBOT,ITTOP
         N = 0
         LM2 = 0
         DO L2 = 0,(NLFP-1)
            DO M2 = -L2,L2
               LM2 = LM2 + 1
               IF ( KFP_LMT_OLD(LM2,IT).EQ.1 ) N = N + 1
            END DO
         END DO
         NFPTMAX_OLD = MAX(NFPTMAX_OLD,N)
      END DO
C
C=======================================================================
C                     rotate the potentials
C=======================================================================
C
      NFPTMAX_NEW = 0
C
      DO IT = ITBOT,ITTOP
         IQ = IQAT(1,IT)
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for REAL spherical harmonics
C ----------------------------------------------------------------------
C
         CALL ROTMAT_RYLM(NLFP,NLMFP,QMPHI(IQ),QMTET(IQ),QMGAM(IQ),
     &                    DROT_FP)
C
         IF ( IPRINT.GE.1 ) THEN
            WRITE (6,99001) IT,IQ,QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
            CALL RMATSTRUCT('DROT_FP',DROT_FP,NLMFP,NLMFP,1,1,0,1D-8,6)
         END IF
C
C-----------------------------------------------------------------------
C
         DO LM1 = 1,NLMFP
C
            KLMFP_NE_0 = .FALSE.
C
            DO LM2 = 1,NLMFP
C
               D_21 = DROT_FP(LM2,LM1)
               IF ( ABS(D_21).GT.1D-8 .AND. KFP_LMT_OLD(LM2,IT).EQ.1 )
     &              THEN
C
                  VNST(JRNSMIN:NRMAX,LM1,IT)
     &               = VNST(JRNSMIN:NRMAX,LM1,IT)
     &               + VNST_OLD(JRNSMIN:NRMAX,LM2,IT)*D_21
                  BNST(JRNSMIN:NRMAX,LM1,IT)
     &               = BNST(JRNSMIN:NRMAX,LM1,IT)
     &               + BNST_OLD(JRNSMIN:NRMAX,LM2,IT)*D_21
                  RHO2NS(1:NRMAX,LM1,IT,1:3)
     &               = RHO2NS(1:NRMAX,LM1,IT,1:3)
     &               + RHO2NS_OLD(1:NRMAX,LM2,IT,1:3)*D_21
C
                  IF ( ORBPOL(1:6).EQ.'SIGMA ' .OR. DMFT ) THEN
                     SEVNST(JRNSMIN:NRMAX,LM1,IT)
     &                  = SEVNST(JRNSMIN:NRMAX,LM1,IT)
     &                  + SEVNST_OLD(JRNSMIN:NRMAX,LM2,IT)*D_21
                     SEBNST(JRNSMIN:NRMAX,LM1,IT)
     &                  = SEBNST(JRNSMIN:NRMAX,LM1,IT)
     &                  + SEBNST_OLD(JRNSMIN:NRMAX,LM2,IT)*D_21
                  END IF
C
                  KLMFP_NE_0 = .TRUE.
               END IF
C
            END DO
C
            IF ( KLMFP_NE_0 ) THEN
               NFPT(IT) = NFPT(IT) + 1
               IFP = NFPT(IT)
               LMIFP(IFP,IT) = LM1
               KLMFP(LM1,IT) = 1
            END IF
C
         END DO
C
         NFPTMAX_NEW = MAX(NFPTMAX_NEW,NFPT(IT))
C
C-----------------------------------------------------------------------
C
      END DO
C
C=======================================================================
C
      WRITE (6,99004) NFPTMAX_NEW,NFPTMAX_OLD
      DO IT = ITBOT,ITTOP
         IQ = IQAT(1,IT)
C
         WRITE (6,99005) 'type',IT
C
         IF ( IPRINT.GE.-1 ) THEN
            WRITE (6,99008) 'old settings '
            WRITE (6,99006) 'NFP ',NFPT_OLD(IT),QMPHI_POTFIL(IQ),
     &                      QMTET_POTFIL(IQ),QMGAM_POTFIL(IQ)
            WRITE (6,99007) (LM_FPT_OLD(IFP,IT),L_LM(LM_FPT_OLD(IFP,IT))
     &                      ,M_LM(LM_FPT_OLD(IFP,IT)),IFP=1,NFPT_OLD(IT)
     &                      )
            WRITE (6,99008) 'new settings '
         END IF
C
         WRITE (6,99006) 'NFP ',NFPT(IT),QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
         WRITE (6,99007) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                   M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
      END DO
C
      WRITE (6,99005)
C
C-----------------------------------------------------------------------
99001 FORMAT (5X,'rotation matrix for IT =',I3,2X,'IQ =',I3,2X,'ALF =',
     &        F6.1,2X,'BET =',F6.1,2X,'GAM =',F6.1,/)
99002 FORMAT (/,1X,79('*'),/,33X,'<FP_POT_LOCAL>',/,1X,79('*'),//,10X,
     &        'angles of canted magnetic moments in potential file',/,
     &        10X,'differ from angle settings in input file',/,10X,
     &        'the non-spherical potentials will be transferred ',/,10X,
     &        'to the new local frames aligned to magnetic moments',//,
     &        10X,'   IQ         TET       PHI       GAM (pot)',11X,
     &        'TET       PHI       GAM (inp)')
99003 FORMAT (I15,3X,3F10.3,10X,3F10.3)
99004 FORMAT (/,10X,'rotated non-spherical potentials generated',//,10X,
     &        'NFPTMAX =',I6,'   new number of potential terms ',/,I25,
     &        '   setting prior to rotation',/)
99005 FORMAT (/,:,10X,A,I4,/)
99006 FORMAT (/,:,10X,A,I4,5X,'PHI =',F7.2,5X,'TET =',F7.2,5X,'GAM =',
     &        F7.2,/)
99007 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99008 FORMAT (/,10X,A,/)
      END
