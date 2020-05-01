C*==thermal_init_ufmat.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE THERMAL_INIT_UFMAT(ERYD)
C   ********************************************************************
C   *                                                                  *
C   * in case of thermal lattice vibrations and/or spin fluctuations   *
C   *                                                                  *
C   * - init displacement matrices                 UMAT_VT             *
C   * - init fluctuation  matrices                 FMAT_FT             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_THERMAL,ONLY:SVEC_VT,UMAT_VT,NVIBRA,NFLUCT,NVT,NFT,
     &    FTET_FT,FPHI_FT,FMAT_FT,MFMAT_FT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,MROT_T,DROT_T
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,NLM,WKM1,WLM1,NK
      USE MOD_CALCMODE,ONLY:IREL,THERMAL_VIBRA_FLUCT,MOMENTS_ROTATED
      USE MOD_CONSTANTS,ONLY:C0,C1,DEG_ARC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='THERMAL_INIT_UFMAT')
      REAL*8 ZER0_ANGLE
      PARAMETER (ZER0_ANGLE=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
C
C Local variables
C
      INTEGER I,IFLUCT,IFT,IT,IVIBRA,IVT,M
      REAL*8 MAUX(3,3),XPHI,XTET
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      IF ( .NOT.ALLOCATED(UMAT_VT) ) ALLOCATE (UMAT_VT(M,M,NVT))
      IF ( .NOT.ALLOCATED(FMAT_FT) ) THEN
         ALLOCATE (FMAT_FT(M,M,NFT))
         ALLOCATE (MFMAT_FT(3,3,NFT))
      END IF
C
C-----------------------------------------------------------------------
C                  standard mode: NO thermal effects
C-----------------------------------------------------------------------
      IF ( .NOT.THERMAL_VIBRA_FLUCT ) THEN
C
         UMAT_VT(:,:,1:NVT) = C0
         FMAT_FT(:,:,1:NFT) = C0
         MFMAT_FT(:,:,1:NFT) = 0D0
C
         DO I = 1,M
            UMAT_VT(I,I,1:NVT) = C1
            FMAT_FT(I,I,1:NFT) = C1
         END DO
         DO I = 1,3
            MFMAT_FT(I,I,1:NFT) = 1D0
         END DO
C
         RETURN
C
      END IF
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C     initialize type-dependent rotation matrices MROT_T and DROT_T
C          if the average moment is rotated from the global z-axis
C-----------------------------------------------------------------------
C
      CALL INIT_MOD_TYPES_ROTMAG
C
      DO IT = ITBOT,ITTOP
C
C ----------------------------------------------------------------------
C    initialize transformation matrix UMAT for displacements
C-----------------------------------------------------------------------
C
         IF ( NVIBRA.EQ.1 ) THEN
C
            UMAT_VT(:,:,1:NVT) = C0
C
            DO I = 1,M
               UMAT_VT(I,I,1:NVT) = C1
            END DO
C
         ELSE
C
            IVT = (IT-1)*NVIBRA
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               IF ( IREL.LE.2 ) THEN
C
                  CALL INIT_UMAT(SVEC_VT(1,IVT),WLM1,ERYD)
C
                  UMAT_VT(1:NLM,1:NLM,IVT) = WLM1(1:NLM,1:NLM)
C
               ELSE
C
                  CALL INIT_UMAT(SVEC_VT(1,IVT),WLM1,ERYD)
C
                  WKM1(1:NKMMAX,1:NKMMAX) = C0
C
                  WKM1(1:NLM,1:NLM) = WLM1(1:NLM,1:NLM)
                  WKM1(NLM+1:NLM+NLM,NLM+1:NLM+NLM) = WLM1(1:NLM,1:NLM)
C
                  CALL CHANGEREP(NKM,NKMMAX,WKM1(1,1),'RLM>REL',
     &                           UMAT_VT(1,1,IVT))
C
               END IF
C
            END DO
C
         END IF
C
C ----------------------------------------------------------------------
C    initialize transformation matrix FMAT for fluctuations
C
C         NOTE: setting of FMAT_FT and MFMAT_FT for (TET,PHI)
C           as for DROTQ and MROTQ in <INIT_MOD_SITES_ROTMAG>
C-----------------------------------------------------------------------
C
         IFT = (IT-1)*NFLUCT
         DO IFLUCT = 1,NFLUCT
            IFT = IFT + 1
C
            IF ( NFLUCT.EQ.1 .OR. 
     &           (ABS(FTET_FT(IFT)).LT.ZER0_ANGLE .AND. ABS(FPHI_FT(IFT)
     &           ).LT.ZER0_ANGLE) ) THEN
C
               FMAT_FT(:,:,IFT) = C0
               MFMAT_FT(:,:,IFT) = 0D0
C
               DO I = 1,M
                  FMAT_FT(I,I,IFT) = C1
               END DO
               DO I = 1,3
                  MFMAT_FT(I,I,IFT) = 1D0
               END DO
C
            ELSE IF ( IREL.LE.2 ) THEN
C
               CALL STOP_MESSAGE(ROUTINE,'IREL <= 2')
C
            ELSE
C
               CALL ROTMAT(NK,3,FPHI_FT(IFT),FTET_FT(IFT),0D0,
     &                     FMAT_FT(1,1,IFT),NKMMAX)
C
               CALL GETMROT(0.0D0,-FTET_FT(IFT),-FPHI_FT(IFT),
     &                      MFMAT_FT(1,1,IFT))
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C               check NEW set up of rotation matrices
C               this section can later on be removed
C               use extended version based on GETMROT instead
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
               MAUX(:,:) = MFMAT_FT(:,:,IFT)
C
               XPHI = FPHI_FT(IFT)*DEG_ARC
               XTET = FTET_FT(IFT)*DEG_ARC
C
               MFMAT_FT(1,1,IFT) = COS(XPHI)*COS(XTET)
               MFMAT_FT(1,2,IFT) = SIN(XPHI)*COS(XTET)
               MFMAT_FT(1,3,IFT) = -SIN(XTET)
               MFMAT_FT(2,1,IFT) = -SIN(XPHI)
               MFMAT_FT(2,2,IFT) = COS(XPHI)
               MFMAT_FT(2,3,IFT) = 0D0
               MFMAT_FT(3,1,IFT) = COS(XPHI)*SIN(XTET)
               MFMAT_FT(3,2,IFT) = SIN(XPHI)*SIN(XTET)
               MFMAT_FT(3,3,IFT) = COS(XTET)
C
               IF ( .NOT.RVEC_SAME(9,MAUX,MFMAT_FT(1,1,IFT),1D-8) ) THEN
                  WRITE (6,*) '  for IFT = ',IFT
                  WRITE (6,'(''MAUX  = '',/,(3F12.8))') MAUX
                  WRITE (6,'(''MFMAT = '',/,(3F12.8))')
     &                   MFMAT_FT(:,:,IFT)
                  CALL STOP_MESSAGE(ROUTINE,'MROT matrices DIFFER')
               END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C ----------------------------------------------------------------------
C    multiply with type-dependent rotation matrices MROT_T and DROT_T
C          if the average moment is rotated from the global z-axis
C-----------------------------------------------------------------------
               IF ( MOMENTS_ROTATED ) THEN
C
                  WKM1(:,:) = FMAT_FT(:,:,IFT)
                  FMAT_FT(:,:,IFT) = MATMUL(DROT_T(:,:,IT),WKM1(:,:))
C
                  MAUX(:,:) = MFMAT_FT(:,:,IFT)
                  MFMAT_FT(:,:,IFT) = MATMUL(MAUX(:,:),MROT_T(:,:,IT))
C
               END IF
C
            END IF
C
         END DO
C
      END DO
C
      END
