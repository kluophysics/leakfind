C*==gfunmat_drive.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GFUNMAT_DRIVE(ERYD,IE,TAUT,MSST,GFMAT)
C
C   ********************************************************************
C   *                                                                  *
C   *  Driver to calculate Greens function matrix                      *
C   *                                                                  *
C   *  Set of localised basis functions is stored in IFILGFWF          *
C   *  Choise of basis functions:                                      *
C   *  In order to obtain full set of fucntions we set BT,BNST,VNST=0  *
C   *  BASIS=1: This is now default: Basis function only within MT     *
C   *                                                                  *
C   *  - if  LHS_SOL_EQ_RHS_SOL = .FALSE.                              *
C   *    ZGA1,ZFA1,JGA1,JFA1 correspond to the LHS solution            *
C   *    ZGA2,ZFA2,JGA2,JFA2 correspond to the RHS solution            *
C   *                                                                  *
C   *    otherwise LHS = RHS solution with ZGA1 = ZGA2                 *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:IFILGFWF,IFILCBWF,IPRINT,IFILCBWF_LHS
      USE MOD_TYPES,ONLY:BT,VNST,BNST,NTMAX,ITBOT,ITTOP,NT,NLMFPMAX,
     &    NCPLWFMAX,VT,CTL
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRNSMIN
      USE MOD_SITES,ONLY:IQBOT,IQTOP
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NLMAX
      USE MOD_CALCMODE,ONLY:DMFT,LDAU,LHS_SOL_EQ_RHS_SOL
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,EREFDMFT,EREFLDAU,KSELF,IBASIS
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_RMESH,ONLY:FLMSF
      IMPLICIT NONE
C*--GFUNMAT_DRIVE33
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IE
      COMPLEX*16 GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BNST0(:,:,:),BT0(:,:),CTL0(:,:),SFTMP(:,:,:),VNST0(:,:,:),
     &       VT0(:,:)
      LOGICAL CALCINT
      COMPLEX*16 EREFTMP,JFA1(:,:),JFA2(:,:),JGA1(:,:),JGA2(:,:),
     &           MEZJ0(:,:,:,:),MEZZ0(:,:,:,:),MSST0(:,:,:),P,
     &           SSST(:,:,:),TMAT(:,:,:),TSST(:,:,:),ZFA1(:,:),ZFA2(:,:)
     &           ,ZGA1(:,:),ZGA2(:,:)
      INTEGER I,IA_ERR,IQBOTSAV,IQTOPSAV,IT,ITBOTSAV,ITTOPSAV,IWRIRRWF,
     &        IWRREGWF,J,K,MC
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BT0,VNST0,BNST0,TSST,MSST0,SSST
      ALLOCATABLE VT0,CTL0
      ALLOCATABLE MEZZ0,MEZJ0,TMAT
      ALLOCATABLE SFTMP
      ALLOCATABLE ZGA1,ZFA1,ZGA2,ZFA2
      ALLOCATABLE JGA1,JFA1,JGA2,JFA2
C
C
      ALLOCATE (BT0(NRMAX,NTMAX))
      ALLOCATE (VT0(NRMAX,NTMAX))
      IF ( FULLPOT ) THEN
         ALLOCATE (VNST0(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (BNST0(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
      END IF
      ALLOCATE (CTL0(NTMAX,NLMAX))
      ALLOCATE (MEZJ0(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (MEZZ0(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (TSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (MSST0(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (SSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TMAT(NKMMAX,NKMMAX,NTMAX))
C
      IWRREGWF = 1
      IWRIRRWF = 1
      CALCINT = .TRUE.
      CALL RVECCOP(NRMAX*NT,BT,BT0)
      CALL RINIT(NRMAX*NT,BT)
      CALL RVECCOP(NRMAX*NT,VT,VT0)
      CTL0 = CTL
      CTL = CTL/((1D-6)**0.5D0)
      IF ( FULLPOT ) THEN
         VNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &      = VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
         BNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &      = BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
C
         VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0.0D0
         BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0.0D0
      END IF
C
      DMFTSIG(1:NKM,1:NKM,1:NT) = C0
C
C
C=======================================================================
C                         calculate projectors for local orbitals
C=======================================================================
C
      ITBOTSAV = ITBOT
      ITTOPSAV = ITTOP
      IQBOTSAV = IQBOT
      IQTOPSAV = IQTOP
C
      DO IT = ITBOT,ITTOP
         CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT(1,1,IT))
      END DO
C
      TSST = C0
      SSST = C0
      MSST0 = C0
      MEZZ0 = C0
      MEZJ0 = C0
C
      CALL DMFT_SSITE_ROTATE(TMAT)
C
      I = SIZE(FLMSF,1)
      J = SIZE(FLMSF,2)
      K = SIZE(FLMSF,3)
      ALLOCATE (SFTMP(I,J,K))
      SFTMP = FLMSF
      IF ( IBASIS.EQ.1 ) FLMSF(:,1:J,:) = 0.0D0
      IF ( IBASIS.EQ.2 ) FLMSF(:,2:J,:) = 0.0D0
      IF ( IBASIS.EQ.3 ) FLMSF(:,1:J,:) = 0.0D0
      IF ( IBASIS.EQ.4 ) FLMSF(:,2:J,:) = 0.0D0
C      IF (.NOT. DMFT_FIX_BAS .OR. ITRSCF.LE.2) THEN
      I = ITBOTSAV
      J = ITTOPSAV
      DO IT = I,J
         IF ( KSELF(IT).EQ.1 ) THEN
            ITBOT = IT
            ITTOP = IT
            IF ( DMFT ) THEN
               EREFTMP = EREFDMFT(IT)
            ELSE IF ( LDAU ) THEN
               EREFTMP = EREFLDAU(IT)
            END IF
c            WRITE(*,*) "XXX",IT
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILGFWF,.TRUE.,
     &                    EREFTMP,P,IPRINT,TSST,MSST0,SSST,MEZZ0,MEZJ0,
     &                    'NONE      ')
         END IF
      END DO
C
C      ELSE
C         WRITE(6,*) "DMFT LDAU BASIS FIXED"
C      END IF
      ITBOT = ITBOTSAV
      ITTOP = ITTOPSAV
      IQBOT = IQBOTSAV
      IQTOP = IQTOPSAV
      VT = VT0
      CTL = CTL0
      CALL RVECCOP(NRMAX*NT,BT0,BT)
      IF ( FULLPOT ) THEN
         VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &      = VNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
         BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &      = BNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
      END IF
C
      IF ( IBASIS.GT.2 ) FLMSF = SFTMP
C
C=======================================================================
C                           calculate and store GFMAT for present E-mesh
C=======================================================================
C
C
      MC = NCPLWFMAX
      IF ( .NOT.FULLPOT ) MC = 2
C
      ALLOCATE (JFA2(NRMAX,MC),ZFA2(NRMAX,MC))
      ALLOCATE (JGA2(NRMAX,MC))
      ALLOCATE (ZGA2(NRMAX,MC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:gfunmat->ZGA2'
      JFA2(:,:)=C0
      ZFA2(:,:)=C0
      JGA2(:,:)=C0
      ZGA2(:,:)=C0

C
      IF ( LHS_SOL_EQ_RHS_SOL ) THEN
         CALL GFUNMAT(IFILGFWF,IFILCBWF,IFILCBWF_LHS,IPRINT,ERYD,MEZZ0,
     &                TMAT,ZGA2,ZGA2,ZFA2,ZFA2,JGA2,JGA2,JFA2,JFA2,MC,
     &                IE,GFMAT)
C
      ELSE
         ALLOCATE (JFA1(NRMAX,MC),ZFA1(NRMAX,MC))
         ALLOCATE (JGA1(NRMAX,MC))
         ALLOCATE (ZGA1(NRMAX,MC),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:gfunmat->ZGA1'
C
         JFA2(:,:)=C0
         ZFA2(:,:)=C0
         JGA2(:,:)=C0
         ZGA2(:,:)=C0
         CALL GFUNMAT(IFILGFWF,IFILCBWF,IFILCBWF_LHS,IPRINT,ERYD,MEZZ0,
     &                TMAT,ZGA1,ZGA2,ZFA1,ZFA2,JGA1,JGA2,JFA1,JFA2,MC,
     &                IE,GFMAT)
      END IF
      FLMSF = SFTMP
      END
