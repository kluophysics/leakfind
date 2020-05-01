C*==fpnrssite.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPNRSSITE(IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,ERYD,P,
     &                     IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOLSS)
C   ********************************************************************
C   *                                                                  *
C   *    call <FPRWFBS> to solve the dirac equation without            *
C   *    spin-orbit coupling and set up the radial integrals           *
C   *                                                                  *
C   *  IREL = 0 para-magnetic  non-relativistic     }                  *
C   *                                               } (l,m_l)-repres.  *
C   *  IREL = 1 para-magnetic  scalar-relativistic  }                  *
C   *                                                                  *
C   *  IREL = 2 spin-polarized  scalar-relativistic                    *
C   *           in (l,m_l,m_s)-representation                          *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   *  NOTE: the minor component F is seen as an auxilary function     *
C   *        without meaning - it enters only the Wronskian            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NLMAX,NCPLWF,NKM,NL,L_LM,NLM,
     &    NSPIN,WKM2,IPIVKM
      USE MOD_RMESH,ONLY:FINITE_NUCLEUS,JRNSMIN,NMMAX,NRMAX,NPAN,JRNS1,
     &    JRCUT,JRCRI,DRDI,R,NM,Y_LEBGRID,W_LEBGRID,D_LGM,NRSFMAX
      USE MOD_TYPES,ONLY:CTL,NTMAX,NCPLWFMAX,NLMFPMAX,NLMFP,VAMEG,
     &    IKMSOLBLK,LMIFP,NBLK,NFPT,NSOLBLK,IKMCPLWF,NLMFPT,KLMFP,LOPT,
     &    BNST,VNST,Z,BT,VT,IMT,ITBOT,ITTOP,NT,NLT
      USE MOD_CALCMODE,ONLY:IREL,SOLVER_FP,GF_CONV_RH,BREAKPOINT,
     &    L_PROJECTED_ME,FP_USE_DIRECTIONS
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,SEVT,SEBT,SEBNST,SEVNST
      USE MOD_FILES,ONLY:FILNAM,LFILNAM,IOTMP,DEBUG,IFILCBWF_SPH,
     &    IFILGFWF,IFILMEZZL
      USE MOD_CONSTANTS,ONLY:C0,CI,SQRT_4PI
      IMPLICIT NONE
C*--FPNRSSITE38
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPNRSSITE')
      LOGICAL WRONSKI
      PARAMETER (WRONSKI=.FALSE.)
      REAL*8 TOL
      PARAMETER (TOL=1D-6)
      INTEGER N_LEBGRID
      PARAMETER (N_LEBGRID=434)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      LOGICAL GETIRRSOL
      INTEGER IFILSS,IPRINT,IWRIRRWF,IWRREGWF
      CHARACTER*10 ORBPOLSS
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BAUX_LGD(:),BAUX_LM(:),BISLM(:,:),BIS_LGD(:,:),C,CSQR,
     &       GAMMA(:,:),IGAMMA(:,:),RCRI_LAST(:),RGAM(:,:,:),
     &       RIGAM(:,:,:),SUMW,VAUX_LGD(:),VAUX_LM(:),VISLM(:,:),
     &       VIS_LGD(:,:)
      LOGICAL CALC_GFWF,DUMP_SSITE,GET_ALL_SPHER_SOL,INITIALIZE,LDAU,
     &        MODIFYV,RMESH_CHANGED
      COMPLEX*16 CBLK(:,:),CBNST(:,:),CBT(:),CFAC,CHIS(:,:,:),CHL,
     &           CHLM(:),CHLP(:),CJL,CJLP1,CJL_J,CNL,CNLP1,CRSQ,CSUM,
     &           CV0(:),CVNST(:,:),CVT(:),DCHLM(:),DCHLP(:),DGR(:,:),
     &           DPRIRTOP(:,:),GAMS(:,:,:),GBLK(:,:),GR(:,:),JF(:,:,:),
     &           JG(:,:,:),MEZJL(:,:,:,:),MEZZL(:,:,:,:),PH0LS(:,:,:),
     &           PIM(:,:,:),PJL(:,:),PR0LS(:,:,:),PRM(:,:,:),QHL(:,:),
     &           QIM(:,:,:),QRL(:,:),QRM(:,:,:),RSCLP,RSCLQ,RTOP,
     &           SBLK(:,:),SMT0L(:),SPNWGT,TBLK(:,:),TMT0L(:),W(:,:),
     &           WFAC,XBLK(:,:),YBLK(:,:),YBLKINV(:,:),ZBOT,ZF(:,:,:),
     &           ZG(:,:,:),ZTOP
      COMPLEX*16 CJLZ,CNLZ
      REAL*8 DDOT
      INTEGER I,I1,I2,IA_ERR,IBLK,IL,ILM,ILMBOT,ILMS,ILMTOP,IM,IMP,INFO,
     &        IPIV(:),IPOT,IR,IRCRIT,IRMTIN,IRP,IRSF,IRTOP,IS,ISOL,IT,
     &        ITP,I_LEBGRID,J,JLM,JLMS,K,L,LI,LJ,LM,LMP,LMSOFF,
     &        NRSPHERSOL,NSOL,NT_LAST_CALL
      CHARACTER*10 SOLVER_LAST_CALL
      SAVE GAMMA,IGAMMA,RCRI_LAST,RGAM,RIGAM
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
      DATA NT_LAST_CALL/0/,SOLVER_LAST_CALL/'**********'/
      DATA DUMP_SSITE/.FALSE./
C
      ALLOCATABLE QRL,QHL,PR0LS,PH0LS,CBLK,GBLK,CHIS,GAMS
      ALLOCATABLE JF,JG,GR,ZF,ZG,CBT,DGR,SBLK,XBLK,YBLK,PIM,QIM,CVT,CV0
      ALLOCATABLE PRM,QRM,TBLK,CHLM,CHLP,YBLKINV,IPIV,DCHLM,DCHLP,CBNST
      ALLOCATABLE CVNST,DPRIRTOP,TMT0L,SMT0L
      ALLOCATABLE RIGAM,RGAM,GAMMA,IGAMMA
      ALLOCATABLE MEZJL,MEZZL,W,RCRI_LAST
      ALLOCATABLE PJL
      ALLOCATABLE VISLM,BISLM,VIS_LGD,BIS_LGD
      ALLOCATABLE VAUX_LGD,BAUX_LGD,VAUX_LM,BAUX_LM
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX),JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX),ZG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (GR(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (CBT(NRMAX),DGR(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (GBLK(NCPLWFMAX,NCPLWFMAX),GAMS(NLM,NLM,NSPIN))
      ALLOCATE (CBLK(NCPLWFMAX,NCPLWFMAX),CHIS(NLM,NLM,NSPIN))
      ALLOCATE (XBLK(NCPLWFMAX,NCPLWFMAX),YBLK(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (PIM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      ALLOCATE (PRM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      ALLOCATE (CVT(NRMAX),CV0(NRMAX))
      ALLOCATE (SBLK(NCPLWFMAX,NCPLWFMAX),SMT0L(NLMAX))
      ALLOCATE (TBLK(NCPLWFMAX,NCPLWFMAX),TMT0L(NLMAX),CHLM(0:NLMAX))
      ALLOCATE (CHLP(0:NLMAX),YBLKINV(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (IPIV(NCPLWFMAX),DCHLM(0:NLMAX),DCHLP(0:NLMAX))
      ALLOCATE (CBNST(JRNSMIN:NRMAX,NLMFPMAX))
      ALLOCATE (CVNST(JRNSMIN:NRMAX,NLMFPMAX))
      ALLOCATE (DPRIRTOP(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (W(NCPLWFMAX,NCPLWFMAX))
      IF ( SOLVER_FP.EQ.'BS' ) THEN
         ALLOCATE (QRM(NCPLWFMAX,NCPLWFMAX,NRMAX))
         ALLOCATE (QIM(NCPLWFMAX,NCPLWFMAX,NRMAX))
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) THEN
         IPRINT = 5
         DUMP_SSITE = .TRUE.
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C
C-----------------------------------------------------------------------
      IF ( IFILSS.EQ.IFILGFWF ) THEN
         CALC_GFWF = .TRUE.
      ELSE
         CALC_GFWF = .FALSE.
      END IF
C-----------------------------------------------------------------------
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      IF ( FP_USE_DIRECTIONS ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'inter_stitial_pot.dat')
         ALLOCATE (VISLM(NRSFMAX,NLMFPMAX),BISLM(NRSFMAX,NLMFPMAX))
         ALLOCATE (VIS_LGD(NRSFMAX,N_LEBGRID),VAUX_LGD(N_LEBGRID))
         ALLOCATE (BIS_LGD(NRSFMAX,N_LEBGRID),BAUX_LGD(N_LEBGRID))
         ALLOCATE (VAUX_LM(NLMFPMAX))
         ALLOCATE (BAUX_LM(NLMFPMAX))
      END IF
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      IF ( DEBUG ) WRITE (6,99003) IWRREGWF,IWRIRRWF,IFILSS,GETIRRSOL,
     &                             ORBPOLSS,ERYD
C
C--------------------------------------------------------- dummy setting
      CVT(1) = 999999D0
      CBT(1) = 999999D0
C
      IF ( SOLVER_FP.EQ.'BS' ) THEN
         GET_ALL_SPHER_SOL = .FALSE.
      ELSE
         SOLVER_FP = 'BORN'
         GET_ALL_SPHER_SOL = .TRUE.
      END IF
C
      ALLOCATE (QRL(NRMAX,NL),QHL(NRMAX,NL))
      ALLOCATE (PR0LS(NRMAX,NL,NSPIN))
      ALLOCATE (PH0LS(NRMAX,NL,NSPIN))
      ALLOCATE (PJL(NRMAX,NL))
C
      IF ( WRONSKI ) IPRINT = 5
C
      C = CTL(1,1)
      CSQR = C*C
C
C=======================================================================
      IF ( INITIALIZE ) THEN
         ALLOCATE (RCRI_LAST(NMMAX))
         RCRI_LAST(1:NMMAX) = 0D0
         INITIALIZE = .FALSE.
      END IF
C
      RMESH_CHANGED = .FALSE.
      DO IM = 1,NM
         IF ( ABS(RCRI_LAST(IM)-R(JRCRI(IM),IM)).GT.1D-12 )
     &        RMESH_CHANGED = .TRUE.
      END DO
C
      IF ( NT_LAST_CALL.LT.NT .OR. SOLVER_LAST_CALL.NE.SOLVER_FP .OR. 
     &     RMESH_CHANGED ) THEN
C
         DO IM = 1,NM
            RCRI_LAST(IM) = R(JRCRI(IM),IM)
         END DO
C
         NT_LAST_CALL = NT
         SOLVER_LAST_CALL = SOLVER_FP
C
         IF ( ALLOCATED(GAMMA) ) DEALLOCATE (GAMMA,IGAMMA,RGAM,RIGAM)
C
C-----------------------------------------------------------------------
C       exponent GAMMA  and  radial weight functions  RGAM
C-----------------------------------------------------------------------
C
         ALLOCATE (RIGAM(1:NRMAX,0:NLMAX,1:NTMAX))
         ALLOCATE (RGAM(1:NRMAX,0:NLMAX,1:NTMAX))
         ALLOCATE (GAMMA(0:NLMAX,1:NTMAX))
         ALLOCATE (IGAMMA(0:NLMAX,1:NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IGAMMA')
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IRTOP = JRCRI(IM)
            DO IL = 1,NLMAX
               L = IL - 1
               IF ( IREL.EQ.0 ) THEN
                  GAMMA(L,IT) = L + 1.0D0
                  IGAMMA(L,IT) = -L
               ELSE IF ( Z(IT).EQ.0 ) THEN
                  GAMMA(L,IT) = L + 1.0D0
                  IGAMMA(L,IT) = -L
               ELSE
                  GAMMA(L,IT) = SQRT(DBLE(L*L+L+1)-4.0D0*Z(IT)**2/C**2)
                  IGAMMA(L,IT) = -GAMMA(L,IT)
               END IF
               DO IR = 1,IRTOP
                  RGAM(IR,L,IT) = R(IR,IM)**GAMMA(L,IT)
                  RIGAM(IR,L,IT) = R(IR,IM)**IGAMMA(L,IT)
               END DO
            END DO
         END DO
C
      END IF
C=======================================================================
C           prepare calculation of l-projected matrix elements
C=======================================================================
C
      IF ( L_PROJECTED_ME ) THEN
         ALLOCATE (MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
         ALLOCATE (MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
      ELSE
         ALLOCATE (MEZJL(1,1,1,1),MEZZL(1,1,1,1))
      END IF
C
      IF ( L_PROJECTED_ME .AND. .NOT.CALC_GFWF ) REWIND IFILMEZZL
C
C=======================================================================
C
C --------------------------------------------------- calculate momentum
C
      CALL GET_MOMENTUM(IREL,C,ERYD,P)
C
      IF ( IREL.EQ.0 ) THEN
         WFAC = C
      ELSE
         WFAC = C*(1.0D0+ERYD/CSQR)
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         MEZZ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
         MEZJ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
C
         IM = IMT(IT)
         IRTOP = JRCRI(IM)
         RTOP = R(IRTOP,IM)
C
         IF ( FINITE_NUCLEUS ) CALL STOP_MESSAGE(ROUTINE,'RNUC NOT SET')
C
         CALL CINIT(NRMAX*NCPLWFMAX*NKMMAX,ZG)
         CALL CINIT(NRMAX*NCPLWFMAX*NKMMAX,ZF)
C
         DO J = 1,NKMMAX
            DO I = 1,NKMMAX
               TSST(I,J,IT) = C0
               MSST(I,J,IT) = C0
               SSST(I,J,IT) = C0
            END DO
C----------- fill diagonal with 1 to allow matrix inversion for NLT < NL
            TSST(J,J,IT) = 1D0
            MSST(J,J,IT) = 1D0
            SSST(J,J,IT) = 1D0
         END DO
         CHIS(:,:,:) = C0
         GAMS(:,:,:) = C0
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
         IF ( FP_USE_DIRECTIONS ) THEN
            IRMTIN = JRCUT(1,IM)
            IRCRIT = JRCRI(IM)
C
C-----------------------------------------------------------------------
C  read ALL potential terms V_lm  without shape functions for IS region
C-----------------------------------------------------------------------
            WRITE (IOTMP,99001) ITP,IMP
            DO LM = 1,NLMFP
               DO IR = IRMTIN + 1,IRCRIT
                  IRSF = IR - IRMTIN
                  READ (IOTMP,99001) ITP,IMP,LMP,IRP,VISLM(IRSF,LM),
     &                               BISLM(IRSF,LM)
               END DO
               IF ( IT.NE.ITP .OR. IM.NE.IMP .OR. LM.NE.LMP ) THEN
                  WRITE (*,*) '***** <FPNRSSITE>  IT,IM,LM,IR:',IT,IM,
     &                        LM,IR,'   ITP,IMP,LMP,IRP:',ITP,IMP,LMP,
     &                        IRP
                  CALL STOP_MESSAGE(ROUTINE,'reading IS potential')
               END IF
            END DO
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
            DO IR = IRMTIN + 1,IRCRIT
               IRSF = IR - IRMTIN
C
C-----------------------------------------------------------------------
C  construct V(->r) for radial mesh point IR and directions I_LEBGRID
C    set V(->r) = 0 for ->r outside the atomic cell i.e. r > D_LGM
C-----------------------------------------------------------------------
               VAUX_LGD(:) = 0D0
               BAUX_LGD(:) = 0D0
               DO LM = 1,NLMFP
                  CALL DAXPY(N_LEBGRID,VISLM(IRSF,LM),Y_LEBGRID(1,LM),1,
     &                       VAUX_LGD,1)
                  CALL DAXPY(N_LEBGRID,BISLM(IRSF,LM),Y_LEBGRID(1,LM),1,
     &                       BAUX_LGD,1)
               END DO
C
               DO I_LEBGRID = 1,N_LEBGRID
                  IF ( R(IR,IM).LT.D_LGM(I_LEBGRID,IM)*ALAT ) THEN
                     VIS_LGD(IRSF,I_LEBGRID) = VAUX_LGD(I_LEBGRID)
                     BIS_LGD(IRSF,I_LEBGRID) = BAUX_LGD(I_LEBGRID)
                  ELSE
                     VIS_LGD(IRSF,I_LEBGRID) = 0D0
                     BIS_LGD(IRSF,I_LEBGRID) = 0D0
                  END IF
               END DO
C
C-----------------------------------------------------------------------
C  construct V(->r) for radial mesh point IR and directions I_LEBGRID
C  using the interstitial potential   VNST  including shape functions
C-----------------------------------------------------------------------
               VAUX_LGD(:) = VT(IR,IT)
               BAUX_LGD(:) = BT(IR,IT)
               DO LM = 2,NLMFP
                  CALL DAXPY(N_LEBGRID,VNST(IR,LM,IT),Y_LEBGRID(1,LM),1,
     &                       VAUX_LGD,1)
                  CALL DAXPY(N_LEBGRID,BNST(IR,LM,IT),Y_LEBGRID(1,LM),1,
     &                       BAUX_LGD,1)
               END DO
C
C-----------------------------------------------------------------------
C  write cut and convoluted potential  V(->r) for each direction
C-----------------------------------------------------------------------
               DO I_LEBGRID = 1,N_LEBGRID
                  WRITE (1000+I_LEBGRID,'(5e20.12)') R(IR,IM),
     &                   VIS_LGD(IRSF,I_LEBGRID),VAUX_LGD(I_LEBGRID),
     &                   BIS_LGD(IRSF,I_LEBGRID),BAUX_LGD(I_LEBGRID)
               END DO
C
C-----------------------------------------------------------------------
C  write cut and convoluted potential  V(->r) for each direction
C-----------------------------------------------------------------------
               DO I_LEBGRID = 1,N_LEBGRID
                  VAUX_LGD(I_LEBGRID) = VIS_LGD(IRSF,I_LEBGRID)
                  BAUX_LGD(I_LEBGRID) = BIS_LGD(IRSF,I_LEBGRID)
               END DO
C
C------------------------------------ construct cut potential term  V_lm
               DO LM = 1,NLMFP
                  VAUX_LM(LM) = DDOT(N_LEBGRID,VAUX_LGD,1,W_LEBGRID(1,LM
     &                          ),1)
                  BAUX_LM(LM) = DDOT(N_LEBGRID,BAUX_LGD,1,W_LEBGRID(1,LM
     &                          ),1)
               END DO
C
C---------------------------- add spherical term to convoluted potential
               KLMFP(1,IT) = 1
               VNST(IR,1,IT) = VT(IR,IT)*SQRT_4PI
               BNST(IR,1,IT) = BT(IR,IT)*SQRT_4PI
C
C---------------------------------------------- write potentials to file
               DO LM = 1,NLMFPT(IT)
C                  IF ( KLMFP(LM,IT).NE.0 )
                  WRITE (2000+LM,'(5e20.12)') R(IR,IM),VAUX_LM(LM),
     &                   VNST(IR,LM,IT),BAUX_LM(LM),BNST(IR,LM,IT)
               END DO
C
            END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
            CALL STOP_MESSAGE(ROUTINE,'DDDDDDDDDDDDDDirections')
         END IF
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         DO IS = 1,NSPIN
C
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
            LMSOFF = NLM*(IS-1)
C
            ZTOP = P*RTOP
C
            DO L = 0,NLT(IT)
               CJL = CJLZ(L,ZTOP)
               CNL = CNLZ(L,ZTOP)
               CHLP(L) = CJL + CI*CNL
               CHLM(L) = CJL - CI*CNL
C
               CJLP1 = CJLZ(L+1,ZTOP)
               CNLP1 = CNLZ(L+1,ZTOP)
               DCHLP(L) = (DBLE(L)/RTOP)*CHLP(L) - P*(CJLP1+CI*CNLP1)
               DCHLM(L) = (DBLE(L)/RTOP)*CHLM(L) - P*(CJLP1-CI*CNLP1)
            END DO
C
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            DO IBLK = 1,NBLK(IT)
C
               NSOL = NSOLBLK(IBLK,IT)
C
C-----------------------------------------------------------------------
C
               CALL CINIT(NCPLWFMAX*NCPLWFMAX,W)
C
               MODIFYV = .FALSE.
C
               IF ( ORBPOLSS.NE.'NONE      ' ) THEN
                  DO ISOL = 1,NSOLBLK(IBLK,IT)
                     L = L_LM(IKMSOLBLK(ISOL,IBLK,IT))
                     IF ( L.EQ.LOPT(IT) ) THEN
                        IF ( ORBPOLSS(1:5).EQ.'SIGMA' .OR. ORBPOLSS(1:4)
     &                       .EQ.'DMFT' .OR. ORBPOLSS(1:5).EQ.'LDA+U' )
     &                       MODIFYV = .TRUE.
                        EXIT
                     END IF
                  END DO
               END IF
C
C-----------------------------------------------------------------------
               IF ( .NOT.(MODIFYV) ) THEN
                  LDAU = .FALSE.
C
                  DO IR = 1,JRCRI(IM)
                     CVT(IR) = VT(IR,IT)
                     CBT(IR) = BT(IR,IT)
                  END DO
C
                  DO IPOT = 2,NFPT(IT)
                     LM = LMIFP(IPOT,IT)
                     DO IR = JRNS1(IM),JRCRI(IM)
                        CVNST(IR,IPOT) = VNST(IR,LM,IT)
                        CBNST(IR,IPOT) = BNST(IR,LM,IT)
                     END DO
                  END DO
C
               ELSE IF ( ORBPOLSS(1:4).EQ.'DMFT' .OR. ORBPOLSS(1:5)
     &                   .EQ.'LDA+U' ) THEN
                  LDAU = .FALSE.
C -------------------------------------------- copy interaction matrix W
C
                  DO IR = 1,JRCRI(IM)
                     CVT(IR) = VT(IR,IT)
                     CBT(IR) = BT(IR,IT)
                  END DO
C
                  DO IPOT = 2,NFPT(IT)
                     LM = LMIFP(IPOT,IT)
                     DO IR = JRNS1(IM),JRCRI(IM)
                        CVNST(IR,IPOT) = VNST(IR,LM,IT)
                        CBNST(IR,IPOT) = BNST(IR,LM,IT)
                     END DO
                  END DO
C
                  SUMW = 0D0
                  ILMBOT = (LOPT(IT)**2+1)
                  ILMTOP = ILMBOT - 1 + (2*LOPT(IT)+1)
                  W = C0
C
                  DO J = 1,NSOLBLK(IBLK,IT)
                     JLM = IKMSOLBLK(J,IBLK,IT)
                     IF ( JLM.GE.ILMBOT .AND. JLM.LE.ILMTOP ) THEN
                        I1 = JLM + (NL**2)*(IS-1)
                        DO I = 1,NSOLBLK(IBLK,IT)
                           ILM = IKMSOLBLK(I,IBLK,IT)
                           IF ( ILM.GE.ILMBOT .AND. ILM.LE.ILMTOP ) THEN
                              I2 = ILM + (NL**2)*(IS-1)
                              W(I,J) = DMFTSIG(I2,I1,IT)
                              SUMW = SUMW + ABS(W(I,J))
                              IF ( ABS(W(I,J)).LT.TOL ) W(I,J) = C0
                           END IF
                        END DO
                     END IF
                  END DO
                  IF ( ABS(SUMW).GT.TOL ) LDAU = .TRUE.
C
               ELSE IF ( ORBPOLSS(1:5).EQ.'SIGMA' ) THEN
C
                  DO IR = 1,JRCRI(IM)
                     CVT(IR) = VT(IR,IT) + SEVT(IR,IT)
                     CBT(IR) = BT(IR,IT) + SEBT(IR,IT)
                  END DO
                  DO LM = 2,NLMFPT(IT)
                     IF ( KLMFP(LM,IT).NE.0 ) THEN
                        DO IR = JRNS1(IM),JRCRI(IM)
                           CVNST(IR,LM) = VNST(IR,LM,IT)
     &                        + SEVNST(IR,LM,IT)
                           CBNST(IR,LM) = BNST(IR,LM,IT)
     &                        + SEBNST(IR,LM,IT)
                        END DO
                     END IF
                  END DO
C
               ELSE
C
                  CALL STOP_MESSAGE(ROUTINE,
     &                              'MODIFYV = T  AND  ORBPOL != SIGMA')
C
               END IF
C
C-----------------------------------------------------------------------
C              set up spin-dependent spherical potential
C-----------------------------------------------------------------------
C
               DO IR = 1,JRCRI(IM)
                  CV0(IR) = CVT(IR) + SPNWGT*CBT(IR)
               END DO
C
C-----------------------------------------------------------------------
C  get l-dependent wave functions for spherical potential
C  SOLVER = BS:   return only UN-NORMALIZED
C                 regular RWFs (PR0LS(IR,IL,IS),QRL)
C                 without factor  r^gamma
C           BORN  return all NORMALIZED RWFs
C                 including factors  r^gamma  and r^(-gamma)
C-----------------------------------------------------------------------
               IF ( IBLK.EQ.1 ) THEN
C
                  IF ( SOLVER_FP.EQ.'BS' ) THEN
                     NRSPHERSOL = JRNS1(IM)
                  ELSE
                     NRSPHERSOL = JRCRI(IM)
                  END IF
C
                  CALL RWFBS_NORM(IT,NRSPHERSOL,GET_ALL_SPHER_SOL,C,
     &                            ERYD,P,IS,CV0,PR0LS(1,1,IS),QRL,
     &                            PH0LS(1,1,IS),QHL,TMT0L,SMT0L,GAMMA,
     &                            IGAMMA,RGAM,RIGAM,GF_CONV_RH,PJL)
C
               END IF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
               IF ( SOLVER_FP.EQ.'BS' ) THEN
C
                  CALL FPRWFBS(IREL,SPNWGT,GETIRRSOL,C,ERYD,P,CV0,CVNST,
     &                         CBNST,LDAU,W,NFPT(IT),JRNS1(IM),R(1,IM),
     &                         DRDI(1,IM),VAMEG(1,1,1,1,IBLK,IT),IBLK,
     &                         NSOLBLK(1,IT),IKMSOLBLK(1,1,IT),
     &                         JRCUT(0,IM),NPAN(IM),PRM,QRM,PIM,QIM,
     &                         PH0LS(1,1,IS),QHL,GAMMA(0,IT),
     &                         IGAMMA(0,IT),RGAM(1,0,IT),RIGAM(1,0,IT),
     &                         DPRIRTOP,GF_CONV_RH)
C
C --------------------- GENERATE T-MATRIX FROM CALCULATED WAVE FUNCTIONS
C
                  DO J = 1,NSOL
                     LJ = L_LM(IKMSOLBLK(J,IBLK,IT))
                     DO I = 1,NSOL
C
                        GR(I,J) = PRM(I,J,IRTOP)/R(IRTOP,IM)
C
                        DGR(I,J) = (DPRIRTOP(I,J)/DRDI(IRTOP,IM)-GR(I,J)
     &                             )/R(IRTOP,IM)
C
                     END DO
                  END DO
C
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        LI = L_LM(IKMSOLBLK(I,IBLK,IT))
C
                        XBLK(I,J) = CHLM(LI)*DGR(I,J) - DCHLM(LI)
     &                              *GR(I,J)
                        TBLK(I,J) = CHLP(LI)*DGR(I,J) - DCHLP(LI)
     &                              *GR(I,J)
                     END DO
                  END DO
C
                  CALL CMATINV3(NSOL,NCPLWFMAX,IPIVKM,TBLK,WKM2,YBLKINV)
C
                  CFAC = -CI/(2*P)
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        YBLK(I,J) = CFAC*(TBLK(I,J)+XBLK(I,J))
                     END DO
                  END DO
C
                  CALL CMATMUL(NSOL,NCPLWFMAX,YBLK,YBLKINV,TBLK)
C
               ELSE
C
                  PH0LS(:,:,IS) = PJL(:,:)
                  CALL FPRWFBORN(SPNWGT,GETIRRSOL,IT,CVNST,CBNST,LDAU,W,
     &                           VAMEG(1,1,1,1,IBLK,IT),IBLK,PRM,PIM,
     &                           PR0LS(1,1,IS),PH0LS(1,1,IS),SBLK,SMT0L,
     &                           TBLK,TMT0L,CBLK,GBLK,L_LM,DPRIRTOP)
C
                  DO J = 1,NSOL
                     JLMS = IKMSOLBLK(J,IBLK,IT)
                     DO I = 1,NSOL
                        ILMS = IKMSOLBLK(I,IBLK,IT)
                        CHIS(ILMS,JLMS,IS) = CBLK(I,J)
                        GAMS(ILMS,JLMS,IS) = GBLK(I,J)
                     END DO
                  END DO
C
               END IF
C-----------------------------------------------------------------------
C
               CALL CMATINV3(NSOL,NCPLWFMAX,IPIVKM,TBLK,WKM2,YBLKINV)
C
               DO J = 1,NSOL
                  JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                  DO I = 1,NSOL
                     ILMS = LMSOFF + IKMSOLBLK(I,IBLK,IT)
                     MSST(ILMS,JLMS,IT) = YBLKINV(I,J)
                     TSST(ILMS,JLMS,IT) = TBLK(I,J)
                  END DO
               END DO
C
C ---------------------------- remove scaling factor from wave functions
C
               DO IR = 1,IRTOP
                  RSCLP = 1D0/R(IR,IM)
C
                  DO J = 1,NSOL
                     JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                     DO I = 1,NSOL
                        PRM(I,J,IR) = PRM(I,J,IR)*RSCLP
                        JG(IR,I,JLMS) = PIM(I,J,IR)*RSCLP
                     END DO
                  END DO
               END DO
C
               IF ( WRONSKI .AND. SOLVER_FP.EQ.'BS' ) THEN
                  DO IR = 1,IRTOP
                     RSCLP = 1D0/R(IR,IM)
                     RSCLQ = RSCLP/C
C
                     DO J = 1,NSOL
                        JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                        DO I = 1,NSOL
                           QRM(I,J,IR) = QRM(I,J,IR)*RSCLQ
                           JF(IR,I,JLMS) = QIM(I,J,IR)*RSCLQ
                        END DO
                     END DO
                  END DO
               END IF
C
C transformation of wave functions ------- set up matrix of coefficients
C
               XBLK(:,:) = C0
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
               IF ( .NOT.GF_CONV_RH ) THEN
C
                  DO I = 1,NSOL
                     LI = L_LM(IKMSOLBLK(I,IBLK,IT))
                     CJL = 0.5D0*(CHLP(LI)+CHLM(LI))
                     CHL = CHLP(LI)
                     DO J = 1,NSOL
                        JLM = IKMSOLBLK(J,IBLK,IT)
                        YBLK(I,J) = PRM(I,J,IRTOP)
                        XBLK(I,J) = CJL*YBLKINV(I,J)
                        IF ( J.EQ.I ) XBLK(I,J) = XBLK(I,J) - CI*P*CHL
                     END DO
                  END DO
C
C ---------- solve system of linear equations determining transformation
C
                  CALL ZGETRF(NSOL,NSOL,YBLK,NCPLWFMAX,IPIV,INFO)
                  CALL ZGETRS('N',NSOL,NSOL,YBLK,NCPLWFMAX,IPIV,XBLK,
     &                        NCPLWFMAX,INFO)
C
                  DO J = 1,NSOL
                     JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                     NCPLWF(JLMS) = NSOL
                     DO I = 1,NSOL
                        IKMCPLWF(I,JLMS) = IKMSOLBLK(I,IBLK,IT)
                        DO IR = 1,IRTOP
                           DO K = 1,NSOL
                              ZG(IR,I,JLMS) = ZG(IR,I,JLMS)
     &                           + PRM(I,K,IR)*XBLK(K,J)
                           END DO
                        END DO
                     END DO
                  END DO
C
               ELSE
C                                                basis functions R and H
C
                  DO J = 1,NSOL
                     JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                     NCPLWF(JLMS) = NSOL
                     DO I = 1,NSOL
                        IKMCPLWF(I,JLMS) = IKMSOLBLK(I,IBLK,IT)
                        DO IR = 1,IRTOP
                           ZG(IR,I,JLMS) = ZG(IR,I,JLMS) + PRM(I,J,IR)
                        END DO
                     END DO
                  END DO
C
               END IF
C
               IF ( WRONSKI .AND. SOLVER_FP.EQ.'BS' ) THEN
                  DO J = 1,NSOL
                     JLMS = LMSOFF + IKMSOLBLK(J,IBLK,IT)
                     NCPLWF(JLMS) = NSOL
                     DO I = 1,NSOL
                        IKMCPLWF(I,JLMS) = IKMSOLBLK(I,IBLK,IT)
                        DO IR = 1,IRTOP
                           DO K = 1,NSOL
                              ZF(IR,I,JLMS) = ZF(IR,I,JLMS)
     &                           + QRM(I,K,IR)*XBLK(K,J)
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
C
C --------------------- set alpha matrix SSST required for Lloyd formula
C
               DO J = 1,NSOL
                  JLM = IKMSOLBLK(J,IBLK,IT)
                  JLMS = LMSOFF + JLM
                  LJ = L_LM(JLM)
                  ZBOT = P*R(1,IM)
                  CJL_J = CJLZ(LJ,ZBOT)
C
                  DO I = 1,NSOL
                     IKMCPLWF(I,JLMS) = IKMSOLBLK(I,IBLK,IT)
                     ILMS = LMSOFF + IKMSOLBLK(I,IBLK,IT)
C
                     SSST(ILMS,JLMS,IT) = ZG(1,I,JLMS)/CJL_J
C
                  END DO
               END DO
C
            END DO
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C
C ------------------------ test wave functions:  WRONSKI test (1st kind)
C
            IF ( WRONSKI ) THEN
               WRITE (6,99002) IT,JLMS,ERYD
               DO JLMS = LMSOFF + 1,LMSOFF + NLM
                  WRITE (6,*) ' '
                  DO IR = 1,IRTOP
                     CRSQ = WFAC*R(IR,IM)**2
                     CSUM = C0
                     DO I = 1,NCPLWF(JLMS)
                        CSUM = CSUM + 
     &                         (ZF(IR,I,JLMS)*JG(IR,I,JLMS)-ZG(IR,I,
     &                         JLMS)*JF(IR,I,JLMS))*CRSQ
                     END DO
                     WRITE (6,'(3I4,2F20.12)') IT,JLMS,IR,CSUM
                  END DO
               END DO
            END IF
C
C ----------------------------- matrix elements ---  not needed for E(k)
C
            CALL FPNRMATELM(IT,IS,LMSOFF,ZG,JG,MEZZ,MEZJ,MEZZL,MEZJL)
C
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         DO ILMS = 1,NKM
C
C ------------------------------------------- write wavefunctions for IT
C
            IF ( IWRREGWF.NE.0 ) WRITE (IFILSS,REC=ILMS+(IT-1)*NKM) IT,
     &                                  'REG',ILMS,NCPLWF(ILMS),
     &                                  (IKMCPLWF(K,ILMS),K=1,
     &                                  NCPLWF(ILMS)),
     &                                  ((ZG(I,K,ILMS),I=1,IRTOP),K=1,
     &                                  NCPLWFMAX)
C
            IF ( IWRIRRWF.NE.0 ) WRITE (IFILSS,REC=ILMS+(IT-1+NT)*NKM)
     &                                  IT,'IRR',ILMS,
     &                                  ((JG(I,K,ILMS),I=1,IRTOP),K=1,
     &                                  NCPLWFMAX)
C
         END DO
C
         IF ( SOLVER_FP.EQ.'BORN' .AND. IWRREGWF.NE.0 .AND. 
     &        IFILSS.NE.IFILGFWF ) WRITE (IFILCBWF_SPH,REC=IT)
     &        (((CHIS(I,J,IS),I=1,NLM),J=1,NLM),
     &        ((GAMS(I,J,IS),I=1,NLM),J=1,NLM),
     &        ((PR0LS(IR,IL,IS)/R(IR,IM),PH0LS(IR,IL,IS)/R(IR,IM),IR=1,
     &        IRTOP),IL=1,NL),IS=1,NSPIN)
C
C-----------------------------------------------------------------------
C               store l-projected matrix elements
C-----------------------------------------------------------------------
         IF ( L_PROJECTED_ME .AND. .NOT.CALC_GFWF ) WRITE (IFILMEZZL)
     &        IT,MEZZL,MEZJL
C-----------------------------------------------------------------------
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                  CHECK
         IF ( IPRINT.GE.3 ) THEN
            WRITE (6,99005) 'atom type   ',IT,ERYD,CTL(IT,1)
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,IREL,
     &                      IREL,1,1.0D-9,6)
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,6)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,6)
         END IF
C
         IF ( DUMP_SSITE ) THEN
            IF ( IT.EQ.1 ) THEN
               LFILNAM = LEN_TRIM(SOLVER_FP)
               FILNAM = 'zzzzzz_ssite_'//SOLVER_FP(1:LFILNAM)
               LFILNAM = LFILNAM + 13
               IF ( GF_CONV_RH ) THEN
                  FILNAM = FILNAM(1:LFILNAM)//'_RH.dat'
               ELSE
                  FILNAM = FILNAM(1:LFILNAM)//'_ZJ.dat'
               END IF
               CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:(LFILNAM+7)))
            END IF
            WRITE (IOTMP,99005) 'atom type   ',IT,ERYD,CTL(IT,1),
     &                          SOLVER_FP,GF_CONV_RH
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,IREL,
     &                      IREL,1,1.0D-9,IOTMP)
            DO I = 1,NKM
               DO J = 1,I
                  WRITE (IOTMP,'(2I3,2F20.12,2X,A)') I,J,TSST(I,J,IT),
     &                   SOLVER_FP
                  WRITE (IOTMP,'(2I3,2F20.12,2X,A)') J,I,TSST(J,I,IT),
     &                   SOLVER_FP
               END DO
            END DO
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IOTMP)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IOTMP)
            DO I = 1,NKM
               DO J = 1,NKM
                  WRITE (IOTMP,99004) I,J,TSST(I,J,IT),MEZZ(I,J,IT,1),
     &                                MEZJ(I,J,IT,1)
               END DO
            END DO
C
            IF ( IT.EQ.ITTOP ) THEN
               CLOSE (IOTMP)
               WRITE (6,99006) FILNAM(1:(LFILNAM+7))
            END IF
         END IF
C                                                                  CHECK
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C                 symetrize ALL single site matrices
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
CC    CALL SSITE_SYMMETRIZE('TSST',TSST)
CC    CALL SSITE_SYMMETRIZE('MSST',MSST)
Cc    CALL SSITE_SYMMETRIZE('SSST',SSST)
C
C=======================================================================
C
      DEALLOCATE (JF,JG,GR,ZF,ZG,CBT,DGR,XBLK,YBLK,PIM,CVT)
      DEALLOCATE (PRM,TBLK,CHLM,CHLP,YBLKINV,IPIV,DCHLM,DCHLP,CBNST)
      DEALLOCATE (CVNST,DPRIRTOP)
      IF ( SOLVER_FP.EQ.'BS' ) DEALLOCATE (QRM,QIM)
C
      IF ( WRONSKI ) CALL STOP_MESSAGE(ROUTINE,'check of Wronskian done'
     &                                 )
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) CALL STOP_BREAKPOINT(ROUTINE)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C
99001 FORMAT (4I4,2E25.14)
99002 FORMAT (/,' wronski relation for IT =',I3,' JLMS =',I3,' E =',
     &        2F10.5)
99003 FORMAT (/,' <<DEBUG>>',/,' <FPNRSSITE>:  IWRREGWF =',I3,
     &        '  IWRIRRWF =',I3,'  IFILSS =',I3,'  GETIRRSOL =',L3,/,
     &        15X,'ORBPOLSS = ',A10,6X,'ERYD =',2F15.10,/)
99004 FORMAT (2I3,2X,'TSST ',2E25.15,/,8X,'MEZZ ',2E25.15,/,8X,'MEZJ ',
     &        2E25.15,/)
99005 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),//,10X,'energy',12X,
     &        2F12.6,/,10X,'speed of light',F16.6,/,:,10X,
     &        'solver            ',A,/,10X,'GF_CONV_RH        ',L10,/)
99006 FORMAT (//,5X,'>>>>>  DUMP written to file  ',A,//)
      END
C*==fpnrmatelm.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPNRMATELM(IT,IS,LMSOFF,ZG,JG,MEZZ,MEZJ,MEZZL,MEZJL)
C   ********************************************************************
C   *                                                                  *
C   *  calculate all necessary matrix elements                         *
C   *               MEZZ =  < Z_LAM | A | Z_LAM' >                     *
C   *               MEZJ =  < Z_LAM | A | J_LAM' >                     *
C   *  with A = 1, sig_z             (index IME)                       *
C   *                                                                  *
C   *  L_PROJECTED_ME = T:  calculate the l-diagonal part MEZZL, MEZJL *
C   *                     to be used in the FP-DOS calculation         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:L_PROJECTED_ME
      USE MOD_TYPES,ONLY:IMT,NTMAX,NCPLWFMAX,IKMCPLWF
      USE MOD_RMESH,ONLY:R2DRDI_W_RADINT,NRMAX,JRCUT,NPAN,NM,NSF,LMISF,
     &    WINTLM
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,NSPIN,NLM,NLMAX,NMEMAX,NKMMAX,
     &    NCPLWF,NL,IDOS,ISMT
      IMPLICIT NONE
C*--FPNRMATELM949
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPNRMATELM')
C
C Dummy arguments
C
      INTEGER IS,IT,LMSOFF
      COMPLEX*16 JG(NRMAX,NCPLWFMAX,NKMMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),
     &           ZG(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 G123,W,ZZGNT(:,:,:,:)
      REAL*8 GAUNT_RYLM
      INTEGER I1,I2,IA_ERR,ICALL,IL1,IL2,ILM1,ILM2,IM,ISF,J,JLMS1,JLMS2,
     &        JOFF,JSF,L1,L2,LM,LM1,LM2,LMAX,LSF,M1,M2,MSF,NSFZZMAX
      COMPLEX*16 SIRR,SREG,ZGJG,ZGZG
      SAVE NSFZZMAX,ZZGNT
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE ZZGNT
C
      CALL TRACK_INFO(ROUTINE)
C
C=======================================================================
C                           INITIALIZE
C=======================================================================
      IF ( ICALL.EQ.0 ) THEN
C
         ICALL = 1
C
C-----------------------------------------------------------------------
C             Gaunt coefficients for REAL spherical harmonics
C-----------------------------------------------------------------------
C
         LMAX = NL - 1
C
         NSFZZMAX = 0
         DO IM = 1,NM
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               LSF = L_LM(LM)
               IF ( LSF.LE.2*LMAX ) NSFZZMAX = MAX(NSFZZMAX,ISF)
            END DO
         END DO
C
         ALLOCATE (ZZGNT(NLM,NLM,NSFZZMAX,NM),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZZGNT')
C
         DO IM = 1,NM
            DO ISF = 1,MIN(NSFZZMAX,NSF(IM))
               LM = LMISF(ISF,IM)
               LSF = L_LM(LM)
               MSF = M_LM(LM)
               DO LM2 = 1,NLM
                  L2 = L_LM(LM2)
                  M2 = M_LM(LM2)
                  DO LM1 = 1,NLM
                     L1 = L_LM(LM1)
                     M1 = M_LM(LM1)
                     ZZGNT(LM1,LM2,ISF,IM)
     &                  = GAUNT_RYLM(L1,M1,L2,M2,LSF,MSF)
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C=======================================================================
C
      IF ( L_PROJECTED_ME .AND. IS.EQ.1 ) THEN
         MEZZL(:,:,:,:) = C0
         MEZJL(:,:,:,:) = C0
      END IF
C
C=======================================================================
C
      IM = IMT(IT)
C
      DO JLMS2 = LMSOFF + 1,LMSOFF + NLM
         DO I2 = 1,NCPLWF(JLMS2)
            ILM2 = IKMCPLWF(I2,JLMS2)
            IL2 = L_LM(ILM2) + 1
C
            DO JLMS1 = LMSOFF + 1,LMSOFF + NLM
               DO I1 = 1,NCPLWF(JLMS1)
                  ILM1 = IKMCPLWF(I1,JLMS1)
                  IL1 = L_LM(ILM1) + 1
C
C-----------------------------------------------------------------------
C               integrals within muffin-tin regime
C-----------------------------------------------------------------------
C
                  IF ( ILM1.EQ.ILM2 ) THEN
                     SREG = C0
                     IF ( JLMS1.EQ.JLMS2 ) THEN
                        SIRR = C0
                        DO J = 1,JRCUT(1,IM)
                           W = R2DRDI_W_RADINT(J,IM)
                           ZGZG = ZG(J,I2,JLMS2)*ZG(J,I1,JLMS1)
                           ZGJG = ZG(J,I2,JLMS2)*JG(J,I1,JLMS1)
                           SREG = SREG + W*ZGZG
                           SIRR = SIRR + W*ZGJG
                        END DO
                        MEZJ(JLMS1,JLMS2,IT,IDOS)
     &                     = MEZJ(JLMS1,JLMS2,IT,IDOS) + SIRR
                        IF ( L_PROJECTED_ME ) THEN
                           IF ( IL1.EQ.IL2 ) MEZJL(JLMS1,JLMS2,IL1,IDOS)
     &                          = MEZJL(JLMS1,JLMS2,IL1,IDOS) + SIRR
                        END IF
                     ELSE
                        DO J = 1,JRCUT(1,IM)
                           W = R2DRDI_W_RADINT(J,IM)
                           ZGZG = ZG(J,I2,JLMS2)*ZG(J,I1,JLMS1)
                           SREG = SREG + W*ZGZG
                        END DO
                     END IF
                     MEZZ(JLMS1,JLMS2,IT,IDOS)
     &                  = MEZZ(JLMS1,JLMS2,IT,IDOS) + SREG
                     IF ( L_PROJECTED_ME ) THEN
                        IF ( IL1.EQ.IL2 ) MEZZL(JLMS1,JLMS2,IL1,IDOS)
     &                       = MEZZL(JLMS1,JLMS2,IL1,IDOS) + SREG
                     END IF
                  END IF
C
C-----------------------------------------------------------------------
C            integrals within interstitial regime
C-----------------------------------------------------------------------
C
                  DO ISF = 1,MIN(NSFZZMAX,NSF(IM))
                     G123 = ZZGNT(ILM1,ILM2,ISF,IM)
                     IF ( ABS(G123).GT.1D-8 ) THEN
C
                        JOFF = JRCUT(1,IM)
                        SREG = C0
                        IF ( JLMS1.EQ.JLMS2 ) THEN
                           SIRR = C0
                           DO J = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                              JSF = J - JOFF
                              W = WINTLM(JSF,ISF,IM)
                              ZGZG = ZG(J,I2,JLMS2)*ZG(J,I1,JLMS1)
                              ZGJG = ZG(J,I2,JLMS2)*JG(J,I1,JLMS1)
                              SREG = SREG + W*ZGZG
                              SIRR = SIRR + W*ZGJG
                           END DO
                           MEZJ(JLMS1,JLMS2,IT,IDOS)
     &                        = MEZJ(JLMS1,JLMS2,IT,IDOS) + SIRR*G123
                           IF ( L_PROJECTED_ME .AND. ISF.EQ.1 ) THEN
                              IF ( IL1.EQ.IL2 )
     &                             MEZJL(JLMS1,JLMS2,IL1,IDOS)
     &                             = MEZJL(JLMS1,JLMS2,IL1,IDOS)
     &                             + SIRR*G123
                           END IF
                        ELSE
                           DO J = JRCUT(1,IM) + 1,JRCUT(NPAN(IM),IM)
                              JSF = J - JOFF
                              W = WINTLM(JSF,ISF,IM)
                              ZGZG = ZG(J,I2,JLMS2)*ZG(J,I1,JLMS1)
                              SREG = SREG + W*ZGZG
                           END DO
                        END IF
                        MEZZ(JLMS1,JLMS2,IT,IDOS)
     &                     = MEZZ(JLMS1,JLMS2,IT,IDOS) + SREG*G123
                        IF ( L_PROJECTED_ME .AND. ISF.EQ.1 ) THEN
                           IF ( IL1.EQ.IL2 ) MEZZL(JLMS1,JLMS2,IL1,IDOS)
     &                          = MEZZL(JLMS1,JLMS2,IL1,IDOS)
     &                            + SREG*G123
                        END IF
                     END IF
                  END DO
C
               END DO
            END DO
C
         END DO
      END DO
C
      IF ( NSPIN.EQ.1 ) RETURN
C
C=======================================================================
C                     deal with spin polarisation
C=======================================================================
C
      IF ( IS.EQ.1 ) THEN
         MEZZ(1:NLM,1:NLM,IT,ISMT) = -MEZZ(1:NLM,1:NLM,IT,IDOS)
         MEZJ(1:NLM,1:NLM,IT,ISMT) = -MEZJ(1:NLM,1:NLM,IT,IDOS)
      ELSE
         I1 = NLM + 1
         I2 = NLM + NLM
         MEZZ(I1:I2,I1:I2,IT,ISMT) = MEZZ(I1:I2,I1:I2,IT,IDOS)
         MEZJ(I1:I2,I1:I2,IT,ISMT) = MEZJ(I1:I2,I1:I2,IT,IDOS)
      END IF
C
C-----------------------------------------------------------------------
C
      IF ( L_PROJECTED_ME ) THEN
C
         IF ( IS.EQ.1 ) THEN
            MEZZL(1:NLM,1:NLM,1:NL,ISMT) = -MEZZL(1:NLM,1:NLM,1:NL,IDOS)
            MEZJL(1:NLM,1:NLM,1:NL,ISMT) = -MEZJL(1:NLM,1:NLM,1:NL,IDOS)
         ELSE
            I1 = NLM + 1
            I2 = NLM + NLM
            MEZZL(I1:I2,I1:I2,1:NL,ISMT) = MEZZL(I1:I2,I1:I2,1:NL,IDOS)
            MEZJL(I1:I2,I1:I2,1:NL,ISMT) = MEZJL(I1:I2,I1:I2,1:NL,IDOS)
         END IF
C
      END IF
C
      END
