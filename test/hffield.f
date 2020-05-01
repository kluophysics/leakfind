C*==hff_split.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE HFF_SPLIT(NWRLOG,IEPATH,DEFERMI,IPRINT,MSST,TAUT,
     &                     IECURR,OBS_LT,HINT)
C   ********************************************************************
C   *                                                                  *
C   *  calculates matrix elements of several hyperfine interaction     *
C   *  connected quantities like elctric field gradient etc.           *
C   *  the individual quantities are stored in arrays, where           *
C   *  HINT is saved for the next energy point and for energy          *
C   *  integration, respectively. all the related arrays have a        *
C   *  counting index as the last index of the array indicates the     *
C   *  corresponding physical property.                                *
C   *                                                                  *
C   *  index-list for matrix elements                                  *
C   *                                                                  *
C   *  1      expectation value of (1/r)^3                             *
C   *  2      EFG                                                      *
C   *  3      nulcear-spin - electron-spin - contact hyperfine field   *
C   *  4      nuclear-spin - electron-spin  hyperfine field            *
C   *  5      nuclear-spin - electron-orbit hyperfine field            *
C   *                                                                  *
C   *  revised 2014  HE                                                *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,JRCRI,FULLPOT
      USE MOD_ANGMOM,ONLY:IDOS,NL,NMEHFMAX,NOBSMAX,NLMAX,NKMMAX,QMOFF,
     &    QOFF,QMDIA,QDIA,SMOFF,SOFF,SMDIA,SDIA,NKM,TXT_L,WKM1,WKM2
      USE MOD_TYPES,ONLY:NT,NTMAX,Z,IMT,TXT_T,NLT,NCPLWFMAX
      USE MOD_ENERGY,ONLY:NEPATH,WETAB,NETAB
      USE MOD_FILES,ONLY:IFILCBWF,IFILLOG
      USE MOD_CALCMODE,ONLY:IREL,LHS_SOL_EQ_RHS_SOL
      USE MOD_CONSTANTS,ONLY:A0_CGS,MB_CGS,PI,C0,C1
      IMPLICIT NONE
C*--HFF_SPLIT34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='HFF_SPLIT')
      INTEGER IHFF
      PARAMETER (IHFF=4)
      REAL*8 F1,F2
      PARAMETER (F1=1.0D0,F2=2.0D0*MB_CGS/(A0_CGS*A0_CGS*A0_CGS))
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C Dummy arguments
C
      REAL*8 DEFERMI
      INTEGER IECURR,IEPATH,IPRINT,NWRLOG
      COMPLEX*16 HINT(NLMAX,NTMAX,NMEHFMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX)
      REAL*8 OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX)
C
C Local variables
C
      REAL*8 AHF_F(:,:,:),AHF_G(:,:,:),BSUM,BSUML,DOVR(NRMAX),
     &       DROVR1(NRMAX),DROVRN(NRMAX),EVSUM,H(NTMAX,NMEHFMAX),
     &       HI(NTMAX,NMEHFMAX),MJ,RNUC
      COMPLEX*16 CMATTRC
      COMPLEX*16 HL(NMEHFMAX),HTOT(:,:),JF(:,:,:),JFR(:,:,:),JG(:,:,:),
     &           JGR(:,:,:),TMAT(:,:),ZFL(:,:,:),ZFR(:,:,:),ZFZF(:,:),
     &           ZFZF1(:,:),ZFZFC(:,:),ZFZFD(:,:),ZFZFO(:,:),ZGL(:,:,:),
     &           ZGR(:,:,:),ZGZG(:,:),ZGZG1(:,:),ZGZGC(:,:),ZGZGD(:,:),
     &           ZGZGO(:,:),ZHFIZ(:,:,:)
      INTEGER IA_ERR,IBSPLIT,IFIL_LHS,IFIL_RHS,IKM1,IKM2,
     &        IKMCPLWF_T(:,:,:),IL,IM,IME,IR,IRTOP,IR_NUC,IT,IW,IWRLOG,
     &        K1,K2,KAP1,KAP2,L,LAM1,LAM2,M,MUE,MUETOP,N,NCPLWF_T(:,:),
     &        NSOL
      INTEGER IKAPMUE
      LOGICAL INITIALIZE
      REAL*8 RNUCTAB
      SAVE AHF_F,AHF_G,HTOT
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE NCPLWF_T,IKMCPLWF_T,ZHFIZ,TMAT,HTOT
      ALLOCATABLE ZGL,ZFL
      ALLOCATABLE JFR,JGR,ZGR,ZFR,JG,JF
      ALLOCATABLE AHF_G,AHF_F
      ALLOCATABLE ZGZG,ZGZG1,ZGZGC,ZGZGD,ZGZGO
      ALLOCATABLE ZFZF,ZFZF1,ZFZFC,ZFZFD,ZFZFO
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (ZGZG(NCPLWFMAX,NCPLWFMAX),ZFZF(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (ZGZG1(NCPLWFMAX,NCPLWFMAX),ZFZF1(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (ZGZGC(NCPLWFMAX,NCPLWFMAX),ZFZFC(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (ZGZGD(NCPLWFMAX,NCPLWFMAX),ZFZFD(NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (ZGZGO(NCPLWFMAX,NCPLWFMAX),ZFZFO(NCPLWFMAX,NCPLWFMAX))
C
      ALLOCATE (ZGL(NRMAX,NCPLWFMAX,NKM),ZFL(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JFR(NRMAX,NCPLWFMAX,NKM),JGR(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGR(NRMAX,NCPLWFMAX,NKM),ZFR(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKM),JG(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (ZHFIZ(NKMMAX,NKMMAX,5),TMAT(NKMMAX,NKMMAX))
      ALLOCATE (NCPLWF_T(NKMMAX,NTMAX))
      ALLOCATE (IKMCPLWF_T(NCPLWFMAX,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZF')
C
      IF ( IREL.NE.3 ) RETURN
C
      IF ( (IEPATH.EQ.NEPATH) .AND. (IECURR.EQ.NETAB(1)) )
     &     WRITE (6,99001)
C
C     if first energy point then initialize integrated quantities
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         HINT(:,:,:) = C0
         ALLOCATE (HTOT(NTMAX,5))
         HTOT(:,:) = C0
C
         ALLOCATE (AHF_G(NKMMAX,NKMMAX,5),AHF_F(NKMMAX,NKMMAX,5))
         AHF_G(:,:,:) = 0D0
         AHF_F(:,:,:) = 0D0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,NL - 1
            IL = L + 1
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C
            MUETOP = 2*L + 2
            MJ = -(DBLE(L)+0.5D0) - 1D0
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MUE = 1,MUETOP
               MJ = MJ + 1D0
               IF ( ABS(MJ).GE.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C-----------------------------------------------------------------------
               IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C-----------------------------------------------------------------------
C
C <1/R^3> -----------------------------------------------------------111
C
               IME = 1
               AHF_G(IKM1,IKM1,IME) = 1D0
               AHF_G(IKM2,IKM2,IME) = 1D0
               AHF_F(IKM1,IKM1,IME) = 1D0
               AHF_F(IKM2,IKM2,IME) = 1D0
C
C coefficients for the quadrupolar matrix elements ------------------222
C
               IME = 2
               AHF_G(IKM1,IKM1,IME) = QDIA(IKM1)
               AHF_G(IKM2,IKM2,IME) = QDIA(IKM2)
               AHF_G(IKM1,IKM2,IME) = QOFF(IKM1)
               AHF_G(IKM2,IKM1,IME) = QOFF(IKM1)
               AHF_F(IKM1,IKM1,IME) = QMDIA(IKM1)
               AHF_F(IKM2,IKM2,IME) = QMDIA(IKM2)
               AHF_F(IKM1,IKM2,IME) = QMOFF(IKM1)
               AHF_F(IKM2,IKM1,IME) = QMOFF(IKM1)
C
C coefficients to calculate the spin-spin field ---------------------333
C
               IME = 3
               AHF_G(IKM1,IKM1,IME) = -MJ/(KAP1+0.5D0)
               AHF_F(IKM1,IKM1,IME) = -MJ/(-KAP1+0.5D0)
               IF ( NSOL.EQ.2 ) THEN
                  AHF_G(IKM1,IKM2,IME)
     &               = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
                  AHF_G(IKM2,IKM1,IME) = AHF_G(IKM1,IKM2,IME)
                  AHF_G(IKM2,IKM2,IME) = -MJ/(KAP2+0.5D0)
                  AHF_F(IKM2,IKM2,IME) = -MJ/(-KAP2+0.5D0)
               END IF
C
C coefficients for the spin-dipolar matrix elements -----------------444
C
               IME = 4
               AHF_G(IKM1,IKM1,IME) = SDIA(IKM1)
               AHF_F(IKM1,IKM1,IME) = SMDIA(IKM1)
               IF ( NSOL.EQ.2 ) THEN
                  AHF_G(IKM1,IKM2,IME) = SOFF(IKM1)
                  AHF_G(IKM2,IKM1,IME) = AHF_G(IKM1,IKM2,IME)
                  AHF_G(IKM2,IKM2,IME) = SDIA(IKM2)
                  AHF_F(IKM1,IKM2,IME) = SMOFF(IKM1)
                  AHF_F(IKM2,IKM1,IME) = SMOFF(IKM1)
                  AHF_F(IKM2,IKM2,IME) = SMDIA(IKM2)
               END IF
C
C coefficients to calculate the orbital field -----------------------555
C
               IME = 5
               AHF_G(IKM1,IKM1,IME) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               AHF_F(IKM1,IKM1,IME) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               IF ( NSOL.EQ.2 ) THEN
                  AHF_G(IKM1,IKM2,IME)
     &               = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
                  AHF_G(IKM2,IKM1,IME) = AHF_G(IKM1,IKM2,IME)
                  AHF_G(IKM2,IKM2,IME) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
                  AHF_F(IKM2,IKM2,IME) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
               END IF
C
            END DO
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      H(:,:) = 0.0D0
      HI(:,:) = 0.0D0
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
C!!!!!!!!!!!!!!!!!!!!! TO BE CHECKED
C
Cc      CALL SET_IFIL_LHS(IFIL_RHS,IFIL_LHS)
C
      IFIL_RHS = IFILCBWF
      IFIL_LHS = IFILCBWF
C
C=======================================================================
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT:DO IT = 1,NT
C
         IF ( Z(IT).EQ.0 ) CYCLE LOOP_IT
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
         CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT)
C
         IM = IMT(IT)
C
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
C-------------------------------- prepare for finite nucleus calculation
C
         RNUC = RNUCTAB(Z(IT))
         IR = 1
         DO WHILE ( R(IR,IM).LT.RNUC )
            IR = IR + 1
         END DO
         IR_NUC = IR + 2
         IF ( MOD(IR_NUC,2).EQ.0 ) IR_NUC = IR_NUC - 1
C
C---------------------------------- NOTE: the RHS WFs will be convoluted
         DO IR = 1,IRTOP
            DOVR(IR) = 1D0/R(IR,IM)**3
            DROVR1(IR) = (R(IR,IM)/RNUC)**3/R(IR,IM)**2
            DROVRN(IR) = DROVR1(IR)/R(IR,IM)
         END DO
C
C ---------------------------------------- read in wavefunctions for RHS
C                                              use arrays *LHS as buffer
C
         CALL WAVFUN_READ_REL(IFIL_RHS,IT,0,ZGL,ZFL,JG,JF,IRTOP,
     &                        NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
C-------------- convolute the RHS wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
         CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZGL,ZFL,JG,JF,ZGR,ZFR,JGR,JFR,
     &                               NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
C --------------------------- read in wavefunctions for LHS if necessary
C
         IF ( .NOT.LHS_SOL_EQ_RHS_SOL )
     &        CALL WAVFUN_READ_REL(IFIL_LHS,IT,0,ZGL,ZFL,JG,JF,IRTOP,
     &        NCPLWF_T(1,IT),IKMCPLWF_T(1,1,IT))
C
         ZHFIZ(:,:,:) = C0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO LAM1 = 1,NKM
            DO LAM2 = 1,NKM
C
               CALL HFF_SPLIT_RADINT(IM,ZGZG,ZGL(1,1,LAM1),ZGR(1,1,LAM2)
     &                               ,DOVR,R,0.0D0,NCPLWF_T(LAM1,IT),
     &                               NCPLWF_T(LAM2,IT),IRTOP,NRMAX)
C
               CALL HFF_SPLIT_RADINT(IM,ZFZF,ZFL(1,1,LAM1),ZFR(1,1,LAM2)
     &                               ,DOVR,R,0.0D0,NCPLWF_T(LAM1,IT),
     &                               NCPLWF_T(LAM2,IT),IRTOP,NRMAX)
C
               CALL HFF_SPLIT_RADINT(IM,ZGZG1,ZGL(1,1,LAM1),
     &                               ZGR(1,1,LAM2),DOVR,R,RNUC,
     &                               NCPLWF_T(LAM1,IT),NCPLWF_T(LAM2,IT)
     &                               ,IR_NUC,NRMAX)
C
               CALL HFF_SPLIT_RADINT(IM,ZFZF1,ZFL(1,1,LAM1),
     &                               ZFR(1,1,LAM2),DOVR,R,RNUC,
     &                               NCPLWF_T(LAM1,IT),NCPLWF_T(LAM2,IT)
     &                               ,IR_NUC,NRMAX)
C
               CALL HFF_SPLIT_RADINT(IM,ZGZGC,ZGL(1,1,LAM1),
     &                               ZGR(1,1,LAM2),DROVRN,R,RNUC,
     &                               NCPLWF_T(LAM1,IT),NCPLWF_T(LAM2,IT)
     &                               ,IR_NUC,NRMAX)
C
               CALL HFF_SPLIT_RADINT(IM,ZFZFC,ZFL(1,1,LAM1),
     &                               ZFR(1,1,LAM2),DROVRN,R,RNUC,
     &                               NCPLWF_T(LAM1,IT),NCPLWF_T(LAM2,IT)
     &                               ,IR_NUC,NRMAX)
C
               DO K1 = 1,NCPLWF_T(LAM1,IT)
                  DO K2 = 1,NCPLWF_T(LAM2,IT)
C
                     ZGZGD(K1,K2) = ZGZG(K1,K2) - ZGZG1(K1,K2)
                     ZFZFD(K1,K2) = ZFZF(K1,K2) - ZFZF1(K1,K2)
C
                     ZGZGO(K1,K2) = ZGZG(K1,K2) - ZGZG1(K1,K2)
     &                              + ZGZGC(K1,K2)
                     ZFZFO(K1,K2) = ZFZF(K1,K2) - ZFZF1(K1,K2)
     &                              + ZFZFC(K1,K2)
C
                  END DO
               END DO
C
C <1/R^3> -----------------------------------------------------------111
C
               CALL HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,1,F1,ZGZG,F1,ZFZF,
     &                              AHF_F,AHF_G,NCPLWF_T(1,IT),
     &                              IKMCPLWF_T(1,1,IT))
C
C electric field gradient contribution ------------------------------222
C
               CALL HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,2,F1,ZGZG,F1,ZFZF,
     &                              AHF_F,AHF_G,NCPLWF_T(1,IT),
     &                              IKMCPLWF_T(1,1,IT))
C
C modifications for B_SSC which is zero outside the nucleus ---------333
C
               CALL HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,3,F2,ZGZGC,-F2,
     &                              ZFZFC,AHF_F,AHF_G,NCPLWF_T(1,IT),
     &                              IKMCPLWF_T(1,1,IT))
C
C modifications for B_SP which is zero inside the nucleus -----------444
C
               CALL HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,4,F1,ZGZGD,-F1,
     &                              ZFZFD,AHF_F,AHF_G,NCPLWF_T(1,IT),
     &                              IKMCPLWF_T(1,1,IT))
C
C -------------------------------------------------------------------555
C
               CALL HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,5,F2,ZGZGO,-F2,
     &                              ZFZFO,AHF_F,AHF_G,NCPLWF_T(1,IT),
     &                              IKMCPLWF_T(1,1,IT))
C
            END DO
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         N = NKM
         M = NKMMAX
C
CME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME
         DO IME = 1,NMEHFMAX
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C               perform a l-decomposition of the MEs
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO IL = 1,NL
C
               L = IL - 1
               IKM1 = 2*L**2 + 1
               IKM2 = 2*(L+1)**2
C
               WKM1(:,:) = C0
               WKM1(IKM1:IKM2,IKM1:IKM2)
     &            = ZHFIZ(IKM1:IKM2,IKM1:IKM2,IME)
C
               CALL ZGEMM('N','N',N,N,N,CPRE,WKM1,M,TMAT,M,C0,WKM2,M)
C
               HL(IME) = CMATTRC(N,M,WKM2)
C
               HINT(IL,IT,IME) = HINT(IL,IT,IME) + WETAB(IECURR,IEPATH)
     &                           *HL(IME)
               H(IT,IME) = H(IT,IME) + DIMAG(HL(IME))
               HI(IT,IME) = HI(IT,IME) + DIMAG(HINT(IL,IT,IME))
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
            CALL ZGEMM('N','N',N,N,N,CPRE,ZHFIZ(1,1,IME),M,TMAT,M,C0,
     &                 WKM2,M)
C
            HTOT(IT,IME) = HTOT(IT,IME) + WETAB(IECURR,IEPATH)
     &                     *CMATTRC(N,M,WKM2)
C
         END DO
CME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME ME
C
      END DO LOOP_IT
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( ((IEPATH.EQ.NEPATH) .AND. (IECURR.EQ.NETAB(1))) .OR. 
     &     (IPRINT.GT.0) ) THEN
C
         DO IT = 1,NT
C
            IF ( Z(IT).EQ.0 ) CYCLE
C
            DO IME = 1,NMEHFMAX
               HI(IT,IME) = HI(IT,IME) + DEFERMI*H(IT,IME)
            END DO
C
            DO IWRLOG = 1,NWRLOG
               IF ( IWRLOG.EQ.1 ) IW = 6
               IF ( IWRLOG.EQ.2 ) IW = IFILLOG
               IF ( (IEPATH.EQ.NEPATH) .AND. (IECURR.EQ.NETAB(1)) ) THEN
                  WRITE (IW,99002) IT,TXT_T(IT),RNUCTAB(Z(IT))
               ELSE
                  WRITE (IW,99002) IT,TXT_T(IT)
               END IF
C
               BSUM = 0D0
               EVSUM = 0D0
               DO IL = 1,NLT(IT)
                  BSUML = 0D0
                  DO IBSPLIT = 3,5
                     BSUML = BSUML + DIMAG(HINT(IL,IT,IBSPLIT))
                  END DO
                  BSUM = BSUM + BSUML
                  EVSUM = EVSUM + OBS_LT(0,IHFF,IL,IT)
                  IF ( IL.EQ.1 ) THEN
                     WRITE (IW,99003) DIMAG(HINT(IL,IT,2)),
     &                                (DIMAG(HINT(IL,IT,IBSPLIT))/1D3,
     &                                IBSPLIT=3,5),BSUML/1D3,
     &                                (BSUML-OBS_LT(0,IHFF,IL,IT))/1D3
                  ELSE
                     WRITE (IW,99004) TXT_L(IL-1),DIMAG(HINT(IL,IT,1))
     &                                /OBS_LT(0,IDOS,IL,IT),
     &                                DIMAG(HINT(IL,IT,2)),
     &                                (DIMAG(HINT(IL,IT,IBSPLIT))/1D3,
     &                                IBSPLIT=3,5),BSUML/1D3,
     &                                (BSUML-OBS_LT(0,IHFF,IL,IT))/1D3
                  END IF
               END DO
               WRITE (IW,99005) HI(IT,2),
     &                          (HI(IT,IBSPLIT)/1D3,IBSPLIT=3,5),
     &                          BSUM/1D3,(BSUM-EVSUM)/1D3
            END DO
         END DO
      END IF
C
99001 FORMAT (//,1X,79('*'),/,35X,'<HFFIELD>',/,24X,
     &        'hyperfine interaction parameters',/,1X,79('*'),//,10X,
     &        'DEL: deviation of the sum of contributions from the ',/,
     &        10X,
     &        'magnetic hyperfine field obtained from Breit expression',
     &        /)
99002 FORMAT (/,' IT=',I2,2X,A2,:,5X,'r_nuc =',F9.6,' a.u.',/)
99003 FORMAT (12X,' 1/r^3        EFG           ',
     &        'decomposition of hyperfine field B [kG]',/,12X,
     &        '[a.u.]   [10^22 cm^2]',
     &        '       FC       dip      orb      sum      DEL',/,
     &        '  s  ',13X,E15.4,2X,5F16.8)
99004 FORMAT (2X,A,2X,F13.4,E15.4,2X,5F16.8)
99005 FORMAT (' sum ',13X,E15.4,2X,5F16.8)
C99004 FORMAT (2X,A,2X,F13.4,E15.4,2X,5F9.3)
C99005 FORMAT (' sum ',13X,E15.4,2X,5F9.3)
      END
C*==hff_split_radint.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE HFF_SPLIT_RADINT(IM,GG,GA,GB,DR,R,RNUC,NCPLWF1,NCPLWF2,
     &                            IRTOP,NRC)
C   ********************************************************************
C   *                                                                  *
C   *  perform the radial integrals connected with the                 *
C   *  splitting of the relativistic hyperfine field                   *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_TYPES,ONLY:NCPLWFMAX
      IMPLICIT NONE
C*--HFF_SPLIT_RADINT506
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,NCPLWF1,NCPLWF2,NRC
      REAL*8 RNUC
      REAL*8 DR(NRC),R(NRC)
      COMPLEX*16 GA(NRC,NCPLWFMAX),GB(NRC,NCPLWFMAX),
     &           GG(NCPLWFMAX,NCPLWFMAX)
C
C Local variables
C
      COMPLEX*16 CYLAG
      INTEGER I,IR,K1,K2
      REAL*8 X(5)
      COMPLEX*16 Y(5),YI(NRC),ZI(NRC)
C
C*** End of declarations rewritten by SPAG
C
      DO K1 = 1,NCPLWF1
         DO K2 = 1,NCPLWF2
C
            DO IR = 1,IRTOP
               YI(IR) = GA(IR,K1)*GB(IR,K2)*DR(IR)
            END DO
C
            CALL CRADINT_R(IM,YI,ZI)
C
            IF ( ABS(RNUC).GT.1D-10 ) THEN
               DO I = 1,4
                  X(I) = R(IRTOP-4+I)
                  Y(I) = ZI(IRTOP-4+I)
               END DO
               ZI(IRTOP) = CYLAG(RNUC,X,Y,0,3,4)
            END IF
C
            X(1) = 1.0D0
            X(2) = 3.0D0
            X(3) = 5.0D0
            Y(1) = ZI(IRTOP) - ZI(1)
            Y(2) = ZI(IRTOP) - ZI(3)
            Y(3) = ZI(IRTOP) - ZI(5)
C
            GG(K1,K2) = CYLAG(0.0D0,X,Y,0,2,3)
C
         END DO
      END DO
      END
C*==hff_split_sumup.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE HFF_SPLIT_SUMUP(ZHFIZ,LAM1,LAM2,IME,VG,G,VF,F,AHF_F,
     &                           AHF_G,NCPLWF,IKMCPLWF)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NCPLWFMAX
      IMPLICIT NONE
C*--HFF_SPLIT_SUMUP579
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IME,LAM1,LAM2
      REAL*8 VF,VG
      REAL*8 AHF_F(NKMMAX,NKMMAX,5),AHF_G(NKMMAX,NKMMAX,5)
      COMPLEX*16 F(NCPLWFMAX,NCPLWFMAX),G(NCPLWFMAX,NCPLWFMAX),
     &           ZHFIZ(NKMMAX,NKMMAX,5)
      INTEGER IKMCPLWF(NCPLWFMAX,NKMMAX),NCPLWF(NKMMAX)
C
C Local variables
C
      INTEGER IKM1,IKM2,K1,K2
      COMPLEX*16 SUMF,SUMG
C
C*** End of declarations rewritten by SPAG
C
      SUMG = 0.0D0
      SUMF = 0.0D0
      DO K2 = 1,NCPLWF(LAM2)
         IKM2 = IKMCPLWF(K2,LAM2)
         DO K1 = 1,NCPLWF(LAM1)
            IKM1 = IKMCPLWF(K1,LAM1)
C
            SUMG = SUMG + G(K1,K2)*AHF_G(IKM1,IKM2,IME)
            SUMF = SUMF + F(K1,K2)*AHF_F(IKM1,IKM2,IME)
C
         END DO
      END DO
C
      ZHFIZ(LAM1,LAM2,IME) = VG*SUMG + VF*SUMF
C
      END
