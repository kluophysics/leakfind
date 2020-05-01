C*==chiradint.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIRADINT(ZGA,ZFA,JGA,JFA,NOBS,NOBSDNS,NOP,IDOS,ISPN,
     &                     IORB,IHFI,IHVV,ICDIA,IKDIA,IRM3,NSPINOBS,
     &                     NSPINOP,DZL,DZLR,D0ZL,D0XL,D0XR,D1ZL,DIJZL,
     &                     DIJXL,DIJZLTR,DIJXLTR,DIJXR,DOBS,RHORTL,
     &                     TAUTLIN,TKTKTT,DDTAUTAUT,BXCNMM,BXCNNM,
     &                     BXCNMN,BXCNNN,AXCN,CHIPRINT,NTK,NTKMAX,
     &                     NTKTKMAX)
C   ********************************************************************
C   *                                                                  *
C   *  CALCULATION OF TR{ INT(R)INT(R') G(R,R',E)*G(R',R,E) }          *
C   *    NEEDED AS INTEGRAND FOR ENERGY INTEGRATION OF G*G IN CHICALC  *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1      1   IDOS                           *
C   *  spin moment           b s_z  2   ISPN                           *
C   *  orbital moment        b l_z  3   IORB    ____ NOP               *
C   *  hyperfine interaction H_hf,z 4   IHFI                           *
C   *  VanVlecK hyperfine    H_VV,z 5   IHVV    ____ NOBS = NOP+2      *
C   *  diamagnetic suscept   r^2    6   ICDIA                          *
C   *  diamagnetic hyperfine r^-1   7   IKDIA                          *
C   *  hyperfine r-weight    r^-3   8   IRM3    ____ NOBSDNS = NOBS+3  *
C   *                                                                  *
C   *  MD + HF + HE 1996 - 2001                                        *
C   ********************************************************************
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_TYPES,ONLY:NT,NTMAX,IMT,NCPLWFMAX,IKMCPLWF
      USE MOD_RMESH,ONLY:NRMAX,R,DRDI,R2DRDI,JRWS
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CONSTANTS,ONLY:B_AU2CGS
      USE MOD_ANGMOM,ONLY:ISMT,AMEOPO,CGC,NLMAX,NKMMAX,NLQ,AME_G,
     &    IMKM_IKM,NCPLWF
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      IMPLICIT NONE
C*--CHIRADINT36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHIRADINT')
      REAL*8 F1
      PARAMETER (F1=1.0D0)
C
C Dummy arguments
C
      INTEGER CHIPRINT,ICDIA,IDOS,IHFI,IHVV,IKDIA,IORB,IRM3,ISPN,NOBS,
     &        NOBSDNS,NOP,NTK,NTKMAX,NTKTKMAX
      REAL*8 AXCN(NRMAX,2,NLMAX,NTMAX),BXCNMM(NRMAX,NTMAX),
     &       BXCNMN(NRMAX,NTMAX),BXCNNM(NRMAX,NTMAX),BXCNNN(NRMAX,NTMAX)
      COMPLEX*16 D0XL(NLMAX,NTMAX,2,NOBS,2,NOP),
     &           D0XR(NRMAX,NTMAX,2,NOBS,2,NOP),
     &           D0ZL(NLMAX,NTMAX,2,NOBS,2,NOP),
     &           D1ZL(NLMAX,NTMAX,2,NOBS,2,NOP),
     &           DDTAUTAUT(NTKTKMAX,NTMAX),
     &           DIJXL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP),
     &           DIJXLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP),
     &           DIJXR(NRMAX,NTMAX,NTMAX,2,NOBS,2,NOP),
     &           DIJZL(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP),
     &           DIJZLTR(NLMAX,NTMAX,NTMAX,2,NOBS,2,NOP),
     &           DOBS(NLMAX,NTMAX,2,NOBSDNS),
     &           DZL(NLMAX,NTMAX,2,NOBS,2,NOP),
     &           DZLR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP),JFA(NRMAX,2,2),
     &           JGA(NRMAX,2,2),RHORTL(NRMAX,NLMAX,NTMAX,2,NOBS),
     &           TAUTLIN(NTKMAX,NTMAX),TKTKTT(NTKTKMAX,NTMAX,NTMAX),
     &           ZFA(NRMAX,2,2),ZGA(NRMAX,2,2)
      INTEGER NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      REAL*8 AMEF(:,:,:,:),AMEG(:,:,:,:),CFF(2,2),CFG(2,2),CGF(2,2),
     &       CGG(2,2),CH(2,2),COF(2,2,2),COG(2,2,2),MJ,R1M(2,2),WGTF(:),
     &       WGTG(:),WROBS(:,:),XF1,XG1
      COMPLEX*16 AUXDIJXL(:,:,:,:),AUXDIJXLTR(:,:,:,:),AUXDIJXR(:,:,:,:)
     &           ,CDRDI(:),CR1(:),CR2DRDI(:),D1Z(:,:,:,:,:,:),
     &           D1ZR(:,:,:,:,:,:,:),D2X(:,:,:,:,:,:),
     &           D2XR(:,:,:,:,:,:,:),D2Z(:,:,:,:,:,:),
     &           D2ZR(:,:,:,:,:,:,:),D3X(:,:,:,:,:,:),
     &           D3XR(:,:,:,:,:,:,:),D3Z(:,:,:,:,:,:),
     &           D3ZR(:,:,:,:,:,:,:),D4X(:,:,:,:,:,:),
     &           D4XR(:,:,:,:,:,:,:),D4Z(:,:,:,:,:,:),
     &           D4ZR(:,:,:,:,:,:,:),DAAX(:,:,:,:,:),DAAXR(:,:,:,:,:,:),
     &           DAAZ(:,:,:,:),DAAZR(:,:,:,:,:),DBBX(:,:,:,:,:),
     &           DBBXR(:,:,:,:,:,:),DBBZ(:,:,:,:),DBBZR(:,:,:,:,:),
     &           DCDIAZJ,DCDIAZZ,DKDIAZJ,DKDIAZZ,DRM3ZJ,DRM3ZZ,
     &           HAX(:,:,:,:,:,:),HAXT(:,:,:,:,:),HAZ(:,:,:,:,:),
     &           HAZR(:,:,:,:,:),HAZRT(:,:,:,:,:),HAZT(:,:,:,:),
     &           HBX(:,:,:,:,:,:),HBZ(:,:,:,:,:),HBZR(:,:,:,:,:),
     &           HCX(:,:,:,:,:,:),HCZ(:,:,:,:,:),HCZR(:,:,:,:,:),
     &           HDX(:,:,:,:,:,:),HDZ(:,:,:,:,:),HDZR(:,:,:,:,:),
     &           JF(:,:,:),JG(:,:,:),KXR2DRDI(:,:,:,:),
     &           RMEHF(NCPLWFMAX,NCPLWFMAX),TAU12,TAU13,TAU34,TTIJ,
     &           TTIJTR,WFWF(:,:,:),WGTRCDIA(:),WGTRKDIA(:),
     &           WGTROOX(:,:,:,:,:,:),WGTROOZ(:,:,:,:,:),WGTRRM3(:),
     &           WGWG(:,:,:),ZF(:,:,:),ZG(:,:,:),ZJ,ZZ
      LOGICAL CHECKORB,SUPPRESS
      COMPLEX*16 CSUMSSLT
      INTEGER I,IA_ERR,ICHI,IKM1,IKM2,IKMCB(2),IL,IM,IMKM1,IMKM2,IOBS,
     &        IOP,IQ,IR,IRTOP,IS,IS1,IS2,ISPIN,IT,J,JT,K,K1,K2,K3,K4,
     &        KAP1,KAP2,L,LAM,LIN,LIN12,LINKK(2,2),LINTKTK,MJM05,MUEM05,
     &        NCHI(:),NLINL,NSOL
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
      DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
C
      ALLOCATABLE KXR2DRDI,WGTRCDIA,WGTRKDIA,AUXDIJXL,AUXDIJXR,CR1
      ALLOCATABLE D2X,D1Z,D2Z,D3X,D3Z,D4X,D4Z,AUXDIJXLTR,HAX,HAZ,HBX
      ALLOCATABLE HBZ,HCX,HCZ,HDX,HDZ,D2XR,D1ZR,D2ZR,D3XR,D3ZR,D4XR
      ALLOCATABLE D4ZR,AMEF,AMEG,DAAX,DAAZ,DBBX,DBBZ,NCHI,HAXT,HAZR
      ALLOCATABLE WGTF,WGTG,HAZT,HBZR,HCZR,HDZR,WFWF,WGWG,CDRDI
      ALLOCATABLE DAAXR,DAAZR,DBBXR,DBBZR,HAZRT,WROBS
      ALLOCATABLE CR2DRDI,WGTRRM3,WGTROOX,WGTROOZ
      ALLOCATABLE ZF,ZG,JF,JG
C
      ALLOCATE (KXR2DRDI(NRMAX,2,NOP,2),WGTRCDIA(NRMAX))
      ALLOCATE (WGTRKDIA(NRMAX),AUXDIJXL(NLMAX,NOBS,NOP,2))
      ALLOCATE (AUXDIJXR(NRMAX,NOBS,NOP,2),CR1(NRMAX))
      ALLOCATE (D2X(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D1Z(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D2Z(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D3X(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D3Z(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D4X(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D4Z(NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (AUXDIJXLTR(NLMAX,NOBS,NOP,2))
      ALLOCATE (HAX(2,2,NRMAX,2,NOBS,2),HAZ(2,2,NRMAX,2,NOBS))
      ALLOCATE (HBX(2,2,NRMAX,2,NOBS,2),HBZ(2,2,NRMAX,2,NOBS))
      ALLOCATE (HCX(2,2,NRMAX,2,NOBS,2),HCZ(2,2,NRMAX,2,NOBS))
      ALLOCATE (HDX(2,2,NRMAX,2,NOBS,2),HDZ(2,2,NRMAX,2,NOBS))
      ALLOCATE (D2XR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D1ZR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D2ZR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D3XR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D3ZR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D4XR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (D4ZR(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP))
      ALLOCATE (AMEF(2,2,2,NOBS),AMEG(2,2,2,NOBS))
      ALLOCATE (DAAX(2,NOBS,2,NOP,2),DAAZ(2,NOBS,2,NOP))
      ALLOCATE (DBBX(2,NOBS,2,NOP,2),DBBZ(2,NOBS,2,NOP),NCHI(NOBS))
      ALLOCATE (HAXT(NTKMAX,NTMAX,2,NOBS,2),HAZR(2,2,NRMAX,2,NOBS))
      ALLOCATE (WGTF(NOBS),WGTG(NOBS),HAZT(NTKMAX,NTMAX,2,NOBS))
      ALLOCATE (HBZR(2,2,NRMAX,2,NOBS))
      ALLOCATE (HCZR(2,2,NRMAX,2,NOBS),HDZR(2,2,NRMAX,2,NOBS))
      ALLOCATE (WFWF(2,2,NRMAX),WGWG(2,2,NRMAX),CDRDI(NRMAX))
      ALLOCATE (DAAXR(NRMAX,2,NOBS,2,NOP,2))
      ALLOCATE (DAAZR(NRMAX,2,NOBS,2,NOP))
      ALLOCATE (DBBXR(NRMAX,2,NOBS,2,NOP,2))
      ALLOCATE (DBBZR(NRMAX,2,NOBS,2,NOP))
      ALLOCATE (HAZRT(NRMAX,NTKMAX,NTMAX,2,NOBS),WROBS(NRMAX,NOBS))
      ALLOCATE (CR2DRDI(NRMAX))
      ALLOCATE (WGTRRM3(NRMAX),WGTROOX(NRMAX,2,NOBS,2,NOP,2))
      ALLOCATE (WGTROOZ(NRMAX,2,NOBS,2,NOP))
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      CHECKORB = .FALSE.
C
      NCHI(1) = 2
      NCHI(2) = 2
      DO IOP = 3,NOBS
         NCHI(IOP) = 1
      END DO
C
      WGTG(1) = 1D0
      WGTF(1) = 1D0
      WGTG(2) = 1D0
      WGTF(2) = -1D0
      WGTG(3) = 1D0
      WGTF(3) = -1D0
C
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS,RHORTL)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D1Z)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D2Z)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D2X)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D3Z)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D3X)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D4Z)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D4X)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,DZL)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D1ZL)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D0ZL)
      CALL CINIT(NLMAX*NTMAX*2*NOBS*2*NOP,D0XL)
      CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJZL)
      CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXL)
      CALL CINIT(NLMAX*NOBS*NOP*2,AUXDIJXL)
      CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJZLTR)
      CALL CINIT(NLMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXLTR)
      CALL CINIT(NRMAX*NTMAX*NTMAX*2*NOBS*2*NOP,DIJXR)
      CALL CINIT(NLMAX*NOBS*NOP*2,AUXDIJXLTR)
      CALL CINIT(NRMAX*NOBS*NOP*2,AUXDIJXR)
C
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D1ZR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D2ZR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D2XR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D3ZR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D3XR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D4ZR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,D4XR)
      CALL CINIT(NRMAX*NLMAX*NTMAX*2*NOBS*2*NOP,DZLR)
      CALL CINIT(NRMAX*NTMAX*2*NOBS*2*NOP,D0XR)
      CALL CINIT(2*2*NRMAX*2*NOBS,HAZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HBZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HCZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HDZ)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HAX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HBX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HCX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HDX)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
         DO I = 1,IRTOP
            CR1(I) = C1
            CDRDI(I) = DCMPLX(DRDI(I,IM),0.0D0)
            CR2DRDI(I) = DCMPLX(R2DRDI(I,IM),0.0D0)
C
C---------------------------------- terms for 1*(ENN*chi_n + ENM*chi_s)
            ICHI = IDOS
            KXR2DRDI(I,1,IDOS,ICHI) = BXCNNN(I,IT)*CR2DRDI(I)
            KXR2DRDI(I,1,ISPN,ICHI) = BXCNMN(I,IT)*CR2DRDI(I)
C---------------------- terms for  beta*sigma_z*(EMN*chi_n + EMM*chi_s)
            ICHI = ISPN
            KXR2DRDI(I,1,IDOS,ICHI) = BXCNNM(I,IT)*CR2DRDI(I)
            KXR2DRDI(I,1,ISPN,ICHI) = BXCNMM(I,IT)*CR2DRDI(I)
C
         END DO
C
         DO IOBS = 1,NOBS
            DO I = 1,IRTOP
               WROBS(I,IOBS) = R2DRDI(I,IM)
            END DO
         END DO
C
         CALL DCOPY(IRTOP,DRDI(1,IM),1,WROBS(1,IHFI),1)
C
         DO I = 1,IRTOP
            WROBS(I,IHVV) = DRDI(I,IM)/R(I,IM)
         END DO
C
         DO IOBS = 1,NOBSDNS
            DO IS1 = 1,2
               DO IL = 1,NLMAX
                  DOBS(IL,IT,IS1,IOBS) = C0
               END DO
            END DO
         END DO
C
         LIN = 0
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,(NLQ(IQ)-1)
C
            IL = L + 1
C ------------------------------------------- AXCN <> 0 only for L = LOP
            DO IS = 1,2
               DO I = 1,IRTOP
                  KXR2DRDI(I,IS,IORB,1) = -AXCN(I,IS,IL,IT)*CR2DRDI(I)
               END DO
            END DO
C
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MJM05 = -L - 1, + L
               MJ = DBLE(MJM05) + 0.5D0
C
               MUEM05 = NINT(MJ-0.5D0)
C
               IF ( ABS(MJ).GE.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C-----------------------------------------------------------------------
               IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
               IMKM1 = IMKM_IKM(IKM1)
               IMKM2 = IMKM_IKM(IKM2)
               IKMCB(1) = IKM1
               IKMCB(2) = IKM2
C-----------------------------------------------------------------------
C
C=======================================================================
C                copy wave functions to OLD 2x2 indexing scheme
C=======================================================================
               DO K = 1,NSOL
                  LAM = IKMCB(K)
                  DO I = 1,NSOL
                     ZGA(:,I,K) = ZG(:,I,LAM)
                     ZFA(:,I,K) = ZF(:,I,LAM)
                     JGA(:,I,K) = JG(:,I,LAM)
                     JFA(:,I,K) = JF(:,I,LAM)
                  END DO
               END DO
C=======================================================================
C
C
C=======================================================================
C===============================================================  LINKK
C
               DO K1 = 1,NSOL
                  DO K2 = 1,NSOL
                     LIN = LIN + 1
                     LINKK(K1,K2) = LIN
                  END DO
               END DO
C
C=======================================================================
C
C   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
C
               CGG(1,1) = AME_G(IKM1,IKM1,2,ISMT)
               CGG(1,2) = AME_G(IKM1,IKM2,2,ISMT)
               CGG(2,1) = AME_G(IKM2,IKM1,2,ISMT)
               CGG(2,2) = AME_G(IKM2,IKM2,2,ISMT)
               CALL RINIT(4,CGF)
               CGF(1,1) = AME_G(IMKM1,IMKM1,2,ISMT)
               CGF(2,2) = AME_G(IMKM2,IMKM2,2,ISMT)
C
C-----------------------------------------------------------------------
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C-----------------------------------------------------------------------
C   COEFFICIENTS TO CALCULATE THE ORBITAL POLARISATION
C
               CALL RINIT(4*2,COF)
               DO IS = 1,2
                  COG(1,1,IS) = AMEOPO(IKM1,IKM1,IS)
                  COG(2,2,IS) = AMEOPO(IKM2,IKM2,IS)
                  COG(1,2,IS) = AMEOPO(IKM1,IKM2,IS)
                  COG(2,1,IS) = AMEOPO(IKM2,IKM1,IS)
C
                  COF(1,1,IS) = CGC(IMKM1,IS)*CGC(IMKM1,IS)
     &                          *DBLE(MUEM05-IS+2)
                  IF ( NSOL.EQ.2 ) COF(2,2,IS) = CGC(IMKM2,IS)
     &                 *CGC(IMKM2,IS)*DBLE(MUEM05-IS+2)
               END DO
C
               DO I = 1,NSOL
                  DO J = 1,NSOL
                     AMEG(I,J,1,IDOS) = R1M(I,J)
                     AMEF(I,J,1,IDOS) = R1M(I,J)
                     AMEG(I,J,1,ISPN) = CGG(I,J)
                     AMEF(I,J,1,ISPN) = CGF(I,J)
                     DO IS = 1,2
                        AMEG(I,J,IS,IORB) = COG(I,J,IS)
                        AMEF(I,J,IS,IORB) = COF(I,J,IS)
                     END DO
                     XG1 = COG(I,J,1) + COG(I,J,2) - CFG(I,J)
                     XF1 = COF(I,J,1) + COF(I,J,2) - CFF(I,J)
                     IF ( ABS(XG1).GT.1D-14 .OR. ABS(XF1).GT.1D-14 )
     &                    THEN
                        WRITE (6,*) '###',I,J,CFG(I,J),COG(I,J,1)
     &                              + COG(I,J,2) - CFG(I,J)
                        WRITE (6,*) '***',I,J,CFF(I,J),COF(I,J,1)
     &                              + COF(I,J,2) - CFF(I,J)
                     END IF
                  END DO
               END DO
C
C-----------------------------------------------------------------------
C   use FULL ORBITAL angular matrix elements for CHI(SPIN) calculations
C
               IF ( CHECKORB ) THEN
                  CALL DCOPY(4,CFG,1,CGG,1)
                  CALL DCOPY(4,CFF,1,CGF,1)
               END IF
C
C-----------------------------------------------------------------------
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C
C-----------------------------------------------------------------------
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
               CH(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
               CH(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
               IF ( NSOL.EQ.2 ) THEN
                  CH(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
                  CH(2,1) = CH(1,2)
               END IF
C=======================================================================
C
               CALL CHITERM1(HAZR,HBZR,HCZR,HDZR,HAZ,HBZ,HCZ,HDZ,HAX,
     &                       HBX,HCX,HDX,ZGA,JGA,WGTG,WGWG,AMEG,ZFA,JFA,
     &                       WGTF,WFWF,AMEF,CR1,CR2DRDI,KXR2DRDI,NSOL,
     &                       IRTOP,NSPINOBS,NOBS,NOP,NRMAX,NCHI)
C
C----------------------------------------------------------- H_hf Zeeman
C
               CALL CHIINTHFF3(IM,HAZ(1,1,IRTOP,1,IHFI),ZGA,ZGA,ZFA,ZFA,
     &                         B_AU2CGS,RMEHF,CH,CDRDI,NSOL,IRTOP)
C
               CALL CHIINTGHFF3(HAZR(1,1,1,1,IHFI),ZGA,ZGA,ZFA,ZFA,
     &                          B_AU2CGS,CH,CR1,NSOL,IRTOP,NRMAX)
C
               CALL CHIINTGHFF3(HBZR(1,1,1,1,IHFI),JGA,ZGA,JFA,ZFA,
     &                          B_AU2CGS,CH,CR1,NSOL,IRTOP,NRMAX)
               CALL CHIINTGHFF3(HCZR(1,1,1,1,IHFI),ZGA,JGA,ZFA,JFA,
     &                          B_AU2CGS,CH,CR1,NSOL,IRTOP,NRMAX)
               CALL CHIINTGHFF3(HDZR(1,1,1,1,IHFI),JGA,JGA,JFA,JFA,
     &                          B_AU2CGS,CH,CR1,NSOL,IRTOP,NRMAX)
C
C------------------------------------------------------- H_l r^-3 Zeeman
               DO I = 1,IRTOP
                  WGTRRM3(I) = CDRDI(I)/R(I,IM)
               END DO
C
               CALL CHIINTHVV3(IM,HAZ(1,1,1,1,IHVV),ZGA,ZGA,F1,CFG,ZFA,
     &                         ZFA,-F1,CFF,WGTRRM3,NSOL,IRTOP)
C
C ===============================================  RAD. INT. FOR 1.TERM:
C
               DO K1 = 1,NSOL
                  DO K2 = 1,NSOL
                     LIN12 = LINKK(K1,K2)
C
                     DO IOBS = 1,NOBS
                        DO ISPIN = 1,NSPINOBS(IOBS)
                           HAZT(LIN12,IT,ISPIN,IOBS)
     &                        = HAZ(K1,K2,IRTOP,ISPIN,IOBS)
                           DO I = 1,IRTOP
                              HAZRT(I,LIN12,IT,ISPIN,IOBS)
     &                           = HAZR(K1,K2,I,ISPIN,IOBS)
                           END DO
                        END DO
                     END DO
C
                     DO IOP = 1,NOP
                        DO ISPIN = 1,NSPINOP(IOP)
                           HAXT(LIN12,IT,ISPIN,IOP,1)
     &                        = HAX(K1,K2,IRTOP,ISPIN,IOP,1)
                           HAXT(LIN12,IT,ISPIN,IOP,2)
     &                        = HAX(K1,K2,IRTOP,ISPIN,IOP,2)
C
                        END DO
                     END DO
C
                  END DO
               END DO
C
C =============================================================  2.TERM:
               DO K1 = 1,NSOL
                  DO K3 = 1,NSOL
                     DO K4 = 1,NSOL
C
C--------------------------------------------------------------- part AA
C
                        CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HAZ,HAX,
     &                               'out',K1,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTPREP(DAAZR,DAAXR,HCZR,HAZ,HAX,'out',
     &                                  K4,K1,K1,K3,IRTOP,NSPINOBS,NOBS,
     &                                  NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINT(DAAZ,DAAX,ZGA(1,1,K4),JGA(1,1,K1),
     &                              WGTG,WGWG,AMEG,ZFA(1,1,K4),
     &                              JFA(1,1,K1),WGTF,WFWF,AMEF,WGTROOZ,
     &                              WGTROOX,NSOL,IRTOP,NSPINOBS,NOBS,
     &                              NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTHF(IM,DAAZ,DAAX,ZGA(1,1,K4),
     &                                JGA(1,1,K1),ZFA(1,1,K4),
     &                                JFA(1,1,K1),B_AU2CGS,RMEHF,CH,
     &                                WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                                NSPINOBS,NOBS,NSPINOP,NOP,IHFI,
     &                                NCHI)
C
                        DO IS1 = 1,2
                           CALL CHIINTHVV4(IM,DAAZ(1,IHVV,IS1,IORB),
     &                        ZGA(1,1,K4),JGA(1,1,K1),F1,CFG,ZFA(1,1,K4)
     &                        ,JFA(1,1,K1),-F1,CFF,
     &                        WGTROOZ(1,1,IHVV,IS1,IORB),NSOL,IRTOP)
                        END DO
C
C--------------------------------------------------------------- part BB
C
                        CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HBZ,HBX,
     &                               'in ',K1,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTPREP(DBBZR,DBBXR,HAZR,HBZ,HBX,'in ',
     &                                  K4,K1,K1,K3,IRTOP,NSPINOBS,NOBS,
     &                                  NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINT(DBBZ,DBBX,ZGA(1,1,K4),ZGA(1,1,K1),
     &                              WGTG,WGWG,AMEG,ZFA(1,1,K4),
     &                              ZFA(1,1,K1),WGTF,WFWF,AMEF,WGTROOZ,
     &                              WGTROOX,NSOL,IRTOP,NSPINOBS,NOBS,
     &                              NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTHF(IM,DBBZ,DBBX,ZGA(1,1,K4),
     &                                ZGA(1,1,K1),ZFA(1,1,K4),
     &                                ZFA(1,1,K1),B_AU2CGS,RMEHF,CH,
     &                                WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                                NSPINOBS,NOBS,NSPINOP,NOP,IHFI,
     &                                NCHI)
C
                        DO IS1 = 1,2
                           CALL CHIINTHVV4(IM,DBBZ(1,IHVV,IS1,IORB),
     &                        ZGA(1,1,K4),ZGA(1,1,K1),F1,CFG,ZFA(1,1,K4)
     &                        ,ZFA(1,1,K1),-F1,CFF,
     &                        WGTROOZ(1,1,IHVV,IS1,IORB),NSOL,IRTOP)
                        END DO
C
C---------------------------------------------------- sum part AA and BB
C
                        TAU34 = TAUTLIN(LINKK(K3,K4),IT)
C
                        CALL CHIDSUM(D2Z,D2X,D2ZR,D2XR,DAAZ,DBBZ,DAAX,
     &                               DBBX,DAAZR,DBBZR,DAAXR,DBBXR,TAU34,
     &                               IL,IT,IRTOP,NLMAX,NTMAX,NSPINOBS,
     &                               NOBS,NSPINOP,NOP,NRMAX,NCHI)
C
                     END DO
                  END DO
               END DO
C
C =============================================================  3.TERM:
               DO K3 = 1,NSOL
                  DO K1 = 1,NSOL
                     DO K2 = 1,NSOL
C
C--------------------------------------------------------------- part AA
C
                        CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HAZ,HAX,
     &                               'out',K2,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTPREP(DAAZR,DAAXR,HBZR,HAZ,HAX,'out',
     &                                  K3,K1,K2,K3,IRTOP,NSPINOBS,NOBS,
     &                                  NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINT(DAAZ,DAAX,JGA(1,1,K3),ZGA(1,1,K1),
     &                              WGTG,WGWG,AMEG,JFA(1,1,K3),
     &                              ZFA(1,1,K1),WGTF,WFWF,AMEF,WGTROOZ,
     &                              WGTROOX,NSOL,IRTOP,NSPINOBS,NOBS,
     &                              NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTHF(IM,DAAZ,DAAX,JGA(1,1,K3),
     &                                ZGA(1,1,K1),JFA(1,1,K3),
     &                                ZFA(1,1,K1),B_AU2CGS,RMEHF,CH,
     &                                WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                                NSPINOBS,NOBS,NSPINOP,NOP,IHFI,
     &                                NCHI)
C
                        DO IS2 = 1,2
                           CALL CHIINTHVV4(IM,DAAZ(1,IHVV,IS2,IORB),
     &                        JGA(1,1,K3),ZGA(1,1,K1),F1,CFG,JFA(1,1,K3)
     &                        ,ZFA(1,1,K1),-F1,CFF,
     &                        WGTROOZ(1,1,IHVV,IS2,IORB),NSOL,IRTOP)
                        END DO
C
C--------------------------------------------------------------- part BB
C
                        CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HCZ,HCX,
     &                               'in ',K2,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTPREP(DBBZR,DBBXR,HAZR,HCZ,HCX,'in ',
     &                                  K3,K1,K2,K3,IRTOP,NSPINOBS,NOBS,
     &                                  NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINT(DBBZ,DBBX,ZGA(1,1,K3),ZGA(1,1,K1),
     &                              WGTG,WGWG,AMEG,ZFA(1,1,K3),
     &                              ZFA(1,1,K1),WGTF,WFWF,AMEF,WGTROOZ,
     &                              WGTROOX,NSOL,IRTOP,NSPINOBS,NOBS,
     &                              NSPINOP,NOP,NRMAX,NCHI)
C
                        CALL CHIINTHF(IM,DBBZ,DBBX,ZGA(1,1,K3),
     &                                ZGA(1,1,K1),ZFA(1,1,K3),
     &                                ZFA(1,1,K1),B_AU2CGS,RMEHF,CH,
     &                                WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                                NSPINOBS,NOBS,NSPINOP,NOP,IHFI,
     &                                NCHI)
C
                        DO IS2 = 1,2
                           CALL CHIINTHVV4(IM,DBBZ(1,IHVV,IS2,IORB),
     &                        ZGA(1,1,K3),ZGA(1,1,K1),F1,CFG,ZFA(1,1,K3)
     &                        ,ZFA(1,1,K1),-F1,CFF,
     &                        WGTROOZ(1,1,IHVV,IS2,IORB),NSOL,IRTOP)
                        END DO
C
C---------------------------------------------------- sum part AA and BB
C
                        TAU12 = TAUTLIN(LINKK(K1,K2),IT)
C
                        CALL CHIDSUM(D3Z,D3X,D3ZR,D3XR,DAAZ,DBBZ,DAAX,
     &                               DBBX,DAAZR,DBBZR,DAAXR,DBBXR,TAU12,
     &                               IL,IT,IRTOP,NLMAX,NTMAX,NSPINOBS,
     &                               NOBS,NSPINOP,NOP,NRMAX,NCHI)
                     END DO
                  END DO
               END DO
C
C =============================================================  4.TERM:
C
               DO K1 = 1,NSOL
                  DO K3 = 1,NSOL
C
C--------------------------------------------------------------- part AA
C
                     CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HAZ,HAX,'out',
     &                            K1,K3,IRTOP,NSPINOBS,NOBS,NSPINOP,NOP,
     &                            NRMAX,NCHI)
C
                     CALL CHIINTPREP(DAAZR,DAAXR,HDZR,HAZ,HAX,'out',K1,
     &                               K3,K1,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                     CALL CHIINT(DAAZ,DAAX,JGA(1,1,K1),JGA(1,1,K3),WGTG,
     &                           WGWG,AMEG,JFA(1,1,K1),JFA(1,1,K3),WGTF,
     &                           WFWF,AMEF,WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                           NSPINOBS,NOBS,NSPINOP,NOP,NRMAX,NCHI)
C
                     CALL CHIINTHF(IM,DAAZ,DAAX,JGA(1,1,K1),JGA(1,1,K3),
     &                             JFA(1,1,K1),JFA(1,1,K3),B_AU2CGS,
     &                             RMEHF,CH,WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                             NSPINOBS,NOBS,NSPINOP,NOP,IHFI,NCHI)
C
                     DO IS1 = 1,2
                        DO IS2 = 1,2
                           CALL CHIINTABR4(DAAZ(IS1,IORB,IS2,IORB),
     &                        JGA(1,1,K1),JGA(1,1,K3),F1,WGWG,
     &                        COG(1,1,IS1),JFA(1,1,K1),JFA(1,1,K3),-F1,
     &                        WFWF,COF(1,1,IS1),
     &                        WGTROOZ(1,IS1,IORB,IS2,IORB),NSOL,IRTOP,
     &                        NRMAX)
                        END DO
                     END DO
C
                     DO IS2 = 1,2
                        CALL CHIINTHVV4(IM,DAAZ(1,IHVV,IS2,IORB),
     &                                  JGA(1,1,K1),JGA(1,1,K3),F1,CFG,
     &                                  JFA(1,1,K1),JFA(1,1,K3),-F1,CFF,
     &                                  WGTROOZ(1,1,IHVV,IS2,IORB),NSOL,
     &                                  IRTOP)
                     END DO
C
C--------------------------------------------------------------- part BB
C
                     CALL CHIWGTR(WGTROOZ,WGTROOX,WROBS,HDZ,HDX,'out',
     &                            K1,K3,IRTOP,NSPINOBS,NOBS,NSPINOP,NOP,
     &                            NRMAX,NCHI)
C
                     CALL CHIINTPREP(DBBZR,DBBXR,HAZR,HDZ,HDX,'out',K1,
     &                               K3,K1,K3,IRTOP,NSPINOBS,NOBS,
     &                               NSPINOP,NOP,NRMAX,NCHI)
C
                     CALL CHIINT(DBBZ,DBBX,ZGA(1,1,K1),ZGA(1,1,K3),WGTG,
     &                           WGWG,AMEG,ZFA(1,1,K1),ZFA(1,1,K3),WGTF,
     &                           WFWF,AMEF,WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                           NSPINOBS,NOBS,NSPINOP,NOP,NRMAX,NCHI)
C
                     CALL CHIINTHF(IM,DBBZ,DBBX,ZGA(1,1,K1),ZGA(1,1,K3),
     &                             ZFA(1,1,K1),ZFA(1,1,K3),B_AU2CGS,
     &                             RMEHF,CH,WGTROOZ,WGTROOX,NSOL,IRTOP,
     &                             NSPINOBS,NOBS,NSPINOP,NOP,IHFI,NCHI)
C
                     DO IS2 = 1,2
                        CALL CHIINTHVV4(IM,DBBZ(1,IHVV,IS2,IORB),
     &                                  ZGA(1,1,K1),ZGA(1,1,K3),F1,CFG,
     &                                  ZFA(1,1,K1),ZFA(1,1,K3),-F1,CFF,
     &                                  WGTROOZ(1,1,IHVV,IS2,IORB),NSOL,
     &                                  IRTOP)
                     END DO
C
C---------------------------------------------------- sum part AA and BB
C
                     CALL CHIDSUM(D4Z,D4X,D4ZR,D4XR,DAAZ,DBBZ,DAAX,DBBX,
     &                            DAAZR,DBBZR,DAAXR,DBBXR,C1,IL,IT,
     &                            IRTOP,NLMAX,NTMAX,NSPINOBS,NOBS,
     &                            NSPINOP,NOP,NRMAX,NCHI)
C
C-----------------------------------------------------------------------
C         calculate  CHI_dia, K_dia and some auxilary expectation values
C
                     DO I = 1,IRTOP
                        WGTRCDIA(I) = CR2DRDI(I)*R(I,IM)*R(I,IM)
                        WGTRKDIA(I) = CDRDI(I)*R(I,IM)
                     END DO
C
                     CALL CHIINTABR4(DCDIAZZ,ZGA(1,1,K1),ZGA(1,1,K3),F1,
     &                               WGWG,R1M,ZFA(1,1,K1),ZFA(1,1,K3),
     &                               F1,WFWF,R1M,WGTRCDIA,NSOL,IRTOP,
     &                               NRMAX)
                     CALL CHIINTABR4(DKDIAZZ,ZGA(1,1,K1),ZGA(1,1,K3),F1,
     &                               WGWG,R1M,ZFA(1,1,K1),ZFA(1,1,K3),
     &                               F1,WFWF,R1M,WGTRKDIA,NSOL,IRTOP,
     &                               NRMAX)
                     CALL CHIINTHVV4(IM,DRM3ZZ,ZGA(1,1,K1),ZGA(1,1,K3),
     &                               F1,R1M,ZFA(1,1,K1),ZFA(1,1,K3),F1,
     &                               R1M,WGTRRM3,NSOL,IRTOP)
C
                     IF ( K1.EQ.K3 ) THEN
                        CALL CHIINTABR4(DCDIAZJ,ZGA(1,1,K1),JGA(1,1,K3),
     &                                  F1,WGWG,R1M,ZFA(1,1,K1),
     &                                  JFA(1,1,K3),F1,WFWF,R1M,
     &                                  WGTRCDIA,NSOL,IRTOP,NRMAX)
                        CALL CHIINTABR4(DKDIAZJ,ZGA(1,1,K1),JGA(1,1,K3),
     &                                  F1,WGWG,R1M,ZFA(1,1,K1),
     &                                  JFA(1,1,K3),F1,WFWF,R1M,
     &                                  WGTRKDIA,NSOL,IRTOP,NRMAX)
                        CALL CHIINTHVV4(IM,DRM3ZJ,ZGA(1,1,K1),
     &                                  JGA(1,1,K3),F1,R1M,ZFA(1,1,K1),
     &                                  JFA(1,1,K3),F1,R1M,WGTRRM3,NSOL,
     &                                  IRTOP)
                     ELSE
                        DCDIAZJ = C0
                        DKDIAZJ = C0
                        DRM3ZJ = C0
                     END IF
C
                     TAU13 = TAUTLIN(LINKK(K1,K3),IT)
C
                     DO IOBS = 1,NOBS
                        DO ISPIN = 1,NSPINOBS(IOBS)
                           ZZ = HAZ(K3,K1,IRTOP,ISPIN,IOBS)
                           IF ( K1.EQ.K3 ) THEN
                              ZJ = HCZ(K3,K1,IRTOP,ISPIN,IOBS)
                           ELSE
                              ZJ = 0D0
                           END IF
                           DOBS(IL,IT,ISPIN,IOBS)
     &                        = DOBS(IL,IT,ISPIN,IOBS) - (TAU13*ZZ-ZJ)
     &                        /PI
                        END DO
                     END DO
C
                     DOBS(IL,IT,1,ICDIA) = DOBS(IL,IT,1,ICDIA)
     &                  + (TAU13*DCDIAZZ-DCDIAZJ)/PI
C
                     DOBS(IL,IT,1,IKDIA) = DOBS(IL,IT,1,IKDIA)
     &                  + (TAU13*DKDIAZZ-DKDIAZJ)/PI
C
                     DOBS(IL,IT,1,IRM3) = DOBS(IL,IT,1,IRM3)
     &                  - (TAU13*DRM3ZZ-DRM3ZJ)/PI
C
                     IF ( K1.EQ.K3 ) THEN
C
                        DO IOBS = 1,NOBS
                           DO ISPIN = 1,NSPINOBS(IOBS)
                              DO I = 1,IRTOP
                                 ZZ = HAZR(K3,K1,I,ISPIN,IOBS)
                                 ZJ = HCZR(K3,K1,I,ISPIN,IOBS)
                                 RHORTL(I,IL,IT,ISPIN,IOBS)
     &                              = RHORTL(I,IL,IT,ISPIN,IOBS)
     &                              - (TAU13*ZZ-ZJ)/PI
                              END DO
                           END DO
                        END DO
C
                     ELSE
C
                        DO IOBS = 1,NOBS
                           DO ISPIN = 1,NSPINOBS(IOBS)
                              DO I = 1,IRTOP
                                 ZZ = HAZR(K3,K1,I,ISPIN,IOBS)
                                 ZJ = HCZR(K3,K1,I,ISPIN,IOBS)
                                 RHORTL(I,IL,IT,ISPIN,IOBS)
     &                              = RHORTL(I,IL,IT,ISPIN,IOBS)
     &                              - TAU13*ZZ/PI
                              END DO
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
C ======================================================================
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IQ = IQAT(1,IT)
C =============================================================  1.TERM:
C --------- JT ---------------------------------------------------------
         DO JT = 1,NT
C
            SUPPRESS = .FALSE.
C *** interaction of metal with vacuum layers suppressed
C *** if next two lines uncommented
C      IF ( ( Z(IT).GT.0 ) .AND. ( Z(JT).EQ.0 ) )
C     &         SUPPRESS=.TRUE.
C
            IF ( .NOT.SUPPRESS ) THEN
               IL = 1
               NLINL = 2
               LINTKTK = 0
               IF ( LIN.NE.NTK ) THEN
                  WRITE (6,*) 'LIN=',LIN,' <> NTK=',NTK
                  STOP
               END IF
               DO I = 1,LIN
                  DO J = 1,LIN
                     LINTKTK = LINTKTK + 1
C
                     IF ( IT.EQ.JT ) THEN
                        TTIJ = TKTKTT(LINTKTK,IT,JT)
     &                         + DDTAUTAUT(LINTKTK,IT)
                     ELSE
                        TTIJ = TKTKTT(LINTKTK,IT,JT)
                     END IF
                     TTIJTR = TKTKTT(LINTKTK,IT,JT)
C
                     DO IOP = 1,NOP
                        DO IS2 = 1,NSPINOP(IOP)
                           DO IOBS = 1,NOBS
                              DO IS1 = 1,NSPINOBS(IOBS)
C
                                 D1Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                              = D1Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                              + TTIJ*HAZT(I,IT,IS1,IOBS)
     &                              *HAZT(J,JT,IS2,IOP)
C
                                 DIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                              = DIJZL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                              + TTIJ*HAZT(I,IT,IS1,IOBS)
     &                              *HAZT(J,JT,IS2,IOP)
C
                                 DIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                              = DIJZLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                              + TTIJTR*HAZT(I,IT,IS1,IOBS)
     &                              *HAZT(J,JT,IS2,IOP)
C
                                 DO ICHI = 1,NCHI(IOP)
                                    DIJXL(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                                 = DIJXL(IL,IT,JT,IS1,IOBS,IS2,
     &                                 IOP) + TTIJ*HAZT(I,IT,IS1,IOBS)
     &                                 *HAXT(J,JT,IS2,IOP,ICHI)
C
                                    DIJXLTR(IL,IT,JT,IS1,IOBS,IS2,IOP)
     &                                 = DIJXLTR(IL,IT,JT,IS1,IOBS,IS2,
     &                                 IOP) + TTIJTR*HAZT(I,IT,IS1,IOBS)
     &                                 *HAXT(J,JT,IS2,IOP,ICHI)
C
                                 END DO
C
                                 DO IR = 1,IRTOP
C
                                    D1ZR(IR,IL,IT,IS1,IOBS,IS2,IOP)
     &                                 = D1ZR(IR,IL,IT,IS1,IOBS,IS2,IOP)
     &                                 + TTIJ*(HAZRT(IR,I,IT,IS1,IOBS)*
     &                                 HAZT(J,JT,IS2,IOP))
C
                                    DO ICHI = 1,NCHI(IOP)
                                       DIJXR(IR,IT,JT,IS1,IOBS,IS2,IOP)
     &                                    = DIJXR(IR,IT,JT,IS1,IOBS,IS2,
     &                                    IOP)
     &                                    + TTIJ*(HAZRT(IR,I,IT,IS1,
     &                                    IOBS)*HAXT(J,JT,IS2,IOP,ICHI))
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
C   ================================================================== J
                  IF ( I.EQ.NLINL ) THEN
                     IL = IL + 1
                     NLINL = NLINL + 8*(IL-1) + 2
                  END IF
               END DO
C   ================================================================== I
            END IF
C --- .NOT. SUPPRESS ---------------------------------------------------
         END DO
C   ================================================================= JT
C
C =============================================================== SUM UP
         DO IL = 1,NLQ(IQ)
C
            DO IOP = 1,NOP
               DO IS2 = 1,NSPINOP(IOP)
                  DO IOBS = 1,NOBS
                     DO IS1 = 1,NSPINOBS(IOBS)
C
                        DZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = D1Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                     - D2Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                     - D3Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + D4Z(IL,IT,IS1,IOBS,IS2,IOP)
C
                        D1ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = D1Z(IL,IT,IS1,IOBS,IS2,IOP)
C
                        D0ZL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = -D2Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                     - D3Z(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + D4Z(IL,IT,IS1,IOBS,IS2,IOP)
C
                        D0XL(IL,IT,IS1,IOBS,IS2,IOP)
     &                     = -D2X(IL,IT,IS1,IOBS,IS2,IOP)
     &                     - D3X(IL,IT,IS1,IOBS,IS2,IOP)
     &                     + D4X(IL,IT,IS1,IOBS,IS2,IOP)
C
                        DO J = 1,IRTOP
                           DZLR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        = D1ZR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        - D2ZR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        - D3ZR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        + D4ZR(J,IL,IT,IS1,IOBS,IS2,IOP)
C
                           D0XR(J,IT,IS1,IOBS,IS2,IOP)
     &                        = D0XR(J,IT,IS1,IOBS,IS2,IOP)
     &                        - D2XR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        - D3XR(J,IL,IT,IS1,IOBS,IS2,IOP)
     &                        + D4XR(J,IL,IT,IS1,IOBS,IS2,IOP)
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  the DHFLXL(IL,IT,JT) are determined like
C            DIJXL(IL,IT,IT,1,ihfi,IS,iorb) because this
C            is two-sites matrix
C
            DIJXL(IL,IT,IT,1,IHFI,1,ISPN)
     &         = DIJXL(IL,IT,IT,1,IHFI,1,ISPN)
     &         + D0XL(IL,IT,1,IHFI,1,ISPN)
C
            DO IS = 1,2
               DIJXL(IL,IT,IT,1,IHFI,IS,IORB)
     &            = DIJXL(IL,IT,IT,1,IHFI,IS,IORB)
     &            + D0XL(IL,IT,1,IHFI,IS,IORB)
            END DO
C
C ======================================================================
            IF ( CHIPRINT.GE.2 ) THEN
               WRITE (6,'(1X,79(''-''))')
               WRITE (6,*) '--- GG ---'
               WRITE (6,99001) IL,D1Z(IL,IT,1,IDOS,1,IDOS),
     &                         D2Z(IL,IT,1,IDOS,1,IDOS),
     &                         D3Z(IL,IT,1,IDOS,1,IDOS),
     &                         D4Z(IL,IT,1,IDOS,1,IDOS)
               WRITE (6,*) '--- GGTL ---'
               WRITE (6,*) 'IL=',IL,' GGTL=',DZL(IL,IT,1,IDOS,1,IDOS)
               WRITE (6,*) '--- SSZ ---'
               WRITE (6,99001) IL,D1Z(IL,IT,1,ISPN,1,ISPN),
     &                         D2Z(IL,IT,1,ISPN,1,ISPN),
     &                         D3Z(IL,IT,1,ISPN,1,ISPN),
     &                         D4Z(IL,IT,1,ISPN,1,ISPN)
               WRITE (6,*) '--- DSSZL ---'
               WRITE (6,*) 'IL=',IL,' DSSZL=',DZL(IL,IT,1,ISPN,1,ISPN)
               WRITE (6,*) '--- SSX ---'
               WRITE (6,99001) IL,0D0,0D0,D2X(IL,IT,1,ISPN,1,ISPN),
     &                         D3X(IL,IT,1,ISPN,1,ISPN),
     &                         D4X(IL,IT,1,ISPN,1,ISPN)
               WRITE (6,*) '--- LLZ ---'
               WRITE (6,99001) IL,
     &                         CSUMSSLT(D1Z,IL,IT,IORB,IORB,NLMAX,NTMAX,
     &                         2,NOBS,2,NOP),
     &                         CSUMSSLT(D2Z,IL,IT,IORB,IORB,NLMAX,NTMAX,
     &                         2,NOBS,2,NOP),
     &                         CSUMSSLT(D3Z,IL,IT,IORB,IORB,NLMAX,NTMAX,
     &                         2,NOBS,2,NOP),
     &                         CSUMSSLT(D4Z,IL,IT,IORB,IORB,NLMAX,NTMAX,
     &                         2,NOBS,2,NOP)
               WRITE (6,*) '--- D0SSXL ---'
               WRITE (6,*) 'IL=',IL,'D0SSXL1=',D0XL(IL,IT,1,ISPN,1,ISPN)
               WRITE (6,*) 'IL=',IL,'D0SSXL2=',D0XL(IL,IT,1,ISPN,1,ISPN)
               WRITE (6,*) '--- DLLZL ---'
               WRITE (6,*) 'IL=',IL,' DLLZL=',
     &                     CSUMSSLT(DZL,IL,IT,IORB,IORB,NLMAX,NTMAX,2,
     &                     NOBS,2,NOP)
               WRITE (6,*) '--- DCDIAL ---'
               WRITE (6,*) 'IL=',IL,' DCDIAL=',DOBS(IL,IT,1,ICDIA)*PI
               WRITE (6,*) '--- SLZ ---'
               WRITE (6,99001) IL,(D1Z(IL,IT,1,ISPN,1,IORB)+D1Z(IL,IT,1,
     &                         ISPN,2,IORB)),
     &                         (D2Z(IL,IT,1,ISPN,1,IORB)+D2Z(IL,IT,1,
     &                         ISPN,2,IORB)),
     &                         (D3Z(IL,IT,1,ISPN,1,IORB)+D3Z(IL,IT,1,
     &                         ISPN,2,IORB)),
     &                         (D4Z(IL,IT,1,ISPN,1,IORB)+D4Z(IL,IT,1,
     &                         ISPN,2,IORB))
               WRITE (6,*) '--- DSLZL ---'
               WRITE (6,*) 'IL=',IL,' DSLZL=',DZL(IL,IT,1,ISPN,1,IORB)
     &                     + DZL(IL,IT,1,ISPN,2,IORB)
               WRITE (6,*) '--- SLX ---'
               WRITE (6,99001) IL,0D0,0D0,
     &                         (D2X(IL,IT,1,ISPN,1,IORB)+D2X(IL,IT,1,
     &                         ISPN,2,IORB)),
     &                         (D3X(IL,IT,1,ISPN,1,IORB)+D3X(IL,IT,1,
     &                         ISPN,2,IORB)),
     &                         (D4X(IL,IT,1,ISPN,1,IORB)+D4X(IL,IT,1,
     &                         ISPN,2,IORB))
               WRITE (6,*) '--- LSZ ---'
               WRITE (6,99002) (IL,D1Z(IL,IT,IS1,IORB,1,ISPN),D2Z(IL,IT,
     &                         IS1,IORB,1,ISPN),IS1,
     &                         D3Z(IL,IT,IS1,IORB,1,ISPN),
     &                         D4Z(IL,IT,IS1,IORB,1,ISPN),IS1=1,2)
               WRITE (6,*) '--- DLSZSL ---'
               WRITE (6,99003) (IL,IS1,'DLSZSL',DZL(IL,IT,IS1,IORB,1,
     &                         ISPN),IS1=1,2)
               WRITE (6,*) '--- LSX ---'
               WRITE (6,99002) (IL,0D0,0D0,D2X(IL,IT,IS1,IORB,1,ISPN),
     &                         IS1,D3X(IL,IT,IS1,IORB,1,ISPN),
     &                         D4X(IL,IT,IS1,IORB,1,ISPN),IS1=1,2)
               WRITE (6,*) '--- HFSZ ---'
               WRITE (6,99001) IL,D1Z(IL,IT,1,IHFI,1,ISPN),
     &                         D2Z(IL,IT,1,IHFI,1,ISPN),
     &                         D3Z(IL,IT,1,IHFI,1,ISPN),
     &                         D4Z(IL,IT,1,IHFI,1,ISPN)
               WRITE (6,*) '--- DHFSZL ---'
               WRITE (6,*) 'IL=',IL,' DHFSZL=',DZL(IL,IT,1,IHFI,1,ISPN)
               WRITE (6,*) '--- HFSX ---'
               WRITE (6,99001) IL,0D0,0D0,D2X(IL,IT,1,IHFI,1,ISPN),
     &                         D3X(IL,IT,1,IHFI,1,ISPN),
     &                         D4X(IL,IT,1,IHFI,1,ISPN)
               WRITE (6,*) '--- DHFSXL ---'
               WRITE (6,*) 'IL=',IL,' DHFSXL=',
     &                     DIJXL(IL,IT,1,1,IHFI,1,ISPN)
               WRITE (6,*) '--- HFLZ ---'
               WRITE (6,99001) IL,
     &                         CSUMSSLT(D1Z,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D2Z,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D3Z,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D4Z,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP)
               WRITE (6,*) '--- DHFLZL ---'
               WRITE (6,*) 'IL=',IL,' DHFLZL=',DZL(IL,IT,1,IHFI,1,IORB)
     &                     + DZL(IL,IT,1,IHFI,2,IORB)
               WRITE (6,*) '--- HFLZ2 ---'
               WRITE (6,99001) IL,0D0,0D0,
     &                         CSUMSSLT(D2X,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D3X,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D4X,IL,IT,IHFI,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP)
               WRITE (6,*) '--- DHFLXL ---'
               WRITE (6,*) 'IL=',IL,' DHFLXL=',
     &                     DIJXL(IL,IT,1,1,IHFI,1,IORB)
     &                     + DIJXL(IL,IT,1,1,IHFI,2,IORB)
               WRITE (6,*) '--- L3L ---'
               WRITE (6,99001) IL,D1Z(IL,IT,1,IHVV,1,IORB)
     &                         + D1Z(IL,IT,1,IHVV,2,IORB),
     &                         CSUMSSLT(D2Z,IL,IT,IHVV,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D3Z,IL,IT,IHVV,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP),
     &                         CSUMSSLT(D4Z,IL,IT,IHVV,IORB,NLMAX,NTMAX,
     &                         1,NOBS,2,NOP)
               WRITE (6,*) '--- DL3LZL ---'
               WRITE (6,*) 'IL=',IL,' DL3LZL=',DZL(IL,IT,1,IHVV,1,IORB)
     &                     + DZL(IL,IT,1,IHVV,2,IORB)
               WRITE (6,*) '--- DKDIAL ---'
               WRITE (6,*) 'IL=',IL,' DKDIAL=',DOBS(IL,IT,1,IKDIA)*PI
               WRITE (6,*) '--- LLZ ---'
               WRITE (6,99002) (IL,D1Z(IL,IT,IS1,IORB,IS1,IORB),D2Z(IL,
     &                         IT,IS1,IORB,IS1,IORB),IS1,
     &                         D3Z(IL,IT,IS1,IORB,IS1,IORB),
     &                         D4Z(IL,IT,IS1,IORB,IS1,IORB),IS1=1,2)
               WRITE (6,*) '--- DLLZSL ---'
               WRITE (6,99003) (IL,IS1,'DLLZSL',DZL(IL,IT,IS1,IORB,IS1,
     &                         IORB),IS1=1,2)
               WRITE (6,*) '--- D1LLZSL ---'
               WRITE (6,99003) (IL,IS1,'D1LLZSL',D1Z(IL,IT,IS1,IORB,IS1,
     &                         IORB),IS1=1,2)
               WRITE (6,*) '--- D0LLZSL ---'
               WRITE (6,99003) (IL,IS1,'D0LLZSL',D0ZL(IL,IT,IS1,IORB,IS1
     &                         ,IORB),IS1=1,2)
C
               DO IS2 = 1,2
                  WRITE (6,'(A,I2,A)') '--- LLX (IS2=',IS2,') ---'
                  WRITE (6,99002) (IL,0D0,0D0,D2X(IL,IT,IS1,IORB,IS2,
     &                            IORB),IS1,D3X(IL,IT,IS1,IORB,IS2,IORB)
     &                            ,D4X(IL,IT,IS1,IORB,IS2,IORB),IS1=1,2)
                  WRITE (6,'(A,I2,A)') '--- D0LLXSL (IS2=',IS2,') ---'
                  WRITE (6,99003) (IL,IS1,'D0LLXSL',D0XL(IL,IT,IS1,IORB,
     &                            IS2,IORB),IS1=1,2)
                  DO JT = 1,NT
                     WRITE (6,'(A,I2,A,I2,A)') '--- DIJLLXSL (IS2=',IS2,
     &                      ' JT=',JT,') ---'
                     WRITE (6,99003) (IL,IS1,'DIJLLXSL',DIJXL(IL,IT,JT,
     &                               IS1,IORB,IS2,IORB),IS1=1,2)
                  END DO
               END DO
            END IF
C ======================================================================
         END DO
      END DO
99001 FORMAT ('il=',I2,' t1=',2E16.8,'  t2=',2E16.8,/,5X,' t3=',2E16.8,
     &        '  t4=',2E16.8)
99002 FORMAT (' il=',I2,' t1=',2E16.8,'  t2=',2E16.8,/,'is1=',I2,' t3=',
     &        2E16.8,'  t4=',2E16.8)
99003 FORMAT ('il=',I2,' is1=',I2,' ',A,'=',2E16.8)
C
      END
C*==chiterm1.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHITERM1(HAZR,HBZR,HCZR,HDZR,HAZ,HBZ,HCZ,HDZ,HAX,HBX,
     &                    HCX,HDX,ZG,JG,WGTG,WGWG,AMEG,ZF,JF,WGTF,WFWF,
     &                    AMEF,CR1,CR2DRDI,KXR2DRDI,NSOL,IRTOP,NSPINOBS,
     &                    NOBS,NOP,NRMAX,NCHI)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHITERM11164
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,NOBS,NOP,NRMAX,NSOL
      REAL*8 AMEF(2,2,2,NOBS),AMEG(2,2,2,NOBS),WGTF(NOBS),WGTG(NOBS)
      COMPLEX*16 CR1(NRMAX),CR2DRDI(NRMAX),HAX(2,2,NRMAX,2,NOBS,2),
     &           HAZ(2,2,NRMAX,2,NOBS),HAZR(2,2,NRMAX,2,NOBS),
     &           HBX(2,2,NRMAX,2,NOBS,2),HBZ(2,2,NRMAX,2,NOBS),
     &           HBZR(2,2,NRMAX,2,NOBS),HCX(2,2,NRMAX,2,NOBS,2),
     &           HCZ(2,2,NRMAX,2,NOBS),HCZR(2,2,NRMAX,2,NOBS),
     &           HDX(2,2,NRMAX,2,NOBS,2),HDZ(2,2,NRMAX,2,NOBS),
     &           HDZR(2,2,NRMAX,2,NOBS),JF(NRMAX,2,2),JG(NRMAX,2,2),
     &           KXR2DRDI(NRMAX,2,NOP,2),WFWF(2,2,NRMAX),WGWG(2,2,NRMAX)
     &           ,ZF(NRMAX,2,2),ZG(NRMAX,2,2)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS)
C
C Local variables
C
      INTEGER ICHI,IOBS,IOP,ISPIN
C
C*** End of declarations rewritten by SPAG
C
      CALL CINIT(2*2*NRMAX*2*NOBS,HAZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HBZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HCZ)
      CALL CINIT(2*2*NRMAX*2*NOBS,HDZ)
C
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HAX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HBX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HCX)
      CALL CINIT(2*2*NRMAX*2*NOBS*2,HDX)
C
      DO IOBS = 1,NOP
         DO ISPIN = 1,NSPINOBS(IOBS)
C
C------------------------------------------------------- H_op Zeeman (R)
C
            CALL CHIINTGABR3(HAZR(1,1,1,ISPIN,IOBS),ZG,ZG,WGTG(IOBS),
     &                       WGWG,AMEG(1,1,ISPIN,IOBS),ZF,ZF,WGTF(IOBS),
     &                       WFWF,AMEF(1,1,ISPIN,IOBS),CR1,NSOL,IRTOP,
     &                       NRMAX)
C
            CALL CHIINTGABR3(HBZR(1,1,1,ISPIN,IOBS),JG,ZG,WGTG(IOBS),
     &                       WGWG,AMEG(1,1,ISPIN,IOBS),JF,ZF,WGTF(IOBS),
     &                       WFWF,AMEF(1,1,ISPIN,IOBS),CR1,NSOL,IRTOP,
     &                       NRMAX)
C
            CALL CHIINTGABR3(HCZR(1,1,1,ISPIN,IOBS),ZG,JG,WGTG(IOBS),
     &                       WGWG,AMEG(1,1,ISPIN,IOBS),ZF,JF,WGTF(IOBS),
     &                       WFWF,AMEF(1,1,ISPIN,IOBS),CR1,NSOL,IRTOP,
     &                       NRMAX)
C
            CALL CHIINTGABR3(HDZR(1,1,1,ISPIN,IOBS),JG,JG,WGTG(IOBS),
     &                       WGWG,AMEG(1,1,ISPIN,IOBS),JF,JF,WGTF(IOBS),
     &                       WFWF,AMEF(1,1,ISPIN,IOBS),CR1,NSOL,IRTOP,
     &                       NRMAX)
C
C----------------------------------------------------------- H_op Zeeman
C
            CALL CHIINTABR3(HAZ(1,1,1,ISPIN,IOBS),ZG,ZG,WGTG(IOBS),WGWG,
     &                      AMEG(1,1,ISPIN,IOBS),ZF,ZF,WGTF(IOBS),WFWF,
     &                      AMEF(1,1,ISPIN,IOBS),CR2DRDI,NSOL,IRTOP,
     &                      NRMAX)
C
            CALL CHIINTABR3(HBZ(1,1,1,ISPIN,IOBS),JG,ZG,WGTG(IOBS),WGWG,
     &                      AMEG(1,1,ISPIN,IOBS),JF,ZF,WGTF(IOBS),WFWF,
     &                      AMEF(1,1,ISPIN,IOBS),CR2DRDI,NSOL,IRTOP,
     &                      NRMAX)
C
            CALL CHIINTABR3(HCZ(1,1,1,ISPIN,IOBS),ZG,JG,WGTG(IOBS),WGWG,
     &                      AMEG(1,1,ISPIN,IOBS),ZF,JF,WGTF(IOBS),WFWF,
     &                      AMEF(1,1,ISPIN,IOBS),CR2DRDI,NSOL,IRTOP,
     &                      NRMAX)
C
            CALL CHIINTABR3IN(HDZ(1,1,1,ISPIN,IOBS),JG,JG,WGTG(IOBS),
     &                        WGWG,AMEG(1,1,ISPIN,IOBS),JF,JF,WGTF(IOBS)
     &                        ,WFWF,AMEF(1,1,ISPIN,IOBS),CR2DRDI,NSOL,
     &                        IRTOP,NRMAX)
C
C--------------------------------------------------------- H_op eXchange
C
C     ICHI = IDOS  related to the term              1*(ENN*chi_n + ENM*chi_s)
C     ICHI = ISPN  related to the term   beta*sigma_z*(EMN*chi_n + EMM*chi_s)
C
            IOP = IOBS
            IF ( IOP.LT.3 ) THEN
               DO ICHI = 1,NCHI(IOP)
C
                  CALL CHIINTABR3(HAX(1,1,1,ISPIN,IOP,ICHI),ZG,ZG,
     &                            WGTG(ICHI),WGWG,AMEG(1,1,ISPIN,ICHI),
     &                            ZF,ZF,WGTF(ICHI),WFWF,
     &                            AMEF(1,1,ISPIN,ICHI),
     &                            KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                            NRMAX)
C
                  CALL CHIINTABR3(HBX(1,1,1,ISPIN,IOP,ICHI),JG,ZG,
     &                            WGTG(ICHI),WGWG,AMEG(1,1,ISPIN,ICHI),
     &                            JF,ZF,WGTF(ICHI),WFWF,
     &                            AMEF(1,1,ISPIN,ICHI),
     &                            KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                            NRMAX)
C
                  CALL CHIINTABR3(HCX(1,1,1,ISPIN,IOP,ICHI),ZG,JG,
     &                            WGTG(ICHI),WGWG,AMEG(1,1,ISPIN,ICHI),
     &                            ZF,JF,WGTF(ICHI),WFWF,
     &                            AMEF(1,1,ISPIN,ICHI),
     &                            KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                            NRMAX)
C
                  CALL CHIINTABR3IN(HDX(1,1,1,ISPIN,IOP,ICHI),JG,JG,
     &                              WGTG(ICHI),WGWG,AMEG(1,1,ISPIN,ICHI)
     &                              ,JF,JF,WGTF(ICHI),WFWF,
     &                              AMEF(1,1,ISPIN,ICHI),
     &                              KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,
     &                              IRTOP,NRMAX)
               END DO
            ELSE
               ICHI = 1
               CALL CHIINTABR3(HAX(1,1,1,ISPIN,IOP,ICHI),ZG,ZG,WGTG(IOP)
     &                         ,WGWG,AMEG(1,1,ISPIN,IOP),ZF,ZF,WGTF(IOP)
     &                         ,WFWF,AMEF(1,1,ISPIN,IOP),
     &                         KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                         NRMAX)
C
               CALL CHIINTABR3(HBX(1,1,1,ISPIN,IOP,ICHI),JG,ZG,WGTG(IOP)
     &                         ,WGWG,AMEG(1,1,ISPIN,IOP),JF,ZF,WGTF(IOP)
     &                         ,WFWF,AMEF(1,1,ISPIN,IOP),
     &                         KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                         NRMAX)
C
               CALL CHIINTABR3(HCX(1,1,1,ISPIN,IOP,ICHI),ZG,JG,WGTG(IOP)
     &                         ,WGWG,AMEG(1,1,ISPIN,IOP),ZF,JF,WGTF(IOP)
     &                         ,WFWF,AMEF(1,1,ISPIN,IOP),
     &                         KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                         NRMAX)
C
               CALL CHIINTABR3IN(HDX(1,1,1,ISPIN,IOP,ICHI),JG,JG,
     &                           WGTG(IOP),WGWG,AMEG(1,1,ISPIN,IOP),JF,
     &                           JF,WGTF(IOP),WFWF,AMEF(1,1,ISPIN,IOP),
     &                           KXR2DRDI(1,ISPIN,IOP,ICHI),NSOL,IRTOP,
     &                           NRMAX)
            END IF
         END DO
      END DO
C
      END
C*==chiint.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
C
C
      SUBROUTINE CHIINT(DSUMZ,DSUMX,AG,BG,WGTG,WGWG,AMEG,AF,BF,WGTF,
     &                  WFWF,AMEF,WGTROOZ,WGTROOX,NSOL,IRTOP,NSPINOBS,
     &                  NOBS,NSPINOP,NOP,NRMAX,NCHI)
C   ********************************************************************
C   *                                                                  *
C   *      INTEGRAL FROM 0 TO R_IRTOP OF PRODUCT OF                     *
C   *        REGULAR/IRREGULAR SOLUTIONS OF DIRAC EQUATION             *
C   *        (E.G. Z(A)Z(B), Z(A)J(B), J(A)Z(B) UND J(A)J(B))          *
C   *        BY CALCULATING THE INTEGRALS G*G AND F*F                  *
C   *        AND PROPERLY SUMMING OVER INTERNAL COUPLINGS.             *
C   *      USED FOR CALCULATING DOUBLE INTEGRAL IN CHIRADINT WITH      *
C   *      INNER R-DEPENDENT INTEGRAL AS WEIGHT FUNCTION.              *
C   *                                                                  *
C   *      HF/MD: 12/96                                                *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHIINT1344
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,NOBS,NOP,NRMAX,NSOL
      COMPLEX*16 AF(NRMAX,2),AG(NRMAX,2),BF(NRMAX,2),BG(NRMAX,2),
     &           DSUMX(2,NOBS,2,NOP,2),DSUMZ(2,NOBS,2,NOP),
     &           WFWF(2,2,NRMAX),WGTROOX(NRMAX,2,NOBS,2,NOP,2),
     &           WGTROOZ(NRMAX,2,NOBS,2,NOP),WGWG(2,2,NRMAX)
      REAL*8 AMEF(2,2,2,NOBS),AMEG(2,2,2,NOBS),WGTF(NOBS),WGTG(NOBS)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      INTEGER ICHI,IOBS,IOP,IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
      DO IOBS = 1,NOP
         DO IS1 = 1,NSPINOBS(IOBS)
            DO IOP = 1,NOP
               DO IS2 = 1,NSPINOP(IOP)
C
                  CALL CHIINTABR4(DSUMZ(IS1,IOBS,IS2,IOP),AG,BG,
     &                            WGTG(IOBS),WGWG,AMEG(1,1,IS1,IOBS),AF,
     &                            BF,WGTF(IOBS),WFWF,AMEF(1,1,IS1,IOBS),
     &                            WGTROOZ(1,IS1,IOBS,IS2,IOP),NSOL,
     &                            IRTOP,NRMAX)
C
C
C     ICHI = IDOS  related to the term              1*(ENN*chi_n + ENM*chi_s)
C     ICHI = ISPN  related to the term   beta*sigma_z*(EMN*chi_n + EMM*chi_s)
C
                  DO ICHI = 1,NCHI(IOP)
                     CALL CHIINTABR4(DSUMX(IS1,IOBS,IS2,IOP,ICHI),AG,BG,
     &                               WGTG(IOBS),WGWG,AMEG(1,1,IS1,IOBS),
     &                               AF,BF,WGTF(IOBS),WFWF,
     &                               AMEF(1,1,IS1,IOBS),
     &                               WGTROOX(1,IS1,IOBS,IS2,IOP,ICHI),
     &                               NSOL,IRTOP,NRMAX)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      END
C*==chiinthf.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTHF(IM,DSUMZ,DSUMX,AG,BG,AF,BF,WGTHF,RMEHF,AMEHF,
     &                    WGTROOZ,WGTROOX,NSOL,IRTOP,NSPINOBS,NOBS,
     &                    NSPINOP,NOP,IHFI,NCHI)
C   ********************************************************************
C   *                                                                  *
C   *      INTEGRAL FROM 0 TO R_IRTOP OF PRODUCT OF                     *
C   *        REGULAR/IRREGULAR SOLUTIONS OF DIRAC EQUATION             *
C   *        (E.G. Z(A)Z(B), Z(A)J(B), J(A)Z(B) UND J(A)J(B))          *
C   *        BY CALCULATING THE INTEGRALS G*G AND F*F                  *
C   *        AND PROPERLY SUMMING OVER INTERNAL COUPLINGS.             *
C   *      USED FOR CALCULATING DOUBLE INTEGRAL IN CHIRADINT WITH      *
C   *      INNER R-DEPENDENT INTEGRAL AS WEIGHT FUNCTION.              *
C   *                                                                  *
C   *      HF/MD: 12/96                                                *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--CHIINTHF1424
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IHFI,IM,IRTOP,NOBS,NOP,NSOL
      REAL*8 WGTHF
      COMPLEX*16 AF(NRMAX,2),AG(NRMAX,2),BF(NRMAX,2),BG(NRMAX,2),
     &           DSUMX(2,NOBS,2,NOP,2),DSUMZ(2,NOBS,2,NOP),
     &           RMEHF(NCPLWFMAX,NCPLWFMAX),
     &           WGTROOX(NRMAX,2,NOBS,2,NOP,2),
     &           WGTROOZ(NRMAX,2,NOBS,2,NOP)
      REAL*8 AMEHF(2,2)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      INTEGER ICHI,IOP,IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
      DO IS1 = 1,NSPINOBS(IHFI)
         DO IOP = 1,NOP
            DO IS2 = 1,NSPINOP(IOP)
C
               CALL CHIINTHFF4(IM,DSUMZ(IS1,IHFI,IS2,IOP),AG,BG,AF,BF,
     &                         WGTHF,RMEHF,AMEHF,
     &                         WGTROOZ(1,IS1,IHFI,IS2,IOP),NSOL,IRTOP)
C
C ICHI = IDOS  related to the term             1*(ENN*chi_n + ENM*chi_s)
C ICHI = ISPN  related to the term  beta*sigma_z*(EMN*chi_n + EMM*chi_s)
C
               DO ICHI = 1,NCHI(IOP)
                  CALL CHIINTHFF4(IM,DSUMX(IS1,IHFI,IS2,IOP,ICHI),AG,BG,
     &                            AF,BF,WGTHF,RMEHF,AMEHF,
     &                            WGTROOX(1,IS1,IHFI,IS2,IOP,ICHI),NSOL,
     &                            IRTOP)
               END DO
C
            END DO
         END DO
      END DO
      END
C*==chiwgtr.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIWGTR(WGTZ,WGTX,WROBS,HZ,HX,KEY,K1,K2,IRTOP,NSPINOBS,
     &                   NOBS,NSPINOP,NOP,NRMAX,NCHI)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHIWGTR1487
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,K1,K2,NOBS,NOP,NRMAX
      CHARACTER*3 KEY
      COMPLEX*16 HX(2,2,NRMAX,2,NOBS,2),HZ(2,2,NRMAX,2,NOBS),
     &           WGTX(NRMAX,2,NOBS,2,NOP,2),WGTZ(NRMAX,2,NOBS,2,NOP)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS),NSPINOP(NOP)
      REAL*8 WROBS(NRMAX,NOBS)
C
C Local variables
C
      INTEGER I,ICHI,IOBS,IOP,IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
      IF ( KEY.EQ.'out' ) THEN
         DO IOBS = 1,NOBS
            DO IS1 = 1,NSPINOBS(IOBS)
               DO IOP = 1,NOP
                  DO IS2 = 1,NSPINOP(IOP)
                     DO I = 1,IRTOP
                        WGTZ(I,IS1,IOBS,IS2,IOP) = WROBS(I,IOBS)
     &                     *HZ(K1,K2,I,IS2,IOP)
C
                        DO ICHI = 1,NCHI(IOP)
                           WGTX(I,IS1,IOBS,IS2,IOP,ICHI) = WROBS(I,IOBS)
     &                        *HX(K1,K2,I,IS2,IOP,ICHI)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE IF ( KEY.EQ.'in ' ) THEN
         DO IOBS = 1,NOBS
            DO IS1 = 1,NSPINOBS(IOBS)
               DO IOP = 1,NOP
                  DO IS2 = 1,NSPINOP(IOP)
                     DO I = 1,IRTOP
                        WGTZ(I,IS1,IOBS,IS2,IOP) = WROBS(I,IOBS)
     &                     *(HZ(K1,K2,IRTOP,IS2,IOP)-HZ(K1,K2,I,IS2,IOP)
     &                     )
                        DO ICHI = 1,NCHI(IOP)
                           WGTX(I,IS1,IOBS,IS2,IOP,ICHI) = WROBS(I,IOBS)
     &                        *(HX(K1,K2,IRTOP,IS2,IOP,ICHI)
     &                        -HX(K1,K2,I,IS2,IOP,ICHI))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      END
C*==chiintprep.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIINTPREP(DZR,DXR,HZR,HZ,HX,KEY,K1,K2,K3,K4,IRTOP,
     &                      NSPINOBS,NOBS,NSPINOP,NOP,NRMAX,NCHI)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHIINTPREP1565
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,K1,K2,K3,K4,NOBS,NOP,NRMAX
      CHARACTER*3 KEY
      COMPLEX*16 DXR(NRMAX,2,NOBS,2,NOP,2),DZR(NRMAX,2,NOBS,2,NOP),
     &           HX(2,2,NRMAX,2,NOBS,2),HZ(2,2,NRMAX,2,NOBS),
     &           HZR(2,2,NRMAX,2,NOBS)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      INTEGER I,ICHI,IOBS,IOP,IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
      IF ( KEY.EQ.'out' ) THEN
         DO IOBS = 1,NOBS
            DO IS1 = 1,NSPINOBS(IOBS)
               DO IOP = 1,NOP
                  DO IS2 = 1,NSPINOP(IOP)
                     DO I = 1,IRTOP
                        DZR(I,IS1,IOBS,IS2,IOP) = HZR(K1,K2,I,IS1,IOBS)
     &                     *HZ(K3,K4,I,IS2,IOP)
                        DO ICHI = 1,NCHI(IOP)
                           DXR(I,IS1,IOBS,IS2,IOP,ICHI)
     &                        = HZR(K1,K2,I,IS1,IOBS)
     &                        *HX(K3,K4,I,IS2,IOP,ICHI)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE IF ( KEY.EQ.'in ' ) THEN
         DO IOBS = 1,NOBS
            DO IS1 = 1,NSPINOBS(IOBS)
               DO IOP = 1,NOP
                  DO IS2 = 1,NSPINOP(IOP)
                     DO I = 1,IRTOP
                        DZR(I,IS1,IOBS,IS2,IOP) = HZR(K1,K2,I,IS1,IOBS)
     &                     *(HZ(K3,K4,IRTOP,IS2,IOP)-HZ(K3,K4,I,IS2,IOP)
     &                     )
                        DO ICHI = 1,NCHI(IOP)
                           DXR(I,IS1,IOBS,IS2,IOP,ICHI)
     &                        = HZR(K1,K2,I,IS1,IOBS)
     &                        *(HX(K3,K4,IRTOP,IS2,IOP,ICHI)
     &                        -HX(K3,K4,I,IS2,IOP,ICHI))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      END
C*==chidsum.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIDSUM(DZ,DX,DRZ,DRX,DAZ,DBZ,DAX,DBX,DARZ,DBRZ,DARX,
     &                   DBRX,TAU,IL,IT,IRTOP,NLMAX,NTMAX,NSPINOBS,NOBS,
     &                   NSPINOP,NOP,NRMAX,NCHI)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHIDSUM1645
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IL,IRTOP,IT,NLMAX,NOBS,NOP,NRMAX,NTMAX
      COMPLEX*16 TAU
      COMPLEX*16 DARX(NRMAX,2,NOBS,2,NOP,2),DARZ(NRMAX,2,NOBS,2,NOP),
     &           DAX(2,NOBS,2,NOP,2),DAZ(2,NOBS,2,NOP),
     &           DBRX(NRMAX,2,NOBS,2,NOP,2),DBRZ(NRMAX,2,NOBS,2,NOP),
     &           DBX(2,NOBS,2,NOP,2),DBZ(2,NOBS,2,NOP),
     &           DRX(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP),
     &           DRZ(NRMAX,NLMAX,NTMAX,2,NOBS,2,NOP),
     &           DX(NLMAX,NTMAX,2,NOBS,2,NOP),
     &           DZ(NLMAX,NTMAX,2,NOBS,2,NOP)
      INTEGER NCHI(NOBS),NSPINOBS(NOBS),NSPINOP(NOP)
C
C Local variables
C
      INTEGER I,ICHI,IOBS,IOP,IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
      DO IOBS = 1,NOBS
         DO IS1 = 1,NSPINOBS(IOBS)
            DO IOP = 1,NOP
               DO IS2 = 1,NSPINOP(IOP)
C
                  DZ(IL,IT,IS1,IOBS,IS2,IOP)
     &               = DZ(IL,IT,IS1,IOBS,IS2,IOP)
     &               + TAU*(DAZ(IS1,IOBS,IS2,IOP)+DBZ(IS1,IOBS,IS2,IOP))
C
                  DO I = 1,IRTOP
                     DRZ(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                  = DRZ(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                  + TAU*(DARZ(I,IS1,IOBS,IS2,IOP)
     &                  +DBRZ(I,IS1,IOBS,IS2,IOP))
                  END DO
C
                  DO ICHI = 1,NCHI(IOP)
                     DX(IL,IT,IS1,IOBS,IS2,IOP)
     &                  = DX(IL,IT,IS1,IOBS,IS2,IOP)
     &                  + TAU*(DAX(IS1,IOBS,IS2,IOP,ICHI)
     &                  +DBX(IS1,IOBS,IS2,IOP,ICHI))
C
                     DO I = 1,IRTOP
                        DRX(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                     = DRX(I,IL,IT,IS1,IOBS,IS2,IOP)
     &                     + TAU*(DARX(I,IS1,IOBS,IS2,IOP,ICHI)
     &                     +DBRX(I,IS1,IOBS,IS2,IOP,ICHI))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      END
