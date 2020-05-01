C*==comptonme1.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COMPTONME1(NCPLWF_T,IKMCPLWF_T,NP,PFE,PFEMAX,ME_PWZR,
     &                      ME_ZLPW,ME_IRR,RME_PW_ZR_P,RME_ZL_PW_P,
     &                      RME_IRR_P,PFEDIR,EPSP,NPMAX)
C     ******************************************************************
C     *                                                                *
C     *                         PART 1                                 *
C     *                                                                *
C     *   speed up the calculation of the matrix elements by           *
C     *   interpolating the p-dependent part                           *
C     *                                                                *
C     *----------------------------------------------------------------*
C     *                                                                *
C     *   calculate the FULL matrix elements                           *
C     *   for the momentum representation                              *
C     *                                                                *
C     *            ME_PWZR = < p m_s | Z_Lam^R >                       *
C     *                                                                *
C     *            ME_ZLPW = < Z_Lam^L | p m_s >                       *
C     *                                                                *
C     ******************************************************************
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_ANGMOM,ONLY:NL,NKMMAX,CGC,L_IKM,NKM,MUEM05_IKM
      USE MOD_TYPES,ONLY:NT,NTMAX,NCPLWFMAX
      USE MOD_CONSTANTS,ONLY:CI,C0
      IMPLICIT NONE
C*--COMPTONME127
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTONME1')
C
C Dummy arguments
C
      REAL*8 EPSP,PFE,PFEMAX
      INTEGER NP,NPMAX
      INTEGER IKMCPLWF_T(NCPLWFMAX,NKMMAX,NTMAX),NCPLWF_T(NKMMAX,NTMAX)
      COMPLEX*16 ME_IRR(2,NTMAX),ME_PWZR(NKMMAX,2,NTMAX),
     &           ME_ZLPW(NKMMAX,2,NTMAX),
     &           RME_IRR_P(NPMAX,NKMMAX,NKMMAX,NTMAX),
     &           RME_PW_ZR_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX),
     &           RME_ZL_PW_P(NPMAX,NCPLWFMAX,NKMMAX,NTMAX)
      REAL*8 PFEDIR(3)
C
C Local variables
C
      REAL*8 DIR(3),DPFE,DX,DXSQ,MJ,MJ1,MJ2,WLAG(NPMAX),XIP
      COMPLEX*16 FPW1,FPW2,IPWML,IPWML1,IPWML2,PWZG,RME_IRR,YLMP(:),ZLPW
      INTEGER IA_ERR,IKM,IKM1,IKM2,ILM,ILM1,ILM2,IP,IPBOT,IPTOP,IS,IT,
     &        KK,L,L1,L2,LAM,LCHPMAX,ML,ML1,ML2,NLMCHPMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE YLMP
C
      CALL TRACK_INFO(ROUTINE)
C
C-----------------------------------------------------------------------
C               initialize complex spherical harmonics
C-----------------------------------------------------------------------
C
      LCHPMAX = NL - 1
      NLMCHPMAX = (LCHPMAX+1)**2
C
      ALLOCATE (YLMP(NLMCHPMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CHP')
C
      DIR(1:3) = PFEDIR(1:3)
      CALL RVECNORM(3,DIR)
C
      CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),YLMP,LCHPMAX,NLMCHPMAX)
C
C-----------------------------------------------------------------------
C
      ME_PWZR(:,:,:) = C0
      ME_ZLPW(:,:,:) = C0
      ME_IRR(:,:) = C0
C
      DPFE = PFEMAX/DBLE(NP-1)
      XIP = 1D0 + PFE/DPFE
      IP = MIN(INT(XIP),NP-2)
      IPBOT = IP
      IPTOP = IP + 2
C
      DX = XIP - IPBOT
      DXSQ = DX*DX
      WLAG(IPBOT+0) = 0.5D0*(2-3*DX+DXSQ)
      WLAG(IPBOT+1) = 2*DX - DXSQ
      WLAG(IPBOT+2) = 0.5D0*(-DX+DXSQ)
C
C  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
C***********************************************************************
C           REGULAR PART   < PW | Z > and  < Z | PW >
C***********************************************************************
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO LAM = 1,NKM
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
            DO KK = 1,NCPLWF_T(LAM,IT)
C
               PWZG = 0D0
               ZLPW = 0D0
               DO IP = IPBOT,IPTOP
                  PWZG = PWZG + WLAG(IP)*RME_PW_ZR_P(IP,KK,LAM,IT)
                  ZLPW = ZLPW + WLAG(IP)*RME_ZL_PW_P(IP,KK,LAM,IT)
               END DO
C
               IKM = IKMCPLWF_T(KK,LAM,IT)
C
               L = L_IKM(IKM)
               IPWML = 1.0D0/CI**L
               MJ = DBLE(MUEM05_IKM(IKM)) + 0.5D0
C
               DO IS = 1,2
C
                  ML = NINT(MJ-(IS-1.5D0))
                  ILM = L**2 + L + ML + 1
                  IF ( ABS(ML).LE.L ) THEN
C
                     ME_PWZR(LAM,IS,IT) = ME_PWZR(LAM,IS,IT)
     &                  + EPSP*IPWML*CGC(IKM,IS)*YLMP(ILM)*PWZG
C
                     ME_ZLPW(LAM,IS,IT) = ME_ZLPW(LAM,IS,IT)
     &                  + EPSP*CGC(IKM,IS)*DCONJG(IPWML*YLMP(ILM))*ZLPW
C
                  END IF
C
               END DO
            END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C
C***********************************************************************
C           IRRREGULAR PART   < PW | ZJ + JZ | PW >
C***********************************************************************
C           NOTE: a summation is done w.r.t to the index LAM
C
         IF ( .NOT.SPLITSS ) THEN
            DO IKM1 = 1,NKM
C
               L1 = L_IKM(IKM1)
               IPWML1 = 1.0D0/CI**L1
               MJ1 = DBLE(MUEM05_IKM(IKM1)) + 0.5D0
C
               DO IKM2 = 1,NKM
C
                  L2 = L_IKM(IKM2)
                  IPWML2 = 1.0D0/CI**L2
                  MJ2 = DBLE(MUEM05_IKM(IKM2)) + 0.5D0
C
                  RME_IRR = 0D0
                  DO IP = IPBOT,IPTOP
                     RME_IRR = RME_IRR + WLAG(IP)
     &                         *RME_IRR_P(IP,IKM1,IKM2,IT)
                  END DO
C
                  DO IS = 1,2
C
                     ML1 = NINT(MJ1-(IS-1.5D0))
                     ILM1 = L1**2 + L1 + ML1 + 1
C
                     ML2 = NINT(MJ2-(IS-1.5D0))
                     ILM2 = L2**2 + L2 + ML2 + 1
C
                     IF ( ABS(ML1).LE.L1 .AND. ABS(ML2).LE.L2 ) THEN
C
                        FPW1 = EPSP*IPWML1*CGC(IKM1,IS)*YLMP(ILM1)
                        FPW2 = EPSP*CGC(IKM2,IS)
     &                         *DCONJG(IPWML2*YLMP(ILM2))
C
                        ME_IRR(IS,IT) = ME_IRR(IS,IT)
     &                                  + FPW1*FPW2*RME_IRR
C
                     END IF
C
                  END DO
               END DO
            END DO
C
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
