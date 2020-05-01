C*==compton_core1.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COMPTON_CORE1(NCPLWF_CORT,IKMCPLWF_CORT,NP,PFE,PFEMAX,
     &                         RME_COR_P,MDCOR_ST,PFEDIR,EPSP,NPMAX)
C     ******************************************************************
C     *                                                                *
C     *                    CORE CONTRIBUTION                           *
C     *                                                                *
C     *     to Compton or positron or positron annihilation spectrum   *
C     *                                                                *
C     *----------------------------------------------------------------*
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
C     *            ME_COR = < p m_s | Z_Lam >                          *
C     *                                                                *
C     ******************************************************************
      USE MOD_ANGMOM,ONLY:NL,CGC,L_IKM,MUEM05_IKM
      USE MOD_TYPES,ONLY:NT,NTMAX,NCPLWFMAX,NCORMAX,NCORT
      USE MOD_CONSTANTS,ONLY:CI
      IMPLICIT NONE
C*--COMPTON_CORE129
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
      INTEGER IKMCPLWF_CORT(NCPLWFMAX,NCORMAX,NTMAX),
     &        NCPLWF_CORT(NCORMAX,NTMAX)
      REAL*8 MDCOR_ST(2,NTMAX),PFEDIR(3)
      COMPLEX*16 RME_COR_P(NPMAX,NCPLWFMAX,NCORMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CADD,IPWML,PWZG,YLMP(:)
      REAL*8 DIR(3),DPFE,DX,DXSQ,MJ,WLAG(NPMAX),XIP
      INTEGER IA_ERR,ICOR,IKM,ILM,IP,IPBOT,IPTOP,IS,IT,KK,L,LCHPMAX,ML,
     &        NLMCHPMAX
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
      MDCOR_ST(:,:) = 0D0
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
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO ICOR = 1,NCORT(IT)
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
            DO KK = 1,NCPLWF_CORT(ICOR,NTMAX)
C
               PWZG = 0D0
               DO IP = IPBOT,IPTOP
                  PWZG = PWZG + WLAG(IP)*RME_COR_P(IP,KK,ICOR,IT)
               END DO
C
               IKM = IKMCPLWF_CORT(KK,ICOR,IT)
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
                     CADD = EPSP*IPWML*CGC(IKM,IS)*YLMP(ILM)*PWZG
C
                     MDCOR_ST(IS,IT) = MDCOR_ST(IS,IT)
     &                                 + DREAL(CADD*DCONJG(CADD))
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
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
