C*==me_spin_nab_alf.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SPIN_NAB_ALF(IFILA,ERYDA,IFILB,ERYDB,IT,MZAZB,MZBZA,
     &                           MIRR_2,MIRR_3,MIRR_4,C)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *    MZAZB  =   < Z^+(E_a) | H_{lam,lam_spin} | Z(E_b) >           *
C   *                                                                  *
C   *    MZBZA  =   < Z^+(E_b) | H_{lam,lam_spin} | Z(E_a) >           *
C   *                                                                  *
C   *    MIRR_2 =   < Z^+(E_b) | H_{lam,lam_spin} * IZAZB | J(E_a) >   *
C   *             + < Z^+(E_b) | H_{lam,lam_spin} * IJAZB | Z(E_a) >   *
C   *                                                                  *
C   *    MIRR_3 =   < J^+(E_b) | H_{lam,lam_spin} * IZAZB | J(E_a) >   *
C   *             + < Z^+(E_b) | H_{lam,lam_spin} * IZAJB | Z(E_a) >   *
C   *                                                                  *
C   *    MIRR_4 =   < J^+(E_b) | H_{lam,lam_spin} * IZAZB | J(E_a) >   *
C   *             + < Z^+(E_b) | H_{lam,lam_spin} * IZAJB | Z(E_a) >   *
C   *                                                                  *
C   *     H_{lam,lam_spin}   spin current operator in the              *
C   *                        NABLA-ALPHA-form, i.e. the ME-parts       *
C   *                        associated with polarization operator are *
C   *                        transformed to NABLA, MEs associated with *
C   *                        current  operator are in ALPHA-form       *
C   *                                                                  *
C   *                                                                  *
C   *     and the matrix element integral functions                    *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_{lam,lam_spin} | Y(E_b) >_r       *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(0),(+)        *
C   *                                                                  *
C   *                                                                  *
C   *  DK, 2013 Dez                                                    *
C   *  DK, 2015 Jan                                                    *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,IKMCPLWF,NCPLWFMAX
      USE MOD_RMESH,ONLY:RWS,JRWS,JRCRI,FULLPOT,NRMAX,R
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL
      IMPLICIT NONE
C*--ME_SPIN_NAB_ALF47
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SPIN_NAB_ALF')
      INTEGER NPOL
      PARAMETER (NPOL=3)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFILA,IFILB,IT
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3,3),MIRR_3(NKMMAX,NKMMAX,3,3,3)
     &           ,MIRR_4(NKMMAX,NKMMAX,3,3,3),MZAZB(NKMMAX,NKMMAX,3,3),
     &           MZBZA(NKMMAX,NKMMAX,3,3)
C
C Local variables
C
      INTEGER IA1,IA_ERR,IB4,IFIL_LHSA,IFIL_LHSB,IFIL_RHSA,IFIL_RHSB,IM,
     &        IRTOP,JPOL,LAMA1,LAMA2,LAMB3,LAMB4
      COMPLEX*16 IJAJB(:),IJAZB(:),IMC,IZAJB(:),IZAZB(:),JFLA(:,:,:),
     &           JFLB(:,:,:),JFLBP(:,:,:),JFRA(:,:,:),JFRAP(:,:,:),
     &           JFRB(:,:,:),JGLA(:,:,:),JGLB(:,:,:),JGLBP(:,:,:),
     &           JGRA(:,:,:),JGRAP(:,:,:),JGRB(:,:,:),MZBZANAB(:,:,:,:),
     &           PAB,PBA,PF_NAB_FORM,PF_NAB_FORM_ALPHA,SJAJB,SJAZB,
     &           SZAJB,SZAZB,ZFLA(:,:,:),ZFLB(:,:,:),ZFLBP(:,:,:),
     &           ZFRA(:,:,:),ZFRAP(:,:,:),ZFRB(:,:,:),ZGLA(:,:,:),
     &           ZGLB(:,:,:),ZGLBP(:,:,:),ZGRA(:,:,:),ZGRAP(:,:,:),
     &           ZGRB(:,:,:)
      LOGICAL K_SELECTION_RULES
      REAL*8 RTOP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZGLA,ZFLA,ZGRA,ZFRA,ZGRAP,ZFRAP
      ALLOCATABLE JGLA,JFLA,JGRA,JFRA,JGRAP,JFRAP
      ALLOCATABLE ZGRB,ZFRB,ZGLB,ZFLB,ZGLBP,ZFLBP
      ALLOCATABLE JGRB,JFRB,JGLB,JFLB,JGLBP,JFLBP
      ALLOCATABLE IZAZB,IZAJB,IJAZB,IJAJB
      ALLOCATABLE MZBZANAB
C
      ALLOCATE (JFLA(NRMAX,NCPLWFMAX,NKM),JFRB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGLA(NRMAX,NCPLWFMAX,NKM),JGRB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFLA(NRMAX,NCPLWFMAX,NKM),ZFRB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,NKM),ZGRB(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKM),ZGRAP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFRA(NRMAX,NCPLWFMAX,NKM),ZFRAP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGRA(NRMAX,NCPLWFMAX,NKM),JGRAP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKM),JFRAP(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (ZGLB(NRMAX,NCPLWFMAX,NKM),ZGLBP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,NKM),ZFLBP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGLB(NRMAX,NCPLWFMAX,NKM),JGLBP(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JFLB(NRMAX,NCPLWFMAX,NKM),JFLBP(NRMAX,NCPLWFMAX,NKM),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ZGLB')
C
      ALLOCATE (MZBZANAB(NKMMAX,NKMMAX,3,3))
C
      ALLOCATE (IZAZB(NRMAX))
      ALLOCATE (IZAJB(NRMAX))
      ALLOCATE (IJAZB(NRMAX))
      ALLOCATE (IJAJB(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IZAZB')
C
C=======================================================================
C ------------------------- set the chanel number for the wave functions
C
      IFIL_RHSA = IFILA
      IFIL_RHSB = IFILB
C
      IF ( LHS_SOL_EQ_RHS_SOL ) THEN
         IFIL_LHSA = IFILA
         IFIL_LHSB = IFILB
      ELSE
         IFIL_LHSA = IFIL_RHSA + 1
         IFIL_LHSB = IFIL_RHSB + 1
      END IF
C
C=======================================================================
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
         RTOP = R(IRTOP,IM)
      ELSE
         IRTOP = JRWS(IM)
         RTOP = RWS(IM)
      END IF
C
      IMC = CI*0.5D0*C
C
      PAB = (C/(C**2+ERYDA+ERYDB))*RTOP**2
      PBA = (C/(C**2+ERYDA+ERYDB))*RTOP**2
C
C     TODO : check whether DCONJG at ERYDB is needed t!!!
      PF_NAB_FORM = 1.0D0/(1.0D0+(ERYDB+ERYDA)/(2*0.5D0*C**2))
C      PF_NAB_FORM = 1.0D0/(1.0D0+(DCONJG(ERYDB)+ERYDA)/(2*0.5D0*C**2))
C      PF_NAB_FORM = C0
C
      PF_NAB_FORM_ALPHA = PF_NAB_FORM/(0.5D0*C**2)
C
C=======================================================================
C
      MZAZB(1:NKMMAX,1:NKMMAX,1:3,1:3) = C0
      MZBZA(1:NKMMAX,1:NKMMAX,1:3,1:3) = C0
      MIRR_2(1:NKMMAX,1:NKMMAX,1:3,1:3,1:3) = C0
      MIRR_3(1:NKMMAX,1:NKMMAX,1:3,1:3,1:3) = C0
      MIRR_4(1:NKMMAX,1:NKMMAX,1:3,1:3,1:3) = C0
C
      MZBZANAB(1:NKMMAX,1:NKMMAX,1:3,1:3) = C0
C
C
C=======================================================================
C
C      evaluate the matrix elements of the type
C
C               MZBZA = < Z^+(E_b) | H_{lam,lam_spin} | Z(E_a) >
C
C=======================================================================
C
C ------------------------------------------ read in wavefunctions for B
C
      CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGLB,ZFLB,JGLB,JFLB,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
      DO LAMB4 = 1,NKM
         DO IB4 = 1,NCPLWF(LAMB4)
            CALL CDIFFER(IM,ZGLB(1,IB4,LAMB4),ZGLBP(1,IB4,LAMB4))
            CALL CDIFFER(IM,ZFLB(1,IB4,LAMB4),ZFLBP(1,IB4,LAMB4))
            CALL CDIFFER(IM,JGLB(1,IB4,LAMB4),JGLBP(1,IB4,LAMB4))
            CALL CDIFFER(IM,JFLB(1,IB4,LAMB4),JFLBP(1,IB4,LAMB4))
         END DO
      END DO
C
C ------------------------------------------ read in wavefunctions for A
C
      CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
      DO LAMA1 = 1,NKM
         DO IA1 = 1,NCPLWF(LAMA1)
            CALL CDIFFER(IM,ZGRA(1,IA1,LAMA1),ZGRAP(1,IA1,LAMA1))
            CALL CDIFFER(IM,ZFRA(1,IA1,LAMA1),ZFRAP(1,IA1,LAMA1))
            CALL CDIFFER(IM,JGRA(1,IA1,LAMA1),JGRAP(1,IA1,LAMA1))
            CALL CDIFFER(IM,JFRA(1,IA1,LAMA1),JFRAP(1,IA1,LAMA1))
         END DO
      END DO
C
      CALL ME_SPIN_CURR_NAB_B4A1(IT,IM,IMC,PF_NAB_FORM,
     &                           PF_NAB_FORM_ALPHA,PBA,ZGLB,ZFLB,ZGRA,
     &                           ZFRA,ZGLBP,ZFLBP,ZGRAP,ZFRAP,MZBZA,
     &                           MZBZANAB)
C
CC---------------------------------------- for nabla related surface term
CC-------------------------------------------------------------------  ZZ
C      CALL ME_SPIN_CURR_NAB_BXAY_NAB(IM,ZGLB,ZFLB,ZGRA,ZFRA,ZGLBP,ZFLBP,
C     &                               ZGRAP,ZFRAP,MZBZANAB)
CC
C=======================================================================
C
C      evaluate the matrix elements involving irregular solutions
C
C=======================================================================
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C ------------------------------------------ read in wavefunctions for A
C                                                use arrays *B as buffer
C
      CALL WAVFUN_READ_REL(IFIL_LHSA,IT,1,ZGRB,ZFRB,JGRB,JFRB,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
      CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZGRB,ZFRB,JGRB,JFRB,ZGLA,ZFLA,JGLA,
     &                            JFLA,NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for B
C
      CALL WAVFUN_READ_REL(IFIL_RHSB,IT,1,ZGRB,ZFRB,JGRB,JFRB,IRTOP,
     &                     NCPLWF,IKMCPLWF)
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
C
      IF ( NKM.GE.0 ) THEN
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
         DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
            DO LAMA2 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
               DO JPOL = 1,NPOL
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
                  CALL ME_ALF_A2B3(LAMA2,LAMB3,JPOL,IM,IMC,PAB,ZGLA,
     &                             ZFLA,JGLA,JFLA,ZGRB,ZFRB,JGRB,JFRB,
     &                             IZAZB,IZAJB,IJAZB,IJAJB,SZAZB,SZAJB,
     &                             SJAZB,SJAJB,MZAZB,K_SELECTION_RULES)
C
C#######################################################################
C
                  IF ( K_SELECTION_RULES ) THEN
C
C***********************************************************************
C***********************************************************************
C     TERM 2    - TAU(3,4,b) * MIRR_2(3,4)                 LAMA1 = LAMA2
C***********************************************************************
C***********************************************************************
C
                     LAMA1 = LAMA2
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
                     DO LAMB4 = 1,NKM
C
                        CALL ME_SPIN_CURR_NAB_BXAXW(IT,IM,IMC,
     &                     PF_NAB_FORM,PF_NAB_FORM_ALPHA,PBA,LAMB4,ZGLB,
     &                     ZFLB,ZGLB,ZFLB,LAMA1,ZGRA,ZFRA,JGRA,JFRA,
     &                     ZGLBP,ZFLBP,ZGLBP,ZFLBP,ZGRAP,ZFRAP,JGRAP,
     &                     JFRAP,IZAZB,SZAZB,IJAZB,SJAZB,MIRR_2,LAMB4,
     &                     LAMB3,JPOL,MZBZANAB)
C
                     END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB LAMB4
C
C
C***********************************************************************
C***********************************************************************
C     TERM 3    - TAU(1,2,a) * MIRR_3(1,2)                 LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
                     DO LAMA1 = 1,NKM
C
                        CALL ME_SPIN_CURR_NAB_BXAXW(IT,IM,IMC,
     &                     PF_NAB_FORM,PF_NAB_FORM_ALPHA,PBA,LAMB3,ZGLB,
     &                     ZFLB,JGLB,JFLB,LAMA1,ZGRA,ZFRA,ZGRA,ZFRA,
     &                     ZGLBP,ZFLBP,JGLBP,JFLBP,ZGRAP,ZFRAP,ZGRAP,
     &                     ZFRAP,IZAZB,SZAZB,IZAJB,SZAJB,MIRR_3,LAMA2,
     &                     LAMA1,JPOL,MZBZANAB)
C
                     END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA LAMA1
C
C
C***********************************************************************
C***********************************************************************
C     TERM 4    MIRR_4(1,3)           LAMA1 = LAMA2  and   LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
                     LAMA1 = LAMA2
C
                     CALL ME_SPIN_CURR_NAB_BXAXW(IT,IM,IMC,PF_NAB_FORM,
     &                  PF_NAB_FORM_ALPHA,PBA,LAMB3,ZGLB,ZFLB,JGLB,JFLB,
     &                  LAMA1,ZGRA,ZFRA,JGRA,JFRA,ZGLBP,ZFLBP,JGLBP,
     &                  JFLBP,ZGRAP,ZFRAP,JGRAP,JFRAP,IZAZB,SZAZB,IJAJB,
     &                  SJAJB,MIRR_4,LAMA1,LAMB3,JPOL,MZBZANAB)
C
C***********************************************************************
C
                  END IF
C################################################################ IFLAG1
               END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
            END DO
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
         END DO
      END IF
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C-----------------------------------------------------------------------
C    MZAZB refer to the standard electric current operator
C    the spin index IPOL_SPIN is therefore a dummy index
C    use the result for IPOL_SPIN=1 throughout
C-----------------------------------------------------------------------
C
      MZAZB(:,:,:,2) = MZAZB(:,:,:,1)
      MZAZB(:,:,:,3) = MZAZB(:,:,:,1)
C
      END
C*==me_spin_curr_nab_b4a1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SPIN_CURR_NAB_B4A1(IT,IM,IMC,PF_NAB_FORM,
     &                                 PF_NAB_FORM_ALPHA,PBA,ZGLB,ZFLB,
     &                                 ZGRA,ZFRA,ZGLBP,ZFLBP,ZGRAP,
     &                                 ZFRAP,MZBZA,MZBZANAB)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements of the type                        *
C   *                                                                  *
C   *           MZBZA = < Z^+(E_b) | H_{lam,lam_spin} | Z(E_a) >       *
C   *                                                                  *
C   *     H_{lam,lam_spin}   spin current operator in the ALFA-form    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF,NCPLWFMAX,VT
      USE MOD_RMESH,ONLY:R,R2DRDI,JRWS,JRCRI,NRMAX,FULLPOT
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_2
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,A1_SPIN_CURR_ALF,A2_SPIN_CURR_ALF,
     &    NCPLWF,L_IKM,LB_IKM,AF1_SPIN_CURR_NAB,AF2_SPIN_CURR_NAB,
     &    AG1_SPIN_CURR_NAB,AG2_SPIN_CURR_NAB
      USE MOD_SIG,ONLY:LSP_NO_ALPHA,LSP_NO_NABLA
      IMPLICIT NONE
C*--ME_SPIN_CURR_NAB_B4A1386
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SPIN_CURR_NAB_B4A1')
      INTEGER NPOL
      PARAMETER (NPOL=3)
C
C Dummy arguments
C
      INTEGER IM,IT
      COMPLEX*16 IMC,PBA,PF_NAB_FORM,PF_NAB_FORM_ALPHA
      COMPLEX*16 MZBZA(NKMMAX,NKMMAX,3,3),MZBZANAB(NKMMAX,NKMMAX,3,3),
     &           ZFLB(NRMAX,NCPLWFMAX,NKM),ZFLBP(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZFRAP(NRMAX,NCPLWFMAX,NKM),
     &           ZGLB(NRMAX,NCPLWFMAX,NKM),ZGLBP(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM),ZGRAP(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 AF1_B4A1,AF2_B4A1,AF3_A1B4,AF4_A1B4,AG1_B4A1,AG2_B4A1,
     &           AG3_A1B4,AG4_A1B4,CADD,CINT1(:),CINT2(:),CINT3(:),
     &           CINT4(:),FAG,FZBZA(:),GAF,IMCA1,IMCA2,ME_D(3,3),
     &           ME_PM0(3,3),ME_P_XYZ(3,3),PA,PFF,PFG,PRE_NABLA,R2_ZFLB,
     &           R2_ZFRA,R2_ZGLB,R2_ZGRA,RMEF1(2,2),RMEF2(2,2),
     &           RMEF3(2,2),RMEF4(2,2),RMEG1(2,2),RMEG2(2,2),RMEG3(2,2),
     &           RMEG4(2,2),S2BAR,TERM12,TERM34,V2BAR
      LOGICAL CNON0
      INTEGER IA1,IA_ERR,IB4,IFLAG1,IKMA1,IKMB4,IPOL,IPOL_SPIN,IR,IRTOP,
     &        LA1,LAMA1,LAMB4,LB4,LBA1,LBB4
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZBZA,CINT1,CINT2,CINT3,CINT4
C
      ALLOCATE (CINT3(NRMAX),CINT1(NRMAX),FZBZA(NRMAX))
      ALLOCATE (CINT4(NRMAX),CINT2(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate FZBZA')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      PRE_NABLA = 1.0D0/(2.0D0*0.5D0)*(1.0D0/CI)*0.5D0
C
C---------------------combine (multiply) prefactor (PF) and (\beta-1/PF)
C
      PFG = 1D0 - PF_NAB_FORM
      PFF = 1D0 + PF_NAB_FORM
C
C      IF ( LSP_NO_NABLA ) PRE_NABLA = 0D0
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
      DO LAMB4 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
         DO LAMA1 = 1,NKM
C
C-----------------------------------------------------------------------
C       BEGIN                alpha - term transformed to nabla form
C                                    contribution in
C                                                     ******************
C =================================================== *  V/c alfa . a  *
C                                                     ******************
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
            DO IPOL_SPIN = 1,NPOL
C               IF ( IPOL_SPIN.NE.3 ) CYCLE
C
               DO IPOL = 1,NPOL
C
                  DO IR = 1,IRTOP
                     FZBZA(IR) = C0
                  END DO
                  S2BAR = C0
C
                  IFLAG1 = 0
C
                  DO IB4 = 1,NCPLWF(LAMB4)
                     IKMB4 = IKMCPLWF(IB4,LAMB4)
C
                     DO IA1 = 1,NCPLWF(LAMA1)
                        IKMA1 = IKMCPLWF(IA1,LAMA1)
C
                        IF ( CNON0(A1_SPIN_CURR_ALF(IKMB4,IKMA1,IPOL,
     &                       IPOL_SPIN)) .OR. 
     &                       CNON0(A2_SPIN_CURR_ALF(IKMB4,IKMA1,IPOL,
     &                       IPOL_SPIN)) ) THEN
C
                           IFLAG1 = 1
C
                           IMCA1 = IMC*A1_SPIN_CURR_ALF(IKMB4,IKMA1,
     &                             IPOL,IPOL_SPIN)
                           IMCA2 = IMC*A2_SPIN_CURR_ALF(IKMB4,IKMA1,
     &                             IPOL,IPOL_SPIN)
C
                           DO IR = 1,IRTOP
C
                              GAF = IMCA1*ZGLB(IR,IB4,LAMB4)
     &                              *ZFRA(IR,IA1,LAMA1)
                              FAG = IMCA2*ZFLB(IR,IB4,LAMB4)
     &                              *ZGRA(IR,IA1,LAMA1)
C
CDK switch
C                                              NB!!!!  additional factor
                              FZBZA(IR) = FZBZA(IR) + (GAF-FAG)
     &                           *R2DRDI(IR,IM)*VT(IR,IT)
     &                           *PF_NAB_FORM_ALPHA
C
                           END DO
C
CDK switch surface term
C                                              NB!!!!  additional factor
                           IR = IRTOP
                           PA = -0.5D0*PBA*(IMCA1-IMCA2)*VT(IR,IT)
     &                          *PF_NAB_FORM_ALPHA
                           S2BAR = S2BAR + 
     &                             (ZGLB(IR,IB4,LAMB4)*ZGRA(IR,IA1,
     &                             LAMA1)-ZFLB(IR,IB4,LAMB4)
     &                             *ZFRA(IR,IA1,LAMA1))*PA
C                           s2bar = c0
C
                        END IF
                     END DO
                  END DO
C
C#######################################################################
C
                  IF ( IFLAG1.EQ.1 ) THEN
C
                     CALL CRADINT(IM,FZBZA,V2BAR)
C
                     IF ( .NOT.LSP_NO_ALPHA )
     &                    MZBZA(LAMB4,LAMA1,IPOL,IPOL_SPIN)
     &                    = MZBZA(LAMB4,LAMA1,IPOL,IPOL_SPIN) + V2BAR + 
     &                    S2BAR
C
                  END IF
C################################################################ IFLAG1
               END DO
            END DO
C-----------------------------------------------------------------------
C       END           Alpha - term transformed to nabla form (V/c alpha)
C-----------------------------------------------------------------------
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
C-----------------------------------------------------------------------
C       BEGIN     small Nabla - related term of BW polarisation operator
C-----------------------------------------------------------------------
C
C ------------------------------------- calculate radial matrix elements
C
            DO IB4 = 1,NCPLWF(LAMB4)
               IKMB4 = IKMCPLWF(IB4,LAMB4)
               LB4 = L_IKM(IKMB4)
               LBB4 = LB_IKM(IKMB4)
C
               DO IA1 = 1,NCPLWF(LAMA1)
                  IKMA1 = IKMCPLWF(IA1,LAMA1)
                  LA1 = L_IKM(IKMA1)
                  LBA1 = LB_IKM(IKMA1)
C
                  DO IR = 1,IRTOP
C
                     R2_ZGRA = ZGRA(IR,IA1,LAMA1)*R2DRDI(IR,IM)
                     CINT3(IR) = R2_ZGRA*(ZGLBP(IR,IB4,LAMB4)-ZGLB(IR,
     &                           IB4,LAMB4)*LB4/R(IR,IM))
                     CINT4(IR) = R2_ZGRA*(ZGLBP(IR,IB4,LAMB4)+ZGLB(IR,
     &                           IB4,LAMB4)*(LB4+1)/R(IR,IM))
C
                     R2_ZGLB = ZGLB(IR,IB4,LAMB4)*R2DRDI(IR,IM)
                     CINT1(IR) = R2_ZGLB*(ZGRAP(IR,IA1,LAMA1)-ZGRA(IR,
     &                           IA1,LAMA1)*LA1/R(IR,IM))
                     CINT2(IR) = R2_ZGLB*(ZGRAP(IR,IA1,LAMA1)+ZGRA(IR,
     &                           IA1,LAMA1)*(LA1+1)/R(IR,IM))
C
                  END DO
C
                  CALL CRADINT(IM,CINT1,RMEG1(IA1,IB4))
                  CALL CRADINT(IM,CINT2,RMEG2(IA1,IB4))
                  CALL CRADINT(IM,CINT3,RMEG3(IA1,IB4))
                  CALL CRADINT(IM,CINT4,RMEG4(IA1,IB4))
C
C
                  DO IR = 1,IRTOP
C
                     R2_ZFRA = ZFRA(IR,IA1,LAMA1)*R2DRDI(IR,IM)
                     CINT3(IR) = R2_ZFRA*(ZFLBP(IR,IB4,LAMB4)-ZFLB(IR,
     &                           IB4,LAMB4)*LBB4/R(IR,IM))
                     CINT4(IR) = R2_ZFRA*(ZFLBP(IR,IB4,LAMB4)+ZFLB(IR,
     &                           IB4,LAMB4)*(LBB4+1)/R(IR,IM))
C
                     R2_ZFLB = ZFLB(IR,IB4,LAMB4)*R2DRDI(IR,IM)
                     CINT1(IR) = R2_ZFLB*(ZFRAP(IR,IA1,LAMA1)-ZFRA(IR,
     &                           IA1,LAMA1)*LBA1/R(IR,IM))
                     CINT2(IR) = R2_ZFLB*(ZFRAP(IR,IA1,LAMA1)+ZFRA(IR,
     &                           IA1,LAMA1)*(LBA1+1)/R(IR,IM))
C
                  END DO
C
                  CALL CRADINT(IM,CINT1,RMEF1(IA1,IB4))
                  CALL CRADINT(IM,CINT2,RMEF2(IA1,IB4))
                  CALL CRADINT(IM,CINT3,RMEF3(IA1,IB4))
                  CALL CRADINT(IM,CINT4,RMEF4(IA1,IB4))
C
C -------------------------------------- calculate total matrix elements
C
                  DO IPOL_SPIN = 1,NPOL
C
                     DO IPOL = 1,NPOL
C
                        AG1_B4A1 = AG1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFG
                        AG2_B4A1 = AG2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFG
                        AF1_B4A1 = AF1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFF
                        AF2_B4A1 = AF2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFF
C
C                           AF1_B4A1 = 0d0
C                           AF2_B4A1 = 0d0
C
                        TERM12 = RMEG1(IA1,IB4)*AG1_B4A1 - 
     &                           RMEG2(IA1,IB4)*AG2_B4A1 + 
     &                           RMEF1(IA1,IB4)*AF1_B4A1 - 
     &                           RMEF2(IA1,IB4)*AF2_B4A1
C
                        AG3_A1B4 = AG2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFG
                        AG4_A1B4 = AG1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFG
                        AF3_A1B4 = AF2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFF
                        AF4_A1B4 = AF1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)*PFF
C
                        TERM34 = RMEG3(IA1,IB4)*AG3_A1B4 - 
     &                           RMEG4(IA1,IB4)*AG4_A1B4 + 
     &                           RMEF3(IA1,IB4)*AF3_A1B4 - 
     &                           RMEF4(IA1,IB4)*AF4_A1B4
C
                        CADD = TERM12 + TERM34
CDEBUG
C                       cadd = c0
C
                        ME_PM0(IPOL,IPOL_SPIN) = PRE_NABLA*CADD
C
                     END DO
C
                  END DO
C
C------------ transform to cartesian basis within the spinpol-component
C
                  ME_P_XYZ(1:3,1) = (-ME_PM0(1:3,3)+ME_PM0(1:3,1))
     &                              /SQRT_2
C
                  ME_P_XYZ(1:3,2) = (ME_PM0(1:3,3)+ME_PM0(1:3,1))
     &                              *CI/SQRT_2
C
                  ME_P_XYZ(1:3,3) = ME_PM0(1:3,2)
C
C
                  IF ( .NOT.LSP_NO_NABLA ) THEN
                     MZBZA(LAMB4,LAMA1,1:3,1:3)
     &                  = MZBZA(LAMB4,LAMA1,1:3,1:3) + ME_P_XYZ(1:3,1:3)
C
C
                     MZBZANAB(LAMB4,LAMA1,1:3,1:3)
     &                  = MZBZANAB(LAMB4,LAMA1,1:3,1:3)
     &                  + ME_P_XYZ(1:3,1:3)
                  END IF
C
C
               END DO
            END DO
C-----------------------------------------------------------------------
C       END     small  Nabla - related term of BW polarisatioon operator
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C       BEGIN                Alpha - term transformed to nabla form
C                                    contribution in
C                                                     ******************
C---------------------------------------------------- *      p . a     *
C                                                     ******************
C
C ------------------------------------- calculate radial matrix elements
C
            DO IB4 = 1,NCPLWF(LAMB4)
               IKMB4 = IKMCPLWF(IB4,LAMB4)
               LB4 = L_IKM(IKMB4)
               LBB4 = LB_IKM(IKMB4)
C
               DO IA1 = 1,NCPLWF(LAMA1)
                  IKMA1 = IKMCPLWF(IA1,LAMA1)
                  LA1 = L_IKM(IKMA1)
                  LBA1 = LB_IKM(IKMA1)
C
                  DO IR = 1,IRTOP
C
                     R2_ZGRA = ZGRA(IR,IA1,LAMA1)*R2DRDI(IR,IM)
                     CINT3(IR) = R2_ZGRA*(ZGLBP(IR,IB4,LAMB4)-ZGLB(IR,
     &                           IB4,LAMB4)*LB4/R(IR,IM))
                     CINT4(IR) = R2_ZGRA*(ZGLBP(IR,IB4,LAMB4)+ZGLB(IR,
     &                           IB4,LAMB4)*(LB4+1)/R(IR,IM))
C
                     R2_ZGLB = ZGLB(IR,IB4,LAMB4)*R2DRDI(IR,IM)
                     CINT1(IR) = R2_ZGLB*(ZGRAP(IR,IA1,LAMA1)-ZGRA(IR,
     &                           IA1,LAMA1)*LA1/R(IR,IM))
                     CINT2(IR) = R2_ZGLB*(ZGRAP(IR,IA1,LAMA1)+ZGRA(IR,
     &                           IA1,LAMA1)*(LA1+1)/R(IR,IM))
C
                  END DO
C
                  CALL CRADINT(IM,CINT1,RMEG1(IA1,IB4))
                  CALL CRADINT(IM,CINT2,RMEG2(IA1,IB4))
                  CALL CRADINT(IM,CINT3,RMEG3(IA1,IB4))
                  CALL CRADINT(IM,CINT4,RMEG4(IA1,IB4))
C
                  DO IR = 1,IRTOP
C
                     R2_ZFRA = ZFRA(IR,IA1,LAMA1)*R2DRDI(IR,IM)
                     CINT3(IR) = R2_ZFRA*(ZFLBP(IR,IB4,LAMB4)-ZFLB(IR,
     &                           IB4,LAMB4)*LBB4/R(IR,IM))
                     CINT4(IR) = R2_ZFRA*(ZFLBP(IR,IB4,LAMB4)+ZFLB(IR,
     &                           IB4,LAMB4)*(LBB4+1)/R(IR,IM))
C
                     R2_ZFLB = ZFLB(IR,IB4,LAMB4)*R2DRDI(IR,IM)
                     CINT1(IR) = R2_ZFLB*(ZFRAP(IR,IA1,LAMA1)-ZFRA(IR,
     &                           IA1,LAMA1)*LBA1/R(IR,IM))
                     CINT2(IR) = R2_ZFLB*(ZFRAP(IR,IA1,LAMA1)+ZFRA(IR,
     &                           IA1,LAMA1)*(LBA1+1)/R(IR,IM))
C
                  END DO
C
                  CALL CRADINT(IM,CINT1,RMEF1(IA1,IB4))
                  CALL CRADINT(IM,CINT2,RMEF2(IA1,IB4))
                  CALL CRADINT(IM,CINT3,RMEF3(IA1,IB4))
                  CALL CRADINT(IM,CINT4,RMEF4(IA1,IB4))
C
C
C -------------------------------------- calculate total matrix elements
C
                  DO IPOL_SPIN = 1,NPOL
C
                     DO IPOL = 1,NPOL
C
                        AG1_B4A1 = AG1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
                        AG2_B4A1 = AG2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
C                                                NB!!!!!! the minus sign
                        AF1_B4A1 = -AF1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
                        AF2_B4A1 = -AF2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
C
                        TERM12 = RMEG1(IA1,IB4)*AG1_B4A1 - 
     &                           RMEG2(IA1,IB4)*AG2_B4A1 + 
     &                           RMEF1(IA1,IB4)*AF1_B4A1 - 
     &                           RMEF2(IA1,IB4)*AF2_B4A1
C
                        AG3_A1B4 = AG2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
                        AG4_A1B4 = AG1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
C                                                NB!!!!!! the minus sign
                        AF3_A1B4 = -AF2_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
                        AF4_A1B4 = -AF1_SPIN_CURR_NAB(IKMB4,IKMA1,IPOL,
     &                             IPOL_SPIN)
C
                        TERM34 = RMEG3(IA1,IB4)*AG3_A1B4 - 
     &                           RMEG4(IA1,IB4)*AG4_A1B4 + 
     &                           RMEF3(IA1,IB4)*AF3_A1B4 - 
     &                           RMEF4(IA1,IB4)*AF4_A1B4
C
                        CADD = TERM12 + TERM34
C
C                                             NB!!!!!! additional factor
                        ME_PM0(IPOL,IPOL_SPIN)
     &                     = PRE_NABLA*PF_NAB_FORM*CADD
C
                     END DO
C
                  END DO
C
C------------ transform to cartesian basis within the spinpol-component
C
                  ME_P_XYZ(1:3,1) = (-ME_PM0(1:3,3)+ME_PM0(1:3,1))
     &                              /SQRT_2
C
                  ME_P_XYZ(1:3,2) = (ME_PM0(1:3,3)+ME_PM0(1:3,1))
     &                              *CI/SQRT_2
C
                  ME_P_XYZ(1:3,3) = ME_PM0(1:3,2)
C                  CALL  CMATSTR('ME_P_XYZ',ME_P_XYZ,3,3,0,0,0,1D-8,6)
C
C-------------------- transform to cartesian basis within ipol component
                  ME_D(1,1:3) = (-ME_P_XYZ(3,1:3)+ME_P_XYZ(1,1:3))
     &                          /SQRT_2
                  ME_D(2,1:3) = (ME_P_XYZ(3,1:3)+ME_P_XYZ(1,1:3))
     &                          *CI/SQRT_2
                  ME_D(3,1:3) = ME_P_XYZ(2,1:3)
C----------------------------------------------------swap ipol, ipolspin
                  ME_D = TRANSPOSE(ME_D)
C----------------back transform to spherical basis within ipol component
                  ME_P_XYZ(3,1:3) = -(ME_D(1,1:3)+CI*ME_D(2,1:3))/SQRT_2
                  ME_P_XYZ(1,1:3) = (ME_D(1,1:3)-CI*ME_D(2,1:3))/SQRT_2
                  ME_P_XYZ(2,1:3) = ME_D(3,1:3)
C                  CALL  CMATSTR('ME_P_XYZ',ME_P_XYZ,3,3,0,0,0,1D-8,6)
C
                  MZBZA(LAMB4,LAMA1,1:3,1:3)
     &               = MZBZA(LAMB4,LAMA1,1:3,1:3) + ME_P_XYZ(1:3,1:3)
C
C
                  MZBZANAB(LAMB4,LAMA1,1:3,1:3)
     &               = MZBZANAB(LAMB4,LAMA1,1:3,1:3) + ME_P_XYZ(1:3,1:3)
C
               END DO
            END DO
C-----------------------------------------------------------------------
C       END                   Alpha - term transformed to nabla form p.a
C-----------------------------------------------------------------------
C
         END DO
C                                                      energy A -- LAMA1
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB4
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
C*==me_spin_curr_nab_bxaxw.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SPIN_CURR_NAB_BXAXW(IT,IM,IMC,PF_NAB_FORM,
     &                                  PF_NAB_FORM_ALPHA,PBA,LAMB,ZGLB,
     &                                  ZFLB,JGLB,JFLB,LAMA,ZGRA,ZFRA,
     &                                  JGRA,JFRA,ZGLBP,ZFLBP,JGLBP,
     &                                  JFLBP,ZGRAP,ZFRAP,JGRAP,JFRAP,
     &                                  IZAZB,SZAZB,IJAJB,SJAJB,MIRR,I,
     &                                  J,JPOL,MZBZANAB)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements of the type                        *
C   *                                                                  *
C   *           MIRR = < Z^+(E_b) | H_{lam,lam_spin} | Z(E_a) >        *
C   *                                                                  *
C   *     H_{lam,lam_spin}   spin current operator in the ALFA-form    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF,NCPLWFMAX,VT
      USE MOD_RMESH,ONLY:R,R2DRDI,JRWS,JRCRI,NRMAX,FULLPOT
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_2
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,A1_SPIN_CURR_ALF,A2_SPIN_CURR_ALF,
     &    NCPLWF,L_IKM,LB_IKM,AF1_SPIN_CURR_NAB,AF2_SPIN_CURR_NAB,
     &    AG1_SPIN_CURR_NAB,AG2_SPIN_CURR_NAB
      USE MOD_SIG,ONLY:LSP_NO_ALPHA,LSP_NO_NABLA
      IMPLICIT NONE
C*--ME_SPIN_CURR_NAB_BXAXW867
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SPIN_CURR_NAB_BXAXW')
      INTEGER NPOL
      PARAMETER (NPOL=3)
C
C Dummy arguments
C
      INTEGER I,IM,IT,J,JPOL,LAMA,LAMB
      COMPLEX*16 IMC,PBA,PF_NAB_FORM,PF_NAB_FORM_ALPHA,SJAJB,SZAZB
      COMPLEX*16 IJAJB(NRMAX),IZAZB(NRMAX),JFLB(NRMAX,NCPLWFMAX,NKM),
     &           JFLBP(NRMAX,NCPLWFMAX,NKM),JFRA(NRMAX,NCPLWFMAX,NKM),
     &           JFRAP(NRMAX,NCPLWFMAX,NKM),JGLB(NRMAX,NCPLWFMAX,NKM),
     &           JGLBP(NRMAX,NCPLWFMAX,NKM),JGRA(NRMAX,NCPLWFMAX,NKM),
     &           JGRAP(NRMAX,NCPLWFMAX,NKM),MIRR(NKMMAX,NKMMAX,3,3,3),
     &           MZBZANAB(NKMMAX,NKMMAX,3,3),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFLBP(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           ZFRAP(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGLBP(NRMAX,NCPLWFMAX,NKM),ZGRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGRAP(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 A1JFRA,A1ZFRA,A2JGRA,A2ZGRA,AF1_B4A1,AF2_B4A1,AF3_A1B4,
     &           AF4_A1B4,AG1_B4A1,AG2_B4A1,AG3_A1B4,AG4_A1B4,CADD,
     &           CINT1(:),CINT2(:),CINT3(:),CINT4(:),FS(:),FV(:),IMCA1,
     &           IMCA2,ME_D(3,3),ME_PM0(3,3),ME_P_XYZ(3,3),PA,PFF,PFG,
     &           PRE_NABLA,R2_IJAJB,R2_IZAZB,R2_JFLB,R2_JFRA,R2_JGLB,
     &           R2_JGRA,R2_ZFLB,R2_ZFRA,R2_ZGLB,R2_ZGRA,RMEF1(2,2),
     &           RMEF2(2,2),RMEF3(2,2),RMEF4(2,2),RMEG1(2,2),RMEG2(2,2),
     &           RMEG3(2,2),RMEG4(2,2),S2REG,SCNTR,SCNTRNAB(3),TERM12,
     &           TERM34,VCNTR,VIRR,WIRR,WREG
      LOGICAL CNON0
      INTEGER IA,IA_ERR,IB,IFLAG2,IKMA,IKMB,IPOL,IPOL_SPIN,IR,IRTOP,LA,
     &        LB,LBA,LBB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FV,FS,CINT1,CINT2,CINT3,CINT4
C
      ALLOCATE (CINT3(NRMAX),CINT1(NRMAX),FV(NRMAX),FS(NRMAX))
      ALLOCATE (CINT4(NRMAX),CINT2(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate FV')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
         CALL STOP_MESSAGE(ROUTINE,'FULLPOT not yet available')
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      PRE_NABLA = 1.0D0/(2.0D0*0.5D0)*(1.0D0/CI)*0.5D0
C
      PFG = 1D0 - PF_NAB_FORM
      PFF = 1D0 + PF_NAB_FORM
C
C      IF ( LSP_NO_NABLA ) PRE_NABLA = 0D0
C
C-----------------------------------------------------------------------
C       BEGIN                Alpha - term transformed to nabla form
C                                    contribution in
C                                                     ******************
C =================================================== *  V/c alfa . a  *
C                                                     ******************
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      LOOP_IPOL_SPIN_1:DO IPOL_SPIN = 1,NPOL
C
         LOOP_IPOL_1:DO IPOL = 1,NPOL
C
            DO IR = 1,IRTOP
               FV(IR) = C0
               FS(IR) = C0
            END DO
C
            S2REG = C0
C
            IFLAG2 = 0
C
C-----------------------------------------------------------------------
C                                                                IA IB
            DO IB = 1,NCPLWF(LAMB)
               IKMB = IKMCPLWF(IB,LAMB)
C
               DO IA = 1,NCPLWF(LAMA)
                  IKMA = IKMCPLWF(IA,LAMA)
C
                  IF ( CNON0(A1_SPIN_CURR_ALF(IKMB,IKMA,IPOL,IPOL_SPIN))
     &                 .OR. 
     &                 CNON0(A2_SPIN_CURR_ALF(IKMB,IKMA,IPOL,IPOL_SPIN))
     &                 ) THEN
C
                     IFLAG2 = 1
C
                     IMCA1 = IMC*A1_SPIN_CURR_ALF(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
                     IMCA2 = IMC*A2_SPIN_CURR_ALF(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
C
                     DO IR = 1,IRTOP
C
                        A2ZGRA = IMCA2*ZGRA(IR,IA,LAMA)
                        A1ZFRA = IMCA1*ZFRA(IR,IA,LAMA)
                        A2JGRA = IMCA2*JGRA(IR,IA,LAMA)
                        A1JFRA = IMCA1*JFRA(IR,IA,LAMA)
C
                        WREG = JGLB(IR,IB,LAMB)
     &                         *A1JFRA - JFLB(IR,IB,LAMB)*A2JGRA
                        WIRR = ZGLB(IR,IB,LAMB)
     &                         *A1ZFRA - ZFLB(IR,IB,LAMB)*A2ZGRA
C
CDK switch
C                                              NB!!!!  additional factor
                        WREG = WREG*R2DRDI(IR,IM)*VT(IR,IT)
     &                         *PF_NAB_FORM_ALPHA
                        WIRR = WIRR*R2DRDI(IR,IM)*VT(IR,IT)
     &                         *PF_NAB_FORM_ALPHA
C
                        FV(IR) = FV(IR) + WREG*IZAZB(IR)
     &                           + WIRR*IJAJB(IR)
                        FS(IR) = FS(IR) + WIRR
C
                     END DO
C
                     IR = IRTOP
C                                              NB!!!!  additional factor
                     PA = -0.5D0*PBA*(IMCA1-IMCA2)*VT(IR,IT)
     &                    *PF_NAB_FORM_ALPHA
                     S2REG = S2REG + (JGLB(IR,IB,LAMB)*JGRA(IR,IA,LAMA)
     &                       -JFLB(IR,IB,LAMB)*JFRA(IR,IA,LAMA))*PA
C
                  END IF
               END DO
            END DO
C                                                                IA IB
C-----------------------------------------------------------------------
C
C***********************************************************************
            IF ( IFLAG2.EQ.1 ) THEN
C
C------------------------------------------------------ pure volume term
C
               CALL CRADINT(IM,FV,VCNTR)
C
C---------------------------------------- alpha - related   surface term
C
               CALL CRADINT(IM,FS,VIRR)
C
CDK switch surface term
               SCNTR = VIRR*SJAJB + S2REG*(IZAZB(IRTOP)+SZAZB)
C               SCNTR = C0
C
               IF ( .NOT.LSP_NO_ALPHA ) MIRR(I,J,IPOL,JPOL,IPOL_SPIN)
     &              = MIRR(I,J,IPOL,JPOL,IPOL_SPIN) + VCNTR + SCNTR
C
            END IF
C***********************************************************************
C
         END DO LOOP_IPOL_1
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      END DO LOOP_IPOL_SPIN_1
C-----------------------------------------------------------------------
C       END           alpha - term transformed to nabla form (V/c alpha)
C-----------------------------------------------------------------------C
C-----------------------------------------------------------------------
C       BEGIN   small   Nabla - related term of BW polarisation operator
C-----------------------------------------------------------------------
C
C ------------------------------------- calculate radial matrix elements
C
      LOOP_IB_2:DO IB = 1,NCPLWF(LAMB)
         IKMB = IKMCPLWF(IB,LAMB)
         LB = L_IKM(IKMB)
         LBB = LB_IKM(IKMB)
C
         LOOP_IA_2:DO IA = 1,NCPLWF(LAMA)
            IKMA = IKMCPLWF(IA,LAMA)
            LA = L_IKM(IKMA)
            LBA = LB_IKM(IKMA)
C
            DO IR = 1,IRTOP
C
               R2_IZAZB = R2DRDI(IR,IM)*IZAZB(IR)
               R2_IJAJB = R2DRDI(IR,IM)*IJAJB(IR)
C
               R2_JGRA = JGRA(IR,IA,LAMA)*R2_IZAZB
               CINT3(IR) = R2_JGRA*(JGLBP(IR,IB,LAMB)-JGLB(IR,IB,LAMB)
     &                     *LB/R(IR,IM))
               CINT4(IR) = R2_JGRA*(JGLBP(IR,IB,LAMB)+JGLB(IR,IB,LAMB)
     &                     *(LB+1)/R(IR,IM))
C
               R2_JGLB = JGLB(IR,IB,LAMB)*R2_IZAZB
               CINT1(IR) = R2_JGLB*(JGRAP(IR,IA,LAMA)-JGRA(IR,IA,LAMA)
     &                     *LA/R(IR,IM))
               CINT2(IR) = R2_JGLB*(JGRAP(IR,IA,LAMA)+JGRA(IR,IA,LAMA)
     &                     *(LA+1)/R(IR,IM))
C
               R2_ZGRA = ZGRA(IR,IA,LAMA)*R2_IJAJB
               CINT3(IR) = CINT3(IR)
     &                     + R2_ZGRA*(ZGLBP(IR,IB,LAMB)-ZGLB(IR,IB,LAMB)
     &                     *LB/R(IR,IM))
               CINT4(IR) = CINT4(IR)
     &                     + R2_ZGRA*(ZGLBP(IR,IB,LAMB)+ZGLB(IR,IB,LAMB)
     &                     *(LB+1)/R(IR,IM))
C
               R2_ZGLB = ZGLB(IR,IB,LAMB)*R2_IJAJB
               CINT1(IR) = CINT1(IR)
     &                     + R2_ZGLB*(ZGRAP(IR,IA,LAMA)-ZGRA(IR,IA,LAMA)
     &                     *LA/R(IR,IM))
               CINT2(IR) = CINT2(IR)
     &                     + R2_ZGLB*(ZGRAP(IR,IA,LAMA)+ZGRA(IR,IA,LAMA)
     &                     *(LA+1)/R(IR,IM))
C
            END DO
C
            CALL CRADINT(IM,CINT1,RMEG1(IA,IB))
            CALL CRADINT(IM,CINT2,RMEG2(IA,IB))
            CALL CRADINT(IM,CINT3,RMEG3(IA,IB))
            CALL CRADINT(IM,CINT4,RMEG4(IA,IB))
C
            DO IR = 1,IRTOP
C
               R2_IZAZB = R2DRDI(IR,IM)*IZAZB(IR)
               R2_IJAJB = R2DRDI(IR,IM)*IJAJB(IR)
C
               R2_JFRA = JFRA(IR,IA,LAMA)*R2_IZAZB
               CINT3(IR) = R2_JFRA*(JFLBP(IR,IB,LAMB)-JFLB(IR,IB,LAMB)
     &                     *LBB/R(IR,IM))
               CINT4(IR) = R2_JFRA*(JFLBP(IR,IB,LAMB)+JFLB(IR,IB,LAMB)
     &                     *(LBB+1)/R(IR,IM))
C
               R2_JFLB = JFLB(IR,IB,LAMB)*R2_IZAZB
               CINT1(IR) = R2_JFLB*(JFRAP(IR,IA,LAMA)-JFRA(IR,IA,LAMA)
     &                     *LBA/R(IR,IM))
               CINT2(IR) = R2_JFLB*(JFRAP(IR,IA,LAMA)+JFRA(IR,IA,LAMA)
     &                     *(LBA+1)/R(IR,IM))
C
               R2_ZFRA = ZFRA(IR,IA,LAMA)*R2_IJAJB
               CINT3(IR) = CINT3(IR)
     &                     + R2_ZFRA*(ZFLBP(IR,IB,LAMB)-ZFLB(IR,IB,LAMB)
     &                     *LBB/R(IR,IM))
               CINT4(IR) = CINT4(IR)
     &                     + R2_ZFRA*(ZFLBP(IR,IB,LAMB)+ZFLB(IR,IB,LAMB)
     &                     *(LBB+1)/R(IR,IM))
C
               R2_ZFLB = ZFLB(IR,IB,LAMB)*R2_IJAJB
               CINT1(IR) = CINT1(IR)
     &                     + R2_ZFLB*(ZFRAP(IR,IA,LAMA)-ZFRA(IR,IA,LAMA)
     &                     *LBA/R(IR,IM))
               CINT2(IR) = CINT2(IR)
     &                     + R2_ZFLB*(ZFRAP(IR,IA,LAMA)+ZFRA(IR,IA,LAMA)
     &                     *(LBA+1)/R(IR,IM))
C
            END DO
C
            CALL CRADINT(IM,CINT1,RMEF1(IA,IB))
            CALL CRADINT(IM,CINT2,RMEF2(IA,IB))
            CALL CRADINT(IM,CINT3,RMEF3(IA,IB))
            CALL CRADINT(IM,CINT4,RMEF4(IA,IB))
C
C -------------------------------------- calculate total matrix elements
C
            LOOP_IPOL_SPIN_A:DO IPOL_SPIN = 1,NPOL
C
               LOOP_IPOL_A:DO IPOL = 1,NPOL
C
                  AG1_B4A1 = AG1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFG
                  AG2_B4A1 = AG2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFG
                  AF1_B4A1 = AF1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFF
                  AF2_B4A1 = AF2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFF
C
                  TERM12 = RMEG1(IA,IB)*AG1_B4A1 - RMEG2(IA,IB)
     &                     *AG2_B4A1 + RMEF1(IA,IB)*AF1_B4A1 - 
     &                     RMEF2(IA,IB)*AF2_B4A1
C
                  AG3_A1B4 = AG2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFG
                  AG4_A1B4 = AG1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFG
                  AF3_A1B4 = AF2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFF
                  AF4_A1B4 = AF1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
     &                       *PFF
C
                  TERM34 = RMEG3(IA,IB)*AG3_A1B4 - RMEG4(IA,IB)
     &                     *AG4_A1B4 + RMEF3(IA,IB)*AF3_A1B4 - 
     &                     RMEF4(IA,IB)*AF4_A1B4
C
                  CADD = TERM12 + TERM34
C
                  ME_PM0(IPOL,IPOL_SPIN) = PRE_NABLA*CADD
C
               END DO LOOP_IPOL_A
            END DO LOOP_IPOL_SPIN_A
C
C------------ transform to cartesian basis within the spinpol-component
C
            ME_P_XYZ(1:3,1) = (-ME_PM0(1:3,3)+ME_PM0(1:3,1))/SQRT_2
C
            ME_P_XYZ(1:3,2) = (ME_PM0(1:3,3)+ME_PM0(1:3,1))*CI/SQRT_2
C
            ME_P_XYZ(1:3,3) = ME_PM0(1:3,2)
C
            IF ( .NOT.LSP_NO_NABLA ) MIRR(I,J,1:3,JPOL,1:3)
     &           = MIRR(I,J,1:3,JPOL,1:3) + ME_P_XYZ(1:3,1:3)
C
         END DO LOOP_IA_2
      END DO LOOP_IB_2
C-----------------------------------------------------------------------
C       END  small     Nabla - related term of BW polarisatioon operator
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C       BEGIN                Alpha - term transformed to nabla form
C                                    contribution in
C                                                     ******************
C---------------------------------------------------- *      p . a     *
C                                                     ******************
C
      LOOP_IB_3:DO IB = 1,NCPLWF(LAMB)
         IKMB = IKMCPLWF(IB,LAMB)
         LB = L_IKM(IKMB)
         LBB = LB_IKM(IKMB)
C
         LOOP_IA_3:DO IA = 1,NCPLWF(LAMA)
            IKMA = IKMCPLWF(IA,LAMA)
            LA = L_IKM(IKMA)
            LBA = LB_IKM(IKMA)
C
            DO IR = 1,IRTOP
C
               R2_IZAZB = R2DRDI(IR,IM)*IZAZB(IR)
               R2_IJAJB = R2DRDI(IR,IM)*IJAJB(IR)
C
               R2_JGRA = JGRA(IR,IA,LAMA)*R2_IZAZB
               CINT3(IR) = R2_JGRA*(JGLBP(IR,IB,LAMB)-JGLB(IR,IB,LAMB)
     &                     *LB/R(IR,IM))
               CINT4(IR) = R2_JGRA*(JGLBP(IR,IB,LAMB)+JGLB(IR,IB,LAMB)
     &                     *(LB+1)/R(IR,IM))
C
               R2_JGLB = JGLB(IR,IB,LAMB)*R2_IZAZB
               CINT1(IR) = R2_JGLB*(JGRAP(IR,IA,LAMA)-JGRA(IR,IA,LAMA)
     &                     *LA/R(IR,IM))
               CINT2(IR) = R2_JGLB*(JGRAP(IR,IA,LAMA)+JGRA(IR,IA,LAMA)
     &                     *(LA+1)/R(IR,IM))
C
               R2_ZGRA = ZGRA(IR,IA,LAMA)*R2_IJAJB
               CINT3(IR) = CINT3(IR)
     &                     + R2_ZGRA*(ZGLBP(IR,IB,LAMB)-ZGLB(IR,IB,LAMB)
     &                     *LB/R(IR,IM))
               CINT4(IR) = CINT4(IR)
     &                     + R2_ZGRA*(ZGLBP(IR,IB,LAMB)+ZGLB(IR,IB,LAMB)
     &                     *(LB+1)/R(IR,IM))
C
               R2_ZGLB = ZGLB(IR,IB,LAMB)*R2_IJAJB
               CINT1(IR) = CINT1(IR)
     &                     + R2_ZGLB*(ZGRAP(IR,IA,LAMA)-ZGRA(IR,IA,LAMA)
     &                     *LA/R(IR,IM))
               CINT2(IR) = CINT2(IR)
     &                     + R2_ZGLB*(ZGRAP(IR,IA,LAMA)+ZGRA(IR,IA,LAMA)
     &                     *(LA+1)/R(IR,IM))
C
            END DO
C
            CALL CRADINT(IM,CINT1,RMEG1(IA,IB))
            CALL CRADINT(IM,CINT2,RMEG2(IA,IB))
            CALL CRADINT(IM,CINT3,RMEG3(IA,IB))
            CALL CRADINT(IM,CINT4,RMEG4(IA,IB))
C
            DO IR = 1,IRTOP
C
               R2_IZAZB = R2DRDI(IR,IM)*IZAZB(IR)
               R2_IJAJB = R2DRDI(IR,IM)*IJAJB(IR)
C
               R2_JFRA = JFRA(IR,IA,LAMA)*R2_IZAZB
               CINT3(IR) = R2_JFRA*(JFLBP(IR,IB,LAMB)-JFLB(IR,IB,LAMB)
     &                     *LBB/R(IR,IM))
               CINT4(IR) = R2_JFRA*(JFLBP(IR,IB,LAMB)+JFLB(IR,IB,LAMB)
     &                     *(LBB+1)/R(IR,IM))
C
               R2_JFLB = JFLB(IR,IB,LAMB)*R2_IZAZB
               CINT1(IR) = R2_JFLB*(JFRAP(IR,IA,LAMA)-JFRA(IR,IA,LAMA)
     &                     *LBA/R(IR,IM))
               CINT2(IR) = R2_JFLB*(JFRAP(IR,IA,LAMA)+JFRA(IR,IA,LAMA)
     &                     *(LBA+1)/R(IR,IM))
C
               R2_ZFRA = ZFRA(IR,IA,LAMA)*R2_IJAJB
               CINT3(IR) = CINT3(IR)
     &                     + R2_ZFRA*(ZFLBP(IR,IB,LAMB)-ZFLB(IR,IB,LAMB)
     &                     *LBB/R(IR,IM))
               CINT4(IR) = CINT4(IR)
     &                     + R2_ZFRA*(ZFLBP(IR,IB,LAMB)+ZFLB(IR,IB,LAMB)
     &                     *(LBB+1)/R(IR,IM))
C
               R2_ZFLB = ZFLB(IR,IB,LAMB)*R2_IJAJB
               CINT1(IR) = CINT1(IR)
     &                     + R2_ZFLB*(ZFRAP(IR,IA,LAMA)-ZFRA(IR,IA,LAMA)
     &                     *LBA/R(IR,IM))
               CINT2(IR) = CINT2(IR)
     &                     + R2_ZFLB*(ZFRAP(IR,IA,LAMA)+ZFRA(IR,IA,LAMA)
     &                     *(LBA+1)/R(IR,IM))
C
            END DO
C
            CALL CRADINT(IM,CINT1,RMEF1(IA,IB))
            CALL CRADINT(IM,CINT2,RMEF2(IA,IB))
            CALL CRADINT(IM,CINT3,RMEF3(IA,IB))
            CALL CRADINT(IM,CINT4,RMEF4(IA,IB))
C
C -------------------------------------- calculate total matrix elements
C
            LOOP_IPOL_SPIN_B:DO IPOL_SPIN = 1,NPOL
C
               LOOP_IPOL_B:DO IPOL = 1,NPOL
C
                  AG1_B4A1 = AG1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
                  AG2_B4A1 = AG2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
C                                                NB!!!!!! the minus sign
                  AF1_B4A1 = -AF1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
                  AF2_B4A1 = -AF2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
C
                  TERM12 = RMEG1(IA,IB)*AG1_B4A1 - RMEG2(IA,IB)
     &                     *AG2_B4A1 + RMEF1(IA,IB)*AF1_B4A1 - 
     &                     RMEF2(IA,IB)*AF2_B4A1
C
                  AG3_A1B4 = AG2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
                  AG4_A1B4 = AG1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,IPOL_SPIN)
C                                                NB!!!!!! the minus sign
                  AF3_A1B4 = -AF2_SPIN_CURR_NAB(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
                  AF4_A1B4 = -AF1_SPIN_CURR_NAB(IKMB,IKMA,IPOL,
     &                       IPOL_SPIN)
C
                  TERM34 = RMEG3(IA,IB)*AG3_A1B4 - RMEG4(IA,IB)
     &                     *AG4_A1B4 + RMEF3(IA,IB)*AF3_A1B4 - 
     &                     RMEF4(IA,IB)*AF4_A1B4
C
                  CADD = TERM12 + TERM34
C
C                                             NB!!!!!! additional factor
                  ME_PM0(IPOL,IPOL_SPIN) = PRE_NABLA*PF_NAB_FORM*CADD
C
               END DO LOOP_IPOL_B
            END DO LOOP_IPOL_SPIN_B
C
C------------ transform to cartesian basis within the spinpol-component
C
            ME_P_XYZ(1:3,1) = (-ME_PM0(1:3,3)+ME_PM0(1:3,1))/SQRT_2
C
            ME_P_XYZ(1:3,2) = (ME_PM0(1:3,3)+ME_PM0(1:3,1))*CI/SQRT_2
C
            ME_P_XYZ(1:3,3) = ME_PM0(1:3,2)
C
C-------------------- transform to cartesian basis within ipol component
            ME_D(1,1:3) = (-ME_P_XYZ(3,1:3)+ME_P_XYZ(1,1:3))/SQRT_2
            ME_D(2,1:3) = (ME_P_XYZ(3,1:3)+ME_P_XYZ(1,1:3))*CI/SQRT_2
            ME_D(3,1:3) = ME_P_XYZ(2,1:3)
C----------------------------------------------------swap ipol, ipolspin
            ME_D = TRANSPOSE(ME_D)
C----------------back transform to spherical basis within ipol component
            ME_P_XYZ(3,1:3) = -(ME_D(1,1:3)+CI*ME_D(2,1:3))/SQRT_2
            ME_P_XYZ(1,1:3) = (ME_D(1,1:3)-CI*ME_D(2,1:3))/SQRT_2
            ME_P_XYZ(2,1:3) = ME_D(3,1:3)
C
C
            MIRR(I,J,1:3,JPOL,1:3) = MIRR(I,J,1:3,JPOL,1:3)
     &                               + ME_P_XYZ(1:3,1:3)
C
         END DO LOOP_IA_3
      END DO LOOP_IB_3
C
C-----------------------------------------------------------------------
C       END                   Alpha - term transformed to nabla form p.a
C-----------------------------------------------------------------------
C
C---------------------------------------- add nabla related surface term
C
      DO IPOL_SPIN = 1,NPOL
C
         SCNTRNAB(1:3) = MZBZANAB(LAMB,LAMA,1:3,IPOL_SPIN)*SJAJB
C         SCNTRNAB(1:3) = C0
C
         MIRR(I,J,1:3,JPOL,IPOL_SPIN) = MIRR(I,J,1:3,JPOL,IPOL_SPIN)
     &                                  + SCNTRNAB(1:3)
      END DO
C
      END
