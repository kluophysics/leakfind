C*==me_sig_alf.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SIG_ALF(IFILA,ERYDA,IFILB,ERYDB,IT,MZAZB,MZBZA,
     &                      MIRR_2,MIRR_3,MIRR_4,C)
C   ********************************************************************
C   *                                                                  *
C   *    calculate matrix elements for the spin operator               *
C   *                                                                  *
C   *                       s = sigma                                  *
C   *                                                                  *
C   *    for the calculations related to Edelstein effect              *
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *         MZAZB  =   < Z^+(E_a) | j_lam | Z(E_b) >                 *
C   *                                                                  *
C   *         MZBZA  =   < Z^+(E_b) | s_lam | Z(E_a) >                 *
C   *                                                                  *
C   *         MIRR_2 =   < Z^+(E_b) | s_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAZB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_3 =   < J^+(E_b) | s_lam * IZAZB | Z(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_4 =   < J^+(E_b) | s_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAJB | Z(E_a) >         *
C   *                                                                  *
C   *     with    j_lam = mc alfa_lam = mc alfa * A_lam                *
C   *                                                                  *
C   *     and the matrix element integral functions                    *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | j_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(z),(+)        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,IKMCPLWF,NCPLWFMAX
      USE MOD_CALCMODE,ONLY:LHS_SOL_EQ_RHS_SOL,IREL
      USE MOD_RMESH,ONLY:R,RWS,JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF,NPOL
      IMPLICIT NONE
C*--ME_SIG_ALF45
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SIG_ALF')
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFILA,IFILB,IT
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3),MIRR_3(NKMMAX,NKMMAX,3,3),
     &           MIRR_4(NKMMAX,NKMMAX,3,3),MZAZB(NKMMAX,NKMMAX,3),
     &           MZBZA(NKMMAX,NKMMAX,3)
C
C Local variables
C
      INTEGER IA_ERR,IFIL_LHSA,IFIL_LHSB,IFIL_RHSA,IFIL_RHSB,IM,IRTOP,
     &        JPOL,LAMA1,LAMA2,LAMB3,LAMB4,N
      COMPLEX*16 IJAJB(:),IJAZB(:),IMC,IZAJB(:),IZAZB(:),JFLA(:,:,:),
     &           JFLB(:,:,:),JFRA(:,:,:),JFRB(:,:,:),JGLA(:,:,:),
     &           JGLB(:,:,:),JGRA(:,:,:),JGRB(:,:,:),PAB,SJAJB,SJAZB,
     &           SZAJB,SZAZB,ZFLA(:,:,:),ZFLB(:,:,:),ZFRA(:,:,:),
     &           ZFRB(:,:,:),ZGLA(:,:,:),ZGLB(:,:,:),ZGRA(:,:,:),
     &           ZGRB(:,:,:)
      LOGICAL K_SELECTION_RULES
      REAL*8 RTOP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IZAZB,IZAJB,IJAZB,IJAJB
      ALLOCATABLE JFLA,JGLA,ZGLA,ZFLA,JGRB,JFRB,ZFRB,ZGRB
      ALLOCATABLE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
C
      N = NKM
      ALLOCATE (JFLA(NRMAX,NCPLWFMAX,N),JGLA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,N),ZFLA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (JGRB(NRMAX,NCPLWFMAX,N),JFRB(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZFRB(NRMAX,NCPLWFMAX,N),ZGRB(NRMAX,NCPLWFMAX,N))
C
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,N),JGRA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,N),ZFRA(NRMAX,NCPLWFMAX,N))
      ALLOCATE (JGLB(NRMAX,NCPLWFMAX,N),JFLB(NRMAX,NCPLWFMAX,N))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,N),ZGLB(NRMAX,NCPLWFMAX,N))
C
      ALLOCATE (IZAZB(NRMAX))
      ALLOCATE (IZAJB(NRMAX))
      ALLOCATE (IJAZB(NRMAX))
      ALLOCATE (IJAJB(NRMAX),STAT=IA_ERR)
C
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IAUX')
      IF ( IREL.LE.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL <= 2')
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
C
C=======================================================================
C
      MZAZB(:,:,:) = C0
      MZBZA(:,:,:) = C0
      MIRR_2(:,:,:,:) = C0
      MIRR_3(:,:,:,:) = C0
      MIRR_4(:,:,:,:) = C0
C
C=======================================================================
C
C      evaluate matrix element integral functions of the type
C
C          IXAYB(r) = < X^+(E_a) | j_lam | Y(E_b) >_r
C
C=======================================================================
C
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
C
C=======================================================================
C
C      evaluate the matrix elements of the type
C
C               MZBZA = < Z^+(E_b) | sigma_lam | Z(E_a) >
C
C=======================================================================
C
C ------------------------------------------ read in wavefunctions for B
C                                                use arrays *A as buffer
C
      CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
      CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZGRA,ZFRA,JGRA,JFRA,ZGLB,ZFLB,JGLB,
     &                            JFLB,NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for A
C
      CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF,IKMCPLWF)
C
      CALL ME_SIG_ALF_B4A1(IM,ZFLB,ZGLB,ZGRA,ZFRA,MZBZA)
C
C=======================================================================
C
C      evaluate the matrix elements involving irregular solutions
C
C=======================================================================
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
               CALL ME_ALF_A2B3(LAMA2,LAMB3,JPOL,IM,IMC,PAB,ZGLA,ZFLA,
     &                          JGLA,JFLA,ZGRB,ZFRB,JGRB,JFRB,IZAZB,
     &                          IZAJB,IJAZB,IJAJB,SZAZB,SZAJB,SJAZB,
     &                          SJAJB,MZAZB,K_SELECTION_RULES)
C
C#######################################################################
C
               IF ( K_SELECTION_RULES ) THEN
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
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
                     CALL ME_SIG_ALF_BXAXW(IM,LAMB4,ZFLB,ZGLB,ZGLB,ZFLB,
     &                  LAMA1,ZGRA,ZFRA,JFRA,JGRA,IZAZB,SZAZB,IJAZB,
     &                  SJAZB,MIRR_2,LAMB4,LAMB3,JPOL)
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
                     CALL ME_SIG_ALF_BXAXW(IM,LAMB3,ZFLB,ZGLB,JGLB,JFLB,
     &                  LAMA1,ZGRA,ZFRA,ZFRA,ZGRA,IZAZB,SZAZB,IZAJB,
     &                  SZAJB,MIRR_3,LAMA2,LAMA1,JPOL)
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
                  CALL ME_SIG_ALF_BXAXW(IM,LAMB3,ZFLB,ZGLB,JGLB,JFLB,
     &                                  LAMA1,ZGRA,ZFRA,JFRA,JGRA,IZAZB,
     &                                  SZAZB,IJAJB,SJAJB,MIRR_4,LAMA1,
     &                                  LAMB3,JPOL)
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
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
C*==me_sig_alf_bxaxw.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SIG_ALF_BXAXW(IM,LAMB,ZFLB,ZGLB,JGLB,JFLB,LAMA,ZGRA,
     &                            ZFRA,JFRA,JGRA,IZAZB,SZAZB,IJAJB,
     &                            SJAJB,MIRR,I,J,JPOL)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the IRREGULAR matrix elements of the type              *
C   *                                                                  *
C   *         MIRR_2 =   < Z^+(E_b) | s_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAZB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_3 =   < J^+(E_b) | s_lam * IZAZB | Z(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_4 =   < J^+(E_b) | s_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAJB | Z(E_a) >         *
C   *                                                                  *
C   *     with    s   the spin operator sigma                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF,NCPLWFMAX
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF,IMKM_IKM,AME_G,ISMT,NPOL
      IMPLICIT NONE
C*--ME_SIG_ALF_BXAXW329
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SIG_ALF_BXAXW')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER I,IM,J,JPOL,LAMA,LAMB
      COMPLEX*16 SJAJB,SZAZB
      COMPLEX*16 IJAJB(NRMAX),IZAZB(NRMAX),JFLB(NRMAX,NCPLWFMAX,NKM),
     &           JFRA(NRMAX,NCPLWFMAX,NKM),JGLB(NRMAX,NCPLWFMAX,NKM),
     &           JGRA(NRMAX,NCPLWFMAX,NKM),MIRR(NKMMAX,NKMMAX,3,3),
     &           ZFLB(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGLB(NRMAX,NCPLWFMAX,NKM),ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      REAL*8 AFF,AGG,RARG
      COMPLEX*16 FS(:),FV(:),IAUX(:),SCNTR,VCNTR,VIRR,WIRR_F,WIRR_G,
     &           WREG_F,WREG_G
      INTEGER IA,IA_ERR,IB,IFLAG2,IKMA,IKMB,IMKMA,IMKMB,IPOL,IR,IRTOP
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FV,FS,IAUX
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FV(NRMAX),FS(NRMAX),IAUX(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IAUX')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      DO IPOL = 1,NPOL
C
         DO IR = 1,IRTOP
            FV(IR) = C0
            FS(IR) = C0
         END DO
C
         IFLAG2 = 0
C
C-----------------------------------------------------------------------
C                                                                IA IB
         DO IA = 1,NCPLWF(LAMA)
            IKMA = IKMCPLWF(IA,LAMA)
            IMKMA = IMKM_IKM(IKMA)
C
            DO IB = 1,NCPLWF(LAMB)
               IKMB = IKMCPLWF(IB,LAMB)
               IMKMB = IMKM_IKM(IKMB)
C
               AGG = AME_G(IKMB,IKMA,IPOL,ISMT)
C
               IF ( RNON0(AGG) ) THEN
C
                  IFLAG2 = 1
C
                  IF ( (IMKMA.LE.NKM) .AND. (IMKMB.LE.NKM) ) THEN
C
                     AFF = AME_G(IMKMB,IMKMA,IPOL,ISMT)
C
                     DO IR = 1,IRTOP
C
                        WREG_G = JGLB(IR,IB,LAMB)*AGG*JGRA(IR,IA,LAMA)
                        WIRR_G = ZGLB(IR,IB,LAMB)*AGG*ZGRA(IR,IA,LAMA)
C
                        WREG_F = JFLB(IR,IB,LAMB)*AFF*JFRA(IR,IA,LAMA)
                        WIRR_F = ZFLB(IR,IB,LAMB)*AFF*ZFRA(IR,IA,LAMA)
C
                        FV(IR) = FV(IR) + (WREG_G-WREG_F)*IZAZB(IR)
     &                           + (WIRR_G-WIRR_F)*IJAJB(IR)
                        FS(IR) = FS(IR) + (WIRR_G-WIRR_F)
C
                     END DO
C
                  ELSE
C
                     DO IR = 1,IRTOP
C
                        WREG_G = JGLB(IR,IB,LAMB)*AGG*JGRA(IR,IA,LAMA)
                        WIRR_G = ZGLB(IR,IB,LAMB)*AGG*ZGRA(IR,IA,LAMA)
C
                        FV(IR) = FV(IR) + WREG_G*IZAZB(IR)
     &                           + WIRR_G*IJAJB(IR)
                        FS(IR) = FS(IR) + WIRR_G
C
                     END DO
C
                  END IF
C
               END IF
C
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
            CALL CRADINT_R(IM,FV,IAUX)
            VCNTR = IAUX(IRTOP)
C
C---------------------------------------- alpha - related   surface term
C
            CALL CRADINT_R(IM,FS,IAUX)
            VIRR = IAUX(IRTOP)
C
            SCNTR = SJAJB*VIRR
C
            MIRR(I,J,IPOL,JPOL) = MIRR(I,J,IPOL,JPOL) + VCNTR + SCNTR
C
         END IF
C***********************************************************************
C
      END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
      END
C*==me_sig_alf_b4a1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_SIG_ALF_B4A1(IM,ZFLB,ZGLB,ZGRA,ZFRA,MZBZA)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the REGULAR matrix elements of the type                *
C   *                                                                  *
C   *           MZBZA = < Z^+(E_b) | sigma_lam | Z(E_a) >              *
C   *                                                                  *
C   *     with    sigma_lam   the spin operator sigma                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF,NCPLWFMAX
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NPOL,ISMT,NKM,NKMMAX,NCPLWF,AME_G,IMKM_IKM
      IMPLICIT NONE
C*--ME_SIG_ALF_B4A1494
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_SIG_ALF_B4A1')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER IM
      COMPLEX*16 MZBZA(NKMMAX,NKMMAX,3),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      REAL*8 AFF,AGG,RARG
      COMPLEX*16 FZBZA(:),IAUX(:),V2BAR
      INTEGER IA1,IA_ERR,IB4,IFLAG1,IKMA1,IKMB4,IMKMA1,IMKMB4,IPOL,IR,
     &        IRTOP,LAMA1,LAMB4
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZBZA,IAUX
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FZBZA(NRMAX),IAUX(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IAUX')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
      DO LAMB4 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
         DO LAMA1 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
            DO IPOL = 1,NPOL
C
               DO IR = 1,IRTOP
                  FZBZA(IR) = C0
               END DO
C
               IFLAG1 = 0
C
C-----------------------------------------------------------------------
C                                                                IA1 IB4
               DO IA1 = 1,NCPLWF(LAMA1)
                  IKMA1 = IKMCPLWF(IA1,LAMA1)
                  IMKMA1 = IMKM_IKM(IKMA1)
C
                  DO IB4 = 1,NCPLWF(LAMB4)
                     IKMB4 = IKMCPLWF(IB4,LAMB4)
                     IMKMB4 = IMKM_IKM(IKMB4)
C
                     AGG = AME_G(IKMB4,IKMA1,IPOL,ISMT)
C
                     IF ( RNON0(AGG) ) THEN
C
                        IFLAG1 = 1
C
                        IF ( (IMKMA1.LE.NKM) .AND. (IMKMB4.LE.NKM) )
     &                       THEN
C
                           AFF = AME_G(IMKMB4,IMKMA1,IPOL,ISMT)
C
                           DO IR = 1,IRTOP
C
                              FZBZA(IR) = FZBZA(IR)
     &                           + AGG*ZGLB(IR,IB4,LAMB4)
     &                           *ZGRA(IR,IA1,LAMA1)
     &                           - AFF*ZFLB(IR,IB4,LAMB4)
     &                           *ZFRA(IR,IA1,LAMA1)
C
                           END DO
C
                        ELSE
C
                           DO IR = 1,IRTOP
C
                              FZBZA(IR) = FZBZA(IR)
     &                           + AGG*ZGLB(IR,IB4,LAMB4)
     &                           *ZGRA(IR,IA1,LAMA1)
C
                           END DO
C
C
                        END IF
C
                     END IF
                  END DO
               END DO
C                                                                IA1 IB4
C-----------------------------------------------------------------------
C
C#######################################################################
C
               IF ( IFLAG1.EQ.1 ) THEN
C
                  CALL CRADINT_R(IM,FZBZA,IAUX)
                  V2BAR = IAUX(IRTOP)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
                  MZBZA(LAMB4,LAMA1,IPOL) = MZBZA(LAMB4,LAMA1,IPOL)
     &               + V2BAR
C
               END IF
C################################################################ IFLAG1
            END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
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
