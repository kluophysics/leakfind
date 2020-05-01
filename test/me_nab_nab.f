C*==me_nab_nab.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_NAB_NAB(LAMBOTA,LAMTOPA,IFIL_RHSA,ERYDA,LAMBOTB,
     &                      LAMTOPB,IFIL_RHSB,ERYDB,ME_CC_BRA_RWF,IT,
     &                      MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4,C,
     &                      K_CALC_ME)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate matrix elements          *
C   *                                                                  *
C   *         MZAZB  =   < Z^+(E_a) | j_lam | Z(E_b) >                 *
C   *                                                                  *
C   *         MZBZA  =   < Z^+(E_b) | j_lam | Z(E_a) >                 *
C   *                                                                  *
C   *         MIRR_2 =   < Z^+(E_b) | j_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAZB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_3 =   < J^+(E_b) | j_lam * IZAZB | Z(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *         MIRR_4 =   < J^+(E_b) | j_lam * IZAZB | J(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IJAJB | Z(E_a) >         *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator                 *
C   *                     in the NABLA form                            *
C   *                                                                  *
C   *     and the matrix element integral functions                    *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | j_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(z),(+)        *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *    K_CALC_ME = 1:    calculate MZBZA                             *
C   *              = 2:  + calculate MZAZB                             *
C   *              = 3:  + calculate MIRR_2,MIRR_3,MIRR_4              *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  ME_CC_BRA_RWF = .T.   take complex conjugate                    *
C   *                        <BRA| radial wave function                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,CI,Y00
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF_LA,IKMCPLWF_LB,
     &    IKMCPLWF_RA,IKMCPLWF_RB,BT,VT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,AN1_NAB,AN2_NAB,AB1_NAB,AB2_NAB,
     &    NPOL,NCPLWF_LA,NCPLWF_LB,NCPLWF_RA,NCPLWF_RB
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRMT,JRCRI,FULLPOT,R2DRDI,FLMSF
      IMPLICIT NONE
C*--ME_NAB_NAB54
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_NAB_NAB')
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFIL_RHSA,IFIL_RHSB,IT,K_CALC_ME,LAMBOTA,LAMBOTB,LAMTOPA,
     &        LAMTOPB
      LOGICAL ME_CC_BRA_RWF
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,3,3),MIRR_3(NKMMAX,NKMMAX,3,3),
     &           MIRR_4(NKMMAX,NKMMAX,3,3),MZAZB(NKMMAX,NKMMAX,3),
     &           MZBZA(NKMMAX,NKMMAX,3)
C
C Local variables
C
      INTEGER IA,IA_ERR,IB,IFIL_LHSA,IFIL_LHSB,IKM,IM,IR,IRSF,IRTOP,
     &        IRTOP_SPHERE,ISF,JA,JB,JPOL,LAMA1,LAMA2,LAMB3,LAMB4
      COMPLEX*16 IJAJB(:),IJAZB(:),IZAJB(:),IZAZB(:),JFLA(:,:,:),
     &           JFLB(:,:,:),JFPLA(:,:,:),JFPLB(:,:,:),JFPRA(:,:,:),
     &           JFPRB(:,:,:),JFRA(:,:,:),JFRB(:,:,:),JGLA(:,:,:),
     &           JGLB(:,:,:),JGPLA(:,:,:),JGPLB(:,:,:),JGPRA(:,:,:),
     &           JGPRB(:,:,:),JGRA(:,:,:),JGRB(:,:,:),PF,PFC,PFDUM,
     &           RMEF1(:,:,:),RMEF2(:,:,:),RMEF3(:,:,:),RMEF4(:,:,:),
     &           RMEG1(:,:,:),RMEG2(:,:,:),RMEG3(:,:,:),RMEG4(:,:,:),
     &           WB(:),WBX(:),WP(:),WPX(:),WV(:),WVX(:),ZFLA(:,:,:),
     &           ZFLB(:,:,:),ZFPLA(:,:,:),ZFPLB(:,:,:),ZFPRA(:,:,:),
     &           ZFPRB(:,:,:),ZFRA(:,:,:),ZFRB(:,:,:),ZGLA(:,:,:),
     &           ZGLB(:,:,:),ZGPLA(:,:,:),ZGPLB(:,:,:),ZGPRA(:,:,:),
     &           ZGPRB(:,:,:),ZGRA(:,:,:),ZGRB(:,:,:)
      LOGICAL K_SELECTION_RULES
      REAL*8 WSF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JFLA,JGLA,ZGLA,ZFLA,JGRB,JFRB,ZFRB,ZGRB
      ALLOCATABLE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
      ALLOCATABLE JFPLA,JGPLA,ZGPLA,ZFPLA,JGPRB,JFPRB,ZFPRB,ZGPRB
      ALLOCATABLE JFPRA,JGPRA,ZGPRA,ZFPRA,JGPLB,JFPLB,ZFPLB,ZGPLB
C
      ALLOCATABLE IZAZB,IZAJB,IJAZB,IJAJB
      ALLOCATABLE RMEF1,RMEF2
      ALLOCATABLE RMEF3,RMEF4,RMEG1,RMEG2,RMEG3,RMEG4
      ALLOCATABLE WP,WV,WB,WPX,WVX,WBX
C
      ALLOCATE (WP(NRMAX),WV(NRMAX),WB(NRMAX))
      ALLOCATE (WPX(NRMAX),WVX(NRMAX),WBX(NRMAX))
C
      ALLOCATE (RMEG1(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEG2(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEG3(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEG4(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEF1(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEF2(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEF3(NRMAX,NCPLWFMAX,NCPLWFMAX))
      ALLOCATE (RMEF4(NRMAX,NCPLWFMAX,NCPLWFMAX))
C
      ALLOCATE (JFLA(NRMAX,NCPLWFMAX,NKM),JGLA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,NKM),ZFLA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGRB(NRMAX,NCPLWFMAX,NKM),JFRB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFRB(NRMAX,NCPLWFMAX,NKM),ZGRB(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKM),JGRA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGLB(NRMAX,NCPLWFMAX,NKM),JFLB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFLB(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JFPLA(NRMAX,NCPLWFMAX,NKM),JGPLA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGPLA(NRMAX,NCPLWFMAX,NKM),ZFPLA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGPRB(NRMAX,NCPLWFMAX,NKM),JFPRB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFPRB(NRMAX,NCPLWFMAX,NKM),ZGPRB(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JFPRA(NRMAX,NCPLWFMAX,NKM),JGPRA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZGPRA(NRMAX,NCPLWFMAX,NKM),ZFPRA(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (JGPLB(NRMAX,NCPLWFMAX,NKM),JFPLB(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZFPLB(NRMAX,NCPLWFMAX,NKM),ZGPLB(NRMAX,NCPLWFMAX,NKM),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZFIP')
C
      ALLOCATE (IZAZB(NRMAX),IZAJB(NRMAX))
      ALLOCATE (IJAZB(NRMAX),IJAJB(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IZAZB')
C
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      CALL SET_IFIL_LHS(IFIL_RHSA,IFIL_LHSA)
C
      CALL SET_IFIL_LHS(IFIL_RHSB,IFIL_LHSB)
C
C=======================================================================
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      PF = 1.0D0/(1.0D0+(ERYDA+ERYDB)/(2*0.5D0*C**2))
      IF ( ME_CC_BRA_RWF ) THEN
         PFC = 1.0D0/(1.0D0+(DCONJG(ERYDA)+ERYDB)/(2*0.5D0*C**2))
C         PFSURFP = PF*(CI*C/(DCONJG(ERYDF)-ERYDI))*RWS(IM)**2
      ELSE
         PFC = PF
C         PFSURFP = PF*(CI*C/(ERYDF-ERYDI))*RWS(IM)**2
      END IF
C
      DO IR = 1,IRTOP
         WP(IR) = PF*(1.0D0/CI)*R2DRDI(IR,IM)*0.5D0
         WV(IR) = PF*(1.0D0/C)*CI*VT(IR,IT)*R2DRDI(IR,IM)
         WB(IR) = PF*(-CI/C)*CI*CI*BT(IR,IT)*R2DRDI(IR,IM)
C
         PFDUM = PF
         IF ( ME_CC_BRA_RWF ) PFDUM = PFC
C
         WPX(IR) = PFDUM*(1.0D0/CI)*R2DRDI(IR,IM)*0.5D0
         WVX(IR) = PFDUM*(1.0D0/C)*CI*VT(IR,IT)*R2DRDI(IR,IM)
         WBX(IR) = PFDUM*(-CI/C)*CI*CI*BT(IR,IT)*R2DRDI(IR,IM)
      END DO
C
C=======================================================================
C ------------------------------------------------ select specific terms
C
      IF ( 1.EQ.0 ) THEN
         WRITE (6,*) '<MENABIRR> Warning checking  ada vs v/c a'
         WRITE (6,*) '<MENABIRR> NB: switch off surface terms'
         WRITE (6,*) '<MENABIRR>     in MEADAIRR'
         CALL RINIT(NKMMAX*NKMMAX*NPOL,AN1_NAB)
         CALL RINIT(NKMMAX*NKMMAX*NPOL,AN2_NAB)
         CALL RINIT(NKMMAX*NKMMAX*NPOL,AB1_NAB)
         CALL RINIT(NKMMAX*NKMMAX*NPOL,AB2_NAB)
      END IF
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
C          IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r
C
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for LA
C
      CALL WAVFUN_READ_REL(IFIL_LHSA,IT,1,ZGLA,ZFLA,JGLA,JFLA,IRTOP,
     &                     NCPLWF_LA,IKMCPLWF_LA)
C
C weight the wave function with the radial factor R2DRDI for integration
C
Cccc      CALL WAVFUN_RAD_WGT_REL(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF_LA,'R')
C
C--------- initialisation needed for ME check routines incl. core states
C
      ZGRB = C0
      ZFRB = C0
      JGRB = C0
      JFRB = C0
C
C ----------------------------------------- read in wavefunctions for RB
C
      CALL WAVFUN_READ_REL(IFIL_RHSB,IT,1,ZGRB,ZFRB,JGRB,JFRB,IRTOP,
     &                     NCPLWF_RB,IKMCPLWF_RB)
C
C------------- take complex conjugate of wave function ZGLA if requested
C
      IF ( ME_CC_BRA_RWF ) THEN
         DO IKM = 1,NKM
            ZGLA(:,:,IKM) = DCONJG(ZGLA(:,:,IKM))
            ZFLA(:,:,IKM) = DCONJG(ZFLA(:,:,IKM))
            JGLA(:,:,IKM) = DCONJG(JGLA(:,:,IKM))
            JFLA(:,:,IKM) = DCONJG(JFLA(:,:,IKM))
         END DO
      END IF
C
C-----------------------------------------------------------------------
C                     differentiate wave functions
C-----------------------------------------------------------------------
C
      DO JA = LAMBOTA,LAMTOPA
         DO IA = 1,NCPLWF_LA(JA)
            CALL CDIFFER(IM,ZGLA(1,IA,JA),ZGPLA(1,IA,JA))
            CALL CDIFFER(IM,ZFLA(1,IA,JA),ZFPLA(1,IA,JA))
            CALL CDIFFER(IM,JGLA(1,IA,JA),JGPLA(1,IA,JA))
            CALL CDIFFER(IM,JFLA(1,IA,JA),JFPLA(1,IA,JA))
         END DO
      END DO
C
      DO JB = LAMBOTB,LAMTOPB
         DO IB = 1,NCPLWF_RB(JB)
            CALL CDIFFER(IM,ZGRB(1,IB,JB),ZGPRB(1,IB,JB))
            CALL CDIFFER(IM,ZFRB(1,IB,JB),ZFPRB(1,IB,JB))
            CALL CDIFFER(IM,JGRB(1,IB,JB),JGPRB(1,IB,JB))
            CALL CDIFFER(IM,JFRB(1,IB,JB),JFPRB(1,IB,JB))
         END DO
      END DO
C
C=======================================================================
C
C      evaluate the matrix elements of the type
C
C               MZBZA = < Z^+(E_b) | H_lam | Z(E_a) >
C
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for LB
C
      CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGLB,ZFLB,JGLB,JFLB,IRTOP,
     &                     NCPLWF_LB,IKMCPLWF_LB)
C
C ----------------------------------------- read in wavefunctions for RA
C
      CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF_RA,IKMCPLWF_RA)
C
C-------------- take complex conjugate of wave function ZGRA if requested
C
      IF ( ME_CC_BRA_RWF ) THEN
         DO IKM = 1,NKM
            ZGLB(:,:,IKM) = DCONJG(ZGLB(:,:,IKM))
            ZFLB(:,:,IKM) = DCONJG(ZFLB(:,:,IKM))
            JGLB(:,:,IKM) = DCONJG(JGLB(:,:,IKM))
            JFLB(:,:,IKM) = DCONJG(JFLB(:,:,IKM))
         END DO
      END IF
C
C-----------------------------------------------------------------------
C                     differentiate wave functions
C-----------------------------------------------------------------------
C
      DO JA = LAMBOTA,LAMTOPA
         DO IA = 1,NCPLWF_RA(JA)
            CALL CDIFFER(IM,ZGRA(1,IA,JA),ZGPRA(1,IA,JA))
            CALL CDIFFER(IM,ZFRA(1,IA,JA),ZFPRA(1,IA,JA))
            CALL CDIFFER(IM,JGRA(1,IA,JA),JGPRA(1,IA,JA))
            CALL CDIFFER(IM,JFRA(1,IA,JA),JFPRA(1,IA,JA))
         END DO
      END DO
C
      DO JB = LAMBOTB,LAMTOPB
         DO IB = 1,NCPLWF_LB(JB)
            CALL CDIFFER(IM,ZGLB(1,IB,JB),ZGPLB(1,IB,JB))
            CALL CDIFFER(IM,ZFLB(1,IB,JB),ZFPLB(1,IB,JB))
            CALL CDIFFER(IM,JGLB(1,IB,JB),JGPLB(1,IB,JB))
            CALL CDIFFER(IM,JFLB(1,IB,JB),JFPLB(1,IB,JB))
         END DO
      END DO
C
C=======================================================================
C
C    weight all LHS shape functions with shape function f_00 * Y00
C
C=======================================================================
C
      IF ( FULLPOT ) THEN
C
         IRTOP_SPHERE = JRMT(IM)
         ISF = 1
C
         IF ( ABS(1D0-FLMSF(1,ISF,IM)*Y00).GT.1D-12 )
     &        CALL STOP_MESSAGE(ROUTINE,'f00*Y00 <> 1 for JRMT')
C
         DO IR = IRTOP_SPHERE + 1,JRCRI(IM)
C
            IRSF = IR - IRTOP_SPHERE
            WSF = FLMSF(IRSF,ISF,IM)*Y00
C
            ZGLA(IR,:,:) = WSF*ZGLA(IR,:,:)
            ZFLA(IR,:,:) = WSF*ZFLA(IR,:,:)
            JGLA(IR,:,:) = WSF*JGLA(IR,:,:)
            JFLA(IR,:,:) = WSF*JFLA(IR,:,:)
C
            ZGPLA(IR,:,:) = WSF*ZGPLA(IR,:,:)
            ZFPLA(IR,:,:) = WSF*ZFPLA(IR,:,:)
            JGPLA(IR,:,:) = WSF*JGPLA(IR,:,:)
            JFPLA(IR,:,:) = WSF*JFPLA(IR,:,:)
C
            ZGLB(IR,:,:) = WSF*ZGLB(IR,:,:)
            ZFLB(IR,:,:) = WSF*ZFLB(IR,:,:)
            JGLB(IR,:,:) = WSF*JGLB(IR,:,:)
            JFLB(IR,:,:) = WSF*JFLB(IR,:,:)
C
            ZGPLB(IR,:,:) = WSF*ZGPLB(IR,:,:)
            ZFPLB(IR,:,:) = WSF*ZFPLB(IR,:,:)
            JGPLB(IR,:,:) = WSF*JGPLB(IR,:,:)
            JFPLB(IR,:,:) = WSF*JFPLB(IR,:,:)
C
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C
      CALL MENAB_B4A1(IM,LAMBOTB,LAMTOPB,ZGLB,ZFLB,ZGPLB,ZFPLB,LAMBOTA,
     &                LAMTOPA,ZGRA,ZFRA,ZGPRA,ZFPRA,RMEG1,RMEG2,RMEG3,
     &                RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WPX,WVX,WBX,MZBZA)
C
      IF ( K_CALC_ME.EQ.1 ) RETURN
C
C=======================================================================
C
C      evaluate the matrix elements involving irregular solutions
C
C=======================================================================
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
      LOOP_LAMB3:DO LAMB3 = LAMBOTB,LAMTOPB
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
         LOOP_LAMA2:DO LAMA2 = LAMBOTA,LAMTOPA
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
            LOOP_JPOL:DO JPOL = 1,NPOL
C
               CALL MENAB_A2B3(LAMA2,LAMB3,JPOL,IM,ZGLA,ZFLA,ZGPLA,
     &                         ZFPLA,JGLA,JFLA,JGPLA,JFPLA,ZGRB,ZFRB,
     &                         ZGPRB,ZFPRB,JGRB,JFRB,JGPRB,JFPRB,IZAZB,
     &                         IZAJB,IJAZB,IJAJB,RMEG1,RMEG2,RMEG3,
     &                         RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WP,WV,WB,
     &                         MZAZB,K_SELECTION_RULES)
C
C#######################################################################
C
               IF ( .NOT.K_SELECTION_RULES ) CYCLE LOOP_JPOL
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
               DO LAMB4 = LAMBOTB,LAMTOPB
C
                  CALL MENAB_BXAXW(IM,LAMB4,ZGLB,ZFLB,ZGPLB,ZFPLB,ZGLB,
     &                             ZFLB,ZGPLB,ZFPLB,LAMA1,ZGRA,ZFRA,
     &                             ZGPRA,ZFPRA,JGRA,JFRA,JGPRA,JFPRA,
     &                             RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,RMEF2,
     &                             RMEF3,RMEF4,WPX,WVX,WBX,IZAZB,IJAZB,
     &                             MIRR_2,LAMB4,LAMB3,JPOL)
C
               END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB LAMB4
C
               IF ( K_CALC_ME.EQ.2 ) CYCLE LOOP_JPOL
C
C***********************************************************************
C***********************************************************************
C     TERM 3    - TAU(1,2,a) * MIRR_3(1,2)                 LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
               DO LAMA1 = LAMBOTA,LAMTOPA
C
                  CALL MENAB_BXAXW(IM,LAMB3,ZGLB,ZFLB,ZGPLB,ZFPLB,JGLB,
     &                             JFLB,JGPLB,JFPLB,LAMA1,ZGRA,ZFRA,
     &                             ZGPRA,ZFPRA,ZGRA,ZFRA,ZGPRA,ZFPRA,
     &                             RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,RMEF2,
     &                             RMEF3,RMEF4,WPX,WVX,WBX,IZAZB,IZAJB,
     &                             MIRR_3,LAMA2,LAMA1,JPOL)
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
               CALL MENAB_BXAXW(IM,LAMB3,ZGLB,ZFLB,ZGPLB,ZFPLB,JGLB,
     &                          JFLB,JGPLB,JFPLB,LAMA1,ZGRA,ZFRA,ZGPRA,
     &                          ZFPRA,JGRA,JFRA,JGPRA,JFPRA,RMEG1,RMEG2,
     &                          RMEG3,RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WPX,
     &                          WVX,WBX,IZAZB,IJAJB,MIRR_4,LAMA1,LAMB3,
     &                          JPOL)
C
C***********************************************************************
C
            END DO LOOP_JPOL
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
         END DO LOOP_LAMA2
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO LOOP_LAMB3
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( IT.LE.0 ) THEN
         IM = 100
         WRITE (IM,*) 'MZAZB'
         WRITE (100,'(2e22.14)') MZAZB
         WRITE (IM,*) 'MZBZA'
         WRITE (100,'(2e22.14)') MZBZA
         WRITE (IM,*) 'MIRR_2'
         WRITE (100,'(2e22.14)') MIRR_2
         WRITE (IM,*) 'MIRR_3'
         WRITE (100,'(2e22.14)') MIRR_3
         WRITE (IM,*) 'MIRR_4'
         WRITE (100,'(2e22.14)') MIRR_4
         STOP
      END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      END
C*==menab_b4a1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MENAB_B4A1(IM,LAMBOTB,LAMTOPB,ZGLB,ZFLB,ZGPLB,ZFPLB,
     &                      LAMBOTA,LAMTOPA,ZGRA,ZFRA,ZGPRA,ZFPRA,RMEG1,
     &                      RMEG2,RMEG3,RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,
     &                      WPX,WVX,WBX,MZBZA)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements of the type                        *
C   *                                                                  *
C   *           MZBZA = < Z^+(E_b) | H_lam | Z(E_a) >                  *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator                 *
C   *                                                                  *
C   *             in the NABLA form                                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX,IKMCPLWF_LB,IKMCPLWF_RA
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRCRI,JRWS
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,
     &    AB1_NAB,AB2_NAB,NPOL,NCPLWF_LB,NCPLWF_RA
      IMPLICIT NONE
C*--MENAB_B4A1528
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MENAB_B4A1')
C
C Dummy arguments
C
      INTEGER IM,LAMBOTA,LAMBOTB,LAMTOPA,LAMTOPB
      COMPLEX*16 MZBZA(NKMMAX,NKMMAX,3),RMEF1(NRMAX,NCPLWFMAX,NCPLWFMAX)
     &           ,RMEF2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF4(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG4(NRMAX,NCPLWFMAX,NCPLWFMAX),WBX(NRMAX),WPX(NRMAX),
     &           WVX(NRMAX),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFPLB(NRMAX,NCPLWFMAX,NKM),ZFPRA(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGPLB(NRMAX,NCPLWFMAX,NKM),ZGPRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 FZBZA(:),V2BAR
      INTEGER IA1,IA_ERR,IB4,IKMA1,IKMB4,IPOL,IRTOP,LAMA1,LAMB4
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZBZA
C
      ALLOCATE (FZBZA(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FZBZA')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB4
      DO LAMB4 = LAMBOTB,LAMTOPB
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA1
         DO LAMA1 = LAMBOTA,LAMTOPA
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
            LOOP_IPOL:DO IPOL = 1,NPOL
C
               FZBZA(:) = C0
C
               DO IB4 = 1,NCPLWF_LB(LAMB4)
                  IKMB4 = IKMCPLWF_LB(IB4,LAMB4)
                  DO IA1 = 1,NCPLWF_RA(LAMA1)
                     IKMA1 = IKMCPLWF_RA(IA1,LAMA1)
                     IF ( RNON0(AN1_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                     IF ( RNON0(AN2_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                     IF ( RNON0(AV1_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                     IF ( RNON0(AV2_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                     IF ( RNON0(AB1_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                     IF ( RNON0(AB2_NAB(IKMB4,IKMA1,IPOL)) ) GOTO 10
                  END DO
               END DO
C     ---------------------------------- all angular matrix elements = 0
               CYCLE LOOP_IPOL
C
C     ------------------------------ non-0 angular matrix elements found
C
C
 10            CONTINUE
               CALL MENAB_MEAUX(IPOL,IM,IRTOP,ZGLB(1,1,LAMB4),
     &                          ZFLB(1,1,LAMB4),ZGPLB(1,1,LAMB4),
     &                          ZFPLB(1,1,LAMB4),IKMCPLWF_LB(1,LAMB4),
     &                          NCPLWF_LB(LAMB4),ZGRA(1,1,LAMA1),
     &                          ZFRA(1,1,LAMA1),ZGPRA(1,1,LAMA1),
     &                          ZFPRA(1,1,LAMA1),IKMCPLWF_RA(1,LAMA1),
     &                          NCPLWF_RA(LAMA1),RMEG1,RMEG2,RMEG3,
     &                          RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WPX,WVX,
     &                          WBX,FZBZA)
C
               CALL CRADINT(IM,FZBZA,V2BAR)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
               MZBZA(LAMB4,LAMA1,IPOL) = MZBZA(LAMB4,LAMA1,IPOL) + V2BAR
C
            END DO LOOP_IPOL
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
C*==menab_a2b3.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MENAB_A2B3(LAMA2,LAMB3,JPOL,IM,ZGLA,ZFLA,ZGPLA,ZFPLA,
     &                      JGLA,JFLA,JGPLA,JFPLA,ZGRB,ZFRB,ZGPRB,ZFPRB,
     &                      JGRB,JFRB,JGPRB,JFPRB,IZAZB,IZAJB,IJAZB,
     &                      IJAJB,RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,RMEF2,
     &                      RMEF3,RMEF4,WP,WV,WB,MZAZB,
     &                      K_SELECTION_RULES)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements                                    *
C   *                                                                  *
C   *              M = < Z^+(E_a) | H_lam | Z(E_b) >                   *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator                 *
C   *                     in the NABLA form                            *
C   *                                                                  *
C   *  driver routine for case A2 - B3                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX,IKMCPLWF_LA,IKMCPLWF_RB
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRCRI,JRWS
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,
     &    AB1_NAB,AB2_NAB,NCPLWF_LA,NCPLWF_RB
      IMPLICIT NONE
C*--MENAB_A2B3676
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MENAB_A2B3')
C
C Dummy arguments
C
      INTEGER IM,JPOL,LAMA2,LAMB3
      LOGICAL K_SELECTION_RULES
      COMPLEX*16 IJAJB(NRMAX),IJAZB(NRMAX),IZAJB(NRMAX),IZAZB(NRMAX),
     &           JFLA(NRMAX,NCPLWFMAX,NKM),JFPLA(NRMAX,NCPLWFMAX,NKM),
     &           JFPRB(NRMAX,NCPLWFMAX,NKM),JFRB(NRMAX,NCPLWFMAX,NKM),
     &           JGLA(NRMAX,NCPLWFMAX,NKM),JGPLA(NRMAX,NCPLWFMAX,NKM),
     &           JGPRB(NRMAX,NCPLWFMAX,NKM),JGRB(NRMAX,NCPLWFMAX,NKM),
     &           MZAZB(NKMMAX,NKMMAX,3),RMEF1(NRMAX,NCPLWFMAX,NCPLWFMAX)
     &           ,RMEF2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF4(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG4(NRMAX,NCPLWFMAX,NCPLWFMAX),WB(NRMAX),WP(NRMAX),
     &           WV(NRMAX),ZFLA(NRMAX,NCPLWFMAX,NKM),
     &           ZFPLA(NRMAX,NCPLWFMAX,NKM),ZFPRB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRB(NRMAX,NCPLWFMAX,NKM),ZGLA(NRMAX,NCPLWFMAX,NKM),
     &           ZGPLA(NRMAX,NCPLWFMAX,NKM),ZGPRB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRB(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 FJAJB(:),FJAZB(:),FZAJB(:),FZAZB(:)
      INTEGER IA2,IA_ERR,IB3,IRTOP,JF,JI
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZAZB,FZAJB,FJAZB,FJAJB
C
      ALLOCATE (FZAZB(NRMAX),FZAJB(NRMAX),FJAZB(NRMAX),FJAJB(NRMAX),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FZAZB')
C
      K_SELECTION_RULES = .FALSE.
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      DO IB3 = 1,NCPLWF_RB(LAMB3)
         JI = IKMCPLWF_RB(IB3,LAMB3)
         DO IA2 = 1,NCPLWF_LA(LAMA2)
            JF = IKMCPLWF_LA(IA2,LAMA2)
            IF ( RNON0(AN1_NAB(JF,JI,JPOL)) ) GOTO 100
            IF ( RNON0(AN2_NAB(JF,JI,JPOL)) ) GOTO 100
            IF ( RNON0(AV1_NAB(JF,JI,JPOL)) ) GOTO 100
            IF ( RNON0(AV2_NAB(JF,JI,JPOL)) ) GOTO 100
            IF ( RNON0(AB1_NAB(JF,JI,JPOL)) ) GOTO 100
            IF ( RNON0(AB2_NAB(JF,JI,JPOL)) ) GOTO 100
         END DO
      END DO
C
C     ---------------------------------- all angular matrix elements = 0
      RETURN
C
C     ------------------------------ non-0 angular matrix elements found
C
 100  CONTINUE
      K_SELECTION_RULES = .TRUE.
C
      CALL MENAB_MEAUX(JPOL,IM,IRTOP,ZGLA(1,1,LAMA2),ZFLA(1,1,LAMA2),
     &                 ZGPLA(1,1,LAMA2),ZFPLA(1,1,LAMA2),
     &                 IKMCPLWF_LA(1,LAMA2),NCPLWF_LA(LAMA2),
     &                 ZGRB(1,1,LAMB3),ZFRB(1,1,LAMB3),ZGPRB(1,1,LAMB3),
     &                 ZFPRB(1,1,LAMB3),IKMCPLWF_RB(1,LAMB3),
     &                 NCPLWF_RB(LAMB3),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                 RMEF2,RMEF3,RMEF4,WP,WV,WB,FZAZB)
C
      CALL MENAB_MEAUX(JPOL,IM,IRTOP,ZGLA(1,1,LAMA2),ZFLA(1,1,LAMA2),
     &                 ZGPLA(1,1,LAMA2),ZFPLA(1,1,LAMA2),
     &                 IKMCPLWF_LA(1,LAMA2),NCPLWF_LA(LAMA2),
     &                 JGRB(1,1,LAMB3),JFRB(1,1,LAMB3),JGPRB(1,1,LAMB3),
     &                 JFPRB(1,1,LAMB3),IKMCPLWF_RB(1,LAMB3),
     &                 NCPLWF_RB(LAMB3),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                 RMEF2,RMEF3,RMEF4,WP,WV,WB,FZAJB)
C
      CALL MENAB_MEAUX(JPOL,IM,IRTOP,JGLA(1,1,LAMA2),JFLA(1,1,LAMA2),
     &                 JGPLA(1,1,LAMA2),JFPLA(1,1,LAMA2),
     &                 IKMCPLWF_LA(1,LAMA2),NCPLWF_LA(LAMA2),
     &                 ZGRB(1,1,LAMB3),ZFRB(1,1,LAMB3),ZGPRB(1,1,LAMB3),
     &                 ZFPRB(1,1,LAMB3),IKMCPLWF_RB(1,LAMB3),
     &                 NCPLWF_RB(LAMB3),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                 RMEF2,RMEF3,RMEF4,WP,WV,WB,FJAZB)
C
      CALL MENAB_MEAUX(JPOL,IM,IRTOP,JGLA(1,1,LAMA2),JFLA(1,1,LAMA2),
     &                 JGPLA(1,1,LAMA2),JFPLA(1,1,LAMA2),
     &                 IKMCPLWF_LA(1,LAMA2),NCPLWF_LA(LAMA2),
     &                 JGRB(1,1,LAMB3),JFRB(1,1,LAMB3),JGPRB(1,1,LAMB3),
     &                 JFPRB(1,1,LAMB3),IKMCPLWF_RB(1,LAMB3),
     &                 NCPLWF_RB(LAMB3),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                 RMEF2,RMEF3,RMEF4,WP,WV,WB,FJAJB)
C
      CALL CRADINT_R(IM,FZAZB,IZAZB)
C
      CALL CRADINT_INW_R(IM,FZAJB,IZAJB)
C
      CALL CRADINT_INW_R(IM,FJAZB,IJAZB)
C
      CALL CRADINT_INW_R(IM,FJAJB,IJAJB)
C
C-----------------------------------------------------------------------
C
      MZAZB(LAMA2,LAMB3,JPOL) = MZAZB(LAMA2,LAMB3,JPOL) + IZAZB(IRTOP)
C
C-----------------------------------------------------------------------
C
      END
C*==menab_bxaxw.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MENAB_BXAXW(IM,LAMB,ZGLB,ZFLB,ZGPLB,ZFPLB,JGLB,JFLB,
     &                       JGPLB,JFPLB,LAMA,ZGRA,ZFRA,ZGPRA,ZFPRA,
     &                       JGRA,JFRA,JGPRA,JFPRA,RMEG1,RMEG2,RMEG3,
     &                       RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WPX,WVX,WBX,
     &                       IZAZB,IJAJB,MIRR,I,J,JPOL)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate matrix element integral functions of the type          *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_TYPES,ONLY:NCPLWFMAX,IKMCPLWF_LB,IKMCPLWF_RA
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,
     &    AB1_NAB,AB2_NAB,NPOL,NCPLWF_LB,NCPLWF_RA
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRCRI,JRWS
      IMPLICIT NONE
C*--MENAB_BXAXW833
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MENAB_BXAXW')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER I,IM,J,JPOL,LAMA,LAMB
      COMPLEX*16 IJAJB(NRMAX),IZAZB(NRMAX),JFLB(NRMAX,NCPLWFMAX,NKM),
     &           JFPLB(NRMAX,NCPLWFMAX,NKM),JFPRA(NRMAX,NCPLWFMAX,NKM),
     &           JFRA(NRMAX,NCPLWFMAX,NKM),JGLB(NRMAX,NCPLWFMAX,NKM),
     &           JGPLB(NRMAX,NCPLWFMAX,NKM),JGPRA(NRMAX,NCPLWFMAX,NKM),
     &           JGRA(NRMAX,NCPLWFMAX,NKM),MIRR(NKMMAX,NKMMAX,3,3),
     &           RMEF1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF4(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG4(NRMAX,NCPLWFMAX,NCPLWFMAX),WBX(NRMAX),WPX(NRMAX),
     &           WVX(NRMAX),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFPLB(NRMAX,NCPLWFMAX,NKM),ZFPRA(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGPLB(NRMAX,NCPLWFMAX,NKM),ZGPRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 FV(:),VCNTR,WIRR(:),WREG(:)
      INTEGER IA,IA_ERR,IB,IKMA,IKMB,IPOL,IR,IRTOP
      REAL*8 RARG
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FV,WREG,WIRR
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (WREG(NRMAX),WIRR(NRMAX),FV(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WREG')
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      LOOP_IPOL:DO IPOL = 1,NPOL
C
         WREG(:) = C0
         WIRR(:) = C0
C
         DO IB = 1,NCPLWF_LB(LAMB)
            IKMB = IKMCPLWF_LB(IB,LAMB)
            DO IA = 1,NCPLWF_RA(LAMA)
               IKMA = IKMCPLWF_RA(IA,LAMA)
               IF ( RNON0(AN1_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
               IF ( RNON0(AN2_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
               IF ( RNON0(AV1_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
               IF ( RNON0(AV2_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
               IF ( RNON0(AB1_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
               IF ( RNON0(AB2_NAB(IKMB,IKMA,IPOL)) ) GOTO 50
            END DO
         END DO
C     ---------------------------------- all angular matrix elements = 0
         CYCLE LOOP_IPOL
C
C     ------------------------------ non-0 angular matrix elements found
C
C
 50      CONTINUE
         CALL MENAB_MEAUX(IPOL,IM,IRTOP,JGLB(1,1,LAMB),JFLB(1,1,LAMB),
     &                    JGPLB(1,1,LAMB),JFPLB(1,1,LAMB),
     &                    IKMCPLWF_LB(1,LAMB),NCPLWF_LB(LAMB),
     &                    JGRA(1,1,LAMA),JFRA(1,1,LAMA),JGPRA(1,1,LAMA),
     &                    JFPRA(1,1,LAMA),IKMCPLWF_RA(1,LAMA),
     &                    NCPLWF_RA(LAMA),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                    RMEF2,RMEF3,RMEF4,WPX,WVX,WBX,WREG)
C
         CALL MENAB_MEAUX(IPOL,IM,IRTOP,ZGLB(1,1,LAMB),ZFLB(1,1,LAMB),
     &                    ZGPLB(1,1,LAMB),ZFPLB(1,1,LAMB),
     &                    IKMCPLWF_LB(1,LAMB),NCPLWF_LB(LAMB),
     &                    ZGRA(1,1,LAMA),ZFRA(1,1,LAMA),ZGPRA(1,1,LAMA),
     &                    ZFPRA(1,1,LAMA),IKMCPLWF_RA(1,LAMA),
     &                    NCPLWF_RA(LAMA),RMEG1,RMEG2,RMEG3,RMEG4,RMEF1,
     &                    RMEF2,RMEF3,RMEF4,WPX,WVX,WBX,WIRR)
C
         DO IR = 1,IRTOP
            FV(IR) = WREG(IR)*IZAZB(IR) + WIRR(IR)*IJAJB(IR)
         END DO
C
         CALL CRADINT(IM,FV,VCNTR)
C
         MIRR(I,J,IPOL,JPOL) = MIRR(I,J,IPOL,JPOL) + VCNTR
C
      END DO LOOP_IPOL
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
      END
C*==menab_meaux.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MENAB_MEAUX(JPOL,IM,IRTOP,ZGLA,ZFLA,ZGPLA,ZFPLA,
     &                       IKMCPLWF_LA,NCPLWF_LA,ZGRB,ZFRB,ZGPRB,
     &                       ZFPRB,IKMCPLWF_RB,NCPLWF_RB,RMEG1,RMEG2,
     &                       RMEG3,RMEG4,RMEF1,RMEF2,RMEF3,RMEF4,WP,WV,
     &                       WB,FZAZB)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements                                    *
C   *                                                                  *
C   *              M = < Z^+(E_a) | H_lam | Z(E_b) >                   *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator                 *
C   *                     in the NABLA form                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX
      USE MOD_RMESH,ONLY:NRMAX,R
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,AB1_NAB,
     &    AB2_NAB,L_IKM,LB_IKM,IMKM_IKM
      IMPLICIT NONE
C*--MENAB_MEAUX978
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IM,IRTOP,JPOL,NCPLWF_LA,NCPLWF_RB
      COMPLEX*16 FZAZB(NRMAX),RMEF1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEF4(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG1(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG2(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG3(NRMAX,NCPLWFMAX,NCPLWFMAX),
     &           RMEG4(NRMAX,NCPLWFMAX,NCPLWFMAX),WB(NRMAX),WP(NRMAX),
     &           WV(NRMAX),ZFLA(NRMAX,NCPLWFMAX),ZFPLA(NRMAX,NCPLWFMAX),
     &           ZFPRB(NRMAX,NCPLWFMAX),ZFRB(NRMAX,NCPLWFMAX),
     &           ZGLA(NRMAX,NCPLWFMAX),ZGPLA(NRMAX,NCPLWFMAX),
     &           ZGPRB(NRMAX,NCPLWFMAX),ZGRB(NRMAX,NCPLWFMAX)
      INTEGER IKMCPLWF_LA(NCPLWFMAX),IKMCPLWF_RB(NCPLWFMAX)
C
C Local variables
C
      REAL*8 AF1,AF2,AF3,AF4
      COMPLEX*16 FFA12,FFA34,FG,FGA,GF,GFA,GGA12,GGA34
      INTEGER IA,IB,IKMA,IKMB,IMKMA,IMKMB,IR,LB_A,LB_B,L_A,L_B
C
C*** End of declarations rewritten by SPAG
C
      FZAZB(:) = C0
C
C-----------------------------------------------------------------------
C                                                                  IA IB
      LOOP_IA_A:DO IA = 1,NCPLWF_LA
         L_A = L_IKM(IKMCPLWF_LA(IA))
         LB_A = LB_IKM(IKMCPLWF_LA(IA))
C
         LOOP_IB_A:DO IB = 1,NCPLWF_RB
            L_B = L_IKM(IKMCPLWF_RB(IB))
            LB_B = LB_IKM(IKMCPLWF_RB(IB))
C
C ----------------------------------------------------------------------
C                radial matrix elements for p . a - term
C ----------------------------------------------------------------------
C
            DO IR = 1,IRTOP
C
               RMEG1(IR,IA,IB) = ZGLA(IR,IA)
     &                           *(ZGPRB(IR,IB)-ZGRB(IR,IB)*L_B/R(IR,IM)
     &                           )
               RMEG2(IR,IA,IB) = ZGLA(IR,IA)
     &                           *(ZGPRB(IR,IB)+ZGRB(IR,IB)*(L_B+1)
     &                           /R(IR,IM))
               RMEG3(IR,IA,IB) = ZGRB(IR,IB)
     &                           *(ZGPLA(IR,IA)-ZGLA(IR,IA)*L_A/R(IR,IM)
     &                           )
               RMEG4(IR,IA,IB) = ZGRB(IR,IB)
     &                           *(ZGPLA(IR,IA)+ZGLA(IR,IA)*(L_A+1)
     &                           /R(IR,IM))
C
               RMEF1(IR,IA,IB) = ZFLA(IR,IA)
     &                           *(ZFPRB(IR,IB)-ZFRB(IR,IB)*LB_B/R(IR,
     &                           IM))
               RMEF2(IR,IA,IB) = ZFLA(IR,IA)
     &                           *(ZFPRB(IR,IB)+ZFRB(IR,IB)*(LB_B+1)
     &                           /R(IR,IM))
               RMEF3(IR,IA,IB) = ZFRB(IR,IB)
     &                           *(ZFPLA(IR,IA)-ZFLA(IR,IA)*LB_A/R(IR,
     &                           IM))
               RMEF4(IR,IA,IB) = ZFRB(IR,IB)
     &                           *(ZFPLA(IR,IA)+ZFLA(IR,IA)*(LB_A+1)
     &                           /R(IR,IM))
            END DO
C
         END DO LOOP_IB_A
C
      END DO LOOP_IA_A
C
C ----------------------------------------------------------------------
C                   calculate  TOTAL  matrix elements
C ----------------------------------------------------------------------
C
      LOOP_IA_B:DO IA = 1,NCPLWF_LA
         IKMA = IKMCPLWF_LA(IA)
         IMKMA = IMKM_IKM(IKMA)
C
         LOOP_IB_B:DO IB = 1,NCPLWF_RB
            IKMB = IKMCPLWF_RB(IB)
            IMKMB = IMKM_IKM(IKMB)
C
            IF ( (IMKMA.LE.NKM) .AND. (IMKMB.LE.NKM) ) THEN
               AF1 = AN1_NAB(IMKMA,IMKMB,JPOL)
               AF2 = AN2_NAB(IMKMA,IMKMB,JPOL)
               AF3 = AN2_NAB(IMKMA,IMKMB,JPOL)
               AF4 = AN1_NAB(IMKMA,IMKMB,JPOL)
            ELSE
               AF1 = 0.0D0
               AF2 = 0.0D0
               AF3 = 0.0D0
               AF4 = 0.0D0
            END IF
C
            LOOP_IR:DO IR = 1,IRTOP
C                                                     ******************
C =================================================== *      p . a     *
C                                                     ******************
C prefactor: PF * hbar/i (momentum operator)
C
C              WP(IR) = PF*(1.0D0/CI)*R2DRDI(IR,IM)*0.5D0
C
               GGA12 = RMEG1(IR,IA,IB)*AN1_NAB(IKMA,IKMB,JPOL)
     &                 - RMEG2(IR,IA,IB)*AN2_NAB(IKMA,IKMB,JPOL)
C
               FFA12 = RMEF1(IR,IA,IB)*AF1 - RMEF2(IR,IA,IB)*AF2
C
               GGA34 = RMEG3(IR,IA,IB)*AN2_NAB(IKMA,IKMB,JPOL)
     &                 - RMEG4(IR,IA,IB)*AN1_NAB(IKMA,IKMB,JPOL)
C
               FFA34 = RMEF3(IR,IA,IB)*AF3 - RMEF4(IR,IA,IB)*AF4
C
               FZAZB(IR) = FZAZB(IR) + WP(IR)*(GGA12+FFA12+GGA34+FFA34)
C
C ======================================================================
C
               GF = ZGLA(IR,IA)*ZFRB(IR,IB)
               FG = ZFLA(IR,IA)*ZGRB(IR,IB)
C                                                     ******************
C =================================================== *  V/c alfa . a  *
C                                                     ******************
C prefactor: PF * V/c * i  (NB: factor i because of mixing g and f)
C
C         WV(IR) = PF*(1.0D0/C)*CI*VT(IR,IT)*R2DRDI(IR,IM)
C
               GFA = GF*AV1_NAB(IKMA,IKMB,JPOL)
               FGA = FG*AV2_NAB(IKMA,IKMB,JPOL)
C
               FZAZB(IR) = FZAZB(IR) + WV(IR)*(GFA-FGA)
C
C                                                ***********************
C ============================================== * -ibB/c (alfa x a)_z *
C                                                ***********************
C prefactor: PF *   (-i/c)          *         i        *       i
C                     ^                       ^                ^
C                     |                       |                |
C                  operator        alpha (mix g/f)     compensate 1/i
C                                                        in ABx_NAB
C
C                        BR3 = PF*(1d0/C)*CI*BT(IR,IT)*R2DRDI(IR,IM)
C
C          WB(IR) = PF*(-CI/C)*CI*CI*BT(IR,IT)*R2DRDI(IR,IM)
C
               GFA = GF*AB1_NAB(IKMA,IKMB,JPOL)
               FGA = FG*AB2_NAB(IKMA,IKMB,JPOL)
C
               FZAZB(IR) = FZAZB(IR) + WB(IR)*(GFA+FGA)
C
            END DO LOOP_IR
C
         END DO LOOP_IB_B
C
      END DO LOOP_IA_B
C
      END
