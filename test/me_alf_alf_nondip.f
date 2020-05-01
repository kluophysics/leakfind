C*==me_alf_alf_nondip.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_ALF_NONDIP(LAMBOTA,LAMTOPA,IFIL_RHSA,ERYDA_KM,
     &                             LAMBOTB,LAMTOPB,IFIL_RHSB,ERYDB,
     &                             ME_CC_BRA_RWF,IT,MZAZB,MZBZA,MIRR_2,
     &                             MIRR_3,MIRR_4,C,K_CALC_ME)
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
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF_LA,IKMCPLWF_LB,
     &    IKMCPLWF_RA,IKMCPLWF_RB
      USE MOD_RMESH,ONLY:JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF_LA,NCPLWF_LB,NCPLWF_RA,
     &    NCPLWF_RB,AME_G,AG_RGNT,NKM_EXT,NLM_AME_RLM_EXT,ISMT,
     &    A_SIG_RLM,NPOL
      USE MOD_CALCMODE,ONLY:PUBLIC_VERSION
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_ALF_NONDIP')
      LOGICAL CHECK_BXAXW
      PARAMETER (CHECK_BXAXW=.FALSE.)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDB
      INTEGER IFIL_RHSA,IFIL_RHSB,IT,K_CALC_ME,LAMBOTA,LAMBOTB,LAMTOPA,
     &        LAMTOPB
      LOGICAL ME_CC_BRA_RWF
      COMPLEX*16 ERYDA_KM(NKMMAX),MIRR_2(NKMMAX,NKMMAX,3,3,3),
     &           MIRR_3(NKMMAX,NKMMAX,3,3,3),MIRR_4(NKMMAX,NKMMAX,3,3,3)
     &           ,MZAZB(NKMMAX,NKMMAX,3,3),MZBZA(NKMMAX,NKMMAX,3,3)
C
C Local variables
C
      INTEGER I,IA_ERR,IFIL_LHSA,IFIL_LHSB,IKM,IM,IPOL,IRTOP,
     &        IRUN_CHECK_BXAXW,J,JPOL,LAMA1,LAMA2,LAMB3,LAMB4,LAMQ,LM,N
      COMPLEX*16 IJAJB(:),IJAZB(:),IMC,IZAJB(:),IZAZB(:),JFLA(:,:,:),
     &           JFLB(:,:,:),JFRA(:,:,:),JFRB(:,:,:),JGLA(:,:,:),
     &           JGLB(:,:,:),JGRA(:,:,:),JGRB(:,:,:),MJBJA(:,:,:,:),
     &           MJBZA(:,:,:,:),MZBJA(:,:,:,:),PAB,PBA,SJAJB,SJAZB,
     &           SZAJB,SZAZB,ZFLA(:,:,:),ZFLB(:,:,:),ZFRA(:,:,:),
     &           ZFRB(:,:,:),ZGLA(:,:,:),ZGLB(:,:,:),ZGRA(:,:,:),
     &           ZGRB(:,:,:)
      LOGICAL K_SELECTION_RULES
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IZAZB,IZAJB,IJAZB,IJAJB
      ALLOCATABLE JFLA,JGLA,ZGLA,ZFLA,JGRB,JFRB,ZFRB,ZGRB
      ALLOCATABLE JFRA,JGRA,ZGRA,ZFRA,JGLB,JFLB,ZFLB,ZGLB
      ALLOCATABLE MZBJA,MJBJA,MJBZA
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
C        RTOP = R(IRTOP,IM)
C
         IF ( .NOT.ALLOCATED(A_SIG_RLM) ) THEN
C
            ALLOCATE (A_SIG_RLM(NKM_EXT,NKM_EXT,3,NLM_AME_RLM_EXT))
C
            N = NKM_EXT
C
            DO IPOL = 1,3
C
               DO LM = 1,NLM_AME_RLM_EXT
C
                  A_SIG_RLM(1:N,1:N,IPOL,LM)
     &               = MATMUL(AME_G(1:N,1:N,IPOL,ISMT),
     &               AG_RGNT(1:N,1:N,LM))
C
               END DO
C
            END DO
C
         END IF
      ELSE
         IRTOP = JRWS(IM)
C        RTOP = RWS(IM)
      END IF
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCC       IMC = CI*0.5D0*C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( .NOT.PUBLIC_VERSION ) WRITE (*,*) 
     &                          'XXXXXXXXXX STATEMENT TO BE REMOVED in '
     &                          ,ROUTINE
      IMC = 1D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C      IF ( ME_CC_BRA_RWF ) THEN
C         PAB = (C/(C**2+DCONJG(ERYDA)+ERYDB))*RTOP**2
C         PBA = (C/(C**2+ERYDA+DCONJG(ERYDB)))*RTOP**2
C      ELSE
C         PAB = (C/(C**2+ERYDA+ERYDB))*RTOP**2
C         PBA = (C/(C**2+ERYDA+ERYDB))*RTOP**2
C      END IF
C
      PAB = 0D0
      PBA = 0D0
C
C=======================================================================
C
      MZAZB(:,:,:,:) = C0
      MZBZA(:,:,:,:) = C0
      MIRR_2(:,:,:,:,:) = C0
      MIRR_3(:,:,:,:,:) = C0
      MIRR_4(:,:,:,:,:) = C0
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
      CALL WAVFUN_RAD_WGT_REL(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF_LA,'R')
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
C weight the wave function with the radial factor R2DRDI for integration
C
      CALL WAVFUN_RAD_WGT_REL(IT,ZGLB,ZFLB,JGLB,JFLB,NCPLWF_LB,'R')
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
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IRUN_CHECK_BXAXW = 0
 100  CONTINUE
      IRUN_CHECK_BXAXW = IRUN_CHECK_BXAXW + 1
C
      IF ( CHECK_BXAXW ) THEN
         PAB = 0D0
         PBA = 0D0
         IF ( .NOT.ALLOCATED(MZBJA) ) THEN
            ALLOCATE (MZBJA(NKMMAX,NKMMAX,3,3),MJBJA(NKMMAX,NKMMAX,3,3))
            ALLOCATE (MJBZA(NKMMAX,NKMMAX,3,3),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBJA')
         END IF
         MZBZA(:,:,:,:) = C0
         MJBZA(:,:,:,:) = C0
         MZBJA(:,:,:,:) = C0
         MJBJA(:,:,:,:) = C0
         CALL ME_ALF_B4A1_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMBOTB,
     &                           LAMTOPB,JGLB,JFLB,LAMBOTA,LAMTOPA,ZGRA,
     &                           ZFRA,MJBZA)
         CALL ME_ALF_B4A1_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMBOTB,
     &                           LAMTOPB,ZGLB,ZFLB,LAMBOTA,LAMTOPA,JGRA,
     &                           JFRA,MZBJA)
         CALL ME_ALF_B4A1_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMBOTB,
     &                           LAMTOPB,JGLB,JFLB,LAMBOTA,LAMTOPA,JGRA,
     &                           JFRA,MJBJA)
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      CALL ME_ALF_B4A1_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMBOTB,
     &                        LAMTOPB,ZGLB,ZFLB,LAMBOTA,LAMTOPA,ZGRA,
     &                        ZFRA,MZBZA)
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
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ LAMQ
            LOOP_LAMQ:DO LAMQ = 1,3
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
               LOOP_JPOL:DO JPOL = 1,NPOL
C
                  CALL ME_ALF_A2B3_NONDIP(LAMA2,LAMB3,JPOL,LAMQ,IM,IMC,
     &               C,ERYDA_KM,ERYDB,PAB,ZGLA,ZFLA,JGLA,JFLA,ZGRB,ZFRB,
     &               JGRB,JFRB,IZAZB,IZAJB,IJAZB,IJAJB,SZAZB,SZAJB,
     &               SJAZB,SJAJB,MZAZB,K_SELECTION_RULES)
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
                     IF ( CHECK_BXAXW )
     &                    CALL ME_CHECK_BXAXW(2,IRUN_CHECK_BXAXW,1,
     &                    MZBJA,MZBZA,LAMB4,LAMA1,IZAZB,SZAZB,IZAZB,
     &                    SZAZB,IJAZB,SJAZB,IJAZB,SJAZB,MIRR_2,LAMB4,
     &                    LAMB3,JPOL)
C
                     CALL ME_ALF_BXAXW_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,
     &                  PBA,LAMB4,ZGLB,ZFLB,ZGLB,ZFLB,LAMA1,ZGRA,ZFRA,
     &                  JFRA,JGRA,IZAZB,SZAZB,IJAZB,SJAZB,MIRR_2,LAMB4,
     &                  LAMB3,JPOL,LAMQ)
C
                     IF ( CHECK_BXAXW )
     &                    CALL ME_CHECK_BXAXW(2,IRUN_CHECK_BXAXW,2,
     &                    MZBJA,MZBZA,LAMB4,LAMA1,IZAZB,SZAZB,IZAZB,
     &                    SZAZB,IJAZB,SJAZB,IJAZB,SJAZB,MIRR_2,LAMB4,
     &                    LAMB3,JPOL)
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
                     IF ( CHECK_BXAXW )
     &                    CALL ME_CHECK_BXAXW(3,IRUN_CHECK_BXAXW,1,
     &                    MJBZA,MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,
     &                    SZAZB,IZAJB,SZAJB,IZAJB,SZAJB,MIRR_3,LAMA2,
     &                    LAMA1,JPOL)
C
                     CALL ME_ALF_BXAXW_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,
     &                  PBA,LAMB3,ZGLB,ZFLB,JGLB,JFLB,LAMA1,ZGRA,ZFRA,
     &                  ZFRA,ZGRA,IZAZB,SZAZB,IZAJB,SZAJB,MIRR_3,LAMA2,
     &                  LAMA1,JPOL,LAMQ)
C
                     IF ( CHECK_BXAXW )
     &                    CALL ME_CHECK_BXAXW(3,IRUN_CHECK_BXAXW,2,
     &                    MJBZA,MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,
     &                    SZAZB,IZAJB,SZAJB,IZAJB,SZAJB,MIRR_3,LAMA2,
     &                    LAMA1,JPOL)
C
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
                  IF ( CHECK_BXAXW )
     &                 CALL ME_CHECK_BXAXW(4,IRUN_CHECK_BXAXW,1,MJBJA,
     &                 MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAJB,
     &                 SJAJB,IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,JPOL)
C
                  CALL ME_ALF_BXAXW_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,
     &               LAMB3,ZGLB,ZFLB,JGLB,JFLB,LAMA1,ZGRA,ZFRA,JFRA,
     &               JGRA,IZAZB,SZAZB,IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,
     &               JPOL,LAMQ)
C
                  IF ( CHECK_BXAXW )
     &                 CALL ME_CHECK_BXAXW(4,IRUN_CHECK_BXAXW,2,MJBJA,
     &                 MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAJB,
     &                 SJAJB,IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,JPOL)
C
C***********************************************************************
C
               END DO LOOP_JPOL
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
            END DO LOOP_LAMQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ LAMQ
C
         END DO LOOP_LAMA2
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO LOOP_LAMB3
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHECK_BXAXW ) THEN
         IF ( IRUN_CHECK_BXAXW.EQ.1 ) GOTO 100
         STOP
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( .NOT.PUBLIC_VERSION ) WRITE (*,*) 
     &                          'XXXXXXXXXX STATEMENT TO BE REMOVED in '
     &                          ,ROUTINE
      DO LAMQ = 1,3
         DO I = 1,NKM
            DO J = 1,NKM
               MZBZA(I,J,2,LAMQ) = -MZAZB(J,I,2,LAMQ)
            END DO
         END DO
      END DO
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( IT.LE.-110 ) THEN
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
C*==me_alf_bxaxw_nondip.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_BXAXW_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMB,
     &                               ZGLB,ZFLB,JGLB,JFLB,LAMA,ZGRA,ZFRA,
     &                               JFRA,JGRA,IZAZB,SZAZB,IJAJB,SJAJB,
     &                               MIRR,I,J,JPOL,LAMQ)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate matrix element integral functions of the type          *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF_LB,IKMCPLWF_RA,NCPLWFMAX
      USE MOD_RMESH,ONLY:R2DRDI,JRWS,JRMT,JRCRI,FULLPOT,NRMAX,NSF,LMISF,
     &    FLMSF,KLMSF,R
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:B1_ADA,B2_ADA,NKM,NKMMAX,NCPLWF_LB,NCPLWF_RA,
     &    A_SIG_RLM,IMKM_IKM,NLM_AME_RLM_EXT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_BXAXW_NONDIP')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      INTEGER NPOL
      PARAMETER (NPOL=3)
      LOGICAL SURFACE_TERM
      PARAMETER (SURFACE_TERM=.FALSE.)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDB,IMC,PBA,SJAJB,SZAZB
      INTEGER I,IM,J,JPOL,LAMA,LAMB,LAMQ
      COMPLEX*16 ERYDA_KM(NKMMAX),IJAJB(NRMAX),IZAZB(NRMAX),
     &           JFLB(NRMAX,NCPLWFMAX,NKM),JFRA(NRMAX,NCPLWFMAX,NKM),
     &           JGLB(NRMAX,NCPLWFMAX,NKM),JGRA(NRMAX,NCPLWFMAX,NKM),
     &           MIRR(NKMMAX,NKMMAX,3,3,3),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 B1,B1_JFRA,B1_ZFRA,B2,B2_JGRA,B2_ZGRA,CAUX,FS(:),FV(:),
     &           IMC_B1,IMC_B2,PA,S2REG,SCNTR,VCNTR,VIRR,WIRR,WREG
      LOGICAL C_LT_EPS
      INTEGER IA,IA_ERR,IB,IFLAG2,IKMA,IKMB,IMKMA,IMKMB,IPOL,IR,IRCRIT,
     &        IRSF,IRTOP_SPHERE,ISF,LMSF
      REAL*8 QLNG,RARG
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FV,FS
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FV(NRMAX),FS(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FV')
C
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRCRIT = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRCRIT = JRWS(IM)
      END IF
C
      QLNG = ABS(DREAL(ERYDA_KM(LAMA)-ERYDB))/C
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
      DO IPOL = 1,NPOL
C
         DO IR = 1,IRCRIT
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
         DO IA = 1,NCPLWF_RA(LAMA)
            IKMA = IKMCPLWF_RA(IA,LAMA)
            IMKMA = IMKM_IKM(IKMA)
C
            DO IB = 1,NCPLWF_LB(LAMB)
               IKMB = IKMCPLWF_LB(IB,LAMB)
               IMKMB = IMKM_IKM(IKMB)
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
               IF ( RNON0(B1_ADA(IKMB,IKMA,IPOL,LAMQ)) .OR. 
     &              RNON0(B2_ADA(IKMB,IKMA,IPOL,LAMQ)) ) THEN
C
                  IFLAG2 = 1
C
                  IMC_B1 = IMC*B1_ADA(IKMB,IKMA,IPOL,LAMQ)
                  IMC_B2 = IMC*B2_ADA(IKMB,IKMA,IPOL,LAMQ)
C
                  DO IR = 1,IRTOP_SPHERE
C
                     B2_ZGRA = IMC_B2*ZGRA(IR,IA,LAMA)
                     B1_ZFRA = IMC_B1*ZFRA(IR,IA,LAMA)
                     B2_JGRA = IMC_B2*JGRA(IR,IA,LAMA)
                     B1_JFRA = IMC_B1*JFRA(IR,IA,LAMA)
C
                     WREG = JGLB(IR,IB,LAMB)*B1_JFRA - JFLB(IR,IB,LAMB)
     &                      *B2_JGRA
                     WREG = WREG*R(IR,IM)
                     WIRR = ZGLB(IR,IB,LAMB)*B1_ZFRA - ZFLB(IR,IB,LAMB)
     &                      *B2_ZGRA
                     WIRR = WIRR*R(IR,IM)
C
                     FV(IR) = FV(IR) + WREG*IZAZB(IR) + WIRR*IJAJB(IR)
                     FS(IR) = FS(IR) + WIRR
C
                  END DO
C
C------------------------------------------------------ ASA surface term
                  IF ( .NOT.FULLPOT .AND. SURFACE_TERM ) THEN
                     IR = IRTOP_SPHERE
                     PA = -0.5D0*PBA*(IMC_B1-IMC_B2)/R2DRDI(IR,IM)
                     S2REG = S2REG + (JGLB(IR,IB,LAMB)*JGRA(IR,IA,LAMA)
     &                       -JFLB(IR,IB,LAMB)*JFRA(IR,IA,LAMA))*PA
                  END IF
C-----------------------------------------------------------------------
C
               END IF
C
C-----------------------------------------------------------------------
C                interstitial regime in case of  FULLPOT
C                there should be NO surface contribution
C-----------------------------------------------------------------------
C
               LOOP_ISF:DO ISF = 1,NSF(IM)
C
                  IF ( .NOT.FULLPOT ) CYCLE LOOP_ISF
C
                  LMSF = LMISF(ISF,IM)
C
                  IF ( LMSF.GT.NLM_AME_RLM_EXT ) CYCLE LOOP_ISF
                  IF ( KLMSF(LMSF,IM).NE.1 ) CYCLE LOOP_ISF
C
                  B1 = A_SIG_RLM(IKMB,IMKMA,IPOL,LMSF)
                  B2 = A_SIG_RLM(IMKMB,IKMA,IPOL,LMSF)
C
                  IF ( C_LT_EPS(B1,NON0_TOL) .AND. C_LT_EPS(B2,NON0_TOL)
     &                 ) CYCLE LOOP_ISF
C
                  IFLAG2 = 1
C
                  IMC_B1 = IMC*B1
                  IMC_B2 = IMC*B2
C
                  DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                     IRSF = IR - IRTOP_SPHERE
C
                     B2_ZGRA = IMC_B2*ZGRA(IR,IA,LAMA)
                     B1_ZFRA = IMC_B1*ZFRA(IR,IA,LAMA)
                     B2_JGRA = IMC_B2*JGRA(IR,IA,LAMA)
                     B1_JFRA = IMC_B1*JFRA(IR,IA,LAMA)
C
                     WREG = JGLB(IR,IB,LAMB)*B1_JFRA - JFLB(IR,IB,LAMB)
     &                      *B2_JGRA
                     WREG = WREG*R(IR,IM)
                     WIRR = ZGLB(IR,IB,LAMB)*B1_ZFRA - ZFLB(IR,IB,LAMB)
     &                      *B2_ZGRA
                     WIRR = WIRR*R(IR,IM)
C
                     CAUX = WREG*IZAZB(IR) + WIRR*IJAJB(IR)
C
                     FV(IR) = FV(IR) + CAUX*FLMSF(IRSF,ISF,IM)
C
                  END DO
C
               END DO LOOP_ISF
C-----------------------------------------------------------------------
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
            CALL CRADINT(IM,FV,VCNTR)
C
C---------------------------------------- alpha - related   surface term
C
            IF ( FULLPOT ) THEN
C
               SCNTR = 0D0
C
            ELSE
C
               CALL CRADINT(IM,FS,VIRR)
C
               SCNTR = SJAJB*VIRR + (IZAZB(IRCRIT)+SZAZB)*S2REG
C
            END IF
C
            MIRR(I,J,IPOL,JPOL,LAMQ) = MIRR(I,J,IPOL,JPOL,LAMQ)
     &                                 + (VCNTR+SCNTR)*QLNG*QLNG
C
         END IF
C***********************************************************************
C
      END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
      END
C*==me_alf_b4a1_nondip.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_B4A1_NONDIP(IM,IMC,C,ERYDB,ERYDA_KM,PBA,LAMBOTB,
     &                              LAMTOPB,ZGLB,ZFLB,LAMBOTA,LAMTOPA,
     &                              ZGRA,ZFRA,MZBZA)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the matrix elements of the type                        *
C   *                                                                  *
C   *           MZBZA = < Z^+(E_b) | H_lam | Z(E_a) >                  *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator in the form     *
C   *                                                                  *
C   *             H_lam = mc alfa_lam = mc alfa * A_lam                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF_LB,IKMCPLWF_RA,NCPLWFMAX
      USE MOD_RMESH,ONLY:R2DRDI,JRWS,JRMT,JRCRI,FULLPOT,NRMAX,NSF,LMISF,
     &    FLMSF,KLMSF,R
      USE MOD_CONSTANTS,ONLY:C0,Y00
      USE MOD_ANGMOM,ONLY:AME_G,B1_ADA,B2_ADA,ISMT,NKM,NKMMAX,NCPLWF_LB,
     &    NCPLWF_RA,A_SIG_RLM,IMKM_IKM,NLM_AME_RLM_EXT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_B4A1_NONDIP')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      INTEGER NPOL
      PARAMETER (NPOL=3)
      LOGICAL CHECK,SURFACE_TERM
      PARAMETER (CHECK=.FALSE.,SURFACE_TERM=.FALSE.)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDB,IMC,PBA
      INTEGER IM,LAMBOTA,LAMBOTB,LAMTOPA,LAMTOPB
      COMPLEX*16 ERYDA_KM(NKMMAX),MZBZA(NKMMAX,NKMMAX,3,3),
     &           ZFLB(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGLB(NRMAX,NCPLWFMAX,NKM),ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 B1,B100,B2,B200,CAUX,FZBZA(:),IMC_B1,IMC_B2,PA,S2BAR,
     &           V2BAR,ZBZA
      LOGICAL C_LT_EPS
      INTEGER IA1,IA_ERR,IB4,IFLAG1,IKMA1,IKMB4,IMKMA1,IMKMB4,IPOL,IR,
     &        IRCRIT,IRSF,IRTOP_SPHERE,ISF,LAMA1,LAMB4,LAMQ,LMSF
      REAL*8 QLNG,RARG
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZBZA
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FZBZA(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FZBZA')
C
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRCRIT = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRCRIT = JRWS(IM)
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
            QLNG = ABS(DREAL(ERYDA_KM(LAMA1)-ERYDB))/C
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ LAMQ
            DO LAMQ = 1,3
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
               DO IPOL = 1,NPOL
C
                  DO IR = 1,IRCRIT
                     FZBZA(IR) = C0
                  END DO
                  S2BAR = C0
C
                  IFLAG1 = 0
C
C-----------------------------------------------------------------------
C                                                                IA1 IB4
                  DO IA1 = 1,NCPLWF_RA(LAMA1)
                     IKMA1 = IKMCPLWF_RA(IA1,LAMA1)
                     IMKMA1 = IMKM_IKM(IKMA1)
C
                     DO IB4 = 1,NCPLWF_LB(LAMB4)
                        IKMB4 = IKMCPLWF_LB(IB4,LAMB4)
                        IMKMB4 = IMKM_IKM(IKMB4)
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
                        IF ( RNON0(B1_ADA(IKMB4,IKMA1,IPOL,LAMQ)) .OR. 
     &                       RNON0(B2_ADA(IKMB4,IKMA1,IPOL,LAMQ)) ) THEN
C
                           IFLAG1 = 1
C
                           IMC_B1 = IMC*B1_ADA(IKMB4,IKMA1,IPOL,LAMQ)
                           IMC_B2 = IMC*B2_ADA(IKMB4,IKMA1,IPOL,LAMQ)
C
                           DO IR = 1,IRTOP_SPHERE
C
                              ZBZA = +IMC_B1*ZGLB(IR,IB4,LAMB4)
     &                               *ZFRA(IR,IA1,LAMA1)
     &                               - IMC_B2*ZFLB(IR,IB4,LAMB4)
     &                               *ZGRA(IR,IA1,LAMA1)
C
                              FZBZA(IR) = FZBZA(IR) + ZBZA*R(IR,IM)
C
                           END DO
C
C------------------------------------------------------ ASA surface term
                           IF ( .NOT.FULLPOT .AND. SURFACE_TERM ) THEN
                              IR = IRTOP_SPHERE
                              PA = -0.5D0*PBA*(IMC_B1-IMC_B2)
     &                             /R2DRDI(IR,IM)
                              S2BAR = S2BAR + 
     &                                (ZGRA(IR,IA1,LAMA1)*ZGLB(IR,IB4,
     &                                LAMB4)-ZFRA(IR,IA1,LAMA1)
     &                                *ZFLB(IR,IB4,LAMB4))*PA
                           END IF
C-----------------------------------------------------------------------
C
                        END IF
C
C-----------------------------------------------------------------------
C                interstitial regime in case of  FULLPOT
C                there should be NO surface contribution
C-----------------------------------------------------------------------
C
                        LOOP_ISF:DO ISF = 1,NSF(IM)
C
                           IF ( .NOT.FULLPOT ) CYCLE LOOP_ISF
C
                           LMSF = LMISF(ISF,IM)
C
                           IF ( LMSF.GT.NLM_AME_RLM_EXT ) CYCLE LOOP_ISF
                           IF ( KLMSF(LMSF,IM).NE.1 ) CYCLE LOOP_ISF
C
                           B1 = A_SIG_RLM(IKMB4,IMKMA1,IPOL,LMSF)
                           B2 = A_SIG_RLM(IMKMB4,IKMA1,IPOL,LMSF)
C
C???????????????????????????????????????????????????????????????????????
                           IF ( CHECK .AND. LMSF.EQ.1 ) THEN
                              B100 = B1/Y00
                              B200 = B2/Y00
                              IF ( ABS(AME_G(IKMB4,IMKMA1,IPOL,ISMT)-
     &                             B100).GT.1D-10 .OR. 
     &                             ABS(AME_G(IMKMB4,IKMA1,IPOL,ISMT)
     &                             -B200).GT.1D-10 ) THEN
                                 WRITE (6,*) 'check failed in routine ',
     &                                  ROUTINE
                                 WRITE (6,*) 'IKMB4,IKMA1,IPOL,LMSF:',
     &                                  IKMB4,IKMA1,IPOL,LMSF
                                 WRITE (6,*) 'A_SIG_RLM ',DREAL(B100),
     &                                  DREAL(B200)
                                 WRITE (6,*) 'B1 B2  ',
     &                                  AME_G(IKMB4,IMKMA1,IPOL,ISMT),
     &                                  AME_G(IMKMB4,IKMA1,IPOL,ISMT)
                              END IF
                           END IF
C???????????????????????????????????????????????????????????????????????
C
                           IF ( C_LT_EPS(B1,NON0_TOL) .AND. 
     &                          C_LT_EPS(B2,NON0_TOL) ) CYCLE LOOP_ISF
C
                           IFLAG1 = 1
C
                           IMC_B1 = IMC*B1
                           IMC_B2 = IMC*B2
C
                           DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                              IRSF = IR - IRTOP_SPHERE
C
                              CAUX = IMC_B1*ZGLB(IR,IB4,LAMB4)
     &                               *ZFRA(IR,IA1,LAMA1)
     &                               - IMC_B2*ZFLB(IR,IB4,LAMB4)
     &                               *ZGRA(IR,IA1,LAMA1)
C
                              FZBZA(IR) = FZBZA(IR)
     &                           + CAUX*FLMSF(IRSF,ISF,IM)*R(IR,IM)
C
                           END DO
C
                        END DO LOOP_ISF
C-----------------------------------------------------------------------
C
                     END DO
                  END DO
C                                                                IA1 IB4
C-----------------------------------------------------------------------
C
C#######################################################################
C
                  IF ( IFLAG1.EQ.1 ) THEN
C
                     CALL CRADINT(IM,FZBZA,V2BAR)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
                     IF ( FULLPOT ) S2BAR = 0D0
C
                     MZBZA(LAMB4,LAMA1,IPOL,LAMQ)
     &                  = MZBZA(LAMB4,LAMA1,IPOL,LAMQ) + (V2BAR+S2BAR)
     &                  *QLNG
C
                  END IF
C################################################################ IFLAG1
C
               END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
            END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ LAMQ
C
         END DO
C                                                      energy A -- LAMA1
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB4
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
Ccc      stop
      END
C*==me_alf_a2b3_nondip.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_A2B3_NONDIP(LAMA2,LAMB3,JPOL,LAMQ,IM,IMC,C,
     &                              ERYDA_KM,ERYDB,PAB,ZGLA,ZFLA,JGLA,
     &                              JFLA,ZGRB,ZFRB,JGRB,JFRB,IZAZB,
     &                              IZAJB,IJAZB,IJAJB,SZAZB,SZAJB,SJAZB,
     &                              SJAJB,MZAZB,K_SELECTION_RULES)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate matrix element integral functions of the type          *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *      SXAYB : corresponding surface contributions                 *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *     with    H_lam   the electric dipole operator in the form     *
C   *                                                                  *
C   *             H_lam = mc alfa_lam = mc alfa * A_lam                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IKMCPLWF_LA,IKMCPLWF_RB,NCPLWFMAX
      USE MOD_RMESH,ONLY:R2DRDI,JRWS,JRMT,JRCRI,FULLPOT,NRMAX,NSF,LMISF,
     &    FLMSF,KLMSF,R
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF_LA,NCPLWF_RB,A_SIG_RLM,
     &    IMKM_IKM,NLM_AME_RLM_EXT,B1_ADA,B2_ADA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_A2B3_NONDIP')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      LOGICAL SURFACE_TERM
      PARAMETER (SURFACE_TERM=.FALSE.)
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDB,IMC,PAB,SJAJB,SJAZB,SZAJB,SZAZB
      INTEGER IM,JPOL,LAMA2,LAMB3,LAMQ
      LOGICAL K_SELECTION_RULES
      COMPLEX*16 ERYDA_KM(NKMMAX),IJAJB(NRMAX),IJAZB(NRMAX),IZAJB(NRMAX)
     &           ,IZAZB(NRMAX),JFLA(NRMAX,NCPLWFMAX,NKM),
     &           JFRB(NRMAX,NCPLWFMAX,NKM),JGLA(NRMAX,NCPLWFMAX,NKM),
     &           JGRB(NRMAX,NCPLWFMAX,NKM),MZAZB(NKMMAX,NKMMAX,3,3),
     &           ZFLA(NRMAX,NCPLWFMAX,NKM),ZFRB(NRMAX,NCPLWFMAX,NKM),
     &           ZGLA(NRMAX,NCPLWFMAX,NKM),ZGRB(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 B1,B1_JGLA,B1_ZGLA,B2,B2_JFLA,B2_ZFLA,CAUX,FJAJB(:),
     &           FJAZB(:),FZAJB(:),FZAZB(:),IMC_B1,IMC_B2,PA
      LOGICAL C_LT_EPS
      REAL*8 FLMSF_R,QLNG,RARG
      INTEGER IA2,IA_ERR,IB3,IFLAG1,IKMA2,IKMB3,IMKMA2,IMKMB3,IR,IRCRIT,
     &        IRSF,IRTOP_SPHERE,ISF,LMSF
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZAZB,FZAJB,FJAZB,FJAJB
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (FZAZB(NRMAX),FZAJB(NRMAX),FJAZB(NRMAX),FJAJB(NRMAX),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FZAZB')
C
      IF ( FULLPOT ) THEN
         IRTOP_SPHERE = JRMT(IM)
         IRCRIT = JRCRI(IM)
      ELSE
         IRTOP_SPHERE = JRWS(IM)
         IRCRIT = JRWS(IM)
      END IF
C
      K_SELECTION_RULES = .FALSE.
C
      SZAZB = C0
      SZAJB = C0
      SJAZB = C0
      SJAJB = C0
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C
      DO IR = 1,IRCRIT
         FZAZB(IR) = C0
         FZAJB(IR) = C0
         FJAZB(IR) = C0
         FJAJB(IR) = C0
      END DO
C
      IFLAG1 = 0
C
C-----------------------------------------------------------------------
C                                                                IA2 IB3
      DO IA2 = 1,NCPLWF_LA(LAMA2)
         IKMA2 = IKMCPLWF_LA(IA2,LAMA2)
         IMKMA2 = IMKM_IKM(IKMA2)
C
         DO IB3 = 1,NCPLWF_RB(LAMB3)
            IKMB3 = IKMCPLWF_RB(IB3,LAMB3)
            IMKMB3 = IMKM_IKM(IKMB3)
C
            QLNG = ABS(DREAL(ERYDA_KM(LAMA2)-ERYDB))/C
C
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
            IF ( RNON0(B1_ADA(IKMA2,IKMB3,JPOL,LAMQ)) .OR. 
     &           RNON0(B2_ADA(IKMA2,IKMB3,JPOL,LAMQ)) ) THEN
C
               IFLAG1 = 1
C
               IMC_B1 = IMC*B1_ADA(IKMA2,IKMB3,JPOL,LAMQ)
               IMC_B2 = IMC*B2_ADA(IKMA2,IKMB3,JPOL,LAMQ)
C
               DO IR = 1,IRTOP_SPHERE
C
                  B1_ZGLA = IMC_B1*ZGLA(IR,IA2,LAMA2)*R(IR,IM)
                  B2_ZFLA = IMC_B2*ZFLA(IR,IA2,LAMA2)*R(IR,IM)
                  B1_JGLA = IMC_B1*JGLA(IR,IA2,LAMA2)*R(IR,IM)
                  B2_JFLA = IMC_B2*JFLA(IR,IA2,LAMA2)*R(IR,IM)
C
                  FZAZB(IR) = FZAZB(IR) + B1_ZGLA*ZFRB(IR,IB3,LAMB3)
     &                        - B2_ZFLA*ZGRB(IR,IB3,LAMB3)
                  FZAJB(IR) = FZAJB(IR) + B1_ZGLA*JFRB(IR,IB3,LAMB3)
     &                        - B2_ZFLA*JGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + B1_JGLA*ZFRB(IR,IB3,LAMB3)
     &                        - B2_JFLA*ZGRB(IR,IB3,LAMB3)
                  FJAJB(IR) = FJAJB(IR) + B1_JGLA*JFRB(IR,IB3,LAMB3)
     &                        - B2_JFLA*JGRB(IR,IB3,LAMB3)
C
               END DO
C
C------------------------------------------------------ ASA surface term
               IF ( .NOT.FULLPOT .AND. SURFACE_TERM ) THEN
                  IR = IRTOP_SPHERE
                  PA = -0.5D0*PAB*(IMC_B1-IMC_B2)/R2DRDI(IR,IM)
C
                  SZAZB = SZAZB + (ZGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                    -ZFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3))*PA
C
                  SZAJB = SZAJB + (ZGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                    -ZFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3))*PA
C
                  SJAZB = SJAZB + (JGLA(IR,IA2,LAMA2)*ZGRB(IR,IB3,LAMB3)
     &                    -JFLA(IR,IA2,LAMA2)*ZFRB(IR,IB3,LAMB3))*PA
C
                  SJAJB = SJAJB + (JGLA(IR,IA2,LAMA2)*JGRB(IR,IB3,LAMB3)
     &                    -JFLA(IR,IA2,LAMA2)*JFRB(IR,IB3,LAMB3))*PA
               END IF
C-----------------------------------------------------------------------
C
            END IF
C
C-----------------------------------------------------------------------
C                interstitial regime in case of  FULLPOT
C                there should be NO surface contribution
C-----------------------------------------------------------------------
C
            LOOP_ISF:DO ISF = 1,NSF(IM)
C
               IF ( .NOT.FULLPOT ) CYCLE LOOP_ISF
C
               LMSF = LMISF(ISF,IM)
C
               IF ( LMSF.GT.NLM_AME_RLM_EXT ) CYCLE LOOP_ISF
               IF ( KLMSF(LMSF,IM).NE.1 ) CYCLE LOOP_ISF
C
               B1 = A_SIG_RLM(IKMA2,IMKMB3,JPOL,LMSF)
               B2 = A_SIG_RLM(IMKMA2,IKMB3,JPOL,LMSF)
C
               IF ( C_LT_EPS(B1,NON0_TOL) .AND. C_LT_EPS(B2,NON0_TOL) )
     &              CYCLE LOOP_ISF
C
               IFLAG1 = 1
C
               IMC_B1 = IMC*B1
               IMC_B2 = IMC*B2
C
               DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                  IRSF = IR - IRTOP_SPHERE
C
                  B1_ZGLA = IMC_B1*ZGLA(IR,IA2,LAMA2)
                  B2_ZFLA = IMC_B2*ZFLA(IR,IA2,LAMA2)
                  B1_JGLA = IMC_B1*JGLA(IR,IA2,LAMA2)
                  B2_JFLA = IMC_B2*JFLA(IR,IA2,LAMA2)
C
                  CAUX = B1_ZGLA*ZFRB(IR,IB3,LAMB3)
     &                   - B2_ZFLA*ZGRB(IR,IB3,LAMB3)
C
                  FLMSF_R = FLMSF(IRSF,ISF,IM)*R(IR,IM)
C
                  FZAZB(IR) = FZAZB(IR) + CAUX*FLMSF_R
C
                  CAUX = B1_ZGLA*JFRB(IR,IB3,LAMB3)
     &                   - B2_ZFLA*JGRB(IR,IB3,LAMB3)
C
                  FZAJB(IR) = FZAJB(IR) + CAUX*FLMSF_R
C
                  CAUX = B1_JGLA*ZFRB(IR,IB3,LAMB3)
     &                   - B2_JFLA*ZGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + CAUX*FLMSF_R
C
                  CAUX = B1_JGLA*JFRB(IR,IB3,LAMB3)
     &                   - B2_JFLA*JGRB(IR,IB3,LAMB3)
C
                  FJAJB(IR) = FJAJB(IR) + CAUX*FLMSF_R
C
               END DO
C
            END DO LOOP_ISF
C-----------------------------------------------------------------------
C
         END DO
      END DO
C                                                                IA2 IB3
C-----------------------------------------------------------------------
C
C#######################################################################
C
      IF ( IFLAG1.EQ.1 ) THEN
C
         CALL CRADINT_R(IM,FZAZB,IZAZB)
C
         CALL CRADINT_INW_R(IM,FZAJB,IZAJB)
C
         CALL CRADINT_INW_R(IM,FJAZB,IJAZB)
C
         CALL CRADINT_INW_R(IM,FJAJB,IJAJB)
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
         IF ( FULLPOT ) THEN
            SZAZB = 0D0
            SZAJB = 0D0
            SJAZB = 0D0
            SJAJB = 0D0
         END IF
C
         MZAZB(LAMA2,LAMB3,JPOL,LAMQ) = MZAZB(LAMA2,LAMB3,JPOL,LAMQ)
     &                                  + (IZAZB(IRCRIT)+SZAZB)*QLNG
C
         K_SELECTION_RULES = .TRUE.
C
      END IF
C################################################################ IFLAG1
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP JPOL
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
