C*==me_alf_alf.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_ALF(LAMBOTA,LAMTOPA,IFIL_RHSA,ERYDA,LAMBOTB,
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
      USE MOD_RMESH,ONLY:R,RWS,JRWS,NRMAX,JRCRI,FULLPOT
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF_LA,NCPLWF_LB,NCPLWF_RA,
     &    NCPLWF_RB,AME_G,AG_RGNT,NKM_EXT,NLM_AME_RLM_EXT,ISMT,
     &    A_SIG_RLM,NPOL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_ALF')
      LOGICAL CHECK_BXAXW
      PARAMETER (CHECK_BXAXW=.FALSE.)
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
      INTEGER IA_ERR,IFIL_LHSA,IFIL_LHSB,IKM,IM,IPOL,IRTOP,
     &        IRUN_CHECK_BXAXW,JPOL,KIRR,LAMA1,LAMA2,LAMB3,LAMB4,LM,N
      COMPLEX*16 IJAJB(:),IJAZB(:),IMC,IZAJB(:),IZAZB(:),JFLA(:,:,:),
     &           JFLB(:,:,:),JFRA(:,:,:),JFRB(:,:,:),JGLA(:,:,:),
     &           JGLB(:,:,:),JGRA(:,:,:),JGRB(:,:,:),MJBJA(:,:,:),
     &           MJBZA(:,:,:),MZBJA(:,:,:),PAB,PBA,SJAJB,SJAZB,SZAJB,
     &           SZAZB,ZFLA(:,:,:),ZFLB(:,:,:),ZFRA(:,:,:),ZFRB(:,:,:),
     &           ZGLA(:,:,:),ZGLB(:,:,:),ZGRA(:,:,:),ZGRB(:,:,:)
      LOGICAL K_SELECTION_RULES
      REAL*8 RTOP
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
C         read irregular wave functions only if necessary
C
      IF ( K_CALC_ME.GE.3 ) THEN
         KIRR = 1
      ELSE
         KIRR = 0
      END IF
C
C=======================================================================
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
         RTOP = R(IRTOP,IM)
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
         RTOP = RWS(IM)
      END IF
C
      IMC = CI*0.5D0*C
C
      IF ( ME_CC_BRA_RWF ) THEN
         PAB = (C/(C**2+DCONJG(ERYDA)+ERYDB))*RTOP**2
         PBA = (C/(C**2+ERYDA+DCONJG(ERYDB)))*RTOP**2
      ELSE
         PAB = (C/(C**2+ERYDA+ERYDB))*RTOP**2
         PBA = (C/(C**2+ERYDA+ERYDB))*RTOP**2
      END IF
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
C          IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r
C
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for LA
C
      CALL WAVFUN_READ_REL(IFIL_LHSA,IT,KIRR,ZGLA,ZFLA,JGLA,JFLA,IRTOP,
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
      CALL WAVFUN_READ_REL(IFIL_RHSB,IT,KIRR,ZGRB,ZFRB,JGRB,JFRB,IRTOP,
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
      CALL WAVFUN_READ_REL(IFIL_LHSB,IT,KIRR,ZGLB,ZFLB,JGLB,JFLB,IRTOP,
     &                     NCPLWF_LB,IKMCPLWF_LB)
C
C weight the wave function with the radial factor R2DRDI for integration
C
      CALL WAVFUN_RAD_WGT_REL(IT,ZGLB,ZFLB,JGLB,JFLB,NCPLWF_LB,'R')
C
C ----------------------------------------- read in wavefunctions for RA
C
      CALL WAVFUN_READ_REL(IFIL_RHSA,IT,KIRR,ZGRA,ZFRA,JGRA,JFRA,IRTOP,
     &                     NCPLWF_RA,IKMCPLWF_RA)
C
C-------------- take complex conjugate of wave function ZGRB if requested
C
Cc      IF ( ME_CC_BRA_RWF ) THEN
Cc         DO IKM = 1,NKM
Cc            ZGLB(:,:,IKM) = DCONJG(ZGLB(:,:,IKM))
Cc            ZFLB(:,:,IKM) = DCONJG(ZFLB(:,:,IKM))
Cc            JGLB(:,:,IKM) = DCONJG(JGLB(:,:,IKM))
Cc            JFLB(:,:,IKM) = DCONJG(JFLB(:,:,IKM))
Cc         END DO
Cc      END IF
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
            ALLOCATE (MZBJA(NKMMAX,NKMMAX,3),MJBJA(NKMMAX,NKMMAX,3))
            ALLOCATE (MJBZA(NKMMAX,NKMMAX,3),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBJA')
         END IF
         MZBZA(:,:,:) = C0
         MJBZA(:,:,:) = C0
         MZBJA(:,:,:) = C0
         MJBJA(:,:,:) = C0
         CALL ME_ALF_B4A1(IM,IMC,PBA,LAMBOTB,LAMTOPB,JGLB,JFLB,LAMBOTA,
     &                    LAMTOPA,ZGRA,ZFRA,MJBZA)
         CALL ME_ALF_B4A1(IM,IMC,PBA,LAMBOTB,LAMTOPB,ZGLB,ZFLB,LAMBOTA,
     &                    LAMTOPA,JGRA,JFRA,MZBJA)
         CALL ME_ALF_B4A1(IM,IMC,PBA,LAMBOTB,LAMTOPB,JGLB,JFLB,LAMBOTA,
     &                    LAMTOPA,JGRA,JFRA,MJBJA)
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      CALL ME_ALF_B4A1(IM,IMC,PBA,LAMBOTB,LAMTOPB,ZGLB,ZFLB,LAMBOTA,
     &                 LAMTOPA,ZGRA,ZFRA,MZBZA)
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
               CALL ME_ALF_A2B3(LAMA2,LAMB3,JPOL,IM,IMC,PAB,ZGLA,ZFLA,
     &                          JGLA,JFLA,ZGRB,ZFRB,JGRB,JFRB,IZAZB,
     &                          IZAJB,IJAZB,IJAJB,SZAZB,SZAJB,SJAZB,
     &                          SJAJB,MZAZB,K_SELECTION_RULES)
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
     &                 CALL ME_CHECK_BXAXW(2,IRUN_CHECK_BXAXW,1,MZBJA,
     &                 MZBZA,LAMB4,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAZB,
     &                 SJAZB,IJAZB,SJAZB,MIRR_2,LAMB4,LAMB3,JPOL)
C
                  CALL ME_ALF_BXAXW(IM,IMC,PBA,LAMB4,ZGLB,ZFLB,ZGLB,
     &                              ZFLB,LAMA1,ZGRA,ZFRA,JFRA,JGRA,
     &                              IZAZB,SZAZB,IJAZB,SJAZB,MIRR_2,
     &                              LAMB4,LAMB3,JPOL)
C
                  IF ( CHECK_BXAXW )
     &                 CALL ME_CHECK_BXAXW(2,IRUN_CHECK_BXAXW,2,MZBJA,
     &                 MZBZA,LAMB4,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAZB,
     &                 SJAZB,IJAZB,SJAZB,MIRR_2,LAMB4,LAMB3,JPOL)
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
     &                 CALL ME_CHECK_BXAXW(3,IRUN_CHECK_BXAXW,1,MJBZA,
     &                 MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IZAJB,
     &                 SZAJB,IZAJB,SZAJB,MIRR_3,LAMA2,LAMA1,JPOL)
C
                  CALL ME_ALF_BXAXW(IM,IMC,PBA,LAMB3,ZGLB,ZFLB,JGLB,
     &                              JFLB,LAMA1,ZGRA,ZFRA,ZFRA,ZGRA,
     &                              IZAZB,SZAZB,IZAJB,SZAJB,MIRR_3,
     &                              LAMA2,LAMA1,JPOL)
C
                  IF ( CHECK_BXAXW )
     &                 CALL ME_CHECK_BXAXW(3,IRUN_CHECK_BXAXW,2,MJBZA,
     &                 MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IZAJB,
     &                 SZAJB,IZAJB,SZAJB,MIRR_3,LAMA2,LAMA1,JPOL)
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
     &              CALL ME_CHECK_BXAXW(4,IRUN_CHECK_BXAXW,1,MJBJA,
     &              MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAJB,
     &              SJAJB,IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,JPOL)
C
               CALL ME_ALF_BXAXW(IM,IMC,PBA,LAMB3,ZGLB,ZFLB,JGLB,JFLB,
     &                           LAMA1,ZGRA,ZFRA,JFRA,JGRA,IZAZB,SZAZB,
     &                           IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,JPOL)
C
               IF ( CHECK_BXAXW )
     &              CALL ME_CHECK_BXAXW(4,IRUN_CHECK_BXAXW,2,MJBJA,
     &              MZBZA,LAMB3,LAMA1,IZAZB,SZAZB,IZAZB,SZAZB,IJAJB,
     &              SJAJB,IJAJB,SJAJB,MIRR_4,LAMA1,LAMB3,JPOL)
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
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHECK_BXAXW ) THEN
         IF ( IRUN_CHECK_BXAXW.EQ.1 ) GOTO 100
         STOP
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
C*==me_alf_bxaxw.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_BXAXW(IM,IMC,PBA,LAMB,ZGLB,ZFLB,JGLB,JFLB,LAMA,
     &                        ZGRA,ZFRA,JFRA,JGRA,IZAZB,SZAZB,IJAJB,
     &                        SJAJB,MIRR,I,J,JPOL)
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
     &    FLMSF,KLMSF
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:AME_G,ISMT,NKM,NKMMAX,NCPLWF_LB,NCPLWF_RA,
     &    A_SIG_RLM,IMKM_IKM,NLM_AME_RLM_EXT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_BXAXW')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      INTEGER NPOL
      PARAMETER (NPOL=3)
C
C Dummy arguments
C
      INTEGER I,IM,J,JPOL,LAMA,LAMB
      COMPLEX*16 IMC,PBA,SJAJB,SZAZB
      COMPLEX*16 IJAJB(NRMAX),IZAZB(NRMAX),JFLB(NRMAX,NCPLWFMAX,NKM),
     &           JFRA(NRMAX,NCPLWFMAX,NKM),JGLB(NRMAX,NCPLWFMAX,NKM),
     &           JGRA(NRMAX,NCPLWFMAX,NKM),MIRR(NKMMAX,NKMMAX,3,3),
     &           ZFLB(NRMAX,NCPLWFMAX,NKM),ZFRA(NRMAX,NCPLWFMAX,NKM),
     &           ZGLB(NRMAX,NCPLWFMAX,NKM),ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 A1,A1_JFRA,A1_ZFRA,A2,A2_JGRA,A2_ZGRA,CAUX,FS(:),FV(:),
     &           IMC_A1,IMC_A2,PA,S2REG,SCNTR,VCNTR,VIRR,WIRR,WREG
      LOGICAL C_LT_EPS
      INTEGER IA,IA_ERR,IB,IFLAG2,IKMA,IKMB,IMKMA,IMKMB,IPOL,IR,IRCRIT,
     &        IRSF,IRTOP_SPHERE,ISF,LMSF
      REAL*8 RARG
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
               IF ( RNON0(AME_G(IKMB,IMKMA,IPOL,ISMT)) .OR. 
     &              RNON0(AME_G(IMKMB,IKMA,IPOL,ISMT)) ) THEN
C
                  IFLAG2 = 1
C
                  IMC_A1 = IMC*AME_G(IKMB,IMKMA,IPOL,ISMT)
                  IMC_A2 = IMC*AME_G(IMKMB,IKMA,IPOL,ISMT)
C
                  DO IR = 1,IRTOP_SPHERE
C
                     A2_ZGRA = IMC_A2*ZGRA(IR,IA,LAMA)
                     A1_ZFRA = IMC_A1*ZFRA(IR,IA,LAMA)
                     A2_JGRA = IMC_A2*JGRA(IR,IA,LAMA)
                     A1_JFRA = IMC_A1*JFRA(IR,IA,LAMA)
C
                     WREG = JGLB(IR,IB,LAMB)*A1_JFRA - JFLB(IR,IB,LAMB)
     &                      *A2_JGRA
                     WIRR = ZGLB(IR,IB,LAMB)*A1_ZFRA - ZFLB(IR,IB,LAMB)
     &                      *A2_ZGRA
C
                     FV(IR) = FV(IR) + WREG*IZAZB(IR) + WIRR*IJAJB(IR)
                     FS(IR) = FS(IR) + WIRR
C
                  END DO
C
C------------------------------------------------------ ASA surface term
                  IF ( .NOT.FULLPOT ) THEN
                     IR = IRTOP_SPHERE
                     PA = -0.5D0*PBA*(IMC_A1-IMC_A2)/R2DRDI(IR,IM)
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
                  A1 = A_SIG_RLM(IKMB,IMKMA,IPOL,LMSF)
                  A2 = A_SIG_RLM(IMKMB,IKMA,IPOL,LMSF)
C
                  IF ( C_LT_EPS(A1,NON0_TOL) .AND. C_LT_EPS(A2,NON0_TOL)
     &                 ) CYCLE LOOP_ISF
C
                  IFLAG2 = 1
C
                  IMC_A1 = IMC*A1
                  IMC_A2 = IMC*A2
C
                  DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                     IRSF = IR - IRTOP_SPHERE
C
                     A2_ZGRA = IMC_A2*ZGRA(IR,IA,LAMA)
                     A1_ZFRA = IMC_A1*ZFRA(IR,IA,LAMA)
                     A2_JGRA = IMC_A2*JGRA(IR,IA,LAMA)
                     A1_JFRA = IMC_A1*JFRA(IR,IA,LAMA)
C
                     WREG = JGLB(IR,IB,LAMB)*A1_JFRA - JFLB(IR,IB,LAMB)
     &                      *A2_JGRA
                     WIRR = ZGLB(IR,IB,LAMB)*A1_ZFRA - ZFLB(IR,IB,LAMB)
     &                      *A2_ZGRA
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
            MIRR(I,J,IPOL,JPOL) = MIRR(I,J,IPOL,JPOL) + VCNTR + SCNTR
C
         END IF
C***********************************************************************
C
      END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPOL
C
      END
C*==me_alf_b4a1.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_B4A1(IM,IMC,PBA,LAMBOTB,LAMTOPB,ZGLB,ZFLB,
     &                       LAMBOTA,LAMTOPA,ZGRA,ZFRA,MZBZA)
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
     &    FLMSF,KLMSF
      USE MOD_CONSTANTS,ONLY:C0,Y00
      USE MOD_ANGMOM,ONLY:AME_G,ISMT,NKM,NKMMAX,NCPLWF_LB,NCPLWF_RA,
     &    A_SIG_RLM,IMKM_IKM,NLM_AME_RLM_EXT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_B4A1')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      INTEGER NPOL
      PARAMETER (NPOL=3)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      INTEGER IM,LAMBOTA,LAMBOTB,LAMTOPA,LAMTOPB
      COMPLEX*16 IMC,PBA
      COMPLEX*16 MZBZA(NKMMAX,NKMMAX,3),ZFLB(NRMAX,NCPLWFMAX,NKM),
     &           ZFRA(NRMAX,NCPLWFMAX,NKM),ZGLB(NRMAX,NCPLWFMAX,NKM),
     &           ZGRA(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 A1,A100,A2,A200,CAUX,FZBZA(:),IMC_A1,IMC_A2,PA,S2BAR,
     &           V2BAR
      LOGICAL C_LT_EPS
      INTEGER IA1,IA_ERR,IB4,IFLAG1,IKMA1,IKMB4,IMKMA1,IMKMB4,IPOL,IR,
     &        IRCRIT,IRSF,IRTOP_SPHERE,ISF,LAMA1,LAMB4,LMSF
      REAL*8 RARG
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
                     IF ( RNON0(AME_G(IKMB4,IMKMA1,IPOL,ISMT)) .OR. 
     &                    RNON0(AME_G(IMKMB4,IKMA1,IPOL,ISMT)) ) THEN
C
                        IFLAG1 = 1
C
                        IMC_A1 = IMC*AME_G(IKMB4,IMKMA1,IPOL,ISMT)
                        IMC_A2 = IMC*AME_G(IMKMB4,IKMA1,IPOL,ISMT)
C
                        DO IR = 1,IRTOP_SPHERE
C
                           FZBZA(IR) = FZBZA(IR)
     &                                 + IMC_A1*ZGLB(IR,IB4,LAMB4)
     &                                 *ZFRA(IR,IA1,LAMA1)
     &                                 - IMC_A2*ZFLB(IR,IB4,LAMB4)
     &                                 *ZGRA(IR,IA1,LAMA1)
C
                        END DO
C
C------------------------------------------------------ ASA surface term
                        IF ( .NOT.FULLPOT ) THEN
                           IR = IRTOP_SPHERE
                           PA = -0.5D0*PBA*(IMC_A1-IMC_A2)/R2DRDI(IR,IM)
                           S2BAR = S2BAR + 
     &                             (ZGRA(IR,IA1,LAMA1)*ZGLB(IR,IB4,
     &                             LAMB4)-ZFRA(IR,IA1,LAMA1)
     &                             *ZFLB(IR,IB4,LAMB4))*PA
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
                        A1 = A_SIG_RLM(IKMB4,IMKMA1,IPOL,LMSF)
                        A2 = A_SIG_RLM(IMKMB4,IKMA1,IPOL,LMSF)
C
C???????????????????????????????????????????????????????????????????????
                        IF ( CHECK .AND. LMSF.EQ.1 ) THEN
                           A100 = A1/Y00
                           A200 = A2/Y00
                           IF ( ABS(AME_G(IKMB4,IMKMA1,IPOL,ISMT)-A100)
     &                          .GT.1D-10 .OR. 
     &                          ABS(AME_G(IMKMB4,IKMA1,IPOL,ISMT)-A200)
     &                          .GT.1D-10 ) THEN
                              WRITE (6,*) 'check failed in routine ',
     &                               ROUTINE
                              WRITE (6,*) 'IKMB4,IKMA1,IPOL,LMSF:',
     &                               IKMB4,IKMA1,IPOL,LMSF
                              WRITE (6,*) 'A_SIG_RLM ',DREAL(A100),
     &                               DREAL(A200)
                              WRITE (6,*) 'A1 A2  ',
     &                               AME_G(IKMB4,IMKMA1,IPOL,ISMT),
     &                               AME_G(IMKMB4,IKMA1,IPOL,ISMT)
                           END IF
                        END IF
C???????????????????????????????????????????????????????????????????????
C
                        IF ( C_LT_EPS(A1,NON0_TOL) .AND. 
     &                       C_LT_EPS(A2,NON0_TOL) ) CYCLE LOOP_ISF
C
                        IFLAG1 = 1
C
                        IMC_A1 = IMC*A1
                        IMC_A2 = IMC*A2
C
                        DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                           IRSF = IR - IRTOP_SPHERE
C
                           CAUX = IMC_A1*ZGLB(IR,IB4,LAMB4)
     &                            *ZFRA(IR,IA1,LAMA1)
     &                            - IMC_A2*ZFLB(IR,IB4,LAMB4)
     &                            *ZGRA(IR,IA1,LAMA1)
C
                           FZBZA(IR) = FZBZA(IR)
     &                                 + CAUX*FLMSF(IRSF,ISF,IM)
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
                  MZBZA(LAMB4,LAMA1,IPOL) = MZBZA(LAMB4,LAMA1,IPOL)
     &               + V2BAR + S2BAR
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
C*==me_alf_a2b3.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_ALF_A2B3(LAMA2,LAMB3,JPOL,IM,IMC,PAB,ZGLA,ZFLA,JGLA,
     &                       JFLA,ZGRB,ZFRB,JGRB,JFRB,IZAZB,IZAJB,IJAZB,
     &                       IJAJB,SZAZB,SZAJB,SJAZB,SJAJB,MZAZB,
     &                       K_SELECTION_RULES)
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
     &    FLMSF,KLMSF
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF_LA,NCPLWF_RB,A_SIG_RLM,
     &    IMKM_IKM,NLM_AME_RLM_EXT,AME_G,ISMT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_ALF_A2B3')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER IM,JPOL,LAMA2,LAMB3
      COMPLEX*16 IMC,PAB,SJAJB,SJAZB,SZAJB,SZAZB
      LOGICAL K_SELECTION_RULES
      COMPLEX*16 IJAJB(NRMAX),IJAZB(NRMAX),IZAJB(NRMAX),IZAZB(NRMAX),
     &           JFLA(NRMAX,NCPLWFMAX,NKM),JFRB(NRMAX,NCPLWFMAX,NKM),
     &           JGLA(NRMAX,NCPLWFMAX,NKM),JGRB(NRMAX,NCPLWFMAX,NKM),
     &           MZAZB(NKMMAX,NKMMAX,3),ZFLA(NRMAX,NCPLWFMAX,NKM),
     &           ZFRB(NRMAX,NCPLWFMAX,NKM),ZGLA(NRMAX,NCPLWFMAX,NKM),
     &           ZGRB(NRMAX,NCPLWFMAX,NKM)
C
C Local variables
C
      COMPLEX*16 A1,A1_JGLA,A1_ZGLA,A2,A2_JFLA,A2_ZFLA,CAUX,FJAJB(:),
     &           FJAZB(:),FZAJB(:),FZAZB(:),IMC_A1,IMC_A2,PA
      LOGICAL C_LT_EPS
      INTEGER IA2,IA_ERR,IB3,IFLAG1,IKMA2,IKMB3,IMKMA2,IMKMB3,IR,IRCRIT,
     &        IRSF,IRTOP_SPHERE,ISF,LMSF
      REAL*8 RARG
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
C-----------------------------------------------------------------------
C                    regime of muffin tin or ASA sphere
C-----------------------------------------------------------------------
C
            IF ( RNON0(AME_G(IKMA2,IMKMB3,JPOL,ISMT)) .OR. 
     &           RNON0(AME_G(IMKMA2,IKMB3,JPOL,ISMT)) ) THEN
C
               IFLAG1 = 1
C
               IMC_A1 = IMC*AME_G(IKMA2,IMKMB3,JPOL,ISMT)
               IMC_A2 = IMC*AME_G(IMKMA2,IKMB3,JPOL,ISMT)
C
               DO IR = 1,IRTOP_SPHERE
C
                  A1_ZGLA = IMC_A1*ZGLA(IR,IA2,LAMA2)
                  A2_ZFLA = IMC_A2*ZFLA(IR,IA2,LAMA2)
                  A1_JGLA = IMC_A1*JGLA(IR,IA2,LAMA2)
                  A2_JFLA = IMC_A2*JFLA(IR,IA2,LAMA2)
C
                  FZAZB(IR) = FZAZB(IR) + A1_ZGLA*ZFRB(IR,IB3,LAMB3)
     &                        - A2_ZFLA*ZGRB(IR,IB3,LAMB3)
                  FZAJB(IR) = FZAJB(IR) + A1_ZGLA*JFRB(IR,IB3,LAMB3)
     &                        - A2_ZFLA*JGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + A1_JGLA*ZFRB(IR,IB3,LAMB3)
     &                        - A2_JFLA*ZGRB(IR,IB3,LAMB3)
                  FJAJB(IR) = FJAJB(IR) + A1_JGLA*JFRB(IR,IB3,LAMB3)
     &                        - A2_JFLA*JGRB(IR,IB3,LAMB3)
C
               END DO
C
C------------------------------------------------------ ASA surface term
               IF ( .NOT.FULLPOT ) THEN
                  IR = IRTOP_SPHERE
                  PA = -0.5D0*PAB*(IMC_A1-IMC_A2)/R2DRDI(IR,IM)
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
               A1 = A_SIG_RLM(IKMA2,IMKMB3,JPOL,LMSF)
               A2 = A_SIG_RLM(IMKMA2,IKMB3,JPOL,LMSF)
C
               IF ( C_LT_EPS(A1,NON0_TOL) .AND. C_LT_EPS(A2,NON0_TOL) )
     &              CYCLE LOOP_ISF
C
               IFLAG1 = 1
C
               IMC_A1 = IMC*A1
               IMC_A2 = IMC*A2
C
               DO IR = IRTOP_SPHERE + 1,IRCRIT
C
                  IRSF = IR - IRTOP_SPHERE
C
                  A1_ZGLA = IMC_A1*ZGLA(IR,IA2,LAMA2)
                  A2_ZFLA = IMC_A2*ZFLA(IR,IA2,LAMA2)
                  A1_JGLA = IMC_A1*JGLA(IR,IA2,LAMA2)
                  A2_JFLA = IMC_A2*JFLA(IR,IA2,LAMA2)
C
                  CAUX = A1_ZGLA*ZFRB(IR,IB3,LAMB3)
     &                   - A2_ZFLA*ZGRB(IR,IB3,LAMB3)
C
                  FZAZB(IR) = FZAZB(IR) + CAUX*FLMSF(IRSF,ISF,IM)
C
                  CAUX = A1_ZGLA*JFRB(IR,IB3,LAMB3)
     &                   - A2_ZFLA*JGRB(IR,IB3,LAMB3)
C
                  FZAJB(IR) = FZAJB(IR) + CAUX*FLMSF(IRSF,ISF,IM)
C
                  CAUX = A1_JGLA*ZFRB(IR,IB3,LAMB3)
     &                   - A2_JFLA*ZGRB(IR,IB3,LAMB3)
C
                  FJAZB(IR) = FJAZB(IR) + CAUX*FLMSF(IRSF,ISF,IM)
C
                  CAUX = A1_JGLA*JFRB(IR,IB3,LAMB3)
     &                   - A2_JFLA*JGRB(IR,IB3,LAMB3)
C
                  FJAJB(IR) = FJAJB(IR) + CAUX*FLMSF(IRSF,ISF,IM)
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
         MZAZB(LAMA2,LAMB3,JPOL) = MZAZB(LAMA2,LAMB3,JPOL)
     &                             + IZAZB(IRCRIT) + SZAZB
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
C*==me_check_bxaxw.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE ME_CHECK_BXAXW(IIRR,IRUN_CHECK,ISTEP,MREF1,MREF2,IREF,
     &                          JREF,IZAZB_1,SZAZB_1,IZAZB_2,SZAZB_2,
     &                          IJAJB_1,SJAJB_1,IJAJB_2,SJAJB_2,MIRR,I,
     &                          J,JPOL)
C   ********************************************************************
C   *                                                                  *
C   *  check result of the *_BXAXW routine, as for example             *
C   *                                                                  *
C   *         MIRR_* =   < J^+(E_b) | j_lam * IZAZB | Z(E_a) >         *
C   *                  + < Z^+(E_b) | j_lam * IZAJB | Z(E_a) >         *
C   *                                                                  *
C   *                                                                  *
C   *  for  IRUN_CHECK = 1 (2) the first (second) term is checked      *
C   *  by setting   IZAZB=1, IZAJB=0        (IZAZB=0, IZAJB=1)         *
C   *  leading to  <J^+(E_b)|j_lam|Z(E_a)>  (<Z^+(E_b)|j_lam|Z(E_a)>)  *
C   *  that is equal to  MREF1              (MREF2)                    *
C   *                                                                  *
C   *  the integral functions I* may consists of 2 terms I*_1 and I*_2 *
C   *  the surface terms S* are not dealt with so far                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_ANGMOM,ONLY:NKMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL_CMP,THRESH_CMP
      PARAMETER (TOL_CMP=1D-14,THRESH_CMP=1D-14)
C
C Dummy arguments
C
      INTEGER I,IIRR,IREF,IRUN_CHECK,ISTEP,J,JPOL,JREF
      COMPLEX*16 SJAJB_1,SJAJB_2,SZAZB_1,SZAZB_2
      COMPLEX*16 IJAJB_1(NRMAX),IJAJB_2(NRMAX),IZAZB_1(NRMAX),
     &           IZAZB_2(NRMAX),MIRR(NKMMAX,NKMMAX,3,3),
     &           MREF1(NKMMAX,NKMMAX,3),MREF2(NKMMAX,NKMMAX,3)
C
C Local variables
C
      REAL*8 ADIFF,RDIFF
      INTEGER IPOL
      COMPLEX*16 MCHK,MREF
C
C*** End of declarations rewritten by SPAG
C
      IF ( ISTEP.EQ.1 ) THEN
C
         MIRR(:,:,:,:) = 0D0
C
         IF ( IRUN_CHECK.EQ.1 ) THEN
            IZAZB_1(:) = 1D0
            SZAZB_1 = 1D0
            IZAZB_2(:) = 1D0
            SZAZB_2 = 1D0
C
            IJAJB_1(:) = 0D0
            SJAJB_1 = 0D0
            IJAJB_2(:) = 0D0
            SJAJB_2 = 0D0
         ELSE
            IZAZB_1(:) = 0D0
            SZAZB_1 = 0D0
            IZAZB_2(:) = 0D0
            SZAZB_2 = 0D0
C
            IJAJB_1(:) = 1D0
            SJAJB_1 = 0D0
            IJAJB_2(:) = 1D0
            SJAJB_2 = 0D0
         END IF
C
         RETURN
C
      ELSE
C
CHECK ------------------------------------------ >   MIRR(I,J,IPOL,JPOL)
C
         DO IPOL = 1,3
C
            MCHK = MIRR(I,J,IPOL,JPOL)
C
            IF ( IRUN_CHECK.EQ.1 ) THEN
               MREF = MREF1(IREF,JREF,IPOL)
            ELSE
               MREF = MREF2(IREF,JREF,IPOL)
            END IF
C
            ADIFF = 0D0
            RDIFF = 0D0
            IF ( ABS(MCHK).GT.THRESH_CMP ) THEN
               RDIFF = ABS(1D0-MREF/MCHK)
            ELSE IF ( ABS(MREF).GT.THRESH_CMP ) THEN
               RDIFF = ABS(1D0-MCHK/MREF)
            ELSE
               ADIFF = ABS(MCHK-MREF)
            END IF
C
            IF ( (RDIFF.GT.TOL_CMP) .OR. (ADIFF.GT.THRESH_CMP) ) THEN
C
               WRITE (6,99001) IIRR,IRUN_CHECK,IREF,JREF,I,J,IPOL
               WRITE (6,99002) TOL_CMP,MCHK,MREF
               WRITE (6,99003) ADIFF,RDIFF
C
            END IF
C
         END DO
C
      END IF
C
99001 FORMAT (/,10X,'ME check for IIRR=',I2,'  IRUN=',I2,'  IREF=',I2,
     &        '  JREF=',I2,'  I=',I2,'  J=',I2,'  IPOL=',I2)
99002 FORMAT (/,10X,'MIRR and MREF deviate more than ',E10.1,/,10X,
     &        'MIRR: ',2E22.14,/,10X,'MREF: ',2E22.14)
99003 FORMAT (/,10X,'ADIFFMAX =',F20.14,/,10X,'RDIFFMAX =',F20.14,/)
C
      END
