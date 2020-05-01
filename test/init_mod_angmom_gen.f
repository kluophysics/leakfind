C*==init_mod_angmom_gen.f    processed by SPAG 6.70Rc at 21:09 on 19 Dec 2016
      SUBROUTINE INIT_MOD_ANGMOM_GEN(NL)
C   ********************************************************************
C   *                                                                  *
C   *  module to store general tables dependent ONLY on     NLMAX      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,L_LMS,M_LMS,A_SIGMA,IKM_KAP_MUEM05,
     &    KAPPA_IKM,MUEM05_IKM,L_IKM,NLMAX,NSPIN,ILMBOT_LM,ILMTOP_LM,
     &    IKMLLIM1,IKMLLIM2,IMKM_IKM,FACT,DBLFACT_L,NLMMADMAX,JP05_IKM,
     &    LB_IKM,WKM1,WKM2,WKM3,WKM4,WLM1,WLM2,WLM3,WKM_LPK,IPIVKM,
     &    NKMMAX,NLMMAX,CGC,CREL,RC,RREL,DELTA_L_EXT,NL_EXT,NK_EXT,
     &    NLM_EXT,NKM_EXT,NKMP_EXT,NL_AME_RLM_EXT,NLM_AME_RLM_EXT
      USE MOD_CONSTANTS,ONLY:SQRT_2
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_ANGMOM_GEN')
C
C Dummy arguments
C
      INTEGER NL
C
C Local variables
C
      INTEGER I,IA_ERR,IFLAG,IKM,ILJM,IMKM,JP05,K,KAPPA,L,LB,LLIM,LM,
     &        LMLIM,M,MUEM05,NLM
      REAL*8 XJ,XJM,XJP
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------------------ FACT
C
      ALLOCATE (FACT(0:100))
C
      FACT(0) = 1.0D0
      DO I = 1,100
         FACT(I) = FACT(I-1)*DBLE(I)
      END DO
C
C------------------------------------------------------------- DBLFACT_L
C
      ALLOCATE (DBLFACT_L(0:(8*NLMAX)))
C
      DBLFACT_L(0) = 1.0D0
      DO L = 1,8*NLMAX
         DBLFACT_L(L) = DBLFACT_L(L-1)*DBLE(2*L+1)
      END DO
C
C----------------------------------------------------------------- NSPIN
      IF ( IREL.EQ.2 ) THEN
         NSPIN = 2
      ELSE
         NSPIN = 1
      END IF
C ----------------------------------------------------------------------
C
C---------------------- dimension of charge multipoles - ASA and FULLPOT
C
      NLMMADMAX = (2*(NLMAX-1)+1)**2
C
C----------------------------------------------- tables  l(LM) and m(LM)
C   with linear index  LM = L*L + L + M + 1
C   allow for highest LM to occur for shape functions: l(SF) = 4 * l_max
C
      LLIM = 4*(NL-1)
      LMLIM = (LLIM+1)**2
C
      ALLOCATE (L_LM(LMLIM),M_LM(LMLIM))
      ALLOCATE (ILMBOT_LM(LMLIM),ILMTOP_LM(LMLIM))
C
      LM = 0
      DO L = 0,LLIM
         DO M = -L, + L
            LM = LM + 1
            L_LM(LM) = L
            M_LM(LM) = M
            ILMBOT_LM(LM) = L**2 + 1
            ILMTOP_LM(LM) = L**2 + 2*L + 1
         END DO
      END DO
C
C--------------------------------------------- tables  l(LMS) and m(LMS)
C   with linear index  LMS = NLM*(IS-1) + L*L + L + M + 1
C   LMS <= 2*NLM = 2*NL**2
C
      LLIM = NL - 1
      NLM = NL**2
C
      ALLOCATE (L_LMS(2*NLM),M_LMS(2*NLM))
C
      LM = 0
      DO L = 0,LLIM
         DO M = -L, + L
            LM = LM + 1
            L_LMS(LM) = L
            M_LMS(LM) = M
            L_LMS(NLM+LM) = L
            M_LMS(NLM+LM) = M
         END DO
      END DO
C-----------------------------------------------------------------------
C                   set extended l-dependent parameters
C-----------------------------------------------------------------------
C
      IF ( DELTA_L_EXT.LT.0 ) CALL STOP_MESSAGE(ROUTINE,'DELTA_L_EXT<0')
      NL_EXT = NL + DELTA_L_EXT
      NLM_EXT = NL_EXT*NL_EXT
      NK_EXT = 2*NL_EXT - 1
      NKM_EXT = 2*NLM_EXT
      NKMP_EXT = NKM_EXT + 2*NL_EXT
C
C----------------------------------------- array size for <LAM|RLM|LAM'>
      NL_AME_RLM_EXT = 2*(NL_EXT-1) + 1
      NLM_AME_RLM_EXT = NL_AME_RLM_EXT**2
C
C---------------------------------------------------- relativistic table
C   set up the tables for IKM = 1,...,IKMPMAX
C   with  IKM = 2*L*(J+1/2) + J    + MUE    + 1
C             = 2*L* JP05   + JP05 + MUEM05 + 1
C   IKMPMAX corresponds to KAPPA = l_max + 1, i.e. f_5/2 for l_max=2
C   to deal properly with the minor components for l_max
C   IKM_KAP_MUEM05:  linear index IKM for (KAPPA,MUE)
C   KAPPA_IKM:    KAPPA(IKM)
C   MUEM05_IKM:   MUEM05(IKM)
C   L_IKM:        L(IKM)
C
      ALLOCATE (IKM_KAP_MUEM05(-NL_EXT:NL_EXT,-NL_EXT:NL_EXT))
      ALLOCATE (L_IKM(NKMP_EXT),LB_IKM(NKMP_EXT),JP05_IKM(NKMP_EXT))
      ALLOCATE (KAPPA_IKM(NKMP_EXT),MUEM05_IKM(NKMP_EXT))
      ALLOCATE (IMKM_IKM(NKMP_EXT))
      ALLOCATE (IKMLLIM1(NKMP_EXT),IKMLLIM2(NKMP_EXT))
C
      IKM_KAP_MUEM05(:,:) = 0
      KAPPA_IKM(:) = 0
      MUEM05_IKM(:) = 0
C
      IFLAG = 0
      IKM = 0
      I = 0
      DO K = 1,(NK_EXT+1)
         L = K/2
         IF ( K.EQ.2*L ) THEN
            KAPPA = L
         ELSE
            KAPPA = -L - 1
         END IF
         LB = L - SIGN(1,KAPPA)
         JP05 = ABS(KAPPA)
         IF ( JP05.GT.NL_EXT ) IFLAG = 1
C
         DO MUEM05 = -JP05,JP05 - 1
            IKM = IKM + 1
            IF ( IKM.LE.NKMP_EXT ) THEN
               L_IKM(IKM) = L
               LB_IKM(IKM) = L - SIGN(1,KAPPA)
               JP05_IKM(IKM) = JP05
               KAPPA_IKM(IKM) = KAPPA
               MUEM05_IKM(IKM) = MUEM05
               IKM_KAP_MUEM05(KAPPA,MUEM05) = IKM
               IMKM = LB*2*JP05 + JP05 + MUEM05 + 1
               IMKM_IKM(IKM) = IMKM
            END IF
            ILJM = 2*L*JP05 + JP05 + MUEM05 + 1
            IF ( IKM.NE.ILJM ) IFLAG = 1
         END DO
C
C ----------------------------------------------------------------------
C     find the bounding indices  IKMLLIM1  and  IKMLLIM2  for IKM-loops
C     assuming that the matrix elements are diagonal with respect to l
C     this does not hold for B_hf for which there are l - (l+/-2)-terms
C ----------------------------------------------------------------------
C
         XJM = L - 0.5D0
         XJP = L + 0.5D0
         IF ( MOD(K,2).EQ.1 ) THEN
            XJ = L + 0.5D0
         ELSE
            XJ = L - 0.5D0
         END IF
         DO M = 1,NINT(2*XJ) + 1
            I = I + 1
            IKMLLIM1(I) = NINT(L*2*(XJM+0.5D0)+1)
            IKMLLIM2(I) = NINT(L*2*(XJP+0.5D0)+2*XJP+1)
         END DO
C
      END DO
C
      IF ( L.NE.NL_EXT ) IFLAG = 1
      IF ( IKM.NE.NKMP_EXT ) IFLAG = 1
C
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99001) L,NL_EXT,IKM,NKMP_EXT,JP05,NL_EXT
         CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
      END IF
C
C------------ calculate the matrix elements for the Pauli spin operators
C
C   A_SIGMA(IS,LMABDA,IS') = <chi(m_s)|SIG(lambda)|chi(m_s')>
C   IS = 1,2 == m_s = -1/2, +1/2   lambda = -1, 0, +1
C
      ALLOCATE (A_SIGMA(2,-1:+1,2))
C
      A_SIGMA(1:2,-1:+1,1:2) = 0D0
C
      A_SIGMA(1,-1,2) = SQRT_2
      A_SIGMA(1,0,1) = -1D0
C
      A_SIGMA(2,0,2) = +1D0
      A_SIGMA(2,+1,1) = SQRT_2
C
C ----------------------------------- workspace for matrix inversion etc
C
      ALLOCATE (IPIVKM(NKMMAX),WKM_LPK(NKMMAX,NKMMAX))
      ALLOCATE (WKM1(NKMMAX,NKMMAX),WLM1(NLMMAX,NLMMAX))
      ALLOCATE (WKM2(NKMMAX,NKMMAX),WLM2(NLMMAX,NLMMAX))
      ALLOCATE (WKM3(NKMMAX,NKMMAX),WLM3(NLMMAX,NLMMAX))
      ALLOCATE (WKM4(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WKM1')
C
C---------------------------------------- variables depending only on NL
C
      ALLOCATE (CREL(NKMMAX,NKMMAX),CGC(NKMP_EXT,2))
      ALLOCATE (RC(NKMMAX,NKMMAX),RREL(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CREL')
C
      CALL CALCCGC(CGC,NK_EXT,NKMP_EXT)
C
C------------------------------ supply matrices for basis transformation
C
      CALL BASTRMAT(NL-1,RC,CREL,RREL,NKMMAX)
C
99001 FORMAT (//,'ERROR in <INIT_MOD_ANGMOM_GEN> ',/,10X,'L =',I3,
     &        '   NL_EXT   =',I3,/,10X,'IKM =',I3,'   NKMP_EXT =',I3,/,
     &        10X,'JP05=',I3,'   NL_EXT   =',I3)
C
      END
