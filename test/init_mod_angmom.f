C*==init_mod_angmom.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE INIT_MOD_ANGMOM(IREL,ISMQHFI,IPRINT,NLQ0,NKKR,NLFPMAX)
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tables dependent ONLY on                    *
C   *                                                                  *
C   *      NLMAX  and  NQMAX   and derived array dimensions            *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CALCMODE,ONLY:ORBPOL,BLCOUPL,BREITINT,TASK
      USE MOD_ANGMOM,ONLY:IMKM_IKM,IKM1LIN,IKM2LIN,NSOLLM,KAPTAB,LTAB,
     &    LBTAB,NMUETAB,SDIA,SMDIA,SOFF,SMOFF,QDIA,QMDIA,QOFF,QMOFF,NKQ,
     &    NLINQ,NKMQ,IND0Q,NLQ,TXT_L,NCPLWF,NCPLWF_RA,NCPLWF_LA,
     &    NCPLWF_RB,NCPLWF_LB,CGCFP,AMEOPC,AMEOPO,NMVECMAX,NL,NK,NLM,
     &    NKM,NLIN,NKMMAX,NMUEMAX,LINMAX,NLMAX,NLMQ,AMEBI1,AMEBI2,NXM,
     &    NXM_Q,TXT_J,IKMCPL_KM,NKMCPL_KM,NLMMAX,LMAT2,LMAT3,A_RGNT,
     &    AG_RGNT,AF_RGNT,NL_EXT,NLM_EXT,NKM_EXT,TAUT,TSST,MSST,SSST,
     &    TAUQ,TSSQ,MSSQ,MEZJ,MEZZ,NMEMAX,WXM1,WXM2,WXM3,WXM4,NXMMAX,
     &    AME_RLM,AME_G,NL_AME_RLM_EXT,NLM_AME_RLM_EXT,U_CS,U_SC,ISMT,
     &    IOMT
      USE MOD_SITES,ONLY:NQ,NQMAX
      USE MOD_TYPES,ONLY:NTMAX,NLMFP
      USE MOD_CONSTANTS,ONLY:C0,C1,CI,SQRT_2,SQRT_4PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_ANGMOM')
C
C Dummy arguments
C
      INTEGER IPRINT,IREL,ISMQHFI,NKKR,NLFPMAX
      INTEGER NLQ0(NQMAX)
C
C Local variables
C
      COMPLEX*16 AII,AIJ,AJIX,A_F,A_G,RRELX(:,:),UTST(3,3),WKMX1(:,:),
     &           WKMX2(:,:)
      LOGICAL C_GT_EPS
      REAL*8 DEL_A,FMOFF,FOFF,GMOFF,GNTTAB(:,:,:),GOFF,M1PWL,W
      REAL*8 GAUNT_RYLM
      INTEGER I,IA_ERR,IFLAG,IKM1,IKM2,IL,IMKM1,IMKM2,IMUE,IMVEC,IPOL,
     &        IQ,J,J2,JCPLP05,JP05,K,KAP,KAPCPL,L,L1,L2,L3,LIN,LM,LM1,
     &        LM2,LM3,M,M1,M2,M3,MMAX,MUE2,MUEM05,N
      INTEGER IKAPMUE
      CHARACTER*25 STR25_A,STR25_G
      CHARACTER*5 STR5
      CHARACTER*1 TXT_LAUX(0:3)
C
C*** End of declarations rewritten by SPAG
C
C--------------------------------- temporary table of Gaunt coefficients
C
      ALLOCATABLE GNTTAB,WKMX1,WKMX2,RRELX
C
      DATA TXT_LAUX/'s','p','d','f'/
      DATA STR25_G/'GNT(RLM) for lm = '/,STR25_A/'GNT(rel) for lm = '/
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( IPRINT.GE.5 ) WRITE (6,99004) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                          'NQMAX',NQMAX
C
C-----------------------------------------------------------------------
C  create transformation matri   U  cartesian/spherical ccordinates
C-----------------------------------------------------------------------
C
C  ( a_x )        ( a_- )                   1    (  1     0     -1   )
C  ( a_y ) = U_CS ( a_0 )   with  U_CS = ------- (  i     0      i   )
C  ( a_z )        ( a_+ )                sqrt(2) (  0  sqrt(2)   0   )
C
C-- convert polarisation: SPHERICAL (-),(0),(+) to CARTESIAN (x),(y),(z)
C
      W = 1.0D0/SQRT_2
C
      U_SC(1,1) = W
      U_SC(1,2) = -CI*W
      U_SC(1,3) = 0.0D0
      U_SC(2,1) = 0.0D0
      U_SC(2,2) = 0.0D0
      U_SC(2,3) = 1.0D0
      U_SC(3,1) = -W
      U_SC(3,2) = -CI*W
      U_SC(3,3) = 0.0D0
C
      U_CS(:,:) = TRANSPOSE(U_SC(:,:))
      U_CS(:,:) = DCONJG(U_CS(:,:))
C
      UTST = MATMUL(U_CS,U_SC)
C
      IFLAG = 0
      DO I = 1,3
         IF ( ABS(UTST(I,I)-1D0).GT.1D-12 ) IFLAG = 1
         UTST(I,I) = 0D0
         DO J = 1,3
            IF ( ABS(UTST(I,J)).GT.1D-12 ) IFLAG = 1
         END DO
      END DO
C
      IF ( IPRINT.GE.5 .OR. IFLAG.NE.0 ) THEN
         WRITE (6,*) ' '
         WRITE (6,*) ' transformation matrices U_CS and U_SC'
         WRITE (6,*) ' U_CS'
         WRITE (6,'(6e15.6)') ((U_CS(I,J),J=1,3),I=1,3)
         WRITE (6,*) ' '
         WRITE (6,*) ' U_SC'
         WRITE (6,'(6e15.6)') ((U_SC(I,J),J=1,3),I=1,3)
         WRITE (6,*) ' '
         WRITE (6,*) ' U_CS*U_SC'
         WRITE (6,'(6e15.6)') UTST
         IF ( IFLAG.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'U_CS*U_SC <> 1')
      END IF
C
C---------------------------------------- variables depending only on NQ
C
      ALLOCATE (NKQ(NQMAX),NLINQ(NQMAX),NKMQ(NQMAX),IND0Q(NQMAX))
      ALLOCATE (NLQ(NQMAX),NLMQ(NQMAX),NXM_Q(NQMAX))
C
C---------------------------------------- variables depending only on NL
C
      M = NKMMAX
      ALLOCATE (NCPLWF(NKMMAX))
      ALLOCATE (NCPLWF_RA(NKMMAX),NCPLWF_LA(NKMMAX))
      ALLOCATE (NCPLWF_RB(NKMMAX),NCPLWF_LB(NKMMAX))
      ALLOCATE (QDIA(M),QOFF(M),QMDIA(M),QMOFF(M))
      ALLOCATE (SDIA(M),SOFF(M),SMDIA(M),SMOFF(M))
      ALLOCATE (IKMCPL_KM(M),NKMCPL_KM(M))
      ALLOCATE (IKM1LIN(LINMAX),IKM2LIN(LINMAX))
      ALLOCATE (KAPTAB(NMUEMAX),LBTAB(NMUEMAX))
      ALLOCATE (LTAB(NMUEMAX),NMUETAB(NMUEMAX))
      ALLOCATE (NSOLLM(NLMAX,NMUEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc NSOLLM')
      QDIA = 999999D0
      QMDIA = 999999D0
      QOFF = 999999D0
      NCPLWF = 999999
C
C----------------------------------------------------------------- TXT_L
C
      ALLOCATE (TXT_L(0:NLMAX))
C
      DO L = 0,MIN(3,NLMAX)
         TXT_L(L) = TXT_LAUX(L)
      END DO
      DO L = 4,NLMAX
         TXT_L(L) = CHAR(ICHAR('f')+L-4)
      END DO
C
C----------------------------------------------------------------- TXT_J
C
      M = MAX(NLMAX+1,5)
      ALLOCATE (TXT_J(M))
C
      DO JP05 = 1,M
         I = 2*JP05 - 1
         WRITE (TXT_J(JP05),FMT='(I2,''/2'')') I
      END DO
C
C=======================================================================
C
      DO IQ = 1,NQ
         NLQ(IQ) = NLQ0(IQ)
      END DO
C
      NLM = NL**2
C
      IF ( IREL.LT.2 ) THEN
C
C-------------------------------------- non-relativistic KKR-calculation
         NK = NL
         NKM = NLM
C
         NKKR = 0
         DO IQ = 1,NQ
            IF ( IQ.EQ.1 ) THEN
               IND0Q(IQ) = 0
            ELSE
               IND0Q(IQ) = IND0Q(IQ-1) + NLQ(IQ-1)**2
            END IF
            NKKR = NKKR + NLQ(IQ)**2
            NKQ(IQ) = NLQ(IQ)
            NKMQ(IQ) = NLQ(IQ)**2
            NLMQ(IQ) = NLQ(IQ)**2
            NLINQ(IQ) = NLQ(IQ)**2
C
            LIN = 0
            DO IL = 1,NLMAX
               LTAB(IL) = IL - 1
               NMUETAB(IL) = 2*LTAB(IL) + 1
C
               DO IMUE = 1,NMUETAB(IL)
                  NSOLLM(IL,IMUE) = 1
                  LIN = LIN + 1
                  IKM1LIN(LIN) = LIN
                  IKM2LIN(LIN) = LIN
               END DO
            END DO
         END DO
C
         NLIN = NKKR
C
      ELSE
C------------------------------------------------ scalar-relativistic OR
C------------------------------------------ relativistic KKR-calculation
C                                                              IREL >= 2
         NK = 2*NL - 1
         NKM = 2*NLM
C
C                     s   p   p   d   d   f   f   g   g   h
C     DATA LTAB    /  0,  1,  1,  2,  2,  3,  3,  4,  4,  5 /
C     DATA LBTAB   /  1,  0,  2,  1,  3,  2,  4,  3,  5,  4 /
C     DATA KAPTAB  / -1,  1, -2,  2, -3,  3, -4,  4, -5,  5 /
C     DATA NMUETAB /  2,  2,  4,  4,  6,  6,  8,  8, 10, 10 /
C
         DO I = 1,NMUEMAX
            LTAB(I) = I/2
            IF ( 2*LTAB(I).EQ.I ) THEN
               LBTAB(I) = I
               LBTAB(I) = LTAB(I) - 1
               KAPTAB(I) = LTAB(I)
            ELSE
               LBTAB(I) = LTAB(I) + 1
               KAPTAB(I) = -LTAB(I) - 1
            END IF
            NMUETAB(I) = 2*ABS(KAPTAB(I))
         END DO
C
         DO IL = 1,NLMAX
            MMAX = 2*IL
            DO IMUE = 1,MMAX
               IF ( (IMUE.EQ.1) .OR. (IMUE.EQ.MMAX) ) THEN
                  NSOLLM(IL,IMUE) = 1
               ELSE
                  NSOLLM(IL,IMUE) = 2
               END IF
            END DO
         END DO
C
         I = 0
         DO K = 1,NK
            L = K/2
            IF ( MOD(K,2).EQ.1 ) THEN
               KAP = -L - 1
            ELSE
               KAP = L
            END IF
            JP05 = ABS(KAP)
            KAPCPL = -KAP - 1
            JCPLP05 = ABS(KAPCPL)
            DO MUEM05 = -JP05,JP05 - 1
               I = I + 1
               IF ( MUEM05.EQ.-JP05 .OR. MUEM05.EQ.(JP05-1) ) THEN
                  NKMCPL_KM(I) = 1
                  IKMCPL_KM(I) = 0
               ELSE
                  NKMCPL_KM(I) = 2
                  IKMCPL_KM(I) = L*2*JCPLP05 + JCPLP05 + MUEM05 + 1
               END IF
            END DO
         END DO
C
         NLIN = 2*NL*(2*NL-1)
C
         CALL IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,LINMAX,
     &               NLMAX)
C
         IF ( ISMQHFI.NE.0 ) CALL CALCSMQ(SDIA,SMDIA,SOFF,SMOFF,QDIA,
     &        QMDIA,QOFF,QMOFF)
C
         CALL CALCAME
C
         IF ( IPRINT.GE.2 ) WRITE (6,99002)
         I = 0
         NKKR = 0
         DO IQ = 1,NQ
            NKQ(IQ) = 2*NLQ(IQ) - 1
            NKMQ(IQ) = 2*NLQ(IQ)**2
            NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
            NLMQ(IQ) = NLQ(IQ)**2
C
            NKKR = NKKR + 2*NLQ(IQ)**2
            IF ( IQ.EQ.1 ) THEN
               IND0Q(IQ) = 0
            ELSE
               IND0Q(IQ) = IND0Q(IQ-1) + 2*NLQ(IQ-1)**2
            END IF
C
            DO K = 1,NKQ(IQ)
               L = LTAB(K)
               KAP = KAPTAB(K)
               J2 = 2*ABS(KAP) - 1
               MUE2 = -J2 - 2
C
               DO M = 1,NMUETAB(K)
                  IF ( IPRINT.GE.2 ) THEN
                     I = I + 1
                     MUE2 = MUE2 + 2
                     MUEM05 = (MUE2-1)/2
                     IKM1 = IKAPMUE(KAP,MUEM05)
                     IF ( KAP.NE.-1 ) THEN
                        IKM2 = IKAPMUE(-KAP-1,MUEM05)
                     ELSE
                        IKM2 = IKM1
                     END IF
                     IMKM1 = IMKM_IKM(IKM1)
                     IMKM2 = IMKM_IKM(IKM2)
C
                     IF ( IQ.EQ.1 ) THEN
                        IF ( ABS(MUE2).GT.2*L ) THEN
                           GOFF = 0D0
                           GMOFF = 0D0
                           FOFF = 0D0
                           FMOFF = 0D0
                        ELSE
                           GOFF = AME_G(IKM1,IKM2,2,ISMT)
                           GMOFF = AME_G(IMKM1,IMKM2,2,ISMT)
                           FOFF = AME_G(IKM1,IKM2,2,IOMT)
                           FMOFF = AME_G(IMKM1,IMKM2,2,IOMT)
                        END IF
                        WRITE (6,99003) I,L,KAP,J2,MUE2,I,
     &                                  AME_G(IKM1,IKM1,2,ISMT),GOFF,
     &                                  AME_G(IMKM1,IMKM1,2,ISMT),GMOFF,
     &                                  AME_G(IKM1,IKM1,2,IOMT),FOFF,
     &                                  AME_G(IMKM1,IMKM1,2,IOMT),FMOFF
                     ELSE
                        WRITE (6,99003) I,L,KAP,J2,MUE2,I
                     END IF
                  END IF
               END DO
            END DO
         END DO
C
         IF ( IREL.EQ.2 ) THEN
            NKKR = NKKR/2
            IND0Q(1:NQ) = IND0Q(1:NQ)/2
         END IF
C
      END IF
C
C------------------------------------------------------------------- NXM
C
      IF ( IREL.LE.1 ) THEN
         NXM = NLM
         NXM_Q(1:NQ) = NLMQ(1:NQ)
C ----------------------------------------------------------------------
      ELSE IF ( IREL.EQ.2 ) THEN
         NXM = NLM
         NXM_Q(1:NQ) = NLMQ(1:NQ)
C ----------------------------------------------------------------------
      ELSE IF ( IREL.EQ.3 ) THEN
         NXM = NKM
         NXM_Q(1:NQ) = NKMQ(1:NQ)
C ----------------------------------------------------------------------
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'IREL out of bounds')
      END IF
C
      NXMMAX = NXM
C
C ----------------------------------- workspace for matrix inversion etc
C
      ALLOCATE (WXM1(NXMMAX,NXMMAX),WXM2(NXMMAX,NXMMAX))
      ALLOCATE (WXM3(NXMMAX,NXMMAX),WXM4(NXMMAX,NXMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WXM1')
C
C------------------------------------ variables for FULLPOT calculations
C
      ALLOCATE (CGCFP(2*NLFPMAX**2+2*NLFPMAX,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc CGCFP')
      CGCFP(:,:) = 999999D0
C
      CALL CALCCGC(CGCFP,2*NLFPMAX-1,2*NLFPMAX**2+2*NLFPMAX)
C
C ------------------------- angular matrix elements for Brooks formalism
C
      IF ( IREL.EQ.3 .AND. (ORBPOL(1:6).EQ.'BROOKS' .OR. BLCOUPL .OR. 
     &     TASK(1:3).EQ.'CHI') ) THEN
C
         CALL AMEOP
C
      ELSE
C
         ALLOCATE (AMEOPC(1,1,1),AMEOPO(1,1,1))
C
      END IF
C
C ------------------------ angular matrix elements for Breit interaction
C                   to suppress retardation term: set GBIL and GBIG to 0
C
      IF ( BREITINT ) THEN
C
         CALL AMEBI
C
      ELSE
C
         ALLOCATE (AMEBI1(1,1,1,1),AMEBI2(1,1,1,1))
C
      END IF
C
C=======================================================================
C   it will be assumed that  IREL< 2:  NLM = NLMMAX = NKM = NKMMAX / 2
C                            IREL>=2:  NLM = NLMMAX   NKM = NKMMAX
C                                            NLMMAX =       NKMMAX / 2
C
      IF ( IREL.LT.2 ) THEN
         IF ( NLM.NE.NLMMAX .OR. NLM.NE.NKM .OR. NLM.NE.NKMMAX/2 )
     &        CALL STOP_MESSAGE(ROUTINE,'error occurred for IREL<2')
      ELSE IF ( NLM.NE.NLMMAX .OR. NLM.NE.NKM/2 .OR. NKM.NE.NKMMAX )
     &          THEN
         CALL STOP_MESSAGE(ROUTINE,'error occurred for IREL>=2')
      END IF
C
C=======================================================================
C           set up the L-matrix   L[i,j] = delta[i,j]*(-1)^l
C           LMAT2:  (l,m,s)-rep.  LMAT3:  (kappa,mue)-rep.
C=======================================================================
C
      ALLOCATE (LMAT2(NKMMAX,NKMMAX),LMAT3(NKMMAX,NKMMAX))
C
      LMAT2(:,:) = C0
      I = 0
      M1PWL = -1D0
      DO IL = 1,NLMAX
         M1PWL = -M1PWL
         DO M = 1,(2*IL-1)
            I = I + 1
            LMAT2(I,I) = M1PWL
            LMAT2(NLM+I,NLM+I) = M1PWL
         END DO
      END DO
C
      LMAT3(:,:) = C0
      I = 0
      M1PWL = -1D0
      DO IL = 1,NLMAX
         M1PWL = -M1PWL
         DO M = 1,2*(2*IL-1)
            I = I + 1
            LMAT3(I,I) = M1PWL
         END DO
      END DO
C
C=======================================================================
C  matrix of Gaunt coefficients  A_RGNT  A^L_{LAM,LAM'} = <LAM|L|LAM'>
C=======================================================================
C
C-----------------------------------------------------------------------
C             Gaunt coefficients for REAL spherical harmonics
C-----------------------------------------------------------------------
C
      ALLOCATE (GNTTAB(NLM_EXT,NLM_EXT,NLM_AME_RLM_EXT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc: GNTTAB')
      GNTTAB(:,:,:) = 0D0
C
      LM1 = 0
      DO L1 = 0,(NL_EXT-1)
         DO M1 = -L1, + L1
            LM1 = LM1 + 1
C
            LM2 = 0
            DO L2 = 0,(NL_EXT-1)
               DO M2 = -L2, + L2
                  LM2 = LM2 + 1
C
                  LM3 = 0
                  DO L3 = 0,(NL_AME_RLM_EXT-1)
                     DO M3 = -L3, + L3
                        LM3 = LM3 + 1
C
                        GNTTAB(LM1,LM2,LM3)
     &                     = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C
                     END DO
                  END DO
C
               END DO
            END DO
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      IF ( IREL.LE.2 ) THEN
C
         ALLOCATE (A_RGNT(NKMMAX,NKMMAX,NLMFP))
C
         A_RGNT(:,:,:) = C0
C
         DO LM = 1,NLMFP
            DO LM1 = 1,NLM
               DO LM2 = 1,NLM
                  A_RGNT(LM1,LM2,LM) = DCMPLX(GNTTAB(LM1,LM2,LM),0D0)
               END DO
            END DO
         END DO
C
      ELSE
C
         N = NKM_EXT
         M = NKM_EXT
C
         ALLOCATE (WKMX1(NKM_EXT,NKM_EXT),WKMX2(NKM_EXT,NKM_EXT))
         ALLOCATE (RRELX(NKM_EXT,NKM_EXT))
C
C------------------------------ supply matrices for basis transformation
C
         CALL BASTRMAT(NL_EXT-1,WKMX1,WKMX2,RRELX,NKM_EXT)
C
         ALLOCATE (AG_RGNT(NKM_EXT,NKM_EXT,NLM_AME_RLM_EXT))
         ALLOCATE (AF_RGNT(NKM_EXT,NKM_EXT,NLM_AME_RLM_EXT))
C
         AG_RGNT(:,:,:) = C0
         AF_RGNT(:,:,:) = C0
C
         IFLAG = 0
C
         DO LM = 1,NLM_AME_RLM_EXT
C
            WKMX1(:,:) = C0
            DO LM1 = 1,NLM_EXT
               DO LM2 = 1,NLM_EXT
                  WKMX1(LM1,LM2) = DCMPLX(GNTTAB(LM1,LM2,LM),0D0)
                  WKMX1(NLM_EXT+LM1,NLM_EXT+LM2) = WKMX1(LM1,LM2)
               END DO
            END DO
C
            CALL ZGEMM('C','N',N,N,N,C1,RRELX,M,WKMX1,M,C0,WKMX2,M)
            CALL ZGEMM('N','N',N,N,N,C1,WKMX2,M,RRELX,M,C0,
     &                 AG_RGNT(1,1,LM),M)
C
C           CALL CHANGEREP(NKM,NKMMAX,WKM1,'RLM>REL',AG_RGNT(1,1,LM))
C
C-----------------------------------------------------------------------
C     set up the angular matrix element for the small component
C-----------------------------------------------------------------------
C
            DO IKM1 = 1,NKM_EXT
               IMKM1 = IMKM_IKM(IKM1)
               DO IKM2 = 1,NKM_EXT
                  IMKM2 = IMKM_IKM(IKM2)
C
                  IF ( IMKM1.LE.NKM_EXT .AND. IMKM2.LE.NKM_EXT )
     &                 AF_RGNT(IKM1,IKM2,LM) = AG_RGNT(IMKM1,IMKM2,LM)
C
               END DO
            END DO
C
C----------------------------- print Gaunt in LMS and IKM representation
            IF ( IPRINT.GE.2 ) THEN
               WRITE (STR5,'(I5)') LM
               CALL STRING_TRIM_LEFT(STR5)
               STR25_G = STR25_G(1:18)//STR5
               CALL CMATSTRUCT(STR25_G,WKMX1,NKM,NKM_EXT,2,2,0,1D-8,6)
               STR25_A = STR25_A(1:18)//STR5
               CALL CMATSTRUCT(STR25_A,AG_RGNT(1,1,LM),NKM,NKM_EXT,3,3,
     &                         0,1D-8,6)
            END IF
C
C---------------------------- check that  <LAM|L|LAM'> =  <-LAM|L|-LAM'>
C
            DO IKM1 = 1,NKM_EXT
               IMKM1 = IMKM_IKM(IKM1)
               DO IKM2 = 1,NKM_EXT
                  IMKM2 = IMKM_IKM(IKM2)
C
                  IF ( IMKM1.LE.NKM_EXT .AND. IMKM2.LE.NKM_EXT ) THEN
C
                     A_G = AG_RGNT(IKM1,IKM2,LM)
                     A_F = AG_RGNT(IMKM1,IMKM2,LM)
                     DEL_A = ABS(A_G-A_F)
C
                     IF ( DEL_A.GT.1D-12 ) THEN
                        WRITE (6,99001) IKM1,IKM2,LM,A_G,IMKM1,IMKM2,LM,
     &                                  A_F
                        IFLAG = 1
                     END IF
C
                  END IF
               END DO
            END DO
C
         END DO
C
         IF ( IFLAG.EQ.1 ) CALL STOP_MESSAGE(ROUTINE,'A(g) !=A(f)')
C
C-----------------------------------------------------------------------
C                all purpose angular matrix elements
C-----------------------------------------------------------------------
C
         ALLOCATE (AME_RLM(NKM_EXT,NKM_EXT,3,NLM_AME_RLM_EXT,NMVECMAX))
         AME_RLM(:,:,:,:,:) = 0D0
C
         N = NKM_EXT
         DO IMVEC = 1,NMVECMAX
            DO IPOL = 1,3
               DO LM = 1,NLM_AME_RLM_EXT
C
                  AME_RLM(1:N,1:N,IPOL,LM,IMVEC)
     &               = MATMUL(AG_RGNT(1:N,1:N,LM),
     &               AME_G(1:N,1:N,IPOL,IMVEC))
C
C-----------------------------------------------------------------------
C     check and enforce hermiticity of the angular matrix elements
C-----------------------------------------------------------------------
C
                  IF ( IPOL.NE.2 ) CYCLE
                  IF ( IMVEC.GT.2 ) CYCLE
C
                  DO J = 1,N
                     I = J
                     AII = AME_RLM(I,J,IPOL,LM,IMVEC)
                     IF ( DIMAG(AII).GT.1D-14 ) THEN
                        WRITE (6,99005) IMVEC,IPOL,LM,I,J,AII
                        CALL STOP_MESSAGE(ROUTINE,'AME_RLM not REAL')
                     END IF
                     AME_RLM(I,J,IPOL,LM,IMVEC) = DREAL(AII)
                     DO I = J + 1,N
                        AIJ = AME_RLM(I,J,IPOL,LM,IMVEC)
                        AJIX = DCONJG(AME_RLM(J,I,IPOL,LM,IMVEC))
                        IF ( C_GT_EPS((AIJ-AJIX),1D-14) ) THEN
                           WRITE (6,99005) IMVEC,IPOL,LM,I,J,AIJ,AJIX
                           CALL STOP_MESSAGE(ROUTINE,
     &                        'AME_RLM not HERMITEAN')
                        END IF
                        AIJ = (AIJ+AJIX)/2
                        AME_RLM(I,J,IPOL,LM,IMVEC) = AIJ
                        AME_RLM(J,I,IPOL,LM,IMVEC) = DCONJG(AIJ)
                     END DO
                  END DO
C-----------------------------------------------------------------------
C
               END DO
            END DO
         END DO
C
      END IF
C
      DEALLOCATE (GNTTAB)
C
C=======================================================================
C        type and site dependent scattering and overlap matrices
C=======================================================================
C
      M = NKMMAX
      ALLOCATE (MEZJ(M,M,NTMAX,NMEMAX),MEZZ(M,M,NTMAX,NMEMAX))
      ALLOCATE (TSST(M,M,NTMAX),MSST(M,M,NTMAX),TAUT(M,M,NTMAX))
      ALLOCATE (TSSQ(M,M,NQMAX),MSSQ(M,M,NQMAX),TAUQ(M,M,NQMAX))
C
      ALLOCATE (SSST(M,M,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEZZ')
C
C=======================================================================
99001 FORMAT (/,3I4,'  A(g) ',2F20.12,/,3I4,'  A(f) ',2F20.12)
99002 FORMAT (' ',/,10X,'RELATIVISTIC QUANTUM NUMBERS',/,10X,
     &        '============================',//,10X,
     &        ' I  L  KAP  J   MUE    I','  GDIA   GOFF   GMDIA  GMOFF',
     &        '  FDIA   FOFF   FMDIA  FMOFF ')
99003 FORMAT (8X,I4,I3,I4,I3,'/2',I3,'/2',I5,8F7.3)
99004 FORMAT ('<',A,'>',('  ',A,'=',I4))
99005 FORMAT (/,10X,'setting up angular matrix elements AME_RLM',/,10X,
     &        'for IMVEC =',I3,'  IPOL =',I3,'  LM =',I3,'  I =',I3,
     &        '  J =',I3,/,10X,'AIJ = ',2E25.15,/,:,10X,'AJI = ',
     &        2E25.15,/)
      END
