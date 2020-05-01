C*==pshift.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE PSHIFT(IPRINT,TSST,MSST,SSST,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of the  phase shift for checking purposes      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NETAB
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,DRDI
      USE MOD_ANGMOM,ONLY:TXT_L,NMEMAX,NKMMAX,NLMAX,NKM,NK,NL,NLQ
      USE MOD_TYPES,ONLY:CTL,NTMAX,IMT,Z,BT,VT,NT,LTXT_T,TXT_T
      USE MOD_FILES,ONLY:IFILCBWF,LSYSTEM,SYSTEM,LDATSET,DATSET,
     &    IFILBUILDBOT,WRBUILDBOT
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:PI,RY_EV
      IMPLICIT NONE
C*--PSHIFT19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='PSHIFT')
C
C Dummy arguments
C
      INTEGER IPRINT
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ARG,CB(:,:),CJL,CJLB1,CNL,CNLB1,CV(:),E,F11,G11,P,PFAC,
     &           PSQR,WRCA(:,:,:),WRCB(:,:,:),WRCC(:,:,:),WRCD(:,:,:),
     &           WRCE(:,:,:),WRCF(:,:,:)
      CHARACTER*1 C2J
      REAL*8 CG1,CG2,CG4,CG5,CG8,COSDELTA,COTDELTA,D,DESO,DEXC,DOVR(:),
     &       DVC,DVCSQR,E1,E2,EPLOT(:),ERES(:,:),IMTSS(:,:,:),LDTAB(:,:)
     &       ,MJ,PS(:,:),PS0(:,:),PSTAB(:,:,:),RETSS(:,:,:),SINDELTA,SK,
     &       TMAX,TMIN,WGT,YMAX(2),YMIN(2)
      COMPLEX*16 CJLZ,CNLZ
      CHARACTER*80 FILNAM,SYS,YTXT1,YTXT2
      INTEGER I,IA_ERR,IC,IE,IFIL,IL,IM,IMJ,IS,ISK,ISPIN,IT,IWRIRRWF,
     &        IWRREGWF,JRCUT(0:1),JTOP,K,KAPPA,L,L1,LB1,LFILNAM,LMS,
     &        LMSOFF,LS,LSYS,LYTXT1,LYTXT2,N,NCURVES,NE,NLM,NMJ,NPAN,
     &        NSPIN
      CHARACTER*20 LEG(:),STR20
      LOGICAL RESONANCE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WRCA,WRCB,WRCC,WRCD,WRCE,WRCF,ERES,DOVR,CV,CB,LDTAB
      ALLOCATABLE PSTAB,EPLOT,PS,PS0,LEG,RETSS,IMTSS
C
      ALLOCATE (WRCA(2,2,NRMAX),WRCB(2,2,NRMAX),PS0(NK,3))
      ALLOCATE (WRCC(2,2,NRMAX),WRCD(2,2,NRMAX),ERES(NK,3))
      ALLOCATE (WRCE(2,2,NRMAX),WRCF(2,2,NRMAX),PS(NK,3))
      ALLOCATE (DOVR(NRMAX),CB(NRMAX,3),LEG(2*NLMAX),LDTAB(NEMAX,NK))
      ALLOCATE (RETSS(NEMAX,NL,2),IMTSS(NEMAX,NL,2))
      ALLOCATE (PSTAB(NEMAX,NK,3),EPLOT(NEMAX),CV(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: EPLOT')
C
      WRITE (6,99001)
C
      NE = NETAB(1)
      SYS = SYSTEM
      LSYS = LSYSTEM
      IWRIRRWF = 0
      IWRREGWF = 0
C
      CALL XMGRSUBSCRIPTS(SYS,LSYS,80)
C
      IFIL = 7
C
      DO IE = 1,NE
         EPLOT(IE) = DREAL(ETAB(IE,1))*RY_EV
      END DO
C
      IF ( IREL.LT.3 ) THEN
C=======================================================================
C                         NON-RELATIVISTIC
C=======================================================================
C
         NLM = NKM/2
         IF ( IREL.EQ.2 ) THEN
            NSPIN = 2
         ELSE
            NSPIN = 1
         END IF
C
         DO IT = 1,NT
C
            WRITE (6,99002) IT,TXT_T(IT)
C
            DVC = CTL(IT,1)
            DVCSQR = DVC*DVC
C
            IM = IMT(IT)
            JTOP = JRWS(IM)
C
            DO I = 1,JTOP
               DOVR(I) = DRDI(I,IM)/R(I,IM)
               CB(I,1) = 0D0
               CB(I,2) = BT(I,IT)
               CB(I,3) = BT(I,IT)
               CV(I) = VT(I,IT)
            END DO
C
            DO ISPIN = 1,2
               DO IL = 1,NL
                  PS0(IL,ISPIN) = 0.0D0
               END DO
            END DO
C
            YMIN(1:2) = +9999.0D0
            YMAX(1:2) = -9999.0D0
            TMIN = +9999.0D0
            TMAX = -9999.0D0
            DO IE = 1,NE
               E = DREAL(ETAB(IE,1))
C
               CALL NRSSITE(IWRREGWF,IWRIRRWF,IFILCBWF,E,P,IPRINT,TSST,
     &                      MSST,SSST,MEZZ,MEZJ)
C
               DO ISPIN = 1,NSPIN
                  LMSOFF = NLM*(ISPIN-1)
C
                  DO IL = 1,NLQ(IQAT(1,IT))
C
                     LMS = LMSOFF + IL*IL
C
                     RETSS(IE,IL,ISPIN) = -DREAL(TSST(LMS,LMS,IT))
                     IMTSS(IE,IL,ISPIN) = -DIMAG(TSST(LMS,LMS,IT))
                     COTDELTA = -DREAL(MSST(LMS,LMS,IT)/P)
C
                     D = ATAN(1D0/COTDELTA) - PS0(IL,ISPIN)
C
                     PS(IL,ISPIN) = D + PS0(IL,ISPIN)
                     IF ( ABS(D).GT.ABS(D-PI) ) PS(IL,ISPIN) = D - PI + 
     &                    PS0(IL,ISPIN)
                     IF ( ABS(D).GT.ABS(D+PI) ) PS(IL,ISPIN) = D + PI + 
     &                    PS0(IL,ISPIN)
                     PS0(IL,ISPIN) = PS(IL,ISPIN)
                     PSTAB(IE,IL,ISPIN) = PS(IL,ISPIN)
C
                     YMIN(ISPIN) = MIN(YMIN(ISPIN),PS(IL,ISPIN))
                     YMAX(ISPIN) = MAX(YMAX(ISPIN),PS(IL,ISPIN))
                     TMIN = MIN(TMIN,RETSS(IE,IL,ISPIN),
     &                      IMTSS(IE,IL,ISPIN))
                     TMAX = MAX(TMAX,RETSS(IE,IL,ISPIN),
     &                      IMTSS(IE,IL,ISPIN))
C
                  END DO
C
               END DO
C
            END DO
C
C ----------------------------------------------------------------------
C                              PHASE SHIFT
C ----------------------------------------------------------------------
C
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(E)'
            LS = LTXT_T(IT) + 5
            IF ( NSPIN.EQ.1 ) THEN
               YTXT1 = '!xd!F!sl '//STR20(1:LS)
               LYTXT1 = 9 + LS
               YTXT2 = ' '
               LYTXT2 = 1
            ELSE
               YTXT1 = '!xd!F!m{1}!S!UP!M{1}!N!sl '//STR20(1:LS)
               LYTXT1 = 27 + LS
               YTXT2 = '!xd!F!m{1}!S!DN!M{1}!N!sl '//STR20(1:LS)
               LYTXT2 = 27 + LS
            END IF
C
            CALL XMGRHEAD(DATSET,LDATSET,'pshift',6,TXT_T(IT),LTXT_T(IT)
     &                    ,FILNAM,80,LFILNAM,IFIL,NSPIN,EPLOT(1),1,
     &                    EPLOT(NE),1,YMIN(1),1,YMAX(1),1,YMIN(2),1,
     &                    YMAX(2),1,'energy (eV)',11,YTXT1,LYTXT1,YTXT2,
     &                    LYTXT2,'SPR-KKR calculations for '//
     &                    SYS(1:LSYS),25+LSYS,'Phase shift of '//
     &                    TXT_T(IT)(1:LTXT_T(IT))
     &                    //' in '//SYSTEM(1:LSYSTEM),
     &                    (15+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
            LEG(1) = 's'
            LEG(2) = 'p'
            LEG(3) = 'd'
            DO IC = 4,NL
               LEG(IC) = CHAR(ICHAR('f')+IC-4)
            END DO
C
            CALL XMGRLEGEND(IFIL,NSPIN,NL,NL,LEG,LEG)
C
            CALL XMGRCURVES(IFIL,NSPIN,NL,NL,2,1,0)
C
            DO ISPIN = 1,NSPIN
C
               WGT = 1.0D0
               DO IL = 1,NL
                  CALL XMGRTABLE((ISPIN-1),(IL-1),EPLOT,
     &                           PSTAB(1,IL,ISPIN),WGT,NE,IFIL)
               END DO
C
            END DO
C
            WRITE (6,*) '  '
            WRITE (6,99003) 'phase shift ',FILNAM(1:LFILNAM)
            CLOSE (IFIL)
C
C ----------------------------------------------------------------------
C        look for the first resonance in the valence band region
C ----------------------------------------------------------------------
C
            IF ( NL.GE.3 .AND. NSPIN.GT.1 ) THEN
               DO IL = 3,NL
                  RESONANCE = .FALSE.
C
                  DO ISPIN = 1,NSPIN
                     ERES(IL,ISPIN) = -9999.90D0
                     DO IE = 2,NE
                        E1 = DREAL(ETAB(IE-1,1))
                        E2 = DREAL(ETAB(IE,1))
                        IF ( E2.GT.1D0 ) EXIT
                        IF ( PSTAB(IE,IL,ISPIN).GT.PI/2D0 ) THEN
                           ERES(IL,ISPIN) = E1 + (E2-E1)
     &                        *(PI/2D0-PSTAB(IE-1,IL,ISPIN))
     &                        /(PSTAB(IE,IL,ISPIN)-PSTAB(IE-1,IL,ISPIN))
                           RESONANCE = .TRUE.
                           EXIT
                        END IF
                     END DO
                  END DO
C
                  IF ( RESONANCE ) THEN
                     DEXC = ERES(IL,2) - ERES(IL,1)
                     WRITE (6,99004) IL - 1,0D0,0D0,DEXC,DEXC*RY_EV
                  ELSE
                     WRITE (6,99005) IL - 1
                  END IF
C
               END DO
C
            END IF
C
C ----------------------------------------------------------------------
C                           t-matrix
C ----------------------------------------------------------------------
C
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(E)'
            LS = LTXT_T(IT) + 5
            IF ( NSPIN.EQ.1 ) THEN
               YTXT1 = '- t!sl '//STR20(1:LS)
               LYTXT1 = 7 + LS
               YTXT2 = ' '
               LYTXT2 = 1
            ELSE
               YTXT1 = '- t!m{1}!S!UP!M{1}!N!sl '//STR20(1:LS)
               LYTXT1 = 24 + LS
               YTXT2 = '- t!m{1}!S!DN!M{1}!N!sl '//STR20(1:LS)
               LYTXT2 = 24 + LS
            END IF
C
            CALL XMGRHEAD(DATSET,LDATSET,'t-matrix',8,TXT_T(IT),
     &                    LTXT_T(IT),FILNAM,80,LFILNAM,IFIL,NSPIN,
     &                    EPLOT(1),1,EPLOT(NE),1,TMIN,1,TMAX,1,TMIN,1,
     &                    TMAX,1,'energy (eV)',11,YTXT1,LYTXT1,YTXT2,
     &                    LYTXT2,'SPR-KKR calculations for '//
     &                    SYS(1:LSYS),25+LSYS,'t-matrix of '//TXT_T(IT)
     &                    (1:LTXT_T(IT))//' in '//SYSTEM(1:LSYSTEM),
     &                    (12+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
            LEG(1) = 's'
            LEG(2) = 'p'
            LEG(3) = 'd'
            DO IC = 4,NL
               LEG(IC) = CHAR(ICHAR('f')+IC-4)
            END DO
C
            CALL XMGRLEGEND(IFIL,NSPIN,NL,NL,LEG,LEG)
C
            CALL XMGRCURVES(IFIL,NSPIN,NL,NL,2,1,0)
C
            WGT = 1.0D0
            DO ISPIN = 1,NSPIN
               DO IL = 1,NL
                  CALL XMGRTABLE((ISPIN-1),(IL-1),EPLOT,
     &                           IMTSS(1,IL,ISPIN),WGT,NE,IFIL)
               END DO
               DO IL = 1,NL
                  CALL XMGRTABLE((ISPIN-1),NL+(IL-1),EPLOT,
     &                           RETSS(1,IL,ISPIN),WGT,NE,IFIL)
               END DO
            END DO
C
            WRITE (6,*) '  '
            WRITE (6,99003) 't-matrix    ',FILNAM(1:LFILNAM)
            CLOSE (IFIL)
C
C ----------------------------------------------------------------------
C                       LOGARITHMIC DERIVATIVE
C ----------------------------------------------------------------------
C***************************************************************** TO DO
            IF ( NL.LT.0 ) THEN
C
               CALL XMGRHEAD(DATSET,LDATSET,'logdrv',6,TXT_T(IT),
     &                       LTXT_T(IT),FILNAM,80,LFILNAM,IFIL,1,
     &                       EPLOT(1),1,EPLOT(NE),1,-5D0,0,+5D0,0,0D0,0,
     &                       0D0,0,'energy (eV)',11,
     &                       'D!s!xk!N!F(E)  (B!sxc!N=0)',26,' ',0,
     &                       'SPR-KKR calculations for '//SYS(1:LSYS),
     &                       25+LSYS,'Logarithmic derivative of '//
     &                       TXT_T(IT)(1:LTXT_T(IT))
     &                       //' in '//SYSTEM(1:LSYSTEM),
     &                       (26+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
               DO IL = 1,NL
                  L = IL/2
                  WRITE (C2J,'(I1)') 2*L - 1 + 2*(IL-2*L)
                  LEG(IL) = TXT_L(L)//'!s'//C2J//'/2!N'
               END DO
C
               CALL XMGRLEG1(IFIL,0,NL,LEG,0.18D0,0.83D0)
C
               CALL XMGRCURVES(IFIL,1,NL,0,2,1,0)
C
               DO IL = 1,NL
                  CALL XMGRTABLE(0,(IL-1),EPLOT,LDTAB(1,IL),1.0D0,NE,
     &                           IFIL)
               END DO
C
               WRITE (6,99003) 'log. deriv. ',FILNAM(1:LFILNAM)
               WRITE (6,*) ' '
C
               CLOSE (IFIL)
C
C ----------------------------------------------------------------------
            END IF
C***************************************************************** TO DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99006)
     &                               ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                               (((PSTAB(IE,IL,ISPIN),IE=1,
     &                               MIN(3,NE)),IL=1,NL),ISPIN=1,NSPIN)
C               WRITE (IFILBUILDBOT,99007) ROUTINE(1:LEN_TRIM(ROUTINE)),
C     &                IT,((LDTAB(IE,IL),IE=1,min(3,ne)),IL = 1,NL)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         END DO
C
      ELSE
C=======================================================================
C                           RELATIVISTIC
C=======================================================================
C
         DO IT = 1,NT
C
            WRITE (6,99002) IT,TXT_T(IT)
C
            DVC = CTL(IT,1)
            DVCSQR = DVC*DVC
C
            IM = IMT(IT)
            JTOP = JRWS(IM)
            NPAN = 1
            JRCUT(0) = 0
            JRCUT(1) = JTOP
C
            DO I = 1,JTOP
               DOVR(I) = DRDI(I,IM)/R(I,IM)
               CB(I,1) = 0D0
               CB(I,2) = BT(I,IT)
               CB(I,3) = BT(I,IT)
               CV(I) = VT(I,IT)
            END DO
C
            DO IMJ = 1,3
               DO K = 1,NK
                  PS0(K,IMJ) = 0.0D0
               END DO
            END DO
C
            YMIN(1:2) = +9999.0D0
            YMAX(1:2) = -9999.0D0
            DO IE = 1,NE
               E = DREAL(ETAB(IE,1))
C
               PSQR = E*(1.0D0+E/DVCSQR)
               P = SQRT(PSQR)
               PFAC = P/(1.0D0+E/DVCSQR)
C
               DO K = 1,NK
                  L = K/2
                  IF ( MOD(K,2).EQ.1 ) THEN
                     KAPPA = -L - 1
                     NMJ = 3
                  ELSE
                     KAPPA = +L
                     NMJ = 1
                  END IF
C
                  DO IMJ = 1,NMJ
C
                     IF ( IMJ.EQ.1 ) MJ = -0.5D0
                     IF ( IMJ.EQ.2 ) MJ = -DBLE(L) - 0.5D0
                     IF ( IMJ.EQ.3 ) MJ = +DBLE(L) + 0.5D0
C
                     CALL DIRBS(.FALSE.,CTL(IT,1),E,L,MJ,KAPPA,KAPPA,P,
     &                          CG1,CG2,CG4,CG5,CG8,CV,CB(1,IMJ),Z(IT),
     &                          R(1,IM),DRDI(1,IM),DOVR,JTOP,JRCUT,NPAN,
     &                          WRCC,WRCD,WRCE,WRCF,WRCA,WRCB,NRMAX,
     &                          .FALSE.)
C
C ------------------------------- wavefunctions at the muffin-tin-radius
C
                     N = JTOP
                     ISK = ISIGN(1,KAPPA)
                     SK = DBLE(ISK)
                     L1 = L
                     LB1 = L - ISK
                     ARG = P*R(JTOP,IM)
                     CJL = CJLZ(L1,ARG)
                     CJLB1 = CJLZ(LB1,ARG)
                     CNL = CNLZ(L1,ARG)
                     CNLB1 = CNLZ(LB1,ARG)
C
                     G11 = WRCC(1,1,N)/R(N,IM)
                     F11 = WRCD(1,1,N)/(R(N,IM)*DVC)
C
                     COSDELTA = DBLE(CNL*DVC*F11-PFAC*SK*CNLB1*G11)
                     SINDELTA = DBLE(CJL*DVC*F11-PFAC*SK*CJLB1*G11)
C
                     IF ( ABS(COSDELTA).GT.1D-8 ) THEN
                        D = ATAN(SINDELTA/COSDELTA) - PS0(K,IMJ)
                     ELSE
                        D = PI/2.0D0
                     END IF
C
                     PS(K,IMJ) = D + PS0(K,IMJ)
                     IF ( ABS(D).GT.ABS(D-PI) ) PS(K,IMJ) = D - PI + 
     &                    PS0(K,IMJ)
                     IF ( ABS(D).GT.ABS(D+PI) ) PS(K,IMJ) = D + PI + 
     &                    PS0(K,IMJ)
                     PS0(K,IMJ) = PS(K,IMJ)
                     PSTAB(IE,K,IMJ) = PS(K,IMJ)
C
                     IF ( IMJ.EQ.1 ) THEN
                        YMIN(1) = MIN(YMIN(1),PS(K,IMJ))
                        YMAX(1) = MAX(YMAX(1),PS(K,IMJ))
                     ELSE
                        YMIN(2) = MIN(YMIN(2),PS(K,IMJ))
                        YMAX(2) = MAX(YMAX(2),PS(K,IMJ))
                     END IF
C
C
                     IF ( IMJ.EQ.1 ) LDTAB(IE,K) = DVC*DREAL(F11/G11)
     &                    - 1D0
C
                  END DO
C
               END DO
C
            END DO
C
C ----------------------------------------------------------------------
C                              PHASE SHIFT
C ----------------------------------------------------------------------
C
            CALL XMGRHEAD(DATSET,LDATSET,'pshift',6,TXT_T(IT),LTXT_T(IT)
     &                    ,FILNAM,80,LFILNAM,IFIL,2,EPLOT(1),1,EPLOT(NE)
     &                    ,1,YMIN(1),1,YMAX(1),1,YMIN(2),1,YMAX(2),1,
     &                    'energy (eV)',11,'!xd!sk!N!F(E)  (B!sxc!N=0)',
     &                    26,'!xd!sk!N!F(E) (j=l!x+1/2;|m|=!Fj)',33,
     &                    'SPR-KKR calculations for '//SYS(1:LSYS),
     &                    25+LSYS,'Phase shift of '//TXT_T(IT)
     &                    (1:LTXT_T(IT))//' in '//SYSTEM(1:LSYSTEM),
     &                    (15+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
C---------------------------------------------------- lower panel pshift
C
            DO K = 1,NK
               L = K/2
               WRITE (C2J,'(I1)') 2*L - 1 + 2*(K-2*L)
               LEG(K) = TXT_L(L)//'!s'//C2J//'/2!N'
            END DO
C
            CALL XMGRLEG1(IFIL,0,NK,LEG,0.18D0,0.49D0)
C
            CALL XMGRCURVES(IFIL,2,NK,((NK/2)+1)*2,2,-1,-1)
C
            DO K = 1,NK
               CALL XMGRTABLE(0,(K-1),EPLOT,PSTAB(1,K,1),1.0D0,NE,IFIL)
            END DO
C
C---------------------------------------------------- upper panel pshift
C
            IS = 0
            DO K = 1,NK,2
               DO IMJ = 2,3
                  IS = IS + 1
                  L = K/2
                  IF ( IMJ.EQ.2 ) THEN
                     LEG(IS) = TXT_L(L)
                  ELSE
                     LEG(IS) = ' '
                  END IF
               END DO
            END DO
            NCURVES = IS
C
            CALL XMGRLEG1(IFIL,1,NCURVES,LEG,0.18D0,0.84D0)
C
            IS = -1
            DO K = 1,NK,2
               DO IMJ = 2,3
                  IS = IS + 1
                  CALL XMGRTABLE(1,IS,EPLOT,PSTAB(1,K,IMJ),1.0D0,NE,
     &                           IFIL)
               END DO
            END DO
C
            WRITE (6,99003) 'phase shift ',FILNAM(1:LFILNAM)
C
            CLOSE (IFIL)
C
C ----------------------------------------------------------------------
C        look for the first resonance in the valence band region
C ----------------------------------------------------------------------
            IF ( NK.GE.4 ) THEN
               DO K = 4,NK
                  L = K/2
                  RESONANCE = .FALSE.
                  IF ( MOD(K,2).EQ.0 ) THEN
                     NMJ = 1
                  ELSE
                     NMJ = 3
                  END IF
C
                  DO IMJ = 1,NMJ
                     ERES(K,IMJ) = -9999.90D0
                     DO IE = 2,NE
                        E1 = DREAL(ETAB(IE-1,1))
                        E2 = DREAL(ETAB(IE,1))
                        IF ( E2.GT.1D0 ) EXIT
                        IF ( PSTAB(IE,K,IMJ).GT.PI/2D0 ) THEN
                           ERES(K,IMJ) = E1 + (E2-E1)
     &                        *(PI/2D0-PSTAB(IE-1,K,IMJ))
     &                        /(PSTAB(IE,K,IMJ)-PSTAB(IE-1,K,IMJ))
                           RESONANCE = .TRUE.
                           EXIT
                        END IF
                     END DO
                  END DO
C
                  IF ( MOD(K,2).EQ.1 ) THEN
                     IF ( RESONANCE ) THEN
                        DESO = ERES(K,1) - ERES(K-1,1)
                        DEXC = ERES(K,2) - ERES(K,3)
                        WRITE (6,99004) L,DESO,DESO*RY_EV,DEXC,
     &                                  DEXC*RY_EV
                     ELSE
                        WRITE (6,99005) L
                     END IF
C
                  END IF
C
               END DO
C
            END IF
C
C ----------------------------------------------------------------------
C                       LOGARITHMIC DERIVATIVE
C ----------------------------------------------------------------------
C
            CALL XMGRHEAD(DATSET,LDATSET,'logdrv',6,TXT_T(IT),LTXT_T(IT)
     &                    ,FILNAM,80,LFILNAM,IFIL,1,EPLOT(1),1,EPLOT(NE)
     &                    ,1,-5D0,0,+5D0,0,0D0,0,0D0,0,'energy (eV)',11,
     &                    'D!s!xk!N!F(E)  (B!sxc!N=0)',26,' ',0,
     &                    'SPR-KKR calculations for '//SYS(1:LSYS),
     &                    25+LSYS,'Logarithmic derivative of '//
     &                    TXT_T(IT)(1:LTXT_T(IT))
     &                    //' in '//SYSTEM(1:LSYSTEM),
     &                    (26+LTXT_T(IT)+4+LSYSTEM),.FALSE.)
C
            DO K = 1,NK
               L = K/2
               WRITE (C2J,'(I1)') 2*L - 1 + 2*(K-2*L)
               LEG(K) = TXT_L(L)//'!s'//C2J//'/2!N'
            END DO
C
            CALL XMGRLEG1(IFIL,0,NK,LEG,0.18D0,0.83D0)
C
            CALL XMGRCURVES(IFIL,1,NK,0,2,1,0)
C
            DO K = 1,NK
               CALL XMGRTABLE(0,(K-1),EPLOT,LDTAB(1,K),1.0D0,NE,IFIL)
            END DO
C
            WRITE (6,99003) 'log. deriv. ',FILNAM(1:LFILNAM)
            WRITE (6,*) ' '
C
            CLOSE (IFIL)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT ) THEN
               WRITE (IFILBUILDBOT,99006) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                IT,(((PSTAB(IE,K,IMJ),IE=1,MIN(3,NE)),K=1,NK,2),
     &                IMJ=2,3)
               WRITE (IFILBUILDBOT,99007) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                IT,((LDTAB(IE,K),IE=1,MIN(3,NE)),K=1,NK)
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
C ----------------------------------------------------------------------
         END DO
C
C=======================================================================
      END IF
C
      DEALLOCATE (WRCA,WRCB,WRCC,WRCD,WRCE,WRCF,ERES,DOVR,CV,CB,LDTAB)
      DEALLOCATE (PSTAB,EPLOT,PS,PS0,LEG)
C
      WRITE (6,*) '          <PSHIFT> finished'
      STOP
C
99001 FORMAT (' '//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      *****           ****   *    *  *  ******  *******     *'
     &  ,/,10X,
     &  '*      *    *         *    *  *    *  *  *          *        *'
     &  ,/,10X,
     &  '*      *    *         *       *    *  *  *          *        *'
     &  ,/,10X,
     &  '*      *****    ***    ****   ******  *  *****      *        *'
     &  ,/,10X,
     &  '*      *                   *  *    *  *  *          *        *'
     &  ,/,10X,
     &  '*      *              *    *  *    *  *  *          *        *'
     &  ,/,10X,
     &  '*      *               ****   *    *  *  *          *        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,10X,'IT=',I2,2X,A)
99003 FORMAT (/,10X,A,' written to the file ',A)
99004 FORMAT (/,10X,'deduced from ',
     &        'the resonance with E_res < 1 Ry for l =',i2,/,10X,
     &        'SO-splitting ',F10.5,' Ry',F10.5,' eV',/,10X,
     &        'XC-splitting ',F10.5,' Ry',F10.5,' eV')
99005 FORMAT (/,10X,'no resonance with E_res < 1 Ry found for l =',i2)
99006 FORMAT ('# BUILDBOT: ',A,':  phase shifts          ','  for IT =',
     &        I5,/,(1PE22.14))
99007 FORMAT ('# BUILDBOT: ',A,':  logarithmic derivative','  for IT =',
     &        I5,/,(1PE22.14))
      END
