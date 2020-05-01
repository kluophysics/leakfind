C*==fpcoupl.f    processed by SPAG 6.70Rc at 08:21 on 26 Apr 2017
      SUBROUTINE FPCOUPL(IPRINT_LOC,NKM_LOC,NL_LOC)
C   ********************************************************************
C   *                                                                  *
C   *  find out the coupling of the wave functions to the              *
C   *  non-spehrical potential terms                                   *
C   *                                                                  *
C   *  in:   IREL      relativistic level                              *
C   *        NFPT      number of non- zero LM-potential terms          *
C   *        NLMFPT    highest LM to be considered for potential       *
C   *        LMIFP     LM-character of term IFP                        *
C   *        KLMFP     key = 1/0 for LM-character occurs/not           *
C   *                                                                  *
C   *  out:  NBLK      number of sub-blocks in t-matrix                *
C   *        NSOLBLK   number of solutions in block IBLK               *
C   *        IKMSOLBLK IKM-index for solution ISOL in block IBLK       *
C   *        ISOLIKM   number of solution (in block IBLK)              *
C   *                  with dominating character  IKM                  *
C   *                                                                  *
C   *        VAMEG     potential matrix element <LAM| V_L |LAM'> for g *
C   *        VAMEF     potential matrix element <LAM| V_L |LAM'> for f *
C   *                                                                  *
C   *                                                                  *
C   *  IREL <= 2:  VAMEG = <Y_LM1|Y_LM3|Y_(LM2>                        *
C   *              Gaunt coefficients for real spherical harmonics     *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  the array size  NCPLWFMAX  will be fixed and the dependent      *
C   *  arrays VAMEG, VAMEF are passed via file  IOTMP  to <KKRMAIN>    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQCLU
      USE MOD_CONSTANTS,ONLY:PI,C0
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_CALCMODE,ONLY:IREL,LHS_SOL_EQ_RHS_SOL
      USE MOD_TYPES,ONLY:TXT_T,LTXT_T,ITBOT,ITTOP,NTMAX,NLMFPMAX,NLMFPT,
     &    KLMFP,ISOLIKM,NSOLBLK,IKMSOLBLK,LMIFP,NBLK,NFPT,IKMCPLWF,NT,
     &    VAMEF,VAMEG
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,SYSTEM_DIMENSION
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,NKMMAX,WKM1,WKM2,NLM,L_IKM,AME_RLM,
     &    IDOS,ISMT,IMKM_IKM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPCOUPL')
C
C Dummy arguments
C
      INTEGER IPRINT_LOC,NKM_LOC,NL_LOC
C
C Local variables
C
      REAL*8 ALAT,RGNT(:),SQRTFP
      COMPLEX*16 AMEF(:,:,:),AMEG(:,:,:),CIPWL(:),CLEBLIN(:,:),
     &           CMAT(NKM_LOC,NKM_LOC),VAMEF_TMP(:,:,:,:,:,:),
     &           VAMEG_TMP(:,:,:,:,:,:)
      LOGICAL C_GT_EPS
      INTEGER I,IA_ERR,IBLK,IBLKIKM(:),IBLKTOP,ICH0,IERR,IFMT,IFMT_NREL,
     &        IFP,IG123,IKM,IKM1,IKM2,IKMSTACK(:),IKMTOP,IKMX,IMKM1,
     &        IMKM2,IPOT,IPOTNS,IRGNT(:),ISOL,IT,ITAB,ITABMAT(:,:),
     &        ITBOTSAV,ITTOPSAV,J,J1,J1TOP,JKM,K,KHIT,L1,L2,LCTTOT,LIN,
     &        LINCNT,LINTAB(:,:),LINTABMAX,LM1,LM1LM2,LM2,LM3,LMAX,LMFP,
     &        LMFPNS(NLMFPMAX),LRGNT12,LRGNT123,LTBLOC(:,:),MC,
     &        NCPLWFMAX,NKAME,NKMAME,NKMSTACK,NKM_PRINT,NPOTNS,NRGNT(:),
     &        NRGNT123TAB(20),NSOL,NSOLMAX,SWAP
C
C*** End of declarations rewritten by SPAG
C
      DATA NRGNT123TAB/1,15,96,388,1181,2917,6342,12452,22525,38289,
     &     61912,95914,143531,208371,294744,407644,552931,736829,966544,
     &     1250346/
C
      ALLOCATABLE CLEBLIN,AMEG,AMEF,LINTAB,LTBLOC,VAMEF_TMP,VAMEG_TMP
      ALLOCATABLE IRGNT,NRGNT,RGNT,CIPWL,IBLKIKM,IKMSTACK,ITABMAT
C
C------------------------------------------- NKAME should be > NK=2*NL-1
      NKAME = 2*NL_LOC + 1
      NKMAME = 2*((NKAME+1)/2)**2
      IF ( NKM_LOC.GT.NKMAME ) STOP 'in <FPCOUPL>  NKM > NKAME'
C
      LINTABMAX = NKMAME**2
C
      ALLOCATE (IBLKIKM(NKMMAX),IKMSTACK(NKMMAX),ITABMAT(NKMMAX,NKMMAX))
      ALLOCATE (CLEBLIN(LINTABMAX,2))
      ALLOCATE (LINTAB(LINTABMAX,2),LTBLOC(LINTABMAX,2))
      ALLOCATE (AMEG(LINTABMAX,NLMFPMAX,2))
      ALLOCATE (AMEF(LINTABMAX,NLMFPMAX,2),STAT=IA_ERR)
C
      IF ( IA_ERR.NE.0 ) STOP 'alloc:scfinitpot -> LTBLOC'
C
      SQRTFP = DSQRT(4.0D0*PI)
C
      NCPLWFMAX = NKMMAX
      MC = NCPLWFMAX
C
      ALLOCATE (VAMEF_TMP(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEF_TMP'
      ALLOCATE (VAMEG_TMP(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEG_TMP'
C
      VAMEG_TMP(:,:,:,:,:,:) = C0
      VAMEF_TMP(:,:,:,:,:,:) = C0
C
      NSOLMAX = 0
C
      WRITE (6,99004)
C
      IF ( IREL.LE.2 .OR. CHECK ) THEN
C-----------------------------------------------------------------------
C                   NON- and SCALAR RELATIVISTIC CASE
C-----------------------------------------------------------------------
C  NOTE: the non- and scalar relativistic case have been obtained by
C        modifying the relat. case without change of variable names
C
         IKMTOP = NL_LOC**2
         IFMT = 1
         J1TOP = 1
C
C ----------------------------------------- calculate Gaunt coefficients
C
         LRGNT123 = NRGNT123TAB(NL_LOC)
         LRGNT12 = (NL_LOC**2*(NL_LOC**2+1)/2)
C
         ALLOCATE (IRGNT(LRGNT123),NRGNT(LRGNT12),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:STRINIT->IRGNT'
         ALLOCATE (RGNT(LRGNT123),CIPWL((2*NL_LOC)**2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:STRINIT->RGNT'
C
         LMAX = NL_LOC - 1
C
C     compensate the prefactor  PRE=4*PI*2*PI/ALAT  in <STRGAUNT>
C
         ALAT = 4*PI*2*PI
C
         CALL STRGAUNT(LMAX,ALAT,RGNT,NRGNT,IRGNT,CIPWL,NL_LOC,IG123,
     &                 LRGNT12,LRGNT123)
C
      END IF
C
C-----------------------------------------------------------------------
C                       FULLY RELATIVISTIC CASE
C-----------------------------------------------------------------------
C
      IF ( IREL.GE.3 ) THEN
C
         IKMTOP = NKM_LOC
         IFMT = 3
         J1TOP = 2
C
      END IF
C-----------------------------------------------------------------------
C
C
C==================================================================== IT
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SYSTEM_TYPE(1:3)
     &     .NE.'VIV' ) THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         ITBOT = 1
         ITTOP = NT
      END IF
C
      DO IT = ITBOT,ITTOP
C
         CALL CINIT(LINTABMAX*NLMFPMAX*2,AMEG)
         CALL CINIT(LINTABMAX*NLMFPMAX*2,AMEF)
C
         DO IKM1 = 1,NKMMAX
            DO IKM2 = 1,NKMMAX
               ITABMAT(IKM1,IKM2) = 0
            END DO
         END DO
C
         DO LIN = 1,LINTABMAX
            DO I = 1,2
               LINTAB(LIN,I) = 0
            END DO
         END DO
C
         LCTTOT = 0
         NPOTNS = 0
C
C------------------------------------------ loop over LM-potential terms
C
         DO LMFP = 1,NLMFPT(IT)
C
            IF ( (KLMFP(LMFP,IT).GT.0) .OR. (LMFP.EQ.1) ) THEN
C
               IF ( LMFP.NE.1 ) THEN
                  NPOTNS = NPOTNS + 1
                  LMFPNS(NPOTNS) = LMFP
               END IF
C
C-----------------------------------------------------------------------
C                           IREL <= 2
C-----------------------------------------------------------------------
C
               IF ( IREL.LE.2 ) THEN
C
                  LINCNT = 0
                  IG123 = 0
                  LM1LM2 = 0
                  DO LM1 = 1,NL_LOC**2
                     DO LM2 = 1,LM1
                        LM1LM2 = LM1LM2 + 1
                        DO I = 1,NRGNT(LM1LM2)
                           IG123 = IG123 + 1
                           LM3 = IRGNT(IG123)
                           IF ( LM3.EQ.LMFP ) THEN
                              LINCNT = LINCNT + 1
                              LTBLOC(LINCNT,1) = LM1
                              LTBLOC(LINCNT,2) = LM2
                              CLEBLIN(LINCNT,1) = RGNT(IG123)
                              IF ( LM1.NE.LM2 ) THEN
                                 LINCNT = LINCNT + 1
                                 LTBLOC(LINCNT,1) = LM2
                                 LTBLOC(LINCNT,2) = LM1
                                 CLEBLIN(LINCNT,1) = RGNT(IG123)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
C
               ELSE
C
C-----------------------------------------------------------------------
C                           IREL = 3
C-----------------------------------------------------------------------
C
                  LIN = 0
                  DO IKM1 = 1,NKMAME
                     DO IKM2 = 1,NKMAME
C
                        IF ( C_GT_EPS(AME_RLM(IKM1,IKM2,2,LMFP,IDOS),
     &                       1D-12) .OR. 
     &                       C_GT_EPS(AME_RLM(IKM1,IKM2,2,LMFP,ISMT),
     &                       1D-12) ) THEN
C
                           LIN = LIN + 1
C
                           LTBLOC(LIN,1) = IKM1
                           LTBLOC(LIN,2) = IKM2
C
                           CLEBLIN(LIN,1)
     &                        = AME_RLM(IKM1,IKM2,2,LMFP,IDOS)
                           CLEBLIN(LIN,2)
     &                        = AME_RLM(IKM1,IKM2,2,LMFP,ISMT)
C
                        END IF
                     END DO
                  END DO
                  LINCNT = LIN
C
                  IF ( M_LM(LMFP).LT.0 ) LHS_SOL_EQ_RHS_SOL = .FALSE.
C
               END IF
C-----------------------------------------------------------------------
C
C--------------------------------------- coupling due to large component
C
               DO LIN = 1,LINCNT
C
                  IKM1 = LTBLOC(LIN,1)
                  IKM2 = LTBLOC(LIN,2)
C
                  IF ( IREL.LE.2 ) THEN
                     L1 = INT(SQRT(DBLE(IKM1-1)))
                     L2 = INT(SQRT(DBLE(IKM2-1)))
                  ELSE
                     L1 = L_IKM(IKM1)
                     L2 = L_IKM(IKM2)
                  END IF
C
                  IF ( IKM1.LE.IKMTOP .AND. IKM2.LE.IKMTOP ) THEN
C
                     IF ( ITABMAT(IKM1,IKM2).EQ.0 ) THEN
C
                        LCTTOT = LCTTOT + 1
                        ITABMAT(IKM1,IKM2) = LCTTOT
                        ITAB = LCTTOT
                        LINTAB(ITAB,1) = IKM1
                        LINTAB(ITAB,2) = IKM2
C
                     ELSE
                        ITAB = ITABMAT(IKM1,IKM2)
                     END IF
C
                     IF ( LMFP.EQ.1 ) THEN
                        AMEG(ITAB,LMFP,1) = CLEBLIN(LIN,1)*SQRTFP
                        AMEG(ITAB,LMFP,2) = CLEBLIN(LIN,2)*SQRTFP
                     ELSE
                        AMEG(ITAB,LMFP,1) = CLEBLIN(LIN,1)
                        AMEG(ITAB,LMFP,2) = CLEBLIN(LIN,2)
                     END IF
C
                  END IF
               END DO
C
               IF ( IREL.GE.3 ) THEN
C--------------------------------------- coupling due to minor component
C
                  DO LIN = 1,LINCNT
C
                     IMKM1 = LTBLOC(LIN,1)
                     IMKM2 = LTBLOC(LIN,2)
C
                     IKM1 = IMKM_IKM(IMKM1)
                     IKM2 = IMKM_IKM(IMKM2)
C
                     L1 = L_IKM(IKM1)
                     L2 = L_IKM(IKM2)
C
                     IF ( IKM1.LE.IKMTOP .AND. IKM2.LE.IKMTOP ) THEN
C
C ------------------------------- omit couplings between Delta l = +/- 2
C                          in spherically symmetric potential components
C
                        IF ( (LMFP.NE.1) .OR. (L1.EQ.L2) ) THEN
C
                           IF ( ITABMAT(IKM1,IKM2).EQ.0 ) THEN
                              LCTTOT = LCTTOT + 1
                              ITABMAT(IKM1,IKM2) = LCTTOT
                              ITAB = LCTTOT
                              LINTAB(ITAB,1) = IKM1
                              LINTAB(ITAB,2) = IKM2
                           ELSE
                              ITAB = ITABMAT(IKM1,IKM2)
                           END IF
C
                           IF ( LMFP.EQ.1 ) THEN
                              AMEF(ITAB,LMFP,1) = CLEBLIN(LIN,1)*SQRTFP
                              AMEF(ITAB,LMFP,2) = CLEBLIN(LIN,2)*SQRTFP
                           ELSE
                              AMEF(ITAB,LMFP,1) = CLEBLIN(LIN,1)
                              AMEF(ITAB,LMFP,2) = CLEBLIN(LIN,2)
                           END IF
C
                        END IF
C
                     END IF
                  END DO
               END IF
C-----------------------------------------------------------------------
            END IF
         END DO
C
C=======================================================================
C                     set up of    COUPLING TABLES
C=======================================================================
C-------------------------------------------- find number of blocks NBLK
C------------------------------ block number IBLKIKM(IKM)  for index IKM
C
         DO IKM = 1,NKMMAX
            IBLKIKM(IKM) = 0
         END DO
C
         NBLK(IT) = 0
         DO IKM = 1,IKMTOP
            IF ( IBLKIKM(IKM).EQ.0 ) THEN
               NBLK(IT) = NBLK(IT) + 1
               IBLKIKM(IKM) = NBLK(IT)
               IKMSTACK(1) = IKM
               NKMSTACK = 1
 10            CONTINUE
               DO I = 1,LCTTOT
                  IKM1 = LINTAB(I,1)
                  IKM2 = LINTAB(I,2)
                  DO J = 1,NKMSTACK
                     IKMX = IKMSTACK(J)
                     IF ( IKMX.EQ.IKM1 ) THEN
                        IBLKIKM(IKM2) = NBLK(IT)
                        KHIT = 0
                        DO K = 1,NKMSTACK
                           IF ( IKM2.EQ.IKMSTACK(K) ) KHIT = 1
                        END DO
                        IF ( KHIT.EQ.0 ) THEN
                           NKMSTACK = NKMSTACK + 1
                           IKMSTACK(NKMSTACK) = IKM2
                           GOTO 10
                        END IF
                     ELSE IF ( IKMX.EQ.IKM2 ) THEN
                        IBLKIKM(IKM1) = NBLK(IT)
                        KHIT = 0
                        DO K = 1,NKMSTACK
                           IF ( IKM1.EQ.IKMSTACK(K) ) KHIT = 1
                        END DO
                        IF ( KHIT.EQ.0 ) THEN
                           NKMSTACK = NKMSTACK + 1
                           IKMSTACK(NKMSTACK) = IKM1
                           GOTO 10
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END DO
C
C--------------------------- find number of solutions in block:  NSOLBLK
C-------------------------- ISOLIKM   number of solution (in block IBLK)
C---------------------------------------- with dominating character  IKM
C
         DO IKM = 1,NKMMAX
            ISOLIKM(IKM,IT) = 0
            NSOLBLK(IKM,IT) = 0
            DO IKM1 = 1,NKMMAX
               IKMSOLBLK(IKM,IKM1,IT) = 0
            END DO
         END DO
C
         DO IKM = 1,IKMTOP
            NSOLBLK(IBLKIKM(IKM),IT) = NSOLBLK(IBLKIKM(IKM),IT) + 1
            ISOLIKM(IKM,IT) = NSOLBLK(IBLKIKM(IKM),IT)
            IKMSOLBLK(ISOLIKM(IKM,IT),IBLKIKM(IKM),IT) = IKM
            NSOLMAX = MAX(NSOLMAX,ISOLIKM(IKM,IT))
         END DO
C
         IF ( NSOLMAX.GT.NCPLWFMAX ) THEN
            WRITE (6,99001) NSOLMAX,NCPLWFMAX
            STOP
         END IF
C
C
C swap IKMSOLBLK ISOLIKM into 'right' order  --------- highest ikm first
C
         DO I = 1,NBLK(IT)
            DO J = 1,NSOLBLK(I,IT) - 1
               DO J1 = J + 1,NSOLBLK(I,IT)
                  IF ( IKMSOLBLK(J1,I,IT).GT.IKMSOLBLK(J,I,IT) ) THEN
                     IKM1 = IKMSOLBLK(J1,I,IT)
                     IKM = IKMSOLBLK(J,I,IT)
                     SWAP = IKMSOLBLK(J,I,IT)
                     IKMSOLBLK(J,I,IT) = IKMSOLBLK(J1,I,IT)
                     IKMSOLBLK(J1,I,IT) = SWAP
                     SWAP = ISOLIKM(IKM,IT)
                     ISOLIKM(IKM,IT) = ISOLIKM(IKM1,IT)
                     ISOLIKM(IKM1,IT) = SWAP
                  END IF
               END DO
            END DO
         END DO
C
C=======================================================================
C                set up  <LAM|V LM|LAM'>  and  <LAM|B LM|LAM'>
C=======================================================================
C
         DO LIN = 1,LCTTOT
            IKM1 = LINTAB(LIN,1)
            IKM2 = LINTAB(LIN,2)
C
            IF ( IKM1.LE.IKMTOP .AND. IKM2.LE.IKMTOP ) THEN
               I = ISOLIKM(IKM1,IT)
               J = ISOLIKM(IKM2,IT)
               IBLK = IBLKIKM(IKM1)
               DO IPOT = 1,NPOTNS + 1
                  IF ( IPOT.GT.1 ) THEN
                     LMFP = LMFPNS(IPOT-1)
                  ELSE
                     LMFP = 1
                  END IF
C
                  DO J1 = 1,J1TOP
                     VAMEG_TMP(I,J,IPOT,J1,IBLK,IT) = AMEG(LIN,LMFP,J1)
                     VAMEG_TMP(J,I,IPOT,J1,IBLK,IT)
     &                  = DCONJG(VAMEG_TMP(I,J,IPOT,J1,IBLK,IT))
C
                     VAMEF_TMP(I,J,IPOT,J1,IBLK,IT) = AMEF(LIN,LMFP,J1)
                     VAMEF_TMP(J,I,IPOT,J1,IBLK,IT)
     &                  = DCONJG(VAMEF_TMP(I,J,IPOT,J1,IBLK,IT))
                  END DO
               END DO
            END IF
         END DO
C
C
C
C=======================================================================
C                     PRINT OUT OF COUPLING TABLES
C=======================================================================
C
         WRITE (6,99005) IT,TXT_T(IT)(1:LTXT_T(IT)),
     &                   (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                   M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
C
         WRITE (6,99006) NBLK(IT),(IBLKIKM(I),I=1,IKMTOP)
         WRITE (6,99007) (ISOLIKM(I,IT),I=1,IKMTOP)
         WRITE (6,*) ' '
         ICH0 = ICHAR('A') - 1
         DO IBLK = 1,NBLK(IT)
            NSOL = NSOLBLK(IBLK,IT)
            WRITE (6,99008) IBLK,CHAR(ICH0+IBLK),NSOL,
     &                      (IKMSOLBLK(I,IBLK,IT),I=1,NSOL)
         END DO
C
         CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
         DO I = 1,LCTTOT
            CMAT(LINTAB(I,1),LINTAB(I,2)) = IBLKIKM(LINTAB(I,1))*10
            CMAT(LINTAB(I,2),LINTAB(I,1)) = IBLKIKM(LINTAB(I,1))*10
         END DO
         CALL CMATSTRUCT('  wave function coupling scheme',CMAT,IKMTOP,
     &                   NKM_LOC,IFMT,IFMT,1,1.0D-9,-6)
C
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
C
C-----------------------------------------------------------------------
C                       scan the solution blocks
C-----------------------------------------------------------------------
C
            WKM1(:,:) = 0D0
            WKM2(:,:) = 0D0
C
            DO IBLK = 1,NBLK(IT)
C
               WRITE (6,'(///,5X,''WAVE FUNCTION COUPLING BLOCK'',I5)')
     &                IBLK
C
               CMAT(:,:) = 0D0
               DO I = 1,LCTTOT
                  IF ( IBLKIKM(LINTAB(I,1)).EQ.IBLK ) THEN
                     CMAT(LINTAB(I,1),LINTAB(I,2))
     &                  = IBLKIKM(LINTAB(I,1))*10
                     CMAT(LINTAB(I,2),LINTAB(I,1))
     &                  = IBLKIKM(LINTAB(I,1))*10
                  END IF
               END DO
C
               IF ( IREL.EQ.3 ) THEN
                  CALL CMATSTRUCT('  coupling scheme   (kappa,mu) ',
     &                            CMAT,IKMTOP,NKM_LOC,IFMT,IFMT,1,
     &                            1.0D-9,6)
C
                  CALL CHANGEREP(NKM_LOC,NKM_LOC,CMAT,'REL>CLM',WKM2)
C
                  DO I = 1,NKM_LOC
                     DO J = 1,NKM_LOC
                        IF ( ABS(WKM2(I,J)).GT.1D-5 ) WKM2(I,J)
     &                       = IBLK*10
                     END DO
                  END DO
                  CMAT(1:NKM_LOC,1:NKM_LOC) = WKM2(1:NKM_LOC,1:NKM_LOC)
                  IFMT_NREL = 2
               ELSE
                  IFMT_NREL = IREL
               END IF
C
               CALL CMATSTRUCT('  coupling scheme   (l,ml,ms)  ',CMAT,
     &                         IKMTOP,NKM_LOC,IFMT_NREL,IFMT_NREL,1,
     &                         1.0D-9,6)
C
               WKM1(1:NKM_LOC,1:NKM_LOC) = WKM1(1:NKM_LOC,1:NKM_LOC)
     &            + CMAT(1:NKM_LOC,1:NKM_LOC)
C
            END DO
C
            CMAT(1:NKM_LOC,1:NKM_LOC) = WKM1(1:NKM_LOC,1:NKM_LOC)
C
            CALL CMATSTRUCT('  ALL couplings     (l,ml,ms)  ',CMAT,
     &                      IKMTOP,NKM_LOC,IFMT_NREL,IFMT_NREL,1,1.0D-9,
     &                      -6)
C
C-----------------------------------------------------------------------
C                         scan potential terms
C-----------------------------------------------------------------------
C
            DO IPOT = 1,NPOTNS + 1
C
               IF ( IPOT.GT.1 ) THEN
                  LMFP = LMFPNS(IPOT-1)
               ELSE
                  LMFP = 1
               END IF
               WRITE (6,'(///,5X,''POTENTIAL COUPLING MATRIX'',4I5)')
     &                IPOT,LMFP,L_LM(LMFP),M_LM(LMFP)
C
               CMAT(:,:) = 0D0
C
               DO LIN = 1,LCTTOT
                  IKM1 = LINTAB(LIN,1)
                  IKM2 = LINTAB(LIN,2)
C
                  IF ( IKM1.LE.IKMTOP .AND. IKM2.LE.IKMTOP ) THEN
                     I = ISOLIKM(IKM1,IT)
                     J = ISOLIKM(IKM2,IT)
                     IBLK = IBLKIKM(IKM1)
C
                     J1 = 1
                     CMAT(IKM1,IKM2) = VAMEG_TMP(I,J,IPOT,J1,IBLK,IT)
                     CMAT(IKM2,IKM1) = DCONJG(CMAT(IKM1,IKM2))
C
                  END IF
               END DO
C
               IF ( IREL.EQ.3 ) THEN
C
                  CALL CMATSTRUCT('  potential terms   (kappa,mu) ',
     &                            CMAT,IKMTOP,NKM_LOC,IFMT,IFMT,1,
     &                            1.0D-9,6)
C
                  CALL CHANGEREP(NKM_LOC,NKM_LOC,CMAT,'REL>RLM',WKM2)
C
                  CMAT(1:NKM_LOC,1:NKM_LOC) = WKM2(1:NKM_LOC,1:NKM_LOC)
                  IFMT_NREL = 2
                  NKM_PRINT = 2*NLM
                  IFMT_NREL = 1
                  NKM_PRINT = NLM
               ELSE
                  IFMT_NREL = IREL
                  NKM_PRINT = NLM
               END IF
C
               CALL CMATSTRUCT('  potential terms   (l,ml,ms)  ',CMAT,
     &                         NKM_PRINT,NKM_LOC,IFMT_NREL,IFMT_NREL,1,
     &                         1.0D-9,6)
C
            END DO
C
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C=======================================================================
C                PRINT OUT OF  <LAM|V LM|LAM'>  IF REQUESTED
C=======================================================================
C
         IF ( IPRINT_LOC.GT.0 ) THEN
C
            CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
            DO J = 1,IKMTOP
               DO I = 1,IKMTOP
                  IF ( ITABMAT(I,J).EQ.0 ) THEN
                     CMAT(I,J) = 0
                  ELSE
                     CMAT(I,J) = 2
                  END IF
               END DO
            END DO
            CALL CMATSTRUCT('coupling matrix  ITAB',CMAT,IKMTOP,NKM_LOC,
     &                      IFMT,IFMT,1,1.0D-9,-6)
C
            DO IPOTNS = 1,NPOTNS
               IPOT = IPOTNS + 1
               LMFP = LMFPNS(IPOTNS)
C
C ------------------------------------------------------ large component
C
               CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
               DO IBLK = 1,NBLK(IT)
                  NSOL = NSOLBLK(IBLK,IT)
                  DO J = 1,NSOL
                     JKM = IKMSOLBLK(J,IBLK,IT)
                     DO I = 1,NSOL
                        IKM = IKMSOLBLK(I,IBLK,IT)
                        CMAT(IKM,JKM) = VAMEG_TMP(I,J,IPOT,1,IBLK,IT)
                     END DO
                  END DO
               END DO
               WRITE (6,99011) 'V Clebsch Gordon matrix for LM = ',LMFP
               CALL CMATSTRUCT('<LAM|V LM|LAM''>',CMAT,IKMTOP,NKM_LOC,
     &                         IFMT,IFMT,1,1.0D-9,6)
C
               IF ( IREL.GE.3 ) THEN
                  CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
                  DO IBLK = 1,NBLK(IT)
                     NSOL = NSOLBLK(IBLK,IT)
                     DO J = 1,NSOL
                        JKM = IKMSOLBLK(J,IBLK,IT)
                        DO I = 1,NSOL
                           IKM = IKMSOLBLK(I,IBLK,IT)
                           CMAT(IKM,JKM) = VAMEG_TMP(I,J,IPOT,2,IBLK,IT)
                        END DO
                     END DO
                  END DO
                  WRITE (6,99011) 'B Clebsch Gordon matrix for LM = ',
     &                            LMFP
                  CALL CMATSTRUCT('<LAM|B LM|LAM''>',CMAT,IKMTOP,
     &                            NKM_LOC,3,3,1,1.0D-9,6)
C
C ------------------------------------------------------ small component
C
                  CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
                  DO IBLK = 1,NBLK(IT)
                     NSOL = NSOLBLK(IBLK,IT)
                     DO J = 1,NSOL
                        JKM = IKMSOLBLK(J,IBLK,IT)
                        DO I = 1,NSOL
                           IKM = IKMSOLBLK(I,IBLK,IT)
                           CMAT(IKM,JKM) = VAMEF_TMP(I,J,IPOT,1,IBLK,IT)
                        END DO
                     END DO
                  END DO
                  WRITE (6,99011) 'V Clebsch Gordon matrix for LM = ',
     &                            LMFP
                  CALL CMATSTRUCT('<-LAM|V LM|-LAM''>',CMAT,IKMTOP,
     &                            NKM_LOC,IFMT,IFMT,1,1.0D-9,6)
                  CALL CINIT(NKM_LOC*NKM_LOC,CMAT)
                  DO IBLK = 1,NBLK(IT)
                     NSOL = NSOLBLK(IBLK,IT)
                     DO J = 1,NSOL
                        JKM = IKMSOLBLK(J,IBLK,IT)
                        DO I = 1,NSOL
                           IKM = IKMSOLBLK(I,IBLK,IT)
                           CMAT(IKM,JKM) = VAMEF_TMP(I,J,IPOT,2,IBLK,IT)
                        END DO
                     END DO
                  END DO
                  WRITE (6,99011) 'B Clebsch Gordon matrix for LM = ',
     &                            LMFP
                  CALL CMATSTRUCT('<-LAM|B LM|-LAM''>',CMAT,IKMTOP,
     &                            NKM_LOC,IFMT,IFMT,1,1.0D-9,6)
               END IF
C
            END DO
C
         END IF
C
C=======================================================================
C                        CHECK CONSISTENCY
C=======================================================================
C
         IERR = 0
C
         DO IKM = 1,IKMTOP
            IBLK = IBLKIKM(IKM)
            ISOL = ISOLIKM(IKM,IT)
            JKM = IKMSOLBLK(ISOL,IBLK,IT)
            IF ( IKM.NE.JKM ) THEN
               WRITE (6,99002) I,IBLKIKM(LINTAB(I,1)),LINTAB(I,1),
     &                         IBLKIKM(LINTAB(I,2)),LINTAB(I,2)
               IERR = 1
            END IF
         END DO
C
         IF ( NPOTNS+1.NE.NFPT(IT) ) IERR = 1
         DO IFP = 1,NPOTNS
            IF ( LMFPNS(IFP).NE.LMIFP(IFP+1,IT) ) IERR = 1
         END DO
         DO I = 1,LCTTOT
            IF ( IBLKIKM(LINTAB(I,1)).NE.IBLKIKM(LINTAB(I,2)) ) THEN
               WRITE (6,99003) I,IBLKIKM(LINTAB(I,1)),LINTAB(I,1),
     &                         IBLKIKM(LINTAB(I,2)),LINTAB(I,2)
               IERR = 1
            END IF
         END DO
         IF ( IERR.EQ.1 ) THEN
            WRITE (6,99010) 'IT        ',IT
            WRITE (6,99010) 'NLMFPT(IT)',NLMFPT(IT)
            WRITE (6,99010) 'NKM       ',NKM_LOC
            WRITE (6,99009) 'KLMFP(IT) ',(KLMFP(I,IT),I=1,NLMFPT(IT))
            WRITE (6,99010) 'LCTTOT    ',LCTTOT
            WRITE (6,99010) 'NPOTNS    ',NPOTNS
            WRITE (6,99010) 'NFP(IT)   ',NFPT(IT)
            WRITE (6,99010) 'LMFPNS    ',(LMFPNS(I),I=1,NPOTNS)
            WRITE (6,99010) 'LMIFP(IT) ',(LMIFP(I+1,IT),I=1,NPOTNS)
            WRITE (6,99010) 'LINTAB 1  ',(LINTAB(I,1),I=1,LCTTOT)
            WRITE (6,99010) 'LINTAB 2  ',(LINTAB(I,2),I=1,LCTTOT)
            STOP 'in <FPCOUPL> because of inconsistencies'
         END IF
C=======================================================================
C
      END DO
C==================================================================== IT
C
      IF ( IREL.EQ.3 ) WRITE (6,99012) LHS_SOL_EQ_RHS_SOL
C
C ======================================================================
C    copy temporary arrays  VAMEG_TMP, VAMEF_TMP  onto  VAMEG, VAMEF
C ======================================================================
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C              NEW VERSION   -----   REPLACES TRANSFER VIA FILE
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IF ( NKMMAX.GT.999999 ) THEN
         IF ( NQCLU.EQ.0 ) THEN
            NCPLWFMAX = NSOLMAX
         ELSE
            NCPLWFMAX = NKMMAX
         END IF
         IBLKTOP = IKMTOP
C
         IF ( ALLOCATED(VAMEG) ) DEALLOCATE (VAMEG,VAMEF,IKMCPLWF)
C
         MC = NCPLWFMAX
         ALLOCATE (VAMEG(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEG'
         ALLOCATE (IKMCPLWF(NCPLWFMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> IKMCPLWF'
         IKMCPLWF = 0
C
         VAMEG(:,:,:,:,:,:) = C0
C
         DO IT = ITBOT,ITTOP
            DO IBLK = 1,IBLKTOP
               DO J1 = 1,J1TOP
                  DO IPOT = 1,NLMFPMAX
                     DO J = 1,NCPLWFMAX
                        DO I = 1,NCPLWFMAX
                           VAMEG(I,J,IPOT,J1,IBLK,IT)
     &                        = VAMEG_TMP(I,J,IPOT,J1,IBLK,IT)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
         IF ( IREL.LE.2 ) THEN
C
            ALLOCATE (VAMEF(1,1,1,1,1,1),STAT=IA_ERR)
            VAMEF = 999999D0
C
         ELSE
C
            ALLOCATE (VAMEF(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEF'
C
            VAMEF(:,:,:,:,:,:) = C0
C
            DO IT = ITBOT,ITTOP
               DO IBLK = 1,IBLKTOP
                  DO J1 = 1,J1TOP
                     DO IPOT = 1,NLMFPMAX
                        DO J = 1,NCPLWFMAX
                           DO I = 1,NCPLWFMAX
                              VAMEF(I,J,IPOT,J1,IBLK,IT)
     &                           = VAMEF_TMP(I,J,IPOT,J1,IBLK,IT)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
         END IF
C
      END IF
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C>>>>>>>>>>>>>>>>>>>>>>> WILL GET OBSOLETE
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ======================================================================
C       write tables to temporary transfer file - read in <KKRMAIN>
C ======================================================================
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      NCPLWFMAX = NSOLMAX
      IBLKTOP = IKMTOP
      WRITE (IOTMP) NCPLWFMAX,IBLKTOP,J1TOP
      WRITE (IOTMP) ((((((VAMEG_TMP(I,J,IPOT,J1,IBLK,IT),I=1,NCPLWFMAX),
     &              J=1,NCPLWFMAX),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,
     &              IBLKTOP),IT=ITBOT,ITTOP)
      IF ( IREL.GE.3 ) WRITE (IOTMP) ((((((VAMEF_TMP(I,J,IPOT,J1,IBLK,IT
     &                               ),I=1,NCPLWFMAX),J=1,NCPLWFMAX),
     &                               IPOT=1,NLMFPMAX),J1=1,J1TOP),
     &                               IBLK=1,IBLKTOP),IT=ITBOT,ITTOP)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                        WILL GET OBSOLETE  <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C-----------------------------------------------------------------------
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SYSTEM_TYPE(1:3)
     &     .NE.'VIV' ) THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
      END IF
C
      DEALLOCATE (CLEBLIN,AMEG,AMEF,LINTAB,LTBLOC,VAMEF_TMP,VAMEG_TMP)
      DEALLOCATE (IBLKIKM,IKMSTACK,ITABMAT)
C
      IF ( ALLOCATED(IRGNT) ) DEALLOCATE (IRGNT,NRGNT,RGNT,CIPWL)
C
C-----------------------------------------------------------------------
C
99001 FORMAT (/,1X,79('*'),/,10X,'STOP in <FPCOUPL>',/,10X,
     &        'maximum number of coupled wave functions  NSOLMAX   =',
     &        I3,/,10X,'available array size',22X,'NCPLWFMAX =',I3,/)
99002 FORMAT ('>>>>>>> JKM =',I3,' <>  IKM =',I3,' IBLK=',I3,
     &        '    ISOL=',I3)
99003 FORMAT ('>>>>>>> IBLK1 <> IBLK2:  I =',I4,' IBLK1=',I2,
     &        '  LINTAB1=',I3,' IBLK2=',I2,'  LINTAB2=',I3)
99004 FORMAT (/,1X,79('*'),/,35X,'<FPCOUPL>',/,1X,79('*'),/)
99005 FORMAT (/,10X,'IT =',I3,3X,A,//,10X,'LM  ',
     &        5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99006 FORMAT (/,10X,'NBLK      ',I3,/,10X,'IBLKIKM   ',18I3,:,/,
     &        (20X,18I3))
99007 FORMAT (10X,'ISOLIKM   ',18I3,:,/,(20X,18I3))
99008 FORMAT (10X,'IBLK =',I3,2X,A,4X,'NSOL =',I3,4X,'IKM: ',10I3,:,/,
     &        (44X,10I3))
99009 FORMAT (2X,A,25I2,:,/,(2X,10X,25I2))
99010 FORMAT (2X,A,18I3,:,/,(2X,10X,18I3))
99011 FORMAT (/,2X,A,18I3,:,/,(2X,10X,18I3))
99012 FORMAT (/,10X,'LHS_SOL_EQ_RHS_SOL = ',L1,//)
      END
