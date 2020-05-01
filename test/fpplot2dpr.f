C*==fpplot2dpr.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPPLOT2DPR
C   ********************************************************************
C   *                                                                  *
C   *   create the full potential and charge densities                 *
C   *   for 2D plotting on a given rectangular plot window             *
C   *                                                                  *
C   *   AVEC, BVEC  specify rectangular window                         *
C   *   NA,   NB    specify grid aling AVEC and BVEC                   *
C   *   DVEC        gives shift of window away from origin             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NL
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,ABAS,ADAINV,A5BAS,ALAT,
     &    SYSTEM_DIMENSION
      USE MOD_SITES,ONLY:QBAS,ITOQ,NQCLU,NQHOST,N5VEC_QCLU,IQ_QCLU,NQ_L,
     &    NQ_I,IMQ
      USE MOD_CALCMODE,ONLY:NONMAG
      USE MOD_TYPES,ONLY:BNST,VNST,BT,VT,RHO2NS,LMIFP,NFPT,IMT,NLFPMAX,
     &    NLMFPMAX,NTMAX,NT,JALF_LMCT,ITBOT,ITTOP
      USE MOD_FILES,ONLY:IOTMP,PLOT2DPR,LDATSET0,DATSET0,PLOT_JALF,
     &    PLOT_JORB,FOUND_INTEGER,N_FOUND
      USE MOD_RMESH,ONLY:FULLPOT,NPAN,JRCUT,JRNS1,R,JRNSMIN,NRMAX,JRCRI,
     &    JRWS
      USE MOD_SYMMETRY,ONLY:IQREPQ,MROTR,ISYMGENQ
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      IMPLICIT NONE
C*--FPPLOT2DPR29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPPLOT2DPR')
C
C Local variables
C
      REAL*8 ADOTA,ADOTB,AHAT(3,3),ALNG(3),APLOT(:),AVEC(3),BPLOT(:),
     &       BVEC(3),D2VEC,DAVEC(3),DBVEC(3),DELSQ,DELSQMIN,DVEC(3),
     &       JALF3,JPLOT(:,:,:),NORM,NORMA,NORMB,NORMDA,NORMDB,PC0VEC(3)
     &       ,PC1VEC(3),PC50VEC(5),PC51VEC(5),POT3,POTPLOT(:,:,:),
     &       PR0VEC(3),PR1VEC(3),PR2VEC(3),PR3,PR3HAT(3),PR3VEC(3),
     &       R12VEC(3),R3VEC(3),RCRIT,RHO3,RHOPLOT(:,:,:),RIVEC(3),RR,
     &       RYLM(:),TMPVEC(3),W3A,W3B,WR,X(3),XX,YY
      CHARACTER*80 AVECSTR,BVECSTR,FILNAM
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER I,I1,I2,I3,I5,IAVEC,IBVEC,IC,IFP,IM,IM3,IPLOT,IQ,IQ2,IQ3,
     &        IQA,IQB,IQCLU,IR,IR3A,IR3B,IRTOP,IRTOPIN,IT,IT3,ITIN,
     &        LFILNAM,LM,LMAX,LMR,LMR_MAX,LMR_MAXIN,LPLOT_MAX,LR_MAX,
     &        LR_MAXIN,NAVEC,NBVEC,NLMPLOT_MAX,NPLOT,NPTS,NPTS_INSIDE,
     &        PN0VEC(3),PN50VEC(5),PN53VEC(5)
      LOGICAL PLOT_J,PLOT_V_RHO,PR_INSIDE,SAME
C
C*** End of declarations rewritten by SPAG
C
      DATA AVEC/2D0,0D0,0D0/,BVEC/0D0,2D0,0D0/,DVEC/ - 1D0, - 1D0,0D0/
      DATA NAVEC/50/,NBVEC/50/
C
      ALLOCATABLE RYLM,RHOPLOT,POTPLOT,JPLOT,APLOT,BPLOT
C
      PLOT_V_RHO = PLOT2DPR(1) .OR. PLOT2DPR(1)
      PLOT_J = PLOT_JALF .OR. PLOT_JORB
C
      IF ( PLOT_V_RHO .AND. PLOT_J ) CALL STOP_ERROR(ROUTINE,
     &     'PLOT_V_RHO .and. PLOT_J not allowed')
C
      IF ( .NOT.FULLPOT .AND. PLOT_V_RHO ) RETURN
C
C     copy spherical parts of variables to nonspherical variables
C
      IF ( FULLPOT .AND. PLOT_V_RHO ) THEN
         VNST(JRNSMIN:NRMAX,1,1:NTMAX) = VT(JRNSMIN:NRMAX,1:NTMAX)
     &      *SQRT_4PI
         BNST(JRNSMIN:NRMAX,1,1:NTMAX) = BT(JRNSMIN:NRMAX,1:NTMAX)
     &      *SQRT_4PI
      END IF
C
C----------------------------------- set NPLOT=2 for spin-polarized case
      IF ( NONMAG ) THEN
         NPLOT = 1
      ELSE
         NPLOT = 2
      END IF
C=======================================================================
C                            set plot window
C=======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',1)
C
      CALL SECTION_SET_REAL_ARRAY('AVEC',AVEC,N_FOUND,3,0,9999D0,1)
C
      CALL SECTION_SET_REAL_ARRAY('BVEC',BVEC,N_FOUND,3,0,9999D0,1)
C
      ADOTB = DDOT(3,AVEC,1,BVEC,1)
      IF ( ABS(ADOTB).GT.1D-12 ) THEN
         NORMB = DNRM2(3,BVEC,1)
         ADOTA = DDOT(3,AVEC,1,AVEC,1)
         BVEC(1:3) = BVEC(1:3) - (ADOTB/ADOTA)*AVEC(1:3)
         NORM = DNRM2(3,BVEC,1)
         BVEC(1:3) = BVEC(1:3)*(NORMB/NORM)
      END IF
C
      CALL SECTION_SET_REAL_ARRAY('DVEC',DVEC,N_FOUND,3,0,9999D0,0)
C
      CALL SECTION_SET_INTEGER('NA',NAVEC,9999,0)
      CALL SECTION_SET_INTEGER('NB',NBVEC,9999,0)
C
      NORMA = DNRM2(3,AVEC,1)
      NORMB = DNRM2(3,BVEC,1)
C
      IF ( .NOT.FOUND_INTEGER ) NBVEC = NINT(NAVEC*NORMB/NORMA)
C
      DAVEC(1:3) = AVEC(1:3)/DBLE(NAVEC)
      DBVEC(1:3) = BVEC(1:3)/DBLE(NBVEC)
C
      NORMDA = DNRM2(3,DAVEC,1)
      NORMDB = DNRM2(3,DBVEC,1)
C
      WRITE (6,99001) AVEC,NORMA,NAVEC,NORMDA,BVEC,NORMB,NBVEC,NORMDB,
     &                DVEC
C
      ALLOCATE (APLOT(0:NAVEC),BPLOT(0:NBVEC))
C
      APLOT(:) = 0D0
      BPLOT(:) = 0D0
C
C VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      IF ( PLOT_V_RHO ) THEN
C
         LPLOT_MAX = NLFPMAX - 1
         NLMPLOT_MAX = NLMFPMAX
C
         ALLOCATE (POTPLOT(0:NAVEC,0:NBVEC,NPLOT))
         ALLOCATE (RHOPLOT(0:NAVEC,0:NBVEC,NPLOT))
C
         POTPLOT(:,:,:) = 0D0
         RHOPLOT(:,:,:) = 0D0
C
C*********************************************************************
C*********************************************************************
C*********************************************************************
         IF ( NT.LT.0 ) THEN
C
            DO I = 1,3
               ALNG(I) = SQRT(DDOT(3,ABAS(1,I),1,ABAS(1,I),1))
               AHAT(1:3,I) = ABAS(1:3,I)/ALNG(I)
            END DO
C
            WRITE (*,'(/,A,/,5(3f8.3,/))') ' ABAS   ',ABAS
            WRITE (*,'(/,A,/,5(3f8.3,/))') ' ADAINV ',ADAINV
            WRITE (*,'(/,A,/,5(3f8.3,/))') ' A5BAS ',A5BAS
            WRITE (*,'(/,A,/,5(3f8.3,/))') ' AHAT   ',AHAT
            WRITE (*,'(/,A,/,5(f8.3,/))') ' ALNG   ',ALNG
C
            DO IT = 1,NT
               DO IFP = 1,NFPT(IT)
                  LM = LMIFP(IFP,IT)
                  IF ( LM.GT.11111 ) THEN
                     RHO2NS(1:NRMAX,LM,IT,1:2) = 0D0
                     VNST(JRNSMIN:NRMAX,LM,IT) = 0D0
                     BNST(JRNSMIN:NRMAX,LM,IT) = 0D0
                  END IF
               END DO
            END DO
         END IF
C*********************************************************************
C*********************************************************************
C*********************************************************************
C
C
C VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      ELSE
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
         IF ( PLOT_JALF ) THEN
            FILNAM = DATSET0(1:LDATSET0)//'_JALF.dat'
         ELSE
            FILNAM = DATSET0(1:LDATSET0)//'_JORB.dat'
         END IF
C
         LMAX = NL - 1
         LR_MAX = 2*LMAX + 1
         LMR_MAX = (LR_MAX+1)**2
C
         LPLOT_MAX = LR_MAX
         NLMPLOT_MAX = LMR_MAX
C
         ALLOCATE (JPLOT(0:NAVEC,0:NBVEC,3))
         JPLOT(:,:,:) = 0D0
         ALLOCATE (JALF_LMCT(NRMAX,LMR_MAX,3,NTMAX))
C
         OPEN (IOTMP,FILE=FILNAM)
         CALL RDHEAD(IOTMP)
         READ (IOTMP,*)
C
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
C
            READ (IOTMP,*)
            READ (IOTMP,*) ITIN,IRTOPIN,LR_MAXIN,LMR_MAXIN
            IF ( IT.NE.ITIN .OR. IRTOP.NE.IRTOPIN .OR. 
     &           LR_MAX.NE.LR_MAXIN .OR. LMR_MAX.NE.LMR_MAXIN )
     &           CALL STOP_ERROR(ROUTINE,'data file not consistent')
            WRITE (*,*) ITIN,IRTOPIN,LR_MAXIN,LMR_MAXIN
C
            DO LMR = 1,LMR_MAX
               DO IR = 1,IRTOP
                  READ (IOTMP,'(3E25.10)')
     &                  (JALF_LMCT(IR,LMR,IC,IT),IC=1,3)
               END DO
            END DO
C
         END DO
C
         CLOSE (IOTMP)
      END IF
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
C
C
      ALLOCATE (RYLM(NLMPLOT_MAX))
C
C
      NPTS = 0
      NPTS_INSIDE = 0
C=======================================================================
C                     scan the plot grid
C=======================================================================
      DO IAVEC = 0,NAVEC
         APLOT(IAVEC) = IAVEC*NORMDA
C
         DO IBVEC = 0,NBVEC
            BPLOT(IBVEC) = IBVEC*NORMDB
C
            NPTS = NPTS + 1
C
C------------------------------------------------------------- PR_INSIDE
C----- indicates PR0VEC is inside regime for which V and RHO are defined
            PR_INSIDE = .TRUE.
C
            PR0VEC(1:3) = IAVEC*DAVEC(1:3) + IBVEC*DBVEC(1:3)
     &                    + DVEC(1:3)
C
            XX = PR0VEC(1)
            YY = PR0VEC(2)
            RR = SQRT(XX*XX+YY*YY)
            IF ( RR.GT.1D-8 ) THEN
               XX = -PR0VEC(2)/RR
               YY = +PR0VEC(1)/RR
            ELSE
               XX = 0D0
               YY = 0D0
            END IF
C
C-----------------------------------------------------------------------
C      expand plot vector PR0VEC w.r.t basisvectors ABAS  and
C      modify expansion coefficients PC0VEC to get
C      equivalent plot vector PR1VEC in unit cell
C-----------------------------------------------------------------------
C
C=======================================================================
            IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C=======================================================================
C
               CALL RVECEXPAND(PR0VEC,ABAS,ADAINV,PC0VEC)
C
C-----------------------------------------------------------------------
C      get vector PR1VEC within the unit cell with  0 <= C1 < 1
C-----------------------------------------------------------------------
C
               DO I = 1,3
                  PN0VEC(I) = INT(PC0VEC(I)+1000D0) - 1000
                  PC1VEC(I) = PC0VEC(I) - PN0VEC(I)
               END DO
C
               CALL RVECLCRVB(3,PC1VEC,ABAS,PR1VEC)
C
C-----------------------------------------------------------------------
C      scan all sites IQ within central and neighboring unit cells
C      to find closest atomic site IQ2 with distance vector PR2VEC
C      PN53VEC  gives unit cell for site IQ2 w.r.t. origin
C      used to check for occupation by embedded cluster
C-----------------------------------------------------------------------
               PN53VEC(4:5) = 0
               DELSQMIN = 1D20
               IQ2 = 0
               DO I1 = -1,1
                  DO I2 = -1,1
                     DO I3 = -1,1
C
                        CALL RVECLCIB(I1,I2,I3,ABAS,RIVEC)
C
                        DO IQ = 1,NQHOST
C
                           DELSQ = 0D0
                           DO I = 1,3
                              TMPVEC(I) = PR1VEC(I)
     &                           - (RIVEC(I)+QBAS(I,IQ))
                              DELSQ = DELSQ + TMPVEC(I)*TMPVEC(I)
                           END DO
C
                           IF ( DELSQ.LT.DELSQMIN ) THEN
                              DELSQMIN = DELSQ
                              IQ2 = IQ
                              PR2VEC(1:3) = TMPVEC(1:3)
                              PN53VEC(1) = PN0VEC(1) + I1
                              PN53VEC(2) = PN0VEC(2) + I2
                              PN53VEC(3) = PN0VEC(3) + I3
                           END IF
C
                        END DO
C
                     END DO
                  END DO
               END DO
C
               IF ( IQ2.EQ.0 ) CALL STOP_TRACE_BACK(ROUTINE,'IQ2 = 0')
C
C=======================================================================
            ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C=======================================================================
C
               CALL RVEC5EXPAND(PR0VEC,ABAS,ADAINV,A5BAS,PN50VEC,
     &                          PC50VEC)
C
C........................................................................
               X(1:3) = 0D0
               DO I = 1,5
                  X(1:3) = X(1:3) + (PN50VEC(I)+PC50VEC(I))*A5BAS(1:3,I)
               END DO
               DO I = 1,3
                  IF ( ABS(X(I)-PR0VEC(I)).GT.1D-8 ) THEN
                     WRITE (*,'(a)') 'RVEC5EXPAND A'
                     WRITE (*,'(a,5f8.4)') '######## c5   ',PC50VEC(1:5)
                     WRITE (*,'(a,3f8.4)') '######## r(0) ',PR0VEC(1:3)
                     WRITE (*,'(a,3f8.4)') '######## x(0) ',X(1:3)
                     CALL STOP_MESSAGE(ROUTINE,'RVEC5EXPAND A')
                  END IF
               END DO
C........................................................................
C-----------------------------------------------------------------------
C      get vector PR1VEC within the unit cell
C-----------------------------------------------------------------------
C
               PC51VEC(1:5) = PC50VEC(1:5)
C
               CALL RVECLCRVB(5,PC51VEC,A5BAS,PR1VEC)
C
               CALL RVECEXPAND(PR1VEC,ABAS,ADAINV,X)
               DO I = 1,3
                  IF ( (X(I).LT.-1D-8) .OR. X(I).GE.1D0 ) THEN
                     WRITE (*,'(a)') 'RVEC5EXPAND B'
                     WRITE (*,'(a,i8)') '######## i    ',I
                     WRITE (*,'(a,5f12.8)') '######## x    ',X
                     WRITE (*,'(a,3f8.4)') '######## r(1) ',PR1VEC(1:3)
                     WRITE (*,'(a,3f8.4)') '######## r(0) ',PR0VEC(1:3)
                     WRITE (*,'(a,5f8.4)') '######## c5(0)',PC50VEC(1:5)
                     WRITE (*,'(a,5i8)') '######## n5(0)',PN50VEC(1:5)
                     CALL STOP_MESSAGE(ROUTINE,'RVEC5EXPAND B')
                  END IF
               END DO
C
C-----------------------------------------------------------------------
C      scan all sites IQ within central and neighboring unit cells
C      to find closest atomic site IQ2 with distance vector PR2VEC
C      PN53VEC  gives unit cell for site IQ2 w.r.t. origin
C      used to check for occupation by embedded cluster
C      2D: cluster sites allowed only in I-zone, i.e. PN53VEC(3)=0
C-----------------------------------------------------------------------
               PN53VEC(4:5) = 0
               DELSQMIN = 1D20
               IQ2 = 0
               DO I1 = -1,1
                  DO I2 = -1,1
C
                     CALL RVECLCIB(I1,I2,0,ABAS,R12VEC)
C
                     DO I3 = -1,1
C
                        IF ( I3.EQ.-1 ) THEN
                           IQA = 1
                           IQB = NQ_L
                           R3VEC(1:3) = -A5BAS(1:3,4)
                        ELSE IF ( I3.EQ.0 ) THEN
                           IQA = 1
                           IQB = NQHOST
                           R3VEC(1:3) = 0D0
                        ELSE
                           IQA = NQ_L + NQ_I + 1
                           IQB = NQHOST
                           R3VEC(1:3) = +A5BAS(1:3,5)
                        END IF
C
                        DO IQ = IQA,IQB
C
                           DELSQ = 0D0
                           DO I = 1,3
                              TMPVEC(I) = PR1VEC(I)
     &                           - (R12VEC(I)+R3VEC(I)+QBAS(I,IQ))
                              DELSQ = DELSQ + TMPVEC(I)*TMPVEC(I)
                           END DO
C
                           IF ( DELSQ.LT.DELSQMIN ) THEN
                              DELSQMIN = DELSQ
                              IQ2 = IQ
                              PR2VEC(1:3) = TMPVEC(1:3)
                              PN53VEC(1) = PN0VEC(1) + I1
                              PN53VEC(2) = PN0VEC(2) + I2
                              PN53VEC(3) = 0
                           END IF
C
                        END DO
C
                     END DO
                  END DO
               END DO
C
C--------------------------- check whether site is inside defined regime
               IF ( SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
                  IF ( (IQ2.LE.NQ_L) .OR. (IQ2.GE.(NQ_L+NQ_I+1)) )
     &                 PR_INSIDE = .FALSE.
C
               END IF
C
               IF ( PR_INSIDE ) THEN
                  IM = IMQ(IQ2)
                  RCRIT = R(JRCRI(IM),IM)
                  D2VEC = DNRM2(3,PR2VEC,1)*ALAT
                  IF ( D2VEC.GT.RCRIT ) THEN
                     WRITE (*,'(a)') 'RVEC5EXPAND D'
                     WRITE (*,'(a,i8)') '######## iq2 ',IQ2
                     WRITE (*,'(a,5f12.8)') '######## d2vec',D2VEC
                     WRITE (*,'(a,5f12.8)') '######## rcrit',RCRIT
                     CALL STOP_MESSAGE(ROUTINE,'RVEC5EXPAND D')
                  END IF
               END IF
C=======================================================================
C           SYSTEM_DIMENSION(1:2) = '0D'  not implemented yet
C=======================================================================
            ELSE
               CALL STOP_ERROR(ROUTINE,
     &                        'Plot not implemented for 0D calculations'
     &                        )
            END IF
C
C=======================================================================
C                          embedded cluster
C=======================================================================
            IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
               DO IQCLU = 1,NQCLU
                  SAME = .TRUE.
                  DO I5 = 1,3
                     IF ( PN53VEC(I5).NE.N5VEC_QCLU(I5,IQCLU) )
     &                    SAME = .FALSE.
                     IF ( IQ2.NE.IQ_QCLU(IQCLU) ) SAME = .FALSE.
                  END DO
                  IF ( SAME ) THEN
                     IQ2 = NQHOST + IQCLU
                     EXIT
                  END IF
               END DO
            END IF
C=======================================================================
C
C
C***********************************************************************
            IF ( PR_INSIDE ) THEN
C
               NPTS_INSIDE = NPTS_INSIDE + 1
C
C-----------------------------------------------------------------------
C      rotate vector PR2VEC back by MROTR**(-1) to vector PR3VEC
C      within representative atomic cell IQ3 = IQREPQ(IQ2)
C-----------------------------------------------------------------------
C
               IQ3 = IQREPQ(IQ2)
C
               CALL DGEMV('T',3,3,1D0,MROTR(1,1,ISYMGENQ(IQ2)),3,PR2VEC,
     &                    1,0D0,PR3VEC,1)
C
               PR3 = DNRM2(3,PR3VEC,1)
               IF ( ABS(PR3).GT.1D-12 ) THEN
                  PR3HAT(1:3) = PR3VEC(1:3)/PR3
C
                  CALL CALC_RHPLM(PR3HAT(1),PR3HAT(2),PR3HAT(3),RYLM,
     &                            LPLOT_MAX,NLMPLOT_MAX)
C
               ELSE
                  PR3HAT(:) = 0D0
                  RYLM(:) = 0D0
               END IF
C
C------------------------------------------------------------------------
C
               IT3 = ITOQ(1,IQ3)
               IM3 = IMT(IT3)
C
               IR3B = JRCUT(NPAN(IM3),IM3)
               PR3 = PR3*ALAT
C
               DO IR = 2,JRCUT(NPAN(IM3),IM3)
                  IF ( R(IR,IM3).GT.PR3 ) THEN
C               write(*,*) 'EXIT'
                     IR3B = IR
                     EXIT
                  END IF
               END DO
               IF ( IR3B.EQ.0 ) CALL STOP_TRACE_BACK(ROUTINE,'IR3B = 0')
               IR3A = IR3B - 1
C               write(*,'(3F7.4,5i4)')APLOT(IAVEC), BPLOT(IBVEC),
C     &                       PR3,iavec,ibvec,IR3B,JRCUT(NPAN(IM3),IM3)
C
               WR = 1D0/(R(IR3B,IM3)-R(IR3A,IM3))
               W3A = (PR3-R(IR3A,IM3))*WR
               W3B = (R(IR3B,IM3)-PR3)*WR
C
C VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
               IF ( PLOT_V_RHO ) THEN
                  DO IPLOT = 1,NPLOT
                     DO IFP = 1,NFPT(IT3)
                        LM = LMIFP(IFP,IT3)
C
C--------------------------- potential arrays restricted to r > R(JRNS1)
                        IF ( IR3A.GE.JRNS1(IM3) ) THEN
C
                           IF ( IPLOT.EQ.1 ) THEN
                              POT3 = W3A*VNST(IR3A,LM,IT3)
     &                               + W3B*VNST(IR3B,LM,IT3)
                           ELSE
                              POT3 = W3A*BNST(IR3A,LM,IT3)
     &                               + W3B*BNST(IR3B,LM,IT3)
                           END IF
C
                           POTPLOT(IAVEC,IBVEC,IPLOT)
     &                        = POTPLOT(IAVEC,IBVEC,IPLOT)
     &                        + POT3*RYLM(LM)
C
                        END IF
C
                        RHO3 = W3A*RHO2NS(IR3A,LM,IT3,IPLOT)
     &                         /(R(IR3A,IM3)*R(IR3A,IM3))
     &                         + W3B*RHO2NS(IR3B,LM,IT3,IPLOT)
     &                         /(R(IR3B,IM3)*R(IR3B,IM3))
C
                        RHOPLOT(IAVEC,IBVEC,IPLOT)
     &                     = RHOPLOT(IAVEC,IBVEC,IPLOT) + RHO3*RYLM(LM)
C
                     END DO
                  END DO
C VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
               ELSE
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
                  DO IC = 1,3
                     DO LM = 1,LMR_MAX
C                     if( abs(m_lm(lm)).eq. 1 ) then
C
                        JALF3 = W3A*JALF_LMCT(IR3A,LM,IC,IT3)
     &                          + W3B*JALF_LMCT(IR3B,LM,IC,IT3)
C
                        JALF3 = JALF3*PR3*PR3
C
C                  if( ic.eq.1 ) then
C                      jalf3 = xx
C                  else if ( ic.eq.2 ) then
C                      jalf3 = yy
C                  else
C                      jalf3 = 0
C                  end if
C
C
                        JPLOT(IAVEC,IBVEC,IC) = JPLOT(IAVEC,IBVEC,IC)
     &                     + JALF3*RYLM(LM)
C
C                     end if
                     END DO
                  END DO
C
               END IF
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
C***********************************************************************
C
            ELSE IF ( PLOT_V_RHO ) THEN
               POTPLOT(IAVEC,IBVEC,1:NPLOT) = 0D0
               RHOPLOT(IAVEC,IBVEC,1:NPLOT) = 0D0
            ELSE
C
               WRITE (*,*) '###########################'
               JPLOT(IAVEC,IBVEC,1:3) = 0D0
C
            END IF
C                                                              PR_INSIDE
C***********************************************************************
C
            IF ( PR0VEC(3).GT.6.9D0 ) THEN
               IF ( ABS(POTPLOT(IAVEC,IBVEC,1)).GT.1D-7 )
     &               WRITE (*,'(a,i3,2f12.8)') '****',IQ3,PR0VEC(3),
     &              POTPLOT(IAVEC,IBVEC,1)
            END IF
C
         END DO
      END DO
C=======================================================================
C
C
C
C=======================================================================
C                            POTENTIAL
C=======================================================================
      IF ( PLOT2DPR(1) ) THEN
C
         DO IPLOT = 1,NPLOT
C
            LFILNAM = LDATSET0
C
            IF ( IPLOT.EQ.1 ) THEN
               FILNAM = DATSET0(1:LDATSET0)//'_V'
            ELSE
               FILNAM = DATSET0(1:LDATSET0)//'_B'
            END IF
            LFILNAM = LDATSET0 + 2
            CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFILNAM))
C
            DO IAVEC = 0,NAVEC
               DO IBVEC = 0,NBVEC
C
                  WRITE (IOTMP,'(3f25.10)') APLOT(IAVEC),BPLOT(IBVEC),
     &                   POTPLOT(IAVEC,IBVEC,IPLOT)
C
               END DO
               WRITE (IOTMP,*)
            END DO
C
            IF ( IPLOT.EQ.1 ) THEN
               WRITE (6,99002) 'NON-spherical spin-averaged  potential '
     &                         ,FILNAM(1:LFILNAM)
            ELSE
               WRITE (6,99002) 'NON-spherical spin-dependent potential '
     &                         ,FILNAM(1:LFILNAM)
            END IF
C
            CLOSE (IOTMP)
         END DO
C
      END IF
C
C=======================================================================
C                            DENSITIES
C=======================================================================
      IF ( PLOT2DPR(2) ) THEN
C
         DO IPLOT = 1,NPLOT
C
            LFILNAM = LDATSET0
C
            IF ( IPLOT.EQ.1 ) THEN
               FILNAM = DATSET0(1:LDATSET0)//'_RHOCHR'
            ELSE
               FILNAM = DATSET0(1:LDATSET0)//'_RHOSPN'
            END IF
            LFILNAM = LDATSET0 + 7
C
            CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFILNAM))
C
            DO IAVEC = 0,NAVEC
               DO IBVEC = 0,NBVEC
                  WRITE (IOTMP,'(3f25.10)') APLOT(IAVEC),BPLOT(IBVEC),
     &                   RHOPLOT(IAVEC,IBVEC,IPLOT)
C
               END DO
               WRITE (IOTMP,*)
            END DO
C
            IF ( IPLOT.EQ.1 ) THEN
               WRITE (6,99002) 'NON-spherical spin-averaged  densities '
     &                         ,FILNAM(1:LFILNAM)
            ELSE
               WRITE (6,99002) 'NON-spherical spin-dependent densities '
     &                         ,FILNAM(1:LFILNAM)
            END IF
C
            CLOSE (IOTMP)
         END DO
C
      END IF
C
C=======================================================================
C                     create Gnuplot macro file .gpl
C=======================================================================
C
      WRITE (AVECSTR,'(F5.1,F5.1,F5.1)') AVEC(1),AVEC(2),AVEC(3)
      WRITE (BVECSTR,'(F5.1,F5.1,F5.1)') BVEC(1),BVEC(2),BVEC(3)
C
      IF ( PLOT2DPR(1) .OR. PLOT2DPR(2) ) THEN
         WRITE (6,*)
         WRITE (6,*) '     generating gnuplot macro files, ',
     &               'use command(s)'
C
         IF ( PLOT2DPR(1) ) THEN
C------------------------------------------------------plot V potential
C
            CALL GNUPLOT2D(DATSET0(1:LDATSET0)//'_V',DATSET0(1:LDATSET0)
     &                     //'_V.gpl',APLOT(0),APLOT(NAVEC),BPLOT(0),
     &                     BPLOT(NBVEC),'SPR-KKR calculation for '//
     &                     DATSET0(1:LDATSET0),
     &                     'Direction ('//AVECSTR(1:15)//' )',
     &                     'Direction ('//BVECSTR(1:15)//' )',
     &                     'Spin-averaged potential V')
C
            WRITE (6,*) '     gnuplot ',DATSET0(1:LDATSET0)//'_V.gpl'
C
C------------------------------------------------------plot B potential
            IF ( NPLOT.EQ.2 ) THEN
C
               CALL GNUPLOT2D(DATSET0(1:LDATSET0)//'_B',
     &                        DATSET0(1:LDATSET0)//'_B.gpl',APLOT(0),
     &                        APLOT(NAVEC),BPLOT(0),BPLOT(NBVEC),
     &                        'SPR-KKR calculation for '//
     &                        DATSET0(1:LDATSET0),
     &                        'Direction ('//AVECSTR(1:15)//' )',
     &                        'Direction ('//BVECSTR(1:15)//' )',
     &                        'Spin-dependent potential B')
C
               WRITE (6,*) '     gnuplot ',DATSET0(1:LDATSET0)//'_B.gpl'
C
            END IF
C
         END IF
C
C--------------------------------------------------plot charge density
         IF ( PLOT2DPR(2) ) THEN
C
            CALL GNUPLOT2D(DATSET0(1:LDATSET0)//'_RHOCHR',
     &                     DATSET0(1:LDATSET0)//'_RHOCHR.gpl',APLOT(0),
     &                     APLOT(NAVEC),BPLOT(0),BPLOT(NBVEC),
     &                     'SPR-KKR calculation for '//
     &                     DATSET0(1:LDATSET0),
     &                     'Direction ('//AVECSTR(1:15)//' )',
     &                     'Direction ('//BVECSTR(1:15)//' )',
     &                     'Charge density')
C
            WRITE (6,*) '     gnuplot ',DATSET0(1:LDATSET0),
     &                  '_RHOCHR.gpl'
C
C---------------------------------------------------plot spin density
            IF ( NPLOT.EQ.2 ) THEN
C
               CALL GNUPLOT2D(DATSET0(1:LDATSET0)//'_RHOSPN',
     &                        DATSET0(1:LDATSET0)//'_RHOSPN.gpl',
     &                        APLOT(0),APLOT(NAVEC),BPLOT(0),
     &                        BPLOT(NBVEC),'SPR-KKR calculation for '//
     &                        DATSET0(1:LDATSET0),
     &                        'Direction ('//AVECSTR(1:15)//' )',
     &                        'Direction ('//BVECSTR(1:15)//' )',
     &                        'Spin density')
C
               WRITE (6,*) '     gnuplot ',DATSET0(1:LDATSET0),
     &                     '_RHOSPN.gpl'
C
            END IF
C
         END IF
C
         WRITE (6,*) '     to view graphs '
C
      END IF
C
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
      IF ( PLOT_J ) THEN
C
         LFILNAM = LDATSET0
C
         FILNAM = DATSET0(1:LDATSET0)//'_J'
         LFILNAM = LDATSET0 + 2
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFILNAM))
C
         DO IAVEC = 0,NAVEC
            DO IBVEC = 0,NBVEC
C
               WRITE (IOTMP,'(5f25.10)') APLOT(IAVEC),BPLOT(IBVEC),
     &                (JPLOT(IAVEC,IBVEC,IC),IC=1,3)
C
            END DO
            WRITE (IOTMP,*)
         END DO
C
         WRITE (6,99002) 'current density data file',FILNAM(1:LFILNAM)
C
         CLOSE (IOTMP)
C
         CALL GNUPLOT2D_VEC(FILNAM,FILNAM(1:LEN_TRIM(FILNAM))//'.gpl',
     &                      APLOT(0),APLOT(NAVEC),BPLOT(0),BPLOT(NBVEC),
     &                      'SPR-KKR calculation for '//
     &                      DATSET0(1:LDATSET0),
     &                      'Direction ('//AVECSTR(1:15)//' )',
     &                      'Direction ('//BVECSTR(1:15)//' )',
     &                      'current density')
C
      END IF
C JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
      WRITE (6,99003) NPTS,NPTS_INSIDE
C
99001 FORMAT (//,1X,79('*'),/,34X,'<FPPLOT2DPR>',/,1X,79('*'),//,10X,
     &        'creating 2D plots using the plot window parameters ',//,
     &        10X,'->A = (',F6.2,',',F6.2,',',F6.2,' )  ','A =',F7.3,
     &        '   NA =',I4,'  DA =',F6.3,/,10X,'->B = (',F6.2,',',F6.2,
     &        ',',F6.2,' )  ','B =',F7.3,'   NB =',I4,'  DB =',F6.3,/,
     &        10X,'->D = (',F6.2,',',F6.2,',',F6.2,' )     ',//)
99002 FORMAT (/,10X,A,' written to file ',A)
99003 FORMAT (/,10X,'number of plots points     ',I4,/,10X,
     &        'points inside plot regime  ',I4/)
      END
C*==rvec5expand.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE RVEC5EXPAND(A,B,BBINV,B5,N5,C5)
C   ********************************************************************
C   *                                                                  *
C   *   expand A with respect to basis vectors B5_i                    *
C   *                                                                  *
C   *                   A = sum(i) (N5_i + C5_i) * B5_i                *
C   *                                                                  *
C   *                   with  0 <= C5_i < 1  for i=1,2                 *
C   *                                                                  *
C   *   L-bulk:               0 <= C5_4 < 1   C5_3=C5_5=0              *
C   *   I-zone:               0 <= C5_3 < 1   C5_4=C5_5=0              *
C   *   R-bulk:               0 <= C5_5       C5_3=C5_4=0              *
C   *                                                                  *
C   *   for a 2D layered system                                        *
C   *                                                                  *
C   *   A     real vector of dimension 3                               *
C   *   B     real 3x3-matrix containing the basis vectors  B1, B2, B3 *
C   *   BBINV real 3x3-matrix with                                     *
C   *         BBINV = BB^(-1)  and  BB(i,j) = (B_i*B_j)                *
C   *   B5    real 3x5-matrix containing ALL basis vectors             *
C   *   N5    integer vector of dimension 5 with expansion coeff.      *
C   *   C5    real vector of dimension 5 with expansion coefficients   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--RVEC5EXPAND861
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='RVEC5EXPAND')
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      REAL*8 A(3),B(3,3),B5(3,5),BBINV(3,3),C5(5)
      INTEGER N5(5)
C
C Local variables
C
      LOGICAL CHECK
      DOUBLE PRECISION DNRM2
      INTEGER I
      REAL*8 LC,LCA,LCI,LCL,LCR,X(3)
C
C*** End of declarations rewritten by SPAG
C
      CHECK = .TRUE.
C
      C5(1:5) = 0D0
      N5(1:5) = 0
C
C-------------- N5(i) + C5(i) give coefficients w.r.t. B(1:3,i)  i=1,2,3
C
      CALL RVECEXPAND(A,B,BBINV,C5)
C
C........................................................................
      IF ( CHECK ) THEN
         X = 0D0
         DO I = 1,3
            X = X + C5(I)*B(1:3,I)
         END DO
         DO I = 1,3
            IF ( ABS(X(I)-A(I)).GT.1D-8 ) THEN
               WRITE (*,'(a)') 'RVEC5EXPAND A'
               WRITE (*,'(a,3f8.4)') '######## c5 ',C5(1:3)
               WRITE (*,'(a,3f8.4)') '######## a  ',A
               WRITE (*,'(a,3f8.4)') '######## x  ',X
               STOP 'RVEC5EXPAND A'
            END IF
         END DO
         X(1:2) = C5(1:2)
      END IF
C........................................................................
C
      DO I = 1,2
         N5(I) = INT(C5(I)+1000D0) - 1000
         C5(I) = C5(I) - N5(I)
C........................................................................
         IF ( CHECK ) THEN
            IF ( ABS(X(I)-N5(I)-C5(I)).GT.1D-8 .OR. N5(I).LT.0 ) THEN
               WRITE (*,'(a)') 'RVEC5EXPAND B'
               WRITE (*,'(a,3f8.4)') '######## x  ',X(1:2)
               WRITE (*,'(a,3f8.4)') '######## c5 ',C5(1:2)
               WRITE (*,'(a,3i8)') '######## n5 ',N5(1:2)
               STOP 'RVEC5EXPAND B'
            END IF
         END IF
C........................................................................
      END DO
C
C
C
C---------- find remaining coefficients C5(i)  w.r.t. B(1:3,i)  i=3,..,5
C
      LC = DNRM2(3,B5(1,3),1)
      LCL = DNRM2(3,B5(1,4),1)
      LCR = DNRM2(3,B5(1,5),1)
      LCI = LC - LCL - LCR
C
      LCA = C5(3)*LC
C
C----------------------------------------------------------- A in L bulk
      IF ( LCA.LT.LCL ) THEN
         N5(3) = 0
         C5(3) = 0D0
         N5(5) = 0
         C5(5) = 0D0
         N5(4) = INT(LCA/LCL+1000D0) - 1000
         C5(4) = LCA/LCL - N5(4)
         IF ( C5(4).LT.(0D0-TOL) .OR. C5(4).GE.1D0+TOL )
     &         STOP '<RVEC5EXPAND> case A'
C         if(check) write(*,'(a)') 'A in L bulk'
C----------------------------------------------------------- A in I zone
      ELSE IF ( LCA.LT.LCL+LCI ) THEN
         C5(3) = LCA/LC
C         if(check) write(*,'(a)') 'A in I zone'
C----------------------------------------------------------- A in R bulk
      ELSE
         N5(3:4) = 0
         C5(3:4) = 0D0
         IF ( LCA.LT.LC ) THEN
            N5(5) = 0
            C5(5) = LCA/LCR
         ELSE
            N5(5) = INT((LCA-LC)/LCR) + 1
            C5(5) = LCA/LCR - N5(5)
         END IF
         IF ( C5(5).LT.(0D0-TOL) .OR. C5(5).GE.LC/LCR )
     &         STOP '<RVEC5EXPAND> case B'
         IF ( N5(5).LT.0 ) STOP '<RVEC5EXPAND> case C'
C         if(check) write(*,'(a)') 'A in R bulk'
      END IF
C
C........................................................................
      IF ( CHECK ) THEN
         X = 0D0
         DO I = 1,5
            X = X + (N5(I)+C5(I))*B5(1:3,I)
         END DO
         DO I = 1,3
            IF ( ABS(X(I)-A(I)).GT.1D-8 ) THEN
               WRITE (*,'(a)') 'RVEC5EXPAND C'
               WRITE (*,'(a,5f8.4)') '######## c5 ',C5(1:5)
               WRITE (*,'(a,3f8.4)') '######## a  ',A
               WRITE (*,'(a,3f8.4)') '######## x  ',X
               CALL STOP_MESSAGE(ROUTINE,'RVEC5EXPAND C')
            END IF
         END DO
      END IF
C........................................................................
C
      END
C*==gnuplot2d.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GNUPLOT2D(DFILNAM,MFILNAM,GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &                     GPTITLE,XATITLE,YATITLE,ZATITLE)
C **********************************************************************
C * Subroutine to generate gnuplot macro file for 2-dimensional plot   *
C * DFILNAM  LDFILNAM       Name of Data file, and length of string    *
C * MFILNAM  LMFILNAM       Name of macro file, and length of string   *
C * GPXMAX   GPXMIN   and  Y accordingly                               *
C *              lowest and highest X / Y values for exact plot area   *
C *              NOTE: range of Z values is not explicitly set and     *
C *              thus determined by gnuplot                            *
C * GPTITLE LGPTITLE         Title of Plot, and length of string       *
C * X/Y/ZATITLE LX/Y/ZATITLE Title of x/y/z axis, and length of string *
C *                                                                    *
C * MKO dec. 2009                                                      *
C **********************************************************************
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C*--GNUPLOT2D1024
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GNUPLOT2D')
C
C Dummy arguments
C
      CHARACTER*(*) DFILNAM,GPTITLE,MFILNAM,XATITLE,YATITLE,ZATITLE
      DOUBLE PRECISION GPXMAX,GPXMIN,GPYMAX,GPYMIN
C
C*** End of declarations rewritten by SPAG
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,MFILNAM(1:LEN_TRIM(MFILNAM)))
C
      WRITE (IOTMP,'(5A)') 'set contour both'
      WRITE (IOTMP,99001) 'x',GPXMIN,GPXMAX
      WRITE (IOTMP,99001) 'y',GPYMIN,GPYMAX
      WRITE (IOTMP,'(5A)') 'set title "',GPTITLE(1:LEN_TRIM(GPTITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set xlabel "',XATITLE(1:LEN_TRIM(XATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set ylabel "',YATITLE(1:LEN_TRIM(YATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set zlabel "',ZATITLE(1:LEN_TRIM(ZATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') '#set hidden3d'
      WRITE (IOTMP,'(5A)') 'unset key'
      WRITE (IOTMP,'(5A)') 'splot "',DFILNAM(1:LEN_TRIM(DFILNAM)),
     &                     '" u 1:2:3 w l'
      WRITE (IOTMP,'(5A)') '#print "rotate view at will, then "'
      WRITE (IOTMP,'(5A)') '#print "type Ctrl-C to exit, or type',
     &                     ' Enter to generate Postscript"'
      WRITE (IOTMP,'(5A)') 'pause -1'
      WRITE (IOTMP,'(5A)') '#set terminal postscript color'
      WRITE (IOTMP,'(5A)') '#set output "',DFILNAM(1:LEN_TRIM(DFILNAM)),
     &                     '.ps"'
      WRITE (IOTMP,'(5A)') '#replot'
      CLOSE (IOTMP)
C
      WRITE (6,99002) MFILNAM(1:LEN_TRIM(MFILNAM))
C
99001 FORMAT ('set ',A,'range [',E15.5,':',E15.5,']')
99002 FORMAT (/,10X,'gnuplot call: gnuplot  ',A,/)
      END
C*==gnuplot2d_vec.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE GNUPLOT2D_VEC(DFILNAM,MFILNAM,GPXMIN,GPXMAX,GPYMIN,
     &                         GPYMAX,GPTITLE,XATITLE,YATITLE,ZATITLE)
C **********************************************************************
C * Subroutine to generate gnuplot macro file for 2-dimensional plot   *
C * DFILNAM  LDFILNAM       Name of Data file, and length of string    *
C * MFILNAM  LMFILNAM       Name of macro file, and length of string   *
C * GPXMAX   GPXMIN   and  Y accordingly                               *
C *              lowest and highest X / Y values for exact plot area   *
C *              NOTE: range of Z values is not explicitly set and     *
C *              thus determined by gnuplot                            *
C * GPTITLE LGPTITLE         Title of Plot, and length of string       *
C * X/Y/ZATITLE LX/Y/ZATITLE Title of x/y/z axis, and length of string *
C *                                                                    *
C * modified to plot vector field                                      *
C *                                                                    *
C **********************************************************************
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C*--GNUPLOT2D_VEC1102
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GNUPLOT2D_VEC')
C
C Dummy arguments
C
      CHARACTER*(*) DFILNAM,GPTITLE,MFILNAM,XATITLE,YATITLE,ZATITLE
      DOUBLE PRECISION GPXMAX,GPXMIN,GPYMAX,GPYMIN
C
C*** End of declarations rewritten by SPAG
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,MFILNAM(1:LEN_TRIM(MFILNAM)))
C
      WRITE (IOTMP,'(5A)') 'set terminal wxt size 350,262 enhanced font'
     &                     ,' "Verdana,10" persist'
      WRITE (IOTMP,'(5A)') 'set title "',GPTITLE(1:LEN_TRIM(GPTITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set xlabel "',XATITLE(1:LEN_TRIM(XATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set ylabel "',YATITLE(1:LEN_TRIM(YATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'set zlabel "',ZATITLE(1:LEN_TRIM(ZATITLE)),
     &                     '"'
      WRITE (IOTMP,'(5A)') 'unset key'
      WRITE (IOTMP,'(5A)') 'unset tics'
      WRITE (IOTMP,'(5A)') 'unset colorbox'
      WRITE (IOTMP,'(5A)') 'set border 0'
      WRITE (IOTMP,'(5A)') 'set palette defined ( 0 "#ffffff", '//
     &                     '1 "#ffee00", 2 "#ff7000", '//
     &                     '3 "#ee0000", 4 "#7f0000")'
      WRITE (IOTMP,99001) 'x',GPXMIN,GPXMAX
      WRITE (IOTMP,99001) 'y',GPYMIN,GPYMAX
      WRITE (IOTMP,'(5A)') 'set cbrange [0:40]'
      WRITE (IOTMP,'(5A)') '# vector size'
      WRITE (IOTMP,'(5A)') 'h = 0.12'
      WRITE (IOTMP,'(5A)') ''
      WRITE (IOTMP,'(5A)') 'plot "',DFILNAM(1:LEN_TRIM(DFILNAM)),
     &                     '"  u ($1):($2):(h*$3):(h*$4) with vectors'
      WRITE (IOTMP,'(5A)') '#print "type Ctrl-C to exit, or type',
     &                     ' Enter to generate Postscript"'
      WRITE (IOTMP,'(5A)') 'pause -1'
      WRITE (IOTMP,'(5A)') '#set terminal postscript color'
      WRITE (IOTMP,'(5A)') '#set output "',DFILNAM(1:LEN_TRIM(DFILNAM)),
     &                     '.ps"'
      WRITE (IOTMP,'(5A)') '#replot'
      CLOSE (IOTMP)
C
      WRITE (6,99002) MFILNAM(1:LEN_TRIM(MFILNAM))
C
99001 FORMAT ('set ',A,'range [',E15.5,':',E15.5,']')
99002 FORMAT (/,10X,'gnuplot call: gnuplot  ',A,/)
      END
