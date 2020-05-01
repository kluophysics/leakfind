C*==cluscfinitpot.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE CLUSCFINITPOT
C   ********************************************************************
C   *                                                                  *
C   *   prepare the cluster SCF-cycle starting from  scratch           *
C   *                                                                  *
C   *   - set up a guess for the charge and spin density               *
C   *     RHOCHR  and  RHOSPN, respectively, using atomic data         *
C   *     and the Mattheiss construction                               *
C   *   - set up the potential functions  VT and BT                    *
C   *     default:    use the ionicity from the Mattheiss construction *
C   *                 QION may be specified in the input instead       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL,NONMAG
      USE MOD_ANGMOM,ONLY:L_LM,M_LM,NKMMAX
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_RMESH,ONLY:R,NRMAX,FULLPOT,JRWS,R2DRDI,DRDI,NM,NMMAX,RWS,
     &    JRCRI,NSFMAX,JRCUT,JRNSMIN,FLMSF,KLMSF,ISFLM,RMT,NRSFTOT
      USE MOD_TYPES,ONLY:NT,NTMAX,TXT_T,RHOSPN,RHOCHR,CONC,Z,NAT,IMT,
     &    QEL,NTCLU,NTHOST,VMTZ,VT,BT,ITBOT,ITTOP,BEXT,VNST,BNST,LTXT_T,
     &    LMIFP,NFPT,NLMFPT,NLMFPMAX,NCPLWFMAX,NBLK,VAMEG,VAMEF,NLT,
     &    KLMFP,ISOLIKM,NSOLBLK,IKMSOLBLK,CMNTT
      USE MOD_SITES,ONLY:NQCLU,NOQ,IMQ,NQHOST,IQ_QCLU,ITOQ,NQMAX,IQAT,
     &    VLMMAD_BACK,AVMAD,CMNTQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUSCFINITPOT')
      LOGICAL USE_HOST_POTENTIAL
      PARAMETER (USE_HOST_POTENTIAL=.TRUE.)
      REAL*8 DPAS
      PARAMETER (DPAS=0.05D0)
      INTEGER NRAT
      PARAMETER (NRAT=251)
C
C Local variables
C
      REAL*8 AUX,BTOP_TMP,DELQ,DQQ(:),DQT,DRDI_TMP(:,:),DVCT,DX_TMP,
     &       ECORTAB(:,:),ETOT,EXC(:),EXPDX,FLTMP(:),MSPIN(:),QEL0(:),
     &       QES,QION(:),QION0,QIONINP(:),QMAG,QTEST,R1AT(:),R1_TMP,
     &       R2DRDI_TMP(:,:),RAT(:,:),RHO4PI(:,:),RHOAT(:,:),RHOSAT(:,:)
     &       ,RHOVAT(:,:),RINT(:),RTMP(:),RTOP_TMP(:),R_TMP(:,:),
     &       SUMAT(:,:),SUMMSPIN,SUMQEL,SUMQION,SUMQIONINP,SUMQIONMC,
     &       SUMZ,U1(:),U2(:),V00MAD,VAUX(:,:),VMAD(:),VR2AT(:,:),
     &       VTOP_TMP,W00_TMP(:,:),WA(:),WB(:),WC(:),WEXC(:,:),WRAT(:,:)
     &       ,W_RADINT_TMP(:,:),X_EXT,X_HOST,ZZ(:)
      REAL*8 DDOT,TABESTMOM,YLAG
      LOGICAL ESPRESENT,NOSSITER,USEIONMATT,USEQION
      INTEGER I,IA_ERR,IFLAG,IFP,IM,IMEXT,IMHOST,IO,IQ,IQCLU,IQEXT,
     &        IQHOST,IR,IRCUT_TMP(0:1),IRMTIN,IRSF,IRTOP,IRTOP_TMP(:),
     &        ISF,IT,ITEXT,ITHOST,ITMP,JQ,JRAT(:),LM,MC,NPAN_TMP,Z_EXT,
     &        Z_HOST
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE VR2AT,ECORTAB,QIONINP,R1AT,MSPIN,QION
      ALLOCATABLE RINT,JRAT,ZZ,SUMAT,QEL0,RAT,WRAT
      ALLOCATABLE WA,WB,WC,RHOAT,RHOSAT,RHOVAT,W_RADINT_TMP
      ALLOCATABLE R_TMP,R2DRDI_TMP,DRDI_TMP,IRTOP_TMP,RTOP_TMP
      ALLOCATABLE VMAD,DQQ,WEXC,EXC,VAUX,RHO4PI,U1,U2,W00_TMP,FLTMP,RTMP
C
      ALLOCATE (R1AT(NTMAX),ECORTAB(120,NTMAX),QIONINP(NTMAX))
      ALLOCATE (VR2AT(NRAT,NTMAX),MSPIN(NTMAX),QION(NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: VR2AT')
      ALLOCATE (JRAT(NTMAX),QEL0(NTMAX),RINT(NRMAX))
      ALLOCATE (RHOVAT(NRAT,NTMAX),WA(NRMAX),ZZ(NTMAX))
      ALLOCATE (RHOAT(NRAT,NTMAX),WC(NRMAX),WB(NRMAX))
      ALLOCATE (RHOSAT(NRAT,NTMAX),SUMAT(NRAT,NTMAX))
      ALLOCATE (RAT(NRAT,NTMAX),WRAT(NRAT,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RAT')
C
      R1AT(1:NTMAX) = 999999D0
      VR2AT(1:NRAT,1:NTMAX) = 999999D0
C
      ALLOCATE (R_TMP(NRMAX,NMMAX),R2DRDI_TMP(NRMAX,NMMAX))
      ALLOCATE (W_RADINT_TMP(NRMAX,NMMAX))
      ALLOCATE (DRDI_TMP(NRMAX,NMMAX),IRTOP_TMP(NMMAX),RTOP_TMP(NMMAX))
      ALLOCATE (W00_TMP(NRMAX,NMMAX))
      ALLOCATE (FLTMP(NRMAX),RTMP(NRMAX))
C
C
C=======================================================================
C   FULLPOT:      create a temporary  ASA  logarithmic mesh
C=======================================================================
      IF ( FULLPOT ) THEN
C
         DO IM = 1,NM
C
            IRTOP_TMP(IM) = JRCRI(IM)
            RTOP_TMP(IM) = R(JRCRI(IM),IM)
            R1_TMP = R(1,IM)
            DX_TMP = LOG(RTOP_TMP(IM)/R1_TMP)/DBLE(IRTOP_TMP(IM)-1)
            EXPDX = EXP(DX_TMP)
C
            R_TMP(1,IM) = R1_TMP
            DO IR = 2,NRMAX
               R_TMP(IR,IM) = R_TMP(IR-1,IM)*EXPDX
            END DO
            DO IR = 1,NRMAX
               DRDI_TMP(IR,IM) = DX_TMP*R_TMP(IR,IM)
               R2DRDI_TMP(IR,IM) = DRDI_TMP(IR,IM)*R_TMP(IR,IM)**2
            END DO
C
            ITMP = 0
            DO IRSF = 1,NRSFTOT(IM)
               IR = JRCUT(1,IM) + IRSF
               IFLAG = 0
               IF ( IRSF.EQ.NRSFTOT(IM) ) THEN
                  IFLAG = 1
               ELSE
                  IF ( ABS(R(IR,IM)-R(IR+1,IM)).GT.1D-12 ) IFLAG = 1
               END IF
               IF ( IFLAG.EQ.1 ) THEN
                  ITMP = ITMP + 1
                  RTMP(ITMP) = R(IR,IM)
                  FLTMP(ITMP) = FLMSF(IRSF,1,IM)
               END IF
            END DO
C
            DO IR = 1,NRMAX
               IF ( R_TMP(IR,IM).LE.RMT(IM) .OR. R_TMP(IR,IM)
     &              .GT.RTOP_TMP(IM) ) THEN
                  W00_TMP(IR,IM) = 1D0
               ELSE
                  W00_TMP(IR,IM) = YLAG(R_TMP(IR,IM),RTMP,FLTMP,0,3,
     &                             ITMP)/SQRT_4PI
               END IF
            END DO
C
            IRCUT_TMP(0) = 0
            IRCUT_TMP(1) = JRCRI(IM)
            NPAN_TMP = 1
            CALL GET_INTEGRATION_WEIGHTS(NPAN_TMP,IRCUT_TMP(0),JRCRI(IM)
     &         ,W_RADINT_TMP(1,IM))
C
         END DO
C
      ELSE
C
C=======================================================================
C  ASA:              use the current ASA radial mesh
C=======================================================================
C
         R_TMP(1:NRMAX,1:NM) = R(1:NRMAX,1:NM)
         R2DRDI_TMP(1:NRMAX,1:NM) = R2DRDI(1:NRMAX,1:NM)
         DRDI_TMP(1:NRMAX,1:NM) = DRDI(1:NRMAX,1:NM)
         W00_TMP(1:NRMAX,1:NM) = 1D0
         IRTOP_TMP(1:NM) = JRWS(1:NM)
         RTOP_TMP(1:NM) = RWS(1:NM)
C
         DO IM = 1,NM
            IRCUT_TMP(0) = 0
            IRCUT_TMP(1) = JRWS(IM)
            NPAN_TMP = 1
            CALL GET_INTEGRATION_WEIGHTS(NPAN_TMP,IRCUT_TMP(0),JRWS(IM),
     &         W_RADINT_TMP(1,IM))
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C
      CALL RINIT(NRAT*NTMAX,RHOAT)
      CALL RINIT(NRAT*NTMAX,RHOSAT)
      CALL RINIT(NRAT*NTMAX,RHOVAT)
C
      WRITE (6,99013)
C
      CALL RINIT(120*NTMAX,ECORTAB)
C
      DO IT = 1,NT
         QIONINP(IT) = 0D0
         MSPIN(IT) = TABESTMOM(Z(IT))
      END DO
C
      USEQION = .FALSE.
      USEIONMATT = .TRUE.
      NOSSITER = .FALSE.
C
      CALL INPUT_FIND_SECTION('SCF',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_REAL_ARRAY('QION',QIONINP,N_FOUND,NTMAX,0,
     &                               9999D0,0)
         USEQION = FOUND_REAL_ARRAY
         IF ( USEQION ) USEIONMATT = .FALSE.
         CALL SECTION_SET_REAL_ARRAY('MSPIN',MSPIN,N_FOUND,NTMAX,0,
     &                               9999D0,0)
         CALL SECTION_FIND_KEYWORD('NOSSITER',NOSSITER)
      END IF
      IF ( NONMAG ) MSPIN(1:NTMAX) = 0D0
      IF ( NT.EQ.1 ) NOSSITER = .TRUE.
C
      ESPRESENT = .FALSE.
C
C ======================================================================
C               calculate the starting charge and spin density
C ======================================================================
C
      DO IT = 1,NT
         ZZ(IT) = DBLE(Z(IT))
      END DO
C
      CALL SCF0ATOM(IPRINT,NT,ZZ,R1AT,TXT_T,RHOAT,RHOSAT,RHOVAT,VR2AT,
     &              NRAT)
C
      CMNTT(1:NLMFPMAX,1:NT) = 0D0
C
C---------------------------------------------------- set up atomic mesh
C
      EXPDX = EXP(DPAS)
      DO IT = 1,NT
         RAT(1,IT) = R1AT(IT)
         DO IR = 2,NRAT
            RAT(IR,IT) = RAT(IR-1,IT)*EXPDX
         END DO
C
         DO IR = 1,NRAT
            WRAT(IR,IT) = DPAS*RAT(IR,IT)
         END DO
      END DO
C
      IF ( IPRINT.GT.0 ) THEN
         DO IT = 1,NT
            IF ( IT.EQ.1 ) WRITE (6,*)
            WRITE (6,99001) IT,DDOT(NRAT,RHOAT(1,IT),1,WRAT(1,IT),1),
     &                      DDOT(NRAT,RHOSAT(1,IT),1,WRAT(1,IT),1)
         END DO
      END IF
C
C=======================================================================
C
      CALL CLUSCF0MATT(DPAS,SUMAT,RHOAT,R1AT,JRAT,NRAT,WA,RTOP_TMP)
C
      IF ( IPRINT.GT.0 ) THEN
         DO IT = 1,NT
            IF ( IT.EQ.1 ) WRITE (6,*)
            WRITE (6,99002) IT,DDOT(NRAT,SUMAT(1,IT),1,WRAT(1,IT),1)
         END DO
      END IF
C
      SUMQION = 0D0
      SUMQIONINP = 0D0
      SUMQIONMC = 0D0
      SUMMSPIN = 0D0
      QES = 0D0
C
C=======================================================================
C                                           loop over cluster atom types
      DO IT = NTHOST + 1,NTHOST + NTCLU
         IM = IMT(IT)
         IRTOP = IRTOP_TMP(IM)
C
C---------------------------------------------- renormalize spin density
C
         DO IR = 1,IRTOP
            WB(IR) = YLAG(R_TMP(IR,IM),RAT(1,IT),RHOSAT(1,IT),0,3,
     &               JRAT(IT))*DRDI_TMP(IR,IM)*W00_TMP(IR,IM)
         END DO
C
         AUX = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
         IF ( AUX.LE.0.0D0 ) THEN
            AUX = 1D0
         ELSE
            AUX = MSPIN(IT)/AUX
         END IF
C
         DO IR = 1,IRTOP
            WB(IR) = AUX*WB(IR)
            RHOSPN(IR,IT) = WB(IR)/R2DRDI_TMP(IR,IM)
         END DO
C
         QMAG = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
C-----------------------------------------------------------------------
C                use RHOVAT to set up magnetisation density
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG).LT.1D-6 .AND. ABS(MSPIN(IT)).GT.1D-6 ) THEN
C
            DO IR = 1,IRTOP
               WB(IR) = YLAG(R_TMP(IR,IM),RAT(1,IT),RHOVAT(1,IT),0,3,
     &                  JRAT(IT))*DRDI_TMP(IR,IM)*W00_TMP(IR,IM)
            END DO
C
            AUX = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
            IF ( AUX.LE.0.0D0 ) THEN
               AUX = 1D0
            ELSE
               AUX = MSPIN(IT)/AUX
            END IF
C
            DO IR = 1,IRTOP
               WB(IR) = AUX*WB(IR)
               RHOSPN(IR,IT) = WB(IR)/R2DRDI_TMP(IR,IM)
            END DO
C
            QMAG = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
         END IF
C-----------------------------------------------------------------------
C
         IF ( ABS(QMAG-MSPIN(IT)).GT.1D-6 )
     &         CALL STOP_MESSAGE(ROUTINE,'|QMAG-MSPIN(IT)| > 1D-6')
C
C---------------------------------------------- deal with charge density
C
         DO IR = 1,IRTOP
            WB(IR) = YLAG(R_TMP(IR,IM),RAT(1,IT),SUMAT(1,IT),0,3,
     &               JRAT(IT))*DRDI_TMP(IR,IM)*W00_TMP(IR,IM)
            RHOCHR(IR,IT) = WB(IR)/R2DRDI_TMP(IR,IM)
         END DO
C
         QEL0(IT) = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
         IF ( Z(IT).LE.0 ) THEN
            QES = QES + QEL0(IT)*CONC(IT)*NAT(IT)
            ESPRESENT = .TRUE.
         END IF
C
         IF ( USEIONMATT ) THEN
            QION(IT) = ZZ(IT) - QEL0(IT)
         ELSE
            QION(IT) = QIONINP(IT)
         END IF
C
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
         SUMQIONINP = SUMQIONINP + QIONINP(IT)*CONC(IT)*NAT(IT)
         SUMQIONMC = SUMQIONMC + (ZZ(IT)-QEL0(IT))*CONC(IT)*NAT(IT)
         SUMMSPIN = SUMMSPIN + MSPIN(IT)*CONC(IT)*NAT(IT)
C
         IF ( IPRINT.GT.0 ) THEN
            IF ( IT.EQ.NTHOST+1 ) WRITE (6,*)
            WRITE (6,99003) IT,QEL0(IT),QION(IT)
         END IF
C
      END DO
C                                           loop over cluster atom types
C=======================================================================
C
      DELQ = SUMQION/DBLE(NQCLU)
      SUMQION = 0.0D0
      DO IT = NTHOST + 1,NTHOST + NTCLU
         QION(IT) = QION(IT) - DELQ
         SUMQION = SUMQION + QION(IT)*CONC(IT)*NAT(IT)
      END DO
C
      IF ( USEQION ) THEN
         WRITE (6,99015) ' Q (inp)     '
      ELSE
         WRITE (6,99015) ' '
      END IF
      DO IT = NTHOST + 1,NTHOST + NTCLU
         IF ( USEQION ) THEN
            WRITE (6,99016) IT,TXT_T(IT),QIONINP(IT),ZZ(IT) - QEL0(IT),
     &                      QION(IT),MSPIN(IT)
         ELSE
            WRITE (6,99016) IT,TXT_T(IT),ZZ(IT) - QEL0(IT),QION(IT),
     &                      MSPIN(IT)
         END IF
      END DO
      IF ( USEQION ) THEN
         WRITE (6,99017) SUMQIONINP,SUMQIONMC,SUMQION,SUMMSPIN
         IF ( ESPRESENT ) WRITE (6,99014) '             ',QES
      ELSE
         WRITE (6,99017) SUMQIONMC,SUMQION,SUMMSPIN
         IF ( ESPRESENT ) WRITE (6,99014) ' ',QES
      END IF
C
C-----------  Renormalization of the charge density according to SUMQION
C
      SUMQEL = 0D0
      SUMZ = 0D0
C
      DO IT = NTHOST + 1,NTHOST + NTCLU
         IM = IMT(IT)
         IRTOP = IRTOP_TMP(IM)
C
         QION0 = Z(IT) - QEL0(IT)
C
         IF ( Z(IT).GT.0 ) THEN
C----------- normalize valence orbital density and add to charge density
            DO IR = 1,IRTOP
               AUX = YLAG(R_TMP(IR,IM),RAT(1,IT),RHOVAT(1,IT),0,3,
     &               JRAT(IT))*W00_TMP(IR,IM)
               WC(IR) = AUX/R_TMP(IR,IM)**2
               WB(IR) = AUX*DRDI_TMP(IR,IM)
            END DO
C
            AUX = DDOT(IRTOP,WB,1,W_RADINT_TMP(1,IM),1)
C
            AUX = (QION(IT)-QION0)/AUX
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = RHOCHR(IR,IT) - AUX*WC(IR)
            END DO
         ELSE
C----------------------------- normalize charge density for empty sphere
            AUX = QION(IT)/QION0
            DO IR = 1,IRTOP
               RHOCHR(IR,IT) = AUX*RHOCHR(IR,IT)
            END DO
         END IF
C
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*R2DRDI_TMP(IR,IM)
         END DO
C
         QEL(IT) = DDOT(IRTOP,RINT,1,W_RADINT_TMP(1,IM),1)
C
         SUMQEL = SUMQEL + QEL(IT)*CONC(IT)*NAT(IT)
         SUMZ = SUMZ + ZZ(IT)*CONC(IT)*NAT(IT)
C
      END DO
C
      IF ( ABS(SUMQEL-SUMZ).GT.1D-5 ) THEN
         WRITE (6,99018) SUMQEL,SUMZ
         CALL STOP_MESSAGE(ROUTINE,'|SUMQEL-SUMZ| > 1D-5')
      END IF
C
      IF ( USEIONMATT ) THEN
         WRITE (6,99011)
      ELSE
         WRITE (6,99012)
      END IF
C
C=======================================================================
C         calculate guess   ASA   potential and muffintinize
C=======================================================================
C
      IF ( FULLPOT ) THEN
C
         ALLOCATE (VMAD(NQMAX),DQQ(NQMAX),WEXC(NRMAX,2),EXC(NRMAX))
         ALLOCATE (RHO4PI(NRMAX,2),VAUX(NRMAX,2),U1(NRMAX),U2(NRMAX))
C
         VMAD(1:NQMAX) = 999999D0
         DQQ(1:NQMAX) = 999999D0
C
C-----------------------------------------------------------------------
C               net charge for  ALL cluster sites IQ
C-----------------------------------------------------------------------
C
         DO IQ = NQHOST + 1,NQHOST + NQCLU
            DQQ(IQ) = 0.0D0
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               DQQ(IQ) = DQQ(IQ) + CONC(IT)*(QEL(IT)-DBLE(Z(IT)))
            END DO
            CMNTQ(1,IQ) = DQQ(IQ)/SQRT_4PI
         END DO
C
C-----------------------------------------------------------------------
C              Madelung contribution to potential
C-----------------------------------------------------------------------
C
         DO IQ = NQHOST + 1,NQHOST + NQCLU
            V00MAD = 0D0
            DO JQ = NQHOST + 1,NQHOST + NQCLU
               V00MAD = V00MAD + AVMAD(IQ,JQ,1,1)*CMNTQ(1,JQ)
            END DO
            VMAD(IQ) = VLMMAD_BACK(1,IQ)/SQRT_4PI + V00MAD/SQRT_4PI
         END DO
C
C-----------------------------------------------------------------------
C                      loop over types IT
C-----------------------------------------------------------------------
         IF ( IPRINT.GE.0 ) WRITE (6,99004)
C
         DO IT = ITBOT,ITTOP
            IQ = IQAT(1,IT)
            IM = IMT(IT)
            IRTOP = IRTOP_TMP(IM)
C
            DQT = QEL(IT) - DBLE(Z(IT))
C
            CALL SCFPOISSON(RHOCHR(1,IT),VAUX,Z(IT),DRDI_TMP(1,IM),
     &                      R_TMP(1,IM),IRTOP)
C
            QTEST = VAUX(IRTOP,1)*R_TMP(IRTOP,IM)/2.0D0
            IF ( ABS(QTEST-DQT).GT.1D-5 ) WRITE (6,99010) IT,DQT,QTEST
C
            DVCT = -VAUX(IRTOP,1) + 2.0D0*DQT/RTOP_TMP(IM)
C
            IF ( IPRINT.GE.0 ) WRITE (6,99009) IT,IRTOP,Z(IT),QEL(IT),
     &                                DQT,DVCT,IQ,DQQ(IQ),VMAD(IQ)
C
            DO I = 1,IRTOP
               VAUX(I,1) = VAUX(I,1) + DVCT + VMAD(IQ)
               VAUX(I,2) = VAUX(I,1)
               RHO4PI(I,1) = RHOCHR(I,IT)
               RHO4PI(I,2) = RHOSPN(I,IT)
            END DO
C
C --------------------------------------------- skip or add XC-potential
C
            WEXC(1:IRTOP,1) = 0D0
C
            CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
C ---------------------- convention VAUX(I,IS) IS=1,2 == up, down (AKAI)
C
            DO I = 1,IRTOP
               VT(I,IT) = (VAUX(I,1)+VAUX(I,2))/2.0D0 - VMTZ
               BT(I,IT) = (VAUX(I,1)-VAUX(I,2))/2.0D0 - BEXT
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C                      initialize potentials
C-----------------------------------------------------------------------
C
         DO IT = ITBOT,ITTOP
C
            IM = IMT(IT)
C
            DO IR = 1,NLMFPMAX
               DO I = JRNSMIN,NRMAX
                  VNST(I,IR,IT) = 0.0D0
                  BNST(I,IR,IT) = 0.0D0
               END DO
            END DO
C
C-----------------------------------------------------------------------
C        interpolate potential functions V and B from ASA to FP mesh
C-----------------------------------------------------------------------
C
            DO IR = 1,IRTOP_TMP(IM)
               U1(IR) = VT(IR,IT) + 2D0*Z(IT)/R_TMP(IR,IM)
               U2(IR) = BT(IR,IT)
            END DO
C
            IR = IRTOP_TMP(IM)
            VTOP_TMP = VT(IR,IT)
            BTOP_TMP = BT(IR,IT)
            DO IR = 1,JRCRI(IM)
               IF ( R(IR,IM).LE.RTOP_TMP(IM) ) THEN
                  VT(IR,IT) = YLAG(R(IR,IM),R_TMP(1,IM),U1,0,3,
     &                        IRTOP_TMP(IM)) - 2D0*Z(IT)/R(IR,IM)
                  BT(IR,IT) = YLAG(R(IR,IM),R_TMP(1,IM),U2,0,3,
     &                        IRTOP_TMP(IM))
               ELSE
                  VT(IR,IT) = VTOP_TMP
                  BT(IR,IT) = BTOP_TMP
               END IF
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C        convolute spherical  ASA  potential with shape function
C-----------------------------------------------------------------------
C
         WRITE (6,99008)
         DO IT = ITBOT,ITTOP
C
C-----------------------------------------------------------------------
C            IQ = IQAT(1,IT)
C            IQHOST = IQ_QCLU(IQ-NQHOST)
C            ITHOST = ITOQ(1,IQHOST)
CC
C            NFPT(IT) = NFPT(ITHOST)
C            NLMFPT(IT) = NLMFPT(ITHOST)
C            DO IFP = 1,NFPT(IT)
C               LMIFP(IFP,IT) = LMIFP(IFP,ITHOST)
C            END DO
C-----------------------------------------------------------------------
C
            WRITE (6,99005) ' '
            WRITE (6,99006) 'atom type',IT,TXT_T(IT)(1:LTXT_T(IT))
            WRITE (6,99006) 'NFP ',NFPT(IT),'mesh IM ',IM
            WRITE (6,99007) (LMIFP(IFP,IT),L_LM(LMIFP(IFP,IT)),
     &                      M_LM(LMIFP(IFP,IT)),IFP=1,NFPT(IT))
C
            IRMTIN = JRCUT(1,IM)
C
            DO LM = 2,NLMFPT(IT)
               IF ( KLMSF(LM,IM).EQ.1 ) THEN
                  ISF = ISFLM(LM,IM)
C
                  IF ( ISF.GT.NSFMAX .OR. ISF.LT.1 )
     &                 CALL STOP_MESSAGE(ROUTINE,'ISF out of range !!')
C
                  DO IR = IRMTIN + 1,JRCRI(IM)
                     IRSF = IR - IRMTIN
                     VNST(IR,LM,IT) = VT(IR,IT)*FLMSF(IRSF,ISF,IM)
                     BNST(IR,LM,IT) = BT(IR,IT)*FLMSF(IRSF,ISF,IM)
                  END DO
C
               END IF
            END DO
C
            DO IR = IRMTIN + 1,JRCRI(IM)
               IRSF = IR - IRMTIN
               VT(IR,IT) = VT(IR,IT)*FLMSF(IRSF,1,IM)/SQRT_4PI
               BT(IR,IT) = BT(IR,IT)*FLMSF(IRSF,1,IM)/SQRT_4PI
            END DO
C
            VNST(:,:,IT) = 0D0
            BNST(:,:,IT) = 0D0
         END DO
C
      ELSE
C------------------------------------------------------------------- ASA
C
         CALL SCFNEWPOT(.FALSE.,ETOT,ECORTAB)
C
         CALL VMUFTIN(VMTZ)
C
      END IF
C
C=======================================================================
C         use potential of unperturbed host as guess potential
C=======================================================================
C
      IF ( USE_HOST_POTENTIAL ) THEN
C
         DO IQCLU = 1,NQCLU
C
            IQEXT = NQHOST + IQCLU
            IMEXT = IMQ(IQEXT)
C
            IQHOST = IQ_QCLU(IQCLU)
            IMHOST = IMQ(IQHOST)
C
C--------------------------------------------- check for same occupation
            IF ( NOQ(IQEXT).EQ.NOQ(IQEXT) .AND. IMEXT.EQ.IMHOST ) THEN
C
               DO IO = 1,NOQ(IQEXT)
C
                  ITEXT = ITOQ(IO,IQEXT)
                  Z_EXT = Z(ITEXT)
                  X_EXT = CONC(ITEXT)
C
                  ITHOST = ITOQ(IO,IQHOST)
                  Z_HOST = Z(ITHOST)
                  X_HOST = CONC(ITHOST)
C
                  IF ( Z_EXT.NE.Z_HOST .OR. ABS(X_EXT-X_HOST).GT.1D-6 )
     &                 EXIT
C
                  VT(1:NRMAX,ITEXT) = VT(1:NRMAX,ITHOST)
                  BT(1:NRMAX,ITEXT) = BT(1:NRMAX,ITHOST)
C
C                  IF ( FULLPOT ) THEN
C                     VNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITEXT)
C     &                  = VNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITHOST)
C                     BNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITEXT)
C     &                  = BNST(JRNSMIN:NRMAX,1:NLMFPMAX,ITHOST)
C                  END IF
C
               END DO
C
            END IF
C ----------------------------------------------------------------------
C
         END DO
C
      END IF
C ======================================================================
C
      DEALLOCATE (RINT,JRAT,RHOAT,SUMAT,QEL0)
      DEALLOCATE (RAT,WRAT,WA,WB,WC,RHOSAT,ZZ,RHOVAT)
C
C ======================================================================
C                         ONLY FOR TESTING
C ======================================================================
      IF ( FULLPOT .AND. NT.LT.-100 ) THEN
         WRITE (*,*) '###########################################'
         WRITE (*,*) '###########################################'
         WRITE (*,*) '##    BULK SETTINGS COPIED TO CLUSTER    ##'
         WRITE (*,*) '###########################################'
         WRITE (*,*) '###########################################'
C
         MC = NCPLWFMAX
C
         DO IT = 2,NT
            NLT(IT) = NLT(1)
            NFPT(IT) = NFPT(1)
            NLMFPT(IT) = NLMFPT(1)
            NBLK(IT) = NBLK(1)
            ISOLIKM(1:NKMMAX,IT) = ISOLIKM(1:NKMMAX,1)
            NSOLBLK(1:NKMMAX,IT) = NSOLBLK(1:NKMMAX,1)
            IKMSOLBLK(1:NKMMAX,1:NKMMAX,IT)
     &         = IKMSOLBLK(1:NKMMAX,1:NKMMAX,1)
            KLMFP(1:NLMFPMAX,IT) = KLMFP(1:NLMFPMAX,1)
            LMIFP(1:NLMFPMAX,IT) = LMIFP(1:NLMFPMAX,1)
C
            VT(1:NRMAX,IT) = VT(1:NRMAX,1)
            BT(1:NRMAX,IT) = BT(1:NRMAX,1)
            VNST(JRNSMIN:NRMAX,1:NLMFPMAX,IT)
     &         = VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1)
            BNST(JRNSMIN:NRMAX,1:NLMFPMAX,IT)
     &         = BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1)
            VAMEG(1:MC,1:MC,1:NLMFPMAX,1:2,1:NKMMAX,IT)
     &         = VAMEG(1:MC,1:MC,1:NLMFPMAX,1:2,1:NKMMAX,1)
            IF ( IREL.GE.3 ) VAMEF(1:MC,1:MC,1:NLMFPMAX,1:2,1:NKMMAX,IT)
     &           = VAMEF(1:MC,1:MC,1:NLMFPMAX,1:2,1:NKMMAX,1)
         END DO
C
      END IF
C
C ======================================================================
C
99001 FORMAT (10X,'IT =',I3,3X,'atomic check sums  Q =',F10.6,3X,'M =',
     &        F10.6)
99002 FORMAT (10X,'IT =',I3,3X,'Mattheiss sum      Q =',F10.6)
99003 FORMAT (10X,'IT =',I3,3X,'charges         QEL0 =',F10.6,3X,
     &        'QION =',F10.6)
99004 FORMAT (/,1X,79('*'),/,21X,
     &        'setting up potential from Mattheis charge densities',/,
     &        1X,79('*'),//,5X,
     &        ' IT IRTOP   Z    QEL       DQT       DVCT  ',
     &        '   IQ    DQQ       VMAD  ')
99005 FORMAT (10X,A,I4,:,4X,'JRCUT   ',I5,5X,'R =',F12.8,5X,F12.8)
99006 FORMAT (10X,A,I4,4X,A,I5)
99007 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99008 FORMAT (/,10X,'full potential parameters:')
99009 FORMAT (5X,I3,I6,I4,3F10.5,I5,3F10.5)
99010 FORMAT (//,5X,'WARNING from <CLUSCFINITPOT>',/,5X,
     &        'effective charge   DQT   not consistent for IT=',I3,/,5X,
     &        'QEL=',F12.6,' from <SCFPOISSON>:',F12.6,/,1X,79('!'),/)
99011 FORMAT (/,1X,79('*'),/,10X,'start potential set up using',
     &        ' Mattheiss construction for  RHO',/,20X,
     &        'the corresponding ionicity has been used',/,1X,79('*'),/)
99012 FORMAT (/,1X,79('*'),/,10X,'start potential set up using',
     &        ' Mattheiss construction for  RHO',/,1X,79('*'),/)
99013 FORMAT (//,1X,79('*'),/,32X,'<CLUSCFINITPOT>',/,11X,
     &        'setting up a cluster potential from scratch',/,1X,79('*')
     &        ,/)
99014 FORMAT (/,5X,'charge in empty spheres QES',A,F10.4)
99015 FORMAT (5X,'charge and spin moment from input and Mattheiss ',
     &        'construction',//,27X,A,' Z-Q_el     Q(corr)      m_spin')
99016 FORMAT (5X,'type',I3,' ',A,' :',4F12.4)
99017 FORMAT (5X,'sum ',12x,' :',4F12.4)
99018 FORMAT (/,5X,'number of electrons    ',F12.6,/,5X,
     &        'weighted sum over  Z   ',F12.6,/)
      END
C*==cluscf0matt.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE CLUSCF0MATT(DPAS,SUMYAT,YAT,R1AT,JRAT,NRAT,WA,RTOP_TMP)
C   ********************************************************************
C   *                                                                  *
C   *   use the Mattheis prescription to sum up contributions          *
C   *   of neighboring sites for a central site                        *
C   *   summation may be made for   charge   or   potential            *
C   *   depending on supplied input  YAT                               *
C   *   store results in variable    SUMYAT                            *
C   *                                                                  *
C   *   MOL = .FALSE.  -->  solid state calculation                    *
C   *                       contributions from neighboring unit cells  *
C   *         .TRUE.   -->  cluster calculation                        *
C   *                       NO neighbor contributions                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:MOL
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC,IMT,NTHOST,NTCLU
      USE MOD_SITES,ONLY:NQCLU,NQMAX,QBAS,NQ_L,NQ_R,NOQ,NO_QCLU,ITOQ,
     &    IQ_QCLU,IT_OQCLU,N5VEC_QCLU,NQHOST
      USE MOD_LATTICE,ONLY:ALAT,ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,
     &    ADAINV_R,ADAINV_I,SYSTEM_DIMENSION,SWS
      USE MOD_FILES,ONLY:IPRINT,IOTMP,IDUMMY
      USE MOD_RMESH,ONLY:NRMAX,NMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUSCF0MATT')
      INTEGER NTR,NTRM
      PARAMETER (NTR=4,NTRM=(2*NTR+1)*(2*NTR+1)*(2*NTR+1))
C
C Dummy arguments
C
      REAL*8 DPAS
      INTEGER NRAT
      INTEGER JRAT(NTMAX)
      REAL*8 R1AT(*),RTOP_TMP(NMMAX),SUMYAT(NRAT,*),WA(NRMAX),
     &       YAT(NRAT,NTMAX)
C
C Local variables
C
      LOGICAL DONE,KDONET(NT)
      REAL*8 DQMATT(:),MATTRAD,RAD,RQMATT(:,:)
      INTEGER I,IA_ERR,IO,IPRINTLOC,IQCLU,IQCNTR,IQMATT,IQ_IQMATT(:),IT,
     &        J,JO,JQ,JQCLU,JT,N5VEC_QMATT(:,:),NQMATT,NQMATT_I,
     &        NQMATT_L,NQMATT_R,NSHLMATT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DQMATT,IQ_IQMATT,RQMATT,N5VEC_QMATT
C
      IF ( NTRM.LE.0 ) CALL STOP_MESSAGE(ROUTINE,'NTRM.LE.0')
C
      WRITE (6,99001) SWS/ALAT,SWS
C
C-----------------------------------------------------------------------
C                     add on-site contribution
C-----------------------------------------------------------------------
      DO IT = 1,NTHOST + NTCLU
         JRAT(IT) = MIN(NRAT,NINT(LOG(RTOP_TMP(IMT(IT))/R1AT(IT))/DPAS)
     &              +2)
         DO I = 1,JRAT(IT)
            SUMYAT(I,IT) = YAT(I,IT)
         END DO
         KDONET(IT) = .FALSE.
      END DO
C
C=======================================================================
C
      DO IQCLU = 1,NQCLU
         DONE = .TRUE.
         DO IO = 1,NO_QCLU(IQCLU)
            DONE = DONE .AND. KDONET(IT_OQCLU(IO,IQCLU))
         END DO
C
         IF ( .NOT.DONE ) THEN
C
C-----------------------------------------------------------------------
C                 generate cluster around site IQ
C-----------------------------------------------------------------------
            IQCNTR = IQ_QCLU(IQCLU)
            NSHLMATT = 0
            MATTRAD = 5.0D0*(SWS/ALAT)
C
            IPRINTLOC = 0
C
            CALL CLUSSITES(IOTMP,IPRINTLOC,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,MATTRAD,IQCNTR,NQMATT,NQMATT_L,
     &                     NQMATT_I,NQMATT_R,NSHLMATT,NQHOST,NQ_L,NQ_R,
     &                     NQMAX)
C
            IF ( IQCLU.GT.1 ) DEALLOCATE (RQMATT,DQMATT,IQ_IQMATT,
     &           N5VEC_QMATT)
            ALLOCATE (RQMATT(3,NQMATT),DQMATT(NQMATT))
            ALLOCATE (N5VEC_QMATT(5,NQMATT))
            ALLOCATE (IQ_IQMATT(NQMATT),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQMATTS')
C
            READ (IOTMP) ((RQMATT(J,I),J=1,3),DQMATT(I),IQ_IQMATT(I),
     &                   I=1,NQMATT),(IDUMMY,I=1,NSHLMATT)
            READ (IOTMP) ((N5VEC_QMATT(J,I),J=1,5),I=1,NQMATT)
C
            CLOSE (IOTMP)
            IF ( IPRINT.GE.5 ) THEN
               WRITE (6,'(I5)') IDUMMY
               WRITE (6,'(I5,4F10.6,I5)')
     &                (I,(RQMATT(J,I),J=1,3),DQMATT(I),IQ_IQMATT(I),I=1,
     &                NQMATT)
               WRITE (6,'(5X,5I3)') ((N5VEC_QMATT(J,I),J=1,5),I=1,NQMATT
     &                              )
            END IF
C
            DO IO = 1,NO_QCLU(IQCLU)
               IT = IT_OQCLU(IO,IQCLU)
               KDONET(IT) = .TRUE.
C
C-----------------------------------------------------------------------
C        add contributions from surrounding ---> IQMATT = 2,NQMATT
C-----------------------------------------------------------------------
C
               DO IQMATT = 2,NQMATT
                  JQ = IQ_IQMATT(IQMATT)
                  RAD = DQMATT(IQMATT)*ALAT
C
                  CALL SEARCH_N5QTAB(N5VEC_QMATT(1,IQMATT),JQ,
     &                               N5VEC_QCLU,IQ_QCLU,NQCLU,JQCLU)
C
                  IF ( JQCLU.EQ.0 ) THEN
C----------------------------------------------- site is outside cluster
                     DO JO = 1,NOQ(JQ)
                        JT = ITOQ(JO,JQ)
C
                        CALL SCF0SUMUP(RAD,R1AT(JT),R1AT(IT),DPAS,
     &                                 YAT(1,JT),SUMYAT(1,IT),JRAT(IT),
     &                                 NRAT,CONC(JT),WA,NRMAX)
                     END DO
                  ELSE
C------------------------------------------------ site is inside cluster
                     DO JO = 1,NO_QCLU(JQCLU)
                        JT = IT_OQCLU(JO,JQCLU)
C
                        CALL SCF0SUMUP(RAD,R1AT(JT),R1AT(IT),DPAS,
     &                                 YAT(1,JT),SUMYAT(1,IT),JRAT(IT),
     &                                 NRAT,CONC(JT),WA,NRMAX)
                     END DO
                  END IF
C
               END DO
C
            END DO
C
         END IF
C
      END DO
C
C----------------------------------------------------------------------
99001 FORMAT (/,10X,'average Wigner-Seitz radius ',F12.6,' a  = ',F12.6,
     &        ' a.u.',/)
      END
C*==search_n5qtab.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SEARCH_N5QTAB(N5VEC,JQ,N5VEC_TAB,IQ_TAB,NTAB,JTAB)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JQ,JTAB,NTAB
      INTEGER IQ_TAB(NTAB),N5VEC(5),N5VEC_TAB(5,NTAB)
C
C Local variables
C
      INTEGER I,ITAB
C
C*** End of declarations rewritten by SPAG
C
      JTAB = 0
C
      DO ITAB = 1,NTAB
C
         IF ( JQ.EQ.IQ_TAB(ITAB) ) THEN
            DO I = 1,5
               IF ( N5VEC(I).NE.N5VEC_TAB(I,ITAB) ) GOTO 100
            END DO
C
            JTAB = ITAB
            RETURN
         END IF
C
 100  END DO
      END
