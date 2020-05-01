C*==scfnewpot.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFNEWPOT(CALCETOT,ETOT,ECORTAB)
C   ********************************************************************
C   *                                                                  *
C   *   this routine sets up the new potential in  VT and BT           *
C   *   using the new charge and spin densities  RHOCHR and RHOSPN     *
C   *   all is done using   <SCFPOISON>  to solve the poisson          *
C   *   equation and the XC-routines according to  SCFVXC              *
C   *                                                                  *
C   *   concerning the calculation of the total energy ETOT on the     *
C   *   basis of the ASA, see e.g.: Abriksosov PRB 56, 9319 (1997)     *
C   *                                                                  *
C   *   AVMAD: Madelung matrix using Drittler's definition             *
C   *                                                                  *
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *        sdf parametrization by vbh          VBH                   *
C   *                            by mjw          MJW                   *
C   *                            by vwn          VWN   DEFAULT         *
C   *    Perdew, Burke, Ernzendorfer GGA         PBE                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:IBND
      USE MOD_SYMMETRY,ONLY:NSYM,ISYMGENQ,SYMEULANG,IQREPQ
      USE MOD_SCF,ONLY:SCFVXC,SCFSIM
      USE MOD_RMESH,ONLY:DRDI,R2DRDI,R,RWS,JRWS,NRMAX
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,PI,RY_EV
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,CONC,NAT,NT,IMT,VT,BT,BEXT,RHOCHR,
     &    RHOSPN,QEL,Z,NTMAX,CMNTT,KLMFP,NLMFPMAX,VMTZ,OBS_T
      USE MOD_SITES,ONLY:IQBOT,IQTOP,NQ_L,QBAS,NQHOST,NQCLU,NQ,NQMAX,
     &    IQAT,RNNQ,AVMAD,CMNTQ,NLQMAD,NLMQMAD,VLMMAD_HOST,VLMMAD_BACK,
     &    VMAD2D_A,VMAD2D_B,NLMMAD,MAGROT_Q,QMGAM,QMTET,QMPHI
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE,SUB_SYSTEM
      USE MOD_ENERGY,ONLY:EFERMI,EWORK
      IMPLICIT NONE
C*--SCFNEWPOT39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFNEWPOT')
C
C Dummy arguments
C
      LOGICAL CALCETOT
      REAL*8 ETOT
      REAL*8 ECORTAB(120,NTMAX)
C
C Local variables
C
      REAL*8 AGRDRHO(:),AGRDRHOD(:),AGRDRHOU(:),CMNTISQ(:,:),
     &       CMNTIST(:,:),CMNTMTQ(:,:),CMNTMTT(:,:),DQT,DROT_QLM(:,:,:),
     &       DSYM_RLM(:,:,:),DVCT,DVSIM,D_INV(:,:),D_ONE(:,:),D_WRK(:,:)
     &       ,ECORT(NT),EEXCT(NT),EMADT(NT),ESPT(NT),ET,EVCBT(NT),
     &       EVT(NT),EXC(:),GDGAG(:),GDGAGD(:),GDGAGU(:),LAPRHO(:),
     &       LAPRHOD(:),LAPRHOU(:),QCORT,QNETT(:),QTEST,RHO(:),
     &       RHO4PI(:,:),RHOD(:),RHODN,RHOU(:),RHOUP,RINT(:),V00MAD,
     &       VAUX(:,:),VDN,VMAD(:),VUP,WEXC(:,:),X
      LOGICAL GGA,INITIALIZE,Q1M_OCCURS
      INTEGER I,IC,IFLAG_NEUTRALITY,IINV,IM,IPRINT_CMNT,IQ,IR,IRTOP,
     &        ISYM,ISYMP,IT,JQ,JQBOT,JQTOP,L,LM,M,NSYMH
      SAVE DROT_QLM,DSYM_RLM
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE DSYM_RLM,D_INV,D_ONE,D_WRK,DROT_QLM
      ALLOCATABLE CMNTIST,CMNTISQ,CMNTMTT,CMNTMTQ,QNETT
      ALLOCATABLE RHO4PI,WEXC,EXC,RHO,RHOD,RHOU
      ALLOCATABLE RINT,VAUX,VMAD
      ALLOCATABLE AGRDRHO,AGRDRHOD,AGRDRHOU
      ALLOCATABLE LAPRHO,LAPRHOD,LAPRHOU
      ALLOCATABLE GDGAG,GDGAGD,GDGAGU
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (RHO4PI(NRMAX,2),WEXC(NRMAX,2),EXC(NRMAX))
      ALLOCATE (RHO(NRMAX),RHOD(NRMAX),RHOU(NRMAX))
      ALLOCATE (RINT(NRMAX),VAUX(NRMAX,2),VMAD(NQMAX))
C
C   ********************************************************************
C                       INITIALISATION - START
C   ********************************************************************
C
      IF ( INITIALIZE ) THEN
C
C=======================================================================
C              rotation matrices for REAL spherical harmonics
C=======================================================================
C NOTE: the following section is taken over from <FPSCFNEWPOT>
C
         ALLOCATE (DROT_QLM(NLMQMAD,NLMQMAD,IQBOT:IQTOP))
         ALLOCATE (D_ONE(NLMQMAD,NLMQMAD))
         ALLOCATE (D_WRK(NLMQMAD,NLMQMAD),D_INV(NLMQMAD,NLMQMAD))
C
C-----------------------------------------------------------------------
C                create 1-matrix and matrix for inversion
C-----------------------------------------------------------------------
C
         D_INV(1:NLMQMAD,1:NLMQMAD) = 0D0
         D_ONE(1:NLMQMAD,1:NLMQMAD) = 0D0
         I = 0
         DO L = 0,(NLQMAD-1)
            DO M = 1,2*L + 1
               I = I + 1
               D_INV(I,I) = (-1.0D0)**L
               D_ONE(I,I) = 1D0
            END DO
         END DO
C
         NSYMH = NSYM/2
C
         DO IQ = IQBOT,IQTOP
C
            IF ( IQREPQ(IQ).NE.IQ ) THEN
C
C---- connection of NON-representative to associated representative site
C---------------------------------------------------------- IQ <> IQREPQ
C
               IF ( IQREPQ(IQ).GE.IQ )
     &               CALL STOP_MESSAGE(ROUTINE,'IQREP >= IQ')
               ISYM = ISYMGENQ(IQ)
               ISYMP = ISYM - NSYMH
C
               IF ( ISYM.LE.NSYMH ) THEN
C
                  CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYM),
     &                             SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                             DROT_QLM(1,1,IQ))
C
               ELSE
C
                  CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYMP),
     &                             SYMEULANG(2,ISYMP),SYMEULANG(3,ISYMP)
     &                             ,D_WRK)
C
                  DROT_QLM(1:NLMQMAD,1:NLMQMAD,IQ)
     &               = MATMUL(D_INV(1:NLMQMAD,1:NLMQMAD),
     &               D_WRK(1:NLMQMAD,1:NLMQMAD))
C
               END IF
C--------------------------------------------- representative site IQREP
C----------------------------------------------------------- IQ = IQREPQ
C
            ELSE IF ( MAGROT_Q(IQ) ) THEN
C
               CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,QMPHI(IQ),QMTET(IQ),
     &                          QMGAM(IQ),DROT_QLM(1,1,IQ))
C
            ELSE
C
               DROT_QLM(1:NLMQMAD,1:NLMQMAD,IQ)
     &            = D_ONE(1:NLMQMAD,1:NLMQMAD)
C
            END IF
C---------------------------------------------------------------- IQREPQ
C
         END DO
C
C=======================================================================
C
         ALLOCATE (DSYM_RLM(NLMQMAD,NLMQMAD,NSYM))
C
         DSYM_RLM(:,:,:) = 0.0D0
C
         NSYMH = NSYM/2
C
         IINV = NSYMH + 1
C
         DSYM_RLM(1:NLMQMAD,1:NLMQMAD,1) = D_ONE(1:NLMQMAD,1:NLMQMAD)
C
         DSYM_RLM(1:NLMQMAD,1:NLMQMAD,IINV) = D_INV(1:NLMQMAD,1:NLMQMAD)
C
         DO ISYM = 2,NSYMH
C
            CALL ROTMAT_RYLM(NLQMAD,NLMQMAD,SYMEULANG(1,ISYM),
     &                       SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                       DSYM_RLM(1,1,ISYM))
C
C-----------------------------------------------------------------------
C                          include inversion
C-----------------------------------------------------------------------
C
            ISYMP = NSYMH + ISYM
C
            CALL DGEMM('N','N',NLMQMAD,NLMQMAD,NLMQMAD,1D0,
     &                 DSYM_RLM(1,1,ISYM),NLMQMAD,DSYM_RLM(1,1,IINV),
     &                 NLMQMAD,0D0,DSYM_RLM(1,1,ISYMP),NLMQMAD)
C
         END DO
C
         DEALLOCATE (D_ONE,D_WRK,D_INV)
C
         INITIALIZE = .FALSE.
C
      END IF
C
C **********************************************************************
C                       INITALISATION - END
C **********************************************************************
C
C=======================================================================
C                      GGA - parametrisations
C=======================================================================
      IF ( SCFVXC(1:3).EQ.'PBE' .OR. SCFVXC(1:6).EQ.'EV-GGA' ) THEN
C
         GGA = .TRUE.
C
         ALLOCATE (AGRDRHO(NRMAX),AGRDRHOD(NRMAX),AGRDRHOU(NRMAX))
         ALLOCATE (LAPRHO(NRMAX),LAPRHOD(NRMAX),LAPRHOU(NRMAX))
         ALLOCATE (GDGAG(NRMAX),GDGAGD(NRMAX),GDGAGU(NRMAX))
C
      ELSE
C
         GGA = .FALSE.
C
      END IF
C=======================================================================
C
      VMAD(1:NQMAX) = 999999D0
C
      IFLAG_NEUTRALITY = 0
C
C=======================================================================
C                set charge moments for  ALL  sites IQ
C=======================================================================
C
      ALLOCATE (CMNTIST(NLMFPMAX,NTMAX),CMNTISQ(NLMFPMAX,NQMAX))
      ALLOCATE (CMNTMTT(NLMFPMAX,NTMAX),CMNTMTQ(NLMFPMAX,NQMAX))
      ALLOCATE (QNETT(NTMAX))
C
      IPRINT_CMNT = -1
      Q1M_OCCURS = .FALSE.
      DO IT = ITBOT,ITTOP
         CMNTIST(1:NLMFPMAX,IT) = 0D0
         CMNTMTT(2:NLMFPMAX,IT) = CMNTT(2:NLMFPMAX,IT)
         CMNTMTT(1,IT) = QEL(IT)/SQRT_4PI
         IF ( NLMMAD.GE.4 ) THEN
            DO LM = 2,4
               IF ( ABS(CMNTMTT(LM,IT)).GT.1D-6 ) THEN
                  IPRINT_CMNT = IPRINT
                  Q1M_OCCURS = .TRUE.
                  KLMFP(LM,IT) = 1
               END IF
            END DO
         END IF
      END DO
C
      CALL SCFCMNT(IPRINT_CMNT,Q1M_OCCURS,IFLAG_NEUTRALITY,CMNTIST,
     &             CMNTISQ,CMNTMTT,CMNTMTQ,CMNTQ,QNETT,DROT_QLM,
     &             DSYM_RLM)
C
C------------------------------------- suppress dipole terms for testing
CCCCC CMNTQ(2:NLMFPMAX,1:NQ) = 0D0
C
C-----------------------------------------------------------------------
C                         host calculation
C-----------------------------------------------------------------------
      IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
         VLMMAD_HOST(1:NLMFPMAX,IQBOT:IQTOP) = 0D0
C
C----------------------------------------------------------- bulk system
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' .OR. SUB_SYSTEM(2:6)
     &        .EQ.'-BULK' .OR. SYSTEM_DIMENSION(1:2).EQ.'0D' ) THEN
C
            JQBOT = IQBOT
            JQTOP = IQTOP
C
C------------------------------------------------------------- 2D system
         ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SUB_SYSTEM(1:6)
     &             .EQ.'I-ZONE' ) THEN
C
            JQBOT = 1
            JQTOP = NQ
C
            CALL SCFMAD2D_AB
C
         END IF
C
C-----------------------------------------------------------------------
C                          embedded cluster
C-----------------------------------------------------------------------
C
      ELSE
C
         JQBOT = NQHOST + 1
         JQTOP = NQHOST + NQCLU
C
      END IF
C-----------------------------------------------------------------------
C
C=======================================================================
C       Madelung contribution to potential for  SELECTED   sites IQ
C=======================================================================
C
      DO IQ = IQBOT,IQTOP
         V00MAD = 0D0
         DO JQ = JQBOT,JQTOP
            DO LM = 1,NLMQMAD
               V00MAD = V00MAD + AVMAD(IQ,JQ,1,LM)*CMNTQ(LM,JQ)
            END DO
         END DO
         VMAD(IQ) = V00MAD/SQRT_4PI
      END DO
C
C-----------------------------------------------------------------------
C                         host calculation
C-----------------------------------------------------------------------
      IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SUB_SYSTEM(1:6)
     &        .EQ.'I-ZONE' ) THEN
C
            DO IQ = IQBOT,IQTOP
C
               VMAD(IQ) = VMAD(IQ) + VMAD2D_A*(QBAS(3,IQ)-QBAS(3,NQ_L))
     &                    + VMAD2D_B
            END DO
C
         END IF
C
         VLMMAD_HOST(1,IQBOT:IQTOP) = VMAD(IQBOT:IQTOP)*SQRT_4PI
C
C-----------------------------------------------------------------------
C                          embedded cluster
C-----------------------------------------------------------------------
C
      ELSE
C
         DO IQ = IQBOT,IQTOP
C
            VMAD(IQ) = VLMMAD_BACK(1,IQ)/SQRT_4PI + VMAD(IQ)
C
         END DO
C
      END IF
C
C=======================================================================
C
      IF ( IPRINT.GE.0 ) WRITE (6,99002)
C
      ETOT = 0.0D0
C
C-----------------------------------------------------------------------
C                      loop over types IT
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C              Test whether all densities are non-negative
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( IPRINT.GE.1 ) THEN
            DO IR = 1,IRTOP
               IF ( RHOCHR(IR,IT).LT.0.0D0 )
     &               WRITE (6,'(a,a,2i6,1pe15.7)')
     &               'Negative density in <scfnewpot>:    ',
     &              'IT, IR, RHOCHR(IR,IT) =',IT,IR,RHOCHR(IR,IT)
               IF ( RHOCHR(IR,IT).LT.RHOSPN(IR,IT) )
     &               WRITE (6,'(a,a,2i6,2(1pe20.10))')
     &               'Spin larger than charge in <scfnewpot>:    ',
     &              'IT, IR, RHOCHR(IR,IT), RHOSPN(IR,IT) =',IT,IR,
     &              RHOCHR(IR,IT),RHOSPN(IR,IT)
            END DO
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         DQT = QEL(IT) - DBLE(Z(IT))
C
         EMADT(IT) = DQT*VMAD(IQ)/2.0D0
C
         CALL SCFPOISSON(RHOCHR(1,IT),VAUX,Z(IT),DRDI(1,IM),R(1,IM),
     &                   IRTOP)
C
         QTEST = VAUX(IRTOP,1)*R(IRTOP,IM)/2.0D0
         IF ( ABS(QTEST-DQT).GT.1D-5 ) WRITE (6,99007) IT,DQT,QTEST
C
         DO IR = 1,IRTOP
            RINT(IR) = RHOCHR(IR,IT)*(0.5D0*VAUX(IR,1)-Z(IT)/R(IR,IM))
     &                 *R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,RINT,EVCBT(IT))
C
         DVCT = -VAUX(IRTOP,1) + 2.0D0*DQT/RWS(IM)
C
         IF ( IPRINT.GE.0 ) WRITE (6,99003) IT,IRTOP,Z(IT),QEL(IT),DQT,
     &                             DVCT,IQ,CMNTQ(1,IQ)*SQRT_4PI,VMAD(IQ)
C
         DVSIM = -SCFSIM*2.0D0*DQT/RNNQ(IQ)
C
         DO IR = 1,IRTOP
            VAUX(IR,1) = VAUX(IR,1) + DVCT + VMAD(IQ) + DVSIM
            VAUX(IR,2) = VAUX(IR,1)
            RHO4PI(IR,1) = RHOCHR(IR,IT)
            RHO4PI(IR,2) = RHOSPN(IR,IT)
         END DO
C
C --------------------------------------------- skip or add XC-potential
         DO IR = 1,IRTOP
            WEXC(IR,1) = 0D0
         END DO
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
         IF ( SCFVXC.EQ.'NONE      ' ) THEN
            DO IR = 1,IRTOP
               WEXC(IR,1) = 0D0
            END DO
C
         ELSE IF ( SCFVXC(1:3).EQ.'VBH' ) THEN
C
            CALL EXCVBH(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
         ELSE IF ( SCFVXC(1:3).EQ.'MJW' ) THEN
C
            CALL EXCMJW(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX,1)
C
         ELSE IF ( SCFVXC(1:3).EQ.'VWN' ) THEN
C
            CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
C
C= LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA = LDA
         ELSE IF ( GGA ) THEN
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
C
C-------------------- get spin-resolved densities and remove factor 4 PI
C
            X = 4.0D0*PI
            DO IR = 1,IRTOP
               RHO(IR) = RHOCHR(IR,IT)/X
               RHOU(IR) = 0.5D0*(RHOCHR(IR,IT)+RHOSPN(IR,IT))/X
               RHOD(IR) = 0.5D0*(RHOCHR(IR,IT)-RHOSPN(IR,IT))/X
            END DO
C
C---------- calculate  AGRD: n' = d rho/dr  and  LAP: n'' = d^2 rho/dr^2
C
            CALL CALCDFDR(RHO,AGRDRHO,DRDI(1,IM),IRTOP)
            CALL CALCDFDR(AGRDRHO,LAPRHO,DRDI(1,IM),IRTOP)
            CALL CALCDFDR(RHOU,AGRDRHOU,DRDI(1,IM),IRTOP)
            CALL CALCDFDR(AGRDRHOU,LAPRHOU,DRDI(1,IM),IRTOP)
            CALL CALCDFDR(RHOD,AGRDRHOD,DRDI(1,IM),IRTOP)
            CALL CALCDFDR(AGRDRHOD,LAPRHOD,DRDI(1,IM),IRTOP)
C
            DO IR = 1,IRTOP
C
C------------------------------------------------- check whether n' <= 0
C
               IF ( Z(IT).GT.0 ) THEN
                  IF ( AGRDRHO(IR).GT.0.0 )
     &                  WRITE (6,'(a,f13.9,a,1pe12.5,2(a,i3))')
     &                  'WARNING IN <SCFNEWPOT> D RHO/DR > 0.0 FOR R=',
     &                 R(IR,IM),'   RHO(IR) =',RHO(IR),'   IR=',IR,
     &                 '   IT =',IT
                  IF ( AGRDRHOU(IR).GT.0.0 )
     &                  WRITE (6,'(a,f13.9,a,1pe12.5,2(a,i3))') 
     &                'WARNING IN <SCFNEWPOT> D RHO(UP)/DR > 0.0 FOR R='
     &                ,R(IR,IM),'   RHOU(IR) =',RHOU(IR),'   IR=',IR,
     &                '   IT =',IT
                  IF ( AGRDRHOD(IR).GT.0.0 )
     &                  WRITE (6,'(a,f13.9,a,1pe12.5,2(a,i3))') 
     &                'WARNING IN <SCFNEWPOT> D RHO(DN)/DR > 0.0 FOR R='
     &                ,R(IR,IM),'   RHOD(IR) =',RHOD(IR),'   IR=',IR,
     &                '   IT =',IT
               END IF
C
C ---------------------- GDGAG: grad RHO . grad | grad RHO | = n' (-n" )
C
               GDGAG(IR) = -AGRDRHO(IR)*LAPRHO(IR)
               GDGAGU(IR) = -AGRDRHOU(IR)*LAPRHOU(IR)
               GDGAGD(IR) = -AGRDRHOD(IR)*LAPRHOD(IR)
C
C -------------------------------------- LAP: grad^2 RHO = n" + 2 n' / r
C
               LAPRHO(IR) = LAPRHO(IR) + 2D0*AGRDRHO(IR)/R(IR,IM)
               LAPRHOU(IR) = LAPRHOU(IR) + 2D0*AGRDRHOU(IR)/R(IR,IM)
               LAPRHOD(IR) = LAPRHOD(IR) + 2D0*AGRDRHOD(IR)/R(IR,IM)
C
C ------------------------------------------ AGRD: | grad RHO | = | n' |
C
               AGRDRHO(IR) = ABS(AGRDRHO(IR))
               AGRDRHOU(IR) = ABS(AGRDRHOU(IR))
               AGRDRHOD(IR) = ABS(AGRDRHOD(IR))
C
            END DO
C
            WEXC(1:NRMAX,1:2) = 0D0
C
            IF ( SCFVXC(1:3).EQ.'PBE' .OR. SCFVXC(1:6).EQ.'EV-GGA' )
     &           THEN
C
               CALL EXCPBE(VAUX,EXC,WEXC,IRTOP,NRMAX,RHO,AGRDRHO,LAPRHO,
     &                     GDGAG,RHOU,AGRDRHOU,LAPRHOU,GDGAGU,RHOD,
     &                     AGRDRHOD,LAPRHOD,GDGAGD)
C
            ELSE
               STOP '<SCFNEWPOT>:  no available GGA parametrisation'
            END IF
C
            X = 4.0D0*PI
            WEXC(1:IRTOP,1:2) = X*WEXC(1:IRTOP,1:2)
C
C= GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA = GGA
         ELSE
            WRITE (6,99001) SCFVXC
            STOP
         END IF
C-----------------------------------------------------------------------
C
         DO IR = 1,IRTOP
            RINT(IR) = WEXC(IR,1)*R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,RINT,EEXCT(IT))
C
C ---------------------- convention VAUX(I,IS) IS=1,2 == up, down (AKAI)
         DO IR = 1,IRTOP
C
            RHOUP = 0.5D0*(RHOCHR(IR,IT)+RHOSPN(IR,IT))
            RHODN = 0.5D0*(RHOCHR(IR,IT)-RHOSPN(IR,IT))
C
            VUP = VT(IR,IT) + BT(IR,IT)
            VDN = VT(IR,IT) - BT(IR,IT)
C
            RINT(IR) = (RHOUP*VUP+RHODN*VDN)*R2DRDI(IR,IM)
C
         END DO
C
         CALL RRADINT(IM,RINT,EVT(IT))
C
         DO IR = 1,IRTOP
            VT(IR,IT) = (VAUX(IR,1)+VAUX(IR,2))/2.0D0
            BT(IR,IT) = (VAUX(IR,1)-VAUX(IR,2))/2.0D0 - BEXT
         END DO
C
      END DO
      IF ( .NOT.CALCETOT ) RETURN
C
C=======================================================================
C---------------------------------------------- single particle energies
C
      DO IT = ITBOT,ITTOP
         QCORT = 0.0D0
         ECORT(IT) = 0D0
         DO IC = 1,120
            IF ( ABS(ECORTAB(IC,IT)).LE.1D-6 ) EXIT
            ECORT(IT) = ECORT(IT) + ECORTAB(IC,IT)
            QCORT = QCORT + 1D0
         END DO
         ECORT(IT) = ECORT(IT)
         ESPT(IT) = ECORT(IT) + OBS_T(0,IBND,IT)
      END DO
C
      IF ( IPRINT.GT.0 ) WRITE (6,99004)
C------------------------------------------------------------------- <V>
      DO IT = ITBOT,ITTOP
C
         ET = ESPT(IT) + EVCBT(IT) + EEXCT(IT) - EVT(IT) + EMADT(IT)
C
         ETOT = ETOT + NAT(IT)*CONC(IT)*ET
C
         IF ( IPRINT.GT.0 ) WRITE (6,99005) IT,ESPT(IT),EVCBT(IT),
     &                             EEXCT(IT),EVT(IT),EMADT(IT)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99010)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                            ESPT(IT),EVCBT(IT),EEXCT(IT),EVT(IT),
     &                            EMADT(IT),ETOT
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO
C
C=======================================================================
      WRITE (6,99006) ETOT
C
      IF ( SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
         EWORK = VMAD(IQTOP) - VMTZ - EFERMI
         WRITE (6,99009) EWORK,EWORK*RY_EV
      END IF
C
      IF ( IFLAG_NEUTRALITY.GT.0 .AND. SYSTEM_TYPE(1:16)
     &     .NE.'EMBEDDED-CLUSTER' .AND. SYSTEM_TYPE(1:10)
     &     .NE.'LIV       ' .AND. SYSTEM_TYPE(1:10).NE.'LIR       ' )
     &     THEN
         WRITE (6,99008)
         CALL STOP_MESSAGE(ROUTINE,'IFLAG_NEUTRALITY > 0')
      END IF
C=======================================================================
C
99001 FORMAT (//,1X,79('#'),/,10X,'ERROR in <SCFNEWPOT> ',/,10X,
     &        'SCFVXC = ',A,'  not known ',/,1X,79('#'),/)
99002 FORMAT (/,1X,79('*'),/,21X,
     &        'setting up new potential in <SCFNEWPOT>',/,1X,79('*'),//,
     &        5X,' IT  JTOP   Z    QEL       DQT       DVCT  ',
     &        '   IQ    DQQ       VMAD  ')
99003 FORMAT (5X,I3,I6,I4,3F10.5,I5,3F10.5)
99004 FORMAT (/,1X,'type resolved contributions to the total energy',/,
     &        5X,'IT        ESP           VCB           EXC',
     &        '           V             MAD')
99005 FORMAT (5X,I2,1X,5F14.6)
99006 FORMAT (/,5X,'total energy ',F20.8)
99007 FORMAT (//,5X,'WARNING from <SCFNEWPOT>',/,5X,
     &        'effective charge   DQT   not consistent for IT=',I3,/,5X,
     &        'QEL=',F12.6,' from <SCFPOISSON>:',F12.6,/,1X,79('!'),/)
99008 FORMAT (2(/,1X,79('#')),/,10X,
     &        'inconsitencies detected by <SCFNEWPOT>',/,10X,
     &        'see output for information - execution stopped',
     &        2(/,1X,79('#')),/)
99009 FORMAT (/,5X,'work function',F20.8,' Ry',/,5X,'work function',
     &        F20.8,' eV')
99010 FORMAT ('# BUILDBOT: ',A,':  ESPT EVCBT EEXCT EVT EMADT (+ETOT)',
     &        ' for IT =',I5,/,(1PE22.14))
C
      END
