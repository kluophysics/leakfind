C*==finalstate.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     sumat.f
C
C     contains the subroutines:
C         finalstate     =>
C         -  splee       =>
C         -  coeffh      =>
C         -  ujgs        =>
C         -  sumat       =>
C            -  akmj     =>
C            -  vmxv     =>
C
      SUBROUTINE FINALSTATE(GANZ,LANZ1,NSPIN,NREL,EHIGH,PHI,THETA,N,VPR,
     &                      TEST,NSPLEE,IBLOCH,IREL,CLIGHT,CELM,TYP,
     &                      ISTR,POL0,POL0L,IPOL,NATL,POS,IRUMP,LAYS,
     &                      NOD,NEV,NTOT,BRUSEP,BULKX,NCLM,LAYB,TSE,TSO,
     &                      TSEO,TSOE,BXY,ELOW,EPSX,BK,ADU,ADD,IGP,
     &                      IGPRO,ZMESH,Q1,Q2,Q3,Q4,LAYER1,MAXL,K1,AWR,
     &                      LAYP,IFM,PSI2G,IGVAL,IOD,IED,IOED,IEOD,VIH,
     &                      ZPARU,ZPARD,IEDLAYER,IODLAYER,IEODLAYER,
     &                      IOEDLAYER,JT,NT,SPOLX)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,NPM,MLS,LL,LG,MLQNAT,LH,LB,MHDIM1,
     &    MHDIM2,CZERO
      USE MOD_SPEC_KPARA,ONLY:K0
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BXY,GANZ,IBLOCH,IED,IEOD,IFM,IGP,IOD,IOED,IPOL,IREL,JT,
     &        LANZ1,LAYB,LAYER1,LAYP,LAYS,MAXL,N,NCLM,NREL,NSPIN,NSPLEE,
     &        NT,SPOLX,TYP
      REAL*8 CLIGHT,EPSX,PHI,TEST,THETA,VIH,VPR
      COMPLEX*16 EHIGH,ELOW,K1,Q1,Q2,Q3,Q4
      REAL*8 ADD(3),ADU(3),BRUSEP(3),CELM(36124),POL0(3),POL0L(3),
     &       POS(3,NATLM,LAYSM),ZMESH(6),ZPARD(3),ZPARU(3)
      COMPLEX*16 AWR(MLS,LL,NATLM,2),BK(LH,2),BULKX(LH,LH,LAYSM),
     &           PSI2G(MHDIM2,LH,4),TSE(MLQNAT,MLQNAT,LAYSM),
     &           TSEO(MLQNAT,MLQNAT,LAYSM),TSO(MLQNAT,MLQNAT,LAYSM),
     &           TSOE(MLQNAT,MLQNAT,LAYSM)
      INTEGER IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IGPRO(LH),IGVAL(LH),
     &        IODLAYER(0:LAYSM),IOEDLAYER(0:LAYSM),IRUMP(LAYSM),ISTR(2),
     &        NATL(LAYSM),NEV(LAYSM),NOD(LAYSM),NTOT(LAYSM)
C
C Local variables
C
      COMPLEX*16 C1(:,:),C2(:,:),C3(:,:),C4(:,:),C5(:,:),C6(:,:),CMG(:),
     &           E1,EIGA(:),EIGB(:),EPOSHM(:,:,:),EPOSHP(:,:,:),
     &           EPOSLM(:,:,:),EPOSLP(:,:,:),KGZBLK(:),KGZH(:),PS(:,:),
     &           QH3(:,:,:),QH4(:,:,:),QMM(:,:),QMP(:,:),QPM(:,:),
     &           QPP(:,:),QT1(:,:,:),QT2(:,:,:),QT3(:,:,:),QT4(:,:,:),
     &           TMM(:),TMP(:),TPM(:),TPP(:),U0GM(:),U0GP(:),U1GP(:),
     &           UGSM(:,:),UGSP(:,:),WA(:,:),WB(:,:),WMHF(:,:),WMHG(:,:)
     &           ,WMSTF(:,:),WMSTG(:,:),WPHF(:,:),WPHG(:,:),XFEE(:,:,:),
     &           XFEO(:,:,:),XFOE(:,:,:),XFOO(:,:,:)
      REAL*8 EZRY,EZRYH(:),GN(:,:),KGT(:,:),RLVS(:),RLVX(:),RLVY(:),
     &       VS(:),VX(:),VY(:),ZAD,ZAU,ZBD,ZBU,ZSD,ZSU
      INTEGER I,IJ,ILAYER,IZBG,J,K,MAXG,NNSK,NNSR,NSK(:),NSR(:),QD,QE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C1,C2,C3,C4,C5,C6,GN,WA,WB
      ALLOCATABLE PS,VS,VX,VY,QH3,QH4,QT1,QT2,QT3,QT4,CMG,KGT,QMM
      ALLOCATABLE QMP,QPM,TMM,QPP,TMP,TPM,TPP,NSK,NSR
      ALLOCATABLE U0GM,U0GP,U1GP,EIGA,EIGB,XFEE,WMHF,KGZH
      ALLOCATABLE WMHG,WPHF,UGSM,UGSP,WPHG,XFEO,XFOE,XFOO
      ALLOCATABLE RLVS,RLVX,RLVY,WMSTF,WMSTG,EZRYH
      ALLOCATABLE KGZBLK,EPOSHM,EPOSHP,EPOSLM,EPOSLP
      ALLOCATE (C1(LH,LH),C2(LH,LH))
      ALLOCATE (C3(LH,LH),C4(LH,LH),C5(LH,LH),C6(LH,LH),GN(LG,2))
      ALLOCATE (WA(LB,LB),WB(LB,LB))
      ALLOCATE (PS(LH,2),VS(NPM),VX(NPM))
      ALLOCATE (VY(NPM),QH3(LH,LH,LAYSM),QH4(LH,LH,LAYSM))
      ALLOCATE (QT1(LH,LH,LAYSM),QT2(LH,LH,LAYSM),QT3(LH,LH,LAYSM))
      ALLOCATE (QT4(LH,LH,LAYSM),CMG(LH),KGT(2,LG),QMM(LH,LH))
      ALLOCATE (QMP(LH,LH),QPM(LH,LH),TMM(LH),QPP(LH,LH),TMP(LH))
      ALLOCATE (TPM(LH),TPP(LH))
      ALLOCATE (NSK(NPM),NSR(NPM))
      ALLOCATE (U0GM(LH),U0GP(LH),U1GP(LH))
      ALLOCATE (EIGA(LB),EIGB(LB))
      ALLOCATE (XFEE(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (WMHF(MHDIM1,LH),KGZH(LH),WMHG(MHDIM1,LH))
      ALLOCATE (WPHF(MHDIM1,LH),UGSM(LH,LL),UGSP(LH,LL))
      ALLOCATE (WPHG(MHDIM1,LH))
      ALLOCATE (XFEO(MLQNAT,MLQNAT,LAYSM),XFOE(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (XFOO(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (RLVS(NPM),RLVX(NPM),RLVY(NPM))
      ALLOCATE (WMSTF(6,LH),WMSTG(6,LH),EZRYH(LH),KGZBLK(LH))
      ALLOCATE (EPOSHM(NATLM,LH,LAYSM),EPOSHP(NATLM,LH,LAYSM))
      ALLOCATE (EPOSLM(NATLM,LH,LAYSM),EPOSLP(NATLM,LH,LAYSM))
C
C*** End of declarations rewritten by SPAG
C
C     only for final state
C
      MAXG = LG
      DO I = 1,MLQNAT
         DO J = 1,MLQNAT
            DO K = 1,LAYSM
               XFEE(I,J,K) = CZERO
               XFEO(I,J,K) = CZERO
               XFOE(I,J,K) = CZERO
               XFOO(I,J,K) = CZERO
            END DO
         END DO
      END DO
C
      IF ( IBLOCH.EQ.0 ) TEST = REAL(ELOW)
      CALL SPLEE(MAXG,GANZ,LANZ1,NSPIN,NREL,EHIGH,GN,PHI,THETA,K0,QH3,
     &           QH4,QD,1,N,VPR,TEST,NSPLEE,IBLOCH,IREL,CLIGHT,NT,CELM,
     &           TYP,ISTR,POL0,POL0L,IPOL,NATL,POS,NNSK,NNSR,NSK,NSR,
     &           RLVS,RLVX,RLVY,VS,VX,VY,IRUMP,LAYS,EPOSLM,EPOSLP,
     &           EPOSHM,EPOSHP,NOD,NEV,NTOT,BRUSEP,QPM,QMP,QPP,QMM,QT1,
     &           QT2,QT3,QT4,BULKX,NCLM,KGZBLK,KGT,LAYB,TSE,TSO,TSEO,
     &           TSOE,XFEE,XFOO,XFEO,XFOE,IFM,VIH,EPSX,ZPARU,ZPARD)
C
      IF ( IBLOCH.EQ.0 ) THEN
         QE = 2*QD
         E1 = EHIGH
         CALL BLOCH(QPM,QPP,QMM,QMP,K0(1),K0(2),QD,QE,BRUSEP,E1,THETA,
     &              JT,NT)
      END IF
C
C     end for bandstructure and spleed calculation
C
      IF ( IBLOCH.LE.1 ) RETURN
C
      CALL ESENK(GANZ,QD,VPR,ELOW,EHIGH,1,IREL,CLIGHT,KGZH,EZRYH,KGT,BK)
C
      CALL POTIN(VPR,VIH,1,ZPARU,ZPARD,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
C
      CALL RSPROP(PS,GANZ,QD,VPR,ELOW,EHIGH,1,IREL,CLIGHT,BK,ADU,ADD,
     &            EPSX,TPM,TPP,TMP,TMM,IZBG,EZRY,ZAU,ZSU,ZBU,ZAD,ZSD,
     &            ZBD,IBLOCH)
C
      CALL CWH2(KGZH,GANZ,EZRYH,WMHF,WMHG,WPHF,WPHG,IGP,IGVAL,IGPRO,
     &          WMSTF,WMSTG,ZMESH,IZBG)
C
      QE = 2*IGP
      IF ( IREL.EQ.2 ) QE = 4*IGP
C
      DO I = 1,LH
         DO J = 1,LH
            C1(I,J) = CZERO
            C2(I,J) = CZERO
            C3(I,J) = CZERO
            C4(I,J) = CZERO
            C5(I,J) = CZERO
            C6(I,J) = CZERO
         END DO
      END DO
      CALL COEFFH(EIGA,KGZBLK,QD,U0GP,CMG,WMHF,WMHG,WPHF,WPHG,WA,WB,QE,
     &            IREL,IGP,IGVAL,IGPRO,EIGB,GANZ,TPM,TPP,TMP,TMM,PS,
     &            U1GP,U0GM,CLIGHT,EHIGH,K0,Q1,Q2,Q3,Q4,IPOL,SPOLX,
     &            BULKX,C1,C2,C3,C5,C6,IBLOCH)
C
      IF ( SPOLX.EQ.4 ) THEN
         DO IJ = 1,GANZ
            IF ( IPOL.EQ.2 ) THEN
               U0GP(IJ) = CZERO
               U0GM(IJ) = CZERO
            ELSE IF ( IPOL.EQ.1 ) THEN
               U0GP(IJ+GANZ) = CZERO
               U0GM(IJ+GANZ) = CZERO
            END IF
         END DO
      END IF
C
      IF ( IP.GE.3 ) THEN
         IF ( IREL.EQ.1 ) THEN
            WRITE (NOUT1,99001) 0
         ELSE
            WRITE (NOUT1,99002) 0
         END IF
         WRITE (NOUT1,99003) (U0GP(IJ),IJ=1,QD)
         WRITE (NOUT1,99003) (U0GM(IJ),IJ=1,QD)
      END IF
C
      CALL CPSI2G(CMG,PSI2G,WMHF,WPHF,WMHG,WPHG,IREL,GANZ,EZRYH,QD,
     &            WMSTF,WMSTG,ZMESH,IGP,IGVAL,IGPRO,U0GP,U0GM,KGZBLK,
     &            IPOL)
C
      DO I = 1,LH
         DO J = 1,LH
            C1(I,J) = CZERO
            C2(I,J) = CZERO
            C3(I,J) = CZERO
            C4(I,J) = CZERO
            C5(I,J) = CZERO
         END DO
      END DO
      CALL UJGS(QH3,QH4,QT2,QT4,UGSM,UGSP,LAYER1,QD,LAYS,BULKX,U1GP,
     &          LAYB,C1,C2,C3,C4,C5)
C
      IF ( SPOLX.EQ.4 ) THEN
         DO ILAYER = 1,LAYP
            DO IJ = 1,GANZ
               IF ( IPOL.EQ.2 ) THEN
                  UGSP(IJ,ILAYER) = CZERO
                  UGSM(IJ,ILAYER) = CZERO
               ELSE IF ( IPOL.EQ.1 ) THEN
                  UGSP(IJ+GANZ,ILAYER) = CZERO
                  UGSM(IJ+GANZ,ILAYER) = CZERO
               END IF
            END DO
         END DO
      END IF
C
      IF ( IP.GE.3 ) THEN
         DO ILAYER = 1,LAYP
            IF ( IREL.EQ.1 ) THEN
               WRITE (NOUT1,99001) ILAYER
            ELSE
               WRITE (NOUT1,99002) ILAYER
            END IF
            WRITE (NOUT1,99003) (UGSP(IJ,ILAYER),IJ=1,QD)
            WRITE (NOUT1,99003) (UGSM(IJ,ILAYER),IJ=1,QD)
         END DO
      END IF
C
      IF ( IREL.EQ.2 ) THEN
         K1 = CDSQRT(2.0*ELOW+ELOW*ELOW/(CLIGHT*CLIGHT))
      ELSE
         K1 = CDSQRT(2.0*ELOW)
      END IF
C
      CALL SUMAT(MAXL,IREL,CLIGHT,EHIGH,GANZ,XFEE,XFOO,UGSP,UGSM,LAYER1,
     &           AWR,IPOL,NATL,LAYS,EPOSHM,EPOSHP,XFEO,XFOE,IRUMP,NOD,
     &           NEV,NTOT,IOD,IED,IOED,IEOD,KGT,LAYB,IEDLAYER,IODLAYER,
     &           IEODLAYER,IOEDLAYER,BXY)
C
      RETURN
C
99001 FORMAT (1x,'ujg+ and ujg-        for layer',4x,i3)
99002 FORMAT (1x,'ujgs+(+,-) and ujgs-(+,-) for layer',4x,i3)
99003 FORMAT (5(e12.5,e12.5))
      END
C*==splee.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPLEE(MAXG,GANZ,LANZ,NSPIN,NREL,EEG,GN,PHI,THETA,K0,
     &                 QH3,QH4,QD,LHL,NTEST,VPR,TEST,NSPLEE,IBLOCH,IREL,
     &                 CLIGHT,NT,CELM,TYP,ISTR,POL0,POL0L,IPOL,NATL,POS,
     &                 NNSK,NNSR,NSK,NSR,RLVS,RLVX,RLVY,VS,VX,VY,IRUMP,
     &                 LAYS,EPOSLM,EPOSLP,EPOSHM,EPOSHP,NOD,NEV,NTOT,
     &                 BRUSEP,QPM1,QMP1,QPP1,QMM1,QT1,QT2,QT3,QT4,BULKX,
     &                 NCLM,KGZBLK,KGT,LAYB,TSE,TSO,TSEO,TSOE,XFEE,XFOO,
     &                 XFEO,XFOE,IFM,VIX,EPSX,ZPARU,ZPARD)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQ,NPM,LG,MLQNAT,LH,MLZP,PI,
     &    CZERO
      USE MOD_SPEC_LMKMS,ONLY:KME,KMO,SW,TAUW
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_POSI,ONLY:PXE,PXO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLR,MLT,MLF,MLA,MLB,MAXRFX,NFM,NDLMM,MLRNAT,MLTNAT,MLFNAT
      PARAMETER (MLR=(ML+1)*ML/2,MLT=ML*(ML-1)/2,MLF=2*MLQ,MLA=2*ML-1,
     &           MLB=MLA*MLA,MAXRFX=100,NFM=NATLM*NATLM-NATLM+1,
     &           NDLMM=(2*(ML-1)+1)**2,MLRNAT=MLR*NATLM,
     &           MLTNAT=MLT*NATLM,MLFNAT=MLF*NATLM)
C
C Dummy arguments
C
      REAL*8 CLIGHT,EPSX,PHI,TEST,THETA,VIX,VPR
      COMPLEX*16 EEG
      INTEGER GANZ,IBLOCH,IFM,IPOL,IREL,LANZ,LAYB,LAYS,LHL,MAXG,NCLM,
     &        NNSK,NNSR,NREL,NSPIN,NSPLEE,NT,NTEST,QD,TYP
      REAL*8 BRUSEP(3),CELM(MLZP),GN(LG,2),K0(2),KGT(2,LG),POL0(3),
     &       POL0L(3),POS(3,NATLM,LAYSM),RLVS(NPM),RLVX(NPM),RLVY(NPM),
     &       VS(NPM),VX(NPM),VY(NPM),ZPARD(3),ZPARU(3)
      COMPLEX*16 BULKX(LH,LH,LAYSM),EPOSHM(NATLM,LH,LAYSM),
     &           EPOSHP(NATLM,LH,LAYSM),EPOSLM(NATLM,LH,LAYSM),
     &           EPOSLP(NATLM,LH,LAYSM),KGZBLK(LH),QH3(LH,LH,LAYSM),
     &           QH4(LH,LH,LAYSM),QMM1(LH,LH),QMP1(LH,LH),QPM1(LH,LH),
     &           QPP1(LH,LH),QT1(LH,LH,LAYSM),QT2(LH,LH,LAYSM),
     &           QT3(LH,LH,LAYSM),QT4(LH,LH,LAYSM),
     &           TSE(MLQNAT,MLQNAT,LAYSM),TSEO(MLQNAT,MLQNAT,LAYSM),
     &           TSO(MLQNAT,MLQNAT,LAYSM),TSOE(MLQNAT,MLQNAT,LAYSM),
     &           XFEE(MLQNAT,MLQNAT,LAYSM),XFEO(MLQNAT,MLQNAT,LAYSM),
     &           XFOE(MLQNAT,MLQNAT,LAYSM),XFOO(MLQNAT,MLQNAT,LAYSM)
      INTEGER IRUMP(LAYSM),ISTR(2),NATL(LAYSM),NEV(LAYSM),NOD(LAYSM),
     &        NSK(NPM),NSR(NPM),NTOT(LAYSM)
C
C Local variables
C
      REAL*8 ADD(3),ADU(3),EEGR,EKIN,EMACH,EMACH1,EZRY,EZRYH(LH),FAC1(:)
     &       ,FAC2(:),KGS(:,:,:),PHIS(:),PK(:),PKGM(:),PKGMN(:),PKN(:),
     &       PN(:),POL0G(:,:),THETAS(:),V0I,YLM(:),ZAD,ZAU,ZBD,ZBU,ZSD,
     &       ZSU
      COMPLEX*16 BK(LH,2),BPM(LH,LH),BPP(LH,LH),C1(:,:),C1ONE(:,:),
     &           C2(:,:),C3(:,:),C4(:,:),C5(:,:),C6(:,:),C7(:,:),C8(:,:)
     &           ,CLE(:,:),CLO(:,:),CLXE(:,:,:),CLXEO(:,:,:),CLXO(:,:,:)
     &           ,CLXOE(:,:,:),D1(:,:),D2(:,:),D3(:,:),D4(:,:),DLM1(:,:)
     &           ,GEVEN(:,:),GEVOD(:,:),GODD(:,:),GODEV(:,:),HALM(:,:),
     &           INTFAK(:),INTG(:),INTGPM(:),K03,K0BTRS,KGFAK1(:,:,:),
     &           KGFAK2(:,:,:),KGM(:),KGZH(LH),MCON(:),ME(:,:,:,:),
     &           MMM(LH,LH),MMP(LH,LH),MO(:,:,:,:),MPM(LH,LH),MPP(LH,LH)
     &           ,MSE(:,:,:,:),MSO(:,:,:,:),POLG(:,:),PS(LH,2),QMM2(:,:)
     &           ,QMP2(:,:),QPM2(:,:),QPP2(:,:),SG(:,:,:),SUME(:),
     &           SUMEO(:),SUMO(:),SUMOE(:),TE(:),TMM(LH),TMP(LH),TO(:),
     &           TPM(LH),TPP(LH),VALM(:,:),XALM(:,:),XEEH(:,:),XEOH(:,:)
     &           ,XOEH(:,:),XOOH(:,:)
      INTEGER ESCAPE(:),G,GINT(:,:),GRFX(MAXRFX,2),I,IAN,IG,IL,IM,INAT,
     &        INATS,IOE,IT,IXX,IZBG,J,JHELP(:),JL,KKK,L,LAYER,LXX,N1,N2,
     &        NALM,NAT,PQ(:,:),RFXANZ,SD,XED,XOD
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
      ALLOCATABLE C1,C2,C3,C4,C5,C6,C7,C8,D1,D2,D3,D4,ME,SG,TE,MO,PK
      ALLOCATABLE PN,TO,PQ,CLE,CLO,KGM,MSE,KGS,MSO,PKN,YLM,FAC1,FAC2
      ALLOCATABLE DLM1,QMM2,QMP2,QPM2,QPP2,GODD,HALM,XEEH,CLXE
      ALLOCATABLE CLXO,INTG,MCON,POLG,SUME,SUMO,VALM,XALM,XEOH,XOEH
      ALLOCATABLE XOOH,PHIS,PKGM,GINT,C1ONE,POL0G,JHELP,GEVEN,GEVOD
      ALLOCATABLE GODEV,CLXEO,CLXOE,PKGMN,SUMEO,SUMOE,KGFAK1,KGFAK2
      ALLOCATABLE ESCAPE,INTFAK,THETAS,INTGPM
      ALLOCATE (C1(LH,LH),C2(LH,LH),C3(LH,LH),C4(LH,LH),C5(LH,LH))
      ALLOCATE (C6(LH,LH),C7(LH,LH),C8(LH,LH),D1(MLQNAT,MLQNAT))
      ALLOCATE (D2(MLQNAT,MLQNAT),D3(MLQNAT,MLQNAT))
      ALLOCATE (D4(MLQNAT,MLQNAT),ME(LG,2,MLQNAT,2),SG(LG,2,2))
      ALLOCATE (TE(MLQNAT),MO(LG,2,MLQNAT,2),PK(LG),PN(LG))
      ALLOCATE (TO(MLQNAT),PQ(LG,2),CLE(MLQNAT,2),CLO(MLQNAT,2))
      ALLOCATE (KGM(LG),MSE(LG,2,MLQNAT,2),KGS(LG,2,3))
      ALLOCATE (MSO(LG,2,MLQNAT,2),PKN(LG),YLM(MLB),FAC1(MLA))
      ALLOCATE (FAC2(MLB),DLM1(NDLMM,NFM),QMM2(LH,LH),QMP2(LH,LH))
      ALLOCATE (QPM2(LH,LH),QPP2(LH,LH),GODD(MLTNAT,MLTNAT))
      ALLOCATE (HALM(MLFNAT,MLFNAT))
      ALLOCATE (XEEH(MLQNAT,MLQNAT),CLXE(MLQNAT,MLQNAT,2))
      ALLOCATE (CLXO(MLQNAT,MLQNAT,2),INTG(LG),MCON(LG),POLG(LG,3))
      ALLOCATE (SUME(MLQNAT),SUMO(MLQNAT),VALM(MLFNAT,MLFNAT))
      ALLOCATE (XALM(MLFNAT,MLFNAT),XEOH(MLQNAT,MLQNAT))
      ALLOCATE (XOEH(MLQNAT,MLQNAT),XOOH(MLQNAT,MLQNAT),PHIS(LG))
      ALLOCATE (PKGM(LG),GINT(NPM,2),C1ONE(LH,LH),POL0G(LG,3))
      ALLOCATE (JHELP(LAYSM),GEVEN(MLRNAT,MLRNAT))
      ALLOCATE (GEVOD(MLRNAT,MLTNAT),GODEV(MLTNAT,MLRNAT))
      ALLOCATE (CLXEO(MLQNAT,MLQNAT,2),CLXOE(MLQNAT,MLQNAT,2))
      ALLOCATE (PKGMN(LG),SUMEO(MLQNAT),SUMOE(MLQNAT))
      ALLOCATE (KGFAK1(LG,2,LAYSM),KGFAK2(LG,2,LAYSM),ESCAPE(LG))
      ALLOCATE (INTFAK(LG),THETAS(LG),INTGPM(LG))
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /****************************************************************/
C     # purpose:                                                       *
C                                                                      *
C     # subroutines and functions called from this routine :           *
C       ilmkms    iposi     celmg     ixcon     check1    pqset        *
C       k0set     atphas    streu     gikam     gikamb    xmat         *
C       xpen      xrel      cmatpr    onemx     mmat      tcopy        *
C       scopy     prop      mult      one       sdouble   copy         *
C       reflex    kgtset    sgset     polkmp    splout                 *
C                                                                      *
C     Note:   uncomment call streu to calculate to, te                 *
C     /****************************************************************/
C
C
C
C
C
C
C     *** common blocks ***
C
C
C
C
C
C
      CALL CHECK1(NSPIN,NREL,GANZ,MAXG)
C
      EEGR = DBLE(EEG)
      EMACH = 1.D-9
      EMACH1 = 1.D-9
      RFXANZ = 0
      N1 = 1
      N2 = 2
C
      DO I = 1,MLFNAT
         DO J = 1,MLFNAT
            XALM(I,J) = CZERO
            HALM(I,J) = CZERO
            VALM(I,J) = CZERO
         END DO
      END DO
      DO I = 1,MLTNAT
         DO J = 1,MLTNAT
            GODD(I,J) = CZERO
         END DO
      END DO
      DO I = 1,MLRNAT
         DO J = 1,MLRNAT
            GEVEN(I,J) = CZERO
         END DO
      END DO
      DO I = 1,MLTNAT
         DO J = 1,MLRNAT
            GODEV(I,J) = CZERO
         END DO
      END DO
      DO I = 1,MLRNAT
         DO J = 1,MLTNAT
            GEVOD(I,J) = CZERO
         END DO
      END DO
C
      IF ( IREL.EQ.1 ) THEN
         QD = GANZ
         SD = 1
      ELSE IF ( IREL.EQ.2 ) THEN
         QD = 2*GANZ
         SD = 2
      END IF
C
      V0I = DIMAG(EEG)
C
      CALL PQSET(PQ,GANZ)
C
      DO IT = 1,LAYS
         IF ( IREL.EQ.1 ) THEN
            XED = NEV(IT)
            XOD = NOD(IT)
            NALM = NTOT(IT)
         ELSE IF ( IREL.EQ.2 ) THEN
            XED = NTOT(IT)
            XOD = NTOT(IT)
            NALM = 2*NTOT(IT)
         END IF
C
         CALL ILMKMS()
         CALL IPOSI
         CALL IXCON(CLE,CLXE,KME,SW)
         CALL IXCON(CLO,CLXO,KMO,SW)
         CALL IXCONA(CLE,CLO,CLXEO,KME,SW,KMO)
         CALL IXCONA(CLO,CLE,CLXOE,KMO,SW,KME)
         CALL K0SET(K0,K03,K0BTRS,ME,MO,MSE,MSO,MCON,KGFAK1,GN,GINT,EEG,
     &              THETA,PHI,CLE,CLO,GANZ,NREL,LHL,TEST,XED,XOD,INTFAK,
     &              KGM,CLIGHT,IBLOCH,VX,VY,VS,NNSR,NSR,RLVX,RLVY,RLVS,
     &              NNSK,NSK,KGFAK2,BRUSEP,NEV(IT),NOD(IT),IT,KGZBLK,
     &              KGT)
         IF ( LHL.EQ.1 .AND. IP.GT.1 ) WRITE (NOUT1,99001) K0(1),K0(2)
         IF ( LHL.EQ.1 .AND. NTEST.EQ.1 .AND. IP.GT.1 .AND. IPOL.EQ.1 )
     &        THEN
            IF ( NT.LE.1 ) THEN
               WRITE (NOUT1,99002)
               WRITE (NOUT1,99003) (GINT(I,1),GINT(I,2),I=1,GANZ)
               WRITE (NOUT1,99004) (GN(I,1),GN(I,2),I=1,GANZ)
            END IF
         END IF
C
         CALL ATPHAS(POS,EPOSLM,EPOSLP,EPOSHM,EPOSHP,EEG,GANZ,LHL,IREL,
     &               CLIGHT,IT,NATL,KGT)
         NAT = NATL(IT)
         DO IAN = 1,NAT
            IF ( IP.GE.4 ) WRITE (NOUT1,99005) IT,IAN
            DO IM = 1,IREL
               DO L = 1,GANZ
                  IG = GANZ*(IM-1) + L
                  IF ( LHL.EQ.1 ) THEN
                     IF ( IP.GT.2 ) WRITE (NOUT1,99006)
     &                    EPOSHM(IAN,IG,IT),EPOSHP(IAN,IG,IT)
                  ELSE IF ( LHL.EQ.2 ) THEN
                     IF ( IP.GT.2 ) WRITE (NOUT1,99006)
     &                    EPOSLM(IAN,IG,IT),EPOSLP(IAN,IG,IT)
                  END IF
               END DO
            END DO
         END DO
C
         CALL COPYSTREU(TE,TO,TSE,TSO,IT)
C
C         if irel is equal to 1 then write data from tse/tso to te/to:
C
         IF ( IREL.EQ.1 ) THEN
            DO I = 1,LAYS
               DO J = 1,NATL(I)
                  IF ( J.GT.1 .OR. ML.NE.5 ) THEN
                     WRITE (*,*) 'restrictions apply for this kind ',
     &                           'of calculation. check subroutine ',
     &                           'splee for more detail.'
                     STOP
                  END IF
               END DO
            END DO
            DO I = 1,MLQ
               IF ( I.LE.15 ) THEN
                  TE(16-I) = TSE(I,I,IT)
               ELSE
                  TO(I-15) = TSO(I,I,IT)
               END IF
            END DO
         END IF
C
         IF ( IRUMP(IT).EQ.0 ) THEN
            IF ( REAL(EEG).LT.20.0D0 .OR. IBLOCH.EQ.1 ) THEN
C
               IF ( NAT.EQ.1 ) THEN
                  CALL GIKAM(GODD,GEVEN,0.D0,0.D0,EEGR,V0I,K0(1),K0(2),
     &                       CELM,EMACH1,CLIGHT,NREL,NEV(IT),NOD(IT))
C
               ELSE IF ( NAT.GT.1 ) THEN
                  CALL GIKAMM(IT,EEGR,V0I,CLIGHT,NATL,POS,NNSK,NNSR,NSK,
     &                        NSR,RLVS,RLVX,RLVY,VS,VX,VY,K0(1),K0(2),
     &                        NREL,DLM1)
                  IOE = 0
                  CALL XMAT(GODD,NOD(IT),NOD(IT),IT,IOE,MLT,NATL,DLM1,
     &                      CELM,NCLM)
                  IOE = 1
                  CALL XMAT(GEVEN,NEV(IT),NEV(IT),IT,IOE,MLR,NATL,DLM1,
     &                      CELM,NCLM)
C
               END IF
            END IF
         ELSE IF ( IRUMP(IT).GT.0 ) THEN
            IF ( REAL(EEG).LT.20.0D0 .OR. IBLOCH.EQ.1 ) THEN
C
               CALL GIKAMM(IT,EEGR,V0I,CLIGHT,NATL,POS,NNSK,NNSR,NSK,
     &                     NSR,RLVS,RLVX,RLVY,VS,VX,VY,K0(1),K0(2),NREL,
     &                     DLM1)
               IOE = 0
               CALL XMAT(GODD,NOD(IT),NOD(IT),IT,IOE,MLT,NATL,DLM1,CELM,
     &                   NCLM)
               IOE = 1
               CALL XMAT(GEVEN,NEV(IT),NEV(IT),IT,IOE,MLR,NATL,DLM1,
     &                   CELM,NCLM)
C
               CALL CELMR(CELM,YLM,FAC2,FAC1,2)
C
               CALL XMATR(POS,GEVEN,GODD,EEGR,CELM,V0I,EMACH,K0(1),
     &                    K0(2),0.D0,0.D0,NAT,NTOT(IT),NREL,CLIGHT,IT,
     &                    GODEV,GEVOD,NEV(IT),NOD(IT))
            END IF
         END IF
C
         DO I = 1,MLQNAT
            DO J = 1,MLQNAT
               XEEH(I,J) = CZERO
               XOOH(I,J) = CZERO
               XEOH(I,J) = CZERO
               XOEH(I,J) = CZERO
               D1(I,J) = CZERO
               D2(I,J) = CZERO
               D3(I,J) = CZERO
               D4(I,J) = CZERO
            END DO
         END DO
C
         IF ( IREL.EQ.1 ) THEN
            IF ( IRUMP(IT).GT.0 ) THEN
C
               CALL XPEN(XEEH,TE,GEVEN,NEV(IT))
               CALL XPENR(XOEH,TO,GODEV,NOD(IT),NEV(IT))
               CALL XPEN(XOOH,TO,GODD,NOD(IT))
               CALL XPENR(XEOH,TE,GEVOD,NEV(IT),NOD(IT))
               CALL XCOPR(XALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,1)
C
            ELSE
C
               CALL XPEN(XEEH,TE,GEVEN,NEV(IT))
               CALL XPEN(XOOH,TO,GODD,NOD(IT))
C
            END IF
C
         ELSE IF ( IREL.EQ.2 ) THEN
            IF ( IRUMP(IT).NE.1 .AND. IFM.NE.1 ) THEN
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELP(XEEH,TE,CLXE,NTOT(IT),GODD,PXO,NOD(IT),
     &                          GEVEN,PXE,NEV(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELP(XOOH,TO,CLXO,NTOT(IT),GEVEN,PXE,NEV(IT),
     &                          GODD,PXO,NOD(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
            ELSE IF ( IRUMP(IT).EQ.1 .AND. IFM.NE.1 ) THEN
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELP(XEEH,TE,CLXE,NTOT(IT),GODD,PXO,NOD(IT),
     &                          GEVEN,PXE,NEV(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELP(XOOH,TO,CLXO,NTOT(IT),GEVEN,PXE,NEV(IT),
     &                          GODD,PXO,NOD(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELPR(XEOH,TE,CLXEO,NTOT(IT),GEVOD,PXO,
     &                           NOD(IT),GODEV,PXE,NEV(IT),INAT,INATS,
     &                           NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELPR(XOEH,TO,CLXOE,NTOT(IT),GODEV,PXE,
     &                           NEV(IT),GEVOD,PXO,NOD(IT),INAT,INATS,
     &                           NAT)
                  END DO
               END DO
C
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,1)
               DO I = 1,MLQNAT
                  D1(I,I) = TE(I)
                  D2(I,I) = TO(I)
               END DO
               CALL XCOPR(HALM,D1,D2,D3,D4,NALM,XED,XOD,1)
C
            ELSE IF ( IRUMP(IT).NE.1 .AND. IFM.EQ.1 ) THEN
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
C
                     CALL XRELFE(XEEH,CLXE,NTOT(IT),GODD,PXO,NOD(IT),
     &                           GEVEN,PXE,NEV(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELFE(XOOH,CLXO,NTOT(IT),GEVEN,PXE,NEV(IT),
     &                           GODD,PXO,NOD(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO I = 1,MLQNAT
                  DO J = 1,MLQNAT
                     XEOH(I,J) = CZERO
                     XOEH(I,J) = CZERO
                  END DO
               END DO
C
               CALL XCOPR(XALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,1)
C
               CALL COPYTS(TSE,D1,XED,IT)
               CALL COPYTS(TSO,D2,XOD,IT)
               CALL COPYTSR(TSEO,D3,XED,XOD,IT)
               CALL COPYTSR(TSOE,D4,XOD,XED,IT)
C
               CALL XCOPR(HALM,D1,D2,D3,D4,NALM,XED,XOD,1)
               CALL MULT(HALM,XALM,VALM,NALM)
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
C
            ELSE IF ( IRUMP(IT).EQ.1 .AND. IFM.EQ.1 ) THEN
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELFE(XEEH,CLXE,NTOT(IT),GODD,PXO,NOD(IT),
     &                           GEVEN,PXE,NEV(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELFE(XOOH,CLXO,NTOT(IT),GEVEN,PXE,NEV(IT),
     &                           GODD,PXO,NOD(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELFER(XEOH,CLXEO,NTOT(IT),GEVOD,PXO,NOD(IT),
     &                            GODEV,PXE,NEV(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               DO INAT = 1,NAT
                  DO INATS = 1,NAT
                     CALL XRELFER(XOEH,CLXOE,NTOT(IT),GODEV,PXE,NEV(IT),
     &                            GEVOD,PXO,NOD(IT),INAT,INATS,NAT)
                  END DO
               END DO
C
               CALL XCOPR(XALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,1)
C
               CALL COPYTS(TSE,D1,XED,IT)
               CALL COPYTS(TSO,D2,XOD,IT)
               CALL COPYTSR(TSEO,D3,XED,XOD,IT)
               CALL COPYTSR(TSOE,D4,XOD,XED,IT)
C
               CALL XCOPR(HALM,D1,D2,D3,D4,NALM,XED,XOD,1)
               CALL MULT(HALM,XALM,VALM,NALM)
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
C
            END IF
         END IF
C
         IF ( IP.GT.4 ) THEN
            CALL CMATPR('xe      ',XEEH,XED,NOUT1)
            CALL CMATPR('xo      ',XOOH,XOD,NOUT1)
            CALL CMATPR('xeo     ',XEOH,XED,NOUT1)
            CALL CMATPR('xoe     ',XOEH,XOD,NOUT1)
         END IF
C
C         *** calculation of the inverse matrix of (1-x) ***
C
         IF ( IFM.EQ.1 .OR. IRUMP(IT).EQ.1 ) THEN
            IF ( REAL(EEG).LT.10.0D0 ) THEN
               CALL ONEMX(VALM,XALM,NALM)
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
            ELSE IF ( REAL(EEG).GT.10.0D0 .AND. REAL(EEG).LT.20.0D0 )
     &                THEN
               CALL ONEMX1(VALM,XALM,NALM)
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
            ELSE
               CALL ONEMX(VALM,XALM,NALM)
               CALL XCOPR(VALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
            END IF
C
            CALL TCOPY(XEEH,XFEE,XED,IT)
            CALL TCOPY(XOOH,XFOO,XOD,IT)
            CALL TCOPYR(XEOH,XFEO,XED,XOD,IT)
            CALL TCOPYR(XOEH,XFOE,XOD,XED,IT)
C
            IF ( IP.GT.4 ) THEN
               CALL CMATPR('i(1-ee))',XEEH,XED,NOUT1)
               CALL CMATPR('i(1-oo))',XOOH,XOD,NOUT1)
               CALL CMATPR('i(1-eo))',XEOH,XED,NOUT1)
               CALL CMATPR('i(1-oe))',XOEH,XOD,NOUT1)
            END IF
         ELSE IF ( IFM.NE.1 .AND. IRUMP(IT).NE.1 ) THEN
            IF ( REAL(EEG).LT.10.0D0 ) THEN
               CALL ONEMX(XEEH,XALM,XED)
               CALL TCOPY(XEEH,XFEE,XED,IT)
            ELSE IF ( REAL(EEG).GT.10.0D0 .AND. REAL(EEG).LT.20.0D0 )
     &                THEN
               CALL ONEMX1(XEEH,XALM,XED)
               CALL TCOPY(XEEH,XFEE,XED,IT)
            ELSE
               CALL ONEMX(XEEH,XALM,XED)
               CALL TCOPY(XEEH,XFEE,XED,IT)
            END IF
            IF ( REAL(EEG).LT.10.0D0 ) THEN
               CALL ONEMX(XOOH,XALM,XOD)
               CALL TCOPY(XOOH,XFOO,XOD,IT)
            ELSE IF ( REAL(EEG).GT.10.0D0 .AND. REAL(EEG).LT.20.0D0 )
     &                THEN
               CALL ONEMX1(XOOH,XALM,XOD)
               CALL TCOPY(XOOH,XFOO,XOD,IT)
            ELSE
               CALL ONEMX(XOOH,XALM,XOD)
               CALL TCOPY(XOOH,XFOO,XOD,IT)
            END IF
C
            IF ( IP.GT.4 ) THEN
               CALL CMATPR('i(1-ee))',XEEH,XED,NOUT1)
               CALL CMATPR('i(1-oo))',XOOH,XOD,NOUT1)
            END IF
         END IF
C
C         *** calculation of the inverse matrix of (1-x)*T ***
C
         IF ( IFM.EQ.1 .OR. IRUMP(IT).EQ.1 ) THEN
            CALL MULT(VALM,HALM,XALM,NALM)
            CALL XCOPR(XALM,XEEH,XOOH,XEOH,XOEH,NALM,XED,XOD,2)
         ELSE
            CALL MULTP(XEEH,TE,D1,XED,XED)
            CALL COPY(D1,XEEH,XED)
            CALL MULTP(XOOH,TO,D2,XOD,XOD)
            CALL COPY(D2,XOOH,XOD)
            IF ( IP.GT.4 ) THEN
               CALL CMATPR('i(1-ee)T',XEEH,XED,NOUT1)
               CALL CMATPR('i(1-oo)T',XOOH,XOD,NOUT1)
            END IF
         END IF
C
C         *** calculation of layer scattering matrices ***
C         *** calculation of mpm ***
C
         IF ( LHL.EQ.1 ) THEN
C
            CALL MMAT(C1,2,1,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSHM,EPOSHP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT1,QD,IT)
C
         ELSE IF ( LHL.EQ.2 ) THEN
C
            CALL MMAT(C1,2,1,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSLM,EPOSLP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT1,QD,IT)
C
         END IF
         CALL COPYQ(QT1,C1,QD,IT)
C
         IF ( IP.GT.3 ) WRITE (NOUT1,99007) IT
         IF ( IP.GT.3 ) CALL CMATPR(' mpm    ',C1,QD,NOUT1)
C
C         *** calculation of mpp ***
C
         IF ( LHL.EQ.1 ) THEN
C
            CALL MMAT(C2,2,2,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSHM,EPOSHP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C2,QT2,QD,IT)
C
         ELSE IF ( LHL.EQ.2 ) THEN
C
            CALL MMAT(C2,2,2,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSLM,EPOSLP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C2,QT2,QD,IT)
C
         END IF
         CALL COPYQ(QT2,C2,QD,IT)
C
         IF ( IP.GT.3 ) WRITE (NOUT1,99007) IT
         IF ( IP.GT.3 ) CALL CMATPR(' mpp    ',C2,QD,NOUT1)
C
C         *** calculation of mmm ***
C
         IF ( LHL.EQ.1 ) THEN
C
            CALL MMAT(C1,1,1,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSHM,EPOSHP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT3,QD,IT)
C
         ELSE IF ( LHL.EQ.2 ) THEN
C
            CALL MMAT(C1,1,1,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSLM,EPOSLP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT3,QD,IT)
C
         END IF
         CALL COPYQ(QT3,C1,QD,IT)
C
         IF ( IP.GT.3 ) WRITE (NOUT1,99007) IT
         IF ( IP.GT.3 ) CALL CMATPR(' mmm    ',C1,QD,NOUT1)
C
C         *** calculation of mmp ***
C
         IF ( LHL.EQ.1 ) THEN
C
            CALL MMAT(C1,1,2,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSHM,EPOSHP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT4,QD,IT)
C
         ELSE IF ( LHL.EQ.2 ) THEN
C
            CALL MMAT(C1,1,2,ME,MO,MSE,MSO,MCON,SUME,SUMO,PQ,XED,XOD,QD,
     &                SD,GANZ,EPOSLM,EPOSLP,IT,NATL,IRUMP(IT),SUMEO,
     &                SUMOE,IFM,XEEH,XEOH,XOEH,XOOH)
            CALL COPYBQ(C1,QT4,QD,IT)
C
         END IF
         CALL COPYQ(QT4,C1,QD,IT)
C
         IF ( IP.GT.3 ) WRITE (NOUT1,99007) IT
         IF ( IP.GT.3 ) CALL CMATPR(' mmp    ',C1,QD,NOUT1)
C
         CALL PROP(C1,2,KGFAK1,PQ,QD,SD,GANZ,LAYS,IT)
         CALL COPYBQ(C1,QH3,QD,IT)
         CALL PROP(C1,1,KGFAK1,PQ,QD,SD,GANZ,LAYS,IT)
         CALL COPYBQ(C1,QH4,QD,IT)
         CALL COPYQ(QH3,C1,QD,IT)
         IF ( IP.GT.3 ) WRITE (NOUT1,99008) IT
         IF ( IP.GT.3 ) CALL CMATPR(' bpp    ',C1,QD,NOUT1)
         CALL COPYQ(QH4,C1,QD,IT)
         IF ( IP.GT.3 ) WRITE (NOUT1,99008) IT
         IF ( IP.GT.3 ) CALL CMATPR(' bpm    ',C1,QD,NOUT1)
      END DO
C
      DO IXX = LAYB,LAYS
         JHELP(IXX-LAYB+1) = IXX
      END DO
C
      DO IL = LAYB,LAYS
         CALL COPYQ(QT1,QPM1,QD,JHELP(1))
         CALL COPYQ(QH3,C1,QD,JHELP(1))
         CALL COPYQ(QT2,C2,QD,JHELP(1))
         CALL MULT(C1,C2,QPP1,QD)
         CALL COPYQ(QT3,C1,QD,JHELP(1))
         CALL COPYQ(QH4,C2,QD,JHELP(1))
         CALL MULT(C1,C2,QMM1,QD)
         CALL COPYQ(QH3,C1,QD,JHELP(1))
         CALL COPYQ(QT4,C2,QD,JHELP(1))
         CALL MULT(C1,C2,C3,QD)
         CALL COPYQ(QH4,C1,QD,JHELP(1))
         CALL MULT(C3,C1,QMP1,QD)
         DO JL = 2,LAYS - LAYB + 1
            CALL COPYQ(QT1,QPM2,QD,JHELP(JL))
            CALL COPYQ(QH3,C1,QD,JHELP(JL))
            CALL COPYQ(QT2,C2,QD,JHELP(JL))
            CALL MULT(C1,C2,QPP2,QD)
            CALL COPYQ(QT3,C1,QD,JHELP(JL))
            CALL COPYQ(QH4,C2,QD,JHELP(JL))
            CALL MULT(C1,C2,QMM2,QD)
            CALL COPYQ(QH3,C1,QD,JHELP(JL))
            CALL COPYQ(QT4,C2,QD,JHELP(JL))
            CALL MULT(C1,C2,C3,QD)
            CALL COPYQ(QH4,C1,QD,JHELP(JL))
            CALL MULT(C3,C1,QMP2,QD)
            CALL ONE(C1ONE,QD)
            CALL SDOUBLE(QPP1,QPM1,QMP1,QMM1,QPP2,QPM2,QMP2,QMM2,C1ONE,
     &                   C1,C2,C3,C4,QD)
         END DO
C
         IF ( IBLOCH.NE.0 ) THEN
C             bulkrepeat-unit doubling
            IF ( IP.GT.3 ) CALL CMATPR(' qpm    ',QPM1,QD,NOUT1)
            IF ( IP.GT.3 ) CALL CMATPR(' qpp1   ',QPP1,QD,NOUT1)
            IF ( IP.GT.3 ) CALL CMATPR(' qmp1   ',QMP1,QD,NOUT1)
            IF ( IP.GT.3 ) CALL CMATPR(' qmm1   ',QMM1,QD,NOUT1)
            DO LAYER = 1,LANZ
               CALL ONE(C1ONE,QD)
               CALL COPY(QPP1,QPP2,QD)
               CALL COPY(QPM1,QPM2,QD)
               CALL COPY(QMP1,QMP2,QD)
               CALL COPY(QMM1,QMM2,QD)
               CALL SDOUBLE(QPP1,QPM1,QMP1,QMM1,QPP2,QPM2,QMP2,QMM2,
     &                      C1ONE,C1,C2,C3,C4,QD)
C
            END DO
            IF ( IL.EQ.LAYB ) THEN
               CALL COPY(QPP1,C5,QD)
               CALL COPY(QPM1,C6,QD)
               CALL COPY(QMP1,C7,QD)
               CALL COPY(QMM1,C8,QD)
            END IF
            IF ( IP.GT.3 ) CALL CMATPR(' btrans ',QPP1,QD,NOUT1)
            CALL COPYBQ(QPM1,BULKX,QD,IL)
            CALL COPYQ(BULKX,C1,QD,IL)
            IF ( IP.GT.3 ) WRITE (NOUT1,99009) JHELP(1)
            IF ( IP.GT.3 ) CALL CMATPR('rbulk1  ',C1,QD,NOUT1)
            KKK = JHELP(1)
            DO LXX = 2,LAYS - LAYB + 1
               JHELP(LXX-1) = JHELP(LXX)
            END DO
            JHELP(LAYS-LAYB+1) = KKK
         END IF
      END DO
C
      IF ( LAYB.GT.1 .AND. IBLOCH.NE.0 ) THEN
         DO IL = LAYB - 1,1, - 1
            CALL COPYQ(QT1,QPM1,QD,IL)
            CALL COPYQ(QH3,C1,QD,IL)
            CALL COPYQ(QT2,C2,QD,IL)
            CALL MULT(C1,C2,QPP1,QD)
            CALL COPYQ(QT3,C1,QD,IL)
            CALL COPYQ(QH4,C2,QD,IL)
            CALL MULT(C1,C2,QMM1,QD)
            CALL COPYQ(QH3,C1,QD,IL)
            CALL COPYQ(QT4,C2,QD,IL)
            CALL MULT(C1,C2,C3,QD)
            CALL COPYQ(QH4,C1,QD,IL)
            CALL MULT(C3,C1,QMP1,QD)
            CALL SDOUBLE(QPP1,QPM1,QMP1,QMM1,C5,C6,C7,C8,C1ONE,C1,C2,C3,
     &                   C4,QD)
            CALL COPYBQ(QPM1,BULKX,QD,IL)
            CALL COPYQ(BULKX,C1,QD,IL)
            IF ( IP.GT.3 ) WRITE (NOUT1,99009) IL
            IF ( IP.GT.3 ) CALL CMATPR('rbulk1  ',C1,QD,NOUT1)
            CALL COPY(QPP1,C5,QD)
            CALL COPY(QPM1,C6,QD)
            CALL COPY(QMP1,C7,QD)
            CALL COPY(QMM1,C8,QD)
         END DO
C
      END IF
C
      IF ( NSPLEE.EQ.0 ) THEN
         CALL ESENK(GANZ,QD,VPR,EEG,EEG,1,IREL,CLIGHT,KGZH,EZRYH,KGT,BK)
C
         CALL POTIN(VPR,VIX,1,ZPARU,ZPARD,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
C
         CALL RSPROP(PS,GANZ,QD,VPR,EEG,EEG,1,IREL,CLIGHT,BK,ADU,ADD,
     &               EPSX,TPM,TPP,TMP,TMM,IZBG,EZRY,ZAU,ZSU,ZBU,ZAD,ZSD,
     &               ZBD,IBLOCH)
C
         CALL BMAT(BPP,BPM,MPM,MPP,MMP,MMM,TPM,TPP,TMP,TMM,PS,QD,2,
     &             IBLOCH)
C
C        IF ( IP.GT.3 ) CALL CMATPR(' bppbar ',BPP,QD,NOUT1)
C        IF ( IP.GT.3 ) CALL CMATPR(' bpmbar ',BPM,QD,NOUT1)
C        IF ( IP.GT.3 ) CALL CMATPR(' mpmbar ',MPM,QD,NOUT1)
C        IF ( IP.GT.3 ) CALL CMATPR(' mppbar ',MPP,QD,NOUT1)
C        IF ( IP.GT.3 ) CALL CMATPR(' mmmbar ',MMM,QD,NOUT1)
C        IF ( IP.GT.3 ) CALL CMATPR(' mmpbar ',MMP,QD,NOUT1)
C
         CALL COPY(MPM,QPM1,QD)
         CALL MULT(BPP,MPP,QPP1,QD)
         CALL MULT(MMM,BPM,QMM1,QD)
         CALL MULT(BPP,MMP,C1,QD)
         CALL MULT(C1,BPM,QMP1,QD)
C
         IF ( IP.GT.2 ) CALL CMATPR(' qpmbar ',QPM1,QD,NOUT1)
         IF ( IP.GT.2 ) CALL CMATPR(' qppbar ',QPP1,QD,NOUT1)
         IF ( IP.GT.2 ) CALL CMATPR(' qmmbar ',QMM1,QD,NOUT1)
         IF ( IP.GT.2 ) CALL CMATPR(' qmpbar ',QMP1,QD,NOUT1)
C
         CALL ONE(C1ONE,QD)
         CALL SDOUBLE(QPP1,QPM1,QMP1,QMM1,C5,C6,C7,C8,C1ONE,C1,C2,C3,C4,
     &                QD)
C
C        calculate possible spleed beams and their directions
C
         EKIN = DBLE(EEG) - VPR
         CALL COPY(QPM1,C1,QD)
         IF ( IP.GT.2 ) CALL CMATPR('RBULKB  ',C1,QD,NOUT1)
         CALL REFLEX(EKIN,K0,RFXANZ,GRFX)
         CALL KGTSET(K0,GN,GANZ,KGS,EKIN,TAUW,THETAS,PHIS,ESCAPE)
         IF ( IP.GE.1 ) THEN
            DO G = 1,LG
               WRITE (NOUT1,99010) G,(GINT(G,I),I=1,2),(GN(G,I),I=1,2),
     &                             (KGS(G,1,I),I=1,3),THETAS(G)
     &                             *180.D0/PI,PHIS(G)*180.D0/PI
            END DO
         END IF
         CALL SGSET(SG,C1,NREL,PQ,GANZ,QD)
         CALL POLKMP(POLG,INTG,SG,INTFAK,POL0,POL0G,POL0L,GANZ,TYP,
     &               THETA,PHI,THETAS,PHIS,INTGPM,PK,PKGM,PN,PKGMN,PKN)
C
         CALL SPLOUT(INTG,INTGPM,POLG,TYP,GANZ,ISTR,GINT,N1,N2,EKIN,
     &               NREL,NTEST,THETA)
      END IF
C
      RETURN
C
99001 FORMAT (1x,'kx ky',4x,2E12.5)
99002 FORMAT (1x,'lattice vectors')
99003 FORMAT (20I3)
99004 FORMAT (10E12.5)
99005 FORMAT (1x,'layertyp ',i3,' atom in position ',i3)
99006 FORMAT (1x,2E12.5,4x,2E12.5)
99007 FORMAT (1x,'scattering matrix for layertyp',i3)
99008 FORMAT (1x,'bulk propagator for layertyp',i3)
99009 FORMAT (1x,'bulk reflection matrix of type',8I3)
99010 FORMAT (i3,2x,i3,1x,i3,2x,f5.2,1x,f5.2,2x,f5.2,1x,f5.2,1x,f5.2,2x,
     &        f7.2,1x,f7.2)
      END
C*==coeffh.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE COEFFH(XV,KGZBLK,QD,U0GP,CMG,WMHF,WMHG,WPHF,WPHG,WA,WB,
     &                  QE,IREL,IGP,IGVAL,IGPRO,BV,GANZ,TPM,TPP,TMP,TMM,
     &                  PS,U1GP,U0GM,CLIGHT,EMV1,K0,Q1,Q2,Q3,Q4,IPOL,
     &                  SPOL,BULKX,C1,C2,C3,C5,C6,IBLOCH)
C
      USE MOD_SPEC,ONLY:LAYSM,LH,MHDIM1,CZERO,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT
      COMPLEX*16 EMV1,Q1,Q2,Q3,Q4
      INTEGER GANZ,IBLOCH,IGP,IPOL,IREL,QD,QE,SPOL
      COMPLEX*16 BULKX(LH,LH,LAYSM),BV(QE),C1(QD,QD),C2(QD,QD),C3(QD,QD)
     &           ,C5(QD,QD),C6(QD,QD),CMG(LH),KGZBLK(LH),PS(LH,2),
     &           TMM(LH),TMP(LH),TPM(LH),TPP(LH),U0GM(LH),U0GP(LH),
     &           U1GP(LH),WA(QE,QE),WB(QE,QE),WMHF(MHDIM1,LH),
     &           WMHG(MHDIM1,LH),WPHF(MHDIM1,LH),WPHG(MHDIM1,LH),XV(QE)
      INTEGER IGPRO(LH),IGVAL(LH)
      REAL*8 K0(2)
C
C Local variables
C
      COMPLEX*16 BPM(:,:),BPP(:,:),KZ,MMM(:,:),MMP(:,:),MPM(:,:),
     &           MPP(:,:),U0
      REAL*8 CROOT2,KXY
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
      ALLOCATABLE BPM,BPP,MMM,MMP,MPM,MPP
      ALLOCATE (BPM(LH,LH),BPP(LH,LH),MMM(LH,LH),MMP(LH,LH))
      ALLOCATE (MPM(LH,LH),MPP(LH,LH))
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     # purpose:       calculates coefficients cg-, u0g+ which are     *
C                      used to couple whittaker functions to bulk      *
C                      solutions for energy e + omega                  *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       bmat      mult      onemx     add       inv       mvmult       *
C     /****************************************************************/
C     /* inputs */
C     /* outputs */
C     /* local */
C
C     /* working arrays */
C     /* common */
C
      CROOT2 = REAL(CDSQRT(DCMPLX(2.D0,0.D0)))
C
      KXY = K0(1)**2 + K0(2)**2
      IF ( IREL.EQ.1 ) THEN
         KZ = CDSQRT(2.D0*EMV1-KXY)
      ELSE
         KZ = CDSQRT(2.D0*EMV1+(EMV1/CLIGHT)**2-KXY)
      END IF
      U0 = -CIMAG/KZ
C
      DO I = 1,LH
         CMG(I) = CZERO
         U0GP(I) = CZERO
         U0GM(I) = CZERO
         U1GP(I) = CZERO
      END DO
C
      DO I = 1,QE
         DO J = 1,QE
            WA(I,J) = CZERO
            WB(I,J) = CZERO
         END DO
      END DO
C
      IF ( IREL.EQ.1 ) THEN
         DO I = 1,IGP
            WA(I,I) = -WMHF(MHDIM1,IGVAL(I))
            WA(I+IGP,I) = -WMHG(MHDIM1,IGVAL(I))
         END DO
      ELSE
         DO I = 1,IGP
            WA(I,I) = -WMHF(MHDIM1,IGVAL(I))
            WA(I+IGP,I+IGP) = -WMHF(MHDIM1,IGVAL(I))
            WA(I+2*IGP,I) = -WMHG(MHDIM1,IGVAL(I)+GANZ)
            WA(I+3*IGP,I+IGP) = -WMHG(MHDIM1,IGVAL(I)+GANZ)
         END DO
      END IF
C
C     calculating multiple scattering between surface-barrier and bulk
C
      CALL BMAT(BPP,BPM,MPM,MPP,MMP,MMM,TPM,TPP,TMP,TMM,PS,QD,2,IBLOCH)
C
      CALL COPYQ(BULKX,C6,QD,1)
C     CALL MULT(BPP,MPM,C1,QD)
      CALL MULT(BPP,MMP,C1,QD)
      CALL MULT(C1,BPM,C2,QD)
      CALL MULT(C2,C6,C1,QD)
      CALL ONEMX(C1,C2,QD)
      CALL MULT(BPP,MPP,C3,QD)
      CALL MULT(C1,C3,C2,QD)
      CALL COPY(C2,C5,QD)
      CALL MULT(C6,C2,C1,QD)
      CALL MULT(BPM,C1,C3,QD)
C
      IF ( IREL.EQ.1 ) THEN
         DO I = 1,IGP
            DO J = 1,IGP
               IF ( I.EQ.J ) THEN
                  WA(I,J+IGP) = 1.D0 + C3(IGVAL(I),IGVAL(J))
                  WA(I+IGP,J+IGP) = CIMAG*KGZBLK(IGVAL(I))
     &                              *(1.D0-C3(IGVAL(I),IGVAL(J)))
               ELSE
                  WA(I,J+IGP) = C3(IGVAL(I),IGVAL(J))
                  WA(I+IGP,J+IGP) = -CIMAG*KGZBLK(IGVAL(I))
     &                              *C3(IGVAL(I),IGVAL(J))
               END IF
            END DO
         END DO
      ELSE
         DO I = 1,IGP
            DO J = 1,IGP
               IF ( I.EQ.J ) THEN
                  WA(I,J+2*IGP) = 1.D0 + C3(IGVAL(I),IGVAL(J))
                  WA(I+IGP,J+3*IGP) = 1.D0 + C3(IGVAL(I)+GANZ,IGVAL(J)+
     &                                GANZ)
                  WA(I+2*IGP,J+2*IGP) = CIMAG*KGZBLK(IGVAL(I))
     &                                  *(1.D0-C3(IGVAL(I),IGVAL(J)))
                  WA(I+3*IGP,J+3*IGP) = CIMAG*KGZBLK(IGVAL(I)+GANZ)
     &                                  *(1.D0-
     &                                  C3(IGVAL(I)+GANZ,IGVAL(J)+GANZ))
               ELSE
                  WA(I,J+2*IGP) = C3(IGVAL(I),IGVAL(J))
                  WA(I+IGP,J+3*IGP) = C3(IGVAL(I)+GANZ,IGVAL(J)+GANZ)
                  WA(I,J+3*IGP) = C3(IGVAL(I),IGVAL(J)+GANZ)
                  WA(I+IGP,J+2*IGP) = C3(IGVAL(I)+GANZ,IGVAL(J))
                  WA(I+2*IGP,J+2*IGP) = -CIMAG*KGZBLK(IGVAL(I))
     &                                  *C3(IGVAL(I),IGVAL(J))
                  WA(I+3*IGP,J+3*IGP) = -CIMAG*KGZBLK(IGVAL(I)+GANZ)
     &                                  *C3(IGVAL(I)+GANZ,IGVAL(J)+GANZ)
                  WA(I+2*IGP,J+3*IGP) = -CIMAG*KGZBLK(IGVAL(I))
     &                                  *C3(IGVAL(I),IGVAL(J)+GANZ)
                  WA(I+3*IGP,J+2*IGP) = -CIMAG*KGZBLK(IGVAL(I)+GANZ)
     &                                  *C3(IGVAL(I)+GANZ,IGVAL(J))
               END IF
            END DO
         END DO
      END IF
C
      CALL INV(WA,WB,QE)
C
      DO I = 1,QE
         BV(I) = CZERO
      END DO
C
      IF ( IREL.EQ.1 ) THEN
         DO I = 1,1
            IF ( IGPRO(I).EQ.1 ) BV(I) = U0
            IF ( IGPRO(I).EQ.1 ) BV(I+IGP) = 1
         END DO
      ELSE
         DO I = 1,1
            IF ( SPOL.EQ.1 .AND. IPOL.EQ.1 ) THEN
               IF ( IGPRO(I).EQ.1 ) THEN
                  BV(I) = U0/CROOT2
                  BV(I+IGP) = U0/CROOT2
                  BV(I+2*IGP) = 1.D0/CROOT2
                  BV(I+3*IGP) = 1.D0/CROOT2
               END IF
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.1 ) THEN
               IF ( IGPRO(I).EQ.1 ) THEN
                  BV(I) = U0*Q1
                  BV(I+IGP) = U0*Q2
                  BV(I+2*IGP) = Q1
                  BV(I+3*IGP) = Q2
               END IF
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.2 ) THEN
               IF ( IGPRO(I).EQ.1 ) THEN
                  BV(I) = U0*Q3
                  BV(I+IGP) = U0*Q4
                  BV(I+2*IGP) = Q3
                  BV(I+3*IGP) = Q4
               END IF
            ELSE IF ( SPOL.EQ.4 ) THEN
               IF ( IGPRO(I).EQ.1 ) THEN
                  BV(I) = U0/CROOT2
                  BV(I+IGP) = U0/CROOT2
                  BV(I+2*IGP) = 1.D0/CROOT2
                  BV(I+3*IGP) = 1.D0/CROOT2
               END IF
            END IF
         END DO
      END IF
C
      CALL MVMULT(WB,BV,XV,QE)
C
      IF ( IREL.EQ.1 ) THEN
         DO I = 1,IGP
            CMG(IGVAL(I)) = XV(I)
            U0GP(IGVAL(I)) = XV(I+IGP)
         END DO
      ELSE
         DO I = 1,IGP
            IF ( SPOL.EQ.1 .AND. IPOL.EQ.1 ) THEN
               CMG(IGVAL(I)) = XV(I)
               CMG(IGVAL(I)+GANZ) = XV(I+IGP)
               U0GP(IGVAL(I)) = XV(I+2*IGP)
               U0GP(IGVAL(I)+GANZ) = XV(I+3*IGP)
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.1 ) THEN
               CMG(IGVAL(I)) = XV(I)
               CMG(IGVAL(I)+GANZ) = XV(I+IGP)
               U0GP(IGVAL(I)) = XV(I+2*IGP)
C              U0GP(IGVAL(I)+GANZ) = XV(I+3*IGP)
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.2 ) THEN
               CMG(IGVAL(I)) = XV(I)
               CMG(IGVAL(I)+GANZ) = XV(I+IGP)
C              U0GP(IGVAL(I)) = XV(I+2*IGP)
               U0GP(IGVAL(I)+GANZ) = XV(I+3*IGP)
            ELSE IF ( SPOL.EQ.4 ) THEN
               CMG(IGVAL(I)) = XV(I)
               CMG(IGVAL(I)+GANZ) = XV(I+IGP)
               U0GP(IGVAL(I)) = XV(I+2*IGP)
               U0GP(IGVAL(I)+GANZ) = XV(I+3*IGP)
            END IF
         END DO
      END IF
C
      DO I = 1,QD
         DO J = 1,QD
            U1GP(I) = U1GP(I) + C5(I,J)*U0GP(J)
         END DO
      END DO
C
      CALL MULT(BPM,C6,C1,QD)
      DO I = 1,QD
         DO J = 1,QD
            U0GM(I) = U0GM(I) + C1(I,J)*U1GP(J)
         END DO
      END DO
C
      IF ( IREL.EQ.1 ) THEN
         DO I = 1,IGP
            DO J = 1,MHDIM1
               WPHF(J,IGVAL(I)) = (U0/WPHF(MHDIM1,IGVAL(I)))
     &                            *WPHF(J,IGVAL(I))
               WPHG(J,IGVAL(I)) = (1.D0/WPHG(MHDIM1,IGVAL(I)))
     &                            *WPHG(J,IGVAL(I))
            END DO
         END DO
      ELSE
         DO I = 1,IGP
            IF ( SPOL.EQ.1 .AND. IPOL.EQ.1 ) THEN
               DO J = 1,MHDIM1
                  WPHF(J,IGVAL(I)) = (U0/WPHF(MHDIM1,IGVAL(I)))
     &                               *WPHF(J,IGVAL(I))/CROOT2
                  WPHG(J,IGVAL(I)) = (1.D0/WPHG(MHDIM1,IGVAL(I)))
     &                               *WPHG(J,IGVAL(I))/CROOT2
                  WPHF(J,IGVAL(I)+GANZ) = WPHF(J,IGVAL(I))
                  WPHG(J,IGVAL(I)+GANZ) = WPHG(J,IGVAL(I))
               END DO
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.1 ) THEN
               DO J = 1,MHDIM1
                  WPHF(J,IGVAL(I)) = Q1*U0*WPHF(J,IGVAL(I))
     &                               /WPHF(MHDIM1,IGVAL(I))
                  WPHG(J,IGVAL(I)) = Q1*WPHG(J,IGVAL(I))
     &                               /WPHG(MHDIM1,IGVAL(I))
                  WPHF(J,IGVAL(I)+GANZ) = Q2*U0*WPHF(J,IGVAL(I)+GANZ)
     &               /WPHF(MHDIM1,IGVAL(I)+GANZ)
                  WPHG(J,IGVAL(I)+GANZ) = Q2*WPHG(J,IGVAL(I)+GANZ)
     &               /WPHG(MHDIM1,IGVAL(I)+GANZ)
               END DO
            ELSE IF ( SPOL.EQ.2 .AND. IPOL.EQ.2 ) THEN
               DO J = 1,MHDIM1
                  WPHF(J,IGVAL(I)) = Q3*U0*WPHF(J,IGVAL(I))
     &                               /WPHF(MHDIM1,IGVAL(I))
                  WPHG(J,IGVAL(I)) = Q3*WPHG(J,IGVAL(I))
     &                               /WPHG(MHDIM1,IGVAL(I))
                  WPHF(J,IGVAL(I)+GANZ) = Q4*U0*WPHF(J,IGVAL(I)+GANZ)
     &               /WPHF(MHDIM1,IGVAL(I)+GANZ)
                  WPHG(J,IGVAL(I)+GANZ) = Q4*WPHG(J,IGVAL(I)+GANZ)
     &               /WPHG(MHDIM1,IGVAL(I)+GANZ)
               END DO
            ELSE IF ( SPOL.EQ.4 ) THEN
               DO J = 1,MHDIM1
                  WPHF(J,IGVAL(I)) = (U0/WPHF(MHDIM1,IGVAL(I)))
     &                               *WPHF(J,IGVAL(I))/CROOT2
                  WPHG(J,IGVAL(I)) = (1.D0/WPHG(MHDIM1,IGVAL(I)))
     &                               *WPHG(J,IGVAL(I))/CROOT2
                  WPHF(J,IGVAL(I)+GANZ) = WPHF(J,IGVAL(I))
                  WPHG(J,IGVAL(I)+GANZ) = WPHG(J,IGVAL(I))
               END DO
            END IF
         END DO
      END IF
C
      END
C*==ujgs.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE UJGS(QH3,QH4,QT2,QT4,UGSM,UGSP,LAYER1,QD,LAYS,BULKX,
     &                U1GP,LAYB,C1,C2,C3,C4,C5)
C
      USE MOD_SPEC,ONLY:LAYSM,LL,LH,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYB,LAYER1,LAYS,QD
      COMPLEX*16 BULKX(LH,LH,LAYSM),C1(QD,QD),C2(QD,QD),C3(QD,QD),
     &           C4(QD,QD),C5(QD,QD),QH3(LH,LH,LAYSM),QH4(LH,LH,LAYSM),
     &           QT2(LH,LH,LAYSM),QT4(LH,LH,LAYSM),U1GP(LH),UGSM(LH,LL),
     &           UGSP(LH,LL)
C
C Local variables
C
      INTEGER I,IBULK,ISEQ,J,L,LN,LSEQ
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     # purpose      :                                                 *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       copyq      mult      onemx                                     *
C     /****************************************************************/
C     /* input */
C     /* output */
C     /* local */
C     /* working arrays */
C
      DO I = 1,LH
         DO J = 1,LL
            UGSP(I,J) = CZERO
            UGSM(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,QD
         UGSP(I,1) = U1GP(I)
      END DO
C
      L = 1
      DO LN = 1,LAYER1
         IF ( LN.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LN.GT.1 ) THEN
            IBULK = LAYB
         END IF
         DO ISEQ = IBULK,LAYS
            L = L + 1
            LSEQ = ISEQ + 1
            IF ( ISEQ.EQ.LAYS ) LSEQ = LAYB
            CALL COPYQ(QH3,C1,QD,ISEQ)
            CALL COPYQ(QT4,C2,QD,ISEQ)
            CALL MULT(C1,C2,C3,QD)
            CALL COPYQ(QH4,C2,QD,ISEQ)
            CALL MULT(C3,C2,C1,QD)
            CALL COPYQ(BULKX,C5,QD,LSEQ)
            CALL MULT(C1,C5,C2,QD)
            CALL ONEMX(C2,C3,QD)
            CALL COPYQ(QH3,C1,QD,ISEQ)
            CALL COPYQ(QT2,C3,QD,ISEQ)
            CALL MULT(C1,C3,C4,QD)
            CALL MULT(C2,C4,C1,QD)
            DO I = 1,QD
               DO J = 1,QD
                  UGSP(I,L) = UGSP(I,L) + C1(I,J)*UGSP(J,L-1)
               END DO
            END DO
            CALL COPYQ(QH4,C2,QD,ISEQ)
            CALL MULT(C2,C5,C3,QD)
            DO I = 1,QD
               DO J = 1,QD
                  UGSM(I,L-1) = UGSM(I,L-1) + C3(I,J)*UGSP(J,L)
               END DO
            END DO
         END DO
      END DO
C
      END
C*==sumat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE SUMAT(MAXL,IREL,CLIGHT,EHIGH,GANZ,XE,XO,UGSP,UGSM,
     &                 LAYER1,AWR,IPOL,NATL,LAYS,EPOSHM,EPOSHP,XEO,XOE,
     &                 IRUMP,NOD,NEV,NTOT,IOD,IED,IOED,IEOD,KGT,LAYB,
     &                 IEDLAYER,IODLAYER,IEODLAYER,IOEDLAYER,BXY)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLS,LL,LG,MLQNAT,LH,LX,PI,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLY
      PARAMETER (MLY=2*ML+1)
C
C Dummy arguments
C
      INTEGER BXY,GANZ,IED,IEOD,IOD,IOED,IPOL,IREL,LAYB,LAYER1,LAYS,MAXL
      REAL*8 CLIGHT
      COMPLEX*16 EHIGH
      COMPLEX*16 AWR(MLS,LL,NATLM,2),EPOSHM(NATLM,LH,LAYSM),
     &           EPOSHP(NATLM,LH,LAYSM),UGSM(LH,LL),UGSP(LH,LL),
     &           XE(MLQNAT,MLQNAT,LAYSM),XEO(MLQNAT,MLQNAT,LAYSM),
     &           XO(MLQNAT,MLQNAT,LAYSM),XOE(MLQNAT,MLQNAT,LAYSM)
      INTEGER IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IODLAYER(0:LAYSM),
     &        IOEDLAYER(0:LAYSM),IRUMP(LAYSM),NATL(LAYSM),NEV(LAYSM),
     &        NOD(LAYSM),NTOT(LAYSM)
      REAL*8 KGT(2,LG)
C
C Local variables
C
      COMPLEX*16 AKM(:,:,:,:),XEEH(:,:),XEOH(:,:),XOEH(:,:),XOOH(:,:),
     &           YP1(:,:,:,:),YP2(:,:,:,:),YS1(:,:,:,:),YS2(:,:,:,:)
      INTEGER IAN,IBULK,IKMU,IL,ILAUF,IM,IRO,IRUM,ISI,IT,IWRITE,JLA,K,
     &        LK,LN,MAXLK,NAT,NSI,XED,XOD
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
      ALLOCATABLE YP1,YP2,YS1,YS2,AKM,XEEH,XEOH,XOEH,XOOH
      ALLOCATE (YP1(ML,MLY,LX,NATLM),YP2(ML,MLY,LX,NATLM))
      ALLOCATE (YS1(ML,MLY,LX,NATLM),YS2(ML,MLY,LX,NATLM))
      ALLOCATE (AKM(MLS,2,LL,NATLM),XEEH(MLQNAT,MLQNAT))
      ALLOCATE (XEOH(MLQNAT,MLQNAT),XOEH(MLQNAT,MLQNAT))
      ALLOCATE (XOOH(MLQNAT,MLQNAT))
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /****************************************************************/
C     # purpose:                                                       *
C                                                                      *
C     # subroutines and functions called from this routine:            *
C       akmj0     akmj      mumati    copyq    akmji    vmxv *
C     /****************************************************************/
C
C
C
C
C
C
C
      YP1(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YP2(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YS1(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YS2(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
C
      MAXLK = MAXL + 1
      IEDLAYER(0) = 0
      IODLAYER(0) = 0
      IEODLAYER(0) = 0
      IOEDLAYER(0) = 0
      IL = 1
      DO LN = 1,LAYER1
         IF ( LN.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LN.GT.1 ) THEN
            IBULK = LAYB
         END IF
         DO IT = IBULK,LAYS
            NAT = NATL(IT)
            DO IAN = 1,NAT
               CALL AKMJ0(MAXL,EHIGH,GANZ,CLIGHT,IREL,UGSP,UGSM,IPOL,1,
     &                    IAN,YP1,YP2,YS1,YS2,EPOSHM,EPOSHP,IT,IL,KGT,
     &                    LH,LL,LX)
            END DO
            IL = IL + 1
         END DO
      END DO
      IL = 1
      IED = 1
      IOD = 1
      IOED = 1
      IEOD = 1
      DO IT = 1,LAYS
         IF ( IREL.EQ.1 ) THEN
            XED = NEV(IT)
            XOD = NOD(IT)
         ELSE IF ( IREL.EQ.2 ) THEN
            XED = NTOT(IT)
            XOD = NTOT(IT)
         END IF
         NAT = NATL(IT)
         DO IAN = 1,NAT
            JLA = 0
            DO LK = 1,MAXLK
               DO NSI = 1,IREL
                  ISI = (-1)*(-1)**NSI
                  IF ( IREL.NE.2 ) THEN
                     K = -LK
                     ILAUF = 2*IDINT(DABS(DBLE(K))) - 1
                  ELSE
                     K = ISI*(LK-1) + (ISI-1)/2
                     ILAUF = 2*IDINT(DABS(DBLE(K)))
                  END IF
                  IF ( K.NE.0 ) THEN
                     DO IRO = 1,IREL
                        CALL AKMJI(K,MAXL,IRO,IREL,XED,XOD,IAN,NAT,IED,
     &                             IOD,IOED,IEOD,JLA)
                     END DO
                     JLA = JLA + ILAUF
                  END IF
               END DO
            END DO
         END DO
         IL = IL + 1
         IEDLAYER(IT) = IED - 1
         IODLAYER(IT) = IOD - 1
         IOEDLAYER(IT) = IOED - 1
         IEODLAYER(IT) = IEOD - 1
      END DO
      IED = IED - 1
      IOD = IOD - 1
      IOED = IOED - 1
      IEOD = IEOD - 1
      IL = 1
      DO LN = 1,LAYER1
         IF ( LN.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LN.GT.1 ) THEN
            IBULK = LAYB
         END IF
         DO IT = IBULK,LAYS
            IF ( IREL.EQ.1 ) THEN
               XED = NEV(IT)
               XOD = NOD(IT)
            ELSE IF ( IREL.EQ.2 ) THEN
               XED = NTOT(IT)
               XOD = NTOT(IT)
            END IF
            IRUM = IRUMP(IT)
            CALL COPYTS(XE,XEEH,XED,IT)
            CALL COPYTS(XO,XOOH,XOD,IT)
            CALL COPYTSR(XEO,XEOH,XED,XOD,IT)
            CALL COPYTSR(XOE,XOEH,XOD,XED,IT)
C
            CALL AKMJ(AKM,MAXL,XEEH,XOOH,IPOL,YP1,YP2,YS1,YS2,XEOH,XOEH,
     &                IRUM,IT,IL,XED,XOD,IEDLAYER,IODLAYER,IOEDLAYER,
     &                IEODLAYER,LX,BXY)
            IL = IL + 1
         END DO
      END DO
      IL = 1
      DO LN = 1,LAYER1
         IF ( LN.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LN.GT.1 ) THEN
            IBULK = LAYB
         END IF
         DO IT = IBULK,LAYS
            NAT = NATL(IT)
            DO IAN = 1,NAT
               IF ( IREL.EQ.2 ) THEN
                  DO IKMU = 1,2*MAXLK*MAXLK
                     AWR(IKMU,IL,IAN,IPOL)
     &                  = 4.0*PI*(AKM(IKMU,1,IL,IAN)+AKM(IKMU,2,IL,IAN))
                  END DO
               ELSE
                  DO IKMU = 1,MAXLK*MAXLK
                     AWR(IKMU,IL,IAN,IPOL) = 4.0*PI*AKM(IKMU,1,IL,IAN)
                  END DO
               END IF
C
               IF ( IP.GT.2 ) THEN
                  IF ( IREL.EQ.1 ) THEN
                     WRITE (NOUT1,99001) IAN,IL
                  ELSE
                     WRITE (NOUT1,99002) IAN,IL
                  END IF
                  IF ( IREL.EQ.1 ) THEN
                     IWRITE = MAXLK**2
                  ELSE
                     IWRITE = 2*MAXLK**2
                  END IF
                  WRITE (NOUT1,99003) (AWR(IM,IL,IAN,IPOL),IM=1,IWRITE)
               END IF
            END DO
            IL = IL + 1
         END DO
      END DO
C
      RETURN
C
99001 FORMAT (1x,'almj atom in position ',3x,i3,5x,'for layer',i3)
99002 FORMAT (1x,'akmj atom in position ',3x,i3,5x,'for layer',i3)
99003 FORMAT (5(e12.5,e12.5))
      END
C*==akmj.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     /****************************************************************/
C
C
C
      SUBROUTINE AKMJ(AKM,MAXL,XE,XO,IPOL,YP1,YP2,YS1,YS2,XEO,XOE,IRUM,
     &                IT,IL,XED,XOD,IEDLAYER,IODLAYER,IOEDLAYER,
     &                IEODLAYER,DIMA,BXY)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,LL,CZERO,CONE,CIMAG
      USE MOD_SPEC_COE1,ONLY:LSPOO,IIOO,POOO,POSOO,IGZOO,IROOO,IAOO
      USE MOD_SPEC_COE2,ONLY:IAOOS,LSPEE,IIEE,PEEE,PESEE,IGZEE,IROEE,
     &    IAEE
      USE MOD_SPEC_COE3,ONLY:IAEES,LSPOE,IIOE,POOE,PESOE,IGZOE,IROOE,
     &    IAOE
      USE MOD_SPEC_COE4,ONLY:IAOES,LSPEO,IIEO,PEEO,POSEO,IGZEO,IROEO,
     &    IAEO
      USE MOD_SPEC_COE5,ONLY:IAEOS,COO,CEE,COE,CEO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MLY
      PARAMETER (MLY=2*ML+1)
C
C Dummy arguments
C
      INTEGER BXY,DIMA,IL,IPOL,IRUM,IT,MAXL,XED,XOD
      COMPLEX*16 AKM(MQD,2,LL,NATLM),XE(XED,XED),XEO(XED,XOD),
     &           XO(XOD,XOD),XOE(XOD,XED),YP1(ML,MLY,DIMA,NATLM),
     &           YP2(ML,MLY,DIMA,NATLM),YS1(ML,MLY,DIMA,NATLM),
     &           YS2(ML,MLY,DIMA,NATLM)
      INTEGER IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IODLAYER(0:LAYSM),
     &        IOEDLAYER(0:LAYSM)
C
C Local variables
C
      INTEGER I1,I2,I3,IA,IEDA,IEDE,IEODA,IEODE,IODA,IODE,IOEDA,IOEDE,IX
      COMPLEX*16 IC(ML)
C
C*** End of declarations rewritten by SPAG
C
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /****************************************************************/
C     # purpose:                                                       *
C     /****************************************************************/
C
C
C
      IC(1) = CONE
      DO IA = 1,MAXL
         IC(IA+1) = IC(IA)*CIMAG
      END DO
      DO I1 = 1,MQD
         DO I2 = 1,2
            DO I3 = 1,NATLM
               AKM(I1,I2,IL,I3) = CZERO
            END DO
         END DO
      END DO
C      IEDA = 1 + (IT-1)*IED/LAYS
C      IEDE = IED/LAYS + (IT-1)*IED/LAYS
C      IODA = 1 + (IT-1)*IOD/LAYS
C      IODE = IOD/LAYS + (IT-1)*IOD/LAYS
C      IEODA = 1 + (IT-1)*IEOD/LAYS
C      IEODE = IEOD/LAYS + (IT-1)*IEOD/LAYS
C      IOEDA = 1 + (IT-1)*IOED/LAYS
C      IOEDE = IOED/LAYS + (IT-1)*IOED/LAYS
C
      IEDA = IEDLAYER(IT-1) + 1
      IEDE = IEDLAYER(IT)
      IODA = IODLAYER(IT-1) + 1
      IODE = IODLAYER(IT)
      IEODA = IEODLAYER(IT-1) + 1
      IEODE = IEODLAYER(IT)
      IOEDA = IOEDLAYER(IT-1) + 1
      IOEDE = IOEDLAYER(IT)
C
      IF ( IPOL.EQ.1 ) THEN
         DO IX = IEDA,IEDE
            AKM(IGZEE(IX),IROEE(IX),IL,IAEE(IX))
     &         = AKM(IGZEE(IX),IROEE(IX),IL,IAEE(IX)) + CEE(IX)
     &         *IC(LSPEE(IX))
     &         *(YP1(LSPEE(IX),IIEE(IX),IL+LL*(IROEE(IX)-1),IAEES(IX))
     &         +YS1(LSPEE(IX),IIEE(IX),IL+LL*(IROEE(IX)-1),IAEES(IX)))
     &         *XE(PESEE(IX),PEEE(IX))
         END DO
         DO IX = IODA,IODE
            AKM(IGZOO(IX),IROOO(IX),IL,IAOO(IX))
     &         = AKM(IGZOO(IX),IROOO(IX),IL,IAOO(IX)) + COO(IX)
     &         *IC(LSPOO(IX))
     &         *(YP1(LSPOO(IX),IIOO(IX),IL+LL*(IROOO(IX)-1),IAOOS(IX))
     &         +YS1(LSPOO(IX),IIOO(IX),IL+LL*(IROOO(IX)-1),IAOOS(IX)))
     &         *XO(POSOO(IX),POOO(IX))
         END DO
         IF ( IRUM.GT.0 .OR. BXY.EQ.1 ) THEN
            DO IX = IEODA,IEODE
               AKM(IGZEO(IX),IROEO(IX),IL,IAEO(IX))
     &            = AKM(IGZEO(IX),IROEO(IX),IL,IAEO(IX)) + CEO(IX)
     &            *IC(LSPEO(IX))
     &            *(YP1(LSPEO(IX),IIEO(IX),IL+LL*(IROEO(IX)-1),IAEOS(IX)
     &            )+YS1(LSPEO(IX),IIEO(IX),IL+LL*(IROEO(IX)-1),IAEOS(IX)
     &            ))*XOE(POSEO(IX),PEEO(IX))
            END DO
            DO IX = IOEDA,IOEDE
               AKM(IGZOE(IX),IROOE(IX),IL,IAOE(IX))
     &            = AKM(IGZOE(IX),IROOE(IX),IL,IAOE(IX)) + COE(IX)
     &            *IC(LSPOE(IX))
     &            *(YP1(LSPOE(IX),IIOE(IX),IL+LL*(IROOE(IX)-1),IAOES(IX)
     &            )+YS1(LSPOE(IX),IIOE(IX),IL+LL*(IROOE(IX)-1),IAOES(IX)
     &            ))*XEO(PESOE(IX),POOE(IX))
            END DO
         END IF
      ELSE IF ( IPOL.EQ.2 ) THEN
         DO IX = IEDA,IEDE
            AKM(IGZEE(IX),IROEE(IX),IL,IAEE(IX))
     &         = AKM(IGZEE(IX),IROEE(IX),IL,IAEE(IX)) + CEE(IX)
     &         *IC(LSPEE(IX))
     &         *(YP2(LSPEE(IX),IIEE(IX),IL+LL*(IROEE(IX)-1),IAEES(IX))
     &         +YS2(LSPEE(IX),IIEE(IX),IL+LL*(IROEE(IX)-1),IAEES(IX)))
     &         *XE(PESEE(IX),PEEE(IX))
         END DO
         DO IX = IODA,IODE
            AKM(IGZOO(IX),IROOO(IX),IL,IAOO(IX))
     &         = AKM(IGZOO(IX),IROOO(IX),IL,IAOO(IX)) + COO(IX)
     &         *IC(LSPOO(IX))
     &         *(YP2(LSPOO(IX),IIOO(IX),IL+LL*(IROOO(IX)-1),IAOOS(IX))
     &         +YS2(LSPOO(IX),IIOO(IX),IL+LL*(IROOO(IX)-1),IAOOS(IX)))
     &         *XO(POSOO(IX),POOO(IX))
         END DO
         IF ( IRUM.GT.0 .OR. BXY.EQ.1 ) THEN
            DO IX = IEODA,IEODE
               AKM(IGZEO(IX),IROEO(IX),IL,IAEO(IX))
     &            = AKM(IGZEO(IX),IROEO(IX),IL,IAEO(IX)) + CEO(IX)
     &            *IC(LSPEO(IX))
     &            *(YP2(LSPEO(IX),IIEO(IX),IL+LL*(IROEO(IX)-1),IAEOS(IX)
     &            )+YS2(LSPEO(IX),IIEO(IX),IL+LL*(IROEO(IX)-1),IAEOS(IX)
     &            ))*XOE(POSEO(IX),PEEO(IX))
            END DO
            DO IX = IOEDA,IOEDE
               AKM(IGZOE(IX),IROOE(IX),IL,IAOE(IX))
     &            = AKM(IGZOE(IX),IROOE(IX),IL,IAOE(IX)) + COE(IX)
     &            *IC(LSPOE(IX))
     &            *(YP2(LSPOE(IX),IIOE(IX),IL+LL*(IROOE(IX)-1),IAOES(IX)
     &            )+YS2(LSPOE(IX),IIOE(IX),IL+LL*(IROOE(IX)-1),IAOES(IX)
     &            ))*XEO(PESOE(IX),POOE(IX))
            END DO
         END IF
      END IF
C
      END
C*==copystreu.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     /****************************************************************/
C
      SUBROUTINE COPYSTREU(TE,TO,TSE,TSO,LAY)
C
      USE MOD_SPEC,ONLY:LAYSM,MLQNAT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAY
      COMPLEX*16 TE(MLQNAT),TO(MLQNAT),TSE(MLQNAT,MLQNAT,LAYSM),
     &           TSO(MLQNAT,MLQNAT,LAYSM)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /****************************************************************/
C     # purpose:       copy t-matrix from tse/tso to te/to             *
C                                                                      *
C     /****************************************************************/
C     /* input */
C     /* output */
C     /* local */
C
      DO I = 1,MLQNAT
         TE(I) = TSE(I,I,LAY)
         TO(I) = TSO(I,I,LAY)
      END DO
C
      END
C*==transpo.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C     /****************************************************************/
C
      SUBROUTINE TRANSPO(A,B,NALM)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NALM
      COMPLEX*16 A(NALM,NALM),B(NALM,NALM)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /* input */
C     /* output */
C     /* local */
C
      DO I = 1,NALM
         DO J = 1,NALM
            B(I,J) = A(J,I)
         END DO
      END DO
C
      END
