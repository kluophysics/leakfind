C*==initialstate.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     vbrun_subs.f
C
C     contains subroutines used only in vbrun:
C             initialstate
C             blmj0
C             vmxg
C             xalmm1
C             bkmu
C             intraf
C             ajgs
C             sphrm4
C             bjgs
C             cwl2
C             coeffl
C             cphi1g
C             djgs
C             gkmu
C             sortbkm
C
C     /****************************************************************/
C
C
      SUBROUTINE INITIALSTATE(GANZ,LANZ2,NSPIN,NREL,ELOW,PHI,THETA,N,
     &                        VPR,TEST,NSPLEE,IBLOCH,IREL,CLIGHT,CELM,
     &                        TYP,ISTR,POL0,POL0L,IPOL,NATL,POS,IRUMP,
     &                        LAYS,NOD,NEV,NTOT,BRUSEP,BULKX,NCLM,LAYB,
     &                        TSE,TSO,TSEO,TSOE,BXY,EHIGH,LAYER1,K1,AWR,
     &                        IAN,IT,MAXL,IRUM,IFM,RMATO,RMATS,LAYP,
     &                        OMHAR,ZMESH,AAZ,ADU,ADD,BK,EPSX,IGVAL,
     &                        IGPRO,PSI2G,IGP,VIL,ZPARU,ZPARD,IEDLAYER,
     &                        IODLAYER,IEODLAYER,IOEDLAYER,ROSUR,ROSULA,
     &                        ROTER,ROTELA,ROTRA,ROTRLA,DETR,DETI,NT,
     &                        SPOL,ROTRASD,ROTERSD,ROSURSD,ROTRLASD,
     &                        ROTELASD,ROSULASD)
C
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,NPM,MLS,LL,LG,MLQNAT,LH,MHDIM1,
     &    MHDIM2,EPS12,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MQDNAT
      PARAMETER (MQDNAT=2*MLQNAT)
C
C Dummy arguments
C
      REAL*8 AAZ,CLIGHT,DETI,DETR,EPSX,OMHAR,PHI,ROSUR,ROTER,ROTRA,TEST,
     &       THETA,VIL,VPR
      INTEGER BXY,GANZ,IAN,IBLOCH,IFM,IGP,IPOL,IREL,IRUM,IT,LANZ2,LAYB,
     &        LAYER1,LAYP,LAYS,MAXL,N,NCLM,NREL,NSPIN,NSPLEE,NT,SPOL,TYP
      COMPLEX*16 EHIGH,ELOW,K1
      REAL*8 ADD(3),ADU(3),BRUSEP(3),CELM(36124),POL0(3),POL0L(3),
     &       POS(3,NATLM,LAYSM),ROSULA(LL),ROTELA(LL),ROTRLA(LL),
     &       ZMESH(6),ZPARD(3),ZPARU(3)
      COMPLEX*16 AWR(MLS,LL,NATLM,2),BK(LH,2),BULKX(LH,LH,LAYSM),
     &           PSI2G(MHDIM2,LH,4),RMATO(MQD,MQD,LAYSM,NATLM),
     &           RMATS(MQD,MQD,LAYSM,NATLM),ROSULASD(LL),ROSURSD(4),
     &           ROTELASD(LL),ROTERSD(4),ROTRASD(4),ROTRLASD(LL),
     &           TSE(MLQNAT,MLQNAT,LAYSM),TSEO(MLQNAT,MLQNAT,LAYSM),
     &           TSO(MLQNAT,MLQNAT,LAYSM),TSOE(MLQNAT,MLQNAT,LAYSM)
      INTEGER IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IGPRO(LH),IGVAL(LH),
     &        IODLAYER(0:LAYSM),IOEDLAYER(0:LAYSM),IRUMP(LAYSM),ISTR(2),
     &        NATL(LAYSM),NEV(LAYSM),NOD(LAYSM),NTOT(LAYSM)
C
C Local variables
C
      COMPLEX*16 A1GP(:),A1GPS(:),A1GZ(:,:),A1GZS(:,:),AGSM(:,:),
     &           AGSP(:,:),BGSM(:,:),BKM0(:,:,:),BKM1(:,:,:),BKM2(:,:,:)
     &           ,C1(:,:),C2(:,:),C3(:,:),C4(:,:),C5(:,:),C6(:,:),
     &           C7(:,:),CMG(:),CMGS(:),D0GM(:),D0GP(:),D1GP(:),
     &           DGSM(:,:),DGSP(:,:),EPOSHM(:,:,:),EPOSHP(:,:,:),
     &           EPOSLM(:,:,:),EPOSLP(:,:,:),GKM2(:,:,:),H1(:,:),H2(:,:)
     &           ,KGZBLK(:),KGZL(:),PHI1G(:,:),PHI1GS(:,:),POTDER(:),
     &           PS(:,:),QH3(:,:,:),QH4(:,:,:),QMM(:,:),QMP(:,:),
     &           QPM(:,:),QPP(:,:),QT1(:,:,:),QT2(:,:,:),QT3(:,:,:),
     &           QT4(:,:,:),TMM(:),TMP(:),TPM(:),TPP(:),WMHF(:,:),
     &           WMHG(:,:),WMSTF(:,:),WMSTG(:,:),WPHF(:,:),WPHG(:,:),
     &           XE1H(:,:),XEEH(:,:),XEO1H(:,:),XEOH(:,:),XFEE(:,:,:),
     &           XFEO(:,:,:),XFOE(:,:,:),XFOO(:,:,:),XO1H(:,:),
     &           XOE1H(:,:),XOEH(:,:),XOOH(:,:)
      INTEGER DIMVAR,I,IBULK,IJ,IL,ILAYER,IZBG,J,K,LU,MAXG,NAT,NNSK,
     &        NNSR,NSK(:),NSR(:),QD,XED,XOD
      REAL*8 EZRY,EZRYL(:),GN(:,:),K0(2),KGT(:,:),RLVS(:),RLVX(:),
     &       RLVY(:),VS(:),VX(:),VY(:),ZAD,ZAU,ZBD,ZBU,ZSD,ZSU
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C1,C2,C3,C4,C5,C6,C7,H1,H2,GN
      ALLOCATABLE PS,VS,VX,VY,QH3,QH4,QT1,QT2,QT3,QT4,CMG
      ALLOCATABLE KGT,QMM,QMP,QPM,QPP,TMM,TMP,TPM,TPP
      ALLOCATABLE NSK,NSR,D0GM,A1GP,BKM0,BKM1,BKM2,D0GP,D1GP
      ALLOCATABLE GKM2,A1GZ,XE1H,XO1H,AGSM,AGSP,BGSM,CMGS,DGSM,DGSP
      ALLOCATABLE KGZL,WMHF,WMHG,WPHF,WPHG,XEEH,XEOH,XOEH,XOOH,XFEE
      ALLOCATABLE XFEO,XFOE,XFOO,RLVS,RLVX,RLVY,PHI1G,A1GPS,A1GZS
      ALLOCATABLE XEO1H,XOE1H,WMSTF,WMSTG,EZRYL,PHI1GS,KGZBLK,EPOSHM
      ALLOCATABLE POTDER,EPOSHP,EPOSLM,EPOSLP
C
      ALLOCATE (C1(LH,LH),C2(LH,LH))
      ALLOCATE (C3(LH,LH),C4(LH,LH),C5(LH,LH),C6(LH,LH),C7(LH,LH))
      ALLOCATE (H1(MQDNAT,MQDNAT),H2(MQDNAT,MQDNAT),GN(LG,2))
      ALLOCATE (PS(LH,2))
      ALLOCATE (VS(NPM),VX(NPM),VY(NPM))
      ALLOCATE (QH3(LH,LH,LAYSM),QH4(LH,LH,LAYSM),QT1(LH,LH,LAYSM))
      ALLOCATE (QT2(LH,LH,LAYSM),QT3(LH,LH,LAYSM),QT4(LH,LH,LAYSM))
      ALLOCATE (CMG(LH))
      ALLOCATE (KGT(2,LG))
      ALLOCATE (QMM(LH,LH),QMP(LH,LH),QPM(LH,LH),QPP(LH,LH),TMM(LH))
      ALLOCATE (TMP(LH),TPM(LH),TPP(LH))
      ALLOCATE (NSK(NPM),NSR(NPM))
      ALLOCATE (D0GM(LH),A1GP(LH),BKM0(MLS,LL,NATLM))
      ALLOCATE (BKM1(MLS,LL,NATLM),BKM2(MLS,LL,NATLM),D0GP(LH))
      ALLOCATE (D1GP(LH),GKM2(MLS,LL,NATLM),A1GZ(MHDIM1,LH))
      ALLOCATE (XE1H(MLQNAT,MLQNAT),XO1H(MLQNAT,MLQNAT),AGSM(LH,LL))
      ALLOCATE (AGSP(LH,LL),BGSM(LH,LL),CMGS(LH),DGSM(LH,LL))
      ALLOCATE (DGSP(LH,LL),KGZL(LH),WMHF(MHDIM1,LH))
      ALLOCATE (WMHG(MHDIM1,LH),WPHF(MHDIM1,LH),WPHG(MHDIM1,LH))
      ALLOCATE (XEEH(MLQNAT,MLQNAT),XEOH(MLQNAT,MLQNAT))
      ALLOCATE (XOEH(MLQNAT,MLQNAT),XOOH(MLQNAT,MLQNAT))
      ALLOCATE (XFEE(MLQNAT,MLQNAT,LAYSM),XFEO(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (XFOE(MLQNAT,MLQNAT,LAYSM),XFOO(MLQNAT,MLQNAT,LAYSM))
      ALLOCATE (RLVS(NPM),RLVX(NPM),RLVY(NPM),PHI1G(MHDIM2,LH))
      ALLOCATE (A1GPS(LH),A1GZS(MHDIM1,LH),XEO1H(MLQNAT,MLQNAT))
      ALLOCATE (XOE1H(MLQNAT,MLQNAT),WMSTF(6,LH),WMSTG(6,LH))
      ALLOCATE (EZRYL(LH),PHI1GS(MHDIM2,LH),KGZBLK(LH))
      ALLOCATE (EPOSHM(NATLM,LH,LAYSM),POTDER(MHDIM2))
      ALLOCATE (EPOSHP(NATLM,LH,LAYSM),EPOSLM(NATLM,LH,LAYSM))
      ALLOCATE (EPOSLP(NATLM,LH,LAYSM))
C
C*** End of declarations rewritten by SPAG
C
C     /****************************************************************/
C     subroutines and functions called from this routine               *
C     splee       esenk       blmj0       copyq                 *
C     vmxg        bkmu        intraf      ajgs        potin            *
C     rsprop      ca1gp       bjgs        cwl2        coeffl           *
C     cphi1g      djgs        gkmu        rhosva                       *
C     /****************************************************************/
C
      MAXG = LG
      DO I = 1,MLQNAT
         DO J = 1,MLQNAT
            XE1H(I,J) = CZERO
            XEO1H(I,J) = CZERO
            XOE1H(I,J) = CZERO
            XO1H(I,J) = CZERO
            XEEH(I,J) = CZERO
            XEOH(I,J) = CZERO
            XOEH(I,J) = CZERO
            XOOH(I,J) = CZERO
            DO K = 1,LAYSM
               XFEE(I,J,K) = CZERO
               XFEO(I,J,K) = CZERO
               XFOE(I,J,K) = CZERO
               XFOO(I,J,K) = CZERO
            END DO
         END DO
      END DO
C
      CALL SPLEE(MAXG,GANZ,LANZ2,NSPIN,NREL,ELOW,GN,PHI,THETA,K0,QH3,
     &           QH4,QD,2,N,VPR,TEST,NSPLEE,IBLOCH,IREL,CLIGHT,NT,CELM,
     &           TYP,ISTR,POL0,POL0L,IPOL,NATL,POS,NNSK,NNSR,NSK,NSR,
     &           RLVS,RLVX,RLVY,VS,VX,VY,IRUMP,LAYS,EPOSLM,EPOSLP,
     &           EPOSHM,EPOSHP,NOD,NEV,NTOT,BRUSEP,QPM,QMP,QPP,QMM,QT1,
     &           QT2,QT3,QT4,BULKX,NCLM,KGZBLK,KGT,LAYB,TSE,TSO,TSEO,
     &           TSOE,XFEE,XFOO,XFEO,XFOE,IFM,VIL,EPSX,ZPARU,ZPARD)
C
      IF ( QD.GT.LH ) THEN
         WRITE (*,*) 'stopped: qd=',QD,' > lh=',LH
         STOP
      END IF
C
      CALL ESENK(GANZ,QD,VPR,ELOW,EHIGH,2,IREL,CLIGHT,KGZL,EZRYL,KGT,BK)
C
C     *** calculating intra-layer contribution ***
C
      IL = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
            IBULK = LAYB
         END IF
C
         DO IT = IBULK,LAYS
            NAT = NATL(IT)
            DO IAN = 1,NAT
               CALL BLMJ0(RMATO,K1,AWR,BKM0,IT,IAN,IL,IPOL)
            END DO
            IL = IL + 1
         END DO
      END DO
C
      DO I = 1,MLS
         DO J = 1,LL
            DO K = 1,NATLM
               BKM1(I,J,K) = CZERO
               BKM2(I,J,K) = CZERO
            END DO
         END DO
      END DO
C
      IL = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
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
            NAT = NATL(IT)
C
            CALL COPYTS(XFEE,XEEH,XED,IT)
            CALL COPYTS(XFOO,XOOH,XOD,IT)
            CALL COPYTSR(XFEO,XEOH,XED,XOD,IT)
            CALL COPYTSR(XFOE,XOEH,XOD,XED,IT)
C
C       *** overwrite bkm0 with g^(-1)*((1-X)^(-1))-1) to bkm2 ***
C       *** overwrite bkm0 with g^(-1)*(1-X)^(-1)*g    to bkm1 ***
C
            DIMVAR = 2*XED
            CALL VMXG(IT,XEEH,XOOH,XEOH,XOEH,XE1H,XO1H,XEO1H,XOE1H,XED,
     &                XOD,DIMVAR,H1,H2,TSE,TSEO,TSOE,TSO)
C
            DO IAN = 1,NAT
               CALL BKMU(BKM0,IREL,MAXL,XEEH,XOOH,BKM1,BKM2,XED,XOD,IAN,
     &                   NAT,XEOH,XOEH,IRUM,IL,XE1H,XO1H,XEO1H,XOE1H,
     &                   BXY)
            END DO
C
            IL = IL + 1
         END DO
      END DO
C
      CALL INTRAF(BKM2,RMATS,AWR,NATL,LAYS,LAYP,ROTRA,ROTRLA,LAYB,
     &            LAYER1,IPOL,SPOL,ROTRASD,ROTRLASD)
C
C     ** calculate surface and inter-layer contributions
C
      CALL AJGS(IREL,BKM1,GANZ,AGSP,AGSM,LAYER1,K1,NATL,LAYS,EPOSLM,
     &          EPOSLP,KGT,LAYB)
C
      IF ( IP.GE.3 ) THEN
         DO ILAYER = 1,LAYP
            IF ( IREL.EQ.1 ) THEN
               WRITE (NOUT1,99001) ILAYER
            ELSE
               WRITE (NOUT1,99002) ILAYER
            END IF
            WRITE (NOUT1,99007) (AGSP(IJ,ILAYER),IJ=1,QD)
            WRITE (NOUT1,99007) (AGSM(IJ,ILAYER),IJ=1,QD)
         END DO
      END IF
C
      CALL POTIN(VPR,VIL,2,ZPARU,ZPARD,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
C
      CALL RSPROP(PS,GANZ,QD,VPR,ELOW,EHIGH,2,IREL,CLIGHT,BK,ADU,ADD,
     &            EPSX,TPM,TPP,TMP,TMM,IZBG,EZRY,ZAU,ZSU,ZBU,ZAD,ZSD,
     &            ZBD,IBLOCH)
C
      CALL CA1GP(KGZBLK,GANZ,PSI2G,A1GP,QD,IREL,AAZ,EZRY,A1GZ,OMHAR,TPM,
     &           A1GPS,POTDER,A1GZS,VPR,IGP,IGVAL,IPOL)
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
      CALL BJGS(QD,QH3,QH4,QT3,QT4,AGSP,AGSM,BGSM,BULKX,LAYS,LAYP,
     &          LAYER1,LAYB,C1,C2,C3,C4,C5)
C
      IF ( IP.GE.3 ) THEN
         DO ILAYER = 1,LAYP
            IF ( IREL.EQ.1 ) THEN
               WRITE (NOUT1,99003) ILAYER
            ELSE
               WRITE (NOUT1,99004) ILAYER
            END IF
            WRITE (NOUT1,99007) (BGSM(IJ,ILAYER),IJ=1,QD)
         END DO
      END IF
C
      CALL CWL2(KGZL,GANZ,EZRYL,WMHF,WMHG,WPHF,WPHG,IGP,IGVAL,IGPRO,
     &          WMSTF,WMSTG,ZMESH,IZBG)
C
      DO I = 1,LH
         DO J = 1,LH
            C1(I,J) = CZERO
            C2(I,J) = CZERO
            C3(I,J) = CZERO
            C4(I,J) = CZERO
            C5(I,J) = CZERO
            C6(I,J) = CZERO
            C7(I,J) = CZERO
         END DO
      END DO
      CALL COEFFL(KGZBLK,QD,WMHF,WMHG,IGVAL,IREL,IGP,GANZ,A1GP,A1GPS,
     &            BGSM,TPM,TPP,TMP,TMM,PS,D0GM,D1GP,BULKX,CMG,CMGS,D0GP,
     &            DETR,DETI,C1,C2,C3,C4,C5,C6,C7,IBLOCH)
C
      IF ( IP.GE.3 ) THEN
         IF ( IREL.EQ.1 ) THEN
            WRITE (NOUT1,99005) 0
         ELSE
            WRITE (NOUT1,99006) 0
         END IF
         WRITE (NOUT1,99007) (D0GP(IJ),IJ=1,QD)
         WRITE (NOUT1,99007) (D0GM(IJ),IJ=1,QD)
      END IF
C
      CALL CPHI1G(IREL,GANZ,IGP,QD,IGVAL,ZMESH,EZRYL,A1GZ,A1GZS,CMG,
     &            CMGS,WMHF,WMHG,D0GM,D0GP,KGZBLK,PHI1G,PHI1GS)
C
C     *** calculate surface contribution ***
C
C     rosur = 0.d0
      DO I = 1,LL
         ROSULA(I) = 0.D0
      END DO
C
      IF ( ABS(AAZ).GT.EPS12 ) CALL RHOSVA(OMHAR,AAZ,EZRY,PSI2G,PHI1G,
     &     POTDER,ROSUR,ROSULA,GANZ,IPOL,SPOL,ROSURSD,ROSULASD)
C
C     *** calculate inter-layer contribution ***
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
      CALL DJGS(QT2,QT4,BGSM,AGSP,QD,LAYER1,QH3,QH4,LAYS,DGSM,DGSP,
     &          BULKX,D1GP,LAYB,C1,C2,C3,C4,C5,C6)
C
      IF ( IP.GE.3 ) THEN
         DO ILAYER = 1,LAYP
            IF ( IREL.EQ.1 ) THEN
               WRITE (NOUT1,99005) ILAYER
            ELSE
               WRITE (NOUT1,99006) ILAYER
            END IF
            WRITE (NOUT1,99007) (DGSP(IJ,ILAYER),IJ=1,QD)
            WRITE (NOUT1,99007) (DGSM(IJ,ILAYER),IJ=1,QD)
         END DO
      END IF
C
      IF ( IPOL.EQ.1 ) THEN
         CALL GKMU(IREL,MAXL,ELOW,GANZ,CLIGHT,XFEE,XFOO,GKM2,DGSM,DGSP,
     &             LAYER1,NEV,NOD,1,NATL,LAYS,EPOSLM,EPOSLP,XFEO,XFOE,
     &             IRUMP,NTOT,KGT,LAYB,IEDLAYER,IODLAYER,IEODLAYER,
     &             IOEDLAYER,BXY)
      ELSE IF ( IPOL.EQ.2 ) THEN
         CALL GKMU(IREL,MAXL,ELOW,GANZ,CLIGHT,XFEE,XFOO,GKM2,DGSM,DGSP,
     &             LAYER1,NEV,NOD,2,NATL,LAYS,EPOSLM,EPOSLP,XFEO,XFOE,
     &             IRUMP,NTOT,KGT,LAYB,IEDLAYER,IODLAYER,IEODLAYER,
     &             IOEDLAYER,BXY)
      ELSE IF ( IPOL.EQ.3 ) THEN
         CALL GKMU(IREL,MAXL,ELOW,GANZ,CLIGHT,XFEE,XFOO,GKM2,DGSM,DGSP,
     &             LAYER1,NEV,NOD,2,NATL,LAYS,EPOSLM,EPOSLP,XFEO,XFOE,
     &             IRUMP,NTOT,KGT,LAYB,IEDLAYER,IODLAYER,IEODLAYER,
     &             IOEDLAYER,BXY)
      ELSE IF ( IPOL.EQ.4 ) THEN
         CALL GKMU(IREL,MAXL,ELOW,GANZ,CLIGHT,XFEE,XFOO,GKM2,DGSM,DGSP,
     &             LAYER1,NEV,NOD,1,NATL,LAYS,EPOSLM,EPOSLP,XFEO,XFOE,
     &             IRUMP,NTOT,KGT,LAYB,IEDLAYER,IODLAYER,IEODLAYER,
     &             IOEDLAYER,BXY)
      END IF
C
      CALL INTRAF(GKM2,RMATS,AWR,NATL,LAYS,LAYP,ROTER,ROTELA,LAYB,
     &            LAYER1,IPOL,SPOL,ROTERSD,ROTELASD)
C
      RETURN
C
99001 FORMAT (1x,'ajg+ and ajg-        for layer',4x,i3)
99002 FORMAT (1x,'ajgs+(+,-) and ajgs-(+,-) for layer',4x,i3)
99003 FORMAT (1x,'bjg- for        layer',4x,i3)
99004 FORMAT (1x,'bjgs-(+,-) for layer',4x,i3)
99005 FORMAT (1x,'djg+ and djg-        for layer',4x,i3)
99006 FORMAT (1x,'djgs+(+,-) and djgs-(+,-) for layer',4x,i3)
99007 FORMAT (5(e12.5,e12.5))
      END
C*==blmj0.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
      SUBROUTINE BLMJ0(RMATO,K1,AWR,BKM0,IT,IAN,J,IPOL)
C     /****************************************************************/
C     calls:  rind        sortbkm
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,LL,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IAN,IPOL,IT,J
      COMPLEX*16 K1
      COMPLEX*16 AWR(MQD,LL,NATLM,2),BKM0(MQD,LL,NATLM),
     &           RMATO(MQD,MQD,LAYSM,NATLM)
C
C Local variables
C
      INTEGER DREL(:),I,I1,I11,I2,I21,NB(:),SB(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NB,SB,DREL
      ALLOCATE (NB(MQD),SB(MQD),DREL(MQD))
C
C
      CALL SORTBKM(NB,SB,DREL)
C
      DO I = 1,MQD
         BKM0(I,J,IAN) = CZERO
      END DO
C
      IF ( IPOL.EQ.1 ) THEN
         DO I1 = 1,MQD
            I11 = REL(NB(I1),SB(I1))
            DO I2 = 1,MQD
               I21 = REL(NB(I2),SB(I2))
               BKM0(I1,J,IAN) = BKM0(I1,J,IAN) + RMATO(I11,I21,IT,IAN)
     &                          *DCONJG(AWR(I2,J,IAN,1))
            END DO
         END DO
      ELSE IF ( IPOL.EQ.2 ) THEN
         DO I1 = 1,MQD
            I11 = REL(NB(I1),SB(I1))
            DO I2 = 1,MQD
               I21 = REL(NB(I2),SB(I2))
               BKM0(I1,J,IAN) = BKM0(I1,J,IAN) + RMATO(I11,I21,IT,IAN)
     &                          *DCONJG(AWR(I2,J,IAN,2))
            END DO
         END DO
      ELSE IF ( IPOL.EQ.3 ) THEN
         DO I1 = 1,MQD
            I11 = REL(NB(I1),SB(I1))
            DO I2 = 1,MQD
               I21 = REL(NB(I2),SB(I2))
               BKM0(I1,J,IAN) = BKM0(I1,J,IAN) + RMATO(I11,I21,IT,IAN)
     &                          *DCONJG(AWR(I2,J,IAN,2))
            END DO
         END DO
      ELSE IF ( IPOL.EQ.4 ) THEN
         DO I1 = 1,MQD
            I11 = REL(NB(I1),SB(I1))
            DO I2 = 1,MQD
               I21 = REL(NB(I2),SB(I2))
               BKM0(I1,J,IAN) = BKM0(I1,J,IAN) + RMATO(I11,I21,IT,IAN)
     &                          *DCONJG(AWR(I2,J,IAN,1))
            END DO
         END DO
      END IF
C
      DO I1 = 1,MQD
         BKM0(I1,J,IAN) = K1*BKM0(I1,J,IAN)
      END DO
C
      IF ( IP.GE.3 ) THEN
         WRITE (NOUT1,99001) IAN
         WRITE (NOUT1,99002) J
         WRITE (NOUT1,99003) (BKM0(I,J,IAN),I=1,MQD)
      END IF
      RETURN
C
99001 FORMAT (1x,'atom in position ',i3)
99002 FORMAT (1x,'blmj0 for              layer',4x,i3)
99003 FORMAT (5(1x,2E13.5))
      END
C*==vmxg.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
      SUBROUTINE VMXG(LAY,XEEH,XOOH,XEOH,XOEH,XE1H,XO1H,XEO1H,XOE1H,XED,
     &                XOD,DIMVAR,H1,H2,TSE,TSEO,TSOE,TSO)
C     /****************************************************************/
C     calls:
C         rind        sortbkm     xalmm1
C         xcopr       mult        inv
C     /****************************************************************/
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,MLQNAT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MQDNAT
      PARAMETER (MQDNAT=MQD*NATLM)
C
C Dummy arguments
C
      INTEGER DIMVAR,LAY,XED,XOD
      COMPLEX*16 H1(MQDNAT,MQDNAT),H2(MQDNAT,MQDNAT),
     &           TSE(MLQNAT,MLQNAT,LAYSM),TSEO(MLQNAT,MLQNAT,LAYSM),
     &           TSO(MLQNAT,MLQNAT,LAYSM),TSOE(MLQNAT,MLQNAT,LAYSM),
     &           XE1H(XED,XED),XEEH(XED,XED),XEO1H(XED,XOD),
     &           XEOH(XED,XOD),XO1H(XOD,XOD),XOE1H(XOD,XED),
     &           XOEH(XOD,XED),XOOH(XOD,XOD)
C
C Local variables
C
      COMPLEX*16 D1(:,:),D2(:,:),D3(:,:),D4(:,:),GALM(:,:),HALM(:,:),
     &           XALM(:,:)
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE D1,D2,D3,D4,GALM,HALM,XALM
      ALLOCATE (D1(MLQNAT,MLQNAT),D2(MLQNAT,MLQNAT))
      ALLOCATE (D3(MLQNAT,MLQNAT),D4(MLQNAT,MLQNAT))
      ALLOCATE (GALM(MQDNAT,MQDNAT),HALM(MQDNAT,MQDNAT))
      ALLOCATE (XALM(MQDNAT,MQDNAT))
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MLQNAT
         DO J = 1,MLQNAT
            D1(I,J) = CZERO
            D2(I,J) = CZERO
            D3(I,J) = CZERO
            D4(I,J) = CZERO
         END DO
      END DO
C
      CALL COPYTS(TSE,D1,XED,LAY)
      CALL COPYTS(TSO,D2,XOD,LAY)
      CALL COPYTS(TSEO,D3,XED,LAY)
      CALL COPYTS(TSOE,D4,XOD,LAY)
C
      DO I = 1,MQDNAT
         DO J = 1,MQDNAT
            XALM(I,J) = CZERO
            GALM(I,J) = CZERO
            HALM(I,J) = CZERO
            H1(I,J) = CZERO
            H2(I,J) = CZERO
         END DO
      END DO
C
      CALL XCOPR(XALM,XEEH,XOOH,XEOH,XOEH,DIMVAR,XED,XOD,1)
      CALL XCOPR(GALM,D1,D2,D3,D4,DIMVAR,XED,XOD,1)
C
      CALL INV(GALM,HALM,DIMVAR)
C
C
C     ---------bkmu1--for--ajgs-inter-layer-contribution-------
C
C
      CALL MULT(HALM,XALM,H1,DIMVAR)
      CALL MULT(H1,GALM,H2,DIMVAR)
      CALL XCOPR(H2,XEEH,XOOH,XEOH,XOEH,DIMVAR,XED,XOD,2)
C
C     ---------bkmu2--for--intraf-intra-layer-contribution-------
C
C     IF ( NCPA.EQ.0 ) THEN
      CALL XALMM1(XALM,H2,DIMVAR)
      CALL MULT(HALM,H2,H1,DIMVAR)
      CALL XCOPR(H1,XE1H,XO1H,XEO1H,XOE1H,DIMVAR,XED,XOD,2)
C     ELSE IF ( ATA.EQ.1 ) THEN
C        CALL XALMM1(XALM,H2,DIM)
C        CALL MULT(HALM,H2,H1,DIM)
C        CALL XCOPR(H1,XE1H,XO1H,XEO1H,XOE1H,DIM,XED,XOD,2)
C     ELSE IF ( ATA.EQ.0 ) THEN
C        CALL MULT(HALM,XALM,H1,DIM)
C        CALL XCOPR(H1,XE1H,XO1H,XEO1H,XOE1H,DIM,XED,XOD,2)
C     END IF
C
      END
C*==xalmm1.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C
C
      SUBROUTINE XALMM1(XALM,H2,DIMVAR)
C
      USE MOD_SPEC,ONLY:CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 H2(DIMVAR,DIMVAR),XALM(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            H2(I,J) = XALM(I,J)
            IF ( I.EQ.J ) H2(I,J) = XALM(I,J) - CONE
         END DO
      END DO
C
      END
C*==bkmu.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
      SUBROUTINE BKMU(BKM0,IREL,MAXL,XE,XO,BKM1,BKM2,XED,XOD,IAN,NAT,
     &                XEO,XOE,IRUM,IL,XE1,XO1,XEO1,XOE1,BXY)
C
      USE MOD_SPEC,ONLY:NATLM,ML,MLS,LL,CZERO
      USE MOD_SPEC_LMKMS,ONLY:LME,LMO,KME,KMO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BXY,IAN,IL,IREL,IRUM,MAXL,NAT,XED,XOD
      COMPLEX*16 BKM0(MLS,LL,NATLM),BKM1(MLS,LL,NATLM),
     &           BKM2(MLS,LL,NATLM),XE(XED,XED),XE1(XED,XED),
     &           XEO(XED,XOD),XEO1(XED,XOD),XO(XOD,XOD),XO1(XOD,XOD),
     &           XOE(XOD,XED),XOE1(XOD,XED)
C
C Local variables
C
      INTEGER IANS,IFE1,IFE2,II1,II2,INDE1,INDE2,INDES1,INDES2,INDO1,
     &        INDO2,INDOS1,INDOS2,IO,ISI,ISIS,IU,IWRITE,JU,K,KAP1,KAP2,
     &        KS,LK,LKS,MAXLK,MI,MIS,MU,MUS,NSI,NSIS,PE,PES,PO,POS,XED1,
     &        XOD1
      REAL*8 M,MS
      COMPLEX*16 SUMB,SUME
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
C
C
C
C
C
C
      MAXL = ML - 1
      MAXLK = ML
      DO LK = 1,MAXLK
         DO NSI = 1,IREL
            ISI = (-1)*(-1)**NSI
            K = ISI*(LK-1) + (ISI-1)/2
            IF ( K.NE.0 .OR. IREL.NE.2 ) THEN
               IFE1 = 2*K*K - IDINT(ABS(DBLE(K))) + K
               IF ( IREL.EQ.1 ) IFE1 = K*K
C
               IF ( IREL.EQ.1 ) THEN
                  KAP1 = 2*K + 1
               ELSE IF ( IREL.EQ.2 ) THEN
                  KAP1 = 2*ABS(K)
               END IF
C
               MU = 1
               DO II1 = 1,KAP1
                  IF ( IREL.EQ.1 ) THEN
                     M = (-ABS(K)+II1-1)*1.0D0
                  ELSE IF ( IREL.EQ.2 ) THEN
                     M = -ABS(K) + 0.5D0 + (II1*1.D0-1.0D0)
                  END IF
C
                  PE = 0
                  PO = 0
                  SUMB = CZERO
                  SUME = CZERO
                  XED1 = XED/NAT
                  INDE1 = 1 + (IAN-1)*XED1
                  INDE2 = INDE1 + XED1 - 1
                  XOD1 = XOD/NAT
                  INDO1 = 1 + (IAN-1)*XOD1
                  INDO2 = INDO1 + XOD1 - 1
                  MI = IDINT(M)
                  IF ( IREL.EQ.2 ) THEN
                     DO IU = INDE1,INDE2
                        IF ( K-IDINT(KME(IU,1)).EQ.0 .AND. 
     &                       IDINT(M-KME(IU,2)).EQ.0 ) PE = IU
                        IF ( K-IDINT(KMO(IU,1)).EQ.0 .AND. 
     &                       IDINT(M-KMO(IU,2)).EQ.0 ) PO = IU
                     END DO
                  ELSE
                     DO IU = INDE1,INDE2
                        IF ( K-LME(IU,1).EQ.0 .AND. MI.EQ.LME(IU,2) )
     &                       PE = IU
                     END DO
                     DO IU = INDO1,INDO2
                        IF ( K-LMO(IU,1).EQ.0 .AND. MI.EQ.LMO(IU,2) )
     &                       PO = IU
                     END DO
                  END IF
                  DO IANS = 1,NAT
                     INDES1 = 1 + (IANS-1)*XED1
                     INDES2 = INDES1 + XED1 - 1
                     INDOS1 = 1 + (IANS-1)*XOD1
                     INDOS2 = INDOS1 + XOD1 - 1
                     DO LKS = 1,MAXLK
                        DO NSIS = 1,IREL
                           ISIS = (-1)*(-1)**NSIS
                           KS = ISIS*(LKS-1) + (ISIS-1)/2
                           IFE2 = 2*KS*KS - IDINT(ABS(DBLE(KS))) + KS
                           IF ( IREL.EQ.1 ) IFE2 = KS*KS
                           IF ( KS.NE.0 .OR. IREL.NE.2 ) THEN
C
                              IF ( IREL.EQ.1 ) THEN
                                 KAP2 = 2*KS + 1
                              ELSE IF ( IREL.EQ.2 ) THEN
                                 KAP2 = 2*ABS(KS)
                              END IF
C
                              MUS = 1
                              DO II2 = 1,KAP2
                                 IF ( IREL.EQ.1 ) THEN
                                    MS = (-ABS(KS)+II2-1)*1.0D0
                                 ELSE IF ( IREL.EQ.2 ) THEN
                                    MS = -ABS(KS) + 0.5D0 + 
     &                                 (II2*1.D0-1.0D0)
                                 END IF
C
                                 PES = 0
                                 POS = 0
                                 MIS = IDINT(MS)
                                 IF ( IREL.EQ.2 ) THEN
                                    DO JU = INDES1,INDES2
                                       IF ( KS-IDINT(KME(JU,1))
     &                                    .EQ.0 .AND. 
     &                                    IDINT(MS-KME(JU,2)).EQ.0 )
     &                                    PES = JU
                                       IF ( KS-IDINT(KMO(JU,1))
     &                                    .EQ.0 .AND. 
     &                                    IDINT(MS-KMO(JU,2)).EQ.0 )
     &                                    POS = JU
                                    END DO
                                 ELSE IF ( IREL.EQ.1 ) THEN
                                    DO JU = INDES1,INDES2
                                       IF ( KS-LME(JU,1).EQ.0 .AND. 
     &                                    MIS.EQ.LME(JU,2) ) PES = JU
                                    END DO
                                    DO JU = INDOS1,INDOS2
                                       IF ( KS-LMO(JU,1).EQ.0 .AND. 
     &                                    MIS.EQ.LMO(JU,2) ) POS = JU
                                    END DO
                                 END IF
C
                                 IF ( PES.NE.0 .AND. PE.NE.0 ) THEN
                                    SUMB = SUMB + BKM0(MUS+IFE2,IL,IANS)
     &                                 *XE(PES,PE)
                                    SUME = SUME + BKM0(MUS+IFE2,IL,IANS)
     &                                 *XE1(PES,PE)
                                 ELSE IF ( POS.NE.0 .AND. PO.NE.0 ) THEN
                                    SUMB = SUMB + BKM0(MUS+IFE2,IL,IANS)
     &                                 *XO(POS,PO)
                                    SUME = SUME + BKM0(MUS+IFE2,IL,IANS)
     &                                 *XO1(POS,PO)
                                 END IF
C
                                 IF ( IRUM.GT.0 .OR. BXY.GT.0 ) THEN
                                    IF ( POS.NE.0 .AND. PE.NE.0 ) THEN
                                       SUMB = SUMB + 
     &                                    BKM0(MUS+IFE2,IL,IANS)
     &                                    *XOE(POS,PE)
                                       SUME = SUME + 
     &                                    BKM0(MUS+IFE2,IL,IANS)
     &                                    *XOE1(POS,PE)
                                    ELSE IF ( PES.NE.0 .AND. PO.NE.0 )
     &                                 THEN
                                       SUMB = SUMB + 
     &                                    BKM0(MUS+IFE2,IL,IANS)
     &                                    *XEO(PES,PO)
                                       SUME = SUME + 
     &                                    BKM0(MUS+IFE2,IL,IANS)
     &                                    *XEO1(PES,PO)
                                    END IF
                                 END IF
C
                                 MUS = MUS + 1
                              END DO
                           END IF
                        END DO
                     END DO
                  END DO
                  BKM1(MU+IFE1,IL,IAN) = SUMB
                  BKM2(MU+IFE1,IL,IAN) = SUME
                  MU = MU + 1
               END DO
            END IF
         END DO
      END DO
C
      IF ( IP.GT.2 ) THEN
         WRITE (NOUT1,99001)
         IF ( IREL.EQ.1 ) THEN
            IWRITE = MAXLK**2
         ELSE
            IWRITE = 2*MAXLK**2
         END IF
         WRITE (NOUT1,99002) IAN,IL
         WRITE (NOUT1,99004) (BKM1(IO,IL,IAN),IO=1,IWRITE)
C
         WRITE (NOUT1,99003) IAN,IL
         WRITE (NOUT1,99004) (BKM2(IO,IL,IAN),IO=1,IWRITE)
      END IF
C
      RETURN
C
99001 FORMAT (1x,'bkmj',2x,'****')
99002 FORMAT (1x,'bkmj1 atom in position ',2x,i3,5x,'for layer',i3)
99003 FORMAT (1x,'bkmj2 atom in position ',2x,i3,5x,'for layer',i3)
99004 FORMAT (5(2x,e12.5,1x,e12.5))
      END
C*==intraf.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
      SUBROUTINE INTRAF(BKM2,RMATS,AWR,NATL,LAYS,LAYP,RO,ROL,LAYB,
     &                  LAYER1,IPOL,SPOL,ROSD,ROLSD)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,LL,PI,CZERO,CIMAG
      USE MOD_SPEC_GEOM,ONLY:TV
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPOL,LAYB,LAYER1,LAYP,LAYS,SPOL
      REAL*8 RO
      COMPLEX*16 AWR(MQD,LL,NATLM,2),BKM2(MQD,LL,NATLM),
     &           RMATS(MQD,MQD,LAYSM,NATLM),ROLSD(LL),ROSD(4)
      INTEGER NATL(LAYSM)
      REAL*8 ROL(LL)
C
C Local variables
C
      COMPLEX*16 ALAGE(:)
      INTEGER DREL(:),I,I1,I11,I2,I21,IAN,IBULK,IT,J,LU,NA,NB(:),SB(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NB,SB,DREL,ALAGE
      ALLOCATE (NB(MQD),SB(MQD),DREL(MQD),ALAGE(LL))
C
C*** End of declarations rewritten by SPAG
C
      CALL SORTBKM(NB,SB,DREL)
C
      RO = 0.0D0
      ROSD = CZERO
      DO I = 1,LL
         ALAGE(I) = CZERO
      END DO
C
      J = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
            IBULK = LAYB
         END IF
C
         DO IT = IBULK,LAYS
            NA = NATL(IT)
            DO IAN = 1,NA
               IF ( IPOL.EQ.1 ) THEN
                  DO I1 = 1,MQD
                     I11 = REL(NB(I1),SB(I1))
                     DO I2 = 1,MQD
                        I21 = REL(NB(I2),SB(I2))
                        ALAGE(J) = ALAGE(J) + AWR(I1,J,IAN,1)
     &                             *RMATS(I11,I21,IT,IAN)*BKM2(I2,J,IAN)
                     END DO
                  END DO
               ELSE IF ( IPOL.EQ.2 ) THEN
                  DO I1 = 1,MQD
                     I11 = REL(NB(I1),SB(I1))
                     DO I2 = 1,MQD
                        I21 = REL(NB(I2),SB(I2))
                        ALAGE(J) = ALAGE(J) + AWR(I1,J,IAN,2)
     &                             *RMATS(I11,I21,IT,IAN)*BKM2(I2,J,IAN)
                     END DO
                  END DO
               ELSE IF ( IPOL.EQ.3 ) THEN
                  DO I1 = 1,MQD
                     I11 = REL(NB(I1),SB(I1))
                     DO I2 = 1,MQD
                        I21 = REL(NB(I2),SB(I2))
                        ALAGE(J) = ALAGE(J) + AWR(I1,J,IAN,1)
     &                             *RMATS(I11,I21,IT,IAN)*BKM2(I2,J,IAN)
                     END DO
                  END DO
               ELSE IF ( IPOL.EQ.4 ) THEN
                  DO I1 = 1,MQD
                     I11 = REL(NB(I1),SB(I1))
                     DO I2 = 1,MQD
                        I21 = REL(NB(I2),SB(I2))
                        ALAGE(J) = ALAGE(J) + AWR(I1,J,IAN,2)
     &                             *RMATS(I11,I21,IT,IAN)*BKM2(I2,J,IAN)
                     END DO
                  END DO
               END IF
            END DO
            J = J + 1
         END DO
      END DO
C
      IF ( SPOL.LE.2 ) THEN
         DO J = 1,LAYP
            ROL(J) = DIMAG(CIMAG*ALAGE(J)/TV)/PI
            RO = RO + ROL(J)
         END DO
      ELSE IF ( SPOL.EQ.4 ) THEN
         DO J = 1,LAYP
            ROLSD(J) = CIMAG*ALAGE(J)/(TV*PI)
            ROSD = ROSD + ROLSD(J)
         END DO
      END IF
C
      END
C*==ajgs.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
C
      SUBROUTINE AJGS(IREL,BKM1,GANZ,AGSP,AGSM,LAYER1,K1,NATL,LAYS,
     &                EPOSLM,EPOSLP,KGT,LAYB)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MLQ,MLS,LL,LG,LH,PI,CZERO,CIMAG
      USE MOD_SPEC_GEOM,ONLY:TV
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
      INTEGER GANZ,IREL,LAYB,LAYER1,LAYS
      COMPLEX*16 K1
      COMPLEX*16 AGSM(LH,LL),AGSP(LH,LL),BKM1(MLS,LL,NATLM),
     &           EPOSLM(NATLM,LH,LAYSM),EPOSLP(NATLM,LH,LAYSM)
      REAL*8 KGT(2,LG)
      INTEGER NATL(LAYSM)
C
C Local variables
C
      REAL*8 CLEGOR
      REAL*8 G,KXX,KYY,M,S1
      INTEGER I,IAN,IBULK,IFE,IG,II1,IKZ,IL,IS,ISI,IT,IX,IY,J,K,KAP1,L,
     &        L1,LII,LK,LU,LV1,M1,M2,MAXL,MAXLK,MU,MU1,NAT,NSI
      COMPLEX*16 IC(:),KZ,KZZ(2),P,Q1,Q2,YLM(:),YP(:,:),YS(:,:),ZPI
      EXTERNAL CLEGOR
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
      ALLOCATABLE IC,YP,YS,YLM
      ALLOCATE (IC(ML),YP(ML,MLY),YS(ML,MLY),YLM(MLQ))
C
C*** End of declarations rewritten by SPAG
C
C
C
C     /****************************************************************/
C     subroutines and functions called from this routine               *
C     sphrm4    clegor                                                 *
C     /****************************************************************/
C     /* input */
C     /* output */
C     /* local */
C     /* common */
C     /* function */
C
      DO I = 1,LH
         DO J = 1,LL
            AGSP(I,J) = CZERO
            AGSM(I,J) = CZERO
         END DO
      END DO
      DO I = 1,ML
         IC(I) = CZERO
      END DO
C
      ZPI = DCMPLX(0.D0,2.D0*PI)/(K1*TV)
      MAXL = ML - 1
      MAXLK = ML
C
      IC(1) = CIMAG
      DO LV1 = 1,MAXL
         IC(LV1+1) = IC(LV1)*IC(1)
      END DO
C
      IL = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
            IBULK = LAYB
         END IF
         DO IT = IBULK,LAYS
            NAT = NATL(IT)
            IY = 1
            DO IG = 1,GANZ
               KXX = KGT(1,IG)
               KYY = KGT(2,IG)
               KZZ(1) = DCMPLX(-1.D0,0.D0)*CDSQRT(K1*K1-KXX*KXX-KYY*KYY)
               KZZ(2) = -KZZ(1)
               DO IKZ = 1,2
                  KZ = KZZ(IKZ)
                  CALL SPHRM4(KXX,KYY,KZ,YLM)
                  J = 1
                  DO L1 = 1,ML
                     M1 = 2*L1 - 1
                     DO M2 = 1,M1
                        IF ( IKZ.EQ.1 ) YS(L1,M2) = YLM(J)
                        IF ( IKZ.EQ.2 ) YP(L1,M2) = YLM(J)
                        J = J + 1
                     END DO
                  END DO
               END DO
               DO IS = 1,IREL
C
C        is=1,2 : spin = 1/2,-1/2
C
                  IX = IG
                  IF ( IS.EQ.2 ) IX = IG + GANZ
                  AGSP(IX,IL) = CZERO
                  AGSM(IX,IL) = CZERO
                  S1 = -0.5D0*(-1.D0)**IS
                  DO IAN = 1,NAT
                     DO LK = 1,MAXLK
                        DO NSI = 1,IREL
                           ISI = (-1)*(-1)**NSI
                           IF ( IREL.EQ.2 ) THEN
                              K = ISI*(LK-1) + (ISI-1)/2
                              IFE = 2*K*K - IABS(K) + K
                           ELSE
                              K = LK - 1
                              IFE = K*K
                           END IF
                           IF ( K.NE.0 .OR. IREL.NE.2 ) THEN
                              L = K
                              IF ( K.LT.0 ) L = ((-K)) - 1
                              P = -CIMAG
                              IF ( L.NE.0 ) P = P/IC(L)
C
                              IF ( IREL.EQ.1 ) THEN
                                 KAP1 = 2*K + 1
                              ELSE IF ( IREL.EQ.2 ) THEN
                                 KAP1 = 2*ABS(K)
                              END IF
C
                              MU = 1
                              DO II1 = 1,KAP1
                                 IF ( IREL.EQ.1 ) THEN
                                    M = (-ABS(K)+II1-1)*1.0D0
                                 ELSE IF ( IREL.EQ.2 ) THEN
                                    M = -ABS(K) + 0.5D0 + 
     &                                  (II1*1.D0-1.0D0)
                                 END IF
C
                                 G = 1.0D0
                                 IF ( IREL.EQ.2 )
     &                                G = CLEGOR(DBLE(K),M,S1)
                                 IF ( ABS(G).GT.0.D0 ) THEN
                                    LII = MU + IFE
                                    MU1 = MU
                                    IF ( IREL.EQ.2 ) THEN
                                       IF ( K.LT.0 .AND. IS.EQ.1 )
     &                                    MU1 = MU - 1
                                       IF ( K.GT.0 .AND. IS.EQ.2 )
     &                                    MU1 = MU + 1
                                    END IF
                                    Q1 = P*G*YP(L+1,MU1)*ZPI/KZZ(2)
                                    Q2 = P*G*YS(L+1,MU1)*ZPI/KZZ(2)
                                    AGSP(IX,IL) = AGSP(IX,IL)
     &                                 + Q1*BKM1(LII,IL,IAN)
     &                                 /EPOSLP(IAN,IY,IT)
                                    AGSM(IX,IL) = AGSM(IX,IL)
     &                                 + Q2*BKM1(LII,IL,IAN)
     &                                 /EPOSLM(IAN,IY,IT)
                                 END IF
                                 MU = MU + 1
                              END DO
                           END IF
                        END DO
                     END DO
                  END DO
                  IY = IY + 1
               END DO
            END DO
            IL = IL + 1
         END DO
      END DO
C
      END
C*==sphrm4.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
C
      SUBROUTINE SPHRM4(KXX,KYY,KZ,YLM)
C
      USE MOD_SPEC,ONLY:ML,MLQ,PI,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KXX,KYY
      COMPLEX*16 KZ
      COMPLEX*16 YLM(MLQ)
C
C Local variables
C
      REAL*8 A,ASG,B,CL,CM,FAC1(:),FAC2(:),FAC3(:),P
      COMPLEX*16 CF,CT,R,SA,SF,ST,XY
      INTEGER I,L,L1,LM,LM2,LM3,LMAX,LN,LO,LP,LQ,M
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE FAC1,FAC2,FAC3
      ALLOCATE (FAC1(ML),FAC2(MLQ),FAC3(ML))
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C
C
      DO I = 1,MLQ
         YLM(I) = CZERO
      END DO
C
      LMAX = ML - 1
      XY = DCMPLX(KXX*KXX+KYY*KYY,0.D0)
      R = CDSQRT(XY+KZ*KZ)
      P = SQRT(KXX*KXX+KYY*KYY)
      CT = KZ/R
      ST = DCMPLX(P,0.D0)/R
C
      IF ( P.LT.1.D-8 ) THEN
         CF = CONE
      ELSE
         CF = DCMPLX(KXX/P,KYY/P)
      END IF
C
      LM = 0
      CL = 0.D0
      A = 1.D0
      B = 1.D0
      ASG = 1.D0
      L1 = LMAX + 1
      DO L = 1,L1
         FAC1(L) = ASG*SQRT((2.D0*CL+1.D0)*A/(4.D0*PI*B*B))
         FAC3(L) = SQRT(2.0*CL)
         CM = -CL
         LN = L + L - 1
         DO M = 1,LN
            LO = LM + M
            FAC2(LO) = SQRT((CL+1.D0+CM)*(CL+1.D0-CM)/((2.D0*CL+3.D0)*(
     &                 2.D0*CL+1.D0)))
            CM = CM + 1.0D0
         END DO
         CL = CL + 1.D0
         A = A*2.0*CL*(2.0*CL-1.D0)/4.D0
         B = B*CL
         ASG = -ASG
         LM = LM + LN
      END DO
C
      LM = 1
      CL = 1.D0
      ASG = -1.D0
      SF = CF
      SA = CONE
      YLM(1) = DCMPLX(FAC1(1),0.D0)
      DO L = 1,LMAX
         LN = LM + L + L + 1
         YLM(LN) = FAC1(L+1)*SA*SF*ST
         YLM(LM+1) = ASG*FAC1(L+1)*SA*ST/SF
         YLM(LN-1) = -FAC3(L+1)*FAC1(L+1)*SA*SF*CT/CF
         YLM(LM+2) = ASG*FAC3(L+1)*FAC1(L+1)*SA*CT*CF/SF
         SA = ST*SA
         SF = SF*CF
         CL = CL + 1.D0
         ASG = -ASG
         LM = LN
      END DO
C
      LM = 1
      L1 = LMAX - 1
      DO L = 1,L1
         LN = L + L - 1
         LM2 = LM + LN + 4
         LM3 = LM - LN
         DO M = 1,LN
            LO = LM2 + M
            LP = LM3 + M
            LQ = LM + M + 1
            YLM(LO) = -(FAC2(LP)*YLM(LP)-CT*YLM(LQ))/FAC2(LQ)
         END DO
         LM = LM + L + L + 1
      END DO
C
      END
C*==bjgs.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
C
      SUBROUTINE BJGS(QD,QH3,QH4,QT3,QT4,AGSP,AGSM,BGSM,BULKX,LAYS,LAYP,
     &                LAYER1,LAYB,C1,C2,C3,C4,C5)
C
      USE MOD_SPEC,ONLY:LAYSM,LL,LH,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYB,LAYER1,LAYP,LAYS,QD
      COMPLEX*16 AGSM(LH,LL),AGSP(LH,LL),BGSM(LH,LL),BULKX(LH,LH,LAYSM),
     &           C1(QD,QD),C2(QD,QD),C3(QD,QD),C4(QD,QD),C5(QD,QD),
     &           QH3(LH,LH,LAYSM),QH4(LH,LH,LAYSM),QT3(LH,LH,LAYSM),
     &           QT4(LH,LH,LAYSM)
C
C Local variables
C
      INTEGER I,IBULK,ISEQ,J,K,L,LN,LSEQ
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     subroutines and functions called from this routine               *
C     copyq     mult      onemx                                        *
C     /****************************************************************/
C     /* input */
C
C     /* output */
C     /* local */
C     /* working arrays */
C
      DO I = 1,LH
         DO J = 1,LL
            BGSM(I,J) = CZERO
         END DO
      END DO
C
      CALL COPYQ(QH3,C1,QD,LAYS)
      CALL COPYQ(QT4,C2,QD,LAYS)
      CALL MULT(C1,C2,C3,QD)
      CALL COPYQ(QH4,C1,QD,LAYS)
      CALL MULT(C3,C1,C2,QD)
      CALL COPYQ(BULKX,C5,QD,LAYB)
      CALL MULT(C5,C2,C1,QD)
      CALL ONEMX(C1,C3,QD)
      CALL COPYQ(QT3,C3,QD,LAYS)
      CALL COPYQ(QH4,C2,QD,LAYS)
      CALL MULT(C3,C2,C4,QD)
      CALL MULT(C4,C1,C2,QD)
      CALL COPYQ(QH3,C1,QD,LAYS)
      CALL MULT(C5,C1,C3,QD)
      CALL MULT(C2,C3,C4,QD)
      DO J = 1,QD
         DO K = 1,QD
            BGSM(J,LAYP) = BGSM(J,LAYP) + C4(J,K)*AGSP(K,LAYP)
         END DO
         BGSM(J,LAYP) = BGSM(J,LAYP) + AGSM(J,LAYP)
      END DO
C
      L = 1
      DO LN = 1,LAYER1
         IF ( LN.LT.LAYER1 ) THEN
            IBULK = LAYB
         ELSE IF ( LN.EQ.LAYER1 ) THEN
            IBULK = 2
         END IF
         DO ISEQ = LAYS,IBULK, - 1
            LSEQ = ISEQ - 1
            IF ( LN.LT.LAYER1 .AND. ISEQ.EQ.IBULK ) LSEQ = LAYS
            CALL COPYQ(QH3,C1,QD,LSEQ)
            CALL COPYQ(QT4,C2,QD,LSEQ)
            CALL MULT(C1,C2,C3,QD)
            CALL COPYQ(QH4,C1,QD,LSEQ)
            CALL MULT(C3,C1,C2,QD)
            CALL COPYQ(BULKX,C5,QD,ISEQ)
            CALL MULT(C5,C2,C1,QD)
            CALL ONEMX(C1,C3,QD)
            CALL COPYQ(QH4,C2,QD,LSEQ)
            CALL COPYQ(QT3,C3,QD,LSEQ)
            CALL MULT(C3,C2,C4,QD)
            CALL MULT(C4,C1,C2,QD)
            CALL COPYQ(QH3,C1,QD,LSEQ)
            CALL MULT(C5,C1,C3,QD)
            CALL MULT(C2,C3,C4,QD)
            DO J = 1,QD
               DO K = 1,QD
                  BGSM(J,LAYP-L) = BGSM(J,LAYP-L) + C4(J,K)
     &                             *AGSP(K,LAYP-L) + C2(J,K)
     &                             *BGSM(K,LAYP+1-L)
               END DO
               BGSM(J,LAYP-L) = BGSM(J,LAYP-L) + AGSM(J,LAYP-L)
            END DO
            L = L + 1
         END DO
      END DO
C
      END
C*==coeffl.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
C
      SUBROUTINE COEFFL(KGZBLK,QD,WMHF,WMHG,IGVAL,IREL,IGP,GANZ,A1GP,
     &                  A1GPS,BGSM,TPM,TPP,TMP,TMM,PS,D0GM,D1GP,BULKX,
     &                  CMG,CMGS,D0GP,DETR,DETI,C1,C2,C3,C4,C5,C6,C7,
     &                  IBLOCH)
C
      USE MOD_SPEC,ONLY:LAYSM,LL,LH,MHDIM1,CZERO,CONE,CIMAG
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DETI,DETR
      INTEGER GANZ,IBLOCH,IGP,IREL,QD
      COMPLEX*16 A1GP(LH),A1GPS(LH),BGSM(LH,LL),BULKX(LH,LH,LAYSM),
     &           C1(QD,QD),C2(QD,QD),C3(QD,QD),C4(QD,QD),C5(QD,QD),
     &           C6(QD,QD),C7(QD,QD),CMG(LH),CMGS(LH),D0GM(LH),D0GP(LH),
     &           D1GP(LH),KGZBLK(LH),PS(LH,2),TMM(LH),TMP(LH),TPM(LH),
     &           TPP(LH),WMHF(MHDIM1,LH),WMHG(MHDIM1,LH)
      INTEGER IGVAL(LH)
C
C Local variables
C
      COMPLEX*16 BPM(:,:),BPP(:,:),DET,MMM(:,:),MMP(:,:),MPM(:,:),
     &           MPP(:,:)
      INTEGER I,IH,IK,J
C
C*** End of declarations rewritten by SPAG
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
C     purpose        calculates coefficients cg-, d0g+ which are used  *
C                    to couple whittaker functions to bulk solutions   *
C                    for energy e                                      *
C                                                                      *
C     subroutines and functions called from this routine               *
C     bmat      mult      onemx     copy      add       inv            *
C     mvmult                                                           *
C     /****************************************************************/
C     /* input */
C
C
C     /* output */
C     /* local */
C     /* working arrays */
C     /* common */
C
      DO I = 1,LH
         CMG(I) = CZERO
         CMGS(I) = CZERO
         D0GM(I) = CZERO
         D0GP(I) = CZERO
         D1GP(I) = CZERO
C        A1GP(I) = CZERO
C        A1GPS(I) = CZERO
      END DO
C
      IF ( IP.GE.3 ) THEN
         IF ( IREL.EQ.1 ) THEN
            WRITE (NOUT1,99001)
         ELSE
            WRITE (NOUT1,99002)
         END IF
         WRITE (NOUT1,99003) (A1GP(I),I=1,QD)
      END IF
      IF ( IGP.GT.0 ) THEN
C
C     calculating multiple scattering between surface-barrier and bulk
C
         CALL BMAT(BPP,BPM,MPM,MPP,MMP,MMM,TPM,TPP,TMP,TMM,PS,QD,2,
     &             IBLOCH)
C
         CALL COPYQ(BULKX,C6,QD,1)
C
         DO I = 1,QD
            DO J = 1,QD
               IF ( I.NE.J ) C6(I,J) = 1.00*C6(I,J)
            END DO
         END DO
C
         CALL MULT(BPM,C6,C1,QD)
         CALL MULT(C1,BPP,C2,QD)
         CALL MULT(C2,MPM,C1,QD)
         CALL ONEMX(C1,C2,QD)
         CALL MULT(C1,BPM,C2,QD)
         CALL MULT(C2,C6,C3,QD)
         CALL MULT(C3,BPP,C1,QD)
C
         DO I = 1,QD
            DO J = 1,QD
               C3(I,J) = CZERO
               C7(I,J) = CZERO
            END DO
         END DO
         DO I = 1,QD
            C3(I,I) = CONE
         END DO
C
         CALL MULT(BPP,MPM,C4,QD)
         CALL MULT(C4,BPM,C5,QD)
         CALL MULT(C5,C6,C4,QD)
         DO I = 1,QD
            DO J = 1,QD
               C5(I,J) = C3(I,J) - C4(I,J)
            END DO
         END DO
C
         CALL DETERMINANTE(C5,C3,C7,QD,DET)
         DETR = DREAL(DET)
         DETI = DIMAG(DET)
C
         DO I = 1,QD
            DO J = 1,QD
               D0GM(I) = D0GM(I) + C2(I,J)*BGSM(J,1) + C1(I,J)*A1GP(J)
            END DO
         END DO
C
         CALL COPY(BPP,C1,QD)
         CALL MULT(BPP,MPM,C2,QD)
         DO I = 1,QD
            DO J = 1,QD
               D1GP(I) = D1GP(I) + C2(I,J)*D0GM(J) + C1(I,J)*A1GP(J)
            END DO
         END DO
C
         DO I = 1,QD
            D0GP(I) = A1GP(I) + TPM(I)*D0GM(I)
         END DO
C
         IF ( IREL.EQ.1 ) THEN
            DO I = 1,IGP
               IH = IGVAL(I)
               CMGS(IH) = ((D0GP(IH)+D0GM(IH))*(WMHG(MHDIM1,IH)+A1GPS(IH
     &                    )))/(CIMAG*KGZBLK(IH)*(D0GP(IH)-D0GM(IH)))
     &                    - A1GP(IH) - WMHF(MHDIM1,IH)
               CMG(IH) = (D0GP(IH)+D0GM(IH))
     &                   /(A1GP(IH)+WMHF(MHDIM1,IH)+CMGS(IH))
            END DO
         ELSE IF ( IREL.EQ.2 ) THEN
            DO I = 1,IGP
               IH = IGVAL(I)
               IK = IH + GANZ
               CMGS(IH) = ((D0GP(IH)+D0GM(IH))*(WMHG(MHDIM1,IH)+A1GPS(IH
     &                    )))/(CIMAG*KGZBLK(IH)*(D0GP(IH)-D0GM(IH)))
     &                    - A1GP(IH) - WMHF(MHDIM1,IH)
               CMG(IH) = (D0GP(IH)+D0GM(IH))
     &                   /(A1GP(IH)+WMHF(MHDIM1,IH)+CMGS(IH))
               CMGS(IK) = ((D0GP(IK)+D0GM(IK))*(WMHG(MHDIM1,IK)+A1GPS(IK
     &                    )))/(CIMAG*KGZBLK(IK)*(D0GP(IK)-D0GM(IK)))
     &                    - A1GP(IK) - WMHF(MHDIM1,IK)
               CMG(IK) = (D0GP(IK)+D0GM(IK))
     &                   /(A1GP(IK)+WMHF(MHDIM1,IK)+CMGS(IK))
            END DO
         END IF
      ELSE IF ( IGP.EQ.0 ) THEN
      END IF
      RETURN
C
99001 FORMAT (1x,'a1g+')
99002 FORMAT (1x,'a1g+(+,-)')
99003 FORMAT (5(e12.5,e12.5))
      END
C*==djgs.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
C
C
      SUBROUTINE DJGS(QT2,QT4,BGSM,AGSP,QD,LAYER1,QH3,QH4,LAYS,DGSM,
     &                DGSP,BULKX,D1GP,LAYB,C1,C2,C3,C4,C5,C6)
C
      USE MOD_SPEC,ONLY:LAYSM,LL,LH,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYB,LAYER1,LAYS,QD
      COMPLEX*16 AGSP(LH,LL),BGSM(LH,LL),BULKX(LH,LH,LAYSM),C1(QD,QD),
     &           C2(QD,QD),C3(QD,QD),C4(QD,QD),C5(QD,QD),C6(QD,QD),
     &           D1GP(LH),DGSM(LH,LL),DGSP(LH,LL),QH3(LH,LH,LAYS),
     &           QH4(LH,LH,LAYSM),QT2(LH,LH,LAYSM),QT4(LH,LH,LAYSM)
C
C Local variables
C
      INTEGER I,IBULK,ISEQ,J,K,L,LN,LSEQ
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     subroutines and functions called from this routine               *
C     mult      copyq     onemx     copy                               *
C     /****************************************************************/
C     /* working arrays */
C
      DO I = 1,LH
         DO J = 1,LL
            DGSP(I,J) = CZERO
            DGSM(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,QD
         DGSP(I,1) = D1GP(I)
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
            CALL MULT(C3,C2,C4,QD)
            CALL COPYQ(QT2,C2,QD,ISEQ)
            CALL MULT(C1,C2,C3,QD)
            CALL COPYQ(BULKX,C6,QD,LSEQ)
            CALL MULT(C4,C6,C2,QD)
            CALL ONEMX(C2,C5,QD)
            CALL COPYQ(QH3,C1,QD,ISEQ)
            CALL MULT(C2,C1,C5,QD)
            CALL MULT(C2,C3,C1,QD)
            CALL MULT(C2,C4,C3,QD)
            DO J = 1,QD
               DO K = 1,QD
                  DGSP(J,L) = DGSP(J,L) + C5(J,K)*AGSP(K,L-1) + C1(J,K)
     &                        *DGSP(K,L-1) + C3(J,K)*BGSM(K,L)
               END DO
            END DO
            CALL COPYQ(QH4,C4,QD,ISEQ)
            CALL MULT(C4,C6,C3,QD)
            DO J = 1,QD
               DO K = 1,QD
                  DGSM(J,L-1) = DGSM(J,L-1) + C4(J,K)*BGSM(K,L)
     &                          + C3(J,K)*DGSP(K,L)
               END DO
            END DO
         END DO
      END DO
C
      END
C*==gkmu.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
      SUBROUTINE GKMU(IREL,MAXL,EMV,GANZ,CLIGHT,XE,XO,GKM2,DGSM,DGSP,
     &                LAYER1,NEV,NOD,IPOL,NATL,LAYS,EPOSLM,EPOSLP,XEO,
     &                XOE,IRUMP,NTOT,KGT,LAYB,IEDLAYER,IODLAYER,
     &                IEODLAYER,IOEDLAYER,BXY)
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
      INTEGER BXY,GANZ,IPOL,IREL,LAYB,LAYER1,LAYS,MAXL
      REAL*8 CLIGHT
      COMPLEX*16 EMV
      COMPLEX*16 DGSM(LH,LL),DGSP(LH,LL),EPOSLM(NATLM,LH,LAYSM),
     &           EPOSLP(NATLM,LH,LAYSM),GKM2(MLS,LL,NATLM),
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
      INTEGER IAN,IBULK,IIF,IKMU,IL,IRUM,IT,IW,LN,MAXLK,NAT,XED,XOD
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
C     subroutines and functions called from this routine               *
C     akmj0     akmj                                                   *
C     /****************************************************************/
C
C
C
      GKM2(1:MLS,1:LL,1:NATLM) = CZERO
      YP1(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YP2(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YS1(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
      YS2(1:ML,1:MLY,1:LX,1:NATLM) = CZERO
C
      MAXLK = MAXL + 1
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
               CALL AKMJ0(MAXL,EMV,GANZ,CLIGHT,IREL,DGSP,DGSM,IPOL,2,
     &                    IAN,YP1,YP2,YS1,YS2,EPOSLM,EPOSLP,IT,IL,KGT,
     &                    LH,LL,LX)
            END DO
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
C
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
               IF ( IREL.EQ.1 ) THEN
                  DO IKMU = 1,MAXLK*MAXLK
                     GKM2(IKMU,IL,IAN) = 4.D0*PI*AKM(IKMU,1,IL,IAN)
                  END DO
               ELSE IF ( IREL.EQ.2 ) THEN
                  DO IKMU = 1,2*MAXLK*MAXLK
                     GKM2(IKMU,IL,IAN)
     &                  = 4.D0*PI*(AKM(IKMU,1,IL,IAN)+AKM(IKMU,2,IL,IAN)
     &                  )
                  END DO
               END IF
               IF ( IP.GT.2 ) THEN
                  IF ( IREL.EQ.1 ) WRITE (NOUT1,99001) IAN,IL
                  IF ( IREL.EQ.2 ) WRITE (NOUT1,99002) IAN,IL
                  IF ( IREL.EQ.1 ) IW = MAXLK**2
                  IF ( IREL.EQ.2 ) IW = 2*MAXLK**2
                  WRITE (NOUT1,99003) (GKM2(IIF,IL,IAN),IIF=1,IW)
               END IF
            END DO
            IL = IL + 1
         END DO
      END DO
C
      RETURN
99001 FORMAT (1x,'glmj atom in position ',3x,i3,5x,'for layer',i3)
99002 FORMAT (1x,'gkmj atom in position ',3x,i3,5x,'for layer',i3)
99003 FORMAT (5(e12.5,e12.5))
      END
C*==sortbkm.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
C     /****************************************************************/
      SUBROUTINE SORTBKM(NB,SB,DREL)
C
      USE MOD_SPEC,ONLY:ML,MQD,MLS
      USE MOD_SPEC_RINDC,ONLY:REL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DREL(MQD),NB(MQD),SB(MQD)
C
C Local variables
C
      INTEGER A,C,D,E,F,I,IREL,J,S
C
C*** End of declarations rewritten by SPAG
C
C     /****************************************************************/
C     # purpose      :                                                 *
C                                                                      *
C     # subroutines called from this routine:                          *
C     /****************************************************************/
C
      E = 1
      C = 0
      A = 1
      S = 2
      D = 0
      F = 2*ML
      DO I = 1,F
         S = 3 - S
         DO J = 1,A
            D = D + 1
            IF ( D.GT.MLS ) GOTO 99999
            IF ( S.EQ.2 ) THEN
               C = C + 1
               NB(D) = C
            ELSE IF ( S.EQ.1 ) THEN
               IF ( J.LE.INT(A/2) ) NB(D) = E + INT(A/2+1)
               IF ( J.GT.INT(A/2) ) NB(D) = E - INT(A/2)
               E = E + 1
            END IF
            SB(D) = 3 - S
C           sb(d) = s
            IREL = REL(NB(D),SB(D))
            DREL(IREL) = D
         END DO
         A = A + 1
      END DO
C
99999 CONTINUE
      END
C*==determinante.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
      SUBROUTINE DETERMINANTE(A,B,C,QD,DET)
C     ==================================================================
C
C     calculates the determinant of a general complex matrix A
C
C     ==================================================================
      USE MOD_SPEC,ONLY:CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DET
      INTEGER QD
      COMPLEX*16 A(QD,QD),B(QD,QD),C(QD,QD)
C
C Local variables
C
      COMPLEX*16 EIG,EIGA(:),EIGB(:)
      INTEGER ITER(:),K2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE EIGA,EIGB,ITER
      ALLOCATE (EIGA(QD),EIGB(QD),ITER(QD))
C
      CALL LZHES(A,B,C,QD)
      CALL LZIT(A,B,C,QD,EIGA,EIGB,ITER)
C
      DET = CONE
C
      DO K2 = 1,QD
C
         IF ( CDABS(EIGB(K2)).LE.1D-20 ) THEN
            DET = DCMPLX(999999D0,9999999D0)
         ELSE IF ( CDABS(EIGA(K2)).LE.1D-20 ) THEN
            DET = DCMPLX(0D0,0D0)
         ELSE
            IF ( CDABS(EIGA(K2)**2).LE.1.D38*CDABS(EIGB(K2)**2) )
     &           EIG = EIGA(K2)/EIGB(K2)
            DET = DET*EIG
         END IF
C
      END DO
C
C
      END
