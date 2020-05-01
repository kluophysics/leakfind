C*==posanipot.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANIPOT
C   ********************************************************************
C   *                                                                  *
C   *  set up the pure Coulomb potential seen by a positron            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_RMESH,ONLY:NM,NRMAX,NSFMAX,NRSFMAX,FLMSF,R,JRNS1,ISFLM,
     &    DRDI,LMISF,KLMSF,NSF,JRCRI,JRCUT,JRWS,NPAN,FULLPOT,
     &    NMAX_LEBGRID,RHAT_LEBGRID,N_LEBGRID,W_LEBGRID,Y_LEBGRID
      USE MOD_TYPES,ONLY:NLMFPMAX,NLFPMAX,NTMAX,BNST,VNST,BT,VT,NLMFPT,
     &    KLMFP,IMT,VMTZ,RHO2NS,RHOCHR,RHOSPN,ITBOT,ITTOP,NTCLU,NTHOST,
     &    NLFP,NLMFP
      USE MOD_SITES,ONLY:NQHOST,ITOQ,NOQ
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,CONST_4PI,PI
      USE MOD_SCF,ONLY:SCFVXC,POSVC
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SYSTEM_TYPE
      IMPLICIT NONE
C*--POSANIPOT22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POSANIPOT')
C
C Local variables
C
      REAL*8 AGRDRHO(:),AGRDRHOD(:),AGRDRHOU(:),BLMTMP(:,:),BNEW(:,:,:),
     &       EXC(:),EXCIJ(:),FPIPR2,GDGAG(:),GDGAGD(:),GDGAGU(:),
     &       LAPRHO(:),LAPRHOD(:),LAPRHOU(:),RGNT_VSF(:),RHO(:),
     &       RHO4PI(:,:),RHOD(:),RHOU(:),RWGT,VAUX(:,:),VLMTMP(:,:),
     &       VLMXC,VLMXCAVR,VLMXCDIF,VLMXCDN,VLMXCUP,VNEW(:,:,:),
     &       VXC(:,:),WEXC(:,:),WEXCIJ(:,:),X,Y00
      REAL*8 DDOT
      LOGICAL GGA
      INTEGER I,IA_ERR,ICALL,IM,IO,IPAN1,IQ,IR,IRCRIT,IRMTIN,IRSF,IRTOP,
     &        ISF,ISPIN,IT,ITBOT_SAV,ITTOP_SAV,J,LM,LM1,LM2,LM3,
     &        LMRGNT_VSF(:,:),NRGNT_VSF,NRGNT_VSF_LM(:),NSPIN
      CHARACTER*30 SYSTEM_TYPE_SAV
      EXTERNAL DAXPY,DDOT
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
C------------------------------------------------------------------- ASA
      ALLOCATABLE RHO4PI,WEXC,EXC,RHO,RHOD,RHOU
      ALLOCATABLE VAUX
C--------------------------------------------------------------- FULLPOT
      ALLOCATABLE EXCIJ,WEXCIJ,VXC
      ALLOCATABLE LMRGNT_VSF,NRGNT_VSF_LM,RGNT_VSF
      ALLOCATABLE VLMTMP,BLMTMP,BNEW,VNEW
C------------------------------------------------------------------- GGA
      ALLOCATABLE AGRDRHO,AGRDRHOD,AGRDRHOU
      ALLOCATABLE LAPRHO,LAPRHOD,LAPRHOU
      ALLOCATABLE GDGAG,GDGAGD,GDGAGU
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C                      GGA - parametrisations
C=======================================================================
      IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
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
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
         WRITE (6,99002)
      ELSE
         WRITE (6,99003)
         RETURN
      END IF
C
      IF ( SYSTEM_DIMENSION(1:2).NE.'3D' )
     &      STOP 'in <POSANIPOT>: SYSTEM_DIMENSION <> 3D'
C
      ITBOT_SAV = ITBOT
      ITTOP_SAV = ITTOP
      SYSTEM_TYPE_SAV = SYSTEM_TYPE
C
      IF ( NTCLU.GT.0 ) THEN
C
         ITBOT = NTHOST
         ITTOP = 0
         DO IQ = 1,NQHOST
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               ITBOT = MIN(ITBOT,IT)
               ITTOP = MAX(ITTOP,IT)
            END DO
         END DO
         IF ( ITBOT.NE.1 ) STOP 'in <POSANIPOT>:  ITBOT <> 1'
C
         SYSTEM_TYPE = 'HOST      '
C
      END IF
C
C=======================================================================
C                                  ASA
C=======================================================================
C
      IF ( .NOT.FULLPOT ) THEN
C
         ALLOCATE (RHO4PI(NRMAX,2),WEXC(NRMAX,2),EXC(NRMAX))
         ALLOCATE (RHO(NRMAX),RHOD(NRMAX),RHOU(NRMAX))
         ALLOCATE (VAUX(NRMAX,2))
C
C---------------------------- reverse sign of TOTAL electronic potential
C
         VT(1:NRMAX,1:(NTHOST+NTCLU)) = -VT(1:NRMAX,1:(NTHOST+NTCLU))
         BT(1:NRMAX,1:(NTHOST+NTCLU)) = -BT(1:NRMAX,1:(NTHOST+NTCLU))
C
         DO IT = 1,(NTHOST+NTCLU)
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            DO IR = 1,IRTOP
               VAUX(IR,1) = VT(IR,IT) + BT(IR,IT)
               VAUX(IR,2) = VT(IR,IT) - BT(IR,IT)
               RHO4PI(IR,1) = RHOCHR(IR,IT)
               RHO4PI(IR,2) = RHOSPN(IR,IT)
               WEXC(IR,1) = 0D0
            END DO
C
C---------------------- remove electronic exchange-correlation potential
C
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
                  IF ( AGRDRHO(IR).GT.0.0 ) WRITE (6,*)
     &                  ' WARNING IN <SCFNEWPOT> D RHO/DR > 0.0 FOR R=',
     &                 RHO(IR),' IR=',IR
                  IF ( AGRDRHOU(IR).GT.0.0 ) WRITE (6,*) 
     &               ' WARNING IN <SCFNEWPOT> D RHO(UP)/DR > 0.0 FOR R='
     &               ,RHOU(IR),' IR=',IR
                  IF ( AGRDRHOD(IR).GT.0.0 ) WRITE (6,*) 
     &               ' WARNING IN <SCFNEWPOT> D RHO(DN)/DR > 0.0 FOR R='
     &               ,RHOD(IR),' IR=',IR
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
               IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
                  CALL EXCPBE(VAUX,EXC,WEXC,IRTOP,NRMAX,RHO,AGRDRHO,
     &                        LAPRHO,GDGAG,RHOU,AGRDRHOU,LAPRHOU,GDGAGU,
     &                        RHOD,AGRDRHOD,LAPRHOD,GDGAGD)
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
C------------------------- add electronic-positron correlation potential
C
            IF ( POSVC.EQ.'BM        ' )
     &           CALL VCEPBM(RHO4PI,VAUX,IRTOP,NRMAX)
C
C ---------------------- convention VAUX(I,IS) IS=1,2 == up, down (AKAI)
            DO IR = 1,IRTOP
               VT(IR,IT) = (VAUX(IR,1)+VAUX(IR,2))/2.0D0
C              BT(IR,IT) = (VAUX(IR,1)-VAUX(IR,2))/2.0D0
               BT(IR,IT) = 0D0
            END DO
C
         END DO
C
C---------------------------------------- muffintinize Coulomb potential
C
         CALL VMUFTIN(VMTZ)
C
C---------------------------- shift potential for embedded cluster atoms
         IF ( NTCLU.GT.0 ) THEN
C
            DO IT = ITTOP + 1,(NTHOST+NTCLU)
               IM = IMT(IT)
               IRTOP = JRWS(IM)
               VT(1:IRTOP,IT) = VT(1:IRTOP,IT) - VMTZ
               VT((IRTOP+1):NRMAX,IT) = 0.0D0
            END DO
C
         END IF
C
      ELSE
C=======================================================================
C                        FULL POTENTIAL
C=======================================================================
C
         ALLOCATE (EXCIJ(NMAX_LEBGRID),RHO4PI(NMAX_LEBGRID,2))
         ALLOCATE (WEXCIJ(NMAX_LEBGRID,2),VXC(NMAX_LEBGRID,2))
C
         IF ( IREL.GE.2 ) THEN
            NSPIN = 2
         ELSE
            NSPIN = 1
         END IF
C
C ------------------- NOTE: the same  LM-expansion is used for ALL sites
C
         IF ( NLFP.GT.NLFPMAX ) STOP ' in <POSANIPOT>:  NLFP > NLFPMAX'
C
         ALLOCATE (VNEW(NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (BNEW(NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (VLMTMP(NRSFMAX,NLMFPMAX),BLMTMP(NRSFMAX,NLMFPMAX))
C
         VNEW(1:NRMAX,1:NLMFPMAX,1:NTMAX) = 0D0
         BNEW(1:NRMAX,1:NLMFPMAX,1:NTMAX) = 0D0
C
C ----------------------------- construct angular mesh for xc-quantities
C
         CALL FPSPHERE
C
         DEALLOCATE (RHAT_LEBGRID)
C
C ------- initialize local variables  RGNT_VSF, NRGNT_VSF_LM, LMRGNT_VSF
C
         ALLOCATE (NRGNT_VSF_LM(0:NLMFPMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: POSANIPOT -> RGNT_VSF'
C
         CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
Cccccccccccccccccccc         CALL FPRGNT_VSF(NLFP,NRGNT_VSF_LM,NRGNT_VSF)
         CALL RGNTSF_SETUP(NLFP-1,NLFP-1,NRGNT_VSF_LM,NRGNT_VSF)
C
         ALLOCATE (RGNT_VSF(NRGNT_VSF),LMRGNT_VSF(NRGNT_VSF,3))
C
         REWIND IOTMP
         DO I = 1,NRGNT_VSF
            READ (IOTMP) (LMRGNT_VSF(I,J),J=1,3),RGNT_VSF(I)
         END DO
         CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C                check consistency of KLMSF and LMISF
C-----------------------------------------------------------------------
C
         DO IM = 1,NM
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               IF ( KLMSF(LM,IM).NE.1 ) THEN
                  WRITE (6,*) '   KLMSF <> 1 for IM=',IM,' LM =',LM
                  STOP ' in <POSANIPOT> '
               END IF
            END DO
         END DO
C
C ----------------------------------------------------------------------
C       construct the exchange correlation potential to  V  and  B
C ----------------------------------------------------------------------
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
         DO IT = 1,(NTHOST+NTCLU)
C
            IM = IMT(IT)
            IPAN1 = NPAN(IM)
            IRCRIT = JRCUT(IPAN1,IM)
            IRMTIN = JRCUT(1,IM)
C
C---> loop over radial mesh
C
            DO IR = 1,IRCRIT
C
C---> generate the densities on an angular mesh
C
               RHO4PI(1:N_LEBGRID,1:2) = 0.0D0
C
               FPIPR2 = CONST_4PI/R(IR,IM)**2
               DO ISPIN = 1,NSPIN
                  DO LM = 1,NLMFP
                     CALL DAXPY(N_LEBGRID,RHO2NS(IR,LM,IT,ISPIN)*FPIPR2,
     &                          Y_LEBGRID(1,LM),1,RHO4PI(1,ISPIN),1)
                  END DO
               END DO
C
C------------------------------------- to be consistent with EXCVWN etc.
               RHO4PI(1:N_LEBGRID,2) = -RHO4PI(1:N_LEBGRID,2)
C
               VXC(1:NMAX_LEBGRID,1:2) = 0.0D0
               EXCIJ(1:NMAX_LEBGRID) = 0.0D0
               WEXCIJ(1:NMAX_LEBGRID,1:2) = 0.0D0
C
C---> calculate the ex.-cor. potential
C
               IF ( SCFVXC.EQ.'NONE      ' ) THEN
C
C
               ELSE IF ( SCFVXC(1:3).EQ.'VBH' ) THEN
C
                  CALL EXCVBH(RHO4PI,VXC,EXCIJ,WEXCIJ,N_LEBGRID,
     &                        NMAX_LEBGRID)
C
               ELSE IF ( SCFVXC(1:3).EQ.'MJW' ) THEN
C
                  CALL EXCMJW(RHO4PI,VXC,EXCIJ,WEXCIJ,N_LEBGRID,
     &                        NMAX_LEBGRID,1)
C
               ELSE IF ( SCFVXC(1:3).EQ.'VWN' ) THEN
C
                  CALL EXCVWN(RHO4PI,VXC,EXCIJ,WEXCIJ,N_LEBGRID,
     &                        NMAX_LEBGRID)
C
               ELSE IF ( SCFVXC(1:3).EQ.'PBE' ) THEN
C
                  STOP
C
C            CALL EXCPBE(RHOAUX,VAUX,EXC,DRDI(1,IM),R(1,IM),NRMAX,JTOP,
C     &                  RHO,RHO1,RHO2,RHOU,RHOU1,RHOU2,RHOD,RHOD1,RHOD2)
C
               ELSE
                  WRITE (6,99001) SCFVXC
                  STOP
               END IF
C
C---> expand the ex.-cor. potential into spherical harmonics ,
C       using the orthogonality
C
               IF ( NSPIN.EQ.1 ) THEN
C
                  DO LM = 1,NLMFP
                     VLMXC = DDOT(N_LEBGRID,VXC(1,1),1,W_LEBGRID(1,LM),
     &                       1)
                     VNEW(IR,LM,IT) = VNEW(IR,LM,IT) + VLMXC
                  END DO
C
               ELSE
C
                  DO LM = 1,NLMFP
C
                     VLMXCDN = DDOT(N_LEBGRID,VXC(1,1),1,W_LEBGRID(1,LM)
     &                         ,1)
                     VLMXCUP = DDOT(N_LEBGRID,VXC(1,2),1,W_LEBGRID(1,LM)
     &                         ,1)
C
                     VLMXCAVR = 0.5D0*(VLMXCUP+VLMXCDN)
                     VLMXCDIF = 0.5D0*(VLMXCUP-VLMXCDN)
C
                     VNEW(IR,LM,IT) = VNEW(IR,LM,IT) + VLMXCAVR
                     BNEW(IR,LM,IT) = BNEW(IR,LM,IT) + VLMXCDIF
C
                  END DO
C
               END IF
C
            END DO
C
         END DO
         DO I = 1,JRCRI(1)
            WRITE (101,*) R(I,1),VNEW(I,1,1)
         END DO
C
C        Muffin Tin Shift done here: zero near cell boundary
C        next step: the electron potential is added, which has been
C        shifted in the same way, so afterwards no additional shift
C        should be needed.
C
         CALL FPVMUFTIN(NLMFP,VNEW,VMTZ)
C
C        CALL VMTZ_RCRI(VNEW,NRMAX,NLMFPMAX,VMTZ)
C
C-----------------------------------------------------------------------
C   convolute potential functions  B  and  V  with shape function
C-----------------------------------------------------------------------
C
         DO IT = 1,(NTHOST+NTCLU)
            IM = IMT(IT)
            IRMTIN = JRCUT(1,IM)
            IRCRIT = JRCRI(IM)
C
            DO LM = 1,NLMFP
               DO IR = 1,IRCRIT - IRMTIN
                  VLMTMP(IR,LM) = 0.0D0
                  BLMTMP(IR,LM) = 0.0D0
               END DO
            END DO
C
            DO J = 1,NRGNT_VSF_LM(NLMFP)
               LM1 = LMRGNT_VSF(J,1)
               IF ( LM1.GT.NLMFPMAX )
     &               STOP ' in <POSANIPOT>:  LMRGNT_VSF > NLMFPMAX'
               LM2 = LMRGNT_VSF(J,2)
               LM3 = LMRGNT_VSF(J,3)
               IF ( KLMSF(LM3,IM).EQ.1 ) THEN
                  ISF = ISFLM(LM3,IM)
C ----------------------------------------------------------------------
                  IF ( ISF.LE.NSFMAX .AND. ISF.GE.1 ) THEN
C ----------------------------------------------------------------------
                     DO IR = IRMTIN + 1,IRCRIT
                        IRSF = IR - IRMTIN
                        RWGT = RGNT_VSF(J)*FLMSF(IRSF,ISF,IM)
C
                        VLMTMP(IRSF,LM1) = VLMTMP(IRSF,LM1)
     &                     + RWGT*VNEW(IR,LM2,IT)
                        BLMTMP(IRSF,LM1) = BLMTMP(IRSF,LM1)
     &                     + RWGT*BNEW(IR,LM2,IT)
                     END DO
C ----------------------------------------------------------------------
                  ELSE
                     WRITE (6,*) ' ISF: ',ISF,J,LM1,LM2,LM3,IM,NSFMAX
                     STOP ' in <POSANIPOT>:  ISF out of range  !!'
                  END IF
C ----------------------------------------------------------------------
               END IF
            END DO
C ----------------------------------------------------------------------
C
            DO LM = 1,NLMFP
               DO IR = IRMTIN + 1,IRCRIT
                  IRSF = IR - IRMTIN
                  VNEW(IR,LM,IT) = VLMTMP(IRSF,LM)
                  BNEW(IR,LM,IT) = BLMTMP(IRSF,LM)
               END DO
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C    subtract original potential  (V,B)  from XC-potential (VNEW,BNEW)
C    to get pure Coulomb potential with sign reversed as seen by e^(+)
C-----------------------------------------------------------------------
C
         Y00 = 1D0/SQRT_4PI
C
         DO IT = 1,(NTHOST+NTCLU)
C
            IM = IMT(IT)
            DO IR = 1,JRCRI(IM)
               VNEW(IR,1,IT) = VNEW(IR,1,IT) - VT(IR,IT)/Y00
               BNEW(IR,1,IT) = BNEW(IR,1,IT) - BT(IR,IT)/Y00
            END DO
C
            DO LM = 2,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 .AND. LM.LE.NLMFPMAX ) THEN
                  DO IR = JRNS1(IM),JRCRI(IM)
                     VNEW(IR,LM,IT) = VNEW(IR,LM,IT) - VNST(IR,LM,IT)
                     BNEW(IR,LM,IT) = BNEW(IR,LM,IT) - BNST(IR,LM,IT)
                  END DO
               END IF
            END DO
C
         END DO
C
C-----------------------------------------------------------------------
C     determine muffin tin zero and shift potential  V  accordingly
C-----------------------------------------------------------------------
C
         !write(6,*) 'vmtz was: ', vmtz
         !write(6,*) 'fpvmuftin suppressed in posanipot.f'
         !vmtz=0d0
C
C---------------------------- shift potential for embedded cluster atoms
C --------------------------MKO bugfix: better don't do that...
C         IF ( NTCLU.GT.0 ) THEN
CC
C            VCORR = SQRT_4PI*VMTZ
C            DO IT = ITTOP + 1,(NTHOST+NTCLU)
C               IM = IMT(IT)
C               IRTOP = JRCUT(NPAN(IM),IM)
C               VNEW(1:IRTOP,1,IT) = VNEW(1:IRTOP,1,IT) - VCORR
C            END DO
CC
C         END IF
C
C-----------------------------------------------------------------------
C    convert to the standard SPR-KKR potential functions
C    VT    BT     spherical part INCLUDING Y00
C    VNST  BNST   NON-spherical parts for LM>1
C-----------------------------------------------------------------------
C
         Y00 = 1D0/SQRT_4PI
C
         DO IT = 1,(NTHOST+NTCLU)
C
            IM = IMT(IT)
            DO IR = 1,JRCRI(IM)
               VT(IR,IT) = VNEW(IR,1,IT)*Y00
               BT(IR,IT) = BNEW(IR,1,IT)*Y00
            END DO
C
            DO LM = 2,NLMFPT(IT)
               IF ( KLMFP(LM,IT).NE.0 .AND. LM.LE.NLMFPMAX ) THEN
                  DO IR = JRNS1(IM),JRCRI(IM)
                     VNST(IR,LM,IT) = VNEW(IR,LM,IT)
                     BNST(IR,LM,IT) = BNEW(IR,LM,IT)
                  END DO
               END IF
            END DO
C
         END DO
C
         DEALLOCATE (VLMTMP,BLMTMP)
C
      END IF
C=======================================================================
C
C------------------------------------------------- reset ITBOT and ITTOP
      ITBOT = ITBOT_SAV
      ITTOP = ITTOP_SAV
      SYSTEM_TYPE = SYSTEM_TYPE_SAV
C
C=======================================================================
99001 FORMAT (//,1X,79('#'),/,10X,'ERROR in <POSANIPOT> ',/,10X,
     &        'SCFVXC = ',A,'  not known ',/,1X,79('#'),/)
99002 FORMAT (/,1X,79('*'),/,18X,
     &        'setting up positron potential in <POSANIPOT>',/,1X,
     &        79('*'),/)
99003 FORMAT (/,1X,79('*'),/,13X,'positron potential already set up',
     &        ' -- leaving <POSANIPOT>',/,1X,79('*'),/)
      END
