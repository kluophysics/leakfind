C*==core_part.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     eigen.f
C
C     new subroutine savecore has to be made
C            writecore has a bug if the states are rotatet
C     calcore has to be changed for arbitrary B-field
C            such that rotation is principally not longer
C            necesary for in-plane magnetization
C
C     contains the subroutines:
C         core_part
C         -   calcore
C             -  mnewt
C                -   usrfun
C                    -   inward
C                    -   outward
C                -   ludcmp
C                -   lubksb
C                -  make_cwf
C             -  writecore
C             -  savecore NOT YET DONE
C         -   readcore
C             -  detnzrowcol
C             -  multsparse
C
C     and functions:
C         ctp
C         frelp
C
C     /****************************************************************/
C
C
      SUBROUTINE CORE_PART(Z,CLIGHT,RN,VLM,ALPHA,LAYS,NATL,EB,ESTART,
     &                     EEND,EVAL,CWFX,CWFY,CWF,CWFSTATES,CWFENERGY,
     &                     CWFKAPMUE,ROT,ROTM1,CORE)
C     /****************************************************************/
C     # purpose      : calculate core wavefunctions                    *
C                      / read wavefunctions from file                  *
C                                                                      *
C     # parameter    :                                                 *
C       input:                                                         *
C       ======                                                         *
C       z -- atomic number information                                 *
C       clight -- speed of light                                       *
C       rn -- atomic sphere radii                                      *
C       vlm -- atomic potentials                                       *
C       alpha -- parameter for the exponential radial mesh             *
C       lays -- number of layers in material                           *
C       natl -- number of atoms per layer                              *
C       eb -- magnetic field vector                                    *
C       estart -- starting point for determination of core states      *
C       eend -- ending point for determination of core states          *
C       eua, eub, euc -- euler angles which specify the rotation of    *
C                        the magnetic field vector                     *
C       rot -- rotation matrix to perform the rotation of wavefunctions*
C              according to the magnetic field                         *
C       rotm1 -- the inverse of the roatation matrix                   *
C       core -- flag to specify action:                                *
C               core = 1 -- calculate wavefunctions                    *
C               core = 2 -- read wavefunctions from file               *
C       output:                                                        *
C       =======                                                        *
C       eval -- eigenvalues of core-states                             *
C       cwfx -- indexfield for core-wavefunction array                 *
C       cwfy -- indexfield for core-wavefunction array                 *
C       cwf -- core-wavefunction array                                 *
C       cwfstates -- number of core-states for each atom               *
C       cwfenergy -- core energies for every atom                      *
C                       (same as eval but different order)             *
C       cwfkapmue -- kappa, mue information for every core-state       *
C     # calls the functions:                                           *
C       calcore       readcore                                         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MAXEIG,MAXN,MAXCORE,
     &    MAXCSTATES,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CLIGHT,EEND,ESTART
      INTEGER CORE,LAYS
      COMPLEX*16 CWF(RSTEP,MAXCORE,2),ROT(MQD,MQD),ROTM1(MQD,MQD),
     &           VLM(LAYSM,NATLM,MQD,RSTEP)
      REAL*8 CWFENERGY(MAXCSTATES,NATLM,LAYSM),
     &       CWFKAPMUE(MAXCSTATES,NATLM,LAYSM,2),EB(3),
     &       EVAL(MAXEIG,MAXN,NATLM,LAYSM),RN(NATLM,LAYSM),
     &       Z(NATLM,LAYSM)
      INTEGER CWFSTATES(NATLM,LAYSM),CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),
     &        CWFY(MAXCORE),NATL(LAYSM)
C
C Local variables
C
      INTEGER CHANNEL,I,J,K,L
      CHARACTER*7 CHNAME
      LOGICAL EXISTS
C
C*** End of declarations rewritten by SPAG
C
C
C     initialize arrays
      DO I = 1,MQD
         DO J = 1,MAXCSTATES
            DO K = 1,NATLM
               DO L = 1,LAYSM
                  CWFX(I,J,K,L,1) = 0
                  CWFX(I,J,K,L,2) = -1
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,MQD
         CWFY(I) = 0
      END DO
C
      DO I = 1,NATLM
         DO J = 1,LAYSM
            CWFSTATES(I,J) = 0
         END DO
      END DO
C
      DO K = 1,MAXCSTATES
         DO I = 1,NATLM
            DO J = 1,LAYSM
               CWFENERGY(K,I,J) = 0.D0
               CWFKAPMUE(K,I,J,1) = 0.D0
               CWFKAPMUE(K,I,J,2) = 0.D0
            END DO
         END DO
      END DO
C
      DO I = 1,RSTEP
         DO J = 1,MAXCORE
            CWF(I,J,1) = CZERO
            CWF(I,J,2) = CZERO
         END DO
      END DO
C
      CHANNEL = 34
      CHNAME = 'fort.34'
C
      IF ( CORE.EQ.1 ) THEN
C       calculate core wavefunctions:
C       =============================
C         open file to write the wave functions
         INQUIRE (FILE=CHNAME,EXIST=EXISTS)
         IF ( .NOT.EXISTS ) THEN
            OPEN (UNIT=CHANNEL,FILE=CHNAME,STATUS='new')
         ELSE
            OPEN (UNIT=CHANNEL,FILE=CHNAME,STATUS='replace')
         END IF
C
         IF ( IP.GT.0 ) WRITE (*,99001)
         CALL CALCORE(Z,CLIGHT,RN,VLM,ALPHA,LAYS,NATL,EB,ESTART,EEND,
     &                EVAL,CHANNEL)
         CLOSE (CHANNEL)
      END IF
      IF ( IP.GT.1 ) WRITE (*,99002)
C
C     read core wavefunctions and energies from file:
C     ===============================================
C     open file to write the wave functions
      INQUIRE (FILE=CHNAME,EXIST=EXISTS)
      IF ( .NOT.EXISTS ) THEN
         WRITE (*,99004)
         STOP
      ELSE
         OPEN (UNIT=CHANNEL,FILE=CHNAME,STATUS='old')
      END IF
      IF ( IP.GT.1 ) WRITE (*,99003)
      CALL READCORE(CWFX,CWFY,CWF,EVAL,LAYS,NATL,CWFSTATES,CWFENERGY,
     &              CWFKAPMUE,ROT,ROTM1,EB,CHANNEL)
      CLOSE (CHANNEL)
C
      RETURN
C
99001 FORMAT ('calculating core wavefunctions')
99002 FORMAT ('core wavefunctions calculated')
99003 FORMAT ('reading core wavefunctions from file')
99004 FORMAT ('stop in corepart: fort.34 does not exist !')
      END
C*==calcore.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CALCORE(Z,CLIGHT,RN,VLM,ALPHA,LAYS,NATL,EB,ESTART,EEND,
     &                   EVAL,CHANNEL)
C     /****************************************************************/
C     # purpose      : determine core wavefunctions and energies       *
C                      for the given potentials vlm                    *
C                                                                      *
C     # parameter:                                                     *
C     # input    :                                                     *
C       z        -- atomic number information                          *
C       clight   -- speed of light                                     *
C       rn       -- atomic sphere radii                                *
C       vlm      -- atomic potentials                                  *
C       alpha    -- parameter for the exponential radial mesh          *
C       lays     -- number of layers in material                       *
C       natl     -- number of atoms per layer                          *
C       eb       -- magnetic field vector                              *
C       estart   -- starting point for determination of core states    *
C       eend     -- ending point for determination of core states      *
C       channel  -- i/o-channel for wave functions                     *
C     # on return:                                                     *
C       ==========                                                     *
C       no return values are specified. the wavefunctions are written  *
C       to file. from this they are re-read in subroutine readcore.    *
C                                                                      *
C       strategy :  the eigenvalue problem is solved using the newton  *
C       ==========  raphson algorithm as specified by ebert            *
C                   (j. phys. cond. matter 1, 9111, 1989).             *
C                   the algotrithm can be found in "numerical recipes".*
C                   the eigenvalues and wavefunctions will be written  *
C                   to the file specified by channel.                  *
C                                                                      *
C     # calls the subroutines and functions:                           *
C       rmesh     rmeshda     mnewt     writecore                      *
C       kapmue2k  kapmue2l   frelp      ctp                            *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,RSTEPP5,MAXEIG,MLXPS,
     &    MAXN,HARTRE,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CLIGHT,EEND,ESTART
      INTEGER CHANNEL,LAYS
      REAL*8 EB(3),EVAL(MAXEIG,MAXN,NATLM,LAYSM),RN(NATLM,LAYSM),
     &       Z(NATLM,LAYSM)
      INTEGER NATL(LAYSM)
      COMPLEX*16 VLM(LAYSM,NATLM,MQD,RSTEP)
C
C Local variables
C
      REAL*8 ASTART,B0,BR(:),BZ,DE,EA,EMAG(2),EREL,LAST,MU,NEXT,POT(:),
     &       RAD(:),RADP4(:),V00,VAL,VR(:),VR1(:,:,:),Y00,ZZ
      INTEGER ATOM,FITAT,FLAG,I,IE,IKM,IKP,IMU,KK,LAY,LLC,MUSTART,
     &        MUSTOP,NE,NN,PTP
      INTEGER CTP,KAPMUE2K,KAPMUE2L
      COMPLEX*16 CWFFS(:,:,:),CWFGS(:,:,:)
      REAL*8 FRELP
      EXTERNAL CTP,FRELP,KAPMUE2K,KAPMUE2L
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE BR,VR,VR1,RAD,POT,RADP4,CWFFS,CWFGS
      ALLOCATE (BR(RSTEP),VR(RSTEP),VR1(RSTEP,NATLM,LAYSM))
      ALLOCATE (RAD(RSTEP),POT(RSTEP),RADP4(RSTEPP5))
      ALLOCATE (CWFFS(RSTEP,2,2),CWFGS(RSTEP,2,2))
C
C*** End of declarations rewritten by SPAG
C
      Y00 = 1.D0/DSQRT(4.D0*PI)
      BZ = EB(3)
C
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
C         prepare data for atom atom in layer lay:
C         ========================================
            ZZ = Z(ATOM,LAY)
C            IF ( ZZ.NE.0.D0 ) THEN
            IF ( ABS(ZZ).GT.1.0D-16 ) THEN
C
               CALL RMESH(RAD,RN,ALPHA,RSTEP,ATOM,LAY)
               CALL RMESHDA(RADP4,RN,ALPHA,RSTEP,4,ATOM,LAY)
C
               IF ( ML.GT.1 ) THEN
                  IKM = KAPMUE2K(-1,-0.5D0,ML)
                  IKP = KAPMUE2K(-1,+0.5D0,ML)
               ELSE
                  WRITE (*,*) 'stop because of ml .le.1 !!!'
                  STOP
               END IF
C
               V00 = DBLE(VLM(LAY,ATOM,IKM,1)+VLM(LAY,ATOM,IKP,1))/2.D0
               B0 = BZ*DBLE(VLM(LAY,ATOM,IKM,1)-VLM(LAY,ATOM,IKP,1))
     &              /2.D0
               DO I = 1,RSTEP
                  VAL = RAD(I)*Y00
                  VR(I) = DBLE(VLM(LAY,ATOM,IKM,I)+VLM(LAY,ATOM,IKP,I))
     &                    *VAL/2.D0
                  BR(I) = BZ*DBLE(VLM(LAY,ATOM,IKM,I)-VLM(LAY,ATOM,IKP,I
     &                    ))*VAL/2.D0
                  POT(I) = Y00*DBLE(VLM(LAY,ATOM,IKP,I)+VLM(LAY,ATOM,IKM
     &                     ,I))/2.D0
                  VR1(I,ATOM,LAY) = POT(I)*RAD(I) + ZZ
               END DO
C
               DO LLC = 0,MLXPS - 1    ! usually s to f, see parms.h
                  DO KK = LLC, - LLC - 1, - (2*LLC+1)
                     IF ( KK.NE.0 ) THEN
                        NN = 0
                        MU = ABS(KK) - 0.5D0
C
                        FITAT = CTP(VLM,ESTART,LAY,ATOM)
C                        IF ( FITAT.NE.0 .AND. ESTART.NE.0.D0 ) THEN
                        IF ( FITAT.NE.0 .AND. ABS(ESTART).GT.1.0D-16 )
     &                       THEN
C
                           LAST = FRELP(ESTART,KK,VR1,RAD,CLIGHT,RSTEP,
     &                            ALPHA,Z,ATOM,LAY,FITAT)
C
                           DE = 0.5D0/HARTRE
                           ASTART = MIN(ESTART,EEND)
                           NE = 1 + INT(ABS((EEND-ESTART)/DE))
                           DO IE = 0,NE
                              EA = ASTART + IE*DE
                              FITAT = CTP(VLM,EA,LAY,ATOM)
C                              IF ( FITAT.EQ.0 .OR. EA.EQ.0.D0 ) THEN
                              IF ( FITAT.EQ.0 .OR. ABS(EA).LE.1.0D-16 )
     &                             THEN
                                 LAST = NEXT
                                 CYCLE
                              END IF
C
                              NEXT = FRELP(EA,KK,VR1,RAD,CLIGHT,RSTEP,
     &                               ALPHA,Z,ATOM,LAY,FITAT)
C
                              IF ( LAST*NEXT.LT.0.D0 ) THEN
                                 NN = NN + 1
                                 MUSTART = -(2*ABS(KK)-1)
                                 MUSTOP = 2*ABS(KK) - 1
                                 DO IMU = MUSTART,MUSTOP,2
                                    MU = DBLE(IMU)*0.5D0
                                    CALL MNEWT(EA,ZZ,CLIGHT,KK,MU,V00,
     &                                 B0,VR,BR,RAD,RADP4,LAY,ATOM,VLM,
     &                                 ALPHA,FLAG,EREL,EMAG,CWFFS,CWFGS)
C
C                                    IF ( EREL.NE.0.D0 ) THEN
                                    IF ( ABS(EREL).GT.1.0D-16 ) THEN
                                       PTP = KAPMUE2L(KK,DBLE(MU))
                                       EVAL(PTP,NN,ATOM,LAY) = EREL
C                                 write core values to file:
C                                 ==========================
                                       CALL WRITECORE(EREL,EMAG,KK,MU,
     &                                    NN,ATOM,LAY,CHANNEL)
                                    END IF
                                 END DO
                              END IF
                              LAST = NEXT
                           END DO
                        END IF
                     END IF
C             end kk:
C             =======
                  END DO
C         end llc:
C         ========
               END DO
            END IF
C     end atom:
C     =========
         END DO
C     end lay:
C     ========
      END DO
C
      END
C*==ctp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION CTP(VLM,EREL,LAY,ATOM)
C     /****************************************************************/
C     # purpose      : determine classical turning point for erel      *
C                      and potential vlm                               *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,LAY
      REAL*8 EREL
      COMPLEX*16 VLM(LAYSM,NATLM,MQD,RSTEP)
C
C Local variables
C
      INTEGER COUNTVAR,IKM
      REAL*8 ETEST
      INTEGER KAPMUE2K
      EXTERNAL KAPMUE2K
C
C*** End of declarations rewritten by SPAG
C
      IKM = KAPMUE2K(-1,-0.5D0,ML)
      CTP = RSTEP
      DO COUNTVAR = RSTEP,1, - 1
         ETEST = DBLE(VLM(LAY,ATOM,IKM,COUNTVAR)) - EREL
         IF ( ETEST.LE.0.D0 ) THEN
            CTP = COUNTVAR + 1
            EXIT
         END IF
      END DO
C
      IF ( CTP.GE.RSTEP ) THEN
         CTP = RSTEP - 5
         WRITE (*,99001) CTP
         WRITE (NOUT1,99001) CTP
      END IF
C
      IF ( CTP.LE.5 ) THEN
         WRITE (*,99002) CTP
         WRITE (NOUT1,99002) CTP
      END IF
C
      RETURN
C
99001 FORMAT ('warning in ctp: ctp >= rstep = max(rm) set ctp=',i2)
99002 FORMAT ('warning in ctp: ctp=',i3,' <= 5 = min(rm)')
      END
C*==frelp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      REAL*8 FUNCTION FRELP(E,KAPPA,VR1,RAD,CLIGHT,NMESH,ALPHA,Z,ATOM,
     &                      LAY,FITAT)
C     /****************************************************************/
C     # purpose      : calculate coefficient f1                        *
C                      for relativistic paramagnetic case              *
C                                                                      *
C     # calls the subroutines:                                         *
C       start0    hank1                                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,RSTEP,MLP,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CLIGHT,E
      INTEGER ATOM,FITAT,KAPPA,LAY,NMESH
      REAL*8 RAD(RSTEP),VR1(RSTEP,NATLM,LAYSM),Z(NATLM,LAYSM)
C
C Local variables
C
      COMPLEX*16 AI1(RSTEP),BI1(RSTEP),DF,DK,DKR,EC,WFF,WFG,XH1(MLP),
     &           XJ(MLP),XN(MLP)
      REAL*8 DG,DKAPPA,F1,G,KK,P1,P2,R(2),VO,XRL
      INTEGER I,IREL
C
C*** End of declarations rewritten by SPAG
C
C     irel is dummy for start0, to stay compatible to solvedirac
      IREL = 2
      EC = CONE*E
      DK = CDSQRT(2.D0*EC+(EC/CLIGHT)**2)
      XRL = 1.D0
      DKR = DK*RAD(NMESH)
C
      CALL HANK1(DKR,XJ,XN,XH1)
C
      VO = VR1(1,ATOM,LAY)/RAD(1)
      KK = SQRT(-2.D0*E)
C      F1 = 0.D0
      DKAPPA = DBLE(KAPPA)
C
      CALL START0(Z,NMESH,ALPHA,EC,VO,CLIGHT,XRL,DKAPPA,RAD,VR1,AI1,BI1,
     &            IREL,ATOM,LAY)
C
      G = SQRT(DABS(DKAPPA**2-(Z(ATOM,LAY)/CLIGHT)**2))
C
      DO I = 1,RSTEP
         WFF = AI1(FITAT)*RAD(FITAT)**G
         WFG = BI1(FITAT)*RAD(FITAT)**G
      END DO
C
C     determine derivative of ffs and fgs at rad(fitat):
C     ==================================================
      P1 = (E+2.D0*CLIGHT**2-(VR1(FITAT,ATOM,LAY)-Z(ATOM,LAY))
     &     /RAD(FITAT))/CLIGHT
      P2 = (E-(VR1(FITAT,ATOM,LAY)-Z(ATOM,LAY))/RAD(FITAT))/CLIGHT
C
C     f = 1/r * ffs:
C     ==============
      DF = (-KAPPA)/RAD(FITAT)*WFF + P1*WFG
      DF = -WFF/RAD(FITAT)**2 + DF/RAD(FITAT)
C
C     g = 1/r * fgs:
C     ==============
      DG = REAL(KAPPA/RAD(FITAT)*WFG-P2*WFF)
      DG = REAL(-WFG/RAD(FITAT)**2+DG/RAD(FITAT))
C
C      IF ( DBLE(WFF).NE.0.D0 ) THEN
      IF ( ABS(DBLE(WFF)).GT.1.0D-16 ) THEN
         R(1) = DBLE(WFF/RAD(FITAT))
         R(2) = DBLE(DF)
      ELSE
         R(1) = DIMAG(WFF/RAD(FITAT))
         R(2) = DIMAG(DF)
      END IF
C
      F1 = (R(1)*KK+R(2))/(2.D0*EXP(KK*RAD(FITAT))*KK)
C
      FRELP = F1
C
      END
C*==mnewt.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE MNEWT(EIN,ZZ,CLIGHT,KK,MU,V00,B0,VR,BR,RAD,RADP4,LAY,
     &                 ATOM,VLM,ALPHA,FLAG,EREL,EMAG,CWFFS,CWFGS)
C     /****************************************************************/
C     # purpose      : solve eigenvalue problem using                  *
C                      newton raphson algorithm                        *
C                                                                      *
C     # calls the subroutines:                                         *
C       usrfun    ludcmp      lubksb                                   *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,RSTEPP5,EPS12,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NTRIAL,N,NP,NSOL
      PARAMETER (NTRIAL=150,N=4,NP=4,NSOL=2)
      REAL*8 TOLX,TOLF
      PARAMETER (TOLX=EPS12,TOLF=EPS12)
C
C Dummy arguments
C
      REAL*8 ALPHA,B0,CLIGHT,EIN,EREL,MU,V00,ZZ
      INTEGER ATOM,FLAG,KK,LAY
      REAL*8 BR(RSTEP),EMAG(NSOL),RAD(RSTEP),RADP4(RSTEPP5),VR(RSTEP)
      COMPLEX*16 CWFFS(RSTEP,2,2),CWFGS(RSTEP,2,2),
     &           VLM(LAYSM,NATLM,MQD,RSTEP)
C
C Local variables
C
      REAL*8 CHIFF(2),CHIFG,CHIGG(2),D,DX,DYDX(N,N),ERRF,ERRX,FY(N),
     &       INDX(N),KAPPA(2),S,SOLF(:,:,:),SOLG(:,:,:),X(N)
      LOGICAL COUPLE,OFFDIAG
      INTEGER I,ISOL,J,K,KB,RM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SOLF,SOLG
      ALLOCATE (SOLF(RSTEP,2,2),SOLG(RSTEP,2,2))
C
C     clear output fields
      DO I = 1,RSTEP
         DO J = 1,2
            DO K = 1,2
               CWFFS(I,J,K) = CZERO
               CWFGS(I,J,K) = CZERO
            END DO
         END DO
      END DO
C
C     energy parameter
      X(1) = EIN
      DX = MIN(0.01D0,EIN*1.D-4)
C
C     check how many solutions have to be calculated
Ccghf may go partially one subroutine higher
C      ISOL = 1
C     no b-field -> no coupling !!
      COUPLE = .FALSE.
      OFFDIAG = .FALSE.
      IF ( ABS(B0).GT.EPS12 ) THEN
         COUPLE = .TRUE.
         OFFDIAG = .TRUE.
C         no coupling for j=l+s (k<0) and |mj|=j !!
         IF ( (KK.LT.0 .AND. ABS(MU)+0.5D0.GE.ABS(KK)) )
     &        OFFDIAG = .FALSE.
         KB = -KK - 1
         IF ( (KB.LT.0 .AND. ABS(MU)+0.5D0.GE.ABS(KB)) )
     &        OFFDIAG = .FALSE.
      END IF
Ccghf the correct selection of the no of solution has to be implemented
Ccghf later together with the bxy case                !!!!!!!!!!!!!!!!!
C     if(offdiag) isol = 2
      ISOL = 2
C
      KAPPA(1) = DBLE(KK)
      KAPPA(2) = -DBLE(KK) - 1.D0
      S = 0.5D0
      CHIGG(1) = 0.D0
      CHIGG(2) = 0.D0
      CHIFF(1) = 0.D0
      CHIFF(2) = 0.D0
      CHIFG = 0.D0
C
C     ***** chi : matrix elements of sigma_z induce a coupling *****
C     *****       of radial wave  functions with the same l    *****
C     *****       for spinpolarized dirac equation             *****
C     ***** magnetization along z                              *****
C     ***** diagonal    :  chi_ii        <k||sz||k>            *****
C     ***** off-diagonal:  chi_fg        <k||sz||-k-1>         *****
C
      IF ( COUPLE ) THEN
         DO I = 1,ISOL
            CHIGG(I) = -MU/(-KAPPA(I)+S)
            CHIFF(I) = -MU/(KAPPA(I)+S)
         END DO
C
C             note: result should not change
C                   if kappa(1) is replaced by kappa(2)
         IF ( OFFDIAG ) CHIFG = -SQRT(1.D0-(MU/(KAPPA(1)+S))**2)
      END IF
C
      DO K = 1,NTRIAL
C         initialize in first loop (set flag=1)
         IF ( K.EQ.1 ) FLAG = 1
C         usrfun returns flag=0 (ok) or 3 (fault)
         CALL USRFUN(ISOL,FLAG,KAPPA,CHIGG,CHIFF,CHIFG,ZZ,CLIGHT,V00,B0,
     &               VR,BR,RAD,RADP4,LAY,ATOM,VLM,ALPHA,RM,X,DYDX,FY,DX,
     &               SOLF,SOLG)
         IF ( FLAG.EQ.3 ) THEN
C             something was wrong in usrfun
            EREL = 0.D0
            EXIT
         END IF
C
C         start newton-raphson (flag=0 see usrfun)
         ERRF = 0.D0
         DO I = 1,N
            ERRF = ERRF + ABS(FY(I))
         END DO
         IF ( ERRF.LE.TOLF ) THEN
C             binding energy was found
            FLAG = 2
            EXIT
         END IF
C
         CALL LUDCMP(DYDX,N,NP,INDX,D)
         CALL LUBKSB(DYDX,N,NP,INDX,FY)
         ERRX = 0.D0
         DO I = 1,N
            ERRX = ERRX + ABS(FY(I))
            X(I) = X(I) + FY(I)
         END DO
         IF ( ERRX.LE.TOLX ) THEN
C             binding energy was found
            FLAG = 2
            EXIT
         END IF
C         repeat untill ntrial is exceeded
      END DO
C
      IF ( FLAG.EQ.2 ) THEN
C         binding energy was found, calculate final wave-functions
         EREL = X(1)
         CALL MAKE_CWF(X,SOLF,SOLG,COUPLE,OFFDIAG,1,CHIGG,CHIFF,CHIFG,
     &                 BR,RAD,RM,ALPHA,EMAG,CWFFS,CWFGS)
      ELSE
C         no energy found, not converged, or an fault occured
         DX = 0.1D0
         CALL USRFUN(ISOL,FLAG,KAPPA,CHIGG,CHIFF,CHIFG,ZZ,CLIGHT,V00,B0,
     &               VR,BR,RAD,RADP4,LAY,ATOM,VLM,ALPHA,RM,X,DYDX,FY,DX,
     &               SOLF,SOLG)
C
         WRITE (NOUT1,99001)
         WRITE (NOUT1,99002) KK,MU,(FY(I),I=1,4)
         WRITE (*,99001)
         WRITE (*,99002) KK,MU,(FY(I),I=1,4)
         EREL = 0.D0
      END IF
C
      RETURN
C
99001 FORMAT ('no core-energy found in mnewt for')
99002 FORMAT ('state: k,m:',i3,f4.1,1x,'fy:',4(1x,e11.4))
      END
C*==usrfun.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE USRFUN(ISOL,FLAG,KAPPA,CHIGG,CHIFF,CHIFG,ZZ,CLIGHT,V00,
     &                  B0,VR,BR,RAD,RADP4,LAY,ATOM,VLM,ALPHA,RM,X,DYDX,
     &                  FY,DX,SOLF,SOLG)
C     /****************************************************************/
C     # purpose      : defines how to integrate the dirac equation     *
C                                                                      *
C     # note:        : always call first inward before outward         *
C                      avoids overwriting of the correct value at rm   *
C                                                                      *
C     # calls the subroutines and function:                            *
C       ctp   inward      outward                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,RSTEPP5
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,B0,CHIFG,CLIGHT,DX,V00,ZZ
      INTEGER ATOM,FLAG,ISOL,LAY,RM
      REAL*8 BR(RSTEP),CHIFF(2),CHIGG(2),DYDX(4,4),FY(4),KAPPA(2),
     &       RAD(RSTEP),RADP4(RSTEPP5),SOLF(RSTEP,2,2),SOLG(RSTEP,2,2),
     &       VR(RSTEP),X(4)
      COMPLEX*16 VLM(LAYSM,NATLM,MQD,RSTEP)
C
C Local variables
C
      REAL*8 A,CIN,COUT,DE,DFY(4),E,FIN(2,2),FOUT(2,2),GIN(2,2),
     &       GOUT(2,2)
      INTEGER CTP
      INTEGER I
      EXTERNAL CTP
C
C*** End of declarations rewritten by SPAG
C
      E = X(1)
      IF ( FLAG.EQ.1 ) RM = CTP(VLM,E,LAY,ATOM)
C      IF ( RM.EQ.0 .OR. X(1).EQ.0.D0 ) THEN
      IF ( RM.EQ.0 .OR. ABS(X(1)).LE.1.0D-16 ) THEN
         FLAG = 3
         RETURN
      END IF
C
C     calculate solution at energy e
      CALL INWARD(E,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RADP4,ALPHA,
     &            VR,BR,SOLF,SOLG,FIN,GIN)
      CALL OUTWARD(E,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RAD,ALPHA,
     &             ZZ,V00,B0,VR,BR,SOLF,SOLG,FOUT,GOUT)
C
      IF ( FLAG.EQ.1 ) THEN
C         init parameters for first iteration of newton raphson:
C         ======================================================
         FLAG = 0
         X(2) = FOUT(1,1)/FIN(1,1)
         X(3) = 0.0D0
         X(4) = 0.0D0
      END IF
C
C     used as test if not converged
C      IF ( DX.EQ.0.1D0 ) X(2) = FOUT(1,1)/FIN(1,1)
      IF ( ABS(DX-0.1D0).LE.1.0D-16 ) X(2) = FOUT(1,1)/FIN(1,1)
C
      E = X(1)
      A = X(2)
      COUT = X(3)
      CIN = X(4)
C
      FY(1) = -((FOUT(1,1)+COUT*FOUT(1,2))-A*(FIN(1,1)+CIN*FIN(1,2)))
      FY(2) = -((GOUT(1,1)+COUT*GOUT(1,2))-A*(GIN(1,1)+CIN*GIN(1,2)))
      FY(3) = -((FOUT(2,1)+COUT*FOUT(2,2))-A*(FIN(2,1)+CIN*FIN(2,2)))
      FY(4) = -((GOUT(2,1)+COUT*GOUT(2,2))-A*(GIN(2,1)+CIN*GIN(2,2)))
C
C     calculate deltaf_alph/deltax_beta:
C     ===================================
      DYDX(1,2) = -(FIN(1,1)+CIN*FIN(1,2))
      DYDX(2,2) = -(GIN(1,1)+CIN*GIN(1,2))
      DYDX(3,2) = -(FIN(2,1)+CIN*FIN(2,2))
      DYDX(4,2) = -(GIN(2,1)+CIN*GIN(2,2))
C
      DYDX(1,3) = FOUT(1,2)
      DYDX(2,3) = GOUT(1,2)
      DYDX(3,3) = FOUT(2,2)
      DYDX(4,3) = GOUT(2,2)
C
      DYDX(1,4) = -A*FIN(1,2)
      DYDX(2,4) = -A*GIN(1,2)
      DYDX(3,4) = -A*FIN(2,2)
      DYDX(4,4) = -A*GIN(2,2)
C
C     calculate solutions at energy e + dx for derivatives
      DE = E + DX
      CALL INWARD(DE,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RADP4,ALPHA,
     &            VR,BR,SOLF,SOLG,FIN,GIN)
      CALL OUTWARD(DE,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RAD,ALPHA,
     &             ZZ,V00,B0,VR,BR,SOLF,SOLG,FOUT,GOUT)
C
      DFY(1) = -((FOUT(1,1)+COUT*FOUT(1,2))-A*(FIN(1,1)+CIN*FIN(1,2)))
      DFY(2) = -((GOUT(1,1)+COUT*GOUT(1,2))-A*(GIN(1,1)+CIN*GIN(1,2)))
      DFY(3) = -((FOUT(2,1)+COUT*FOUT(2,2))-A*(FIN(2,1)+CIN*FIN(2,2)))
      DFY(4) = -((GOUT(2,1)+COUT*GOUT(2,2))-A*(GIN(2,1)+CIN*GIN(2,2)))
C
      DO I = 1,4
         DYDX(I,1) = -(DFY(I)-FY(I))/DX
      END DO
C
      END
C*==make_cwf.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE MAKE_CWF(X,SOLF,SOLG,COUPLE,OFFDIAG,JSOL,CHIGG,CHIFF,
     &                    CHIFG,BR,RAD,RM,ALPHA,EMAG,CWFFS,CWFGS)
C     /****************************************************************/
C     # purpose      : normalization of core wavefunctions             *
C                      calculate off-diagonal energy for coupled       *
C                      wave functions                                  *
C                                                                      *
C     # calls the following subroutine:                                *
C       simps                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:RSTEP,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CHIFG
      LOGICAL COUPLE,OFFDIAG
      INTEGER JSOL,RM
      REAL*8 BR(RSTEP),CHIFF(2),CHIGG(2),EMAG(2),RAD(RSTEP),
     &       SOLF(RSTEP,2,2),SOLG(RSTEP,2,2),X(4)
      COMPLEX*16 CWFFS(RSTEP,2,2),CWFGS(RSTEP,2,2)
C
C Local variables
C
      REAL*8 A,CIN,COUT
      COMPLEX*16 F1(:),F2(:),G1(:),G2(:),INTSQA,SOL(:,:),SQA(:),SQNORM
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F1,F2,G1,G2,SQA,SOL
      ALLOCATE (F1(RSTEP),F2(RSTEP),G1(RSTEP),G2(RSTEP),SQA(RSTEP))
      ALLOCATE (SOL(RSTEP,4))
C
      A = X(2)
      COUT = X(3)
      CIN = X(4)
C
C          write(*,'(i4, 4(1x, e12.5))') rm, e*hartre, cout, a, cin
C
      DO I = 1,RSTEP
         DO J = 1,4
            SOL(I,J) = CZERO
         END DO
      END DO
C     outward [0,rm]
      DO I = 1,RM
         SOL(I,1) = CONE*(SOLF(I,1,1)+COUT*SOLF(I,1,2))
         SOL(I,2) = CONE*(SOLG(I,1,1)+COUT*SOLG(I,1,2))
      END DO
C     matched inward [rm+1,infinity=rstep]
      DO I = RM + 1,RSTEP
         SOL(I,1) = CONE*A*(SOLF(I,1,1)+CIN*SOLF(I,1,2))
         SOL(I,2) = CONE*A*(SOLG(I,1,1)+CIN*SOLG(I,1,2))
      END DO
C     same if offdiagonal waves are present
      IF ( OFFDIAG ) THEN
         DO I = 1,RM
            SOL(I,3) = CONE*(SOLF(I,2,1)+COUT*SOLF(I,2,2))
            SOL(I,4) = CONE*(SOLG(I,2,1)+COUT*SOLG(I,2,2))
         END DO
         DO I = RM + 1,RSTEP
            SOL(I,3) = CONE*A*(SOLF(I,2,1)+CIN*SOLF(I,2,2))
            SOL(I,4) = CONE*A*(SOLG(I,2,1)+CIN*SOLG(I,2,2))
         END DO
      END IF
C
C     determine norm of wavefunction
      DO I = 1,RSTEP
         F1(I) = SOL(I,1)
         G1(I) = SOL(I,2)
         F2(I) = SOL(I,3)
         G2(I) = SOL(I,4)
         SQA(I) = (F1(I)*DCONJG(F1(I))+G1(I)*DCONJG(G1(I))+F2(I)
     &            *DCONJG(F2(I))+G2(I)*DCONJG(G2(I)))*RAD(I)
      END DO
C
      CALL SIMPS(ALPHA,SQA,INTSQA)
      SQNORM = 1.D0/CDSQRT(INTSQA)
C
C     normalize wavefunction and write into cwf array
      DO I = 1,RSTEP
         CWFFS(I,1,1) = F1(I)*SQNORM
         CWFFS(I,2,2) = F2(I)*SQNORM
         CWFGS(I,1,1) = G1(I)*SQNORM
         CWFGS(I,2,2) = G2(I)*SQNORM
      END DO
C
C     calculate magnetic energy terms in case of magnetic field
      EMAG(1) = 0.D0
      EMAG(2) = 0.D0
C
C     diagonal magnetic energy terms
      IF ( COUPLE ) THEN
         DO I = 1,RSTEP
            SQA(I) = (DCONJG(F1(I))*F1(I)*CHIFF(JSOL)-DCONJG(G1(I))
     &               *G1(I)*CHIGG(JSOL))*BR(I)*SQNORM
         END DO
         CALL SIMPS(ALPHA,SQA,INTSQA)
         EMAG(1) = DBLE(INTSQA)
      END IF
C
C     off-diagonal energy terms
      IF ( OFFDIAG ) THEN
         DO I = 1,RSTEP
            SQA(I) = DCONJG(F1(I))*F2(I)*CHIFG*BR(I)*SQNORM
         END DO
         CALL SIMPS(ALPHA,SQA,INTSQA)
         EMAG(2) = DBLE(INTSQA)
      END IF
C
      END
C*==inward.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE INWARD(E,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RADP4,
     &                  ALPHA,VR,BR,SOLF,SOLG,FIN,GIN)
C     /****************************************************************/
C     # purpose      : integrate dirac equation from the outside to    *
C                      a point rm in the extended mesh radp4           *
C                                                                      *
C     # parameters   :                                                 *
C     # input  :                                                       *
C        e     : energy for which to solve dirac equation              *
C        isol  : no of coupled equations (1 or 2)                      *
C        kappa : relativistic quantum number                           *
C        chiij : coupling oefficients                                  *
C        rm    : matchig point                                         *
C        clight: speed of light                                        *
C        radp4 : mesh on which to solve the equation                   *
C        alpha : parameter for exponential mesh                        *
C        vr    : potential times r                                     *
C        br    : magnetic field times r                                *
C     # output     :                                                   *
C        solf,solg : inward solution of dirac equation in [rm,rstep]   *
C        fin, gin  : inward solution of dirac equation at rm           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:RSTEP,RSTEPP5,EPS12
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,CHIFG,CLIGHT,E
      INTEGER ISOL,RM
      REAL*8 BR(RSTEP),CHIFF(2),CHIGG(2),FIN(2,2),GIN(2,2),KAPPA(2),
     &       RADP4(RSTEPP5),SOLF(RSTEP,2,2),SOLG(RSTEP,2,2),VR(RSTEP)
C
C Local variables
C
      REAL*8 DF(2,2,4),DG(2,2,4),EF(2),EFF,EG(2),EGG,F(:,:,:),G(:,:,:),
     &       GAMMA,H,H3,KF(2),KG(2),OF(2,2),OG(2,2),PK,QK,R,RCP,V1,V2
      INTEGER I,ITER,J,K,M,NITER
      LOGICAL ICONVERGE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F,G
      ALLOCATE (F(2,2,RSTEPP5),G(2,2,RSTEPP5))
C
C*** End of declarations rewritten by SPAG
C
C     ***** maximum number of iterations in pece algorithm *****
      NITER = 30
C
      H = ALPHA
      H3 = -H/3.D0
C
      DO I = 1,ISOL
         KF(I) = -KAPPA(I)
         KG(I) = KAPPA(I)
      END DO
C
      GAMMA = SQRT(ABS(-E-(E/CLIGHT)**2))
      PK = 1.D0
      QK = -PK*GAMMA/(1.D0+E/(CLIGHT**2))
C
C     ***** inward wkb solution for large r (r > rstep -> vr,br=0) *****
      DO K = RSTEP + 4,RSTEP + 1, - 1
         M = RSTEP + 4 - K + 1
         R = RADP4(K)
         DO J = 1,ISOL
            F(J,J,K) = PK*EXP(-GAMMA*RADP4(K))*R
            G(J,J,K) = QK*EXP(-GAMMA*RADP4(K))*R
            DF(J,J,M) = H3*(-GAMMA*F(J,J,K)+F(J,J,K)/R)
            DG(J,J,M) = H3*(-GAMMA*G(J,J,K)+G(J,J,K)/R)
C
            I = 3 - J
            F(I,J,K) = 0.0D0
            G(I,J,K) = 0.0D0
            DF(I,J,M) = 0.0D0
            DG(I,J,M) = 0.0D0
         END DO
      END DO
C
C     ***** inward integration to muffin-tin-radius *****
      DO K = RSTEP,RM, - 1
C         ***** predictor *****
         DO J = 1,ISOL
            DO I = 1,ISOL
               F(I,J,K) = F(I,J,K+4)
     &                    + 8.D0*(DF(I,J,3)-0.5D0*DF(I,J,2)+DF(I,J,1))
               G(I,J,K) = G(I,J,K+4)
     &                    + 8.D0*(DG(I,J,3)-0.5D0*DG(I,J,2)+DG(I,J,1))
            END DO
         END DO
C         ***** corrector *****
         R = RADP4(K)
         EGG = (E*R-VR(K))/CLIGHT + 2.D0*CLIGHT*R
         EFF = -(E*R-VR(K))/CLIGHT
         RCP = BR(K)*CHIFG/CLIGHT
C
         DO I = 1,ISOL
            EG(I) = EGG - BR(K)*CHIGG(I)/CLIGHT
            EF(I) = EFF + BR(K)*CHIFF(I)/CLIGHT
         END DO
C
         DO ITER = 1,NITER
            ICONVERGE = .FALSE.
            DO J = 1,ISOL
               DO I = 1,ISOL
                  OF(I,J) = F(I,J,K)
                  OG(I,J) = G(I,J,K)
C                     ***** derivate at the new point  *****
                  DF(I,J,4) = H3*(KF(I)*F(I,J,K)+EG(I)*G(I,J,K))
                  DG(I,J,4) = H3*(KG(I)*G(I,J,K)+EF(I)*F(I,J,K)+RCP*F(3-
     &                        I,J,K))
C                     ***** functions at the new point *****
                  F(I,J,K) = F(I,J,K+1) + 1.125*DF(I,J,4)
     &                       + 2.375*DF(I,J,3) - 0.625*DF(I,J,2)
     &                       + 0.125*DF(I,J,1)
                  G(I,J,K) = G(I,J,K+1) + 1.125*DG(I,J,4)
     &                       + 2.375*DG(I,J,3) - 0.625*DG(I,J,2)
     &                       + 0.125*DG(I,J,1)
C                     ***** test if converged *****
                  V1 = ABS(F(I,J,K)-OF(I,J)) - ABS(F(I,J,K)*EPS12)
                  V2 = ABS(G(I,J,K)-OG(I,J)) - ABS(G(I,J,K)*EPS12)
                  IF ( V1.LT.0. .OR. V2.LT.0. ) ICONVERGE = .TRUE.
               END DO
            END DO
            IF ( ICONVERGE ) EXIT
         END DO
C          if (iter.ge.niter) write(*,*) 'poor convergence in inward'
C
C         ***** restore derivatives for calculation of the next point *****
         DO J = 1,ISOL
            DO I = 1,ISOL
               DF(I,J,1) = DF(I,J,2)
               DG(I,J,1) = DG(I,J,2)
               DF(I,J,2) = DF(I,J,3)
               DG(I,J,2) = DG(I,J,3)
C                 ***** final derivative at the new point *****
               DF(I,J,3) = H3*(KF(I)*F(I,J,K)+EG(I)*G(I,J,K))
               DG(I,J,3) = H3*(KG(I)*G(I,J,K)+EF(I)*F(I,J,K)+RCP*F(3-I,J
     &                     ,K))
            END DO
         END DO
      END DO
C
C     ***** save inward solution *****
      DO K = RM,RSTEP
         R = 1.D0/RADP4(K)
         DO J = 1,ISOL
            DO I = 1,ISOL
               SOLF(K,I,J) = F(I,J,K)*R
               SOLG(K,I,J) = G(I,J,K)*R
            END DO
         END DO
      END DO
C     ***** inward solution at rm *****
      DO J = 1,ISOL
         DO I = 1,ISOL
            FIN(I,J) = F(I,J,RM)
            GIN(I,J) = CLIGHT*G(I,J,RM)
         END DO
      END DO
C
      END
C*==outward.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE OUTWARD(E,ISOL,KAPPA,CHIGG,CHIFF,CHIFG,RM,CLIGHT,RAD,
     &                   ALPHA,ZZ,V00,B0,VR,BR,SOLF,SOLG,FOUT,GOUT)
C     /****************************************************************/
C     # purpose      : integrate dirac equation from the origin to     *
C                      a point rm in the mesh rad                      *
C     # parameters   :                                                 *
C     # input  :                                                       *
C        e     : energy for which to solve dirac equation              *
C        isol  : no of coupled equations (1 or 2)                      *
C        kappa : relativistic quantum number                           *
C        chiij : coupling oefficients                                  *
C        rm    : matchig point                                         *
C        clight: speed of light                                        *
C        rad   : mesh on which to solve the equation                   *
C        alpha : parameter for exponential mesh                        *
C        zz    : atomic number                                         *
C        v00,b0: potential and magnetic field  at r=0                  *
C        vr    : potential times r                                     *
C        br    : magnetic field times r                                *
C     # output      :                                                  *
C        solf, solg : outward solution of dirac equation in [1,rm]     *
C        fout, gout : outward solution of dirac equation at rm         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:RSTEP,EPS12
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,B0,CHIFG,CLIGHT,E,V00,ZZ
      INTEGER ISOL,RM
      REAL*8 BR(RSTEP),CHIFF(2),CHIGG(2),FOUT(2,2),GOUT(2,2),KAPPA(2),
     &       RAD(RSTEP),SOLF(RSTEP,2,2),SOLG(RSTEP,2,2),VR(RSTEP)
C
C Local variables
C
      REAL*8 CP,D11,D22,DELTA,DF(2,2,4),DG(2,2,4),EF(2),EF0(2),EFF,EFF0,
     &       EG(2),EG0(2),EGG,EGG0,F(:,:,:),F0(2,2,50),FIP,G(:,:,:),
     &       G0(2,2,50),GAM(2),H,H3,KF(2),KG(2),OF(2,2),OG(2,2),Q12,Q21,
     &       R,RCP,RG,V1,V2,WFF,WGG
      INTEGER I,IP,ITER,J,K,NITER,NPOW
      LOGICAL ICONVERGE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE F,G
      ALLOCATE (F(2,2,RSTEP),G(2,2,RSTEP))
C
C     ***** maximum number of coefficients           *****
C     ***** for power series expansion close to zero *****
      NPOW = 25
C
C     maximum number of iterations in in pece algorithm *****
      NITER = 30
C
      H = ALPHA
      H3 = H/3.D0
C
      Q12 = -ZZ/CLIGHT
      Q21 = ZZ/CLIGHT
C
      DO I = 1,ISOL
         KF(I) = -KAPPA(I)
         KG(I) = KAPPA(I)
         GAM(I) = SQRT(ABS(KAPPA(I)**2+Q12*Q21))
      END DO
C
C     ***** power series expansion of the radial wave functions *****
C     ***** construct coefficients for the expansion            *****
C
      EGG0 = (E-V00)/CLIGHT + 2.D0*CLIGHT
      EFF0 = -(E-V00)/CLIGHT
      CP = B0*CHIFG/CLIGHT
C
      DO I = 1,ISOL
         EG0(I) = EGG0 - B0*CHIGG(I)/CLIGHT
         EF0(I) = EFF0 + B0*CHIFF(I)/CLIGHT
      END DO
C
C     **** boundary conditions ****
      DO J = 1,ISOL
         I = 3 - J
         F0(J,J,1) = SQRT(ABS(KAPPA(J)-GAM(J)))
         G0(J,J,1) = SIGN(1.D0,KAPPA(J))*SQRT(ABS(KAPPA(J)+GAM(J)))
C
         F0(I,J,1) = 0.0D0
         G0(I,J,1) = 0.0D0
      END DO
C
C     **** coefficients ****
      FIP = 0.
      DO IP = 2,NPOW
         FIP = FIP + 1.D0
         DO J = 1,ISOL
            DO I = 1,ISOL
               WFF = EG0(I)*G0(I,J,IP-1)
               WGG = EF0(I)*F0(I,J,IP-1) + CP*F0(3-I,J,IP-1)
C
               D11 = GAM(J) + KAPPA(I) + FIP
               D22 = GAM(J) - KAPPA(I) + FIP
               DELTA = D11*D22 - Q21*Q12
C
               F0(I,J,IP) = (WFF*D22-Q21*WGG)/DELTA
               G0(I,J,IP) = (D11*WGG-WFF*Q12)/DELTA
            END DO
         END DO
      END DO
C
C     ***** determine values of the radial wave functions for *****
C     ***** the first four points of the radial mesh and  the *****
C     ***** derivatives for the first 3 integration points    *****
C
      DO I = 1,2
         DO J = 1,2
            DO K = 1,RM
               F(I,J,K) = 0.D0
               G(I,J,K) = 0.D0
            END DO
         END DO
      END DO
C
      DO I = 1,2
         DO J = 1,2
            DO K = 1,4
               DF(I,J,K) = 0.D0
               DG(I,J,K) = 0.D0
            END DO
         END DO
      END DO
C
      DO K = 1,4
         R = RAD(K)
         DO J = 1,ISOL
            RG = R**(GAM(J))
            DO I = 1,ISOL
               F(I,J,K) = F0(I,J,1)*RG
               G(I,J,K) = G0(I,J,1)*RG
               DF(I,J,K) = H3*(GAM(J)+0)*F0(I,J,1)*RG
               DG(I,J,K) = H3*(GAM(J)+0)*G0(I,J,1)*RG
            END DO
C
            DO IP = 2,NPOW
               RG = RG*R
               DO I = 1,ISOL
                  F(I,J,K) = F(I,J,K) + F0(I,J,IP)*RG
                  G(I,J,K) = G(I,J,K) + G0(I,J,IP)*RG
                  DF(I,J,K) = DF(I,J,K) + H3*(GAM(J)+IP-1)*F0(I,J,IP)*RG
                  DG(I,J,K) = DG(I,J,K) + H3*(GAM(J)+IP-1)*G0(I,J,IP)*RG
               END DO
            END DO
         END DO
C
C         ***** shift index of derivatives back by 1 *****
         DO J = 1,ISOL
            DO I = 1,ISOL
               DF(I,J,1) = DF(I,J,2)
               DG(I,J,1) = DG(I,J,2)
               DF(I,J,2) = DF(I,J,3)
               DG(I,J,2) = DG(I,J,3)
               DF(I,J,3) = DF(I,J,4)
               DG(I,J,3) = DG(I,J,4)
            END DO
         END DO
      END DO
C
C     ***** outward integration to muffin-tin-radius *****
      DO K = 5,RM
C         ***** predictor *****
         DO J = 1,ISOL
            DO I = 1,ISOL
               F(I,J,K) = F(I,J,K-4)
     &                    + 8.D0*(DF(I,J,3)-0.5D0*DF(I,J,2)+DF(I,J,1))
               G(I,J,K) = G(I,J,K-4)
     &                    + 8.D0*(DG(I,J,3)-0.5D0*DG(I,J,2)+DG(I,J,1))
            END DO
         END DO
C         ***** corrector *****
         R = RAD(K)
         EGG = (E*R-VR(K))/CLIGHT + 2.D0*CLIGHT*R
         EFF = -(E*R-VR(K))/CLIGHT
         RCP = BR(K)*CHIFG/CLIGHT
C
         DO I = 1,ISOL
            EG(I) = EGG - BR(K)*CHIGG(I)/CLIGHT
            EF(I) = EFF + BR(K)*CHIFF(I)/CLIGHT
         END DO
C
         DO ITER = 1,NITER
            ICONVERGE = .FALSE.
            DO J = 1,ISOL
               DO I = 1,ISOL
                  OF(I,J) = F(I,J,K)
                  OG(I,J) = G(I,J,K)
C                     ***** derivate at the new point *****
                  DF(I,J,4) = H3*(KF(I)*F(I,J,K)+EG(I)*G(I,J,K))
                  DG(I,J,4) = H3*(KG(I)*G(I,J,K)+EF(I)*F(I,J,K)+RCP*F(3-
     &                        I,J,K))
C                     ***** functions at the new point *****
                  F(I,J,K) = F(I,J,K-1) + 1.125*DF(I,J,4)
     &                       + 2.375*DF(I,J,3) - 0.625*DF(I,J,2)
     &                       + 0.125*DF(I,J,1)
                  G(I,J,K) = G(I,J,K-1) + 1.125*DG(I,J,4)
     &                       + 2.375*DG(I,J,3) - 0.625*DG(I,J,2)
     &                       + 0.125*DG(I,J,1)
C                     ***** test if converged *****
                  V1 = ABS(F(I,J,K)-OF(I,J)) - ABS(F(I,J,K)*EPS12)
                  V2 = ABS(G(I,J,K)-OG(I,J)) - ABS(G(I,J,K)*EPS12)
                  IF ( V1.LT.0. .OR. V2.LT.0. ) ICONVERGE = .TRUE.
               END DO
            END DO
            IF ( ICONVERGE ) EXIT
         END DO
C          if (iter.ge.niter) write(*,*) 'poor convergence in inward'
C
C         ***** restore derivatives for calculation of the next point *****
         DO J = 1,ISOL
            DO I = 1,ISOL
               DF(I,J,1) = DF(I,J,2)
               DG(I,J,1) = DG(I,J,2)
               DF(I,J,2) = DF(I,J,3)
               DG(I,J,2) = DG(I,J,3)
C                 ***** final derivative at the new point *****
               DF(I,J,3) = H3*(KF(I)*F(I,J,K)+EG(I)*G(I,J,K))
               DG(I,J,3) = H3*(KG(I)*G(I,J,K)+EF(I)*F(I,J,K)+RCP*F(3-I,J
     &                     ,K))
            END DO
         END DO
      END DO
C
C     ***** save outward solution *****
      DO K = 1,RM
         R = 1.D0/RAD(K)
         DO J = 1,ISOL
            DO I = 1,ISOL
               SOLF(K,I,J) = F(I,J,K)*R
               SOLG(K,I,J) = G(I,J,K)*R
            END DO
         END DO
      END DO
C     ***** outward solution at rm *****
      DO J = 1,ISOL
         DO I = 1,ISOL
            FOUT(I,J) = F(I,J,RM)
            GOUT(I,J) = CLIGHT*G(I,J,RM)
         END DO
      END DO
C
      END
C*==ludcmp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C     /****************************************************************/
C     from numerical recipes                                           *
C     /****************************************************************/
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NMAX
      PARAMETER (NMAX=100)
      REAL*8 TINYVAR
      PARAMETER (TINYVAR=1.0D-20)
C
C Dummy arguments
C
      REAL*8 D
      INTEGER N,NP
      REAL*8 A(NP,NP),INDX(N)
C
C Local variables
C
      REAL*8 AAMAX,DUM,SUMVAR,VV(NMAX)
      INTEGER I,IMAX,J,K
C
C*** End of declarations rewritten by SPAG
C
C
      D = 1.D0
      DO I = 1,N
         AAMAX = 0.D0
         DO J = 1,N
            IF ( ABS(A(I,J)).GT.AAMAX ) AAMAX = ABS(A(I,J))
         END DO
C         IF ( AAMAX.EQ.0. ) WRITE (*,*) 'singular matrix.'
         IF ( ABS(AAMAX).LE.1.0D-16 ) WRITE (*,*) 'singular matrix.'
         VV(I) = 1./AAMAX
      END DO
C
      DO J = 1,N
         IF ( J.GT.1 ) THEN
            DO I = 1,J - 1
               SUMVAR = A(I,J)
               IF ( I.GT.1 ) THEN
                  DO K = 1,I - 1
                     SUMVAR = SUMVAR - A(I,K)*A(K,J)
                  END DO
                  A(I,J) = SUMVAR
               END IF
            END DO
         END IF
C
         AAMAX = 0.D0
         DO I = J,N
            SUMVAR = A(I,J)
            IF ( J.GT.1 ) THEN
               DO K = 1,J - 1
                  SUMVAR = SUMVAR - A(I,K)*A(K,J)
               END DO
               A(I,J) = SUMVAR
            END IF
            DUM = VV(I)*ABS(SUMVAR)
            IF ( DUM.GE.AAMAX ) THEN
               IMAX = I
               AAMAX = DUM
            END IF
         END DO
         IF ( J.NE.IMAX ) THEN
            DO K = 1,N
               DUM = A(IMAX,K)
               A(IMAX,K) = A(J,K)
               A(J,K) = DUM
            END DO
            D = -D
            VV(IMAX) = VV(J)
         END IF
         INDX(J) = IMAX
         IF ( J.NE.N ) THEN
C            IF ( A(J,J).EQ.0.D0 ) A(J,J) = TINYVAR
            IF ( ABS(A(J,J)).LE.1.0D-16 ) A(J,J) = TINYVAR
            DUM = 1./A(J,J)
            DO I = J + 1,N
               A(I,J) = A(I,J)*DUM
            END DO
         END IF
      END DO
C
C      IF ( A(N,N).EQ.0. ) A(N,N) = TINYVAR
      IF ( ABS(A(N,N)).LE.1.0D-16 ) A(N,N) = TINYVAR
C
      END
C*==lubksb.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C     /****************************************************************/
C     from numerical recipes                                           *
C     /****************************************************************/
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,NP
      REAL*8 A(NP,NP),B(N),INDX(N)
C
C Local variables
C
      INTEGER I,II,J
      REAL*8 LL,SUMVAR
C
C*** End of declarations rewritten by SPAG
C
      II = 0
      DO I = 1,N
         LL = INDX(I)
         SUMVAR = B(INT(LL))
         B(INT(LL)) = B(I)
         IF ( II.NE.0 ) THEN
            DO J = II,I - 1
               SUMVAR = SUMVAR - A(I,J)*B(J)
            END DO
C         ELSE IF ( SUMVAR.NE.0. ) THEN
         ELSE IF ( ABS(SUMVAR).GT.1.0D-16 ) THEN
            II = I
         END IF
         B(I) = SUMVAR
      END DO
C
      DO I = N,1, - 1
         SUMVAR = B(I)
         IF ( I.LT.N ) THEN
            DO J = I + 1,N
               SUMVAR = SUMVAR - A(I,J)*B(J)
            END DO
         END IF
         B(I) = SUMVAR/A(I,I)
      END DO
C
      END
C*==writecore.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE WRITECORE(EVAL,EMAG,KAPPA,MUE,N,ATOM,LAY,CHANNEL)
C     /****************************************************************/
C     # purpose      : write core wavefunctions and energies to file   *
C                                                                      *
C     # calls the following functions:                                 *
C       simps                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:RSTEP,HARTRE
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,CHANNEL,KAPPA,LAY,N
      REAL*8 EVAL,MUE
      REAL*8 EMAG(2)
C
C Local variables
C
      INTEGER I,KL,KNOTS,NQN
      REAL*8 KJ
C
C*** End of declarations rewritten by SPAG
C
      KNOTS = 0
C      IF ( ABS(EVAL).NE.0.D0 ) THEN
      IF ( ABS(EVAL).GT.1.0D-16 ) THEN
         WRITE (CHANNEL,99001) KAPPA,MUE,N,ATOM,LAY,EVAL*HARTRE
         DO I = 1,RSTEP
C             write (channel,1000) rad(i),
C    1                            ((cwffs(i,j,k),j=1,2),k=1,2)
C             if (i.gt.1
C    1        .and.dble(cwffs(i,1,1)*cwffs(i-1,1,1)).lt.0.d0)
C    2        knots=knots+1
         END DO
         WRITE (CHANNEL,99002) KAPPA,MUE,N,ATOM,LAY,EVAL*HARTRE
         DO I = 1,RSTEP
C             write (channel,1000) rad(i),
C    1                             ((cwfgs(i,j,k),j=1,2),k=1,2)
         END DO
      END IF
C
C     write file with binding energies
C     (maybe moved to outputs making use of readcore)
C
      IF ( KAPPA.GT.0 ) THEN
         KL = KAPPA
         KJ = KAPPA - 0.5
      ELSE
         KL = IABS(KAPPA+1)
         KJ = ABS(KAPPA+0.5)
      END IF
      NQN = 1 + KNOTS + KL
      WRITE (NOUT1,99003) NQN,KAPPA,KL,KJ,MUE,EVAL*HARTRE,EMAG(1)
     &                    *HARTRE*1.D3,EMAG(2)*HARTRE*1.D3
Ccghfend
      RETURN
C
99001 FORMAT ('#f: kappa=',i3,1x,'mu=',f5.1,1x,'n=',i3,1x,'atom=',i3,1x,
     &        'layer=',i3,1x,'energy=',f15.8)
99002 FORMAT ('#g: kappa=',i3,1x,'mu=',f5.1,1x,'n=',i3,1x,'atom=',i3,1x,
     &        'layer=',i3,1x,'energy=',f15.8)
C
99003 FORMAT (i3,1x,i3,2x,i2,2x,f4.1,2x,f4.1,2x,3F15.8)
      END
C*==readcore.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE READCORE(CWFX,CWFY,CWF,EVAL,LAYS,NATL,CWFSTATES,
     &                    CWFENERGY,CWFKAPMUE,ROT,ROTM1,EB,CHANNEL)
C     /****************************************************************/
C     # purpose      : read core wavefunctions and energies from file. *
C                      the wavefunctions are read from the file in     *
C                      which they are stored for the case of           *
C                      magnetic field in z-direction.                  *
C                      if the real magnetic field differs from that,   *
C                      the read wavefunctions are rotated before they  *
C                      are stored  in the cwf-field                    *
C                      (core-wavefunctions).                           *
C                                                                      *
C     # parameters   :                                                 *
C         cwfx,cwfy,cwf : pointer fields and field for storage of      *
C                         core-wavefunctions for every atom and layer  *
C         eval : storage for the energy-eigenvalues                    *
C         lays,natl: number of layers and atoms in layers              *
C         cwfstates: how many cwf-states are present for atom in layer *
C         cwfenergy: all eigenvalues for every atom in layer           *
C         cwfkapmue: kappa and mue values for the above eigenvalues    *
C                    (not really needed, just for information)         *
C         eua,eub,euc: euler angles for rotation (only needed for      *
C         rot,rotm1: rotation matrix and its inverse                   *
C     # calls the following functions:                                 *
C       simps        kapmue2k    kapmue2l                              *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,MAXEIG,MAXN,MAXCORE,
     &    MAXCSTATES,HARTRE,EPS12,CZERO,CONE
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER CHANNEL,LAYS
      COMPLEX*16 CWF(RSTEP,MAXCORE,2),ROT(MQD,MQD),ROTM1(MQD,MQD)
      REAL*8 CWFENERGY(MAXCSTATES,NATLM,LAYSM),
     &       CWFKAPMUE(MAXCSTATES,NATLM,LAYSM,2),EB(3),
     &       EVAL(MAXEIG,MAXN,NATLM,LAYSM)
      INTEGER CWFSTATES(NATLM,LAYSM),CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),
     &        CWFY(MAXCORE),NATL(LAYSM)
C
C Local variables
C
      INTEGER ACSTATE,ATOM,COLNZ(:),COUNTVAR,CT,CWFMAP(:,:,:,:,:),I,IN,
     &        IN1,IN2,J,K,KAPPA,L,LAY,LAYER,LINECOUNT,M,N1,N2,NEW,NZN,
     &        PT1,PT2,PTP,ROWNZ(:),RR,STATES,TATOM,TKAPPA,TLAYER
      COMPLEX*16 DUMWFF(:,:),DUMWFG(:,:),HA(:),HB(:),HC(:),WF1(:,:),
     &           WF2(:,:)
      REAL*8 EREL,MUE,SHIFT,TEREL,TMUE
      LOGICAL IB
      INTEGER KAPMUE2K,KAPMUE2L
      CHARACTER READFLAG
      EXTERNAL KAPMUE2K,KAPMUE2L
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE HA,HB,HC,WF1,WF2,COLNZ,ROWNZ,CWFMAP,DUMWFF,DUMWFG
      ALLOCATE (HA(MQD*MQD),HB(MQD*MQD),HC(MQD*MQD),WF1(MQD,MQD))
      ALLOCATE (WF2(MQD,MQD),COLNZ(MQD),ROWNZ(MQD))
      ALLOCATE (CWFMAP(MQD,MQD,MAXCSTATES,NATLM,LAYSM))
      ALLOCATE (DUMWFF(RSTEP,4),DUMWFG(RSTEP,4))
C
C     shift the energylevels by shift ev:
C     ===================================
Ccghf should be read from input, and added in calcore or writecore
C     shift=24.5
      SHIFT = 0.D0
      REWIND (CHANNEL)
C
      IF ( ABS(EB(1)).GT.EPS12 .OR. ABS(EB(2)).GT.EPS12 .OR. ABS(EB(3))
     &     .GT.EPS12 ) THEN
         IB = .TRUE.
      ELSE
         IB = .FALSE.
      END IF
C
C     init arrays:
C     ============
C
      DO I = 1,MQD
         DO J = 1,MQD
            DO K = 1,MAXCSTATES
               DO L = 1,NATLM
                  DO M = 1,LAYSM
                     CWFMAP(I,J,K,L,M) = 0
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,MAXEIG
         DO J = 1,MAXN
            DO K = 1,NATLM
               DO L = 1,LAYSM
                  EVAL(I,J,K,L) = 0.D0
               END DO
            END DO
         END DO
      END DO
C
      DO J = 1,MQD
         DO K = 1,MAXCSTATES
            DO L = 1,NATLM
               DO M = 1,LAYSM
                  CWFX(J,K,L,M,1) = 0
                  CWFX(J,K,L,M,2) = -1
               END DO
            END DO
         END DO
      END DO
C
      DO I = 1,RSTEP
         DO J = 1,MAXCORE
            CWF(I,J,1) = CZERO
            CWF(I,J,2) = CZERO
         END DO
      END DO
C
      DO J = 1,MAXCORE
         CWFY(J) = 0
      END DO
C
      DO I = 1,NATLM
         DO J = 1,LAYSM
            CWFSTATES(I,J) = 0
         END DO
      END DO
C
      LINECOUNT = 1
      STATES = 0
      COUNTVAR = 0
      ACSTATE = 0
C     start reading file:
C     ===================
 100  CONTINUE
      READ (CHANNEL,99002,END=300,ERR=200) READFLAG,KAPPA,MUE,N1,ATOM,
     &      LAYER,EREL
C
      LINECOUNT = LINECOUNT + 1
      IF ( ATOM.GT.NATLM ) THEN
         WRITE (*,*) 'in the core-wavefunction file an atom',
     &               ' number was specified'
         WRITE (*,*) 'which is larger than natlm!! aborting ...'
      END IF
C
      IF ( LAYER.GT.LAYSM ) THEN
         WRITE (*,*) 'in the core-wavefunction file a layer',
     &               ' number was specified'
         WRITE (*,*) 'which is larger than laysm!! aborting ...'
      END IF
C
      IF ( READFLAG.NE.'!' ) THEN
C         read f-core-wavefunction:
C     =========================
         DO RR = 1,RSTEP
            READ (CHANNEL,99004,END=300,ERR=200) (DUMWFF(RR,J),J=1,4)
         END DO
         READ (CHANNEL,99002,END=300,ERR=200) READFLAG,TKAPPA,TMUE,N2,
     &         TATOM,TLAYER,TEREL
CC
C         IF ( KAPPA.NE.TKAPPA .OR. MUE.NE.TMUE .OR. N1.NE.N2 .OR.
C     &        ATOM.NE.TATOM .OR. LAYER.NE.TLAYER .OR. EREL.NE.TEREL )
C     &        THEN
         IF ( KAPPA.NE.TKAPPA .OR. ABS(MUE-TMUE).GT.1.0D-16 .OR. 
     &        N1.NE.N2 .OR. ATOM.NE.TATOM .OR. LAYER.NE.TLAYER .OR. 
     &        ABS(EREL-TEREL).GT.1.0D-16 ) THEN
            WRITE (*,99001)
            STOP
         END IF
C         read g-core-wavefunction:
C         =========================
         DO RR = 1,RSTEP
            READ (CHANNEL,99004,END=300,ERR=200) (DUMWFG(RR,J),J=1,4)
         END DO
         CWFSTATES(ATOM,LAYER) = CWFSTATES(ATOM,LAYER) + 1
         IF ( CWFSTATES(ATOM,LAYER).GT.MAXCSTATES ) THEN
            WRITE (NOUT1,*) 'error in readcore; arrays to small to',
     &                      ' store all elements; increase maxcstates'
            WRITE (*,*) 'error in readcore; arrays to small to',
     &                  ' store all elements; increase maxcstates'
            WRITE (*,*) 'aborting ...'
            STOP
         END IF
         CWFENERGY(CWFSTATES(ATOM,LAYER),ATOM,LAYER) = (EREL-SHIFT)
     &      /HARTRE
         CWFKAPMUE(CWFSTATES(ATOM,LAYER),ATOM,LAYER,1) = KAPPA
         CWFKAPMUE(CWFSTATES(ATOM,LAYER),ATOM,LAYER,2) = MUE
C     wavefunction read, now rotate it:
C     =================================
         ACSTATE = ACSTATE + 1
         ACSTATE = KAPMUE2K(KAPPA,MUE,ML)
         PTP = KAPMUE2L(KAPPA,MUE)
         PT1 = KAPMUE2K(KAPPA,MUE,ML)
         PT2 = PT1
         IF ( IB .AND. -KAPPA-1.NE.0 .AND. ABS(MUE).LT.ABS(-KAPPA-1) )
     &        PT2 = KAPMUE2K(-KAPPA-1,MUE,ML)
C
         IF ( PTP.GT.MAXEIG ) THEN
            WRITE (NOUT1,*) 'error in readcore: maxeig not large',
     &                      ' enough to store core wavefunctions!'
            WRITE (*,*) 'error in readcore: maxeig not large',
     &                  ' enough to store core wavefunctions!'
            STOP
         END IF
C
         EVAL(PTP,N1,ATOM,LAYER) = (EREL-SHIFT)/HARTRE
C
C         first determine resulting wavefunction
C         components (that means: what components will
C         be present after rotation; therefore
C         fill testmatrix with 1 at the position,
C         where the wavefunction exists. then
C         look what other positions will be effected
C         during the rotation. this positions can be
C         possibly non zero in the resulting wavefunction).
C         this procedure is neccessary because the core-
C         wavefunction is stored as a vector and not as
C         a matrix. so the information on what and how many
C         elements have to be stored must be present when
C         preparing the storage-vectors:
C         =================================================
C
         DO I = 1,MQD
            DO J = 1,MQD
               WF1(I,J) = CZERO
               WF2(I,J) = CZERO
            END DO
         END DO
C
         WF1(PT1,PT1) = CONE
C         IF ( PT1.NE.PT2 .AND. CDABS(DUMWFF(200,4)).NE.0.D0 )
C     &        WF1(PT2,PT2) = CONE
         IF ( PT1.NE.PT2 .AND. CDABS(DUMWFF(200,4)).GT.1.0D-16 )
     &        WF1(PT2,PT2) = CONE
C         test what components are  effected in wf2=(wf1 x rot):
C         ======================================================
         DO I = 1,MQD
            DO J = 1,MQD
               DO K = 1,MQD
C                  IF ( CDABS(WF1(I,K)*ROT(K,J)).NE.0.D0 ) WF2(I,J)
C     &                 = CONE
                  IF ( CDABS(WF1(I,K)*ROT(K,J)).GT.0.D0 ) WF2(I,J)
     &                 = CONE
               END DO
            END DO
         END DO
C
         DO I = 1,MQD
            DO J = 1,MQD
               WF1(I,J) = 0.D0
            END DO
         END DO
C
C         test what components are  effected in wf1=(rotm1 x wf2):
C         ========================================================
         DO I = 1,MQD
            DO J = 1,MQD
               DO K = 1,MQD
C                  IF ( CDABS(ROTM1(I,K)*WF2(K,J)).NE.0.D0 ) WF1(I,J)
C     &                 = CONE
                  IF ( CDABS(ROTM1(I,K)*WF2(K,J)).GT.1.0D-16 ) WF1(I,J)
     &                 = CONE
               END DO
            END DO
         END DO
C         now prepare map with information where the matrix-elements
C         are stored in the vector; such a map is constructed for
C         every core-state of every atom in every layer:
C         ==========================================================
         DO IN1 = 1,MQD
            DO IN2 = 1,MQD
               CWFMAP(IN1,IN2,CWFSTATES(ATOM,LAYER),ATOM,LAYER) = 0
               IF ( CDABS(WF1(IN1,IN2)).GT.EPS12 ) THEN
                  STATES = STATES + 1
                  IF ( STATES.GT.MAXCORE ) THEN
                     WRITE (NOUT1,*) 'error while reading core file;',
     &                            ' too many wavefunctions; aborting...'
                     WRITE (*,*) 'error while reading core file;',
     &                           ' too many wavefunctions; aborting...'
                     STOP
                  END IF
                  CWFMAP(IN1,IN2,CWFSTATES(ATOM,LAYER),ATOM,LAYER)
     &               = STATES
               END IF
            END DO
         END DO
C         determine nonzero_row/columns for
C         multsparse from the cwfmap:
C         =================================
         CALL DETNZROWCOL(NZN,COLNZ,ROWNZ,
     &                    CWFMAP(1,1,CWFSTATES(ATOM,LAYER),ATOM,LAYER))
C
C         now determine rotated wavefunction:
C         ===================================
         DO RR = 1,RSTEP
C             rotate large dirac component:
C             =============================
            DO I = 1,MQD
               DO J = 1,MQD
                  WF1(I,J) = CZERO
                  WF2(I,J) = CZERO
               END DO
            END DO
C
            WF1(PT1,PT1) = DUMWFF(RR,1)
            IF ( PT1.NE.PT2 ) WF1(PT2,PT2) = DUMWFF(RR,4)
C
            CALL MULTSPARSE(NZN,COLNZ,ROWNZ,WF1,ROT,WF2,HA,HB,HC)
            CALL MULTSPARSE(NZN,COLNZ,ROWNZ,ROTM1,WF2,WF1,HA,HB,HC)
C
            DO IN1 = 1,MQD
               DO IN2 = 1,MQD
                  IN = CWFMAP(IN1,IN2,CWFSTATES(ATOM,LAYER),ATOM,LAYER)
                  IF ( IN.NE.0 ) CWF(RR,IN,1) = WF1(IN1,IN2)
                  WF1(IN1,IN2) = CZERO
               END DO
            END DO
C             rotate small dirac component:
C             =============================
            DO I = 1,MQD
               DO J = 1,MQD
                  WF1(I,J) = CZERO
                  WF2(I,J) = CZERO
               END DO
            END DO
C
            WF1(PT1,PT1) = DUMWFG(RR,1)
            IF ( PT1.NE.PT2 ) WF1(PT2,PT2) = DUMWFG(RR,4)
C
            CALL MULTSPARSE(NZN,COLNZ,ROWNZ,WF1,ROT,WF2,HA,HB,HC)
            CALL MULTSPARSE(NZN,COLNZ,ROWNZ,ROTM1,WF2,WF1,HA,HB,HC)
C
            DO IN1 = 1,MQD
               DO IN2 = 1,MQD
                  IN = CWFMAP(IN1,IN2,CWFSTATES(ATOM,LAYER),ATOM,LAYER)
                  IF ( IN.NE.0 ) CWF(RR,IN,2) = WF1(IN1,IN2)
                  WF1(IN1,IN2) = CZERO
               END DO
            END DO
         END DO
      ELSE
C         case of readflag='!'; skip this wavefunction:
C         =============================================
C         read dummy:
C         ===========
         DO I = 1,RSTEP
C            READ (CHANNEL,99006,END=300,ERR=200) DUM,DUM,DUM,DUM,DUM,
C     &            DUM,DUM,DUM,DUM
            LINECOUNT = LINECOUNT + 1
         END DO
         READ (CHANNEL,99002,END=300,ERR=200) READFLAG,KAPPA,MUE,N1,
     &         ATOM,LAYER,EREL
         LINECOUNT = LINECOUNT + 1
         DO I = 1,RSTEP
C            READ (CHANNEL,99006,END=300,ERR=200) DUM,DUM,DUM,DUM,DUM,
C     &           DUM,DUM,DUM,DUM
            LINECOUNT = LINECOUNT + 1
         END DO
      END IF
C
      GOTO 100
 200  CONTINUE
      WRITE (*,*) 'error in line ',LINECOUNT,' attempting to ',
     &            'read core data. aborting ...'
      STOP
C
 300  CONTINUE
      WRITE (NOUT1,*) '# read or rotated states:',STATES
C
C     now fill the other core index fields cwfx and cwfy:
C     ===================================================
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            DO CT = 1,CWFSTATES(ATOM,LAY)
               DO IN1 = 1,MQD
                  NEW = -1
                  DO IN2 = 1,MQD
                     IN = CWFMAP(IN1,IN2,CT,ATOM,LAY)
                     IF ( IN.NE.0 ) THEN
                        IF ( NEW.EQ.-1 ) THEN
                           CWFX(IN1,CT,ATOM,LAY,1) = IN
                           NEW = 0
                        END IF
                        NEW = NEW + 1
                        CWFY(IN) = IN2
                     END IF
                  END DO
                  CWFX(IN1,CT,ATOM,LAY,2) = CWFX(IN1,CT,ATOM,LAY,1)
     &               + NEW - 1
               END DO
            END DO
         END DO
      END DO
C
C     output to file:
C     ===============
      IF ( IP.GT.4 ) THEN
         WRITE (NOUT1,99005)
         J = 0
         DO I = 1,CWFSTATES(1,1)
            DO IN1 = 1,MQD
               DO IN = CWFX(IN1,I,1,1,1),CWFX(IN1,I,1,1,2)
                  IF ( I.NE.J ) THEN
                     J = I
                     COUNTVAR = 1
                  ELSE
                     COUNTVAR = COUNTVAR + 1
                  END IF
                  IN2 = CWFY(IN)
                  WRITE (NOUT1,99003) IN,I,COUNTVAR,IN1,IN2,
     &                                CWF(200,IN,1),CWF(200,IN,2)
               END DO
            END DO
            WRITE (NOUT1,*) '#'
         END DO
      END IF
C
      RETURN
C
99001 FORMAT ('in readcore: expecting to read g component of wf')
99002 FORMAT (1x,1a,8x,i3,4x,f5.1,3x,i3,6x,i3,7x,i3,8x,f15.8)
99003 FORMAT (5I4,2(1x,e15.8),1x,2(1x,e15.8))
99004 FORMAT (9(1x,e15.8))
99005 FORMAT ('core-wavefunctions after reading and probably rotation:')
      END
C*==detnzrowcol.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETNZROWCOL(NZN,COLNZ,ROWNZ,CWFMAP)
C     /****************************************************************/
C     # purpose      : determine which cols / rows are non zero        *
C                                                                      *
C     # parameters   :                                                 *
C          cwfmap: gives the non zero elements which will occur        *
C                  in the whole multiplication procedure               *
C                  ( [wf'=(r^-1 wf r)] has been calculated before )    *
C          nzn:    non zero number of rows (non zero number of         *
C                  cols should be the same)                            *
C          colnz/rownz: fields specify the non zero elements           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NZN
      INTEGER COLNZ(MQD),CWFMAP(MQD,MQD),ROWNZ(MQD)
C
C Local variables
C
      INTEGER FOUND,I,J,NZ_COL,NZ_ROW
C
C*** End of declarations rewritten by SPAG
C
      NZ_COL = 0
      NZ_ROW = 0
      DO I = 1,MQD
         COLNZ(I) = 0
         ROWNZ(I) = 0
      END DO
C
      DO I = 1,MQD
         FOUND = 0
         DO J = 1,MQD
C            IF ( ABS(CWFMAP(I,J)).NE.0.D0 ) FOUND = 1
            IF ( ABS(CWFMAP(I,J)).GT.1.0D-16 ) FOUND = 1
         END DO
         IF ( FOUND.EQ.1 ) THEN
            NZ_COL = NZ_COL + 1
            COLNZ(NZ_COL) = I
         END IF
      END DO
C
      DO I = 1,MQD
         FOUND = 0
         DO J = 1,MQD
C            IF ( ABS(CWFMAP(J,I)).NE.0.D0 ) FOUND = 1
            IF ( ABS(CWFMAP(J,I)).GT.1.0D-16 ) FOUND = 1
         END DO
         IF ( FOUND.EQ.1 ) THEN
            NZ_ROW = NZ_ROW + 1
            ROWNZ(NZ_ROW) = I
         END IF
      END DO
C
      IF ( NZ_COL.NE.NZ_ROW ) THEN
         WRITE (*,99001)
         WRITE (*,99002)
         STOP
      END IF
      NZN = NZ_COL
C
      RETURN
C
99001 FORMAT ('stopped: something is wrong in detnzrowcol,')
99002 FORMAT ('expecting quadratic matrices...')
      END
C*==multsparse.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE MULTSPARSE(NZN,COLNZ,ROWNZ,A,B,C,HA,HB,HC)
C     /****************************************************************/
C     # purpose      : multiply the core-wavefunction matrix with the  *
C                      rotation matrix. since the core-wavefunction    *
C                      matrix has just one block unequal to zero,      *
C                      this block is written to a smaller matrix to.   *
C                      increase speed the non zero blcok is specified  *
C                      by the fields colnz (columns not zero)          *
C                      and rownz (rows not zero).                      *
C                      this fields are calculated in detnzrowcol.      *
C                                                                      *
C     # parameters   :                                                 *
C          nzn: non zero number of rows (non zero number of            *
C               columns should be equal)                               *
C          colnz,rownz: specifiy which cols/rows are non zero          *
C          a,b,c: a x b = c  (matrix multiplication)                   *
C          ha,hb,hc: ha x hb = hc (reduced matrix multiplication)      *
C     # calls the subroutine:                                          *
C       mult                                                           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NZN
      COMPLEX*16 A(MQD,MQD),B(MQD,MQD),C(MQD,MQD),HA(MQD*MQD),
     &           HB(MQD*MQD),HC(MQD*MQD)
      INTEGER COLNZ(MQD),ROWNZ(MQD)
C
C Local variables
C
      INTEGER COUNTVAR,I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MQD*MQD
         HA(I) = CZERO
         HB(I) = CZERO
         HC(I) = CZERO
      END DO
C
C   fill reduced matrix (stored as a vector for implementation reasons):
C   ====================================================================
      COUNTVAR = 0
      DO J = 1,NZN
         DO I = 1,NZN
            COUNTVAR = COUNTVAR + 1
            HA(COUNTVAR) = A(COLNZ(I),ROWNZ(J))
            HB(COUNTVAR) = B(COLNZ(I),ROWNZ(J))
         END DO
      END DO
C
      CALL MULT(HA,HB,HC,NZN)
C
      DO I = 1,MQD
         DO J = 1,MQD
            C(I,J) = CZERO
         END DO
      END DO
      COUNTVAR = 0
      DO J = 1,NZN
         DO I = 1,NZN
            COUNTVAR = COUNTVAR + 1
            C(COLNZ(I),ROWNZ(J)) = HC(COUNTVAR)
         END DO
      END DO
C
      END
