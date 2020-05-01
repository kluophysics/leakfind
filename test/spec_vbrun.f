C*==spec_upsrun.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
      SUBROUTINE SPEC_UPSRUN(IBLOCH,NOUT1,IP,NREL,IREL,NSPIN,SPOL,
     &                       NSPLEE,EFEV,VPR,VIH,VIL,OMEV,OMHAR,THQ,FIQ,
     &                       ALQ,DELQ,ICIRC,IDREH,MCD,NATL,LAYS,LAYB,
     &                       LAYP,LAYER1,GANZ,LANZ1,LANZ2,POL0,POL0L,
     &                       TYP,ISTR,Q1,Q2,Q3,Q4,NEW,PKSCAN,PKEMIN,
     &                       PKEMAX,ESTEP,THETA,PHI,NT,NP,DELTAT,DELTAP,
     &                       TMINA,PMINA,IFM,POS,IRUMP,BRUSEP,NCLM,CELM,
     &                       CLIGHT,MAXL,NOD,NEV,NTOT,USEEULER,IBXY,Z,
     &                       ZPARU,ZPARD,EPSX,EBIND,PHOTOCSD,PHOTOC,
     &                       ROATLA,RMATO,RMATS,AA,AWR,ZMAT0,AMAT1X,
     &                       AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,AMAT2V,
     &                       AMAT2,PSI2G,BULKX,IGPRO,IGVAL,TSEL,TSOL,
     &                       TSEOL,TSOEL,TSEH,TSOH,TSEOH,TSOEH,IEDLAYER,
     &                       IODLAYER,IEODLAYER,IOEDLAYER,ZMESH,ADU,ADD,
     &                       BK,IFSP,EFERM,RADMESHM,VLM,MEZJ,MEZZ,MSSQ,
     &                       SSST,TAUQ,TAUT,PHASK,FIXTHQ,PSPIN,
     &                       POL0_VEC_TILT,POL0_INITIAL)
C
C   ********************************************************************
C   *                                                                  *
C   *  Driver to calculate:                                            *
C   *  IBLOCH = 0: band structure calculations                         *
C   *           1: SP-LEED                                             *
C   *           2: Angle integrated PES                                *
C   *           4: Angle resolved PES                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CPA,ONLY:NCPA
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NQMAX
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_MPI,ONLY:MPI,MPI_ID,NPROCS
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,MLS,LL,MAXSPEC,XMAXE,NFULLPOT,
     &    MLQNAT,MEW,MPW,LH,MHDIM2,MLZP,HARTRE,PI,NTPHOMAX,NVFTPHOMAX
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT,IBLOCHTMP
      USE MOD_CALCMODE,ONLY:DMFT,LDAU
      USE MOD_SPEC_ESIGMA,ONLY:EREAL
      USE MOD_THERMAL,ONLY:NVFO_Q
      USE MOD_FILES,ONLY:IFILSPECOU1,IFILSPECOU2
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALQ,CLIGHT,DELQ,DELTAP,DELTAT,EFERM,EFEV,EPSX,ESTEP,FIQ,
     &       FIXTHQ,OMEV,OMHAR,PHI,PKEMAX,PKEMIN,PMINA,THETA,THQ,TMINA,
     &       VIH,VIL,VPR
      INTEGER GANZ,IBLOCH,IBXY,ICIRC,IDREH,IFM,IFSP,IP,IREL,LANZ1,LANZ2,
     &        LAYB,LAYER1,LAYP,LAYS,MAXL,MCD,NCLM,NEW,NOUT1,NP,NREL,
     &        NSPIN,NSPLEE,NT,PKSCAN,SPOL,TYP,USEEULER
      LOGICAL POL0_VEC_TILT
      COMPLEX*16 Q1,Q2,Q3,Q4
      COMPLEX*16 AA(3),AMAT1(XMAXE,4,NFULLPOT,2),AMAT2(XMAXE,4,NFULLPOT)
     &           ,AWR(MLS,LL,NATLM,2),BK(LH,2),BULKX(LH,LH,LAYSM),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),PHASK(NEMAX),
     &           PHOTOCSD(NEW,NT,NP,4,2),PSI2G(MHDIM2,LH,4),
     &           RMATO(MQD,MQD,LAYSM,NATLM),RMATS(MQD,MQD,LAYSM,NATLM),
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TSEH(MLQNAT,MLQNAT,LAYSM),
     &           TSEL(MLQNAT,MLQNAT,LAYSM),TSEOH(MLQNAT,MLQNAT,LAYSM),
     &           TSEOL(MLQNAT,MLQNAT,LAYSM),TSOEH(MLQNAT,MLQNAT,LAYSM),
     &           TSOEL(MLQNAT,MLQNAT,LAYSM),TSOH(MLQNAT,MLQNAT,LAYSM),
     &           TSOL(MLQNAT,MLQNAT,LAYSM),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX),
     &           ZMAT0(MQD,MQD,LAYSM,NATLM)
      REAL*8 ADD(3),ADU(3),BRUSEP(3),CELM(MLZP),EBIND(MEW),
     &       PHOTOC(NEW,NT,NP,4,2),POL0(3),POL0L(3),POL0_INITIAL(3),
     &       POS(3,NATLM,LAYSM),PSPIN(3),RADMESHM(NRMAX,NATLM,LAYSM),
     &       ROATLA(LL),Z(NATLM,LAYSM),ZMESH(6),ZPARD(3),ZPARU(3)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IGPRO(LH),IGVAL(LH),
     &        IODLAYER(0:LAYSM),IOEDLAYER(0:LAYSM),IRUMP(LAYSM),ISTR(2),
     &        NATL(LAYSM),NEV(LAYSM),NOD(LAYSM),NTOT(LAYSM)
C
C Local variables
C
      REAL*8 AAZ,ALLQ,CALC_E(:),DET,DETF(:,:,:),DETI,DETR,DMFTMIX,EKHAR,
     &       KPARAF(:,:),KXF(:,:),KYF(:,:),PHID,PHIQ,PHIQ_INP,PKESTEP,
     &       PSPIN_INP(3),ROALL(2),ROALLF(:,:,:),ROATLAINC(:,:),ROINC(:)
     &       ,ROLAYRES(:,:,:,:,:),ROSULA(:),ROSUR,ROTELA(:),ROTER,ROTRA,
     &       ROTRLA(:),RWORK(:,:,:),RWORK1(:),TESTVAC,THEQ,THETAD,
     &       THETAF(:),THQ_INP,TMP
      COMPLEX*16 AMAT1T(:,:,:,:),AMAT2T(:,:,:),CPAPROJ(:,:,:,:,:,:),
     &           EHIGH,ELOW,ESTAT,ETAB(:),GAMMA1(:,:,:,:,:),K1,
     &           MSSTS(:,:,:,:),P(2),PF,PFI,ROALLSD(:),ROATLAINCSD(:,:),
     &           ROATLASD(LL),ROATSD(:),ROINCSD(:,:),
     &           ROLAYRESSD(:,:,:,:,:,:),ROSULASD(:),ROSURSD(:),
     &           ROTELASD(:),ROTERSD(:),ROTRASD(:),ROTRLASD(:),
     &           TAUMT(:,:,:,:,:),TSSTS(:,:,:,:),UMAT_VT_PES(:,:,:,:,:),
     &           ZMAT0INC(:,:,:,:,:)
      INTEGER ATA,ATOM,CPAATOM(:,:),FINAL,I,I2,IBPOL,IDRQ,IE,IED,IEOD,
     &        IERR,IFILE,IGP,INC,INITIAL,IOD,IOED,IPHPOL,IPL,IPOL,IPROC,
     &        IPROCA(:),IPROCE(:),IQ,IRSTATE,IRUM,ISDML,ISTATE,J,JA,JE,
     &        JP,JT,K,LAY,NTMP,PHSTATE
      LOGICAL CALCDOS,TILT
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
      ALLOCATABLE ROATSD,ROTRASD,ROTERSD,ROSURSD,ROALLSD
      ALLOCATABLE IPROCE,THETAF,ROALLF,RWORK,IPROCA
      ALLOCATABLE ROATLAINC,ROATLAINCSD,ROINC,ROINCSD
      ALLOCATABLE CPAPROJ
      ALLOCATABLE TAUMT,ZMAT0INC,ETAB
      ALLOCATABLE RWORK1
      ALLOCATABLE UMAT_VT_PES
      ALLOCATABLE GAMMA1,AMAT1T,AMAT2T,CALC_E
      ALLOCATABLE CPAATOM
      ALLOCATABLE ROSULA,ROTELA,ROTRLA
      ALLOCATABLE KPARAF,DETF
      ALLOCATABLE ROLAYRES,ROTRLASD,ROTELASD,ROSULASD,ROLAYRESSD
      ALLOCATABLE TSSTS,MSSTS,KXF,KYF
C ======================================================================
      ALLOCATE (ROSULASD(LL),ROTELASD(LL),ROTRLASD(LL),ROALLSD(4))
      ALLOCATE (ROSULA(LL),ROTELA(LL),ROTRLA(LL))
      ALLOCATE (ROATSD(4),ROTRASD(4),ROTERSD(4),ROSURSD(4))
C
      ALLOCATE (GAMMA1(LAYSM,NATLM,3,MQD,MQD))
      ALLOCATE (AMAT1T(XMAXE,4,NFULLPOT,2),AMAT2T(XMAXE,4,NFULLPOT))
      ALLOCATE (CALC_E(MAXSPEC),CPAATOM(LAYSM,NATLM))
C
      ALLOCATE (TSSTS(NKMMAX,NKMMAX,NTMAX,2))
      ALLOCATE (MSSTS(NKMMAX,NKMMAX,NTMAX,2))
C
      ALLOCATE (ETAB(MEW),ROLAYRESSD(1,NT,NP,4,2,LAYP))
      ALLOCATE (ROLAYRES(NT,NP,2,2,LAYP))
C ======================================================================
C
C
      ROTRLA(:) = 0.0D0
      ROTELA(:) = 0.0D0
      ROSULA(:) = 0.0D0
      ROTRLASD(:) = C0
      ROTELASD(:) = C0
      ROSULASD(:) = C0
C
      DETR = 0.0D0
      DETI = 0.0D0
      ROLAYRES = 0.0D0
      ROLAYRESSD = C0
      PHOTOC = 0.0D0
      PHOTOCSD = C0
C
      GAMMA1(:,:,:,:,:) = C0
      CPAATOM(:,:) = 0
C
      TILT = .FALSE.
      THETAD = 0.0D0
      PHID = 0.0D0
C
      INITIAL = 1
      FINAL = 2
      IRSTATE = 3
      CALCDOS = .FALSE.
C
      IBPOL = 1
      IPHPOL = 1
      ALLQ = ALQ
      IDRQ = IDREH
      PSPIN_INP(1:3) = PSPIN(1:3)
      THQ_INP = THQ
      PHIQ_INP = FIQ
C
      IF ( IBLOCH.EQ.2 ) WRITE (6,99002) 'ANGLE INTEGRATED PES'
C
C ======================================================================
C
C        ATA - parameter; if ATA=1 we have the ATA-approximation for the
C                                  initial state of the photocurrent
C
C        ATA = 1
C        ATAVIB - parameter; if ATAVIB=1 we calculate the Debye-Waller
C                                        like photocurrent
C        ATAVIB - parameter; if ATAVIB=2 we calculate the direct type
C                                        photocurrent
C ======================================================================
      ATA = 0
C      ATAVIB = 0
C      ATAVIB = 1
      IF ( IBLOCH.EQ.0 ) ATA = 1
C ==========================================  Init of tilt calculations
      CALL SPEC_TILT(THETAD,PHID,TILT,0,NP,POL0,POL0_VEC_TILT,
     &               POL0_INITIAL,IP)
C
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C                                            initialize MPI - parameters
C
C
      ALLOCATE (IPROCE(MEW))
      ALLOCATE (IPROCA(NT*NP))
      ALLOCATE (THETAF(MPW),ROALLF(MPW,2,MEW),RWORK(MPW,2,MEW))
      ALLOCATE (KPARAF(MPW,MEW),DETF(MPW,2,MEW))
      ALLOCATE (KXF(MPW,MPW),KYF(MPW,MPW))
C
      IPROCE(:) = 0
      IPROCA(:) = 0
      KXF = 0.0D0
      KYF = 0.0D0
C
      IF ( NEW.LE.1 ) THEN
C
C     MPI over angular loops
C
         IPROCE(:) = MPI_ID
         IF ( NPROCS.GT.1 ) THEN
            IPROC = 0
            INC = 1
            DO JE = 1,NT*NP - 1
               IE = JE
               IPROC = IPROC + INC
               IF ( IPROC.EQ.NPROCS ) THEN
                  IPROC = NPROCS - 1
                  INC = -1
               ELSE IF ( IPROC.EQ.-1 ) THEN
                  IPROC = 0
                  INC = 1
               END IF
               IPROCA(IE) = IPROC
            END DO
         END IF
      ELSE IF ( NPROCS.GT.1 ) THEN
C
C     MPI over energy loop
C
         IPROCA(:) = MPI_ID
         IPROC = 0
         INC = 1
         DO JE = 1,NEW - 1
            IE = JE
            IPROC = IPROC + INC
            IF ( IPROC.EQ.NPROCS ) THEN
               IPROC = NPROCS - 1
               INC = -1
            ELSE IF ( IPROC.EQ.-1 ) THEN
               IPROC = 0
               INC = 1
            END IF
            IPROCE(IE) = IPROC
         END DO
      END IF
      IF ( MPI ) THEN
         THETAF(:) = 0.0D0
         ROALLF(:,:,:) = 0.0D0
         DETF(1:MPW,1:2,1:MEW) = 0.0D0
         KPARAF(1:MPW,1:MEW) = 0.0D0
      END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C
C ======================================================================
C     calculate energy vector, that is the vector which specifies
C     where to calculate the spectra:
C ======================================================================
      SELECT CASE (PKSCAN)
      CASE (0)
         DO JE = 1,NEW
            CALC_E(JE) = DBLE(JE-1)*ESTEP + EFEV
         END DO
      CASE DEFAULT
         IF ( NEW.GT.1 ) THEN
            PKESTEP = (PKEMAX-PKEMIN)/DBLE(NEW-1)
            DO JE = 1,NEW
               CALC_E(JE) = PKEMIN + DBLE(JE-1)*PKESTEP + EFEV
            END DO
         ELSE
            CALC_E(1) = PKEMIN + EFEV
         END IF
      END SELECT
C
C ======================================================================
      IF ( NSPLEE.NE.0 ) THEN
         THQ = THQ_INP*PI/180.0D0
         FIQ = PHIQ_INP*PI/180.0D0
         CALL SPEC_TILT(THQ,FIQ,TILT,1,NP,POL0,POL0_VEC_TILT,
     &                  POL0_INITIAL,IP)
         THQ = THQ*180.0D0/PI
         FIQ = FIQ*180.0D0/PI
         CALL PHOTONVEC(AA,THQ,FIQ,ALLQ,DELQ,OMHAR,ICIRC,IDRQ)
      END IF
      AAZ = DBLE(AA(3))
C ======================================================================
C
C
C=======================================================================
C     DMFT-calculational mode
C=======================================================================
C
      IBLOCHTMP = IBLOCH
      IF ( DMFT .OR. LDAU ) THEN
         DO JE = 1,NEW
            EREAL = CALC_E(JE)
            CALL SPEC_SETVIL(VIL,EREAL,EFERM)
            ELOW = DCMPLX(EREAL,VIL)/HARTRE
C
            IF ( IBLOCH.EQ.1 ) THEN
               EHIGH = DCMPLX(EREAL+VPR*HARTRE,VIH)/HARTRE
               ELOW = EHIGH
            END IF
C
            ETAB(JE) = 2.0D0*ELOW
         END DO
         CALL INIT_MOD_DMFT_LDAU(NEW,ETAB,DMFTMIX,NEW)
      END IF
C=======================================================================
C     Initialise temperature and CPA
C=======================================================================
C
      TMP = 0D0
      NTMP = 1
C
      CALL THERMAL_INIT(0,NTMP,TMP)
      CALL THERMAL_INIT(1,NTMP,TMP)
      CALL INIT_MOD_CPA(NCPA)
C
      IF ( NTMP.GT.1 ) STOP '<SPEC_VBRUN>: only one temperature allowed'
C
C     Calculate number of types needed for
C     Photoemission with THERMAL_VIBRA_FLUCT
C
      NVFTPHOMAX = 0
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            IQ = IQ_SPR_LAYAT(ATOM,LAY)
            NVFTPHOMAX = MAX(NVFTPHOMAX,NVFO_Q(IQ))
         END DO
      END DO
      ALLOCATE (UMAT_VT_PES(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX))
      ALLOCATE (ROATLAINC(LL,NVFTPHOMAX+1),ROINC(NVFTPHOMAX+1))
      ALLOCATE (ROATLAINCSD(LL,NVFTPHOMAX+1))
      ALLOCATE (ROINCSD(NVFTPHOMAX+1,4))
      ALLOCATE (CPAPROJ(LAYSM,NATLM,MQD,MQD,2,NVFTPHOMAX))
      ALLOCATE (TAUMT(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX+1))
      ALLOCATE (ZMAT0INC(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX+1))
C
      UMAT_VT_PES = C0
      CPAPROJ = C0
      TAUMT = C0
C
C=======================================================================
C
C*********************************************************** ENERGY LOOP
      LOOP_ENERGY:DO JE = 1,NEW
         EREAL = CALC_E(JE)
C
         CALL SPEC_SETVIL(VIL,EREAL,EFERM)
C
         ELOW = DCMPLX(EREAL,VIL)/HARTRE
         EHIGH = DCMPLX(EREAL+OMEV,VIH)/HARTRE
         IF ( IBLOCH.EQ.1 ) EHIGH = DCMPLX(EREAL+VPR*HARTRE,VIH)/HARTRE
         TESTVAC = 1.0D0
         IF ( IBLOCH.GE.2 ) THEN
            TESTVAC = REAL(EHIGH) - VPR
            IF ( TESTVAC.LT.0.0D0 ) WRITE (6,*) 'ELECTRON CANNOT ESCAPE'
         END IF
         EBIND(JE) = EREAL - EFERM
         EKHAR = DBLE(EHIGH) - VPR
         IF ( IBLOCH.LE.1 ) EKHAR = DBLE(ELOW)
         IF ( IP.GT.0 ) WRITE (6,*) 'VPR,EF,EKIN',VPR*HARTRE,DBLE(EHIGH)
     &                              *HARTRE,EKHAR*HARTRE
C
C******************************************************IF MPI ENERGY
         IF ( IPROCE(JE).EQ.MPI_ID ) THEN
C
            IFM = 0
            CALL POTFULLM2(LAYS,NATL,Z,IFM,VLM)
C
            IF ( IBLOCH.LE.1 ) THEN
               I2 = 1
            ELSE
               I2 = 2
            END IF
C
            DO ISTATE = 1,I2
               IF ( TESTVAC.GT.0.0D0 ) THEN
C
                  IF ( ISTATE.EQ.1 .AND. IBLOCH.NE.1 ) THEN
                     ESTAT = ELOW
                     PF = SQRT(2.0D0*ELOW+2.0D0*ELOW*2.0D0*ELOW/(CLIGHT*
     &                    CLIGHT))
                     PFI = PF
                  ELSE IF ( ISTATE.EQ.1 .AND. IBLOCH.EQ.1 ) THEN
                     ESTAT = EHIGH
                     PF = SQRT
     &                    (2.0D0*EHIGH+2.0D0*EHIGH*2.0D0*EHIGH/(CLIGHT*
     &                    CLIGHT))
                  ELSE IF ( ISTATE.EQ.2 ) THEN
                     ESTAT = EHIGH
                     PF = SQRT
     &                    (2.0D0*EHIGH+2.0D0*EHIGH*2.0D0*EHIGH/(CLIGHT*
     &                    CLIGHT))
                  END IF
C
                  IF ( IBLOCH.EQ.2 .AND. ISTATE.EQ.1 ) THEN
                     CALCDOS = .TRUE.
                  ELSE
                     CALCDOS = .FALSE.
                  END IF
C
                  CALL SPEC_TAU_PREPARE(ISTATE,JE,NEW,LAYS,NATL,ESTAT,
     &                                  UMAT_VT_PES,TAUMT,TSSTS,MSSTS,
     &                                  GAMMA1,CPAPROJ,CPAATOM,P,MEZZ,
     &                                  MEZJ,MSSQ,SSST,TAUQ,TAUT,PHASK,
     &                                  NOUT1,CALCDOS,ATA,IP,PF,IBLOCH)
               END IF
            END DO
C
C
            IF ( IBLOCH.LE.1 ) GAMMA1(:,:,2,:,:) = GAMMA1(:,:,1,:,:)
C******************************************************   ANGULAR LOOPS
            JA = 0
            LOOP_PHI:DO JP = 1,NP
               LOOP_THETA:DO JT = 1,NT
                  JA = JA + 1
C******************************************************IF MPI ANGLES
                  IF ( IPROCA(JA).EQ.MPI_ID ) THEN
C
                     THETAD = (TMINA+DELTAT*DBLE(JT-1))
                     PHID = (PMINA+DELTAP*DBLE(JP-1))
C
                     CALL SPEC_TILT(THETAD,PHID,TILT,JP,NP,PSPIN,.TRUE.,
     &                              PSPIN_INP,IP)
C=======================================================================
C     PSPIN IS PER DEFAULT NOT ROTATED (IF NEEDED COMMENT OUT)
C=======================================================================
C
                     PSPIN = PSPIN_INP
C
                     THETA = THETAD
                     PHI = PHID
                     IF ( TESTVAC.GT.0.0D0 ) THEN
C
                        IF ( IFSP.EQ.1 ) THEN
                           THEQ = THETA*180.0/PI - FIXTHQ
                           PHIQ = PHIQ_INP*PI/180.0D0
                           CALL SPEC_TILT(THEQ,PHIQ,TILT,JP,NP,POL0,
     &                        POL0_VEC_TILT,POL0_INITIAL,IP)
                           THEQ = THEQ*180.0D0/PI
                           PHIQ = PHIQ*180.0D0/PI
                           CALL PHOTONVEC(AA,THEQ,PHIQ,ALLQ,DELQ,OMHAR,
     &                        ICIRC,IDRQ)
                        ELSE
                           THQ = THQ_INP*PI/180.0D0
                           FIQ = PHIQ_INP*PI/180.0D0
                           CALL SPEC_TILT(THQ,FIQ,TILT,JP,NP,POL0,
     &                        POL0_VEC_TILT,POL0_INITIAL,IP)
                           THQ = THQ*180.0D0/PI
                           FIQ = FIQ*180.0D0/PI
                           CALL PHOTONVEC(AA,THQ,FIQ,ALLQ,DELQ,OMHAR,
     &                        ICIRC,IDRQ)
                        END IF
C
C=======================================================================
C     Calculation of matrix elements
C=======================================================================
C
                        IF ( IBLOCH.GT.1 )
     &                       CALL SPEC_MECALC(RMATO,RMATS,ZMAT0,
     &                       ZMAT0INC,TSSTS,MSSTS,P,GAMMA1,TAUMT,
     &                       UMAT_VT_PES,CPAPROJ,CPAATOM,VLM,AA,AMAT1X,
     &                       AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,AMAT2V,
     &                       AMAT2,AMAT1T,AMAT2T,PFI,CLIGHT,OMHAR,IP,
     &                       IBLOCH,EHIGH,ELOW,USEEULER,JE,ATA,LAYS,
     &                       NATL,INITIAL,FINAL,IRSTATE,NOUT1,RADMESHM)
C=======================================================================
                     END IF
C
C=======================================================================
C     convert scattering phases to even/odd form
C=======================================================================
C
                     DO LAY = 1,LAYS
                        DO ATOM = 1,NATL(LAY)
                           DO PHSTATE = 1,2
                              CALL CKMULM(MLQNAT,TSEL,TSOL,TSEOL,TSOEL,
     &                           TSEH,TSOH,TSEOH,TSOEH,PHSTATE,LAY,ATOM,
     &                           GAMMA1)
                           END DO
                        END DO
                     END DO
C
                     IF ( SPOL.EQ.1 ) THEN
                        ISDML = 1
                     ELSE IF ( SPOL.GE.2 ) THEN
                        ISDML = 2
                     END IF
                     DO IPL = 1,ISDML
                        IF ( IBLOCH.EQ.0 ) EHIGH = ELOW
                        CALL FINALSTATE(GANZ,LANZ1,NSPIN,NREL,EHIGH,PHI,
     &                                  THETA,JE,VPR,EKHAR,NSPLEE,
     &                                  IBLOCH,IREL,CLIGHT,CELM,TYP,
     &                                  ISTR,POL0,POL0L,IPL,NATL,POS,
     &                                  IRUMP,LAYS,NOD,NEV,NTOT,BRUSEP,
     &                                  BULKX,NCLM,LAYB,TSEH,TSOH,TSEOH,
     &                                  TSOEH,IBXY,ELOW,EPSX,BK,ADU,ADD,
     &                                  IGP,IGPRO,ZMESH,Q1,Q2,Q3,Q4,
     &                                  LAYER1,MAXL,K1,AWR,LAYP,IFM,
     &                                  PSI2G,IGVAL,IOD,IED,IOED,IEOD,
     &                                  VIH,ZPARU,ZPARD,IEDLAYER,
     &                                  IODLAYER,IEODLAYER,IOEDLAYER,JT,
     &                                  NT,SPOL)
                     END DO
C
                     DO J = 1,MHDIM2
                        DO K = 1,LH
                           PSI2G(J,K,4) = PSI2G(J,K,1)
                           PSI2G(J,K,3) = PSI2G(J,K,2)
                        END DO
                     END DO
C
C******************************************************   ANGULAR LOOPS
                     LOOP_SPOL:DO IPOL = 1,SPOL
                        IF ( TESTVAC.GT.0.0D0 ) THEN
                           IF ( IBLOCH.EQ.0 ) EHIGH = ELOW
                           IF ( IBLOCH.GT.1 ) THEN
                              IF ( IBLOCH.EQ.2 .OR. IBLOCH.EQ.4 )
     &                             CALL INITIALSTATE(GANZ,LANZ2,NSPIN,
     &                             NREL,ELOW,PHI,THETA,JE,VPR,EKHAR,
     &                             NSPLEE,IBLOCH,IREL,CLIGHT,CELM,TYP,
     &                             ISTR,POL0,POL0L,IPOL,NATL,POS,IRUMP,
     &                             LAYS,NOD,NEV,NTOT,BRUSEP,BULKX,NCLM,
     &                             LAYB,TSEL,TSOL,TSEOL,TSOEL,IBXY,
     &                             EHIGH,LAYER1,K1,AWR,ATOM,LAY,MAXL,
     &                             IRUM,IFM,RMATO,RMATS,LAYP,OMHAR,
     &                             ZMESH,AAZ,ADU,ADD,BK,EPSX,IGVAL,
     &                             IGPRO,PSI2G,IGP,VIL,ZPARU,ZPARD,
     &                             IEDLAYER,IODLAYER,IEODLAYER,
     &                             IOEDLAYER,ROSUR,ROSULA,ROTER,ROTELA,
     &                             ROTRA,ROTRLA,DETR,DETI,JT,SPOL,
     &                             ROTRASD,ROTERSD,ROSURSD,ROTRLASD,
     &                             ROTELASD,ROSULASD)
C
                              IF ( IBLOCH.GT.1 ) THEN
                                 CALL XPSPHOTO(CLIGHT,OMHAR,ELOW,AWR,
     &                              ZMAT0,NATL,LAYS,LAYER1,ROATLA,1,
     &                              LAYB,IPOL,SPOL,ROATLASD)
C
                                 IF ( ATA.LE.1 .OR. IBLOCH.EQ.2 ) THEN
                                    DO I = 1,LAYP
                                       DO J = 1,NVFTPHOMAX + 1
                                         ROATLAINC(I,J) = 0.0D0
                                         ROATLAINCSD(I,J) = C0
                                       END DO
                                    END DO
                                    CALL XPSPHOTOINC(CLIGHT,ELOW,AWR,
     &                                 ZMAT0INC,NATL,LAYS,LAYER1,
     &                                 ROATLAINC,1,LAYB,IPOL,SPOL,
     &                                 ROATLAINCSD)
                                 END IF
                              END IF
                              CALL SPEC_SUMUP_INT(NEW,JE,JP,JT,EBIND,
     &                           EHIGH,IBLOCH,LAYP,IPOL,IBPOL,IPHPOL,
     &                           MCD,NT,NP,LL,NOUT1,IP,SPOL,ROTRA,ROTER,
     &                           ROSUR,ROATLA,ROINC,ROATSD,ROTRASD,
     &                           ROTERSD,ROSURSD,ROINCSD,ROTRLA,ROTELA,
     &                           ROSULA,ROTRLASD,ROTELASD,ROSULASD,
     &                           ROLAYRES,ROLAYRESSD,ROATLASD,ROATLAINC,
     &                           ROATLAINCSD,ROALL,ROALLSD,NCPA,ATA,
     &                           PSPIN,THETA,PHI,KXF,KYF,THETAF,KPARAF,
     &                           DETF,ROALLF,VPR,PHOTOC,TESTVAC,DETR,
     &                           DETI,TMINA,DELTAT,1)
                           END IF
                        ELSE IF ( TESTVAC.LT.0.0D0 ) THEN
                           ROALL(IPOL) = 0.0D0
                           ROALLSD(IPOL) = C0
                        END IF
C
                     END DO LOOP_SPOL
C******************************************************   LOOP SPOL
                     CALL SPEC_SUMUP_INT(NEW,JE,JP,JT,EBIND,EHIGH,
     &                  IBLOCH,LAYP,IPOL,IBPOL,IPHPOL,MCD,NT,NP,LL,
     &                  NOUT1,IP,SPOL,ROTRA,ROTER,ROSUR,ROATLA,ROINC,
     &                  ROATSD,ROTRASD,ROTERSD,ROSURSD,ROINCSD,ROTRLA,
     &                  ROTELA,ROSULA,ROTRLASD,ROTELASD,ROSULASD,
     &                  ROLAYRES,ROLAYRESSD,ROATLASD,ROATLAINC,
     &                  ROATLAINCSD,ROALL,ROALLSD,NCPA,ATA,PSPIN,THETA,
     &                  PHI,KXF,KYF,THETAF,KPARAF,DETF,ROALLF,VPR,
     &                  PHOTOC,TESTVAC,DETR,DETI,TMINA,DELTAT,2)
                  END IF
C******************************************************IF MPI ANGLES
               END DO LOOP_THETA
            END DO LOOP_PHI
C******************************************************ANGULAR LOOPS
         END IF
C******************************************************IF MPI ENERGY
      END DO LOOP_ENERGY
C******************************************************ENERGY LOOP
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C=======================================================================
C     Collect all results and write out into spc files
C=======================================================================
C
      IF ( IBLOCH.GT.1 ) THEN
         IF ( MPI ) THEN
            INC = MEW*MPW*2
            DO JE = 1,MEW
               DO JT = 1,MPW
                  RWORK(JT,1,JE) = 0.0D0
                  RWORK(JT,2,JE) = 0.0D0
               END DO
            END DO
C
            IF ( NEW.GT.1 ) THEN
               CALL DRV_MPI_REDUCE_R(ROALLF,RWORK(1,1,1),INC)
               RWORK = 0.0D0
               CALL DRV_MPI_REDUCE_R(DETF,RWORK(1,1,1),INC)
               RWORK = 0.0D0
               INC = MEW*MPW*1
               CALL DRV_MPI_REDUCE_R(KPARAF,RWORK(1,1,1),INC)
            ELSE
               INC = MPW*MPW
               ALLOCATE (RWORK1(INC))
               RWORK1 = 0.0D0
               CALL DRV_MPI_REDUCE_R(KXF,RWORK1(1),INC)
               RWORK1 = 0.0D0
               CALL DRV_MPI_REDUCE_R(KYF,RWORK1(1),INC)
               DEALLOCATE (RWORK1)
               INC = NEW*NT*NP*2*2
               ALLOCATE (RWORK1(INC))
               RWORK1 = 0.0D0
               CALL DRV_MPI_REDUCE_R(PHOTOC,RWORK1(1),INC)
               DEALLOCATE (RWORK1)
               INC = NT*NP*2*2*LAYP
               ALLOCATE (RWORK1(INC))
               RWORK1 = 0.0D0
               CALL DRV_MPI_REDUCE_R(ROLAYRES,RWORK1(1),INC)
            END IF
C
         END IF
         IF ( MPI_ID.EQ.0 ) THEN
            IF ( IPHPOL.EQ.1 ) THEN
               IFILE = IFILSPECOU1
            ELSE IF ( IPHPOL.EQ.2 ) THEN
               IFILE = IFILSPECOU2
            END IF
            DO JE = 1,NEW
               DO JP = 1,NP
                  DO JT = 1,NT
                     IF ( NEW.GT.1 ) THEN
                        DET = DETF(JT,1,JE)**2 + DETF(JT,2,JE)**2
                        WRITE (IFILE,99001) THETAF(JT),EBIND(JE),
     &                         (ROALLF(JT,1,JE)+ROALLF(JT,2,JE))/2.D0,
     &                         0.5*ROALLF(JT,1,JE),0.5*ROALLF(JT,2,JE),
     &                         (ROALLF(JT,1,JE)-ROALLF(JT,2,JE))
     &                         /(ROALLF(JT,1,JE)+ROALLF(JT,2,JE)),
     &                         KPARAF(JT,JE),DET
                     ELSE
                        WRITE (IFILE,99001) KXF(JT,JP),KYF(JT,JP),
     &                         (PHOTOC(JE,JT,JP,1,IBPOL)
     &                         +PHOTOC(JE,JT,JP,2,IBPOL))/2.D0,
     &                         0.5*PHOTOC(JE,JT,JP,1,IBPOL),
     &                         0.5*PHOTOC(JE,JT,JP,2,IBPOL),
     &                         (PHOTOC(JE,JT,JP,1,IBPOL)
     &                         -PHOTOC(JE,JT,JP,2,IBPOL))
     &                         /(PHOTOC(JE,JT,JP,1,IBPOL)
     &                         +PHOTOC(JE,JT,JP,2,IBPOL)),THETAF(JT)
                     END IF
                  END DO
               END DO
            END DO
            CALL FLUSH(IFILE)
         END IF
      END IF
C=======================================================================
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      IF ( MPI ) CALL MPI_FINALIZE(IERR)
C
      RETURN
C
99001 FORMAT (2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,
     &        e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5)
99002 FORMAT (//,1X,79('#'),/,10X,A,/,1X,79('#'),/)
C
      END
C*==spec_setvil.f    processed by SPAG 6.70Rc at 13:13 on 10 Mar 2017
C
C
      SUBROUTINE SPEC_SETVIL(VIL,EREAL,EFERM)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EFERM,EREAL,VIL
C
C*** End of declarations rewritten by SPAG
C
C
C       vil=0.05+0.002*(ereal-eferm)**2
C       vil=atan(0.05+0.03*(ereal-eferm)**2+0.5*(ereal-eferm)**4) ! Co(0001)
C       vil=atan(0.05+0.02*(ereal-eferm)**2+0.4*(ereal-eferm)**4) !  Fe(110) normal emission
C       vil=atan(0.15+0.02*(ereal-eferm)**2+0.004*(ereal-eferm)**4) !  Fe(110) winkelserie
C       vil=atan(0.15+0.02*(ereal-eferm)**2+0.1*(ereal-eferm)**4)
C       vil=atan(0.2+0.01*ereal**1.3) !W(100) JB
C
      VIL = VIL
C
      END
