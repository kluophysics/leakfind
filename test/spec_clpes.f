C*==spec_clpes.f    processed by SPAG 6.70Rc at 16:35 on 28 Feb 2017
      SUBROUTINE SPEC_CLPES(IBLOCH,NOUT1,IP,NREL,IREL,NSPIN,SPOL,NSPLEE,
     &                      EFEV,VPR,VIH,VIL,OMHAR,THQ,FIQ,ALQ,DELQ,
     &                      NPOL,ICIRC,IDREH,NATL,LAYS,LAYB,LAYP,LAYER1,
     &                      GANZ,LANZ1,POL0,POL0L,TYP,ISTR,Q1,Q2,Q3,Q4,
     &                      NEW,ESTEP,THETA,PHI,NT,NP,DELTAT,DELTAP,
     &                      TMINA,PMINA,IFM,POS,IRUMP,BRUSEP,NCLM,CELM,
     &                      CLIGHT,MAXL,NOD,NEV,NTOT,USEEULER,IBXY,ROT,
     &                      ROTM1,Z,ZPARU,ZPARD,EPSX,EBIND,PHOTOC,
     &                      ROATLA,AA,AWR,ZMAT0,AMAT1X,AMAT1Y,AMAT1V,
     &                      AMAT1,AMAT2X,AMAT2Y,AMAT2V,AMAT2,PSI2G,
     &                      BULKX,IGPRO,IGVAL,TSEL,TSOL,TSEOL,TSOEL,
     &                      TSEH,TSOH,TSEOH,TSOEH,IEDLAYER,IODLAYER,
     &                      IEODLAYER,IOEDLAYER,ZMESH,ADU,ADD,BK,IFSP,
     &                      VLM,MEZJ,MEZZ,MSSQ,SSST,TAUQ,TAUT,PHASK,
     &                      FIXTHQ,ITXRAY,LCXRAY,NCXRAY,ASIG,FROMSIGMA,
     &                      TOSIGMA,EMIN,EMAX,ASYM,IE_ICST,
     &                      POL0_VEC_TILT,POL0_INITIAL)
C
C   ********************************************************************
C   *                                                                  *
C   *  driver routien to calculate angle resolved CLPES and XPD        *
C   *                                                                  *
C   ********************************************************************
C
C
C
C
      USE MOD_CPA,ONLY:NCPA
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R2DRDI,R
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_TYPES,ONLY:NTMAX,IMT
      USE MOD_SITES,ONLY:NQMAX,NOQ
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_MPI,ONLY:MPI,MPI_ID,NPROCS
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,MQK,MLS,LL,XMAXE,NFULLPOT,
     &    MLQNAT,MEW,MPW,LH,MHDIM2,MLZP,HARTRE,PI,NTPHOMAX,NVFTPHOMAX,
     &    MAXCSTATES,MAXCORE,MAXEIG,MAXN,RSTEPP5
      USE MOD_SPEC_WAVE,ONLY:FF_S2,FF_S2M
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC_POTLM,ONLY:EB,BB,ALPHA
      USE MOD_FILES,ONLY:IFILSPECOU3,IFILSPECOU1,IPRINT,DATSET,LDATSET,
     &    IFILSPECIO3
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALQ,ASYM,CLIGHT,DELQ,DELTAP,DELTAT,EFEV,EMAX,EMIN,EPSX,
     &       ESTEP,FIQ,FIXTHQ,FROMSIGMA,OMHAR,PHI,PMINA,THETA,THQ,TMINA,
     &       TOSIGMA,VIH,VIL,VPR
      INTEGER ASIG,GANZ,IBLOCH,IBXY,ICIRC,IDREH,IE_ICST,IFM,IFSP,IP,
     &        IREL,ITXRAY,LANZ1,LAYB,LAYER1,LAYP,LAYS,MAXL,NCLM,NEW,
     &        NOUT1,NP,NPOL,NREL,NSPIN,NSPLEE,NT,SPOL,TYP,USEEULER
      LOGICAL POL0_VEC_TILT
      COMPLEX*16 Q1,Q2,Q3,Q4
      COMPLEX*16 AA(3),AMAT1(XMAXE,4,NFULLPOT,2),AMAT2(XMAXE,4,NFULLPOT)
     &           ,AWR(MLS,LL,NATLM,2),BK(LH,2),BULKX(LH,LH,LAYSM),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),PHASK(NEMAX),
     &           PSI2G(MHDIM2,LH,4),ROT(MQD,MQD),ROTM1(MQD,MQD),
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
     &       POS(3,NATLM,LAYSM),ROATLA(LL),Z(NATLM,LAYSM),ZMESH(6),
     &       ZPARD(3),ZPARU(3)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        IEDLAYER(0:LAYSM),IEODLAYER(0:LAYSM),IGPRO(LH),IGVAL(LH),
     &        IODLAYER(0:LAYSM),IOEDLAYER(0:LAYSM),IRUMP(LAYSM),ISTR(2),
     &        LCXRAY(NTMAX),NATL(LAYSM),NCXRAY(NTMAX),NEV(LAYSM),
     &        NOD(LAYSM),NTOT(LAYSM)
C
C Local variables
C
      REAL*8 ALLQ,BCOR(:),BCORS(:),CWFENERGY(:,:,:),CWFG(:,:,:),
     &       CWFKAPMUE(:,:,:,:),E1,ECOR(:),EMAXCORE,EMINCORE,
     &       EVAL(:,:,:,:),FCOR(:,:,:),GCOR(:,:,:),KPARAF(:,:),KPX,KPY,
     &       KXF(:,:),KYF(:,:),MUE,PHID,PHIQ,RD,RINT(:),ROAT,
     &       RWORK(:,:,:),RWORK1(:),SZCOR(:),TESTVAC,THEQ,THETAD,
     &       XNORM(2),XPHO(2)
      COMPLEX*16 AMAT1T(:,:,:,:),AMAT2T(:,:,:),C1(:,:),C2(:,:),
     &           CWF(:,:,:),DMATTG(:,:,:),DOS(:,:,:),DTILTG(:,:,:),
     &           EHIGH,ELOW,GAMMA1(:,:,:,:,:),K1,MSSTS(:,:,:,:),P(2),PF,
     &           RDIP1M(:,:,:,:,:),RDIP2M(:,:,:,:),ROATLASD(:),
     &           TAUT_GLO(:,:,:),TSSQ(:,:,:),TSSTS(:,:,:,:),
     &           TSST_GLO(:,:,:),ZMATD(:,:,:,:)
      INTEGER ATOM,CWFSTATES(:,:),CWFX(:,:,:,:,:),CWFY(:),G1(:),G2(:),I,
     &        I1,IA_ERR,ICST,IDRQ,IE,IED,IEOD,IERR,IFILE,IGP,IKMCOR(:,:)
     &        ,IM,INC,IO,IOD,IOED,IPOL,IPROC,IPROCA(:),IPROCE(:),IQ,
     &        ISTATE,IT,IZERO(:),J,JA,JE,JP,JT,K,KAPCOR(:),LAY,LOUT,
     &        MAXG,MM05COR(:),N,NCST,NCSTMAX,NKPCOR(:)
      LOGICAL CALCDOS,TILT
      CHARACTER*80 DATAFILE
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE FCOR,GCOR,ECOR,SZCOR,KAPCOR,RWORK,RWORK1
      ALLOCATABLE MM05COR,IKMCOR,NKPCOR,IZERO,KPARAF,KXF,KYF
      ALLOCATABLE BCOR,BCORS,RINT,GAMMA1
      ALLOCATABLE AMAT1T,AMAT2T,CWFSTATES
      ALLOCATABLE CWFENERGY,CWFKAPMUE
      ALLOCATABLE cwfg,g1,g2,zmatd
      ALLOCATABLE     cwfx,cwfy
      ALLOCATABLE     cwf,eval
      ALLOCATABLE RDIP1M,RDIP2M,TSSTS,MSSTS
      ALLOCATABLE tssq
      ALLOCATABLE DMATTG,DTILTG,c1,c2
      ALLOCATABLE TSST_GLO,TAUT_GLO,DOS,ROATLASD,iproca,iproce
C
C
      NCSTMAX = 14
      ALLOCATE (GAMMA1(LAYSM,NATLM,3,MQD,MQD))
      ALLOCATE (TSSTS(NKMMAX,NKMMAX,NTMAX,2))
      ALLOCATE (MSSTS(NKMMAX,NKMMAX,NTMAX,2))
      ALLOCATE (CWF(NRMAX,MAXCORE,2))
      ALLOCATE (CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),CWFY(MAXCORE))
      ALLOCATE (EVAL(MAXEIG,MAXN,NATLM,LAYSM))
      ALLOCATE (RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX))
      ALLOCATE (RDIP2M(NFULLPOT,NRMAX,2,NTPHOMAX))
      ALLOCATE (G1(MQK),G2(MQK))
      ALLOCATE (CWFG(MAXCSTATES,NATLM,LAYSM))
      ALLOCATE (CWFSTATES(NATLM,LAYSM))
      ALLOCATE (CWFENERGY(MAXCSTATES,NATLM,LAYSM))
      ALLOCATE (CWFKAPMUE(MAXCSTATES,NATLM,LAYSM,2))
      ALLOCATE (ECOR(NCSTMAX),SZCOR(NCSTMAX))
      ALLOCATE (KAPCOR(NCSTMAX),MM05COR(NCSTMAX),IZERO(NCSTMAX))
      ALLOCATE (IKMCOR(NCSTMAX,2),NKPCOR(NCSTMAX))
      ALLOCATE (BCOR(NTMAX),BCORS(NTMAX),RINT(NRMAX))
      ALLOCATE (AMAT1T(XMAXE,4,NFULLPOT,2),AMAT2T(XMAXE,4,NFULLPOT))
      ALLOCATE (ZMATD(MQD,MQD,LAYSM,NATLM))
      ALLOCATE (ROATLASD(LL))
      ALLOCATE (TSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (DMATTG(NKMMAX,NKMMAX,NTMAX),DTILTG(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TSST_GLO(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TAUT_GLO(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DOS(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (C1(MQD,MQD),C2(MQD,MQD))
      ALLOCATE (FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: SPEC_CLPES -> GCOR'
C
      IF ( .NOT.ALLOCATED(FF_S2M) ) THEN
         ALLOCATE (FF_S2M(NRMAX,MQD,MQD))
         ALLOCATE (FF_S2(RSTEPP5,MQD,MQD))
         FF_S2 = DCMPLX(0.0D0,0.0D0)
         FF_S2M = DCMPLX(0.0D0,0.0D0)
      END IF
C
      IF ( NCPA.GT.0 ) STOP '<SPEC_CLPES>: Not yet in CPA'
C
C
      IT = ITXRAY
      IM = IMT(IT)
      RD = PI/180.0D0
      PHOTOC(:,:,:,:,:) = 0.0D0
      TILT = .FALSE.
      THETAD = 0.0D0
      PHID = 0.0D0
      CALL SPEC_TILT(THETAD,PHID,TILT,0,NP,POL0,POL0_VEC_TILT,
     &               POL0_INITIAL,IP)
C
      DATAFILE = DATSET(1:LDATSET)//'_corewf'
      LOUT = LDATSET + 7
C
      IF ( MPI_ID.EQ.0 ) THEN
         WRITE (6,99005) DATAFILE
         OPEN (UNIT=IFILSPECIO3,FILE=DATAFILE(1:LOUT))
      END IF
      EMINCORE = 1D10
      EMAXCORE = -1D10
C
C
C     for loop over photon polarisation (linear or circular dichroism)
C     (always inner loop)
      IF ( NPOL.EQ.0 .OR. NPOL.EQ.3 )
     &      STOP '<SPEC_CLPES>: LOOP OVER NPHPON NOT YET'
CC
C=======================================================================
C                   calculate core levels and wave functions
C     NB: this part needs to be generalise in case of overlapping core
C     levels
C=======================================================================
C
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            WRITE (6,99001) ITXRAY,NCXRAY(ITXRAY),LCXRAY(ITXRAY)
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
            NCST = 4*LCXRAY(IT) + 2
C
            WRITE (6,99004)
C
            DO ICST = 1,NCST
C
               DO K = 1,NKPCOR(ICST)
                  DO N = 1,JRWS(IM)
                     RINT(N) = R2DRDI(N,IM)
     &                         *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
                  END DO
                  CALL RRADINT(IM,RINT,XNORM(K))
               END DO
               WRITE (6,99002) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                         (2*MM05COR(ICST)+1),IKMCOR(ICST,1),
     &                         XNORM(1),ECOR(ICST),ECOR(ICST)*RY_EV,
     &                         SZCOR(ICST),IZERO(ICST)
               IF ( NKPCOR(ICST).EQ.2 ) WRITE (6,99003) IKMCOR(ICST,2),
     &              XNORM(2)
C
               MUE = (2.0D0*FLOAT(MM05COR(ICST))+1.0D0)/2.0D0
C
               EMINCORE = MIN(EMINCORE,ECOR(ICST)*RY_EV)
               EMAXCORE = MAX(EMAXCORE,ECOR(ICST)*RY_EV)
C
C
               IF ( MPI_ID.EQ.0 )
     &              CALL SPEC_CORE_WRWF(FCOR,GCOR,JRWS(IM),R(1,IM),LAY,
     &              ATOM,KAPCOR(ICST),MUE,ECOR(ICST),NCSTMAX,
     &              IFILSPECIO3,ICST,NRMAX)
C
            END DO
         END DO
      END DO
C
      IF ( MPI_ID.EQ.0 ) CLOSE (IFILSPECIO3)
C
      IF ( MPI ) THEN
         CALL SLEEP(2)
         CALL DRV_MPI_BARRIER
      END IF
      OPEN (UNIT=IFILSPECIO3,FILE=DATAFILE(1:LOUT))
C
      CALCDOS = .FALSE.
C
C=======================================================================
C     read core wavefunctions and energies from file:
C=======================================================================
C
      CALL READCORE(CWFX,CWFY,CWF,EVAL,LAYS,NATL,CWFSTATES,CWFENERGY,
     &              CWFKAPMUE,ROT,ROTM1,EB,IFILSPECIO3)
C
C
      IFM = 0
      CALL POTFULLM2(LAYS,NATL,Z,IFM,VLM)
      ALLQ = ALQ
      IDRQ = IDREH
      THQ = THQ*PI/180.0D0
      FIQ = FIQ*PI/180.0D0
      CALL SPEC_TILT(THQ,FIQ,TILT,1,NP,POL0,POL0_VEC_TILT,POL0_INITIAL,
     &               IP)
      THQ = THQ*180.0D0/PI
      FIQ = FIQ*180.0D0/PI
      CALL PHOTONVEC(AA,THQ,FIQ,ALLQ,DELQ,OMHAR,ICIRC,IDRQ)
C
      CALL DETMAT(AA,BB,VLM,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,
     &            AMAT2V,AMAT2,AMAT1T,AMAT2T,NATL,LAYS)
C
C
C=======================================================================
C     Add 5 eV to energy window
C=======================================================================
C
      EMINCORE = EMINCORE - 5.0D0
      EMAXCORE = EMAXCORE + 5.0D0
      IF ( NEW.NE.1 ) THEN
         ESTEP = (EMAXCORE-EMINCORE)/(NEW-1.)
         DO I1 = 1,NEW
            EBIND(I1) = EMINCORE + (I1-1)*ESTEP
         END DO
      ELSE
         EBIND(1) = EFEV
         ICST = IE_ICST
         EBIND(1) = ECOR(ICST)*RY_EV
         WRITE (6,99012)
         WRITE (6,99002) ICST,NCXRAY(IT),LCXRAY(IT),KAPCOR(ICST),
     &                   (2*MM05COR(ICST)+1),IKMCOR(ICST,1),XNORM(1),
     &                   ECOR(ICST),ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                   IZERO(ICST)
C
         EMINCORE = EBIND(1)
         EMAXCORE = EBIND(1)
      END IF
C     =================================================================
      ASIG = 0
      FROMSIGMA = 1.0D0
      TOSIGMA = 0.5D0
      EMIN = EMINCORE
      EMAX = -EMAXCORE
      NVFTPHOMAX = NTPHOMAX
      ISTATE = 2
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C                                            initialize MPI - parameters
C
C
      ALLOCATE (IPROCE(MEW))
      ALLOCATE (IPROCA(NT*NP))
      ALLOCATE (RWORK(MPW,2,MEW))
      ALLOCATE (KPARAF(MPW,MEW))
      ALLOCATE (KXF(MPW,MPW),KYF(MPW,MPW))
C
      IPROCE(:) = 0
      IPROCA(:) = 0
      KXF = 0.0D0
      KYF = 0.0D0
      RWORK(:,:,:) = 0.0D0
      KPARAF(:,:) = 0.0D0
C
      IF ( NEW.LE.1 ) THEN
C
C     Parallel set for the angular loops
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
C
C        Parallel set of for the energy loop
C
      ELSE IF ( NPROCS.GT.1 ) THEN
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
C
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C
      LOOP_I1_NEW:DO I1 = 1,NEW
CMPIMPIMPI
         IF ( IPROCE(I1).EQ.MPI_ID ) THEN
CMPIMPIMPI
            E1 = EBIND(I1)
            ELOW = CMPLX(REAL(E1),VIL)/HARTRE
            EHIGH = CMPLX(REAL(E1)+OMHAR*HARTRE,VIH)/HARTRE
            PF = SQRT
     &           (2.0D0*EHIGH+2.0D0*EHIGH*2.0D0*EHIGH/(CLIGHT*CLIGHT))
            TESTVAC = REAL(EHIGH) - VPR
            WRITE (6,*) 'Kinetic energy of:',REAL(EHIGH)*HARTRE
            WRITE (6,*) 'Initial energy of: ',REAL(ELOW)*HARTRE
            IF ( TESTVAC.LT.0.0D0 ) THEN
               WRITE (6,*) 'ELECTRON CANNOT ESCAPE'
               STOP
            END IF
C
C
            CALL DETG(ELOW,LAYS,NATL,ASIG,FROMSIGMA,TOSIGMA,EMIN,EMAX,
     &                CWFSTATES,CWFENERGY,CWFKAPMUE,CWFG,ASYM)
C
            DO ISTATE = 1,2
               CALL SPEC_TAU_DRIVE(MEZJ,MEZZ,MSSQ,TSSQ,
     &                             MSSTS(1,1,1,ISTATE),SSST,TAUQ,TAUT,
     &                             TSSTS(1,1,1,ISTATE),PHASK,DMATTG,
     &                             DTILTG,2.0D0*EHIGH,ISTATE,JE,
     &                             TSST_GLO,TAUT_GLO,DOS,CALCDOS,NEW,
     &                             P(ISTATE))
            END DO
            ISTATE = 2
C
C
            DO LAY = 1,LAYS
               DO ATOM = 1,NATL(IT)
                  IQ = IQ_SPR_LAYAT(ATOM,LAY)
                  IF ( NCPA.EQ.0 ) THEN
C
                     IQ = IQ_SPR_LAYAT(ATOM,LAY)
                     C1(:,:) = TSSQ(:,:,IQ)
                     CALL SPEC_TRANSPHO(C1,C2,NKMMAX,1,PF)
C
C
                     DO I = 1,MQD
                        DO J = 1,MQD
                           GAMMA1(LAY,ATOM,ISTATE,I,J) = C2(I,J)
                        END DO
                     END DO
                  ELSE
                     STOP 'CPA NOT YET IN XPS'
C
                  END IF
C
                  DO IO = 1,NOQ(IQ)
C
                     CALL RDIPOLEM(LAY,ATOM,VLM,RDIP1M,RDIP2M,CLIGHT,
     &                             OMHAR,ALPHA(LAY,ATOM),IO)
C
                     DO ISTATE = 1,2
C
                        CALL FULLF(LAY,ATOM,IO,ISTATE,EHIGH,CLIGHT,
     &                             USEEULER,IBLOCH,JE,VLM,MAXG,G1,G2,
     &                             TSSTS(1,1,1,ISTATE),
     &                             MSSTS(1,1,1,ISTATE),P(ISTATE))
                     END DO
                     ISTATE = 2
C
                     CALL CKMULM(MLQNAT,TSEL,TSOL,TSEOL,TSOEL,TSEH,TSOH,
     &                           TSEOH,TSOEH,ISTATE,LAY,ATOM,GAMMA1)
C
                     CALL ZMATXPS(CWFSTATES,CWFG,CWFX,CWFY,CWF,RDIP1M,
     &                            OMHAR,ZMAT0,ZMATD,ALPHA(LAY,ATOM),
     &                            ATOM,LAY,AMAT1X,AMAT1Y,AMAT1V,AMAT1,
     &                            AMAT2X,AMAT2Y,AMAT2V,IBLOCH,DOS,IO,
     &                            AMAT1T,AMAT2T)
                  END DO
               END DO
            END DO
            JA = 0
            DO JP = 1,NP
               DO JT = 1,NT
                  JA = JA + 1
C
CMPIMPIMPI
                  IF ( IPROCA(JA).EQ.MPI_ID ) THEN
CMPIMPIMPI
C
                     THETAD = (TMINA+DELTAT*DBLE(JT-1))
                     PHID = (PMINA+DELTAP*DBLE(JP-1))
                     CALL SPEC_TILT(THETAD,PHID,TILT,JP,NP,POL0,
     &                              POL0_VEC_TILT,POL0_INITIAL,IP)
                     THETA = THETAD
                     PHI = PHID
                     IF ( IFSP.EQ.1 ) THEN
                        THEQ = THETA*180.0/PI - FIXTHQ
                        THEQ = THEQ*PI/180.0D0
                        PHIQ = PHIQ*PI/180.0D0
                        CALL SPEC_TILT(THEQ,PHIQ,TILT,JP,NP,POL0,
     &                                 POL0_VEC_TILT,POL0_INITIAL,IP)
                        THEQ = THEQ*180.0D0/PI
                        PHIQ = PHIQ*180.0D0/PI
                        CALL PHOTONVEC(AA,THEQ,PHIQ,ALLQ,DELQ,OMHAR,
     &                                 ICIRC,IDRQ)
                        CALL DETMAT(AA,BB,VLM,AMAT1X,AMAT1Y,AMAT1V,
     &                              AMAT1,AMAT2X,AMAT2Y,AMAT2V,AMAT2,
     &                              AMAT1T,AMAT2T,NATL,LAYS)
                     END IF
C
                     DO IPOL = 1,SPOL
C
                        CALL FINALSTATE(GANZ,LANZ1,NSPIN,NREL,EHIGH,PHI,
     &                                  THETA,JE,VPR,DREAL(EHIGH),
     &                                  NSPLEE,IBLOCH,IREL,CLIGHT,CELM,
     &                                  TYP,ISTR,POL0,POL0L,IPOL,NATL,
     &                                  POS,IRUMP,LAYS,NOD,NEV,NTOT,
     &                                  BRUSEP,BULKX,NCLM,LAYB,TSEH,
     &                                  TSOH,TSEOH,TSOEH,IBXY,ELOW,EPSX,
     &                                  BK,ADU,ADD,IGP,IGPRO,ZMESH,Q1,
     &                                  Q2,Q3,Q4,LAYER1,MAXL,K1,AWR,
     &                                  LAYP,IFM,PSI2G,IGVAL,IOD,IED,
     &                                  IOED,IEOD,VIH,ZPARU,ZPARD,
     &                                  IEDLAYER,IODLAYER,IEODLAYER,
     &                                  IOEDLAYER,JT,NT,SPOL)
C
                        CALL XPSPHOTO(CLIGHT,OMHAR,ELOW,AWR,ZMAT0,NATL,
     &                                LAYS,LAYER1,ROATLA,2,LAYB,IPOL,
     &                                SPOL,ROATLASD)
C
                        IF ( IPOL.EQ.1 .AND. I1.EQ.1 .AND. JA.EQ.1 )
     &                       WRITE (6,99010)
                        ROAT = 0.0
                        XPHO(IPOL) = 0.0
                        DO I = 1,LAYP
                           ROAT = ROAT + ROATLA(I)
                           IF ( IP.GT.0 ) WRITE (NOUT1,99009) ROATLA(I)
                        END DO
                        IF ( IP.GT.0 ) WRITE (NOUT1,99006) REAL(ELOW)
     &                       *27.2116 - EMINCORE,ROAT
                        XPHO(IPOL) = ROAT
                     END DO
C
                     KPX = 0.3622605*SIN(THETA)
     &                     *SQRT(2.0*HARTRE*(CDABS(EHIGH)-VPR))*COS(PHI)
                     KPY = 0.3622605*SIN(THETA)
     &                     *SQRT(2.0*HARTRE*(CDABS(EHIGH)-VPR))*SIN(PHI)
                     KXF(JT,JP) = KPX
                     KYF(JT,JP) = KPY
C
C
                     DO IPOL = 1,SPOL
                        PHOTOC(I1,JT,JP,IPOL,1) = XPHO(IPOL)
                     END DO
C
                     IF ( NT.EQ.1 .AND. NP.EQ.1 ) THEN
                        WRITE (6,99007) REAL(ELOW)*27.2116,
     &                                  (XPHO(1)+XPHO(2))/2.0,XPHO(1),
     &                                  XPHO(2),(XPHO(1)-XPHO(2))
     &                                  /(XPHO(1)+XPHO(2))
                        WRITE (IFILSPECOU3,99007) REAL(ELOW)*27.2116,
     &                         (XPHO(1)+XPHO(2))/2.0,XPHO(1),XPHO(2),
     &                         (XPHO(1)-XPHO(2))/(XPHO(1)+XPHO(2))
                     ELSE
                        WRITE (6,99008) REAL(ELOW)*27.2116,THETA/RD,
     &                                  PHI/RD,(XPHO(1)+XPHO(2))/2.0,
     &                                  XPHO(1),XPHO(2),
     &                                  (XPHO(1)-XPHO(2))
     &                                  /(XPHO(1)+XPHO(2))
                        WRITE (IFILSPECOU3,99008) REAL(ELOW)*27.2116,
     &                         THETA/RD,PHI/RD,(XPHO(1)+XPHO(2))/2.0,
     &                         XPHO(1),XPHO(2),(XPHO(1)-XPHO(2))
     &                         /(XPHO(1)+XPHO(2))
                     END IF
C
                  END IF
               END DO
            END DO
C
         END IF
      END DO LOOP_I1_NEW
      CLOSE (IFILSPECIO3)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
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
            RWORK = 0.0D0
            INC = MEW*MPW*1
            CALL DRV_MPI_REDUCE_R(KPARAF,RWORK(1,1,1),INC)
            INC = NEW*NT*NP*2*2
            ALLOCATE (RWORK1(INC))
            RWORK1 = 0.0D0
            CALL DRV_MPI_REDUCE_R(PHOTOC,RWORK1(1),INC)
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
         END IF
C
      END IF
      IF ( MPI_ID.EQ.0 ) THEN
         IFILE = IFILSPECOU1
         DO JE = 1,NEW
            JA = 0
            DO JP = 1,NP
               DO JT = 1,NT
                  JA = JA + 1
                  THETAD = (TMINA+DELTAT*DBLE(JT-1))
                  PHID = (PMINA+DELTAP*DBLE(JP-1))
                  CALL SPEC_TILT(THETAD,PHID,TILT,JP,NP,POL0,
     &                           POL0_VEC_TILT,POL0_INITIAL,IP)
                  THETA = THETAD*180.D0/PI
                  PHI = PHID*180.D0/PI
                  IF ( NEW.GT.1 ) THEN
C
                     WRITE (IFILE,99011) THETA,EBIND(JE),
     &                      (PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1))
     &                      /2.D0,0.5*PHOTOC(JE,JT,JP,1,1),
     &                      0.5*PHOTOC(JE,JT,JP,2,1),
     &                      (PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1))
     &                      /(PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1))
     &                      ,KPARAF(JT,JE)
                  ELSE
                     WRITE (IFILE,99011) KXF(JT,JP),KYF(JT,JP),
     &                      (PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1))
     &                      /2.D0,0.5*PHOTOC(JE,JT,JP,1,1),
     &                      0.5*PHOTOC(JE,JT,JP,2,1),
     &                      (PHOTOC(JE,JT,JP,1,1)-PHOTOC(JE,JT,JP,2,1))
     &                      /(PHOTOC(JE,JT,JP,1,1)+PHOTOC(JE,JT,JP,2,1))
     &                      ,THETA,PHI
                  END IF
               END DO
            END DO
         END DO
         CALL FLUSH(IFILE)
      END IF
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
      IF ( MPI ) CALL MPI_FINALIZE(IERR)
C
      STOP 'CL_PES DONE'
C
99001 FORMAT (/,4X,' CORE QUANTUM-NUMBERS  FOR  IT=',I2,':   N=',I2,
     &        '  L=',I2,/)
99002 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5,/)
99003 FORMAT (22X,I4,F12.6)
99004 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
C
99005 FORMAT (/,10X,'*************************',/,10X,
     &        'Core WFs writen in file: ',/,10X,
     &        '*************************',/,10X,' ',/,10X,A,/)
99006 FORMAT (6(1x,e14.7))
99007 FORMAT (5(1x,e14.7))
99008 FORMAT (7(1x,e14.7))
99009 FORMAT (4(1x,e14.7))
99010 FORMAT (1x,'    INT    ','     UP    ','     DN    ','     P    ')
99011 FORMAT (2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,
     &        e14.5,2x,e14.5,2x,e14.5)
99012 FORMAT (/,4X,'Constant energy mode: ',/,4X,
     &        'Energy fixed to core state:')
C
      END
C*==spec_core_wrwf.f    processed by SPAG 6.70Rc at 16:35 on 28 Feb 2017
C
      SUBROUTINE SPEC_CORE_WRWF(FCOR,GCOR,JRWS,R,LAY,ATOM,KAPCOR,MUE,
     &                          ECOR,NCSTMAX,IFILSPECIO3,ICST,NRMAX)
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,ICST,IFILSPECIO3,JRWS,KAPCOR,LAY,NCSTMAX,NRMAX
      REAL*8 ECOR,MUE
      REAL*8 FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),R(NRMAX)
C
C Local variables
C
      INTEGER IR,J
C
C*** End of declarations rewritten by SPAG
C
C      write (IFILSPECIO3,1100) kapcor,mue,nkpcor,atom,lay,ECOR*RY_EV
      WRITE (IFILSPECIO3,99002) KAPCOR,MUE,1,ATOM,LAY,ECOR*RY_EV
      DO IR = 1,JRWS
         WRITE (IFILSPECIO3,99001) R(IR),
     &                             (DCMPLX(GCOR(IR,J,ICST),0.0D0),J=1,1)
     &                             ,(DCMPLX(0.0D0,0.0D0),J=1,2),
     &                             (DCMPLX(GCOR(IR,J,ICST),0.0D0),J=2,2)
C
      END DO
C      write (IFILSPECIO3,1200) kapcor,mue,nkpcor,atom,lay,ECOR*RY_EV
      WRITE (IFILSPECIO3,99003) KAPCOR,MUE,1,ATOM,LAY,ECOR*RY_EV
      DO IR = 1,JRWS
         WRITE (IFILSPECIO3,99001) R(IR),
     &                             (DCMPLX(FCOR(IR,J,ICST),0.0D0),J=1,1)
     &                             ,(DCMPLX(0.0D0,0.0D0),J=1,2),
     &                             (DCMPLX(FCOR(IR,J,ICST),0.0D0),J=2,2)
C
      END DO
C
99001 FORMAT (9(1x,e15.8))
99002 FORMAT ('#f: kappa=',i3,1x,'mu=',f5.1,1x,'n=',i3,1x,'atom=',i3,1x,
     &        'layer=',i3,1x,'energy=',f15.8)
99003 FORMAT ('#g: kappa=',i3,1x,'mu=',f5.1,1x,'n=',i3,1x,'atom=',i3,1x,
     &        'layer=',i3,1x,'energy=',f15.8)
      END
C
