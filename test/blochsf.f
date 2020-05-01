C*==blochsf.f    processed by SPAG 6.70Rc at 09:36 on 30 Mar 2017
      SUBROUTINE BLOCHSF
C   ********************************************************************
C   *                                                                  *
C   *   Bloch spectral function  A(k,E)                                *
C   *                                                                  *
C   *          NREL  SREL  SP-SREL   SPR                               *
C   * IREL      0     1      2       3                                 *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   * ASA       -     -      -       x                                 *
C   *                                                                  *
C   * FP        x     x      x       x                                 *
C   *                                                                  *
C   *                                                                  *
C   *   NBSFOP: 1: total bsf                                           *
C   *           2: Sigma_x                                             *
C   *           3: Sigma_y                                             *
C   *           4: Sigma_z                                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:KTAB,NKTABMAX
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT
      USE MOD_CALCMODE,ONLY:KMROT,ORBPOL,DMFT,LDAU,KKRMODE,IREL,
     &    GF_CONV_RH
      USE MOD_ENERGY,ONLY:EFERMI,NEMAX,NETAB,ETAB,EMAX,EMIN
      USE MOD_LATTICE,ONLY:BBAS,BRAVAIS,SUB_SYSTEM,SYSTEM_TYPE,
     &    SYSTEM_DIMENSION
      USE MOD_ANGMOM,ONLY:NKKR,IND0Q,NKM,NKMMAX,NKMQ,NLMMAX,NLMQ,MEZJ,
     &    MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,NQMAX,DROTQ,IQAT,NOQ,ITOQ,IQBOT,IQTOP
      USE MOD_TYPES,ONLY:NT,NCPLWFMAX,NTMAX,CTL,CONC,ITBOT,ITTOP
      USE MOD_FILES,ONLY:DATSET,IFILCBWF,IPRINT,LDATSET,LRECREAL8,
     &    LSYSTEM,SYSTEM,IFILBUILDBOT,WRBUILDBOT,FOUND_INTEGER,
     &    FOUND_REAL_ARRAY
      USE MOD_CPA,ONLY:NCPA
      USE MOD_CONSTANTS,ONLY:C0,C1,PI,RY_EV
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='BLOCHSF')
      INTEGER NKDIRMAX,NBSFOP
      PARAMETER (NKDIRMAX=20,NBSFOP=4)
C
C Local variables
C
      REAL*8 AKE(:,:,:,:,:,:),AMAX,AMIN,APLOT(:,:,:),ASUM,D11,D21,D22,
     &       DEL(3),EPHOT,EPHOTEV,EPLOT(:),KA(:,:),KE(:,:),KP(:),KSHIFT,
     &       KVEC(3),RWKT(:),SCL,SCL1,SCL2,TRC,VECK1(3),VECK2(3),
     &       VECK3(3),VECKA(3),WGT
      COMPLEX*16 CWGT,DMAMC(:,:),DMATT(:,:,:),DTILT(:,:,:),
     &           DZJAB(:,:,:,:),DZZ(1,1,1,1),DZZAB(:,:,:,:,:),ERYD,
     &           ESHIFT,FC(:,:),FCCQL(:,:,:,:),FFTAUQL(:,:,:,:),
     &           FJQL(:,:,:,:),MAUX(:,:),MQS(:,:,:,:),MSSQL(:,:),P,RAT,
     &           RSUM,RSUMD,TAUK(:,:),TAUKQQL(:,:),TAUKSP(:,:,:),
     &           TAUQL(:,:,:),TMAT(:,:),W1(:,:),W2(:,:)
      DOUBLE PRECISION DDOT,DNRM2
      CHARACTER*80 FILNAM,YTXT1,YTXT2
      LOGICAL GETIRRSOL
      INTEGER I,I1,IA_ERR,IBSFOP,ID,IE,IELOOP,IERR,IFIL,IFILBSF,IFILSP,
     &        IK,IK1,IK2,IKD,INC,INDKDIR(:),INFO,IO,IOA,IOB,IPIV(:),
     &        IPROC,IPROCE(:),IQ,IQBOTSAV,IQTOPSAV,IS,IT,ITA,ITB,
     &        ITBOTSAV,ITTOPSAV,IWRIRRWF,IWRREGWF,J,J1,KPATH,LFILNAM,
     &        LYTXT1,LYTXT2,M,N,NE,NK1,NK2,NKDIR,NKTAB0,NKTABD(:),
     &        NLM_LOC,NXYZ,RECLNGWF
      CHARACTER*8 LBLKDIR(:)
      CHARACTER*20 LEG(5)
      CHARACTER*10 MODE
      LOGICAL RVEC_SAME
      CHARACTER*4 STR4
C
C*** End of declarations rewritten by SPAG
C
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE FJQL,TAUK,MAUX,DMAMC,FCCQL,EPLOT,APLOT
      ALLOCATABLE DZZAB,DTILT,TAUQL,MSSQL,W1,W2,FC,KP,AKE
      ALLOCATABLE FFTAUQL,TAUKQQL,DMATT,DZJAB,KA,KE,INDKDIR
      ALLOCATABLE NKTABD,LBLKDIR,TMAT,IPIV,IPROCE,RWKT
      ALLOCATABLE MQS,TAUKSP
C
      ALLOCATE (TAUK(NKKR,NKKR),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      ALLOCATE (DZZAB(NKMMAX,NKMMAX,NTMAX,NTMAX,NBSFOP))
      ALLOCATE (DZJAB(NKMMAX,NKMMAX,NTMAX,NBSFOP),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DZJAB')
C
C
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (TAUKQQL(NKMMAX,NKMMAX))
      ALLOCATE (W2(NKMMAX,NKMMAX),FC(NKMMAX,NKMMAX))
      ALLOCATE (MSSQL(NKMMAX,NKMMAX),W1(NKMMAX,NKMMAX))
C
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NTMAX),STAT=IA_ERR)
C
      ALLOCATE (KA(3,NKDIRMAX),KE(3,NKDIRMAX))
      ALLOCATE (INDKDIR(NKDIRMAX),NKTABD(NKDIRMAX))
      ALLOCATE (LBLKDIR(NKDIRMAX))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FT')
C
      ALLOCATE (FFTAUQL(NKMMAX,NKMMAX,NQMAX,NBSFOP))
      ALLOCATE (FCCQL(NKMMAX,NKMMAX,NQMAX,NBSFOP))
      ALLOCATE (FJQL(NKMMAX,NKMMAX,NQMAX,NBSFOP))
      ALLOCATE (TMAT(NKMMAX,NKMMAX),IPIV(NKMMAX))
      ALLOCATE (TAUQL(NKMMAX,NKMMAX,NQMAX),KP(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FFTAUQL')
      IF ( IREL.EQ.2 ) THEN
         ALLOCATE (MQS(NLMMAX,NLMMAX,NQ,2))
         ALLOCATE (TAUKSP(NKKR,NKKR,2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUKSP')
      END IF
      ALLOCATE (IPROCE(NEMAX))
C
      IF ( IREL.LT.3 .AND. .NOT.FULLPOT )
     &      CALL STOP_MESSAGE(ROUTINE,'SP-SREL not available in ASA')
C
      TAUQ = C0
      TMAT = C0
      TAUQL = C0
      FJQL = C0
      FCCQL = C0
      FFTAUQL = C0
C
      GETIRRSOL = .TRUE.
C
      WRITE (6,99016)
C
      IFILBSF = 11
      IFILSP = 12
C
      FILNAM = DATSET(1:LDATSET)//'_spol.bsf'
      CALL WRHEAD(IFILSP,FILNAM,'BSF-SPOL  ',NETAB(1))
C
      NKDIR = 0
      DO ID = 1,NKDIRMAX
         LBLKDIR(ID) = '        '
      END DO
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      IF ( NCPA.NE.0 ) THEN
         CALL READTAU(9,ERYD,0,NE,W1,0,NT,.FALSE.,TAUQ,MSSQ,0,NQ,0,
     &                NKMMAX,NKMMAX,IPRINT)
      ELSE
         NE = NETAB(1)
      END IF
C
      IF ( NE.GT.NEMAX ) CALL STOP_MESSAGE(ROUTINE,'NE > NEMAX ')
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      EPHOTEV = 0D0
      CALL SECTION_SET_REAL('EPHOT',EPHOTEV,9999D0,0)
      EPHOTEV = MAX(0D0,EPHOTEV)
      EPHOT = EPHOTEV/RY_EV
C
C-----------------------------------------------------------------------
C                                                               E-k-plot
      IF ( NE.GT.1 ) THEN
C
         MODE = 'E-k-plot'
         NXYZ = 3
C
         WRITE (6,99010)
         NK1 = 1
         CALL SECTION_SET_INTEGER('NK',NK2,800,0)
         CALL SECTION_SET_INTEGER('NKDIR',NKDIR,9999,0)
         IF ( .NOT.FOUND_INTEGER .AND. NK2.LE.2 ) NKDIR = 1
         NKDIR = MIN(NKDIR,9)
C
         IF ( NKDIR.NE.0 ) THEN
C
            DO ID = 1,NKDIR
               STR4 = 'KA'//CHAR(ICHAR('1')-1+ID)
               CALL SECTION_SET_REAL_ARRAY(STR4,KA(1,ID),NXYZ,3,1,
     &            9999D0,0)
               STR4 = 'KE'//CHAR(ICHAR('1')-1+ID)
               CALL SECTION_SET_REAL_ARRAY(STR4,KE(1,ID),NXYZ,3,1,
     &            9999D0,0)
               IF ( .NOT.FOUND_REAL_ARRAY .AND. NKDIR.EQ.1 .AND. 
     &              NK2.EQ.1 ) KE(:,1) = KA(:,1)
               WRITE (6,99017) ID,(KA(I,ID),I=1,3),(KE(I,ID),I=1,3)
               LBLKDIR(ID) = '        '
            END DO
C
            WRITE (6,*) ' '
C
         ELSE
C
            CALL SECTION_SET_INTEGER('KPATH',KPATH,9999,1)
C
            CALL KDIRTAB(BRAVAIS,KPATH,NKDIR,LBLKDIR,KA,KE,BBAS)
C
         END IF
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
            DO ID = 1,NKDIR
               KA(3,ID) = 0.0D0
               KE(3,ID) = 0.0D0
            END DO
            WRITE (6,*) 
     &        'WARNING from <BLOCHSF>: KA_Z and KE_Z set to zero for 2D'
         END IF
C
         IF ( NK2.GT.NKTABMAX ) THEN
            NK2 = NKTABMAX
            WRITE (6,*) 
     &               'WARNING from <BLOCHSF>:  NK reduced to array size'
     &               ,NKTABMAX
         END IF
C
C-----------------------------------------------------------------------
         IF ( NKDIR.EQ.0 ) THEN
            NKDIR = 1
            DO I = 1,3
               KA(I,1) = 0.0D0
               KE(I,1) = 0.0D0
            END DO
            KE(1,1) = 1.0D0
            NKTABD(1) = NK2
            LBLKDIR(1) = 'GG-x -1 '
         ELSE
            RSUMD = 0.0D0
            DO ID = 1,NKDIR
               DO I = 1,3
                  DEL(I) = KE(I,ID) - KA(I,ID)
               END DO
               RSUM = DNRM2(3,DEL,1)
               NKTABD(ID) = INT(1000000*RSUM)
               RSUMD = RSUMD + RSUM
            END DO
            NKTAB0 = 0
            DO ID = 1,NKDIR
               NKTABD(ID) = INT(NK2*NKTABD(ID)/(RSUMD*1000000))
               NKTABD(ID) = MAX(2,NKTABD(ID))
               NKTAB0 = NKTAB0 + NKTABD(ID)
            END DO
            IF ( NKTAB0.GT.NKTABMAX ) THEN
               RAT = 0.95D0*DBLE(NKTABMAX)/DBLE(NKTAB0)
               NKTAB0 = 0
               DO ID = 1,NKDIR
                  NKTABD(ID) = MAX(2,INT(RAT*NKTABD(ID)))
                  NKTAB0 = NKTAB0 + NKTABD(ID)
               END DO
            END IF
C
            NK2 = NKTAB0
            IF ( NKDIR.EQ.1 ) THEN
               IF ( RVEC_SAME(3,KA(1,1),KE(1,1),1D-6) ) NK2 = 1
            END IF
C
         END IF
C
         IF ( NK1.EQ.1 .AND. NK2.LE.2 ) MODE = 'E-plots   '
C
         INDKDIR(1) = NKTABD(1)
         DO ID = 2,NKDIR
            INDKDIR(ID) = INDKDIR(ID-1) + NKTABD(ID)
         END DO
C
C-----------------------------------------------------------------------
C----------------------------------- set up k-points and linear array KP
C
         IK = 0
         KP(1) = 0D0
         DO ID = 1,NKDIR
C
            DO I = 1,3
               DEL(I) = (KE(I,ID)-KA(I,ID))/DBLE(NKTABD(ID)-1)
            END DO
C
            DO IKD = 1,NKTABD(ID)
               IK = IK + 1
               DO I = 1,3
                  KTAB(I,IK) = KA(I,ID) + DEL(I)*DBLE(IKD-1)
               END DO
               IF ( IKD.NE.1 ) THEN
                  KP(IK) = KP(IK-1) + DNRM2(3,DEL,1)
               ELSE IF ( IK.EQ.1 ) THEN
                  KP(IK) = 0.0D0
               ELSE
                  KP(IK) = KP(IK-1)
               END IF
            END DO
         END DO
C
         IF ( IPRINT.GE.1 ) THEN
            IK = 0
            DO ID = 1,NKDIR
               DO IKD = 1,NKTABD(ID)
                  IK = IK + 1
                  WRITE (6,'(A,4I5,3F8.3,F10.3)') ' ->K ',ID,NKTABD(ID),
     &                   IKD,IK,(KTAB(I,IK),I=1,3),KP(IK)
               END DO
            END DO
         END IF
C
         WRITE (6,99011) EFERMI,EMIN,EMAX,NE,NK2
         IF ( EPHOT.GT.1D-3 ) WRITE (6,99001) 'E_ph      ',EPHOTEV,
     &                               ' eV    ',EPHOT,' Ry'
         DO IFIL = IFILBSF,IFILSP
            WRITE (IFIL,99008)
            WRITE (IFIL,99009) DATSET(1:LDATSET),NE,NK2
            WRITE (IFIL,99008)
            WRITE (IFIL,99007) 'EMIN EMAX ',EMIN,EMAX
            WRITE (IFIL,99008)
            WRITE (IFIL,99004) 'NKDIR     ',NKDIR
            DO I = 1,NKDIR
               WRITE (IFIL,99005) LBLKDIR(I)
            END DO
            WRITE (IFIL,99008)
            DO I = 1,NKDIR
               WRITE (IFIL,99004) 'INDKDIR   ',INDKDIR(I)
            END DO
            WRITE (IFIL,99008)
            WRITE (IFIL,99004) 'NK         ',NK2
            DO IK = 1,NK2
               WRITE (IFIL,99006) KP(IK)
            END DO
            WRITE (IFIL,99008)
         END DO
C
C-----------------------------------------------------------------------
C                                                               k-k-plot
      ELSE
C
         MODE = 'k-k-plot'
         NXYZ = 3
C
         WRITE (6,99012)
C
         CALL SECTION_SET_INTEGER('NK1',NK1,9999,1)
         CALL SECTION_SET_REAL_ARRAY('K1',VECK1,NXYZ,3,0,9999D0,1)
         CALL SECTION_SET_INTEGER('NK2',NK2,9999,1)
         CALL SECTION_SET_REAL_ARRAY('K2',VECK2,NXYZ,3,0,9999D0,1)
         CALL SECTION_SET_REAL_ARRAY('KA',VECKA,NXYZ,3,0,9999D0,0)
C
         IF ( .NOT.FOUND_REAL_ARRAY ) THEN
            WRITE (6,*) 'WARNING: in <BLOCHSF>   ->KA not set'
            WRITE (6,*) 'K1 and K2 are supposed to start at 0,0,0'
            VECKA(1:3) = 0.0D0
         END IF
C
         IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
            VECK1(3) = 0.0D0
            VECK2(3) = 0.0D0
            WRITE (6,*) 
     &        'WARNING from <BLOCHSF>: K1_Z and K2_Z set to zero for 2D'
         END IF
C
         VECK1 = VECK1 - VECKA
         VECK2 = VECK2 - VECKA
         D11 = DDOT(3,VECK1,1,VECK1,1)
         D22 = DDOT(3,VECK2,1,VECK2,1)
C
         WRITE (6,99013) EFERMI,ETAB(1,1),NK1,VECK1 + VECKA,SQRT(D11),
     &                   NK2,VECK2 + VECKA,SQRT(D22)
         WRITE (6,*) 'VECKA',VECKA
         IF ( EPHOT.GT.1D-3 ) WRITE (6,99001) 'E_ph      ',EPHOTEV,
     &                               ' eV    ',EPHOT,' Ry'
         D21 = DDOT(3,VECK2,1,VECK1,1)
         IF ( ABS(D21).GT.1D-6 ) THEN
            SCL = -D21/D11
            CALL DAXPY(3,SCL,VECK1,1,VECK2,1)
            VECK3 = VECK2 + VECKA
            WRITE (6,99014) NK2,VECK2 + VECKA,
     &                      SQRT(DDOT(3,VECK3,1,VECK3,1))
         END IF
         DO IFIL = IFILBSF,IFILSP
            WRITE (IFIL,99008)
            WRITE (IFIL,99015) DATSET(1:LDATSET),NK1,NK2,ETAB(1,1)
            WRITE (IFIL,99004) 'NK1        ',NK1
            WRITE (IFIL,99007) 'VECK1      ',VECK1 + VECKA
            WRITE (IFIL,99004) 'NK2        ',NK2
            WRITE (IFIL,99007) 'VECK2      ',VECK2 + VECKA
            WRITE (IFIL,99008)
         END DO
      END IF
C-----------------------------------------------------------------------
C
      CLOSE (IFILCBWF)
C
      IF ( FULLPOT ) THEN
         RECLNGWF = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         RECLNGWF = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      KSHIFT = 1D-4
C
      ALLOCATE (AKE(NK1,NK2,NE,NQ,3,NBSFOP),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: AKE')
C
      AKE(:,:,:,:,:,:) = 0.0D0
C
C=======================================================================
C               loop over the energy windows separated by   EPHOT
C=======================================================================
      IELOOP = 1
      ESHIFT = 0D0
 100  CONTINUE
      IF ( EPHOT.GT.1D-3 ) WRITE (6,99002) IELOOP
C
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      IF ( DMFT .OR. LDAU ) CALL DMFT_READSIG(NETAB(1),ETAB(1,1),IPRINT)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IPROCE(1:NEMAX) = 0
C
      IF ( NPROCS.GT.1 ) THEN
         IPROC = 0
         INC = 1
         DO IE = 1,NE - 1
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
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      DO IE = 1,NE
C
C======= read site-diagonal TAU- and inverse eff. single site t-matrix m
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCE(IE) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            IF ( NCPA.NE.0 ) THEN
               DO IQ = 1,NQ
                  IF ( MPI ) THEN
                     ERYD = ETAB(IE,1)
                     CALL READTAU(9,ERYD,IE,NE,W1,0,NT,.FALSE.,
     &                            TAUQ(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,2,
     &                            NKMQ(IQ),NKMMAX,IPRINT)
                  ELSE
                     CALL READTAU(9,ERYD,IE,NE,W1,0,NT,.FALSE.,
     &                            TAUQ(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,1,
     &                            NKMQ(IQ),NKMMAX,IPRINT)
                  END IF
               END DO
C
               ETAB(IE,1) = ERYD
            ELSE
               ERYD = ETAB(IE,1) + ESHIFT
            END IF
C
C
            IF ( DIMAG(ERYD).GT.1D-8 ) GETIRRSOL = .TRUE.
C
C ------------------------------ calculate WF and DOS - radial integrals
C
            IF ( DMFT .OR. LDAU ) DMFTSIG(1:NKM,1:NKM,1:NT)
     &           = DMFTSIGMA(1:NKM,1:NKM,1:NT,IE)
C
            IWRREGWF = 1
            IWRIRRWF = 1
C
C
            IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6)
     &           .EQ.'I-ZONE' ) THEN
               ITBOTSAV = ITBOT
               ITTOPSAV = ITTOP
               IQBOTSAV = IQBOT
               IQTOPSAV = IQTOP
               ITBOT = 1
               ITTOP = NT
               IQBOT = 1
               IQTOP = NQ
            END IF
C
            CALL RUNSSITE(.TRUE.,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,
     &                    ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            IF ( IREL.EQ.3 ) THEN
               CALL ME_OBS_REL(IFILCBWF,IFILCBWF,.FALSE.,DZZ,DZJAB,
     &                         DZZAB,NBSFOP)
            ELSE
               CALL BLOCHINT(MEZZ,MEZJ,DZZAB,DZJAB,NBSFOP)
            END IF
C
C --------------------------- initialize MSSQ in case of ordered systems
C
            IF ( NCPA.EQ.0 ) CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
            IF ( SYSTEM_TYPE(1:3).NE.'VIV' .AND. SUB_SYSTEM(1:6)
     &           .EQ.'I-ZONE' ) THEN
               ITBOT = ITBOTSAV
               ITTOP = ITTOPSAV
               IQBOT = IQBOTSAV
               IQTOP = IQTOPSAV
            END IF
C
C
C--------------------------------------- calculate relativistic momentum
C
            P = SQRT(ERYD*(1.0D0+ERYD/CTL(1,1)**2))
C
C-------- set up the matrices TAUQL and MSSQL referred to the local axis
C----- the resulting projection matrices DMATT and DTILT, as well as the
C---------------------------- matrices TAUT also refer to the local axis
C-------------------------  TAUQ(IQ) for equivalent sites are the same !
C
            M = NKMMAX
            N = NKM
C
            DO IT = ITBOT,ITTOP
               IQ = IQAT(1,IT)
C
               IF ( NKMQ(IQ).NE.N ) THEN
                  WRITE (6,*) 'NKM     = N: ',N
                  WRITE (6,*) 'NKMQ(IQ):    ',NKMQ(IQ)
                  WRITE (6,*) 'IQ  IT:      ',IT,IQ
                  CALL STOP_MESSAGE(ROUTINE,
     &                         'the same NKMQ required for all sites IQ'
     &                         )
               END IF
C
               IF ( KMROT.NE.0 ) THEN
C
                  CALL ROTATE(MSSQ(1,1,IQ),'G->L',MSSQL,N,DROTQ(1,1,IQ),
     &                        M)
C
                  CALL ROTATE(TAUQ(1,1,IQ),'G->L',TAUQL(1,1,IQ),N,
     &                        DROTQ(1,1,IQ),M)
C
               ELSE
C
                  DO J = 1,N
                     CALL ZCOPY(N,MSSQ(1,J,IQ),1,MSSQL(1,J),1)
                     CALL ZCOPY(N,TAUQ(1,J,IQ),1,TAUQL(1,J,IQ),1)
                  END DO
C
               END IF
C
               CALL TAUGFCONV(MSSQL(1,1),TAUQL(1,1,IQ),TMAT)
               TAUQL(1:NKMMAX,1:NKMMAX,IQ) = TMAT(1:NKMMAX,1:NKMMAX)
C
               IF ( GF_CONV_RH ) THEN
C
C                  ii     i       i    (-1)
C      D  = [ 1 + G   * (t    -  t  ) ]
C       a          CPA    CPA     a
C
C      _            i       i      ii   (-1)
C      D  = [ 1 + (t    -  t  ) * G    ]
C       a           CPA     a      CPA
C
                  TMAT(:,:) = MSSQL(:,:)
                  CALL ZGETRF(N,N,TMAT,NKMMAX,IPIV,INFO)
                  CALL ZGETRI(N,TMAT,NKMMAX,IPIV,MAUX,NKMMAX*NKMMAX,
     &                        INFO)
C
                  CALL GETDMAT(TAUQL(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT)
     &                         ,DMAMC,N,TSST(1,1,IT),TMAT,M)
               ELSE
                  CALL GETDMAT(TAUQL(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT)
     &                         ,DMAMC,N,MSSQL,MSST(1,1,IT),M)
               END IF
C
               CALL ZGEMM('N','N',N,N,N,C1,DMATT(1,1,IT),M,TAUQL(1,1,IQ)
     &                    ,M,C0,TAUT(1,1,IT),M)
            END DO
C
C=======================================================================
C
C
C
C
C=======================================================================
C
            CALL CINIT(NKMMAX*NKMMAX*NQMAX*NBSFOP,FJQL)
            CALL CINIT(NKMMAX*NKMMAX*NQMAX*NBSFOP,FCCQL)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
            DO IQ = IQBOT,IQTOP
C
               DO IBSFOP = 1,NBSFOP
C
                  N = NKMQ(IQ)
                  M = NKMMAX
C
                  IF ( GETIRRSOL ) THEN
                     DO IO = 1,NOQ(IQ)
                        IT = ITOQ(IO,IQ)
                        CWGT = -CONC(IT)/PI
C
C
                        DO I = 1,N
                           FJQL(I,I,IQ,IBSFOP) = FJQL(I,I,IQ,IBSFOP)
     &                        + CWGT*DZJAB(I,I,IT,IBSFOP)
                        END DO
C
                     END DO
                  END IF
C
C=======================================================================
C
                  CALL CINIT(NKMMAX*NKMMAX,FC)
C
                  DO IO = 1,NOQ(IQ)
                     IT = ITOQ(IO,IQ)
                     CWGT = -CONC(IT)/PI
C
                     CALL ZGEMM('N','N',N,N,N,C1,DZZAB(1,1,IT,IT,IBSFOP)
     &                          ,M,DMATT(1,1,IT),M,C0,W2,M)
C
                     DO J = 1,N
                        CALL ZAXPY(N,CWGT,W2(1,J),1,FC(1,J),1)
                     END DO
                  END DO
C
C
C=======================================================================
C
                  DO IOA = 1,NOQ(IQ)
                     ITA = ITOQ(IOA,IQ)
                     DO IOB = 1,NOQ(IQ)
                        ITB = ITOQ(IOB,IQ)
C
                        CALL ZGEMM('N','N',N,N,N,C1,
     &                             DZZAB(1,1,ITA,ITB,IBSFOP),M,
     &                             DMATT(1,1,ITB),M,C0,W1,M)
C
                        CALL ZGEMM('N','N',N,N,N,C1,DTILT(1,1,ITA),M,W1,
     &                             M,C0,W2,M)
C
                        CWGT = -CONC(ITA)*CONC(ITB)/PI
C
                        DO J = 1,N
C
                           CALL ZAXPY(N,CWGT,W2(1,J),1,
     &                                FCCQL(1,J,IQ,IBSFOP),1)
C
                        END DO
                     END DO
                  END DO
C
C
C=======================================================================
C
                  DO J = 1,N
                     DO I = 1,N
                        W1(I,J) = FC(I,J) - FCCQL(I,J,IQ,IBSFOP)
                     END DO
                  END DO
C
                  CALL ZGEMM('N','N',N,N,N,C1,W1,M,TAUQL(1,1,IQ),M,C0,
     &                       FFTAUQL(1,1,IQ,IBSFOP),M)
C
               END DO
            END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
            WRITE (6,99003) IE,ERYD
C
C -------------- calculate energy - dependent terms of str.constants
C
            IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
               CALL STRCC(ERYD,.FALSE.)
C
            ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C
               CALL TAU_TB_1K(ERYD,P,MSSQ,TAUQ,.TRUE.,KVEC)
C
            END IF
C
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
            DO IK1 = 1,NK1
               IF ( MOD(IK1,100).EQ.0 ) WRITE (6,'(10X,A,I5,A)')
     &               'outer k-loop: ',IK1,' k-points done'
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
               DO IK2 = 1,NK2
                  IF ( MOD(IK1,100).EQ.0 ) WRITE (6,'(10X,A,I5,A)')
     &                  'inner k-loop: ',IK2,' k-points done'
C
C---------------------------------------- CREATE ->K IN  UNITS OF  2PI/A
C
                  IF ( MODE.EQ.'E-k-plot  ' .OR. MODE.EQ.'E-plots   ' )
     &                 THEN
                     KVEC(1:3) = KTAB(1:3,IK2) + KSHIFT
                  ELSE
                     SCL1 = DBLE(IK1-1)/DBLE(NK1-1)
                     SCL2 = DBLE(IK2-1)/DBLE(NK2-1)
                     KVEC(1:3) = SCL1*VECK1(1:3) + SCL2*VECK2(1:3)
     &                           + KSHIFT + VECKA(1:3)
                  END IF
C
C
                  IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C
                     IF ( IREL.NE.2 ) THEN
C
                        CALL STRSET(1,KVEC,TAUK,MAUX,P)
C
                        CALL SETKKR(NQ,NKMQ,IND0Q,MAUX,MSSQ,NQMAX,NKKR,
     &                              NKMMAX)
C
                        CALL CINVLU(MAUX,TAUK,NKKR,NKKR)
C
                     ELSE
C
                        DO IS = 1,2
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  combine IQ-loops and replace use of   MQS
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                           DO IQ = 1,NQ
                              NLM_LOC = NLMQ(IQ)
                              DO J = 1,NLM_LOC
                                 CALL ZCOPY(NLM_LOC,MSSQ(1,J,IQ),1,
     &                              MQS(1,J,IQ,1),1)
                                 CALL ZCOPY(NLM_LOC,
     &                              MSSQ(NLM_LOC+1,NLM_LOC+J,IQ),1,
     &                              MQS(1,J,IQ,2),1)
                              END DO
                           END DO
C
                           CALL STRSET(1,KVEC,TAUK,MAUX,P)
C
                           DO IQ = 1,NQ
                              I1 = IND0Q(IQ) + 1
                              N = NLMQ(IQ)
                              DO J = 1,N
                                 J1 = IND0Q(IQ) + J
                                 CALL ZAXPY(N,C1,MQS(1,J,IQ,IS),1,
     &                              MAUX(I1,J1),1)
                              END DO
                           END DO
C
                           CALL CINVLU(MAUX,TAUK,NKKR,NKKR)
C
                           TAUKSP(1:NKKR,1:NKKR,IS)
     &                        = TAUK(1:NKKR,1:NKKR)
C
                        END DO
C
                     END IF
C
                  ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C
                     CALL TAU_TB_1K(ERYD,P,MSSQ,TAUQ,.FALSE.,KVEC)
C
                  ELSE
C
                     CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
                  END IF
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
                  IERR = 0
                  DO IQ = IQBOT,IQTOP
C
                     N = NKMQ(IQ)
                     M = NKMMAX
                     NLM_LOC = NKMQ(IQ)/2
                     I1 = IND0Q(IQ) + 1
C
C -------------------------------- get the matrix TAUKQQ = TAU(k)[IQ,IQ]
C ------- this has to refer to the local axis ------ rotate if necessary
C ------------------ calling <ROTATE> it is implicitely assumed that the
C -------------------- local axis has the same orientation for all sites
C
                     IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
                        IF ( IREL.NE.2 ) THEN
                           IF ( KMROT.NE.0 ) THEN
C
                              DO J = 1,N
                                 J1 = IND0Q(IQ) + J
                                 CALL ZCOPY(N,TAUK(I1,J1),1,W1(1,J),1)
                              END DO
C
                              CALL ROTATE(W1,'G->L',TAUKQQL,N,
     &                           DROTQ(1,1,IQ),M)
C
                           ELSE
C
                              DO J = 1,N
                                 J1 = IND0Q(IQ) + J
                                 CALL ZCOPY(N,TAUK(I1,J1),1,TAUKQQL(1,J)
     &                              ,1)
                              END DO
C
                           END IF
                        ELSE IF ( IREL.EQ.2 ) THEN
                           DO J = 1,NLM_LOC
                              J1 = IND0Q(IQ) + J
                              CALL ZCOPY(N,TAUKSP(I1,J1,1),1,
     &                           TAUKQQL(1,J),1)
                              CALL ZCOPY(N,TAUKSP(I1,J1,2),1,
     &                           TAUKQQL(NLM_LOC+1,NLM_LOC+J),1)
                           END DO
                        END IF
C
                     ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C
                        IF ( KMROT.NE.0 ) THEN
C
                           CALL ROTATE(TAUQ(1,1,IQ),'G->L',TAUKQQL,N,
     &                                 DROTQ(1,1,IQ),M)
C
                        ELSE
                           TAUKQQL(1:N,1:N) = TAUQ(1:N,1:N,IQ)
C
                        END IF
                     END IF
C
                     CALL TAUGFCONV(MSSQ(1,1,IQ),TAUKQQL(1,1),TMAT)
C
                     TAUKQQL(1:NKMMAX,1:NKMMAX)
     &                  = TMAT(1:NKMMAX,1:NKMMAX)
C
                     DO IBSFOP = 1,NBSFOP
C
                        CALL ZGEMM('N','N',N,N,N,C1,FCCQL(1,1,IQ,IBSFOP)
     &                             ,M,TAUKQQL,M,C0,W1,M)
C
                        IF ( NCPA.NE.0 ) THEN
                           CWGT = 1D0
                           DO J = 1,N
                              CALL ZAXPY(N,CWGT,FFTAUQL(1,J,IQ,IBSFOP),
     &                           1,W1(1,J),1)
                           END DO
                        END IF
C
                        CWGT = -1D0
                        DO J = 1,N
                           CALL ZAXPY(N,CWGT,FJQL(1,J,IQ,IBSFOP),1,
     &                                W1(1,J),1)
                        END DO
C
                        IF ( IBSFOP.EQ.1 ) THEN
                           TRC = 0.0D0
                           DO I = 1,N
                              TRC = TRC + DIMAG(W1(I,I))
                           END DO
                           IF ( TRC.LT.0.0D0 ) THEN
                              IERR = IERR + 1
                              W1(1:N,1:N) = C0
                           END IF
                           TRC = 0.0D0
                           DO I = 1,N
                              TRC = TRC + DIMAG(W1(I,I))
                           END DO
C
                           AKE(IK1,IK2,IE,IQ,3,IBSFOP)
     &                        = AKE(IK1,IK2,IE,IQ,3,IBSFOP) + TRC
C
                           IF ( IREL.GE.3 ) THEN
C
                              CALL CHANGEREP(NKM,NKMMAX,W1,'REL>CLM',W2)
C
                              DO IS = 1,2
                                 TRC = 0.0D0
                                 DO I = (IS-1)*NLM_LOC + 1,IS*NLM_LOC
                                    TRC = TRC + DIMAG(W2(I,I))
                                 END DO
                                 AKE(IK1,IK2,IE,IQ,IS,IBSFOP)
     &                              = AKE(IK1,IK2,IE,IQ,IS,IBSFOP) + TRC
                              END DO
C
                           ELSE IF ( IREL.EQ.2 ) THEN
                              DO IS = 1,2
                                 TRC = 0.0D0
                                 DO I = (IS-1)*NLM_LOC + 1,IS*NLM_LOC
                                    TRC = TRC + DIMAG(W1(I,I))
                                 END DO
                                 AKE(IK1,IK2,IE,IQ,IS,IBSFOP)
     &                              = AKE(IK1,IK2,IE,IQ,IS,IBSFOP) + TRC
                              END DO
                           END IF
                        ELSE
                           TRC = 0.0D0
                           DO I = 1,N
                              TRC = TRC + DIMAG(W1(I,I))
                           END DO
                           AKE(IK1,IK2,IE,IQ,3,IBSFOP)
     &                        = AKE(IK1,IK2,IE,IQ,3,IBSFOP) + TRC
                        END IF
                     END DO
                  END DO
                  IF ( IERR.GT.0 ) WRITE (6,99018)
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
               END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
            END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
         END IF
      END DO
CEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C=======================================================================
C               loop over the energy windows separated by   EPHOT
C=======================================================================
      IF ( EPHOT.GT.1D-3 .AND. IELOOP.EQ.1 ) THEN
         IELOOP = IELOOP + 1
         ESHIFT = EPHOT
         GOTO 100
      END IF
C
C=======================================================================
      IF ( MPI ) THEN
         CALL DRV_MPI_BARRIER
         ALLOCATE (RWKT(NK1))
         DO IK2 = 1,NK2
            DO IE = 1,NE
               DO IQ = 1,NQ
                  DO IS = 1,3
                     DO IBSFOP = 1,NBSFOP
                        CALL DRV_MPI_REDUCE_R(AKE(1,IK2,IE,IQ,IS,IBSFOP)
     &                     ,RWKT(1),NK1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      IF ( MPI_ID.EQ.0 ) THEN
C
C=======================================================================
C                     E-k-plot  and  k-k-plot
C=======================================================================
C
         IF ( MODE.EQ.'E-k-plot  ' .OR. MODE.EQ.'k-k-plot  ' ) THEN
C
C-----------------------------------------------------------------------
C                        write    AKE: original version
C-----------------------------------------------------------------------
C
            DO IE = 1,NE
               DO IK1 = 1,NK1
                  WRITE (IFILBSF,'(E12.5)')
     &                   (((AKE(IK1,IK2,IE,IQ,IS,1),IK2=1,NK2),IQ=IQBOT,
     &                   IQTOP),IS=1,2)
               END DO
            END DO
C
            DO IE = 1,NE
               DO IK1 = 1,NK1
                  WRITE (IFILBSF,'(E12.5)')
     &                   (((AKE(IK1,IK2,IE,IQ,IS,1),IK2=1,NK2),IQ=IQBOT,
     &                   IQTOP),IS=3,3)
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                        write    AKE: including spin polarisation
C-----------------------------------------------------------------------
C
            DO IE = 1,NE
               DO IBSFOP = 1,NBSFOP
                  DO IK1 = 1,NK1
                     WRITE (IFILSP,'(E12.5)')
     &                      ((AKE(IK1,IK2,IE,IQ,3,IBSFOP),IK2=1,NK2),
     &                      IQ=IQBOT,IQTOP)
                  END DO
               END DO
            END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT ) THEN
               DO IQ = IQBOT,IQTOP
                  DO IS = 1,3
                     DO IBSFOP = 1,NBSFOP
                        WRITE (IFILBUILDBOT,99020)
     &                         ROUTINE(1:LEN_TRIM(ROUTINE)),IQ,IS,
     &                         IBSFOP,
     &                         (((AKE(IK1,IK2,IE,IQ,3,IBSFOP),IK1=1,
     &                         MIN(3,NK1)),IK2=1,MIN(3,NK2)),IE=1,
     &                         MIN(3,NE))
                     END DO
                  END DO
               END DO
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            CLOSE (IFILBSF)
            CLOSE (IFILSP)
C
C=======================================================================
C                          E-plots
C=======================================================================
C
         ELSE IF ( MODE.EQ.'E-plots   ' ) THEN
C
C ----------------------------------------------------------------------
C                           E-plots
C ----------------------------------------------------------------------
C
            ALLOCATE (EPLOT(NEMAX),APLOT(NEMAX,2,0:2))
            AMAX = 0D0
            AMIN = 0D0
            DO IE = 1,NE
               EPLOT(IE) = DREAL(ETAB(IE,1))
               DO IK = 1,NK2
                  DO IS = 1,2
                     ASUM = 0D0
                     DO IQ = IQBOT,IQTOP
                        ASUM = ASUM + AKE(1,IK,IE,IQ,IS,1)
                     END DO
                     APLOT(IE,IK,IS) = ASUM
                  END DO
                  APLOT(IE,IK,0) = APLOT(IE,IK,1) + APLOT(IE,IK,2)
                  AMAX = MAX(AMAX,APLOT(IE,IK,0))
               END DO
            END DO
C
            YTXT1 = 'A!sB!N(!m{1}k!m{2}!M{1}!S!AR!N!M{2},E)'
            LYTXT1 = 38
            YTXT2 = ' '
            LYTXT2 = 1
C
            CALL XMGRHEAD(DATSET,LDATSET,'BLOCHSF',7,' ',0,FILNAM,80,
     &                    LFILNAM,IFIL,1,EPLOT(1),1,EPLOT(NE),1,AMIN,1,
     &                    AMAX,1,AMIN,1,AMAX,1,'energy (eV)',11,YTXT1,
     &                    LYTXT1,YTXT2,LYTXT2,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'Bloch spectral function '//
     &                    'A!sB!N(!m{1}k!m{2}!M{1}!S!AR!N!M{2},E)',62,
     &                    .FALSE.)
C
            LEG(1) = 'tot'
            LEG(2) = 'up'
            LEG(3) = 'down'
C
            CALL XMGRLEGEND(IFIL,1,3,3,LEG,LEG)
C
            WGT = 1.0D0
            DO IK = 1,NK2
               DO IS = 0,2
                  CALL XMGRTABLE(0,(IK-1)*3+IS,EPLOT,APLOT(1,IK,IS),WGT,
     &                           NE,IFIL)
               END DO
            END DO
C
            WRITE (6,*) '  '
            WRITE (6,99019) 'Bloch spectral function   ',
     &                      FILNAM(1:LFILNAM)
            CLOSE (IFIL)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT ) THEN
               DO IK = 1,MIN(3,NK2)
                  DO IS = 0,2
                     WRITE (IFILBUILDBOT,99021)
     &                      ROUTINE(1:LEN_TRIM(ROUTINE)),IK,IS,
     &                      (APLOT(IE,IK,IS),IE=1,MIN(3,NE))
                  END DO
               END DO
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         ELSE
C
            CALL STOP_MESSAGE(ROUTINE,'MODE not known')
C
         END IF
C
      END IF
      IF ( MPI ) CALL DRV_MPI_BARRIER
      DEALLOCATE (FJQL,TAUK,MAUX,DMAMC,FCCQL)
      DEALLOCATE (DTILT,TAUQL,MSSQL,W1,W2,FC,KP,AKE)
      DEALLOCATE (FFTAUQL,TAUKQQL,DMATT)
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
C=======================================================================
99001 FORMAT (10X,A10,4(F10.5,A))
99002 FORMAT (/,10X,'energy loop ',I3,/)
99003 FORMAT (5X,'IE=',I4,'  E=',2F10.4,
     &        '    all matrices set up  -  enter k-loop ')
99004 FORMAT (A10,2I10)
99005 FORMAT (10X,A)
99006 FORMAT (8F10.4)
99007 FORMAT (A10,8F10.4)
99008 FORMAT (80('#'))
99009 FORMAT ('#DATASET  ',A,/,'MODE      EK-REL    ',/,'NE        ',
     &        I10,/,'NK        ',I10)
99010 FORMAT (/,1X,79('*'),/,24X,'Bloch spectral function  A(k,E)',/,1X,
     &        79('*'),//,10X,'E(k)-relation mode for: ',/)
99011 FORMAT (10X,'Fermi energy       ',F10.6,' Ry',/,10X,
     &        'energy range       ',F10.6,'   --- ',F10.6,' Ry',/,10X,
     &        'with  NE =',I5,'  NK =',I5,/)
99012 FORMAT (/,1X,79('*'),/,24X,'Bloch spectral function  A(k,E)',/,1X,
     &        79('*'),//,10X,'constant E mode for: ',/)
99013 FORMAT (10X,'Fermi energy       ',F10.6,' Ry',/,10X,
     &        'fixed energy       ',2F10.6,//,10X,'k-space parameters:',
     &        /,10X,'NK1 =',I5,'    K1 =',2X,'(',3F8.3,' )   |K1| =',
     &        F7.3,/,10X,'NK2 =',I5,'    K2 =',2X,'(',3F8.3,
     &        ' )   |K2| =',F7.3,/)
99014 FORMAT (10X,'vector 2 had to be orthogonalized to vector 1 ',/,
     &        10X,'NK2 =',I5,'    K2 =',2X,'(',3F8.3,' )   |K2| =',F7.3,
     &        /)
99015 FORMAT ('#DATASET  ',A,/,'MODE      CONST-E   ',/,'NK1       ',
     &        I10,/,'NK2       ',I10/,'ERYD      ',2F10.5,/,80('#'))
99016 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*   *****  *      ****   ****  *    *        ****   ******   *'
     &  ,/,10X,
     &  '*   *    * *     *    * *    * *    *       *    *  *        *'
     &  ,/,10X,
     &  '*   *    * *     *    * *      *    *       *       *        *'
     &  ,/,10X,
     &  '*   *****  *     *    * *      ******  ***   ****   ****     *'
     &  ,/,10X,
     &  '*   *    * *     *    * *      *    *            *  *        *'
     &  ,/,10X,
     &  '*   *    * *     *    * *    * *    *       *    *  *        *'
     &  ,/,10X,
     &  '*   *****  *****  ****   ****  *    *        ****   *        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99017 FORMAT (10X,'k-dir',I3,':   KA = (',3F5.2,' )    KE = (',3F5.2,
     &        ' )')
99018 FORMAT ('WARNING <BLOCHBSF>: Spectral function negative')
99019 FORMAT (/,10X,A,/,10X,' written to the file ',A)
99020 FORMAT ('# BUILDBOT: ',A,':  Bloch spectral function A',
     &        '  for IQ =',I5,'  IS =',I5,'  IBSFOP =',I5,/,(1PE22.14))
99021 FORMAT ('# BUILDBOT: ',A,':  Bloch spectral function A',
     &        '  for IK =',I5,'  IS =',I5,/,(1PE22.14))
      END
