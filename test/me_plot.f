C*==me_plot.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_PLOT(NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to plot     DIPOLE MATRIX ELEMENTS                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_SITES,ONLY:IQAT
      USE MOD_FILES,ONLY:LDATSET,LRECREAL8,LSYSTEM,SYSTEM,DATSET,IPRINT,
     &    IOTMP,FOUND_SECTION,FOUND_REAL,FOUND_INTEGER
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EFERMI,NETAB
      USE MOD_TYPES,ONLY:NT,NTMAX,NCPLWFMAX,RHOSPN,RHOCHR,BT,VT,IMT,
     &    LTXT_T,TXT_T,CTL
      USE MOD_ANGMOM,ONLY:NKMMAX,NLQ,NKM,NLM,MZAZB,MZBZA,MIRR_2,MEZJ,
     &    MEZZ,SSST,MSST,TSST,MIRR_3,MIRR_4,NPOL,NPOLMAX
      USE MOD_RMESH,ONLY:FULLPOT,NMMAX,NRMAX,R,JRWS,R2DRDI
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C*--ME_PLOT21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_PLOT')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      INTEGER NCSTMAX
C
C Local variables
C
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),BFINT(:,:),BT0(:,:),C,DEPHOT,
     &       ECOR(NCSTMAX),EMAX,EMIN,EPHOT,EPHOTMAX,EPHOTMIN,
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),JF,JI,
     &       JIKM(NKMMAX),JQN,MIKM(NKMMAX),MJ,MJF,MJI,RARG,SCL,
     &       SEFINNORM(NTMAX),SZCOR(NCSTMAX),VFINT(:,:),VT0(:,:),XMAX,
     &       XMIN,XX(:),YMAX,YMAX0,YMAX1,YMIN,YMIN0,YMIN1,YY(:,:,:)
      LOGICAL CALCINT,GETIRRSOL,ME_CC_BRA_RWF,NONMAG,TRANS,UDT2
      CHARACTER*1 CHPOL(5)
      COMPLEX*16 CINTB(NRMAX),CINTV(NRMAX),ERYD,ERYDF,ERYDI,
     &           NORMF(NKMMAX),NORMI(NKMMAX),P,SEBFIN(NEMAX),
     &           SEBTFIN(:,:),SEVFIN(NEMAX),SEVTFIN(:,:)
      CHARACTER*2 CL
      CHARACTER*80 FILNAM
      CHARACTER*4 FUNTXTJ
      CHARACTER*1 FUNTXTL
      INTEGER I,IA_ERR,IC0,ICST,IE,IECURR,IEPHOT,IFIL,IFILF,IFILI,IG,
     &        IKM,IKMCOR(NCSTMAX,2),IKMCPLWFCOR(NCPLWFMAX),IKMF,IKMF1,
     &        IKMF2,IKMI,IKMI1,IKMI2,ILEG,IM,IOL,IPOL,IRTOP,IS,IT,ITSEL,
     &        ITXRAY,IWME,IY,IZERO(NCSTMAX),J,K,KAP,KAPCOR(NCSTMAX),
     &        KIKM(NKMMAX),KSEFIN,L,LCXRAY(NTMAX),LF,LFILNAM,LI,
     &        LIKM(NKMMAX),LISEL,LL,LSTRF,LSTRI,LSTRQN,LXTXT,LYTXT,
     &        MJM05,MM05COR(NCSTMAX),MPOL(3),NC,NCST,NCXRAY(NTMAX),NE,
     &        NEPHOT,NGRAPH,NK,NKPCOR(NCSTMAX),NL,NSOL,NXX,NYY
      CHARACTER*4 INITIALSTATE,STRQN
      CHARACTER*40 LEG(:,:),XTXT,YTXT0,YTXT1
      CHARACTER*3 MEFORM
      LOGICAL RNON0
      CHARACTER*20 STR20,STRF,STRI
C
C*** End of declarations rewritten by SPAG
C
      DATA CHPOL/'+','-','z','x','y'/,KSEFIN/0/
      DATA MPOL/ + 1, - 1,0/,CALCINT/.TRUE./
C
      ALLOCATABLE XX,YY,LEG,VFINT,BFINT,VT0,BT0
      ALLOCATABLE SEVTFIN,SEBTFIN
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (VT0(NRMAX,NTMAX),BT0(NRMAX,NTMAX))
      ALLOCATE (SEVTFIN(NRMAX,NTMAX),SEBTFIN(NRMAX,NTMAX))
      ALLOCATE (VFINT(NRMAX,NTMAX),BFINT(NRMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: VFINT')
C
      ALLOCATE (MZAZB(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MZBZA(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MIRR_2(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_3(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_4(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZAZB')
C
      IF ( FULLPOT ) THEN
         NC = NCPLWFMAX
      ELSE
         NC = 2
      END IF
      IF ( NC.LE.0 ) CALL STOP_MESSAGE(ROUTINE,'array sizes NC ???')
C
      WRITE (6,99007)
C
      NL = NINT(SQRT(DBLE(NKM/2)))
      NK = 2*NL - 1
      IKM = 0
      DO K = 1,NK
         L = K/2
         IF ( MOD(K,2).EQ.1 ) THEN
            KAP = -L - 1
         ELSE
            KAP = +L
         END IF
         JQN = ABS(KAP) - 0.5D0
C
         DO MJM05 = NINT(-JQN-0.5D0),NINT(JQN-0.5D0)
            MJ = DBLE(MJM05) + 0.5D0
            IKM = IKM + 1
            LIKM(IKM) = L
            KIKM(IKM) = KAP
            JIKM(IKM) = JQN
            MIKM(IKM) = MJ
C
            IF ( IPRINT.GT.0 ) WRITE (6,'(a,3i4,2f5.1)') '#########',
     &                                IKM,LIKM(IKM),KIKM(IKM),JIKM(IKM),
     &                                MIKM(IKM)
         END DO
      END DO
C
      MJI = 0.5D0
      IPOL = 3
C
C     GETIRRSOL = .TRUE.
      GETIRRSOL = .FALSE.
C
      IFIL = IOTMP
C
      NONMAG = .FALSE.
      LISEL = -1
      ITSEL = 1
      INITIALSTATE = 'BAND'
C      INITIALSTATE = 'CORE'
      ERYDI = 0D0
      ERYDF = 0D0
      EPHOT = 100D0
      NEPHOT = 1
C
      YMIN = 0D0
      YMAX = 0D0
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_FIND_KEYWORD('NONMAG',NONMAG)
         CALL SECTION_SET_INTEGER('ITSEL',ITSEL,9999,0)
         CALL SECTION_SET_INTEGER('IPRINT',IPRINT,9999,0)
         CALL SECTION_SET_INTEGER('SEFIN',KSEFIN,9999,0)
         CALL SECTION_SET_STRING('INITIALSTATE',INITIALSTATE,'9999',0)
         IT = ITSEL
C
         IF ( INITIALSTATE.EQ.'CORE' ) THEN
C
            CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
            STRQN = CL
            LSTRQN = 2
C
         ELSE IF ( INITIALSTATE.EQ.'BAND' ) THEN
C
            CALL SECTION_SET_INTEGER('L',LISEL,NLQ(IQAT(1,IT))-1,0)
C
            IF ( FOUND_INTEGER ) THEN
               STRQN = FUNTXTL(LISEL)
               LSTRQN = 1
            ELSE
               STRQN = 'VB'
               LSTRQN = 2
               LISEL = -1
            END IF
C
C --------------------------------- read photon energy information in eV
C
            CALL SECTION_SET_INTEGER('NEPHOT',NEPHOT,9999,0)
            CALL SECTION_SET_REAL('EPHOTMIN',EPHOTMIN,9999D0,0)
            UDT2 = FOUND_REAL
            CALL SECTION_SET_REAL('EPHOTMAX',EPHOTMAX,9999D0,0)
            IF ( FOUND_INTEGER .OR. UDT2 .OR. FOUND_REAL ) THEN
               IF ( FOUND_INTEGER .AND. UDT2 .AND. .NOT.FOUND_REAL )
     &              THEN
                  WRITE (6,99001)
                  CALL STOP_MESSAGE(ROUTINE,' ')
               END IF
               WRITE (6,99002) EPHOTMIN,EPHOTMAX,NEPHOT
            ELSE
C
               CALL SECTION_SET_REAL('EPHOT',EPHOT,9999D0,0)
               IF ( .NOT.FOUND_REAL ) THEN
                  CALL SECTION_SET_REAL('EPHOTMIN',EPHOT,9999D0,0)
                  IF ( .NOT.FOUND_REAL )
     &                  CALL SECTION_SET_REAL('EPHOTMAX',EPHOT,9999D0,0)
               END IF
C
               WRITE (6,99003) EPHOT
C
            END IF
C
            IF ( KSEFIN.NE.0 .AND. NEPHOT.GT.1 ) THEN
               KSEFIN = 0
               WRITE (6,99016)
            END IF
C
         ELSE
            WRITE (6,99004) INITIALSTATE
            CALL STOP_MESSAGE(ROUTINE,'INITIALSTATE')
         END IF
C
      END IF
C
      MEFORM = 'ADA'
C
      IWME = 0
      CALL AME_INIT(MEFORM,IWME)
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
      C = CTL(IT,1)
C
C=======================================================================
C
      IFILI = 87
      IFILF = 88
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFILI,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
      OPEN (UNIT=IFILF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
C-----------------------------------------------------------------------
C
      WRITE (6,99011) 'NPOL      ',NPOL,(CHPOL(I),I=1,NPOL)
      WRITE (6,99010) 'NLM       ',NLM
      WRITE (6,99010) 'NKM       ',NKM
      WRITE (6,99014) ITSEL
C
C ======================================================================
C            find array sizes for plotting functions and allocate
C ======================================================================
      IF ( INITIALSTATE.EQ.'CORE' ) THEN
         IKMI1 = 2*LCXRAY(IT)**2 + 1
         IKMI2 = 2*(LCXRAY(IT)+1)**2
      ELSE IF ( LISEL.GE.0 ) THEN
         IKMI1 = 2*LISEL**2 + 1
         IKMI2 = 2*(LISEL+1)**2
      ELSE
         IKMI1 = 1
         IKMI2 = NKM
      END IF
      IKMF1 = 1
      IKMF2 = NKM
C
      IS = -1
      DO IKMI = IKMI1,IKMI2
C----------------- select initial state with magnetic quantum number MJI
         IF ( .NOT.RNON0(MIKM(IKMI)-MJI) ) THEN
            LI = LIKM(IKMI)
            JI = JIKM(IKMI)
C
            DO IKMF = IKMF1,IKMF2
               LF = LIKM(IKMF)
               JF = JIKM(IKMF)
               MJF = MIKM(IKMF)
C
C------------------------------------------ apply dipole selection rules
               IF ( ABS(LI-LF).NE.1 ) THEN
                  TRANS = .FALSE.
               ELSE IF ( NINT(MJF-MJI-MPOL(IPOL)).NE.0 ) THEN
                  TRANS = .FALSE.
               ELSE IF ( NONMAG .AND. ABS(NINT(JF-JI)).GT.1 ) THEN
                  TRANS = .FALSE.
               ELSE
                  TRANS = .TRUE.
               END IF
               IF ( TRANS ) THEN
C
                  IS = IS + 1
                  NYY = IS + 1
C
               END IF
C
            END DO
         END IF
      END DO
C
      IF ( NEPHOT.EQ.1 ) THEN
         NE = NEMAX
      ELSE
         NE = NEPHOT
      END IF
C
      NXX = MAX(NE,NRMAX)
C
      ALLOCATE (XX(NXX),YY(NXX,0:NYY-1,0:1),LEG(NYY,0:1),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: XX')
C
C ======================================================================
C                     suppress magnetic field - if requested
C ======================================================================
C
      IF ( NONMAG ) THEN
         WRITE (6,99015)
         DO I = 1,IRTOP
            BT(I,IT) = 0D0
         END DO
      END IF
C
C ======================================================================
C                     copy final state potential
C ======================================================================
C
      DO IT = 1,NT
C
         CALL DCOPY(NRMAX,VT(1,IT),1,VT0(1,IT),1)
         CALL DCOPY(NRMAX,BT(1,IT),1,BT0(1,IT),1)
         CALL DCOPY(NRMAX,VT(1,IT),1,VFINT(1,IT),1)
         CALL DCOPY(NRMAX,BT(1,IT),1,BFINT(1,IT),1)
C
         IF ( KSEFIN.NE.0 ) THEN
            IM = IMT(IT)
            CALL RRADINT(IM,R2DRDI(1,IM),SEFINNORM(IT))
         END IF
      END DO
C
      IT = ITSEL
C
C ======================================================================
C                       band energy and photon loops
C ======================================================================
C
      IF ( INITIALSTATE.EQ.'CORE' ) THEN
         IF ( NEPHOT.NE.1 ) THEN
            WRITE (6,99006) NEPHOT
            NEPHOT = 1
         END IF
      END IF
C
      DO IE = 1,NETAB(1)
         ERYD = ETAB(IE,1)
C
         DO IEPHOT = 1,NEPHOT
C
            IF ( NEPHOT.NE.1 ) THEN
               DEPHOT = (EPHOTMAX-EPHOTMIN)/DBLE(NEPHOT-1)
               EPHOT = EPHOTMIN + DEPHOT*(IEPHOT-1)
            END IF
C
C ======================================================================
C                        initial BAND state
C ======================================================================
C
            IF ( INITIALSTATE.NE.'CORE' ) THEN
C
               ERYDI = ERYD
               IF ( LISEL.GE.0 ) THEN
                  IKMI1 = 2*LISEL**2 + 1
                  IKMI2 = 2*(LISEL+1)**2
               ELSE
                  IKMI1 = 1
                  IKMI2 = NKM
               END IF
C
               IF ( (IE.EQ.1 .AND. IEPHOT.EQ.1) .OR. IPRINT.GT.0 ) THEN
                  WRITE (6,99012) INITIALSTATE
                  WRITE (6,99013) 'initial',ERYDI,IKMI1,IKMI2
               END IF
C
               IF ( FULLPOT ) THEN
C
                  CALL FPSSITE(1,1,IFILI,GETIRRSOL,ERYDI,P,IPRINT,TSST,
     &                         MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
               ELSE
C
                  CALL SSITE(1,1,IFILI,CALCINT,GETIRRSOL,ERYDI,P,IPRINT,
     &                       NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
               END IF
C
               DO IKMI = IKMI1,IKMI2
                  NORMI(IKMI) = SQRT(MEZZ(IKMI,IKMI,IT,1))
               END DO
C
C ======================================================================
C                        initial CORE state
C ======================================================================
C
            ELSE IF ( IE.EQ.1 ) THEN
C
               ITXRAY = IT
               NCST = 4*LCXRAY(IT) + 2
C
               WRITE (6,99012) INITIALSTATE,NCXRAY(IT),LCXRAY(IT)
C
               IF ( FULLPOT ) THEN
C
                  CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,
     &                      NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,BCORS,
     &                      NCSTMAX)
C
                  DO ICST = 1,NCST
                     IF ( ABS(MM05COR(ICST)+0.5D0).GT.LCXRAY(IT) ) THEN
                        NSOL = 1
                     ELSE
                        NSOL = 2
                     END IF
                     DO J = 1,NCPLWFMAX
                        IKMCPLWFCOR(J) = 0
                     END DO
                     DO J = 1,NSOL
                        IKMCPLWFCOR(J) = IKMCOR(ICST,J)
                     END DO
C
                     WRITE (IFILI,REC=IKMCOR(ICST,1)+(IT-1)*NKM) IT,
     &                      'COR',IKMCOR(ICST,1),JRWS(IM),
     &                      (IKMCPLWFCOR(J),J=1,NCPLWFMAX),NSOL,
     &                      ((DCMPLX(GCOR(I,K,ICST),0.0D0),I=1,JRWS(IM))
     &                      ,K=1,NSOL),
     &                      ((DCMPLX(FCOR(I,K,ICST),0.0D0),I=1,JRWS(IM))
     &                      ,K=1,NSOL)
C
                  END DO
C
               ELSE
C
                  CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,
     &                      NKPCOR,IKMCOR,IZERO,ITXRAY,BCOR,BCORS,
     &                      NCSTMAX)
C
                  DO ICST = 1,NCST
                     WRITE (IFILI,REC=IKMCOR(ICST,1)+(IT-1)*NKM) IT,
     &                      LCXRAY(IT),DBLE(MM05COR(ICST)+0.5D0),
     &                      NKPCOR(ICST),'COR',KAPCOR(ICST),
     &                      IKMCOR(ICST,1),
     &                      (DCMPLX(GCOR(I,1,ICST),0.0D0),
     &                      DCMPLX(FCOR(I,1,ICST),0.0D0),I=1,IRTOP),
     &                      - KAPCOR(ICST) - 1,IKMCOR(ICST,2),
     &                      (DCMPLX(GCOR(I,2,ICST),0.0D0),
     &                      DCMPLX(FCOR(I,2,ICST),0.0D0),I=1,IRTOP)
                  END DO
C
               END IF
C
               IKMI1 = IKMCOR(1,1)
               IKMI2 = IKMCOR(NCST,1)
               ERYDI = DCMPLX(ECOR(1))
C
               WRITE (6,99013) 'initial',ERYDI,IKMI1,IKMI2
               DO IKMI = IKMI1,IKMI2
                  NORMI(IKMI) = 1D0
               END DO
C
            END IF
C
C ======================================================================
C                         final BAND state
C ======================================================================
C
            IF ( INITIALSTATE.EQ.'CORE' ) THEN
               ERYDF = ERYD
            ELSE
               ERYDF = ERYD + EPHOT/RY_EV
            END IF
C
            IF ( KSEFIN.NE.0 ) THEN
C
               IECURR = IE
C
               CALL SEFIN(IECURR,ERYD,EFERMI,NT,RHOCHR,RHOSPN,KSEFIN,
     &                    SEBTFIN,SEVTFIN,JRWS,IMT,NRMAX,NMMAX,NTMAX)
C
               DO I = 1,NRMAX
                  VFINT(I,IT) = VT(I,IT) + DREAL(SEVTFIN(I,IT))
                  BFINT(I,IT) = BT(I,IT) + DREAL(SEBTFIN(I,IT))
C
                  CINTV(I) = SEVTFIN(I,IT)*R2DRDI(I,IM)
                  CINTB(I) = SEBTFIN(I,IT)*R2DRDI(I,IM)
               END DO
C
               CALL CRADINT(IM,CINTV,SEVFIN(IECURR))
               CALL CRADINT(IM,CINTB,SEBFIN(IECURR))
C
            END IF
C
            IKMF1 = 1
            IKMF2 = NKM
            IF ( (IE.EQ.1 .AND. IEPHOT.EQ.1) .OR. IPRINT.GT.0 )
     &           WRITE (6,99013) 'final  ',ERYDF,IKMF1,IKMF2
C
            CALL RVECCOP(NRMAX*NT,VFINT,VT)
            CALL RVECCOP(NRMAX*NT,BFINT,BT)
C
            IF ( FULLPOT ) THEN
C
               CALL FPSSITE(1,1,IFILF,GETIRRSOL,ERYDF,P,IPRINT,TSST,
     &                      MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            ELSE
C
               CALL SSITE(1,1,IFILF,CALCINT,GETIRRSOL,ERYDF,P,IPRINT,
     &                    NKM,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            END IF
C
            CALL RVECCOP(NRMAX*NT,VT0,VT)
            CALL RVECCOP(NRMAX*NT,BT0,BT)
C
            DO IKMF = IKMF1,IKMF2
               NORMF(IKMF) = SQRT(MEZZ(IKMF,IKMF,IT,1))
            END DO
C
C ======================================================================
C                            matrix elements
C ======================================================================
C
            ME_CC_BRA_RWF = .TRUE.
            IT = ITSEL
            C = CTL(IT,1)
C
            CALL ME_DRIVE(IKMF1,IKMF2,IFILF,ERYDF,IKMI1,IKMI2,IFILI,
     &                    ERYDI,ME_CC_BRA_RWF,IT,MZAZB,MZBZA,MIRR_2,
     &                    MIRR_3,MIRR_4,C,K_CALC_ME,MEFORM,IWME)
C
C ======================================================================
C                         copy matrix elements
C ======================================================================
C
            IS = -1
            IG = 0
            DO IKMI = IKMI1,IKMI2
C----------------- select initial state with magnetic quantum number MJI
               IF ( .NOT.RNON0(MIKM(IKMI)-MJI) ) THEN
                  LI = LIKM(IKMI)
                  JI = JIKM(IKMI)
C
                  DO IKMF = IKMF1,IKMF2
                     LF = LIKM(IKMF)
                     JF = JIKM(IKMF)
                     MJF = MIKM(IKMF)
C
C------------------------------------------ apply dipole selection rules
                     IF ( ABS(LI-LF).NE.1 ) THEN
                        TRANS = .FALSE.
                     ELSE IF ( NINT(MJF-MJI-MPOL(IPOL)).NE.0 ) THEN
                        TRANS = .FALSE.
                     ELSE IF ( NONMAG .AND. ABS(NINT(JF-JI)).GT.1 ) THEN
                        TRANS = .FALSE.
                     ELSE
                        TRANS = .TRUE.
                     END IF
                     IF ( TRANS ) THEN
                        IS = IS + 1
                        IF ( (IE.EQ.1 .AND. IEPHOT.EQ.1) .OR. 
     &                       IPRINT.GT.0 ) WRITE (6,99005) IG,IS,IKMF,
     &                       LF,NINT(2*JIKM(IKMF)),NINT(2*MJF),IKMI,LI,
     &                       NINT(2*JIKM(IKMI)),NINT(2*MJI),IPOL
C
                        IF ( NEPHOT.EQ.1 ) THEN
                           XX(IE) = DREAL(ERYD-EFERMI)*RY_EV
                           IY = IE
                        ELSE
                           XX(IEPHOT) = EPHOT
                           IY = IEPHOT
                        END IF
                        YY(IY,IS,IG) = ABS(MZAZB(IKMF,IKMI,IPOL)/(NORMI(
     &                                 IKMI)*NORMF(IKMF)))
                        YMIN = MIN(YMIN,YY(IY,IS,IG))
                        YMAX = MAX(YMAX,YY(IY,IS,IG))
                     END IF
C
                  END DO
               END IF
            END DO
C ======================================================================
C
         END DO
C ==============================================================  IEPHOT
C
      END DO
C ==================================================================  IE
C
C
C ======================================================================
C                    set up legends for curves
C ======================================================================
      IS = -1
      IG = 0
      DO IKMI = IKMI1,IKMI2
C----------------- select initial state with magnetic quantum number MJI
         IF ( .NOT.RNON0(MIKM(IKMI)-MJI) ) THEN
            LI = LIKM(IKMI)
            JI = JIKM(IKMI)
C
            IF ( INITIALSTATE.EQ.'CORE' ) THEN
               STRI = CL//'!s'//FUNTXTJ(JIKM(IKMI))//'!N'
            ELSE
               STRI = FUNTXTL(LI)//'!s'//FUNTXTJ(JIKM(IKMI))//'!N'
            END IF
            LSTRI = LEN_TRIM(STRI)
C
            DO IKMF = IKMF1,IKMF2
               LF = LIKM(IKMF)
               JF = JIKM(IKMF)
               MJF = MIKM(IKMF)
C
C------------------------------------------ apply dipole selection rules
               IF ( ABS(LI-LF).NE.1 ) THEN
                  TRANS = .FALSE.
               ELSE IF ( NINT(MJF-MJI-MPOL(IPOL)).NE.0 ) THEN
                  TRANS = .FALSE.
               ELSE IF ( NONMAG .AND. ABS(NINT(JF-JI)).GT.1 ) THEN
                  TRANS = .FALSE.
               ELSE
                  TRANS = .TRUE.
               END IF
               IF ( TRANS ) THEN
C
                  IS = IS + 1
                  NYY = IS + 1
                  ILEG = IS + 1
C
                  STRF = FUNTXTL(LF)//'!s'//FUNTXTJ(JIKM(IKMF))//'!N'
                  LSTRF = LEN_TRIM(STRF)
                  LSTRI = MIN(40-LSTRF-3,LSTRI)
C
                  LEG(ILEG,IG) = STRF(1:LSTRF)//'!AL'//STRI(1:(LSTRI-2))
C
                  WRITE (6,99005) IG,IS,IKMF,LF,NINT(2*JIKM(IKMF)),
     &                            NINT(2*MJF),IKMI,LI,NINT(2*JIKM(IKMI))
     &                            ,NINT(2*MJI),IPOL
C
               END IF
C
            END DO
         END IF
      END DO
C ======================================================================
C
      NGRAPH = 1
C
      YTXT0 = '|M!sfi!N(E)| (a.u.)'
      YTXT1 = '                   '
      LYTXT = 19
      IF ( NEPHOT.LE.1 ) THEN
         XTXT = 'energy (eV)'
         LXTXT = 11
      ELSE
         XTXT = 'photon energy (eV)'
         LXTXT = 18
      END IF
C
      IF ( NEPHOT.EQ.1 ) THEN
         NXX = NETAB(1)
         EMIN = XX(1)
         EMAX = XX(NXX)
      ELSE
         NXX = NEPHOT
         EMIN = EPHOTMIN
         EMAX = EPHOTMAX
      END IF
C
      IC0 = ICHAR('0')
      CALL STRING_CONVERT_TO_LC(INITIALSTATE)
      IF ( KSEFIN.EQ.0 ) THEN
         STR20 = STRQN(1:LSTRQN)
         LL = LSTRQN
      ELSE
         STR20 = STRQN(1:LSTRQN)//'_SE'//CHAR(IC0+KSEFIN)
         LL = LSTRQN + 4
      END IF
C
      CALL XMGRHEAD(DATSET,LDATSET,STR20,LL,TXT_T(IT),LTXT_T(IT),FILNAM,
     &              80,LFILNAM,IFIL,NGRAPH,EMIN,0,EMAX,0,YMIN,0,YMAX,1,
     &              YMIN,0,YMAX,0,XTXT,LXTXT,YTXT0,LYTXT,YTXT1,LYTXT,
     &              'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &              25+LSYSTEM,INITIALSTATE//
     &              '-band dipole transition matrix elements for '//
     &              TXT_T(IT)(1:LTXT_T(IT)),(48+LTXT_T(IT)),.FALSE.)
C
      CALL XMGRLEG1(IFIL,0,NYY,LEG,0.708D0,0.83D0)
C
      CALL XMGRCURVES(IFIL,NGRAPH,NYY,NYY,2,1,0)
C
      DO IG = 0,(NGRAPH-1)
         DO IS = 0,(NYY-1)
            CALL XMGRTABLE(IG,IS,XX,YY(1,IS,IG),1.0D0,NXX,IFIL)
         END DO
      END DO
C
      WRITE (6,99008) INITIALSTATE,FILNAM(1:LFILNAM)
      CLOSE (IFIL)
C
      IF ( KSEFIN.EQ.0 ) CALL STOP_REGULAR(ROUTINE,'done')
C=======================================================================
C                  Self energy for final state  ---  SEFIN
C=======================================================================
C
      NGRAPH = 2
      NYY = 2
C
      YTXT0 = '!xS!0!sV,fin!N(E) (Ry)'
      YTXT1 = '!xS!0!sB,fin!N(E) (Ry)'
      LYTXT = 22
      XTXT = 'energy (eV)'
      LXTXT = 11
C
      NXX = NETAB(1)
      EMIN = XX(1)
      EMAX = XX(NXX)
C
      YMIN0 = 0D0
      YMAX0 = 0D0
      YMIN1 = 0D0
      YMAX1 = 0D0
      DO I = 1,NXX
         YY(I,0,0) = DREAL(SEVFIN(I))
         YY(I,1,0) = DIMAG(SEVFIN(I))
         YMIN0 = MIN(YMIN0,YY(I,0,0),YY(I,1,0))
         YMAX0 = MAX(YMAX0,YY(I,0,0),YY(I,1,0))
C
         YY(I,0,1) = DREAL(SEBFIN(I))
         YY(I,1,1) = DIMAG(SEBFIN(I))
         YMIN1 = MIN(YMIN1,YY(I,0,1),YY(I,1,1))
         YMAX1 = MAX(YMAX1,YY(I,0,1),YY(I,1,1))
      END DO
C
      STR20 = 'SEfin_E_SE'//CHAR(IC0+KSEFIN)
C
      CALL XMGRHEAD(DATSET,LDATSET,STR20,11,TXT_T(IT),LTXT_T(IT),FILNAM,
     &              80,LFILNAM,IFIL,NGRAPH,EMIN,0,EMAX,0,YMIN0,1,YMAX0,
     &              1,YMIN1,1,YMAX1,1,XTXT,LXTXT,YTXT0,LYTXT,YTXT1,
     &              LYTXT,'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &              ,25+LSYSTEM,
     &              'self energy  !xS!0!sfin!N(E)  for final state for '
     &              //TXT_T(IT)(1:LTXT_T(IT)),(50+LTXT_T(IT)),.FALSE.)
C
      LEG(1,0) = 'Re !xS!0!sV,fin!N'
      LEG(2,0) = 'Im !xS!0!sV,fin!N'
      LEG(1,1) = 'Re !xS!0!sB,fin!N'
      LEG(2,1) = 'Im !xS!0!sB,fin!N'
C
      CALL XMGRLEG1(IFIL,0,NYY,LEG(1,0),0.62D0,0.48D0)
      CALL XMGRLEG1(IFIL,1,NYY,LEG(1,1),0.62D0,0.83D0)
C
      CALL XMGRCURVES(IFIL,NGRAPH,NYY,NYY,2,1,0)
C
      DO IG = 0,(NGRAPH-1)
         DO IS = 0,1
            CALL XMGRTABLE(IG,IS,XX,YY(1,IS,IG),1.0D0,NXX,IFIL)
         END DO
      END DO
C
      WRITE (6,99009) 'E',FILNAM(1:LFILNAM)
      CLOSE (IFIL)
C
C-----------------------------------------------------------------------
C
      NGRAPH = 2
      NXX = IRTOP
C
      YTXT0 = '!xS!0!sV,fin!N(r) (Ry)'
      YTXT1 = '!xS!0!sB,fin!N(r) (Ry)'
      LYTXT = 22
      XTXT = 'radius r (a.u.)'
      LXTXT = 16
C
      XMIN = 0D0
      XMAX = R(IRTOP,IM)
C
      YMIN0 = 0D0
      YMAX0 = 0D0
      YMIN1 = 0D0
      YMAX1 = 0D0
      DO I = 1,IRTOP
         XX(I) = R(I,IM)
         SCL = XX(I)
         SCL = 1D0
C
         YY(I,0,0) = SCL*DREAL(SEVTFIN(I,IT))
         YY(I,1,0) = SCL*DIMAG(SEVTFIN(I,IT))
         YMIN0 = MIN(YMIN0,YY(I,0,0),YY(I,1,0))
         YMAX0 = MAX(YMAX0,YY(I,0,0),YY(I,1,0))
C
         YY(I,0,1) = SCL*DREAL(SEBTFIN(I,IT))
         YY(I,1,1) = SCL*DIMAG(SEBTFIN(I,IT))
         YMIN1 = MIN(YMIN1,YY(I,0,1),YY(I,1,1))
         YMAX1 = MAX(YMAX1,YY(I,0,1),YY(I,1,1))
      END DO
C
      STR20 = 'SEfin_r_SE'//CHAR(IC0+KSEFIN)
C
      CALL XMGRHEAD(DATSET,LDATSET,STR20,11,TXT_T(IT),LTXT_T(IT),FILNAM,
     &              80,LFILNAM,IFIL,NGRAPH,XMIN,0,XMAX,0,YMIN0,1,YMAX0,
     &              1,YMIN1,1,YMAX1,1,XTXT,LXTXT,YTXT0,LYTXT,YTXT1,
     &              LYTXT,'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &              ,25+LSYSTEM,
     &              'self energy  !xS!0!sfin!N(E)  for final state for '
     &              //TXT_T(IT)(1:LTXT_T(IT)),(50+LTXT_T(IT)),.FALSE.)
C
      CALL XMGRLEG1(IFIL,0,NYY,LEG(1,0),0.62D0,0.48D0)
      CALL XMGRLEG1(IFIL,1,NYY,LEG(1,1),0.62D0,0.83D0)
C
      CALL XMGRCURVES(IFIL,NGRAPH,NYY,NYY,2,1,0)
C
      DO IG = 0,(NGRAPH-1)
         DO IS = 0,1
            CALL XMGRTABLE(IG,IS,XX,YY(1,IS,IG),1.0D0,NXX,IFIL)
         END DO
      END DO
C
      WRITE (6,99009) 'r',FILNAM(1:LFILNAM)
      CLOSE (IFIL)
C
C=======================================================================
C
      CALL STOP_REGULAR(ROUTINE,'plot of  DIPOLE MATRIX ELEMENTS  done')
C
C=======================================================================
99001 FORMAT (/,1x,79('#'),/,10X,
     &        'information on photon energy mesh incomplete',/,10X,
     &        'input format:   EPHOTMIN=100 EPHOTMAX=600 NEPHOT=20')
99002 FORMAT (/,10X,'matrix elements for varying photon energy',/,10X,
     &        'E_phot =',F9.2,' ...',F9.2,'  eV   in ',I4,' steps',/)
99003 FORMAT (/,10X,'fixed photon energy for matrix elements',
     &        '  E_phot =',F9.2,' eV   ',/)
99004 FORMAT (/,1x,79('#'),/,10X,'from input:  INITIALSTATE = ',A,/,10X,
     &        'should be   CORE   or   BAND ')
99005 FORMAT (10x,'G=',I1,'  S=',I2,3x,'trans ',2I3,2(i2,'/2'),' <-',
     &        2I3,2(i2,'/2'),3x,'POL',i3)
99006 FORMAT (/,10X,'############## WARNING from  <MEPLOT>',/,10X,
     &        '############## initial core state  AND  NEPHOT =',I4,/,
     &        10X,'##############  NEPHOT  set to 1 ',/)
99007 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      *     *  ******  *****   *        ****   *******      *'
     &  ,/,10X,
     &  '*      **   **  *       *    *  *       *    *     *         *'
     &  ,/,10X,
     &  '*      * * * *  *       *    *  *       *    *     *         *'
     &  ,/,10X,
     &  '*      *  *  *  ***     *****   *       *    *     *         *'
     &  ,/,10X,
     &  '*      *     *  *       *       *       *    *     *         *'
     &  ,/,10X,
     &  '*      *     *  *       *       *       *    *     *         *'
     &  ,/,10X,
     &  '*      *     *  ******  *       ******   ****      *         *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99008 FORMAT (/,10X,A,'-band matrix elements written to: ',A,/)
99009 FORMAT (/,10X,'self energy S_fin(',A,') written to: ',A,/)
99010 FORMAT (10X,A10,I10)
99011 FORMAT (10X,A10,I10,5(:,'  (',A,')'))
99012 FORMAT (/,10X,'initial state type: ',A,/,:,/,10X,
     &        'quantum numbers:     n = ',I1,'  l = ',I1,/)
99013 FORMAT (/,10X,A,' state    energy  E = ',2F12.6,//,33X,'IKM =',I9,
     &        '   ...   ',I3,/)
99014 FORMAT (//,10X,'selected atom type   IT =',I3,/)
99015 FORMAT (/,10X,'spin-dependent part B of potential suppressed',/)
99016 FORMAT (/,10X,'############## WARNING from  <MEPLOT>',/,10X,
     &        '############## SEFIN set in input   AND  NEPHOT <> 0 ',/,
     &        10X,'##############  SEFIN setting ignored  ',/)
      END
