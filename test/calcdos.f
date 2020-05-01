C*==calcdos.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCDOS(NWRLOG,INTEFERMI,SPLITSS,IEPATH,NCPAFAIL,
     &                   IPRINT,ERYD,MEZZ,MEZJ,TSST,MSST,TAUT,MSSQ,TAUQ,
     &                   IECURR,WE,BCOR,BCORS,DOBS_LTX,DOBS_TX,OBS_LTX,
     &                   OBS_TX,DOBS_TX_GLO,OBS_TX_GLO,CTOTDOS)
C   ********************************************************************
C   *                                                                  *
C   * SUBROUTINE TO CALCULATE THE DOS AND MAGNETIC SPIN AND ORBITAL    *
C   * MOMENT AND HYPERFINE FIELD WITHIN AN ATOMIC CELL                 *
C   *                                                                  *
C   *  N(EF) = -1/PI * Im    Int    TRACE G(R,R,E) dE                  *
C   *                       E=0..EF       =                            *
C   *                                                                  *
C   *  M(EF) = -1/PI * Im    Int    TRACE BET*SIG*G(R,R,E) dE          *
C   *                       E=0..EF        =   =  =                    *
C   *                                                                  *
C   * THE BACKSCATTERING CONTRIBUTION AND THE FULL CRYSTAL             *
C   * CONTRIBUTIONS ARE CALCULATED FOR ALL QUANTITIES                  *
C   * CRYSTAL:          G = SUM Z*TAU*Z - SUM Z*J                      *
C   * BACKSCATTERING:   G = SUM Z*(TAU-T)*Z                            *
C   *                                                                  *
C   *  FP:   the matrix elements  MEZZ and MEZJ  are recalculated      *
C   *        in a l-projected way                                      *
C   *        the print out of the mj- and kappa-projected DOS          *
C   *        has been removed                                          *
C   *                                                                  *
C   * 20/09/15 HE major revision                                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_ID
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS,JRCRI
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_SITES,ONLY:IQAT,NQMAX,DROTQ
      USE MOD_THERMAL,ONLY:X_FT,FTET_FT,NFLUCT,TEMP_LAT
      USE MOD_ANGMOM,ONLY:IDOS,ISMT,IOMT,IHFF,IBND,ISDM,IKDS,NLM,NMEMAX,
     &    NOBSMAX,NKMMAX,NKM,NMUEMAX,NLMAX,WKM1,NCPLWF
      USE MOD_TYPES,ONLY:NTMAX,TXT_T,NVALTOT,CONC,NAT,NLT,ITBOT,ITTOP,
     &    DOBS_BS_EF_TX,DOBS_BS_EF_LTX,TOTDOS_BS_EF,JFRA,JGRA,ZGRA,ZFRA,
     &    JGLA,JFLA,ZFLA,ZGLA,IMT,NCPLWFMAX,IKMCPLWF,DOBS_LTEX
      USE MOD_CALCMODE,ONLY:IREL,ICORE,UPDATE_EFERMI,LHS_SOL_EQ_RHS_SOL,
     &    THERMAL_VIBRA_FLUCT,L_PROJECTED_ME,KMROT,MOMENTS_ROTATED,
     &    WAVE_FUNCTIONS_AVAILABLE,PROGNAME
      USE MOD_KSPACE,ONLY:IBZINT,NKTAB,NPTMAX,NPTMIN
      USE MOD_SYMMETRY,ONLY:NWEDGE,IWEDGEROT
      USE MOD_ENERGY,ONLY:NEMAX,IGRID,NETAB,EFERMI,NEPATH
      USE MOD_FILES,ONLY:SYSTEM,LSYSTEM,IFILLOG,IFILMEZZL,IFILBUILDBOT,
     &    WRBUILDBOT,IDUMMY,IFILCBWF
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      IMPLICIT NONE
C*--CALCDOS51
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CALCDOS')
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
      INTEGER NME,NOBS
      PARAMETER (NME=7,NOBS=3)
      LOGICAL CHECK_OBS_ME
      PARAMETER (CHECK_OBS_ME=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,WE
      INTEGER IECURR,IEPATH,IPRINT,NCPAFAIL,NWRLOG
      LOGICAL INTEFERMI,SPLITSS
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),OBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX),
     &       OBS_TX(0:3,NOBSMAX,NTMAX),OBS_TX_GLO(0:3,NOBSMAX,NTMAX)
      COMPLEX*16 CTOTDOS(NEMAX),DOBS_LTX(0:3,NOBSMAX,NLMAX,NTMAX),
     &           DOBS_TX(0:3,NOBSMAX,NTMAX),
     &           DOBS_TX_GLO(0:3,NOBSMAX,NTMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BDUM(3),DQ,E_BAND,MIX_DL(:),MIX_IL(:),MUEORB,MUESPN,
     &       OBS_SUM_L_TX(:,:,:),SHFTEF,SPF,TOTDOS,TOTNOS,WANG_FT,WSPIN,
     &       WSUM
      COMPLEX*16 CTOTALDOS,DOBS_LMX(:,:,:,:),DOBS_SUM_L_TX(:,:,:),
     &           MEZJL(:,:,:,:),MEZZL(:,:,:,:),MZBJA(:,:,:,:),
     &           MZBZA(:,:,:,:),TMAT(:,:),TMAT_LOC(:,:)
      INTEGER I,I0,I2,IA_ERR,IEF,IFIL_LHSB,IFIL_RHSA,IFLUCT,IFT,IL,IM,
     &        IME,IQ,IRTOP,IT,ITP,IW,IWRLOG,J,L,LMAX,M,MX,MXTOP,N
      LOGICAL INITIALIZE
      CHARACTER*20 STR20
      CHARACTER*9 STR9
      CHARACTER*10 TXTMJ
      SAVE DOBS_SUM_L_TX,OBS_SUM_L_TX
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE MZBJA,MZBZA
      ALLOCATABLE DOBS_SUM_L_TX,OBS_SUM_L_TX
      ALLOCATABLE MEZJL,MEZZL,MIX_DL,MIX_IL
      ALLOCATABLE TMAT,TMAT_LOC,DOBS_LMX
C
      ALLOCATE (TMAT(NKMMAX,NKMMAX),TMAT_LOC(NKMMAX,NKMMAX))
C
      ALLOCATE (DOBS_LMX(0:3,NOBSMAX,NLMAX,NMUEMAX))
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (OBS_SUM_L_TX(0:3,NOBSMAX,NTMAX))
         ALLOCATE (DOBS_SUM_L_TX(0:3,NOBSMAX,NTMAX))
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
C         set the chanel number for the LHS wave functions
C
      IFIL_RHSA = IFILCBWF
C
      CALL SET_IFIL_LHS(IFIL_RHSA,IFIL_LHSB)
C
C=======================================================================
C
C-----------------------------------------------------------------------
C         read l-projected matrix elements MEZZL from IFILMEZZL
C-----------------------------------------------------------------------
      IF ( FULLPOT ) THEN
C
         IF ( .NOT.L_PROJECTED_ME )
     &         CALL STOP_MESSAGE(ROUTINE,'L_PROJECTED_ME = .FALSE.')
C
         ALLOCATE (MIX_DL(NMEMAX))
         ALLOCATE (MIX_IL(NMEMAX))
C
         ALLOCATE (MEZJL(NKMMAX,NKMMAX,NLMAX,NMEMAX))
         ALLOCATE (MEZZL(NKMMAX,NKMMAX,NLMAX,NMEMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEZZL')
C
         REWIND IFILMEZZL
C
C--------------------------------------------------------- position file
         DO ITP = 1,ITBOT - 1
C
            READ (IFILMEZZL) IDUMMY,MEZZL,MEZJL
C
            IF ( IDUMMY.EQ.ITBOT-1 ) EXIT
            IF ( IDUMMY.EQ.ITBOT .AND. ITP.EQ.1 ) THEN
               REWIND IFILMEZZL
               EXIT
            END IF
            IF ( IDUMMY.GE.ITBOT ) CALL STOP_MESSAGE(ROUTINE,
     &           'reading IFILMEZZL: IDUMMY >= ITBOT')
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C      calculate matrix elements  MZBZA and MZBJA in the LOCAL frame
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
      IF ( THERMAL_VIBRA_FLUCT .AND. .NOT.WAVE_FUNCTIONS_AVAILABLE )
     &     CALL STOP_MESSAGE(ROUTINE,
     &        'THERMAL_VIBRA_FLUCT .and. .NOT. WAVE_FUNCTIONS_AVAILABLE'
     &        )
C
      IF ( WAVE_FUNCTIONS_AVAILABLE .AND. 
     &     (THERMAL_VIBRA_FLUCT .OR. (MOMENTS_ROTATED .AND. IREL.EQ.3))
     &     ) THEN
C
         ALLOCATE (ZGRA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFRA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JGRA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFRA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZGLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JGLA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFLA(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGLA')
C
         ALLOCATE (MZBJA(NKMMAX,NKMMAX,3,NOBS))
         ALLOCATE (MZBZA(NKMMAX,NKMMAX,3,NOBS),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZBZA')
C
      END IF
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C
      IF ( IBND.LT.IDOS ) CALL STOP_MESSAGE(ROUTINE,'IBND < IDOS')
      IF ( NME.GT.NMEMAX ) CALL STOP_MESSAGE(ROUTINE,'NME > NMEMAX')
      IF ( IREL.GT.1 ) THEN
         WSPIN = 1D0
      ELSE
         WSPIN = 2D0
      END IF
C
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
      IF ( IECURR.EQ.1 .AND. IEPATH.EQ.1 ) THEN
C
         IF ( ICORE.EQ.0 ) THEN
            BCOR(ITBOT:ITTOP) = 0.0D0
            BCORS(ITBOT:ITTOP) = 0.0D0
         END IF
C
      END IF
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
C
      IF ( .NOT.SPLITSS ) THEN
         STR20 = 'CRYSTAL TERMS       '
         STR9 = 'crystal  '
      ELSE IF ( IEPATH.EQ.1 ) THEN
         STR20 = 'BACKSCATTERING TERMS'
         STR9 = 'backscat.'
      ELSE
         STR20 = 'SINGLE SITE TERMS   '
         STR9 = 'sgl. site'
      END IF
C
      CTOTALDOS = 0.0D0
      TOTNOS = 0.0D0
      TOTDOS = 0.0D0
      MUESPN = 0.0D0
      MUEORB = 0.0D0
      E_BAND = 0.0D0
C
      DOBS_LTX(:,:,:,:) = C0
      DOBS_TX(:,:,:) = C0
      DOBS_SUM_L_TX(:,:,:) = C0
C
      OBS_SUM_L_TX(:,:,:) = 0.0D0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT_A:DO IT = ITBOT,ITTOP
         N = NKM
         M = NKMMAX
C
         IM = IMT(IT)
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C      calculate matrix elements  MZBZA and MZBJA in the LOCAL frame
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
         IF ( WAVE_FUNCTIONS_AVAILABLE .AND. 
     &        (THERMAL_VIBRA_FLUCT .OR. (MOMENTS_ROTATED .AND. 
     &        IREL.EQ.3)) ) THEN
C
C=======================================================================
C                   read valence band wave functions
C                   NCPLWF  =   NCPLWF_RA =   NCPLWF_LB
C                 IKMCPLWF  = IKMCPLWF_RA = IKMCPLWF_LB
C=======================================================================
C
C ----------------------------------------- read in wavefunctions for RA
C
            CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGRA,ZFRA,JGRA,JFRA,
     &                           IRTOP,NCPLWF,IKMCPLWF)
C
            IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
               ZGLA(:,:,:) = ZGRA(:,:,:)
               ZFLA(:,:,:) = ZFRA(:,:,:)
               JGLA(:,:,:) = JGRA(:,:,:)
               JFLA(:,:,:) = JFRA(:,:,:)
C
            ELSE
C ----------------------------------------- read in wavefunctions for LA
C
               CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGLA,ZFLA,JGLA,JFLA,
     &                              IRTOP,NCPLWF,IKMCPLWF)
C
            END IF
C=======================================================================
C
C--------------------------------------------------------------IREL = 3
C                                                              LHS = RHS
            IF ( LHS_SOL_EQ_RHS_SOL ) THEN
C
               CALL CALC_OBS_ME(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA,MZBJA,NOBS,'S>C',
     &                          CHECK_OBS_ME)
C
C-------------------------------------------------------------- IREL = 3
C                                                             LHS != RHS
            ELSE
C
               CALL CALC_OBS_ME(IT,ZGLA,ZFLA,JGLA,JFLA,NCPLWF,IKMCPLWF,
     &                          ZGRA,ZFRA,JGRA,JFRA,NCPLWF,IKMCPLWF,
     &                          MEZZ,MEZJ,MZBZA,MZBJA,NOBS,'S>C',
     &                          CHECK_OBS_ME)
C
            END IF
C
         END IF
CRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTRTR
C
C----------------------------------------------------------- TEMPERATURE
         IF ( .NOT.THERMAL_VIBRA_FLUCT ) THEN
C
C--------------------- site diagonal multiple scattering matrix G or TAU
C
            CALL TAUGFCONV(MSST(1,1,IT),TAUT(1,1,IT),TMAT)
C
C--- provide site diagonal multiple scattering matrix in the LOCAL frame
C
            IF ( KMROT.EQ.0 ) THEN
C
               TMAT_LOC(:,:) = TMAT(:,:)
C
            ELSE
C
               IQ = IQAT(1,IT)
C
               CALL ROTATE(TMAT,'G->L',TMAT_LOC,N,DROTQ(1,1,IQ),M)
C
            END IF
C
         ELSE
C---------------------------------------------------- finite temperature
C
            CALL THERMAL_PROPERTIES(IT,IECURR,IEPATH,WE,TSST,MSST,MSSQ,
     &                              TAUQ,TMAT,MZBZA,MZBJA,DOBS_TX_GLO,
     &                              OBS_TX_GLO,NOBS)
C
            CALL TAUGFCONV(MSST(1,1,IT),TMAT,TMAT_LOC)
C
         END IF
C----------------------------------------------------------- TEMPERATURE
C
C--- provide site diagonal multiple scattering matrix in the LOCAL frame
C
         IF ( KMROT.EQ.0 ) THEN
C
            TMAT_LOC(:,:) = TMAT(:,:)
C
         ELSE
C
            IQ = IQAT(1,IT)
C
            CALL ROTATE(TMAT,'G->L',TMAT_LOC,N,DROTQ(1,1,IQ),M)
C
         END IF
C
         LMAX = NLT(IT) - 1
C
C-----------------------------------------------------------------------
C         read l-projected matrix elements MEZZL from IFILMEZZL
C-----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            READ (IFILMEZZL) IDUMMY,MEZZL,MEZJL
C
            IF ( IDUMMY.NE.IT )
     &            CALL STOP_MESSAGE(ROUTINE,'reading IFILMEZZL')
         END IF
C-----------------------------------------------------------------------
C
C
C MEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEME
         DO IME = 1,NME
C
            IF ( IME.NE.IBND ) THEN
C
               CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ(1,1,IT,IME),M,
     &                    TMAT_LOC,M,C0,WKM1,M)
C
               WKM1(:,:) = WKM1(:,:) - CPRE*MEZJ(:,:,IT,IME)
C
C---------------------------- compress to first (l,ml)-block if IREL = 2
C
               IF ( IREL.EQ.2 ) THEN
                  DO J = 1,NLM
                     DO I = 1,NLM
                        WKM1(I,J) = WKM1(I,J) + WKM1(NLM+I,NLM+J)
                        WKM1(NLM+I,NLM+J) = C0
                     END DO
                  END DO
               END IF
C
               IF ( IREL.LE.2 ) THEN
                  I2 = NLT(IT)*NLT(IT)
               ELSE
                  I2 = 2*NLT(IT)*NLT(IT)
               END IF
C
               DO I = 1,I2
                  DOBS_TX(0,IME,IT) = DOBS_TX(0,IME,IT) + WKM1(I,I)
               END DO
C
            ELSE
C
               DOBS_TX(0,IBND,IT) = WSPIN*DOBS_TX(0,IDOS,IT)*ERYD
C
            END IF
C
            OBS_TX(0,IME,IT) = OBS_TX(0,IME,IT)
     &                         + DIMAG(WE*DOBS_TX(0,IME,IT))
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO L = 0,LMAX
               IL = L + 1
C
C---------------------------------------- ASA and FULLPOT and IME = IBND
               IF ( IME.EQ.IBND ) THEN
                  DOBS_LTX(0,IBND,IL,IT) = WSPIN*DOBS_LTX(0,IDOS,IL,IT)
     &               *ERYD
C
C----------------------------------------------- FULLPOT and IME != IBND
               ELSE IF ( FULLPOT ) THEN
C
                  CALL ZGEMM('N','N',N,N,N,CPRE,MEZZL(1,1,IL,IME),M,
     &                       TMAT_LOC,M,C0,WKM1,M)
C
                  DO J = 1,N
                     CALL ZAXPY(N,-CPRE,MEZJL(1,J,IL,IME),1,WKM1(1,J),1)
                  END DO
C
C---------------------------- compress to first (l,ml)-block if IREL = 2
C
                  IF ( IREL.EQ.2 ) THEN
                     DO J = 1,NLM
                        DO I = 1,NLM
                           WKM1(I,J) = WKM1(I,J) + WKM1(NLM+I,NLM+J)
                           WKM1(NLM+I,NLM+J) = C0
                        END DO
                     END DO
                  END IF
C
                  IF ( IREL.LE.2 ) THEN
                     I2 = NLT(IT)*NLT(IT)
                  ELSE
                     I2 = 2*NLT(IT)*NLT(IT)
                  END IF
C
                  DOBS_LTX(0,IME,IL,IT) = 0.0D0
                  DO I = 1,I2
                     DOBS_LTX(0,IME,IL,IT) = DOBS_LTX(0,IME,IL,IT)
     &                  + WKM1(I,I)
                  END DO
C
               ELSE
C-------------------------------------------------- ASA and IME != IBND
C
                  IF ( IREL.LE.2 ) THEN
                     I0 = L*L
                     DO MX = 1,2*L + 1
                        DOBS_LMX(0,IME,IL,MX) = WKM1(I0+MX,I0+MX)
                     END DO
                  ELSE
                     I0 = 2*L*L + L + L
                     MXTOP = 2*L + 2
                     DO MX = 1,MXTOP
                        DOBS_LMX(0,IME,IL,MX) = WKM1(I0+MX,I0+MX)
                     END DO
                     I0 = 2*(L-1)*L + L + L
                     DO MX = 2,MXTOP - 1
                        DOBS_LMX(0,IME,IL,MX) = DOBS_LMX(0,IME,IL,MX)
     &                     + WKM1(I0+MX-1,I0+MX-1)
                     END DO
                  END IF
C
                  IF ( IREL.GT.2 ) THEN
                     MXTOP = 2*L + 2
                  ELSE
                     MXTOP = 2*L + 1
                  END IF
C
                  DO MX = 1,MXTOP
                     DOBS_LTX(0,IME,IL,IT) = DOBS_LTX(0,IME,IL,IT)
     &                  + DOBS_LMX(0,IME,IL,MX)
C
                  END DO
C
               END IF
C-------------------------------------------------- setting DOBS_LTX END
C
               OBS_LTX(0,IME,IL,IT) = OBS_LTX(0,IME,IL,IT)
     &                                + DIMAG(WE*DOBS_LTX(0,IME,IL,IT))
C
               DOBS_SUM_L_TX(:,IME,IT) = DOBS_SUM_L_TX(:,IME,IT)
     &            + DOBS_LTX(:,IME,IL,IT)
C
               OBS_SUM_L_TX(:,IME,IT) = OBS_SUM_L_TX(:,IME,IT)
     &                                  + OBS_LTX(:,IME,IL,IT)
C
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
               IF ( ((IPRINT.GT.0) .OR. (IECURR.EQ.NETAB(IEPATH))) .AND. 
     &              IME.EQ.NME ) THEN
C
                  IF ( FULLPOT ) THEN
                     TXTMJ = '          '
                  ELSE
                     TXTMJ = ' SUM(MJ)  '
                  END IF
                  WRITE (6,99006) IECURR,ERYD,L,IT,TXT_T(IT),STR20,
     &                            OBS_LTX(0,IDOS,IL,IT),
     &                            OBS_LTX(0,ISMT,IL,IT),
     &                            OBS_LTX(0,IOMT,IL,IT),
     &                            (OBS_LTX(0,IHFF,IL,IT)*1D-3),TXTMJ,
     &                            DOBS_LTX(0,IDOS,IL,IT),
     &                            DOBS_LTX(0,ISMT,IL,IT),
     &                            DOBS_LTX(0,IOMT,IL,IT),
     &                            (DOBS_LTX(0,IHFF,IL,IT)*1D-6)
C
                  IF ( .NOT.FULLPOT ) THEN
                     IF ( IREL.GT.2 ) THEN
                        MXTOP = 2*L + 2
C
                        WRITE (6,99007) ((-MXTOP-1+2*MX),DOBS_LMX(0,IDOS
     &                                  ,IL,MX),DOBS_LMX(0,ISMT,IL,MX),
     &                                  DOBS_LMX(0,IOMT,IL,MX),
     &                                  (DOBS_LMX(0,IHFF,IL,MX)*1D-6),
     &                                  MX=1,MXTOP)
C
                     ELSE
                        MXTOP = 2*L + 1
C
                        WRITE (6,99019) ((-L-1+MX),DOBS_LMX(0,IDOS,IL,MX
     &                                  ),DOBS_LMX(0,ISMT,IL,MX),
     &                                  DOBS_LMX(0,IOMT,IL,MX),
     &                                  (DOBS_LMX(0,IHFF,IL,MX)*1D-6),
     &                                  MX=1,MXTOP)
C
                     END IF
                  END IF
C
               END IF
C wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
C
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         END DO
C MEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEMEME
C
         IF ( NWRLOG.NE.0 ) WRITE (6,99008) IECURR,ERYD,IT,TXT_T(IT),
     &                             STR9,OBS_TX(0,IDOS,IT),
     &                             OBS_TX(0,ISMT,IT),OBS_TX(0,IOMT,IT),
     &                             (OBS_TX(0,IHFF,IT)*1D-3),STR9,
     &                             DIMAG(DOBS_TX(0,IDOS,IT)),
     &                             DIMAG(DOBS_TX(0,ISMT,IT)),
     &                             DIMAG(DOBS_TX(0,IOMT,IT)),
     &                             (DIMAG(DOBS_TX(0,IHFF,IT))*1D-3)
C
         CTOTALDOS = CTOTALDOS + DOBS_TX(0,IDOS,IT)*CONC(IT)*NAT(IT)
C
         TOTDOS = TOTDOS + DIMAG(DOBS_TX(0,IDOS,IT))*CONC(IT)*NAT(IT)
         TOTNOS = TOTNOS + OBS_TX(0,IDOS,IT)*CONC(IT)*NAT(IT)
         MUESPN = MUESPN + OBS_TX(0,ISMT,IT)*CONC(IT)*NAT(IT)
         MUEORB = MUEORB + OBS_TX(0,IOMT,IT)*CONC(IT)*NAT(IT)
C
         IF ( NWRLOG.NE.0 ) THEN
            IF ( IT.LT.ITTOP ) THEN
               WRITE (6,'(1X,79(''-''))')
            ELSE IF ( (IPRINT.GT.0) .OR. (IECURR.EQ.NETAB(IEPATH)) )
     &                THEN
               WRITE (6,99014) TOTDOS,TOTNOS,MUESPN,MUEORB
            ELSE
               WRITE (6,'('' '',79(''=''))')
            END IF
         END IF
C
      END DO LOOP_IT_A
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( IECURR.LT.NEMAX ) CTOTDOS(IECURR) = CTOTALDOS
C
C------------------------------------- store backscattering contribution
C
      IF ( SPLITSS .AND. (IEPATH.EQ.1) .AND. (IECURR.EQ.NETAB(IEPATH)) )
     &     THEN
C
         TOTDOS_BS_EF = TOTDOS
C
         DOBS_BS_EF_TX(:,:,ITBOT:ITTOP) = DOBS_TX(:,:,ITBOT:ITTOP)
         DOBS_BS_EF_LTX(:,:,:,ITBOT:ITTOP) = DOBS_LTX(:,:,:,ITBOT:ITTOP)
C
      END IF
C
C
C ======================================================================
C ======================================================================
C     EF reached for path NEPATH  ---  correction at the Fermi energy
C ======================================================================
C ======================================================================
C
C
      IF ( (IGRID(IEPATH).LT.6 .OR. IGRID(IEPATH).EQ.8) .AND. 
     &     INTEFERMI .AND. (IECURR.EQ.NETAB(IEPATH)) .AND. 
     &     (IEPATH.EQ.NEPATH) ) THEN
C
C---------------- add the backscattering contribution that has been kept
         IF ( SPLITSS ) THEN
C
            TOTDOS = TOTDOS + TOTDOS_BS_EF
C
            DOBS_TX(:,:,ITBOT:ITTOP) = DOBS_TX(:,:,ITBOT:ITTOP)
     &                                 + DOBS_BS_EF_TX(:,:,ITBOT:ITTOP)
C
            DOBS_LTX(:,:,:,ITBOT:ITTOP) = DOBS_LTX(:,:,:,ITBOT:ITTOP)
     &         + DOBS_BS_EF_LTX(:,:,:,ITBOT:ITTOP)
C
         END IF
C
         DQ = TOTNOS - NVALTOT/WSPIN
C
         IF ( (SUB_SYSTEM(1:6).EQ.'I-ZONE' .AND. SYSTEM_TYPE(1:3)
     &        .NE.'VIV') .OR. SUB_SYSTEM(1:6).EQ.'R-BULK' .OR. 
     &        .NOT.UPDATE_EFERMI ) THEN
C
            SHFTEF = 0.0D0
C
         ELSE
C
            SHFTEF = -DQ/TOTDOS
C
         END IF
C
         IF ( .NOT.FULLPOT ) EFERMI = EFERMI + SHFTEF !!!!!!!!!!!!!!!!!!
         IEF = IECURR + 1
C
         IF ( NWRLOG.NE.0 ) WRITE (6,'(/)')
         DO IWRLOG = 1,NWRLOG
            IF ( IWRLOG.EQ.1 ) IW = 6
            IF ( IWRLOG.EQ.2 ) IW = IFILLOG
            WRITE (IW,99016) SYSTEM(1:LSYSTEM)
            IF ( UPDATE_EFERMI ) THEN
               WRITE (IW,99017) DQ,SHFTEF,EFERMI
            ELSE
               WRITE (IW,99018) DQ,EFERMI
            END IF
            WRITE (IW,99001) NWEDGE,(IWEDGEROT(I),I=1,NWEDGE)
            IF ( IBZINT.EQ.1 ) WRITE (IW,99002) NPTMIN,NPTMAX
            IF ( IBZINT.EQ.2 ) WRITE (IW,99003) NKTAB
            IF ( IBZINT.EQ.3 ) WRITE (IW,99004) NKTAB
            WRITE (IW,99005) NCPAFAIL
C
         END DO
C
C=========================================================== TEMPERATURE
C             project m_z in local frame onto the z-axis of global frame
C                                   in case of thermal spin fluctuations
C
         IF ( THERMAL_VIBRA_FLUCT .AND. NFLUCT.GT.1 ) THEN
C
            DO IT = ITBOT,ITTOP
C
               IFT = (IT-1)*NFLUCT
               WSUM = 0D0
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
                  WANG_FT = COS(FTET_FT(IFT)*PI/180D0)
                  WSUM = WSUM + X_FT(IFT)*WANG_FT
               END DO
C
               OBS_TX(0,ISMT,IT) = OBS_TX(0,ISMT,IT)*WSUM
               OBS_TX(0,IOMT,IT) = OBS_TX(0,IOMT,IT)*WSUM
            END DO
C
            WRITE (111,'(15f12.6)') TEMP_LAT,OBS_TX(0,ISMT,IT)/WSUM,
     &                              OBS_TX(0,ISMT,IT),OBS_TX(0,IOMT,IT)
     &                              /WSUM,OBS_TX(0,IOMT,IT),WSUM
         END IF
C ======================================================================
C
         TOTNOS = 0.0D0
         MUESPN = 0.0D0
         MUEORB = 0.0D0
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         LOOP_IT_B:DO IT = ITBOT,ITTOP
C
            DO IME = 1,NME
C
               OBS_TX(0,IME,IT) = OBS_TX(0,IME,IT)
     &                            + SHFTEF*DIMAG(DOBS_TX(0,IME,IT))
C
               OBS_SUM_L_TX(0,IME,IT) = OBS_SUM_L_TX(0,IME,IT)
     &                                  + SHFTEF*DIMAG
     &                                  (DOBS_SUM_L_TX(0,IME,IT))
C
               DO IL = 1,NLT(IT)
                  OBS_LTX(0,IME,IL,IT) = OBS_LTX(0,IME,IL,IT)
     &               + SHFTEF*DIMAG(DOBS_LTX(0,IME,IL,IT))
               END DO
C
            END DO
C
            TOTNOS = TOTNOS + OBS_TX(0,IDOS,IT)*CONC(IT)*NAT(IT)
            MUESPN = MUESPN + OBS_TX(0,ISMT,IT)*CONC(IT)*NAT(IT)
            MUEORB = MUEORB + OBS_TX(0,IOMT,IT)*CONC(IT)*NAT(IT)
            E_BAND = E_BAND + OBS_TX(0,IBND,IT)*CONC(IT)*NAT(IT)
C
            DO IWRLOG = 1,NWRLOG
               IF ( IWRLOG.EQ.1 ) IW = 6
               IF ( IWRLOG.EQ.2 ) IW = IFILLOG
               WRITE (IW,99008) IEF,EFERMI,0.0D0,IT,TXT_T(IT)
C
               BDUM(1) = BCORS(IT)*1D-3
               BDUM(2) = (BCOR(IT)-BCORS(IT))*1D-3
               BDUM(3) = BCOR(IT)*1D-3
C
               WRITE (IW,99009) (DIMAG(DOBS_LTX(0,IDOS,IL,IT)),OBS_LTX(0
     &                          ,IDOS,IL,IT),
     &                          DIMAG(DOBS_LTX(0,ISMT,IL,IT)),
     &                          OBS_LTX(0,ISMT,IL,IT),
     &                          DIMAG(DOBS_LTX(0,IOMT,IL,IT)),
     &                          OBS_LTX(0,IOMT,IL,IT),
     &                          OBS_LTX(0,IHFF,IL,IT)*1D-3,BDUM(IL),
     &                          IL=1,MIN(3,NLT(IT)))
C
               IF ( NLT(IT).GT.3 ) WRITE (IW,99010)
     &              (DIMAG(DOBS_LTX(0,IDOS,IL,IT)),OBS_LTX(0,IDOS,IL,IT)
     &              ,DIMAG(DOBS_LTX(0,ISMT,IL,IT)),OBS_LTX(0,ISMT,IL,IT)
     &              ,DIMAG(DOBS_LTX(0,IOMT,IL,IT)),OBS_LTX(0,IOMT,IL,IT)
     &              ,OBS_LTX(0,IHFF,IL,IT)*1D-3,IL=4,NLT(IT))
C
C-----------------------------------------------------------------------
C                           mixed contributions
C-----------------------------------------------------------------------
               IF ( FULLPOT ) THEN
C
                  DO IME = 1,NME
                     MIX_DL(IME) = DIMAG(DOBS_TX(0,IME,IT)-DOBS_SUM_L_TX
     &                             (0,IME,IT))
                     MIX_IL(IME) = OBS_TX(0,IME,IT)
     &                             - OBS_SUM_L_TX(0,IME,IT)
                  END DO
C
                  WRITE (IW,99011) 'mix',MIX_DL(IDOS),MIX_IL(IDOS),
     &                             MIX_DL(ISMT),MIX_IL(ISMT),
     &                             MIX_DL(IOMT),MIX_IL(IOMT),
     &                             MIX_IL(IHFF)*1D-3
C
               END IF
C-----------------------------------------------------------------------
C
               WRITE (IW,99011) 'sum',DIMAG(DOBS_TX(0,IDOS,IT)),
     &                          OBS_TX(0,IDOS,IT),
     &                          DIMAG(DOBS_TX(0,ISMT,IT)),
     &                          OBS_TX(0,ISMT,IT),
     &                          DIMAG(DOBS_TX(0,IOMT,IT)),
     &                          OBS_TX(0,IOMT,IT),
     &                          (OBS_TX(0,IHFF,IT)*1D-3),
     &                          ((OBS_TX(0,IHFF,IT)+BCOR(IT))*1D-3)
C
               WRITE (IW,99012) (DIMAG(DOBS_LTX(0,ISDM,IL,IT)),OBS_LTX(0
     &                          ,ISDM,IL,IT),IL=1,NLT(IT))
               WRITE (IW,99013) DIMAG(DOBS_TX(0,ISDM,IT)),
     &                          OBS_TX(0,ISDM,IT)
C
               IF ( IT.LT.ITTOP ) THEN
                  WRITE (IW,'(1X,79(''-''))')
               ELSE
                  WRITE (IW,99015) TOTDOS,TOTNOS,MUESPN,MUEORB,E_BAND
               END IF
            END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C NOTE: HFF taken out for the moment as it is too sensitive
            IF ( WRBUILDBOT .AND. MPI_ID.EQ.0 )
     &           WRITE (IFILBUILDBOT,99021) ROUTINE(1:LEN_TRIM(ROUTINE))
     &           ,IT,DIMAG(DOBS_TX(0,IDOS,IT)),OBS_TX(0,IDOS,IT),
     &           OBS_TX(0,ISMT,IT),OBS_TX(0,IOMT,IT),OBS_TX(0,IBND,IT)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         END DO LOOP_IT_B
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END IF
C============================================ EF reached for path NEPATH
C=======================================================================
C
      IF ( ALLOCATED(ZGRA) ) DEALLOCATE (ZGRA,ZFRA,JGRA,JFRA)
      IF ( ALLOCATED(JFLA) ) DEALLOCATE (JGLA,JFLA,ZFLA,ZGLA)
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C                  store DOS for later write out
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      IF ( PROGNAME(4:6).EQ.'GEN' ) DOBS_LTEX(:,:,:,:,IECURR)
     &     = DIMAG(DOBS_LTX(:,:,:,:))
C
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT .AND. MPI_ID.EQ.0 .AND. IECURR.LE.3 ) THEN
         IF ( IREL.GT.1 ) THEN
            SPF = 0.5D0
         ELSE
            SPF = 1.0D0
         END IF
         WRITE (IFILBUILDBOT,99020) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'DOS       ',IECURR,DREAL(ERYD),
     &                              (((SPF*DIMAG(DOBS_LTX(0,IDOS,IL,IT)
     &                              +DOBS_LTX(0,ISMT,IL,IT))),IL=1,
     &                              NLT(IT)),
     &                              ((SPF*DIMAG(DOBS_LTX(0,IDOS,IL,IT)
     &                              -DOBS_LTX(0,ISMT,IL,IT))),IL=1,
     &                              NLT(IT)),IT=ITBOT,ITTOP)
         IF ( FULLPOT ) RETURN
         WRITE (IFILBUILDBOT,99020) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'DOS-kappa ',IECURR,DREAL(ERYD),
     &                              ((DIMAG(DOBS_LTX(0,IDOS,IL,IT)),
     &                              IL=1,NLT(IT)),
     &                              DIMAG(DOBS_LTX(0,IDOS,1,IT)),
     &                              (0.5D0*DIMAG(DOBS_LTX(0,IDOS,IL,IT)
     &                              -DOBS_LTX(0,IKDS,IL,IT)),
     &                              0.5D0*DIMAG(DOBS_LTX(0,IDOS,IL,IT)
     &                              +DOBS_LTX(0,IKDS,IL,IT)),IL=2,
     &                              NLT(IT)),IT=ITBOT,ITTOP)
         WRITE (IFILBUILDBOT,99020) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                              'T_z - term',IECURR,DREAL(ERYD),
     &                              ((DIMAG(DOBS_LTX(0,ISMT,IL,IT)),
     &                              IL=1,NLT(IT)),
     &                              (DIMAG(DOBS_LTX(0,IOMT,IL,IT)),IL=1,
     &                              NLT(IT)),
     &                              (DIMAG(DOBS_LTX(0,ISDM,IL,IT)),IL=1,
     &                              NLT(IT)),IT=ITBOT,ITTOP)
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
99001 FORMAT (' parameters for k-space integration:','  NWEDGE:',I3,
     &        '  wedges:',25I3)
99002 FORMAT (' WEYL - point sampling  NPT:',I5,' -',I5,
     &        '  according to Im(E)')
99003 FORMAT (' regular k-mesh point sampling  NKTAB:',I5)
99004 FORMAT (' tetrahedron method used with mesh  NF:',I5)
99005 FORMAT (' number of energies CPA failed ',I5)
99006 FORMAT (/,I4,' E=',2F7.4,3X,'L=',I2,3X,'IT=',I4,2X,A,2X,A20,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE)  ',8X,F8.3,10X,F8.3,10X,
     &        F8.3,10X,F8.1,/,A,2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1)
99007 FORMAT (' MJ= ',I2,'/2 ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1)
99008 FORMAT (/,I4,' E=',2F7.4,10X,'IT=',I4,2X,A,:,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE) ',A9,F8.3,4X,F14.3,4X,F14.3,
     &        4X,F14.1,/,' TOTAL   ',A9,F8.3,4X,F14.3,4X,F14.3,4X,F14.1)
99009 FORMAT ('         DOS      NOS     P_spin   m_spin',
     &        '    P_orb    m_orb    B_val      B_core',/,'  s ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' s  ',F8.2,:,/,'  p ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' ns ',F8.2,:,/,'  d ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' cor',F8.2)
99010 FORMAT ('  f ',2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,:,/,'  g ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,:,/,'  h ',2F9.4,F10.4,F9.4,
     &        F10.5,F9.5,F8.2,:,/,'  i ',2F9.4,F10.4,F9.4,F10.5,F9.5,
     &        F8.2)
99011 FORMAT (1X,A3,2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,:,' v+c',F8.2)
99012 FORMAT ('                          P(T_z)   T_z   ',/,'  s ',18X,
     &        F10.4,F9.4,:,/,'  p ',18X,F10.4,F9.4,:,/,'  d ',18X,F10.4,
     &        F9.4,:,/,'  f ',18X,F10.4,F9.4,:,/,'  g ',18X,F10.4,F9.4,
     &        :,/,'  h ',18X,F10.4,F9.4,:,/,'  i ',18X,F10.4,F9.4)
99013 FORMAT (' sum',18X,F10.4,F9.4)
99014 FORMAT (1X,79('-'),/,' TOT',2F9.4,10X,F9.4,10X,F9.5,/,1X,79('='))
99015 FORMAT (1X,79('-'),/,' TOT',2F9.4,10X,F9.4,10X,F9.5,/,' E_band',
     &        F19.8,' [Ry]',/,' ',79('='))
99016 FORMAT ((' ',79('*'),/),/,' SPRKKR-run for: ',A)
99017 FORMAT (/,' results extrapolated to corrected Fermi energy:',/,
     &        ' CHARGE MISFIT     ',F12.8,' els.',/,
     &        ' E_F CORRECTION    ',F12.8,/,' NEW FERMI ENERGY  ',F12.8,
     &        /)
99018 FORMAT (/,' results integrated up to fixed Fermi energy:',/,
     &        ' CHARGE MISFIT     ',F12.8,' els.',/,
     &        ' FERMI ENERGY      ',F12.8,/)
99019 FORMAT (' ML= ',I2,'   ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1)
99020 FORMAT ('# BUILDBOT: ',A,':  ',A,'  for IE, E =',I5,F12.5,/,
     &        (1PE22.14))
99021 FORMAT ('# BUILDBOT: ',A,':  DOS NOS SMT OMT BND for IT =',I5,/,
     &        (1PE22.14))
      END
