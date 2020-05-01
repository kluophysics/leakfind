C*==cludrive.f    processed by SPAG 6.70Rc at 08:55 on  8 Mar 2017
      SUBROUTINE CLUDRIVE(DOSREP,DOSCORSYS)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine for embedded cluster calculations                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C_AU
      USE MOD_TABLES,ONLY:TAB_CHSYM
      USE MOD_SCF,ONLY:SCFSTATUS_CLU
      USE MOD_RMESH,ONLY:RWS,RMT,FULLPOT,NM,NMHOST,NMCLU
      USE MOD_CALCMODE,ONLY:KKRMODE,PROGNAME,CLUTYPE,KMROT,TXTKMROT,
     &    TASK,UPDATE_EFERMI,RELAX_CLU,U_CALCULATION,
     &    IT_CENTER_U_CALCULATION,SOLVER,IREL,MOMENTS_ROTATED
      USE MOD_ENERGY,ONLY:NETAB
      USE MOD_TYPES,ONLY:NTCLU,NTHOST,NTMAX,ITBOT,ITTOP,SOCTL,CTL,NLT,
     &    NLIN_T,NKM_T,NAT,NVALTOT,CONC,NVALT,TXT_T,LTXT_T,IMT,Z,
     &    NA_TCLU,Z_TCLU,X_TCLU,NTCLUMAX,NCORT,NSEMCORSHLT,TXT_TCLU,
     &    NCOR_TCLU,NVAL_TCLU,NSEMCORSHL_TCLU,NT,LOPT
      USE MOD_SITES,ONLY:NQMAX,NQCLU,IQ_QCLU,IQBOT,IQTOP,NQHOST,ITOQ,
     &    NO_QCLU,NOQ,IQBOT_CLU,IQTOP_CLU,QMGAM,QMTET,QMPHI,QMVEC,ICPA,
     &    IMQ,IQAT,NQ,IT_OQCLU,ITCLU_OQCLU,IQCLU_ATCLU,NQCLUMAX,
     &    VLMMAD_HOST,N5VEC_QCLU,DQBAS_QCLU,QBAS_QCLU,SHIFT_CLU,
     &    QBAS0_QCLU,IQBOT_HOST,IQTOP_HOST,FORCE_QCLU,MAGROT_Q
      USE MOD_ANGMOM,ONLY:IXM0_QCLU,NXM_QCLU,NKMQ,NLMQ,NKKR_CLU,NLMAX,
     &    NLQ,NLINQ,NXM_Q,NL
      USE MOD_STR,ONLY:RGNT,IRGNT,NRGNT,IJQ,NIJQ,GGJLRS,PWP,SRREL,NRREL,
     &    IRREL,G1,G2,G3,R1,R2,R3,QQPZ,QQPY,QQPX,HP,EXPGNQ,INDR,SMAX,
     &    QQMLRS,CILMAT,D1TERM3,IILERS,DLLMMKE,CHP,SHP,WK,PWEX2K1,
     &    PWEX2K2,PWEX2K3,PWEXIK1,PWEXIK2,PWEXIK3
      USE MOD_FILES,ONLY:DOSFIL,DATSET,LDATSET,POTFIL,LPOTFIL,IFILDOS,
     &    IFILCBWF,IFILCBWF_SPH,IFILCBWF_LHS,RECLNGWF,RECLNGWF_SPH,
     &    FOUND_SECTION,FOUND_REAL_ARRAY,N_FOUND
      USE MOD_SYMMETRY,ONLY:ISYMGENQ,IQREPQ,IQEQ,NQEQ,NQCLU_EQCLU,
     &    IQCLU_EQCLU
      USE MOD_LATTICE,ONLY:SWS,SYSTEM_TYPE
      USE MOD_MPI,ONLY:MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUDRIVE')
C
C Dummy arguments
C
      CHARACTER*3 DOSCORSYS,DOSREP
C
C Local variables
C
      REAL*8 AUXL(:),NVALTOT_CLU,SCL
      REAL*8 DNRM2
      LOGICAL FULLPOT5,INIT_RELAXATION,MANPULC,MANPULSOC,RUNELOOP
      INTEGER IA,IEQ,IFLAG,IL,IM,IMCLU,IO,IQ,IQCLU,IQCLU_CENTER,IQHOST,
     &        IQ_ATCLU(:,:),IT,ITCLU,IT_CENTER,KBZI,KMOL
      CHARACTER*2 SOCII(-2:-1)
      CHARACTER*10 STR10
      INTEGER TABNVAL
C
C*** End of declarations rewritten by SPAG
C
C ======================================================================
      DATA INIT_RELAXATION/.FALSE./
      DATA SOCII/'xy','zz'/
C ======================================================================
C
      ALLOCATABLE AUXL,IQ_ATCLU
C
      CALL TRACK_INFO(ROUTINE)
C
      WRITE (6,99017)
C
      UPDATE_EFERMI = .FALSE.
C
      CALL CLUPOTRD(2)
C
C=======================================================================
C shift atomic positions along direction of forces in case of relaxation
C=======================================================================
C
      IF ( RELAX_CLU .AND. ABS(SHIFT_CLU).GT.1D-10 ) THEN
C
         INIT_RELAXATION = .TRUE.
C
         DO IQCLU = 1,NQCLU
C
            IF ( SHIFT_CLU.LT.0D0 ) THEN
               SCL = DNRM2(3,FORCE_QCLU(1:3,IQCLU),1)
               IF ( ABS(SCL).GT.1D-10 ) THEN
                  SCL = SHIFT_CLU/SCL
               ELSE
                  SCL = 0D0
               END IF
C
               DQBAS_QCLU(1:3,IQCLU) = SCL*FORCE_QCLU(1:3,IQCLU)
            END IF
C
            DQBAS_QCLU(1:3,IQCLU) = SHIFT_CLU*QBAS0_QCLU(1:3,IQCLU)
C           dqbas_qclu(1:3,iqclu) = 0d0
C
            QBAS_QCLU(1:3,IQCLU) = QBAS0_QCLU(1:3,IQCLU)
     &                             + DQBAS_QCLU(1:3,IQCLU)
C
         END DO
C
      END IF
C
C=======================================================================
C              initialize all tables for cluster calculation
C=======================================================================
C
      ALLOCATE (IT_OQCLU(NTCLUMAX,NQCLUMAX))
C
C-----------------------------------------------------------------------
C                       site specific settings
C-----------------------------------------------------------------------
C
      IQBOT_HOST = IQBOT
      IQTOP_HOST = IQTOP
C
      ALLOCATE (IXM0_QCLU(NQCLU),NXM_QCLU(NQCLU))
C
      NAT(NTHOST+1:NTHOST+NTCLU) = 0
C
      DO IQCLU = 1,NQCLU
C
         IQHOST = IQ_QCLU(IQCLU)
C
C--------------------------------------- settings for cluster KKR matrix
C
         NXM_QCLU(IQCLU) = NXM_Q(IQHOST)
C
         IF ( IQCLU.EQ.1 ) THEN
            IXM0_QCLU(IQCLU) = 0
         ELSE
            IXM0_QCLU(IQCLU) = IXM0_QCLU(IQCLU-1) + NXM_QCLU(IQCLU-1)
         END IF
C
C-------- settings for site dependent angular momentum expansion and CPA
C
         IQ = NQHOST + IQCLU
C
         NLQ(IQ) = NLQ(IQHOST)
         NLMQ(IQ) = NLMQ(IQHOST)
         NKMQ(IQ) = NKMQ(IQHOST)
         NLINQ(IQ) = NLINQ(IQHOST)
         IMQ(IQ) = IMQ(IQHOST)
C
         NOQ(IQ) = NO_QCLU(IQCLU)
         IF ( NOQ(IQ).GT.1 ) THEN
            ICPA(IQ) = 1
         ELSE
            ICPA(IQ) = 0
         END IF
C
         NQEQ(IQ) = NQCLU_EQCLU(IQCLU)
         DO IEQ = 1,NQEQ(IQ)
            IQEQ(IEQ,IQ) = NQHOST + IQCLU_EQCLU(IEQ,IQCLU)
         END DO
C
         DO IO = 1,NO_QCLU(IQCLU)
            IT = NTHOST + ITCLU_OQCLU(IO,IQCLU)
            ITOQ(IO,IQ) = IT
            IT_OQCLU(IO,IQCLU) = IT
         END DO
C
      END DO
C-----------------------------------------------------------------------
C                       type specific settings
C-----------------------------------------------------------------------
C
      DO ITCLU = 1,NTCLU
C
         IT = NTHOST + ITCLU
C
         NAT(IT) = NA_TCLU(ITCLU)
         DO IA = 1,NA_TCLU(ITCLU)
            IQ = NQHOST + IQCLU_ATCLU(IA,ITCLU)
            IQAT(IA,IT) = IQ
         END DO
C
         IQ = IQAT(1,IT)
C
         NLT(IT) = NLQ(IQ)
         NKM_T(IT) = NKMQ(IQ)
         NLIN_T(IT) = NLINQ(IQ)
         IMT(IT) = IMQ(IQ)
C
         SOCTL(IT,1:NLMAX) = 1.0D0
         CTL(IT,1:NLMAX) = C_AU
C
         Z(IT) = Z_TCLU(ITCLU)
         CONC(IT) = X_TCLU(ITCLU)
C
         IF ( SCFSTATUS_CLU(1:5).EQ.'START' ) THEN
            TXT_T(IT) = TAB_CHSYM(Z(IT))
            LTXT_T(IT) = LEN_TRIM(TXT_T(IT))
            NVALT(IT) = TABNVAL(Z(IT))
            NCORT(IT) = Z(IT) - NVALT(IT)
            NSEMCORSHLT(IT) = 0
         END IF
      END DO
C
      NKKR_CLU = IXM0_QCLU(NQCLU) + NXM_QCLU(NQCLU)
C
C=======================================================================
C  specify whether the SOC should be manipulated l- and type-dependent
C    MANPULC:      via CTL     scale the speed of light
C    MANPULSOC:    via SOCTL   scale strength of SOC
C                              or select between xy- and zz-SOC term
C=======================================================================
      IF ( IREL.EQ.3 ) THEN
C
         MANPULC = .FALSE.
         MANPULSOC = .FALSE.
C
         ALLOCATE (AUXL(NLMAX))
C
C------------------------------ apply the same manipulation to ALL types
C
         CALL INPUT_FIND_SECTION('MODE',0)
C
         IF ( FOUND_SECTION ) THEN
C
            STR10 = 'C'
            CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,N_FOUND,NLMAX,0,
     &                                  9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
               DO ITCLU = 1,NTCLU
                  IT = NTHOST + ITCLU
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        CTL(IT,IL) = CTL(IT,IL)/AUXL(IL)**0.5D0
                        MANPULC = .TRUE.
                     END IF
                  END DO
               END DO
            END IF
C
            STR10 = 'SOC'
            CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,N_FOUND,NLMAX,0,
     &                                  9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
               DO ITCLU = 1,NTCLU
                  IT = NTHOST + ITCLU
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-I  '
                     END IF
                     IF ( AUXL(IL).LT.0D0 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-II '
                     END IF
                  END DO
               END DO
            END IF
C
C------------------ apply the manipulation to individually for the types
C
            DO ITCLU = 1,NTCLU
               IT = NTHOST + ITCLU
C
               STR10 = 'C'
               CALL STRING_ADD_N(STR10,IT)
               CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,N_FOUND,NLMAX,0,
     &            9999D0,0)
C
               IF ( FOUND_REAL_ARRAY ) THEN
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        CTL(IT,IL) = CTL(IT,IL)/AUXL(IL)**0.5D0
                        MANPULC = .TRUE.
                     END IF
                  END DO
               END IF
C
               STR10 = 'SOC'
               CALL STRING_ADD_N(STR10,IT)
C
               CALL SECTION_SET_REAL_ARRAY(STR10,AUXL,N_FOUND,NLMAX,0,
     &            9999D0,0)
C
               IF ( FOUND_REAL_ARRAY ) THEN
                  DO IL = 1,NLMAX
                     IF ( ABS(AUXL(IL)-1.0D0).GT.1D-6 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-I  '
                     END IF
                     IF ( AUXL(IL).LT.0D0 ) THEN
                        SOCTL(IT,IL) = AUXL(IL)
                        MANPULSOC = .TRUE.
                        SOLVER = 'BS-SOC-II '
                     END IF
                  END DO
               END IF
            END DO
C
         END IF
C
         IF ( MANPULC ) THEN
            WRITE (6,'(/,10X,''the speed of light is scaled with:'')')
            DO ITCLU = 1,NTCLU
               IT = NTHOST + ITCLU
               WRITE (6,99001) IT,(NL-1),((C_AU/CTL(IT,IL))**2,IL=1,NL)
            END DO
         END IF
         IF ( MANPULSOC ) THEN
            IF ( SOLVER.EQ.'BS-SOC-I  ' ) THEN
               WRITE (6,'(/,10X,''the SOC is scaled with:'')')
               DO ITCLU = 1,NTCLU
                  IT = NTHOST + ITCLU
                  WRITE (6,99001) IT,(NL-1),(SOCTL(IT,IL),IL=1,NL)
                  DO IL = 1,NL
                     IF ( SOCTL(IT,IL).LT.0D0 )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                  'all SOCTL should be >= 0 >>>> check input file'
     &                  )
                  END DO
               END DO
               SOLVER = 'BS-SOC    '
            ELSE
               WRITE (6,99002)
               DO ITCLU = 1,NTCLU
                  IT = NTHOST + ITCLU
                  WRITE (6,99003) IT,NL - 1,
     &                            (SOCII(NINT(SOCTL(IT,IL))),IL=1,NL)
                  DO IL = 1,NL
                     IF ( SOCTL(IT,IL).GT.-0.1D0 )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                   'all SOCTL should be < 0 >>>> check input file'
     &                   )
                  END DO
               END DO
               SOLVER = 'BS-SOC    '
            END IF
         END IF
C
      END IF
C
C=======================================================================
C       determine set of equivalent cluster site pairs within TAU_ij
C=======================================================================
C
      CALL TAUIJ_EQUIV_PAIRS
C
C=======================================================================
C                 open or create requested TAUIJ file
C                perform host calculation in later case
C=======================================================================
C
      NQ = NQHOST
      NT = NTHOST
C
      CALL CLUTAUIJ_DRIVE(DOSREP,DOSCORSYS)
C
C=======================================================================
C
      IF ( CLUTYPE.EQ.'embedded' ) THEN
         KKRMODE = 'EMBEDDED-CLUSTER'
         SYSTEM_TYPE = 'EMBEDDED-CLUSTER'
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'CLUTYPE <> embedded')
      END IF
C
C-----------------------------------------------------------------------
C in case of ASA calculations no additional radial meshes ar needed
C
      IF ( .NOT.FULLPOT ) NMCLU = 0
C
      NQ = NQHOST + NQCLU
      NT = NTHOST + NTCLU
      NM = NMHOST + NMCLU
C
      IQBOT_CLU = NQHOST + 1
      IQTOP_CLU = NQHOST + NQCLU
C
      IQBOT = NQHOST + 1
      IQTOP = NQHOST + NQCLU
C
      ITBOT = NTHOST + 1
      ITTOP = NTHOST + NTCLU
C
      NVALTOT_CLU = 0.0D0
      DO IT = NTHOST + 1,NTHOST + NTCLU
         NVALTOT_CLU = NVALTOT_CLU + CONC(IT)*NAT(IT)*NVALT(IT)
      END DO
      NVALTOT = NVALTOT_CLU
C
      CALL EXTEND_TXT_T
C
      DO IT = NTHOST + 1,NTHOST + NTCLU
         ITCLU = IT - NTHOST
C
         TXT_TCLU(ITCLU) = TXT_T(IT)
         Z_TCLU(ITCLU) = Z(IT)
         X_TCLU(ITCLU) = CONC(IT)
         NCOR_TCLU(ITCLU) = NCORT(IT)
         NVAL_TCLU(ITCLU) = NVALT(IT)
         NSEMCORSHL_TCLU(ITCLU) = NSEMCORSHLT(IT)
      END DO
C
C=======================================================================
C      calculate rotation matrices  MROTQ  and  DROTQ
C      connecting LOCAL and GLOBAL frame individually for
C      all sites IQ according to the flag  KMROT
C=======================================================================
C
      CALL INIT_MOD_SITES_ROTMAG
C
C-----------------------------------------------------------------------
C                prepare FULLPOT calculations
C-----------------------------------------------------------------------
C
      IF ( FULLPOT ) CALL CLUINITFULLPOT(INIT_RELAXATION)
C
C-----------------------------------------------------------------------
C
      WRITE (6,99005) NQCLU
C
      IF ( (UBOUND(MAGROT_Q,1).LT.(NQHOST+NQCLU)) .OR. 
     &     (UBOUND(QMTET,1).LT.(NQHOST+NQCLU)) )
     &     CALL STOP_MESSAGE(ROUTINE,'UBOUND(MAGROT_Q) < NQ')
C
      DO IQ = NQHOST + 1,NQHOST + NQCLU
C
         IF ( ABS(QMTET(IQ)).LT.1D-6 .AND. ABS(QMPHI(IQ)).LT.1D-6 .AND. 
     &        ABS(QMGAM(IQ)).LT.1D-6 ) THEN
            MAGROT_Q(IQ) = .FALSE.
         ELSE
            MAGROT_Q(IQ) = .TRUE.
         END IF
C
         WRITE (6,99006) IQ,NLQ(IQ),IMQ(IQ),NINT(QMGAM(IQ)),
     &                   NINT(QMTET(IQ)),NINT(QMPHI(IQ)),RWS(IMQ(IQ)),
     &                   ICPA(IQ),NOQ(IQ),IQREPQ(IQ),ISYMGENQ(IQ),
     &                   (IQEQ(IEQ,IQ),IEQ=1,NQEQ(IQ))
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            WRITE (6,'(43X,I5,1X,A,2X,I4,F5.2)') IT,TAB_CHSYM(Z(IT)),
     &             NVALT(IT),CONC(IT)
         END DO
      END DO
C
      WRITE (6,99007) SWS,NVALTOT_CLU
      IF ( KMROT.EQ.1 .OR. KMROT.EQ.2 ) MOMENTS_ROTATED = .TRUE.
      WRITE (6,99008) KMROT,TXTKMROT(KMROT),MOMENTS_ROTATED,QMVEC
C
      WRITE (6,99009)
      DO IQCLU = 1,NQCLU
         IQ = NQHOST + IQCLU
         WRITE (6,99010) IQ,IQCLU,N5VEC_QCLU(1:5,IQCLU),IQ_QCLU(IQCLU),
     &                   DQBAS_QCLU(1:3,IQCLU),QBAS_QCLU(1:3,IQCLU)
      END DO
C
      WRITE (6,99011) NTCLU
      DO ITCLU = 1,NTCLU
         IT = NTHOST + ITCLU
         IQCLU = IQCLU_ATCLU(1,ITCLU)
         WRITE (6,99012) IT,ITCLU,TXT_T(IT)(1:7),NLQ(IQAT(1,IT)),IMT(IT)
     &                   ,NAT(IT),CONC(IT),IQ_QCLU(IQCLU),
     &                   NQHOST + IQCLU,IQCLU_ATCLU(1,ITCLU),
     &                   QBAS_QCLU(1:3,IQCLU)
         DO IA = 2,NAT(IT)
            IQCLU = IQCLU_ATCLU(IA,ITCLU)
            WRITE (6,99013) IQ_QCLU(IQCLU),NQHOST + IQCLU,
     &                      IQCLU_ATCLU(1,ITCLU),QBAS_QCLU(1:3,IQCLU)
         END DO
      END DO
C
      WRITE (6,99014) NMCLU
      IF ( NMCLU.GT.0 ) WRITE (6,99015)
      DO IMCLU = 1,NMCLU
         IM = NMHOST + IMCLU
         WRITE (6,99016) IM,IMCLU,RMT(IMT(IT)),RWS(IMT(IT))
      END DO
C
      WRITE (6,'(//)')
C
CUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                   calculation of U parameter
CUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C
      IF ( U_CALCULATION ) THEN
C
         IQCLU_CENTER = 0
         DO IQCLU = 1,NQCLU
            IF ( ABS(QBAS_QCLU(1,IQCLU)).LT.1D-5 .AND. 
     &           ABS(QBAS_QCLU(2,IQCLU)).LT.1D-5 .AND. 
     &           ABS(QBAS_QCLU(3,IQCLU)).LT.1D-5 ) THEN
               IQCLU_CENTER = IQCLU
               EXIT
            END IF
         END DO
         IF ( IQCLU_CENTER.EQ.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'IQCLU_CENTER = 0')
C
         IT_CENTER = 0
         DO IT = NTHOST + 1,NTHOST + NTCLU
            IF ( IQAT(1,IT).EQ.NQHOST+IQCLU_CENTER ) THEN
               IT_CENTER = IT
               EXIT
            END IF
         END DO
         IF ( IT_CENTER.EQ.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ITCLU_CENTER = 0')
         IF ( NAT(IT_CENTER).NE.1 )
     &         CALL STOP_MESSAGE(ROUTINE,'NAT(IT_CENTER) <> 1')
C
         LOPT(IT_CENTER) = 2
C
         WRITE (6,99019) IQCLU_CENTER,IQCLU_CENTER + NQHOST,IT_CENTER,
     &                   LOPT(IT_CENTER)
C
         IT_CENTER_U_CALCULATION = IT_CENTER
C
         CALL SCF_UCALC(IT_CENTER,IQCLU_CENTER)
C
      END IF
C
C-----------------------------------------------------------------------
C       deallocate all storage connected with host calculation
C-----------------------------------------------------------------------
      IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
         DEALLOCATE (RGNT,IRGNT,NRGNT,IJQ,NIJQ,GGJLRS,PWP)
         DEALLOCATE (SRREL,NRREL,IRREL)
         DEALLOCATE (G1,G2,G3,R1,R2,R3,QQPZ,QQPY,QQPX,HP)
         DEALLOCATE (EXPGNQ,INDR,SMAX,QQMLRS,CILMAT)
         DEALLOCATE (D1TERM3,IILERS,DLLMMKE,CHP,SHP,WK)
         DEALLOCATE (PWEX2K1,PWEX2K2,PWEX2K3,PWEXIK1,PWEXIK2,PWEXIK3)
      END IF
C
C-----------------------------------------------------------------------
      KBZI = 1
      KMOL = 0
      FULLPOT5 = FULLPOT
C
C-----------------------------------------------------------------------
C              reopen RHS and LHS wave function files
C           as LHS_SOL_EQ_RHS_SOL may be .F. in contrast to host
C-----------------------------------------------------------------------
C
      OPEN (IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      CALL SET_IFIL_LHS(IFILCBWF,IFILCBWF_LHS)
      OPEN (IFILCBWF_LHS,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF)
C
      CALL SET_IFIL_SPH(IFILCBWF,IFILCBWF_SPH)
      OPEN (IFILCBWF_SPH,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=RECLNGWF_SPH)
C
C-----------------------------------------------------------------------
C
      IF ( PROGNAME(4:6).EQ.'SCF' ) THEN
C
C------ make sure that host Madelung potential and charges are available
         IFLAG = 0
         IF ( VLMMAD_HOST(1,1).GT.900000D0 ) IFLAG = 1
         IF ( IFLAG.NE.0 ) THEN
            WRITE (6,99018) POTFIL(1:LPOTFIL)
            CALL STOP_MESSAGE(ROUTINE,'host Madelung potential')
         END IF
C
         CALL SCF(KBZI,KMOL,FULLPOT5)
C
      ELSE IF ( PROGNAME(4:6).EQ.'GEN' ) THEN
C
         RUNELOOP = .TRUE.
C        DOSREP = 'XXX'
C
         CLOSE (IFILDOS)
         DOSFIL = DATSET(1:LDATSET)//'.dos'
C
         IF ( MPI_ID.EQ.0 ) THEN
            OPEN (UNIT=IFILDOS,FILE=DOSFIL)
         ELSE
            OPEN (UNIT=IFILDOS,STATUS='SCRATCH')
         END IF
C
         ALLOCATE (IQ_ATCLU(NQMAX,NTMAX))
         DO IT = NTHOST + 1,NTHOST + NTCLU
            DO IA = 1,NAT(IT)
               IQ_ATCLU(IA,IT-NTHOST) = IQAT(IA,IT) - NQHOST
            END DO
         END DO
C
         CALL WRHEAD(IFILDOS,DOSFIL,'DOS       ',NETAB(1))
C
         WRITE (IFILDOS,99004) 'DOS-FMT:  ','OLD-SPRKKR'
C
         IF ( TASK.NE.'CHI       ' ) THEN
C
            CALL GEN(RUNELOOP,TASK,DOSREP,DOSCORSYS)
C
         ELSE
C
            CALL CLU_CHI_DRIVE
C
         END IF
C
      END IF
C
      CALL STOP_REGULAR(ROUTINE,'cluster calculation completed')
C
C=======================================================================
99001 FORMAT (10X,'IT=',I2,' L=0,..,',I1,':',20F7.3)
99002 FORMAT (/,10X,'the SOC is manipulated -- part of the SOC kept: ')
99003 FORMAT (10X,'IT=',I2,' L=0,..,',I1,':',20(2X,A))
99004 FORMAT (A10,A,A)
99005 FORMAT (//,'   number of cluster sites (NQCLU):',I4,//,
     &        '   site NL mesh MGAM MTET MPHI   RWS  CPA NOQ',
     &        ' IT ELMT VAL CONC QREP SYMG eq.sit.')
99006 FORMAT (I6,I4,I4,1X,3I5,F7.4,I3,I4,18X,I4,I5,1X,3I3,:,/,(67X,3I3))
99007 FORMAT (3X,'average',19X,F8.4,10X,'total',F8.2,/)
99008 FORMAT (3X,'KMROT: ',I5,3X,A,3X,'MOMENTS_ROTATED: ',L2,/,3X,
     &        'QMVEC: ',3F7.3,/)
99009 FORMAT (//,5X,'positions of cluster sites',//,
     &        '    IQ IQCLU     N5VEC_QCLU   (IQ)    ',
     &        'DQX    DQY    DQZ       QX     QY     QZ')
99010 FORMAT (2I6,I5,4I3,I4,2(2X,3F7.3))
99011 FORMAT (//,'   number of atomic types in cluster (NTCLU):',I4,
     &        //8X,'types',32X,'on sites',/,
     &        '    IT ITCLU  TXTT   NL mesh NAT',
     &        '  CONC   (IQ)   IQ IQCLU      QX     QY     QZ ')
99012 FORMAT (2I6,2X,A,I2,I5,I4,F6.3,3I6,2X,3F7.3)
99013 FORMAT (38X,3I6,2X,3F7.3)
99014 FORMAT (//,'   number of additional mesh',
     &        ' types in cluster (NMCLU):',I4,/)
99015 FORMAT ('    IM IMCLU    RMT    RWS ')
99016 FORMAT (2I6,2X,3F7.3)
99017 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*  *****  *   *  ****   *****  ****    ****   *****  ****    *'
     &  ,/,10X,
     &  '*  *      ** **  *   *  *      *   *   *   *  *      *   *   *'
     &  ,/,10X,
     &  '*  ****   * * *  ****   ****   *   *   *   *  ****   *   *   *'
     &  ,/,10X,
     &  '*  *      *   *  *   *  *      *   *   *   *  *      *   *   *'
     &  ,/,10X,
     &  '*  *****  *   *  ****   *****  ****    ****   *****  ****    *'
     &  ,/,10X,'*',60X,'*',/,10X,
     &  '*  *****   *       *     *   ****   *******  *******  *****  *'
     &  ,/,10X,
     &  '* *     *  *       *     *  *    *     *     *        *    * *'
     &  ,/,10X,
     &  '* *        *       *     *  *          *     *        *    * *'
     &  ,/,10X,
     &  '* *        *       *     *   ****      *     ****     *****  *'
     &  ,/,10X,
     &  '* *        *       *     *       *     *     *        * *    *'
     &  ,/,10X,
     &  '* *     *  *       *     *  *    *     *     *        *   *  *'
     &  ,/,10X,
     &  '*  *****   ******   *****    ****      *     *******  *    * *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99018 FORMAT (2(/,1X,79('#')),/,10X,'TROUBLE in <CLUDRIVE>',/,10X,
     &        'host Madelung potential not found in potfile',2X,A,/,10X,
     &        'the program will stop !',2(/,1X,79('#')))
99019 FORMAT (//,2(/,1X,79('U')),/,10X,
     &        'calculation of U-parameter for ',/,10X,
     &        'central site  IQCLU =',I3,'    IQ =',I3,'    type IT =',
     &        I3,'    l =',I2,/,2(/,1X,79('U')),//)
      END
C*==clutauij_drive.f    processed by SPAG 6.70Rc at 08:55 on  8 Mar 2017
      SUBROUTINE CLUTAUIJ_DRIVE(DOSREP,DOSCORSYS)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine to supply the host TAUIJ data                    *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SITES,ONLY:NQCLU,IQBOT_HOST,IQTOP_HOST
      USE MOD_ANGMOM,ONLY:NXM,NKM
      USE MOD_ENERGY,ONLY:IGRID,EMIN,EMAX,EFERMI,NETAB
      USE MOD_FILES,ONLY:IFILTAUIJ,TAUIJFIL,DATSET,LDATSET,LTAUIJFIL,
     &    WRTAUIJ,RECLNG_TAUIJ,LRECREAL8
      USE MOD_TAUIJ,ONLY:NTAUIJ,IQ_QTAB1_TAUIJ,JQ_QTAB2_TAUIJ,
     &    NQTAB2_TAUIJ,JQTAB_TAUIJMAX,IQ_TAUIJ,JQ_TAUIJ,N5VEC_TAUIJ,
     &    N123TAUIJMAX,ITAUIJ_LOOP,NLOOP_TAUIJ,NQTAB1_TAUIJ,
     &    NTAUIJ_QTAB2_TAUIJ
      USE MOD_FILES,ONLY:FOUND_STRING
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUTAUIJ_DRIVE')
C
C Dummy arguments
C
      CHARACTER*3 DOSCORSYS,DOSREP
C
C Local variables
C
      INTEGER I,J
      LOGICAL RUNELOOP,TAUIJ_AVAILABLE
      CHARACTER*10 TASK
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      CALL INPUT_FIND_SECTION('CONTROL',1)
C
      CALL SECTION_SET_STRING('TAUIJFIL',TAUIJFIL,'9999',0)
C
      IF ( .NOT.FOUND_STRING ) TAUIJFIL = DATSET(1:LDATSET)//'.tauij'
C
      LTAUIJFIL = LEN_TRIM(TAUIJFIL)
C
      CALL CLUTAUIJ_READHEAD(TAUIJ_AVAILABLE)
C
C=======================================================================
C                create host TAUIJ data and write to disk
C=======================================================================
C
      IF ( .NOT.TAUIJ_AVAILABLE ) THEN
C
         CALL INIT_MOD_TAUIJ_CLUSTER
C
C-----------------------------------------------------------------------
C
         RECLNG_TAUIJ = 1 + 3 + 2*(NKM**2)*(IQTOP_HOST-IQBOT_HOST+1)
         RECLNG_TAUIJ = RECLNG_TAUIJ + (NXM*NQCLU)**2
         RECLNG_TAUIJ = RECLNG_TAUIJ*2*LRECREAL8
C
C *****  RECLNG_TAUIJ = RECLNG_TAUIJ + (NXM**2)*NTAUIJ
C
         OPEN (UNIT=IFILTAUIJ,FILE=TAUIJFIL(1:LTAUIJFIL),
     &         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLNG_TAUIJ,
     &         ERR=100)
C
         WRITE (IFILTAUIJ,REC=1) RECLNG_TAUIJ,0,0D0,0D0,0D0,0
C
         WRTAUIJ = .TRUE.
C
         RUNELOOP = .TRUE.
Cabc         DOSREP = 'XXX'
         TASK = 'NONE'
C
         CALL GEN(RUNELOOP,TASK,DOSREP,DOSCORSYS)
C
         WRITE (IFILTAUIJ,REC=1) RECLNG_TAUIJ,IGRID(1),EMIN,EMAX,EFERMI,
     &                           NETAB(1),NTAUIJ,N123TAUIJMAX,
     &                           NLOOP_TAUIJ,NQTAB1_TAUIJ,
     &                           JQTAB_TAUIJMAX,
     &                           (IQ_TAUIJ(I),JQ_TAUIJ(I),
     &                           (N5VEC_TAUIJ(J,I),J=1,5),I=1,NTAUIJ),
     &                           (ITAUIJ_LOOP(I),I=1,NLOOP_TAUIJ),
     &                           (IQ_QTAB1_TAUIJ(I),NQTAB2_TAUIJ(I),
     &                           (JQ_QTAB2_TAUIJ(I,J),
     &                           NTAUIJ_QTAB2_TAUIJ(I,J),J=1,
     &                           JQTAB_TAUIJMAX),I=1,NQTAB1_TAUIJ)
C
         WRITE (6,99001) TAUIJFIL(1:LTAUIJFIL)
C
      END IF
C
C=======================================================================
C
      RETURN
 100  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,
     &                  'TAUIJ file could not be opened for writing')
C-----------------------------------------------------------------------
99001 FORMAT (2(/,1X,79('*')),/,10X,'INFO from   <CLUTAUIJ_DRIVE>',/,
     &        10X,'the TAUIJ - data for the host have been created',/,
     &        10X,'the TAUIJ-file   ',A,'   is available now',
     &        2(/,1X,79('*')),/)
      END
C*==clutauij_readhead.f    processed by SPAG 6.70Rc at 08:55 on  8 Mar 2017
      SUBROUTINE CLUTAUIJ_READHEAD(TAUIJ_AVAILABLE)
C   ********************************************************************
C   *                                                                  *
C   *   read header of the TAUIJ file                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILTAUIJ,TAUIJFIL,LTAUIJFIL,IPRINT,
     &    RECLNG_TAUIJ
      USE MOD_TAUIJ,ONLY:NTAUIJ,IQ_QTAB1_TAUIJ,JQ_QTAB2_TAUIJ,
     &    NQTAB2_TAUIJ,JQTAB_TAUIJMAX,IQ_TAUIJ,JQ_TAUIJ,N5VEC_TAUIJ,
     &    N123TAUIJMAX,ITAUIJ_LOOP,NLOOP_TAUIJ,NQTAB1_TAUIJ,
     &    NTAUIJ_QTAB2_TAUIJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUTAUIJ_READHEAD')
C
C Dummy arguments
C
      LOGICAL TAUIJ_AVAILABLE
C
C Local variables
C
      REAL*8 EFERMI_IN,EMAX_IN,EMIN_IN
      INTEGER I,IGRID_IN,J,NE_IN,RECLNG_TMP
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
C   ====================================================================
C                  check consistency of TAU-file
C   ====================================================================
C
      RECLNG_TMP = 100
C
      OPEN (UNIT=IFILTAUIJ,FILE=TAUIJFIL(1:LTAUIJFIL),
     &      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLNG_TMP,ERR=100)
C
      READ (IFILTAUIJ,REC=1,ERR=100) RECLNG_TAUIJ,IGRID_IN,EMIN_IN,
     &                               EMAX_IN,EFERMI_IN,NE_IN
C
      CLOSE (IFILTAUIJ)
C
      IF ( NE_IN.EQ.0 ) THEN
         WRITE (6,99003) TAUIJFIL(1:LTAUIJFIL)
         RETURN
      END IF
C
C-----------------------------------------------------------------------
      OPEN (UNIT=IFILTAUIJ,FILE=TAUIJFIL(1:LTAUIJFIL),
     &      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLNG_TAUIJ,
     &      ERR=100)
C
C-----------------------------------------------------------------------
      READ (IFILTAUIJ,REC=1,ERR=100) RECLNG_TAUIJ,IGRID_IN,EMIN_IN,
     &                               EMAX_IN,EFERMI_IN,NE_IN,NTAUIJ,
     &                               N123TAUIJMAX,NLOOP_TAUIJ,
     &                               NQTAB1_TAUIJ,JQTAB_TAUIJMAX
C
C-----------------------------------------------------------------------
      ALLOCATE (IQ_TAUIJ(NTAUIJ),JQ_TAUIJ(NTAUIJ))
      ALLOCATE (N5VEC_TAUIJ(5,NTAUIJ))
C
      ALLOCATE (ITAUIJ_LOOP(NLOOP_TAUIJ))
C
      ALLOCATE (IQ_QTAB1_TAUIJ(NQTAB1_TAUIJ))
      ALLOCATE (NQTAB2_TAUIJ(NQTAB1_TAUIJ))
      ALLOCATE (JQ_QTAB2_TAUIJ(NQTAB1_TAUIJ,JQTAB_TAUIJMAX))
      ALLOCATE (NTAUIJ_QTAB2_TAUIJ(NQTAB1_TAUIJ,JQTAB_TAUIJMAX))
C
C-----------------------------------------------------------------------
      READ (IFILTAUIJ,REC=1,ERR=100) RECLNG_TAUIJ,IGRID_IN,EMIN_IN,
     &                               EMAX_IN,EFERMI_IN,NE_IN,NTAUIJ,
     &                               N123TAUIJMAX,NLOOP_TAUIJ,
     &                               NQTAB1_TAUIJ,JQTAB_TAUIJMAX,
     &                               (IQ_TAUIJ(I),JQ_TAUIJ(I),
     &                               (N5VEC_TAUIJ(J,I),J=1,5),I=1,
     &                               NTAUIJ),
     &                               (ITAUIJ_LOOP(I),I=1,NLOOP_TAUIJ),
     &                               (IQ_QTAB1_TAUIJ(I),NQTAB2_TAUIJ(I),
     &                               (JQ_QTAB2_TAUIJ(I,J),
     &                               NTAUIJ_QTAB2_TAUIJ(I,J),J=1,
     &                               JQTAB_TAUIJMAX),I=1,NQTAB1_TAUIJ)
C
C-----------------------------------------------------------------------
C
      IF ( IPRINT.GT.-10 ) WRITE (6,99001) TAUIJFIL(1:LTAUIJFIL),
     &                            IGRID_IN,EMIN_IN,EMAX_IN,EFERMI_IN,
     &                            NE_IN,NTAUIJ,NLOOP_TAUIJ,N123TAUIJMAX,
     &                            NQTAB1_TAUIJ,JQTAB_TAUIJMAX
C
      TAUIJ_AVAILABLE = .TRUE.
      RETURN
C
C-----------------------------------------------------------------------
 100  CONTINUE
      CLOSE (IFILTAUIJ)
      TAUIJ_AVAILABLE = .FALSE.
      WRITE (6,99002) TAUIJFIL(1:LTAUIJFIL)
C
C-----------------------------------------------------------------------
99001 FORMAT (/,1X,79('*'),/,32X,'<CLUTAUIJ_READHEAD>',/,1X,79('*'),//,
     &        10X,'information from TAUIJ-file: ',A,//,10X,
     &        'IGRID          =',I12,/,10X,'EMIN           =',F12.6,/,
     &        10X,'EMAX           =',F12.6,/,10X,'EFERMI         =',
     &        F12.6,/,10X,'NE             =',I12,/,10X,
     &        'NTAUIJ         =',I12,/,10X,'NLOOP_TAUIJ    =',I12,/,10X,
     &        'N123TAUIJMAX   =',I12,/,10X,'NQTAB1_TAUIJ   =',I12,/,10X,
     &        'JQTAB_TAUIJMAX =',I12,/)
99002 FORMAT (2(/,1X,79('*')),/,10X,'INFO from   <CLUTAUIJ_READHEAD>',/,
     &        10X,'the TAUIJ-file   ',A,'   is empty',/,10X,
     &        'the TAUIJ - data for the host will be created',
     &        2(/,1X,79('*')),/)
99003 FORMAT (2(/,1X,79('*')),/,10X,'INFO from   <CLUTAUIJ_READHEAD>',/,
     &        10X,'the TAUIJ-file   ',A,'   is incomplete',/,10X,
     &        'NE=0: i.e. the number of energies has not yet been set',
     &        /,10X,'the TAUIJ - data for the host will be created',
     &        2(/,1X,79('*')),/)
      END
C*==cluinitfullpot.f    processed by SPAG 6.70Rc at 08:55 on  8 Mar 2017
      SUBROUTINE CLUINITFULLPOT(INIT_RELAXATION)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SCF,ONLY:SCFSTATUS_CLU
      USE MOD_FILES,ONLY:IOTMP,IPRINT
      USE MOD_TYPES,ONLY:VAMEF,VAMEG,TXT_T,LTXT_T,NLMFPT,NFPT,LMIFP,
     &    KLMFP,NLMFPMAX,NCPLWFMAX,NAT,ITBOT,ITTOP
      USE MOD_RMESH,ONLY:ISFLM,KLMSF,LMISF,NSF,NMMAX,NLMSFMAX
      USE MOD_ANGMOM,ONLY:NL,NKM
      USE MOD_SITES,ONLY:IQAT,KFP_LMQ,IMQ
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUINITFULLPOT')
C
C Dummy arguments
C
      LOGICAL INIT_RELAXATION
C
C Local variables
C
      INTEGER I,IBLK,IBLKTOP,IPOT,IT,J,J1,J1TOP,NCPLWF
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      WRITE (6,99001)
C
C=======================================================================
      IF ( INIT_RELAXATION ) THEN
C
C ----------------------------------------------------------------------
C      create shape functions, update array sizes
C      keeping NMMAX that might be too large to avoid reallocation
C      and  set up radial mesh
C ----------------------------------------------------------------------
C
         WRITE (6,99002)
C
         CALL SFNCREATECLU
C
      END IF
C=======================================================================
C
      IF ( SCFSTATUS_CLU(1:5).NE.'START' ) THEN
C
         CALL FPPICKRULES('CHECK',NMMAX,TXT_T,LTXT_T,IMQ,IQAT,NAT,
     &                    KFP_LMQ,KLMFP,NLMFPT,NFPT,LMIFP,ISFLM,KLMSF,
     &                    NSF,LMISF,NLMFPMAX,NLMSFMAX)
C
      ELSE
C
         CALL FPPICKRULES('SETUP',NMMAX,TXT_T,LTXT_T,IMQ,IQAT,NAT,
     &                    KFP_LMQ,KLMFP,NLMFPT,NFPT,LMIFP,ISFLM,KLMSF,
     &                    NSF,LMISF,NLMFPMAX,NLMSFMAX)
C
      END IF
C
      CALL FPCOUPL(IPRINT,NKM,NL)
C
      REWIND (IOTMP)
      READ (IOTMP) NCPLWF,IBLKTOP,J1TOP
C
      IF ( NCPLWF.GT.NCPLWFMAX )
     &      CALL STOP_MESSAGE(ROUTINE,'NCPLWF > NCPLWFMAX')
C
      READ (IOTMP) ((((((VAMEG(I,J,IPOT,J1,IBLK,IT),I=1,NCPLWF),J=1,
     &             NCPLWF),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,IBLKTOP),
     &             IT=ITBOT,ITTOP)
      IF ( IREL.LE.2 ) THEN
         VAMEF = 999999D0
      ELSE
         READ (IOTMP) ((((((VAMEF(I,J,IPOT,J1,IBLK,IT),I=1,NCPLWF),J=1,
     &                NCPLWF),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,
     &                IBLKTOP),IT=ITBOT,ITTOP)
      END IF
      CLOSE (IOTMP)
C
99001 FORMAT (/,1X,79('*'),/,32X,'<CLUINITFULLPOT>',/,1X,79('*'),//,10X,
     &        'initialize full potential cluster calculations',/)
99002 FORMAT (10X,'set up cluster specific meshes for relaxation',/)
      END
