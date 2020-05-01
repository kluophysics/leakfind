C*==gil.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GIL
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the Gilbert damping parameter alpha     *
C   *                                                                  *
C   *  NOTE: TAUQ is always recalculated                               *
C   *        never taken from a previous run                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,CHIZ
      USE MOD_THERMAL,ONLY:NFTMAX,UMAT_VT,I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT
      USE MOD_CPA,ONLY:CPATOL,ITCPAMAX
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,
     &    TSST,TAUT
      USE MOD_SITES,ONLY:NQMAX,NOMAX,NOQ,ITOQ,IQBOT,IQTOP,NQ
      USE MOD_TYPES,ONLY:NTMAX,CONC,ITBOT,ITTOP,NT
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY_LINRESP
      USE MOD_CALCMODE,ONLY:ORBPOL,KKRMODE,KMROT
      USE MOD_ENERGY,ONLY:EFERMI,PHASK
      USE MOD_FILES,ONLY:IFILCBWF,IPRINT,WRMAT,IOTMP
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,MPI_ID
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--GIL27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GIL')
C
C Local variables
C
      COMPLEX*16 ALF0Q(:,:,:),ALF0QO(:,:,:,:),ALF1QOQ_NV(:,:,:,:,:),
     &           ALF1QOQ_VC(:,:,:,:,:),ALF1QQ_NV(:,:,:,:),
     &           ALF1QQ_VC(:,:,:,:),ERYD,MAQAB(:,:,:,:),
     &           MAQOAB(:,:,:,:,:),MBQAB(:,:,:,:),MBQOAB(:,:,:,:,:),
     &           MCQAB(:,:,:,:),MCQOAB(:,:,:,:,:),MDQAB(:,:,:,:),
     &           MDQOAB(:,:,:,:,:),MREG_FT(:,:,:,:),P,TAUQZ(:,:,:,:)
      LOGICAL CALCINT,GETIRRSOL,INITIALIZE,VERTEX
      REAL*8 CPACHNG,CPACHNGMAX,TIME,TIME0
      INTEGER IA_ERR,ICPACONV,ICPAFLAG,IERR,IFIL_DUMP_TAU,ITCPA,
     &        IWRIRRWF,IWRREGWF,M,NCPAFAIL
C
C*** End of declarations rewritten by SPAG
C
      DATA CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/,IFIL_DUMP_TAU/6/
      DATA VERTEX/.TRUE./
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE MAQAB,MBQAB,MCQAB,MDQAB,MREG_FT
      ALLOCATABLE MAQOAB,MBQOAB,MCQOAB,MDQOAB
      ALLOCATABLE ALF0Q,ALF0QO,ALF1QQ_NV,ALF1QQ_VC,ALF1QOQ_NV,ALF1QOQ_VC
      ALLOCATABLE TAUQZ
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
      ALLOCATE (MAQAB(M,M,3,NQMAX))
      ALLOCATE (MBQAB(M,M,3,NQMAX))
      ALLOCATE (MCQAB(M,M,3,NQMAX))
      ALLOCATE (MDQAB(M,M,3,NQMAX))
      ALLOCATE (MAQOAB(M,M,3,NQMAX,NOMAX))
      ALLOCATE (MBQOAB(M,M,3,NQMAX,NOMAX))
      ALLOCATE (MCQOAB(M,M,3,NQMAX,NOMAX))
      ALLOCATE (MDQOAB(M,M,3,NQMAX,NOMAX))
      ALLOCATE (ALF0Q(3,3,NQMAX),ALF0QO(3,3,NQMAX,NOMAX))
      ALLOCATE (ALF1QQ_NV(3,3,NQMAX,NQMAX))
      ALLOCATE (ALF1QQ_VC(3,3,NQMAX,NQMAX))
      ALLOCATE (ALF1QOQ_NV(3,3,NQMAX,NOMAX,NQMAX))
      ALLOCATE (ALF1QOQ_VC(3,3,NQMAX,NOMAX,NQMAX))
      ALLOCATE (TAUQZ(NKMMAX,NKMMAX,NQMAX,2),STAT=IA_ERR)
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAQAB')
C
C ======================================================================
      WRITE (6,99002)
C ======================================================================
      CALL CPU_TIME(TIME0)
C ======================================================================
C
      IF ( KMROT.NE.0 )
     &      CALL STOP_MESSAGE(ROUTINE,'not working for rotated moments')
C
      IF ( IBZINT.NE.2 .AND. IBZINT.NE.6 )
     &      CALL STOP_MESSAGE(ROUTINE,'IBZINT <> 2 AND IBZINT <> 6')
C
C
      NCPAFAIL = 0
      ITCPAMAX = 30
C
      MPI_KLOOP = MPI
C
C ======================================================================
C   set the range of sites IQ and types IT to maximum
C   to allow for TB-calculations on 2D systems
C ======================================================================
C
      ITBOT = 1
      ITTOP = NT
      IQBOT = 1
      IQTOP = NQ
C
C***********************************************************************
C***********************************************************************
C****************   START OF ALPHA CALCULATION   ***********************
C***********************************************************************
C***********************************************************************
C
      ICPAFLAG = 0
      CPACHNG = 0.0D0
      CPACHNGMAX = 0D0
C
      ERYD = EFERMI
C
      WRITE (6,99001) EFERMI
C
      ERYDA_EQ_ERYDB = .TRUE.
C
C=======================================================================
C===============         LOOP OVER TEMPERATURE          ========== START
C=======================================================================
C
      TEMP_LAT = 0D0
      N_TEMP_LAT = 1
C
      CALL THERMAL_INIT(0,N_TEMP_LAT,TEMP_LAT)
C
      DO I_TEMP_LAT = 1,N_TEMP_LAT
C
C ======================================================================
C    set up  NVIBRA  vectors SVEC_VT for displacements
C            NFLUCT  vectors DVEC_FT for fluctuations
C            and occupation tables X_VFT, NVFO_Q and IVFT_VFOQ
C
         CALL THERMAL_INIT(I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
         IF ( INITIALIZE ) THEN
C
            ALLOCATE (MREG_FT(NKMMAX,NKMMAX,3,NFTMAX),STAT=IA_ERR)
            INITIALIZE = .FALSE.
C
         END IF
C=======================================================================
C
C ===================================== solve SS - differential equation
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,
     &                 ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C=========================================================== TEMPERATURE
C         in case of thermal lattice vibrations and/or spin fluctuations
C                     init E-dependent displacement matrices     UMAT_VT
C
         CALL THERMAL_INIT_UFMAT(ERYD)
C
C ======================================================================
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C***********************************************************************
         IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
            CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                           ICPACONV,CONC,NOQ,ITOQ,PHASK,1,NTMAX,
     &                           TSST,MSST,TSSQ,MSSQ,TAUQ)
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
            CALL TAU_TB(ERYD,P,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,TSST,
     &                  MSST,TSSQ,MSSQ,TAUQ)
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
            TAUQ(:,:,:) = TSSQ(:,:,:)
C
C***********************************************************************
         ELSE
C
            CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
         END IF
C***********************************************************************
C
         IF ( ICPAFLAG.NE.0 ) THEN
            NCPAFAIL = NCPAFAIL + 1
            CPACHNGMAX = MAX(CPACHNGMAX,ABS(CPACHNG))
         END IF
C ??????????????????????????????????????????????????????????????????????
         IF ( IPRINT.GE.5 .OR. WRMAT ) THEN
C
            CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
            CALL DUMPTAU(1,ERYD,IFIL_DUMP_TAU,MSST,MSSQ,TAUT,TAUQ)
C
         END IF
C ??????????????????????????????????????????????????????????????????????
C
C***********************************************************************
         IF ( KKRMODE(1:12).EQ.'STANDARD-KKR' ) THEN
C***********************************************************************
C
            IF ( IBZINT.EQ.4 ) THEN
C
               CALL STOP_MESSAGE(ROUTINE,
     &                   'only <SIGKLOOPS> available for BZ integration'
     &                   )
C
            ELSE IF ( NO_SYMMETRY_LINRESP ) THEN
C
               CALL SIGKLOOPS_NOSYM(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,
     &                              MSSQ)
C
            ELSE
C
               CALL SIGKLOOPS(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQ,MSSQ,MSSQ)
C
            END IF
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:6).EQ.'TB-KKR' ) THEN
C***********************************************************************
C
            CALL SIGKLOOPSDRV_TB(ERYD,P,MSSQ,TAUQZ)
C
C***********************************************************************
         ELSE IF ( KKRMODE(1:11).EQ.'SINGLE-SITE' ) THEN
C***********************************************************************
C
            CHIZ(:,:,:) = C0
C
C***********************************************************************
         ELSE
C
            CALL STOP_MESSAGE(ROUTINE,'KKRMODE incorrect')
C
         END IF
C***********************************************************************
C
C ??????????????????????????????????????????????????????????????????????
         IF ( WRMAT ) THEN
C
            OPEN (IOTMP,FILE='CHIZ_gil.dump')
            WRITE (IOTMP,'(2D25.12)') CHIZ
C
            CALL STOP_REGULAR(ROUTINE,'CHIZ dumped to CHIZ_gil.dump')
C
         END IF
C ??????????????????????????????????????????????????????????????????????
C
C======================================================================
C                                                             MPI_ID = 0
         IF ( MPI_ID.EQ.0 ) THEN
C
            CALL GILME(MSST,TSST,MSSQ,TAUQ,MAQAB,MBQAB,MCQAB,MDQAB,
     &                 MAQOAB,MBQOAB,MCQOAB,MDQOAB,MREG_FT)
C
            CALL GIL0(MREG_FT,TSST,MSSQ,TAUQ,ALF0Q,ALF0QO)
C
            CALL GIL1(MAQAB,MBQAB,MCQOAB,MDQOAB,ALF1QQ_NV,ALF1QOQ_NV,
     &                'N',MAQOAB,MBQAB,MCQOAB,MDQAB)
C
            IF ( VERTEX ) THEN
C
               CALL LINRESP_VERTEX(TAUQ,TAUQ,TSST,TSST,UMAT_VT,UMAT_VT,
     &                             MSSQ,MSSQ)
C
               CALL GIL1(MAQAB,MBQAB,MCQOAB,MDQOAB,ALF1QQ_VC,ALF1QOQ_VC,
     &                   'V',MAQOAB,MBQAB,MCQOAB,MDQAB)
C
            ELSE
C
               ALF1QQ_VC(:,:,:,:) = ALF1QQ_NV(:,:,:,:)
               ALF1QOQ_VC(:,:,:,:,:) = ALF1QOQ_NV(:,:,:,:,:)
C
            END IF
C
C ======================================================================
C
            CALL GILSUM(ALF0Q,ALF1QQ_NV,ALF1QQ_VC,ALF0QO,ALF1QOQ_NV,
     &                  ALF1QOQ_VC,I_TEMP_LAT,N_TEMP_LAT,TEMP_LAT)
C
            IF ( NCPAFAIL.NE.0 ) WRITE (6,99003) NCPAFAIL,CPATOL,
     &                                  CPACHNGMAX
C
            CALL CPU_TIME(TIME)
            WRITE (6,99004) TIME - TIME0
C
         END IF
C                                                             MPI_ID = 0
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI ) CALL DRV_MPI_BARRIER
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C
C===============  LOOP OVER TEMPERATURE - IF REQUESTED  ========== END =
C=======================================================================
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI ) THEN
         CALL DRV_MPI_BARRIER
         CALL MPI_FINALIZE(IERR)
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      WRITE (6,*) '          ALPHA-job done      FEIERABEND '
      STOP
C
C   ====================================================================
99001 FORMAT (/,79('*'),/,79('*'),/,30X,'energy loop',/,79('*'),/,
     &        79('*'),//,25X,'EFERMI =',F10.5,/)
99002 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      ****   ***  *      *****   *****  *****   *******     *'
     &  ,/,10X,
     &  '*     *    *   *   *      *    *  *      *    *     *        *'
     &  ,/,10X,
     &  '*     *        *   *      *    *  *      *    *     *        *'
     &  ,/,10X,
     &  '*     *  ***   *   *      *****   ****   *****      *        *'
     &  ,/,10X,
     &  '*     *    *   *   *      *    *  *      *  *       *        *'
     &  ,/,10X,
     &  '*     *    *   *   *      *    *  *      *   *      *        *'
     &  ,/,10X,
     &  '*      ****   ***  *****  *****   *****  *    *     *        *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99003 FORMAT (/,1X,79('*'),/,10X,'CPA not converged for',I3,
     &        ' energies:',/,10X,'tolerance for CPA-cycle:',F15.7,/,10X,
     &        'maximum deviation:      ',F15.7,/,1X,79('*'),/)
99004 FORMAT (/,5X,'execution time for <GIL>:',F14.3,' secs',/)
C
      END
