C*==eloop.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE ELOOP(IPRINT,PHASK,TAUQ,TAUT,CALCINT,GETIRRSOL)
C   ********************************************************************
C   *                                                                  *
C   * this routine runs the energy loop in a similar way as in main.   *
C   * it is primarily meant to create the TAU-files for the various    *
C   * spectroscopy applications. for that reasons some restrictions    *
C   * have been imposed:                                               *
C   *                                                                  *
C   * - no run of CORE routines                                        *
C   * - no search for FERMI energy                                     *
C   * - NEPATH = 1 is assumed ===  SPLITSS = FALSE                     *
C   * - no output of LOG in CALCDOS to LOGFILE                         *
C   * - no HYPERFINE calculations                                      *
C   *                                                                  *
C   * 04/01/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:SYSTEM,LSYSTEM,WRTAU,WRTAUMQ,IFILTAU
      USE MOD_ENERGY,ONLY:NETAB,NEMAX,WETAB,ETAB,SPLITSS
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_TYPES,ONLY:NT,NTMAX,TXT_T,DOBS_LTX,DOBS_TX,OBS_LTX,OBS_TX,
     &    DOBS_TX_GLO,OBS_TX_GLO
      USE MOD_SITES,ONLY:NQ,NQMAX,IQAT
      USE MOD_FILES,ONLY:IFILCBWF
      IMPLICIT NONE
C*--ELOOP30
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER IPRINT
      COMPLEX*16 PHASK(NEMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BCOR(:),BCORS(:),CPACHNG,CPACHTAB(:)
      COMPLEX*16 CTOTDOS(NEMAX),ERYD,MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),MSSQ(:,:,:),
     &           MSST(:,:,:),P,SSST(NKMMAX,NKMMAX,NTMAX),TSSQ(:,:,:),
     &           TSST(:,:,:)
      INTEGER I,IA_ERR,ICPACONV,ICPAFLAG,IE,IECPAFAIL(:),IECURR,IEPATH,
     &        IPRINTBAND,IQ,IT,ITCPA,IWRIRRWF,IWRREGWF,NCPAFAIL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BCOR,CPACHTAB,TSSQ,MSSQ,MSST,TSST,IECPAFAIL,BCORS
C
      ALLOCATE (BCOR(NTMAX),CPACHTAB(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:eloop -> BCOR'
C
      ALLOCATE (TSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSSQ(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:eloop -> MSSQ'
C
      ALLOCATE (MSST(NKMMAX,NKMMAX,NTMAX),BCORS(NTMAX))
      ALLOCATE (TSST(NKMMAX,NKMMAX,NTMAX),IECPAFAIL(NEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:eloop -> IECPAFAIL'
C
      IPRINTBAND = IPRINT
      IEPATH = 1
C ======================================================================
C
      DO IT = 1,NTMAX
         BCOR(IT) = 0D0
         BCORS(IT) = 0D0
      END DO
C
      IWRREGWF = 0
      IWRIRRWF = 0
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
C
C----------- NO calculation of matrix elements, irregular wave function,
C--------------------------- and NO output of wave functions in <SSITE>,
C------------------ nor of TAU-matrix in <PROJTAU> for CLUSTER ELOOP-run
      IF ( IBZINT.EQ.0 ) THEN
         IWRREGWF = 0
         IWRIRRWF = 0
C
         CALCINT = .TRUE.
         GETIRRSOL = .TRUE.
C
         WRTAU = .FALSE.
         WRTAUMQ = .FALSE.
      END IF
C
C=======================================================================
C
      IF ( WRTAU .OR. (IBZINT.EQ.0) ) THEN
         REWIND IFILTAU
         WRITE (IFILTAU,99001) NT,NQ
         WRITE (IFILTAU,99002) (IT,IQAT(1,IT),TXT_T(IT),IT=1,NT)
         DO IQ = 1,NQ
            WRITE (IFILTAU,'(''site'',I3,'' of '',A)') IQ,
     &             SYSTEM(1:LSYSTEM)
         END DO
      ELSE IF ( WRTAUMQ ) THEN
         REWIND IFILTAU
         WRITE (IFILTAU,99003) NQ
         DO IQ = 1,NQ
            WRITE (IFILTAU,'(''site'',I3,'' of '',A)') IQ,
     &             SYSTEM(1:LSYSTEM)
         END DO
      END IF
C
      NCPAFAIL = 0
C
      DO IE = 1,NETAB(1)
         IECURR = IE
         ICPAFLAG = 0
         CPACHNG = 0.0D0
C
         ERYD = ETAB(IE,1)
C
C====================================== solve SS - differential equation
C
         IWRREGWF = 1
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,
     &                 ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C                                                  START:  CALCULATE TAU
C
         CALL TAU_DRIVE(IECURR,IPRINTBAND,ERYD,P,TSSQ,MSSQ,TSST,MSST,
     &                  TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
         IF ( ICPAFLAG.NE.0 ) THEN
            NCPAFAIL = NCPAFAIL + 1
            CPACHTAB(NCPAFAIL) = CPACHNG
            IECPAFAIL(NCPAFAIL) = IECURR
         END IF
C                                                    END:  CALCULATE TAU
C ======================================================================
C
C
         CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
C ======================================================================
C
         CALL CALCDOS(1,.FALSE.,SPLITSS,IEPATH,NCPAFAIL,IPRINT,ERYD,
     &                MEZZ,MEZJ,TSST,MSST,TAUT,MSSQ,TAUQ,IECURR,
     &                WETAB(IECURR,IEPATH),BCOR,BCORS,DOBS_LTX,DOBS_TX,
     &                OBS_LTX,OBS_TX,DOBS_TX_GLO,OBS_TX_GLO,CTOTDOS)
C
      END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      CLOSE (IFILCBWF)
      IF ( NCPAFAIL.NE.0 ) THEN
         WRITE (6,99004) CPATOL,NCPAFAIL,
     &                   (IECPAFAIL(I),DREAL(ETAB(IECPAFAIL(I),1)),
     &                   CPACHTAB(I),I=1,NCPAFAIL)
         WRITE (6,99006)
      ELSE IF ( NCPA.NE.0 ) THEN
         WRITE (6,99005)
      END IF
C
      DEALLOCATE (BCOR,CPACHTAB,MSSQ,MSST,TSST,IECPAFAIL,BCORS)
C
99001 FORMAT (/,70('*'),/,10X,' TAU(LAM,LAM'')',/,70('*'),/,5X,' NT =',
     &        I3,' NQ =',I3,' FMT = 2')
99002 FORMAT (5X,' IT =',I3,' IQ =',I3,': ',A)
99003 FORMAT (/,70('*'),/,10X,' TAU(LAM,LAM'') and M(LAM,LAM'')',/,
     &        70('*'),/,13X,' NQ =',I3,' FMT = 2')
99004 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99005 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99006 FORMAT (1X,79('*'),/)
      END
