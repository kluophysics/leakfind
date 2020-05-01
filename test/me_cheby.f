C*==me_cheby.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_CHEBY(TSST,MSST,SSST,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_FILES,ONLY:LRECREAL8,IPRINT
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,CTL,IMT,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NLM,NPOL,NCPLWF
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS,JRCRI
      USE MOD_ENERGY,ONLY:ECHEBY1,ECHEBY2,TCHEBY1,TCHEBY2,NCHEBY1,
     &    NCHEBY2,ZGCHEBY1,ZFCHEBY1,JGCHEBY1,JFCHEBY1,ZGCHEBY2,ZFCHEBY2,
     &    JGCHEBY2,JFCHEBY2,ECHEBY1_BOT,ECHEBY1_TOP,ECHEBY2_BOT,
     &    ECHEBY2_TOP
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--ME_CHEBY19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_CHEBY')
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 ADEL,C,CHKSUM(:,:),MAXDEV,RELDIFF,RELDIFFMAX,SUMDEV,TJE1(:)
     &       ,TJE2(:),W12
      LOGICAL CALCINT,GETIRRSOL,ME_CC_BRA_RWF
      CHARACTER*1 CHPOL(NPOL)
      COMPLEX*16 ERYD1,ERYD2,ETAB1(:),ETAB2(:),JF1(:,:,:),JF2(:,:,:),
     &           JG1(:,:,:),JG2(:,:,:),MCHB,MIRR_2(:,:,:,:),
     &           MIRR_3(:,:,:,:),MIRR_4(:,:,:,:),MREF,MZAZB(:,:,:),
     &           MZBZA(:,:,:),MZBZA_CHEBY(:,:,:,:,:),
     &           MZBZA_REF(:,:,:,:,:),MZBZA_X(:,:,:),P,ZF1(:,:,:),
     &           ZF2(:,:,:),ZG1(:,:,:),ZG2(:,:,:)
      INTEGER I,IA_ERR,IE1,IE2,IEC1,IEC2,IFIL1,IFIL10,IFIL2,IFIL20,
     &        IFILCHEB1,IFILCHEB10,IFILCHEB2,IFILCHEB20,IFLAG,IM,IOL,
     &        IPOL,IR,IRTOP,IT,IWRIRRWF,IWRREGWF,J,J1,J2,KIRR,LAM,M,NSUM
C
C*** End of declarations rewritten by SPAG
C
      DATA CHPOL/'-','z','+'/
      DATA CALCINT/.TRUE./
C
      ALLOCATABLE MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4,MZBZA_X
      ALLOCATABLE MZBZA_CHEBY,MZBZA_REF
      ALLOCATABLE JF1,JG1,ZG1,ZF1,JF2,JG2,ZG2,ZF2
      ALLOCATABLE TJE1,TJE2,ETAB2,ETAB1,CHKSUM
C
      WRITE (6,99001)
C
      M = NKMMAX
C
      ALLOCATE (MIRR_2(M,M,3,3),MIRR_3(M,M,3,3),MIRR_4(M,M,3,3))
      ALLOCATE (MZAZB(M,M,3),MZBZA(M,M,3),MZBZA_X(M,M,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MIRR')
C
      ALLOCATE (JF1(NRMAX,NCPLWFMAX,NKM),JG1(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZG1(NRMAX,NCPLWFMAX,NKM),ZF1(NRMAX,NCPLWFMAX,NKM))
C
      ALLOCATE (JF2(NRMAX,NCPLWFMAX,NKM),JG2(NRMAX,NCPLWFMAX,NKM))
      ALLOCATE (ZG2(NRMAX,NCPLWFMAX,NKM),ZF2(NRMAX,NCPLWFMAX,NKM))
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZF2')
C
      IFIL10 = 100
      IFILCHEB10 = 1000
      NCHEBY1 = 10
C      NCHEBY1 = 1
      ALLOCATE (ECHEBY1(NCHEBY1),TJE1(NCHEBY1))
      ALLOCATE (TCHEBY1(NCHEBY1,NCHEBY1))
      ALLOCATE (ZGCHEBY1(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY1))
      ALLOCATE (ZFCHEBY1(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY1))
      ALLOCATE (JGCHEBY1(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY1))
      ALLOCATE (JFCHEBY1(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY1))
C
      IFIL20 = 200
      IFILCHEB20 = 2000
      NCHEBY2 = 10
C      NCHEBY2 = 1
      ALLOCATE (ECHEBY2(NCHEBY2),TJE2(NCHEBY2))
      ALLOCATE (TCHEBY2(NCHEBY2,NCHEBY2))
      ALLOCATE (ZGCHEBY2(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY2))
      ALLOCATE (ZFCHEBY2(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY2))
      ALLOCATE (JGCHEBY2(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY2))
      ALLOCATE (JFCHEBY2(NRMAX,NCPLWFMAX,NKMMAX,NCHEBY2))
C
      M = NKMMAX
      ALLOCATE (MZBZA_CHEBY(M,M,3,NCHEBY2,NCHEBY1))
      ALLOCATE (MZBZA_REF(M,M,3,NCHEBY2,NCHEBY1))
C
      ALLOCATE (ETAB2(NCHEBY2),ETAB1(NCHEBY1))
      ALLOCATE (CHKSUM(NCHEBY2,NCHEBY1))
C
      IT = 1
      C = CTL(IT,1)
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
C=======================================================================
      CALCINT = .FALSE.
      ME_CC_BRA_RWF = .FALSE.
      IWRREGWF = 1
C
C     IWRIRRWF = 0
      IWRIRRWF = 1
C
      GETIRRSOL = IWRREGWF.EQ.1
      KIRR = IWRIRRWF
C=======================================================================
C
C      ECHEBY1_BOT = (0.1D0,0.01D0)
C      ECHEBY1_TOP = (0.7D0,0.01D0)
C
C      ECHEBY2_BOT = (1.1D0,0.02D0)
C      ECHEBY2_TOP = (1.7D0,0.02D0)
C
      ECHEBY1_BOT = (0.1D0,0.01D0)
      ECHEBY1_TOP = (0.7D0,0.01D0)
C
      ECHEBY2_BOT = (1.1D0,0.02D0)
      ECHEBY2_TOP = (1.7D0,0.02D0)
C
C=======================================================================
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      DO IE1 = 1,NCHEBY1
         IFIL1 = IFIL10 + IE1
         OPEN (UNIT=IFIL1,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
      END DO
C
      DO IE2 = 1,NCHEBY2
         IFIL2 = IFIL20 + IE2
         OPEN (UNIT=IFIL2,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
      END DO
C
      DO IEC1 = 1,NCHEBY1
         IFILCHEB1 = IFILCHEB10 + IEC1
         OPEN (UNIT=IFILCHEB1,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
      END DO
C
      DO IEC2 = 1,NCHEBY2
         IFILCHEB2 = IFILCHEB20 + IEC2
         OPEN (UNIT=IFILCHEB2,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
      END DO
C
C-----------------------------------------------------------------------
C
      WRITE (6,99003) 'NPOL      ',NPOL,(CHPOL(I),I=1,NPOL)
      WRITE (6,99002) 'NLM       ',NLM
      WRITE (6,99002) 'NKM       ',NKM
      WRITE (6,99005) IT
C
C ======================================================================
C           initial BAND state on Chebychev E-grid ECHEBY1
C         supply reference values for Chebychev approximation
C         and calculate corresponding Chebychev wave functions
C ======================================================================
C
      CALL CHEBY_INIT_ENERGY(NCHEBY1,ECHEBY1_BOT,ECHEBY1_TOP,ECHEBY1,
     &                       TCHEBY1)
C
      DO IE1 = 1,NCHEBY1
C
         IFIL1 = IFIL10 + IE1
C
         ERYD1 = ECHEBY1(IE1)
C
         WRITE (6,99004) 'initial',ERYD1,1,NKM
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL1,GETIRRSOL,ERYD1,
     &                 P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      END DO
C
      CALL ME_CHEBY_WAVFUN(IT,KIRR,IFIL10,JF1,JG1,ZG1,ZF1,IFILCHEB10,
     &                     NCHEBY1,TCHEBY1,ZGCHEBY1,ZFCHEBY1,JGCHEBY1,
     &                     JFCHEBY1)
C
C ======================================================================
C             final BAND state on Chebychev E-grid ECHEBY2
C         supply reference values for Chebychev approximation
C         and calculate corresponding Chebychev wave functions
C ======================================================================
C
      CALL CHEBY_INIT_ENERGY(NCHEBY2,ECHEBY2_BOT,ECHEBY2_TOP,ECHEBY2,
     &                       TCHEBY2)
C
      DO IE2 = 1,NCHEBY2
C
         IFIL2 = IFIL20 + IE2
C
         ERYD2 = ECHEBY2(IE2)
C
         WRITE (6,99004) 'final  ',ERYD2,1,NKM
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL2,GETIRRSOL,ERYD2,
     &                 P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      END DO
C
      CALL ME_CHEBY_WAVFUN(IT,KIRR,IFIL20,JF2,JG2,ZG2,ZF2,IFILCHEB20,
     &                     NCHEBY2,TCHEBY2,ZGCHEBY2,ZFCHEBY2,JGCHEBY2,
     &                     JFCHEBY2)
C
C ======================================================================
C   check fit of wave function for Chebychev energies  ECHEBY1
C ======================================================================
C
      WRITE (6,'(//,''  self-check for fit of wave function  '',/)')
C
      IFLAG = 0
C
      DO IE1 = 1,NCHEBY1
C
         ERYD1 = ECHEBY1(IE1)
C
         CALL CHEBY_POLYNOM_TJE(ERYD1,ECHEBY1_BOT,ECHEBY1_TOP,NCHEBY1,
     &                          TJE1)
C
         IFIL1 = IFIL10 + IE1
C
         CALL WAVFUN_READ_REL(IFIL1,IT,KIRR,ZG1,ZF1,JG1,JF1,IRTOP,
     &                        NCPLWF,IKMCPLWF)
C
         ZG2(:,:,:) = 0D0
         DO J1 = 1,NCHEBY1
            IFILCHEB1 = IFILCHEB10 + J1
            CALL WAVFUN_READ_REL(IFILCHEB1,IT,KIRR,ZGCHEBY1(1,1,1,J1),
     &                           ZFCHEBY1(1,1,1,J1),JGCHEBY1(1,1,1,J1),
     &                           JFCHEBY1(1,1,1,J1),IRTOP,NCPLWF,
     &                           IKMCPLWF)
            ZG2(:,:,:) = ZG2(:,:,:) + ZGCHEBY1(:,:,:,J1)*TJE1(J1)
         END DO
C
         RELDIFFMAX = 0D0
         DO LAM = 1,NKM
            DO I = 1,NCPLWF(LAM)
               DO IR = 1,IRTOP,20
                  RELDIFF = ABS(DREAL(1D0-ZG2(IR,I,LAM)/ZG1(IR,I,LAM)))
                  RELDIFFMAX = MAX(RELDIFF,RELDIFFMAX)
                  IF ( ABS(RELDIFF).GT.1D-12 ) THEN
                     WRITE (6,*) '######################'
                     WRITE (6,'(2i5,2e16.5,e22.12)') LAM,IR,
     &                      DREAL(ZG1(IR,I,LAM)),DREAL(ZG2(IR,I,LAM)),
     &                      RELDIFF
                     IFLAG = 1
                  END IF
               END DO
            END DO
         END DO
C
         WRITE (6,'(/,'' IE1 ='',i3,''  max. rel. diff.:'',E20.6)') IE1,
     &          RELDIFFMAX
C
      END DO
C
      IF ( IFLAG.NE.0 ) STOP
C
C ======================================================================
C                      calculate matrix elements
C       for proper wave functions on E-grids  ECHEBY1  and  ECHEBY2
C ======================================================================
C
      DO IE1 = 1,NCHEBY1
C
         IFIL1 = IFIL10 + IE1
C
         ERYD1 = ECHEBY1(IE1)
C
         DO IE2 = 1,NCHEBY2
C
            IFIL2 = IFIL20 + IE2
C
            ERYD2 = ECHEBY2(IE2)
C
            CALL ME_ALF_ALF(1,NKM,IFIL2,ERYD2,1,NKM,IFIL1,ERYD1,
     &                      ME_CC_BRA_RWF,IT,MZAZB,
     &                      MZBZA_REF(1,1,1,IE2,IE1),MIRR_2,MIRR_3,
     &                      MIRR_4,C,K_CALC_ME)
C
            WRITE (6,'(A,2I6,6F8.3)') '##### ME ENERGIE',IE1,IE2,ERYD1,
     &                                ERYD2
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,1,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,2,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,3,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C
C
         END DO
C
      END DO
C
C ======================================================================
C                      calculate matrix elements
C     for Chebychev wave functions on E-grids  ECHEBY1  and  ECHEBY2
C ======================================================================
C
      DO IE1 = 1,NCHEBY1
C
         IFIL1 = IFILCHEB10 + IE1
C
         ERYD1 = ECHEBY1(IE1)
C
         DO IE2 = 1,NCHEBY2
C
            IFIL2 = IFILCHEB20 + IE2
C
            ERYD2 = ECHEBY2(IE2)
C
            CALL ME_ALF_ALF(1,NKM,IFIL2,ERYD2,1,NKM,IFIL1,ERYD1,
     &                      ME_CC_BRA_RWF,IT,MZAZB,
     &                      MZBZA_CHEBY(1,1,1,IE2,IE1),MIRR_2,MIRR_3,
     &                      MIRR_4,C,K_CALC_ME)
C
            WRITE (6,'(A,2I6,6F8.3)') '##### ME ENERGIE',IE1,IE2,ERYD1,
     &                                ERYD2
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,1,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,2,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C            CALL  CMATSTR('MZBZA C',MZBZA_CHEBY(1,1,3,IE2,IE1),NKM,
C     &                   NKMMAX,3,3,0,1D-8,6)
C
C
         END DO
C
      END DO
C
C ======================================================================
C    self-check of MEs for Chebychev energies  ECHEBY1  and  ECHEBY2
C ======================================================================
C
      DO IE1 = 1,NCHEBY1
C
         ERYD1 = ECHEBY1(IE1)
C
         CALL CHEBY_POLYNOM_TJE(ERYD1,ECHEBY1_BOT,ECHEBY1_TOP,NCHEBY1,
     &                          TJE1)
C
         DO IE2 = 1,NCHEBY2
C
            ERYD2 = ECHEBY2(IE2)
C
            CALL CHEBY_POLYNOM_TJE(ERYD2,ECHEBY2_BOT,ECHEBY2_TOP,
     &                             NCHEBY2,TJE2)
C
C-----------------------------------------------------------------------
C                    calculate Chebychev matrix element
C-----------------------------------------------------------------------
C
            MZBZA_X(:,:,:) = C0
C
            CHKSUM(:,:) = 0D0
C
            DO J1 = 1,NCHEBY1
               DO J2 = 1,NCHEBY2
C
                  W12 = TJE1(J1)*TJE2(J2)
C
                  WRITE (*,*) IE1,IE2,J1,J2,TJE1(J1),TJE2(J2),W12
C
                  MZBZA_X(:,:,:) = MZBZA_X(:,:,:)
     &                             + W12*MZBZA_CHEBY(:,:,:,J2,J1)
C
                  CHKSUM(1,J1) = CHKSUM(1,J1) + W12
C
               END DO
            END DO
C
            DO J1 = 1,NCHEBY1
               WRITE (6,'('' CHECKSUM '',I4,20F12.6)') J1,CHKSUM(1,J1)
            END DO
C-----------------------------------------------------------------------
C
            SUMDEV = 0D0
            MAXDEV = 0D0
            NSUM = 0
C
            DO IPOL = 1,NPOL
               IF ( IPOL.NE.1 ) CYCLE
               DO J = 1,NKM
                  DO I = 1,NKM
C
                     MREF = MZBZA_REF(I,J,IPOL,IE2,IE1)
                     MCHB = MZBZA_X(I,J,IPOL)
C
                     IF ( ABS(MREF).GT.1D-6 ) THEN
                        WRITE (6,99007) CHPOL(IPOL),I,J,'REF',MREF
C
                        NSUM = NSUM + 1
                        ADEL = ABS(1D0-MCHB/MREF)
                        SUMDEV = SUMDEV + ADEL
                        MAXDEV = MAX(MAXDEV,ADEL)
                        WRITE (6,99007) CHPOL(IPOL),I,J,'CHB',MCHB,ADEL
                        WRITE (6,*) ' '
C
                     END IF
C
                  END DO
               END DO
            END DO
C
            SUMDEV = SUMDEV/DBLE(NSUM)
            WRITE (6,99006) 'ADA',NSUM,NSUM,SUMDEV,MAXDEV
C
         END DO
C
         EXIT
      END DO
C
C ======================================================================
C                   initial BAND state - regular E-grid
C ======================================================================
C
      DO IE1 = 1,NCHEBY1
C
         IFIL1 = IFIL10 + IE1
C
         ETAB1(IE1) = ECHEBY1_BOT + (ECHEBY1_TOP-ECHEBY1_BOT)*(IE1-1)
     &                /DBLE(NCHEBY1-1)
C
         ERYD1 = ETAB1(IE1)
C
         WRITE (6,99004) 'initial',ERYD1,1,NKM
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL1,GETIRRSOL,ERYD1,
     &                 P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      END DO
C
C ======================================================================
C                   final BAND state - regular E-grid
C ======================================================================
C
      DO IE2 = 1,NCHEBY2
C
         IFIL2 = IFIL20 + IE2
C
         ETAB2(IE2) = ECHEBY2_BOT + (ECHEBY2_TOP-ECHEBY2_BOT)*(IE2-1)
     &                /DBLE(NCHEBY2-1)
C
         ERYD2 = ETAB2(IE2)
C
         WRITE (6,99004) 'initial',ERYD2,1,NKM
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL2,GETIRRSOL,ERYD2,
     &                 P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      END DO
C
C ======================================================================
C                 check of MEs for regular regular E-grid
C ======================================================================
C
      DO IE1 = 1,NCHEBY1
C
         IFIL1 = IFIL10 + IE1
C
         ERYD1 = ETAB1(IE1)
C
         CALL CHEBY_POLYNOM_TJE(ERYD1,ECHEBY1_BOT,ECHEBY1_TOP,NCHEBY1,
     &                          TJE1)
C
         DO IE2 = 1,NCHEBY2
C
            IFIL2 = IFIL20 + IE2
C
            ERYD2 = ETAB2(IE2)
C
            CALL CHEBY_POLYNOM_TJE(ERYD2,ECHEBY2_BOT,ECHEBY2_TOP,
     &                             NCHEBY2,TJE2)
C
C-----------------------------------------------------------------------
C                    calculate Chebychev matrix element
C-----------------------------------------------------------------------
C
            MZBZA_X(:,:,:) = C0
C
            DO J1 = 1,NCHEBY1
               DO J2 = 1,NCHEBY2
C
                  W12 = TJE1(J1)*TJE2(J2)
C
                  MZBZA_X(:,:,:) = MZBZA_X(:,:,:)
     &                             + W12*MZBZA_CHEBY(:,:,:,J2,J1)
C
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                    calculate STANDARD matrix element
C-----------------------------------------------------------------------
C
            CALL ME_ALF_ALF(1,NKM,IFIL2,ERYD2,1,NKM,IFIL1,ERYD1,
     &                      ME_CC_BRA_RWF,IT,MZAZB,MZBZA,MIRR_2,MIRR_3,
     &                      MIRR_4,C,K_CALC_ME)
C
C-----------------------------------------------------------------------
C
            SUMDEV = 0D0
            MAXDEV = 0D0
            NSUM = 0
C
            DO IPOL = 1,NPOL
               DO J = 1,NKM
                  DO I = 1,NKM
C
                     MREF = MZBZA(I,J,IPOL)
                     MCHB = MZBZA_X(I,J,IPOL)
C
                     IF ( ABS(MREF).GT.1D-6 ) THEN
                        WRITE (6,99007) CHPOL(IPOL),I,J,'REF',MREF
C
                        NSUM = NSUM + 1
                        ADEL = ABS(1D0-MCHB/MREF)
                        SUMDEV = SUMDEV + ADEL
                        MAXDEV = MAX(MAXDEV,ADEL)
                        WRITE (6,99007) CHPOL(IPOL),I,J,'CHB',MCHB,ADEL
                        WRITE (6,*) ' '
C
                     END IF
C
                  END DO
               END DO
            END DO
C
            SUMDEV = SUMDEV/DBLE(NSUM)
            WRITE (6,99006) 'ADA',NSUM,NSUM,SUMDEV,MAXDEV
C
         END DO
C
      END DO
C
      STOP
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '* *     *  ******   ****   *     *  ******   ****   *     *  *'
     &  ,/,10X,
     &  '* **   **  *       *    *  *     *  *       *    *  *    *   *'
     &  ,/,10X,
     &  '* * * * *  *       *       *     *  *       *       *  *     *'
     &  ,/,10X,
     &  '* *  *  *  *****   *       *******  *****   *       * *      *'
     &  ,/,10X,
     &  '* *     *  *       *       *     *  *       *       *  *     *'
     &  ,/,10X,
     &  '* *     *  *       *    *  *     *  *       *    *  *    *   *'
     &  ,/,10X,
     &  '* *     *  ******   ****   *     *  ******   ****   *     *  *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (5X,A10,I10)
99003 FORMAT (5X,A10,I10,5(:,'  (',A,')'))
99004 FORMAT (/,5X,A,' state    energy  E = ',2F12.6,//,33X,'IKM =',I9,
     &        '   ...   ',I3,/)
99005 FORMAT (//,5X,'selected atom type   IT =',I3,/)
99006 FORMAT (/,5X,A,'  number of non-0 terms:       ',I12,/,:,13X,
     &        '  number of deviations:        ',I12,/,13X,
     &        '  average relative deviation:  ',F12.8,/,:,13X,
     &        '  maximum relative deviation:  ',F12.8,/)
99007 FORMAT (5X,'(',A,')',2I3,2X,A,2F20.12,F16.12)
      END
C*==cheby_init_energy.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE CHEBY_INIT_ENERGY(N,A,B,X,T)
C   ********************************************************************
C   *                                                                  *
C   *  initialize Chebychev grid points X_k for the energy             *
C   *                                                                  *
C   *           with A <= X <= B and -1 <= Y <= +1                     *
C   *                                                                  *
C   *  an energy grid with a constant imaginary part is assumed        *
C   *                                                                  *
C   *                                                                  *
C   *  tabulate Chebychev polynoms T(j,k)                              *
C   *  of order j for the grid points k                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CI,PI
      IMPLICIT NONE
C*--CHEBY_INIT_ENERGY622
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHEBY_INIT_ENERGY')
C
C Dummy arguments
C
      COMPLEX*16 A,B
      INTEGER N
      REAL*8 T(N,N)
      COMPLEX*16 X(N)
C
C Local variables
C
      INTEGER I,J,K
      REAL*8 IMX,S,TOL,XAVG,XDEL,Y
C
C*** End of declarations rewritten by SPAG
C
      IMX = DIMAG(A)
      IF ( ABS(IMX-DIMAG(B)).GT.1D-8 )
     &      CALL STOP_MESSAGE(ROUTINE,'inconsistent imaginary part of E'
     &     )
C
      XDEL = DREAL(B-A)/2.0D0
      XAVG = DREAL(B+A)/2.0D0
      WRITE (*,'(A,6F8.3)') '############### ENERGY A    ',A
      WRITE (*,'(A,6F8.3)') '############### ENERGY B    ',B
      WRITE (*,'(A,6F8.3)') '############### ENERGY XDEL ',XDEL
      WRITE (*,'(A,6F8.3)') '############### ENERGY XAVG ',XAVG
C
C-------------------------------- set up energy grid with  A <= X_k <= B
      DO K = 1,N
         Y = COS(PI*(K-0.5D0)/DBLE(N))
         X(K) = XAVG + Y*XDEL + IMX*CI
         WRITE (*,'(A,I6,6F8.3)') '############### ENERGY X(K)',K,Y,X(K)
      END DO
C
C---------------------------------------- tabulate Chebychev polynomials
      DO J = 1,N
         DO K = 1,N
            T(J,K) = COS(PI*(J-1)*(K-0.5D0)/DBLE(N))
         END DO
      END DO
C
C=======================================================================
      WRITE (6,'(//,''  check the sum rules  '',/)')
C=======================================================================
C
      DO J = 1,N
         WRITE (6,'(''    T(J,K): '',20f8.3)') (T(J,K),K=1,N)
      END DO
C
      TOL = 1D-12
      WRITE (6,'(/,''  tolerance '',1pe20.1,/)') TOL
C
      DO I = 1,N
         DO K = 1,N
            S = T(1,K)*T(1,I)/2D0
            DO J = 2,N
               S = S + T(J,K)*T(J,I)
            END DO
            S = S*2/DBLE(N)
            IF ( (I.EQ.K .AND. ABS(S-1D0).GT.TOL) .OR. 
     &           (I.NE.K .AND. ABS(S).GT.TOL) ) WRITE (6,'(2I3,F18.12)')
     &           I,K,S
         END DO
      END DO
C
      END
C*==cheby_polynom_tje.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE CHEBY_POLYNOM_TJE(X,A,B,N,T)
C   ********************************************************************
C   *                                                                  *
C   *  evaluate the Chebychev polynoms  T(J,X)  J=1,..,N               *
C   *                                                                  *
C   *  for X  with A <= X <= B and -1 <= Y <= +1                       *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHEBY_POLYNOM_TJE721
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHEBY_POLYNOM_TJE')
C
C Dummy arguments
C
      COMPLEX*16 A,B,X
      INTEGER N
      REAL*8 T(N)
C
C Local variables
C
      REAL*8 IMX,XAVG,XDEL,Y,Y2
      INTEGER K
C
C*** End of declarations rewritten by SPAG
C
      IMX = DIMAG(A)
C
      IF ( ABS(IMX-DIMAG(B)).GT.1D-8 )
     &      CALL STOP_MESSAGE(ROUTINE,'Im E inconsistent for A and B')
C
      XDEL = DREAL(B-A)/2.0D0
      XAVG = DREAL(B+A)/2.0D0
C
      Y = DREAL(X-XAVG)/XDEL
C
      T(1) = 1.0D0
C
      IF ( N.NE.1 ) THEN
C
         T(2) = Y
C
         IF ( N.NE.2 ) THEN
C
            Y2 = Y + Y
C
            DO K = 3,N
C
               T(K) = Y2*T(K-1) - T(K-2)
C
            END DO
C
         END IF
C
      END IF
C
      END
C*==me_cheby_wavfun.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_CHEBY_WAVFUN(IT,KIRR,IFIL0,JF,JG,ZG,ZF,IFILCHEB0,N,
     &                           TC,ZGC,ZFC,JGC,JFC)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the Chebychev wave functions                          *
C   *                                                                  *
C   *           CZ(j) = 2/N SUM_k=1..N Z(k) * T(j,k)                   *
C   *                                                                  *
C   *  j=1: divide weight by 2 (exploiting T(1,k) = 1)                 *
C   *                                                                  *
C   *  k    numbers the zeros of the Chebychev polynom of order N      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX,IKMCPLWF,IMT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,NCPLWF
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,JRWS,JRCRI
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--ME_CHEBY_WAVFUN808
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFIL0,IFILCHEB0,IT,KIRR,N
      COMPLEX*16 JF(NRMAX,NCPLWFMAX,NKM),JFC(NRMAX,NCPLWFMAX,NKMMAX,N),
     &           JG(NRMAX,NCPLWFMAX,NKM),JGC(NRMAX,NCPLWFMAX,NKMMAX,N),
     &           ZF(NRMAX,NCPLWFMAX,NKM),ZFC(NRMAX,NCPLWFMAX,NKMMAX,N),
     &           ZG(NRMAX,NCPLWFMAX,NKM),ZGC(NRMAX,NCPLWFMAX,NKMMAX,N)
      REAL*8 TC(N,N)
C
C Local variables
C
      INTEGER I,IFIL,IM,IR,IRTOP,J,K,LAM
      REAL*8 W
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      ZGC(:,:,:,:) = C0
      ZFC(:,:,:,:) = C0
C
      IF ( KIRR.EQ.1 ) THEN
         JGC(:,:,:,:) = C0
         JFC(:,:,:,:) = C0
      END IF
C
C=======================================================================
C              determine the Chebychev wave functions
C=======================================================================
C
      DO K = 1,N
C
         IFIL = IFIL0 + K
C
         WRITE (*,*) 'CCCCCCCC',N,IFIL
         WRITE (6,*) 'FIT ',K,IFIL0,IFIL,IT,KIRR,IRTOP
         CALL WAVFUN_READ_REL(IFIL,IT,KIRR,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
C-----------------------------------------------------------------------
         DO J = 1,N
C
            W = TC(J,K)*2/DBLE(N)
C
C----------------------- deal with Chebychev polynomial of order 0 (J=1)
C
            IF ( J.EQ.1 ) W = W/2D0
            WRITE (6,*) 'XXXXXXXXXXXXXXXXXXXXXXXX',K,J,W
C
            DO LAM = 1,NKM
               DO I = 1,NCPLWF(LAM)
C
                  IF ( KIRR.EQ.1 ) THEN
                     DO IR = 1,IRTOP
                        ZGC(IR,I,LAM,J) = ZGC(IR,I,LAM,J)
     &                     + W*ZG(IR,I,LAM)
                        ZFC(IR,I,LAM,J) = ZFC(IR,I,LAM,J)
     &                     + W*ZF(IR,I,LAM)
                        JGC(IR,I,LAM,J) = JGC(IR,I,LAM,J)
     &                     + W*JG(IR,I,LAM)
                        JFC(IR,I,LAM,J) = JFC(IR,I,LAM,J)
     &                     + W*JF(IR,I,LAM)
                     END DO
                  ELSE
                     DO IR = 1,IRTOP
                        ZGC(IR,I,LAM,J) = ZGC(IR,I,LAM,J)
     &                     + W*ZG(IR,I,LAM)
                        ZFC(IR,I,LAM,J) = ZFC(IR,I,LAM,J)
     &                     + W*ZF(IR,I,LAM)
                     END DO
                  END IF
C
               END DO
            END DO
C
         END DO
C-----------------------------------------------------------------------
C
      END DO
C=======================================================================
C
C
C=======================================================================
C                 write the Chebychev wave functions
C=======================================================================
      DO J = 1,N
C
         IFIL = IFILCHEB0 + J
C
         CALL WAVFUN_WRITE_REL(IFIL,IT,KIRR,ZGC(1,1,1,J),ZFC(1,1,1,J),
     &                         JGC(1,1,1,J),JFC(1,1,1,J),IRTOP,NCPLWF,
     &                         IKMCPLWF)
C
      END DO
C=======================================================================
C
      END
