C*==cpabroyden.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CPABROYDEN(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,
     &                      NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUQ,KMROT,
     &                      DROTQ,NTMAX,NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme applied to the CPA - cycle            *
C   *                                                                  *
C   * -----------------------------------------------------------------*
C   *                                                                  *
C   * imix=(0,..,6) : key for using different mixing schemes           *
C   *                          0 means straight mixing                 *
C   *                          3 broyden's first method used           *
C   *                          4 broyden's second method used          *
C   *                          5 strmix with given jacobian            *
C   *                            of broyden's first method             *
C   *                          6 strmix with given jacobian            *
C   *                            of broyden's second method            *
C   *                                                                  *
C   * temporary data is stored on scratch file IFILBROY_CPA            *
C   * modify record length   lngr8   according to machine              *
C   *                                                                  *
C   * _________________________________________________________________*
C   *                                                                  *
C   * on entry: itcpa=1 MSSQ       OLD     initialize tables           *
C   *           else    MSSQ       NEW                                 *
C   *                                                                  *
C   * on exit:          MSSQ     = MSSQ[mix]                           *
C   * _________________________________________________________________*
C   *                                                                  *
C   * recommended            imix         4                            *
C   *                        istbry       1                            *
C   *                        itdept      40                            *
C   *                        mixing    5 .. 20 %                       *
C   ********************************************************************
      USE MOD_FILES,ONLY:IFILBROY_CPA
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--CPABROYDEN39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CPABROYDEN')
      INTEGER ITDEPTMAX,ITDCMMAX,IMIX,ISTBRY
      PARAMETER (ITDEPTMAX=40,ITDCMMAX=3,IMIX=4,ISTBRY=1)
      REAL*8 TOL,MIXING
      PARAMETER (TOL=1D-6,MIXING=0.05D0)
      INTEGER ITDEPT
      PARAMETER (ITDEPT=40)
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPACORR,CPAERR
      INTEGER IPRINT,ITCPA,KMROT,NKMMAX,NQ,NQMAX,NTMAX
      REAL*8 CONC(NTMAX)
      COMPLEX*16 DROTQ(NKMMAX,NKMMAX,NQMAX),MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER ICPA(NQMAX),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)
C
C Local variables
C
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),CHK,ERR,FM(:),FM1(:),
     &       IMM,REM,SM(:),SM1(:),SN,UI(:,:),VI(:,:),WN,WO,XINP(:),
     &       XOUT(:)
      COMPLEX*16 CJ(:,:),CSUM,MNEW,WK(:,:,:)
      INTEGER I,I1,I2,IA_ERR,ICALL,IO,IPF,IQ,IT,IVV,J,M,MIT,MM,N,NMAP,
     &        NMAPMAX
      SAVE AM,BM,FM,FM1,IPF,MIT,NMAP,NMAPMAX,SM,SM1,WN,WO,XINP,XOUT
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE CJ,WK,UI,VI,FM,SM,FM1,SM1,XINP,XOUT
C
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
C
         NMAP = 0
         DO IQ = 1,NQ
            IF ( ICPA(IQ).NE.0 ) NMAP = NMAP + NKMQ(IQ)**2
         END DO
         NMAP = 2*NMAP
C
         ALLOCATE (FM(NMAP),FM1(NMAP))
         ALLOCATE (SM(NMAP),SM1(NMAP))
         ALLOCATE (XINP(NMAP),XOUT(NMAP),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: XINP')
         NMAPMAX = NMAP
      END IF
C
      ALLOCATE (UI(NMAP,2:ITDCMMAX),VI(NMAP,2:ITDCMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: VI')
C
      ALLOCATE (CJ(NKMMAX,NKMMAX),WK(NKMMAX,NKMMAX,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: CJ')
C
C   ====================================================================
C
C------------ array set up and definition of input parameter -----------
C
      IF ( ITCPA.EQ.1 ) THEN
         MIT = 1
         IPF = 6
C
         IF ( ITDCMMAX.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'ITDCMMAX.NE.3')
C
         IF ( IPRINT.GT.2 ) THEN
C
            WRITE (6,FMT=99001)
C
            IF ( IMIX.EQ.3 .OR. IMIX.EQ.4 ) WRITE (6,FMT=99003) (IMIX-2)
     &           ,ITDEPT,ISTBRY
C
            IF ( IMIX.GE.5 ) WRITE (6,FMT=99004) (IMIX-4),ITDEPT - 1,
     &                              ISTBRY
C
            WRITE (6,FMT=99002) MIXING
C
         END IF
C
C-------------------------------------------------------- store old MSSQ
C
         IVV = 0
         DO IQ = 1,NQ
            IF ( ICPA(IQ).NE.0 ) THEN
               DO J = 1,NKMQ(IQ)
                  DO I = 1,NKMQ(IQ)
                     REM = DREAL(MSSQ(I,J,IQ))
                     IMM = DIMAG(MSSQ(I,J,IQ))
                     IVV = IVV + 1
                     XINP(IVV) = REM
                     IVV = IVV + 1
                     XINP(IVV) = IMM
                  END DO
               END DO
            END IF
         END DO
C
         IF ( IVV.NE.NMAP ) CALL STOP_MESSAGE(ROUTINE,'IVV.NE.NMAP')
C
         OPEN (IFILBROY_CPA,STATUS='SCRATCH',FORM='UNFORMATTED')
C
      END IF
C
C   ====================================================================
C
      CPACORR = 0D0
      CPACHNG = 0D0
      CPAERR = 0D0
      WO = 1.0D0 - MIXING
      WN = MIXING
C     WO = 1.0D0
C     WN = 0d0
C
      M = NKMMAX
      MM = NKMMAX*NKMMAX
C
      IVV = 0
      DO IQ = 1,NQ
         IF ( ICPA(IQ).NE.0 ) THEN
C
            N = NKMQ(IQ)
C
            CALL CMATCOP(N,M,TAUQ(1,1,IQ),WK)
            CALL CMATINV(N,M,WK,CJ)
            DO I2 = 1,N
               DO I1 = 1,N
                  CJ(I1,I2) = MSSQ(I1,I2,IQ) - CJ(I1,I2)
               END DO
            END DO
C
            CALL CINIT(MM,WK)
C
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               DO I2 = 1,N
                  DO I1 = 1,N
                     WK(I1,I2,2) = -CJ(I1,I2)
                  END DO
               END DO
C
C ------------------------- rotate the single site m-matrix if necessary
               IF ( KMROT.NE.0 ) THEN
C
                  CALL ROTATE(MSST(1,1,IT),'L->G',WK(1,1,3),N,
     &                        DROTQ(1,1,IQ),M)
C
                  DO I2 = 1,N
                     DO I1 = 1,N
                        WK(I1,I2,2) = WK(I1,I2,2) + WK(I1,I2,3)
                     END DO
                  END DO
C
               ELSE
C
                  DO I2 = 1,N
                     DO I1 = 1,N
                        WK(I1,I2,2) = WK(I1,I2,2) + MSST(I1,I2,IT)
                     END DO
                  END DO
C
               END IF
C
               CALL CMATINV(N,M,WK(1,1,2),WK(1,1,3))
C
               DO I2 = 1,N
                  DO I1 = 1,N
                     WK(I1,I2,1) = WK(I1,I2,1) + CONC(IT)*WK(I1,I2,3)
                  END DO
               END DO
            END DO
C
            CALL CMATINV(N,M,WK(1,1,1),WK(1,1,3))
C
            CHK = 0D0
            DO I2 = 1,N
               DO I1 = 1,N
                  WK(I1,I2,1) = CJ(I1,I2) + WK(I1,I2,3) - MSSQ(I1,I2,IQ)
                  SN = ABS(MSSQ(I1,I2,IQ))
                  ERR = ABS(WK(I1,I2,1))*SN/(1D0+SN**2)
                  IF ( ERR.GT.CHK ) CHK = ERR
               END DO
            END DO
C
            IF ( CHK.GT.TOL ) THEN
               DO J = 1,N
                  DO I = 1,N
C ------------------------------------------------------------- new MSSQ
                     MNEW = MSSQ(I,J,IQ) + WK(I,J,1)
                     REM = DREAL(MNEW)
                     IMM = DIMAG(MNEW)
C ------------------------------------------------------------- mix MSSQ
                     IVV = IVV + 1
                     XOUT(IVV) = REM
                     XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
                     IVV = IVV + 1
                     XOUT(IVV) = IMM
                     XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
C
C                    MSSQ(i,j,IQ) = MSSQ(i,j,IQ) + CPADMP*WK(i,j,1)
                  END DO
               END DO
            END IF
C
            CPAERR = MAX(CHK,CPAERR)
C
            IF ( IPRINT.GE.2 ) THEN
               CSUM = C0
               DO I = 1,N
                  CSUM = CSUM + MSSQ(I,I,IQ)
               END DO
               WRITE (6,99005) IQ,CPAERR,CPACORR,CSUM
            END IF
C
         END IF
      END DO
C
C   ====================================================================
C
C----> broyden updating schemes
C
      IF ( IMIX.GE.3 .AND. ITCPA.GE.ISTBRY )
     &     CALL SCFBROYPT2(IPRINT,XINP,XOUT,FM,FM1,SM,SM1,UI,VI,AM,BM,
     &     MIXING,ITDEPT,ITDEPTMAX,IMIX,IFILBROY_CPA,IPF,MIT,NMAP,
     &     NMAPMAX)
C
      IVV = 0
      DO IQ = 1,NQ
         IF ( ICPA(IQ).NE.0 ) THEN
            DO J = 1,NKMQ(IQ)
               DO I = 1,NKMQ(IQ)
                  IVV = IVV + 1
                  XINP(IVV) = XOUT(IVV)
                  REM = XOUT(IVV)
                  IVV = IVV + 1
                  XINP(IVV) = XOUT(IVV)
                  IMM = XOUT(IVV)
                  MSSQ(I,J,IQ) = DCMPLX(REM,IMM)
               END DO
            END DO
         END IF
      END DO
C
C   ====================================================================
C
      DEALLOCATE (CJ,WK,UI,VI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (/,1X,79('*'),/,28X,'BROYDEN - mixing scheme')
99002 FORMAT (/,10X,'parameter for broyden-update:   ',f12.7,/,1X,
     &        79('*'),/)
99003 FORMAT (/,10X,'broyden''s method # ',i1,' used',/,10X,
     &        'iterationdepth to accumulate the jacobian:  ',i3,/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99004 FORMAT (/,10X,'broyden''s method # ',i3,
     &        ' is used up to iteration depth: ',i3,/,10X,
     &        'then jacobian is fixed and potential',
     &        ' is updated using that jacobian',/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99005 FORMAT (' CPA:  IQ',I3,'  ERR',F12.5,'  CORR',F13.5,'  M',
     &        18(1X,2(1PE14.6)))
      END
