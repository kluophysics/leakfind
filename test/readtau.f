C*==readtau.f    processed by SPAG 6.70Rc at 08:46 on  8 Mar 2017
      SUBROUTINE READTAU(IFIL,ERYD,IE,NETAU,TAUT,IT,NT,RDTAU,TAUQ,MSSQ,
     &                   IQ,NQ,KEY,N,M,IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *   read TAU from file IFIL accoring to the specific situation     *
C   *                                                                  *
C   *   RDTAU = .TRUE.    TAUT(IT) is read                             *
C   *           .FALSE.   TAUQ(IQ)  and  MSSQ(IQ) is read              *
C   *                     (implies that RDTAUMQ = .TRUE.)              *
C   *                                                                  *
C   *   KEY=0    just read and check the header and find the           *
C   *              number of tabulated energy mesh-points              *
C   *   KEY=1    read TAU for current energy                           *
C   *                                                                  *
C   *   KEY=2    read TAU only for ERYDINP in MPI ELOOP                *
C   *                                                                  *
C   *   the TAU-matrices have been written by <PROJTAU> or <CLUSTER>   *
C   *   the end of an entry is indicated by an empty line !            *
C   *   old input formats are accepted indicated by  IFMT              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:CMD00,CMDUC
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='READTAU')
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IE,IFIL,IPRINT,IQ,IT,KEY,M,N,NETAU,NQ,NT
      LOGICAL RDTAU
      COMPLEX*16 MSSQ(M,M),TAUQ(M,M),TAUT(M,M)
C
C Local variables
C
      COMPLEX*16 ERYDINP,ERYDLAST,MSSINP,TAUINP
      CHARACTER*40 FMTE,FMTT
      INTEGER I,IFMT,IOS,IQP,ITP,ITP1,J,NLINES,NQP,NTP
      CHARACTER*30 STR30
      CHARACTER*5 STR5
      SAVE FMTT
C
C*** End of declarations rewritten by SPAG
C
      DATA IFMT/0/,FMTE(1:1)/' '/
C
      IOS = 0
      ERYDINP = ERYD
      IF ( IPRINT.GE.1 ) WRITE (6,*) 'ERYDINP = ',ERYDINP
C
C
C   ********************************************************************
C   ********************************************************************
C                        RDTAU - read  TAUT(IT)
C   ********************************************************************
C   ********************************************************************
C
      IF ( RDTAU ) THEN
C
         IF ( KEY.EQ.0 ) THEN
C   ====================================================================
C                  check consistency of TAU-file
C   ====================================================================
C
            REWIND IFIL
            NETAU = 0
            NLINES = 0
            DO
               READ (IFIL,*,END=20)
               NLINES = NLINES + 1
            END DO
 20         CONTINUE
            REWIND IFIL
C
            READ (IFIL,99001,END=800,ERR=900) STR30,CMD00
            CMDUC = CMD00
            CALL SECTION_SET_INTEGER('NT',NTP,9999,1)
            CALL SECTION_SET_INTEGER('NQ',NQP,9999,1)
            CALL SECTION_SET_INTEGER('FMT',IFMT,9999,1)
C
            IF ( NTP.NE.NT ) CALL STOP_MESSAGE(ROUTINE,
     &           'inconsistent header in  TAU-file    for RDTAU')
            IF ( STR30.NE.' TAU(LAM,LAM'')                ' )
     &           CALL STOP_MESSAGE(ROUTINE,
     &                   'inconsistent header in  TAU-file    for RDTAU'
     &                   )
            IF ( NQP.NE.NQ )
     &            CALL STOP_MESSAGE(ROUTINE,'NQ   <>   NQ (TAU-file)')
C
            READ (IFIL,'(A5)',END=40,ERR=40) (STR5,I=1,(NTP+NQP))
C
            FMTE = '(5X,2F8.4,27X,I2,5X,I2)'
            FMTT = '(2I5,2E15.7)'
            IF ( IFMT.EQ.1 ) THEN
               FMTE = '(2F15.9,18X,I2,5X,I2,6X,I2,F15.6)'
               FMTT = '(2I5,1P4E17.9)'
            ELSE IF ( IFMT.EQ.2 ) THEN
               FMTE = '(2F21.15,18X,I2,5X,I2,6X,I2,F15.6)'
               FMTT = '(2I5,1P4E22.14)'
            ELSE IF ( IFMT.EQ.3 ) THEN
               FMTE = '(2F21.15,18X,I5,5X,I5,6X,I2,F15.8)'
               FMTT = '(2I5,1P4E22.14)'
            END IF
C
C-----------------------------------------------------------------------
C
            ERYDLAST = DCMPLX(-99D0,-99D0)
            ITP1 = 0
            DO I = 1,NLINES
 30            CONTINUE
               READ (IFIL,'(A5)',END=40,ERR=40) STR5
               IF ( STR5.NE.'*****' ) GOTO 30
               READ (IFIL,FMT=FMTE,END=40,ERR=40) ERYD,ITP,IQP
               IF ( NETAU.EQ.0 ) ITP1 = ITP
               IF ( ITP.EQ.ITP1 ) THEN
                  IF ( ABS(ERYD-ERYDLAST).GT.1D-6 ) THEN
                     NETAU = NETAU + 1
                     ERYDLAST = ERYD
                  END IF
               END IF
            END DO
C
 40         CONTINUE
            REWIND IFIL
            READ (IFIL,'(A5)') (STR5,I=1,(4+NTP+NQP))
C
            WRITE (6,99003) '        TAUT(IT)',NETAU,IFMT
C
         ELSE IF ( KEY.EQ.1 ) THEN
C   ====================================================================
C                  read   TAUT(IT)
C   ====================================================================
C
            TAUT(:,:) = C0
C
 60         CONTINUE
            READ (IFIL,'(A5)',END=600) STR5
            IF ( STR5.NE.'*****' ) GOTO 60
C
            READ (IFIL,FMT=FMTE,END=600) ERYD,ITP,IQP
            WRITE (6,99004) IE,ERYD
            IF ( IPRINT.GT.0 ) WRITE (6,99004) IE,ERYD
            IF ( IT.NE.ITP ) GOTO 60
C
 80         CONTINUE
            READ (IFIL,FMT=FMTT,END=500,IOSTAT=IOS) I,J,TAUINP
            IF ( (I+J).NE.0 ) THEN
               TAUT(I,J) = TAUINP
               IF ( (I+J).LT.2*N ) GOTO 80
            END IF
C
         ELSE IF ( KEY.EQ.2 ) THEN
C   ====================================================================
C                  read TAUT(IT) only for ERYD = ERYDINP
C   ====================================================================
C
            TAUT(:,:) = C0
C
            REWIND IFIL
            READ (IFIL,'(A5)') (STR5,I=1,(4+NTP+NQP))
C
 100        CONTINUE
            READ (IFIL,'(A5)',END=600) STR5
            IF ( STR5.NE.'*****' ) GOTO 100
C
            READ (IFIL,FMT=FMTE,END=700) ERYD,ITP,IQP
            IF ( IPRINT.GT.0 ) WRITE (6,99004) IE,ERYD
            IF ( IT.NE.ITP .OR. ABS(ERYDINP-ERYD).GT.1D-10 ) GOTO 100
C
 120        CONTINUE
            READ (IFIL,FMT=FMTT,END=500,IOSTAT=IOS) I,J,TAUINP
            IF ( (I+J).NE.0 ) THEN
               TAUT(I,J) = TAUINP
               IF ( (I+J).LT.2*N ) GOTO 120
            END IF
C
         END IF
C
C   ********************************************************************
C   ********************************************************************
C                  RDTAUMQ - read TAUQ(IQ)  and  MSSQ(IQ)
C   ********************************************************************
C   ********************************************************************
C
      ELSE IF ( KEY.EQ.0 ) THEN
C   ====================================================================
C                  check consistency of TAU-file
C   ====================================================================
C
         REWIND IFIL
         NETAU = 0
         NLINES = 0
         DO
            READ (IFIL,*,END=150)
            NLINES = NLINES + 1
         END DO
 150     CONTINUE
         REWIND IFIL
C
         READ (IFIL,99002,END=800,ERR=1000) STR30,CMD00
         CMDUC = CMD00
         CALL SECTION_SET_INTEGER('NQ',NQP,9999,1)
         CALL SECTION_SET_INTEGER('FMT',IFMT,9999,1)
C
         NTP = 0
         IF ( STR30.NE.' TAU(LAM,LAM'') and M(LAM,LAM'')' ) THEN
            WRITE (6,*) 
     &               'STR30.NE.'' TAU(LAM,LAM'''') and M(LAM,LAM'''')'''
            CALL STOP_MESSAGE(ROUTINE,
     &                 'inconsistent header in  TAU-file    for RDTAUMQ'
     &                 )
         END IF
         IF ( NQP.NE.NQ )
     &         CALL STOP_MESSAGE(ROUTINE,'NQ   <>   NQ (TAU-file)')
C
         READ (IFIL,'(A5)',ERR=200) (STR5,I=1,(NTP+NQP))
C
         FMTE = '(5X,2F8.4,27X,I2,5X,I2)'
         FMTT = '(2I5,2E15.7)'
         IF ( IFMT.EQ.1 ) THEN
            FMTE = '(2F15.9,18X,I2,5X,I2,6X,I2,F15.6)'
            FMTT = '(2I5,1P4E17.9)'
         ELSE IF ( IFMT.EQ.2 ) THEN
            FMTE = '(2F21.15,25X,I2,9X,I2,F15.6)'
            FMTT = '(2I5,1P,4E22.14)'
         ELSE IF ( IFMT.EQ.3 ) THEN
            FMTE = '(2F21.15,25X,I5,9X,I5,F15.6)'
            FMTT = '(2I5,1P4E22.14)'
         END IF
C
C-----------------------------------------------------------------------
C
         DO I = 1,NLINES
 160        CONTINUE
            READ (IFIL,'(A5)',END=200,ERR=200) STR5
            IF ( STR5.NE.'*****' ) GOTO 160
            READ (IFIL,FMT=FMTE,ERR=200) ERYD,IQP
            IF ( IQP.EQ.1 ) NETAU = NETAU + 1
         END DO
C
 200     CONTINUE
         REWIND IFIL
         READ (IFIL,'(A5)') (STR5,I=1,(4+NTP+NQP))
C
         WRITE (6,99003) '   TAUQ(IQ)  and    mss(IQ)',NETAU,IFMT
C
      ELSE IF ( KEY.EQ.1 ) THEN
C   ====================================================================
C                  read   TAUQ(IQ)  and    mss(IQ)
C   ====================================================================
         IF ( FMTE(1:1).EQ.' ' ) THEN
            IFMT = 2
            FMTE = '(2F21.15,25X,I2,9X,I2,F15.6)'
            FMTT = '(2I5,1P,4E22.14)'
         END IF
C
         TAUQ(:,:) = C0
         MSSQ(:,:) = C0
C
 250     CONTINUE
         READ (IFIL,'(A5)',END=600) STR5
         IF ( STR5.NE.'*****' ) GOTO 250
C
         READ (IFIL,FMT=FMTE,END=600) ERYD,IQP
         IF ( IPRINT.GT.0 ) WRITE (6,99004) IE,ERYD
         IF ( IQ.NE.IQP ) GOTO 250
C
 300     CONTINUE
         READ (IFIL,FMT=FMTT,END=500,IOSTAT=IOS) I,J,TAUINP,MSSINP
         IF ( (I+J).NE.0 ) THEN
            TAUQ(I,J) = TAUINP
            MSSQ(I,J) = MSSINP
            IF ( (I+J).LT.2*N ) GOTO 300
         END IF
C
      ELSE IF ( KEY.EQ.2 ) THEN
C   ====================================================================
C          read  TAUQ(IQ)  and  mss(IQ)  only for ERYD = ERYDINP
C   ====================================================================
         IF ( FMTE(1:1).EQ.' ' ) THEN
            IFMT = 2
            FMTE = '(2F21.15,25X,I2,9X,I2,F15.6)'
            FMTT = '(2I5,1P,4E22.14)'
         END IF
C
         TAUQ(:,:) = C0
         MSSQ(:,:) = C0
C
CSW       WARNING: POSSIBLY FATAL! BUT WHY #TYPES + #SITES (SEE BELOW)?
         NTP = 0
CSW       WARNING: POSSIBLY FATAL! BUT WHY #TYPES + #SITES (SEE BELOW)?
         NQP = NQ
C
         REWIND IFIL
         READ (IFIL,'(A5)') (STR5,I=1,(4+NTP+NQP))
C
 350     CONTINUE
         READ (IFIL,'(A5)',END=600) STR5
         IF ( STR5.NE.'*****' ) GOTO 350
C
         READ (IFIL,FMT=FMTE,END=700) ERYD,IQP
         IF ( IPRINT.GT.0 ) WRITE (6,99004) IE,ERYD,IQP
         IF ( IQ.NE.IQP .OR. ABS(ERYDINP-ERYD).GT.1D-10 ) GOTO 350
C
 400     CONTINUE
         READ (IFIL,FMT=FMTT,END=500,IOSTAT=IOS) I,J,TAUINP,MSSINP
         IF ( (I+J).NE.0 ) THEN
            TAUQ(I,J) = TAUINP
            MSSQ(I,J) = MSSINP
            IF ( (I+J).LT.2*N ) GOTO 400
         END IF
C
      END IF
C
C   ====================================================================
 500  CONTINUE
      IF ( IOS.NE.0 ) WRITE (6,99005) IOS,ERYD
C   ====================================================================
      RETURN
 600  CONTINUE
      STOP
 700  CONTINUE
      STOP 'ERYDINP NOT FOUND'
 800  CONTINUE
      NETAU = 0
      RETURN
 900  CONTINUE
      CALL STOP_MESSAGE(ROUTINE,
     &                  'inconsistent header in  TAU-file    for RDTAU')
 1000 CONTINUE
      CALL STOP_MESSAGE(ROUTINE,
     &                 'inconsistent header in  TAU-file    for RDTAUMQ'
     &                 )
C
C-----------------------------------------------------------------------
99001 FORMAT (//,10X,A30,//,A)
99002 FORMAT (//,10X,A30,//,A)
99003 FORMAT (/,10X,'reading TAU-file for:',A,/,10X,
     &        'number of tabulated energies ',I5,/,10X,
     &        'format type                  ',I5,/)
99004 FORMAT (/,10X,'energy   ',I10,2F10.5,'  IQ ',I3)
99005 FORMAT (/,'<READTAU>   EOF  (IOSTAT=',I3,')   for ERYD =',2F21.15,
     &        /)
      END
