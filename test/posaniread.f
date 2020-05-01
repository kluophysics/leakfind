C*==posaniread.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANIREAD(DATSET,LDATSET,NT,IPRINT,JRWS,IMT,TAUAPOS,
     &                      DMATTPOS,DTILTPOS,TAUABPOS,TAUQQPOS,ZGPOS,
     &                      NQMAX,NTMAX,NRMAX,NMMAX)
C    *******************************************************************
C    *                                                                 *
C    *  READ      the information on the positron state at the         *
C    *  bottom of the positron band                                    *
C    *                                                                 *
C    *  NOTE:  the imaginary part of the positron wave function g      *
C    *         is suppressed                                           *
C    *                                                                 *
C    *******************************************************************
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C*--POSANIREAD16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POSANIREAD')
C
C Dummy arguments
C
      CHARACTER*80 DATSET
      INTEGER IPRINT,LDATSET,NMMAX,NQMAX,NRMAX,NT,NTMAX
      COMPLEX*16 DMATTPOS(NTMAX),DTILTPOS(NTMAX),TAUABPOS(NTMAX,NTMAX),
     &           TAUAPOS(NTMAX),TAUQQPOS(NQMAX,NQMAX),ZGPOS(NRMAX,NTMAX)
      INTEGER IMT(NTMAX),JRWS(NMMAX)
C
C Local variables
C
      CHARACTER*80 FILNAM
      INTEGER I,IM,IQA,IQB,IT,ITA,ITB,JTOP,LFILNAM
      COMPLEX*16 MEZJ(2,2,NTMAX),MEZZ(2,2,NTMAX),ZF
      REAL*8 RDUM
C
C*** End of declarations rewritten by SPAG
C
      LFILNAM = LDATSET + 12
      FILNAM = DATSET(1:LDATSET)//'_pos_dat.pan'
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFILNAM))
C
      DO I = 1,3
         READ (IOTMP,99001,ERR=100,END=100)
      END DO
C
C-----------------------------------------------------------------------
      DO IT = 1,NT
         IM = IMT(IT)
         JTOP = JRWS(IM)
C
         READ (IOTMP,99001,ERR=100,END=100)
         READ (IOTMP,99001,ERR=100,END=100)
         READ (IOTMP,99002,ERR=100,END=100) TAUAPOS(IT)
         READ (IOTMP,99002,ERR=100,END=100) DMATTPOS(IT)
         READ (IOTMP,99002,ERR=100,END=100) DTILTPOS(IT)
         READ (IOTMP,99002,ERR=100,END=100) (MEZZ(I,I,IT),I=1,2)
         READ (IOTMP,99002,ERR=100,END=100) (MEZJ(I,I,IT),I=1,2)
         READ (IOTMP,99001,ERR=100,END=100)
         DO I = 1,JTOP
            READ (IOTMP,99003,ERR=100,END=100) RDUM,ZGPOS(I,IT),ZF
         END DO
C--- suppress imaginary part of wave function
         DO I = 1,JTOP
            ZGPOS(I,IT) = DREAL(ZGPOS(I,IT))
         END DO
C
         IF ( IPRINT.GT.3 ) THEN
            WRITE (6,99005) IT
            WRITE (6,99002) TAUAPOS(IT)
            WRITE (6,99002) (MEZZ(I,I,IT),I=1,2)
            WRITE (6,99002) (MEZJ(I,I,IT),I=1,2)
            WRITE (6,99003) RDUM,ZGPOS(1,IT),ZF
         END IF
      END DO
C
C----------------------------------------------------------------- TAUAB
C
      READ (IOTMP,99001,ERR=100,END=100)
      DO ITA = 1,NT
         DO ITB = 1,NT
            READ (IOTMP,99006) I,I,TAUABPOS(ITA,ITB)
         END DO
      END DO
C
C----------------------------------------------------------------- TAUQQ
C
      READ (IOTMP,99001,ERR=100,END=100)
      DO IQA = 1,NT
         DO IQB = 1,NT
            READ (IOTMP,99006) I,I,TAUQQPOS(IQA,IQB)
         END DO
      END DO
C
      WRITE (6,99007) FILNAM(1:LFILNAM)
C
      RETURN
C-----------------------------------------------------------------------
 100  CONTINUE
      WRITE (6,99004) FILNAM(1:LFILNAM)
      STOP
C
99001 FORMAT (A)
99002 FORMAT (/,:,(2(2E14.6,:,3x)))
99003 FORMAT (5E14.6)
99004 FORMAT (//,60('*'),/,10X,' STOP in <POSANIREAD>',/,10X,
     &        'trouble reading positron data file ',A)
99005 FORMAT ('<POSANIREAD> read for IT=',I3)
99006 FORMAT (2I3,2(2E14.6,3x))
99007 FORMAT (/,10X,'<POSANIREAD> read positron data from file ',A,/)
      END
C
