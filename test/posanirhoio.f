C*==posanirhoio.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANIRHOIO(MODE,RHOPOS_ON_DISK)
C   ********************************************************************
C   *                                                                  *
C   *  read / write charge density for the positron                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRWS,FULLPOT,JRCRI
      USE MOD_TYPES,ONLY:NLMFPT,KLMFP,IMT,RHO2NS,RHOCHR,NT
      USE MOD_FILES,ONLY:LDATSET,DATSET,IOTMP
      IMPLICIT NONE
C*--POSANIRHOIO13
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*5 MODE
      LOGICAL RHOPOS_ON_DISK
C
C Local variables
C
      CHARACTER*80 FILNAM
      CHARACTER*40 FMT07
      INTEGER IM,IR,IRTOP,IT,LFILNAM,LM
      CHARACTER*10 SDUM
C
C*** End of declarations rewritten by SPAG
C
      FMT07 = '(1P,5E16.9)'
C
      WRITE (6,99001)
C
      FILNAM = DATSET(1:LDATSET)//'_positron.rho'
      LFILNAM = LDATSET + 13
      RHOPOS_ON_DISK = .FALSE.
C
C***********************************************************************
C                                READ
C***********************************************************************
C
      IF ( MODE.EQ.'READ ' ) THEN
C
         OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),ERR=100)
C
         DO IT = 1,NT
            IM = IMT(IT)
C
            IF ( .NOT.FULLPOT ) THEN
               IRTOP = JRWS(IM)
               READ (IOTMP,FMT=FMT07,ERR=100,END=100)
     &               (RHOCHR(IR,IT),IR=1,IRTOP)
            ELSE
               IRTOP = JRCRI(IM)
               DO LM = 1,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     READ (IOTMP,*,ERR=100,END=100) SDUM
                     READ (IOTMP,FMT=FMT07,ERR=100,END=100)
     &                     (RHO2NS(IR,LM,IT,1),IR=1,IRTOP)
                  END IF
               END DO
            END IF
C
         END DO
C
         WRITE (6,99002) 'read from: ',FILNAM(1:LFILNAM)
         RHOPOS_ON_DISK = .TRUE.
         IF ( NT.LT.0 ) WRITE (*,*) SDUM
C
C***********************************************************************
C                                WRITE
C***********************************************************************
C
      ELSE IF ( MODE.EQ.'WRITE' ) THEN
C
         OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),ERR=200)
C
         DO IT = 1,NT
            IM = IMT(IT)
C
            IF ( .NOT.FULLPOT ) THEN
               IRTOP = JRWS(IM)
               WRITE (IOTMP,FMT=FMT07,ERR=200)
     &                (RHOCHR(IR,IT),IR=1,IRTOP)
            ELSE
               IRTOP = JRCRI(IM)
               DO LM = 1,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     WRITE (IOTMP,'(A,2I3)',ERR=200) 'IT LM:',IT,LM
                     WRITE (IOTMP,FMT=FMT07,ERR=200)
     &                      (RHO2NS(IR,LM,IT,1),IR=1,IRTOP)
                  END IF
               END DO
            END IF
C
         END DO
C
         WRITE (6,99002) 'written to: ',FILNAM(1:LFILNAM)
         RHOPOS_ON_DISK = .TRUE.
C
      ELSE
         STOP 'in <POSANIRHOIO>:  MODE not found'
      END IF
C***********************************************************************
C
      CLOSE (IOTMP)
      RETURN
C
C=======================================================================
 100  CONTINUE
      CLOSE (IOTMP)
      WRITE (6,99003) 'read'
      RETURN
C
C=======================================================================
 200  CONTINUE
      CLOSE (IOTMP)
      WRITE (6,99003) 'write'
      RETURN
C
99001 FORMAT (/,1X,79('*'),/,33X,'<POSANIRHOIO>',/,79('*'),/)
99002 FORMAT (/,10X,'positron charge density ',A,/)
99003 FORMAT (/,10X,'<POSANIRHOIO> failed to ',A,
     &        ' positron charge density',//)
      END
