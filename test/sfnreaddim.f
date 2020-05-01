C*==sfnreaddim.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNREADDIM(NPANLIM,NRSFLIM,SFN_AVAILABLE)
C   ********************************************************************
C   *                                                                  *
C   *      read   SHAPE FUNCTION  data  to fix array sizes             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:STORE_SFN
      USE MOD_RMESH,ONLY:NM,NMMAX,NSFMAX
      USE MOD_FILES,ONLY:SFNFIL,LSFNFIL,IOTMP,DATSET0,RDUMMY,IDUMMY,
     &    FOUND_SECTION,FOUND_STRING
      IMPLICIT NONE
C*--SFNREADDIM13
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NPANLIM,NRSFLIM
      LOGICAL SFN_AVAILABLE
C
C Local variables
C
      INTEGER IM,IPAN1,ISF,J,N,NMSFN,NPAN,NRSFTOT,NSF,NSFLIM
C
C*** End of declarations rewritten by SPAG
C
      LSFNFIL = 0
      CALL INPUT_FIND_SECTION('CONTROL',0)
      IF ( FOUND_SECTION ) CALL SECTION_SET_STRING('SFNFIL',SFNFIL,
     &     '9999',0)
      IF ( .NOT.FOUND_STRING ) THEN
         LSFNFIL = LEN_TRIM(DATSET0)
         SFNFIL = DATSET0(1:LSFNFIL)//'.sfn'
         LSFNFIL = LSFNFIL + 4
      END IF
C
C-----------------------------------------------------------------------
      IF ( .NOT.STORE_SFN ) THEN
         SFN_AVAILABLE = .FALSE.
         WRITE (6,99007)
         RETURN
      END IF
C
C-----------------------------------------------------------------------
C
      IF ( SFN_AVAILABLE ) STOP
C
      INQUIRE (FILE=SFNFIL(1:LSFNFIL),EXIST=SFN_AVAILABLE)
C------------------------------------- shape function file not available
      IF ( .NOT.SFN_AVAILABLE ) THEN
         WRITE (6,99006) SFNFIL(1:LSFNFIL)
         RETURN
      END IF
C
      OPEN (IOTMP,FILE=SFNFIL(1:LSFNFIL),ERR=100)
C
      READ (IOTMP,FMT=99003,END=100) NMSFN
C
      DO IM = 1,NM
         READ (IOTMP,FMT=99004,END=100) RDUMMY
      END DO
C
      NRSFLIM = 0
      NPANLIM = 0
      NSFLIM = 0
      DO IM = 1,NM
C
         READ (IOTMP,FMT=99003) NPAN,NRSFTOT
C
C ---- the 1st panel has to be added, for it contains no shape functions
C
         NPAN = NPAN + 1
C
         READ (IOTMP,FMT=99003) (IDUMMY,IPAN1=2,NPAN)
         READ (IOTMP,FMT=99004) (RDUMMY,RDUMMY,J=1,NRSFTOT)
         READ (IOTMP,FMT=99003) NSF
C
         DO ISF = 1,NSF
            READ (IOTMP,FMT=99003) IDUMMY
            READ (IOTMP,FMT=99004) (RDUMMY,N=1,NRSFTOT)
         END DO
C
         NRSFLIM = MAX(NRSFLIM,NRSFTOT)
         NPANLIM = MAX(NPANLIM,NPAN)
         NSFLIM = MAX(NSFLIM,NSF)
      END DO
C
      CLOSE (IOTMP)
C
C ======================================================================
C
      IF ( NM.NE.NMSFN .OR. NMSFN.GT.NMMAX .OR. NSFMAX.LT.NSFLIM ) THEN
         SFN_AVAILABLE = .FALSE.
         WRITE (6,99005) SFNFIL(1:LSFNFIL)
         WRITE (6,99001) NMMAX,NM,NMSFN,NSFMAX,NSFLIM
         RETURN
      END IF
C
      WRITE (6,99002) NPANLIM,NRSFLIM
      RETURN
C
C----------------------------------- trouble reading shape function file
 100  CONTINUE
      WRITE (6,99008) SFNFIL(1:LSFNFIL)
      STOP
99001 FORMAT (/,10X,'number of radial meshes ',/,10X,
     &        'dimension fixed               NMMAX: ',I5,/,10X,
     &        'expected                      NM:    ',I5,/,10X,
     &        'found in shape function file  NMSFN: ',I5,//,10X,
     &        'number of shape functions ',/,10X,
     &        'dimension fixed               NSFMAX:',I5,/,10X,
     &        'found in shape function file  NSFLIM:',I5,/)
99002 FORMAT (/,1X,79('*'),/,35X,'<SFNREADDIM>',/,22X,
     &        'array sizes for the shape functions',/,1X,79('*'),//,10X,
     &        'NPANLIM =',I4,5X,'NRSFLIM =',I4,//)
99003 FORMAT (16I5)
99004 FORMAT (4D20.12)
99005 FORMAT (2(/,1X,79('#')),//,34X,'<SFNREADDIM>',//,10X,
     &        'scanning the shape function file    ',A,/,10X,
     &        'inconsistencies found with internal settings',/,10X,
     &        'the shape functions will be recalculated '//,
     &        2(1X,79('#'),/))
99006 FORMAT (2(/,1X,79('#')),//,34X,'<SFNREADDIM>',/,14X,
     &        'shape function file not available:  ',A,/,2(1X,79('#'),/)
     &        )
99007 FORMAT (1(/,1X,79('I')),//,34X,'<SFNREADDIM>',/,10X,
     &        'STORE_SFN = .FALSE.  ',/,10X,
     &        'shape function will be created and not stored',/,
     &        1(1X,79('I'),/))
99008 FORMAT (2(/,1X,79('#')),//,34X,'<SFNREADDIM>',/,14X,
     &        'trouble reading shape function file ',A,/,2(1X,79('#'),/)
     &        )
      END
