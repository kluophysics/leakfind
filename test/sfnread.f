C*==sfnread.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNREAD(NMLOC,NRPAN,XRSF,DRSF,SCLM)
C   ********************************************************************
C   *                                                                  *
C   *             read   SHAPE FUNCTION  data                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NMMAX,NPANMAX,NSFMAX,NRSFMAX,NLMSFMAX,NRSFTOT,
     &    NPAN,FLMSF,ISFLM,KLMSF,NSF,LMISF
      USE MOD_FILES,ONLY:SFNFIL,LSFNFIL,IOTMP
      IMPLICIT NONE
C*--SFNREAD12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NMLOC
      REAL*8 DRSF(NRSFMAX,NMMAX),SCLM(NMMAX),XRSF(NRSFMAX,NMMAX)
      INTEGER NRPAN(NPANMAX,NMMAX)
C
C Local variables
C
      INTEGER IM,IPAN1,ISF,ISUM,J,LM,N,NFUN,NMSFN
      LOGICAL SFN_AVAILABLE
C
C*** End of declarations rewritten by SPAG
C
      INQUIRE (FILE=SFNFIL,EXIST=SFN_AVAILABLE)
C
      IF ( SFN_AVAILABLE ) THEN
C
         OPEN (IOTMP,FILE=SFNFIL,ERR=100)
C
         READ (IOTMP,FMT=99002,END=100) NMSFN
C
         IF ( NMLOC.NE.NMSFN .OR. NMSFN.GT.NMMAX ) THEN
            WRITE (6,*) '  NMLOC: ',NMLOC,' NMSFN: ',NMSFN,' NMMAX: ',
     &                  NMMAX
            STOP 'in <SFNREAD>:  NMLOC != NMSFN'
         END IF
C
         DO IM = 1,NMLOC
            READ (IOTMP,FMT=99003,END=100) SCLM(IM)
         END DO
C
         DO IM = 1,NMLOC
            READ (IOTMP,FMT=99002) NPAN(IM),NRSFTOT(IM)
            IF ( NPAN(IM).GT.NPANMAX ) THEN
               WRITE (6,99001) IM,'NPAN(IM)',NPAN(IM),'  > NPANMAX',
     &                         NPANMAX
               STOP
            END IF
C
C ---- the 1st panel has to be added, for it contains no shape functions
C
            NPAN(IM) = NPAN(IM) + 1
C
            READ (IOTMP,FMT=99002) (NRPAN(IPAN1,IM),IPAN1=2,NPAN(IM))
            READ (IOTMP,FMT=99003) (XRSF(J,IM),DRSF(J,IM),J=1,NRSFTOT(IM
     &                             ))
            ISUM = 0
            DO IPAN1 = 2,NPAN(IM)
               ISUM = ISUM + NRPAN(IPAN1,IM)
            END DO
            IF ( ISUM.NE.NRSFTOT(IM) .OR. NRSFTOT(IM).GT.NRSFMAX ) THEN
               WRITE (6,99001) IM,'SUM(NRPAN)',ISUM,'  =?  NRSFTOT',
     &                         NRSFTOT(IM),'NRSFTOT',NRSFTOT(IM),
     &                         '  =?  NRSFMAX: ',NRSFMAX
               STOP
            END IF
C
            READ (IOTMP,FMT=99002) NSF(IM)
            NFUN = NSF(IM)
C
            IF ( NFUN.GT.NSFMAX ) THEN
               WRITE (6,99001) IM,'NSF',NSF(IM),'  >  NSFMAX',NSFMAX
               STOP
            END IF
C
            DO LM = 1,NLMSFMAX
               KLMSF(LM,IM) = 0
            END DO
C
            DO ISF = 1,NFUN
               READ (IOTMP,FMT=99002) LM
               LMISF(ISF,IM) = LM
               KLMSF(LM,IM) = 1
               ISFLM(LM,IM) = ISF
               READ (IOTMP,FMT=99003) (FLMSF(N,ISF,IM),N=1,NRSFTOT(IM))
            END DO
         END DO
C
         CLOSE (IOTMP)
C
         WRITE (6,99004) SFNFIL(1:LSFNFIL),1,NMLOC
C
C ======================================================================
         RETURN
      END IF
 100  CONTINUE
      WRITE (6,99005) SFNFIL(1:LSFNFIL)
      STOP
C
99001 FORMAT (/,' ##### TROUBLE in <SFNREAD> ',52('#'),/,10X,'MESH  ',
     &        I4,:,(/,2(A,I4,3X)))
99002 FORMAT (16I5)
99003 FORMAT (4D20.12)
99004 FORMAT (/,10X,'the shape functions have been read successfully',/,
     &        10X,'from file: ',A,/,10X,'for meshes IM =',I4,'  --',I3)
99005 FORMAT (2(/,1X,79('#')),//,34X,'<SFNREAD>',/,14X,
     &        'trouble reading shape function file ',A,/,2(1X,79('#'),/)
     &        )
      END
C*==sfndump.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNDUMP
C   ********************************************************************
C   *                                                                  *
C   *             DUMP   SHAPE FUNCTION  data                          *
C   *                                                                  *
C   *             in case of BREAKPOINT = 1                            *
C   *                                                                  *
C   *             called in <SCF>                                      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NLMSFMAX,NRSFTOT,FLMSF,ISFLM,KLMSF,NSF,LMISF,NM
      USE MOD_FILES,ONLY:SFNFIL,LSFNFIL,IOTMP
      USE MOD_CALCMODE,ONLY:ROTATE_LATTICE
      IMPLICIT NONE
C*--SFNDUMP141
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNDUMP')
C
C Local variables
C
      INTEGER IM,ISF,NRSF
C
C*** End of declarations rewritten by SPAG
C
      IF ( ROTATE_LATTICE ) THEN
         SFNFIL = 'sfn_global.dat'
         LSFNFIL = 14
      ELSE
         SFNFIL = 'sfn_local.dat'
         LSFNFIL = 13
      END IF
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,SFNFIL(1:LSFNFIL))
C
      DO IM = 1,NM
C
         NRSF = NRSFTOT(IM)
C
         WRITE (IOTMP,99001) 'NRSF      ',NRSF
         WRITE (IOTMP,99001) 'NSF       ',NSF(IM)
         WRITE (IOTMP,99001) 'KLMSF     ',KLMSF(1:NLMSFMAX,IM)
         WRITE (IOTMP,99001) 'ISFLM     ',ISFLM(1:NLMSFMAX,IM)
C
         DO ISF = 1,NSF(IM)
            WRITE (IOTMP,99001) 'LMISF     ',LMISF(ISF,IM)
C
            WRITE (IOTMP,99002) FLMSF(1:NRSF,ISF,IM)
         END DO
      END DO
C
      WRITE (6,*) ' SHAPE FUNCTIONS DUMPED TO ',SFNFIL(1:LSFNFIL)
C
      CALL STOP_BREAKPOINT(ROUTINE)
C
99001 FORMAT (A,12I5,:,/,(:,10X,12I5))
99002 FORMAT (4D20.12)
      END
