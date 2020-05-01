C*==stop_trace_back.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE STOP_TRACE_BACK(ROUTINE,MESSAGE)
C   ********************************************************************
C   *                                                                  *
C   *   stop program execution and print   MESSAGE                     *
C   *   ROUTINE is the name of the calling routine                     *
C   *                                                                  *
C   *   NOTE: a division by 0 is done to force a trace back            *
C   *         provided the program was compiled with test flags on     *
C   *         this should help to locate the problem unambigously      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      IMPLICIT NONE
C*--STOP_TRACE_BACK198
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) MESSAGE,ROUTINE
C
C Local variables
C
      INTEGER LM,LR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI ) WRITE (6,99003) MPI_ID,NPROCS
C
      LR = LEN_TRIM(ROUTINE)
      LM = LEN_TRIM(MESSAGE)
C
      WRITE (6,99001) ROUTINE(1:LR),MESSAGE(1:LM)
C
      WRITE (6,99002)
      LR = 0
      WRITE (6,*) 'division by 0:',NPROCS/LR
C
      CALL FLUSH(6)
C
      STOP ' via <STOP_TRACE_BACK>'
C
99001 FORMAT (2(/,1X,79('#')),//,5X,'STOP in subroutine <',A,'>',//,10X,
     &        'with error message:  ',A,//,2(1X,79('#'),/))
99002 FORMAT (//,2(/,1X,79('T')),//,5X,
     &        'forcing TRACE BACK by dividing by 0',//,2(1X,79('T'),/),
     &        //)
99003 FORMAT (2(/,20(' MPI')),/,5X,'MPI process number',I3,'  out of ',
     &        I3,' processes',//,2(20(' MPI'),/))
      END
C*==stop_error.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE STOP_ERROR(ROUTINE,MESSAGE)
C   ********************************************************************
C   *                                                                  *
C   *   stop program execution and print   MESSAGE                     *
C   *   ROUTINE is the name of the calling routine                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      IMPLICIT NONE
C*--STOP_ERROR256
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) MESSAGE,ROUTINE
C
C Local variables
C
      INTEGER LM,LR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI ) WRITE (6,99002) MPI_ID,NPROCS
C
      LR = LEN_TRIM(ROUTINE)
      LM = LEN_TRIM(MESSAGE)
C
      WRITE (6,99001) ROUTINE(1:LR),MESSAGE(1:LM)
C
      CALL FLUSH(6)
      STOP ' via <STOP_ERROR>'
C
99001 FORMAT (2(/,1X,79('#')),//,5X,'STOP in subroutine <',A,'>',//,10X,
     &        'with error message:  ',A,//,2(1X,79('#'),/))
99002 FORMAT (2(/,20(' MPI')),/,5X,'MPI process number',I3,'  out of ',
     &        I3,' processes',//,2(20(' MPI'),/))
      END
C*==stop_regular.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE STOP_REGULAR(ROUTINE,MESSAGE)
C   ********************************************************************
C   *                                                                  *
C   *   stop program execution and print   MESSAGE                     *
C   *   ROUTINE is the name of the calling routine                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:CPU_TIME_PROGRAM_START,WALL_TIME_PROGRAM_START
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      USE MOD_CALCMODE,ONLY:TASK
      IMPLICIT NONE
C*--STOP_REGULAR308
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) MESSAGE,ROUTINE
C
C Local variables
C
      INTEGER IERR,LM,LR,RATE,WALL_TIME
      REAL*8 TIME
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI ) WRITE (6,99001) MPI_ID,NPROCS
C
      LR = LEN_TRIM(ROUTINE)
      LM = LEN_TRIM(MESSAGE)
C
      IF ( LM.EQ.0 ) THEN
         WRITE (6,99002) TASK,ROUTINE(1:LR)
      ELSE
         WRITE (6,99003) TASK,ROUTINE(1:LR),MESSAGE(1:LM)
      END IF
C
      CALL CPU_TIME(TIME)
      CALL SYSTEM_CLOCK(WALL_TIME,RATE)
C
      CALL FLUSH(6)
C
      WRITE (6,99004) TIME - CPU_TIME_PROGRAM_START,
     &                DBLE(WALL_TIME-WALL_TIME_PROGRAM_START)/DBLE(RATE)
C
      IF ( MPI ) CALL MPI_FINALIZE(IERR)
C
      CALL FLUSH(6)
      STOP
C
99001 FORMAT (2(/,20(' MPI')),//,5X,'MPI process number',I3,'  out of ',
     &        I3,' processes',//,2(20(' MPI'),/))
99002 FORMAT (2(/,1X,79('*')),//,10X,'calculation for task = ',A,//,10X,
     &        'finished in subroutine <',A,'>')
99003 FORMAT (2(/,1X,79('*')),//,10X,'calculation for task = ',A,//,10X,
     &        'finished in subroutine <',A,'>',//,10X,A)
99004 FORMAT (/,10X,'run time info  CPU ',F10.3,' sec',/,25X,'WALL',
     &        F10.3,' sec',//,2(1X,79('*'),/),/,10X,
     &        'program stopped regularly via STOP_REGULAR',/)
      END
C*==track_info.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE TRACK_INFO(ROUTINE)
C   ********************************************************************
C   *                                                                  *
C   *  print info to track run of the program                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IPRINT,WR_TRACK_INFO,CPU_TIME_PROGRAM_START,
     &    CPU_TIME_LAST_CALL,WALL_TIME_PROGRAM_START,
     &    WALL_TIME_LAST_CALL,IOTMP
      USE MOD_TYPES,ONLY:NT,NTMAX,ITBOT,ITTOP
      USE MOD_SITES,ONLY:NQ,NQMAX,NQHOST,NQCLU,NQ_L,NQ_I,NQ_R,IQBOT,
     &    IQTOP
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_SCF,ONLY:SCFSTATUS,SCFSTATUS_CLU
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_CALCMODE,ONLY:IREL,LLOYD,DMFT,LDAU,NONMAG,KMROT
      USE MOD_ANGMOM,ONLY:NL,NLMAX,NLM,NKM
      USE MOD_MPI,ONLY:MPI,MPI_ID,NPROCS
      IMPLICIT NONE
C*--TRACK_INFO388
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) ROUTINE
C
C Local variables
C
      LOGICAL FILE_IS_OPEN
      INTEGER RATE,WALL_TIME
      REAL*8 TIME
C
C*** End of declarations rewritten by SPAG
C
      IF ( IPRINT.LE.2 .AND. .NOT.WR_TRACK_INFO ) RETURN
C
      CALL CPU_TIME(TIME)
      CALL SYSTEM_CLOCK(WALL_TIME,RATE)
C
      WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),SYSTEM_TYPE,
     &                SYSTEM_DIMENSION,SUB_SYSTEM,SCFSTATUS,
     &                SCFSTATUS_CLU
      WRITE (6,99002) IREL,KMROT,FULLPOT,LLOYD,NONMAG,DMFT,LDAU
      WRITE (6,99003) NQ_L,NQ_I,NQ_R,NQHOST,NQCLU,NQ,NQMAX,IQBOT,IQTOP
      WRITE (6,99004) NT,NTMAX,ITBOT,ITTOP,NL,NLMAX,NLM,NKM
C
C-------------------------------------------- check whether file is open
C
      INQUIRE (IOTMP,OPENED=FILE_IS_OPEN)
C
      WRITE (6,99005) IOTMP,FILE_IS_OPEN
C
C-------------------------------------------------- report on MPI status
C
      IF ( MPI ) WRITE (6,99006) MPI_ID,NPROCS
C
      IF ( IPRINT.GT.5 ) WRITE (6,99007) TIME - CPU_TIME_PROGRAM_START,
     &                          TIME - CPU_TIME_LAST_CALL,
     &                          DBLE(WALL_TIME-WALL_TIME_PROGRAM_START)
     &                          /DBLE(RATE),
     &                          DBLE(WALL_TIME-WALL_TIME_LAST_CALL)
     &                          /DBLE(RATE)
C
      CALL CPU_TIME(CPU_TIME_LAST_CALL)
      CALL SYSTEM_CLOCK(WALL_TIME_LAST_CALL,RATE)
C
      CALL FLUSH(6)
C
99001 FORMAT (/,'TRACK INFO ',/,'TRACK INFO ',69('i'),/'TRACK INFO ',/,
     &        'TRACK INFO    entering subroutine  ',A,/,'TRACK INFO ',/,
     &        'TRACK INFO ',69('i'),/,'TRACK INFO',/,
     &        'TRACK INFO     SYSTYPE:    ',A,'  SYSDIM: ',A,
     &        '  SUBSYSTEM: ',A,/,'TRACK INFO     SCFSTATUS:  ',A,
     &        '  SCFSTATUS_CLU: ',A)
99002 FORMAT ('TRACK INFO     IREL: ',I4,'  KMROT: ',I3,'  FULLPOT: ',
     &        L1,'  LLOYD:   ',L1,/,'TRACK INFO     NONMAG:  ',L1,
     &        '  DMFT:    ',L1,'  LDAU:    ',L1)
99003 FORMAT ('TRACK INFO     NQ_L: ',I4,'  NQ_I: ',I4,'  NQ_R: ',I4,
     &        '  NQHOST:',I3,'  NQCLU:',I4,/,'TRACK INFO     NQ:   ',I4,
     &        '  NQMAX:',I4,'  IQBOT:',I4,'  IQTOP:',I4)
99004 FORMAT ('TRACK INFO     NT:   ',I4,'  NTMAX:',I4,'  ITBOT:',I4,
     &        '  ITTOP:',I4,/,'TRACK INFO     NL:   ',I4,'  NLMAX:',I4,
     &        '  NLM:  ',I4,'  NKM:  ',I4,/,'TRACK INFO    ')
99005 FORMAT ('TRACK INFO     temporary file IOTMP =',I3,'  is open:',
     &        L3,/,'TRACK INFO    ')
99006 FORMAT ('TRACK INFO     MPI process number',I3,' out of ',I3,
     &        ' processes',/,'TRACK INFO    ')
99007 FORMAT ('TRACK INFO     time',5X,'CPU_TOT',3X,'CPU_LAST',4X,
     &        'WALL_TOT',3X,'WALL_LAST',/,'TRACK INFO     sec',3X,F10.3,
     &        1X,F10.3,2X,F10.3,2X,F10.3,/,'TRACK INFO    ')
C
      END
C*==stop_message.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE STOP_MESSAGE(ROUTINE,MESSAGE)
C   ********************************************************************
C   *                                                                  *
C   *  to be replaced by *ERROR or *TRACE_BACK                         *
C   *                                                                  *
C   ********************************************************************
      USE MOD_MPI,ONLY:MPI,NPROCS,MPI_ID
      IMPLICIT NONE
C*--STOP_MESSAGE482
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) MESSAGE,ROUTINE
C
C Local variables
C
      INTEGER LM,LR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI ) WRITE (6,99003) MPI_ID,NPROCS
C
      LR = LEN_TRIM(ROUTINE)
      LM = LEN_TRIM(MESSAGE)
C
      WRITE (6,99001) ROUTINE(1:LR),MESSAGE(1:LM)
C
      WRITE (6,99002)
      LR = 0
      WRITE (6,*) 'division by 0:',NPROCS/LR
C
      CALL FLUSH(6)
      STOP ' via <STOP_MESSAGE>'
C
99001 FORMAT (2(/,1X,79('#')),//,5X,'STOP in subroutine <',A,'>',//,10X,
     &        'with error message:  ',A,//,2(1X,79('#'),/))
99002 FORMAT (//,2(/,1X,79('T')),//,5X,
     &        'forcing TRACE BACK by dividing by 0',//,2(1X,79('T'),/),
     &        //)
99003 FORMAT (2(/,20(' MPI')),/,5X,'MPI process number',I3,'  out of ',
     &        I3,' processes',//,2(20(' MPI'),/))
      END
C*==info_message.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE INFO_MESSAGE(ROUTINE,MESSAGE)
C   ********************************************************************
C   *                                                                  *
C   *  write important INFO to standard OUTPUT                         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--INFO_MESSAGE538
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) MESSAGE,ROUTINE
C
C Local variables
C
      INTEGER LM,LR
C
C*** End of declarations rewritten by SPAG
C
      LR = LEN_TRIM(ROUTINE)
      LM = LEN_TRIM(MESSAGE)
C
      WRITE (6,99001) ROUTINE(1:LR),MESSAGE(1:LM)
C
      RETURN
C
99001 FORMAT (/,3(/,1X,79('*')),//,5X,'INFO from subroutine <',A,'>',//,
     &        5X,A,//,3(1X,79('*'),/),/)
      END
C*==under_construction.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE UNDER_CONSTRUCTION(ROUTINE)
C   ********************************************************************
C   *                                                                  *
C   *   routine to give a WARNING from sections of the code            *
C   *   that are under construction and not fully checked              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--UNDER_CONSTRUCTION583
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) ROUTINE
C
C Local variables
C
      INTEGER L
C
C*** End of declarations rewritten by SPAG
C
      L = LEN_TRIM(ROUTINE)
C
      WRITE (6,99001) ROUTINE(1:L)
C
99001 FORMAT (2(/,1X,79('#')),/,10X,'warning from subroutine <',A,'>',/,
     &        2(1X,79('#'),/),/,10X,
     &        'this routine is under construction ',/,10X,
     &        'the results may be unreliable ',///,5X,31X,'******',/,5X,
     &        29X,'**********',/,5X,27X,'***',8X,'***',/,5X,25X,'***',
     &        4X,'*****',3X,'***',/,5X,23X,'***',5X,'*******',4X,'***',
     &        /,5X,21X,'***',6X,'*********',5X,'***',/,5X,19X,'***',8X,
     &        '*********',7X,'***',/,5X,17X,'***',11X,'*******',10X,
     &        '***',/,5X,15X,'***',3X,'****',8X,'****',13X,'***',/,5X,
     &        13X,'***',3X,'***************',18X,'***',/,5X,11X,'***',
     &        3X,'*****************',20X,'***',/,5X,9X,'***',3X,'**',4X,
     &        '***********',24X,'***',/,5X,7X,'***',7X,'**',2X,
     &        '************',25X,'***',/,5X,5X,'***',11X,'**********',
     &        3X,'**',26X,'***',/,5X,3X,'***',13X,'*********',4X,'**',
     &        28X,'***',/,5X,' ***',14X,'*********',5X,'**',30X,'***',/,
     &        5X,'***',14X,'***********',4X,'**',31X,'***',/,5X,'***',
     &        14X,'*********',2X,'**',2X,'**',31X,'***',/,5X,'***',14X,
     &        '**********',3X,'****',31X,'***',/,5X,'***',14X,
     &        '***** *****',4X,'**',31X,'***',/,5X,' ***',13X,'*****',
     &        2X,'*****',5X,'**',28X,'***',/,5X,3X,'***',11X,'*****',3X,
     &        '*****',7X,'**',23X,'***',/,5X,5X,'***',9X,'*****',3X,
     &        '*****',8X,'***',19X,'***',/,5X,7X,'***',7X,'*****',3X,
     &        '*****',7X,'*******',14X,'***',/,5X,9X,'***',5X,'*****',
     &        3X,'*****',6X,'********',12X,'***',/,5X,11X,'***',3X,
     &        '*****',3X,'**',9X,'**********',8X,'***',/,5X,13X,'***',
     &        3X,'***',3X,'*****',4X,'**************',4X,'***',/,5X,15X,
     &        '***',3X,'*',3X,'*****',3X,'**************** ***',/,5X,
     &        17X,'***',28X,'***',/,5X,19X,'***',24X,'***',/,5X,21X,
     &        '***',20X,'***',/,5X,23X,'***',16X,'***',/,5X,25X,'***',
     &        12X,'***',/,5X,27X,'***',8X,'***',/,5X,29X,'**********',/,
     &        5X,31X,'******',///)
      END
