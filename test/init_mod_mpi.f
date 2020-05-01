C*==init_mod_mpi.f    processed by SPAG 6.70Rc at 08:58 on 12 Dec 2016
      SUBROUTINE INIT_MOD_MPI
C   ********************************************************************
C   *                                                                  *
C   *  initialize all variable dealing with  MPI  mode                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI,MPI_ID,MPI_OWN_CSUM,MPI_OWN_RSUM,NPROCS,
     &    MPI_COMM_WORLD,MPI_CHARACTER,MPI_STATUS_SIZE,MPI_STATUS,
     &    MPI_TAG,MPI_OWN_ISUM
      USE MOD_FILES,ONLY:IPRINT,IFILINP,FOUND_SECTION
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IERR,LL,PROCID
      CHARACTER*80 INPFILE,OUTFILE
      CHARACTER*10 IOUT
      EXTERNAL OWN_CSUM,OWN_ISUM,OWN_RSUM
C
C*** End of declarations rewritten by SPAG
C
      CALL MPI_INIT(IERR)
C
      IF ( IERR.NE.0 ) THEN
C
C-------------------------------------------- SEQUENTIAL MODE --- NO MPI
C
         MPI = .FALSE.
         MPI_ID = 0
         NPROCS = 1
C
      ELSE
C------------------------------------------------- PARALLEL MODE --- MPI
C
         MPI = .TRUE.
C
         ALLOCATE (MPI_STATUS(MPI_STATUS_SIZE))
C
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPI_ID,IERR)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
C
         CALL MPI_OP_CREATE(OWN_ISUM,.TRUE.,MPI_OWN_ISUM,IERR)
         CALL MPI_OP_CREATE(OWN_RSUM,.TRUE.,MPI_OWN_RSUM,IERR)
         CALL MPI_OP_CREATE(OWN_CSUM,.TRUE.,MPI_OWN_CSUM,IERR)
C
         IF ( MPI_ID.EQ.0 ) THEN
C
            CALL GETARG(1,INPFILE)
C
            OPEN (IFILINP,FILE=INPFILE,STATUS='OLD')
C
            CALL INPUT_FIND_SECTION('CONTROL',0)
C
          IF ( FOUND_SECTION ) CALL SECTION_SET_INTEGER('PRINT',IPRINT,
     &           9999,0)
C
            DO PROCID = 1,NPROCS - 1
               CALL MPI_SEND(INPFILE,80,MPI_CHARACTER,PROCID,0,
     &                       MPI_COMM_WORLD,IERR)
C
               CALL DRV_MPI_SEND_IX(IPRINT,1,PROCID,MPI_TAG)
            END DO
C
         ELSE
C
            CALL MPI_RECV(INPFILE,80,MPI_CHARACTER,0,0,MPI_COMM_WORLD,
     &                    MPI_STATUS,IERR)
C
            OPEN (IFILINP,FILE=INPFILE,STATUS='OLD')
C
            CALL DRV_MPI_RECV_IX(IPRINT,1,0,MPI_TAG)
C
            IF ( IPRINT.GT.2 ) THEN
C
               IOUT = '0000'
               CALL STRING_ADD_N(IOUT,MPI_ID)
               LL = LEN_TRIM(IOUT)
               OUTFILE = 'SCRATCH.out'//IOUT((LL-3):LL)
C
               OPEN (UNIT=6,FILE=OUTFILE(1:15))
            ELSE
               OUTFILE = '/dev/null'
               OPEN (UNIT=6,FILE=OUTFILE(1:9))
            END IF
C
         END IF
C
      END IF
c modified by XJQ: nprocs==1
      if(nprocs==1) then
        mpi = .false.
        mpi_id = 0
      endif
c end-mod-xjq
C
      END
