C*==mpi_init.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_INIT(IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_INIT9
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 9999
      END
C*==mpi_comm_rank.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_COMM_RANK37
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,MPI_COMM_WORLD,MYID
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,*) MPI_COMM_WORLD,MYID,IERR
      STOP 'in MPI_COMM_RANK'
      END
C*==mpi_comm_size.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_COMM_SIZE66
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,MPI_COMM_WORLD,NPROCS
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,*) MPI_COMM_WORLD,NPROCS,IERR
      STOP 'in MPI_COMM_SIZE'
      END
C*==mpi_op_create.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_OP_CREATE(OWN_CSUM,FLAG,MPI_OWN_CSUM,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_OP_CREATE95
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FLAG
      INTEGER IERR,MPI_OWN_CSUM
      REAL OWN_CSUM
C
C Local variables
C
      COMPLEX*16 A(1),B(1)
      INTEGER N,TYPE
      EXTERNAL OWN_CSUM
C
C*** End of declarations rewritten by SPAG
C
      MPI_OWN_CSUM = 0
C
      WRITE (6,*) FLAG,MPI_OWN_CSUM,IERR
      IF ( FLAG ) WRITE (6,*) OWN_CSUM(A,B,N,TYPE)
      STOP 'in MPI_OP_CREATE'
      END
C*==mpi_reduce.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_REDUCE(A,B,I,MPI_DOUBLE_PRECISION,MPI_OWN_RSUM,J,
     &                      MPI_COMM_WORLD,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_REDUCE139
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,B
      INTEGER I,IERR,J,MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,MPI_OWN_RSUM
C
C*** End of declarations rewritten by SPAG
C
      IERR = 9999
      B = 0D0
      MPI_OWN_RSUM = 0
      WRITE (6,*) A,B,I,MPI_DOUBLE_PRECISION,MPI_OWN_RSUM,J,
     &            MPI_COMM_WORLD,IERR
      STOP 'in MPI_REDUCE'
      END
C*==mpi_send.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_SEND(STR,LSTR,MPI_CHARACTER,I,J,MPI_COMM_WORLD,
     &                    IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_SEND174
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IERR,J,LSTR,MPI_CHARACTER,MPI_COMM_WORLD
      CHARACTER*(*) STR
C
C*** End of declarations rewritten by SPAG
C
C
      WRITE (6,*) STR,LSTR,MPI_CHARACTER,I,J,MPI_COMM_WORLD,IERR
      STOP 'in MPI_SEND'
      END
C*==mpi_barrier.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_BARRIER(MPI_COMM_WORLD,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_BARRIER205
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR,MPI_COMM_WORLD
C
C*** End of declarations rewritten by SPAG
C
      IERR = 9999
      WRITE (6,*) MPI_COMM_WORLD,IERR
      STOP 'in MPI_BARRIER'
      END
C*==mpi_bcast.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_BCAST(A,I,MPI_DOUBLE_PRECISION,J,MPI_COMM_WORLD,
     &                     IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_BCAST236
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A
      INTEGER I,IERR,J,MPI_COMM_WORLD,MPI_DOUBLE_PRECISION
C
C*** End of declarations rewritten by SPAG
C
      IERR = 9999
      WRITE (6,*) A,I,MPI_DOUBLE_PRECISION,J,MPI_COMM_WORLD,IERR
      STOP 'in MPI_BCAST'
      END
C*==mpi_recv.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_RECV(STR,LSTR,MPI_CHARACTER,I,J,MPI_COMM_WORLD,
     &                    STATUS,IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_RECV268
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,IERR,J,LSTR,MPI_CHARACTER,MPI_COMM_WORLD
      CHARACTER*(*) STR
      INTEGER STATUS(1)
C
C*** End of declarations rewritten by SPAG
C
      STR = '*'
      STATUS = 0
      WRITE (6,*) STR,LSTR,MPI_CHARACTER,I,J,MPI_COMM_WORLD,STATUS,IERR
      STOP 'in MPI_RECV'
      END
C*==mpi_finalize.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_FINALIZE(IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_FINALIZE301
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
      WRITE (6,*) IERR
      STOP 'in MPI_FINALIZE'
      END
C*==mpi_gather.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MPI_GATHER(IERR)
C   ********************************************************************
C   *                                                                  *
C   *   MPI  dummy subroutine                                          *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--MPI_GATHER331
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,*) IERR
      STOP 'in MPI_GATHER'
      END
