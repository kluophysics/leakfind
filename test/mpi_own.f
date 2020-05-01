C*==mpi_distribute.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE MPI_DISTRIBUTE(IPROC_XLOOP,NLOOP,MPI_XLOOP,KEY)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for  MPI system                               *
C   *                                                                  *
C   *   distribute a DO-loop with  I=1,NLOOP  over  NPROCS  processes  *
C   *   loop  X  will be done by process  IPROC = IPROC_XLOOP(I)       *
C   *   loop  X = NLOOP  will be done by  IPROC = 0 (master)           *
C   *                                                                  *
C   *   this supplies the distribution if                              *
C   *   EITHER  the  E-, the K-  OR  the M-loop is run parallel        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:NPROCS,MPI,MPI_ID
      USE MOD_CALCMODE,ONLY:TASK
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*1 KEY
      LOGICAL MPI_XLOOP
      INTEGER NLOOP
      INTEGER IPROC_XLOOP(NLOOP)
C
C Local variables
C
      INTEGER I,IPROC
C
C*** End of declarations rewritten by SPAG
C
      IPROC_XLOOP(1:NLOOP) = MPI_ID
C
C------------------------------------------------------------------------
C  NOTE: for k-space parallelisation for SIGMA_SPIN=.T. has the k-tables
C            chopped and already distributed to the individual processes
C
      IF ( (TASK.EQ.'SIGMA' .OR. TASK.EQ.'PHONONS') .AND. KEY.EQ.'K' )
     &     RETURN
C
C
C------------------------------------------------------------------------
C
      IF ( MPI_XLOOP ) THEN
C
         IF ( .NOT.MPI ) STOP '<MPI_DISTRIBUTE>  MPI = .FALSE. '
C
         IPROC = -1
         DO I = NLOOP,1, - 1
            IPROC = IPROC + 1
            IF ( IPROC.EQ.NPROCS ) IPROC = 0
            IPROC_XLOOP(I) = IPROC
         END DO
C
         CALL DRV_MPI_BARRIER
C
      END IF
C
      END
C*==own_isum.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE OWN_ISUM(A,B,N,TYPE)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for  MPI system                               *
C   *   add   INTEGER   vector  A  of length  N  to vector  B          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,TYPE
      INTEGER A(N),B(N)
C
C*** End of declarations rewritten by SPAG
C
      IF ( TYPE.NE.MPI_INTEGER ) STOP 'TYPE in <OWN_ISUM>'
C
      B(1:N) = B(1:N) + A(1:N)
C
      END
C*==own_csum.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE OWN_CSUM(A,B,N,TYPE)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for  MPI system                               *
C   *   add   COMPLEX   vector  A  of length  N  to vector  B          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_COMPLEX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER N,TYPE
      COMPLEX*16 A(N),B(N)
C
C*** End of declarations rewritten by SPAG
C
      IF ( TYPE.NE.MPI_DOUBLE_COMPLEX ) STOP 'TYPE in <OWN_CSUM>'
C
      CALL ZAXPY(N,C1,A(1),1,B(1),1)
C
      END
C*==own_rsum.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE OWN_RSUM(A,B,N,TYPE)
C   ********************************************************************
C   *                                                                  *
C   *   auxilary routine for  MPI system                               *
C   *   add    REAL     vector  A  of length  N  to vector  B          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 R1
      PARAMETER (R1=1.0D0)
C
C Dummy arguments
C
      INTEGER N,TYPE
      REAL*8 A(N),B(N)
C
C*** End of declarations rewritten by SPAG
C
      IF ( TYPE.NE.MPI_DOUBLE_PRECISION ) STOP 'TYPE in <OWN_RSUM>'
C
      CALL DAXPY(N,R1,A(1),1,B(1),1)
C
      END
C*==drv_mpi_barrier.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_BARRIER
C   ********************************************************************
C   *                                                                  *
C   *   driver to call MPI routine  MPI_BARRIER                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.MPI ) RETURN
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_BARRIER>'
C
      END
C*==drv_mpi_bcast_i.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_BCAST_I(ID_SENDER,A,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_BCAST          *
C   *   for  INTEGER  array  A  of length N                            *
C   *   sent by process  ID_SENDER  to all other processes             *
C   *                                                                  *
C   *   NOTE: A declared as SCALAR to avoid SPAG error message         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER A,ID_SENDER,N
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      CALL MPI_BCAST(A,N,MPI_INTEGER,ID_SENDER,MPI_COMM_WORLD,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_BCASTX_I>'
C
      END
C*==drv_mpi_bcast_r.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_BCAST_R(ID_SENDER,A,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_BCAST          *
C   *   for  REAL  array  A  of length N                               *
C   *   sent by process  ID_SENDER  to all other processes             *
C   *                                                                  *
C   *   NOTE: A declared as SCALAR to avoid SPAG error message         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A
      INTEGER ID_SENDER,N
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      CALL MPI_BCAST(A,N,MPI_DOUBLE_PRECISION,ID_SENDER,MPI_COMM_WORLD,
     &               IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_BCASTX_R>'
C
      END
C*==drv_mpi_bcast_c.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_BCAST_C(ID_SENDER,A,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_BCAST          *
C   *   for  COMPLEX  array  A  of length N                            *
C   *   sent by process  ID_SENDER  to all other processes             *
C   *                                                                  *
C   *   NOTE: A declared as SCALAR to avoid SPAG error message         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 A
      INTEGER ID_SENDER,N
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      CALL MPI_BCAST(A,N,MPI_DOUBLE_COMPLEX,ID_SENDER,MPI_COMM_WORLD,
     &               IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_BCASTX_CA>'
C
      END
C*==drv_mpi_bcast_l.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_BCAST_L(ID_SENDER,A,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_BCAST          *
C   *   for  LOGICAL  array  A  of length N                            *
C   *   sent by process  ID_SENDER  to all other processes             *
C   *                                                                  *
C   *   NOTE: A declared as SCALAR to avoid SPAG error message         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_LOGICAL,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL A
      INTEGER ID_SENDER,N
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      CALL MPI_BCAST(A,N,MPI_LOGICAL,ID_SENDER,MPI_COMM_WORLD,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_BCASTX_L>'
C
      END
C*==drv_mpi_send_i.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_SEND_I(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   INTEGER   array  A  of length  N                     *
C   *   from process  IPROC <> 0  to master process  IPROC = 0         *
C   *                                                                  *
C   *   MPI_RECV for IPROC = 0  and  MPI_BARRIER  are called           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER,MPI_COMM_WORLD,MPI_ID,MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      INTEGER A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      IF ( IPROC.EQ.0 ) RETURN
C
      IF ( MPI_ID.EQ.0 ) THEN
C
         CALL MPI_RECV(A,N,MPI_INTEGER,IPROC,MPI_TAG,MPI_COMM_WORLD,
     &                 MPI_STATUS,IERR)
C
      ELSE IF ( MPI_ID.EQ.IPROC ) THEN
C
         CALL MPI_SEND(A,N,MPI_INTEGER,0,MPI_TAG,MPI_COMM_WORLD,IERR)
C
      END IF
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_I>'
C
      END
C*==drv_mpi_send_r.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_SEND_R(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   REAL  array  A  of length  N                         *
C   *   from process  IPROC <> 0  to master process  IPROC = 0         *
C   *                                                                  *
C   *   MPI_RECV for IPROC = 0  and  MPI_BARRIER  are called           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,MPI_ID,
     &    MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      REAL*8 A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      IF ( IPROC.EQ.0 ) RETURN
C
      IF ( MPI_ID.EQ.0 ) THEN
C
         CALL MPI_RECV(A,N,MPI_DOUBLE_PRECISION,IPROC,MPI_TAG,
     &                 MPI_COMM_WORLD,MPI_STATUS,IERR)
C
      ELSE IF ( MPI_ID.EQ.IPROC ) THEN
C
         CALL MPI_SEND(A,N,MPI_DOUBLE_PRECISION,0,MPI_TAG,
     &                 MPI_COMM_WORLD,IERR)
C
      END IF
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_reduce_i.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_REDUCE_I(A,B,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_REDUCE         *
C   *   summing   INTEGER  array  A  of length  N  on array  B         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER,MPI_COMM_WORLD,MPI_OWN_ISUM,MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      INTEGER A(N),B(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI_ID.EQ.0 ) B(1:N) = 0
C
      CALL MPI_REDUCE(A,B,N,MPI_INTEGER,MPI_OWN_ISUM,0,MPI_COMM_WORLD,
     &                IERR)
C
      IF ( MPI_ID.EQ.0 ) A(1:N) = B(1:N)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_REDUCE_I>'
C
      END
C*==drv_mpi_reduce_r.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_REDUCE_R(A,B,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_REDUCE         *
C   *   summing   REAL  array  A  of length  N  on array  B            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,MPI_OWN_RSUM,
     &    MPI_ID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(N),B(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI_ID.EQ.0 ) B(1:N) = 0D0
C
      CALL MPI_REDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_OWN_RSUM,0,
     &                MPI_COMM_WORLD,IERR)
C
      IF ( MPI_ID.EQ.0 ) A(1:N) = B(1:N)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_REDUCE_R>'
C
      END
C*==drv_mpi_reduce_c.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_REDUCE_C(A,B,N)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_REDUCE         *
C   *   summing   COMPLEX  array  A  of length  N  on array  B         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,MPI_OWN_CSUM,
     &    MPI_ID
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N),B(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IF ( MPI_ID.EQ.0 ) B(1:N) = C0
C
      CALL MPI_REDUCE(A,B,N,MPI_DOUBLE_COMPLEX,MPI_OWN_CSUM,0,
     &                MPI_COMM_WORLD,IERR)
C
      IF ( MPI_ID.EQ.0 ) A(1:N) = B(1:N)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_REDUCE_C>'
C
      END
C*==drv_mpi_send_rx.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_SEND_RX(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   REAL  array  A  of length  N                         *
C   *   from calling process to process IPROC                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      REAL*8 A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      CALL MPI_SEND(A,N,MPI_DOUBLE_PRECISION,IPROC,MPI_TAG,
     &              MPI_COMM_WORLD,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_recv_rx.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_RECV_RX(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_RECV           *
C   *   copying   REAL  array  A  of length  N                         *
C   *   from process IPROC  to calling process                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      REAL*8 A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      CALL MPI_RECV(A,N,MPI_DOUBLE_PRECISION,IPROC,MPI_TAG,
     &              MPI_COMM_WORLD,MPI_STATUS,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_send_ix.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_SEND_IX(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   INTEGER  array  A  of length  N                      *
C   *   from calling process to process IPROC                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_COMM_WORLD,MPI_INTEGER
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      INTEGER A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      CALL MPI_SEND(A,N,MPI_INTEGER,IPROC,MPI_TAG,MPI_COMM_WORLD,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_recv_ix.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_RECV_IX(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_RECV           *
C   *   copying  INTEGER array  A  of length  N                        *
C   *   from process IPROC  to calling process                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER,MPI_COMM_WORLD,MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      INTEGER A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      CALL MPI_RECV(A,N,MPI_INTEGER,IPROC,MPI_TAG,MPI_COMM_WORLD,
     &              MPI_STATUS,IERR)
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_recv_r.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_RECV_R(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   REAL  array  A  of length  N                         *
C   *   from master process  IPROC = 0  to process  IPROC <> 0         *
C   *                                                                  *
C   *   MPI_SEND for IPROC = 0  called                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,MPI_ID,
     &    MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      REAL*8 A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      IF ( IPROC.EQ.0 ) RETURN
C
      IF ( MPI_ID.EQ.IPROC ) THEN
C
         CALL MPI_RECV(A,N,MPI_DOUBLE_PRECISION,0,MPI_TAG,
     &                 MPI_COMM_WORLD,MPI_STATUS,IERR)
C
      ELSE IF ( MPI_ID.EQ.0 ) THEN
C
         CALL MPI_SEND(A,N,MPI_DOUBLE_PRECISION,IPROC,MPI_TAG,
     &                 MPI_COMM_WORLD,IERR)
C
      END IF
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_R>'
C
      END
C*==drv_mpi_recv_i.f    processed by SPAG 6.70Rc at 18:58 on 29 Dec 2014
      SUBROUTINE DRV_MPI_RECV_I(A,N,IPROC,MPI_TAG)
C   ********************************************************************
C   *                                                                  *
C   *   driver to call general purpose MPI routine  MPI_SEND           *
C   *   copying   INTEGER   array  A  of length  N                     *
C   *   from master process  IPROC = 0  to process  IPROC <> 0         *
C   *                                                                  *
C   *   MPI_SEND for IPROC = 0  called                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_MPI,ONLY:MPI_INTEGER,MPI_COMM_WORLD,MPI_ID,MPI_STATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPROC,MPI_TAG,N
      INTEGER A(N)
C
C Local variables
C
      INTEGER IERR
C
C*** End of declarations rewritten by SPAG
C
      IERR = 0
C
      IF ( IPROC.EQ.0 ) RETURN
C
      IF ( MPI_ID.EQ.IPROC ) THEN
C
         CALL MPI_RECV(A,N,MPI_INTEGER,IPROC,MPI_TAG,MPI_COMM_WORLD,
     &                 MPI_STATUS,IERR)
C
      ELSE IF ( MPI_ID.EQ.0 ) THEN
C
         CALL MPI_SEND(A,N,MPI_INTEGER,0,MPI_TAG,MPI_COMM_WORLD,IERR)
C
      END IF
C
      IF ( IERR.NE.0 ) STOP 'in <DRV_MPI_SEND_I>'
C
      END
