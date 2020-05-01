C*==mod_mpi.f    processed by SPAG 6.70Rc at 07:51 on  2 Jul 2012
      MODULE MOD_MPI
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables connected with   MPI  mode        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
      INCLUDE 'mpif.h'
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      LOGICAL MPI,MPI_ELOOP,MPI_KLOOP
      INTEGER MPI_ID,MPI_OWN_CSUM,MPI_OWN_ISUM,MPI_OWN_RSUM,
     &        MPI_STATUS(:),NPROCS,IPROC_K(:),IPROC_E(:)
      SAVE MPI,MPI_STATUS,NPROCS
C
C*** End of declarations rewritten by SPAG
C
      DATA MPI_ID/0/
      DATA MPI_OWN_ISUM/999999/,MPI_OWN_RSUM/999999/
      DATA MPI_OWN_CSUM/999999/
      DATA MPI_KLOOP/.FALSE./,MPI_ELOOP/.FALSE./
C
      ALLOCATABLE MPI_STATUS,IPROC_K,IPROC_E
C
      SAVE IPROC_K,IPROC_E

      END
