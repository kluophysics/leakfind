C*==fpnrcore_tau.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPNRCORE_TAU(IPRINT,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &                        NMAX_LEBGRID)
C   ********************************************************************
C   *                                                                  *
C   *  DUMMY ---> 7.4.1                                                *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TYPES,ONLY:NLMFPMAX
      IMPLICIT NONE
C*--FPNRCORE_TAU11
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NMAX_LEBGRID,N_LEBGRID
      REAL*8 W_LEBGRID(NMAX_LEBGRID,NLMFPMAX),
     &       Y_LEBGRID(NMAX_LEBGRID,NLMFPMAX)
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,*) 'FPNRCORE_TAU ---> 7.4.1'
      IF ( IPRINT.LT.999 ) WRITE (6,*) IPRINT,Y_LEBGRID,W_LEBGRID,
     &                                 N_LEBGRID,NMAX_LEBGRID
      STOP
      END
