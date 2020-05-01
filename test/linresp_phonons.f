C*==linresp_phonons.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_PHONONS
C   ********************************************************************
C   *                                                                  *
C   *  main subroutine for phonon spectra                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:PHASK
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_FILES,ONLY:WRTAU,RDTAU
      USE MOD_ANGMOM,ONLY:MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX
      USE MOD_MPI,ONLY:MPI
      IMPLICIT NONE
C*--LINRESP_PHONONS16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_PHONONS')
C
C Local variables
C
      REAL*8 RHO2NSX(:,:,:,:),TIME1,TIME2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RHO2NSX
C
      ALLOCATE (RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3))
C
      CALL UNDER_CONSTRUCTION(ROUTINE)
C
C ======================================================================
C
      CALL LINRESP_INIT('PHONONS   ')
C
C ======================================================================
C
      WRITE (6,99001)
C
C=======================================================================
C
      CALL CPU_TIME(TIME1)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
      WRITE (6,99002) RDTAU,WRTAU
C
C ======================================================================
C
      CALL LINRESP_ELOOP(MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST,
     &                   RHO2NSX,PHASK)
C
      CALL CPU_TIME(TIME2)
C
      WRITE (6,99003) TIME2 - TIME1
C
      CALL UNDER_CONSTRUCTION(ROUTINE)
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
99001 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*  *****   *    *   ****   *     *   ****   *     *   ****   *'
     &  ,/,10X,
     &  '*  *    *  *    *  *    *  **    *  *    *  **    *  *    *  *'
     &  ,/,10X,
     &  '*  *    *  *    *  *    *  * *   *  *    *  * *   *  *       *'
     &  ,/,10X,
     &  '*  *****   ******  *    *  *  *  *  *    *  *  *  *   ****   *'
     &  ,/,10X,
     &  '*  *       *    *  *    *  *   * *  *    *  *   * *       *  *'
     &  ,/,10X,
     &  '*  *       *    *  *    *  *    **  *    *  *    **  *    *  *'
     &  ,/,10X,
     &  '*  *       *    *   ****   *     *   ****   *     *   ****   *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,10X,'RDTAU     =',L2,'  WRTAU     =',L2)
99003 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
C
      END
