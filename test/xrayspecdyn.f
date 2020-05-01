C*==xrayspecdyn.f    processed by SPAG 6.70Rc at 15:53 on 19 Dec 2016
      SUBROUTINE XRAYSPECDYN(CALCINT,GETIRRSOL,ETAB0,TAUQ,TAUT,MSSQ,
     &                       TSST,MSST,SSST,MEZZ,MEZJ,GCOR,FCOR,ECOR,
     &                       SZCOR,KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,
     &                       ITXRAY,BCOR,BCORS,NCSTMAX,NPOLMAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of  x-ray-absorption spectra for               *
C   *       circularly and linearly polarized radiation                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_SITES,ONLY:NQMAX
      USE MOD_TYPES,ONLY:NTMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      INTEGER ITXRAY,NCSTMAX,NPOLMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      COMPLEX*16 ETAB0(NEMAX),MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
C
C*** End of declarations rewritten by SPAG
C
C
C
C=======================================================================
      CALL STOP_MESSAGE('XRAYSPECDYN',
     &                  'routine not available at the moment')
C=======================================================================
C
      WRITE (*,*) CALCINT,GETIRRSOL,ETAB0,TAUQ,TAUT,MSSQ,TSST,MSST,SSST,
     &            MEZZ,MEZJ,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &            IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX,NPOLMAX
      END
