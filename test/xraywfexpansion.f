C*==xrayexpansion.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XRAYEXPANSION(MEZJ,MEZZ,MSST,SSST,TAUT,TSST,IXRAYCORWF)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      IMPLICIT NONE
C*--XRAYEXPANSION10
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='XRAYEXPANSION')
C
C Dummy arguments
C
      INTEGER IXRAYCORWF
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
      CALL STOP_MESSAGE(ROUTINE,'routine not available at the moment')
C=======================================================================
C
      WRITE (*,*) MEZJ,MEZZ,MSST,SSST,TAUT,TSST,IXRAYCORWF
C
      END
