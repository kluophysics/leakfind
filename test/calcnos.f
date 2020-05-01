C*==calcnos.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE CALCNOS(CTOTDOS,PHASK,PHAST,PHASA)
C   ********************************************************************
C   *                                                                  *
C   *  calculate NOS using the Lloyd formula and write result          *
C   *  the necessary data are created by  <LLOYDPT1>  and  <LLOYDPT2>  *
C   *                                                                  *
C   *  the corresp. DOS is calculated via numerical differantiation    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILNOS
      USE MOD_SITES,ONLY:IQAT
      USE MOD_TYPES,ONLY:NAT,CONC,NTMAX,NT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ
      USE MOD_ENERGY,ONLY:NETAB,ETAB,NEMAX
      USE MOD_CONSTANTS,ONLY:C0,PI
      IMPLICIT NONE
C*--CALCNOS19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CTOTDOS(NEMAX),PHASA(NKMMAX,NTMAX,NEMAX),PHASK(NEMAX),
     &           PHAST(NKMMAX,NTMAX,NEMAX)
C
C Local variables
C
      COMPLEX*16 CSSNOS(:),CSUMA,CSUMSS,CTOTNOS(:),WT
      INTEGER I,IA_ERR,IE,IT,N,NELD,NELDMAX
      REAL*8 RETAB(:),RTOTDOS(:),RTOTNOS(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CTOTNOS,CSSNOS,RTOTNOS,RTOTDOS,RETAB
C
      NELD = NETAB(1)
C
      NELDMAX = NELD
C
      ALLOCATE (RTOTNOS(NELDMAX),RTOTDOS(NELDMAX),RETAB(NELDMAX))
      ALLOCATE (CSSNOS(NELDMAX),CTOTNOS(NELDMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: PCOEF -> LLOYDPT2'
C
C------------------------------------- remove jumps from PHAST and PHASA
C                                   this is already done in   <LLOYDPT2>
C
C     CALL SBRANCH(PHAST,NKM,NT,NKMMAX,NTMAX,NELD)
C     CALL SBRANCH(PHASA,NKM,NT,NKMMAX,NTMAX,NELD)
C
C------------------------------------------------- sum up to get CTOTNOS
C
      DO IE = 1,NELD
         CSUMA = C0
         CSUMSS = C0
         DO IT = 1,NT
            N = NKMQ(IQAT(1,IT))
            WT = NAT(IT)*CONC(IT)
            DO I = 1,N
               CSUMSS = CSUMSS - WT*(PHAST(I,IT,IE)-PHASA(I,IT,IE))
               CSUMA = CSUMA + WT*PHASA(I,IT,IE)
            END DO
         END DO
C
         CTOTNOS(IE) = (PHASK(IE)+CSUMSS)/PI
         CSSNOS(IE) = CSUMA/PI
C
         RTOTNOS(IE) = DIMAG(CTOTNOS(IE))
         RETAB(IE) = DREAL(ETAB(IE,1))
C
C
      END DO
C
C ======================================================================
C
      CALL RDFDX(RTOTNOS,RETAB,RTOTDOS,NELD)
C
C ======================================================================
C
      DO IE = 1,NELD
         WRITE (IFILNOS,99001) RETAB(IE),RTOTNOS(IE),DIMAG(CSSNOS(IE)),
     &                         RTOTDOS(IE),DIMAG(CTOTDOS(IE))
      END DO
C
C ======================================================================
C
      DEALLOCATE (CTOTNOS,CSSNOS)
C
C ----------------------------------------------------------------------
99001 FORMAT (12E20.6)
      END
