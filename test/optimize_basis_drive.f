C*==optimize_basis_drive.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE OPTIMIZE_BASIS_DRIVE(ABAS,QBAS,NQ,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *           optimize basis geometry if requested                   *
C   *                                                                  *
C   * this is done only in case of an SCF start for a host system      *
C   * or the key OPTBASIS is set in section CONTROL of the input file  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      USE MOD_FILES,ONLY:IOTMP,FOUND_SECTION
      USE MOD_SITES,ONLY:OPTIMIZE_BASIS_GEOMETRY
      USE MOD_SCF,ONLY:SCFSTATUS_HOST
      USE MOD_CALCMODE,ONLY:PROGNAME
      IMPLICIT NONE
C*--OPTIMIZE_BASIS_DRIVE18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='OPTIMIZE_BASIS_DRIVE')
C
C Dummy arguments
C
      INTEGER NQ,NQMAX
      REAL*8 ABAS(3,3),QBAS(3,NQMAX)
C
C Local variables
C
      INTEGER IFLAG
C
C*** End of declarations rewritten by SPAG
C
      OPTIMIZE_BASIS_GEOMETRY = .FALSE.
      IFLAG = 0
C
      IF ( SYSTEM_DIMENSION(1:2).NE.'2D' .AND. SYSTEM_DIMENSION(1:2)
     &     .NE.'3D' ) RETURN
C
      IF ( (SCFSTATUS_HOST(1:5).EQ.'START' .AND. PROGNAME(4:6).EQ.'SCF')
     &     ) THEN
C
         OPTIMIZE_BASIS_GEOMETRY = .TRUE.
C
      ELSE
C
         CALL INPUT_FIND_SECTION('CONTROL',0)
C
         IF ( FOUND_SECTION ) CALL SECTION_FIND_KEYWORD('OPTBASIS',
     &        OPTIMIZE_BASIS_GEOMETRY)
C
      END IF
C
      IF ( .NOT.OPTIMIZE_BASIS_GEOMETRY ) RETURN
C
C      IF ( NQ.GT.0 ) RETURN
C
C ======================================================================
C                      optimize basis geometry
C ======================================================================
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D        ' ) THEN
C
         CALL OPTIMIZE_BASIS_2D(ABAS,QBAS,NQ,NQMAX,IFLAG)
C
      ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'3D        ' ) THEN
C
         CALL OPTIMIZE_BASIS_3D(ABAS,QBAS,NQ,NQMAX,IFLAG,IOTMP)
C
      END IF
C
C      IF ( IFLAG.LT.0 )
C     &      CALL STOP_MESSAGE(ROUTINE,'optimisation of basis failed')
C
C ======================================================================
      END
