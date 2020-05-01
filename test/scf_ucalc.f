C*==scf_ucalc.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCF_UCALC(IT_CENTER,IQCLU_CENTER)
C   ********************************************************************
C   *                                                                  *
C   *  driver routine for the calculation of the U parameter           *
C   *  using a real space cluster approach                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_FILES,ONLY:IFILUCALC,FILNAM,LFILNAM,LSYSTEM,SYSTEM
      USE MOD_SITES,ONLY:NQHOST,NQCLU,QBAS_QCLU,IT_OQCLU,NO_QCLU
      USE MOD_TYPES,ONLY:LOPT,ITBOT,ITTOP
      IMPLICIT NONE
C*--SCF_UCALC15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF_UCALC')
C
C Dummy arguments
C
      INTEGER IQCLU_CENTER,IT_CENTER
C
C Local variables
C
      INTEGER I,IQCLU,IT,J
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      FILNAM = 'UCALC.dat'
      LFILNAM = LEN_TRIM(FILNAM)
C
      OPEN (UNIT=IFILUCALC,FILE=FILNAM(1:LFILNAM))
C
      WRITE (IFILUCALC,99003) SYSTEM(1:LSYSTEM),NQCLU,IQCLU_CENTER,
     &                        IQCLU_CENTER + NQHOST,IT_CENTER,
     &                        LOPT(IT_CENTER),ITBOT,ITTOP
C
      WRITE (IFILUCALC,99004)
      DO IQCLU = 1,NQCLU
         IF ( NO_QCLU(IQCLU).GT.1 ) STOP '<SCF_UCALC> NO_QCLU > 1'
         IT = IT_OQCLU(1,IQCLU)
         WRITE (IFILUCALC,99005) IQCLU,IT,QBAS_QCLU(1:3,IQCLU)
      END DO
C
      WRITE (IFILUCALC,99002)
      WRITE (IFILUCALC,99001) ((BBAS(I,J),I=1,3),J=1,3)
C
99001 FORMAT (3D20.12)
99002 FORMAT (/,10X,'primitive vectors of reciprocal space',
     &        ' Bravais lattice',/,10X,'       in units of 2*pi/a',/)
99003 FORMAT (//,2(/,1X,79('U')),/,10X,
     &        'calculation of U-parameter for ',A,/,10X,
     &        'cluster size NQCLU =',I5,/,10X,'central site  IQCLU =',
     &        I3,'    IQ =',I3,'    type IT =',I3,'    l =',I2,/,10X,
     &        'ITBOT =',I3,3X,'ITTOP =',I3,/)
99004 FORMAT ('     IQCLU   IT       QBAS(X)',13X,'QBAS(Y)',13X,
     &        'QBAS(Z)')
99005 FORMAT (I10,I5,3D20.12)
C
      END
C*==scf_ucalc_dndv.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCF_UCALC_DNDV(KEY,OBS_LT)
C   ********************************************************************
C   *                                                                  *
C   *  auxilary routine for the calculation of the U parameter         *
C   *                                                                  *
C   *  calculation and writing of the difference quotient              *
C   *                                                                  *
C   *                    d N_l / d V_shift                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS
      USE MOD_SCF,ONLY:ITRSCF
      USE MOD_CALCMODE,ONLY:U_CALCULATION,U_POT_SHIFT,
     &    IT_CENTER_U_CALCULATION
      USE MOD_ANGMOM,ONLY:NLMAX,NL
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX,LOPT
      USE MOD_FILES,ONLY:IFILUCALC,IOTMP,IDUMMY
      IMPLICIT NONE
C*--SCF_UCALC_DNDV102
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF_UCALC_DNDV')
C
C Dummy arguments
C
      INTEGER KEY
      REAL*8 OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX)
C
C Local variables
C
      REAL*8 DNDV,DNDV1,NOS0_LT(:,:),NOS1_LT(:,:)
      INTEGER IL,IT
      LOGICAL INITIALIZE
      SAVE NOS0_LT,NOS1_LT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NOS0_LT,NOS1_LT
C
      DATA INITIALIZE/.TRUE./
C
C=======================================================================
C
      IF ( .NOT.U_CALCULATION ) RETURN
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (NOS0_LT(NLMAX,NTMAX),NOS1_LT(NLMAX,NTMAX))
C
         INITIALIZE = .FALSE.
      END IF
C=======================================================================
C
C
      IF ( ITRSCF.EQ.1 ) THEN
C
         NOS1_LT(1:NL,ITBOT:ITTOP) = OBS_LT(0,IDOS,1:NL,ITBOT:ITTOP)
C
         RETURN
C
      END IF
C
C=======================================================================
C
      IF ( KEY.EQ.0 ) RETURN
C
C
      IL = LOPT(IT_CENTER_U_CALCULATION) + 1
C
      WRITE (6,99001) U_POT_SHIFT
C
      IF ( ABS(U_POT_SHIFT).LT.1D-5 ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'NOS0.dat')
         DO IT = ITBOT,ITTOP
            WRITE (IOTMP,99002) IT,OBS_LT(0,IDOS,IL,IT)
         END DO
         CLOSE (IOTMP)
         STOP
      END IF
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'NOS0.dat')
C
      DO IT = ITBOT,ITTOP
         READ (IOTMP,*) IDUMMY,NOS0_LT(IL,IT)
C
         DNDV1 = (NOS1_LT(IL,IT)-NOS0_LT(IL,IT))/U_POT_SHIFT
         DNDV = (OBS_LT(0,IDOS,IL,IT)-NOS0_LT(IL,IT))/U_POT_SHIFT
C
C
         WRITE (6,99002) IT,NOS0_LT(IL,IT),NOS1_LT(IL,IT),
     &                   OBS_LT(0,IDOS,IL,IT),DNDV1,DNDV
C
         WRITE (IFILUCALC,99002) IT,NOS0_LT(IL,IT),NOS1_LT(IL,IT),
     &                           OBS_LT(0,IDOS,IL,IT),DNDV1,DNDV
      END DO
C
      RETURN
C=======================================================================
C
99001 FORMAT (//,2(/,1X,79('U')),/,10X,
     &        'calculation of U-parameter for ',/,10X,'U_POT_SHIFT = ',
     &        F12.5,/,2(/,1X,79('U')),//)
99002 FORMAT (i10,3F20.10,2X,2F20.10)
      END
