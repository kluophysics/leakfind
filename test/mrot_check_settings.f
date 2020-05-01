C*==mrot_check_settings.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MROT_CHECK_SETTINGS(NQ,IREL,QMPHI,QMTET,QMGAM,KMROT,
     &                               QMVEC,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *  check and correct settings for calculations                     *
C   *  with rotated magnetic moments                                   *
C   *                                                                  *
C   *  QMVEC <> 0-vector  >>>  spin spiral calculation                 *
C   *                                                                  *
C   *       ==>  check Theta to fix KMROT = 3 or 4                     *
C   *                                                                  *
C   *  ELSE                                                            *
C   *                                                                  *
C   *       ==>  non-0 angles TET and PHI                              *
C   *                                                                  *
C   *            all angles the same:    KMROT = 1                     *
C   *            all angles individual:  KMROT = 2                     *
C   *                                                                  *
C   *      ELSE                                                        *
C   *                                                                  *
C   *            no rotation: KMROT = 0                                *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  KMROT                                                           *
C   *                                                                  *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site IQ   *
C   *  2: global COMMON rotation of the magnetisation                  *
C   *  3: spin spiral    Theta =  90 }    QMVEC <> 0-vector            *
C   *  4: spin spiral    Theta <> 90 }    IREL = 2                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:MOMENTS_ROTATED
      IMPLICIT NONE
C*--MROT_CHECK_SETTINGS37
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MROT_CHECK_SETTINGS')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-6)
C
C Dummy arguments
C
      INTEGER IREL,KMROT,NQ,NQMAX
      REAL*8 QMGAM(NQMAX),QMPHI(NQMAX),QMTET(NQMAX),QMVEC(3)
C
C Local variables
C
      INTEGER I,IQ,KMROT_ORIGINAL
      REAL*8 RARG
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      KMROT_ORIGINAL = KMROT
C
      KMROT = 0
C
      DO I = 1,3
         IF ( RNON0(QMVEC(I)) ) KMROT = 4
      END DO
C
C-----------------------------------------------------------------------
C                       spin spiral calculation
C-----------------------------------------------------------------------
C
      IF ( KMROT.EQ.4 ) THEN
C
         IF ( IREL.NE.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'KMROT = 3 or 4 but IREL <> 2')
C
         DO IQ = 1,NQ
            IF ( RNON0(90D0-QMTET(IQ)) ) THEN
               KMROT = 3
               EXIT
            END IF
         END DO
C
C-----------------------------------------------------------------------
C           collinear  or non-collinear spin structure
C-----------------------------------------------------------------------
      ELSE
C
         KMROT = 0
         DO IQ = 1,NQ
            IF ( RNON0(QMGAM(IQ)) ) THEN
               KMROT = 1
               EXIT
            END IF
            IF ( RNON0(QMTET(IQ)) ) THEN
               KMROT = 2
               EXIT
            END IF
            IF ( RNON0(QMPHI(IQ)) ) THEN
               KMROT = 2
               EXIT
            END IF
         END DO
C
         IF ( KMROT.EQ.2 ) THEN
C
            DO IQ = 2,NQ
               IF ( RNON0(QMTET(IQ)-QMTET(1)) ) THEN
                  KMROT = 1
                  EXIT
               END IF
               IF ( RNON0(QMPHI(IQ)-QMPHI(1)) ) THEN
                  KMROT = 1
                  EXIT
               END IF
            END DO
C
         END IF
C
      END IF
C
C-----------------------------------------------------------------------
C
      IF ( KMROT.NE.KMROT_ORIGINAL ) WRITE (6,99001)
     &     ROUTINE(1:LEN_TRIM(ROUTINE)),KMROT_ORIGINAL,KMROT
C
      IF ( KMROT.EQ.1 .OR. KMROT.EQ.2 ) MOMENTS_ROTATED = .TRUE.
C
C=======================================================================
99001 FORMAT (//,70('*'),/,10X,'WARNING from  <<',A,'>> ',/,10X,
     &        'the flag KMROT has been adjusted from ',I3,'  to ',i3,/,
     &        /,70('*'),/)
C
      END
