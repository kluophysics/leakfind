C*==init_spiral.f    processed by SPAG 6.70Rc at 16:59 on 31 Jan 2017
      SUBROUTINE INIT_SPIRAL(KMROT,QMVEC,MDIRQ,QMTET,QMPHI,NQ,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SUBROUTINE TO INITIALIZE KMROT                                 *
C   *   AND PREPARE VALUES FOR SPINSPIRALS                             *
C   *                                                                  *
C   *   HERE KMROT =                                                   *
C   *        3: SPIN SPIRAL         THETA=90                           *
C   *        4: SPIN SPIRAL         THETA<>90                          *
C   *        5: SPIN SPIRAL with    THETA/PHI from MQDIR               *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:DEG_ARC,ARC_DEG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KMROT,NQ,NQMAX
      REAL*8 MDIRQ(3,NQMAX),QMPHI(NQMAX),QMTET(NQMAX),QMVEC(3)
C
C Local variables
C
      REAL*8 DDOT
      INTEGER I,IQ
      REAL*8 P,PHIM,PHIV,QMVECN(3),STET,TETV,VX,VXY,VXYZ,VY,VZ
      LOGICAL PTEST,TTEST
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
CCGHF NOTES:
C     ROTATION PHI OF MAGNETIZATION ABOUT QVEC IS NOT YET PERFECT !
C                  AND MAY NEED IMPROVEMENT.
C     DIFFERENT SIGN OF THETA AT DIFFERENT SITES RESULTS IN AFM STATE !
C
C     Q IS GIVEN IN CARTHESIANS THAT IS IN HEXAGONAL OR TETRAGONAL
C       SYSTEMS [0 0 1] SHOULD BE GIVEN IN THE INPUT AS [0 0 a/c]
C-----------------------------------------------------------------------
C
      VXYZ = SQRT(QMVEC(1)**2+QMVEC(2)**2+QMVEC(3)**2)
C
      IF ( VXYZ.GT.1.D-6 ) THEN
C
C------------------------------------------------- MAKE UNIT VECTOR OF Q
         VX = QMVEC(1)/VXYZ
         VY = QMVEC(2)/VXYZ
         VZ = QMVEC(3)/VXYZ
         VXY = SQRT(VX**2+VY**2)
C
C--------------------------------------------------- NORMALISED Q-VECTOR
         QMVECN(1) = VX
         QMVECN(2) = VY
         QMVECN(3) = VZ
C
         PTEST = .FALSE.
C
         IF ( KMROT.NE.5 ) THEN
C
C------------------------ SET M PERPENDICULAR TO Q IF QMTET = 90 degrees
C
            DO IQ = 1,NQ
               TETV = ACOS(VZ)*ARC_DEG + QMTET(IQ)
               PHIV = ATAN2(VY,VX)*ARC_DEG + QMPHI(IQ)
               STET = SIN(TETV*DEG_ARC)
               MDIRQ(1,IQ) = STET*COS(PHIV*DEG_ARC)
               MDIRQ(2,IQ) = STET*SIN(PHIV*DEG_ARC)
               MDIRQ(3,IQ) = COS(TETV*DEG_ARC)
            END DO
            PTEST = .TRUE.
C
         ELSE
C
C--------------------------------------- CALCULATE ANGLE BETWEEN M AND Q
C
            DO IQ = 1,NQ
               P = DDOT(3,MDIRQ(1,IQ),1,QMVECN,1)
               IF ( ABS(P).LE.1D-6 ) THEN
                  PTEST = .TRUE.
                  QMTET(IQ) = ACOS(P)*ARC_DEG
C                   QMPHI IS TAKEN AS DIFFERENCE OF THE PROJECTIONS OF
C                         M AND Q ONTO THE X-Y PLANE
                  IF ( VXY.GT.1.D-6 ) THEN
                     PHIV = ATAN2(VY,VX)
                     PHIM = ATAN2(MDIRQ(2,IQ),MDIRQ(1,IQ))
                     QMPHI(IQ) = (PHIM-PHIV)*ARC_DEG
                  ELSE
                     QMPHI(IQ) = ATAN2(MDIRQ(2,IQ),MDIRQ(1,IQ))
                  END IF
               ELSE
                  QMTET(IQ) = 0.D0
                  QMPHI(IQ) = 0.D0
               END IF
            END DO
C
         END IF
C
C------------ TEST THAT AT LEAST ONE COMPONENT OF M IS NOT PARALLEL TO Q
C
         IF ( PTEST ) THEN
            KMROT = 4
         ELSE
            STOP 'NO M COMPONENT PERPENDICULAR TO Q FOUND'
         END IF
C
C---------------------------- TEST IF QMTET = +-90 degrees  AT ALL SITES
C
         TTEST = .TRUE.
         DO IQ = 1,NQ
            IF ( ABS(COS(QMTET(IQ))).GT.1.D-6 ) TTEST = .FALSE.
         END DO
         IF ( TTEST ) KMROT = 3
C
         WRITE (77,99001) QMVEC
         DO IQ = 1,NQ
            WRITE (77,99002) IQ,(MDIRQ(I,IQ),I=1,3)
         END DO
         WRITE (77,99003) (QMTET(IQ),IQ=1,NQ)
         WRITE (77,99004) (QMPHI(IQ),IQ=1,NQ)
      END IF
C
      RETURN
99001 FORMAT ('*  QVEC           =',2X,3F10.5)
99002 FORMAT ('*  MVEC(',I3,')      =',2X,3F10.5)
99003 FORMAT ('*  QMTET(IQ)      =',2X,9F8.2)
99004 FORMAT ('*  QMPHI(IQ)      =',2X,9F8.2)
      END
