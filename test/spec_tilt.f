C*==spec_tilt.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_TILT(THETAINP,PHIINP,TILT,JT,NT,POL0,
     &                     POL0_VEC_TILT,POL0_INITIAL,IP)
C   ********************************************************************
C   *                                                                  *
C   *  Rotates frame of reference for THETAINP and PHIINP arround      *
C   *  axis ROTAXIS with angle BETA                                    *
C   *  e.g.: PES Calculations are always done in the local frame       *
C   *  of reference  (with respect to the surface normal) and input    *
C   *  parameters for theta and phi are to be done in the global       *
C   *  laboratory frame of reference (with respect to detector)        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FOUND_SECTION,FOUND_REAL,N_FOUND
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IP,JT,NT
      REAL*8 PHIINP,THETAINP
      LOGICAL POL0_VEC_TILT,TILT
      REAL*8 POL0(3),POL0_INITIAL(3)
C
C Local variables
C
      REAL*8 BETA,BETA1,BETA2,MXY,NEWVEC(3),NEWVEC_POL(3),PHI,PI,QXY,
     &       ROTAXIS(3),ROTMAT(3,3),ROTMATX(3,3),TET,THETA,VEC(3),X,Y,Z
      REAL*8 DDOT
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      BETA = 99999.0D0
      ROTAXIS(1:3) = 0.0D0
      BETA1 = -30.0D0
      BETA2 = 30.0D0
      PI = 3.141592653589793238462643D0
      CALL INPUT_FIND_SECTION('SPEC_EL',0)
      IF ( FOUND_SECTION ) THEN
         IF ( NT.GT.1 ) THEN
            CALL SECTION_SET_REAL('BETA1',BETA1,9999D0,0)
            CALL SECTION_SET_REAL('BETA2',BETA2,9999D0,0)
         ELSE
            CALL SECTION_SET_REAL('BETA',BETA,9999D0,0)
            BETA1 = BETA
            BETA2 = BETA
         END IF
         IF ( FOUND_REAL ) THEN
            TILT = .TRUE.
            CALL SECTION_SET_REAL_ARRAY('ROTAXIS',ROTAXIS,N_FOUND,3,0,
     &                                  9999D0,1)
         END IF
      END IF
C
C     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      IF ( TILT .AND. JT.EQ.0 ) THEN
         WRITE (6,*) '===================================='
         WRITE (6,*) 'SAMPLE TILT WILL BE CONSIDERED:'
         WRITE (6,*) 'NP (Number of phi angles from input):',NT
         WRITE (6,*) 'BETA1=',BETA1
         WRITE (6,*) 'BETA2=',BETA2
         WRITE (6,*) 'ROTIONAL AXIS:',ROTAXIS
         WRITE (6,*) '===================================='
         RETURN
      ELSE IF ( .NOT.TILT ) THEN
         RETURN
      END IF
C
      IF ( NT.GT.1 ) BETA = BETA1 + (BETA2-BETA1)*DBLE(JT-1)/DBLE(NT-1)
C
      IF ( IP.GT.0 ) THEN
         WRITE (6,*) '======================================='
         WRITE (6,*) 'SAMPLE TILT JP=:',JT
         WRITE (6,*) 'BETA=:',BETA
         WRITE (6,*) 'THETAINP=:',THETAINP*180D0/PI
         WRITE (6,*) 'PHIINP=:',PHIINP*180D0/PI
         WRITE (6,*) '======================================='
      END IF
C
      THETA = BETA*PI/180D0
C
      X = ROTAXIS(1)
      Y = ROTAXIS(2)
      Z = ROTAXIS(3)
C
      TET = THETAINP
      PHI = PHIINP
C
      VEC(3) = COS(TET)
      QXY = SIN(TET)
      VEC(1) = COS(PHI)*QXY
      VEC(2) = SIN(PHI)*QXY
C
      ROTMAT(1,1) = COS(THETA) + (1.0D0-COS(THETA))*X**2
      ROTMAT(1,2) = (1.0D0-COS(THETA))*X*Y - SIN(THETA)*Z
      ROTMAT(1,3) = (1.0D0-COS(THETA))*X*Z - SIN(THETA)*Y
C
      ROTMAT(2,1) = (1.0D0-COS(THETA))*Y*X + SIN(THETA)*Z
      ROTMAT(2,2) = COS(THETA) + (1.0D0-COS(THETA))*Y**2
      ROTMAT(2,3) = (1.0D0-COS(THETA))*Y*Z - SIN(THETA)*X
C
      ROTMAT(3,1) = (1.0D0-COS(THETA))*Z*X - SIN(THETA)*Y
      ROTMAT(3,2) = (1.0D0-COS(THETA))*Z*Y + SIN(THETA)*X
      ROTMAT(3,3) = COS(THETA) + (1.0D0-COS(THETA))*Z**2
C
      ROTMATX = ROTMAT
      DO I = 1,3
         DO J = 1,3
            ROTMAT(I,J) = ROTMATX(J,I)
         END DO
      END DO
C
      IF ( .NOT.POL0_VEC_TILT ) THEN
         NEWVEC = MATMUL(ROTMAT,VEC)
         MXY = SQRT(NEWVEC(1)**2+NEWVEC(2)**2)
      ELSE
         IF ( IP.GT.0 ) THEN
            WRITE (*,*) ''
            WRITE (*,*) 'OLD_VEC'
            WRITE (*,*) POL0_INITIAL
            WRITE (*,*) ''
         END IF
         NEWVEC = MATMUL(ROTMAT,VEC)
         MXY = SQRT(NEWVEC(1)**2+NEWVEC(2)**2)
C
         NEWVEC_POL = MATMUL(ROTMAT,POL0_INITIAL)
         POL0 = NEWVEC_POL
C
         IF ( IP.GT.0 ) THEN
            WRITE (*,*) THETA
            WRITE (*,*) 'NEW_VEC'
            WRITE (*,*) POL0
            WRITE (*,*) ''
         END IF
      END IF
C
      IF ( ABS(MXY).LT.1D-8 ) THEN
         PHI = 0D0
      ELSE IF ( NEWVEC(2).GE.0D0 ) THEN
         PHI = ACOS(NEWVEC(1)/MXY)
      ELSE IF ( NEWVEC(1).LT.0D0 ) THEN
         PHI = PI + ACOS(-NEWVEC(1)/MXY)
      ELSE
         PHI = 2*PI - ACOS(NEWVEC(1)/MXY)
      END IF
      TET = ACOS(NEWVEC(3)/SQRT(DDOT(3,NEWVEC(1),1,NEWVEC(1),1)))
      THETAINP = TET
      PHIINP = PHI
      END
