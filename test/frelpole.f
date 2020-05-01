C*==frelpole.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FRELPOLE(ERYDTOP,IPRINT)
C   ********************************************************************
C   *                                                                  *
C   *   generate set of recip. vectors  GN  that might cause a         *
C   *   free electron pole   (k+GN)**2 - E = 0                         *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_LATTICE,ONLY:ALAT,BBAS
      USE MOD_KSPACE,ONLY:GFEP,NGFEPMAX,NKTAB,KTAB,NGFEP,IBZINT,KTET,
     &    NKTET
      USE MOD_STR,ONLY:G3,G2,G1,NGRL
      IMPLICIT NONE
C*--FRELPOLE15
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ERYDTOP
      INTEGER IPRINT
C
C Local variables
C
      REAL*8 EDUTOP,GX,GY,GZ,KNSQ
      INTEGER I,IK,IX
C
C*** End of declarations rewritten by SPAG
C
      EDUTOP = ERYDTOP/(2*PI/ALAT)**2
      NGFEP = 0
C
C-------------------------------  scan all reciprocal lattice  k-vectors
      DO I = 1,NGRL
C
         GX = G1(I)*BBAS(1,1) + G2(I)*BBAS(1,2) + G3(I)*BBAS(1,3)
         GY = G1(I)*BBAS(2,1) + G2(I)*BBAS(2,2) + G3(I)*BBAS(2,3)
         GZ = G1(I)*BBAS(3,1) + G2(I)*BBAS(3,2) + G3(I)*BBAS(3,3)
C
C-----------------------------------  scan all tabulated k-vectors in BZ
C
         IF ( IBZINT.NE.4 ) THEN
C
            DO IK = 1,NKTAB
               KNSQ = (KTAB(1,IK)+GX)**2 + (KTAB(2,IK)+GY)
     &                **2 + (KTAB(3,IK)+GZ)**2
               IF ( KNSQ.LT.EDUTOP ) THEN
                  NGFEP = NGFEP + 1
                  IF ( NGFEP.LE.NGFEPMAX ) THEN
                     GFEP(1,NGFEP) = GX
                     GFEP(2,NGFEP) = GY
                     GFEP(3,NGFEP) = GZ
                  END IF
                  EXIT
               END IF
            END DO
C
         ELSE
C
            DO IK = 1,NKTET
               KNSQ = (KTET(1,IK)+GX)**2 + (KTET(2,IK)+GY)
     &                **2 + (KTET(3,IK)+GZ)**2
               IF ( KNSQ.LT.EDUTOP ) THEN
                  NGFEP = NGFEP + 1
                  IF ( NGFEP.LE.NGFEPMAX ) THEN
                     GFEP(1,NGFEP) = GX
                     GFEP(2,NGFEP) = GY
                     GFEP(3,NGFEP) = GZ
                  END IF
                  EXIT
               END IF
            END DO
C
         END IF
C
      END DO
C-----------------------------------------------------------------------
C
      IF ( NGFEP.GT.NGFEPMAX ) THEN
         WRITE (6,*) ' STOP in <FRELPOLE>  NGFEP    =',NGFEP
         WRITE (6,*) ' array size:         NGFEPMAX =',NGFEPMAX
         STOP
      END IF
C
      IF ( IPRINT.LT.1 ) RETURN
C-----------------------------------------------------------------------
      WRITE (6,99001)
      DO I = 1,NGFEP
         WRITE (6,99002) I,(GFEP(IX,I),IX=1,3)
      END DO
      WRITE (6,*) ' '
C-----------------------------------------------------------------------
99001 FORMAT (/,1X,79('*'),/,35X,'<FRELPOLE>',/,1X,79('*'),//,10X,
     &    'reciprocal lattice vectors with possible free electron poles'
     &    ,/)
99002 FORMAT (10X,I3,3F10.4)
      END
