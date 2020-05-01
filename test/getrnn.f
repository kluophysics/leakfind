C*==getrnn.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GETRNN
C   ********************************************************************
C   *                                                                  *
C   *   find the number and distance of nearest neighbor atoms         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ,QBAS,RNNQ,NNNQ
      USE MOD_LATTICE,ONLY:ABAS,ALAT
      IMPLICIT NONE
C*--GETRNN12
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 DNRM2
      REAL*8 DQVEC(3),DR,DRNN,DRVEC(3)
      INTEGER I1,I2,I3,IQ,JQ,LOOP
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,99001)
C
      DO IQ = 1,NQ
         DRNN = 1D9
         NNNQ(IQ) = 0
         DO LOOP = 1,2
            DO JQ = 1,NQ
               DQVEC(1:3) = QBAS(1:3,JQ) - QBAS(1:3,IQ)
               DO I1 = -1, + 1
                  DO I2 = -1, + 1
                     DO I3 = -1, + 1
C
C
                        CALL RVECLCIB(I1,I2,I3,ABAS,DRVEC)
C
                        DRVEC(1:3) = DRVEC(1:3) + DQVEC(1:3)
C
                        DR = DNRM2(3,DRVEC,1)
C
                        IF ( DR.GT.1.0D-8 ) DRNN = MIN(DRNN,DR)
                        IF ( LOOP.EQ.2 ) THEN
                           IF ( ABS(DR-DRNN).LT.1D-6 ) NNNQ(IQ)
     &                          = NNNQ(IQ) + 1
                        END IF
                     END DO
                  END DO
               END DO
            END DO
            IF ( LOOP.EQ.2 ) THEN
               RNNQ(IQ) = DRNN*ALAT
               WRITE (6,99002) IQ,NNNQ(IQ),RNNQ(IQ),RNNQ(IQ)/ALAT
            END IF
         END DO
C
      END DO
C
99001 FORMAT (/,10X,'<GETRNN> number and distance of',
     &        ' nearest neighbor atoms')
99002 FORMAT (10X,'IQ=',I2,I4,' NN''s with ','r_NN=',F12.5,' a.u. ',
     &        F12.5,' a ')
      END
