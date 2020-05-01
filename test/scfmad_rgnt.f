C*==scfmad_rgnt.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFMAD_RGNT(NLVMAD_LOC,NLQMAD_LOC,RGNTMAD,IRGNTMAD,
     &                       NRGNTMAD,NRGNTMADMAX)
C   ********************************************************************
C   *                                                                  *
C   *  set up of the Gaunt coefficients with an index field            *
C   *  recognize that they are needed here only for L3=LV+LQ           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SCFMAD_RGNT12
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLQMAD_LOC,NLVMAD_LOC,NRGNTMAD,NRGNTMADMAX
      INTEGER IRGNTMAD(NRGNTMADMAX,3)
      REAL*8 RGNTMAD(NRGNTMADMAX)
C
C Local variables
C
      REAL*8 GAUNT_RYLM
      INTEGER I,L3,LQ,LV,M3,MQ,MV
      REAL*8 RGAUNT
C
C*** End of declarations rewritten by SPAG
C
      I = 1
      DO LV = 0,(NLVMAD_LOC-1)
         DO LQ = 0,(NLQMAD_LOC-1)
            L3 = LV + LQ
            DO MV = -LV,LV
               DO MQ = -LQ,LQ
                  DO M3 = -L3,L3
C
                     RGAUNT = GAUNT_RYLM(LV,MV,LQ,MQ,L3,M3)
C
                     IF ( ABS(RGAUNT).GT.1.D-10 ) THEN
                        IF ( I.GT.NRGNTMADMAX ) THEN
                           WRITE (6,99001) I,NRGNTMADMAX
                           STOP 'in <SCFMAD_RGNT> '
                        END IF
C
                        RGNTMAD(I) = RGAUNT
                        IRGNTMAD(I,1) = LV*(LV+1) + MV + 1
                        IRGNTMAD(I,2) = LQ*(LQ+1) + MQ + 1
                        IRGNTMAD(I,3) = L3*(L3+1) + M3 + 1
                        I = I + 1
                     END IF
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      NRGNTMAD = I - 1
C
      RETURN
99001 FORMAT ('<SCFMAD_RGNT>: I=',I5,' > NRGNTMADMAX =',I5)
      END
