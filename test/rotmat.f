C*==rotmat.f    processed by SPAG 6.70Rc at 08:13 on  8 Mar 2017
      SUBROUTINE ROTMAT(NK,IREL,ALFDEG,BETDEG,GAMDEG,ROT,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
C   *           ( ALFDEG, BETDEG, GAMDEG )                             *
C   *                                                                  *
C   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
C   *            EQS. (4.8), (4.12) AND (4.13)                         *
C   *                                                                  *
C   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
C   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
C   *                                                                  *
C   *   27/12/02  HE  deal with beta = 0, 180                          *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:FACT
      USE MOD_CONSTANTS,ONLY:CI,C0,DEG_ARC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ROTMAT')
C
C Dummy arguments
C
      REAL*8 ALFDEG,BETDEG,GAMDEG
      INTEGER IREL,NK,NKMMAX
      COMPLEX*16 ROT(NKMMAX,NKMMAX)
C
C Local variables
C
      REAL*8 CB05,CB05PW,CB05SQ,DOM,J,M1,M2,MSB05,MSB05PW,MSB05SQ,NUM,
     &       RSUM,X
      COMPLEX*16 EMIM1G,EMIM2A
      INTEGER IM1,IM2,K,L,NMUE,OFF,S,SHIGH,SLOW
      REAL*8 RFAC
C
C*** End of declarations rewritten by SPAG
C
C inline function    factorial for real argument
C
      RFAC(X) = FACT(NINT(X))
C
      IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,
     &     'routine called for IREL = 2   deal with that case outside')
      IF ( IREL.EQ.3 .AND. MOD(NK,2).EQ.0 ) CALL STOP_MESSAGE(ROUTINE,
     &     'routine called for IREL = 3   and   even  NK')
C
      ROT(1:NKMMAX,1:NKMMAX) = C0
C
      CB05 = DCOS(BETDEG*0.5D0*DEG_ARC)
      CB05SQ = CB05*CB05
      MSB05 = -DSIN(BETDEG*0.5D0*DEG_ARC)
      MSB05SQ = MSB05*MSB05
C
      OFF = 0
      DO K = 1,NK
         IF ( IREL.LT.2 ) THEN
            L = K - 1
            J = L
         ELSE
            L = K/2
            IF ( L*2.EQ.K ) THEN
               J = L - 0.5D0
            ELSE
               J = L + 0.5D0
            END IF
         END IF
C
         NMUE = NINT(2*J+1)
C
         DO IM2 = 1,NMUE
            M2 = -J + (IM2-1.0D0)
            EMIM2A = CDEXP(-CI*M2*ALFDEG*DEG_ARC)
C
            DO IM1 = 1,NMUE
               M1 = -J + (IM1-1.0D0)
               EMIM1G = CDEXP(-CI*M1*GAMDEG*DEG_ARC)
C
               IF ( DABS(BETDEG).LT.1D-8 ) THEN
C
                  IF ( IM1.EQ.IM2 ) THEN
                     RSUM = 1.0D0
                  ELSE
                     RSUM = 0.0D0
                  END IF
C
               ELSE IF ( DABS(BETDEG-180).LT.1D-8 ) THEN
C
                  IF ( ABS(M1+M2).LT.1D-8 ) THEN
                     RSUM = (-1.0D0)**NINT(-J+M1)
                  ELSE
                     RSUM = 0.0D0
                  END IF
C
               ELSE
C
                  SLOW = MAX(0,NINT(M1-M2))
                  SHIGH = MIN(NINT(J-M2),NINT(J+M1))
                  CB05PW = CB05**NINT(2*J+M1-M2-2*SLOW+2)
                  MSB05PW = MSB05**NINT(M2-M1+2*SLOW-2)
                  DOM = (-1.0D0)**(SLOW-1)
     &                  *DSQRT(RFAC(J+M1)*RFAC(J-M1)*RFAC(J+M2)
     &                  *RFAC(J-M2))
                  RSUM = 0.0D0
C
                  DO S = SLOW,SHIGH
                     DOM = -DOM
                     NUM = FACT(S)*RFAC(J-M2-S)*RFAC(J+M1-S)
     &                     *RFAC(M2-M1+S)
                     CB05PW = CB05PW/CB05SQ
                     MSB05PW = MSB05PW*MSB05SQ
                     RSUM = RSUM + (DOM/NUM)*CB05PW*MSB05PW
                  END DO
               END IF
C
               ROT(OFF+IM2,OFF+IM1) = EMIM1G*RSUM*EMIM2A
            END DO
C
         END DO
C
         OFF = OFF + NMUE
      END DO
C
      END
