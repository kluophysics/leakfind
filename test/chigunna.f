C*==chigunna.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIGUNNA(IOTMP,IPRINT,DATSET,LDATSET,TXT_T,LTXT_T,
     &                    RHOCHR,RHOSPN,EMM,ENM,ENN,Z,IT,JTOP,NRMAX,
     &                    NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the curvature of Exc with respect       *
C   *  to the spin and charge densities at s >= 0                      *
C   *  after O GUNNARSSON J. Phys. F 6, 587 (1976)                     *
C   *  the Gunnarsson exchange-correlation potential Eq. (2.4)         *
C   *                                                                  *
C   *                   2       (          1     DELTA*S    )          *
C   *     V(+-) = - -----------*(BETHA +- ---*------------- )          *
C   *               PI*ALPHA*RS (          3   1 +- GAMMA*S )          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--CHIGUNNA20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*80 DATSET
      INTEGER IOTMP,IPRINT,IT,JTOP,LDATSET,NRMAX,NTMAX,Z
      REAL*8 EMM(NRMAX),ENM(NRMAX),ENN(NRMAX),RHOCHR(NRMAX),
     &       RHOSPN(NRMAX)
      INTEGER LTXT_T(NTMAX)
      CHARACTER*8 TXT_T(NTMAX)
C
C Local variables
C
      REAL*8 ALPHA,COEFF,DELTA,GAMMA,GAMMAS2,RHO,RHOS,RS,S,TERM1,TERM2
      CHARACTER*80 FILNAM
      INTEGER I,LFN
      LOGICAL NEGVALS
C
C*** End of declarations rewritten by SPAG
C
      NEGVALS = .FALSE.
C
      IF ( IPRINT.GT.0 ) THEN
         CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_GUN_A-xc.dat',
     &                 13,FILNAM,LFN,1,IOTMP,'A-xc file ',10,NTMAX)
         WRITE (IOTMP,99001)
      END IF
C
      IF ( Z.GT.0 ) THEN
         DO I = 1,JTOP
            RHO = RHOCHR(I)/(4.0D0*PI)
            RHOS = RHOSPN(I)/(4.0D0*PI)
C
C *** IF (RHO .GT. 1.725D-5) THEN EMM(I) < 0.0 ==> EMM(I) = 0.0
C *** BUT FIRST CHECK FOR NUMERICALLY TOO SMALL RHO
C
            IF ( RHO.GE.1D-20 ) THEN
C
               S = RHOS/RHO
               RS = (3D0/(4D0*PI*RHO))**(1D0/3D0)
               ALPHA = (4D0/(9D0*PI))**(1D0/3D0)
               DELTA = 1.0D0 - 0.036D0*RS - 1.36D0*RS/(1.0D0+10.0D0*RS)
               GAMMA = 0.297D0
               GAMMAS2 = (GAMMA*S)**2
               COEFF = (-2D0/3D0)*(1.D0/(PI*ALPHA*RS))
               TERM1 = (1.D0+GAMMAS2)/(1.D0-GAMMAS2)
               EMM(I) = -COEFF*DELTA/RHO*TERM1/(1.D0-GAMMAS2)
               TERM1 = 2.D0*GAMMA*S/(1.D0-GAMMAS2)
               ENM(I) = COEFF*DELTA/RHO*TERM1/(1.D0-GAMMAS2)
               TERM1 = 1.D0 + 0.0545D0*11.4D0/(RS+11.4D0)
               TERM2 = (1.D0/3.D0)*GAMMA*S*S/(1.D0-GAMMAS2)
     &                 - 2.D0*DELTA*GAMMA*S*S/(1.D0-GAMMAS2)**2
               ENN(I) = COEFF/RHO*(TERM1-TERM2)
               IF ( EMM(I).LT.0.0D0 ) THEN
                  EMM(I) = 0.0D0
                  NEGVALS = .TRUE.
               END IF
C
               IF ( IPRINT.GT.0 ) WRITE (IOTMP,'(1I4,1X,6E12.4)') I,
     &              RHOCHR(I),RS,EMM(I),COEFF
            END IF
         END DO
      END IF
C
      IF ( NEGVALS ) WRITE (6,*) 
     &       'WARNING !!! NEGATIVE VALUES FOR STONER KERNEL IN CHIGUNNA'
     &       ,' FOR IT=',IT
C
      CLOSE (IOTMP)
C
99001 FORMAT ('# FILE fort.(520+IT)',/,
     &        '# GUNNARSSON: I,RHOCHR,RS,STONER,COEFF')
      END
