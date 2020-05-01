C*==rgntsf_setup.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE RGNTSF_SETUP(L1MAX,L2MAX,NRGNTSF_LM1,NRGNTSF)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the Gaunt coefficients    G(LM1,LM2,LM3=LM_SF)        *
C   *                                                                  *
C   *   for L1<=L1MAX, L2<=L2MAX, L3=L_SF<=NLSF-1                      *
C   *                                                                  *
C   *  used for folding the rho or potential with the shape functions  *
C   *  accordingly only those are tabulated for which any L_SF of the  *
C   *  shape functions is non-zero                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NL
      USE MOD_TYPES,ONLY:NLMFPMAX
      USE MOD_RMESH,ONLY:NM,KLMSF
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C*--RGNTSF_SETUP20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1.0D-10)
C
C Dummy arguments
C
      INTEGER L1MAX,L2MAX,NRGNTSF
      INTEGER NRGNTSF_LM1(0:NLMFPMAX)
C
C Local variables
C
      REAL*8 GAUNT_RYLM
      INTEGER I,IM,ISUM,L1,L2,L3,LM1,LM2,LM3,LSFMAX,M1,M2,M3
      LOGICAL KLSF(:)
      REAL*8 RGNT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KLSF
C
      LSFMAX = 4*(NL-1)
C
      ALLOCATE (KLSF(0:LSFMAX))
      KLSF(0:LSFMAX) = .FALSE.
      DO IM = 1,NM
         DO L3 = 0,LSFMAX
            DO M3 = -L3,L3
               LM3 = L3*L3 + L3 + M3 + 1
               IF ( KLMSF(LM3,IM).NE.0 ) KLSF(L3) = .TRUE.
            END DO
         END DO
      END DO
C
      I = 1
      DO L1 = 0,L1MAX
         DO M1 = -L1,L1
            LM1 = L1*L1 + L1 + M1 + 1
            NRGNTSF_LM1(LM1-1) = I - 1
C
            DO L3 = 0,LSFMAX
               DO M3 = -L3,L3
                  LM3 = L3*L3 + L3 + M3 + 1
                  ISUM = 0
                  DO IM = 1,NM
                     ISUM = ISUM + KLMSF(LM3,IM)
                  END DO
                  IF ( ISUM.GT.0 ) THEN
C                  IF ( KLSF(L3) ) THEN
C
                     DO L2 = 0,L2MAX
                        DO M2 = -L2,L2
                           LM2 = L2*(L2+1) + M2 + 1
C
                           RGNT = GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
C
                           IF ( ABS(RGNT).GT.TOL ) THEN
                              WRITE (IOTMP) LM1,LM2,LM3,RGNT
                              I = I + 1
                           END IF
C
                        END DO
                     END DO
C
                  END IF
               END DO
            END DO
C
         END DO
      END DO
C
      NRGNTSF_LM1(LM1) = I - 1
      NRGNTSF = NRGNTSF_LM1(LM1)
C
      WRITE (6,99001) L1MAX,L2MAX,LSFMAX,NRGNTSF
99001 FORMAT (/,1X,79('*'),/,33X,'<RGNTSF_SETUP>',/,1X,79('*'),//,10X,
     &        'for L1MAX =',I3,'   L2MAX =',I3,'   LSFMAX =',I3,/,10X,
     &        'number of Gaunt coefficients RGNTSF created:',I8,/)
      END
