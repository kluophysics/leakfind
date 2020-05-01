C*==sig_spin_checkme.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG_SPIN_CHECKME(ERYDA,ERYDB,CCF,IFILA,IFILB)
C   ********************************************************************
C   *  Check MEs                                                       *
C   *                                                                  *
C   *     + make sanity tests on various symmetry properties of MEs    *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SIG,ONLY:NSPR,LIST_ISPR,NSPINPROJ
      USE MOD_TYPES,ONLY:CTL,NT,NTMAX
      USE MOD_MPI,ONLY:MPI
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--SIG_SPIN_CHECKME16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG_SPIN_CHECKME')
      INTEGER K_CALC_ME
      PARAMETER (K_CALC_ME=3)
C
C Dummy arguments
C
      LOGICAL CCF
      COMPLEX*16 ERYDA,ERYDB
      INTEGER IFILA,IFILB
C
C Local variables
C
      REAL*8 C
      COMPLEX*16 DELMSST(NKMMAX,NKMMAX,3,NTMAX),DUMMY(:,:,:,:),
     &           DUMMY1(:,:,:),DUMMY2(:,:,:,:),YBARAB(:,:,:,:),
     &           YBARBA(:,:,:,:),YIRR2AB(:,:,:,:,:),YIRR2BA(:,:,:,:,:),
     &           YIRR3AB(:,:,:,:,:),YIRR3BA(:,:,:,:,:),
     &           YIRR4AB(:,:,:,:,:),YIRR4BA(:,:,:,:,:),YREGAB(:,:,:,:),
     &           YREGBA(:,:,:,:)
      INTEGER I,IA_ERR,INDPOLREV(3),IPOL,IPOL_SPIN,ISP,ISPR,IT,IXX,J,
     &        JPOL,M,SIGNPOL(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA INDPOLREV/3,2,1/,SIGNPOL/ + 1, - 1, + 1/
C
      ALLOCATABLE YBARAB,YIRR4BA,YREGAB,YREGBA
      ALLOCATABLE YBARBA,YIRR2AB,YIRR2BA,YIRR3AB,YIRR3BA,YIRR4AB
      ALLOCATABLE DUMMY,DUMMY1,DUMMY2
C ---------------------------------------------------------------- TESTS
C
      M = NKMMAX
      ALLOCATE (YBARAB(M,M,3,NSPINPROJ),YBARBA(M,M,3,NSPINPROJ))
      ALLOCATE (YIRR2AB(M,M,3,3,NSPINPROJ),YIRR2BA(M,M,3,3,NSPINPROJ))
      ALLOCATE (YIRR3AB(M,M,3,3,NSPINPROJ),YIRR3BA(M,M,3,3,NSPINPROJ))
      ALLOCATE (YIRR4AB(M,M,3,3,NSPINPROJ),YIRR4BA(M,M,3,3,NSPINPROJ))
      ALLOCATE (YREGAB(M,M,3,NSPINPROJ),YREGBA(M,M,3,NSPINPROJ))
      ALLOCATE (DUMMY(NKMMAX,NKMMAX,3,3),DUMMY1(NKMMAX,NKMMAX,3))
      ALLOCATE (DUMMY2(NKMMAX,NKMMAX,3,3))
C
C=======================================================================
C                     calculate angular matrix elements
C=======================================================================
C
      CALL AME_ALF_SIG
C
      CALL AME_NAB_SIG
C
      CALL AME_NAB
C
      CALL AME_INIT('ADA',0)
C
C***********************************************************************
C symmetry properties of angular matrix-elements coming from ME_ALF_ALF
C***********************************************************************
C
      WRITE (6,'(60("-"),/,5X,A,2F12.6)') 'ERYDA',ERYDA
      WRITE (6,'(5X,A,2F12.6)') 'ERYDB',ERYDB
      WRITE (6,'(5X,A,2F12.6,/,60("-"))') 'DELE ',ERYDB - ERYDA
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99005)
     &                         ROUTINE(1:LEN_TRIM(ROUTINE))
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      DO IT = 1,NT
C
         WRITE (6,99001) IT
C
         C = CTL(IT,1)
C
C-------------------------------------------------- energies E_b and E_a
C
         YREGBA = C0
         YBARBA = C0
C
C-------------------------------------------------- energies E_a and E_b
C
         YREGAB = C0
         YBARAB = C0
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR ME_ALF_ALF
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
         IXX = 1
C
C-------------------------------------------------- energies E_b and E_a
C
         CALL ME_ALF_ALF(1,NKM,IFILB,ERYDB,1,NKM,IFILA,ERYDA,CCF,IT,
     &                   YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                   YIRR2BA(1,1,1,1,IXX),YIRR3BA(1,1,1,1,IXX),
     &                   YIRR4BA(1,1,1,1,IXX),C,K_CALC_ME)
C
C-------------------------------------------------- energies E_a and E_b
C
         CALL ME_ALF_ALF(1,NKM,IFILA,ERYDA,1,NKM,IFILB,ERYDB,CCF,IT,
     &                   YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                   YIRR2AB(1,1,1,1,IXX),YIRR3AB(1,1,1,1,IXX),
     &                   YIRR4AB(1,1,1,1,IXX),C,K_CALC_ME)
C
         CALL ME_CHECK_STEP1(IT,IXX,'ME_ALF_ALF   ',YBARBA,YIRR2AB,
     &                       YIRR3BA,YIRR4AB,YIRR4BA,YREGAB,YREGBA,
     &                       DUMMY,DUMMY1)
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR MENAB
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
         DO ISPR = 2,NSPR
            ISP = LIST_ISPR(ISPR)
            IF ( ISP.EQ.9 ) THEN
               IXX = 9
C
               CALL ME_NAB_NAB(1,NKM,IFILB,ERYDB,1,NKM,IFILA,ERYDA,CCF,
     &                         IT,YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                         YIRR2BA(1,1,1,1,IXX),YIRR3BA(1,1,1,1,IXX)
     &                         ,YIRR4BA(1,1,1,1,IXX),C,K_CALC_ME)
C
C-------------------------------------------------- energies E_a and E_b
C
               CALL ME_NAB_NAB(1,NKM,IFILA,ERYDA,1,NKM,IFILB,ERYDB,CCF,
     &                         IT,YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                         YIRR2AB(1,1,1,1,IXX),YIRR3AB(1,1,1,1,IXX)
     &                         ,YIRR4AB(1,1,1,1,IXX),C,K_CALC_ME)
C
               WRITE (6,99002) ROUTINE,'MENAB'
C
C***********************************************************************
C  TEST 1     MZBZA_m(YBARBA) vs MZAZB_mbar(YREGBA)
C***********************************************************************
C
               DO IPOL = 1,3
                  DUMMY1(1:NKM,1:NKM,IPOL) = SIGNPOL(IPOL)
     &               *TRANSPOSE(YREGBA(1:NKM,1:NKM,INDPOLREV(IPOL),9))
               END DO
C
               CALL CMPMAT3(3,YBARBA(1,1,1,9),DUMMY1,'TEST1',IXX,IT)
C
C***********************************************************************
C  TEST 2   ( MZBZA(YREGAB) vs MZAZB(YBARBA) )
C***********************************************************************
C
               CALL CMPMAT3(3,YBARBA,YREGAB,'TEST2',IXX,IT)
C
C***********************************************************************
C  TEST 3    MIRR_3_mn(YIRR3BA) vs (MIRR_2_nbar_mbar(YIRR2AB))^T
C***********************************************************************
C
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &                  *SIGNPOL(JPOL)
     &                  *TRANSPOSE(YIRR2AB(1:NKM,1:NKM,INDPOLREV(IPOL),
     &                  INDPOLREV(JPOL),9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,9),'TEST3',IXX,IT)
C
C***********************************************************************
C  TEST 6    MIRR_4_mn(YIRR4AB) =(MIRR_4_mbar_nbar(YIRR4BA))^T
C***********************************************************************
C
               DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &                  *SIGNPOL(JPOL)
     &                  *TRANSPOSE(YIRR4BA(1:NKM,1:NKM,INDPOLREV(IPOL),
     &                  INDPOLREV(JPOL),9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR4AB(1,1,1,1,9),'TEST6',IXX,IT)
            END IF
         END DO
C
C=======================================================================
C                   overwrite Js withs Zs
C=======================================================================
C
         WRITE (6,99004) ROUTINE
         CALL WF_J_TO_Z(IFILA,1,NKM,IT)
         CALL WF_J_TO_Z(IFILB,1,NKM,IT)
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR ME_ALF_ALF
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
         IXX = 1
C
         CALL ME_ALF_ALF(1,NKM,IFILB,ERYDB,1,NKM,IFILA,ERYDA,CCF,IT,
     &                   YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                   YIRR2BA(1,1,1,1,IXX),YIRR3BA(1,1,1,1,IXX),
     &                   YIRR4BA(1,1,1,1,IXX),C,K_CALC_ME)
C
C-------------------------------------------------- energies E_a and E_b
C
         CALL ME_ALF_ALF(1,NKM,IFILA,ERYDA,1,NKM,IFILB,ERYDB,CCF,IT,
     &                   YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                   YIRR2AB(1,1,1,1,IXX),YIRR3AB(1,1,1,1,IXX),
     &                   YIRR4AB(1,1,1,1,IXX),C,K_CALC_ME)
C
         WRITE (6,99002) ROUTINE,'ME_ALF_ALF'
C
C***********************************************************************
C  TEST 1     MZBZA_m(YBARBA) vs MZAZB_mbar(YREGBA)
C***********************************************************************
C
         DO IPOL = 1,3
            DUMMY1(1:NKM,1:NKM,IPOL) = SIGNPOL(IPOL)
     &                                 *TRANSPOSE(YREGBA(1:NKM,1:NKM,
     &                                 INDPOLREV(IPOL),1))
         END DO
C
C
         CALL CMPMAT3(3,YBARBA(1,1,1,1),DUMMY1,'TEST1',IXX,IT)
C
C***********************************************************************
C  TEST 2   ( MZBZA(YREGAB) vs MZAZB(YBARBA) )
C***********************************************************************
C
         CALL CMPMAT3(3,YBARBA,YREGAB,'TEST2',IXX,IT)
C
C***********************************************************************
C  TEST 3    MIRR_3_mn(YIRR3BA) vs (MIRR_2_nbar_mbar(YIRR2AB))^T
C***********************************************************************
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &            *SIGNPOL(JPOL)
     &            *TRANSPOSE(YIRR2AB(1:NKM,1:NKM,INDPOLREV(IPOL),
     &            INDPOLREV(JPOL),1))
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,1),'TEST3',IXX,IT)
C
C***********************************************************************
C  TEST 4    MIRR_2_mn(YIRR2BA) = MZBZA_m*MZAZB_n (replace J with Z)
C***********************************************************************
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &            = MATMUL(YBARBA(1:NKM,1:NKM,IPOL,1),
     &            YREGBA(1:NKM,1:NKM,JPOL,1))
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,YIRR2BA(1,1,1,1,1),'TEST4',IXX,IT)
C
C***********************************************************************
C  TEST 5    MIRR_3_mn(YIRR3BA) = MZAZB_n*MZBZA_m (replace J with Z)
C***********************************************************************
C
         DO IPOL = 1,3
            DO JPOL = 1,3
               DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &            = MATMUL(YREGBA(1:NKM,1:NKM,JPOL,1),
     &            YBARBA(1:NKM,1:NKM,IPOL,1))
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,1),'TEST5',IXX,IT)
C
C***********************************************************************
C  TEST 6    MIRR_4_mn(YIRR4AB) =(MIRR_4_mbar_nbar(YIRR4BA))^T
C***********************************************************************
C
         DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
         DO IPOL = 1,3
            DO JPOL = 1,3
               DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &            *SIGNPOL(JPOL)
     &            *TRANSPOSE(YIRR4BA(1:NKM,1:NKM,INDPOLREV(IPOL),
     &            INDPOLREV(JPOL),1))
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,YIRR4AB(1,1,1,1,1),'TEST6',IXX,IT)
C
C***********************************************************************
C  TEST 7    \sum_lam'    MIRR_4_mn(lam',lam) = MIRR_2_mn(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
         DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
         DO I = 1,NKM
            DO J = 1,NKM
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(I,I,IPOL,JPOL) = DUMMY(I,I,IPOL,JPOL)
     &                  + YIRR4BA(J,I,IPOL,JPOL,1)
                  END DO
               END DO
            END DO
         END DO
C
         DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
         DO I = 1,NKM
            DO IPOL = 1,3
               DO JPOL = 1,3
                  DUMMY2(I,I,IPOL,JPOL) = YIRR2BA(I,I,IPOL,JPOL,1)
               END DO
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,DUMMY2,'TEST7',IXX,IT)
C
C***********************************************************************
C  TEST 8  \sum_lam' MIRR_4_mn(lam',lam) = MIRR_3_mbar_nbar(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
         DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
         DO I = 1,NKM
            DO J = 1,NKM
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(I,I,IPOL,JPOL) = DUMMY(I,I,IPOL,JPOL)
     &                  + YIRR4BA(J,I,IPOL,JPOL,1)
                  END DO
               END DO
            END DO
         END DO
C
         DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
         DO I = 1,NKM
            DO IPOL = 1,3
               DO JPOL = 1,3
                  DUMMY2(I,I,IPOL,JPOL) = YIRR3AB(I,I,JPOL,IPOL,1)
               END DO
            END DO
         END DO
C
         CALL CMPMAT9(DUMMY,DUMMY2,'TEST8',IXX,IT)
C
C-----------------------------------------------------------------------
         DO ISPR = 2,NSPR
            ISP = LIST_ISPR(ISPR)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR ME_SPIN_CURR_ALF
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
            IF ( ISP.EQ.2 ) THEN
               IXX = 2
C
               CALL ME_SPIN_ALF_ALF(IFILB,ERYDB,IFILA,ERYDA,IT,
     &                              YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                              YIRR2BA(1,1,1,1,IXX),
     &                              YIRR3BA(1,1,1,1,IXX),
     &                              YIRR4BA(1,1,1,1,IXX),C)
C
C-------------------------------------------------- energies E_a and E_b
C
               CALL ME_SPIN_ALF_ALF(IFILA,ERYDA,IFILB,ERYDB,IT,
     &                              YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                              YIRR2AB(1,1,1,1,IXX),
     &                              YIRR3AB(1,1,1,1,IXX),
     &                              YIRR4AB(1,1,1,1,IXX),C)
C
               WRITE (6,99002) ROUTINE,'ME_SPIN_CURR_ALF'
C
C***********************************************************************
C  TEST 1    MIRR_2_mn(YIRR2BA) = MZBZA_m*MZAZB_n     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 2,4
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN),
     &                     YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN))
                     END DO
                  END DO
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY,YIRR2BA(1,1,1,1,IPOL_SPIN),'TEST1',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 2    MIRR_3_mn(YIRR3BA) = MZAZB_n*MZBZA_m (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 2,4
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN),
     &                     YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN))
                     END DO
                  END DO
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,IPOL_SPIN),'TEST2',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 3    \sum_lam' MIRR_4_mn(lam',lam) = MIRR_2_mn(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 2,4
                  DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO J = 1,NKM
                        DO IPOL = 1,3
                           DO JPOL = 1,3
                              DUMMY(I,I,IPOL,JPOL)
     &                           = DUMMY(I,I,IPOL,JPOL)
     &                           + YIRR4BA(J,I,IPOL,JPOL,IPOL_SPIN)
                           END DO
                        END DO
                     END DO
                  END DO
C
                  DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO IPOL = 1,3
                        DO JPOL = 1,3
                           DUMMY2(I,I,IPOL,JPOL)
     &                        = YIRR2BA(I,I,IPOL,JPOL,IPOL_SPIN)
                        END DO
                     END DO
                  END DO
C
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY2,DUMMY,'TEST3',IPOL_SPIN,IT)
C
               END DO
C
            END IF
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR ME_SPIN_CURR_NAB
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
            IF ( ISP.EQ.5 ) THEN
               IXX = 5
C
               CALL ME_SPIN_NAB_ALF(IFILB,ERYDB,IFILA,ERYDA,IT,
     &                              YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                              YIRR2BA(1,1,1,1,IXX),
     &                              YIRR3BA(1,1,1,1,IXX),
     &                              YIRR4BA(1,1,1,1,IXX),C)
C
C-------------------------------------------------- energies E_a and E_b
C
               CALL ME_SPIN_NAB_ALF(IFILA,ERYDA,IFILB,ERYDB,IT,
     &                              YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                              YIRR2AB(1,1,1,1,IXX),
     &                              YIRR3AB(1,1,1,1,IXX),
     &                              YIRR4AB(1,1,1,1,IXX),C)
C
               WRITE (6,99002) ROUTINE,'ME_SPIN_CURR_NAB'
C
C***********************************************************************
C  TEST 1    MIRR_2_mn(YIRR2BA) = MZBZA_m*MZAZB_n     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 5,7
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN),
     &                     YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN))
                     END DO
                  END DO
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY,YIRR2BA(1,1,1,1,IPOL_SPIN),'TEST1',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 2    MIRR_3_mn(YIRR3BA) = MZAZB_n*MZBZA_m (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 5,7
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN),
     &                     YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN))
                     END DO
                  END DO
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,IPOL_SPIN),'TEST2',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 3    \sum_lam' MIRR_4_mn(lam',lam) = MIRR_2_mn(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 5,7
                  DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO J = 1,NKM
                        DO IPOL = 1,3
                           DO JPOL = 1,3
                              DUMMY(I,I,IPOL,JPOL)
     &                           = DUMMY(I,I,IPOL,JPOL)
     &                           + YIRR4BA(J,I,IPOL,JPOL,IPOL_SPIN)
                           END DO
                        END DO
                     END DO
                  END DO
C
                  DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO IPOL = 1,3
                        DO JPOL = 1,3
                           DUMMY2(I,I,IPOL,JPOL)
     &                        = YIRR2BA(I,I,IPOL,JPOL,IPOL_SPIN)
                        END DO
                     END DO
                  END DO
C
                  WRITE (6,99003) IPOL_SPIN
                  CALL CMPMAT9(DUMMY2,DUMMY,'TEST3',IPOL_SPIN,IT)
               END DO
C
            END IF
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR ME_SOT
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
            IF ( ISP.EQ.8 ) THEN
               IXX = 8
C
               CALL ME_SOT_ALF(IFILB,ERYDB,IFILA,ERYDA,IT,
     &                         YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                         YIRR2BA(1,1,1,1,IXX),YIRR3BA(1,1,1,1,IXX)
     &                         ,YIRR4BA(1,1,1,1,IXX),C)
C
C-------------------------------------------------- energies E_a and E_b
C
               CALL ME_SOT_ALF(IFILA,ERYDA,IFILB,ERYDB,IT,
     &                         YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                         YIRR2AB(1,1,1,1,IXX),YIRR3AB(1,1,1,1,IXX)
     &                         ,YIRR4AB(1,1,1,1,IXX),C)
C
               WRITE (6,99002) ROUTINE,'ME_SOT'
C
               IF ( ERYDA.EQ.ERYDB .AND. ERYDA.EQ.DREAL(ERYDA) ) THEN
                  CALL XCPLTENME(DELMSST)
                  DUMMY1(1:NKM,1:NKM,1:3) = DELMSST(1:NKM,1:NKM,1:3,IT)
                  CALL CMPMAT3(3,YBARBA(1,1,1,IXX),DUMMY1,'TEST0',IXX,
     &                         IT)
               ELSE
                  WRITE (6,*)
                  WRITE (6,*) 
     &             'Comparison to MEs of exchange tensor cannot be done'
               END IF
C
C***********************************************************************
C  TEST 1    MIRR_2_mn(YIRR2BA) = MZBZA_m*MZAZB_n     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 8,8
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN),
     &                     YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN))
                     END DO
                  END DO
C
                  CALL CMPMAT9(DUMMY,YIRR2BA(1,1,1,1,IPOL_SPIN),'TEST1',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 2    MIRR_3_mn(YIRR3BA) = MZAZB_n*MZBZA_m (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 8,8
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                     = MATMUL(YREGBA(1:NKM,1:NKM,JPOL,IPOL_SPIN),
     &                     YBARBA(1:NKM,1:NKM,IPOL,IPOL_SPIN))
                     END DO
                  END DO
C
                  CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,IPOL_SPIN),'TEST2',
     &                         IPOL_SPIN,IT)
               END DO
C
C***********************************************************************
C  TEST 3    \sum_lam' MIRR_4_mn(lam',lam) = MIRR_2_mn(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
               DO IPOL_SPIN = 8,8
                  DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO J = 1,NKM
                        DO IPOL = 1,3
                           DO JPOL = 1,3
                              DUMMY(I,I,IPOL,JPOL)
     &                           = DUMMY(I,I,IPOL,JPOL)
     &                           + YIRR4BA(J,I,IPOL,JPOL,IPOL_SPIN)
                           END DO
                        END DO
                     END DO
                  END DO
C
                  DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
                  DO I = 1,NKM
                     DO IPOL = 1,3
                        DO JPOL = 1,3
                           DUMMY2(I,I,IPOL,JPOL)
     &                        = YIRR2BA(I,I,IPOL,JPOL,IPOL_SPIN)
                        END DO
                     END DO
                  END DO
C
                  CALL CMPMAT9(DUMMY2,DUMMY,'TEST3',IPOL_SPIN,IT)
               END DO
            END IF
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                  INTERNAL TESTS FOR MENAB
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C-------------------------------------------------- energies E_b and E_a
C
            IF ( ISP.EQ.9 ) THEN
               IXX = 9
C
               CALL ME_NAB_NAB(1,NKM,IFILB,ERYDB,1,NKM,IFILA,ERYDA,CCF,
     &                         IT,YREGBA(1,1,1,IXX),YBARBA(1,1,1,IXX),
     &                         YIRR2BA(1,1,1,1,IXX),YIRR3BA(1,1,1,1,IXX)
     &                         ,YIRR4BA(1,1,1,1,IXX),C,K_CALC_ME)
C
C-------------------------------------------------- energies E_a and E_b
C
               CALL ME_NAB_NAB(1,NKM,IFILA,ERYDA,1,NKM,IFILB,ERYDB,CCF,
     &                         IT,YREGAB(1,1,1,IXX),YBARAB(1,1,1,IXX),
     &                         YIRR2AB(1,1,1,1,IXX),YIRR3AB(1,1,1,1,IXX)
     &                         ,YIRR4AB(1,1,1,1,IXX),C,K_CALC_ME)
C
               WRITE (6,99002) ROUTINE,'MENAB'
C
C***********************************************************************
C  TEST 1     MZBZA_m(YBARBA) vs MZAZB_mbar(YREGBA)
C***********************************************************************
C
               DO IPOL = 1,3
                  DUMMY1(1:NKM,1:NKM,IPOL) = SIGNPOL(IPOL)
     &               *TRANSPOSE(YREGBA(1:NKM,1:NKM,INDPOLREV(IPOL),9))
               END DO
C
               CALL CMPMAT3(3,YBARBA(1,1,1,9),DUMMY1,'TEST1',IXX,IT)
C
C***********************************************************************
C  TEST 2   ( MZBZA(YREGAB) vs MZAZB(YBARBA) )
C***********************************************************************
C
               CALL CMPMAT3(3,YBARBA,YREGAB,'TEST2',IXX,IT)
C
C***********************************************************************
C  TEST 3    MIRR_3_mn(YIRR3BA) vs (MIRR_2_nbar_mbar(YIRR2AB))^T
C***********************************************************************
C
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &                  *SIGNPOL(JPOL)
     &                  *TRANSPOSE(YIRR2AB(1:NKM,1:NKM,INDPOLREV(IPOL),
     &                  INDPOLREV(JPOL),9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,9),'TEST3',IXX,IT)
C
C***********************************************************************
C  TEST 4    MIRR_2_mn(YIRR2BA) = MZBZA_m*MZAZB_n (replace J with Z)
C***********************************************************************
C
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                  = MATMUL(YBARBA(1:NKM,1:NKM,IPOL,9),
     &                  YREGBA(1:NKM,:NKM,JPOL,9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR2BA(1,1,1,1,9),'TEST4',IXX,IT)
C
C***********************************************************************
C  TEST 5    MIRR_3_mn(YIRR3BA) = MZAZB_n*MZBZA_m (replace J with Z)
C***********************************************************************
C
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL)
     &                  = MATMUL(YREGBA(1:NKM,1:NKM,JPOL,9),
     &                  YBARBA(1:NKM,1:NKM,IPOL,9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,9),'TEST5',IXX,IT)
C
C***********************************************************************
C  TEST 6    MIRR_4_mn(YIRR4AB) =(MIRR_4_mbar_nbar(YIRR4BA))^T
C***********************************************************************
C
               DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
               DO IPOL = 1,3
                  DO JPOL = 1,3
                     DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)
     &                  *SIGNPOL(JPOL)
     &                  *TRANSPOSE(YIRR4BA(1:NKM,1:NKM,INDPOLREV(IPOL),
     &                  INDPOLREV(JPOL),9))
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,YIRR4AB(1,1,1,1,9),'TEST6',IXX,IT)
C
C***********************************************************************
C  TEST 7    \sum_lam'    MIRR_4_mn(lam',lam) = MIRR_2_mn(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
               DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
               DO I = 1,NKM
                  DO J = 1,NKM
                     DO IPOL = 1,3
                        DO JPOL = 1,3
                           DUMMY(I,I,IPOL,JPOL) = DUMMY(I,I,IPOL,JPOL)
     &                        + YIRR4BA(J,I,IPOL,JPOL,9)
                        END DO
                     END DO
                  END DO
               END DO
C
               DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
               DO I = 1,NKM
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY2(I,I,IPOL,JPOL) = YIRR2BA(I,I,IPOL,JPOL,9)
                     END DO
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,DUMMY2,'TEST7',IXX,IT)
C
C***********************************************************************
C  TEST 8  \sum_lam' MIRR_4_mn(lam',lam) = MIRR_3_mbar_nbar(lam,lam)
C                                                     (replace J with Z)
C***********************************************************************
C
               DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
               DO I = 1,NKM
                  DO J = 1,NKM
                     DO IPOL = 1,3
                        DO JPOL = 1,3
                           DUMMY(I,I,IPOL,JPOL) = DUMMY(I,I,IPOL,JPOL)
     &                        + YIRR4BA(J,I,IPOL,JPOL,9)
                        END DO
                     END DO
                  END DO
               END DO
C
               DUMMY2(1:NKM,1:NKM,1:3,1:3) = C0
               DO I = 1,NKM
                  DO IPOL = 1,3
                     DO JPOL = 1,3
                        DUMMY2(I,I,IPOL,JPOL) = YIRR3AB(I,I,JPOL,IPOL,9)
                     END DO
                  END DO
               END DO
C
               CALL CMPMAT9(DUMMY,DUMMY2,'TEST8',IXX,IT)
C
C-----------------------------------------------------------------------
            END IF
         END DO
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END DO ! IT
C
      IF ( MPI ) CALL MPI_FINALIZE(IA_ERR)
C
      WRITE (6,*) '          ALL DONE'
      STOP
C
99001 FORMAT (///,' J - matrix elements for atom type   IT =',I4)
99002 FORMAT (//,60('-'),/,A18,' Internal tests: ',A,/,60('-'))
99003 FORMAT ('IPOL_SPIN=',I1)
99004 FORMAT (//,60('-'),/,A18,' will overwrite Js with Zs',/,60('-'))
99005 FORMAT ('# BUILDBOT: ',A)
      END
C*==cmpmat3.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE CMPMAT3(N,A,B,TXT,IYY,IT)
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--CMPMAT3854
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT,IYY,N
      CHARACTER*5 TXT
      COMPLEX*16 A(NKMMAX,NKMMAX,N),B(NKMMAX,NKMMAX,N)
C
C Local variables
C
      REAL*8 ADEL,ADELMAX,FLAG,RDEL,RDELMAX
      INTEGER I,J,K
C
C*** End of declarations rewritten by SPAG
C
      RDELMAX = 0D0
      ADELMAX = 0D0
C
      DO K = 1,N
C
         DO I = 1,NKM
            DO J = 1,NKM
C
               ADEL = ABS(A(I,J,K)-B(I,J,K))
               ADELMAX = MAX(ADEL,ADELMAX)
C
               IF ( ABS(A(I,J,K)).GT.1D-12 ) THEN
                  RDEL = ABS(1D0-B(I,J,K)/A(I,J,K))
                  RDELMAX = MAX(RDEL,RDELMAX)
               END IF
            END DO
         END DO
      END DO
C
      IF ( RDELMAX.GT.1D-9 ) WRITE (6,99002)
C
      IF ( RDELMAX.GT.1D-9 ) THEN
         FLAG = 1D0
      ELSE
         FLAG = 0D0
      END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99003) IYY,IT,TXT,FLAG
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      WRITE (6,99001) TXT,ADELMAX,RDELMAX
99001 FORMAT (5X,A,'    abs: ',E15.6,'    rel: ',E15.6)
99002 FORMAT (5X,'Test fails ')
99003 FORMAT ('# BUILDBOT: ISPINPROJ=',I2,2X,'IT=',I2,2X,A,/,(1PE22.14))
      END
C*==cmpmat9.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
C
C
      SUBROUTINE CMPMAT9(A,B,TXT,IYY,IT)
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--CMPMAT9926
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT,IYY
      CHARACTER*5 TXT
      COMPLEX*16 A(NKMMAX,NKMMAX,3,3),B(NKMMAX,NKMMAX,3,3)
C
C Local variables
C
      REAL*8 ADEL,ADELMAX,FLAG,RDEL,RDELMAX
      INTEGER I,J,K,L
C
C*** End of declarations rewritten by SPAG
C
      RDELMAX = 0D0
      ADELMAX = 0D0
C
      DO K = 1,3
         DO L = 1,3
C
            DO I = 1,NKM
               DO J = 1,NKM
C
                  ADEL = ABS(A(I,J,K,L)-B(I,J,K,L))
C
                  IF ( ADEL.GT.1D-8 .AND. 1.EQ.0 ) THEN
                     WRITE (6,*) I,J,K,L
                     WRITE (6,*) A(I,J,K,L)
                     WRITE (6,*) B(I,J,K,L)
                     WRITE (6,*) A(I,J,K,L)/B(I,J,K,L)
                     WRITE (6,*)
                  END IF
C
                  ADELMAX = MAX(ADEL,ADELMAX)
C
                  IF ( ABS(A(I,J,K,L)).GT.1D-12 ) THEN
                     RDEL = ABS(1D0-B(I,J,K,L)/A(I,J,K,L))
                     RDELMAX = MAX(RDEL,RDELMAX)
                  END IF
               END DO
            END DO
C
         END DO
      END DO
C
      IF ( RDELMAX.GT.1D-9 ) WRITE (6,99002)
C
      IF ( RDELMAX.GT.1D-9 ) THEN
         FLAG = 1D0
      ELSE
         FLAG = 0D0
      END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99003) IYY,IT,TXT,FLAG
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      WRITE (6,99001) TXT,ADELMAX,RDELMAX
99001 FORMAT (5X,A,'    abs: ',E15.6,'    rel: ',E15.6)
99002 FORMAT (5X,'Test fails')
99003 FORMAT ('# BUILDBOT: ISPINPROJ=',I2,2X,'IT=',I2,2X,A,/,(1PE22.14))
      END
C*==wf_j_to_z.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE WF_J_TO_Z(IFIL,LAM1,LAM2,IT)
C   ********************************************************************
C   *                                                                  *
C   *     replace    J   by   Z  in the files on disc                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NT,NCPLWFMAX,IMT,IKMCPLWF
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,NCPLWF
      USE MOD_RMESH,ONLY:NRMAX,JRWS,JRCRI,FULLPOT
      IMPLICIT NONE
C*--WF_J_TO_Z1014
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WF_J_TO_Z')
C
C Dummy arguments
C
      INTEGER IFIL,IT,LAM1,LAM2
C
C Local variables
C
      INTEGER I,IA_ERR,IM,IRTOP,K,LAM
      COMPLEX*16 JF(:,:,:),JG(:,:,:),ZF(:,:,:),ZG(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZF,ZG,JF,JG
C
      WRITE (6,99001) IFIL,LAM1,LAM2,IT
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
C ======================================================================
C         replace    J   by   Z  in the files on disc
C ======================================================================
C
      IM = IMT(IT)
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      CALL WAVFUN_READ_REL(IFIL,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,IKMCPLWF)
C
      DO LAM = LAM1,LAM2
C
         WRITE (IFIL,REC=LAM+(IT-1+NT)*NKM) IT,'IRR',LAM,IRTOP,
     &          ((ZG(I,K,LAM),I=1,IRTOP),K=1,NCPLWF(LAM)),
     &          ((ZF(I,K,LAM),I=1,IRTOP),K=1,NCPLWF(LAM))
C
      END DO
C
99001 FORMAT (/,10x,'<WF_J_TO_Z> overwriting Js with Zs in',/,
     &        'in file number',i5,' ikm1',i5,' ikm2',i5,' type',i5)
      END
C*==me_check_step1.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE ME_CHECK_STEP1(IT,IXX,METYPE,YBARBA,YIRR2AB,YIRR3BA,
     &                          YIRR4AB,YIRR4BA,YREGAB,YREGBA,DUMMY,
     &                          DUMMY1)
C   ********************************************************************
C   *  Check MEs                                                       *
C   *                                                                  *
C   *     + make sanity tests on various symmetry properties of MEs    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SIG,ONLY:NSPINPROJ
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--ME_CHECK_STEP11097
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ME_CHECK_STEP1')
C
C Dummy arguments
C
      INTEGER IT,IXX
      CHARACTER*(*) METYPE
      COMPLEX*16 DUMMY(NKMMAX,NKMMAX,3,3),DUMMY1(NKMMAX,NKMMAX,3),
     &           YBARBA(NKMMAX,NKMMAX,3,NSPINPROJ),
     &           YIRR2AB(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           YIRR3BA(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           YIRR4AB(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           YIRR4BA(NKMMAX,NKMMAX,3,3,NSPINPROJ),
     &           YREGAB(NKMMAX,NKMMAX,3,NSPINPROJ),
     &           YREGBA(NKMMAX,NKMMAX,3,NSPINPROJ)
C
C Local variables
C
      INTEGER INDPOLREV(3),IPOL,JPOL,SIGNPOL(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA INDPOLREV/3,2,1/,SIGNPOL/ + 1, - 1, + 1/
C
      WRITE (6,99001) ROUTINE,METYPE
C
C***********************************************************************
C  TEST 1     MZBZA_m(YBARBA) vs MZAZB_mbar(YREGBA)
C***********************************************************************
C
      DO IPOL = 1,3
         DUMMY1(1:NKM,1:NKM,IPOL) = SIGNPOL(IPOL)
     &                              *TRANSPOSE(YREGBA(1:NKM,1:NKM,
     &                              INDPOLREV(IPOL),IXX))
      END DO
C
      CALL CMPMAT3(3,YBARBA(1,1,1,IXX),DUMMY1,'TEST1',IXX,IT)
C
C***********************************************************************
C  TEST 2   ( MZBZA(YREGAB) vs MZAZB(YBARBA) )
C***********************************************************************
C
      CALL CMPMAT3(3,YBARBA,YREGAB,'TEST2',IXX,IT)
C
C***********************************************************************
C  TEST 3    MIRR_3_mn(YIRR3BA) vs (MIRR_2_nbar_mbar(YIRR2AB))^T
C***********************************************************************
C
      DO IPOL = 1,3
         DO JPOL = 1,3
            DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)*SIGNPOL(JPOL)
     &         *TRANSPOSE
     &         (YIRR2AB(1:NKM,1:NKM,INDPOLREV(IPOL),INDPOLREV(JPOL),IXX)
     &         )
         END DO
      END DO
C
      CALL CMPMAT9(DUMMY,YIRR3BA(1,1,1,1,IXX),'TEST3',IXX,IT)
C
C***********************************************************************
C  TEST 6    MIRR_4_mn(YIRR4AB) =(MIRR_4_mbar_nbar(YIRR4BA))^T
C***********************************************************************
C
      DUMMY(1:NKM,1:NKM,1:3,1:3) = C0
      DO IPOL = 1,3
         DO JPOL = 1,3
            DUMMY(1:NKM,1:NKM,IPOL,JPOL) = SIGNPOL(IPOL)*SIGNPOL(JPOL)
     &         *TRANSPOSE
     &         (YIRR4BA(1:NKM,1:NKM,INDPOLREV(IPOL),INDPOLREV(JPOL),IXX)
     &         )
         END DO
      END DO
C
      CALL CMPMAT9(DUMMY,YIRR4AB(1,1,1,1,1),'TEST6',IXX,IT)
C
99001 FORMAT (//,60('-'),/,A18,' Internal tests: ',A,/,60('-'))
      END
