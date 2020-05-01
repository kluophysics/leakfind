C*==init_umat.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE INIT_UMAT(SVEC,UMAT,ERYD)
C   ********************************************************************
C   *                                                                  *
C   *   This subroutine calculates the transformation matrix Ull'(E)   *
C   *   for the Green's function in the shifted position               *
C   *   accoring to N. Stefanou et al PRB 36 6372 (1987) eq.(6).       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NL,NLMMAX,L_LM,NLMAX,NRGNT123TAB,NLM
      USE MOD_CONSTANTS,ONLY:C1,CI,C0,PI
      IMPLICIT NONE
C*--INIT_UMAT14
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
      COMPLEX*16 ERYD
      REAL*8 SVEC(3)
      COMPLEX*16 UMAT(NLMMAX,NLMMAX)
C
C Local variables
C
      REAL*8 ALATDUM,M1_PW_LL(:,:),RGNT(:),RYLM(:),SABS,X,Y,Z
      COMPLEX*16 CHL(:),CIPWL(:),CJL(:),CLLL,CNL(:),ZZ
      REAL*8 DNRM2
      INTEGER I,IA_ERR,IG123,IRGNT(:),L,L1,L2,L3,L3MAX,LM,LM1,LM1LM2,
     &        LM2,LM3,LMAX,LRGNT12,LRGNT123,NLM3MAX,NRGNT(:)
      LOGICAL INITIALIZE
      SAVE CIPWL,IRGNT,LRGNT12,LRGNT123,M1_PW_LL,NLM3MAX,NRGNT,RGNT,RYLM
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./,IA_ERR/0/
C
      ALLOCATABLE IRGNT,NRGNT,RGNT,CIPWL,M1_PW_LL,RYLM
      ALLOCATABLE CJL,CNL,CHL
C
      LMAX = NL - 1
      L3MAX = 2*(NL-1)
C
      ALLOCATE (CJL(0:L3MAX),CNL(0:L3MAX),CHL(0:L3MAX))
C
C***********************************************************************
C                        INITIALZE
C***********************************************************************
      IF ( INITIALIZE ) THEN
C
         WRITE (6,99001) NL
C
C-----------------------------------------------------------------------
C                 - allocate storage for GAUNT coefficients
C                 - set up modified GAUNT's
C-----------------------------------------------------------------------
C
         LRGNT123 = NRGNT123TAB(NL)
         LRGNT12 = (NL**2*(NL**2+1)/2)
C
         ALLOCATE (IRGNT(LRGNT123),NRGNT(LRGNT12),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:STRINIT->IRGNT'
         ALLOCATE (RGNT(LRGNT123),CIPWL((2*NLMAX)**2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:STRINIT->RGNT'
C
         LMAX = NL - 1
C
         ALATDUM = 2*PI
C
         CALL STRGAUNT(LMAX,ALATDUM,RGNT,NRGNT,IRGNT,CIPWL,NLMAX,IG123,
     &                 LRGNT12,LRGNT123)
C
C.......................................................................
C
         DEALLOCATE (CIPWL)
         ALLOCATE (CIPWL(-(2*NL):(2*NL)))
C
         CIPWL(0) = 1D0
         DO L = 1,2*NL
            CIPWL(L) = CIPWL(L-1)*CI
            CIPWL(-L) = 1D0/CIPWL(L)
         END DO
C
         ALLOCATE (M1_PW_LL(NLMMAX,NLMMAX))
C
         DO LM1 = 1,NLMMAX
            DO LM2 = 1,LM1
               L1 = L_LM(LM1)
               L2 = L_LM(LM2)
               M1_PW_LL(LM1,LM2) = DBLE((-1)**(L1+L2))
               M1_PW_LL(LM2,LM1) = M1_PW_LL(LM1,LM2)
            END DO
         END DO
C
C-----------------------------------------------------------------------
C              calculation of the harmonic polynomials
C-----------------------------------------------------------------------
C
         NLM3MAX = (L3MAX+1)*(L3MAX+1)
C
         ALLOCATE (RYLM(NLM3MAX))
C
         INITIALIZE = .FALSE.
C
      END IF
C***********************************************************************
C
      UMAT(1:NLMMAX,1:NLMMAX) = C0
C
      SABS = DNRM2(3,SVEC,1)
C
      IF ( SABS.LT.TOL ) THEN
         DO LM = 1,NLMMAX
            UMAT(LM,LM) = C1
         END DO
         RETURN
      END IF
C
C--------------------------- normalize vector to get spherical harmonics
C
      X = SVEC(1)/SABS
      Y = SVEC(2)/SABS
      Z = SVEC(3)/SABS
C
      CALL CALC_RHPLM(X,Y,Z,RYLM,L3MAX,NLM3MAX)
C
      IF ( DIMAG(ERYD).GE.-1D-15 ) THEN
         ZZ = SQRT(ERYD)*SABS
      ELSE
         ZZ = -DCONJG(SQRT(DCONJG(ERYD)))*SABS
      END IF
C
C
C      ZZ = SQRT(ERYD)*SABS
C
      CALL BESHAN(CHL,CJL,CNL,ZZ,L3MAX)
C
C=======================================================================
C
      IG123 = 0
      LM1LM2 = 0
      DO LM1 = 1,NLM
         L1 = L_LM(LM1)
         DO LM2 = 1,LM1
            L2 = L_LM(LM2)
            LM1LM2 = LM1LM2 + 1
C
            DO I = 1,NRGNT(LM1LM2)
               IG123 = IG123 + 1
               LM3 = IRGNT(IG123)
               L3 = L_LM(LM3)
C
               CLLL = CIPWL(L1+L3-L2)*RGNT(IG123)
C
               UMAT(LM1,LM2) = UMAT(LM1,LM2) + CLLL*CJL(L3)*RYLM(LM3)
C
            END DO
C
            UMAT(LM2,LM1) = UMAT(LM1,LM2)*M1_PW_LL(LM2,LM1)
C
         END DO
      END DO
C
99001 FORMAT (//,1X,79('*'),/,34X,'<INIT_UMAT>',/,1X,79('*'),//,10X,
     &        'initialize the U-transformation matrices for NL =',I4,/)
      END
