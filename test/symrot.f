C*==symrot.f    processed by SPAG 6.70Rc at 08:10 on  8 Mar 2017
      SUBROUTINE SYMROT(IPRINT,NL,MROTR,NSYM,SYMACCEPTED,SYMUNITARY,
     &                  SYMEULANG,IREL,NKM,NK,DROT,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  - set up the rotation matrices  DROT                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1,PI
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMROT')
      LOGICAL CHECK
      PARAMETER (CHECK=.TRUE.)
C
C Dummy arguments
C
      INTEGER IPRINT,IREL,NK,NKM,NKMMAX,NL,NSYM
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX)
      REAL*8 MROTR(3,3,NSYMMAX),SYMEULANG(3,NSYMMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      REAL*8 CTET,DEVMAX,DEVTOL,DIR(3),MINV(3,3),MJ,PHI,RJ,RV(3),RVP(3),
     &       SK,STET,WGAU(:),XGAU(:)
      COMPLEX*16 D1,D2,DNUM(:,:,:),DTIM(:,:),DYLM(:),W2(:,:),YLM(:)
      INTEGER I,I1,I2,IA_ERR,IINV,IRELEFF,ISYM,ISYMP,ITOP,IXY,IZ,J,K,L,
     &        LCHPMAX,LM1,LM2,M,MJM05,NERR,NGAU,NKEFF,NLM,NLMCHPMAX,
     &        NSBLOCK,NSYMH,NYLM
      CHARACTER*2 STR2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DTIM,W2
      ALLOCATABLE XGAU,WGAU,YLM,DYLM,DNUM
C
C-----------------------------------------------------------------------
C               initialize complex spherical harmonics
C-----------------------------------------------------------------------
      LCHPMAX = NL - 1
      NLMCHPMAX = (LCHPMAX+1)**2
C
      ALLOCATE (YLM(NLMCHPMAX),DYLM(NLMCHPMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: Q_CHP')
C
C-----------------------------------------------------------------------
C
      ALLOCATE (DTIM(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DTIM')
C
      WRITE (6,99002) IREL
C
      NSYMH = NSYM/2
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C            calculate the rotation matrices   DROT  for  COMPLEX
C      spherical harmonics analytically and numerically and compare
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IF ( CHECK ) THEN
C
         CALL CINIT(NKMMAX*NKMMAX*NSYM,DROT)
C
         DO ISYM = 1,NSYMH
            CALL ROTMAT(NL,0,SYMEULANG(1,ISYM),SYMEULANG(2,ISYM),
     &                  SYMEULANG(3,ISYM),DROT(1,1,ISYM),NKMMAX)
         END DO
C
         NYLM = NL*NL
         NGAU = 20
         DEVTOL = 1D-8
C
         M = NKMMAX
         ALLOCATE (XGAU(NGAU),WGAU(NGAU))
         ALLOCATE (DNUM(M,M,NSYMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DNUM')
C
         CALL GAULEG(-1D0,+1D0,XGAU,WGAU,NGAU)
C
         DNUM(:,:,:) = C0
C
         DO IZ = 1,NGAU
            CTET = XGAU(IZ)
            STET = SQRT(1D0-CTET*CTET)
            DO IXY = 1,NGAU
               PHI = PI*(1D0+XGAU(IXY))
               RV(1) = STET*COS(PHI)
               RV(2) = STET*SIN(PHI)
               RV(3) = CTET
C
               DO ISYM = 1,NSYMH
                  DO I = 1,3
                     DO J = 1,3
                        MINV(I,J) = MROTR(J,I,ISYM)
                     END DO
                  END DO
                  CALL DGEMV('N',3,3,1D0,MINV,3,RV,1,0D0,RVP,1)
C
                  DIR(1:3) = RV(1:3)
                  CALL RVECNORM(3,DIR)
C
                  CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),YLM,LCHPMAX,
     &                            NLMCHPMAX)
C
                  DIR(1:3) = RVP(1:3)
                  CALL RVECNORM(3,DIR)
C
                  CALL CALC_CHPLM(DIR(1),DIR(2),DIR(3),DYLM,LCHPMAX,
     &                            NLMCHPMAX)
C
                  DO LM2 = 1,NYLM
                     DO LM1 = 1,NYLM
                        DNUM(LM1,LM2,ISYM) = DNUM(LM1,LM2,ISYM)
     &                     + WGAU(IZ)*WGAU(IXY)*DCONJG(YLM(LM1))
     &                     *DYLM(LM2)*PI
                     END DO
                  END DO
               END DO
C
            END DO
         END DO
C
         NERR = 0
         DEVMAX = 0D0
         DO ISYM = 1,NSYMH
            DO LM1 = 1,NYLM
               DO LM2 = 1,NYLM
                  D1 = DROT(LM1,LM2,ISYM)
                  D2 = DNUM(LM1,LM2,ISYM)
                  DEVMAX = MAX(DEVMAX,ABS(D1-D2))
                  IF ( ABS(D1-D2).GT.DEVTOL .AND. IPRINT.GT.2 ) THEN
                     WRITE (6,99001) ISYM,LM1,LM2,D1,D2
                     NERR = NERR + 1
                  END IF
               END DO
            END DO
         END DO
C
         WRITE (6,99003) NSYMH*NYLM*NYLM,NERR,DEVTOL,DEVMAX
C
         DEALLOCATE (YLM,DYLM,XGAU,WGAU,DNUM)
C
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C-----------------------------------------------------------------------
C                    initialize all rotation matrices
C-----------------------------------------------------------------------
C
      CALL CINIT(NKMMAX*NKMMAX*NSYM,DROT)
C
C-----------------------------------------------------------------------
C                       create rotation matrices
C-----------------------------------------------------------------------
C
      IF ( IREL.LE.2 ) THEN
         IRELEFF = 0
         NKEFF = NL
      ELSE
         IRELEFF = 3
         NKEFF = NK
      END IF
C
      DO ISYM = 1,NSYMH
C
         CALL ROTMAT(NKEFF,IRELEFF,SYMEULANG(1,ISYM),SYMEULANG(2,ISYM),
     &               SYMEULANG(3,ISYM),DROT(1,1,ISYM),NKMMAX)
C
      END DO
C-----------------------------------------------------------------------
C                     create matrix for inversion
C-----------------------------------------------------------------------
      IINV = NSYMH + 1
      CALL CINIT(NKMMAX*NKMMAX,DROT(1,1,IINV))
C
      I = 0
      IF ( IREL.GT.2 ) THEN
         NSBLOCK = 2
      ELSE
         NSBLOCK = 1
      END IF
      DO L = 0,(NL-1)
         DO M = 1,NSBLOCK*(2*L+1)
            I = I + 1
            DROT(I,I,IINV) = (-1.0D0)**L
         END DO
      END DO
      ITOP = I
C
      CALL DGEMM('N','N',3,3,3,-1D0,MROTR(1,1,1),3,MROTR(1,1,1),3,0D0,
     &           MROTR(1,1,IINV),3)
C
C-----------------------------------------------------------------------
C                         include inversion
C-----------------------------------------------------------------------
      DO ISYM = 2,NSYMH
         ISYMP = NSYMH + ISYM
C
         CALL ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM),NKMMAX,
     &              DROT(1,1,IINV),NKMMAX,C0,DROT(1,1,ISYMP),NKMMAX)
C
         CALL DGEMM('N','N',3,3,3,1D0,MROTR(1,1,ISYM),3,MROTR(1,1,IINV),
     &              3,0D0,MROTR(1,1,ISYMP),3)
C
      END DO
C
C-----------------------------------------------------------------------
C            add second spin-diagonal block for  IREL=2
C            spin off-diagonal blocks have been initialized before
C-----------------------------------------------------------------------
      IF ( IREL.EQ.2 ) THEN
C
         NLM = NKM/2
         IF ( ITOP.NE.NLM )
     &         CALL STOP_MESSAGE(ROUTINE,'ITOP <> NLM for IREL = 2')
C
         DO ISYM = 1,NSYM
            DO J = 1,NLM
               CALL ZCOPY(NLM,DROT(1,J,ISYM),1,DROT(NLM+1,NLM+J,ISYM),1)
            END DO
         END DO
C
      END IF
C-----------------------------------------------------------------------
C                     create matrix for time reversal
C-----------------------------------------------------------------------
      IF ( IREL.GT.2 ) THEN
C
         CALL CINIT(NKMMAX*NKMMAX,DTIM)
C
         I = 0
         DO K = 1,NK
            L = K/2
            IF ( L*2.EQ.K ) THEN
               SK = -1D0
            ELSE
               SK = +1D0
            END IF
            RJ = L + SK*0.5D0
            DO MJM05 = NINT(-RJ-0.5D0),NINT(RJ-0.5D0)
               MJ = DBLE(MJM05) + 0.5D0
               I1 = NINT(2*L*(RJ+0.5D0)+RJ+MJ+1)
               I2 = NINT(2*L*(RJ+0.5D0)+RJ-MJ+1)
C               DTIM(I1,I2) = CI*SK*(-1)**NINT(MJ+0.5D0)
               DTIM(I1,I2) = SK*(-1)**NINT(MJ+0.5D0)
            END DO
         END DO
C
      END IF
C=======================================================================
C            set up of transformation matrices completed
C=======================================================================
C
C=======================================================================
C   include time reversal operation for anti-unitary transformations
C   this can occur only for IREL > 2 i.e. fully relativistic case
C=======================================================================
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
            IF ( .NOT.SYMUNITARY(ISYM) ) THEN
C
               IF ( IREL.LE.2 ) CALL STOP_MESSAGE(ROUTINE,
     &             'anti-unitary symmetry matrices created for IREL = 2'
     &             )
C
               CALL ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM),NKMMAX,
     &                    DTIM,NKMMAX,C0,W2,NKMMAX)
C
               DO J = 1,NKM
                  CALL ZCOPY(NKM,W2(1,J),1,DROT(1,J,ISYM),1)
               END DO
C
            END IF
         END IF
      END DO
C
C=======================================================================
C                 print transformation matrices on request
C=======================================================================
C
      IF ( IPRINT.GT.0 ) THEN
C
         DO ISYM = 1,NSYM
            IF ( SYMACCEPTED(ISYM) ) THEN
               WRITE (STR2,FMT='(I2)') ISYM
C
               CALL CMATSTRUCT('symmetry operation '//STR2,
     &                         DROT(1,1,ISYM),NKM,NKMMAX,IREL,IREL,0,
     &                         1D-8,6)
            END IF
         END DO
C
         CALL CMATSTRUCT('Inversion     MATRIX',DROT(1,1,IINV),NKM,
     &                   NKMMAX,IREL,IREL,0,1D-8,6)
         IF ( IREL.GT.2 ) CALL CMATSTRUCT('Time reversal MATRIX',DTIM,
     &        NKM,NKMMAX,3,3,0,1D-8,6)
      END IF
C
C=======================================================================
C
      DEALLOCATE (DTIM,W2,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'DEALLOC'
C
99001 FORMAT (2X,'TROUBLE with  ROT-matrix  D  for SYM',I3,3X,2I3,3X,
     &        2E12.5,/,53X,2E12.5,/)
99002 FORMAT (//,1X,79('*'),/,36X,'<SYMROT>',/,1X,79('*'),//,10X,
     &        'set up of rotation matrices  DROT  for IREL =',I2)
99003 FORMAT (/,10X,'CHECK of rotation matrices   DROT  for COMPLEX ',
     &        'spherical harmonics',/,10X,'comparison of ',
     &        'analytical and numerical results',/,10X,'out of',I6,
     &        ' elements',I6,' deviate more than TOL =',1PE10.1,/,10X,
     &        'maximum deviation found: ',1PE10.1,/)
      END
C
