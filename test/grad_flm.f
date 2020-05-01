C*==grad_flm.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GRAD_FLM(FLM,KFLM,GFLM,KGFLM,NLMFP,IM)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the gradient of    SUM_LM  FLM(IR,LM) * Y(LM)         *
C   *                                                                  *
C   *  grad SUM_LM FLM(IR,LM) * Y(LM)                                  *
C   *                                                                  *
C   *     = SUM_LM' GFLM(IR,LM',IC) * Y(LM') * ->E(IC)                 *
C   *                                                                  *
C   *                   with IC = 1,2,3 == x,y,z                       *
C   *                   and REAL spherical harmonics Y(LM)             *
C   *                                                                  *
C   *  making use of gradient formula for COMPLEX spherical harmonics  *
C   *                                                                  *
C   *  KFLM(LM) indicates non-0 terms  FLM(IR,LM)                      *
C   *                                                                  *
C   *  NOTE: the conversion factors are initialized in first call      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NLFP
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:L_LM
      USE MOD_RMESH,ONLY:NRMAX,R,NPAN,JRCUT,DRDI,JRCRI,FULLPOT,JRWS
      USE MOD_CONSTANTS,ONLY:C0,C1,CI,SQRT_2
      IMPLICIT NONE
C*--GRAD_FLM28
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GRAD_FLM')
C
C Dummy arguments
C
      INTEGER IM,NLMFP
      REAL*8 FLM(NRMAX,NLMFP),GFLM(NRMAX,NLMFP,3)
      INTEGER KFLM(NLMFP),KGFLM(NLMFP,3)
C
C Local variables
C
      REAL*8 CGC_RACAH
      REAL*8 DFDR(:),FAC_LMIN1,FAC_LPLS1,FOVR,FRMIN(:),FRPLS(:),
     &       SGMIN(:,:,:),SGPLS(:,:,:),SMIN,SPLS
      COMPLEX*16 F,FNAB,FNRC,RC_FP(:,:),SNAB(-1:+1,3),TGMIN(:,:),
     &           TGPLS(:,:)
      INTEGER IC,IPAN,IR,IR1,IRTOP,L,L2,L3,LM,LM1,LM2,LM3,M,M1,M1LIM,
     &        M1STEP,M2,M3,M3LIM,M3STEP,MM,MMLIM,NR
      LOGICAL INITIALIZE
      SAVE SGMIN,SGPLS
      EXTERNAL DAXPY
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RC_FP,SGPLS,SGMIN,TGPLS,TGMIN
      ALLOCATABLE FRPLS,FRMIN,DFDR
C
      DATA INITIALIZE/.TRUE./
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
C ----------------------------------------------------------------------
C     allow to deal with full potential charge         NLFP = 2 * NL - 1
C ----------------------------------------------------------------------
C
         IF ( NLMFP.NE.NLFP**2 )
     &         CALL STOP_MESSAGE(ROUTINE,' NLMFP<> NLFP**2')
C
C-----------------------------------------------------------------------
C   transformation matrix  RC_FP = RC:  complex -> real spher. harm.
C-----------------------------------------------------------------------
C
         ALLOCATE (RC_FP(NLMFP,NLMFP))
C
         CALL CALC_RC(NLFP-1,RC_FP,NLMFP,1)
C
C-----------------------------------------------------------------------
C            transformation matrices  SGPLS and SGMIN
C-----------------------------------------------------------------------
C
         SNAB(:,:) = C0
         SNAB(-1,1) = +C1/SQRT_2
         SNAB(+1,1) = -C1/SQRT_2
         SNAB(-1,2) = +CI/SQRT_2
         SNAB(+1,2) = +CI/SQRT_2
         SNAB(0,3) = C1
C
         ALLOCATE (SGPLS(NLMFP,NLMFP,3),SGMIN(NLMFP,NLMFP,3))
         ALLOCATE (TGPLS(NLMFP,NLMFP),TGMIN(NLMFP,NLMFP))
         SGPLS(:,:,:) = 0D0
         SGMIN(:,:,:) = 0D0
C
         LM = 0
         DO L = 0,NLFP - 1
            FAC_LPLS1 = SQRT(DBLE(L+1)/DBLE(2*L+3))
            FAC_LMIN1 = SQRT(DBLE(L)/DBLE(2*L-1))
C
            DO M = -L, + L
               LM = LM + 1
C
               DO IC = 1,3
C
                  TGPLS(:,:) = 0D0
                  TGMIN(:,:) = 0D0
C
                  IF ( IC.LE.2 ) THEN
                     MMLIM = 1
                  ELSE
                     MMLIM = 0
                  END IF
C
                  DO MM = -MMLIM, + MMLIM,2
C
                     FNAB = SNAB(MM,IC)
C
                     M1LIM = ABS(M)
                     M1STEP = MAX(1,2*M1LIM)
                     DO M1 = -M1LIM, + M1LIM,M1STEP
                        LM1 = L*(L+1) + M1 + 1
                        FNRC = FNAB*DCONJG(RC_FP(LM,LM1))
C
C-------------------------------------------------------------- L2 = L+1
                        L2 = L + 1
                        M2 = M1 + MM
                        LM2 = L2*(L2+1) + M2 + 1
C
                        IF ( L2.LE.NLFP-1 .AND. ABS(M2).LE.L2 ) THEN
C
                           F = FNRC*FAC_LPLS1*CGC_RACAH(DBLE(L),1D0,
     &                         DBLE(L2),DBLE(M1),DBLE(MM),DBLE(M2))
C
                           L3 = L2
                           M3LIM = ABS(M2)
                           M3STEP = MAX(1,2*M3LIM)
                           DO M3 = -M3LIM, + M3LIM,M3STEP
                              LM3 = L3*(L3+1) + M3 + 1
                              TGPLS(LM,LM3) = TGPLS(LM,LM3)
     &                           + F*RC_FP(LM3,LM2)
                           END DO
                        END IF
C
C-------------------------------------------------------------- L2 = L-1
                        L2 = L - 1
                        M2 = M1 + MM
                        LM2 = L2*(L2+1) + M2 + 1
C
                        IF ( L2.GE.0 .AND. ABS(M2).LE.L2 ) THEN
C
                           F = FNRC*FAC_LMIN1*CGC_RACAH(DBLE(L),1D0,
     &                         DBLE(L2),DBLE(M1),DBLE(MM),DBLE(M2))
C
                           L3 = L2
                           M3LIM = ABS(M2)
                           M3STEP = MAX(1,2*M3LIM)
                           DO M3 = -M3LIM, + M3LIM,M3STEP
                              LM3 = L3*(L3+1) + M3 + 1
                              TGMIN(LM,LM3) = TGMIN(LM,LM3)
     &                           + F*RC_FP(LM3,LM2)
                           END DO
                        END IF
C-----------------------------------------------------------------------
C
                     END DO
C
                  END DO
C
C-----------------------------------------------------------------------
                  DO LM3 = 1,NLMFP
C
                     SGMIN(LM,LM3,IC) = DREAL(TGMIN(LM,LM3))
                     SGPLS(LM,LM3,IC) = DREAL(TGPLS(LM,LM3))
C
                     IF ( ABS(DIMAG(TGMIN(LM,LM3))).GT.1D-10 )
     &                    WRITE (6,*) 'TGMIN ',LM,LM3,TGMIN(LM,LM3)
                     IF ( ABS(DIMAG(TGPLS(LM,LM3))).GT.1D-10 )
     &                    WRITE (6,*) 'TGPLS ',LM,LM3,TGPLS(LM,LM3)
                  END DO
C-----------------------------------------------------------------------
C
               END DO
C
            END DO
         END DO
C
         DEALLOCATE (RC_FP)
C
         IF ( IPRINT.GE.1 ) THEN
            M = NLMFP
            CALL RMATSTRUCT('SGMIN1',SGMIN(1,1,1),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGMIN2',SGMIN(1,1,2),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGMIN3',SGMIN(1,1,3),M,M,1,1,0,1D-8,6)
C
            CALL RMATSTRUCT('SGPLS1',SGPLS(1,1,1),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGPLS2',SGPLS(1,1,2),M,M,1,1,0,1D-8,6)
            CALL RMATSTRUCT('SGPLS3',SGPLS(1,1,3),M,M,1,1,0,1D-8,6)
         END IF
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      ALLOCATE (FRPLS(NRMAX),FRMIN(NRMAX),DFDR(NRMAX))
C
C---------------------------------------- ensure proper settings for ASA
      IF ( .NOT.FULLPOT ) THEN
         NPAN(IM) = 1
         JRCUT(0,IM) = 0
         JRCUT(1,IM) = JRWS(IM)
         JRCRI(IM) = JRWS(IM)
      END IF
C
      IRTOP = JRCRI(IM)
      GFLM(:,:,:) = 0D0
      KGFLM(:,:) = 0
C
      DO LM = 1,NLMFP
         IF ( KFLM(LM).NE.0 ) THEN
C
            DO IPAN = 1,NPAN(IM)
               IR1 = JRCUT(IPAN-1,IM) + 1
               NR = JRCUT(IPAN,IM) - JRCUT(IPAN-1,IM)
C
               CALL CALCDFDR(FLM(IR1,LM),DFDR(IR1),DRDI(IR1,IM),NR)
C
            END DO
C
            L = L_LM(LM)
C
            DO IR = 1,IRTOP
               FOVR = FLM(IR,LM)/R(IR,IM)
               FRMIN(IR) = -(DFDR(IR)+(L+1)*FOVR)
               FRPLS(IR) = DFDR(IR) - L*FOVR
            END DO
C
C-----------------------------------------------------------------------
            DO IC = 1,3
               DO LM3 = 1,NLMFP
C
                  SMIN = SGMIN(LM,LM3,IC)
                  IF ( ABS(SMIN).GT.1D-8 ) THEN
                     KGFLM(LM3,IC) = 1
                     CALL DAXPY(IRTOP,SMIN,FRMIN,1,GFLM(1,LM3,IC),1)
                  END IF
C
                  SPLS = SGPLS(LM,LM3,IC)
                  IF ( ABS(SPLS).GT.1D-8 ) THEN
                     KGFLM(LM3,IC) = 1
                     CALL DAXPY(IRTOP,SPLS,FRPLS,1,GFLM(1,LM3,IC),1)
                  END IF
C
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END IF
      END DO
C
      END
