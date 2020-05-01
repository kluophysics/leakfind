C*==scfbimad.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBIMAD
C   ********************************************************************
C   *                                                                  *
C   *   calculate the Madelung matrix connected with the               *
C   *   Breit Interaction                                              *
C   *                                                                  *
C   *   - it is assumed that all moments are collinearly               *
C   *     aligned a common local axis specified by the standard        *
C   *     Euler angles (ALFDEG,BETDEG,GAMDEG)                          *
C   *   - the Madelung sum is expressed with respect to the local      *
C   *     frame of reference, as a consequence the lattice vectors     *
C   *     ABAS, and QBAS have to be transferred to the local frame     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI,C0,CI
      USE MOD_SITES,ONLY:ABIMAD,NQ,NQMAX,QBAS,ALFDEG,BETDEG,GAMDEG,
     &    IQBOT,IQTOP
      USE MOD_ANGMOM,ONLY:NL,NLABIMAX
      USE MOD_LATTICE,ONLY:ALAT,ABAS,SYSTEM_DIMENSION
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--SCFBIMAD24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBIMAD')
      LOGICAL CHECK
      PARAMETER (CHECK=.TRUE.)
C
C Local variables
C
      REAL*8 ABAS_LOC(:,:),DROT_RLM(:,:),M1PWABSM2,M1PWL2,MROT(3,3),
     &       PREFAC,PRE_GNT,QBAS_LOC(:,:),RSUM,SMAT(:,:,:),
     &       SMAT_TMP(:,:,:),V(3),VP(3),WZ05
      REAL*8 GAUNT_CYLM
      INTEGER I,IA_ERR,IQ,JQ,L,L1,L2,L2SQ,LM1,LM2,LMRM,LMRP,M,M1,M2,MRM,
     &        MRP,NLAMAX,NLM2MAX,NLMAMAX,NLM_ROT,NLX,NL_ROT,NQ_LOC
      COMPLEX*16 SMATC,WM,WP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SMAT,ABAS_LOC,QBAS_LOC,DROT_RLM,SMAT_TMP
C
      CALL TRACK_INFO(ROUTINE)
C
      WRITE (6,99001) ALFDEG,BETDEG,GAMDEG
C
C-----------------------------------------------------------------------
C
      NLAMAX = (2*(NL-1)+1) + 1
      NLMAMAX = NLAMAX*NLAMAX
      NQ_LOC = NQ
C
C-----------------------------------------------------------------------
C       calculate lattice sum    SMAT = sum(R) 1/R^(l+2) Y_L(R)
C                      in the GLOBAL frame
C-----------------------------------------------------------------------
C
      NLX = NL + 1
      NLM2MAX = (2*(NLX-1)+1)**2
C
      ALLOCATE (SMAT(NQ_LOC,NQ_LOC,NLM2MAX))
      ALLOCATE (SMAT_TMP(NQ_LOC,NQ_LOC,NLM2MAX))
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C
         CALL SCFMAD3D_SMAT(ALAT,ABAS,QBAS,NQ_LOC,IPRINT,NLX,SMAT_TMP,
     &                      NLAMAX,NLM2MAX)
C
      ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
         CALL SCFBI2DSTRMAT(SMAT_TMP,NLX,NLM2MAX)
C
      END IF
C
C-----------------------------------------------------------------------
C       rotate     SMAT    to the LOCAL frame
C-----------------------------------------------------------------------
C
      NL_ROT = 2*(NLX-1) + 1
      NLM_ROT = NL_ROT**2
      IF ( NLM_ROT.LT.NLM2MAX ) STOP 'NLM_ROT < NLM2MAX'
      ALLOCATE (DROT_RLM(NLM_ROT,NLM_ROT))
C
      CALL ROTMAT_RYLM(NL_ROT,NLM_ROT,-GAMDEG,-BETDEG,-ALFDEG,DROT_RLM)
C
      DO IQ = IQBOT,IQTOP
         DO JQ = IQBOT,IQTOP
            DO LM1 = 1,NLM2MAX
               RSUM = 0D0
               DO LM2 = 1,NLM2MAX
                  RSUM = RSUM + SMAT_TMP(IQ,JQ,LM2)*DROT_RLM(LM2,LM1)
               END DO
               SMAT(IQ,JQ,LM1) = RSUM
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C     CHECK: rotate basis and primitive vectors to LOCAL frame
C-----------------------------------------------------------------------
C
      IF ( CHECK .AND. SYSTEM_DIMENSION(1:2).EQ.'3D' ) THEN
C
         CALL GETMROT(ALFDEG,BETDEG,GAMDEG,MROT)
C
         ALLOCATE (ABAS_LOC(3,3),QBAS_LOC(3,NQ_LOC))
C
         WRITE (6,99002)
         DO I = 1,3
            V(1:3) = ABAS(1:3,I)
            CALL DGEMV('N',3,3,1D0,MROT,3,V,1,0D0,VP,1)
            ABAS_LOC(1:3,I) = VP(1:3)
            WRITE (6,99004) I,V,VP
         END DO
C
         WRITE (6,99003)
         DO IQ = 1,NQ_LOC
            V(1:3) = QBAS(1:3,IQ)
            CALL DGEMV('N',3,3,1D0,MROT,3,V,1,0D0,VP,1)
            QBAS_LOC(1:3,IQ) = VP(1:3)
            WRITE (6,99004) IQ,V,VP
         END DO
C
CHECK --------------- calculate SMAT_TMP for rotated vectors and compare
C
         CALL SCFMAD3D_SMAT(ALAT,ABAS_LOC,QBAS_LOC,NQ_LOC,IPRINT,NLX,
     &                      SMAT_TMP,NLAMAX,NLM2MAX)
C
         DO IQ = IQBOT,IQTOP
            DO JQ = IQBOT,IQTOP
               DO LM1 = 1,NLM2MAX
                  IF ( ABS(SMAT_TMP(IQ,JQ,LM1)-SMAT(IQ,JQ,LM1))
     &                 .GT.1D-6 ) WRITE (*,*) 'IQ JQ ',IQ,JQ,' LM1 ',
     &                 LM1,'  SMAT:',SMAT_TMP(IQ,JQ,LM1),SMAT(IQ,JQ,LM1)
               END DO
            END DO
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C
      PREFAC = -(4D0*PI/SQRT(2D0))*SQRT(8D0*PI/3D0)
C
      WZ05 = SQRT(0.5D0)
C
C-----------------------------------------------------------------------
C                 initialize Madelung constants
C-----------------------------------------------------------------------
C
      IF ( .NOT.ALLOCATED(ABIMAD) )
     &     ALLOCATE (ABIMAD(NQMAX,NQMAX,-1:1,NLMAMAX))
C
      ABIMAD(1:NQMAX,1:NQMAX,-1:1,1:NLMAMAX) = C0
C
C-----------------------------------------------------------------------
C                 calculate Madelung constants
C-----------------------------------------------------------------------
C
      L = 1
      DO M = -1, + 1,2
C
         LM1 = 0
         DO L1 = 0,(2*NL-1)
            L2 = L1 + L
            L2SQ = L2*L2
            M1PWL2 = (-1D0)**L2
C
            DO M1 = -L1,L1
               LM1 = LM1 + 1
               M2 = M - M1
C
               PRE_GNT = M*PREFAC*GAUNT_CYLM(L,M,L1,M1,L2,M2)
C
C--------------------------- convert from real Y_L(-R) to complex Y_L(R)
C
               MRM = -ABS(M2)
               MRP = -MRM
               LMRM = L2SQ + L2 + MRM + 1
               LMRP = L2SQ + L2 + MRP + 1
C
               IF ( M2.LT.0 ) THEN
                  WM = -WZ05*CI
                  WP = +WZ05
               ELSE IF ( M2.EQ.0 ) THEN
                  WM = 0D0
                  WP = 1D0
               ELSE
                  M1PWABSM2 = (-1D0)**ABS(M2)
                  WM = M1PWABSM2*WZ05*CI
                  WP = M1PWABSM2*WZ05
               END IF
C
               DO JQ = IQBOT,IQTOP
                  DO IQ = IQBOT,IQTOP
C
                     SMATC = M1PWL2*(WM*SMAT(IQ,JQ,LMRM)
     &                       +WP*SMAT(IQ,JQ,LMRP))
C
                     ABIMAD(IQ,JQ,M,LM1) = PRE_GNT*SMATC
C
                  END DO
               END DO
C
            END DO
         END DO
C
      END DO
C
      IF ( IPRINT.GE.0 ) THEN
         WRITE (6,*) ' Madelung-matrix (A(11) term)'
C
         DO IQ = IQBOT,IQTOP
            WRITE (6,99005) (ABIMAD(IQ,JQ,1,1),JQ=1,NQ)
         END DO
         DO M = -1,1,2
            DO LM1 = 1,NLABIMAX
               WRITE (6,*) ' Madelung-matrix (A term)',M,LM1
               DO IQ = IQBOT,IQTOP
                  WRITE (6,99005) (ABIMAD(IQ,JQ,M,LM1),JQ=IQBOT,IQTOP)
               END DO
            END DO
         END DO
         WRITE (6,*) ' '
      END IF
C
      DEALLOCATE (SMAT,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFBIMAD'
C
99001 FORMAT (//,1X,79('*'),/,34X,'<SCFBIMAD>',/,1X,79('*'),//,10X,
     &        'rotate the spatial lattice vectors ',
     &        'according to the  Euler-angles: ',3F6.1)
99002 FORMAT (/,10X,'primitive lattice vectors  ->A_i ',/)
99003 FORMAT (/,10X,'basis vectors  Q',/)
99004 FORMAT (10X,I3,3X,2('(',F5.2,',',F5.2,',',F5.2,')',:,'  ->  '))
99005 FORMAT ((5X,7F10.6))
C
      END
