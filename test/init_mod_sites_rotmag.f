C*==init_mod_sites_rotmag.f    processed by SPAG 6.70Rc at 09:36 on 30 Mar 2017
      SUBROUTINE INIT_MOD_SITES_ROTMAG
C   ********************************************************************
C   *                                                                  *
C   *  calculate rotation matrices  MROTQ  and  DROTQ                  *
C   *  connecting LOCAL and GLOBAL frame individually for              *
C   *  all sites IQ according to the flag  KMROT                       *
C   *                                                                  *
C   *  KMROT                                                           *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site IQ   *
C   *  2: global COMMON rotation of the magnetisation                  *
C   *  3: spin spiral    Theta =  90                                   *
C   *  4: spin spiral    Theta <> 90                                   *
C   *                                                                  *
C   *  NOTE: theta and phi are the standard polar angles               *
C   *        of the magnetic moment vector w.r.t the GLOBAL frame      *
C   *        the angle gamma allows to perform an additional           *
C   *        reorientation of the LOCAL frame by rotating around the   *
C   *        magnetic moment vector to make full use of symmetry       *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:FOUND_SECTION
      USE MOD_SITES,ONLY:QMPHI,QMTET,QMGAM,MROTQ,DROTQ,NQ,NQMAX
      USE MOD_ANGMOM,ONLY:NK,NLM,NKMMAX,WKM1,WKM2
      USE MOD_CONSTANTS,ONLY:C1,C0,CI,DEG_ARC
      USE MOD_CALCMODE,ONLY:KMROT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_SITES_ROTMAG')
C
C Local variables
C
      LOGICAL CHANGE_FRAME
      REAL*8 CTET05,MCHK(3,3),STET05,XPHI,XTET
      COMPLEX*16 EMPHI05,EPPHI05
      INTEGER I,IA_ERR,IQ,N
      LOGICAL RVEC_SAME
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( ALLOCATED(MROTQ) ) DEALLOCATE (MROTQ,DROTQ)
C
C=======================================================================
C  Finds whether the representation is supposed to be changed
C     (no matter whether magnetization is rotated or not)
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( .NOT.FOUND_SECTION ) THEN
C
         CHANGE_FRAME = .FALSE.
C
      ELSE
C
         CALL SECTION_FIND_KEYWORD('FRAMETET',CHANGE_FRAME)
C
         IF ( .NOT.CHANGE_FRAME ) THEN
C
            CALL SECTION_FIND_KEYWORD('FRAMEPHI',CHANGE_FRAME)
C
            IF ( .NOT.CHANGE_FRAME )
     &            CALL SECTION_FIND_KEYWORD('FRAMEGAM',CHANGE_FRAME)
C
         END IF
C
      END IF
C
C=======================================================================
C
      IF ( KMROT.LE.0 .AND. .NOT.CHANGE_FRAME ) THEN
         ALLOCATE (MROTQ(1,1,1),DROTQ(1,1,1))
         MROTQ = 999999D0
         DROTQ = 999999D0
         RETURN
      END IF
C
C=======================================================================
C
      ALLOCATE (MROTQ(3,3,NQMAX))
      ALLOCATE (DROTQ(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROTQ')
C
C-----------------------------------------------------------------------
C                  non-collinear spin structure
C-----------------------------------------------------------------------
C
      IF ( KMROT.LE.2 ) THEN
C
         DO IQ = 1,NQ
C
            CALL ROTMAT(NK,3,QMPHI(IQ),QMTET(IQ),0.0D0,DROTQ(1,1,IQ),
     &                  NKMMAX)
C
            CALL GETMROT(-QMGAM(IQ),-QMTET(IQ),-QMPHI(IQ),MROTQ(1,1,IQ))
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C               check NEW set up of rotation matrices
C               this section can later on be removed
C               use extended version based on GETMROT instead
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            XPHI = QMPHI(IQ)*DEG_ARC
            XTET = QMTET(IQ)*DEG_ARC
C
            MCHK(1,1) = COS(XPHI)*COS(XTET)
            MCHK(1,2) = SIN(XPHI)*COS(XTET)
            MCHK(1,3) = -SIN(XTET)
            MCHK(2,1) = -SIN(XPHI)
            MCHK(2,2) = COS(XPHI)
            MCHK(2,3) = 0D0
            MCHK(3,1) = COS(XPHI)*SIN(XTET)
            MCHK(3,2) = SIN(XPHI)*SIN(XTET)
            MCHK(3,3) = COS(XTET)
C
            IF ( .NOT.RVEC_SAME(9,MCHK,MROTQ(1,1,IQ),1D-8) ) THEN
               WRITE (6,*) '  for site IQ = ',IQ
               WRITE (6,'(''MCHK = '',/,(3F12.8))') MCHK
               WRITE (6,'(''MROT = '',/,(3F12.8))') MROTQ(:,:,IQ)
               CALL STOP_MESSAGE(ROUTINE,'MROT matrices DIFFER')
            END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
         END DO
C
      ELSE
C
C-----------------------------------------------------------------------
C                    spin spiral - implies IREL=2
C                    set up the rotation matrices
C-----------------------------------------------------------------------
C
         N = 2*NLM
C
         DO IQ = 1,NQ
C
            EPPHI05 = EXP(CI*0.5D0*QMPHI(IQ)*DEG_ARC)
            EMPHI05 = C1/EPPHI05
            CTET05 = COS(0.5D0*QMTET(IQ)*DEG_ARC)
            STET05 = SIN(0.5D0*QMTET(IQ)*DEG_ARC)
C
            WKM1(1:N,1:N) = C0
            WKM2(1:N,1:N) = C0
C
            DO I = 1,NLM
               WKM1(I,I) = CTET05
               WKM1(NLM+I,NLM+I) = CTET05
               WKM1(I,NLM+I) = STET05
               WKM1(NLM+I,I) = -STET05
               WKM2(I,I) = EPPHI05
               WKM2(NLM+I,NLM+I) = EMPHI05
            END DO
C
            DROTQ(1:N,1:N,IQ) = MATMUL(WKM1(1:N,1:N),WKM2(1:N,1:N))
C
         END DO
C
         MROTQ(1:3,1:3,1:NQ) = 0.0D0
C
      END IF
C
      END
C*==init_mod_types_rotmag.f    processed by SPAG 6.70Rc at 09:36 on 30 Mar 2017
      SUBROUTINE INIT_MOD_TYPES_ROTMAG
C   ********************************************************************
C   *                                                                  *
C   *  initialize  MROT_T  and  DROT_T according to  MROTQ  and  DROTQ *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT,MROTQ,DROTQ
      USE MOD_TYPES,ONLY:MROT_T,DROT_T,ITBOT,ITTOP,NTMAX
      USE MOD_ANGMOM,ONLY:NKMMAX,NK
      USE MOD_CALCMODE,ONLY:MOMENTS_ROTATED
      USE MOD_SCF,ONLY:SCF_THETA_DEPENDENT_POT
      USE MOD_THERMAL,ONLY:NTET_FLUCT_POT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_TYPES_ROTMAG')
C
C Local variables
C
      COMPLEX*16 DROT_TET(:,:,:)
      REAL*8 DTET,MGAM,MPHI,MROT_TET(:,:,:),MTET
      INTEGER IA_ERR,IQ,IT,ITET,IT_CHEM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MROT_TET,DROT_TET
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( .NOT.MOMENTS_ROTATED ) RETURN
C
C=======================================================================
C
      IF ( .NOT.ALLOCATED(MROT_T) ) THEN
         ALLOCATE (MROT_T(3,3,NTMAX))
         ALLOCATE (DROT_T(NKMMAX,NKMMAX,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROT_T')
      END IF
C
C-----------------------------------------------------------------------
C     theta-dependent potential for type IT - i.e. pseudo types used
C-----------------------------------------------------------------------
C
      IF ( SCF_THETA_DEPENDENT_POT ) THEN
C
         ALLOCATE (MROT_TET(3,3,NTET_FLUCT_POT))
         ALLOCATE (DROT_TET(NKMMAX,NKMMAX,NTET_FLUCT_POT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DROT_TET')
C
         DTET = 180.D0/(NTET_FLUCT_POT-1)
         MPHI = 0D0
         MGAM = 0D0
         DO ITET = 1,NTET_FLUCT_POT
            MTET = (ITET-1)*DTET
C
            CALL ROTMAT(NK,3,MPHI,MTET,MGAM,DROT_TET(1,1,ITET),NKMMAX)
C
            CALL GETMROT(-MGAM,-MTET,-MPHI,MROT_TET(1,1,ITET))
C
         END DO
C
C------------------------------------------------ loop over pseudo types
         DO IT = ITBOT,ITTOP
C
            IT_CHEM = (IT-1)/NTET_FLUCT_POT + 1
            ITET = IT - (IT_CHEM-1)*NTET_FLUCT_POT
C
            MROT_T(:,:,IT) = MROT_TET(:,:,ITET)
            DROT_T(:,:,IT) = DROT_TET(:,:,ITET)
C
         END DO
C
      ELSE
C-----------------------------------------------------------------------
C  common average orientation of moments for all types IT on site IQ
C-----------------------------------------------------------------------
C
         DO IT = ITBOT,ITTOP
C
            IQ = IQAT(1,IT)
C
            MROT_T(:,:,IT) = MROTQ(:,:,IQ)
            DROT_T(:,:,IT) = DROTQ(:,:,IQ)
C
         END DO
C
C-----------------------------------------------------------------------
      END IF
C
      END
