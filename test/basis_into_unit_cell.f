C*==basis_into_unit_cell.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE BASIS_INTO_UNIT_CELL(DIMENSION,ABAS,ADAINV,ABAS_I,
     &                                ADAINV_I,QBAS,QBAS_UC,NQ,
     &                                OVERWRITE)
C   ********************************************************************
C   *                                                                  *
C   *   shift NQ basis vectors into the unit cell                      *
C   *                                                                  *
C   *                     QBAS -> QBAS_UC                              *
C   *                                                                  *
C   *   using a set prinitive basis vectors  ABAS_DIM                  *
C   *   depending on the depending on the dimension of the system      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--BASIS_INTO_UNIT_CELL17
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-6)
C
C Dummy arguments
C
      CHARACTER*10 DIMENSION
      INTEGER NQ
      LOGICAL OVERWRITE
      REAL*8 ABAS(3,3),ABAS_I(3,3),ADAINV(3,3),ADAINV_I(3,3),QBAS(3,NQ),
     &       QBAS_UC(3,NQ)
C
C Local variables
C
      REAL*8 ABAS_DIM(3,3),ADAINV_DIM(3,3),CF(3),QBAS_TMP(:,:),
     &       QVEC_SHIFT(3)
      INTEGER I,IQ,NI
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE QBAS_TMP
C
      ALLOCATE (QBAS_TMP(3,NQ))
C
C--------------------------- use a set prinitive basis vectors  ABAS_DIM
C              depending on the depending on the dimension of the system
C
      IF ( DIMENSION(1:2).EQ.'3D' ) THEN
C
         ABAS_DIM(:,:) = ABAS(:,:)
         ADAINV_DIM(:,:) = ADAINV(:,:)
C
      ELSE
C
         ABAS_DIM(:,:) = ABAS_I(:,:)
         ADAINV_DIM(:,:) = ADAINV_I(:,:)
C
      END IF
C
C=======================================================================
      DO IQ = 1,NQ
C
         CALL RVECEXPAND(QBAS(1,IQ),ABAS_DIM,ADAINV_DIM,CF)
C
C-------- shift ->q to have ->q = sum_i ->a(i) + c(i) with 0 <= c(i) < 1
C
         DO I = 1,3
            IF ( ABS(CF(I)).GT.(0D0-TOL) ) THEN
               NI = INT(CF(I))
            ELSE
               NI = -INT(-CF(I)) - 1
            END IF
C
            CF(I) = CF(I) - NI
C
            IF ( ABS(1D0-CF(I)).LT.TOL ) CF(I) = 0D0
         END DO
C
         CALL RVECLCRVB(3,CF,ABAS_DIM,QBAS_TMP(1,IQ))
C
C---------------------------------------------------------- check result
C
         QVEC_SHIFT(:) = QBAS_TMP(:,IQ) - QBAS(:,IQ)
C
         CALL RVECEXPAND(QVEC_SHIFT,ABAS_DIM,ADAINV_DIM,CF)
C
         DO I = 1,3
            IF ( ABS(CF(I)-NINT(CF(I))).GT.TOL ) THEN
               WRITE (6,*) ' STOP in <BASIS_INTO_UNIT_CELL>'
               WRITE (6,*) ' TROUBLE:   for IQ=',IQ,'   CF=',CF
               STOP
            END IF
         END DO
C
      END DO
C=======================================================================
C
      WRITE (6,99001)
      DO IQ = 1,NQ
         WRITE (6,99002) (QBAS(I,IQ),QBAS_TMP(I,IQ),I=1,3)
      END DO
C
C---------------------------------------------------------- store result
C
      IF ( OVERWRITE ) THEN
         QBAS(:,1:NQ) = QBAS_TMP(:,1:NQ)
      ELSE
         QBAS_UC(:,1:NQ) = QBAS_TMP(:,1:NQ)
      END IF
C
99001 FORMAT (//,2(1X,79('*'),/),29X,'<BASIS_INTO_UNIT_CELL>',/,
     &        2(1X,79('*')),//,10X,'basis vectors in units of a',/,10X,
     &        'original',6X,'witin unit cell')
99002 FORMAT (5X,'(',F10.5,',',F10.5,',',F10.5,' )',3X,'(',F10.5,',',
     &        F10.5,',',F10.5,' )')
      END
