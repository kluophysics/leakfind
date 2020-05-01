C*==symrotlat.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMROTLAT
C   ********************************************************************
C   *                                                                  *
C   *  rotate the lattice vectors    ->BR  and  ->Q   ABAS and QBAS    *
C   *  corresponding to the Euler angles  LALF, LBET, LGAM (in degree) *
C   *  according to the ACTIVE-TEMPORARY AXIS convention  (Rose)       *
C   *                                                                  *
C   *  NOTE: if the angles are not found the corresponding             *
C   *        information supplied by  MALF, MBET, MGAM  or  MDIR       *
C   *        is used instead if available in the input                 *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:ROTATE_LATTICE
      USE MOD_SITES,ONLY:QBAS,NQ
      USE MOD_LATTICE,ONLY:ABAS,SYSTEM_DIMENSION
      USE MOD_FILES,ONLY:IPRINT,RDUMMY,FOUND_SECTION,FOUND_REAL,
     &    FOUND_REAL_ARRAY,N_FOUND
      IMPLICIT NONE
C*--SYMROTLAT20
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ALFDEG,BETDEG,GAMDEG,LALF,LBET,LGAM,MDIR(3),MROT(3,3),V(3)
      LOGICAL ANGLEFOUND
      INTEGER I,IQ
C
C*** End of declarations rewritten by SPAG
C
      DATA ALFDEG/0D0/,BETDEG/0D0/,GAMDEG/0D0/
      DATA LALF/0D0/,LBET/0D0/,LGAM/0D0/
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( .NOT.FOUND_SECTION ) RETURN
C
      CALL SECTION_FIND_KEYWORD('ROTLAT',ROTATE_LATTICE)
C
C=======================================================================
C
      IF ( .NOT.ROTATE_LATTICE ) RETURN
C
C=======================================================================
C
      CALL SECTION_SET_REAL('LALF',LALF,9999D0,0)
      ANGLEFOUND = FOUND_REAL
      CALL SECTION_SET_REAL('LBET',LBET,9999D0,0)
      ANGLEFOUND = ANGLEFOUND .OR. FOUND_REAL
      CALL SECTION_SET_REAL('LGAM',LGAM,9999D0,0)
      ANGLEFOUND = ANGLEFOUND .OR. FOUND_REAL
C
C=======================================================================
C             try to get orientation of magnetic moment
C=======================================================================
      IF ( .NOT.ANGLEFOUND ) THEN
C
         CALL SECTION_SET_REAL('MALF',ALFDEG,9999D0,0)
         ANGLEFOUND = ANGLEFOUND .OR. FOUND_REAL
         CALL SECTION_SET_REAL('MBET',BETDEG,9999D0,0)
         ANGLEFOUND = ANGLEFOUND .OR. FOUND_REAL
         CALL SECTION_SET_REAL('MGAM',GAMDEG,9999D0,0)
         ANGLEFOUND = ANGLEFOUND .OR. FOUND_REAL
C
         IF ( .NOT.ANGLEFOUND ) THEN
            CALL SECTION_SET_REAL_ARRAY('MDIR',MDIR,N_FOUND,3,0,9999D0,
     &                                  0)
            ANGLEFOUND = FOUND_REAL_ARRAY
C
            IF ( ANGLEFOUND ) CALL CONVERT_CART_TO_SPHER(MDIR(1),MDIR(2)
     &           ,MDIR(3),RDUMMY,BETDEG,ALFDEG,.TRUE.)
         END IF
C
         LALF = -GAMDEG
         LBET = -BETDEG
         LGAM = -ALFDEG
C
      END IF
C
C=======================================================================
      IF ( .NOT.ANGLEFOUND ) THEN
         WRITE (6,99005)
         STOP
      END IF
C=======================================================================
C
      WRITE (6,99001) LALF,LBET,LGAM
C
      CALL GETMROT(LALF,LBET,LGAM,MROT)
C
      WRITE (6,99002)
      DO I = 1,3
         V(1:3) = ABAS(1:3,I)
         ABAS(1:3,I) = MATMUL(MROT(1:3,1:3),V(1:3))
         WRITE (6,99004) I,V,ABAS(1:3,I)
      END DO
C
      WRITE (6,99003)
      DO IQ = 1,NQ
         V(1:3) = QBAS(1:3,IQ)
         QBAS(1:3,IQ) = MATMUL(MROT(1:3,1:3),V(1:3))
         WRITE (6,99004) IQ,V,QBAS(1:3,IQ)
      END DO
C
C=======================================================================
C           find a suitable Bravais lattice w.r.t. the symmetry
C=======================================================================
C
      CALL INIT_MOD_LATTICE(IPRINT,SYSTEM_DIMENSION)
C
      CALL SYMBRAVAIS
C
      RETURN
99001 FORMAT (//,1X,79('*'),/,34X,'<SYMROTLAT>',/,1X,79('*'),//,10X,
     &        'rotate the spatial lattice',/10X,
     &        'according to the  Euler-angles: ',3F6.1)
99002 FORMAT (/,10X,'primitive lattice vectors  ->A_i ',/)
99003 FORMAT (/,10X,'basis vectors  Q',/)
99004 FORMAT (10X,I3,3X,2('(',F5.2,',',F5.2,',',F5.2,')',:,'  ->  '))
99005 FORMAT (//,1X,79('#'),/,34X,'<SYMROTLAT>',/,1X,79('#'),//,10X,
     &        'ROTLAT  set but no Euler angles given in section MODE',/,
     &        10X,
     &        'specify angles:  LALF=xxx LBET=xxx LGAM=xxx in degrees',
     &        /,1X,79('#'),/)
      END
C*==convert_cart_to_spher.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE CONVERT_CART_TO_SPHER(X,Y,Z,R,TET,PHI,DEG)
C   ********************************************************************
C   *                                                                  *
C   *  convert cartesian coordinates to spherical ones                 *
C   *  give angles in degrees if requested by  DEG                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--CONVERT_CART_TO_SPHER145
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL DEG
      REAL*8 PHI,R,TET,X,Y,Z
C
C Local variables
C
      REAL*8 D
C
C*** End of declarations rewritten by SPAG
C
      R = SQRT(X*X+Y*Y+Z*Z)
C
      D = SQRT(X*X+Y*Y)
C
C------------------------------------------------------------- x = y = 0
      IF ( ABS(D).LT.1D-8 ) THEN
C
         PHI = 0D0
C
C---------------------------------------------------------------- x >= 0
      ELSE IF ( X.GE.0 ) THEN
C
         IF ( Y.GE.0 ) THEN
            PHI = ACOS(X/D)
         ELSE
            PHI = 2*PI - ACOS(X/D)
         END IF
C
C----------------------------------------------------------------- x < 0
      ELSE IF ( Y.GE.0 ) THEN
         PHI = PI - ACOS(-X/D)
      ELSE
         PHI = PI + ACOS(-X/D)
      END IF
C-----------------------------------------------------------------------
C
      TET = ACOS(Z/R)
C
      IF ( DEG ) THEN
         TET = TET*180D0/PI
         PHI = PHI*180D0/PI
      END IF
C
      END
