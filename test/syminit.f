C*==syminit.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMINIT(IPRINT,BRAVAIS,NSYM,NSYMCRYSYS,SYMCRYSYS,
     &                   SYMSYMBL,SYMEULANG)
C   ********************************************************************
C   *                                                                  *
C   *   The main task of the routine is to set up the symmetry         *
C   *   operations compatible with the selected Bravais lattice        *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *         1 triclinic   primitive      -1     C_i                  *
C   *         2 monoclinic  primitive      2/m    C_2h                 *
C   *         3 monoclinic  base centered  2/m    C_2h                 *
C   *         4 orthorombic primitive      mmm    D_2h                 *
C   *         5 orthorombic base-centered  mmm    D_2h                 *
C   *         6 orthorombic body-centered  mmm    D_2h                 *
C   *         7 orthorombic face-centered  mmm    D_2h                 *
C   *         8 tetragonal  primitive      4/mmm  D_4h                 *
C   *         9 tetragonal  body-centered  4/mmm  D_4h                 *
C   *        10 trigonal    primitive      -3m    D_3d                 *
C   *        11 hexagonal   primitive      6/mmm  D_6h                 *
C   *        12 cubic       primitive      m3m    O_h                  *
C   *        13 cubic       face-centered  m3m    O_h                  *
C   *        14 cubic       body-centered  m3m    O_h                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NO_SYMMETRY,NSYMMAX
      USE MOD_LATTICE,ONLY:TXTBRAVAIS
      IMPLICIT NONE
C*--SYMINIT31
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER BRAVAIS,IPRINT,NSYM,NSYMCRYSYS
      LOGICAL SYMCRYSYS(NSYMMAX)
      REAL*8 SYMEULANG(3,NSYMMAX)
      CHARACTER*4 SYMSYMBL(NSYMMAX)
C
C Local variables
C
      CHARACTER*3 CRYSYSNAM(7)
      REAL*8 EULANGCUB(3,24),EULANGHEX(3,12)
      INTEGER ICRYSYS,ICRYSYSBRAVAIS(14),ISYM,ISYMCRYSYS(24,7),NSYMH,
     &        NSYMTOP(7)
      CHARACTER*4 SYMSYMCUB(NSYMMAX),SYMSYMHEX(24)
C
C*** End of declarations rewritten by SPAG
C
      DATA CRYSYSNAM/'TCL','MCL','ORB','TET','TRI','HEX','CUB'/
C
C--------------------- gives the crystal system for each bravais lattice
      DATA ICRYSYSBRAVAIS/1,2,2,3,3,3,3,4,4,5,6,7,7,7/
C--------------------------- number of symmetry operations to be scanned
C---------- indicates also that the point group is subgroup of Oh or D6h
      DATA NSYMTOP/48,48,48,48,24,24,48/
C--------------- activates symmetry operations for a given crystal class
C-------------- the operations including the inversion are excluded here
C
C triclinic
C monoclinic
C orthorhombic
C tetragonal
C trigonal
C hexagonal
C cubic
      DATA ISYMCRYSYS/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
     &     0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,
     &     0,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,
     &     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
C
C=======================================================================
C the following tables give the Euler angles corresponding to the
C proper rotations of the cubic group Oh and the hexagonal group D6h
C taken from Bradley & Cracknell's table 1.4.
C Note: These authors use active fixed rotation axes
C       Rose's convention used here is to use ACTIVE TEMPORARY AXES
C       accordingly alpha and gamme had to be interchanged
C
      DATA EULANGCUB/0.0D0,0.0D0,0.0D0,0.0D0,180.0D0,180.0D0,0.0D0,
     &     180.0D0,0.0D0,0.0D0,0.0D0,180.0D0,0.0D0,90.0D0,90.0D0,
     &     180.0D0,90.0D0,270.0D0,180.0D0,90.0D0,90.0D0,0.0D0,90.0D0,
     &     270.0D0,90.0D0,90.0D0,180.0D0,270.0D0,90.0D0,0.0D0,90.0D0,
     &     90.0D0,0.0D0,270.0D0,90.0D0,180.0D0,270.0D0,90.0D0,90.0D0,
     &     0.0D0,90.0D0,0.0D0,0.0D0,0.0D0,90.0D0,90.0D0,90.0D0,270.0D0,
     &     180.0D0,90.0D0,180.0D0,0.0D0,0.0D0,270.0D0,0.0D0,180.0D0,
     &     90.0D0,0.0D0,180.0D0,270.0D0,0.0D0,90.0D0,180.0D0,90.0D0,
     &     90.0D0,90.0D0,180.0D0,90.0D0,0.0D0,270.0D0,90.0D0,270.0D0/
C
      DATA SYMSYMCUB/'E   ','C2x ','C2y ','C2z ','C+31','C+32','C+33',
     &     'C+34','C-31','C-32','C-33','C-34','C+4x','C+4y','C+4z',
     &     'C-4x','C-4y','C-4z','C2a ','C2b ','C2c ','C2d ','C2e ',
     &     'C2f ','I   ','sx  ','sy  ','sz  ','S-61','S-62','S-63',
     &     'S-64','S+61','S+62','S+63','S+64','S-4x','S-4y','S-4z',
     &     'S+4x','S+4y','S+4z','sda ','sdb ','sdc ','sdd ','sde ',
     &     'sdf '/
C
      DATA EULANGHEX/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,60.0D0,0.0D0,0.0D0,
     &     120.0D0,0.0D0,0.0D0,180.0D0,0.0D0,0.0D0,240.0D0,0.0D0,0.0D0,
     &     300.0D0,0.0D0,180.0D0,180.0D0,0.0D0,180.0D0,300.0D0,0.0D0,
     &     180.0D0,60.0D0,0.0D0,180.0D0,0.0D0,0.0D0,180.0D0,120.0D0,
     &     0.0D0,180.0D0,240.0D0/
C
      DATA SYMSYMHEX/'E   ','C+6 ','C+3 ','C2  ','C-3 ','C-6 ','C''21',
     &     'C''22','C''23','C"21','C"22','C"23','I   ','S-3 ','S-6 ',
     &     'sh  ','S+6 ','S+3 ','sd1 ','sd2 ','sd3 ','sv1 ','sv2 ',
     &     'sv3 '/
C
      IF ( BRAVAIS.LT.1 .OR. BRAVAIS.GT.14 ) THEN
         WRITE (6,*) ' ********************************* '
         WRITE (6,*) ' STOP IN  <SYMINIT> '
         WRITE (6,*) ' BRAVAIS     = ',BRAVAIS,TXTBRAVAIS(BRAVAIS)
         STOP
      END IF
C
      IF ( IPRINT.GE.0 ) WRITE (6,99001) TXTBRAVAIS(BRAVAIS)
C
C=======================================================================
C        set the symmetry operation flag  SYMCRYSYS
C          according to the crystal system  ICRYSYS
C        copy the Euler angles and the symmetry operation symbols
C=======================================================================
C
      NSYMCRYSYS = 0
      ICRYSYS = ICRYSYSBRAVAIS(BRAVAIS)
      NSYM = NSYMTOP(ICRYSYS)
      NSYMH = NSYM/2
C
      DO ISYM = 1,NSYMH
C
         IF ( ISYMCRYSYS(ISYM,ICRYSYS).EQ.1 ) THEN
            NSYMCRYSYS = NSYMCRYSYS + 1
            SYMCRYSYS(ISYM) = .TRUE.
            SYMCRYSYS(ISYM+NSYMH) = .TRUE.
         ELSE
            SYMCRYSYS(ISYM) = .FALSE.
            SYMCRYSYS(ISYM+NSYMH) = .FALSE.
         END IF
C
         IF ( NSYM.EQ.48 ) THEN
            CALL DCOPY(3,EULANGCUB(1,ISYM),1,SYMEULANG(1,ISYM),1)
            SYMSYMBL(ISYM) = SYMSYMCUB(ISYM)
            SYMSYMBL(ISYM+NSYMH) = SYMSYMCUB(ISYM+NSYMH)
         ELSE
            CALL DCOPY(3,EULANGHEX(1,ISYM),1,SYMEULANG(1,ISYM),1)
            SYMSYMBL(ISYM) = SYMSYMHEX(ISYM)
            SYMSYMBL(ISYM+NSYMH) = SYMSYMHEX(ISYM+NSYMH)
         END IF
C
      END DO
C
      NSYMCRYSYS = NSYMCRYSYS*2
C
      IF ( NO_SYMMETRY ) THEN
         SYMCRYSYS(2:NSYM) = .FALSE.
         NSYMCRYSYS = 1
         WRITE (6,99003)
      END IF
C
      IF ( IPRINT.GE.0 ) WRITE (6,99002) ICRYSYS,CRYSYSNAM(ICRYSYS),
     &                          NSYM,NSYMCRYSYS,
     &                          (SYMCRYSYS(ISYM),ISYM=1,NSYM)
C
      RETURN
99001 FORMAT (//,1X,79('*'),/,35X,'<SYMINIT>',/,1X,79('*'),//,10X,
     &        'set up possible symmetry operations for the ',/,10X,
     &        'real space Bravais lattice: ',A)
99002 FORMAT (/,10X,'crystal system ',I3,5X,A,//,10X,
     &        'symmetry operations to scan       NSYM ',I3,/10X,
     &        'maximum possible symmetry operations   ',I3,//,(10X,24L2)
     &        )
99003 FORMAT (2(/,1X,79('*')),/,10X,
     &        'use of symmetry suppressed for testing',/,2(/,1X,79('*'))
     &        ,/)
      END
