C*==mod_lattice.f    processed by SPAG 6.70Rc at 09:20 on  3 Aug 2011
      MODULE MOD_LATTICE
C   ********************************************************************
C   *                                                                  *
C   *  module to store all variables connected with files              *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--MOD_LATTICE11
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 A5BAS(3,5),ABAS(3,3),ABAS_2D(3,3),ABAS_I(3,3),ABAS_L(3,3),
     &       ABAS_R(3,3),ADAINV(3,3),ADAINV_I(3,3),ADAINV_L(3,3),
     &       ADAINV_R(3,3),ADAMAT(3,3),ADAMAT_I(3,3),ADAMAT_L(3,3),
     &       ADAMAT_R(3,3),ALAT,BBAS(3,3),BBAS_2D(3,3),BBAS_I(3,3),
     &       BBAS_L(3,3),BBAS_R(3,3),BDBINV(3,3),BDBMAT(3,3),BOA,COA,
     &       SWS,VOLUC,VOLUC_2D
      INTEGER BRAVAIS
      CHARACTER*10 SYSTEM_DIMENSION
      CHARACTER*30 SUB_SYSTEM,SYSTEM_TYPE
      CHARACTER*38 TXTBRAVAIS(14)
      SAVE A5BAS,ABAS,ABAS_2D,ABAS_I,ABAS_L,ABAS_R,ADAINV,ADAINV_I,
     &     ADAINV_L,ADAINV_R,ADAMAT,ADAMAT_I,ADAMAT_L,ADAMAT_R,ALAT,
     &     BBAS,BBAS_2D,BBAS_I,BBAS_L,BBAS_R,BDBINV,BDBMAT,BRAVAIS,SWS,
     &     VOLUC,VOLUC_2D
C
C*** End of declarations rewritten by SPAG
C
      DATA SYSTEM_DIMENSION/'3D        '/,SYSTEM_TYPE/'BULK      '/
      DATA BOA/999999D0/,COA/999999D0/,SUB_SYSTEM/'BULK      '/
      DATA TXTBRAVAIS/'triclinic   primitive      -1     C_i ',
     &     'monoclinic  primitive      2/m    C_2h',
     &     'monoclinic  base centered  2/m    C_2h',
     &     'orthorombic primitive      mmm    D_2h',
     &     'orthorombic base-centered  mmm    D_2h',
     &     'orthorombic body-centered  mmm    D_2h',
     &     'orthorombic face-centered  mmm    D_2h',
     &     'tetragonal  primitive      4/mmm  D_4h',
     &     'tetragonal  body-centered  4/mmm  D_4h',
     &     'trigonal    primitive      -3m    D_3d',
     &     'hexagonal   primitive      6/mmm  D_6h',
     &     'cubic       primitive      m3m    O_h ',
     &     'cubic       face-centered  m3m    O_h ',
     &     'cubic       body-centered  m3m    O_h '/
C
      END
