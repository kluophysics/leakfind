C*==kdirtab.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE KDIRTAB(BRAVAIS,KPATH,NKDIR,LBLKDIR,KA,KE,BBAS)
C **********************************************************************
C *                                                                    *
C *   generate a path in the reciprocal lattice according to           *
C *   BRAVAIS    and    KPATH                                          *
C *   the path may consists of several segments which do not have      *
C *   to join continously                                              *
C *   NKDIR   number of k-path segments created                        *
C *   LBLKDIR name of each segment - contains name of first, last      *
C *           and intermediate k-point (symmetry point in BZ)          *
C *                                                                    *
C *   in multiples of  2 * PI/A                                        *
C *                                                                    *
C *   KPATH = 10  collect all edges of the irreducible edge of the BZ  *
C *                                                                    *
C *--------------------------------------------------------------------*
C *  BRAVAIS                             KPATH    NKDIR   LBLKDIR      *
C *                                                                    *
C *      1 triclinic   primitive                                       *
C *      2 monoclinic  primitive                                       *
C *      3 monoclinic  base centered                                   *
C *      4 orthorombic primitive                                       *
C *      5 orthorombic base-centered                                   *
C *      6 orthorombic body-centered                                   *
C *      7 orthorombic face-centered                                   *
C *      8 tetragonal  primitive                                       *
C *      9 tetragonal  body-centered                                   *
C *     10 trigonal    primitive                                       *
C *     11 hexagonal   primitive                                       *
C *     12 cubic       primitive                                       *
C *     13 cubic       face-centered                                   *
C *     14 cubic       body-centered                                   *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C*--KDIRTAB37
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 R0,R1,R2,R3,R4,R5,R8
      PARAMETER (R0=0.0D0,R1=1.0D0,R2=2.0D0,R3=3.0D0,R4=4.0D0,R5=5.0D0,
     &           R8=8.0D0)
C
C Dummy arguments
C
      INTEGER BRAVAIS,KPATH,NKDIR
      REAL*8 BBAS(3,3),KA(3,20),KE(3,20)
      CHARACTER*8 LBLKDIR(20)
C
C Local variables
C
      REAL*8 AVEC(3),B1VEC(3),B2VEC(3),B3VEC(3),BVEC(3),CVEC(3),DVEC(3),
     &       EVEC(3),GAMV(3),HVEC(3),KVEC(3),LVEC(3),MVEC(3),NVEC(3),
     &       PVEC(3),RVEC(3),SVEC(3),TVEC(3),UVEC(3),WVEC(3),XVEC(3),
     &       YVEC(3),ZVEC(3)
      INTEGER I,ID
C
C*** End of declarations rewritten by SPAG
C
      B1VEC(1:3) = BBAS(1:3,1)
      B2VEC(1:3) = BBAS(1:3,2)
      B3VEC(1:3) = BBAS(1:3,3)
C
      CALL LC3RVEC(GAMV,R0,R0,R0,B1VEC,B2VEC,B3VEC,3)
C
      NKDIR = 0
      ID = 0
C
C----------------------------------------------- 1 triclinic   primitive
      IF ( BRAVAIS.EQ.1 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C----------------------------------------------- 2 monoclinic  primitive
      ELSE IF ( BRAVAIS.EQ.2 ) THEN
         CALL LC3RVEC(BVEC,-R1/R2,R0,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(YVEC,R0,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(ZVEC,R0,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(CVEC,R0,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(DVEC,-R1/R2,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(AVEC,-R1/R2,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(EVEC,-R1/R2,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
C
         STOP '<KDIRTAB>: Bravais lattice not treated '
C------------------------------------------- 3 monoclinic  base centered
      ELSE IF ( BRAVAIS.EQ.3 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C----------------------------------------------- 4 orthorombic primitive
      ELSE IF ( BRAVAIS.EQ.4 ) THEN
         CALL LC3RVEC(YVEC,R0,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(XVEC,R1/R2,R0,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(ZVEC,R0,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(UVEC,R1/R2,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(TVEC,R0,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(SVEC,R1/R2,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(RVEC,R1/R2,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
C
         IF ( KPATH.EQ.1 ) NKDIR = 12
         IF ( KPATH.EQ.2 ) NKDIR = 7
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) THEN
            NKDIR = 3
            GOTO 50
         END IF
         IF ( KPATH.EQ.5 ) NKDIR = 1
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 1
            GOTO 50
         END IF
         IF ( KPATH.EQ.7 ) THEN
            NKDIR = 2
            GOTO 50
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-X '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -G -U '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -A -Z '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,ZVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'Z -GL-GG'
         CALL DCOPY(3,ZVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
 50      CONTINUE
         IF ( KPATH.EQ.7 ) THEN
            ID = ID + 1
            LBLKDIR(ID) = 'X -GS-GG'
            CALL DCOPY(3,XVEC,1,KA(1,ID),1)
            CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-Y '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,YVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'Y -H -T '
         CALL DCOPY(3,YVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'T -B -Z '
         CALL DCOPY(3,TVEC,1,KA(1,ID),1)
         CALL DCOPY(3,ZVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -D -S '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,SVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'S -C -Y '
         CALL DCOPY(3,SVEC,1,KA(1,ID),1)
         CALL DCOPY(3,YVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -P -R '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,RVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'R -E -T '
         CALL DCOPY(3,RVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'S -Q -T '
         CALL DCOPY(3,SVEC,1,KA(1,ID),1)
         CALL DCOPY(3,TVEC,1,KE(1,ID),1)
C------------------------------------------- 5 orthorombic base-centered
      ELSE IF ( BRAVAIS.EQ.5 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C------------------------------------------- 6 orthorombic body-centered
      ELSE IF ( BRAVAIS.EQ.6 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C------------------------------------------- 7 orthorombic face-centered
      ELSE IF ( BRAVAIS.EQ.7 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C----------------------------------------------- 8 tetragonal  primitive
      ELSE IF ( BRAVAIS.EQ.8 ) THEN
         CALL LC3RVEC(MVEC,R1/R2,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(ZVEC,R0,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(AVEC,R1/R2,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(RVEC,R0,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(XVEC,R0,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
C
         STOP '<KDIRTAB>: Bravais lattice not treated '
C------------------------------------------- 9 tetragonal  body-centered
      ELSE IF ( BRAVAIS.EQ.9 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C---------------------------------------------- 10 trigonal    primitive
      ELSE IF ( BRAVAIS.EQ.10 ) THEN
         STOP '<KDIRTAB>: Bravais lattice not treated '
C---------------------------------------------- 11 hexagonal   primitive
      ELSE IF ( BRAVAIS.EQ.11 ) THEN
         CALL LC3RVEC(MVEC,R0,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(AVEC,R0,R0,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(LVEC,R0,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(KVEC,-R1/R3,R2/R3,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(HVEC,-R1/R3,R2/R3,R1/R2,B1VEC,B2VEC,B3VEC,3)
C
         IF ( KPATH.EQ.1 ) NKDIR = 9
         IF ( KPATH.EQ.2 ) NKDIR = 7
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) NKDIR = 1
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 1
            GOTO 100
         END IF
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-M '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -T''-K '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,KVEC,1,KE(1,ID),1)
 100     CONTINUE
         ID = ID + 1
         LBLKDIR(ID) = 'K -T -GG'
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-A '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,AVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'A -R -L '
         CALL DCOPY(3,AVEC,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -S''-H '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'H -S -A '
         CALL DCOPY(3,HVEC,1,KA(1,ID),1)
         CALL DCOPY(3,AVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -U -L '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'K -P -H '
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
C---------------------------------------------- 12 cubic       primitive
      ELSE IF ( BRAVAIS.EQ.12 ) THEN
         CALL LC3RVEC(XVEC,R0,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(MVEC,R1/R2,R1/R2,R0,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(RVEC,R1/R2,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
C
         NKDIR = 1
         IF ( KPATH.EQ.1 ) NKDIR = 5
         IF ( KPATH.EQ.2 ) NKDIR = 4
         IF ( KPATH.EQ.3 ) NKDIR = 3
         IF ( KPATH.EQ.4 ) NKDIR = 2
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 3
            DO ID = 1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL RINIT(3,KE(1,ID))
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 0.5D0
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 0.5D0
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 0.5D0
            ID = NKDIR
         END IF
C
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-X '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -Y -M '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'M -V -R '
         CALL DCOPY(3,MVEC,1,KA(1,ID),1)
         CALL DCOPY(3,RVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'R -GL-GG'
         CALL DCOPY(3,RVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GS-M '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,MVEC,1,KE(1,ID),1)
C------------------------------------------ 13 cubic       face-centered
      ELSE IF ( BRAVAIS.EQ.13 ) THEN
         CALL LC3RVEC(XVEC,R0,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(LVEC,R1/R2,R1/R2,R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(WVEC,R1/R4,R1/R2,R3/R4,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(KVEC,R3/R8,R3/R8,R3/R4,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(UVEC,R1/R4,R5/R8,R5/R8,B1VEC,B2VEC,B3VEC,3)
C
         IF ( KPATH.EQ.0 ) KPATH = 4
         IF ( KPATH.EQ.1 ) NKDIR = 9
         IF ( KPATH.EQ.2 ) NKDIR = 5
         IF ( KPATH.EQ.3 ) NKDIR = 2
         IF ( KPATH.EQ.4 ) THEN
            NKDIR = 1
            ID = ID + 1
            LBLKDIR(ID) = 'GG-GD-X '
            CALL DCOPY(3,GAMV,1,KA(1,ID),1)
            CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         END IF
         IF ( KPATH.EQ.5 ) THEN
            NKDIR = 1
            GOTO 150
         END IF
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 3
            DO ID = 1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL RINIT(3,KE(1,ID))
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 1D0
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 1D0
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 1D0
            ID = NKDIR
         END IF
C
         ID = ID + 1
         LBLKDIR(ID) = 'X -GD-GG'
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
 150     CONTINUE
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GL-L '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,LVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -Q -W '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,WVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'W -N -K '
         CALL DCOPY(3,WVEC,1,KA(1,ID),1)
         CALL DCOPY(3,KVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'K -GS-GG'
         CALL DCOPY(3,KVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'L -M -U '
         CALL DCOPY(3,LVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'U -S -X '
         CALL DCOPY(3,UVEC,1,KA(1,ID),1)
         CALL DCOPY(3,XVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'X -Z -W '
         CALL DCOPY(3,XVEC,1,KA(1,ID),1)
         CALL DCOPY(3,WVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'W -B -U '
         CALL DCOPY(3,WVEC,1,KA(1,ID),1)
         CALL DCOPY(3,UVEC,1,KE(1,ID),1)
         IF ( KPATH.EQ.10 ) THEN
            NKDIR = 9
            ID = 3
            LBLKDIR(ID) = 'L -? -K '
            CALL DCOPY(3,LVEC,1,KA(1,ID),1)
            CALL DCOPY(3,KVEC,1,KE(1,ID),1)
         END IF
C------------------------------------------ 14 cubic       body-centered
      ELSE IF ( BRAVAIS.EQ.14 ) THEN
         CALL LC3RVEC(HVEC,R1/R2,R1/R2,-R1/R2,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(PVEC,R1/R4,R1/R4,R1/R4,B1VEC,B2VEC,B3VEC,3)
         CALL LC3RVEC(NVEC,R1/R2,R0,R0,B1VEC,B2VEC,B3VEC,3)
C
         IF ( KPATH.EQ.0 ) KPATH = 5
         IF ( KPATH.EQ.1 ) NKDIR = 6
         IF ( KPATH.EQ.2 ) NKDIR = 5
         IF ( KPATH.EQ.3 ) NKDIR = 4
         IF ( KPATH.EQ.4 ) NKDIR = 3
         IF ( KPATH.EQ.5 ) NKDIR = 1
         IF ( KPATH.EQ.6 ) THEN
            NKDIR = 3
            DO ID = 1,NKDIR
               CALL DCOPY(3,GAMV,1,KA(1,ID),1)
               CALL RINIT(3,KE(1,ID))
            END DO
            LBLKDIR(1) = 'GG-GD-X '
            KE(1,1) = 1D0
            LBLKDIR(2) = 'GG-GD-Y '
            KE(2,2) = 1D0
            LBLKDIR(3) = 'GG-GD-Z '
            KE(3,3) = 1D0
            ID = NKDIR
         END IF
C
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GD-H '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'H -G -N '
         CALL DCOPY(3,HVEC,1,KA(1,ID),1)
         CALL DCOPY(3,NVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'N -GS-GG'
         CALL DCOPY(3,NVEC,1,KA(1,ID),1)
         CALL DCOPY(3,GAMV,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'GG-GL-P '
         CALL DCOPY(3,GAMV,1,KA(1,ID),1)
         CALL DCOPY(3,PVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'P -F -H '
         CALL DCOPY(3,PVEC,1,KA(1,ID),1)
         CALL DCOPY(3,HVEC,1,KE(1,ID),1)
         ID = ID + 1
         LBLKDIR(ID) = 'N -D -P '
         CALL DCOPY(3,NVEC,1,KA(1,ID),1)
         CALL DCOPY(3,PVEC,1,KE(1,ID),1)
      ELSE
         WRITE (6,*) '<KDIRTAB> called for BRAVAIS=',BRAVAIS
         STOP
      END IF
C ......................................................................
      WRITE (6,99001) KPATH,NKDIR
      WRITE (6,99002) '->B(1): ',B1VEC
      WRITE (6,99002) '->B(2): ',B2VEC
      WRITE (6,99002) '->B(3): ',B3VEC
      WRITE (6,*) ' '
      DO ID = 1,NKDIR
         WRITE (6,99002) LBLKDIR(ID),(KA(I,ID),I=1,3),(KE(I,ID),I=1,3)
      END DO
      WRITE (6,*) ' '
99001 FORMAT (/,1X,79('*'),/,35X,'<KDIRTAB>',/,1X,79('*'),//,5X,
     &        'for KPATH =',I2,I10,' k-directions created',/)
99002 FORMAT (5X,A,2X,'(',3F8.4,' )',:,'   ...   (',3F8.4,' )')
      END
C*==lc3rvec.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LC3RVEC(V,C1,C2,C3,V1,V2,V3,L)
C **********************************************************************
C *                                                                    *
C *     ->V = C1 * ->V1 + C2 * ->V2 + C3 * ->V3                        *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT NONE
C*--LC3RVEC467
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 C1,C2,C3
      INTEGER L
      REAL*8 V(L),V1(L),V2(L),V3(L)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,L
         V(I) = C1*V1(I) + C2*V2(I) + C3*V3(I)
      END DO
      END
