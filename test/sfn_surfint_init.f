C*==sfn_surfint_init.f    processed by SPAG 6.70Rc at 13:55 on 27 Jan 2017
      SUBROUTINE SFN_SURFINT_INIT(NFACE_M,NEDGE_FCM,RVERT_VFM,NVERTMAX,
     &                            NFACEMAX)
C   ********************************************************************
C   *                                                                  *
C   *  find the intersections of the direction vectors                 *
C   *  connected with a triangular integration grid                    *
C   *  with the faces of an atomic cell                                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,DATSET0,LDATSET0,IOTMP
      USE MOD_RMESH,ONLY:NM,NMMAX,RHAT_VORONOI_SURF_GRID,
     &    W_VORONOI_SURF_GRID,D_VORONOI_SURF_GRID,N_VORONOI_SURF_GRID,
     &    NMAX_VORONOI_SURF_GRID
      USE MOD_CONSTANTS,ONLY:A0_ANG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFN_SURFINT_INIT')
      INTEGER N_GRID
      PARAMETER (N_GRID=6)
C
C Dummy arguments
C
      INTEGER NFACEMAX,NVERTMAX
      INTEGER NEDGE_FCM(NFACEMAX,NMMAX),NFACE_M(NMMAX)
      REAL*8 RVERT_VFM(3,NVERTMAX,NFACEMAX,NMMAX)
C
C Local variables
C
      REAL*8 AVEC(3),A_TRI,BVEC(3),COLOR,D,DVEC_VORONOI_SURF_GRID(:,:),
     &       RV(3,3),RVEC(3),S,WGHT(N_GRID),ZETA(3,N_GRID)
      REAL*8 DNRM2
      CHARACTER*80 FILPDB,FILRAS
      INTEGER I,IC,IC0,IC1,ICP,IFACE,IM,IVERT,I_GRID,
     &        I_VORONOI_SURF_GRID,J,LFIL,LL,N,NC
      CHARACTER*20 STR20
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DVEC_VORONOI_SURF_GRID
C
C----- use larger scaling factor to avoid spurious lines drawn by rasmol
      S = 3*A0_ANG*8D0
      S = 6*S
C
C=======================================================================
C        weights and grid positions for integration over a Triangle
C=======================================================================
C
      WGHT(1:3) = 0.109951743655322D0
C
      ZETA(1:3,1:3) = 0.091576213509771D0
      ZETA(1,1) = 0.816847572980459D0
      ZETA(2,2) = 0.816847572980459D0
      ZETA(3,3) = 0.816847572980459D0
C
      WGHT(4:6) = 0.223381589678011D0
C
      ZETA(1:3,4:6) = 0.445948490915965D0
      ZETA(1,4) = 0.108103018168070D0
      ZETA(2,5) = 0.108103018168070D0
      ZETA(3,6) = 0.108103018168070D0
C
C=======================================================================
C     get max number NMAX_VORONOI_SURF_GRID of gridpoints per mesh
C=======================================================================
C
      NMAX_VORONOI_SURF_GRID = 0
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      DO IM = 1,NM
         N = 0
         DO IFACE = 1,NFACE_M(IM)
            N = N + NEDGE_FCM(IFACE,IM)*N_GRID
         END DO
         NMAX_VORONOI_SURF_GRID = MAX(NMAX_VORONOI_SURF_GRID,N)
      END DO
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
C
      ALLOCATE (N_VORONOI_SURF_GRID(NMMAX))
      ALLOCATE (D_VORONOI_SURF_GRID(NMAX_VORONOI_SURF_GRID,NMMAX))
      ALLOCATE (W_VORONOI_SURF_GRID(NMAX_VORONOI_SURF_GRID,NMMAX))
      ALLOCATE (RHAT_VORONOI_SURF_GRID(3,NMAX_VORONOI_SURF_GRID,NMMAX))
      ALLOCATE (DVEC_VORONOI_SURF_GRID(3,NMAX_VORONOI_SURF_GRID))
C
C
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      DO IM = 1,NM
C
         I_VORONOI_SURF_GRID = 0
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         DO IFACE = 1,NFACE_M(IM)
C
            IF ( IPRINT.GT.0 ) WRITE (6,99001) IM,IFACE
C
C-----------------------------------------------------------------------
C              center of gravity
C-----------------------------------------------------------------------
C
            RVEC(1:3) = 0D0
            DO IVERT = 1,NEDGE_FCM(IFACE,IM)
               RVEC(1:3) = RVEC(1:3) + RVERT_VFM(1:3,IVERT,IFACE,IM)
            END DO
            RVEC(1:3) = RVEC(1:3)/DBLE(NEDGE_FCM(IFACE,IM))
C
            RV(1:3,1) = RVEC(1:3)
C
Cttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
C       run over all triangles with a corner at ->r_cg
Cttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
            DO IVERT = 1,NEDGE_FCM(IFACE,IM)
               RV(1:3,2) = RVERT_VFM(1:3,IVERT,IFACE,IM)
               IF ( IVERT.LT.NEDGE_FCM(IFACE,IM) ) THEN
                  RV(1:3,3) = RVERT_VFM(1:3,IVERT+1,IFACE,IM)
               ELSE
                  RV(1:3,3) = RVERT_VFM(1:3,1,IFACE,IM)
               END IF
C
               AVEC(1:3) = RV(1:3,1) - RV(1:3,2)
               BVEC(1:3) = RV(1:3,3) - RV(1:3,1)
               CALL RVECXPRO(AVEC,BVEC,RVEC)
C
               A_TRI = DABS(DNRM2(3,RVEC,1))/2D0
C
               DO I_GRID = 1,N_GRID
                  I_VORONOI_SURF_GRID = I_VORONOI_SURF_GRID + 1
                  I = I_VORONOI_SURF_GRID
                  RVEC(1:3) = 0D0
                  DO J = 1,3
                     RVEC(1:3) = RVEC(1:3) + ZETA(J,I_GRID)*RV(1:3,J)
                  END DO
                  D = DNRM2(3,RVEC,1)
C
                  W_VORONOI_SURF_GRID(I,IM) = WGHT(I_GRID)*A_TRI
                  D_VORONOI_SURF_GRID(I,IM) = D
                  RHAT_VORONOI_SURF_GRID(1:3,I,IM) = RVEC(1:3)/D
                  DVEC_VORONOI_SURF_GRID(1:3,I) = RVEC(1:3)
               END DO
C
            END DO
Cttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
C
         END DO
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
         N_VORONOI_SURF_GRID(IM) = I_VORONOI_SURF_GRID
C
C=======================================================================
C           plot central polyhedron + intersection points
C=======================================================================
C
         FILPDB = DATSET0(1:LDATSET0)//'_VSI_M'
         CALL STRING_ADD_N(FILPDB,IM)
         LFIL = LEN_TRIM(FILPDB)
         FILPDB = FILPDB(1:LFIL)//'.pdb'
         LFIL = LFIL + 4
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILPDB(1:LFIL))
         WRITE (IOTMP,99007) 
     &                      'polyhedron + surface integration grid for '
     &                      //DATSET0(1:LDATSET0)
C
C ----------------------------------------------------------------------
C                  Plotting of the intersection points of
C            surface integration grid and polyhedron surface
C ----------------------------------------------------------------------
C
         DO I_GRID = 1,N_VORONOI_SURF_GRID(IM)
C
            COLOR = 1.0D0
            WRITE (IOTMP,FMT=99003) I_GRID,I_GRID,
     &                              DVEC_VORONOI_SURF_GRID(1:3,I_GRID)
     &                              *S,COLOR
         END DO
C
C ----------------------------------------------------------------------
C                  Plotting of the polyhedron
C ----------------------------------------------------------------------
C
         COLOR = 2.0D0
         IC0 = N_VORONOI_SURF_GRID(IM)
         IC = IC0
         DO I = 1,NFACE_M(IM)
            DO J = 1,NEDGE_FCM(I,IM)
               IC = IC + 1
               WRITE (IOTMP,FMT=99004) IC,IC,RVERT_VFM(1:3,J,I,IM)*S,
     &                                 COLOR
            END DO
         END DO
C
         NC = IC
C
         IC = IC0
         DO I = 1,NFACE_M(IM)
            IC1 = IC + 1
            DO J = 1,NEDGE_FCM(I,IM)
               IC = IC + 1
               IF ( J.LT.NEDGE_FCM(I,IM) ) THEN
                  ICP = IC + 1
               ELSE
                  ICP = IC1
               END IF
               WRITE (IOTMP,99006) IC,ICP
            END DO
         END DO
C
         WRITE (IOTMP,99005)
C
         CLOSE (IOTMP)
C
C ----------------------------------------------------------------------
C                  write script file
C ----------------------------------------------------------------------
C
         FILRAS = FILPDB(1:(LFIL-4))//'.ras'
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILRAS(1:LFIL))
C
         WRITE (IOTMP,*) 'load '''//FILPDB(1:LFIL)//'''  '
         WRITE (IOTMP,*) 'set background white'
         WRITE (IOTMP,*) 'color temperature'
         WRITE (IOTMP,*) 'set fontsize 20'
         STR20 = 'select 1-'
         CALL STRING_ADD_N(STR20,N_VORONOI_SURF_GRID(IM))
         WRITE (IOTMP,*) STR20
         WRITE (IOTMP,*) 'cpk  150'
         STR20 = 'select '
         CALL STRING_ADD_N(STR20,N_VORONOI_SURF_GRID(IM)+1)
         LL = LEN_TRIM(STR20)
         STR20 = STR20(1:LL)//'-'
         CALL STRING_ADD_N(STR20,NC)
         WRITE (IOTMP,*) STR20
         WRITE (IOTMP,*) 'cpk  10'
         WRITE (IOTMP,*) 'select all '
         WRITE (IOTMP,*) 'set axes on '
C
         CLOSE (IOTMP)
C
      END DO
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
C
      WRITE (6,99002) ROUTINE(1:LEN_TRIM(ROUTINE))
C
99001 FORMAT (/,'MESH ',I3,'   FACE ',I3)
99002 FORMAT (/,1X,79('*'),/,36X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'for all directions crossing with surface found',/)
99003 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99004 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99005 FORMAT ('END  ')
99006 FORMAT ('CONECT',2I5)
99007 FORMAT ('HEADER    ',A,/,'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
      END
