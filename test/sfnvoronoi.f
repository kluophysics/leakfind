C*==sfnvoronoi.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNVORONOI(NVEC,RVEC,NEDGEMAX,WEIGHT0,WEIGHT,RMT,ROUT,
     &                      VOLUME,NFACE,A3,B3,C3,D3,NEDGE,XEDGE,YEDGE,
     &                      ZEDGE,NFACEMAX,IFLAG_VERTEX)
C***********************************************************************
C
C    this file contains the voronoi program as well as all auxilary
C    subroutines and functions. the following modifications were made:
C
C    - all routines renamed to  SFN...
C    - all COMMON blocks removed
C    - all include statements removed
C    - all PARAMETER statements removed exept in  SFNVORONOI
C      the argument lists have been extended accordingly
C    - all temporary arrays are accessed using   ALLOCATE
C
C    HE 2006
C----------------------------------------------------------------------
C Given a cluster of atomic positions at RVEC(3,NVEC), this subroutine
C returns information about the Voronoi cell around the origin. It is
C supposed, of course, that the origin corresponds to an atomic position
C which is not included in the RVEC(3,N). The information returned is:
C
C RMT: Muffin-tin radius (radius of inscribed sphere centered at the
C      origin).
C
C ROUT: Radius of circumscribed sphere centered at the origin.
C
C VOLUME: Volume of the Voronoi cell.
C
C NFACE: Number of faces of the cell.
C
C A3(NFACE),B3(NFACE),C3(NFACE),D3(NFACE): Coefficients defining the
C faces via A3*x+B3*y+C3*z=D3. The arrays are filled in the first
C NFACE positions with the information, the rest can be garbage.
C
C NEDGE(NFACE): Number of edges of each face.
C
C XEDGE(NEDGEMAX,NFACE),YEDGE(NEDGEMAX,NFACE),ZEDGE(NEDGEMAX,NFACE):
C Coordinates of these edges.
C
C The Voronoi construction performed here allows for different than
C 50%/50% bisections. For this, each atomic site, positioned say at
C vector r(i), is assigned a weight w(i). Then a point r in space
C belongs to the cell i if, for all j,
C
C        |r-r(i)|**2 - w(i) < |r-r(j)|**2 - w(j).    (1)
C
C The points for which the unequality becomes an equality is the
C bisector, and it can be easily shown that it is a plane perpendicular
C to the vector r(j)-r(i).
C One can parametrize the segment connecting r(i) and r(j) using a
C parameter t, 0 < t < 1. For t=0 we are at r(i), for t=1 at r(j), and
C in-between we have a linear dependence of the distance from r(i) with
C t. I.e., the position vector is r = r(i)*(1-t) + r(j)*t.
C The point where the bisector cuts this segment will then be at
C
C        t=(1/2)*( (w(i)-w(j))/dist**2 + 1 )         (2)
C
C where dist is the distance of the two points. As a special case we see
C that, for w(i)=w(j), the bisector will be in the middle, else the
C atomic site with the bigger weight gains more space.
C
C The above procedure is guaranteed to divide space into tesselating
C (:=space-filling) convex polyhedra (one of their sides could be at
C infinity in special cases). The convexity of the polyhedra is
C guaranteed by the constuction, i.e. as mutual cut of half-spaces, and
C the property of tesselation is guaranteed by the fact that every point
C in space is assigned to some atom.
C
C However, this procedure might place the polyhedron in such a way that
C its own atomic site is not contained in it. This can happen if there
C is a big enough difference in the weights, whence, from eq.(2) above,
C t can become <0 or >1 (then the bisector has passed one of the two
C points). This cannot be allowed in a KKR calculation, and therefore
C each time it is checked that this is not the case. If such a case
C occurs, a different set of weights must be chosen.
C
C Uses subroutines NORMALPLANE, POLYHEDRON, and function SFNDISTPLANE.
C
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--SFNVORONOI83
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNVORONOI')
      INTEGER NPOIMAX,NNEIMAX,NPLANEMAX
      PARAMETER (NPOIMAX=1000,NNEIMAX=600,NPLANEMAX=600)
C
C Dummy arguments
C
      INTEGER IFLAG_VERTEX,NEDGEMAX,NFACE,NFACEMAX,NVEC
      REAL*8 RMT,ROUT,VOLUME,WEIGHT0
      REAL*8 A3(NFACEMAX),B3(NFACEMAX),C3(NFACEMAX),D3(NFACEMAX),
     &       RVEC(3,NVEC),WEIGHT(NVEC),XEDGE(NEDGEMAX,NFACEMAX),
     &       YEDGE(NEDGEMAX,NFACEMAX),ZEDGE(NEDGEMAX,NFACEMAX)
      INTEGER NEDGE(NFACEMAX)
C
C Local variables
C
      INTEGER I,IEDGE,IFACE,IVEC,NPLANE,NVERTEX
      REAL*8 RSQ,TAU,TEMP,TETRVOL,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3
      REAL*8 SFNDISTPLANE
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------------------------------
C Check that the origin is not included in RVEC.
C
      DO IVEC = 1,NVEC
         IF ( DABS(RVEC(1,IVEC))+DABS(RVEC(2,IVEC))+DABS(RVEC(3,IVEC))
     &        .LT.1.D-12 ) THEN
            WRITE (6,'(/,10X,A,I8,3X,A,/)') 'vector',IVEC,'is zero.'
            CALL STOP_MESSAGE(ROUTINE,'INCONSISTENCIES')
         END IF
      END DO
C
C---------------------------------------------------------------
C Define the planes as normal to the vectors RVEC, passing from t as
C in eq. (2) above:
C
      DO IVEC = 1,NVEC
C
         RSQ = RVEC(1,IVEC)*RVEC(1,IVEC) + RVEC(2,IVEC)*RVEC(2,IVEC)
     &         + RVEC(3,IVEC)*RVEC(3,IVEC)
         TAU = 0.5D0*((WEIGHT0-WEIGHT(IVEC))/RSQ+1.D0)
C
         CALL SFNNORMALPLANE(0.D0,0.D0,0.D0,RVEC(1,IVEC),RVEC(2,IVEC),
     &                       RVEC(3,IVEC),TAU,A3(IVEC),B3(IVEC),C3(IVEC)
     &                       ,D3(IVEC))
      END DO
C
C---------------------------------------------------------------
C Find the Voronoi polyhedron.
      NPLANE = NVEC
C
      CALL SFNPOLYHEDRON(NPLANE,A3,B3,C3,D3,NEDGEMAX,NFACE,NEDGE,XEDGE,
     &                   YEDGE,ZEDGE,NVERTEX,NPOIMAX,NNEIMAX,NPLANEMAX,
     &                   NFACEMAX,IFLAG_VERTEX)
C
      IF ( IFLAG_VERTEX.NE.0 ) RETURN
C
C---------------------------------------------------------------
C Calculate the volume as sum of the volumes of all tetrahedra
C connecting the origin to the faces. Use for each tetrahedron
C volume = det((r0-r1),(r0-r2),(r0-r3))/6, where r0 is here the
C origin and r1,r2,r3 the vectors of the 3 other edges.
      VOLUME = 0.D0
C
      DO IFACE = 1,NFACE
         X1 = XEDGE(1,IFACE)
         Y1 = YEDGE(1,IFACE)
         Z1 = ZEDGE(1,IFACE)
C
         DO IEDGE = 2,NEDGE(IFACE) - 1
            X2 = XEDGE(IEDGE,IFACE)
            Y2 = YEDGE(IEDGE,IFACE)
            Z2 = ZEDGE(IEDGE,IFACE)
            X3 = XEDGE(IEDGE+1,IFACE)
            Y3 = YEDGE(IEDGE+1,IFACE)
            Z3 = ZEDGE(IEDGE+1,IFACE)
            TETRVOL = X1*(Y2*Z3-Y3*Z2) + X2*(Y3*Z1-Y1*Z3)
     &                + X3*(Y1*Z2-Y2*Z1)
C
            VOLUME = VOLUME + DABS(TETRVOL)
         END DO
      END DO
      VOLUME = VOLUME/6.D0
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,99002) NFACE
         DO I = 1,NFACE
            WRITE (6,99001) I,NEDGE(I)
         END DO
         WRITE (6,99003) NVERTEX,VOLUME
      END IF
C-----------------------------------------------------------------------
C Find RMT:
      RMT = SFNDISTPLANE(A3(1),B3(1),C3(1),D3(1))
      DO IFACE = 2,NFACE
         TEMP = SFNDISTPLANE(A3(IFACE),B3(IFACE),C3(IFACE),D3(IFACE))
         IF ( TEMP.LT.RMT ) RMT = TEMP
      END DO
C Fint ROUT:
      ROUT = 0.D0
      DO IFACE = 1,NFACE
         DO IEDGE = 1,NEDGE(IFACE)
            TEMP = XEDGE(IEDGE,IFACE)*XEDGE(IEDGE,IFACE)
     &             + YEDGE(IEDGE,IFACE)*YEDGE(IEDGE,IFACE)
     &             + ZEDGE(IEDGE,IFACE)*ZEDGE(IEDGE,IFACE)
            IF ( TEMP.GT.ROUT ) ROUT = TEMP
         END DO
      END DO
      ROUT = DSQRT(ROUT)
99001 FORMAT (10X,'Face ',I4,'   has ',I4,' vertices ')
99002 FORMAT (10X,'Polyhedron properties ',/,10X,'Number of faces :',I5)
99003 FORMAT (10X,'Total number of vertices  :',I10,/,10X,
     &        'The Volume is             :',F10.5)
      END
C*==sfnnormalplane.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNNORMALPLANE(X1,Y1,Z1,X2,Y2,Z2,TAU,A,B,C,D)
C***********************************************************************
C Given two points in space, r1=(X1,Y1,Z1) and r2=(X2,Y2,Z2), this
C subroutine returns the coefficients defining a plane through the
C equation A*x+B*y+C*z=D, which is normal to the vector r2-r1 and passes
C through the point (1.-TAU)*r1 + TAU*r2 (TAU thus being a parameter
C defining how close the plane is to each of the two points).
C***********************************************************************
C
      IMPLICIT NONE
C*--SFNNORMALPLANE228
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,B,C,D,TAU,X1,X2,Y1,Y2,Z1,Z2
C
C Local variables
C
      REAL*8 ONEMTAU
C
C*** End of declarations rewritten by SPAG
C
C The plane is defined as
C (A,B,C)*(X-X1,Y-Y1,Z-Z1)=const=
C                         =(distance from r1 to (1.-TAU)*r1 + TAU*r2)**2
C so A,B,C are the coords. of a vector connecting the point r1 to
C the point (1.-TAU)*r1 + TAU*r2.
      ONEMTAU = 1.D0 - TAU
C
      A = ONEMTAU*X1 + TAU*X2
      B = ONEMTAU*Y1 + TAU*Y2
      C = ONEMTAU*Z1 + TAU*Z2
      D = A*(A+X1) + B*(B+Y1) + C*(C+Z1)
C
      END
C*==sfnpolyhedron.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNPOLYHEDRON(NPLANE,A3,B3,C3,D3,NEDGEMAX,NFACE,NEDGE,
     &                         XEDGE,YEDGE,ZEDGE,NVERTEX,NPOIMAX,
     &                         NNEIMAX,NPLANEMAX,NFACEMAX,IFLAG_VERTEX)
C
C Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
C a convex part of space (the minimal one containing the origin,
C usually a WS-polyhedron), this subroutine returns the actual faces of
C the polyhedron, discarding the planes that do not contain faces. Also,
C the coordinates of the edges of the faces XEDGE,YEDGE,ZEDGE and their
C number NEDGE per face are returned. The coefficients of the actual
C faces are returned in the same arrays A3,B3,C3, and D3.
C
      IMPLICIT NONE
C*--SFNPOLYHEDRON280
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFLAG_VERTEX,NEDGEMAX,NFACE,NFACEMAX,NNEIMAX,NPLANE,
     &        NPLANEMAX,NPOIMAX,NVERTEX
      REAL*8 A3(NFACEMAX),B3(NFACEMAX),C3(NFACEMAX),D3(NFACEMAX),
     &       XEDGE(NEDGEMAX,NFACEMAX),YEDGE(NEDGEMAX,NFACEMAX),
     &       ZEDGE(NEDGEMAX,NFACEMAX)
      INTEGER NEDGE(NFACEMAX)
C
C Local variables
C
      INTEGER ID,IEDGE,INDEXLOW
C
C*** End of declarations rewritten by SPAG
C
C Input:
C                  ! Initial number of planes.
C Input and Output
C                  ! Max. number of edges per plane.
C                  ! dimensioned >= NPLANE.
C Output:
C                  ! Coefs. defining the planes,
C                  ! Number of edges found for each face
C                  ! Cartesian coords. of edges for each plane
C                  ! (2nd index is for planes).
C Inside:
C---------------------------------------------------------------
C Find the faces and edges of the polyhedron.
C
      CALL SFNFINDVERTEX(NPLANE,A3,B3,C3,D3,NEDGEMAX,NFACE,NEDGE,XEDGE,
     &                   YEDGE,ZEDGE,NVERTEX,NPOIMAX,NNEIMAX,NPLANEMAX,
     &                   NFACEMAX,IFLAG_VERTEX)
C
      IF ( IFLAG_VERTEX.NE.0 ) RETURN
C
C---------------------------------------------------------------
C Pack the planes that contain faces at the beginning of the arrays
C A3, B3, C3, D3, and do the same for NEDGE,XEDGE,YEDGE,ZEDGE. The
C order is changed.
C
C You have to fill up the arrays up to NFACE, so...
      INDEXLOW = 1
      DO WHILE ( INDEXLOW.LE.NFACE )
C
         IF ( NEDGE(INDEXLOW).EQ.0 ) THEN
C     promote all planes by one
            DO ID = INDEXLOW + 1,NPLANE
               A3(ID-1) = A3(ID)
               B3(ID-1) = B3(ID)
               C3(ID-1) = C3(ID)
               D3(ID-1) = D3(ID)
               NEDGE(ID-1) = NEDGE(ID)
               DO IEDGE = 1,NEDGE(ID-1)
                  XEDGE(IEDGE,ID-1) = XEDGE(IEDGE,ID)
                  YEDGE(IEDGE,ID-1) = YEDGE(IEDGE,ID)
                  ZEDGE(IEDGE,ID-1) = ZEDGE(IEDGE,ID)
               END DO
            END DO
         ELSE
            INDEXLOW = INDEXLOW + 1
         END IF
C
      END DO
C
      END
C*==sfnfindvertex.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNFINDVERTEX(NPLANE,A3,B3,C3,D3,NEDGEMAX,NFACE,NEDGE,
     &                         XEDGE,YEDGE,ZEDGE,NVERTEX,NPOIMAX,
     &                         NNEIMAX,NPLANEMAX,NFACEMAX,IFLAG_VERTEX)
C Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
C a convex part of space (the minimal one containing the origin,
C usually a WS-polyhedron), this subroutine returns the edges
C of this polyhedron in cartesian coordinates. For the planes that
C are not faces of the polyhedron, a value NEDGE(IPLANE)=0 is returned.
C The total no. of faces found is returned as NFACE.
C
C Updated on 3.9.2001 Algorithm is changed
C Present Algorithm:
C   Define a cube around the  point.
C   Define: POI2PLANE  (Point belongs to planes... )
C           POI2POI    (Point is connected with other points...)
C
C   Now suppose a new canditate plane, first check if all points are
C   on the correct side of space (same side as origin)
C   find how many are in wrong side, and look at all connections
C   between points in wrong side to points in correct side (use POI2POI)
C
C   Define new points and also arrays POI2POI,PO2PLANE, for these points
C   Bookkeeping...
C   If you found 3 new points then accept the plane and the new points
C   a few points may already exist, this is taken care later!
C   goto next plane..
C
C Uses logical function SFNHALFSPACE
C
      IMPLICIT NONE
C*--SFNFINDVERTEX391
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNFINDVERTEX')
C
C Dummy arguments
C
      INTEGER IFLAG_VERTEX,NEDGEMAX,NFACE,NFACEMAX,NNEIMAX,NPLANE,
     &        NPLANEMAX,NPOIMAX,NVERTEX
      REAL*8 A3(NFACEMAX),B3(NFACEMAX),C3(NFACEMAX),D3(NFACEMAX),
     &       XEDGE(NEDGEMAX,NFACEMAX),YEDGE(NEDGEMAX,NFACEMAX),
     &       ZEDGE(NEDGEMAX,NFACEMAX)
      INTEGER NEDGE(NFACEMAX)
C
C Local variables
C
      REAL*8 ALFA,POLYPLANEA3(:),POLYPLANEB3(:),POLYPLANEC3(:),
     &       POLYPLANED3(:),POLYPOINTS(:,:),R,SMALL,V1(3)
      INTEGER I,I1,IEDGE,INDTAB(:,:),IPLANE,ITEST,IV,K,NPOI,NPOLYPLAN,
     &        NPOLYPOI,NTRUEEDGE,ORDER(:),POI2PLANE(:,:),POI2POI(:,:),
     &        POIOFPLANE(:,:)
      LOGICAL LKEEP(:),LTAKEN(:)
C
C*** End of declarations rewritten by SPAG
C
C -------------------------------------------------------
C Input:
C             ! Number of planes. (changed on out)
C             ! Max. number of edges per plane.
C             ! dimensioned >= NPLANE.
C             ! changed on output
C Output:
C             ! Coefs. defining the planes,
C             ! Number of edges found for each face
C             ! Number of faces found (with nedge>0).
C
C             ! Cartesian coords. of edges for each plane
C             ! (2nd index is for planes).
C Inside:
C             ! number of vertices
C             ! Plane indices
C The following are for sorting the edges of each face:
C             ! Edge index
C             ! the polyhedron.
      DATA SMALL/1.D-5/
C
      ALLOCATABLE POI2PLANE,POIOFPLANE,POLYPOINTS,POLYPLANEA3
      ALLOCATABLE POLYPLANEB3,POLYPLANEC3,POLYPLANED3,LKEEP,INDTAB
      ALLOCATABLE ORDER,LTAKEN,POI2POI
C
      ALLOCATE (POI2PLANE(0:NPLANEMAX,NPOIMAX))
      ALLOCATE (POIOFPLANE(NPLANEMAX,NPOIMAX),POLYPOINTS(3,NPOIMAX))
      ALLOCATE (POLYPLANEA3(NPLANEMAX),POLYPLANEB3(NPLANEMAX))
      ALLOCATE (POLYPLANEC3(NPLANEMAX),POLYPLANED3(NPLANEMAX))
      ALLOCATE (LKEEP(NEDGEMAX),INDTAB(NEDGEMAX,NPLANEMAX))
      ALLOCATE (ORDER(NEDGEMAX),LTAKEN(NPOIMAX))
      ALLOCATE (POI2POI(0:NNEIMAX,NPOIMAX))
C---------------------------------------------------------------
C Check & initialize
C
      IF ( NPLANE.LT.4 ) WRITE (6,*) 'EDGE3D: NPLANE was only',NPLANE
C
      DO IPLANE = 1,NPLANE
         NEDGE(IPLANE) = 0
      END DO
      DO IPLANE = 1,NPLANEMAX
         DO I1 = 1,NPOIMAX
            POIOFPLANE(IPLANE,I1) = 0
         END DO
      END DO
C
C Define bounding cube this is the initial polyhedron
C
      ALFA = 30.D0
      CALL SFNDEFCUBE(ALFA,NPOLYPOI,NPOLYPLAN,POLYPOINTS,POI2POI,
     &                POI2PLANE,POLYPLANEA3,POLYPLANEB3,POLYPLANEC3,
     &                POLYPLANED3,NPOIMAX,NNEIMAX,NPLANEMAX)
C
C
C now cut original cube with all planes availiable
C
      DO IPLANE = 1,NPLANE
C
         CALL SFNCUTPLANE(NPOLYPOI,NPOLYPLAN,POLYPOINTS,POI2POI,
     &                    POI2PLANE,POLYPLANEA3,POLYPLANEB3,POLYPLANEC3,
     &                    POLYPLANED3,A3(IPLANE),B3(IPLANE),C3(IPLANE),
     &                    D3(IPLANE),NEDGE,POIOFPLANE,NPOIMAX,NNEIMAX,
     &                    NPLANEMAX)
      END DO
C
C canditate plane a3,b3,c3,d3 is checked if it belongs to polyhedron :
C  1. Is copied in polyplaneA3, polyplaneB3 ...
C     this includes also the original cube planes this is checked
C     later.
C  2. The arrays poi2poi, poi2plane are updated this is used to find
C     all points in a face.
C
C Now all info is availiable, just prepare arrays neaded.
C
      DO I = 1,6
         IF ( NEDGE(I).GT.0 ) WRITE (6,*) 
     &                         ' C A U T I O N : Polyhedron maybe open!'
      END DO
C the first 6 planes are the original cube...
C
      DO IPLANE = 1,NPLANE
         DO IV = 1,NEDGE(IPLANE)
C
            INDTAB(IV,IPLANE) = POIOFPLANE(IPLANE,IV)
            XEDGE(IV,IPLANE) = POLYPOINTS(1,POIOFPLANE(IPLANE,IV))
            YEDGE(IV,IPLANE) = POLYPOINTS(2,POIOFPLANE(IPLANE,IV))
            ZEDGE(IV,IPLANE) = POLYPOINTS(3,POIOFPLANE(IPLANE,IV))
C
         END DO
      END DO
C
C Now go on to find vertices of each face
C
C
      NFACE = 0
      DO IPLANE = 1,NPLANE
         IF ( NEDGE(IPLANE).GE.3 ) THEN
            NFACE = NFACE + 1
C
C     now order using previous info about neighbours
C     flag all atoms exept atoms in plane
C
            DO I = 1,NPOIMAX
               LTAKEN(I) = .TRUE.
            END DO
            DO I = 1,NEDGE(IPLANE)
               LTAKEN(INDTAB(I,IPLANE)) = .FALSE.
               ORDER(I) = 0
            END DO
C
            K = 1
            ORDER(K) = INDTAB(1,IPLANE)
            LTAKEN(ORDER(K)) = .TRUE.
            NPOI = ORDER(K)
            K = K + 1
C
C
            DO WHILE ( K.LE.NEDGE(IPLANE) )
C
               IF ( K.GE.4 ) LTAKEN(ORDER(1)) = .FALSE.
C              allow to return back
C
               I = 0
               DO WHILE ( I.LT.POI2POI(0,NPOI) )
                  I = I + 1
                  ITEST = POI2POI(I,NPOI)
                  IF ( .NOT.LTAKEN(ITEST) ) THEN
                     LTAKEN(ITEST) = .TRUE.
                     ORDER(K) = ITEST
                     I = POI2POI(0,NPOI)
C                    go out of loop
                  END IF
               END DO
               IF ( ORDER(K).EQ.0 ) THEN
                  WRITE (6,*) ' Could not find neighbours to connect '
                  WRITE (6,*) ' Edge not closing  STOPING ...        '
                  WRITE (6,*) ' atom :',ITEST,K,IPLANE
                  IF ( IFLAG_VERTEX.LT.0 )
     &                  CALL STOP_MESSAGE(ROUTINE,'2nd run also failed')
                  IFLAG_VERTEX = 1
                  RETURN
               END IF
C
               NPOI = ORDER(K)
               K = K + 1
            END DO
C
            DO IEDGE = 1,NEDGE(IPLANE)
               IV = ORDER(IEDGE)
               XEDGE(IEDGE,IPLANE) = POLYPOINTS(1,IV)
               YEDGE(IEDGE,IPLANE) = POLYPOINTS(2,IV)
               ZEDGE(IEDGE,IPLANE) = POLYPOINTS(3,IV)
            END DO
         END IF
C------- (NEDGE(IPLANE).GE.3)
C
C Now ordering is done
C
      END DO
C---- IPLANE = 1,NPLANE
C
C Now reject vertices that occur more than once
C
      DO IPLANE = 1,NPLANE
         DO IEDGE = 1,NEDGE(IPLANE)
            LKEEP(IEDGE) = .TRUE.
         END DO
C
         DO IEDGE = 1,NEDGE(IPLANE)
            DO ITEST = IEDGE + 1,NEDGE(IPLANE)
               V1(1) = XEDGE(ITEST,IPLANE) - XEDGE(IEDGE,IPLANE)
               V1(2) = YEDGE(ITEST,IPLANE) - YEDGE(IEDGE,IPLANE)
               V1(3) = ZEDGE(ITEST,IPLANE) - ZEDGE(IEDGE,IPLANE)
               R = SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
               IF ( R.LT.SMALL ) THEN
                  LKEEP(ITEST) = .FALSE.
                  WRITE (6,*) 'Vertex ',ITEST,' in plane ',IPLANE,
     &                        ' found again'
               END IF
            END DO
         END DO
C    now get rid of double points
         NTRUEEDGE = 0
         DO IEDGE = 1,NEDGE(IPLANE)
            IF ( LKEEP(IEDGE) ) NTRUEEDGE = NTRUEEDGE + 1
         END DO
         IF ( NTRUEEDGE.NE.NEDGE(IPLANE) ) THEN
C    next plane
C
            IEDGE = 1
            DO WHILE ( IEDGE.LE.NTRUEEDGE )
C
               IF ( .NOT.LKEEP(IEDGE) ) THEN
C    promote by one all verices
                  DO ITEST = IEDGE + 1,NEDGE(IPLANE)
                     XEDGE(ITEST-1,IPLANE) = XEDGE(ITEST,IPLANE)
                     YEDGE(ITEST-1,IPLANE) = YEDGE(ITEST,IPLANE)
                     ZEDGE(ITEST-1,IPLANE) = ZEDGE(ITEST,IPLANE)
                     INDTAB(ITEST-1,IPLANE) = INDTAB(ITEST,IPLANE)
                     LKEEP(ITEST-1) = LKEEP(ITEST)
                  END DO
               ELSE
                  IEDGE = IEDGE + 1
               END IF
            END DO
C
            NEDGE(IPLANE) = NTRUEEDGE
         END IF
C
      END DO
C
C
C Now update the plane coordinates
C
      NPLANE = NPOLYPLAN
      DO IPLANE = 1,NPLANE
         A3(IPLANE) = POLYPLANEA3(IPLANE)
         B3(IPLANE) = POLYPLANEB3(IPLANE)
         C3(IPLANE) = POLYPLANEC3(IPLANE)
         D3(IPLANE) = POLYPLANED3(IPLANE)
      END DO
      NVERTEX = NPOLYPOI
C
      IFLAG_VERTEX = 0
C
      DEALLOCATE (POI2PLANE,POIOFPLANE,POLYPOINTS,POLYPLANEA3)
      DEALLOCATE (POLYPLANEB3,POLYPLANEC3,POLYPLANED3,LKEEP,INDTAB)
      DEALLOCATE (ORDER,LTAKEN,POI2POI)
C
      END
C*==sfndefcube.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNDEFCUBE(ALFA,NPOINTS,NPLANES,POLYPOINTS,POI2POI,
     &                      POI2PLANE,POLYPLANEA3,POLYPLANEB3,
     &                      POLYPLANEC3,POLYPLANED3,NPOIMAX,NNEIMAX,
     &                      NPLANEMAX)
      IMPLICIT NONE
C*--SFNDEFCUBE670
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALFA
      INTEGER NNEIMAX,NPLANEMAX,NPLANES,NPOIMAX,NPOINTS
      INTEGER POI2PLANE(0:NPLANEMAX,NPOIMAX),POI2POI(0:NNEIMAX,NPOIMAX)
      REAL*8 POLYPLANEA3(NPLANEMAX),POLYPLANEB3(NPLANEMAX),
     &       POLYPLANEC3(NPLANEMAX),POLYPLANED3(NPLANEMAX),
     &       POLYPOINTS(3,NPOIMAX)
C
C Local variables
C
      REAL*8 AO2
      INTEGER I,IPLANE,IPOINT1
C
C*** End of declarations rewritten by SPAG
C
C This just defines all arrays for a cube of edge alfa
C
      DO IPOINT1 = 0,NNEIMAX
         DO I = 1,NPOIMAX
            POI2POI(IPOINT1,I) = 0
         END DO
      END DO
      DO IPOINT1 = 0,NPLANEMAX
         DO I = 1,NPOIMAX
            POI2PLANE(IPOINT1,I) = 0
         END DO
      END DO
C
      DO IPOINT1 = 1,NPOIMAX
         DO I = 1,3
            POLYPOINTS(I,IPOINT1) = 0.D0
         END DO
      END DO
      NPOINTS = 8
      AO2 = ALFA/2.D0
      POLYPOINTS(1,1) = AO2
      POLYPOINTS(2,1) = AO2
      POLYPOINTS(3,1) = AO2
C
      POLYPOINTS(1,2) = AO2
      POLYPOINTS(2,2) = -AO2
      POLYPOINTS(3,2) = AO2
C
      POLYPOINTS(1,3) = -AO2
      POLYPOINTS(2,3) = -AO2
      POLYPOINTS(3,3) = AO2
C
      POLYPOINTS(1,4) = -AO2
      POLYPOINTS(2,4) = AO2
      POLYPOINTS(3,4) = AO2
C
      POLYPOINTS(1,5) = AO2
      POLYPOINTS(2,5) = AO2
      POLYPOINTS(3,5) = -AO2
C
      POLYPOINTS(1,6) = AO2
      POLYPOINTS(2,6) = -AO2
      POLYPOINTS(3,6) = -AO2
C
      POLYPOINTS(1,7) = -AO2
      POLYPOINTS(2,7) = -AO2
      POLYPOINTS(3,7) = -AO2
C
      POLYPOINTS(1,8) = -AO2
      POLYPOINTS(2,8) = AO2
      POLYPOINTS(3,8) = -AO2
C ************************************************
      POI2POI(0,1) = 3 ! number of neighbors of point 1
      POI2POI(1,1) = 2 ! point indeces
      POI2POI(2,1) = 4
      POI2POI(3,1) = 5
      POI2PLANE(0,1) = 3  ! number of planes passing from point 1
      POI2PLANE(1,1) = 1  ! plane indeces
      POI2PLANE(2,1) = 3
      POI2PLANE(3,1) = 5
C
      POI2POI(0,2) = 3 ! number of neighbors of point
      POI2POI(1,2) = 1 ! point indeces
      POI2POI(2,2) = 3
      POI2POI(3,2) = 6
      POI2PLANE(0,2) = 3  ! number of planes passing from point
      POI2PLANE(1,2) = 1  ! plane indeces
      POI2PLANE(2,2) = 4
      POI2PLANE(3,2) = 5
C
      POI2POI(0,3) = 3 ! number of neighbors of point
      POI2POI(1,3) = 2 ! point indeces
      POI2POI(2,3) = 4
      POI2POI(3,3) = 7
      POI2PLANE(0,3) = 3  ! number of planes passing from point
      POI2PLANE(1,3) = 2  ! plane indeces
      POI2PLANE(2,3) = 4
      POI2PLANE(3,3) = 5
C
      POI2POI(0,4) = 3 ! number of neighbors of point
      POI2POI(1,4) = 1 ! point indeces
      POI2POI(2,4) = 3
      POI2POI(3,4) = 8
      POI2PLANE(0,4) = 3  ! number of planes passing from point
      POI2PLANE(1,4) = 2  ! plane indeces
      POI2PLANE(2,4) = 3
      POI2PLANE(3,4) = 5
C
      POI2POI(0,5) = 3 ! number of neighbors of point
      POI2POI(1,5) = 1 ! point indeces
      POI2POI(2,5) = 6
      POI2POI(3,5) = 8
      POI2PLANE(0,5) = 3  ! number of planes passing from point
      POI2PLANE(1,5) = 1  ! plane indeces
      POI2PLANE(2,5) = 3
      POI2PLANE(3,5) = 6
C
      POI2POI(0,6) = 3 ! number of neighbors of point
      POI2POI(1,6) = 2 ! point indeces
      POI2POI(2,6) = 5
      POI2POI(3,6) = 7
      POI2PLANE(0,6) = 3  ! number of planes passing from point
      POI2PLANE(1,6) = 1  ! plane indeces
      POI2PLANE(2,6) = 4
      POI2PLANE(3,6) = 6
C
      POI2POI(0,7) = 3 ! number of neighbors of point
      POI2POI(1,7) = 3 ! point indeces
      POI2POI(2,7) = 6
      POI2POI(3,7) = 8
      POI2PLANE(0,7) = 3  ! number of planes passing from point
      POI2PLANE(1,7) = 2  ! plane indeces
      POI2PLANE(2,7) = 4
      POI2PLANE(3,7) = 6
C
      POI2POI(0,8) = 3 ! number of neighbors of point
      POI2POI(1,8) = 4 ! point indeces
      POI2POI(2,8) = 5
      POI2POI(3,8) = 7
      POI2PLANE(0,8) = 3  ! number of planes passing from point
      POI2PLANE(1,8) = 2  ! plane indeces
      POI2PLANE(2,8) = 3
      POI2PLANE(3,8) = 6
C ---------------------------
      NPLANES = 6
      DO IPLANE = 1,NPLANES
         POLYPLANEA3(IPLANE) = 0.D0
         POLYPLANEB3(IPLANE) = 0.D0
         POLYPLANEC3(IPLANE) = 0.D0
         POLYPLANED3(IPLANE) = 0.D0
      END DO
      POLYPLANEA3(1) = 1.D0
      POLYPLANED3(1) = ALFA/2.D0
C
      POLYPLANEA3(2) = 1.D0
      POLYPLANED3(2) = -ALFA/2.D0
C
      POLYPLANEB3(3) = 1.D0
      POLYPLANED3(3) = ALFA/2.D0
C
      POLYPLANEB3(4) = 1.D0
      POLYPLANED3(4) = -ALFA/2.D0
C
      POLYPLANEC3(5) = 1.D0
      POLYPLANED3(5) = ALFA/2.D0
C
      POLYPLANEC3(6) = 1.D0
      POLYPLANED3(6) = -ALFA/2.D0
C
      END
C*==sfncutplane.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNCUTPLANE(NPOLYPOI,NPOLYPLAN,POLYPOINTS,POI2POI,
     &                       POI2PLANE,POLYPLANEA3,POLYPLANEB3,
     &                       POLYPLANEC3,POLYPLANED3,A3,B3,C3,D3,NEDGE,
     &                       POIOFPLANE,NPOIMAX,NNEIMAX,NPLANEMAX)
C
C * This sub takes a polyhedron and the equation of a plane
C   A3*x+B3*y+C3*z=D3       and produces  new polyhedron cuted
C   by this plane.
C
      IMPLICIT NONE
C*--SFNCUTPLANE862
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNCUTPLANE')
C
C Dummy arguments
C
      REAL*8 A3,B3,C3,D3
      INTEGER NNEIMAX,NPLANEMAX,NPOIMAX,NPOLYPLAN,NPOLYPOI
      INTEGER NEDGE(*),POI2PLANE(0:NPLANEMAX,NPOIMAX),
     &        POI2POI(0:NNEIMAX,NPOIMAX),POIOFPLANE(NPLANEMAX,NPOIMAX)
      REAL*8 POLYPLANEA3(NPLANEMAX),POLYPLANEB3(NPLANEMAX),
     &       POLYPLANEC3(NPLANEMAX),POLYPLANED3(NPLANEMAX),
     &       POLYPOINTS(3,NPOIMAX)
C
C Local variables
C
      REAL*8 A1,NEWPOINT(:,:),X1,X2,XCUT,Y1,Y2,YCUT,Z1,Z2,ZCUT
      INTEGER I,I1,IN,IP,IP0,IP1,IP2,ISIMILAR(:),ISIMILAR1(:),
     &        LOCALPLANES(:),NEWPLANES(:,:),NNEW,NUM
      LOGICAL KEEPIT(:)
      LOGICAL SFNHALFSPACE
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE ISIMILAR,NEWPOINT,ISIMILAR1,NEWPLANES,LOCALPLANES
      ALLOCATABLE KEEPIT
C
      ALLOCATE (ISIMILAR(NNEIMAX),NEWPOINT(3,NNEIMAX))
      ALLOCATE (ISIMILAR1(NNEIMAX),NEWPLANES(0:NNEIMAX,NPOIMAX))
      ALLOCATE (LOCALPLANES(0:NPLANEMAX),KEEPIT(NPOIMAX))
C
      DO IP = 1,NPOLYPOI
         XCUT = POLYPOINTS(1,IP)
         YCUT = POLYPOINTS(2,IP)
         ZCUT = POLYPOINTS(3,IP)
         KEEPIT(IP) = (SFNHALFSPACE(A3,B3,C3,D3,XCUT,YCUT,ZCUT))
      END DO
      NNEW = 0
C
      DO I = 1,NNEIMAX
         ISIMILAR(I) = 0
         ISIMILAR1(I) = 0
      END DO
C
      DO IP = 1,NPOLYPOI
         IF ( .NOT.KEEPIT(IP) ) THEN
C     one point is left outside, have to redefine polyhedron
C     Introduce new points
C
            X1 = POLYPOINTS(1,IP)
            Y1 = POLYPOINTS(2,IP)
            Z1 = POLYPOINTS(3,IP)
C
C     loop in all neighbours of ip
            DO IP0 = 1,POI2POI(0,IP)
               IP2 = POI2POI(IP0,IP)
C
C     Now find crossing of line connecting point ip with all neighbors ip2
C     and plane a3,b3,c3,d3
C
               IF ( KEEPIT(IP2) ) THEN
                  X2 = POLYPOINTS(1,IP2)
                  Y2 = POLYPOINTS(2,IP2)
                  Z2 = POLYPOINTS(3,IP2)
C
                  CALL SFNCROSSPOIPLANE(X1,Y1,Z1,X2,Y2,Z2,A3,B3,C3,D3,
     &                                  XCUT,YCUT,ZCUT,A1)
C     new point introduced, keep it until all points are found
                  NNEW = NNEW + 1
                  NEWPOINT(1,NNEW) = XCUT
                  NEWPOINT(2,NNEW) = YCUT
                  NEWPOINT(3,NNEW) = ZCUT
C
                  IF ( ABS(A1).LT.1.D-8 ) ISIMILAR(NNEW) = IP
                  IF ( ABS(A1-1.D0).LT.1.D-8 ) ISIMILAR1(NNEW) = IP2
C     find planes that go through ip,ip2 and keep
                  CALL SFNCOMMONPLANE(IP,IP2,LOCALPLANES,POI2PLANE,
     &                                NPOIMAX,NNEIMAX,NPLANEMAX)
C
                  NEWPLANES(0,NNEW) = LOCALPLANES(0)
                  DO I = 1,LOCALPLANES(0)
                     NEWPLANES(I,NNEW) = LOCALPLANES(I)
                  END DO
C
               END IF
            END DO
C   --------- all neigbours of excluded point are done now decide if you
C     will keep the plane, and if you will keep the excluded point...
C
C
         END IF
C------  if (.NOT.keepit(ip)) THEN
      END DO
C---- ip=1,npolypoi
C
C New points are now stored see if we have a cut in the polyhedron
C
      IF ( NNEW.GT.2 ) THEN
C     add plane in the list of planes and update points
         NPOLYPLAN = NPOLYPLAN + 1
         IF ( NPOLYPLAN.GT.NPLANEMAX )
     &         CALL STOP_MESSAGE(ROUTINE,'increase NPLANEMAX')
C
         POLYPLANEA3(NPOLYPLAN) = A3
         POLYPLANEB3(NPOLYPLAN) = B3
         POLYPLANEC3(NPOLYPLAN) = C3
         POLYPLANED3(NPOLYPLAN) = D3
C
         DO IN = 1,NNEW
C
            IP1 = ISIMILAR(IN)
            IP2 = ISIMILAR1(IN)
C if all ok they can never be both zero
            IP = 0
            IF ( IP1.NE.0 ) IP = IP1
            IF ( IP2.NE.0 ) IP = IP2
            IF ( IP1.NE.0 .AND. IP2.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'check IP1, IP2')
            IF ( IP.NE.0 ) THEN
C     point already exists add plane in list of planes and find
C     neighbors in this plane --  Now updating existing point
               POI2PLANE(0,IP) = POI2PLANE(0,IP) + 1
C
               NUM = POI2PLANE(0,IP)
               POI2PLANE(NUM,IP) = NPOLYPLAN
               KEEPIT(IP) = .TRUE.
            ELSE
C     point is new, just append it to the list
C
               NPOLYPOI = NPOLYPOI + 1
               IF ( NPOLYPOI.GT.NPOIMAX )
     &               CALL STOP_MESSAGE(ROUTINE,'increase NPOIMAX')
C
               DO I = 1,3
                  POLYPOINTS(I,NPOLYPOI) = NEWPOINT(I,IN)
               END DO
C
               POI2PLANE(0,NPOLYPOI) = NEWPLANES(0,IN) + 1
               POI2PLANE(1,NPOLYPOI) = NPOLYPLAN
               DO I1 = 2,POI2PLANE(0,NPOLYPOI)
                  POI2PLANE(I1,NPOLYPOI) = NEWPLANES(I1-1,IN)
               END DO
C
               KEEPIT(NPOLYPOI) = .TRUE.
C
C
            END IF
C
C     Now the plane list is updated, nead to update
C     the poi2poi array (neighbors list)
C
         END DO
C
C  Points and Planes updated -- now erase points outside polyhedron
C
         CALL SFNERASEPOINTS(NPOLYPOI,KEEPIT,POI2PLANE,POI2POI,
     &                       POLYPOINTS,NNEIMAX,NPOIMAX,NPLANEMAX)
C
         CALL SFNPOI2POIUPDATE(POI2POI,POI2PLANE,NPOLYPLAN,NPOLYPOI,
     &                         NEDGE,POIOFPLANE,NPOIMAX,NNEIMAX,
     &                         NPLANEMAX)
C
      END IF
C
      DEALLOCATE (ISIMILAR,NEWPOINT,ISIMILAR1,NEWPLANES,LOCALPLANES)
      DEALLOCATE (KEEPIT)
C
      END
C*==sfncrosspoiplane.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNCROSSPOIPLANE(X1,Y1,Z1,X2,Y2,Z2,A3,B3,C3,D3,XCUT,
     &                            YCUT,ZCUT,A)
C
C      R = Ro + a V      (eq of line )
C      -   -      -
C
C      C.R = D           (eq of plane)
C      - -
C
C
C  Cross point if we find  a
C
C       a = (D - C.Ro ) / C.V
C                - -      - -
C
C   Ro = (x1,y1,z1), V = (x2-x1,y2-y1,z2-z1), C = (A3,B3,C3), D = D3
C
      IMPLICIT NONE
C*--SFNCROSSPOIPLANE1069
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,A3,B3,C3,D3,X1,X2,XCUT,Y1,Y2,YCUT,Z1,Z2,ZCUT
C
C*** End of declarations rewritten by SPAG
C
      A = (D3-(A3*X1+B3*Y1+C3*Z1))/(A3*(X2-X1)+B3*(Y2-Y1)+C3*(Z2-Z1))
      XCUT = X1 + A*(X2-X1)
      YCUT = Y1 + A*(Y2-Y1)
      ZCUT = Z1 + A*(Z2-Z1)
      END
C*==sfncommonplane.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNCOMMONPLANE(IP1,IP2,LOCALPLANES,POI2PLANE,NPOIMAX,
     &                          NNEIMAX,NPLANEMAX)
C
C This looks at the line connecting 2 points and returns the index
C of common planes
C
      IMPLICIT NONE
C*--SFNCOMMONPLANE1100
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IP1,IP2,NNEIMAX,NPLANEMAX,NPOIMAX
      INTEGER LOCALPLANES(0:NPLANEMAX),POI2PLANE(0:NPLANEMAX,NPOIMAX)
C
C Local variables
C
      INTEGER I,I1,I2,IC,INDTAB(NPLANEMAX),NL
C
C*** End of declarations rewritten by SPAG
C
      IF ( NNEIMAX.LE.0 ) WRITE (6,*) 'DUMMY',NNEIMAX
      NL = 0
      DO I1 = 1,POI2PLANE(0,IP1)
C
         DO I2 = 1,POI2PLANE(0,IP2)
C
            IC = POI2PLANE(I1,IP1) - POI2PLANE(I2,IP2)
C
            IF ( IC.EQ.0 ) THEN
C     common plane found
               NL = NL + 1
               INDTAB(NL) = POI2PLANE(I1,IP1)
            END IF
         END DO
      END DO
C
      LOCALPLANES(0) = NL
C
      DO I = 1,NL
         LOCALPLANES(I) = INDTAB(I)
      END DO
      END
C*==sfnerasepoints.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNERASEPOINTS(NPOLYPOI,KEEPIT,POI2PLANE,POI2POI,
     &                          POLYPOINTS,NNEIMAX,NPOIMAX,NPLANEMAX)
      IMPLICIT NONE
C*--SFNERASEPOINTS1152
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NNEIMAX,NPLANEMAX,NPOIMAX,NPOLYPOI
      LOGICAL KEEPIT(NPOIMAX)
      INTEGER POI2PLANE(0:NPLANEMAX,NPOIMAX),POI2POI(0:NNEIMAX,NPOIMAX)
      REAL*8 POLYPOINTS(3,NPOIMAX)
C
C Local variables
C
      INTEGER I,I0,IP,N,T_POI2PLANE(:),T_POI2POI(:)
      REAL*8 TR(3)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE T_POI2POI,T_POI2PLANE
      ALLOCATE (T_POI2POI(0:NNEIMAX),T_POI2PLANE(0:NPLANEMAX))
C
      N = 0
      I0 = 0
C
      DO IP = 1,NPOLYPOI
C
         IF ( .NOT.KEEPIT(IP) ) THEN
C
            I0 = I0 + 1
C           ! !!!!!!!! erase(i0) = ip
         ELSE
C
            N = N + 1
C     copy to tmp
            T_POI2POI(0) = POI2POI(0,IP)
C
            DO I = 1,T_POI2POI(0)
               T_POI2POI(I) = POI2POI(I,IP)
            END DO
            T_POI2PLANE(0) = POI2PLANE(0,IP)
            DO I = 1,T_POI2PLANE(0)
               T_POI2PLANE(I) = POI2PLANE(I,IP)
            END DO
            DO I = 1,3
               TR(I) = POLYPOINTS(I,IP)
            END DO
C     now map back
C
            POI2POI(0,N) = T_POI2POI(0)
            DO I = 1,T_POI2POI(0)
               POI2POI(I,N) = T_POI2POI(I)
            END DO
            POI2PLANE(0,N) = T_POI2PLANE(0)
            DO I = 1,T_POI2PLANE(0)
               POI2PLANE(I,N) = T_POI2PLANE(I)
C
            END DO
            DO I = 1,3
               POLYPOINTS(I,N) = TR(I)
            END DO
C
         END IF
C
      END DO
C  now we know which have to be errased...
      NPOLYPOI = N
C
      DEALLOCATE (T_POI2POI,T_POI2PLANE)
C
      END
C*==sfnpoi2poiupdate.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SFNPOI2POIUPDATE(POI2POI,POI2PLANE,NPOLYPLAN,NPOLYPOI,
     &                            NEDGE,POIOFPLANE,NPOIMAX,NNEIMAX,
     &                            NPLANEMAX)
C
C    This sub updates the neighbours array poi2poi, it checks all atoms
C    some are newly defined, and looks if a pair of atoms belong
C    to the same 2 planes
C
      IMPLICIT NONE
C*--SFNPOI2POIUPDATE1243
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NNEIMAX,NPLANEMAX,NPOIMAX,NPOLYPLAN,NPOLYPOI
      INTEGER NEDGE(*),POI2PLANE(0:NPLANEMAX,NPOIMAX),
     &        POI2POI(0:NNEIMAX,NPOIMAX),POIOFPLANE(NPLANEMAX,NPOIMAX)
C
C Local variables
C
      INTEGER IND,IP1,IP2,IPLANE,LOCALPLANES(:),NNEI,NP,POINTS(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE LOCALPLANES,POINTS
      ALLOCATE (LOCALPLANES(0:NPLANEMAX),POINTS(NPOIMAX))
C
      DO IP1 = 1,NPOLYPOI
         NNEI = 0
         DO IP2 = 1,NPOLYPOI
            IF ( IP1.NE.IP2 ) THEN
C
               CALL SFNCOMMONPLANE(IP1,IP2,LOCALPLANES,POI2PLANE,
     &                             NPOIMAX,NNEIMAX,NPLANEMAX)
C
               IF ( LOCALPLANES(0).EQ.2 ) THEN
C     this pair are neighbours, add it to the list of poi2poi
                  NNEI = NNEI + 1
                  POI2POI(NNEI,IP1) = IP2
               END IF
            END IF
         END DO
         POI2POI(0,IP1) = NNEI
      END DO
C
C    Poi2poi updated now look if all planes have 3 or more atoms in them
C
      DO IP1 = 1,NPOIMAX
         POINTS(IP1) = 0
      END DO
      DO IP1 = 1,NPOLYPOI
         NP = POI2PLANE(0,IP1)
         DO IPLANE = 1,NP
            IND = POI2PLANE(IPLANE,IP1)
            IF ( IND.GT.0 ) THEN
               POINTS(IND) = POINTS(IND) + 1
               POIOFPLANE(IND,POINTS(IND)) = IP1
            END IF
         END DO
      END DO
C
      NP = 0
      DO IPLANE = 1,NPOLYPLAN
         NEDGE(IPLANE) = POINTS(IPLANE)
      END DO
C
      DEALLOCATE (LOCALPLANES,POINTS)
C
      END
C*==sfndistplane.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
C***********************************************************************
      REAL*8 FUNCTION SFNDISTPLANE(A,B,C,D)
C Returns the distance of a plane A*x+B*y+C*z=D to the origin.
      IMPLICIT NONE
C*--SFNDISTPLANE1320
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNDISTPLANE')
C
C Dummy arguments
C
      REAL*8 A,B,C,D
C
C Local variables
C
      REAL*8 ABCSQ
C
C*** End of declarations rewritten by SPAG
C
      ABCSQ = A*A + B*B + C*C
C
      IF ( ABCSQ.LT.1.D-100 )
     &      CALL STOP_MESSAGE(ROUTINE,'ABCSQ < 1.D-100')
C
      SFNDISTPLANE = DABS(D)/DSQRT(ABCSQ)
C
      END
C*==sfnhalfspace.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      LOGICAL FUNCTION SFNHALFSPACE(A,B,C,D,X,Y,Z)
C   ********************************************************************
C   *                                                                  *
C   * Given a plane A*x+B*y+C*z=D, and a point (X,Y,Z) in space, this  *
C   * function takes the value TRUE if (X,Y,Z) lies in the half-space  *
C   * defined by the plane and the origin (0,0,0) (including the plane *
C   * itself). Else, the value FALSE is returned.                      *
C   *                                                                  *
C   * criterion used:  the inner product of the vector (X,Y,Z)  with   *
C   * the vector d connecting the origin to the plane vertically be    *
C   * less than or equal to d**2:  (d_x,d_y,d_z)*(X,Y,Z) =< d**2.      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--SFNHALFSPACE1377
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNHALFSPACE')
C
C Dummy arguments
C
      REAL*8 A,B,C,D,X,Y,Z
C
C*** End of declarations rewritten by SPAG
C
      IF ( DABS(A)+DABS(B)+DABS(C).LT.1.D-80 )
     &     CALL STOP_MESSAGE(ROUTINE,'A,B,C too small')
C
      SFNHALFSPACE = .FALSE.
C
      IF ( D*(A*X+B*Y+C*Z).LE.D*D+1.D-7 ) SFNHALFSPACE = .TRUE.
C
      END
