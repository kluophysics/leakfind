C*==sfnshape.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNSHAPE(NRSF0,AFACE8,BFACE8,CFACE8,DFACE8,NVERTICES8,
     &                    XEDGE8,YEDGE8,ZEDGE8,NFACE,LMAX,DLT,KEYPAN,
     &                    NM8,NCELL_S,SCALE_S,NPAN_S,MESHN_S,NM_S,XRN_S,
     &                    DRN_S,NFUN_S,LMIFUN_S,FLMSF_S,NPANMAX,NRSFMAX,
     &                    LSFMAX,NLMSFMAX,NFACEMAX,NVERTMAX,
     &                    IFLAG_POLYHEDRON)
C   ********************************************************************
C   *                                                                  *
C   *     set of subroutines based on                                  *
C   *                                                                  *
C   *                    S H A P E   P R O G R A M                     *
C   *     F O R  A R B I T R A R Y  V O R O N O I  P O L Y H E D R A   *
C   *                                                                  *
C   *                                                by   N.Stefanou   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   this file contains the shape program as well as all            *
C   *   auxilary subroutines and functions. the following              *
C   *   modifications were made:                                       *
C   *                                                                  *
C   *   - all routines renamed to  SFN...                              *
C   *   - all COMMON blocks removed                                    *
C   *   - all include statements removed                               *
C   *   - all PARAMETER statements removed exept in  SFNSHAPE          *
C   *     the argument lists have been extended accordingly            *
C   *   - all temporary arrays are accessed using   ALLOCATE           *
C   *                                                                  *
C   *   10/05/14 output to scratch files suppressed                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   this program calculates the angular momentum components of the *
C   *   shape function for an arbitrary voronoi polyhedron.            *
C   *   a real spherical harmonic basis is used for the decomposition. *
C   *   on input we give :                                             *
C   *                                                                  *
C   *      LMAX          :  maximum angular momentum                   *
C   *                                                                  *
C   *      DLT           :  defines the step for Gauss-Legendre calc.  *
C   *                       to be at the save side use 0.05 or 0.1     *
C   *                                                                  *
C   *      NFACE         :  number of faces of the polyhedron          *
C   *      KEYPAN        :  key to define  the  radial  mesh.          *
C   *                       if KEYPAN = X  the default radial division *
C   *                       of pannels given in data statement is used.*
C   *                    ** IN THIS VERSION THE MESH IS DETERMINED     *
C   *                       BY SUBROUTINE MESH0.                       *
C   *      Z(I)          :  coefficients of the equation of a face     *
C   *                       Z(1)*X + Z(2)*Y + Z(3)*Z  =  1             *
C   *      NVERT         :  number of vertices of a face               *
C   *      V(I,IVERT)    :  coordinates of the vertices of a face      *
C   *      NEWSCH(IFACE) :  integer   parameter to calculate   (=1)    *
C   *                       the contribution of the corresponding      *
C   *                       pyramid to the shape functions  .  if      *
C   *                       NEWSCH.NE.1 the contribution is taken      *
C   *                       equal to that of the previous pyramid      *
C   *                                                                  *
C   *                                                                  *
C   *   the definition of real spherical harmonics is not the standard *
C   *   one refered in the paper:                                      *
C   *   N. Stefanou, H. Akai and R. Zeller,                            *
C   *   Computer Phys. Commun. 60 (1990) 231                           *
C   *   if you want to have angular momentum components in the         *
C   *   standard basis change the following statements in the routines:*
C   *   IN CCOEF      ISI=1                        ->  ISI=1-2*MOD(M,  *
C   *   IN DROTREAL   IF(MOD(M+MP),2).EQ.0) D=-D   ->  DELETE THE LIN  *
C   *                                                                  *
C   *                                                                  *
C   * The equation        A1*x + A2*y + A3*z = A4                      *
C   * is transformed to   A1*x'+ A2*y'+ A3*z'= A4 -(A1*dx+A2*dy+A3dz)  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNSHAPE')
      INTEGER ICD,NDIM
      PARAMETER (ICD=13000,NDIM=20000)
C
C Dummy arguments
C
      REAL*8 DLT,SCALE_S
      INTEGER IFLAG_POLYHEDRON,KEYPAN,LMAX,LSFMAX,MESHN_S,NCELL_S,NFACE,
     &        NFACEMAX,NFUN_S,NLMSFMAX,NPANMAX,NPAN_S,NRSF0,NRSFMAX,
     &        NVERTMAX
      REAL*8 AFACE8(NFACEMAX),BFACE8(NFACEMAX),CFACE8(NFACEMAX),
     &       DFACE8(NFACEMAX),DRN_S(NRSFMAX),FLMSF_S(NRSFMAX,NLMSFMAX),
     &       XEDGE8(NVERTMAX,NFACEMAX),XRN_S(NRSFMAX),
     &       YEDGE8(NVERTMAX,NFACEMAX),ZEDGE8(NVERTMAX,NFACEMAX)
      INTEGER LMIFUN_S(NLMSFMAX),NM8(NPANMAX),NM_S(NPANMAX),
     &        NVERTICES8(NFACEMAX)
C
C Local variables
C
      REAL*8 A1,A2,A3,A4,AFACE(:),ALPHA(:),ARG1,ARG2,B(:),BETA(:),
     &       BFACE(:),C(:),CFACE(:),CL(:),CRT(:),DFACE(:),DMATL(:,:),
     &       DRN(:),FA(:),FB(:),FD(:),FK,FL,FLMSF_LOC(:,:),GAMMA(:),R,
     &       R0(:),RAP,RD(:),RDOWN,RSUM(:,:),RUPSQ(:),S(:,:),S1(:,:),
     &       S2(:,:),S3(:,:),SCL,V(:,:),XEDGE(:,:),XRN(:),YEDGE(:,:),
     &       Z(3),ZEDGE(:,:)
      INTEGER I,I1,IB,IBM,IC,ICE,ICED,ICEDIM,ICOUNT,IFACE,IMAX,IP,IPAN,
     &        IPMAX,IS,IS0,ISIGNU(:),ISU,ISUM,ISUMD,ISW(:),ITESTSFN,
     &        ITET,IV,IVERT,IVTOT,K,K0,L,LM0,LMA2D,LOFM(:),M,MESHN,
     &        MESHNMAX,MO,MOFM(:),MP,N,NCELL,NEWSCH(:),NFUN,NM(:),NPAN,
     &        NTET,NTT(:),NVERT,NVERTICES(:),NVTOT,NVTOTD
      EXTERNAL SFNCCOEF,SFNCRIT,SFNDROTREAL,SFNMESH,SFNPINTG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE B,C,S,V,R0,S1,S2,S3,NVERTICES,FA,FB,FD,CL,RD,NM
      ALLOCATABLE DRN,CRT,ISW,RSUM,NTT,XRN,BETA,LOFM,MOFM,AFACE,BFACE
      ALLOCATABLE CFACE,DFACE,GAMMA,ALPHA,XEDGE,YEDGE,ZEDGE,DMATL
      ALLOCATABLE RUPSQ,NEWSCH,ISIGNU,FLMSF_LOC
C
      ICED = ((LSFMAX+1)*(LSFMAX+2))/2
      MESHNMAX = NRSFMAX
      NVTOTD = NFACEMAX*NVERTMAX
      IF ( MESHNMAX.EQ.0 ) CALL STOP_MESSAGE(ROUTINE,'MESHNMAX')
C
      ALLOCATE (B(NLMSFMAX),C(ICED),S(-LSFMAX:LSFMAX,0:LSFMAX))
      ALLOCATE (V(3,NVERTMAX),R0(NFACEMAX),S1(-LSFMAX:LSFMAX,0:LSFMAX))
      ALLOCATE (S2(-LSFMAX:LSFMAX,0:LSFMAX),NTT(NFACEMAX))
      ALLOCATE (S3(-LSFMAX:LSFMAX,0:LSFMAX),NVERTICES(NFACEMAX))
      ALLOCATE (FA(NVTOTD),FB(NVTOTD),FD(NVTOTD),CL(ICD))
      ALLOCATE (RD(NVTOTD),NM(NPANMAX),DRN(MESHNMAX),CRT(NPANMAX))
      ALLOCATE (ISW(NLMSFMAX),RSUM(0:LSFMAX,2),XRN(MESHNMAX))
      ALLOCATE (BETA(NFACEMAX),LOFM(NLMSFMAX),MOFM(NLMSFMAX))
      ALLOCATE (AFACE(NFACEMAX),BFACE(NFACEMAX),CFACE(NFACEMAX))
      ALLOCATE (DFACE(NFACEMAX),GAMMA(NFACEMAX),ALPHA(NFACEMAX))
      ALLOCATE (XEDGE(NVERTMAX,NFACEMAX),YEDGE(NVERTMAX,NFACEMAX))
      ALLOCATE (ZEDGE(NVERTMAX,NFACEMAX),RUPSQ(NVTOTD))
      ALLOCATE (NEWSCH(NFACEMAX),ISIGNU(NVTOTD))
      ALLOCATE (FLMSF_LOC(NRSFMAX,NLMSFMAX))
C
C----------------------------------------------------------------------
C
      ITESTSFN = 0
C
      IF ( ITESTSFN.NE.0 ) WRITE (6,*) 'ICED ',ICED
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      DO I = 1,NFACEMAX
         NVERTICES(I) = NVERTICES8(I)
      END DO
C
      DO I = 1,NPANMAX
         NM(I) = NM8(I)
      END DO
C
      DO I = 1,NFACEMAX
         AFACE(I) = AFACE8(I)
         BFACE(I) = BFACE8(I)
         CFACE(I) = CFACE8(I)
         DFACE(I) = DFACE8(I)
         DO I1 = 1,NVERTMAX
C
            XEDGE(I1,I) = XEDGE8(I1,I)
            YEDGE(I1,I) = YEDGE8(I1,I)
            ZEDGE(I1,I) = ZEDGE8(I1,I)
         END DO
      END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  THIS CALL DOES SOME GEOMETRICAL TESTS (N.STEFANOU 98)
C
      CALL SFNPOLCHK(NFACE,NVERTICES,XEDGE,YEDGE,ZEDGE,NVERTMAX,
     &               NFACEMAX,IFLAG_POLYHEDRON)
C
      IF ( IFLAG_POLYHEDRON.EQ.1 ) RETURN
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      ISUM = 0
      DO L = 0,LMAX
         IS0 = (2*L+1)*(2*L+1)
         ISUM = ISUM + IS0
      END DO
C
      ISUMD = ISUM
      ALLOCATE (DMATL(ISUMD,NFACE))
C
      IF ( LMAX.GT.LSFMAX ) THEN
         WRITE (6,99004) ISUM,ISUMD,LMAX,LSFMAX
         CALL STOP_MESSAGE(ROUTINE,'LMAX > LSFMAX')
      END IF
C
      IPAN = 0
      IVTOT = 0
      CRT(1) = 9999.9999D0
C.......................................................................
C     C A L C U L A T I O N    O F    R O T A T I O N    M A T R I C E S
C.......................................................................
      DO IFACE = 1,NFACE
C
         A1 = AFACE(IFACE)
         A2 = BFACE(IFACE)
         A3 = CFACE(IFACE)
         A4 = DFACE(IFACE)
C
         NVERT = NVERTICES(IFACE)
         NEWSCH(IFACE) = 1
C        THIS is ALWAYS ONE!
C
         Z(1) = A1/A4
         Z(2) = A2/A4
         Z(3) = A3/A4
C
         DO IVERT = 1,NVERT
C
            V(1,IVERT) = XEDGE(IVERT,IFACE)
            V(2,IVERT) = YEDGE(IVERT,IFACE)
            V(3,IVERT) = ZEDGE(IVERT,IFACE)
C
         END DO
C
         CALL SFNCRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,CRT,FA,FB,FD,R0,RD,
     &                ISIGNU,NTT,ALPHA,BETA,GAMMA,NPANMAX,NVERTMAX,
     &                NFACEMAX,NVTOTD)
C
         IF ( IPRINT.GT.0 ) WRITE (6,99003) IFACE,NTT(IFACE)
         IF ( NEWSCH(IFACE).EQ.0 ) WRITE (6,99005)
C
         CALL SFNDROTREAL(LMAX,ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE),
     &                    LSFMAX,ISUMD,DMATL(1,IFACE))
C
      END DO
C
C.......................................................................
C     D E F I N I T I O N    O F    T H E    S U I T A B L E    M E S H
C.......................................................................
      NVTOT = IVTOT
      NPAN = IPAN
C
      CALL SFNMESH(CRT,NPAN,NM,XRN,DRN,MESHN,NRSF0,KEYPAN,NPANMAX,
     &             NRSFMAX)
C
      IF ( ITESTSFN.NE.0 ) THEN
         WRITE (6,99002)
         DO IV = 1,NVTOT
            WRITE (6,99001) IV,FA(IV)/PI,FB(IV)/PI,FD(IV)/PI,RD(IV),
     &                      ISIGNU(IV)
         END DO
      END IF
C.......................................................................
C     E X P A N S I O N    C O E F F I C I E N T S
C.......................................................................
C
      ICEDIM = ((LSFMAX+1)*(LSFMAX+2))/2
      LMA2D = LSFMAX/2 + 1
      CALL SFNCCOEF(LMAX,CL,C,LSFMAX,LMA2D,ITESTSFN,ICEDIM,ICD)
      IVTOT = 0
      DO IFACE = 1,NFACE
         NTET = NTT(IFACE)
         DO ITET = 1,NTET
            IVTOT = IVTOT + 1
            RUPSQ(IVTOT) = SQRT((RD(IVTOT)-R0(IFACE))*(RD(IVTOT)+R0(
     &                     IFACE)))
         END DO
      END DO
      DO IBM = 1,NLMSFMAX
         ISW(IBM) = 0
      END DO
C
C.......................................................................
C     L O O P    O V E R    R A D I A L    M E S H    P O I N T S
C.......................................................................
      DO N = 1,MESHN
         R = XRN(N)
         DO IBM = 1,NLMSFMAX
            B(IBM) = 0.0D0
         END DO
         IVTOT = 0
C.......................................................................
C     L O O P    O V E R    P Y R A M I D S
C.......................................................................
         DO IFACE = 1,NFACE
            NTET = NTT(IFACE)
C
            IF ( R.GT.R0(IFACE) ) THEN
               IF ( NEWSCH(IFACE).EQ.1 ) THEN
                  ARG1 = R0(IFACE)/R
                  RDOWN = SQRT((R-R0(IFACE))*(R+R0(IFACE)))
                  DO I = 0,LMAX
                     S(0,I) = 0.0D0
                  END DO
                  DO M = 1,LMAX
                     DO I = 0,LMAX - M
                        S(-M,I) = 0.0D0
                        S(M,I) = 0.0D0
                     END DO
                  END DO
C.......................................................................
C     L O O P     O V E R     T E T R A H E D R A
C.......................................................................
                  DO ITET = 1,NTET
                     IVTOT = IVTOT + 1
                     IF ( R.LE.RD(IVTOT) ) THEN
                        CALL SFNPINTG(FA(IVTOT),FB(IVTOT),DLT,S1,LMAX,
     &                                ISIGNU(IVTOT),ARG1,FD(IVTOT),0,
     &                                LSFMAX,NDIM)
                        DO I = 0,LMAX
                           S(0,I) = S(0,I) + S1(0,I)
                        END DO
                        DO M = 1,LMAX
                           DO I = 0,LMAX - M
                              S(-M,I) = S(-M,I) + S1(-M,I)
                              S(M,I) = S(M,I) + S1(M,I)
                           END DO
                        END DO
                     ELSE
                        RAP = RUPSQ(IVTOT)/RDOWN
                        ARG2 = RUPSQ(IVTOT)/R0(IFACE)
                        FK = FD(IVTOT) - ACOS(RAP)
                        FL = FD(IVTOT) + ACOS(RAP)
C
                        FK = DMAX1(FA(IVTOT),FK)
                        FL = DMAX1(FA(IVTOT),FL)
                        FK = DMIN1(FB(IVTOT),FK)
                        FL = DMIN1(FB(IVTOT),FL)
                        CALL SFNPINTG(FA(IVTOT),FK,DLT,S1,LMAX,
     &                                ISIGNU(IVTOT),ARG1,FD(IVTOT),0,
     &                                LSFMAX,NDIM)
                        CALL SFNPINTG(FK,FL,DLT,S2,LMAX,ISIGNU(IVTOT),
     &                                ARG2,FD(IVTOT),1,LSFMAX,NDIM)
                        CALL SFNPINTG(FL,FB(IVTOT),DLT,S3,LMAX,
     &                                ISIGNU(IVTOT),ARG1,FD(IVTOT),0,
     &                                LSFMAX,NDIM)
                        DO I = 0,LMAX
                           S(0,I) = S(0,I) + S1(0,I) + S2(0,I) + S3(0,I)
                        END DO
                        DO M = 1,LMAX
                           DO I = 0,LMAX - M
                              S(-M,I) = S(-M,I) + S1(-M,I) + S2(-M,I)
     &                                  + S3(-M,I)
                              S(M,I) = S(M,I) + S1(M,I) + S2(M,I)
     &                                 + S3(M,I)
                           END DO
                        END DO
                     END IF
                  END DO
               ELSE
                  IVTOT = IVTOT + NTET
               END IF
C.......................................................................
C    I N T E G R A L   E X P A N S I O N    B A C K - R O T A T I O N
C.......................................................................
               IB = 0
               IC = 0
               ICE = 0
               ISU = 0
               DO L = 0,LMAX
                  IB = IB + L + 1
                  DO MP = L,1, - 1
                     RSUM(MP,1) = 0.0D0
                     RSUM(MP,2) = 0.0D0
                     ICE = ICE + 1
                     K0 = (L+MP+1)/2
                     DO K = L,K0, - 1
                        IS = 2*K - L - MP
                        IC = IC + 1
C
                        RSUM(MP,2) = RSUM(MP,2) + CL(IC)*S(-MP,IS)
                        RSUM(MP,1) = RSUM(MP,1) + CL(IC)*S(MP,IS)
                     END DO
                     RSUM(MP,2) = RSUM(MP,2)*C(ICE)
                     RSUM(MP,1) = RSUM(MP,1)*C(ICE)
                  END DO
                  RSUM(0,1) = 0.0D0
                  ICE = ICE + 1
                  K0 = (L+1)/2
                  DO K = L,K0, - 1
                     IS = 2*K - L
                     IC = IC + 1
C
                     RSUM(0,1) = RSUM(0,1) + CL(IC)*S(0,IS)
                  END DO
                  RSUM(0,1) = RSUM(0,1)*C(ICE)
                  IMAX = 1
                  M = 0
 5                CONTINUE
                  DO I = 1,IMAX
                     MO = (3-2*I)*M
                     IBM = IB + MO
                     LOFM(IBM) = L
                     MOFM(IBM) = MO
                     IPMAX = 1
                     MP = 0
 6                   CONTINUE
                     DO IP = 1,IPMAX
                        ISU = ISU + 1
                        B(IBM) = B(IBM) + RSUM(MP,IP)*DMATL(ISU,IFACE)
                     END DO
                     IPMAX = 2
                     MP = MP + 1
                     IF ( MP.LE.L ) GOTO 6
                  END DO
                  IMAX = 2
                  M = M + 1
                  IF ( M.LE.L ) GOTO 5
                  IB = IB + L
               END DO
            ELSE
               IVTOT = IVTOT + NTET
               DO I = 0,LMAX
                  S(0,I) = 0.0D0
               END DO
               DO M = 1,LMAX
                  DO I = 0,LMAX - M
                     S(-M,I) = 0.0D0
                     S(M,I) = 0.0D0
                  END DO
               END DO
            END IF
         END DO
C.......................................................................
C     D E F I N E S   A N D    S A V E S   S H A P E    F U N C T I O N
C.......................................................................
         B(1) = SQRT_4PI - B(1)/SQRT_4PI
         DO IBM = 2,NLMSFMAX
            B(IBM) = -B(IBM)/SQRT_4PI
         END DO
         DO IBM = 1,NLMSFMAX
            IF ( ABS(B(IBM)).GT.1.D-6 ) ISW(IBM) = 1
C
            FLMSF_LOC(N,IBM) = B(IBM)
C
         END DO
      END DO
      NFUN = 0
      DO IBM = 1,NLMSFMAX
         IF ( ISW(IBM).EQ.1 ) NFUN = NFUN + 1
      END DO
C THIS IS FOR DIFFERENT SHELLS...
      NCELL = 1
      SCL = 1.D0
C
      NCELL_S = NCELL
      SCALE_S = SCL
      NPAN_S = NPAN - 1
      MESHN_S = MESHN
      DO IPAN = 1,NPAN - 1
         NM_S(IPAN) = NM(IPAN)
      END DO
      DO N = 1,MESHN
         XRN_S(N) = XRN(N)
         DRN_S(N) = DRN(N)
      END DO
      NFUN_S = NFUN
C
      ICOUNT = 0
      DO IBM = 1,NLMSFMAX
C
         IF ( ISW(IBM).NE.0 ) THEN
C
            ICOUNT = ICOUNT + 1
            LM0 = LOFM(IBM)*LOFM(IBM) + LOFM(IBM) + MOFM(IBM) + 1
            LMIFUN_S(ICOUNT) = LM0
            DO N = 1,MESHN
C
               FLMSF_S(N,ICOUNT) = FLMSF_LOC(N,IBM)
C
            END DO
C
         END IF
      END DO
C
      DEALLOCATE (B,C,S,V,R0,S1,S2,S3,NVERTICES,FA,FB,FD,CL,RD,NM)
      DEALLOCATE (DRN,CRT,ISW,RSUM,NTT,XRN,BETA,LOFM,MOFM,AFACE,BFACE)
      DEALLOCATE (CFACE,DFACE,GAMMA,ALPHA,XEDGE,YEDGE,ZEDGE,DMATL)
      DEALLOCATE (RUPSQ,NEWSCH,ISIGNU)
C
      RETURN
C
99001 FORMAT (I10,4F10.4,I10)
99002 FORMAT (//15X,'FA/PI',5X,'FB/PI',5X,'FD/PI',6X,'RD',8X,'ISIGNU'/)
99003 FORMAT (/,10X,I3,'-TH PYRAMID SUBDIVIDED IN ',I3,' TETRAHEDRA')
99004 FORMAT (23X,'FROM MAIN : ISUM=',I7,'  GREATER THAN DIMENSIONED',
     &        I7/23X,'       OR   LMAX=',I7,
     &        '  GREATER THAN DIMENSIONED',I7)
99005 FORMAT (11X,'BUT IS IDENTICAL TO A PREVIOUS ONE.')
      END
C*==sfncrit.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNCRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,CRT,FA,FB,FD,R0,RD,
     &                   ISIGNU,NTT,ALPHA,BETA,GAMMA,NPANMAX,NVERTMAX,
     &                   NFACEMAX,NVTOTD)
C   ********************************************************************
C   *                                                                  *
C   * THIS ROUTINE CALCULATES THE CRITICAL POINTS 'CRT' OF THE SHAPE   *
C   * FUNCTIONS DUE TO THE FACE: Z(1)*X + Z(2)*Y + Z(3)*Z = 1          *
C   * THE FACE IS ROTATED THROUGH THE APPROPRIATE EULER ANGLES TO BE   *
C   * PERPENDICULAR TO THE Z-AXIS. A FURTHER SUBDIVISION OF THE CEN-   *
C   * TRAL PYRAMID INTO ELEMENTARY TETRAHEDRA IS PERFORMED.            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNCRIT')
C
C Dummy arguments
C
      INTEGER IFACE,IPAN,IVTOT,NFACEMAX,NPANMAX,NVERT,NVERTMAX,NVTOTD
      REAL*8 ALPHA(NFACEMAX),BETA(NFACEMAX),CRT(NPANMAX),FA(NVTOTD),
     &       FB(NVTOTD),FD(NVTOTD),GAMMA(NFACEMAX),R0(NFACEMAX),
     &       RD(NVTOTD),V(3,NVERTMAX),Z(3)
      INTEGER ISIGNU(NVTOTD),NTT(NFACEMAX)
C
C Local variables
C
      REAL*8 A1,A2,A3,ARG,CF1,CF2,CF3,CO,CRRT,D1,D2,DD,DOWN,F1,F2,FF,
     &       OMEGA,ORIGIN(3),RDD,RDV(3),S,SF1,SF2,SF3,UP,VZ(:,:),XJ,YJ,
     &       ZMOD2,ZVMOD
      INTEGER I,IBACK,ICORN,IN(:),INEW,IP,IVERT,IVERT1,IVERTP,IX
      LOGICAL INSIDE
      EXTERNAL SFNEULER,SFNPERP
C
C*** End of declarations rewritten by SPAG
C
      DATA ORIGIN/3*0.0D0/
C-----------------------------------------------------------------------
      ALLOCATABLE IN,VZ
      ALLOCATE (IN(NVERTMAX),VZ(3,NVERTMAX))
C
      NTT(IFACE) = 0
      IF ( IPRINT.GT.0 ) WRITE (6,99004) IFACE,(Z(I),I=1,3)
      ZMOD2 = Z(1)*Z(1) + Z(2)*Z(2) + Z(3)*Z(3)
C
      IF ( ZMOD2.LE.1.0D-6 ) THEN
C
         WRITE (6,99001) IFACE,(Z(I),I=1,3)
         CALL STOP_MESSAGE(ROUTINE,'ZMOD2 <= 1.0D-6')
C
      ELSE
C
         S = 2.0D0*PI
         DO I = 1,3
            Z(I) = Z(I)/ZMOD2
         END DO
         IX = 1
         ZVMOD = SQRT((V(1,1)-Z(1))**2+(V(2,1)-Z(2))**2+(V(3,1)-Z(3))
     &           **2)
         IF ( ZVMOD.LT.1.0D-6 ) IX = 2
C
         CALL SFNEULER(Z,V(1,IX),IFACE,ALPHA,BETA,GAMMA,NFACEMAX)
C
         IF ( IPRINT.GT.0 ) WRITE (6,99005) ALPHA(IFACE)/PI,BETA(IFACE)
     &                             /PI,GAMMA(IFACE)/PI
C
         CALL SFNROTATE(V,VZ,IFACE,NVERT,ALPHA,BETA,GAMMA,NFACEMAX,
     &                  NVERTMAX)
C
         R0(IFACE) = 1.0D0/SQRT(ZMOD2)
         ICORN = 0
         IF ( IPRINT.GT.0 ) WRITE (6,99008)
C
         DO IVERT = 1,NVERT
            IF ( ABS(R0(IFACE)-VZ(3,IVERT)).GT.1.0D-6 ) GOTO 100
C.......................................................................
C     D I S T A N C E S   O F   V E R T I C E S   F R O M   C E N T E R
C.......................................................................
            CRRT = SQRT(VZ(1,IVERT)**2+VZ(2,IVERT)**2+VZ(3,IVERT)**2)
            INEW = 1
            DO IP = 1,IPAN
               IF ( ABS(CRRT-CRT(IP)).LT.1.0D-6 ) INEW = 0
            END DO
            IF ( INEW.EQ.1 ) THEN
               IPAN = IPAN + 1
               IF ( IPAN.GT.NPANMAX ) GOTO 200
               CRT(IPAN) = CRRT
            END IF
            IVERTP = IVERT + 1
            IF ( IVERT.EQ.NVERT ) IVERTP = 1
C.......................................................................
C     D I S T A N C E S   O F   E D G E S   F R O M   C E N T E R
C.......................................................................
            CALL SFNPERP(ORIGIN,VZ(1,IVERT),VZ(1,IVERTP),RDV,INSIDE)
            RDD = SQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2)+RDV(3)*RDV(3))
            IF ( INSIDE ) THEN
               INEW = 1
               DO IP = 1,IPAN
                  IF ( ABS(RDD-CRT(IP)).LT.1.0D-4 ) INEW = 0
               END DO
               IF ( INEW.EQ.1 ) THEN
                  IPAN = IPAN + 1
                  IF ( IPAN.GT.NPANMAX ) GOTO 200
                  CRT(IPAN) = RDD
               END IF
            END IF
            A1 = SQRT(VZ(1,IVERT)*VZ(1,IVERT)+VZ(2,IVERT)*VZ(2,IVERT))
            A2 = SQRT(VZ(1,IVERTP)*VZ(1,IVERTP)+VZ(2,IVERTP)
     &           *VZ(2,IVERTP))
            DOWN = A1*A2
            UP = VZ(1,IVERT)*VZ(1,IVERTP) + VZ(2,IVERT)*VZ(2,IVERTP)
            IF ( DOWN.GT.1.0D-06 ) THEN
               ARG = UP/DOWN
               IF ( ABS(ARG).GE.1.0D0 ) ARG = SIGN(1D0,ARG)
               OMEGA = ACOS(ARG)
               S = S - OMEGA
               IF ( ABS(OMEGA-PI).GT.1.0D-06 ) THEN
C.......................................................................
C     S U B D I V I S I O N    I N T O    T E T R A H E D R A
C.......................................................................
                  NTT(IFACE) = NTT(IFACE) + 1
                  IVTOT = IVTOT + 1
                  IF ( IPRINT.GT.0 ) THEN
                     WRITE (6,99006) IVTOT,IVERT,(VZ(I,IVERT),I=1,3)
                     WRITE (6,99007) IVERTP,(VZ(I,IVERTP),I=1,3)
                  END IF
                  A3 = SQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2))
                  RD(IVTOT) = RDD
                  ISIGNU(IVTOT) = 1
                  CF1 = VZ(1,IVERT)/A1
                  CF2 = VZ(1,IVERTP)/A2
                  SF1 = VZ(2,IVERT)/A1
                  SF2 = VZ(2,IVERTP)/A2
                  CF3 = RDV(1)/A3
                  SF3 = RDV(2)/A3
                  IF ( ABS(SF1).LT.1.0D-5 .AND. ABS(CF1+1.0D0)
     &                 .LT.1.0D-05 ) THEN
                     F1 = PI
                  ELSE
                     F1 = 2.0D0*ATAN2(SF1,CF1+1.0D0)
                  END IF
                  IF ( ABS(SF2).LT.1.0D-5 .AND. ABS(CF2+1.0D0)
     &                 .LT.1.0D-05 ) THEN
                     F2 = PI
                  ELSE
                     F2 = 2.0D0*ATAN2(SF2,CF2+1.0D0)
                  END IF
                  IF ( ABS(SF3).LT.1.0D-5 .AND. ABS(CF3+1.0D0)
     &                 .LT.1.0D-05 ) THEN
                     FD(IVTOT) = PI
                  ELSE
                     FD(IVTOT) = 2.0D0*ATAN2(SF3,CF3+1.0D0)
                  END IF
C
                  FA(IVTOT) = DMIN1(F1,F2)
                  FB(IVTOT) = DMAX1(F1,F2)
                  IF ( (FB(IVTOT)-FA(IVTOT)).GT.PI ) THEN
                     FF = FA(IVTOT) + 2.0D0*PI
                     FA(IVTOT) = FB(IVTOT)
                     FB(IVTOT) = FF
                  END IF
                  IF ( (FA(IVTOT)-FD(IVTOT)).GT.PI ) FD(IVTOT)
     &                 = 2.0D0*PI + FD(IVTOT)
                  IF ( (FD(IVTOT)-FA(IVTOT)).GT.PI ) FD(IVTOT)
     &                 = -2.0D0*PI + FD(IVTOT)
               END IF
            ELSE
               ICORN = 1
            END IF
         END DO
C.......................................................................
C     F O O T   O F   T H E    P E R P E N D I C U L A R   TO    T H E
C     F A C E   O U T S I D E   O R  I N S I D E   T H E   P O L Y G O N
C.......................................................................
         IF ( S.LT.1.0D-06 .OR. ICORN.EQ.1 ) THEN
            INEW = 1
            DO IP = 1,IPAN
               IF ( ABS(R0(IFACE)-CRT(IP)).LT.1.0D-4 ) INEW = 0
            END DO
            IF ( INEW.EQ.1 ) THEN
               IPAN = IPAN + 1
               IF ( IPAN.GT.NPANMAX ) GOTO 200
               CRT(IPAN) = R0(IFACE)
            END IF
         ELSE
            DO IVERT1 = 1,NVERT
               IN(IVERT1) = 0
               DO IVERT = 1,NVERT
                  IVERTP = IVERT + 1
                  IF ( IVERT.EQ.NVERT ) IVERTP = 1
                  IF ( IVERT.NE.IVERT1 .AND. IVERTP.NE.IVERT1 ) THEN
                     DOWN = VZ(2,IVERT1)*(VZ(1,IVERTP)-VZ(1,IVERT))
     &                      - VZ(1,IVERT1)*(VZ(2,IVERTP)-VZ(2,IVERT))
                     IF ( ABS(DOWN).GT.1.0D-06 ) THEN
                        UP = VZ(1,IVERT1)
     &                       *(VZ(2,IVERT)*(VZ(1,IVERTP)+VZ(1,IVERT))
     &                       -VZ(1,IVERT)*(VZ(2,IVERTP)+VZ(2,IVERT)))
                        XJ = UP/DOWN
                        YJ = XJ*VZ(2,IVERT1)/VZ(1,IVERT1)
                        DD = (VZ(1,IVERTP)-VZ(1,IVERT))
     &                       **2 + (VZ(2,IVERTP)-VZ(2,IVERT))**2
                        D1 = (XJ-VZ(1,IVERT))**2 + (YJ-VZ(2,IVERT))**2
                        D2 = (XJ-VZ(1,IVERTP))**2 + (YJ-VZ(2,IVERTP))**2
C
                        CO = DD - DMAX1(D1,D2)
                        IF ( CO.GT.1.0D-06 ) THEN
                           IN(IVERT1) = 1
                           EXIT
                        END IF
                     END IF
                  END IF
               END DO
            END DO
            IBACK = IVTOT - NVERT
            DO IVERT = 1,NVERT
               IBACK = IBACK + 1
               IVERTP = IVERT + 1
               IF ( IVERT.EQ.NVERT ) IVERTP = 1
               IF ( IN(IVERT).EQ.0 .AND. IN(IVERTP).EQ.0 ) ISIGNU(IBACK)
     &              = -1
            END DO
         END IF
C
         DEALLOCATE (IN,VZ)
         RETURN
      END IF
C
 100  CONTINUE
      WRITE (6,99002) IFACE,R0(IFACE),
     &                ((VZ(I,IVERT),I=1,3),IVERT=1,NVERT)
      CALL STOP_MESSAGE(ROUTINE,' XXX ')
C
 200  CONTINUE
      WRITE (6,99003) IPAN,NPANMAX
      CALL STOP_MESSAGE(ROUTINE,' XXX ')
C
99001 FORMAT (//,13X,'FATAL ERROR FROM CRIT: THE',I3,
     &        '-TH FACE OF THE POLYHEDRON PASSES THROUGH THE CENTER',/,
     &        13X,'(',3E14.7,' )')
99002 FORMAT (//,13X,'FATAL ERROR FROM CRIT: THE VERTICES OF THE',I3,
     &        '-TH ROTATED POLYGON DO NOT LIE ON THE PLANE:',E13.6,
     &        ' *Z = 1'/30(/13X,3E13.6))
99003 FORMAT (//,13X,'ERROR FROM CRIT: NUMBER OF PANNELS=',I5,
     &        ' GREATER THAN DIMENSIONED=',I5)
99004 FORMAT (//,80('*'),/,10X,'FACE:',I3,' EQUATION:',F10.4,'*X +',
     &        F10.4,'*Y +',F10.4,'*Z  =  1')
99005 FORMAT (10X,'ROTATION ANGLES  :',3(F10.4,4X)/)
99006 FORMAT (10X,I5,'       VZ(',I2,')  =  (',3F10.4,' )')
99007 FORMAT (10X,5X,'       VZ(',I2,')  =  (',3F10.4,' )')
99008 FORMAT (/,10X,'TETRAHEDRON',22X,'COORDINATES',/,10X,11('*'),22X,
     &        11('*')/)
      END
C*==sfnrotate.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNROTATE(V,VZ,IFACE,NVERT,ALPHA,BETA,GAMMA,NFACEMAX,
     &                     NVERTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  THIS ROUTINE PERFORMS THE ROTATION OF NVERT VECTORS THROUGH     *
C   *  THE EULER ANGLES: ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE).        *
C   *   V (I,IVERT) : INPUT   VECTORS                                  *
C   *   VZ(I,IVERT) : ROTATED VECTORS                                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IFACE,NFACEMAX,NVERT,NVERTMAX
      REAL*8 ALPHA(NFACEMAX),BETA(NFACEMAX),GAMMA(NFACEMAX),
     &       V(3,NVERTMAX),VZ(3,NVERTMAX)
C
C Local variables
C
      REAL*8 A(3,3),CA,CB,CG,SA,SB,SG
      INTEGER I,IVERT,J
C
C*** End of declarations rewritten by SPAG
C
      CA = COS(ALPHA(IFACE))
      SA = SIN(ALPHA(IFACE))
      CB = COS(BETA(IFACE))
      SB = SIN(BETA(IFACE))
      CG = COS(GAMMA(IFACE))
      SG = SIN(GAMMA(IFACE))
      A(1,1) = CA*CB*CG - SA*SG
      A(2,1) = SA*CB*CG + CA*SG
      A(3,1) = -SB*CG
      A(1,2) = -CA*CB*SG - SA*CG
      A(2,2) = -SA*CB*SG + CA*CG
      A(3,2) = SB*SG
      A(1,3) = CA*SB
      A(2,3) = SA*SB
      A(3,3) = CB
      DO IVERT = 1,NVERT
         DO I = 1,3
            VZ(I,IVERT) = 0.0D0
            DO J = 1,3
               VZ(I,IVERT) = VZ(I,IVERT) + A(I,J)*V(J,IVERT)
            END DO
         END DO
      END DO
      END
C*==sfneuler.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNEULER(Z,XX,IFACE,ALPHA,BETA,GAMMA,NFACEMAX)
C-----------------------------------------------------------------------
C     GIVEN TWO DISTINCT POINTS (Z(1),Z(2),Z(3)) AND (XX(1),XX(2),XX(3))
C     THIS ROUTINE DEFINES  A LOCAL COORDINATE  SYSTEM WITH THE  Z- AXIS
C     PASSING THROUGH (Z(1),Z(2),Z(3))  AND THE X- AXIS PARALLEL TO  THE
C     VECTOR : (XX(1)-Z(1),XX(2)-Z(2),XX(3)-Z(3)).
C     THE EULER ANGLES ROTATING THIS LOCAL COORDINATE SYSTEM BACK TO THE
C     ORIGINAL FRAME OF REFERENCE ARE CALCULATED
C-----------------------------------------------------------------------
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNEULER')
C
C Dummy arguments
C
      INTEGER IFACE,NFACEMAX
      REAL*8 ALPHA(NFACEMAX),BETA(NFACEMAX),GAMMA(NFACEMAX),XX(3),Z(3)
C
C Local variables
C
      REAL*8 CA,CG,P,RX,RZ,RZP,S,SA,SG,X(3),Y(3)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      IF ( IFACE.GT.NFACEMAX ) THEN
         WRITE (6,99001) IFACE,NFACEMAX
         CALL STOP_MESSAGE(ROUTINE,'IFACE > NFACE')
      ELSE
         DO I = 1,3
            X(I) = XX(I) - Z(I)
         END DO
         RX = SQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
         RZ = SQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
         IF ( RX.GE.1D-6 .AND. RZ.GE.1D-6 ) THEN
            S = X(1)*Z(1) + X(2)*Z(2) + X(3)*Z(3)
            IF ( S.LE.1D-6 ) THEN
               P = SQRT(Z(1)*Z(1)+Z(2)*Z(2))
               DO I = 1,3
                  X(I) = X(I)/RX
                  Z(I) = Z(I)/RZ
               END DO
               ALPHA(IFACE) = 0.0D0
               BETA(IFACE) = ACOS(Z(3))
               IF ( P.LT.1D-10 ) THEN
                  SG = -Z(3)*X(2)
                  CG = Z(3)*X(1)
                  IF ( ABS(SG).LT.1.0D-5 .AND. ABS(CG+1.0D0)
     &                 .LT.1.0D-05 ) THEN
                     GAMMA(IFACE) = PI
                  ELSE
                     GAMMA(IFACE) = 2.0D0*ATAN2(SG,CG+1.0D0)
                  END IF
                  DO I = 1,3
                     Z(I) = Z(I)*RZ
                  END DO
                  RETURN
               ELSE
                  RZP = RZ/P
                  Y(1) = Z(2)*X(3) - Z(3)*X(2)
                  Y(2) = Z(3)*X(1) - Z(1)*X(3)
                  Y(3) = Z(1)*X(2) - Z(2)*X(1)
                  SA = Y(3)*RZP
                  CA = X(3)*RZP
                  IF ( ABS(SA).LT.1.0D-5 .AND. ABS(CA+1.D0).LT.1.D-05 )
     &                 THEN
                     ALPHA(IFACE) = PI
                  ELSE
                     ALPHA(IFACE) = 2.0D0*ATAN2(SA,CA+1.0D0)
                  END IF
                  SG = Z(2)*RZP
                  CG = -Z(1)*RZP
                  IF ( ABS(SG).LT.1.0D-5 .AND. ABS(CG+1.0D0)
     &                 .LT.1.0D-05 ) THEN
                     GAMMA(IFACE) = PI
                  ELSE
                     GAMMA(IFACE) = 2.0D0*ATAN2(SG,CG+1.0D0)
                  END IF
                  DO I = 1,3
                     Z(I) = Z(I)*RZ
                  END DO
                  RETURN
               END IF
            END IF
         END IF
         WRITE (6,99002) (X(I),I=1,3),(Z(I),I=1,3)
         CALL STOP_MESSAGE(ROUTINE,' XXX ')
      END IF
C
99001 FORMAT (//13X,'NUMBER OF FACES:',I5,' GREATER THAN DIMENSIONED',
     &        I5)
99002 FORMAT (/13X,'FROM EULER,ILLEGAL VECTORS:'/13X,2(' (',3E13.6,' )')
     &        )
      END
C*==sfnperp.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNPERP(R0,R1,R2,RD,INSIDE)
C-----------------------------------------------------------------------
C     GIVEN  TWO  DISTINCT  POINTS   R1 , R2, THIS  ROUTINE CALCULATES
C     THE COORDINATES  OF THE FOOT OF  THE  PERPENDICULAR FROM A POINT
C     R0 TO THE LINE JOINING R1   AND  R2. THE LOGICAL VARIABLE INSIDE
C     GIVES THE ADDITIONAL INFORMATION WHETHER THE FOOT OF THE PERPEN-
C     DICULAR LIES WITHIN THE SEGMENT OR NOT.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNPERP')
C
C Dummy arguments
C
      LOGICAL INSIDE
      REAL*8 R0(3),R1(3),R2(3),RD(3)
C
C Local variables
C
      REAL*8 CO,D,D1,D2,DA,DB,DC,DX,DY,DZ,S
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DX = R2(1) - R1(1)
      DY = R2(2) - R1(2)
      DZ = R2(3) - R1(3)
      S = R0(1)*DX + R0(2)*DY + R0(3)*DZ
      D = DX*DX + DY*DY + DZ*DZ
      IF ( SQRT(D).LT.1.0D-6 ) THEN
         WRITE (6,99001) (R1(I),I=1,3),(R2(I),I=1,3)
         CALL STOP_MESSAGE(ROUTINE,'IDENTICAL POINTS')
      END IF
      DA = S*DX + DY*(R1(1)*R2(2)-R1(2)*R2(1))
     &     + DZ*(R1(1)*R2(3)-R1(3)*R2(1))
      DB = S*DY + DZ*(R1(2)*R2(3)-R1(3)*R2(2))
     &     + DX*(R1(2)*R2(1)-R1(1)*R2(2))
      DC = S*DZ + DX*(R1(3)*R2(1)-R1(1)*R2(3))
     &     + DY*(R1(3)*R2(2)-R1(2)*R2(3))
      RD(1) = DA/D
      RD(2) = DB/D
      RD(3) = DC/D
      D1 = (RD(1)-R1(1))**2 + (RD(2)-R1(2))**2 + (RD(3)-R1(3))**2
      D2 = (RD(1)-R2(1))**2 + (RD(2)-R2(2))**2 + (RD(3)-R2(3))**2
C
      CO = D - DMAX1(D1,D2)
      INSIDE = .FALSE.
      IF ( CO.GT.1.0D-6 ) INSIDE = .TRUE.
      RETURN
99001 FORMAT (///,10X,2('(',3E14.6,')',3X))
      END
C*==sfnmesh.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNMESH(CRT,NPAN,NM,XRN,DRN,MESHN,NRSF0,KEYPAN,NPANMAX,
     &                   NRSFMAX)
C-----------------------------------------------------------------------
C     THIS ROUTINE DEFINES A UNIQUE SUITABLE RADIAL MESH 'XRN,DRN' OF
C     'MESHN' POINTS,DISTRIBUTED INTO 'NPAN' PANNELS  DEFINED  BY THE
C     CRITICAL POINTS 'CRT'
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNMESH')
C
C Dummy arguments
C
      INTEGER KEYPAN,MESHN,NPAN,NPANMAX,NRSF0,NRSFMAX
      REAL*8 CRT(NPANMAX),DRN(NRSFMAX),XRN(NRSFMAX)
      INTEGER NM(NPANMAX)
C
C Local variables
C
      REAL*8 C,D
      INTEGER IORD,IPAN,K,N1,N2
C
C*** End of declarations rewritten by SPAG
C
      DO IORD = 1,NPAN
         C = CRT(IORD)
         DO IPAN = NPAN,IORD, - 1
            IF ( CRT(IPAN).LE.C ) THEN
               C = CRT(IPAN)
               CRT(IPAN) = CRT(IORD)
               CRT(IORD) = C
            END IF
         END DO
      END DO
C
C     CALCULATE AN APPROPRIATE MESH
C
      CALL SFNMESH0(KEYPAN,CRT,NM,NPAN,NRSF0,NPANMAX)
C
      WRITE (6,99002)
      N2 = 0
      DO IPAN = 1,NPAN - 1
         WRITE (6,99003) IPAN,CRT(IPAN),CRT(IPAN+1),NM(IPAN)
         N1 = N2 + 1
         N2 = N2 + NM(IPAN)
         IF ( NRSFMAX.LT.N2 ) GOTO 100
C
         C = (CRT(IPAN+1)-CRT(IPAN))/DBLE(N2-N1)
         D = CRT(IPAN) - C*DBLE(N1)
         DO K = N1,N2
            XRN(K) = C*DBLE(K) + D
            DRN(K) = C
         END DO
      END DO
      MESHN = N2
      WRITE (6,99004) MESHN
C      WRITE(6,101)(K,DRN(K),XRN(K),K=1,MESHN)
      RETURN
 100  CONTINUE
      WRITE (6,99001) NRSFMAX
      CALL STOP_MESSAGE(ROUTINE,'NXR IS TOO SMALL')
99001 FORMAT (/,10X,'NXR =',I6)
99002 FORMAT (/,24X,'suitable radial mesh',//,15X,'IPAN',6X,'from',10X,
     &        'to',10X,'POINTS',/)
99003 FORMAT (13X,I5,2F14.7,I10)
99004 FORMAT (13X,'total',28X,I10,/)
      END
C*==sfnpintg.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNPINTG(X1,X2,DLT,S,LMAX,ISI,ARG,FD,ITYPE,LSFMAX,NDIM)
C-----------------------------------------------------------------------
C     THIS ROUTINE  ACCOMPLISHES THE  FI-INTEGRATION  OF REAL  SPHERICAL
C     HARMONICS BY THE REPEATED SIMPSON'S METHOD , OR ANALYTICALLY ACCOR
C     DING TO THE VALUE OF ITYPE. THE OBTAINED RESULTS HAVE TO BE MULTI-
C     PLIED BY THE APPROPRIATE EXPANSION COEFFICIENTS.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNPINTG')
C
C Dummy arguments
C
      REAL*8 ARG,DLT,FD,X1,X2
      INTEGER ISI,ITYPE,LMAX,LSFMAX,NDIM
      REAL*8 S(-LSFMAX:LSFMAX,0:LSFMAX)
C
C Local variables
C
      INTEGER I,K,M,N
      REAL*8 THETA,W,WW(:),X,XX(:)
      EXTERNAL SFNRECUR,SFNRECUR0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE WW,XX
      ALLOCATE (WW(NDIM),XX(NDIM))
C
      IF ( LMAX.LE.LSFMAX ) THEN
         DO I = 0,LMAX
            S(0,I) = 0.0D0
         END DO
         DO M = 1,LMAX
            DO I = 0,LMAX - M
               S(-M,I) = 0.0D0
               S(M,I) = 0.0D0
            END DO
         END DO
         IF ( ITYPE.NE.0 ) THEN
C
            N = INT((X2-X1)/DLT) + 3
            IF ( N.GT.NDIM ) CALL STOP_MESSAGE(ROUTINE,'INCREASE NDIM')
            CALL GAULEG(X1,X2,XX,WW,N)
            DO K = 1,N
               X = XX(K)
               W = DBLE(ISI)*WW(K)
               THETA = ATAN(ARG/COS(X-FD))
               CALL SFNRECUR(LMAX,X,THETA,W,S,LSFMAX)
            END DO
            GOTO 99999
         END IF
      ELSE
         WRITE (6,99001) LMAX,LSFMAX
         CALL STOP_MESSAGE(ROUTINE,'LMAX > LSFMAX')
      END IF
      THETA = ACOS(ARG)
      CALL SFNRECUR0(LMAX,X1,THETA,-DBLE(ISI),S,LSFMAX)
      CALL SFNRECUR0(LMAX,X2,THETA,DBLE(ISI),S,LSFMAX)
C
      DEALLOCATE (WW,XX)
C
      RETURN
99001 FORMAT (3X,'LMAX=',I4,' greater than dimensioned',I4)
99999 CONTINUE
      END
C*==sfnrecur.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNRECUR(LMAX,X,THETA,FAC,S,LSFMAX)
C-----------------------------------------------------------------------
C     THIS ROUTINE IS USED TO PERFORM THE FI-INTEGRATION OF REAL SPHE-
C     RICAL HARMONICS .THE THETA-INTEGRATION IS PERFORMED ANALYTICALLY
C     USING RECURRENCE RELATIONS.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FAC,THETA,X
      INTEGER LMAX,LSFMAX
      REAL*8 S(-LSFMAX:LSFMAX,0:LSFMAX)
C
C Local variables
C
      REAL*8 C01(:),C02(:),C1,C2,CC,CCA(:),EL,EL0,OL,OL0,SS,SSA(:)
      INTEGER I,M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C01,C02,CCA,SSA
      ALLOCATE (C01(LSFMAX),C02(LSFMAX),CCA(LSFMAX+2),SSA(LSFMAX+2))
C
      SS = SIN(THETA)
      CC = COS(THETA)
      DO I = 1,LMAX
         C01(I) = FAC*SIN(DBLE(I)*X)
         C02(I) = FAC*COS(DBLE(I)*X)
      END DO
      DO I = 1,LMAX + 2
         SSA(I) = SS**I
         CCA(I) = CC**I
      END DO
      OL0 = (THETA-SS*CC)/2.0D0
      EL0 = 0.0D0
      DO M = 1,LMAX,2
         OL = OL0
         EL = EL0
         C1 = C01(M)
         C2 = C02(M)
         DO I = 0,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
            OL = (DBLE(I+1)*OL+SSA(M+2)*CCA(I+1))/DBLE(I+M+3)
         END DO
         EL0 = DBLE(M+2)*EL0/DBLE(M+3)
         OL0 = (DBLE(M+2)*OL0-SSA(M+2)*CC)/DBLE(M+3)
      END DO
      OL0 = SSA(3)/3.0D0
      EL0 = 0.0D0
      DO M = 1,LMAX,2
         OL = OL0
         EL = EL0
         C1 = C01(M)
         C2 = C02(M)
         DO I = 1,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+SSA(M+2)*CCA(I+1))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-SSA(M+2)*CCA(2))/DBLE(M+4)
         EL0 = DBLE(M+2)*EL0/DBLE(M+4)
      END DO
      OL0 = -CC
      EL0 = -1.0D0
      OL = OL0
      EL = EL0
      C2 = FAC
      DO I = 0,LMAX,2
         S(0,I) = S(0,I) + C2*(OL-EL)
         OL = (DBLE(I+1)*OL+SSA(2)*CCA(I+1))/DBLE(I+3)
         EL = DBLE(I+1)*EL/DBLE(I+3)
      END DO
      OL0 = (2.0D0*OL0-SSA(2)*CC)/3.0D0
      EL0 = 2.0D0*EL0/3.0D0
      DO M = 2,LMAX,2
         OL = OL0
         EL = EL0
         C1 = C01(M)
         C2 = C02(M)
         DO I = 0,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+SSA(M+2)*CCA(I+1))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-SSA(M+2)*CC)/DBLE(M+3)
         EL0 = DBLE(M+2)*EL0/DBLE(M+3)
      END DO
      OL0 = -CCA(2)/2.0D0
      EL0 = -0.5D0
      OL = OL0
      EL = EL0
      C2 = FAC
      DO I = 1,LMAX,2
         S(0,I) = S(0,I) + C2*(OL-EL)
         OL = (DBLE(I+1)*OL+SS*SS*CCA(I+1))/DBLE(I+3)
         EL = DBLE(I+1)*EL/DBLE(I+3)
      END DO
      OL0 = (2.0D0*OL0-SSA(2)*CCA(2))/4.0D0
      EL0 = 2.0D0*EL0/4.0D0
      DO M = 2,LMAX,2
         OL = OL0
         EL = EL0
         C1 = C01(M)
         C2 = C02(M)
         DO I = 1,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+SSA(M+2)*CCA(I+1))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-SSA(M+2)*CCA(2))/DBLE(M+4)
         EL0 = DBLE(M+2)*EL0/DBLE(M+4)
      END DO
C
      DEALLOCATE (C01,C02,CCA,SSA)
C
      END
C*==sfnrecur0.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNRECUR0(LMAX,X,THETA,FAC,S,LSFMAX)
C-----------------------------------------------------------------------
C     THIS ROUTINE IS USED TO PERFORM  THE  FI - INTEGRATION  OF REAL SP
C     RICAL HARMONICS ANALYTICALLY.  THE  THETA-INTEGRATION  IS   PERFOR
C     ALSO ANALYTICALLY USING RECURRENCE RELATIONS.(THETA IS FI-INDEPEND
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FAC,THETA,X
      INTEGER LMAX,LSFMAX
      REAL*8 S(-LSFMAX:LSFMAX,0:LSFMAX)
C
C Local variables
C
      REAL*8 C1,C2,CC,EL,EL0,OL,OL0,SS
      INTEGER I,M
C
C*** End of declarations rewritten by SPAG
C
      SS = SIN(THETA)
      CC = COS(THETA)
      OL0 = (THETA-SS*CC)/2.0D0
      EL0 = 0.0D0
      DO M = 1,LMAX,2
         OL = OL0
         EL = EL0
         C1 = -FAC*COS(DBLE(M)*X)/DBLE(M)
         C2 = FAC*SIN(DBLE(M)*X)/DBLE(M)
         DO I = 0,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
            OL = (DBLE(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DBLE(I+M+3)
         END DO
         EL0 = DBLE(M+2)*EL0/DBLE(M+3)
         OL0 = (DBLE(M+2)*OL0-(SS**(M+2))*CC)/DBLE(M+3)
      END DO
      OL0 = SS*SS*SS/3.0D0
      EL0 = 0.0D0
      DO M = 1,LMAX,2
         OL = OL0
         EL = EL0
         C1 = -FAC*COS(DBLE(M)*X)/DBLE(M)
         C2 = FAC*SIN(DBLE(M)*X)/DBLE(M)
         DO I = 1,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-(SS**(M+2))*CC*CC)/DBLE(M+4)
         EL0 = DBLE(M+2)*EL0/DBLE(M+4)
      END DO
      OL0 = -CC
      EL0 = -1.0D0
      OL = OL0
      EL = EL0
      C2 = FAC*X
      DO I = 0,LMAX,2
         S(0,I) = S(0,I) + C2*(OL-EL)
         OL = (DBLE(I+1)*OL+SS*SS*(CC**(I+1)))/DBLE(I+3)
         EL = DBLE(I+1)*EL/DBLE(I+3)
      END DO
      OL0 = (2.0D0*OL0-SS*SS*CC)/3.0D0
      EL0 = 2.0D0*EL0/3.0D0
      DO M = 2,LMAX,2
         OL = OL0
         EL = EL0
         C1 = -FAC*COS(DBLE(M)*X)/DBLE(M)
         C2 = FAC*SIN(DBLE(M)*X)/DBLE(M)
         DO I = 0,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-(SS**(M+2))*CC)/DBLE(M+3)
         EL0 = DBLE(M+2)*EL0/DBLE(M+3)
      END DO
      OL0 = -CC*CC/2.0D0
      EL0 = -0.5D0
      OL = OL0
      EL = EL0
      C2 = FAC*X
      DO I = 1,LMAX,2
         S(0,I) = S(0,I) + C2*(OL-EL)
         OL = (DBLE(I+1)*OL+SS*SS*(CC**(I+1)))/DBLE(I+3)
         EL = DBLE(I+1)*EL/DBLE(I+3)
      END DO
      OL0 = (2.0D0*OL0-SS*SS*CC*CC)/4.0D0
      EL0 = 2.0D0*EL0/4.0D0
      DO M = 2,LMAX,2
         OL = OL0
         EL = EL0
         C1 = -FAC*COS(DBLE(M)*X)/DBLE(M)
         C2 = FAC*SIN(DBLE(M)*X)/DBLE(M)
         DO I = 1,LMAX - M,2
            S(-M,I) = S(-M,I) + C1*(OL-EL)
            S(M,I) = S(M,I) + C2*(OL-EL)
            OL = (DBLE(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DBLE(I+M+3)
            EL = DBLE(I+1)*EL/DBLE(I+M+3)
         END DO
         OL0 = (DBLE(M+2)*OL0-(SS**(M+2))*CC*CC)/DBLE(M+4)
         EL0 = DBLE(M+2)*EL0/DBLE(M+4)
      END DO
      END
C*==sfnreduce.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNREDUCE(NMBR,IFMX,IFI,IEXP)
C-----------------------------------------------------------------------
C     THIS ROUTINE REDUCES A POSITIVE INTEGER   INPUT NUMBER 'NMBR'
C     TO A PRODUCT  OF  FIRST  NUMBERS 'IFI' , AT POWERS  'IEXP'.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNREDUCE')
C
C Dummy arguments
C
      INTEGER IFMX,NMBR
      INTEGER IEXP(IFMX),IFI(IFMX)
C
C Local variables
C
      INTEGER I,NMB
C
C*** End of declarations rewritten by SPAG
C
      IF ( NMBR.LE.0 ) THEN
         WRITE (6,99002) NMBR
         CALL STOP_MESSAGE(ROUTINE,'NON POSITIVE NUMBER')
      END IF
      DO I = 1,IFMX
         IEXP(I) = 0
      END DO
      IF ( NMBR.EQ.1 ) RETURN
      NMB = NMBR
      DO I = 1,IFMX
         IEXP(I) = 0
 50      CONTINUE
         IF ( MOD(NMB,IFI(I)).EQ.0 ) THEN
            NMB = NMB/IFI(I)
            IEXP(I) = IEXP(I) + 1
            GOTO 50
         ELSE IF ( NMB.EQ.1 ) THEN
            RETURN
         END IF
      END DO
      WRITE (6,99001) NMBR
      CALL STOP_MESSAGE(ROUTINE,'NMBR inconsistent')
99001 FORMAT (3X,I15,'  CANNOT BE REDUCED IN THE BASIS OF FIRST NUMBERS'
     &        ,' GIVEN'/20X,'INCREASE THE BASIS OF FIRST NUMBERS')
99002 FORMAT (3X,I15,'  NON POSITIVE NUMBER')
      END
C*==sfnccoef.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNCCOEF(LMAX,CL,COE,LSFMAX,LMA2D,ITESTSFN,ICED,ICD)
C-----------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE COEFFICIENTS OF A POLYNOMIAL EXPANSION
C     OF RENORMALIZED LEGENDRE FUNCTIONS IN POWERS OF COSINES.
C     THE POSSIBILITY OF OVERFLOW (HIGH LMAX) IS AVOIDED BY USING FACTO-
C     RIZED FORMS FOR THE NUMBERS.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNCCOEF')
      INTEGER IFMX
      PARAMETER (IFMX=25)
C
C Dummy arguments
C
      INTEGER ICD,ICED,ITESTSFN,LMA2D,LMAX,LSFMAX
      REAL*8 CL(ICD),COE(ICED)
C
C Local variables
C
      REAL*8 DOWN,UP,UPSQ
      INTEGER I1,IC,IC1,IC2,ICE,ICMAX,IE(IFMX,LMA2D),IEA(IFMX),IEB(IFMX)
     &        ,IED(IFMX),IEINT,IEMOD,IEUPSQ,IFI(IFMX),IL2P(IFMX),IR,IRE,
     &        ISI,JM0(IFMX),K,K0,L,L1(IFMX),L1ST(IFMX),L2(IFMX),L2P,
     &        L2ST(IFMX),LA,LB,LI,M
      EXTERNAL SFNREDUCE
C
C*** End of declarations rewritten by SPAG
C
      DATA IFI/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
     &     73,79,83,89,97/
C-----------------------------------------------------------------------
      ICMAX = 0
      DO L = 0,LMAX
         LI = L/2 + 1
         ICMAX = ICMAX + (LMAX+1-L)*LI
      END DO
      IF ( LMAX.GT.LSFMAX .OR. ICMAX.GT.ICD ) THEN
         WRITE (6,99004) LMAX,LSFMAX,ICMAX,ICD
         CALL STOP_MESSAGE(ROUTINE,'LMAX > LSFMAX .OR. ICMAX > ICD')
      END IF
      IF ( ITESTSFN.NE.0 ) WRITE (6,99003) ICMAX
      ICE = 0
      IC = 1
      L = 0
      DO I1 = 1,IFMX
         L1ST(I1) = 0
         L2ST(I1) = 0
      END DO
C
 100  CONTINUE
      L2P = 2*L + 1
      CALL SFNREDUCE(L2P,IFMX,IFI,IL2P)
      IL2P(1) = IL2P(1) + 1
      M = L
      DO I1 = 1,IFMX
         L1(I1) = L1ST(I1)
         L2(I1) = L2ST(I1)
         JM0(I1) = 0
      END DO
 200  CONTINUE
      ICE = ICE + 1
      ISI = 1
C  THIS IS CHANGED
C     ISI=1-2*MOD(M,2)
C
      K0 = (L+M+1)/2
      K = L
      IRE = 1
      IC1 = IC
      DO I1 = 1,IFMX
         IE(I1,IRE) = JM0(I1)
      END DO
 300  CONTINUE
      IF ( (K-1).LT.K0 ) THEN
         IC2 = IC
         DO I1 = 1,IFMX
            IED(I1) = IE(I1,1)
            DO IR = 2,IRE
               IF ( IE(I1,IR).LT.IED(I1) ) IED(I1) = IE(I1,IR)
            END DO
            DO IR = 1,IRE
               IE(I1,IR) = IE(I1,IR) - IED(I1)
            END DO
         END DO
         IR = 0
         DO IC = IC1,IC2
            IR = IR + 1
            CL(IC) = 1.D0
            DO I1 = 1,IFMX
               CL(IC) = CL(IC)*IFI(I1)**IE(I1,IR)
            END DO
            CL(IC) = ISI*CL(IC)
            ISI = -ISI
         END DO
         IF ( M.EQ.0 ) IL2P(1) = IL2P(1) - 1
         UP = 1.D0
         UPSQ = 1.D0
         DOWN = 1.D0
         DO I1 = 1,IFMX
            IEUPSQ = 2*IED(I1) + IL2P(I1) + L1(I1)
            IEINT = IEUPSQ/2 - L2(I1)
            IEMOD = MOD(IEUPSQ,2)
            UPSQ = UPSQ*IFI(I1)**IEMOD
            IF ( IEINT.GE.0 ) THEN
               UP = UP*IFI(I1)**IEINT
            ELSE
               DOWN = DOWN*IFI(I1)**(-IEINT)
            END IF
         END DO
         COE(ICE) = SQRT(UPSQ)*UP/DOWN
         IF ( ITESTSFN.NE.0 ) WRITE (6,99001) L,M,UP,UPSQ,DOWN,
     &                               (CL(IC),IC=IC1,IC2)
         IF ( M.EQ.0 ) THEN
            IF ( ITESTSFN.NE.0 ) WRITE (6,99002)
            IF ( L.NE.LMAX ) THEN
               LA = (2*L+1)*(2*L+2)
               LB = (L+1)*2
               CALL SFNREDUCE(LA,IFMX,IFI,IEA)
               CALL SFNREDUCE(LB,IFMX,IFI,IEB)
               DO I1 = 1,IFMX
                  L1ST(I1) = L1ST(I1) + IEA(I1)
                  L2ST(I1) = L2ST(I1) + IEB(I1)
               END DO
               L = L + 1
               GOTO 100
            END IF
         ELSE
            LA = L + M
            LB = L - M + 1
            CALL SFNREDUCE(LA,IFMX,IFI,IEA)
            CALL SFNREDUCE(LB,IFMX,IFI,IEB)
            DO I1 = 1,IFMX
               JM0(I1) = JM0(I1) + IEA(I1) - IEB(I1)
               L1(I1) = L1(I1) - IEA(I1) + IEB(I1)
            END DO
            M = M - 1
            GOTO 200
         END IF
      ELSE
         IRE = IRE + 1
         IC = IC + 1
         LA = (2*K-L-M)*(2*K-L-M-1)
         LB = 2*(2*K-1)*(L-K+1)
         CALL SFNREDUCE(LA,IFMX,IFI,IEA)
         CALL SFNREDUCE(LB,IFMX,IFI,IEB)
         DO I1 = 1,IFMX
            IE(I1,IRE) = IE(I1,IRE-1) + IEA(I1) - IEB(I1)
         END DO
         K = K - 1
         GOTO 300
      END IF
      RETURN
99001 FORMAT (2X,'L=',I2,' M=',I2,F10.3,' *SQRT(',F16.2,')/',F10.3/2X,
     &        'CL  :',6F14.2)
99002 FORMAT (80('*'))
99003 FORMAT (13X,'THERE ARE',I5,'  COEFFICIENTS'/)
99004 FORMAT (13X,'FROM CCOEF: INCONSISTENCY DATA-DIMENSION'/14X,
     &        'LMAX:',2I5/13X,'ICMAX:',2I5)
      END
C*==sfndrotreal.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNDROTREAL(LMAX,ALPHA,BETA,GAMMA,LSFMAX,ISUMD,DMATL)
C------------------------------------------------------------------
C     THIS ROUTINE COMPUTES TRANSFORMATION MATRICES ASSOCIATED TO
C     THE ROTATION THROUGH THE EULER ANGLES ALPHA,BETA,GAMMA  FOR
C     REAL SPHERICAL HARMONICS UP TO QUANTUM NUMBER LMAX.
C------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,BETA,GAMMA
      INTEGER ISUMD,LMAX,LSFMAX
      REAL*8 DMATL(ISUMD)
C
C Local variables
C
      REAL*8 D,D1,D2,DMN(:,:),DPL(:,:),FAC,FAC1,FAC2,SQR2
      INTEGER I,IMAX,IP,IPMAX,ISU,L,M,MP
      REAL*8 SFNDROT
      EXTERNAL SFNDROT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMN,DPL
C
      ALLOCATE (DMN(LSFMAX+1,LSFMAX+1),DPL(LSFMAX+1,LSFMAX+1))
C
      IF ( LSFMAX.LE.0 ) WRITE (6,*) 'DUMMY',LSFMAX,ISUMD
C
      SQR2 = SQRT(2.D0)
      ISU = 0
      DO L = 0,LMAX
         FAC2 = 1.0D0
         M = 0
 50      CONTINUE
         FAC1 = 1.0D0
         MP = 0
 100     CONTINUE
         FAC = FAC1*FAC2/2.0D0
         D1 = SFNDROT(L,MP,M,BETA)
         D2 = SFNDROT(L,MP,-M,BETA)
         IF ( MOD(M,2).NE.0 ) D2 = -D2
         DPL(MP+1,M+1) = (D1+D2)*FAC
         DMN(MP+1,M+1) = (D1-D2)*FAC
         IF ( MOD(M+MP,2).NE.0 ) THEN
            DMN(M+1,MP+1) = -DMN(MP+1,M+1)
            DPL(M+1,MP+1) = -DPL(MP+1,M+1)
         ELSE
            DPL(M+1,MP+1) = DPL(MP+1,M+1)
            DMN(M+1,MP+1) = DMN(MP+1,M+1)
         END IF
         FAC1 = SQR2
         MP = 1 + MP
         IF ( MP.LE.M ) GOTO 100
         FAC2 = SQR2
         M = 1 + M
         IF ( M.LE.L ) GOTO 50
         IMAX = 1
         M = 0
 150     CONTINUE
         DO I = 1,IMAX
            IPMAX = 1
            MP = 0
 160        CONTINUE
            DO IP = 1,IPMAX
               IF ( I.EQ.2 ) THEN
                  IF ( IP.EQ.2 ) THEN
                     D = -SIN(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     &                   + COS(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                  ELSE
                     D = -COS(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     &                   - SIN(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                  END IF
               ELSE IF ( IP.EQ.2 ) THEN
                  D = SIN(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     &                + COS(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
               ELSE
                  D = COS(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     &                - SIN(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
               END IF
C THIS IS CHANGED
               IF ( MOD(M+MP,2).NE.0 ) D = -D
C
               ISU = ISU + 1
               DMATL(ISU) = D
            END DO
            IPMAX = 2
            MP = MP + 1
            IF ( MP.LE.L ) GOTO 160
         END DO
         IMAX = 2
         M = M + 1
         IF ( M.LE.L ) GOTO 150
      END DO
C
      DEALLOCATE (DMN,DPL)
C
      END
C*==sfndrot.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      FUNCTION SFNDROT(L,MP,M,BETA)
C-----------------------------------------------------------------------
C     CALCULATION OF D COEFFICIENT ACCORDING TO ROSE, ELEMENTARY THEORY
C     ANGULAR MOMENTUM,J.WILEY & SONS ,1957 , EQ. (4.13).
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BETA
      INTEGER L,M,MP
      REAL*8 SFNDROT
C
C Local variables
C
      REAL*8 BETA2,COSB,FF,SINB,TERM
      INTEGER I,K,KMAX,KMIN,L0,LTRM,M0,MP0,N,N1,N2,N3,N4,NF(4),NN
C
C*** End of declarations rewritten by SPAG
C
      EQUIVALENCE (N1,NF(1))
      EQUIVALENCE (N2,NF(2))
      EQUIVALENCE (N3,NF(3))
      EQUIVALENCE (N4,NF(4))
C
      DATA L0,M0,MP0/ - 1,1,1/
      L0 = -1
      M0 = 1
      MP0 = 1
      IF ( L.EQ.L0 ) THEN
         WRITE (6,99001) L0,M0,MP0
         IF ( (IABS(M).EQ.IABS(M0) .AND. IABS(MP).EQ.IABS(MP0)) .OR. 
     &        (IABS(MP).EQ.IABS(M0) .AND. IABS(M).EQ.IABS(MP0)) )
     &        GOTO 100
      END IF
      FF = 1.0D0
      IF ( IABS(M).LE.L .AND. IABS(MP).LE.L ) THEN
         N1 = L + M
         N2 = L - M
         N3 = L + MP
         N4 = L - MP
         L0 = L
         M0 = M
         MP0 = MP
         DO N = 1,4
            NN = NF(N)
            IF ( NN.NE.0 ) THEN
               DO I = 1,NN
                  FF = FF*I
               END DO
            END IF
         END DO
         FF = SQRT(FF)
      ELSE
         WRITE (6,99002) L,M,MP
         RETURN
      END IF
 100  CONTINUE
      BETA2 = BETA/2.0D0
      COSB = COS(BETA2)
      SINB = -SIN(BETA2)
      IF ( ABS(COSB).LT.1.0D-4 ) THEN
         LTRM = L
         TERM = FF
         IF ( SINB.LT.0.0D0 .AND. MOD(MP-M,2).NE.0 ) TERM = -TERM
      ELSE IF ( ABS(SINB).LT.1.0D-4 ) THEN
         LTRM = 0
         TERM = FF
         IF ( COSB.LT.0.0D0 .AND. MOD(MP-M,2).NE.0 ) TERM = -TERM
      ELSE
         KMAX = MIN0(L-MP,L+M)
         KMIN = MAX0(M-MP,0)
         TERM = COSB**(2*L+M-MP-2*KMIN)*SINB**(MP-M+2*KMIN)*FF
         GOTO 200
      END IF
      KMAX = M - MP
      IF ( MOD(KMAX,2).NE.0 ) THEN
         SFNDROT = 0.0D0
         GOTO 99999
      ELSE
         KMAX = LTRM + KMAX/2
         IF ( KMAX.LT.MAX0(M-MP,0) ) THEN
            SFNDROT = 0.0D0
            GOTO 99999
         ELSE IF ( KMAX.GT.MIN0(L-MP,L+M) ) THEN
            SFNDROT = 0.0D0
            GOTO 99999
         ELSE
            KMIN = KMAX
         END IF
      END IF
 200  CONTINUE
      IF ( MOD(KMIN,2).NE.0 ) TERM = -TERM
      N1 = L - MP - KMIN
      N2 = L + M - KMIN
      N3 = KMIN + MP - M
      N4 = KMIN
      DO N = 1,4
         NN = NF(N)
         IF ( NN.NE.0 ) THEN
            DO I = 1,NN
               TERM = TERM/I
            END DO
         END IF
      END DO
      SFNDROT = TERM
      IF ( KMIN.EQ.KMAX ) RETURN
      KMIN = KMIN + 1
      COSB = COSB**2
      SINB = SINB**2
      N3 = N3 + 1
      DO K = KMIN,KMAX
         TERM = -N1*N2*TERM*SINB/(COSB*K*N3)
         SFNDROT = SFNDROT + TERM
         N1 = N1 - 1
         N2 = N2 - 1
         N3 = N3 + 1
      END DO
      RETURN
99001 FORMAT (3X,'L0,M0,MP0=',3I3)
99002 FORMAT ('     L=',I5,'    M=',I5,'    MP =',I5)
99999 CONTINUE
      END
C*==sfnintsim.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNINTSIM(FD,RATIO,X1,X2,DLT,XEL)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DLT,FD,RATIO,X1,X2
      REAL*8 XEL(5)
C
C Local variables
C
      REAL*8 FFF
      REAL*8 H,X
      INTEGER I,INB,INC,K,N
      EXTERNAL FFF
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,5
         XEL(I) = 0.0D0
      END DO
      N = INT((X2-X1)/DLT) + 1
      INC = 2*N
      INB = INC - 1
      H = (X2-X1)/DBLE(INC)
      DO K = 1,INB,2
         X = DBLE(K)*H + X1
         DO I = 1,5
            XEL(I) = XEL(I) + 4.0D0*FFF(I,X,FD,RATIO)
         END DO
      END DO
      DO K = 2,INB,2
         X = DBLE(K)*H + X1
         DO I = 1,5
            XEL(I) = XEL(I) + 2.0D0*FFF(I,X,FD,RATIO)
         END DO
      END DO
      DO I = 1,5
         XEL(I) = XEL(I) + FFF(I,X1,FD,RATIO) + FFF(I,X2,FD,RATIO)
      END DO
      DO I = 1,5
         XEL(I) = H*XEL(I)/3.0D0
      END DO
      END
C*==fff.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      REAL*8 FUNCTION FFF(I,X,FD,RATIO)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 FD,RATIO,X
      INTEGER I
C
C Local variables
C
      REAL*8 A,A1,A2,AEXP,B,B1,C,C1,C2
C
C*** End of declarations rewritten by SPAG
C
      AEXP = 3.0D0/2.0D0
      A1 = RATIO*COS(X-FD)*(1.0D0-RATIO*RATIO)
      A2 = (1.0D0-RATIO*RATIO*SIN(X-FD)*SIN(X-FD))**AEXP
      A = A1/A2
      C1 = RATIO*COS(X-FD)
      C2 = RATIO*RATIO*(1+2.0D0*SIN(X-FD)*SIN(X-FD)) - 3.0D0
      C = 2.0D0 + C1*C2/A2
      B1 = (1.0D0-RATIO*RATIO)**AEXP
      B = B1/A2
      IF ( I.EQ.1 ) FFF = SQRT(15.0D0)*SIN(2.0D0*X)*C/6.0D0
      IF ( I.EQ.2 ) FFF = -SQRT(5.0D0/3.0D0)*SIN(X)*B
      IF ( I.EQ.3 ) FFF = SQRT(5.0D0)*A/2.0D0
      IF ( I.EQ.4 ) FFF = -SQRT(5.0D0/3.0D0)*COS(X)*B
      IF ( I.EQ.5 ) FFF = SQRT(15.0D0)*COS(2.0D0*X)*C/6.0D0
      END
C*==sfnmesh0.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNMESH0(KEYPAN,CRT,NM,NPAN,NRSF0,NPANMAX)
C ***********************************************************
C *  THIS SUBROUTINE CALCULATES AN APPROPRIATE MESH FOR
C *  THE SHAPE FUNCTIONS. MORE THAN NMIN POINTS BETWEEN TWO
C *  CRITICAL POINTS
C *  In case of more dense mesh increase NMIN
C *
C *  force the number of mesh points to be odd
C *  to apply Simpson rule without exceptions
C *
C *
C *  12/02/20  reset to original values
C *            NMIN = 7    KEYPAN = 0
C *
C ***********************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNMESH0')
      INTEGER NMIN_TOP,NMIN_BOT
      PARAMETER (NMIN_TOP=7,NMIN_BOT=3)
C
C Dummy arguments
C
      INTEGER KEYPAN,NPAN,NPANMAX,NRSF0
      REAL*8 CRT(NPANMAX)
      INTEGER NM(NPANMAX)
C
C Local variables
C
      REAL*8 D1,DIST,DPAN,DX0,R1,RA,RB
      INTEGER I,IRA,N,NA,NMIN,NTOT
C
C*** End of declarations rewritten by SPAG
C
      DO NMIN = NMIN_TOP,NMIN_BOT, - 2
         IF ( (NPAN-1)*NMIN.LE.NRSF0 ) GOTO 100
      END DO
C
      WRITE (6,*) '<SFNMESH0>  issue with radial grid setup:'
      WRITE (6,*) 'NPAN, NMIN_TOP, NMIN_BOT =',NPAN,NMIN_TOP,NMIN_BOT
      WRITE (6,*) '(NPAN-1)*NMIN_BOT =',(NPAN-1)*NMIN_BOT
      WRITE (6,*) '            NRSF0 =',NRSF0
      CALL STOP_MESSAGE(ROUTINE,'increase number of points  NRSF0')
C
 100  CONTINUE
      KEYPAN = 0
C
C-----------------------------------------------------------------------
C                           original settings
C-----------------------------------------------------------------------
      IF ( KEYPAN.EQ.0 ) THEN
C
         DIST = ABS(CRT(1)-CRT(NPAN))
         DO I = 1,NPAN - 1
            D1 = ABS(CRT(I)-CRT(I+1))
            NM(I) = INT(NRSF0*D1/DIST)
Ccc      IF ( MOD(NM(I),2) .EQ. 0 ) NM(I) = NM(I) + 1
            IF ( NM(I).LT.NMIN ) NM(I) = NMIN
         END DO
         N = NMIN*(NPAN-1)
         N = NRSF0 - N
         IF ( N.LE.0 ) CALL STOP_MESSAGE(ROUTINE,
     &        'INCREACE NUMBER OF MESH POINTS')
         D1 = DBLE(N)/DBLE(NRSF0)
         NTOT = N
         DO I = 1,NPAN - 1
            NA = NINT(D1*DBLE(NM(I)))
            IF ( NM(I).GT.NMIN .AND. (NTOT-NA).GT.0 ) THEN
               NM(I) = NMIN + NA
               NTOT = NTOT - NA
            END IF
         END DO
         NM(1) = NM(1) + NTOT
C
C-----------------------------------------------------------------------
C                    allow increasing step size
C-----------------------------------------------------------------------
      ELSE
C
         R1 = 1D-6
         DX0 = 0.021D0
C
         NTOT = 0
         DO I = 1,NPAN - 1
C
            DPAN = (CRT(I+1)-CRT(I))
            RA = CRT(I) + 0.25D0*DPAN
            IRA = INT(LOG(RA/R1)/DX0) + 1
            RA = R1*EXP(DX0*(IRA-1))
            RB = R1*EXP(DX0*IRA)
C
            NM(I) = INT(DPAN/(RB-RA)) + 1
            NM(I) = MAX(NMIN,NM(I))
C
            NTOT = NTOT + NM(I)
         END DO
C
         IF ( NTOT.GT.NRSF0 ) THEN
            WRITE (6,*) NPAN,NTOT,NRSF0
            CALL STOP_MESSAGE(ROUTINE,'increase number of points')
         END IF
C
      END IF
C
      END
C*==sfnpolchk.f    processed by SPAG 6.70Rc at 09:32 on 30 Mar 2017
      SUBROUTINE SFNPOLCHK(NFACE,NVERTICES,XEDGE,YEDGE,ZEDGE,NVERTMAX,
     &                     NFACEMAX,IFLAG_POLYHEDRON)
C     ----------------------------------------------------------------
C     THIS SUBROUTINE READS THE COORDINATES OF THE VERTICES OF EACH
C     (POLYGON)  FACE OF  A CONVEX POLYHEDRON AND  CHECKS  IF THESE
C     VERTICES ARRANGED  CONSECUTIVELY DEFINE A  POLYGON. THEN  THE
C     SUBROUTINE  DETERMINES  THE  VERTICES  AND  THE  EDGES OF THE
C     POLYHEDRON AND CHECKS IF  THE  NUMBER  OF  VERTICES  PLUS THE
C     NUMBER OF FACES EQUALS THE NUMBER OF EDGES PLUS 2.
C     ----------------------------------------------------------------
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNPOLCHK')
C
C Dummy arguments
C
      INTEGER IFLAG_POLYHEDRON,NFACE,NFACEMAX,NVERTMAX
      INTEGER NVERTICES(NFACEMAX)
      REAL*8 XEDGE(NVERTMAX,NFACEMAX),YEDGE(NVERTMAX,NFACEMAX),
     &       ZEDGE(NVERTMAX,NFACEMAX)
C
C Local variables
C
      REAL*8 A1,A2,ARG,DOWN,DV1(3),DV2(3),DVM(3),DVP(3),FISUM,UP,V(:,:),
     &       V1(:,:),V2(:,:),VRT(:,:)
      REAL*8 DDOT,DNRM2
      INTEGER IEDGE,IFACE,INEW,IVERT,IVERTM,IVERTP,IVRT,NEDGE,NVERT,
     &        NVRT,NVRTMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE V,V1,V2,VRT
C
      NVRTMAX = SUM(NVERTICES(1:NFACE))
      ALLOCATE (V(3,NVERTMAX),V1(3,NVERTMAX),V2(3,NVERTMAX))
      ALLOCATE (VRT(3,NVRTMAX))
C
      V1(:,:) = 0D0
      V2(:,:) = 0D0
      VRT(:,:) = 0D0
C
      REWIND (100)
      NVRT = 0
      NEDGE = 0
      DO IFACE = 1,NFACE
         NVERT = NVERTICES(IFACE)
         FISUM = (NVERT-2)*PI
C
         DO IVERT = 1,NVERT
            V(1,IVERT) = XEDGE(IVERT,IFACE)
            V(2,IVERT) = YEDGE(IVERT,IFACE)
            V(3,IVERT) = ZEDGE(IVERT,IFACE)
         END DO
C
C------> T R E A T M E N T   O F   V E R T I C E S
C
         WRITE (100,'(2I4,3F20.8)') IFACE
         DO IVERT = 1,NVERT
C
            WRITE (100,'(2I4,3F20.8)') IFACE,IVERT,V(1:3,IVERT)
C
            INEW = 1
C           Save all different vertices
            DO IVRT = 1,NVRT
               DV1(1:3) = V(1:3,IVERT) - VRT(1:3,IVRT)
               IF ( DDOT(3,DV1,1,DV1,1).LT.1.D-10 ) THEN
                  INEW = 0
                  EXIT
               END IF
            END DO
            IF ( INEW.EQ.1 ) THEN
               NVRT = NVRT + 1
               VRT(1:3,NVRT) = V(1:3,IVERT)
            END IF
C
            IVERTP = IVERT + 1
            IF ( IVERT.EQ.NVERT ) IVERTP = 1
            DVP(1:3) = V(1:3,IVERTP) - V(1:3,IVERT)
C
            IVERTM = IVERT - 1
            IF ( IVERT.EQ.1 ) IVERTM = NVERT
            DVM(1:3) = V(1:3,IVERTM) - V(1:3,IVERT)
C
C     Check if the  consecutive vertices define a polygon
C
            A1 = DNRM2(3,DVP,1)
            A2 = DNRM2(3,DVM,1)
            DOWN = A1*A2
            UP = DDOT(3,DVP,1,DVM,1)
C
C           IF ( DOWN.GT.1.D-6 ) THEN
            IF ( DOWN.GT.1.D-10 ) THEN
C
               ARG = UP/DOWN
               IF ( ABS(ARG).GE.1.D0 ) ARG = SIGN(1.D0,ARG)
               FISUM = FISUM - ACOS(ARG)
C
            ELSE
               WRITE (6,*) ' '
               WRITE (6,*) 'IFACE ',IFACE
               WRITE (6,*) 'IVERT ',IVERT
               WRITE (6,'(A,3F14.8)') 'V - ',V(1:3,IVERTM)
               WRITE (6,'(A,3F14.8)') 'V   ',V(1:3,IVERT)
               WRITE (6,'(A,3F14.8)') 'V + ',V(1:3,IVERTP)
               WRITE (6,'(A,3F16.10)') 'A1  ',A1
               WRITE (6,'(A,3F16.10)') 'A2  ',A2
               WRITE (6,'(A,3F16.10)') 'DOWN',DOWN
               WRITE (6,'(A,3F16.10)') 'UP  ',UP
C
               CALL STOP_MESSAGE(ROUTINE,
     &                           'IDENTICAL CONSECUTIVE VERTICES')
            END IF
C
C------> T R E A T M E N T   O F   E D G E S
C
C---------------------------------------------- Save all different edges
C
            INEW = 1
            DO IEDGE = 1,NEDGE
               DV1(1:3) = V(1:3,IVERT) - V1(1:3,IEDGE)
               IF ( DDOT(3,DV1,1,DV1,1).LT.1.D-10 ) THEN
                  DV2(1:3) = V(1:3,IVERTP) - V2(1:3,IEDGE)
                  IF ( DDOT(3,DV2,1,DV2,1).LT.1.D-10 ) THEN
                     INEW = 0
                     EXIT
                  END IF
               ELSE
                  DV1(1:3) = V(1:3,IVERT) - V2(1:3,IEDGE)
                  IF ( DDOT(3,DV1,1,DV1,1).LT.1.D-10 ) THEN
                     DV2(1:3) = V(1:3,IVERTP) - V1(1:3,IEDGE)
                     IF ( DDOT(3,DV2,1,DV2,1).LT.1.D-10 ) THEN
                        INEW = 0
                        EXIT
                     END IF
                  END IF
               END IF
            END DO
C
            IF ( INEW.EQ.1 ) THEN
               NEDGE = NEDGE + 1
               IF ( NEDGE.GT.NVERTMAX )
     &               CALL STOP_MESSAGE(ROUTINE,'INSUFFICIENT NEDGED')
               V1(1:3,NEDGE) = V(1:3,IVERT)
               V2(1:3,NEDGE) = V(1:3,IVERTP)
            END IF
C
         END DO
C
         IF ( FISUM.GT.1.D-6 ) THEN
            WRITE (6,*) 'FISUM =',FISUM
            CALL STOP_MESSAGE(ROUTINE,
     &                        'NOT CONSECUTIVE VERTICES OF A POLYGON')
         END IF
      END DO
      IF ( (NVRT+NFACE).NE.(NEDGE+2) ) THEN
         WRITE (6,99001) NVRT,NFACE,NEDGE
         IFLAG_POLYHEDRON = 1
      END IF
C
      DEALLOCATE (V,V1,V2,VRT)
99001 FORMAT (/,10X,'WARNING from <SFNPOLCHK>: illegal polyhedron',/,
     &        10X,'NVRT+NFACE .NE. NEDGE+2',/,10X,'NVRT  =',I4,/,10X,
     &        'NFACE =',I4,/,10X,'NEDGE =',I4,/)
C
      END
