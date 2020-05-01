C*==fpvmuftin.f    processed by SPAG 6.70Rc at 21:36 on 19 Dec 2016
      SUBROUTINE FPVMUFTIN(NLMFP,V,VMTZ)
C   ********************************************************************
C   *                                                                  *
C   * determine muffin tin zero and shift potential to muffin tin zero *
C   *                                                                  *
C   *       empty spheres are excluded form the AVERAGING              *
C   *                                                                  *
C   * 08/07/14  VMTZ: sign convention as in ASA                        *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NLMFPMAX,ITBOT,ITTOP,CONC,NAT,NTMAX,IMT,Z,
     &    KLMFP,NLMFPT
      USE MOD_RMESH,ONLY:R,JRMT,JRCUT,NPAN,FLMSF,KLMSF,ISFLM,NRMAX,
     &    DRDI_W_RADINT,N_VORONOI_SURF_GRID,D_VORONOI_SURF_GRID,
     &    Y_VORONOI_SURF_GRID,IRLAG_VORONOI_SURF_GRID,
     &    W_VORONOI_SURF_GRID
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION,SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_CALCMODE,ONLY:VMTZ_AVERAGE_SURFACE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPVMUFTIN')
C
C Dummy arguments
C
      INTEGER NLMFP
      REAL*8 VMTZ
      REAL*8 V(NRMAX,NLMFPMAX,NTMAX)
C
C Local variables
C
      LOGICAL CALCVMTZ,INITIALIZE
      REAL*8 DDOT
      REAL*8 DR,DRSQ,V_GRID,V_LM_LAG,WLAG(0:2),Y_LM
      INTEGER IM,IR,IR1,IRCRIT,IRLAG,IRMTIN,IRSF,ISF,IT,I_GRID,J,LM,NR
      DOUBLE PRECISION RWGT,SUM_V_WGHT,SUM_WGHT,V1(NRMAX),V2(NRMAX),VAV,
     &                 VCORR,VOL
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
C     VMTZ_AVERAGE_SURFACE=.TRUE.
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE .AND. VMTZ_AVERAGE_SURFACE ) THEN
C
         CALL VORONOI_SURF_INTEG_INIT
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' .OR. SYSTEM_TYPE(1:3)
     &     .EQ.'VIV' ) THEN
         CALCVMTZ = .TRUE.
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
         CALCVMTZ = .FALSE.
      ELSE IF ( SUB_SYSTEM(1:6).EQ.'L-BULK' ) THEN
         CALCVMTZ = .TRUE.
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'R-BULK not yet implemented')
      END IF
C
      IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) CALCVMTZ = .FALSE.
C
C=======================================================================
      IF ( CALCVMTZ ) THEN
C
         SUM_V_WGHT = 0.D0
         SUM_WGHT = 0.D0
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         DO IT = ITBOT,ITTOP
C
C---------------------------------- EXCLUDE empty spheres from averaging
            IF ( Z(IT).EQ.0 ) CYCLE
C
            IM = IMT(IT)
C
CVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C      take the average over the VOLUME of the interstitial regime
CVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            IF ( .NOT.VMTZ_AVERAGE_SURFACE ) THEN
C
               IRMTIN = JRMT(IM)
               IRCRIT = JRCUT(NPAN(IM),IM)
C
               DO IR = IRMTIN + 1,IRCRIT
                  RWGT = R(IR,IM)**2*FLMSF(IR-IRMTIN,1,IM)
                  V1(IR) = V(IR,1,IT)*RWGT
                  V2(IR) = RWGT*SQRT_4PI
               END DO
C
               DO LM = 2,NLMFP
                  IF ( KLMSF(LM,IM).GT.0 ) THEN
                     ISF = ISFLM(LM,IM)
                     DO IR = IRMTIN + 1,IRCRIT
                        IRSF = IR - IRMTIN
                        V1(IR) = V1(IR) + R(IR,IM)**2*V(IR,LM,IT)
     &                           *FLMSF(IRSF,ISF,IM)
                     END DO
                  END IF
C
               END DO
C
               NR = IRCRIT - IRMTIN
               IR1 = IRMTIN + 1
               VAV = DDOT(NR,V1(IR1),1,DRDI_W_RADINT(IR1,IM),1)
               VOL = DDOT(NR,V2(IR1),1,DRDI_W_RADINT(IR1,IM),1)
C
CVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            ELSE
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C           take the average over the SURFACE of the polyhedron
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
               VOL = 0.0D0
               VAV = 0.0D0
C
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
               DO I_GRID = 1,N_VORONOI_SURF_GRID(IM)
C
                  IRLAG = IRLAG_VORONOI_SURF_GRID(I_GRID,IM)
C
                  DR = D_VORONOI_SURF_GRID(I_GRID,IM) - R(IRLAG,IM)
                  DRSQ = DR*DR
                  WLAG(0) = 0.5D0*(2-3*DR+DRSQ)
                  WLAG(1) = 2*DR - DRSQ
                  WLAG(2) = 0.5D0*(-DR+DRSQ)
C
                  IF ( ABS(1D0-SUM(WLAG(0:2))).GT.1D-10 )
     &                 CALL STOP_MESSAGE(ROUTINE,'SUM(f) w(f) != 1')
C
                  V_GRID = 0.0D0
C
                  DO LM = 1,NLMFPT(IT)
                     IF ( KLMFP(LM,IT).NE.0 .AND. LM.LE.NLMFPMAX ) THEN
C
                        V_LM_LAG = 0D0
                        DO J = 0,2
                           V_LM_LAG = V_LM_LAG + WLAG(J)
     &                                *V(IRLAG+J,LM,IT)
                        END DO
C
                        Y_LM = Y_VORONOI_SURF_GRID(LM,I_GRID,IM)
C
                        V_GRID = V_GRID + V_LM_LAG*Y_LM
C
                     END IF
                  END DO
C
                  VOL = VOL + W_VORONOI_SURF_GRID(I_GRID,IM)
                  VAV = VAV + W_VORONOI_SURF_GRID(I_GRID,IM)*V_GRID
C
               END DO
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
            END IF
C
            SUM_V_WGHT = SUM_V_WGHT + VAV*CONC(IT)*NAT(IT)
            SUM_WGHT = SUM_WGHT + VOL*CONC(IT)*NAT(IT)
C
         END DO
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
         VMTZ = SUM_V_WGHT/SUM_WGHT
C
         WRITE (6,99001) VMTZ
C
      ELSE
C
         WRITE (6,99002) VMTZ
C
      END IF
C=======================================================================
C
C-----------------------------------------------------------------------
C     shift potential to muffin tin zero for ALL sites
C-----------------------------------------------------------------------
      VCORR = SQRT_4PI*VMTZ
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         DO IR = 1,JRCUT(NPAN(IM),IM)
            V(IR,1,IT) = V(IR,1,IT) - VCORR
         END DO
      END DO
C
99001 FORMAT (/,1X,79('*'),/,10X,'shift of muffin-tin zero  VMTZ ',
     &        F10.6,/,1X,79('*'),/)
99002 FORMAT (/,1X,79('*'),/,10X,'shift of muffin-tin zero  VMTZ ',
     &        F10.6,'     fixed by host'/,1X,79('*'),/)
      END
C*==voronoi_surf_integ_init.f    processed by SPAG 6.70Rc at 21:36 on 19 Dec 2016
      SUBROUTINE VORONOI_SURF_INTEG_INIT
C   ********************************************************************
C   *                                                                  *
C   *  find the grid points next to the intersection of the            *
C   *  directions and the polyhedron surface used for interpolation    *
C   *                                                                  *
C   *  tabulate real spherical harmonics for directions                *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NLMFPMAX,NLFPMAX
      USE MOD_RMESH,ONLY:NM,NMMAX,JRCUT,NPAN,R,NMAX_VORONOI_SURF_GRID,
     &    N_VORONOI_SURF_GRID,D_VORONOI_SURF_GRID,Y_VORONOI_SURF_GRID,
     &    RHAT_VORONOI_SURF_GRID,IRLAG_VORONOI_SURF_GRID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='VORONOI_SURF_INTEG_INIT')
      INTEGER NLAG
      PARAMETER (NLAG=3)
C
C Local variables
C
      REAL*8 D,RHAT(3)
      INTEGER IA_ERR,IM,IPAN,IR,IRBOT_PAN,IRLAG,I_GRID,N
C
C*** End of declarations rewritten by SPAG
C
      N = NMAX_VORONOI_SURF_GRID
      ALLOCATE (IRLAG_VORONOI_SURF_GRID(N,NMMAX))
      ALLOCATE (Y_VORONOI_SURF_GRID(NLMFPMAX,N,NMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 )
     &      CALL STOP_MESSAGE(ROUTINE,'ALLOC: Y_VORONOI_SURF_GRID')
C-----------------------------------------------------------------------
C
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      DO IM = 1,NM
C
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
         DO I_GRID = 1,N_VORONOI_SURF_GRID(IM)
C
            D = D_VORONOI_SURF_GRID(I_GRID,IM)
            IRLAG = 0
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
            DO IPAN = 2,NPAN(IM)
C
               IRBOT_PAN = JRCUT(IPAN-1,IM) + 1
C
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
               DO IR = IRBOT_PAN,JRCUT(IPAN,IM)
C
                  IF ( R(IR,IM).GE.D .AND. IRLAG.EQ.0 ) THEN
                     IRLAG = MAX(IRBOT_PAN,IR-(NLAG-1))
                     IRLAG_VORONOI_SURF_GRID(I_GRID,IM) = IRLAG
                     IF ( IRLAG+(NLAG-1).GT.JRCUT(IPAN,IM) )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                    'IRLAG+(NLAG-1) > JRCUT(IPAN,IM)')
                     EXIT
                  END IF
C
               END DO
Crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
               IF ( IRLAG.NE.9999 ) EXIT
            END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C
            RHAT(1:3) = RHAT_VORONOI_SURF_GRID(1:3,I_GRID,IM)
C
            CALL CALC_RHPLM(RHAT(1),RHAT(2),RHAT(3),
     &                      Y_VORONOI_SURF_GRID(1,I_GRID,IM),NLFPMAX-1,
     &                      NLMFPMAX)
C
         END DO
Cddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
C
      END DO
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
C
      WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),NLFPMAX - 1,
     &                NMAX_VORONOI_SURF_GRID
C
      RETURN
C-----------------------------------------------------------------------
99001 FORMAT (/,1X,79('*'),/,35X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'spherical harmonics created for l_max =',I3,/,10X,
     &        'on a triangular integration mesh with max. ',I5,
     &        '  points',/,10X,'on the surface of an atomic polyhedron',
     &        /)
      END
