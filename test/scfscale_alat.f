C*==scfscale_alat.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFSCALE_ALAT(RSCL,QAVAILABLE)
C   ********************************************************************
C   *                                                                  *
C   *   scale the radial mesh by the factor    RSCL                    *
C   *                                                                  *
C   *   assuming ASA and exponential mesh (so far)                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NM,DX,R,RMESHTYPE,RWS,RMT,JRWS,JRMT,NRMAX,NMMAX
      USE MOD_TYPES,ONLY:NT,VT,BT,CONC,IMT,NAT
      USE MOD_SITES,ONLY:NQ
      USE MOD_LATTICE,ONLY:ALAT,ABAS,SWS
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--SCFSCALE_ALAT16
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL QAVAILABLE
      REAL*8 RSCL
C
C Local variables
C
      REAL*8 ALAT_REF,B0(NRMAX),DX_REF(:),R0(NRMAX,NMMAX),RI,RMT_REF(:),
     &       RWS_REF(:),R_REF(:,:),SWS_REF,VR0(NRMAX),VUC,VWS
      INTEGER I,IERR,IM,IT,JRMT_REF(:),JRWS_REF(:),JTOP
      LOGICAL INITIALIZE
      REAL*8 YLAG
      SAVE ALAT_REF,DX_REF,JRMT_REF,JRWS_REF,RMT_REF,RWS_REF,R_REF,
     &     SWS_REF
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE R_REF,RWS_REF,RMT_REF,JRWS_REF,JRMT_REF,DX_REF
C
      DATA INITIALIZE/.TRUE./
C*** End of declarations rewritten by SPAG
C
C============================================================ INITIALIZE
      IF ( INITIALIZE ) THEN
C
         ALAT_REF = ALAT
         SWS_REF = SWS
C
         ALLOCATE (RWS_REF(NMMAX),RMT_REF(NMMAX),R_REF(1,NMMAX))
         ALLOCATE (JRWS_REF(NMMAX),JRMT_REF(NMMAX),DX_REF(NMMAX))
C
         DX_REF(1:NM) = DX(1:NM)
         R_REF(1,1:NM) = R(1,1:NM)
         RWS_REF(1:NM) = RWS(1:NM)
         RMT_REF(1:NM) = RMT(1:NM)
         JRWS_REF(1:NM) = JRWS(1:NM)
         JRMT_REF(1:NM) = JRMT(1:NM)
C
         INITIALIZE = .FALSE.
      END IF
C============================================================ INITIALIZE
C
      WRITE (6,99001) RSCL
C
      IERR = 0
C
C ----------------------------------------------------------------------
C
      IF ( RMESHTYPE.NE.'EXPONENTIAL ' ) THEN
         WRITE (6,*) 'radial mesh should be EXPONENTIAL !!'
         WRITE (6,*) 'no scaling done '
         RETURN
      END IF
C
      ALAT = ALAT_REF*RSCL
C
      DO IM = 1,NM
         DO I = 1,NRMAX
            R0(I,IM) = R(I,IM)
         END DO
C
C------------------------------------------------  scale ALL mesh points
         DX(IM) = DX_REF(IM)
         R(1,IM) = R_REF(1,IM)*RSCL
         RWS(IM) = RWS_REF(IM)*RSCL
         RMT(IM) = RMT_REF(IM)*RSCL
         JRWS(IM) = JRWS_REF(IM)
         JRMT(IM) = JRMT_REF(IM)
C
C--------------------------------------------  keep 1st mesh point fixed
         R(1,IM) = R_REF(1,IM)
         RWS(IM) = RWS_REF(IM)*RSCL
         RMT(IM) = RMT_REF(IM)*RSCL
         JRWS(IM) = JRWS_REF(IM)
C
         DX(IM) = LOG(RWS(IM)/R(1,IM))/DBLE(JRWS(IM)-1)
         JRMT(IM) = INT(LOG(RMT(IM)/R(1,IM))/DX(IM)) + 2
C
      END DO
C
C --------------------------------------------------------- check volume
C
      SWS = 0.0D0
      DO IT = 1,NT
         SWS = SWS + CONC(IT)*NAT(IT)*RWS(IMT(IT))**3
      END DO
      SWS = (SWS/DBLE(NQ))**(1.0D0/3.0D0)
      VWS = SWS**3*4D0*PI/3D0
C
      CALL RVECSPAT(ABAS(1,1),ABAS(1,2),ABAS(1,3),VUC,1)
C
      VUC = ABS(VUC)*ALAT**3
C
      IF ( ABS(1D0-NQ*VWS/VUC).GT.1D-8 ) THEN
         WRITE (6,99002) VUC,NQ*VWS,ABS(VUC-NQ*VWS)
         IERR = 1
      END IF
C
      WRITE (6,99005) ALAT_REF,ALAT,SWS_REF,SWS
C
C ---------------------------- print mesh info specific for atomic types
      WRITE (6,99003)
      DO IM = 1,NM
         WRITE (6,99004) IM,JRMT(IM),RMT_REF(IM),JRWS_REF(IM),
     &                   RWS_REF(IM),DX_REF(IM)
      END DO
C
C-----------------------------------------------------------------------
C                        set up new radial mesh
C-----------------------------------------------------------------------
C
      CALL RMESHASA
C
C------------------------------------------------- interpolate potential
C
      DO IT = 1, - NT
         IM = IMT(IT)
         JTOP = JRWS(IM)
C
         DO I = 1,JTOP
            VR0(I) = VT(I,IT)*R0(I,IM)
            B0(I) = BT(I,IT)
         END DO
C
         DO I = 1,JTOP
            RI = R(I,IM)
            VT(I,IT) = YLAG(RI,R0,VR0,0,3,JTOP)/RI
            BT(I,IT) = YLAG(RI,R0,B0,0,3,JTOP)
         END DO
C
      END DO
C
C=======================================================================
C
      QAVAILABLE = .FALSE.
C
      WRITE (6,99006)
C
C ======================================================================
      IF ( IERR.GT.0 ) THEN
         WRITE (6,99007)
         WRITE (6,99006)
         STOP
      END IF
C ======================================================================
99001 FORMAT (//,1X,79('*'),/,32X,'<SCFSCALE_ALAT>',/,1X,79('*'),//,10X,
     &        'scale lattice parameters by   RSCL = ',F10.6,/)
99002 FORMAT (2(/,1X,79('#')),/,10X,'WARNING from <SCFSCALE_ALAT>',/,
     &        10X,'radial mesh parameters in potfile  inconsistent '/,
     &        10X,'a1*(a2xa3) * a^3     ',F15.8,/,10X,
     &        'NQ * SWS^3 * 4*PI/3  ',F15.8,/,10X,
     &        'DEVIATION            ',F15.8,/,10X,2(/,1X,79('#')),/)
99003 FORMAT (/,10X,'radial mesh parameters: (reference values)',//,10X,
     &        'IM    JRMT              RMT      JRWS    RWS       DX')
99004 FORMAT (8X,I4,I7,10X,F10.5,I7,F10.5,F14.9)
99005 FORMAT (10X,'lattice constant  ALAT      ',F12.5,3X,'>>>',F10.5,/,
     &        10X,'average Wigner-Seitz radius ',F12.5,3X,'>>>',F10.5)
99006 FORMAT (/,1X,79('*'),/)
99007 FORMAT (10X,'inconsistency due to scaling -- program stops')
      END
