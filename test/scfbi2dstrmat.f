C*==scfbi2dstrmat.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBI2DSTRMAT(SMAT,NLMAD,NLM2MAX)
C **********************************************************************
C *                                                                    *
C *                                                                    *
C **********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_SITES,ONLY:NQMAX,QBAS,NQ_I,NQ_L
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE,ALAT,ABAS_2D,BBAS_2D,VOLUC_2D
      IMPLICIT NONE
C*--SCFBI2DSTRMAT13
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBI2DSTRMAT')
C
C Dummy arguments
C
      INTEGER NLM2MAX,NLMAD
      REAL*8 SMAT(NQMAX,NQMAX,NLM2MAX)
C
C Local variables
C
      REAL*8 DNRM2
      REAL*8 GMAX,GN2(:,:),RGNTMAD(:),RM2(:,:),RMAX,SUM2D(:),
     &       VOLUC_2D_AU,X1,X2
      INTEGER IPRINT,IQ,IRGNTMAD(:,:),ISHLD,JQ,LMAX_RGNTMAD,NGMAX,NMAXD,
     &        NRGNT123TAB(20),NRGNTMAD,NRGNTMADMAX,NRMAX,NSG(:),NSHLG,
     &        NSHLR,NSR(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA NRGNT123TAB/1,15,96,388,1181,2917,6342,12452,22525,38289,
     &     61912,95914,143531,208371,294744,407644,552931,736829,966544,
     &     1250346/
C
      ALLOCATABLE RGNTMAD,IRGNTMAD,SUM2D,GN2,RM2,NSG,NSR
C
      CALL TRACK_INFO(ROUTINE)
C
      IPRINT = 0
C
      LMAX_RGNTMAD = 2*(NLMAD-1)
C
      NRGNTMADMAX = NRGNT123TAB(NLMAD)
C
      RMAX = 7.0D0
      GMAX = 65.0D0
      VOLUC_2D_AU = VOLUC_2D*ALAT**2
C-------------------------------- fix array sizes for lattice summations
C
      X1 = DNRM2(3,ABAS_2D(1,1),1)
      X2 = DNRM2(3,ABAS_2D(1,2),1)
      ISHLD = INT(RMAX/MIN(X1,X2)+1)
C
      X1 = DNRM2(3,BBAS_2D(1,1),1)
      X2 = DNRM2(3,BBAS_2D(1,2),1)
      ISHLD = MAX(ISHLD,INT(GMAX/MIN(X1,X2)+1))
C
      ISHLD = ISHLD**2
      NMAXD = 4*ISHLD
C
      ALLOCATE (RGNTMAD(NRGNTMADMAX),IRGNTMAD(NRGNTMADMAX,3))
      ALLOCATE (SUM2D(NLM2MAX))
      ALLOCATE (GN2(2,NMAXD),RM2(2,NMAXD),NSG(ISHLD),NSR(ISHLD))
C
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(79(''=''))')
      WRITE (6,'(18X,A)') 'SCFMAD2D: setting 2D Madelung coefficients'
      WRITE (6,'(79(''=''))')
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ======================================================================
      CALL LATTICE2D(RMAX,GMAX,ALAT,ABAS_2D,BBAS_2D,NGMAX,NRMAX,NSHLG,
     &               NSHLR,NSG,NSR,GN2,RM2,IPRINT,NMAXD,ISHLD)
C ======================================================================
C
C --> calculate the gaunt coefs
C
      CALL MADELGAUNT(NLMAD,RGNTMAD,IRGNTMAD,NRGNTMAD,NRGNTMADMAX)
C
      IF ( SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV' )
     &     STOP 'NOT YET IN SCFBI2DSTRMAT>'
C
C **********************************************************************
C                          IQ   within  I-ZONE
C                          JQ   within  I-ZONE
C **********************************************************************
      DO IQ = NQ_L + 1,NQ_L + NQ_I
         DO JQ = NQ_L + 1,NQ_L + NQ_I
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
            IF ( IQ.EQ.NQ_L+1 .AND. JQ.EQ.NQ_L+1 ) THEN
               WRITE (6,'(5X,2A,/)') 
     &                '< SCFMAD2D_SMAT > : calculating 2D-lattice sums '
     &                ,'inside the slab'
               IF ( IPRINT.GE.2 ) WRITE (6,99001)
            END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C make ewald sumation in plane and inverse space
C sum if rz<>0 (out of plane)
C
            CALL SCFMAD2D_SMAT(NLMAD,ALAT,QBAS(1,IQ),QBAS(1,JQ),RM2,
     &                         NRMAX,NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,
     &                         SUM2D,VOLUC_2D_AU,LMAX_RGNTMAD,NLM2MAX)
C
            SUM2D(1) = SUM2D(1)*SQRT(4.0D0*PI)
            SMAT(IQ,JQ,1:NLM2MAX) = SUM2D(1:NLM2MAX)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
            IF ( IPRINT.GE.2 ) THEN
               WRITE (6,99002) IQ,JQ,SUM2D(1:4)
               IF ( JQ.EQ.NQ_L+NQ_I .AND. IQ.NE.NQ_L+NQ_I )
     &               WRITE (6,'(20X,20(''-''))')
            END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C
         END DO
      END DO
C
C
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC ######################################################################
CcccCIF ( SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV' ) THEN
CcccC
Cccc      NREPBASALL_L = NREPBAS_L*NQ_L
Cccc      NREPBASALL_R = NREPBAS_R*NQ_R
CcccC
CcccC **********************************************************************
CcccC                      IQ      within  I-ZONE
CcccC                      JQ      within  L-BULK basis
CcccC                      JQ_XL   within  extended L-BULK
CcccC **********************************************************************
Cccc      DO IQ = NQ_L + 1,NQ_L + NQ_I
CcccC
Cccc         JQ_XL = 0
CcccC ++++++++++++++++++++++++++++++++ loop over all sites in the left host
Cccc         DO IREPBAS_L = 1,NREPBAS_L
CcccC
Cccc            DO JQ = 1,NQ_L
CcccC
Cccc               VEC_JQ(1:3) = QBAS(1:3,JQ) - (IREPBAS_L-1)*ABAS_L(1:3,3)
CcccC
Cccc               JQ_XL = JQ_XL + 1
CcccC
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
Cccc               IF ( IQ.EQ.NQ_L+1 .AND. JQ_XL.EQ.1 ) THEN
Cccc                  WRITE (6,'(5X,2A,/)')
Cccc     &                      '< SCFMAD2D_SMAT > : calculating 2D-lattice sums '
Cccc     &                      ,'slab - left host'
Cccc                  IF ( IPRINT.GE.2 ) WRITE (6,99001)
Cccc               END IF
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC
CcccC-->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
CcccC     Inverse space sum if rz<>0 (out of plane)
CcccC
Cccc               SUM2D = 0.0D0
Cccc               CALL SCFMAD2D_SMAT(NLMAD,ALAT,QBAS(1,IQ),VEC_JQ,RM2,NRMAX,
Cccc     &                      NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SUM2D,
Cccc     &                      VOLUC_2D_AU,LMAX_RGNTMAD,NLM2MAX)
CcccC
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
Cccc               IF ( IPRINT.GE.2 ) THEN
Cccc                  WRITE (6,99002) IQ,JQ_XL,SUM2D(1)
Cccc                  IF ( JQ_XL.EQ.NREPBASALL_L .AND. IQ.NE.NQ_I )
Cccc     &                 WRITE (6,'(20X,20(''-''))')
Cccc               END IF
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC
Cccc            END DO         ! JQ loop in left host basis
Cccc         END DO            ! IREPBAS_L loop in layers to get convergence
CcccC
Cccc         IF ( JQ_XL.NE.NREPBASALL_L ) THEN
Cccc            WRITE (6,*) ' < scfmad2d > : index error ',
Cccc     &                  'JQ_XL <> NREPBAS_L*NQ_L   '
Cccc            STOP
Cccc         END IF
Cccc      END DO                    ! ILAY1 loop
CcccC **********************************************************************
CcccC
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
Cccc      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC **********************************************************************
CcccC                      IQ      within  I-ZONE
CcccC                      JQ      within  R-BULK basis
CcccC                      JQ_XR   within  extended R-BULK
CcccC **********************************************************************
Cccc      DO IQ = NQ_L + 1,NQ_L + NQ_I
CcccC
Cccc         JQ_XR = 0
Cccc         DO IREPBAS_R = 1,NREPBAS_R
Cccc            DO JQ = NQ - NQ_R + 1,NQ
CcccC
Cccc               VEC_JQ(1:3) = QBAS(1:3,JQ) + (IREPBAS_R-1)*ABAS_R(1:3,3)
CcccC
Cccc               JQ_XR = JQ_XR + 1
CcccC
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
Cccc               IF ( IQ.EQ.NQ_L+1 .AND. JQ_XR.EQ.1 ) THEN
Cccc                  WRITE (6,'(5X,2A,/)')
Cccc     &                      '< SCFMAD2D_SMAT > : calculating 2D-lattice sums '
Cccc     &                      ,'slab - right host'
Cccc                  IF ( IPRINT.GE.2 ) WRITE (6,99001)
Cccc               END IF
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC-->  make ewald sumation (in plane) and
CcccC     Inverse space sum if rz<>0 (out of plane)
CcccC
Cccc               CALL SCFMAD2D_SMAT(NLMAD,ALAT,QBAS(1,IQ),VEC_JQ,RM2,NRMAX,
Cccc     &                      NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SUM2D,
Cccc     &                      VOLUC_2D_AU,LMAX_RGNTMAD,NLM2MAX)
CcccC
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
Cccc               IF ( IPRINT.GE.2 ) THEN
Cccc                  WRITE (6,99002) IQ,JQ_XR,SUM2D(1)
Cccc                  IF ( JQ_XR.EQ.NREPBASALL_R .AND. IQ.NE.NQ_I )
Cccc     &                 WRITE (6,'(20X,20(''-''))')
Cccc               END IF
CcccC OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
CcccC
CcccC
Cccc            END DO         ! JQ loop in right host basis
Cccc         END DO            ! IREPBAS_R loop in layers to get convergence
CcccC
CcccC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Cccc         IF ( JQ_XR.NE.NREPBASALL_R ) THEN
Cccc            WRITE (6,*) ' < scfmad2d > : index error ',
Cccc     &                  'JQ_XR <> NREPBAS_R*NQ_R'
Cccc            STOP
Cccc         END IF
CcccC
Cccc      END DO                    ! IQ loop
C **********************************************************************
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C      END IF
C ######################################################################
C
99001 FORMAT (8X,'2D Lattice sum (LMXSP = 1)',/,18X,'  IQ  JQ  SUM',/,
     &        18X,23('-'))
99002 FORMAT (18X,2I5,1000D12.4)
      END
