C*==scfmad2d.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD2D(IQ_LREF,IQ_RREF,NREPBAS_L,NREPBAS_R)
C **********************************************************************
C *                                                                    *
C * calculates the Madelung potential coefficients in the 2D case      *
C * For each layer the summation for a representative site IQ_LREF     *
C * in L-bulk is split into three parts                                *
C * - within the slab                                                  *
C * - NREPBAS_L*NQ_L left  host sites                                  *
C * - NREPBAS_R*NQ_R right host sites                                  *
C *                                                                    *
C *   IQ_LREF      is the index of reference site in L-BULK            *
C *   NREPBAS_L    number of repetitions of L-BULK basis               *
C *                                                                    *
C *   analogously for R-bulk and I-zone                                *
C *                                                                    *
C * all positions must be scaled with ALAT to get them correct         *
C * (done in SCFMAD2D_SMAT)                                            *
C *                                                                    *
C **********************************************************************
      USE MOD_LATTICE,ONLY:ALAT,ABAS_L,ABAS_R,ABAS_2D,BBAS_2D,VOLUC_2D
      USE MOD_SITES,ONLY:QBAS,NQ_L,NQ_R,NQ_I,NQ,AVMAD_LI,AVMAD_LR,
     &    AVMAD_LL,AVMAD_RI,AVMAD_RR,AVMAD_RL,NLMAD,NLMMAD,AVMAD
      USE MOD_ANGMOM,ONLY:NRGNT123TAB,NLMMADMAX
      IMPLICIT NONE
C*--SCFMAD2D26
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD2D')
C
C Dummy arguments
C
      INTEGER IQ_LREF,IQ_RREF,NREPBAS_L,NREPBAS_R
C
C Local variables
C
      REAL*8 AM(:,:),BM(:),GMAX,GN2(:,:),RGNTMAD(:),RM2(:,:),RMAX,
     &       SMAT(:),VEC_JQ(3),VEC_QL(3),VEC_QR(3),VOLUC_2D_AU,X1,X2
      REAL*8 DNRM2
      INTEGER IPRINT,IQ,IQ1_I,IQ1_L,IQ1_R,IQ2_I,IQ2_L,IQ2_R,IREPBAS_L,
     &        IREPBAS_R,IRGNTMAD(:,:),ISHLD,JQ,JQ_XL,JQ_XR,LMAX_RGNTMAD,
     &        NGMAX,NLMAD_LOC,NLMMAD_LOC,NLMMAX_RGNTMAD,NMAXD,
     &        NREPBASALL_L,NREPBASALL_R,NRGNTMAD,NRGNTMADMAX,NRMAX,
     &        NSG(:),NSHLG,NSHLR,NSR(:)
      EXTERNAL MADELGAUNT,SCFMAD2D_COEF,SCFMAD2D_SMAT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RGNTMAD,IRGNTMAD,SMAT,GN2,RM2,NSG,NSR,AM,BM
C
      CALL TRACK_INFO(ROUTINE)
C
      IPRINT = 0
C
      NLMAD_LOC = MAX(4,NLMAD)
      NLMMAD_LOC = NLMAD_LOC**2
C
      WRITE (6,99003) ROUTINE(1:LEN_TRIM(ROUTINE)),NLMAD,NLMMAD,
     &                NLMAD_LOC,NLMMAD_LOC
C
      NRGNTMADMAX = 2*NRGNT123TAB(NLMAD_LOC)
C
      LMAX_RGNTMAD = 2*(NLMAD_LOC-1)
      NLMMAX_RGNTMAD = (LMAX_RGNTMAD+1)**2
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
      ALLOCATE (SMAT(NLMMAX_RGNTMAD),AM(NLMMAD_LOC,NLMMAD_LOC))
      ALLOCATE (BM(NLMMAD_LOC))
      ALLOCATE (GN2(2,NMAXD),RM2(2,NMAXD),NSG(ISHLD),NSR(ISHLD))
C
      IQ1_L = 1
      IQ2_L = NQ_L
      IQ1_I = NQ_L + 1
      IQ2_I = NQ_L + NQ_I
      IQ1_R = NQ_L + NQ_I + 1
      IQ2_R = NQ_L + NQ_I + NQ_R
C
      AVMAD(IQ1_I:IQ2_I,IQ1_I:IQ2_I,NLMMADMAX,NLMMADMAX) = 0D0
C
      ALLOCATE (AVMAD_LL(IQ1_L:IQ2_L,NLMMADMAX,NLMMADMAX))
      ALLOCATE (AVMAD_LI(IQ1_I:IQ2_I,NLMMADMAX,NLMMADMAX))
      ALLOCATE (AVMAD_LR(IQ1_R:IQ2_R,NLMMADMAX,NLMMADMAX))
C
      ALLOCATE (AVMAD_RL(IQ1_L:IQ2_L,NLMMADMAX,NLMMADMAX))
      ALLOCATE (AVMAD_RI(IQ1_I:IQ2_I,NLMMADMAX,NLMMADMAX))
      ALLOCATE (AVMAD_RR(IQ1_R:IQ2_R,NLMMADMAX,NLMMADMAX))
C
C ======================================================================
      CALL LATTICE2D(RMAX,GMAX,ALAT,ABAS_2D,BBAS_2D,NGMAX,NRMAX,NSHLG,
     &               NSHLR,NSG,NSR,GN2,RM2,IPRINT,NMAXD,ISHLD)
C ======================================================================
C
C --> calculate the Gaunt coefs
C
      CALL MADELGAUNT(NLMAD_LOC,RGNTMAD,IRGNTMAD,NRGNTMAD,NRGNTMADMAX)
C
C --> calculate the Madelung coefficients to be used for VMAD
C
C **********************************************************************
C                          calculate V_ql
C
C     IQ = IQ_LREF in L
C     AVMAD_LI: JQ in I   AVMAD_LL: JQ in L   AVMAD_LR: JQ in R
C **********************************************************************
C
      VEC_QL(1:3) = QBAS(1:3,IQ_LREF)
C
      MAD_LI:DO JQ = NQ_L + 1,NQ_L + NQ_I
C
         CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QL,QBAS(1,JQ),RM2,NRMAX,
     &                      NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                      VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
         CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                      IRGNTMAD,NRGNTMAD,NLMMAD_LOC,NLMMAX_RGNTMAD,
     &                      NRGNTMADMAX)
C
         AVMAD_LI(JQ,1:NLMMAD,1:NLMMAD) = AM(1:NLMMAD,1:NLMMAD)
C
      END DO MAD_LI
C
C
      AVMAD_LL(IQ1_L:IQ2_L,1:NLMMAD,1:NLMMAD) = 0.0D0
C
      MAD_LL:DO IREPBAS_L = 1,NREPBAS_L
         DO JQ = 1,NQ_L
C
            VEC_JQ(1:3) = QBAS(1:3,JQ) - (IREPBAS_L-1)*ABAS_L(1:3,3)
C
            CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QL,VEC_JQ,RM2,NRMAX,
     &                         NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                         VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
            CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
C
            AVMAD_LL(JQ,1:NLMMAD,1:NLMMAD)
     &         = AVMAD_LL(JQ,1:NLMMAD,1:NLMMAD) + AM(1:NLMMAD,1:NLMMAD)
C
         END DO
      END DO MAD_LL
C
C
      AVMAD_LR(IQ1_R:IQ2_R,1:NLMMAD,1:NLMMAD) = 0.0D0
C
      MAD_LR:DO IREPBAS_R = 1,NREPBAS_R
         DO JQ = NQ - NQ_R + 1,NQ
C
            VEC_JQ(1:3) = QBAS(1:3,JQ) + (IREPBAS_R-1)*ABAS_R(1:3,3)
C
            CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QL,VEC_JQ,RM2,NRMAX,
     &                         NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                         VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
            CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
C
            AVMAD_LR(JQ,1:NLMMAD,1:NLMMAD)
     &         = AVMAD_LR(JQ,1:NLMMAD,1:NLMMAD) + AM(1:NLMMAD,1:NLMMAD)
C
         END DO
      END DO MAD_LR
C
C **********************************************************************
C                          calculate V_qr
C
C     IQ = IQ_RREF in R
C     AVMAD_RI: JQ in I   AVMAD_RL: JQ in L   AVMAD_RR: JQ in R
C **********************************************************************
C
      VEC_QR(1:3) = QBAS(1:3,IQ_RREF)
C
      MAD_RI:DO JQ = NQ_L + 1,NQ_L + NQ_I
C
         CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QR,QBAS(1,JQ),RM2,NRMAX,
     &                      NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                      VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
         CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                      IRGNTMAD,NRGNTMAD,NLMMAD_LOC,NLMMAX_RGNTMAD,
     &                      NRGNTMADMAX)
C
         AVMAD_RI(JQ,1:NLMMAD,1:NLMMAD) = AM(1:NLMMAD,1:NLMMAD)
C
      END DO MAD_RI
C
C
      AVMAD_RL(IQ1_L:IQ2_L,1:NLMMAD,1:NLMMAD) = 0.0D0
C
      MAD_RL:DO IREPBAS_L = 1,NREPBAS_L
         DO JQ = 1,NQ_L
C
            VEC_JQ(1:3) = QBAS(1:3,JQ) - (IREPBAS_L-1)*ABAS_L(1:3,3)
C
            CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QR,VEC_JQ,RM2,NRMAX,
     &                         NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                         VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
            CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
C
            AVMAD_RL(JQ,1:NLMMAD,1:NLMMAD)
     &         = AVMAD_RL(JQ,1:NLMMAD,1:NLMMAD) + AM(1:NLMMAD,1:NLMMAD)
C
         END DO
      END DO MAD_RL
C
C
      AVMAD_RR(IQ1_R:IQ2_R,1:NLMMAD,1:NLMMAD) = 0.0D0
C
      MAD_RR:DO IREPBAS_R = 1,NREPBAS_R
         DO JQ = NQ - NQ_R + 1,NQ
C
            VEC_JQ(1:3) = QBAS(1:3,JQ) + (IREPBAS_R-1)*ABAS_R(1:3,3)
C
            CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,VEC_QR,VEC_JQ,RM2,NRMAX,
     &                         NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                         VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
            CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
C
            AVMAD_RR(JQ,1:NLMMAD,1:NLMMAD)
     &         = AVMAD_RR(JQ,1:NLMMAD,1:NLMMAD) + AM(1:NLMMAD,1:NLMMAD)
C
         END DO
      END DO MAD_RR
C
C **********************************************************************
C                          IQ   within  I-ZONE
C                          JQ   within  I-ZONE
C **********************************************************************
      WRITE (6,'(5X,2A,/)') '<'//ROUTINE(1:LEN_TRIM(ROUTINE))//'>:',
     &                    ' calculating 2D-lattice sums inside the slab'
      IF ( IPRINT.GE.2 ) WRITE (6,99001)
C
      MAD_II:DO IQ = NQ_L + 1,NQ_L + NQ_I
         DO JQ = NQ_L + 1,NQ_L + NQ_I
C
C make ewald sumation in plane and inverse space
C sum if rz<>0 (out of plane)
C
            CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,QBAS(1,IQ),QBAS(1,JQ),RM2,
     &                         NRMAX,NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,SMAT,
     &                         VOLUC_2D_AU,LMAX_RGNTMAD,NLMMAX_RGNTMAD)
C
C ----------------------------------------------------------------------
            IF ( IPRINT.GE.2 ) THEN
               WRITE (6,99002) IQ,JQ,SMAT(1)
               IF ( JQ.EQ.NQ_L+NQ_I .AND. IQ.NE.NQ_L+NQ_I )
     &               WRITE (6,'(20X,20(''-''))')
            END IF
C ----------------------------------------------------------------------
C
            CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                         IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                         NLMMAX_RGNTMAD,NRGNTMADMAX)
C
            AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD)
     &         = AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD) + AM(1:NLMMAD,1:NLMMAD)
C
C
         END DO
      END DO MAD_II
C
      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
C
      NREPBASALL_L = NREPBAS_L*NQ_L
      NREPBASALL_R = NREPBAS_R*NQ_R
C
C **********************************************************************
C                      IQ      within  I-ZONE
C                      JQ      within  L-BULK basis
C                      JQ_XL   within  extended L-BULK
C **********************************************************************
C
      WRITE (6,'(5X,2A,/)') '<'//ROUTINE(1:LEN_TRIM(ROUTINE))//'>:',
     &                   ' calculating 2D-lattice sums slab - left host'
      IF ( IPRINT.GE.2 ) WRITE (6,99001)
C
      MAD_IL:DO IQ = NQ_L + 1,NQ_L + NQ_I
C
         JQ_XL = 0
         DO IREPBAS_L = 1,NREPBAS_L
            DO JQ = 1,NQ_L
C
               VEC_JQ(1:3) = QBAS(1:3,JQ) - (IREPBAS_L-1)*ABAS_L(1:3,3)
C
               JQ_XL = JQ_XL + 1
C
C-->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
C     Inverse space sum if rz<>0 (out of plane)
C
               CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,QBAS(1,IQ),VEC_JQ,RM2,
     &                            NRMAX,NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,
     &                            SMAT,VOLUC_2D_AU,LMAX_RGNTMAD,
     &                            NLMMAX_RGNTMAD)
C
C ----------------------------------------------------------------------
               IF ( IPRINT.GE.2 ) THEN
                  WRITE (6,99002) IQ,JQ_XL,SMAT(1)
                  IF ( JQ_XL.EQ.NREPBASALL_L .AND. IQ.NE.NQ_I )
     &                 WRITE (6,'(20X,20(''-''))')
               END IF
C ----------------------------------------------------------------------
C
               CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                            IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                            NLMMAX_RGNTMAD,NRGNTMADMAX)
C
               AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD)
     &            = AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD)
     &            + AM(1:NLMMAD,1:NLMMAD)
C
            END DO         ! JQ loop in left host basis
         END DO            ! IREPBAS_L loop in layers to get convergence
C
      END DO MAD_IL
C **********************************************************************
C
      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
C
C **********************************************************************
C                      IQ      within  I-ZONE
C                      JQ      within  R-BULK basis
C                      JQ_XR   within  extended R-BULK
C **********************************************************************
      WRITE (6,'(5X,2A,/)') '<'//ROUTINE(1:LEN_TRIM(ROUTINE))//'>:',
     &                  ' calculating 2D-lattice sums slab - right host'
      IF ( IPRINT.GE.2 ) WRITE (6,99001)
C
      MAD_IR:DO IQ = NQ_L + 1,NQ_L + NQ_I
C
         JQ_XR = 0
         DO IREPBAS_R = 1,NREPBAS_R
            DO JQ = NQ - NQ_R + 1,NQ
C
               VEC_JQ(1:3) = QBAS(1:3,JQ) + (IREPBAS_R-1)*ABAS_R(1:3,3)
C
               JQ_XR = JQ_XR + 1
C
C-->  make ewald sumation (in plane) and
C     Inverse space sum if rz<>0 (out of plane)
C
               CALL SCFMAD2D_SMAT(NLMAD_LOC,ALAT,QBAS(1,IQ),VEC_JQ,RM2,
     &                            NRMAX,NSHLR,NSR,GN2,NGMAX,NSHLG,NSG,
     &                            SMAT,VOLUC_2D_AU,LMAX_RGNTMAD,
     &                            NLMMAX_RGNTMAD)
C
C ----------------------------------------------------------------------
               IF ( IPRINT.GE.2 ) THEN
                  WRITE (6,99002) IQ,JQ_XR,SMAT(1)
                  IF ( JQ_XR.EQ.NREPBASALL_R .AND. IQ.NE.NQ_I )
     &                 WRITE (6,'(20X,20(''-''))')
               END IF
C ----------------------------------------------------------------------
C
               CALL SCFMAD2D_COEF(.TRUE.,NLMAD_LOC,AM,BM,SMAT,RGNTMAD,
     &                            IRGNTMAD,NRGNTMAD,NLMMAD_LOC,
     &                            NLMMAX_RGNTMAD,NRGNTMADMAX)
C
               AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD)
     &            = AVMAD(IQ,JQ,1:NLMMAD,1:NLMMAD)
     &            + AM(1:NLMMAD,1:NLMMAD)
C
            END DO         ! JQ loop in right host basis
         END DO            ! IREPBAS_R loop in layers to get convergence
C
      END DO MAD_IR
C **********************************************************************
C
C ----------------------------------------------------------------------
      IF ( IPRINT.GE.2 ) WRITE (6,'(18X,22(''-''),/)')
      IF ( IPRINT.GT.0 ) CALL RMATSTRUCT('AVMAD Madelung matrix',AVMAD,
     &     NQ,NQ,1,0,0,1D-8,6)
C      END IF
C ######################################################################
C
99001 FORMAT (8X,'2D Lattice sum (LMXSP = 1)',/,18X,'  IQ  JQ  SUM',/,
     &        18X,23('-'))
99002 FORMAT (18X,2I5,D12.4)
99003 FORMAT (//,1X,79('*'),/,34X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'set up Madelung matrix for 2D-bulk calculations',//,10X,
     &        'NLMAD =',I4,5X,'NLMMAD =',I4,'   global setting',/,10X,
     &        'NLMAD =',I4,5X,'NLMMAD =',I4,
     &        '   local temporary setting',/)
      END
