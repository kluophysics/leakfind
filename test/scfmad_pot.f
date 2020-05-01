C*==scfmad_pot.f    processed by SPAG 6.70Rc at 07:51 on 26 Apr 2017
      SUBROUTINE SCFMAD_POT(USE_KLMFP,NLFP,NLMFP,VNEW,DROT_QLM)
C   ********************************************************************
C   *                                                                  *
C   * calculate the spin-independent inter-cell Madelung-potentials    *
C   * and add these to the type dependent spin-averaged potential V    *
C   * using the Madelung matrices   AVMAD                              *
C   * and site dependent charge-moments CMNTQ                          *
C   *                                                                  *
C   *   V(r,lm,IQ) =  (-r)**l * SUM_{JQ,l'm'}                          *
C   *                             A(IQ,JQ,lm,l'm')*CMNTQ(JQ,l'm')      *
C   *                                                                  *
C   * -----------------------------------------------------------------*
C   *  NOTE: a) The Madelung-potential at a specific lattice-site IQ   *
C   *           depends only on this lattice-site, but not on the      *
C   *           atom-type, it is occupied with.                        *
C   *        b) the Madelung-potential of equivalent sites my differ   *
C   *           but are connected by a symmetry operation              *
C   *        c) if a type  IT  occupies more than 1 site (NAT>1) its   *
C   *           potential may depend which site IQ is considered due   *
C   *           due to b)                                              *
C   *        d) only the representative sites IQ  (IQREPQ(IQ)=IQ) are  *
C   *           dealt with below. This way  V_Mad(IT) is the potential *
C   *           for the first site IQ it occurs                        *
C   *        e) CMNTQ includes the nuclear contribution -Z/SQRT(4*pi)  *
C   *                                                                  *
C   *                 based on   B. Drittler's  routines               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_CALCMODE,ONLY:BREAKPOINT,ROTATE_LATTICE
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_TYPES,ONLY:IMT,KLMFP,NTMAX,NLMFPMAX,ITBOT,ITTOP
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,NQHOST,CMNTQ,IQBOT,AVMAD,
     &    IQTOP,MAGROT_Q,KFP_LMQ,NQCLU,VLMMAD_HOST,VLMMAD_BACK,VMAD2D_A,
     &    VMAD2D_B,NQ_L,QBAS
      USE MOD_RMESH,ONLY:JRCUT,NPAN,R,NRMAX,FULLPOT
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_CONSTANTS,ONLY:PI
      USE MOD_FILES,ONLY:IOTMP,IPRINT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD_POT')
C
C Dummy arguments
C
      INTEGER NLFP,NLMFP
      LOGICAL USE_KLMFP
      REAL*8 DROT_QLM(NLMFP,NLMFP,IQBOT:IQTOP),
     &       VNEW(NRMAX,NLMFPMAX,NTMAX)
C
C Local variables
C
      REAL*8 AC,DLMLMP,RPWML(:),SQRT_1OV4PI,SQRT_3OV4PI,VIC(:)
      INTEGER IFLAG,IM,IO,IQ,IR,IRTOP,IT,IT1,JQ,JQBOT,JQTOP,L,LM,LM2,
     &        LMP,M
      LOGICAL KDONE_T(:),KFP_Q,KFP_T
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RPWML,VIC,KDONE_T
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (VIC(NRMAX),RPWML(NRMAX),KDONE_T(NTMAX))
C
      KDONE_T(1:NTMAX) = .FALSE.
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 3
      IF ( FULLPOT .AND. BREAKPOINT.EQ.3 ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'VLM_check_sum.dat')
         WRITE (IOTMP,99005) ROTATE_LATTICE
         DO IT = 1,ITBOT,ITTOP
            WRITE (IOTMP,99003) 'V_intra',IT
            IM = IMT(IT)
            IRTOP = JRCUT(NPAN(IM),IM)
            DO LM = 1,NLMFP
               WRITE (IOTMP,99004) LM,L_LM(LM),M_LM(LM),
     &                             SUM(VNEW(1:IRTOP,LM,IT))
            END DO
         END DO
         VNEW(:,:,:) = 0D0
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 3
C
C-----------------------------------------------------------------------
C                         host calculation
C-----------------------------------------------------------------------
      IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
         VLMMAD_HOST(1:NLMFPMAX,1:NQMAX) = 0D0
C
C----------------------------------------------------------- bulk system
         IF ( SYSTEM_DIMENSION(1:2).EQ.'3D' .OR. SUB_SYSTEM(2:6)
     &        .EQ.'-BULK' .OR. SYSTEM_DIMENSION(1:2).EQ.'0D' ) THEN
C
            JQBOT = IQBOT
            JQTOP = IQTOP
C
C------------------------------------------------------------- 2D system
         ELSE IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SUB_SYSTEM(1:6)
     &             .EQ.'I-ZONE' ) THEN
C
            JQBOT = 1
            JQTOP = NQ
C
            CALL SCFMAD2D_AB
C
         END IF
C
C-----------------------------------------------------------------------
C                          embedded cluster
C-----------------------------------------------------------------------
C
      ELSE
C
         JQBOT = NQHOST + 1
         JQTOP = NQHOST + NQCLU
C
      END IF
C-----------------------------------------------------------------------
C
      SQRT_1OV4PI = SQRT(1D0/(4D0*PI))
      SQRT_3OV4PI = SQRT(3D0/(4D0*PI))
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C--------------------------------------------------- loop over all sites
C--------- restrict to those representative for a symmetry-related class
C
      DO IQ = IQBOT,IQTOP
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            IT1 = ITOQ(1,IQ)
            IM = IMT(IT1)
            IRTOP = JRCUT(NPAN(IM),IM)
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO L = 0,NLFP - 1
C
               IF ( L.EQ.0 ) THEN
                  RPWML(1:NRMAX) = 1D0
               ELSE
                  DO IR = 1,IRTOP
                     RPWML(IR) = -RPWML(IR)*R(IR,IM)
                  END DO
               END IF
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
               DO M = -L,L
                  LM = L*L + L + M + 1
C
C=======================================================================
C
                  IF ( .NOT.USE_KLMFP .OR. KFP_LMQ(LM,IQ).NE.0 ) THEN
C
                     AC = 0.0D0
C
C-----------------------------------------------------------------------
C    NQ=1: the LM = 1 component disappears for pure and CPA case
C    CMNTQ(0,1) contains the nuclear charge contribution
C-----------------------------------------------------------------------
                     IF ( NQ.EQ.1 ) THEN
C
                        DO LM2 = 2,NLMFP
                           AC = AC + AVMAD(1,1,LM,LM2)*CMNTQ(LM2,1)
                        END DO
C
                     ELSE
C
                        DO JQ = JQBOT,JQTOP
                           DO LM2 = 1,NLMFP
                              AC = AC + AVMAD(IQ,JQ,LM,LM2)
     &                             *CMNTQ(LM2,JQ)
                           END DO
                        END DO
C
                     END IF
C
                     IF ( LM.EQ.1 .AND. IPRINT.GT.0 ) WRITE (6,99001)
     &                    IQ,(AC/SQRT(4.0D0*PI))
C
C-----------------------------------------------------------------------
C                         host calculation
C-----------------------------------------------------------------------
                     IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
                        IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &                       SUB_SYSTEM(1:6).EQ.'I-ZONE' ) THEN
C
                           IF ( LM.EQ.1 ) THEN
C
                              AC = AC + 
     &                             (VMAD2D_A*(QBAS(3,IQ)-QBAS(3,NQ_L))
     &                             +VMAD2D_B)/SQRT_1OV4PI
C
                           ELSE IF ( LM.EQ.3 ) THEN
C
C -------------------- inclusion of a L=(1,0) term seems to make problems
C                             AC = AC + VMAD2D_A/SQRT_3OV4PI
                              AC = AC*1D0
C
                           END IF
C
                        END IF
C
                        VLMMAD_HOST(LM,IQ) = AC
C
C-----------------------------------------------------------------------
C                          embedded cluster
C-----------------------------------------------------------------------
C
                     ELSE
C
                        AC = VLMMAD_BACK(LM,IQ) + AC
C
                     END IF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  add the type independent intercell-potential VIC (see NOTE above)
C-----------------------------------------------------------------------
C
                     VIC(1:IRTOP) = RPWML(1:IRTOP)*AC
C
C-----------------------------------------------------------------------
                     IF ( MAGROT_Q(IQ) ) THEN
C
                        DO LMP = 1,NLMFP
C
                           DLMLMP = DROT_QLM(LM,LMP,IQ)
C
                           KFP_T = KLMFP(LMP,IT1).EQ.1
                           KFP_Q = KFP_LMQ(LM,IQ).EQ.1 .AND. ABS(DLMLMP)
     &                             .GT.1D-6
C
C                           IF ( (KFP_T .AND. .NOT.KFP_Q) .OR.
C     &                          (.NOT.KFP_T .AND. KFP_Q) )
C     &                          WRITE (6,99002) IQ,LM,KFP_Q,
C     &                          KFP_LMQ(LM,IQ),IT1,LMP,KFP_T,
C     &                          KLMFP(LMP,IT1)
C99002 FORMAT (/,1X,79('W'),/,10X,'TROUBLE in <SCFMAD_POT>',/,10X,
C     &        '  IQ   LM ',2I6,'KFP (LM)',2I6,/,10X,'  IT1  LMP',2I6,
C     &        'KFP (LM)',2I6,/,1X,79('W'),/)
C
                           IF ( KFP_T ) THEN
C
                              DO IO = 1,NOQ(IQ)
                                 IT = ITOQ(IO,IQ)
                                 KDONE_T(IT) = .TRUE.
C
                                 DO IR = 1,IRTOP
                                    VNEW(IR,LMP,IT) = VNEW(IR,LMP,IT)
     &                                 + VIC(IR)*DLMLMP
                                 END DO
C
                              END DO
C
                           END IF
C
                        END DO
C
                     ELSE
C-----------------------------------------------------------------------
C
                        DO IO = 1,NOQ(IQ)
                           IT = ITOQ(IO,IQ)
                           KDONE_T(IT) = .TRUE.
C
                           DO IR = 1,IRTOP
                              VNEW(IR,LM,IT) = VNEW(IR,LM,IT) + VIC(IR)
                           END DO
C
                        END DO
C
                     END IF
C-----------------------------------------------------------------------
C
                  END IF
C=======================================================================
C
               END DO
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
            END DO
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         END IF
      END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      IFLAG = 0
      DO IT = ITBOT,ITTOP
         IF ( .NOT.KDONE_T(IT) ) THEN
            WRITE (6,99002) IT,ITBOT,ITTOP
            IFLAG = 1
         END IF
      END DO
      IF ( IFLAG.EQ.1 ) STOP '<SCFMAD_POT> not all types treated'
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 3
      IF ( FULLPOT .AND. BREAKPOINT.EQ.3 ) THEN
         DO IT = 1,ITBOT,ITTOP
            WRITE (IOTMP,99003) 'V_madel',IT
            IM = IMT(IT1)
            IRTOP = JRCUT(NPAN(IM),IM)
            DO LM = 1,NLMFP
               WRITE (IOTMP,99004) LM,L_LM(LM),M_LM(LM),
     &                             SUM(VNEW(1:IRTOP,LM,IT))
            END DO
         END DO
C
         CALL STOP_BREAKPOINT(ROUTINE)
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 3
C
C=======================================================================
99001 FORMAT (5x,'spherically averaged Madelung-potential for site',
     &        ' IQ=',I2,': ',1P,D14.6)
99002 FORMAT (10X,'within range ',I3,' - ',I3,' type IT =',I3,
     &        ' not treated')
99003 FORMAT (/,10X,'V(l,m) check sum for ',A,' IT=',I4,/)
99004 FORMAT (3I3,2F20.10)
99005 FORMAT (10X,'BREAKPOINT in <SCFMAD_POT>',//,10X,
     &        'ROTATE_LATTICE = ',L3,/)
      END
