C*==sfncreate.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SFNCREATE(NT,IQAT,ITOQ,NAT,RMTRED0,NSFLIM,NPANLIM,
     &                     NRSFLIM)
C   ********************************************************************
C   *                                                                  *
C   *  create the shape functions for full potential calculations      *
C   *  - write results to file DATASET0.sfn                            *
C   *  - create data and script files for  rasmol                      *
C   *  - return limits NSFLIM, NPANLIM, NRSFLIM for storage allocation *
C   *                                                                  *
C   *  NTMAX, IQAT, ITOQ, NAT are already fixed                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:SQRT_4PI,PI
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_TYPES,ONLY:NTMAX,NLMFPMAX
      USE MOD_ANGMOM,ONLY:NL,L_LM,M_LM
      USE MOD_RMESH,ONLY:NM,NMMAX,NRMAX,NLMSFMAX,JRCUT_SAV,NPAN_SAV
      USE MOD_FILES,ONLY:DATSET0,LDATSET0,IPRINT,SFNFIL,LSFNFIL,IOTMP,
     &    IDUMMY,RDUMMY
      USE MOD_SITES,ONLY:QBAS,NQ,NQMAX,NQ_L,NQ_R,IMQ
      USE MOD_LATTICE,ONLY:ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,ADAINV_I,
     &    ADAINV_R,SYSTEM_DIMENSION,VOLUC,ALAT
      USE MOD_SYMMETRY,ONLY:NSFTSYMQ,SYMTVEC,MROTR,ISYMGENQ,IQREPQ
      USE MOD_CALCMODE,ONLY:FP_USE_DIRECTIONS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNCREATE')
      INTEGER NVERTMAX,NFACEMAX,NPANMAX_S,NRSFMAX_S
      PARAMETER (NVERTMAX=300,NFACEMAX=1600,NPANMAX_S=600,
     &           NRSFMAX_S=2000)
      LOGICAL CHECK_SYM,CHECK_SCL
      PARAMETER (CHECK_SYM=.FALSE.,CHECK_SCL=.FALSE.)
C
C Dummy arguments
C
      INTEGER NPANLIM,NRSFLIM,NSFLIM,NT
      INTEGER IQAT(NQMAX,NTMAX),ITOQ(NTMAX,NQMAX),NAT(NTMAX)
      REAL*8 RMTRED0(NMMAX)
C
C Local variables
C
      REAL*8 A3(:),AUX,B3(:),C3(:),CLURAD_SFN,D3(:),DEDW,DLT,DQBAS(:,:),
     &       DRN_S(:,:),DXH,ERR1,ERRA,ERRB,ERRI,FLMSF_PR(:,:),
     &       FLMSF_S(:,:,:),RINT0(:),RMTFILL_M(:),RMTREDFILL,ROUT,
     &       RQCLU_SFN(:,:),RVERT_VFM(:,:,:,:),SCALE_S,SIZEFAC(:),
     &       SUMVOL_NUM,VOLREDIS,VOLREDIS_NUM,VOLUC_AU,VOLUME,VOL_M(:),
     &       VOL_NUM_M(:),WEIGHT0,WEIGHTCLU_SFN(:),WH,WH2,WH2A,WH2B,
     &       W_RADINT_S(:,:),XEDGE(:,:),XH,XH1,XH2,XRN_S(:,:),YEDGE(:,:)
     &       ,ZEDGE(:,:)
      LOGICAL ACCEPT,KDONE(:),KEEP_SFN_S,MOL
      INTEGER I,IA_ERR,IED,IFC,IFLAG,IFLAG_POLYHEDRON,IFLAG_VERTEX,IM,
     &        IMCURR,IPAN,IPRINT_LOW,IPROC,IPROCM(:),IQ,IQCLU_SFN,
     &        IQCNTR,IQ_MAUX(:),IQ_QCLU_SFN(:),IRCUT(:),IRSF,ISF,
     &        ISFLM_PR(:,:),ISF_S,IT,ITER,IWSF,J,KEY_PAN,KFP_LMQ_PR(:,:)
     &        ,KLMFP_PR(:,:),KLMSF_PR(:,:),LM,LMIFP_PR(:,:),
     &        LMISF_PR(:,:),LMISF_S(:,:),LSF,LSFMAX,LTXT_PR(:),
     &        MESHN_S(:),NCELL_S,NEDGE(:),NEDGE_FCM(:,:),NFACE,
     &        NFACE_M(:),NFPT_PR(:),NFUN_S(:),NLMFPT_PR(:),NMSF(:),
     &        NM_S(:,:),NPAN_S(:),NQCLU_I_SFN,NQCLU_L_SFN,NQCLU_R_SFN,
     &        NQCLU_SFN,NRSF0,NSF_PR(:),NSF_PRMAX,NSHLCLU_SFN,NWARN
      CHARACTER*8 TXT_PR(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA DLT/0.05D0/,NRSF0/200/
C
      ALLOCATABLE SIZEFAC,VOL_M,DQBAS,NPAN_S,NFUN_S,IPROCM,KDONE
      ALLOCATABLE RQCLU_SFN,IQ_QCLU_SFN,WEIGHTCLU_SFN,VOL_NUM_M
      ALLOCATABLE A3,B3,C3,D3,XEDGE,YEDGE,ZEDGE,NEDGE
      ALLOCATABLE RVERT_VFM,NEDGE_FCM,W_RADINT_S
      ALLOCATABLE NMSF,NM_S,MESHN_S,LMISF_S,XRN_S,RMTFILL_M,IRCUT
      ALLOCATABLE DRN_S,FLMSF_S,NFACE_M,FLMSF_PR,KFP_LMQ_PR,RINT0
      ALLOCATABLE KLMFP_PR,NLMFPT_PR,NFPT_PR,LMIFP_PR,ISFLM_PR
      ALLOCATABLE KLMSF_PR,NSF_PR,LMISF_PR,TXT_PR,LTXT_PR,IQ_MAUX
C
      LSF = 4*(NL-1)
      LSFMAX = LSF
C
      ALLOCATE (SIZEFAC(NTMAX),NFACE_M(NMMAX),VOL_M(NMMAX))
      ALLOCATE (A3(NFACEMAX),B3(NFACEMAX),IPROCM(NMMAX))
      ALLOCATE (C3(NFACEMAX),D3(NFACEMAX),KDONE(NMMAX))
      ALLOCATE (RVERT_VFM(3,NVERTMAX,NFACEMAX,NMMAX))
      ALLOCATE (XEDGE(NVERTMAX,NFACEMAX))
      ALLOCATE (YEDGE(NVERTMAX,NFACEMAX),VOL_NUM_M(NMMAX))
      ALLOCATE (ZEDGE(NVERTMAX,NFACEMAX),IRCUT(0:NPANMAX_S))
      ALLOCATE (NEDGE_FCM(NFACEMAX,NMMAX),NMSF(NPANMAX_S))
      ALLOCATE (NEDGE(NFACEMAX),NPAN_S(NMMAX),NFUN_S(NMMAX))
      ALLOCATE (NM_S(NPANMAX_S,NMMAX),MESHN_S(NMMAX))
      ALLOCATE (LMISF_S(NLMSFMAX,NMMAX),RMTFILL_M(NMMAX))
      ALLOCATE (DRN_S(NRSFMAX_S,NMMAX),XRN_S(NRSFMAX_S,NMMAX))
      ALLOCATE (RINT0(NRMAX))
      ALLOCATE (W_RADINT_S(NRSFMAX_S,NMMAX))
      ALLOCATE (FLMSF_S(NRSFMAX_S,NLMSFMAX,NMMAX))
      ALLOCATE (FLMSF_PR(NRSFMAX_S,NLMSFMAX),DQBAS(3,NQMAX))
C
      DQBAS(1:3,1:NQMAX) = 0D0
      NPAN_S(1:NMMAX) = 0
      NEDGE_FCM(1:NFACEMAX,1:NMMAX) = 0
      LMISF_S(1:NLMSFMAX,1:NMMAX) = 0
      NFUN_S(1:NMMAX) = 0
      LMISF_S(1:NLMSFMAX,1:NMMAX) = 0
      XRN_S(1:NRSFMAX_S,1:NMMAX) = 0D0
      DRN_S(1:NRSFMAX_S,1:NMMAX) = 0D0
      FLMSF_S(1:NRSFMAX_S,1:NLMSFMAX,1:NMMAX) = 0D0
      NM_S(1:NPANMAX_S,1:NMMAX) = 0
C
      MOL = .FALSE.
      IF ( IPRINT.LE.0 ) THEN
         IPRINT_LOW = -1
      ELSE
         IPRINT_LOW = IPRINT
      END IF
C
C-------------------------------------------------------- dummy settings
      DO IT = 1,NTMAX
         SIZEFAC(IT) = 1D0
      END DO
C
C=======================================================================
C             find the picking rules for the shape functions
C=======================================================================
      NSF_PRMAX = NLMSFMAX
C
      ALLOCATE (KFP_LMQ_PR(NLMFPMAX,NQMAX))
      ALLOCATE (KLMFP_PR(NLMFPMAX,NTMAX))
      ALLOCATE (ISFLM_PR(NLMSFMAX,NMMAX))
      ALLOCATE (NLMFPT_PR(NTMAX),NFPT_PR(NTMAX))
      ALLOCATE (LMIFP_PR(NLMFPMAX,NTMAX))
      ALLOCATE (KLMSF_PR(NLMSFMAX,NMMAX))
      ALLOCATE (LMISF_PR(NSF_PRMAX,NMMAX))
      ALLOCATE (NSF_PR(NMMAX),IQ_MAUX(NMMAX))
      ALLOCATE (TXT_PR(NTMAX),LTXT_PR(NTMAX))
C
      NSF_PR(1:NMMAX) = 999999
      NFPT_PR(1:NTMAX) = 999999
C
      NLMFPT_PR(:) = 0
      KLMFP_PR(:,:) = 0
      LMIFP_PR(:,:) = 0
      KLMSF_PR(:,:) = 0
      LMISF_PR(:,:) = 0
C
      DO IT = 1,NTMAX
         TXT_PR(IT) = '  '
         LTXT_PR(IT) = 2
      END DO
      DO IM = 1,NM
         IQ_MAUX(IM) = 0
         DO IT = 1,NT
            IQ = IQAT(1,IT)
            IF ( IMQ(IQ).EQ.IM ) THEN
               IQ_MAUX(IM) = IQ
               EXIT
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C              determine non-0 shape functions via symmetry
C-----------------------------------------------------------------------
C
      IF ( CHECK_SYM ) CALL FPPICKRULES('SHAPE',NMMAX,TXT_PR,LTXT_PR,
     &                                  IMQ,IQAT,NAT,KFP_LMQ_PR,
     &                                  KLMFP_PR,NLMFPT_PR,NFPT_PR,
     &                                  LMIFP_PR,ISFLM_PR,KLMSF_PR,
     &                                  NSF_PR,LMISF_PR,NLMFPMAX,
     &                                  NLMSFMAX)
C
C=======================================================================
C
      WRITE (6,99012) NM
C
      KEY_PAN = 0
C
      IFLAG_POLYHEDRON = 0
      IFLAG_VERTEX = 0
      IFLAG = 0
      NWARN = 0
C
      CALL MPI_DISTRIBUTE(IPROCM,NM,MPI,'M')
C
      NSHLCLU_SFN = 25
C            NSHLCLU_SFN = 32
C=======================================================================
      DO IM = 1,NM
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCM(IM) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            WRITE (6,99017) IM
C
 20         CONTINUE
            NQCLU_SFN = 0
            NQCLU_L_SFN = 0
            NQCLU_I_SFN = 0
            NQCLU_R_SFN = 0
C
            IQ = IQ_MAUX(IM)
C
            IQCNTR = IQ
            CLURAD_SFN = 0D0
C
            IF ( IPRINT.GT.0 ) WRITE (6,99018) IM,'CLUSSITES'
C
            CALL CLUSSITES(IOTMP,IPRINT_LOW,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,CLURAD_SFN,IQCNTR,NQCLU_SFN,
     &                     NQCLU_L_SFN,NQCLU_I_SFN,NQCLU_R_SFN,
     &                     NSHLCLU_SFN,NQ,NQ_L,NQ_R,NQMAX)
C
            IF ( ALLOCATED(RQCLU_SFN) )
     &           DEALLOCATE (RQCLU_SFN,IQ_QCLU_SFN,WEIGHTCLU_SFN)
            ALLOCATE (RQCLU_SFN(3,NQCLU_SFN),WEIGHTCLU_SFN(NQCLU_SFN))
            ALLOCATE (IQ_QCLU_SFN(NQCLU_SFN),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQ_QCLU_SFN')
            WEIGHTCLU_SFN(1:NQCLU_SFN) = 0D0
C
            READ (IOTMP) ((RQCLU_SFN(J,I),J=1,3),RDUMMY,IQ_QCLU_SFN(I),
     &                   I=1,NQCLU_SFN),(IDUMMY,I=1,NSHLCLU_SFN)
            CLOSE (IOTMP)
C
            DO IQCLU_SFN = 2,NQCLU_SFN
               DO J = 1,3
                  RQCLU_SFN(J,IQCLU_SFN-1) = RQCLU_SFN(J,IQCLU_SFN)
               END DO
               IQ_QCLU_SFN(IQCLU_SFN-1) = IQ_QCLU_SFN(IQCLU_SFN)
            END DO
            NQCLU_SFN = NQCLU_SFN - 1
C
            IT = ITOQ(1,IQ)
            WEIGHT0 = SIZEFAC(IT)
C
            DO IQCLU_SFN = 1,NQCLU_SFN
               IT = ITOQ(1,IQ_QCLU_SFN(IQCLU_SFN))
               WEIGHTCLU_SFN(IQCLU_SFN) = SIZEFAC(IT)
            END DO
C
            A3(1:NFACEMAX) = 0D0
            B3(1:NFACEMAX) = 0D0
            C3(1:NFACEMAX) = 0D0
            D3(1:NFACEMAX) = 0D0
            XEDGE(1:NVERTMAX,1:NFACEMAX) = 0D0
            YEDGE(1:NVERTMAX,1:NFACEMAX) = 0D0
            ZEDGE(1:NVERTMAX,1:NFACEMAX) = 0D0
C
            IF ( IPRINT.GT.0 ) WRITE (6,99018) IM,'SFNVORONOI'
C
            CALL SFNVORONOI(NQCLU_SFN,RQCLU_SFN,NVERTMAX,WEIGHT0,
     &                      WEIGHTCLU_SFN,RMTREDFILL,ROUT,VOLUME,NFACE,
     &                      A3,B3,C3,D3,NEDGE,XEDGE,YEDGE,ZEDGE,
     &                      NFACEMAX,IFLAG_VERTEX)
C
C----------------------------------- make 2nd attempt in case of trouble
            IF ( IFLAG_VERTEX.NE.0 ) THEN
               IFLAG_VERTEX = -1
               NSHLCLU_SFN = 32
               WRITE (6,99016) NSHLCLU_SFN
               GOTO 20
            END IF
C
            VOL_M(IM) = VOLUME*ALAT**3
C
            NFACE_M(IM) = NFACE
            DO IFC = 1,NFACE
               NEDGE_FCM(IFC,IM) = NEDGE(IFC)
               DO IED = 1,NEDGE(IFC)
                  RVERT_VFM(1,IED,IFC,IM) = XEDGE(IED,IFC)
                  RVERT_VFM(2,IED,IFC,IM) = YEDGE(IED,IFC)
                  RVERT_VFM(3,IED,IFC,IM) = ZEDGE(IED,IFC)
               END DO
            END DO
C
            IMCURR = IM
C
            IF ( IPRINT.GT.0 ) THEN
               WRITE (6,99018) IM,'SFNRASMOL'
C
               CALL SFNRASMOL(1,IMCURR,NQ,QBAS,DQBAS,DATSET0,LDATSET0,
     &                        IMQ,NQCLU_SFN,RQCLU_SFN,IOTMP,NFACE_M,
     &                        NEDGE_FCM,RVERT_VFM,NSFTSYMQ,SYMTVEC,
     &                        MROTR,ISYMGENQ,IQREPQ,NVERTMAX,NFACEMAX,
     &                        NMMAX,NQMAX)
            END IF
C
C---------------------- Calculate shape functions for each Voronoi shape
C
            NMSF(1:NPANMAX_S) = 0
C
            IF ( IPRINT.GT.0 ) WRITE (6,99018) IM,'SFNSHAPE'
C
            CALL SFNSHAPE(NRSF0,A3,B3,C3,D3,NEDGE,XEDGE,YEDGE,ZEDGE,
     &                    NFACE,LSF,DLT,KEY_PAN,NMSF,NCELL_S,SCALE_S,
     &                    NPAN_S(IM),MESHN_S(IM),NM_S(1,IM),XRN_S(1,IM),
     &                    DRN_S(1,IM),NFUN_S(IM),LMISF_S(1,IM),
     &                    FLMSF_S(1,1,IM),NPANMAX_S,NRSFMAX_S,LSFMAX,
     &                    NLMSFMAX,NFACEMAX,NVERTMAX,IFLAG_POLYHEDRON)
C
            IF ( IFLAG_POLYHEDRON.EQ.1 ) THEN
C
               IFLAG_POLYHEDRON = 0
               IF ( NSHLCLU_SFN.LT.50 ) THEN
                  NSHLCLU_SFN = NSHLCLU_SFN + 4
                  WRITE (6,99001) NSHLCLU_SFN
                  GOTO 20
               ELSE
                  CALL STOP_MESSAGE(ROUTINE,'NSHLCLU_SFN got to large')
               END IF
            END IF
C
C--CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM
            IF ( CHECK_SYM ) THEN
C
C--------------------------------------------------- check picking rules
C
               IF ( NSF_PR(IM).NE.NFUN_S(IM) ) THEN
                  WRITE (6,99008) 'number of non-0 terms inconsistent'
                  WRITE (6,99010) 'IM               ',IM
                  WRITE (6,99010) 'NSF picking rules',NSF_PR(IM)
                  WRITE (6,99010) 'NSF in <SFNSHAPE>',NFUN_S(IM)
                  IFLAG = MAX(1,IFLAG)
                  NWARN = NWARN + 1
                  IF ( NFUN_S(IM).GT.NSF_PR(IM) ) IFLAG = MAX(2,IFLAG)
               END IF
C
               ACCEPT = NSF_PR(IM).GE.NFUN_S(IM)
               DO ISF = 1,NFUN_S(IM)
                  LM = LMISF_S(ISF,IM)
                  IF ( KLMSF_PR(LM,IM).NE.1 ) THEN
                     WRITE (6,99009) 
     &                             'set of non-0 terms inconsistent for'
                     WRITE (6,99010) 'IM       LM      ',IM,LM
                     IFLAG = MAX(2,IFLAG)
                     ACCEPT = .FALSE.
                  END IF
               END DO
               IF ( IFLAG.GE.1 ) THEN
                  WRITE (6,99010) 'KLMSF  from picking rules'
                  WRITE (6,'(10X,20I5)')
     &                   (LMISF_PR(ISF,IM),ISF=1,NSF_PR(IM))
                  WRITE (6,99010) 
     &                          'KLMSF  from shape function calculation'
                  WRITE (6,'(10X,20I5)')
     &                   (LMISF_S(ISF,IM),ISF=1,NFUN_S(IM))
C
C-----------------------------------------------------------------------
C   if the picking rules allow for more non-zero shape functions:
C   fill up table with "missing terms" -- this case may happen if
C   not all restricted symmetry operations are used or recognized
C
C   KEEP_SFN_S = .TRUE. -->  DON'T insert dummy SFNs
C
                  KEEP_SFN_S = .TRUE.
C
                  IF ( ACCEPT .AND. NSF_PR(IM).GT.NFUN_S(IM) .AND. 
     &                 .NOT.KEEP_SFN_S ) THEN
C
                     FLMSF_PR(1:NRSFMAX_S,1:NSF_PR(IM)) = 0D0
C
                     DO ISF = 1,NSF_PR(IM)
                        LM = LMISF_PR(ISF,IM)
                        DO ISF_S = 1,NFUN_S(IM)
                           IF ( LMISF_S(ISF_S,IM).EQ.LM )
     &                          FLMSF_PR(1:NRSFMAX_S,ISF)
     &                          = FLMSF_S(1:NRSFMAX_S,ISF_S,IM)
                        END DO
                     END DO
C
                     NFUN_S(IM) = NSF_PR(IM)
                     FLMSF_S(1:NRSFMAX_S,1:NFUN_S(IM),IM)
     &                  = FLMSF_PR(1:NRSFMAX_S,1:NFUN_S(IM))
                     LMISF_S(1:NFUN_S(IM),IM)
     &                  = LMISF_PR(1:NFUN_S(IM),IM)
                     WRITE (6,99013) IM
                  END IF
C
               END IF
C
            END IF
C--CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM
C
C
C--CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL
C
            VOLREDIS = VOL_M(IM) - (4D0*PI/3D0)*(RMTREDFILL*ALAT)**3
C
            IRCUT(0) = 0
            DO IPAN = 1,NPAN_S(IM)
               IRCUT(IPAN) = IRCUT(IPAN-1) + NM_S(IPAN,IM)
            END DO
C
            DO IRSF = 1,MESHN_S(IM)
               RINT0(IRSF) = XRN_S(IRSF,IM)**2*FLMSF_S(IRSF,1,IM)
            END DO
C
            IF ( IPRINT.GT.0 ) WRITE (6,99018) IM,
     &                                'GET_INTEGRATION_WEIGHTS'
C
            CALL GET_INTEGRATION_WEIGHTS(NPAN_S(IM),IRCUT(0),MESHN_S(IM)
     &         ,W_RADINT_S(1,IM))
C
            XH1 = XRN_S(1,IM)
            XH2 = XRN_S(MESHN_S(IM),IM)
            DXH = XH2 - XH1
            WH2 = 1D0
C
            ITER = 0
C
            WH2A = 0D0
            ERRA = 0D0
 40         CONTINUE
            ITER = ITER + 1
C
            AUX = 0.0D0
            DO IRSF = 1,MESHN_S(IM)
               XH = XRN_S(IRSF,IM)
               WH = ((XH2-XH)+WH2*(XH-XH1))/DXH
               AUX = AUX + RINT0(IRSF)*WH*DRN_S(IRSF,IM)
     &               *W_RADINT_S(IRSF,IM)
            END DO
C
            VOLREDIS_NUM = AUX*SQRT_4PI*ALAT**3
C
            ERRI = (1D0-VOLREDIS_NUM/VOLREDIS)
C
            IF ( ITER.EQ.1 ) ERR1 = ERRI
C
            IF ( CHECK_SCL ) WRITE (6,*) 'ITERATION ',ITER,WH2,ERRI
C
            IF ( ABS(ERRI).GT.1D-10 ) THEN
C
               ERRB = ERRA
               WH2B = WH2A
               ERRA = ERRI
               WH2A = WH2
               IF ( ITER.EQ.1 ) THEN
                  WH2 = 1D0 + SIGN(0.01D0,ERRA)
               ELSE
                  DEDW = (ERRA-ERRB)/(WH2A-WH2B)
                  WH2 = WH2 - ERRA/DEDW
               END IF
C
               IF ( ITER.LE.20 ) GOTO 40
C
            END IF
            IF ( CHECK_SCL ) WRITE (6,*) 'ITERATION ',ITER,WH2
C
C--------------------------------------------- scale all shape functions
            DO IRSF = 1,MESHN_S(IM)
               XH = XRN_S(IRSF,IM)
               WH = ((XH2-XH)+WH2*(XH-XH1))/DXH
               DO ISF = 1,NFUN_S(IM)
                  FLMSF_S(IRSF,ISF,IM) = FLMSF_S(IRSF,ISF,IM)*WH
               END DO
            END DO
C
            WRITE (6,99006) ERR1,WH2
C
            IF ( ABS(ERRI).GT.1D-10 ) THEN
               WRITE (6,99005) IM,VOLREDIS,VOLREDIS_NUM,
     &                         VOLREDIS/VOLREDIS_NUM
               CALL STOP_MESSAGE(ROUTINE,'volume scaling not converged')
            END IF
C
C--CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL-CHECK_VOL
C
            RMTFILL_M(IM) = RMTREDFILL*ALAT
            IF ( RMTRED0(IM).LT.1D-6 ) RMTRED0(IM) = RMTREDFILL
C
            IF ( IPRINT.GT.0 ) WRITE (6,99018) IM,'SFNMTMESH'
C
            CALL SFNMTMESH(IM,NPAN_S(IM),MESHN_S(IM),NM_S(1,IM),
     &                     XRN_S(1,IM),DRN_S(1,IM),NFUN_S(IM),
     &                     FLMSF_S(1,1,IM),RMTREDFILL,RMTRED0(IM),
     &                     NRSFMAX_S,NPANMAX_S,NLMSFMAX)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C=======================================================================
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C collect: NPAN_S  MESHN_S  NM_S    XRN_S     DRN_S   NFUN_S
C collect: FLMSF_S LMISF_S  NFACE_M NEDGE_FCM RVERT_VFM  RMT RMTFILL_M
C
      IF ( MPI ) THEN
C=======================================================================
         CALL DRV_MPI_BARRIER
C
         DO IM = 1,NM
            IPROC = IPROCM(IM)
C
            CALL DRV_MPI_SEND_I(NPAN_S(IM),1,IPROC,1)
            CALL DRV_MPI_SEND_I(MESHN_S(IM),1,IPROC,2)
            CALL DRV_MPI_SEND_I(NM_S(1,IM),NPANMAX_S,IPROC,3)
            CALL DRV_MPI_SEND_I(NFUN_S(IM),1,IPROC,4)
            CALL DRV_MPI_SEND_I(LMISF_S(1,IM),NLMSFMAX,IPROC,5)
            CALL DRV_MPI_SEND_I(NFACE_M(IM),1,IPROC,6)
            CALL DRV_MPI_SEND_I(NEDGE_FCM(1,IM),NFACEMAX,IPROC,7)
C
            CALL DRV_MPI_SEND_R(XRN_S(1,IM),NRSFMAX_S,IPROC,8)
            CALL DRV_MPI_SEND_R(DRN_S(1,IM),NRSFMAX_S,IPROC,9)
            CALL DRV_MPI_SEND_R(FLMSF_S(1,1,IM),(NRSFMAX_S*NLMSFMAX),
     &                          IPROC,10)
            CALL DRV_MPI_SEND_R(RVERT_VFM(1,1,1,IM),
     &                          (3*NVERTMAX*NFACEMAX),IPROC,11)
            CALL DRV_MPI_SEND_R(RMTFILL_M(IM),1,IPROC,12)
            CALL DRV_MPI_SEND_R(VOL_M(IM),1,IPROC,13)
            CALL DRV_MPI_SEND_R(W_RADINT_S(1,IM),NRSFMAX_S,IPROC,14)
C
         END DO
C=======================================================================
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_ID.EQ.0 ) THEN
C
C=======================================================================
         IWSF = 20
C
         OPEN (IWSF,FILE=SFNFIL(1:LSFNFIL))
C
         WRITE (IWSF,FMT='(I5)') NM
         DO IM = 1,NM
            WRITE (IWSF,FMT='(D20.12)') 1D0
         END DO
C
         NPANLIM = 0
         NRSFLIM = 0
         NSFLIM = 0
C
         DO IM = 1,NM
C
            WRITE (6,99014) IM,'NSF ',NFUN_S(IM),'NRSFTOT ',MESHN_S(IM)
            WRITE (6,99015) (LMISF_S(ISF,IM),L_LM(LMISF_S(ISF,IM)),
     &                      M_LM(LMISF_S(ISF,IM)),ISF=1,NFUN_S(IM))
C
            NPANLIM = MAX(NPANLIM,NPAN_S(IM)+1)
            NRSFLIM = MAX(NRSFLIM,MESHN_S(IM))
            NSFLIM = MAX(NSFLIM,NFUN_S(IM))
C
            CALL SFNWRITE(IWSF,NPAN_S(IM),MESHN_S(IM),NM_S(1,IM),
     &                    XRN_S(1,IM),DRN_S(1,IM),NFUN_S(IM),
     &                    FLMSF_S(1,1,IM),LMISF_S(1,IM),NRSFMAX_S,
     &                    NPANMAX_S,NLMSFMAX)
C
         END DO
C
         CLOSE (IWSF)
C
C=======================================================================
C    keep NPAN and JRCUT in case it is overwritten during ASA start
C=======================================================================
         ALLOCATE (JRCUT_SAV(0:NPANLIM,NMMAX),NPAN_SAV(NMMAX))
         NPAN_SAV(:) = 0
         JRCUT_SAV(:,:) = 0
C
         DO IM = 1,NM
            NPAN_SAV(IM) = NPAN_S(IM)
            JRCUT_SAV(0,IM) = 0
            DO IPAN = 1,NPAN_S(IM)
               JRCUT_SAV(IPAN,IM) = JRCUT_SAV(IPAN-1,IM) + NM_S(IPAN,IM)
            END DO
         END DO
C=======================================================================
C
         WRITE (6,99003)
         VOLUME = 0D0
         SUMVOL_NUM = 0D0
         KDONE(1:NMMAX) = .FALSE.
         DO IQ = 1,NQ
            IM = IMQ(IQ)
C
            IF ( .NOT.KDONE(IM) ) THEN
C
               IRCUT(0) = 0
               DO IPAN = 1,NPAN_S(IM)
                  IRCUT(IPAN) = IRCUT(IPAN-1) + NM_S(IPAN,IM)
               END DO
C
               CALL GET_INTEGRATION_WEIGHTS(NPAN_S(IM),IRCUT(0),
     &            MESHN_S(IM),W_RADINT_S(1,IM))
C
               AUX = 0.0D0
               DO IRSF = 1,MESHN_S(IM)
                  AUX = AUX + XRN_S(IRSF,IM)**2*FLMSF_S(IRSF,1,IM)
     &                  *W_RADINT_S(IRSF,IM)*DRN_S(IRSF,IM)
               END DO
C
               VOL_NUM_M(IM) = (4D0*PI/3D0)*(RMTRED0(IM)*ALAT)
     &                         **3 + AUX*ALAT**3*SQRT_4PI
C
               WRITE (6,99004) IM,VOL_M(IM),VOL_NUM_M(IM),RMTRED0(IM)
     &                         *ALAT,RMTFILL_M(IM)
C
               KDONE(IM) = .TRUE.
C
            END IF
C
            VOLUME = VOLUME + VOL_M(IM)
            SUMVOL_NUM = SUMVOL_NUM + VOL_NUM_M(IM)
         END DO
         VOLUC_AU = VOLUC*ALAT**3
         WRITE (6,99007) VOLUME,SUMVOL_NUM,VOLUC_AU
C
C
C=======================================================================
C
         IF ( IPRINT.GT.0 ) CALL SFNRASMOL(2,IMCURR,NQ,QBAS,DQBAS,
     &        DATSET0,LDATSET0,IMQ,NQCLU_SFN,RQCLU_SFN,IOTMP,NFACE_M,
     &        NEDGE_FCM,RVERT_VFM,NSFTSYMQ,SYMTVEC,MROTR,ISYMGENQ,
     &        IQREPQ,NVERTMAX,NFACEMAX,NMMAX,NQMAX)
C
C=======================================================================
C          create Lebedev grid on the surface of the polyhedra
C=======================================================================
C
         IF ( FP_USE_DIRECTIONS )
     &        CALL SFNLEBEDEV(NFACE_M,NEDGE_FCM,RVERT_VFM,NVERTMAX,
     &        NFACEMAX)
C
C=======================================================================
C  create grid on the surface of the Voronoi polyhedra for integration
C=======================================================================
C
         CALL SFN_SURFINT_INIT(NFACE_M,NEDGE_FCM,RVERT_VFM,NVERTMAX,
     &                         NFACEMAX)
C
C=======================================================================
C
         WRITE (6,99011)
         WRITE (6,99010) 'running <SCFNCREATE> '
         WRITE (6,99010) 'mumber of warnings   NWARN =',NWARN
         WRITE (6,99010) 'error flag           IFLAG =',IFLAG
         WRITE (6,99010) 'polyhedron flag      IFLAG =',IFLAG_POLYHEDRON
         WRITE (6,99011)
C
         WRITE (6,99011)
         IF ( IFLAG.EQ.2 .OR. IFLAG_POLYHEDRON.NE.0 )
     &        CALL STOP_MESSAGE(ROUTINE,'INCONSISTENCIES')
C
         IF ( ABS(VOLUME-SUMVOL_NUM).GT.1D-5 ) THEN
            WRITE (6,99002)
            CALL STOP_MESSAGE(ROUTINE,'INCONSISTENCIES')
         END IF
         IF ( ABS(VOLUME-VOLUC_AU).GT.1D-6*NQ )
     &        CALL STOP_MESSAGE(ROUTINE,'VOLUME ????')
C
C=======================================================================
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI ) THEN
         CALL DRV_MPI_BARRIER
         CALL DRV_MPI_BCAST_I(0,NPANLIM,1)
         CALL DRV_MPI_BCAST_I(0,NRSFLIM,1)
      END IF
C
99001 FORMAT (10X,'INFO from <SFNCREATE>:  NSHLCLU_SFN set to ',I4)
99002 FORMAT (/,10X,'<SCFNCREATE>: volume of unit cell not properly',/,
     &        10X,'reproduced by integration over shape functions',/)
99003 FORMAT (/,10X,'volume of the atomic cells',//,10X,'IM',10X,
     &        'vol (poly)    vol (intd)',6X,'r_mt (used)   r_mt (fill)')
99004 FORMAT (I12,3X,2(3X,2F14.8))
99005 FORMAT (10X,'for IM =',I3,3X,'volume polyhedron:',F14.8,
     &        '  numerical:',F14.8,4X,F14.8)
99006 FORMAT (10X,'initial relative volume mismatch     ERR =',F12.8,/,
     &        10X,'shape functions scaled linearly with WH2 =',F12.8)
99007 FORMAT (10X,'SUM',F19.8,F14.8,/,10X,'VUC',F19.8,/)
99008 FORMAT (/,' ##### WARNING from <SFNCREATE> ',48('#'),/,10X,A)
99009 FORMAT (/,' ##### TROUBLE in <SFNCREATE> ',50('#'),/,10X,A)
99010 FORMAT (10X,A,10I5)
99011 FORMAT (/,1X,79('*'),/,:,/,10X,'mesh IM =',I3,/)
99012 FORMAT (//,1X,79('*'),/,34X,'<SCFNCREATE>',/,1X,79('*'),//,10X,
     &        'create the shape functions for ',I4,' meshes',/)
99013 FORMAT (/,' ***** INFO from <SFNCREATE> ',51('*'),/,10X,
     &        'for mesh IM=',I3,'   dummy shape function added',/)
99014 FORMAT (/,:,10X,'mesh type IM =',I4,10X,A,I4,4X,A,I5)
99015 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99016 FORMAT (/,10X,'trouble in <SFNVORONOI>: ',' trying NSHLCLU_SFN =',
     &        I3,/)
99017 FORMAT (3(/,1X,79('m')),/,:,/,30X,'mesh IM =',I4,/,3(/,1X,79('m'))
     &        )
99018 FORMAT (1(/,1X,79('x')),/,:,/,20X,'IM =',I4,4X,
     &        'calling subroutine ',A,/,1(/,1X,79('x')))
      END
C*==sfnrasmol.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SFNRASMOL(KEY_SFN_RASMOL,IMCURR,NQ,QBAS0,DQBAS,DATSET0,
     &                     LDATSET0,IMQ,NQCLU_SFN,RQCLU_SFN,IOTMP,
     &                     NFACE_M,NEDGE_FCM,RVERT_VFM,NSFTSYMQ,SYMTVEC,
     &                     MROTR,ISYMGENQ,IQREPQ,NVERTMAX,NFACEMAX,
     &                     NMMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *    create data and script files for  rasmol                      *
C   *                                                                  *
C   *    - cluster and polyhedron for each mesh                        *
C   *    - unit cell and polyhedra                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:A0_ANG
      USE MOD_LATTICE,ONLY:ABAS,ADAINV
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNRASMOL')
C
C Dummy arguments
C
      CHARACTER*80 DATSET0
      INTEGER IMCURR,IOTMP,KEY_SFN_RASMOL,LDATSET0,NFACEMAX,NMMAX,NQ,
     &        NQCLU_SFN,NQMAX,NVERTMAX
      REAL*8 DQBAS(3,NQMAX),MROTR(3,3,NSYMMAX),QBAS0(3,NQMAX),
     &       RQCLU_SFN(3,NQCLU_SFN),RVERT_VFM(3,NVERTMAX,NFACEMAX,NMMAX)
     &       ,SYMTVEC(3,NSYMMAX)
      INTEGER IMQ(NQMAX),IQREPQ(NQMAX),ISYMGENQ(NQMAX),
     &        NEDGE_FCM(NFACEMAX,NMMAX),NFACE_M(NMMAX),
     &        NSFTSYMQ(3,NSYMMAX,NQMAX)
C
C Local variables
C
      REAL*8 COLOR,PCVEC(3),QBASPI(3),QBASPII(3),S,SFT,VEC(3),VECP(3)
      CHARACTER*80 FILPDB,FILRAS
      INTEGER I,IC,IC0,IC1,ICP,IM,IPV,IQ,IQCLU_SFN,IQREP,ISYM,IX,J,LFIL,
     &        LL,NC,NPOINTS_SITES,NPOINTS_TOT,NPOINTS_UC
      CHARACTER*20 STR20
C
C*** End of declarations rewritten by SPAG
C
C----- use larger scaling factor to avoid spurious lines drawn by rasmol
      S = 3*A0_ANG*8D0
      S = 6*S
C
      IF ( KEY_SFN_RASMOL.EQ.1 ) THEN
C=======================================================================
C           plot 'shape function' cluster and central polyhedron
C=======================================================================
C
         IM = IMCURR
C
         FILPDB = DATSET0(1:LDATSET0)//'_SFN_M'
         CALL STRING_ADD_N(FILPDB,IM)
         LFIL = LEN_TRIM(FILPDB)
         FILPDB = FILPDB(1:LFIL)//'.pdb'
         LFIL = LFIL + 4
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILPDB(1:LFIL))
         WRITE (IOTMP,99005) 'cluster and polyhedron for '//
     &                       DATSET0(1:LDATSET0)
C
C ----------------------------------------------------------------------
C                  Plotting of the cluster
C ----------------------------------------------------------------------
         DO IQCLU_SFN = 1,NQCLU_SFN
C
            COLOR = 1.0D0
            WRITE (IOTMP,FMT=99001) IQCLU_SFN,IQCLU_SFN,
     &                              (RQCLU_SFN(IX,IQCLU_SFN)*S,IX=1,3),
     &                              COLOR
         END DO
C
C ----------------------------------------------------------------------
C                  Plotting of the polyhedron
C ----------------------------------------------------------------------
C
         COLOR = 2.0D0
         IC0 = NQCLU_SFN
         IC = IC0
         DO I = 1,NFACE_M(IM)
            DO J = 1,NEDGE_FCM(I,IM)
               IC = IC + 1
               WRITE (IOTMP,FMT=99002) IC,IC,
     &                                 (RVERT_VFM(IX,J,I,IM)*S,IX=1,3),
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
               WRITE (IOTMP,99004) IC,ICP
            END DO
         END DO
C
         WRITE (IOTMP,99003)
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
         CALL STRING_ADD_N(STR20,NQCLU_SFN)
         WRITE (IOTMP,*) STR20
         WRITE (IOTMP,*) 'cpk  150'
         STR20 = 'select '
         CALL STRING_ADD_N(STR20,NQCLU_SFN+1)
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
      ELSE
C=======================================================================
C KEY_SFN_RASMOL = 2:  plot unit cell and polyhedra around sites
C KEY_SFN_RASMOL = 3:  plot embedded cluster and polyhedra around sites
C                      no unit cell is plotted
C=======================================================================
C
         IF ( KEY_SFN_RASMOL.EQ.2 ) THEN
            FILPDB = DATSET0(1:LDATSET0)//'_SFN_UC'
            LFIL = LDATSET0 + 7
         ELSE
            FILPDB = DATSET0(1:LDATSET0)//'_SFN_CLU'
            LFIL = LDATSET0 + 8
         END IF
         FILPDB = FILPDB(1:LFIL)//'.pdb'
         LFIL = LFIL + 4
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILPDB(1:LFIL))
         WRITE (IOTMP,99005) 'unit cell and polyhedron for '//
     &                       DATSET0(1:LDATSET0)
C
C------------------------------------------ specify corners of unit cell
C
         IF ( KEY_SFN_RASMOL.EQ.2 ) THEN
C
            COLOR = 3D0
            WRITE (IOTMP,FMT=99002) 1,1,0D0,0D0,0D0,COLOR
            DO I = 1,3
               WRITE (IOTMP,FMT=99002) (I+1),(I+1),(ABAS(IX,I)*S,IX=1,3)
     &                                 ,COLOR
            END DO
C
            WRITE (IOTMP,FMT=99002) 5,5,S*(ABAS(1,1)+ABAS(1,2)),
     &                              S*(ABAS(2,1)+ABAS(2,2)),
     &                              S*(ABAS(3,1)+ABAS(3,2)),COLOR
            WRITE (IOTMP,FMT=99002) 6,6,S*(ABAS(1,1)+ABAS(1,3)),
     &                              S*(ABAS(2,1)+ABAS(2,3)),
     &                              S*(ABAS(3,1)+ABAS(3,3)),COLOR
            WRITE (IOTMP,FMT=99002) 7,7,S*(ABAS(1,3)+ABAS(1,2)),
     &                              S*(ABAS(2,3)+ABAS(2,2)),
     &                              S*(ABAS(3,3)+ABAS(3,2)),COLOR
            WRITE (IOTMP,FMT=99002) 8,8,
     &                              S*(ABAS(1,1)+ABAS(1,2)+ABAS(1,3)),
     &                              S*(ABAS(2,1)+ABAS(2,2)+ABAS(2,3)),
     &                              S*(ABAS(3,1)+ABAS(3,2)+ABAS(3,3)),
     &                              COLOR
C
            NPOINTS_UC = 8
            NPOINTS_SITES = 0
C
         ELSE
C
C-------------------------------------------- unshifted atomic positions
C
            COLOR = 3.0D0
            I = 0
            DO IQ = 1,NQ
               I = I + 1
               WRITE (IOTMP,99001) I,I,QBAS0(1:3,IQ)*S,COLOR
            END DO
C
            NPOINTS_UC = 0
            NPOINTS_SITES = NQ
C
         END IF
C
C---------------------------------------------- specify atomic positions
C
         COLOR = 1.0D0
         I = NPOINTS_UC + NPOINTS_SITES
         DO IQ = 1,NQ
            I = I + 1
            WRITE (IOTMP,99001) I,I,
     &                          ((QBAS0(IX,IQ)+DQBAS(IX,IQ))*S,IX=1,3),
     &                          COLOR
         END DO
         NPOINTS_SITES = NPOINTS_SITES + NQ
C
C-------------------------------------------------------- plot polyhedra
C
         COLOR = 2.0D0
         IC0 = NPOINTS_UC + NPOINTS_SITES
C
         DO IQ = 1,NQ
            ISYM = ISYMGENQ(IQ)
            IQREP = IQREPQ(IQ)
C
C- find the primitive shift of the atomic site due to symmetry operation
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,QBAS0(1,IQREP),1,
     &                 0D0,QBASPI,1)
            QBASPI(1:3) = QBASPI(1:3) + SYMTVEC(1:3,ISYM)
C
            QBASPII(1:3) = QBASPI(1:3) - QBAS0(1:3,IQ)
C
            CALL RVECEXPAND(QBASPII,ABAS,ADAINV,PCVEC)
            DO I = 1,3
               IF ( ABS(NINT(PCVEC(I))-PCVEC(I)).GT.1D-9 ) THEN
                  WRITE (6,*) '<SFNRASMOL> trouble with exp. coeff',I
                  WRITE (6,*) 'PCVEC  ',PCVEC
               END IF
            END DO
C
            IM = IMQ(IQ)
            IC = IC0
            DO I = 1,NFACE_M(IM)
               DO J = 1,NEDGE_FCM(I,IM)
                  IC = IC + 1
C
                  VEC(1:3) = QBAS0(1:3,IQREP) + RVERT_VFM(1:3,J,I,IM)
C
                  CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,VEC,1,0D0,
     &                       VECP,1)
C
                  VECP(1:3) = VECP(1:3) + SYMTVEC(1:3,ISYM)
C
                  DO IX = 1,3
                     SFT = SYMTVEC(IX,ISYM)
                     DO IPV = 1,3
                        SFT = SFT - ABAS(IX,IPV)
     &                        *NSFTSYMQ(IPV,ISYM,IQREP)
                     END DO
C                     VECP(IX) = VECP(IX) + SFT
                  END DO
                  VECP(1:3) = VECP(1:3) - QBASPII(1:3)
C
                  WRITE (IOTMP,FMT=99002) IC,IC,(VECP(IX)*S,IX=1,3),
     &                   COLOR
               END DO
            END DO
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
                  WRITE (IOTMP,99004) IC,ICP
               END DO
            END DO
C
            IC0 = IC
         END DO
         NPOINTS_TOT = IC
C
         IF ( KEY_SFN_RASMOL.EQ.2 ) THEN
            WRITE (IOTMP,FMT=99007)
         ELSE
            WRITE (IOTMP,'(''END'')')
         END IF
         CLOSE (IOTMP)
C
C--------------------------------------------------- write RASMOL script
C
         FILRAS = FILPDB(1:(LFIL-4))//'.ras'
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILRAS(1:LFIL))
C
         WRITE (IOTMP,*) 'load '''//FILPDB(1:LFIL)//'''  '
         WRITE (IOTMP,*) 'set background white'
         WRITE (IOTMP,*) 'color temperature'
         WRITE (IOTMP,*) 'set axes on '
C
         IF ( KEY_SFN_RASMOL.EQ.2 ) THEN
            WRITE (IOTMP,99008) 1,NPOINTS_UC
            WRITE (IOTMP,*) 'cpk  0'
         END IF
         WRITE (IOTMP,99008) (NPOINTS_UC+1),(NPOINTS_UC+NPOINTS_SITES)
         WRITE (IOTMP,*) 'cpk 150'
         WRITE (IOTMP,99008) (NPOINTS_UC+NPOINTS_SITES+1),NPOINTS_TOT
         WRITE (IOTMP,*) 'cpk  0'
C
         CLOSE (IOTMP)
C
         WRITE (6,99006) FILPDB(1:LFIL),FILRAS(1:LFIL)
C
      END IF
C=======================================================================
C
99001 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99002 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99003 FORMAT ('END  ')
99004 FORMAT ('CONECT',2I5)
99005 FORMAT ('HEADER    ',A,/,'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
99006 FORMAT (/,5X,'unit cell data stored in rasmol data-file ',A,/,5X,
     &        'view via:   rasmol  -script ',A)
99007 FORMAT ('CONECT    1    2',/,'CONECT    1    3',/,
     &        'CONECT    1    4',/,'CONECT    5    2',/,
     &        'CONECT    5    3',/,'CONECT    5    8',/,
     &        'CONECT    6    2',/,'CONECT    6    4',/,
     &        'CONECT    6    8',/,'CONECT    7    3',/,
     &        'CONECT    7    4',/,'CONECT    7    8',/,'END')
99008 FORMAT ('SELECT ',I5,'-',I5)
      END
C*==sfnlebedev.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SFNLEBEDEV(NFACE_M,NEDGE_FCM,RVERT_VFM,NVERTMAX,
     &                      NFACEMAX)
C   ********************************************************************
C   *                                                                  *
C   *  find the intersections of the direction vectors connected       *
C   *  with the Lebedev grid with the faces of an atomic cell          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:A0_ANG
      USE MOD_FILES,ONLY:IPRINT,DATSET0,LDATSET0,IOTMP
      USE MOD_RMESH,ONLY:NM,NMMAX,RHAT_LEBGRID,D_LGM,N_LEBGRID
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNLEBEDEV')
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      INTEGER NFACEMAX,NVERTMAX
      INTEGER NEDGE_FCM(NFACEMAX,NMMAX),NFACE_M(NMMAX)
      REAL*8 RVERT_VFM(3,NVERTMAX,NFACEMAX,NMMAX)
C
C Local variables
C
      REAL*8 AA,AB,ALFA,AVEC(3),BB,BETA,BVEC(3),COLOR,D,DET,
     &       DVEC_LEBGRID(:,:),H,HI,NHAT(3),N_DOT_D,RA,RB,RSUM,RVEC(3),
     &       S,VEC(3),XHIGH,XLOW
      LOGICAL ALL_DONE,DONE_LEBGRID(:)
      REAL*8 DDOT,DNRM2
      CHARACTER*80 FILPDB,FILRAS
      INTEGER I,IC,IC0,IC1,ICP,IFACE,IM,IVERT,I_LEBGRID,J,LFIL,LL,NC
      CHARACTER*20 STR20
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DONE_LEBGRID,DVEC_LEBGRID
C
      ALLOCATE (D_LGM(N_LEBGRID,NMMAX))
      ALLOCATE (DONE_LEBGRID(N_LEBGRID),DVEC_LEBGRID(3,N_LEBGRID))
C
      DONE_LEBGRID(:) = .FALSE.
C----- use larger scaling factor to avoid spurious lines drawn by rasmol
      S = 3*A0_ANG*8D0
      S = 6*S
C
Cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      DO IM = 1,NM
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
         DO IFACE = 1,NFACE_M(IM)
C
            IF ( IPRINT.GT.0 ) WRITE (6,99001) IM,IFACE
C
C-----------------------------------------------------------------------
C              find normal vector NHAT for surface IFACE
C-----------------------------------------------------------------------
C
            AVEC(1:3) = RVERT_VFM(1:3,1,IFACE,IM)
     &                  - RVERT_VFM(1:3,2,IFACE,IM)
            BVEC(1:3) = RVERT_VFM(1:3,3,IFACE,IM)
     &                  - RVERT_VFM(1:3,1,IFACE,IM)
C
            CALL RVECXPRO(AVEC,BVEC,RVEC)
C
            NHAT(1:3) = RVEC(1:3)/DNRM2(3,RVEC,1)
            H = DDOT(3,NHAT,1,RVERT_VFM(1,1,IFACE,IM),1)
C
            IF ( H.LT.0D0 ) THEN
               NHAT(1:3) = -NHAT(1:3)
               H = -H
            END IF
C
C--------------- check that all corner vectors fullfil   n^ * ->c_i = h
            DO IVERT = 1,NEDGE_FCM(IFACE,IM)
C
               HI = DDOT(3,NHAT,1,RVERT_VFM(1,IVERT,IFACE,IM),1)
C
               IF ( ABS(1D0-HI/H).GT.1D-10 ) WRITE (6,99002) IVERT,H,
     &              (1D0-HI/H)
C
            END DO
C
Clllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
            DO I_LEBGRID = 1,N_LEBGRID
               IF ( .NOT.DONE_LEBGRID(I_LEBGRID) ) THEN
C
                  N_DOT_D = DDOT(3,NHAT,1,RHAT_LEBGRID(1,I_LEBGRID),1)
C
C-----------------------------------------------------------------------
                  IF ( N_DOT_D.GT.0D0 ) THEN
C
                     D = H/N_DOT_D
C
                     RVEC(1:3) = D*RHAT_LEBGRID(1:3,I_LEBGRID)
     &                           - RVERT_VFM(1:3,1,IFACE,IM)
C
C--------------------------- run over all trangles with corner at ->c_1
C
                     DO IVERT = 2,NEDGE_FCM(IFACE,IM) - 1
C
C---------- expand  (->d - ->c_1)  in terms of edge vectors ->a and ->b
C
                        AVEC(1:3) = RVERT_VFM(1:3,IVERT,IFACE,IM)
     &                              - RVERT_VFM(1:3,1,IFACE,IM)
                        BVEC(1:3) = RVERT_VFM(1:3,IVERT+1,IFACE,IM)
     &                              - RVERT_VFM(1:3,1,IFACE,IM)
C
                        AA = DDOT(3,AVEC,1,AVEC,1)
                        AB = DDOT(3,AVEC,1,BVEC,1)
                        BB = DDOT(3,BVEC,1,BVEC,1)
                        RA = DDOT(3,RVEC,1,AVEC,1)
                        RB = DDOT(3,RVEC,1,BVEC,1)
                        DET = AA*BB - AB*AB
                        ALFA = (BB*RA-AB*RB)/DET
                        BETA = (-AB*RA+AA*RB)/DET
C
C-------------------------------- check the expansion for (->d - ->c_1)
C
                        VEC(1:3) = ALFA*AVEC(1:3) + BETA*BVEC(1:3)
C
                        RSUM = 0D0
                        DO I = 1,3
                           RSUM = RSUM + ABS(RVEC(I)-VEC(I))
                        END DO
                        IF ( ABS(RSUM).GT.TOL ) WRITE (6,*) IVERT,
     &                       RVEC(1:3),VEC(1:3)
C
C check whether (->d - ->c_1) is in the triangle spanned by ->a and ->b
C
                        XLOW = -TOL
                        XHIGH = 1D0 + TOL
                        IF ( XLOW.LT.ALFA .AND. ALFA.LT.XHIGH .AND. 
     &                       XLOW.LE.BETA .AND. BETA.LE.XHIGH .AND. 
     &                       (ALFA+BETA).LE.XHIGH ) THEN
C
                           DONE_LEBGRID(I_LEBGRID) = .TRUE.
                           D_LGM(I_LEBGRID,IM) = D
C
                           DVEC_LEBGRID(1:3,I_LEBGRID)
     &                        = D*RHAT_LEBGRID(1:3,I_LEBGRID)
C
                           IF ( IPRINT.GT.0 ) WRITE (6,*)
     &                           'LEB direction ',I_LEBGRID,' crosses ',
     &                          'surface ',IFACE
C
                           EXIT
C
                        END IF
C
                     END DO
C
                  END IF
C-----------------------------------------------------------------------
C
               END IF
            END DO
Clllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
C
         END DO
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
         ALL_DONE = .TRUE.
         DO I_LEBGRID = 1,N_LEBGRID
            IF ( .NOT.DONE_LEBGRID(I_LEBGRID) ) THEN
               WRITE (6,*) 'NO CROSSING FOUND FOR LEB direction ',
     &                     I_LEBGRID
               ALL_DONE = .FALSE.
            END IF
         END DO
         IF ( .NOT.ALL_DONE ) CALL STOP_MESSAGE(ROUTINE,
     &        'crossings not found for ALL directions')
C
C=======================================================================
C           plot central polyhedron + intersection points
C=======================================================================
C
         FILPDB = DATSET0(1:LDATSET0)//'_LGD_M'
         CALL STRING_ADD_N(FILPDB,IM)
         LFIL = LEN_TRIM(FILPDB)
         FILPDB = FILPDB(1:LFIL)//'.pdb'
         LFIL = LFIL + 4
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILPDB(1:LFIL))
         WRITE (IOTMP,99008) 'polyhedron + Lebedev grid for '//
     &                       DATSET0(1:LDATSET0)
C
C ----------------------------------------------------------------------
C                  Plotting of the intersection points of
C            Lebedev grid and polyhedron surface
C ----------------------------------------------------------------------
C
         DO I_LEBGRID = 1,N_LEBGRID
C
            COLOR = 1.0D0
            WRITE (IOTMP,FMT=99004) I_LEBGRID,I_LEBGRID,
     &                              DVEC_LEBGRID(1:3,I_LEBGRID)*S,COLOR
         END DO
C
C ----------------------------------------------------------------------
C                  Plotting of the polyhedron
C ----------------------------------------------------------------------
C
         COLOR = 2.0D0
         IC0 = N_LEBGRID
         IC = IC0
         DO I = 1,NFACE_M(IM)
            DO J = 1,NEDGE_FCM(I,IM)
               IC = IC + 1
               WRITE (IOTMP,FMT=99005) IC,IC,RVERT_VFM(1:3,J,I,IM)*S,
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
               WRITE (IOTMP,99007) IC,ICP
            END DO
         END DO
C
         WRITE (IOTMP,99006)
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
         CALL STRING_ADD_N(STR20,N_LEBGRID)
         WRITE (IOTMP,*) STR20
         WRITE (IOTMP,*) 'cpk  150'
         STR20 = 'select '
         CALL STRING_ADD_N(STR20,N_LEBGRID+1)
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
      WRITE (6,99003) ROUTINE(1:LEN_TRIM(ROUTINE))
C
99001 FORMAT (/,'MESH ',I3,'   FACE ',I3)
99002 FORMAT ('VERTEX ',I3,'   DISTANCE and rel error:',2F20.14)
99003 FORMAT (/,1X,79('*'),/,36X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'for all Lebedev directions crossing with surface found',
     &        /)
99004 FORMAT ('ATOM  ',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99005 FORMAT ('HETATM',I5,'          ',I5,'    ',3F8.3,'  0.00',F8.3)
99006 FORMAT ('END  ')
99007 FORMAT ('CONECT',2I5)
99008 FORMAT ('HEADER    ',A,/,'SOURCE    SPRKKR - program       ',/,
     &        'AUTHOR    H. Ebert               ',/,
     &        'REMARK    None                   ')
      END
