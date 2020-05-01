C*==sfncreateclu.f    processed by SPAG 6.70Rc at 14:57 on 28 Apr 2017
      SUBROUTINE SFNCREATECLU
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
C OS: parameters increased  16/05/31
C
C     PARAMETER (NVERTMAX=150,NFACEMAX=200,NPANMAX_S=300,NRSFMAX_S=335)
C
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_TYPES,ONLY:NTMAX,NAT,ITBOT,ITTOP,NLMFPMAX
      USE MOD_ANGMOM,ONLY:NL
      USE MOD_RMESH,ONLY:NRSFMAX,NLMSFMAX,NSFMAX,NPANMAX,NMMAX,LMISF,
     &    NSF,KLMSF,ISFLM,FLMSF,NPAN,NRSFTOT,NM,NMHOST,NMCLU,JRCUT,
     &    WINTLM,R,DX,JRNS1,RMT,JRMT
      USE MOD_FILES,ONLY:DATSET0,LDATSET0,IPRINT,SFNFIL,LSFNFIL,IOTMP,
     &    RDUMMY,IDUMMY
      USE MOD_SITES,ONLY:IQAT,ITOQ,QBAS,NQ,NQMAX,NQ_L,NQ_R,IMQ,IQBOT,
     &    IQTOP,NQHOST,IQ_QCLU,QBAS0_QCLU,DQBAS_QCLU,NQCLU
      USE MOD_LATTICE,ONLY:ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV_L,ADAINV_I,
     &    ADAINV_R,SYSTEM_DIMENSION,VOLUC,ALAT
      USE MOD_SYMMETRY,ONLY:NSFTSYMQ,SYMTVEC,MROTR,ISYMGENQ,IQREPQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SFNCREATECLU')
      INTEGER NVERTMAX,NFACEMAX,NPANMAX_S,NRSFMAX_S
      PARAMETER (NVERTMAX=300,NFACEMAX=1600,NPANMAX_S=600,
     &           NRSFMAX_S=2000)
      LOGICAL CHECK_SYM
      PARAMETER (CHECK_SYM=.FALSE.)
C
C Local variables
C
      REAL*8 A3(:),B3(:),C3(:),CLURAD_SFN,D3(:),DLT,DQBAS(:,:),
     &       DRN_S(:,:),DRSF(:,:),FLMSF_PR(:,:),FLMSF_S(:,:,:),RMTRED0,
     &       RMTREDFILL,ROUT,RQCLU_SFN(:,:),RVERT_VFM(:,:,:,:),SCALE_S,
     &       SCLM(:),SIZEFAC(:),VOLUME,VOL_M(:),WEIGHT0,WEIGHTCLU_SFN(:)
     &       ,WINTLM_HOST(:,:,:),XEDGE(:,:),XRN_S(:,:),XRSF(:,:),
     &       YEDGE(:,:),ZEDGE(:,:)
      LOGICAL ACCEPT,KEEP_SFN_S,MOL,TEST
      INTEGER I,IA_ERR,IED,IFC,IFLAG,IFLAG_POLYHEDRON,IFLAG_VERTEX,IM,
     &        IMBOT,IMCURR,IMHOST,IMHOST_M(:),IMTOP,IPROC,IPROCM(:),IQ,
     &        IQ1,IQCLU,IQCLU_SFN,IQCNTR,IQHOST,IQIMP,IQ_MAUX(:),
     &        IQ_QCLU_SFN(:),ISF,ISFLM_PR(:,:),ISF_S,IT,IWSF,J,
     &        JRCUT_HOST(:,:),KEY_PAN,KFP_LMQ_PR(:,:),KLMFP_PR(:,:),
     &        KLMSF_HOST(:,:),KLMSF_PR(:,:),LM,LMIFP_PR(:,:),
     &        LMISF_PR(:,:),LMISF_S(:,:),LSF,LSFMAX,LTXT_PR(:),
     &        MESHN_S(:),NCELL_S,NEDGE(:),NEDGE_FCM(:,:),NFACE,
     &        NFACE_M(:),NFPT_PR(:),NFUN_S(:),NLMFPT_PR(:),NLMSFMAXHOST,
     &        NMSF(:),NM_S(:,:),NPANLIM,NPANMAXHOST,NPAN_S(:),
     &        NQCLU_I_SFN,NQCLU_L_SFN,NQCLU_R_SFN,NQCLU_SFN,NRPAN(:,:),
     &        NRSF0,NRSFLIM,NRSFMAXHOST,NSFLIM,NSFMAXHOST,NSF_PR(:),
     &        NSF_PRMAX,NSHLCLU_SFN,NWARN
      CHARACTER*8 TXT_PR(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA DLT/0.05D0/,NRSF0/200/,TEST/.TRUE./
      DATA NQCLU_SFN/0/,NQCLU_L_SFN/0/,NQCLU_I_SFN/0/,NQCLU_R_SFN/0/
C     ------------------------------------------------------------------
C
      ALLOCATABLE JRCUT_HOST,KLMSF_HOST,WINTLM_HOST,IMHOST_M,IPROCM
      ALLOCATABLE DRSF,SCLM,XRSF,NRPAN,NPAN_S,NFUN_S
      ALLOCATABLE SIZEFAC,VOL_M,DQBAS
      ALLOCATABLE RQCLU_SFN,IQ_QCLU_SFN,WEIGHTCLU_SFN
      ALLOCATABLE A3,B3,C3,D3,XEDGE,YEDGE,ZEDGE,NEDGE
      ALLOCATABLE RVERT_VFM,NEDGE_FCM
      ALLOCATABLE NMSF,NM_S,MESHN_S,LMISF_S,XRN_S
      ALLOCATABLE DRN_S,FLMSF_S,NFACE_M,FLMSF_PR,KFP_LMQ_PR
      ALLOCATABLE KLMFP_PR,NLMFPT_PR,NFPT_PR,LMIFP_PR,ISFLM_PR
      ALLOCATABLE KLMSF_PR,NSF_PR,LMISF_PR,TXT_PR,LTXT_PR,IQ_MAUX
C
      LSF = 4*(NL-1)
      LSFMAX = LSF
      IF ( NLMSFMAX.NE.(LSFMAX+1)*(LSFMAX+1) )
     &     CALL STOP_MESSAGE(ROUTINE,'LSFMAX  inconsistent')
C
      ALLOCATE (IMHOST_M(NMMAX))
      IMHOST_M(1:NMMAX) = 0
      ALLOCATE (DRSF(NRSFMAX,NMMAX),SCLM(NMMAX))
      ALLOCATE (XRSF(NRSFMAX,NMMAX),NRPAN(NPANMAX,NMMAX))
      ALLOCATE (SIZEFAC(NTMAX),NFACE_M(NMMAX),VOL_M(NMMAX))
      ALLOCATE (A3(NFACEMAX),B3(NFACEMAX),IPROCM(NMMAX))
      ALLOCATE (C3(NFACEMAX),D3(NFACEMAX))
      ALLOCATE (RVERT_VFM(3,NVERTMAX,NFACEMAX,NMMAX))
      ALLOCATE (XEDGE(NVERTMAX,NFACEMAX))
      ALLOCATE (YEDGE(NVERTMAX,NFACEMAX))
      ALLOCATE (ZEDGE(NVERTMAX,NFACEMAX))
      ALLOCATE (NEDGE_FCM(NFACEMAX,NMMAX),NMSF(NPANMAX_S))
      ALLOCATE (NEDGE(NFACEMAX),NPAN_S(NMMAX),NFUN_S(NMMAX))
      ALLOCATE (NM_S(NPANMAX_S,NMMAX),MESHN_S(NMMAX))
      ALLOCATE (LMISF_S(NLMSFMAX,NMMAX))
      ALLOCATE (DRN_S(NRSFMAX_S,NMMAX),XRN_S(NRSFMAX_S,NMMAX))
      ALLOCATE (FLMSF_S(NRSFMAX_S,NLMSFMAX,NMMAX))
      ALLOCATE (FLMSF_PR(NRSFMAX_S,NLMSFMAX),DQBAS(3,NQMAX))
      DQBAS(1:3,1:NQMAX) = 0D0
C
      MOL = .FALSE.
C
C-------------------------------------------------------- dummy settings
      SIZEFAC(1:NTMAX) = 1D0
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
      IM = NMHOST
      DO IQ = IQBOT,IQTOP
         IF ( IQ.EQ.IQREPQ(IQ) ) THEN
            IM = IM + 1
            IMQ(IQ) = IM
         END IF
         IMQ(IQ) = IMQ(IQREPQ(IQ))
      END DO
      IF ( (IM-NMHOST).NE.NMCLU )
     &      CALL STOP_MESSAGE(ROUTINE,'(IM-NMHOST) .NE. NMCLU')
C
      IMBOT = NMMAX
      IMTOP = 0
      DO IQ = IQBOT,IQTOP
         IMBOT = MIN(IMBOT,IMQ(IQ))
         IMTOP = MAX(IMTOP,IMQ(IQ))
      END DO
      IF ( IMBOT.NE.NMHOST+1 )
     &      CALL STOP_MESSAGE(ROUTINE,'IMBOT.NE.NMHOST+1')
      IF ( IMTOP.NE.NMHOST+NMCLU )
     &      CALL STOP_MESSAGE(ROUTINE,'IMBOT.NE.NMHOST+NMCLU')
C
      DO IT = 1,NTMAX
         TXT_PR(IT) = '  '
         LTXT_PR(IT) = 2
      END DO
      DO IM = IMBOT,IMTOP
         IQ_MAUX(IM) = 0
         DO IT = ITBOT,ITTOP
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
C-----------------------------------------------------------------------
C                    reread shape functions for HOST
C-----------------------------------------------------------------------
C
      CALL SFNREAD(NMHOST,NRPAN,XRSF,DRSF,SCLM)
C
C-----------------------------------------------------------------------
C
      WRITE (6,99009) IMBOT,IMTOP,SFNFIL(1:LSFNFIL)
C
      KEY_PAN = 0
      NPANLIM = 0
      NRSFLIM = 0
      NSFLIM = 0
C
      IFLAG = 0
      NWARN = 0
C
C-----------------------------------------------------------------------
C   now create additional shape functions for the CLUSTER
C        and add to the shape functions file
C-----------------------------------------------------------------------
C
      CALL MPI_DISTRIBUTE(IPROCM(IMBOT),(IMTOP-IMBOT+1),MPI,'M')
C
C=======================================================================
      DO IM = IMBOT,IMTOP
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCM(IM) ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            NQCLU_SFN = 0
            NQCLU_L_SFN = 0
            NQCLU_I_SFN = 0
            NQCLU_R_SFN = 0
            LMISF_S(1:NLMSFMAX,IM) = 0
            XRN_S(1:NRSFMAX_S,IM) = 0D0
            DRN_S(1:NRSFMAX_S,IM) = 0D0
            FLMSF_S(1:NRSFMAX_S,1:NLMSFMAX,IM) = 0D0
            NM_S(1:NPANMAX_S,IM) = 0
C
            IQIMP = IQ_MAUX(IM)
            IQCLU = IQ - NQHOST
            IQHOST = IQ_QCLU(IQCLU)
C
            IMHOST = IMQ(IQHOST)
            IMHOST_M(IM) = IMHOST
C
            IQCNTR = IQHOST
            NSHLCLU_SFN = 10
            CLURAD_SFN = 0D0
C
C???????????????????????????????????????????????????????????????????????
            IF ( TEST ) THEN
               WRITE (*,*) '************* IM        ',IM
               WRITE (*,*) '************* IMHOST    ',IMHOST
               WRITE (*,*) '************* IQIMP     ',IQIMP
               WRITE (*,*) '************* IQCLU     ',IQCLU
               WRITE (*,*) '************* IQHOST    ',IQHOST
               WRITE (*,*) '************* IQCNTR    ',IQCNTR
            END IF
C???????????????????????????????????????????????????????????????????????
C
            CALL CLUSSITES(IOTMP,IPRINT,MOL,SYSTEM_DIMENSION,ABAS,
     &                     ABAS_L,ABAS_I,ABAS_R,ADAINV_L,ADAINV_I,
     &                     ADAINV_R,QBAS,CLURAD_SFN,IQCNTR,NQCLU_SFN,
     &                     NQCLU_L_SFN,NQCLU_I_SFN,NQCLU_R_SFN,
     &                     NSHLCLU_SFN,NQHOST,NQ_L,NQ_R,NQMAX)
C
            IF ( IM.GT.IMBOT ) THEN
               DEALLOCATE (RQCLU_SFN,IQ_QCLU_SFN)
               DEALLOCATE (WEIGHTCLU_SFN)
            END IF
            ALLOCATE (RQCLU_SFN(3,NQCLU_SFN),WEIGHTCLU_SFN(NQCLU_SFN))
            ALLOCATE (IQ_QCLU_SFN(NQCLU_SFN),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &           CALL STOP_MESSAGE(ROUTINE,'ALLOC: IQCLUS')
            WEIGHTCLU_SFN(1:NQCLU_SFN) = 0D0
C
            READ (IOTMP) ((RQCLU_SFN(J,I),J=1,3),RDUMMY,IQ_QCLU_SFN(I),
     &                   I=1,NQCLU_SFN),(IDUMMY,I=1,NSHLCLU_SFN)
            CLOSE (IOTMP)
C
            WRITE (6,99008)
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
            IT = 1
            WEIGHT0 = SIZEFAC(IT)
C
            DO IQCLU_SFN = 1,NQCLU_SFN
               IT = ITOQ(1,IQ_QCLU_SFN(IQCLU_SFN))
               WEIGHTCLU_SFN(IQCLU_SFN) = SIZEFAC(IT)
            END DO
C
            CALL SFNVORONOI(NQCLU_SFN,RQCLU_SFN,NVERTMAX,WEIGHT0,
     &                      WEIGHTCLU_SFN,RMTREDFILL,ROUT,VOLUME,NFACE,
     &                      A3,B3,C3,D3,NEDGE,XEDGE,YEDGE,ZEDGE,
     &                      NFACEMAX,IFLAG_VERTEX)
C
C-----------------------------------------------------------------------
            IF ( IFLAG_VERTEX.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'IFLAG_VERTEX <> 0')
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
            CALL SFNRASMOL(1,IMCURR,NQ,QBAS,DQBAS,DATSET0,LDATSET0,IMQ,
     &                     NQCLU_SFN,RQCLU_SFN,IOTMP,NFACE_M,NEDGE_FCM,
     &                     RVERT_VFM,NSFTSYMQ,SYMTVEC,MROTR,ISYMGENQ,
     &                     IQREPQ,NVERTMAX,NFACEMAX,NMMAX,NQMAX)
C
C---------------------- Calculate shape functions for each Voronoi shape
C
            NMSF(1:NPANMAX_S) = 0
C
            CALL SFNSHAPE(NRSF0,A3,B3,C3,D3,NEDGE,XEDGE,YEDGE,ZEDGE,
     &                    NFACE,LSF,DLT,KEY_PAN,NMSF,NCELL_S,SCALE_S,
     &                    NPAN_S(IM),MESHN_S(IM),NM_S(1,IM),XRN_S(1,IM),
     &                    DRN_S(1,IM),NFUN_S(IM),LMISF_S(1,IM),
     &                    FLMSF_S(1,1,IM),NPANMAX_S,NRSFMAX_S,LSFMAX,
     &                    NLMSFMAX,NFACEMAX,NVERTMAX,IFLAG_POLYHEDRON)
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( TEST ) THEN
               NSF_PR(IM) = NFUN_S(IM)
               LMISF_PR(1:NFUN_S(IM),IM) = LMISF_S(1:NFUN_S(IM),IM)
               WRITE (6,99010) IM
            END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C--CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM
            IF ( CHECK_SYM ) THEN
C
C--------------------------------------------------- check picking rules
C
               IF ( NSF_PR(IM).NE.NFUN_S(IM) ) THEN
                  WRITE (6,99005) 'number of non-0 terms inconsistent'
                  WRITE (6,99007) 'IM               ',IM
                  WRITE (6,99007) 'NSF picking rules',NSF_PR(IM)
                  WRITE (6,99007) 'NSF in <SFNSHAPE>',NFUN_S(IM)
                  IFLAG = MAX(1,IFLAG)
                  NWARN = NWARN + 1
                  IF ( NFUN_S(IM).GT.NSF_PR(IM) ) IFLAG = MAX(2,IFLAG)
               END IF
C
               ACCEPT = NSF_PR(IM).GE.NFUN_S(IM)
               DO ISF = 1,NFUN_S(IM)
                  LM = LMISF_S(ISF,IM)
                  IF ( KLMSF_PR(LM,IM).NE.1 ) THEN
                     WRITE (6,99006) 
     &                             'set of non-0 terms inconsistent for'
                     WRITE (6,99007) 'IM       LM      ',IM,LM
                     IFLAG = MAX(2,IFLAG)
                     ACCEPT = .FALSE.
                  END IF
               END DO
               IF ( IFLAG.GE.1 ) THEN
                  WRITE (6,99007) 'KLMSF  from picking rules'
                  WRITE (6,'(10X,20I5)')
     &                   (LMISF_PR(ISF,IM),ISF=1,NSF_PR(IM))
                  WRITE (6,99007) 
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
                     WRITE (6,99010) IM
                  END IF
C
               END IF
            END IF
C--CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM-CHECK_SYM
C
            RMTRED0 = MIN(RMT(IMHOST)/ALAT,RMTREDFILL)
C
            CALL SFNMTMESH(IM,NPAN_S(IM),MESHN_S(IM),NM_S(1,IM),
     &                     XRN_S(1,IM),DRN_S(1,IM),NFUN_S(IM),
     &                     FLMSF_S(1,1,IM),RMTREDFILL,RMTRED0,NRSFMAX_S,
     &                     NPANMAX_S,NLMSFMAX)
C
            RMT(IM) = RMTRED0*ALAT
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
C collect: FLMSF_S LMISF_S  NFACE_M NEDGE_FCM RVERT_VFM  RMT
C
      IF ( MPI ) THEN
C=======================================================================
         DO IM = IMBOT,IMTOP
C
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
            CALL DRV_MPI_SEND_R(RMT(IM),1,IPROC,12)
C
         END DO
C=======================================================================
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_ID.EQ.0 ) THEN
C
C=======================================================================
         NPANLIM = 0
         NRSFLIM = 0
         NSFLIM = 0
C
         IWSF = 20
C
         OPEN (IWSF,FILE=SFNFIL(1:LSFNFIL))
C
         WRITE (IWSF,FMT='(I5)') NM
         DO IM = 1,NM
            WRITE (IWSF,FMT='(D20.12)') 1D0
         END DO
C
C-----------------------------------------------------------------------
C   rewrite shape functions for the HOST -- note shift in panel index
C-----------------------------------------------------------------------
C
         DO IM = 1,NMHOST
C
            CALL SFNWRITE(IWSF,(NPAN(IM)-1),NRSFTOT(IM),NRPAN(2,IM),
     &                    XRSF(1,IM),DRSF(1,IM),NSF(IM),FLMSF(1,1,IM),
     &                    LMISF(1,IM),NRSFMAX,NPANMAX,NLMSFMAX)
C
         END DO
C
         DO IM = IMBOT,IMTOP
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
C=======================================================================
C
         WRITE (6,99001)
         VOLUME = 0D0
         DO IQ = IQBOT,IQTOP
            IM = IMQ(IQ)
            WRITE (6,99002) IQ,IM,VOL_M(IM),RMTRED0*ALAT
            VOLUME = VOLUME + VOL_M(IM)
         END DO
         WRITE (6,99003) VOLUME,VOLUC*ALAT**3
C=======================================================================
C
         IQ1 = NQHOST + 1
C
         CALL SFNRASMOL(3,IMCURR,NQCLU,QBAS0_QCLU,DQBAS_QCLU,DATSET0,
     &                  LDATSET0,IMQ(IQ1),NQCLU_SFN,RQCLU_SFN,IOTMP,
     &                  NFACE_M,NEDGE_FCM,RVERT_VFM,NSFTSYMQ(1,1,IQ1),
     &                  SYMTVEC,MROTR,ISYMGENQ(IQ1),IQREPQ(IQ1),
     &                  NVERTMAX,NFACEMAX,NMMAX,NQMAX)
C
         WRITE (6,99008)
         WRITE (6,99007) 'running <SCFNCREATECLU> '
         WRITE (6,99007) 'mumber of warnings   NWARN =',NWARN
         WRITE (6,99007) 'error flag           IFLAG =',IFLAG
         WRITE (6,99008)
         IF ( IFLAG.EQ.2 )
     &         CALL STOP_MESSAGE(ROUTINE,'IFLAG.EQ.2: INCONSISTENCIES')
C
C-----------------------------------------------------------------------
C             reallocate storage and reread shape functions
C-----------------------------------------------------------------------
C
         DEALLOCATE (DRSF,SCLM,XRSF,NRPAN)
C
         ALLOCATE (DRSF(NRSFMAX,NMMAX),SCLM(NMMAX))
         ALLOCATE (XRSF(NRSFMAX,NMMAX),NRPAN(NPANMAX,NMMAX))
C
C
C-------------------------------------- variables depending on NM and SF
C
         NPANMAXHOST = NPANMAX
         NLMSFMAXHOST = NLMSFMAX
         NRSFMAXHOST = NRSFMAX
         NSFMAXHOST = NSFMAX
C
         ALLOCATE (JRCUT_HOST(0:NPANMAXHOST,NMHOST))
         JRCUT_HOST(0:NPANMAXHOST,1:NMHOST)
     &      = JRCUT(0:NPANMAXHOST,1:NMHOST)
         ALLOCATE (KLMSF_HOST(NLMSFMAXHOST,NMHOST))
         KLMSF_HOST(1:NLMSFMAXHOST,1:NMHOST)
     &      = KLMSF(1:NLMSFMAXHOST,1:NMHOST)
         ALLOCATE (WINTLM_HOST(1:NRSFMAXHOST,1:NSFMAXHOST,1:NMHOST))
         WINTLM_HOST(1:NRSFMAXHOST,1:NSFMAXHOST,1:NMHOST)
     &      = WINTLM(1:NRSFMAXHOST,1:NSFMAXHOST,1:NMHOST)
C
         NPANMAX = MAX(NPANLIM,NPANMAX)
         NRSFMAX = MAX(NRSFLIM,NRSFMAX)
         NSFMAX = MAX(NSFLIM,NSFMAX)
C----------------- keep NMMAX to avoid reallocation and swapping of data
C     NMMAX = NM
C
         WRITE (6,99004) NPANMAX,NRSFMAX,NSFMAX
C
         WRITE (*,*) '       NMhost ',NMHOST
         WRITE (*,*) '       NMclu  ',NMCLU
         WRITE (*,*) '       NM     ',NM
         WRITE (*,*) '       NMmax  ',NMMAX
C
         DEALLOCATE (LMISF,FLMSF)
         DEALLOCATE (JRCUT,WINTLM,ISFLM,KLMSF)
C
         ALLOCATE (FLMSF(NRSFMAX,NSFMAX,NMMAX),JRCUT(0:NPANMAX,NMMAX))
         ALLOCATE (WINTLM(NRSFMAX,NSFMAX,NMMAX),ISFLM(NLMSFMAX,NMMAX))
         ALLOCATE (KLMSF(NLMSFMAX,NMMAX),LMISF(NLMSFMAX,NMMAX))
C
         JRCUT(0:NPANMAXHOST,1:NMHOST)
     &      = JRCUT_HOST(0:NPANMAXHOST,1:NMHOST)
         KLMSF(1:NLMSFMAXHOST,1:NMHOST)
     &      = KLMSF_HOST(1:NLMSFMAXHOST,1:NMHOST)
         WINTLM(1:NRSFMAXHOST,1:NSFMAXHOST,1:NMHOST)
     &      = WINTLM_HOST(1:NRSFMAXHOST,1:NSFMAXHOST,1:NMHOST)
C
         DEALLOCATE (JRCUT_HOST,KLMSF_HOST,WINTLM_HOST)
C
         CALL SFNREAD(NM,NRPAN,XRSF,DRSF,SCLM)
C
C ----------------------------------------------------------------------
C
         DO IM = NMHOST + 1,NMHOST + NMCLU
            IMHOST = IMHOST_M(IM)
            R(1,IM) = R(1,IMHOST)
            JRNS1(IM) = JRNS1(IMHOST)
            JRMT(IM) = NINT(LOG(RMT(IM)/R(1,IM))/DX(IMHOST)+1)
            DX(IM) = LOG(RMT(IM)/R(1,IM))/DBLE(JRMT(IM)-1)
         END DO
C
         CALL RMESHFP(NRPAN,XRSF,DRSF,SCLM)
C
         DEALLOCATE (SIZEFAC,RQCLU_SFN,IQ_QCLU_SFN,WEIGHTCLU_SFN)
         DEALLOCATE (A3,B3,C3,D3,XEDGE,YEDGE,ZEDGE,NEDGE,RVERT_VFM,
     &               NEDGE_FCM)
         DEALLOCATE (NMSF,NM_S,MESHN_S,LMISF_S,XRN_S)
         DEALLOCATE (DRN_S,FLMSF_S,NFACE_M)
         DEALLOCATE (KLMFP_PR,NLMFPT_PR,NFPT_PR,LMIFP_PR,ISFLM_PR)
         DEALLOCATE (KLMSF_PR,NSF_PR,LMISF_PR,TXT_PR,LTXT_PR,IQ_MAUX)
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C
C=======================================================================
99001 FORMAT (/,10X,'volume of the atomic cells',//,10X,'IQ   IM',16X,
     &        'Volume',11X,'r_mt')
99002 FORMAT (I12,I5,7X,4F17.8)
99003 FORMAT (10X,'SUM',11X,F17.8,/,10X,'VUC',11X,F17.8,/)
99004 FORMAT (/,10X,'array sizes for the shape functions',/,10X,
     &        'NPANMAX =',I4,5X,'NRSFMAX =',I4,5X,'NSFMAX =',I4,//)
99005 FORMAT (/,' ##### WARNING from <SFNCREATECLU> ',48('#'),/,10X,A)
99006 FORMAT (/,' ##### TROUBLE in <SFNCREATECLU> ',50('#'),/,10X,A)
99007 FORMAT (10X,A,10I5)
99008 FORMAT (/,1X,79('*'),/)
99009 FORMAT (//,1X,79('*'),/,34X,'<SCFNCREATECLU>',/,1X,79('*'),//,10X,
     &        'create the shape functions for meshes IM =',I4,'  --',I3,
     &        /,10X,'write results to file ',A,//)
99010 FORMAT (/,' ***** INFO from <SFNCREATECLU> ',51('*'),/,10X,
     &        'for mesh IM=',I3,'   dummy shape function added',/)
      END
