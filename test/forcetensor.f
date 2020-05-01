C*==forcetensor.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FORCETENSOR(IECURR,P,TAUQ,MSSQ,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the site-off force tensor   PHI_ij                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NCL,IQ_MBCL
      USE MOD_SITES,ONLY:NQMAX,IQAT,ITOQ,NOQ,NQ,QBAS
      USE MOD_ANGMOM,ONLY:NL,NKMQ,NLM,NKMMAX,NKM,NXM,NSPIN
      USE MOD_CALCMODE,ONLY:KMROT,IREL
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,CONC,NT,LTXT_T
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:CI,C0,PI,RY_EV,EV_ERG,A0_CGS
      USE MOD_TAUIJ,ONLY:TAUIJ,DQCLU_TAUIJ_CL,NTAUIJ,ITAUIJ_LOOP,
     &    ITAUJI_LOOP,IQ_TAUIJ_CL,NQCLU_TAUIJ_CL,N5VEC_TAUIJ_CL,
     &    RQCLU_TAUIJ_CL,SELECTED_TAUIJ_CL
      IMPLICIT NONE
C*--FORCETENSOR21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FORCETENSOR')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 P
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      REAL*8 AUTOCGS,DITJT(:,:,:),PHI_IJ(3,3),PHI_ITJT(:,:,:,:,:),
     &       PREFAC,XMAX,XMIN,YMAX,YMIN
      COMPLEX*16 CFAC,CSUM,DELMSST(:,:,:,:,:),DMAMC(:,:),DMATTS(:,:,:,:)
     &           ,DTILTS(:,:,:,:),PHI_IJINT(:,:,:,:,:),UBAR(:,:,:),
     &           W1(:,:),W2(:,:),WBAR(:,:,:),WES,WSA(:,:),WSB(:,:)
      LOGICAL CHECK
      CHARACTER*80 FILNAM
      REAL*8 GAUNT_RYLM
      INTEGER I,IA_ERR,IC,ICL,ILOOP,IO,IQ,IS,IT,ITAUIJ,ITAUJI,IY,J,JC,
     &        JO,JQ,JQCLU,JT,L1,L2,L3,LFILNAM,LFN,LM1,LM2,M,M1,M2,M3,
     &        M_CART(3),N,NIJ,NITJT(:,:)
      CHARACTER*20 LEG(NT)
      SAVE PHI_IJINT,WBAR
C
C*** End of declarations rewritten by SPAG
C
      DATA IA_ERR/0/,M_CART/ + 1, - 1,0/,CHECK/.TRUE./
C
      ALLOCATABLE WBAR,UBAR
      ALLOCATABLE W1,W2,PHI_IJINT,WSA,WSB
      ALLOCATABLE DMAMC,DELMSST,DTILTS,DMATTS,PHI_ITJT,DITJT,NITJT
C
      IF ( NEPATH.NE.1 ) STOP 'in <FORCETENSOR>:  NEPATH <> 0'
      IF ( IGRID(1).NE.5 ) STOP 'in <FORCETENSOR>:  IGRID <> 5'
      IF ( IREL.EQ.2 .AND. KMROT.NE.0 )
     &      STOP 'in <FORCETENSOR>:  KMROT <> 0 for IREL = 2'
C
C-----------------------------------------------------------------------
C                           initialize
C-----------------------------------------------------------------------
C
      IF ( IREL.EQ.1 ) THEN
C
         WES = WETAB(IECURR,1)*2D0
C
      ELSE
C
         WES = WETAB(IECURR,1)
C
      END IF
C
      ALLOCATE (WSA(NXM,NXM),WSB(NXM,NXM))
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FORCETENSOR->W1'
C
      ALLOCATE (DMAMC(NKMMAX,NKMMAX),UBAR(NKMMAX,NKMMAX,3))
      ALLOCATE (DELMSST(NXM,NXM,NT,NSPIN,3))
      ALLOCATE (DTILTS(NXM,NXM,NT,NSPIN))
      ALLOCATE (DMATTS(NXM,NXM,NT,NSPIN),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:FORCETENSOR->DMATT'
C
      AUTOCGS = RY_EV*EV_ERG/(A0_CGS*A0_CGS)
C
C***********************************************************************
C             initialize PHI_IJINT and WBAR if  IECURR = 1
C***********************************************************************
C
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (PHI_IJINT(NTAUIJ,NT,NT,3,3),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc: FORCETENSOR -> PHI_IJINT'
C
         CALL CINIT(NTAUIJ*NT*NT*3*3,PHI_IJINT)
C
C-----------------------------------------------------------------------
C               auxilary transformation matrix  WBAR
C-----------------------------------------------------------------------
         ALLOCATE (WBAR(NKM,NKM,3),STAT=IA_ERR)
C
         CALL CINIT(NKM*NKM*3,WBAR)
C
C-------------------------------------- initiate WBAR for 2nd spin index
C
         PREFAC = SQRT(4D0*PI/3D0)
C
         LM1 = 0
         DO L1 = 0,(NL-1)
            DO M1 = -L1,L1
               LM1 = LM1 + 1
C
               LM2 = 0
               DO L2 = 0,(NL-1)
                  CFAC = PREFAC*CI**(L1-L2+1)
                  DO M2 = -L2,L2
                     LM2 = LM2 + 1
C
                     L3 = 1
                     DO IC = 1,3
                        M3 = M_CART(IC)
                        WBAR(LM1,LM2,IC)
     &                     = CFAC*GAUNT_RYLM(L1,M1,L2,M2,L3,M3)
                     END DO
C
                  END DO
               END DO
            END DO
         END DO
C
C-------------------------------------- initiate WBAR for 2nd spin index
C
         IF ( IREL.GE.2 ) WBAR(NLM+1:NKM,NLM+1:NKM,1:3)
     &        = WBAR(1:NLM,1:NLM,1:3)
C
C--------------------------------- convert to (kappa,mue) representation
         IF ( IREL.EQ.3 ) THEN
            DO IC = 1,3
C
               CALL CHANGEREP(NKM,NKMMAX,WBAR(1,1,IC),'RLM>REL',W1)
C
               WBAR(1:NKM,1:NKM,IC) = W1(1:NKM,1:NKM)
C
            END DO
         END IF
C
C.......................................................................
         IF ( CHECK ) THEN
            DO IC = 1,3
               CALL CMATSTRUCT('WBAR  '//(CHAR(ICHAR('1')+IC-1)),
     &                         WBAR(1,1,IC),NKM,NKM,IREL,IREL,0,1D-8,6)
            END DO
         END IF
C.......................................................................
C
      END IF
C
C***********************************************************************
C
C-----------------------------------------------------------------------
C                  reduced transformation matrix  UBAR
C-----------------------------------------------------------------------
      DO IC = 1,3
         DO J = 1,NKM
            DO I = 1,NKM
               UBAR(I,J,IC) = P*WBAR(I,J,IC)
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                      DELMSST  =  U * m - m U
C-----------------------------------------------------------------------
C
      DO IT = 1,NT
C
         DO IC = 1,3
C
            CALL CMATMUL(NKM,NKMMAX,UBAR(1,1,IC),MSST(1,1,IT),W1)
            CALL CMATMUL(NKM,NKMMAX,MSST(1,1,IT),UBAR(1,1,IC),W2)
C
            IF ( IREL.NE.2 ) THEN
C
               DELMSST(1:NKM,1:NKM,IT,1,IC) = W1(1:NKM,1:NKM)
     &            - W2(1:NKM,1:NKM)
C
            ELSE
C
               DELMSST(1:NLM,1:NLM,IT,1,IC) = W1(1:NLM,1:NLM)
     &            - W2(1:NLM,1:NLM)
C
               DELMSST(1:NLM,1:NLM,IT,2,IC) = W1(NLM+1:NKM,NLM+1:NKM)
     &            - W2(NLM+1:NKM,NLM+1:NKM)
C
            END IF
         END DO
C
C.......................................................................
         IF ( CHECK ) THEN
            IF ( IREL.LE.2 ) THEN
               IY = 1
            ELSE
               IY = 3
            END IF
            DO IC = 1,3
               DO IS = 1,NSPIN
                  CALL CMATSTRUCT('DELM '//(CHAR(ICHAR('1')+IS-1))
     &                            //' '//(CHAR(ICHAR('1')+IC-1)),
     &                            DELMSST(1,1,IT,IS,IC),NXM,NXM,IY,IY,0,
     &                            1D-8,6)
               END DO
            END DO
         END IF
C.......................................................................
C
C--------------------------- calculate projection matrices DMAT and DTIL
C
         IQ = IQAT(1,IT)
         N = NKMQ(IQ)
         M = NKMMAX
C
         CALL GETDMAT(TAUQ(1,1,IQ),W1,W2,DMAMC,N,MSSQ(1,1,IQ),
     &                MSST(1,1,IT),M)
C
         IF ( IREL.NE.2 ) THEN
C
            DMATTS(1:NXM,1:NXM,IT,1) = W1(1:NXM,1:NXM)
            DTILTS(1:NXM,1:NXM,IT,1) = W2(1:NXM,1:NXM)
C
         ELSE
C
            DMATTS(1:NLM,1:NLM,IT,1) = W1(1:NLM,1:NLM)
            DMATTS(1:NLM,1:NLM,IT,2) = W1(NLM+1:NKM,NLM+1:NKM)
C
            DTILTS(1:NLM,1:NLM,IT,1) = W2(1:NLM,1:NLM)
            DTILTS(1:NLM,1:NLM,IT,2) = W2(NLM+1:NKM,NLM+1:NKM)
C
         END IF
      END DO
C.......................................................................
      IF ( CHECK ) THEN
         IF ( IREL.LE.2 ) THEN
            IY = 1
         ELSE
            IY = 3
         END IF
         WRITE (*,*) '####',IREL,NSPIN
         IT = 1
         DO IS = 1,NSPIN
            CALL CMATSTRUCT('DMAT   '//(CHAR(ICHAR('1')+IS-1)),
     &                      DMATTS(1,1,IT,IS),NXM,NXM,IY,IY,0,1D-8,6)
         END DO
         DO IS = 1,NSPIN
            CALL CMATSTRUCT('DTIL   '//(CHAR(ICHAR('1')+IS-1)),
     &                      DTILTS(1,1,IT,IS),NXM,NXM,IY,IY,0,1D-8,6)
         END DO
      END IF
C.......................................................................
C
C=======================================================================
C                calculate site-off force tensor   PHI_ij
C       in case of alloy system: perform projection on atom types
C=======================================================================
C
      N = NXM
      M = NXM
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               ITAUJI = ITAUJI_LOOP(ILOOP)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
C
                  DO JO = 1,NOQ(JQ)
                     JT = ITOQ(JO,JQ)
C
                     DO IC = 1,3
                        DO JC = 1,3
                           DO IS = 1,NSPIN
C
C.......................................................................
                              IF ( CHECK .AND. ITAUIJ.LT.10 .AND. 
     &                             IC.EQ.1 .AND. JC.EQ.1 .AND. 
     &                             IS.EQ.NSPIN ) THEN
                                 IF ( IREL.LE.2 ) THEN
                                    IY = 1
                                 ELSE
                                    IY = 3
                                 END IF
                                 CALL CMATSTRUCT
     &                              ('TAUIJ '//(CHAR(ICHAR('1')
     &                              +ITAUIJ-1)),TAUIJ(1,1,1,1),NXM,NXM,
     &                              IY,IY,0,1D-8,6)
                              END IF
C.......................................................................
C
                              CALL CMATMUL(N,M,TAUIJ(1,1,IS,ITAUJI),
     &                           DTILTS(1,1,IT,IS),WSA)
C
                              CALL CMATMUL(N,M,DMATTS(1,1,JT,IS),WSA,
     &                           WSB)
C
                              CALL CMATMUL(N,M,DELMSST(1,1,JT,IS,JC),
     &                           WSB,WSA)
C
                              CALL CMATMUL(N,M,DTILTS(1,1,JT,IS),WSA,
     &                           WSB)
C
                              CALL CMATMUL(N,M,TAUIJ(1,1,IS,ITAUIJ),WSB,
     &                           WSA)
C
                              CALL CMATMUL(N,M,DMATTS(1,1,IT,IS),WSA,
     &                           WSB)
C
                              CALL CMATMUL(N,M,DELMSST(1,1,IT,IS,IC),
     &                           WSB,WSA)
C
                              CSUM = C0
                              DO I = 1,NXM
                                 CSUM = CSUM + WSA(I,I)
                              END DO
C
                              PHI_IJINT(ITAUIJ,IT,JT,IC,JC)
     &                           = PHI_IJINT(ITAUIJ,IT,JT,IC,JC)
     &                           + WES*CSUM
C
                           END DO
                        END DO
                     END DO
C
                  END DO
               END DO
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      IF ( IECURR.NE.NETAB(1) ) THEN
         DEALLOCATE (TAUIJ,W1,W2,DMAMC,DTILTS,DMATTS)
         RETURN
      END IF
C
C-----------------------------------------------------------------------
C                              write results
C-----------------------------------------------------------------------
C
      FILNAM = DATSET(1:LDATSET)//'_F_ij.dat'
      LFN = LDATSET + 9
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
C
      WRITE (*,99003)
      WRITE (IOTMP,99003)
C
      WRITE (IOTMP,99005)
      DO IQ = 1,NQ
         WRITE (IOTMP,99006) IQ,NOQ(IQ),
     &                       (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,NOQ(IQ)
     &                       )
      END DO
C
      NIJ = SUM(NQCLU_TAUIJ_CL(1:NCL))
C
      ALLOCATE (PHI_ITJT(NIJ,NT,NT,3,3),DITJT(NIJ,NT,NT),STAT=IA_ERR)
      ALLOCATE (NITJT(NT,NT),STAT=IA_ERR)
C
      NITJT(1:NT,1:NT) = 0
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
C-----------------------------------------------------------------------
            IF ( SELECTED_TAUIJ_CL(JQCLU,ICL) ) THEN
C
               JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
               ITAUJI = ITAUJI_LOOP(ILOOP)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
C
                  DO JO = 1,NOQ(JQ)
                     JT = ITOQ(JO,JQ)
C
                     WRITE (*,99001) IQ,IT,JQ,JT,(QBAS(I,IQ),I=1,3),
     &                               (QBAS(I,JQ),I=1,3)
C
                     WRITE (IOTMP,99001) IQ,IT,JQ,JT,(QBAS(I,IQ),I=1,3),
     &                      (QBAS(I,JQ),I=1,3)
C
                     NITJT(IT,JT) = NITJT(IT,JT) + 1
                     I = NITJT(IT,JT)
                     DO JC = 1,3
                        DO IC = 1,3
C
                           PHI_IJ(IC,JC)
     &                        = DIMAG(PHI_IJINT(ITAUIJ,IT,JT,IC,JC))/PI
                           PHI_ITJT(I,IT,JT,IC,JC) = PHI_IJ(IC,JC)
     &                        *AUTOCGS*1D-3
                           DITJT(I,IT,JT) = DQCLU_TAUIJ_CL(JQCLU,ICL)
                        END DO
                     END DO
C
                     WRITE (6,99002) ITAUIJ,ITAUJI,
     &                               N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                               RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                               DQCLU_TAUIJ_CL(JQCLU,ICL),
     &                               ((PHI_IJ(IC,JC),PHI_IJ(IC,JC)
     &                               *AUTOCGS,JC=1,3),IC=1,3)
C
                     WRITE (IOTMP,99002) ITAUIJ,ITAUJI,
     &                      N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                      RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                      DQCLU_TAUIJ_CL(JQCLU,ICL),
     &                      ((PHI_IJ(IC,JC)*AUTOCGS,JC=1,3),IC=1,3)
C
                  END DO
               END DO
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
      WRITE (*,99004) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C
      XMIN = 1D10
      XMAX = 0D0
      DO IT = 1,NT
         LEG(IT) = TXT_T(IT)
         DO JT = 1,NT
            DO I = 1,NITJT(IT,JT)
               XMIN = MIN(XMIN,DITJT(I,IT,JT))
               XMAX = MAX(XMAX,DITJT(I,IT,JT))
            END DO
         END DO
      END DO
C
      DO IT = 1,NT
C
         YMIN = 0D0
         YMAX = 0D0
C
         DO JC = 1,3
            DO IC = 1,3
               DO JT = 1,NT
                  DO I = 1,NITJT(IT,JT)
                     YMIN = MIN(YMIN,PHI_ITJT(I,IT,JT,IC,JC))
                     YMAX = MAX(YMAX,PHI_ITJT(I,IT,JT,IC,JC))
                  END DO
               END DO
            END DO
         END DO
C
         CALL XMGRHEAD(DATSET,LDATSET,'Fij',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,'R!sij!N/a',9,
     &                 '!xF!0!sij!N (10!S3!N erg/cm!S2!N)',33,' ',0,
     &                 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &                 25+LSYSTEM,
     &                 'force tensor elements    !xF!0!sij!N for '//
     &                 TXT_T(IT)(1:LTXT_T(IT))//' at center',
     &                 (41+LTXT_T(IT)+10),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,1,NT,0,LEG,LEG)
C
         CALL XMGRPOINTS(IOTMP,1,NT,0,2,1,1)
C
         IY = -1
         DO JT = 1,NT
            DO JC = 1,3
               DO IC = 1,3
C
                  IY = IY + 1
                  CALL XMGRTABLE(0,IY,DITJT(1,IT,JT),
     &                           PHI_ITJT(1,IT,JT,IC,JC),1.0D0,
     &                           NITJT(IT,JT),IOTMP)
C
               END DO
            END DO
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '   force tensor elements written to file ',
     &               FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IOTMP)
C
      END DO
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (PHI_IJINT,TAUIJ,W1,W2)
      DEALLOCATE (DMAMC,DTILTS,DMATTS,PHI_ITJT,DITJT,NITJT)
C
C-----------------------------------------------------------------------
99001 FORMAT (/,5X,'IQ =',I3,' IT =',I3,20X,' JQ =',I3,' JT =',I3,/,5X,
     &        2(3X,'->Q = (',F7.3,',',F7.3,',',F7.3,')',5X),/,6X,
     &        'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR ',
     &        '    PHI_ij [Ry]  PHI_ij [eV]')
99002 FORMAT (5X,2I7,2X,3I3,1X,3F7.3,1X,F7.3,1X,6F14.6,2(/,61X,6F14.6))
99003 FORMAT (/,1X,79('*'),/,35X,'<FORCETENSOR>',/,28X,
     &        'Force tensor       PHI_ij',/,1X,79('*'),/)
99004 FORMAT (/,10X,'results written to file:',A,/)
99005 FORMAT (/,10X,'number of sites   NQ = ',I3,/,10X,
     &        'number of types   NT = ',I3,/,10X,'site occupation:')
99006 FORMAT (10X,2I4,10(I3,F6.3))
      END
