C*==cluxcpltensor.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUXCPLTENSOR(IECURR,TAUQ,MSSQ,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the site-off diagonal exchange tensor  J_ij           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQMAX,IQAT,IT_OQCLU,ITCLU_OQCLU,NQCLU,NQHOST,
     &    QBAS_QCLU,NO_QCLU
      USE MOD_ANGMOM,ONLY:NKMQ,NKMMAX,NKM
      USE MOD_CALCMODE,ONLY:NONMAG,KMROT,IREL
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,X_TCLU,NT,LTXT_T,Z,NTCLU,NTHOST
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      USE MOD_TAUIJ,ONLY:TAUQQ_CLU
      USE MOD_MPI,ONLY:MPI_ELOOP,MPI_KLOOP,MPI_ID
      IMPLICIT NONE
C*--CLUXCPLTENSOR20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUXCPLTENSOR')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,CWK(:),DELMSST(:,:,:,:),DMAMC(:,:),DMATT(:,:,:),
     &           DTILT(:,:,:),J_IJINT(:,:,:,:,:,:),WES,WSA(:,:),WSB(:,:)
      REAL*8 DIJ,DITJT(:,:,:),J_IJ(3,3),J_ITJT(:,:,:,:,:),RIJVEC(3),
     &       XMAX,XMIN,YMAX,YMIN
      REAL*8 DNRM2
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IC,IO,IQ,IQCLU,IT,ITCLU,IY,JC,JO,JQCLU,JT,JTCLU,
     &        LFILNAM,LFN,LWK,M,N,NIJ,NITJT(:,:)
      LOGICAL INITIALIZE
      CHARACTER*20 LEG(NT)
      CHARACTER*1 TXT_DM_PROJ(3)
      SAVE J_IJINT
C
C*** End of declarations rewritten by SPAG
C
      DATA TXT_DM_PROJ/'x','y','z'/,INITIALIZE/.TRUE./
C
      ALLOCATABLE J_IJINT,WSA,WSB,CWK
      ALLOCATABLE DMAMC,DELMSST,DTILT,DMATT,J_ITJT,DITJT,NITJT
C
      IF ( NONMAG ) CALL STOP_MESSAGE(ROUTINE,'non-magnetic system')
      IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KMROT <> 0 ')
      IF ( NEPATH.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NEPATH <> 1')
      IF ( IGRID(1).NE.5 ) CALL STOP_MESSAGE(ROUTINE,'IGRID <> 5')
      IF ( IREL.LE.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL < 3')
      IF ( NKM.NE.NKMMAX ) CALL STOP_MESSAGE(ROUTINE,'NKM <> NKMMAX ')
C
C-----------------------------------------------------------------------
C                           initialize
C-----------------------------------------------------------------------
C
      WES = WETAB(IECURR,1)
C
      ALLOCATE (WSA(NKMMAX,NKMMAX),WSB(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: W1')
C
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (DELMSST(NKMMAX,NKMMAX,3,NT))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NT))
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DMATT')
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (J_IJINT(NQCLU,NQCLU,NT,NT,3,3),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: J_IJINT')
C
         J_IJINT(:,:,:,:,:,:) = 0D0
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ELOOP .AND. IECURR.GT.NETAB(1) ) THEN
C
C-----------------------------------  collect results from all processes
C
         CALL DRV_MPI_BARRIER
C
         LWK = NQCLU*NQCLU*NT*NT*3*3
         ALLOCATE (CWK(LWK))
C
         CALL DRV_MPI_REDUCE_C(J_IJINT(1,1,1,1,1,1),CWK(1),LWK)
C
         IF ( MPI_ID.EQ.0 ) GOTO 100
         RETURN
C
      END IF
C
      IF ( MPI_KLOOP .AND. MPI_ID.NE.0 ) RETURN
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C***********************************************************************
C             calculate the matrix elements  DELMSST
C***********************************************************************
C
      CALL XCPLTENME(DELMSST)
C
C***********************************************************************
C             calculate projection matrices DMAT and DTIL
C***********************************************************************
C
      DO IT = 1,NT
C
         IQ = IQAT(1,IT)
C
         IF ( Z(IT).EQ.0 ) CYCLE
         IF ( IQ.LE.NQHOST ) CYCLE
C
         N = NKMQ(IQ)
         M = NKMMAX
C
         CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT),DMAMC,N,
     &                MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
      END DO
C
C=======================================================================
C        calculate elements of site-off exchange tensor  J_ij
C       in case of alloy system: perform projection on atom types
C=======================================================================
C
      N = NKM
      M = NKMMAX
C
      DO IQCLU = 1,NQCLU
C
         IQ = NQHOST + IQCLU
C
         DO JQCLU = 1,NQCLU
C
            IF ( IQCLU.EQ.JQCLU ) CYCLE
C
C-----------------------------------------------------------------------
C
            DO IO = 1,NO_QCLU(IQCLU)
               IT = IT_OQCLU(IO,IQCLU)
               ITCLU = ITCLU_OQCLU(IO,IQCLU)
               IF ( Z(IT).EQ.0 ) CYCLE
C
               DO JO = 1,NO_QCLU(JQCLU)
                  JT = IT_OQCLU(JO,JQCLU)
                  JTCLU = ITCLU_OQCLU(JO,JQCLU)
                  IF ( Z(JT).EQ.0 ) CYCLE
C
                  DO IC = 1,3
                     DO JC = 1,3
C
                        CALL CMATMUL(N,M,TAUQQ_CLU(1,1,JQCLU,IQCLU),
     &                               DTILT(1,1,IT),WSA)
C
                        CALL CMATMUL(N,M,DMATT(1,1,JT),WSA,WSB)
C
                        CALL CMATMUL(N,M,DELMSST(1,1,JC,JT),WSB,WSA)
C
                        CALL CMATMUL(N,M,DTILT(1,1,JT),WSA,WSB)
C
                        CALL CMATMUL(N,M,TAUQQ_CLU(1,1,IQCLU,JQCLU),WSB,
     &                               WSA)
C
                        CALL CMATMUL(N,M,DMATT(1,1,IT),WSA,WSB)
C
                        CALL CMATMUL(N,M,DELMSST(1,1,IC,IT),WSB,WSA)
C
                        CSUM = C0
                        DO I = 1,NKM
                           CSUM = CSUM + WSA(I,I)
                        END DO
C
                        J_IJINT(IQCLU,JQCLU,ITCLU,JTCLU,IC,JC)
     &                     = J_IJINT(IQCLU,JQCLU,ITCLU,JTCLU,IC,JC)
     &                     + WES*CSUM*0.5D0
C
C-----------------------------------------------------------------------
C              NOTE:      the factor 0.5 here - to keep the definition
C                         of Heisenberg Hamitonian in the same form
C                         as used in non-relativistic J_ij calculations
C-----------------------------------------------------------------------
C
                     END DO
                  END DO
C
               END DO
            END DO
C-----------------------------------------------------------------------
         END DO
      END DO
C-----------------------------------------------------------------------
C
      IF ( IECURR.NE.NETAB(1) ) RETURN
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C process MPI_ID = 0 continues here in case of MPI and IECURR > NETAB(1)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C-----------------------------------------------------------------------
C                 write results for isotropic exchange
C-----------------------------------------------------------------------
C
 100  CONTINUE
      FILNAM = DATSET(1:LDATSET)//'_JJij.dat'
      LFN = LDATSET + 9
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
      WRITE (6,99008)
      WRITE (6,99001)
C
      WRITE (IOTMP,99008)
C
      WRITE (IOTMP,99006) NQCLU,NTCLU
      DO IQCLU = 1,NQCLU
         WRITE (IOTMP,99007) IQCLU,NO_QCLU(IQCLU),
     &                       (ITCLU_OQCLU(IO,IQCLU),
     &                       X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                       NO_QCLU(IQCLU))
      END DO
C
      WRITE (IOTMP,99001)
C
C-------------------------------------------------to be on the safe side
      NIJ = MAX(NTCLU*NTCLU,NQCLU*NQCLU)
C
      ALLOCATE (J_ITJT(NIJ,NTCLU,NTCLU,3,3))
      ALLOCATE (DITJT(NIJ,NTCLU,NTCLU))
      ALLOCATE (NITJT(NTCLU,NTCLU),STAT=IA_ERR)
C
      NITJT(1:NTCLU,1:NTCLU) = 0
C
      DO IQCLU = 1,NQCLU
C
         IQ = NQHOST + IQCLU
C
         DO JQCLU = 1,NQCLU
C
            IF ( IQCLU.EQ.JQCLU ) CYCLE
C
            RIJVEC(:) = QBAS_QCLU(1:3,IQCLU) - QBAS_QCLU(1:3,JQCLU)
            DIJ = DNRM2(3,RIJVEC,1)
C
C-----------------------------------------------------------------------
C
            DO IO = 1,NO_QCLU(IQCLU)
               IT = IT_OQCLU(IO,IQCLU)
               ITCLU = ITCLU_OQCLU(IO,IQCLU)
               IF ( Z(IT).EQ.0 ) CYCLE
C
               DO JO = 1,NO_QCLU(JQCLU)
                  JT = IT_OQCLU(JO,JQCLU)
                  JTCLU = ITCLU_OQCLU(JO,JQCLU)
                  IF ( Z(JT).EQ.0 ) CYCLE
C
C-------------------------------------------------- exclude on-site term
C
                  WRITE (6,99002) IQCLU,ITCLU,JQCLU,JTCLU,
     &                            (QBAS_QCLU(I,IQCLU),I=1,3),
     &                            (QBAS_QCLU(I,JQCLU),I=1,3)
C
                  NITJT(ITCLU,JTCLU) = NITJT(ITCLU,JTCLU) + 1
                  I = NITJT(ITCLU,JTCLU)
C
                  DO JC = 1,3
                     DO IC = 1,3
C
                        J_IJ(IC,JC) = DIMAG(J_IJINT(IQCLU,JQCLU,ITCLU,
     &                                JTCLU,IC,JC))/PI
                        J_ITJT(I,ITCLU,JTCLU,IC,JC) = J_IJ(IC,JC)
     &                     *RY_EV*1.D3
                        DITJT(I,ITCLU,JTCLU) = DIJ
                     END DO
                  END DO
C
                  IF ( (I.EQ.1) .AND. (IQCLU.EQ.JQCLU) )
     &                 J_ITJT(1,ITCLU,JTCLU,1,1) = 0.D0
C
                  WRITE (6,99003) ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                            RIJVEC(1:3),DIJ,
     &                            J_ITJT(I,ITCLU,JTCLU,1,1)/RY_EV,
     &                            J_ITJT(I,ITCLU,JTCLU,1,1)
C
                  WRITE (IOTMP,99003) ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                                RIJVEC(1:3),DIJ,
     &                                J_ITJT(I,ITCLU,JTCLU,1,1)/RY_EV,
     &                                J_ITJT(I,ITCLU,JTCLU,1,1)
C
               END DO
            END DO
C
C-----------------------------------------------------------------------
         END DO
      END DO
C
      WRITE (6,99005) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C                      find suitable plot windows
C-----------------------------------------------------------------------
C
      XMIN = 1D10
      XMAX = 0D0
      DO ITCLU = 1,NTCLU
         LEG(ITCLU) = TXT_T(ITCLU)
         DO JTCLU = 1,NTCLU
            DO I = 1,NITJT(ITCLU,JTCLU)
               XMIN = MIN(XMIN,DITJT(I,ITCLU,JTCLU))
               XMAX = MAX(XMAX,DITJT(I,ITCLU,JTCLU))
            END DO
         END DO
      END DO
C
C-------------------------------------------- start plot at X = R_ij = 0
      XMIN = 0D0
C
      DO ITCLU = 1,NTCLU
C
         IT = NTHOST + ITCLU
C
         YMIN = 0D0
         YMAX = 0D0
C
         DO JTCLU = 1,NTCLU
            DO I = 1,NITJT(ITCLU,JTCLU)
               YMIN = MIN(YMIN,J_ITJT(I,ITCLU,JTCLU,1,1))
               YMAX = MAX(YMAX,J_ITJT(I,ITCLU,JTCLU,1,1))
            END DO
         END DO
C
         CALL XMGRHEAD(DATSET,LDATSET,'Jij',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,'R!sij!N/a',9,
     &                 'J!0!sij!N (meV)',15,' ',0,
     &                 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &                 25+LSYSTEM,
     &                 'Isotropic exchange coupling   J!0!sij!N for '//
     &                 TXT_T(IT)(1:LTXT_T(IT))//' at center',
     &                 (42+LTXT_T(IT)+10),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,1,NTCLU,0,LEG,LEG)
C
         CALL XMGRPOINTS(IOTMP,1,NTCLU,0,2,1,1)
C
         IY = -1
         DO JTCLU = 1,NTCLU
            IY = IY + 1
            CALL XMGRTABLE(0,IY,DITJT(1,ITCLU,JTCLU),
     &                     J_ITJT(1,ITCLU,JTCLU,1,1),1.0D0,
     &                     NITJT(ITCLU,JTCLU),IOTMP)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '   exchange couplings are written to file ',
     &               FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IOTMP)
C
      END DO
C
C-----------------------------------------------------------------------
C         write results for Dzyaloshinski-Moriya interaction
C-----------------------------------------------------------------------
C
      FILNAM = DATSET(1:LDATSET)//'_Dij_x.dat'
      LFN = LDATSET + 10
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP+1,FILNAM(1:LFN))
      FILNAM = DATSET(1:LDATSET)//'_Dij_y.dat'
      LFN = LDATSET + 10
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP+2,FILNAM(1:LFN))
      FILNAM = DATSET(1:LDATSET)//'_Dij_z.dat'
      LFN = LDATSET + 10
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP+3,FILNAM(1:LFN))
C
      WRITE (6,99010)
      WRITE (IOTMP+1,99004) TXT_DM_PROJ(1)
      WRITE (IOTMP+1,99006) NQCLU,NTCLU
      DO IQCLU = 1,NQCLU
         WRITE (IOTMP+1,99007) IQCLU,NO_QCLU(IQCLU),
     &                         (ITCLU_OQCLU(IO,IQCLU),
     &                         X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                         NO_QCLU(IQCLU))
      END DO
C
      WRITE (IOTMP+2,99004) TXT_DM_PROJ(2)
      WRITE (IOTMP+2,99006) NQCLU,NTCLU
      DO IQCLU = 1,NQCLU
         WRITE (IOTMP+2,99007) IQCLU,NO_QCLU(IQCLU),
     &                         (ITCLU_OQCLU(IO,IQCLU),
     &                         X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                         NO_QCLU(IQCLU))
      END DO
C
      WRITE (IOTMP+3,99004) TXT_DM_PROJ(3)
      WRITE (IOTMP+3,99006) NQCLU,NTCLU
      DO IQCLU = 1,NQCLU
         WRITE (IOTMP+3,99007) IQCLU,NO_QCLU(IQCLU),
     &                         (ITCLU_OQCLU(IO,IQCLU),
     &                         X_TCLU(ITCLU_OQCLU(IO,IQCLU)),IO=1,
     &                         NO_QCLU(IQCLU))
      END DO
C
C-------------------------------------------------to be on the safe side
      NIJ = MAX(NTCLU*NTCLU,NQCLU*NQCLU)
      NITJT(1:NTCLU,1:NTCLU) = 0
C
      WRITE (6,99001)
      WRITE (IOTMP+1,99001)
      WRITE (IOTMP+2,99001)
      WRITE (IOTMP+3,99001)
C
      DO IQCLU = 1,NQCLU
C
         IQ = NQHOST + IQCLU
C
         DO JQCLU = 1,NQCLU
C
            IF ( IQCLU.EQ.JQCLU ) CYCLE
C
            RIJVEC(:) = QBAS_QCLU(1:3,IQCLU) - QBAS_QCLU(1:3,JQCLU)
            DIJ = DNRM2(3,RIJVEC,1)
C
C-----------------------------------------------------------------------
C
            DO IO = 1,NO_QCLU(IQCLU)
               IT = IT_OQCLU(IO,IQCLU)
               ITCLU = ITCLU_OQCLU(IO,IQCLU)
               IF ( Z(IT).EQ.0 ) CYCLE
C
               DO JO = 1,NO_QCLU(JQCLU)
                  JT = IT_OQCLU(JO,JQCLU)
                  JTCLU = ITCLU_OQCLU(JO,JQCLU)
                  IF ( Z(JT).EQ.0 ) CYCLE
C
                  WRITE (6,99002) IQCLU,ITCLU,JQCLU,JTCLU,
     &                            (QBAS_QCLU(I,IQCLU),I=1,3),
     &                            (QBAS_QCLU(I,JQCLU),I=1,3)
C
C
                  NITJT(ITCLU,JTCLU) = NITJT(ITCLU,JTCLU) + 1
                  I = NITJT(ITCLU,JTCLU)
C
C
                  J_ITJT(I,ITCLU,JTCLU,3,3)
     &               = (J_ITJT(I,ITCLU,JTCLU,1,2)-
     &               J_ITJT(I,ITCLU,JTCLU,2,1))*0.5D0
                  J_ITJT(I,ITCLU,JTCLU,2,2)
     &               = -(J_ITJT(I,ITCLU,JTCLU,1,3)
     &               -J_ITJT(I,ITCLU,JTCLU,3,1))*0.5D0
                  J_ITJT(I,ITCLU,JTCLU,1,1)
     &               = (J_ITJT(I,ITCLU,JTCLU,2,3)-
     &               J_ITJT(I,ITCLU,JTCLU,3,2))*0.5D0
C
                  WRITE (6,99009) 'x',ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                            RIJVEC(1:3),DIJ,
     &                            J_ITJT(I,ITCLU,JTCLU,1,1)/RY_EV,
     &                            J_ITJT(I,ITCLU,JTCLU,1,1)
C
                  WRITE (IOTMP+1,99003) ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                                  RIJVEC(1:3),DIJ,
     &                                  J_ITJT(I,ITCLU,JTCLU,1,1)/RY_EV,
     &                                  J_ITJT(I,ITCLU,JTCLU,1,1)
C
                  WRITE (6,99009) 'y',ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                            RIJVEC(1:3),DIJ,
     &                            J_ITJT(I,ITCLU,JTCLU,2,2)/RY_EV,
     &                            J_ITJT(I,ITCLU,JTCLU,2,2)
C
                  WRITE (IOTMP+2,99003) ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                                  RIJVEC(1:3),DIJ,
     &                                  J_ITJT(I,ITCLU,JTCLU,2,2)/RY_EV,
     &                                  J_ITJT(I,ITCLU,JTCLU,2,2)
C
                  WRITE (6,99009) 'z',ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                            RIJVEC(1:3),DIJ,
     &                            J_ITJT(I,ITCLU,JTCLU,3,3)/RY_EV,
     &                            J_ITJT(I,ITCLU,JTCLU,3,3)
C
                  WRITE (IOTMP+3,99003) ITCLU,IQCLU,JTCLU,JQCLU,0,0,0,
     &                                  RIJVEC(1:3),DIJ,
     &                                  J_ITJT(I,ITCLU,JTCLU,3,3)/RY_EV,
     &                                  J_ITJT(I,ITCLU,JTCLU,3,3)
C
               END DO
            END DO
C-----------------------------------------------------------------------
         END DO
      END DO
C
      WRITE (6,99005) DATSET(1:LDATSET)//'_Dij_x.dat'
      WRITE (6,99005) DATSET(1:LDATSET)//'_Dij_y.dat'
      WRITE (6,99005) DATSET(1:LDATSET)//'_Dij_z.dat'
C-----------------------------------------------------------------------
C
      XMIN = 1D10
      XMAX = 0D0
      DO ITCLU = 1,NTCLU
         LEG(ITCLU) = TXT_T(ITCLU)
         DO JTCLU = 1,NTCLU
            DO I = 1,NITJT(ITCLU,JTCLU)
               XMIN = MIN(XMIN,DITJT(I,ITCLU,JTCLU))
               XMAX = MAX(XMAX,DITJT(I,ITCLU,JTCLU))
            END DO
         END DO
      END DO
C
      DO ITCLU = 1,NTCLU
         IT = NTHOST + ITCLU
         DO IC = 1,3
C
            YMIN = 0D0
            YMAX = 0D0
C
            DO JTCLU = 1,NTCLU
               DO I = 1,NITJT(ITCLU,JTCLU)
                  YMIN = MIN(YMIN,J_ITJT(I,ITCLU,JTCLU,IC,IC))
                  YMAX = MAX(YMAX,J_ITJT(I,ITCLU,JTCLU,IC,IC))
               END DO
            END DO
C
            CALL XMGRHEAD(DATSET,LDATSET,'Dij_'//TXT_DM_PROJ(IC),5,
     &                    TXT_T(IT),LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,
     &                    1,XMIN,1,XMAX,1,YMIN,0,YMAX,1,YMIN,0,YMAX,0,
     &                    'R!sij!N/a',9,'D!0!S'//TXT_DM_PROJ(IC)
     &                    //'!N!sij!N (meV)',20,' ',0,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'Components of DM vector D!S'//
     &                    TXT_DM_PROJ(IC)//'!N!sij!N for '//TXT_T(IT)
     &                    (1:LTXT_T(IT))//' at center',
     &                    (28+13+LTXT_T(IT)+10),.FALSE.)
C
            CALL XMGRLEGEND(IOTMP,1,NTCLU,0,LEG,LEG)
C
            CALL XMGRPOINTS(IOTMP,1,NTCLU,0,2,1,1)
C
            IY = -1
            DO JTCLU = 1,NTCLU
               IY = IY + 1
               CALL XMGRTABLE(0,IY,DITJT(1,ITCLU,JTCLU),
     &                        J_ITJT(1,ITCLU,JTCLU,IC,IC),1.0D0,
     &                        NITJT(ITCLU,JTCLU),IOTMP)
            END DO
C
            WRITE (6,*) ' '
            WRITE (6,*) 
     &                '   Exchange tensor elements are written to file '
     &                ,FILNAM(1:LFILNAM)
            WRITE (6,*) ' '
            CLOSE (IOTMP)
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
99001 FORMAT (6X,
     &     '  ITCLU IQCLU JTCLU JQCLU DUMMY   DRX    DRY    DRZ     DR '
     &     ,'        J_ij [mRy]     J_ij [meV]')
99002 FORMAT (/,5X,'ITCLU =',I3,' IQCLU =',I3,20X,' JTCLU =',I3,
     &        ' JQCLU =',I3,/,5X,
     &        2(3X,'->Q = (',F7.3,',',F7.3,',',F7.3,')'))
99003 FORMAT (5X,4I5,2X,3I3,1X,3F7.3,1X,F7.3,1X,6F16.9)
99004 FORMAT (/,1X,79('*'),/,35X,'<XCPLTENSOR>',/,28X,
     &        'Dzyaloshinski-Moriya couplings  Dij_',A1/,1X,79('*'),/)
99005 FORMAT (/,10X,'results written to file:',A,/)
99006 FORMAT (/,10X,'number of sites   NQCLU = ',I3,/,10X,
     &        'number of types   NTCLU = ',I3,/,10X,'site occupation:')
99007 FORMAT (10X,2I4,10(I3,F6.3))
99008 FORMAT (/,1X,79('*'),/,35X,'<XCPLTENSOR>',/,28X,
     &        'Isotropic exchange couplings  Jij',A1/,1X,79('*'),/)
99009 FORMAT (3X,A1,3X,4I5,2X,3I3,1X,3F7.3,1X,F7.3,1X,6F16.9)
99010 FORMAT (/,1X,79('*'),/,35X,'<XCPLTENSOR>',/,10X,
     &        'Components of Dzyaloshinski-Moriya vector Dij',A1/,1X,
     &        79('*'),/)
      END
