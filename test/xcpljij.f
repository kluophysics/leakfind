C*==xcpljij.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XCPLJIJ(IECURR,TAUQ,MSSQ,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the site-off diagonal  XC-coupling parameters   J_ij  *
C   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
C   *                                                                  *
C   *  NOTE: the site off-diagonal scattering path operator  TAUIJ     *
C   *        is imported via the module  MOD_TAUIJ                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NCL,IQ_MBCL
      USE MOD_SITES,ONLY:NQMAX,IQAT,ITOQ,NOQ,NQ,QBAS,IQBOT,IQTOP
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,NLM
      USE MOD_CALCMODE,ONLY:NONMAG,KMROT,IREL
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,CONC,NT,LTXT_T,ITBOT,ITTOP
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET,
     &    IFILBUILDBOT,WRBUILDBOT
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV,EV_J,KB_SI
      USE MOD_TAUIJ,ONLY:TAUIJ,DQCLU_TAUIJ_CL,NTAUIJ,ITAUIJ_LOOP,
     &    ITAUJI_LOOP,IQ_TAUIJ_CL,NQCLU_TAUIJ_CL,N5VEC_TAUIJ_CL,
     &    RQCLU_TAUIJ_CL,SELECTED_TAUIJ_CL
      USE MOD_MPI,ONLY:MPI_ELOOP,MPI_KLOOP,MPI_ID
      IMPLICIT NONE
C*--XCPLJIJ27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='XCPLJIJ')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,CWK(:),DELMSST(NLM,NLM,NT),DMAMC(:,:),
     &           DMATTS(:,:,:,:),DTILTS(:,:,:,:),J_IJINT(:,:,:),
     &           MQS(:,:,:,:),W1(:,:),W2(:,:),WSA(NLM,NLM),WSB(NLM,NLM)
      REAL*8 DITJT(:,:,:),JXC0(NT),JXCIJ,JXCITJT(:,:,:),JXCITJT_0(NT,NT)
     &       ,TCURIE,T_C(NT,NT),W4(:),WI(NT),WR(NT),WW1(NT,NT),
     &       WW2(NT,NT),XMAX,XMIN,YMAX,YMIN
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,ICL,ILOOP,INFO,IO,IOP,IQ,IT,ITAUIJ,ITAUJI,J,JO,
     &        JQ,JQCLU,JT,LFILNAM,LFN,LWK,LWMAX,M,N,NIJ,NITJT(:,:),NTC
      LOGICAL INITIALIZE
      CHARACTER*20 LEG(NT)
      SAVE J_IJINT
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE W1,W2,W4,MQS,J_IJINT,CWK
      ALLOCATABLE DMAMC,DTILTS,DMATTS,JXCITJT,DITJT,NITJT
C
      IF ( NONMAG ) CALL STOP_MESSAGE(ROUTINE,'non-magnetic system')
      IF ( IREL.NE.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL <> 2 ')
      IF ( KMROT.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'KMROT <> 0 ')
      IF ( NEPATH.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NEPATH <> 0')
      IF ( IGRID(1).NE.5 ) CALL STOP_MESSAGE(ROUTINE,'IGRID <> 5')
C
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX))
      ALLOCATE (DMAMC(NKMMAX,NKMMAX),MQS(NLM,NLM,NQ,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      ALLOCATE (DMATTS(NLM,NLM,NT,2),DTILTS(NLM,NLM,NT,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DMATT')
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (J_IJINT(NTAUIJ,NT,NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: J_IJINT')
C
         J_IJINT(:,:,:) = 0D0
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
         LWK = NTAUIJ*NT*NT
         ALLOCATE (CWK(LWK))
C
         CALL DRV_MPI_REDUCE_C(J_IJINT(1,1,1),CWK(1),LWK)
C
         DEALLOCATE (CWK)
         CALL DRV_MPI_BARRIER
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
      DO IQ = IQBOT,IQTOP
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,MSSQ(1,J,IQ),1,MQS(1,J,IQ,1),1)
            CALL ZCOPY(NLM,MSSQ(NLM+1,NLM+J,IQ),1,MQS(1,J,IQ,2),1)
         END DO
C
      END DO
C
      DO IT = ITBOT,ITTOP
C
         DO J = 1,NLM
            DO I = 1,NLM
               DELMSST(I,J,IT) = MSST(NLM+I,NLM+J,IT) - MSST(I,J,IT)
            END DO
         END DO
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
         DO J = 1,NLM
            CALL ZCOPY(NLM,W1(1,J),1,DMATTS(1,J,IT,1),1)
            CALL ZCOPY(NLM,W1(NLM+1,NLM+J),1,DMATTS(1,J,IT,2),1)
         END DO
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,W2(1,J),1,DTILTS(1,J,IT,1),1)
            CALL ZCOPY(NLM,W2(NLM+1,NLM+J),1,DTILTS(1,J,IT,2),1)
         END DO
C
      END DO
C
C=======================================================================
C                   calculate J_ij via Eq. (19)
C       in case of alloy system: perform projection on atom types
C=======================================================================
C
      N = NLM
      M = NLM
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
                     CALL CMATMUL(N,M,TAUIJ(1,1,1,ITAUJI),
     &                            DTILTS(1,1,IT,1),WSA)
C
                     CALL CMATMUL(N,M,DMATTS(1,1,JT,1),WSA,WSB)
C
                     CALL CMATMUL(N,M,DELMSST(1,1,JT),WSB,WSA)
C
                     CALL CMATMUL(N,M,DTILTS(1,1,JT,2),WSA,WSB)
C
                     CALL CMATMUL(N,M,TAUIJ(1,1,2,ITAUIJ),WSB,WSA)
C
                     CALL CMATMUL(N,M,DMATTS(1,1,IT,2),WSA,WSB)
C
                     CALL CMATMUL(N,M,DELMSST(1,1,IT),WSB,WSA)
C
                     CSUM = C0
                     DO I = 1,NLM
                        CSUM = CSUM + WSA(I,I)
                     END DO
C
                     J_IJINT(ITAUIJ,IT,JT) = J_IJINT(ITAUIJ,IT,JT)
     &                  + CSUM*WETAB(IECURR,1)
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
         DEALLOCATE (TAUIJ,W1,W2,MQS)
         DEALLOCATE (DMAMC,DTILTS,DMATTS)
         RETURN
      END IF
C
      IF ( IECURR.NE.NETAB(1) .OR. MPI_ELOOP ) RETURN
C
C-----------------------------------------------------------------------
C  print out results only by MPI_ID = 0 for:
C  sequential:  IECURR = NETAB
C  MPI_KLOOP:   IECURR = NETAB
C  MPI_ELOOP:   IECURR.GT.NETAB(1)  add. dummy call, get here via  GOTO
C-----------------------------------------------------------------------
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C process MPI_ID = 0 continues here in case of MPI and IECURR > NETAB(1)
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
 100  CONTINUE
      FILNAM = DATSET(1:LDATSET)//'_J_ij.dat'
      LFN = LDATSET + 9
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
      WRITE (6,99003)
      WRITE (IOTMP,99003)
C
      WRITE (IOTMP,99005)
      DO IQ = IQBOT,IQTOP
         WRITE (IOTMP,99006) IQ,NOQ(IQ),
     &                       (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,NOQ(IQ)
     &                       )
      END DO
C
      NIJ = SUM(NQCLU_TAUIJ_CL(1:NCL))
C
      ALLOCATE (JXCITJT(NIJ,NT,NT),DITJT(NIJ,NT,NT),NITJT(NT,NT))
C
      NITJT(1:NT,1:NT) = 0
      JXCITJT_0(1:NT,1:NT) = 0.D0
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
C-------------------------------------------------- exclude on-site term
C
                     IF ( DQCLU_TAUIJ_CL(JQCLU,ICL).LT.1D-6 ) CYCLE
C
                     WRITE (6,99001) IQ,IT,JQ,JT,(QBAS(I,IQ),I=1,3),
     &                               (QBAS(I,JQ),I=1,3)
C
                     WRITE (IOTMP,99001) IQ,IT,JQ,JT,(QBAS(I,IQ),I=1,3),
     &                      (QBAS(I,JQ),I=1,3)
C
                     JXCIJ = DIMAG(J_IJINT(ITAUIJ,IT,JT))/(4*PI)
C
                     NITJT(IT,JT) = NITJT(IT,JT) + 1
                     I = NITJT(IT,JT)
                     JXCITJT(I,IT,JT) = JXCIJ*RY_EV*1D3
                     DITJT(I,IT,JT) = DQCLU_TAUIJ_CL(JQCLU,ICL)
C
                     JXCITJT_0(IT,JT) = JXCITJT_0(IT,JT)
     &                                  + JXCITJT(I,IT,JT)
C
C
                     WRITE (6,99002) ITAUIJ,ITAUJI,
     &                               N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                               RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                               DQCLU_TAUIJ_CL(JQCLU,ICL),JXCIJ,
     &                               JXCIJ*RY_EV
                     WRITE (IOTMP,99002) ITAUIJ,ITAUJI,
     &                      N5VEC_TAUIJ_CL(1:3,JQCLU,ICL),
     &                      RQCLU_TAUIJ_CL(1:3,JQCLU,ICL),
     &                      DQCLU_TAUIJ_CL(JQCLU,ICL),JXCIJ,JXCIJ*RY_EV
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
                     IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99011)
     &                    ROUTINE(1:LEN_TRIM(ROUTINE)),ITAUIJ,ITAUJI,
     &                    JXCIJ
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
                  END DO
               END DO
C
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
C
C------------------------------------------- Evaluate mean-field T_Curie
C
      LWMAX = 20*NT*NT
      ALLOCATE (W4(LWMAX))
C
      T_C(:,:) = 0.D0
      DO IQ = 1,NQ
C
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
C
            DO JQ = 1,NQ
               JXC0(IT) = 0D0
               DO IOP = 1,NOQ(JQ)
                  JT = ITOQ(IOP,JQ)
                  JXC0(IT) = CONC(JT)*JXCITJT_0(IT,JT)*1.D-3
                  T_C(IT,JT) = CONC(IT)*JXC0(IT)
               END DO
            END DO
         END DO
      END DO
C
      NTC = NT
C
      CALL DGEEV('N','N',NTC,T_C,NTC,WR,WI,WW1,NTC,WW2,NTC,W4,-1,INFO)
      LWK = MIN(LWMAX,INT(W4(1)))
      CALL DGEEV('N','N',NTC,T_C,NTC,WR,WI,WW1,NTC,WW2,NTC,W4,LWK,INFO)
C
      IF ( INFO.EQ.0 ) THEN
         TCURIE = 0
         DO IT = 1,NT
C
            TCURIE = MAX(TCURIE,T_C(IT,IT))
C
         END DO
         WRITE (6,99007) TCURIE*(2D0/3D0)*EV_J/KB_SI
         WRITE (6,99008)
         WRITE (6,99009)
         WRITE (6,99010)
C
      ELSE
         WRITE (6,*) 'WARNING:  <SGEGV> INFO = ',INFO
      END IF
C
C ----------------------------------------------------------------------
C
      WRITE (6,99004) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C                      find suitable plot windows
C-----------------------------------------------------------------------
C
      XMIN = 1D10
      XMAX = 0D0
      DO IT = ITBOT,ITTOP
         LEG(IT) = TXT_T(IT)
         DO JT = ITBOT,ITTOP
            DO I = 1,NITJT(IT,JT)
               XMIN = MIN(XMIN,DITJT(I,IT,JT))
               XMAX = MAX(XMAX,DITJT(I,IT,JT))
            END DO
         END DO
      END DO
C
C-------------------------------------------- start plot at X = R_ij = 0
      XMIN = 0D0
C
      DO IT = ITBOT,ITTOP
C
         YMIN = 0D0
         YMAX = 0D0
C
         DO JT = ITBOT,ITTOP
            DO I = 1,NITJT(IT,JT)
               YMIN = MIN(YMIN,JXCITJT(I,IT,JT))
               YMAX = MAX(YMAX,JXCITJT(I,IT,JT))
            END DO
         END DO
C
         CALL XMGRHEAD(DATSET,LDATSET,'Jij',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IOTMP,1,XMIN,1,XMAX,1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,
     &                 'distance (units of lattice parameter)',37,
     &                 'J!sij!N (meV)',13,' ',0,
     &                 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM),
     &                 25+LSYSTEM,
     &                 'exchange coupling parameters J!sij!N for '//
     &                 TXT_T(IT)(1:LTXT_T(IT))//' at center',
     &                 (41+LTXT_T(IT)+10),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,1,(ITTOP-ITBOT)+1,0,LEG,LEG)
C
         CALL XMGRPOINTS(IOTMP,1,(ITTOP-ITBOT)+1,0,2,1,1)
C
         DO JT = ITBOT,ITTOP
            CALL XMGRTABLE(0,(ITBOT-1),DITJT(1,IT,JT),JXCITJT(1,IT,JT),
     &                     1.0D0,NITJT(IT,JT),IOTMP)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '   exchange-coupling parameter written to file ',
     &               FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IOTMP)
C
      END DO
C
C-----------------------------------------------------------------------
C
99001 FORMAT (/,5X,'IQ =',I3,' IT =',I3,20X,' JQ =',I3,' JT =',I3,/,5X,
     &        2(3X,'->Q = (',F7.3,',',F7.3,',',F7.3,')',5X),/,5X,
     &        'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR ',
     &        '    J_ij [Ry]  J_ij [eV]')
99002 FORMAT (5X,2I7,2X,3I3,1X,3F7.3,1X,F7.3,1X,2F11.6)
99003 FORMAT (/,1X,79('*'),/,35X,'<XCPLJIJ>',/,28X,
     &        'XC-coupling constants J_ij',/,1X,79('*'),/)
99004 FORMAT (/,10X,'results written to file:',A,/)
99005 FORMAT (/,10X,'number of sites   NQ = ',I3,/,10X,
     &        'number of types   NT = ',I3,/,10X,'site occupation:')
99006 FORMAT (10X,2I4,10(I3,F6.3))
99007 FORMAT (/,5X,'Curie temperature within mean field',
     &        ' approximation  T_C =',F9.1,' K',/)
99008 FORMAT (/,5X,'NOTE 1:',/,5X,
     &        'For calculations of mean field T_C one has to use',/,5X,
     &        'the cluster size, CLURAD, as large as possible',/)
99009 FORMAT (/,5X,'NOTE 2:',/,5X,
     &        'The calculated temperature should be treated ',/,5X,
     &        'as Neel temperature, when the exchange parameters',/,5X,
     &        'between sub-lattices are negative: J_0^{\mu\nu} < 0',/)
99010 FORMAT (5X,
     &        'Mean-field Curie temperature for multicomponent system',
     &        /,5X,'is evaluated according to  P. W. Anderson [1]:',/,
     &        5X,
     &        '[1] P.W.Anderson, Theory of magn. exchange interactions:'
     &        ,/,5X,'Exchange in insulators and semiconductors.',/,5X,
     &        'Solid State Physics, edited by F. Seitz and D. Turnbull,'
     &        ,/,5X,'(Academic Press, New York), Vol. 14 pp. 99-214.',/)
99011 FORMAT ('# BUILDBOT: ',A,':  XC-coupling constants J_ij',
     &        '  for ITAUIJ,ITAUJI =',2I5,/,(1PE22.14))
      END
