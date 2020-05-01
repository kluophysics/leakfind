C*==xcplj0.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XCPLJ0(TAUT,MSST,IECURR)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the site-diagonal  XC-coupling parameter   J0         *
C   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_SYMMETRY,ONLY:NQEQ,IQEQ
      USE MOD_CALCMODE,ONLY:NONMAG,KMROT
      USE MOD_SITES,ONLY:NQ,ITOQ,NOQ
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,CONC,NT
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM
      USE MOD_CONSTANTS,ONLY:C0,PI,EV_J,RY_EV,KB_SI
      USE MOD_TAUIJ,ONLY:NQTAB1_TAUIJ,IQ_QTAB1_TAUIJ
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      USE MOD_MPI,ONLY:MPI_ELOOP,MPI_KLOOP,MPI_ID
      IMPLICIT NONE
C*--XCPLJ021
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='XCPLJ0')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,CWK(:),DELMSS(NLM,NLM),DELTAU(NLM,NLM),JXC0INTT(:)
     &           ,MSSTS(NLM,NLM,2),TAUTS(NLM,NLM,2),WS1(NLM,NLM),
     &           WS2(NLM,NLM),WSA(NLM,NLM),WSB(NLM,NLM)
      INTEGER I,IA_ERR,IEQ,IO,IQ,IQTAB,IT,J,LWK,M,N
      LOGICAL INITIALIZE
      REAL*8 JXC0,JXC0T(NT),TCURIE
      SAVE JXC0INTT
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE JXC0INTT,CWK
C
      IF ( NONMAG ) RETURN
      IF ( KMROT.NE.0 ) RETURN
      IF ( NEPATH.NE.1 ) RETURN
      IF ( IGRID(1).NE.5 ) RETURN
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         ALLOCATE (JXC0INTT(NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: J_IJINT')
C
         JXC0INTT(:) = 0D0
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
         LWK = NT
         ALLOCATE (CWK(LWK))
C
         CALL DRV_MPI_REDUCE_C(JXC0INTT(1),CWK(1),LWK)
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
C-----------------------------------------------------------------------
C                       evaluate Eq. (18)
C-----------------------------------------------------------------------
C
      DO IT = 1,NT
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,MSST(1,J,IT),1,MSSTS(1,J,1),1)
            CALL ZCOPY(NLM,MSST(NLM+1,NLM+J,IT),1,MSSTS(1,J,2),1)
         END DO
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,TAUT(1,J,IT),1,TAUTS(1,J,1),1)
            CALL ZCOPY(NLM,TAUT(NLM+1,NLM+J,IT),1,TAUTS(1,J,2),1)
         END DO
C
         DO J = 1,NLM
            DO I = 1,NLM
               DELMSS(I,J) = MSSTS(I,J,2) - MSSTS(I,J,1)
               DELTAU(I,J) = TAUTS(I,J,2) - TAUTS(I,J,1)
            END DO
         END DO
C
         N = NLM
         M = NLM
C
         CALL CMATMUL(N,M,DELMSS,DELTAU,WSA)
         CALL CMATMUL(N,M,DELMSS,TAUTS(1,1,1),WS1)
         CALL CMATMUL(N,M,TAUTS(1,1,2),WS1,WS2)
         CALL CMATMUL(N,M,DELMSS,WS2,WSB)
C
         CSUM = C0
         DO I = 1,NLM
            CSUM = CSUM + WSA(I,I) + WSB(I,I)
         END DO
C
         JXC0INTT(IT) = JXC0INTT(IT) + WETAB(IECURR,1)*CSUM
C
      END DO
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
      WRITE (6,99001)
C
      DO IT = 1,NT
         JXC0T(IT) = -DIMAG(JXC0INTT(IT))/(4*PI)
         WRITE (6,99002) IT,TXT_T(IT),JXC0T(IT),JXC0T(IT)*RY_EV
      END DO
C
      WRITE (6,*)
C
C------------------------------------- average over types IT -- Eq. (36)
C
      DO IQTAB = 1,NQTAB1_TAUIJ
         IQ = IQ_QTAB1_TAUIJ(IQTAB)
C
         JXC0 = 0D0
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            JXC0 = JXC0 + CONC(IT)*JXC0T(IT)
         END DO
C
         WRITE (6,99003) IQ,JXC0,JXC0*RY_EV
         WRITE (6,99004) (ITOQ(IO,IQ),IO=1,NOQ(IQ))
         WRITE (6,99005) (IQEQ(IEQ,IQ),IEQ=1,NQEQ(IQ))
C
C----------------------------------------- Curie temperature -- Eq. (33)
C
         IF ( NQ.EQ.1 ) THEN
            TCURIE = (2D0/3D0)*JXC0*RY_EV*EV_J/KB_SI
            WRITE (6,99006) TCURIE
         END IF
C
      END DO
C
      WRITE (6,*)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
      IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99007)
     &                         ROUTINE(1:LEN_TRIM(ROUTINE)),IQ,JXC0
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
99001 FORMAT (/,1X,79('*'),/,36X,'<XCPLJ0>',/,28X,
     &        'XC-coupling constants J_0',/,1X,79('*'),/)
99002 FORMAT (5X,'IT = ',I3,2X,A,'  J0 = ',F10.6,' Ry  ',F10.6,' eV')
99003 FORMAT (5X,'IQ = ',I3,6X,'  J0 = ',F10.6,' Ry  ',F10.6,' eV')
99004 FORMAT (5X,'occupation     IT = ',20I3)
99005 FORMAT (5X,'equiv. sites   IQ = ',20I3)
99006 FORMAT (/,5X,'Curie temperature within mean field',
     &        ' approximation  T_C =',F9.1,' K',/)
99007 FORMAT ('# BUILDBOT: ',A,':  XC-coupling constants J_0',
     &        '  for IQ =',I5,/,(1PE22.14))
      END
