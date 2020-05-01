C*==tauij_host_invwri.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUIJ_HOST_INVWRI(IECURR,ERYD,TAUQ,MSSQ)
C   ********************************************************************
C   *                                                                  *
C   *   -  set up the cluster matrix TAU for the host                  *
C   *   -  invert TAU(HOST)                                            *
C   *   -  write  1 / TAU(HOST)                                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQMAX,NQCLU,IQBOT_HOST,IQTOP_HOST,QBAS_QCLU
      USE MOD_ANGMOM,ONLY:NKMMAX,IXM0_QCLU,NKKR_CLU,NSPIN,NXM,NXM_QCLU,
     &    NKM,L_LM,L_IKM
      USE MOD_SYMMETRY,ONLY:NSYM,SYMACCEPTED
      USE MOD_FILES,ONLY:IFILTAUIJ
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TAUIJ,ONLY:NTAUIJ,TAUIJ,ITAUIJ_LOOP
      IMPLICIT NONE
C*--TAUIJ_HOST_INVWRI19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAUIJ_HOST_INVWRI')
      LOGICAL CHECK_TAUIJ
      PARAMETER (CHECK_TAUIJ=.FALSE.)
      REAL*8 TOL_TAUIJ
      PARAMETER (TOL_TAUIJ=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IECURR
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      REAL*8 DVECIJ(3)
      INTEGER I,I0,IA_ERR,IINV,ILOOP,INFO,IPIV(:),IQ,IQCLU,IREC,IRELEFF,
     &        ISPIN,ITAUIJ,J,J0,JJ,JQCLU,LI,LJ,LWKOPT,NB,NI,NJ
      INTEGER ILAENV
      COMPLEX*16 TAU(:,:),TAUIJCLU(:,:,:,:),TIJ,TJI,WORK(:),XIJ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IPIV,WORK,TAU,TAUIJCLU
C
      CALL TRACK_INFO(ROUTINE)
C
      NB = ILAENV(1,'ZGETRI',' ',NKKR_CLU,-1,-1,-1)
      LWKOPT = NKKR_CLU*NB
C
      ALLOCATE (IPIV(NKKR_CLU))
      ALLOCATE (WORK(LWKOPT),TAU(NKKR_CLU,NKKR_CLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUTAUIJ_CREATE->WORK'
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( CHECK_TAUIJ ) THEN
         ALLOCATE (TAUIJCLU(NXM,NXM,NQCLU,NQCLU))
         IF ( IREL.LE.2 ) THEN
            IRELEFF = 1
         ELSE
            IRELEFF = 3
         END IF
         IINV = NSYM/2 + 1
      END IF
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DO ISPIN = 1,NSPIN
C
C-----------------------------------------------------------------------
C                          set up TAU-matrix
C-----------------------------------------------------------------------
C
         ILOOP = 0
         DO JQCLU = 1,NQCLU
            J0 = IXM0_QCLU(JQCLU)
            NJ = NXM_QCLU(JQCLU)
C
            DO IQCLU = 1,NQCLU
               I0 = IXM0_QCLU(IQCLU)
               NI = NXM_QCLU(IQCLU)
C
               ILOOP = ILOOP + 1
               ITAUIJ = ITAUIJ_LOOP(ILOOP)
C
               DO J = 1,NJ
                  JJ = J0 + J
                  DO I = 1,NI
                     TAU(I0+I,JJ) = TAUIJ(I,J,ISPIN,ITAUIJ)
                  END DO
               END DO
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               IF ( CHECK_TAUIJ ) THEN
C
                  TAUIJCLU(:,:,IQCLU,JQCLU) = TAUIJ(:,:,ISPIN,ITAUIJ)
C
                  DVECIJ(:) = QBAS_QCLU(:,IQCLU) - QBAS_QCLU(:,JQCLU)
C
                  WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),ISPIN,
     &                            ITAUIJ,ILOOP,IQCLU,JQCLU,DVECIJ
C
                  CALL CMATSTRUCT('TAUIJ-matrix',TAUIJ(1,1,ISPIN,ITAUIJ)
     &                            ,NXM,NXM,IRELEFF,IRELEFF,1,TOL_TAUIJ,
     &                            6)
C
               END IF
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
            END DO
         END DO
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK_TAUIJ .AND. IREL.LE.2 ) THEN
C
            DO JQCLU = 1,NQCLU
               NJ = NXM_QCLU(JQCLU)
C
               DO IQCLU = 1,JQCLU
                  NI = NXM_QCLU(IQCLU)
C
                  DO I = 1,NI
                     IF ( IREL.LE.2 ) THEN
                        LI = L_LM(I)
                     ELSE
                        LI = L_IKM(I)
                     END IF
                     DO J = 1,NJ
                        IF ( IREL.LE.2 ) THEN
                           LJ = L_LM(J)
                        ELSE
                           LJ = L_IKM(J)
                        END IF
C
                        TIJ = TAUIJCLU(I,J,IQCLU,JQCLU)
                        XIJ = TAUIJCLU(J,I,IQCLU,JQCLU)*((-1)**(LI-LJ))
                        TJI = TAUIJCLU(J,I,JQCLU,IQCLU)
C
C-------------------------------------- check for time reversal symmetry
C
                        IF ( ABS(TIJ-TJI).GT.TOL_TAUIJ ) WRITE (6,99002)
     &                       ROUTINE(1:LEN_TRIM(ROUTINE)),ISPIN,IQCLU,
     &                       JQCLU,I,J,TIJ,'TJI',TJI
C
C------------------------------------------ check for inversion symmetry
C
                        IF ( SYMACCEPTED(IINV) ) THEN
                           IF ( ABS(TIJ-XIJ).GT.TOL_TAUIJ )
     &                          WRITE (6,99002)
     &                          ROUTINE(1:LEN_TRIM(ROUTINE)),ISPIN,
     &                          IQCLU,JQCLU,I,J,TIJ,'XIJ',XIJ
                        END IF
C
                     END DO
                  END DO
C
               END DO
            END DO
         END IF
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C-----------------------------------------------------------------------
C                           symmetrize TAU-matrix
C-----------------------------------------------------------------------
C
         CALL TAUIJ_SYMMETRIZE(TAU)
C
C-----------------------------------------------------------------------
C                           invert TAU-matrix
C-----------------------------------------------------------------------
C
         CALL ZGETRF(NKKR_CLU,NKKR_CLU,TAU,NKKR_CLU,IPIV,INFO)
C
         CALL ZGETRI(NKKR_CLU,TAU,NKKR_CLU,IPIV,WORK,LWKOPT,INFO)
C
C-----------------------------------------------------------------------
C                       write  inverse TAU-matrix
C-----------------------------------------------------------------------
C
         IREC = 1 + NSPIN*(IECURR-1) + ISPIN
C
         WRITE (IFILTAUIJ,REC=IREC) ERYD,IQBOT_HOST,IQTOP_HOST,NTAUIJ,
     &                              NKM,NXM,ISPIN,
     &                              (TAUQ(1:NKM,1:NKM,IQ),MSSQ(1:NKM,
     &                              1:NKM,IQ),IQ=IQBOT_HOST,IQTOP_HOST),
     &                              TAU(1:NKKR_CLU,1:NKKR_CLU)
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
99001 FORMAT (//,10X,'<',A,'>  TAU[i,j] matrices for ',/,10X,'ISPIN =',
     &        I4,'  ITAUIJ =',I4,'  ILOOP =',I4,/,10X,'IQCLU =',I4,
     &        '  JQCLU =',I4,/,10X,'DVEC[i,j] = ',3F12.6,/)
99002 FORMAT (//,10X,'<',A,'>  TAU[i,j] matrices for ',/,10X,'ISPIN =',
     &        I4,/,10X,'IQCLU =',I4,'  JQCLU =',I4,'  I =',I4,'  J =',
     &        I4,'    DEVIATE !!!!!!!!!!!!!!!!!!!',/,10X,'TIJ = ',
     &        2F15.10,/10X,A,' = ',2F15.10,/)
      END
