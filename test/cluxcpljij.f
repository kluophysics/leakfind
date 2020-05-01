C*==cluxcpljij.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUXCPLJIJ(IECURR,MSSQ,MSST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the site-off diagonal  XC-coupling parameters   J_ij  *
C   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:NONMAG,KMROT
      USE MOD_ANGMOM,ONLY:NLM,NKMMAX,NKMQ,NKM
      USE MOD_SITES,ONLY:NQMAX,IQAT,ITOQ,NOQ,NQ,QBAS
      USE MOD_ENERGY,ONLY:NETAB,NEPATH,WETAB,IGRID
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,CONC,NT,Z
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      IMPLICIT NONE
C*--CLUXCPLJIJ18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUXCPLJIJ')
C
C Dummy arguments
C
      INTEGER IECURR
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUM,DELMSST(NLM,NLM,NT),DMAMC(:,:),DMATTS(:,:,:,:),
     &           DTILTS(:,:,:,:),JXCIJINT(:,:,:,:),MQS(:,:,:),
     &           MSSTS(NLM,NLM,2),MTS(:,:,:),TAUQ(:,:,:),TAUQQ1(:,:,:,:)
     &           ,TAUQQ2(:,:,:,:),W1(:,:),W2(:,:),WSA(NLM,NLM),
     &           WSB(NLM,NLM)
      REAL*8 DNRM2
      REAL*8 DQVEC(3),DR,JXCIJ
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IFIL,IO,IQ,IT,J,JO,JQ,JT,LFN,M,N,NI,NJ
      SAVE JXCIJINT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1,W2,JXCIJINT
      ALLOCATABLE DMAMC,DTILTS,DMATTS
      ALLOCATABLE TAUQ,TAUQQ1,TAUQQ2,MQS,MTS
C
      IF ( NONMAG ) RETURN
      IF ( KMROT.NE.0 ) RETURN
      IF ( NEPATH.NE.1 ) RETURN
      IF ( IGRID(1).NE.5 ) RETURN
C
C
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX))
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (DMATTS(NLM,NLM,NT,2))
      ALLOCATE (DTILTS(NLM,NLM,NT,2))
      ALLOCATE (MTS(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (MQS(NKMMAX,NKMMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUXCPLJIJ->DMATT'
C
C----------------------------------------  read off-diagonal  TAU-matrix
C
      ALLOCATE (TAUQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TAUQQ1(NLM,NLM,NQ,NQ))
      ALLOCATE (TAUQQ2(NLM,NLM,NQ,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUXCPLJIJ->TAUQQ'
C
C
C
      IFIL = 19
      REWIND IFIL
C
      READ (IFIL) TAUQ
C
      DO JQ = 1,NQ
         NJ = NLM
         DO IQ = 1,NQ
            NI = NLM
C
            READ (IFIL) ((TAUQQ1(I,J,IQ,JQ),I=1,NI),J=1,NJ)
            READ (IFIL) ((TAUQQ2(I,J,IQ,JQ),I=1,NI),J=1,NJ)
C
         END DO
      END DO
      CLOSE (IFIL)
C
C-----------------------------------------------------------------------
C                 initialize all if  IECURR = 1
C-----------------------------------------------------------------------
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (JXCIJINT(NQ,NQ,NT,NT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc: CLUXCPLJIJ -> JXCIJINT'
C
         CALL CINIT(NQ*NQ*NT*NT,JXCIJINT)
C
      END IF
C-----------------------------------------------------------------------
C
      DO IT = 1,NT
C
         CALL CHANGEREP(NKM,NKMMAX,MSST(1,1,IT),'REL>RLM',MTS(1,1,IT))
C
         DO J = 1,NLM
            CALL ZCOPY(NLM,MTS(1,J,IT),1,MSSTS(1,J,1),1)
            CALL ZCOPY(NLM,MTS(NLM+1,NLM+J,IT),1,MSSTS(1,J,2),1)
         END DO
C
         DO J = 1,NLM
            DO I = 1,NLM
               DELMSST(I,J,IT) = MSSTS(I,J,2) - MSSTS(I,J,1)
            END DO
         END DO
C
C--------------------------- calculate projection matrices DMAT and DTIL
C---------------------------------- all matrices in  RLM  representation
C
         IQ = IQAT(1,IT)
         N = NKMQ(IQ)
         M = NKMMAX
C
         CALL CHANGEREP(NKM,NKMMAX,MSSQ(1,1,IQ),'REL>RLM',MQS(1,1,IQ))
C
         CALL GETDMAT(TAUQ(1,1,IQ),W1,W2,DMAMC,N,MQS(1,1,IQ),MTS(1,1,IT)
     &                ,M)
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
      DO IQ = 1,NQ
C
         DO JQ = 1,NQ
C
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               DO JO = 1,NOQ(JQ)
                  JT = ITOQ(JO,JQ)
C
                  CALL CMATMUL(N,M,TAUQQ1(1,1,JQ,IQ),DTILTS(1,1,IT,1),
     &                         WSA)
C
                  CALL CMATMUL(N,M,DMATTS(1,1,JT,1),WSA,WSB)
C
                  CALL CMATMUL(N,M,DELMSST(1,1,JT),WSB,WSA)
C
                  CALL CMATMUL(N,M,DTILTS(1,1,JT,2),WSA,WSB)
C
                  CALL CMATMUL(N,M,TAUQQ2(1,1,IQ,JQ),WSB,WSA)
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
                  JXCIJINT(IQ,JQ,IT,JT) = JXCIJINT(IQ,JQ,IT,JT)
     &               + WETAB(IECURR,1)*CSUM
C
               END DO
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      IF ( IECURR.NE.NETAB(1) ) THEN
         DEALLOCATE (W1,W2,DMAMC,DTILTS,DMATTS)
         DEALLOCATE (TAUQ,TAUQQ1,TAUQQ2,MQS,MTS)
         RETURN
      END IF
C
C-----------------------------------------------------------------------
C                              write results
C-----------------------------------------------------------------------
C
      FILNAM = DATSET(1:LDATSET)//'_J_ij.dat'
      LFN = LDATSET + 9
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
      WRITE (6,99003) (' ',I=1,6)
      WRITE (IOTMP,99003) ('#',I=1,6)
C
      WRITE (IOTMP,99005) NQ,NT
      DO IQ = 1,NQ
         WRITE (IOTMP,99006) IQ,NOQ(IQ),
     &                       (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,NOQ(IQ)
     &                       )
      END DO
      WRITE (IOTMP,99007)
C
      DO IQ = 1,NQ
         DO JQ = 1,NQ
            IF ( IQ.NE.JQ ) THEN
C
               DQVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
               DR = DNRM2(3,DQVEC,1)
C
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
                  IF ( Z(IT).NE.0 ) THEN
C
                     DO JO = 1,NOQ(JQ)
                        JT = ITOQ(JO,JQ)
                        IF ( Z(JT).NE.0 ) THEN
C
                           JXCIJ = DIMAG(JXCIJINT(IQ,JQ,IT,JT))/(4*PI)
C
                           WRITE (6,99001) IQ,IT,JQ,JT,
     &                            (QBAS(I,IQ),I=1,3),(QBAS(I,JQ),I=1,3)
                           WRITE (6,99002) DQVEC,DR,JXCIJ,JXCIJ*RY_EV
C
                           WRITE (IOTMP,99008) DR,JXCIJ*RY_EV*1D3,
     &                            JXCIJ*1D3,IQ,IT,TXT_T(IT),JQ,JT,
     &                            TXT_T(JT),(QBAS(I,IQ),I=1,3),
     &                            (QBAS(I,JQ),I=1,3)
C
                        END IF
                     END DO
C
                  END IF
               END DO
C
            END IF
         END DO
      END DO
C
      WRITE (6,99004) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (JXCIJINT,W1,W2,DMAMC,DTILTS,DMATTS)
      DEALLOCATE (TAUQ,TAUQQ1,TAUQQ2,MQS,MTS,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:CLUXCPLJIJ->JXCIJINT'
C
C-----------------------------------------------------------------------
99001 FORMAT (/,5X,'IQ =',I3,' IT =',I3,20X,' JQ =',I3,' JT =',I3,/,5X,
     &        2(3X,'->Q = (',F7.3,',',F7.3,',',F7.3,')',5X),/,5X,
     &        '      DRX    DRY    DRZ     DR ',
     &        '    J_ij [Ry]  J_ij [eV]')
99002 FORMAT (5X,10F11.6)
99003 FORMAT (1A,/,1A,79('*'),/,1A,32X,'<CLUXCPLJIJ>',/,1A,27X,
     &        'XC-coupling constants J_ij',/,1A,79('*'),/,1A)
99004 FORMAT (/,10X,'results written to file:',A,/)
99005 FORMAT ('#',/,'#',9X,'number of sites   NQ = ',I3,/,'#',9X,
     &        'number of types   NT = ',I3,/,'#',9X,'site occupation:')
99006 FORMAT ('#',9X,2I4,10(I3,F6.3))
99007 FORMAT ('#',/,'# DR   J (meV)   J (mRy)    ',
     &        ' IQ  IT       JQ  JT       QX(i)  QY(i)  QZ(i)',
     &        '  QX(j)  QY(j)  QZ(j)')
99008 FORMAT (F6.3,2F10.5,2(2I4,1X,A),3F7.2,3F7.2)
      END
