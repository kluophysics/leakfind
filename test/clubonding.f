C*==clubonding.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUBONDING(IECURR,ERYD,MEZZ,MEZJ,MSSQ,MSST,TSST)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NT,NTMAX
      USE MOD_CONSTANTS,ONLY:PI,C0,C1
      USE MOD_ENERGY,ONLY:NEPATH
      USE MOD_ANGMOM,ONLY:NL,NMEMAX,NKMMAX,NKMQ,NKM,NLM
      USE MOD_SITES,ONLY:QBAS,NQMAX,IQAT,ITOQ,NOQ,NQ
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      IMPLICIT NONE
C*--CLUBONDING15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IECURR
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      LOGICAL CHECKTAUDIA
      COMPLEX*16 DMAMC(:,:),DMATT(:,:,:),DTILT(:,:,:),G0QQ(:,:,:,:),
     &           MQS(:,:,:),MTS(:,:,:),TAUQ(:,:,:),TAUQQDD(:,:,:,:),
     &           TAUQQDU(:,:,:,:),TAUQQUD(:,:,:,:),TAUQQUU(:,:,:,:),
     &           TTS(:,:,:),W1(:,:),W2(:,:),WGG(:,:),WGIJ(:,:),WGJI(:,:)
     &           ,WTG(:,:),WTIJ(:,:),WZJ(:,:),WZZ(:,:)
      REAL*8 DNRM2
      REAL*8 DQVEC(3),DR,DSH(:),N00(0:NL,2),NBS(0:NL,2),NGG(0:NL,2),
     &       NSS(0:NL,2)
      INTEGER FILG0,FILTAU,I,IA_ERR,IDN,IL,IML,IO,IQ,IS,ISH,ISHQ(:),IT,
     &        IUP,J,JDN,JO,JQ,JSH,JT,JUP,M,N,NI,NJ,NSH
      SAVE DSH,ISHQ,NSH
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1,W2,WGG,WGIJ,WTIJ,WGJI,WZJ,WZZ,DSH,ISHQ,WTG
      ALLOCATABLE DMAMC,DTILT,DMATT,G0QQ
      ALLOCATABLE TAUQ,TAUQQDD,TAUQQUU,TAUQQDU,TAUQQUD,MQS,MTS,TTS
C
      IF ( KMROT.NE.0 ) RETURN
      IF ( NEPATH.NE.1 ) RETURN
C
      M = NKMMAX
      ALLOCATE (WTIJ(M,M),WGIJ(M,M),WGJI(M,M),WGG(M,M))
      ALLOCATE (WZZ(M,M),WZJ(M,M),W1(M,M),W2(M,M),WTG(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUBONDING->W1'
C
      ALLOCATE (DMAMC(M,M),DMATT(M,M,NT),DTILT(M,M,NT))
      ALLOCATE (TTS(M,M,NT),MTS(M,M,NT),MQS(M,M,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUBONDING->DMATT'
C
      CHECKTAUDIA = .TRUE.
C
C---------------------------------  read off-diagonal G0- and TAU-matrix
C
      ALLOCATE (TAUQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (G0QQ(NLM,NLM,NQ,NQ))
      ALLOCATE (TAUQQDD(NLM,NLM,NQ,NQ))
      ALLOCATE (TAUQQUU(NLM,NLM,NQ,NQ))
      ALLOCATE (TAUQQDU(NLM,NLM,NQ,NQ))
      ALLOCATE (TAUQQUD(NLM,NLM,NQ,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:CLUBONDING->TAUQQ'
C
      FILG0 = 18
      FILTAU = 19
C
      REWIND FILG0
      REWIND FILTAU
C
      READ (FILTAU) TAUQ
C
      DO JQ = 1,NQ
         NJ = NLM
         DO IQ = 1,NQ
            NI = NLM
C
            READ (FILG0) ((G0QQ(I,J,IQ,JQ),I=1,NI),J=1,NJ)
C
            IF ( IREL.LE.2 ) THEN
               READ (FILTAU) ((TAUQQDD(I,J,IQ,JQ),I=1,NI),J=1,NJ)
               READ (FILTAU) ((TAUQQUU(I,J,IQ,JQ),I=1,NI),J=1,NJ)
               TAUQQDU(:,:,IQ,JQ) = C0
               TAUQQUD(:,:,IQ,JQ) = C0
            ELSE
               READ (FILTAU) ((TAUQQDD(I,J,IQ,JQ),I=1,NI),J=1,NJ)
               READ (FILTAU) ((TAUQQDU(I,J,IQ,JQ),I=1,NI),J=1,NJ)
               READ (FILTAU) ((TAUQQUD(I,J,IQ,JQ),I=1,NI),J=1,NJ)
               READ (FILTAU) ((TAUQQUU(I,J,IQ,JQ),I=1,NI),J=1,NJ)
            END IF
C
         END DO
      END DO
      CLOSE (FILG0)
      CLOSE (FILTAU)
C
C-----------------------------------------------------------------------
      IF ( CHECKTAUDIA ) THEN
C
         DO IQ = 1,NQ
C
            DO J = 1,NLM
               JDN = J
               JUP = NLM + J
               DO I = 1,NLM
                  IDN = I
                  IUP = NLM + I
                  W1(IDN,JDN) = TAUQQDD(I,J,IQ,IQ)
                  W1(IDN,JUP) = TAUQQDU(I,J,IQ,IQ)
                  W1(IUP,JDN) = TAUQQUD(I,J,IQ,IQ)
                  W1(IUP,JUP) = TAUQQUU(I,J,IQ,IQ)
               END DO
            END DO
C
            DO J = 1,2*NLM
               DO I = 1,2*NLM
                  IF ( ABS(TAUQ(I,J,IQ)-W1(I,J)).GT.1D-8 ) THEN
                     WRITE (6,*) 'TROUBLE..... IQ I J ',IQ,I,J
                     WRITE (6,*) 'TAU ',TAUQ(I,J,IQ)
                     WRITE (6,*) 'W1  ',W1(I,J)
                  END IF
               END DO
            END DO
C
         END DO
C
      END IF
C
C-----------------------------------------------------------------------
C                 initialize all if  IECURR = 1
C-----------------------------------------------------------------------
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (DSH(NQ),ISHQ(NQ),STAT=IA_ERR)
C
         NSH = 0
         DSH(1) = 0D0
         IQ = 1
C
         DO JQ = 2,NQ
C
            DQVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
C
            DR = DNRM2(3,DQVEC,1)
C
            DO ISH = 1,NSH
               IF ( ABS(DR-DSH(ISH)).LT.1D-6 ) THEN
                  ISHQ(JQ) = ISH
                  GOTO 20
               END IF
C
               IF ( DSH(ISH).GT.DR ) THEN
                  DO JSH = NSH,ISH, - 1
                     DSH(JSH+1) = DSH(JSH)
                  END DO
                  DSH(ISH) = DR
                  ISHQ(JQ) = ISH
                  NSH = NSH + 1
                  GOTO 20
               END IF
            END DO
C
            NSH = NSH + 1
            DSH(NSH) = DR
            ISHQ(JQ) = NSH
C
 20         CONTINUE
            WRITE (6,*) '***************  jq - shell ',JQ,ISHQ(JQ)
         END DO
C
         DO ISH = 1,NSH
            WRITE (6,*) ' shell ',ISH,DSH(ISH)
         END DO
         DO IQ = 2,NQ
            WRITE (6,*) ' iq - shell ',IQ,ISHQ(IQ)
         END DO
C
      END IF
C-----------------------------------------------------------------------
C
      N = NKM
      M = NKMMAX
C
C--------------------------- calculate projection matrices DMAT and DTIL
C---------------------------------- all matrices in  CLM  representation
C
      DO IQ = 1,NQ
C
         CALL CHANGEREP(NKM,NKMMAX,MSSQ(1,1,IQ),'REL>CLM',MQS(1,1,IQ))
C
      END DO
C
      DO IT = 1,NT
C
         CALL CHANGEREP(NKM,NKMMAX,TSST(1,1,IT),'REL>CLM',TTS(1,1,IT))
C
         CALL CHANGEREP(NKM,NKMMAX,MSST(1,1,IT),'REL>CLM',MTS(1,1,IT))
C
         IQ = IQAT(1,IT)
         N = NKMQ(IQ)
C
         CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT),DMAMC,N,
     &                MQS(1,1,IQ),MTS(1,1,IT),M)
C
C            CALL  CMATSTR('DMATT A ',DMATT(1,1,IT),NKM,NKMMAX,3,3,0,
C     &                   1D-8,6)
      END DO
C
C=======================================================================
C=======================================================================
C
      IQ = 1
C
      DO IO = 1,NOQ(IQ)
         IT = ITOQ(IO,IQ)
C
         DO J = 1,NKM
            DO I = J + 1,NKM
               MEZJ(I,J,IT,1) = C0
               MEZJ(J,I,IT,1) = C0
            END DO
         END DO
C
         CALL CHANGEREP(NKM,NKMMAX,MEZZ(1,1,IT,1),'REL>CLM',WZZ)
C
         CALL CHANGEREP(NKM,NKMMAX,MEZJ(1,1,IT,1),'REL>CLM',WZJ)
C
         CALL CMATMUL(N,M,TAUQ(1,1,IQ),DTILT(1,1,IT),W1)
C
         CALL ZGEMM('N','N',N,N,N,CPRE,WZZ,M,W1,M,C0,W2,M)
C
         DO J = 1,N
            W2(J,J) = W2(J,J) - CPRE*WZJ(J,J)
         END DO
C
         I = 0
         DO IS = 1,2
            N00(0,IS) = 0D0
            DO IL = 1,NL
               N00(IL,IS) = 0D0
               DO IML = 1,2*IL - 1
                  I = I + 1
                  N00(IL,IS) = N00(IL,IS) + DIMAG(W2(I,I))
               END DO
               N00(0,IS) = N00(0,IS) + N00(IL,IS)
               WRITE (6,'(2i3,2f10.4)') IL,IS,N00(0,IS)
            END DO
         END DO
C
         CALL ZGEMM('N','N',N,N,N,CPRE,WZZ,M,TTS(1,1,IT),M,C0,W2,M)
C
         DO J = 1,N
            W2(J,J) = W2(J,J) - CPRE*WZJ(J,J)
         END DO
C
         I = 0
         DO IS = 1,2
            NSS(0,IS) = 0D0
            DO IL = 1,NL
               NSS(IL,IS) = 0D0
               DO IML = 1,2*IL - 1
                  I = I + 1
                  NSS(IL,IS) = NSS(IL,IS) + DIMAG(W2(I,I))
               END DO
               NSS(0,IS) = NSS(0,IS) + NSS(IL,IS)
               WRITE (6,'(2i3,2f10.4)') IL,IS,NSS(0,IS)
            END DO
         END DO
C
C
         CALL CINIT(NKMMAX*NKMMAX,WTG)
         CALL CINIT(NKMMAX*NKMMAX,WGG)
C
         DO JQ = 2,NQ
C
            DO J = 1,NLM
               JDN = J
               JUP = NLM + J
               DO I = 1,NLM
                  IDN = I
                  IUP = NLM + I
                  WTIJ(IDN,JDN) = TAUQQDD(I,J,IQ,JQ)
                  WTIJ(IDN,JUP) = TAUQQDU(I,J,IQ,JQ)
                  WTIJ(IUP,JDN) = TAUQQUD(I,J,IQ,JQ)
                  WTIJ(IUP,JUP) = TAUQQUU(I,J,IQ,JQ)
C
                  WGIJ(IDN,JDN) = G0QQ(I,J,IQ,JQ)
                  WGIJ(IDN,JUP) = 0D0
                  WGIJ(IUP,JDN) = 0D0
                  WGIJ(IUP,JUP) = G0QQ(I,J,IQ,JQ)
C
                  WGJI(IDN,JDN) = G0QQ(I,J,JQ,IQ)
                  WGJI(IDN,JUP) = 0D0
                  WGJI(IUP,JDN) = 0D0
                  WGJI(IUP,JUP) = G0QQ(I,J,JQ,IQ)
               END DO
            END DO
C
            CALL ZGEMM('N','N',N,N,N,C1,WTIJ,M,WGJI,M,C1,WTG,M)
C
            JT = ITOQ(1,JQ)
C
            CALL CMATMUL(N,M,WGIJ,TTS(1,1,JT),W2)
            CALL CMATMUL(N,M,TTS(1,1,IT),W2,WGIJ)
            CALL ZGEMM('N','N',N,N,N,C1,WGIJ,M,WGJI,M,C1,WGG,M)
C
         END DO
C
         CALL CMATMUL(N,M,WTG,TTS(1,1,IT),W1)
         CALL ZGEMM('N','N',N,N,N,CPRE,WZZ,M,W1,M,C0,WTG,M)
C
         CALL CMATMUL(N,M,WGG,TTS(1,1,IT),W1)
         CALL ZGEMM('N','N',N,N,N,CPRE,WZZ,M,W1,M,C0,WGG,M)
C
         I = 0
         DO IS = 1,2
            NBS(0,IS) = 0D0
            NGG(0,IS) = 0D0
            DO IL = 1,NL
               NBS(IL,IS) = 0D0
               NGG(IL,IS) = 0D0
               DO IML = 1,2*IL - 1
                  I = I + 1
                  NBS(IL,IS) = NBS(IL,IS) + DIMAG(WTG(I,I))
                  NGG(IL,IS) = NGG(IL,IS) + DIMAG(WGG(I,I))
               END DO
               NBS(0,IS) = NBS(0,IS) + NBS(IL,IS)
               NGG(0,IS) = NGG(0,IS) + NGG(IL,IS)
               WRITE (6,'(2i3,2f10.4)') IL,IS,NBS(0,IS),NGG(0,IS)
            END DO
         END DO
C
C
         WRITE (6,'(2i3,10f10.4)')
         DO IS = 1,2
            DO IL = 1,NL
               WRITE (6,'(2i3,10f10.4)') IL,IS,N00(IL,IS),NSS(IL,IS)
     &                + NBS(IL,IS),NSS(IL,IS),NBS(IL,IS),NGG(IL,IS)
            END DO
         END DO
C
         WRITE (22,'(20f12.5)') DREAL(ERYD),
     &                          ((N00(IL,IS),IL=0,0),(NSS(IL,IS),IL=0,0)
     &                          ,(NBS(IL,IS),IL=0,0),(NGG(IL,IS),IL=0,0)
     &                          ,IS=1,2)
C
      END DO
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
Cc      FILNAM = DATSET(1:LDATSET)//'_J_ij.dat'
Cc      LFN = LDATSET + 9
Cc      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
C
      IQ = 1
      DO IO = 1,NOQ(IQ)
         IT = ITOQ(IO,IQ)
C
         DO JQ = 2,NQ
C
            DQVEC(1:3) = QBAS(1:3,IQ) - QBAS(1:3,JQ)
C
            DR = DNRM2(3,DQVEC,1)
C
            DO JO = 1,NOQ(JQ)
               JT = ITOQ(JO,JQ)
C
            END DO
C
         END DO
C
      END DO
C
Ccccccccccc      WRITE (6,99001) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (W1,W2,WGG,WGIJ,WTIJ,WGJI,WZJ,WZZ,WTG)
      DEALLOCATE (DMAMC,DTILT,DMATT,G0QQ)
      DEALLOCATE (TAUQ,TAUQQDD,TAUQQUU,TAUQQDU,TAUQQUD,MQS,MTS,TTS)
C
C99001 FORMAT (/,10X,'results written to file:',A,/)
      END
