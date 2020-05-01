C*==symtau.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMTAU(ITEST)
C   ********************************************************************
C   *                                                                  *
C   *  - find the table of coefficients that represent all symmetry    *
C   *    operations - only for IBZINT=3 (tetrahedron method)           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYM,SYMACCEPTED,NSYMACCEPTED,SYMUNITARY,
     &    DROT,IQORGQP
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKKR,NKMQ,NLMQ
      USE MOD_SITES,ONLY:NQ,NQMAX
      USE MOD_FILES,ONLY:IPRINT,IOTMP
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SYMMETRY,ONLY:SYMTAUTET_AVAILABLE
      USE MOD_KSPACE,ONLY:NELMT
      IMPLICIT NONE
C*--SYMTAU19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMTAU')
C
C Dummy arguments
C
      INTEGER ITEST
C
C Local variables
C
      REAL*8 ABST
      INTEGER I,I0,IA_ERR,IQ,ISYM,ITBZ(:),IW,IWR,J,JTBZ(:),K,L,LIN,N,
     &        NELMTMAX,NKMTOP,NLIN,NON0(:),NTAUUV(:),QTBZ(:),SNTAUUVMAX,
     &        UTAUUV(:),VTAUUV(:)
      COMPLEX*16 ST(:,:),TAUK(:,:),WTAUUV(:),X
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,ST,NON0
      ALLOCATABLE ITBZ,JTBZ,NTAUUV,QTBZ
      ALLOCATABLE WTAUUV,UTAUUV,VTAUUV
C
      ALLOCATE (NON0(NQMAX))
      ALLOCATE (TAUK(NKKR,NKKR),ST(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:symtau -> NON0'
C
C----------------------------------------------------- assume worst case
C                                       if SNTAUUVMAX is still too small
C                          a RESTART with then proper value will be done
C
      NELMTMAX = NKKR**2
      SNTAUUVMAX = 10*NELMTMAX
C
      ALLOCATE (ITBZ(NELMTMAX),JTBZ(NELMTMAX))
      ALLOCATE (QTBZ(NELMTMAX),NTAUUV(NELMTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: symtau -> QTBZ'
C
 100  CONTINUE
      ALLOCATE (WTAUUV(SNTAUUVMAX))
      ALLOCATE (UTAUUV(SNTAUUVMAX),VTAUUV(SNTAUUVMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: symtau -> WTAUUV'
C
      WRITE (6,99001)
C
      SYMTAUTET_AVAILABLE = .TRUE.
C
C=======================================================================
C       find the structure of the site-diagonal TAU - matrices  TAUQ
C=======================================================================
C
      NELMT = 0
      NLIN = 0
      IW = 6
C
      DO IQ = 1,NQ
         NON0(IQ) = 0
         IF ( IREL.EQ.2 ) THEN
            NKMTOP = NLMQ(IQ)
         ELSE
            NKMTOP = NKMQ(IQ)
         END IF
C
c         IF ( IPRINT.GT.0 ) WRITE (6,99004) IQ
         WRITE (6,99004) IQ
         DO I = 1,NKMTOP
            DO J = 1,NKMTOP
               ST(I,J) = 0.0D0
C
               CALL CINIT(NKKR*NKKR,TAUK)
C
               DO ISYM = 1,NSYM
                  IF ( SYMACCEPTED(ISYM) ) THEN
                     I0 = IND0Q(IQORGQP(ISYM,IQ))
C
                     IF ( SYMUNITARY(ISYM) ) THEN
                        DO L = 1,NKMTOP
                           DO K = 1,NKMTOP
                              TAUK(I0+K,I0+L) = TAUK(I0+K,I0+L)
     &                           + DROT(I,K,ISYM)*DCONJG(DROT(J,L,ISYM))
                           END DO
                        END DO
                     ELSE
                        DO L = 1,NKMTOP
                           DO K = 1,NKMTOP
                              TAUK(I0+L,I0+K) = TAUK(I0+L,I0+K)
     &                           + DROT(I,K,ISYM)*DCONJG(DROT(J,L,ISYM))
                           END DO
                        END DO
                     END IF
                  END IF
               END DO
C
               LIN = 0
               IWR = 0
c modified by XJQ: speed up in the case with many atoms and low symmetry
               DO ISYM = 1,NSYM
                  IF ( SYMACCEPTED(ISYM) ) THEN
                     I0 = IND0Q(IQORGQP(ISYM,IQ))
                     DO K=I0+1,I0+NKMTOP
                       DO L=I0+1,I0+NKMTOP
c               DO K = 1,NKKR
c                  DO L = 1,NKKR
                     ABST = ABS(TAUK(K,L))
                     ST(I,J) = ST(I,J) + ABST
                     IF ( ABST.GT.1D-8 ) THEN
                        X = TAUK(K,L)/DBLE(NSYMACCEPTED)
C
                        IF ( IPRINT.GT.1 ) THEN
                           IF ( IWR.EQ.0 ) THEN
                              IWR = 1
                              WRITE (IW,99002) I,J,IQ,X,K,L
                           ELSE
                              WRITE (IW,99003) X,K,L
                           END IF
                        END IF
                        LIN = LIN + 1
                        IF ( NLIN+LIN.LE.SNTAUUVMAX ) THEN
                           UTAUUV(NLIN+LIN) = K
                           VTAUUV(NLIN+LIN) = L
                           WTAUUV(NLIN+LIN) = X
                        END IF
                     END IF
C
c                  END DO
c               END DO
                       END DO
                     END DO
                  END IF
               END DO
c end-mod-xjq
C
               IF ( LIN.GT.0 ) THEN
                  NLIN = NLIN + LIN
                  NELMT = NELMT + 1
                  NON0(IQ) = NON0(IQ) + 1
                  IF ( NELMT.LE.NELMTMAX ) THEN
                     NTAUUV(NELMT) = LIN
                     ITBZ(NELMT) = I
                     JTBZ(NELMT) = J
                     QTBZ(NELMT) = IQ
                  END IF
               END IF
C
               IF ( ABS(ST(I,J)).GT.1D-5 ) ST(I,J) = 2
C
            END DO
         END DO
C
         IF ( IPRINT.GT.1 ) CALL CMATSTRUCT('TAU-MAT',ST,NKMTOP,NKMMAX,
     &        3,3,0,1D-8,-6)
      END DO
C
C-----------------------------------------------------------------------
      WRITE (6,99005) NELMT,(NON0(IQ),IQ=1,MIN(9,NQ))
      IF ( NQ.GT.9 ) WRITE (6,99006) (NON0(IQ),IQ=10,NQ)
      WRITE (6,99007) NELMTMAX,NLIN,SNTAUUVMAX
C      IF ( (NELMT.GT.NELMTMAX) .OR. (NLIN.GT.SNTAUUVMAX) ) THEN
      IF ( NLIN.GT.SNTAUUVMAX ) THEN
         WRITE (6,99009)
         IF ( ITEST.EQ.2 ) THEN
            NELMT = 0
         ELSE
            DEALLOCATE (UTAUUV,VTAUUV,WTAUUV,STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) STOP 'alloc: symtau -> WTAUUV'
            SNTAUUVMAX = NLIN
            GOTO 100
C            STOP
         END IF
      END IF
C
      IF ( (IPRINT.GT.1) .OR. (ITEST.EQ.2) ) THEN
         N = 0
         DO I = 1,NELMT
            WRITE (6,99008) ITBZ(I),JTBZ(I),QTBZ(I),NTAUUV(I),
     &                      (WTAUUV(N+J),UTAUUV(N+J),VTAUUV(N+J),J=1,
     &                      NTAUUV(I))
            N = N + NTAUUV(I)
         END DO
      END IF
C
C=======================================================================
C                   write results for rereading by main program
C=======================================================================
      N = 0
      DO I = 1,NELMT
         N = N + NTAUUV(I)
      END DO
      IF ( N.NE.NLIN ) STOP 'in <SYMTAU>   N <> NLIN'
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
      WRITE (IOTMP) NELMT,NLIN
      WRITE (IOTMP) (ITBZ(I),JTBZ(I),NTAUUV(I),QTBZ(I),I=1,NELMT)
      WRITE (IOTMP) (UTAUUV(I),VTAUUV(I),WTAUUV(I),I=1,NLIN)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (TAUK,ST,NON0,ITBZ,JTBZ,NTAUUV,QTBZ)
      DEALLOCATE (WTAUUV,UTAUUV,VTAUUV)
C
C-----------------------------------------------------------------------
99001 FORMAT (//,1X,79('*'),/,36X,'<SYMTAU>',/,1X,79('*'),//,10X,
     &        'BZ-integral   Int d^3k  TAU(k)[L1,L2] ')
99002 FORMAT ('     TAUQ(',I2,',',I2,',',I2,') =  ','   ',2F10.4,' * <',
     &        I3,'|T(K)|',I3,'>')
99003 FORMAT (23X,' + ',2F10.4,' * <',I3,'|T(K)|',I3,'>')
99004 FORMAT (//,1X,58('='),/,
     &        '   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=',I3,
     &        /,1X,58('='),/)
99005 FORMAT (/,10X,'non-0 TAU-elements          ',I5,'   Q:',9I4)
99006 FORMAT (43X,9I4)
99007 FORMAT (10X,'array size (NELMTMAX)  ',I10,/,10X,
     &        'terms to sum up        ',I10,/,10X,
     &        'array size (SNTAUUVMAX)',I10,/)
99008 FORMAT (' TAU(',2(I2,','),I2,')',I3,1X,5(2F6.3,2I3),:,/,
     &        (18X,5(2F6.3,2I3),:))
99009 FORMAT (10X,'WARNING: array size exceeded  >>>  ',
     &        'RESTART  with increased size ',/)
      END
