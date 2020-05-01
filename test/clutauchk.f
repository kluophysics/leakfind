C*==clutauchk.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CLUTAUCHK(IPRINT,IREL,NKKR,NLM,IQCLU1,IQCLU2,NLMQCLU,
     &                     NKMQCLU,CONC,WA,WB,WC,WD,IPIV,ITCPAMAX,
     &                     NOQCLU,ITOQCLU,CPATOL,NCPA,ICPACLU,CPACHNG,
     &                     ICPAFLAG,NSYM,DROT,IQCLUORG,SYMACCEPTED,
     &                     SYMUNITARY,SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                     TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,TAUQCLU)
C   ********************************************************************
C   *                                                                  *
C   *   check     cluster TAU-matrix                                   *
C   *                                                                  *
C   *   - calculate TAU inverting the full KKR-matrix                  *
C   *   - make use of simplifications due to matrix structure          *
C   *   - compare results                                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--CLUTAUCHK22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CLUTAUCHK')
      REAL*8 TOLR,TOL
      PARAMETER (TOLR=1.0D-8,TOL=1.0D-8)
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPATOL
      INTEGER ICPAFLAG,IPRINT,IQCLU1,IQCLU2,IREL,ITCPAMAX,NCPA,NKKR,
     &        NKMMAX,NLM,NQCLU,NSYM,NTMAX
      LOGICAL SYMMETRIZE
      REAL*8 CONC(NTMAX)
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),
     &           MSSQCLU(NKMMAX,NKMMAX,NQCLU),
     &           MSSTCLU(NKMMAX,NKMMAX,NTMAX),
     &           TAUQCLU(NKMMAX,NKMMAX,NQCLU),
     &           TSSQCLU(NKMMAX,NKMMAX,NQCLU),
     &           TSSTCLU(NKMMAX,NKMMAX,NTMAX),WA(NKKR,NKKR),
     &           WB(NKKR,NKKR),WC(NKKR,NKKR),WD(NKKR,NKKR)
      INTEGER ICPACLU(NQCLU),IPIV(NKKR),IQCLUORG(NSYMMAX,NQCLU),
     &        ITOQCLU(NTMAX,NQCLU),NKMQCLU(NQCLU),NLMQCLU(NQCLU),
     &        NOQCLU(NQCLU)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      INTEGER I,IA_ERR,II0,II1,IINV,IIPIV(:),INFO,IQA,IQB,IQCLU,IQCLUP,
     &        ISYM,IUP1,J,JDN,JJ,JUP,KIJ,M,N,NERR5,NERROR,NKM,NN2,NN2SQ,
     &        NSYMACPTLOC
      COMPLEX*16 MG0MAT(:,:),T1,T2,TT(:,:,:),W1(NKMMAX,NKMMAX),
     &           W2(NKMMAX,NKMMAX),WWA(:,:),WWB(:,:)
      LOGICAL SYMACPTLOC(NSYMMAX)
      REAL*8 TA1,TA2,TI1,TI2,TIME,TIME0,TR1,TR2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MG0MAT,WWA,WWB,IIPIV,TT
C
      SYMMETRIZE = .TRUE.
C      symmetrize = .false.
      NKM = 2*NLM
      IINV = NSYM/2 + 1
C
      NSYMACPTLOC = 0
      DO ISYM = 1,NSYM
C         IF ( SYMACCEPTED(ISYM) .AND. SYMUNITARY(ISYM) ) THEN
         IF ( SYMACCEPTED(ISYM) ) THEN
            SYMACPTLOC(ISYM) = .TRUE.
            NSYMACPTLOC = NSYMACPTLOC + 1
C
            IF ( IINV.EQ.0 ) THEN
               IF ( .NOT.SYMUNITARY(ISYM) ) THEN
C     &       call zscal(nkmmax*nkmmax,CI,DROT(1,1,isym) )
C
                  WRITE (6,*) '############## ISYM',ISYM
                  CALL CMATSTRUCT('DROT  ANTIUNITARY   ',DROT(1,1,ISYM),
     &                            NKM,NKMMAX,2,2,0,1D-8,6)
C
                  CALL ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM),
     &                       NKMMAX,DROT(1,1,IINV),NKMMAX,C0,W1,NKMMAX)
C             call cmatcop(nkm,nkmmax,w1,DROT(1,1,ISYM))
C
                  WRITE (6,*) '############## ISYM',ISYM
                  CALL CMATSTRUCT('DROT  ANTIUNITARY   ',DROT(1,1,ISYM),
     &                            NKM,NKMMAX,2,2,0,1D-8,6)
               END IF
            END IF
C
         ELSE
            SYMACPTLOC(ISYM) = .FALSE.
         END IF
      END DO
      CALL CMATSTRUCT('Inversion     MATRIX',DROT(1,1,IINV),NKM,NKMMAX,
     &                2,2,0,1D-8,6)
C     stop
C-----------------------------------------------------------------------
C
      NCPA = 0
C
      WRITE (6,99003) SYMMETRIZE
C
C-------------------------------------------------------- store G-matrix
C
      IF ( IREL.GE.2 ) THEN
C
         ALLOCATE (MG0MAT(NKKR,NKKR),STAT=IA_ERR)
C
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MG0MAT')
C
         CALL ZCOPY(NKKR*NKKR,WA,1,MG0MAT,1)
C
      END IF
C
C=======================================================================
C         non-/scalar-relativistic AND non-spin-polarized calculation
C=======================================================================
      IF ( IREL.LE.1 ) THEN
C
         CALL STOP_MESSAGE(ROUTINE,'IREL.LE.1')
C
C=======================================================================
C         non-/scalar-relativistic AND spin-polarized calculation
C=======================================================================
C             fully relativistic / spin-polarized calculation
C=======================================================================
      ELSE IF ( IREL.GE.2 ) THEN
C
         NN2 = 2*NKKR
         ALLOCATE (TT(NKMMAX,NKMMAX,NQCLU))
         ALLOCATE (WWA(NN2,NN2),WWB(NN2,NN2),IIPIV(NN2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: WWA')
C
         NN2SQ = NN2*NN2
C
C-----------------------------------------------------------------------
         CALL CPU_TIME(TIME0)
C
         CALL CINIT(NN2SQ,WWA)
C
         DO J = 1,NKKR
            DO I = 1,NKKR
               WWA(I,J) = WA(I,J)
               WWA(NKKR+I,NKKR+J) = WA(I,J)
            END DO
         END DO
C
         II0 = 0
         DO IQCLU = 1,NQCLU
            N = NLMQCLU(IQCLU)
            IUP1 = NLM + 1
            II1 = II0 + 1
            DO JDN = 1,N
               JUP = NLM + JDN
               JJ = II0 + JDN
C                                                                  dn-dn
               CALL ZAXPY(N,C1,MSSQCLU(1,JDN,IQCLU),1,WWA(II1,JJ),1)
C                                                                  dn-up
               CALL ZAXPY(N,C1,MSSQCLU(1,JUP,IQCLU),1,WWA(II1,NKKR+JJ),
     &                    1)
C                                                                  up-dn
               CALL ZAXPY(N,C1,MSSQCLU(IUP1,JDN,IQCLU),1,
     &                    WWA(NKKR+II1,JJ),1)
C                                                                  up-up
               CALL ZAXPY(N,C1,MSSQCLU(IUP1,JUP,IQCLU),1,
     &                    WWA(NKKR+II1,NKKR+JJ),1)
            END DO
            II0 = II0 + N
         END DO
C
         CALL ZGETRF(NN2,NN2,WWA,NN2,IIPIV,INFO)
         CALL ZGETRI(NN2,WWA,NN2,IIPIV,WWB,NN2SQ,INFO)
C
         CALL CINIT(NKMMAX*NKMMAX*NQCLU,TT)
C
         II0 = 0
         DO IQCLU = 1,NQCLU
            N = NLMQCLU(IQCLU)
            IUP1 = NLM + 1
            II1 = II0 + 1
            DO JDN = 1,N
               JUP = NLM + JDN
               JJ = II0 + JDN
C                                                                  dn-dn
               CALL ZCOPY(N,WWA(II1,JJ),1,TT(1,JDN,IQCLU),1)
C                                                                  dn-up
               CALL ZCOPY(N,WWA(II1,NKKR+JJ),1,TT(1,JUP,IQCLU),1)
C                                                                  up-dn
               CALL ZCOPY(N,WWA(NKKR+II1,JJ),1,TT(IUP1,JDN,IQCLU),1)
C                                                                  up-up
               CALL ZCOPY(N,WWA(NKKR+II1,NKKR+JJ),1,TT(IUP1,JUP,IQCLU),
     &                    1)
            END DO
            II0 = II0 + N
         END DO
C
         CALL CPU_TIME(TIME)
         WRITE (6,99004) 'LAPACK             ',TIME - TIME0
         DEALLOCATE (WWB)
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'clutauchk_LAPACK')
         DO IQCLU = 1,NQCLU
            WRITE (IOTMP,99001) IQCLU
            CALL CMATSTRUCT('LAPACK',TT(1,1,IQCLU),NKMQCLU(IQCLU),
     &                      NKMMAX,2,2,0,1D-8,IOTMP)
         END DO
C
C-----------------------------------------------------------------------
         IF ( IREL.EQ.3 ) THEN
C
            CALL CPU_TIME(TIME0)
            KIJ = 0
C
            CALL CLUTAUMAT(KIJ,IPRINT,IREL,NKKR,NLM,IQCLU1,IQCLU2,
     &                     NLMQCLU,NKMQCLU,CONC,WA,WB,WC,WD,IPIV,
     &                     ITCPAMAX,NOQCLU,ITOQCLU,CPATOL,NCPA,ICPACLU,
     &                     CPACHNG,ICPAFLAG,NSYM,NSYMACPTLOC,DROT,
     &                     IQCLUORG,SYMACPTLOC,SYMUNITARY,SYMMETRIZE,
     &                     NKMMAX,NQCLU,NTMAX,TSSTCLU,MSSTCLU,TSSQCLU,
     &                     MSSQCLU,TAUQCLU)
C
            CALL CPU_TIME(TIME)
            WRITE (6,99004) '<CLUTAUMAT>  KIJ=0 ',TIME - TIME0
C
            NERROR = 0
C
            DO IQCLU = 1,NQCLU
               N = NLMQCLU(IQCLU)
               DO J = 1,2*N
                  DO I = 1,2*N
                     T1 = TT(I,J,IQCLU)
                     T2 = TAUQCLU(I,J,IQCLU)
                     TR1 = DREAL(T1)
                     TI1 = DIMAG(T1)
                     TA1 = ABS(TR1) + ABS(TI1)
                     TR2 = DREAL(T2)
                     TI2 = DIMAG(T2)
                     TA2 = ABS(TR2) + ABS(TI2)
C
                     IF ( TA1.GT.1D-6 ) THEN
                        IF ( ABS(1D0-TA2/TA1).GT.TOL ) THEN
                           WRITE (6,99005) IQCLU,I,J,T1,T2,'A'
                           NERROR = NERROR + 1
                        END IF
                     ELSE IF ( ABS(TA2-TA1).GT.TOL ) THEN
                        WRITE (6,99005) IQCLU,I,J,T1,T2,'B'
                        NERROR = NERROR + 1
                     END IF
                  END DO
               END DO
            END DO
C
            WRITE (6,99006) KIJ,NERROR,TOL
C
            CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'clutauchk_KIJ0')
            DO IQCLU = 1,NQCLU
               WRITE (IOTMP,99001) IQCLU
               CALL CMATSTRUCT('KIJ = 0',TAUQCLU(1,1,IQCLU),
     &                         NKMQCLU(IQCLU),NKMMAX,2,2,0,1D-8,IOTMP)
            END DO
C
            CALL ZCOPY(NKKR*NKKR,MG0MAT,1,WA,1)
C
         END IF
C-----------------------------------------------------------------------
C
         CALL CPU_TIME(TIME0)
         KIJ = 1
C
         CALL CLUTAUMAT(KIJ,IPRINT,IREL,NKKR,NLM,IQCLU1,IQCLU2,NLMQCLU,
     &                  NKMQCLU,CONC,WA,WB,WC,WD,IPIV,ITCPAMAX,NOQCLU,
     &                  ITOQCLU,CPATOL,NCPA,ICPACLU,CPACHNG,ICPAFLAG,
     &                  NSYM,NSYMACPTLOC,DROT,IQCLUORG,SYMACPTLOC,
     &                  SYMUNITARY,SYMMETRIZE,NKMMAX,NQCLU,NTMAX,
     &                  TSSTCLU,MSSTCLU,TSSQCLU,MSSQCLU,TAUQCLU)
C
         CALL CPU_TIME(TIME)
         WRITE (6,99004) '<CLUTAUMAT>  KIJ=1 ',TIME - TIME0
C
         NERROR = 0
C
         DO IQCLU = 1,NQCLU
            N = NLMQCLU(IQCLU)
            DO J = 1,2*N
               DO I = 1,2*N
                  T1 = TT(I,J,IQCLU)
                  T2 = TAUQCLU(I,J,IQCLU)
                  TR1 = DREAL(T1)
                  TI1 = DIMAG(T1)
                  TA1 = ABS(TR1) + ABS(TI1)
                  TR2 = DREAL(T2)
                  TI2 = DIMAG(T2)
                  TA2 = ABS(TR2) + ABS(TI2)
C
                  IF ( TA1.GT.1D-6 ) THEN
                     IF ( ABS(1D0-TA2/TA1).GT.TOL ) THEN
                        WRITE (6,99005) IQCLU,I,J,T1,T2,'A'
                        NERROR = NERROR + 1
                     END IF
                  ELSE IF ( ABS(TA2-TA1).GT.TOL ) THEN
                     WRITE (6,99005) IQCLU,I,J,T1,T2,'B'
                     NERROR = NERROR + 1
                  END IF
C
C IF ( NERROR.GT.4 ) CALL STOP_MESSAGE(ROUTINE,'NERROR > 4 ')
               END DO
            END DO
         END DO
C
         WRITE (6,99006) KIJ,NERROR,TOL
C
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'clutauchk_KIJ1')
         DO IQCLU = 1,NQCLU
            WRITE (IOTMP,99001) IQCLU
            CALL CMATSTRUCT('KIJ = 1',TAUQCLU(1,1,IQCLU),NKMQCLU(IQCLU),
     &                      NKMMAX,2,2,0,1D-8,IOTMP)
         END DO
C
C-----------------------------------------------------------------------
         M = NKMMAX
         NERR5 = 0
C
         DO ISYM = 1,NSYM
            IF ( SYMACPTLOC(ISYM) ) THEN
C
               DO IQCLU = 1,NQCLU
C
                  IQCLUP = IQCLUORG(ISYM,IQCLU)
                  IQB = IQCLU
                  IQA = IQCLUP
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
                  IF ( SYMUNITARY(ISYM) ) THEN
                     CNT = 'N'
                  ELSE
                     CNT = 'T'
                  END IF
C-----------------------------------------------------------------------
                  CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                       TT(1,1,IQA),M,C0,W2,M)
C
                  CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),
     &                       M,C0,W1,M)
C-----------------------------------------------------------------------
C
                  DO J = 1,NKM
                     DO I = 1,NKM
                        IF ( ABS(W1(I,J)-TT(I,J,IQB)).GT.TOLR ) THEN
                           NERR5 = NERR5 + 1
                           WRITE (6,99002) ISYM,CNT,IQA,I,J,TT(I,J,IQB),
     &                            W1(I,J),W1(I,J) - TT(I,J,IQB)
                        END IF
                     END DO
                  END DO
C
               END DO
C
            END IF
         END DO
C
C-----------------------------------------------------------------------
C
      END IF
C
      WRITE (6,99007)
C
      IF ( ALLOCATED(MG0MAT) ) DEALLOCATE (MG0MAT)
      IF ( ALLOCATED(WWA) ) DEALLOCATE (WWA)
      IF ( ALLOCATED(WWB) ) DEALLOCATE (WWB)
      IF ( ALLOCATED(IIPIV) ) DEALLOCATE (IIPIV)
      IF ( ALLOCATED(TT) ) DEALLOCATE (TT)
C
      STOP ' in <CLUTAUCHK> --- test of <CLUTAUMAT> completed'
99001 FORMAT (/,5X,'TAUQQ for IQCLU =',I4)
99002 FORMAT (/,I3,1X,A,' Q',I3,2X,2I3,3X,'TAU(q'',q'')    Q',2E13.5,/,
     &        21X,'R TAU(q,q) R+  ',2E13.5,2F20.12)
99003 FORMAT (/,1X,79('*'),/,35X,'<CLUTAUCHK>',/,20X,
     &        'checking matrix inversion by <CLUTAUMAT>',/,1X,79('*'),
     &        //,10X,'SYMMETRIZE = ',L1,/)
99004 FORMAT (/,5X,'execution time for ',A,F14.3,' secs',/)
99005 FORMAT (2X,'TROUBLE for  IQ',I3,3X,2I3,3X,2E22.14,/,32X,2E22.14,
     &        3X,A,/)
99006 FORMAT (5X,'number of deviations for run with  KIJ =',I2,/,5X,
     &        'NERROR =',I7,'           for TOL = ',E12.3,/)
99007 FORMAT (/,5X,'site diagonal blocks of TAU written to ',
     &        'files:  clutauchk_*',/)
      END
