C*==chikloops.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHIKLOOPS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,NQ,NKM,NKKR,
     &                     NKMQ,IND0Q,KMROT,DROTQ,CPATOL,NCPA,ICPA,
     &                     ICPAALG,ITCPA,ITCPAMAX,ICPACONV,DROT,IQORGQP,
     &                     SYMUNITARY,SYMACCEPTED,NOQ,ITOQ,CONC,WKTAB,
     &                     KTAB,NKTAB,NSYM,NSYMACCEPTED,NTMAX,NQMAX,
     &                     NKMMAX,NKTABMAX,TSST,MSST,TSSQ,MSSQ,TAUQ,
     &                     TAUQZ)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array                    *
C   *------------------------------------------------------------------*
C   * 01/09/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:ERYDA_EQ_ERYDB,CHIZ,NZ12,NZ12MAX,ITTA,ITTB,
     &    ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,JTTX,WTTJ,NTKTKLIN
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--CHIKLOOPS33
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPATOL
      COMPLEX*16 ERYD,P
      INTEGER ICPAALG,ICPACONV,ICPAFLAG,IPRINT,ITCPA,ITCPAMAX,KMROT,
     &        NCPA,NKKR,NKM,NKMMAX,NKTAB,NKTABMAX,NQ,NQMAX,NSYM,
     &        NSYMACCEPTED,NTMAX
      REAL*8 CONC(NTMAX),KTAB(3,NKTABMAX),WKTAB(NKTABMAX)
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),DROTQ(NKMMAX,NKMMAX,NQMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      INTEGER ICPA(NQMAX),IND0Q(NQMAX),IQORGQP(NSYMMAX,NQMAX),
     &        ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      REAL*8 CPACORR,CPAERR,CPAERRL,WK,WKSUM
      COMPLEX*16 CSCL,CSUM1,CSUMX,CWK,MAUX(:,:),
     &           SUMQ(NKMMAX,NKMMAX,NQMAX),TAUKLIN(:),TKTKQQ(:,:),
     &           W1(NKMMAX,NKMMAX),WT
      INTEGER I,I1,IA,IA_ERR,IB,IC,ID,IK,INFO,IPIV(NKKR),IQ,IQP,ISYM,J,
     &        JQ,JTT,KADI,KBCJ,KCBJ,KDAI,LV(32),N,NKMSQ,Z2
C
C*** End of declarations rewritten by SPAG
C
      DATA LV/0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,
     &     3,3,3/
C
      ALLOCATABLE MAUX,TAUKLIN,TKTKQQ
C
      ALLOCATE (MAUX(NKKR,NKKR),TKTKQQ(NTKTKLIN,NZ12))
      ALLOCATE (TAUKLIN(NKKR*NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:blochsf -> TAUKLIN'
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      TKTKQQ(1:NTKTKLIN,1:NZ12) = C0
C
C------------------------- assume the same L-expansion for every site IQ
      IF ( NZ12.GT.1 ) THEN
         NKM = NKMQ(1)
         DO IQ = 2,NQ
            IF ( NKMQ(IQ).NE.NKM ) STOP '<CHIKLOOPS>:  NKMQ(IQ)<>NKM'
         END DO
         IF ( NKM.GT.32 ) STOP '<CHIKLOOPS>:  NKM > 32'
      END IF
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ICPACONV = 0
      CPAERRL = 1.0D+6
      ITCPA = 0
 100  CONTINUE
      ITCPA = ITCPA + 1
C
      CALL CINIT(NKMMAX*NKMMAX*NQ,TAUQ)
      CALL CINIT(NKMMAX*NKMMAX*NQ,SUMQ)
C
      WKSUM = 0D0
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
         CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLIN,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLIN,MSSQ,NQMAX,NKKR,NKMMAX)
C
         CALL ZGETRF(NKKR,NKKR,TAUKLIN,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUKLIN,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C------------------------------------------------------------ store TAUQ
         WK = WKTAB(IK)
         WKSUM = WKSUM + WK
         CWK = DCMPLX(WK,0D0)
         DO IQ = 1,NQ
            I1 = IND0Q(IQ)*NKKR + IND0Q(IQ) + 1 - NKKR
            N = NKMQ(IQ)
            DO J = 1,N
               I1 = I1 + NKKR
               CALL ZAXPY(N,CWK,TAUKLIN(I1),1,SUMQ(1,J,IQ),1)
            END DO
         END DO
C--------------------------------------------------------- store TAU*TAU
         IF ( (NCPA.EQ.0) .OR. (ICPACONV.EQ.1) ) THEN
            IF ( NZ12.EQ.1 ) THEN
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     CSUM1 = CSUM1 + TAUKLIN(JTT1(JTT))
     &                       *TAUKLIN(JTT2(JTT))*WTTJ(JTT)
                  END DO
                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
               END DO
C.......................................................................
            ELSE
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  CSUMX = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
                     CSUM1 = CSUM1 + WT*TAUKLIN(JTT2(JTT))
                     CSUMX = CSUMX + WT*DCONJG(TAUKLIN(JTTX(JTT)))
                  END DO
                  TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
                  TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
               END DO
            END IF
C.......................................................................
         END IF
C-----------------------------------------------------------------------
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C-----------------------------------------------------------------------
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
            IF ( SYMUNITARY(ISYM) ) THEN
               CNT = 'N'
            ELSE
               CNT = 'T'
            END IF
C
            DO IQ = 1,NQ
C
               N = NKMQ(IQ)
               IQP = IQORGQP(ISYM,IQ)
C-----------------------------------------------------------------------
               CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),NKMMAX,
     &                    SUMQ(1,1,IQP),NKMMAX,C0,W1,NKMMAX)
C
               CALL ZGEMM('N','C',N,N,N,C1,W1,NKMMAX,DROT(1,1,ISYM),
     &                    NKMMAX,C1,TAUQ(1,1,IQ),NKMMAX)
C.......................................................................
            END DO
         END IF
      END DO
C
      CSCL = 1D0/DBLE(NSYMACCEPTED*WKSUM)
C
      DO IQ = 1,NQ
         DO J = 1,NKMQ(IQ)
            CALL ZSCAL(NKMQ(IQ),CSCL,TAUQ(1,J,IQ),1)
         END DO
      END DO
C
      IF ( (NCPA.EQ.0) .OR. (ICPACONV.EQ.1) ) THEN
         CSCL = 1D0/WKSUM
         DO Z2 = 1,NZ12
            CALL ZSCAL(NTKTKLIN,CSCL,TKTKQQ(1,Z2),1)
         END DO
      END IF
C
C cpa-loop--------------------------------------------------------------
      IF ( (NCPA.GT.0) .AND. (ICPACONV.EQ.0) ) THEN
C
         IF ( ICPAALG.EQ.1 ) THEN
            CALL CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                    NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,NTMAX,
     &                    NQMAX,NKMMAX)
C
         ELSE
            CALL CPANESBET(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NKMQ,
     &                     NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,KMROT,
     &                     DROTQ,NTMAX,NQMAX,NKMMAX)
         END IF
C
         IF ( IPRINT.EQ.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
         IF ( CPAERR.LE.CPATOL ) THEN
            ICPACONV = 1
            WRITE (6,99001) ITCPA,CPAERR,CPACORR,CPACHNG
            GOTO 100
         ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
            WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 1
         ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
            WRITE (6,99003) ITCPA
            WRITE (6,99004) CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 2
         ELSE
            CPAERRL = CPAERR
            GOTO 100
         END IF
C
      END IF
C
      IF ( IPRINT.GE.2 .AND. NCPA.GT.0 )
     &      WRITE (12,'(''E '',2F10.5,'' CPA '',I5,3E12.5)') ERYD,
     &     ICPAFLAG,CPAERR,CPACORR,CPACHNG
C
C-----------------------------------------------------------------------
C                    set up the matrix   TAUZ
C-----------------------------------------------------------------------
      DO IQ = 1,NQ
C
         DO J = 1,NKMQ(IQ)
            CALL ZCOPY(NKMQ(IQ),TAUQ(1,J,IQ),1,TAUQZ(1,J,IQ,1),1)
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C                    set up the super-matrix  CHI(K1,K2)
C  if possible: make use of symmetry relations   for   z1 = z2
C-----------------------------------------------------------------------
      NKMSQ = NKM*NKM
C
      DO Z2 = 1,NZ12
         CALL CINIT(NKMMAX**4*NQMAX**2,CHIZ(1,1,Z2))
      END DO
C
      DO I = 1,NTKTKLIN
C
         IA = ITTA(I)
         IB = ITTB(I)
         IC = ITTC(I)
         ID = ITTD(I)
         IQ = ITTQ1(I)
         JQ = ITTQ2(I)
C
         KADI = (IA-1)*NKM + ID + (IQ-1)*NKMSQ
         KBCJ = (IB-1)*NKM + IC + (JQ-1)*NKMSQ
C
         DO Z2 = 1,NZ12
C
            IF ( Z2.EQ.2 ) TKTKQQ(I,Z2) = TKTKQQ(I,Z2)*(-1)
     &           **(LV(IC)-LV(ID))
C
            CHIZ(KADI,KBCJ,Z2) = TKTKQQ(I,Z2)
C
            IF ( ERYDA_EQ_ERYDB ) THEN
               KCBJ = (IC-1)*NKM + IB + (JQ-1)*NKMSQ
               KDAI = (ID-1)*NKM + IA + (IQ-1)*NKMSQ
               CHIZ(KCBJ,KDAI,Z2) = TKTKQQ(I,Z2)
C
               IF ( IQ.EQ.JQ ) THEN
                  CHIZ(KADI,KBCJ,Z2) = CHIZ(KADI,KBCJ,Z2)
     &                                 - TAUQZ(IA,IB,IQ,1)
     &                                 *TAUQZ(IC,ID,IQ,Z2)
                  IF ( IA.NE.IC .OR. IB.NE.ID ) CHIZ(KCBJ,KDAI,Z2)
     &                 = CHIZ(KCBJ,KDAI,Z2) - TAUQZ(IA,IB,IQ,1)
     &                 *TAUQZ(IC,ID,IQ,Z2)
               END IF
            ELSE IF ( IQ.EQ.JQ ) THEN
               CHIZ(KADI,KBCJ,Z2) = CHIZ(KADI,KBCJ,Z2)
     &                              - TAUQZ(IA,IB,IQ,1)
     &                              *TAUQZ(IC,ID,IQ,Z2)
            END IF
         END DO
      END DO
C
      DEALLOCATE (MAUX,TAUKLIN)
C
99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations: ERROR ',F12.8,' CORRECTION ',F15.8,
     &        ' CHANGE ',F15.8,10('!'))
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F17.6,
     &        '    CHANGE ',F15.8)
      END
