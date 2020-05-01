C*==kloopsij.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE KLOOPSIJ(ERYD,P,NQ,NKKR,NKMQ,IND0Q,TAUQQ,DROT,
     &                    SYMUNITARY,SYMACCEPTED,MSSQ,WKTAB,KTAB,NKTAB,
     &                    NSYM,NSYMACCEPTED,NQMAX,NKMMAX,NKTABMAX)
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
C   *    NOTE  **************************************************      *
C   *                                                                  *
C   *    ALL IQ-JQ - elements of TAU are calculated                    *
C   *    NO  CPA iteration is done -- has to be done before !          *
C   *                                                                  *
C   * 01/09/2000 HE                                                    *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--KLOOPSIJ29
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      INTEGER NKKR,NKMMAX,NKTAB,NKTABMAX,NQ,NQMAX,NSYM,NSYMACCEPTED
      COMPLEX*16 DROT(NKMMAX,NKMMAX,NSYMMAX),MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQQ(NKMMAX,NKMMAX,NQMAX,NQMAX)
      INTEGER IND0Q(NQMAX),NKMQ(NQMAX)
      REAL*8 KTAB(3,NKTABMAX),WKTAB(NKTABMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      COMPLEX*16 CSCL,CWK,MAUX(:,:),SUMQQ(:,:),TAUK(:,:),W1(:,:)
      INTEGER I1,IA_ERR,IK,INFO,IPIV(:),IQ,ISYM,J,J1,JJ,JQ,N,NI,NJ
      REAL*8 WK,WKSUM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,IPIV,MAUX,SUMQQ,W1
C
      ALLOCATE (W1(NKMMAX,NKMMAX))
      ALLOCATE (TAUK(NKKR,NKKR),IPIV(NKKR))
      ALLOCATE (MAUX(NKKR,NKKR),SUMQQ(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:kloops -> SUMQQ'
C
      IF ( NSYM.GT.1 ) THEN
         WRITE (6,*) ' ################################################'
         WRITE (6,*) ' using symmetry not yet included in <KLOOPSIJ>'
         WRITE (6,*) ' calling <KLOOPSIJ> only allowed for NSYM = 1'
         WRITE (6,*) ' ################################################'
         STOP
      END IF
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL CINIT(NKMMAX*NKMMAX*NQMAX*NQMAX,TAUQQ)
      CALL CINIT(NKKR*NKKR,SUMQQ)
      CALL CINIT(NKKR*NKKR,TAUK)
C
      WKSUM = 0D0
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
         CALL STRSET(IK,KTAB(1,IK),TAUK,MAUX,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,NKMMAX)
C
         CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C----------------------------------------------------------- store TAUQQ
         WK = WKTAB(IK)
         WKSUM = WKSUM + WK
         CWK = DCMPLX(WK,0D0)
         DO J = 1,NKKR
            CALL ZAXPY(NKKR,CWK,TAUK(1,J),1,SUMQQ(1,J),1)
         END DO
C-----------------------------------------------------------------------
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      CSCL = 1D0/DBLE(NSYMACCEPTED*WKSUM)
C
      DO J = 1,NKKR
         CALL ZSCAL(NKKR,CSCL,SUMQQ(1,J),1)
      END DO
C
C=======================================================================
      IF ( NSYM.EQ.1 ) THEN
         DO JQ = 1,NQ
            J1 = IND0Q(JQ) + 1
            NJ = NKMQ(JQ)
            DO IQ = 1,NQ
               I1 = IND0Q(IQ) + 1
               NI = NKMQ(IQ)
               JJ = 0
               DO J = J1,J1 + NJ - 1
                  JJ = JJ + 1
                  CALL ZCOPY(NI,SUMQQ(I1,J),1,TAUQQ(1,JJ,IQ,JQ),1)
               END DO
            END DO
         END DO
C=======================================================================
      ELSE
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
C-----------------------------------------------------------------------
                  CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),NKMMAX,
     &                       SUMQQ(1,1),NKMMAX,C0,W1,NKMMAX)
C
                  CALL ZGEMM('N','C',N,N,N,C1,W1,NKMMAX,DROT(1,1,ISYM),
     &                       NKMMAX,C1,TAUQQ(1,1,IQ,JQ),NKMMAX)
C.......................................................................
               END DO
            END IF
         END DO
      END IF
C=======================================================================
C
      DEALLOCATE (TAUK,IPIV,MAUX,SUMQQ,W1,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:kloops -> W1'
      END
