C*==symsumrtr.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMSUMRTR(LOCAL_SYM_ONLY,SWAP,CHECK_SYMMETRISATION,
     &                     WKSUM,T1,T2,W1,NQ,NKMQ,DSYM,IQORGQP,
     &                     SYMUNITARY,SYMACCEPTED,NSYM,NSYMACCEPTED,
     &                     NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the symmetric average                                 *
C   *                                                                  *
C   *                  T2 = Sum(R)  R+ T1 R                            *
C   *                                                                  *
C   *  for:                                                            *
C   *                                                                  *
C   *  - site-diagonal TAU-matrix after Brillouin-zone integration     *
C   *                                                                  *
C   *    SWAP = .FALSE. =>  no swap needed                             *
C   *    WKSUM    accounts for scaling from BZ-integral                *
C   *                                                                  *
C   *  - single-site m-matrix in <MSSINIT>                             *
C   *  - single-site t-matrix or m-matrix  after CPA-update            *
C   *  - site-diagonal TAU-matrix for embedded cluster                 *
C   *                                                                  *
C   *    SWAP = .TRUE.  =>  T1 and T2 data have to be swapped          *
C   *                       i.e. T1 is used as work space              *
C   *    WKSUM == 1  dummy weight                                      *
C   *                                                                  *
C   *    LOCAL_SYM_ONLY: restrict symmetry operations to those         *
C   *                    that leave the site postion i.e.  IQ = IQP    *
C   *                                                                  *
C   *    NOTE: for embedded clusters IQ, ...  correspond to IQCLU, ... *
C   *                                                                  *
C   *    NOTE: CHECK_SYMMETRISATION = .T. >> the matrices before and   *
C   *          after symmetrisation are compared with each other.      *
C   *          if the change is beyond the threshold the program stops *
C   *          this test is sensible only if the symmetrisation is     *
C   *          used to remove numerical noise ! In the case of a       *
C   *          BZ-integration over the irreducible wedge the symmetr.  *
C   *          the output TAU refers to the FULL BZ, thus the input    *
C   *          and output TAU will differ in general appreciably.      *
C   *          If input is TAU from an embedded cluster calculation    *
C   *          or TSS after calculation in SSITE, the test makes sense *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--SYMSUMRTR50
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMSUMRTR')
      REAL*8 TOL_SYMSUMRTR,THRESH_SYMSUMRTR
      PARAMETER (TOL_SYMSUMRTR=1D-8,THRESH_SYMSUMRTR=1D-8)
C
C Dummy arguments
C
      LOGICAL CHECK_SYMMETRISATION,LOCAL_SYM_ONLY,SWAP
      INTEGER NKMMAX,NQ,NQMAX,NSYM,NSYMACCEPTED
      REAL*8 WKSUM
      COMPLEX*16 DSYM(NKMMAX,NKMMAX,NSYMMAX),T1(NKMMAX,NKMMAX,NQMAX),
     &           T2(NKMMAX,NKMMAX,NQMAX),W1(NKMMAX,NKMMAX)
      INTEGER IQORGQP(NSYMMAX,NQMAX),NKMQ(NQMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      COMPLEX*16 CSCL
      LOGICAL ERROR_FLAG_NSYMACCEPTED,SAME,SAME_IQ
      INTEGER IQ,IQP,IRELEFF,ISYM,J,N,NSYMACCEPTED_LOCAL(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NSYMACCEPTED_LOCAL
C
      ALLOCATE (NSYMACCEPTED_LOCAL(NQ))
      NSYMACCEPTED_LOCAL(:) = 0
      ERROR_FLAG_NSYMACCEPTED = .FALSE.
C
C-----------------------------------------------------------------------
C                           swap T2 to T1
C-----------------------------------------------------------------------
C
      IF ( SWAP ) THEN
         DO IQ = 1,NQ
            DO J = 1,NKMQ(IQ)
               CALL ZCOPY(NKMQ(IQ),T2(1,J,IQ),1,T1(1,J,IQ),1)
            END DO
         END DO
      END IF
C
C-----------------------------------------------------------------------
C
      DO IQ = 1,NQ
         DO J = 1,NKMQ(IQ)
            CALL CINIT(NKMQ(IQ),T2(1,J,IQ))
         END DO
      END DO
C
C-----------------------------------------------------------------------
C                      symmetrize matrices
C-----------------------------------------------------------------------
      LOOP_ISYM:DO ISYM = 1,NSYM
         IF ( .NOT.SYMACCEPTED(ISYM) ) CYCLE LOOP_ISYM
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
         IF ( SYMUNITARY(ISYM) ) THEN
            CNT = 'N'
         ELSE
            CNT = 'T'
         END IF
C
         LOOP_IQ:DO IQ = 1,NQ
C
            N = NKMQ(IQ)
            IQP = IQORGQP(ISYM,IQ)
C
            IF ( IQ.EQ.IQP .OR. (.NOT.LOCAL_SYM_ONLY) ) THEN
               NSYMACCEPTED_LOCAL(IQ) = NSYMACCEPTED_LOCAL(IQ) + 1
C
C.......................................................................
               CALL ZGEMM('N',CNT,N,N,N,C1,DSYM(1,1,ISYM),NKMMAX,
     &                    T1(1,1,IQP),NKMMAX,C0,W1,NKMMAX)
C
               CALL ZGEMM('N','C',N,N,N,C1,W1,NKMMAX,DSYM(1,1,ISYM),
     &                    NKMMAX,C1,T2(1,1,IQ),NKMMAX)
C.......................................................................
            END IF
C
         END DO LOOP_IQ
C
      END DO LOOP_ISYM
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C                     normalize symmetrized matrices
C-----------------------------------------------------------------------
      DO IQ = 1,NQ
C
         IF ( (.NOT.LOCAL_SYM_ONLY) .AND. 
     &        (NSYMACCEPTED_LOCAL(IQ).NE.NSYMACCEPTED) )
     &        ERROR_FLAG_NSYMACCEPTED = .TRUE.
C
         CSCL = 1D0/DBLE(NSYMACCEPTED_LOCAL(IQ)*WKSUM)
C
         DO J = 1,NKMQ(IQ)
            CALL ZSCAL(NKMQ(IQ),CSCL,T2(1,J,IQ),1)
         END DO
C
      END DO
C-----------------------------------------------------------------------
C
      IF ( ERROR_FLAG_NSYMACCEPTED ) THEN
         WRITE (6,99002) NSYMACCEPTED,NQ,NSYMACCEPTED_LOCAL(1:NQ)
         CALL STOP_MESSAGE(ROUTINE,'NSYMACCEPTED inconsistent')
      END IF
C
C=======================================================================
      IF ( .NOT.CHECK_SYMMETRISATION ) RETURN
C=======================================================================
C
      IF ( IREL.LE.2 ) THEN
         IRELEFF = 1
      ELSE
         IRELEFF = 3
      END IF
C
C-----------------------------------------------------------------------
C               compare original and symmetrized matrices
C-----------------------------------------------------------------------
C
      SAME = .TRUE.
      DO IQ = 1,NQ
C
         WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),IQ,NKMQ(IQ),
     &                   NKMMAX,IRELEFF
C
         CALL CMATCMP(NKMQ(IQ),NKMMAX,IRELEFF,
     &                'T1 - before symmetrisation',T1(1,1,IQ),
     &                'T2 - after  symmetrisation',T2(1,1,IQ),
     &                THRESH_SYMSUMRTR,TOL_SYMSUMRTR,SAME_IQ)
C
         SAME = SAME .AND. SAME_IQ
C
      END DO
C
      IF ( .NOT.SAME ) CALL STOP_MESSAGE(ROUTINE,'tolerance exceeded')
C
99001 FORMAT (//,10X,'<',A,'>  calling <CMATCMP> for IQ =',I4,/,10X,
     &        'NKMQ(IQ) =',I4,'  NKMMAX =',I4,'  IRELEFF =',I4,/,10X,
     &        'comparing original and symmetrized matrices')
99002 FORMAT (/,10X,'expected number of symmetry operations ',
     &        '  NSYMACCEPTED =',I3,/,10X,'found for IQ = 1  ...',I4,/,
     &        (10X,20I3))
C
      END
