C*==ssite_symmetrize.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SSITE_SYMMETRIZE(TXT_SS_MATRIX,XSS)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the symmetric average                                 *
C   *                                                                  *
C   *                  Y = Sum(R)  R+ X R                              *
C   *                                                                  *
C   *  for site-diagonal single-site  t,  m  or s                      *
C   *                                                                  *
C   *    restrict symmetry operations to those                         *
C   *    that leave the site postion i.e.  IQ = IQP                    *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:NXM,NKMMAX,NLM,NSPIN,WKM1,WKM2,WKM3
      USE MOD_SITES,ONLY:IQAT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SYMMETRY,ONLY:DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM
      IMPLICIT NONE
C*--SSITE_SYMMETRIZE24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SSITE_SYMMETRIZE')
      LOGICAL CHECK_SYMMETRISATION
      PARAMETER (CHECK_SYMMETRISATION=.TRUE.)
      REAL*8 TOL_SYMSUMRTR,THRESH_SYMSUMRTR
      PARAMETER (TOL_SYMSUMRTR=1D-6,THRESH_SYMSUMRTR=1D-6)
C
C Dummy arguments
C
      CHARACTER*(*) TXT_SS_MATRIX
      COMPLEX*16 XSS(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      COMPLEX*16 CSCL
      INTEGER IQ,IQP,IRELEFF,IS,ISYM,IT,LMSOFF,M,N,NSYMACCEPTED_LOCAL
      LOGICAL SAME,SAME_IT
C
C*** End of declarations rewritten by SPAG
C
      N = NXM
      M = NKMMAX
C
      IF ( IREL.LE.2 ) THEN
         IRELEFF = 1
      ELSE
         IRELEFF = 3
      END IF
C
      SAME = .TRUE.
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         LOOP_IS:DO IS = 1,NSPIN
C
            LMSOFF = NLM*(IS-1)
C
C-----------------------------------------------------------------------
C
            WKM1(1:N,1:N) = XSS(LMSOFF+1:LMSOFF+N,LMSOFF+1:LMSOFF+N,IT)
C
            WKM3(1:N,1:N) = C0
C
C-----------------------------------------------------------------------
C                      symmetrize SS-matrix
C-----------------------------------------------------------------------
            NSYMACCEPTED_LOCAL = 0
C
            LOOP_ISYM:DO ISYM = 1,NSYM
               IF ( .NOT.SYMACCEPTED(ISYM) ) CYCLE LOOP_ISYM
               IQP = IQORGQP(ISYM,IQ)
               IF ( IQ.NE.IQP ) CYCLE LOOP_ISYM
C
               NSYMACCEPTED_LOCAL = NSYMACCEPTED_LOCAL + 1
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C
               CALL ZGEMM('N',CNT,N,N,N,C1,DROT(1,1,ISYM),M,WKM1,M,C0,
     &                    WKM2,M)
C
               CALL ZGEMM('N','C',N,N,N,C1,WKM2,M,DROT(1,1,ISYM),M,C1,
     &                    WKM3,M)
C.......................................................................
C
            END DO LOOP_ISYM
C-----------------------------------------------------------------------
C                     normalize symmetrized SS-matrix
C-----------------------------------------------------------------------
C
            CSCL = 1D0/DBLE(NSYMACCEPTED_LOCAL)
C
            WKM3(1:N,1:N) = CSCL*WKM3(1:N,1:N)
C
            XSS(LMSOFF+1:LMSOFF+N,LMSOFF+1:LMSOFF+N,IT) = WKM3(1:N,1:N)
C
C-----------------------------------------------------------------------
C
C=======================================================================
            IF ( .NOT.CHECK_SYMMETRISATION ) CYCLE LOOP_IS
C=======================================================================
C
C-----------------------------------------------------------------------
C               compare original and symmetrized matrices
C-----------------------------------------------------------------------
C
            WRITE (6,99001) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,IS,
     &                      TXT_SS_MATRIX(1:LEN_TRIM(TXT_SS_MATRIX))
C
            CALL CMATCMP(N,NKMMAX,IRELEFF,'XSS - before symmetrisation',
     &                   WKM1,'YSS - after  symmetrisation',WKM3,
     &                   THRESH_SYMSUMRTR,TOL_SYMSUMRTR,SAME_IT)
C
            SAME = SAME .AND. SAME_IT
C
         END DO LOOP_IS
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      END DO LOOP_IT
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( .NOT.SAME ) CALL STOP_MESSAGE(ROUTINE,'tolerance exceeded')
C
99001 FORMAT (//,10X,'<',A,'>  calling <CMATCMP> for IT =',I4,4X,'IS =',
     &        I4,/,10X,'comparing original and symmetrized   ',A,
     &        '   matrices')
C
      END
