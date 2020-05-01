C*==sig_drv_kmesh.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG_DRV_KMESH(IEA,NEA,SIG_MODE)
      USE MOD_CALCMODE,ONLY:ITEST,MOL,KMROT
      USE MOD_KSPACE,ONLY:NKTABMAX,NKTABINP
      USE MOD_SYMMETRY,ONLY:SNTAUUVMAX,NELMTMAX
      USE MOD_SIG,ONLY:NKTABSIG
      USE MOD_CPA,ONLY:NCPA
      IMPLICIT NONE
C*--SIG_DRV_KMESH9
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IEA,NEA
      CHARACTER*10 SIG_MODE
C
C Local variables
C
      INTEGER DIV,DIV_VIC_EF
      LOGICAL INITIALIZE
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
C##------------------------------------------------------- old behaviour
C##      nktabsig = nktabinp
C##      CALL INIT_MOD_KSPACE(.TRUE.,MOL,KMROT,ITEST,NKTABMAX,NELMTMAX,
C##  &                        SNTAUUVMAX)
C##      return
C
      WRITE (6,99005) 'calling SIG_DRV_KMESH'
C
      NKTABSIG = NKTABINP
C
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C       adjustment of k-mesh with imaginary part of energy E_a
C       in case of DC conductivity via Kubo-bastin formula
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      IF ( SIG_MODE.NE.'OPTICS    ' ) THEN
C
         DIV = 1000
C
         IF ( NKTABINP.GE.DIV ) THEN
C
C----------------------differentiate between CPA and non-CPA calculation
            IF ( NCPA.EQ.0 ) THEN
               DIV_VIC_EF = 10
            ELSE
               DIV_VIC_EF = 100
            END IF
C
            IF ( INITIALIZE ) THEN
               IF ( NCPA.EQ.1 ) WRITE (6,99001)
               WRITE (6,99002) NKTABINP
               IF ( NEA.GE.2 ) WRITE (6,99003) NKTABINP/DIV_VIC_EF
               IF ( NEA.GT.2 ) WRITE (6,99004) NKTABINP/DIV
               INITIALIZE = .FALSE.
            END IF
C
            IF ( IEA.LT.(NEA-1) ) THEN
               NKTABSIG = NKTABINP/DIV
            ELSE IF ( IEA.EQ.(NEA-1) ) THEN
               NKTABSIG = NKTABINP/DIV_VIC_EF
            ELSE
               NKTABSIG = NKTABINP
            END IF
C
         ELSE
            WRITE (6,99006) NKTABINP,DIV
         END IF
C
      END IF
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      CALL INIT_MOD_KSPACE(.TRUE.,MOL,KMROT,ITEST,NKTABMAX,NELMTMAX,
     &                     SNTAUUVMAX)
C
99001 FORMAT (/,10X,'Doing CPA calculation ',/)
99002 FORMAT (10X,'NKTAB for point at      EF is',i15)
99003 FORMAT (10X,'NKTAB for point next to EF is',i15)
99004 FORMAT (10X,'NKTAB for other points     is',i15)
99005 FORMAT (/,79('*'),/,10X,A,/,79('*'),/)
99006 FORMAT (/,10X,'NKTAB',i8,' lesser than DIV:',i8,/,10X,
     &        ' ---> no adjustment of kpoints done')
      END
