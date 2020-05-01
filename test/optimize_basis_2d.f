C*==optimize_basis_2d.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE OPTIMIZE_BASIS_2D(ABAS,QBAS,NQ,NQMAX,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  optimize the list of 2D basis vectors   ->q_i                   *
C   *  to have all basisvectors as close as possible to the z-axis     *
C   *                                                                  *
C   *  IFLAG = 0:   no changes done                                    *
C   *          1:   changes done                                       *
C   *         -1:   changes done - result questionable                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--OPTIMIZE_BASIS_2D15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      INTEGER IFLAG,NQ,NQMAX
      REAL*8 ABAS(3,3),QBAS(3,NQMAX)
C
C Local variables
C
      REAL*8 AVEC(3),DZ,DZA1,DZA2,DZMAX_NEW,DZMAX_OLD,DZMAX_Q,
     &       QBAS_OLD(:,:),QIPVEC(1:3),ZHAT(3)
      LOGICAL CHANGES_MADE
      REAL*8 DNRM2
      INTEGER I1,I123TOP,I1MIN,I2,I2MIN,IQ
C
C*** End of declarations rewritten by SPAG
C
      DATA ZHAT/0D0,0D0,1D0/
C
      ALLOCATABLE QBAS_OLD
C
      IFLAG = 0
C
      IF ( NQ.EQ.1 ) RETURN
C
      ALLOCATE (QBAS_OLD(3,NQ))
C
      QBAS_OLD(1:3,1:NQ) = QBAS(1:3,1:NQ)
C
C-----------------------------------------------------------------------
C          find max. distance from z^ for old set of basis vectors
C-----------------------------------------------------------------------
C
      DZMAX_OLD = 0D0
      DO IQ = 1,NQ
         CALL RVECXPRO(QBAS(1,IQ),ZHAT,AVEC)
         DZ = DNRM2(3,AVEC,1)
         DZMAX_OLD = MAX(DZ,DZMAX_OLD)
      END DO
C
C-----------------------------------------------------------------------
C          find range of neighboring unit cells to scan
C-----------------------------------------------------------------------
C
      DZA1 = DNRM2(3,ABAS(1,1),1)
      DZA2 = DNRM2(3,ABAS(1,2),1)
C
      I123TOP = INT(DZMAX_OLD/MIN(DZA1,DZA2)) + 2
C
C-----------------------------------------------------------------------
C
      DZMAX_NEW = 0D0
      CHANGES_MADE = .FALSE.
C----------------------------------------------------- scan all sites IQ
      DO IQ = 1,NQ
C
         DZMAX_Q = 1D+20
C
C------------------------------------------- scan neighboring unit cells
         DO I1 = I123TOP, - I123TOP, - 1
            DO I2 = I123TOP, - I123TOP, - 1
C
               QIPVEC(1:3) = QBAS(1:3,IQ) + I1*ABAS(1:3,1)
     &                       + I2*ABAS(1:3,2)
C
               CALL RVECXPRO(QIPVEC,ZHAT,AVEC)
               DZ = DNRM2(3,AVEC,1)
C
               IF ( DZ.LT.(DZMAX_Q-TOL) ) THEN
                  DZMAX_Q = DZ
                  I1MIN = I1
                  I2MIN = I2
               END IF
C
            END DO
         END DO
C-----------------------------------------------------------------------
C
         IF ( I1MIN.NE.0 .OR. I2MIN.NE.0 ) THEN
C
            QBAS(1:3,IQ) = QBAS(1:3,IQ) + I1MIN*ABAS(1:3,1)
     &                     + I2MIN*ABAS(1:3,2)
C
            CHANGES_MADE = .TRUE.
         END IF
C
         DZMAX_NEW = MAX(DZMAX_Q,DZMAX_NEW)
C
      END DO
C----------------------------------------------------- scan all sites IQ
C
      IF ( CHANGES_MADE ) THEN
C
         WRITE (6,99001)
         DO IQ = 1,NQ
            WRITE (6,99002) QBAS_OLD(1:3,IQ),QBAS(1:3,IQ)
         END DO
         WRITE (6,99003) DZMAX_OLD,DZMAX_NEW
         IFLAG = 1
C
         IF ( DZMAX_NEW.GT.(DZMAX_OLD+TOL) ) THEN
            WRITE (6,99005)
            IFLAG = -1
         END IF
C
      ELSE
C
         WRITE (6,99004) DZMAX_NEW
C
      END IF
C
99001 FORMAT (//,2(1X,79('*'),/),30X,'<OPTIMIZE_BASIS_2D>',/,
     &        2(1X,79('*'),/),//,15X,'OLD basis vectors',21X,
     &        'NEW basis vectors '/)
99002 FORMAT (2X,2(3X,'(',F10.5,',',F10.5,',',F10.5,' )'))
99003 FORMAT (/,10X,'largest distance from z-axis',/,10X,
     &        'originally:      ',F10.5,' a.u.',/,10X,
     &        'optimized:       ',F10.5,' a.u.',/)
99004 FORMAT (//,2(1X,79('*'),/),30X,'<OPTIMIZE_BASIS_2D>',/,
     &        2(1X,79('*'),/),//,10X,'no changes done on basis vectors',
     &        /,10X,'largest distance from z-axis:',F10.5,' a.u.',/)
99005 FORMAT (//,2(1X,79('#'),/),10X,'TROUBLE finding optimzed basis',/,
     &        2(1X,79('#'),/),//)
      END
