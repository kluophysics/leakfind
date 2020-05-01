C*==init_mod_lattice.f    processed by SPAG 6.70Rc at 17:35 on 30 Oct 2016
      SUBROUTINE INIT_MOD_LATTICE(IPRINT,SYSTEM_DIMENSION)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LATTICE,ONLY:ABAS,ABAS_L,ABAS_R,ABAS_I,ADAINV,ADAINV_L,
     &    ADAINV_R,ADAINV_I,ADAMAT,ADAMAT_L,ADAMAT_R,ADAMAT_I,BBAS,
     &    BDBINV,BDBMAT,VOLUC,SWS,TXTBRAVAIS,BRAVAIS,BBAS_L,SUB_SYSTEM,
     &    BBAS_R,BBAS_I,ABAS_2D,BBAS_2D,VOLUC_2D,A5BAS
      USE MOD_SITES,ONLY:NQ,QBAS,NQ_I,NQ_L,NQ_R
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT
      CHARACTER*10 SYSTEM_DIMENSION
C
C Local variables
C
      REAL*8 DDOT
      REAL*8 DET,VOLUC_I,VOLUC_L,VOLUC_R,WRK3X3(3,3)
      INTEGER I,I1,I2,IQ,J
C
C*** End of declarations rewritten by SPAG
C
      WRITE (6,99001) SYSTEM_DIMENSION(1:2)
C
C=======================================================================
C                 basis vectors in 3D - real space
C=======================================================================
C
      DO I = 1,3
         DO J = 1,3
            ADAMAT(I,J) = DDOT(3,ABAS(1,I),1,ABAS(1,J),1)
         END DO
      END DO
C
      WRK3X3(1:3,1:3) = ADAMAT(1:3,1:3)
C
      CALL RMATINV(3,3,WRK3X3,ADAINV)
C
      CALL RVECSPAT(ABAS(1,1),ABAS(1,2),ABAS(1,3),VOLUC,1)
C
      SWS = (VOLUC/(NQ*4.D0*PI/3.D0))**(1.D0/3.D0)
C
C
      IF ( IPRINT.GE.0 ) THEN
C
         WRITE (6,99002) TXTBRAVAIS(BRAVAIS)
C
         DO J = 1,3
            WRITE (6,99003) (ABAS(I,J),I=1,3)
         END DO
C
         WRITE (6,'(/,10X,''basis vectors in units of a''/)')
         DO IQ = 1,NQ
            WRITE (6,99003) (QBAS(I,IQ),I=1,3)
         END DO
C
      END IF
C
      A5BAS(1:3,1:3) = ABAS(1:3,1:3)
      A5BAS(1:3,4:5) = 0D0
C
C=======================================================================
C             basis vectors of 3D - reciprocal space
C=======================================================================
C
      DO I = 1,3
         I1 = 1 + MOD(I,3)
         I2 = 1 + MOD(I1,3)
         CALL RVECXPRO(ABAS(1,I1),ABAS(1,I2),BBAS(1,I))
      END DO
C
      CALL DSCAL(9,1D0/VOLUC,BBAS,1)
C
C-----------------------------------------------------------------------
C
C-------------------------- primitive vectors (BBAS) of reciprocal space
C
      DO I = 1,3
         I1 = 1 + MOD(I,3)
         I2 = 1 + MOD(I1,3)
         BBAS(1,I) = ABAS(2,I1)*ABAS(3,I2) - ABAS(3,I1)*ABAS(2,I2)
         BBAS(2,I) = ABAS(3,I1)*ABAS(1,I2) - ABAS(1,I1)*ABAS(3,I2)
         BBAS(3,I) = ABAS(1,I1)*ABAS(2,I2) - ABAS(2,I1)*ABAS(1,I2)
      END DO
      VOLUC = DABS(ABAS(1,1)*BBAS(1,1)+ABAS(2,1)*BBAS(2,1)+ABAS(3,1)
     &        *BBAS(3,1))
C
      DO I = 1,3
         BBAS(1,I) = BBAS(1,I)/VOLUC
         BBAS(2,I) = BBAS(2,I)/VOLUC
         BBAS(3,I) = BBAS(3,I)/VOLUC
      END DO
C
      DO I = 1,3
         DO J = 1,3
            BDBMAT(I,J) = DDOT(3,BBAS(1,I),1,BBAS(1,J),1)
         END DO
      END DO
C
      WRK3X3(1:3,1:3) = BDBMAT(1:3,1:3)
C
      CALL RMATINV(3,3,WRK3X3,BDBINV)
C
      IF ( IPRINT.GE.0 ) THEN
         WRITE (6,99004)
         DO J = 1,3
            WRITE (6,99003) (BBAS(I,J),I=1,3)
         END DO
      END IF
C
C=======================================================================
C                 basis vectors in 2D - real space
C=======================================================================
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' ) THEN
C
         NQ_I = NQ - NQ_L - NQ_R
C
         SUB_SYSTEM = 'I-ZONE'
C
C-------------------------- primitive vectors (BBAS) of reciprocal space
C                                                           3D LEFT HOST
C
         CALL RVECSPAT(ABAS_L(1,1),ABAS_L(1,2),ABAS_L(1,3),VOLUC_L,1)
C
         DO I = 1,3
            I1 = 1 + MOD(I,3)
            I2 = 1 + MOD(I1,3)
            CALL RVECXPRO(ABAS_L(1,I1),ABAS_L(1,I2),BBAS_L(1,I))
         END DO
C
         CALL DSCAL(9,1D0/VOLUC_L,BBAS_L,1)
C
         IF ( IPRINT.GE.0 ) THEN
            WRITE (6,99005) 'L-HOST'
            DO J = 1,3
               WRITE (6,99003) (BBAS_L(I,J),I=1,3)
            END DO
         END IF
C
         A5BAS(1:3,4) = ABAS_L(1:3,3)
C
C-------------------------- primitive vectors (BBAS) of reciprocal space
C                                                          3D RIGHT HOST
C
         CALL RVECSPAT(ABAS_R(1,1),ABAS_R(1,2),ABAS_R(1,3),VOLUC_R,1)
C
         DO I = 1,3
            I1 = 1 + MOD(I,3)
            I2 = 1 + MOD(I1,3)
            CALL RVECXPRO(ABAS_R(1,I1),ABAS_R(1,I2),BBAS_R(1,I))
         END DO
C
         CALL DSCAL(9,1D0/VOLUC_R,BBAS_R,1)
C
         IF ( IPRINT.GE.0 ) THEN
            WRITE (6,99005) 'R-HOST'
            DO J = 1,3
               WRITE (6,99003) (BBAS_R(I,J),I=1,3)
            END DO
         END IF
C
         A5BAS(1:3,5) = ABAS_R(1:3,3)
C
C-------------------------- primitive vectors (BBAS) of reciprocal space
C                                                              2D I-ZONE
C
         ABAS_I(1:3,1:2) = ABAS(1:3,1:2)
         ABAS_I(1:3,3) = ABAS(1:3,3) - ABAS_L(1:3,3) - ABAS_R(1:3,3)
C
         DET = ABAS_I(1,1)*ABAS_I(2,2) - ABAS_I(1,2)*ABAS_I(2,1)
C
         IF ( ABS(DET).LT.1D-8 ) STOP 
     &        ' ERROR: <INIT_MOD_LATTICE> ABAS_I are linearly dependent'
C
         BBAS_I(1,1) = ABAS_I(2,2)/DET
         BBAS_I(2,1) = -ABAS_I(1,2)/DET
         BBAS_I(1,2) = -ABAS_I(2,1)/DET
         BBAS_I(2,2) = ABAS_I(1,1)/DET
C
         VOLUC_I = ABS(DET)
C
         IF ( IPRINT.GE.0 ) THEN
            WRITE (6,99005) 'I-ZONE'
            DO J = 1,3
               WRITE (6,99003) BBAS_I(1:3,J)
            END DO
         END IF
C
C.......................................................................
C
         ABAS_2D(1:3,1:2) = ABAS(1:3,1:2)
         ABAS_2D(1:3,3) = 0D0
         ABAS_2D(3,3) = 1D0
C
         CALL RVECSPAT(ABAS_2D(1,1),ABAS_2D(1,2),ABAS_2D(1,3),VOLUC_2D,
     &                 1)
C
         DO I = 1,3
            I1 = 1 + MOD(I,3)
            I2 = 1 + MOD(I1,3)
            CALL RVECXPRO(ABAS_2D(1,I1),ABAS_2D(1,I2),BBAS_2D(1,I))
         END DO
C
         CALL DSCAL(9,1D0/VOLUC_2D,BBAS_2D,1)
C
         IF ( IPRINT.GE.0 ) THEN
            WRITE (6,99005) '2D-LATTICE'
            DO J = 1,3
               WRITE (6,99003) (BBAS_2D(I,J),I=1,3)
            END DO
         END IF
C
         IF ( ABS(VOLUC_2D-VOLUC_I).GT.1D-8 ) WRITE (*,*)
     &         '**** voluc_2d ',I,J,VOLUC_2D,VOLUC_I
         DO I = 1,3
            DO J = 1,2
               IF ( ABS(BBAS_2D(I,J)-BBAS_I(I,J)).GT.1D-8 ) WRITE (*,*)
     &               '**** BBAS_2D ',I,J,BBAS_2D(I,J),BBAS_I(I,J)
            END DO
         END DO
C-----------------------------------------------------------------------
         DO I = 1,3
            DO J = 1,3
               ADAMAT_L(I,J) = DDOT(3,ABAS_L(1,I),1,ABAS_L(1,J),1)
            END DO
         END DO
C
         WRK3X3(1:3,1:3) = ADAMAT_L(1:3,1:3)
C
         CALL RMATINV(3,3,WRK3X3,ADAINV_L)
C
C
         DO I = 1,3
            DO J = 1,3
               ADAMAT_I(I,J) = DDOT(3,ABAS_I(1,I),1,ABAS_I(1,J),1)
            END DO
         END DO
C
         WRK3X3(1:3,1:3) = ADAMAT_I(1:3,1:3)
C
         CALL RMATINV(3,3,WRK3X3,ADAINV_I)
C
         DO I = 1,3
            DO J = 1,3
               ADAMAT_R(I,J) = DDOT(3,ABAS_R(1,I),1,ABAS_R(1,J),1)
            END DO
         END DO
C
         WRK3X3(1:3,1:3) = ADAMAT_R(1:3,1:3)
C
         CALL RMATINV(3,3,WRK3X3,ADAINV_R)
C
      ELSE
C
         ADAMAT_L = 999999D0
         ADAMAT_I = 999999D0
         ADAMAT_R = 999999D0
         ADAINV_L = 999999D0
         ADAINV_I = 999999D0
         ADAINV_R = 999999D0
C
      END IF
C
99001 FORMAT (//,1X,79('*'),/,31X,'<INIT_MOD_LATTICE>',/,1X,79('*'),//,
     &        10X,'initialize the module <MOD_LATTICE>',
     &        ' according to SYSDIM = ',A,//,10X,
     &        'set up basis vectors for reciprocal space ',/)
99002 FORMAT (/,10X,'real space Bravais lattice: ',A,//,10X,
     &        'primitive vectors for real space Bravais lattice',/)
99003 FORMAT (12X,'(',F10.5,',',F10.5,',',F10.5,' )')
99004 FORMAT (/,10X,'primitive vectors of reciprocal space',
     &        ' Bravais lattice',/,10X,'       in units of 2*pi/a',/)
99005 FORMAT (/,10X,'primitive vectors of reciprocal space',/,10X,
     &        '       in units of 2*pi/a',/,10X,'       for subsystem:',
     &        A,/)
      END
