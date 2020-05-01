C*==strtest.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE STRTEST(ETOP)
C   ********************************************************************
C   *                                                                  *
C   *  check the calculation of the standard KKR structure constants   *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ
      USE MOD_ENERGY,ONLY:EMIN,EIMAG,NE
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_ANGMOM,ONLY:NKKR,NLM
      USE MOD_FILES,ONLY:LDATSET,DATSET
      IMPLICIT NONE
C*--STRTEST16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='STRTEST')
C
C Dummy arguments
C
      REAL*8 ETOP
C
C Local variables
C
      COMPLEX*16 ERYD,MAUX(:,:),P,TAUK(:,:)
      CHARACTER*80 FILNAM
      INTEGER I,I1,I2,I3,IA_ERR,IE,IK,ILM,ILOOP,IQ,J,JLM,JQ,LFILNAM,
     &        NKDIV,NLOOP
      REAL*8 KVEC(3),SUM_TIME,TIME,TIME0
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,MAUX
C
      ALLOCATE (TAUK(NKKR,NKKR),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MAUX')
C
      EMIN = -0.2D0
      EIMAG = 0.2D0
      NE = 2
      NKDIV = 1
C
      FILNAM = DATSET(1:LDATSET)//'_str_test.dat'
      LFILNAM = LDATSET + 13
C
      OPEN (8,FILE=FILNAM(1:LFILNAM))
C
      SUM_TIME = 0D0
      NLOOP = 100
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IE = 1,NE
C
         ERYD = DCMPLX(EMIN+(IE-1)*(ETOP-EMIN)/DBLE(NE-1),EIMAG)
         P = SQRT(ERYD)
C
C ------------------ calculate energy - dependent terms of str.constants
C
         CALL STRCC(ERYD,.FALSE.)
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
         IK = 0
         DO I1 = 0,NKDIV
            DO I2 = 0,NKDIV
               DO I3 = 0,NKDIV
C
                  IK = IK + 1
C
                  KVEC(1:3) = (I1*BBAS(1:3,1)+I2*BBAS(1:3,2)+I3*BBAS(1:3
     &                        ,3))/(NKDIV*2D0)
C
                  CALL CPU_TIME(TIME0)
C
                  DO ILOOP = 1,NLOOP
                     CALL STRSET(IK,KVEC,TAUK,MAUX,P)
                  END DO
                  CALL CPU_TIME(TIME)
                  SUM_TIME = SUM_TIME + (TIME-TIME0)/DBLE(NLOOP)
C
                  WRITE (6,99001) ERYD,KVEC
                  WRITE (8,99001) ERYD,KVEC
C
                  I = 0
                  DO IQ = 1,NQ
                     DO ILM = 1,NLM
                        I = I + 1
C
                        J = 0
                        DO JQ = 1,NQ
                           DO JLM = 1,NLM
                              J = J + 1
C
                              WRITE (8,99002) I,J,IQ,ILM,JQ,JLM,
     &                               MAUX(I,J)
C
                           END DO
                        END DO
C
                     END DO
                  END DO
C
               END DO
            END DO
         END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      WRITE (6,99003) FILNAM(1:LFILNAM)
C
      WRITE (6,99004) SUM_TIME
C
      STOP
C
99001 FORMAT (' E = ',2F10.5,'     k = ',3F10.5)
99002 FORMAT (2I4,'   I: ',2I3,'   J: ',2I3,2F22.14)
99003 FORMAT (/,10X,'data for structure constant cehck written to ',A,/)
99004 FORMAT (/,5X,'execution time ',F14.8,' secs',/)
      END
