C*==getsystem.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GETSYSTEM(LSYSMAX)
C   ******************************************************************
C   *                                                                *
C   *   create name  SYSTEM  from the ordering numbers   Z  of the   *
C   *   components                                                   *
C   *                                                                *
C   ******************************************************************
C
C*** Start of declarations rewritten by SPAG
C
      USE MOD_CPA,ONLY:NCPA
      USE MOD_TYPES,ONLY:Z,CONC,NAT,NT
      USE MOD_FILES,ONLY:SYSTEM,LSYSTEM
      USE MOD_TABLES,ONLY:TAB_CHSYM
      IMPLICIT NONE
C*--GETSYSTEM17
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LSYSMAX
C
C Local variables
C
      CHARACTER*2 CS
      INTEGER I01,I10,IA_ERR,IC0,IFLAG,IT,JT,LOOP,N,NA_TAUX(:),Z_TAUX(:)
      CHARACTER*4 STR4
      REAL*8 X_TAUX(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE NA_TAUX,X_TAUX,Z_TAUX
C
      ALLOCATE (NA_TAUX(NT),X_TAUX(NT),Z_TAUX(NT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:getsystem -> Z_TAUX'
C
      DO IT = 1,NT
         NA_TAUX(IT) = NAT(IT)
         X_TAUX(IT) = CONC(IT)
         Z_TAUX(IT) = Z(IT)
      END DO
C
C------------------------------------ combine types with the same name
C
      IF ( NCPA.EQ.0 ) THEN
         DO IT = 1,NT
            DO JT = IT + 1,NT
               IF ( NA_TAUX(JT).NE.0 ) THEN
                  IF ( Z_TAUX(IT).EQ.Z_TAUX(JT) ) THEN
                     NA_TAUX(IT) = NA_TAUX(IT) + NA_TAUX(JT)
                     NA_TAUX(JT) = 0
                  END IF
               END IF
            END DO
         END DO
      END IF
C
C--------------------------- remove common factor in the stoichiometry
      N = 7
      DO LOOP = 1,5
         N = N - 1
         IFLAG = 1
         DO IT = 1,NT
            IF ( NA_TAUX(IT).NE.(N*(NA_TAUX(IT)/N)) ) IFLAG = 0
         END DO
         IF ( IFLAG.EQ.1 ) THEN
            DO IT = 1,NT
               NA_TAUX(IT) = NA_TAUX(IT)/N
            END DO
         END IF
      END DO
C
C-------------------------------------------------- create system name
      LSYSTEM = 0
      IC0 = ICHAR('0')
C
      DO IT = 1,NT
         IF ( NA_TAUX(IT).GT.0 ) THEN
            CS = TAB_CHSYM(Z_TAUX(IT))
            IF ( LSYSTEM.GE.LSYSMAX-2 ) GOTO 100
            IF ( IT.EQ.1 ) THEN
               SYSTEM = CS
            ELSE
               SYSTEM = SYSTEM(1:LSYSTEM)//CS
            END IF
            LSYSTEM = LSYSTEM + 2 - INDEX(CS,' ')/2
C
            IF ( NA_TAUX(IT).GT.1 ) THEN
               I10 = NA_TAUX(IT)/10 + IC0
               I01 = NA_TAUX(IT) - 10*(NA_TAUX(IT)/10) + IC0
               IF ( NA_TAUX(IT).GE.10 ) THEN
                  IF ( LSYSTEM.GE.LSYSMAX-2 ) GOTO 100
                  SYSTEM = SYSTEM(1:LSYSTEM)//CHAR(I10)
                  LSYSTEM = LSYSTEM + 1
               END IF
               IF ( LSYSTEM.GE.LSYSMAX-2 ) GOTO 100
               SYSTEM = SYSTEM(1:LSYSTEM)//CHAR(I01)
               LSYSTEM = LSYSTEM + 1
            END IF
C
            IF ( ABS(X_TAUX(IT)-1.0D0).GT.0.001D0 ) THEN
               IF ( NA_TAUX(IT).GT.1 ) THEN
                  IF ( LSYSTEM.GE.LSYSMAX-2 ) GOTO 100
                  SYSTEM = SYSTEM(1:LSYSTEM)//'_'
                  LSYSTEM = LSYSTEM + 1
               END IF
               OPEN (99,STATUS='SCRATCH')
               WRITE (99,'(F5.2)') X_TAUX(IT)
               REWIND 99
               READ (99,'(A)') STR4
               CLOSE (99)
               IF ( LSYSTEM.GE.LSYSMAX-4 ) GOTO 100
               SYSTEM = SYSTEM(1:LSYSTEM)//STR4
               IF ( STR4(4:4).EQ.'0' ) THEN
                  LSYSTEM = LSYSTEM + 3
               ELSE
                  LSYSTEM = LSYSTEM + 4
               END IF
            END IF
         END IF
      END DO
      IF ( LSYSTEM.GT.LSYSMAX ) STOP '<GETSYSTEM> LSYSTEM > LSYSMAX'
C
      DEALLOCATE (NA_TAUX,X_TAUX,Z_TAUX,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:getsystem -> Z_TAUX'
      RETURN
 100  CONTINUE
      WRITE (6,*) ' WARNING from  <GETSYSTEM>:  LSYSTEM > LSYSMAX'
      END
