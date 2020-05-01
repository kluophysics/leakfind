C*==symbravais.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE SYMBRAVAIS
C   ********************************************************************
C   *                                                                  *
C   *     find   BRAVAIS   requiring that all rotations connected      *
C   *     with a Bravais lattice have to map the primitive vectors     *
C   *     on themselves. If more than 1 acceptable Bravais type is     *
C   *     found that with the highest symmetry is taken                *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LATTICE,ONLY:ABAS,ADAINV,TXTBRAVAIS,BRAVAIS
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--SYMBRAVAIS15
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ABASP(3),BV(3),CF(3),MROTR(3,3,NSYMMAX),
     &       SYMEULANG(3,NSYMMAX)
      REAL*8 DDOT
      INTEGER I,IBRAACC(14),IBRAMAX,II,IINV,IPRINTLOC,IROT,IROTP,IRUN,J,
     &        NACCEPTED,NSYM,NSYMACC(14),NSYMACCEPTED,NSYMACCEPTEDMAX,
     &        NSYMCRYSYS,NSYMH
      LOGICAL SYMCRYSYS(NSYMMAX)
      CHARACTER*4 SYMSYMBL(NSYMMAX)
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
C                LOOP OVER ALL POSSIBLE BRAVAIS LATTICES
C=======================================================================
C
      NACCEPTED = 0
      NSYMACCEPTEDMAX = 0
      IBRAMAX = 0
C
      DO BRAVAIS = 1,14
C
         IPRINTLOC = IPRINT - 2
         CALL SYMINIT(IPRINTLOC,BRAVAIS,NSYM,NSYMCRYSYS,SYMCRYSYS,
     &                SYMSYMBL,SYMEULANG)
C
C=======================================================================
C             create the rotation matrices    MROTR
C=======================================================================
C
         NSYMH = NSYM/2
         DO IROT = 1,NSYMH
C
            CALL GETMROT(SYMEULANG(1,IROT),SYMEULANG(2,IROT),
     &                   SYMEULANG(3,IROT),MROTR(1,1,IROT))
C
         END DO
C
C-----------------------------------------------------------------------
C                     create matrix for inversion
C-----------------------------------------------------------------------
         IINV = NSYMH + 1
C
         CALL DGEMM('N','N',3,3,3,-1D0,MROTR(1,1,1),3,MROTR(1,1,1),3,
     &              0D0,MROTR(1,1,IINV),3)
C
C-----------------------------------------------------------------------
C                         include inversion
C-----------------------------------------------------------------------
C
         DO IROT = 2,NSYMH
            IROTP = NSYMH + IROT
C
            CALL DGEMM('N','N',3,3,3,1D0,MROTR(1,1,IROT),3,
     &                 MROTR(1,1,IINV),3,0D0,MROTR(1,1,IROTP),3)
C
         END DO
C
C-----------------------------------------------------------------------
C         check whether ALL rotations are proper
C         symmetry operations of the Bravais lattice
C-----------------------------------------------------------------------
C
         IRUN = 1
C
 50      CONTINUE
         NSYMACCEPTED = 0
         DO IROT = 1,NSYM
            IF ( SYMCRYSYS(IROT) ) THEN
               NSYMACCEPTED = NSYMACCEPTED + 1
C
               DO I = 1,3
C
                  CALL DGEMV('N',3,3,1D0,MROTR(1,1,IROT),3,ABAS(1,I),1,
     &                       0D0,ABASP,1)
C
                  DO J = 1,3
                     BV(J) = DDOT(3,ABAS(1,J),1,ABASP,1)
                  END DO
                  CALL DGEMV('N',3,3,1D0,ADAINV,3,BV,1,0D0,CF,1)
C
                  DO J = 1,3
                     IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) THEN
C
C-------------------------------------------------------------- trigonal
                        IF ( BRAVAIS.NE.10 ) GOTO 100
                        IF ( IRUN.NE.1 ) GOTO 100
                        IRUN = IRUN + 1
                        DO II = 1,3
                           SYMCRYSYS(6+II) = .FALSE.
                           SYMCRYSYS(9+II) = .TRUE.
                           SYMCRYSYS(18+II) = .FALSE.
                           SYMCRYSYS(21+II) = .TRUE.
                        END DO
                        GOTO 50
                     END IF
                  END DO
C
               END DO
C
            END IF
         END DO
C
         NACCEPTED = NACCEPTED + 1
C
         IBRAACC(NACCEPTED) = BRAVAIS
         NSYMACC(NACCEPTED) = NSYMACCEPTED
         IF ( NSYMACCEPTED.GT.NSYMACCEPTEDMAX ) THEN
            NSYMACCEPTEDMAX = NSYMACCEPTED
            IBRAMAX = BRAVAIS
         END IF
C
 100  END DO
C=======================================================================
C
      WRITE (6,99001)
C
      IF ( NACCEPTED.EQ.0 ) THEN
         WRITE (6,99002)
         STOP 'in <SYMBRAVAIS>'
      END IF
C
      WRITE (6,99003) NACCEPTED
      WRITE (6,99004) (IBRAACC(I),I=1,NACCEPTED)
      WRITE (6,99005) (NSYMACC(I),I=1,NACCEPTED)
C
      BRAVAIS = IBRAMAX
C
      WRITE (6,99006) BRAVAIS,TXTBRAVAIS(BRAVAIS),NSYMACCEPTEDMAX
C
99001 FORMAT (//,1X,79('*'),/,34X,'<SYMBRAVAIS>',/,1X,79('*'),//,10X,
     &        'find out Bravais type of the lattice ',/)
99002 FORMAT (/,' ##### TROUBLE in <SYMBRAVAIS> ',47('#'),:,/,10X,
     &        'NO BRAVAIS LATTICE WAS FOUND ',/,1X,79('#'))
99003 FORMAT (10X,I5,' possible Bravais lattices found',/)
99004 FORMAT (10X,'BRAVAIS   ',14I4)
99005 FORMAT (10X,'NSYM      ',14I4)
99006 FORMAT (/,10X,'selected Bravais lattice:',I3,' = ',A,/,10X,'with',
     &        I3,' symmetry operations',/)
      END
