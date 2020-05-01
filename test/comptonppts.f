C*==comptonppts.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE COMPTONPPTS(IPP1BOT,IPP1TOP,IPP2BOT,IPP2TOP,IPNBOT,NPN,
     &                       PPVEC1,PPVEC2,PNVEC,PEQUIV,NSYM,
     &                       SYMACCEPTED,MROTK,IPRINT,NPPTS)
C   ********************************************************************
C   *                                                                  *
C   *  scan the p-mesh used for the subroutine  <COMPTON>  and search  *
C   *  for equivalent p-points                                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--COMPTONPPTS14
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='COMPTONPPTS')
C
C Dummy arguments
C
      INTEGER IPNBOT,IPP1BOT,IPP1TOP,IPP2BOT,IPP2TOP,IPRINT,NPN,NPPTS,
     &        NSYM
      REAL*8 MROTK(3,3,NSYMMAX),PNVEC(3),PPVEC1(3),PPVEC2(3)
      INTEGER PEQUIV(NPPTS)
      LOGICAL SYMACCEPTED(NSYMMAX)
C
C Local variables
C
      REAL*8 BGINV(3,3),BGMAT(3,3),BGP(3,3),BV(3),CF(3),PSPAN(3,3)
      REAL*8 DDOT
      INTEGER I,IA_ERR,IND2(3),IPN,IPP1,IPP2,IPTAB(3),IROT,IS,J,LIN,
     &        LIN0,LIN2,NBGL(3),NBGP(3,3,NSYMMAX),NPSKIP,NPTAB
      LOGICAL PNOTDONE(:),SYMOK(NSYMMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE PNOTDONE
C
      ALLOCATE (PNOTDONE(NPPTS),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: PNOTDONE')
C
      NBGL(1) = NPN - IPNBOT + 1
      NBGL(2) = IPP2TOP - IPP2BOT + 1
      NBGL(3) = IPP1TOP - IPP1BOT + 1
      IF ( NPPTS.NE.NBGL(1)*NBGL(2)*NBGL(3) )
     &      CALL STOP_MESSAGE(ROUTINE,'NPPTS != NBGL(1)*NBGL(2)*NBGL(3)'
     &     )
C
      LIN0 = IPNBOT*NBGL(3)*NBGL(2) + IPP2BOT*NBGL(3) + IPP1BOT
      LIN = 0
      DO IPN = IPNBOT,NPN
         DO IPP2 = IPP2BOT,IPP2TOP
            DO IPP1 = IPP1BOT,IPP1TOP
               LIN = LIN + 1
               PNOTDONE(LIN) = .TRUE.
               LIN2 = IPN*NBGL(3)*NBGL(2) + IPP2*NBGL(3) + IPP1 + 1 - 
     &                LIN0
               IF ( LIN.NE.LIN2 ) WRITE (6,*) '###LIN-INDEX#######',IPN,
     &              IPP2,IPP1,LIN,LIN2
            END DO
         END DO
      END DO
      IF ( NPPTS.NE.LIN ) CALL STOP_MESSAGE(ROUTINE,'NPPTS != LIN(max)')
C
C-----------------------------------------------------------------------
C                                                       spanning vectors
      PSPAN(1:3,1) = PNVEC(1:3)
      PSPAN(1:3,2) = PPVEC2(1:3)
      PSPAN(1:3,3) = PPVEC1(1:3)
C
      DO J = 1,3
         DO I = 1,3
            BGMAT(I,J) = DDOT(3,PSPAN(1,I),1,PSPAN(1,J),1)
         END DO
      END DO
C
      CALL RINVGJ(BGINV,BGMAT,3,3)
C
C-----------------------------------------------------------------------
C     rotate the 3 spanning vectors   PSPAN  --  if it is possible
C     to express ALL rotated vectors  B(j) in terms of the old ones
C     B(i) = SUM(j) PSPAN(j) * n(j,i)  with integer coefficients n(j,i)
C     the symmetry operation creates a new point compatible with the
C     p-mesh spanned by the spanning vectors
C
      DO IROT = 1,NSYM
         IF ( SYMACCEPTED(IROT) ) THEN
            IF ( IPRINT.GT.0 ) WRITE (6,99003) IROT
            SYMOK(IROT) = .TRUE.
C
            DO I = 1,3
               CALL DGEMV('N',3,3,1D0,MROTK(1,1,IROT),3,PSPAN(1,I),1,
     &                    0D0,BGP(1,I),1)
C
               DO J = 1,3
                  BV(J) = DDOT(3,PSPAN(1,J),1,BGP(1,I),1)
               END DO
               CALL DGEMV('N',3,3,1D0,BGINV,3,BV,1,0D0,CF,1)
C
               DO J = 1,3
                  IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) SYMOK(IROT)
     &                 = .FALSE.
                  NBGP(J,I,IROT) = NINT(CF(J))
               END DO
C
            END DO
C
         ELSE
            SYMOK(IROT) = .FALSE.
         END IF
         IF ( IPRINT.GT.0 ) WRITE (6,99001) IROT,SYMOK(IROT)
      END DO
C
C-----------------------------------------------------------------------
C                           scan the p-mesh and remove equivalent points
C
      NPTAB = 0
      NPSKIP = 0
      LIN = 0
      DO IPN = IPNBOT,NPN
         IPTAB(1) = IPN
         DO IPP2 = IPP2BOT,IPP2TOP
            IPTAB(2) = IPP2
            DO IPP1 = IPP1BOT,IPP1TOP
               IPTAB(3) = IPP1
               LIN = LIN + 1
               IF ( PNOTDONE(LIN) ) THEN
                  NPTAB = NPTAB + 1
C
                  DO IROT = 1,NSYM
                     IF ( SYMOK(IROT) ) THEN
C
C               rotate p-vector  LIN
C          ROT p = SUM(i) m(i) ROT pspan(i)
C                = SUM(i) m(i) SUM(j) n(i,j) pspan(j)
C                = SUM(j) [SUM(i) m(i) n(i,j)] pspan(j)
C
                        DO J = 1,3
                           IS = 0
                           DO I = 1,3
                              IS = IS + IPTAB(I)*NBGP(J,I,IROT)
                           END DO
                           IS = MOD(IS,NBGL(J))
                           IF ( IS.LT.0 ) IS = IS + NBGL(J)
                           IND2(J) = IS
                        END DO
C
                        LIN2 = IND2(1)*NBGL(3)*NBGL(2) + IND2(2)*NBGL(3)
     &                         + IND2(3) + 1 - LIN0
                        IF ( PNOTDONE(LIN2) ) THEN
                           PNOTDONE(LIN2) = .FALSE.
                           PEQUIV(LIN2) = LIN
                           NPSKIP = NPSKIP + 1
                        END IF
C
                     END IF
                  END DO
C
                  PEQUIV(LIN) = LIN
C
               END IF
            END DO
         END DO
      END DO
C
      WRITE (6,99004) NPPTS,NPTAB,NPSKIP
C
      IF ( NPPTS.NE.NPTAB+NPSKIP ) THEN
         WRITE (6,99002) NPPTS,NPTAB + NPSKIP,NPTAB,NPSKIP
         STOP
      END IF
C
      DEALLOCATE (PNOTDONE,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
      RETURN
99001 FORMAT (10X,'IROT=',I3,'  SYMOK=',L2)
99002 FORMAT (/,10X,'trouble in <COMPTONPPTS>:',/,10X,
     &        'NPPTS != NPTAB+NPSKIP',/,10X,'NPPTS         ',I6,/,10X,
     &        'NPTAB+NPSKIP  ',I6,/,10X,'NPTAB         ',I6,/,10X,
     &        'NPSKIP        ',I6)
99003 FORMAT (10X,'rotated PSPAN  for IROT=',I3)
99004 FORMAT (/,10X,'<COMPTONPPTS>:',/,10X,
     &        'number of vectors in p-mesh ',I8,/,10X,
     &        'p-points to be treated      ',I8,/,10X,
     &        'p-points to be skipped      ',I8,/)
      END
