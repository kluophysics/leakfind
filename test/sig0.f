C*==sig0.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG0(INTEG,MTAB,MTBA,MIRRTAB,TAUTA,TAUTB,SIG0Q,
     &                SIG_MODE)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *                J_m Im TAU00(t) J_n Im TAU00(t)                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,WKM1,WKM2
      USE MOD_SITES,ONLY:NQMAX,ITOQ,NOQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_CALCMODE,ONLY:TASK
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--SIG019
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG0')
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      LOGICAL INTEG
      CHARACTER*10 SIG_MODE
      COMPLEX*16 MIRRTAB(NKMMAX,1,NKMMAX,3,3,NTMAX),
     &           MTAB(NKMMAX,NKMMAX,3,NTMAX),MTBA(NKMMAX,NKMMAX,3,NTMAX)
     &           ,SIG0Q(3,3,NQMAX),TAUTA(NKMMAX,NKMMAX,NTMAX),
     &           TAUTB(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER I,IA_ERR,IO,IQ,IT,J,M,MUE,N,NUE
      COMPLEX*16 IMTAUTA(:,:),IMTAUTB(:,:),RETAUTB(:,:),SIG0II(3,3),
     &           SIG0IR(3,3),SIG0IR_IRR(3,3),TRACE
      REAL*8 PRE_SIGOFF
      CHARACTER*5 TXTVAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IMTAUTA,IMTAUTB,RETAUTB
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      ALLOCATE (IMTAUTA(M,M),IMTAUTB(M,M),RETAUTB(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc IMTAUT')
C
      IF ( TASK.EQ.'GILBERT   ' ) THEN
         TXTVAR = 'alpha'
      ELSE
         TXTVAR = 'sigma'
      END IF
C
      IF ( SIG_MODE.EQ.'INTEGRAL  ' ) THEN
         PRE_SIGOFF = 1D0
      ELSE
         PRE_SIGOFF = 0.5D0
      END IF
C
      WRITE (6,99001)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IF ( IPRINT.GT.0 ) WRITE (6,99002) IQ
C
         CALL CINIT(3*3,SIG0II)
         CALL CINIT(3*3,SIG0IR)
         CALL CINIT(3*3,SIG0IR_IRR)
C         CALL CINIT(NKMMAX*NKMMAX*3*3*NTMAX,MIRRTAB(1,1,1,1,1,1))
C
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            N = NKMQ(IQ)
C
            DO J = 1,N
               DO I = 1,N
                  IMTAUTA(I,J) = (TAUTA(I,J,IT)-DCONJG(TAUTA(J,I,IT)))
     &                           /(2D0*CI)
                  IMTAUTB(I,J) = (TAUTB(I,J,IT)-DCONJG(TAUTB(J,I,IT)))
     &                           /(2D0*CI)
                  RETAUTB(I,J) = (TAUTB(I,J,IT)+DCONJG(TAUTB(J,I,IT)))
     &                           /2D0
C
               END DO
            END DO
C
            DO MUE = 1,3
               DO NUE = 1,3
C
C----------------------------------------------------- j_m Im G j_n Im G
C
                  IF ( .NOT.INTEG ) THEN
C
                     CALL CMATMUL(N,M,MTAB(1,1,NUE,IT),IMTAUTB,WKM1)
C
                     CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
                     CALL CMATMUL(N,M,MTBA(1,1,MUE,IT),WKM2,WKM1)
C
                     DO I = 1,NKM
                        SIG0II(MUE,NUE) = SIG0II(MUE,NUE) + CONC(IT)
     &                     *WKM1(I,I)
                     END DO
C
                  END IF
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C
                  CALL CMATMUL(N,M,MTAB(1,1,NUE,IT),RETAUTB,WKM1)
C
                  CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
                  CALL CMATMUL(N,M,MTBA(1,1,MUE,IT),WKM2,WKM1)
C
                  TRACE = C0
                  DO I = 1,NKM
                     TRACE = TRACE + WKM1(I,I)
                  END DO
C
                  CALL CMATMUL(N,M,MTAB(1,1,MUE,IT),RETAUTB,WKM1)
C
                  CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
                  CALL CMATMUL(N,M,MTBA(1,1,NUE,IT),WKM2,WKM1)
C
                  DO I = 1,NKM
                     TRACE = TRACE - WKM1(I,I)
                  END DO
C
                  SIG0IR(MUE,NUE) = SIG0IR(MUE,NUE) + CONC(IT)*CI*TRACE
C
C
C
C
C---------------------------------------------------------irregular part
C
                  CALL CMATMUL(N,M,IMTAUTA,MIRRTAB(1,1,1,NUE,MUE,IT),
     &                         WKM1)
C
                  TRACE = C0
                  DO I = 1,NKM
                     TRACE = TRACE + WKM1(I,I)
                  END DO
C
                  CALL CMATMUL(N,M,IMTAUTA,MIRRTAB(1,1,1,MUE,NUE,IT),
     &                         WKM1)
C
                  DO I = 1,NKM
                     TRACE = TRACE - WKM1(I,I)
                  END DO
C
                  SIG0IR_IRR(MUE,NUE) = SIG0IR_IRR(MUE,NUE) + CONC(IT)
     &                                  *TRACE
C
               END DO
            END DO
C
         END DO
C
C ----------------------------------- suppress small elements and sum up
C
         DO MUE = 1,3
            DO NUE = 1,3
               IF ( ABS(SIG0II(MUE,NUE)).LT.TOL ) SIG0II(MUE,NUE) = C0
C
               IF ( ABS(DIMAG(SIG0II(MUE,NUE))).LT.TOL ) SIG0II(MUE,NUE)
     &              = DREAL(SIG0II(MUE,NUE))
C
               IF ( ABS(SIG0IR(MUE,NUE)).LT.TOL ) SIG0IR(MUE,NUE) = C0
C
               IF ( ABS(DIMAG(SIG0IR(MUE,NUE))).LT.TOL ) SIG0IR(MUE,NUE)
     &              = DREAL(SIG0IR(MUE,NUE))
C
               IF ( ABS(SIG0IR_IRR(MUE,NUE)).LT.TOL )
     &              SIG0IR_IRR(MUE,NUE) = C0
C
               IF ( ABS(DIMAG(SIG0IR_IRR(MUE,NUE))).LT.TOL )
     &              SIG0IR_IRR(MUE,NUE) = DREAL(SIG0IR_IRR(MUE,NUE))
C
C- only 1/2 of SIG0IR and SIG0IR_IRR --> Crepieux formula PRB 64, 014416
C
C
               IF ( TASK.EQ.'GILBERT   ' ) THEN
C
                  SIG0Q(MUE,NUE,IQ) = SIG0II(MUE,NUE)
                  SIG0IR(MUE,NUE) = 0D0
                  SIG0IR_IRR(MUE,NUE) = 0D0
C
               ELSE
C
                  SIG0Q(MUE,NUE,IQ) = SIG0II(MUE,NUE)
     &                                + PRE_SIGOFF*SIG0IR(MUE,NUE)
     &                                - PRE_SIGOFF*SIG0IR_IRR(MUE,NUE)
C
               END IF
C
C----------------------------------------------------------------
C
            END DO
         END DO
C
C --------------------------------- write out site-diagonal term sigma_0
C
         IF ( IPRINT.GE.3 .AND. TASK.NE.'GILBERT   ' ) THEN
            WRITE (6,99003) TXTVAR
            DO MUE = 1,3
               WRITE (6,99007) (SIG0II(MUE,NUE),NUE=1,3)
            END DO
            WRITE (6,99004) TXTVAR
            DO MUE = 1,3
               WRITE (6,99007) (SIG0IR(MUE,NUE),NUE=1,3)
            END DO
            WRITE (6,99009) TXTVAR
            DO MUE = 1,3
               WRITE (6,99007) (SIG0IR_IRR(MUE,NUE),NUE=1,3)
            END DO
         END IF
C
         IF ( IPRINT.GE.2 ) THEN
            WRITE (6,99005) TXTVAR
            DO MUE = 1,3
               WRITE (6,99008) (DREAL(SIG0Q(MUE,NUE,IQ)),NUE=1,3)
            END DO
C
         END IF
C
C ----------------------------- check if imaginary part of SIG0Q is zero
C
C         IF ( .NOT.INTEG ) THEN
         DO MUE = 1,3
            DO NUE = 1,3
               IF ( ABS(DIMAG(SIG0Q(MUE,NUE,IQ))).GT.1D-5 )
     &              WRITE (6,99006) TXTVAR,DIMAG(SIG0Q(MUE,NUE,IQ)),MUE,
     &                              NUE,IQ
            END DO
         END DO
C         END IF
C
      END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (IMTAUTA,IMTAUTB,RETAUTB,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC IMTAUTA')
C
99001 FORMAT (//,1X,79('*'),/,37X,'<SIG0>',/,1X,79('*'),/)
99002 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99003 FORMAT (/,10X,'site-diagonal term ',A,'_0    j_m Im G j_n Im G',/)
99004 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',/)
99005 FORMAT (/,10X,'site-diagonal term ',A,'_0',/)
99006 FORMAT (' WARNING!! Im(',A,') =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99007 FORMAT (3(F14.6,F12.6))
99008 FORMAT (10X,3F14.6)
99009 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',/,
     &        13X,'(contribution from irregular solution)',/)
      END
