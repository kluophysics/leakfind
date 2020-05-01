C*==signlc0.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLC0(INTEG,MTAB,MTBA,TAUTA,TAUTB,CONC,NOQ,ITOQ,
     &                   SIG0Q,IPRINT,NCFG,PCFG,IQCPA,NQNLCPA,W1HAT,
     &                   MSSTA,MSSTB,OMEGAHATA,OMEGAHATB,NDIMCLU,
     &                   IND0QCLU,NKMQ,NKMMAX,NTMAX,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *                J_m Im TAU00(t) J_n Im TAU00(t)                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:TASK
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_SITES,ONLY:IQBOT_CHI,IQTOP_CHI
      IMPLICIT NONE
C*--SIGNLC018
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLC0')
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      LOGICAL INTEG
      INTEGER IPRINT,IQCPA,NCFG,NDIMCLU,NKMMAX,NQMAX,NQNLCPA,NTMAX
      REAL*8 CONC(NTMAX),PCFG(NCFG)
      INTEGER IND0QCLU(NQNLCPA),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)
      COMPLEX*16 MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           MTAB(NKMMAX,NKMMAX,3,NTMAX),MTBA(NKMMAX,NKMMAX,3,NTMAX)
     &           ,OMEGAHATA(NDIMCLU,NDIMCLU),OMEGAHATB(NDIMCLU,NDIMCLU),
     &           SIG0Q(3,3,NQMAX),TAUTA(NKMMAX,NKMMAX,NTMAX),
     &           TAUTB(NKMMAX,NKMMAX,NTMAX),W1HAT(NDIMCLU,NDIMCLU)
C
C Local variables
C
      INTEGER I,I1_1,I1_N,I2_1,I2_N,IA_ERR,ICFG,IO,IOCC1,IOCC2,IQ,
     &        IQCLU1,IQCLU2,IQCLUREP,IT,IT1,IT2,J,MUE,N,NUE
      COMPLEX*16 IMTAUTA(:,:),IMTAUTB(:,:),RETAUTB(:,:),SIG0II(3,3),
     &           SIG0IR(3,3),TAU12A(:,:),TAU12B(:,:),TAU21A(:,:),
     &           TAU21B(:,:),TAUIMPA(:,:),TAUIMPB(:,:)
      INTEGER NLCPACONF
      CHARACTER*5 TXTVAR
      REAL*8 WGT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAU12A,TAU12B,TAU21A,TAU21B
      ALLOCATABLE IMTAUTA,IMTAUTB,RETAUTB,TAUIMPA,TAUIMPB
C
      ALLOCATE (TAU12A(NKMMAX,NKMMAX),TAU21A(NKMMAX,NKMMAX))
      ALLOCATE (TAU12B(NKMMAX,NKMMAX),TAU21B(NKMMAX,NKMMAX))
      ALLOCATE (RETAUTB(NKMMAX,NKMMAX))
      ALLOCATE (IMTAUTA(NKMMAX,NKMMAX),IMTAUTB(NKMMAX,NKMMAX))
      ALLOCATE (TAUIMPB(NDIMCLU,NDIMCLU))
      ALLOCATE (TAUIMPA(NDIMCLU,NDIMCLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate IMTAUT')
C
      IF ( TASK.EQ.'GILBERT   ' ) THEN
         TXTVAR = 'alpha'
      ELSE
         TXTVAR = 'sigma'
      END IF
C
      WRITE (6,99001)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ = IQBOT_CHI,IQTOP_CHI
         WRITE (6,99002) IQ
C
         CALL CINIT(3*3,SIG0II)
         CALL CINIT(3*3,SIG0IR)
C
         N = NKMQ(IQ)
C
C-----------------------------------------------------------------------
         IF ( IQ.EQ.IQCPA ) THEN
C
            DO ICFG = 1,NCFG
C
               CALL NLCPAPROJ(ICFG,IQCPA,MSSTA,OMEGAHATA,TAUIMPA,W1HAT,
     &                        NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,
     &                        NTMAX,NKMMAX)
C
               CALL NLCPAPROJ(ICFG,IQCPA,MSSTB,OMEGAHATB,TAUIMPB,W1HAT,
     &                        NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,
     &                        NTMAX,NKMMAX)
C
               DO IQCLUREP = 1,NQNLCPA
C
                  IQCLU1 = IQCLUREP
C
                  IOCC1 = NLCPACONF(ICFG,NQNLCPA,IQCLU1,1,2)
                  IT1 = ITOQ(IOCC1,IQCPA)
                  I1_1 = IND0QCLU(IQCLU1) + 1
                  I1_N = IND0QCLU(IQCLU1) + N
C
                  DO IQCLU2 = 1,NQNLCPA
C
                     IOCC2 = NLCPACONF(ICFG,NQNLCPA,IQCLU2,1,2)
                     IT2 = ITOQ(IOCC2,IQCPA)
                     I2_1 = IND0QCLU(IQCLU2) + 1
                     I2_N = IND0QCLU(IQCLU2) + N
                     TAU12A(1:N,1:N) = TAUIMPA(I1_1:I1_N,I2_1:I2_N)
                     TAU21A(1:N,1:N) = TAUIMPA(I2_1:I2_N,I1_1:I1_N)
                     TAU12B(1:N,1:N) = TAUIMPB(I1_1:I1_N,I2_1:I2_N)
                     TAU21B(1:N,1:N) = TAUIMPB(I2_1:I2_N,I1_1:I1_N)
C
                     DO J = 1,N
                        DO I = 1,N
                           IMTAUTA(I,J)
     &                        = (TAU12A(I,J)-DCONJG(TAU21A(J,I)))
     &                        /(2D0*CI)
                           IMTAUTB(I,J)
     &                        = (TAU12B(I,J)-DCONJG(TAU21B(J,I)))
     &                        /(2D0*CI)
                           RETAUTB(I,J)
     &                        = (TAU12B(I,J)+DCONJG(TAU21B(J,I)))/2D0
                        END DO
                     END DO
C
                     WGT = PCFG(ICFG)/DBLE(NQNLCPA)
C
                     CALL SIGNLC0AUX(INTEG,IT1,IT2,N,SIG0II,SIG0IR,
     &                               IMTAUTA,IMTAUTB,RETAUTB,MTAB,MTBA,
     &                               WGT,NKMMAX,NTMAX)
C
                  END DO
C
               END DO
C
            END DO
C
         ELSE
C-----------------------------------------------------------------------
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               DO J = 1,N
                  DO I = 1,N
                     IMTAUTA(I,J) = (TAUTA(I,J,IT)-DCONJG(TAUTA(J,I,IT))
     &                              )/(2D0*CI)
                     IMTAUTB(I,J) = (TAUTB(I,J,IT)-DCONJG(TAUTB(J,I,IT))
     &                              )/(2D0*CI)
                     RETAUTB(I,J) = (TAUTB(I,J,IT)+DCONJG(TAUTB(J,I,IT))
     &                              )/2D0
                  END DO
               END DO
C
               WGT = CONC(IT)
C
               CALL SIGNLC0AUX(INTEG,IT,IT,N,SIG0II,SIG0IR,IMTAUTA,
     &                         IMTAUTB,RETAUTB,MTAB,MTBA,WGT,NKMMAX,
     &                         NTMAX)
C
            END DO
C
         END IF
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
C
               IF ( TASK.EQ.'GILBERT   ' ) THEN
C
                  SIG0Q(MUE,NUE,IQ) = SIG0II(MUE,NUE)
                  SIG0IR(MUE,NUE) = 0D0
C
               ELSE
C
                  SIG0Q(MUE,NUE,IQ) = SIG0II(MUE,NUE) + SIG0IR(MUE,NUE)
C
               END IF
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
         IF ( .NOT.INTEG ) THEN
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(DIMAG(SIG0Q(MUE,NUE,IQ))).GT.1D-5 )
     &                 WRITE (6,99006) TXTVAR,DIMAG(SIG0Q(MUE,NUE,IQ)),
     &                                 MUE,NUE,IQ
               END DO
            END DO
         END IF
C
      END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (IMTAUTA,IMTAUTB,RETAUTB)
      DEALLOCATE (TAU12A,TAU12B,TAU21A,TAU21B,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC TAU12A')
C
99001 FORMAT (//,1X,79('*'),/,36X,'<SIGNLC0>',/,1X,79('*'),/)
99002 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99003 FORMAT (/,10X,'site-diagonal term ',A,'_0    j_m Im G j_n Im G',/)
99004 FORMAT (/,10X,'site-diagonal term ',A,'_0  i j_m Im G j_n Re G',/)
99005 FORMAT (/,10X,'site-diagonal term ',A,'_0',/)
99006 FORMAT (' WARNING!! Im(',A,') =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99007 FORMAT (3(F14.6,F12.6))
99008 FORMAT (10X,3F14.6)
      END
C*==signlc0aux.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLC0AUX(INTEG,IT1,IT2,N,SIG0II,SIG0IR,IMTAUTA,
     &                      IMTAUTB,RETAUTB,MTAB,MTBA,WGT,NKMMAX,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  auxilary routine to <SIGNLC0>                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:WKM1,WKM2
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--SIGNLC0AUX271
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL INTEG
      INTEGER IT1,IT2,N,NKMMAX,NTMAX
      REAL*8 WGT
      COMPLEX*16 IMTAUTA(NKMMAX,NKMMAX),IMTAUTB(NKMMAX,NKMMAX),
     &           MTAB(NKMMAX,NKMMAX,3,NTMAX),MTBA(NKMMAX,NKMMAX,3,NTMAX)
     &           ,RETAUTB(NKMMAX,NKMMAX),SIG0II(3,3),SIG0IR(3,3)
C
C Local variables
C
      INTEGER I,M,MUE,NUE
      COMPLEX*16 TRACE
C
C*** End of declarations rewritten by SPAG
C
      M = NKMMAX
C
      DO MUE = 1,3
         DO NUE = 1,3
C
C----------------------------------------------------- j_m Im G j_n Im G
C
            IF ( .NOT.INTEG ) THEN
C
               CALL CMATMUL(N,M,MTAB(1,1,NUE,IT1),IMTAUTB,WKM1)
C
               CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
               CALL CMATMUL(N,M,MTBA(1,1,MUE,IT2),WKM2,WKM1)
C
               DO I = 1,N
                  SIG0II(MUE,NUE) = SIG0II(MUE,NUE) + WGT*WKM1(I,I)
               END DO
C
            END IF
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C
            CALL CMATMUL(N,M,MTAB(1,1,NUE,IT1),RETAUTB,WKM1)
C
            CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
            CALL CMATMUL(N,M,MTBA(1,1,MUE,IT2),WKM2,WKM1)
C
            TRACE = C0
            DO I = 1,N
               TRACE = TRACE + WKM1(I,I)
            END DO
C
C
            CALL CMATMUL(N,M,MTAB(1,1,MUE,IT1),RETAUTB,WKM1)
C
            CALL CMATMUL(N,M,IMTAUTA,WKM1,WKM2)
C
            CALL CMATMUL(N,M,MTBA(1,1,NUE,IT2),WKM2,WKM1)
C
            DO I = 1,N
               TRACE = TRACE - WKM1(I,I)
            END DO
C
            SIG0IR(MUE,NUE) = SIG0IR(MUE,NUE) + WGT*CI*TRACE
C
         END DO
      END DO
C
      END
