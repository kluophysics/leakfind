C*==init_mod_strlag.f    processed by SPAG 6.55Rc at 09:09 on  9 Jul 2007
      SUBROUTINE INIT_MOD_STRLAG(GFEP,NGFEP_ARG,NGFEPMAX,IGRID,EMIN,
     &                           EMAX,EIMAG,ETAB,NETAB,NEMAX,KTAB,
     &                           NKTAB_ARG,NKTABMAX)
C   ********************************************************************
C   *                                                                  *
C   *  initialize the Lagrangian interpolation of the structure        *
C   *  constants for a given                                           *
C   *                                                                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_STR,ONLY:NGRL,G1,G2,G3,EXPGNQ,NQQP_STR,ETA,NLLMMMAX
      USE MOD_STRLAG,ONLY:EXPGNQ_RED,NSTRLAG,ESTRLAG,WSTRLAG,
     &                    NGFEP,DLMLAG
      USE MOD_KSPACE,ONLY:NKTAB
      USE MOD_CONSTANTS,ONLY:PI

      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EIMAG,EMAX,EMIN
      INTEGER IGRID,NEMAX,NGFEPMAX,NGFEP_ARG,NKTABMAX,NKTAB_ARG
      COMPLEX*16 ETAB(NEMAX,2)
      REAL*8 GFEP(3,NGFEPMAX),KTAB(3,NKTABMAX)
      INTEGER NETAB(2)
C
C Local variables
C
      COMPLEX*16 CSWAP,DE,F
      DOUBLE PRECISION DBLE
      REAL*8 GX,GY,GZ,IME,PHI,REE
      INTEGER I,I1,I2,I3,IE,IGFEP,IK,IQQP,IX,J,KGFEP(:),N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KGFEP
C
      ALLOCATE (KGFEP(NGRL))
C
      NGFEP = NGFEP_ARG
      IF( NKTAB .NE. NKTAB_ARG ) THEN 
         WRITE(6,*) 'INCONSISTENCY in <INIT_MOD_STRLAG>'
         WRITE(6,*) 'NKTAB     = ',NKTAB
         WRITE(6,*) 'NKTAB_ARG = ',NKTAB_ARG
         STOP
      END IF 
C
C  =====================================================================
C  indicate those rec. lattice vectors that give rise to free el. poles
C  =====================================================================
C
      DO I = 1,NGRL
C
         KGFEP(I) = 0
C
         GX = G1(I)*BBAS(1,1) + G2(I)*BBAS(1,2) + G3(I)*BBAS(1,3)
         GY = G1(I)*BBAS(2,1) + G2(I)*BBAS(2,2) + G3(I)*BBAS(2,3)
         GZ = G1(I)*BBAS(3,1) + G2(I)*BBAS(3,2) + G3(I)*BBAS(3,3)
C
         DO IGFEP = 1,NGFEP
C
            IF ( ABS(GX-GFEP(1,IGFEP)).LT.1D-8 .AND. 
     &           ABS(GY-GFEP(2,IGFEP)).LT.1D-8 .AND. 
     &           ABS(GZ-GFEP(3,IGFEP)).LT.1D-8 ) KGFEP(I) = 1
C
         END DO
C
      END DO
C
C  =====================================================================
C         sort rec. lattice vectors to have the poles at the begining
C  =====================================================================
C
      DO I = 1,NGRL
C
         IF ( KGFEP(I).NE.1 ) THEN
            DO J = I + 1,NGRL
               IF ( KGFEP(J).EQ.1 ) THEN
                  I1 = G1(I)
                  G1(I) = G1(J)
                  G1(J) = I1
                  I2 = G2(I)
                  G2(I) = G2(J)
                  G2(J) = I2
                  I3 = G3(I)
                  G3(I) = G3(J)
                  G3(J) = I3
                  DO IQQP = 1,NQQP_STR
                     CSWAP = EXPGNQ(I,IQQP)
                     EXPGNQ(I,IQQP) = EXPGNQ(J,IQQP)
                     EXPGNQ(J,IQQP) = CSWAP
                  END DO
                  KGFEP(I) = 1
                  KGFEP(J) = 0
                  GOTO 100
               END IF
            END DO
         END IF
 100  END DO
C
      WRITE (6,99001)
      DO I = 1,NGFEP
         WRITE (6,99002) I,(GFEP(IX,I),IX=1,3),G1(I),G2(I),G3(I)
      END DO
      WRITE (6,*) ' '
C
C
C  =====================================================================
C      remove factor exp(-K_n^2/eta) EXPGNQ and store in EXPGNQ_RED
C  =====================================================================
C
      ALLOCATE (EXPGNQ_RED(NGRL,NQQP_STR))
C
      DO N = 1,NGRL
C
         GX = G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)*BBAS(1,3)
         GY = G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)*BBAS(2,3)
         GZ = G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)*BBAS(3,3)
C
         F = EXP(-(GX**2+GY**2+GZ**2)/ETA)
C
         DO IQQP = 1,NQQP_STR
            EXPGNQ_RED(N,IQQP) = EXPGNQ(N,IQQP)/F
         END DO
      END DO
C
C  =====================================================================
C                  prepare array for interpolation
C  =====================================================================
C
      NSTRLAG = MIN(5,NETAB(1))
      NSTRLAG = 11
C
      ALLOCATE (ESTRLAG(NSTRLAG),WSTRLAG(NSTRLAG))
C
      IF ( NSTRLAG.EQ.NETAB(1) ) THEN
C
         DO IE = 1,NSTRLAG
            ESTRLAG(IE) = ETAB(IE,1)
         END DO
C
      ELSE IF ( IGRID.EQ.3 ) THEN
C
         DE = (EMAX-EMIN)/DBLE(NSTRLAG-1)
         DO IE = 1,NSTRLAG
            ESTRLAG(IE) = DCMPLX(EMIN,EIMAG) + DE*(IE-1)
         END DO
C
      ELSE IF ( IGRID.EQ.5 ) THEN
C
         DO IE = 1,NSTRLAG
            PHI = -PI*COS(PI*DBLE(2*IE-1)/DBLE(2*NSTRLAG))/2D0
            IME = COS(PHI)*(EMAX-EMIN)/2D0
            REE = SIN(PHI)*(EMAX-EMIN)/2D0 + (EMAX+EMIN)/2D0
            ESTRLAG(IE) = DCMPLX(REE,IME)
            WRITE (*,'(A,I3,f18.12,2X,2f18.12)') '###',IE,PHI,
     &             ESTRLAG(IE)
         END DO
C
      ELSE
         WRITE (*,*) 'IGRID = ',IGRID
         STOP 'in <STRINITLAG>'
C
      END IF
C
      ALLOCATE (DLMLAG(NLLMMMAX,NQQP_STR,NSTRLAG,NKTAB))
C
      DO IE = 1,NSTRLAG
C
         CALL STRCC(ESTRLAG(IE),.FALSE.)
C
         DO IK = 1,NKTAB
C
            CALL STRBBDDRED(DLMLAG(1,1,IE,IK),KTAB(1,IK),KTAB(2,IK),
     &                      KTAB(3,IK))
C
         END DO
C
      END DO
C
      WRITE (*,99003) NKTAB,NSTRLAG
C
C-----------------------------------------------------------------------
99001 FORMAT (/,1X,79('*'),/,34X,'<STRINITLAG>',/,1X,79('*'),//,10X,
     &    'reciprocal lattice vectors with possible free electron poles'
     &    ,/)
99002 FORMAT (10X,I3,3F10.4,5X,3I4)
99003 FORMAT (/,10X,'interpolation of structure constants initiated for'
     &        ,/,10X,'NKTAB =',I6,' k-vectors',5X,'NSTRLAG =',I3,
     &        ' energy grid points',/)
      END
