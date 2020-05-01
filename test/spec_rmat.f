C*==zmatxps.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C     contains the subroutines
C     zmatxps  calzmat  detg  romat  cind  cpcorewff  simpsnn
C     romatm
C
      SUBROUTINE ZMATXPS(CWFSTATES,CWFG,CWFX,CWFY,CWF,RDIP1M,OMHAR,
     &                   ZMAT0,ZMATD,ALPHA,IAN,IT,AMAT1X,AMAT1Y,AMAT1V,
     &                   AMAT1,AMAT2X,AMAT2Y,AMAT2V,IBLOCH,DOS,IO,
     &                   AMAT1T,AMAT2T)
C
C     purpose : prepare data for the call to subroutine calzmat
C
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,XMAXE,NFULLPOT,MAXCORE,
     &    MAXCSTATES,EPS12,CZERO,NTPHOMAX,NVFTPHOMAX
      USE MOD_SPEC_WAVE,ONLY:WFFM,WFFXM,WFFYM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,OMHAR
      INTEGER IAN,IBLOCH,IO,IT
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),CWF(RSTEP,MAXCORE,2),
     &           DOS(MQD,MQD,NATLM,LAYSM),
     &           RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           ZMAT0(MQD,MQD,LAYSM,NATLM),ZMATD(MQD,MQD,LAYSM,NATLM)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        CWFSTATES(NATLM,LAYSM),CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),
     &        CWFY(MAXCORE)
      REAL*8 CWFG(MAXCSTATES,NATLM,LAYSM)
C
C Local variables
C
      INTEGER CONTRIB,I,J
      COMPLEX*16 ZMAT_ATOM(:,:),ZMAT_DOS(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZMAT_DOS,ZMAT_ATOM
      ALLOCATE (ZMAT_DOS(MQD,MQD),ZMAT_ATOM(MQD,MQD))
C
C*** End of declarations rewritten by SPAG
C
C     check if atom gives a contribution to photocurrent:
C     ===================================================
      CONTRIB = 0
      DO I = 1,CWFSTATES(IAN,IT)
         IF ( ABS(CWFG(I,IAN,IT)).GT.EPS12 ) CONTRIB = 1
C         write (*,*) 'contrib',ian,it,i,contrib,cwfg(i,ian,it)
      END DO
C
      DO I = 1,MQD
         DO J = 1,MQD
            ZMAT0(I,J,IT,IAN) = CZERO
            ZMATD(I,J,IT,IAN) = CZERO
         END DO
      END DO
      IF ( CONTRIB.NE.0 ) THEN
C
C         calculate z-matrix:
C         calculate zmatrix for atom it,ian and specified polarization:
C         ==============================================================
         CALL CALZMAT(IT,IAN,CWFSTATES,CWFG,CWFX,CWFY,CWF,AMAT1,RDIP1M,
     &                ZMAT_ATOM,ZMAT_DOS,AMAT1X,AMAT1Y,AMAT1V,AMAT2X,
     &                AMAT2Y,AMAT2V,ALPHA,WFFM,WFFXM,WFFYM,OMHAR,IBLOCH,
     &                DOS,IO,AMAT1T,AMAT2T,NRMAX,NTPHOMAX,NVFTPHOMAX)
C
C         set up z-matrix:
C         ================
         DO I = 1,MQD
            DO J = 1,MQD
               ZMAT0(I,J,IT,IAN) = ZMAT0(I,J,IT,IAN) + ZMAT_ATOM(I,J)
               IF ( IBLOCH.NE.2 ) ZMATD(I,J,IT,IAN) = ZMATD(I,J,IT,IAN)
     &              + ZMAT_DOS(I,J)
            END DO
         END DO
      END IF
C
      END
C*==calzmat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CALZMAT(LAY,ATOM,CWFSTATES,CWFG,CWFX,CWFY,CWF,AMAT1,
     &                   RDIP1M,ZMAT_A,ZMAT_D,AMAT1X,AMAT1Y,AMAT1V,
     &                   AMAT2X,AMAT2Y,AMAT2V,ALPHA,WFFM,WFFXM,WFFYM,
     &                   OMHAR,IBLOCH,DOS,IO,AMAT1T,AMAT2T,NRMAX,
     &                   NTPHOMAX,NVFTPHOMAX)
C
C     purpose       : calculate the zmatrix
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MAXWF,NSTATES,XMAXE,
     &    NFULLPOT,MAXCORE,MAXCSTATES,EPS12,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ALPHA,OMHAR
      INTEGER ATOM,IBLOCH,IO,LAY,NRMAX,NTPHOMAX,NVFTPHOMAX
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),CWF(RSTEP,MAXCORE,2),
     &           DOS(MQD,MQD,NATLM,LAYSM),
     &           RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           WFFM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),ZMAT_A(MQD,MQD),
     &           ZMAT_D(MQD,MQD)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        CWFSTATES(NATLM,LAYSM),CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),
     &        CWFY(MAXCORE),WFFXM(MQD,NSTATES,2,NTPHOMAX),
     &        WFFYM(MQD,NSTATES,2,NTPHOMAX)
      REAL*8 CWFG(MAXCSTATES,NATLM,LAYSM)
C
C Local variables
C
      INTEGER CWFS,FINAL,I,INITIAL,J,K1,K2,K3,K4,STATE1,STATE2
      COMPLEX*16 PFI,RMATO(:,:,:,:,:),RMATS(:,:,:,:,:),WFFZ(:,:,:,:,:),
     &           ZMIJ,ZMIJ0
      REAL*8 ZMAX
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE rmato,rmats,wffz
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (RMATO(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (RMATS(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX))
      ALLOCATE (WFFZ(NRMAX,MAXWF,NSTATES,2,NATLM))
      RMATO = CZERO
      RMATS = CZERO
      WFFZ = CZERO
      PFI = CZERO
      DO I = 1,MQD
         DO J = 1,MQD
            ZMAT_A(I,J) = CZERO
            ZMAT_D(I,J) = CZERO
         END DO
      END DO
      STATE1 = 1
      STATE2 = 2
C
C     define initial and final state:
C     ===============================
      IF ( IBLOCH.NE.5 ) THEN
C     xps or xas
         FINAL = STATE1
         INITIAL = STATE2
      ELSE
C     xes
         FINAL = STATE2
         INITIAL = STATE1
      END IF
C
      DO CWFS = 1,CWFSTATES(ATOM,LAY)
C         IF ( CWFG(CWFS,ATOM,LAY).NE.0.D0 ) THEN
         IF ( ABS(CWFG(CWFS,ATOM,LAY)).GT.1.0D-16 ) THEN
C             copy core wavefunction into wff array:
C             ======================================
            CALL CPCOREWFF(WFFM,WFFXM,WFFYM,CWF,CWFX,CWFY,CWFS,STATE2,
     &                     ATOM,LAY,NTPHOMAX,NRMAX)
C
            CALL ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,AMAT2V,
     &                  RDIP1M,WFFM,WFFXM,WFFYM,LAY,ATOM,RMATO,RMATS,
     &                  ALPHA,OMHAR,FINAL,IO,INITIAL,AMAT1T,AMAT2T,WFFZ,
     &                  PFI,0,1,IO)
C
C             calculate z-matrix:
C             ===================
            IF ( IBLOCH.NE.5 ) THEN
C             xps or xas
C             ==========
C                 loop over initial states
               DO K3 = 1,MQD
C                 2 loops over final states
                  DO K1 = 1,MQD
                     IF ( CDABS(RMATO(K3,K1,LAY,ATOM,1)).GE.EPS12 ) THEN
                        DO K2 = 1,MQD
                           IF ( CDABS(RMATO(K3,K2,LAY,ATOM,1))
     &                          .GE.EPS12 ) THEN
C
C                 xps
                              ZMIJ0 = RMATO(K3,K1,LAY,ATOM,1)
     &                                *DCONJG(RMATO(K3,K2,LAY,ATOM,1))
     &                                *CWFG(CWFS,ATOM,LAY)
                              ZMAT_A(K1,K2) = ZMAT_A(K1,K2) + ZMIJ0
                              IF ( IBLOCH.NE.2 ) THEN
C                 xas (dos belongs to final state)
                                 IF ( CDABS(DOS(K1,K2,ATOM,LAY))
     &                                .GE.EPS12 ) THEN
                                    ZMIJ = ZMIJ0*DOS(K3,K3,ATOM,LAY)
                                    ZMAT_D(K1,K2) = ZMAT_D(K1,K2) + ZMIJ
                                 END IF
                              END IF
                           END IF
                        END DO
                     END IF
                  END DO
               END DO
C
            ELSE
C             xes
C             ===
C                 loop over final state
               DO K1 = 1,MQD
C                 2 loops over initial states
                  DO K3 = 1,MQD
                     IF ( CDABS(RMATO(K3,K1,LAY,ATOM,1)).GE.EPS12 ) THEN
                        DO K4 = 1,MQD
                           IF ( CDABS(RMATO(K4,K1,LAY,ATOM,1))
     &                          .GE.EPS12 .AND. 
     &                          CDABS(DOS(K3,K4,ATOM,LAY)).GE.EPS12 )
     &                          THEN
C                     xes (dos belongs to initial state)
                              ZMIJ0 = RMATO(K3,K1,LAY,ATOM,1)
     &                                *DCONJG(RMATO(K4,K1,LAY,ATOM,1))
     &                                *CWFG(CWFS,ATOM,LAY)
                              ZMIJ = ZMIJ0*DOS(K1,K1,ATOM,LAY)
                              ZMAT_A(K3,K4) = ZMAT_A(K3,K4) + ZMIJ0
                              ZMAT_D(K3,K4) = ZMAT_D(K3,K4) + ZMIJ
                           END IF
                        END DO
                     END IF
                  END DO
               END DO
            END IF
         END IF
      END DO
C
C     output:
C     =======
      IF ( IP.GE.4 ) THEN
         WRITE (NOUT1,*) 'zmatrix of atom ',ATOM,' in layer ',LAY
         ZMAX = REAL(ZMAT_A(1,1))
         DO I = 1,MQD
            DO J = 1,MQD
               IF ( CDABS(ZMAT_A(I,J)).GT.ZMAX )
     &              ZMAX = REAL(ZMAT_A(I,J))
            END DO
         END DO
         DO I = 1,MQD
            DO J = 1,MQD
               IF ( CDABS(ZMAT_A(I,J)).GT.ZMAX*0.1D-6 )
     &              WRITE (NOUT1,99001) I,J,ZMAT_A(I,J)
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2I6,2(1x,e15.8))
C
      END
C*==detg.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETG(ENERGY,LAYS,NATL,ASIG,FROMSIGMA,TOSIGMA,EMIN,EMAX,
     &                CWFSTATES,CWFENERGY,CWFKAPMUE,CWFG,ASYM)
C
C     # purpose       : determine the greensfunction energy factor *
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MAXCSTATES,HARTRE,PI,EPS12
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ASIG,LAYS
      REAL*8 ASYM,EMAX,EMIN,FROMSIGMA,TOSIGMA
      COMPLEX*16 ENERGY
      REAL*8 CWFENERGY(MAXCSTATES,NATLM,LAYSM),
     &       CWFG(MAXCSTATES,NATLM,LAYSM),
     &       CWFKAPMUE(MAXCSTATES,NATLM,LAYSM,2)
      INTEGER CWFSTATES(NATLM,LAYSM),NATL(LAYSM)
C
C Local variables
C
      INTEGER ATOM,I,LAY
      COMPLEX*16 CDGAM
      REAL*8 DELTA,DELTA1,EE,G,GMAX,RE
      COMPLEX*16 GAMMA
      EXTERNAL CDGAM
C
C*** End of declarations rewritten by SPAG
C
C     set selfenergy:
C     ===============
      IF ( ASIG.EQ.1 ) THEN
C          write (*,*) fromsigma,tosigma,emax,emin
         DELTA = FROMSIGMA + (TOSIGMA-FROMSIGMA)*(DBLE(ENERGY)-EMIN)
     &           /(EMAX-EMIN)
         DELTA = DELTA/HARTRE
      ELSE
         DELTA = DIMAG(ENERGY)
      END IF
C
      RE = DBLE(ENERGY)
      GMAX = 0.D0
      GAMMA = DCMPLX(1.D0-ASYM,0.D0)
C      IF ( DELTA.EQ.0.D0 ) THEN
      IF ( ABS(DELTA).LE.1.0D-16 ) THEN
         DELTA = EPS12
         WRITE (NOUT1,99001) DELTA
         WRITE (*,99001) DELTA
      END IF
C
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            DO I = 1,CWFSTATES(ATOM,LAY)
               IF ( ASIG.NE.2 ) THEN
                  DELTA1 = DELTA
               ELSE IF ( CWFKAPMUE(I,ATOM,LAY,1).LT.0 ) THEN
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.-2.5D0 )
C     &                 DELTA1 = 3.5D0*DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)+2.5D0).LE.1.0D-16 )
     &                 DELTA1 = 3.5D0*DELTA
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.-1.5D0 )
C     &                 DELTA1 = 3.D0*DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)+1.5D0).LE.1.0D-16 )
     &                 DELTA1 = 3.D0*DELTA
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.-0.5D0 )
C     &                 DELTA1 = 2.5D0*DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)+0.5D0).LE.1.0D-16 )
     &                 DELTA1 = 2.5D0*DELTA
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.0.5D0 )
C     &                 DELTA1 = 2.D0*DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)-0.5D0).LE.1.0D-16 )
     &                 DELTA1 = 2.D0*DELTA
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.1.5D0 )
C     &                 DELTA1 = 1.5D0*DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)-1.5D0).LE.1.0D-16 )
     &                 DELTA1 = 1.5D0*DELTA
C                  IF ( CWFKAPMUE(I,ATOM,LAY,2).EQ.2.5D0 ) DELTA1 = DELTA
                  IF ( ABS(CWFKAPMUE(I,ATOM,LAY,2)-2.5D0).LE.1.0D-16 )
     &                 DELTA1 = DELTA
               ELSE
                  DELTA1 = 5.D0*DELTA
               END IF
C
               EE = RE - CWFENERGY(I,ATOM,LAY)
               G = REAL(CDGAM(GAMMA)
     &             *COS(PI*ASYM/2.D0+(1-ASYM)*ATAN(EE/DELTA1))
     &             /(EE**2+DELTA1**2)**(DBLE(1-ASYM)/2.D0))
C
C               write (*,99002)  cwfkapmue(i,atom,lay,1),
C     1              cwfkapmue(i,atom,lay,2),
C     2              cwfenergy(i,atom,lay),delta1,g,gmax
               IF ( G.GT.GMAX ) GMAX = G
               CWFG(I,ATOM,LAY) = G
            END DO
         END DO
      END DO
C
C      IF ( GMAX.EQ.0.D0 ) THEN
      IF ( ABS(GMAX).LE.1.0D-16 ) THEN
         WRITE (NOUT1,99002)
         WRITE (*,99002)
         STOP
      END IF
C
C     set all cwfg()'s which are 10**(4) times less than the
C     maximum gmax to zero:
C     ======================================================
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            DO I = 1,CWFSTATES(ATOM,LAY)
               IF ( CWFG(I,ATOM,LAY).LT.GMAX*1.D-4 ) CWFG(I,ATOM,LAY)
     &              = 0.D0
            END DO
         END DO
      END DO
C
      IF ( IP.GE.4 ) THEN
         WRITE (NOUT1,99003) RE
         WRITE (NOUT1,99004) DELTA
         WRITE (NOUT1,99005)
         DO LAY = 1,LAYS
            DO ATOM = 1,NATL(LAY)
               DO I = 1,CWFSTATES(ATOM,LAY)
C                  IF ( CWFG(I,ATOM,LAY).NE.0.D0 ) THEN
                  IF ( ABS(CWFG(I,ATOM,LAY)).GT.1.0D-16 ) THEN
                     WRITE (NOUT1,99006) ATOM,LAY,I,CWFG(I,ATOM,LAY),
     &                      CWFENERGY(I,ATOM,LAY)*HARTRE
C                  ELSE IF ( CWFENERGY(I,ATOM,LAY).NE.0.D0 ) THEN
                  ELSE IF ( ABS(CWFENERGY(I,ATOM,LAY)).GT.1.0D-16 ) THEN
                     WRITE (NOUT1,99007) ATOM,LAY,I,CWFG(I,ATOM,LAY),
     &                      CWFENERGY(I,ATOM,LAY)*HARTRE
                  END IF
               END DO
            END DO
         END DO
      END IF
      RETURN
C
99001 FORMAT ('warning in detg: imaginary part of core energy was zero',
     &        ' reset to',e8.1)
99002 FORMAT ('stop in detg :  gmax is zero! -> nothing to do.')
99003 FORMAT ('energy for greensfunction factor:',e13.5)
99004 FORMAT ('                           sigma:',e13.5)
99005 FORMAT ('g(atom,layer,state)=...')
99006 FORMAT ('g(',i3,',',i3,',',i3,') = ',e15.8,15x,'(e=',f12.6,')')
99007 FORMAT ('g(',i3,',',i3,',',i3,') = ',e15.8,' [set to zero]',
     &        ' (e=',f12.6,')')
      END
C*==romatm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ROMATM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,AMAT2V,
     &                  RDIP1M,WFF,WFFX,WFFY,LAY,ATOM,RMATO,RMATS,HM,
     &                  OMHAR,FINAL,IO,INITIAL,AMAT1T,AMAT2T,WFFZ,PFI,
     &                  NCPA,ATA,ICV)
C
C     purpose       : calculates the single radial matrix elements
C                     from munich wave functions
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MAXWF,NSTATES,XMAXE,
     &    NFULLPOT,EPS12,CZERO,CIMAG,NTPHOMAX,NVFTPHOMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATA,ATOM,FINAL,ICV,INITIAL,IO,LAY,NCPA
      REAL*8 HM,OMHAR
      COMPLEX*16 PFI
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),
     &           RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           RMATO(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           RMATS(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           WFF(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),
     &           WFFZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        WFFX(MQD,NSTATES,2,NTPHOMAX),WFFY(MAXWF,NSTATES,NTPHOMAX)
C
C Local variables
C
      INTEGER AA,AMX(:,:,:),AMY(:,:),ANGULAR_TYPE,I1,I2,I3,K1,K2,K3,K4,
     &        LMP,RR
      COMPLEX*16 AMV(:,:),AMVT(:,:),ANG,FF(:),RD1(:),RES
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FF,RD1,AMV,AMX,AMY,AMVT
      ALLOCATE (FF(NRMAX),RD1(NRMAX),AMV(XMAXE,8),AMX(MQD,2,8))
      ALLOCATE (AMY(XMAXE,8),AMVT(XMAXE,8))
C
C*** End of declarations rewritten by SPAG
C
      DO I1 = 1,MQD
         DO I2 = 1,MQD
            RMATO(I1,I2,LAY,ATOM,ICV) = CZERO
            RMATS(I1,I2,LAY,ATOM,ICV) = CZERO
         END DO
      END DO
C
C     next 2 loops are loops over the potential parts:
C
      DO LMP = 1,NFULLPOT
         DO AA = 1,2
            IF ( CDABS(RDIP1M(LMP,AA,RSTEP-30,1,IO)).GE.EPS12 ) THEN
               DO RR = 1,RSTEP
                  RD1(RR) = RDIP1M(LMP,AA,RR,1,IO)
               END DO
C
C         extract angular matrixelements for this l,m component:
C
               CALL CIND(LMP,AA,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,
     &                   AMAT2Y,AMAT2V,AMX,AMY,AMV,AMAT1T,AMAT2T,AMVT)
C
C         calculate radial matrix elements for d type angular:
C         loop over initial state phi(k2,k1):
C
               ANGULAR_TYPE = 1
               DO K2 = 1,MQD
                  DO I1 = WFFX(K2,INITIAL,1,IO),WFFX(K2,INITIAL,2,IO)
                     K1 = WFFY(I1,INITIAL,IO)
                     DO I2 = AMX(K2,1,ANGULAR_TYPE),
     &                  AMX(K2,2,ANGULAR_TYPE)
                        K3 = AMY(I2,ANGULAR_TYPE)
                        ANG = AMV(I2,ANGULAR_TYPE)
                        DO I3 = WFFX(K3,FINAL,1,IO),WFFX(K3,FINAL,2,IO)
                           K4 = WFFY(I3,FINAL,IO)
C
                           IF ( NCPA.GT.0 ) THEN
                              DO RR = 1,RSTEP
                                 FF(RR)
     &                              = (WFFZ(RR,I1,INITIAL,1,IO)*DCONJG
     &                              (WFF(RR,I3,FINAL,1,IO))
     &                              +WFFZ(RR,I1,INITIAL,2,IO)
     &                              *DCONJG(WFF(RR,I3,FINAL,2,IO)))
     &                              *RD1(RR)
                              END DO
                           ELSE
                              DO RR = 1,RSTEP
                                 FF(RR)
     &                              = (WFF(RR,I1,INITIAL,1,IO)*DCONJG
     &                              (WFF(RR,I3,FINAL,1,IO))
     &                              +WFF(RR,I1,INITIAL,2,IO)
     &                              *DCONJG(WFF(RR,I3,FINAL,2,IO)))
     &                              *RD1(RR)
                              END DO
                           END IF
C
                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C
                           IF ( NCPA.EQ.0 ) THEN
                              RMATO(K1,K4,LAY,ATOM,IO)
     &                           = RMATO(K1,K4,LAY,ATOM,IO)
     &                           + RES*DCONJG(ANG)/OMHAR
                           ELSE IF ( ATA.EQ.0 ) THEN
                              RMATO(K1,K4,LAY,ATOM,ICV)
     &                           = RMATO(K1,K4,LAY,ATOM,ICV)
     &                           + RES*DCONJG(ANG)/OMHAR/(-CIMAG*PFI)
                           ELSE IF ( ATA.EQ.1 ) THEN
                              RMATO(K1,K4,LAY,ATOM,ICV)
     &                           = RMATO(K1,K4,LAY,ATOM,ICV)
     &                           + RES*DCONJG(ANG)/OMHAR
                           END IF
C
                           IF ( NCPA.GT.0 ) THEN
                              DO RR = 1,RSTEP
                                 FF(RR)
     &                              = (WFFZ(RR,I3,INITIAL,1,IO)*WFF(RR,
     &                              I1,FINAL,1,IO)
     &                              +WFFZ(RR,I3,INITIAL,2,IO)
     &                              *WFF(RR,I1,FINAL,2,IO))*RD1(RR)
                              END DO
                           ELSE
                              DO RR = 1,RSTEP
                                 FF(RR)
     &                              = (WFF(RR,I3,INITIAL,1,IO)*WFF(RR,
     &                              I1,FINAL,1,IO)
     &                              +WFF(RR,I3,INITIAL,2,IO)
     &                              *WFF(RR,I1,FINAL,2,IO))*RD1(RR)
                              END DO
                           END IF
C
                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C
                           IF ( NCPA.EQ.0 ) THEN
                              RMATS(K1,K4,LAY,ATOM,IO)
     &                           = RMATS(K1,K4,LAY,ATOM,IO)
     &                           - RES*ANG/OMHAR
                           ELSE IF ( ATA.EQ.0 ) THEN
                              RMATS(K1,K4,LAY,ATOM,ICV)
     &                           = RMATS(K1,K4,LAY,ATOM,ICV)
     &                           - RES*ANG/OMHAR/(-CIMAG*PFI)
                           ELSE IF ( ATA.EQ.1 ) THEN
                              RMATS(K1,K4,LAY,ATOM,ICV)
     &                           = RMATS(K1,K4,LAY,ATOM,ICV)
     &                           - RES*ANG/OMHAR
                           END IF
C
                        END DO
                     END DO
                  END DO
C
C         skip other matrixelements (they are very small):
C
               END DO
C
C         calculate radial matrix elements for a+- type angular:
C         loop over final state phi(k2,k1):
C
C               ANGULAR_TYPE = 5
C               DO K2 = 1,MQD
C                  DO I1 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                     K1 = WFFY(I1,FINAL,IO)
C
C               loop over a+-(k2,k3):
C
C                     DO I2 = AMX(K2,1,ANGULAR_TYPE),
C     &                  AMX(K2,2,ANGULAR_TYPE)
C                        K3 = AMY(I2,ANGULAR_TYPE)
C                        ANG = AMV(I2,ANGULAR_TYPE)
C
C                  loop over initial state phi(k3,k4):
C
C                        DO I3 = WFFX(K3,INITIAL,1,IO),
C     &                     WFFX(K3,INITIAL,2,IO)
C                           K4 = WFFY(I3,INITIAL,IO)
C
C                           DO RR = 1,RSTEP
C                              FF(RR) = DCONJG(WFF(RR,I1,FINAL,1,IO))
C     &                                 *WFF(RR,I3,INITIAL,2,IO)
C     &                                 *DCONJG(RDIP2M(LMP,RR,1,IO))
C                           END DO
C
C                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C                           RMATO(K1,K4,LAY,ATOM,IO)
C     &                        = RMATO(K1,K4,LAY,ATOM,IO)
C     &                        - RES*DCONJG(ANG)/OMHAR
C
C                           IF ( FLAG.NE.1 ) THEN
C
C                              DO RR = 1,RSTEP
C                                 FF(RR) = WFF(RR,I1,FINAL,1,IO)
C     &                              *WFF(RR,I3,INITIAL,2,IO)
C     &                              *RDIP2M(LMP,RR,1,IO)
C                              END DO
C                              CALL SIMPSNN(HM,FF,RES,RSTEP)
C                              RMATS(K4,K1,LAY,ATOM,IO)
C     &                           = RMATS(K4,K1,LAY,ATOM,IO)
C     &                           - RES*ANG/OMHAR
C                           END IF
C                        END DO
C                     END DO
C                  END DO
C               END DO
C
C         calculate radial matrix elements for a-+ type angular:
C         loop over final state phi(k2,k1):
C
C               ANGULAR_TYPE = 6
C               DO K2 = 1,MQD
C                  DO I1 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                     K1 = WFFY(I1,FINAL,IO)
C
C               loop over a-+(k2,k3):
C
C                     DO I2 = AMX(K2,1,ANGULAR_TYPE),
C     &                  AMX(K2,2,ANGULAR_TYPE)
C                        K3 = AMY(I2,ANGULAR_TYPE)
C                        ANG = AMV(I2,ANGULAR_TYPE)
C
C                  loop over initial state phi(k3,k4):
C
C                        DO I3 = WFFX(K3,INITIAL,1,IO),
C     &                     WFFX(K3,INITIAL,2,IO)
C                           K4 = WFFY(I3,INITIAL,IO)
C
C                           DO RR = 1,RSTEP
C                              FF(RR) = DCONJG(WFF(RR,I1,FINAL,2,IO))
C     &                                 *WFF(RR,I3,INITIAL,1,IO)
C     &                                 *DCONJG(RDIP2M(LMP,RR,1,IO))
C                           END DO
C
C                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C                           RMATO(K1,K4,LAY,ATOM,IO)
C     &                        = RMATO(K1,K4,LAY,ATOM,IO)
C     &                        - RES*DCONJG(ANG)/OMHAR
C                           IF ( FLAG.NE.1 ) THEN
C
C                              DO RR = 1,RSTEP
C                                 FF(RR) = WFF(RR,I1,FINAL,2,IO)
C     &                              *WFF(RR,I3,INITIAL,1,IO)
C     &                              *RDIP2M(LMP,RR,1,IO)
C                              END DO
C
C                              CALL SIMPSNN(HM,FF,RES,RSTEP)
C                              RMATS(K4,K1,LAY,ATOM,IO)
C     &                           = RMATS(K4,K1,LAY,ATOM,IO)
C     &                           - RES*ANG/OMHAR
C                           END IF
C                        END DO
C                     END DO
C                  END DO
C               END DO
            END IF
         END DO
      END DO
C
      IF ( IP.GE.1 ) THEN
         WRITE (NOUT1,99002) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATO(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATO(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
         WRITE (NOUT1,99003) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATS(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATS(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2I3,2(1x,e15.7))
99002 FORMAT ('radial matrix elements rmato for atom',2x,i3,2x,'layer',
     &        2x,i3,2x,'type',2x,i3)
99003 FORMAT ('radial matrix elements rmats for atom',2x,i3,2x,'layer',
     &        2x,i3,2x,'type',2x,i3)
C
      END
C*==romatd.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ROMATD(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,AMAT2V,
     &                  WFFX,WFFY,LAY,ATOM,RMATO,RMATS,FINAL,IO,INITIAL,
     &                  AMAT1T,AMAT2T,ICV,IC5)
C
C     purpose       : calculates the single radial matrix elements
C                     from munich wave functions
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,MAXWF,NSTATES,XMAXE,NFULLPOT,
     &    EPS12,CZERO,NTPHOMAX,NVFTPHOMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,FINAL,ICV,INITIAL,IO,LAY
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),
     &           RMATO(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           RMATS(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        IC5(MQD,MQD),WFFX(MQD,NSTATES,2,NTPHOMAX),
     &        WFFY(MAXWF,NSTATES,NTPHOMAX)
C
C Local variables
C
      INTEGER AA,AMX(:,:,:),AMY(:,:),ANGULAR_TYPE,I1,I2,I3,K1,K2,K3,K4,
     &        LMP
      COMPLEX*16 AMV(:,:),AMVT(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AMV,AMX,AMY,AMVT
      ALLOCATE (AMV(XMAXE,8),AMX(MQD,2,8))
      ALLOCATE (AMY(XMAXE,8),AMVT(XMAXE,8))
C
C*** End of declarations rewritten by SPAG
C
      DO I1 = 1,MQD
         DO I2 = 1,MQD
            RMATO(I1,I2,LAY,ATOM,ICV) = CZERO
            RMATS(I1,I2,LAY,ATOM,ICV) = CZERO
            IF ( ICV.EQ.1 ) IC5(I1,I2) = 0
         END DO
      END DO
C
C     next 2 loops are loops over the potential parts:
C
      DO LMP = 1,NFULLPOT
         DO AA = 1,2
C
            CALL CIND(LMP,AA,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,
     &                AMAT2V,AMX,AMY,AMV,AMAT1T,AMAT2T,AMVT)
C
C         calculate radial matrix elements for d type angular:
C         loop over initial state phi(k2,k1):
C
            ANGULAR_TYPE = 1
            DO K2 = 1,MQD
               DO I1 = WFFX(K2,INITIAL,1,IO),WFFX(K2,INITIAL,2,IO)
                  K1 = WFFY(I1,INITIAL,IO)
                  DO I2 = AMX(K2,1,ANGULAR_TYPE),AMX(K2,2,ANGULAR_TYPE)
                     K3 = AMY(I2,ANGULAR_TYPE)
                     DO I3 = WFFX(K3,FINAL,1,IO),WFFX(K3,FINAL,2,IO)
                        K4 = WFFY(I3,FINAL,IO)
                        IC5(K1,K4) = 1
                     END DO
                  END DO
               END DO
C
C         skip other matrixelements (they are very small):
C
            END DO
C
C         calculate radial matrix elements for a+- type angular:
C         loop over final state phi(k2,k1):
C
C            ANGULAR_TYPE = 5
C            DO K2 = 1,MQD
C               DO I1 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                  K1 = WFFY(I1,FINAL,IO)
C
C               loop over a+-(k2,k3):
C
C                  DO I2 = AMX(K2,1,ANGULAR_TYPE),AMX(K2,2,ANGULAR_TYPE)
C                     K3 = AMY(I2,ANGULAR_TYPE)
C                     ANG = AMV(I2,ANGULAR_TYPE)
C
C                  loop over initial state phi(k3,k4):
C
C                     DO I3 = WFFX(K3,INITIAL,1,IO),WFFX(K3,INITIAL,2,IO)
C                        K4 = WFFY(I3,INITIAL,IO)
C
C                        DO RR = 1,RSTEP
C                           FF(RR) = DCONJG(WFF(RR,I1,FINAL,1,IO))
C     &                              *WFF(RR,I3,INITIAL,2,IO)
C     &                              *DCONJG(RDIP2M(LMP,RR,1,IO))
C                        END DO
C
C                        CALL SIMPSNN(HM,FF,RES,RSTEP)
C                        RMATO(K1,K4,LAY,ATOM,IO)
C     &                     = RMATO(K1,K4,LAY,ATOM,IO) - RES*DCONJG(ANG)
C     &                     /OMHAR
C
C                        IF ( FLAG.NE.1 ) THEN
C
C                           DO RR = 1,RSTEP
C                              FF(RR) = WFF(RR,I1,FINAL,1,IO)
C     &                                 *WFF(RR,I3,INITIAL,2,IO)
C     &                                 *RDIP2M(LMP,RR,1,IO)
C                           END DO
C                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C                           RMATS(K4,K1,LAY,ATOM,IO)
C     &                        = RMATS(K4,K1,LAY,ATOM,IO) - RES*ANG/OMHAR
C                        END IF
C                     END DO
C                  END DO
C               END DO
C            END DO
C
C         calculate radial matrix elements for a-+ type angular:
C         loop over final state phi(k2,k1):
C
C            ANGULAR_TYPE = 6
C            DO K2 = 1,MQD
C               DO I1 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                  K1 = WFFY(I1,FINAL,IO)
C
C               loop over a-+(k2,k3):
C
C                  DO I2 = AMX(K2,1,ANGULAR_TYPE),AMX(K2,2,ANGULAR_TYPE)
C                     K3 = AMY(I2,ANGULAR_TYPE)
C                     ANG = AMV(I2,ANGULAR_TYPE)
C
C                  loop over initial state phi(k3,k4):
C
C                     DO I3 = WFFX(K3,INITIAL,1,IO),WFFX(K3,INITIAL,2,IO)
C                        K4 = WFFY(I3,INITIAL,IO)
C
C                        DO RR = 1,RSTEP
C                           FF(RR) = DCONJG(WFF(RR,I1,FINAL,2,IO))
C     &                              *WFF(RR,I3,INITIAL,1,IO)
C     &                              *DCONJG(RDIP2M(LMP,RR,1,IO))
C                        END DO
C
C                        CALL SIMPSNN(HM,FF,RES,RSTEP)
C                        RMATO(K1,K4,LAY,ATOM,IO)
C     &                     = RMATO(K1,K4,LAY,ATOM,IO) - RES*DCONJG(ANG)
C     &                     /OMHAR
C                        IF ( FLAG.NE.1 ) THEN
C
C                           DO RR = 1,RSTEP
C                              FF(RR) = WFF(RR,I1,FINAL,2,IO)
C     &                                 *WFF(RR,I3,INITIAL,1,IO)
C     &                                 *RDIP2M(LMP,RR,1,IO)
C                           END DO
C
C                           CALL SIMPSNN(HM,FF,RES,RSTEP)
C                           RMATS(K4,K1,LAY,ATOM,IO)
C     &                        = RMATS(K4,K1,LAY,ATOM,IO) - RES*ANG/OMHAR
C                        END IF
C                     END DO
C                  END DO
C               END DO
C            END DO
         END DO
      END DO
C
      IF ( IP.GE.1 ) THEN
         WRITE (NOUT1,99002) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATO(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATO(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
         WRITE (NOUT1,99003) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATS(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATS(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2I3,2(1x,e15.7))
99002 FORMAT ('radial matrix elements rmato for atom',2x,i3,2x,'layer',
     &        2x,i3,2x,'type',2x,i3)
99003 FORMAT ('radial matrix elements rmats for atom',2x,i3,2x,'layer',
     &        2x,i3,2x,'type',2x,i3)
C
      END
C*==romatmcpa.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ROMATMCPA(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,
     &                     AMAT2V,RDIP1M,WFF,WFFX,WFFY,LAY,ATOM,RMATO,
     &                     RMATS,HM,OMHAR,FINAL,IO,INITIAL,AMAT1T,
     &                     AMAT2T,C1,MAXG,WFFZ,G2,PMS,PMS1,PMS2,ICV)
C
C     purpose       : calculates the single radial matrix elements
C                     from munich wave functions
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MQK,MAXWF,NSTATES,XMAXE,
     &    NFULLPOT,EPS12,CZERO,NTPHOMAX,NVFTPHOMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,FINAL,ICV,INITIAL,IO,LAY,MAXG,PMS
      REAL*8 HM,OMHAR
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),C1(MQD,MQD),
     &           RDIP1M(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           RMATO(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           RMATS(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           WFF(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),
     &           WFFZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),G2(MQK),
     &        PMS1(MQK),PMS2(MQK),WFFX(MQD,NSTATES,2,NTPHOMAX),
     &        WFFY(MAXWF,NSTATES,NTPHOMAX)
C
C Local variables
C
      INTEGER AA,AAMAX,AAMIN,AMX(:,:,:),AMY(:,:),ANGULAR_TYPE,I,I1,I2,
     &        I3,J,K1,K2,K3,K4,LMP,RR
      COMPLEX*16 AMV(:,:),AMVT(:,:),ANG,C2(:,:),FF(:),RES,RESMATO(:,:),
     &           RESMATS(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C2,FF,AMV,AMX,AMY,AMVT,RESMATO,RESMATS
      ALLOCATE (C2(MQD,MQD),FF(NRMAX),AMV(XMAXE,8),AMX(MQD,2,8))
      ALLOCATE (AMY(XMAXE,8),AMVT(XMAXE,8),RESMATO(MAXWF,MAXWF))
      ALLOCATE (RESMATS(MAXWF,MAXWF))
C
      DO I1 = 1,MAXWF
         DO I2 = 1,MAXWF
            RESMATO(I1,I2) = CZERO
            RESMATS(I1,I2) = CZERO
         END DO
      END DO
C
      DO I1 = 1,MQD
         DO I2 = 1,MQD
            C2(I1,I2) = CZERO
            RMATO(I1,I2,LAY,ATOM,ICV) = CZERO
            RMATS(I1,I2,LAY,ATOM,ICV) = CZERO
         END DO
      END DO
C
      DO I = 1,MQD
         DO J = 1,MQD
            C2(I,J) = C1(J,I)
         END DO
      END DO
C
      DO LMP = 1,NFULLPOT
         IF ( NFULLPOT.EQ.1 ) THEN
            AAMIN = 2
            AAMAX = 2
         ELSE
            AAMIN = 1
            AAMAX = 2
         END IF
C
         DO AA = AAMIN,AAMAX
C
            IF ( IP.GE.2 ) WRITE (*,*) 'MAXG,PMS',MAXG,PMS
            DO I1 = 1,MAXG
               DO I2 = 1,MAXG
C
                  DO RR = 1,RSTEP
                     FF(RR) = (WFFZ(RR,I1,INITIAL,1,IO)*DCONJG(WFF(RR,I2
     &                        ,FINAL,1,IO))+WFFZ(RR,I1,INITIAL,2,IO)
     &                        *DCONJG(WFF(RR,I2,FINAL,2,IO)))
     &                        *RDIP1M(LMP,AA,RR,1,IO)
                  END DO
C
                  CALL SIMPSNN(HM,FF,RES,RSTEP)
                  RESMATO(I1,I2) = RESMATO(I1,I2) + RES
C
                  DO RR = 1,RSTEP
                     FF(RR) = (WFFZ(RR,I2,INITIAL,1,IO)*WFF(RR,I1,FINAL,
     &                        1,IO)+WFFZ(RR,I2,INITIAL,2,IO)
     &                        *WFF(RR,I1,FINAL,2,IO))
     &                        *RDIP1M(LMP,AA,RR,1,IO)
                  END DO
C
                  CALL SIMPSNN(HM,FF,RES,RSTEP)
                  RESMATS(I1,I2) = RESMATS(I1,I2) - RES
C
               END DO
            END DO
C
         END DO
      END DO
C
      DO LMP = 1,NFULLPOT
         DO AA = 1,2
            IF ( CDABS(RDIP1M(LMP,AA,RSTEP-30,1,IO)).GE.EPS12 ) THEN
C
               CALL CIND(LMP,AA,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,
     &                   AMAT2Y,AMAT2V,AMX,AMY,AMV,AMAT1T,AMAT2T,AMVT)
C
               ANGULAR_TYPE = 1
               DO K2 = 1,MQD
                  DO I1 = WFFX(K2,INITIAL,1,IO),WFFX(K2,INITIAL,2,IO)
                     K1 = WFFY(I1,INITIAL,IO)
                     DO I2 = AMX(K2,1,ANGULAR_TYPE),
     &                  AMX(K2,2,ANGULAR_TYPE)
                        K3 = AMY(I2,ANGULAR_TYPE)
                        ANG = AMV(I2,ANGULAR_TYPE)
                        DO I3 = WFFX(K3,FINAL,1,IO),WFFX(K3,FINAL,2,IO)
                           K4 = WFFY(I3,FINAL,IO)
C
                           DO I = 1,PMS
                              IF ( K1.EQ.PMS1(I) .AND. K4.EQ.G2(I3) )
     &                             RMATO(K1,K4,LAY,ATOM,ICV)
     &                             = RMATO(K1,K4,LAY,ATOM,ICV)
     &                             + C2(PMS1(I),PMS2(I))*RESMATO(I1,I3)
     &                             *DCONJG(ANG)/OMHAR
                           END DO
C
                           DO I = 1,PMS
                              IF ( K4.EQ.PMS2(I) .AND. K1.EQ.G2(I1) )
     &                             RMATS(K1,K4,LAY,ATOM,ICV)
     &                             = RMATS(K1,K4,LAY,ATOM,ICV)
     &                             + RESMATS(I1,I3)*C1(PMS1(I),PMS2(I))
     &                             *ANG/OMHAR
                           END DO
C
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
C
      IF ( IP.GE.1 ) THEN
         WRITE (NOUT1,99002) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATO(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATO(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
         WRITE (NOUT1,99003) ATOM,LAY,ICV
         DO I1 = 1,MQD
            DO I2 = 1,MQD
               IF ( CDABS(RMATS(I1,I2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99001) I1,I2,RMATS(I1,I2,LAY,ATOM,ICV)
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2I3,2(1x,e15.7))
99002 FORMAT ('CPA matrix elements rmato for atom',2x,i3,2x,'layer',2x,
     &        i3,2x,'type',2x,i3)
99003 FORMAT ('CPA matrix elements rmats for atom',2x,i3,2x,'layer',2x,
     &        i3,2x,'type',2x,i3)
C
      END
C*==cind.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CIND(LMP,A,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,
     &                AMAT2V,AMX,AMY,AMV,AMAT1T,AMAT2T,AMVT)
C
C     purpose:    determine some index fields
C
C     full potential angular matrix elements:
C
      USE MOD_SPEC,ONLY:MQD,XMAXE,NFULLPOT,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER A,LMP
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),AMV(XMAXE,8),AMVT(XMAXE,8)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        AMX(MQD,2,8),AMY(XMAXE,8)
C
C Local variables
C
      INTEGER I,INDEXVAR,MCT
C
C*** End of declarations rewritten by SPAG
C
      DO MCT = 1,8
         DO I = 1,MQD
            AMX(I,1,MCT) = 0
            AMX(I,2,MCT) = -1
         END DO
      END DO
      DO MCT = 1,8
         DO I = 1,XMAXE
            AMY(I,MCT) = 0
            AMV(I,MCT) = CZERO
         END DO
      END DO
C     d++, d--, e++, e--:
C     ===================
      DO MCT = 1,4
         INDEXVAR = 1
         DO I = 1,MQD
            IF ( I.EQ.AMAT1X(INDEXVAR,MCT,LMP,A) ) THEN
               AMX(I,1,MCT) = AMAT1Y(INDEXVAR,MCT,LMP,A)
               AMX(I,2,MCT) = AMAT1Y(INDEXVAR+1,MCT,LMP,A) - 1
               INDEXVAR = INDEXVAR + 1
            ELSE
               AMX(I,1,MCT) = 0
               AMX(I,2,MCT) = -1
            END IF
         END DO
         DO I = 1,XMAXE
            AMY(I,MCT) = AMAT1V(I,MCT,LMP,A)
            AMV(I,MCT) = AMAT1(I,MCT,LMP,A)
            AMVT(I,MCT) = AMAT1T(I,MCT,LMP,A)
         END DO
      END DO
C
C     a+-, a-+, f++, f--:
C     ===================
      DO MCT = 5,8
         INDEXVAR = 1
         DO I = 1,MQD
            IF ( I.EQ.AMAT2X(INDEXVAR,MCT-4,LMP) ) THEN
               AMX(I,1,MCT) = AMAT2Y(INDEXVAR,MCT-4,LMP)
               AMX(I,2,MCT) = AMAT2Y(INDEXVAR+1,MCT-4,LMP) - 1
               INDEXVAR = INDEXVAR + 1
            ELSE
               AMX(I,1,MCT) = 0
               AMX(I,2,MCT) = -1
            END IF
         END DO
         DO I = 1,XMAXE
            AMY(I,MCT) = AMAT2V(I,MCT-4,LMP)
            AMVT(I,MCT) = AMAT2T(I,MCT-4,LMP)
         END DO
      END DO
C
      END
C*==cpcorewff.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CPCOREWFF(WFF,WFFX,WFFY,CWF,CWFX,CWFY,CWFS,STATE,ATOM,
     &                     LAY,NTPHOMAX,NRMAX)
C
C     /****************************************************************/
C     purpose       : copy core wavefunction data from cwf to wff
C     parameters:
C     ===========
C      on input:
C     =========
C      cwf   -- core-wavefunction array
C      cwfx  -- indexfield for core-wavefunction array
C      cwfy  -- indexfield for core-wavefunction array
C      cwfs  -- pointer into cwf array,
C               specifies the core wavefunction to copy
C      atom  -- current atom
C      lay   -- current layer
C      state -- pointer into wff,
C               specify the state where to store the wavefunction
C
C      on return:
C      ==========
C      wff -- array for the wavefunctions
C      wffx,wffy -- index vectors for wavefunction
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MAXWF,NSTATES,MAXCORE,
     &    MAXCSTATES,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,CWFS,LAY,NRMAX,NTPHOMAX,STATE
      COMPLEX*16 CWF(RSTEP,MAXCORE,2),
     &           WFF(NRMAX,MAXWF,NSTATES,2,NTPHOMAX)
      INTEGER CWFX(MQD,MAXCSTATES,NATLM,LAYSM,2),CWFY(MAXCORE),
     &        WFFX(MQD,NSTATES,2,NTPHOMAX),WFFY(MAXWF,NSTATES,NTPHOMAX)
C
C Local variables
C
      INTEGER I,I1,J,K,K1,PT,RR
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,RSTEP
         DO J = 1,MAXWF
            DO K = 1,2
               WFF(I,J,STATE,K,1) = CZERO
            END DO
         END DO
      END DO
      DO K1 = 1,MQD
         WFFX(K1,STATE,1,1) = CWFX(K1,CWFS,ATOM,LAY,1)
     &                        - CWFX(1,CWFS,ATOM,LAY,1) + 1
         WFFX(K1,STATE,2,1) = CWFX(K1,CWFS,ATOM,LAY,2)
     &                        - CWFX(1,CWFS,ATOM,LAY,1) + 1
         DO I1 = CWFX(K1,CWFS,ATOM,LAY,1),CWFX(K1,CWFS,ATOM,LAY,2)
            PT = I1 - CWFX(1,CWFS,ATOM,LAY,1) + 1
            WFFY(PT,STATE,1) = CWFY(I1)
C
            DO RR = 1,RSTEP
               WFF(RR,PT,STATE,1,1) = CWF(RR,I1,1)
               WFF(RR,PT,STATE,2,1) = CWF(RR,I1,2)
            END DO
         END DO
      END DO
C
      END
C*==simpsnn.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SIMPSNN(H,F,QMT,NN)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H
      INTEGER NN
      COMPLEX*16 QMT
      COMPLEX*16 F(NN)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      IF ( MOD(NN,2).EQ.0 ) THEN
         WRITE (*,99001)
         STOP
      END IF
C
      QMT = F(1)
      DO I = 1,NN - 2,2
         QMT = QMT + 4.D0*F(I+1) + 2.D0*F(I+2)
      END DO
      QMT = (H/3.D0)*(QMT-F(NN))
C
      RETURN
C
99001 FORMAT (1x,'** stop in simpsnn: only odd number of points allowed'
     &        )
      END
