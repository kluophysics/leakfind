C*==spec_sumup_int.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_SUMUP_INT(NEW,JE,JP,JT,EBIND,EHIGH,IBLOCH,LAYP,
     &                          IPOL,IBPOL,IPHPOL,MCD,NT,NP,LL,NOUT1,IP,
     &                          SPOL,ROTRA,ROTER,ROSUR,ROATLA,ROINC,
     &                          ROATSD,ROTRASD,ROTERSD,ROSURSD,ROINCSD,
     &                          ROTRLA,ROTELA,ROSULA,ROTRLASD,ROTELASD,
     &                          ROSULASD,ROLAYRES,ROLAYRESSD,ROATLASD,
     &                          ROATLAINC,ROATLAINCSD,ROALL,ROALLSD,
     &                          NCPA,ATA,PSPIN,THETA,PHI,KXF,KYF,THETAF,
     &                          KPARAF,DETF,ROALLF,VPR,PHOTOC,TESTVAC,
     &                          DETR,DETI,TMINA,DELTAT,ICALL)
C
C   ********************************************************************
C   *                                                                  *
C   *  Summing up diffrent contributions to photocurrent               *
C   *  ICALL = 1: Call from SPOL LOOP                                  *
C   *  ICALL = 2: Call after SPOL LOOP                                 *
C   ********************************************************************
C
      USE MOD_SPEC,ONLY:NVFTPHOMAX,MEW,HARTRE,MPW,PI
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_FILES,ONLY:IFILSPECOU3,IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='ARPES')
C
C Dummy arguments
C
      INTEGER ATA,IBLOCH,IBPOL,ICALL,IP,IPHPOL,IPOL,JE,JP,JT,LAYP,LL,
     &        MCD,NCPA,NEW,NOUT1,NP,NT,SPOL
      REAL*8 DELTAT,DETI,DETR,PHI,ROSUR,ROTER,ROTRA,TESTVAC,THETA,TMINA,
     &       VPR
      COMPLEX*16 EHIGH
      REAL*8 DETF(MPW,2,MEW),EBIND(MEW),KPARAF(MPW,MEW),KXF(MPW,MPW),
     &       KYF(MPW,MPW),PHOTOC(NEW,NT,NP,4,2),PSPIN(3),ROALL(2),
     &       ROALLF(MPW,2,MEW),ROATLA(LL),ROATLAINC(LL,NVFTPHOMAX+1),
     &       ROINC(NVFTPHOMAX+1),ROLAYRES(NT,NP,2,2,LAYP),ROSULA(LL),
     &       ROTELA(LL),ROTRLA(LL),THETAF(MPW)
      COMPLEX*16 ROALLSD(4),ROATLAINCSD(LL,NVFTPHOMAX+1),ROATLASD(LL),
     &           ROATSD(4),ROINCSD(NVFTPHOMAX+1,4),
     &           ROLAYRESSD(1,NT,NP,4,2,LAYP),ROSULASD(LL),ROSURSD(4),
     &           ROTELASD(LL),ROTERSD(4),ROTRASD(4),ROTRLASD(LL)
C
C Local variables
C
      INTEGER I,J
      COMPLEX*16 ISS(2,2),ROATLAINCSSD(4)
      REAL*8 KPARA,KPX,KPY,POL(2),POLV(3),ROAT,ROATLAINCS,SPOLP,TOT
C
C*** End of declarations rewritten by SPAG
C
C
C
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( SPOL.LE.2 ) THEN
            IF ( IP.GT.0 ) WRITE (NOUT1,99008)
         ELSE IF ( SPOL.EQ.4 ) THEN
            IF ( IP.GT.0 ) WRITE (NOUT1,99009)
         END IF
         ROAT = 0.D0
         ROTRA = 0.D0
         ROTER = 0.D0
         ROSUR = 0.D0
         ROINC(:) = 0.0D0
C
         ROATSD(:) = C0
         ROTRASD(:) = C0
         ROTERSD(:) = C0
         ROSURSD(:) = C0
         ROINCSD(:,:) = C0
C
         IF ( IBLOCH.EQ.2 ) THEN
            DO I = 1,LAYP
               ROTRLA(I) = 0.0D0
               ROTELA(I) = 0.0D0
C               ROSULA(I) = 0.0D0
C               ROTRLASD(I) = C0
C               ROTELASD(I) = C0
C               ROSULASD(I) = C0
            END DO
         END IF
C
         IF ( IP.GT.0 ) WRITE (6,*)
         DO I = 1,LAYP
            IF ( SPOL.LE.2 ) THEN
               ROATLA(I) = DABS(ROATLA(I))
               ROAT = ROAT + ROATLA(I)
               ROTRA = ROTRA + ROTRLA(I)
               ROTER = ROTER + ROTELA(I)
               ROSUR = ROSUR + ROSULA(I)
            ELSE IF ( SPOL.EQ.4 ) THEN
               ROATSD(IPOL) = ROATSD(IPOL) + ROATLASD(I)
               ROTRASD(IPOL) = ROTRASD(IPOL) + ROTRLASD(I)
               ROTERSD(IPOL) = ROTERSD(IPOL) + ROTELASD(I)
               ROSURSD(IPOL) = ROSURSD(IPOL) + ROSULASD(I)
            END IF
C
            IF ( SPOL.LE.2 ) THEN
               ROLAYRES(JT,JP,IPOL,IBPOL,I) = ROATLA(I) + ROTRLA(I)
     &            + ROTELA(I) + ROSULA(I)
            ELSE IF ( SPOL.EQ.4 ) THEN
               ROLAYRESSD(1,JT,JP,IPOL,IBPOL,I) = ROATLASD(I)
     &            + ROTRLASD(I) + ROTELASD(I) + ROSULASD(I)
            END IF
C
C
            IF ( IBLOCH.EQ.2 ) THEN
               DO J = 1,NVFTPHOMAX
                  IF ( SPOL.LE.2 ) THEN
                     ROINC(J) = ROINC(J) + ROATLAINC(I,J)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     ROINCSD(J,IPOL) = ROINCSD(J,IPOL)
     &                                 + ROATLAINCSD(I,J)
                  END IF
               END DO
            ELSE IF ( IBLOCH.EQ.4 ) THEN
               IF ( NCPA.GT.0 .AND. ATA.LE.1 ) THEN
                  DO J = 1,NVFTPHOMAX + 1
                     IF ( SPOL.LE.2 ) THEN
                        ROINC(J) = ROINC(J) + ROATLAINC(I,J)
                     ELSE IF ( SPOL.EQ.4 ) THEN
                        ROINCSD(J,IPOL) = ROINCSD(J,IPOL)
     &                     + ROATLAINCSD(I,J)
                     END IF
                  END DO
C
                  DO J = 1,NVFTPHOMAX
                     IF ( SPOL.LE.2 ) THEN
                        ROLAYRES(JT,JP,IPOL,IBPOL,I)
     &                     = ROLAYRES(JT,JP,IPOL,IBPOL,I)
     &                     + ROATLAINC(I,J)
                     ELSE IF ( SPOL.EQ.4 ) THEN
                        ROLAYRESSD(1,JT,JP,IPOL,IBPOL,I)
     &                     = ROLAYRESSD(1,JT,JP,IPOL,IBPOL,I)
     &                     + ROATLAINCSD(I,J)
                     END IF
                  END DO
                  IF ( SPOL.LE.2 ) THEN
                     ROLAYRES(JT,JP,IPOL,IBPOL,I)
     &                  = ROLAYRES(JT,JP,IPOL,IBPOL,I)
     &                  - ROATLAINC(I,NVFTPHOMAX+1)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     ROLAYRESSD(1,JT,JP,IPOL,IBPOL,I)
     &                  = ROLAYRESSD(1,JT,JP,IPOL,IBPOL,I)
     &                  - ROATLAINCSD(I,NVFTPHOMAX+1)
                  END IF
               END IF
            END IF
C
            IF ( NCPA.GT.0 ) THEN
               ROATLAINCS = 0.0D0
               DO J = 1,NVFTPHOMAX
                  IF ( SPOL.LE.2 ) THEN
                     ROATLAINCS = ROATLAINCS + ROATLAINC(I,J)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     ROATLAINCSSD(IPOL) = ROATLAINCSSD(IPOL)
     &                  + ROATLAINCSD(I,J)
                  END IF
               END DO
               IF ( IBLOCH.EQ.2 ) THEN
                  IF ( SPOL.LE.2 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99007) ROATLA(I),
     &                    ROTRLA(I),ROTELA(I),ROSULA(I),
     &                    (ROATLAINC(I,J),J=1,NVFTPHOMAX)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99011) ROATLASD(I),
     &                    ROTRLASD(I),ROTELASD(I),ROSULASD(I),
     &                    (ROATLAINCSD(I,J),J=1,NVFTPHOMAX)
                  END IF
               ELSE IF ( IBLOCH.EQ.4 ) THEN
                  IF ( SPOL.LE.2 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99007) ROATLA(I),
     &                    ROTRLA(I),ROTELA(I),ROSULA(I),
     &                    (ROATLAINC(I,J),J=1,NVFTPHOMAX+1)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99011) ROATLASD(I),
     &                    ROTRLASD(I),ROTELASD(I),ROSULASD(I),
     &                    (ROATLAINCSD(I,J),J=1,NVFTPHOMAX+1)
                  END IF
               END IF
            ELSE IF ( NCPA.EQ.0 ) THEN
               IF ( IBLOCH.EQ.2 ) THEN
                  IF ( IP.GE.0 ) THEN
                  END IF
                  IF ( SPOL.LE.2 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99007) ROATLA(I),
     &                    ROTRLA(I),ROTELA(I),ROSULA(I),
     &                    (ROATLAINC(I,J),J=1,NVFTPHOMAX+1)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99011) ROATLASD(I),
     &                    ROTRLASD(I),ROTELASD(I),ROSULASD(I),
     &                    (ROATLAINCSD(I,J),J=1,NVFTPHOMAX+1)
                  END IF
               ELSE IF ( IBLOCH.EQ.4 ) THEN
                  IF ( SPOL.LE.2 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99007) ROATLA(I),
     &                    ROTRLA(I),ROTELA(I),ROSULA(I)
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     IF ( IP.GT.0 ) WRITE (NOUT1,99011) ROATLASD(I),
     &                    ROTRLASD(I),ROTELASD(I),ROSULASD(I)
                  END IF
               END IF
            END IF
         END DO
C
         IF ( IBLOCH.EQ.4 ) THEN
            IF ( SPOL.LE.2 ) THEN
               ROALL(IPOL) = ROAT + ROTRA + ROTER + ROSUR
            ELSE IF ( SPOL.EQ.4 ) THEN
               ROALLSD(IPOL) = ROATSD(IPOL) + ROTRASD(IPOL)
     &                         + ROTERSD(IPOL) + ROSURSD(IPOL)
            END IF
         ELSE IF ( IBLOCH.EQ.2 ) THEN
            IF ( SPOL.LE.2 ) THEN
               ROALL(IPOL) = ROAT + ROSUR
            ELSE IF ( SPOL.EQ.4 ) THEN
               ROALLSD(IPOL) = ROATSD(IPOL) + ROSURSD(IPOL)
            END IF
         END IF
C
         IF ( IBLOCH.EQ.4 ) THEN
            IF ( NCPA.GT.0 .AND. ATA.LE.1 ) THEN
               IF ( SPOL.LE.2 ) THEN
                  DO I = 1,NVFTPHOMAX
                     ROALL(IPOL) = ROALL(IPOL) + ROINC(I)
                  END DO
                  ROALL(IPOL) = ROALL(IPOL) - ROINC(NVFTPHOMAX+1)
               ELSE IF ( SPOL.EQ.4 ) THEN
                  DO I = 1,NVFTPHOMAX
                     ROALLSD(IPOL) = ROALLSD(IPOL) + ROINCSD(I,IPOL)
                  END DO
                  ROALLSD(IPOL) = ROALLSD(IPOL)
     &                            - ROINCSD(NVFTPHOMAX+1,IPOL)
               END IF
            END IF
         ELSE IF ( IBLOCH.EQ.2 ) THEN
            IF ( SPOL.LE.2 ) THEN
               DO I = 1,NVFTPHOMAX
                  ROALL(IPOL) = ROALL(IPOL) + ROINC(I)
               END DO
            ELSE IF ( SPOL.EQ.4 ) THEN
               DO I = 1,NVFTPHOMAX
                  ROALLSD(IPOL) = ROALLSD(IPOL) + ROINCSD(I,IPOL)
               END DO
            END IF
         END IF
C
         IF ( SPOL.LE.2 ) THEN
            IF ( IP.GT.0 ) THEN
               IF ( NCPA.GT.0 ) THEN
                  WRITE (NOUT1,99006) JE,IPOL,EBIND(JE),ROALL(IPOL)
                  WRITE (NOUT1,99007) ROAT,ROTRA,ROTER,ROSUR,
     &                                (ROINC(I),I=1,NVFTPHOMAX+1)
               ELSE IF ( NCPA.EQ.0 ) THEN
                  WRITE (NOUT1,99006) JE,IPOL,EBIND(JE),ROALL(IPOL)
                  WRITE (NOUT1,99007) ROAT,ROTRA,ROTER,ROSUR
               END IF
            END IF
         ELSE IF ( SPOL.EQ.4 ) THEN
            IF ( IP.GT.0 ) THEN
               IF ( NCPA.GT.0 ) THEN
                  WRITE (NOUT1,99010) JE,IPOL,EBIND(JE),ROALLSD(IPOL)
                  WRITE (NOUT1,99011) ROATSD(IPOL),ROTRASD(IPOL),
     &                                ROTERSD(IPOL),ROSURSD(IPOL),
     &                                (ROINCSD(I,IPOL),I=1,NVFTPHOMAX+1)
               ELSE IF ( NCPA.EQ.0 ) THEN
                  WRITE (NOUT1,99010) JE,IPOL,EBIND(JE),ROALLSD(IPOL)
                  WRITE (NOUT1,99011) ROATSD(IPOL),ROTRASD(IPOL),
     &                                ROTERSD(IPOL),ROSURSD(IPOL)
               END IF
            END IF
         END IF
C
      ELSE IF ( ICALL.EQ.2 ) THEN
         IF ( SPOL.EQ.4 ) THEN
            ISS(1,1) = ROALLSD(1)
            ISS(2,2) = ROALLSD(2)
            ISS(1,2) = ROALLSD(3)
            ISS(2,1) = ROALLSD(4)
C     Rotation in Rashba-Konfiguration with PHI_E
C                             PSPIN(1) = COS(PHID+0.5d0*PI)
C                             PSPIN(2) = SIN(PHID+0.5d0*PI)
C                             PSPIN(3) = 0.0d0
C                             IF ( IPOL.EQ.1 ) THEN
C                                 WRITE (*,*) 'PHIR',PHID*180.0/PI+90.0
C                                 WRITE (*,*) 'PSPINX','PSPINY',
C    &                                  PSPIN(1),PSPIN(2)
C                             END IF
C     in-plane senkrecht zu GM fcc111
C                             PSPIN(1) = 1.0d0
C                             PSPIN(2) = 0.0d0
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GX fcc100
C                             PSPIN(1) = 0.5D0*DSQRT(2.0d0)
C                             PSPIN(2) =-0.5D0*DSQRT(2.0d0)
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GK fcc111
C                             PSPIN(1) = 0.0d0
C                             PSPIN(2) = 1.0d0
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GK hcp111
C            PSPIN(1) = 1.0D0
C            PSPIN(2) = 0.0D0
C            PSPIN(3) = 0.0D0
C     in-plane senkrecht zu GM hcp111
C                             PSPIN(1) = 0.0d0
C                             PSPIN(2) = 1.0d0
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GS bcc110
C                             PSPIN(1) = 0.25D0*DSQRT(2.0d0)
C                             PSPIN(2) =-0.5d0
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GN bcc110
C                             PSPIN(1) = 1.0d0
C                             PSPIN(2) = 0.0d0
C                             PSPIN(3) = 0.0d0
C     in-plane senkrecht zu GH bcc110
C                             PSPIN(1) = 0.0d0
C                             PSPIN(2) = 1.0d0
C                             PSPIN(3) = 0.0d0
C     out of plane
C                             PSPIN(1) = 0.0d0
C                             PSPIN(2) = 0.0d0
C                             PSPIN(3) = 1.0d0
C
            CALL SPEC_VBPESINTTR(ISS,TOT,POL,PSPIN,POLV)
C
            ROALL(1) = POL(1)
            ROALL(2) = POL(2)
         END IF
C
         IF ( IBLOCH.GT.1 ) THEN
            IF ( SPOL.EQ.1 ) ROALL(2) = ROALL(1)
            KPARA = 0.3622605*SIN(THETA)
     &              *SQRT(2.0*HARTRE*(CDABS(EHIGH)-VPR))
            KPX = 0.3622605*SIN(THETA)
     &            *SQRT(2.0*HARTRE*(CDABS(EHIGH)-VPR))*COS(PHI)
            KPY = 0.3622605*SIN(THETA)
     &            *SQRT(2.0*HARTRE*(CDABS(EHIGH)-VPR))*SIN(PHI)
            KXF(JT,JP) = KPX
            KYF(JT,JP) = KPY
C
            IF ( IP.GT.0 .AND. SPOL.EQ.4 ) WRITE (6,*) KPX,KPY,POLV(1),
     &           POLV(2),POLV(3)
            IF ( IP.GT.0 ) WRITE (6,*) THETA,EHIGH,KPARA
C
            IF ( MCD.EQ.0 ) THEN
               IF ( SPOL.LE.2 ) THEN
                  DO IPOL = 1,SPOL
                     PHOTOC(JE,JT,JP,IPOL,IPHPOL)
     &                  = PHOTOC(JE,JT,JP,IPOL,IPHPOL) + ROALL(IPOL)
                  END DO
               ELSE IF ( SPOL.EQ.4 ) THEN
                  DO IPOL = 1,SPOL - 2
                     PHOTOC(JE,JT,JP,IPOL,IPHPOL)
     &                  = PHOTOC(JE,JT,JP,IPOL,IPHPOL) + ROALL(IPOL)
                  END DO
               END IF
            ELSE IF ( SPOL.LE.2 ) THEN
               DO IPOL = 1,SPOL
                  PHOTOC(JE,JT,JP,IPOL,IBPOL)
     &               = PHOTOC(JE,JT,JP,IPOL,IBPOL) + ROALL(IPOL)
               END DO
            ELSE IF ( SPOL.EQ.4 ) THEN
               DO IPOL = 1,SPOL - 2
                  PHOTOC(JE,JT,JP,IPOL,IBPOL)
     &               = PHOTOC(JE,JT,JP,IPOL,IBPOL) + ROALL(IPOL)
               END DO
            END IF
C
            IF ( TESTVAC.LT.0.0D0 ) THEN
               SPOLP = 0.0D0
               DETR = 0.0D0
               DETI = 0.0D0
            ELSE
               SPOLP = (ROALL(1)-ROALL(2))/(ROALL(1)+ROALL(2))
            END IF
C
            IF ( NT.EQ.1 .AND. NP.EQ.1 ) THEN
               WRITE (6,99002) EBIND(JE),(ROALL(1)+ROALL(2))/2.D0,SPOLP,
     &                         0.5*ROALL(1),0.5*ROALL(2),
     &                         DETR**2 + DETI**2,KPARA
               WRITE (IFILSPECOU3,99001) EBIND(JE),(ROALL(1)+ROALL(2))
     &                /2.D0,SPOLP,0.5*ROALL(1),0.5*ROALL(2),
     &                DETR**2 + DETI**2,KPARA
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT .AND. JE.LE.4 )
     &              WRITE (IFILBUILDBOT,99012)
     &              ROUTINE(1:LEN_TRIM(ROUTINE)),JT,JP,JE,EBIND(JE),
     &              ROALL(1),ROALL(2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            ELSE IF ( NT.GT.1 .AND. NP.EQ.1 ) THEN
               WRITE (6,99003) KPARA,EBIND(JE),(ROALL(1)+ROALL(2))/2.D0,
     &                         SPOLP,0.5*ROALL(1),0.5*ROALL(2)
               WRITE (IFILSPECOU3,99005) KPARA,EBIND(JE),
     &                (ROALL(1)+ROALL(2))/2.D0,SPOLP,0.5*ROALL(1),
     &                0.5*ROALL(2)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT .AND. JE.LE.4 .AND. JT.LE.2 )
     &              WRITE (IFILBUILDBOT,99012)
     &              ROUTINE(1:LEN_TRIM(ROUTINE)),JT,JP,JE,EBIND(JE),
     &              ROALL(1),ROALL(2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            ELSE
               WRITE (6,99004) KPX,KPY,(ROALL(1)+ROALL(2))/2.D0,SPOLP,
     &                         0.5*ROALL(1),0.5*ROALL(2)
               WRITE (IFILSPECOU3,99005) KPX,KPY,(ROALL(1)+ROALL(2))
     &                /2.D0,SPOLP,0.5*ROALL(1),0.5*ROALL(2)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
               IF ( WRBUILDBOT .AND. JE.LE.4 .AND. JT.LE.2 .AND. 
     &              JP.LE.2 ) WRITE (IFILBUILDBOT,99012)
     &                               ROUTINE(1:LEN_TRIM(ROUTINE)),JT,JP,
     &                               JE,EBIND(JE),ROALL(1),ROALL(2)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
            END IF
C
            THETAF(JT) = (TMINA+DELTAT*DBLE(JT-1))*180.0/PI
            KPARAF(JT,JE) = KPARA
C             thetaf(jt)=kpara
            ROALLF(JT,1,JE) = ROALL(1)
            ROALLF(JT,2,JE) = ROALL(2)
            DETF(JT,1,JE) = DETR
            DETF(JT,2,JE) = DETI
C
         END IF
      END IF
C
99001 FORMAT (2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,
     &        e14.5,2x,e14.5)
99002 FORMAT (2x,'e=',f10.3,2x,'i=',e14.5,2x,'p=',e14.5,2x,'iup=',e14.5,
     &        2x,'idn=',e14.5,2x,'det=',e14.5,2x,'kp=',e14.5)
99003 FORMAT (2x,'kp=',f10.3,2x,'e=',e14.5,2x,'i=',e14.5,2x,'p=',e14.5,
     &        2x,'iun=',e14.5,2x,'idn=',e14.5)
99004 FORMAT (2x,'kx=',e14.5,2x,'ky=',e14.5,2x,'i=',e14.5,2x,'p=',e14.5,
     &        2x,'iup=',e14.5,2x,'idn=',e14.5)
99005 FORMAT (2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,
     &        e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5,2x,e14.5)
99006 FORMAT (1x,i4,1x,'s=',i1,2(2x,e14.7))
99007 FORMAT (10(2x,e12.4))
99008 FORMAT (3x,'roatom',10x,'rointra',9x,'rointer',9x,'rosurf')
99009 FORMAT (15x,'roatom',20x,'rointra',20x,'rointer',20x,'rosurf')
99010 FORMAT (1x,i4,1x,'s=',i1,4(2x,e14.7))
99011 FORMAT (20(2x,e12.4))
99012 FORMAT ('# BUILDBOT: ',A,': (ROALL(ISPIN))',
     &        ',ISPIN=1,2)  for JT=',I5,4x,'JP=',I5,/,'#',10X,
     &        'energy  JE =',I5,' EBIND = ',F10.6,/,(1PE22.14))
C
      END
C*==spec_vbpesinttr.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_VBPESINTTR(ISS,TOT,POL,NSPIN,P)
C   ********************************************************************
C   *                                                                  *
C   *             get polarisation and total intensity                 *
C   *             (copy of VBPESINTTR)                                 *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0,CI,C1
      PARAMETER (C0=(0.0D0,0.0D0),CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      REAL*8 TOT
      COMPLEX*16 ISS(2,2)
      REAL*8 NSPIN(3),P(3),POL(2)
C
C Local variables
C
      REAL*8 DDOT
      INTEGER I,IMS,IMSP
      COMPLEX*16 PAULI(2,2,1:3),POLMAT(2,2),RHO(2,2)
      REAL*8 PN
C
C*** End of declarations rewritten by SPAG
C
C
      PAULI(1,1,1) = C0
      PAULI(2,1,1) = C1
      PAULI(1,2,1) = C1
      PAULI(2,2,1) = C0
C
      PAULI(1,1,2) = C0
      PAULI(2,1,2) = CI
      PAULI(1,2,2) = -CI
      PAULI(2,2,2) = C0
C
      PAULI(1,1,3) = C1
      PAULI(2,1,3) = C0
      PAULI(1,2,3) = C0
      PAULI(2,2,3) = -C1
C
C--------------------------------------------------- spin density matrix
C
      DO IMS = 1,2
         DO IMSP = 1,2
            RHO(IMS,IMSP) = (ISS(IMS,IMSP)-DCONJG(ISS(IMSP,IMS)))
     &                      /(2.0D0*CI)
         END DO
      END DO
C
C------------------------------------------------------- total intensity
C
      TOT = DREAL(RHO(1,1)+RHO(2,2))
C
C--------------------------------------------------- polarisation vector
C
      DO I = 1,3
         CALL CMATMUL(2,2,PAULI(1,1,I),RHO,POLMAT(1,1))
         IF ( ABS(TOT).GT.1D-12 ) THEN
            P(I) = DREAL(POLMAT(1,1)+POLMAT(2,2))/TOT
         ELSE
            P(I) = 0.0D0
         END IF
      END DO
C
C
C------------------------------------------ spin-projected photo-current
C
      PN = DDOT(3,NSPIN,1,P(1),1)
C
      POL(1) = (1D0-PN)*TOT/2D0
      POL(2) = (1D0+PN)*TOT/2D0
C
C
      END
