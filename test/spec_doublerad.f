C*==doubleradm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C     contains subroutines:
C     doublerad   doubleint  doubleradm   doubleintm
C
      SUBROUTINE DOUBLERADM(AMAT1,AMAT1X,AMAT2X,AMAT1Y,AMAT2Y,AMAT1V,
     &                      AMAT2V,RDIP1,WFF,WFFX,WFFY,LAY,ATOM,ZMAT,HM,
     &                      OMHAR,FINAL,INITIAL,IRSTATE,RADMESH,AMAT1T,
     &                      AMAT2T,IO,ICV)
C
C     # purpose   :calculates double matrix elements for atomic part
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,MAXWF,NSTATES,XMAXE,
     &    NFULLPOT,EPS12,CZERO,NTPHOMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_THERMAL,ONLY:NVFTMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,FINAL,ICV,INITIAL,IO,IRSTATE,LAY
      REAL*8 HM,OMHAR
      COMPLEX*16 AMAT1(XMAXE,4,NFULLPOT,2),AMAT1T(XMAXE,4,NFULLPOT,2),
     &           AMAT2T(XMAXE,4,NFULLPOT),
     &           RDIP1(NFULLPOT,2,NRMAX,2,NTPHOMAX),
     &           WFF(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),
     &           ZMAT(MQD,MQD,LAYSM,NATLM,NVFTMAX)
      INTEGER AMAT1V(XMAXE,4,NFULLPOT,2),AMAT1X(MQD+1,4,NFULLPOT,2),
     &        AMAT1Y(MQD+1,4,NFULLPOT,2),AMAT2V(XMAXE,4,NFULLPOT),
     &        AMAT2X(MQD+1,4,NFULLPOT),AMAT2Y(MQD+1,4,NFULLPOT),
     &        WFFX(MQD,NSTATES,2,NTPHOMAX),WFFY(MAXWF,NSTATES,NTPHOMAX)
      REAL*8 RADMESH(NRMAX,NATLM,LAYSM)
C
C Local variables
C
      INTEGER AA,AMX(:,:,:),AMY(:,:),ANGULAR_TYPE,I,I1,I2,I3,I4,I5,I6,J,
     &        K,K1,K2,K3,K4,K5,K6,K7,LMP,RMAXM,RR,V1(:,:),V2(:,:)
      COMPLEX*16 AMV(:,:),AMVT(:,:),ANG,DUM1,DUM2,FF_G1(:,:,:),
     &           FF_G2(:,:,:),FF_S1(:,:,:),FF_S2(:,:,:),FR1G(:),FR1S(:),
     &           FR2G(:),FR2S(:),RD1(:),RES
      REAL*8 DRE,DRE1,HH3,RM1,RM2,RM3,RM4
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE V1,V2,RD1,AMV,AMX,AMY,FR1G,FR2G,FR1S,FR2S
      ALLOCATABLE AMVT,FF_G1,FF_G2,FF_S1,FF_S2
      ALLOCATE (V1(MQD,MQD),V2(MQD,MQD),RD1(NRMAX))
      ALLOCATE (AMV(XMAXE,8),AMX(MQD,2,8),AMY(XMAXE,8))
      ALLOCATE (FR1G(NRMAX),FR2G(NRMAX),FR1S(NRMAX),FR2S(NRMAX))
      ALLOCATE (AMVT(XMAXE,8),FF_G1(NRMAX,MQD,MQD))
      ALLOCATE (FF_G2(NRMAX,MQD,MQD),FF_S1(NRMAX,MQD,MQD))
      ALLOCATE (FF_S2(NRMAX,MQD,MQD))
C
C*** End of declarations rewritten by SPAG
C
      RMAXM = RSTEP - 1
      HH3 = HM/3.D0
      RM1 = RADMESH(RMAXM,ATOM,LAY)
      RM2 = RADMESH(RSTEP,ATOM,LAY)
      RM3 = RADMESH(1,ATOM,LAY)
      RM4 = RADMESH(2,ATOM,LAY)
C
      DRE = (RADMESH(RSTEP,ATOM,LAY)-RADMESH(RMAXM,ATOM,LAY))/2.D0
      DRE1 = (RADMESH(2,ATOM,LAY)-RADMESH(1,ATOM,LAY))/2.D0
C
      DO I = 1,MQD
         DO J = 1,MQD
            V1(I,J) = 0
            V2(I,J) = 0
         END DO
      END DO
C
      DO I = 1,MQD
         DO J = 1,MQD
            ZMAT(I,J,LAY,ATOM,ICV) = CZERO
         END DO
      END DO
C
      DO I = 1,MQD
         DO J = 1,MQD
            DO K = 1,NRMAX
               FF_S1(K,I,J) = CZERO
               FF_S2(K,I,J) = CZERO
               FF_G1(K,I,J) = CZERO
               FF_G2(K,I,J) = CZERO
            END DO
         END DO
      END DO
C
      DO LMP = 1,NFULLPOT
         DO AA = 1,2
            CALL CIND(LMP,AA,AMAT1X,AMAT1Y,AMAT1V,AMAT1,AMAT2X,AMAT2Y,
     &                AMAT2V,AMX,AMY,AMV,AMAT1T,AMAT2T,AMVT)
C
C         this part is for d-type matrixelements:
C
C            IF ( CDABS(RDIP1(LMP,AA,RSTEP-30,1,IO)).NE.0.0 ) THEN
            IF ( CDABS(RDIP1(LMP,AA,RSTEP-30,1,IO)).GT.1.0D-16 ) THEN
               DO RR = 1,RSTEP
                  RD1(RR) = RDIP1(LMP,AA,RR,1,IO)
               END DO
C
               ANGULAR_TYPE = 1
               DO K2 = 1,MQD
                  DO I1 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
                     K1 = WFFY(I1,FINAL,IO)
                     DO I2 = AMX(K1,1,ANGULAR_TYPE),
     &                  AMX(K1,2,ANGULAR_TYPE)
                        K4 = AMY(I2,ANGULAR_TYPE)
                        ANG = AMV(I2,ANGULAR_TYPE)
                        DO K3 = 1,MQD
                           DO I3 = WFFX(K3,INITIAL,1,IO),
     &                        WFFX(K3,INITIAL,2,IO)
                              IF ( WFFY(I3,INITIAL,IO).EQ.K4 ) THEN
                                 V1(K2,K3) = 1
C
                                 DO RR = 1,RSTEP
                                    DUM1 = WFF(RR,I1,FINAL,1,IO)
     &                                 *ANG*RD1(RR)
                                    DUM2 = WFF(RR,I1,FINAL,2,IO)
     &                                 *ANG*RD1(RR)
C
                                    FF_S1(RR,K2,K3) = FF_S1(RR,K2,K3)
     &                                 + DUM1*WFF(RR,I3,INITIAL,1,IO)
     &                                 + DUM2*WFF(RR,I3,INITIAL,2,IO)
                                    FF_G1(RR,K2,K3) = FF_G1(RR,K2,K3)
     &                                 + DUM1*WFF(RR,I3,IRSTATE,1,IO)
     &                                 + DUM2*WFF(RR,I3,IRSTATE,2,IO)
                                 END DO
C
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
C
               DO K3 = 1,MQD
                  DO I4 = WFFX(K3,INITIAL,1,IO),WFFX(K3,INITIAL,2,IO)
                     K5 = WFFY(I4,INITIAL,IO)
                     DO I5 = AMX(K5,1,ANGULAR_TYPE),
     &                  AMX(K5,2,ANGULAR_TYPE)
                        K7 = AMY(I5,ANGULAR_TYPE)
                        ANG = AMV(I5,ANGULAR_TYPE)
                        DO K6 = 1,MQD
                           DO I6 = WFFX(K6,FINAL,1,IO),
     &                        WFFX(K6,FINAL,2,IO)
                              IF ( WFFY(I6,FINAL,IO).EQ.K7 ) THEN
                                 V2(K3,K6) = 1
C
                                 DO RR = 1,RSTEP
                                    DUM1 = DCONJG(WFF(RR,I6,FINAL,1,IO)
     &                                 *ANG)*RD1(RR)
                                    DUM2 = DCONJG(WFF(RR,I6,FINAL,2,IO)
     &                                 *ANG)*RD1(RR)
C
                                    FF_S2(RR,K3,K6) = FF_S2(RR,K3,K6)
     &                                 + DUM1*WFF(RR,I4,INITIAL,1,IO)
     &                                 + DUM2*WFF(RR,I4,INITIAL,2,IO)
                                    FF_G2(RR,K3,K6) = FF_G2(RR,K3,K6)
     &                                 + DUM1*WFF(RR,I4,IRSTATE,1,IO)
     &                                 + DUM2*WFF(RR,I4,IRSTATE,2,IO)
                                 END DO
C
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
C
               END DO
C
C         this part is for a-type matrixelements:
C
C               IF ( CDABS(RDIP2(LMP,100,1,IO)).NE.0.0 .AND. AA.NE.1 )
C     &              THEN
C                  DO RR = 1,RSTEP
C                     RD2(RR) = RDIP2(LMP,RR,1,IO)
C                  END DO
C
C         calculate radial matrix elements of a type:
C
C                  ANGULAR_TYPE = 5
C                  DO K2 = 1,MQD
C                     DO I1 = AMX(K2,1,ANGULAR_TYPE),
C     &                  AMX(K2,2,ANGULAR_TYPE)
C                        K3 = AMY(I1,ANGULAR_TYPE)
C                        ANG = AMV(I1,ANGULAR_TYPE)
C                        DO I2 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                           K1 = WFFY(I2,FINAL,IO)
C                           DO I3 = WFFX(K3,INITIAL,1,IO),
C     &                        WFFX(K3,INITIAL,2,IO)
C                              K7 = WFFY(I3,INITIAL,IO)
C                              V1(K2,K7) = 1
C
C                              DO RR = 1,RSTEP
C
C                          a+- contribution:
C
C                                 DUM1 = WFF(RR,I1,FINAL,1,IO)*RD2(RR)
C                                 FF_S1(RR,K2,K7) = FF_S1(RR,K2,K7)
C     &                              + DUM1*WFF(RR,I3,INITIAL,2,IO)*ANG
C                                 FF_G1(RR,K2,K7) = FF_G1(RR,K2,K7)
C     &                              + DUM1*WFF(RR,I3,IRSTATE,2,IO)*ANG
C                              END DO
C                           END DO
C                        END DO
C                     END DO
C                  END DO
C
C                  DO K7 = 1,MQD
C                     DO I4 = WFFX(K7,INITIAL,1,IO),WFFX(K7,INITIAL,2,IO)
C                        K4 = WFFY(I4,INITIAL,IO)
C                        DO I5 = AMX(K4,1,ANGULAR_TYPE),
C     &                     AMX(K4,2,ANGULAR_TYPE)
C                           K5 = AMY(I5,ANGULAR_TYPE)
C                           ANG = DCONJG(AMV(I5,ANGULAR_TYPE))
C                           DO I6 = WFFX(K5,FINAL,1,IO),
C     &                        WFFX(K5,FINAL,2,IO)
C                              K6 = WFFY(I6,FINAL,IO)
C                              V2(K7,K6) = 1
C
C                              DO RR = 1,RSTEP
C
C                          a+- contribution:
C
C                                 DUM1 = DCONJG(WFF(RR,I4,FINAL,2,IO)
C     &                                  *RD2(RR))
C                                 FF_S2(RR,K7,K6) = FF_S2(RR,K7,K6)
C     &                              + DUM1*WFF(RR,I6,INITIAL,1,IO)*ANG
C                                 FF_G2(RR,K7,K6) = FF_G2(RR,K7,K6)
C     &                              + DUM1*WFF(RR,I6,IRSTATE,1,IO)*ANG
C                              END DO
C                           END DO
C                        END DO
C                     END DO
C                  END DO
C
C                  ANGULAR_TYPE = 6
C                  DO K2 = 1,MQD
C                     DO I1 = AMX(K2,1,ANGULAR_TYPE),
C     &                  AMX(K2,2,ANGULAR_TYPE)
C                        K3 = AMY(I1,ANGULAR_TYPE)
C                        ANG = AMV(I1,ANGULAR_TYPE)
C                        DO I2 = WFFX(K2,FINAL,1,IO),WFFX(K2,FINAL,2,IO)
C                           K1 = WFFY(I2,FINAL,IO)
C                           DO I3 = WFFX(K3,INITIAL,1,IO),
C     &                        WFFX(K3,INITIAL,2,IO)
C                              K7 = WFFY(I3,INITIAL,IO)
C                              V1(K2,K7) = 1
C
C                              DO RR = 1,RSTEP
C
C                          a+- contribution:
C
C                                 DUM1 = WFF(RR,I1,FINAL,2,IO)*RD2(RR)
C                                 FF_S1(RR,K2,K7) = FF_S1(RR,K2,K7)
C     &                              + DUM1*WFF(RR,I3,INITIAL,1,IO)*ANG
C                                 FF_G1(RR,K2,K7) = FF_G1(RR,K2,K7)
C     &                              + DUM1*WFF(RR,I3,IRSTATE,1,IO)*ANG
C                              END DO
C                           END DO
C                        END DO
C                     END DO
C                  END DO
C
C                  DO K7 = 1,MQD
C                     DO I4 = WFFX(K7,INITIAL,1,IO),WFFX(K7,INITIAL,2,IO)
C                        K4 = WFFY(I4,INITIAL,IO)
C                        DO I5 = AMX(K4,1,ANGULAR_TYPE),
C     &                     AMX(K4,2,ANGULAR_TYPE)
C                           K5 = AMY(I5,ANGULAR_TYPE)
C                           ANG = DCONJG(AMV(I5,ANGULAR_TYPE))
C                           DO I6 = WFFX(K5,FINAL,1,IO),
C     &                        WFFX(K5,FINAL,2,IO)
C                              K6 = WFFY(I6,FINAL,IO)
C                              V2(K7,K6) = 1
C
C                              DO RR = 1,RSTEP
C
C                          a-+ contribution:
C
C                                 DUM1 = DCONJG(WFF(RR,I4,FINAL,1,IO)
C     &                                  *RD2(RR))
C                                 FF_S2(RR,K7,K6) = FF_S2(RR,K7,K6)
C     &                              + DUM1*WFF(RR,I6,INITIAL,2,IO)*ANG
C                                 FF_G2(RR,K7,K6) = FF_G2(RR,K7,K6)
C     &                              + DUM1*WFF(RR,I6,IRSTATE,2,IO)*ANG
C                              END DO
C                           END DO
C                        END DO
C                     END DO
C                  END DO
C               END IF
            END IF
         END DO
      END DO
C
C     integration:
C
      DO K2 = 1,MQD
         DO K4 = 1,MQD
            DO K7 = 1,MQD
               IF ( V1(K2,K7)*V2(K7,K4).EQ.1 ) THEN
C
                  DO RR = 1,RSTEP
                     FR1S(RR) = FF_S1(RR,K2,K7)
                     FR1G(RR) = FF_G1(RR,K2,K7)
                     FR2S(RR) = FF_S2(RR,K7,K4)
                     FR2G(RR) = FF_G2(RR,K7,K4)
                  END DO
C
                  CALL DOUBLEINTM(RES,RSTEP,FR1S,FR1G,FR2S,FR2G,RM1,RM2,
     &                            RM3,RM4,HH3,DRE,DRE1)
C
                  ZMAT(K2,K4,LAY,ATOM,ICV) = ZMAT(K2,K4,LAY,ATOM,ICV)
     &               - RES/(OMHAR*OMHAR)
C
               END IF
            END DO
         END DO
      END DO
C
      IF ( IP.GE.2 ) THEN
         WRITE (NOUT1,99001) ATOM,LAY,ICV
         DO K1 = 1,MQD
            DO K2 = 1,MQD
               IF ( CDABS(ZMAT(K1,K2,LAY,ATOM,ICV)).GT.EPS12 )
     &              WRITE (NOUT1,99002) K1,K2,ZMAT(K1,K2,LAY,ATOM,ICV)
            END DO
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2x,'zmatrix from doubleradm for atom',2x,i3,2x,'layer',2x,
     &        i3,2x,'type',2x,i3)
99002 FORMAT (2x,i3,2x,i3,2x,2E14.7)
C
      END
C*==doubleintm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DOUBLEINTM(RES,RMAX,FR1S,FR1G,FR2S,FR2G,RM1,RM2,RM3,
     &                      RM4,HH3,DRE,DRE1)
C
C     # purpose      : calculate radial integrals in doublerad
C     # calulate     : int(0-r) dr1 {int(0-r1) dr2 {fr1g(r1)*fr2s(r2)}}
C                    + int(0-r) dr1 {int(r1-r) dr2 {fr1s(r1)*fr2g(r2)}}
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DRE,DRE1,HH3,RM1,RM2,RM3,RM4
      COMPLEX*16 RES
      INTEGER RMAX
      COMPLEX*16 FR1G(NRMAX),FR1S(NRMAX),FR2G(NRMAX),FR2S(NRMAX)
C
C Local variables
C
      COMPLEX*16 FF1(:),FF2(:),QMT1,QMT2,TRVAL
      INTEGER I,J,RMAXM,RMAXM2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FF1,FF2
      ALLOCATE (FF1(NRMAX),FF2(NRMAX))
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NRMAX
         FF1(I) = CZERO
         FF2(I) = CZERO
      END DO
C
      RMAXM = RMAX - 1
      RMAXM2 = RMAX - 2
C
C     calculate integral over r2(r1) (odd numbers of r1):
C
      QMT1 = FR2S(1)
      DO J = 1,RMAXM2,2
         QMT1 = QMT1 + 4.D0*FR2S(J+1) + 2.D0*FR2S(J+2)
         FF1(J+2) = (QMT1-FR2S(J+2))
      END DO
C
      QMT2 = FR2G(RMAX)
      DO J = RMAX,3, - 2
         QMT2 = QMT2 + 4.D0*FR2G(J-1) + 2.D0*FR2G(J-2)
         FF2(J-2) = (QMT2-FR2G(J-2))
      END DO
C
C     calculate integral over r2(r1) (even numbers of r1):
C
      TRVAL = (FR2S(2)/RM4+FR2S(1)/RM3)*DRE1/HH3
      FF1(2) = TRVAL
      QMT1 = FR2S(2)
      DO J = 2,RMAXM2,2
         QMT1 = QMT1 + 4.D0*FR2S(J+1) + 2.D0*FR2S(J+2)
         FF1(J+2) = (QMT1-FR2S(J+2)+TRVAL)
      END DO
C
      TRVAL = (FR2G(RMAXM)/RM1+FR2G(RMAX)/RM2)*DRE/HH3
      FF2(RMAX-1) = TRVAL
      QMT2 = FR2G(RMAXM)
      DO J = RMAXM,4, - 2
         QMT2 = QMT2 + 4.D0*FR2G(J-1) + 2.D0*FR2G(J-2)
         FF2(J-2) = (QMT2-FR2G(J-2)+TRVAL)
      END DO
C
      DO I = 1,NRMAX
         FF1(I) = FF1(I)*FR1G(I)
         FF2(I) = FF2(I)*FR1S(I)
      END DO
C
C     determine integrals over r1:
C
      QMT1 = FF1(1)
      QMT2 = FF2(1)
      DO I = 1,RMAXM2,2
         QMT1 = QMT1 + 4.D0*FF1(I+1) + 2.D0*FF1(I+2)
         QMT2 = QMT2 + 4.D0*FF2(I+1) + 2.D0*FF2(I+2)
      END DO
      QMT1 = (QMT1-FF1(RMAX))
      QMT2 = (QMT2-FF2(RMAX))
C
      RES = HH3*HH3*(QMT1+QMT2)
C
      END
C
