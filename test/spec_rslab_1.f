C*==checkpot.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     contains the subroutines:
C         checkpot
C         detradmesh
C         xpsphoto
C         xpsphotoinc
C         cparr
C         cpradmat
C         start0
C         k2kapmue
C         l2kapmue
C         calc_e_vec
C     and functions:
C         k2l2
C         kapmue2l2
C         kapmue2l
C         k2l
C         l2k
C         kapmue2k
C
      SUBROUTINE CHECKPOT(LAYS,NATL,PTYPE)
C     /****************************************************************/
C     # purpose      : check if identical potentials exist;            *
C                      if so determinewhich potentials are equal and   *
C                      store information in ptype(atom,lay,n):         *
C                      n=1 :atomspecifier of equal atom                *
C                      n=2 :layerspecifier of equal atom               *
C                      n=3 :ptype(atom,lay,3)th number                 *
C                           of atom of this type                       *
C                      e.g.:   ptype(2,2,1)=1                          *
C                              ptype(2,2,2)=1                          *
C                              ptype(2,2,3)=4                          *
C                      => means: atom 2 in layer 2 has the same        *
C                         potential as atom 1 in layer 1; it's the     *
C                         4. atom of this potential type               *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SITES,ONLY:NOQ
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYS
      INTEGER NATL(LAYSM),PTYPE(NATLM,LAYSM,3)
C
C Local variables
C
      INTEGER ATOM1,ATOM2,CC1,CC2,I,IQ,J,KC,LAY1,LAY2
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,NATLM
         DO J = 1,LAYSM
            PTYPE(I,J,1) = 0
            PTYPE(I,J,2) = 0
            PTYPE(I,J,3) = 0
         END DO
      END DO
C
      CC1 = 0
      DO LAY1 = 1,LAYS
         DO ATOM1 = 1,NATL(LAY1)
            CC1 = CC1 + 1
            CC2 = 0
            IF ( PTYPE(ATOM1,LAY1,1).EQ.0 ) THEN
               PTYPE(ATOM1,LAY1,1) = ATOM1
               PTYPE(ATOM1,LAY1,2) = LAY1
               PTYPE(ATOM1,LAY1,3) = 1
            END IF
            DO LAY2 = 1,LAYS
               DO ATOM2 = 1,NATL(LAY2)
                  IQ = IQ_SPR_LAYAT(ATOM2,LAY2)
                  DO KC = 1,NOQ(IQ)
                     CC2 = CC2 + 1
C                     IF ( CC2.GT.CC1 .AND. PTYPE(ATOM2,LAY2,1).EQ.0 )
C     &                    THEN
C                        DO RR = 1,RSTEP
C                           DO I = 1,MQD
C                              IF ( DBLE(VLM(LAY1,ATOM1,I,RR,KC))
C     &                             .NE.DBLE(VLM(LAY2,ATOM2,I,RR,KC)) )
C     &                             NEQUAL = 1
C                           END DO
C                        END DO
C                         if (nequal.eq.0) then
C                         potential type atom1,lay1 is equal to
C                         potential type atom2, lay2:
C                         =====================================
C                         write (*,*) 'equal',lay1,atom1,lay2,atom2
C                         ptype(atom2,lay2,1)=atom1
C                         ptype(atom2,lay2,2)=lay1
C                         nptype(atom1,lay1)=nptype(atom1,lay1)+1
C                         ptype(atom2,lay2,3)=nptype(atom1,lay1)
C                         endif
C                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      IF ( IP.GT.0 ) THEN
         WRITE (NOUT1,*) ' '
         WRITE (NOUT1,99001) NATLM,LAYSM
         WRITE (NOUT1,99002)
         DO I = 1,LAYS
            DO J = 1,NATL(I)
               WRITE (NOUT1,99003) I,J,PTYPE(J,I,2),PTYPE(J,I,1),
     &                             PTYPE(J,I,3)
               WRITE (*,99003) I,J,PTYPE(J,I,2),PTYPE(J,I,1),
     &                         PTYPE(J,I,3)
            END DO
         END DO
         WRITE (NOUT1,*) ' '
      END IF
      RETURN
C
99001 FORMAT ('checkpot  natlm,laysm=',2I4)
99002 FORMAT ('the following equal potentials have',' been found:')
99003 FORMAT ('layer',i3,' atom ',i3,' the same as layer',i3,' atom',i3,
     &        ' [',i3,'. atom of this type]')
      END
C*==detradmesh.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETRADMESH(RN,RADMESH,RADMESH3,H)
C     /****************************************************************/
C     # purpose      : determine radial mesh and radial mesh ** 3      *
C                      for all atoms in layers                         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,RSTEP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H(LAYSM,NATLM),RADMESH(RSTEP,NATLM,LAYSM),
     &       RADMESH3(RSTEP,NATLM,LAYSM),RN(NATLM,LAYSM)
C
C Local variables
C
      INTEGER ATOM,I,LAY
C
C*** End of declarations rewritten by SPAG
C
      DO LAY = 1,LAYSM
         DO ATOM = 1,NATLM
            DO I = 1,RSTEP
               RADMESH(I,ATOM,LAY) = RN(ATOM,LAY)
     &                               *EXP(H(LAY,ATOM)*DBLE(I-RSTEP))
               RADMESH3(I,ATOM,LAY) = RADMESH(I,ATOM,LAY)**3
            END DO
         END DO
      END DO
C
      END
C*==detradmeshm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DETRADMESHM(RADMESHM,RADMESH3M,NATL,LAYS,ZM)
C
C     determine muenchner radial mesh and
C     radial mesh**3 for all atoms in layers
C
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R
      USE MOD_TYPES,ONLY:NTMAX,IMT,Z
      USE MOD_SITES,ONLY:ITOQ,NOQ
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC,ONLY:LAYSM,NATLM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LAYS
      INTEGER NATL(LAYSM)
      REAL*8 RADMESH3M(NRMAX,NATLM,LAYSM),RADMESHM(NRMAX,NATLM,LAYSM),
     &       ZM(NATLM,LAYSM,NTMAX)
C
C Local variables
C
      INTEGER ATOM,I,IM,IO,IQ,ITS,JTOP,LAY
C
C*** End of declarations rewritten by SPAG
C
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            IQ = IQ_SPR_LAYAT(ATOM,LAY)
            ITS = ITOQ(1,IQ)
            IM = IMT(ITS)
            JTOP = JRWS(IM)
            DO IO = 1,NOQ(IQ)
               ITS = ITOQ(IO,IQ)
               ZM(ATOM,LAY,IO) = Z(ITS)
            END DO
            DO I = 1,JTOP
               RADMESHM(I,ATOM,LAY) = R(I,IM)
               RADMESH3M(I,ATOM,LAY) = RADMESHM(I,ATOM,LAY)**3
            END DO
         END DO
      END DO
C
      END
C*==xpsphoto.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE XPSPHOTO(CLIGHT,OMHAR,EMV,AWR,ZMAT,NATL,LAYS,LAYER1,
     &                    ROATLA,FLAG,LAYB,IPOL,SPOL,ROATLASD)
C     /****************************************************************/
C     # purpose      :  calculate atomic contribution to the           *
C                       photocurrent                                   *
C                                                                      *
C     # parameter    :                                                 *
C                      flag    =>  type of calculation                 *
C                      =0: xes, xas, or atomic part for xps            *
C                          or 1st ionization in 2ppe                   *
C                          or secondary emission spectrum              *
C                          awr is set to 1                             *
C                      =1: ups                                         *
C                      =2: xps                                         *
C                                                                      *
C     # calls the following function:                                  *
C       k2l2                                                           *
C     /****************************************************************/
C
      USE MOD_SPEC
      USE MOD_SPEC_GEOM,ONLY:TV
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT,OMHAR
      COMPLEX*16 EMV
      INTEGER FLAG,IPOL,LAYB,LAYER1,LAYS,SPOL
      COMPLEX*16 AWR(MQD,LL,NATLM,2),ROATLASD(LL),
     &           ZMAT(MQD,MQD,LAYSM,NATLM)
      INTEGER NATL(LAYSM)
      REAL*8 ROATLA(LL)
C
C Local variables
C
      COMPLEX*16 AC,FAC,FAC1,FAC2,VEC(:),VECC(:),Z(:,:)
      INTEGER ATOM,I,IBULK,IL,IND(:),IT,J,LU
      REAL*8 CLIGHT2
      INTEGER K2L2
      EXTERNAL K2L2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE Z,IND,VEC,VECC
      ALLOCATE (Z(MQD,MQD),IND(MQD),VEC(MQD),VECC(MQD))
C
C*** End of declarations rewritten by SPAG
C
C     determine prefactor of dipole operator
C
      CLIGHT2 = CLIGHT**2
      FAC = 2.D0*CLIGHT2/(2.D0*(EMV+CLIGHT2)+OMHAR)
      FAC2 = FAC*DCONJG(FAC)
      FAC1 = CIMAG*CDSQRT(2.D0*EMV+EMV*EMV/CLIGHT2)/(PI*TV)
C
C     determine field for indexconversion (akm and zmat notation)
C
      IF ( FLAG.NE.0 ) THEN
C         e.g.: ups, xps
         DO I = 1,MQD
            IND(I) = K2L2(I,ML)
         END DO
      END IF
C
      DO I = 1,LL
         ROATLA(I) = 0.0D0
         ROATLASD(I) = CZERO
      END DO
C
      IL = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
            IBULK = LAYB
         END IF
C
         DO IT = IBULK,LAYS
            AC = CZERO
            DO ATOM = 1,NATL(IT)
               DO I = 1,MQD
                  IF ( FLAG.EQ.0 ) THEN
                     VEC(I) = CONE
                     VECC(I) = CONE
                  ELSE IF ( IPOL.EQ.1 ) THEN
                     VEC(I) = AWR(IND(I),IL,ATOM,1)
                     VECC(I) = AWR(IND(I),IL,ATOM,1)
                  ELSE IF ( IPOL.EQ.2 ) THEN
                     VEC(I) = AWR(IND(I),IL,ATOM,2)
                     VECC(I) = AWR(IND(I),IL,ATOM,2)
                  ELSE IF ( IPOL.EQ.3 ) THEN
                     VEC(I) = AWR(IND(I),IL,ATOM,1)
                     VECC(I) = AWR(IND(I),IL,ATOM,2)
                  ELSE IF ( IPOL.EQ.4 ) THEN
                     VEC(I) = AWR(IND(I),IL,ATOM,2)
                     VECC(I) = AWR(IND(I),IL,ATOM,1)
                  END IF
                  DO J = 1,MQD
                     Z(I,J) = ZMAT(I,J,IT,ATOM)
                  END DO
               END DO
C
               DO I = 1,MQD
                  DO J = 1,MQD
                     AC = AC + VEC(I)*Z(I,J)*DCONJG(VECC(J))
C                    ac = ac + vec(i)*dconjg(vec(j))
C                    ac = ac + z(i,j)
                  END DO
               END DO
            END DO
C
            IF ( FLAG.EQ.2 ) THEN
C               xes,xas,xps
               ROATLA(IL) = DBLE(AC*FAC2)
C               ups, xps
            ELSE IF ( SPOL.LE.2 ) THEN
               ROATLA(IL) = DIMAG(AC*FAC1)
            ELSE IF ( SPOL.GT.2 ) THEN
               ROATLASD(IL) = AC*FAC1
            END IF
            IL = IL + 1
         END DO
      END DO
C
      END
C*==xpsphotoinc.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE XPSPHOTOINC(CLIGHT,EMV,AWR,ZMAT,NATL,LAYS,LAYER1,
     &                       ROATLA,FLAG,LAYB,IPOL,SPOL,ROATLASD)
C     /****************************************************************/
C     # purpose      :  calculate incoherent contribution to the       *
C                       photocurrent in CPA case                       *
C                                                                      *
C     # parameter    :                                                 *
C                      flag    =>  type of calculation                 *
C                      =0: xes, xas, or atomic part for xps            *
C                          or 1st ionization in 2ppe                   *
C                          or secondary emission spectrum              *
C                          awr is set to 1                             *
C                      =1: ups                                         *
C                      =2: xps                                         *
C                                                                      *
C     # calls the following function:                                  *
C       k2l2                                                           *
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,LL,PI,CZERO,CONE,CIMAG,
     &    NVFTPHOMAX
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_SPEC_GEOM,ONLY:TV
      USE MOD_THERMAL,ONLY:NVFO_Q
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT
      COMPLEX*16 EMV
      INTEGER FLAG,IPOL,LAYB,LAYER1,LAYS,SPOL
      COMPLEX*16 AWR(MQD,LL,NATLM,2),ROATLASD(LL,NVFTPHOMAX+1),
     &           ZMAT(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX+1)
      INTEGER NATL(LAYSM)
      REAL*8 ROATLA(LL,NVFTPHOMAX+1)
C
C Local variables
C
      INTEGER ATOM,I,IBULK,IL,IND(:),IO,IQ,IT,J,K,LU
      REAL*8 CLIGHT2
      COMPLEX*16 FAC1,VEC(:),VECC(:),Z(:,:)
      INTEGER K2L2
      EXTERNAL K2L2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE Z,IND,VEC,VECC
      ALLOCATE (Z(MQD,MQD),IND(MQD),VEC(MQD),VECC(MQD))
C
C     determine prefactor of dipole operator
C
      CLIGHT2 = CLIGHT**2
      FAC1 = CIMAG*CDSQRT(2.D0*EMV+EMV*EMV/CLIGHT2)/(PI*TV)
C
C     determine field for indexconversion (akm and zmat notation)
C
      IF ( FLAG.NE.0 ) THEN
         DO I = 1,MQD
            IND(I) = K2L2(I,ML)
         END DO
      END IF
C
      DO K = 1,NVFTPHOMAX + 1
         DO I = 1,LL
            ROATLA(I,K) = 0.0D0
            ROATLASD(I,K) = CZERO
         END DO
      END DO
C
      IL = 1
      DO LU = 1,LAYER1
         IF ( LU.EQ.1 ) THEN
            IBULK = 1
         ELSE IF ( LU.GT.1 ) THEN
            IBULK = LAYB
         END IF
C
         DO IT = IBULK,LAYS
            DO ATOM = 1,NATL(IT)
C
               IQ = IQ_SPR_LAYAT(ATOM,IT)
               DO IO = 1,NVFO_Q(IQ) + 1
C
                  DO I = 1,MQD
                     IF ( FLAG.EQ.0 ) THEN
                        VEC(I) = CONE
                        VECC(I) = CONE
                     ELSE IF ( IPOL.EQ.1 ) THEN
                        VEC(I) = AWR(IND(I),IL,ATOM,1)
                        VECC(I) = AWR(IND(I),IL,ATOM,1)
                     ELSE IF ( IPOL.EQ.2 ) THEN
                        VEC(I) = AWR(IND(I),IL,ATOM,2)
                        VECC(I) = AWR(IND(I),IL,ATOM,2)
                     ELSE IF ( IPOL.EQ.3 ) THEN
                        VEC(I) = AWR(IND(I),IL,ATOM,1)
                        VECC(I) = AWR(IND(I),IL,ATOM,2)
                     ELSE IF ( IPOL.EQ.4 ) THEN
                        VEC(I) = AWR(IND(I),IL,ATOM,2)
                        VECC(I) = AWR(IND(I),IL,ATOM,1)
                     END IF
C
                     DO J = 1,MQD
                        Z(I,J) = ZMAT(I,J,IT,ATOM,IO)
                     END DO
                  END DO
C
                  IF ( SPOL.LE.2 ) THEN
                     DO I = 1,MQD
                        DO J = 1,MQD
                           ROATLA(IL,IO) = ROATLA(IL,IO)
     &                        + DIMAG(VEC(I)*Z(I,J)*DCONJG(VECC(J))
     &                        *FAC1)
                        END DO
                     END DO
                  ELSE IF ( SPOL.EQ.4 ) THEN
                     DO I = 1,MQD
                        DO J = 1,MQD
                           ROATLASD(IL,IO) = ROATLASD(IL,IO) - VEC(I)
     &                        *Z(I,J)*DCONJG(VECC(J))*CIMAG*FAC1
                        END DO
                     END DO
                  END IF
C
               END DO
C
            END DO
C
            IL = IL + 1
         END DO
      END DO
C
      END
C*==zmatinco.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ZMATINCO(RMATOTP,RMATSTP,RMATO,RMATS,ZMAT0INC,TAUMT,
     &                    LAY,ATOM,UMAT_VT_PES)
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,CZERO,NVFTPHOMAX
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_THERMAL,ONLY:X_VFT,NVFO_Q,IVFT_VFOQ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,LAY
      COMPLEX*16 RMATO(MQD,MQD,LAYSM,NATLM),
     &           RMATOTP(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           RMATS(MQD,MQD,LAYSM,NATLM),
     &           RMATSTP(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX),
     &           TAUMT(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX+1),
     &           UMAT_VT_PES(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX),
     &           ZMAT0INC(MQD,MQD,LAYSM,NATLM,NVFTPHOMAX+1)
C
C Local variables
C
      COMPLEX*16 C1(:,:),C2(:,:),C3(:,:),C4(:,:)
      REAL*8 CONZ
      INTEGER I,IO,IQ,ITS,J
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE C1,C2,C3,C4
      ALLOCATE (C1(MQD,MQD),C2(MQD,MQD),C3(MQD,MQD),C4(MQD,MQD))
C
C*** End of declarations rewritten by SPAG
C
      IQ = IQ_SPR_LAYAT(ATOM,LAY)
      DO IO = 1,NVFO_Q(IQ)
         ITS = IVFT_VFOQ(IO,IQ)
         CONZ = X_VFT(ITS)
         DO I = 1,MQD
            DO J = 1,MQD
               ZMAT0INC(I,J,LAY,ATOM,IO) = CZERO
            END DO
         END DO
         DO I = 1,MQD
            DO J = 1,MQD
               C1(I,J) = RMATSTP(I,J,LAY,ATOM,IO)
               C2(I,J) = TAUMT(LAY,ATOM,I,J,IO)
               C3(I,J) = RMATOTP(I,J,LAY,ATOM,IO)
            END DO
         END DO
         CALL MULT(C1,C2,C4,MQD)
         CALL MULT(C4,C3,C1,MQD)
C
         DO I = 1,MQD
            DO J = 1,MQD
               C2(I,J) = UMAT_VT_PES(LAY,ATOM,I,J,IO)
            END DO
         END DO
         CALL MULT(C2,C1,C3,MQD)
         DO I = 1,MQD
            DO J = 1,MQD
               C4(I,J) = UMAT_VT_PES(LAY,ATOM,J,I,IO)
            END DO
         END DO
         CALL MULT(C3,C4,C1,MQD)
C
         DO I = 1,MQD
            DO J = 1,MQD
               ZMAT0INC(I,J,LAY,ATOM,IO) = CONZ*C1(I,J)
            END DO
         END DO
      END DO
C
      IF ( NVFO_Q(IQ).GT.1 ) THEN
         DO I = 1,MQD
            DO J = 1,MQD
               C1(I,J) = RMATS(I,J,LAY,ATOM)
               C2(I,J) = TAUMT(LAY,ATOM,I,J,NVFO_Q(IQ)+1)
               C3(I,J) = RMATO(I,J,LAY,ATOM)
            END DO
         END DO
         CALL MULT(C1,C2,C4,MQD)
         CALL MULT(C4,C3,C1,MQD)
C
         DO I = 1,MQD
            DO J = 1,MQD
               ZMAT0INC(I,J,LAY,ATOM,NVFO_Q(IQ)+1) = C1(I,J)
            END DO
         END DO
      ELSE
         DO I = 1,MQD
            DO J = 1,MQD
               ZMAT0INC(I,J,LAY,ATOM,NVFO_Q(IQ)+1) = CZERO
            END DO
         END DO
C
      END IF
C
      END
C*==cparr.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CPARR(ATOM2,LAY2,ATOM1,LAY1,ZMAT0,ZMATD)
C     /****************************************************************/
C     # last purpose :  copies parameters from atom1 to atom2          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM1,ATOM2,LAY1,LAY2
      COMPLEX*16 ZMAT0(MQD,MQD,LAYSM,NATLM),ZMATD(MQD,MQD,LAYSM,NATLM)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MQD
         DO J = 1,MQD
            ZMAT0(I,J,LAY2,ATOM2) = ZMAT0(I,J,LAY1,ATOM1)
            ZMATD(I,J,LAY2,ATOM2) = ZMATD(I,J,LAY1,ATOM1)
         END DO
      END DO
C
      END
C*==cpradmat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CPRADMAT(ATOM2,LAY2,ATOM1,LAY1,RMATO,RMATS,ZMAT)
C     /****************************************************************/
C     # last purpose :  copies parameters from atom1 to atom2          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM1,ATOM2,LAY1,LAY2
      COMPLEX*16 RMATO(MQD,MQD,LAYSM,NATLM),RMATS(MQD,MQD,LAYSM,NATLM),
     &           ZMAT(MQD,MQD,LAYSM,NATLM)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MQD
         DO J = 1,MQD
            RMATO(I,J,LAY2,ATOM2) = RMATO(I,J,LAY1,ATOM1)
            RMATS(I,J,LAY2,ATOM2) = RMATS(I,J,LAY1,ATOM1)
            ZMAT(I,J,LAY2,ATOM2) = ZMAT(I,J,LAY1,ATOM1)
         END DO
      END DO
C
      END
C*==start0.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE START0(Z,NMESH,H,E,VO,CLIGHT,XRL,DKAPPA,RAD,VR,AI1,BI1,
     &                  IREL,ATOM,LAY)
C     /****************************************************************/
C     # purpose      : calculate regular solutions of the              *
C                      radial equation                                 *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,RSTEP,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IREL,LAY,NMESH
      REAL*8 CLIGHT,DKAPPA,H,VO,XRL
      COMPLEX*16 E
      COMPLEX*16 AI1(RSTEP),BI1(RSTEP)
      REAL*8 RAD(RSTEP),VR(RSTEP,NATLM,LAYSM),Z(NATLM,LAYSM)
C
C Local variables
C
      COMPLEX*16 A0X(5),B0X(5),DA(4),DB(4),P12,P21,RP12,RP21,W1,W2
      REAL*8 C,FI,G,H3,Q11,Q12,Q21,Q22,R
      INTEGER I,K
C
C*** End of declarations rewritten by SPAG
C
      C = CLIGHT
      H3 = H/3.D0
      P21 = (((-E))+VO)/C
      P12 = 2.D0*C - P21*XRL
      Q21 = -Z(ATOM,LAY)/C
      Q12 = -Q21*XRL
      IF ( IREL.NE.2 ) THEN
         G = DKAPPA
         Q11 = 0.D0
         Q22 = -2.D0*G
         A0X(1) = CONE
         B0X(1) = DCMPLX(Z(ATOM,LAY)/Q22,0.D0)
      ELSE
         G = SQRT(ABS(DKAPPA**2-Q12**2))
         Q11 = (-DKAPPA) - G
         Q22 = DKAPPA - G
         A0X(1) = DCMPLX(SQRT(ABS(Q22)),0.D0)
         B0X(1) = DCMPLX(SIGN(1.D0,DKAPPA)*SQRT(ABS(Q11)),0.D0)
      END IF
C
      FI = 0.D0
      DO I = 2,5
         FI = FI + 1.D0
         W1 = P12*B0X(I-1)/FI
         W2 = P21*A0X(I-1)/FI
         A0X(I) = W1 + (Q11*W1+Q12*W2)/(FI+2.D0*G)
         B0X(I) = W2 + (Q21*W1+Q22*W2)/(FI+2.D0*G)
      END DO
C
      DO K = 1,4
         R = RAD(K)
         AI1(K) = A0X(1) + R*(A0X(2)+R*(A0X(3)+R*(A0X(4)+R*A0X(5))))
         BI1(K) = B0X(1) + R*(B0X(2)+R*(B0X(3)+R*(B0X(4)+R*B0X(5))))
      END DO
      DO K = 2,4
         DA(K-1) = H3*(Q11*AI1(K)+(RAD(K)*P12+Q12)*BI1(K))
         DB(K-1) = H3*(Q22*BI1(K)+(RAD(K)*P21+Q21)*AI1(K))
      END DO
C
      DO K = 5,NMESH
         R = RAD(K)
         RP21 = -(E*R-VR(K,ATOM,LAY)+Z(ATOM,LAY))/C
         RP12 = 2.D0*C*R - RP21*XRL
         AI1(K) = AI1(K-4) + 8.D0*(DA(3)-0.5D0*DA(2)+DA(1))
         BI1(K) = BI1(K-4) + 8.D0*(DB(3)-0.5D0*DB(2)+DB(1))
         DA(4) = H3*Q11*AI1(K) + H3*RP12*BI1(K)
         DB(4) = H3*Q22*BI1(K) + H3*RP21*AI1(K)
         AI1(K) = AI1(K-1) + 1.125D0*DA(4) + 2.375D0*DA(3)
     &            - 0.625D0*DA(2) + 0.125D0*DA(1)
         BI1(K) = BI1(K-1) + 1.125D0*DB(4) + 2.375D0*DB(3)
     &            - 0.625D0*DB(2) + 0.125D0*DB(1)
         DA(1) = DA(2)
         DB(1) = DB(2)
         DA(2) = DA(3)
         DB(2) = DB(3)
         DA(3) = H3*Q11*AI1(K) + H3*RP12*BI1(K)
         DB(3) = H3*Q22*BI1(K) + H3*RP21*AI1(K)
      END DO
C
      END
C*==k2kapmue.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE K2KAPMUE(K,KAPPA,MUE,ML)
C     /****************************************************************/
C     # purpose      : determine kappa and mue for a given index k     *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K,KAPPA,ML
      REAL*8 MUE
C
C Local variables
C
      INTEGER DUM,KDIV,L,P
C
C*** End of declarations rewritten by SPAG
C
      KDIV = 2*ML + ML*(ML-1)
      IF ( K.LE.KDIV ) THEN
C         kappa is negative:
C         ==================
         P = KDIV - K + 1
C         L = AINT(SQRT(DBLE(P)+0.25D0)-1.5D0+0.999D0)
         L = INT(SQRT(DBLE(P)+0.25D0)-1.5D0+0.999D0)
         KAPPA = -L - 1
         DUM = P - (2*(L+1)+L*(L+1))
         MUE = -ABS(KAPPA) + 0.5D0 - DBLE(DUM)
      ELSE
C         kappa is positive:
C         ==================
         P = K - KDIV
C         L = AINT(SQRT(DBLE(P)+0.25D0)-1.5D0+0.999D0) + 1
         L = INT(SQRT(DBLE(P)+0.25D0)-1.5D0+0.999D0) + 1
         KAPPA = L
         DUM = P - (2*(L)+L*(L-1))
         MUE = ABS(KAPPA) - 0.5D0 + DBLE(DUM)
      END IF
C
      END
C*==l2kapmue.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE L2KAPMUE(P,KAPPA,MUE)
C     /****************************************************************/
C     # purpose      : determine kappa and mue for a given p           *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAPPA,P
      REAL*8 MUE
C
C Local variables
C
      INTEGER DUM,L
C
C*** End of declarations rewritten by SPAG
C
C      L = AINT(SQRT(DBLE(P-1)/2.D0))
      L = INT(SQRT(DBLE(P-1)/2.D0))
      DUM = P - (2*L+2*L*(L-1))
C
      IF ( DUM.LE.2*L ) THEN
         KAPPA = L
      ELSE
         KAPPA = -L - 1
         DUM = DUM - 2*L
      END IF
      MUE = ABS(KAPPA) - 0.5D0 - DBLE(DUM+1)
C
      END
C*==k2l2.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION K2L2(K,ML)
C     /****************************************************************/
C     # purpose      : determine pointer from 'kappa-storage'          *
C                      to 'l-storage'                                  *
C                                                                      *
C     # calls the subroutine and function:                             *
C       k2kapmue      kapmue2l2                                        *
C     /****************************************************************/
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K,ML
C
C Local variables
C
      INTEGER KAPMUE2L2
      INTEGER KAPPA
      REAL*8 MUE
      EXTERNAL KAPMUE2L2
C
C*** End of declarations rewritten by SPAG
C
      CALL K2KAPMUE(K,KAPPA,MUE,ML)
C
      K2L2 = KAPMUE2L2(KAPPA,MUE)
C
      END
C*==kapmue2l2.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION KAPMUE2L2(KAPPA,MUE)
C     /****************************************************************/
C     # purpose       : returns pointer into matrizes                  *
C     /****************************************************************/
C
C     this is for calculation of core wavefunctions;
C     in case of paramagnetic calulation use lowest mue value:
C     ========================================================
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAPPA
      REAL*8 MUE
C
C*** End of declarations rewritten by SPAG
C
C      IF ( MUE.EQ.0.D0 ) MUE = -(ABS(KAPPA)-0.5D0)
      IF ( ABS(MUE).LE.1.0D-16 ) MUE = -(ABS(KAPPA)-0.5D0)
C
      KAPMUE2L2 = 4*(ABS(KAPPA)-1)*(ABS(KAPPA))/2
      IF ( KAPPA.GT.0 ) KAPMUE2L2 = KAPMUE2L2 + 2*ABS(KAPPA)
C      KAPMUE2L2 = KAPMUE2L2 + ABS(KAPPA) - (-MUE-0.5D0)
      KAPMUE2L2 = KAPMUE2L2 + INT(ABS(KAPPA)-(-MUE-0.5D0))
C
      END
C*==kapmue2l.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION KAPMUE2L(KAPPA,MUE)
C     /****************************************************************/
C     # purpose       : returns pointer into matrizes                  *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAPPA
      REAL*8 MUE
C
C*** End of declarations rewritten by SPAG
C
C     this is for calculation of core wavefunctions;
C     in case of paramagnetic calulation use lowest mue value:
C     ========================================================
C      IF ( MUE.EQ.0.D0 ) MUE = -(ABS(KAPPA)-0.5D0)
      IF ( ABS(MUE).LE.1.0D-16 ) MUE = -(ABS(KAPPA)-0.5D0)
C
      KAPMUE2L = 4*(ABS(KAPPA)-1)*(ABS(KAPPA))/2
      IF ( KAPPA.GT.0 ) KAPMUE2L = KAPMUE2L + 2*ABS(KAPPA)
C      KAPMUE2L = KAPMUE2L + ABS(KAPPA) + (-MUE+0.5D0)
      KAPMUE2L = KAPMUE2L + INT(ABS(KAPPA)+(-MUE+0.5D0))
C
      END
C*==k2l.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION K2L(K,ML)
CC     /****************************************************************/
C     # purpose       : determine pointer from 'kappa-storage'         *
C                       to 'l-storage'                                 *
C     # calls the subroutine and function:                             *
C                                                                      *
C       k2kapmue      kapmue2l                                         *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER K,ML
C
C Local variables
C
      INTEGER KAPMUE2L
      INTEGER KAPPA
      REAL*8 MUE
      EXTERNAL KAPMUE2L
C
C*** End of declarations rewritten by SPAG
C
      CALL K2KAPMUE(K,KAPPA,MUE,ML)
C
      K2L = KAPMUE2L(KAPPA,MUE)
C
      END
C*==l2k.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION L2K(L)
C     /****************************************************************/
C     # purpose       : determine pointer from 'l-storage'             *
C                       to 'kappa-storage'                             *
C                                                                      *
C     # calls the subroutine and function:                             *
C       l2kapmue      kapmue2k                                         *
C     /****************************************************************/
      USE MOD_SPEC,ONLY:ML
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L
C
C Local variables
C
      INTEGER KAPMUE2K
      INTEGER KAPPA
      REAL*8 MUE
      EXTERNAL KAPMUE2K
C
C*** End of declarations rewritten by SPAG
C
      CALL L2KAPMUE(L,KAPPA,MUE)
C
      L2K = KAPMUE2K(KAPPA,MUE,ML)
C
      END
C*==kapmue2k.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      INTEGER FUNCTION KAPMUE2K(KAPPA,MUE,MAXL)
C     /****************************************************************/
C     # purpose       : return field index for kappa and mue           *
C                       in case of paramagnetic calculation use mue=0.0*
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KAPPA,MAXL
      REAL*8 MUE
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
C     in case of paramagnetic calulation use lowest mue value:
C     ========================================================
C      IF ( MUE.EQ.0.D0 ) MUE = -(ABS(KAPPA)-0.5D0)
      IF ( ABS(MUE).LE.1.0D-16 ) MUE = -(ABS(KAPPA)-0.5D0)
C
      I = -MAXL
      KAPMUE2K = 1
      DO WHILE ( I.NE.KAPPA )
         KAPMUE2K = KAPMUE2K + ABS(I)*2
         I = I + 1
      END DO
      KAPMUE2K = KAPMUE2K + ABS(KAPPA) + INT(MUE-0.5D0)
C
      END
C*==calc_e_vec.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CALC_E_VEC(PKSCAN,CALC_E,ESTART,EEND,DELTA,NPTS,EVAL)
C     /****************************************************************/
C     # purpose       : calculate energy vector in case of core-level  *
C                       spectra calculation; the mesh of energypoints  *
C                       in the peakregion is calculated finer than in  *
C                       the valley regions.                            *
C                                                                      *
C     # parameter    :                                                 *
C         pkscan : =1: calculate mesh of energypoints in the peakregion*
C                      finer than in the valley region.                *
C                  =2: calculate intensities only for the exact        *
C                             eigenvalues within the energy-interval   *
C                             estart-eend.                             *
C         calc_e :     vector which holds the energy-values where to   *
C                      calculate spectrum.                             *
C         estart,eend: interval in which to calculate the spectrum     *
C         delta:       imaginary part of self energy (parameter)       *
C         npts:        number of points in spectrum (if pkscan=2 then  *
C                      npts is set to the number of core-states within *
C                      he energyinterval estart-eend)                  *
C         eval:        energy-eigenvalues of the core-states           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MAXSPEC,MAXEIG,MAXN,HARTRE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MAXV
      PARAMETER (MAXV=500)
C
C Dummy arguments
C
      REAL*8 DELTA,EEND,ESTART
      INTEGER NPTS,PKSCAN
      REAL*8 CALC_E(MAXSPEC),EVAL(MAXEIG,MAXN,NATLM,LAYSM)
C
C Local variables
C
      REAL*8 ADD,DUM,EMAX,ENERGY,ENERGYVALS(MAXV),EPS,G,INTERVAL,STEP,
     &       SUMVAR,VAL(MAXV)
      INTEGER COUNTVAR,I,INDEXVAR,J,K,L,M
C
C*** End of declarations rewritten by SPAG
C
      INTERVAL = ABS(ESTART-EEND)
      STEP = INTERVAL/(MAXV-1)
C
      IF ( ESTART.GT.EEND ) THEN
         DUM = ESTART
         ESTART = EEND
         EEND = DUM
      END IF
C
      DO I = 1,MAXV
         VAL(I) = 0.D0
      END DO
      DO I = 1,MAXSPEC
         CALC_E(I) = 0.D0
      END DO
C
      IF ( PKSCAN.EQ.2 ) THEN
C         calculate spectrum only for exact core-level evergies:
C         ======================================================
         COUNTVAR = 0
         EMAX = 0
C         find core level energies in energy interval:
C         ============================================
         DO I = 1,MAXEIG
            DO J = 1,MAXN
               DO K = 1,NATLM
                  DO L = 1,LAYSM
                     IF ( EVAL(I,J,K,L)*HARTRE.GT.ESTART .AND. 
     &                    EVAL(I,J,K,L)*HARTRE.LT.EEND ) THEN
                        COUNTVAR = COUNTVAR + 1
                        CALC_E(COUNTVAR) = EVAL(I,J,K,L)*HARTRE
                        EMAX = MAX(EMAX,ABS(CALC_E(COUNTVAR)))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         NPTS = COUNTVAR
         EPS = EMAX*1E-6
C         remove double countings:
C         ========================
         DO I = 1,NPTS - 1
            DO J = I + 1,NPTS
C               IF ( ABS(CALC_E(I)-CALC_E(J)).LT.EPS .AND. CALC_E(J)
C     &              .NE.0.D0 ) THEN
               IF ( ABS(CALC_E(I)-CALC_E(J)).LT.EPS .AND. ABS(CALC_E(J))
     &              .GT.1.0D-16 ) THEN
                  COUNTVAR = COUNTVAR - 1
                  CALC_E(J) = 0.D0
               END IF
            END DO
         END DO
C         sort energy values:
C         ===================
         DO I = 1,NPTS
            DO J = 1,NPTS - 1
               IF ( CALC_E(J).GT.CALC_E(J+1) ) THEN
                  DUM = CALC_E(J)
                  CALC_E(J) = CALC_E(J+1)
                  CALC_E(J+1) = DUM
               END IF
            END DO
         END DO
         NPTS = COUNTVAR
      ELSE
C         put a gaussian to every energy eigenvalue (broaden by delta)
C         this is a test-spectrum which gives an idea about the shape
C         of the real spectrum:
C         =============================================================
         INDEXVAR = 0
         DO I = 1,MAXV
            ENERGY = ESTART + (I-1)*STEP
            INDEXVAR = INDEXVAR + 1
            ENERGYVALS(INDEXVAR) = ENERGY
            DO J = 1,MAXEIG
               DO K = 1,MAXN
                  DO L = 1,NATLM
                     DO M = 1,LAYSM
                        G = DELTA/((ENERGY-EVAL(J,K,L,M)*HARTRE)
     &                      **2+DELTA**2)
                        VAL(INDEXVAR) = VAL(INDEXVAR) + G
                     END DO
                  END DO
               END DO
            END DO
            IF ( INDEXVAR.EQ.MAXV ) EXIT
         END DO
C200      continue
C
C    now spread energypoints where to calculate spectrum over interval;
C    idea for the spreading the points: the integral between two energy-
C    points should be the same, the steeper the test-spectrum is in one
C    region, the more points will be calculated there in comparison to
C    a flat region
C    ===================================================================
C    first integrate the test-spectrum (val):
C    ========================================
         DO I = 1,MAXV - 1
            SUMVAR = SUMVAR + (ENERGYVALS(I+1)-ENERGYVALS(I))
     &               *(VAL(I)+VAL(I+1))/2.D0
         END DO
C         add on:
C         =======
         ADD = 0.5D0*SUMVAR
         SUMVAR = SUMVAR + ADD
         ADD = ADD/MAXV
C
C         calculate integral value which defines step between two points:
C         ===============================================================
         STEP = SUMVAR/NPTS
C         start spreading:
C         ================
         INDEXVAR = 0
         DO I = 1,MAXV - 1
            SUMVAR = SUMVAR + (ENERGYVALS(I+1)-ENERGYVALS(I))
     &               *(VAL(I)+VAL(I+1))/2.D0
            SUMVAR = SUMVAR + ADD
            IF ( SUMVAR.GT.STEP ) THEN
               INDEXVAR = INDEXVAR + 1
               SUMVAR = 0.D0
               IF ( INDEXVAR.GT.MAXSPEC ) THEN
                  WRITE (*,99001)
                  WRITE (*,99002)
                  STOP
               ELSE
                  CALC_E(INDEXVAR) = ENERGYVALS(I)
               END IF
            END IF
         END DO
         NPTS = INDEXVAR
      END IF
C
      RETURN
C
99001 FORMAT ('stop in calc_e_vec')
99002 FORMAT ('array index out of bounds for maxspec!')
      END
