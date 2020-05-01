C*==potbar.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE POTBAR(IBAR,ZPARU,ZPARD,BPAR,EPSX,BARABU,BARABD,INPVER)
C     /****************************************************************/
C     # purpose      : read parameter for surface potential            *
C     # note         : file with lu=10 is opened in openfiles          *
C     /****************************************************************/
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      USE MOD_FILES,ONLY:IFILSPECBAR,FOUND_SECTION,FOUND_REAL,
     &    FOUND_REAL_ARRAY,N_FOUND
      USE MOD_SPEC,ONLY:TRANSP_BAR
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 BARABD,BARABU,EPSX
      INTEGER IBAR,INPVER
      REAL*8 BPAR(3),ZPARD(3),ZPARU(3)
C
C Local variables
C
      LOGICAL FOUND
      INTEGER I,LU_POTBAR
      REAL*8 RINP(2),RINP3(3)
C
C*** End of declarations rewritten by SPAG
C
C
C     ibar = controls type of surface barrier
C     ibar : = 0 step, = 1 rm-barrier, = 2 rm + gauss
C     epsx  = tol for barrier integration
C     barabu = spin up:
C              distance of surface layer from the first bulk-layer
C     barabd = spin down:
C              distance of surface layer from the first bulk-layer
C     zpar(1:3) = rm potential parameter
C     bpar(1:3) = gaussian baqrrier parameter (distance, width, height)
C
C     TRANSP_BAR = transparent barier
      TRANSP_BAR = .FALSE.
C
      IF ( INPVER.EQ.1 ) THEN
         IBAR = 1
         EPSX = 0.5D-02
         BARABU = 0.25D0
         BARABD = 0.25D0
         ZPARU(1) = 0.0D0
         ZPARU(2) = -1.0D0
         ZPARU(3) = 1.0D0
         ZPARD(1) = 0.0D0
         ZPARD(2) = -1.0D0
         ZPARD(3) = 1.0D0
C
         BPAR(1) = 1.34D0
         BPAR(2) = -2.0D0
         BPAR(3) = 1.018D0
C
         CALL INPUT_FIND_SECTION('SPEC_STR',0)
C
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_REAL_ARRAY('SURF_BAR',RINP,N_FOUND,2,0,
     &                                  9999D0,0)
C
            IF ( FOUND_REAL_ARRAY ) THEN
               BARABU = RINP(1)
               BARABD = RINP(2)
            ELSE
               CALL SECTION_SET_REAL('BARABU',BARABU,9999D0,0)
               CALL SECTION_SET_REAL('BARABD',BARABD,9999D0,0)
            END IF
            CALL SECTION_SET_INTEGER('IBAR',IBAR,9999,0)
            CALL SECTION_SET_REAL('EPSX',EPSX,9999D0,0)
            IF ( FOUND_REAL ) EPSX = MAX(1.D-6,EPSX)
            CALL SECTION_SET_REAL_ARRAY('ZPARU',RINP3,N_FOUND,3,0,
     &                                  9999D0,0)
            IF ( FOUND_REAL_ARRAY ) ZPARU(1:3) = RINP3(1:3)
            CALL SECTION_SET_REAL_ARRAY('ZPARD',RINP3,N_FOUND,3,0,
     &                                  9999D0,0)
            IF ( FOUND_REAL_ARRAY ) ZPARD(1:3) = RINP3(1:3)
            CALL SECTION_SET_REAL_ARRAY('BPAR',RINP3,N_FOUND,3,0,9999D0,
     &                                  0)
            IF ( FOUND_REAL_ARRAY ) BPAR(1:3) = RINP3(1:3)
C
            CALL SECTION_FIND_KEYWORD('TRANSP_BAR',FOUND)
            IF ( FOUND ) TRANSP_BAR = .TRUE.
         END IF
         WRITE (6,99011) 'Surface barier:'
         IF ( TRANSP_BAR ) THEN
            WRITE (6,99011) '      Transparent barier will be used'
         ELSE
            WRITE (6,99012) 'Distance from the first layer Spin UP:',
     &                      BARABU
            WRITE (6,99012) 'Distance from the first layer Spin DN:',
     &                      BARABD
         END IF
         WRITE (6,99013)
      ELSE
         LU_POTBAR = IFILSPECBAR
         REWIND LU_POTBAR
C
         READ (LU_POTBAR,99001) IBAR
         READ (LU_POTBAR,99002) EPSX,BARABU
         READ (LU_POTBAR,99003) (ZPARU(I),I=1,3)
C
         READ (LU_POTBAR,99002) EPSX,BARABD
         EPSX = MAX(1.D-6,EPSX)
C
C     zpar(1:3) = rm potential parameter
         READ (LU_POTBAR,99003) (ZPARD(I),I=1,3)
C
         READ (LU_POTBAR,99003) (BPAR(I),I=1,3)
C
      END IF
      WRITE (NOUT1,99004)
      WRITE (NOUT1,99005) IBAR
      WRITE (NOUT1,99006) EPSX
      WRITE (NOUT1,99007) (ZPARU(I),I=1,3)
      WRITE (NOUT1,99009) (ZPARD(I),I=1,3)
      WRITE (NOUT1,99008) (BPAR(I),I=1,3)
      WRITE (NOUT1,99010)
C
      RETURN
C
99001 FORMAT (1I3)
99002 FORMAT (2E14.7)
99003 FORMAT (3E14.7)
C
99004 FORMAT (1x,'parameter for potential barrier:')
99005 FORMAT (1x,'# ibar:',i3)
99006 FORMAT (1x,'# epsx:',e14.7)
99007 FORMAT (1x,'# zparup:',3E14.7)
99008 FORMAT (1x,'# bparp:',3E14.7)
99009 FORMAT (1x,'# zpardn:',3E14.7)
99010 FORMAT (1x,'=======================')
99011 FORMAT (//,10X,A)
99012 FORMAT (10X,A,5F10.6)
99013 FORMAT (/,1X,79('-'),//)
      END
C*==esenk.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ESENK(GANZ,QD,VPR,EMV,EMV1,ILH,IREL,CLIGHT,KGZ,EZRY,
     &                 KGT,BK)
C
C     /****************************************************************/
C     # purpose      : determine perpendicular part of e-momentum      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LG,LH,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT,VPR
      COMPLEX*16 EMV,EMV1
      INTEGER GANZ,ILH,IREL,QD
      COMPLEX*16 BK(LH,2),KGZ(LH)
      REAL*8 EZRY(LH),KGT(2,LG)
C
C Local variables
C
      REAL*8 AE,AEH
      COMPLEX*16 BX,E
      INTEGER I,IS
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,LH
         EZRY(I) = 0.0
         KGZ(I) = CZERO
         BK(I,1) = CZERO
         BK(I,2) = CZERO
      END DO
C
      E = EMV1
      IF ( ILH.EQ.2 ) E = EMV
C
C     hr = (30.d0 - zz(3,1))/mhmax
      DO IS = 1,GANZ
         BK(IS,1) = DCMPLX(KGT(1,IS),0.D0)
         BK(IS,2) = DCMPLX(KGT(2,IS),0.D0)
      END DO
C
C     hr = (30.d0 - zz(3,2))/mhmax
      DO IS = GANZ + 1,QD
         BK(IS,1) = DCMPLX(KGT(1,IS-GANZ),0.D0)
         BK(IS,2) = DCMPLX(KGT(2,IS-GANZ),0.D0)
      END DO
C
      IF ( IP.GT.1 ) WRITE (NOUT1,99001)
C
      DO IS = 1,QD
         BX = BK(IS,1)**2 + BK(IS,2)**2
         IF ( IREL.EQ.2 ) THEN
            AE = DBLE(2.D0*E-BX) + DBLE(E/CLIGHT)**2
         ELSE
            AE = DBLE(2.D0*E-BX)
         END IF
         AEH = AE/2.D0
         EZRY(IS) = AEH - VPR
         KGZ(IS) = CDSQRT(DCMPLX(2.D0*EZRY(IS),0.D0))
         IF ( IP.GT.1 ) WRITE (NOUT1,99002) ILH,IS,EZRY(IS),KGZ(IS)
      END DO
      RETURN
C
99001 FORMAT (1x,'energy ezrg calculated for rm-potential')
99002 FORMAT (1x,'h=',i3,2x,'rv=',i3,2x,'ez=',e12.5,2x,'kz=',e12.5,
     &        e12.5)
      END
C*==potin.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE POTIN(VPR,VI,ILH,U,D,ZAU,ZSU,ZBU,ZAD,ZSD,ZBD)
C
C     /****************************************************************/
C     # purpose      : potin initializes values in common block compot *
C                      and common block poly that are later used by    *
C                      subroutine cbpot.                               *
C                      zz(i), i=1,6    are the position parameters     *
C                      vv(i) are the potential parameters              *
C                                                                      *
C     # note         : the potential barrier is given in ry            *
C                      1 hartree = 2 ry                                *
C                                                                      *
C     # subroutines called from this routine:                          *
C       herpol                                                         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:HARTRE
      USE MOD_SPEC_COMPOT,ONLY:ZZ,A,B,VV,KTRB
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ILH
      REAL*8 VI,VPR,ZAD,ZAU,ZBD,ZBU,ZSD,ZSU
      REAL*8 D(3),U(3)
C
C Local variables
C
      REAL*8 BOHRI,H1,H2,H3,RYDBI,V(5,2,2),Z(7,2),Z113,Z12,Z123,Z13
      INTEGER I,IMM
C
C*** End of declarations rewritten by SPAG
C
C     bohr = 0.52917706d0
C     rydb = 13.605804d0  = 1/2 hartree
C
      BOHRI = 1.D0/0.52917706D0
      RYDBI = 1.D0/13.605804D0
C
      KTRB = 0
C
      DO I = 1,5
         VV(I,ILH,1) = 0.D0
         VV(I,ILH,2) = 0.D0
      END DO
C
      DO I = 1,4
         A(I,1) = 0.D0
         A(I,2) = 0.D0
         B(I,1) = 0.D0
         B(I,2) = 0.D0
      END DO
C
      Z(1,1) = U(1)
      Z(2,1) = U(2)
      Z(3,1) = U(3)
C
      Z(1,2) = D(1)
      Z(2,2) = D(2)
      Z(3,2) = D(3)
C
      DO IMM = 1,2
         Z(4,IMM) = Z(2,IMM)
         Z(5,IMM) = Z(3,IMM)
         Z(6,IMM) = 0.D0
         Z(7,IMM) = 0.D0
C
         V(1,ILH,IMM) = 0.D0
         V(2,ILH,IMM) = 0.D0
         V(3,ILH,IMM) = -VPR*HARTRE
         V(4,ILH,IMM) = -VI
         V(5,ILH,IMM) = 0.D0
C
         IF ( KTRB.NE.1 ) THEN
C
C           rm potential
C
            DO I = 1,7
               ZZ(I,IMM) = Z(I,IMM)*BOHRI
            END DO
C
            IF ( ZZ(2,IMM).GT.ZZ(3,IMM) ) THEN
               WRITE (NOUT1,99001)
               STOP
            END IF
            IF ( ZZ(4,IMM).GT.ZZ(5,IMM) ) THEN
               WRITE (NOUT1,99002)
               STOP
            END IF
C
            DO I = 1,2
               VV(I,ILH,IMM) = V(I,ILH,IMM)*BOHRI*BOHRI*RYDBI
            END DO
            VV(3,ILH,IMM) = V(3,ILH,IMM)*RYDBI
            VV(4,ILH,IMM) = V(4,ILH,IMM)*RYDBI
C
C           find za and zb, limits of the barrier
C
            IF ( IMM.EQ.1 ) THEN
               ZSU = ZZ(1,IMM)
               ZAU = ZZ(2,IMM)
               ZBU = ZZ(2,IMM)
               DO I = 3,5
                  ZAU = MIN(ZAU,ZZ(I,IMM))
                  ZBU = MAX(ZBU,ZZ(I,IMM))
               END DO
               IF ( ZSU-ZAU.LT..01D0 ) ZAU = ZSU - .01D0
            ELSE IF ( IMM.EQ.2 ) THEN
               ZSD = ZZ(1,IMM)
               ZAD = ZZ(2,IMM)
               ZBD = ZZ(2,IMM)
               DO I = 3,5
                  ZAD = MIN(ZAD,ZZ(I,IMM))
                  ZBD = MAX(ZBD,ZZ(I,IMM))
               END DO
               IF ( ZSD-ZAD.LT..01D0 ) ZAD = ZSD - .01D0
            END IF
C
C           polynomial fit for rm potential
C
C            IF ( ZZ(2,IMM).NE.ZZ(3,IMM) ) THEN
            IF ( ABS(ZZ(2,IMM)-ZZ(3,IMM)).GT.1.0D-16 ) THEN
C
C              real part polynomial fit, for z(2,imm).lt.z.lt.z(3,imm)
C
               Z12 = 1./(ZZ(2,IMM)-ZZ(1,IMM))
               Z13 = Z12**2
               Z123 = 0.5D0*Z12 + VV(1,ILH,IMM)*Z13
               Z113 = -0.5D0*Z13 - 2.*VV(1,ILH,IMM)*Z13*Z12
               H1 = ZZ(2,IMM)
               H2 = ZZ(3,IMM)
               H3 = VV(3,ILH,IMM)
               CALL HERPOL(H1,Z123,Z113,H2,H3,0.D0,A,IMM)
            END IF
C
C            IF ( ZZ(4,IMM).NE.ZZ(5,IMM) ) THEN
            IF ( ABS(ZZ(4,IMM)-ZZ(5,IMM)).GT.1.0D-16 ) THEN
C               IF ( V(4,ILH,IMM).NE.0.D0 .OR. V(2,ILH,IMM).NE.0.D0 )
C     &              THEN
               IF ( ABS(V(4,ILH,IMM)).GT.1.0D-16 .OR. ABS(V(2,ILH,IMM))
     &              .GT.1.0D-16 ) THEN
C
C                   imag part polynomial fit, for z(4,imm).lt.z.lt.z(5,imm)
C
                  Z12 = 1.D0/(ZZ(4,IMM)-ZZ(1,IMM))
                  Z13 = VV(2,ILH,IMM)*Z12*Z12
                  Z123 = -2.D0*Z13*Z12
                  H1 = ZZ(4,IMM)
                  H2 = ZZ(5,IMM)
                  H3 = VV(4,ILH,IMM)
                  CALL HERPOL(H1,Z13,Z123,H2,H3,0.D0,B,IMM)
               END IF
            END IF
C
C           end rm potential
C
         ELSE
C
C           some other potential (looks like step)
C
            DO I = 1,5
               ZZ(I,IMM) = Z(I,IMM)*BOHRI
               VV(I,ILH,IMM) = V(I,ILH,IMM)*RYDBI
            END DO
C
C           find za and zb, the limits of the barrier
C
            IF ( IMM.EQ.1 ) THEN
               ZSU = ZZ(1,IMM)
               ZAU = ZZ(2,IMM)
               ZBU = ZZ(2,IMM)
               DO I = 3,5
                  ZAU = MIN(ZAU,ZZ(I,IMM))
                  ZBU = MAX(ZBU,ZZ(I,IMM))
               END DO
               IF ( ZSU-ZAU.LT..01D0 ) ZAU = ZSU - .01D0
            ELSE IF ( IMM.EQ.2 ) THEN
               ZSD = ZZ(1,IMM)
               ZAD = ZZ(2,IMM)
               ZBD = ZZ(2,IMM)
               DO I = 3,5
                  ZAD = MIN(ZAD,ZZ(I,IMM))
                  ZBD = MAX(ZBD,ZZ(I,IMM))
               END DO
               IF ( ZSD-ZAD.LT..01D0 ) ZAD = ZSD - .01D0
            END IF
C
            IF ( ZZ(2,IMM).GT.ZZ(1,IMM) ) THEN
C
C              real part polynomial fits
C
               H1 = ZZ(2,IMM)
               H2 = ZZ(3,IMM)
               H3 = VV(3,ILH,IMM)
               CALL HERPOL(H1,Z123,Z113,H2,H3,0.D0,A,IMM)
            END IF
C
            IF ( ZZ(4,IMM).GT.ZZ(3,IMM) ) THEN
C
C              imaginary part polynomial fits
C
               H1 = ZZ(4,IMM)
               H2 = ZZ(5,IMM)
               H3 = VV(4,ILH,IMM)
               CALL HERPOL(H1,Z13,Z123,H2,H3,0.D0,B,IMM)
            END IF
         END IF
      END DO
C
      RETURN
C
99001 FORMAT (' ***** potin-z(2).gt.z(3)')
99002 FORMAT (' ***** potin-z(4).gt.z(5)')
      END
C*==herpol.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE HERPOL(X0,F0,F01,X1,F1,F11,C,IMM)
C     /****************************************************************/
C     # purpose      : cubic hermite interpolation                     *
C                      in x0.lt.x.lt.x1 f0,f01/f1,f11 are ?????        *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 F0,F01,F1,F11,X0,X1
      INTEGER IMM
      REAL*8 C(4,2)
C
C Local variables
C
      REAL*8 A,B,S,SA
C
C*** End of declarations rewritten by SPAG
C
      S = X1 - X0
      C(1,IMM) = F0
      C(2,IMM) = F01
      B = F1 - F0 - F01*S
      A = 3.D0*B - S*(F11-F01)
      SA = 1.D0/S
      S = SA*SA
      C(3,IMM) = A*S
      C(4,IMM) = (B-A)*S*SA
C
      END
C*==cbpot.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CBPOT(EZ,Z,ILH,CBX,IMM)
C
C     /****************************************************************/
C     # purpose      : cbpot(z,ilh) gives the potential-energy         *
C                      in ry (!) at position z in bohr (!)             *
C                                                                      *
C     # function called from this routine:                             *
C       gauss                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC_COMPOT,ONLY:ZZ,A,B,VV,KTRB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CBX
      REAL*8 EZ,Z
      INTEGER ILH,IMM
C
C Local variables
C
      REAL*8 D,DS,DZS,DZT,VI,VR
C
C*** End of declarations rewritten by SPAG
C
      D = Z
      IF ( KTRB.NE.1 ) THEN
C
C         rm barrier
C         real part
C
         IF ( D.LE.ZZ(2,IMM) ) THEN
            DZS = 1./(D-ZZ(1,IMM))
            DZT = DZS*DZS
            VR = 0.5*DZS
C            IF ( VV(1,ILH,IMM).NE.0. ) VR = VR + VV(1,ILH,IMM)*DZT
            IF ( ABS(VV(1,ILH,IMM)).GT.1.0D-16 ) VR = VR + VV(1,ILH,IMM)
     &           *DZT
         ELSE IF ( D.LE.ZZ(3,IMM) ) THEN
            DS = D - ZZ(2,IMM)
            VR = A(1,IMM) + DS*(A(2,IMM)+DS*(A(3,IMM)+DS*A(4,IMM)))
         ELSE
            VR = VV(3,ILH,IMM)
         END IF
C         imaginary part
         IF ( D.LE.ZZ(4,IMM) ) THEN
C            IF ( VV(2,ILH,IMM).EQ.0. ) THEN
            IF ( ABS(VV(2,ILH,IMM)).LE.1.0D-16 ) THEN
               VI = 0.0
            ELSE IF ( D.GT.ZZ(2,IMM) ) THEN
               DZT = D - ZZ(1,IMM)
               DZT = 1./(DZT*DZT)
            ELSE
               VI = VV(2,ILH,IMM)*DZT
            END IF
         ELSE IF ( D.LE.ZZ(5,IMM) ) THEN
C            IF ( VV(2,ILH,IMM).EQ.0. .AND. VV(4,ILH,IMM).EQ.0. ) THEN
            IF ( ABS(VV(2,ILH,IMM)).LE.1.0D-16 .AND. ABS(VV(4,ILH,IMM))
     &           .LE.1.0D-16 ) THEN
               VI = 0.0
            ELSE
               DS = D - ZZ(4,IMM)
               VI = B(1,IMM) + DS*(B(2,IMM)+DS*(B(3,IMM)+DS*B(4,IMM)))
            END IF
         ELSE
            VI = VV(4,ILH,IMM)
         END IF
      ELSE
C         step
C         real part
         IF ( D.GE.ZZ(2,IMM) ) THEN
            VR = VV(3,ILH,IMM)
         ELSE IF ( D.GT.ZZ(1,IMM) ) THEN
            DS = D - ZZ(1,IMM)
            VR = A(1,IMM) + DS*(A(2,IMM)+DS*(A(3,IMM)+DS*A(4,IMM)))
         ELSE
            VR = VV(1,ILH,IMM)
         END IF
C         imaginary part
         IF ( D.GE.ZZ(4,IMM) ) THEN
            VI = VV(4,ILH,IMM)
         ELSE IF ( D.GT.ZZ(3,IMM) ) THEN
            DS = D - ZZ(3,IMM)
            VI = B(1,IMM) + DS*(B(2,IMM)+DS*(B(3,IMM)+DS*B(4,IMM)))
         ELSE
            VI = VV(2,ILH,IMM)
         END IF
      END IF
      CBX = DCMPLX(VR-EZ,VI)
C
Ccghf add a gaussian barrier
C      if (ktrb .eq. 2 .and. ilh.eq.1) then
C          barab = -9.0
C          barwidth = 2.0
C          barheight = 1.0
C          cg = barheight * gauss(d,barab,barwidth)
C          cbx = dcmplx(vr - ez + cg,vi)
C      else
C          cbx = dcmplx(vr - ez, vi)
C      endif
C
      END
C*==rsprop.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE RSPROP(PS,GANZ,QD,VPR,EMV,EMV1,ILH,IREL,CLIGHT,BK,ADU,
     &                  ADD,EPSX,TPM,TPP,TMP,TMM,IZBG,EZRY,ZAU,ZSU,ZBU,
     &                  ZAD,ZSD,ZBD,IBLOCH)
C     /****************************************************************/
C     # purpose      : calculate reflection and transmission           *
C                      coefficients for the surface potential          *
C                                                                      *
C     # subroutines called from this routine:                          *
C       reflrm                                                         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,TRANSP_BAR,CONE,CZERO
      USE MOD_SPEC_GEOM,ONLY:SEP
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 CLIGHT,EPSX,EZRY,VPR,ZAD,ZAU,ZBD,ZBU,ZSD,ZSU
      COMPLEX*16 EMV,EMV1
      INTEGER GANZ,IBLOCH,ILH,IREL,IZBG,QD
      REAL*8 ADD(3),ADU(3)
      COMPLEX*16 BK(LH,2),PS(LH,2),TMM(LH),TMP(LH),TPM(LH),TPP(LH)
C
C Local variables
C
      REAL*8 A,AE,B,C,RE,ZA,ZB,ZS
      COMPLEX*16 BX,CRT(4,LH),E,EH,EK,EL,EM
      INTEGER I,IDN,IMM,IN,IS,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,LH
         PS(I,1) = CZERO
         PS(I,2) = CZERO
         TPP(I) = CZERO
         TPM(I) = CZERO
         TMP(I) = CZERO
         TMM(I) = CZERO
C
         DO J = 1,4
            CRT(J,I) = CZERO
         END DO
      END DO
C
      IF ( ILH.EQ.1 ) THEN
         E = EMV1
         IZBG = 1
      ELSE IF ( ILH.EQ.2 ) THEN
         E = EMV
         IZBG = 1
      END IF
C
      IF ( IBLOCH.EQ.1 ) THEN
         E = EMV
         IZBG = GANZ
      END IF
C
C     izbg controls the number of reciprocal lattice vectors, for
C     which the rm-barrier potential is used
C
      IF ( IZBG.GT.GANZ ) THEN
         WRITE (NOUT1,99002)
         STOP
      END IF
C
C     ad = distance of spin-dependent surface barrier from the first bulk layer
C
      ADU(1) = SEP(1,0)
      ADU(2) = 0.0
      ADU(3) = 0.0
C
      ADD(1) = SEP(2,0)
      ADD(2) = 0.0
      ADD(3) = 0.0
C
      DO IMM = 1,2
         EH = E
         EL = E - DCMPLX(VPR,0.D0)
         DO IS = 1,IZBG
            IF ( IMM.EQ.1 ) IDN = IS
            IF ( IMM.EQ.2 ) IDN = IS + GANZ
            BX = BK(IDN,1)**2 + BK(IDN,2)**2
            IF ( IREL.EQ.2 ) THEN
               EM = CDSQRT(2.D0*EH+(EH/CLIGHT)**2-BX)
               EK = CDSQRT(2.D0*EL+(EL/CLIGHT)**2-BX)
            ELSE
               EM = CDSQRT(2.D0*EH-BX)
               EK = CDSQRT(2.D0*EL-BX)
            END IF
            RE = DBLE(E)
            IF ( IREL.EQ.2 ) THEN
               AE = 2.D0*RE - DBLE(BX) + (RE/CLIGHT)**2
            ELSE
               AE = 2.D0*RE - DBLE(BX)
            END IF
            EZRY = AE - VPR*2.D0
            IF ( IMM.EQ.1 ) THEN
               ZA = ZAU
               ZS = ZSU
               ZB = ZBU
            ELSE IF ( IMM.EQ.2 ) THEN
               ZA = ZAD
               ZS = ZSD
               ZB = ZBD
            END IF
            CALL REFLRM(EZRY,EPSX,ILH,CRT,IDN,ZA,ZS,ZB,IMM)
            IF ( DABS(EZRY).GT.20.0D0 .AND. IBLOCH.NE.1 ) THEN
               CRT(3,IDN) = (EK-EM)/(EK+EM)
               CRT(2,IDN) = 2.D0*EM/(EK+EM)
               CRT(1,IDN) = (EK-EM)/(EK+EM)
               CRT(4,IDN) = 2.D0*EM/(EK+EM)
            END IF
         END DO
C
         IF ( IMM.EQ.1 ) THEN
            DO IS = 1,GANZ
               IF ( IS.LE.IZBG ) THEN
                  TPM(IS) = CRT(3,IS)
                  TPP(IS) = CRT(2,IS)
                  TMP(IS) = CRT(1,IS)
                  TMM(IS) = CRT(4,IS)
               ELSE
                  BX = BK(IS,1)*BK(IS,1) + BK(IS,2)*BK(IS,2)
                  EH = E - DCMPLX(VPR,0.D0)
                  EL = E
                  IF ( IREL.EQ.2 ) THEN
                     EM = CDSQRT(2.D0*EH+(EH/CLIGHT)**2-BX)
                     EK = CDSQRT(2.D0*EL+(EL/CLIGHT)**2-BX)
                  ELSE
                     EM = CDSQRT(2.D0*EH-BX)
                     EK = CDSQRT(2.D0*EL-BX)
                  END IF
                  TMP(IS) = (EK-EM)/(EK+EM)
                  TPP(IS) = 2.D0*EM/(EK+EM)
                  TPM(IS) = (EK-EM)/(EK+EM)
                  TMM(IS) = 2.D0*EM/(EK+EM)
               END IF
            END DO
         ELSE IF ( IMM.EQ.2 ) THEN
            DO IS = GANZ + 1,QD
               IF ( IS.LE.GANZ+IZBG ) THEN
                  TPM(IS) = CRT(3,IS)
                  TPP(IS) = CRT(2,IS)
                  TMP(IS) = CRT(1,IS)
                  TMM(IS) = CRT(4,IS)
               ELSE
                  BX = BK(IS,1)*BK(IS,1) + BK(IS,2)*BK(IS,2)
                  EH = E - DCMPLX(VPR,0.D0)
                  EL = E
                  IF ( IREL.EQ.2 ) THEN
                     EM = CDSQRT(2.D0*EH+(EH/CLIGHT)**2-BX)
                     EK = CDSQRT(2.D0*EL+(EL/CLIGHT)**2-BX)
                  ELSE
                     EM = CDSQRT(2.D0*EH-BX)
                     EK = CDSQRT(2.D0*EL-BX)
                  END IF
                  TMP(IS) = (EK-EM)/(EK+EM)
                  TPP(IS) = 2.D0*EM/(EK+EM)
                  TPM(IS) = (EK-EM)/(EK+EM)
                  TMM(IS) = 2.D0*EM/(EK+EM)
               END IF
            END DO
         END IF
      END DO
C
C     transparent barrier potential
C
      IF ( TRANSP_BAR ) THEN
         IF ( ILH.EQ.2 ) THEN
            DO IS = 1,QD
               TMP(IS) = CZERO
               TPP(IS) = CONE
               TPM(IS) = CZERO
               TMM(IS) = CONE
            END DO
         END IF
      END IF
C
      DO IMM = 1,2
         IF ( IMM.EQ.1 ) THEN
            DO IS = 1,GANZ
               BX = BK(IS,1)*BK(IS,1) + BK(IS,2)*BK(IS,2)
               EH = E - DCMPLX(VPR,0.D0)
               IF ( IREL.EQ.2 ) THEN
                  EK = CDSQRT(2.D0*E+(E/CLIGHT)**2-BX)
               ELSE
                  EK = CDSQRT(2.D0*E-BX)
               END IF
               A = DBLE(BK(IS,1))*ADU(2) + DBLE(BK(IS,2))*ADU(3)
               B = -DIMAG(EK)*ADU(1)
               C = DBLE(EK)*ADU(1)
               PS(IS,1) = CDEXP(DCMPLX(B,C+A))
               PS(IS,2) = CDEXP(DCMPLX(B,C-A))
            END DO
         ELSE IF ( IMM.EQ.2 ) THEN
            DO IS = 1 + GANZ,QD
               BX = BK(IS,1)*BK(IS,1) + BK(IS,2)*BK(IS,2)
               EH = E - DCMPLX(VPR,0.D0)
               IF ( IREL.EQ.2 ) THEN
                  EK = CDSQRT(2.D0*E+(E/CLIGHT)**2-BX)
               ELSE
                  EK = CDSQRT(2.D0*E-BX)
               END IF
               A = DBLE(BK(IS,1))*ADD(2) + DBLE(BK(IS,2))*ADD(3)
               B = -DIMAG(EK)*ADD(1)
               C = DBLE(EK)*ADD(1)
               PS(IS,1) = CDEXP(DCMPLX(B,C+A))
               PS(IS,2) = CDEXP(DCMPLX(B,C-A))
            END DO
         END IF
      END DO
C
      IF ( IP.GE.3 ) THEN
         WRITE (NOUT1,99004)
         WRITE (NOUT1,99003) (TPM(J),J=1,QD)
         WRITE (NOUT1,99003) (TPP(J),J=1,QD)
         WRITE (NOUT1,99003) (TMP(J),J=1,QD)
         WRITE (NOUT1,99003) (TMM(J),J=1,QD)
         IN = 0
         WRITE (NOUT1,99001) IN,(ADU(J),J=1,3)
         WRITE (NOUT1,99001) IN,(ADD(J),J=1,3)
         WRITE (NOUT1,99003) (PS(J,1),J=1,QD)
         WRITE (NOUT1,99003) (PS(J,2),J=1,QD)
      END IF
      RETURN
C
99001 FORMAT (1x,'rsprop: layer',i3,15x,'displ= ',3F12.5)
99002 FORMAT (1x,'error in rsprop(izbg.gt.ganz)')
99003 FORMAT (5(2x,e12.5,1x,e12.5))
99004 FORMAT (1x,'rsprop: r+-,r++,r-+,r--')
      END
C*==reflrm.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE REFLRM(ENERGY,EPSX,ILH,CRT,IS,ZA,ZS,ZB,IMM)
C     /****************************************************************/
C     # purpose      :  calculate reflection and transmission          *
C                       coefficients                                   *
C                 crt(1)       reflection vacuum to vacuum     r-+     *
C                 crt(2)       transmission vacuum to metal    t++     *
C                 crt(3)       reflection metal to metal       r+-     *
C                 crt(4)       transmission metal to vacuum    t--     *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       cbpot     rkg       cdwc                                       *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,PI,CZERO,CONE,CIMAG
      USE MOD_SPEC_COMPOT,ONLY:VV,KTRB
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 ENERGY,EPSX,ZA,ZB,ZS
      INTEGER ILH,IMM,IS
      COMPLEX*16 CRT(4,LH)
C
C Local variables
C
      REAL*8 AK,AS,ERR,EZ,FN,X
      COMPLEX*16 C,CBX,CD,CFM,CFM1,CFN,CFP,CFP1,CFPFM,CGM,CGM1,CGP,CGP1,
     &           CGPGM,CIQ,CIS,CM,CM1,CP,CP1,CQ,CT,CX
C
C*** End of declarations rewritten by SPAG
C
      CGP = CZERO
      CGP1 = CZERO
      CGM = CZERO
      CGM1 = CZERO
C
C     g-fundamental solutions at zb
C
      EZ = ENERGY
      CALL CBPOT(EZ,ZB,ILH,CBX,IMM)
      CQ = CDSQRT((-CBX))
      CIQ = CIMAG*CQ
      C = CDEXP(CIQ*ZB)
      CP = C
      CM = CONE/C
      CP1 = CIQ*CP
      CM1 = -CIQ*CM
C
C     write (nout1,*) 'test', ' is=',is, ' ktrb=',ktrb
C     write (nout1,*) 'cp=',cp,'cm=',cm
C
C     integrate through the barrier from zb to za. cgp,cgm are
C     wavefunctions at za, cgp1,cgm1 derivatives.
C
      CALL RKG(EZ,ZB,CP,CP1,ZA,CGP,CGP1,20,EPSX,ERR,ILH,IMM)
      IF ( ERR.GE.EPSX ) THEN
         WRITE (NOUT1,99001) ENERGY
         STOP
      ELSE
         CALL RKG(EZ,ZB,CM,CM1,ZA,CGM,CGM1,20,EPSX,ERR,ILH,IMM)
         IF ( ERR.GT.EPSX ) THEN
            WRITE (NOUT1,99002) ENERGY
            STOP
         END IF
C
C         find f-fundamental solutions at zb
C
         IF ( KTRB.EQ.1 ) THEN
C
C             cfp,cfm are normalized to 1 at z=0
C
            CALL CBPOT(EZ,ZA,ILH,CBX,IMM)
            CQ = CDSQRT((-CBX))
            CIQ = CIMAG*CQ
            CFP = CDEXP(CIQ*ZA)
C             ca = cone
            CFM = CONE/CFP
            CFP1 = CIQ*CFP
            CFM1 = -CIQ*CFM
         ELSE
C
C             find parameters and arguments of the whittaker functions
C
            AK = 2.D0*SQRT(ABS(ENERGY))
            AS = 1.D0/(2.D0*AK)
            FN = PI*AS/2.D0
            X = AK*(ZS-ZA)
            IF ( ENERGY.LE.0.D0 ) THEN
               CIS = DCMPLX((-AS),0.D0)
               CX = DCMPLX((-X),1.D-38)
               CD = DCMPLX(AK,0.D0)
               CFN = DCMPLX(COS(FN),SIN(FN))
            ELSE
               CIS = DCMPLX(0.D0,(-AS))
               CX = DCMPLX(0.D0,X)
               CD = DCMPLX(0.D0,(-AK))
               CFN = DCMPLX(EXP((-FN)),0.D0)
            END IF
C
C     evaluation of the whittaker functions
C
C            IF ( VV(1,ILH,IMM).EQ.0.D0 .AND. VV(2,ILH,IMM).EQ.0.D0 )
C     &           THEN
            IF ( ABS(VV(1,ILH,IMM)).LE.1.0D-16 .AND. ABS(VV(2,ILH,IMM))
     &           .LE.1.0D-16 ) THEN
               CT = DCMPLX(0.5D0,0.0D0)
            ELSE
               CT = CDSQRT(DCMPLX(0.25D0+VV(1,ILH,IMM),VV(2,ILH,IMM)))
            END IF
            CALL CDWC(CIS,CT,CX,CFP,CFP1,EPSX)
C
            CFP = CFN*CFP
            CFP1 = CFN*CFP1*CD
C            IF ( ENERGY.LE.0.D0 .OR. VV(2,ILH,IMM).NE.0.D0 ) THEN
            IF ( ENERGY.LE.0.D0 .OR. ABS(VV(2,ILH,IMM)).GT.1.0D-16 )
     &           THEN
               CALL CDWC((-CIS),CT,(-CX),CFM,CFM1,EPSX)
               CFM = CFN*CFM
               CFM1 = CFN*CFM1*(-CD)
            ELSE
               CFM = DCONJG(CFP)
               CFM1 = DCONJG(CFP1)
            END IF
         END IF
         CGPGM = CGP*CGM1 - CGM*CGP1
         CFPFM = CFP*CFM1 - CFM*CFP1
         C = CONE/(CFM*CGP1-CGP*CFM1)
         CRT(1,IS) = -(CFP*CGP1-CGP*CFP1)*C
         CRT(2,IS) = -CFPFM*C
         CRT(3,IS) = -(CFM*CGM1-CGM*CFM1)*C
         CRT(4,IS) = -CGPGM*C
      END IF
C
      RETURN
C
99001 FORMAT (' ***** reflex-energy  ',e15.6,' ry',2x,
     &        'integration of g+ could not be done to desired accuracy')
99002 FORMAT (' ***** reflex-energy  ',e15.6,' ry',2x,
     &        'integration of g- could not be done to desired accuracy')
      END
C*==rkg.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE RKG(EZ,DX1,CY1,CZ1,DX2,CY2,CZ2,NSTEP,EPS,ERR,ILH,IMM)
C     /****************************************************************/
C     # purpose      : integration of equation cy''(x)=cy(x)*cbpot(x)  *
C                      from x=x1 to x=x2.                              *
C     # method       : runge-kutta integration as modified by gill.    *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:CZERO
      USE MOD_FILES,ONLY:RDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CY1,CY2,CZ1,CZ2
      REAL*8 DX1,DX2,EPS,ERR,EZ
      INTEGER ILH,IMM,NSTEP
C
C Local variables
C
      COMPLEX*16 C,CBX,CCY1,CCY2,CCZ1,CCZ2,CK(3),CQ(3),CY(3)
      REAL*8 DCY,DD,DERP,DH,ERP,R1,R2,R3,R4,U1,U2
      INTEGER I,IND,INDMAX,J,N
C
C*** End of declarations rewritten by SPAG
C
                                                  !,CCY2,CCZ2
C
C*** End of declarations rewritten by SPAG
C
      R1 = 0.5D0                      !,U1,U2,DERP
      R2 = 0.29289321881345D0
      R3 = 1.70710678118650D0
      R4 = 0.166666666666667D0
      U1 = 0.0D0
      U2 = 0.0D0
      DERP = 0.0D0
      RDUMMY = DERP
      CCY2 = CZERO
      CCZ2 = CZERO
      INDMAX = 10
      N = NSTEP
      IND = 0
      CCY1 = CY1
      CCZ1 = CZ1
C      IF ( DX1.EQ.DX2 ) THEN
      IF ( ABS(DX1-DX2).LE.1.0D-16 ) THEN
         CY2 = CY1
         CZ2 = CZ1
         ERR = 0.D0
      ELSE
         DD = (DBLE(DX2)-DBLE(DX1))/DBLE(N)
         DH = DD
 50      CONTINUE
         CY(1) = DCMPLX(DX1,0.D0)
         CY(2) = CCY1
         CY(3) = CCZ1
         CK(1) = DCMPLX(DD,0.D0)
C
         DO I = 1,3
            CQ(I) = CZERO
         END DO
C
         DO J = 1,N
            IF ( J.GE.N ) THEN
               DD = DX2 - DBLE(CY(1))
               CK(1) = DCMPLX(DD,0.D0)
            END IF
            CK(2) = DD*CY(3)
            DCY = DBLE(CY(1))
            CALL CBPOT(EZ,DCY,ILH,CBX,IMM)
            CK(3) = DD*CY(2)*CBX
            DO I = 1,3
               C = (CK(I)-2.D0*CQ(I))*R1
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R1
            END DO
            CK(2) = DD*CY(3)
            DCY = DBLE(CY(1))
            CALL CBPOT(EZ,DCY,ILH,CBX,IMM)
            CK(3) = DD*CY(2)*CBX
            DO I = 1,3
               C = R2*(CK(I)-CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R2
            END DO
            CK(2) = DD*CY(3)
            DCY = DBLE(CY(1))
            CALL CBPOT(EZ,DCY,ILH,CBX,IMM)
            CK(3) = DD*CY(2)*CBX
            DO I = 1,3
               C = R3*(CK(I)-CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R3
            END DO
            CK(2) = DD*CY(3)
            DCY = DBLE(CY(1))
            CALL CBPOT(EZ,DCY,ILH,CBX,IMM)
            CK(3) = DD*CY(2)*CBX
            DO I = 1,3
               C = R4*(CK(I)-2.D0*CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R1
            END DO
         END DO
         DD = DH
         IND = IND + 1
         IF ( IND.GT.INDMAX ) STOP
         ERP = 1.D0
         IF ( IND.GT.1 ) THEN
            U1 = CDABS((CY(2)-CCY2)/CCY2)
            U2 = CDABS((CY(3)-CCZ2)/CCZ2)
            DERP = (U1+U2)*0.5D0
            ERP = DERP
         END IF
         IF ( ERP.LT.EPS ) THEN
            CY2 = CY(2)
            CZ2 = CY(3)
            ERR = ERP
         ELSE
            DD = DD*0.5D0
            DH = DD
            N = N*2
            CCY2 = CY(2)
            CCZ2 = CY(3)
            GOTO 50
         END IF
      END IF
C
      END
C*==cdwc.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CDWC(CK,CM,CZ,CW,CW1,EPS)
C     /****************************************************************/
C     # purpose      : whittaker w function of complex parameters      *
C                      ck and cm, and argument cz.                     *
C                      cw=function value, cw1=value of derivative      *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       cdwha     cdwhb     cdwhc     crkg    cvcp                     *
C     /****************************************************************/
C
      USE MOD_SPEC_COMCP,ONLY:CKC,CMC
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CK,CM,CW,CW1,CZ
      REAL*8 EPS
C
C Local variables
C
      REAL*8 A,ALFA,B,ERQ,ERR,X,XL
      COMPLEX*16 CA,CB,CWW,CWW1,CZZ
C
C*** End of declarations rewritten by SPAG
C
      ALFA = 3.D0
      CKC = CK
      CMC = CM
      X = CDABS(CZ)
      CA = 0.5D0 + CM - CK
      A = CDABS(CA)
      CB = 1.D0 + 2.D0*CM
      B = CDABS(CB)
      XL = ALFA*B/A
      IF ( X.LE.XL ) THEN
C         evaluation by convergent series
C         IF ( DBLE(CM).NE.0.5D0 .AND. DIMAG(CM).NE.0.D0 )
C     &        CALL CDWHA(CK,CM,CZ,CW,CW1,EPS)
         IF ( ABS(DBLE(CM)-0.5D0).GT.1.0D-16 .AND. ABS(DIMAG(CM))
     &        .GT.1.0D-16 ) CALL CDWHA(CK,CM,CZ,CW,CW1,EPS)
Ccghf ???? seems the if's are not working correctly if cm is calculated
C         from root(0.25) check whether an else construct is possible !!!!!
C         else
C         IF ( DBLE(CM).EQ.0.5D0 .AND. DIMAG(CM).EQ.0.D0 )
C     &        CALL CDWHB(CW,CW1,CK,CZ)
         IF ( ABS(DBLE(CM)-0.5D0).LE.1.0D-16 .AND. ABS(DIMAG(CM))
     &        .LE.1.0D-16 ) CALL CDWHB(CW,CW1,CK,CZ)
         RETURN
      END IF
C
      IF ( SQRT(X)*EXP((-X)).LE.EPS ) THEN
         CALL CDWHC(CK,CM,CZ,CW,CW1,EPS,ERR)
         IF ( ERR.LE.EPS ) RETURN
      END IF
C
      CZZ = CZ*XL/X
C      IF ( DBLE(CM).NE.0.5D0 .AND. DIMAG(CM).NE.0.D0 )
C     &     CALL CDWHA(CK,CM,CZZ,CW,CW1,EPS)
C      IF ( DBLE(CM).EQ.0.5D0 .AND. DIMAG(CM).EQ.0.D0 )
C     &     CALL CDWHB(CW,CW1,CK,CZZ)
      IF ( ABS(DBLE(CM)-0.5D0).GT.1.0D-16 .AND. ABS(DIMAG(CM))
     &     .GT.1.0D-16 ) CALL CDWHA(CK,CM,CZZ,CW,CW1,EPS)
      IF ( ABS(DBLE(CM)-0.5D0).LE.1.0D-16 .AND. ABS(DIMAG(CM))
     &     .LE.1.0D-16 ) CALL CDWHB(CW,CW1,CK,CZZ)
C
      CALL CRKG(CZZ,CW,CW1,CZ,CWW,CWW1,20,EPS,ERQ)         !, cvcp)
      IF ( ERQ.GT.EPS .AND. IP.GT.(-1) ) WRITE (NOUT1,99001)
      CW = CWW
      CW1 = CWW1
      RETURN
C
99001 FORMAT ('** cdwc-warning desired accuracy in crkg not attained')
      END
C*==cdwha.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CDWHA(CK,CM,CZ,CW,CW1,EPS)
C     /****************************************************************/
C     # purpose      : whittaker w function of complex parameters      *
C                      ck and cm, and argument cz. cw=function value,  *
C                      cw1=value of derivative of function.            *
C                      2*cm must not be an integer, positive,          *
C                      negative or zero !                              *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       cdgam                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:PI,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CK,CM,CW,CW1,CZ
      REAL*8 EPS
C
C Local variables
C
      COMPLEX*16 CA,CAA,CB,CBB,CL,CMA,CMAA,CMB,CMBB,CS1,CS2,CS3,CS4,CT1,
     &           CT2,CXA,CXX,CXZ
      COMPLEX*16 CDGAM
      REAL*8 DI,EPSS
      INTEGER I
      EXTERNAL CDGAM
C
C*** End of declarations rewritten by SPAG
C
C
      EPSS = EPS*1.D-2
      CA = 0.5D0 + CM - CK
      CB = 1.D0 + 2.D0*CM
      CAA = 0.5D0 - CM - CK
      CBB = 1.D0 - 2.D0*CM
      CS1 = CONE
      CS3 = CONE
      CS2 = CZERO
      CS4 = CZERO
      CT1 = CA/CB
      CT2 = CAA/CBB
      DO I = 1,40
         CS2 = CS2 + CT1
         CS4 = CS4 + CT2
         DI = DBLE(I)
         CL = DCMPLX(DI,0.D0)
         CXZ = CZ*DCMPLX(1.D0/DI,0.D0)
         CT1 = CT1*CXZ
         CT2 = CT2*CXZ
         CS1 = CS1 + CT1
         CS3 = CS3 + CT2
         IF ( CDABS(CT1).LT.EPSS ) EXIT
         CT1 = CT1*(CA+CL)/(CB+CL)
         CT2 = CT2*(CAA+CL)/(CBB+CL)
      END DO
C
      CXA = 1.D0/CZ
      CL = CDLOG(CZ)
      CXX = CDEXP((0.5D0+CM)*CL-CZ/2.D0)
      CMA = CS1*CXX
      CMAA = ((0.5D0+CM)*CXA-0.5D0)*CMA + CXX*CS2
      CXX = CDEXP((0.5D0-CM)*CL-CZ/2.D0)
      CMB = CS3*CXX
      CMBB = ((0.5D0-CM)*CXA-0.5D0)*CMB + CXX*CS4
      CT1 = 2.D0*CM
      CT2 = CDGAM(CT1)
      CS2 = CT2/CDGAM(CA)
      CS1 = -PI/(CT1*CT2*CDSIN(PI*CT1)*CDGAM(CAA))
      CW = CMA*CS1 + CMB*CS2
      CW1 = CMAA*CS1 + CMBB*CS2
C
      END
C*==cdwhb.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CDWHB(CW,CW1,CK,CZ)
C     /****************************************************************/
C     # purpose      : whittaker w function of parameters              *
C                      ck, 1/2 and argument cz. value function = cw,   *
C                      value of derivative of function = cw1.          *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       cdgam     cdpsi                                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:PI,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CK,CW,CW1,CZ
C
C Local variables
C
      COMPLEX*16 C,CAB,CABQ,CAX,CC,CG1,CGK,CKZ,CLG,CM,CM1,CPIK,CSM1,
     &           CSM2,CSUM1,CSUM2,CV1(40),CXP
      COMPLEX*16 CDGAM,CDPSI
      REAL*8 DPSIR(40),DT,FLI
      INTEGER I,JJ,N,NJ,NN,NTERM
      EXTERNAL CDGAM,CDPSI
C
C*** End of declarations rewritten by SPAG
C
C     dpsir(i) contains psi(i)+psi(1+i) for i=1,40 where psi is the
C     digamma function
C
      DATA DPSIR/ - 1.5443132980307D-01,1.3455686701969D0,
     &     2.1789020035303D0,2.7622353368636D0,3.2122353368636D0,
     &     3.5789020035303D0,3.8884258130541D0,4.1562829559112D0,
     &     4.3923067022300D0,4.6035051781334D0,4.7414269042500D0,
     &     4.9686566932850D0,5.1289131035414D0,5.2772647518930D0,
     &     5.4153599899883D0,5.5445266566549D0,5.6658501860667D0,
     &     5.7802292710340D0,5.8884164055369D0,5.9910479844843D0,
     &     6.0886670321033D0,6.1817406251769D0,6.2706734315010D0,
     &     6.3558183590373D0,6.4374850257039D0,6.5156564165500D0,
     &     6.5914451396641D0,6.6641964624154D0,6.7343935067504D0,
     &     6.8022095987044D0,6.8678009965538D0,6.9313090610700D0,
     &     6.9928620913730D0,7.0525768863819D0,7.1105600796592D0,
     &     7.1669092860084D0,7.2217140908132D0,7.2750569073140D0,
     &     7.3270137224287D0,7.3776547480697D0/
C
      NTERM = 16
C
      C = CK
C      IF ( DIMAG(C).NE.0.D0 .OR. DBLE(C).LE.0.D0 ) THEN
      IF ( ABS(DIMAG(C)).GT.1.0D-16 .OR. DBLE(C).LE.0.D0 ) THEN
         CGK = 1.D0/CDGAM((-CK))
      ELSE
C         N = DBLE(C)
         N = INT(DBLE(C))
         IF ( DBLE(N).NE.DBLE(C) ) THEN
            CGK = 1.D0/CDGAM((-CK))
         ELSE
            CGK = CZERO
         END IF
      END IF
      NN = NTERM + 1
      CPIK = CK*PI
C      CC = CONE - CK
      CG1 = CDGAM(CONE+CK)
      DO I = 1,NN
         IF ( DBLE(C).LE.0.5D0 ) THEN
            CV1(I) = CGK*CDPSI(CONE-C)
         ELSE
            CV1(I) = -CG1*(CDCOS(CPIK)+CDSIN(CPIK)*CDPSI(C)/PI)
         END IF
         C = C - CONE
      END DO
      CAX = -CK
      CC = CONE
      DT = 1.D0
      DO I = 1,NN
         CV1(I) = CC*(CV1(I)-DPSIR(I)*CGK)/DT
         FLI = DBLE(I)
         CC = CC*(CAX+DCMPLX(FLI,0.D0))
         DT = DT*FLI*(FLI+1.D0)
      END DO
      CSUM1 = CZERO
      CSUM2 = CSUM1
      NJ = NN + 2
      DO I = 2,NN
         JJ = NJ - I
         CSUM1 = CSUM1*CZ + CV1(JJ)
         CSUM2 = CSUM2*CZ + CV1(JJ)*DCMPLX(DBLE(JJ-1),0.D0)
      END DO
      CSUM1 = CSUM1*CZ + CV1(1)
      CKZ = CGK/(CK*CZ)
      CSUM1 = CSUM1 - CKZ
      CSUM2 = CSUM2 + CKZ/CZ
      CSM1 = CONE
      CSM2 = CZERO
      CABQ = CSM1 - CK
      CAB = CSM1
      DO I = 1,40
         FLI = DBLE(I)
         CAB = CAB*CABQ*DCMPLX(1.D0/(FLI+1.D0),0.D0)
         CSM2 = CSM2 + CAB
         CAB = CAB*CZ/FLI
         CSM1 = CSM1 + CAB
         IF ( CDABS(CAB).LT.1.D-10 ) EXIT
         CABQ = CABQ + CONE
      END DO
C
      CXP = CDEXP((-CZ/2.D0))
      CM = CZ*CXP*CSM1
      CM1 = CXP*((CONE-CZ/2.D0)*CSM1+CZ*CSM2)
      CLG = CDLOG(CZ)
      CW = CM*CLG*CGK + CZ*CXP*CSUM1
      CW1 = CGK*(CM1*CLG+CM/CZ) + CXP*((1.D0-CZ/2.D0)*CSUM1+CZ*CSUM2)
C
      END
C*==cdwhc.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CDWHC(CK,CM,CZ,CW,CW1,EPS,ERR)
C     /****************************************************************/
C     # purpose      :                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CK,CM,CW,CW1,CZ
      REAL*8 EPS,ERR
C
C Local variables
C
      COMPLEX*16 CA,CAN,CB,CBN,CS,CS1,CS2,CT,CX
      REAL*8 E1,E2,EEPS,FN
      INTEGER N,NT
C
C*** End of declarations rewritten by SPAG
C
      NT = 25
      CA = 0.5D0 + CM - CK
      CB = 0.5D0 - CM - CK
      CS = -1.D0/CZ
      CT = CA*CB
      EEPS = EPS*CDABS(CT)
      CS2 = CT
      CT = CT*CS
      CS1 = CONE + CT
      CAN = CA
      CBN = CB
      FN = 1.D0
      E2 = 1.D10
      N = 1
 100  CONTINUE
      N = N + 1
      CAN = CAN + CONE
      CBN = CBN + CONE
      CT = CT*CAN*CBN
      E1 = CDABS(CT)
      IF ( E1.LE.E2 ) THEN
         IF ( E1.GE.EEPS ) THEN
            IF ( N.LT.NT ) THEN
               E2 = E1
               CS2 = CS2 + CT
               FN = FN + 1.D0
               CT = CT*CS*DCMPLX(1.D0/FN,0.D0)
               CS1 = CS1 + CT
               GOTO 100
            END IF
         END IF
      END IF
      CX = CDEXP((-CZ/2.D0)+CK*CDLOG(CZ))
      CW = CX*CS1
      CW1 = (CS2*CS*CS-CS1*(0.5D0+CK*CS))*CX
      ERR = E1
      END
C*==crkg.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CRKG(CX1,CY1,CZ1,CX2,CY2,CZ2,NSTEP,EPS,ERR)
C     /****************************************************************/
C     # purpose      : integration of equation cy''(x)=cy(x)*cfct(x)   *
C                      from cx1 to cx2                                 *
C     # method       : runge-kutta integration as modified by gill.    *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:CZERO
      USE MOD_FILES,ONLY:RDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CX1,CX2,CY1,CY2,CZ1,CZ2
      REAL*8 EPS,ERR
      INTEGER NSTEP
C
C Local variables
C
      COMPLEX*16 C,CCD,CCH,CCX1,CCX2,CCY1,CCY2,CCZ1,CCZ2,CK(3),CQ(3),
     &           CY(3),DCY
      COMPLEX*16 CVCP
      REAL*8 DERP,ERP,R1,R2,R3,R4,U1,U2
      INTEGER I,IND,INDMAX,J,N
      EXTERNAL CVCP
C
C*** End of declarations rewritten by SPAG
C
                          !,CCY2,CCZ2
C
C*** End of declarations rewritten by SPAG
C
C     /* from cdwc to cvcp */
C     /* function */
C
      R1 = 0.5D0            !,U1,U2
      R2 = 0.29289321881345D0
      R3 = 1.70710678118650D0
      R4 = 1.D0/6.D0
      U1 = 0.0D0
      U2 = 0.0D0
      DERP = 0.0D0
      RDUMMY = DERP
      CCY2 = CZERO
      CCZ2 = CZERO
C
      INDMAX = 10
      N = NSTEP
      IND = 0
C
      CCY1 = CY1
      CCZ1 = CZ1
      IF ( CX1.EQ.CX2 ) THEN
         CY2 = CY1
         CZ2 = CZ1
         ERR = 0.D0
      ELSE
         CCX1 = CX1
         CCX2 = CX2
         CCD = (CCX2-CCX1)/DBLE(N)
         CCH = CCD
 50      CONTINUE
         CY(1) = CCX1
         CY(2) = CCY1
         CY(3) = CCZ1
         CK(1) = CCD
         DO I = 1,3
            CQ(I) = CZERO
         END DO
         DO J = 1,N
            IF ( J.GE.N ) THEN
               CCD = CCX2 - CY(1)
               CK(1) = CCD
            END IF
            CK(2) = CCD*CY(3)
            DCY = CY(1)
            CK(3) = CCD*CY(2)*CVCP(DCY)
            DO I = 1,3
               C = (CK(I)-2.D0*CQ(I))*R1
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R1
            END DO
            CK(2) = CCD*CY(3)
            DCY = CY(1)
            CK(3) = CCD*CY(2)*CVCP(DCY)
            DO I = 1,3
               C = R2*(CK(I)-CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R2
            END DO
            CK(2) = CCD*CY(3)
            DCY = CY(1)
            CK(3) = CCD*CY(2)*CVCP(DCY)
            DO I = 1,3
               C = R3*(CK(I)-CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R3
            END DO
            CK(2) = CCD*CY(3)
            DCY = CY(1)
            CK(3) = CCD*CY(2)*CVCP(DCY)
            DO I = 1,3
               C = R4*(CK(I)-2.D0*CQ(I))
               CY(I) = CY(I) + C
               CQ(I) = CQ(I) + 3.D0*C - CK(I)*R1
            END DO
         END DO
         CCD = CCH
         IND = IND + 1
         IF ( IND.GT.INDMAX ) STOP
         ERP = 1.D0
         IF ( IND.GT.1 ) THEN
            U1 = CDABS((CY(2)-CCY2)/CCY2)
            U2 = CDABS((CY(3)-CCZ2)/CCZ2)
            ERP = (U1+U2)*0.5D0
         END IF
         IF ( ERP.LT.EPS ) THEN
            CY2 = CY(2)
            CZ2 = CY(3)
            ERR = ERP
         ELSE
            CCD = CCD*0.5D0
            CCH = CCD
            N = N*2
            CCY2 = CY(2)
            CCZ2 = CY(3)
            GOTO 50
         END IF
      END IF
C
      END
C*==cvcp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      COMPLEX*16 FUNCTION CVCP(C)
C     /****************************************************************/
C     # purpose      : for use in crkg when integrating                *
C                      the whittaker function                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:CONE
      USE MOD_SPEC_COMCP,ONLY:CKC,CMC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 C
C
C Local variables
C
      COMPLEX*16 CA,CC
C
C*** End of declarations rewritten by SPAG
C
      CA = DCMPLX(0.25D0,0.D0)
      CC = CONE/C
      CVCP = -(((CA-CMC*CMC)*CC+CKC)*CC-CA)
C
      END
C*==cdgam.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      COMPLEX*16 FUNCTION CDGAM(C)
C     /****************************************************************/
C     # purpose      : cdgam(cz) = gamma function of argument cz       *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:PI,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 C
C
C Local variables
C
      COMPLEX*16 CC,CCC,CZ
      REAL*8 DA,DB(9),DP
      INTEGER I,M
C
C*** End of declarations rewritten by SPAG
C
      DATA DP,DB/.91893853320467D0,8.3333333333333D-2,
     &     - 2.7777777777778D-3,7.9365079365079D-4,
     &     - 5.9523809523810D-4,8.4175084175084D-4,
     &     - 1.9175269175269D-3,6.4102564102564D-3, - 2.95506535771D-2,
     &     1.7964437236883D-1/
C
      CDGAM = CONE
      CC = C
      DA = DBLE(CC) - 0.5D0
      IF ( DA.LT.0. ) CC = 1.D0 - CC
      M = INT(12-DBLE(CC))
      IF ( M.GE.1 ) THEN
         DO I = 1,M
            CDGAM = CDGAM*CC
            CC = CC + CONE
         END DO
      END IF
      CZ = CZERO
      CCC = CONE/CC**2
      DO I = 1,9
         CZ = (CZ+DB(10-I))*CCC
      END DO
      CZ = CZ*CC
      CZ = CZ + (CC-0.5D0)*CDLOG(CC) - CC + DP
      CDGAM = CDEXP(CZ)/CDGAM
      IF ( DA.LT.0. ) CDGAM = PI/(CDSIN(PI*C)*CDGAM)
C
      END
C*==cdpsi.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      COMPLEX*16 FUNCTION CDPSI(CZ)
C     /****************************************************************/
C     # purpose      : cdpsi = digamma function of argument cz         *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:PI,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CZ
C
C Local variables
C
      COMPLEX*16 C,CC,CP
      REAL*8 DB(9)
      INTEGER I,M,N
C
C*** End of declarations rewritten by SPAG
C
      DATA DB/8.3333333333333D-2, - 8.3333333333333D-3,
     &     3.9682539682540D-3, - 4.1666666666667D-3,7.5757575757576D-3,
     &     - 2.1092796092796D-2,8.3333333333333D-2,
     &     - 4.4325980392157D-1,3.0539543302701D0/
C
      C = CZ
      CDPSI = CZERO
      IF ( DBLE(C).LT.0. ) THEN
         CP = C*PI
         CDPSI = -PI*CDCOS(CP)/CDSIN(CP)
         C = CONE - C
      END IF
      N = IDINT(DBLE(C))
      M = 12 - N
      IF ( M.GE.1 ) THEN
         CC = C
         DO I = 1,M
            CDPSI = CDPSI - CONE/CC
            CC = CC + CONE
         END DO
         C = C + DCMPLX(DBLE(M),0.D0)
      END IF
      CDPSI = CDPSI + CDLOG(C)
      C = CONE/C
      CDPSI = CDPSI - C/2.D0
      C = C*C
      CC = CZERO
      DO I = 1,9
         CC = (CC+DB(10-I))*C
      END DO
      CDPSI = CDPSI - CC
C
      END
C*==cpsi2g.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CPSI2G(CMG,PSI2G,WMHF,WPHF,WMHG,WPHG,IREL,GANZ,EZRYH,
     &                  QD,WMSTF,WMSTG,ZMESH,IGP,IGVAL,IGPRO,U0GP,U0GM,
     &                  KGZBLK,IPOL)
C     /****************************************************************/
C     # purpose      : calculates the final state wavefield            *
C                      at the first bulk-layer                         *
C                                                                      *
C     # subroutine called from this routine                            *
C       cbpot                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM,MHDIM1,MHDIM2,CZERO,CIMAG
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,IGP,IPOL,IREL,QD
      COMPLEX*16 CMG(LH),KGZBLK(LH),PSI2G(MHDIM2,LH,4),U0GM(LH),U0GP(LH)
     &           ,WMHF(MHDIM1,LH),WMHG(MHDIM1,LH),WMSTF(6,LH),
     &           WMSTG(6,LH),WPHF(MHDIM1,LH),WPHG(MHDIM1,LH)
      REAL*8 EZRYH(LH),ZMESH(6)
      INTEGER IGPRO(LH),IGVAL(LH)
C
C Local variables
C
      COMPLEX*16 DELTA,DELTAS,F1(7),F1S(7),F2(7),F2S(7),V
      REAL*8 EZ,HB,ZZB
      INTEGER I,ILK,IMM,INDM,J,K,L,MHDIMH
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,LH
         DO J = 1,MHDIM2
            PSI2G(J,I,IPOL) = CZERO
         END DO
      END DO
C
C     calculate psi2g in [1,mhdim1]
C
      DO IMM = 1,2
         DO K = 1,IGP
            IF ( IMM.EQ.1 ) INDM = IGVAL(K)
            IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
            IF ( IREL.EQ.1 ) THEN
               IF ( IGPRO(K).EQ.1 ) THEN
                  DO I = 1,MHDIM1
                     PSI2G(I,INDM,IPOL) = CMG(INDM)*WMHF(I,INDM)
     &                  + WPHF(I,INDM)
                  END DO
               ELSE
                  DO I = 1,MHDIM1
                     PSI2G(I,INDM,IPOL) = CMG(INDM)*WMHF(I,INDM)
                  END DO
               END IF
            ELSE IF ( IREL.EQ.2 ) THEN
               IF ( IGPRO(K).EQ.1 ) THEN
                  DO I = 1,MHDIM1
                     PSI2G(I,INDM,IPOL) = CMG(INDM)*WMHF(I,INDM)
     &                  + WPHF(I,INDM)
                  END DO
               ELSE
                  DO I = 1,MHDIM1
                     PSI2G(I,INDM,IPOL) = CMG(INDM)*WMHF(I,INDM)
                  END DO
               END IF
            END IF
         END DO
C
         MHDIMH = MHDIM/2
         IF ( IMM.EQ.1 ) THEN
            DO K = 1,GANZ
               IF ( IP.GT.1 ) WRITE (NOUT1,99001)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - 30.D0,
     &                               PSI2G(1,K,IPOL)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - 15.D0,
     &                               PSI2G(MHDIMH,K,IPOL)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - ZZ(3,IMM),
     &                               PSI2G(MHDIM1,K,IPOL)
            END DO
         ELSE IF ( IMM.EQ.2 ) THEN
            DO K = 1 + GANZ,QD
               IF ( IP.GT.1 ) WRITE (NOUT1,99001)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - 30.D0,
     &                               PSI2G(1,K,IPOL)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - 15.D0,
     &                               PSI2G(MHDIMH,K,IPOL)
               IF ( IP.GT.1 ) WRITE (NOUT1,99003) K, - ZZ(3,IMM),
     &                               PSI2G(MHDIM1,K,IPOL)
            END DO
         END IF
C
         IF ( IP.GT.1 ) THEN
            WRITE (NOUT1,99002)
C
            IF ( IREL.EQ.1 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) INDM = IGVAL(K)
                  IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
                  DELTA = U0GP(INDM) + U0GM(INDM)
     &                    - PSI2G(MHDIM1,INDM,IPOL)
                  DELTAS = CMG(INDM)*WMHG(MHDIM1,INDM)
     &                     + WPHG(MHDIM1,INDM) + CIMAG*KGZBLK(INDM)
     &                     *(U0GM(INDM)-U0GP(INDM))
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTA
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTAS
               END DO
            ELSE IF ( IREL.EQ.2 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) INDM = IGVAL(K)
                  IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
                  DELTA = U0GP(INDM) + U0GM(INDM)
     &                    - PSI2G(MHDIM1,INDM,IPOL)
                  DELTAS = CMG(INDM)*WMHG(MHDIM1,INDM)
     &                     + WPHG(MHDIM1,INDM) + CIMAG*KGZBLK(INDM)
     &                     *(U0GM(INDM)-U0GP(INDM))
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTA
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTAS
               END DO
            END IF
         END IF
C
C        calculate psi2g in [mhdim1+1,mhdim2]
C
         DO ILK = 1,IREL
            DO K = 1,IGP
               IF ( IMM.EQ.1 ) INDM = IGVAL(K)
               IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
C              ind = igval(k) + (ilk - 1)*ganz
               EZ = 2.D0*EZRYH(INDM)
               IF ( IGPRO(K).EQ.1 ) THEN
                  DO I = 1,6
                     ZZB = ZMESH(I)
                     F1(I) = CMG(INDM)*WMSTF(I,INDM)
     &                       + DCONJG(WMSTF(I,INDM))
                     F2(I) = CMG(INDM)*WMSTG(I,INDM)
     &                       + DCONJG(WMSTG(I,INDM))
                     F1S(I) = F2(I)
                     CALL CBPOT(EZ,ZZB,1,V,IMM)
                     F2S(I) = V*F1(I)
                  END DO
               ELSE
                  DO I = 1,6
                     ZZB = ZMESH(I)
                     F1(I) = CMG(INDM)*WMSTF(I,INDM)
                     F2(I) = CMG(INDM)*WMSTG(I,INDM)
                     F1S(I) = F2(I)
                     CALL CBPOT(EZ,ZZB,1,V,IMM)
                     F2S(I) = V*F1(I)
                  END DO
               END IF
               HB = ZZ(3,IMM)/MHDIM
C
C              calculating predictor for psi2g-solutions
C
               DO I = 1,MHDIM
                  F1(7) = F1(6) + (HB/1440)
     &                    *(4277*F1S(6)-7923*F1S(5)+9982*F1S(4)
     &                    -7298*F1S(3)+2877*F1S(2)-475*F1S(1))
                  F2(7) = F2(6) + (HB/1440)
     &                    *(4277*F2S(6)-7923*F2S(5)+9982*F2S(4)
     &                    -7298*F2S(3)+2877*F2S(2)-475*F2S(1))
                  F1S(7) = F2(7)
                  ZZB = I*HB
                  CALL CBPOT(EZ,ZZB,1,V,IMM)
                  F2S(7) = V*F1(7)
C
C                 corrector
C
                  F1(7) = F1(6) + (HB/60480)
     &                    *(19087*F1S(7)+65112*F1S(6)-46461*F1S(5)
     &                    +37504*F1S(4)-20211*F1S(3)+6312*F1S(2)
     &                    -863*F1S(1))
                  F2(7) = F2(6) + (HB/60480)
     &                    *(19087*F2S(7)+65112*F2S(6)-46461*F2S(5)
     &                    +37504*F2S(4)-20211*F2S(3)+6312*F2S(2)
     &                    -863*F2S(1))
                  F1S(7) = F2(7)
                  F2S(7) = V*F1(7)
                  PSI2G(I+MHDIM1,INDM,IPOL) = F1(7)
C
                  DO L = 1,6
                     F1(L) = F1(L+1)
                     F2(L) = F2(L+1)
                     F1S(L) = F1S(L+1)
                     F2S(L) = F2S(L+1)
                  END DO
C
               END DO
            END DO
         END DO
      END DO
C
      IF ( IP.GT.1 ) THEN
         DO K = 1,QD
            WRITE (NOUT1,99001)
            WRITE (NOUT1,99003) K,0.D0,PSI2G(MHDIM2,K,IPOL)
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (1x,'rvg',7x,'zr',15x,'psi2g')
99002 FORMAT (1x,'rvg',7x,'zr',6x,'deltapsi and deltapsis')
99003 FORMAT (1x,i3,2x,3E12.5)
      END
C*==cwh2.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CWH2(KGZH,GANZ,EZRYH,WMHF,WMHG,WPHF,WPHG,IGP,IGVAL,
     &                IGPRO,WMSTF,WMSTG,ZMESH,IZBG)
C     /****************************************************************/
C     # purpose      : calculates wave functions to energy ehigh       *
C                      in a rm-potential                               *
C                                                                      *
C     # subroutines called from this routine                           *
C       wh        cbpot                                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM1,MHSTEP,MHMAX,MHMAX1,CZERO
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,IGP,IZBG
      REAL*8 EZRYH(LH),ZMESH(6)
      INTEGER IGPRO(LH),IGVAL(LH)
      COMPLEX*16 KGZH(LH),WMHF(MHDIM1,LH),WMHG(MHDIM1,LH),WMSTF(6,LH),
     &           WMSTG(6,LH),WPHF(MHDIM1,LH),WPHG(MHDIM1,LH)
C
C Local variables
C
      REAL*8 EMACH,EZ,HR,TOLHF,TOLHG,ZR,ZZB
      COMPLEX*16 F1(7),F1S(7),F2(7),F2S(7),V,WFCT,WFCTS
      INTEGER G,I,IMM,INDM,J,JTEST,K
C
C*** End of declarations rewritten by SPAG
C
      EMACH = 1.D-6
C
      DO I = 1,MHDIM1
         DO J = 1,LH
            WMHF(I,J) = CZERO
            WMHG(I,J) = CZERO
            WPHF(I,J) = CZERO
            WPHG(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,6
         DO J = 1,LH
            WMSTF(I,J) = CZERO
            WMSTG(I,J) = CZERO
         END DO
      END DO
      DO IMM = 1,2
         HR = (30.D0-ZZ(3,IMM))/MHMAX
C
         DO G = 1,IZBG
            IF ( IMM.EQ.1 ) INDM = G
            IF ( IMM.EQ.2 ) INDM = G + GANZ
C
            IF ( IP.GT.1 ) WRITE (NOUT1,99001)
C
            EZ = 2.D0*EZRYH(INDM)
            DO I = 1,6
               ZZB = (-30.D0) + DBLE(I-1)*HR + ZZ(3,IMM)
               ZR = ZZB - ZZ(3,IMM)
               CALL WH(KGZH(INDM),ZR,1.D-4,WFCT,WFCTS,EZ)
               F1(I) = WFCT
               F2(I) = WFCTS
               F1S(I) = WFCTS
C
C             v = v(z) - vpr - ez  (rydberg)
C
               CALL CBPOT(EZ,ZZB,1,V,IMM)
               F2S(I) = V*WFCT
            END DO
            WMHF(1,INDM) = F1(1)
            WMHG(1,INDM) = F2(1)
            IF ( IP.GT.1 ) WRITE (NOUT1,99003) INDM, - 30.D0,
     &                            WMHF(1,INDM),WMHG(1,INDM)
C
C         calculating predictor for wgm-soltutions
C
            DO I = 7,MHMAX1
               F1(7) = F1(6) + (HR/1440)
     &                 *(4277*F1S(6)-7923*F1S(5)+9982*F1S(4)-7298*F1S(3)
     &                 +2877*F1S(2)-475*F1S(1))
               F2(7) = F2(6) + (HR/1440)
     &                 *(4277*F2S(6)-7923*F2S(5)+9982*F2S(4)-7298*F2S(3)
     &                 +2877*F2S(2)-475*F2S(1))
               F1S(7) = F2(7)
               ZZB = (-30.D0) + DBLE(I-1)*HR + ZZ(3,IMM)
               CALL CBPOT(EZ,ZZB,1,V,IMM)
               F2S(7) = V*F1(7)
C
C            corrector
C
               F1(7) = F1(6) + (HR/60480)
     &                 *(19087*F1S(7)+65112*F1S(6)-46461*F1S(5)
     &                 +37504*F1S(4)-20211*F1S(3)+6312*F1S(2)-863*F1S(1)
     &                 )
               F2(7) = F2(6) + (HR/60480)
     &                 *(19087*F2S(7)+65112*F2S(6)-46461*F2S(5)
     &                 +37504*F2S(4)-20211*F2S(3)+6312*F2S(2)-863*F2S(1)
     &                 )
               F1S(7) = F2(7)
               F2S(7) = V*F1(7)
               IF ( MOD(I-1,MHSTEP).EQ.0 ) THEN
                  J = 1 + INT((I-1)/MHSTEP)
                  WMHF(J,INDM) = F1(7)
                  WMHG(J,INDM) = F2(7)
               END IF
               IF ( I.GT.MHMAX1-6 ) THEN
                  WMSTF(I+6-MHMAX1,INDM) = F1(7)
                  WMSTG(I+6-MHMAX1,INDM) = F2(7)
                  ZMESH(I+6-MHMAX1) = ZZB
               END IF
               DO K = 1,6
                  F1(K) = F1(K+1)
                  F2(K) = F2(K+1)
                  F1S(K) = F1S(K+1)
                  F2S(K) = F2S(K+1)
               END DO
            END DO
C
            IF ( IP.GT.1 ) WRITE (NOUT1,99003) INDM,(-ZZ(3,IMM)),
     &                            WMHF(MHDIM1,INDM),WMHG(MHDIM1,INDM)
         END DO
C
         EZ = 2.D0*EZRYH(INDM)
         IF ( EZ.GT.0.D0 ) THEN
            IF ( IP.GT.1 ) WRITE (NOUT1,99002)
            DO I = 1,MHDIM1
               WPHF(I,INDM) = DCONJG(WMHF(I,INDM))
               WPHG(I,INDM) = DCONJG(WMHG(I,INDM))
               IF ( IP.GT.1 .AND. I.EQ.1 ) WRITE (NOUT1,99003) INDM,
     &              - 30.D0,WPHF(I,INDM),WPHG(I,INDM)
            END DO
            IF ( IP.GT.1 ) WRITE (NOUT1,99003) INDM,(-ZZ(3,IMM)),
     &                            WPHF(MHDIM1,INDM),WPHG(MHDIM1,INDM)
         END IF
C
         IGP = 0
         DO I = 1,IZBG
            TOLHF = CDABS(WMHF(MHDIM1,INDM))
            TOLHG = CDABS(WMHG(MHDIM1,INDM))
            IF ( TOLHF.GT.EMACH .AND. TOLHG.GT.EMACH ) THEN
               IGP = IGP + 1
               IGVAL(IGP) = I
            END IF
         END DO
C
         DO I = 1,IGP
            IGPRO(INDM) = 0
            EZ = 2.D0*EZRYH(INDM)
            IF ( EZ.GT.0.0D0 ) IGPRO(INDM) = 1
         END DO
C
         JTEST = 0
         DO I = 1,IGP
            JTEST = JTEST + IGPRO(INDM)
         END DO
C         stop
         IF ( JTEST.EQ.0 ) WRITE (NOUT1,99004)
      END DO
C
      RETURN
99001 FORMAT (1x,'rvg',7x,'zr',15x,'wfmg',20x,'wfmgs')
99002 FORMAT (1x,'rvg',7x,'zr',15x,'wfpg',20x,'wfpgs')
99003 FORMAT (8x,i3,2x,5E12.5)
99004 FORMAT (1x,'cdabs(propagating solution) .lt. emach')
      END
C*==wh.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE WH(KZ,Z,WTOL,W,WS,EZ)
C     /****************************************************************/
C     # purpose      : calculate asymptotic form of                    *
C                      whittaker functions                             *
C                                                                      *
C     # function called from this routine                              *
C       cdgam                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MHDIM,PI,CZERO,CONE,CIMAG
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EZ,WTOL,Z
      COMPLEX*16 KZ,W,WS
C
C Local variables
C
      COMPLEX*16 B,CDG,CDGE,CE,F1,F2,KZI,SUM1,SUM2
      COMPLEX*16 CDGAM
      INTEGER I,ILAST
      REAL*8 R1,R2
      EXTERNAL CDGAM
C
C*** End of declarations rewritten by SPAG
C
      W = CZERO
      WS = CZERO
      ILAST = 0
      KZI = CONE/KZ
C
      IF ( EZ.GE.(-0.05D0) .AND. EZ.LE.0.05D0 ) THEN
         IF ( EZ.LE.0.D0 ) THEN
            EZ = -0.05D0
            KZ = CDSQRT(DCMPLX(2.D0*EZ,0.D0))
            KZI = CONE/KZ
         ELSE
            EZ = 0.05D0
            KZ = CDSQRT(DCMPLX(2.D0*EZ,0.D0))
            KZI = CONE/KZ
         END IF
      END IF
C
      CE = -CIMAG*(KZ*Z-0.25D0*KZI*CDLOG(2.D0*KZ*Z))
      F1 = 2.D0*CDEXP((-3.D0*PI*KZI/8.D0))
      CDG = CONE + CIMAG*0.25D0*KZI
C
      CDGE = CDGAM(CDG)
C
      CDG = CONE/CDGE
      F1 = F1*CDEXP(CE)*CDG
      F2 = 0.25D0*CIMAG*KZI/Z - CIMAG*KZ
      B = CONE
      SUM1 = B
      SUM2 = CZERO
      DO I = 1,MHDIM
         R1 = DBLE(I)
         R2 = DBLE(I-1)
         B = B*(R1-0.25D0*CIMAG*KZI)
         B = B*(R2-0.25D0*CIMAG*KZI)
         B = B/(-2.D0*CIMAG*KZ*R1*Z)
         IF ( CDABS(B).LT.WTOL ) EXIT
         ILAST = I
         SUM1 = SUM1 + B
         SUM2 = SUM2 - R1*B/Z
      END DO
C
      IF ( ILAST.GE.MHDIM ) THEN
         WRITE (*,99001)
         WRITE (NOUT1,99001)
         WRITE (NOUT1,99002) MHDIM
         STOP
      ELSE
         W = F1*SUM1
         WS = F1*(SUM2+F2*SUM1)
      END IF
C
      RETURN
99001 FORMAT (2x,'asymptotic series in wh not converged !')
99002 FORMAT (2x,'ilast = mhdim =',i4)
      END
C*==derivat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DERIVAT(POTDER,H,FLAG,EZRY,IMM)
C     /****************************************************************/
C     # purpose      : calculate derivative of cbpot                   *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       cbpot                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MHDIM,MHDIM1,MHDIM2,MHSTEP,MHMAX1
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EZRY,H
      INTEGER FLAG,IMM
      COMPLEX*16 POTDER(MHDIM2)
C
C Local variables
C
      REAL*8 EZ,F1I(3),F1R(3),F2I(3),F2R(3),H2,IMV(2),REV(2),ZWEIH(3),
     &       ZZB
      INTEGER I,IDER,J,L
      COMPLEX*16 V
C
C*** End of declarations rewritten by SPAG
C
      IF ( FLAG.EQ.1 ) THEN
         IDER = MHMAX1
C        ider = mhdim1
      ELSE IF ( FLAG.EQ.2 ) THEN
         IDER = MHDIM
      END IF
C
      H2 = H/2.D0
      ZWEIH(1) = 2.D0
      ZWEIH(2) = 4.D0
      ZWEIH(3) = 8.D0
      EZ = EZRY
      DO I = 1,IDER
         DO J = 1,3
            IF ( FLAG.EQ.1 ) THEN
               ZZB = ZZ(3,IMM) - 30.D0 + DBLE(I-1)*H + H2*ZWEIH(J)
            ELSE IF ( FLAG.EQ.2 ) THEN
               ZZB = DBLE(I)*H + H2*ZWEIH(J)
            END IF
            CALL CBPOT(EZ,ZZB,1,V,IMM)
            REV(1) = DBLE(V)
            IMV(1) = DIMAG(V)
C
            IF ( FLAG.EQ.1 ) THEN
               ZZB = ZZ(3,IMM) - 30.D0 + DBLE(I-1)*H - H2*ZWEIH(J)
            ELSE IF ( FLAG.EQ.2 ) THEN
               ZZB = DBLE(I)*H - H2*ZWEIH(J)
            END IF
            CALL CBPOT(EZ,ZZB,1,V,IMM)
            REV(2) = DBLE(V)
            IMV(2) = DIMAG(V)
C
            F1R(J) = (REV(1)-REV(2))/(H*ZWEIH(J))
            F1I(J) = (IMV(1)-IMV(2))/(H*ZWEIH(J))
         END DO
         F2R(1) = F1R(1) + (F1R(1)-F1R(2))/3.D0
         F2I(1) = F1I(1) + (F1I(1)-F1I(2))/3.D0
         F2R(2) = F1R(2) + (F1R(2)-F1R(3))/3.D0
         F2I(2) = F1I(2) + (F1I(2)-F1I(3))/3.D0
         F2R(3) = F2R(1) + (F2R(1)-F2R(2))/15.D0
         F2I(3) = F2I(1) + (F2I(1)-F2I(2))/15.D0
         IF ( FLAG.EQ.1 ) THEN
            IF ( MOD(I-1,MHSTEP).EQ.0 ) THEN
               L = 1 + INT(DBLE(I-1)/MHSTEP)
               POTDER(L) = DCMPLX(F2R(3),F2I(3))
C                 write (8,*) 'FLAG=1',l,potder(l)
            END IF
C              potder(i) = dcmplx(f2r(3),f2i(3))
         ELSE IF ( FLAG.EQ.2 ) THEN
            L = I + MHDIM1
            POTDER(L) = DCMPLX(F2R(3),F2I(3))
C                 write (9,*) 'FLAG=2',l,potder(l)
         END IF
      END DO
C
      END
C*==cwl2.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CWL2(KGZL,GANZ,EZRYL,WMHF,WMHG,WPHF,WPHG,IGP,IGVAL,
     &                IGPRO,WMSTF,WMSTG,ZMESH,IZBG)
C     /****************************************************************/
C     # purpose      : calculates wave functions to energy elow        *
C                      in a rm-potential                               *
C                                                                      *
C     # subroutines called from this routine                           *
C       wh        cbpot                                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM1,MHSTEP,MHMAX,MHMAX1,CZERO
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,IGP,IZBG
      REAL*8 EZRYL(LH),ZMESH(6)
      INTEGER IGPRO(LH),IGVAL(LH)
      COMPLEX*16 KGZL(LH),WMHF(MHDIM1,LH),WMHG(MHDIM1,LH),WMSTF(6,LH),
     &           WMSTG(6,LH),WPHF(MHDIM1,LH),WPHG(MHDIM1,LH)
C
C Local variables
C
      REAL*8 EMACH,EZ,HR,TOLHF,TOLHG,ZR,ZZB
      COMPLEX*16 F1(7),F1S(7),F2(7),F2S(7),V,WFCT,WFCTS
      INTEGER G,I,IMM,INDM,J,K
C
C*** End of declarations rewritten by SPAG
C
      EMACH = 1.0D-6
C
      DO I = 1,MHDIM1
         DO J = 1,LH
            WMHF(I,J) = CZERO
            WMHG(I,J) = CZERO
            WPHF(I,J) = CZERO
            WPHG(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,6
         DO J = 1,LH
            WMSTF(I,J) = CZERO
            WMSTG(I,J) = CZERO
         END DO
      END DO
C
      DO IMM = 1,2
         HR = (30.D0-ZZ(3,IMM))/MHMAX
         DO G = 1,IZBG
            IF ( IMM.EQ.1 ) INDM = G
            IF ( IMM.EQ.2 ) INDM = G + GANZ
            IF ( IP.GT.1 ) WRITE (NOUT1,99001)
            EZ = 2.0*EZRYL(INDM)
            DO I = 1,6
               ZZB = (-30.D0) + DBLE(I-1)*HR + ZZ(3,IMM)
               ZR = ZZB - ZZ(3,IMM)
               CALL WH(KGZL(INDM),ZR,1.0D-4,WFCT,WFCTS,EZ)
               F1(I) = WFCT
               F2(I) = WFCTS
               F1S(I) = WFCTS
               CALL CBPOT(EZ,ZZB,2,V,IMM)
               F2S(I) = V*WFCT
            END DO
            WMHF(1,INDM) = F1(1)
            WMHG(1,INDM) = F2(1)
            IF ( IP.GT.1 ) WRITE (NOUT1,99002) INDM, - 30.D0,
     &                            WMHF(1,INDM),WMHG(1,INDM)
C
C         calculating predictor for wgm-soltutions
C
            DO I = 7,MHMAX1
               F1(7) = F1(6) + (HR/1440)
     &                 *(4277*F1S(6)-7923*F1S(5)+9982*F1S(4)-7298*F1S(3)
     &                 +2877*F1S(2)-475*F1S(1))
               F2(7) = F2(6) + (HR/1440)
     &                 *(4277*F2S(6)-7923*F2S(5)+9982*F2S(4)-7298*F2S(3)
     &                 +2877*F2S(2)-475*F2S(1))
               F1S(7) = F2(7)
               ZZB = (-30.D0) + DBLE(I-1)*HR + ZZ(3,IMM)
               CALL CBPOT(EZ,ZZB,2,V,IMM)
               F2S(7) = V*F1(7)
C
C              corrector
C
               F1(7) = F1(6) + (HR/60480)
     &                 *(19087*F1S(7)+65112*F1S(6)-46461*F1S(5)
     &                 +37504*F1S(4)-20211*F1S(3)+6312*F1S(2)-863*F1S(1)
     &                 )
               F2(7) = F2(6) + (HR/60480)
     &                 *(19087*F2S(7)+65112*F2S(6)-46461*F2S(5)
     &                 +37504*F2S(4)-20211*F2S(3)+6312*F2S(2)-863*F2S(1)
     &                 )
               F1S(7) = F2(7)
               F2S(7) = V*F1(7)
               IF ( MOD(I-1,MHSTEP).EQ.0 ) THEN
                  J = 1 + INT((I-1)/MHSTEP)
                  WMHF(J,INDM) = F1(7)
                  WMHG(J,INDM) = F2(7)
               END IF
               IF ( I.GT.MHMAX1-6 ) THEN
                  WMSTF(I+6-MHMAX1,INDM) = F1(7)
                  WMSTG(I+6-MHMAX1,INDM) = F2(7)
                  ZMESH(I+6-MHMAX1) = ZZB
               END IF
               DO K = 1,6
                  F1(K) = F1(K+1)
                  F2(K) = F2(K+1)
                  F1S(K) = F1S(K+1)
                  F2S(K) = F2S(K+1)
               END DO
            END DO
            IF ( IP.GT.1 ) WRITE (NOUT1,99002) INDM,(-ZZ(3,IMM)),
     &                            WMHF(MHDIM1,INDM),WMHG(MHDIM1,INDM)
         END DO
C
         IGP = 0
         DO I = 1,IZBG
            IF ( IMM.EQ.1 ) INDM = I
            IF ( IMM.EQ.2 ) INDM = I + GANZ
            TOLHF = CDABS(WMHF(MHDIM1,INDM))
            TOLHG = CDABS(WMHG(MHDIM1,INDM))
            IF ( TOLHF.GT.EMACH .AND. TOLHG.GT.EMACH ) THEN
               IGP = IGP + 1
               IGVAL(IGP) = I
            END IF
         END DO
         DO I = 1,IGP
            IGPRO(I) = 0
         END DO
      END DO
C
      RETURN
C
99001 FORMAT (1x,'rvg',7x,'zr',15x,'wfmg',20x,'wfmgs')
99002 FORMAT (1x,i3,2x,5E12.5)
      END
C*==cphi1g.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CPHI1G(IREL,GANZ,IGP,QD,IGVAL,ZMESH,EZRYL,A1GZ,A1GZS,
     &                  CMG,CMGS,WMHF,WMHG,D0GM,D0GP,KGZBLK,PHI1G,
     &                  PHI1GS)
C     /****************************************************************/
C     # purpose     : calculates the final state wavefield             *
C                     at the first bulk-layer                          *
C                                                                      *
C     # subroutines called from this routine                           *
C       cbpot                                                          *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM,MHDIM1,MHDIM2,CZERO,CIMAG
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER GANZ,IGP,IREL,QD
      COMPLEX*16 A1GZ(MHDIM1,LH),A1GZS(MHDIM1,LH),CMG(LH),CMGS(LH),
     &           D0GM(LH),D0GP(LH),KGZBLK(LH),PHI1G(MHDIM2,LH),
     &           PHI1GS(MHDIM2,LH),WMHF(MHDIM1,LH),WMHG(MHDIM1,LH)
      REAL*8 EZRYL(LH),ZMESH(6)
      INTEGER IGVAL(LH)
C
C Local variables
C
      COMPLEX*16 DELTA,DELTAS,F1(7),F1S(7),F2(7),F2S(7),V
      REAL*8 EZ,HB,ZZB
      INTEGER I,ILK,IMM,INDM,IX,J,K,L,MHDIMH
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,MHDIM2
         DO J = 1,LH
            PHI1G(I,J) = CZERO
            PHI1GS(I,J) = CZERO
         END DO
      END DO
C
C     calclate phi1g in [1,mhdim1]
C
      DO IMM = 1,2
         IF ( IGP.GT.0 ) THEN
            DO K = 1,IGP
               IF ( IMM.EQ.1 ) INDM = IGVAL(K)
               IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
               IF ( IREL.EQ.1 ) THEN
                  DO I = 1,MHDIM1
                     PHI1G(I,INDM) = CMG(INDM)
     &                               *(WMHF(I,INDM)+A1GZ(I,INDM)
     &                               +CMGS(INDM))
                     PHI1GS(I,INDM) = CMG(INDM)
     &                                *(WMHG(I,INDM)+A1GZS(I,INDM))
                  END DO
               ELSE IF ( IREL.EQ.2 ) THEN
                  DO I = 1,MHDIM1
                     PHI1G(I,INDM) = CMG(INDM)
     &                               *(WMHF(I,INDM)+A1GZ(I,INDM)
     &                               +CMGS(INDM))
                     PHI1GS(I,INDM) = CMG(INDM)
     &                                *(WMHG(I,INDM)+A1GZS(I,INDM))
                  END DO
               END IF
            END DO
         END IF
C
         MHDIMH = MHDIM/2
         IF ( IP.GT.1 ) THEN
            IF ( IMM.EQ.1 ) THEN
               DO K = 1,GANZ
                  WRITE (NOUT1,99001)
                  WRITE (NOUT1,99003) K, - 30.0,PHI1G(1,K)
                  WRITE (NOUT1,99003) K, - 15.0,PHI1G(MHDIMH,K)
                  WRITE (NOUT1,99003) K, - ZZ(3,IMM),PHI1G(MHDIM1,K)
               END DO
            ELSE IF ( IMM.EQ.2 ) THEN
               DO K = 1 + GANZ,QD
                  WRITE (NOUT1,99001)
                  WRITE (NOUT1,99003) K, - 30.0,PHI1G(1,K)
                  WRITE (NOUT1,99003) K, - 15.0,PHI1G(MHDIMH,K)
                  WRITE (NOUT1,99003) K, - ZZ(3,IMM),PHI1G(MHDIM1,K)
               END DO
            END IF
         END IF
C
         IF ( IP.GT.1 ) THEN
            WRITE (NOUT1,99002)
            IF ( IREL.EQ.1 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) INDM = IGVAL(K)
                  IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
                  DELTA = D0GP(INDM) + D0GM(INDM) - PHI1G(MHDIM1,INDM)
                  DELTAS = PHI1GS(MHDIM1,INDM) + CIMAG*KGZBLK(INDM)
     &                     *(D0GM(INDM)-D0GP(INDM))
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTA
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTAS
               END DO
            ELSE IF ( IREL.EQ.2 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) INDM = IGVAL(K)
                  IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
                  DELTA = D0GP(INDM) + D0GM(INDM) - PHI1G(MHDIM1,INDM)
                  DELTAS = PHI1GS(MHDIM1,INDM) + CIMAG*KGZBLK(INDM)
     &                     *(D0GM(INDM)-D0GP(INDM))
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTA
                  IF ( IP.GT.1 ) WRITE (NOUT1,99003) K,(-ZZ(3,IMM)),
     &                                  DELTAS
               END DO
            END IF
         END IF
C
C        calculate phi1g in [mhdim1+1,mhdim2]
C
         IF ( IGP.GT.0 ) THEN
            DO ILK = 1,IREL
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) INDM = IGVAL(K)
                  IF ( IMM.EQ.2 ) INDM = IGVAL(K) + GANZ
C                  ind = igval(k) + (ilk - 1)*ganz
                  EZ = 2.D0*EZRYL(INDM)
                  DO I = 1,6
                     IX = MHDIM1 - 6 + I
                     ZZB = ZMESH(I)
                     F1(I) = PHI1G(IX,INDM)
                     F2(I) = PHI1GS(IX,INDM)
                     F1S(I) = F2(I)
                     CALL CBPOT(EZ,ZZB,1,V,IMM)
                     F2S(I) = V*F1(I)
                  END DO
                  HB = ZZ(3,IMM)/MHDIM
C
C                  calculating predictor for phi1g-soltutions
C
                  DO I = 1,MHDIM
                     F1(7) = F1(6) + (HB/1440)
     &                       *(4277*F1S(6)-7923*F1S(5)+9982*F1S(4)
     &                       -7298*F1S(3)+2877*F1S(2)-475*F1S(1))
                     F2(7) = F2(6) + (HB/1440)
     &                       *(4277*F2S(6)-7923*F2S(5)+9982*F2S(4)
     &                       -7298*F2S(3)+2877*F2S(2)-475*F2S(1))
                     F1S(7) = F2(7)
                     ZZB = DBLE(I)*HB
                     CALL CBPOT(EZ,ZZB,1,V,IMM)
                     F2S(7) = V*F1(7)
C
C                     corrector
C
                     F1(7) = F1(6) + (HB/60480)
     &                       *(19087*F1S(7)+65112*F1S(6)-46461*F1S(5)
     &                       +37504*F1S(4)-20211*F1S(3)+6312*F1S(2)
     &                       -863*F1S(1))
                     F2(7) = F2(6) + (HB/60480)
     &                       *(19087*F2S(7)+65112*F2S(6)-46461*F2S(5)
     &                       +37504*F2S(4)-20211*F2S(3)+6312*F2S(2)
     &                       -863*F2S(1))
                     F1S(7) = F2(7)
                     F2S(7) = V*F1(7)
                     PHI1G(I+MHDIM1,INDM) = F1(7)
                     PHI1GS(I+MHDIM1,INDM) = F2(7)
                     DO L = 1,6
                        F1(L) = F1(L+1)
                        F2(L) = F2(L+1)
                        F1S(L) = F1S(L+1)
                        F2S(L) = F2S(L+1)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
C
      MHDIMH = MHDIM/2
      IF ( IP.GT.1 ) THEN
         DO K = 1,GANZ
            WRITE (NOUT1,99001)
            WRITE (NOUT1,99003) K,0.0,PHI1G(MHDIM2,K)
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (1x,'rvg',7x,'zr',15x,'phi1g')
99002 FORMAT (1x,'rvg',7x,'zr',6x,'deltapsi and deltapsis')
99003 FORMAT (1x,i3,2x,3E12.5)
      END
C*==ca1gp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE CA1GP(K1GZ,GANZ,PSI2G,A1GP,QD,IREL,AAZ,EZRY,A1GZ,OMEGA,
     &                 RS,A1GPS,POTDER,A1GZS,VPR,IGP,IGVAL,IPOL)
C     /****************************************************************/
C     # purpose      :                                                 *
C                                                                      *
C     # subroutines and functions called from this routine             *
C       derivat   whsimp                                               *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM1,MHDIM2,MHSTEP,MHMAX,CZERO,CONE
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 AAZ,EZRY,OMEGA,VPR
      INTEGER GANZ,IGP,IPOL,IREL,QD
      COMPLEX*16 A1GP(LH),A1GPS(LH),A1GZ(MHDIM1,LH),A1GZS(MHDIM1,LH),
     &           K1GZ(LH),POTDER(MHDIM2),PSI2G(MHDIM2,LH,4),RS(LH)
      INTEGER IGVAL(LH)
C
C Local variables
C
      COMPLEX*16 CONST,FACTOR,NORM(LH)
      REAL*8 EMACH,HR,HRM
      INTEGER I,IH,IMM,J,K
C
C*** End of declarations rewritten by SPAG
C
      EMACH = 1.D-6
C
      DO I = 1,LH
         A1GP(I) = CZERO
         A1GPS(I) = CZERO
      END DO
C
      DO I = 1,LH
         DO J = 1,MHDIM1
            A1GZ(J,I) = CZERO
            A1GZS(J,I) = CZERO
         END DO
      END DO
C
      CONST = DCMPLX(1.D0/(2.D0*OMEGA*137.036),0.D0)
C      IF ( AAZ.NE.0.D0 ) THEN
      IF ( ABS(AAZ).GT.1.0D-16 ) THEN
         DO IMM = 1,2
            HR = (30.D0-ZZ(3,IMM))/MHMAX
            HRM = HR*MHSTEP
C
C           calculating the derivative of the potential by a richardson
C           interpolation procedure up to o(h**6)
C
            CALL DERIVAT(POTDER,HR,1,EZRY,IMM)
C
            IF ( IMM.EQ.1 ) THEN
               DO K = 1,GANZ
                  FACTOR = AAZ*CONST/K1GZ(K)
                  NORM(K) = FACTOR*VPR*(CONE+RS(K))
     &                      *DCONJG(PSI2G(MHDIM1,K,IPOL))
                  DO I = 1,MHDIM1
                     A1GZS(I,K) = FACTOR*POTDER(I)
     &                            *DCONJG(PSI2G(I,K,IPOL))*(CONE+RS(K))
                  END DO
                  CALL WHSIMP(A1GZS,HRM,A1GZ,K)
               END DO
            ELSE
               DO K = 1 + GANZ,QD
                  FACTOR = AAZ*CONST/K1GZ(K)
                  NORM(K) = FACTOR*VPR*(CONE+RS(K))
     &                      *DCONJG(PSI2G(MHDIM1,K,IPOL))
                  DO I = 1,MHDIM1
                     A1GZS(I,K) = FACTOR*POTDER(I)
     &                            *DCONJG(PSI2G(I,K,IPOL))*(CONE+RS(K))
                  END DO
                  CALL WHSIMP(A1GZS,HRM,A1GZ,K)
               END DO
            END IF
C
            IF ( IREL.EQ.1 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) IH = IGVAL(K)
                  IF ( IMM.EQ.2 ) IH = IGVAL(K) + GANZ
                  DO I = 1,MHDIM1
                     A1GZ(I,IH) = NORM(IH)*A1GZ(I,IH)/A1GZ(MHDIM1,IH)
                     A1GZS(I,IH) = NORM(IH)*A1GZS(I,IH)/A1GZ(MHDIM1,IH)
                  END DO
               END DO
            ELSE IF ( IREL.EQ.2 ) THEN
               DO K = 1,IGP
                  IF ( IMM.EQ.1 ) IH = IGVAL(K)
                  IF ( IMM.EQ.2 ) IH = IGVAL(K) + GANZ
                  DO I = 1,MHDIM1
                     IF ( CDABS(A1GZ(MHDIM1,IH)).GT.EMACH ) THEN
                        A1GZ(I,IH) = NORM(IH)*A1GZ(I,IH)/A1GZ(MHDIM1,IH)
                        A1GZS(I,IH) = NORM(IH)*A1GZS(I,IH)
     &                                /A1GZ(MHDIM1,IH)
                     ELSE
                        A1GZ(I,IH) = CZERO
                        A1GZS(I,IH) = CZERO
                     END IF
                  END DO
               END DO
            END IF
         END DO
C
         DO K = 1,QD
            A1GP(K) = A1GZ(MHDIM1,K)
            A1GPS(K) = A1GZS(MHDIM1,K)
         END DO
      END IF
C
      END
C*==whsimp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE WHSIMP(A1GZS,H,A1GZ,K)
C     /****************************************************************/
C     # purpose      :                                                 *
C                                                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,MHDIM1,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 H
      INTEGER K
      COMPLEX*16 A1GZ(MHDIM1,LH),A1GZS(MHDIM1,LH)
C
C Local variables
C
      REAL*8 H3
      INTEGER I,J
      COMPLEX*16 SUMVAR
C
C*** End of declarations rewritten by SPAG
C
      H3 = H/3.D0
      A1GZ(1,K) = CZERO
      J = 1
      SUMVAR = CZERO
      DO I = 1,MHDIM1 - 2,2
         SUMVAR = SUMVAR + H3*(A1GZS(I,K)+4.D0*A1GZS(I+1,K)+A1GZS(I+2,K)
     &            )
         A1GZ(2*J+1,K) = SUMVAR
         J = J + 1
      END DO
C
      A1GZ(2,K) = (H/2.D0)*(A1GZS(1,K)+A1GZS(2,K))
      J = 1
      SUMVAR = CZERO
      DO I = 1,MHDIM1 - 3,2
         SUMVAR = SUMVAR + H3*(A1GZS(I+1,K)+4.D0*A1GZS(I+2,K)
     &            +A1GZS(I+3,K))
         A1GZ(2*J+2,K) = SUMVAR
         J = J + 1
      END DO
C
      END
C*==rhosva.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE RHOSVA(OMEGA,AAZ,EZRY,PSI2G,PHI1G,POTDER,ROSUR,ROSULA,
     &                  GANZ,IPOL,SPOL,ROSURSD,ROSULASD)
C     /****************************************************************/
C     # purpose      : calculate surface contribution                  *
C                                                                      *
C     # calls:  derivat                                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LL,LH,MHDIM,MHDIM1,MHDIM2,PI,CZERO,CIMAG
      USE MOD_SPEC_COMPOT,ONLY:ZZ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 AAZ,EZRY,OMEGA,ROSUR
      INTEGER GANZ,IPOL,SPOL
      COMPLEX*16 PHI1G(MHDIM2,LH),POTDER(MHDIM2),PSI2G(MHDIM2,LH,4),
     &           ROSULASD(LL),ROSURSD(4)
      REAL*8 ROSULA(LL)
C
C Local variables
C
      REAL*8 HB
      INTEGER IMM,INDM,J
      COMPLEX*16 SUMA1,SUMO
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,LL
         ROSULA(J) = 0.0D0
         ROSULASD(J) = CZERO
      END DO
C
      SUMA1 = CZERO
      DO IMM = 1,2
         HB = ZZ(3,IMM)/MHDIM
         CALL DERIVAT(POTDER,HB,2,EZRY,IMM)
C
         DO J = 1,GANZ
            IF ( IMM.EQ.1 ) INDM = J
            IF ( IMM.EQ.2 ) INDM = J + GANZ
            SUMA1 = SUMA1 + PSI2G(MHDIM1,INDM,IPOL)*PHI1G(MHDIM1,INDM)
     &              *POTDER(MHDIM1)
         END DO
      END DO
C
      SUMO = AAZ*SUMA1*DCMPLX(2.D0/(137.036*OMEGA),0.D0)
      IF ( SPOL.EQ.2 ) THEN
         ROSUR = -DIMAG(SUMO*CIMAG)/PI
         ROSULA(1) = ROSUR
      ELSE IF ( SPOL.EQ.4 ) THEN
         ROSURSD(IPOL) = -CIMAG*SUMO/PI
         ROSULASD(1) = ROSURSD(IPOL)
      END IF
C
      END
C*==bmat.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE BMAT(BPP,BPM,MPM,MPP,MMP,MMM,TPM,TPP,TMP,TMM,PS,QD,
     &                ISTREU,IBLOCH)
C     /****************************************************************/
C     purpose      : copy reflection coefficiens                       *
C                    or switch them off                                *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LH,CZERO,CONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IBLOCH,ISTREU,QD
      COMPLEX*16 BPM(QD,QD),BPP(QD,QD),MMM(QD,QD),MMP(QD,QD),MPM(QD,QD),
     &           MPP(QD,QD),PS(LH,2),TMM(QD),TMP(QD),TPM(QD),TPP(QD)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,QD
         DO J = 1,QD
            BPP(I,J) = CZERO
            BPM(I,J) = CZERO
            MPM(I,J) = CZERO
            MPP(I,J) = CZERO
            MMP(I,J) = CZERO
            MMM(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,QD
         BPP(I,I) = PS(I,1)
         BPM(I,I) = PS(I,2)
C        BPP(I,I) = PS(I,2)
C        BPM(I,I) = PS(I,1)
      END DO
      IF ( ISTREU.EQ.2 ) THEN
         DO I = 1,QD
            IF ( IBLOCH.EQ.1 ) THEN
               MPM(I,I) = TMP(I)
               MPP(I,I) = TMM(I)
               MMP(I,I) = TPM(I)
               MMM(I,I) = TPP(I)
            ELSE
               MPM(I,I) = TPM(I)
               MPP(I,I) = TPP(I)
               MMP(I,I) = TMP(I)
               MMM(I,I) = TMM(I)
            END IF
C
         END DO
      ELSE
         DO I = 1,QD
            MPM(I,I) = CZERO
            MPP(I,I) = CONE
            MMP(I,I) = CZERO
            MMM(I,I) = CZERO
         END DO
      END IF
C
      END
C
