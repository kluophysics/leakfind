C*==dmft_vertex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_VERTEX(UU,UJ,UC,UD,UM,US,UT,NLM,L,IPRINT)
C
      IMPLICIT NONE
C*--DMFT_VERTEX5
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,L,NLM
      REAL*8 UJ,UU
      COMPLEX*16 UC(NLM,NLM,NLM,NLM),UD(NLM,NLM,NLM,NLM),
     &           UM(NLM,NLM,NLM,NLM),US(NLM,NLM,NLM,NLM),
     &           UT(NLM,NLM,NLM,NLM)
C
C Local variables
C
      INTEGER IAUM,K,M1,M2,M3,M4
      REAL*8 RCL(4,1),UCL(7,7,7,7),UR(7,7,7,7)
      COMPLEX*16 U(NLM,NLM,NLM,NLM)
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C----------- For a moment - d-orbital and 1 atom per cell !!!!!!!!
      IAUM = 1
      K = 1
      RCL(1,K) = UU
      IF ( L.EQ.2 ) THEN
         RCL(2,K) = UJ*14D0/(1.D0+0.63D0)
         RCL(3,K) = 0.63D0*RCL(2,K)
      ELSE IF ( L.EQ.3 ) THEN
         RCL(2,K) = 6435.D0*UJ/(286.D0+195.D0*451.D0/675.D0+250.D0*
     &              1001.D0/2025.D0)
         RCL(3,K) = 451.D0*RCL(2,K)/675.D0
         RCL(4,K) = 1001.D0*RCL(2,K)/2025.D0
      END IF
C
      IF ( IPRINT.GT.0 ) WRITE (6,'(A,1x,10E16.6)') 'F0, F2, F4 =',
     &                          RCL(1,K),RCL(2,K),RCL(3,K)
C
      CALL DMFT_U4IND(UR,UCL,RCL,L,IAUM)
C--------- Only for Fe with s,d (6-orbitals!!!!!!!!!!)
      DO M1 = 1,NLM
         DO M2 = 1,NLM
            DO M3 = 1,NLM
               DO M4 = 1,NLM
                  U(M1,M2,M3,M4) = DCMPLX(UCL(M1,M2,M3,M4),0.0D0)
               END DO
            END DO
         END DO
      END DO
C-------------- Vertex for different channel
      DO M1 = 1,NLM
         DO M2 = 1,NLM
            DO M3 = 1,NLM
               DO M4 = 1,NLM
C------------- Uc - just COMPLEX U
                  UC(M1,M2,M3,M4) = U(M1,M2,M3,M4)
                  UD(M1,M2,M3,M4) = 2.D0*U(M1,M3,M2,M4) - U(M1,M3,M4,M2)
                  UM(M1,M2,M3,M4) = -U(M1,M3,M4,M2)
                  US(M1,M2,M3,M4) = 0.5D0*(U(M1,M2,M3,M4)+U(M1,M2,M4,M3)
     &                              )
                  UT(M1,M2,M3,M4) = 0.5D0*(U(M1,M2,M3,M4)-U(M1,M2,M4,M3)
     &                              )
               END DO
            END DO
         END DO
      END DO
C
      END
C*==dmft_u4ind.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
      SUBROUTINE DMFT_U4IND(U,UC,RCL,L,IAUM)
C..................................................................u4ind
C----> calculation of <m1m2|1/r12|m3m4> coulomb integrals
C
C
      IMPLICIT NONE
C*--DMFT_U4IND94
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IAUM,L
      DOUBLE PRECISION RCL(4,*),U(7,7,7,7),UC(7,7,7,7)
C
C Local variables
C
      DOUBLE PRECISION AM,CGK0,CGK1,CGK2,XK,XL,XM,XM1,XM2,XM3,XM4,
     &                 YOI(7,7),YOR(7,7)
      DOUBLE COMPLEX AM1,AM2,AM3,AM4
      DOUBLE PRECISION CGK
      INTEGER K,K2P1,MMAX,MS1,MS2,MS3,MS4,MS5,MS6,MS7,MS8
C
C*** End of declarations rewritten by SPAG
C
C
C Passed variables
C
C Local variables
C
C
      MMAX = 2*L + 1
      XL = DFLOAT(L)
C
      UC(1:7,1:7,1:7,1:7) = 0.0D0
C
      DO K = 0,2*L,2
         K2P1 = K/2 + 1
         XK = DFLOAT(K)
         CGK0 = CGK(XL,0.D0,XK,0.D0,XL,0.D0)
         DO MS1 = 1,MMAX
            XM1 = DFLOAT(MS1-L-1)
            DO MS2 = 1,MMAX
               XM2 = DFLOAT(MS2-L-1)
               DO MS3 = 1,MMAX
                  XM3 = DFLOAT(MS3-L-1)
                  XM = XM1 - XM3
                  DO MS4 = 1,MMAX
                     IF ( (MS1+MS2-MS3-MS4).EQ.0 ) THEN
                        XM4 = DFLOAT(MS4-L-1)
                        CGK1 = CGK(XL,XM3,XK,XM,XL,XM1)
                        CGK2 = CGK(XL,XM2,XK,XM,XL,XM4)
                        UC(MS1,MS2,MS3,MS4) = UC(MS1,MS2,MS3,MS4)
     &                     + RCL(K2P1,IAUM)*CGK0*CGK0*CGK1*CGK2
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL DMFT_CTORMT(YOR,YOI,L)
      U(1:7,1:7,1:7,1:7) = 0.0D0
C
      DO MS1 = 1,MMAX
         DO MS2 = 1,MMAX
            DO MS3 = 1,MMAX
               DO MS4 = 1,MMAX
                  DO MS5 = 1,MMAX
                     AM1 = DCMPLX(YOR(MS1,MS5),-YOI(MS1,MS5))
                     DO MS6 = 1,MMAX
                        AM2 = DCMPLX(YOR(MS2,MS6),-YOI(MS2,MS6))
                        DO MS7 = 1,MMAX
                           AM3 = DCMPLX(YOR(MS3,MS7),YOI(MS3,MS7))
                           DO MS8 = 1,MMAX
                              AM4 = DCMPLX(YOR(MS4,MS8),YOI(MS4,MS8))
                              AM = DREAL(AM1*AM2*AM3*AM4)
                              U(MS1,MS2,MS3,MS4) = U(MS1,MS2,MS3,MS4)
     &                           + AM*UC(MS5,MS6,MS7,MS8)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      END
C*==cgk.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C....................................................................cgk
      DOUBLE PRECISION FUNCTION CGK(A,AL,B,BE,C,GA)
      IMPLICIT NONE
C*--CGK192
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION A,AL,B,BE,C,GA
C
C Local variables
C
      DOUBLE PRECISION FA(0:20)
      INTEGER I1,I10,I11,I12,I13,I2,I3,I4,I5,I6,I7,I8,I9,Z,ZMAX,ZMIN
C
C*** End of declarations rewritten by SPAG
C
C
C Passed variables
C
C Local variables
C
      DATA FA/1.D0,1.D0,2.D0,6.D0,24.D0,12.D1,72.D1,504.D1,4032.D1,
     &     36288.D1,36288.D2,399168.D2,4790016.D2,62270208.D2,
     &     871782912.D2,1307674368.D3,20922789888.D3,355687428096.D3,
     &     6402373705728.D3,121645100408832.D3,243290200817664.D4/
C
C
      I1 = 0
      I2 = IDINT(A+B-C)
      I3 = IDINT(A-AL)
      I4 = IDINT(B+BE)
      I5 = IDINT(C-B+AL)
      I6 = IDINT(C-A-BE)
      ZMIN = MAX0(I1,-I5,-I6)
      ZMAX = MIN0(I2,I3,I4)
      CGK = 0.D0
      IF ( DABS(AL).GT.A ) RETURN
      IF ( DABS(BE).GT.B ) RETURN
      IF ( DABS(GA).GT.C ) RETURN
      IF ( ZMIN.GT.ZMAX ) RETURN
      IF ( DABS(AL+BE-GA).GT.1.D-10 ) RETURN
      I7 = IDINT(A-B+C)
      I8 = IDINT(C+B-A)
      I9 = IDINT(C+B+A)
      I10 = IDINT(A+AL)
      I11 = IDINT(B-BE)
      I12 = IDINT(C+GA)
      I13 = IDINT(C-GA)
      DO Z = ZMIN,ZMAX
         CGK = CGK + (-1.D0)
     &         **Z/(FA(Z)*FA(I2-Z)*FA(I3-Z)*FA(I4-Z)*FA(I5+Z)*FA(I6+Z))
      END DO
      CGK = CGK*DSQRT(FA(I2)*FA(I7)*FA(I8)*FA(I10)*FA(I3)*FA(I4)*FA(I11)
     &      *FA(I12)*FA(I13)*(2.D0*C+1.D0)/FA(I9+1))
C
      END
C*==dmft_ctormt.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
      SUBROUTINE DMFT_CTORMT(YOR,YOI,L)
C.................................................................ctormt
C
C---->    transformation from (ms) to real harmonics basis set
C
      IMPLICIT NONE
C*--DMFT_CTORMT266
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L
      DOUBLE PRECISION YOI(7,7),YOR(7,7)
C
C Local variables
C
      DOUBLE PRECISION SQTWO
C
C*** End of declarations rewritten by SPAG
C
      YOR(1:7,1:7) = 0.0D0
      YOI(1:7,1:7) = 0.0D0
      SQTWO = 1.D0/DSQRT(2.D0)
      IF ( L.EQ.0 ) THEN
         YOR(1,1) = 1.D0
C------- AIL - change to my for p
      ELSE IF ( L.EQ.1 ) THEN
         YOI(1,1) = -SQTWO
         YOI(1,3) = SQTWO
         YOR(2,2) = 1.D0
         YOR(3,1) = SQTWO
         YOR(3,3) = SQTWO
      ELSE IF ( L.EQ.2 ) THEN
         YOI(1,1) = SQTWO
         YOI(1,5) = -SQTWO
         YOI(2,2) = SQTWO
         YOI(2,4) = SQTWO
         YOR(3,3) = 1.D0
         YOR(4,2) = SQTWO
         YOR(4,4) = -SQTWO
         YOR(5,1) = SQTWO
         YOR(5,5) = SQTWO
      ELSE IF ( L.EQ.3 ) THEN
         YOI(1,1) = SQTWO
         YOI(1,7) = SQTWO
         YOI(2,2) = SQTWO
         YOI(2,6) = -SQTWO
         YOI(3,3) = SQTWO
         YOI(3,5) = SQTWO
         YOR(4,4) = 1.D0
         YOR(5,3) = SQTWO
         YOR(5,5) = -SQTWO
         YOR(6,2) = SQTWO
         YOR(6,6) = SQTWO
         YOR(7,1) = SQTWO
         YOR(7,7) = -SQTWO
      END IF
      END
