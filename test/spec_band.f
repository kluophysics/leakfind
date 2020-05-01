C*==bloch.f    processed by SPAG 6.70Rc at 16:35 on 28 Feb 2017
      SUBROUTINE BLOCH(MPM,MPP,MMM,MMP,AK2,AK3,QD,QE,BRUSEP,E2,THETA,JT,
     &                 NT)
C     /****************************************************************/
C     # purpose      : calculates bandstructure: e(kz), e(kparallel)   *
C                                                                      *
C     # subroutines called:                                            *
C       eigen                                                          *
C     /****************************************************************/
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 AK2,AK3,THETA
      COMPLEX*16 E2
      INTEGER JT,NT,QD,QE
      REAL*8 BRUSEP(3)
      COMPLEX*16 MMM(QD,QD),MMP(QD,QD),MPM(QD,QD),MPP(QD,QD)
C
C Local variables
C
      COMPLEX*16 A(:,:),B(:,:),C(:,:),EIG,EIGA(:),EIGB(:),EIGENVAL(:),
     &           EIGKC,EIGKCZ
      REAL*8 AKCP,EOUT,KPARA
      INTEGER I1,I2,ITER(:),J1,J2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE EIGENVAL,A,B,C,EIGA,EIGB,ITER
      ALLOCATE (EIGENVAL(QE),A(QE,QE),B(QE,QE),C(QE,QE),EIGA(QE))
      ALLOCATE (EIGB(QE),ITER(QE))
C
C*** End of declarations rewritten by SPAG
C
      EOUT = DREAL(E2)*27.2116
      AKCP = (BRUSEP(2)*AK2+BRUSEP(3)*AK3)/BRUSEP(1)
      KPARA = 0.3622605*SIN(THETA)*DSQRT(2.0*EOUT)
      WRITE (6,*) 'EOUT,KPARA'
C
      DO I1 = 1,QD
         I2 = I1 + QD
         DO J1 = 1,QD
            J2 = J1 + QD
            A(J1,I1) = C0
            A(J2,I1) = C0
            A(J1,I2) = -MPM(J1,I1)
            A(J2,I2) = MPP(J1,I1)
            B(J1,I1) = MMM(J1,I1)
            B(J2,I1) = -MMP(J1,I1)
            B(J1,I2) = C0
            B(J2,I2) = C0
         END DO
         A(I1,I1) = C1
         B(I2,I2) = C1
      END DO
C
      CALL LZHES(A,B,C,QE)
      CALL LZIT(A,B,C,QE,EIGA,EIGB,ITER)
C
      DO I1 = 1,QE
         IF ( CDABS(EIGA(I1)**2).GT.1.D50*CDABS(EIGB(I1)**2) ) THEN
            WRITE (6,*) 'not converged'
         ELSE
            EIGENVAL(I1) = EIGB(I1)/EIGA(I1)
         END IF
         EIG = EIGENVAL(I1)
         EIGKC = -DCMPLX(0.0,1.0)*CDLOG(EIG)/BRUSEP(1)
         EIGKCZ = EIGKC - DCMPLX(AKCP,0.0D0)
         IF ( DABS(DIMAG(EIGKCZ)).LT.0.001D0 ) THEN
            IF ( NT.EQ.1 ) THEN
               WRITE (NOUT1,99001) DREAL(EIGKCZ),EOUT,KPARA
               WRITE (6,99003) EIGKCZ,EOUT
            ELSE
               WRITE (NOUT1,99002) KPARA,EOUT,JT
               WRITE (6,*) KPARA,EOUT,JT
            END IF
         END IF
      END DO
C
      RETURN
C
99001 FORMAT (3(1x,d12.5))
99002 FORMAT (2(1x,d12.5),2x,i3)
99003 FORMAT (5(1x,d12.5))
      END
