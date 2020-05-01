C*==rvec_same.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION RVEC_SAME(N,A,B,TOL)
C   ********************************************************************
C   *                                                                  *
C   *   perform  comparison of 2 vectors  A  and  B                    *
C   *                                                                  *
C   *   N       dimension of A and B                                   *
C   *   A,B     real vectors                                           *
C   *   TOL     max. difference between components A(i) and B(i)       *
C   *                                                                  *
C   *   FALSE   if TOL is exceeded for at least 1 coordinate           *
C   *   TRUE    else                                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 TOL
      REAL*8 A(N),B(N)
      LOGICAL RVEC_SAME
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      RVEC_SAME = .FALSE.
      DO I = 1,N
         IF ( ABS(A(I)-B(I)).GT.TOL ) RETURN
      END DO
      RVEC_SAME = .TRUE.
      END
C*==rvecwgt.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECWGT(N,A,W,B)
C   ********************************************************************
C   *                                                                  *
C   *   weight a  REAL vector A with REAL weights W    B_i = A_i * W_i *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(N),B(N),W(N)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         B(I) = A(I)*W(I)
      END DO
C
      END
C*==rvecnorm.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECNORM(N,A)
C   ********************************************************************
C   *                                                                  *
C   *   normalize REAL vector to 1                     A = A / |A|     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(N)
C
C Local variables
C
      REAL*8 DNRM2
      REAL*8 NORM
C
C*** End of declarations rewritten by SPAG
C
      NORM = DNRM2(N,A,1)
C
      IF ( ABS(NORM).GT.1D-20 ) A(1:N) = A(1:N)/NORM
C
      END
C*==rvecortho.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECORTHO(A,B)
C   ********************************************************************
C   *                                                                  *
C   *   A     real 3x3-matrix containing the original basis vectors    *
C   *         not overwritten                                          *
C   *   B     real 3x3-matrix containing the   NEW    basis vectors    *
C   *         ortho-normal vectors obained via Schmidt-s method        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(3,3),B(3,3)
C
C Local variables
C
      REAL*8 DDOT,DNRM2
      REAL*8 DOT21,DOT31,DOT32,NORM
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CALL RVECCOP(9,A,B)
C
      NORM = 1D0/DNRM2(3,B(1,1),1)
      CALL DSCAL(3,NORM,B(1,1),1)
C
C
      DOT21 = DDOT(3,B(1,2),1,B(1,1),1)
      DO I = 1,3
         B(I,2) = B(I,2) - DOT21*B(I,1)
      END DO
C
      NORM = 1D0/DNRM2(3,B(1,2),1)
      CALL DSCAL(3,NORM,B(1,2),1)
C
      DOT31 = DDOT(3,B(1,3),1,B(1,1),1)
      DOT32 = DDOT(3,B(1,3),1,B(1,2),1)
      DO I = 1,3
         B(I,3) = B(I,3) - DOT31*B(I,1) - DOT32*B(I,2)
      END DO
C
      NORM = 1D0/DNRM2(3,B(1,3),1)
      CALL DSCAL(3,NORM,B(1,3),1)
C
      END
C*==rveclcib.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECLCIB(I,J,K,B,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the scalar-vector operation    A = I*B1 + J*B2 +K*B3  *
C   *                                                                  *
C   *   A       real vector of dimension 3                             *
C   *   I,J,K   integer weights                                        *
C   *   B       3 real basis vectors of dimension 3                    *
C   *           stored as 3x3 matrix                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,J,K
      REAL*8 A(3),B(3,3)
C
C Local variables
C
      INTEGER IC
C
C*** End of declarations rewritten by SPAG
C
      DO IC = 1,3
         A(IC) = I*B(IC,1) + J*B(IC,2) + K*B(IC,3)
      END DO
C
      END
C*==rveclcivb.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECLCIVB(N,IV,B,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the scalar-vector operation    A = SUM(i) I_i * B_i   *
C   *                                                                  *
C   *   N       number of coefficients in vector IV                    *
C   *   A       real vector of dimension 3                             *
C   *   IV      integer weights                                        *
C   *   B       N real basis vectors of dimension 3                    *
C   *           stored as 3xN matrix                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(3),B(3,N)
      INTEGER IV(N)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      A(1:3) = 0D0
C
      DO J = 1,N
         DO I = 1,3
            A(I) = A(I) + B(I,J)*IV(J)
         END DO
      END DO
C
      END
C*==rveclcrb.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECLCRB(U,V,W,B,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the scalar-vector operation    A = U*B1 + V*B2 +W*B3  *
C   *                                                                  *
C   *   A       real vector of dimension 3                             *
C   *   U,V,W   real weights                                           *
C   *   B       3 real basis vectors of dimension 3                    *
C   *           stored as 3x3 matrix                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 U,V,W
      REAL*8 A(3),B(3,3)
C
C Local variables
C
      INTEGER IC
C
C*** End of declarations rewritten by SPAG
C
      DO IC = 1,3
         A(IC) = U*B(IC,1) + V*B(IC,2) + W*B(IC,3)
      END DO
C
      END
C*==rveclcrvb.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECLCRVB(N,RV,B,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the scalar-vector operation    A = SUM(i) R_i * B_i   *
C   *                                                                  *
C   *   N       number of coefficients in vector RV                    *
C   *   A       real vector of dimension 3                             *
C   *   RV      real weights                                           *
C   *   B       N real basis vectors of dimension 3                    *
C   *           stored as 3xN matrix                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(3),B(3,N),RV(N)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      A(1:3) = 0D0
C
      DO J = 1,N
         DO I = 1,3
            A(I) = A(I) + B(I,J)*RV(J)
         END DO
      END DO
C
      END
C*==rvecxpro.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECXPRO(A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the vector-vector operation           C = A x B       *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(3),B(3),C(3)
C
C*** End of declarations rewritten by SPAG
C
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
C
      END
C*==rvecspro.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECSPRO(A,B,C,V)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the vector-vector operation         V = (A x B) * C   *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *   V       spatial product of A, B and C with sign kept           *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 V
      REAL*8 A(3),B(3),C(3)
C
C Local variables
C
      REAL*8 D(3)
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CALL RVECXPRO(A,B,D)
C
      V = 0D0
      DO I = 1,3
         V = V + D(I)*C(I)
      END DO
C
      END
C*==rvecspat.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECSPAT(A,B,C,V,KEY)
C   ********************************************************************
C   *                                                                  *
C   *   spatial product                          V = A * ( B x C )     *
C   *                                                                  *
C   *   A,B,C   real vectors of dimension 3                            *
C   *                                                                  *
C   *   for  KEY = 1:                            V ->  |V|             *
C   *                                                                  *
C   *                                                                  *
C   *   NOTE:  the input is copied to local arrays to avoid            *
C   *          warnings when checking with ftnchek                     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER KEY
      REAL*8 V
      REAL*8 A(3),B(3),C(3)
C
C Local variables
C
      REAL*8 AA(3),BB(3),CC(3),DD(3)
      REAL*8 DDOT
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,3
         AA(I) = A(I)
         BB(I) = B(I)
         CC(I) = C(I)
      END DO
C
      CALL RVECXPRO(BB,CC,DD)
C
      V = DDOT(3,AA,1,DD,1)
C
      IF ( KEY.EQ.1 ) V = ABS(V)
C
      END
C*==rvecexpand.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECEXPAND(A,B,BBINV,C)
C   ********************************************************************
C   *                                                                  *
C   *   expand A with respect to basis vectors B_i                     *
C   *                                                                  *
C   *                   A = sum(i) c_i * B_i                           *
C   *                                                                  *
C   *   with          c_i =  sum(j)  BBINV(i,j) * (A*B_j)              *
C   *                                                                  *
C   *   A     real vector of dimension 3                               *
C   *   B     real 3x3-matrix containing the basis vectors             *
C   *   BBINV real 3x3-matrix with                                     *
C   *         BBINV = BB^(-1)  and  BB(i,j) = (B_i*B_j)                *
C   *   C     real vector of dimension 3 with expansion coefficients   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(3),B(3,3),BBINV(3,3),C(3)
C
C Local variables
C
      REAL*8 ADOTB(3)
      REAL*8 DDOT
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,3
         ADOTB(I) = DDOT(3,A,1,B(1,I),1)
      END DO
C
      DO I = 1,3
         C(I) = 0D0
         DO J = 1,3
            C(I) = C(I) + BBINV(I,J)*ADOTB(J)
         END DO
      END DO
C
      END
C*==rmatmul.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RMATMUL(N,M,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = A * B       *
C   *                                                                  *
C   *   A,B,C   real     SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 R0,EPS
      PARAMETER (R0=0.0D0,EPS=1D-16)
C
C Dummy arguments
C
      INTEGER M,N
      REAL*8 A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      REAL*8 BLJ
      INTEGER I,J,L
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = R0
         END DO
      END DO
C
      DO J = 1,N
         DO L = 1,N
            BLJ = B(L,J)
            IF ( ABS(BLJ).GT.EPS ) THEN
               DO I = 1,N
                  C(I,J) = C(I,J) + A(I,L)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      END
C*==rmatinv.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RMATINV(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           B = 1 / A       *
C   *                                                                  *
C   *   A,B     real  SQUARE  N x N - matrices                        *
C   *   N       dimension of A and B                                   *
C   *   M       array size of A, B    with M >= N                      *
C   *                                                                  *
C   *   on exit A contains the LU-decomposition of A                   *
C   *   NO PIVOTING is used                                            *
C   *                                                                  *
C   *   based on a subroutine by H. Akai                               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      REAL*8 A(M,M),B(M,M)
C
C Local variables
C
      INTEGER I,J,L
      REAL*8 S
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N - 1
         DO J = I + 1,N
            S = -A(J,I)/A(I,I)
C
            DO L = 1,N
               A(J,L) = A(J,L) + S*A(I,L)
            END DO
            A(J,I) = S
         END DO
      END DO
C
      DO I = N,1, - 1
C
         DO J = 1,I - 1
            B(I,J) = A(I,J)
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         DO J = I + 1,N
            B(I,J) = 0D0
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         B(I,I) = 1D0
         DO L = I + 1,N
            B(I,I) = B(I,I) - A(I,L)*B(L,I)
         END DO
C
         S = 1D0/A(I,I)
         DO J = 1,N
            B(I,J) = B(I,J)*S
         END DO
C
      END DO
C
      END
C*==rmat3x3det.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION RMAT3X3DET(A)
C   ********************************************************************
C   *                                                                  *
C   * returns the determinant of a real 3x3 matrix  A                  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A(3,3)
      REAL*8 RMAT3X3DET
C
C*** End of declarations rewritten by SPAG
C
      RMAT3X3DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) + A(1,2)
     &             *(A(2,3)*A(3,1)-A(2,1)*A(3,3)) + A(1,3)
     &             *(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      END
C*==cmatdagger.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATDAGGER(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           B = CONJG(A^T)  *
C   *                                                                  *
C   *   A,B     complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A and B                                   *
C   *   M       array size of A, B with M >= N                         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            B(I,J) = DCONJG(A(J,I))
         END DO
      END DO
C
      END
C*==cmatconjg.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATCONJG(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           B = CONJG(A)    *
C   *                                                                  *
C   *   A,B     complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A and B                                   *
C   *   M       array size of A, B with M >= N                         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            B(I,J) = DCONJG(A(I,J))
         END DO
      END DO
C
      END
C*==cmatmul.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATMUL(N,M,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = A * B       *
C   *                                                                  *
C   *   A,B,C   complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
      REAL*8 EPS
      PARAMETER (EPS=1D-16)
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      COMPLEX*16 BLJ
      INTEGER I,J,L
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = C0
         END DO
      END DO
C
      DO J = 1,N
         DO L = 1,N
            BLJ = B(L,J)
            IF ( ABS(DREAL(BLJ))+ABS(DIMAG(BLJ)).GT.EPS ) THEN
               DO I = 1,N
                  C(I,J) = C(I,J) + A(I,L)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      END
C*==cmatadd.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATADD(N,M,FA,A,FB,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation      C = FA*A + FB* B     *
C   *                                                                  *
C   *   A,B,C   complex  SQUARE  N x N - matrices                      *
C   *   FA,FB   complex weight factors                                 *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 FA,FB
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = FA*A(I,J) + FB*B(I,J)
         END DO
      END DO
C
      END
C*==rinit.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RINIT(N,A)
C   ********************************************************************
C   *                                                                  *
C   *   initialize  REAL   array  A                       A = 0        *
C   *                                                                  *
C   *   A       real array                                             *
C   *   N       dimension of A                                         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 R0
      PARAMETER (R0=0.0D0)
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(N)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         A(I) = R0
      END DO
C
      END
C*==cinit.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CINIT(N,A)
C   ********************************************************************
C   *                                                                  *
C   *   initialize COMPLEX array  A                       A = 0        *
C   *                                                                  *
C   *   A       complex array                                          *
C   *   N       dimension of A                                         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         A(I) = C0
      END DO
C
      END
C*==rveccop.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RVECCOP(N,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   copy real array  A  onto  B                       B = A        *
C   *                                                                  *
C   *   A, B    real arrays                                            *
C   *   N       dimension of A                                         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 A(N),B(N)
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         B(I) = A(I)
      END DO
C
      END
C*==cvecdot.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION CVECDOT(N,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   dot-product of complex vectors   A  and  B       C =  A . B    *
C   *                                                                  *
C   *   A, B    complex arrays                                         *
C   *   N       dimension of A and B                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N),B(N)
      COMPLEX*16 CVECDOT
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CSUM = C0
      DO I = 1,N
         CSUM = CSUM + A(I)*B(I)
      END DO
      CVECDOT = CSUM
C
      END
C*==cvecdotcr.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION CVECDOTCR(N,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   dot-product of vectors   A  and  B               C =  A . B    *
C   *                                                                  *
C   *   A       complex array                                          *
C   *   B       real array                                             *
C   *   N       dimension of A and B                                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N)
      REAL*8 B(N)
      COMPLEX*16 CVECDOTCR
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CSUM = C0
      DO I = 1,N
         CSUM = CSUM + A(I)*B(I)
      END DO
      CVECDOTCR = CSUM
C
      END
C*==cmatcmp.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATCMP(N,M,ISTRUCT,TXTA,A,TXTB,B,THRESH,TOL,SAME)
C   ********************************************************************
C   *                                                                  *
C   *   compare the 2 complex matrices  A  and  B                      *
C   *                                                                  *
C   *   A,B     complex  SQUARE  N x N - matrices                      *
C   *   TXTA/B  text to identify matrix A/B                            *
C   *   N       dimension of A, B                                      *
C   *   M       array size of A, B    with M >= N                      *
C   *   ISTRUCT structure of the matrix corresponding to IREL          *
C   *           0: NONE, 1: NREL/SRA, 2: SP-SREL, 3: REL               *
C   *   THRESH  compare only for  A(i,j) and B(i,j) > THRESH           *
C   *   TOL     accept differences smaller than TOL                    *
C   *   IPRINT  0: no print  1: print differences                      *
C   *   SAME    TRUE:  |(A(i,j)-B(i,j))/A(i,j)|<TOL   for all i,j      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CMATCMP')
C
C Dummy arguments
C
      INTEGER ISTRUCT,M,N
      LOGICAL SAME
      REAL*8 THRESH,TOL
      CHARACTER*(*) TXTA,TXTB
      COMPLEX*16 A(M,M),B(M,M)
C
C Local variables
C
      REAL*8 ADIFFMAX,IM1,IM2,RDIFFMAX,RE1,RE2
      COMPLEX*16 C(:,:)
      INTEGER I,J,LROUT,LTXTA,LTXTB
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C
C
      LROUT = LEN_TRIM(ROUTINE)
      LTXTA = LEN_TRIM(TXTA)
      LTXTB = LEN_TRIM(TXTB)
C
      WRITE (6,99001) ROUTINE(1:LROUT),TXTA(1:LTXTA),TXTB(1:LTXTB),
     &                THRESH,TOL
C
      SAME = .TRUE.
      ADIFFMAX = 0D0
      RDIFFMAX = 0D0
C
      DO I = 1,N
         DO J = 1,N
C
            RE1 = DREAL(A(I,J))
            IM1 = DIMAG(A(I,J))
            RE2 = DREAL(B(I,J))
            IM2 = DIMAG(B(I,J))
C
            IF ( ABS(RE1-RE2).GT.TOL .OR. ABS(IM1-IM2).GT.TOL ) THEN
               WRITE (6,99002) TOL,I,J,A(I,J),B(I,J)
               SAME = .FALSE.
            END IF
C
            ADIFFMAX = MAX(ADIFFMAX,ABS(RE1-RE2),ABS(IM1-IM2))
C
            IF ( ABS(A(I,J)).GT.THRESH ) THEN
               RDIFFMAX = MAX(RDIFFMAX,ABS(1D0-B(I,J)/A(I,J)))
            ELSE IF ( ABS(B(I,J)).GT.THRESH ) THEN
               RDIFFMAX = MAX(RDIFFMAX,ABS(1D0-A(I,J)/B(I,J)))
            END IF
C
         END DO
      END DO
C
      WRITE (6,99003) ADIFFMAX,RDIFFMAX
C
      IF ( .NOT.SAME ) THEN
C
         ALLOCATE (C(M,M))
         C(1:N,1:N) = A(1:N,1:N) - B(1:N,1:N)
C
         CALL CMATSTRUCT('A matrix ',A,N,M,ISTRUCT,ISTRUCT,1,THRESH,6)
C
         CALL CMATSTRUCT('B matrix ',B,N,M,ISTRUCT,ISTRUCT,1,THRESH,6)
C
         CALL CMATSTRUCT('C = A - B',C,N,M,ISTRUCT,ISTRUCT,1,THRESH,6)
C
      END IF
C
99001 FORMAT (/,10X,'<',A,'>   comparison of the matrices ',//,10X,
     &        'A = ',A,/,10X,'B = ',A,//,10X,'for   THRESH =',1PE8.1,
     &        10X,'TOL = ',1PE8.1)
99002 FORMAT (/,5X,'A and B deviate more than ',E10.1,' for    I=',I3,
     &        '  J=',I3,/,5X,' A: ',2E20.12,/,5X,' B: ',2E20.12)
99003 FORMAT (/,10X,'ADIFFMAX =',F20.10,/,10X,'RDIFFMAX =',F20.10,/)
      END
C*==cmatcop.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATCOP(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           B = A           *
C   *                                                                  *
C   *   A,B     complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B                                      *
C   *   M       array size of A, B    with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO J = 1,N
         DO I = 1,N
            B(I,J) = A(I,J)
         END DO
      END DO
C
      END
C*==cmatinv.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATINV(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           B = 1 / A       *
C   *                                                                  *
C   *   A,B     complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A and B                                   *
C   *   M       array size of A, B    with M >= N                      *
C   *                                                                  *
C   *   on exit A contains the LU-decomposition of A                   *
C   *   NO PIVOTING is used                                            *
C   *                                                                  *
C   *   based on a subroutine by H. Akai                               *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0,C1
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M)
C
C Local variables
C
      INTEGER I,J,L
      COMPLEX*16 S
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N - 1
         DO J = I + 1,N
            S = -A(J,I)/A(I,I)
C
            DO L = 1,N
               A(J,L) = A(J,L) + S*A(I,L)
            END DO
            A(J,I) = S
         END DO
      END DO
C
      DO I = N,1, - 1
C
         DO J = 1,I - 1
            B(I,J) = A(I,J)
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         DO J = I + 1,N
            B(I,J) = C0
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         B(I,I) = C1
         DO L = I + 1,N
            B(I,I) = B(I,I) - A(I,L)*B(L,I)
         END DO
C
         S = C1/A(I,I)
         DO J = 1,N
            B(I,J) = B(I,J)*S
         END DO
C
      END DO
C
      END
C*==cmatinv2.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATINV2(N,M,IPIV,W,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           A = 1 / A       *
C   *                                                                  *
C   *   A       matrix to be inverted - overwritten by  inv( A )       *
C   *   W       workspace - contains LU-decoposition of A on exit      *
C   *                                                                  *
C   *   A,W     COMPLEX  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, W                                      *
C   *   M       array size of A, W    with M >= N                      *
C   *   IPIV    workspace - integer array at least of size N           *
C   *                                                                  *
C   *   driver routine using LAPACK routines with PIVOTING             *
C   *                                                                  *
C   *   NOTE: in contrast to <CMATINV3>                                *
C   *         the input matrix A is     OVERWRITTEN                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CMATINV2')
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),W(M,M)
      INTEGER IPIV(N)
C
C Local variables
C
      INTEGER INFO
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------- invert matrix A
C
      CALL ZGETRF(N,N,A,M,IPIV,INFO)
C
      IF ( INFO.NE.0 ) CALL STOP_MESSAGE(ROUTINE,' --> ZGETRF')
C
      CALL ZGETRI(N,A,M,IPIV,W,M*M,INFO)
C
      IF ( INFO.NE.0 ) CALL STOP_MESSAGE(ROUTINE,' --> ZGETRI')
C
      END
C*==cmatinv3.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE CMATINV3(N,M,IPIV,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = 1 / A       *
C   *                                                                  *
C   *   A       matrix to be inverted - unchanged on exit              *
C   *   B       workspace - contains LU-decoposition of A on exit      *
C   *   C       inv( A )                                               *
C   *                                                                  *
C   *   A,B,C   COMPLEX  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *   IPIV    workspace - integer array at least of size N           *
C   *                                                                  *
C   *   driver routine using LAPACK routines with PIVOTING             *
C   *                                                                  *
C   *   NOTE: in contrast to <CMATINV2>                                *
C   *         the input matrix A is NOT OVERWRITTEN                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CMATINV3')
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M),B(M,M),C(M,M)
      INTEGER IPIV(N)
C
C Local variables
C
      INTEGER INFO
C
C*** End of declarations rewritten by SPAG
C
C--------------------------------------------- copy matrix A to matrix C
C
      C(1:N,1:N) = A(1:N,1:N)
C
C------------------------------------------------------- invert matrix C
C
      CALL ZGETRF(N,N,C,M,IPIV,INFO)
C
      IF ( INFO.NE.0 ) CALL STOP_MESSAGE(ROUTINE,' --> ZGETRF')
C
      CALL ZGETRI(N,C,M,IPIV,B,M*M,INFO)
C
      IF ( INFO.NE.0 ) CALL STOP_MESSAGE(ROUTINE,' --> ZGETRI')
C
      END
C*==cmatdet.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      COMPLEX*16 FUNCTION CMATDET(N,M,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the operation           det (A)                       *
C   *                                                                  *
C   *   A       COMPLEX  SQUARE  N x N - matrix                        *
C   *                                                                  *
C   *   N       dimension of A                                         *
C   *   M       array size of A M >= N                                 *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_CONSTANTS,ONLY:C1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M)
C
C Local variables
C
      INTEGER I,INFO,IPIV(:)
      REAL*8 SGN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IPIV
C
      ALLOCATE (IPIV(M))
C
      IPIV = 0
C
      CALL ZGETRF(N,N,A,M,IPIV,INFO)
C
      CMATDET = C1
      DO I = 1,N
         CMATDET = CMATDET*A(I,I)
      END DO
C
      SGN = 1.D0
      DO I = 1,N
         IF ( IPIV(I).NE.I ) SGN = -SGN
      END DO
C
      CMATDET = SGN*CMATDET
C
      END
C*==clogdet.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION CLOGDET(N,M,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the operation     log ( det (A) )  with A = L * U     *
C   *                                                                  *
C   *   A       LU decomposed COMPLEX  SQUARE  N x N - matrix          *
C   *                                                                  *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M)
      COMPLEX*16 CLOGDET
C
C Local variables
C
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CLOGDET = 0D0
      DO I = 1,N
         CLOGDET = CLOGDET + LOG(A(I,I))
      END DO
C
      END
C*==rmatinv2.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RMATINV2(N,M,IPIV,W,A)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           A = 1 / A       *
C   *                                                                  *
C   *   A       matrix to be inverted - overwritten by  inv( A )       *
C   *   W       workspace - contains LU-decoposition of A on exit      *
C   *                                                                  *
C   *   A,W     REAL SQUARE  N x N - matrices                          *
C   *   N       dimension of A, W                                      *
C   *   M       array size of A, W    with M >= N                      *
C   *   IPIV    workspace - integer array at least of size N           *
C   *                                                                  *
C   *   driver routine using LAPACK routines with PIVOTING             *
C   *                                                                  *
C   *   NOTE: in contrast to <RMATINV3>                                *
C   *         the input matrix A is     OVERWRITTEN                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      REAL*8 A(M,M),W(M,M)
      INTEGER IPIV(N)
C
C Local variables
C
      INTEGER INFO
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------- invert matrix A
C
      CALL DGETRF(N,N,A,M,IPIV,INFO)
C
      IF ( INFO.NE.0 ) STOP '<RMATINV2> --> DGETRF'
C
      CALL DGETRI(N,A,M,IPIV,W,M*M,INFO)
C
      IF ( INFO.NE.0 ) STOP '<RMATINV2> --> DGETRI'
C
      END
C*==rmatinv3.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RMATINV3(N,M,IPIV,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = 1 / A       *
C   *                                                                  *
C   *   A       matrix to be inverted - unchanged on exit              *
C   *   B       workspace - contains LU-decoposition of A on exit      *
C   *   C       inv( A )                                               *
C   *                                                                  *
C   *   A,B,C   REAL SQUARE  N x N - matrices                          *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *   IPIV    workspace - integer array at least of size N           *
C   *                                                                  *
C   *   driver routine using LAPACK routines with PIVOTING             *
C   *                                                                  *
C   *   NOTE: in contrast to <RMATINV2>                                *
C   *         the input matrix A is NOT OVERWRITTEN                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      REAL*8 A(M,M),B(M,M),C(M,M)
      INTEGER IPIV(N)
C
C Local variables
C
      INTEGER INFO
C
C*** End of declarations rewritten by SPAG
C
C--------------------------------------------- copy matrix A to matrix C
C
      C(1:N,1:N) = A(1:N,1:N)
C
C------------------------------------------------------- invert matrix C
C
      CALL DGETRF(N,N,C,M,IPIV,INFO)
C
      IF ( INFO.NE.0 ) STOP '<RMATINV3> --> DGETRF'
C
      CALL DGETRI(N,C,M,IPIV,B,M*M,INFO)
C
      IF ( INFO.NE.0 ) STOP '<RMATINV3> --> DGETRI'
C
      END
C*==cmattrc.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION CMATTRC(N,M,A)
C   ********************************************************************
C   *                                                                  *
C   *   takes the trace of a matrix                                    *
C   *                                                                  *
C   *   A       complex  SQUARE  N x N - matrix                        *
C   *   N       dimension of A                                         *
C   *   M       array size of A       with M >= N                      *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER M,N
      COMPLEX*16 A(M,M)
      COMPLEX*16 CMATTRC
C
C Local variables
C
      COMPLEX*16 CSUM
      INTEGER I
C
C*** End of declarations rewritten by SPAG
C
      CSUM = C0
      DO I = 1,N
         CSUM = CSUM + A(I,I)
      END DO
      CMATTRC = CSUM
C
      END
C*==rnon0.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION RNON0(X)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  REAL  argument  X  is numerically <> 0          *
C   *   i.e. whether  ABS(X) > EPS  holds                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EPS
      PARAMETER (EPS=1D-14)
C
C Dummy arguments
C
      REAL*8 X
      LOGICAL RNON0
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(X).GT.EPS ) THEN
         RNON0 = .TRUE.
      ELSE
         RNON0 = .FALSE.
      END IF
C
      END
C*==requ0.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION REQU0(X)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  REAL  argument  X  is numerically = 0           *
C   *   i.e. whether  ABS(X) < EPS  holds                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EPS
      PARAMETER (EPS=1D-14)
C
C Dummy arguments
C
      REAL*8 X
      LOGICAL REQU0
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(X).LT.EPS ) THEN
         REQU0 = .TRUE.
      ELSE
         REQU0 = .FALSE.
      END IF
C
      END
C*==r_lt_eps.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION R_LT_EPS(X,EPS)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  REAL     argument  X  < EPS                     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPS,X
      LOGICAL R_LT_EPS
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( ABS(X).LT.EPS ) THEN
         R_LT_EPS = .TRUE.
      ELSE
         R_LT_EPS = .FALSE.
      END IF
C
      END
C*==r_gt_eps.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION R_GT_EPS(X,EPS)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  REAL     argument  X  > EPS                     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPS,X
      LOGICAL R_GT_EPS
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( ABS(X).GT.EPS ) THEN
         R_GT_EPS = .TRUE.
      ELSE
         R_GT_EPS = .FALSE.
      END IF
C
      END
C*==cnon0.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION CNON0(X)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  COMPLEX  argument  X  is numerically <> 0       *
C   *   i.e. whether  ABS(REAL(X)) > EPS                               *
C   *             OR  ABS(IMAG(X)) > EPS   holds                       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 EPS
      PARAMETER (EPS=1D-14)
C
C Dummy arguments
C
      COMPLEX*16 X
      LOGICAL CNON0
C
C*** End of declarations rewritten by SPAG
C
      IF ( ABS(DREAL(X)).GT.EPS ) THEN
         CNON0 = .TRUE.
      ELSE IF ( ABS(DIMAG(X)).GT.EPS ) THEN
         CNON0 = .TRUE.
      ELSE
         CNON0 = .FALSE.
      END IF
C
      END
C*==c_lt_eps.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION C_LT_EPS(X,EPS)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  COMPLEX  argument  X  < EPS                     *
C   *   i.e. whether  ABS(REAL(X)) < EPS                               *
C   *            AND  ABS(IMAG(X)) < EPS   holds                       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPS
      COMPLEX*16 X
      LOGICAL C_LT_EPS
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( ABS(DREAL(X)).LT.EPS .AND. ABS(DIMAG(X)).LT.EPS ) THEN
         C_LT_EPS = .TRUE.
      ELSE
         C_LT_EPS = .FALSE.
      END IF
C
      END
C*==c_gt_eps.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      FUNCTION C_GT_EPS(X,EPS)
C   ********************************************************************
C   *                                                                  *
C   *   tests whether  COMPLEX  argument  X  > EPS                     *
C   *   i.e. whether  ABS(REAL(X)) > EPS                               *
C   *             OR  ABS(IMAG(X)) > EPS   holds                       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EPS
      COMPLEX*16 X
      LOGICAL C_GT_EPS
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( ABS(DREAL(X)).GT.EPS .OR. ABS(DIMAG(X)).GT.EPS ) THEN
         C_GT_EPS = .TRUE.
      ELSE
         C_GT_EPS = .FALSE.
      END IF
C
      END
C*==rdfdx.f    processed by SPAG 6.70Rc at 08:22 on 26 Apr 2017
      SUBROUTINE RDFDX(F,X,DFDX,N)
C   *******************************************************************
C   *                                                                 *
C   *   differentiate  REAL   function   F(x)   to get  dF / dx       *
C   *   for arbitrary mesh  x(i)  i=1,..,N                            *
C   *   using a 3-point formula                                       *
C   *                                                                 *
C   *******************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      REAL*8 DFDX(N),F(N),X(N)
C
C Local variables
C
      INTEGER I,I0,J
      REAL*8 T0,T1,T2
C
C*** End of declarations rewritten by SPAG
C
      T0(I0,J) = F(I0)*(2*X(J)-X(I0+1)-X(I0+2))
     &           /((X(I0)-X(I0+1))*(X(I0)-X(I0+2)))
C
      T1(I0,J) = F(I0+1)*(2*X(J)-X(I0)-X(I0+2))
     &           /((X(I0+1)-X(I0))*(X(I0+1)-X(I0+2)))
C
      T2(I0,J) = F(I0+2)*(2*X(J)-X(I0)-X(I0+1))
     &           /((X(I0+2)-X(I0))*(X(I0+2)-X(I0+1)))
C
C     FORWARD DIFFERENCE AT THE BEGINNING OF THE TABLE
C
      DFDX(1) = T0(1,1) + T1(1,1) + T2(1,1)
C
C     CENTRAL DIFFERENCE AT THE INTERIOR OF THE TABLE
C
      DO I = 2,N - 1
         DFDX(I) = T0(I-1,I) + T1(I-1,I) + T2(I-1,I)
      END DO
C
C     BACKWARD DIFFERENCE AT THE END OF THE TABLE
C
      DFDX(1) = T0(N-2,N) + T1(N-2,N) + T2(N-2,N)
C
      END
C-----------------------------------------------------------------------
C   ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
C
C   C := alpha*op( A )*op( B ) + beta*C
C
C   for complex  ALPHA  A(M<=LDA,K)  B(K<=LDB,N)  BETA  C(M<=LDC,*)
C
C-----------------------------------------------------------------------
C  ZGETRF( M, N, A, LDA, IPIV, INFO )
C
C  performs LU decomposition   A = P * L * U
C
C  for complex  A(M<=LDA,N)
C  IPIV    N pivot indices
C  INFO    = 0:  successful exit
C
C-----------------------------------------------------------------------
C  ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
C  solves     A * X = B,  A**T * X = B,  or  A**H * X = B
C  TRANS   = 'N'          'T'                'C'
C
C  for complex  A(LDA,N)  B(LDB,NRHS)
C  IPIV    N pivot indices
C          On exit:  the solution matrix X in B
C  INFO    = 0:  successful exit
C
C-----------------------------------------------------------------------
