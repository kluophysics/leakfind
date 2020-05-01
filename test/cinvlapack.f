C*==cinvdiablk.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVDIABLK(N,A,UDIA,TAUQ,NKMMAX,NQ,IND0Q,NKMQ,ADD,WGTK)
C   ********************************************************************
C   *                                                                  *
C   *    CINVDIABLK  calculates just the diagonal blocks of a matrix   *
C   *                                                                  *
C   *    it computes the inverse of a COMPLEX matrix  A  using the     *
C   *    LU factorization. This method inverts U and then gets inv(A)  *
C   *    by solving the system  inv(A) = inv(U) * inv(L) for inv(A)    *
C   *    using NO PIVOTING.                                            *
C   *                                                                  *
C   *    The factorization has the form                                *
C   *                               A = L * U                          *
C   *    where L is lower triangular with unit diagonal elements       *
C   *    and   U is upper triangular.                                  *
C   *                                                                  *
C   *                                                                  *
C   *  N       (input) INTEGER                                         *
C   *          The order of the matrix A.  N >= 0.                     *
C   *                                                                  *
C   *  A       (input) COMPLEX*16 array, dimension NxN                 *
C   *          in general this is the KKR-matrix                       *
C   *                                                                  *
C   *  UDIA    (output) COMPLEX*16 array, dimension (N)                *
C   *          as an indermediate result one has   A = L*U             *
C   *          On exit, the diagnal elements of U                      *
C   *                                                                  *
C   *  TAUQ    (output) COMPLEX*16 array, dimension NKMMAXxNKMMAXxNQ   *
C   *          site-diagonal blocks of the original matrix A.          *
C   *                                                                  *
C   *  NKMMAX  (input) INTEGER  dimension of TAUQ                      *
C   *                                                                  *
C   *  NQ      (input) INTEGER  dimension of TAUQ                      *
C   *          specifies the number of sites IQ represented by matrx A *
C   *                                                                  *
C   *  IND0Q   (input) INTEGER  array of dimension NQ                  *
C   *          specifies reference index to start block IQ = 1,..,NQ   *
C   *          with IND0Q(1)=0,  IND0Q(2)=IND0Q(1)+NKMQ(1), etc        *
C   *          allowing for different l-expansion foe each site        *
C   *                                                                  *
C   *  NKMQ    (input) INTEGER  array of dimension NQ                  *
C   *          specifies the size of each block IQ = 1,..,NQ           *
C   *                                                                  *
C   *  ADD     (input) LOGICAL   if TRUE; accumulate TAUQ              *
C   *                                                                  *
C   *  WGTK    (input) REAL      weight for TAUQ when accumulating     *
C   *          e.g. the weight of a k-point within a BZ-integration    *
C   *                                                                  *
C   *  INFO    = 0:  successful exit                                   *
C   *          < 0   = -i, the i-th argument had an illegal value      *
C   *          > 0   =  i, U(i,i) is exactly zero; the matrix          *
C   *                is singular and the program will stop             *
C   *                                                                  *
C   *  NOTE    no test performed on input parameters                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CINVDIABLK59
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D+0,0.0D+0))
C
C Dummy arguments
C
      LOGICAL ADD
      INTEGER N,NKMMAX,NQ
      REAL*8 WGTK
      COMPLEX*16 A(N,N),TAUQ(NKMMAX,NKMMAX,NQ),UDIA(N)
      INTEGER IND0Q(N),NKMQ(N)
C
C Local variables
C
      INTEGER I,II,INFO,IQ,J,JB,JJ,KK,KK0,NB,NKM
      INTEGER ILAENV
      COMPLEX*16 TII,TIJ,TJI
C
C*** End of declarations rewritten by SPAG
C
C     Determine the block size for this environment.
C
      NB = ILAENV(1,'ZGETRF',' ',N,N,-1,-1)
C
      IF ( NB.LE.1 .OR. NB.GE.N ) THEN
C
C        Use unblocked code.
C
         CALL ZGETF2_NO_PIV(N,N,A,N)
C
      ELSE
C
C        Use blocked code.
C
         DO J = 1,N,NB
            JB = MIN(N-J+1,NB)
C
C           Factor diagonal and subdiagonal blocks and test for exact
C           singularity.
C
            CALL ZGETF2_NO_PIV(N-J+1,JB,A(J,J),N)
C
            IF ( J+JB.LE.N ) THEN
C
C              Compute block row of U.
C
               CALL ZTRSM('Left','Lower','No transpose','Unit',JB,
     &                    N-J-JB+1,C1,A(J,J),N,A(J,J+JB),N)
C
C                 Update trailing submatrix.
C
               CALL ZGEMM('No transpose','No transpose',N-J-JB+1,
     &                    N-J-JB+1,JB,-C1,A(J+JB,J),N,A(J,J+JB),N,C1,
     &                    A(J+JB,J+JB),N)
            END IF
         END DO
      END IF
C
C=======================================================================
C       store diagonal elements of matrix  U  for use of Lloyd formula
C
      DO J = 1,N
         UDIA(J) = A(J,J)
      END DO
C
C=======================================================================
C                calculate  inv(U) and inv(L)
C=======================================================================
C  the upper triangular Non-unit structure of inv(U) is same as for U
C  ths applies also to the lower triangular Unit matrices
C  inv(L) and L for which one has implicitly L(i,i) = inv(L)(i,i) = 1
C  accordingly all operations can be done using only array A
C
      CALL ZTRTRI('Upper','Non-unit',N,A,N,INFO)
C
      CALL ZTRTRI('Lower','Unit',N,A,N,INFO)
C
C   calculate   TAUQ = inv(A) = inv(U) * inv(L)
C   and store only diagonal blocks specified by NQ, IND0Q, NKMQ in TAUQ
C   accumulate results if ADD = .TRUE. :   TAUQ = sum WGTK * inv(A)
C
      IF ( ADD ) THEN
C
         DO IQ = 1,NQ
            KK0 = IND0Q(IQ)
            NKM = NKMQ(IQ)
            DO I = 1,NKM
               II = KK0 + I
C
               TII = A(II,II)
               DO KK = II + 1,N
                  TII = TII + A(II,KK)*A(KK,II)
               END DO
               TAUQ(I,I,IQ) = TAUQ(I,I,IQ) + WGTK*TII
C
               DO J = 1,I - 1
                  JJ = KK0 + J
                  TIJ = A(II,II)*A(II,JJ)
                  TJI = A(JJ,II)
                  DO KK = II + 1,N
                     TIJ = TIJ + A(II,KK)*A(KK,JJ)
                     TJI = TJI + A(JJ,KK)*A(KK,II)
                  END DO
                  TAUQ(I,J,IQ) = TAUQ(I,J,IQ) + WGTK*TIJ
                  TAUQ(J,I,IQ) = TAUQ(J,I,IQ) + WGTK*TJI
               END DO
C
            END DO
         END DO
C
      ELSE
C
         DO IQ = 1,NQ
            KK0 = IND0Q(IQ)
            NKM = NKMQ(IQ)
            DO I = 1,NKM
               II = KK0 + I
C
               TII = A(II,II)
               DO KK = II + 1,N
                  TII = TII + A(II,KK)*A(KK,II)
               END DO
               TAUQ(I,I,IQ) = TII
C
               DO J = 1,I - 1
                  JJ = KK0 + J
                  TIJ = A(II,II)*A(II,JJ)
                  TJI = A(JJ,II)
                  DO KK = II + 1,N
                     TIJ = TIJ + A(II,KK)*A(KK,JJ)
                     TJI = TJI + A(JJ,KK)*A(KK,II)
                  END DO
                  TAUQ(I,J,IQ) = TIJ
                  TAUQ(J,I,IQ) = TJI
               END DO
C
            END DO
         END DO
C
      END IF
C
      END
C*==cmatlufact.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CMATLUFACT(N,M,A,B)
C   ********************************************************************
C   *                                                                  *
C   *    CMATLUFACT  calculates just the diagonal blocks of a matrix   *
C   *                                                                  *
C   *    it computes the LU factorization of a COMPLEX matrix  A       *
C   *    using NO PIVOTING.                                            *
C   *                                                                  *
C   *    The factorization has the form                                *
C   *                               B = L * U = A                      *
C   *    where L is lower triangular with unit diagonal elements       *
C   *    and   U is upper triangular.                                  *
C   *                                                                  *
C   *  N       (input) INTEGER                                         *
C   *          The order of the matrix A.  N >= 0.                     *
C   *                                                                  *
C   *  M       (input) INTEGER                                         *
C   *          the array size of A         M >= N                      *
C   *                                                                  *
C   *  A       (input)  COMPLEX*16 array, dimension NxN                *
C   *          unchanged on exit                                       *
C   *                                                                  *
C   *  B       (output) COMPLEX*16 array, dimension NxN                *
C   *          On exit, B contains  L  and  U  with  L(i,i) = 1        *
C   *                                                                  *
C   *  NOTE    no test performed on input parameters                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1
      IMPLICIT NONE
C*--CMATLUFACT252
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
      INTEGER ILAENV
      INTEGER J,JB,NB
C
C*** End of declarations rewritten by SPAG
C
C------------------------------------------------------------ copy B = A
C
      B(1:N,1:N) = A(1:N,1:N)
C
C     Determine the block size for this environment.
C
      NB = ILAENV(1,'ZGETRF',' ',N,N,-1,-1)
C
      IF ( NB.LE.1 .OR. NB.GE.N ) THEN
C
C        Use unblocked code.
C
         CALL ZGETF2_NO_PIV(N,N,B,N)
C
      ELSE
C
C        Use blocked code.
C
         DO J = 1,N,NB
            JB = MIN(N-J+1,NB)
C
C           Factor diagonal and subdiagonal blocks and test for exact
C           singularity.
C
            CALL ZGETF2_NO_PIV(N-J+1,JB,B(J,J),N)
C
            IF ( J+JB.LE.N ) THEN
C
C              Compute block row of U.
C
               CALL ZTRSM('Left','Lower','No transpose','Unit',JB,
     &                    N-J-JB+1,C1,B(J,J),N,B(J,J+JB),N)
C
C                 Update trailing submatrix.
C
               CALL ZGEMM('No transpose','No transpose',N-J-JB+1,
     &                    N-J-JB+1,JB,-C1,B(J+JB,J),N,B(J,J+JB),N,C1,
     &                    B(J+JB,J+JB),N)
            END IF
         END DO
      END IF
C
      END
C*==zgetf2_no_piv.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE ZGETF2_NO_PIV(M,N,A,LDA)
C   ********************************************************************
C   *                                                                  *
C   *  ZGETF2_NO_PIV computes an LU factorization of a                 *
C   *  general m-by-n matrix A using NO PIVOTING.                      *
C   *                                                                  *
C   *  The factorization has the form                                  *
C   *     A = L * U                                                    *
C   *  where L is lower triangular with unit diagonal elements         *
C   *  and   U is upper triangular.                                    *
C   *                                                                  *
C   *  This is the right-looking version of the algorithm.             *
C   *                                                                  *
C   *  Arguments                                                       *
C   *  =========                                                       *
C   *                                                                  *
C   *  M       (input) INTEGER                                         *
C   *          The number of rows of the matrix A.  M >= 0.            *
C   *                                                                  *
C   *  N       (input) INTEGER                                         *
C   *          The number of columns of the matrix A.  N >= 0.         *
C   *                                                                  *
C   *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)      *
C   *          On entry, the m by n matrix to be factored.             *
C   *          On exit, A = (L\U)  with  A(i,i) = U(i,i)               *
C   *                                                                  *
C   *  LDA     (input) INTEGER                                         *
C   *          The leading dimension of the array A.  LDA >= max(1,M). *
C   *                                                                  *
C   *  NOTE    no test performed on input parameters                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C1,C0
      IMPLICIT NONE
C*--ZGETF2_NO_PIV358
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LDA,M,N
      COMPLEX*16 A(LDA,*)
C
C Local variables
C
      INTEGER J,JTOP
C
C*** End of declarations rewritten by SPAG
C
      JTOP = MIN(M,N)
C
      DO J = 1,JTOP
C
C        test for singularity.
C
         IF ( A(J,J).NE.C0 ) THEN
C
C           Compute elements J+1:M of J-th column.
C
            IF ( J.LT.M ) CALL ZSCAL(M-J,C1/A(J,J),A(J+1,J),1)
C
         ELSE
C
            STOP 'in <ZGETF2_NO_PIV>: matrix singular '
C
         END IF
C
C           Update trailing submatrix.
C
C------ use library version of ZGERU instead of reduced version ZGERU_X
C
C        IF ( J.LT.JTOP ) CALL ZGERU_X(M-J,N-J,A(J+1,J),A(J,J+1),LDA,
C    &                                 A(J+1,J+1),LDA)
C
         IF ( J.LT.JTOP ) CALL ZGERU(M-J,N-J,-C1,A(J+1,J),1,A(J,J+1),
     &                               LDA,A(J+1,J+1),LDA)
      END DO
C
C     End of ZGETF2_NO_PIV
C
      END
