C*==cinvlu.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVLU(A,B,N,M)
C   ********************************************************************
C   *                                                                  *
C   *  Given a complex matrix a, this returns its inverse matrix.      *
C   *  After calling of this program a is replaced by the lu           *
C   *  decomposition. No pivoting is employed.                         *
C   *  In this version the operations running over the first index of  *
C   *  the matrix come up to the outermost loop, suitable for a big    *
C   *  matrix with rather small size regardint the first index.        *
C   *  coded by H.Akai., Feb. 4, 1992, Osaka                           *
C   *  modified HE 2003                                                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CINVLU17
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
            B(I,J) = (0D0,0D0)
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         B(I,I) = (1D0,0D0)
         DO L = I + 1,N
            B(I,I) = B(I,I) - A(I,L)*B(L,I)
         END DO
C
         S = (1D0,0D0)/A(I,I)
         DO J = 1,N
            B(I,J) = B(I,J)*S
         END DO
C
      END DO
C
      END
C*==cinvlu_spsrel.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVLU_SPSREL(A,B,N,M)
C   ********************************************************************
C   *                                                                  *
C   *  Given a complex matrix a, this returns its inverse matrix.      *
C   *  After calling of this program a is replaced by the lu           *
C   *  decomposition. No pivoting is employed.                         *
C   *  In this version the operations running over the first index of  *
C   *  the matrix come up to the outermost loop, suitable for a big    *
C   *  matrix with rather small size regardint the first index.        *
C   *                                                                  *
C   *  special version of CINVLU to handle case of  NKM < NKMMAX       *
C   *  for spin-polarized relativistic case (IREL=2)                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NLM
      IMPLICIT NONE
C*--CINVLU_SPSREL102
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
      INTEGER I,J,L,NH
      COMPLEX*16 S
C
C*** End of declarations rewritten by SPAG
C
      IF ( M.NE.NKMMAX ) STOP '<CINVLU_SPSREL>   M <> NKMMAX'
      IF ( N.GE.NKMMAX ) STOP '<CINVLU_SPSREL>   N >= NKMMAX'
C
      NH = N/2
      IF ( NH*2.NE.N ) STOP '<CINVLU_SPSREL>   N is ODD'
C
C------------------------------------ remove empty stripes from matrix A
C
      DO I = 1,NH
         DO J = 1,NH
            A(NH+I,J) = A(NLM+I,J)
            A(I,NH+J) = A(I,NLM+J)
            A(NH+I,NH+J) = A(NLM+I,NLM+J)
         END DO
      END DO
C
C---------------------------------------------- invert shrunken matrix A
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
            B(I,J) = (0D0,0D0)
            DO L = I + 1,N
               B(I,J) = B(I,J) - A(I,L)*B(L,J)
            END DO
         END DO
C
         B(I,I) = (1D0,0D0)
         DO L = I + 1,N
            B(I,I) = B(I,I) - A(I,L)*B(L,I)
         END DO
C
         S = (1D0,0D0)/A(I,I)
         DO J = 1,N
            B(I,J) = B(I,J)*S
         END DO
C
      END DO
C
      END
