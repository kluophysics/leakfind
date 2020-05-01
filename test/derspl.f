C*==derspl.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DERSPL(N,X,F,D)
C   ********************************************************************
C   *                                                                  *
C   *   F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N      *
C   *   AND THE SPLINE DERIVATIVES D(I) ARE FOUND.                     *
C   *   THE DIMENSION OF A MUST NOT BE LESS THAN 3*N.                  *
C   *                                                                  *
C   *   sepp redinger 1985                                             *
C   *                                                                  *
C   *   NOTE: the work array A is no more imported as argument         *
C   *         but temporarily allocated                                *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--DERSPL16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DERSPL')
C
C Dummy arguments
C
      INTEGER N
      REAL*8 D(N),F(N),X(N)
C
C Local variables
C
      REAL*8 A(:),H1,H2,P
      INTEGER I,IA_ERR,J,K
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE A
C
      ALLOCATE (A(N*3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: A')
C
      DO I = 2,N
         IF ( X(I).LE.X(I-1) ) THEN
            WRITE (6,99001) I
            CALL STOP_MESSAGE(ROUTINE,'X(I).LE.X(I-1)')
         END IF
      END DO
C
      DO I = 1,N
         J = 2
         IF ( I.NE.1 ) THEN
            J = N - 1
            IF ( I.NE.N ) THEN
               H1 = 1.D0/(X(I)-X(I-1))
               H2 = 1.D0/(X(I+1)-X(I))
               A(3*I-2) = H1
               A(3*I-1) = 2.D0*(H1+H2)
               A(3*I) = H2
               D(I) = 3*(F(I+1)*H2*H2+F(I)*(H1*H1-H2*H2)-F(I-1)*H1*H1)
               CYCLE
            END IF
         END IF
         H1 = 1.D0/(X(J)-X(J-1))
         H2 = 1.D0/(X(J+1)-X(J))
         A(3*I-2) = H1*H1
         A(3*I-1) = H1*H1 - H2*H2
         A(3*I) = -H2*H2
         D(I) = 2.D0*(F(J)*(H2*H2*H2+H1*H1*H1)-F(J+1)*H2*H2*H2-F(J-1)
     &          *H1*H1*H1)
      END DO
C
      P = A(4)/A(1)
      A(5) = A(5) - P*A(2)
      A(6) = A(6) - P*A(3)
      D(2) = D(2) - P*D(1)
      DO I = 3,N
         K = 3*I - 4
         P = A(K+2)/A(K)
         A(K+3) = A(K+3) - P*A(K+1)
         D(I) = D(I) - P*D(I-1)
         IF ( I.EQ.N-1 ) THEN
            P = A(K+5)/A(K)
            A(K+5) = A(K+6) - P*A(K+1)
            A(K+6) = A(K+7)
            D(N) = D(N) - P*D(N-2)
         END IF
      END DO
      D(N) = D(N)/A(3*N-1)
      DO I = 3,N
         J = N + 2 - I
         D(J) = (D(J)-A(3*J)*D(J+1))/A(3*J-1)
      END DO
      D(1) = (D(1)-D(2)*A(2)-D(3)*A(3))/A(1)
      A(1) = 0.D0
C
      DEALLOCATE (A,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
99001 FORMAT ('STOP in <DERSPL>:  mesh point ',I3,' out of order',/)
      END
