C*==one.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C     linalg_nrc.f
C
C     linear algebra subroutines for matrix manipulation
C
C     here with subroutines from numerical recipes
C
C     contains the following subroutines:
C         - one
C         - add, sub
C         - onemx
C         - mult, multp
C         - mvmult
C         - dot, kreuz
C         - inv, inverrzero
C         - zge
C         - zsu
C         - eigendeter
C         - eigen
C         - lzit
C         - lzhes
C
      SUBROUTINE ONE(A,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 A(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            A(I,J) = DCMPLX(0.D0,0.D0)
         END DO
         A(I,I) = DCMPLX(1.D0,0.D0)
      END DO
C
      END
C*==onemx.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE ONEMX(X,XH1,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),XH1(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 CONE
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      CONE = DCMPLX(1.D0,0.D0)
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            XH1(I,J) = -X(I,J)
         END DO
         XH1(I,I) = CONE - X(I,I)
      END DO
C
      CALL INV(XH1,X,DIMVAR)
C
      END
C*==onemx1.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE ONEMX1(X,XH1,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),XH1(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 CONE
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
C
      CONE = DCMPLX(1.D0,0.D0)
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            XH1(I,J) = -X(I,J)
         END DO
         XH1(I,I) = CONE - X(I,I)
      END DO
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            X(I,J) = XH1(I,J)
         END DO
      END DO
C
      END
C*==onepx.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ONEPX(X,XH1,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),XH1(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 CONE
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      CONE = DCMPLX(1.D0,0.D0)
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            XH1(I,J) = X(I,J)
         END DO
         XH1(I,I) = CONE + X(I,I)
      END DO
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            X(I,J) = XH1(I,J)
         END DO
      END DO
C
      END
C*==add.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ADD(X,Y,Z,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),Y(DIMVAR,DIMVAR),Z(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Z(I,J) = X(I,J) + Y(I,J)
         END DO
      END DO
      END
C*==sub.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SUB(X,Y,Z,DIMVAR)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),Y(DIMVAR,DIMVAR),Z(DIMVAR,DIMVAR)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            Z(I,J) = X(I,J) - Y(I,J)
         END DO
      END DO
C
      END
C*==mult.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE MULT(X,Y,Z,DIMVAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 X(DIMVAR,DIMVAR),Y(DIMVAR,DIMVAR),Z(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 C0,C1
C
C*** End of declarations rewritten by SPAG
C
      C1 = DCMPLX(1.0D0,0.0D0)
      C0 = DCMPLX(0.0D0,0.0D0)
      CALL ZGEMM('N','N',DIMVAR,DIMVAR,DIMVAR,C1,X,DIMVAR,Y,DIMVAR,C0,Z,
     &           DIMVAR)
CCC      DO I = 1,DIMVAR
CCC         DO J = 1,DIMVAR
CCC            Z(I,J) = DCMPLX(0.D0,0.D0)
CCC         END DO
CCC      END DO
CCCC
CCC      DO K = 1,DIMVAR
CCC         DO J = 1,DIMVAR
CCC            DO I = 1,DIMVAR
CCC               Z(K,J) = Z(K,J) + X(K,I)*Y(I,J)
CCC            END DO
CCC         END DO
CCC      END DO
C
      END
C*==multp.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE MULTP(X,Y,Z,DIM1,DIM2)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIM1,DIM2
      COMPLEX*16 X(DIM1,DIM2),Y(DIM2),Z(DIM1,DIM2)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
C
C
      DO I = 1,DIM1
         DO J = 1,DIM2
            Z(I,J) = DCMPLX(0.D0,0.D0)
         END DO
      END DO
C
      DO I = 1,DIM1
         DO J = 1,DIM2
            Z(I,J) = X(I,J)*Y(J)
         END DO
      END DO
C
      END
C*==mvmult.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE MVMULT(A,X,B,DIMVAR)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 A(DIMVAR,DIMVAR),B(DIMVAR),X(DIMVAR)
C
C Local variables
C
      INTEGER I,J
      COMPLEX*16 SUMVAR
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,DIMVAR
         SUMVAR = DCMPLX(0.D0,0.D0)
         DO J = 1,DIMVAR
            SUMVAR = SUMVAR + A(I,J)*X(J)
         END DO
         B(I) = SUMVAR
      END DO
C
      END
C*==dot.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DOT(VA,VB,P,G)
C
      USE MOD_SPEC,ONLY:LG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER G
      REAL*8 P(LG),VA(3),VB(LG,3)
C
C*** End of declarations rewritten by SPAG
C
      P(G) = VA(1)*VB(G,1) + VA(2)*VB(G,2) + VA(3)*VB(G,3)
C
      END
C*==kreuz.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE KREUZ(A,B,C,G)
C
      USE MOD_SPEC,ONLY:LG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER G
      REAL*8 A(LG,3),B(LG,3),C(LG,3)
C
C*** End of declarations rewritten by SPAG
C
C
C
      C(G,1) = A(G,2)*B(G,3) - A(G,3)*B(G,2)
      C(G,2) = ((-A(G,1)*B(G,3))) + A(G,3)*B(G,1)
      C(G,3) = A(G,1)*B(G,2) - A(G,2)*B(G,1)
C
      END
C*==inv_woz.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE INV_WOZ(XIN,XOUT,DIMVAR)
C
      USE MOD_SPEC,ONLY:EPS12
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 XIN(DIMVAR,DIMVAR),XOUT(DIMVAR,DIMVAR)
C
C Local variables
C
      REAL*8 EMACH
      INTEGER I,INTA(:),J
      COMPLEX*16 X(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE X,INTA
      ALLOCATE (X(DIMVAR,DIMVAR),INTA(DIMVAR))
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     # purpose      : matrix inversion                                *
C     # uses nrc subroutines:    zge, zsu                              *
C     /****************************************************************/
C     /* input */
C     /* output */
C     /* local */
C
      EMACH = EPS12
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            X(I,J) = XIN(I,J)
         END DO
      END DO
C
      CALL ZGE(X,INTA,DIMVAR,DIMVAR,EMACH)
      CALL ONE(XOUT,DIMVAR)
      DO I = 1,DIMVAR
         CALL ZSU(X,INTA,XOUT(1,I),DIMVAR,DIMVAR,EMACH)
      END DO
C
      END
C*==inv.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE INV(XIN,XOUT,DIMVAR)
C
      USE MOD_SPEC,ONLY:EPS12,EPS16,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 XIN(DIMVAR,DIMVAR),XOUT(DIMVAR,DIMVAR)
C
C Local variables
C
      REAL*8 EMACH
      INTEGER I,INTA(:),J
      COMPLEX*16 X(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE X,INTA
      ALLOCATE (X(DIMVAR,DIMVAR),INTA(DIMVAR))
C
C*** End of declarations rewritten by SPAG
C
C
C     /****************************************************************/
C     # purpose      : matrix inversion                                *
C                      with error message and                          *
C                      removal of small components                     *
C     # uses nrc subroutines:    zge, zsu (ludcmp, lubksb)             *
C     /****************************************************************/
C     /* input */
C     /* output */
C     /* local */
C
C     dummy to avoid warning in nrc version
C
      EMACH = EPS12
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            X(I,J) = XIN(I,J)
         END DO
      END DO
C
      CALL ZGE(X,INTA,DIMVAR,DIMVAR,EMACH)
      CALL ONE(XOUT,DIMVAR)
      DO I = 1,DIMVAR
         CALL ZSU(X,INTA,XOUT(1,I),DIMVAR,DIMVAR,EMACH)
      END DO
C
C     remove small things
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            IF ( CDABS(XOUT(I,J)).LT.EPS16 ) XOUT(I,J) = CZERO
         END DO
      END DO
C
      END
C*==zge.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
C
C
      SUBROUTINE ZGE(A,INTA,NR,NC,EMACH)
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EMACH
      INTEGER NC,NR
      COMPLEX*16 A(NR,NC)
      INTEGER INTA(NC)
C
C Local variables
C
      COMPLEX*16 DUM,YR
      INTEGER I,II,IN,J,K,N
C
C*** End of declarations rewritten by SPAG
C
C     /****************************************************************/
C	# purpose       : lu decomposition                               *
C                                                                      *
C     # parameter   : a = input (complex*16) matrix,                   *
C                     will be overwritten by a triangularized          *
C                     matrix, to be transmitted to subroutine zsu.     *
C                     inta = storage for permutation of matrix columns *
C                     (and rows), to be transmitted to subroutine zsu. *
C                     nr = first dimension of a (.ge.nc).              *
C                     nc = order of a.                                 *
C                     emach = machine accuracy.                        *
C                                                                      *
C	# note:         corresponds to nrc 2.0 ludcmp                    *
C     /****************************************************************/
C
C
      N = NC
      DO II = 2,N
         I = II - 1
         YR = A(I,I)
         IN = I
         DO J = II,N
            IF ( CDABS(YR)-CDABS(A(J,I)).LT.0. ) THEN
               YR = A(J,I)
               IN = J
            END IF
         END DO
         INTA(I) = IN
         IF ( IN.NE.I ) THEN
            DO J = I,N
               DUM = A(I,J)
               A(I,J) = A(IN,J)
               A(IN,J) = DUM
            END DO
         END IF
         IF ( CDABS(YR)-EMACH.GT.0. ) THEN
            DO J = II,N
               IF ( CDABS(A(J,I))-EMACH.GT.0. ) THEN
                  A(J,I) = A(J,I)/YR
                  DO K = II,N
                     A(J,K) = A(J,K) - A(I,K)*A(J,I)
                  END DO
               END IF
            END DO
         END IF
      END DO
C
      END
C*==zsu.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE ZSU(A,INTA,X,NR,NC,EMACH)
C     /****************************************************************/
C     # purpose       : lu backsubstitution                            *
C                                                                      *
C     # parameter   : a = input matrix, prepared by subroutine zge.    *
C                     inta = input permutation from subroutine zge.    *
C                     x = input constan and output resulting vector.   *
C                     nr, nc, emach: like zge.                         *
C                                                                      *
C	# note:         corresponds to nrc 2.0 lubksb                  *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 EMACH
      INTEGER NC,NR
      COMPLEX*16 A(NR,NC),X(NC)
      INTEGER INTA(NC)
C
C Local variables
C
      COMPLEX*16 DUM
      INTEGER I,II,IJ,IN,J,N
C
C*** End of declarations rewritten by SPAG
C
      N = NC
      DO II = 2,N
         I = II - 1
         IF ( INTA(I).NE.I ) THEN
            IN = INTA(I)
            DUM = X(IN)
            X(IN) = X(I)
            X(I) = DUM
         END IF
         DO J = II,N
            IF ( CDABS(A(J,I)).GT.EMACH ) X(J) = X(J) - A(J,I)*X(I)
         END DO
      END DO
      DO II = 1,N
         I = N - II + 1
         IJ = I + 1
         IF ( I.NE.N ) THEN
            DO J = IJ,N
               X(I) = X(I) - A(I,J)*X(J)
            END DO
         END IF
         IF ( CDABS(A(I,I))-EMACH*1.0D-5.LT.0. ) A(I,I)
     &        = EMACH*1.0D-5*(1.D0,1.D0)
         X(I) = X(I)/A(I,I)
      END DO
C
      END
C*==diagonal.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE DIAGONAL(A,B,C,D,DIMVAR)
C     /****************************************************************/
C     # purpose       : prepares matrix diagonalisation                *
C                       calculates the eigenvalues and eigenvectors of *
C                       a complex matrix                               *
C                       the eigenvalues build a diagonal matrix and    *
C                       the eigenvectors provide the unitary           *
C                       transformation matrix                          *
C                                                                      *
C     # parameter     : a   input matrix                               *
C                       b   eigenvalues in b(i,i)                      *
C                       c   left eigenvector v of                      *
C                           a * v  = lambda * v                        *
C                       d   c^-1 inverse of c                          *
C                                                                      *
C                                                                      *
C     uses lapack routine : zgeev  (parameter see linalg.f)            *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR
      COMPLEX*16 A(DIMVAR,DIMVAR),B(DIMVAR,DIMVAR),C(DIMVAR,DIMVAR),
     &           D(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 AIN(:,:),VL(:,:),VR(:,:),W(:),WORK(:)
      INTEGER DIM2,I,INFO,J
      REAL*8 RWORK(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W,VL,VR,AIN,WORK,RWORK
      ALLOCATE (W(DIMVAR),VL(DIMVAR,DIMVAR),VR(DIMVAR,DIMVAR),
     &          AIN(DIMVAR,DIMVAR))
      ALLOCATE (WORK(2*DIMVAR),RWORK(2*DIMVAR))
C
      INFO = 1
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            AIN(I,J) = A(I,J)
            B(I,J) = DCMPLX(0.D0,0.D0)
         END DO
      END DO
C
      DIM2 = 2*DIMVAR
      CALL ZGEEV('v','n',DIMVAR,AIN,DIMVAR,W,VL,DIMVAR,VR,DIMVAR,WORK,
     &           DIM2,RWORK,INFO)
C                 jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr,
C                 work, lwork, rwork, info )
C
      IF ( INFO.NE.0 ) THEN
         WRITE (NOUT1,99001) INFO
         WRITE (*,99001) INFO
         STOP
      END IF
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            B(I,J) = CZERO
         END DO
      END DO
C
      DO I = 1,MQD
         B(I,I) = CDSQRT(W(I))
      END DO
C
      DO I = 1,DIMVAR
         DO J = 1,DIMVAR
            C(I,J) = VL(I,J)
            D(I,J) = DCONJG(VL(J,I))
         END DO
      END DO
C
      RETURN
99001 FORMAT (1x,'error in diagonal zgeev: info=',i3)
      END
C*==eigendeter.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE EIGENDETER(A,B,DIMVAR,DET)
C     /****************************************************************/
C     # purpose       : calculates the complexdeterminant of a general *
C                       complex*16 matrix a                            *
C                                                                      *
C     # calls:        eigen                                            *
C     /****************************************************************/
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DET
      INTEGER DIMVAR
      COMPLEX*16 A(DIMVAR,DIMVAR),B(DIMVAR,DIMVAR)
C
C Local variables
C
      COMPLEX*16 EIGA(:),EIGB(:),EIGENVAL(:)
      INTEGER I,NEIG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE EIGENVAL,EIGA,EIGB
      ALLOCATE (EIGENVAL(DIMVAR),EIGA(DIMVAR),EIGB(DIMVAR))
C
C*** End of declarations rewritten by SPAG
C
      CALL EIGEN(A,B,DIMVAR,EIGA,EIGB,EIGENVAL,NEIG)
C
      DET = DCMPLX(1.D0,0.D0)
      DO I = 1,NEIG
         IF ( CDABS(EIGENVAL(I)).GT.0.D0 ) DET = DET*EIGENVAL(I)
      END DO
C
      END
C*==eigen.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE EIGEN(A,B,DIMVAR,EIGA,EIGB,EIGENVAL,NEIG)
C     /****************************************************************/
C     # purpose       : calculates the eigenvalues and, if required,   *
C                       the eigenvectors of the complex generalized    *
C                       eigenproblem:      a*x = lambda*b*x            *
C                                                                      *
C     # using the nrc routines:       lzit, lzhes                      *
C     /****************************************************************/
C
      USE MOD_SPEC_OUTPUT,ONLY:NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DIMVAR,NEIG
      COMPLEX*16 A(DIMVAR,DIMVAR),B(DIMVAR,DIMVAR),EIGA(DIMVAR),
     &           EIGB(DIMVAR),EIGENVAL(DIMVAR)
C
C Local variables
C
      COMPLEX*16 AIN(:,:),BIN(:,:),C(:,:)
      INTEGER I,INFO,IPVT(:),J
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE C,AIN,BIN,IPVT
      ALLOCATE (C(DIMVAR,DIMVAR),AIN(DIMVAR,DIMVAR),BIN(DIMVAR,DIMVAR),
     &          IPVT(DIMVAR))
C
      INFO = 0
C
      DO I = 1,DIMVAR
         EIGENVAL(I) = DCMPLX(0.D0,0.D0)
         DO J = 1,DIMVAR
            AIN(I,J) = A(I,J)
            BIN(I,J) = B(I,J)
         END DO
      END DO
C
      CALL LZHES(AIN,BIN,C,DIMVAR)
      CALL LZIT(AIN,BIN,C,DIMVAR,EIGA,EIGB,IPVT)
C
      IF ( INFO.NE.0 ) THEN
         WRITE (NOUT1,99001) INFO,DIMVAR
C          write (*,1000) info, dim
         IF ( INFO.LT.0 .OR. INFO.GE.DIMVAR ) THEN
            WRITE (NOUT1,99002)
            WRITE (*,99002)
            STOP
         END IF
      END IF
C
C     eigenvalues are at least from info+1 partially correct
      NEIG = 0
      DO I = INFO + 1,DIMVAR
         IF ( CDABS(EIGB(I)).GT.1.D38*CDABS(EIGA(I)) ) THEN
            WRITE (NOUT1,99003)
            EIGENVAL(NEIG) = DCMPLX(0.D0,0.D0)
         ELSE
            NEIG = NEIG + 1
            EIGENVAL(NEIG) = EIGA(I)/EIGB(I)
         END IF
      END DO
C
      RETURN
C
99001 FORMAT (1x,'warning in eigen: zgegv info:',i3,' dim=',i3)
99002 FORMAT (1x,'fatal error in eigen: --> stop')
99003 FORMAT (1x,'betr(eigb) > 1.d38 * betr(eiga)')
      END
C*==lzit.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE LZIT(A,B,C,N,EIGA,EIGB,ITER)
C
      USE MOD_SPEC,ONLY:CZERO
      USE MOD_FILES,ONLY:CDUMMY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N,N),B(N,N),C(N,N),EIGA(N),EIGB(N)
      INTEGER ITER(N)
C
C Local variables
C
      COMPLEX*16 ALFM,ANM1M1,ANNM1,BETM,CTMP,D,DEN,NUM,S,SL,W,Y,Z
      REAL*8 ANI,ANORM,BNI,BNORM,D0,D1,D2,E0,E1,EPSA,EPSB,F,R,SS
      INTEGER I,ITS,J,K,L,L1,LOR1,LZ,M,MB,NL,NM1,NN,NNORN
C
C*** End of declarations rewritten by SPAG
C
      NN = N
      ANORM = 0.
      BNORM = 0.
      ANNM1 = CZERO
      CDUMMY = ANNM1
      DO I = 1,N
         ANI = 0.0D0
         IF ( I.NE.1 ) THEN
            Y = A(I,I-1)
            ANI = ANI + DABS(DREAL(Y)) + DABS(DIMAG(Y))
         END IF
         BNI = 0.
         DO J = I,N
            ANI = ANI + DABS(DREAL(A(I,J))) + DABS(DIMAG(A(I,J)))
            BNI = BNI + DABS(DREAL(B(I,J))) + DABS(DIMAG(B(I,J)))
         END DO
         ANORM = DMAX1(ANI,ANORM)
         BNORM = DMAX1(BNI,BNORM)
      END DO
      IF ( ANORM.EQ.0.0D0 ) ANORM = 1.0D0
      IF ( BNORM.EQ.0.0D0 ) BNORM = 1.0D0
C      IF ( ABS(ANORM).LE.1.0D-16 ) ANORM = 1.0D0
C      IF ( ABS(BNORM).LE.1.0D-16 ) BNORM = 1.0D0
      EPSB = BNORM
      EPSA = ANORM
 100  CONTINUE
      EPSA = EPSA/2.0D0
      EPSB = EPSB/2.0D0
      F = ANORM + EPSA
      IF ( F.GT.ANORM ) GOTO 100
      IF ( N.LE.1 ) GOTO 500
 200  CONTINUE
      ITS = 0
      NM1 = NN - 1
 300  CONTINUE
      D2 = DABS(DREAL(A(NN,NN))) + DABS(DIMAG(A(NN,NN)))
      DO LZ = 2,NN
         L = NN + 2 - LZ
         SS = D2
         Y = A(L-1,L-1)
         D2 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
         SS = SS + D2
         Y = A(L,L-1)
         R = SS + DABS(DREAL(Y)) + DABS(DIMAG(Y))
         IF ( R.EQ.SS ) GOTO 400
C         IF ( ABS(R-SS).LE.1.0D-16 ) GOTO 400
      END DO
      L = 1
 400  CONTINUE
      IF ( L.NE.NN ) THEN
         IF ( ITS.GE.30 ) THEN
            ITER(NN) = -1
            IF ( DABS(DREAL(A(NN,NM1)))+DABS(DIMAG(A(NN,NM1)))
     &           .GT.0.8D0*DABS(DREAL(ANNM1))+DABS(DIMAG(ANNM1)) )
     &           RETURN
         END IF
         IF ( ITS.NE.10 .AND. ITS.NE.20 ) THEN
            ANNM1 = A(NN,NM1)
            ANM1M1 = A(NM1,NM1)
            S = A(NN,NN)*B(NM1,NM1) - ANNM1*B(NM1,NN)
            W = ANNM1*B(NN,NN)*(A(NM1,NN)*B(NM1,NM1)-B(NM1,NN)*ANM1M1)
            Y = (ANM1M1*B(NN,NN)-S)/2.
            Z = CDSQRT(Y*Y+W)
            IF ( DREAL(Z).NE.0.0D0 .OR. DIMAG(Z).NE.0.0D0 ) THEN
               CTMP = Y/Z
               D0 = DREAL(CTMP)
               IF ( D0.LT.0.0D0 ) Z = -Z
            END IF
            DEN = (Y+Z)*B(NM1,NM1)*B(NN,NN)
            IF ( DREAL(DEN).EQ.0.0D0 .AND. DIMAG(DEN).EQ.0.0D0 )
     &           DEN = EPSA
C            IF ( ABS(DREAL(DEN)).LE.1.0D-16 .AND. ABS(DIMAG(DEN))
C     &           .LE.1.0D-16 ) DEN = EPSA
            NUM = (Y+Z)*S - W
         ELSE
            Y = A(NM1,NN-2)
            NUM = DCMPLX(DABS(DREAL(ANNM1))+DABS(DIMAG(ANNM1)),
     &            DABS(DREAL(Y))+DABS(DIMAG(Y)))
            DEN = (1.0D0,0.0D0)
         END IF
         IF ( NN.NE.L+1 ) THEN
            D2 = DABS(DREAL(A(NM1,NM1))) + DABS(DIMAG(A(NM1,NM1)))
            E1 = DABS(DREAL(ANNM1)) + DABS(DIMAG(ANNM1))
            D1 = DABS(DREAL(A(NN,NN))) + DABS(DIMAG(A(NN,NN)))
            NL = NN - (L+1)
            DO MB = 1,NL
               M = NN - MB
               E0 = E1
               Y = A(M,M-1)
               E1 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
               D0 = D1
               D1 = D2
               Y = A(M-1,M-1)
               D2 = DABS(DREAL(Y)) + DABS(DIMAG(Y))
               Y = A(M,M)*DEN - B(M,M)*NUM
               D0 = (D0+D1+D2)*(DABS(DREAL(Y))+DABS(DIMAG(Y)))
               E0 = E0*E1*(DABS(DREAL(DEN))+DABS(DIMAG(DEN))) + D0
               IF ( E0.EQ.D0 ) GOTO 450
C               IF ( ABS(E0-D0).LE.1.0D-16 ) GOTO 450
            END DO
         END IF
         M = L
 450     CONTINUE
         ITS = ITS + 1
         W = A(M,M)*DEN - B(M,M)*NUM
         Z = A(M+1,M)*DEN
         D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
         D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
         LOR1 = L
         NNORN = NN
         LOR1 = 1
         NNORN = N
         DO I = M,NM1
            J = I + 1
            IF ( I.NE.M ) THEN
               W = A(I,I-1)
               Z = A(J,I-1)
               D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
               D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
               IF ( D1.EQ.0.0D0 ) EXIT
C               IF ( ABS(D1).LE.1.0D-16 ) EXIT
            END IF
            IF ( D2.LE.D1 ) THEN
               DO K = I,NNORN
                  Y = A(I,K)
                  A(I,K) = A(J,K)
                  A(J,K) = Y
                  Y = B(I,K)
                  B(I,K) = B(J,K)
                  B(J,K) = Y
               END DO
               IF ( I.GT.M ) A(I,I-1) = A(J,I-1)
               IF ( D2.EQ.0.0D0 ) GOTO 460
C               IF ( ABS(D2).LE.1.0D-16 ) GOTO 460
               Y = DCMPLX(DREAL(W)/D1,DIMAG(W)/D1)
     &             /DCMPLX(DREAL(Z)/D1,DIMAG(Z)/D1)
            ELSE
               Y = DCMPLX(DREAL(Z)/D2,DIMAG(Z)/D2)
     &             /DCMPLX(DREAL(W)/D2,DIMAG(W)/D2)
            END IF
            DO K = I,NNORN
               A(J,K) = A(J,K) - Y*A(I,K)
               B(J,K) = B(J,K) - Y*B(I,K)
            END DO
            IF ( I.GT.M ) A(J,I-1) = (0.0D0,0.0D0)
 460        CONTINUE
            Z = B(J,I)
            W = B(J,J)
            D2 = DABS(DREAL(W)) + DABS(DIMAG(W))
            D1 = DABS(DREAL(Z)) + DABS(DIMAG(Z))
            IF ( D1.EQ.0.0D0 ) EXIT
C            IF ( ABS(D1).LE.1.0D-16 ) EXIT
            IF ( D2.LE.D1 ) THEN
               DO K = LOR1,J
                  Y = A(K,J)
                  A(K,J) = A(K,I)
                  A(K,I) = Y
                  Y = B(K,J)
                  B(K,J) = B(K,I)
                  B(K,I) = Y
               END DO
               IF ( I.NE.NM1 ) THEN
                  Y = A(J+1,J)
                  A(J+1,J) = A(J+1,I)
                  A(J+1,I) = Y
               END IF
               DO K = 1,N
                  Y = C(K,J)
                  C(K,J) = C(K,I)
                  C(K,I) = Y
               END DO
               IF ( D2.EQ.0.0D0 ) CYCLE
C               IF ( ABS(D2).LE.1.0D-16 ) CYCLE
               Z = DCMPLX(DREAL(W)/D1,DIMAG(W)/D1)
     &             /DCMPLX(REAL(Z)/D1,DIMAG(Z)/D1)
            ELSE
               Z = DCMPLX(DREAL(Z)/D2,DIMAG(Z)/D2)
     &             /DCMPLX(DREAL(W)/D2,DIMAG(W)/D2)
            END IF
            DO K = LOR1,J
               A(K,I) = A(K,I) - Z*A(K,J)
               B(K,I) = B(K,I) - Z*B(K,J)
            END DO
            B(J,I) = (0.0D0,0.0D0)
            IF ( I.LT.NM1 ) A(I+2,I) = A(I+2,I) - Z*A(I+2,J)
            DO K = 1,N
               C(K,I) = C(K,I) - Z*C(K,J)
            END DO
         END DO
         GOTO 300
      END IF
 500  CONTINUE
      EIGA(NN) = A(NN,NN)
      EIGB(NN) = B(NN,NN)
      IF ( NN.EQ.1 ) THEN
         M = N
      ELSE
         ITER(NN) = ITS
         NN = NM1
         IF ( NN.GT.1 ) GOTO 200
         ITER(1) = 0
         GOTO 500
      END IF
 600  CONTINUE
      ALFM = A(M,M)
      BETM = B(M,M)
      B(M,M) = (1.0D0,0.0D0)
      L = M - 1
      IF ( L.EQ.0 ) GOTO 800
 700  CONTINUE
      L1 = L + 1
      SL = 0.
      DO J = L1,M
         SL = SL + (BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
      END DO
      Y = BETM*A(L,L) - ALFM*B(L,L)
      IF ( DREAL(Y).EQ.0.0D0 .AND. DIMAG(Y).EQ.0.0D0 ) Y = (EPSA+EPSB)
     &     /2.0D0
C      IF ( ABS(DREAL(Y)).LE.1.0D-16 .AND. ABS(DIMAG(Y)).LE.1.0D-16 )
C     &     Y = (EPSA+EPSB)/2.0D0
      B(L,M) = -SL/Y
      L = L - 1
 800  CONTINUE
      IF ( L.GT.0 ) GOTO 700
      M = M - 1
      IF ( M.GT.0 ) GOTO 600
      M = N
 900  CONTINUE
      DO I = 1,N
         S = 0.
         DO J = 1,M
            S = S + C(I,J)*B(J,M)
         END DO
         C(I,M) = S
      END DO
      M = M - 1
      IF ( M.GT.0 ) GOTO 900
      M = N
 1000 CONTINUE
      SS = 0.
      DO I = 1,N
         R = DABS(DREAL(C(I,M))) + DABS(DIMAG(C(I,M)))
         IF ( R.GE.SS ) THEN
            SS = R
            D = C(I,M)
         END IF
      END DO
      IF ( SS.NE.0.0D0 ) THEN
C      IF ( ABS(SS).GT.1.0D-16 ) THEN
         DO I = 1,N
            C(I,M) = C(I,M)/D
         END DO
      END IF
      M = M - 1
      IF ( M.GT.0 ) GOTO 1000
      END
C*==lzhes.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE LZHES(A,B,C,N)
C
      USE MOD_SPEC,ONLY:CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 A(N,N),B(N,N),C(N,N)
C
C Local variables
C
      REAL*8 D,F
      INTEGER I,II,IM1,IMJ,IP1,J,JM2,JP1,K,NM1,NM2
      COMPLEX*16 W,Y,Z
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,N
         DO J = 1,N
            C(I,J) = CZERO
         END DO
      END DO
C
      NM1 = N - 1
      DO I = 1,NM1
         D = 0.0D0
         IP1 = I + 1
         DO K = IP1,N
            Y = B(K,I)
            F = DABS(DREAL(Y)) + DABS(DIMAG(Y))
            IF ( F.GT.D ) THEN
               D = F
               II = K
            END IF
         END DO
C         IF ( D.NE.0.0D0 ) THEN
         IF ( ABS(D).GT.1.0D-16 ) THEN
            Y = B(I,I)
            IF ( D.GT.DABS(DREAL(Y))+DABS(DIMAG(Y)) ) THEN
               DO J = 1,N
                  Y = A(I,J)
                  A(I,J) = A(II,J)
                  A(II,J) = Y
               END DO
               DO J = I,N
                  Y = B(I,J)
                  B(I,J) = B(II,J)
                  B(II,J) = Y
               END DO
            END IF
            DO J = IP1,N
               Y = B(J,I)/B(I,I)
C               IF ( DREAL(Y).NE.0.0D0 .OR. DIMAG(Y).NE.0.0D0 ) THEN
               IF ( ABS(DREAL(Y)).GT.1.0D-16 .OR. ABS(DIMAG(Y))
     &              .GT.1.0D-16 ) THEN
                  DO K = 1,N
                     A(J,K) = A(J,K) - Y*A(I,K)
                  END DO
                  DO K = IP1,N
                     B(J,K) = B(J,K) - Y*B(I,K)
                  END DO
               END IF
            END DO
            B(IP1,I) = (0.0D0,0.0D0)
         END IF
      END DO
      DO I = 1,N
         DO J = 1,N
            C(I,J) = (0.0D0,0.0D0)
         END DO
         C(I,I) = (1.0D0,0.0D0)
      END DO
      NM2 = N - 2
      IF ( NM2.GE.1 ) THEN
         DO J = 1,NM2
            JM2 = NM1 - J
            JP1 = J + 1
            DO II = 1,JM2
               I = N + 1 - II
               IM1 = I - 1
               IMJ = I - J
               W = A(I,J)
               Z = A(IM1,J)
               IF ( DABS(DREAL(W))+DABS(DIMAG(W)).GT.DABS(DREAL(Z))
     &              +DABS(DIMAG(Z)) ) THEN
                  DO K = J,N
                     Y = A(I,K)
                     A(I,K) = A(IM1,K)
                     A(IM1,K) = Y
                  END DO
                  DO K = IM1,N
                     Y = B(I,K)
                     B(I,K) = B(IM1,K)
                     B(IM1,K) = Y
                  END DO
               END IF
               Z = A(I,J)
C               IF ( DREAL(Z).NE.0.0D0 .OR. DIMAG(Z).NE.0.0D0 ) THEN
               IF ( ABS(DREAL(Z)).GT.1.0D-16 .OR. ABS(DIMAG(Z))
     &              .GT.1.0D-16 ) THEN
                  Y = Z/A(IM1,J)
                  DO K = JP1,N
                     A(I,K) = A(I,K) - Y*A(IM1,K)
                  END DO
                  DO K = IM1,N
                     B(I,K) = B(I,K) - Y*B(IM1,K)
                  END DO
               END IF
               W = B(I,IM1)
               Z = B(I,I)
               IF ( DABS(DREAL(W))+DABS(DIMAG(W)).GT.DABS(DREAL(Z))
     &              +DABS(DIMAG(Z)) ) THEN
                  DO K = 1,I
                     Y = B(K,I)
                     B(K,I) = B(K,IM1)
                     B(K,IM1) = Y
                  END DO
                  DO K = 1,N
                     Y = A(K,I)
                     A(K,I) = A(K,IM1)
                     A(K,IM1) = Y
                  END DO
                  DO K = IMJ,N
                     Y = C(K,I)
                     C(K,I) = C(K,IM1)
                     C(K,IM1) = Y
                  END DO
               END IF
               Z = B(I,IM1)
C               IF ( DREAL(Z).NE.0.0D0 .OR. DIMAG(Z).NE.0.0D0 ) THEN
               IF ( ABS(DREAL(Z)).GT.1.0D-16 .OR. ABS(DIMAG(Z))
     &              .GT.1.0D-16 ) THEN
                  Y = Z/B(I,I)
                  DO K = 1,IM1
                     B(K,IM1) = B(K,IM1) - Y*B(K,I)
                  END DO
                  B(I,IM1) = (0.0D0,0.0D0)
                  DO K = 1,N
                     A(K,IM1) = A(K,IM1) - Y*A(K,I)
                  END DO
                  DO K = IMJ,N
                     C(K,IM1) = C(K,IM1) - Y*C(K,I)
                  END DO
               END IF
            END DO
            A(JP1+1,J) = (0.0D0,0.0D0)
         END DO
      END IF
      END
