C*==taugfconv.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE TAUGFCONV(MSS,TAU,TMAT)
C   ********************************************************************
C   *                                                                  *
C   *  set the site-projected multiple scattering matrix  TMAT         *
C   *  that includes the back scattering according to convention used  *
C   *                                                                  *
C   *  GF_CONV_RH = .T.      TMAT = G[i,i] = m[i] TAU[i,i] m[i] - m[i] *
C   *                                                                  *
C   *  GF_CONV_RH = .F.      TMAT =               TAU[i,i]             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_CALCMODE,ONLY:GF_CONV_RH
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--TAUGFCONV18
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 MSS(NKMMAX,NKMMAX),TAU(NKMMAX,NKMMAX),
     &           TMAT(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER N
      COMPLEX*16 W1(:,:),W2(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1,W2
C
      ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX))
C
      TMAT(1:NKMMAX,1:NKMMAX) = C0
C
C----------------------------------------- convention for Green function
      IF ( GF_CONV_RH ) THEN
C                                                basis functions R and H
C
         N = NKM
         W1(1:N,1:N) = MATMUL(MSS(1:N,1:N),TAU(1:N,1:N))
         W2(1:N,1:N) = MATMUL(W1(1:N,1:N),MSS(1:N,1:N))
C
         TMAT(1:N,1:N) = W2(1:N,1:N) - MSS(1:N,1:N)
C
      ELSE
C                                                basis functions Z and J
C
         TMAT(1:NKMMAX,1:NKMMAX) = TAU(1:NKMMAX,1:NKMMAX)
C
      END IF
C
      END
