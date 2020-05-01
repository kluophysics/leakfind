C*==cinvrij1.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVRIJ1(KIJ,A,C,MQ,TAUQ,IPIV,NQ,NLMQ,NN,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  inversion of the NON-relativistic KKR-matrix by factorisation   *
C   *  using the (l,m,s)-representation                                *
C   *                                                                  *
C   *                  (P | Q)           P = m(dn,dn) - G              *
C   *             M =  (-----)   with    S = 0                         *
C   *                  (R | S)           Q = 0                         *
C   *                                    R = 0                         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   * this version is used for:                                        *
C   * R   Real space (cluster) calculations                            *
C   * IJ  supplying the IJ (site-off-diagonal) blocks                  *
C   *     of TAU - matrix if requested (KIJ=1)                         *
C   *     on exit:     WA = TAUIJ(dn,dn)                               *
C   *                  WB = TAUIJ(up,up)                               *
C   * 1   implies IREL=0,1 i.e. NON-relativistic AND                   *
C   *                           NON-spin-polarized case                *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CINVRIJ126
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CINVRII3')
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER KIJ,NKMMAX,NN,NQ
      COMPLEX*16 A(NN,NN),C(NN,NN),MQ(NKMMAX,NKMMAX,NQ),
     &           TAUQ(NKMMAX,NKMMAX,NQ)
      INTEGER IPIV(NN),NLMQ(NQ)
C
C Local variables
C
      INTEGER I,II,II0,II1,INFO,IQ,J,JJ,N,NNSQ
C
C*** End of declarations rewritten by SPAG
C
      IF ( KIJ.LT.0 ) CALL STOP_MESSAGE(ROUTINE,'KIJ < 0 ')
      NNSQ = NN*NN
C
C------------------------------------------------- KKR-matrix: M = m - G
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         II1 = II0 + 1
         DO J = 1,N
            JJ = II0 + J
C
            CALL ZAXPY(N,C1,MQ(1,J,IQ),1,A(II1,JJ),1)
C
         END DO
         II0 = II0 + N
      END DO
C
C-------------------------------------------------------------- invert M
C
      CALL ZGETRF(NN,NN,A,NN,IPIV,INFO)
      CALL ZGETRI(NN,A,NN,IPIV,C,NNSQ,INFO)
C
C-------------------------------------------- store site-diagonal blocks
C
      II = II0
C
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         JJ = II0
         DO J = 1,N
            JJ = JJ + 1
            II = II0
            DO I = 1,N
               II = II + 1
               TAUQ(I,I,IQ) = A(II,JJ)
            END DO
         END DO
         II0 = II0 + N
      END DO
C
      END
C*==cinvrij2.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVRIJ2(KIJ,A,B,C,MQ,TAUQ,IPIV,NQ,NLM,NLMQ,NN,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  inversion of the NON-relativistic KKR-matrix by factorisation   *
C   *  using the (l,m,s)-representation                                *
C   *                                                                  *
C   *                  (P | Q)           P = m(dn,dn) - G              *
C   *             M =  (-----)   with    S = m(up,up) - G              *
C   *                  (R | S)           Q = 0                         *
C   *                                    R = 0                         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   * this version is used for:                                        *
C   * R   Real space (cluster) calculations                            *
C   * IJ  supplying the IJ (site-off-diagonal) blocks                  *
C   *     of TAU - matrix if requested (KIJ=1)                         *
C   *     on exit:     WA = TAUIJ(dn,dn)                               *
C   *                  WB = TAUIJ(up,up)                               *
C   * 2   implies IREL=2 i.e. NON-relativistic spin-polarized case     *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CINVRIJ2130
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CINVRIJ2')
      COMPLEX*16 C1
      PARAMETER (C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER KIJ,NKMMAX,NLM,NN,NQ
      COMPLEX*16 A(NN,NN),B(NN,NN),C(NN,NN),MQ(NKMMAX,NKMMAX,NQ),
     &           TAUQ(NKMMAX,NKMMAX,NQ)
      INTEGER IPIV(NN),NLMQ(NQ)
C
C Local variables
C
      INTEGER IDN,II,II0,II1,INFO,IQ,IUP,IUP1,JDN,JJ,JUP,N,NNSQ
C
C*** End of declarations rewritten by SPAG
C
      IF ( KIJ.LT.0 ) CALL STOP_MESSAGE(ROUTINE,'KIJ < 0 ')
      NNSQ = NN*NN
C
C---------------------------------------------------- initialize P and S
C
      CALL ZCOPY(NNSQ,A,1,B,1)
C
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         IUP1 = NLM + 1
         II1 = II0 + 1
         DO JDN = 1,N
            JUP = NLM + JDN
            JJ = II0 + JDN
C                                                               P: dn-dn
            CALL ZAXPY(N,C1,MQ(1,JDN,IQ),1,A(II1,JJ),1)
C                                                               S: up-up
            CALL ZAXPY(N,C1,MQ(IUP1,JUP,IQ),1,B(II1,JJ),1)
         END DO
         II0 = II0 + N
      END DO
C
C-------------------------------------------------------------- invert P
C
      CALL ZGETRF(NN,NN,A,NN,IPIV,INFO)
      CALL ZGETRI(NN,A,NN,IPIV,C,NNSQ,INFO)
C
C-------------------------------------------------------------- invert S
C
      CALL ZGETRF(NN,NN,B,NN,IPIV,INFO)
      CALL ZGETRI(NN,B,NN,IPIV,C,NNSQ,INFO)
C
C-------------------------------------------- store site-diagonal blocks
C
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         JJ = II0
         DO JDN = 1,N
            JUP = NLM + JDN
            JJ = JJ + 1
            II = II0
            DO IDN = 1,N
               IUP = NLM + IDN
               II = II + 1
C
               TAUQ(IDN,JDN,IQ) = A(II,JJ)
               TAUQ(IUP,JUP,IQ) = B(II,JJ)
            END DO
         END DO
         II0 = II0 + N
      END DO
C
      END
C*==cinvrii3.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVRII3(A,B,C,D,MQ,TAUQ,IPIV,IQ1,IQ2,NQ,NLM,NLMQ,NN,
     &                    NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  inversion of the relativistic KKR-matrix by factorisation       *
C   *  using the (l,m,s)-representation                                *
C   *                                                                  *
C   *                  (P | Q)           P = m(dn,dn) - G              *
C   *             M =  (-----)   with    S = m(up,up) - G              *
C   *                  (R | S)           Q = m(dn,up)                  *
C   *                                    R = m(up,dn)                  *
C   *                                                                  *
C   *   use is made of the spareness of matrices  Q  and  R            *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   * this version is used for:                                        *
C   * R   Real space (cluster) calculations                            *
C   * II  supplying the II (site-diagonal) blocks of the TAU - matrix  *
C   * 3   for  IQ = IQ1, IQ2  and the fully relativistic case (IREL=3) *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--CINVRII3250
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ1,IQ2,NKMMAX,NLM,NN,NQ
      COMPLEX*16 A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),
     &           MQ(NKMMAX,NKMMAX,NQ),TAUQ(NKMMAX,NKMMAX,NQ)
      INTEGER IPIV(NN),NLMQ(NQ)
C
C Local variables
C
      INTEGER IDN,II,II0,II1,INFO,IQ,IUP,IUP1,JDN,JJ,JJ0,JUP,KK,KK1,KK2,
     &        KUP,KUP1,KUP2,N,NNSQ
C
C*** End of declarations rewritten by SPAG
C
      NNSQ = NN*NN
C
C---------------------------------------------------- initialize P and S
C
      CALL ZCOPY(NNSQ,A,1,B,1)
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         IUP1 = NLM + 1
         II1 = II0 + 1
         DO JDN = 1,N
            JUP = NLM + JDN
            JJ = II0 + JDN
C                                                               P: dn-dn
            CALL ZAXPY(N,C1,MQ(1,JDN,IQ),1,A(II1,JJ),1)
C                                                               S: up-up
            CALL ZAXPY(N,C1,MQ(IUP1,JUP,IQ),1,B(II1,JJ),1)
         END DO
         II0 = II0 + N
      END DO
C
C-------------------------------------------------------------- invert S
C
      CALL ZGETRF(NN,NN,B,NN,IPIV,INFO)
      CALL ZGETRI(NN,B,NN,IPIV,C,NNSQ,INFO)
C
C---------------------------------------------------- store (1/S)  -> S~
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         IUP1 = NLM + 1
         II1 = II0 + 1
         JJ = II0
         DO JUP = IUP1,NLM + N
            JJ = JJ + 1
C                                                               S: up-up
            CALL ZCOPY(N,B(II1,JJ),1,TAUQ(IUP1,JUP,IQ),1)
         END DO
         II0 = II0 + N
      END DO
C
C--------------------------------------------------------- U = Q * (1/S)
C
      II = 0
      DO IQ = 1,NQ
         II0 = II
         N = NLMQ(IQ)
         KK1 = II0 + 1
         KK2 = II0 + N
         DO IDN = 1,N
            II = II + 1
            DO JJ = 1,NN
               C(II,JJ) = C0
               KUP = NLM
               DO KK = KK1,KK2
                  KUP = KUP + 1
                  C(II,JJ) = C(II,JJ) + MQ(IDN,KUP,IQ)*B(KK,JJ)
               END DO
            END DO
         END DO
      END DO
C
C----------------------------------------------------------- W = P - U*R
C----------------------------------------------------------- V = (1/S)*R
C
      JJ = 0
      DO IQ = 1,NQ
         JJ0 = JJ
         N = NLMQ(IQ)
         KK1 = JJ0 + 1
         KK2 = JJ0 + N
         KUP1 = NLM + 1
         KUP2 = NLM + N
         DO JDN = 1,N
            JJ = JJ + 1
            DO II = 1,NN
               D(II,JJ) = C0
               KK = JJ0
               DO KUP = KUP1,KUP2
                  KK = KK + 1
                  A(II,JJ) = A(II,JJ) - C(II,KK)*MQ(KUP,JDN,IQ)
                  D(II,JJ) = D(II,JJ) + B(II,KK)*MQ(KUP,JDN,IQ)
               END DO
            END DO
         END DO
      END DO
C
C------------------------------------------------------- invert W  -> P~
C
      CALL ZGETRF(NN,NN,A,NN,IPIV,INFO)
      CALL ZGETRI(NN,A,NN,IPIV,B,NNSQ,INFO)
C
C----------------------------------------------------------- Q_ = P~ * U
C
      CALL ZGEMM('N','N',NN,NN,NN,C1,A,NN,C,NN,C0,B,NN)
C
C-------------------------------------------------------------- store P~
C-------------------------------------------------------------- store Q~
C--------------------------------------------------------- R~ = - V * P~
C--------------------------------------------------- S~ = V * Q_ + (1/S)
C
      II0 = 0
      IF ( IQ1.NE.1 ) THEN
         DO IQ = IQ1,IQ1 - 1
            II0 = II0 + NLMQ(IQ)
         END DO
      END IF
C
      II = II0
C
      DO IQ = IQ1,IQ2
         N = NLMQ(IQ)
         DO IDN = 1,N
            IUP = NLM + IDN
            II = II + 1
            JJ = II0
            DO JDN = 1,N
               JUP = NLM + JDN
               JJ = JJ + 1
               TAUQ(IDN,JDN,IQ) = A(II,JJ)
               TAUQ(IDN,JUP,IQ) = -B(II,JJ)
               TAUQ(IUP,JDN,IQ) = C0
               DO KK = 1,NN
                  TAUQ(IUP,JDN,IQ) = TAUQ(IUP,JDN,IQ) - D(II,KK)
     &                               *A(KK,JJ)
                  TAUQ(IUP,JUP,IQ) = TAUQ(IUP,JUP,IQ) + D(II,KK)
     &                               *B(KK,JJ)
               END DO
            END DO
         END DO
         II0 = II0 + N
      END DO
C
      END
C*==cinvrij3.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CINVRIJ3(A,B,C,D,MQ,TAUQ,IPIV,NQ,NLM,NLMQ,NN,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *  inversion of the relativistic KKR-matrix by factorisation       *
C   *  using the (l,m,s)-representation                                *
C   *                                                                  *
C   *                  (P | Q)           P = m(dn,dn) - G              *
C   *             M =  (-----)   with    S = m(up,up) - G              *
C   *                  (R | S)           Q = m(dn,up)                  *
C   *                                    R = m(up,dn)                  *
C   *                                                                  *
C   *   use is made of the spareness of matrices  Q  and  R            *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   * this version is used for:                                        *
C   * R   Real space (cluster) calculations                            *
C   * II  supplying the II (site-diagonal) blocks of the TAU - matrix  *
C   * 3   for  IQ = IQ1, IQ2  and the fully relativistic case (IREL=3) *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--CINVRIJ3439
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CINVRIJ3')
C
C Dummy arguments
C
      INTEGER NKMMAX,NLM,NN,NQ
      COMPLEX*16 A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),
     &           MQ(NKMMAX,NKMMAX,NQ),TAUQ(NKMMAX,NKMMAX,NQ)
      INTEGER IPIV(NN),NLMQ(NQ)
C
C Local variables
C
      COMPLEX*16 AKJ,BKJ,E(:,:)
      INTEGER I,IA_ERR,IDN,II,II0,II1,INFO,IQ,IUP,IUP1,J,JDN,JJ,JJ0,JUP,
     &        K,KK,KK1,KK2,KUP,KUP1,KUP2,N,NNSQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE E
C
      ALLOCATE (E(NN,NN),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: E')
C
      NNSQ = NN*NN
C
C---------------------------------------------------- initialize P and S
C
      CALL ZCOPY(NNSQ,A,1,D,1)
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         IUP1 = NLM + 1
         II1 = II0 + 1
         DO JDN = 1,N
            JUP = NLM + JDN
            JJ = II0 + JDN
C                                                               P: dn-dn
            CALL ZAXPY(N,C1,MQ(1,JDN,IQ),1,A(II1,JJ),1)
C                                                               S: up-up
            CALL ZAXPY(N,C1,MQ(IUP1,JUP,IQ),1,D(II1,JJ),1)
         END DO
         II0 = II0 + N
      END DO
C
C-------------------------------------------------------------- invert S
C
      CALL ZGETRF(NN,NN,D,NN,IPIV,INFO)
      CALL ZGETRI(NN,D,NN,IPIV,C,NNSQ,INFO)
C
C--------------------------------------------------------- U = Q * (1/S)
C
      II = 0
      DO IQ = 1,NQ
         II0 = II
         N = NLMQ(IQ)
         KK1 = II0 + 1
         KK2 = II0 + N
         DO IDN = 1,N
            II = II + 1
            DO JJ = 1,NN
               C(II,JJ) = C0
               KUP = NLM
               DO KK = KK1,KK2
                  KUP = KUP + 1
                  C(II,JJ) = C(II,JJ) + MQ(IDN,KUP,IQ)*D(KK,JJ)
               END DO
            END DO
         END DO
      END DO
C
C----------------------------------------------------------- W = P - U*R
C----------------------------------------------------------- V = (1/S)*R
C
      JJ = 0
      DO IQ = 1,NQ
         JJ0 = JJ
         N = NLMQ(IQ)
         KK1 = JJ0 + 1
         KK2 = JJ0 + N
         KUP1 = NLM + 1
         KUP2 = NLM + N
         DO JDN = 1,N
            JJ = JJ + 1
            DO II = 1,NN
               E(II,JJ) = C0
               KK = JJ0
               DO KUP = KUP1,KUP2
                  KK = KK + 1
                  A(II,JJ) = A(II,JJ) - C(II,KK)*MQ(KUP,JDN,IQ)
                  E(II,JJ) = E(II,JJ) + D(II,KK)*MQ(KUP,JDN,IQ)
               END DO
            END DO
         END DO
      END DO
C
C------------------------------------------------------- invert W  -> P~
C
      CALL ZGETRF(NN,NN,A,NN,IPIV,INFO)
      CALL ZGETRI(NN,A,NN,IPIV,B,NNSQ,INFO)
C
C----------------------------------------------------------- Q~ = P~ * U
C
      CALL ZGEMM('N','N',NN,NN,NN,-C1,A,NN,C,NN,C0,B,NN)
C
C--------------------------------------------------------- R~ = - V * P~
C------------------------------------------------- S~ = - V * Q~ + (1/S)
C
      DO J = 1,NN
         DO I = 1,NN
            C(I,J) = C0
         END DO
      END DO
C
      DO J = 1,NN
         DO K = 1,NN
            AKJ = A(K,J)
            BKJ = B(K,J)
            DO I = 1,NN
               C(I,J) = C(I,J) - E(I,K)*AKJ
               D(I,J) = D(I,J) - E(I,K)*BKJ
            END DO
         END DO
      END DO
C
C-------------------------- store site-diagonal blocks of P~, Q~, R~, S~
C
      II0 = 0
      DO IQ = 1,NQ
         N = NLMQ(IQ)
         JJ = II0
         DO JDN = 1,N
            JUP = NLM + JDN
            JJ = JJ + 1
            II = II0
            DO IDN = 1,N
               IUP = NLM + IDN
               II = II + 1
C                                                                  dn-dn
               TAUQ(IDN,JDN,IQ) = A(II,JJ)
C                                                                  dn-up
               TAUQ(IDN,JUP,IQ) = B(II,JJ)
C                                                                  up-dn
               TAUQ(IUP,JDN,IQ) = C(II,JJ)
C                                                                  up-up
               TAUQ(IUP,JUP,IQ) = D(II,JJ)
C
            END DO
         END DO
         II0 = II0 + N
      END DO
C
      DEALLOCATE (E)
C
      END
