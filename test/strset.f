C*==strset.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRSET(IK,KVEC,GNR,TAUKINV,P)
C   ********************************************************************
C   *                                                                  *
C   *   SETS UP NON-RELATIVISTIC STRUCTURE CONSTANTS  GNR              *
C   *   FOR  REAL SPHERICAL HARMONICS AND CONVERTS TO  RELAT.  G'S     *
C   *                                                                  *
C   *   NOTE:  - GNR  IS STORED IN  SUPPLIED WORKSPACE                 *
C   *          - ONLY THE ELEMENTS  G(I,J) FOR I<=J ARE EVALUATED      *
C   *            THOSE FOR J>I ARE OBTAINED FROM THE FACT THAT         *
C   *            G(I,J) IS HERMITIAN:        G(I,J) = G(J,I)*          *
C   *          - FOR IREL<2:                 G(I,J) = GNR(I,J)         *
C   *          - G(I,J) IS STORED AS -G IN TAUKINV  !!!!!              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_ANGMOM,ONLY:NL,NLM,NKMQ,NLMQ,IND0Q,NKKR
      USE MOD_STR,ONLY:STRMODE,NQQP_STR,RGNT,NRGNT,IRGNT,SRREL,NRREL,
     &    IRREL,NIJQ,IJQ,DLLMMKE,WK,CILMAT
      USE MOD_CONSTANTS,ONLY:CI,C0
c modified by XJQ: allow a supercell containing > 100 atoms
      use mod_sites,only:nqhost
c end-mod-xjq
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IK
      COMPLEX*16 P
      COMPLEX*16 GNR(NKKR,NKKR),TAUKINV(NKKR,NKKR)
      REAL*8 KVEC(3)
C
C Local variables
C
      COMPLEX*16 CIP,CSUM,CSUM1,CSUM2
      INTEGER I,I1,I2,IG123,IKM1,IKM2,ILM1,ILM2,IQ1,IQ2,IQQP,IS,J1,J2,
     &        JKM1,JKM10,JKM2,JLM1,JLM10,JLM2,JQ1,JQ2,JQQP,LM,LM1,
     &        LM1LM2,LM2,LM3,N1,N2
      REAL*8 KX,KY,KZ
C
C*** End of declarations rewritten by SPAG
C
      KX = KVEC(1)
      KY = KVEC(2)
      KZ = KVEC(3)
C  =====================================================================
C                calculate auxilary function  DLM(K,E)
C  =====================================================================
C
      IF ( STRMODE.EQ.1 ) THEN
C-------------------------------------------------------------- STANDARD
C
         CALL STRBBDD(DLLMMKE,KX,KY,KZ)
C
      ELSE IF ( STRMODE.EQ.2 ) THEN
C--------------------------------------------------------- INTERPOLATION
C
         CALL STRBBDDLAG(DLLMMKE,IK,KX,KY,KZ)
C
      ELSE
C------------------------------------ for test purposes SPLITTING SCHEME
C
         CALL STRBBDDSPLIT(DLLMMKE,KX,KY,KZ)
C
      END IF

c     do i=1, NQQP_STR
c       write(*, "(i, *('('sf6.2xspf6.2x'i)':x))") i,DLLMMKE(1,i)
c     end do
C
C  =====================================================================
C               set up   KKR strcuture constants matrix  G(k,E)
C  =====================================================================
C
      CIP = CI*P
C
      IF ( IREL.LE.2 ) THEN
C
C     NON-RELATIVISTIC CALCULATION
C     ============================
C
         DO IQQP = 1,NQQP_STR
c modified by XJQ: allow a supercell containing > 99 atoms
            IQ1 = IJQ(1,IQQP)/max(nqhost+1,100)
            IQ2 = IJQ(1,IQQP) - max(nqhost+1,100)*IQ1
c end-mod-xjq
C
            IG123 = 0
            LM1LM2 = 0
            DO LM1 = 1,NL**2
               DO LM2 = 1,LM1
                  LM1LM2 = LM1LM2 + 1
                  CSUM = C0
                  DO I = 1,NRGNT(LM1LM2)
                     IG123 = IG123 + 1
                     LM3 = IRGNT(IG123)
                     CSUM = CSUM + RGNT(IG123)*DLLMMKE(LM3,IQQP)
                  END DO
                  WK(LM1,LM2) = -CSUM*CILMAT(LM1,LM2)
                  WK(LM2,LM1) = -CSUM*CILMAT(LM2,LM1)
               END DO
            END DO
C
            IF ( IQQP.EQ.1 ) THEN
               DO LM = 1,NLM
                  WK(LM,LM) = WK(LM,LM) - CIP
               END DO
            END IF
C
C     ---------------------------------------------------
C                          COPY
C     ---------------------------------------------------
C
            DO JQQP = 1,NIJQ(IQQP)
c modified by XJQ: allow a supercell containing > 99 atoms
               JQ1 = IJQ(JQQP,IQQP)/max(nqhost+1,100)
               JQ2 = IJQ(JQQP,IQQP) - max(nqhost+1,100)*JQ1
c end-mod-xjq
               JLM2 = IND0Q(JQ2)
               JLM10 = IND0Q(JQ1)
               DO ILM2 = 1,NLMQ(JQ2)
                  JLM2 = JLM2 + 1
                  JLM1 = JLM10
                  DO ILM1 = 1,NLMQ(JQ1)
                     JLM1 = JLM1 + 1
                     TAUKINV(JLM1,JLM2) = WK(ILM1,ILM2)
                  END DO
               END DO
            END DO
C
         END DO
         RETURN
C
      ELSE
C
C     RELATIVISTIC CALCULATION
C     ========================
C
C     ---------------------------------------------------
C     set up non-relativistic str.const.  GNR  for every
C     representative q-qp-block  IQQP, then transfrom to
C     relativ. representation and finally copy result
C     to the  NIJQ(IQQP)-1  equivalent q-qp-blocks  JQQP
C     ---------------------------------------------------
C
         DO IQQP = 1,NQQP_STR
c modified by XJQ: allow a supercell containing > 99 atoms
            IQ1 = IJQ(1,IQQP)/max(nqhost+1,100)
            IQ2 = IJQ(1,IQQP) - max(nqhost+1,100)*IQ1
c end-mod-xjq
C
            IG123 = 0
            LM1LM2 = 0
            DO LM1 = 1,NL**2
               DO LM2 = 1,LM1
                  LM1LM2 = LM1LM2 + 1
                  CSUM = 0.0D0
C
                  DO I = 1,NRGNT(LM1LM2)
                     IG123 = IG123 + 1
                     LM3 = IRGNT(IG123)
                     CSUM = CSUM + RGNT(IG123)*DLLMMKE(LM3,IQQP)
                  END DO
C
                  GNR(LM1,LM2) = CSUM*CILMAT(LM1,LM2)
                  GNR(LM2,LM1) = CSUM*CILMAT(LM2,LM1)
C
               END DO
            END DO
C
C                                 SITE-DIAGONAL BLOCK
            IF ( IQQP.EQ.1 ) THEN
               DO LM = 1,NLM
                  GNR(LM,LM) = GNR(LM,LM) + CIP
               END DO
            END IF
C
C     ---------------------------------------------------
C                          TRANSFORM
C     ---------------------------------------------------
C
            DO IKM2 = 1,NKMQ(IQ2)
               DO IKM1 = 1,NKMQ(IQ1)
C
                  CSUM1 = C0
                  DO IS = 1,2
                     N1 = NRREL(IS,IKM1)
                     N2 = NRREL(IS,IKM2)
                     DO I1 = 1,N1
                        J1 = IRREL(I1,IS,IKM1)
C
                        CSUM2 = C0
                        DO I2 = 1,N2
                           J2 = IRREL(I2,IS,IKM2)
                           CSUM2 = CSUM2 + GNR(J1,J2)*SRREL(I2,IS,IKM2)
                        END DO
C
                        CSUM1 = CSUM1 + DCONJG(SRREL(I1,IS,IKM1))*CSUM2
                     END DO
                  END DO
                  WK(IKM1,IKM2) = -CSUM1
               END DO
            END DO
C
C     ---------------------------------------------------
C                          COPY
C     ---------------------------------------------------
C
            DO JQQP = 1,NIJQ(IQQP)
c modified by XJQ: allow a supercell containing > 99 atoms
               JQ1 = IJQ(JQQP,IQQP)/max(nqhost+1,100)
               JQ2 = IJQ(JQQP,IQQP) - max(nqhost+1,100)*JQ1
c end-mod-xjq
               JKM2 = IND0Q(JQ2)
               JKM10 = IND0Q(JQ1)
               DO IKM2 = 1,NKMQ(JQ2)
                  JKM2 = JKM2 + 1
                  JKM1 = JKM10
                  DO IKM1 = 1,NKMQ(JQ1)
                     JKM1 = JKM1 + 1
                     TAUKINV(JKM1,JKM2) = WK(IKM1,IKM2)
                  END DO
               END DO
C
            END DO
C
         END DO
C
      END IF
C
C=======================================================================
      END
C*==strbbdd.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRBBDD(DLLMMKE,KX,KY,KZ)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the missing   ->k - dependent quantities              *
C   *  then set up   DLM(K,E)   by performing the lattice sum          *
C   *  the DLM's for the various (IQ,IQ')-blocks of  G  are            *
C   *  evaluated in parallel - as far as possible                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ABAS,BBAS
      USE MOD_STR,ONLY:NLLMMMAX,LLMAX,MMLLMAX,NQQP_STR,NGRL,GMAXSQ,ETA,
     &    HP,G123MAX,R123MAX,QQPX,QQPY,QQPZ,SMAX,R1,R2,R3,G1,G2,G3,
     &    EXPGNQ,INDR,QQMLRS,D1TERM3,PWEX2K1,PWEX2K2,PWEX2K3,PWEXIK1,
     &    PWEXIK2,PWEXIK3,D300,EDU,NQQP_STR_RED,M1PWL,
     &    USE_NEW_BBDD_VERSION
      USE MOD_CONSTANTS,ONLY:CI2PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQP_STR)
C
C Local variables
C
      COMPLEX*16 DENOM,EXIK1,EXIK2,EXIK3,EXIKQPQ,EXIKQQP,EXIKRS,EXIKRSX,
     &           F1,F2,F2QPQ,F2QQP,F3,F3QPQ,F3QQP
      REAL*8 EX2K1,EX2K2,EX2K3,EX2KGN,F0,KNX,KNY,KNZ
      INTEGER I,IPW,IQPQ,IQQP,LL,MM,MMLL,N,S
C
C*** End of declarations rewritten by SPAG
C
C  =====================================================================
C                calculate auxilary function  DLM(K,E)
C  =====================================================================
C
      CALL CINIT(NLLMMMAX*NQQP_STR,DLLMMKE)
C
C  =====================================================================
C                               ********
C                               * DLM1 *
C                               ********
C
C     TABLE TO CALCULATE         EXP( -2 * ->KK[N] * ->K )
C
      EX2K1 = EXP(-2.0D0*(BBAS(1,1)*KX+BBAS(2,1)*KY+BBAS(3,1)*KZ)/ETA)
      EX2K2 = EXP(-2.0D0*(BBAS(1,2)*KX+BBAS(2,2)*KY+BBAS(3,2)*KZ)/ETA)
      EX2K3 = EXP(-2.0D0*(BBAS(1,3)*KX+BBAS(2,3)*KY+BBAS(3,3)*KZ)/ETA)
C
      PWEX2K1(0) = 1.0D0
      PWEX2K2(0) = 1.0D0
      PWEX2K3(0) = 1.0D0
C
      DO IPW = 1,G123MAX
         PWEX2K1(IPW) = PWEX2K1(IPW-1)*EX2K1
         PWEX2K1(-IPW) = 1.0D0/PWEX2K1(IPW)
C
         PWEX2K2(IPW) = PWEX2K2(IPW-1)*EX2K2
         PWEX2K2(-IPW) = 1.0D0/PWEX2K2(IPW)
C
         PWEX2K3(IPW) = PWEX2K3(IPW-1)*EX2K3
         PWEX2K3(-IPW) = 1.0D0/PWEX2K3(IPW)
      END DO
C
C  ---------------------------------------------------------------------
C                                                 RECIPROCAL LATTICE SUM
      F0 = EXP(-(KX*KX+KY*KY+KZ*KZ)/ETA)
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( USE_NEW_BBDD_VERSION ) THEN
C
         DO N = 1,NGRL
C
C     ->K(N) = ->K + ->KK(N)
C
            KNX = KX + G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)
     &            *BBAS(1,3)
            KNY = KY + G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)
     &            *BBAS(2,3)
            KNZ = KZ + G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)
     &            *BBAS(3,3)
C
            DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
C
C
            IF ( DREAL(DENOM).LE.GMAXSQ ) THEN
C  .....................................................................
C
               EX2KGN = PWEX2K1(G1(N))*PWEX2K2(G2(N))*PWEX2K3(G3(N))
C
               F1 = F0*EX2KGN/DENOM
C
C-----------------------------------------------------------------------
C
               CALL CALC_RHPLM(KNX,KNY,KNZ,HP,LLMAX,NLLMMMAX)
C
C-----------------------------------------------------------------------
C
               F2 = F1*EXPGNQ(N,1)
C
               MMLL = 0
               DO LL = 0,LLMAX
                  F3 = F2*D1TERM3(LL)
C
                  DO MM = -LL, + LL
                     MMLL = MMLL + 1
                     DLLMMKE(MMLL,1) = DLLMMKE(MMLL,1) + F3*HP(MMLL)
                  END DO
C
               END DO
C
               DO IQQP = 2,NQQP_STR_RED
                  IQPQ = IQQP + NQQP_STR_RED - 1
C
                  F2QQP = F1*EXPGNQ(N,IQQP)
C                 F2QPQ = F1*EXPGNQ(N,IQPQ)
                  F2QPQ = F1*DCONJG(EXPGNQ(N,IQQP))
C
                  MMLL = 0
                  DO LL = 0,LLMAX
                     F3QQP = F2QQP*D1TERM3(LL)
                     F3QPQ = F2QPQ*D1TERM3(LL)
C
                     DO MM = -LL, + LL
                        MMLL = MMLL + 1
                        DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                     + F3QQP*HP(MMLL)
                        DLLMMKE(MMLL,IQPQ) = DLLMMKE(MMLL,IQPQ)
     &                     + F3QPQ*HP(MMLL)
                     END DO
C
                  END DO
C
               END DO
C
            END IF
C  .....................................................................
C
         END DO
C
C  MULTIPLY THE MISSING PHASE FACTOR    EXP( I ->K * (->Q[I] - ->Q[J]) )
C
         DO IQQP = 2,NQQP_STR_RED
            IQPQ = IQQP + NQQP_STR_RED - 1
            EXIKQQP = CDEXP
     &                (CI2PI*(KX*QQPX(IQQP)+KY*QQPY(IQQP)+KZ*QQPZ(IQQP))
     &                )
            EXIKQPQ = DCONJG(EXIKQQP)
C
C           DO MMLL = 1,MMLLMAX
C              DLLMMKE(MMLL,IQQP) = EXIKQQP*DLLMMKE(MMLL,IQQP)
C              DLLMMKE(MMLL,IQPQ) = EXIKQPQ*DLLMMKE(MMLL,IQPQ)
C           END DO
            CALL ZSCAL(MMLLMAX,EXIKQQP,DLLMMKE(1,IQQP),1)
            CALL ZSCAL(MMLLMAX,EXIKQPQ,DLLMMKE(1,IQPQ),1)
C
         END DO
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
C
         DO N = 1,NGRL
C
C     ->K(N) = ->K + ->KK(N)
C
            KNX = KX + G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)
     &            *BBAS(1,3)
            KNY = KY + G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)
     &            *BBAS(2,3)
            KNZ = KZ + G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)
     &            *BBAS(3,3)
C
            DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
C
C
            IF ( DREAL(DENOM).LE.GMAXSQ ) THEN
C  .....................................................................
C
               EX2KGN = PWEX2K1(G1(N))*PWEX2K2(G2(N))*PWEX2K3(G3(N))
C
               F1 = F0*EX2KGN/DENOM
C
C-----------------------------------------------------------------------
C
               CALL CALC_RHPLM(KNX,KNY,KNZ,HP,LLMAX,NLLMMMAX)
C
C-----------------------------------------------------------------------
C
               DO IQQP = 1,NQQP_STR
C
                  F2 = F1*EXPGNQ(N,IQQP)
C*
C*                MMLL = 1
C*                DO LL = 0,LLMAX
C*                N2LLP1 = LL + LL + 1
C*                CALL ZAXPY(N2LLP1,F3,HP(MMLL),1,DLLMMKE(MMLL,IQQP),1)
C*                MMLL = MMLL + N2LLP1
C
                  MMLL = 0
                  DO LL = 0,LLMAX
                     F3 = F2*D1TERM3(LL)
C
                     DO MM = -LL, + LL
                        MMLL = MMLL + 1
                        DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                     + F3*HP(MMLL)
                     END DO
C
                  END DO
C
               END DO
            END IF
C  .....................................................................
C
         END DO
C
C  MULTIPLY THE MISSING PHASE FACTOR    EXP( I ->K * (->Q[I] - ->Q[J]) )
C
         DO IQQP = 2,NQQP_STR
            EXIKQQP = CDEXP
     &                (CI2PI*(KX*QQPX(IQQP)+KY*QQPY(IQQP)+KZ*QQPZ(IQQP))
     &                )
C
            DO MMLL = 1,MMLLMAX
               DLLMMKE(MMLL,IQQP) = EXIKQQP*DLLMMKE(MMLL,IQQP)
            END DO
C
C           CALL ZSCAL(MMLLMAX,EXIKQQP,DLLMMKE(1,IQQP),1)
C
         END DO
C
      END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  ---------------------------------------------------------------------
C
C
C  =====================================================================
C                               ********
C                               * DLM2 *
C                               ********
C
C     TABLE TO CALCULATE         EXP( I ->K * ->R[S] )
C
      EXIK1 = CDEXP(CI2PI*(ABAS(1,1)*KX+ABAS(2,1)*KY+ABAS(3,1)*KZ))
      EXIK2 = CDEXP(CI2PI*(ABAS(1,2)*KX+ABAS(2,2)*KY+ABAS(3,2)*KZ))
      EXIK3 = CDEXP(CI2PI*(ABAS(1,3)*KX+ABAS(2,3)*KY+ABAS(3,3)*KZ))
C
      PWEXIK1(0) = 1.0D0
      PWEXIK2(0) = 1.0D0
      PWEXIK3(0) = 1.0D0
C
      DO IPW = 1,R123MAX
         PWEXIK1(IPW) = PWEXIK1(IPW-1)*EXIK1
         PWEXIK1(-IPW) = 1.0D0/PWEXIK1(IPW)
C
         PWEXIK2(IPW) = PWEXIK2(IPW-1)*EXIK2
         PWEXIK2(-IPW) = 1.0D0/PWEXIK2(IPW)
C
         PWEXIK3(IPW) = PWEXIK3(IPW-1)*EXIK3
         PWEXIK3(-IPW) = 1.0D0/PWEXIK3(IPW)
      END DO
C
C  ---------------------------------------------------------------------
C                                                     DIRECT LATTICE SUM
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( USE_NEW_BBDD_VERSION ) THEN
C
         DO S = 1,SMAX(1)
C
            I = INDR(S,1)
            EXIKRS = PWEXIK1(R1(I))*PWEXIK2(R2(I))*PWEXIK3(R3(I))
C
            CALL ZAXPY(MMLLMAX,EXIKRS,QQMLRS(1,S,1),1,DLLMMKE(1,1),1)
C
         END DO
C
         DO IQQP = 2,NQQP_STR_RED
            IQPQ = IQQP + NQQP_STR_RED - 1
C
            DO S = 1,SMAX(IQQP)
C
               I = INDR(S,IQQP)
               EXIKRS = PWEXIK1(R1(I))*PWEXIK2(R2(I))*PWEXIK3(R3(I))
               EXIKRSX = DCONJG(EXIKRS)
C
               DO MMLL = 1,MMLLMAX
                  DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                                 + EXIKRS*QQMLRS(MMLL,S,IQQP)
                  DLLMMKE(MMLL,IQPQ) = DLLMMKE(MMLL,IQPQ)
     &                                 + EXIKRSX*QQMLRS(MMLL,S,IQQP)
     &                                 *M1PWL(MMLL)
               END DO
C
            END DO
C
         END DO
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
C
         DO IQQP = 1,NQQP_STR
            DO S = 1,SMAX(IQQP)
C
               I = INDR(S,IQQP)
               EXIKRS = PWEXIK1(R1(I))*PWEXIK2(R2(I))*PWEXIK3(R3(I))
C
               CALL ZAXPY(MMLLMAX,EXIKRS,QQMLRS(1,S,IQQP),1,
     &                    DLLMMKE(1,IQQP),1)
C
            END DO
         END DO
C
      END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  =====================================================================
C                               ********
C                               * DLM3 *
C                               ********
C
C             ADD THE MISSING TERM   D300
C
      DLLMMKE(1,1) = DLLMMKE(1,1) + D300
C
C
      END
C*==strbbddsplit.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRBBDDSPLIT(DLLMMKE,KX,KY,KZ)
C   ********************************************************************
C   *                                                                  *
C   *  for testing: call <STRBBDDRED> and <STRBBDDPOLE>                *
C   *  and add up results                                              *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_STR,ONLY:NLLMMMAX,NQQP_STR
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQP_STR)
C
C*** End of declarations rewritten by SPAG
C
      CALL STRBBDDRED(DLLMMKE,KX,KY,KZ)
C
      CALL STRBBDDPOLE(DLLMMKE,KX,KY,KZ)
C
      END
C*==strbbddred.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRBBDDRED(DLLMMKE,KX,KY,KZ)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the missing   ->k - dependent quantities              *
C   *  then set up   DLM(K,E)   by performing the lattice sum          *
C   *  the DLM's for the various (IQ,IQ')-blocks of  G  are            *
C   *  evaluated in parallel - as far as possible                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:NGFEP
      USE MOD_LATTICE,ONLY:ABAS,BBAS
      USE MOD_STRLAG,ONLY:EXPGNQ_RED
      USE MOD_STR,ONLY:NLLMMMAX,LLMAX,MMLLMAX,NQQP_STR,NGRL,GMAXSQ,ETA,
     &    HP,G123MAX,R123MAX,QQPX,QQPY,QQPZ,SMAX,R1,R2,R3,G1,G2,G3,
     &    EXPGNQ,INDR,QQMLRS,D1TERM3,PWEX2K1,PWEX2K2,PWEX2K3,PWEXIK1,
     &    PWEXIK2,PWEXIK3,D300,EDU,PWP
      USE MOD_CONSTANTS,ONLY:CI2PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQP_STR)
C
C Local variables
C
      COMPLEX*16 DENOM,EXIK1,EXIK2,EXIK3,EXIKQQP,EXIKRS,F1,F2,F3,H0,H1,
     &           H2
      REAL*8 EX2K1,EX2K2,EX2K3,EX2KGN,F0,KNX,KNY,KNZ
      INTEGER I,IPW,IQQP,LL,MM,MMLL,N,S
C
C*** End of declarations rewritten by SPAG
C
C  =====================================================================
C                calculate auxilary function  DLM(K,E)
C  =====================================================================
C
      CALL CINIT(NLLMMMAX*NQQP_STR,DLLMMKE)
C
C  =====================================================================
C                               ********
C                               * DLM1 *
C                               ********
C
C     TABLE TO CALCULATE         EXP( -2 * ->KK[N] * ->K )
C
      EX2K1 = EXP(-2.0D0*(BBAS(1,1)*KX+BBAS(2,1)*KY+BBAS(3,1)*KZ)/ETA)
      EX2K2 = EXP(-2.0D0*(BBAS(1,2)*KX+BBAS(2,2)*KY+BBAS(3,2)*KZ)/ETA)
      EX2K3 = EXP(-2.0D0*(BBAS(1,3)*KX+BBAS(2,3)*KY+BBAS(3,3)*KZ)/ETA)
C
      PWEX2K1(0) = 1.0D0
      PWEX2K2(0) = 1.0D0
      PWEX2K3(0) = 1.0D0
C
      DO IPW = 1,G123MAX
         PWEX2K1(IPW) = PWEX2K1(IPW-1)*EX2K1
         PWEX2K1(-IPW) = 1.0D0/PWEX2K1(IPW)
C
         PWEX2K2(IPW) = PWEX2K2(IPW-1)*EX2K2
         PWEX2K2(-IPW) = 1.0D0/PWEX2K2(IPW)
C
         PWEX2K3(IPW) = PWEX2K3(IPW-1)*EX2K3
         PWEX2K3(-IPW) = 1.0D0/PWEX2K3(IPW)
      END DO
C
C  ---------------------------------------------------------------------
C                                                 RECIPROCAL LATTICE SUM
      F0 = EXP(-(KX*KX+KY*KY+KZ*KZ)/ETA)
C
      DO N = NGFEP + 1,NGRL
C
C     ->K(N) = ->K + ->KK(N)
C
         KNX = KX + G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)*BBAS(1,3)
         KNY = KY + G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)*BBAS(2,3)
         KNZ = KZ + G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)*BBAS(3,3)
C
         DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
C
         IF ( DREAL(DENOM).LE.GMAXSQ ) THEN
C  .....................................................................
C
            EX2KGN = PWEX2K1(G1(N))*PWEX2K2(G2(N))*PWEX2K3(G3(N))
C
            F1 = F0*EX2KGN/DENOM
C
C-----------------------------------------------------------------------
            CALL CALC_RHPLM(KNX,KNY,KNZ,HP,LLMAX,NLLMMMAX)
C-----------------------------------------------------------------------
C
            DO IQQP = 1,NQQP_STR
C
               F2 = F1*EXPGNQ(N,IQQP)
C
               MMLL = 0
               DO LL = 0,LLMAX
                  F3 = F2*D1TERM3(LL)
C
                  DO MM = -LL, + LL
                     MMLL = MMLL + 1
                     DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                  + F3*HP(MMLL)
                  END DO
               END DO
C
            END DO
         END IF
C  .....................................................................
C
      END DO
C
C-----------------------------------------------------------------------
C                   loop over possible free electron poles
C-----------------------------------------------------------------------
C
      H0 = EXP(-EDU/ETA)
C
      DO N = 1,NGFEP
C
C     ->K(N) = ->K + ->KK(N)
C
         KNX = KX + G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)*BBAS(1,3)
         KNY = KY + G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)*BBAS(2,3)
         KNZ = KZ + G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)*BBAS(3,3)
C
         DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
C
         EX2KGN = PWEX2K1(G1(N))*PWEX2K2(G2(N))*PWEX2K3(G3(N))
C
         F1 = F0*EX2KGN/DENOM
         H1 = H0/DENOM
C
         CALL CALC_RHPLM(KNX,KNY,KNZ,HP,LLMAX,NLLMMMAX)
C
         DO IQQP = 1,NQQP_STR
C
            F2 = F1*EXPGNQ(N,IQQP)
            H2 = H1*EXPGNQ_RED(N,IQQP)
C
            MMLL = 0
            DO LL = 0,LLMAX
               F3 = (F2-H2)*D1TERM3(LL)
C
               DO MM = -LL, + LL
                  MMLL = MMLL + 1
                  DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP) + F3*HP(MMLL)
               END DO
            END DO
C
         END DO
C
      END DO
C
C  MULTIPLY THE MISSING PHASE FACTOR    EXP( I ->K * (->Q[I] - ->Q[J]) )
C
      DO IQQP = 2,NQQP_STR
         EXIKQQP = CDEXP
     &             (CI2PI*(KX*QQPX(IQQP)+KY*QQPY(IQQP)+KZ*QQPZ(IQQP)))
C
         DO MMLL = 1,MMLLMAX
            DLLMMKE(MMLL,IQQP) = EXIKQQP*DLLMMKE(MMLL,IQQP)
         END DO
      END DO
C
C  =====================================================================
C                               ********
C                               * DLM2 *
C                               ********
C
C     TABLE TO CALCULATE         EXP( I ->K * ->R[S] )
C
      EXIK1 = CDEXP(CI2PI*(ABAS(1,1)*KX+ABAS(2,1)*KY+ABAS(3,1)*KZ))
      EXIK2 = CDEXP(CI2PI*(ABAS(1,2)*KX+ABAS(2,2)*KY+ABAS(3,2)*KZ))
      EXIK3 = CDEXP(CI2PI*(ABAS(1,3)*KX+ABAS(2,3)*KY+ABAS(3,3)*KZ))
C
      PWEXIK1(0) = 1.0D0
      PWEXIK2(0) = 1.0D0
      PWEXIK3(0) = 1.0D0
C
      DO IPW = 1,R123MAX
         PWEXIK1(IPW) = PWEXIK1(IPW-1)*EXIK1
         PWEXIK1(-IPW) = 1.0D0/PWEXIK1(IPW)
C
         PWEXIK2(IPW) = PWEXIK2(IPW-1)*EXIK2
         PWEXIK2(-IPW) = 1.0D0/PWEXIK2(IPW)
C
         PWEXIK3(IPW) = PWEXIK3(IPW-1)*EXIK3
         PWEXIK3(-IPW) = 1.0D0/PWEXIK3(IPW)
      END DO
C
C  ---------------------------------------------------------------------
C                                                     DIRECT LATTICE SUM
C
      DO IQQP = 1,NQQP_STR
         DO S = 1,SMAX(IQQP)
C
            I = INDR(S,IQQP)
            EXIKRS = PWEXIK1(R1(I))*PWEXIK2(R2(I))*PWEXIK3(R3(I))
C
            CALL ZAXPY(MMLLMAX,EXIKRS,QQMLRS(1,S,IQQP),1,DLLMMKE(1,IQQP)
     &                 ,1)
C
         END DO
      END DO
C
C  =====================================================================
C                               ********
C                               * DLM3 *
C                               ********
C
C             ADD THE MISSING TERM   D300
C
      DLLMMKE(1,1) = DLLMMKE(1,1) + D300
C
C  =====================================================================
C         remove the factor 1/p**l to get smoother variation with E
C  =====================================================================
C
      DO IQQP = 1,NQQP_STR
         DO MMLL = 1,MMLLMAX
            DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)/PWP(MMLL)
         END DO
      END DO
C
      END
C*==strbbddpole.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRBBDDPOLE(DLLMMKE,KX,KY,KZ)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the missing   ->k - dependent quantities              *
C   *  then set up   DLM(K,E)   by performing the lattice sum          *
C   *  the DLM's for the various (IQ,IQ')-blocks of  G  are            *
C   *  evaluated in parallel - as far as possible                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:NGFEP
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_STRLAG,ONLY:EXPGNQ_RED
      USE MOD_STR,ONLY:NLLMMMAX,LLMAX,MMLLMAX,NQQP_STR,HP,QQPX,QQPY,
     &    QQPZ,G1,G2,G3,EDU,PWP
      USE MOD_CONSTANTS,ONLY:C1,CI2PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQP_STR)
C
C Local variables
C
      COMPLEX*16 DENOM,DLM1(:,:),EXIKQQP,H1,H2
      INTEGER IQQP,MMLL,N
      REAL*8 KNX,KNY,KNZ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DLM1
C
      ALLOCATE (DLM1(NLLMMMAX,NQQP_STR))
C
C  =====================================================================
C                     add pole contributions  to DLM1
C  =====================================================================
C
      CALL CINIT(NLLMMMAX*NQQP_STR,DLM1)
C
      DO N = 1,NGFEP
C
C     ->K(N) = ->K + ->KK(N)
C
         KNX = KX + G1(N)*BBAS(1,1) + G2(N)*BBAS(1,2) + G3(N)*BBAS(1,3)
         KNY = KY + G1(N)*BBAS(2,1) + G2(N)*BBAS(2,2) + G3(N)*BBAS(2,3)
         KNZ = KZ + G1(N)*BBAS(3,1) + G2(N)*BBAS(3,2) + G3(N)*BBAS(3,3)
C
         DENOM = (KNX*KNX+KNY*KNY+KNZ*KNZ) - EDU
C
         H1 = C1/DENOM
C
         CALL CALC_RHPLM(KNX,KNY,KNZ,HP,LLMAX,NLLMMMAX)
C
         DO IQQP = 1,NQQP_STR
C
            H2 = H1*EXPGNQ_RED(N,IQQP)
C
            DO MMLL = 1,MMLLMAX
               DLM1(MMLL,IQQP) = DLM1(MMLL,IQQP) + H2*HP(MMLL)
            END DO
C
         END DO
C
      END DO
C
C  MULTIPLY THE MISSING PHASE FACTOR    EXP( I ->K * (->Q[I] - ->Q[J]) )
C
      DO MMLL = 1,MMLLMAX
         DLLMMKE(MMLL,1) = DLLMMKE(MMLL,1) + DLM1(MMLL,1)
      END DO
C
      DO IQQP = 2,NQQP_STR
         EXIKQQP = CDEXP
     &             (CI2PI*(KX*QQPX(IQQP)+KY*QQPY(IQQP)+KZ*QQPZ(IQQP)))
C
         DO MMLL = 1,MMLLMAX
            DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                           + EXIKQQP*DLM1(MMLL,IQQP)
         END DO
      END DO
C
C  ---------------------------------------------------------------------
C
C
C  =====================================================================
C                  add the factor 1/p**l
C  =====================================================================
C
      DO IQQP = 1,NQQP_STR
         DO MMLL = 1,MMLLMAX
            DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)*PWP(MMLL)
         END DO
      END DO
C
      END
C*==strcclag.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRCCLAG(ERYD)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_STR,ONLY:PWP,LLMAX,EDU
      USE MOD_STRLAG,ONLY:NSTRLAG,ESTRLAG,WSTRLAG
      USE MOD_CONSTANTS,ONLY:C1,PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
C
C Local variables
C
      COMPLEX*16 CAUX,DE,PDU
      INTEGER I,J,LL,MM,MMLL
C
C*** End of declarations rewritten by SPAG
C
      EDU = ERYD/(2*PI/ALAT)**2
      PDU = CDSQRT(EDU)
C
      CAUX = PDU
      MMLL = 0
      DO LL = 0,LLMAX
         CAUX = CAUX/PDU
         DO MM = -LL, + LL
            MMLL = MMLL + 1
            PWP(MMLL) = CAUX
         END DO
      END DO
C
      WRITE (*,'(a,3x,2f18.12)') '############## EEE ',ERYD
      CAUX = 0D0
      DO I = 1,NSTRLAG
         WSTRLAG(I) = C1
         DO J = 1,NSTRLAG
            IF ( I.NE.J ) THEN
               DE = (ERYD-ESTRLAG(J))/(ESTRLAG(I)-ESTRLAG(J))
               WSTRLAG(I) = WSTRLAG(I)*DE
            END IF
         END DO
         WRITE (*,'(a,i3,2f18.12)') '############## LAG ',I,WSTRLAG(I)
         CAUX = CAUX + WSTRLAG(I)
      END DO
      WRITE (*,'(a,3x,2f18.12)') '##### SUM #### LAG ',CAUX
C
      END
C*==strbbddlag.f    processed by SPAG 6.70Rc at 22:16 on 20 Dec 2016
      SUBROUTINE STRBBDDLAG(DLLMMKE,IK,KX,KY,KZ)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_STR,ONLY:NLLMMMAX,MMLLMAX,NQQP_STR
      USE MOD_STRLAG,ONLY:NSTRLAG,WSTRLAG,DLMLAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IK
      REAL*8 KX,KY,KZ
      COMPLEX*16 DLLMMKE(NLLMMMAX,NQQP_STR)
C
C Local variables
C
      INTEGER IE,IQQP,MMLL
      COMPLEX*16 WE
C
C*** End of declarations rewritten by SPAG
C
      CALL CINIT(NLLMMMAX*NQQP_STR,DLLMMKE)
C
      DO IE = 1,NSTRLAG
         WE = WSTRLAG(IE)
         DO IQQP = 1,NQQP_STR
            DO MMLL = 1,MMLLMAX
               DLLMMKE(MMLL,IQQP) = DLLMMKE(MMLL,IQQP)
     &                              + WE*DLMLAG(MMLL,IQQP,IE,IK)
            END DO
         END DO
      END DO
C
      CALL STRBBDDPOLE(DLLMMKE,KX,KY,KZ)
C
      END
