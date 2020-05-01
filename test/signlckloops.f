C*==signlckloops.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLCKLOOPS(P,ERYD,TAUQ,TAUQZ,DROT,IQORGQP,SYMUNITARY,
     &                        MSSQ,NSYM,CHIHATZ,IPRINTL,MUEHAT,OMEGAHAT,
     &                        TAUHAT,EIKNRIJ,KTABKN,WKTABKN,NKTABKN,
     &                        SYMACCEPTEDKN,NSYMACCEPTEDKN,LNLCPAAVG,
     &                        SYM_IRR,IQCPA,NQNLCPA,NDIMCLU,IND0QCLU,
     &                        NKN,NKNIRMU,MAXNKNIRMU,NKTABKNMAX,NDIMCHI)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array                    *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *    THIS ROUTINE IS RESTRICTED TO NQMAX = 1 AT THE MOMENT !!!!!   *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_SITES,ONLY:IQBOT_CHI,IQTOP_CHI,NQ,NQMAX
      USE MOD_LINRESP,ONLY:CHIZ,NZ12,NZ12MAX,ITTA,ITTB,ITTC,ITTD,ITTQ1,
     &    ITTQ2,NTKTKLIN
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKKR,NKMQ
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--SIGNLCKLOOPS37
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLCKLOOPS')
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      INTEGER IPRINTL,IQCPA,MAXNKNIRMU,NDIMCHI,NDIMCLU,NKN,NKTABKNMAX,
     &        NQNLCPA,NSYM
      LOGICAL LNLCPAAVG
      COMPLEX*16 CHIHATZ(NDIMCHI,NDIMCHI,NZ12MAX),
     &           DROT(NKMMAX,NKMMAX,NSYMMAX),
     &           EIKNRIJ(NQNLCPA,MAXNKNIRMU,NQNLCPA,NQNLCPA),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MUEHAT(NDIMCLU,NDIMCLU),
     &           OMEGAHAT(NDIMCLU,NDIMCLU),TAUHAT(NDIMCLU,NDIMCLU),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX)
      INTEGER IND0QCLU(NQNLCPA),IQORGQP(NSYMMAX,NQMAX),NKNIRMU(NKN),
     &        NKTABKN(NKN),NSYMACCEPTEDKN(NKN)
      REAL*8 KTABKN(3,NKTABKNMAX,NKN),WKTABKN(NKTABKNMAX,NKN)
      LOGICAL SYMACCEPTEDKN(NSYMMAX,NKN),SYMUNITARY(NSYMMAX),
     &        SYM_IRR(NKN,MAXNKNIRMU,NSYM)
C
C Local variables
C
      COMPLEX*16 CHIKN(:,:,:,:,:),CW11,CW12,CWGT,MAUX(:,:),MSSQKN(:,:,:)
     &           ,MUEHATKN(:,:),SUMQ(:,:,:),TAUHATZ(:,:,:,:,:),
     &           TAUKLIN(:),TAUKN(:,:),TAUKND(:,:),TKTKQQ(:,:),W1(:,:),
     &           W1HAT(:,:),W2HAT(:,:)
      INTEGER I,I0,IA_ERR,II,III,IKN,IKND,IKNM,IQ,IQCLU,ISYM,J,J0,JJ,
     &        JJJ,JQ,JQCLU,K1,K2,KK,L1,L2,L3,L4,LL,M,NDIMCLUSQ,NKMSQ,
     &        NQKMSQ
      LOGICAL LNLCPASYMAVG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUKLIN,SUMQ,TKTKQQ,W1,TAUHATZ
      ALLOCATABLE W1HAT,W2HAT,TAUKND,MUEHATKN,TAUKN,MSSQKN,CHIKN
C
      ALLOCATE (SUMQ(NKMMAX,NKMMAX,NQMAX),W1(NKMMAX,NKMMAX))
      ALLOCATE (TKTKQQ(NTKTKLIN,NZ12))
      ALLOCATE (MAUX(NKKR,NKKR))
      ALLOCATE (TAUKLIN(NKKR*NKKR))
      ALLOCATE (MUEHATKN(NKMMAX,NKMMAX))
      ALLOCATE (TAUHATZ(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA,2))
      ALLOCATE (CHIKN(NKMMAX,NKMMAX,NKMMAX,NKMMAX,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUKLIN')
C
      IF ( NQMAX.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NQMAX <> 1')
      IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL = 2')
      IF ( NZ12.LT.1 .OR. NZ12.GT.2 )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12 not set properly')
      IF ( NZ12.NE.NZ12MAX )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12MAX not set properly')
C
C
C ------------ calculate energy - dependent terms of structure constants
C
      CALL STRCC(ERYD,.FALSE.)
C
C------------------------- assume the same L-expansion for every site IQ
      NKM = NKMQ(1)
      DO IQ = 2,NQ
         IF ( NKMQ(IQ).NE.NKM )
     &         CALL STOP_MESSAGE(ROUTINE,'NKMQ(IQ) <> NKM')
      END DO
C
      M = NKMMAX
      ALLOCATE (MSSQKN(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TAUKN(NKMMAX,NKMMAX),TAUKND(NKMMAX,NKMMAX))
      CALL CINIT(NKMMAX*NKMMAX,TAUKN)
C
      ALLOCATE (W1HAT(NDIMCLU,NDIMCLU))
      ALLOCATE (W2HAT(NDIMCLU,NDIMCLU),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC W1HAT')
C
      NDIMCLUSQ = NDIMCLU*NDIMCLU
      NKMSQ = NKM*NKM
      NQKMSQ = NQNLCPA*NKMSQ
C
      CALL CINIT(NDIMCHI*NDIMCHI*NZ12MAX,CHIHATZ)
      CALL CINIT(NDIMCLUSQ,TAUHAT)
C
C========================================================= KN-LOOP START
      IKND = 0
      DO IKN = 1,NKN
         IKND = IKND + 1
C
C---------------------------------------------------------- MUEHATKN (3)
C
         CALL CINIT(NKMMAX*NKMMAX,MUEHATKN)
C
         DO JQCLU = 1,NQNLCPA
            J0 = IND0QCLU(JQCLU)
            DO IQCLU = 1,NQNLCPA
               I0 = IND0QCLU(IQCLU)
               CWGT = 1D0/EIKNRIJ(IKN,1,IQCLU,JQCLU)/DBLE(NQNLCPA)
               DO J = 1,NKM
                  CALL ZAXPY(NKM,CWGT,MUEHAT(I0+1,J0+J),1,MUEHATKN(1,J),
     &                       1)
               END DO
            END DO
         END DO
C
         IF ( IPRINTL.GT.4 ) THEN
            WRITE (6,'("kn = ",4i5)') IKN
            CALL CMATSTRUCT('MUEHATKN',MUEHATKN,NKM,NKMMAX,0,0,0,1D-8,6)
         END IF
C-------------------------------- add MUEHATKN to MSSQ at the NLCPA site
         DO IQ = IQBOT_CHI,IQTOP_CHI
            M = NKMQ(IQ)
            IF ( IQ.EQ.IQCPA ) THEN
               MSSQKN(1:M,1:M,IQ) = MUEHATKN(1:M,1:M)
            ELSE
               MSSQKN(1:M,1:M,IQ) = MSSQ(1:M,1:M,IQ)
            END IF
         END DO
C
C--------------------------------- do BZ integration over tile for TAUKN
C
         CALL SIGNLBZINT(P,TAUQ,SUMQ,TKTKQQ,TAUKLIN,MAUX,W1,DROT,
     &                   IQORGQP,SYMUNITARY,SYMACCEPTEDKN(1,IKN),MSSQKN,
     &                   WKTABKN(1,IKN),KTABKN(1,1,IKN),NKTABKN(IKN),
     &                   NSYM,NSYMACCEPTEDKN(IKN),CHIKN,NKTABKNMAX)
C
C-----------------------------------------------------------------------
C              set up the matrix   TAUZ and correct
C-----------------------------------------------------------------------
C                      ONLY FOR TESTING
C-----------------------------------------------------------------------
         DO IQ = IQBOT_CHI,IQTOP_CHI
C
            DO J = 1,NKMQ(IQ)
               DO I = 1,NKMQ(IQ)
                  TAUQZ(I,J,IQ,1) = TAUQ(I,J,IQ)
                  TAUQZ(I,J,IQ,2) = DCONJG(TAUQ(J,I,IQ))
               END DO
            END DO
C
         END DO
C
         CALL CINIT(NKMMAX**4*NQMAX**2*NZ12,CHIZ)
C
         DO I = 1,NTKTKLIN
C
            IQ = ITTQ1(I)
            JQ = ITTQ2(I)
            L1 = ITTA(I)
            L2 = ITTB(I)
            L3 = ITTC(I)
            L4 = ITTD(I)
            K1 = (L1-1)*NKM + L4 + (IQ-1)*NKMSQ
            K2 = (L2-1)*NKM + L3 + (JQ-1)*NKMSQ
C
            CHIZ(K1,K2,1) = CHIKN(L1,L2,L3,L4,1)
            CHIZ(K1,K2,2) = CHIKN(L1,L2,L3,L4,2)
C
         END DO
C-----------------------------------------------------------------------
C
         IF ( IPRINTL.GT.3 ) THEN
            WRITE (6,'(70("-"))')
            WRITE (6,'("nlcpa:: TAU for ikn ",9i5)') IKN
            CALL CMATSTRUCT('TAUQ         ',TAUQ(1,1,1),NKM,NKMMAX,IREL,
     &                      IREL,0,1D-8,6)
         END IF
C        ------------------------ obtain tau's for symmetry related kn's
C
         TAUKND(1:NKM,1:NKM) = TAUQ(1:NKM,1:NKM,IQCPA)
C
         DO IKNM = 1,NKNIRMU(IKN)
            IF ( IKNM.GT.1 ) THEN
               IKND = IKND + 1
C        ------------------------------------------- do rotation of tauq
               DO ISYM = 1,NSYM
                  IF ( SYM_IRR(IKN,IKNM-1,ISYM) ) THEN
                     CALL ROTATETAUKN(TAUKND,TAUKN,W1,DROT,NKM,NKMMAX,
     &                                SYMUNITARY,ISYM)
                     IF ( IPRINTL.GT.3 ) THEN
                        WRITE (6,'("nlcpa:: TAU for rotated ikn ",9i5)')
     &                         IKN
                        CALL CMATSTRUCT('TAUQ         ',TAUKN(1,1),NKM,
     &                                  NKMMAX,IREL,IREL,0,1D-8,6)
                     END IF
                     EXIT
                  END IF
               END DO
               TAUQ(1:NKM,1:NKM,IQCPA) = TAUKN(1:NKM,1:NKM)
            END IF
C
            DO JQCLU = 1,NQNLCPA
               J0 = IND0QCLU(JQCLU)
               DO IQCLU = 1,NQNLCPA
                  I0 = IND0QCLU(IQCLU)
C
                  CWGT = EIKNRIJ(IKN,IKNM,IQCLU,JQCLU)/DBLE(NQNLCPA)
                  DO J = 1,NKM
                     CALL ZAXPY(NKM,CWGT,TAUQ(1,J,IQCPA),1,
     &                          TAUHAT(I0+1,J0+J),1)
                  END DO
               END DO
            END DO
C
C***********************************************************************
C     NOTE: indexing of CHI in <SIGNLCKLOOPS>, <LINRESP_VERTEX_NLC> and
C     <SIGNLC1> as well as of  w  in <LINRESP_VERTEX_NLC> have to be
C     consistent with loop sequence (II,LL,JJ,KK,L1,L4,L2,L3)
C***********************************************************************
C                                loops over cluster sites II, LL, JJ, KK
            DO II = 1,NQNLCPA
               DO LL = 1,NQNLCPA
                  DO JJ = 1,NQNLCPA
                     DO KK = 1,NQNLCPA
C
                        CW11 = EIKNRIJ(IKN,IKNM,II,JJ)
     &                         *EIKNRIJ(IKN,IKNM,KK,LL)/DBLE(NQNLCPA)
                        CW12 = EIKNRIJ(IKN,IKNM,II,JJ)
     &                         /EIKNRIJ(IKN,IKNM,KK,LL)/DBLE(NQNLCPA)
C
                        DO I = 1,NTKTKLIN
C
                           IQ = ITTQ1(I)
                           JQ = ITTQ2(I)
                           L1 = ITTA(I)
                           L2 = ITTB(I)
                           L3 = ITTC(I)
                           L4 = ITTD(I)
                           K1 = L4 + (L1-1)*NKM + (LL-1)*NKMSQ + (II-1)
     &                          *NQKMSQ
                           K2 = L3 + (L2-1)*NKM + (KK-1)*NKMSQ + (JJ-1)
     &                          *NQKMSQ
C
                           CHIHATZ(K1,K2,1) = CHIHATZ(K1,K2,1)
     &                        + CW11*CHIKN(L1,L2,L3,L4,1)
C
                           CHIHATZ(K1,K2,2) = CHIHATZ(K1,K2,2)
     &                        + CW12*CHIKN(L1,L2,L3,L4,2)
C
                        END DO
C
                     END DO
                  END DO
               END DO
            END DO
C                                loops over cluster sites II, LL, JJ, KK
C***********************************************************************
C
         END DO
      END DO
C=========================================================== KN-LOOP END
C
C----------------------- enforce symmetry, equalize site diagonal blocks
      IF ( LNLCPAAVG ) THEN
C        --- symetrise the average of site-diagonal blocks
         LNLCPASYMAVG = .TRUE.
C
         CALL NLCPAAVG(TAUHAT,NQNLCPA,IND0QCLU,DROT,SYMUNITARY,NSYM,
     &                 NKMMAX,NKM,IQCPA,NDIMCLU,NQMAX,NKMQ,IQORGQP,
     &                 LNLCPASYMAVG)
      END IF
C-----------------------------------------------------------------------
      IF ( IPRINTL.GT.3 ) THEN
         WRITE (6,'(70("-"))')
         WRITE (6,'("nlcpa:: tauhat (I,J)")')
         CALL CMATSTRUCT('TAUHAT',TAUHAT,NDIMCLU,NDIMCLU,0,0,0,1D-8,6)
      END IF
C
C-------------------------------------------------------------- OMEGAHAT
C
      M = NDIMCLU
C
      W1HAT(1:M,1:M) = TAUHAT(1:M,1:M)
C
      CALL CMATINV(M,M,W1HAT,W2HAT)
C
      OMEGAHAT(1:M,1:M) = MUEHAT(1:M,1:M) - W2HAT(1:M,1:M)
C
C=======================================================================
C
      CALL CINIT(NKMMAX*NKMMAX*NQNLCPA*NQNLCPA*2,TAUHATZ)
C
      JJJ = 0
      DO JJ = 1,NQNLCPA
         DO J = 1,NKM
            JJJ = JJJ + 1
            III = 0
            DO II = 1,NQNLCPA
               DO I = 1,NKM
                  III = III + 1
                  TAUHATZ(I,J,II,JJ,1) = TAUHAT(III,JJJ)
                  TAUHATZ(J,I,JJ,II,2) = DCONJG(TAUHAT(III,JJJ))
               END DO
            END DO
         END DO
      END DO
C
C***********************************************************************
C                                loops over cluster sites II, LL, JJ, KK
      DO II = 1,NQNLCPA
         DO LL = 1,NQNLCPA
            DO JJ = 1,NQNLCPA
               DO KK = 1,NQNLCPA
C
                  DO I = 1,NTKTKLIN
C
                     IQ = ITTQ1(I)
                     JQ = ITTQ2(I)
                     L1 = ITTA(I)
                     L2 = ITTB(I)
                     L3 = ITTC(I)
                     L4 = ITTD(I)
                     K1 = L4 + (L1-1)*NKM + (LL-1)*NKMSQ + (II-1)*NQKMSQ
                     K2 = L3 + (L2-1)*NKM + (KK-1)*NKMSQ + (JJ-1)*NQKMSQ
C
                     CHIHATZ(K1,K2,1) = CHIHATZ(K1,K2,1) - DBLE(NQNLCPA)
     &                                  *TAUHATZ(L1,L2,II,JJ,1)
     &                                  *TAUHATZ(L3,L4,KK,LL,1)
C
                     CHIHATZ(K1,K2,2) = CHIHATZ(K1,K2,2) - DBLE(NQNLCPA)
     &                                  *TAUHATZ(L1,L2,II,JJ,1)
     &                                  *TAUHATZ(L3,L4,KK,LL,2)
C
                  END DO
C
               END DO
            END DO
         END DO
      END DO
C                                loops over cluster sites II, LL, JJ, KK
C***********************************************************************
C
      DEALLOCATE (MAUX,TAUKLIN,SUMQ,TKTKQQ,W1,TAUHATZ)
      DEALLOCATE (W1HAT,W2HAT,TAUKND,MUEHATKN,TAUKN,MSSQKN,CHIKN)
C
      END
C*==signlbzint.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLBZINT(P,TAUQ,SUMQ,TKTKQQ,TAUKLIN,MAUX,W1,DROT,
     &                      IQORGQP,SYMUNITARY,SYMACCEPTED,MSSQ,WKTAB,
     &                      KTAB,NKTAB,NSYM,NSYMACCEPTED,CHIKN,NKTABMAX)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array                    *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *    THIS ROUTINE IS RESTRICTED TO NQMAX = 1 AT THE MOMENT !!!!!   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:NZ12,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,
     &    JTT1,JTT2,JTTX,WTTJ,NTKTKLIN
      USE MOD_ANGMOM,ONLY:NKKR,NKMMAX,IND0Q,NKMQ
      USE MOD_SITES,ONLY:IQBOT_CHI,IQTOP_CHI,NQ,NQMAX
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_SYMMETRY,ONLY:NSYMMAX
      IMPLICIT NONE
C*--SIGNLBZINT423
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLBZINT')
C
C Dummy arguments
C
      INTEGER NKTAB,NKTABMAX,NSYM,NSYMACCEPTED
      COMPLEX*16 P
      COMPLEX*16 CHIKN(NKMMAX,NKMMAX,NKMMAX,NKMMAX,2),
     &           DROT(NKMMAX,NKMMAX,NSYMMAX),MAUX(NKKR,NKKR),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),SUMQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUKLIN(NKKR*NKKR),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TKTKQQ(NTKTKLIN,NZ12),W1(NKMMAX,NKMMAX)
      INTEGER IQORGQP(NSYMMAX,NQMAX)
      REAL*8 KTAB(3,NKTABMAX),WKTAB(NKTABMAX)
      LOGICAL SYMACCEPTED(NSYMMAX),SYMUNITARY(NSYMMAX)
C
C Local variables
C
      COMPLEX*16 CSUM1,CSUMX,CWK,WT
      INTEGER I,I1,IK,INFO,IPIV(NKKR),IQ,J,JQ,JTT,L1,L2,L3,L4,N
      REAL*8 RSCL,WK,WKSUM
C
C*** End of declarations rewritten by SPAG
C
      TAUQ(:,:,:) = C0
      SUMQ(:,:,:) = C0
      TKTKQQ(:,:) = C0
C
      WKSUM = 0D0
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTAB
C
         CALL STRSET(IK,KTAB(1,IK),MAUX,TAUKLIN,P)
C
         CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLIN,MSSQ,NQMAX,NKKR,NKMMAX)
C
         CALL ZGETRF(NKKR,NKKR,TAUKLIN,NKKR,IPIV,INFO)
         CALL ZGETRI(NKKR,TAUKLIN,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C------------------------------------------------------------ store TAUQ
         WK = WKTAB(IK)
         WKSUM = WKSUM + WK
         CWK = DCMPLX(WK,0D0)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            I1 = IND0Q(IQ)*NKKR + IND0Q(IQ) + 1 - NKKR
            N = NKMQ(IQ)
            DO J = 1,N
               I1 = I1 + NKKR
               CALL ZAXPY(N,CWK,TAUKLIN(I1),1,SUMQ(1,J,IQ),1)
            END DO
         END DO
C--------------------------------------------------------- store TAU*TAU
         JTT = 0
         DO I = 1,NTKTKLIN
            CSUM1 = C0
            CSUMX = C0
            DO J = 1,NTTJ(I)
               JTT = JTT + 1
               WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
               CSUM1 = CSUM1 + WT*TAUKLIN(JTT2(JTT))
               CSUMX = CSUMX + WT*DCONJG(TAUKLIN(JTTX(JTT)))
            END DO
            TKTKQQ(I,1) = TKTKQQ(I,1) + CSUM1*WK
            TKTKQQ(I,NZ12) = TKTKQQ(I,NZ12) + CSUMX*WK
         END DO
C-----------------------------------------------------------------------
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      WRITE (6,*) 'calling SYMSUMRTR'
      CALL SYMSUMRTR(.FALSE.,.FALSE.,.FALSE.,WKSUM,SUMQ,TAUQ,W1,NQ,NKMQ,
     &               DROT,IQORGQP,SYMUNITARY,SYMACCEPTED,NSYM,
     &               NSYMACCEPTED,NQMAX,NKMMAX)
C
C-----------------------------------------------------------------------
C            set up the super-matrix  CHIKN(L1,L2,L3,L4)
C-----------------------------------------------------------------------
C
      CALL CINIT(NKMMAX**4*2,CHIKN)
C
      RSCL = 1D0/WKSUM
C
      DO I = 1,NTKTKLIN
C
         L1 = ITTA(I)
         L2 = ITTB(I)
         L3 = ITTC(I)
         L4 = ITTD(I)
         IQ = ITTQ1(I)
         JQ = ITTQ2(I)
C
         CHIKN(L1,L2,L3,L4,1) = RSCL*TKTKQQ(I,1)
         CHIKN(L1,L2,L3,L4,2) = RSCL*TKTKQQ(I,2)
      END DO
C
      IF ( IQ.NE.1 .OR. JQ.NE.1 )
     &      CALL STOP_MESSAGE(ROUTINE,'IQ or JQ <> 1')
C
      END
