C*==linresp_vertex_nlc.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_VERTEX_NLC(CHIHATZ,MUEHATA,MUEHATB,OMEGAHATA,
     &                              OMEGAHATB,MSSTA,MSSTB,NQNLCPA,
     &                              IND0QCLU,ITOQ,NCFG,PCFG,NZ12,IQCPA,
     &                              NKM,NKMQ,NKMMAX,NQMAX,NZ12MAX,
     &                              NDIMCLU,NTMAX,NDIMCHI)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the vertex-corrected CHI matrix                      *
C   *                                                                  *
C   *           chi(vertex) = [ 1 - chi * w ]^-1 * chi                 *
C   *                                                                  *
C   *                                                                  *
C   *   NOTE: the matrix  CHIHATZ  is overwritten                      *
C   *   THIS ROUTINE IS RESTRICTED TO NQMAX = 1 AT THE MOMENT !!!!!    *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C*--LINRESP_VERTEX_NLC19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_VERTEX_NLC')
      LOGICAL CHECK
      PARAMETER (CHECK=.TRUE.)
      REAL*8 T
      PARAMETER (T=1D-8)
C
C Dummy arguments
C
      INTEGER IQCPA,NCFG,NDIMCHI,NDIMCLU,NKM,NKMMAX,NQMAX,NQNLCPA,NTMAX,
     &        NZ12,NZ12MAX
      COMPLEX*16 CHIHATZ(NDIMCHI,NDIMCHI,NZ12MAX),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           MUEHATA(NDIMCLU,NDIMCLU),MUEHATB(NDIMCLU,NDIMCLU),
     &           OMEGAHATA(NDIMCLU,NDIMCLU),OMEGAHATB(NDIMCLU,NDIMCLU)
      INTEGER IND0QCLU(NQNLCPA),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX)
      REAL*8 PCFG(NCFG)
C
C Local variables
C
      INTEGER AA,I,I0,I1,IA_ERR,ICFG,II,INFO,IOCC,IPIV(:),IQCLU,IT,J,J0,
     &        JJ,JQCLU,K,K1,K2,KK,L,L1,L2,L3,L4,LL,LWORK,M,MM,MMQ,N,
     &        NKMSQ,NN,NNQ,NQKMSQ,PP,Z2
      COMPLEX*16 AMAT(:,:),BMAT(:,:),CWGT,CXAXB,TAUIMP(:,:),
     &           W1HAT(NDIMCLU,NDIMCLU),W2HAT(:,:),W3HAT(:,:),
     &           XMATA(:,:,:,:),XMATB(:,:,:,:,:)
      INTEGER NLCPACONF
      CHARACTER*2 ZZ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AMAT,BMAT,XMATA,XMATB,TAUIMP,W2HAT,W3HAT,IPIV
C
      ALLOCATE (W2HAT(NDIMCLU,NDIMCLU),W3HAT(NDIMCLU,NDIMCLU))
      ALLOCATE (IPIV(NDIMCHI))
      ALLOCATE (AMAT(NDIMCHI,NDIMCHI),BMAT(NDIMCHI,NDIMCHI))
      ALLOCATE (TAUIMP(NDIMCLU,NDIMCLU))
      ALLOCATE (XMATA(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA))
      ALLOCATE (XMATB(NKMMAX,NKMMAX,NQNLCPA,NQNLCPA,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate  AMAT')
C
      NKMSQ = NKM*NKM
      NQKMSQ = NQNLCPA*NKMSQ
C
C???????????????????????????????????????????????????????????????????????
      IF ( CHECK ) THEN
         NN = NKMMAX
         MM = NKMMAX
         NNQ = NKMMAX*NQNLCPA
         MMQ = NKMMAX*NQNLCPA
         CALL CMATSTRUCT('MSSTA 1',MSSTA(1,1,1),NN,MM,3,3,0,T,6)
         CALL CMATSTRUCT('MSSTA 2',MSSTA(1,1,2),NN,MM,3,3,0,T,6)
         CALL CMATSTRUCT('MSSTB 1',MSSTB(1,1,1),NN,MM,3,3,0,T,6)
         CALL CMATSTRUCT('MSSTB 2',MSSTB(1,1,2),NN,MM,3,3,0,T,6)
         CALL CMATSTRUCT('MUEHATA',MUEHATA,NNQ,MMQ,0,0,0,T,6)
         CALL CMATSTRUCT('MUEHATB',MUEHATB,NNQ,MMQ,0,0,0,T,6)
         CALL CMATSTRUCT('OMEGAHATA',OMEGAHATA,NNQ,MMQ,0,0,0,T,6)
         CALL CMATSTRUCT('OMEGAHATB',OMEGAHATB,NNQ,MMQ,0,0,0,T,6)
      END IF
C???????????????????????????????????????????????????????????????????????
C
C=======================================================================
C                loop over complex energies
C=======================================================================
C
      DO Z2 = 1,NZ12
C
         AA = 0
C
         IF ( Z2.EQ.1 ) THEN
            ZZ = '  '
         ELSE
            ZZ = ' *'
         END IF
C
         N = NKMQ(IQCPA)
         M = NDIMCLU
C
C=======================================================================
C    set up effective interaction matrix  w  and store in AMAT
C
C                w = sum(ICFG) w(ICFG) * x(ICFG) * x(ICFG)
C=======================================================================
         CALL CINIT(NDIMCHI*NDIMCHI,AMAT)
C
         PP = 1
C
         DO ICFG = 1,NCFG
C
            OPEN (777,FILE='weights')
C
            WRITE (777,*) PP,PCFG(ICFG)
            PP = PP + 1
            CWGT = PCFG(ICFG)/(DBLE(NQNLCPA))
C
C-----------------------------------------------------------------------
C    set up auxilary x-matrix  XMAT  for each configuration  ICFG
C
C              x(ICFG) = [ m(ICFG) - m(NLCPA) ] * D(ICFG)
C-----------------------------------------------------------------------
C----------------------------------------------------- calculate TAU_imp
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               CALL CMATSTRUCT('MSSTA 1',MSSTA(1,1,1),NN,MM,3,3,0,T,6)
               CALL CMATSTRUCT('MSSTA 2',MSSTA(1,1,2),NN,MM,3,3,0,T,6)
               CALL CMATSTRUCT('MUEHATA',MUEHATA,NNQ,MMQ,0,0,0,T,6)
               CALL CMATSTRUCT('OMEGAHATA',OMEGAHATA,NNQ,MMQ,0,0,0,T,6)
            END IF
C???????????????????????????????????????????????????????????????????????
C
            CALL NLCPAPROJ(ICFG,IQCPA,MSSTA,OMEGAHATA,TAUIMP,W1HAT,
     &                     NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,NTMAX,
     &                     NKMMAX)
C
C--- calculate  D from   TAU_imp = D TAU^  and  Omega^ = mue^ - TAU^**-1
            W1HAT(1:M,1:M) = MUEHATA(1:M,1:M) - OMEGAHATA(1:M,1:M)
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               CALL CMATSTRUCT('TAUIMP  A',TAUIMP,NNQ,MMQ,0,0,0,T,6)
               CALL CMATSTRUCT('W1HAT   A',W1HAT,NNQ,MMQ,0,0,0,T,6)
            END IF
C???????????????????????????????????????????????????????????????????????
C
            CALL CMATMUL(M,M,TAUIMP,W1HAT,W2HAT)
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               WRITE (6,'(//,''  CLUSTER CONFIGURATION '',2I4)')
     &                (NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2),IQCLU=1,
     &                NQNLCPA)
               IF ( ICFG.EQ.1 ) THEN
                  CALL CMATSTRUCT('DMATTA 1'//ZZ,W2HAT(1,1),NNQ,MMQ,0,0,
     &                            0,T,6)
               ELSE
                  CALL CMATSTRUCT('DMATTA 2'//ZZ,W2HAT(1,1),NNQ,MMQ,0,0,
     &                            0,T,6)
               END IF
            END IF
C???????????????????????????????????????????????????????????????????????
C
C----------------------------------------------------- calculate  m-mue^
            W1HAT(1:M,1:M) = -MUEHATA(1:M,1:M)
C
            DO IQCLU = 1,NQNLCPA
C
               IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
C
               IT = ITOQ(IOCC,IQCPA)
               I0 = IND0QCLU(IQCLU)
C
               W1HAT(I0+1:I0+N,I0+1:I0+N) = MSSTA(1:N,1:N,IT)
     &            + W1HAT(I0+1:I0+N,I0+1:I0+N)
C
            END DO
C
C------------------------------------------- calculate  x = (m-mue^) * D
            CALL CMATMUL(M,M,W1HAT,W2HAT,W3HAT)
C
            DO IQCLU = 1,NQNLCPA
               I1 = IND0QCLU(IQCLU) + 1
               DO JQCLU = 1,NQNLCPA
                  J0 = IND0QCLU(JQCLU)
                  DO J = 1,N
                     CALL ZCOPY(N,W3HAT(I1,J0+J),1,
     &                          XMATA(1,J,IQCLU,JQCLU),1)
                  END DO
               END DO
            END DO
C
            AA = AA + 1
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               IF ( ICFG.EQ.1 ) THEN
                  DO II = 1,NQNLCPA
                     DO JJ = 1,NQNLCPA
                        WRITE (6,'(//,''  CLUSTER SITES '',2I4)') II,JJ
                        CALL CMATSTRUCT('XMATTA 1'//ZZ,XMATA(1,1,II,JJ),
     &                                  NN,MM,3,3,0,T,6)
                     END DO
                  END DO
               ELSE
                  DO II = 1,NQNLCPA
                     DO JJ = 1,NQNLCPA
                        WRITE (6,'(//,''  CLUSTER SITES '',2I4)') II,JJ
                        CALL CMATSTRUCT('XMATTA 2'//ZZ,XMATA(1,1,II,JJ),
     &                                  NN,MM,3,3,0,T,6)
                     END DO
                  END DO
               END IF
            END IF
C???????????????????????????????????????????????????????????????????????
C
C----------------------------------------------------- calculate TAU_imp
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               CALL CMATSTRUCT('MSSTB 1',MSSTB(1,1,1),NN,MM,3,3,0,T,6)
               CALL CMATSTRUCT('MSSTB 2',MSSTB(1,1,2),NN,MM,3,3,0,T,6)
               CALL CMATSTRUCT('MUEHATB',MUEHATB,NNQ,MMQ,0,0,0,T,6)
               CALL CMATSTRUCT('OMEGAHATB',OMEGAHATB,NNQ,MMQ,0,0,0,T,6)
            END IF
C???????????????????????????????????????????????????????????????????????
C
            CALL NLCPAPROJ(ICFG,IQCPA,MSSTB,OMEGAHATB,TAUIMP,W1HAT,
     &                     NQNLCPA,N,IND0QCLU,ITOQ,NDIMCLU,NQMAX,NTMAX,
     &                     NKMMAX)
C
C--- calculate  D from   TAU_imp = D TAU^  and  Omega^ = mue^ - TAU^**-1
            W1HAT(1:M,1:M) = MUEHATB(1:M,1:M) - OMEGAHATB(1:M,1:M)
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               CALL CMATSTRUCT('TAUIMP  B',TAUIMP,NNQ,MMQ,0,0,0,T,6)
               CALL CMATSTRUCT('W1HAT   B',W1HAT,NNQ,MMQ,0,0,0,T,6)
            END IF
C???????????????????????????????????????????????????????????????????????
C
            CALL CMATMUL(M,M,TAUIMP,W1HAT,W2HAT)
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               WRITE (6,'(//,''  CLUSTER CONFIGURATION '',2I4)')
     &                (NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2),IQCLU=1,
     &                NQNLCPA)
               IF ( ICFG.EQ.1 ) THEN
                  CALL CMATSTRUCT('DMATTB 1'//ZZ,W2HAT(1,1),NNQ,MMQ,0,0,
     &                            0,T,6)
               ELSE
                  CALL CMATSTRUCT('DMATTB 2'//ZZ,W2HAT(1,1),NNQ,MMQ,0,0,
     &                            0,T,6)
               END IF
            END IF
C???????????????????????????????????????????????????????????????????????
C
C----------------------------------------------------- calculate  m-mue^
            W1HAT(1:M,1:M) = -MUEHATB(1:M,1:M)
C
            DO IQCLU = 1,NQNLCPA
C
               IOCC = NLCPACONF(ICFG,NQNLCPA,IQCLU,1,2)
C
               IT = ITOQ(IOCC,IQCPA)
               I0 = IND0QCLU(IQCLU)
C
               W1HAT(I0+1:I0+N,I0+1:I0+N) = MSSTB(1:N,1:N,IT)
     &            + W1HAT(I0+1:I0+N,I0+1:I0+N)
C
            END DO
C
C------------------------------------------- calculate  x = (m-mue^) * D
            CALL CMATMUL(M,M,W1HAT,W2HAT,W3HAT)
C
            DO IQCLU = 1,NQNLCPA
               I1 = IND0QCLU(IQCLU) + 1
               DO JQCLU = 1,NQNLCPA
                  J0 = IND0QCLU(JQCLU)
                  DO J = 1,N
                     CALL ZCOPY(N,W3HAT(I1,J0+J),1,
     &                          XMATB(1,J,IQCLU,JQCLU,1),1)
                  END DO
               END DO
            END DO
C
            IF ( Z2.EQ.2 ) THEN
               DO IQCLU = 1,NQNLCPA
                  DO JQCLU = 1,NQNLCPA
                     DO I = 1,N
                        DO J = 1,N
                           XMATB(I,J,IQCLU,JQCLU,2)
     &                        = DCONJG(XMATB(J,I,JQCLU,IQCLU,1))
                        END DO
                     END DO
                  END DO
               END DO
            END IF
C
C???????????????????????????????????????????????????????????????????????
            IF ( CHECK ) THEN
               IF ( ICFG.EQ.1 ) THEN
                  DO II = 1,NQNLCPA
                     DO JJ = 1,NQNLCPA
                        WRITE (6,'(//,''  CLUSTER SITES '',2I4)') II,JJ
                        CALL CMATSTRUCT('XMATTB 1'//ZZ,
     &                                  XMATB(1,1,II,JJ,Z2),NN,MM,3,3,0,
     &                                  T,6)
                     END DO
                  END DO
               ELSE
                  DO II = 1,NQNLCPA
                     DO JJ = 1,NQNLCPA
                        WRITE (6,'(//,''  CLUSTER SITES '',2I4)') II,JJ
                        CALL CMATSTRUCT('XMATTB 2'//ZZ,
     &                                  XMATB(1,1,II,JJ,Z2),NN,MM,3,3,0,
     &                                  T,6)
                     END DO
                  END DO
               END IF
            END IF
C???????????????????????????????????????????????????????????????????????
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
                        K1 = (LL-1)*NKMSQ + (II-1)*NQKMSQ
                        DO L1 = 1,NKM
                           DO L4 = 1,NKM
                              K1 = K1 + 1
C
                              K2 = (KK-1)*NKMSQ + (JJ-1)*NQKMSQ
                              DO L2 = 1,NKM
                                 DO L3 = 1,NKM
                                    K2 = K2 + 1
C
                                    CXAXB = CWGT*XMATA(L1,L2,II,JJ)
     &                                 *XMATB(L3,L4,KK,LL,Z2)
C
                                    AMAT(K1,K2) = AMAT(K1,K2) + CXAXB
C
                                 END DO
                              END DO
                           END DO
                        END DO
C
                     END DO
                  END DO
               END DO
            END DO
C                                loops over cluster sites II, JJ, KK, LL
C***********************************************************************
C
         END DO
C=======================================================================
C
         IF ( CHECK ) WRITE (6,*) 'XXXXXX  INVERTING BIG MATRIX  XXXXXX'
C
C
         IF ( Z2.EQ.1 ) THEN
C
            OPEN (999,FILE='chihatSTART')
C
C
            K = 0
            L = 0
C
            DO I = 1,NDIMCHI
               DO J = 1,NDIMCHI
                  IF ( ABS(CHIHATZ(I,J,1)).GT.1D-20 ) K = K + 1
               END DO
            END DO
C
C
            DO I = 1,NDIMCHI
               DO J = 1,NDIMCHI
                  IF ( ABS(CHIHATZ(I,J,2)).GT.1D-20 ) L = L + 1
               END DO
            END DO
C
C
            WRITE (999,*) 'Anzahl der Elemente > 1D-20 (z=1)',K
            WRITE (999,*) 'Anzahl der Elemente > 1D-20 (z=2)',L
            WRITE (999,*) 'Anzahl aller Elemente zu einer Energie',
     &                    NDIMCHI*NDIMCHI
C
C
C          DO I = 1,NDIMCHI
C             DO J = 1,NDIMCHI
C                WRITE (999,'(2I6,4e17.10)')i,j,
C     &                                CHIHATZ(I,J,1),CHIHATZ(I,J,2)
C             END DO
C          END DO
C
C
         END IF
C
C
         CALL CMATMUL(NDIMCHI,NDIMCHI,CHIHATZ(1,1,Z2),AMAT,BMAT)
C
         DO I = 1,NDIMCHI
            BMAT(I,I) = -1D0 + BMAT(I,I)
         END DO
C
C        CALL CINVLU(BMAT,AMAT,NDIMCHI,NDIMCHI)
C
         LWORK = NDIMCHI*NDIMCHI
C
         CALL ZGETF2(NDIMCHI,NDIMCHI,BMAT,NDIMCHI,IPIV,INFO)
         CALL ZGETRI(NDIMCHI,BMAT,NDIMCHI,IPIV,AMAT,LWORK,INFO)
C
         CALL CMATMUL(NDIMCHI,NDIMCHI,BMAT,CHIHATZ(1,1,Z2),AMAT)
C
         DO J = 1,NDIMCHI
            DO I = 1,NDIMCHI
               CHIHATZ(I,J,Z2) = -AMAT(I,J)
            END DO
         END DO
C
C
         IF ( Z2.EQ.2 ) THEN
C
            OPEN (888,FILE='chihatEND')
C
C
            K = 0
            L = 0
C
            DO I = 1,NDIMCHI
               DO J = 1,NDIMCHI
                  IF ( ABS(CHIHATZ(I,J,1)).GT.1D-20 ) K = K + 1
               END DO
            END DO
C
C
            DO I = 1,NDIMCHI
               DO J = 1,NDIMCHI
                  IF ( ABS(CHIHATZ(I,J,2)).GT.1D-20 ) L = L + 1
               END DO
            END DO
C
C
            WRITE (888,*) 'Anzahl der Elemente > 1D-20 (z=1)',K
            WRITE (888,*) 'Anzahl der Elemente > 1D-20 (z=2)',L
            WRITE (888,*) 'Anzahl aller Elemente zu einer Energie',
     &                    NDIMCHI*NDIMCHI
C
C
C          DO I = 1,NDIMCHI
C             DO J = 1,NDIMCHI
C                WRITE (888,'(2I6,4e17.10)')i,j,
C     &                                CHIHATZ(I,J,1),CHIHATZ(I,J,2)
C             END DO
C          END DO
C
C
         END IF
C
C
      END DO
C
C
C
C99001 FORMAT ('GAMMMA    I   J   L1  L2 ')
C
      END
