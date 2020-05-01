C*==lloydpt1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LLOYDPT1(PHAST,PHASA,TAUQ,MSSQ,MSST,SSST,ERYD,IE)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the atom specific terms of the Lloyd formula         *
C   *   PHAST   connected with TAU- and t-matrices                     *
C   *   PHASA   connected with the alpha-matrix                        *
C   *                                                                  *
C   *   files opened for test purposes  (ITEST=5)                      *
C   *            OPEN (60,FILE='lloyd_a+b.dat')                        *
C   *            OPEN (61,FILE='lloyd_c.dat')                          *
C   *            OPEN (62,FILE='lloyd_d.dat')                          *
C   *            OPEN (63,FILE='lloyd_a+b+c+d.dat')                    *
C   *            OPEN (64,FILE='lloyd_nos.dat')                        *
C   *            OPEN (77,FILE='lloyd_A+B.dat')                        *
C   *            OPEN (78,FILE='lloyd_A.dat')                          *
C   *            OPEN (79,FILE='lloyd_B.dat')                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:SPLITSS
      USE MOD_ANGMOM,ONLY:NKMQ,NKMMAX
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ
      USE MOD_TYPES,ONLY:NT,NTMAX
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CALCMODE,ONLY:ITEST,LLOYD,GF_CONV_RH,IREL
      IMPLICIT NONE
C*--LLOYDPT128
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IE
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASA(NKMMAX,NTMAX,IE),PHAST(NKMMAX,NTMAX,IE),
     &           SSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 CSUMAB,CSUMC,CSUMD,MSSQLU(NKMMAX,NKMMAX),PHASQ(NKMMAX),
     &           TAUQINV(NKMMAX,NKMMAX),TAUQLU(NKMMAX,NKMMAX),
     &           W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
      INTEGER I,IO,IQ,IT,M,N
      LOGICAL KDONET(NTMAX)
C
C*** End of declarations rewritten by SPAG
C
C55555555555555555555555555555555555555555555555555555555555555555555555
C                          ITEST = 5
C55555555555555555555555555555555555555555555555555555555555555555555555
      IF ( IE.EQ.0 ) THEN
         LLOYD = .TRUE.
         IPRINT = 1
C
         OPEN (60,FILE='lloyd_a+b.dat')
         OPEN (61,FILE='lloyd_c.dat')
         OPEN (62,FILE='lloyd_d.dat')
         OPEN (63,FILE='lloyd_a+b+c+d.dat')
         OPEN (64,FILE='lloyd_nos.dat')
         OPEN (77,FILE='lloyd_A+B.dat')
         OPEN (78,FILE='lloyd_A.dat')
         OPEN (79,FILE='lloyd_B.dat')
C
         WRITE (77,'(A)') '# lloyd_A+B      '
         WRITE (77,'(A)') '# lloyd_A+B    - int(k) det[ TAU(k,E) ]'
         WRITE (77,'(A)') '# lloyd_A+B    sum term A B'
         WRITE (77,'(A)') '# lloyd_A+B      '
C
         WRITE (78,'(A)') '# lloyd_A      '
         WRITE (78,'(A)') '# lloyd_A    free electron term A'
         WRITE (78,'(A)') '# lloyd_A      '
C
         WRITE (79,'(A)') '# lloyd_B      '
         WRITE (79,'(A)') '# lloyd_B    gt-term   term B '
         WRITE (79,'(A)') '# lloyd_B      '
C
         RETURN
      END IF
C55555555555555555555555555555555555555555555555555555555555555555555555
C
      PHASA(1:NKMMAX,1:NTMAX,IE) = 0D0
      PHAST(1:NKMMAX,1:NTMAX,IE) = 0D0
C
      DO IT = 1,NT
         KDONET(IT) = .FALSE.
      END DO
C
      DO IQ = 1,NQ
         IF ( .NOT.KDONET(ITOQ(1,IQ)) ) THEN
C
            N = NKMQ(IQ)
            M = NKMMAX
C
            TAUQLU(1:M,1:M) = TAUQ(1:M,1:M,IQ)
            MSSQLU(1:M,1:M) = MSSQ(1:M,1:M,IQ)
C
            IF ( IREL.EQ.2 .AND. N.LT.NKMMAX ) THEN
               CALL CINVLU_SPSREL(TAUQLU,TAUQINV,N,NKMMAX)
               CALL CINVLU_SPSREL(MSSQLU,W1,N,NKMMAX)
            ELSE
               CALL CINVLU(TAUQLU,TAUQINV,N,NKMMAX)
               CALL CINVLU(MSSQLU,W1,N,NKMMAX)
            END IF
C
            CSUMAB = 0D0
            DO I = 1,N
               PHASQ(I) = LOG(-TAUQLU(I,I)) + LOG(MSSQLU(I,I))
               CSUMAB = CSUMAB + PHASQ(I)
            END DO
            IF ( ITEST.EQ.5 ) WRITE (60,'(2f10.5,3x,2f10.5)')
     &                               DREAL(ERYD),DIMAG(CSUMAB)
C
            DO IO = 1,NOQ(IQ)
C
C----------------------------------------------------------- ln |1/TAUT|
C
               IT = ITOQ(IO,IQ)
C
               IF ( SPLITSS ) THEN
                  PHASQ(1:M) = 0D0
                  W1(1:M,1:M) = MSST(1:M,1:M,IT)
               ELSE
                  W1(1:M,1:M) = TAUQINV(1:M,1:M) - MSSQ(1:M,1:M,IQ)
     &                          + MSST(1:M,1:M,IT)
               END IF
C
               IF ( IREL.EQ.2 .AND. N.LT.NKMMAX ) THEN
                  CALL CINVLU_SPSREL(W1,W2,N,NKMMAX)
               ELSE
                  CALL CINVLU(W1,W2,N,NKMMAX)
               END IF
C
               CSUMC = 0D0
               DO I = 1,N
                  PHAST(I,IT,IE) = PHASQ(I) + LOG(-W1(I,I))
                  CSUMC = CSUMC + LOG(-W1(I,I))
               END DO
C
C----------------------------------------------------- alpha matrix SSST
               W1(1:M,1:M) = SSST(1:M,1:M,IT)
C
               IF ( IREL.EQ.2 .AND. N.LT.NKMMAX ) THEN
                  CALL CINVLU_SPSREL(W1,W2,N,NKMMAX)
               ELSE
                  CALL CINVLU(W1,W2,N,NKMMAX)
               END IF
C
               CSUMD = 0D0
               DO I = 1,N
                  PHASA(I,IT,IE) = LOG(W1(I,I))
                  CSUMD = CSUMD + PHASA(I,IT,IE)
               END DO
C
C----------------------------------------- convention for Green function
C                                                basis functions R and H
               IF ( GF_CONV_RH ) THEN
C
                  W1(1:M,1:M) = MSST(1:M,1:M,IT)
C
                  IF ( IREL.EQ.2 .AND. N.LT.NKMMAX ) THEN
                     CALL CINVLU_SPSREL(W1,W2,N,NKMMAX)
                  ELSE
                     CALL CINVLU(W1,W2,N,NKMMAX)
                  END IF
C
                  DO I = 1,N
                     PHASA(I,IT,IE) = PHASA(I,IT,IE) + LOG(W1(I,I))
                     CSUMD = CSUMD + PHASA(I,IT,IE)
                  END DO
C
               END IF
C
               IF ( ITEST.EQ.5 ) WRITE (61,'(2f10.5,3x,2f10.5)')
     &                                  DREAL(ERYD),DIMAG(CSUMC),
     &                                  DIMAG(CSUMD)
               KDONET(IT) = .TRUE.
            END DO
C
         END IF
      END DO
C
      END
C*==sbranch.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE SBRANCH(A,N,NT,NKMMAX,NTMAX,NE)
C   ********************************************************************
C   *                                                                  *
C   *    Select a proper branch so as to keep continuity.              *
C   *    coded by H.Akai, 1986, Osaka                                  *
C   *    adapted to SPRKKR-package                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--SBRANCH208
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER N,NE,NKMMAX,NT,NTMAX
      COMPLEX*16 A(NKMMAX,NTMAX,NE)
C
C Local variables
C
      INTEGER I,IE,IT,JUMP
C
C*** End of declarations rewritten by SPAG
C
      DO IT = 1,NT
         DO I = 1,N
            JUMP = INT((DIMAG(A(I,IT,1))+PI)/PI/2D0+1.0000001D5)
            JUMP = JUMP - 100000
            A(I,IT,1) = A(I,IT,1) - DBLE(JUMP)*2D0*PI*(0D0,1D0)
         END DO
      END DO
C
      DO IE = 2,NE
         DO IT = 1,NT
            DO I = 1,N
               JUMP = INT((DIMAG(A(I,IT,IE)-A(I,IT,IE-1))+PI)
     &                /PI/2D0+1.0000001D5)
               JUMP = JUMP - 100000
               A(I,IT,IE) = A(I,IT,IE) - DBLE(JUMP)*2D0*PI*(0D0,1D0)
            END DO
         END DO
      END DO
      END
