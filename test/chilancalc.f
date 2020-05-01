C*==chilancalc.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHILANCALC(IECURR,ERYD,P,CHILANLT,TSSQ,MSSQ,TSST,MSST,
     &                      TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the LANDAU susceptibility                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:NETAB,IGRID,WETAB,NEPATH
      USE MOD_CALCMODE,ONLY:KMROT
      USE MOD_TAUIJ,ONLY:TAUIJ,N5VEC_TAUIJ,ITAUIJ_LOOP,ITAUJI_LOOP,
     &    IQ_TAUIJ_CL,NQCLU_TAUIJ_CL
      USE MOD_LATTICE,ONLY:ABAS,ALAT
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_TYPES,ONLY:NT,NLT,NAT,NTMAX,TXT_T,CONC
      USE MOD_SITES,ONLY:QBAS,NQMAX,NOQ,ITOQ,IQAT,ICPA
      USE MOD_CONSTANTS,ONLY:C1,C0,CHI_AU2CGS
      USE MOD_ANGMOM,ONLY:NL,NKM,NLMAX,NKMQ,NKMMAX,TXT_L
      USE MOD_SYMMETRY,ONLY:NCL,IQ_MBCL
      IMPLICIT NONE
C*--CHILANCALC22
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHILANCALC')
      REAL*8 FLAN,CU
      PARAMETER (FLAN=1D0,CU=CHI_AU2CGS*1D+6)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      INTEGER IECURR
      REAL*8 CHILANLT(0:NLMAX,NTMAX,2)
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 AMEA1(:,:,:),AMEA2(:,:,:),AMEB1(:,:,:,:),AMEB1C(:,:,:,:),
     &       AMEB2(:,:,:,:),AMEB2C(:,:,:,:),CHILAN,CPACHNG,DQY,DRY,FA,FB
      COMPLEX*16 CSUM1,CSUM2,DMAMC(:,:),DMATT(:,:),DTILT(:,:),ME(:,:,:),
     &           MEA(:,:,:,:),MEB(:,:,:,:,:),TMT(:,:,:),W1(:,:),WE
      INTEGER I,IA_ERR,ICL,ICPACONV,ICPAFLAG,IL,ILM,ILOOP,IO,IPOL,IPOL1,
     &        IPOL2,IQ,IT,ITAUIJ,ITAUJI,ITCPA,J,JO,JQ,JQCLU,JT,M,MSQ,N,
     &        N1,N2,N3
      SAVE AMEA1,AMEA2,AMEB1,AMEB1C,AMEB2,AMEB2C
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1,ME,TMT,MEA,MEB,DMAMC,DTILT,DMATT
      ALLOCATABLE AMEA1,AMEA2,AMEB1,AMEB2,AMEB1C,AMEB2C
C
      IF ( KMROT.NE.0 ) RETURN
      IF ( NEPATH.NE.1 ) RETURN
      IF ( IGRID(1).NE.5 ) RETURN
C
      M = NKMMAX
C
      ALLOCATE (W1(M,M),ME(M,M,2),TMT(M,M,2))
      ALLOCATE (DMAMC(M,M),DMATT(M,M),DTILT(M,M))
      ALLOCATE (MEA(M,M,3,NT),MEB(M,M,3,3,NT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEA')
C
C-----------------------------------------------------------------------
C          initialize angular matrix elements    if  IECURR = 1
C-----------------------------------------------------------------------
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (AMEA1(M,M,3),AMEB1(M,M,3,3))
         ALLOCATE (AMEA2(M,M,3),AMEB2(M,M,3,3))
         ALLOCATE (AMEB1C(M,M,3,3),AMEB2C(M,M,3,3),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEB1')
C
         CALL CHILANAME(AMEA1,AMEA2,AMEB1,AMEB2,AMEB1C,AMEB2C)
C
         CALL RINIT((1+NLMAX)*NTMAX*2,CHILANLT)
C
      END IF
C
C-----------------------------------------------------------------------
C                calculate site off-diagonal TAUIJ
C-----------------------------------------------------------------------
C
      CALL TAUIJ_DRIVE(IECURR,ERYD,P,TSSQ,MSSQ,TSST,MSST,TAUQ,ICPAFLAG,
     &                 CPACHNG,ITCPA,ICPACONV)
C
C-----------------------------------------------------------------------
C
      DO IT = 1,NT
C
         IQ = IQAT(1,IT)
         N = NKMQ(IQ)
         M = NKMMAX
C
C--------------------------------------------- calculate matrix elements
C
         CALL CHILANME(N,IFILCBWF,IT,AMEA1,AMEA2,AMEB1,AMEB2,AMEB1C,
     &                 AMEB2C,MEA,MEB)
C
C--------------------------- calculate projection matrices DMAT and DTIL
C------------------------ and apply projection matrices in case of alloy
C
         IF ( ICPA(IQ).NE.0 ) THEN
C
            CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,MSSQ(1,1,IQ),
     &                   MSST(1,1,IT),M)
C
            DO IPOL = 1,3
               CALL CMATMUL(N,M,MEA(1,1,IPOL,IT),DMATT,W1)
               CALL CMATMUL(N,M,DTILT,W1,MEA(1,1,IPOL,IT))
            END DO
C
            DO IPOL1 = 1,3
               DO IPOL2 = 1,3
                  CALL CMATMUL(N,M,MEB(1,1,IPOL1,IPOL2,IT),DMATT,W1)
                  CALL CMATMUL(N,M,DTILT,W1,MEB(1,1,IPOL1,IPOL2,IT))
               END DO
            END DO
         END IF
C
      END DO
C
      WE = WETAB(IECURR,1)*FLAN
C
C=======================================================================
C=======================================================================
C
      N = NKM
      M = NKMMAX
      MSQ = M*M
C
      ILOOP = 0
C
      DO ICL = 1,NCL
C
         IQ = IQ_MBCL(1,ICL)
C
         DO JQCLU = 1,NQCLU_TAUIJ_CL(ICL)
C
            JQ = IQ_TAUIJ_CL(JQCLU,ICL)
C
            ILOOP = ILOOP + 1
            ITAUIJ = ITAUIJ_LOOP(ILOOP)
            ITAUJI = ITAUJI_LOOP(ILOOP)
C
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
               CALL CINIT(MSQ*2,TMT)
C
               DQY = QBAS(2,JQ) - QBAS(2,IQ)
C
               N1 = N5VEC_TAUIJ(1,ITAUIJ)
               N2 = N5VEC_TAUIJ(2,ITAUIJ)
               N3 = N5VEC_TAUIJ(3,ITAUIJ)
C
               DRY = N1*ABAS(2,1) + N2*ABAS(2,2) + N3*ABAS(2,3) + DQY
C
               CALL CINIT(MSQ*2,ME)
C
               DO JO = 1,NOQ(JQ)
                  JT = ITOQ(JO,JQ)
C
                  FA = CONC(JT)*0.5D0*DRY*DRY*ALAT*ALAT
                  FB = CONC(JT)*DRY*ALAT
C
                  DO J = 1,N
                     DO I = 1,N
                        ME(I,J,1) = ME(I,J,1) + FA*MEA(I,J,1,JT)
                        ME(I,J,2) = ME(I,J,2) - FB*MEB(I,J,2,1,JT)
                     END DO
                  END DO
               END DO
C
               CALL CMATMUL(N,M,ME(1,1,1),TAUIJ(1,1,1,ITAUJI),W1)
C
               CALL ZGEMM('N','N',N,N,N,C1,TAUIJ(1,1,1,ITAUIJ),M,W1,M,
     &                    C1,TMT(1,1,1),M)
C
               CALL CMATMUL(N,M,ME(1,1,2),TAUIJ(1,1,1,ITAUJI),W1)
C
               CALL ZGEMM('N','N',N,N,N,C1,TAUIJ(1,1,1,ITAUIJ),M,W1,M,
     &                    C1,TMT(1,1,2),M)
C
            END DO
C
            I = 0
            DO IL = 1,NLT(IT)
C
               CSUM1 = C0
               CSUM2 = C0
               DO ILM = 1,2*(2*IL-1)
                  I = I + 1
                  DO J = 1,N
                     CSUM1 = CSUM1 + MEA(I,J,1,IT)*TMT(J,I,1)
                     CSUM2 = CSUM2 + MEA(I,J,1,IT)*TMT(J,I,2)
                  END DO
               END DO
C
               CHILANLT(IL,IT,1) = CHILANLT(IL,IT,1) + DIMAG(WE*CSUM1)
               CHILANLT(IL,IT,2) = CHILANLT(IL,IT,2) + DIMAG(WE*CSUM2)
C
            END DO
C
         END DO
C
      END DO
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (MEA,MEB,W1,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 )
     &     CALL STOP_MESSAGE(ROUTINE,'DEALLOC: MEA,MEB,W1')
C
      IF ( IECURR.LT.-NETAB(1) ) RETURN
C
C ======================================================================
C                        print Landau susceptibility
C ======================================================================
C
      CHILAN = 0D0
      DO IT = 1,NT
         CHILANLT(0,IT,1) = 0D0
         CHILANLT(0,IT,2) = 0D0
C
         DO IL = 1,NLT(IT)
            CHILANLT(0,IT,1) = CHILANLT(0,IT,1) + CHILANLT(IL,IT,1)
            CHILANLT(0,IT,2) = CHILANLT(0,IT,2) + CHILANLT(IL,IT,2)
         END DO
C
         CHILAN = CHILAN + NAT(IT)*CONC(IT)
     &            *(CHILANLT(0,IT,1)+CHILANLT(0,IT,2))
      END DO
C
      WRITE (6,99001) (TXT_L(IL),IL=1,NL)
      DO IT = 1,NT
         WRITE (6,99002) IT,TXT_T(IT),
     &                   (CHILANLT(IL,IT,1)*CU,IL=0,NLT(IT))
      END DO
      DO IT = 1,NT
         WRITE (6,99002) IT,TXT_T(IT),
     &                   (CHILANLT(IL,IT,2)*CU,IL=0,NLT(IT))
      END DO
      IF ( NT.GT.1 ) WRITE (6,99003) CHILAN*CU
      WRITE (6,99004)
C
99001 FORMAT (/,1X,79('='),//,10X,
     &        'Landau magnetic susceptibility  in [10^(-6) cm^3/mol]',
     &        //,1X,79('='),//,10X,'type',9X,'sum         ',5(A,:,10X))
99002 FORMAT (10X,I2,1X,A,6F20.10)
99003 FORMAT (10X,'total  ',6F11.4)
99004 FORMAT (/,1X,79('='),/)
      END
