C*==chitktkq2t.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHITKTKQ2T(MSSQ,MSST,TAUQ,DMATT,DTILT)
C   ********************************************************************
C   *                                                                  *
C   *   convert the data set CHIZ(IQ,IQ') to TKTKTT(IT,IT')            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:CHIZ,DDTAUTAUT,TKTKTT,NTKTKLIN,NTKTKMAX,
     &    IKM1_CHI_LIN,IKM2_CHI_LIN,IKM3_CHI_LIN,IKM4_CHI_LIN,LAMCOLROW,
     &    NCOLROW,NLIN23_CHI,NLIN41_CHI
      USE MOD_SITES,ONLY:NQ,NQMAX,IQAT,ICPA
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC,NAT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,NKMQ
      USE MOD_CPA,ONLY:NCPA
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:C1
      IMPLICIT NONE
C*--CHITKTKQ2T19
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 DMATT(NKMMAX,NKMMAX,NTMAX),DTILT(NKMMAX,NKMMAX,NTMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 DA,DAB,DABC,DABCD,DMAMC(NKMMAX,NKMMAX),EA,EAB,EABC,
     &           EABCD,MUNIT(NKMMAX,NKMMAX)
      INTEGER I,IA,IAT,IB,IC,ID,IDA,IDB,IDC,IDD,INDIQ,INDJQ,IQ,IT,JAT,
     &        JQ,JT,K1,K2,K3,K4,KADI,KBCJ,LIN23,LIN41,LINTKTK,M,N,NKMSQ
C
C*** End of declarations rewritten by SPAG
C
      NKMSQ = NKM*NKM
C ---------------------------------------------------------- UNIT MATRIX
      CALL CINIT(NKMMAX*NKMMAX,MUNIT)
      DO I = 1,NKMMAX
         MUNIT(I,I) = C1
      END DO
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,99001) NTKTKMAX
         DO IQ = 1,NQ
            WRITE (6,99002) NTKTKLIN
         END DO
         WRITE (6,99003) NCPA,(ICPA(IQ),IQ=1,NQ)
      END IF
C
C tktkqq -> tktktt================================================ IT ==
      DO IT = 1,NT
C
         IQ = IQAT(1,IT)
C
         IF ( ICPA(IQ).EQ.0 ) THEN
C
            CALL CINIT(NKMMAX*NKMMAX,DTILT(1,1,IT))
            CALL CINIT(NKMMAX*NKMMAX,DMATT(1,1,IT))
C
            DO I = 1,NKMQ(IQ)
               DTILT(I,I,IT) = C1
               DMATT(I,I,IT) = C1
            END DO
C
         ELSE
C
            M = NKMMAX
            N = NKMQ(IQ)
C
            CALL GETDELM(DMAMC,IT,IQ,NKMQ,MSSQ,MSST,NTMAX,NQMAX,NKMMAX)
C
            CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT),DMAMC,
     &                   N,MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
         END IF
C
      END DO
C================================================================= IT ==
C
      CALL CINIT(NTKTKMAX*NTMAX*NTMAX,TKTKTT)
      CALL CINIT(NTKTKMAX*NTMAX,DDTAUTAUT)
C
C================================================================= IT ==
      DO IT = 1,NT
C--------------------------------------------------------------------IAT
         DO IAT = 1,NAT(IT)
            IQ = IQAT(IAT,IT)
            INDIQ = (IQ-1)*NKMSQ
C
C================================================================= JT ==
            DO JT = 1,NT
C--------------------------------------------------------------------jAT
               DO JAT = 1,NAT(JT)
C
                  JQ = IQAT(JAT,JT)
                  INDJQ = (JQ-1)*NKMSQ
C
                  LINTKTK = 0
C
C---------------------------------------------------------------------41
                  DO LIN41 = 1,NLIN41_CHI
C
                     K4 = IKM4_CHI_LIN(LIN41)
                     K1 = IKM1_CHI_LIN(LIN41)
C
C---------------------------------------------------------------------23
                     DO LIN23 = 1,NLIN23_CHI
C
                        K2 = IKM2_CHI_LIN(LIN23)
                        K3 = IKM3_CHI_LIN(LIN23)
C
                        LINTKTK = LINTKTK + 1
C
C===================================================================== A
                        DO IDA = 1,NCOLROW(K1,IQ)
                           IA = LAMCOLROW(IDA,K1,IQ)
                           DA = CONC(JT)*DMATT(K1,IA,IT)
                           EA = DMATT(K1,IA,IT)
C--------------------------------------------------------------------- B
                           DO IDB = 1,NCOLROW(K2,JQ)
                              IB = LAMCOLROW(IDB,K2,JQ)
                              DAB = DA*DTILT(IB,K2,JT)
                              EAB = EA*MUNIT(IB,K2)
C--------------------------------------------------------------------- C
                              DO IDC = 1,NCOLROW(K3,JQ)
                                 IC = LAMCOLROW(IDC,K3,JQ)
                                 DABC = DAB*DMATT(K3,IC,JT)
                                 EABC = EAB*MUNIT(K3,IC)
C--------------------------------------------------------------------- D
                                 DO IDD = 1,NCOLROW(K4,IQ)
                                    ID = LAMCOLROW(IDD,K4,IQ)
                                    DABCD = DABC*DTILT(ID,K4,IT)
                                    EABCD = EABC*DTILT(ID,K4,IT)
C
                                    KADI = (IA-1)*NKM + ID + INDIQ
                                    KBCJ = (IB-1)*NKM + IC + INDJQ
C
                                    TKTKTT(LINTKTK,IT,JT)
     &                                 = TKTKTT(LINTKTK,IT,JT)
     &                                 + DABCD*CHIZ(KADI,KBCJ,1)
     &                                 /DBLE(NAT(IT))
C
                                    IF ( (IT.EQ.JT) .AND. (IQ.EQ.JQ) )
     &                                 DDTAUTAUT(LINTKTK,IT)
     &                                 = DDTAUTAUT(LINTKTK,IT)
     &                                 + EABCD*TAUQ(IA,IB,IQ)
     &                                 *TAUQ(IC,ID,IQ)/DBLE(NAT(IT))
C
                                 END DO
                              END DO
                           END DO
                        END DO
C========================================================= A = B = C = D
C
                     END DO
C---------------------------------------------------------------------23
                  END DO
C---------------------------------------------------------------------41
                  IF ( (IPRINT.GT.0) ) WRITE (6,99004) IQ,JQ,NTKTKLIN
               END DO
C---------------------------------------------------------------------JQ
            END DO
C================================================================= JT ==
C
         END DO
C--------------------------------------------------------------------IAT
C
C   normalisation not consistent for many atoms/unit cell
C   NOW: normalisation done within the loops above (SM)
C         IF ( NAT(IT).GT.1 ) THEN
C            NORM = 1.0D0/DBLE(NAT(IT))
C           DO JT = 1,NT
C               DO I = 1,LINTKTK
C                  TKTKTT(I,IT,JT) = TKTKTT(I,IT,JT)*NORM
C               END DO
C            END DO
C         END IF
C
      END DO
C================================================================= IT ==
C
99001 FORMAT (/,10X,'<CHITKTKQ2T>',/,10X,'NTKTKMAX ',I12)
99002 FORMAT (10X,'NTKTKLIN ',I12,/)
99003 FORMAT (10X,'NCPA     ',I12,/,10X,'ICPA(IQ) ',9X,20I3,/)
99004 FORMAT (10X,'sites  IQ  JQ                     ',2I6,/,10X,
     &        'number of different TT-terms NTKTK',I12,I5,/)
C
      END
