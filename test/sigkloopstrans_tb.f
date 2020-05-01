C*==sigkloopstrans_tb.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGKLOOPSTRANS_TB(TAUDELTAK,TAUK,USSQX,VSSQX,WA,WB,
     &                             NKKR_TB,EXIKDQ_QTB)
C   ********************************************************************
C   *                                                                  *
C   *   transfer from G- to TAU conventions                            *
C   *                                                                  *
C   *   TAU[q,q'](k) = U[q] * TAU_Delta[q,q'](k) * V[q']               *
C   *                                                                  *
C   *                         * exp( k * (q-q') )                      *
C   *                                                                  *
C   *                                                                  *
C   *  NOTE: this is correct for q <> q'                               *
C   *        for q = q' the single site terms are not accounted for    *
C   *                                                                  *
C   *  The phase factor                                                *
C   *                                                                  *
C   *       EXIKDQ_QTB(IQTB,JQTB) = exp{-2*PI*k*(q[JQ] - q[IQ])}       *
C   *                             = exp( k * (q-q') )                  *
C   *                                                                  *
C   *  accounts that the TB-scheme uses real space G-functions         *
C   *  referring to the individual central site q. Doing the           *
C   *  Fourier transformation this is not accounted for. For site      *
C   *  diagonal TAU- or G-matrices this does not matters, but for      *
C   *  site off-diagonal matrices TAU[q,q'](k) this leads to a         *
C   *  phase factor that has to be added                               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_ANGMOM,ONLY:NKM,NXM
      USE MOD_SITES,ONLY:NQTB
      IMPLICIT NONE
C*--SIGKLOOPSTRANS_TB34
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKKR_TB
      COMPLEX*16 EXIKDQ_QTB(NQTB,NQTB),TAUDELTAK(NKKR_TB,NKKR_TB),
     &           TAUK(NKKR_TB,NKKR_TB),USSQX(NXM,NXM,NQTB),
     &           VSSQX(NXM,NXM,NQTB),WA(NXM,NXM),WB(NXM,NXM)
C
C Local variables
C
      INTEGER IQTB,IX0,IX1,IX2,JQTB,JX0,JX1,JX2
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------------------------------------------
C
      JX0 = -NKM
      DO JQTB = 1,NQTB
         JX0 = JX0 + NXM
         JX1 = JX0 + 1
         JX2 = JX0 + NXM
C
         IX0 = -NKM
         DO IQTB = 1,NQTB
            IX0 = IX0 + NXM
            IX1 = IX0 + 1
            IX2 = IX0 + NXM
C
            WA(1:NXM,1:NXM) = TAUDELTAK(IX1:IX2,JX1:JX2)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,WA,NXM,VSSQX(1,1,JQTB),
     &                 NXM,C0,WB,NXM)
C
            CALL ZGEMM('N','N',NXM,NXM,NXM,C1,USSQX(1,1,IQTB),NXM,WB,
     &                 NXM,C0,WA,NXM)
C
            WA(1:NXM,1:NXM) = WA(1:NXM,1:NXM)
C
            IF ( IQTB.NE.JQTB ) THEN
               TAUK(IX1:IX2,JX1:JX2) = WA(1:NXM,1:NXM)
     &                                 *EXIKDQ_QTB(IQTB,JQTB)
            ELSE
               TAUK(IX1:IX2,JX1:JX2) = WA(1:NXM,1:NXM)
            END IF
C
         END DO
      END DO
C
      END
C*==sigkloopschi.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGKLOOPSCHI(WK,WKSUM,TAUK,CHIZ,NKKR_TB,NZ12MAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SITES,ONLY:NQ_CHI,IQBOT_CHI,IQTOP_CHI
      IMPLICIT NONE
C*--SIGKLOOPSCHI107
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NKKR_TB,NZ12MAX
      REAL*8 WK,WKSUM
      COMPLEX*16 CHIZ(NKMMAX*NKMMAX*NQ_CHI,NKMMAX*NKMMAX*NQ_CHI,NZ12MAX)
     &           ,TAUK(NKKR_TB,NKKR_TB)
C
C Local variables
C
      INTEGER I1,I2,I3,I4,IQ,IQCHI,JQ,JQCHI,K1,K2,L1,L2,L3,L4,NKMSQ
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------------------------------------------
C
      NKMSQ = NKM*NKM
C
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IQCHI = IQ - IQBOT_CHI + 1
C
         DO JQ = IQBOT_CHI,IQTOP_CHI
            JQCHI = JQ - IQBOT_CHI + 1
C
            K1 = (IQCHI-1)*NKMSQ
            DO L1 = 1,NKM
               DO L4 = 1,NKM
                  K1 = K1 + 1
C
                  K2 = (JQCHI-1)*NKMSQ
                  DO L2 = 1,NKM
                     DO L3 = 1,NKM
                        K2 = K2 + 1
C
                        I1 = (IQCHI-1)*NKM + L1
                        I2 = (JQCHI-1)*NKM + L2
C
                        I3 = (JQCHI-1)*NKM + L3
                        I4 = (IQCHI-1)*NKM + L4
C
                        CHIZ(K1,K2,1) = CHIZ(K1,K2,1) + TAUK(I1,I2)
     &                                  *TAUK(I3,I4)*WK/WKSUM
                        CHIZ(K1,K2,2) = CHIZ(K1,K2,2) + TAUK(I1,I2)
     &                                  *DCONJG(TAUK(I4,I3))*WK/WKSUM
C
                     END DO
                  END DO
               END DO
            END DO
C
         END DO
      END DO
C
      END
