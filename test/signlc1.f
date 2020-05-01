C*==signlc1.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGNLC1(CHIHATZ,MAQQAB,MBQQAB,MCQQAB,MDQQAB,SIG1Q,C,
     &                   IPRINT,NQNLCPA,NKM,NKMMAX,NQMAX,NZ12MAX,
     &                   NDIMCHI)
C   ********************************************************************
C   *                                                                  *
C   *   calculate  TRACE jbar(mue,z2,z1)*mat*jbar(nue,z1,z2)           *
C   *   (including Vertex-corrections)                                 *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *    THIS ROUTINE IS RESTRICTED TO NQMAX = 1 AT THE MOMENT !!!!!   *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_CALCMODE,ONLY:TASK
      IMPLICIT NONE
C*--SIGNLC119
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGNLC1')
C
C Dummy arguments
C
      CHARACTER*1 C
      INTEGER IPRINT,NDIMCHI,NKM,NKMMAX,NQMAX,NQNLCPA,NZ12MAX
      COMPLEX*16 CHIHATZ(NDIMCHI,NDIMCHI,NZ12MAX),
     &           MAQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MBQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MCQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           MDQQAB(NKMMAX,NKMMAX,3,NQNLCPA,NQNLCPA),
     &           SIG1Q(3,3,NQMAX,NQMAX)
C
C Local variables
C
      INTEGER IA_ERR,II,IQ1,IQ2,JJ,K1,K2,KK,L1,L2,L3,L4,LL,MUE,NKMSQ,
     &        NQKMSQ,NUE,Z1,Z2
      COMPLEX*16 S11,S12,S21,S22,SIG1IIQ(:,:,:,:),SIG1IRQ(:,:,:,:),
     &           SUMX(3,3,2,2),TMP(3,3)
      CHARACTER*5 TXTVAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIG1IIQ,SIG1IRQ
C
      ALLOCATE (SIG1IIQ(3,3,NQMAX,NQMAX))
      ALLOCATE (SIG1IRQ(3,3,NQMAX,NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate SIG1IIQ')
C
      CALL CINIT(3*3*NQMAX*NQMAX,SIG1Q)
C
      IF ( TASK.EQ.'GILBERT   ' ) THEN
         TXTVAR = 'alpha'
      ELSE
         TXTVAR = 'sigma'
      END IF
C
      NKMSQ = NKM*NKM
      NQKMSQ = NQNLCPA*NKMSQ
C
      IF ( C.EQ.'N' ) WRITE (6,99005) 'without vertex-corrections'
      IF ( C.EQ.'V' ) WRITE (6,99005) 'including vertex-corrections'
C
      IQ1 = 1
      IQ2 = 1
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      WRITE (6,99001) IQ1,IQ2
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
                  CALL CINIT(3*3*2*2,SUMX)
C
                  DO NUE = 1,3
                     DO MUE = 1,3
C
                        S11 = C0
                        S12 = C0
                        S21 = C0
                        S22 = C0
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
C                              sum up for eq. (74) or (38'):
C
                                    S11 = S11 + MAQQAB(L4,L1,MUE,LL,II)
     &                                 *CHIHATZ(K1,K2,1)
     &                                 *MAQQAB(L2,L3,NUE,JJ,KK)
C
                                    S12 = S12 + MCQQAB(L4,L1,MUE,LL,II)
     &                                 *CHIHATZ(K1,K2,2)
     &                                 *MBQQAB(L2,L3,NUE,JJ,KK)
C
                                    S21 = S21 + 
     &                                 DCONJG(MCQQAB(L1,L4,NUE,LL,II))
     &                                 *CHIHATZ(K1,K2,2)
     &                                 *DCONJG(MBQQAB(L3,L2,MUE,JJ,KK))
C
                                    S22 = S22 + 
     &                                 DCONJG(MDQQAB(L1,L4,NUE,LL,II))
     &                                 *CHIHATZ(K1,K2,1)
     &                                 *DCONJG(MDQQAB(L3,L2,MUE,JJ,KK))
C
                                 END DO
                              END DO
                           END DO
                        END DO
C
                        SUMX(MUE,NUE,1,1) = S11
                        SUMX(MUE,NUE,1,2) = S12
                        SUMX(MUE,NUE,2,1) = DCONJG(S21)
                        SUMX(MUE,NUE,2,2) = DCONJG(S22)
C
                     END DO
                  END DO
C
C -----    suppress small elements and complete missing elements of SUMX
C
                  DO Z1 = 1,2
                     DO Z2 = 1,2
                        DO MUE = 1,3
                           DO NUE = 1,3
                              IF ( ABS(SUMX(MUE,NUE,Z1,Z2)).LT.1D-12 )
     &                             SUMX(MUE,NUE,Z1,Z2) = C0
                              IF ( ABS(DIMAG(SUMX(MUE,NUE,Z1,Z2)))
     &                             .LT.1D-12 ) SUMX(MUE,NUE,Z1,Z2)
     &                             = DREAL(SUMX(MUE,NUE,Z1,Z2))
                           END DO
                        END DO
                     END DO
                  END DO
C
C------------------------------------ write out k-integral parts (eq.74)
C
                  IF ( IPRINT.GE.3 ) THEN
                     DO Z1 = 1,2
                        DO Z2 = 1,2
                           WRITE (6,99006) Z1,Z2
                           DO MUE = 1,3
                              WRITE (6,99004)
     &                               (SUMX(MUE,NUE,Z1,Z2),NUE=1,3)
                           END DO
C
                        END DO
                     END DO
                  END IF
C
C----------------------------------------------------- j_m Im G j_n Im G
C
                  DO MUE = 1,3
                     DO NUE = 1,3
                        SIG1IIQ(MUE,NUE,IQ1,IQ2)
     &                     = -0.25D0*(SUMX(MUE,NUE,1,1)
     &                     -SUMX(MUE,NUE,1,2)-SUMX(MUE,NUE,2,1)
     &                     +SUMX(MUE,NUE,2,2))
C
                        IF ( ABS(SIG1IIQ(MUE,NUE,IQ1,IQ2)).LT.1D-12 )
     &                       SIG1IIQ(MUE,NUE,IQ1,IQ2) = C0
                        IF ( ABS(DIMAG(SIG1IIQ(MUE,NUE,IQ1,IQ2)))
     &                       .LT.1D-12 ) SIG1IIQ(MUE,NUE,IQ1,IQ2)
     &                       = DREAL(SIG1IIQ(MUE,NUE,IQ1,IQ2))
                     END DO
                  END DO
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C
                  DO MUE = 1,3
                     DO NUE = 1,3
                        TMP(MUE,NUE) = CI*(0.25D0/CI)
     &                                 *(SUMX(MUE,NUE,1,1)+
     &                                 SUMX(MUE,NUE,1,2)
     &                                 -SUMX(MUE,NUE,2,1)
     &                                 -SUMX(MUE,NUE,2,2))
                     END DO
                  END DO
                  DO MUE = 1,3
                     DO NUE = 1,3
                        SIG1IRQ(MUE,NUE,IQ1,IQ2) = TMP(MUE,NUE)
     &                     - TMP(NUE,MUE)
C
                        IF ( ABS(SIG1IRQ(MUE,NUE,IQ1,IQ2)).LT.1D-12 )
     &                       SIG1IRQ(MUE,NUE,IQ1,IQ2) = C0
                        IF ( ABS(DIMAG(SIG1IRQ(MUE,NUE,IQ1,IQ2)))
     &                       .LT.1D-12 ) SIG1IRQ(MUE,NUE,IQ1,IQ2)
     &                       = DREAL(SIG1IRQ(MUE,NUE,IQ1,IQ2))
                     END DO
                  END DO
C
                  IF ( TASK.EQ.'GILBERT   ' ) THEN
C
                     DO MUE = 1,3
                        DO NUE = 1,3
                           SIG1Q(MUE,NUE,IQ1,IQ2)
     &                        = SIG1Q(MUE,NUE,IQ1,IQ2)
     &                        + SIG1IIQ(MUE,NUE,IQ1,IQ2)
                           SIG1IRQ(MUE,NUE,IQ1,IQ2) = 0D0
                        END DO
                     END DO
C
                  ELSE
C
                     DO MUE = 1,3
C
                        DO NUE = 1,3
                           SIG1Q(MUE,NUE,IQ1,IQ2)
     &                        = SIG1Q(MUE,NUE,IQ1,IQ2)
     &                        + SIG1IIQ(MUE,NUE,IQ1,IQ2)
     &                        + SIG1IRQ(MUE,NUE,IQ1,IQ2)
                        END DO
                     END DO
C
                  END IF
C
               END DO
            END DO
         END DO
      END DO
C                                loops over cluster sites II, JJ, KK, LL
C***********************************************************************
C
C-----------------------------------------------------------------------
      IF ( IPRINT.GE.2 ) THEN
C
         IF ( TASK.NE.'GILBERT   ' ) THEN
C
            WRITE (6,99002) TXTVAR,'  j_m Im G j_n Im G'
            DO MUE = 1,3
               WRITE (6,99007) (DREAL(SIG1IIQ(MUE,NUE,IQ1,IQ2)),NUE=1,3)
            END DO
C
            WRITE (6,99002) TXTVAR,
     &                      '  i ( j_m Im G j_n - j_n Im G j_m ) Re G'
            DO MUE = 1,3
               WRITE (6,99007) (DREAL(SIG1IRQ(MUE,NUE,IQ1,IQ2)),NUE=1,3)
            END DO
C
         END IF
C
         WRITE (6,99002) TXTVAR,'  total'
         DO MUE = 1,3
            WRITE (6,99007) (DREAL(SIG1Q(MUE,NUE,IQ1,IQ2)),NUE=1,3)
         END DO
C
         WRITE (6,*)
      END IF
C-----------------------------------------------------------------------
C
C
C--------------------------------- check if imaginary part of SIG1Q is 0
C
      DO MUE = 1,3
         DO NUE = 1,3
            IF ( ABS(DIMAG(SIG1Q(MUE,NUE,IQ1,IQ2))).GT.1D-5 )
     &           WRITE (6,99003) TXTVAR,DIMAG(SIG1Q(MUE,NUE,IQ1,IQ2)),
     &                           MUE,NUE
         END DO
      END DO
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (SIG1IIQ,SIG1IRQ)
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ1, IQ2 = ',2I3,/,10X,21('='),/)
99002 FORMAT (/,10X,A,' 1 ',5X,A,/)
99003 FORMAT (' WARNING!! Im(',A,') =',e13.5,' for mue,nue=',2I2)
C99004 FORMAT (3(F14.6,F12.6))
99004 FORMAT (3(E16.6,E14.6))
99005 FORMAT (//,1X,79('*'),/,37X,'<SIG1>',/,1X,79('*'),//,10X,A,/)
99006 FORMAT (/,10X,'(z1,z2) = ',2I2,/)
99007 FORMAT (10X,3E15.6)
      END
