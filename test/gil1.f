C*==gil1.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GIL1(MAQAB,MBQAB,MCQOAB,MDQOAB,ALF1QQ,ALF1QOQ,KEY,
     &                MAQOBA,MBQBA,MCQOBA,MDQBA)
C   ********************************************************************
C   *                                                                  *
C   *   calculate  TRACE jbar(mue,z2,z1)*mat*jbar(nue,z1,z2)           *
C   *   (including Vertex-corrections)                                 *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:CHIZ
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SITES,ONLY:NQMAX,NOMAX,IQBOT_CHI,IQTOP_CHI,NOQ,ITOQ
      USE MOD_TYPES,ONLY:CONC
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--GIL121
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GIL1')
C
C Dummy arguments
C
      CHARACTER*1 KEY
      COMPLEX*16 ALF1QOQ(3,3,NQMAX,NOMAX,NQMAX),ALF1QQ(3,3,NQMAX,NQMAX),
     &           MAQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MAQOBA(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NQMAX),
     &           MCQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MCQOBA(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NQMAX),
     &           MDQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX)
C
C Local variables
C
      COMPLEX*16 ALF1II(:,:),S11,S12,S21,S22,SUMX(3,3,2,2)
      INTEGER IO,IQ,IQCHI,IT,JQ,JQCHI,K1,K2,L1,L2,L3,L4,MUE,NKMSQ,NUE,
     &        Z1,Z2
      CHARACTER*5 TXTVAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ALF1II
C
      CALL TRACK_INFO(ROUTINE)
C
      TXTVAR = 'alpha'
C
      NKMSQ = NKM*NKM
C
      ALLOCATE (ALF1II(3,3))
C
      ALF1QQ(:,:,:,:) = C0
C
      IF ( KEY.EQ.'N' ) WRITE (6,99005) 'without vertex-corrections'
      IF ( KEY.EQ.'V' ) WRITE (6,99005) 'including vertex-corrections'
C
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <GIL1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IQCHI = IQ - IQBOT_CHI + 1
         DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            DO JQ = IQBOT_CHI,IQTOP_CHI
               JQCHI = JQ - IQBOT_CHI + 1
C
               IF ( IPRINT.GT.0 ) WRITE (6,99001) IQ,JQ,IO
C
               SUMX(:,:,:,:) = C0
C
               DO NUE = 1,3
C
                  DO MUE = 1,3
C
                     S11 = C0
                     S12 = C0
                     S21 = C0
                     S22 = C0
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
C                              sum up for eq. (74) or (38'):
C
                                 S11 = S11 + MAQOBA(L4,L1,MUE,IQ,IO)
     &                                 *CHIZ(K1,K2,1)
     &                                 *MAQAB(L2,L3,NUE,JQ)
C
                                 S12 = S12 + MCQOBA(L4,L1,MUE,IQ,IO)
     &                                 *CHIZ(K1,K2,2)
     &                                 *MBQAB(L2,L3,NUE,JQ)
C
                                 S21 = S21 + 
     &                                 DCONJG(MCQOAB(L1,L4,NUE,IQ,IO))
     &                                 *CHIZ(K1,K2,2)
     &                                 *DCONJG(MBQBA(L3,L2,MUE,JQ))
C
                                 S22 = S22 + 
     &                                 DCONJG(MDQOAB(L1,L4,NUE,IQ,IO))
     &                                 *CHIZ(K1,K2,1)
     &                                 *DCONJG(MDQBA(L3,L2,MUE,JQ))
C
                              END DO
C
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
                           IF ( ABS(DIMAG(SUMX(MUE,NUE,Z1,Z2)))
     &                          .LT.1D-12 ) SUMX(MUE,NUE,Z1,Z2)
     &                          = DREAL(SUMX(MUE,NUE,Z1,Z2))
                           IF ( ABS(SUMX(MUE,NUE,Z1,Z2)).LT.1D-12 )
     &                          SUMX(MUE,NUE,Z1,Z2) = C0
                        END DO
                     END DO
                  END DO
               END DO
C
C------------------------------------ write out k-integral parts (eq.74)
C
               IF ( IPRINT.GT.0 ) THEN
                  DO Z1 = 1,2
                     DO Z2 = 1,2
                        WRITE (6,99006) Z1,Z2
                        DO MUE = 1,3
                           WRITE (6,99004) (SUMX(MUE,NUE,Z1,Z2),NUE=1,3)
                        END DO
                     END DO
                  END DO
               END IF
C
C----------------------------------------------------- j_m Im G j_n Im G
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     ALF1II(MUE,NUE) = -0.25D0*(SUMX(MUE,NUE,1,1)-SUMX(
     &                                 MUE,NUE,1,2)-SUMX(MUE,NUE,2,1)
     &                                 +SUMX(MUE,NUE,2,2))
C
                     IF ( ABS(ALF1II(MUE,NUE)).LT.1D-12 )
     &                    ALF1II(MUE,NUE) = C0
                     IF ( ABS(DIMAG(ALF1II(MUE,NUE))).LT.1D-12 )
     &                    ALF1II(MUE,NUE) = DREAL(ALF1II(MUE,NUE))
                  END DO
               END DO
C
               ALF1QOQ(:,:,IQ,IO,JQ) = ALF1II(:,:)
C
               ALF1QQ(:,:,IQ,JQ) = ALF1QQ(:,:,IQ,JQ) + CONC(IT)
     &                             *ALF1II(:,:)
C
C-----------------------------------------------------------------------
               IF ( IPRINT.GT.0 ) THEN
C
                  WRITE (6,99002) TXTVAR,'  total'
                  DO MUE = 1,3
                     WRITE (6,99007) (DREAL(ALF1QQ(MUE,NUE,IQ,JQ)),NUE=1
     &                               ,3)
                  END DO
C
                  WRITE (6,*)
               END IF
C-----------------------------------------------------------------------
C
C
C-------------------------------- check if imaginary part of ALF1QQ is 0
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     IF ( ABS(DIMAG(ALF1QQ(MUE,NUE,IQ,JQ))).GT.1D-5 )
     &                    WRITE (6,99003) TXTVAR,
     &                           DIMAG(ALF1QQ(MUE,NUE,IQ,JQ)),MUE,NUE
                  END DO
               END DO
C
            END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         END DO
      END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (ALF1II)
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ, JQ = ',2I3,3X,'IO =',I4,/,10X,
     &        21('='),/)
99002 FORMAT (/,10X,A,' 1 ',5X,A,/)
99003 FORMAT (' WARNING!! Im(',A,') =',e13.5,' for mue,nue=',2I2)
99004 FORMAT (3(F14.6,F12.6))
99005 FORMAT (//,1X,79('*'),/,37X,'<ALF1>',/,1X,79('*'),//,10X,A,/)
99006 FORMAT (/,10X,'(z1,z2) = ',2I2,/)
99007 FORMAT (10X,3F14.6)
      END
