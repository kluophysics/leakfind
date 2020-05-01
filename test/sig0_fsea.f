C*==sig0_fsea.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG0_FSEA(WINTEG,S10AQAB,S10DQAB,S2AQAB,S2DQAB,S3AQAB,
     &                     S3DQAB,S4AQAB,S4DQAB,S10AQBA,S10DQBA,S2AQBA,
     &                     S2DQBA,S3AQBA,S3DQBA,S4AQBA,S4DQBA,NSPINPROJ,
     &                     SIG0Q_FSEA)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the Fermi-sea contribution of conductivity tensor    *
C   *  for the SIGMA0 term according to the Bastin formula             *
C   *                                                                  *
C   *         J_m dG(-)/dE j_n G(-)  -  J_m G(-) j_n dG(-)/dE          *
C   *                                                                  *
C   *       + J_m G(+) j_n dG(+)/dE  -  J_m dG(+)/dE j_n G(+)          *
C   *                                                                  *
C   * - the derivative d / dE is realized by pairs of energies (A,B)   *
C   *   and also the weight for taking the derivative                  *
C   * - the weight WINTEG includes the Gaussian weight for integration *
C   * - energy pairs below (-) and above (+) the real energy axis      *
C   *   are represented by the matrix elements of type  D  and  A      *
C   *                                                                  *
C   * the prefactor 1/4 is included here in the S-terms (Bastin)       *
C   * constants are added in <SIG_SUM> when converting to SI units     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_SITES,ONLY:NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR
      IMPLICIT NONE
C*--SIG0_FSEA30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG0_FSEA')
C
C Dummy arguments
C
      INTEGER NSPINPROJ
      COMPLEX*16 WINTEG
      COMPLEX*16 S10AQAB(3,3,NSPINPROJ,NQMAX),
     &           S10AQBA(3,3,NSPINPROJ,NQMAX),
     &           S10DQAB(3,3,NSPINPROJ,NQMAX),
     &           S10DQBA(3,3,NSPINPROJ,NQMAX),
     &           S2AQAB(3,3,NSPINPROJ,NQMAX),S2AQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S2DQAB(3,3,NSPINPROJ,NQMAX),
     &           S2DQBA(3,3,NSPINPROJ,NQMAX),S3AQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S3AQBA(3,3,NSPINPROJ,NQMAX),
     &           S3DQAB(3,3,NSPINPROJ,NQMAX),S3DQBA(3,3,NSPINPROJ,NQMAX)
     &           ,S4AQAB(3,3,NSPINPROJ,NQMAX),
     &           S4AQBA(3,3,NSPINPROJ,NQMAX),S4DQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4DQBA(3,3,NSPINPROJ,NQMAX),
     &           SIG0Q_FSEA(3,3,NSPINPROJ,NQMAX)
C
C Local variables
C
      COMPLEX*16 DELTA,D_SIG0_FSEA_DE(3,3),SUM_A(3,3),SUM_AD(3,3),
     &           SUM_D(3,3)
      INTEGER IQ,ISP,ISPR,MUE,NUE
C
C*** End of declarations rewritten by SPAG
C
      CALL TRACK_INFO(ROUTINE)
C
      WRITE (6,99001)
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO ISPR = 1,NSPR
         ISP = LIST_ISPR(ISPR)
         WRITE (6,99004)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(ISP)
         WRITE (6,99004)
C
         D_SIG0_FSEA_DE(:,:) = C0
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         LOOP_IQ:DO IQ = IQBOT_CHI,IQTOP_CHI
C
            WRITE (6,99002) IQ
C
C   *         J_m dG(-)/dE j_n G(-)  -  J_m G(-) j_n dG(-)/dE          *
C
C
            SUM_D(1:3,1:3) = S10DQBA(1:3,1:3,ISP,IQ)
     &                       + S2DQBA(1:3,1:3,ISP,IQ)
     &                       + S3DQBA(1:3,1:3,ISP,IQ)
     &                       + S4DQBA(1:3,1:3,ISP,IQ)
     &                       - S10DQAB(1:3,1:3,ISP,IQ)
     &                       - S2DQAB(1:3,1:3,ISP,IQ)
     &                       - S3DQAB(1:3,1:3,ISP,IQ)
     &                       - S4DQAB(1:3,1:3,ISP,IQ)
C
            SUM_D(1:3,1:3) = 0.25D0*SUM_D(1:3,1:3)
C
C
C   *       + J_m G(+) j_n dG(+)/dE  -  J_m dG(+)/dE j_n G(+)          *
C
C
            SUM_A(1:3,1:3) = S10AQAB(1:3,1:3,ISP,IQ)
     &                       + S2AQAB(1:3,1:3,ISP,IQ)
     &                       + S3AQAB(1:3,1:3,ISP,IQ)
     &                       + S4AQAB(1:3,1:3,ISP,IQ)
     &                       - S10AQBA(1:3,1:3,ISP,IQ)
     &                       - S2AQBA(1:3,1:3,ISP,IQ)
     &                       - S3AQBA(1:3,1:3,ISP,IQ)
     &                       - S4AQBA(1:3,1:3,ISP,IQ)
C
            SUM_A(1:3,1:3) = 0.25D0*SUM_A(1:3,1:3)
C
C----------------------------------------------------------------- CHECK
C
            DO MUE = 1,3
               IF ( ABS(SUM_A(MUE,MUE)).GT.1D-6 ) WRITE (6,99006) 'A',
     &              IQ,MUE,MUE,SUM_A(MUE,MUE)
               DO NUE = 1,3
                  DELTA = SUM_A(MUE,NUE) - DCONJG(SUM_D(MUE,NUE))
                  IF ( ABS(DELTA).GT.1D-6 ) WRITE (6,99006) 'B',IQ,MUE,
     &                 NUE,DELTA
                  DELTA = SUM_A(MUE,NUE) + SUM_D(MUE,NUE)
     &                    + SUM_A(NUE,MUE) + SUM_D(NUE,MUE)
                  IF ( ABS(DELTA).GT.1D-6 ) WRITE (6,99006) 'C',IQ,MUE,
     &                 NUE,DELTA
               END DO
            END DO
C
C --------------- site IQ resolved E-dependent contribution to SIG0_FSEA
C
            SUM_AD(1:3,1:3) = DCONJG(WINTEG)*SUM_D(1:3,1:3)
     &                        + WINTEG*SUM_A(1:3,1:3)
C
            D_SIG0_FSEA_DE(1:3,1:3) = D_SIG0_FSEA_DE(1:3,1:3)
     &                                + SUM_AD(1:3,1:3)
C
            SIG0Q_FSEA(1:3,1:3,ISP,IQ) = SIG0Q_FSEA(1:3,1:3,ISP,IQ)
     &         + SUM_AD(1:3,1:3)
C
C--------------------------------- check if imaginary part of SIG0Q is 0
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(DIMAG(SIG0Q_FSEA(MUE,NUE,ISP,IQ))).GT.1D-5 )
     &                 WRITE (6,99003) DIMAG(SIG0Q_FSEA(MUE,NUE,ISP,IQ))
     &                                 ,MUE,NUE,IQ
               END DO
            END DO
C-----------------------------------------------------------------------
         END DO LOOP_IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
         WRITE (6,99005)
C         WRITE (6,99006) ((D_SIG0_FSEA_DE(MUE,NUE),NUE=1,3),MUE=1,3)
         CALL PR_COND_TENSOR(D_SIG0_FSEA_DE,1D0)
C
      END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
99001 FORMAT (//,1X,79('*'),/,34X,'<SIG0_FSEA>',/,1X,79('*'),/)
99002 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99003 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99004 FORMAT (/,37('*'),' SIG0 ',37('*'))
99005 FORMAT (/,10X,'site-diagonal term SIGMA_0   ',
     &        '(Fermi-sea contribution)',/)
C99006 FORMAT (3(2X,2E12.5))
99006 FORMAT (/,'#### <SIG0_FSEA> TEST ',A,' failed  IQ:',I3,3X,2I2,
     &        2F20.12)
      END
