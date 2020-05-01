C*==sig1_fsea.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_FSEA(WINTEG,MAQAB,MAQBA,MDQAB,MDQBA,SIG1Q,KEY,
     &                     NSPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the Fermi-sea contribution of conductivity tensor    *
C   *  for the SIGMA1 term according to the Bastin formula             *
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
C   *   calculate  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)           *
C   *                                                                  *
C   *   where the auxilary quantity                                    *
C   *                                                                  *
C   *           CHI(K1,K2,IZ12) = Int d3k Tau(k,E_a) * Tau(k,E_b)      *
C   *                                                                  *
C   *   may contain the vertex corrections                             *
C   *                                                                  *
C   *  NOTE: IZ12 = 1    energies   Tau(k,E_a) * Tau(k,E_b)            *
C   *               2    energies   Tau(k,E_a) * Tau(k,E_b^*)          *
C   *        here we need only  IZ12 = 1                               *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   * the prefactor 1/4 is included here in the S-terms (Bastin)       *
C   * constants are added in <SIG_SUM> when converting to SI units     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:CHIZ
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,LMAT3
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SOTSI,EESI,
     &    IRESPONSE_SOT,IRESPONSE_EDELSTEIN
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--SIG1_FSEA49
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG1_FSEA')
C
C Dummy arguments
C
      CHARACTER*1 KEY
      INTEGER NSPINPROJ
      COMPLEX*16 WINTEG
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           SIG1Q(3,3,NSPINPROJ,NQMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 DELTA,D_SIG1_FSEA_DE(3,3),MDQABX(:,:,:,:,:),
     &           MDQBAX(:,:,:,:,:),SIG1_MTMTQ(:,:,:,:),SUM_A(3,3),
     &           SUM_A1(3,3),SUM_A2(3,3),SUM_D(3,3),SUM_D1(3,3),
     &           SUM_D2(3,3),S_A1,S_A2,S_D1,S_D2
      INTEGER I,IA_ERR,IQ,IQCHI,ISPINPROJ,ISPR,J,JQ,JQCHI,K1,K2,L1,L2,
     &        L3,L4,MUE,NKMSQ,NUE
      REAL*8 PRESI
      CHARACTER*40 SIUNITS
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIG1_MTMTQ,MDQABX,MDQBAX
C
      CALL TRACK_INFO(ROUTINE)
C
      NKMSQ = NKM*NKM
C
      ALLOCATE (MDQABX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQBAX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (SIG1_MTMTQ(3,3,NQ,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:sig1_fsea -> SIG1_MTMTQ'
C
      CALL CINIT(3*3*NQ*NQ,SIG1_MTMTQ)
C
      IF ( KEY.EQ.'N' ) WRITE (6,99005) 'without vertex-corrections'
      IF ( KEY.EQ.'V' ) WRITE (6,99005) 'including vertex-corrections'
C
C-----------------------------------------------------------------------
C    supply auxilary matrix elements for type D corresponding to (-,-)
C-----------------------------------------------------------------------
C
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               DO I = 1,NKM
                  DO J = 1,NKM
                     MDQABX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQAB(J,I,MUE,ISPINPROJ,IQ))
C
                     MDQBAX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQBA(J,I,MUE,ISPINPROJ,IQ))
                  END DO
               END DO
C
C-------------- multiply the averaged MEs with LMAT to take into account
C-------------- that CHIZ(+,-)  does not include the LMATs: L M L
C
               MDQBAX(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQBAX(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
C
               MDQABX(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQABX(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
C
            END DO
         END DO
      END DO
C
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
C
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
C
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
C---------- torkance prefactor for SI output (C * m)
            PRESI = SOTSI
            SIUNITS = '10**-30 C*m'
         ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
C---------- Edelstein prefactor for SI output (m/V)
            PRESI = EESI
            SIUNITS = 'm/V'
         ELSE
C---------- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
            PRESI = 1D-8*CONSI
            SIUNITS = '1/(muOhm*cm)'
         END IF
C
         WRITE (6,99007)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(ISPINPROJ)
         WRITE (6,99007)
C
         D_SIG1_FSEA_DE(:,:) = C0
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
            IQCHI = IQ - IQBOT_CHI + 1
C
            DO JQ = IQBOT_CHI,IQTOP_CHI
               JQCHI = JQ - IQBOT_CHI + 1
C
               WRITE (6,99001) IQ,JQ
C
               DO NUE = 1,3
C
                  DO MUE = 1,3
C
                     S_A1 = C0
                     S_A2 = C0
                     S_D1 = C0
                     S_D2 = C0
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
                                 S_A1 = S_A1 + 
     &                                  MAQBA(L4,L1,MUE,ISPINPROJ,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MAQAB(L2,L3,NUE,1,JQ)
C
                                 S_A2 = S_A2 + MAQBA(L4,L1,NUE,1,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MAQAB(L2,L3,MUE,ISPINPROJ,JQ)
C
                                 S_D1 = S_D1 + 
     &                                  MDQABX(L4,L1,MUE,ISPINPROJ,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MDQBAX(L2,L3,NUE,1,JQ)
C
                                 S_D2 = S_D2 + MDQABX(L4,L1,NUE,1,IQ)
     &                                  *CHIZ(K1,K2,1)
     &                                  *MDQBAX(L2,L3,MUE,ISPINPROJ,JQ)
C
                              END DO
C
                           END DO
                        END DO
                     END DO
C
C ---------------------------------------------- suppress small elements
                     IF ( ABS(S_A1).LT.1D-12 ) S_A1 = C0
                     IF ( ABS(S_A2).LT.1D-12 ) S_A2 = C0
                     IF ( ABS(S_D1).LT.1D-12 ) S_D1 = C0
                     IF ( ABS(S_D2).LT.1D-12 ) S_D2 = C0
C ----------------------------------------------------------------------
C
                     SUM_A1(MUE,NUE) = 0.25D0*S_A1
                     SUM_A2(MUE,NUE) = 0.25D0*S_A2
C
                     SUM_D1(MUE,NUE) = 0.25D0*DCONJG(S_D1)
                     SUM_D2(MUE,NUE) = 0.25D0*DCONJG(S_D2)
C
                  END DO
               END DO
C
C-----------------------------------------------------------------------
               IF ( IPRINT.GE.3 ) THEN
                  WRITE (6,*) 'SUM_A1'
                  WRITE (6,99004) ((SUM_A1(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM_A2'
                  WRITE (6,99004) ((SUM_A2(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM_D1'
                  WRITE (6,99004) ((SUM_D1(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,*) 'SUM_D2'
                  WRITE (6,99004) ((SUM_D2(MUE,NUE),NUE=1,3),MUE=1,3)
               END IF
C-----------------------------------------------------------------------
C
C w(E_a,+)   * [ j_m G(E_a,+) j_n G(E_b,+) - j_n G(E_b,+) j_m G(E_a,+) ]
C
               SUM_A(1:3,1:3) = SUM_A1(1:3,1:3) - SUM_A2(1:3,1:3)
C
C w(E_a,-)^* * [ j_m G(E_a,-) j_n G(E_b,-) - j_n G(E_b,-) j_m G(E_a,-) ]
C
               SUM_D(1:3,1:3) = SUM_D1(1:3,1:3) - SUM_D2(1:3,1:3)
C
C----------------------------------------------------------------- CHECK
C
               DO MUE = 1,3
                  IF ( ABS(SUM_A(MUE,MUE)).GT.1D-6 ) WRITE (6,99009)
     &                 'A',IQ,MUE,0,SUM_A(MUE,MUE)
                  DO NUE = 1,3
                     DELTA = SUM_A(MUE,NUE) - DCONJG(SUM_D(MUE,NUE))
                     IF ( ABS(DELTA).GT.1D-6 ) WRITE (6,99009) 'B',IQ,
     &                    MUE,NUE,DELTA
                     DELTA = SUM_A(MUE,NUE) + SUM_D(MUE,NUE)
     &                       + SUM_A(NUE,MUE) + SUM_D(NUE,MUE)
                     IF ( ABS(DELTA).GT.1D-6 ) WRITE (6,99009) 'C',IQ,
     &                    MUE,NUE,DELTA
                  END DO
               END DO
C
C------------------------------------------------------------- integrate
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)
     &                  = SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)
     &                  + WINTEG*SUM_A(MUE,NUE) + DCONJG(WINTEG)
     &                  *SUM_D(MUE,NUE)
C
                  END DO
               END DO
C
               D_SIG1_FSEA_DE(1:3,1:3) = D_SIG1_FSEA_DE(1:3,1:3)
     &            + WINTEG*SUM_A(1:3,1:3) + DCONJG(WINTEG)
     &            *SUM_D(1:3,1:3)
C
C-----------------------------------------------------------------------
               IF ( IPRINT.GE.3 ) THEN
C
                  WRITE (6,99002) '  sum A'
                  WRITE (6,99008) ((SUM_A(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,99002) '  sum D'
                  WRITE (6,99008) ((SUM_D(MUE,NUE),NUE=1,3),MUE=1,3)
                  WRITE (6,99002) '  sum total A + D'
                  WRITE (6,99008) ((SUM_A(MUE,NUE)+SUM_D(MUE,NUE),NUE=1,
     &                            3),MUE=1,3)
C
                  WRITE (6,99006) TRIM(SIUNITS)
                  WRITE (6,99002) '  total'
                  DO MUE = 1,3
                     WRITE (6,99008) (PRESI*(SUM_A(MUE,NUE)+SUM_D(MUE,
     &                               NUE)),NUE=1,3)
                  END DO
C
                  WRITE (6,*)
               END IF
C-----------------------------------------------------------------------
C
C--------------------------------- check if imaginary part of SIG1Q is 0
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     IF ( ABS(DIMAG(SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)))
     &                    .GT.1D-5 ) WRITE (6,99003)
     &                    DIMAG(SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)),MUE,NUE
                  END DO
               END DO
            END DO ! JQ
C   q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2
C
C
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
         WRITE (6,99010)
         CALL PR_COND_TENSOR(D_SIG1_FSEA_DE,1D0)
C
      END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DEALLOCATE (SIG1_MTMTQ)
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ, JQ = ',2I3,/,10X,21('='),/)
99002 FORMAT (/,10X,'sigma 1 ',5X,A,/)
99003 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2)
99004 FORMAT (3(2X,2E12.5))
99005 FORMAT (//,1X,79('*'),/,34X,'<SIG1_FSEA>',/,1X,79('*'),//,10X,A,/)
99006 FORMAT (/,'    SI units [',A,']')
99007 FORMAT (/,40('*'),/,40('*'))
C99008 FORMAT (3('(',F14.6,',',F12.6,')'))
99008 FORMAT (3(2X,2E12.5))
99009 FORMAT (/,'#### <SIG1_FSEA> TEST ',A,' failed  IQ:',I3,3X,2I2,
     &        2F20.12)
99010 FORMAT (/,5X,'site-off diag. contr. to SIG1, including integ. ',
     &        'weight (Fermi-sea contribution)',/)
      END
