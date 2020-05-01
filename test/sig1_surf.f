C*==sig1_surf.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_SURF(MAQAB,MBQAB,MCQAB,MDQAB,SIG1Q,KEY,MAQBA,
     &                     MCQBA,MBQBA,MDQBA,NSPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *   calculate  TRACE jbar(mue,z2,z1)*mat*jbar(nue,z1,z2)           *
C   *   (including Vertex-corrections)                                 *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:CHIZ
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,LMAT3
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SOTSI,EESI,
     &    IRESPONSE_SOT,IRESPONSE_EDELSTEIN
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_CALCMODE,ONLY:PUBLIC_VERSION
      IMPLICIT NONE
C*--SIG1_SURF27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG1_SURF')
C
C Dummy arguments
C
      CHARACTER*1 KEY
      INTEGER NSPINPROJ
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           SIG1Q(3,3,NSPINPROJ,NQMAX,NQMAX)
C
C Local variables
C
      INTEGER I,IA_ERR,IQ,IQCHI,ISPINPROJ,ISPR,J,JQ,JQCHI,K1,K2,L1,L2,
     &        L3,L4,MUE,NKMSQ,NUE,Z1,Z2
      COMPLEX*16 MBQAB_L(:,:,:,:,:),MBQBA_L(:,:,:,:,:),
     &           MCQAB_L(:,:,:,:,:),MCQBA_L(:,:,:,:,:),MDQABX(:,:,:,:,:)
     &           ,MDQAB_L(:,:,:,:,:),MDQBAX(:,:,:,:,:),
     &           MDQBA_L(:,:,:,:,:),S11_1,S11_2,S12_1,S12_2,S21_1,S21_2,
     &           S22_1,S22_2,SIG1IIQ(:,:,:,:),SIG1IRQ(:,:,:,:),
     &           SUMX_1(3,3,2,2),SUMX_2(3,3,2,2),SUM_A(3,3),SUM_ALL(3,3)
     &           ,SUM_CB(3,3),SUM_D(3,3),S_A,S_CB,S_D,TMP_1(3,3),
     &           TMP_2(3,3)
      REAL*8 PRESI
      CHARACTER*40 SIUNITS
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIG1IIQ,SIG1IRQ,MDQBAX,MDQABX
      ALLOCATABLE MBQAB_L,MBQBA_L,MCQAB_L,MCQBA_L,MDQAB_L,MDQBA_L
C
      CALL TRACK_INFO(ROUTINE)
C
      NKMSQ = NKM*NKM
C
      ALLOCATE (SIG1IIQ(3,3,NQ,NQ),SIG1IRQ(3,3,NQ,NQ),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:sig1 -> SIG1IIQ'
C
      ALLOCATE (MDQBAX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQABX(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
      ALLOCATE (MBQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MBQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MCQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MCQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQAB_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
      ALLOCATE (MDQBA_L(NKMMAX,NKMMAX,3,NSPINPROJ,IQBOT_CHI:IQTOP_CHI))
C
      CALL CINIT(3*3*NQ*NQ,SIG1IIQ)
      CALL CINIT(3*3*NQ*NQ,SIG1IRQ)
C
      IF ( KEY.EQ.'N' ) WRITE (6,99005) ROUTINE,
     &                                  'without vertex-corrections'
      IF ( KEY.EQ.'V' ) WRITE (6,99005) ROUTINE,
     &                                  'including vertex-corrections'
C
C-----------------------------------------------------------------------
C  multiply the averaged MEs with proper LMAT to take into account
C  that CHIZ does not incude the LMATs
C-----------------------------------------------------------------------
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               MBQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(MBQAB(:,:,MUE,ISPINPROJ,IQ),LMAT3)
               MBQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MBQBA(:,:,MUE,ISPINPROJ,IQ))
C
               MCQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(MCQAB(:,:,MUE,ISPINPROJ,IQ),LMAT3)
               MCQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MCQBA(:,:,MUE,ISPINPROJ,IQ))
C
               MDQAB_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQAB(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
               MDQBA_L(:,:,MUE,ISPINPROJ,IQ)
     &            = MATMUL(LMAT3,MATMUL(MDQBA(:,:,MUE,ISPINPROJ,IQ),
     &            LMAT3))
            END DO
         END DO
      END DO
C
      DO ISPR = 1,NSPR
         ISPINPROJ = LIST_ISPR(ISPR)
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
C
               DO I = 1,NKM
                  DO J = 1,NKM
                     MDQBAX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQBA_L(J,I,MUE,ISPINPROJ,IQ))
C
                     MDQABX(I,J,MUE,ISPINPROJ,IQ)
     &                  = DCONJG(MDQAB_L(J,I,MUE,ISPINPROJ,IQ))
                  END DO
               END DO
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
C
         WRITE (6,99009)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(ISPINPROJ)
         WRITE (6,99009)
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
            IQCHI = IQ - IQBOT_CHI + 1
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            DO JQ = IQBOT_CHI,IQTOP_CHI
               JQCHI = JQ - IQBOT_CHI + 1
C
               WRITE (6,99001) IQ,JQ
C
               CALL CINIT(3*3*2*2,SUMX_1)
               CALL CINIT(3*3*2*2,SUMX_2)
C
               DO NUE = 1,3
C
                  DO MUE = 1,3
C
                     S11_1 = C0
                     S12_1 = C0
                     S21_1 = C0
                     S22_1 = C0
C
                     S11_2 = C0
                     S12_2 = C0
                     S21_2 = C0
                     S22_2 = C0
C
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
                                 S11_1 = S11_1 + 
     &                              MAQBA(L4,L1,MUE,ISPINPROJ,IQ)
     &                              *CHIZ(K1,K2,1)
     &                              *MAQAB(L2,L3,NUE,ISPINPROJ,JQ)
C
                                 S12_1 = S12_1 + 
     &                              MCQBA_L(L4,L1,MUE,ISPINPROJ,IQ)
     &                              *CHIZ(K1,K2,2)
     &                              *MBQAB_L(L2,L3,NUE,ISPINPROJ,JQ)
C
                                 S21_1 = S21_1 + 
     &                              DCONJG(MCQAB_L(L1,L4,NUE,ISPINPROJ,
     &                              IQ))*CHIZ(K1,K2,2)
     &                              *DCONJG(MBQBA_L(L3,L2,MUE,ISPINPROJ,
     &                              JQ))
C
                                 S22_1 = S22_1 + 
     &                              DCONJG(MDQAB_L(L1,L4,NUE,ISPINPROJ,
     &                              IQ))*CHIZ(K1,K2,1)
     &                              *DCONJG(MDQBA_L(L3,L2,MUE,ISPINPROJ,
     &                              JQ))
C
C
                                 S11_2 = S11_2 + 
     &                              MAQBA(L4,L1,MUE,ISPINPROJ,IQ)
     &                              *CHIZ(K1,K2,1)
     &                              *MAQAB(L2,L3,NUE,ISPINPROJ,JQ)
C
                                 S12_2 = S12_2 + 
     &                              MCQBA_L(L4,L1,MUE,ISPINPROJ,IQ)
     &                              *CHIZ(K1,K2,2)
     &                              *MBQAB_L(L2,L3,NUE,ISPINPROJ,JQ)
C
                                 S21_2 = S21_2 + 
     &                              DCONJG(MCQAB_L(L1,L4,NUE,ISPINPROJ,
     &                              IQ))*CHIZ(K1,K2,2)
     &                              *DCONJG(MBQBA_L(L3,L2,MUE,ISPINPROJ,
     &                              JQ))
C
                                 S22_2 = S22_2 + 
     &                              DCONJG(MDQAB_L(L1,L4,NUE,ISPINPROJ,
     &                              IQ))*CHIZ(K1,K2,1)
     &                              *DCONJG(MDQBA_L(L3,L2,MUE,ISPINPROJ,
     &                              JQ))
C
                              END DO
C
                           END DO
                        END DO
                     END DO
C
                     SUMX_1(MUE,NUE,1,1) = S11_1
                     SUMX_1(MUE,NUE,1,2) = S12_1
                     SUMX_1(MUE,NUE,2,1) = DCONJG(S21_1)
                     SUMX_1(MUE,NUE,2,2) = DCONJG(S22_1)
C
                     SUMX_2(MUE,NUE,1,1) = S11_2
                     SUMX_2(MUE,NUE,1,2) = S12_2
                     SUMX_2(MUE,NUE,2,1) = DCONJG(S21_2)
                     SUMX_2(MUE,NUE,2,2) = DCONJG(S22_2)
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
                           IF ( ABS(DIMAG(SUMX_1(MUE,NUE,Z1,Z2)))
     &                          .LT.1D-12 ) SUMX_1(MUE,NUE,Z1,Z2)
     &                          = DREAL(SUMX_1(MUE,NUE,Z1,Z2))
                           IF ( ABS(SUMX_1(MUE,NUE,Z1,Z2)).LT.1D-12 )
     &                          SUMX_1(MUE,NUE,Z1,Z2) = C0
C
                           IF ( ABS(DIMAG(SUMX_2(MUE,NUE,Z1,Z2)))
     &                          .LT.1D-12 ) SUMX_2(MUE,NUE,Z1,Z2)
     &                          = DREAL(SUMX_2(MUE,NUE,Z1,Z2))
                           IF ( ABS(SUMX_2(MUE,NUE,Z1,Z2)).LT.1D-12 )
     &                          SUMX_2(MUE,NUE,Z1,Z2) = C0
C
                        END DO
                     END DO
                  END DO
               END DO
C
C------------------------------------ write out k-integral parts (eq.74)
C
               IF ( IPRINT.GE.3 ) THEN
                  WRITE (6,*) 'SUMX_1'
                  DO Z1 = 1,2
                     DO Z2 = 1,2
                        WRITE (6,99006) Z1,Z2
                        DO MUE = 1,3
                           WRITE (6,99004)
     &                            (SUMX_1(MUE,NUE,Z1,Z2),NUE=1,3)
                        END DO
C
                     END DO
                  END DO
               END IF
C
C
               IF ( IPRINT.GE.3 ) THEN
                  WRITE (6,*) 'SUMX_2'
                  DO Z1 = 1,2
                     DO Z2 = 1,2
                        WRITE (6,99006) Z1,Z2
                        DO MUE = 1,3
                           WRITE (6,99004)
     &                            (SUMX_2(MUE,NUE,Z1,Z2),NUE=1,3)
                        END DO
C
                     END DO
                  END DO
               END IF
C
C----------------------------------------------------- j_m Im G j_n Im G
C----------------------------- symmetric part of tensor (Kubo-Greenwood)
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     SIG1IIQ(MUE,NUE,IQ,JQ)
     &                  = -0.25D0*(SUMX_1(MUE,NUE,1,1)
     &                  -SUMX_1(MUE,NUE,1,2)-SUMX_1(MUE,NUE,2,1)
     &                  +SUMX_1(MUE,NUE,2,2))
C
                     IF ( ABS(SIG1IIQ(MUE,NUE,IQ,JQ)).LT.1D-12 )
     &                    SIG1IIQ(MUE,NUE,IQ,JQ) = C0
                     IF ( ABS(DIMAG(SIG1IIQ(MUE,NUE,IQ,JQ))).LT.1D-12 )
     &                    SIG1IIQ(MUE,NUE,IQ,JQ)
     &                    = DREAL(SIG1IIQ(MUE,NUE,IQ,JQ))
                  END DO
               END DO
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
               IF ( .NOT.PUBLIC_VERSION ) THEN
C
                  IF ( MUE.LT.0 ) CALL SIG1_SURF_NONPUB(SUMX_1,SUMX_2,
     &                 SIG1IRQ,IQ,JQ)
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C------------------------------------------ antisymmetric part of tensor
C
                  DO MUE = 1,3
                     DO NUE = 1,3
                        TMP_1(MUE,NUE) = CI*(0.25D0/CI)
     &                     *(SUMX_1(MUE,NUE,1,1)+SUMX_1(MUE,NUE,1,2)
     &                     -SUMX_1(MUE,NUE,2,1)-SUMX_1(MUE,NUE,2,2))
C
                        TMP_2(MUE,NUE) = CI*(0.25D0/CI)
     &                     *(SUMX_2(MUE,NUE,1,1)+SUMX_2(MUE,NUE,1,2)
     &                     -SUMX_2(MUE,NUE,2,1)-SUMX_2(MUE,NUE,2,2))
C
                     END DO
                  END DO
                  DO MUE = 1,3
                     DO NUE = 1,3
C
                        SIG1IRQ(MUE,NUE,IQ,JQ) = TMP_1(MUE,NUE)
     &                     - TMP_2(NUE,MUE)
C
                        IF ( ABS(SIG1IRQ(MUE,NUE,IQ,JQ)).LT.1D-12 )
     &                       SIG1IRQ(MUE,NUE,IQ,JQ) = C0
                        IF ( ABS(DIMAG(SIG1IRQ(MUE,NUE,IQ,JQ)))
     &                       .LT.1D-12 ) SIG1IRQ(MUE,NUE,IQ,JQ)
     &                       = DREAL(SIG1IRQ(MUE,NUE,IQ,JQ))
                     END DO
                  END DO
C
C
C--------------- only 1/2 of SIG1IRQ --> Crepieux formula PRB 64, 014416
C--------------- also in the Bastin framework this factor has to be
C--------------- included
                  SIG1IRQ(:,:,IQ,JQ) = 0.5D0*SIG1IRQ(:,:,IQ,JQ)
C
               END IF
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
               DO MUE = 1,3
                  DO NUE = 1,3
                     SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)
     &                  = SIG1IIQ(MUE,NUE,IQ,JQ)
     &                  + SIG1IRQ(MUE,NUE,IQ,JQ)
                  END DO
               END DO
C
C-----------------------------------------------------------------------
               IF ( IPRINT.GE.3 ) THEN
C
                  WRITE (6,99002) '  j_m Im G j_n Im G'
                  DO MUE = 1,3
                     WRITE (6,99007) (DREAL(SIG1IIQ(MUE,NUE,IQ,JQ)),NUE=
     &                               1,3)
                  END DO
C
                  WRITE (6,99002) 
     &                        '  i ( j_m Im G j_n - j_n Im G j_m ) Re G'
                  DO MUE = 1,3
                     WRITE (6,99007) (DREAL(SIG1IRQ(MUE,NUE,IQ,JQ)),NUE=
     &                               1,3)
                  END DO
C
                  WRITE (6,99002) '  total'
                  DO MUE = 1,3
                     WRITE (6,99010) (SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ),NUE
     &                               =1,3)
                  END DO
C
C
                  WRITE (6,99008) TRIM(SIUNITS)
C
C
                  WRITE (6,99002) '  j_m Im G j_n Im G'
                  DO MUE = 1,3
                     WRITE (6,99007) (PRESI*DREAL(SIG1IIQ(MUE,NUE,IQ,JQ)
     &                               ),NUE=1,3)
                  END DO
C
                  WRITE (6,99002) 
     &                        '  i ( j_m Im G j_n - j_n Im G j_m ) Re G'
                  DO MUE = 1,3
                     WRITE (6,99007) (PRESI*DREAL(SIG1IRQ(MUE,NUE,IQ,JQ)
     &                               ),NUE=1,3)
                  END DO
C
                  WRITE (6,99002) '  total'
                  DO MUE = 1,3
                     WRITE (6,99010) (PRESI*SIG1Q(MUE,NUE,ISPINPROJ,IQ,
     &                               JQ),NUE=1,3)
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
                     IF ( ABS(DIMAG(SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)))
     &                    .GT.1D-5 ) WRITE (6,99003)
     &                    DIMAG(SIG1Q(MUE,NUE,ISPINPROJ,IQ,JQ)),MUE,NUE
                  END DO
               END DO
            END DO ! JQ
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C   q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2q2
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
                     S_A = C0
                     S_D = C0
                     S_CB = C0
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
                                 S_CB = S_CB + 
     &                                  MCQBA_L(L4,L1,MUE,ISPINPROJ,IQ)
     &                                  *CHIZ(K1,K2,2)
     &                                  *MBQAB_L(L2,L3,NUE,1,JQ)
C
                                 S_D = S_D + MDQBAX(L4,L1,MUE,ISPINPROJ,
     &                                 IQ)*CHIZ(K1,K2,1)
     &                                 *MDQABX(L2,L3,NUE,1,JQ)
C
                                 S_A = S_A + MAQBA(L4,L1,MUE,ISPINPROJ,
     &                                 IQ)*CHIZ(K1,K2,1)
     &                                 *MAQAB(L2,L3,NUE,1,JQ)
C
                              END DO
C
                           END DO
                        END DO
                     END DO
C
                     SUM_A(MUE,NUE) = 0.25D0*S_A
                     SUM_D(MUE,NUE) = 0.25D0*S_D
                     SUM_CB(MUE,NUE) = 0.25D0*S_CB
C
                  END DO
               END DO
C
               SUM_ALL = 2*SUM_CB - DCONJG(SUM_D) - SUM_A
C
C-----------------------------------------------------------------------
C               IF ( IPRINT.GE.3 ) THEN
               WRITE (6,*) 'SUMS'
               WRITE (6,*) 'SUM_ALL'
               WRITE (6,99010) ((SUM_ALL(MUE,NUE),NUE=1,3),MUE=1,3)
               WRITE (6,*) 'SUM_C'
               WRITE (6,99010) ((SUM_CB(MUE,NUE),NUE=1,3),MUE=1,3)
               WRITE (6,*) 'SUM_D'
               WRITE (6,99010) ((SUM_D(MUE,NUE),NUE=1,3),MUE=1,3)
               WRITE (6,*) 'SUM_A'
               WRITE (6,99010) ((SUM_A(MUE,NUE),NUE=1,3),MUE=1,3)
C               END IF
C
               SIG1Q(1:3,1:3,ISPINPROJ,IQ,JQ) = SUM_ALL(1:3,1:3)
C
C-----------------------------------------------------------------------
C
            END DO ! JQ
C
         END DO ! IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DEALLOCATE (SIG1IIQ,SIG1IRQ)
C
99001 FORMAT (/,10X,21('='),/,12X,'IQ, JQ = ',2I3,/,10X,21('='),/)
99002 FORMAT (/,10X,'sigma 1 ',5X,A,/)
99003 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2)
99004 FORMAT (3(F14.6,F12.6))
C99005 FORMAT (//,1X,79('*'),/,34X,'<SIG1_SURF>',/,1X,79('*'),//,10X,A,/)
99005 FORMAT (//,1X,79('*'),/,34X,'<',A9,'>',/,1X,79('*'),//,10X,A,/)
99006 FORMAT (/,10X,'(z1,z2) = ',2I2,/)
99007 FORMAT (10X,3F14.6)
99008 FORMAT (/,'    SI units [',A,']')
99009 FORMAT (/,40('*'),/,40('*'))
99010 FORMAT (3('(',F14.6,',',F12.6,')'))
      END
