C*==linresp_magnet_perturb.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_PERTURB(IFIL_LHSA,IFIL_RHSB,ERYD)
C   ********************************************************************
C   *                                                                  *
C   *   driver routine to                                              *
C   *   set up the perturbation matrix functions for atom type IT      *
C   *                                                                  *
C   *    HA_ZZ(r) = Int_{0}^{r}       Z^x(r',E_a)  Delta H  Z(r',E_b)  *
C   *    HB_JZ(r) = Int_{r}^{r_crit}  J^x(r',E_a)  Delta H  Z(r',E_b)  *
C   *    HC_ZJ(r) = Int_{r}^{r_crit}  Z^x(r',E_a)  Delta H  J(r',E_b)  *
C   *    HD_JJ(r) = Int_{r}^{r_crit}  J^x(r',E_a)  Delta H  J(r',E_b)  *
C   *                                                                  *
C   *   for a scalar angular dependent perturbation                    *
C   *                                                                  *
C   *      Delta H(->r) = Sum_lm  H_lm(r) * Y_lm(^r)                   *
C   *                                                                  *
C   *   with non-vanishing terms indicated by the flag  KDVLMT         *
C   *                                                                  *
C   *   the endings A (B) for the wave functions Z and J indicate the  *
C   *   energies E_A(B) for  Delta G(E_a,E_b) =  G(E_a) Delta H G(E_b) *
C   *                                                                  *
C   *   for IREL = 3 (fully relativistic) it is implied that the       *
C   *   LHS solutions Z^x(r,E_a) and J^x(r,E_a) are read from file     *
C   *   IFIL_LHSA as these may differ from the RHS solutions for E_b   *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:HZ_PERT_LMTP,HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,
     &    HDZ_JJ_T
      USE MOD_ANGMOM,ONLY:NKMMAX,NCPLWF
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:NCPLWFMAX,ITBOT,ITTOP,IMT,IKMCPLWF
      IMPLICIT NONE
C*--LINRESP_MAGNET_PERTURB34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_PERTURB')
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IFIL_LHSA,IFIL_RHSB
C
C Local variables
C
      INTEGER IA_ERR,IM,IRCRIT,IT
      COMPLEX*16 JFA(:,:,:),JFB(:,:,:),JGA(:,:,:),JGB(:,:,:),ZFA(:,:,:),
     &           ZFB(:,:,:),ZGA(:,:,:),ZGB(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JGA,JGB,ZGA,ZGB,JFA,JFB,ZFA,ZFB
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( IREL.EQ.3 ) THEN
         ALLOCATE (JFA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFB(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFB(NRMAX,NCPLWFMAX,NKMMAX))
      END IF
C
      ALLOCATE (JGA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGB(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ZGB')
C
C-----------------------------------------------------------------------
C     initialize the r-dependent perturbation matrix elements
C-----------------------------------------------------------------------
C
      HAZ_ZZ_T(:,:,:,:,:) = 0D0
      HCZ_ZJ_T(:,:,:,:,:) = 0D0
      HBZ_JZ_T(:,:,:,:,:) = 0D0
      HDZ_JJ_T(:,:,:,:,:) = 0D0
C
C-----------------------------------------------------------------------
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRCRIT = JRCRI(IM)
C
C ------------------------------------------ read in wavefunctions for A
C
         CALL WAVFUN_READ_REL(IFIL_LHSA,IT,1,ZGA,ZFA,JGA,JFA,IRCRIT,
     &                        NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for B
C
         CALL WAVFUN_READ_REL(IFIL_RHSB,IT,1,ZGB,ZFB,JGB,JFB,IRCRIT,
     &                        NCPLWF,IKMCPLWF)
C
C---------------------------------- perturbation without enhancement (Z)
C                                     no rescaling of the wave functions
C
         CALL LINRESP_MAGNET_PERT_AUX(IT,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &                                JFB,HZ_PERT_LMTP,HAZ_ZZ_T,
     &                                HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T)
C
C-----------------------------------------------------------------------
C                       matrix elements of NABLA
C-----------------------------------------------------------------------
C                                     no rescaling of the wave functions
C
         CALL LINRESP_MAGNET_ME_NAB(IT,ERYD,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &                              JFB)
C
C-----------------------------------------------------------------------
C                       matrix elements of ALPHA
C-----------------------------------------------------------------------
C                                            rescales the wave functions
C
         CALL LINRESP_MAGNET_ME_ALF(IFIL_LHSA,IFIL_RHSB,ERYD,IT,ZGA,ZFA,
     &                              ZGB,ZFB,JGA,JFA,JGB,JFB)
C
C-----------------------------------------------------------------------
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C FOR TESTING ONLY: OVERWRITE NABLA MATRIX ELEMENTS WITH THOSE OF ALPHA
C
C      do i=17,19
C         HAZ_ZZ_T(:,:,:,:,I+3) = HAZ_ZZ_T(:,:,:,:,I)
C         HBZ_JZ_T(:,:,:,:,I+3) = HBZ_JZ_T(:,:,:,:,I)
C         HCZ_ZJ_T(:,:,:,:,I+3) = HCZ_ZJ_T(:,:,:,:,I)
C         HDZ_JJ_T(:,:,:,:,I+3) = HDZ_JJ_T(:,:,:,:,I)
C      end do
C
      DEALLOCATE (JGA,JGB,ZGA,ZGB)
C
      END
C*==linresp_magnet_pert_aux.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_PERT_AUX(IT,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &   JFB,H_PERT_LMTP,HA_ZZ_T,HB_JZ_T,HC_ZJ_T,HD_JJ_T)
C   ********************************************************************
C   *                                                                  *
C   *   set up the perturbation matrix functions for atom type IT      *
C   *                                                                  *
C   *    HA_ZZ(r) = Int_{0}^{r}       Z^x(r',E_a)  Delta H  Z(r',E_b)  *
C   *    HB_JZ(r) = Int_{r}^{r_crit}  J^x(r',E_a)  Delta H  Z(r',E_b)  *
C   *    HC_ZJ(r) = Int_{r}^{r_crit}  Z^x(r',E_a)  Delta H  J(r',E_b)  *
C   *    HD_JJ(r) = Int_{r}^{r_crit}  J^x(r',E_a)  Delta H  J(r',E_b)  *
C   *                                                                  *
C   *   for a scalar angular dependent perturbation                    *
C   *                                                                  *
C   *      Delta H(->r) = Sum_lm  H_lm(r) * Y_lm(^r)                   *
C   *                                                                  *
C   *   with non-vanishing terms indicated by the flag  KDVLMT         *
C   *                                                                  *
C   *   NOTE: H_lm has been convoluted with the shape functions        *
C   *         and the radial integration weights  R2DRDI               *
C   *                                                                  *
C   *   the endings A (B) for the wave functions Z and J indicate the  *
C   *   energies E_A(B) for  Delta G(E_a,E_b) =  G(E_a) Delta H G(E_b) *
C   *                                                                  *
C   *   for IREL = 3 (fully relativistic) it is implied that the       *
C   *   LHS solutions Z^x(r,E_a) and J^x(r,E_a) are read from file     *
C   *   IFIL_LHSA as these may differ from the RHS solutions for E_b   *
C   *                                                                  *
C   *   fully relativistic version                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LINRESP,ONLY:K_PERT_LMTP,IOPER_PERT,AMEG_LMOP,AMEF_LMOP,
     &    NPERT
      USE MOD_ANGMOM,ONLY:NKMMAX,NCPLWF,NKM,L_IKM
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_CONSTANTS,ONLY:
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF,NTMAX,NLMFPMAX
      IMPLICIT NONE
C*--LINRESP_MAGNET_PERT_AUX192
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IT
      COMPLEX*16 HA_ZZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HB_JZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HC_ZJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HD_JJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),ZFA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
      REAL*8 H_PERT_LMTP(NRMAX,NLMFPMAX,NTMAX,NPERT)
C
C Local variables
C
      COMPLEX*16 AME_F,AME_G,FJAJB(NRMAX),FJAZB(NRMAX),FZAJB(NRMAX),
     &           FZAZB(NRMAX),W_F,W_G,W_JFA,W_JGA,W_ZFA,W_ZGA
      INTEGER IA2,IB3,IM,IOPER,IPERT,IR,IRCRIT,KLMA2,KLMB3,LA,LAMA2,
     &        LAMB3,LB,LH,LMH,MH
C
C*** End of declarations rewritten by SPAG
C
      IM = IMT(IT)
      IRCRIT = JRCRI(IM)
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                     energy B -- LAMB3
      DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                     energy A -- LAMA2
         DO LAMA2 = 1,NKM
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
            DO IPERT = 1,7
C
               FZAZB(:) = 0D0
               FZAJB(:) = 0D0
               FJAZB(:) = 0D0
               FJAJB(:) = 0D0
C
C-----------------------------------------------------------------------
C                                                                IA2 IB3
               DO IB3 = 1,NCPLWF(LAMB3)
                  KLMB3 = IKMCPLWF(IB3,LAMB3)
                  LB = L_IKM(KLMB3)
C
                  DO IA2 = 1,NCPLWF(LAMA2)
                     KLMA2 = IKMCPLWF(IA2,LAMA2)
                     LA = L_IKM(KLMA2)
C
                     DO LH = ABS(LA-LB),(LA+LB),2
                        LMH = LH*LH
                        DO MH = -LH,LH
                           LMH = LMH + 1
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      check angular matrix element  AME_OP  and
C      KPERT flag for contribution LM of perturbation H_PERT
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                           IOPER = IOPER_PERT(IPERT)
                           AME_G = AMEG_LMOP(KLMA2,KLMB3,LMH,IOPER)
                           AME_F = AMEF_LMOP(KLMA2,KLMB3,LMH,IOPER)
C
                           IF ( ABS(AME_G).GT.1D-10 .AND. 
     &                          K_PERT_LMTP(LMH,IT,IPERT) ) THEN
C
                              DO IR = 1,IRCRIT
C
                                 W_G = AME_G*H_PERT_LMTP(IR,LMH,IT,
     &                                 IPERT)
                                 W_F = AME_F*H_PERT_LMTP(IR,LMH,IT,
     &                                 IPERT)
C
                                 W_ZGA = W_G*ZGA(IR,IA2,LAMA2)
                                 W_ZFA = W_F*ZFA(IR,IA2,LAMA2)
                                 FZAZB(IR) = FZAZB(IR)
     &                              + W_ZGA*ZGB(IR,IB3,LAMB3)
     &                              + W_ZFA*ZFB(IR,IB3,LAMB3)
                                 FZAJB(IR) = FZAJB(IR)
     &                              + W_ZGA*JGB(IR,IB3,LAMB3)
     &                              + W_ZFA*JFB(IR,IB3,LAMB3)
C
                                 W_JGA = W_G*JGA(IR,IA2,LAMA2)
                                 W_JFA = W_F*JFA(IR,IA2,LAMA2)
                                 FJAZB(IR) = FJAZB(IR)
     &                              + W_JGA*ZGB(IR,IB3,LAMB3)
     &                              + W_JFA*ZFB(IR,IB3,LAMB3)
                                 FJAJB(IR) = FJAJB(IR)
     &                              + W_JGA*JGB(IR,IB3,LAMB3)
     &                              + W_JFA*JFB(IR,IB3,LAMB3)
                              END DO
C
                           END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                        END DO
                     END DO
C                                                                  LH MH
C-----------------------------------------------------------------------
C
                  END DO
               END DO
C                                                                IA2 IB3
C-----------------------------------------------------------------------
C
               CALL CRADINT_R(IM,FZAZB,HA_ZZ_T(1,LAMA2,LAMB3,IT,IPERT))
C
               CALL CRADINT_INW_R(IM,FZAJB,
     &                            HC_ZJ_T(1,LAMA2,LAMB3,IT,IPERT))
               CALL CRADINT_INW_R(IM,FJAZB,
     &                            HB_JZ_T(1,LAMA2,LAMB3,IT,IPERT))
               CALL CRADINT_INW_R(IM,FJAJB,
     &                            HD_JJ_T(1,LAMA2,LAMB3,IT,IPERT))
C
            END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
         END DO
C                                                     energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                     energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
C*==linresp_magnet_response.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_RESPONSE(IFIL_RHSA,MSSTA,TAUTA,DMATTA,
     &   DTILTA,IFIL_LHSB,MSSTB,TAUTB,DMATTB,DTILTB,MEZZ,MEZJ,WE,ERYD)
C   ********************************************************************
C   *                                                                  *
C   *  calculates r-integrated response to a perurbation               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKMMAX,NMEMAX,NCPLWF,NKM,WKM1,WKM2,WKM3,NKMQ
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_CALCMODE,ONLY:IREL,GF_CONV_RH
      USE MOD_SITES,ONLY:IQAT,IQBOT_CHI,IQTOP_CHI,ITOQ,NOQ
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,ITBOT,ITTOP,NTMAX,IKMCPLWF,CONC
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_LINRESP,ONLY:CHIZ,CHI_DK,HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,
     &    HDZ_JJ_T,NPERT,D0Z,T0Z,D1Z,T1Z,DZ,TZ,DIJZ,TIJZ,THZT_OFF_TP,
     &    THZT_DIA_P,MIRR2_OP,MIRR3_OP,MIRR4_OP,T1Z_DK,T0Z_DK,TZ_DK
      IMPLICIT NONE
C*--LINRESP_MAGNET_RESPONSE356
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_RESPONSE')
      INTEGER NC
      PARAMETER (NC=3)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,WE
      INTEGER IFIL_LHSB,IFIL_RHSA
      COMPLEX*16 DMATTA(NKMMAX,NKMMAX,NTMAX),DMATTB(NKMMAX,NKMMAX,NTMAX)
     &           ,DTILTA(NKMMAX,NKMMAX,NTMAX),
     &           DTILTB(NKMMAX,NKMMAX,NTMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           TAUTA(NKMMAX,NKMMAX,NTMAX),TAUTB(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CSUMZ,CWGT,HZ_TP(:,:,:,:),JFA(:,:,:),JFB(:,:,:),
     &           JGA(:,:,:),JGB(:,:,:),TMATA(:,:),TMATB(:,:),ZFA(:,:,:),
     &           ZFB(:,:,:),ZGA(:,:,:),ZGB(:,:,:)
      INTEGER I,IA_ERR,IC,IM,IPERT,IQ,IRCRIT,IRTOP,IT,J,JO,JQ,JT,K1,K2,
     &        LM1,LM2,LM3,LM4,M,N,NKMSQ
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TMATA,TMATB,HZ_TP
      ALLOCATABLE JGA,JGB,ZGA,ZGB,JFA,JFB,ZFA,ZFB
C
      IF ( IREL.EQ.3 ) THEN
         ALLOCATE (JFA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (JFB(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFA(NRMAX,NCPLWFMAX,NKMMAX))
         ALLOCATE (ZFB(NRMAX,NCPLWFMAX,NKMMAX))
      END IF
C
      ALLOCATE (JGA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGA(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGB(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ZGB')
C
      ALLOCATE (HZ_TP(NKMMAX,NKMMAX,NTMAX,NPERT))
      ALLOCATE (TMATA(NKMMAX,NKMMAX),TMATB(NKMMAX,NKMMAX))
C
      TMATA(1:NKMMAX,1:NKMMAX) = C0
      TMATB(1:NKMMAX,1:NKMMAX) = C0
C
C-----------------------------------------------------------------------
C           combination of TAU and perturbation terms
C-----------------------------------------------------------------------
      M = NKMMAX
      IF ( .NOT.ALLOCATED(THZT_OFF_TP) ) THEN
         ALLOCATE (THZT_OFF_TP(M,M,NTMAX,NPERT))
         ALLOCATE (THZT_DIA_P(M,M,NPERT),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate THZT')
      END IF
C
C=======================================================================
C   H_TP(JT) = x(JT) * D~(JT) [* m(JT)] * H(JT) [* m(JT)] * D(JT)
C
C   renormalize type dependent matrix elements M for perturbation
C   GF_CONV_RH:  change to ZJ convention      M ->   m M m (if needed)
C   CPA       :  include projection matrices  M ->  D~ M D
C
C   H = HAZ_ZZ_T   may be evaluated in RH- or ZJ-convention
C   H_TP(JT)       to be coupled to TAU always in ZJ-convention
C=======================================================================
C
      N = NKM
C
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      DO IPERT = 1,NPERT
C   TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         DO JT = ITBOT,ITTOP
C
            CWGT = CONC(JT)
            IM = IMT(JT)
            IRCRIT = JRCRI(IM)
C
            WKM2(1:N,1:N) = HAZ_ZZ_T(IRCRIT,1:N,1:N,JT,IPERT)
C
            IF ( GF_CONV_RH ) THEN
C
               WKM1(1:N,1:N) = MATMUL(MSSTA(1:N,1:N,JT),WKM2)
               WKM2(1:N,1:N) = MATMUL(WKM1(1:N,1:N),MSSTB(1:N,1:N,JT))
C
               WKM1(1:N,1:N) = MATMUL(MSSTA(1:N,1:N,JT),WKM3)
               WKM3(1:N,1:N) = MATMUL(WKM1(1:N,1:N),MSSTB(1:N,1:N,JT))
C
            END IF
C
            WKM1(1:N,1:N) = MATMUL(DTILTA(1:N,1:N,JT),WKM2)
C
            HZ_TP(1:N,1:N,JT,IPERT) = CWGT*MATMUL(WKM1(1:N,1:N),DMATTB(1
     &                                :N,1:N,JT))
C
         END DO
C   TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      END DO
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C      CALL  CMATSTR('H PERT   1',HZ_TP(1,1,1,1),N,M,3,3,0,1D-8,6)
C      CALL  CMATSTR('H PERT   2',HZ_TP(1,1,1,2),N,M,3,3,0,1D-8,6)
C      CALL  CMATSTR('H PERT   3',HZ_TP(1,1,1,3),N,M,3,3,0,1D-8,6)
C
      NKMSQ = NKM*NKM
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRCRIT = JRCRI(IM)
         IRTOP = JRCRI(IM)
         IQ = IQAT(1,IT)
         M = NKMMAX
         N = NKMQ(IQ)
C
C-----------------------------------------------------------------------
C         TAUT:  type projected TAU containing DMAT in case of CPA
C         TMAT:  site diagonal multiple scattering matrix G or TAU
C-----------------------------------------------------------------------
C
         CALL TAUGFCONV(MSSTA(1,1,IT),TAUTA(1,1,IT),TMATA)
C
         CALL TAUGFCONV(MSSTB(1,1,IT),TAUTB(1,1,IT),TMATB)
C
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                   loop over the various perturbation terms
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C-----------------------------------------------------------------------
C        IPERT
C
C        1        DOS
C        2 -  4   SPN
C        5 -  7   ORB
C        8 - 16   ( (Y_LM1 x sigma_IPOL), LM1=2,4), IPOL=1,3)
C       17 - 19             (sigma_IPOL),           IPOL=1,3)
C       20 - 22             (nabla_IPOL),           IPOL=1,3)
C
C                  LM1 =  2,  3,  4
C                   m  = -1,  0, +1
C                   c  =  y,  z,  x
C-----------------------------------------------------------------------
C
         THZT_OFF_TP(:,:,:,:) = 0D0
C
         DO IPERT = 1,16
C-----------------------------------------------------------------------
C           Delta G_off(1) site off-diagonal contribution       1st step
C-----------------------------------------------------------------------
C  THT_TP(JT) =  SUM{JQ} TAU(IQ,IT;JQ,JT) * H_TP(JT) * TAU(JQ,JT;IQ,IT)
C
C   represents perturbation H_TP(JT) of type JT on all sites JQ
C   leading to a response for type IT on central site IQ
C   with all quantities defined according to ZJ-convention
C-----------------------------------------------------------------------
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            DO JQ = IQBOT_CHI,IQTOP_CHI
               DO JO = 1,NOQ(JQ)
                  JT = ITOQ(JO,JQ)
C
                  K1 = (IQ-1)*NKMSQ
                  DO LM1 = 1,NKM
                     DO LM4 = 1,NKM
                        K1 = K1 + 1
C
                        CSUMZ = 0D0
                        K2 = (JQ-1)*NKMSQ
                        DO LM2 = 1,NKM
                           DO LM3 = 1,NKM
                              K2 = K2 + 1
C
                              CSUMZ = CSUMZ + CHIZ(K1,K2,1)
     &                                *HZ_TP(LM2,LM3,JT,IPERT)
C
                           END DO
                        END DO
C
                        THZT_OFF_TP(LM1,LM4,JT,IPERT)
     &                     = THZT_OFF_TP(LM1,LM4,JT,IPERT) + CSUMZ
C
                     END DO
                  END DO
C
               END DO
            END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C-----------------------------------------------------------------------
C           Delta G_off(1) site off-diagonal contribution       2st step
C-----------------------------------------------------------------------
C  transform THT_TP(JT) to allow use of matrix elements for
C  observables  MZBZA_O  in convention X = RH, ZJ as set by  GF_CONV_RH
C
C     Trace  [ m(IT)] * M(IT) [* m(IT)] * D(IT) * THT_TP(JT) * D~(IT)
C
C   = Trace  M(IT) [* m(IT)] * D(IT) * THT_TP(JT) * D~(IT) [* m(IT)]
C
C   = Trace  M(IT) * THT_OFF_TP(JT)
C-----------------------------------------------------------------------
C
            DO JT = ITBOT,ITTOP
C
               WKM1(1:N,1:N) = MATMUL(DMATTA(1:N,1:N,IT),THZT_OFF_TP(1:N
     &                         ,1:N,JT,IPERT))
               WKM2(1:N,1:N) = MATMUL(WKM1(1:N,1:N),DTILTB(1:N,1:N,IT))
C
               IF ( GF_CONV_RH ) THEN
C
                  WKM1(1:N,1:N) = MATMUL(MSSTA(1:N,1:N,IT),WKM2(1:N,1:N)
     &                            )
C
                  THZT_OFF_TP(1:N,1:N,JT,IPERT)
     &               = MATMUL(WKM1(1:N,1:N),MSSTB(1:N,1:N,IT))
C
               ELSE
C
                  THZT_OFF_TP(1:N,1:N,JT,IPERT) = WKM2(1:N,1:N)
C
               END IF
C
            END DO
C
C-----------------------------------------------------------------------
C            Delta G_dia(1) site diagonal contribution
C-----------------------------------------------------------------------
C   THT_DIA_P(IT ) = TMAT(IT) * HA_ZZ_T(IT) * TMAT(IT)
C
C   all quantities in RH- or ZJ-convention as set by     GF_CONV_RH
C-----------------------------------------------------------------------
C
            DO J = 1,NKM
               DO I = 1,NKM
                  WKM2(I,J) = HAZ_ZZ_T(IRTOP,I,J,IT,IPERT)
               END DO
            END DO
C
            WKM1(1:N,1:N) = MATMUL(TMATA(1:N,1:N),WKM2(1:N,1:N))
C
            THZT_DIA_P(1:N,1:N,IPERT) = MATMUL(WKM1(1:N,1:N),TMATB(1:N,1
     &                                  :N))
C
         END DO
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C-----------------------------------------------------------------------
C ------------------------------------------ read in wavefunctions for A
C                                                use arrays *B as buffer
C
         CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGB,ZFB,JGB,JFB,IRCRIT,
     &                        NCPLWF,IKMCPLWF)
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
         CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZGB,ZFB,JGB,JFB,ZGA,ZFA,JGA,JFA,
     &                               NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for B
C
         CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGB,ZFB,JGB,JFB,IRCRIT,
     &                        NCPLWF,IKMCPLWF)
C
C-------------------------------------- response without enhancement (Z)
C
         CALL LINRESP_MAGNET_RESP_AUX(IT,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &                                JFB,HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,
     &                                HDZ_JJ_T,ERYD,1,16)
C
         CALL LINRESP_MAGNET_RESP_SUM(IT,MEZZ,MEZJ,TMATA,TMATB,
     &                                THZT_DIA_P,THZT_OFF_TP,D0Z,T0Z,
     &                                D1Z,T1Z,DZ,TZ,DIJZ,TIJZ,WE,1,16)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO IC = 1,NC
C
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                   loop over the various perturbation terms
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C-----------------------------------------------------------------------
C        IPERT
C
C        1        DOS
C        2 -  4   SPN
C        5 -  7   ORB
C        8 - 16   ( (Y_LM1 x sigma_IPOL), LM1=2,4), IPOL=1,3)
C       17 - 19             (sigma_IPOL),           IPOL=1,3)
C       20 - 22             (nabla_IPOL),           IPOL=1,3)
C
C                  LM1 =  2,  3,  4
C                   m  = -1,  0, +1
C                   c  =  y,  z,  x
C-----------------------------------------------------------------------
C
            DO IPERT = 17,22
               THZT_OFF_TP(:,:,:,IPERT) = C0
C-----------------------------------------------------------------------
C           Delta G_off(1) site off-diagonal contribution       1st step
C-----------------------------------------------------------------------
C  THT_TP(JT) =  SUM{JQ} TAU(IQ,IT;JQ,JT) * H_TP(JT) * TAU(JQ,JT;IQ,IT)
C
C   represents perturbation H_TP(JT) of type JT on all sites JQ
C   leading to a response for type IT on central site IQ
C   with all quantities defined according to ZJ-convention
C-----------------------------------------------------------------------
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
               DO JQ = IQBOT_CHI,IQTOP_CHI
                  DO JO = 1,NOQ(JQ)
                     JT = ITOQ(JO,JQ)
C
                     K1 = (IQ-1)*NKMSQ
                     DO LM1 = 1,NKM
                        DO LM4 = 1,NKM
                           K1 = K1 + 1
C
                           CSUMZ = 0D0
                           K2 = (JQ-1)*NKMSQ
                           DO LM2 = 1,NKM
                              DO LM3 = 1,NKM
                                 K2 = K2 + 1
C
                                 CSUMZ = CSUMZ + CHI_DK(K1,K2,IC)
     &                              *HZ_TP(LM2,LM3,JT,IPERT)
C
                              END DO
                           END DO
C
                           THZT_OFF_TP(LM1,LM4,JT,IPERT)
     &                        = THZT_OFF_TP(LM1,LM4,JT,IPERT) + CSUMZ
C
                        END DO
                     END DO
C
                  END DO
               END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C-----------------------------------------------------------------------
C           Delta G_off(1) site off-diagonal contribution       2st step
C-----------------------------------------------------------------------
C  transform THT_TP(JT) to allow use of matrix elements for
C  observables  MZBZA_O  in convention X = RH, ZJ as set by  GF_CONV_RH
C
C     Trace  [ m(IT)] * M(IT) [* m(IT)] * D(IT) * THT_TP(JT) * D~(IT)
C
C   = Trace  M(IT) [* m(IT)] * D(IT) * THT_TP(JT) * D~(IT) [* m(IT)]
C
C   = Trace  M(IT) * THT_OFF_TP(JT)
C-----------------------------------------------------------------------
C
               DO JT = ITBOT,ITTOP
C
                  WKM1(1:N,1:N) = MATMUL(DMATTA(1:N,1:N,IT),THZT_OFF_TP(
     &                            1:N,1:N,JT,IPERT))
                  WKM2(1:N,1:N) = MATMUL(WKM1(1:N,1:N),DTILTB(1:N,1:N,IT
     &                            ))
C
                  IF ( GF_CONV_RH ) THEN
C
                     WKM1(1:N,1:N) = MATMUL(MSSTA(1:N,1:N,IT),WKM2(1:N,1
     &                               :N))
C
                     THZT_OFF_TP(1:N,1:N,JT,IPERT)
     &                  = MATMUL(WKM1(1:N,1:N),MSSTB(1:N,1:N,IT))
C
                  ELSE
C
                     THZT_OFF_TP(1:N,1:N,JT,IPERT) = WKM2(1:N,1:N)
C
                  END IF
C
               END DO
C
C-----------------------------------------------------------------------
C            Delta G_dia(1) site diagonal contribution
C-----------------------------------------------------------------------
C   THT_DIA_P(IT ) = TMAT(IT) * HA_ZZ_T(IT) * TMAT(IT)
C
C   all quantities in RH- or ZJ-convention as set by     GF_CONV_RH
C-----------------------------------------------------------------------
C
               DO J = 1,NKM
                  DO I = 1,NKM
                     WKM2(I,J) = HAZ_ZZ_T(IRTOP,I,J,IT,IPERT)
                  END DO
               END DO
C
               WKM1(1:N,1:N) = MATMUL(TMATA(1:N,1:N),WKM2(1:N,1:N))
C
               THZT_DIA_P(1:N,1:N,IPERT)
     &            = MATMUL(WKM1(1:N,1:N),TMATB(1:N,1:N))
C
            END DO
C pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
C
C-----------------------------------------------------------------------
C   FROM       G * G = Delta_G_off(1) + Delta_G_dia(1)
C                    + Delta_G(2) + Delta_G(3) + Delta_G(4)
C-----------------------------------------------------------------------
C   KEEP ONLY
C------------------------- Delta G_off(1) site off-diagonal contribution
C
C-------------------------------------- response without enhancement (Z)
C
            CALL LINRESP_MAGNET_RESP_AUX(IT,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &         JFB,HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T,ERYD,17,22)
C
            DO IPERT = 17,22
               THZT_DIA_P(1:N,1:N,IPERT) = C0
               MIRR2_OP(:,:,:,IPERT) = C0
               MIRR3_OP(:,:,:,IPERT) = C0
               MIRR4_OP(:,:,:,IPERT) = C0
            END DO
C
            CALL LINRESP_MAGNET_RESP_SUM(IT,MEZZ,MEZJ,TMATA,TMATB,
     &         THZT_DIA_P,THZT_OFF_TP,D0Z,T0Z_DK(1,1,1,IC),D1Z,
     &         T1Z_DK(1,1,1,IC),DZ,TZ_DK(1,1,1,IC),DIJZ,TIJZ,WE,17,22)
C
C
         END DO
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
C-----------------------------------------------------------------------
C      read wave functions and calculate induced density   RHO2NS_GG
C-----------------------------------------------------------------------
C
C ------------------------------------------ read in wavefunctions for A
C
Ccc         CALL WAVFUN_READ_REL(IFIL_RHSA,IT,1,ZGA,ZFA,JGA,JFA,IRTOPccc,
Ccc     &                            NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for B
C
Ccc         CALL WAVFUN_READ_REL(IFIL_LHSB,IT,1,ZGB,ZFB,JGB,JFB,IRTOPccc,
Ccc     &                            NCPLWF,IKMCPLWF)
C
C-----------------------------------------------------------------------
C
Ccc         CALL LINRESP_MAGNET_DENSITY_REL(IT,IRTOP,WE,TMATA,TMATB,
Ccc     &      THZT_OFF_TP,THZT_DIA_P,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,JFB)
C
C ======================================================================
C            calculate expectation values via the standard scheme
C                        <A> = -1/PI Im Trace A G
C ======================================================================
C
         CALL LINRESP_MAGNET_STANDARD(IT,WE,TMATA)
C
C ======================================================================
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DEALLOCATE (JGA,JGB,ZGA,ZGB)
C
      END
C*==linresp_magnet_resp_aux.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_RESP_AUX(IT,ZGA,ZFA,ZGB,ZFB,JGA,JFA,JGB,
     &   JFB,HA_ZZ_T,HB_JZ_T,HC_ZJ_T,HD_JJ_T,ERYD,IPERT1,IPERT2)
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
CIE,
      USE MOD_ENERGY,ONLY:EFERMI
      USE MOD_LINRESP,ONLY:NPERT,NOBSE,AMEG_LMOP,AMEF_LMOP,MZBZA_O,
     &    MIRR2_OP,MIRR3_OP,MIRR4_OP,IDOS
      USE MOD_ANGMOM,ONLY:NKMMAX,NCPLWF,NKM
      USE MOD_RMESH,ONLY:NRMAX,JRCRI
      USE MOD_CONSTANTS,ONLY:C0,SQRT_4PI
      USE MOD_TYPES,ONLY:IMT,NTMAX,NCPLWFMAX,IKMCPLWF
      IMPLICIT NONE
C*--LINRESP_MAGNET_RESP_AUX856
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IPERT1,IPERT2,IT
      COMPLEX*16 HA_ZZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HB_JZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HC_ZJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           HD_JJ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),ZFA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 AFF,AGG,AME_F,AME_G,FAUX(:),HA_ZZX(:,:),HB_JZX(:,:),
     &           HC_ZJX(:,:),HD_JJX(:,:),IAUX_TOP,WIRR,WREG,W_ZFA,W_ZFB,
     &           W_ZGA,W_ZGB
      INTEGER IA1,IB3,IB4,IFLAG2,IKM1,IKM4,IKMA1,IKMB3,IKMB4,IM,IOBSE,
     &        IOPER,IPERT,IR,IRCRIT,LAMA1,LAMA2,LAMB3,LAMB4
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FAUX,HA_ZZX,HB_JZX,HC_ZJX,HD_JJX
C
      ALLOCATE (FAUX(NRMAX))
      ALLOCATE (HA_ZZX(NRMAX,NPERT),HB_JZX(NRMAX,NPERT))
      ALLOCATE (HC_ZJX(NRMAX,NPERT),HD_JJX(NRMAX,NPERT))
C
      IM = IMT(IT)
      IRCRIT = JRCRI(IM)
C
      MZBZA_O(:,:,:) = C0
      MIRR2_OP(:,:,:,:) = C0
      MIRR3_OP(:,:,:,:) = C0
      MIRR4_OP(:,:,:,:) = C0
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                     energy B -- LAMB4
      DO LAMB4 = 1,NKM
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                     energy A -- LAMA1
         DO LAMA1 = 1,NKM
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            DO IOBSE = 1,NOBSE
C
               FAUX(1:IRCRIT) = C0
C
C-----------------------------------------------------------------------
               DO IB4 = 1,NCPLWF(LAMB4)
                  IKM4 = IKMCPLWF(IB4,LAMB4)
C
                  DO IA1 = 1,NCPLWF(LAMA1)
                     IKM1 = IKMCPLWF(IA1,LAMA1)
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      check angular matrix element  AME_OP for observable   IOBSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                     IOPER = IDOS
C
                     AME_G = AMEG_LMOP(IKM4,IKM1,1,IOPER)*SQRT_4PI
                     AME_F = AMEF_LMOP(IKM4,IKM1,1,IOPER)*SQRT_4PI
C
                     AME_G = AME_G*(ERYD-EFERMI)
                     AME_F = AME_F*(ERYD-EFERMI)
C
                     IF ( ABS(AME_G).GT.1D-6 ) THEN
C
                        DO IR = 1,IRCRIT
                           AGG = AME_G*ZGB(IR,IB4,LAMB4)
     &                           *ZGA(IR,IA1,LAMA1)
                           AFF = AME_F*ZFB(IR,IB4,LAMB4)
     &                           *ZFA(IR,IA1,LAMA1)
                           FAUX(IR) = FAUX(IR) + (AGG+AFF)
                        END DO
C
                     END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                  END DO
               END DO
C-----------------------------------------------------------------------
C
               CALL CRADINT(IM,FAUX,MZBZA_O(LAMB4,LAMA1,IOBSE))
C
            END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
         END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                     energy B -- LAMB3
      DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                     energy A -- LAMA2
         DO LAMA2 = 1,NKM
C
            DO IPERT = IPERT1,IPERT2
               HA_ZZX(1:IRCRIT,IPERT) = HA_ZZ_T(1:IRCRIT,LAMA2,LAMB3,IT,
     &                                  IPERT)
               HC_ZJX(1:IRCRIT,IPERT) = HC_ZJ_T(1:IRCRIT,LAMA2,LAMB3,IT,
     &                                  IPERT)
               HB_JZX(1:IRCRIT,IPERT) = HB_JZ_T(1:IRCRIT,LAMA2,LAMB3,IT,
     &                                  IPERT)
               HD_JJX(1:IRCRIT,IPERT) = HD_JJ_T(1:IRCRIT,LAMA2,LAMB3,IT,
     &                                  IPERT)
            END DO
C
C***********************************************************************
C***********************************************************************
C     TERM 1    MZBZA(4,1) * TAU(1,2,a) * MZAZB(2,3) * TAU(3,4,b)
C***********************************************************************
C***********************************************************************
C
C
C***********************************************************************
C***********************************************************************
C     TERM 2    - TAU(3,4,b) * MIRR_2(3,4)               LAMA1 = LAMA2
C***********************************************************************
C***********************************************************************
C
            LAMA1 = LAMA2
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                     energy B -- LAMB4
            DO LAMB4 = 1,NKM
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
               DO IOBSE = 1,NOBSE
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
                  DO IPERT = IPERT1,IPERT2
C
                     FAUX(1:IRCRIT) = C0
C
                     IFLAG2 = 0
C
C-----------------------------------------------------------------------
C                                                                IA1 IB4
                     DO IA1 = 1,NCPLWF(LAMA1)
                        IKMA1 = IKMCPLWF(IA1,LAMA1)
C
                        DO IB4 = 1,NCPLWF(LAMB4)
                           IKMB4 = IKMCPLWF(IB4,LAMB4)
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      check angular matrix element  AME_OP for observable   IOBSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                           IOPER = IDOS
C
                           AME_G = AMEG_LMOP(IKMA1,IKMB4,1,IOPER)
     &                             *SQRT_4PI
                           AME_F = AMEF_LMOP(IKMA1,IKMB4,1,IOPER)
     &                             *SQRT_4PI
C
                           AME_G = AME_G*(ERYD-EFERMI)
                           AME_F = AME_F*(ERYD-EFERMI)
C
                           IF ( ABS(AME_G).GT.1D-10 ) THEN
C
                              IFLAG2 = 1
C
                              DO IR = 1,IRCRIT
                                 W_ZGB = AME_G*ZGB(IR,IB4,LAMB4)
                                 W_ZFB = AME_F*ZFB(IR,IB4,LAMB4)
                                 WREG = W_ZGB*JGA(IR,IA1,LAMA1)
     &                                  + W_ZFB*JFA(IR,IA1,LAMA1)
                                 WIRR = W_ZGB*ZGA(IR,IA1,LAMA1)
     &                                  + W_ZFB*ZFA(IR,IA1,LAMA1)
C
                                 FAUX(IR) = FAUX(IR)
     &                              + WREG*HA_ZZX(IR,IPERT)
     &                              + WIRR*HB_JZX(IR,IPERT)
                              END DO
C
                           END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                        END DO
                     END DO
C                                                                IA1 IB4
C-----------------------------------------------------------------------
C
C***********************************************************************
                     IF ( IFLAG2.EQ.1 ) THEN
C
                        CALL CRADINT(IM,FAUX,IAUX_TOP)
C
                        MIRR2_OP(LAMB4,LAMB3,IOBSE,IPERT)
     &                     = MIRR2_OP(LAMB4,LAMB3,IOBSE,IPERT)
     &                     + IAUX_TOP
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
                     END IF
C***********************************************************************
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
                  END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
               END DO
            END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB LAMB4
C
C
C***********************************************************************
C***********************************************************************
C     TERM 3    - TAU(1,2,a) * MIRR_3(1,2)               LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                     energy A -- LAMA1
            DO LAMA1 = 1,NKM
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
               DO IOBSE = 1,NOBSE
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
                  DO IPERT = IPERT1,IPERT2
C
                     FAUX(1:IRCRIT) = C0
C
                     IFLAG2 = 0
C
C-----------------------------------------------------------------------
C                                                                IA1 IB3
                     DO IA1 = 1,NCPLWF(LAMA1)
                        IKMA1 = IKMCPLWF(IA1,LAMA1)
C
                        DO IB3 = 1,NCPLWF(LAMB3)
                           IKMB3 = IKMCPLWF(IB3,LAMB3)
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      check angular matrix element  AME_OP for observable   IOBSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                           IOPER = IDOS
C
                           AME_G = AMEG_LMOP(IKMA1,IKMB3,1,IOPER)
     &                             *SQRT_4PI
                           AME_F = AMEF_LMOP(IKMA1,IKMB3,1,IOPER)
     &                             *SQRT_4PI
C
                           AME_G = AME_G*(ERYD-EFERMI)
                           AME_F = AME_F*(ERYD-EFERMI)
C
                           IF ( ABS(AME_G).GT.1D-10 ) THEN
C
                              IFLAG2 = 1
C
                              DO IR = 1,IRCRIT
                                 W_ZGA = AME_G*ZGA(IR,IA1,LAMA1)
                                 W_ZFA = AME_F*ZFA(IR,IA1,LAMA1)
                                 WREG = JGB(IR,IB3,LAMB3)
     &                                  *W_ZGA + JFB(IR,IB3,LAMB3)*W_ZFA
                                 WIRR = ZGB(IR,IB3,LAMB3)
     &                                  *W_ZGA + ZFB(IR,IB3,LAMB3)*W_ZFA
C
                                 FAUX(IR) = FAUX(IR)
     &                              + WREG*HA_ZZX(IR,IPERT)
     &                              + WIRR*HC_ZJX(IR,IPERT)
                              END DO
C
                           END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                        END DO
                     END DO
C                                                                IA1 IB3
C-----------------------------------------------------------------------
C
C***********************************************************************
                     IF ( IFLAG2.EQ.1 ) THEN
C
                        CALL CRADINT(IM,FAUX,IAUX_TOP)
C
                        MIRR3_OP(LAMA1,LAMA2,IOBSE,IPERT)
     &                     = MIRR3_OP(LAMA1,LAMA2,IOBSE,IPERT)
     &                     + IAUX_TOP
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
                     END IF
C***********************************************************************
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
                  END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
               END DO
            END DO
CcCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA LAMA1
C
C
C***********************************************************************
C***********************************************************************
C     TERM 4    MIRR_4(1,3)       LAMA1 = LAMA2  and   LAMB3 = LAMB4
C***********************************************************************
C***********************************************************************
C
            LAMA1 = LAMA2
C
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            DO IOBSE = 1,NOBSE
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
               DO IPERT = IPERT1,IPERT2
C
                  FAUX(1:IRCRIT) = C0
C
                  IFLAG2 = 0
C
C-----------------------------------------------------------------------
C                                                                IA1 IB3
                  DO IA1 = 1,NCPLWF(LAMA1)
                     IKMA1 = IKMCPLWF(IA1,LAMA1)
C
                     DO IB3 = 1,NCPLWF(LAMB3)
                        IKMB3 = IKMCPLWF(IB3,LAMB3)
C
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C      check angular matrix element  AME_OP for observable   IOBSE
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                        IOPER = IDOS
C
                        AME_G = AMEG_LMOP(IKMA1,IKMB3,1,IOPER)*SQRT_4PI
                        AME_F = AMEF_LMOP(IKMA1,IKMB3,1,IOPER)*SQRT_4PI
C
                        AME_G = AME_G*(ERYD-EFERMI)
                        AME_F = AME_F*(ERYD-EFERMI)
C
                        IF ( ABS(AME_G).GT.1D-10 ) THEN
C
                           IFLAG2 = 1
C
                           DO IR = 1,IRCRIT
C
                              WREG = AME_G*JGB(IR,IB3,LAMB3)
     &                               *JGA(IR,IA1,LAMA1)
     &                               + AME_F*JFB(IR,IB3,LAMB3)
     &                               *JFA(IR,IA1,LAMA1)
                              WIRR = AME_G*ZGB(IR,IB3,LAMB3)
     &                               *ZGA(IR,IA1,LAMA1)
     &                               + AME_F*ZFB(IR,IB3,LAMB3)
     &                               *ZFA(IR,IA1,LAMA1)
C
                              FAUX(IR) = FAUX(IR)
     &                           + WREG*HA_ZZX(IR,IPERT)
     &                           + WIRR*HD_JJX(IR,IPERT)
C
                           END DO
C
                        END IF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
                     END DO
                  END DO
C                                                                IA1 IB3
C-----------------------------------------------------------------------
C
C***********************************************************************
C
C
                  IF ( IFLAG2.EQ.1 )
     &                 CALL CRADINT(IM,FAUX,MIRR4_OP(LAMB3,LAMA1,IOBSE,
     &                 IPERT))
C***********************************************************************
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
               END DO
Cooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            END DO
         END DO
C                                                     energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                     energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
