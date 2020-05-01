C*==fpvintra.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPVINTRA(USE_KLMFP,CMNTMTT,CMNTIST,R2RHO,VNEW,
     &                    LMRGNT_VSF,NRGNT_VSF_LM,RGNT_VSF,NRGNT_VSF)
C   ********************************************************************
C   *                                                                  *
C   * - initialize the spin-averaged potential V                       *
C   * - calculate the spin independent electron-intracell-potentials   *
C   *   contribution and the charge moments of given charge densities  *
C   *                                                                  *
C   * the intracell-potential V is expanded into spherical harmonics   *
C   * LM-term of V of the representive atom i:                         *
C   *                                                                  *
C   *                8pi        r      r'** l                          *
C   *  V(r,LM,IT) =  ----- *  [  S dr' --------   R2RHO(r',LM,IT,1)    *
C   *               2*l+1       0     r **(l+1)                        *
C   *                                                                  *
C   *                           rcut    r ** l                         *
C   *                         +  S dr' ---------   R2RHO(r',LM,IT,1) ] *
C   *                            r     r' **(l+1)                      *
C   *                                                                  *
C   * LM-term of charge moment of the representive atom i:             *
C   *                                                                  *
C   *                            rcut                                  *
C   *          CMNTT(LM,IT) =    S dr r** l R2RHO(r,LM,IT,1)           *
C   *                             0                                    *
C   *                                                                  *
C   * attention:  R2RHO(...,1) is the real charge density times r**2   *
C   *             expanded into spherical harmonics.                   *
C   *                                                                  *
C   * NOTE: all quantities refer to the LOCAL frame                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,DRDI,JRCUT,NPAN,FLMSF,KLMSF,NRMAX,NSFMAX,
     &    ISFLM
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,KLMFP,NLFP,IMT,NTMAX,NLMFPMAX
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--FPVINTRA39
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NRGNT_VSF
      LOGICAL USE_KLMFP
      REAL*8 CMNTIST(NLMFPMAX,NTMAX),CMNTMTT(NLMFPMAX,NTMAX),
     &       R2RHO(NRMAX,NLMFPMAX,NTMAX,3),RGNT_VSF(NRGNT_VSF),
     &       VNEW(NRMAX,NLMFPMAX,NTMAX)
      INTEGER LMRGNT_VSF(NRGNT_VSF,3),NRGNT_VSF_LM(0:NLMFPMAX)
C
C Local variables
C
      REAL*8 FAC,RL,RPWL(:),V1(:),V2(:),VINT1(:),VINT2(:)
      INTEGER IM,IR,IRCRIT,IRMTIN,IRSF,ISF,IT,J,L,LM,LM2,LM3,M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RPWL,V1,V2,VINT1,VINT2
C
      ALLOCATE (V1(NRMAX),V2(NRMAX),RPWL(NRMAX))
      ALLOCATE (VINT1(NRMAX),VINT2(NRMAX))
C
      CMNTIST(1:NLMFPMAX,1:NTMAX) = 0D0
      CMNTMTT(1:NLMFPMAX,1:NTMAX) = 0D0
      VNEW(1:NRMAX,1:NLMFPMAX,1:NTMAX) = 0.0D0
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT:DO IT = ITBOT,ITTOP
C
         IM = IMT(IT)
         IRMTIN = JRCUT(1,IM)
         IRCRIT = JRCUT(NPAN(IM),IM)
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         LOOP_L:DO L = 0,NLFP - 1
C
            IF ( L.EQ.0 ) THEN
               RPWL(1:NRMAX) = 1D0
            ELSE
               DO IR = 1,IRCRIT
                  RPWL(IR) = RPWL(IR)*R(IR,IM)
               END DO
            END IF
C
            FAC = 8.0D0*PI/DBLE(2*L+1)
C
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            LOOP_M:DO M = -L,L
C
               LM = L*L + L + M + 1
C
C=======================================================================
               IF ( .NOT.USE_KLMFP .OR. KLMFP(LM,IT).NE.0 ) THEN
C
C-----------------------------------------------------------------------
C     set up of the integrands V1 and V2
C-----------------------------------------------------------------------
                  V1(1) = 0.0D0
                  V2(1) = 0.0D0
C
                  DO IR = 1,IRMTIN
                     RL = RPWL(IR)
                     V1(IR) = DRDI(IR,IM)*R2RHO(IR,LM,IT,1)*RL
                     V2(IR) = DRDI(IR,IM)*R2RHO(IR,LM,IT,1)/R(IR,IM)/RL
                  END DO
C
C-----------------------------------------------------------------------
C     convolute charge density in interstial with shape function
C-----------------------------------------------------------------------
                  DO IR = IRMTIN + 1,IRCRIT
                     V1(IR) = 0.0D0
                  END DO
C
                  DO J = NRGNT_VSF_LM(LM-1) + 1,NRGNT_VSF_LM(LM)
                     LM2 = LMRGNT_VSF(J,2)
                     LM3 = LMRGNT_VSF(J,3)
                     ISF = ISFLM(LM3,IM)
                     IF ( KLMSF(LM3,IM).EQ.1 .AND. ISF.LE.NSFMAX ) THEN
                        DO IR = IRMTIN + 1,IRCRIT
                           IRSF = IR - IRMTIN
                           V1(IR) = V1(IR) + RGNT_VSF(J)
     &                              *R2RHO(IR,LM2,IT,1)
     &                              *FLMSF(IRSF,ISF,IM)
                        END DO
                     END IF
                  END DO
C
                  DO IR = IRMTIN + 1,IRCRIT
                     RL = RPWL(IR)
                     V2(IR) = DRDI(IR,IM)*V1(IR)/R(IR,IM)/RL
                     V1(IR) = DRDI(IR,IM)*V1(IR)*RL
                  END DO
C
C-----------------------------------------------------------------------
C     now integrate V1 and V2
C-----------------------------------------------------------------------
C
C-- use same integration routine as in <CHRDNS> to avoid inconsistencies
C
                  CALL RRADINT_R(IM,V1,VINT1)
C
                  CALL RRADINT_INW_R(IM,V2,VINT2)
C
C-----------------------------------------------------------------------
C     collect all parts
C-----------------------------------------------------------------------
C
                  DO IR = 1,IRCRIT
                     RL = RPWL(IR)
                     VNEW(IR,LM,IT) = FAC*(VINT1(IR)/R(IR,IM)/RL+VINT2(
     &                                IR)*RL)
                  END DO
C
C-----------------------------------------------------------------------
C     store charge moments   CMNTMTT  w.r.t. muffin tin sphere
C                            CMNTIST  w.r.t. interstitial regime
C-----------------------------------------------------------------------
C
                  CMNTMTT(LM,IT) = VINT1(IRMTIN)
                  CMNTIST(LM,IT) = VINT1(IRCRIT) - VINT1(IRMTIN)
C
               END IF
C=======================================================================
C
            END DO LOOP_M
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
         END DO LOOP_L
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      END DO LOOP_IT
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      END
