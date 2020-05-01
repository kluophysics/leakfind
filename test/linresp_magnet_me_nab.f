C*==linresp_magnet_me_nab.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_ME_NAB(IT,ERYD,ZGA,ZFA,ZGB,ZFB,JGA,JFA,
     &                                 JGB,JFB)
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements for             *
C   *                                                                  *
C   *        <E_f, LAM_f | m c ->alfa * ->A_lam | E_i, LAM_i>          *
C   *                                                                  *
C   * according to ELECTRIC DIPOLE selection rules                     *
C   *                                                                  *
C   * using the        NABLA * A - FORM                                *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NCPLWFMAX,CTL
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_LINRESP,ONLY:HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T
      IMPLICIT NONE
C*--LINRESP_MAGNET_ME_NAB23
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IT
      COMPLEX*16 JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),ZFA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      REAL*8 C
      COMPLEX*16 FAP(:,:,:),FBP(:,:,:),GAP(:,:,:),GBP(:,:,:),PF,
     &           WR_FA(:,:,:),WR_FB(:,:,:),WR_GA(:,:,:),WR_GB(:,:,:)
      INTEGER M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FAP,GAP,FBP,GBP
      ALLOCATABLE WR_FA,WR_FB,WR_GA,WR_GB
C
      M = NKMMAX
      ALLOCATE (GAP(NRMAX,NCPLWFMAX,M),FAP(NRMAX,NCPLWFMAX,M))
      ALLOCATE (GBP(NRMAX,NCPLWFMAX,M),FBP(NRMAX,NCPLWFMAX,M))
C
      ALLOCATE (WR_GA(NRMAX,NCPLWFMAX,M),WR_FA(NRMAX,NCPLWFMAX,M))
      ALLOCATE (WR_GB(NRMAX,NCPLWFMAX,M),WR_FB(NRMAX,NCPLWFMAX,M))
C
      C = CTL(IT,1)
C
      PF = 1.0D0/(1.0D0+(ERYD+ERYD)/(2*0.5D0*C**2))
C
      CALL LINRESP_MAGNET_ME_NAB_AUX(IT,PF,ZGA,ZFA,ZGB,ZFB,GAP,FAP,GBP,
     &                               FBP,WR_GA,WR_FA,WR_GB,WR_FB,
     &                               HAZ_ZZ_T)
C
      CALL LINRESP_MAGNET_ME_NAB_AUX(IT,PF,JGA,JFA,ZGB,ZFB,GAP,FAP,GBP,
     &                               FBP,WR_GA,WR_FA,WR_GB,WR_FB,
     &                               HBZ_JZ_T)
C
      CALL LINRESP_MAGNET_ME_NAB_AUX(IT,PF,ZGA,ZFA,JGB,JFB,GAP,FAP,GBP,
     &                               FBP,WR_GA,WR_FA,WR_GB,WR_FB,
     &                               HCZ_ZJ_T)
C
      CALL LINRESP_MAGNET_ME_NAB_AUX(IT,PF,JGA,JFA,JGB,JFB,GAP,FAP,GBP,
     &                               FBP,WR_GA,WR_FA,WR_GB,WR_FB,
     &                               HDZ_JJ_T)
C
      END
C*==linresp_magnet_me_nab_aux.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_ME_NAB_AUX(IT,PF,GA,FA,GB,FB,GAP,FAP,
     &   GBP,FBP,WR_GA,WR_FA,WR_GB,WR_FB,HAZ_T)
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements for             *
C   *                                                                  *
C   *        <E_f, LAM_f | m c ->alfa * ->A_lam | E_i, LAM_i>          *
C   *                                                                  *
C   * according to ELECTRIC DIPOLE selection rules                     *
C   *                                                                  *
C   * using the        NABLA * A - FORM                                *
C   *                                                                  *
C   *                   -- FULL POTENTIAL VERSION --                   *
C   *                                                                  *
C   ********************************************************************
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
      USE MOD_TYPES,ONLY:NCPLWFMAX,IMT,IKMCPLWF,NTMAX
      USE MOD_ANGMOM,ONLY:FPAN1_NAB,FPAN2_NAB,NKMMAX,NCPLWF,NKM,L_IKM,
     &    LB_IKM,IMKM_IKM,MUEM05_IKM
      USE MOD_RMESH,ONLY:JRCRI,NRMAX,NSFMAX,JRCUT,NSF,R,R2DRDI,FLMSF,
     &    LMISF
      USE MOD_LINRESP,ONLY:NPERT
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_4PI
      IMPLICIT NONE
C*--LINRESP_MAGNET_ME_NAB_AUX128
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER IPERT0
      PARAMETER (IPERT0=19)
C
C Dummy arguments
C
      INTEGER IT
      COMPLEX*16 PF
      COMPLEX*16 FA(NRMAX,NCPLWFMAX,NKMMAX),FAP(NRMAX,NCPLWFMAX,NKMMAX),
     &           FB(NRMAX,NCPLWFMAX,NKMMAX),FBP(NRMAX,NCPLWFMAX,NKMMAX),
     &           GA(NRMAX,NCPLWFMAX,NKMMAX),GAP(NRMAX,NCPLWFMAX,NKMMAX),
     &           GB(NRMAX,NCPLWFMAX,NKMMAX),GBP(NRMAX,NCPLWFMAX,NKMMAX),
     &           HAZ_T(NRMAX,NKMMAX,NKMMAX,NTMAX,NPERT),
     &           WR_FA(NRMAX,NCPLWFMAX,NKMMAX),
     &           WR_FB(NRMAX,NCPLWFMAX,NKMMAX),
     &           WR_GA(NRMAX,NCPLWFMAX,NKMMAX),
     &           WR_GB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 AF1,AF2,AF3,AF4,FF1(:),FF2(:),FF3(:),FF4(:),FRAD(:,:),
     &           GG1(:),GG2(:),GG3(:),GG4(:),PAF1,PAF2,PAF3,PAF4,PAG1,
     &           PAG2,PAG3,PAG4
      LOGICAL CNON0
      INTEGER IA2,IB3,IDXPOL(-1:+1),IKMA2,IKMB3,IM,IMKMA2,IMKMB3,IPERT,
     &        IR,IRSF,IRTOP,ISF,JPOL,LA2,LAMA2,LAMB3,LB3,LBARA2,LBARB3,
     &        LM,M,NPOL
      REAL*8 RSCAL(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RSCAL,FRAD
      ALLOCATABLE FF1,FF2,FF3,FF4
      ALLOCATABLE GG1,GG2,GG3,GG4
C
      DATA IDXPOL/2,3,1/
C
      ALLOCATE (RSCAL(NRMAX,NSFMAX),FRAD(NRMAX,3))
      ALLOCATE (FF1(NRMAX),FF2(NRMAX),FF3(NRMAX),FF4(NRMAX))
      ALLOCATE (GG1(NRMAX),GG2(NRMAX),GG3(NRMAX),GG4(NRMAX))
C
      M = NKMMAX
C
      IM = IMT(IT)
      IRTOP = JRCRI(IM)
      NPOL = 3
C
      CALL RINIT(NRMAX*NSFMAX,RSCAL)
C
      DO ISF = 1,NSF(IM)
         LM = LMISF(ISF,IM)
C
         IF ( LM.EQ.1 ) THEN
            DO IR = 1,JRCUT(1,IM)
               RSCAL(IR,ISF) = R2DRDI(IR,IM)*SQRT_4PI
            END DO
         ELSE
            DO IR = 1,JRCUT(1,IM)
               RSCAL(IR,ISF) = 0D0
            END DO
         END IF
C
C-------------- for remaining panels, all shape functions come into play
C
         DO IR = JRCUT(1,IM) + 1,IRTOP
            IRSF = IR - JRCUT(1,IM)
            RSCAL(IR,ISF) = R2DRDI(IR,IM)*FLMSF(IRSF,ISF,IM)
         END DO
C
      END DO
C=======================================================================
C
      DO LAMB3 = 1,NKM
         DO IB3 = 1,NCPLWF(LAMB3)
C
            DO IR = 1,IRTOP
               WR_GB(IR,IB3,LAMB3) = GB(IR,IB3,LAMB3)*RSCAL(IR,1)
               WR_FB(IR,IB3,LAMB3) = FB(IR,IB3,LAMB3)*RSCAL(IR,1)
            END DO
C
            CALL CDIFFER(IM,GB(1,IB3,LAMB3),GBP(1,IB3,LAMB3))
            CALL CDIFFER(IM,FB(1,IB3,LAMB3),FBP(1,IB3,LAMB3))
C
         END DO
      END DO
C
      DO LAMA2 = 1,NKM
         DO IA2 = 1,NCPLWF(LAMA2)
C
            DO IR = 1,IRTOP
               WR_GA(IR,IA2,LAMA2) = GA(IR,IA2,LAMA2)*RSCAL(IR,1)
               WR_FA(IR,IA2,LAMA2) = FA(IR,IA2,LAMA2)*RSCAL(IR,1)
            END DO
C
            CALL CDIFFER(IM,GA(1,IA2,LAMA2),GAP(1,IA2,LAMA2))
            CALL CDIFFER(IM,FA(1,IA2,LAMA2),FAP(1,IA2,LAMA2))
C
         END DO
      END DO
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C                                                      energy B -- LAMB3
      DO LAMB3 = 1,NKM
C
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C                                                      energy A -- LAMA2
         DO LAMA2 = 1,NKM
C
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPERT
C
            DO IB3 = 1,NCPLWF(LAMB3)
               IKMB3 = IKMCPLWF(IB3,LAMB3)
               DO IA2 = 1,NCPLWF(LAMA2)
                  IKMA2 = IKMCPLWF(IA2,LAMA2)
                  DO JPOL = 1,NPOL
                     IF ( CNON0(FPAN1_NAB(IKMA2,IKMB3,1,JPOL,1)) )
     &                    GOTO 20
                     IF ( CNON0(FPAN2_NAB(IKMA2,IKMB3,1,JPOL,1)) )
     &                    GOTO 20
                  END DO
               END DO
            END DO
C
C     ---------------------------------- all angular matrix elements = 0
            CYCLE
C-----------------------------------------------------------------------
 20         CONTINUE
            FRAD(:,:) = C0
C-----------------------------------------------------------------------
C
C                                                                IA2 IB3
            DO IB3 = 1,NCPLWF(LAMB3)
               IKMB3 = IKMCPLWF(IB3,LAMB3)
               IMKMB3 = IMKM_IKM(IKMB3)
               LB3 = L_IKM(IKMB3)
               LBARB3 = LB_IKM(IKMB3)
C
               DO IA2 = 1,NCPLWF(LAMA2)
                  IKMA2 = IKMCPLWF(IA2,LAMA2)
                  IMKMA2 = IMKM_IKM(IKMA2)
                  LA2 = L_IKM(IKMA2)
                  LBARA2 = LB_IKM(IKMA2)
C
C                  DO ISF = 1,NSF(IM)
                  DO ISF = 1,1
                     LM = LMISF(ISF,IM)
C
                     DO IR = 1,IRTOP
C
                        GG1(IR) = WR_GA(IR,IA2,LAMA2)
     &                            *(GBP(IR,IB3,LAMB3)-GB(IR,IB3,LAMB3)
     &                            *LB3/R(IR,IM))
C
                        GG2(IR) = WR_GA(IR,IA2,LAMA2)
     &                            *(GBP(IR,IB3,LAMB3)+GB(IR,IB3,LAMB3)
     &                            *(LB3+1)/R(IR,IM))
C
                        GG3(IR) = DCONJG(WR_GB(IR,IB3,LAMB3)*(GAP(IR,IA2
     &                            ,LAMA2)-GA(IR,IA2,LAMA2)*LA2/R(IR,IM))
     &                            )
C
                        GG4(IR) = DCONJG(WR_GB(IR,IB3,LAMB3)*(GAP(IR,IA2
     &                            ,LAMA2)+GA(IR,IA2,LAMA2)*(LA2+1)
     &                            /R(IR,IM)))
C
                        FF1(IR) = WR_FA(IR,IA2,LAMA2)
     &                            *(FBP(IR,IB3,LAMB3)-FB(IR,IB3,LAMB3)
     &                            *LBARB3/R(IR,IM))
C
                        FF2(IR) = WR_FA(IR,IA2,LAMA2)
     &                            *(FBP(IR,IB3,LAMB3)+FB(IR,IB3,LAMB3)
     &                            *(LBARB3+1)/R(IR,IM))
C
                        FF3(IR) = DCONJG(WR_FB(IR,IB3,LAMB3)*(FAP(IR,IA2
     &                            ,LAMA2)-FA(IR,IA2,LAMA2)
     &                            *LBARA2/R(IR,IM)))
C
                        FF4(IR) = DCONJG(WR_FB(IR,IB3,LAMB3)*(FAP(IR,IA2
     &                            ,LAMA2)+FA(IR,IA2,LAMA2)*(LBARA2+1)
     &                            /R(IR,IM)))
                     END DO
C
                  END DO
C
C -------------------------------------- calculate total matrix elements
C
                  DO M = -1, + 1
C
                     JPOL = IDXPOL(M)
C
C                     IPERT = IPERT0 + JPOL
C
C-------- store results in the sequence -1, 0, +1 as for the other terms
C         NOTE:  may be the sign convention for m=+1 is incoherent
C                i.e. the sign for AG* and AG* has to be reversed
C
                     IPERT = IPERT0 + 2 + M
C
                     IF ( (MUEM05_IKM(IKMA2)-MUEM05_IKM(IKMB3)).EQ.M )
     &                    THEN
C
                        IF ( (IMKMA2.LE.NKM) .AND. (IMKMB3.LE.NKM) )
     &                       THEN
                           AF1 = FPAN1_NAB(IMKMA2,IMKMB3,1,JPOL,1)
                           AF2 = FPAN2_NAB(IMKMA2,IMKMB3,1,JPOL,1)
                           AF3 = FPAN2_NAB(IMKMA2,IMKMB3,1,JPOL,1)
                           AF4 = FPAN1_NAB(IMKMA2,IMKMB3,1,JPOL,1)
                        ELSE
                           AF1 = 0.0D0
                           AF2 = 0.0D0
                           AF3 = 0.0D0
                           AF4 = 0.0D0
                        END IF
C
                        PAG1 = PF*(1.0D0/CI)
     &                         *0.5D0*FPAN1_NAB(IKMA2,IKMB3,1,JPOL,1)
                        PAG2 = PF*(1.0D0/CI)
     &                         *0.5D0*FPAN2_NAB(IKMA2,IKMB3,1,JPOL,1)
                        PAG3 = DCONJG(PF*(1.0D0/CI))
     &                         *0.5D0*FPAN2_NAB(IKMA2,IKMB3,1,JPOL,1)
                        PAG4 = DCONJG(PF*(1.0D0/CI))
     &                         *0.5D0*FPAN1_NAB(IKMA2,IKMB3,1,JPOL,1)
C
                        PAF1 = PF*(1.0D0/CI)*0.5D0*AF1
                        PAF2 = PF*(1.0D0/CI)*0.5D0*AF2
                        PAF3 = DCONJG(PF*(1.0D0/CI))*0.5D0*AF3
                        PAF4 = DCONJG(PF*(1.0D0/CI))*0.5D0*AF4
C
                        DO IR = 1,IRTOP
                           FRAD(IR,JPOL) = FRAD(IR,JPOL)
     &                        + (GG1(IR)*PAG1-GG2(IR)*PAG2+FF1(IR)
     &                        *PAF1-FF2(IR)
     &                        *PAF2+DCONJG(GG3(IR)*PAG3-GG4(IR)
     &                        *PAG4+FF3(IR)*PAF3-FF4(IR)*PAF4))
                        END DO
C
                     END IF
C
                  END DO
C
               END DO
            END DO
C                                                                IA2 IB3
C-----------------------------------------------------------------------
C
C#######################################################################
            JPOL = 1
            IPERT = IPERT0 + 3
C
            CALL CRADINT_R(IM,FRAD(1,JPOL),HAZ_T(1,LAMA2,LAMB3,IT,IPERT)
     &                     )
C
            DO JPOL = 2,3
               IPERT = IPERT0 + JPOL - 1
C
               CALL CRADINT_R(IM,FRAD(1,JPOL),
     &                        HAZ_T(1,LAMA2,LAMB3,IT,IPERT))
C
            END DO
C#######################################################################
C
         END DO
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      DEALLOCATE (RSCAL)
      END
