C*==linresp_magnet_me_alf.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_MAGNET_ME_ALF(IFIL_LHSA,IFIL_RHSB,ERYD,IT,ZGA,
     &                                 ZFA,ZGB,ZFB,JGA,JFA,JGB,JFB)
C   ********************************************************************
C   *                                                                  *
C   *    read wave function and the calculate                          *
C   *    the matrix element integral functions                         *
C   *                                                                  *
C   *      IXAYB(r) = < X^+(E_a) | H_lam | Y(E_b) >_r                  *
C   *                                                                  *
C   *  X, Y: = Z (J) regular (irregular) solutions to Dirac equation   *
C   *                                                                  *
C   *    polarisation lam:        ipol = 1,2,3  ==  (-),(0),(+)        *
C   *                                                                  *
C   ********************************************************************
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
      USE MOD_TYPES,ONLY:IMT,CTL,IKMCPLWF,NCPLWFMAX
      USE MOD_RMESH,ONLY:R,NRMAX,JRCRI,R2DRDI
      USE MOD_CONSTANTS,ONLY:CI,C0,SQRT_4PI
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NCPLWF
      USE MOD_LINRESP,ONLY:SHAZ_ZZ_T,SHBZ_JZ_T,SHCZ_ZJ_T,SHDZ_JJ_T,
     &    HAZ_ZZ_T,HBZ_JZ_T,HCZ_ZJ_T,HDZ_JJ_T,AMEG_LMOP,AMEF_LMOP
      IMPLICIT NONE
C*--LINRESP_MAGNET_ME_ALF38
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER IPERT1,IPERT2,N_TERMS_A_AND_B
      PARAMETER (IPERT1=8,IPERT2=19,N_TERMS_A_AND_B=9)
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IFIL_LHSA,IFIL_RHSB,IT
      COMPLEX*16 JFA(NRMAX,NCPLWFMAX,NKMMAX),JFB(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,JGA(NRMAX,NCPLWFMAX,NKMMAX),
     &           JGB(NRMAX,NCPLWFMAX,NKMMAX),ZFA(NRMAX,NCPLWFMAX,NKMMAX)
     &           ,ZFB(NRMAX,NCPLWFMAX,NKMMAX),
     &           ZGA(NRMAX,NCPLWFMAX,NKMMAX),ZGB(NRMAX,NCPLWFMAX,NKMMAX)
C
C Local variables
C
      COMPLEX*16 A1,A1JGA,A1ZGA,A2,A2JFA,A2ZFA,CARG,FJAJB(:),FJAZB(:),
     &           FZAJB(:),FZAZB(:),IMC,IMCA1,IMCA2,PA,PAB
      REAL*8 C,RTOP
      LOGICAL CNON0
      INTEGER IA2,IB3,IFLAG1,IKMA2,IKMB3,IM,IPERT,IR,IRTOP,LAMA2,LAMB3
      LOGICAL TERMS_A_AND_B
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FZAZB,FZAJB,FJAZB,FJAJB
C
      CNON0(CARG) = ABS(DREAL(CARG)) + ABS(DIMAG(CARG)).GT.NON0_TOL
C
      ALLOCATE (FZAZB(NRMAX),FZAJB(NRMAX),FJAZB(NRMAX),FJAJB(NRMAX))
C
      SHAZ_ZZ_T(:,:,:,:) = C0
      SHCZ_ZJ_T(:,:,:,:) = C0
      SHBZ_JZ_T(:,:,:,:) = C0
      SHDZ_JJ_T(:,:,:,:) = C0
C
C=======================================================================
C
      C = CTL(IT,1)
      IM = IMT(IT)
      IRTOP = JRCRI(IM)
      RTOP = R(IRTOP,IM)
C
      IMC = CI*0.5D0*C
C
      PAB = (C/(C**2+ERYD+ERYD))*RTOP**2
C
C=======================================================================
C ------------------------------------------ read in wavefunctions for A
C                                                use arrays *B as buffer
C
      CALL WAVFUN_READ_REL(IFIL_LHSA,IT,1,ZGB,ZFB,JGB,JFB,IRTOP,NCPLWF,
     &                     IKMCPLWF)
C
C------------------ convolute the wave function with the shape functions
C------------ and multiply with the radial factor R2DRDI for integration
C
      CALL WAVFUN_WGT_RAD_FLM_REL(IT,ZGB,ZFB,JGB,JFB,ZGA,ZFA,JGA,JFA,
     &                            NCPLWF,IKMCPLWF)
C
C ------------------------------------------ read in wavefunctions for B
C
      CALL WAVFUN_READ_REL(IFIL_RHSB,IT,1,ZGB,ZFB,JGB,JFB,IRTOP,NCPLWF,
     &                     IKMCPLWF)
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
            DO IPERT = IPERT1,IPERT2
C
               IF ( IPERT.LE.(IPERT1+N_TERMS_A_AND_B-1) ) THEN
                  TERMS_A_AND_B = .TRUE.
               ELSE
                  TERMS_A_AND_B = .FALSE.
               END IF
C
               DO IR = 1,IRTOP
                  FZAZB(IR) = C0
                  FZAJB(IR) = C0
                  FJAZB(IR) = C0
                  FJAJB(IR) = C0
               END DO
C
               IFLAG1 = 0
C
C-----------------------------------------------------------------------
C                                                                IA2 IB3
               DO IA2 = 1,NCPLWF(LAMA2)
                  IKMA2 = IKMCPLWF(IA2,LAMA2)
C
                  DO IB3 = 1,NCPLWF(LAMB3)
                     IKMB3 = IKMCPLWF(IB3,LAMB3)
C
                     A1 = AMEG_LMOP(IKMA2,IKMB3,1,IPERT)*SQRT_4PI
                     A2 = AMEF_LMOP(IKMA2,IKMB3,1,IPERT)*SQRT_4PI
C
                     IF ( CNON0(A1) .OR. CNON0(A2) ) THEN
C
                        IFLAG1 = 1
C
                        IMCA1 = IMC*A1
                        IMCA2 = IMC*A2
C
                        DO IR = 1,IRTOP
                           A1ZGA = IMCA1*ZGA(IR,IA2,LAMA2)
                           A2ZFA = IMCA2*ZFA(IR,IA2,LAMA2)
                           A1JGA = IMCA1*JGA(IR,IA2,LAMA2)
                           A2JFA = IMCA2*JFA(IR,IA2,LAMA2)
C
                           FZAZB(IR) = FZAZB(IR)
     &                                 + A1ZGA*ZFB(IR,IB3,LAMB3)
     &                                 - A2ZFA*ZGB(IR,IB3,LAMB3)
                           FZAJB(IR) = FZAJB(IR)
     &                                 + A1ZGA*JFB(IR,IB3,LAMB3)
     &                                 - A2ZFA*JGB(IR,IB3,LAMB3)
C
                           FJAZB(IR) = FJAZB(IR)
     &                                 + A1JGA*ZFB(IR,IB3,LAMB3)
     &                                 - A2JFA*ZGB(IR,IB3,LAMB3)
                           FJAJB(IR) = FJAJB(IR)
     &                                 + A1JGA*JFB(IR,IB3,LAMB3)
     &                                 - A2JFA*JGB(IR,IB3,LAMB3)
C
                        END DO
C
                        IR = IRTOP
                        PA = -0.5D0*PAB*(IMCA1-IMCA2)/R2DRDI(IR,IM)
C
                        SHAZ_ZZ_T(LAMA2,LAMB3,IT,IPERT)
     &                     = SHAZ_ZZ_T(LAMA2,LAMB3,IT,IPERT)
     &                     + (ZGA(IR,IA2,LAMA2)*ZGB(IR,IB3,LAMB3)
     &                     -ZFA(IR,IA2,LAMA2)*ZFB(IR,IB3,LAMB3))*PA
C
                        SHCZ_ZJ_T(LAMA2,LAMB3,IT,IPERT)
     &                     = SHCZ_ZJ_T(LAMA2,LAMB3,IT,IPERT)
     &                     + (ZGA(IR,IA2,LAMA2)*JGB(IR,IB3,LAMB3)
     &                     -ZFA(IR,IA2,LAMA2)*JFB(IR,IB3,LAMB3))*PA
C
                        SHBZ_JZ_T(LAMA2,LAMB3,IT,IPERT)
     &                     = SHBZ_JZ_T(LAMA2,LAMB3,IT,IPERT)
     &                     + (JGA(IR,IA2,LAMA2)*ZGB(IR,IB3,LAMB3)
     &                     -JFA(IR,IA2,LAMA2)*ZFB(IR,IB3,LAMB3))*PA
C
                        SHDZ_JJ_T(LAMA2,LAMB3,IT,IPERT)
     &                     = SHDZ_JJ_T(LAMA2,LAMB3,IT,IPERT)
     &                     + (JGA(IR,IA2,LAMA2)*JGB(IR,IB3,LAMB3)
     &                     -JFA(IR,IA2,LAMA2)*JFB(IR,IB3,LAMB3))*PA
C
C
                     END IF
                  END DO
               END DO
C                                                                IA2 IB3
C-----------------------------------------------------------------------
C
C#######################################################################
C
               IF ( IFLAG1.EQ.1 ) THEN
C
C---------------------- multiply with missing factor R for terms A and B
C
                  IF ( TERMS_A_AND_B ) THEN
                     DO IR = 1,IRTOP
                        FZAZB(IR) = FZAZB(IR)*R(IR,IM)
                        FZAJB(IR) = FZAJB(IR)*R(IR,IM)
                        FJAZB(IR) = FJAZB(IR)*R(IR,IM)
                        FJAJB(IR) = FJAJB(IR)*R(IR,IM)
                     END DO
                  END IF
C-----------------------------------------------------------------------
C
                  CALL CRADINT_R(IM,FZAZB,
     &                           HAZ_ZZ_T(1,LAMA2,LAMB3,IT,IPERT))
C
                  CALL CRADINT_INW_R(IM,FZAJB,
     &                               HCZ_ZJ_T(1,LAMA2,LAMB3,IT,IPERT))
C
                  CALL CRADINT_INW_R(IM,FJAZB,
     &                               HBZ_JZ_T(1,LAMA2,LAMB3,IT,IPERT))
C
                  CALL CRADINT_INW_R(IM,FJAJB,
     &                               HDZ_JJ_T(1,LAMA2,LAMB3,IT,IPERT))
C
               END IF
C################################################################ IFLAG1
C
            END DO
CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP IPERT
C
         END DO
C                                                      energy A -- LAMA2
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      END DO
C                                                      energy B -- LAMB3
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
      END
