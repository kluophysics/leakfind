C*==linresp_susc_stat.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_SUSC_STAT
C   ********************************************************************
C   *                                                                  *
C   *   calculate STATIC susceptibility                                *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1        1 IDOS                           *
C   *  spin moment           b s_z    2 ISPN                           *
C   *  orbital moment        b l_z    3 IORB   ____ NPERT              *
C   *  hyperfine interaction H_hf,z   4 IHFI                           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:PHASK
      USE MOD_RMESH,ONLY:NRMAX,NMMAX,NMAX_LEBGRID,N_LEBGRID,W_LEBGRID,
     &    Y_LEBGRID
      USE MOD_FILES,ONLY:IPRINT,IOTMP,FOUND_SECTION
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,IND0Q,NKKR,MEZJ,MEZZ,TSSQ,MSSQ,
     &    MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQMAX,NQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:NTMAX,NLFP,NLMFPMAX,NT,CONC,NAT,RHO2NS
      USE MOD_MPI,ONLY:MPI
      USE MOD_CALCMODE,ONLY:NONMAG
      USE MOD_LINRESP,ONLY:HX_PERT_LMTP,HZ_PERT_LMTP,RHO2NS_GG,NOBSE,
     &    T0Z,T1Z,TIJZ,TZ,T0X,TIJX,IDOS,ISPN,IORB,IHFI,DOBS_TO,
     &    K_PERT_LMTP,CHI_TO,QVEC_PERT_EQ_0VEC,QVEC_PERT,NZ12MAX,ITTA,
     &    ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,JTTX,WTTJ,ITTMAX,
     &    JTTMAX,NTKTKMAX,NTKTKLIN,CHIZ
      USE MOD_SCF,ONLY:SCFVXC,RMSAVB,RMSAVV,SCFMIX,SCFTOL,NSCFITER
      IMPLICIT NONE
C*--LINRESP_SUSC_STAT33
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_SUSC_STAT')
      LOGICAL HYPERFINE_INTERACTION
      PARAMETER (HYPERFINE_INTERACTION=.FALSE.)
C
C Local variables
C
      REAL*8 AA(:,:),AAHF(:,:),BB(:,:),BBHF_N(:),BBHF_O(:),BBHF_S(:),
     &       BCP(:),BCPS(:),CHI_N(:),CHI_NN(:),CHI_NN0(:),CHI_O(:),
     &       CHI_OO(:),CHI_OO0(:),CHI_OO_TR(:,:),CHI_OS(:),CHI_OS0(:),
     &       CHI_S(:),CHI_SO(:),CHI_SO0(:),CHI_SS(:),CHI_SS0(:),
     &       CHI_SS_TR(:,:),DOBSEF(:,:),DOSEF,DQVEC_PERT,ICHI_O(:),
     &       ICHI_S(:),POBS(:,:),QHAT_PERT(3),QPERT,QVEC_PERT0,
     &       RGNT_VSF(:),RHO2NSX(:,:,:,:),SFT_O(:),SFT_O0(:),
     &       SFT_O_TR(:,:),SFT_S(:),SFT_S0(:),SFT_S_TR(:,:),SUM_WTZ_NN,
     &       SUM_WTZ_NO,SUM_WTZ_NS,TCDIACRI(NTMAX),TCDIARI(NTMAX),TIME1,
     &       TIME2,TKDIACRI(NTMAX),TKDIARI(NTMAX),TOBS(:,:),W1MIX(:,:),
     &       W1NSMIX(:,:,:),W2MIX(:,:),W2NSMIX(:,:,:),WIHF,WIN,WIO,WIS,
     &       WJ,WJP,WTX_NN(NTMAX),WTX_NO(NTMAX),WTX_NS(NTMAX)
      REAL*8 DNRM2
      INTEGER I,IA_ERR,IBOTN,IBOTO,IBOTS,IMIX,INFO,IOBSE,IPIV(:),IQ,
     &        IQVEC_PERT,IQVEC_PERT1,IQVEC_PERT2,ISOCOUPLING,ISTBRY,IT,
     &        ITDEPT,ITRSCF,J,JQ,JRNS1TMP(:),JRNSMINTMP,JT,JTP,L1,L2,L3,
     &        L4,LMRGNT_VSF(:,:),M,NDIM,NOBSEDNS,NRGNT_VSF,
     &        NRGNT_VSF_LM(:)
      LOGICAL SPIN_POLARIZED
      COMPLEX*16 SUMEF,WT0X_NN(NTMAX),WT0X_NO(NTMAX),WT0X_NS(NTMAX),
     &           WTIJX_NN(NTMAX),WTIJX_NO(NTMAX),WTIJX_NS(NTMAX),
     &           WTZ_NN(NTMAX),WTZ_NO(NTMAX),WTZ_NS(NTMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1MIX,W2MIX,W1NSMIX,W2NSMIX,JRNS1TMP
      ALLOCATABLE RHO2NSX
      ALLOCATABLE BCP,BCPS
      ALLOCATABLE LMRGNT_VSF,NRGNT_VSF_LM,RGNT_VSF
C
      ALLOCATABLE AA,BB,AAHF
      ALLOCATABLE CHI_N,CHI_O,CHI_S,POBS,IPIV,SFT_O
      ALLOCATABLE SFT_S,SFT_O0,SFT_S0,BBHF_N,BBHF_O,BBHF_S,ICHI_O,CHI_NN
      ALLOCATABLE CHI_OO,CHI_OS,CHI_SO,CHI_SS,ICHI_S
      ALLOCATABLE CHI_NN0,CHI_OO0
      ALLOCATABLE CHI_OS0,CHI_SO0,CHI_SS0,DOBSEF
      ALLOCATABLE SFT_O_TR,SFT_S_TR
      ALLOCATABLE CHI_OO_TR,CHI_SS_TR
      ALLOCATABLE TOBS
C
C ======================================================================
C
C   CHI vector is defined as follows
C   CHI = ( CHI_SPIN(NT), CHI_ORB(NT), CHI_CHARGE(NT) )
C       = ( CHI(IBOTS+1),...,CHI(IBOTS+NT),
C           CHI(IBOTO+1),...,CHI(IBOTO+NT),
C           CHI(IBOTN+1),...,CHI(IBOTN+NT) )
C
      IBOTS = 0
      IBOTO = NT
      IBOTN = NT*2
C
      CALL LINRESP_INIT('SUSC_STAT ')
C
C ======================================================================
C
      NOBSEDNS = NOBSE + 3
      ALLOCATE (TOBS(NTMAX,NOBSEDNS))
C
C----------------------------------------------------- variables for BXC
C
      ALLOCATE (JRNS1TMP(NMMAX))
      JRNSMINTMP = 1
      JRNS1TMP(1:NMMAX) = 1
      ALLOCATE (W1NSMIX(JRNSMINTMP:NRMAX,NLMFPMAX,NTMAX))
      ALLOCATE (W2NSMIX(JRNSMINTMP:NRMAX,NLMFPMAX,NTMAX))
C
      ALLOCATE (W1MIX(NRMAX,NTMAX),W2MIX(NRMAX,NTMAX))
      SCFVXC = 'VWN'
C
C ------- initialize local variables  RGNT_VSF, NRGNT_VSF_LM, LMRGNT_VSF
C
      ALLOCATE (NRGNT_VSF_LM(0:NLMFPMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate RGNT_VSF')
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
      CALL RGNTSF_SETUP(NLFP-1,NLFP-1,NRGNT_VSF_LM,NRGNT_VSF)
C
      ALLOCATE (RGNT_VSF(NRGNT_VSF),LMRGNT_VSF(NRGNT_VSF,3))
C
      REWIND IOTMP
      DO I = 1,NRGNT_VSF
         READ (IOTMP) (LMRGNT_VSF(I,J),J=1,3),RGNT_VSF(I)
      END DO
      CLOSE (IOTMP)
C
C----------------------------------------------------- variables for CHI
C
      ALLOCATE (AA(NTMAX*3,NTMAX*3),BB(NTMAX*3,1))
      ALLOCATE (AAHF(NTMAX*3,NTMAX*3),BBHF_S(NTMAX*3))
      ALLOCATE (BBHF_N(3*NTMAX),BBHF_O(NTMAX*3))
C
      ALLOCATE (CHI_S(NTMAX),CHI_SS0(NTMAX),CHI_SS(NTMAX))
      ALLOCATE (CHI_O(NTMAX),CHI_OO0(NTMAX),CHI_OO(NTMAX))
      ALLOCATE (CHI_N(NTMAX),CHI_NN0(NTMAX))
      ALLOCATE (CHI_SO(NTMAX),CHI_SO0(NTMAX))
      ALLOCATE (CHI_OS(NTMAX),CHI_OS0(NTMAX))
      ALLOCATE (CHI_OO_TR(NTMAX,NTMAX))
      ALLOCATE (CHI_SS_TR(NTMAX,NTMAX))
      ALLOCATE (SFT_S(NTMAX),SFT_S0(NTMAX))
      ALLOCATE (SFT_O(NTMAX),SFT_O0(NTMAX))
      ALLOCATE (SFT_O_TR(NTMAX,NTMAX))
      ALLOCATE (SFT_S_TR(NTMAX,NTMAX))
      ALLOCATE (ICHI_S(NTMAX),ICHI_O(NTMAX),CHI_NN(NTMAX))
C
      ALLOCATE (POBS(NTMAX,NOBSE),IPIV(3*NTMAX))
      ALLOCATE (DOBSEF(NTMAX,NOBSEDNS))
C
      ALLOCATE (RHO2NSX(NRMAX,NLMFPMAX,NTMAX,3))
      ALLOCATE (BCP(NTMAX),BCPS(NTMAX))
      BCP(:) = 0D0
      BCPS(:) = 0D0
      TCDIACRI(:) = 0D0
      TCDIARI(:) = 0D0
      TKDIACRI(:) = 0D0
      TKDIARI(:) = 0D0
      TOBS(:,:) = 0D0
      CHI_NN(:) = 0D0
      CHI_NN0(:) = 0D0
C
      NONMAG = .FALSE.
C
      IF ( .NOT.NONMAG ) SPIN_POLARIZED = .TRUE.
C
C
      WRITE (6,99001)
C
C=======================================================================
C
      CALL CPU_TIME(TIME1)
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C ======================================================================
C               read variables for CHI SCF-CYCLE
C ======================================================================
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
C------------------------------------------- standard Broyden parameters
      IMIX = 4
C
      CALL SECTION_SET_INTEGER('NITER',NSCFITER,1,0)
      CALL SECTION_SET_INTEGER('ISTBRY',ISTBRY,50,0)
      CALL SECTION_SET_INTEGER('ITDEPT',ITDEPT,40,0)
      CALL SECTION_SET_REAL('MIX',SCFMIX,0.20D0,0)
      CALL SECTION_SET_REAL('TOL',SCFTOL,1D-5,0)
C
      NSCFITER = 2
      IF ( NSCFITER.GT.1 ) THEN
         WRITE (6,99003) SCFTOL,SCFMIX,NSCFITER
      ELSE
         WRITE (6,99004)
      END IF
C
C ======================================================================
C                        q_pert - vectors
C ======================================================================
C
      IQVEC_PERT1 = 0
      IQVEC_PERT2 = 0
      QVEC_PERT0 = 0.0D0
      QHAT_PERT(1) = 1.0D0
      QHAT_PERT(2) = 0.0D0
      QHAT_PERT(3) = 0.0D0
      QPERT = DNRM2(3,QHAT_PERT,1)
      QHAT_PERT(1:3) = QHAT_PERT(1:3)/QPERT
      DQVEC_PERT = QPERT/DBLE(IQVEC_PERT2-1)
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_INTEGER('IPERT1',IQVEC_PERT1,9999,0)
         CALL SECTION_SET_INTEGER('IPERT2',IQVEC_PERT2,9999,0)
         CALL SECTION_SET_REAL('QPERT0',QVEC_PERT0,9999D0,0)
         CALL SECTION_SET_REAL('DQPERT',DQVEC_PERT,9999D0,0)
      END IF
C
C
      IF ( IQVEC_PERT1.EQ.0 .AND. IQVEC_PERT2.EQ.0 ) THEN
C
         QVEC_PERT_EQ_0VEC = .TRUE.
C
         QVEC_PERT0 = 0D0
         DQVEC_PERT = 0D0
         QHAT_PERT(1:3) = 0D0
C
      ELSE
C
         QVEC_PERT_EQ_0VEC = .FALSE.
C
C ======================================================================
C         set up index table for BZ integration for NO SYMMETRY case
C ======================================================================
C
         IF ( ALLOCATED(WTTJ) ) THEN
            DEALLOCATE (WTTJ,JTT1,JTT2,JTTX)
            DEALLOCATE (NTTJ,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2)
         END IF
C
C------------------------------- consider ALL TAU(k)*TAU(k) combinations
         NTKTKMAX = NKM**4
C
         ITTMAX = NTKTKMAX*NQ*NQ
C
         ALLOCATE (NTTJ(ITTMAX))
         ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX))
         ALLOCATE (ITTD(ITTMAX),ITTQ1(ITTMAX),ITTQ2(ITTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ITTD')
C
         JTTMAX = ITTMAX
C
         ALLOCATE (WTTJ(JTTMAX))
         ALLOCATE (JTT1(JTTMAX),JTT2(JTTMAX),JTTX(JTTMAX))
C
         NTKTKLIN = ITTMAX
C
         I = 0
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO JQ = IQBOT_CHI,IQTOP_CHI
               DO L1 = 1,NKM
                  DO L4 = 1,NKM
                     DO L2 = 1,NKM
                        DO L3 = 1,NKM
                           I = I + 1
                           ITTA(I) = L1
                           ITTB(I) = L2
                           ITTC(I) = L3
                           ITTD(I) = L4
                           ITTQ1(I) = IQ
                           ITTQ2(I) = JQ
                           NTTJ(I) = 1
C
                           J = I
                           JTT1(J) = (IND0Q(JQ)+L2-1)*NKKR + IND0Q(IQ)
     &                               + L1
                           JTT2(J) = (IND0Q(IQ)+L4-1)*NKKR + IND0Q(JQ)
     &                               + L3
                           JTTX(J) = (IND0Q(JQ)+L3-1)*NKKR + IND0Q(IQ)
     &                               + L4
                           WTTJ(J) = 1D0
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C ======================================================================
C                         re-allocate CHIZ
C ======================================================================
C
         DEALLOCATE (CHIZ)
C
         NZ12MAX = 2
         M = NKMMAX
C
         ALLOCATE (CHIZ(M*M*NQMAX,M*M*NQMAX,NZ12MAX),STAT=IA_ERR)
C
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate CHIZ')
C
C
C ======================================================================
C                   initialize k-mesh for FULL BZ
C ======================================================================
C
         CALL INIT_MOD_TAUIJ_KMESH
C
      END IF
C
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                        loop over q_pert - vectors
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IQVEC_PERT = IQVEC_PERT1,IQVEC_PERT2
C
         QVEC_PERT(1:3) = QVEC_PERT0 + DQVEC_PERT*IQVEC_PERT*QHAT_PERT
     &                    (1:3)
         WRITE (6,99016) IQVEC_PERT,QVEC_PERT
C
C
C=======================================================================
C                      CHI - SCF - cycle  START
C=======================================================================
C
C------------------------------------------ initialize for 1st iteration
C
         HZ_PERT_LMTP(:,:,:,2) = HZ_PERT_LMTP(:,:,:,1)
         HZ_PERT_LMTP(:,:,:,3) = HZ_PERT_LMTP(:,:,:,1)
         K_PERT_LMTP(:,:,2) = K_PERT_LMTP(:,:,1)
         K_PERT_LMTP(:,:,3) = K_PERT_LMTP(:,:,1)
C
         HX_PERT_LMTP(:,:,:,:) = 0D0
C     HX_PERT_LMTP(:,:,:,:) = HZ_PERT_LMTP(:,:,:,:)
C
         CHI_TO(:,:) = 0D0
         CHI_TO(:,:) = 0D0
         CHI_TO(:,:) = 0D0
C
         DO ITRSCF = 1,NSCFITER
C
            CALL LINRESP_ELOOP(MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,
     &                         TSST,RHO2NSX,PHASK)
C
C ======================================================================
C              terms connected with high-field susceptibility
C ======================================================================
C
C--------------------------------------- last E-mesh point: Fermi energy
            DO IT = 1,NT
               DO IOBSE = 1,NOBSE
                  DOBSEF(IT,IOBSE) = DIMAG(DOBS_TO(IT,IOBSE))
               END DO
            END DO
C
            DOSEF = 0D0
            DO IT = 1,NT
               DOSEF = DOSEF + NAT(IT)*CONC(IT)*DOBSEF(IT,IDOS)
            END DO
C
            DO IOBSE = 1,NOBSE
               DO IT = 1,NT
                  POBS(IT,IOBSE) = DOBSEF(IT,IOBSE)/DOSEF
               END DO
            END DO
C
C-----------------------------------------------------------------------
            SUM_WTZ_NS = 0D0
            SUM_WTZ_NO = 0D0
            SUM_WTZ_NN = 0D0
C
            DO JT = 1,NT
C
               WJ = NAT(JT)*CONC(JT)
C
C---------- There is no charge coupling with external magnetic field
C---------- Therefore WTZ_NN(JT) = WJ*TZ(JT,1,IDOS) = 0
C
               WTZ_NS(JT) = WJ*TZ(JT,IDOS,ISPN)
               WTZ_NO(JT) = WJ*TZ(JT,IDOS,IORB)
               WTZ_NN(JT) = 0.D0
C
               SUM_WTZ_NS = SUM_WTZ_NS + DIMAG(WTZ_NS(JT))
               SUM_WTZ_NO = SUM_WTZ_NO + DIMAG(WTZ_NO(JT))
               SUM_WTZ_NN = SUM_WTZ_NN + DIMAG(WTZ_NN(JT))
C
               WT0X_NS(JT) = WJ*T0X(JT,IDOS,ISPN)
               WT0X_NO(JT) = WJ*T0X(JT,IDOS,IORB)
               WT0X_NN(JT) = WJ*T0X(JT,IDOS,IDOS)
C
               WTIJX_NS(JT) = 0D0
               WTIJX_NO(JT) = 0D0
               WTIJX_NN(JT) = 0D0
C
               DO JTP = 1,NT
C
                  WJP = NAT(JTP)*CONC(JTP)
C
                  WTIJX_NS(JT) = WTIJX_NS(JT)
     &                           + WJP*TIJX(JTP,JT,IDOS,ISPN)
                  WTIJX_NO(JT) = WTIJX_NO(JT)
     &                           + WJP*TIJX(JTP,JT,IDOS,IORB)
                  WTIJX_NN(JT) = WTIJX_NN(JT)
     &                           + WJP*TIJX(JTP,JT,IDOS,IDOS)
C
               END DO
C
               WTX_NS(JT) = DIMAG(WT0X_NS(JT)+WTIJX_NS(JT))
               WTX_NO(JT) = DIMAG(WT0X_NO(JT)+WTIJX_NO(JT))
               WTX_NN(JT) = DIMAG(WT0X_NN(JT)+WTIJX_NN(JT))
C
            END DO
C
C---------- suppress terms connected with particle number N if requested
C
            IF ( .NOT.SPIN_POLARIZED ) THEN
C
               SUM_WTZ_NN = 0D0
               WTX_NN(:) = 0D0
C
               NDIM = 2*NT
C
            ELSE IF ( NT.EQ.1 ) THEN
C
               NDIM = 2
C
            ELSE
C
               NDIM = 3*NT
C
            END IF
C
            AAHF(:,:) = 0D0
C
            DO IT = 1,NT
C
C----contribution to high-field suscept due to shift of the Fermi energy
C
               WIS = POBS(IT,ISPN)
               WIO = POBS(IT,IORB)
               WIN = POBS(IT,IDOS)
C
               BBHF_S(IBOTS+IT) = WIS*SUM_WTZ_NS
               BBHF_O(IBOTS+IT) = WIS*SUM_WTZ_NO
               BBHF_N(IBOTS+IT) = WIS*SUM_WTZ_NN
C
               DO JT = 1,NT
                  AAHF(IBOTS+IT,IBOTS+JT) = -WIS*WTX_NS(JT)
                  AAHF(IBOTS+IT,IBOTO+JT) = -WIS*WTX_NO(JT)
C
                  IF ( NT.EQ.1 ) THEN
                     BBHF_N(IBOTS+IT) = BBHF_N(IBOTS+IT)
     &                                  + WIS*WTX_NN(JT)
                  ELSE
                     AAHF(IBOTS+IT,IBOTN+JT) = -WIS*WTX_NN(JT)
                  END IF
               END DO
C
               BBHF_S(IBOTO+IT) = WIO*SUM_WTZ_NS
               BBHF_O(IBOTO+IT) = WIO*SUM_WTZ_NO
               BBHF_N(IBOTO+IT) = WIO*SUM_WTZ_NN
               DO JT = 1,NT
                  AAHF(IBOTO+IT,IBOTS+JT) = -WIO*WTX_NS(JT)
                  AAHF(IBOTO+IT,IBOTO+JT) = -WIO*WTX_NO(JT)
C
                  IF ( NT.EQ.1 ) THEN
                     BBHF_N(IBOTO+IT) = BBHF_N(IBOTO+IT)
     &                                  + WIO*WTX_NN(JT)
                  ELSE
                     AAHF(IBOTO+IT,IBOTN+JT) = -WIO*WTX_NN(JT)
                  END IF
               END DO
C
               BBHF_S(IBOTN+IT) = WIN*SUM_WTZ_NS
               BBHF_O(IBOTN+IT) = WIN*SUM_WTZ_NO
               BBHF_N(IBOTN+IT) = WIN*SUM_WTZ_NN
               DO JT = 1,NT
                  AAHF(IBOTN+IT,IBOTS+JT) = -WIN*WTX_NS(JT)
                  AAHF(IBOTN+IT,IBOTO+JT) = -WIN*WTX_NO(JT)
C
                  IF ( NT.EQ.1 ) THEN
                     BBHF_N(IBOTN+IT) = BBHF_N(IBOTN+IT)
     &                                  + WIN*WTX_NN(JT)
                  ELSE
                     AAHF(IBOTN+IT,IBOTN+JT) = -WIN*WTX_NN(JT)
                  END IF
               END DO
C
C-------------- the terms due to spin-charge and orbital-charge coupling
C
               BBHF_S(IBOTN+IT) = BBHF_S(IBOTN+IT)
     &                            - DIMAG(TZ(IT,IDOS,ISPN))
               BBHF_O(IBOTN+IT) = BBHF_O(IBOTN+IT)
     &                            - DIMAG(TZ(IT,IDOS,IORB))
               AAHF(IBOTN+IT,IBOTN+IT) = AAHF(IBOTN+IT,IBOTN+IT)
     &            - DIMAG(1D0-T0X(IT,IDOS,IDOS))
               DO JT = 1,NT
                  AAHF(IBOTN+IT,IBOTN+JT) = AAHF(IBOTN+IT,IBOTN+JT)
     &               + DIMAG(TIJX(IT,JT,IDOS,IDOS))
               END DO
C
               IF ( NT.EQ.1 ) THEN
C
                  BBHF_N(IBOTS+IT) = BBHF_N(IBOTS+IT)
     &                               - DIMAG(T0X(IT,ISPN,IDOS))
                  DO JT = 1,NT
                     BBHF_N(IBOTS+IT) = BBHF_N(IBOTS+IT)
     &                                  - DIMAG(TIJX(IT,JT,ISPN,IDOS))
                  END DO
                  BBHF_N(IBOTO+IT) = BBHF_N(IBOTO+IT)
     &                               - DIMAG(T0X(IT,IORB,IDOS))
C
                  DO JT = 1,NT
                     BBHF_N(IBOTO+IT) = BBHF_N(IBOTO+IT)
     &                                  - DIMAG(TIJX(IT,JT,IORB,IDOS))
                  END DO
C
               ELSE
C
                  AAHF(IBOTS+IT,IBOTN+IT) = AAHF(IBOTS+IT,IBOTN+IT)
     &               + DIMAG(T0X(IT,ISPN,IDOS))
                  AAHF(IBOTN+IT,IBOTS+IT) = AAHF(IBOTN+IT,IBOTS+IT)
     &               + DIMAG(T0X(IT,IDOS,ISPN))
                  DO JT = 1,NT
                     AAHF(IBOTS+IT,IBOTN+JT) = AAHF(IBOTS+IT,IBOTN+JT)
     &                  + DIMAG(TIJX(IT,JT,ISPN,IDOS))
                     AAHF(IBOTN+IT,IBOTS+JT) = AAHF(IBOTN+IT,IBOTS+JT)
     &                  + DIMAG(TIJX(IT,JT,IDOS,ISPN))
                  END DO
C
                  AAHF(IBOTO+IT,IBOTN+IT) = AAHF(IBOTO+IT,IBOTN+IT)
     &               + DIMAG(T0X(IT,IORB,IDOS))
                  AAHF(IBOTN+IT,IBOTO+IT) = AAHF(IBOTN+IT,IBOTO+IT)
     &               + DIMAG(T0X(IT,IDOS,IORB))
C
                  DO JT = 1,NT
                     AAHF(IBOTO+IT,IBOTN+JT) = AAHF(IBOTO+IT,IBOTN+JT)
     &                  + DIMAG(TIJX(IT,JT,IORB,IDOS))
                     AAHF(IBOTN+IT,IBOTO+JT) = AAHF(IBOTN+IT,IBOTO+JT)
     &                  + DIMAG(TIJX(IT,JT,IDOS,IORB))
                  END DO
C
               END IF
C
            END DO
C
C ======================================================================
C ======================================================================
C     calculate atom type-resolved  SPIN  and  ORBITAL  susceptibility
C ======================================================================
C ======================================================================
C                                                               S-O-LOOP
            DO ISOCOUPLING = 1,2
C
C--------------------------- set up and solve linear system of equations
C
               AA(:,:) = 0D0
               BB(:,:) = 0D0
C
               DO IT = 1,NT
C
                  BB(IBOTS+IT,1) = DIMAG(TZ(IT,ISPN,ISPN))
                  AA(IBOTS+IT,IBOTS+IT) = 1D0 - DIMAG(T0X(IT,ISPN,ISPN))
                  DO JT = 1,NT
                     AA(IBOTS+IT,IBOTS+JT) = AA(IBOTS+IT,IBOTS+JT)
     &                  - DIMAG(TIJX(IT,JT,ISPN,ISPN))
                  END DO
C
                  BB(IBOTO+IT,1) = DIMAG(TZ(IT,IORB,IORB))
                  AA(IBOTO+IT,IBOTO+IT) = 1D0 - DIMAG(T0X(IT,IORB,IORB))
C
                  DO JT = 1,NT
                     AA(IBOTO+IT,IBOTO+JT) = AA(NT+IT,NT+JT)
     &                  - DIMAG(TIJX(IT,JT,IORB,IORB))
                  END DO
C
               END DO
C
C---------------- exchange coupling of spin and orbital susceptibilities
C-------------- if ignored:  the coefficient matrix AA is block diagonal
C
               IF ( ISOCOUPLING.EQ.1 ) THEN
                  WRITE (6,99012) 'off'
               ELSE
                  WRITE (6,99012) 'on'
C
                  DO IT = 1,NT
                     BB(IBOTS+IT,1) = BB(IBOTS+IT,1)
     &                                + DIMAG(TZ(IT,ISPN,IORB))
                     AA(IBOTS+IT,IBOTO+IT) = AA(IBOTS+IT,IBOTO+IT)
     &                  + DIMAG(T0X(IT,ISPN,IORB))
                     DO JT = 1,NT
                        AA(IBOTS+IT,IBOTO+JT) = AA(IBOTS+IT,IBOTO+JT)
     &                     - DIMAG(TIJX(IT,JT,ISPN,IORB))
                     END DO
C
                     BB(IBOTO+IT,1) = BB(IBOTO+IT,1)
     &                                + DIMAG(TZ(IT,IORB,ISPN))
                     AA(IBOTO+IT,IBOTS+IT) = AA(IBOTO+IT,IBOTS+IT)
     &                  + DIMAG(T0X(IT,IORB,ISPN))
                     DO JT = 1,NT
                        AA(IBOTO+IT,IBOTS+JT) = AA(IBOTO+IT,IBOTS+JT)
     &                     - DIMAG(TIJX(IT,JT,IORB,ISPN))
                     END DO
C
                  END DO
C
               END IF
C
C------------------- add  terms connected with high-field susceptibility
C
               IF ( .NOT.NONMAG ) THEN
C
C------------------------------- the terms created by Fermi energy shift
                  IF ( ISOCOUPLING.EQ.1 ) THEN
C
                     DO I = 1,NT
                        BB(I,1) = BB(I,1) - BBHF_S(I) - BBHF_O(I)
     &                            - BBHF_N(I)
                        DO J = 1,NT
                           AA(I,J) = AA(I,J) - AAHF(I,J)
                        END DO
                     END DO
C
                     DO I = NT + 1,NDIM
                        BB(I,1) = BB(I,1) - BBHF_S(I) - BBHF_O(I)
     &                            - BBHF_N(I)
                        DO J = NT + 1,NDIM
                           AA(I,J) = AA(I,J) - AAHF(I,J)
                        END DO
                     END DO
C
                  ELSE
C
                     DO I = 1,NDIM
                        BB(I,1) = BB(I,1) - BBHF_S(I) - BBHF_O(I)
     &                            - BBHF_N(I)
                        DO J = 1,NDIM
                           AA(I,J) = AA(I,J) - AAHF(I,J)
                        END DO
                     END DO
C
                  END IF
C
               END IF
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
               CALL DGESV(NDIM,1,AA,NTMAX*3,IPIV,BB,NTMAX*3,INFO)
C
               IF ( INFO.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ZGESV: INFO != 0')
C
C-----------------------------------------------------------------------
C----------------------------------- determine enhanced susceptibilities
C
               DO IT = 1,NT
                  CHI_S(IT) = BB(IBOTS+IT,1)
                  CHI_O(IT) = BB(IBOTO+IT,1)
                  CHI_N(IT) = BB(IBOTN+IT,1)
C
                  CHI_TO(IT,ISPN) = CHI_S(IT)
                  CHI_TO(IT,IORB) = CHI_O(IT)
                  CHI_TO(IT,IDOS) = CHI_N(IT)
               END DO
C
C--------------------------------------- determine shift in Fermi energy
C
               IF ( .NOT.NONMAG ) THEN
C
                  SUMEF = SUM_WTZ_NS + SUM_WTZ_NO + SUM_WTZ_NN
C
                  DO JT = 1,NT
C
                     SUMEF = SUMEF + WTX_NS(JT)*CHI_S(JT)
C
                     SUMEF = SUMEF + WTX_NO(JT)*CHI_O(JT)
                     IF ( NT.EQ.1 ) THEN
                        SUMEF = SUMEF + WTX_NN(JT)
                     ELSE
                        SUMEF = SUMEF + WTX_NN(JT)*CHI_N(JT)
                     END IF
C             sumef = sumef + wtxnn(jt)*chin(jt)
                  END DO
C
               END IF
C
               IF ( .NOT.NONMAG ) WRITE (6,99013) ' WITHOUT '
C
C ======================================================================
C             determine remaining susceptibility terms
C ======================================================================
C
               DO IT = 1,NT
C
C-------------------------------------------------------------------- SS
C
                  CHI_SS0(IT) = DIMAG(TZ(IT,ISPN,ISPN))
C
                  CHI_SS(IT) = CHI_SS0(IT) + DIMAG(T0X(IT,ISPN,ISPN))
     &                         *CHI_S(IT)
C
                  DO JT = 1,NT
                     CHI_SS_TR(IT,JT) = DIMAG(TIJX(IT,JT,ISPN,ISPN))
     &                                  *CHI_S(JT)
                     CHI_SS(IT) = CHI_SS(IT) + CHI_SS_TR(IT,JT)
                  END DO
C
                  ICHI_S(IT) = 1.D0 - CHI_SS0(IT)/CHI_SS(IT)
C
C-------------------------------------------------------------------- OO
C
                  CHI_OO0(IT) = DIMAG(TZ(IT,IORB,IORB))
C
                  CHI_OO(IT) = CHI_OO0(IT) + DIMAG(T0X(IT,IORB,IORB))
     &                         *CHI_O(IT)
                  DO JT = 1,NT
                     CHI_OO_TR(IT,JT) = DIMAG(TIJX(IT,JT,IORB,IORB))
                     CHI_OO(IT) = CHI_OO(IT) + CHI_OO_TR(IT,JT)
                  END DO
C
                  ICHI_O(IT) = 1.D0 - CHI_OO0(IT)/CHI_OO(IT)
C
C-------------------------------------------------------------------- SO
C
                  CHI_SO0(IT) = DIMAG(TZ(IT,ISPN,IORB))
C
                  CHI_SO(IT) = CHI_SO0(IT) + DIMAG(T0X(IT,ISPN,IORB))
     &                         *CHI_O(IT)
                  DO JT = 1,NT
C
                     CHI_SO(IT) = CHI_SO(IT)
     &                            + DIMAG(TIJX(IT,JT,ISPN,IORB))
     &                            *CHI_O(JT)
                  END DO
C
C-------------------------------------------------------------------- OS
C
                  CHI_OS0(IT) = DIMAG(TZ(IT,IORB,ISPN))
C
                  CHI_OS(IT) = CHI_OS0(IT) + DIMAG(T0X(IT,IORB,ISPN))
     &                         *CHI_S(IT)
                  DO JT = 1,NT
                     CHI_OS(IT) = CHI_OS(IT)
     &                            + DIMAG(TIJX(IT,JT,IORB,ISPN))
     &                            *CHI_S(JT)
                  END DO
C
               END DO
C
C-----------------------------------------------------------------------
               IF ( IPRINT.GE.1 ) THEN
                  WRITE (6,99005) (IT,DIMAG(T1Z(IT,ISPN,ISPN)),DIMAG(T0Z
     &                            (IT,ISPN,ISPN)),IT=1,NT)
                  WRITE (6,99006) (IT,DIMAG(T0X(IT,ISPN,ISPN)),IT=1,NT)
C
                  DO IT = 1,NT
                     WRITE (6,'('' IT='',I2)') IT
                     WRITE (6,99007) (JT,'Z',DIMAG(TIJZ(IT,JT,ISPN,ISPN)
     &                               ),JT=1,NT)
                     WRITE (6,99007) (JT,'X',DIMAG(TIJX(IT,JT,ISPN,ISPN)
     &                               ),JT=1,NT)
                  END DO
C
                  WRITE (6,99008) (IT,DIMAG(T1Z(IT,IORB,IORB)),DIMAG(T0Z
     &                            (IT,IORB,IORB)),IT=1,NT)
                  WRITE (6,99009) (IT,DIMAG(T0X(IT,IORB,IORB)),IT=1,NT)
C
                  DO IT = 1,NT
                     WRITE (6,'('' IT='',I2)') IT
                     WRITE (6,99010) (JT,DIMAG(TIJZ(IT,JT,IORB,IORB)),JT
     &                               =1,NT)
                     WRITE (6,99011) (JT,DIMAG(TIJX(IT,JT,IORB,IORB)),JT
     &                               =1,NT)
                  END DO
               END IF
C
C ======================================================================
C                    calculate Knight shift
C ======================================================================
C
               IF ( HYPERFINE_INTERACTION ) THEN
C
                  DO IT = 1,NT
C
C------------------------------------------------------------------ HF-S
C
                     SFT_S0(IT) = DIMAG(TZ(IT,IHFI,ISPN))
                     SFT_S(IT) = SFT_S0(IT) + DIMAG(T0X(IT,IHFI,ISPN))
     &                           *CHI_S(IT)
                     DO JT = 1,NT
                        SFT_S_TR(IT,JT) = DIMAG(TIJX(IT,JT,IHFI,ISPN))
     &                     *CHI_S(JT)
                        SFT_S(IT) = SFT_S(IT) + SFT_S_TR(IT,JT)
                     END DO
C
C------------------------------------------------------------------ HF-O
C
                     SFT_O0(IT) = DIMAG(TZ(IT,IHFI,IORB))
                     SFT_O(IT) = SFT_O0(IT) + DIMAG(T0X(IT,IHFI,IORB))
     &                           *CHI_O(IT)
                     DO JT = 1,NT
                        SFT_O_TR(IT,JT) = DIMAG(TIJX(IT,JT,IHFI,IORB))
     &                     *CHI_O(JT)
                        SFT_O(IT) = SFT_O(IT) + SFT_O_TR(IT,JT)
                     END DO
C
                  END DO
C
               END IF
C
C ======================================================================
C            print out of SUSCEPTIBILITY and KNIGHT SHIFT
C                         WITHOUT HF-terms
C ======================================================================
C
C
               CALL LINRESP_SUSC_OUTPUT(BCP,BCPS,CHI_O,CHI_OO,CHI_OO0,
     &                                  CHI_OO_TR,CHI_OS,CHI_OS0,CHI_S,
     &                                  CHI_SO,CHI_SO0,CHI_SS,CHI_SS0,
     &                                  CHI_N,CHI_NN,CHI_NN0,CHI_SS_TR,
     &                                  DOBSEF,ICHI_O,ICHI_S,NONMAG,
     &                                  SFT_O,SFT_O0,SFT_O_TR,SFT_S,
     &                                  SFT_S0,SFT_S_TR,TCDIACRI,
     &                                  TCDIARI,TKDIACRI,TKDIARI,TOBS,
     &                                  NOBSEDNS,HYPERFINE_INTERACTION)
C
C
C ======================================================================
C                                      for FERRO add additional HF-terms
               IF ( .NOT.NONMAG ) THEN
C
                  WRITE (6,99013) 'INCLUDING'
C
C ======================================================================
C             determine remaining susceptibility terms
C ======================================================================
C
                  DO IT = 1,NT
C
C--------------------------------------------------------- high field SS
C
                     WIS = POBS(IT,ISPN)
                     CHI_SS0(IT) = CHI_SS0(IT) - WIS*SUM_WTZ_NS
                     CHI_SS(IT) = CHI_SS(IT) - WIS*SUM_WTZ_NS
C
                     DO JT = 1,NT
                        CHI_SS(IT) = CHI_SS(IT) - WIS*WTX_NS(JT)
     &                               *CHI_S(JT)
                     END DO
C
                     WIS = POBS(IT,ISPN)
                     CHI_SS(IT) = CHI_SS(IT) - WIS*SUM_WTZ_NS
                     DO JT = 1,NT
                        CHI_SS(IT) = CHI_SS(IT) - WIS*WTX_NS(JT)
     &                               *CHI_S(JT)
                     END DO
C
C--------------------------------------------------------- high field OO
C
                     WIO = POBS(IT,IORB)
                     CHI_OO0(IT) = CHI_OO0(IT) - WIO*SUM_WTZ_NO
                     CHI_OO(IT) = CHI_OO(IT) - WIO*SUM_WTZ_NO
                     DO JT = 1,NT
                        CHI_OO(IT) = CHI_OO(IT) - WIO*WTX_NO(JT)
     &                               *CHI_O(JT)
                     END DO
C
                     WIO = POBS(IT,IORB)
                     CHI_OO(IT) = CHI_OO(IT) - WIO*SUM_WTZ_NO
                     DO JT = 1,NT
                        CHI_OO(IT) = CHI_OO(IT) - WIO*WTX_NO(JT)
     &                               *CHI_O(JT)
                     END DO
C
C--------------------------------------------------------- high field SO
C
                     WIS = POBS(IT,ISPN)
C
                     CHI_SO0(IT) = CHI_SO0(IT) - WIS*SUM_WTZ_NO
                     CHI_SO(IT) = CHI_SO(IT) - WIS*SUM_WTZ_NO
C
                     DO JT = 1,NT
                        CHI_SO(IT) = CHI_SO(IT) - WIS*WTX_NO(JT)
     &                               *CHI_O(JT)
                     END DO
C
C--------------------------------------------------------- high field OS
C
                     WIO = POBS(IT,IORB)
C
                     CHI_OS0(IT) = CHI_OS0(IT) - WIO*SUM_WTZ_NS
                     CHI_OS(IT) = CHI_OS(IT) - WIO*SUM_WTZ_NS
C
                     DO JT = 1,NT
                        CHI_OS(IT) = CHI_OS(IT) - WIO*WTX_NS(JT)
     &                               *CHI_S(JT)
                     END DO
C
C--------------------------------------------------------- high field NN
C
                     WIN = POBS(IT,IDOS)
                     CHI_NN0(IT) = -WIN*SUM_WTZ_NN
                     CHI_NN(IT) = CHI_NN0(IT)
                     IF ( NT.EQ.1 ) THEN
                        CHI_NN(IT) = CHI_NN(IT)
     &                               + DIMAG(T0X(IT,IDOS,IDOS))
                        DO JT = 1,NT
                           CHI_NN(IT) = CHI_NN(IT) - WIN*WTX_NN(JT)
                           CHI_NN(IT) = CHI_NN(IT)
     &                                  + DIMAG(TIJX(IT,JT,IDOS,IDOS))
                        END DO
                     ELSE
                        CHI_NN(IT) = CHI_NN(IT)
     &                               + DIMAG(T0X(IT,IDOS,IDOS))
     &                               *CHI_N(IT)
                        DO JT = 1,NT
                           CHI_NN(IT) = CHI_NN(IT) - WIN*WTX_NN(JT)
     &                                  *CHI_N(JT)
                           CHI_NN(IT) = CHI_NN(IT)
     &                                  + DIMAG(TIJX(IT,JT,IDOS,IDOS))
     &                                  *CHI_N(JT)
                        END DO
                     END IF
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C----------------------- Calculation of the Hartree contribution
C
C               DO I = 1,JTOP
C                  VH(I) = 0.D0
C                  DO IP = 1,JTOP
C                     RINT(IP) = CHI_NR(IP,IT)*R2DRDI(IP,IM)
C                  END DO
C                  CALL RRADINT(IM,RINT,CHI_SUM)
C                  VH(I) = VH(I) + CHI_SUM/R(I,IM)
CC
C                  DO IP = 1,JTOP
C                     RINT(IP) = CHI_NR(IP,IT)*R(IP,IM)*DRDI(IP,IM)
C                  END DO
C                  CALL RRADINT(IM,RINT,CHI_SUM)
C                  CALL RRADINT(IM,RINT,CHI_SUM1)
C                  VH(I) = VH(I) + CHI_SUM - CHI_SUM1
C               END DO
C
C-----------------------------------------------------------------------
C
                  END DO
C
C-----------------------------------------------------------------------
                  IF ( IPRINT.GE.1 ) THEN
                     WRITE (6,99005) (IT,DIMAG(T1Z(IT,ISPN,ISPN)),DIMAG(
     &                               T0Z(IT,ISPN,ISPN)),IT=1,NT)
                     WRITE (6,99006) (IT,DIMAG(T0X(IT,ISPN,ISPN)),IT=1,
     &                               NT)
C
                     DO IT = 1,NT
                        WRITE (6,'('' IT='',I2)') IT
                        WRITE (6,99007) (JT,'Z',DIMAG(TIJZ(IT,JT,ISPN,
     &                                  ISPN)),JT=1,NT)
                        WRITE (6,99007) (JT,'X',DIMAG(TIJX(IT,JT,ISPN,
     &                                  ISPN)),JT=1,NT)
                     END DO
C
                     WRITE (6,99008) (IT,DIMAG(T1Z(IT,IORB,IORB)),DIMAG(
     &                               T0Z(IT,IORB,IORB)),IT=1,NT)
                     WRITE (6,99009) (IT,DIMAG(T0X(IT,IORB,IORB)),IT=1,
     &                               NT)
C
                     DO IT = 1,NT
                        WRITE (6,'('' IT='',I2)') IT
                        WRITE (6,99010) (JT,DIMAG(TIJZ(IT,JT,IORB,IORB))
     &                                  ,JT=1,NT)
                        WRITE (6,99011) (JT,DIMAG(TIJX(IT,JT,IORB,IORB))
     &                                  ,JT=1,NT)
                     END DO
                  END IF
C
C ======================================================================
C                    calculate Knight shift
C ======================================================================
C
                  IF ( HYPERFINE_INTERACTION ) THEN
C
                     DO IT = 1,NT
C
C------------------------------------------------------ high field  HF-S
C
C
                        WIHF = POBS(IT,IHFI)
                        SFT_S0(IT) = SFT_S0(IT) + WIHF*SUM_WTZ_NS
                        SFT_S(IT) = SFT_S(IT) + WIHF*SUM_WTZ_NS
                        DO JT = 1,NT
                           SFT_S(IT) = SFT_S(IT) + WIHF*WTX_NS(JT)
     &                                 *CHI_S(JT)
                        END DO
C
C------------------------------------------------------ high field  HF-O
C
                        WIHF = POBS(IT,IHFI)
                        SFT_O0(IT) = SFT_O0(IT) + WIHF*SUM_WTZ_NO
                        SFT_O(IT) = SFT_O(IT) + WIHF*SUM_WTZ_NO
                        DO JT = 1,NT
                           SFT_O(IT) = SFT_O(IT) + WIHF*WTX_NO(JT)
     &                                 *CHI_O(JT)
                        END DO
C
                        WIHF = POBS(IT,IHFI)
C
                        SFT_O0(IT) = SFT_O0(IT) + WIHF*SUM_WTZ_NO
                        SFT_O(IT) = SFT_O(IT) + WIHF*SUM_WTZ_NO
                        DO JT = 1,NT
                           SFT_O(IT) = SFT_O(IT) + WIHF*WTX_NO(JT)
     &                                 *CHI_O(JT)
                        END DO
C
                     END DO
C
                  END IF
C
C ======================================================================
C            print out of SUSCEPTIBILITY and KNIGHT SHIFT
C                         including  HF-terms
C ======================================================================
C
                  CALL LINRESP_SUSC_OUTPUT(BCP,BCPS,CHI_O,CHI_OO,
     &               CHI_OO0,CHI_OO_TR,CHI_OS,CHI_OS0,CHI_S,CHI_SO,
     &               CHI_SO0,CHI_SS,CHI_SS0,CHI_N,CHI_NN,CHI_NN0,
     &               CHI_SS_TR,DOBSEF,ICHI_O,ICHI_S,NONMAG,SFT_O,SFT_O0,
     &               SFT_O_TR,SFT_S,SFT_S0,SFT_S_TR,TCDIACRI,TCDIARI,
     &               TKDIACRI,TKDIARI,TOBS,NOBSEDNS,
     &               HYPERFINE_INTERACTION)
C
               END IF
C                                      for FERRO all terms sumed up now
C ======================================================================
C
            END DO
C                                                               S-O-LOOP
C ======================================================================
C
C
C
C ======================================================================
C                      iterate Delta B_xc
C ======================================================================
C
C ----------------------------------------------------------------------
C                        update the densities
C ----------------------------------------------------------------------
C
            W1MIX(1:NRMAX,1:NTMAX) = RHO2NS_GG(1:NRMAX,1,1:NTMAX,1,1)
            W2MIX(1:NRMAX,1:NTMAX) = RHO2NS_GG(1:NRMAX,1,1:NTMAX,2,2)
C
            W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &         = RHO2NS_GG(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1,1)
            W2NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
     &         = RHO2NS_GG(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,2,2)
C
            IF ( ITRSCF.EQ.1 ) THEN
C
               CALL SCFBROYPT1(IPRINT,0,IMIX,ISTBRY,ITDEPT,SCFMIX,W1MIX,
     &                         W2MIX,W1NSMIX,W2NSMIX,RMSAVV,RMSAVB,
     &                         JRNS1TMP,JRNSMINTMP)
C
            ELSE
C
               CALL SCFBROYPT1(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,SCFMIX,
     &                         W1MIX,W2MIX,W1NSMIX,W2NSMIX,RMSAVV,
     &                         RMSAVB,JRNS1TMP,JRNSMINTMP)
C
               RHO2NS_GG(1:NRMAX,1,1:NTMAX,1,1) = W1MIX(1:NRMAX,1:NTMAX)
C
               RHO2NS_GG(1:NRMAX,1,1:NTMAX,2,2) = W2MIX(1:NRMAX,1:NTMAX)
C
               RHO2NS_GG(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,1,1)
     &            = W1NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
C
               RHO2NS_GG(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX,2,2)
     &            = W2NSMIX(JRNSMINTMP:NRMAX,2:NLMFPMAX,1:NTMAX)
C
            END IF
C
            IF ( (RMSAVV.LE.SCFTOL) .AND. (RMSAVB.LE.SCFTOL) ) EXIT
C
C ----------------------------------------------------------------------
C  determine change to exchange correlation potential to  V  and  B
C ----------------------------------------------------------------------
C
            CALL LINRESP_VXC(RHO2NS,RGNT_VSF,LMRGNT_VSF,NRGNT_VSF_LM,
     &                       NRGNT_VSF,Y_LEBGRID,W_LEBGRID,N_LEBGRID,
     &                       NMAX_LEBGRID)
C
C ----------------------------------------------------------------------
C                 print info on susceptibilities
C ----------------------------------------------------------------------
C
            WRITE (6,99014) IDOS,ISPN,IORB
            DO IT = 1,NT
               WRITE (6,99015) IT,CHI_TO(IT,IDOS),CHI_TO(IT,ISPN),
     &                         CHI_TO(IT,IORB)
            END DO
C ======================================================================
C
         END DO
C=======================================================================
C                      CHI - SCF - cycle  END
C=======================================================================
         CALL CPU_TIME(TIME2)
C
         WRITE (6,99002) TIME2 - TIME1
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                        loop over q_pert - vectors
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      CALL STOP_REGULAR(ROUTINE,' ')
C
99001 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*   ***   *   *   ***    ***      ***   *****   ***   *****  *'
     &  ,/,10X,
     &  '*  *   *  *   *  *   *  *   *    *   *    *    *   *    *    *'
     &  ,/,10X,
     &  '*  *      *   *  *      *        *        *    *   *    *    *'
     &  ,/,10X,
     &  '*   ***   *   *   ***   *     **  ***     *    *****    *    *'
     &  ,/,10X,
     &  '*      *  *   *      *  *            *    *    *   *    *    *'
     &  ,/,10X,
     &  '*  *   *  *   *  *   *  *   *    *   *    *    *   *    *    *'
     &  ,/,10X,
     &  '*   ***    ***    ***    ***      ***     *    *   *    *    *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99003 FORMAT (/,10X,'CHI-SCF-tolerance    ',F16.6,/,10X,
     &        'CHI-mixing-parameter ',F16.6,/,10X,
     &        'max. number of iterations',I12)
99004 FORMAT (/,5X,'no iteration of induced normalized densities',/)
99005 FORMAT (/,8(:,' IT=',I2,' T1SSZ   =',1E13.5,' T0SSZ   =',1E13.5,/)
     &        )
99006 FORMAT (8(:,' IT=',I2,' T0SSX   =',1E13.5,/))
99007 FORMAT (8(:,' JT=',I2,' TIJSS',A,'  =',1E13.5,/))
99008 FORMAT (/,16(:,' IT=',I2,' T1LLZS  =',1E13.5,' T0LLZS  =',1E13.5,
     &        /))
99009 FORMAT (32(:,' IT=',I2,' T0LLXS  =',1E13.5,/))
99010 FORMAT (32(:,' JT=',I2,' TIJLLZS =',1E13.5,/))
99011 FORMAT (32(:,' JT=',I2,' TIJLLXS =',1E13.5,/))
99012 FORMAT (//,1X,79('*'),/,1X,79('*'),/,10X,
     &       'exchange coupling of spin and orbital susceptibilities:  '
     &       ,A,/,1X,79('*'),/,1X,79('*'),/)
99013 FORMAT (//,1X,79('*'),/,1X,16X,'results  ',A,
     &        '  additional high-field terms',/,1X,79('*'),/)
99014 FORMAT (/,10X,'susceptibilites used for next iteration ',/,10X,
     &        ' IT     IDOS =',I2,'     ISPN =',I2,'     IORB =',I2,/)
99015 FORMAT (10X,I3,3F13.5)
99016 FORMAT (I5,'  QVEC_PERT   ',5X,3F10.5)
C
      END
C*==linresp_susc_output.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_SUSC_OUTPUT(BCP,BCPS,CHI_O,CHI_OO,CHI_OO0,
     &                               CHI_OO_TR,CHI_OS,CHI_OS0,CHI_S,
     &                               CHI_SO,CHI_SO0,CHI_SS,CHI_SS0,
     &                               CHI_N,CHI_NN,CHI_NN0,CHI_SS_TR,
     &                               DOBSEF,ICHI_O,ICHI_S,NONMAG,SFT_O,
     &                               SFT_O0,SFT_O_TR,SFT_S,SFT_S0,
     &                               SFT_S_TR,TCDIACRI,TCDIARI,TKDIACRI,
     &                               TKDIARI,TOBS,NOBSEDNS,
     &                               HYPERFINE_INTERACTION)
C   ********************************************************************
C   *                                                                  *
C   *   print results for susceptibility and Knight shift              *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *                                                                  *
C   *  indexing of operators and observables                           *
C   *                                                                  *
C   *  number of particles   1        1 IDOS                           *
C   *  spin moment           b s_z    2 ISPN                           *
C   *  orbital moment        b l_z    3 IORB   ____ NPERT              *
C   *  hyperfine interaction H_hf,z   4 IHFI                           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TYPES,ONLY:NT,NTMAX,TXT_T,NAT,CONC
      USE MOD_CONSTANTS,ONLY:RY_ERG,MB_CGS,ALPHA_FS,CHI_AU2CGS
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_LINRESP,ONLY:ICDIA,IKDIA,IORB
      IMPLICIT NONE
C*--LINRESP_SUSC_OUTPUT1209
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      LOGICAL DIAMAGNETISM
      PARAMETER (DIAMAGNETISM=.FALSE.)
      REAL*8 CU,KU,FKVV,FKDIA
      PARAMETER (CU=CHI_AU2CGS*1D+6,KU=100D0*MB_CGS/RY_ERG,
     &           FKVV=100D0*ALPHA_FS*ALPHA_FS,FKDIA=FKVV/3D0)
C
C Dummy arguments
C
      LOGICAL HYPERFINE_INTERACTION,NONMAG
      INTEGER NOBSEDNS
      REAL*8 BCP(NTMAX),BCPS(NTMAX),CHI_N(NTMAX),CHI_NN(NTMAX),
     &       CHI_NN0(NTMAX),CHI_O(NTMAX),CHI_OO(NTMAX),CHI_OO0(NTMAX),
     &       CHI_OO_TR(NTMAX,NTMAX),CHI_OS(NTMAX),CHI_OS0(NTMAX),
     &       CHI_S(NTMAX),CHI_SO(NTMAX),CHI_SO0(NTMAX),CHI_SS(NTMAX),
     &       CHI_SS0(NTMAX),CHI_SS_TR(NTMAX,NTMAX),
     &       DOBSEF(NTMAX,NOBSEDNS),ICHI_O(NTMAX),ICHI_S(NTMAX),
     &       SFT_O(NTMAX),SFT_O0(NTMAX),SFT_O_TR(NTMAX,NTMAX),
     &       SFT_S(NTMAX),SFT_S0(NTMAX),SFT_S_TR(NTMAX,NTMAX),
     &       TCDIACRI(NTMAX),TCDIARI(NTMAX),TKDIACRI(NTMAX),
     &       TKDIARI(NTMAX),TOBS(NTMAX,NOBSEDNS)
C
C Local variables
C
      REAL*8 CHITOT,CHI_SUM,SFT0
      INTEGER IOBSE,IT,JT
C
C*** End of declarations rewritten by SPAG
C
      CHITOT = 0D0
C
      DO IT = 1,NT
C
         WRITE (6,99006) IT,TXT_T(IT)
C
C----------------------------------------------- polarisation at E_Fermi
C
C
         IF ( .NOT.NONMAG .OR. IPRINT.GT.1 ) WRITE (6,99007)
     &        (DOBSEF(IT,IOBSE),IOBSE=1,IORB)
C
         IF ( .NOT.DIAMAGNETISM ) TCDIARI(IT) = 0D0
         CHI_SUM = CHI_S(IT) + CHI_O(IT) + TCDIARI(IT)/3D0
C
         CHITOT = CHITOT + NAT(IT)*CONC(IT)*CHI_SUM
C
         WRITE (6,99008) CHI_SUM*CU
C
         WRITE (6,99010) 'CHI S      :',CHI_S(IT)*CU
C
         WRITE (6,99010) 'CHI SS     :',CHI_SS(IT)*CU,CHI_SS0(IT)*CU,
     &                   (CHI_SS(IT)-CHI_SS0(IT))*CU
C
         WRITE (6,99010) 'ICHI(S)    :',ICHI_S(IT)
C
         WRITE (6,99010) 'CHI SO     :',CHI_SO(IT)*CU,CHI_SO0(IT)*CU,
     &                   (CHI_SO(IT)-CHI_SO0(IT))*CU
C
         WRITE (6,99010) 'CHI O      :',CHI_O(IT)*CU
C
         WRITE (6,99010) 'CHI OO     :',CHI_OO(IT)*CU,CHI_OO0(IT)*CU,
     &                   (CHI_OO(IT)-CHI_OO0(IT))*CU
C
         WRITE (6,99010) 'ICHI(O)    :',ICHI_O(IT)
C
         WRITE (6,99010) 'CHI OS     :',CHI_OS(IT)*CU,CHI_OS0(IT)*CU,
     &                   (CHI_OS(IT)-CHI_OS0(IT))*CU
C
         IF ( DIAMAGNETISM ) THEN
            WRITE (6,99012) 'CHI DIA    :',TCDIARI(IT)*CU/3D0,'core:',
     &                      TCDIACRI(IT)*CU/3D0,'  cb:',
     &                      (TCDIARI(IT)-TCDIACRI(IT))*CU/3D0
            WRITE (6,99011) TOBS(IT,ICDIA)*CU/3D0
         END IF
C
         WRITE (6,99010) 'CHI N     :',CHI_N(IT)*CU
         WRITE (6,99010) 'CHI NN     :',CHI_NN(IT)*CU,CHI_NN0(IT)*CU,
     &                   (CHI_NN(IT)-CHI_NN0(IT))*CU
C
C
         IF ( HYPERFINE_INTERACTION ) THEN
C -------------------------------------------------- OUTPUT KNIGHT SHIFT
C
            SFT0 = SFT_S(IT) + SFT_O(IT) + BCP(IT)*CHI_S(IT)
C
            WRITE (6,99009) SFT0*KU + TKDIARI(IT)*FKDIA,SFT0*KU
C
            WRITE (6,99010) 'K HF-S  cb :',SFT_S(IT)*KU,SFT_S0(IT)*KU,
     &                      (SFT_S(IT)-SFT_S0(IT))*KU
C
            WRITE (6,99012) 'K HF-S  CP :',BCP(IT)*CHI_S(IT)*KU,'s-el',
     &                      BCPS(IT)*CHI_S(IT)*KU,'non-s',
     &                      (BCP(IT)-BCPS(IT))*CHI_S(IT)*KU
C
            WRITE (6,99010) 'K HF-O  cb :',SFT_O(IT)*KU,SFT_O0(IT)*KU,
     &                      (SFT_O(IT)-SFT_O0(IT))*KU
C
            WRITE (6,99012) 'K DIA      :',TKDIARI(IT)*FKDIA,' core',
     &                      TKDIACRI(IT)*FKDIA,'   cb',
     &                      (TKDIARI(IT)-TKDIACRI(IT))*FKDIA
            WRITE (6,99011) TOBS(IT,IKDIA)*FKDIA
C
            WRITE (6,*)
            WRITE (6,99012) 'BCP [MG]   :',BCP(IT)/1D6,' s-el',BCPS(IT)
     &                      /1D6,'non-s',(BCP(IT)-BCPS(IT))/1D6
C
         END IF
C
C-----------------------------------------------------------------------
C                                                 transferred quantities
C
         IF ( NT.GT.1 ) THEN
C
            WRITE (6,'(/,1X,79(''-''))')
C
            WRITE (6,99002) 'spin susceptibility CHI_SS'
            WRITE (6,99003) CHI_SS(IT)*CU
            DO JT = 1,NT
               WRITE (6,99004) JT,CHI_SS_TR(IT,JT)*CU
            END DO
C
            WRITE (6,99002) 'orbital susceptibility CHI_OO'
            WRITE (6,99003) CHI_OO(IT)*CU
            DO JT = 1,NT
               WRITE (6,99004) JT,CHI_OO_TR(IT,JT)*CU
               IF ( .NOT.NONMAG ) WRITE (6,99005) CHI_OO_TR(IT,JT)*CU
            END DO
C
            IF ( HYPERFINE_INTERACTION ) THEN
               WRITE (6,99002) 'Knight shift K HF-S'
               WRITE (6,99003) SFT_S(IT)*KU
               DO JT = 1,NT
                  WRITE (6,99004) JT,SFT_S_TR(IT,JT)*KU
               END DO
C
               WRITE (6,99002) 'Knight shift K HF-O'
               WRITE (6,99003) SFT_O(IT)*KU
               DO JT = 1,NT
                  WRITE (6,99004) JT,SFT_O_TR(IT,JT)*KU
               END DO
            END IF
C
         END IF
C
         WRITE (6,'(/,1X,79(''-''))')
C
C-----------------------------------------------------------------------
      END DO
C
      IF ( NT.EQ.1 ) THEN
         WRITE (6,'(/,1X,79(''=''))')
      ELSE
         WRITE (6,99001) CHITOT*CU
      END IF
C
99001 FORMAT (/,1X,79('='),//,10X,'total magnetic susceptibility',F15.6,
     &        ' x 10^(-6) cm^3/mol',//,1X,79('='),/)
99002 FORMAT (/,10X,'transferred ',A,//,23X,'sum         ',5(A,:,10X))
99003 FORMAT (10X,'total  ',6F11.4)
99004 FORMAT (10X,'JT =',I3,6F15.4)
99005 FORMAT (15X,A,6F11.4)
99006 FORMAT (/,1X,79('='),//,10X,'results for atom type  IT=',I2,2X,A)
99007 FORMAT (//,10X,'polarisations at the Fermi level',//,30X,
     &        'DOS         P_spin      P_orb',/,22X,F15.6,3F12.6/,
     &        (21X,A,F15.6,3F12.6))
99008 FORMAT (//,10X,'magnetic susceptibility',//,10X,'CHI tot    :',
     &        F15.6,' x 10^(-6) cm^3/mol')
99009 FORMAT (//,10X,'Knight shift',//,10X,'K   tot    :',F15.6,' % '/,
     &        10X,'K   (-dia) :',F15.6,' % ')
99010 FORMAT (/,10X,A12,F15.6,5X,:,' 0:',F12.6,5X,'XC:',F12.6)
99011 FORMAT (21X,A,43X,F12.6,/,(21X,A,43X,F12.6))
99012 FORMAT (/,10X,A12,F15.6,2X,:,A6,F12.6,3X,A5,F12.6)
      END
