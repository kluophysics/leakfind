C*==sigme_aux.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGME_AUX(MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,
     &                     MDQBA,MREG_TAB,MREG_TBA,S10_VFTAB,S2_VFTAB,
     &                     S3_VFTAB,S4_VFTAB,S10AQAB,S10BQAB,S10CQAB,
     &                     S10DQAB,S2AQAB,S2BQAB,S2CQAB,S2DQAB,S3AQAB,
     &                     S3BQAB,S3CQAB,S3DQAB,S4AQAB,S4BQAB,S4CQAB,
     &                     S4DQAB,MSSTA,TSSTA,MSSQA,TAUQA_TMP,MSSTB,
     &                     TSSTB,MSSQB,TAUQB_TMP,TSS_FTA,TSS_FTB,
     &                     TSS_VFTA,TSS_VFTB,MSS_VFTA,MSS_VFTB,UMAT_VTA,
     &                     UMAT_VTB)
C   ********************************************************************
C   *                                                                  *
C   *   determine the site dependent average for the                   *
C   *   the auxilary response terms  S10, S2, S3 and S4                *
C   *   and matrix elements MAQAB, MBQAB, MCQAB, and MDQAB             *
C   *                                                                  *
C   *==================================================================*
C   *                                                                  *
C   *                         IZA       IZB                            *
C   *     METYPE = A:        2  (+)    2  (+)                          *
C   *              B:        2  (+)    1  (-)                          *
C   *              C:        1  (-)    2  (+)                          *
C   *              D:        1  (-)    1  (-)                          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,WKM1,IPIVKM
      USE MOD_SITES,ONLY:NQMAX,IQBOT_CHI,IQTOP_CHI,NOQ,ITOQ
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_THERMAL,ONLY:NVT,NFT,FMAT_FT,MFMAT_FT,NVIBFLU,NVIBRA,
     &    NFLUCT,NVFTMAX,X_VFT
      USE MOD_SIG,ONLY:NSPINPROJ,NSPR,LIST_ISPR
      IMPLICIT NONE
C*--SIGME_AUX38
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGMEAUX')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MAQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MBQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MCQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MDQBA(NKMMAX,NKMMAX,3,NSPINPROJ,NQMAX),
     &           MREG_TAB(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2),
     &           MREG_TBA(NKMMAX,NKMMAX,3,NSPINPROJ,NTMAX,2,2),
     &           MSSQA(NKMMAX,NKMMAX,NQMAX),MSSQB(NKMMAX,NKMMAX,NQMAX),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           MSS_VFTA(NKMMAX,NKMMAX),MSS_VFTB(NKMMAX,NKMMAX),
     &           S10AQAB(3,3,NSPINPROJ,NQMAX),
     &           S10BQAB(3,3,NSPINPROJ,NQMAX),
     &           S10CQAB(3,3,NSPINPROJ,NQMAX),
     &           S10DQAB(3,3,NSPINPROJ,NQMAX),
     &           S10_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S2AQAB(3,3,NSPINPROJ,NQMAX),S2BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S2CQAB(3,3,NSPINPROJ,NQMAX),
     &           S2DQAB(3,3,NSPINPROJ,NQMAX),
     &           S2_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S3AQAB(3,3,NSPINPROJ,NQMAX),S3BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S3CQAB(3,3,NSPINPROJ,NQMAX),
     &           S3DQAB(3,3,NSPINPROJ,NQMAX),
     &           S3_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S4AQAB(3,3,NSPINPROJ,NQMAX),S4BQAB(3,3,NSPINPROJ,NQMAX)
     &           ,S4CQAB(3,3,NSPINPROJ,NQMAX),
     &           S4DQAB(3,3,NSPINPROJ,NQMAX),
     &           S4_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           TAUQA_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TAUQB_TMP(NKMMAX,NKMMAX,NQMAX),
     &           TSSTA(NKMMAX,NKMMAX,NTMAX),TSSTB(NKMMAX,NKMMAX,NTMAX),
     &           TSS_FTA(NKMMAX,NKMMAX),TSS_FTB(NKMMAX,NKMMAX),
     &           TSS_VFTA(NKMMAX,NKMMAX),TSS_VFTB(NKMMAX,NKMMAX),
     &           UMAT_VTA(NKMMAX,NKMMAX,NVT),UMAT_VTB(NKMMAX,NKMMAX,NVT)
C
C Local variables
C
      COMPLEX*16 CWGT,DMAT_VFTA(:,:),DMAT_VFTAZ(:,:,:),DMAT_VFTB(:,:),
     &           DMAT_VFTBZ(:,:,:),DTIL_VFTA(:,:),DTIL_VFTAZ(:,:,:),
     &           DTIL_VFTB(:,:),DTIL_VFTBZ(:,:,:),MREG_FTAB(:,:,:,:,:,:)
     &           ,MREG_FTBA(:,:,:,:,:,:),MREG_VFTAB(:,:,:,:,:,:),
     &           MREG_VFTBA(:,:,:,:,:,:),UMAT_VTAZ(:,:,:),
     &           UMAT_VTBZ(:,:,:)
      CHARACTER*50 FMT_CHECK
      INTEGER I,IFLUCT,IFT,IO,IQ,ISP,IT,IVFT,IVIBRA,IVT,IXX,IZA,IZB,J,M,
     &        MUE,N
      REAL*8 T
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMAT_VFTA,DMAT_VFTB,DTIL_VFTA,DTIL_VFTB
      ALLOCATABLE DMAT_VFTAZ,DMAT_VFTBZ,DTIL_VFTAZ,DTIL_VFTBZ
      ALLOCATABLE MREG_FTAB,MREG_FTBA
      ALLOCATABLE MREG_VFTAB,MREG_VFTBA,UMAT_VTAZ,UMAT_VTBZ
C
      CALL TRACK_INFO(ROUTINE)
C
      FMT_CHECK = '(///,''CHECK '',50(''>''),4X,A,//,(2F26.12))'
C
      M = NKMMAX
      N = NKM
C
      ALLOCATE (UMAT_VTAZ(NKMMAX,NKMMAX,2))
      ALLOCATE (UMAT_VTBZ(NKMMAX,NKMMAX,2))
      ALLOCATE (DMAT_VFTA(M,M),DMAT_VFTB(M,M))
      ALLOCATE (DTIL_VFTA(M,M),DTIL_VFTB(M,M))
      ALLOCATE (DMAT_VFTAZ(M,M,2),DMAT_VFTBZ(M,M,2))
      ALLOCATE (DTIL_VFTAZ(M,M,2),DTIL_VFTBZ(M,M,2))
      ALLOCATE (MREG_FTAB(NKMMAX,NKMMAX,3,NSPINPROJ,2,2))
      ALLOCATE (MREG_FTBA(NKMMAX,NKMMAX,3,NSPINPROJ,2,2))
      ALLOCATE (MREG_VFTAB(NKMMAX,NKMMAX,3,NSPINPROJ,2,2))
      ALLOCATE (MREG_VFTBA(NKMMAX,NKMMAX,3,NSPINPROJ,2,2))
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
C       determine the site dependent average for the
C       the auxilary response terms  S10, S2, S3 and S4
C       and matrix elements MAQAB, MBQAB, MCQAB, and MDQAB
C=======================================================================
C
C                         IZA       IZB
C     METYPE = A:        2  (+)    2  (+)
C              B:        2  (+)    1  (-)
C              C:        1  (-)    2  (+)
C              D:        1  (-)    1  (-)
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      LOOP_IQ:DO IQ = IQBOT_CHI,IQTOP_CHI
C
         IF ( IQREPQ(IQ).NE.IQ .AND. FULLPOT )
     &        CALL STOP_MESSAGE(ROUTINE,
     &        'IQREPQ(IQ) <> IQ .AND. FULLPOT')
C
         IF ( IPRINT.EQ.5 ) WRITE (6,99001) IQ
C
         LOOP_IO:DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
            IFT = (IT-1)*NFLUCT
            LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
               IFT = IFT + 1
C
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
               IF ( NFLUCT.EQ.1 ) THEN
C
                  TSS_FTA(:,:) = TSSTA(:,:,IT)
                  TSS_FTB(:,:) = TSSTB(:,:,IT)
C
                  DO ISP = 1,NSPR
                     IXX = LIST_ISPR(ISP)
                     MREG_FTAB(:,:,:,IXX,:,:)
     &                  = MREG_TAB(:,:,:,IXX,IT,:,:)
                     MREG_FTBA(:,:,:,IXX,:,:)
     &                  = MREG_TBA(:,:,:,IXX,IT,:,:)
                  END DO
C
               ELSE
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           TSSTA(1,1,IT),'SAS+',TSS_FTA)
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           TSSTB(1,1,IT),'SAS+',TSS_FTB)
C
                  DO ISP = 1,NSPR
                     IXX = LIST_ISPR(ISP)
C
                     DO IZA = 1,2
                        DO IZB = 1,2
                           CALL ME_ROTATE_REG(IFT,NFT,FMAT_FT,MFMAT_FT,
     &                        MREG_TAB(1,1,1,IXX,IT,IZA,IZB),
     &                        MREG_TBA(1,1,1,IXX,IT,IZB,IZA),
     &                        MREG_FTAB(1,1,1,IXX,IZA,IZB),
     &                        MREG_FTBA(1,1,1,IXX,IZB,IZA))
                        END DO
                     END DO
C
                  END DO
C
               END IF
C
C============================================================= IVIBRA ==
C                                             perform local displacement
               IVT = (IT-1)*NVIBRA
               LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                  IVT = IVT + 1
C
                  IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                   *NFLUCT + IFLUCT
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                  IF ( NVIBRA.EQ.1 ) THEN
C
                     TSS_VFTA(:,:) = TSS_FTA(:,:)
                     TSS_VFTB(:,:) = TSS_FTB(:,:)
C
                     DO ISP = 1,NSPR
                        IXX = LIST_ISPR(ISP)
                        MREG_VFTAB(:,:,:,IXX,:,:)
     &                     = MREG_FTAB(:,:,:,IXX,:,:)
                        MREG_VFTBA(:,:,:,IXX,:,:)
     &                     = MREG_FTBA(:,:,:,IXX,:,:)
                     END DO
C
                  ELSE
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C---------------------------- NB: umat and t at energies above real axis
C
c modified by XJQ: scfvb
                     call read_scftv(ivt,TSS_FTA)
                     TSS_FTB(:,:)=TSS_FTA(:,:)
c end-mod-xjq
                     CALL CMAT_U_TRANS(UMAT_VTA(1,1,IVT),TSS_FTA,'UAUT',
     &                                 TSS_VFTA)
                     CALL CMAT_U_TRANS(UMAT_VTB(1,1,IVT),TSS_FTB,'UAUT',
     &                                 TSS_VFTB)
C
C------------------------------- MREG_VFT: M_vft = U_v * M_ft * U_v^(-1)
C
                     CALL GET_MATZ_UMAT(UMAT_VTA(1,1,IVT),UMAT_VTAZ)
                     CALL GET_MATZ_UMAT(UMAT_VTB(1,1,IVT),UMAT_VTBZ)
C
                     DO ISP = 1,NSPR
                        IXX = LIST_ISPR(ISP)
C
                        DO IZA = 1,2
                           DO IZB = 1,2
c modified by XJQ: scfvb
                              call read_mreg_scfvb(iza,izb,ivt,
     &                          mreg_ftab(:,:,:,ixx,iza,izb),
     &                          mreg_ftba(:,:,:,ixx,izb,iza))
c end-mod-xjq
                              CALL ME_UTRANS_REG(UMAT_VTAZ(1,1,IZA),
     &                           UMAT_VTBZ(1,1,IZB),
     &                           MREG_FTAB(1,1,1,IXX,IZA,IZB),
     &                           MREG_FTBA(1,1,1,IXX,IZB,IZA),
     &                           MREG_VFTAB(1,1,1,IXX,IZA,IZB),
     &                           MREG_VFTBA(1,1,1,IXX,IZB,IZA))
                           END DO
                        END DO
C
                     END DO
C
                  END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C------ NB: umat and t at energies above real axis --> m above real axis
C
                  IF ( NVIBFLU.EQ.1 ) THEN
                     MSS_VFTA(:,:) = MSSTA(:,:,IT)
                     MSS_VFTB(:,:) = MSSTB(:,:,IT)
                  ELSE
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFTA,WKM1,MSS_VFTA)
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFTB,WKM1,MSS_VFTB)
                  END IF
C
C---------------------------- get projection matrices DMAT_vft, DTIL_vft
C------------------------------------------------------- above real axis
C
                  CALL GETDMAT(TAUQA_TMP(1,1,IQ),DMAT_VFTA,DTIL_VFTA,
     &                         WKM1,NKM,MSSQA(1,1,IQ),MSS_VFTA,NKMMAX)
C
                  CALL GETDMAT(TAUQB_TMP(1,1,IQ),DMAT_VFTB,DTIL_VFTB,
     &                         WKM1,NKM,MSSQB(1,1,IQ),MSS_VFTB,NKMMAX)
C
C-----------------------------------------------------------------------
C        generate projection matrices for energies below real axis
C-----------------------------------------------------------------------
                  CALL GET_MATZ(DMAT_VFTA,DMAT_VFTAZ)
                  CALL GET_MATZ(DMAT_VFTB,DMAT_VFTBZ)
                  CALL GET_MATZ(DTIL_VFTA,DTIL_VFTAZ)
                  CALL GET_MATZ(DTIL_VFTB,DTIL_VFTBZ)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  IF ( CHECK ) THEN
                     IF ( IFLUCT.EQ.1 .AND. IVIBRA.EQ.1 ) THEN
                        WRITE (6,FMT=FMT_CHECK) 'FMAT        ',
     &                         FMAT_FT(:,:,IFT)
                        WRITE (6,FMT=FMT_CHECK) 'MFMAT       ',
     &                         MFMAT_FT(:,:,IFT)
C
                        WRITE (6,FMT=FMT_CHECK) 'MREG_TAB    ',
     &                         MREG_TAB(:,:,:,1,IT,:,:)
                        DO I = 1,NKM
                           DO J = 1,NKM
                              WRITE (6,'(A,2I3)') 'MREG_FTAB   ',I,J
                              WRITE (6,FMT=FMT_CHECK) 'MREG_FTAB   ',
     &                               MREG_FTAB(I,J,:,1,:,:)
                           END DO
                        END DO
                        WRITE (6,FMT=FMT_CHECK) 'MREG_VFTAB  ',
     &                         MREG_VFTAB(:,:,:,1,:,:)
C
                        WRITE (6,FMT=FMT_CHECK) 'DMAT_VFTA   ',DMAT_VFTA
                        WRITE (6,FMT=FMT_CHECK) 'DMAT_VFTB   ',DMAT_VFTB
                        WRITE (6,FMT=FMT_CHECK) 'DTIL_VFTA   ',DTIL_VFTA
                        WRITE (6,FMT=FMT_CHECK) 'DTIL_VFTB   ',DTIL_VFTB
                     END IF
                  END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C-----------------------------------------------------------------------
C
                  CWGT = X_VFT(IVFT)
C
                  LOOP_IXX_B:DO ISP = 1,NSPR
                     IXX = LIST_ISPR(ISP)
C
                     LOOP_MUE:DO MUE = 1,3
C
C----------------------------------------------------- MBAR type A and B
C--------------------------------------------------------------energy AB
C--------------------------------------------------------------------- A
C
                        CALL ZGEMM('N','N',N,N,N,C1,DTIL_VFTAZ(1,1,2),M,
     &                             MREG_VFTAB(1,1,MUE,IXX,2,2),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DMAT_VFTBZ(1,1,2),M,C1,
     &                             MAQAB(1,1,MUE,IXX,IQ),M)
C
C--------------------------------------------------------------------- B
                        CALL ZGEMM('N','N',N,N,N,C1,DTIL_VFTAZ(1,1,2),M,
     &                             MREG_VFTAB(1,1,MUE,IXX,2,1),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DTIL_VFTBZ(1,1,1),M,C1,
     &                             MBQAB(1,1,MUE,IXX,IQ),M)
C
C--------------------------------------------------------------energy BA
C--------------------------------------------------------------------- A
C
                        CALL ZGEMM('N','N',N,N,N,C1,DTIL_VFTBZ(1,1,2),M,
     &                             MREG_VFTBA(1,1,MUE,IXX,2,2),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DMAT_VFTAZ(1,1,2),M,C1,
     &                             MAQBA(1,1,MUE,IXX,IQ),M)
C
C--------------------------------------------------------------------- B
                        CALL ZGEMM('N','N',N,N,N,C1,DTIL_VFTBZ(1,1,2),M,
     &                             MREG_VFTBA(1,1,MUE,IXX,2,1),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DTIL_VFTAZ(1,1,1),M,C1,
     &                             MBQBA(1,1,MUE,IXX,IQ),M)
C
C----------------------------------------------------- MBAR type C and D
C--------------------------------------------------------------energy AB
C--------------------------------------------------------------------- C
C
                        CALL ZGEMM('N','N',N,N,N,C1,DMAT_VFTAZ(1,1,1),M,
     &                             MREG_VFTAB(1,1,MUE,IXX,1,2),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DMAT_VFTBZ(1,1,2),M,C1,
     &                             MCQAB(1,1,MUE,IXX,IQ),M)
C--------------------------------------------------------------------- D
C
                        CALL ZGEMM('N','N',N,N,N,C1,DMAT_VFTAZ(1,1,1),M,
     &                             MREG_VFTAB(1,1,MUE,IXX,1,1),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DTIL_VFTBZ(1,1,1),M,C1,
     &                             MDQAB(1,1,MUE,IXX,IQ),M)
C
C--------------------------------------------------------------energy BA
C--------------------------------------------------------------------- C
C
                        CALL ZGEMM('N','N',N,N,N,C1,DMAT_VFTBZ(1,1,1),M,
     &                             MREG_VFTBA(1,1,MUE,IXX,1,2),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DMAT_VFTAZ(1,1,2),M,C1,
     &                             MCQBA(1,1,MUE,IXX,IQ),M)
C
C--------------------------------------------------------------------- D
                        CALL ZGEMM('N','N',N,N,N,C1,DMAT_VFTBZ(1,1,1),M,
     &                             MREG_VFTBA(1,1,MUE,IXX,1,1),M,C0,
     &                             WKM1,M)
C
                        CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,
     &                             DTIL_VFTAZ(1,1,1),M,C1,
     &                             MDQBA(1,1,MUE,IXX,IQ),M)
C
C----------------------------------------- average terms S10, S2, S3, S4
C
                        CALL SIGME_AVG_SI(MUE,IXX,IVFT,IQ,CWGT,
     &                     S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,
     &                     NSPINPROJ,NVFTMAX,'A',S10AQAB,S2AQAB,S3AQAB,
     &                     S4AQAB)
C
                        CALL SIGME_AVG_SI(MUE,IXX,IVFT,IQ,CWGT,
     &                     S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,
     &                     NSPINPROJ,NVFTMAX,'B',S10BQAB,S2BQAB,S3BQAB,
     &                     S4BQAB)
C
                        CALL SIGME_AVG_SI(MUE,IXX,IVFT,IQ,CWGT,
     &                     S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,
     &                     NSPINPROJ,NVFTMAX,'C',S10CQAB,S2CQAB,S3CQAB,
     &                     S4CQAB)
C
                        CALL SIGME_AVG_SI(MUE,IXX,IVFT,IQ,CWGT,
     &                     S10_VFTAB,S2_VFTAB,S3_VFTAB,S4_VFTAB,
     &                     NSPINPROJ,NVFTMAX,'D',S10DQAB,S2DQAB,S3DQAB,
     &                     S4DQAB)
C
C----------------------------------------- average terms S10, S2, S3, S4
C
                     END DO LOOP_MUE
C
                  END DO LOOP_IXX_B
C
               END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
            END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
         END DO LOOP_IO
C================================================================= IO ==
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
         IF ( IPRINT.EQ.5 ) THEN
            DO MUE = 1,3
               DO ISP = 1,NSPR
                  IXX = LIST_ISPR(ISP)
                  WRITE (6,99002) MUE,IQ
                  T = 1D-8
                  CALL CMATSTRUCT('MAQAB ',MAQAB(1,1,MUE,IXX,IQ),N,M,3,
     &                            3,0,T,6)
                  CALL CMATSTRUCT('MBQAB ',MBQAB(1,1,MUE,IXX,IQ),N,M,3,
     &                            3,0,T,6)
                  CALL CMATSTRUCT('MCQAB ',MCQAB(1,1,MUE,IXX,IQ),N,M,3,
     &                            3,0,T,6)
                  CALL CMATSTRUCT('MDQAB ',MDQAB(1,1,MUE,IXX,IQ),N,M,3,
     &                            3,0,T,6)
               END DO
            END DO
         END IF
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
            WRITE (6,FMT=FMT_CHECK) 'MAQAB     ',MAQAB
            WRITE (6,FMT=FMT_CHECK) 'MBQAB     ',MBQAB
            WRITE (6,FMT=FMT_CHECK) 'MCQAB     ',MCQAB
            WRITE (6,FMT=FMT_CHECK) 'MDQAB     ',MDQAB
            WRITE (6,FMT=FMT_CHECK) 'MAQBA     ',MAQBA
            WRITE (6,FMT=FMT_CHECK) 'MBQBA     ',MBQBA
            WRITE (6,FMT=FMT_CHECK) 'MCQBA     ',MCQBA
            WRITE (6,FMT=FMT_CHECK) 'MDQBA     ',MDQBA
            DO IT = ITBOT,ITTOP
               IFT = (IT-1)*NFLUCT
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
                  IVT = (IT-1)*NVIBRA
                  DO IVIBRA = 1,NVIBRA
                     IVT = IVT + 1
                     IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                      *NFLUCT + IFLUCT
                     IF ( IFLUCT.EQ.1 .AND. IVIBRA.EQ.1 ) THEN
                        WRITE (6,'(/,''VFT '',4I3)') IVIBRA,IFLUCT,IT,
     &                         IVFT
                        WRITE (6,FMT=FMT_CHECK) 'S10_VFTAB     ',
     &                         S10_VFTAB(:,:,1,IVFT,:,:)
                        WRITE (6,FMT=FMT_CHECK) 'S2_VFTAB      ',
     &                         S2_VFTAB(:,:,1,IVFT,:,:)
                        WRITE (6,FMT=FMT_CHECK) 'S3_VFTAB      ',
     &                         S3_VFTAB(:,:,1,IVFT,:,:)
                        WRITE (6,FMT=FMT_CHECK) 'S4_VFTAB      ',
     &                         S4_VFTAB(:,:,1,IVFT,:,:)
                     END IF
                  END DO
               END DO
            END DO
            WRITE (6,FMT=FMT_CHECK) 'S10AQAB       ',S10AQAB
            WRITE (6,FMT=FMT_CHECK) 'S10BQAB       ',S10BQAB
            WRITE (6,FMT=FMT_CHECK) 'S10CQAB       ',S10CQAB
            WRITE (6,FMT=FMT_CHECK) 'S10DQAB       ',S10DQAB
            WRITE (6,FMT=FMT_CHECK) 'S2AQAB        ',S2AQAB
            WRITE (6,FMT=FMT_CHECK) 'S2BQAB        ',S2BQAB
            WRITE (6,FMT=FMT_CHECK) 'S2CQAB        ',S2CQAB
            WRITE (6,FMT=FMT_CHECK) 'S2DQAB        ',S2DQAB
            WRITE (6,FMT=FMT_CHECK) 'S3AQAB        ',S3AQAB
            WRITE (6,FMT=FMT_CHECK) 'S3BQAB        ',S3BQAB
            WRITE (6,FMT=FMT_CHECK) 'S3CQAB        ',S3CQAB
            WRITE (6,FMT=FMT_CHECK) 'S3DQAB        ',S3DQAB
            WRITE (6,FMT=FMT_CHECK) 'S4AQAB        ',S4AQAB
            WRITE (6,FMT=FMT_CHECK) 'S4BQAB        ',S4BQAB
            WRITE (6,FMT=FMT_CHECK) 'S4CQAB        ',S4CQAB
            WRITE (6,FMT=FMT_CHECK) 'S4DQAB        ',S4DQAB
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      END DO LOOP_IQ
C
      IF ( CHECK ) CALL STOP_REGULAR(ROUTINE,'CHECK completed')
C
99001 FORMAT (/,'JBAR-Matrices IQ=',I2,/,44('='))
99002 FORMAT (/,'direction (xyz)=(123)',i2,5x,' IQ=',i2)
      END
C*==sigme_avg_si.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGME_AVG_SI(MUE,IXX,IVFT,IQ,CWGT,S10_VFTAB,S2_VFTAB,
     &                        S3_VFTAB,S4_VFTAB,NSPINPROJ,NVFTMAX,
     &                        METYPE,S10Q,S2Q,S3Q,S4Q)
C   ********************************************************************
C   *                                                                  *
C   * for given site IQ average ALL rerturbation terms S_i             *
C   * over types IVFT with weight CWGT                                 *
C   *                                                                  *
C   *                         IZA       IZB                            *
C   *     METYPE = A:        2  (+)    2  (+)                          *
C   *              B:        2  (+)    1  (-)                          *
C   *              C:        1  (-)    2  (+)                          *
C   *              D:        1  (-)    1  (-)                          *
C   *                                                                  *
C   ********************************************************************
      USE MOD_SITES,ONLY:NQMAX
      IMPLICIT NONE
C*--SIGME_AVG_SI557
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 CWGT
      INTEGER IQ,IVFT,IXX,MUE,NSPINPROJ,NVFTMAX
      CHARACTER*1 METYPE
      COMPLEX*16 S10Q(3,3,NSPINPROJ,NQMAX),
     &           S10_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S2Q(3,3,NSPINPROJ,NQMAX),
     &           S2_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S3Q(3,3,NSPINPROJ,NQMAX),
     &           S3_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2),
     &           S4Q(3,3,NSPINPROJ,NQMAX),
     &           S4_VFTAB(3,3,NSPINPROJ,NVFTMAX,2,2)
C
C Local variables
C
      INTEGER IZBRA,IZKET,NUE
C
C*** End of declarations rewritten by SPAG
C
      IF ( METYPE.EQ.'A' .OR. METYPE.EQ.'B' ) THEN
         IZBRA = 2
      ELSE
         IZBRA = 1
      END IF
C
      IF ( METYPE.EQ.'A' .OR. METYPE.EQ.'C' ) THEN
         IZKET = 2
      ELSE
         IZKET = 1
      END IF
C
      DO NUE = 1,3
C
         S10Q(MUE,NUE,IXX,IQ) = S10Q(MUE,NUE,IXX,IQ)
     &                          + CWGT*S10_VFTAB(MUE,NUE,IXX,IVFT,IZBRA,
     &                          IZKET)
C
         S2Q(MUE,NUE,IXX,IQ) = S2Q(MUE,NUE,IXX,IQ)
     &                         + CWGT*S2_VFTAB(MUE,NUE,IXX,IVFT,IZBRA,
     &                         IZKET)
C
         S3Q(MUE,NUE,IXX,IQ) = S3Q(MUE,NUE,IXX,IQ)
     &                         + CWGT*S3_VFTAB(MUE,NUE,IXX,IVFT,IZBRA,
     &                         IZKET)
C
         S4Q(MUE,NUE,IXX,IQ) = S4Q(MUE,NUE,IXX,IQ)
     &                         + CWGT*S4_VFTAB(MUE,NUE,IXX,IVFT,IZBRA,
     &                         IZKET)
C
      END DO
C
      END
