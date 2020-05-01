C*==gilme.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GILME(MSST,TSST,MSSQ,TAUQ,MAQAB,MBQAB,MCQAB,MDQAB,
     &                 MAQOAB,MBQOAB,MCQOAB,MDQOAB,MREG_FT)
C   ********************************************************************
C   *                                                                  *
C   *   - calculate the ALPHA matrix elements                          *
C   *   - transform the matrix elements for different E-arguments      *
C   *                                                                  *
C   *   the l-expansion is set via   NKM   for all atomic types to NL  *
C   *                                                                  *
C   *   all matrix elements are calculated in the LOCAL and then       *
C   *   then transfered to the GLOBAL   frame of reference             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_FILES,ONLY:IPRINT,CPOL_CART
      USE MOD_ANGMOM,ONLY:NKMMAX,NKM,WKM1,WKM2,IPIVKM
      USE MOD_SITES,ONLY:NQMAX,NOMAX,NOQ,ITOQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_TYPES,ONLY:NTMAX,NT,CONC
      USE MOD_CONSTANTS,ONLY:C1,C0
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_THERMAL,ONLY:NVFT,X_VFT,UMAT_VT,NVIBFLU,NVIBRA,NFLUCT,
     &    FMAT_FT,NFTMAX
      IMPLICIT NONE
C*--GILME27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GILME')
C
C Dummy arguments
C
      COMPLEX*16 MAQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MAQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MBQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MBQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MCQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MCQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MDQAB(NKMMAX,NKMMAX,3,NQMAX),
     &           MDQOAB(NKMMAX,NKMMAX,3,NQMAX,NOMAX),
     &           MREG_FT(NKMMAX,NKMMAX,3,NFTMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 CWGT,DMAMC(:,:),DMATTG(:,:),DTILTG(:,:),MREG_VFT(:,:,:)
     &           ,MSS_VFT(:,:),TSS_FT(:,:),TSS_VFT(:,:)
      INTEGER IA_ERR,IFLUCT,IFT,IO,IQ,IT,IVFT,IVIBRA,IVT,M,MUE,N
      REAL*8 T,X
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMAMC,DMATTG,DTILTG,TSS_FT,TSS_VFT,MSS_VFT,MREG_VFT
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (MSS_VFT(NKMMAX,NKMMAX),MREG_VFT(NKMMAX,NKMMAX,3))
      ALLOCATE (TSS_VFT(NKMMAX,NKMMAX),TSS_FT(NKMMAX,NKMMAX))
      ALLOCATE (DMAMC(NKMMAX,NKMMAX),DMATTG(NKMMAX,NKMMAX))
      ALLOCATE (DTILTG(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ERR')
C
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements                 *
C   *                                                                  *
C   ********************************************************************
C
      M = NKMMAX
      N = NKM
C
      IF ( THERMAL_VIBRA_FLUCT ) THEN
C
C ----------------------------------------------------------------------
C       calculate matrix MREG_VFT in global frame of references
C-----------------------------------------------------------------------
C
         DO IT = 1,NT
C
            IVFT = (IT-1)*NVIBFLU
C
            IVT = (IT-1)*NVIBRA
            DO IVIBRA = 1,NVIBRA
               IVT = IVT + 1
C
               IFT = (IT-1)*NFLUCT
               DO IFLUCT = 1,NFLUCT
                  IFT = IFT + 1
C
                  IVFT = IVFT + 1
C
C------------------------------------------------------------- rotations
C      generate torque matrix elements MREG_FT for tilted
C      magnetic moments represented within the non-rotated displaced
C      frames of references
C-----------------------------------------------------------------------
C
                  CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                           MSST(1,1,IT),'SAS+',WKM2)
C
                  CALL METORQUE(MREG_FT(1,1,1,IFT),WKM2)
C
               END DO
            END DO
C
         END DO
C
         IF ( NVFT.NE.IVFT ) CALL STOP_MESSAGE(ROUTINE,'NVFT <> IVFT')
C
      ELSE
C
         IF ( NVFT.NE.NT ) CALL STOP_MESSAGE(ROUTINE,
     &        'NVFT <> NT for THERMAL_VIBRA_FLUCT=.FALSE.')
C
         DO IT = 1,NT
            CALL METORQUE(MREG_FT(1,1,1,IT),MSST(1,1,IT))
         END DO
C
      END IF
C
C   ********************************************************************
C
C     calculate MBAR (Butler eq. 78)
C     ------------------------------
C
      MAQAB(:,:,:,:) = C0
      MBQAB(:,:,:,:) = C0
      MCQAB(:,:,:,:) = C0
      MDQAB(:,:,:,:) = C0
C
      MAQOAB(:,:,:,:,:) = C0
      MBQOAB(:,:,:,:,:) = C0
      MCQOAB(:,:,:,:,:) = C0
      MDQOAB(:,:,:,:,:) = C0
C
C================================================================= IQ ==
      LOOP_IQ:DO IQ = IQBOT_CHI,IQTOP_CHI
C
         IF ( IQREPQ(IQ).NE.IQ .AND. FULLPOT )
     &        CALL STOP_MESSAGE(ROUTINE,
     &        'IQREPQ(IQ) <> IQ .AND. FULLPOT')
C
         IF ( IPRINT.EQ.5 ) WRITE (6,99001) IQ
C
C================================================================= IO ==
         LOOP_IO:DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
C
C-----------------------------------------------------------------------
C         thermal lattice vibrations and/or spin fluctuations
C-----------------------------------------------------------------------
C
C============================================================= IFLUCT ==
C                                                 perform local rotation
            IFT = (IT-1)*NFLUCT
            LOOP_IFLUCT:DO IFLUCT = 1,NFLUCT
               IFT = IFT + 1
C
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
               CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                        TSST(1,1,IT),'SAS+',TSS_FT)
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
C------------------------------- MREG_VFT: M_vft = U_v * M_ft * U_v^(-1)
C
                  DO MUE = 1,3
                     CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),
     &                                 MREG_FT(1,1,MUE,IFT),'UAUT',
     &                                 MREG_VFT(1,1,MUE))
                  END DO
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                  CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,'UAUT',
     &                              TSS_VFT)
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                  CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WKM1,MSS_VFT)
C
C---------------- get projection matrices DMATTG, DTILTG in global frame
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMATTG,DTILTG,DMAMC,N,
     &                         MSSQ(1,1,IQ),MSS_VFT,M)
C
C-----------------------------------------------------------------------
C
                  CWGT = X_VFT(IVFT)/CONC(IT)
C
C================================================================ MUE ==
                  LOOP_MUE:DO MUE = 1,3
C
C----------------------------------------------------- MBAR type A and B
C
                     CALL ZGEMM('N','N',N,N,N,C1,DTILTG,M,
     &                          MREG_VFT(1,1,MUE),M,C0,WKM1,M)
C
                     CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,DMATTG,M,C1,
     &                          MAQOAB(1,1,MUE,IQ,IO),M)
C
                     CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,DTILTG,M,C1,
     &                          MBQOAB(1,1,MUE,IQ,IO),M)
C
C----------------------------------------------------- MBAR type C and D
C
                     CALL ZGEMM('C','N',N,N,N,C1,DMATTG,M,
     &                          MREG_VFT(1,1,MUE),M,C0,WKM1,M)
C
                     CALL ZGEMM('N','N',N,N,N,CWGT,WKM1,M,DMATTG,M,C1,
     &                          MCQOAB(1,1,MUE,IQ,IO),M)
C
                     CALL ZGEMM('N','C',N,N,N,CWGT,WKM1,M,DTILTG,M,C1,
     &                          MDQOAB(1,1,MUE,IQ,IO),M)
C
                  END DO LOOP_MUE
C================================================================ MUE ==
C
               END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
            END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
            X = CONC(IT)
C
            MAQAB(:,:,:,IQ) = MAQAB(:,:,:,IQ) + X*MAQOAB(:,:,:,IQ,IO)
            MBQAB(:,:,:,IQ) = MBQAB(:,:,:,IQ) + X*MBQOAB(:,:,:,IQ,IO)
            MCQAB(:,:,:,IQ) = MCQAB(:,:,:,IQ) + X*MCQOAB(:,:,:,IQ,IO)
            MDQAB(:,:,:,IQ) = MDQAB(:,:,:,IQ) + X*MDQOAB(:,:,:,IQ,IO)
C
         END DO LOOP_IO
C================================================================= IO ==
C
         IF ( IPRINT.EQ.5 ) THEN
            DO MUE = 1,3
               WRITE (6,99002) IQ,CPOL_CART(MUE)
               T = 1D-8
               CALL CMATSTRUCT('MAQAB ',MAQAB(1,1,MUE,IQ),N,M,3,3,0,T,6)
               CALL CMATSTRUCT('MBQAB ',MBQAB(1,1,MUE,IQ),N,M,3,3,0,T,6)
               CALL CMATSTRUCT('MCQAB ',MCQAB(1,1,MUE,IQ),N,M,3,3,0,T,6)
               CALL CMATSTRUCT('MDQAB ',MDQAB(1,1,MUE,IQ),N,M,3,3,0,T,6)
            END DO
         END IF
C
      END DO LOOP_IQ
C================================================================= IQ ==
C
99001 FORMAT (/,'JBAR-Matrices IQ=',I2,/,44('='))
99002 FORMAT (/,' IQ =',I2,'   POL = ',A)
      END
