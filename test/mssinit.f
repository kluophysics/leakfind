C*==mssinit.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MSSINIT(TSST,MSST,TSSQ,MSSQ)
C   ********************************************************************
C   *                                                                  *
C   *   SET UP THE EFFECTIVE m-MATRIX                                  *
C   *                                                                  *
C   *      CPASTART=1:     m(IQ)  = m(ATA) = 1 / t(ATA)                *
C   *                                                                  *
C   *                      t(ATA) =    SUM(i)  c(i) * t(i)             *
C   *                      =                          =                *
C   *                                                                  *
C   *      CPASTART=2:     unused at the moment                        *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   ALL matrices refer to the GLOBAL frame                         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *   THERMAL_VIBRA_FLUCT = .T.   dealing with thermal effects       *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_TYPES,ONLY:CONC,NTMAX
      USE MOD_CPA,ONLY:CPASTART
      USE MOD_SITES,ONLY:NQMAX,ITOQ,NOQ,IQBOT,IQTOP
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,NLM,WKM1,IPIVKM
      USE MOD_SYMMETRY,ONLY:ISYMGENQ,DROT,SYMUNITARY,IQREPQ,SYMACCEPTED,
     &    NSYM,IQPSYMQ,SYMMETRIZE_MSS
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_CALCMODE,ONLY:IREL,THERMAL_VIBRA_FLUCT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_THERMAL,ONLY:X_VFT,UMAT_VT,NVIBRA,NFLUCT,FMAT_FT
c modified by XJQ: scf of vibrations
      use mod_scfvb_cpa_sigma,only:lscfvb
c end-mod-xjq
      IMPLICIT NONE
C*--MSSINIT36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MSSINIT')
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TSSQ(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
c modified by XJQ: scf of vibrations
      integer ivt_extend
c end-mod-xjq
      INTEGER IA_ERR,IFLUCT,IFT,IO,IQ,IQREP,ISYM,IT,IVFT,IVIBRA,IVT,M,N,
     &        NSYMLOC
      LOGICAL KDONEQ(:)
      COMPLEX*16 TSS_FT(:,:),TSS_VFT(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KDONEQ,TSS_FT,TSS_VFT
C
      ALLOCATE (KDONEQ(NQMAX),TSS_FT(NKMMAX,NKMMAX))
      ALLOCATE (TSS_VFT(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KDONEQ')
C
      IF ( CPASTART.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'CPASTART <> 1')
C
      M = NKMMAX
C
      MSSQ(1:M,1:M,1:NQMAX) = C0
C
      DO IQ = IQBOT,IQTOP
         KDONEQ(IQ) = .FALSE.
      END DO
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ = IQBOT,IQTOP
C
         IF ( CHECK ) WRITE (6,99001) IQ,IREL,NKM,NLM
C
         IF ( IREL.NE.2 ) THEN
            N = NKMQ(IQ)
         ELSE
            N = NKM
         END IF
C
C=======================================================================
C                   representative site
C=======================================================================
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            IF ( NOQ(IQ).EQ.1 ) CONC(ITOQ(1,IQ)) = 1D0
C
C ----------------------------------------------------------------------
C        t(ata) = sum(IT) c(IT) * t(IT)     m(ata) = 1 / t(ata)
C ----------------------------------------------------------------------
C
            TSSQ(:,:,IQ) = C0
C
C================================================================= IO ==
            LOOP_IO:DO IO = 1,NOQ(IQ)
C
               IT = ITOQ(IO,IQ)
C
               IF ( THERMAL_VIBRA_FLUCT ) THEN
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
     &                              TSST(1,1,IT),'SAS+',TSS_FT)
C
C============================================================= IVIBRA ==
C                                             perform local displacement
                     IVT = (IT-1)*NVIBRA
                     LOOP_IVIBRA:DO IVIBRA = 1,NVIBRA
                        IVT = IVT + 1
C
                        IVFT = (IT-1)*NVIBRA*NFLUCT + (IVIBRA-1)
     &                         *NFLUCT + IFLUCT
C
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                        CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,
     &                                    'UAUT',TSS_VFT)
C
                        TSSQ(1:N,1:N,IQ) = TSSQ(1:N,1:N,IQ)
     &                     + X_VFT(IVFT)*TSS_VFT(1:N,1:N)
C
                     END DO LOOP_IVIBRA
C============================================================= IVIBRA ==
C
                  END DO LOOP_IFLUCT
C============================================================= IFLUCT ==
C
               ELSE
C
C-----------------------------------------------------------------------
C                  standard mode: NO thermal effects
C-----------------------------------------------------------------------
C
c modified by XJQ: scf of vibrations
                 if(lscfvb) then
                   if(mod(it,nvibra)==0) then
                     ivt_extend = (it/nvibra-1) * nvibra*nvibra +
     &                            nvibra
                   else
                     ivt_extend = (it/nvibra) * nvibra*nvibra +
     &                            mod(it,nvibra)
                   endif
                   CALL CMAT_U_TRANS(UMAT_VT(1,1,ivt_extend),
     &                               TSST(1,1,IT),'UAUT',TSS_VFT)
                   TSSQ(1:N,1:N,IQ) = TSSQ(1:N,1:N,IQ) + CONC(IT)
     &                                *TSS_VFT(1:N,1:N)
                 else
                   TSSQ(1:N,1:N,IQ) = TSSQ(1:N,1:N,IQ) + CONC(IT)
     &                                *TSST(1:N,1:N,IT)
                 endif
c end-mod-xjq
C
               END IF
C-----------------------------------------------------------------------
C
            END DO LOOP_IO
C================================================================= IO ==
C
            IF ( CHECK ) CALL CMATSTRUCT('<TSS>Q',TSSQ(1,1,IQ),NKM,
     &           NKMMAX,IREL,IREL,0,1D-8,6)
C
C--------------------------------------invert to get MSSQ: m[q] = 1/t[q]
C
            CALL CMATINV3(N,M,IPIVKM,TSSQ(1,1,IQ),WKM1,MSSQ(1,1,IQ))
C
C
C----------------------------------------------------------------------
C                         symmetrize MSSQ
C----------------------------------------------------------------------
C
            IF ( SYMMETRIZE_MSS ) THEN
C
               NSYMLOC = 0
C
               WKM1(1:NKM,1:NKM) = C0
C
               LOOP_ISYM:DO ISYM = 1,NSYM
C
                  IF ( .NOT.SYMACCEPTED(ISYM) ) CYCLE LOOP_ISYM
                  IF ( IQ.NE.IQPSYMQ(ISYM,IQ) ) CYCLE LOOP_ISYM
C
                  NSYMLOC = NSYMLOC + 1
C
                  CALL ROTATEGEN(MSSQ(1,1,IQ),WKM1,NKM,DROT(1,1,ISYM),
     &                           SYMUNITARY(ISYM),NKMMAX)
C
               END DO LOOP_ISYM
C
               MSSQ(1:NKM,1:NKM,IQ) = WKM1(1:NKM,1:NKM)/DBLE(NSYMLOC)
C
            END IF
C
            KDONEQ(IQ) = .TRUE.
C
C=======================================================================
C                NON - representative site
C=======================================================================
C
         ELSE
C
            IQREP = IQREPQ(IQ)
C
            IF ( .NOT.KDONEQ(IQREP) ) THEN
               WRITE (6,*) '###### IQ IQREPQ ',IQ,IQREP
               CALL STOP_MESSAGE(ROUTINE,'.NOT.KDONEQ(IQREP)')
            END IF
C
            IF ( FULLPOT ) THEN
               ISYM = ISYMGENQ(IQ)
               CALL ROTATEGEN(MSSQ(1,1,IQREP),MSSQ(1,1,IQ),N,
     &                        DROT(1,1,ISYM),SYMUNITARY(ISYM),NKMMAX)
            ELSE
               MSSQ(1:N,1:N,IQ) = MSSQ(1:N,1:N,IQREP)
            END IF
C
C--------------------------------------invert to get TSSQ: t[q] = 1/m[q]
C
            CALL CMATINV3(N,M,IPIVKM,MSSQ(1,1,IQ),WKM1,TSSQ(1,1,IQ))
C
C
         END IF
C=======================================================================
C
         IF ( CHECK ) CALL CMATSTRUCT('MSS Q',MSSQ(1,1,IQ),NKM,NKMMAX,
     &                                IREL,IREL,0,1D-8,6)
      END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (KDONEQ,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (//,1X,79('#'),/,10X,'IQ=',I3,'  IREL=',I3,'  NKM=',I3,
     &        '  NLM=',I3,/,1X,79('#'),/)
      END
