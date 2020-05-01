C*==cpamills.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CPAMILLS(CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,NXM_Q,
     &                    NOQ,ITOQ,CONC,TSST,MSST,TSSQ,MSSQ,TAUQ,NTMAX,
     &                    NQMAX,M)
C   ********************************************************************
C   *                                                                  *
C   *   perform  CPA-iteration according    MILLS's  algorithm         *
C   *                                                                  *
C   *   the CPA - iteration step for site IQ is omitted if             *
C   *   ICPA(IQ) = 0 ( set in <INITALL> )                              *
C   *                                                                  *
C   *   only the projection matrix  DTILT is needed                    *
C   *                                                                  *
C   *   allows an atom type IT to have different orientation of        *
C   *   its moment on different but equivalent sites  IQ               *
C   *                                                                  *
C   *   NOTE: all matrices *Q and *T refer to the GLOBAL frame         *
C   *                                                                  *
C   *   NOTE: the Mills algorithm iterates MSSQ                        *
C   *         at the end of the routine TSSQ = 1/MSSQ is also updated  *
C   *                                                                  *
C   *   NOTE:   for IREL <> 2:     N = NKM   M = NKMMAX   NXM_Q = NKMQ *
C   *           for IREL =  2:     N = NLM   M = NLMMAX   NXM_Q = NLMQ *
C   *                           or N = NKM   M = NKMMAX   NXM_Q = NKMQ *
C   *                           depending on calling routine           *
C   *                                                                  *
C   * 13/12/2014 HE                                                    *
C   ********************************************************************
C
C
      USE MOD_ANGMOM,ONLY:IPIVKM
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_THERMAL,ONLY:X_VFT,UMAT_VT,NVIBRA,NFLUCT,FMAT_FT
c modified by XJQ: scf of vibrations
      use mod_scfvb_cpa_sigma,only:lscfvb
c end-mod-xjq
      IMPLICIT NONE
C*--CPAMILLS36
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CPAMILLS')
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPACORR,CPAERR
      INTEGER IPRINT,M,NQ,NQMAX,NTMAX
      REAL*8 CONC(NTMAX)
      INTEGER ICPA(NQMAX),ITOQ(NTMAX,NQMAX),NOQ(NQMAX),NXM_Q(NQMAX)
      COMPLEX*16 MSSQ(M,M,NQMAX),MSST(M,M,NTMAX),TAUQ(M,M,NQMAX),
     &           TSSQ(M,M,NQMAX),TSST(M,M,NTMAX)
C
C Local variables
C
c modified by XJQ: scf of vibrations
      integer ivt_extend
c end-mod-xjq
      COMPLEX*16 CSUM,DMAMC(:,:),DMATT(:,:),DQ(:,:),DTILT(:,:),ERR(:,:),
     &           MSS_VFT(:,:),TSS_FT(:,:),TSS_VFT(:,:),WA(:,:),WB(:,:)
      INTEGER I,IA_ERR,IFLUCT,IFT,IO,IQ,IT,IVFT,IVIBRA,IVT,J,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMAMC,DQ,DMATT,DTILT,ERR,TSS_FT,TSS_VFT,MSS_VFT,WA,WB
C
      ALLOCATE (TSS_FT(M,M))
      ALLOCATE (WA(M,M),WB(M,M))
      ALLOCATE (TSS_VFT(M,M),MSS_VFT(M,M))
      ALLOCATE (DMAMC(M,M))
      ALLOCATE (DQ(M,NQMAX),DMATT(M,M))
      ALLOCATE (DTILT(M,M),ERR(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ERR')
C
      CPAERR = 0.0D0
      CPACORR = 0.0D0
      CPACHNG = 0.0D0
C
C================================================================= IQ ==
      LOOP_IQ:DO IQ = 1,NQ
C
         N = NXM_Q(IQ)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( ICPA(IQ).EQ.0 ) CYCLE LOOP_IQ
C
         DQ(1:N,IQ) = -C1
         ERR(1:N,1:N) = C0
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
C------------------------------------------------ perform local rotation
C----------- keep intermediate t-matrix TSS_FT: t_ft = F_f * t_t * F_f^+
C
                  IF ( NFLUCT.EQ.1 ) THEN
C
                     TSS_FT(:,:) = TSST(:,:,IT)
C
                  ELSE
C
                     CALL CMATTRANS(N,M,WA,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSST(1,1,IT),'SAS+',TSS_FT)
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
     &                      *NFLUCT + IFLUCT
C
C-------------------------------- TSS_VFT: t_vft = U_v * t_ft * U_v^(-1)
C
                     IF ( NVIBRA.EQ.1 ) THEN
C
                        TSS_VFT(:,:) = TSS_FT(:,:)
C
                     ELSE
C
c modified by XJQ: scfvb
                        call read_scftv(IVT,TSS_FT)
c end-mod-xjq
                        CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,
     &                                    'UAUT',TSS_VFT)
C
                     END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                     CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WA,MSS_VFT)
C
C---------------------------------- get projection matrices DMATT, DTILT
C
                     CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,
     &                            MSSQ(1,1,IQ),MSS_VFT,M)
C
                     DO I = 1,N
                        DQ(I,IQ) = DQ(I,IQ) + X_VFT(IVFT)*DTILT(I,I)
                     END DO
C
C     -------------------------------------------
C            - E[a] = D~[a] * ( m[a] - m[c] )
C     -------------------------------------------
                     CALL CMATMUL(N,M,DTILT,DMAMC,WA)
C
C     -------------------------------------------
C            E = SUM[a]  c[a] *  E[a]
C     -------------------------------------------
C
                     ERR(1:N,1:N) = ERR(1:N,1:N) - X_VFT(IVFT)
     &                              *WA(1:N,1:N)
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
     &                          nvibra
                 else
                   ivt_extend = (it/nvibra) * nvibra*nvibra +
     &                          mod(it,nvibra)
                 endif
                 CALL CMAT_U_TRANS(UMAT_VT(1,1,ivt_extend),
     &                             TSST(1,1,IT),'UAUT',TSS_VFT)
                 CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WA,MSS_VFT)
c
                 CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,
     &                        MSSQ(1,1,IQ),MSS_VFT,M)
               else
                 CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,
     &                        MSSQ(1,1,IQ),MSST(1,1,IT),M)
               endif
c end-mod-xjq
C
               DO I = 1,N
                  DQ(I,IQ) = DQ(I,IQ) + CONC(IT)*DTILT(I,I)
               END DO
C
C     -------------------------------------------
C            - E[a] = D~[a] * ( m[a] - m[c] )
C     -------------------------------------------
               CALL CMATMUL(N,M,DTILT,DMAMC,WA)
C
C     -------------------------------------------
C            E = SUM[a]  c[a] *  E[a]
C     -------------------------------------------
C
               ERR(1:N,1:N) = ERR(1:N,1:N) - CONC(IT)*WA(1:N,1:N)
C
            END IF
C================================= standard mode: NO thermal effects END
C
         END DO LOOP_IO
C================================================================= IO ==
C
C     -------------------------------------------
C                   E * TAU
C     -------------------------------------------
C
         CALL CMATMUL(N,M,ERR,TAUQ(1,1,IQ),WA)
C
C     -------------------------------------------
C                1 + E * TAU
C     -------------------------------------------
         DO I = 1,N
            WA(I,I) = C1 + WA(I,I)
         END DO
C
C     -------------------------------------------
C               ( 1 + E * TAU )**(-1)
C     -------------------------------------------
C
         CALL CMATINV(N,M,WA,WB)
C
C     -------------------------------------------
C           ( 1 + E * TAU )**(-1) * E
C     -------------------------------------------
C
         CALL CMATMUL(N,M,WB,ERR,WA)
C
C     -------------------------------------------
C           APPLY CORRECTION  TO  MEFF
C     m{n+1} = m{n} -  ( 1 + E * TAU )**(-1) * E
C     -------------------------------------------
C
         DO J = 1,N
C
            CPAERR = CPAERR + ABS(DREAL(DQ(J,IQ)))
     &               + ABS(DIMAG(DQ(J,IQ)))
            CPACORR = CPACORR + ABS(DREAL(WA(J,J)))
     &                + ABS(DIMAG(WA(J,J)))
            CPACHNG = DBLE(MAX(CPACHNG,ABS(WA(J,J)/MSSQ(J,J,IQ))))
C
            DO I = 1,N
               MSSQ(I,J,IQ) = MSSQ(I,J,IQ) - WA(I,J)
            END DO
C
         END DO
C
         CPAERR = CPACHNG
C
         IF ( IPRINT.GE.2 ) THEN
            CSUM = C0
            DO I = 1,N
               CSUM = CSUM + MSSQ(I,I,IQ)
            END DO
            WRITE (6,99001) IQ,CPAERR,CPACORR,CSUM
         END IF
C
C--------------------------------------invert to get TSSQ: t[q] = 1/m[q]
C
         CALL CMATINV3(N,M,IPIVKM,MSSQ(1,1,IQ),WA,TSSQ(1,1,IQ))
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END DO LOOP_IQ
C================================================================= IQ ==
C
C       stop
      DEALLOCATE (DMAMC,DQ,DMATT,DTILT,ERR,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (' CPA:  IQ',I3,'  ERR',F12.5,'  CORR',F13.5,'  M',
     &        18(1X,2(1PE14.6)))
      END
