C*==gil0.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GIL0(MREG_FT,TSST,MSSQ,TAUQ,ALF0Q,ALF0QO)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *                J_m Im TAU00(t) J_n Im TAU00(t)                   *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_TYPES,ONLY:NTMAX,CONC
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,WKM1,WKM2,IPIVKM
      USE MOD_SITES,ONLY:NQMAX,NOMAX,NOQ,ITOQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_THERMAL,ONLY:X_VFT,UMAT_VT,NVIBRA,NFLUCT,FMAT_FT,NFTMAX
      IMPLICIT NONE
C*--GIL019
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='GIL0')
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      COMPLEX*16 ALF0Q(3,3,NQMAX),ALF0QO(3,3,NQMAX,NOMAX),
     &           MREG_FT(NKMMAX,NKMMAX,3,NFTMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ALF_VFT,IM_TAU_VFT(:,:),MREG_VFT(:,:,:),MSS_VFT(:,:),
     &           TAU_VFT(:,:),TSS_FT(:,:),TSS_VFT(:,:)
      COMPLEX*16 CMATTRC
      INTEGER I,IA_ERR,IFLUCT,IFT,IO,IQ,IT,IVFT,IVIBRA,IVT,J,M,MUE,N,NUE
      CHARACTER*5 TXTVAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAU_VFT,IM_TAU_VFT,MREG_VFT
      ALLOCATABLE TSS_FT,TSS_VFT,MSS_VFT
C
      CALL TRACK_INFO(ROUTINE)
C
      M = NKMMAX
C
      ALLOCATE (MREG_VFT(M,M,3),MSS_VFT(M,M))
      ALLOCATE (TSS_VFT(M,M),TSS_FT(M,M))
      ALLOCATE (TAU_VFT(M,M),IM_TAU_VFT(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: IM_TAU')
C
      TXTVAR = 'alpha'
C
      WRITE (6,99001)
C
      ALF0Q(:,:,:) = 0D0
      ALF0QO(:,:,:,:) = 0D0
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      DO IQ = IQBOT_CHI,IQTOP_CHI
         IF ( IPRINT.GT.0 ) WRITE (6,99002) IQ
C
         N = NKMQ(IQ)
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
C-------------------------------------------------- TAU(t) = TAU * D~(t)
C
                  CALL GET_TAUT(MSS_VFT,MSSQ(1,1,IQ),TAUQ(1,1,IQ),
     &                          TAU_VFT)
C
C-----------------------------------------------------------------------
C
                  DO J = 1,N
                     DO I = 1,N
                        IM_TAU_VFT(I,J)
     &                     = (TAU_VFT(I,J)-DCONJG(TAU_VFT(J,I)))
     &                     /(2D0*CI)
C
                     END DO
                  END DO
C
C----------------------------------------------------- j_m Im G j_n Im G
C
                  DO MUE = 1,3
                     DO NUE = 1,3
C
                        CALL CMATMUL(N,M,MREG_VFT(1,1,NUE),IM_TAU_VFT,
     &                               WKM1)
C
                        CALL CMATMUL(N,M,IM_TAU_VFT,WKM1,WKM2)
C
                        CALL CMATMUL(N,M,MREG_VFT(1,1,MUE),WKM2,WKM1)
C
                        ALF_VFT = CMATTRC(N,M,WKM1)
C
                        ALF0Q(MUE,NUE,IQ) = ALF0Q(MUE,NUE,IQ)
     &                     + X_VFT(IVFT)*ALF_VFT
C
                        ALF0QO(MUE,NUE,IQ,IO) = ALF0QO(MUE,NUE,IQ,IO)
     &                     + (X_VFT(IVFT)/CONC(IT))*ALF_VFT
C
                     END DO
                  END DO
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
C ---------------------------------------------- suppress small elements
C
         DO MUE = 1,3
            DO NUE = 1,3
               IF ( ABS(ALF0Q(MUE,NUE,IQ)).LT.TOL ) ALF0Q(MUE,NUE,IQ)
     &              = C0
C
               IF ( ABS(DIMAG(ALF0Q(MUE,NUE,IQ))).LT.TOL )
     &              ALF0Q(MUE,NUE,IQ) = DREAL(ALF0Q(MUE,NUE,IQ))
C
            END DO
         END DO
C
         IF ( IPRINT.GE.2 ) THEN
            WRITE (6,99003) TXTVAR
            DO MUE = 1,3
               WRITE (6,99005) (DREAL(ALF0Q(MUE,NUE,IQ)),NUE=1,3)
            END DO
         END IF
C
C ----------------------------- check if imaginary part of ALF0Q is zero
C
         DO MUE = 1,3
            DO NUE = 1,3
               IF ( ABS(DIMAG(ALF0Q(MUE,NUE,IQ))).GT.1D-5 )
     &              WRITE (6,99004) TXTVAR,DIMAG(ALF0Q(MUE,NUE,IQ)),MUE,
     &                              NUE,IQ
            END DO
         END DO
C
      END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      DEALLOCATE (IM_TAU_VFT,TAU_VFT,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (//,1X,79('*'),/,37X,'<GIL0>',/,1X,79('*'),/)
99002 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99003 FORMAT (/,10X,'site-diagonal term ',A,'_0',/)
99004 FORMAT (' WARNING!! Im(',A,') =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99005 FORMAT (10X,3F14.6)
      END
