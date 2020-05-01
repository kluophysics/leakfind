C*==linresp_kloop_dk.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_KLOOP_DK(ERYD,P,TAUQ,TAUQBZ,MSSQ)
C   ********************************************************************
C   *                                                                  *
C   *    PERFORM THE K-SPACE INTEGRAL USING SPECIAL POINTS             *
C   *                                                                  *
C   *    - run a loop over the k-points  KTAB and sum TAU(k)           *
C   *      for the irreducible wedges                                  *
C   *    - the k-points  KTAB  have weights  WKTAB  according the      *
C   *      symmetry of the system                                      *
C   *      KTAB and WKTAB are set up in  <KMESHS>                      *
C   *    - the full BZ is accounted for by applying the symmetry       *
C   *      rotations  DROT                                             *
C   *    - using BLAS routines for matrix inversion                    *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  TAUKLIN = TAU(K) is set up as a LINEAR array                    *
C   *------------------------------------------------------------------*
C   *  WKSUM is determined in <INIT_MOD_KSPACE>                        *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  the derivative d / dk is determinded by 2- or 4- point formula  *
C   *  setting the parameter  NP = 2 or 4                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_LINRESP,ONLY:CHIZ,CHI_DK,NZ12,NZ12MAX,ITTA,ITTB,ITTC,ITTD,
     &    ITTQ1,ITTQ2,NTTJ,JTT1,JTT2,WTTJ,NTKTKLIN
      USE MOD_SITES,ONLY:NQ,NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      USE MOD_CONSTANTS,ONLY:CI,C0,PI
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TAUIJ,ONLY:NKTABTAUIJ,WKTABTAUIJ,KTABTAUIJ
      IMPLICIT NONE
C*--LINRESP_KLOOP_DK37
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_KLOOP_DK')
      INTEGER NC,NP
      PARAMETER (NC=3,NP=2)
      REAL*8 K_STEP
      PARAMETER (K_STEP=0.001D0)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUQBZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX)
C
C Local variables
C
      COMPLEX*16 CSCL,CSUM1,CWK,CWORK(:),MAUX(:,:),MAUXLIN_IP(:),
     &           SUMQ(:,:,:),TAUKLIN(:),TAUKLIN_DK(:),TKTKQQ(:),
     &           TKTKQQ_DK(:,:),WT
      REAL*8 DKP(:),DKP2(2),DKP4(4),EHAT(3,3),KPDKVEC(3),WK,WKP(:),
     &       WKP2(2),WKP4(4)
      INTEGER I,I1,IA_ERR,IC,IK,INFO,IP,IPIV(:),IPROCK(:),IQ,J,JQ,JTT,
     &        K14I,K23J,L1,L2,L3,L4,M,N,NKMSQ
C
C*** End of declarations rewritten by SPAG
C
      DATA DKP2/ - 1D0, + 1D0/,WKP2/ - 0.5D0, + 0.5D0/
      DATA DKP4/ - 1.5D0, - 0.5D0, + 0.5D0, + 1.5D0/
      DATA WKP4/ + 1D0, - 27D0, + 27D0, - 1D0/
      DATA EHAT/1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0/
C
      ALLOCATABLE MAUX,MAUXLIN_IP,TAUKLIN,TAUKLIN_DK,SUMQ,TKTKQQ,IPIV
      ALLOCATABLE IPROCK,CWORK,WKP,DKP,TKTKQQ_DK
C
      ALLOCATE (SUMQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TKTKQQ(NTKTKLIN),TKTKQQ_DK(NTKTKLIN,3),IPIV(NKKR))
C
      M = NKKR
      ALLOCATE (TAUKLIN_DK(M*M))
      ALLOCATE (MAUX(M,M),MAUXLIN_IP(M*M),TAUKLIN(M*M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate TAUKLIN')
C
      IF ( IREL.EQ.2 ) CALL STOP_MESSAGE(ROUTINE,'IREL = 2')
      IF ( NZ12.NE.1 )
     &     CALL STOP_MESSAGE(ROUTINE,'NZ12 not set properly')
      IF ( NZ12.NE.1 )
     &      CALL STOP_MESSAGE(ROUTINE,'NZ12MAX not set properly')
C
C=======================================================================
C            parameters for numerical differentiation
C=======================================================================
      IF ( NP.NE.2 .AND. NP.NE.4 )
     &      CALL STOP_MESSAGE(ROUTINE,'NP.NE.2 .AND. NP.NE.4')
      ALLOCATE (WKP(NP),DKP(NP))
C
      IF ( NP.EQ.2 ) THEN
         WKP(1:NP) = WKP2(1:NP)/(K_STEP*2*PI/ALAT)
         DKP(1:NP) = DKP2(1:NP)*K_STEP
      ELSE
         WKP4(1:NP) = WKP4(1:NP)/24D0
         WKP(1:NP) = WKP4(1:NP)/(K_STEP*2*PI/ALAT)
         DKP(1:NP) = DKP4(1:NP)*K_STEP
      END IF
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., NKTAB
C      over processors;   IK=NKTAB  dealt with by IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCK(NKTABTAUIJ))
C
      CALL MPI_DISTRIBUTE(IPROCK,NKTABTAUIJ,MPI_KLOOP,'K')
C
C ------------ calculate energy - dependent terms of structure constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      TKTKQQ(1:NTKTKLIN) = C0
      TKTKQQ_DK(1:NTKTKLIN,:) = C0
C
C------------------------- assume the same L-expansion for every site IQ
      NKM = NKMQ(1)
      DO IQ = 2,NQ
         IF ( NKMQ(IQ).NE.NKM )
     &         CALL STOP_MESSAGE(ROUTINE,'NKMQ(IQ)<>NKM')
      END DO
C
      SUMQ(:,:,:) = C0
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      DO IK = 1,NKTABTAUIJ
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         IF ( MPI_ID.EQ.IPROCK(IK) .OR. MPI_ELOOP ) THEN
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            WK = WKTABTAUIJ(IK)
            CWK = DCMPLX(WK,0D0)
C
C-------------------------------------------------------- calculate TAUA
C
            CALL STRSET(IK,KTABTAUIJ(1,IK),MAUX,TAUKLIN,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUKLIN,MSSQ,NQMAX,NKKR,NKMMAX)
C
            CALL ZGETRF(NKKR,NKKR,TAUKLIN,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUKLIN,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C----------------------------------------------------------- store TAUQ
C
            DO IQ = IQBOT_CHI,IQTOP_CHI
               I1 = IND0Q(IQ)*NKKR + IND0Q(IQ) + 1 - NKKR
               N = NKMQ(IQ)
               DO J = 1,N
                  I1 = I1 + NKKR
                  CALL ZAXPY(N,CWK,TAUKLIN(I1),1,SUMQ(1,J,IQ),1)
               END DO
            END DO
C
C------------------------------------------------------- store TAUA*TAUA
C
            JTT = 0
            DO I = 1,NTKTKLIN
               CSUM1 = C0
               DO J = 1,NTTJ(I)
                  JTT = JTT + 1
                  WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
                  CSUM1 = CSUM1 + WT*TAUKLIN(JTT2(JTT))
               END DO
               TKTKQQ(I) = TKTKQQ(I) + CSUM1*WK
            END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            DO IC = 1,NC
C
               TAUKLIN_DK(:) = C0
C
               DO IP = 1,NP
C
C------------------------------------------------------ calculate TAU_DK
C
                  KPDKVEC(1:3) = KTABTAUIJ(1:3,IK) + EHAT(1:3,IC)
     &                           *DKP(IP)
C
                  CALL STRSET(0,KPDKVEC,MAUX,MAUXLIN_IP,P)
C
                  CALL SETKKR(NQ,NKMQ,IND0Q,MAUXLIN_IP,MSSQ,NQMAX,NKKR,
     &                        NKMMAX)
C
                  CALL ZGETRF(NKKR,NKKR,MAUXLIN_IP,NKKR,IPIV,INFO)
                  CALL ZGETRI(NKKR,MAUXLIN_IP,NKKR,IPIV,MAUX,NKKR*NKKR,
     &                        INFO)
C
                  TAUKLIN_DK(:) = TAUKLIN_DK(:) + WKP(IP)*MAUXLIN_IP(:)
C
               END DO
C
C------------------------------------------------------- store TAUA*TAUB
C
               JTT = 0
               DO I = 1,NTKTKLIN
                  CSUM1 = C0
                  DO J = 1,NTTJ(I)
                     JTT = JTT + 1
                     WT = WTTJ(JTT)*TAUKLIN(JTT1(JTT))
                     CSUM1 = CSUM1 + WT*TAUKLIN_DK(JTT2(JTT))
                  END DO
                  TKTKQQ_DK(I,IC) = TKTKQQ_DK(I,IC) + CSUM1*WK
               END DO
C
            END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
         END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_KLOOP ) THEN
C
         CALL DRV_MPI_BARRIER
C
C        use TAUQ as work space for transfer
         M = NKMMAX*NKMMAX*NQMAX
C
         CALL DRV_MPI_REDUCE_C(SUMQ(1,1,1),TAUQ(1,1,1),M)
C
         M = NTKTKLIN
         ALLOCATE (CWORK(NTKTKLIN))
         CALL DRV_MPI_REDUCE_C(TKTKQQ(1),CWORK(1),M)
         DO IC = 1,NC
            CALL DRV_MPI_REDUCE_C(TKTKQQ_DK(1,IC),CWORK(1),M)
         END DO
         DEALLOCATE (CWORK)
C
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_ID.EQ.0 .OR. MPI_ELOOP ) THEN
C
         TAUQ(:,:,:) = SUMQ(:,:,:)
C
         CSCL = 1D0
C
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO J = 1,NKMQ(IQ)
               CALL ZSCAL(NKMQ(IQ),CSCL,TAUQ(1,J,IQ),1)
            END DO
         END DO
C
         CSCL = 1D0
         CALL ZSCAL(NTKTKLIN,CSCL,TKTKQQ(1),1)
C
         CSCL = 1D0/CI
         TKTKQQ_DK(1:NTKTKLIN,1:3) = CSCL*TKTKQQ_DK(1:NTKTKLIN,1:3)
C
C-----------------------------------------------------------------------
C                   set up the matrix   TAUQBZ
C-----------------------------------------------------------------------
C
         DO IQ = IQBOT_CHI,IQTOP_CHI
            DO J = 1,NKMQ(IQ)
               DO I = 1,NKMQ(IQ)
                  TAUQBZ(I,J,IQ,1) = TAUQ(I,J,IQ)
               END DO
            END DO
         END DO
C
C-----------------------------------------------------------------------
C               set up the super-matrix  CHI(K14I,k2)
C-----------------------------------------------------------------------
C NOTE: indexing of CHI in <SIGKLOOPS>, <LINRESP_VERTEX> and <SIG1>
C       as well as of  w  in <LINRESP_VERTEX>
C       have to be consistent with loop sequence (L1,L4,L2,L3)
C-----------------------------------------------------------------------
         NKMSQ = NKM*NKM
C
         CHIZ(:,:,:) = C0
         CHI_DK(:,:,:) = C0
C
         DO I = 1,NTKTKLIN
C
            L1 = ITTA(I)
            L2 = ITTB(I)
            L3 = ITTC(I)
            L4 = ITTD(I)
            IQ = ITTQ1(I)
            JQ = ITTQ2(I)
C
            K14I = (L1-1)*NKM + L4 + (IQ-1)*NKMSQ
            K23J = (L2-1)*NKM + L3 + (JQ-1)*NKMSQ
C
C-----------------------------------------------------------------------
C
            IF ( IQ.EQ.JQ ) THEN
               CHIZ(K14I,K23J,1) = TKTKQQ(I) - TAUQ(L1,L2,IQ)
     &                             *TAUQBZ(L3,L4,IQ,1)
            ELSE
               CHIZ(K14I,K23J,1) = TKTKQQ(I)
            END IF
C
C-----------------------------------------------------------------------
C
            CHI_DK(K14I,K23J,1:3) = TKTKQQ_DK(I,1:3)
C
         END DO
C
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( MPI_KLOOP ) CALL DRV_MPI_BARRIER
C
      DEALLOCATE (MAUX,TAUKLIN,SUMQ,TKTKQQ,TAUKLIN_DK)
C
      END
