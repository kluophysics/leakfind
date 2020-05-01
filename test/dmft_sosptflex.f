C*==dmft_sosptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_SOSPTFLEX(GF,SIGM,UR,NDEG,NREP,NMSLIN,EGM,TEMP,DC,
     &                          IPRINT,SIGSTAT)
C NDEG spin-orbit index (1..10)
C nrep nrep=ndeg/nspn
C nmslin number of matsubara/2
C temp temperature in units of Ur
C dc integer to specify double counting 0..3
C GF Bath Green's funtion:
C                      on matsubara positive
C SIGM  output self energy
C  U(2l+1,2l+1,2l+1,2l+1)
C***********************************************************
C********** SP-T-FLEX  for  N-band spin-polarised **********
C********** SP-Tmatrix: Hartree and Fock          **********
C********** FLEX-local  ONLY P-H channel          **********
C********** with renormalized U from P-P channel  **********
C********** General complex Hamiltonian Ho        **********
C********** A.Lichtenstein (Nijmegen)             **********
C**********   In collaboration with:              **********
C********** M.Katsnelson (IFM)                    **********
C***********************************************************
C Spin-orbit version (L. Pourovskii)
C GF(0:NMSLIN,2,2,NDEG) is diagonal in orbital quantum number part
C of GF(omega)
C Ur is U-matrix: <1,2|U|3,4> ,where 1,2,3,4 are orbital indexes
C========+=========+=========+=========+=========+=========+=========+=$
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NLM,NN,NNNN,NLMS,NS,NOM,PI,OMEGA,TIME,
     &    DOMEGA,DTIME,MINOM,UC,UW0,IPRT
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C*--DMFT_SOSPTFLEX33
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DC,IPRINT,NDEG,NMSLIN,NREP
      REAL*8 TEMP
      REAL*8 EGM(3),SIGSTAT(NDEG,NDEG),UR(NREP,NREP,NREP,NREP)
      COMPLEX*16 GF(NDEG,NDEG,0:NMSLIN),SIGM(NDEG,NDEG,0:NMSLIN)
C
C Local variables
C
      REAL*8 FNMM
      COMPLEX*16 G(:,:,:),SIG(:,:,:)
      INTEGER IOM,M1
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE G,SIG
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C
C
      FNMM = LOG(DBLE(NMSLIN))/LOG(2D0)
      IF ( ABS(FNMM-NINT(FNMM)).GT.1D-6 ) THEN
         WRITE (6,99001)
         STOP
      END IF
C
      NOM = 2*NMSLIN
      NLM = NDEG/2
      NS = 2
      NLMS = NDEG
      NN = NLMS*NLMS
      NNNN = NN*NN
      IPRT = IPRINT
C
      ALLOCATE (UC(NLMS,NLMS,NLMS,NLMS))
      ALLOCATE (UW0(NLMS,NLMS,NLMS,NLMS))
      ALLOCATE (G(NLMS,NLMS,NOM),SIG(NLMS,NLMS,NOM))
      ALLOCATE (OMEGA(NOM),TIME(NOM),MINOM(NOM))
C
      IF ( IPRT.GT.0 ) THEN
         PRINT *,'SP-T-FLEX+SO:  N=',NLM
C
         PRINT *,'T(EV,K) =',DBLE(TEMP),DBLE(TEMP*11605*RY_EV)
      END IF
      DOMEGA = 2D0*TEMP*PI
      DTIME = 1.D0/NOM/TEMP
C
      CALL DMFT_MATSUB_SOSPTFLEX
C
      CALL DMFT_VERTEXSP_SOSPTFLEX(UR)
      IF ( IPRT.GT.0 ) THEN
         PRINT *,'dT,dOm=',DBLE(DTIME),DBLE(DOMEGA)
         WRITE (6,'(14F9.4)') (DBLE(UC(M1,M1,M1,M1)),M1=1,NLMS)
         PRINT *,' nlm=',NLM,' nn=',NN
      END IF
C
C-------------------------------------------------
      SIG(:,:,:) = DCMPLX(0.0D0,0.0D0)
      G(:,:,:) = DCMPLX(0.0D0,0.0D0)
      DO IOM = 1,NOM/2
         G(:,:,IOM*2) = GF(:,:,IOM-1)
      END DO
      IF ( IPRT.GT.0 ) THEN
C
         DO IOM = 1,NOM/2
            WRITE (46,'(29f12.6)') OMEGA(IOM) + DOMEGA/2D0,
     &                             (G(M1,M1,IOM*2),M1=1,NLMS)
         END DO
      END IF
C
C----------Starting FLEX            --------------
C-------------- Gf
      IF ( IPRT.GT.0 ) PRINT *,'NS,NLM,NOM =',NS,NLM,NOM
C
C------------------------------
C
      CALL DMFT_SIG_SOSPTFLEX(G,SIG,EGM,TEMP,DC,SIGSTAT)
C------------------------------
C---------- Matrix inversion
C        G^-1=Go^-1 - Sig
C-------------- GFF - note that G - was overwrited...
      DO IOM = 1,NOM/2
         SIGM(:,:,IOM-1) = SIG(:,:,IOM*2)
         GF(:,:,IOM-1) = G(:,:,IOM*2)
      END DO
C
      IF ( IPRT.GT.0 ) THEN
C
         DO IOM = 1,NOM/2
            WRITE (41,'(29f12.6)') OMEGA(IOM) + DOMEGA/2D0,
     &                             (SIG(M1,M1,IOM*2),M1=1,NLMS)
         END DO
      END IF
      DEALLOCATE (UC,UW0)
      DEALLOCATE (OMEGA,TIME,MINOM)
      DEALLOCATE (SIG,G)
C
99001 FORMAT ('SPTFLEX: NUMBER OF MATSUBARA POINTS (NMATSUB) ',
     &        'MUST BE POWER OF 2')
      END
