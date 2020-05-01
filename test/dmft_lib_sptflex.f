C*==dmft_sig_sptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_SIG_SPTFLEX(GOM,SOM,TEMP,MINOM,OMEGA,DTIME,DOMEGA,
     &                            UC,DMFTDBLC,UUU,UJJ,NLM,NLMSQ,NLMQU,
     &                            NOM,NSPIN,IPRINT,IFLEX)
C
      IMPLICIT NONE
C*--DMFT_SIG_SPTFLEX7
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*4 DMFTDBLC
      REAL*8 DOMEGA,DTIME,TEMP,UJJ,UUU
      INTEGER IFLEX,IPRINT,NLM,NLMQU,NLMSQ,NOM,NSPIN
      COMPLEX*16 GOM(NLM,NLM,NOM,NSPIN),SOM(NLM,NLM,NOM,NSPIN),
     &           UC(NLMSQ,NLMSQ)
      INTEGER MINOM(NOM)
      REAL*8 OMEGA(NOM)
C
C Local variables
C
      COMPLEX*16 CSUM,PDD0(:,:,:,:,:),PDZ0(:,:,:,:,:),PPM0(:,:,:,:,:),
     &           PPP(:,:,:,:,:,:),PPP0(:,:,:,:,:),TDD(:,:),TDM(:,:),
     &           TMAT(:,:,:,:,:,:),TMM(:,:),TPM(:,:),VDM(:,:,:,:,:,:),
     &           VPM(:,:,:,:,:),VVPM
      REAL*8 DCLDA(:,:,:),OCCM(:,:),RSUM,TIMEE,TIMES
      INTEGER IOM,IS,LM,M1,M2,M3,M4,MIOM
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATABLE OCCM,PDD0,PDZ0,PPM0,PPP0,PPP
      ALLOCATABLE TMAT,TDD,TMM,TDM,TPM
      ALLOCATABLE VDM,VPM,DCLDA
C=================================================================
C     Calculate occupation matrix if LDA+U double counting
C=================================================================
      CALL CPU_TIME(TIMES)
C
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
         ALLOCATE (OCCM(0:NLM,NSPIN))
         OCCM = 0.0D0
         DO IS = 1,NSPIN
            DO LM = 1,NLM
               DO IOM = 1,NOM/2
                  OCCM(LM,IS) = OCCM(LM,IS)
     &                          + TEMP*DREAL(GOM(LM,LM,2*IOM,IS))
               END DO
               OCCM(LM,IS) = OCCM(LM,IS) + 0.5D0
                                !Fermi step
               OCCM(0,IS) = OCCM(0,IS) + OCCM(LM,IS)
            END DO
         END DO
C
         IF ( IPRINT.GT.0 ) THEN
            DO IS = 1,NSPIN
               WRITE (6,'(/,A,I1)') 'Occuparion matrix for Spin',IS
               WRITE (6,'(11f10.5)') (OCCM(LM,IS),LM=1,NLM)
            END DO
            WRITE (6,'(11f10.5)') OCCM(0,1),OCCM(0,2)
         END IF
      END IF
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for OCCM',TIMEE - TIMES
C
C=================================================================
C   obtain g(tau)
C=================================================================
      CALL CPU_TIME(TIMES)
      CALL DMFT_FFT_FTOT2(GOM,DOMEGA,NLMSQ,NOM,NSPIN)
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for G(TAU)',TIMEE - TIMES
C
      CALL CPU_TIME(TIMES)
C=================================================================
C   determine susceptibility Pzz and Pud (up-down)
C   Note that Pud, Ppp not depends on is, but we use both is=1,2
C=================================================================
      ALLOCATE (PDD0(NLM,NLM,NLM,NLM,NOM))
      ALLOCATE (PDZ0(NLM,NLM,NLM,NLM,NOM))
      ALLOCATE (PPM0(NLM,NLM,NLM,NLM,NOM))
C
      DO IOM = 1,NOM
         MIOM = MINOM(IOM)
         DO M1 = 1,NLM
            DO M2 = 1,NLM
               DO M3 = 1,NLM
                  DO M4 = 1,NLM
C-------------density d
                     PDD0(M1,M2,M3,M4,IOM)
     &                  = -0.5D0*(GOM(M4,M1,MIOM,1)*GOM(M2,M3,IOM,1)
     &                  +GOM(M4,M1,MIOM,2)*GOM(M2,M3,IOM,2))
C-------------magnetic -z
                     PDZ0(M1,M2,M3,M4,IOM)
     &                  = -0.5D0*(GOM(M4,M1,MIOM,1)*GOM(M2,M3,IOM,1)
     &                  -GOM(M4,M1,MIOM,2)*GOM(M2,M3,IOM,2))
C-------------magnetic +- = pm. Note spin-ordering!
                     PPM0(M1,M2,M3,M4,IOM) = -GOM(M4,M1,MIOM,2)
     &                  *GOM(M2,M3,IOM,1)
                  END DO
               END DO
            END DO
         END DO
      END DO
C=================================================================
C          Susceptibility Pph
C=================================================================
C
      CALL DMFT_FFT_TTOF4M(PDD0,DTIME,NLMQU,NOM)
      CALL DMFT_FFT_TTOF4M(PDZ0,DTIME,NLMQU,NOM)
      CALL DMFT_FFT_TTOF4M(PPM0,DTIME,NLMQU,NOM)
C
C=================================================================
C-------------  Prticle-Particle - spins!
C=================================================================
      ALLOCATE (PPP0(NLM,NLM,NLM,NLM,NOM))
      ALLOCATE (PPP(NLM,NLM,NLM,NLM,NOM,NSPIN+1))
      DO IS = 1,NSPIN + 1
         DO IOM = 1,NOM
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  DO M3 = 1,NLM
                     DO M4 = 1,NLM
                        IF ( IS.LE.NSPIN ) THEN
                           PPP0(M1,M2,M3,M4,IOM) = GOM(M1,M3,IOM,IS)
     &                        *GOM(M2,M4,IOM,IS)
                        ELSE
                           PPP0(M1,M2,M3,M4,IOM) = GOM(M1,M3,IOM,1)
     &                        *GOM(M2,M4,IOM,2)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C=================================================================
C-------------Ppp(0) - static susc. for different spins:
C=================================================================
C
         CALL DMFT_FFT_TTOF4M(PPP0,DTIME,NLMQU,NOM)
         DO IOM = 1,NOM
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  DO M3 = 1,NLM
                     DO M4 = 1,NLM
                        PPP(M1,M2,M3,M4,IOM,IS) = PPP0(M1,M2,M3,M4,IOM)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for PPP',TIMEE - TIMES
C
C=================================================================
C----------  Screening  U ------------
C=================================================================
      CALL CPU_TIME(TIMES)
      ALLOCATE (TMAT(NLM,NLM,NLM,NLM,NOM,NSPIN+1))
      ALLOCATE (TDD(NLMSQ,NLMSQ),TMM(NLMSQ,NLMSQ),TDM(NLMSQ,NLMSQ))
      ALLOCATE (TPM(NLMSQ,NLMSQ))
      CALL DMFT_TMATR(UC,PPP,TMAT,TDD,TMM,TDM,TPM,NLMSQ,NLM,NOM,NSPIN,
     &                IPRINT)
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for TMAT',TIMEE - TIMES
C
C
C=================================================================
C IF IFLEX.EQ.0 Only TMA-approximation is used
C=================================================================
C
      IF ( IFLEX.EQ.1 ) THEN
C
C=================================================================
C---------- FLEX effective chanal potential: - -------
C---------- Vph spin-pol   -------
C=================================================================
         CALL CPU_TIME(TIMES)
C
         ALLOCATE (VDM(NLM,NLM,NLM,NLM,NOM,NSPIN))
         ALLOCATE (VPM(NLM,NLM,NLM,NLM,NOM))
C
         CALL DMFT_TOTPOT(PDD0,PDZ0,PPM0,VDM,VPM,TDD,TMM,TDM,TPM,NLMSQ,
     &                    NOM,NSPIN,IPRINT)
C
C
C=================================================================
C     Effective potential V(tau)p
C=================================================================
         CALL DMFT_FFT_FTOT4(VDM,DOMEGA,NLMQU,NOM,NSPIN)
         CALL DMFT_FFT_FTOT4M(VPM,DOMEGA,NLMQU,NOM)
         CALL DMFT_FFT_FTOT4S(TMAT,DOMEGA,NLMQU,NOM,NSPIN)
C
         CALL CPU_TIME(TIMEE)
         IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for TOTPOT',TIMEE - TIMES
C=================================================================
C   determine Sigma(tau)= Tmat(HF)+FLEX(p-h)
C=================================================================
         CALL CPU_TIME(TIMES)
C
         DO IS = 1,NSPIN
            DO IOM = 1,NOM
               MIOM = MINOM(IOM)
               DO M1 = 1,NLM
                  DO M2 = 1,NLM
                     CSUM = (0.D0,0.D0)
                     DO M3 = 1,NLM
                        DO M4 = 1,NLM
C-------Since Vmp1342(tau)=Vpm4213(-tau) !
                           IF ( IS.EQ.1 ) THEN
                              VVPM = VPM(M1,M3,M4,M2,IOM)
                           ELSE
                              VVPM = VPM(M4,M2,M1,M3,MIOM)
                           END IF
                           CSUM = CSUM + TMAT(M1,M3,M2,M4,IOM,IS)
     &                            *GOM(M4,M3,MIOM,IS)
     &                            + TMAT(M1,M3,M2,M4,IOM,3)
     &                            *GOM(M4,M3,MIOM,3-IS)
     &                            - TMAT(M1,M4,M3,M2,IOM,IS)
     &                            *GOM(M3,M4,MIOM,IS)
     &                            + VDM(M1,M3,M4,M2,IOM,IS)
     &                            *GOM(M3,M4,IOM,IS)
     &                            + VVPM*GOM(M3,M4,IOM,3-IS)
                        END DO
                     END DO
                     SOM(M1,M2,IOM,IS) = CSUM
                  END DO
               END DO
            END DO
         END DO
         CALL CPU_TIME(TIMEE)
         IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for SOM',TIMEE - TIMES
C=================================================================
      ELSE   !IFLEX
         CALL DMFT_TMAOM(TMAT,GOM,SOM,MINOM,NLM,NSPIN,NOM)
      END IF  !IFLEX
C=================================================================
C   determine Sigma(om)
C=================================================================
      CALL DMFT_FFT_TTOF2(SOM,DTIME,NLMSQ,NOM,NSPIN)
C
C
      RSUM = 0D0
      DO IS = 1,NSPIN
         DO M1 = 1,NLM
            RSUM = RSUM + DREAL(GOM(M1,M1,1,IS))
         END DO
      END DO
C=================================================================
C   Substract double counting
C
C=================================================================
      CALL CPU_TIME(TIMES)
C
      ALLOCATE (DCLDA(NLM,NLM,NSPIN))
      DCLDA = 0.D0
      IF ( DMFTDBLC(1:4).EQ.'INSU' ) THEN
         DO M1 = 1,NLM
            DO M2 = 1,NLM
               DCLDA(M1,M2,1) = DREAL(SOM(M1,M2,2,1)+SOM(M1,M2,2,2))
     &                          /2.D0
               DCLDA(M1,M2,NSPIN) = DCLDA(M1,M2,1)
            END DO
         END DO
      ELSE IF ( DMFTDBLC(1:4).EQ.'META' ) THEN
         DO IS = 1,NSPIN
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  DCLDA(M1,M2,IS) = DREAL(SOM(M1,M2,2,IS))
               END DO
            END DO
         END DO
      ELSE IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
         DO IS = 1,NSPIN
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  DCLDA(M1,M2,IS) = 0.0D0 !DREAL(SOM(M1,M2,2,IS))
C                            dcLDA(m1,m2,is)=dcLDA(m1,m2,is)+
C     &                           (OCCM(M1,IS)-0.5D0)*UUU
C     &                           -(OCCM(M1,IS)-0.5D0)*UJJ
C                  IF ( M1.EQ.M2 ) DCLDA(M1,M2,IS) = DCLDA(M1,M2,IS)
C     &                 + (OCCM(M1,1)+OCCM(M1,2)-0.5D0)
C     &                 *UUU - (OCCM(M1,IS)-0.5D0)*UJJ
                  IF ( M1.EQ.M2 ) DCLDA(M1,M2,IS) = DCLDA(M1,M2,IS)
     &                 + (UUU-UJJ)*(OCCM(M1,IS)-0.5D0)
C
C
               END DO
            END DO
         END DO
      ELSE
         WRITE (6,*) DMFTDBLC(1:4),'Not yet implemented'
         STOP
      END IF
C
      DO IS = 1,NSPIN
         DO IOM = 1,NOM
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  SOM(M2,M1,IOM,IS) = SOM(M2,M1,IOM,IS)
     &                                - DCLDA(M2,M1,IS)
               END DO
            END DO
         END DO
      END DO
      IF ( IPRINT.GT.0 ) THEN
         WRITE (2,'(11f10.5)') (OMEGA(2*IOM),(SOM(M1,M1,2*IOM,1),M1=1,
     &                         NLM),IOM=1,NOM/2)
         WRITE (3,'(11f10.5)') (OMEGA(2*IOM),(SOM(M1,M1,2*IOM,2),M1=1,
     &                         NLM),IOM=1,NOM/2)
      END IF
C
      CALL CPU_TIME(TIMEE)
      IF ( IPRINT.GT.0 ) WRITE (6,*) 'Time for DBLC',TIMEE - TIMES
C
      END
C*==dmft_tmatr.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_TMATR(UC,PPP,TMAT,TDD,TMM,TDM,TPM,NN,NLM,NOM,NS,
     &                      IPRINT)
C
C***********************************************************
C
C     Calculates the T matrix in spirit of Kanamori Approach
C
C     U= U*[1+R0*U]^{-1}
C
C***********************************************************
C
      IMPLICIT NONE
C*--DMFT_TMATR345
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NLM,NN,NOM,NS
      COMPLEX*16 PPP(NN,NN,NOM,NS+1),TDD(NN,NN),TDM(NN,NN),
     &           TMAT(NN,NN,NOM,NS+1),TMM(NN,NN),TPM(NN,NN),UC(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(:,:),R0(:,:),UR(:,:),WORK(:,:)
      INTEGER I,INFO,IPIV(:,:),ISS,M,M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE D,R0,UR,IPIV,WORK
      ALLOCATE (D(NN,NN),R0(NN,NN),UR(NN,NN),IPIV(NN,NN))
      ALLOCATE (WORK(NN,NN))
C
C*** End of declarations rewritten by SPAG
C
C
C
C         COMPLEX(KIND=8), DIMENSION(nn,nn,nom,ns+1) :: Tmat, Ppp
C         COMPLEX(KIND=8), DIMENSION(nn,nn) :: Tdd, Tdm, Tpm, Tmm
C         COMPLEX(KIND=8), DIMENSION(nn,nn) :: R0, D
C         COMPLEX(KIND=8), DIMENSION(nn,nn) :: Uc, Ur
C         COMPLEX(KIND=8), DIMENSION(nn) :: WORK
C         INTEGER, DIMENSION(nn) :: IPIV
C----------- put 0 for Fermion Matsubara
      DO M1 = 1,NN
         DO M2 = 1,NN
            DO I = 1,NOM/2
C           do i=1,nom/4
               DO ISS = 1,NS + 1
                  TMAT(M1,M2,2*I,ISS) = DCMPLX(0.D0,0.D0)
               END DO
            END DO
         END DO
      END DO
C---------- Tmat - screening, spins: iss = uu, dd, ud (ud=du)------
      DO I = 1,NOM,2
C       do i=1,nom/2,2
         DO ISS = 1,NS + 1
            DO M1 = 1,NN
               DO M2 = 1,NN
                  R0(M1,M2) = PPP(M1,M2,I,ISS)
               END DO
            END DO
            IF ( I.EQ.1 .AND. IPRINT.EQ.1 ) THEN
               WRITE (6,*) 'Rpp0(0)diag: is=',ISS
               WRITE (6,'(5f14.5)') (DREAL(R0(M,M)),M=1,NN)
C       print*,'Rpp0(0)=',iss
C       write(6,'(25f4.1)')((dReal(R0(m1,m2)),m2=1,nn),m1=1,nn)
            END IF
C-----------------PRINT----------
C       print*,'Uc:'
C       write(6,'(25f4.1)')((dReal(Uc(m1,m2)),m2=1,nn),m1=1,nn)
C----------  PP-screening: Ur=U*[1+R0*U]^-1 --------
            DO M1 = 1,NN
               D(M1,M1) = DCMPLX(1.D0,0.D0)
               DO M2 = 1,NN
                  IF ( M2.NE.M1 ) D(M1,M2) = DCMPLX(0.D0,0.D0)
               END DO
            END DO
C-----------------PRINT----------
C       print*,'Re Ppp0:'
C       write(6,'(4f10.5)')((dReal(R0(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Im Ppp0:'
C       write(6,'(4f10.5)')((dImag(R0(m1,m2)),m2=1,nn),m1=1,nn)
C------- BLAS multiplication
            CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),R0,NN,UC,NN,
     &                 DCMPLX(1.D0,0.D0),D,NN)
C------- LAPACK inversion of general complex matrix
            CALL ZGETRF(NN,NN,D,NN,IPIV,INFO)
            IF ( INFO.NE.0 ) STOP 'in RPA ZGETRF Info NE 0'
C
            CALL ZGETRI(NN,D,NN,IPIV,WORK,NN,INFO)
            IF ( INFO.NE.0 ) STOP 'in RPA ZGETRI Info NE 0'
C
C------- BLAS multiplication
            CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),UC,NN,D,NN,
     &                 DCMPLX(0.D0,0.D0),UR,NN)
C-----------------PRINT----------
C       print*,'Ur:'
C       write(6,'(25f4.1)')((dReal(Ur(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Im Ur:'
C       write(6,'(4f10.5)')((dImag(Ur(m1,m2)),m2=1,nn),m1=1,nn)
            DO M1 = 1,NN
               DO M2 = 1,NN
                  TMAT(M1,M2,I,ISS) = UR(M1,M2)
               END DO
            END DO
C
         END DO
               ! iss
      END DO  ! iom
      CALL DMFT_SCRSYM(TMAT,TDD,TMM,TDM,TPM,NLM,NOM,NS)
C-----------------PRINT----------
      IF ( IPRINT.EQ.1 ) THEN
         DO ISS = 1,NS + 1
            WRITE (6,*) 'Tmat(ii):is=',ISS
            WRITE (6,'(5f8.3)') (DREAL(TMAT(M,M,1,ISS)),M=1,NN)
         END DO
              ! iss
      END IF
      END
C*==dmft_scrsym.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
      SUBROUTINE DMFT_SCRSYM(TMAT,TDD,TMM,TDM,TPM,NLM,NOM,NS)
      IMPLICIT NONE
C*--DMFT_SCRSYM470
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLM,NOM,NS
      COMPLEX*16 TDD(NLM,NLM,NLM,NLM),TDM(NLM,NLM,NLM,NLM),
     &           TMAT(NLM,NLM,NLM,NLM,NOM,NS+1),TMM(NLM,NLM,NLM,NLM),
     &           TPM(NLM,NLM,NLM,NLM)
C
C Local variables
C
      INTEGER M1,M2,M3,M4
C
C*** End of declarations rewritten by SPAG
C
C-------------- T- Vertex for different spin-channel channel
      DO M1 = 1,NLM
         DO M2 = 1,NLM
            DO M3 = 1,NLM
               DO M4 = 1,NLM
                  TDD(M1,M2,M3,M4) = 0.5D0*(TMAT(M1,M3,M2,M4,1,1)+TMAT(
     &                               M1,M3,M2,M4,1,2)
     &                               +2.D0*TMAT(M1,M3,M2,M4,1,3))
     &                               - 0.5D0*(TMAT(M1,M3,M4,M2,1,1)
     &                               +TMAT(M1,M3,M4,M2,1,2))
                  TMM(M1,M2,M3,M4) = 0.5D0*(TMAT(M1,M3,M2,M4,1,1)+TMAT(
     &                               M1,M3,M2,M4,1,2)
     &                               -2.D0*TMAT(M1,M3,M2,M4,1,3))
     &                               - 0.5D0*(TMAT(M1,M3,M4,M2,1,1)
     &                               +TMAT(M1,M3,M4,M2,1,2))
                  TDM(M1,M2,M3,M4) = 0.5D0*(TMAT(M1,M3,M2,M4,1,1)-TMAT(
     &                               M1,M3,M2,M4,1,2))
     &                               - 0.5D0*(TMAT(M1,M3,M4,M2,1,1)
     &                               -TMAT(M1,M3,M4,M2,1,2))
                  TPM(M1,M2,M3,M4) = -TMAT(M1,M3,M4,M2,1,3)
               END DO
            END DO
         END DO
      END DO
C
      END
C*==dmft_totpot.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C
C
C***********************************************************
      SUBROUTINE DMFT_TOTPOT(PDD0,PDZ0,PPM0,VDM,VPM,TDD,TMM,TDM,TPM,NN,
     &                       NOM,NS,IPRINT)
      IMPLICIT NONE
C*--DMFT_TOTPOT532
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,NN,NOM,NS
      COMPLEX*16 PDD0(NN,NN,NOM),PDZ0(NN,NN,NOM),PPM0(NN,NN,NOM),
     &           TDD(NN,NN),TDM(NN,NN),TMM(NN,NN),TPM(NN,NN),
     &           VDM(NN,NN,NOM,NS),VPM(NN,NN,NOM)
C
C Local variables
C
      INTEGER IOM,IS,M,M1,M2
      COMPLEX*16 RDD0(:,:),RDZ0(:,:),RPM(:,:),RPM0(:,:),UDM(:,:),
     &           VDD(:,:),VDZ(:,:),VZD(:,:),VZZ(:,:),WDM(:,:,:),WPM(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE VDD,UDM,WDM,RPM,VDZ,VZD,WPM,RDD0,VZZ,RPM0,RDZ0
      ALLOCATE (VDD(NN,NN),UDM(2*NN,2*NN),WDM(NN,NN,NS),RPM(NN,NN))
      ALLOCATE (VDZ(NN,NN),VZD(NN,NN),WPM(NN,NN),RDD0(NN,NN))
      ALLOCATE (VZZ(NN,NN),RPM0(NN,NN),RDZ0(NN,NN))
C
C*** End of declarations rewritten by SPAG
C
C-------- Effective interactions in D-Mz channels
      DO M1 = 1,NN
         DO M2 = 1,NN
            UDM(M1,M2) = TDD(M1,M2)
            UDM(M1+NN,M2+NN) = TMM(M1,M2)
            UDM(M1,M2+NN) = TDM(M1,M2)
            UDM(M1+NN,M2) = TDM(M1,M2)
C----------- put 0 for Fermion Matsubara
            DO IOM = 1,NOM/2
C           do iom=1,nom/4
               DO IS = 1,NS
                  VDM(M1,M2,2*IOM,IS) = DCMPLX(0.D0,0.D0)
               END DO
               VPM(M1,M2,2*IOM) = DCMPLX(0.D0,0.D0)
            END DO
         END DO
      END DO
C-----------------PRINT----------
C       print*,'Re Udm:'
C       write(6,'(8f10.4)')((Real(Udm(m1,m2)),m2=1,2*nn),m1=1,2*nn)
C---------- FLEX effective channels potential  --------
C---------- d-m 2x2 chanal  and only Bosonic Matsubara --------
      DO IOM = 1,NOM,2
C       do iom=1,nom/2,2
         DO M1 = 1,NN
            DO M2 = 1,NN
               RDD0(M1,M2) = PDD0(M1,M2,IOM)
               RDZ0(M1,M2) = PDZ0(M1,M2,IOM)
               RPM0(M1,M2) = PPM0(M1,M2,IOM)
            END DO
         END DO
C
C
         IF ( IPRINT.GT.0 ) THEN
            IF ( IOM.EQ.1 ) THEN
               WRITE (6,*) 'Rdd0(0)='
               WRITE (6,'(5f14.6)') (DREAL(RDD0(M,M)),M=1,NN)
C       write(6,'(4f15.5)')((Real(Rdd0(m1,m2)),m2=1,nn),m1=1,nn)
               WRITE (6,*) 'Rdz0(0)='
               WRITE (6,'(5f14.6)') (DREAL(RDZ0(M,M)),M=1,NN)
C       write(6,'(4f15.5)')((Real(Rdz0(m1,m2)),m2=1,nn),m1=1,nn)
               WRITE (6,*) 'Rpm0(0)='
               WRITE (6,'(5f14.6)') (DREAL(RPM0(M,M)),M=1,NN)
C       write(6,'(4f15.5)')((Real(Rpm0(m1,m2)),m2=1,nn),m1=1,nn)
            END IF
         END IF
C----------------------------------------
C----------- Vph M+-  channel  (transverse)
C----------- Vph - D-M-channel (longitudinal)
C
         CALL DMFT_RPASP(UDM,RDD0,RDZ0,VDD,VZZ,VDZ,VZD,NN)
C
C-----------------PRINT----------
C       print*,'Re Rdd:'
C       write(6,'(4f15.5)')((Real(Vdd(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Re Rdz:'
C       write(6,'(4f15.5)')((Real(Vdz(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Re Rzz:'
C       write(6,'(4f15.5)')((Real(Vzz(m1,m2)),m2=1,nn),m1=1,nn)
C----------------------------------------
C
         CALL DMFT_POTSP(VDD,VZZ,VDZ,VZD,WDM,NN,NS)
C----------- Vph M+-  channel  (transverse)
C
         CALL DMFT_RPA(TPM,RPM0,RPM,NN)
C
         CALL DMFT_POT(TPM,RPM0,RPM,WPM,NN)
C------ Total Vph Spin-polarized case
         DO M1 = 1,NN
            DO M2 = 1,NN
               DO IS = 1,NS
                  VDM(M1,M2,IOM,IS) = WDM(M1,M2,IS)
               END DO
               VPM(M1,M2,IOM) = WPM(M1,M2)
            END DO
         END DO
C
      END DO
C
      END
C*==dmft_rpasp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C**********************************************************
      SUBROUTINE DMFT_RPASP(UDM,PDD0,PDZ0,VDD,VZZ,VDZ,VZD,NN)
      IMPLICIT NONE
C*--DMFT_RPASP655
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NN
      COMPLEX*16 PDD0(NN,NN),PDZ0(NN,NN),UDM(2*NN,2*NN),VDD(NN,NN),
     &           VDZ(NN,NN),VZD(NN,NN),VZZ(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(:,:),RDM(:,:),RDM0(:,:),WORK(:,:)
      INTEGER INFO,IPIV(:),M1,M2,N2
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE D,RDM,RDM0,IPIV,WORK
      ALLOCATE (D(2*NN,2*NN),RDM(2*NN,2*NN),RDM0(2*NN,2*NN))
      ALLOCATE (IPIV(2*NN),WORK(2*NN,2*NN))
C
C*** End of declarations rewritten by SPAG
C
C
C
      N2 = 2*NN
C---------- Spin-polarized FLEX "RPA" d-m part 2x2  R=[1+R0*U]^-1*R0 --------
      DO M1 = 1,NN
         DO M2 = 1,NN
            RDM0(M1,M2) = PDD0(M1,M2)
            RDM0(M1+NN,M2+NN) = PDD0(M1,M2)
            RDM0(M1,M2+NN) = PDZ0(M1,M2)
            RDM0(M1+NN,M2) = PDZ0(M1,M2)
         END DO
      END DO
C-----------------PRINT----------
C	print*,'Re Pdm0:'
C       write(6,'(8f10.5)')((Real(Rdm0(m1,m2)),m2=1,n2),m1=1,n2)
C	print*,'Im Pdm0:'
C       write(6,'(8f10.5)')((Imag(Rdm0(m1,m2)),m2=1,n2),m1=1,n2)
C----------------------------------------
C
      DO M1 = 1,N2
         D(M1,M1) = DCMPLX(1.D0,0.D0)
         DO M2 = 1,N2
            IF ( M2.NE.M1 ) D(M1,M2) = DCMPLX(0.D0,0.D0)
         END DO
      END DO
C------- BLAS multiplication
      CALL ZGEMM('N','N',N2,N2,N2,DCMPLX(1.D0,0.D0),RDM0,N2,UDM,N2,
     &           DCMPLX(1.D0,0.D0),D,N2)
C------- LAPACK inversion of general complex matrix
      CALL ZGETRF(N2,N2,D,N2,IPIV,INFO)
      IF ( INFO.NE.0 ) STOP 'in RPA ZGETRF Info NE 0'
C
      CALL ZGETRI(N2,D,N2,IPIV,WORK,N2,INFO)
      IF ( INFO.NE.0 ) STOP 'in RPA ZGETRI Info NE 0'
C------- BLAS multiplication
      CALL ZGEMM('N','N',N2,N2,N2,DCMPLX(1.D0,0.D0),D,N2,RDM0,N2,
     &           DCMPLX(0.D0,0.D0),RDM,N2)
C-----------------PRINT----------
C	print*,'Re P:'
C       write(6,'(8f10.5)')((Real(Rdm(m1,m2)),m2=1,n2),m1=1,n2)
C	print*,'Im P:'
C       write(6,'(8f10.5)')((Imag(Rdm(m1,m2)),m2=1,n2),m1=1,n2)
C----------------------------------------
C--------- Calculate longitudinal potential W=0.5Umd*(Pmd-Po)*Umd
      DO M1 = 1,N2
         DO M2 = 1,N2
            RDM(M1,M2) = RDM(M1,M2) - RDM0(M1,M2)
         END DO
      END DO
C------- BLAS multiplications (factor 0.5 - later!)
      CALL ZGEMM('N','N',N2,N2,N2,DCMPLX(1.D0,0.D0),UDM,N2,RDM,N2,
     &           DCMPLX(0.D0,0.D0),D,N2)
      CALL ZGEMM('N','N',N2,N2,N2,DCMPLX(1.D0,0.D0),D,N2,UDM,N2,
     &           DCMPLX(0.D0,0.D0),RDM,N2)
C---------- Tranfer back to dd and zz
      DO M1 = 1,NN
         DO M2 = 1,NN
            VDD(M1,M2) = RDM(M1,M2)
            VZZ(M1,M2) = RDM(M1+NN,M2+NN)
            VDZ(M1,M2) = RDM(M1,M2+NN)
            VZD(M1,M2) = RDM(M1+NN,M2)
         END DO
      END DO
C
      END
C*==dmft_potsp.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C***********************************************************
      SUBROUTINE DMFT_POTSP(VDD,VZZ,VDZ,VZD,VDM,NN,NS)
      IMPLICIT NONE
C*--DMFT_POTSP760
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NN,NS
      COMPLEX*16 VDD(NN,NN),VDM(NN,NN,NS),VDZ(NN,NN),VZD(NN,NN),
     &           VZZ(NN,NN)
C
C Local variables
C
      INTEGER M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C---------- FLEX effective SP-potential for couple D and M channels -------
C------ Total pot in D-Mz channels
      DO M1 = 1,NN
         DO M2 = 1,NN
            VDM(M1,M2,1) = 0.5D0*(VDD(M1,M2)+VZZ(M1,M2)+VDZ(M1,M2)
     &                     +VZD(M1,M2))
            VDM(M1,M2,2) = 0.5D0*(VDD(M1,M2)+VZZ(M1,M2)-VDZ(M1,M2)
     &                     -VZD(M1,M2))
         END DO
      END DO
C
      END
C*==dmft_rpa.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C***********************************************************
      SUBROUTINE DMFT_RPA(U,R0,R,NN)
      IMPLICIT NONE
C*--DMFT_RPA807
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NN
      COMPLEX*16 R(NN,NN),R0(NN,NN),U(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(:,:),WORK(:,:)
      INTEGER INFO,IPIV(:,:),M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE D,IPIV,WORK
      ALLOCATE (D(NN,NN),IPIV(NN,NN),WORK(NN,NN))
C
C*** End of declarations rewritten by SPAG
C
C----------  FLEX "RPA" part: R=[1+R0*U]^-1*R0 --------
      DO M1 = 1,NN
         D(M1,M1) = DCMPLX(1.D0,0.D0)
         DO M2 = 1,NN
            IF ( M2.NE.M1 ) D(M1,M2) = DCMPLX(0.D0,0.D0)
         END DO
      END DO
C-----------------PRINT----------
C       print*,'Re P0+-:'
C       write(6,'(4f10.5)')((dReal(R0(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Im P0+-:'
C       write(6,'(4f10.5)')((dImag(R0(m1,m2)),m2=1,nn),m1=1,nn)
C------- BLAS multiplication
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),R0,NN,U,NN,
     &           DCMPLX(1.D0,0.D0),D,NN)
C------- LAPACK inversion of general complex matrix
      CALL ZGETRF(NN,NN,D,NN,IPIV,INFO)
      IF ( INFO.NE.0 ) STOP 'in RPA ZGETRF Info NE 0'
C
      CALL ZGETRI(NN,D,NN,IPIV,WORK,NN,INFO)
      IF ( INFO.NE.0 ) STOP 'in RPA ZGETRI Info NE 0'
C
C------- BLAS multiplication
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),D,NN,R0,NN,
     &           DCMPLX(0.D0,0.D0),R,NN)
C-----------------PRINT----------
C       print*,'Re P+-:'
C       write(6,'(4f10.5)')((dReal(R(m1,m2)),m2=1,nn),m1=1,nn)
C       print*,'Im P+-:'
C       write(6,'(4f10.5)')((dImag(R(m1,m2)),m2=1,nn),m1=1,nn)
      END
C*==dmft_pot.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C***********************************************************
      SUBROUTINE DMFT_POT(U,R0,R,V,NN)
      IMPLICIT NONE
C*--DMFT_POT876
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NN
      COMPLEX*16 R(NN,NN),R0(NN,NN),U(NN,NN),V(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(NN,NN)
      INTEGER M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
C
C      COMPLEX(KIND=8), DIMENSION(nn,nn) :: R0, R, D, U, V
C
C---------- FLEX effective chanal potential: Vd=Ud(R-R0)Ud --------
      DO M1 = 1,NN
         DO M2 = 1,NN
            R(M1,M2) = R(M1,M2) - R0(M1,M2)
         END DO
      END DO
C------- BLAS multiplication
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),R,NN,U,NN,
     &           DCMPLX(0.D0,0.D0),D,NN)
C------- BLAS multiplication
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),U,NN,D,NN,
     &           DCMPLX(0.D0,0.D0),V,NN)
C
      END
C*==dmft_tmaom.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C***********************************************************
      SUBROUTINE DMFT_TMAOM(T,GT,SIGT,MINOM,NLM,NS,NOM)
      IMPLICIT NONE
C*--DMFT_TMAOM926
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLM,NOM,NS
      COMPLEX*16 GT(NLM,NLM,NOM,NS),SIGT(NLM,NLM,NOM,NS),
     &           T(NLM,NLM,NLM,NLM,NOM,NS+1)
      INTEGER MINOM(NOM)
C
C Local variables
C
      INTEGER I,IS,M1,M2,M3,M4,M5,M6,M7,M8,MI
      COMPLEX*16 SS
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C---------- Sigma-2:    --------
C-------------- Spin-polarized -FLEX -------- Ns=2
      DO IS = 1,NS
         DO I = 1,NOM
C       do i=1,nom/2
            MI = MINOM(I)
            DO M1 = 1,NLM
               DO M2 = 1,NLM
                  SS = (0.D0,0.D0)
                  DO M3 = 1,NLM
                     DO M4 = 1,NLM
                        DO M5 = 1,NLM
                           DO M6 = 1,NLM
                              DO M7 = 1,NLM
                                 DO M8 = 1,NLM
C----------------- "Hartree" part
                                    SS = SS + T(M1,M5,M3,M6,1,IS)
     &                                 *T(M4,M7,M2,M8,1,IS)
     &                                 *GT(M3,M4,I,IS)*GT(M8,M5,MI,IS)
     &                                 *GT(M6,M7,I,IS)
     &                                 + T(M1,M5,M3,M6,1,3)
     &                                 *T(M4,M7,M2,M8,1,3)
     &                                 *GT(M3,M4,I,IS)*GT(M8,M5,MI,3-IS)
     &                                 *GT(M6,M7,I,3-IS)
C----------------- "Fock" part
                                    SS = SS - T(M1,M5,M3,M6,1,IS)
     &                                 *T(M4,M7,M8,M2,1,IS)
     &                                 *GT(M3,M4,I,IS)*GT(M8,M5,MI,IS)
     &                                 *GT(M6,M7,I,IS)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
C-------------- Total minus sign for Sig-2: -U*U*G*G*G
                  SIGT(M1,M2,I,IS) = SIGT(M1,M2,I,IS) - SS
               END DO
            END DO
         END DO
      END DO
C
      END
C***********************************************************
C
