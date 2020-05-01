C*==dmft_sig_sosptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C***********************************************************
      SUBROUTINE DMFT_SIG_SOSPTFLEX(G,SIG,EGM,TEMP,DC,SIGSTAT)
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NLM,NLMS,NS,NOM,OMEGA,DOMEGA,DTIME,
     &    MINOM,SIGNATURE,ENDSIGMA,TLAST,G_FIL,G_MU,PPP,PPH,UC,UW0,IPRT
      IMPLICIT NONE
C*--DMFT_SIG_SOSPTFLEX9
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER SCREEN
      PARAMETER (SCREEN=1)
C
C Dummy arguments
C
      INTEGER DC
      REAL*8 TEMP
      REAL*8 EGM(3),SIGSTAT(NLMS,NLMS)
      COMPLEX*16 G(NLMS,NLMS,NOM),SIG(NLMS,NLMS,NOM)
C
C Local variables
C
      REAL*8 AVS(:),EGMW,FI(:),RES,SIGAV,SIGDC(:,:),SIGDC2(:,:),
     &       SIGHF(:,:),TAU,UNIT(:,:)
      COMPLEX*16 D(:,:),G1(:,:,:),R0(:,:,:,:),SIGT(:,:,:),SS,WORK(:)
      REAL*8 G_MTAU,G_TAU
      INTEGER I,INFO,IPIV(:),IS,J,M,M1,M2,M3,M4,MI
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C
      ALLOCATABLE IPIV,WORK,UNIT,SIGHF,SIGDC,FI,D,R0,G1,SIGT,AVS
      ALLOCATABLE SIGDC2
C
C
      ALLOCATE (IPIV(NLMS*NLMS),WORK(NLMS))
      ALLOCATE (UNIT(NLMS,NLMS),SIGHF(NLMS,NLMS),SIGDC(NLMS,NLMS))
      ALLOCATE (SIGDC2(NLMS,NLMS))
      ALLOCATE (FI(NOM+1))
      ALLOCATE (D(NLMS,NLMS),R0(NLMS,NLMS,NLMS,NLMS))
      ALLOCATE (G1(NLMS,NLMS,NOM),SIGT(NLMS,NLMS,NOM))
      ALLOCATE (AVS(2))
      IPIV = 0
      WORK = (0.0D0,0.0D0)
      D = (0.0D0,0.0D0)
      UNIT = 0.0D0
      SIGHF = 0.0D0
      SIGDC = 0.0D0
      FI = 0.0D0
      R0 = (0.0D0,0.0D0)
      G1 = (0.0D0,0.0D0)
      SIGT = (0.0D0,0.0D0)
      AVS = 0.0D0
      SIGAV = 0.0D0
C
      ALLOCATE (ENDSIGMA(NLMS,NLMS))
      ALLOCATE (G_MU(NLMS,NLMS))
      ALLOCATE (G_FIL(NLMS,NLMS))
      ALLOCATE (TLAST(NLMS,NLMS,NLMS,NLMS))
C
C
      ENDSIGMA = (0D0,0.0D0)
      G_MU(1:NLMS,1:NLMS) = 0.D0
      DO M1 = 1,NLMS
         UNIT(M1,M1) = 1D0
      END DO
C
C     write(*,*) '---- write green function ----- '
      IF ( IPRT.GT.0 ) THEN
C
         WRITE (6,*) ' DC = ',DC
         DO I = 1,NOM/2
            WRITE (44,'(f12.6,16(2f12.6))') OMEGA(I) + DOMEGA/2D0,
     &             (G(M,M,2*I),M=1,NLMS)
         END DO
C
      END IF
C
C     Save G in omega
      G1(1:NLMS,1:NLMS,1:NOM) = G(1:NLMS,1:NLMS,1:NOM)
C   obtain g(t)
      SIGNATURE = 'G'
      CALL DMFT_FFT_FTOT2_SO(G,TEMP,NLMS,NOM)
C
      ALLOCATE (PPP(NLMS,NLMS,NLMS,NLMS,NOM))
C
C---Particle particle suseptebility in tau
      DO I = 1,NOM
         MI = MINOM(I)
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               DO M3 = 1,NLMS
                  DO M4 = 1,NLMS
C------------- density d
C-------------G(-tau)=-G(-tau+beta)
C pp-chanel Ppp(13,24)
                     PPP(M1,M2,M3,M4,I) = G(M1,M3,I)*G(M2,M4,I)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C Substract analytical part with discontinuty
C
      DO I = 1,NOM
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               PPP(M1,M2,M1,M2,I) = PPP(M1,M2,M1,M2,I)
     &                              - G_TAU(M1,M1,I,TEMP)
     &                              *G_TAU(M2,M2,I,TEMP)
            END DO
         END DO
      END DO
C
C------- Susceptibility Ppp in omega
      CALL DMFT_FFT_TTOF4_SO(PPP,TEMP,NLMS,NOM)
C
C Add analytically transformed part
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
C
            CALL DMFT_PPP_AN(M1,M2,M1,M2,1,TEMP)
C
         END DO
      END DO
C
      IF ( IPRT.GT.0 ) WRITE (6,'(a,14f14.6)') 'Hi(pp)(0):',
     &                        (REAL(PPP(M,M,M,M,1)),M=1,NLMS)
C
C
C----------  Screening  U ------------
      CALL DMFT_TMATR_SO(UC,PPP)
C
C Set the vertex for the PH channel
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
                  IF ( SCREEN.EQ.1 ) THEN
                     UW0(M1,M2,M3,M4) = PPP(M1,M3,M2,M4,1)
     &                                  - PPP(M1,M3,M4,M2,1)
                  ELSE
                     UW0(M1,M2,M3,M4) = UC(M1,M3,M2,M4)
     &                                  - UC(M1,M3,M4,M2)
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF ( IPRT.GT.0 ) WRITE (6,'(a,14f14.6)') 'W(pp)(0):',
     &                        (REAL(PPP(M,M,M,M,1)),M=1,NLMS)
C
C  Remove the first order from TMAT
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
                  DO I = 1,NOM,2
                     PPP(M1,M2,M3,M4,I) = PPP(M1,M2,M3,M4,I)
     &                  - UC(M1,M2,M3,M4)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C Remove the second order (analytical part) from Tmat
C
      DO I = 1,NOM,2
         R0 = (0D0,0D0)
         SS = (0.0D0,0.0D0)
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               CALL DMFT_PPP_AN_I(SS,I,M1,M2,M1,M2,0,TEMP)
               R0(M1,M2,M1,M2) = SS
            END DO
         END DO
         CALL DMFT_TMAT2(R0,UC)
         PPP(:,:,:,:,I) = PPP(:,:,:,:,I) + R0
      END DO
C
C
      CALL DMFT_FFT_FTOT4_SO(PPP,TEMP,NLMS,NOM)
C
C Add analytical part of Tmat in tau
C
C Last point tau=beta
      R0 = (0D0,0D0)
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            R0(M1,M2,M1,M2) = R0(M1,M2,M1,M2)
     &                        - (-UNIT(M1,M1)-G_TAU(M1,M1,1,TEMP))
     &                        *(-UNIT(M2,M2)-G_TAU(M2,M2,1,TEMP))
         END DO
      END DO
      CALL DMFT_TMAT2(R0,UC)
      TLAST(:,:,:,:) = PPP(:,:,:,:,1) + R0
C
C All other points
C
      DO I = 1,NOM
         R0 = (0D0,0D0)
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               R0(M1,M2,M1,M2) = R0(M1,M2,M1,M2) - G_TAU(M1,M1,I,TEMP)
     &                           *G_TAU(M2,M2,I,TEMP)
            END DO
         END DO
         CALL DMFT_TMAT2(R0,UC)
         PPP(:,:,:,:,I) = PPP(:,:,:,:,I) + R0
      END DO
C      do m1=1,nlms
C      do m2=1,nlms
C      do m3=1,nlms
C      do m4=1,nlms
C       IF(ABS(Ppp(m1,m2,m3,m4,1)) > 1d-6)
C     & write(25,'(4i6,2f15.8)')m1,m2,m3,m4,Ppp(m1,m2,m3,m4,1)
C      ENDDO
C      ENDDO
C      ENDDO
C      ENDDO
C      write(25,'(i6,f12.6,2f15.8)')i,dtime*nom,Tlast(1,1,1,1)
C
C
C Calculate the Tmat-contribution to Sig
C
C tau=+0
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            SS = (0.D0,0.D0)
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
                  SS = SS - PPP(M1,M3,M2,M4,1)*(-UNIT(M4,M3)-G(M4,M3,1))
     &                 + PPP(M1,M4,M3,M2,1)*(-UNIT(M3,M4)-G(M3,M4,1))
               END DO
            END DO
            SIG(M1,M2,1) = SS
         END DO
      END DO
C All other points
      DO I = 2,NOM
         MI = MINOM(I)
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               SS = (0.D0,0.D0)
               DO M3 = 1,NLMS
                  DO M4 = 1,NLMS
                     SS = SS - PPP(M1,M3,M2,M4,I)*G(M4,M3,MI)
     &                    + PPP(M1,M4,M3,M2,I)*G(M3,M4,MI)
                  END DO
               END DO
               SIG(M1,M2,I) = SS
            END DO
         END DO
      END DO
C
C
C Add contribution from T-matrix to endsigma
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            SS = (0.D0,0.D0)
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
                  SS = SS - TLAST(M1,M3,M2,M4)*G(M4,M3,1)
     &                 + TLAST(M1,M4,M3,M2)*G(M3,M4,1)
               END DO
            END DO
            ENDSIGMA(M1,M2) = SS
         END DO
      END DO
C      do i=1,nom
C       write(81,'(i6,2f14.7)')i,Sig(1,1,i)
C      enddo
C      write(81,'(i6,2f14.7)')nom+1,endsigma(1,1)
C
      DEALLOCATE (PPP)
C
C Calculate  the particle hole suseptebility
C
      ALLOCATE (PPH(NLMS,NLMS,NLMS,NLMS,NOM))
      DO I = 2,NOM
         MI = MINOM(I)
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               DO M3 = 1,NLMS
                  DO M4 = 1,NLMS
C------------- density d
C-------------G(-tau)=-G(-tau+beta)
                     PPH(M1,M2,M3,M4,I) = G(M4,M1,MI)*G(M2,M3,I)
                  END DO
               END DO
            END DO
         END DO
      END DO
C tau=0+
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
C-------------        -G(-0)=-1-G(+0)
                  PPH(M1,M2,M3,M4,1) = -(UNIT(M4,M1)+G(M4,M1,1))
     &                                 *G(M2,M3,1)
               END DO
            END DO
         END DO
      END DO
C
C Substract analytical part with discontinuity
C
      DO I = 1,NOM
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               PPH(M1,M2,M2,M1,I) = PPH(M1,M2,M2,M1,I)
     &                              + G_MTAU(M1,M1,I,TEMP)
     &                              *G_TAU(M2,M2,I,TEMP)
            END DO
         END DO
      END DO
C
      CALL DMFT_FFT_TTOF4_SO(PPH,TEMP,NLMS,NOM)
C
C Add analytically transformed part
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
C
            CALL DMFT_PPH_AN(M1,M2,M2,M1,1,TEMP)
C
         END DO
      END DO
C
      IF ( IPRT.GT.0 ) WRITE (6,'(a,28f14.6)') 'Hi1(ph)(0):',
     &                        (PPH(M,M,M,M,1),M=1,NLMS)
C
C---------- W   -------
      CALL DMFT_TOTPOT_SO(PPH,UW0)
C
      IF ( IPRT.GT.0 ) WRITE (6,'(a,14f14.6)') 'W(ph)(0):',
     &                        (REAL(PPH(M,M,M,M,1)),M=1,NLMS)
C W in time
      CALL DMFT_FFT_FTOT4_SO(PPH,TEMP,NLMS,NOM)
C
      DO I = 1,NOM
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               SS = (0.D0,0.D0)
               DO M3 = 1,NLMS
                  DO M4 = 1,NLMS
                     SS = SS + PPH(M1,M3,M4,M2,I)*G(M3,M4,I)
                  END DO
               END DO
               SIG(M1,M2,I) = SIG(M1,M2,I) + SS
            END DO
         END DO
      END DO
C
C Last point t=beta
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            SS = (0.D0,0.D0)
            DO M3 = 1,NLMS
               DO M4 = 1,NLMS
                  SS = SS + PPH(M1,M3,M4,M2,1)*(-UNIT(M3,M4)-G(M3,M4,1))
               END DO
            END DO
            ENDSIGMA(M1,M2) = ENDSIGMA(M1,M2) + SS
         END DO
      END DO
C      do i=1,nom
C       write(82,'(i6,2f14.7)')i,Sig(1,1,i)
C      enddo
C      write(82,'(i6,2f14.7)')nom+1,endsigma(1,1)
C
      DEALLOCATE (PPH)
C
C     Save Sigma in tau space
C
      SIGT(1:NLMS,1:NLMS,1:NOM) = SIG(1:NLMS,1:NLMS,1:NOM)
C
C Transform Sig to omega space
      SIGNATURE = 'S'
      CALL DMFT_FFT_TTOF2_SO(SIG,TEMP,NLMS,NOM)
C
C Add real constant contribution from U*delta(tau)
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            DO M3 = 1,NLMS
               DO IS = 1,2
                  DO M4 = 1 + NLM*(IS-1),NLM*IS
                     SIGHF(M1,M2) = SIGHF(M1,M2)
     &                              + DREAL(UC(M1,M3,M2,M4))
     &                              *(UNIT(M4,M3)+DREAL(G(M4,M3,1)))
     &                              - DREAL(UC(M1,M4,M3,M2))
     &                              *(UNIT(M3,M4)+DREAL(G(M3,M4,1)))
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      IF ( IPRT.GT.0 ) CALL RMATSTRUCT('SIGHF',SIGHF,NLMS,NLMS,0,0,0,
     &                                 1D-8,6)
C
C
C
      IF ( DC.EQ.3 ) THEN
C Use LDA+U-AMF DC
         AVS = 0D0
         DO IS = 1,NS
            DO M1 = (IS-1)*NLM + 1,IS*NLM
               AVS(IS) = AVS(IS) + (1D0+DREAL(G(M1,M1,1)))
            END DO
         END DO
C
         AVS = AVS/NLM
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               DO IS = 1,NS
                  DO M3 = (IS-1)*NLM + 1,IS*NLM
                     SIGHF(M1,M2) = SIGHF(M1,M2)
     &                              - DREAL(UC(M1,M3,M2,M3)-
     &                              UC(M1,M3,M3,M2))*AVS(IS)
     &                              *(1D0-1D0/NLM*UNIT(M1,M2))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      IF ( DC.NE.2 ) THEN
C Add HF contribution
         DO I = 2,NOM,2
            SIG(:,:,I) = SIG(:,:,I) + SIGHF
         END DO
      END IF
C
C Print out the self-energy
C
      IF ( IPRT.GT.0 ) THEN
C
         DO I = 1,NOM/2
            WRITE (45,'(f12.6,16(2f12.6))') OMEGA(I) + DOMEGA/2D0,
     &             (SIG(M,M,2*I),M=1,NLMS)
         END DO
      END IF
C
C Substract the double conting
C
      IF ( DC.EQ.2 .OR. DC.EQ.0 ) THEN
         SIGDC = 0
      ELSE IF ( DC.EQ.1 ) THEN
         SIGDC = REAL(SIG(:,:,2),KIND=8)
      ELSE IF ( DC.EQ.4 ) THEN
         SIGDC = 0.0D0
         DO IS = 1,2
            SIGAV = 0.0D0
            DO I = 1,NLM
               J = I + NLM*(IS-1)
               SIGAV = SIGAV + REAL(SIG(J,J,2),KIND=8)
            END DO
            DO I = 1,NLM
               J = I + NLM*(IS-1)
               SIGDC(J,J) = SIGAV/DFLOAT(NLM)
            END DO
         END DO
      ELSE IF ( DC.EQ.5 ) THEN
         SIGDC = 0.0D0
         DO IS = 1,2
            SIGAV = 0.0D0
            DO I = 1,NLM
               J = I + NLM*(IS-1)
               SIGAV = SIGAV + REAL(SIG(J,J,2),KIND=8)
            END DO
            DO I = 1,NLM
               J = I + NLM*(IS-1)
               SIGDC(J,J) = SIGAV/DFLOAT(NLM)
     &                      - 2.0D0*(REAL(SIG(J,J,2),KIND=8)
     &                      -SIGAV/DFLOAT(NLM))
            END DO
         END DO
      END IF
C
C
C************************** IGOR *****************************
C here we should write sigma top
C we change the sign of sigma top since SIGDC is substracted
C and not added to the total self energy
C notice that you can define SIGMATOP both real or complex
      SIGDC(:,:) = SIGDC(:,:) - SIGSTAT(:,:)
C
C re-initialize sigma HF
      SIGHF = 0D0
C************************ IGOR ********************************
C
C
C
      DO I = 2,NOM,2
         SIG(:,:,I) = SIG(:,:,I) - SIGDC(:,:)
      END DO
C
C Print out the self-energy
C
      IF ( IPRT.GT.0 ) THEN
C
         DO I = 1,NOM/2
            WRITE (50,'(f12.6,16(2f12.6))') OMEGA(I) + DOMEGA/2D0,
     &             (SIG(M,M,2*I),M=1,NLMS)
         END DO
      END IF
C
C Get the new local G
C
      DO I = 2,NOM,2
         D = G1(:,:,I)
C------- LAPACK inversion of GF
         CALL ZGETRF(NLMS,NLMS,D,NLMS,IPIV,INFO)
         IF ( INFO.NE.0 ) STOP 'in SPTFLEX ZGETRF Info NE 0'
C
         CALL ZGETRI(NLMS,D,NLMS,IPIV,WORK,NLMS,INFO)
         IF ( INFO.NE.0 ) STOP 'in SPTFLEX ZGETRI Info NE 0'
         D = D - SIG(:,:,I)
C
C------- LAPACK inversion of (GF^-1-Sig)
         CALL ZGETRF(NLMS,NLMS,D,NLMS,IPIV,INFO)
         IF ( INFO.NE.0 ) STOP 'in SPTFLEX ZGETRF Info NE 0'
C
         CALL ZGETRI(NLMS,D,NLMS,IPIV,WORK,NLMS,INFO)
         IF ( INFO.NE.0 ) STOP 'in SPTFLEX ZGETRI Info NE 0'
         G1(:,:,I) = D
      END DO
C
C Save final local G in omega
C
      IF ( IPRT.GT.0 ) THEN
C
         DO I = 1,NOM/2
            WRITE (51,'(f12.6,16(2f12.6))') OMEGA(I) + DOMEGA/2D0,
     &             (G1(M,M,2*I),M=1,NLMS)
         END DO
      END IF
      G = G1
C
C Get local G in tau
C
      SIGNATURE = 'G'
      CALL DMFT_FFT_FTOT2_SO(G1,TEMP,NLMS,NOM)
C
C***********************************************************
C
C Calculate the first contribution to the GM energy
C E=Tr[G(tau)Sig(-tau)]
C
      EGM(1:3) = 0.0D0
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            DO I = 2,NOM
               MI = MINOM(I)
               FI(I) = -DREAL(G1(M1,M2,I)*SIGT(M2,M1,MI))
            END DO
            FI(1) = -DREAL(G1(M1,M2,1)*ENDSIGMA(M2,M1))
            FI(NOM+1) = DREAL(-(-UNIT(M1,M2)-G1(M1,M2,1))*SIGT(M2,M1,1))
            CALL DMFT_SIMPN(FI,DTIME,NOM+1,RES)
            EGM(1) = EGM(1) + RES/2D0
         END DO
      END DO
C
      IF ( IPRT.GT.0 ) PRINT *,'Dynamic contribution',EGM(1)
C
C Calculate the second contribution to GM energy
C   Egm(2)=T/2*Tr[(G(tau=0+)*(SigHF-SigDC)]
C
      IF ( IPRT.GT.0 ) CALL RMATSTRUCT('SIGSTAT',SIGSTAT,NLMS,NLMS,0,0,
     &                                 0,1D-8,6)
C
      DO M1 = 1,NLMS
         DO M2 = 1,NLMS
            EGM(2) = EGM(2) + DREAL(UNIT(M1,M2)+G1(M1,M2,1))
     &               *(SIGSTAT(M2,M1))/2D0
         END DO
      END DO
C
C
C
      EGM(3) = EGM(1) + EGM(2)
C
      IF ( IPRT.GT.0 ) THEN
         PRINT *,'Static contribution',EGM(2)
         PRINT *,'Total',EGM(3)
      END IF
C
C
C*****************************************************
C
C
      IF ( IPRT.GT.0 ) THEN
C Calculate the GM contribution as a sum over the matsubara frequencies
         PRINT *,'Start Matsubara sum for GM contribution to tot E'
         EGMW = 0D0
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               SS = 0D0
               DO I = 1,NOM/2
                  SS = SS + G(M1,M2,2*I)*SIG(M2,M1,2*I)
               END DO
               EGMW = EGMW + TEMP*REAL(SS)
            END DO
         END DO
         PRINT *,'Final Value of the Correction:'
         PRINT *,EGMW
C
C
         PRINT *,'Asymptotics Corrections'
C
         EGMW = 0D0
C
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               IF ( ABS(G(M1,M2,NOM)).GT.1.D-20 ) THEN
                  G_MU(M1,M2) = -REAL(1D0/G(M1,M2,NOM))
               ELSE
                  G_MU(M1,M2) = 0.0D0
               END IF
            END DO
         END DO
         TAU = DOMEGA*(NOM/2) - DOMEGA/2D0
         SIGDC2(:,:) = REAL(SIG(:,:,NOM)) + REAL(G(:,:,NOM))
     &                 /(1D0+(REAL(G(:,:,NOM))**2)*(TAU**2))
C
         DO M1 = 1,NLMS
            DO M2 = 1,NLMS
               SS = 0D0
               IF ( ABS(G(M1,M2,NOM)).GE.1.D-5 ) THEN
                  DO I = 1,NOM/2
                     TAU = DOMEGA*I - DOMEGA/2D0
                     SS = SS + (G_MU(M1,M2)/((G_MU(M1,M2)**2)+(TAU**2)))
                  END DO
                  IF ( G_MU(M1,M2)/TEMP.LT.-70.00D0 ) THEN
                     SS = TEMP*SS
                  ELSE
                     SS = TEMP*SS - 0.5D0/(1+EXP(-G_MU(M1,M2)/TEMP))
                  END IF
               ELSE IF ( REAL(G(M1,M2,NOM)).GE.0 ) THEN
                  SS = 0D0
               ELSE
                  SS = -0.5D0
               END IF
                  !order of magnitude
C
               EGMW = EGMW + REAL(SS)*SIGDC2(M2,M1)
            END DO
         END DO
C
         PRINT *,'Final Value of the Correction:'
         PRINT *,EGMW
      END IF
C
C
      DEALLOCATE (ENDSIGMA,TLAST,G_MU,G_FIL)
      END
C*==dmft_tmatr_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
C
C
C
      SUBROUTINE DMFT_TMATR_SO(UC,PPP)
C***********************************************************
C
C    Calculates the T matrix in spirit of Kanamori Approach
C    U= U*[1+R0*U]^{-1}
C
C***********************************************************
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NN,NOM,IPRT
      IMPLICIT NONE
C*--DMFT_TMATR_SO703
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 PPP(NN,NN,NOM),UC(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(NN,NN),R0(NN,NN),UR(NN,NN),WORK(NN)
      INTEGER I,INFO,IPIV(NN,NN),M,M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C----------- put 0 for Fermion Matsubara
      DO I = 2,NOM,2
         PPP(:,:,I) = (0D0,0D0)
      END DO
C---------- Tmat
      DO I = 1,NOM,2
         DO M1 = 1,NN
            DO M2 = 1,NN
               R0(M1,M2) = PPP(M1,M2,I)
            END DO
         END DO
         IF ( I.EQ.1 .AND. IPRT.GT.0 ) WRITE (6,'(5f14.5)')
     &        (REAL(R0(M,M)),M=1,NN)
C----------  PP-screening: Ur=U*[1+R0*U]^-1 --------
         DO M1 = 1,NN
            D(M1,M1) = DCMPLX(1.D0,0.D0)
            DO M2 = 1,NN
               IF ( M2.NE.M1 ) D(M1,M2) = DCMPLX(0.D0,0.D0)
            END DO
         END DO
C------- BLAS multiplication
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),UC,NN,R0,NN,
     &              DCMPLX(1.D0,0.D0),D,NN)
C------- LAPACK inversion of general complex matrix
         CALL ZGETRF(NN,NN,D,NN,IPIV,INFO)
         IF ( INFO.NE.0 ) STOP 'in RPA ZGETRF Info NE 0'
C
         CALL ZGETRI(NN,D,NN,IPIV,WORK,NN,INFO)
         IF ( INFO.NE.0 ) STOP 'in RPA ZGETRI Info NE 0'
C
C------- BLAS multiplication
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),D,NN,UC,NN,
     &              DCMPLX(0.D0,0.D0),UR,NN)
C
         DO M1 = 1,NN
            DO M2 = 1,NN
               PPP(M1,M2,I) = UR(M1,M2)
            END DO
         END DO
C
      END DO ! iom
      IF ( IPRT.GT.0 ) THEN
C-----------------PRINT----------
         PRINT *,'Tmat(ii)='
         WRITE (6,'(10f8.3)') (DREAL(PPP(M,M,1)),M=1,NN)
C
      END IF
      END
C*==dmft_pph_an.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
      SUBROUTINE DMFT_PPH_AN(M1,M2,M3,M4,ADD,TEMP)
C
C Add to (add=1) PH-susept. analytical part
C or substitute PH-susept. with (add=0) analytical part
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:OMEGA,NOM,CI,PPH,G_FIL,G_MU
      IMPLICIT NONE
C*--DMFT_PPH_AN790
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ADD,M1,M2,M3,M4
      REAL*8 TEMP
C
C Local variables
C
      REAL*8 BETA,MU,MU1,MU2,TMP1
      INTEGER I,JSIGN
      COMPLEX*16 I_W
C
C*** End of declarations rewritten by SPAG
C
C
C     IF(add/=0) add=1
C Set values to 0 if required
      DO I = 1,NOM/2
         PPH(M1,M2,M3,M4,2*I-1) = ADD*PPH(M1,M2,M3,M4,2*I-1)
      END DO
C
      BETA = 1D0/TEMP
C
      MU1 = G_MU(M4,M1)
      MU2 = G_MU(M2,M3)
      IF ( (G_FIL(M4,M1).EQ.1 .AND. G_FIL(M2,M3).EQ.0) .OR. 
     &     (G_FIL(M4,M1).EQ.0 .AND. G_FIL(M2,M3).EQ.1) ) THEN
C One band is completely filled and another one is completely empty
         DO I = 1,NOM/2
            PPH(M1,M2,M3,M4,2*I-1) = PPH(M1,M2,M3,M4,2*I-1) + 1D0
         END DO
C One band is completely filled or empty
      ELSE IF ( G_FIL(M4,M1).NE.2 .AND. G_FIL(M2,M3).EQ.2 .OR. 
     &          G_FIL(M4,M1).EQ.2 .AND. G_FIL(M2,M3).NE.2 ) THEN
         IF ( G_FIL(M2,M3).EQ.2 ) THEN
            MU = G_MU(M2,M3)
            JSIGN = 2*G_FIL(M4,M1) - 1
         ELSE
            MU = G_MU(M4,M1)
            JSIGN = 2*G_FIL(M2,M3) - 1
         END IF
         DO I = 1,NOM/2
            PPH(M1,M2,M3,M4,2*I-1) = PPH(M1,M2,M3,M4,2*I-1)
     &                               + 1D0/(1D0+EXP(JSIGN*MU*BETA))
         END DO
C Some filling in the both bands
      ELSE IF ( G_FIL(M4,M1).EQ.2 .AND. G_FIL(M2,M3).EQ.2 ) THEN
         IF ( ABS(MU1-MU2).GT.1D-6 ) THEN
            TMP1 = (EXP(-MU2*BETA)-EXP(-MU1*BETA))
     &             /((EXP(-MU2*BETA)+1D0)*(EXP(-MU1*BETA)+1D0))
            DO I = 1,NOM/2
               I_W = CI*OMEGA(I)
               PPH(M1,M2,M3,M4,2*I-1) = PPH(M1,M2,M3,M4,2*I-1)
     &                                  + TMP1/(I_W+MU1-MU2)
            END DO
         ELSE
C Add correction only to tau=0 point
            PPH(M1,M2,M3,M4,1) = PPH(M1,M2,M3,M4,1)
     &                           + BETA*EXP(-MU2*BETA)
     &                           /((EXP(-MU2*BETA)+1D0)
     &                           *(EXP(-MU1*BETA)+1D0))
         END IF
      END IF
C
      END
C*==dmft_ppp_an.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
      SUBROUTINE DMFT_PPP_AN(M1,M2,M3,M4,ADD,TEMP)
C
C Add to (add=1) PP-susept. analytical part
C or substitute PP-susept. with (add=0) analytical part
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:OMEGA,NOM,CI,PPP,G_FIL,G_MU
      IMPLICIT NONE
C*--DMFT_PPP_AN882
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ADD,M1,M2,M3,M4
      REAL*8 TEMP
C
C Local variables
C
      REAL*8 BETA,MU,MU1,MU2,TMP
      INTEGER I,JSIGN
      COMPLEX*16 I_W
C
C*** End of declarations rewritten by SPAG
C
C
C Set values to 0 if required
C     IF(add/=0) add=1
      DO I = 1,NOM/2
         PPP(M1,M2,M3,M4,2*I-1) = ADD*PPP(M1,M2,M3,M4,2*I-1)
      END DO
C
      BETA = 1D0/TEMP
C
      IF ( G_FIL(M1,M3).EQ.1 .AND. G_FIL(M2,M4).EQ.1 .OR. G_FIL(M1,M3)
     &     .EQ.0 .AND. G_FIL(M2,M4).EQ.0 ) THEN
C Both bands are completely filled or empty
         DO I = 1,NOM/2
            PPP(M1,M2,M3,M4,2*I-1) = PPP(M1,M2,M3,M4,2*I-1) + 1D0
         END DO
C
      ELSE IF ( G_FIL(M1,M3).NE.2 .AND. G_FIL(M2,M4).EQ.2 .OR. 
     &          G_FIL(M1,M3).NE.2 .AND. G_FIL(M2,M4).EQ.2 ) THEN
C One band is completely filled or empty
         IF ( G_FIL(M1,M3).EQ.2 ) THEN
            MU = G_MU(M1,M3)
            JSIGN = 1 - 2*G_FIL(M2,M4)
         ELSE
            MU = G_MU(M2,M4)
            JSIGN = 1 - 2*G_FIL(M1,M3)
         END IF
         DO I = 1,NOM/2
            PPP(M1,M2,M3,M4,2*I-1) = PPP(M1,M2,M3,M4,2*I-1)
     &                               + 1D0/(1D0+EXP(JSIGN*MU*BETA))
         END DO
C
      ELSE IF ( G_FIL(M1,M3).EQ.2 .AND. G_FIL(M2,M4).EQ.2 ) THEN
C Some filling in the both bands
         MU1 = G_MU(M1,M3)
         MU2 = G_MU(M2,M4)
         IF ( ABS(MU1+MU2).GT.1D-6 ) THEN
C Not half-filled case
            TMP = (EXP(-(MU1+MU2)*BETA)-1D0)/(EXP(-MU1*BETA)+1D0)
     &            /(EXP(-MU2*BETA)+1D0)
            DO I = 1,NOM/2
               I_W = CI*OMEGA(I)
               PPP(M1,M2,M3,M4,2*I-1) = ADD*PPP(M1,M2,M3,M4,2*I-1)
     &                                  + TMP/(I_W-MU1-MU2)
            END DO
         ELSE
C Antisymmetric in mu case, correction only to the first Matsubara
            PPP(M1,M2,M3,M4,1) = ADD*PPP(M1,M2,M3,M4,1)
     &                           + BETA/(EXP(-MU1*BETA)+1D0)
     &                           /(EXP(-MU2*BETA)+1D0)
         END IF
C
      END IF
C
      END
C*==dmft_ppp_an_i.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_PPP_AN_I(RES,I,M1,M2,M3,M4,ADD,TEMP)
C
C Add to (add=1) PP-susept. analytical part
C or substitute PP-susept. with (add=0) analytical part
C
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:OMEGA,CI,G_FIL,G_MU
      IMPLICIT NONE
C*--DMFT_PPP_AN_I976
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ADD,I,M1,M2,M3,M4
      COMPLEX*16 RES
      REAL*8 TEMP
C
C Local variables
C
      REAL*8 BETA,MU,MU1,MU2,TMP
      COMPLEX*16 I_W
      INTEGER JSIGN,K
C
C*** End of declarations rewritten by SPAG
C
C
C Set values to 0 if required
C     IF(add/=0) add=1
      RES = ADD*RES
      K = (I+1)/2
C
      BETA = 1D0/TEMP
C
      IF ( G_FIL(M1,M3).EQ.1 .AND. G_FIL(M2,M4).EQ.1 .OR. G_FIL(M1,M3)
     &     .EQ.0 .AND. G_FIL(M2,M4).EQ.0 ) THEN
C Both bands are completely filled or empty
         RES = RES + 1D0
C
      ELSE IF ( G_FIL(M1,M3).NE.2 .AND. G_FIL(M2,M4).EQ.2 .OR. 
     &          G_FIL(M1,M3).NE.2 .AND. G_FIL(M2,M4).EQ.2 ) THEN
C One band is completely filled or empty
         IF ( G_FIL(M1,M3).EQ.2 ) THEN
            MU = G_MU(M1,M3)
            JSIGN = 1 - 2*G_FIL(M2,M4)
         ELSE
            MU = G_MU(M2,M4)
            JSIGN = 1 - 2*G_FIL(M1,M3)
         END IF
         RES = RES + 1D0/(1D0+EXP(JSIGN*MU*BETA))
C
      ELSE IF ( G_FIL(M1,M3).EQ.2 .AND. G_FIL(M2,M4).EQ.2 ) THEN
C Some filling in the both bands
         MU1 = G_MU(M1,M3)
         MU2 = G_MU(M2,M4)
         IF ( ABS(MU1+MU2).GT.1D-6 ) THEN
C Not half-filled case
            TMP = (EXP(-(MU1+MU2)*BETA)-1D0)/(EXP(-MU1*BETA)+1D0)
     &            /(EXP(-MU2*BETA)+1D0)
            I_W = CI*OMEGA(K)
            RES = RES + TMP/(I_W-MU1-MU2)
         ELSE
C Antisymmetric in mu case, correction only to the first Matsubara
            RES = RES + BETA/(EXP(-MU1*BETA)+1D0)/(EXP(-MU2*BETA)+1D0)
         END IF
C
      END IF
C
      END
C*==dmft_totpot_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C***********************************************************
      SUBROUTINE DMFT_TOTPOT_SO(PPH,UW0)
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NN,NOM
      IMPLICIT NONE
C*--DMFT_TOTPOT_SO1054
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 PPH(NN,NN,NOM),UW0(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(NN,NN),DD(NN,NN),RPP(NN,NN),RW(NN,NN),RW0(NN,NN),
     &           WORK(NN)
      INTEGER INFO,IOM,IPIV(NN,NN),M1,M2
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C----------  only Bosonic Matsubara --------
C Calculate W =  Uw0 * (1 - Pph * Uw0)^-1 * Pph *Uw0
      DO IOM = 1,NOM,2
C
         DO M1 = 1,NN
            D(M1,M1) = DCMPLX(1.D0,0.D0)
            DO M2 = 1,NN
               IF ( M2.NE.M1 ) D(M1,M2) = DCMPLX(0.D0,0.D0)
            END DO
         END DO
C
         DO M1 = 1,NN
            DO M2 = 1,NN
               RPP(M1,M2) = PPH(M1,M2,IOM)
            END DO
         END DO
C
C Perform RPA summation
C
C------- BLAS multiplication
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),RPP,NN,UW0,NN,
     &              DCMPLX(1.D0,0.D0),D,NN)
C------- LAPACK inversion of general complex matrix
         CALL ZGETRF(NN,NN,D,NN,IPIV,INFO)
         IF ( INFO.NE.0 ) STOP 'in RPA ZGETRF Info NE 0'
C
         CALL ZGETRI(NN,D,NN,IPIV,WORK,NN,INFO)
         IF ( INFO.NE.0 ) STOP 'in RPA ZGETRI Info NE 0'
C
C------- BLAS multiplication U*hi_0
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),D,NN,RPP,NN,
     &              DCMPLX(0.D0,0.D0),RW0,NN)
C
C  Substract hi_0
C
         DD = RW0 - RPP
C
C------- BLAS multiplication U*hi_0
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),UW0,NN,DD,NN,
     &              DCMPLX(0.D0,0.D0),D,NN)
C
         CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),D,NN,UW0,NN,
     &              DCMPLX(0.D0,0.D0),RW,NN)
C
         DO M1 = 1,NN
            DO M2 = 1,NN
               PPH(M1,M2,IOM) = RW(M1,M2)
            END DO
         END DO
C
      END DO
      END
C*==dmft_fft_ttof2_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C============================================================================
      SUBROUTINE DMFT_FFT_TTOF2_SO(FT,TEMP,NLMS,NOM)
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:SIGNATURE,ENDSIGMA,ENDSIG
      IMPLICIT NONE
C*--DMFT_FFT_TTOF2_SO1144
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER FERMION
      PARAMETER (FERMION=1)
C
C Dummy arguments
C
      INTEGER NLMS,NOM
      DOUBLE PRECISION TEMP
      COMPLEX*16 FT(NLMS,NLMS,NOM)
C
C Local variables
C
      INTEGER DIAG,IK,IL
      COMPLEX*16 F(NOM),OUT(0:NOM-1)
C
C*** End of declarations rewritten by SPAG
C
C*** Start of declarations inserted by SPAG
C*** End of declarations inserted by SPAG
C
      DO IL = 1,NLMS
         DO IK = 1,NLMS
            DIAG = 0
            IF ( IL.EQ.IK ) DIAG = 1
            IF ( SIGNATURE.EQ.'S' ) THEN
               DIAG = 1
               ENDSIG = ENDSIGMA(IK,IL)
            END IF
            F(1:NOM) = FT(IK,IL,1:NOM)
C      call dfour1(f,nom,-1)
C tau -> omega
            CALL DMFT_NFOURIER2(F,OUT,NOM,NOM-1,DIAG,FERMION,TEMP)
C
            FT(IK,IL,1:NOM) = OUT(0:NOM-1)
         END DO
      END DO
C
      END
C*==dmft_fft_ftot2_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C============================================================================
      SUBROUTINE DMFT_FFT_FTOT2_SO(FT,TEMP,NLMS,NOM)
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:G_FIL,G_F,G_MU,MU_1
      IMPLICIT NONE
C*--DMFT_FFT_FTOT2_SO1209
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER FERMION
      PARAMETER (FERMION=1)
C
C Dummy arguments
C
      INTEGER NLMS,NOM
      DOUBLE PRECISION TEMP
      COMPLEX*16 FT(NLMS,NLMS,NOM)
C
C Local variables
C
      INTEGER DIAG,IK,IL
      COMPLEX*16 F(NOM),OUT(NOM)
C
C*** End of declarations rewritten by SPAG
C
C*** Start of declarations inserted by SPAG
C*** End of declarations inserted by SPAG
C
      DO IL = 1,NLMS
         DO IK = 1,NLMS
            DIAG = 0
            IF ( IL.EQ.IK ) DIAG = 1
            F(1:NOM) = FT(IK,IL,1:NOM)
C      call dfour1(f,nom,-1)
C omega -> tau
            CALL DMFT_INVFOURIER2(F,OUT,NOM,NOM-1,DIAG,FERMION,TEMP)
C
            FT(IK,IL,1:NOM) = OUT(1:NOM)
            IF ( IK.EQ.IL ) THEN
               G_MU(IK,IL) = MU_1
               G_FIL(IK,IL) = G_F
            END IF
         END DO
      END DO
C
      END
C*==dmft_fft_ttof4_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C=============================================================================
C
      SUBROUTINE DMFT_FFT_TTOF4_SO(FT,TEMP,NLMS,NOM)
C
      IMPLICIT NONE
C*--DMFT_FFT_TTOF4_SO1273
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER FERMION
      PARAMETER (FERMION=0)
C
C Dummy arguments
C
      INTEGER NLMS,NOM
      REAL*8 TEMP
      COMPLEX*16 FT(NLMS,NLMS,NLMS,NLMS,NOM)
C
C Local variables
C
      INTEGER DIAG,IK1,IK2,IL1,IL2
      COMPLEX*16 F(NOM),OUT(0:NOM-1)
C
C*** End of declarations rewritten by SPAG
C
C
      DIAG = 0
      DO IL1 = 1,NLMS
         DO IL2 = 1,NLMS
            DO IK1 = 1,NLMS
               DO IK2 = 1,NLMS
                  F(1:NOM) = FT(IK2,IK1,IL2,IL1,1:NOM)
C
C tau -> omega
                  CALL DMFT_NFOURIER2(F,OUT,NOM,NOM-1,DIAG,FERMION,TEMP)
C
                  FT(IK2,IK1,IL2,IL1,1:NOM) = OUT(0:NOM-1)
C
               END DO
            END DO
         END DO
      END DO
C
      END
C*==dmft_fft_ftot4_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C=============================================================================
C
      SUBROUTINE DMFT_FFT_FTOT4_SO(FT,TEMP,NLMS,NOM)
C
      IMPLICIT NONE
C*--DMFT_FFT_FTOT4_SO1335
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER FERMION
      PARAMETER (FERMION=0)
C
C Dummy arguments
C
      INTEGER NLMS,NOM
      REAL*8 TEMP
      COMPLEX*16 FT(NLMS,NLMS,NLMS,NLMS,NOM)
C
C Local variables
C
      INTEGER DIAG,IK1,IK2,IL1,IL2
      COMPLEX*16 F(NOM),OUT(NOM)
C
C*** End of declarations rewritten by SPAG
C
C
      DIAG = 0
      DO IL1 = 1,NLMS
         DO IL2 = 1,NLMS
            DO IK1 = 1,NLMS
               DO IK2 = 1,NLMS
                  F(1:NOM) = FT(IK2,IK1,IL2,IL1,1:NOM)
C
C omega -> tau
                  CALL DMFT_INVFOURIER2(F,OUT,NOM,NOM-1,DIAG,FERMION,
     &                                  TEMP)
C
                  FT(IK2,IK1,IL2,IL1,1:NOM) = OUT(1:NOM)
C
               END DO
            END DO
         END DO
      END DO
C
      END
C*==dmft_nfourier2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************************
C COMMENT: INPUT: rindata(L)  - Greens' function(tau) from 0 to beta
C                 diag        - 1 for diagonal element of GF
C                             - 0(not 1) for non-diagonal element of GF
C                 fermion     - 1 for fermionic GF
C                             - 0(not 1) for bosonic
C                 temp        - temperature in Ry
C
C         OUTPUT: coutdata(0:Iwmax) -Green's function (omega)
C
C========+=========+=========+=========+=========+=========+=========+=$
      SUBROUTINE DMFT_NFOURIER2(RINDATA,COUTDATA,L,IWMAX,DG,FERM,TEMP)
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:SIGNATURE,ENDSIG,CI,PI
      IMPLICIT NONE
C*--DMFT_NFOURIER21409
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DG,FERM,IWMAX,L
      REAL*8 TEMP
      COMPLEX*16 COUTDATA(0:IWMAX),RINDATA(L)
C
C Local variables
C
      REAL*8 BETA,DELTA,DIFF,OM
      COMPLEX*16 CDUMMY,DAT(:)
      INTEGER I,NSIGN
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE DAT
      ALLOCATE (DAT(2*L))
      DIFF = 0D0
      NSIGN = 1 - 2*FERM
      BETA = 1D0/TEMP
C      BOSON = 1 - FERM
      DO I = 1,L
         DAT(I) = RINDATA(I)
         DAT(L+I) = NSIGN*RINDATA(I)
      END DO
      CALL DMFT_FFT_SO(DAT,2*L,1)
      DO I = 1,IWMAX + 1
         COUTDATA(I-1) = DAT(I)*BETA/2D0/REAL(L)
      END DO
C
C put in attenuation factors
C
      IF ( FERM.EQ.0 .AND. DG.EQ.1 ) THEN
         DIFF = 0D0
         COUTDATA(0) = COUTDATA(0) + DG*DIFF*BETA/2D0/REAL(L)
      ELSE IF ( FERM.EQ.1 .AND. DG.EQ.1 ) THEN
         DIFF = 1D0
         IF ( SIGNATURE.EQ.'S' ) DIFF = REAL(-RINDATA(1)-ENDSIG,KIND=8)
      END IF
      DO I = 2 - FERM,IWMAX,2
         OM = I*PI/BETA
         DELTA = BETA/L
         CDUMMY = CI*OM*DELTA
         COUTDATA(I) = (1D0/CI/OM-(EXP(-CDUMMY)-1.D0)/OM**2/DELTA)
     &                 *DG*DIFF + 4D0*(SIN(OM*DELTA/2D0))
     &                 **2/OM**2/DELTA**2*COUTDATA(I)
      END DO
      END
C*==dmft_invfourier2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: invfourier
C TYPE   : subroutine
C PURPOSE: inverse fourier transform F(iw) -> F(tau)
C          Greent, Greenw use physical definition
C          Greent(i)=G((i-1)*deltau) for i=1,...,T
C          Greenw(n)=G(i w_n), for n=0,T/2-1
C                 w_n=(2*n+1)pi/beta
C          Symmetry property:
C          G(iw_(-n)=G(iw_(n-1))*
C          coupled to the impurity
C I/O    :
C VERSION: 11-18-91
C COMMENT:
C========+=========+=========+=========+=========+=========+=========+=$
      SUBROUTINE DMFT_INVFOURIER2(CINDATA,RINDATA,L,IWMAX,DG,FERM,TEMP)
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:G_F,MU_1,CI,PI
      IMPLICIT NONE
C*--DMFT_INVFOURIER21496
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DG,FERM,IWMAX,L
      REAL*8 TEMP
      COMPLEX*16 CINDATA(0:IWMAX),RINDATA(L)
C
C Local variables
C
      REAL*8 ASSIMPT,BETA,CHARGE,DELTA,TAU,W_MAX
      COMPLEX*16 CINPRIME(2*IWMAX+2),G_CORR,I_W
      INTEGER DIAG,I
C
C*** End of declarations rewritten by SPAG
C
C
      BETA = 1D0/TEMP
      DELTA = 1D0*BETA/L
      ASSIMPT = 0D0
      CINPRIME = (0D0,0D0)
      IF ( DG.NE.1 ) DIAG = 0
      IF ( DG.EQ.1 .AND. FERM.EQ.1 ) THEN
C
C  Assimptotics 1/(iw-mu_1) is used. mu_1 is fitted in such way so
C  in tau space  (g-1/(iw-mu_1))(tau=0)=0
C
         DIAG = 1
         W_MAX = IWMAX*PI/BETA
         CHARGE = 1D0 + (2D0*SUM(REAL(CINDATA))/BETA-0.5D0)
C
         IF ( 1D0-CHARGE.LT.1D-5 ) THEN
C Filled band
            G_F = 1
         ELSE IF ( CHARGE.LT.1D-5 ) THEN
C Empty band
            G_F = 0
         ELSE
C Some filling
            G_F = 2
            MU_1 = -TEMP*LOG(-1D0/(CHARGE-1D0)-1D0)
C        CALL newt(cindata,Iwmax+1,mu_1,pi,temp)
            G_CORR = 1D0/(CI*W_MAX-MU_1)
            DO I = 1,IWMAX,2
               I_W = I*PI/BETA*CI
               G_CORR = 1D0/(I_W-MU_1)
               CINDATA(I) = CINDATA(I) - G_CORR
            END DO
         END IF
      ELSE
         DIAG = 0
      END IF
C
      DO I = 1,IWMAX + 1
         CINPRIME(I) = CINDATA(I-1)
      END DO
      IF ( FERM.EQ.0 ) THEN
         CINPRIME(2*IWMAX+2) = (0D0,0D0)
         CINPRIME(IWMAX+2) = REAL(CINDATA(IWMAX-1),KIND=8)
      END IF
      DO I = 1,IWMAX - 1 + FERM
         CINPRIME(2*IWMAX+2+FERM-I) = CONJG(CINDATA(I+1-FERM))
      END DO
C
      CALL DMFT_FFT_SO(CINPRIME,2*IWMAX+2,-1)
C
      DO I = 1,L
         IF ( DIAG.EQ.1 .AND. FERM.EQ.1 ) THEN
            IF ( G_F.EQ.2 ) THEN
C Some filling
               TAU = DELTA*(I-1)
               IF ( BETA*MU_1.GT.0D0 ) THEN
                  ASSIMPT = -EXP(-TAU*MU_1)/(1D0+EXP(-BETA*MU_1))
               ELSE
                  ASSIMPT = -EXP((BETA-TAU)*MU_1)/(1D0+EXP(BETA*MU_1))
               END IF
            ELSE IF ( G_F.EQ.0 ) THEN
C Empty band
               ASSIMPT = 0D0
            ELSE IF ( G_F.EQ.1 ) THEN
C Filled band
               ASSIMPT = 0D0
               IF ( I.EQ.1 ) ASSIMPT = -1D0
            END IF
         END IF
         RINDATA(I) = CINPRIME((I-1)*(IWMAX+1)/L+1)/BETA + ASSIMPT
C        write(78,'(4f14.7)')tau,
C     ,      REAL(cinprime((i-1)*(Iwmax+1)/L+1)/beta),assimpt,
C     ,      REAL(rindata(i))
      END DO
      END
C*==dmft_fft_so.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C
C
C========+=========+=========+=========+=========+=========+=========+=$
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     FFT routine from Numerical recipes slightly changed.
C
C     NOTE, THAT nn MUST BE AN INTEGER POWER OF 2!
C
C     "data" carries the complex input and output vectors of length nn.
C     It does so by being a real array of length 2*nn, where two consecutive
C     components are occupied by one complex number.
C
C     If "isign" is input as 1, the FFT is computed [Eq. (12.1.7) on page
C     497 of the Numerical Recipes]. If "isign" is input as -1, the iFFT
C     is computed [Eq. (12.1.9) with the prefactor 1/N being discarded].
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
      SUBROUTINE DMFT_FFT_SO(CDATAM,NN,JSIGN)
      IMPLICIT NONE
C*--DMFT_FFT_SO1625
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER JSIGN,NN
      COMPLEX*16 CDATAM(NN)
C
C Local variables
C
      INTEGER I,ISTEP,J,M,MMAX,N
      REAL*8 RDATA(:),TEMPI,TEMPR,THETA,WI,WPI,WPR,WR,WTEMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RDATA
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (RDATA(2*NN))
      N = 2*NN
      J = 0
      DO I = 1,NN
         J = J + 1
         RDATA(J) = DREAL(CDATAM(I))
         J = J + 1
         RDATA(J) = DIMAG(CDATAM(I))
      END DO
      J = 1
      DO I = 1,N,2
         IF ( J.GT.I ) THEN
            TEMPR = RDATA(J)
            TEMPI = RDATA(J+1)
            RDATA(J) = RDATA(I)
            RDATA(J+1) = RDATA(I+1)
            RDATA(I) = TEMPR
            RDATA(I+1) = TEMPI
         END IF
         M = N/2
 50      CONTINUE
         IF ( (M.GE.2) .AND. (J.GT.M) ) THEN
            J = J - M
            M = M/2
            GOTO 50
         END IF
         J = J + M
      END DO
      MMAX = 2
 100  CONTINUE
      IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         THETA = 6.28318530717959D0/(JSIGN*MMAX)
         WPR = -2.D0*SIN(0.5D0*THETA)**2
         WPI = SIN(THETA)
         WR = 1.D0
         WI = 0.D0
         DO M = 1,MMAX,2
            DO I = M,N,ISTEP
               J = I + MMAX
               TEMPR = WR*RDATA(J) - WI*RDATA(J+1)
               TEMPI = WR*RDATA(J+1) + WI*RDATA(J)
               RDATA(J) = RDATA(I) - TEMPR
               RDATA(J+1) = RDATA(I+1) - TEMPI
               RDATA(I) = RDATA(I) + TEMPR
               RDATA(I+1) = RDATA(I+1) + TEMPI
            END DO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         END DO
         MMAX = ISTEP
         GOTO 100
      END IF
      J = 1
      DO I = 1,NN
         CDATAM(I) = DCMPLX(RDATA(J),RDATA(J+1))
         J = J + 2
      END DO
      END
C*==dmft_matsub_sosptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************
      SUBROUTINE DMFT_MATSUB_SOSPTFLEX
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NOM,OMEGA,TIME,DOMEGA,DTIME,MINOM
      IMPLICIT NONE
C*--DMFT_MATSUB_SOSPTFLEX1725
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IOM
C
C*** End of declarations rewritten by SPAG
C
C
C
C..... Ordering of the "periodic" bosonic Matsubara frequencies:
C      0,2piT,4piT,..... Nom*2piT
C      For fermionic use omega(i)+domega/2
      DO IOM = 1,NOM
         OMEGA(IOM) = DOMEGA*(IOM-1)
         TIME(IOM) = DTIME*(IOM-1)
      END DO
C------- Define "minus" array for Tau ----
C------- if i=0,...,N-1 then  -i=N-i
C------- if i=1,...,N   then  -i=N+2-i
      MINOM(1) = 1
      DO IOM = 2,NOM
         MINOM(IOM) = NOM + 2 - IOM
C         write(*,'(2i7,1X,f12.8)') iom,minom(iom),omega(iom)
      END DO
      END
C*==dmft_vertexsp_sosptflex.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************
C
      SUBROUTINE DMFT_VERTEXSP_SOSPTFLEX(U)
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:NLM,UC
      IMPLICIT NONE
C*--DMFT_VERTEXSP_SOSPTFLEX1770
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 U(NLM,NLM,NLM,NLM)
C
C Local variables
C
      INTEGER IS1,IS2
C
C*** End of declarations rewritten by SPAG
C
C
C      REAL(KIND=8), DIMENSION(nlms,nlms,nlms,nlms) :: U
C-------------- Vertex
C
C------------- Uc - just COMPLEX U
C If dim[U]=(nlm*ns X nlm*ns)
C       Uc=U
C If dim[U]=(nlm X nlm)
      UC = (0D0,0D0)
      DO IS1 = 1,2
         DO IS2 = 1,2
            UC((IS1-1)*NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM,(IS1-1)
     &         *NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM) = U(:,:,:,:)
         END DO
      END DO
      END
C*==g_tau.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************
C
      REAL(KIND=8) FUNCTION G_TAU(M1,M2,I,TEMP)
C
C Returns analytical value G(tau) = -exp(-t*mu)/(exp(-beta*mu)+1)
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:DTIME,G_FIL,G_MU
      IMPLICIT NONE
C*--G_TAU1823
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,M1,M2
      REAL*8 TEMP
C
C Local variables
C
      REAL*8 BETA,MU,TAU
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( G_FIL(M1,M2).EQ.1 ) THEN
C Filled band
         G_TAU = 0D0
         IF ( I.EQ.1 ) G_TAU = -1D0
         RETURN
      ELSE IF ( G_FIL(M1,M2).EQ.0 ) THEN
C Empty band
         G_TAU = 0D0
         RETURN
      END IF
C
      TAU = (I-1)*DTIME
      MU = G_MU(M1,M2)
      BETA = 1D0/TEMP
      IF ( BETA*MU.GT.0D0 ) THEN
         G_TAU = -EXP(-TAU*MU)/(1D0+EXP(-BETA*MU))
      ELSE
         G_TAU = -EXP((BETA-TAU)*MU)/(1D0+EXP(BETA*MU))
      END IF
      END
C*==g_mtau.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************
C
      REAL(KIND=8) FUNCTION G_MTAU(M1,M2,I,TEMP)
C
C Returns analytical value G(-tau) = exp(-(beta-tau)*mu)/(exp(-beta*mu)+1)
C
C
      USE MOD_DMFT_SOSPTFLEX,ONLY:G_FIL,DTIME,G_MU
      IMPLICIT NONE
C*--G_MTAU1882
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,M1,M2
      REAL*8 TEMP
C
C Local variables
C
      REAL*8 BETA,MU,TAU
C
C*** End of declarations rewritten by SPAG
C
C
      IF ( G_FIL(M1,M2).EQ.1 ) THEN
C Filled band
         G_MTAU = 0D0
         RETURN
      ELSE IF ( G_FIL(M1,M2).EQ.0 ) THEN
C Empty band
         G_MTAU = 0D0
         IF ( I.EQ.1 ) G_MTAU = 1D0
         RETURN
      END IF
C
      TAU = (I-1)*DTIME
      MU = G_MU(M1,M2)
      BETA = 1D0/TEMP
      IF ( BETA*MU.GT.0D0 ) THEN
         G_MTAU = EXP(-(BETA-TAU)*MU)/(1D0+EXP(-BETA*MU))
      ELSE
         G_MTAU = EXP(TAU*MU)/(1D0+EXP(BETA*MU))
      END IF
      END
C*==dmft_tmat2.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
C***********************************************************
C
      SUBROUTINE DMFT_TMAT2(R0,U)
      USE MOD_DMFT_SOSPTFLEX,ONLY:NN
      IMPLICIT NONE
C*--DMFT_TMAT21937
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 R0(NN,NN),U(NN,NN)
C
C Local variables
C
      COMPLEX*16 D(NN,NN)
C
C*** End of declarations rewritten by SPAG
C
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),R0,NN,U,NN,
     &           DCMPLX(0.D0,0.D0),D,NN)
      CALL ZGEMM('N','N',NN,NN,NN,DCMPLX(1.D0,0.D0),U,NN,D,NN,
     &           DCMPLX(0.D0,0.D0),R0,NN)
      END
C*==dmft_newt.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C***********************************************************
C
      SUBROUTINE DMFT_NEWT(G,NOM,MU,PI,TEMP)
C
C Solve equation 2Re SUM_w[G(iw)-1/(iw+mu)] = 0 by the Newton method
C
      IMPLICIT NONE
C*--DMFT_NEWT1976
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER MAXITER
      PARAMETER (MAXITER=100)
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
C
C Dummy arguments
C
      REAL*8 MU,PI,TEMP
      INTEGER NOM
      COMPLEX*16 G(NOM)
C
C Local variables
C
      REAL*8 COFF,G0,GDER,GV,MU_OLD,W
      COMPLEX*16 G_SUM
      INTEGER I,ITER
C
C*** End of declarations rewritten by SPAG
C
C Sum up the real part of G
      G_SUM = (0D0,0D0)
      DO I = 1,NOM/2
         G_SUM = G_SUM + G(2*I)
      END DO
      G0 = REAL(G_SUM,KIND=8)
C
      DO ITER = 1,MAXITER
         MU_OLD = MU
         GV = 0D0
         GDER = 0D0
         DO I = 1,NOM,2
            W = I*PI*TEMP
            COFF = W*W + MU*MU
            GV = GV + MU/COFF
            GDER = GDER + (W*W-MU*MU)/COFF/COFF
         END DO
         MU = MU_OLD - (GV+G0)/GDER
         IF ( ABS(MU-MU_OLD).LT.TOL ) GOTO 99999
      END DO
      WRITE (6,*) 'NEWT: MU is not found'
      STOP
99999 CONTINUE
      END
C*==dmft_simpn.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
C
      SUBROUTINE DMFT_SIMPN(F,D,N,FI)
C     ***************************************************************
C     *                                                             *
C     * Integrates F from X(1) to X(n) by the Simpson's technique   *
C     *                                                             *
C     ***************************************************************
      IMPLICIT NONE
C*--DMFT_SIMPN2049
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 D,FI
      INTEGER N
      REAL*8 F(N)
C
C Local variables
C
      REAL*8 FF
      INTEGER J,N2
C
C*** End of declarations rewritten by SPAG
C
      N2 = N - 1
      IF ( MOD(N,2).EQ.0 ) THEN
         WRITE (6,99001) N
         STOP
      END IF
      FF = 0.D0
      DO J = 2,N2,2
         FF = FF + (F(J-1)+4.D0*F(J)+F(J+1))
      END DO
      FI = D*FF/3.D0
99001 FORMAT (' SIMPN:**  N =',I4,'. Must be even')
      END
C
