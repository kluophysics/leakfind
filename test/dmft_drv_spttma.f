C*==dmft_drv_spttma.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_DRV_SPTTMA(ITRSCF,NKMMAX,NKM,NEGF,DMFTMIX,LOP,
     &                           EGFOLD,EGFNEW,EREFDMFT,WEGF,GFMAT,
     &                           DMFTSIGMA,TXT_T,LTXT_T,EFERMI0,IPRINT,
     &                           UDMFT,NLMAX,DMFTDBLC,TMA_SPINFLIPOFF,
     &                           TMA_REFACTOR,TMA_STATIC,UEFF,JEFF,
     &                           TMA_DBLC,TMA_NOBATH,IT,NTMAX,DOSINT,
     &                           TMA_NOCCSCL,LDAUSIGMA)
C
C   ********************************************************************
C   * SP-TMA DMFT solver working on the real axis
C   * Stanislav Chadov
C   ********************************************************************
C   * Added features :
C   *
C   *  LDAU (set Sigma(EF)=0,  + LDA+U(AMF/AAL))
C   *
C   *  TMA_SPINFLIPOFF:
C   *  if specified, the spin-flip terms of the input GF are neglected
C   *
C   *  TMA_REFACTOR:
C   *  Re[Sigma] scaling factor (applied after Re[Sigma]
C   *  via KKT have been obtained)
C   *  to keep the EF unchanged within META mode
C   *  (at the moment the same for all correlated sorts)
C   *  calculation of occupation numbers OCC,
C   *  used for SIGSTAT (static self-energy)
C   *
C   *  TMA_STATIC:
C   *  if specified, the dynamical part is ignored; added
C   *  only LSDA+U
C   *  We imply that input GF-matrix is symmetric
C   *  ->  Im(GF) is symmetric -> C [0,Im(GF)] -
C   *  is antihermitian component.
C   *
C   *  TMA_DBLC:
C   *  TMA_DBLC=AMF   /LSDA+U "Around Mean-Field"/,
C   *  TMA_DBLC=AAL   /LSDA+U "Around Atomic Limit"/
C   *   Make sense only if DMFTDBLC=LDAU ! By default TMA_DBLC=AMF
C   *
C   *  TMA_NOBATH:  /testing/
C   *  if specified, the effective GF instead of bath is used in the solver
C   *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0,PI
      USE MOD_DMFT_LDAU,ONLY:DMFT_FIX_DYN_SE
      IMPLICIT NONE
C*--DMFT_DRV_SPTTMA48
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1.0D-10)
C
C Dummy arguments
C
      CHARACTER*4 DMFTDBLC
      REAL*8 DMFTMIX,EFERMI0,JEFF,TMA_REFACTOR,UEFF
      COMPLEX*16 EREFDMFT
      INTEGER IPRINT,IT,ITRSCF,LOP,LTXT_T,NEGF,NKM,NKMMAX,NLMAX,NTMAX
      CHARACTER*3 TMA_DBLC
      LOGICAL TMA_NOBATH,TMA_NOCCSCL,TMA_SPINFLIPOFF,TMA_STATIC
      CHARACTER*8 TXT_T
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NEGF),EGFNEW(NEGF),EGFOLD(NEGF)
     &           ,GFMAT(NKMMAX,NKMMAX,NEGF),LDAUSIGMA(NKMMAX,NKMMAX),
     &           WEGF(NEGF)
      REAL*8 DOSINT(NLMAX),UDMFT(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX)
C
C Local variables
C
      COMPLEX*16 CTMP,CWORK1(:),CWORK1P(:),CWORK2(:,:),CWORK2P(:,:),
     &           GG(:,:,:),GR(:,:,:),GZ(:,:,:),SEZ(:,:,:,:),
     &           SEZNEW(:,:,:),SE_AR(:,:,:),SE_AR_SAV(:,:,:)
      REAL*8 DX,EFERMINEW,EFERMIOLD,OCC(:,:),RWORK1(:),RWORK1P(:),
     &       RWORK2(:,:),RWORK4(:,:,:,:),SIGSTAT(:,:),TMP,V(:,:),X(:),
     &       XX(:),XX_SAV(:)
      CHARACTER*80 FILNAM
      INTEGER I,I1,I2,I3,I4,ICALL,IE,IND(:,:),IR,IR1,IR2,IR3,IR4,IRP,
     &        IZERO(:,:),J,LFILNAM,M1,M2,M3,M4,NLM,NR,NREP,NRR
      INTEGER IORB,ISPN
      LOGICAL IPRT,WRSIGKK
      SAVE SEZ,SE_AR_SAV,XX_SAV
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C
C
C
C
C
      DATA ICALL/0/
C
      ALLOCATABLE GZ,SEZ,GR,GG,SEZNEW,SE_AR,V,OCC,SIGSTAT
      ALLOCATABLE X,XX,RWORK1,RWORK1P
      ALLOCATABLE RWORK2,RWORK4
      ALLOCATABLE CWORK1,CWORK1P,CWORK2,CWORK2P
      ALLOCATABLE IZERO,IND,XX_SAV,SE_AR_SAV
C
C
C
C=======================================================
C Set initial paramters:
C=======================================================
      NR = 501
      NRR = 1001
      IPRT = .FALSE.
      WRSIGKK = .TRUE.
      IF ( IPRINT.GT.0 ) IPRT = .TRUE.
C
C
      EFERMIOLD = EFERMI0
      EFERMINEW = DREAL(EGFNEW(NEGF))
C
C      UUU = UEFF/13.605826D0
C      JJJ = JEFF/13.605826D0
      NLM = LOP*2 + 1
      I1 = LOP**2 + 1
      I2 = I1 + NLM - 1
      I3 = NKM/2 + I1
      I4 = I3 + NLM - 1
      NREP = 2*NLM
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Create double index:
C=======================================================
      ALLOCATE (IND(NREP,NREP))
      I = 0
      IND(1:NREP,1:NREP) = 0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            I = I + 1
            IND(IR,IRP) = I
         END DO
      END DO
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C  Set Sigma to 0 if it is the 1st iteration:
C=======================================================
      IF ( DMFTMIX.LE.TOL .AND. ITRSCF.GT.1 ) THEN
         ALLOCATE (SEZNEW(NREP,NREP,NEGF))
         SEZNEW(1:NREP,1:NREP,1:NEGF) = C0
         GOTO 100
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Checking printouts:
C=======================================================
      IF ( IPRT ) THEN
         WRITE (*,'(/,8X,A)') '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
         WRITE (*,'(/,10X,A,/)') 'DMFT-TMA input parameters:'
         WRITE (*,'(/,10X,A,A,/)') 'system suffix:         : ',
     &                             TXT_T(1:LTXT_T)
         WRITE (6,'(/,10X,A,I4)') 'NLMAX                   : ',NLMAX
         WRITE (6,'(/,10X,A,I4)') 'correlated shell        : ',LOP
         WRITE (6,'(/,10X,A,I4)') 'problem dimension       : ',NREP
         WRITE (6,'(/,10X,A,I4)') 'number of energy points : ',NEGF
         WRITE (6,'(/,10X,A,F6.3)') 'E_min [Ry]              : ',
     &                              REAL(EGFNEW(1))
         WRITE (6,'(/,10X,A,F6.3)') 'E_F [Ry]                : ',EFERMI0
         WRITE (6,'(/,10X,A,F6.3)') 'E_max [Ry]              : ',
     &                              REAL(EGFNEW(NEGF))
         IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
            IF ( TMA_DBLC(1:3).EQ.'AMF' ) WRITE (6,'(/,10X,A)')
     &            'LSDA+U with AMF /around mean-field/ double-counting'
            IF ( TMA_DBLC(1:3).EQ.'AAL' ) WRITE (6,'(/,10X,A)') 
     &           'LSDA+U with AAL /around atomic limit/ double-counting'
         END IF
         WRITE (*,'(/,8X,A,/)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Here in GZ later we store the bath GF G_BATH=G^-1+SIGMA:
C=======================================================
      IF ( ICALL.EQ.0 ) THEN
         ALLOCATE (SEZ(NREP,NREP,NEGF,NTMAX))
         SEZ = C0
      END IF
C
      ALLOCATE (GZ(NREP,NREP,NEGF))
      GZ(1:NREP,1:NREP,1:NEGF) = C0
C
      GZ(1:NLM,1:NLM,1:NEGF) = GFMAT(I1:I2,I1:I2,1:NEGF)         !1_1
      GZ(NLM+1:NREP,NLM+1:NREP,1:NEGF) = GFMAT(I3:I4,I3:I4,1:NEGF)
                                                                 !2_2
      GZ(1:NLM,NLM+1:NREP,1:NEGF) = GFMAT(I1:I2,I3:I4,1:NEGF)    !1_2
      GZ(NLM+1:NREP,1:NLM,1:NEGF) = GFMAT(I3:I4,I1:I2,1:NEGF)    !2_1
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C  Check which elements of GZ are nonzero:
C=======================================================
      IF ( TMA_SPINFLIPOFF ) THEN
         DO IR = 1,NREP
            DO IRP = 1,NREP
               IF ( (IR.GT.NREP/2 .AND. IRP.LE.NREP/2) .OR. 
     &              (IRP.GT.NREP/2 .AND. IR.LE.NREP/2) )
     &              GZ(IRP,IR,1:NEGF) = C0
            END DO
         END DO
      END IF
C
      ALLOCATE (IZERO(NREP,NREP))
C
      IZERO(1:NREP,1:NREP) = 0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            CTMP = C0
            DO IE = 1,NEGF
               CTMP = CTMP + GZ(IRP,IR,IE)
            END DO
            IF ( CDABS(CTMP).GE.TOL ) IZERO(IRP,IR) = 1
         END DO
      END DO
C
C
      IF ( IPRT ) THEN
         ALLOCATE (CWORK2(NREP,NREP))
         CWORK2(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,NEGF-1)
         CALL CMATSTRUCT('INPUT GF-MATRIX (AS READ):',CWORK2,NREP,NREP,
     &                   0,0,0,TOL,6)
         DEALLOCATE (CWORK2)
      END IF
C
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Calculate bath GF:
C=======================================================
      IF ( ITRSCF.GT.1 .AND. .NOT.TMA_NOBATH ) THEN
         ALLOCATE (CWORK2(NREP,NREP),CWORK2P(NREP,NREP))
         DO IE = 1,NEGF
            CWORK2(1:NREP,1:NREP) = GZ(1:NREP,1:NREP,IE)
            CALL CMATINV(NREP,NREP,CWORK2,CWORK2P)
            CWORK2(1:NREP,1:NREP) = CWORK2P(1:NREP,1:NREP)
     &                              + SEZ(1:NREP,1:NREP,IE,IT)
            CALL CMATINV(NREP,NREP,CWORK2,CWORK2P)
            GZ(1:NREP,1:NREP,IE) = CWORK2P(1:NREP,1:NREP)
         END DO
         DEALLOCATE (CWORK2,CWORK2P)
      ELSE IF ( TMA_NOBATH ) THEN
         WRITE (*,'(/,10X,A)') 'Without bath GF! (testing case)'
      END IF
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C Static part of the self-energy:
C=======================================================
C
C Calculate occupation numbers:
C
      IF ( DMFTDBLC.EQ.'LDAU' ) THEN
C
         ALLOCATE (OCC(NREP,NREP))
         TMP = 0D0
         OCC(1:NREP,1:NREP) = 0D0
         DO IR = 1,NREP
            DO IRP = 1,NREP
               IF ( IZERO(IR,IRP).EQ.1 ) THEN
                  DO IE = 1,NEGF
                     OCC(IR,IRP) = OCC(IR,IRP)
     &                             - IMAG(GZ(IR,IRP,IE)*WEGF(IE))
                  END DO
               END IF
            END DO
            TMP = TMP + OCC(IR,IR)
         END DO
         OCC(1:NREP,1:NREP) = OCC(1:NREP,1:NREP)/PI
         TMP = TMP/PI
C
         IF ( IPRT ) THEN
            CALL RMATSTRUCT('occupation numbers:',OCC,NREP,NREP,0,0,0,
     &                      TOL,6)
            WRITE (*,'(/,10X,A,F8.4)') 'Trace:',TMP
         END IF
C
         IF ( TMA_NOCCSCL ) THEN
            TMP = DOSINT(LOP+1)/TMP
            DO IR = 1,NREP
               OCC(IR,IR) = OCC(IR,IR)*TMP
            END DO
            IF ( IPRT ) THEN
               CALL RMATSTRUCT('occupation numbers RESCALED:',OCC,NREP,
     &                         NREP,0,0,0,TOL,6)
               TMP = 0.0D0
               DO IR = 1,NREP
                  TMP = TMP + OCC(IR,IR)
               END DO
               WRITE (*,'(/,20X,A,F8.4)') 'NOCC from SCF:',DOSINT(LOP+1)
               WRITE (*,'(/,20X,A,F8.4)') 'Trace: RESCALED',TMP
            END IF
         END IF
C
C
C
      END IF  ! DMFTDBLC='LDAU'
C
C
C Change the index coupling order according to Eschrig:
C     U_m1m2m3m4 = <m1|m3><m2|m4>:
      ALLOCATE (RWORK4(NREP/2,NREP/2,NREP/2,NREP/2))
      DO IR1 = 1,NREP/2
         DO IR2 = 1,NREP/2
            DO IR3 = 1,NREP/2
               DO IR4 = 1,NREP/2
                  RWORK4(IR1,IR2,IR3,IR4) = UDMFT(IR1,IR3,IR2,IR4)
               END DO
            END DO
         END DO
      END DO
      UDMFT(1:NREP/2,1:NREP/2,1:NREP/2,1:NREP/2)
     &   = RWORK4(1:NREP/2,1:NREP/2,1:NREP/2,1:NREP/2)
      DEALLOCATE (RWORK4)
C
C printout UDMFT matrix:
      ALLOCATE (RWORK2((NREP/2)**2,(NREP/2)**2))
C
      RWORK2(1:(NREP/2)**2,1:(NREP/2)**2) = 0D0
      I = 0
      DO M1 = 1,2*LOP + 1
         DO M2 = 1,2*LOP + 1
            I = I + 1
            J = 0
            DO M3 = 1,2*LOP + 1
               DO M4 = 1,2*LOP + 1
                  J = J + 1
                  RWORK2(I,J) = UDMFT(M1,M2,M3,M4)
               END DO
            END DO
         END DO
      END DO
      IF ( IPRT ) CALL RMATSTRUCT('UDMFT matrix (with F_0<>0):',RWORK2,
     &                            (NREP/2)**2,(NREP/2)**2,0,0,0,TOL,6)
C
      DEALLOCATE (RWORK2)
C
C
C  Read in UDMFT matrix for each (s1,s2,s1,s2)-block
C  and comprise the total V-matrix:
      ALLOCATE (V(NREP**2,NREP**2))
      V(1:NREP**2,1:NREP**2) = 0D0
      DO IR1 = 1,NREP
         DO IR2 = 1,NREP
            DO IR3 = 1,NREP
               DO IR4 = 1,NREP
                  IF ( ISPN(IR1,NREP).EQ.ISPN(IR3,NREP) .AND. 
     &                 ISPN(IR2,NREP).EQ.ISPN(IR4,NREP) )
     &                 V(IND(IR1,IR2),IND(IR3,IR4))
     &                 = UDMFT(IORB(IR1,NREP)+LOP+1,IORB(IR2,NREP)
     &                 +LOP+1,IORB(IR3,NREP)+LOP+1,IORB(IR4,NREP)+LOP+1)
               END DO
            END DO
         END DO
      END DO
C
C
C Calculate static self-energy /LSDA+U: Around Mean-Field part:
C  Czyzik and Zawatzky, PRB49(1994)14211: formula (4)/:
C TMP is taken as orbital average for each spin:
C
      IF ( DMFTDBLC(1:4).EQ.'LDAU' ) THEN
         ALLOCATE (SIGSTAT(NREP,NREP))
         SIGSTAT(1:NREP,1:NREP) = 0D0
         I1 = LOP**2 + 1
         I2 = I1 + NLM - 1
         I3 = NKM/2 + I1
         I4 = I3 + NLM - 1
         SIGSTAT(1:NLM,1:NLM) = DREAL(LDAUSIGMA(I1:I2,I1:I2))
         SIGSTAT(NLM+1:NREP,NLM+1:NREP) = DREAL(LDAUSIGMA(I3:I4,I3:I4))
         SIGSTAT(1:NLM,NLM+1:NREP) = DREAL(LDAUSIGMA(I1:I2,I3:I4))
         SIGSTAT(NLM+1:NREP,1:NLM) = DREAL(LDAUSIGMA(I3:I4,I1:I2))
         IF ( IPRINT.GT.0 ) CALL RMATSTRUCT('NEWtic self-energy:',
     &        SIGSTAT,NREP,NREP,0,0,0,TOL,6)
      END IF
C
C      IF ( DMFTDBLC.EQ.'LDAU' ) THEN
CC
C         ALLOCATE (SIGSTAT(NREP,NREP))
C         SIGSTAT(1:NREP,1:NREP) = 0D0
CC
C         TMP1 = 0D0
C         TMP2 = 0D0
C         DO IR = 1,NREP/2
C            TMP1 = TMP1 + OCC(IR,IR)
C            TMP2 = TMP2 + OCC(IR+NREP/2,IR+NREP/2)
C         END DO
C         TMP1 = 2*TMP1/NREP
C         TMP2 = 2*TMP2/NREP
CC
C         DO IR = 1,NREP
C            DO IRP = 1,NREP
CC full static electron interaction:
CC     SIGMA_{12}=\sum_{34}(V_{1324}-V_{1342})*n_{34}
C               DO IR1 = 1,NREP
C                  DO IR2 = 1,NREP
C                     SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP)
C     &                                 + (V(IND(IR,IR1),IND(IRP,IR2))
C     &                                 -V(IND(IR,IR1),IND(IR2,IRP)))
C     &                                 *OCC(IR1,IR2)
C                  END DO
C               END DO
CC subtract double-counting /Tr<SIGMA>/:
C               IF ( IR.EQ.IRP ) THEN
C                  IF ( TMA_DBLC(1:3).EQ.'AMF' ) THEN
CC        AMF case:
CC        SIGMA_{11}=SIGMA_{11}-\sum_{2}(V_{1212}-V_{1221})*<n_{2}>
C                     DO IR1 = 1,NREP/2
C                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
C     &                     - (V(IND(IR,IR1),IND(IR,IR1))
C     &                     -V(IND(IR,IR1),IND(IR1,IR)))*TMP1
C                     END DO
C                     DO IR1 = NREP/2 + 1,NREP
C                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
C     &                     - (V(IND(IR,IR1),IND(IR,IR1))
C     &                     -V(IND(IR,IR1),IND(IR1,IR)))*TMP2
C                     END DO
C                  ELSE IF ( TMA_DBLC(1:3).EQ.'AAL' ) THEN
CC        AAL case:
CC        SIGMA_{11}=SIGMA_{11}-U(N-1/2)+J(N_{1}-1/2)
C                     SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
C     &                                - UUU*((TMP1+TMP2)*NREP/2-0.5D0)
C                     IF ( IR.LE.NREP/2 ) THEN
C                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
C     &                     + JJJ*(TMP1*NREP/2-0.5D0)
C                     ELSE
C                        SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
C     &                     + JJJ*(TMP2*NREP/2-0.5D0)
C                     END IF
C                  ELSE
C                     STOP 'DMFT_DRV_SPTTMA: CHECK TMA_DBLC !'
C                  END IF
C               END IF
C                  ! IR=IRP
C            END DO
C         END DO
CC
C         DEALLOCATE (OCC)
CC
C         IF ( IPRT ) CALL  RMATSTR('Static self-energy:',SIGSTAT,NREP,
C     &                            NREP,0,0,0,TOL,6)
CC
C      END IF                    ! DMFTDBLC='LDAU'
C=======================================================
C
C
C                   *    *    *
C
C
C=======================================================
C     Set energy mesh on real-axis:
C=======================================================
C
C Create Green's function GR(NR,NREP,NREP) on Re-axis X(NR)
C by Pade approximation:
C
      ALLOCATE (X(NR))
C
      X(1) = REAL(EGFOLD(1)+EGFOLD(NEGF))/2
      X(NR) = REAL(2*EGFOLD(NEGF)-X(1))
      DX = (X(NR)-X(1))/(NR-1)
      DO IE = 2,NR - 1
         X(IE) = X(1) + DX*(IE-1)
      END DO
C
      ALLOCATE (GR(NRR,NREP,NREP))
      ALLOCATE (CWORK2(NEGF,NEGF))
      ALLOCATE (CWORK1(NEGF))
C
      GR(1:NRR,1:NREP,1:NREP) = C0
      DO IR = 1,NREP
         DO IRP = IR,NREP
            IF ( IZERO(IRP,IR).EQ.1 ) THEN
               CWORK1(1:NEGF) = GZ(IRP,IR,1:NEGF)
               CALL PADECOEFF(CWORK1,EGFOLD,CWORK2,NEGF,NEGF)
               DO IE = 1,NR
                  CALL PADEAPPROX(CTMP,DCMPLX(X(IE),0D0),EGFOLD,CWORK2,
     &                            NEGF,NEGF)
                  GR(IE,IRP,IR) = CTMP
                  GR(IE,IR,IRP) = CTMP
               END DO
            END IF
            IF ( IRP.EQ.IR ) THEN
               DO IE = 1,NR
                  IF ( DIMAG(GR(IE,IR,IR)).GT.0D0 ) GR(IE,IR,IR)
     &                 = DCMPLX(DREAL(GR(IE,IR,IR)),0D0)
               END DO
            END IF
         END DO
      END DO
      DEALLOCATE (CWORK2,CWORK1,GZ)
C
C
C check printout:
      IF ( IPRT ) THEN
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_GF_PADE_RE.DAT')
         DO IE = 1,NR
            WRITE (1,'(1000E15.5)') X(IE),
     &                              ((DREAL(GR(IE,IR,IRP)),IRP=IR,IR),
     &                              IR=1,NREP)
         END DO
         CLOSE (1)
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_GF_PADE_IM.DAT')
         DO IE = 1,NR
            WRITE (1,'(1000E15.5)') X(IE),
     &                              ((DIMAG(GR(IE,IR,IRP)),IRP=IR,IR),
     &                              IR=1,NREP)
         END DO
         CLOSE (1)
      END IF
C
C interpolation to the fine mesh {X(NR),GR(NR)}-->{XX(NRR),GG(NRR)}:
C new mesh XX is centered at EFERMIOLD
C
      ALLOCATE (XX(NRR),GG(NREP,NREP,NRR))
      IF ( .NOT.ALLOCATED(XX_SAV) ) ALLOCATE (XX_SAV(NRR))
      GG(1:NREP,1:NREP,1:NRR) = C0
C
      XX(1) = X(1) - (EFERMIOLD-X(1))*1.5D0
      XX(NRR) = X(NR) + (EFERMIOLD-X(1))*1.5D0
      DX = (XX(NRR)-XX(1))/(NRR-1D0)
      DO IE = 2,NRR - 1
         XX(IE) = XX(1) + DX*(IE-1D0)
      END DO
C
C interpolation via KKT:
C
C$$$      ALLOCATE(CWORK1P(NRR))
C$$$      CWORK1P(:)=DCMPLX(XX(:),TOL)
C$$$      ALLOCATE(RWORK1(NR),CWORK1(NRR))
C$$$      DO IR=1,NREP
C$$$      DO IRP=1,NREP
C$$$         IF(IZERO(IRP,IR).eq.1) THEN
C$$$            RWORK1(:)=DIMAG(GR(:,IRP,IR))
C$$$            CALL KKT_EXTRAPOL(NR,X(1),RWORK1,NRR,CWORK1P,CWORK1)
C$$$               GG(IRP,IR,:)=dcmplx(0d0,dimag(CWORK1(:)))
C$$$         ENDIF
C$$$      ENDDO
C$$$      ENDDO
C$$$      DEALLOCATE(CWORK1P,RWORK1,CWORK1)
C$$$
C
C quadratic spline interpolation:
C
C$$$      ALLOCATE(RWORK1(NR),RWORK1P(NR),RWORK1PP(NR),RWORK1PPP(NRR))
C$$$      DO IR=1,NREP
C$$$      DO IRP=IR,NREP
C$$$         IF(IZERO(IRP,IR).eq.1) THEN
C$$$            RWORK1(:)=DIMAG(GR(:,IRP,IR))
C$$$            if(IR.eq.IRP) then
C$$$               do ie=1,NR
C$$$                  if(RWORK1(ie).gt.0d0) RWORK1(ie)=0d0
C$$$               enddo
C$$$            endif
C$$$            CALL INTERPOLATE(NR,X,RWORK1,RWORK1P,RWORK1PP,NRR,XX,
C$$$     $           RWORK1PPP)
C$$$            if(IR.eq.IRP) then
C$$$               if(RWORK1PPP(ie).gt.0d0) RWORK1PPP(ie)=0d0
C$$$            endif
C$$$            GG(IRP,IR,:)=DCMPLX(0d0,RWORK1PPP(:))
C$$$            GG(IR,IRP,:)=GG(IRP,IR,:)
C$$$         ENDIF
C$$$      ENDDO
C$$$      ENDDO
C$$$      DEALLOCATE(RWORK1,RWORK1P,RWORK1PP,RWORK1PPP)
C
C
C linear interpolation:
C
      ALLOCATE (RWORK1(NR),RWORK1P(NRR))
C
      RWORK1(1:NR) = 0D0
      RWORK1P(1:NRR) = 0D0
      DO IR = 1,NREP
         DO IRP = IR,NREP
            IF ( IZERO(IRP,IR).EQ.1 ) THEN
               RWORK1(1:NR) = DIMAG(GR(1:NR,IRP,IR))
               IF ( IR.EQ.IRP ) THEN
                  DO IE = 1,NR
                     IF ( RWORK1(IE).GT.0D0 ) RWORK1(IE) = 0D0
                  END DO
               END IF
               DO IE = 1,NRR
                  CALL DMFT_LINNTERP(NR,X,RWORK1,XX(IE),RWORK1P(IE))
                  IF ( IR.EQ.IRP ) THEN
                     IF ( RWORK1P(IE).GT.0D0 ) RWORK1P(IE) = 0D0
                  END IF
               END DO
               GG(IRP,IR,1:NRR) = DCMPLX(0D0,RWORK1P(1:NRR))
               GG(IR,IRP,1:NRR) = GG(IRP,IR,1:NRR)
            END IF
         END DO
      END DO
C
      DEALLOCATE (RWORK1,RWORK1P)
C
      DEALLOCATE (GR,X)
C
C normalization:
C$$$      ALLOCATE(RWORK1(NRR))
C$$$      DO IR=1,NREP
C$$$         RWORK1(:)=DIMAG(GG(IR,IR,:))
C$$$         CALL INTEGRAL(TMP,DX,RWORK1,NRR,1,NRR)
C$$$         GG(IR,IR,:)=GG(IR,IR,:)/DABS(TMP)*PI
C$$$      ENDDO
C$$$      DEALLOCATE(RWORK1)
C
C check printout:
      IF ( IPRT ) THEN
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_GF_PADE_INTPOL_IM.DAT')
         DO IE = 1,NRR
            WRITE (1,'(1000E15.5)') XX(IE),
     &                              ((DIMAG(GG(IR,IRP,IE)),IRP=IR,IR),
     &                              IR=1,NREP)
         END DO
         CLOSE (1)
      END IF
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C I'm not shure if it is necessary, because this DESTROYS THE
C ROTATIONAL INVARIANCE OF SIGMA, but it was given by Drchal (DJK05)
      DO IR1 = 1,NREP
         DO IR2 = 1,NREP
            DO IR3 = 1,NREP
               DO IR4 = 1,NREP
                  TMP = 0D0
                  IF ( (IORB(IR1,NREP).EQ.IORB(IR3,NREP) .AND. 
     &                 IORB(IR2,NREP).EQ.IORB(IR4,NREP)) .OR. 
     &                 (IR1.EQ.IR4 .AND. IR2.EQ.IR3) ) TMP = 1D0
                  V(IND(IR1,IR2),IND(IR3,IR4))
     &               = V(IND(IR1,IR2),IND(IR3,IR4))*TMP
               END DO
            END DO
         END DO
      END DO
C
C printout:
      IF ( IPRT ) CALL RMATSTRUCT(TXT_T(1:LTXT_T)//' 1/2 COULOMB MATRIX'
     &                            ,V,NREP**2/2,NREP**2,0,0,0,1D-8,6)
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Calling DMFT-TMA-solver
C============================================================
      ALLOCATE (SE_AR(NRR,NREP,NREP))
      IF ( .NOT.ALLOCATED(SE_AR_SAV) )
     &     ALLOCATE (SE_AR_SAV(NRR,NREP,NREP))
      XX(1:NRR) = XX(1:NRR) - EFERMIOLD
      SE_AR(1:NRR,1:NREP,1:NREP) = C0
C
      IF ( DMFT_FIX_DYN_SE ) THEN
         XX = XX_SAV
         SE_AR = SE_AR_SAV
      ELSE IF ( .NOT.TMA_STATIC ) THEN
         CALL DMFT_SPTTMA(NREP,NRR,XX,GG,IZERO,SE_AR,V,IPRT,TOL)
         DO IE = 1,NRR
            SE_AR(IE,1:NREP,1:NREP) = SE_AR(IE,1:NREP,1:NREP)
     &                                - SE_AR((NRR/2+1),1:NREP,1:NREP)
         END DO
         XX_SAV = XX
         SE_AR_SAV = SE_AR
      END IF
C
      DEALLOCATE (GG)
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Printout of Im-Sigma on the Re-axis:
C============================================================
      IF ( IPRT ) THEN
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_SE11_IM.DAT')
         OPEN (2,FILE=TXT_T(1:LTXT_T)//'_SE22_IM.DAT')
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            OPEN(3,FILE=TXT_T(1:LTXT_T)//'_SE12_IM.DAT')
C$$$            OPEN(4,FILE=TXT_T(1:LTXT_T)//'_SE21_IM.DAT')
C$$$         ENDIF
         DO IE = 1,NRR
            IF ( DABS(XX(IE)-XX(NRR/2+1)).LE.2D0 ) THEN
C11
               WRITE (1,'(1000E15.5)') XX(IE),
     &                                 ((DIMAG(SE_AR(IE,IR,IRP)),IRP=IR,
     &                                 IR),IR=1,NREP/2)
C22
               WRITE (2,'(1000E15.5)') XX(IE),
     &                                 ((DIMAG(SE_AR(IE,IR,IRP)),IRP=IR,
     &                                 IR),IR=NREP/2+1,NREP)
C$$$               IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$C12
C$$$                  WRITE(3,'(1000E15.5)')
C$$$     $                 XX(IE), ((DIMAG(SE_AR(IE,IR,IRP)),
C$$$     $                 IRP=1,NREP/2),IR=NREP/2+1,NREP)
C$$$C21
C$$$                  WRITE(4,'(1000E15.5)')
C$$$     $                 XX(IE), ((DIMAG(SE_AR(IE,IR,IRP)),
C$$$     $                 IRP=NREP/2+1,NREP),IR=1,NREP/2)
C$$$               ENDIF
            END IF
         END DO
         CLOSE (1)
         CLOSE (2)
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            CLOSE(3)
C$$$            CLOSE(4)
C$$$         ENDIF
      END IF
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Calculate Re-Sigma on the Re-axis:
C============================================================
C
      ALLOCATE (CWORK1(NRR),CWORK1P(NRR),SEZNEW(NREP,NREP,NRR),
     &          CWORK2(NREP,NREP))
C
      CWORK1P(1:NRR) = DCMPLX(XX(1:NRR),TOL)
C
      SEZNEW(1:NREP,1:NREP,1:NRR) = C0
      CWORK1(1:NRR) = C0
      CWORK2(1:NREP,1:NREP) = C0
      DO IR = 1,NREP
         DO IRP = IR,NREP
            IF ( IZERO(IRP,IR).EQ.1 ) THEN
               CALL DMFT_CAUCHYT(NRR,XX,SE_AR(1,IRP,IR),NRR,CWORK1P,
     &                           CWORK1)
               SEZNEW(IRP,IR,1:NRR) = CWORK1(1:NRR)
               SEZNEW(IR,IRP,1:NRR) = CWORK1(1:NRR)
            END IF
         END DO
      END DO
C
      CWORK2(1:NREP,1:NREP) = SEZNEW(1:NREP,1:NREP,NRR/2+1)
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Print out Re-Sigma on the Re-axis:
C============================================================
      IF ( IPRT ) THEN
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_SE11KK_RE.DAT')
                                                    ! up-up SIGMA
         OPEN (2,FILE=TXT_T(1:LTXT_T)//'_SE22KK_RE.DAT')
                                                    ! dn-dn SIGMA
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            OPEN(3,FILE=TXT_T(1:LTXT_T)//'_SE12_RE.DAT') ! up-dn SIGMA
C$$$            OPEN(4,FILE=TXT_T(1:LTXT_T)//'_SE21_RE.DAT') ! dn-up SIGMA
C$$$         ENDIF
         DO IE = 1,NRR
            IF ( DABS(XX(IE)-XX(NRR/2+1)).LE.2D0 ) THEN
C11
               WRITE (1,'(1000E15.5)') XX(IE),
     &                                 ((DREAL(SEZNEW(IR,IRP,IE)),
     &                                 IRP=IR,IR),IR=1,NREP/2)
C22
               WRITE (2,'(1000E15.5)') XX(IE),
     &                                 ((DREAL(SEZNEW(IR,IRP,IE)),
     &                                 IRP=IR,IR),IR=NREP/2+1,NREP)
C$$$               IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$C12
C$$$                  WRITE(3,'(1000E15.5)')
C$$$     $                 XX(IE), ((DREAL(SEZNEW(IR,IRP,IE)),
C$$$     $                 IRP=IR,IR),IR=NREP/2+1,NREP)
C$$$C21
C$$$                  WRITE(4,'(1000E15.5)')
C$$$     $                 XX(IE), ((DREAL(SEZNEW(IR,IRP,IE)),
C$$$     $                 IRP=IR,IR),IR=1,NREP/2)
C$$$               ENDIF
            END IF
         END DO
         CLOSE (1)
         CLOSE (2)
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_SE11KK_IM.DAT')
         OPEN (2,FILE=TXT_T(1:LTXT_T)//'_SE22KK_IM.DAT')
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            OPEN(3,FILE=TXT_T(1:LTXT_T)//'_SE12_IM.DAT')
C$$$            OPEN(4,FILE=TXT_T(1:LTXT_T)//'_SE21_IM.DAT')
C$$$         ENDIF
         DO IE = 1,NRR
            IF ( DABS(XX(IE)-XX(NRR/2+1)).LE.2D0 ) THEN
C11
               WRITE (1,'(1000E15.5)') XX(IE),
     &                                 ((DIMAG(SEZNEW(IR,IRP,IE)),
     &                                 IRP=IR,IR),IR=1,NREP/2)
C22
               WRITE (2,'(1000E15.5)') XX(IE),
     &                                 ((DIMAG(SEZNEW(IR,IRP,IE)),
     &                                 IRP=IR,IR),IR=NREP/2+1,NREP)
C$$$               IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$C12
C$$$                  WRITE(3,'(1000E15.5)')
C$$$     $                 XX(IE), ((DIMAG(SE_AR(IE,IR,IRP)),
C$$$     $                 IRP=1,NREP/2),IR=NREP/2+1,NREP)
C$$$C21
C$$$                  WRITE(4,'(1000E15.5)')
C$$$     $                 XX(IE), ((DIMAG(SE_AR(IE,IR,IRP)),
C$$$     $                 IRP=NREP/2+1,NREP),IR=1,NREP/2)
C$$$               ENDIF
            END IF
         END DO
         CLOSE (1)
         CLOSE (2)
C
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            CLOSE(3)
C$$$            CLOSE(4)
C$$$         ENDIF
      END IF ! IPRT
C
C
      DEALLOCATE (CWORK1,CWORK1P,SEZNEW,CWORK2)
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Return the XX mesh to the initial reference:
C============================================================
      XX(1:NRR) = XX(1:NRR) + EFERMINEW
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C Extrapolate Sigma on to semicircle:
C============================================================
C
C Check zero elements:
      IZERO(1:NREP,1:NREP) = 0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            CTMP = C0
            DO IE = 1,NRR
               CTMP = CTMP + CDABS(SE_AR(IE,IRP,IR))
            END DO
            IF ( CDABS(CTMP).GE.TOL ) IZERO(IRP,IR) = 1
         END DO
      END DO
C
C
C Make KKT of Im-Sigma:
      ALLOCATE (SEZNEW(NREP,NREP,NEGF))
      ALLOCATE (CWORK1(NEGF))
      SEZNEW(1:NREP,1:NREP,1:NEGF) = C0
      CWORK1(1:NEGF) = C0
      DO IR = 1,NREP
         DO IRP = IR,NREP
            IF ( IZERO(IRP,IR).EQ.1 ) THEN
C
               CALL DMFT_CAUCHYT(NRR,XX,SE_AR(1,IRP,IR),NEGF,EGFNEW,
     &                           CWORK1)
C
               SEZNEW(IRP,IR,1:NEGF) = CWORK1(1:NEGF)
               SEZNEW(IR,IRP,1:NEGF) = CWORK1(1:NEGF)
            END IF
         END DO
      END DO
      DEALLOCATE (CWORK1)
C
C
C
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C  Make transformations on the semicircle: DOUBLE COUNTING
C============================================================
C
      ALLOCATE (CWORK2(NREP,NREP))
C
C Sigma(E) = Sigma(E)-Sigma(EF) + LDAU:
      IF ( DMFTDBLC.EQ.'LDAU' ) THEN
         CWORK2(1:NREP,1:NREP) = SEZNEW(1:NREP,1:NREP,NEGF)
         DO IE = 1,NEGF
            SEZNEW(1:NREP,1:NREP,IE) = SEZNEW(1:NREP,1:NREP,IE)
     &                                 - DREAL(CWORK2(1:NREP,1:NREP))
         END DO
         SEZNEW(1:NREP,1:NREP,1:NEGF) = DCMPLX(TMA_REFACTOR*DREAL(SEZNEW
     &                                  (1:NREP,1:NREP,1:NEGF)),
     &                                  DIMAG(SEZNEW(1:NREP,1:NREP,
     &                                  1:NEGF)))
         DO IR = 1,NREP
            DO IRP = 1,NREP
               SEZNEW(IR,IRP,1:NEGF) = SEZNEW(IR,IRP,1:NEGF)
     &                                 + SIGSTAT(IR,IRP)
            END DO
         END DO
      ELSE
C Sigma(EF)=0:
         CWORK2(1:NREP,1:NREP) = SEZNEW(1:NREP,1:NREP,NEGF)
         DO IE = 1,NEGF
            SEZNEW(1:NREP,1:NREP,IE) = SEZNEW(1:NREP,1:NREP,IE)
     &                                 - CWORK2(1:NREP,1:NREP)
         END DO
         SEZNEW(1:NREP,1:NREP,1:NEGF) = DCMPLX(TMA_REFACTOR*DREAL(SEZNEW
     &                                  (1:NREP,1:NREP,1:NEGF)),
     &                                  DIMAG(SEZNEW(1:NREP,1:NREP,
     &                                  1:NEGF)))
         SEZNEW(1:NREP,1:NREP,NEGF) = C0
      END IF                    ! DMFTDBLC
C
C
      DEALLOCATE (CWORK2)
C
C Check zero elements:
      IZERO(1:NREP,1:NREP) = 0
      DO IR = 1,NREP
         DO IRP = 1,NREP
            CTMP = C0
            DO IE = 1,NEGF
               CTMP = CTMP + CDABS(SEZNEW(IRP,IR,IE))
            END DO
            IF ( CDABS(CTMP).GE.TOL ) IZERO(IRP,IR) = 1
         END DO
      END DO
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C write SIGMA in the unformatted file (on the semicircle):
C============================================================
      FILNAM = TXT_T(1:LTXT_T)//'_SIGUNFORM'//'.sig'
      LFILNAM = LTXT_T + 10 + 4
C
      XX(1:NRR) = XX(1:NRR) - EFERMINEW
C
      OPEN (1,FILE=FILNAM(1:LFILNAM),FORM='unformatted')
      WRITE (1) NEGF,NREP,TMA_STATIC
      WRITE (1) EGFNEW,SEZNEW,IZERO
      WRITE (1) EREFDMFT
      WRITE (1) WRSIGKK
      WRITE (1) SIGSTAT
      WRITE (1) NRR
      WRITE (1) XX
      WRITE (1) SE_AR
      CLOSE (1)
C
      DEALLOCATE (IZERO)
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C printout SIGMA on the semicircle:
C============================================================
      IF ( IPRT ) THEN
C
         OPEN (1,FILE=TXT_T(1:LTXT_T)//'_SEZNEW11_IM.DAT')
         OPEN (2,FILE=TXT_T(1:LTXT_T)//'_SEZNEW11_RE.DAT')
         OPEN (3,FILE=TXT_T(1:LTXT_T)//'_SEZNEW22_IM.DAT')
         OPEN (4,FILE=TXT_T(1:LTXT_T)//'_SEZNEW22_RE.DAT')
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            OPEN(5,FILE=TXT_T(1:LTXT_T)//'_SEZNEW12_IM.DAT')
C$$$            OPEN(6,FILE=TXT_T(1:LTXT_T)//'_SEZNEW12_RE.DAT')
C$$$            OPEN(7,FILE=TXT_T(1:LTXT_T)//'_SEZNEW21_IM.DAT')
C$$$            OPEN(8,FILE=TXT_T(1:LTXT_T)//'_SEZNEW21_RE.DAT')
C$$$         ENDIF
C
         DO IE = 1,NEGF
C11
            WRITE (1,'(1000E15.5)') DREAL(EGFNEW(IE)),
     &                              ((DIMAG(SEZNEW(IR,IRP,IE)),IRP=IR,
     &                              IR),IR=1,NREP/2)
            WRITE (2,'(1000E15.5)') DREAL(EGFNEW(IE)),
     &                              ((DREAL(SEZNEW(IR,IRP,IE)),IRP=IR,
     &                              IR),IR=1,NREP/2)
C22
            WRITE (3,'(1000E15.5)') DREAL(EGFNEW(IE)),
     &                              ((DIMAG(SEZNEW(IR,IRP,IE)),IRP=IR,
     &                              IR),IR=NREP/2+1,NREP)
            WRITE (4,'(1000E15.5)') DREAL(EGFNEW(IE)),
     &                              ((DREAL(SEZNEW(IR,IRP,IE)),IRP=IR,
     &                              IR),IR=NREP/2+1,NREP)
C
C$$$            IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$C12
C$$$               WRITE(5,'(1000E15.5)') DREAL(EGFNEW(IE)),
C$$$     $              ((DIMAG(SEZNEW(IR,IRP,IE)), IRP=IR,
C$$$     $              IR),IR=NREP/2+1,NREP)
C$$$               WRITE(6,'(1000E15.5)') DREAL(EGFNEW(IE)),
C$$$     $              ((DREAL(SEZNEW(IR,IRP,IE)),
C$$$     $              IRP=IR,IR),IR=NREP/2+1,NREP)
C$$$C21
C$$$               WRITE(7,'(1000E15.5)') DREAL(EGFNEW(IE)),
C$$$     $              ((DIMAG(SEZNEW(IR,IRP,IE)),
C$$$     $              IRP=IR,IR),IR=1,NREP/2)
C$$$               WRITE(8,'(1000E15.5)') DREAL(EGFNEW(IE)),
C$$$     $              ((DREAL(SEZNEW(IR,IRP,IE)),
C$$$     $              IRP=IR,IR),IR=1,NREP/2)
C$$$            ENDIF
C
         END DO
C
         CLOSE (1)
         CLOSE (2)
         CLOSE (3)
         CLOSE (4)
C$$$         IF(.NOT.TMA_SPINFLIPOFF) THEN
C$$$            CLOSE(5)
C$$$            CLOSE(6)
C$$$            CLOSE(7)
C$$$            CLOSE(8)
C$$$         ENDIF
      END IF
C============================================================
C
C
C                   *    *    *
C
C
C============================================================
C     Mix new self-energy SEZNEW with the old SEZ and
C     store the result in DMFTSIGMA
C============================================================
C
C
 100  CONTINUE
      IF ( ITRSCF.EQ.1 ) THEN
C
         DMFTSIGMA(I1:I2,I1:I2,1:NEGF) = SEZNEW(1:NLM,1:NLM,1:NEGF)
                                         !1_1
         DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
     &      = SEZNEW(NLM+1:NREP,NLM+1:NREP,1:NEGF) !2_2
         DMFTSIGMA(I1:I2,I3:I4,1:NEGF) = SEZNEW(1:NLM,NLM+1:NREP,1:NEGF)
                                              !1_2
         DMFTSIGMA(I3:I4,I1:I2,1:NEGF) = SEZNEW(NLM+1:NREP,1:NLM,1:NEGF)
                                              !2_1
C
         SEZ(1:NREP,1:NREP,1:NEGF,IT) = SEZNEW(1:NREP,1:NREP,1:NEGF)
C
      ELSE
C
         DMFTSIGMA(I1:I2,I1:I2,1:NEGF) = SEZ(1:NLM,1:NLM,1:NEGF,IT)
     &      *(1D0-DMFTMIX) + SEZNEW(1:NLM,1:NLM,1:NEGF)*DMFTMIX
                                                               !1_1
         DMFTSIGMA(I3:I4,I3:I4,1:NEGF)
     &      = SEZ(NLM+1:NREP,NLM+1:NREP,1:NEGF,IT)*(1D0-DMFTMIX)
     &      + SEZNEW(NLM+1:NREP,NLM+1:NREP,1:NEGF)*DMFTMIX !2_2
         DMFTSIGMA(I1:I2,I3:I4,:) = SEZ(1:NLM,NLM+1:NREP,1:NEGF,IT)
     &                              *(1D0-DMFTMIX)
     &                              + SEZNEW(1:NLM,NLM+1:NREP,1:NEGF)
     &                              *DMFTMIX                        !1_2
         DMFTSIGMA(I3:I4,I1:I2,1:NEGF) = SEZ(NLM+1:NREP,1:NLM,1:NEGF,IT)
     &      *(1D0-DMFTMIX) + SEZNEW(NLM+1:NREP,1:NLM,1:NEGF)*DMFTMIX
                                                      !2_1
C
         SEZ(1:NREP,1:NREP,1:NEGF,IT) = SEZ(1:NREP,1:NREP,1:NEGF,IT)
     &                                  *(1D0-DMFTMIX)
     &                                  + SEZNEW(1:NREP,1:NREP,1:NEGF)
     &                                  *DMFTMIX
C
      END IF
      DEALLOCATE (SEZNEW)
C
      ICALL = ICALL + 1
C
      DEALLOCATE (XX,SE_AR)
      DEALLOCATE (V,IND)
      END
C*==ispn.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      FUNCTION ISPN(I,NREP)
      IMPLICIT NONE
C*--ISPN1143
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,NREP
      INTEGER ISPN
C
C*** End of declarations rewritten by SPAG
C
      IF ( I.LE.NREP/2 ) THEN
         ISPN = 1
      ELSE
         ISPN = 2
      END IF
      END
C*==iorb.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
C
      FUNCTION IORB(I,NREP)
      IMPLICIT NONE
C*--IORB1172
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,NREP
      INTEGER IORB
C
C*** End of declarations rewritten by SPAG
C
      IF ( I.LE.NREP/2 ) THEN
         IORB = -NREP/4 - 1 + I
      ELSE
         IORB = -NREP/4 - NREP/2 - 1 + I
      END IF
      END
C
