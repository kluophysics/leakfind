C*==dmft_u_matrix.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_U_MATRIX(EREFLDAU,IPRINT,JEFF,NFILLDAU,UEFF,ULDAU)
C   ********************************************************************
C   * Calculations of U-matrix                                         *
C   *  this is partially taken from FPSCFLDAUINIT                      *
C   *  Radial Slater integrals are calculated                          *
C   *      up to critical/rws radius (FP/ASA)                          *
C   *      and are renormalised w.r.t. input value of U and J          *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,PI,RY_EV
      USE MOD_ANGMOM,ONLY:NKMQ,NSPIN,NKM,NKMMAX,NLMAX,NMEMAX,NCPLWF
      USE MOD_RMESH,ONLY:JRWS,JRCRI,FULLPOT,R2DRDI,NRMAX,R
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NCPLWFMAX,IMT,NTMAX,TXT_T,LOPT,
     &    IKMCPLWF
      USE MOD_SITES,ONLY:IQAT
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_FILES,ONLY:IDUMMY
      IMPLICIT NONE
C*--DMFT_U_MATRIX20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DMFT_U_MATRIX')
      REAL*8 TOL
      PARAMETER (TOL=1.D-10)
C
C Dummy arguments
C
      INTEGER IPRINT,NFILLDAU
      COMPLEX*16 EREFLDAU(NTMAX)
      REAL*8 JEFF(NTMAX),UEFF(NTMAX),
     &       ULDAU(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX,NTMAX)
C
C Local variables
C
      REAL*8 A(:,:,:,:,:),FCLMB(:),G12,G34,QNORM,RAUX,RHO(:),RHORSQ,
     &       RLOP,RPW(:,:),SCL,SG(:),SL(:),SUMFCLMB,TG(:),TL(:),W3J,
     &       WGTFCLMB(:)
      REAL*8 CGC_RACAH,GAUNT_CYLM
      COMPLEX*16 DZJ(:,:,:,:),DZZ(:,:,:,:),DZZBA(1,1,1,1,1),JF(:,:,:),
     &           JG(:,:,:),MEZJ(:,:,:,:),MEZZ(:,:,:,:),TMP,X,ZF(:,:,:),
     &           ZG(:,:,:)
      INTEGER I,IA_ERR,IKM,IKMBOT,IKMTOP,IL,IM,IM1,IM2,IM3,IM4,IQ,IR,
     &        IRTOP,IS,IT,J,JRTOP,K,L,LF,LFMAX,M,M1,M2,M3,M4,NWF,N
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      DATA IA_ERR/0/
C
      ALLOCATABLE WGTFCLMB,A,SG,TG,SL,TL,RHO,RPW
      ALLOCATABLE FCLMB,MEZZ,MEZJ
      ALLOCATABLE DZZ,DZJ,ZF,ZG,JF,JG
C
      ALLOCATE (DZZ(NKMMAX,NKMMAX,NTMAX,1))
      ALLOCATE (DZJ(NKMMAX,NKMMAX,NTMAX,1))
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      ALLOCATE (WGTFCLMB(0:2*NLMAX))
      ALLOCATE (A(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX,0:2*NLMAX))
      ALLOCATE (SG(NRMAX),TG(NRMAX),SL(NRMAX),TL(NRMAX))
      ALLOCATE (RHO(NRMAX),RPW(NRMAX,2*NLMAX))
      ALLOCATE (FCLMB(0:2*NLMAX))
      ALLOCATE (MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
C
      MEZZ = C0
      MEZJ = C0
      DZZ = C0
      DZJ = C0
      RPW = 0.0D0
C
      IF ( IPRINT.GT.0 ) WRITE (6,99001)
C
C     ------------------------------------- calculate overlap matrix for
C                                                  renormalisation of WF
      IF ( IREL.EQ.3 ) THEN
C
         CALL ME_OBS_REL(NFILLDAU,NFILLDAU,.TRUE.,DZZ,DZJ,DZZBA,1)
C
      ELSE
C
         CALL ME_OBS_SRA(NFILLDAU,MEZZ,MEZJ,DZZ,DZJ)
C
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IQ = IQAT(1,IT)
C
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C
C ---------------------------- read, normalize and rewrite wavefunctions
C
         IF ( IREL.EQ.3 ) CALL WAVFUN_READ_REL(NFILLDAU,IT,0,ZG,ZF,JG,
     &        JF,IRTOP,NCPLWF,IKMCPLWF)
C
         DO IKM = 1,NKMQ(IQ)
C
            IF ( IREL.NE.3 ) THEN
               READ (NFILLDAU,REC=IKM+(IT-1)*NKM) IDUMMY,STRP,IDUMMY,
     &               NCPLWF(IKM),(IKMCPLWF(J,IKM),J=1,NCPLWF(IKM)),
     &               ((ZG(I,K,IKM),I=1,IRTOP),K=1,NCPLWF(IKM))
               ZF(:,:,IKM) = C0
            END IF
C

            N= NCPLWF(IKM)
            IF ( ABS(DZZ(IKM,IKM,IT,1)).GT.TOL ) THEN
               X = 1D0/SQRT(DZZ(IKM,IKM,IT,1))
              ZG(1:IRTOP,1:N,IKM) = X*ZG(1:IRTOP,1:N,IKM)
              ZF(1:IRTOP,1:N,IKM) = X*ZF(1:IRTOP,1:N,IKM)
            ELSE
               ZG(1:IRTOP,1:N,IKM) = C0
               ZF(1:IRTOP,1:N,IKM) = C0
            END IF

            IF ( ABS(DZZ(IKM,IKM,IT,1)).GT.TOL ) THEN
               X = 1D0/SQRT(DZZ(IKM,IKM,IT,1))
            ELSE
               X = C0
            END IF
C
            DO K = 1,NCPLWFMAX
               DO I = 1,IRTOP
                  ZG(I,K,IKM) = X*ZG(I,K,IKM)
                  ZF(I,K,IKM) = X*ZF(I,K,IKM)
               END DO
            END DO
         END DO
C
         IF ( LOPT(IT).GE.0 ) THEN
C
            IF ( LOPT(IT).LT.2 .OR. LOPT(IT).GT.3 ) THEN
               WRITE (6,*) ' LOP = ',LOPT(IT),'  for IT = ',IT,TXT_T(IT)
               CALL STOP_MESSAGE(ROUTINE,'LOPT < 2 or LOPT > 3')
            END IF
C
            IF ( IPRINT.GT.0 ) WRITE (6,99004) IT,TXT_T(IT)
C=======================================================================
C     calculate Slater integrals
C=======================================================================
C
C--------------------------------------------Setup upper limit for radial
C                                integral (FULLPOT=critical rad; ASA=RWS)
C
C
            IF ( FULLPOT ) THEN
               JRTOP = JRCRI(IM)
            ELSE
               JRTOP = JRWS(IM)
            END IF
C
            LFMAX = 2*LOPT(IT)
C
            DO IR = 1,IRTOP
               RPW(IR,1) = R(IR,IM)
               DO IL = 2,2*NLMAX
                  RPW(IR,IL) = RPW(IR,IL-1)*R(IR,IM)
               END DO
            END DO
C
            RHO(1:NRMAX) = 0.0D0
            SG(1:NRMAX) = 0.0D0
            TG(1:NRMAX) = 0.0D0
            SL(1:NRMAX) = 0.0D0
            TL(1:NRMAX) = 0.0D0
C
            DO IS = 1,NSPIN
C
               IF ( IREL.EQ.3 ) THEN
                  IKMBOT = 2*LOPT(IT)**2 + 1
                  IKMTOP = IKMBOT - 1 + 2*(2*LOPT(IT)+1)
               ELSE
                  IKMBOT = (NKMQ(IQ)/NSPIN)*(IS-1) + LOPT(IT)**2 + 1
                  IKMTOP = IKMBOT + LOPT(IT)*2
               END IF
C
               IF ( FULLPOT ) THEN
                  NWF = 0
                  DO IKM = IKMBOT,IKMTOP
                     DO J = 1,NCPLWFMAX
                        IF ( IKMCPLWF(J,IKM).EQ.IKM ) THEN
                           NWF = NWF + 1
                           DO IR = 1,IRTOP
                              RHO(IR) = RHO(IR) + DREAL(ZG(IR,J,IKM))**2
                           END DO
                        END IF
                     END DO
                  END DO
C
               ELSE
C
                  DO IKM = IKMBOT,IKMTOP
                     DO J = 1,NCPLWF(IKM)
                        DO IR = 1,IRTOP
                           RHO(IR) = RHO(IR) + DREAL(ZG(IR,J,IKM))**2
                        END DO
                     END DO
                  END DO
               END IF
            END DO
C
            DO IR = 1,IRTOP
               TL(IR) = RHO(IR)*R2DRDI(IR,IM)
            END DO
C
            CALL RRADINT_R(IM,TL,SL)
            QNORM = SL(JRTOP)
C
            DO IR = 1,IRTOP
               RHO(IR) = RHO(IR)/QNORM
            END DO
C
            RLOP = DBLE(LOPT(IT))
            SUMFCLMB = 0D0
C -------------------------------------------------- calculate potential
            DO LF = 2,LFMAX,2
C -------------------------------------- evaluate first radial integrals
               DO IR = 1,IRTOP
                  RHORSQ = 2.0D0*RHO(IR)*R2DRDI(IR,IM)
                  TL(IR) = RHORSQ*RPW(IR,LF)
                  TG(IR) = RHORSQ/RPW(IR,LF+1)
               END DO
C
               CALL RRADINT_R(IM,TL,SL)
               CALL RRADINT_R(IM,TG,SG)
C
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
               DO IR = 1,IRTOP
                  SL(IR) = SL(IR)/RPW(IR,LF+1) + (SG(IRTOP)-SG(IR))
     &                     *RPW(IR,LF)
               END DO
C
C ------------------------------------- calculate Coulomb integrals F^LF
               DO IR = 1,IRTOP
                  SG(IR) = RHO(IR)*R2DRDI(IR,IM)*SL(IR)
               END DO
C
               CALL RRADINT_R(IM,SG,SL)
               FCLMB(LF) = SL(JRTOP)
C
               W3J = (-1)**NINT(RLOP)*(1D0/SQRT(2D0*RLOP+1D0))
     &               *CGC_RACAH(RLOP,DBLE(LF),RLOP,0D0,0D0,0D0)
C
               WGTFCLMB(LF) = ((2*RLOP+1)/(2*RLOP))*W3J**2
               SUMFCLMB = SUMFCLMB + WGTFCLMB(LF)*FCLMB(LF)
            END DO
C
C
            SCL = JEFF(IT)/RY_EV/SUMFCLMB
C
            IF ( IPRINT.GT.0 ) WRITE (6,99002) DREAL(EREFLDAU(IT)),
     &                                UEFF(IT)/RY_EV,UEFF(IT),JEFF(IT)
     &                                /RY_EV,JEFF(IT),SUMFCLMB,
     &                                SUMFCLMB*RY_EV,(UEFF(IT)-JEFF(IT))
     &                                *0.5D0/RY_EV,(UEFF(IT)-JEFF(IT))
     &                                *0.5D0,SCL
            FCLMB(0) = UEFF(IT)/RY_EV
C
            DO LF = 2,LFMAX,2
               IF ( IPRINT.GT.0 ) WRITE (6,99003) LF,WGTFCLMB(LF),
     &              FCLMB(LF),FCLMB(LF)*RY_EV,SCL*FCLMB(LF),
     &              SCL*FCLMB(LF)*RY_EV
               FCLMB(LF) = SCL*FCLMB(LF)
            END DO
C=======================================================================
C                     calculate coefficient matrix
C=======================================================================
            L = LOPT(IT)
C
            CALL RINIT((2*NLMAX)**4*(1+2*NLMAX),A)
C
            DO LF = 0,LFMAX,2
               DO M1 = -L, + L
                  IM1 = L + M1 + 1
                  DO M2 = -L, + L
                     IM2 = L + M2 + 1
                     DO M3 = -L, + L
                        IM3 = L + M3 + 1
                        M4 = M1 - M2 + M3
                        IF ( -L.LE.M4 .AND. M4.LE.L ) THEN
                           IM4 = L + M4 + 1
                           RAUX = 0D0
                           DO K = -LF, + LF
                              G12 = GAUNT_CYLM(L,M1,LF,K,L,M2)
                              G34 = GAUNT_CYLM(L,M3,LF,-K,L,M4)
                              RAUX = RAUX + G12*G34*(-1)**ABS(K)
                           END DO
                           A(IM1,IM2,IM3,IM4,LF)
     &                        = RAUX*4*PI/(2D0*LF+1D0)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
C
            CALL RINIT((NLMAX*2)**4,ULDAU(1,1,1,1,IT))
C
            DO IM1 = 1,2*L + 1
               DO IM2 = 1,2*L + 1
                  DO IM3 = 1,2*L + 1
                     DO IM4 = 1,2*L + 1
                        DO LF = 0,LFMAX,2
                           ULDAU(IM1,IM2,IM3,IM4,IT)
     &                        = ULDAU(IM1,IM2,IM3,IM4,IT)
     &                        + A(IM1,IM2,IM3,IM4,LF)*FCLMB(LF)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
C
C-----------------------------------------------------------------------
C check the sum rules for UDMFT-matrix:
C-----------------------------------------------------------------------
C U=1/(2*l+1)^2*sum_{mm1}U(m,m1,m,m1)
            TMP = 0D0
            DO M = 1,2*L + 1
               DO M1 = 1,2*L + 1
                  TMP = TMP + ULDAU(M,M,M1,M1,IT)
               END DO
            END DO
            TMP = TMP/((2D0*L+1D0)**2)
            IF ( ABS(TMP-UEFF(IT)/RY_EV).GE.1D-10 ) THEN
               PRINT *,'Warning!!'
               WRITE (6,'(10X,A12,2F8.3)') 'sum rule 1: ',UEFF(IT)
     &                /RY_EV,TMP
            END IF
C     1/((2*l+1))*sum_{mm1}U(m,m1,m1,m)=(U+2*l*J):
            TMP = 0D0
            DO M = 1,2*L + 1
               DO M1 = 1,2*L + 1
                  TMP = TMP + ULDAU(M,M1,M1,M,IT)
               END DO
            END DO
            TMP = TMP/(2D0*L+1D0)
            IF ( ABS(TMP-(UEFF(IT)+2*L*JEFF(IT))/RY_EV).GE.1D-10 ) THEN
               PRINT *,'Warning!!'
               WRITE (6,'(10X,A12,2F8.3)') 'sum rule 2 :',
     &                (UEFF(IT)+2*L*JEFF(IT))/RY_EV,TMP
            END IF
C U=1/(2l*(2l+1))*sum_{mm1}[U(m,m1,m,m1)-U(m,m1,m1,m)]:
            TMP = 0D0
            DO M = 1,2*L + 1
               DO M1 = 1,2*L + 1
                  TMP = TMP + ULDAU(M,M,M1,M1,IT) - ULDAU(M,M1,M1,M,IT)
               END DO
            END DO
            TMP = TMP/(2D0*L+1D0)/(2D0*L)
            IF ( ABS(TMP-(UEFF(IT)-JEFF(IT))/RY_EV).GE.1D-10 ) THEN
               PRINT *,'Warning!!'
               WRITE (6,'(10X,A12,2F8.3)') 'sum rule 3:',
     &                (UEFF(IT)-JEFF(IT))/RY_EV,TMP
            END IF
C
C-----------------------------------------------------------------------
C
C=======================================================================
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DEALLOCATE (DZZ,DZJ,ZF,ZG,STAT=IA_ERR)
C=======================================================================
C
99001 FORMAT (//,1X,79('*'),/,32X,'< DMFT_U_MATRIX  >',/,1X,79('*'),/)
99002 FORMAT (5X,'E(ref)  =',F7.4,' Ry  ',/,5X,'U(inp)  =',F7.4,' Ry =',
     &        F7.3,' eV',/,5X,'J(inp)  =',F7.4,' Ry =',F7.3,' eV',/,5X,
     &        'J(calc) =',F7.4,' Ry =',F7.3,' eV',/,5X,'(U-J)/2 =',F7.4,
     &        ' Ry =',F7.3,' eV',//,5X,'scaling factor for F^n: ',F7.4,
     &        //,5X,'n   W(F^n)    F^n  calculated',16X,'F^n  scaled')
99003 FORMAT (5X,I1,F9.4,F10.4,' Ry =',F7.3,' eV   ',' ==>',F9.4,
     &        ' Ry =',F7.3,' eV')
99004 FORMAT (//,5X,'setting up the Coulomb matrix  U  for  IT=',I2,2X,
     &        A,/)
      END
