C*==sigoff.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGOFF(INTEG,MSST,TSST,MSSQ,TAUQ,SIGOFFQ,IFILWF,
     &                  NSPINPROJ,MEOFF_BARG_T,MEOFF_AHE_T,ME_SZ_T,
     &                  X_VFT,NVFO_Q,IVFT_VFOQ,NVFTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-diagonal term of the conductivity tensor    *
C   *                                                                  *
C   *      OFF-diagonal w.r.t. the polarisation                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,DATSET,LDATSET,IOTMP
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ,A_RSIG,NKMPMAX,WKM1,IPIVKM
      USE MOD_SITES,ONLY:NQ,NQMAX,IQAT,IQBOT_CHI,IQTOP_CHI,ITOQ,NOQ
      USE MOD_TYPES,ONLY:NTMAX,NT,ITBOT,ITTOP,CTL
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_THERMAL,ONLY:FMAT_FT,MFMAT_FT,NVIBFLU,NVIBRA,NFLUCT,
     &    UMAT_VT
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,CONSI2,SOTSI2,EESI2,
     &    IRESPONSE_SOT,IRESPONSE_EDELSTEIN,IRESPONSE_ADA_ADA
      IMPLICIT NONE
C*--SIGOFF24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGOFF')
      REAL*8 TOL
      PARAMETER (TOL=1D-10)
C
C Dummy arguments
C
      INTEGER IFILWF,NSPINPROJ,NVFTMAX
      LOGICAL INTEG
      INTEGER IVFT_VFOQ(NVFTMAX,NQMAX),NVFO_Q(NQMAX)
      COMPLEX*16 MEOFF_AHE_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           MEOFF_BARG_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           ME_SZ_T(NKMMAX,NKMMAX,3,3,NTMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           SIGOFFQ(3,3,NSPINPROJ,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
      REAL*8 X_VFT(NVFTMAX)
C
C Local variables
C
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,ICALL,IFLUCT,IFT,IO,IPOL1,IPOL2,IQ,ISPINPROJ,
     &        ISPR,IT,IVFO,IVFT,IVIBRA,IVT,J,JPOL1,JPOL2,M,MUE,N,NUE
      COMPLEX*16 IMTAU_VFT(:,:),IMTSS_VFT(:,:),ME(:,:,:,:),
     &           MET(:,:,:,:,:),ME_FT(:,:,:,:),ME_VFT(:,:,:,:),
     &           MSS_VFT(:,:),SIGTMP(3,3),SIGTMPSS_VFT(:,:,:),
     &           SIGTMP_VFT(:,:,:),TAU_VFT(:,:),TSS_FT(:,:),TSS_VFT(:,:)
     &           ,WRK(:,:,:,:)
      LOGICAL LPRSIG
      CHARACTER*1 MSTR,NSTR
      REAL*8 PRESI2,WZ2
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE IMTSS_VFT,IMTAU_VFT,ME,MET,ME_VFT,SIGTMP_VFT
      ALLOCATABLE SIGTMPSS_VFT
      ALLOCATABLE TAU_VFT,TSS_VFT,ME_FT,TSS_FT,MSS_VFT
      ALLOCATABLE WRK
C
      IF ( INTEG ) RETURN
C
      LPRSIG = .FALSE.
C      LPRSIG = .TRUE.
C
      ICALL = ICALL + 1
C
      ALLOCATE (WRK(NKMMAX,NKMMAX,3,3))
      ALLOCATE (ME_FT(NKMMAX,NKMMAX,3,3))
      ALLOCATE (TSS_FT(NKMMAX,NKMMAX))
      ALLOCATE (MSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (TAU_VFT(NKMMAX,NKMMAX),TSS_VFT(NKMMAX,NKMMAX))
      ALLOCATE (ME(NKMMAX,NKMMAX,-1:+1,-1:+1))
      ALLOCATE (MET(NKMMAX,NKMMAX,3,3,NT))
      ALLOCATE (ME_VFT(NKMMAX,NKMMAX,3,3))
      ALLOCATE (SIGTMP_VFT(3,3,NVFTMAX))
      ALLOCATE (SIGTMPSS_VFT(3,3,NVFTMAX))
      MET(1:NKMMAX,1:NKMMAX,1:3,1:3,1:NT) = C0
C
C=======================================================================
C                               Initialize
C=======================================================================
C
      WZ2 = DSQRT(2.0D0)
C
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( FULLPOT ) THEN
            WRITE (6,*) ' <SIGOFF> Warning'
            WRITE (6,*) ' <SIGOFF> does not work for FULLPOT'
            RETURN
         END IF
C
C --------------------------- allocate array for angular matrix elements
C
         ALLOCATE (A_RSIG(NKMPMAX,NKMPMAX,-1:+1,-1:+1),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 )
     &         CALL STOP_MESSAGE(ROUTINE,'ALLOC AME_RSIG_SPIN')
      END IF
C
C ------------------------------------ calculate angular matrix elements
      CALL AME_RSIG_SPIN(1)
C
C=======================================================================
C
      M = NKMMAX
C
      ALLOCATE (IMTAU_VFT(M,M),IMTSS_VFT(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc IMTAUT')
C
      DO ISPR = 1,1 !NSPR,3
         ISPINPROJ = LIST_ISPR(ISPR)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(ISPINPROJ)
C
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
            PRESI2 = SOTSI2
         ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
            PRESI2 = EESI2
         ELSE
            PRESI2 = 1D-8*CONSI2
         END IF
C
CSW!!!!!!assuming old list of operators!!!!!!
         IF ( ISPINPROJ.EQ.11 ) CALL AME_RSIG_SPIN(11)
CSW!!!!!!assuming old list of operators!!!!!!
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
         LOOP_IT:DO IT = ITBOT,ITTOP
C
            M = NKMMAX
            N = NKM
            IQ = IQAT(1,IT)
C
            WRITE (6,99001) IT
C
            CALL SIGOFFME(IT,IFILWF,ME,ISPINPROJ)
C
C------------ + SIGOFFME calculates matrix elements for r dot alpha
C------------ + velocity operator is  c*alpha --> include c
C
            ME = -CTL(IT,1)*ME
C
C
            DO J = 1,NKM
               DO I = 1,NKM
                  MET(I,J,1,2,IT) = CI*(ME(I,J,-1,+1)-ME(I,J,+1,-1))
                  MET(I,J,2,1,IT) = -MET(I,J,1,2,IT)
                  MET(I,J,1,3,IT) = (ME(I,J,-1,0)-ME(I,J,+1,0)-ME(I,J,0,
     &                              -1)+ME(I,J,0,+1))/WZ2
                  MET(I,J,3,1,IT) = -MET(I,J,1,3,IT)
                  MET(I,J,2,3,IT) = (CI/WZ2)
     &                              *(ME(I,J,-1,0)+ME(I,J,+1,0)-ME(I,J,
     &                              0,-1)-ME(I,J,0,+1))
                  MET(I,J,3,2,IT) = -MET(I,J,2,3,IT)
               END DO
            END DO
C
            IF ( IPRINT.GE.3 .AND. 1.EQ.0 ) THEN
               CALL CMATSTRUCT('MET12',MET(1,1,1,2,IT),NKM,NKMMAX,3,3,0,
     &                         1D-8,6)
               CALL CMATSTRUCT('MET21',MET(1,1,2,1,IT),NKM,NKMMAX,3,3,0,
     &                         1D-8,6)
            END IF
C
CSW!!!!!!assuming old list of operators!!!!!!
            IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
CSW!!!!!!assuming old list of operators!!!!!!
               DO J = 1,NKM
                  DO I = 1,NKM
C
                     MET(J,I,1:3,1:3,IT) = MET(J,I,1:3,1:3,IT)
     &                  + MEOFF_BARG_T(J,I,1:3,1:3,IT)
C
                  END DO
               END DO
            END IF
C
CSW!!!!!!assuming old list of operators!!!!!!
            IF ( ISPINPROJ.EQ.2 ) THEN
CSW!!!!!!assuming old list of operators!!!!!!
C
CDK
C    the term in meoff_ahe_T is switched on here, numerical tests have
C    shown that when muliplying MEs here we have a product of a
C    symmetric and antisymmteric matrix, the trace over such always is
C    zero,  so this term so far gives zero contribution
CDK
               MET(1:NKMMAX,1:NKMMAX,1:3,1:3,IT)
     &            = MEOFF_AHE_T(1:NKMMAX,1:NKMMAX,1:3,1:3,IT)
     &            + ME_SZ_T(1:NKMMAX,1:NKMMAX,1:3,1:3,IT)
C
               IF ( IPRINT.GE.3 ) THEN
                  CALL CMATSTRUCT('meoff_ahe_T',MET(1,1,1,2,IT),NKM,
     &                            NKMMAX,3,3,0,1D-8,6)
                  CALL CMATSTRUCT('meoff_ahe_T',MET(1,1,2,1,IT),NKM,
     &                            NKMMAX,3,3,0,1D-8,6)
               END IF
C
            END IF
C
         END DO LOOP_IT
C                                                                     IT
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
         M = NKMMAX
C
         WRITE (6,99002)
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         DO IQ = IQBOT_CHI,IQTOP_CHI
            WRITE (6,99003) IQ
C
            CALL CINIT(3*3,SIGTMP)
            CALL CINIT(3*3*NVFTMAX,SIGTMP_VFT)
            CALL CINIT(3*3*NVFTMAX,SIGTMPSS_VFT)
C
            LOOP_IO:DO IO = 1,NOQ(IQ)
C
               IT = ITOQ(IO,IQ)
C
C=======================================================================
C         thermal lattice vibrations and/or spin fluctuations
C=======================================================================
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
                     ME_FT(:,:,:,:) = MET(:,:,:,:,IT)
C
                  ELSE
C
                     CALL CMATTRANS(N,M,WKM1,FMAT_FT(1,1,IFT),.TRUE.,
     &                              TSST(1,1,IT),'SAS+',TSS_FT)
C
C--------------------------------------------------------------- ROT MET
C
                     DO IPOL1 = 1,3
                        DO IPOL2 = 1,3
                           CALL ROTATE(MET(1,1,IPOL1,IPOL2,IT),'L->G',
     &                                 WRK(1,1,IPOL1,IPOL2),N,
     &                                 FMAT_FT(1,1,IFT),M)
                        END DO
                     END DO
C
                     ME_FT(:,:,:,:) = C0
C
                     DO IPOL1 = 1,3
                        DO IPOL2 = 1,3
C
                           DO JPOL1 = 1,3
                              DO JPOL2 = 1,3
                                 ME_FT(1:N,1:N,IPOL1,IPOL2)
     &                              = ME_FT(1:N,1:N,IPOL1,IPOL2)
     &                              + MFMAT_FT(JPOL1,IPOL1,IFT)
     &                              *MFMAT_FT(JPOL2,IPOL2,IFT)
     &                              *WRK(1:N,1:N,JPOL1,JPOL2)
                              END DO
                           END DO
C
                        END DO
                     END DO
C--------------------------------------------------------------- ROT MET
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
                        ME_VFT(:,:,:,:) = ME_FT(:,:,:,:)
C
                     ELSE
C
                        CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),TSS_FT,
     &                     'UAUT',TSS_VFT)
C
C------------------------------------------------------------ SHIFTT MET
                        DO IPOL1 = 1,3
                           DO IPOL2 = 1,3
                              CALL CMAT_U_TRANS(UMAT_VT(1,1,IVT),
     &                           ME_FT(1,1,IPOL1,IPOL2),'UAUT',
     &                           ME_VFT(1,1,IPOL1,IPOL2))
                           END DO
                        END DO
C------------------------------------------------------------ SHIFTT MET
C
                     END IF
C
C-------------------------------------------- MSS_VFT: m_vft = 1 / t_vft
C
                     IF ( NVIBFLU.EQ.1 ) THEN
                        MSS_VFT(:,:) = MSST(:,:,IT)
                     ELSE
                        CALL CMATINV3(N,M,IPIVKM,TSS_VFT,WKM1,MSS_VFT)
                     END IF
C
C-----------------------------------------------------------------------
C
                     TSS_VFT(1:N,1:N) = TSST(1:N,1:N,IT)
C
                     CALL GET_TAUT(MSS_VFT,MSSQ(1,1,IQ),TAUQ(1,1,IQ),
     &                             TAU_VFT)
C
C-----------------------------------------------------------------------
C
                     N = NKMQ(IQ)
C
                     DO J = 1,N
                        DO I = 1,N
                           IMTAU_VFT(I,J)
     &                        = (TAU_VFT(I,J)-DCONJG(TAU_VFT(J,I)))
     &                        /(2D0*CI)
C
                           IMTSS_VFT(I,J)
     &                        = (TSS_VFT(I,J)-DCONJG(TSS_VFT(J,I)))
     &                        /(2D0*CI)
                        END DO
                     END DO
C
                     DO MUE = 1,3
                        DO NUE = 1,3
C
C-----------------------------------------------------------------------
C
                           CALL CMATMUL(N,M,IMTAU_VFT,
     &                                  ME_VFT(1,1,MUE,NUE),WKM1)
C
                           DO I = 1,NKM
                              SIGTMP(MUE,NUE) = SIGTMP(MUE,NUE)
     &                           + X_VFT(IVFT)*WKM1(I,I)
C
                              SIGTMP_VFT(MUE,NUE,IVFT)
     &                           = SIGTMP_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*WKM1(I,I)
                           END DO
C
C--------------------------------  get single site t-matrix contribution
C
                           CALL CMATMUL(N,M,IMTSS_VFT,
     &                                  ME_VFT(1,1,MUE,NUE),WKM1)
C
                           DO I = 1,NKM
                              SIGTMPSS_VFT(MUE,NUE,IVFT)
     &                           = SIGTMPSS_VFT(MUE,NUE,IVFT)
     &                           + X_VFT(IVFT)*WKM1(I,I)
                           END DO
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
C ----------------------------------- suppress small elements and sum up
C
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(SIGTMP(MUE,NUE)).LT.TOL ) SIGTMP(MUE,NUE)
     &                 = C0
C
                  IF ( DIMAG(SIGTMP(MUE,NUE)).LT.TOL ) SIGTMP(MUE,NUE)
     &                 = DREAL(SIGTMP(MUE,NUE))
C
                  SIGOFFQ(MUE,NUE,ISPINPROJ,IQ) = SIGTMP(MUE,NUE)
C
               END DO
            END DO
C
C -------------------------------- write out conductivity contriibutions
C
            IF ( IPRINT.GE.3 ) THEN
               WRITE (6,99004)
               DO MUE = 1,3
                  WRITE (6,99007) (PRESI2*SIGTMP(MUE,NUE),NUE=1,3)
               END DO
            END IF
C
            IF ( IPRINT.GT.4 ) THEN
C
               WRITE (6,'(/,10x,"Type resolved (SI - units)")')
C
               DO IVFO = 1,NVFO_Q(IQ)
                  IVFT = IVFT_VFOQ(IVFO,IQ)
C
                  WRITE (6,'(/,10x,"Type IVFT = ",i5)') IVFT
C
                  IF ( .NOT.INTEG ) THEN
                     DO MUE = 1,3
                        WRITE (6,99009) (PRESI2*SIGTMP_VFT(MUE,NUE,IVFT)
     &                                  ,NUE=1,3)
                     END DO
                  END IF
               END DO
C                  ------------------------ write conductivities to file
               IF ( ISPINPROJ.EQ.IRESPONSE_ADA_ADA .AND. 
     &              .NOT.INTEG .AND. LPRSIG .AND. NQ.EQ.1 ) THEN
                  DO MUE = 1,3
                     DO NUE = 1,3
                        WRITE (MSTR(1:1),'(i1)') MUE
                        WRITE (NSTR(1:1),'(i1)') NUE
                        FILNAM = DATSET(1:LDATSET)//MSTR//NSTR//
     &                           '_SIG_Streda.dat'
                        OPEN (IOTMP,FILE=FILNAM,STATUS='replace')
                        WRITE (IOTMP,99011)
C
                        WRITE (IOTMP,99010) X_VFT(1),
     &                         PRESI2*SIGTMP(MUE,NUE),
     &                         (PRESI2*SIGTMP_VFT(MUE,NUE,I),I=1,2),
     &                         (PRESI2*SIGTMPSS_VFT(MUE,NUE,I),I=1,2)
                        CLOSE (IOTMP)
                     END DO
                           ! NUE
                  END DO ! MUE
               END IF
C
            END IF
C
            IF ( IPRINT.GE.2 .AND. 1.EQ.0 ) THEN
               WRITE (6,99005)
               DO MUE = 1,3
                  WRITE (6,99008) (DREAL(SIGOFFQ(MUE,NUE,ISPINPROJ,IQ)),
     &                            NUE=1,3)
               END DO
C
            END IF
C
C --------------------------- check if imaginary part of SIGOFFQ is zero
C
            DO MUE = 1,3
               DO NUE = 1,3
                  IF ( ABS(DIMAG(SIGOFFQ(MUE,NUE,ISPINPROJ,IQ)))
     &                 .GT.1D-5 ) WRITE (6,99006)
     &                 DIMAG(SIGOFFQ(MUE,NUE,ISPINPROJ,IQ)),MUE,NUE,IQ
               END DO
            END DO
C
         END DO
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      END DO
C
      DEALLOCATE (IMTAU_VFT,ME,MET,IMTSS_VFT,SIGTMPSS_VFT,SIGTMP_VFT,
     &            STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC IMTAU_VFT')
C
99001 FORMAT (//,' J - matrix elements for component IT=',I3)
99002 FORMAT (//,1X,79('*'),/,37X,'<SIGOFF>',/,1X,79('*'),/)
99003 FORMAT (/,10X,12('='),/,12X,'IQ = ',I3,/,10X,12('='),/)
99004 FORMAT (/,10X,'Streda Term   (SI - units)',/)
99005 FORMAT (/,10X,'site-diagonal term sigma_0',/)
99006 FORMAT (' WARNING!! Im(sigma) =',e13.5,' for mue,nue=',2I2,
     &        '  IQ=',I2)
99007 FORMAT (3(F14.6,F12.6))
99008 FORMAT (10X,3F14.6)
99009 FORMAT (3(F16.10,F16.10))
99010 FORMAT (99E17.8)
99011 FORMAT ('#',5x,'conc',10x,'sig_streda_sum',19X,'sigs1',29x,
     &        'sigs2',29x,'sigs_ss1',29x,'sigs_ss2')
      END
C*==sigoffme.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGOFFME(IT,IFILWF,ME,ISPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *  read wave function and calculate the off-diagonal               *
C   *  of the conductivity tensor                                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:IMT,NCPLWFMAX,IKMCPLWF
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R2DRDI,R
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,IMKM_IKM,A_RSIG,NCPLWF
      USE MOD_CONSTANTS,ONLY:CI,C0
      USE MOD_SIG,ONLY:IRESPONSE_EDELSTEIN,IRESPONSE_ADA_ADA
      IMPLICIT NONE
C*--SIGOFFME527
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGOFFME')
C
C Dummy arguments
C
      INTEGER IFILWF,ISPINPROJ,IT
      COMPLEX*16 ME(NKMMAX,NKMMAX,-1:+1,-1:+1)
C
C Local variables
C
      REAL*8 AME1,AME2,RWGT
      COMPLEX*16 CINT1(:),CINT2(:),JFX(1,1,1),JGX(1,1,1),R1A1,R2A2,
     &           RME1(2,2),RME2(2,2),ZFF(:,:,:),ZFI(:,:,:),ZGF(:,:,:),
     &           ZGI(:,:,:)
      INTEGER IA_ERR,IKMTF,IKMTI,IKMTOP,IM,IR,IRTOP,JF,JI,JMF,JMI,KF,KI,
     &        LAM1,LAM2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CINT1,CINT2,ZGF,ZFF,ZFI,ZGI
C
      ALLOCATE (CINT1(NRMAX),CINT2(NRMAX))
      ALLOCATE (ZFI(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGI(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGF(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
      IKMTOP = NKM
C
      ME(1:NKMMAX,1:NKMMAX,-1:+1,-1:+1) = C0
C
      CALL WAVFUN_READ_REL(IFILWF,IT,0,ZGF,ZFF,JGX,JFX,IRTOP,NCPLWF,
     &                     IKMCPLWF)
C
      CALL WAVFUN_READ_REL(IFILWF,IT,0,ZGI,ZFI,JGX,JFX,IRTOP,NCPLWF,
     &                     IKMCPLWF)
C
C=======================================================================
      DO IKMTF = 1,IKMTOP
C
C-----------------------------------------------------------------------
C
         DO IKMTI = 1,IKMTOP
C
C------------------------------------------------- check selection rules
C
            DO KI = 1,NCPLWF(IKMTI)
               JI = IKMCPLWF(KI,IKMTI)
               JMI = IMKM_IKM(JI)
               DO KF = 1,NCPLWF(IKMTF)
                  JF = IKMCPLWF(KF,IKMTF)
                  JMF = IMKM_IKM(JF)
C
                  DO LAM1 = -1, + 1
                     DO LAM2 = -1, + 1
                        AME1 = A_RSIG(JF,JMI,LAM1,LAM2)
                        AME2 = A_RSIG(JMF,JI,LAM1,LAM2)
C
                        IF ( ABS(AME1).GT.1D-6 ) GOTO 20
                        IF ( ABS(AME2).GT.1D-6 ) GOTO 20
C
                     END DO
                  END DO
               END DO
            END DO
C
            CYCLE
C
C ------------------------------------- calculate radial matrix elements
C
 20         CONTINUE
            DO KI = 1,NCPLWF(IKMTI)
               DO KF = 1,NCPLWF(IKMTF)
C
C----------------------------------------------------------- volume term
C
                  DO IR = 1,IRTOP
                     RWGT = R(IR,IM)*R2DRDI(IR,IM)
                     CINT1(IR) = ZGF(IR,KF,IKMTF)*ZFI(IR,KI,IKMTI)*RWGT
                     CINT2(IR) = ZFF(IR,KF,IKMTF)*ZGI(IR,KI,IKMTI)*RWGT
                  END DO
C
                  CALL CRADINT(IM,CINT1,RME1(KF,KI))
                  CALL CRADINT(IM,CINT2,RME2(KF,KI))
C
C-------------------------------------------------- alpha - surface term
C
C                       ZFFZFI(KF,KI) = ZFF(IRTOP,KF)*ZFI(IRTOP,KI)
C                       ZGFZGI(KF,KI) = ZGF(IRTOP,KF)*ZGI(IRTOP,KI)
C
               END DO
            END DO
C
C -------------------------------------- calculate total matrix elements
C
            DO KI = 1,NCPLWF(IKMTI)
               JI = IKMCPLWF(KI,IKMTI)
               JMI = IMKM_IKM(JI)
               DO KF = 1,NCPLWF(IKMTF)
                  JF = IKMCPLWF(KF,IKMTF)
                  JMF = IMKM_IKM(JF)
C
                  DO LAM1 = -1, + 1
                     DO LAM2 = -1, + 1
                        AME1 = A_RSIG(JF,JMI,LAM1,LAM2)
                        AME2 = A_RSIG(JMF,JI,LAM1,LAM2)
                        R1A1 = RME1(KF,KI)*AME1
                        R2A2 = RME2(KF,KI)*AME2
C
                        IF ( ISPINPROJ.EQ.IRESPONSE_ADA_ADA ) THEN
C
                           ME(IKMTF,IKMTI,LAM1,LAM2)
     &                        = ME(IKMTF,IKMTI,LAM1,LAM2)
     &                        + CI*(R1A1-R2A2)
C
C----------------------------------------------------------------------up
CSW!!!!!!assuming old list of operators!!!!!!
                        ELSE IF ( ISPINPROJ.EQ.6 ) THEN
CSW!!!!!!assuming old list of operators!!!!!!
C
                           IF ( LAM2.EQ.+1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2) + CI*R1A1
                           ELSE IF ( LAM2.EQ.-1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2) - CI*R2A2
                           END IF
C-----------------------------------------------------------------------
C
C
C--------------------------------------------------------------------down
CSW!!!!!!assuming old list of operators!!!!!!
                        ELSE IF ( ISPINPROJ.EQ.7 ) THEN
CSW!!!!!!assuming old list of operators!!!!!!
C
                           IF ( LAM2.EQ.+1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2) - CI*R2A2
                           ELSE IF ( LAM2.EQ.-1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2) + CI*R1A1
                           END IF
C
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------beta Sigma_z
CSW!!!!!!assuming old list of operators!!!!!!
                        ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN )
     &                            THEN
CSW!!!!!!assuming old list of operators!!!!!!
C
                           IF ( LAM2.EQ.+1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           - CI*(R1A1+R2A2)
                           ELSE IF ( LAM2.EQ.-1 ) THEN
                              ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           = ME(IKMTF,IKMTI,LAM1,LAM2)
     &                           + CI*(R1A1+R2A2)
                           END IF
C
C-----------------------------------------------------------------------
                        END IF
C
                     END DO
                  END DO
               END DO
            END DO
C
         END DO
C
      END DO
C
      DEALLOCATE (CINT1,CINT2,ZGF,ZFF,ZFI,ZGI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC SIGOFFME')
C
      END
C*==ame_rsig_spin.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE AME_RSIG_SPIN(ISPINPROJ)
C   ********************************************************************
C   *                                                                  *
C   *   calculate angular matrix elements occuring for the             *
C   *           full conductivity tensor                               *
C   *                                                                  *
C   *  A_RSIG = < KAP1,MUE1| r_lam * SIG_lam' |KAP2,MUE2>              *
C   *                                                                  *
C   *  14/11/08  HE                                                    *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:PI,SQRT_2
      USE MOD_ANGMOM,ONLY:A_RSIG,A_SIGMA,CGC,NLMAX,NKMPMAX
      USE MOD_SIG,ONLY:IRESPONSE_ADA_ADA
      IMPLICIT NONE
C*--AME_RSIG_SPIN742
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ISPINPROJ
C
C Local variables
C
      REAL*8 AMP,APM,CGNT,F,NORM_Y1M,RSUM,SIG
      REAL*8 GAUNT_CYLM
      INTEGER IKM1,IKM2,IMS1,IMS2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,L2,
     &        LAM1,LAM2,LZ,M1,M2,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      A_RSIG(1:NKMPMAX,1:NKMPMAX,-1:+1,-1:+1) = 0D0
C
      NK = 2*NLMAX
      NORM_Y1M = DSQRT(4D0*PI/3D0)
C
C=======================================================================
      A_SIGMA(1:2,-1:+1,1:2) = 0D0
C
      A_SIGMA(1,-1,2) = SQRT_2
      A_SIGMA(1,0,1) = -1D0
C
      A_SIGMA(2,0,2) = +1D0
      A_SIGMA(2,+1,1) = -SQRT_2
C=======================================================================
C
C=======================================================================
      DO LAM1 = -1, + 1
         DO LAM2 = -1, + 1
C
C ----------------------------------------------------------------------
            IKM1 = 0
            DO K1 = 1,NK
               L1 = K1/2
               IF ( MOD(K1,2).EQ.0 ) THEN
                  KAP1 = L1
               ELSE
                  KAP1 = -L1 - 1
               END IF
               J1P05 = ABS(KAP1)
C
               DO MUE1M05 = -J1P05,J1P05 - 1
                  IKM1 = IKM1 + 1
C
                  IKM2 = 0
                  DO K2 = 1,NK
                     L2 = K2/2
                     IF ( MOD(K2,2).EQ.0 ) THEN
                        KAP2 = L2
                     ELSE
                        KAP2 = -L2 - 1
                     END IF
                     J2P05 = ABS(KAP2)
C
                     DO MUE2M05 = -J2P05,J2P05 - 1
                        IKM2 = IKM2 + 1
C
C ----------------------------------------------------------------------
CSW!!!!!!assuming old list of operators!!!!!!
                        IF ( ISPINPROJ.EQ.11 ) THEN
CSW!!!!!!assuming old list of operators!!!!!!
C
                           RSUM = 0D0
C
                           DO IMS1 = 1,2
                              M1 = MUE1M05 - IMS1 + 2
C
                              DO IMS2 = 1,2
                                 M2 = MUE2M05 - IMS2 + 2
C
C
                                 IF ( IMS2.EQ.1 ) LZ = MUE2M05 + 1
                                 IF ( IMS2.EQ.2 ) LZ = MUE2M05
C
C
                                 CGNT = GAUNT_CYLM(L1,M1,1,LAM1,L2,M2)
                                 SIG = A_SIGMA(IMS1,LAM2,IMS2)*LZ
C
                                 RSUM = RSUM + CGC(IKM1,IMS1)
     &                                  *CGC(IKM2,IMS2)*CGNT*SIG
C
C
                              END DO
                           END DO
C
                           A_RSIG(IKM1,IKM2,LAM1,LAM2) = NORM_Y1M*RSUM
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
                        ELSE
                           RSUM = 0D0
C
                           DO IMS1 = 1,2
                              M1 = MUE1M05 - IMS1 + 2
C
                              DO IMS2 = 1,2
                                 M2 = MUE2M05 - IMS2 + 2
C
                                 CGNT = GAUNT_CYLM(L1,M1,1,LAM1,L2,M2)
                                 SIG = A_SIGMA(IMS1,LAM2,IMS2)
C
                                 RSUM = RSUM + CGC(IKM1,IMS1)
     &                                  *CGC(IKM2,IMS2)*CGNT*SIG
C
                              END DO
                           END DO
C
                           A_RSIG(IKM1,IKM2,LAM1,LAM2) = NORM_Y1M*RSUM
C
                        END IF
C ----------------------------------------------------------------------
C
C
C
                     END DO
                  END DO
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END DO
      END DO
C
C
C
C ======================================================================
C
      F = SQRT_2*NORM_Y1M
C
C ======================================================================
C                       check for (-,+) and (+,-)
C=======================================================================
      IKM1 = 0
      IF ( ISPINPROJ.EQ.IRESPONSE_ADA_ADA ) THEN
         DO K1 = 1,NK
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            J1P05 = ABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  J2P05 = ABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     IKM2 = IKM2 + 1
C
C ----------------------------------------------------------------------
                     IMS1 = 2
                     IMS2 = 1
                     M1 = MUE1M05 - IMS1 + 2
                     M2 = MUE2M05 - IMS2 + 2
                     CGNT = GAUNT_CYLM(L1,M1,1,-1,L2,M2)
C
                     AMP = -F*CGC(IKM1,IMS1)*CGC(IKM2,IMS2)*CGNT
C
                     IMS1 = 1
                     IMS2 = 2
                     M1 = MUE1M05 - IMS1 + 2
                     M2 = MUE2M05 - IMS2 + 2
                     CGNT = GAUNT_CYLM(L1,M1,1,+1,L2,M2)
C
                     APM = F*CGC(IKM1,IMS1)*CGC(IKM2,IMS2)*CGNT
C
                     IF ( ABS(A_RSIG(IKM1,IKM2,-1,+1)-AMP).GT.1D-8 .OR. 
     &                    ABS(A_RSIG(IKM1,IKM2,+1,-1)-APM).GT.1D-8 )
     &                    THEN
                        WRITE (6,*) ' TROUBLE in <AME_RSIG_SPIN>'
                        WRITE (6,*) ' IKM1 = ',IKM1,'    IKM2 = ',IKM2
                        WRITE (6,*) ' A(-,+) ',A_RSIG(IKM1,IKM2,-1,+1),
     &                              AMP
                        WRITE (6,*) ' A(+,-) ',A_RSIG(IKM1,IKM2,+1,-1),
     &                              APM
                     END IF
C ----------------------------------------------------------------------
C
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C ======================================================================
C
      END
