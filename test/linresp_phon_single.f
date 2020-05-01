C*==linresp_phon_single.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LINRESP_PHON_SINGLE
C   ********************************************************************
C   *                                                                  *
C   *  calculate the dynamical matrix  PHI_ij(q)                       *
C   *                                                                  *
C   *  using the single particle energy approximation                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQMAX,IQAT,ITOQ,NOQ,NQ,IQBOT_CHI,IQTOP_CHI
      USE MOD_ANGMOM,ONLY:NL,NKMQ,NLM,NKMMAX,NKM,NXM,WKM1,WKM2,MEZJ,
     &    MEZZ,TSSQ,MSSQ,MSST,SSST,TAUQ,TAUT,TSST
      USE MOD_CALCMODE,ONLY:ORBPOL,KMROT,IREL
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,NETAB,NEPATH,WETAB,IGRID,EFERMI,
     &    EIMAG,EMIN,EILOW,PHASK
      USE MOD_TYPES,ONLY:NTMAX,CONC,NT
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET,IPRINT,IFILCBWF,
     &    FOUND_SECTION
      USE MOD_CONSTANTS,ONLY:CI,C0,PI,RY_EV,EV_ERG,A0_CGS
      USE MOD_TAUIJ,ONLY:TAUIJ
      USE MOD_LINRESP,ONLY:NZ12MAX,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2,NTTJ,
     &    JTT1,JTT2,JTTX,WTTJ,ITTMAX,JTTMAX,NTKTKMAX,CHIZ
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_MPI,ONLY:MPI,MPI_KLOOP,MPI_ELOOP
      IMPLICIT NONE
C*--LINRESP_PHON_SINGLE27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LINRESP_PHON_SINGLE')
C
C Local variables
C
      REAL*8 AUTOCGS,CPACHNG,CPACHTAB(NEMAX),DQVEC_PERT,PREFAC,
     &       QHAT_PERT(3),QPERT,QVEC_PERT(3),QVEC_PERT0,TIME1,TIME2
      LOGICAL CALCINT,CHECK,GETIRRSOL
      COMPLEX*16 CFAC,CSUM,CSUM1,CSUM2,CSUM3,DELMSSQ(:,:,:,:),
     &           DELMSST(:,:,:,:),DMAMC(:,:),DMATT(:,:,:),DTILT(:,:,:),
     &           ERYD,P,PHI_IJ(:,:,:,:),
     &           TAUQZ(NKMMAX,NKMMAX,NQMAX,NZ12MAX),THT_OFF_QC(:,:,:,:),
     &           UBAR(:,:,:),WBAR(:,:,:),WE
      REAL*8 DNRM2,GAUNT_RYLM
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IC,ICPACONV,ICPAFLAG,IE,IECPAFAIL(NEMAX),IO,IQ,
     &        IQVEC_PERT,IQVEC_PERT1,IQVEC_PERT2,IS,IT,ITCPA,IWRI,
     &        IWRIRRWF,IWRREGWF,J,JC,JQ,K1,K2,L1,L2,L3,L4,LFN,LM1,LM2,
     &        LM3,LM4,M,M1,M2,M3,M_CART(3),N,NCPAFAIL,NKMSQ
C
C*** End of declarations rewritten by SPAG
C
      DATA M_CART/ + 1, - 1,0/,CHECK/.TRUE./
C
      ALLOCATABLE DMATT,DTILT
      ALLOCATABLE WBAR,UBAR
      ALLOCATABLE PHI_IJ
      ALLOCATABLE DMAMC,DELMSST,DELMSSQ,THT_OFF_QC
C
      ALLOCATE (THT_OFF_QC(NKMMAX,NKMMAX,NQ,3))
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DTILT(NKMMAX,NKMMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate DTILT')
C
      IF ( NEPATH.NE.1 ) CALL STOP_MESSAGE(ROUTINE,'NEPATH <> 0')
C
      DEALLOCATE (ETAB,WETAB)
      IGRID(1) = 5
      NEMAX = 20
      NETAB(1) = NEMAX
      ALLOCATE (ETAB(NEMAX,2),WETAB(NEMAX,2))
C
      CALL EPATH(IGRID(1),EMIN,EFERMI,EIMAG,NETAB(1),ETAB(1,1),
     &           WETAB(1,1),EILOW,3,NEMAX)
C
      IF ( IGRID(1).NE.5 ) CALL STOP_MESSAGE(ROUTINE,'IGRID <> 5')
      IF ( IREL.EQ.2 .AND. KMROT.NE.0 )
     &      CALL STOP_MESSAGE(ROUTINE,'KMROT <> 0 for IREL = 2')
C
C-----------------------------------------------------------------------
C                           initialize
C-----------------------------------------------------------------------
C
      ALLOCATE (DMAMC(NKMMAX,NKMMAX),UBAR(NKMMAX,NKMMAX,3))
      ALLOCATE (DELMSSQ(NKMMAX,NKMMAX,NQ,3))
      ALLOCATE (DELMSST(NKMMAX,NKMMAX,NT,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc: DELMSST')
C
      AUTOCGS = RY_EV*EV_ERG/(A0_CGS*A0_CGS)
C
C***********************************************************************
C                  initialize PHI_IJ and WBAR
C***********************************************************************
C
      ALLOCATE (PHI_IJ(NQ,NQ,3,3),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc: PHI_IJ')
C
C-----------------------------------------------------------------------
C               auxilary transformation matrix  WBAR
C-----------------------------------------------------------------------
      ALLOCATE (WBAR(NKM,NKM,3))
C
      WBAR(:,:,:) = C0
C
C-------------------------------------- initiate WBAR for 1st spin index
C
      PREFAC = SQRT(4D0*PI/3D0)
C
      LM1 = 0
      DO L1 = 0,(NL-1)
         DO M1 = -L1,L1
            LM1 = LM1 + 1
C
            LM2 = 0
            DO L2 = 0,(NL-1)
               CFAC = PREFAC*CI**(L1-L2+1)
               DO M2 = -L2,L2
                  LM2 = LM2 + 1
C
                  L3 = 1
                  DO IC = 1,3
                     M3 = M_CART(IC)
                     WBAR(LM1,LM2,IC) = CFAC*GAUNT_RYLM(L1,M1,L2,M2,L3,
     &                                  M3)
                  END DO
C
               END DO
            END DO
         END DO
      END DO
C
C-------------------------------------- initiate WBAR for 2nd spin index
C
      IF ( IREL.GE.2 ) WBAR(NLM+1:NKM,NLM+1:NKM,1:3)
     &     = WBAR(1:NLM,1:NLM,1:3)
C
C--------------------------------- convert to (kappa,mue) representation
      IF ( IREL.EQ.3 ) THEN
         DO IC = 1,3
C
            CALL CHANGEREP(NKM,NKMMAX,WBAR(1,1,IC),'RLM>REL',WKM1)
C
            WBAR(1:NKM,1:NKM,IC) = WKM1(1:NKM,1:NKM)
C
         END DO
      END IF
C
C.......................................................................
      IF ( CHECK ) THEN
         DO IC = 1,3
            CALL CMATSTRUCT('WBAR  '//(CHAR(ICHAR('1')+IC-1)),
     &                      WBAR(1,1,IC),NKM,NKM,IREL,IREL,0,1D-8,6)
         END DO
      END IF
C.......................................................................
C
C***********************************************************************
C
C
C
C ======================================================================
C         set up index table for BZ integration for NO SYMMETRY case
C ======================================================================
C
      IF ( ALLOCATED(WTTJ) ) THEN
         DEALLOCATE (WTTJ,JTT1,JTT2,JTTX)
         DEALLOCATE (NTTJ,ITTA,ITTB,ITTC,ITTD,ITTQ1,ITTQ2)
      END IF
C
      JTTMAX = 1
C
      ALLOCATE (WTTJ(JTTMAX))
      ALLOCATE (JTT1(JTTMAX),JTT2(JTTMAX),JTTX(JTTMAX))
C
C------------------------------- consider ALL TAU(k)*TAU(k) combinations
      NTKTKMAX = NKM**4
C
      ITTMAX = NTKTKMAX*NQ*NQ
C
      ALLOCATE (NTTJ(ITTMAX))
      ALLOCATE (ITTA(ITTMAX),ITTB(ITTMAX),ITTC(ITTMAX))
      ALLOCATE (ITTD(ITTMAX),ITTQ1(ITTMAX),ITTQ2(ITTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate ITTD')
C
      I = 0
      DO IQ = IQBOT_CHI,IQTOP_CHI
         DO JQ = IQBOT_CHI,IQTOP_CHI
            DO L1 = 1,NKM
               DO L4 = 1,NKM
                  DO L2 = 1,NKM
                     DO L3 = 1,NKM
                        I = I + 1
                        ITTA(I) = L1
                        ITTB(I) = L2
                        ITTC(I) = L3
                        ITTD(I) = L4
                        ITTQ1(I) = IQ
                        ITTQ2(I) = JQ
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C ======================================================================
C                   initialize k-mesh for FULL BZ
C ======================================================================
C
      IPRINT = 4
      CALL INIT_MOD_TAUIJ_KMESH
C
C ======================================================================
C
      IWRI = 6
C
      IPRINT = 1
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      IWRREGWF = 1
      IWRIRRWF = 1
C
C=======================================================================
C
      NCPAFAIL = 0
C
      CALL CPU_TIME(TIME1)
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      IF ( .NOT.(MPI) ) THEN
         MPI_KLOOP = .FALSE.
         MPI_ELOOP = .FALSE.
      ELSE IF ( MPI ) THEN
         MPI_ELOOP = .FALSE.
         MPI_KLOOP = .TRUE.
      END IF
C
      IF ( MPI ) CALL DRV_MPI_BARRIER
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C
C
C ======================================================================
C                        q_pert - vectors
C ======================================================================
C
      IQVEC_PERT1 = 0
      IQVEC_PERT2 = 11
      QVEC_PERT0 = 0.0D0
      QHAT_PERT(1) = 1.0D0
      QHAT_PERT(2) = 0.0D0
      QHAT_PERT(3) = 0.0D0
      QPERT = DNRM2(3,QHAT_PERT,1)
      QHAT_PERT(1:3) = QHAT_PERT(1:3)/QPERT
      DQVEC_PERT = QPERT/DBLE(IQVEC_PERT2-1)
C
      CALL INPUT_FIND_SECTION('MODE',0)
C
      IF ( FOUND_SECTION ) THEN
         CALL SECTION_SET_INTEGER('IPERT1',IQVEC_PERT1,9999,0)
         CALL SECTION_SET_INTEGER('IPERT2',IQVEC_PERT2,9999,0)
         CALL SECTION_SET_REAL('QPERT0',QVEC_PERT0,9999D0,0)
         CALL SECTION_SET_REAL('DQPERT',DQVEC_PERT,9999D0,0)
      END IF
C
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                        loop over q_pert - vectors
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      DO IQVEC_PERT = IQVEC_PERT1,IQVEC_PERT2
C
         QVEC_PERT(1:3) = QVEC_PERT0 + DQVEC_PERT*IQVEC_PERT*QHAT_PERT
     &                    (1:3)
         WRITE (6,99003) IQVEC_PERT,QVEC_PERT
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         IS = 1
C
         DO IE = 1,NETAB(1)
C
            ERYD = ETAB(IE,1)
C
            IF ( IREL.EQ.1 ) THEN
C
               WE = WETAB(IE,1)*2D0
C
            ELSE
C
               WE = WETAB(IE,1)
C
            END IF
C
            ICPAFLAG = 0
            CPACHNG = 0.0D0
C
            CALL CINIT(NKMMAX*NKMMAX*NTMAX,TAUT)
C
C ===================================== solve SS - differential equation
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,
     &                    ERYD,P,IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
C ======================================================================
C
            IF ( NCPA.GT.0 ) CALL TAU_DRIVE(IE,IPRINT,ERYD,P,TSSQ,MSSQ,
     &           TSST,MSST,TAUQ,ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
            IF ( IREL.EQ.2 ) THEN
C
               CALL SIGKLOOPS_SPSREL(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,
     &                               MSSQ)
C
            ELSE
C
               CALL SIGKLOOPS(ERYD,P,ERYD,P,TAUQ,TAUQ,TAUQZ,MSSQ,MSSQ)
C
            END IF
C
C ======================================================================
C
C-----------------------------------------------------------------------
C                  reduced transformation matrix  UBAR
C-----------------------------------------------------------------------
            DO IC = 1,3
               DO J = 1,NKM
                  DO I = 1,NKM
                     UBAR(I,J,IC) = P*WBAR(I,J,IC)
                  END DO
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                      DELMSST  =  U * m - m U
C-----------------------------------------------------------------------
C
            DO IT = 1,NT
C
               DO IC = 1,3
C
                  CALL CMATMUL(NKM,NKMMAX,UBAR(1,1,IC),MSST(1,1,IT),
     &                         WKM1)
                  CALL CMATMUL(NKM,NKMMAX,MSST(1,1,IT),UBAR(1,1,IC),
     &                         WKM2)
C
                  DELMSST(1:NKM,1:NKM,IT,IC) = WKM1(1:NKM,1:NKM)
     &               - WKM2(1:NKM,1:NKM)
C
               END DO
C
C.......................................................................
               IF ( CHECK ) THEN
                  DO IC = 1,3
                     CALL CMATSTRUCT('DELM '//' '//(CHAR(ICHAR('1')+IC-1
     &                               )),DELMSST(1,1,IT,IC),NXM,NXM,IREL,
     &                               IREL,0,1D-8,6)
                  END DO
               END IF
C.......................................................................
C
C--------------------------- calculate projection matrices DMAT and DTIL
C
               IQ = IQAT(1,IT)
               N = NKMQ(IQ)
               M = NKMMAX
C
               CALL GETDMAT(TAUQ(1,1,IQ),DMATT(1,1,IT),DTILT(1,1,IT),
     &                      DMAMC,N,MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
            END DO
C.......................................................................
            IF ( CHECK ) THEN
               WRITE (*,*) '####',IREL
               IT = 1
               CALL CMATSTRUCT('DMAT   '//(CHAR(ICHAR('1')+IS-1)),
     &                         DMATT(1,1,IT),N,M,IREL,IREL,0,1D-8,6)
               CALL CMATSTRUCT('DTIL   '//(CHAR(ICHAR('1')+IS-1)),
     &                         DTILT(1,1,IT),N,M,IREL,IREL,0,1D-8,6)
            END IF
C.......................................................................
C
C=======================================================================
C                calculate site-off force tensor   PHI_ij
C       in case of alloy system: perform projection on atom types
C=======================================================================
C
            N = NKM
            M = NKMMAX
            NKMSQ = NKM*NKM
C
C-----------------------------------------------------------------------
C             type averaged Delta m matrix for site IQ
C-----------------------------------------------------------------------
C
            DELMSSQ(:,:,:,:) = C0
C
            DO IQ = 1,NQ
               DO IO = 1,NOQ(IQ)
                  IT = ITOQ(IO,IQ)
                  DO IC = 1,3
C
                     WKM1(1:N,1:N) = MATMUL(DTILT(1:N,1:N,IT),DELMSST(1:
     &                               N,1:N,IT,IC))
C
                     WKM2(1:N,1:N) = MATMUL(WKM1(1:N,1:N),DMATT(1:N,1:N,
     &                               IT))
C
                     DELMSSQ(1:N,1:N,IQ,IC) = DELMSSQ(1:N,1:N,IT,IC)
     &                  + CONC(IT)*WKM2(1:N,1:N)
C
                  END DO
               END DO
            END DO
C
C-----------------------------------------------------------------------
C
            PHI_IJ(:,:,:,:) = 0D0
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            DO IQ = 1,NQ
C
               THT_OFF_QC(:,:,:,:) = 0D0
C
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
               DO JQ = IQBOT_CHI,IQTOP_CHI
C
                  K1 = (IQ-1)*NKMSQ
                  DO LM1 = 1,NKM
                     DO LM4 = 1,NKM
                        K1 = K1 + 1
C
                        CSUM1 = 0D0
                        CSUM2 = 0D0
                        CSUM3 = 0D0
C
                        K2 = (JQ-1)*NKMSQ
                        DO LM2 = 1,NKM
                           DO LM3 = 1,NKM
                              K2 = K2 + 1
C
                              CSUM1 = CSUM1 + CHIZ(K1,K2,1)
     &                                *DELMSSQ(LM2,LM3,JQ,1)
C
                              CSUM2 = CSUM2 + CHIZ(K1,K2,1)
     &                                *DELMSSQ(LM2,LM3,JQ,2)
C
                              CSUM3 = CSUM3 + CHIZ(K1,K2,1)
     &                                *DELMSSQ(LM2,LM3,JQ,3)
C
                           END DO
                        END DO
C
                        THT_OFF_QC(LM1,LM4,JQ,1)
     &                     = THT_OFF_QC(LM1,LM4,JQ,1) + CSUM1
C
                        THT_OFF_QC(LM1,LM4,JQ,2)
     &                     = THT_OFF_QC(LM1,LM4,JQ,2) + CSUM2
C
                        THT_OFF_QC(LM1,LM4,JQ,3)
     &                     = THT_OFF_QC(LM1,LM4,JQ,3) + CSUM3
C
                     END DO
                  END DO
C
                  DO IC = 1,3
                     DO JC = 1,3
C
                        WKM1(1:N,1:N) = MATMUL(DELMSSQ(1:N,1:N,IQ,IC),
     &                                  THT_OFF_QC(1:N,1:N,JQ,JC))
C
                        CSUM = C0
                        DO I = 1,NKM
                           CSUM = CSUM + WKM1(I,I)
                        END DO
C
                        PHI_IJ(IQ,JQ,IC,JC) = PHI_IJ(IQ,JQ,IC,JC)
     &                     + WE*CSUM
C
                     END DO
                  END DO
C
               END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
            END DO
C   QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C ======================================================================
C
            IF ( ICPAFLAG.NE.0 ) THEN
               NCPAFAIL = NCPAFAIL + 1
               CPACHTAB(NCPAFAIL) = CPACHNG
               IECPAFAIL(NCPAFAIL) = IE
            END IF
C
C ======================================================================
            IF ( IPRINT.GE.3 ) CALL DUMPTAU(IE,ERYD,IWRI,MSST,MSSQ,TAUT,
     &           TAUQ)
C
         END DO
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
         CALL CPU_TIME(TIME2)
C
         WRITE (IWRI,99004) TIME2 - TIME1
C
         IF ( NCPAFAIL.NE.0 ) THEN
            WRITE (IWRI,99005) CPATOL,NCPAFAIL,
     &                         (IECPAFAIL(IE),DREAL(ETAB(IECPAFAIL(IE),
     &                         1)),CPACHTAB(IE),IE=1,NCPAFAIL)
            WRITE (IWRI,'(1X,79(''*''),/)')
         ELSE IF ( NCPA.NE.0 ) THEN
            WRITE (IWRI,99006)
         END IF
C
C-----------------------------------------------------------------------
C                              write results
C-----------------------------------------------------------------------
C
         IF ( IQVEC_PERT.EQ.1 ) THEN
            FILNAM = DATSET(1:LDATSET)//'_F_ij.dat'
            LFN = LDATSET + 9
            CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFN))
C
            WRITE (*,99008) ROUTINE
            WRITE (IOTMP,99008) ROUTINE
C
            WRITE (IOTMP,99010)
            DO IQ = 1,NQ
               WRITE (IOTMP,99011) IQ,NOQ(IQ),
     &                             (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,
     &                             NOQ(IQ))
            END DO
         END IF
C
         PHI_IJ(:,:,:,:) = PHI_IJ(:,:,:,:)*AUTOCGS*1D-3/PI
C
         WRITE (IOTMP,99007) QVEC_PERT
         DO IQ = 1,NQ
            DO JQ = 1,NQ
               WRITE (IOTMP,99001)
               DO IC = 1,3
                  WRITE (IOTMP,99002) (PHI_IJ(IQ,JQ,IQ,JC),JC=1,3)
               END DO
            END DO
         END DO
C
      END DO
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C                        loop over q_pert - vectors
Cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
C
      WRITE (*,99009) FILNAM(1:LFN)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE (PHI_IJ,TAUIJ,WKM1,WKM2)
      DEALLOCATE (DMAMC,DTILT,DMATT)
C
C-----------------------------------------------------------------------
      CALL STOP_REGULAR(ROUTINE,' ')
C-----------------------------------------------------------------------
C
99001 FORMAT (2I4,'  IQ  JQ')
99002 FORMAT (3(2E15.7,2X))
99003 FORMAT (I5,'  QVEC_PERT   ',5X,3F10.5)
99004 FORMAT (/,' CPU - TIME  USED      ',F10.3,/)
99005 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99006 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99007 FORMAT (3F10.5,'   q_pert')
99008 FORMAT (/,1X,79('*'),/,35X,'<',A,'>',/,28X,
     &        'Force tensor       PHI_ij',/,1X,79('*'),/)
99009 FORMAT (/,10X,'results written to file:',A,/)
99010 FORMAT (/,10X,'number of sites   NQ = ',I3,/,10X,
     &        'number of types   NT = ',I3,/,10X,'site occupation:')
99011 FORMAT (10X,2I4,10(I3,F6.3))
      END
