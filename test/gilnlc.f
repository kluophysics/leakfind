C*==gilnlc.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE GILNLC(IPRINT,WRMAT,NDIMCLUSQ,IND0QCLU,RQNLCPA,KNNLCPA,
     &                  NSYMACCEPTEDKN,KTABKN,WKTABKN,NKTABKN,
     &                  SYMACCEPTEDKN,NKNIRMU,SYM_IRR,PCFG,MAXNKNIRMU,
     &                  EIKNRIJ,NCFG,NDIMCLU,NKN,NKTABKNMAX,NQNLCPA,
     &                  LCPAMILLS,LNLCPASYMAVG,LNLCPAAVG,IPRINTL,IQCPA)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to calculate the Gilbert damping parameter           *
C   *                                                                  *
C   *  on the basis of the   NON-LOCAL CPA                             *
C   *                                                                  *
C   ********************************************************************
      USE MOD_LINRESP,ONLY:CHIZ,NZ12,NZ12MAX
      USE MOD_CPA,ONLY:CPATOL,ALPHASRO,CPAMIX,CPALVL,USENLCPA,ITCPAMAX,
     &    NCPA
      USE MOD_KSPACE,ONLY:IBZINT
      USE MOD_SYMMETRY,ONLY:SYMUNITARY,IQORGQP,DROT,MROTK,NSYM,NSYMMAX
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,NKM,MEZJ,MEZZ,TSSQ,MSSQ,MSST,SSST,
     &    TAUQ,TAUT,TSST
      USE MOD_SITES,ONLY:NQ,NOMAX,NQMAX,ITOQ,NOQ,IQAT
      USE MOD_TYPES,ONLY:NT,NTMAX,CONC
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_ENERGY,ONLY:ETAB,NEPATH,NETAB,EFERMI,PHASK
      USE MOD_FILES,ONLY:IFILCBWF
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--GILNLC28
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPRINT,IPRINTL,IQCPA,MAXNKNIRMU,NCFG,NDIMCLU,NDIMCLUSQ,
     &        NKN,NKTABKNMAX,NQNLCPA
      LOGICAL LCPAMILLS,LNLCPAAVG,LNLCPASYMAVG,WRMAT
      COMPLEX*16 EIKNRIJ(NQNLCPA,MAXNKNIRMU,NQNLCPA,NQNLCPA)
      INTEGER IND0QCLU(NQNLCPA),NKNIRMU(NKN),NKTABKN(NKN),
     &        NSYMACCEPTEDKN(NKN)
      REAL*8 KNNLCPA(3,NQNLCPA),KTABKN(3,NKTABKNMAX,NKN),PCFG(NCFG),
     &       RQNLCPA(3,NQNLCPA),WKTABKN(NKTABKNMAX,NKN)
      LOGICAL SYMACCEPTEDKN(NSYMMAX,NKN),SYM_IRR(NKN,MAXNKNIRMU,NSYM)
C
C Local variables
C
      COMPLEX*16 ALF0Q(3,3,NQMAX),ALF0QO(:,:,:,:),ALF1QOQ_NV(:,:,:,:,:),
     &           ALF1QOQ_VC(:,:,:,:,:),ALF1QQ_NV(3,3,NQMAX,NQMAX),
     &           ALF1QQ_VC(3,3,NQMAX,NQMAX),CHIHATZ(:,:,:),DERYD,
     &           DMAMC(:,:),DMATTA(:,:,:),DMATTB(:,:,:),DTILTA(:,:,:),
     &           DTILTB(:,:,:),ERYD,ERYD9,ERYDA,ERYDB,MAQAB(:,:,:,:),
     &           MAQBA(:,:,:,:),MAQQAB(:,:,:,:,:),MAQQBA(:,:,:,:,:),
     &           MBQAB(:,:,:,:),MBQBA(:,:,:,:),MBQQAB(:,:,:,:,:),
     &           MBQQBA(:,:,:,:,:),MCQAB(:,:,:,:),MCQBA(:,:,:,:),
     &           MCQQAB(:,:,:,:,:),MCQQBA(:,:,:,:,:),MDQAB(:,:,:,:),
     &           MDQBA(:,:,:,:),MDQQAB(:,:,:,:,:),MDQQBA(:,:,:,:,:),
     &           MIRRTAB(:,:,:,:,:,:),MMAT(:,:),
     &           MSSTA(NKMMAX,NKMMAX,NTMAX),MSSTB(NKMMAX,NKMMAX,NTMAX),
     &           MTAB(:,:,:,:),MTBA(:,:,:,:),MUEHAT(:,:),MUEHATA(:,:),
     &           MUEHATB(:,:),OMEGAHAT(:,:),OMEGAHATA(:,:),
     &           OMEGAHATB(:,:),P,PA,PB,TAUHAT(:,:),TAUIMP(:,:),
     &           TAUQA(:,:,:),TAUQB(:,:,:),TAUQZ(NKMMAX,NKMMAX,NQMAX,2),
     &           TAUTA(:,:,:),TAUTB(:,:,:),W1HAT(:,:)
      LOGICAL CALCINT,GETIRRSOL,INTEG,VERTEX,WRTAUMQ
      REAL*8 CPACHNG,CPACHNGMAX,CPAERR,TIME,TIME0,WGIT(:)
      INTEGER I,I0_QU(:,:),I1_Q(:),I1_QU(:,:),I2_Q(:),I2_QU(:,:),IA_ERR,
     &        ICFG,ICPACONV,ICPAFLAG,IE,IEA,IEB,IEPATH,IFILA,IFILB,
     &        IFILTAU,IOCC,IQ,IT,ITCPA,IU,IWRI,IWRIRRWF,IWRREGWF,J,K,M,
     &        N,NAQ,NAQU,NA_Q(:),NCPAFAIL,NDIMCHI,NEA,NEB,NEINTEG,NETAU,
     &        NN,NU,OUT
      INTEGER NLCPACONF
      CHARACTER*10 SIG_MODE
C
C*** End of declarations rewritten by SPAG
C
      DATA OUT/6/,CALCINT/.TRUE./,GETIRRSOL/.TRUE./
      DATA IWRREGWF/1/,IWRIRRWF/1/
      DATA WRTAUMQ/.FALSE./,VERTEX/.TRUE./,SIG_MODE/'STREDA    '/
      DATA ERYD9/(999999D0,999999D0)/
C
      ALLOCATABLE MTAB,MTBA,DTILTA,DMATTA,DTILTB,DMATTB,DMAMC
      ALLOCATABLE MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA
      ALLOCATABLE MAQQAB,MAQQBA,MBQQAB,MBQQBA
      ALLOCATABLE MCQQAB,MCQQBA,MDQQAB,MDQQBA,MIRRTAB
      ALLOCATABLE TAUQA,TAUQB,TAUTA,TAUTB,CHIHATZ
      ALLOCATABLE OMEGAHATA,OMEGAHATB,OMEGAHAT,MUEHATA,MUEHATB,MUEHAT
      ALLOCATABLE TAUHAT,W1HAT,TAUIMP,MMAT,WGIT
      ALLOCATABLE NA_Q,I0_QU,I1_QU,I2_QU,I1_Q,I2_Q
      ALLOCATABLE ALF0QO,ALF1QOQ_NV,ALF1QOQ_VC
C
      M = NKMMAX
      ALLOCATE (NA_Q(NQ),I0_QU(NQ,NQNLCPA),I1_QU(NQ,NQNLCPA))
      ALLOCATE (I1_Q(NQ),I2_Q(NQ),I2_QU(NQ,NQNLCPA))
      ALLOCATE (DMAMC(M,M))
      ALLOCATE (TAUQA(M,M,NQMAX),TAUQB(M,M,NQMAX))
      ALLOCATE (TAUTA(M,M,NTMAX),TAUTB(M,M,NTMAX))
      ALLOCATE (MAQAB(M,M,3,NQMAX),MAQBA(M,M,3,NQMAX))
      ALLOCATE (MBQAB(M,M,3,NQMAX),MBQBA(M,M,3,NQMAX))
      ALLOCATE (MCQAB(M,M,3,NQMAX),MCQBA(M,M,3,NQMAX))
      ALLOCATE (MDQAB(M,M,3,NQMAX),MDQBA(M,M,3,NQMAX))
      ALLOCATE (MTAB(M,M,3,NTMAX),MTBA(M,M,3,NTMAX))
      ALLOCATE (DTILTA(M,M,NTMAX),DMATTA(M,M,NTMAX))
      ALLOCATE (DTILTB(M,M,NTMAX),DMATTB(M,M,NTMAX))
      ALLOCATE (MIRRTAB(NKMMAX,1,NKMMAX,3,3,NTMAX))
      MIRRTAB(1:NKMMAX,1,1:NKMMAX,1:3,1:3,1:NTMAX) = C0
C
      ALLOCATE (ALF0QO(3,3,NQMAX,NOMAX))
      ALLOCATE (ALF1QOQ_NV(3,3,NQMAX,NOMAX,NQMAX))
      ALLOCATE (ALF1QOQ_VC(3,3,NQMAX,NOMAX,NQMAX))
      ALF0QO(:,:,:,:) = 0D0
      ALF1QOQ_NV(:,:,:,:,:) = 0D0
      ALF1QOQ_VC(:,:,:,:,:) = 0D0
C
      M = NDIMCLU
      ALLOCATE (WGIT(NTMAX))
      ALLOCATE (MUEHATA(M,M),OMEGAHATA(M,M))
      ALLOCATE (MUEHATB(M,M),OMEGAHATB(M,M))
      ALLOCATE (MUEHAT(M,M),OMEGAHAT(M,M),TAUHAT(M,M))
      ALLOCATE (W1HAT(M,M),MMAT(M,M),TAUIMP(M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'in <SIGNLC> -> W1HAT'
C
      NN = NKMMAX*NKMMAX*3*NQMAX
      CALL CINIT(NN,MAQAB)
      CALL CINIT(NN,MBQAB)
      CALL CINIT(NN,MCQAB)
      CALL CINIT(NN,MDQAB)
      CALL CINIT(NN,MAQBA)
      CALL CINIT(NN,MBQBA)
      CALL CINIT(NN,MCQBA)
      CALL CINIT(NN,MDQBA)
C
      NN = NQNLCPA
      ALLOCATE (MAQQAB(M,M,3,NN,NN),MAQQBA(M,M,3,NN,NN))
      ALLOCATE (MBQQAB(M,M,3,NN,NN),MBQQBA(M,M,3,NN,NN))
      ALLOCATE (MCQQAB(M,M,3,NN,NN),MCQQBA(M,M,3,NN,NN))
      ALLOCATE (MDQQAB(M,M,3,NN,NN),MDQQBA(M,M,3,NN,NN),STAT=IA_ERR)
C
      IF ( M.LT.NKM ) WRITE (6,*) 'WARNING:   NKMMAX = ',M,NN
      IF ( IA_ERR.NE.0 ) STOP 'alloc:sigma -> MAQAB'
C
C ======================================================================
C
      NU = NQNLCPA
C
      NAQ = 0
      DO IQ = 1,NQ
         NA_Q(IQ) = NKMQ(IQ)
         NAQ = NAQ + NA_Q(IQ)
         IF ( IQ.EQ.1 ) THEN
            I1_Q(IQ) = 1
         ELSE
            I1_Q(IQ) = I2_Q(IQ-1) + 1
         END IF
         I2_Q(IQ) = I1_Q(IQ) + NA_Q(IQ)
      END DO
C
      NAQU = NAQ*NU
C
      DO IU = 1,NU
         DO IQ = 1,NQ
            IF ( IQ.EQ.1 ) THEN
               I0_QU(IQ,IU) = (IU-1)*NAQ
            ELSE
               I0_QU(IQ,IU) = I0_QU(IQ-1,IU) + NA_Q(IQ-1)
            END IF
            I1_QU(IQ,IU) = I0_QU(IQ,IU) + 1
            I2_QU(IQ,IU) = I0_QU(IQ,IU) + NA_Q(IQ)
         END DO
      END DO
C ======================================================================
C
      CALL CPU_TIME(TIME0)
C
      IF ( IBZINT.NE.2 .AND. IBZINT.NE.6 )
     &      STOP ' in <SIGMA>:  IBZINT <> 2 AND IBZINT <> 6'
C
      WRITE (6,99003)
      IF ( USENLCPA ) THEN
         WRITE (6,99002) IQCPA,CPALVL,CPAMIX,ALPHASRO
      ELSE
         WRITE (6,99001)
      END IF
C
      NCPAFAIL = 0
      ITCPAMAX = 30
C
      IFILTAU = 9
      CALL READTAU(IFILTAU,ERYD9,0,NETAU,TAUT,0,NT,.FALSE.,TAUQ(1,1,1),
     &             MSSQ(1,1,1),1,NQ,0,NKM,NKMMAX,IPRINT)
C
      DERYD = 0.001D0
C
      NEINTEG = 0
      NEA = NEINTEG + 1
      NEB = 0
C
      IFILA = IFILCBWF
      IFILB = IFILCBWF
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop 1   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      IEPATH = 1
      CPACHNGMAX = 0D0
C
      IF ( NETAU.EQ.0 ) THEN
Ccc      IF ( NETAU.NE.111 ) THEN
C
         WRTAUMQ = .TRUE.
C
         IE = 0
         DO IEA = 1,NEA
C
            IF ( IEA.LT.NEA ) THEN
               ERYDA = ETAB(IEA,IEPATH)
            ELSE
               ERYDA = EFERMI
            END IF
C
            DO IEB = 0,NEB
               IE = IE + 1
C
               IF ( IEB.EQ.0 ) THEN
                  ERYD = ERYDA
               ELSE IF ( IEB.EQ.1 ) THEN
                  ERYD = ERYDA - DERYD/2D0
               ELSE
                  ERYD = ERYDA + DERYD/2D0
               END IF
C
               ICPAFLAG = 0
               CPACHNG = 0.0D0
C
C ===================================== solve SS - differential equation
C
               CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,
     &                       GETIRRSOL,ERYD,P,IPRINT,TSST,MSST,SSST,
     &                       MEZZ,MEZJ,ORBPOL)
C
               CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
               CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,
     &                              ITCPA,ICPACONV,CONC,NOQ,ITOQ,PHASK,
     &                              IE,NTMAX,TSST,MSST,TSSQ,MSSQ,TAUQ)
C
               IF ( USENLCPA ) THEN
C
C=======================================================================
                  CALL NLCPAITER(CPAERR,P,IPRINTL,NKM,NQ,NKMQ,TAUQ,DROT,
     &                           IQORGQP,SYMUNITARY,MSSQ,MSST,ITOQ,IREL,
     &                           NSYM,NTMAX,NQMAX,NKMMAX,NKTABKNMAX,
     &                           MROTK,MUEHAT,OMEGAHAT,TAUHAT,NDIMCLU,
     &                           IND0QCLU,NQNLCPA,NKTABKN,KTABKN,
     &                           WKTABKN,NKN,NSYMACCEPTEDKN,
     &                           SYMACCEPTEDKN,KNNLCPA,RQNLCPA,IQCPA,
     &                           CPATOL,CPAMIX,NKNIRMU,SYM_IRR,
     &                           MAXNKNIRMU,LCPAMILLS,LNLCPASYMAVG,
     &                           LNLCPAAVG,ITCPA,ITCPAMAX,ICPAFLAG,
     &                           ICPACONV,PCFG,NCFG,NAQ,NA_Q,NU,NAQU,
     &                           I0_QU,I1_QU,I2_QU,I1_Q,I2_Q)
C
C   --------------------------------------------------------------------
C                               determine TAUT
C   --------------------------------------------------------------------
C
                  CALL RINIT(NTMAX,WGIT)
                  NDIMCLUSQ = NDIMCLU*NDIMCLU
                  CALL CINIT(NDIMCLUSQ,MMAT)
C
                  M = NDIMCLU
                  N = NKMQ(IQCPA)
C
C
C     ----------------------------  zero out TAUT blocks for nlcpa IT''s
                  DO IOCC = 1,2
                     IT = ITOQ(IOCC,IQCPA)
                     TAUT(1:N,1:N,IT) = C0
                  END DO
C
C   ---------------------------------------------------------------- CGF
                  DO ICFG = 1,NCFG
C
                     CALL NLCPAPROJ(ICFG,IQCPA,MSST,OMEGAHAT,TAUIMP,
     &                              W1HAT,NQNLCPA,N,IND0QCLU,ITOQ,
     &                              NDIMCLU,NQMAX,NTMAX,NKMMAX)
C
C        -------- get contribution to TAUT(IT) from cluster site IQCLU=1
C
                     IOCC = NLCPACONF(ICFG,NQNLCPA,1,1,2)
                     IT = ITOQ(IOCC,IQCPA)
                     TAUT(1:N,1:N,IT) = TAUT(1:N,1:N,IT) + PCFG(ICFG)
     &                                  *TAUIMP(1:N,1:N)
                     WGIT(IT) = WGIT(IT) + PCFG(ICFG)
C
                  END DO
C
C     --------------------------------------------------- normalise TAUT
C
                  DO IOCC = 1,2
                     IT = ITOQ(IOCC,IQCPA)
                     TAUT(1:N,1:N,IT) = TAUT(1:N,1:N,IT)/WGIT(IT)
                     IF ( IPRINTL.GT.3 ) THEN
                        WRITE (6,'(5x,"<SIGNLC>: TAUT for IT",i5)') IT
                        CALL CMATSTRUCT('TAUT      (KAP,MUE)',
     &                                  TAUT(1,1,IT),N,NKMMAX,IREL,IREL,
     &                                  0,1D-8,6)
                     END IF
                  END DO
C
                  CALL NLCTAUIO(IFILTAU,WRTAUMQ,IPRINT,IE,ERYD,TAUQ,
     &                          MSSQ,OMEGAHAT,MUEHAT,NQ,NKMQ,CPACHNG,
     &                          ICPAFLAG,IQCPA,NDIMCLU,NQMAX,NKMMAX)
C
                  IF ( WRMAT ) THEN
                     IWRI = 6
                     IWRI = 101
                     WRITE (IWRI,*) 'SIGNLC - ITER'
C
                     CALL NLCDUMPTAU(IE,ERYD,IWRI,MUEHAT,OMEGAHAT,
     &                               TAUHAT,IREL,N,NDIMCLU)
                     CALL DUMPTAU(IE,ERYD,6,MSST,MSSQ,TAUT,TAUQ)
                  END IF
C
               ELSE
C ======================================================================
C
                  CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,
     &                                 ITCPA,ICPACONV,CONC,NOQ,ITOQ,
     &                                 PHASK,IE,NTMAX,TSST,MSST,TSSQ,
     &                                 MSSQ,TAUQ)
C
                  IF ( ICPAFLAG.NE.0 ) THEN
                     NCPAFAIL = NCPAFAIL + 1
                     CPACHNGMAX = MAX(CPACHNGMAX,ABS(CPACHNG))
                  END IF
C
                  CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,
     &                         TAUT)
C
C ======================================================================
C
               END IF
C
               IF ( WRMAT .OR. IPRINT.GE.3 )
     &              CALL DUMPTAU(IE,ERYD,6,MSST,MSSQ,TAUT,TAUQ)
C
C-----------------------------------------------------------------------
               IF ( IE.EQ.NETAB(IEPATH) ) THEN
C
                  IF ( NCPAFAIL.NE.0 ) THEN
                     WRITE (OUT,99005) NCPAFAIL,CPATOL,CPACHNGMAX
                  ELSE IF ( NCPA.NE.0 ) THEN
                     WRITE (OUT,99006)
                  END IF
C
               END IF
C
            END DO
         END DO
      END IF
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      IF ( NETAU.EQ.0 ) CALL READTAU(IFILTAU,ERYD9,0,NETAU,TAUT,0,NT,
     &                               .FALSE.,TAUQ(1,1,1),MSSQ(1,1,1),1,
     &                               NQ,0,NKM,NKMMAX,IPRINT)
      WRTAUMQ = .FALSE.
C
C----------------------------------------------------------------------
C            allocate storage for  CHI  calculation
C----------------------------------------------------------------------
C
      NDIMCHI = NQNLCPA*NQNLCPA*NKMMAX*NKMMAX
C
      ALLOCATE (CHIHATZ(NDIMCHI,NDIMCHI,NZ12MAX),STAT=IA_ERR)
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop 2   START
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      NEPATH = 1
      IEPATH = 1
C
      ICPAFLAG = 0
      CPACHNG = 0.0D0
C
      IE = 0
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      DO IEA = 1,NEA
C
         IF ( IEA.LT.NEA ) THEN
            INTEG = .TRUE.
            ERYDA = ETAB(IEA,IEPATH)
            NEB = 2
         ELSE
            INTEG = .FALSE.
            ERYDA = EFERMI
            NEB = 1
         END IF
C
C ===================================== solve SS - differential equation
C
         CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILA,GETIRRSOL,ERYDA,
     &                 PA,IPRINT,TSST,MSSTA,SSST,MEZZ,MEZJ,ORBPOL)
         P = PA
C
         WRITE (6,*) 'setting OMEGAHATA and MUEHATA'
         WRITE (6,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
         WRITE (6,*) 'XXXX    READING TAUQA MSSQ  XXXXXXXXXXX'
         WRITE (6,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
         IF ( USENLCPA ) THEN
            CALL NLCTAUIO(IFILTAU,WRTAUMQ,IPRINT,IEA,ERYD,TAUQA,MSSQ,
     &                    OMEGAHATA,MUEHATA,NQ,NKMQ,CPACHNG,ICPAFLAG,
     &                    IQCPA,NDIMCLU,NQMAX,NKMMAX)
         ELSE
            DO IQ = 1,NQ
               CALL READTAU(IFILTAU,ERYD9,IE,NETAU,TAUT,0,NT,.FALSE.,
     &                      TAUQA(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,1,NKM,
     &                      NKMMAX,IPRINT)
               IF ( ABS(ERYDA-ERYD9).GT.1D-8 )
     &               STOP 'in <SIGMA>: >>> FILE 9'
               TAUQ(:,:,IQ) = TAUQA(:,:,IQ)
            END DO
C
            CALL PROJTAU(ICPAFLAG,CPACHNG,ERYDA,MSSTA,MSSQ,TAUQA,TAUTA)
C
            DO IT = 1,NT
               IQ = IQAT(1,IT)
               CALL GETDMAT(TAUQA(1,1,IQ),DMATTA(1,1,IT),DTILTA(1,1,IT),
     &                      DMAMC,NKM,MSSQ(1,1,IQ),MSSTA(1,1,IT),NKMMAX)
            END DO
         END IF
C
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         DO IEB = 1,NEB
            IE = IE + 1
C
            IF ( INTEG ) THEN
               ERYDB = ERYDA + (IEB-1.5D0)*DERYD
            ELSE
               ERYDB = EFERMI
            END IF
C
C ===================================== solve SS - differential equation
C
            CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILB,GETIRRSOL,
     &                    ERYDB,PB,IPRINT,TSST,MSSTB,SSST,MEZZ,MEZJ,
     &                    ORBPOL)
C
            IF ( .NOT.INTEG ) THEN
               TAUQB(:,:,1:NQ) = TAUQA(:,:,1:NQ)
               IF ( USENLCPA ) THEN
                  WRITE (6,*) 'setting OMEGAHATB and MUEHATB'
                  M = NDIMCLU
                  OMEGAHATB(1:M,1:M) = OMEGAHATA(1:M,1:M)
                  MUEHATB(1:M,1:M) = MUEHATA(1:M,1:M)
                  OMEGAHATB = OMEGAHATA
                  MUEHATB = MUEHATA
C
               END IF
            ELSE IF ( USENLCPA ) THEN
               CALL NLCTAUIO(IFILTAU,WRTAUMQ,IPRINT,IEA,ERYD,TAUQB,MSSQ,
     &                       OMEGAHATB,MUEHATB,NQ,NKMQ,CPACHNG,ICPAFLAG,
     &                       IQCPA,NDIMCLU,NQMAX,NKMMAX)
            ELSE
               DO IQ = 1,NQ
                  CALL READTAU(IFILTAU,ERYD9,IE,NETAU,TAUT,0,NT,.FALSE.,
     &                         TAUQB(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,1,NKM,
     &                         NKMMAX,IPRINT)
                  IF ( ABS(ERYDB-ERYD9).GT.1D-8 )
     &                  STOP 'in <SIGMA>: >>> FILE 9'
               END DO
            END IF
C
            CALL PROJTAU(ICPAFLAG,CPACHNG,ERYDB,MSSTB,MSSQ,TAUQB,TAUTB)
C
            DO IT = 1,NT
               IQ = IQAT(1,IT)
               CALL GETDMAT(TAUQB(1,1,IQ),DMATTB(1,1,IT),DTILTB(1,1,IT),
     &                      DMAMC,NKM,MSSQ(1,1,IQ),MSSTB(1,1,IT),NKMMAX)
            END DO
C
C ======================================================================
C
            CALL SIGNLCKLOOPS(PA,ERYDA,TAUQ,TAUQZ,DROT,IQORGQP,
     &                        SYMUNITARY,MSSQ,NSYM,CHIHATZ,IPRINTL,
     &                        MUEHATA,OMEGAHATA,TAUHAT,EIKNRIJ,KTABKN,
     &                        WKTABKN,NKTABKN,SYMACCEPTEDKN,
     &                        NSYMACCEPTEDKN,LNLCPAAVG,SYM_IRR,IQCPA,
     &                        NQNLCPA,NDIMCLU,IND0QCLU,NKN,NKNIRMU,
     &                        MAXNKNIRMU,NKTABKNMAX,NDIMCHI)
C
            IF ( WRMAT ) THEN
               IWRI = 6
               N = NKMQ(IQCPA)
               IWRI = 102
               WRITE (IWRI,*) 'SIGNLCKLOOPS'
               CALL NLCDUMPTAU(IE,ERYDA,IWRI,MUEHATA,OMEGAHATA,TAUHAT,
     &                         IREL,N,NDIMCLU)
            END IF
C
            IF ( NQ.GT.1 ) THEN
               OPEN (80,FILE='zzzzzz_chiz_NLCPA')
               DO K = 1,NZ12MAX
                  DO J = 1,NKMMAX*NKMMAX*NQMAX
                     DO I = 1,NKMMAX*NKMMAX*NQMAX
                        WRITE (80,'(3i5,2e20.12)') I,J,K,CHIZ(I,J,K)
                     END DO
                  END DO
               END DO
               DO K = 1,NQ
                  DO J = 1,NKMMAX
                     DO I = 1,NKMMAX
                        WRITE (80,'(3i5,2e20.12)') I,J,K,TAUQ(I,J,K)
                     END DO
                  END DO
               END DO
            END IF
C
C ======================================================================
C
            IF ( WRMAT .OR. IPRINT.GE.3 )
     &           CALL DUMPTAU(IE,ERYDB,6,MSSTB,MSSQ,TAUTB,TAUQB)
C
C ======================================================================
C
            CALL GILNLCME(IQCPA,NQNLCPA,NDIMCLU,IND0QCLU,NCFG,PCFG,
     &                    W1HAT,MSSTA,MSSTB,MUEHATA,MUEHATB,OMEGAHATA,
     &                    OMEGAHATB,MAQQAB,MAQQBA,MBQQAB,MBQQBA,MCQQAB,
     &                    MCQQBA,MDQQAB,MDQQBA,MTAB,MTBA)
C
            IPRINT = 3
C
            IF ( USENLCPA ) THEN
               CALL SIGNLC0(INTEG,MTAB,MTBA,TAUTA,TAUTB,CONC,NOQ,ITOQ,
     &                      ALF0Q,IPRINT,NCFG,PCFG,IQCPA,NQNLCPA,W1HAT,
     &                      MSSTA,MSSTB,OMEGAHATA,OMEGAHATB,NDIMCLU,
     &                      IND0QCLU,NKMQ,NKMMAX,NTMAX,NQMAX)
C
            ELSE
C
               CALL SIG0(INTEG,MTAB,MTBA,MIRRTAB,TAUTA,TAUTB,ALF0Q,
     &                   SIG_MODE)
C
            END IF
C
            CALL SIGNLC1(CHIHATZ,MAQQAB,MBQQAB,MCQQAB,MDQQAB,ALF1QQ_NV,
     &                   'N',IPRINT,NQNLCPA,NKM,NKMMAX,NQMAX,NZ12MAX,
     &                   NDIMCHI)
C
            IF ( VERTEX ) THEN
C
               CALL LINRESP_VERTEX_NLC(CHIHATZ,MUEHATA,MUEHATB,
     &                                 OMEGAHATA,OMEGAHATB,MSSTA,MSSTB,
     &                                 NQNLCPA,IND0QCLU,ITOQ,NCFG,PCFG,
     &                                 NZ12,IQCPA,NKM,NKMQ,NKMMAX,NQMAX,
     &                                 NZ12MAX,NDIMCLU,NTMAX,NDIMCHI)
C
               CALL SIGNLC1(CHIHATZ,MAQQAB,MBQQAB,MCQQAB,MDQQAB,
     &                      ALF1QQ_VC,'V',IPRINT,NQNLCPA,NKM,NKMMAX,
     &                      NQMAX,NZ12MAX,NDIMCHI)
C
            ELSE
C
               CALL ZCOPY(3*3*NQMAX*NQMAX,ALF1QQ_NV,1,ALF1QQ_VC,1)
C
            END IF
C
            CALL GILSUM(ALF0Q,ALF1QQ_NV,ALF1QQ_VC,ALF0QO,ALF1QOQ_NV,
     &                  ALF1QOQ_VC,1,1,0D0)
C
         END DO
CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      END DO
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      CALL CPU_TIME(TIME)
      WRITE (6,99004) TIME - TIME0
C
      DEALLOCATE (MTAB,MTBA,DTILTA,DMATTA,DTILTB,DMATTB,DMAMC)
      DEALLOCATE (MAQAB,MAQBA,MBQAB,MBQBA,MCQAB,MCQBA,MDQAB,MDQBA)
      DEALLOCATE (TAUQA,TAUQB,TAUTA,TAUTB)
C
      WRITE (6,*) '          SIGMA-job done      FEIERABEND '
C
      STOP
C
C   ====================================================================
99001 FORMAT (10X,'the standard CPA will be used ',/)
99002 FORMAT (10X,'the NLCPA will be used ',/,10X,
     &        'site with disorder  IQCPA       ',I10,/,10X,
     &        'NLCPA-level (cluster size)      ',I10,/,10X,
     &        'mixing parameter                ',F10.5,/,10X,
     &        'SRO - parameter ALPHA           ',F10.5,/)
99003 FORMAT (///,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*            ****   ***   ****   *     *    ***              *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  **   *    *   *             *'
     &  ,/,10X,
     &  '*           *        *   *       * * * *  *     *            *'
     &  ,/,10X,
     &  '*            ****    *   *  ***  *  *  *  *******            *'
     &  ,/,10X,
     &  '*                *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*           *    *   *   *    *  *     *  *     *            *'
     &  ,/,10X,
     &  '*            ****   ***   ****   *     *  *     *            *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99004 FORMAT (/,5X,'execution time for <SIGNLC>:',F14.3,' secs',/)
99005 FORMAT (/,1X,79('*'),/,10X,'CPA not converged for',I3,
     &        ' energies:',/,10X,'tolerance for CPA-cycle:',F15.7,/,10X,
     &        'maximum deviation:      ',F15.7,/,1X,79('*'),/)
99006 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
      END
