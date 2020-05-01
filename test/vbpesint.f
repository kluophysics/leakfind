C*==vbpesint.f    processed by SPAG 6.70Rc at 15:40 on 19 Dec 2016
      SUBROUTINE VBPESINT(MGEO,NE,MEZZ,MEZJ,NKMINI,NKMFIN,LLAM,MJLAM,
     &                    MSPN,IFILINI,IFILFIN,IFILTAUNN,MEFORM,IPRINT,
     &                    CALCINT,GETIRRSOL,MODE,NPOL,ORBPOLFIN,
     &                    USETAUNN,KGEO,NSPIN,TMBAR,MBAR,MREG,MIRR,
     &                    MEIRR,CCYY,ETABFIN,MEREG,SA,ST,TIE,TM,MEBAR,
     &                    NPOLMAX,TSST,MSST,SSST,TSSQ,MSSQ,TAUT,TSSTFIN,
     &                    MSSTFIN,TSSQFIN,MSSQFIN)
C   ********************************************************************
C   *                                                                  *
C   *   VBPES   angle integrated valence band photo emission spectrum  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ENERGY,ONLY:EFERMI,NEMAX,ETAB
      USE MOD_RMESH,ONLY:JRNSMIN,NRMAX,JRCRI,JRNS1,FULLPOT
      USE MOD_CPA,ONLY:NCPA
      USE MOD_FILES,ONLY:IOTMP,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CALCMODE,ONLY:IREL,DMFT,ORBPOL
      USE MOD_TYPES,ONLY:NT,NCPLWFMAX,NLMFPMAX,NTMAX,VAMEF,VAMEG,CTL,
     &    CONC,NAT,IMT,NLT,NKM_T
      USE MOD_SITES,ONLY:NQ,NQMAX
      USE MOD_ANGMOM,ONLY:CGC,NMEMAX,NKMMAX,NL,NLQ,NKMQ,MZAZB,MZBZA,
     &    MIRR_2,MIRR_3,MIRR_4,WKM1,WKM2,NKM
      USE MOD_DMFT_LDAU,ONLY:SEVT,SEVNST,SEBT,DMFTSIG,DMFTSIGMA,SEBNST,
     &    KSELF
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--VBPESINT30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='VBPESINT')
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL,USETAUNN
      INTEGER IFILFIN,IFILINI,IFILTAUNN,IPRINT,KGEO,MODE,NE,NKMFIN,
     &        NKMINI,NPOL,NPOLMAX
      CHARACTER*3 MEFORM
      CHARACTER*10 ORBPOLFIN
      COMPLEX*16 CCYY(NKMMAX,2,NKMMAX,2),ETABFIN(NEMAX,2),
     &           MBAR(NKMMAX,NKMMAX,3),MEBAR(NKMMAX,NKMMAX,3,NT),
     &           MEIRR(NKMMAX,NKMMAX,3,NT),MEREG(NKMMAX,NKMMAX,3,NT),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),MGEO(3,3),
     &           MIRR(NKMMAX,NKMMAX,3,3),MREG(NKMMAX,NKMMAX,3),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSSQFIN(NKMMAX,NKMMAX,NQMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),MSSTFIN(NKMMAX,NKMMAX,NTMAX),
     &           SA(2,2,NPOLMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           ST(2,2,NPOLMAX,NTMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TIE(NKMMAX,NKMMAX),TM(NKMMAX,NKMMAX),
     &           TMBAR(NKMMAX,NKMMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSSQFIN(NKMMAX,NKMMAX,NQMAX),TSST(NKMMAX,NKMMAX,NTMAX),
     &           TSSTFIN(NKMMAX,NKMMAX,NTMAX)
      INTEGER LLAM(NKMMAX)
      REAL*8 MJLAM(NKMMAX),MSPN(2),NSPIN(3)
C
C Local variables
C
      COMPLEX*16 CSUM1,CSUM2,CYTIT,CYTMTM,EFIN,EFIN9,EINI,ERYD,P,PINI,
     &           TAUTI(NKMMAX,NKMMAX,NTMAX)
      INTEGER I,IA_ERR,IBLK,IBLKTOP,IE,IK2,IPOL,IPOLREV,IPOT,IQ,IT,
     &        ITMP(:),J,J1,J1TOP,LAM,LAM1,LAM2,LAMP,LAMPP,LAMPPP,MC,MS,
     &        MSP,N,NE9,NKMFINT,NKMINIT
      REAL*8 POL(2,NPOLMAX),POLNT(2,NPOLMAX,NTMAX),TOT(NPOLMAX),
     &       TOTNT(NPOLMAX,NTMAX)
C
C*** End of declarations rewritten by SPAG
C
      DATA ERYD/(999999D0,999999D0)/
C
      ALLOCATABLE ITMP
      ALLOCATE (ITMP(NTMAX))
C
      ALLOCATE (MZAZB(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MZBZA(NKMMAX,NKMMAX,NPOLMAX))
      ALLOCATE (MIRR_2(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_3(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX))
      ALLOCATE (MIRR_4(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MZAZB')
C
C=======================================================================
C                         set up prefactor CCYY =  CGC * CGC
C=======================================================================
C
      DO MS = 1,2
         DO MSP = 1,2
            DO LAMPP = 1,NKMFIN
               DO LAM = 1,NKMFIN
                  IF ( (LLAM(LAMPP).EQ.LLAM(LAM)) .AND. 
     &                 (NINT(MJLAM(LAMPP)-MSPN(MSP))
     &                 .EQ.NINT(MJLAM(LAM)-MSPN(MS))) ) THEN
                     CCYY(LAM,MS,LAMPP,MSP) = CGC(LAM,MS)*CGC(LAMPP,MSP)
                  ELSE
                     CCYY(LAM,MS,LAMPP,MSP) = C0
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C ======================================================================
C ========================================= ENERGY LOOP === START ======
C
      WRITE (40,*) NE,EFERMI
      WRITE (40,*) KGEO
      WRITE (40,*) 1
      WRITE (40,*) NT
      NKM_T(1:NT) = NKM
C
      DO IE = 1,NE
         EINI = ETAB(IE,1)
         EFIN = ETABFIN(IE,1)
C
         WRITE (7,'(I5,'' E(ini) = '',2F10.5)') IE,EINI
         WRITE (40,'(I5,'' E(ini) = '',F10.5)') IE,DREAL(EINI)
         WRITE (6,'(I5,'' E(ini) = '',F10.5,10X,'' E(fin) = '',F10.5)')
     &          IE,DREAL(EINI),DREAL(EFIN)
C
C-----------------------------------------------------------------------
C----------------------------------- Calling ssite to prepare t and WF-s
C
         IF ( DMFT ) THEN
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NT)
     &         = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NT,IE)
            DO IT = 1,NT
               IF ( KSELF(IT).EQ.1 ) WRITE (6,
     &              '(/,19X,''SELFENERGY applied for IT:'',I3)') IT
            END DO
         END IF
         IF ( FULLPOT ) THEN
            IF ( ORBPOL(1:5).EQ.'SIGMA' ) THEN
               CALL DMFT_READSELFENE(EINI,EFERMI,KSELF,SEVT,SEBT,SEVNST,
     &                               SEBNST,IMT,JRNS1,JRCRI,JRNSMIN,
     &                               NTMAX,NLMFPMAX,NRMAX,FULLPOT)
               DO IT = 1,NT
                  IF ( KSELF(IT).EQ.1 .AND. IE.EQ.1 ) WRITE (6,
     &                 '(/,19X,''SELFENERGY applied for IT:'',I3)') IT
               END DO
            END IF
C
            DO IQ = 1,NQ
               NLQ(IQ) = NLQ(IQ) - 1
               NL = MAX(NL,NLQ(IQ))
               NKMQ(IQ) = 2*NLQ(IQ)**2
            END DO
            CLOSE (IOTMP)
C
            CALL FPCOUPL(0,NKMINI,NL-1)
C
            REWIND (IOTMP)
C
            READ (IOTMP) MC,IBLKTOP,J1TOP
C
            IF ( MC.GT.NCPLWFMAX ) WRITE (6,*)
     &            'VBPESINT: MC > NCPLWFMAX  A:',MC,NCPLWFMAX
C
            READ (IOTMP) ((((((VAMEG(I,J,IPOT,J1,IBLK,IT),I=1,MC),J=1,MC
     &                   ),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,IBLKTOP),
     &                   IT=1,NT)
C
            IF ( IREL.GE.3 ) READ (IOTMP)
     &                             ((((((VAMEF(I,J,IPOT,J1,IBLK,IT),I=1,
     &                             MC),J=1,MC),IPOT=1,NLMFPMAX),J1=1,
     &                             J1TOP),IBLK=1,IBLKTOP),IT=1,NT)
            CLOSE (IOTMP)
            ITMP(1:NTMAX) = NKM_T(1:NTMAX)
            NKM_T(1:NT) = NKMINI
            CALL FPSSITE(1,1,IFILINI,GETIRRSOL,EINI,P,IPRINT,TSST,MSST,
     &                   SSST,MEZZ,MEZJ,ORBPOL)
C
            NKM_T(1:NT) = ITMP(1:NTMAX)
C
C----------------------------------------------------------- final state
C
            DO IQ = 1,NQ
               NLQ(IQ) = NLQ(IQ) + 1
               NL = MAX(NL,NLQ(IQ))
               NKMQ(IQ) = 2*NLQ(IQ)**2
            END DO
C
            CALL FPCOUPL(IPRINT,NKMFIN,NL)
C
            REWIND (IOTMP)
C
            READ (IOTMP) MC,IBLKTOP,J1TOP
C
            READ (IOTMP) ((((((VAMEG(I,J,IPOT,J1,IBLK,IT),I=1,MC),J=1,MC
     &                   ),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,IBLKTOP),
     &                   IT=1,NT)
C
            IF ( IREL.GE.3 ) READ (IOTMP)
     &                             ((((((VAMEF(I,J,IPOT,J1,IBLK,IT),I=1,
     &                             MC),J=1,MC),IPOT=1,NLMFPMAX),J1=1,
     &                             J1TOP),IBLK=1,IBLKTOP),IT=1,NT)
            CLOSE (IOTMP)
C
            ITMP(1:NTMAX) = NKM_T(1:NTMAX)
            NKM_T(1:NT) = NKMFIN
C
            CALL FPSSITE(1,1,IFILFIN,GETIRRSOL,EFIN,P,IPRINT,TSSTFIN,
     &                   MSST,SSST,MEZZ,MEZJ,ORBPOLFIN)
            NKM_T(1:NT) = ITMP(1:NTMAX)
C
         ELSE
C
            CALL SSITE(1,1,IFILINI,CALCINT,GETIRRSOL,EINI,PINI,IPRINT,
     &                 NKMFIN,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C----------------------------------------------------------- final state
C
            CALL SSITE(1,1,IFILFIN,CALCINT,GETIRRSOL,EFIN,P,IPRINT,
     &                 NKMFIN,TSSTFIN,MSST,SSST,MEZZ,MEZJ,ORBPOLFIN)
C
C --------------------------- initialize MSSQ in case of ordered systems
C
         END IF
C----------------------------------- Calling ssite to prepare t and WF-s
C-----------------------------------------------------------------------
C
C ------------------------------------------------------ initialize MSSQ
C
         IF ( NCPA.EQ.0 ) THEN
C
            CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
            CALL MSSINIT(TSSTFIN,MSSTFIN,TSSQFIN,MSSQFIN)
C
         END IF
C
C=======================================================================
C                            matrix elements
C=======================================================================
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         DO IT = 1,NT
            MREG(:,:,:) = C0
            MBAR(:,:,:) = C0
            MIRR(:,:,:,:) = C0
C
C        MBAR = < Z(ini) | H_lam | Z(fin) >
C
            CALL ME_DRIVE(1,NKMINI,IFILINI,EINI,1,NKMFIN,IFILFIN,EFIN,
     &                    .FALSE.,IT,MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4,
     &                    CTL(IT,1),2,MEFORM,0)
C
            MBAR(:,:,1) = -MZAZB(:,:,3)
            MBAR(:,:,2) = +MZAZB(:,:,1)
            MBAR(:,:,3) = +MZAZB(:,:,2)
C
            MEBAR(:,:,1,IT) = -MZAZB(:,:,3)
            MEBAR(:,:,2,IT) = +MZAZB(:,:,1)
            MEBAR(:,:,3,IT) = +MZAZB(:,:,2)
C
C        MREG = < Z(fin)* | H_lam | Z(ini) >
C        MIRR = < Z(fin)* | H_lam | Z(ini) J(ini) | H_lam+ | Z(fin) >
C
            CALL ME_DRIVE(1,NKMFIN,IFILFIN,EFIN,1,NKMINI,IFILINI,EINI,
     &                    .TRUE.,IT,MZAZB,MZBZA,MIRR_2,MIRR_3,MIRR_4,
     &                    CTL(IT,1),3,MEFORM,0)
C
            MREG(:,:,1) = -MZAZB(:,:,3)
            MREG(:,:,2) = +MZAZB(:,:,1)
            MREG(:,:,3) = +MZAZB(:,:,2)
C
            MBAR(:,:,1) = MEBAR(:,:,1,IT) - MZBZA(:,:,3)
            MBAR(:,:,2) = MEBAR(:,:,2,IT) + MZBZA(:,:,1)
            MBAR(:,:,3) = MEBAR(:,:,3,IT) + MZBZA(:,:,2)
C
            MIRR(:,:,2,1) = -MIRR_3(:,:,3,1)
            MIRR(:,:,1,2) = -MIRR_3(:,:,1,3)
            MIRR(:,:,3,3) = +MIRR_3(:,:,2,2)
C
C --------------------------------------------------------------- IPOL -
C              index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C               M(X) =   [  M(+) + M(-) ] / SQRT(2)
C               M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
            IF ( KGEO.EQ.0 ) THEN
C
               DO IPOL = 1,NPOL
                  IF ( IPOL.EQ.3 ) THEN
                     IPOLREV = IPOL
                  ELSE
                     IPOLREV = 3 - IPOL
                  END IF
C
                  N = NKMMAX*NKMMAX
                  CALL ZCOPY(N,MREG(1,1,IPOL),1,MEREG(1,1,IPOL,IT),1)
                  CALL ZCOPY(N,MBAR(1,1,IPOLREV),1,MEBAR(1,1,IPOL,IT),1)
C
                  N = NKMMAX*NKMMAX
                  CALL ZCOPY(N,MIRR(1,1,IPOL,IPOLREV),1,
     &                       MEIRR(1,1,IPOL,IT),1)
C
               END DO
C
            ELSE
C
               CALL VBPESMETR(MBAR,MEBAR(1,1,1,IT),.TRUE.,.FALSE.,MIRR,
     &                        MEIRR(1,1,1,IT),MGEO,NKMFIN,NKMMAX)
C
               CALL VBPESMETR(MREG,MEREG(1,1,1,IT),.FALSE.,.TRUE.,MIRR,
     &                        MEIRR(1,1,1,IT),MGEO,NKMFIN,NKMMAX)
C
            END IF
C
         END DO
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C-----------------------------------                     matrix elements
C-----------------------------------------------------------------------
C
C
C
C=======================================================================
C                          Calculation of spectra
C=======================================================================
C
         DO IK2 = 1,1
C
            CALL CINIT(2*2*NPOLMAX,SA)
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
            DO IT = 1,NT
               NKMINI = 2*(NLT(IT)-1)**2
               NKMINIT = 2*(NLT(IT)-1)**2
               NKMFINT = 2*NLT(IT)**2
C
CKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C-----------------------------------------------------------------------
C-----------------------------------Reading TAU for ini and final states
C
               CALL CINIT(NKMMAX**2,TAUT(1,1,IT))
               CALL CINIT(NKMMAX**2,TAUTI(1,1,IT))
C
               CALL READTAU(9,ERYD,IE,NE,TAUTI(1,1,IT),IT,NT,.TRUE.,
     &                      WKM1,WKM2,0,NQ,1,NKMINI,NKMMAX,IPRINT)
C
               DO I = 1,NKMMAX
                  DO J = 1,NKMMAX
                     TAUT(I,J,IT) = TAUTI(I,J,IT)
                  END DO
               END DO
C
C     --------------------------------------------------- init t-ss(FIN)
C
               IF ( USETAUNN ) THEN
                  CALL CINIT(NKMMAX*NKMMAX*NTMAX,TSSTFIN)
C
                  CALL READTAU(IFILTAUNN,EFIN9,IE,NE9,TSSTFIN(1,1,IT),
     &                         IT,NT,.TRUE.,WKM1,WKM2,0,NQ,1,NKMFIN,
     &                         NKMMAX,IPRINT)
C
                  IF ( ABS(EFIN9-EFIN).GT.1D-3 ) THEN
                     WRITE (6,*) ' WARNING in <VBPESINT>'
                     WRITE (6,*) ' E_initial    ',ERYD
                     WRITE (6,*) ' E_final      ',EFIN
                     WRITE (6,*) ' E_final(9)   ',EFIN9
                     WRITE (6,*) ' DELTA        ',EFIN9 - EFIN
                  END IF
C
               END IF
C
C IPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIP
C
               DO IPOL = 1,NPOL
C
                  IF ( IPOL.EQ.3 ) THEN
                     IPOLREV = IPOL
                  ELSE
                     IPOLREV = 3 - IPOL
                  END IF
C
C-----------------------------------------------------------------------
C                     set up REGULAR terms
C-----------------------------------------------------------------------
C
                  CALL CINIT(NKMMAX*NKMMAX,TM)
                  CALL CINIT(NKMMAX*NKMMAX,TMBAR)
C
                  DO LAM1 = 1,NKMINIT
                     DO LAM = 1,NKMFINT
                        DO LAMP = 1,NKMFINT
C
                           TM(LAM,LAM1) = TM(LAM,LAM1)
     &                        + DCONJG(TSSTFIN(LAMP,LAM,IT))
     &                        *MEREG(LAMP,LAM1,IPOL,IT)
C
                           TMBAR(LAM,LAM1) = TMBAR(LAM,LAM1)
     &                        + TSSTFIN(LAMP,LAM,IT)
     &                        *MEBAR(LAM1,LAMP,IPOL,IT)*0.5D0
C
                        END DO
                     END DO
                  END DO
C
                  IF ( MODE.EQ.3 ) THEN
                     CALL CINIT(NKMMAX*NKMMAX,TM)
                     CALL CINIT(NKMMAX*NKMMAX,TMBAR)
                  END IF
C
C-----------------------------------------------------------------------
C                     set up IRREGULAR term
C-----------------------------------------------------------------------
C
                  IF ( MODE.GT.1 ) THEN
C
                     WKM1(1:NKMFINT,1:NKMFINT)
     &                  = MEIRR(1:NKMFINT,1:NKMFINT,IPOL,IT)
C
                     DO LAM = 1,NKMFINT
                        DO LAMPP = 1,NKMFINT
                           CSUM2 = C0
                           DO LAMPPP = 1,NKMFINT
                              CSUM1 = C0
C
                              DO LAMP = 1,NKMFINT
                                 CSUM1 = CSUM1 + 
     &                              DCONJG(TSSTFIN(LAMP,LAM,IT))
     &                              *WKM1(LAMP,LAMPPP)
                              END DO
                              CSUM2 = CSUM2 + 
     &                                CSUM1*TSSTFIN(LAMPPP,LAMPP,IT)
C
                           END DO
                           TIE(LAM,LAMPP) = CSUM2
                        END DO
                     END DO
C
                  END IF
C
C ----------------------------------------------------------------------
C                          SUM UP
C ----------------------------------------------------------------------
C
                  DO MS = 1,2
                     DO MSP = 1,2
C
                        ST(MS,MSP,IPOL,IT) = C0
C
                        DO LAM2 = 1,NKMINIT
                           DO LAM1 = 1,NKMINIT
C
                              CYTMTM = C0
                              DO LAMPP = 1,NKMFINT
                                 DO LAM = 1,NKMFINT
C
                                    CYTMTM = CYTMTM + 
     &                                 CCYY(LAM,MS,LAMPP,MSP)
     &                                 *TM(LAM,LAM1)*TMBAR(LAMPP,LAM2)
C
                                 END DO
                              END DO
C
                              ST(MS,MSP,IPOL,IT) = ST(MS,MSP,IPOL,IT)
     &                           - TAUT(LAM1,LAM2,IT)*CYTMTM
C
                           END DO
                        END DO
C
C ----------------------------------------------------------------------
C
                        IF ( MODE.GT.1 ) THEN
                           CYTIT = C0
                           DO LAM = 1,NKMFINT
                              DO LAMPP = 1,NKMFINT
                                 CYTIT = CYTIT + CCYY(LAM,MS,LAMPP,MSP)
     &                              *TIE(LAM,LAMPP)
                              END DO
C
                           END DO
                           ST(MS,MSP,IPOL,IT) = ST(MS,MSP,IPOL,IT)
     &                        + CYTIT
                        END IF
C
C ----------------------------------------------------------------------
C
                        SA(MS,MSP,IPOL) = SA(MS,MSP,IPOL) + CONC(IT)
     &                     *NAT(IT)*ST(MS,MSP,IPOL,IT)
C
                     END DO
                  END DO
               END DO
C IPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIPOLIP
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
            END DO
C
            CALL FLUSH(6)
C
            CALL VBPESINTTR_0(SA,TOT,POL,NSPIN)
C
            WRITE (7,99001) (TOT(IPOL),TOT(IPOL),IPOL=1,NPOL),
     &                      (POL(1,IPOL),POL(2,IPOL),IPOL=1,NPOL)
            WRITE (40,99001) (TOT(IPOL),TOT(IPOL),IPOL=1,NPOL),
     &                       (POL(1,IPOL),POL(2,IPOL),IPOL=1,NPOL)
C
            DO IT = 1,NT
               CALL VBPESINTTR_0(ST(1,1,1,IT),TOTNT(1,IT),POLNT(1,1,IT),
     &                           NSPIN)
            END DO
C
            DO IT = 1,NT
               WRITE (7,99002) IT,(TOTNT(IPOL,IT),TOTNT(IPOL,IT),IPOL=1,
     &                         NPOL),
     &                         (POLNT(1,IPOL,IT),POLNT(2,IPOL,IT),IPOL=1
     &                         ,NPOL)
               WRITE (40,99002) IT,
     &                          (TOTNT(IPOL,IT),TOTNT(IPOL,IT),IPOL=1,
     &                          NPOL),
     &                          (POLNT(1,IPOL,IT),POLNT(2,IPOL,IT),
     &                          IPOL=1,NPOL)
            END DO
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( WRBUILDBOT .AND. IE.LE.MIN(3,NE) ) THEN
               WRITE (IFILBUILDBOT,99003) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                (TOT(IPOL),IPOL=1,NPOL),
     &                (POL(1,IPOL),POL(2,IPOL),IPOL=1,NPOL)
               DO IT = 1,NT
                  WRITE (IFILBUILDBOT,99004)
     &                   ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                   (TOTNT(IPOL,IT),IPOL=1,NPOL),
     &                   (POLNT(1,IPOL,IT),POLNT(2,IPOL,IT),IPOL=1,NPOL)
               END DO
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
CKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
         END DO
C
         CALL FLUSH(7)
         CALL FLUSH(40)
C
      END DO
C =========================================== ENERGY LOOP === END ======
      CLOSE (40)
C
      RETURN
C
C-----------------------------------------------------------------------
99001 FORMAT (' SA     SP ',40E12.5)
99002 FORMAT (' ST(',I2,') SP ',100E12.5)
C------------------------------------------------------------for testing
Cc99001 FORMAT (' SA     SP ',40E25.14)
Cc99002 FORMAT (' ST(',I2,') SP ',100E25.14)
99003 FORMAT ('# BUILDBOT: ',A,':  total   XPS int. and pol.',/,
     &        (1PE22.14))
99004 FORMAT ('# BUILDBOT: ',A,':  partial XPS int. and pol.',
     &        '  for IT =',I5,/,(1PE22.14))
      END
