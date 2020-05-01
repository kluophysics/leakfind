C*==scfoppot.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFOPPOT(EFCORRECT,SHFTEF,IECURR,OBS_T,TAUT,RHOOPC,
     &                    RHOOPO,QLZ,AOPT)
C   ********************************************************************
C   *                                                                  *
C   * calculate spin-resolved charge and orbital momentum densities    *
C   * used in any schemes dealing with the ORBITAL POLARISATION        *
C   *                                                                  *
C   * ORBPOL = BROOKS     M.S.S. Brooks' OP-scheme                     *
C   * l-channel (d or f) selected by LOPT                              *
C   *                                                    x              *
C   * AMEOP*  allocated with dimension  NKMPMAX                        *
C   *         with AMEOP*(IKM) = 0 for IKM > NKM                       *
C   *         this allows NKM < IKM <= NKMPMAX for the small component *
C   *                                                                  *
C   ********************************************************************
C NOTE: IMKM(KA) and IMKM(KB) allowed to get > NKM  (but < NKMPMAX)
C       AMEOP* allocated with dimen. NKMPMAX with AMEOP(I)=0 for I > NKM
C
      USE MOD_TYPES,ONLY:NT,NTMAX,IMT,TXT_T,NCPLWFMAX,IKMCPLWF,LOPT
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,ISMT,IOMT,AMEOPC,AMEOPO,NCPLWF,
     &    NKMMAX,NL,NLMAX
      USE MOD_CONSTANTS,ONLY:PI,RY_EV
      USE MOD_RMESH,ONLY:NRMAX,R,R2DRDI,JRWS
      USE MOD_ENERGY,ONLY:WETAB
      USE MOD_CALCMODE,ONLY:IREL,ORBPOL
      USE MOD_FILES,ONLY:IFILCBWF
      IMPLICIT NONE
C*--SCFOPPOT29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFOPPOT')
      INTEGER NLPREMAX
      PARAMETER (NLPREMAX=6)
C
C Dummy arguments
C
      LOGICAL EFCORRECT
      INTEGER IECURR
      REAL*8 SHFTEF
      REAL*8 AOPT(NRMAX,2,NTMAX),OBS_T(0:3,NOBSMAX,NTMAX),QLZ(NTMAX,2),
     &       RHOOPC(NRMAX,NTMAX,2),RHOOPO(NRMAX,NTMAX,2)
      COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BRACAH,FCLMB,MJ,MJMAX,MJMIN,PRE(NLPREMAX),QNORM,QQC,QQCMS,
     &       QQO,QQOMS,QQS,RHORSQ,RINT(NRMAX),RPW(NRMAX,2*NLMAX),
     &       SG(NRMAX),SL(NRMAX),TG(NRMAX),TL(NRMAX)
      INTEGER I,IA_ERR,IEPATH,IKMCB(2),IL,IM,IMKMCB(2),IMLAST,IR,IRTOP,
     &        IT,K1,K2,KA,KAP1,KAP2,KB,L,LAM1,LAM2,LBOT,LF,LFMAX,LMAX,
     &        LTOP,MS,MUE,MUETOP,NSOL
      INTEGER IKAPMUE
      COMPLEX*16 JF(:,:,:),JG(:,:,:),WCF,WCG,WE,WOF,WOG,WT,ZF(:,:,:),
     &           ZFJF,ZFZF,ZG(:,:,:),ZGJG,ZGZG
      LOGICAL OPF(NTMAX)
      CHARACTER*2 RACPAR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JF,JG,ZF,ZG
C
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IEPATH = 1
      IF ( EFCORRECT ) THEN
         WE = SHFTEF
      ELSE
         WE = WETAB(IECURR,IEPATH)
      END IF
C
      LBOT = 0
      LTOP = NL - 1
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         OPF(IT) = .FALSE.
         IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
            IF ( LOPT(IT).EQ.2 ) THEN
               LBOT = 2
               LTOP = 2
               OPF(IT) = .TRUE.
            ELSE IF ( LOPT(IT).EQ.3 ) THEN
               LBOT = 3
               LTOP = 3
               OPF(IT) = .TRUE.
            END IF
         END IF
C
         IF ( OPF(IT) ) THEN
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C
            LMAX = NL - 1
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO L = 0,LMAX
               IL = L + 1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               KAP1 = -L - 1
               KAP2 = L
               IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
               IF ( IREL.GT.1 ) THEN
                  MJMAX = DBLE(L) + 0.5D0
                  MUETOP = 2*L + 2
               ELSE
                  MJMAX = DBLE(L)
                  MUETOP = 2*L + 1
               END IF
               MJMIN = -MJMAX
               MJ = MJMIN - 1D0
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
               DO MUE = 1,MUETOP
                  MJ = MJ + 1D0
C
                  IF ( IREL.LE.1 ) THEN
                     NSOL = 1
C           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
                  ELSE IF ( ABS(MJ).GT.DBLE(L) ) THEN
                     NSOL = 1
                  ELSE
                     NSOL = 2
                  END IF
C-----------------------------------------------------------------------
                  IKMCB(1) = IKAPMUE(KAP1,NINT(MJ-0.5D0))
                  IKMCB(2) = IKAPMUE(KAP2,NINT(MJ-0.5D0))
                  IMKMCB(1) = IKAPMUE(-KAP1,NINT(MJ-0.5D0))
                  IMKMCB(2) = IKAPMUE(-KAP2,NINT(MJ-0.5D0))
C-----------------------------------------------------------------------
C
                  DO K1 = 1,NSOL
                     LAM1 = IKMCB(K1)
                     DO K2 = 1,NSOL
                        LAM2 = IKMCB(K2)
C
                        WT = -TAUT(LAM2,LAM1,IT)/PI
C
C-----------------------------------------------------------------------
                        IF ( L.GE.LBOT .AND. L.LE.LTOP ) THEN
                           DO MS = 1,2
                              DO KA = 1,NSOL
                                 DO KB = 1,NSOL
                                    WCG = WE*WT*AMEOPC(IKMCB(KA),
     &                                 IKMCB(KB),MS)
                                    WCF = WE*WT*AMEOPC(IMKMCB(KA),
     &                                 IMKMCB(KB),MS)
                                    WOG = WE*WT*AMEOPO(IKMCB(KA),
     &                                 IKMCB(KB),MS)
                                    WOF = WE*WT*AMEOPO(IMKMCB(KA),
     &                                 IMKMCB(KB),MS)
                                    DO I = 1,IRTOP
                                       ZGZG = ZG(I,KA,LAM1)
     &                                    *ZG(I,KB,LAM2)
                                       ZFZF = ZF(I,KA,LAM1)
     &                                    *ZF(I,KB,LAM2)
                                       RHOOPC(I,IT,MS) = RHOOPC(I,IT,MS)
     &                                    + DIMAG(WCG*ZGZG+WCF*ZFZF)
                                       RHOOPO(I,IT,MS) = RHOOPO(I,IT,MS)
     &                                    + DIMAG(WOG*ZGZG-WOF*ZFZF)
                                    END DO
                                 END DO
                              END DO
                           END DO
C
C    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
                           IF ( K1.EQ.K2 ) THEN
C
                              DO MS = 1,2
                                 DO KA = 1,NSOL
                                    DO KB = 1,NSOL
                                       WCG = WE*AMEOPC(IKMCB(KA),
     &                                    IKMCB(KB),MS)/PI
                                       WCF = WE*AMEOPC(IMKMCB(KA),
     &                                    IMKMCB(KB),MS)/PI
                                       WOG = WE*AMEOPO(IKMCB(KA),
     &                                    IKMCB(KB),MS)/PI
                                       WOF = WE*AMEOPO(IMKMCB(KA),
     &                                    IMKMCB(KB),MS)/PI
                                       DO I = 1,IRTOP
                                         ZGJG = ZG(I,KA,LAM1)
     &                                      *JG(I,KB,LAM2)
                                         ZFJF = ZF(I,KA,LAM1)
     &                                      *JF(I,KB,LAM2)
                                         RHOOPC(I,IT,MS)
     &                                      = RHOOPC(I,IT,MS)
     &                                      + DIMAG(WCG*ZGJG+WCF*ZFJF)
                                         RHOOPO(I,IT,MS)
     &                                      = RHOOPO(I,IT,MS)
     &                                      + DIMAG(WOG*ZGJG-WOF*ZFJF)
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END IF
C
                        END IF
C-----------------------------------------------------------------------
                     END DO
                  END DO
C
C
               END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
C
      IF ( .NOT.EFCORRECT ) RETURN
C
C=======================================================================
C
      DO IT = 1,NT
C
         IF ( OPF(IT) ) THEN
C
            QQO = 0.0D0
            QQC = 0.0D0
            QQS = 0.0D0
            WRITE (6,99001)
            DO MS = 1,2
               DO I = 1,IRTOP
                  RINT(I) = RHOOPC(I,IT,MS)*R2DRDI(I,IM)
               END DO
               CALL RRADINT(IM,RINT,QQCMS)
               DO I = 1,IRTOP
                  RINT(I) = RHOOPO(I,IT,MS)*R2DRDI(I,IM)
               END DO
               CALL RRADINT(IM,RINT,QQOMS)
               WRITE (6,99003) IT,' C ',QQCMS,MS
               WRITE (6,99003) IT,' O ',QQOMS,MS
               QQC = QQC + QQCMS
               QQO = QQO + QQOMS
               QQS = QQS + QQCMS*(-1)**MS
            END DO
            WRITE (6,*)
            WRITE (6,99002) IT,' Q ',OBS_T(0,IDOS,IT),QQC,
     &                      OBS_T(0,IDOS,IT)/QQC
            WRITE (6,99002) IT,' S ',OBS_T(0,ISMT,IT),QQS,
     &                      OBS_T(0,ISMT,IT)/QQS
            WRITE (6,99002) IT,' O ',OBS_T(0,IOMT,IT),QQO,
     &                      OBS_T(0,IOMT,IT)/QQO
            WRITE (6,*)
         END IF
C
      END DO
C
      CALL RINIT(NRMAX*2*NTMAX,AOPT)
C
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      IF ( ORBPOL(1:6).EQ.'BROOKS' ) THEN
C
         IMLAST = 0
         DO IT = 1,NT
C
            IF ( LOPT(IT).EQ.2 ) THEN
               RACPAR = 'B '
               PRE(2) = +9.0D0/441.0D0
               PRE(4) = -5.0D0/441.0D0
               LFMAX = 4
            ELSE IF ( LOPT(IT).EQ.3 ) THEN
               RACPAR = 'E3'
               PRE(2) = +5.0D0/(225.0D0*3.0D0)
               PRE(4) = +6.0D0/(1089.0D0*3.0D0)
               PRE(6) = -91.0D0/(7361.64D0*3.0D0)
               LFMAX = 6
            END IF
C
            IF ( OPF(IT) ) THEN
               IM = IMT(IT)
               IRTOP = JRWS(IM)
C
               IF ( IM.NE.IMLAST ) THEN
                  DO IR = 1,NRMAX
                     RPW(IR,1) = R(IR,IM)
                     DO IL = 2,2*NLMAX
                        RPW(IR,IL) = RPW(IR,IL-1)*R(IR,IM)
                     END DO
                  END DO
                  IMLAST = IM
               END IF
C
               DO MS = 1,2
C --------------------------------------------- normalize densities to 1
                  DO I = 1,IRTOP
                     TL(I) = RHOOPC(I,IT,MS)*R2DRDI(I,IM)
                     TG(I) = RHOOPO(I,IT,MS)*R2DRDI(I,IM)
                  END DO
                  CALL RRADINT(IM,TL,QNORM)
                  CALL RRADINT(IM,TG,QLZ(IT,MS))
                  DO I = 1,IRTOP
                     RHOOPC(I,IT,MS) = RHOOPC(I,IT,MS)/QNORM
                  END DO
C
C -------------------------------------------------- calculate potential
                  DO LF = 2,LFMAX,2
C -------------------------------------- evaluate first radial integrals
                     DO I = 1,IRTOP
                        RHORSQ = 2.0D0*RHOOPC(I,IT,MS)*R2DRDI(I,IM)
                        TL(I) = RHORSQ*RPW(I,LF)
                        TG(I) = RHORSQ/RPW(I,LF+1)
                     END DO
C
                     CALL RRADINT_R(IM,TL,SL)
                     CALL RRADINT_R(IM,TG,SG)
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
                     DO I = 1,IRTOP
                        SL(I) = SL(I)/RPW(I,LF+1) + (SG(IRTOP)-SG(I))
     &                          *RPW(I,LF)
                     END DO
C ----------------------------------------- check Coulomb integrals F^LF
                     DO I = 1,IRTOP
                        SG(I) = RHOOPC(I,IT,MS)*R2DRDI(I,IM)*SL(I)
                     END DO
                     CALL RRADINT(IM,SG,FCLMB)
C ----------------------------------------------------- set up potential
                     DO I = 1,IRTOP
                        AOPT(I,MS,IT) = AOPT(I,MS,IT) - PRE(LF)*SL(I)
     &                                  *QLZ(IT,MS)
                     END DO
C
                     WRITE (6,99004) LF,FCLMB*1000D0,FCLMB*RY_EV,IT,
     &                               TXT_T(IT),MS
                  END DO
C
C ----------------------------------------- check Racah parameter B / E3
                  DO I = 1,IRTOP
                     SG(I) = RHOOPC(I,IT,MS)*R2DRDI(I,IM)*AOPT(I,MS,IT)
                  END DO
                  CALL RRADINT(IM,SG,BRACAH)
                  BRACAH = BRACAH/QLZ(IT,MS)
                  WRITE (6,99005) RACPAR,BRACAH*1000D0,BRACAH*RY_EV,IT,
     &                            TXT_T(IT),MS
                  WRITE (6,99006) QLZ(IT,MS),IT,TXT_T(IT),MS
C
               END DO
               WRITE (6,*) ' '
            END IF
C
C ------ suppress vector potential when dealing with a non-magnetic host
C
            IF ( ABS(QLZ(IT,1)+QLZ(IT,2)).LT.1D-6 ) AOPT(1:IRTOP,1:2,IT)
     &           = 0D0
C
         END DO
C
      END IF
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
99001 FORMAT (' integrals in <SCFOPPOT>')
99002 FORMAT (' IT ',I3,A,F20.10,/,10X,4F20.10)
99003 FORMAT (' IT ',I3,A,F20.10,I5)
99004 FORMAT (' COULOMB-intgral F^',I1,' = ',F8.3,' mRy = ',F8.3,' eV',
     &        ' for  IT=',I2,2X,A,'  MS=',I2)
99005 FORMAT (' RACAH-parameter ',A2,1X,' = ',F8.3,' mRy = ',F8.3,' eV',
     &        ' for  IT=',I2,2X,A,'  MS=',I2)
99006 FORMAT (' <l_z>             ',1X,' = ',F8.3,18X,' for  IT=',I2,2X,
     &        A,'  MS=',I2)
      END
