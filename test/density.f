C*==density.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE DENSITY(IECURR,TAUT,OBS_T)
C   ********************************************************************
C   *                                                                  *
C   * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
C   *                  WITHIN AN ATOMIC CELL                           *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ENERGY,ONLY:WETAB,NETAB
      USE MOD_ANGMOM,ONLY:NOBSMAX,IDOS,ISMT,IOMT,NL,NKMMAX,NKM,AME_G,
     &    IMKM_IKM,NCPLWF
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:IMT,NTMAX,LTXT_T,TXT_T,ITBOT,ITTOP,NCPLWFMAX,
     &    IKMCPLWF
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_FILES,ONLY:IFILCBWF,LDATSET,DATSET,LSYSTEM,SYSTEM,
     &    IFILBUILDBOT,WRBUILDBOT
      USE MOD_CONSTANTS,ONLY:A0_ANG,PI
      IMPLICIT NONE
C*--DENSITY20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DENSITY')
      INTEGER NQNEUMAX,NLEGMAX
      PARAMETER (NQNEUMAX=200,NLEGMAX=3)
C
C Dummy arguments
C
      INTEGER IECURR
      REAL*8 OBS_T(0:3,NOBSMAX,NTMAX)
      COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ARG,JF(:,:,:),JG(:,:,:),WDS,WE,WOF,WOG,WSF,WSG,WT,
     &           ZF(:,:,:),ZFJF,ZFZF,ZG(:,:,:),ZGJG,ZGZG
      REAL*8 AUX,BES0,BES2,CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2),CHKO(:),
     &       CHKQ(:),CHKS(:),CHRMAX,CHRMIN,DQNEU,DX,FAVR(NQNEUMAX),
     &       FORB(NQNEUMAX),FSPN(NQNEUMAX),MJ,MJMAX,MJMIN,ORBMAX,ORBMIN,
     &       QNEU,QNEUMAX,R1M(2,2),RHOCHR(:,:),RHOORB(:,:),RHOR2,
     &       RHOSPN(:,:),RINT(:),RINTO(:),RINTS(:),SCAL,SCALCHR,SCALORB,
     &       SPNMAX,SPNMIN,XTAB(:),YMAX,YMIN,YTAB(:,:)
      COMPLEX*16 CJLZ
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IFIL,IKM1,IKM2,IKMCB(2),IL,IM,IMKM1,IMKM2,IQNEU,
     &        IRTOP,IS,IT,K1,K2,KA,KAP1,KAP2,KB,L,LAM1,LAM2,LFILNAM,
     &        LMAX,MUE,MUETOP,NQNEU,NSOL
      INTEGER IKAPMUE
      CHARACTER*20 LEG(NLEGMAX)
      SAVE CHKO,CHKQ,CHKS,RHOCHR,RHOORB,RHOSPN
C
C*** End of declarations rewritten by SPAG
C
      DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
C
      ALLOCATABLE XTAB,YTAB,RINT,RINTO,RINTS,JF,JG,ZF,ZG,CHKQ,CHKS,CHKO
      ALLOCATABLE RHOCHR,RHOSPN,RHOORB
C
      ALLOCATE (XTAB(NRMAX),YTAB(NRMAX,0:2),RINT(NRMAX))
      ALLOCATE (RINTO(NRMAX),RINTS(NRMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IF ( IECURR.EQ.1 ) THEN
C
         ALLOCATE (CHKQ(NTMAX),CHKS(NTMAX),CHKO(NTMAX))
C
         ALLOCATE (RHOORB(NRMAX,NTMAX))
         ALLOCATE (RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX))
         DO IT = ITBOT,ITTOP
            CHKQ(IT) = 0.0D0
            CHKS(IT) = 0.0D0
            CHKO(IT) = 0.0D0
            RHOCHR(:,IT) = 0.0D0
            RHOSPN(:,IT) = 0.0D0
            RHOORB(:,IT) = 0.0D0
         END DO
      END IF
C
C
      WE = WETAB(IECURR,1)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,1,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
         LMAX = NL - 1
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,LMAX
            IL = L + 1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            MUE = 0
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
               IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
               IMKM1 = IMKM_IKM(IKM1)
               IMKM2 = IMKM_IKM(IKM2)
               IKMCB(1) = IKM1
               IKMCB(2) = IKM2
C-----------------------------------------------------------------------
               IF ( IREL.LE.1 ) THEN
                  IKM1 = IL
                  IKM2 = IL
                  IF ( NKM.NE.NL**2 ) WRITE (6,*) 
     &               'warning in <DENSITY>:  IREL<=1 and  NL**2 <> NKM='
     &               ,NKM
               END IF
C
C
C   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
C
               CGG(1,1) = AME_G(IKM1,IKM1,2,ISMT)
               CGG(1,2) = AME_G(IKM1,IKM2,2,ISMT)
               CGG(2,1) = AME_G(IKM2,IKM1,2,ISMT)
               CGG(2,2) = AME_G(IKM2,IKM2,2,ISMT)
               CALL RINIT(4,CGF)
               CGF(1,1) = AME_G(IMKM1,IMKM1,2,ISMT)
               CGF(2,2) = AME_G(IMKM2,IMKM2,2,ISMT)
C
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C-----------------------------------------------------------------------
C
               DO K1 = 1,NSOL
                  LAM1 = IKMCB(K1)
                  DO K2 = 1,NSOL
                     LAM2 = IKMCB(K2)
C
                     WT = -TAUT(LAM2,LAM1,IT)/PI
C
                     DO KA = 1,NSOL
                        DO KB = 1,NSOL
                           WDS = WE*WT*R1M(KA,KB)
                           WSG = WE*WT*CGG(KA,KB)
                           WSF = WE*WT*CGF(KA,KB)
                           WOG = WE*WT*CFG(KA,KB)
                           WOF = WE*WT*CFF(KA,KB)
                           DO I = 1,IRTOP
                              ZGZG = ZG(I,KA,LAM1)*ZG(I,KB,LAM2)
                              ZFZF = ZF(I,KA,LAM1)*ZF(I,KB,LAM2)
                              RHOCHR(I,IT) = RHOCHR(I,IT)
     &                           + DIMAG(WDS*ZGZG+WDS*ZFZF)
                              RHOSPN(I,IT) = RHOSPN(I,IT)
     &                           + DIMAG(WSG*ZGZG-WSF*ZFZF)
                              RHOORB(I,IT) = RHOORB(I,IT)
     &                           + DIMAG(WOG*ZGZG-WOF*ZFZF)
                           END DO
                        END DO
                     END DO
C
C    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
                     IF ( K1.EQ.K2 ) THEN
C
                        DO KA = 1,NSOL
                           DO KB = 1,NSOL
                              WDS = WE*R1M(KA,KB)/PI
                              WSG = WE*CGG(KA,KB)/PI
                              WSF = WE*CGF(KA,KB)/PI
                              WOG = WE*CFG(KA,KB)/PI
                              WOF = WE*CFF(KA,KB)/PI
                              DO I = 1,IRTOP
                                 ZGJG = ZG(I,KA,LAM1)*JG(I,KB,LAM2)
                                 ZFJF = ZF(I,KA,LAM1)*JF(I,KB,LAM2)
                                 RHOCHR(I,IT) = RHOCHR(I,IT)
     &                              + DIMAG(WDS*ZGJG+WDS*ZFJF)
                                 RHOSPN(I,IT) = RHOSPN(I,IT)
     &                              + DIMAG(WSG*ZGJG-WSF*ZFJF)
                                 RHOORB(I,IT) = RHOORB(I,IT)
     &                              + DIMAG(WOG*ZGJG-WOF*ZFJF)
                              END DO
                           END DO
                        END DO
                     END IF
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
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( IECURR.EQ.NETAB(1) ) THEN
C
         WRITE (6,99003)
C
         DO IT = ITBOT,ITTOP
C
            IFIL = 80
C
            WRITE (6,'(//)')
C
C=======================================================================
C                   PLOT DENSITIES
C=======================================================================
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            CHRMIN = 0D0
            CHRMAX = 0D0
            DO I = 1,IRTOP
               RINT(I) = RHOCHR(I,IT)*R2DRDI(I,IM)
               RHOR2 = RHOCHR(I,IT)*R(I,IM)**2
               CHRMIN = MIN(CHRMIN,RHOR2)
               CHRMAX = MAX(CHRMAX,RHOR2)
            END DO
            CALL RRADINT(IM,RINT,AUX)
            CHKQ(IT) = AUX - CHKQ(IT)
C
            SPNMIN = +1D20
            SPNMAX = -1D20
            DO I = 1,IRTOP
               RINT(I) = RHOSPN(I,IT)*R2DRDI(I,IM)
               RHOR2 = RHOSPN(I,IT)*R(I,IM)**2
               SPNMIN = MIN(SPNMIN,RHOR2)
               SPNMAX = MAX(SPNMAX,RHOR2)
            END DO
            CALL RRADINT(IM,RINT,AUX)
            CHKS(IT) = AUX - CHKS(IT)
C
            ORBMIN = +1D20
            ORBMAX = -1D20
            DO I = 1,IRTOP
               RINT(I) = RHOORB(I,IT)*R2DRDI(I,IM)
               RHOR2 = RHOORB(I,IT)*R(I,IM)**2
               ORBMIN = MIN(ORBMIN,RHOR2)
               ORBMAX = MAX(ORBMAX,RHOR2)
            END DO
            CALL RRADINT(IM,RINT,AUX)
            CHKO(IT) = AUX - CHKO(IT)
C
            WRITE (6,99005) ' ',IT,TXT_T(IT),' ',' ',OBS_T(0,IDOS,IT),
     &                      CHKQ(IT),OBS_T(0,IDOS,IT) - CHKQ(IT),' ',
     &                      OBS_T(0,ISMT,IT),CHKS(IT),OBS_T(0,ISMT,IT)
     &                      - CHKS(IT),' ',OBS_T(0,IOMT,IT),CHKO(IT),
     &                      OBS_T(0,IOMT,IT) - CHKO(IT)
            WRITE (6,*) ' '
C
            SCAL = 1D0
            SCALCHR = 0.1D0
            SCALORB = 40D0
            YMIN = MIN(CHRMIN*SCALCHR,SPNMIN,ORBMIN*SCALORB)
            YMAX = MAX(CHRMAX*SCALCHR,SPNMAX,ORBMAX*SCALORB)
            YMIN = DBLE(INT(YMIN/SCAL)-1)*SCAL
            YMAX = DBLE(INT(YMAX/SCAL)+1)*SCAL
            DX = 0.5D0
C
            CALL XMGRHEAD(DATSET,LDATSET,'dens',4,TXT_T(IT),LTXT_T(IT),
     &                    FILNAM,80,LFILNAM,IFIL,1,0.0D0,0,R(IRTOP,IM),
     &                    0,YMIN,0,YMAX,0,YMIN,0,YMAX,0,'r (a!s0!N)',10,
     &                    'r!S2!N !xr!0!s'//TXT_T(IT)(1:LTXT_T(IT))
     &                    //'!N (a.u.)',14+LTXT_T(IT)+9,' ',0,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'density curves for '//TXT_T(IT)
     &                    (1:LTXT_T(IT)),(19+LTXT_T(IT)),.FALSE.)
C
            WRITE (IFIL,99004) '#',DATSET(1:LDATSET)
            WRITE (IFIL,99005) '#',IT,TXT_T(IT),'#','#',OBS_T(0,IDOS,IT)
     &                         ,CHKQ(IT),OBS_T(0,IDOS,IT) - CHKQ(IT),
     &                         '#',OBS_T(0,ISMT,IT),CHKS(IT),
     &                         OBS_T(0,ISMT,IT) - CHKS(IT),'#',
     &                         OBS_T(0,IOMT,IT),CHKO(IT),
     &                         OBS_T(0,IOMT,IT) - CHKO(IT)
C
            LEG(1) = 'charge (x 0.1)'
            LEG(2) = 'spin '
            LEG(3) = 'orbital (x 40)'
C
            CALL XMGRLEG1(IFIL,0,3,LEG,0.5D0,0.8D0)
C
            CALL XMGRCURVES(IFIL,1,3,0,2,1,1)
C
            DO I = 1,IRTOP
               YTAB(I,0) = RHOCHR(I,IT)*R(I,IM)**2*SCALCHR
               YTAB(I,1) = RHOSPN(I,IT)*R(I,IM)**2
               YTAB(I,2) = RHOORB(I,IT)*R(I,IM)**2*SCALORB
            END DO
C
            DO IS = 0,2
               CALL XMGRTABLE(0,IS,R(1,IM),YTAB(1,IS),1.0D0,IRTOP,IFIL)
            END DO
C
            CALL XMGRTABLE(0,3,R(1,IM),YTAB(1,0),0.0D0,IRTOP,IFIL)
C
            WRITE (6,99001) FILNAM(1:LFILNAM)
            CLOSE (IFIL)
C
C=======================================================================
C              CALCULATE AND PLOT SCATTERING AMPLITUDES
C=======================================================================
C
C
            NQNEU = MIN(NQNEUMAX,100)
            QNEUMAX = 4*PI*A0_ANG
            DQNEU = QNEUMAX/DBLE(NQNEU-1)
            DO IQNEU = 1,NQNEU
               QNEU = (IQNEU-1)*DQNEU
               DO I = 1,IRTOP
                  ARG = QNEU*R(I,IM)
                  BES0 = DREAL(CJLZ(0,ARG))
                  BES2 = DREAL(CJLZ(2,ARG))
                  RINTS(I) = BES0*RHOSPN(I,IT)*R2DRDI(I,IM)
                  RINTO(I) = (BES0+BES2)*RHOORB(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RRADINT(IM,RINTS,FSPN(IQNEU))
               CALL RRADINT(IM,RINTO,FORB(IQNEU))
C
               FAVR(IQNEU) = FSPN(IQNEU) + FORB(IQNEU)
            END DO
            DO IQNEU = NQNEU,1, - 1
               FAVR(IQNEU) = FAVR(IQNEU)/FAVR(1)
               FSPN(IQNEU) = FSPN(IQNEU)/FSPN(1)
               FORB(IQNEU) = FORB(IQNEU)/FORB(1)
            END DO
C
            DX = DQNEU/(4*PI*A0_ANG)
C
            CALL XMGRHEAD(DATSET,LDATSET,'fmag',4,TXT_T(IT),LTXT_T(IT),
     &                    FILNAM,80,LFILNAM,IFIL,1,0.0D0,0,(NQNEU-1)*DX,
     &                    0,-0.2D0,0,1.2D0,0,-0.2D0,0,1.2D0,0,
     &                    'sin !xq / l !0(!m{1}A!M{1}!S o !N)!S-1!N',40,
     &                    'f!s'//TXT_T(IT)(1:LTXT_T(IT))//'!N (a.u.)',
     &                    3+LTXT_T(IT)+9,' ',0,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'magnetic form factor for '//
     &                    TXT_T(IT)(1:LTXT_T(IT)),(25+LTXT_T(IT)),
     &                    .FALSE.)
C
            WRITE (IFIL,99004) '#',DATSET(1:LDATSET)
            WRITE (IFIL,99005) '#',IT,TXT_T(IT),'#','#',OBS_T(0,IDOS,IT)
     &                         ,CHKQ(IT),OBS_T(0,IDOS,IT) - CHKQ(IT),
     &                         '#',OBS_T(0,ISMT,IT),CHKS(IT),
     &                         OBS_T(0,ISMT,IT) - CHKS(IT),'#',
     &                         OBS_T(0,IOMT,IT),CHKO(IT),
     &                         OBS_T(0,IOMT,IT) - CHKO(IT)
C
            LEG(1) = 'f!savr  !N'
            LEG(2) = 'f!sspin !N'
            LEG(3) = 'f!sorb  !N'
C
            CALL XMGRLEG1(IFIL,0,3,LEG,0.6D0,0.8D0)
C
            CALL XMGRCURVES(IFIL,1,3,0,2,1,1)
C
            DO I = 1,NQNEU
               XTAB(I) = (I-1)*DX
            END DO
C
            CALL XMGRTABLE(0,0,XTAB,FAVR,1.0D0,NQNEU,IFIL)
            CALL XMGRTABLE(0,1,XTAB,FSPN,1.0D0,NQNEU,IFIL)
            CALL XMGRTABLE(0,2,XTAB,FORB,1.0D0,NQNEU,IFIL)
            CALL XMGRTABLE(0,3,XTAB,FAVR,0.0D0,NQNEU,IFIL)
C
            WRITE (6,99002) FILNAM(1:LFILNAM)
            CLOSE (IFIL)
C
         END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) THEN
            WRITE (IFILBUILDBOT,99006) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'FAVR',IT,
     &                                 (FAVR(IQNEU),IQNEU=1,MIN(5,NQNEU)
     &                                 )
            WRITE (IFILBUILDBOT,99006) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'FSPN',IT,
     &                                 (FSPN(IQNEU),IQNEU=1,MIN(5,NQNEU)
     &                                 )
            WRITE (IFILBUILDBOT,99006) ROUTINE(1:LEN_TRIM(ROUTINE)),
     &                                 'FORB',IT,
     &                                 (FORB(IQNEU),IQNEU=1,MIN(5,NQNEU)
     &                                 )
         END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         CALL STOP_REGULAR(ROUTINE,'all done')
C
      ELSE
         DEALLOCATE (XTAB,YTAB,RINT,RINTO,RINTS,JF,JG,ZF,ZG)
      END IF
C
99001 FORMAT (10X,'densities written to               ',A)
99002 FORMAT (10X,'magnetic form factor written to    ',A)
99003 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*            ******       *     *    **     ****             *'
     &  ,/,10X,
     &  '*            *            **   **   *  *   *    *            *'
     &  ,/,10X,
     &  '*            *            * * * *  *    *  *                 *'
     &  ,/,10X,
     &  '*            *****   ***  *  *  *  ******  *  ***            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *  *    *            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *  *    *            *'
     &  ,/,10X,
     &  '*            *            *     *  *    *   ****             *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99004 FORMAT (A1,9X,'DATASET: ',A)
99005 FORMAT (A1,9X,'IT=',I2,3X,A,/,A1,9X,10X,
     &        '     <CALCDOS>   INT rho d3r        DELTA',/,A1,9X,
     &        'CHARGE    ',3F14.6,/,A1,9X,'m_spin    ',3F14.6,/,A1,9X,
     &        'm_orb     ',3F14.6)
99006 FORMAT ('# BUILDBOT: ',A,':  scattering amnplitude ',A,
     &        '  for IT =',I5,/,(1PE22.14))
      END
