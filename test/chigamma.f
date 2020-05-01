C*==chigamma.f    processed by SPAG 6.70Rc at 16:04 on 16 Jan 2017
      SUBROUTINE CHIGAMMA(IFILCBWF,GAMMAM,GAMMAN,BXCNMM,BXCNNM,BXCNMN,
     &                    BXCNNN,AXCN,NEFTL,ORBSQINT,EMM,ENM,ENN)
C   ********************************************************************
C   *                                                                  *
C   *  Calculate B_xc via GAMMA(r)=|phi(r,E_F)|^2 (normalized)         *
C   *  GAMMA is normalized: int gamma(r)*r^2 dr = 1                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT,IOTMP,DATSET,LDATSET
      USE MOD_TYPES,ONLY:NTMAX,NT,Z,TXT_T,LTXT_T,IMT,RHOCHR,RHOSPN,
     &    NCPLWFMAX,IKMCPLWF
      USE MOD_SITES,ONLY:IQAT
      USE MOD_ANGMOM,ONLY:AMEOPO,AMEOPC,NLMAX,NLQ,NCPLWF,NKMMAX
      USE MOD_CONSTANTS,ONLY:PI,RY_EV
      USE MOD_CALCMODE,ONLY:ORBPOL
      USE MOD_RMESH,ONLY:NRMAX,JRWS,R,R2DRDI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHIGAMMA')
C
C Dummy arguments
C
      INTEGER IFILCBWF
      REAL*8 AXCN(NRMAX,2,NLMAX,NTMAX),BXCNMM(NRMAX,NTMAX),
     &       BXCNMN(NRMAX,NTMAX),BXCNNM(NRMAX,NTMAX),BXCNNN(NRMAX,NTMAX)
     &       ,EMM(NRMAX,NTMAX),ENM(NRMAX,NTMAX),ENN(NRMAX,NTMAX),
     &       GAMMAM(NRMAX,NTMAX),GAMMAN(NRMAX,NTMAX),NEFTL(NTMAX,NLMAX),
     &       ORBSQINT(0:NLMAX,2)
C
C Local variables
C
      COMPLEX*16 AXCINT(2),CINTAXC(NRMAX,2),CSUM,GAMMAINT(0:NLMAX,2),
     &           GAMMAINTG(NRMAX,0:NLMAX,2),GAMMAL(NRMAX,NLMAX,2),
     &           GAMMAS(NRMAX,NTMAX,2),JF(:,:,:),JG(:,:,:),NORM,
     &           ORBINTZGZG(NRMAX),ORBSQINTG(NRMAX,0:NLMAX,2),
     &           ORBSQL(NRMAX,NLMAX,2),ORBZGZG(NRMAX),PHIINTZGZG(NRMAX),
     &           PHIZGZG(NRMAX),ZF(:,:,:),ZG(:,:,:)
      REAL*8 COG(2,2,2),CSG(2,2,2),GAMMA(NRMAX,NTMAX),GAMMAINT0,
     &       GAMMAINTG0(NRMAX),GAMMAOP(NRMAX,2),MJ,NEFT(NTMAX),
     &       RHORINT(NTMAX),RINTRHO(NRMAX),RINTSTO(NRMAX),RPWOP(:,:),
     &       STONERI(NTMAX)
      CHARACTER*80 DATFMT,FILNAM
      INTEGER I5,IA_ERR,IKM1,IKM2,IKMCB(2),IL,IM,IQ,IR,IRTOP,IS,IT,K,K1,
     &        K2,KAP1,KAP2,L,LAM,LDATFMT,LFN,LOP,MJM05,NL,NLOP,NSOL
      INTEGER IKAPMUE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZF,ZG,JF,JG,RPWOP
C
      ALLOCATE (RPWOP(NRMAX,NLMAX*2))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      WRITE (6,99001)
C
      LOP = 2
      IF ( ORBPOL(1:8).EQ.'BROOKS-F' ) LOP = 3
      NLOP = LOP + 1
C
      CALL RINIT(NRMAX*2*NLMAX*NTMAX,AXCN)
C
      I5 = NRMAX*NTMAX
      CALL RINIT(I5,GAMMA)
      CALL CINIT(I5*2,GAMMAS)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         NL = NLQ(IQ)
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         NEFT(IT) = 0.0D0
C
         DO IL = 1,NL
            NEFT(IT) = NEFT(IT) + NEFTL(IT,IL)
         END DO
         IF ( NEFT(IT).LE.1D-8 ) THEN
            WRITE (6,99005) IT
            NEFTL(IT,1:NL) = 1D0
            NEFT(IT) = NL
         END IF
C
         DO IR = 1,IRTOP
            RPWOP(IR,1) = R(IR,IM)
            DO IL = 2,2*NLMAX
               RPWOP(IR,IL) = RPWOP(IR,IL-1)*R(IR,IM)
            END DO
         END DO
C
         I5 = NRMAX*NLMAX*2
         CALL CINIT(I5,GAMMAL)
         CALL CINIT(I5,ORBSQL)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,(NLQ(IQ)-1)
C
            IL = L + 1
C
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MJM05 = -L - 1, + L
               MJ = DBLE(MJM05) + 0.5D0
C
               IF ( ABS(MJ).GE.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C-----------------------------------------------------------------------
               IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
               IKMCB(1) = IKM1
               IKMCB(2) = IKM2
C-----------------------------------------------------------------------
               DO IS = 1,2
                  CSG(1,1,IS) = AMEOPC(IKM1,IKM1,IS)
                  CSG(2,2,IS) = AMEOPC(IKM2,IKM2,IS)
                  CSG(1,2,IS) = AMEOPC(IKM1,IKM2,IS)
                  CSG(2,1,IS) = AMEOPC(IKM2,IKM1,IS)
C
                  COG(1,1,IS) = AMEOPO(IKM1,IKM1,IS)
                  COG(2,2,IS) = AMEOPO(IKM2,IKM2,IS)
                  COG(1,2,IS) = AMEOPO(IKM1,IKM2,IS)
                  COG(2,1,IS) = AMEOPO(IKM2,IKM1,IS)
               END DO
C-----------------------------------------------------------------------
C
C
C========== NORMALIZE ZG WITH |MUE|=J TO 1 AND THEN CALCULATE WGT*|ZG|^2
               DO IS = 1,2
                  DO K = 1,NSOL
                     LAM = IKMCB(K)
C
                     CALL CINIT(NRMAX,PHIZGZG)
                     CALL CINIT(NRMAX,ORBZGZG)
C
                     DO IR = 1,IRTOP
                        DO K1 = 1,NSOL
                           DO K2 = 1,NSOL
                              PHIZGZG(IR) = PHIZGZG(IR) + ZG(IR,K1,LAM)
     &                           *ZG(IR,K2,LAM)*CSG(K1,K2,IS)
                              ORBZGZG(IR) = ORBZGZG(IR) + ZG(IR,K1,LAM)
     &                           *ZG(IR,K2,LAM)*COG(K1,K2,IS)
                           END DO
                        END DO
                        PHIINTZGZG(IR) = PHIZGZG(IR)*R2DRDI(IR,IM)
                        ORBINTZGZG(IR) = ORBZGZG(IR)*R2DRDI(IR,IM)
                     END DO
C
                     CALL CRADINT(IM,PHIINTZGZG,NORM)
C
                     IF ( CDABS(NORM).GT.1.D-15 ) THEN
                        DO IR = 1,IRTOP
                           GAMMAL(IR,IL,IS) = GAMMAL(IR,IL,IS)
     &                        + PHIZGZG(IR)/NORM/DBLE(2*L+1)/DBLE(NSOL)
                        END DO
                     END IF
C
                     CALL CRADINT(IM,ORBINTZGZG,NORM)
C
                     IF ( CDABS(NORM).GT.1.D-15 ) THEN
                        DO IR = 1,IRTOP
                           ORBSQL(IR,IL,IS) = ORBSQL(IR,IL,IS)
     &                        + ORBZGZG(IR)/NORM/DBLE(2*L+1)/DBLE(NSOL)
                        END DO
                     END IF
C
                  END DO
               END DO
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
C ==================================== AVERAGE |ZG|^2 ACCORDING TO JANAK
            DO IS = 1,2
               DO IR = 1,IRTOP
                  GAMMAS(IR,IT,IS) = GAMMAS(IR,IT,IS) + GAMMAL(IR,IL,IS)
     &                               *NEFTL(IT,IL)/NEFT(IT)
                  GAMMAINTG(IR,IL,IS) = GAMMAL(IR,IL,IS)*R2DRDI(IR,IM)
                  ORBSQINTG(IR,IL,IS) = ORBSQL(IR,IL,IS)*R2DRDI(IR,IM)
               END DO
            END DO
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         DO IR = 1,IRTOP
            GAMMAM(IR,IT) = 0.5D0*DREAL(GAMMAS(IR,IT,1)+GAMMAS(IR,IT,2))
            GAMMAN(IR,IT) = 0.D0
            GAMMAINTG0(IR) = GAMMAM(IR,IT)*R2DRDI(IR,IM)
            DO IS = 1,2
               GAMMAINTG(IR,0,IS) = GAMMAS(IR,IT,IS)*R2DRDI(IR,IM)
               GAMMAOP(IR,IS) = DREAL(GAMMAL(IR,NLOP,IS))
            END DO
         END DO
C
         CALL RRADINT(IM,GAMMAINTG0,GAMMAINT0)
C
         DO IS = 1,2
            CALL CRADINT(IM,GAMMAINTG(1,0,IS),GAMMAINT(0,IS))
         END DO
C
         DO IL = 1,NL
            DO IS = 1,2
               CALL CRADINT(IM,GAMMAINTG(1,IL,IS),GAMMAINT(IL,IS))
            END DO
         END DO
C
C ===> GAMMA(L)INT=1 ==> GAMMA(L) IS NORMALIZED TO 1 WITHIN WS
C                        (CONTAINS AVERAGING OVER SPHERE WITH RADIUS R)
         DO IL = 1,NL
            DO IS = 1,2
               CALL CRADINT(IM,ORBSQINTG(1,IL,IS),CSUM)
               ORBSQINT(IL,IS) = DREAL(CSUM)
            END DO
         END DO
C
C=============================== CALCULATE BXC AND TEST STONER-PARAMETER
C
         CALL CHIVWN(RHOCHR(1,IT),RHOSPN(1,IT),EMM(1,IT),ENM(1,IT),
     &               ENN(1,IT),IT,IRTOP)
C
         DO IR = 1,IRTOP
            BXCNMM(IR,IT) = EMM(IR,IT)*GAMMAM(IR,IT)/4D0/PI
            BXCNMN(IR,IT) = ENM(IR,IT)*GAMMAM(IR,IT)/4D0/PI
            BXCNNM(IR,IT) = ENM(IR,IT)*GAMMAN(IR,IT)/4D0/PI
            BXCNNN(IR,IT) = ENN(IR,IT)*GAMMAN(IR,IT)/4D0/PI
C
            RINTSTO(IR) = BXCNMM(IR,IT)*GAMMAM(IR,IT)*R2DRDI(IR,IM)
            RINTRHO(IR) = RHOCHR(IR,IT)*R2DRDI(IR,IM)
         END DO
C
         CALL RRADINT(IM,RINTSTO,STONERI(IT))
         CALL RRADINT(IM,RINTRHO,RHORINT(IT))
C
C=============================== CALCULATE AXC =========================
C
         CALL CHIORBPOL(AXCN(1,1,1,IT),GAMMAOP,ORBSQINT,R2DRDI(1,IM),
     &                  RPWOP,Z(IT),ORBPOL,TXT_T(IT),IT,IRTOP,NLMAX,
     &                  NRMAX)
C
         DO IS = 1,2
            DO IR = 1,IRTOP
               CINTAXC(IR,IS) = AXCN(IR,IS,NLOP,IT)*GAMMAOP(IR,IS)
     &                          *R2DRDI(IR,IM)
            END DO
            CALL CRADINT(IM,CINTAXC(1,IS),AXCINT(IS))
         END DO
C
C------------------------------------------ test output for arbitrary NL
C
         IF ( IPRINT.GT.0 ) THEN
C
            WRITE (6,'(A,I1,A,1E15.7,A)')
     &             ('IS=',IS,' : AXCINT= ',DREAL(AXCINT(IS))*RY_EV,
     &             '[EV/SPIN]',IS=1,2)
C
            WRITE (DATFMT,'(''('',I2,''(1E12.4,1X))'')') 5 + 2*NL
            LDATFMT = 15
C
            CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_g-spin.dat',
     &                    11,FILNAM,LFN,2,IOTMP,'g-spn file',10,NTMAX)
C
            WRITE (IOTMP,99002)
            WRITE (IOTMP,FMT=DATFMT(1:LDATFMT))
     &             (R(IR,IMT(IT)),EMM(IR,IT),GAMMAM(IR,IT),
     &             (DREAL(GAMMAS(IR,IT,IS)),
     &             (DREAL(GAMMAL(IR,IL,IS)),IL=1,NL),IS=1,2),IR=1,IRTOP)
            CLOSE (IOTMP)
C
            CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,'_g-orb.dat',
     &                    10,FILNAM,LFN,2,IOTMP,'g-orb file',10,NTMAX)
C
            WRITE (IOTMP,99003)
            WRITE (IOTMP,'(4(1E15.7,1X))')
     &             (R(IR,IM),GAMMAOP(IR,1),GAMMAOP(IR,2),RHOCHR(IR,IT)
     &             /RHORINT(IT),IR=1,IRTOP)
            CLOSE (IOTMP)
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      WRITE (6,99004) 
     &         '  I [eV]  |  n(E_F) [1/(eV*spin)]  |  I*n(E_F)   |    Q'
     &         ,(IT,2D0*STONERI(IT)*RY_EV,NEFT(IT)/RY_EV/2D0,STONERI(IT)
     &         *NEFT(IT),RHORINT(IT),IT=1,NT)
C
C=======================================================================
99001 FORMAT (//,1X,79('*'),/,35X,'<CHIGAMMA>',/,1X,79('*'),/)
99002 FORMAT ('# FILE fort.(500+IT)',/,
     &        '# INIT: R, STONER, GAMMAM, GAMMAS(1),',
     &        ' GAMMAL(L=0..LMAX,1), GAMMAS(2), GAMMAL(L=0..LMAX,2)')
99003 FORMAT ('# FILE fort.(600+IT)',/,
     &        '# INIT: R, GAMMAOP(1), GAMMAOP(2), RHOCHR')
99004 FORMAT (/,39('ef'),/,/,12X,A,/,1X,79('-'),/,
     &        80(:,' IT= ',I2,2X,F13.7,5X,F13.7,5X,F13.7,5X,F13.7,/))
99005 FORMAT (2(/,1X,79('W')),/,10X,'WARNING from <CHIGAMMA>',/,10X,
     &        'DOS(EF) not supplied in input for IT =',I4,/,10X,
     &        'NEFTL set to 1 !  ---',
     &        '  the Stoner enhancement will be unrealistic',
     &        2(/,1X,79('W')),//)
      END
