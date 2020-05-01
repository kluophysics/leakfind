C*==socpar.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SOCPAR(IPRINT,TSST,MSST,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *   calculation of the SOC - matrix elements                       *
C   *                                                                  *
C   *   see e.g.  Davenport et al. PRB 37 p 9985 (1988)  Eq. (7)       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:ETAB,NEMAX,EFERMI,NETAB
      USE MOD_CALCMODE,ONLY:ORBPOL,SOLVER
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_ANGMOM,ONLY:NL,NMEMAX,NKMMAX,NKM,NCPLWF
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,SYSTEM,LDATSET,LRECREAL8,
     &    IFILCBWF,IFILBUILDBOT,WRBUILDBOT
      USE MOD_TYPES,ONLY:SOCTL,NTMAX,BT,VT,IMT,NT,LTXT_T,TXT_T,CTL,
     &    NCPLWFMAX,IKMCPLWF
      USE MOD_CONSTANTS,ONLY:RY_EV
      IMPLICIT NONE
C*--SOCPAR21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SOCPAR')
      INTEGER NLEGMAX
      PARAMETER (NLEGMAX=50)
C
C Dummy arguments
C
      INTEGER IPRINT
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      REAL*8 BT0(:,:),C,CSQ,DVDR(:,:),EPLOT(:),EXCME(:,:,:),MJ,
     &       SOCME(:,:,:),YMAX,YMAX1T(:),YMAX2T(:),YMIN
      COMPLEX*16 CINT1(:),CINT2(:),ERYD,JF(:,:,:),JG(:,:,:),MASS,NORM,P,
     &           SSST(NKMMAX,NKMMAX,NTMAX),XI,ZF(:,:,:),ZG(:,:,:)
      CHARACTER*80 FILNAM,SYS
      INTEGER IA_ERR,IC,IE,IFIL,IL,IM,IR,IRTOP,IT,KAP1,L,LAM,LFILNAM,
     &        LSYS,NCURVES,NE
      INTEGER IKAPMUE
      CHARACTER*20 LEG(NLEGMAX)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE CINT1,CINT2,DVDR,SOCME,EXCME,EPLOT,YMAX1T,YMAX2T
      ALLOCATABLE BT0
      ALLOCATABLE ZF,ZG,JF,JG
C
      ALLOCATE (YMAX1T(NTMAX),YMAX2T(NTMAX),BT0(NRMAX,NTMAX))
      ALLOCATE (CINT1(NRMAX),CINT2(NRMAX),DVDR(NRMAX,NT),EPLOT(NEMAX))
      ALLOCATE (SOCME(NEMAX,NT,NL),EXCME(NEMAX,NT,NL),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: EPLOT')
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      WRITE (6,99001)
C
      NE = NETAB(1)
      SYS = SYSTEM
      LSYS = LSYSTEM
C
      CALL XMGRSUBSCRIPTS(SYS,LSYS,80)
C
      OPEN (UNIT=IFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=(5+2*(2+2*NRMAX*2))*LRECREAL8)
C
      DO IT = 1,NT
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CALL RDIFFER(IM,VT(1,IT),DVDR(1,IT))
C
         DO IR = 1,IRTOP
            BT(IR,IT) = -BT(IR,IT)
            CINT1(IR) = R2DRDI(IR,IM)*BT(IR,IT)
         END DO
C
         CALL CRADINT(IM,CINT1,NORM)
C
         IF ( ABS(DREAL(NORM)).LT.0 ) THEN
            DO IR = 1,IRTOP
               BT(IR,IT) = -BT(IR,IT)
            END DO
         END IF
         DO IR = 1,IRTOP
            BT0(IR,IT) = BT(IR,IT)
            BT(IR,IT) = 0D0
         END DO
C
         DO IL = 1,NL
            SOCTL(IT,IL) = 1D-3
         END DO
C
         YMAX1T(IT) = 0D0
         YMAX2T(IT) = 0D0
C
      END DO
C
      SOLVER = 'BS-SOC    '
C
C=======================================================================
      DO IE = 1,NE
         ERYD = DREAL(ETAB(IE,1))
C
         CALL SSITE(1,0,IFILCBWF,.FALSE.,.FALSE.,ERYD,P,IPRINT,NKM,TSST,
     &              MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C------------------------------------------------------------------ T --
         DO IT = 1,NT
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            CALL WAVFUN_READ_REL(IFILCBWF,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C
C------------------------------------------------------------------ L --
            DO IL = 2,NL
               L = IL - 1
               C = CTL(IT,L+1)
               CSQ = C*C
C
               MJ = DBLE(L) + 0.5D0
C
               KAP1 = -L - 1
C
               LAM = IKAPMUE(KAP1,NINT(MJ-0.5D0))
C
               CALL CINIT(NRMAX,CINT1)
C
               DO IR = 1,IRTOP
                  CINT1(IR) = CINT1(IR) + ZG(IR,1,LAM)**2*R2DRDI(IR,IM)
               END DO
C
               CALL CRADINT(IM,CINT1,NORM)
               NORM = 1.0D0/SQRT(NORM)
C
               CALL ZSCAL(IRTOP,NORM,ZG(1,1,LAM),1)
C
               DO IR = 1,IRTOP
                  MASS = 0.5D0*(1.0D0+(ERYD-VT(IR,IT))/CSQ)
C
                  CINT1(IR) = ZG(IR,1,LAM)**2*R2DRDI(IR,IM)*DVDR(IR,IT)
     &                        /R(IR,IM)/MASS**2
C
                  CINT2(IR) = ZG(IR,1,LAM)**2*R2DRDI(IR,IM)*BT0(IR,IT)
               END DO
C
               CALL CRADINT(IM,CINT1,XI)
C
               SOCME(IE,IT,IL) = RY_EV*DREAL(XI)/(2.0D0*CSQ)
C
               YMAX1T(IT) = MAX(YMAX1T(IT),SOCME(IE,IT,IL))
C
               CALL CRADINT(IM,CINT2,XI)
C
               EXCME(IE,IT,IL) = 2D0*DREAL(XI)
C
               YMAX2T(IT) = MAX(YMAX2T(IT),EXCME(IE,IT,IL))
C
            END DO
C------------------------------------------------------------------ L --
C
         END DO
C------------------------------------------------------------------ T --
C
      END DO
C=======================================================================
C
      DO IE = 1,NE
         EPLOT(IE) = (DREAL(ETAB(IE,1))-EFERMI)*RY_EV
      END DO
C
      IFIL = 7
C
      NCURVES = NL - 1
      LEG(1) = 'p'
      LEG(2) = 'd'
      DO IC = 3,NCURVES
         LEG(IC) = CHAR(ICHAR('f')+IC-3)
      END DO
C
      YMIN = 0D0
C
      DO IT = 1,NT
C
         CALL XMGRHEAD(DATSET,LDATSET,'soc',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IFIL,1,EPLOT(1),1,EPLOT(NE),1,
     &                 YMIN,0,YMAX1T(IT),1,YMIN,0,YMAX,0,'energy (eV)',
     &                 11,'!xx!F!sl!N(E) (eV)  (B!sxc!N=0)',31,' ',0,
     &                 'SPR-KKR calculations for '//SYS(1:LSYS),25+LSYS,
     &                 'spin-orbit-coupling parameter of '//TXT_T(IT)
     &                 (1:LTXT_T(IT)),(33+LTXT_T(IT)),.FALSE.)
C
         CALL XMGRLEGEND(IFIL,1,NCURVES,0,LEG,LEG)
C
         CALL XMGRCURVES(IFIL,1,(NL-1),0,2,1,1)
C
         DO IL = 2,NL
            CALL XMGRTABLE(0,(IL-2),EPLOT,SOCME(1,IT,IL),1.0D0,NE,IFIL)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '   spin-orbit-coupling parameter written to file '
     &               ,FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IFIL)
C
C-----------------------------------------------------------------------
C
         CALL XMGRHEAD(DATSET,LDATSET,'exc',3,TXT_T(IT),LTXT_T(IT),
     &                 FILNAM,80,LFILNAM,IFIL,1,EPLOT(1),1,EPLOT(NE),1,
     &                 YMIN,0,YMAX2T(IT),1,YMIN,0,YMAX,0,'energy (eV)',
     &                 11,'!xD!FE!sxc!N (eV)',17,' ',0,
     &                 'SPR-KKR calculations for '//SYS(1:LSYS),25+LSYS,
     &                 'exchange splitting parameter of '//TXT_T(IT)
     &                 (1:LTXT_T(IT)),(32+LTXT_T(IT)),.FALSE.)
C
         CALL XMGRLEGEND(IFIL,1,NCURVES,0,LEG,LEG)
C
         CALL XMGRCURVES(IFIL,1,(NL-1),0,2,1,1)
C
         DO IL = 2,NL
            CALL XMGRTABLE(0,(IL-2),EPLOT,EXCME(1,IT,IL),1.0D0,NE,IFIL)
         END DO
C
         WRITE (6,*) ' '
         WRITE (6,*) '   exchange splitting  parameter written to file '
     &               ,FILNAM(1:LFILNAM)
         WRITE (6,*) ' '
         CLOSE (IFIL)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) THEN
            WRITE (IFILBUILDBOT,99002) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 ((SOCME(IE,IT,IL),IE=1,MIN(3,NE))
     &                                 ,IL=2,NL)
            WRITE (IFILBUILDBOT,99003) ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                                 ((EXCME(IE,IT,IL),IE=1,MIN(3,NE))
     &                                 ,IL=2,NL)
         END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO
C
      CALL STOP_REGULAR(ROUTINE,'all done')
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      ****    ****    ****        *****    ****   *****     *'
     &  ,/,10X,
     &  '*     *    *  *    *  *    *       *    *  *    *  *    *    *'
     &  ,/,10X,
     &  '*     *       *    *  *            *    *  *    *  *    *    *'
     &  ,/,10X,
     &  '*      ****   *    *  *       ***  *****   ******  *****     *'
     &  ,/,10X,
     &  '*          *  *    *  *            *       *    *  *  *      *'
     &  ,/,10X,
     &  '*     *    *  *    *  *    *       *       *    *  *   *     *'
     &  ,/,10X,
     &  '*      ****    ****    ****        *       *    *  *    *    *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT ('# BUILDBOT: ',A,':  spin-orbit-coupling parameter',
     &        '  for IT =',I5,/,(1PE22.14))
99003 FORMAT ('# BUILDBOT: ',A,':  exchange splitting parameter ',
     &        '  for IT =',I5,/,(1PE22.14))
C
      END
