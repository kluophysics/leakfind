C*==nonrelt1.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NONRELT1(MEZJ,MEZZ,MSST,SSST,TAUT,ERYD,P,BCORS)
C   ********************************************************************
C   *                                                                  *
C   *    calculate non-relativistic T1-time                            *
C   *    and print results                                             *
C   *                                                                  *
C   *   KB        1.3807D-16        erg/K                              *
C   *   HBAR      1.05459D-27       erg*s                              *
C   *   RY        13.6058D0         erg                                *
C   *   A0        0.529177D-08      cm                                 *
C   *   MB        9.2741D-21        erg/G                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILCBWF,IPRINT
      USE MOD_ANGMOM,ONLY:NL,NLMAX,NKMMAX,NMEMAX,NCPLWF
      USE MOD_RMESH,ONLY:NRMAX,DRDI,R,JRWS
      USE MOD_CALCMODE,ONLY:ORBPOL,IREL
      USE MOD_LATTICE,ONLY:BRAVAIS
      USE MOD_TYPES,ONLY:NTMAX,RHOCHR,CTL,IMT,NT,BT,VT,IKMCPLWF,
     &    NCPLWFMAX
      USE MOD_CONSTANTS,ONLY:A0_CGS,HBAR_CGS,KB_CGS,MB_CGS,PI,RY_ERG
      IMPLICIT NONE
C*--NONRELT125
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NONRELT1')
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      REAL*8 BCORS(NTMAX)
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      INTEGER A,IA_ERR,IL,ILS,IM,IR,IRTOP,IS,IT,IWRIRRWF,IWRREGWF,L,
     &        NSPIN,NZ,Z(:)
      REAL*8 BEXTRA
      REAL*8 BHF(:,:),C,DOSL(:,:),EXC(:),FAC,FSYM,GAMN,HME(:,:),INUC,
     &       MUNUC,NA1G(:),NEG(:),NT1U(:),NT2G(:),QFAC,QNUC,RCPL(:),
     &       RDIP(:),RDIPTST,RFCT(:),RHO4PI(:,:),RORB(:),RORBTST,RQUP(:)
     &       ,RTH,RTOT,VAUX(:,:),WEXC(:,:)
      LOGICAL CALCINT,GETIRRSOL
      CHARACTER*1 CL(5)
      COMPLEX*16 DTR,DZJ0,DZZ0,E,INTR0S,INTY(:),JG(:,:,:),NORM,RPW3,
     &           TSST(NKMMAX,NKMMAX,NTMAX),YXR(:),ZG(:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      DATA CL/'s','p','d','f','g'/
C
      ALLOCATABLE RDIP,RCPL,DOSL,RFCT,RORB,WEXC,exc,VAUX,INTY,RQUP,Z
      ALLOCATABLE RHO4PI,BHF,HME,NEG,JG,ZG,YXR,NA1G
      ALLOCATABLE NT2G,NT1U
C
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (RFCT(NLMAX),RORB(NLMAX))
      ALLOCATE (RDIP(NLMAX),RCPL(NLMAX),DOSL(NLMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DOSL')
      ALLOCATE (NEG(NTMAX),Z(NTMAX),RQUP(NLMAX))
      ALLOCATE (NA1G(NTMAX),NT2G(NTMAX),NT1U(NTMAX))
      ALLOCATE (BHF(NLMAX,NTMAX),HME(NLMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NEG')
      ALLOCATE (VAUX(NRMAX,2),INTY(NRMAX),WEXC(NRMAX,2),EXC(NRMAX))
      ALLOCATE (RHO4PI(NRMAX,2))
      ALLOCATE (YXR(NRMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: NT1U')
C
      NSPIN = 1
C
      E = ERYD
      P = SQRT(E)
C
      IF ( NT.GT.0 ) CALL STOP_MESSAGE(ROUTINE,
     &                                 'not available !!!!!!!!!!!!!')
C
C ======================================================================
C                       solve SS - differential equation
C ======================================================================
      CALCINT = .FALSE.
      GETIRRSOL = .TRUE.
      IWRREGWF = 1
      IWRIRRWF = 1
C
      CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,ERYD,P,
     &              IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IM = IMT(IT)
C
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_SRA(IFILCBWF,IT,1,ZG,JG,IRTOP,NCPLWF,IKMCPLWF)
C
         DO IS = 1,NSPIN
C
            NZ = NINT(-VT(1,IT)*R(1,IM)/2.0D0)
            Z(IT) = NZ
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO IL = 1,NL
C
               ILS = NL*(IS-1) + IL
C
               C = CTL(IT,IL)
               RTH = NZ*2/(0.5D0*C**2)/2.0D0
C
               L = IL - 1
               IF ( L.GT.0 ) CALL STOP_MESSAGE(ROUTINE,'L > 0')
C
               DO IR = 1,IRTOP
                  YXR(IR) = ZG(IR,1,ILS)**2*DRDI(IR,IM)
               END DO
C
               CALL CRADINT(IM,YXR,DZZ0)
C
C irregular part    Z*J
C
               DO IR = 1,IRTOP
                  YXR(IR) = ZG(IR,1,ILS)*JG(IR,1,ILS)*DRDI(IR,IM)
               END DO
C
               CALL CRADINT(IM,YXR,DZJ0)
C
               NORM = 1.0D0/SQRT(DZZ0)
C
               IF ( L.EQ.0 ) THEN
                  IF ( IREL.EQ.0 ) THEN
                     INTR0S = 0.0D0
                     DO IR = 1,IRTOP
                        INTY(IR) = -(ZG(IR,1,ILS)*NORM/R(IR,IM))
     &                             **2/(4*PI)
                     END DO
C
                  ELSE
                     DO IR = 1,IRTOP
                        DTR = (1.D0/(4D0*PI*R(IR,IM)**2))
     &                        *RTH/((1D0+E/C**2)*R(IR,IM)+RTH)**2
                        YXR(IR) = (ZG(IR,1,ILS)*NORM)**2*DTR*DRDI(IR,IM)
                     END DO
C
                     CALL CRADINT(IM,YXR,INTR0S)
C
                  END IF
C
                  HME(IL,IT) = BEXTRA(R(1,IM),INTY,INTR0S)
                  BHF(IL,IT) = (8*PI/3.0D0)*MB_CGS*HME(IL,IT)/A0_CGS**3
C
               ELSE
                  DO IR = 1,IRTOP
                     RPW3 = R(IR,IM)**3*(1D0+IREL*(E-VT(IR,IT))/C**2)
                     YXR(IR) = (ZG(IR,1,ILS)*NORM)**2*DRDI(IR,IM)/RPW3
                  END DO
C
                  CALL CRADINT(IM,YXR,INTR0S)
C
C                 WRITE (6,*)'<r^-3> for component ',IT,': ',INTR0S
C
                  HME(IL,IT) = BEXTRA(R(1,IM),INTY,INTR0S)
                  BHF(IL,IT) = 2*MB_CGS*HME(IL,IT)/A0_CGS**3
C
               END IF
C
               IF ( IL.EQ.1 ) THEN
                  NA1G(IT) = DIMAG(-TAUT(1,1,IT)*DZZ0+DZJ0)/PI
                  DOSL(1,IT) = NA1G(IT)
               END IF
               IF ( IL.EQ.2 ) THEN
                  NT1U(IT) = DIMAG(-TAUT(2,2,IT)*DZZ0+DZJ0)/PI
                  DOSL(2,IT) = 3*NT1U(IT)
               END IF
               IF ( IL.EQ.3 ) THEN
                  NT2G(IT) = DIMAG(-TAUT(5,5,IT)*DZZ0+DZJ0)/PI
                  NEG(IT) = DIMAG(-TAUT(7,7,IT)*DZZ0+DZJ0)/PI
                  DOSL(3,IT) = 3*NT2G(IT) + 2*NEG(IT)
               END IF
C
C cpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcp
               IF ( IL.EQ.3 ) THEN
                  DO IR = 1,IRTOP
                     VAUX(IR,1) = 0.0D0
                     VAUX(IR,2) = 0.0D0
                     RHO4PI(IR,1) = RHOCHR(IR,IT)
                     RHO4PI(IR,2) = (DREAL(ZG(IR,1,ILS)*NORM)/R(IR,IM))
     &                              **2
                     WEXC(IR,1) = 0.0D0
                     WEXC(IR,2) = 0.0D0
                  END DO
                  CALL EXCVWN(RHO4PI,VAUX,EXC,WEXC,IRTOP,NRMAX)
                  DO IR = 1,IRTOP
                     BT(IR,IT) = (VAUX(IR,1)-VAUX(IR,2))/2.0D0
                  END DO
               END IF
C cpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcp
C
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C cpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcp
C      CALL XXXXCORE(IPRINT,NT,NCORT,CTL,VT,BT,Z,FINITE_NUCLEUS,
C     &          R,R2DRDI,DRDI,JRWS,IMT,
C     &          RHOCHR,RHOSPN,ECORTAB,
C     &          GCOR,FCOR,ECOR,SZCOR,
C     &          KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,NCXRAY,
C     &          LCXRAY,0,BCOR,BCORS,
C     &          SDIA,SMDIA,SOFF,SMOFF,
C     &          SDIA,SMDIA,SOFF,SMOFF, NKMMAX,ISMQHFI,
C     &          NTMAX,NRMAX,NMMAX,NCSTMAX,NLMAX)
C      WRITE (6,*) bcor,bcors
C cpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcpcp
C
      DO IT = 1,NT
         DO IL = 1,NL
            RFCT(IL) = 0.0D0
            RDIP(IL) = 0.0D0
            RORB(IL) = 0.0D0
            RCPL(IL) = 0.0D0
            RQUP(IL) = 0.0D0
         END DO
C
         WRITE (6,99001) IT
C
         CALL TABGAMMAN(Z(IT),GAMN)
         FAC = 4*PI*(KB_CGS/HBAR_CGS)*(HBAR_CGS*GAMN/RY_ERG)**2
         CALL TABMUQIN(Z(IT),A,MUNUC,QNUC,INUC,1)
         QFAC = 0.0D0
         IF ( (INUC.GT.0.5D0) .AND. (ABS(MUNUC).GT.1D-8) )
     &        QFAC = 1.81D0*(1D0+4D0/(2D0*INUC-1D0))*(QNUC/MUNUC)**2
C
C--------------------------------------------------------- cubic systems
C
         IF ( BRAVAIS.GE.12 ) THEN
            IF ( ABS(NEG(IT)+NT2G(IT)).LT.1D-8 ) THEN
               FSYM = 0.6D0
            ELSE
               FSYM = 3*NT2G(IT)/(2*NEG(IT)+3*NT2G(IT))
            END IF
C
            RFCT(1) = FAC*(NA1G(IT)*BHF(1,IT))**2
            RCPL(1) = FAC*BCORS(IT)**2*(3D0*NT2G(IT)**2+2D0*NEG(IT)**2)
C
            RORB(2) = 2D0*FAC*(NT1U(IT)*BHF(2,IT))**2
            RDIP(2) = RORB(2)*(3.D0/10.D0)
            RQUP(2) = RDIP(2)*QFAC
C
            RORBTST = FAC*BHF(3,IT)**2*DOSL(3,IT)**2*(2D0/3D0)
     &                *FSYM*(2D0-(5D0/3D0)*FSYM)
            RORB(3) = 2D0*FAC*BHF(3,IT)**2*NT2G(IT)
     &                *(NT2G(IT)+4D0*NEG(IT))
C
            RDIPTST = FAC*BHF(3,IT)**2*DOSL(3,IT)**2*(1D0/49D0)
     &                *((5D0/3D0)*FSYM**2-2D0*FSYM+2D0)
            RDIP(3) = RORB(3)
     &                *(3.D0*(5D0*FSYM**2-6D0*FSYM+6.D0)/(98.D0*(6.D0-
     &                5D0*FSYM)*FSYM))
            RQUP(3) = RDIP(3)*QFAC
            IF ( IPRINT.GT.0 ) WRITE (*,*) 'TEST PARAMETERS  ',RORBTST,
     &                                RDIPTST
         END IF
C
         RTOT = 0.0D0
         DO IL = 1,NL
            RTOT = RTOT + RFCT(IL) + RORB(IL) + RDIP(IL) + RCPL(IL)
         END DO
C
         WRITE (6,99002) FSYM,(CL(IL),IL=1,NL)
         WRITE (6,99003) (DOSL(IL,IT),IL=1,NL)
         WRITE (6,99004) (HME(IL,IT),IL=1,NL)
         WRITE (6,99005) (BHF(IL,IT)/1.0D+6,IL=1,NL)
         WRITE (6,99006) (RFCT(IL),IL=1,1)
         WRITE (6,99007) (RORB(IL),IL=2,NL)
         WRITE (6,99008) (RDIP(IL),IL=2,NL)
         WRITE (6,99009) (RQUP(IL),IL=2,NL)
         WRITE (6,99010) (RCPL(IL),IL=1,1)
         WRITE (6,99011) RTOT
         WRITE (6,99012) BHF(1,IT)*2*MB_CGS*DOSL(1,IT)/RY_ERG*100
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DEALLOCATE (RDIP,RCPL,DOSL,RFCT,RORB,WEXC,VAUX,INTY,RQUP,Z)
      DEALLOCATE (RHO4PI,BHF,HME,NEG,JG,ZG,YXR,NA1G)
      DEALLOCATE (NT2G,NT1U)
C
99001 FORMAT (/,' nuclear spin lattice relaxation rate',
     &        '  R = 1/T_1*T  for component',I3,/,1X,67('='),/)
99002 FORMAT (' symmetry parameter f (eg/t2g):',F12.4,/,29X,
     &        5(:9X,'(',A1,')'))
99003 FORMAT (' l-resolved DOS    [1/Ry*spin]:',7F12.4)
99004 FORMAT (' hyperfine mat.element in [au]:',7F12.4)
99005 FORMAT (' hyperfine field   B   in [MG]:',7F12.4)
99006 FORMAT (' R  Fermi contact   in [1/s*K]:',7F12.4)
99007 FORMAT (' R  orbital         in [1/s*K]:',12X,7F12.4)
99008 FORMAT (' R  dipolar         in [1/s*K]:',12X,7F12.4)
99009 FORMAT (' R  quadrupolar     in [1/s*K]:',12X,7F12.4)
99010 FORMAT (' R  core polaris.   in [1/s*K]:',12X,7F12.4)
99011 FORMAT (' ------------------------------------------',/,
     &        ' R  total           in [1/s*K]:',F12.4,/)
99012 FORMAT (' K  Fermi contact       in [%]:',F12.4,/)
      END
