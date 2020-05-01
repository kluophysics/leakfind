C*==scfchknval.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFCHKNVAL(EMIN,ECTOP)
C   ********************************************************************
C   *                                                                  *
C   *  this routine checks whether there are bound states above the    *
C   *  states declared as core states but below the lower boundary     *
C   *  EMIN   set for the valence band                                 *
C   *  these states are declared to be core states with the  NCOR      *
C   *  and NVAL adjusted accordingly                                   *
C   *  NOTE: the last 'atomic' (n,l)-shell is assumed to be part of    *
C   *        the valence band                                          *
C   *        the sequence of occupied core states  NQNTAB,  LQNTAB     *
C   *        has to be consistent with the routine <CORE>              *
C   *                                                                  *
C   ********************************************************************
      USE MOD_TABLES,ONLY:TAB_CHSYM
      USE MOD_ANGMOM,ONLY:TXT_L,TXT_J
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NT,NCORT,NVALT,NAT,CONC,VT,Z,IMT,
     &    NVALTOT
      USE MOD_RMESH,ONLY:FULLPOT,JRWS,JRCRI,NRCMAX,R_COR,DRDI_COR,
     &    DRDIOVR_COR
      IMPLICIT NONE
C*--SCFCHKNVAL24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFCHKNVAL')
      INTEGER NLSHELLMAX
      PARAMETER (NLSHELLMAX=15)
      REAL*8 EMINDEFAULT,EDIST
      PARAMETER (EMINDEFAULT=-0.15D0,EDIST=0.065D0)
C
C Dummy arguments
C
      REAL*8 ECTOP,EMIN
C
C Local variables
C
      REAL*8 BB(:),DP(:,:,:),DQ(:,:,:),EC,ECOR(2),ECORSHT(:,:,:),ESCBOT,
     &       FC(:,:,:),GC(:,:,:),VV(:),WP(:,:,:),WQ(:,:,:)
      CHARACTER*1 CHANGE
      INTEGER IA_ERR,ILSHELL,IM,IR,IRTOP,ISOL,IT,JSOL,KAPCOR(2),
     &        KCORSHT(:,:,:),LQN,LQNTAB(NLSHELLMAX),NCORT0(NT),
     &        NLSHCORT(NT),NLSHSEMT(NT),NQN,NQNTAB(NLSHELLMAX),NRC,
     &        NSHTAB(:),NSOL,NSOLSHT(:,:),NSTATES
      LOGICAL KSEMI,MODNCOR
      CHARACTER*4 STR4
C
C*** End of declarations rewritten by SPAG
C
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
C
      ALLOCATABLE VV,BB,FC,GC,DP,DQ,WP,WQ,ECORSHT,NSOLSHT
      ALLOCATABLE NSHTAB,KCORSHT
C
      CALL TRACK_INFO(ROUTINE)
C
      IR = NRCMAX
      ALLOCATE (WP(2,2,IR),DP(2,2,IR),FC(2,2,IR),VV(IR))
      ALLOCATE (WQ(2,2,IR),DQ(2,2,IR),GC(2,2,IR),BB(IR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: CORE -> DQ'
      ALLOCATE (KCORSHT(NLSHELLMAX,2,NT))
      ALLOCATE (NSOLSHT(NLSHELLMAX,NT),NSHTAB(NLSHELLMAX))
      ALLOCATE (ECORSHT(NLSHELLMAX,2,NT),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: CORE -> ECORSHT'
C
      IF ( EMIN.GT.999D0 ) EMIN = EMINDEFAULT
      DO IT = ITBOT,ITTOP
         NCORT0(IT) = NCORT(IT)
         NLSHCORT(IT) = 0
         NLSHSEMT(IT) = 0
         DO ILSHELL = 1,NLSHELLMAX
            NSOLSHT(ILSHELL,IT) = 0
         END DO
      END DO
      DO ILSHELL = 1,NLSHELLMAX
         NSHTAB(ILSHELL) = 2*(2*LQNTAB(ILSHELL)+1)
      END DO
C
      IF ( IPRINT.GE.0 ) WRITE (6,99001)
      ESCBOT = 9999D0
C
C=======================================================================
C                    find and tabulate bound states
C=======================================================================
C
      IF ( IPRINT.GT.0 ) WRITE (6,99002)
C
      DO IT = ITBOT,ITTOP
C
         IF ( Z(IT).LE.2 ) CYCLE
C
         IM = IMT(IT)
C
C-----------------------------------------------------------------------
         IF ( FULLPOT ) THEN
            IRTOP = JRCRI(IM)
         ELSE
            IRTOP = JRWS(IM)
         END IF
C-----------------------------------------------------------------------
         NRC = 2*IRTOP
         IF ( NRC.GT.NRCMAX ) CALL STOP_MESSAGE(ROUTINE,'NRC > NRCMAX')
C
         VV(1:IRTOP) = VT(1:IRTOP,IT)
         VV((IRTOP+1):NRC) = 0.0D0
         BB(:) = 0.0D0
C
         NSTATES = 0
C
         DO ILSHELL = 1,NLSHELLMAX
C
            NQN = NQNTAB(ILSHELL)
            LQN = LQNTAB(ILSHELL)
C
            NSTATES = NSTATES + NSHTAB(ILSHELL)
C
            IF ( NSTATES.GE.Z(IT) ) EXIT
C
            CALL SCFCHKCORE(IT,IM,IRTOP,IPRINT,R_COR(1,IM),
     &                      DRDI_COR(1,IM),DRDIOVR_COR(1,IM),VV,BB,FC,
     &                      GC,DP,DQ,WP,WQ,Z(IT),NRC,NSOL,ECOR,KAPCOR,
     &                      NQN,LQN,NRCMAX)
C
            NSOLSHT(ILSHELL,IT) = NSOL
            DO ISOL = 1,NSOL
               ECORSHT(ILSHELL,ISOL,IT) = ECOR(ISOL)
               KCORSHT(ILSHELL,ISOL,IT) = KAPCOR(ISOL)
            END DO
C
            IF ( IPRINT.GT.1 ) THEN
               DO ISOL = NSOL,1, - 1
                  WRITE (6,99003) TAB_CHSYM(Z(IT)),IT,Z(IT),NQN,
     &                            TXT_L(LQN),TXT_J(IABS(KAPCOR(ISOL))),
     &                            ECORSHT(ILSHELL,ISOL,IT)
               END DO
            END IF
C
C--------------------------------- fix number of core / semi-core shells
C
            IF ( NSTATES.LE.NCORT(IT) ) THEN
               NLSHCORT(IT) = NLSHCORT(IT) + 1
            ELSE
               NLSHSEMT(IT) = NLSHSEMT(IT) + 1
               DO ISOL = 1,NSOLSHT(ILSHELL,IT)
                  ESCBOT = MIN(ESCBOT,ECORSHT(ILSHELL,ISOL,IT))
               END DO
            END IF
C
            IF ( NSOL.EQ.0 ) EXIT
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      END DO
C
C=======================================================================
C                  fix lower bound for semi-core states  ESCBOT
C              update number of core / semi-core shells if necessary
C=======================================================================
C
      MODNCOR = .FALSE.
 100  CONTINUE
      KSEMI = .FALSE.
C
      DO IT = ITBOT,ITTOP
C
         IF ( Z(IT).LE.2 ) CYCLE
C
         ILSHELL = NLSHCORT(IT)
C
         DO ISOL = 1,NSOLSHT(ILSHELL,IT)
            EC = ECORSHT(ILSHELL,ISOL,IT)
C
            IF ( EC.GT.ESCBOT ) THEN
               DO JSOL = 1,NSOLSHT(ILSHELL,IT)
                  ESCBOT = MIN(ESCBOT,ECORSHT(ILSHELL,JSOL,IT))
               END DO
               NCORT(IT) = NCORT(IT) - NSHTAB(ILSHELL)
               NLSHCORT(IT) = NLSHCORT(IT) - 1
               NLSHSEMT(IT) = NLSHSEMT(IT) + 1
               KSEMI = .TRUE.
               MODNCOR = .TRUE.
               EXIT
            END IF
         END DO
      END DO
      IF ( KSEMI ) GOTO 100
C
C=======================================================================
C                  fix upper bound for core states  ECTOP
C=======================================================================
C
      ECTOP = -9999D10
      DO IT = ITBOT,ITTOP
         IF ( Z(IT).LE.2 ) CYCLE
         ILSHELL = NLSHCORT(IT)
         DO ISOL = 1,NSOLSHT(ILSHELL,IT)
            ECTOP = MAX(ECTOP,ECORSHT(ILSHELL,ISOL,IT))
         END DO
      END DO
C
C=======================================================================
C            check and update EMIN = bottom of energy contour
C=======================================================================
C
      IF ( ESCBOT.LT.(EMIN+2*EDIST) ) THEN
         IF ( IPRINT.GE.0 ) THEN
            IF ( EMIN.GT.999D0 ) THEN
               WRITE (6,99006)
            ELSE
               WRITE (6,99007) EMIN
            END IF
         END IF
         IR = INT(ABS(ESCBOT)/EDIST) + 4
         EMIN = -IR*EDIST
         IF ( ABS(EMIN-ECTOP).LT.ABS(EMIN-ESCBOT) )
     &        EMIN = 0.5D0*(ECTOP+ESCBOT)
      ELSE
         ESCBOT = EMIN
      END IF
C
      IF ( IPRINT.GE.0 ) THEN
         IF ( ECTOP.LT.-9D10 ) THEN
            WRITE (6,99008) EMIN,ESCBOT
         ELSE
            WRITE (6,99008) EMIN,ESCBOT,ECTOP
         END IF
         IF ( MODNCOR ) WRITE (6,99009)
      END IF
      IF ( ECTOP.GT.ESCBOT ) STOP 'in <SCFCHKNVAL>: ECTOP > ESCBOT'
c      IF ( EMIN.LT.-3.0D0 ) STOP 'in <SCFCHKNVAL>:  EMIN < -3.0 Ry'
C
C=======================================================================
C                         print table
C=======================================================================
C
      IF ( IPRINT.GE.0 ) WRITE (6,99002)
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         IF ( Z(IT).LE.2 ) CYCLE
C
         IF ( IPRINT.GE.0 ) WRITE (6,*)
C
         NSTATES = 0
C
         DO ILSHELL = 1,NLSHELLMAX
C
            NQN = NQNTAB(ILSHELL)
            LQN = LQNTAB(ILSHELL)
C
            NSTATES = NSTATES + NSHTAB(ILSHELL)
C
            IF ( NSTATES.GE.Z(IT) ) EXIT
C
            NSOL = NSOLSHT(ILSHELL,IT)
C
            DO ISOL = NSOL,1, - 1
               KAPCOR(ISOL) = KCORSHT(ILSHELL,ISOL,IT)
               IF ( ISOL.EQ.2 ) THEN
                  IF ( IPRINT.GE.0 ) WRITE (6,99003) TAB_CHSYM(Z(IT)),
     &                 IT,Z(IT),NQN,TXT_L(LQN),TXT_J(IABS(KAPCOR(ISOL)))
     &                 ,ECORSHT(ILSHELL,ISOL,IT)
               ELSE
                  CHANGE = '-'
                  IF ( NSTATES.LE.NCORT(IT) ) THEN
                     STR4 = 'core'
                  ELSE
                     STR4 = 'semi'
                     IF ( NSTATES.LE.NCORT0(IT) ) CHANGE = '*'
                  END IF
C
                  IF ( IPRINT.GE.0 ) WRITE (6,99003) TAB_CHSYM(Z(IT)),
     &                 IT,Z(IT),NQN,TXT_L(LQN),TXT_J(IABS(KAPCOR(ISOL)))
     &                 ,ECORSHT(ILSHELL,ISOL,IT),NSHTAB(ILSHELL),
     &                 NSTATES,STR4,CHANGE
C
               END IF
            END DO
C
            IF ( NSOL.EQ.0 ) EXIT
C
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         IF ( NCORT(IT).NE.NCORT0(IT) ) THEN
            NVALT(IT) = Z(IT) - NCORT(IT)
            IF ( IPRINT.GE.0 ) WRITE (6,99004) NVALT(IT)
         END IF
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( MODNCOR ) THEN
         NVALTOT = 0.0D0
         DO IT = ITBOT,ITTOP
            NVALTOT = NVALTOT + CONC(IT)*NAT(IT)*NVALT(IT)
         END DO
C
         IF ( IPRINT.GE.0 ) WRITE (6,99005) NVALTOT
      ELSE
         IF ( IPRINT.GE.0 ) WRITE (6,'(/)')
      END IF
C
      DEALLOCATE (VV,BB,FC,GC,DP,DQ,WP,WQ,ECORSHT,NSOLSHT)
      DEALLOCATE (NSHTAB,KCORSHT)
C
99001 FORMAT (//,1X,79('*'),/,34X,'<SCFCHKNVAL>',/,1X,79('*'),//,10X,
     &        'fix core and semi-core states and ',
     &        'lower boundary for valence band',/)
99002 FORMAT (10X,'table of detected bound states  ',
     &        '(ignoring exchange splitting)',//,10X,
     &        'elmt  IT   Z   n  l   j',8X,'E (Ry)',6X,
     &        'NSH  NSUM  state  changed')
99003 FORMAT (10X,A2,2X,3I4,2X,A1,1X,A4,F15.6,2I6,4X,A,5X,A)
99004 FORMAT (10X,'*****  NEW number of valence states  NVALT =',I4,
     &        '  *****')
99005 FORMAT (/,1X,79('*'),/,10X,'total number of valence ',
     &        'electrons adjusted to',F10.3,/,1X,79('*'),//)
99006 FORMAT (10X,'bottom of energy contour EMIN fixed by program',/)
99007 FORMAT (10X,'bottom of energy contour EMIN =',F8.3,' too high',/,
     &        10X,'the value will be adjusted by the program',/)
99008 FORMAT (10X,'bottom of energy contour for band states  EMIN   =',
     &        F8.3,/,10X,
     &        'lower limit for semi-core states          ESCBOT =',F8.3,
     &        /,:,10X,
     &        'upper limit for core states               ECTOP  =',F8.3,
     &        /)
99009 FORMAT (1X,79('*'),/,16X,'the number of core states NCORT',
     &        ' has been adjusted',/,1X,79('*'),/)
      END
