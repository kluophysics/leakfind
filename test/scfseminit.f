C*==scfseminit.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFSEMINIT(SEMICORE,NT,NL,NCORT0,NCORT,EMIN,EMINSEMCOR,
     &                      EMAXSEMCOR,NESEMCOR,ECORTAB,NSEMCORSHLT,
     &                      SEMCORSHLT,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  find out which core shells should be treated as semicore        *
C   *                                                                  *
C   *  - pick core shell highest in energy                             *
C   *  - include all core shells belonging to the same (n,l)-shell     *
C   *  - include all core shells that are not separated by more than   *
C   *    the gap ESEGAP                                                *
C   *  - the semicore sates should not include states with energy      *
C   *    below  ESEMCORTHRESHOLD  because of possible problems with    *
C   *    the KKR structure constants                                   *
C   *  - if no promlems were encountered enlarge the energy window     *
C   *    by  ESEMCOROFF  up- and downwards                             *
C   *                                                                  *
C   *  NOTE:       NQNTAB and LQNTAB have to be as used in <CORE>      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:TXT_L
      IMPLICIT NONE
C*--SCFSEMINIT24
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFSEMINIT')
      REAL*8 ESEMCORTHRESHOLD,ESEMCOROFF,ESEMCORGAP
      PARAMETER (ESEMCORTHRESHOLD=-5.5D0,ESEMCOROFF=0.15D0,
     &           ESEMCORGAP=0.5D0)
      INTEGER NCMAX,NLSHELLMAX
      PARAMETER (NCMAX=120,NLSHELLMAX=15)
C
C Dummy arguments
C
      REAL*8 EMAXSEMCOR,EMIN,EMINSEMCOR
      INTEGER NESEMCOR,NL,NT,NTMAX
      LOGICAL SEMICORE
      REAL*8 ECORTAB(120,NTMAX)
      INTEGER NCORT(NTMAX),NCORT0(NTMAX),NSEMCORSHLT(NTMAX)
      CHARACTER*2 SEMCORSHLT(NLSHELLMAX,NTMAX)
C
C Local variables
C
      LOGICAL DONE
      REAL*8 EMINSEMCORL,MJ
      INTEGER I,IA_ERR,IC,ILSHELL,IT,JC,KSEMCORT(:,:),L,LIC(NCMAX),
     &        LQNSEMCOR(NLSHELLMAX),LQNTAB(NLSHELLMAX),MUEM05,N,NCT(:),
     &        NIC(NCMAX),NQN,NQNSEMCOR(NLSHELLMAX),NQNTAB(NLSHELLMAX),
     &        NSOL,S
C
C*** End of declarations rewritten by SPAG
C
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
C
      ALLOCATABLE KSEMCORT,NCT
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (KSEMCORT(NCMAX,NTMAX),NCT(NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'alloc -> NCT')
C
      IC = 0
      DO ILSHELL = 1,NLSHELLMAX
         NQN = NQNTAB(ILSHELL)
         L = LQNTAB(ILSHELL)
         DO MUEM05 = -L - 1, + L
            MJ = MUEM05 + 0.5D0
            IF ( ABS(MJ).GT.L ) THEN
               NSOL = 1
            ELSE
               NSOL = 2
            END IF
            DO S = 1,NSOL
               IC = IC + 1
               IF ( IC.LT.NCMAX ) THEN
                  NIC(IC) = NQN
                  LIC(IC) = L
               END IF
            END DO
         END DO
      END DO
C
      DO IT = 1,NTMAX
         NCT(IT) = 0
         DO IC = 1,NCMAX
            IF ( ECORTAB(IC,IT).LT.-1D-6 ) NCT(IT) = NCT(IT) + 1
            KSEMCORT(IC,IT) = 0
         END DO
      END DO
      NQNSEMCOR(1) = 0
      LQNSEMCOR(1) = 0
C
C   --------------------------------------------------------------------
C                   find upper threshold for semi core states EMAXSEMCOR
C
      EMAXSEMCOR = -1D+20
      DO IT = 1,NT
         DO IC = 1,NCT(IT)
            EMAXSEMCOR = MAX(EMAXSEMCOR,ECORTAB(IC,IT))
         END DO
      END DO
C
C   --------------------------------------------------------------------
C                   find lower threshold for semi core states EMINSEMCOR
C
      EMINSEMCOR = EMAXSEMCOR
      EMINSEMCORL = EMAXSEMCOR
 100  CONTINUE
      DO IT = 1,NT
         DO IC = NCT(IT),1, - 1
C
            IF ( ECORTAB(IC,IT).GT.(EMINSEMCOR-ESEMCORGAP) ) THEN
               EMINSEMCOR = MIN(EMINSEMCOR,ECORTAB(IC,IT))
               KSEMCORT(IC,IT) = 1
C
               DO JC = 1,NCT(IT)
                  IF ( NIC(JC).EQ.NIC(IC) .AND. LIC(JC).EQ.LIC(IC) )
     &                 THEN
                     EMINSEMCOR = MIN(EMINSEMCOR,ECORTAB(JC,IT))
                     KSEMCORT(JC,IT) = 1
                  END IF
               END DO
C
            END IF
C
         END DO
      END DO
C
      IF ( ABS(EMINSEMCORL-EMINSEMCOR).GT.1D-6 ) THEN
         EMINSEMCORL = EMINSEMCOR
         GOTO 100
      END IF
C
      EMINSEMCOR = EMINSEMCOR - ESEMCOROFF
      NESEMCOR = INT((EMAXSEMCOR-EMINSEMCOR)*10D0)
C
      IF ( EMAXSEMCOR.LT.ESEMCORTHRESHOLD ) THEN
         WRITE (6,99001) 'upper',EMAXSEMCOR,ESEMCORTHRESHOLD
         SEMICORE = .FALSE.
         RETURN
      END IF
      IF ( EMINSEMCOR.LT.ESEMCORTHRESHOLD ) THEN
         WRITE (6,99001) 'lower',EMAXSEMCOR,ESEMCORTHRESHOLD
         CALL STOP_MESSAGE(ROUTINE,' --- ')
      END IF
C
C   -------------------------------------------------------------- CHECK
C
      DO IT = 1,NT
         NCORT(IT) = 0
         DO IC = 1,NCT(IT)
            IF ( KSEMCORT(IC,IT).EQ.0 ) THEN
               NCORT(IT) = NCORT(IT) + 1
               IF ( ECORTAB(IC,IT).GT.EMINSEMCOR ) THEN
                  WRITE (6,99004) IT,NIC(IC),LIC(IC),ECORTAB(IC,IT),
     &                            'EMINSEMCOR',EMINSEMCOR
                  CALL STOP_MESSAGE(ROUTINE,' --- ')
               END IF
            END IF
C
            IF ( ECORTAB(IC,IT).GT.EMAXSEMCOR ) THEN
               WRITE (6,99004) IT,NIC(IC),LIC(IC),ECORTAB(IC,IT),
     &                         'EMAXSEMCOR',EMAXSEMCOR
               CALL STOP_MESSAGE(ROUTINE,' --- ')
            END IF
C
            IF ( IC.GT.1 ) THEN
               IF ( KSEMCORT(IC,IT).EQ.0 .AND. KSEMCORT(IC-1,IT).EQ.1 )
     &              THEN
                  WRITE (6,99002) IT,NIC(IC),LIC(IC)
                  CALL STOP_MESSAGE(ROUTINE,' --- ')
               END IF
            END IF
C
            IF ( LIC(IC).GT.(NL-1) ) THEN
               WRITE (6,99003) IT,NIC(IC),LIC(IC),(NL-1)
               CALL STOP_MESSAGE(ROUTINE,' --- ')
            END IF
C
         END DO
      END DO
C
C   --------------------------------------------------------------------
C
      WRITE (6,99005) EMINSEMCOR,EMAXSEMCOR,NESEMCOR
      IF ( EMAXSEMCOR.GT.EMIN ) WRITE (6,99006) EMIN,EMAXSEMCOR
      WRITE (6,99007)
C
      DO IT = 1,NT
         N = 0
         DO IC = 1,NCT(IT)
            IF ( KSEMCORT(IC,IT).EQ.1 ) THEN
               DONE = .FALSE.
               DO I = 1,N
                  IF ( NIC(IC).EQ.NQNSEMCOR(I) .AND. LIC(IC)
     &                 .EQ.LQNSEMCOR(I) ) DONE = .TRUE.
               END DO
               IF ( .NOT.DONE ) THEN
                  N = N + 1
                  NQNSEMCOR(N) = NIC(IC)
                  LQNSEMCOR(N) = LIC(IC)
                  WRITE (SEMCORSHLT(N,IT),FMT='(I1,A1)') NIC(IC),
     &                   TXT_L(LIC(IC))
               END IF
            END IF
         END DO
         NSEMCORSHLT(IT) = N
         WRITE (6,99008) IT,NCORT0(IT),NCORT(IT),NCORT0(IT) - NCORT(IT),
     &                   (SEMCORSHLT(I,IT),I=1,NSEMCORSHLT(IT))
      END DO
C
      WRITE (6,*)
C
      DEALLOCATE (KSEMCORT,NCT,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'dealloc NCT')
C
      RETURN
C
99001 FORMAT (//,1X,79('*'),/,10X,A,' threshold for semi core states',
     &        F12.5,/,10X,'below acceptable limit ESEMCORTHRESHOLD =',
     &        F12.5,/)
99002 FORMAT (//,1X,79('*'),/,10X,'for atom type ',i3,/,10X,'n=',i1,
     &        ' l=',i1,'not recognized as semi core state')
99003 FORMAT (//,1X,79('*'),/,10X,'for atom type ',i3,
     &        ' and core state with  n=',i1,' l=',i1,':',/,10X,
     &        'angular momentum > l_max (valence band)=',I2)
99004 FORMAT (//,1X,79('*'),/,10X,'for atom type ',i3,
     &        ' and core state with  n=',i1,' l=',i1,':',/,10X,'E = ',
     &        f10.5,'  >  ',A,' = ',f10.5)
99005 FORMAT (//,1X,79('*'),/,34X,'<SCFSEMINIT>',/,1X,79('*'),//,10X,
     &        'EMINSEMCOR =',F10.5,'    EMAXSEMCOR =',F10.5,
     &        '     NESEMCOR = ',I3,/)
99006 FORMAT (//,1X,79('#'),/,34X,'WARNING',/,1X,79('#'),//,10X,
     &        'EMIN       =',F10.5,' <  EMAXSEMCOR =',F10.5,/)
99007 FORMAT (/,10X,'IT  NCOR0  NCOR  NSEMCOR   semi core shells ')
99008 FORMAT (10X,I2,3I6,3X,5(4X,A2))
      END
