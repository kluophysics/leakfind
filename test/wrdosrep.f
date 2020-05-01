C*==wrdosrep.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WRDOSREP(DOSREP,DOSCORSYS,MEZZ,MEZJ,TAUT,IECURR)
C   ********************************************************************
C   *                                                                  *
C   *    calculate the l,m_l-projected DOS                             *
C   *                                                                  *
C   *   DOSREP =  RLM     REAL spher.  harmonics - representation      *
C   *             CLM     COMPLEX sph. harmonics - representation      *
C   *             REL     RELATIVISTIC             representation      *
C   *                                                                  *
C   * ---------------------------------------------------------------- *
C   *                                                                  *
C   * CHECK:   a cubic, paramagnetic system is assumed                 *
C   *          and the symmetry imposed degeneracies are checked       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ENERGY,ONLY:EFERMI,ETAB,NETAB
      USE MOD_SITES,ONLY:IQAT,DROTQ
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKMQ,NKM
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,NLT,LTXT_T,ITBOT,ITTOP
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0,C1,PI,RY_EV
      IMPLICIT NONE
C*--WRDOSREP26
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WRDOSREP')
      REAL*8 TOL
      PARAMETER (TOL=1D-6)
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      CHARACTER*3 DOSCORSYS,DOSREP
      INTEGER IECURR
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      LOGICAL DIRECT
      REAL*8 DOS(:,:,:)
      CHARACTER*80 FILNAM
      INTEGER IA_ERR,ICALL,IE,IQ,IT,J,LFN,LL,M,N,NH
      COMPLEX*16 MEZJ_LOC(:,:),MEZZ_LOC(:,:),TAUT_LOC(:,:),
     &           W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
      CHARACTER*12 STR12
      CHARACTER*7 STR7
      SAVE DOS,MEZJ_LOC,MEZZ_LOC,TAUT_LOC
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE DOS,MEZZ_LOC,MEZJ_LOC,TAUT_LOC
C
      IF ( IREL.LE.2 ) RETURN
C
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
         ALLOCATE (DOS(NKM,NETAB(1),NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DOS')
         IF ( NKM.GT.NKMMAX ) CALL STOP_MESSAGE(ROUTINE,'NKM > MKMAX')
C
         ALLOCATE (MEZZ_LOC(NKMMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MEZZ_LOC')
         ALLOCATE (MEZJ_LOC(NKMMAX,NKMMAX))
         ALLOCATE (TAUT_LOC(NKMMAX,NKMMAX))
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         M = NKMMAX
         N = NKMQ(IQAT(1,IT))
C
         IQ = IQAT(1,IT)
C
C If no magnetization rotation was involved in the calculation,
C   calculate the DOS without rotation
C
         DIRECT = .FALSE.
         IF ( .NOT.ALLOCATED(DROTQ) ) THEN
            DIRECT = .TRUE.
         ELSE IF ( SIZE(DROTQ,3).LE.1 ) THEN
            DIRECT = .TRUE.
         END IF
C
C Rotates into a local coordinate system
C
         IF ( DOSCORSYS.EQ.'LOC' .AND. .NOT.DIRECT ) THEN
C
            CALL ROTATE(TAUT(1,1,IT),'G->L',TAUT_LOC,N,DROTQ(1,1,IQ),M)
            CALL ROTATE(MEZZ(1,1,IT,1),'G->L',MEZZ_LOC,N,DROTQ(1,1,IQ),
     &                  M)
            CALL ROTATE(MEZJ(1,1,IT,1),'G->L',MEZJ_LOC,N,DROTQ(1,1,IQ),
     &                  M)
C
         ELSE
            MEZZ_LOC(1:M,1:M) = MEZZ(1:M,1:M,IT,1)
            MEZJ_LOC(1:M,1:M) = MEZJ(1:M,1:M,IT,1)
            TAUT_LOC(1:M,1:M) = TAUT(1:M,1:M,IT)
         END IF
C
         CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ_LOC,M,TAUT_LOC,M,C0,W1,M)
         DO J = 1,N
            W1(J,J) = W1(J,J) - CPRE*MEZJ_LOC(J,J)
         END DO
C
         IF ( DOSREP.EQ.'RLM' ) THEN
C
C-----------------------------------------------------------------------
C   change representation to non-relativistic / REAL spher. harmonics
C-----------------------------------------------------------------------
C
            CALL CHANGEREP(NKM,NKMMAX,W1,'REL>RLM',W2)
C
         ELSE IF ( DOSREP.EQ.'CLM' ) THEN
C
C-----------------------------------------------------------------------
C   change representation to non-relativistic / COMPLEX sph. harmonics
C-----------------------------------------------------------------------
C
            CALL CHANGEREP(NKM,NKMMAX,W1,'REL>CLM',W2)
C
         ELSE IF ( DOSREP.NE.'REL' ) THEN
            CALL STOP_MESSAGE(ROUTINE,'DOSREP.NE.REL')
         END IF
C
         IF ( DOSREP.NE.'REL' ) THEN
            DO J = 1,N
               DOS(J,IECURR,IT) = DIMAG(W2(J,J))
            END DO
         ELSE
            DO J = 1,N
               DOS(J,IECURR,IT) = DIMAG(W1(J,J))
            END DO
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
      IF ( IECURR.LT.NETAB(1) ) RETURN
C=======================================================================
C
C
      IF ( DOSREP.NE.'REL' ) THEN
         IF ( DOSREP.EQ.'RLM' ) THEN
            STR7 = 'REAL   '
            LL = 4
         ELSE IF ( DOSREP.EQ.'CLM' ) THEN
            STR7 = 'COMPLEX'
            LL = 7
         END IF
         WRITE (6,99005) STR7(1:LL)
      ELSE
         WRITE (6,99006)
      END IF
C
      DO IT = ITBOT,ITTOP
C
         N = NKMQ(IQAT(1,IT))
C
         STR12 = '_DOS_'//DOSREP//'.dat'
C
         CALL OPENFILT(DATSET,LDATSET,TXT_T,LTXT_T,IT,STR12,12,FILNAM,
     &                 LFN,2,IOTMP,'(l,ml,ms)-DOS',13,NTMAX)
C
         IF ( DOSREP.NE.'REL' ) THEN
            WRITE (IOTMP,99001) SYSTEM(1:LSYSTEM),TXT_T(IT),IREL,NLT(IT)
     &                          ,STR7,EFERMI
         ELSE
            WRITE (IOTMP,99002) SYSTEM(1:LSYSTEM),TXT_T(IT),IREL,NLT(IT)
     &                          ,EFERMI
         END IF
C
         DO IE = 1,NETAB(1)
            WRITE (IOTMP,99003) DREAL(ETAB(IE,1)-EFERMI)*RY_EV,
     &                          (DOS(J,IE,IT),J=1,N)
         END DO
         CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C
         IF ( CHECK ) THEN
C
            WRITE (6,99004)
            NH = N/2
C
C--- n(up) = n(dn)
            DO IE = 1,NETAB(1)
               DO J = 1,NH
                  IF ( ABS(DOS(J,IE,IT)-DOS(NH+J,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * up-dn: ',IE,J
               END DO
            END DO
C
C--- n(p) degenerated
            DO IE = 1,NETAB(1)
               IF ( ABS(DOS(2,IE,IT)-DOS(3,IE,IT)).GT.TOL ) WRITE (6,*)
     &               '* WRDOSREP * p 2-3: ',IE
               IF ( ABS(DOS(2,IE,IT)-DOS(4,IE,IT)).GT.TOL ) WRITE (6,*)
     &               '* WRDOSREP * p 2-4: ',IE
            END DO
C
            IF ( DOSREP.EQ.'RLM' ) THEN
               DO IE = 1,NETAB(1)
C--- d(t2g) degenerated
                  IF ( ABS(DOS(5,IE,IT)-DOS(6,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * d(t2g) 5-6: ',IE
                  IF ( ABS(DOS(5,IE,IT)-DOS(8,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * d(t2g) 5-8: ',IE
C--- d(eg) degenerated
                  IF ( ABS(DOS(7,IE,IT)-DOS(9,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * d(eg) 7-9: ',IE
               END DO
            ELSE
               DO IE = 1,NETAB(1)
C--- d(ml/-ml) degenerated
                  IF ( ABS(DOS(5,IE,IT)-DOS(9,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * d(ml/-ml) 5-9: ',IE
                  IF ( ABS(DOS(6,IE,IT)-DOS(8,IE,IT)).GT.TOL )
     &                 WRITE (6,*) '* WRDOSREP * d(ml/-ml) 6-8: ',IE
               END DO
            END IF
         END IF
C-------------------------------------------------------------- CHECK --
      END DO
C
99001 FORMAT ('#',/,'#  ',A,/,'#  (l,ml,ms)-resolved DOS for ',A,/,
     &        '#  relativistic mode  IREL =',I3,/,
     &        '#  l-expansion        NL   =',I3,/,
     &        '#  representation           ',A,/,
     &        '#  Fermi energy       [Ry]  ',F10.6,/,'#',
     &        /'#  data structure:   E, n(l,m,s)  ',
     &        'for  m=-l,l; l=0,lmax; s=dn,up',/'#',20X,
     &        'E in [eV] with respect to E_Fermi',/,'#',20X,
     &        'n in [states/Ry*atom*spin]',/,'#')
99002 FORMAT ('#',/,'#  ',A,/,'#  (l,ml,ms)-resolved DOS for ',A,/,
     &        '#  relativistic mode  IREL =',I3,/,
     &        '#  l-expansion        NL   =',I3,/,
     &        '#  Fermi energy       [Ry]  ',F10.6,/,'#',
     &        /'#  data structure:   E, n(kappa,mu)  ',
     &        'for  mu=-l-1/2,+l+1/2; kappa=-1,..,-l_max-1',/'#',20X,
     &        'E in [eV] with respect to E_Fermi',/,'#',20X,
     &        'n in [states/Ry*atom]',/,'#')
99003 FORMAT (25F10.5)
99004 FORMAT (//,1X,79('*'),/,36X,'<WRDOSREP>',/,1X,79('*'),//,10X,
     &        'checking consistency of (l,ml,ms)-resolved DOS',/,10X,
     &        'assuming paramagnetic and cubic system ',/)
99005 FORMAT (//,1X,79('*'),/,36X,'<WRDOSREP>',/,1X,79('*'),//,10X,
     &        'writing (l,ml,ms)-resolved DOS ',/,10X,
     &        'with respect to ',A,' spherical harmonics ',/)
99006 FORMAT (//,1X,79('*'),/,36X,'<WRDOSREP>',/,1X,79('*'),//,10X,
     &        'writing (kappa,mu)-resolved DOS ',/)
      END
