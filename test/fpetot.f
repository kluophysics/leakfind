C*==fpetot.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPETOT(ECOU,EFERMI,EPOTIN,EXC,ETOT)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the total energy                                      *
C   *  gather all energy-parts which are calculated in different       *
C   *  subroutines.                                                    *
C   *  since the program uses group theory only shell-indices          *
C   *  are used instead of atom-indices.                               *
C   *                                                                  *
C   *                            b.drittler   may 1987                 *
C   *                                                                  *
C   *  Note the modified definition of  Espv, Espc                     *
C   *       ( they are no longer spin-resolved )                       *
C   *  Espv(1, ) , Espc(1, ) is the l=0 component, and so on           *
C   *  The core states are assumed to have only  s,p,d,f               *
C   *  but no higher l-components, that is  nlcore=4                   *
C   *  The total energies of the indivdual atom-types are              *
C   *  concentrationally averaged to get the total energy              *
C   *  in case of CPA-calculation                                      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CALCMODE,ONLY:BREAKPOINT,SOLVER_FP
      USE MOD_SITES,ONLY:QMPHI,QMTET,IQAT
      USE MOD_ANGMOM,ONLY:NLMAX,IBND,NLCORE
      USE MOD_TYPES,ONLY:TXT_T,LTXT_T,CONC,NAT,ITBOT,ITTOP,NTMAX,NLFP,
     &    NLFPMAX,OBS_LT,ECOR_LT
      USE MOD_FILES,ONLY:IFILBUILDBOT,WRBUILDBOT
      IMPLICIT NONE
C*--FPETOT30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='FPETOT')
C
C Dummy arguments
C
      REAL*8 EFERMI,ETOT
      REAL*8 ECOU(0:(NLFPMAX-1),NTMAX),EPOTIN(NTMAX),
     &       EXC(0:(NLFPMAX-1),NTMAX)
C
C Local variables
C
      REAL*8 ET,ETSP
      INTEGER IL,IQ,IT,L
      CHARACTER*2 TEXTL(0:5)
C
C*** End of declarations rewritten by SPAG
C
      DATA TEXTL/'s:','p:','d:','f:','g:','h:'/
C
      ETOT = 0.0D0
      WRITE (6,99001) EFERMI
C-----------------------------------------------------------------------
C     loop over atom types
C-----------------------------------------------------------------------
C
      DO IT = ITBOT,ITTOP
C
         IQ = IQAT(1,IT)
C
         WRITE (6,99002) IT,TXT_T(IT)(1:LTXT_T(IT))
         ET = 0.0D0
         WRITE (6,99003) (TEXTL(IL-1),ECOR_LT(IL,IT),IL=1,NLCORE)
         WRITE (6,99004) (TEXTL(IL-1),OBS_LT(0,IBND,IL,IT),IL=1,NLMAX)
         DO IL = 1,NLMAX
            ET = ET + OBS_LT(0,IBND,IL,IT)
         END DO
         DO IL = 1,NLCORE
            ET = ET + ECOR_LT(IL,IT)
         END DO
C
         ETSP = ET
C
         WRITE (6,99005) (L,ECOU(L,IT),L=0,(NLFP-1))
         WRITE (6,99006) (L,EXC(L,IT),L=0,(NLFP-1))
C
         DO L = 0,(NLFP-1)
            ET = ET + ECOU(L,IT) + EXC(L,IT)
         END DO
C
         ET = ET + EPOTIN(IT)
C
         WRITE (6,99007) ETSP,ET
         IF ( BREAKPOINT.NE.0 ) WRITE (6,99008) 'BREAK ',NINT(QMTET(IQ))
     &                                 ,NINT(QMPHI(IQ)),SOLVER_FP,ETSP,
     &                                 'BREAK ',NINT(QMTET(IQ)),
     &                                 NINT(QMPHI(IQ)),SOLVER_FP,ET
C
         ETOT = ETOT + ET*CONC(IT)*NAT(IT)
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99010)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),IT,
     &                            SUM(ECOR_LT(1:NLCORE,IT)),
     &                            SUM(ECOU(0:(NLFP-1),IT)),
     &                            SUM(EXC(0:(NLFP-1),IT)),ETOT
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
      END DO
C
      WRITE (6,99009) ETOT
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
99001 FORMAT (/,1X,79('*'),/,34X,'total energies',/,1X,79('*'),//,1x,
     &        'Fermi energy',F10.6,' Ry',/)
99002 FORMAT (1x,'contribution to total energy of atom type IT=',I3,3X,
     &        A,/)
99003 FORMAT (1x,'core    ',5(1x,a,f14.6))
99004 FORMAT (1x,'valence ',5(1x,a,f14.6),:,/,10X,5(1x,a,f14.6))
99005 FORMAT (1x,'Coulomb ',4(i2,':',f14.6),:,/,(9X,4(i2,':',f14.6)))
99006 FORMAT (1x,'ex.-cor.',4(i2,':',f14.6),:,/,(9X,4(i2,':',f14.6)))
99007 FORMAT (/,1x,'sum of single particle energies:',f18.10,/,1x,
     &        'total contribution of the atom: ',f18.10,/)
99008 FORMAT (/,A,1x,'sum of single particle energies:',2I3,2X,A,F16.10,
     &        /,A,1x,'total contribution of the atom: ',2I3,A,2X,F16.10,
     &        /)
99009 FORMAT (1x,'total energy                    ',f18.10,' Ry')
99010 FORMAT ('# BUILDBOT: ',A,':  ECOR ECOU EXC (+ETOT) for IT =',I5,/,
     &        (1PE22.14))
      END
