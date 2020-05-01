C*==init_mod_sites.f    processed by SPAG 6.70Rc at 10:23 on 31 Oct 2016
      SUBROUTINE INIT_MOD_SITES
C   ********************************************************************
C   *                                                                  *
C   *  initialize all tables dependent ONLY on                         *
C   *                                                                  *
C   *      NQ,  NT  and  NL   and derived array dimensions             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:QMPHI,QMTET,QMGAM,MDIRQ,NOQ,ICPA,IMQ,NQMAX,
     &    QBAS,IQBOT,IQTOP,NQ,RNNQ,NNNQ,QMPHI_POTFIL,QMTET_POTFIL,
     &    QMGAM_POTFIL,MAGROT_Q,KMAD1M_Q
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IA_ERR
C
C*** End of declarations rewritten by SPAG
C
C---------------------------------------- variables depending only on NQ
C
      ALLOCATE (QMPHI(NQMAX),QMPHI_POTFIL(NQMAX))
      ALLOCATE (QMTET(NQMAX),QMTET_POTFIL(NQMAX))
      ALLOCATE (QMGAM(NQMAX),QMGAM_POTFIL(NQMAX))
      QMPHI(:) = 0D0
      QMTET(:) = 0D0
      QMGAM(:) = 0D0
      QMPHI_POTFIL(:) = 0D0
      QMTET_POTFIL(:) = 0D0
      QMGAM_POTFIL(:) = 0D0
C
      ALLOCATE (NOQ(NQMAX),ICPA(NQMAX),IMQ(NQMAX))
      ALLOCATE (MAGROT_Q(NQMAX),MDIRQ(3,NQMAX))
      ALLOCATE (QBAS(3,NQMAX),KMAD1M_Q(0:NQMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: INIT_MOD_SITES -> KMAD1M_Q'
C
      MAGROT_Q(1:NQMAX) = .FALSE.
      KMAD1M_Q(0:NQMAX) = .FALSE.
C
C-----------------------------------------------------------------------
C     loops over sites IQ run from    IQBOT    to    IQTOP
C-----------------------------------------------------------------------
      IQBOT = 1
      IQTOP = NQ
C
C-----------------------------------------------------------------------
C        parameters for screened-impurity-model in case of CPA
C-----------------------------------------------------------------------
      ALLOCATE (RNNQ(NQMAX),NNNQ(NQMAX))
      RNNQ(1:NQMAX) = 999999D0
      NNNQ(1:NQMAX) = 999999
C
      END
