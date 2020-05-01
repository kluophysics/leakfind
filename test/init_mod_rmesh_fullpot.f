C*==init_mod_rmesh_fullpot.f    processed by SPAG 6.70Rc at 17:02 on  5 Feb 2016
      SUBROUTINE INIT_MOD_RMESH_FULLPOT(NTLIM,NT0,NA_TAUX,RMTRED0,
     &                                  NSFLIM,FULLPOT,FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *  initialize all tables dependent ONLY on                         *
C   *                                                                  *
C   *      the radial mesh                                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRCRI,JRNS1,JRNSMIN,NPAN,NRNS,NRSFTOT,NSF,
     &    FLMSF,ISFLM,JRCUT,KLMSF,LMISF,NMMAX,WINTLM,NSFMAX,NPANMAX,
     &    NLMSFMAX,NRSFMAX,SPHERCELL
      USE MOD_SITES,ONLY:NQCLU,IQAT,ITOQ
      USE MOD_ANGMOM,ONLY:NL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FULLPOT,FULLPOT5
      INTEGER NSFLIM,NT0,NTLIM
      INTEGER NA_TAUX(NTLIM)
      REAL*8 RMTRED0(NMMAX)
C
C Local variables
C
      INTEGER NPANLIM,NRSFLIM
      LOGICAL SFN_AVAILABLE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (JRCRI(NMMAX),JRNS1(NMMAX),NPAN(NMMAX))
      ALLOCATE (NRNS(NMMAX),NRSFTOT(NMMAX),NSF(NMMAX))
C
C======================================================================
C    dummy allocation for spherical ASA or muffin-tin geometry case
C======================================================================
C
      IF ( .NOT.(FULLPOT .OR. FULLPOT5) ) THEN
C
         NRSFMAX = 1
         NSFMAX = 1
         NPANMAX = 2
C
         ALLOCATE (FLMSF(NRSFMAX,NSFMAX,NMMAX),JRCUT(0:NPANMAX,NMMAX))
         ALLOCATE (WINTLM(NRSFMAX,NSFMAX,NMMAX),ISFLM(NLMSFMAX,NMMAX))
         ALLOCATE (KLMSF(NLMSFMAX,NMMAX),LMISF(NLMSFMAX,NMMAX))
C
C------------------- set shape function parameters for spherical region
C
         NPAN(:) = 1
         NSF(:) = 1
         LMISF(1,:) = 1
         KLMSF(1,:) = 1
C
         RETURN
C
      END IF
C
C======================================================================
C                     FULL POTENTIAL CASE
C======================================================================
C
      IF ( SPHERCELL ) THEN
C
         NSFMAX = 1
         NPANMAX = 2
         NRSFMAX = 75
C
      ELSE
C
         IF ( NQCLU.EQ.0 ) THEN
            NSFMAX = NSFLIM
         ELSE
            NSFMAX = (4*NL-3)**2
         END IF
C
         CALL SFNREADDIM(NPANLIM,NRSFLIM,SFN_AVAILABLE)
C
         IF ( .NOT.SFN_AVAILABLE ) CALL SFNCREATE(NT0,IQAT,ITOQ,NA_TAUX,
     &        RMTRED0,NSFLIM,NPANLIM,NRSFLIM)
C
         NSFMAX = NSFLIM
         NPANMAX = NPANLIM
         NRSFMAX = NRSFLIM
C
      END IF
C
      JRNSMIN = 1
C
      ALLOCATE (FLMSF(NRSFMAX,NSFMAX,NMMAX),JRCUT(0:NPANMAX,NMMAX))
      ALLOCATE (WINTLM(NRSFMAX,NSFMAX,NMMAX),ISFLM(NLMSFMAX,NMMAX))
      ALLOCATE (KLMSF(NLMSFMAX,NMMAX),LMISF(NLMSFMAX,NMMAX))
C
      JRCUT(:,:) = 0
      ISFLM(:,:) = 0
      KLMSF(:,:) = 0
      LMISF(:,:) = 0
      WINTLM(:,:,:) = 0D0
      FLMSF(:,:,:) = 0D0
C
      END
