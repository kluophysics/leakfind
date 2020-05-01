C*==init_mod_rmesh.f    processed by SPAG 6.70Rc at 08:28 on 14 Jan 2016
      SUBROUTINE INIT_MOD_RMESH
C   ********************************************************************
C   *                                                                  *
C   *  initialize all tables dependent ONLY on                         *
C   *                                                                  *
C   *      the radial mesh                                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:JRMT,JRWS,ARMSH,BRMSH,DX,RWS,RMT,R,DRDI,R2DRDI,
     &     DRDIOVR,NRMAX,NRCMAX,NMMAX,RAUX1,RAUX2,CAUX1,CAUX2,
     &     R_COR,DRDI_COR,DRDIOVR_COR
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (ARMSH(NMMAX),BRMSH(NMMAX),DX(NMMAX))
      ALLOCATE (JRWS(NMMAX),JRMT(NMMAX),RWS(NMMAX),RMT(NMMAX))
C
      ARMSH = 999999D0
      BRMSH = 999999D0
C
C-------------------------------------- variables depending on NM and NR
C
      ALLOCATE (R(NRMAX,NMMAX))
      ALLOCATE (DRDI(NRMAX,NMMAX))
      ALLOCATE (R2DRDI(NRMAX,NMMAX))
      ALLOCATE (DRDIOVR(NRMAX,NMMAX))
      ALLOCATE (R_COR(NRCMAX,NMMAX))
      ALLOCATE (DRDI_COR(NRCMAX,NMMAX))
      ALLOCATE (DRDIOVR_COR(NRCMAX,NMMAX))
      ALLOCATE (RAUX1(NRMAX),RAUX2(NRMAX))
      ALLOCATE (CAUX1(NRMAX),CAUX2(NRMAX))
C
      END
