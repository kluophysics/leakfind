C*==xrayclb.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE XRAYCLB(NCV,KMAT,IT,ZG,ZF,ZNORM_CST,GCOR,FCOR,IRTOP,
     &                   RPWL,AMECIG,AMECIF,IKMCOR,NSOLCOR,NCST,NCSTMAX,
     &                   IFILCBWF_CST)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the Coulomb integrals ( core  - valence )             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NLMAX,NKMMAX
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C*--XRAYCLB14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IRTOP,IT,NCST,NCSTMAX,NCV
      REAL*8 AMECIF(NKMMAX,NKMMAX,2*NLMAX),AMECIG(NKMMAX,NKMMAX,2*NLMAX)
     &       ,FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),
     &       RPWL(NRMAX,0:2*NLMAX),ZNORM_CST(NKMMAX,NCSTMAX)
      INTEGER IFILCBWF_CST(NCSTMAX),IKMCOR(NCSTMAX,2),NSOLCOR(NCSTMAX)
      COMPLEX*16 KMAT(NCV,NCV),ZF(NRMAX,2,NKMMAX,NCSTMAX),
     &           ZG(NRMAX,2,NKMMAX,NCSTMAX)
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
      CALL STOP_MESSAGE('XRAYCLB','routine not available at the moment')
C=======================================================================
C
      WRITE (*,*) NCV,KMAT,IT,ZG,ZF,ZNORM_CST,GCOR,FCOR,IRTOP,RPWL,
     &            AMECIG,AMECIF,IKMCOR,NSOLCOR,NCST,NCSTMAX,IFILCBWF_CST
C
C
      END
