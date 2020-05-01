C*==me_drive.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE ME_DRIVE(IKMFBOT,IKMFTOP,IFILF,ERYDF,IKMIBOT,IKMITOP,
     &                    IFILI,ERYDI,ME_CC_BRA_RWF,IT,MZAZB,MZBZA,
     &                    MIRR_2,MIRR_3,MIRR_4,C,K_CALC_ME,MEFORM,IWME)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the matrix elements of the                            *
C   *                                                                  *
C   *             ELECTRIC DIPOLE INTERACTION OPERATOR                 *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *                                       ASA            FP          *
C   *                                     REG IRR       REG IRR        *
C   *  MEFORM  NAB  Nabla-form             +   -         -   -         *
C   *          ADA  alpha * A - form       +   +         +   +         *
C   *          GRV  gradient V - form      +   -         -   -         *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  IREL < 3  only Nabla form for non-magnetic systems              *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *    K_CALC_ME = 1:    calculate MZBZA                             *
C   *              = 2:  + calculate MZAZB                             *
C   *              = 3:  + calculate MIRR_2,MIRR_3,MIRR_4              *
C   *                                                                  *
C   *                       ALL polarisation combinations lam-lam'     *
C   *                       are treated for MIRR                       *
C   *                                                                  *
C   *                       MIRR(LAM1,LAM2,LAM3) is summed with        *
C   *                       respect to 2nd (initial state) index       *
C   *                                                                  *
C   *  IWME<>0:  write matrix elements MREG to file IWME               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NPOLMAX
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--ME_DRIVE44
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 C
      COMPLEX*16 ERYDF,ERYDI
      INTEGER IFILF,IFILI,IKMFBOT,IKMFTOP,IKMIBOT,IKMITOP,IT,IWME,
     &        K_CALC_ME
      CHARACTER*3 MEFORM
      LOGICAL ME_CC_BRA_RWF
      COMPLEX*16 MIRR_2(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX),
     &           MIRR_3(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX),
     &           MIRR_4(NKMMAX,NKMMAX,NPOLMAX,NPOLMAX),
     &           MZAZB(NKMMAX,NKMMAX,NPOLMAX),
     &           MZBZA(NKMMAX,NKMMAX,NPOLMAX)
C
C Local variables
C
      INTEGER M,N
      REAL*8 T
C
C*** End of declarations rewritten by SPAG
C
      IF ( FULLPOT ) STOP
C
C ======================================================================
C      non-magnetic systems in the non- or scalar relativistic mode
C ======================================================================
      IF ( IREL.LE.2 ) THEN
C
         IF ( IREL.EQ.2 ) STOP 
     &                      '<> matrix elements not avalable for IREL=2'
C
         STOP
C         CALL MENABIRR_SRA(IFILF,IFILI,IT,MREG,MBAR,MIRR       )
C
C         RETURN
C
      END IF
C ======================================================================
C
      IF ( IPRINT.GT.0 ) WRITE (6,99001) MEFORM(1:3)
C ================================================================== ADA
C
      IF ( MEFORM.EQ.'ADA' ) THEN
C
         CALL ME_ALF_ALF(IKMFBOT,IKMFTOP,IFILF,ERYDF,IKMIBOT,IKMITOP,
     &                   IFILI,ERYDI,ME_CC_BRA_RWF,IT,MZAZB,MZBZA,
     &                   MIRR_2,MIRR_3,MIRR_4,C,K_CALC_ME)
C
C ================================================================== NAB
C
      ELSE IF ( MEFORM.EQ.'NAB' ) THEN
C
         CALL ME_NAB_NAB(IKMFBOT,IKMFTOP,IFILF,ERYDF,IKMIBOT,IKMITOP,
     &                   IFILI,ERYDI,ME_CC_BRA_RWF,IT,MZAZB,MZBZA,
     &                   MIRR_2,MIRR_3,MIRR_4,C,K_CALC_ME)
C
C ================================================================== GRV
C
      ELSE IF ( MEFORM.EQ.'GRV' ) THEN
C
Cc         CALL ME_GRV_GRV(IKMFBOT,IKMFTOP,IFILF,ERYDF,IKMIBOT,IKMITOP,IFILI,
Cc     &              ERYDI,ME_CC_BRA_RWF,IT,MZAZB,MZBZA,MIRR_2,MIRR_3,
Cc     &              MIRR_4,C)
         WRITE (6,*) ' MEFORM = ',MEFORM
         STOP 'in <MECALC>:  not implemented'
C
C ================================================================== ???
C
      ELSE
         WRITE (6,*) ' MEFORM = ',MEFORM
         STOP 'in <MECALC>:  not implemented'
      END IF
C
C=======================================================================
C
      IF ( IWME.LE.0 ) RETURN
C
      N = NKM
      M = NKMMAX
      T = 1D-8
C
      CALL CMATSTRUCT('ME (+) '//MEFORM,MZAZB(1,1,1),N,M,3,3,0,T,IWME)
      CALL CMATSTRUCT('ME (-) '//MEFORM,MZAZB(1,1,2),N,M,3,3,0,T,IWME)
      CALL CMATSTRUCT('ME (z) '//MEFORM,MZAZB(1,1,3),N,M,3,3,0,T,IWME)
C
C=======================================================================
99001 FORMAT (//,1X,79('*'),/,36X,'<MECALC>',/,1X,79('*'),//,10X,
     &        'setting up matrix elements for MEFORM = ',A,/)
      END
