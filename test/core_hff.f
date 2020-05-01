C*==core_hff.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CORE_HFF(IM,KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,RC,DRDIC,
     &                    RNUC,IR_ZERO,IRTOP,NRC)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
C   *                  CURRENT  CORE STATE S                           *
C   *                                                                  *
C   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:B_AU2CGS
      IMPLICIT NONE
C*--CORE_HFF15
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CORE_HFF')
C
C Dummy arguments
C
      INTEGER IM,IRTOP,IR_ZERO,KAP1,KAP2,NRC,NSOL,S
      REAL*8 MJ,RNUC
      REAL*8 BHF(2,2),DRDIC(NRC),FCK(2,2,NRC),GCK(2,2,NRC),RC(NRC)
C
C Local variables
C
      REAL*8 AMEHF(2,2),XX(5),YI(:),YY(5),ZI(:)
      INTEGER I,IA_ERR,IR,K1,K2
      LOGICAL RNON0
      REAL*8 YLAG
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE YI,ZI
C
      ALLOCATE (YI(NRC),ZI(NRC),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZI')
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   angular hyperfine matrix elements   see e.g.  E.M.Rose
C        the factor  i  has been omitted
C
      AMEHF(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
      IF ( NSOL.EQ.2 ) THEN
         AMEHF(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
         AMEHF(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
         AMEHF(2,1) = AMEHF(1,2)
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
C
            DO IR = 1,IR_ZERO
               YI(IR) = DRDIC(IR)
     &                  *(GCK(K1,S,IR)*FCK(K2,S,IR)+FCK(K1,S,IR)
     &                  *GCK(K2,S,IR))
            END DO
            YI((IR_ZERO+1):NRC) = 0.0D0
C
            CALL RRADINT_R(IM,YI,ZI)
C
            IF ( RNON0(RNUC) ) THEN
               DO I = 1,5
                  XX(I) = RC(IRTOP-5+I)
                  YY(I) = ZI(IRTOP-5+I)
               END DO
               ZI(IRTOP) = YLAG(RNUC,XX,YY,0,4,5)
            END IF
C
            XX(1) = 1.0D0
C                      !RC( 1)
            XX(2) = 6.0D0
C                      !RC( 6)
            XX(3) = 11.0D0
C                      !RC(11)
            YY(1) = ZI(IRTOP) - ZI(1)
            YY(2) = ZI(IRTOP) - ZI(6)
            YY(3) = ZI(IRTOP) - ZI(11)
C
            BHF(K1,K2) = B_AU2CGS*AMEHF(K1,K2)*YLAG(0.0D0,XX,YY,0,2,3)
C
         END DO
      END DO
C
      END
