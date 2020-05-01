C*==sig1_surf_nonpub.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIG1_SURF_NONPUB(SUMX_1,SUMX_2,SIG1IRQ,IQ,JQ)
C   ********************************************************************
C   *                                                                  *
C   *   calculate  TRACE jbar(mue,z2,z1)*mat*jbar(nue,z1,z2)           *
C   *   (including Vertex-corrections)                                 *
C   *                                                                  *
C   *   or  TRACE jbar(mue,z2,z1)*chi*jbar(nue,z1,z2)                  *
C   *   (neglecting Vertex-corrections)                                *
C   *                                                                  *
C   *  NOTE: CHIZ is defined only for the regime                       *
C   *        IQ = IQBOT_CHI, ... , IQTOP_CHI                           *
C   *        the auxilary site index IQCHI = IQ - IQBOT_CHI + 1        *
C   *        is used to index this regime with IQCHI = 1, ..., NQ_CHI  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:NQ
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--SIG1_SURF_NONPUB21
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IQ,JQ
      COMPLEX*16 SIG1IRQ(3,3,NQ,NQ),SUMX_1(3,3,2,2),SUMX_2(3,3,2,2)
C
C Local variables
C
      INTEGER MUE,NUE
      COMPLEX*16 TMP_1(3,3),TMP_2(3,3)
C
C*** End of declarations rewritten by SPAG
C
C-------------------------------- i ( j_m Im G j_n - j_n Im G j_m ) Re G
C------------------------------------------ antisymmetric part of tensor
C
      DO MUE = 1,3
         DO NUE = 1,3
            TMP_1(MUE,NUE) = CI*(0.25D0/CI)
     &                       *(SUMX_1(MUE,NUE,1,1)+SUMX_1(MUE,NUE,1,2)
     &                       -SUMX_1(MUE,NUE,2,1)-SUMX_1(MUE,NUE,2,2))
C
            TMP_2(MUE,NUE) = CI*(0.25D0/CI)
     &                       *(SUMX_2(MUE,NUE,1,1)+SUMX_2(MUE,NUE,1,2)
     &                       -SUMX_2(MUE,NUE,2,1)-SUMX_2(MUE,NUE,2,2))
C
         END DO
      END DO
      DO MUE = 1,3
         DO NUE = 1,3
C
            SIG1IRQ(MUE,NUE,IQ,JQ) = TMP_1(MUE,NUE) - TMP_2(NUE,MUE)
C
            IF ( ABS(SIG1IRQ(MUE,NUE,IQ,JQ)).LT.1D-12 )
     &           SIG1IRQ(MUE,NUE,IQ,JQ) = C0
            IF ( ABS(DIMAG(SIG1IRQ(MUE,NUE,IQ,JQ))).LT.1D-12 )
     &           SIG1IRQ(MUE,NUE,IQ,JQ) = DREAL(SIG1IRQ(MUE,NUE,IQ,JQ))
         END DO
      END DO
C
C
C--------------- only 1/2 of SIG1IRQ --> Crepieux formula PRB 64, 014416
C--------------- also in the Bastin framework this factor has to be
C--------------- included
C
      SIG1IRQ(:,:,IQ,JQ) = 0.5D0*SIG1IRQ(:,:,IQ,JQ)
C
      END
