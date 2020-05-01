C*==ame_nab.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AME_NAB
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed                   *
C   *                                                                  *
C   *   for the          NABLA * A -  form                             *
C   *                                                                  *
C   *                       ->   ->                                    *
C   *          AN1 = < LAM| NAB . e(POL) | LAM'>     l = l'+1          *
C   *          AN2 = < LAM| NAB . e(POL) | LAM'>     l = l'-1          *
C   *          AV1 = < LAM| SIG . e(POL) |-LAM'>                       *
C   *          AV2 = <-LAM| SIG . e(POL) | LAM'>                       *
C   *          AB1 = < LAM| SIG x e(POL) |-LAM'>/i                     *
C   *          AB2 = <-LAM| SIG x e(POL) | LAM'>/i                     *
C   *                                                                  *
C   *   IPOL  = 1,2,3  ==  (-),(z),(+):                                *
C   *                                                                  *
C   *  ->         ->       ->                                          *
C   *  e(-) =  + [e(x) - i e(y)] / sqrt(2)    RCP / negative helicity  *
C   *                                                                  *
C   *  ->         ->       ->                                          *
C   *  e(+) =  - [e(x) + i e(y)] / sqrt(2)    LCP / positive helicity  *
C   *                                                                  *
C   *  definition of polarisation revised 2015                         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:AN1_NAB,AN2_NAB,AV1_NAB,AV2_NAB,AB1_NAB,
     &    AB2_NAB,CGC,AME_G,L_IKM,LB_IKM,IMKM_IKM,MUEM05_IKM,ISMT,
     &    NKM_EXT
      USE MOD_CONSTANTS,ONLY:SQRT_2
      IMPLICIT NONE
C*--AME_NAB34
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 CGC1
      INTEGER IKM,IKM1,IKM2,IMKM,IMKM1,IMKM2,IPOL,JKM,JMKM,L1,L2,LB1,
     &        LB2,MM,MUE1M05,MUE2M05
C
C*** End of declarations rewritten by SPAG
C
      IF ( .NOT.ALLOCATED(AN1_NAB) ) THEN
         ALLOCATE (AN1_NAB(NKM_EXT,NKM_EXT,3))
         ALLOCATE (AN2_NAB(NKM_EXT,NKM_EXT,3))
         ALLOCATE (AV1_NAB(NKM_EXT,NKM_EXT,3))
         ALLOCATE (AV2_NAB(NKM_EXT,NKM_EXT,3))
         ALLOCATE (AB1_NAB(NKM_EXT,NKM_EXT,3))
         ALLOCATE (AB2_NAB(NKM_EXT,NKM_EXT,3))
      END IF
C
      AN1_NAB(:,:,:) = 0D0
      AN2_NAB(:,:,:) = 0D0
      AV1_NAB(:,:,:) = 0D0
      AV2_NAB(:,:,:) = 0D0
      AB1_NAB(:,:,:) = 0D0
      AB2_NAB(:,:,:) = 0D0
C
C-----------------------------------------------------------------------
C              AV1 = < LAM| SIG . e(POL) |-LAM'>
C              AV2 = <-LAM| SIG . e(POL) | LAM'>
C-----------------------------------------------------------------------
C
      DO IKM = 1,NKM_EXT
         IMKM = IMKM_IKM(IKM)
         DO JKM = 1,NKM_EXT
            JMKM = IMKM_IKM(JKM)
C
            AV1_NAB(IKM,JKM,:) = +AME_G(IKM,JMKM,:,ISMT)
            AV2_NAB(IKM,JKM,:) = +AME_G(IMKM,JKM,:,ISMT)
C
         END DO
      END DO
C
C-----------------------------------------------------------------------
C
      DO IKM1 = 1,NKM_EXT
         L1 = L_IKM(IKM1)
         LB1 = LB_IKM(IKM1)
         MUE1M05 = MUEM05_IKM(IKM1)
         IMKM1 = IMKM_IKM(IKM1)
C
         DO IKM2 = 1,NKM_EXT
            L2 = L_IKM(IKM2)
            LB2 = LB_IKM(IKM2)
            MUE2M05 = MUEM05_IKM(IKM2)
            IMKM2 = IMKM_IKM(IKM2)
C
            MM = MUE1M05 - MUE2M05
C
            IF ( ABS(MM).GT.1 ) CYCLE
C
            IPOL = MM + 2
C
C-----------------------------------------------------------------------
C              AN1 = < LAM| NAB . e(POL) | LAM'>
C              AN2 = < LAM| NAB . e(POL) | LAM'>
C-----------------------------------------------------------------------
C
            IF ( L1.EQ.L2+1 ) THEN
C
               AN1_NAB(IKM1,IKM2,IPOL) = SQRT(DBLE(L2+1)/DBLE(2*L2+3))
     &            *(CGC(IKM1,1)*CGC(IKM2,1)*CGC1(L2,L2+1,(MUE2M05+1),MM)
     &            +CGC(IKM1,2)*CGC(IKM2,2)*CGC1(L2,L2+1,(MUE2M05+0),MM))
C
            ELSE IF ( L1.EQ.L2-1 ) THEN
C
               AN2_NAB(IKM1,IKM2,IPOL) = SQRT(DBLE(L2)/DBLE(2*L2-1))
     &            *(CGC(IKM1,1)*CGC(IKM2,1)*CGC1(L2,L2-1,(MUE2M05+1),MM)
     &            +CGC(IKM1,2)*CGC(IKM2,2)*CGC1(L2,L2-1,(MUE2M05+0),MM))
C
            END IF
C
C-----------------------------------------------------------------------
C              AB1 = < LAM| SIG x e(POL) |-LAM'>/i
C              AB2 = <-LAM| SIG x e(POL) | LAM'>/i
C-----------------------------------------------------------------------
C
            IF ( L1.EQ.LB2 ) THEN
C
               IF ( MM.EQ.+1 ) AB1_NAB(IKM1,IKM2,IPOL)
     &              = -SQRT_2*CGC(IKM1,2)*CGC(IMKM2,1)
C
               IF ( MM.EQ.-1 ) AB1_NAB(IKM1,IKM2,IPOL)
     &              = -SQRT_2*CGC(IKM1,1)*CGC(IMKM2,2)
C
            END IF
C
            IF ( LB1.EQ.L2 ) THEN
C
               IF ( MM.EQ.+1 ) AB2_NAB(IKM1,IKM2,IPOL)
     &              = -SQRT_2*CGC(IMKM1,2)*CGC(IKM2,1)
C
               IF ( MM.EQ.-1 ) AB2_NAB(IKM1,IKM2,IPOL)
     &              = -SQRT_2*CGC(IMKM1,1)*CGC(IKM2,2)
C
            END IF
C-----------------------------------------------------------------------
C
         END DO
      END DO
C
      END
