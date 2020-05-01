C*==ame_init.f    processed by SPAG 6.70Rc at 09:00 on  8 Mar 2017
      SUBROUTINE AME_INIT(MEFORM,IWME)
C   ********************************************************************
C   *                                                                  *
C   *  initialize the angular matrix elements of the                   *
C   *                                                                  *
C   *             ELECTRIC DIPOLE INTERACTION OPERATOR                 *
C   *                                                                  *
C   *  occuring in electron spectroscopy                               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,A1_ADA,A2_ADA,B1_ADA,B2_ADA
      USE MOD_RMESH,ONLY:FULLPOT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='AME_INIT')
C
C Dummy arguments
C
      INTEGER IWME
      CHARACTER*3 MEFORM
C
C Local variables
C
      INTEGER M,N,P
      CHARACTER*8 S
      REAL*8 T
      CHARACTER*3 TXTPOL(3)
C
C*** End of declarations rewritten by SPAG
C
      DATA TXTPOL/'(-)','(z)','(+)'/
C
      WRITE (6,99001) MEFORM
C
      N = NKM
      M = NKMMAX
      T = 1D-8
C
C ------------------------------------ calculate angular matrix elements
C
      CALL MEADA_AMEADA_OLD_VERSION
C
C --------------------- ADA - form of non-dipole angular matrix elements
C
      CALL AME_B_ADA
C
C ================================================================== ADA
C
      IF ( MEFORM.EQ.'ADA' ) THEN
C
         CALL MEADA_AMEADA_OLD_VERSION
C
C ================================================================== NAB
C
      ELSE IF ( MEFORM.EQ.'NAB' ) THEN
C
         CALL AME_NAB
C
C ================================================================== GRV
C
      ELSE IF ( MEFORM.EQ.'GRV' ) THEN
C
         IF ( FULLPOT ) CALL STOP_MESSAGE(ROUTINE,
     &       'the chosen TASK does not YET work for FULL POTENTIAL mode'
     &       )
C
         CALL AME_GRV
C
C ================================================================== ???
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'MEFORM='//MEFORM//
     &                     ' not implemented')
      END IF
C
C=======================================================================
C
      IF ( IWME.EQ.0 ) RETURN
C
C=======================================================================
C
      DO P = 1,3
         S = TXTPOL(P)//'  '//MEFORM
         CALL RMATSTRUCT('AME  A1'//S,A1_ADA(1,1,P),N,M,3,3,0,T,6)
         CALL RMATSTRUCT('AME  A2'//S,A2_ADA(1,1,P),N,M,3,3,0,T,6)
         CALL RMATSTRUCT('AME  B1'//S,B1_ADA(1,1,P,2),N,M,3,3,0,T,6)
         CALL RMATSTRUCT('AME  B2'//S,B2_ADA(1,1,P,2),N,M,3,3,0,T,6)
      END DO
C
99001 FORMAT (//,1X,79('*'),/,36X,'<AMEINIT>',/,1X,79('*'),//,10X,
     &        'setting up angular matrix elements for MEFORM = ',A,/)
      END
C*==ame_b_ada.f    processed by SPAG 6.70Rc at 09:00 on  8 Mar 2017
      SUBROUTINE AME_B_ADA
C   ********************************************************************
C   *                                                                  *
C   *   calculate the NON-dipolar part of ADA angular matrix elements  *
C   *                                                                  *
C   *   B(LAM,LAM';POL,LAM) = <chi(LAM)| sig_POL r_LAM |chi(LAM')>     *
C   *                                                                  *
C   *   IPOL= 1,2,3  ==  (-),(0),(+):  polarisation of light           *
C   *                                                                  *
C   *   LAMQ= 1.2.3  ==  (-),(0),(+):  orientation of q-vector         *
C   *                                                                  *
C   *   given in terms of r_LAM concerning the dot product q.r         *
C   *   using spherical coordinates leading to REAL AMEs               *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:B1_ADA,B2_ADA,NKMMAX,NLMAX,CGC
      USE MOD_CONSTANTS,ONLY:SQRT_2
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 CCCDN,CCCUP,FAC
      REAL*8 CGC1
      INTEGER IKM1,IKM2,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,L2,
     &        LAMQ,LB1,LB2,MQ,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
C
      NK = 2*NLMAX - 1
      IF ( .NOT.ALLOCATED(B1_ADA) ) ALLOCATE (B1_ADA(NKMMAX,NKMMAX,3,3),
     &     B2_ADA(NKMMAX,NKMMAX,3,3))
C
      B1_ADA(1:NKMMAX,1:NKMMAX,1:3,1:3) = 0D0
      B2_ADA(1:NKMMAX,1:NKMMAX,1:3,1:3) = 0D0
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
            LB1 = L1 - 1
         ELSE
            KAP1 = -L1 - 1
            LB1 = L1 + 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
            IMKM1 = LB1*2*J1P05 + J1P05 + MUE1M05 + 1
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
                  LB2 = L2 - 1
               ELSE
                  KAP2 = -L2 - 1
                  LB2 = L2 + 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
                  IMKM2 = LB2*2*J2P05 + J2P05 + MUE2M05 + 1
C
C ----------------------------------------------------------------------
                  DO LAMQ = 1,3
C
                     MQ = LAMQ - 2
C
                     IF ( IABS(L1-LB2).EQ.1 ) THEN
C
                        FAC = DSQRT((2.D0*LB2+1)/(2.D0*L1+1))
     &                        *CGC1(LB2,L1,0,0)
C
C                                                                     (-)
                        IF ( (MUE1M05-MUE2M05).EQ.(-1-MQ) )
     &                       B1_ADA(IKM1,IKM2,1,LAMQ) = FAC*CGC(IKM1,1)
     &                       *CGC(IMKM2,2)*CGC1(LB2,L1,(MUE2M05+0),MQ)
     &                       *SQRT_2
C
C                                                                     (0)
                        IF ( (MUE1M05-MUE2M05).EQ.(0-MQ) ) THEN
                           CCCDN = CGC(IKM1,1)*CGC(IMKM2,1)
     &                             *CGC1(LB2,L1,(MUE2M05+1),MQ)
                           CCCUP = CGC(IKM1,2)*CGC(IMKM2,2)
     &                             *CGC1(LB2,L1,(MUE2M05+0),MQ)
                           B1_ADA(IKM1,IKM2,2,LAMQ) = FAC*(CCCUP-CCCDN)
C
                           B1_ADA(IKM1,IKM2,2,LAMQ) = FAC*CGC(IKM1,1)
     &                        *CGC(IMKM2,2)*CGC1(LB2,L1,(MUE1M05+1),MQ)
     &                        *SQRT_2
C
                        END IF
C
C                                                                     (+)
                        IF ( (MUE1M05-MUE2M05).EQ.(+1-MQ) )
     &                       B1_ADA(IKM1,IKM2,3,LAMQ) = -FAC*CGC(IKM1,2)
     &                       *CGC(IMKM2,1)*CGC1(LB2,L1,(MUE2M05+1),MQ)
     &                       *SQRT_2
                     END IF
C
                     IF ( IABS(LB1-L2).EQ.1 ) THEN
C
                        FAC = DSQRT((2.D0*L2+1)/(2.D0*LB1+1))
     &                        *CGC1(L2,LB1,0,0)
C
C                                                                     (-)
                        IF ( (MUE1M05-MUE2M05).EQ.(-1-MQ) )
     &                       B2_ADA(IKM1,IKM2,1,LAMQ) = FAC*CGC(IMKM1,1)
     &                       *CGC(IKM2,2)*CGC1(L2,LB1,(MUE2M05+0),MQ)
     &                       *SQRT_2
C
C                                                              (0)
                        IF ( (MUE1M05-MUE2M05).EQ.(0-MQ) ) THEN
                           CCCDN = CGC(IMKM1,1)*CGC(IKM2,1)
     &                             *CGC1(LB2,L1,(MUE2M05+1),MQ)
                           CCCUP = CGC(IMKM1,2)*CGC(IKM2,2)
     &                             *CGC1(L2,LB1,(MUE2M05+0),MQ)
                           B2_ADA(IKM1,IKM2,2,LAMQ) = FAC*(CCCUP-CCCDN)
C
                           B2_ADA(IKM1,IKM2,2,LAMQ) = FAC*CGC(IMKM1,1)
     &                        *CGC(IKM2,2)*CGC1(L2,LB1,(MUE1M05+1),MQ)
     &                        *SQRT_2
                        END IF
C
C                                                              (+)
                        IF ( (MUE1M05-MUE2M05).EQ.(+1-MQ) )
     &                       B2_ADA(IKM1,IKM2,3,LAMQ)
     &                       = -FAC*CGC(IMKM1,2)*CGC(IKM2,1)
     &                       *CGC1(L2,LB1,(MUE2M05+1),MQ)*SQRT_2
                     END IF
C
                  END DO
C ----------------------------------------------------------------------
C
               END DO
            END DO
         END DO
      END DO
C
      END
