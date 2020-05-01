C*==ame_nab_sig.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AME_NAB_SIG
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements needed  for the          *
C   *                                                                  *
C   *   Bargmann-Wigner spin-polarization operator                     *
C   *                                                                  *
C   *                ->   ->                                           *
C   *   A1 = < LAM| NAB . e(IPOLSPIN) SIGMA(JPOL) | LAM'>              *
C   *   A2 = < LAM| NAB . e(IPOLSPIN) SIGMA(JPOL) | LAM'>              *
C   *                                                                  *
C   *   A1, A2 are associated with the two terms of the gradient       *
C   *          formula                                                 *
C   *                                                                  *
C   *   ipolspin, jpol = 1,2,3  == (-),(0),(+)                         *
C   *                                                                  *
C   *   NB!: The following convention is used for transforming from    *
C   *        cartesian to (special) spherical coordinates              *
C   *                                                                  *
C   *  ->        ->       ->                                           *
C   *  e(+) = - [e(x) + i e(y)] / sqrt(2) left circ. positive helicity *
C   *                                                                  *
C   *  ->        ->       ->                                           *
C   *  e(-) = + [e(x) - i e(y)] / sqrt(2) right circ. negative helicity*
C   *                                                                  *
C   *  DK, 09/2013,01/2015                                             *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:IMKM_IKM,CGC,NK_EXT,NKM,NKMP_EXT,NKM_EXT,
     &    AF1_SPIN_CURR_NAB,AF2_SPIN_CURR_NAB,AG1_SPIN_CURR_NAB,
     &    AG2_SPIN_CURR_NAB
      USE MOD_CONSTANTS,ONLY:SQRT_2
      IMPLICIT NONE
C*--AME_NAB_SIG35
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 AF1(:,:,:,:),AF2(:,:,:,:),AG1(:,:,:,:),AG2(:,:,:,:),CLCL,
     &       CLCL1,WZ2
      REAL*8 CGC1,SQM,SQP
      INTEGER I1,I2,IDXPOL(-1:+1),IKM1,IKM2,IMKM1,IMKM2,IP,IPOL,J1P05,
     &        J2P05,JP,JPOL,K1,K2,KAP1,KAP2,L1,L2,LB1,LB2,M,MM,MUE1M05,
     &        MUE2M05,N
C
C*** End of declarations rewritten by SPAG
C
      DATA IDXPOL/1,2,3/
C
C-----------------------------------------------------------------------
C
      ALLOCATABLE AG1,AG2,AF1,AF2
C
C=======================================================================
C angular matrix elements for NABLA-related term in spin current density
C=======================================================================
C
      IF ( .NOT.ALLOCATED(AG1_SPIN_CURR_NAB) ) THEN
         N = NKMP_EXT
         ALLOCATE (AG1_SPIN_CURR_NAB(N,N,3,3))
         ALLOCATE (AG2_SPIN_CURR_NAB(N,N,3,3))
         ALLOCATE (AF1_SPIN_CURR_NAB(N,N,3,3))
         ALLOCATE (AF2_SPIN_CURR_NAB(N,N,3,3))
      END IF
C
C=======================================================================
C                          work space
C=======================================================================
C
      N = NKMP_EXT
C
      ALLOCATE (AG1(1:N,1:N,-1:1,-1:1))
      ALLOCATE (AG2(1:N,1:N,-1:1,-1:1))
      ALLOCATE (AF1(1:N,1:N,-1:1,-1:1))
      ALLOCATE (AF2(1:N,1:N,-1:1,-1:1))
C
      AG1(1:N,1:N,-1:1,-1:1) = 0D0
      AG2(1:N,1:N,-1:1,-1:1) = 0D0
      AF1(1:N,1:N,-1:1,-1:1) = 0D0
      AF2(1:N,1:N,-1:1,-1:1) = 0D0
C
      WZ2 = SQRT_2
C
      IKM1 = 0
      LOOP_K1:DO K1 = 1,NK_EXT
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
         LOOP_MUE1M05:DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
            IMKM1 = LB1*2*J1P05 + J1P05 + MUE1M05 + 1
C
            IKM2 = 0
            LOOP_K2:DO K2 = 1,NK_EXT
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
               LOOP_MUE2M05:DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
                  IMKM2 = LB2*2*J2P05 + J2P05 + MUE2M05 + 1
C
                  MM = MUE1M05 - MUE2M05
C-----------------------------------------------------------------------
C
C---------------------------------------NB: IPOL denotes spin-pol index
C---------------------------------------    JPOL denotes alpha component
C
                  LOOP_IPOL:DO IPOL = -1,1
                     LOOP_JPOL:DO JPOL = -1,1
C
                        IF ( JPOL.EQ.0 ) THEN
                           CLCL = CGC(IKM1,2)*CGC(IKM2,2)
                           CLCL1 = CGC(IKM1,1)*CGC(IKM2,1)
                        ELSE
                           I1 = (3+JPOL)/2
                           I2 = (3-JPOL)/2
                           CLCL = CGC(IKM1,I1)*CGC(IKM2,I2)*WZ2
                        END IF
C
                        IF ( MM.EQ.IPOL+JPOL ) THEN
C
                           IF ( L1.EQ.L2+1 ) THEN
C
                              IF ( JPOL.EQ.0 ) THEN
C
                                 AG1(IKM1,IKM2,JPOL,IPOL) = SQP(L2)
     &                              *CGC1(L2,L2+1,MUE2M05,IPOL)
     &                              *CLCL - SQP(L2)
     &                              *CGC1(L2,L2+1,(MUE2M05+1),IPOL)
     &                              *CLCL1
                              ELSE
                                 I1 = MUE2M05
                                 IF ( JPOL.EQ.1 ) I1 = I1 + 1
                                 AG1(IKM1,IKM2,JPOL,IPOL) = SQP(L2)
     &                              *CGC1(L2,L2+1,I1,IPOL)*CLCL*(-JPOL)
                              END IF
                           END IF
C
                           IF ( L1.EQ.L2-1 ) THEN
C
                              IF ( JPOL.EQ.0 ) THEN
C
                                 AG2(IKM1,IKM2,JPOL,IPOL) = SQM(L2)
     &                              *CGC1(L2,L2-1,MUE2M05,IPOL)
     &                              *CLCL - SQM(L2)
     &                              *CGC1(L2,L2-1,(MUE2M05+1),IPOL)
     &                              *CLCL1
                              ELSE
                                 I1 = MUE2M05
                                 IF ( JPOL.EQ.1 ) I1 = I1 + 1
                                 AG2(IKM1,IKM2,JPOL,IPOL) = SQM(L2)
     &                              *CGC1(L2,L2-1,I1,IPOL)*CLCL*(-JPOL)
                              END IF
                           END IF
C
                        END IF
C
                     END DO LOOP_JPOL
                  END DO LOOP_IPOL
C
               END DO LOOP_MUE2M05
            END DO LOOP_K2
         END DO LOOP_MUE1M05
      END DO LOOP_K1
C
      DO IKM1 = 1,NKMP_EXT
         IMKM1 = IMKM_IKM(IKM1)
C
         DO IKM2 = 1,NKMP_EXT
            IMKM2 = IMKM_IKM(IKM2)
C
            DO IPOL = -1,1
               DO JPOL = -1,1
C
                  IF ( (IMKM1.LE.NKM) .AND. (IMKM2.LE.NKM) ) THEN
C
                     AF1(IKM1,IKM2,IPOL,JPOL)
     &                  = AG1(IMKM1,IMKM2,IPOL,JPOL)
                     AF2(IKM1,IKM2,IPOL,JPOL)
     &                  = AG2(IMKM1,IMKM2,IPOL,JPOL)
C
                  END IF
C
               END DO
            END DO
         END DO
      END DO
C
C
C---------------------------- translate from (-),(0),(+)  to (1),(2),(3)
      DO IP = -1,1
         DO JP = -1,1
            AG1_SPIN_CURR_NAB(:,:,IDXPOL(IP),IDXPOL(JP))
     &         = AG1(:,:,IP,JP)
            AG2_SPIN_CURR_NAB(:,:,IDXPOL(IP),IDXPOL(JP))
     &         = AG2(:,:,IP,JP)
            AF1_SPIN_CURR_NAB(:,:,IDXPOL(IP),IDXPOL(JP))
     &         = AF1(:,:,IP,JP)
            AF2_SPIN_CURR_NAB(:,:,IDXPOL(IP),IDXPOL(JP))
     &         = AF2(:,:,IP,JP)
         END DO
      END DO
C
C---------------------------------------------------- AG1 AG2 comparison
      IF ( 1.EQ.0 ) THEN
         N = NKM_EXT
         M = NKMP_EXT
         WRITE (6,*) 'new'
         WRITE (6,'(5E22.14)') AG1_SPIN_CURR_NAB
         WRITE (6,'(5E22.14)') AG2_SPIN_CURR_NAB
C
         CALL RMATSTRUCT('AG1 11',AG1_SPIN_CURR_NAB(:,:,1,1),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 21',AG2_SPIN_CURR_NAB(:,:,2,1),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 21',AG2_SPIN_CURR_NAB(:,:,2,1),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 11',AG1_SPIN_CURR_NAB(:,:,1,1),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 31',AG1_SPIN_CURR_NAB(:,:,3,1),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 31',AG2_SPIN_CURR_NAB(:,:,3,1),N,M,3,3,0,
     &                   1D-8,6)
C
         CALL RMATSTRUCT('AG1 12',AG1_SPIN_CURR_NAB(:,:,1,2),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 22',AG2_SPIN_CURR_NAB(:,:,2,2),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 22',AG2_SPIN_CURR_NAB(:,:,2,2),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 12',AG1_SPIN_CURR_NAB(:,:,1,2),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 32',AG1_SPIN_CURR_NAB(:,:,3,2),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 32',AG2_SPIN_CURR_NAB(:,:,3,2),N,M,3,3,0,
     &                   1D-8,6)
C
         CALL RMATSTRUCT('AG1 13',AG1_SPIN_CURR_NAB(:,:,1,3),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 23',AG2_SPIN_CURR_NAB(:,:,2,3),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 23',AG1_SPIN_CURR_NAB(:,:,2,3),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 13',AG2_SPIN_CURR_NAB(:,:,1,3),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG1 33',AG1_SPIN_CURR_NAB(:,:,3,3),N,M,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AG2 33',AG2_SPIN_CURR_NAB(:,:,3,3),N,M,3,3,0,
     &                   1D-8,6)
C
         STOP
      END IF
C
      END
C*==sqp.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      DOUBLE PRECISION FUNCTION SQP(INTARG)
      IMPLICIT NONE
C*--SQP282
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INTARG
C
C*** End of declarations rewritten by SPAG
C
C
      SQP = DSQRT(DBLE(INTARG+1)/DBLE(2*INTARG+3))
C
      END
C*==sqm.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      DOUBLE PRECISION FUNCTION SQM(INTARG)
      IMPLICIT NONE
C*--SQM307
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INTARG
C
C*** End of declarations rewritten by SPAG
C
      SQM = DSQRT(DBLE(INTARG)/DBLE(2*INTARG-1))
C
      END
C
