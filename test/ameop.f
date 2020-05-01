C*==ameop.f    processed by SPAG 6.70Rc at 15:34 on 19 Dec 2016
      SUBROUTINE AMEOP
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements occuring with the        *
C   *                      ORBITAL POLARISATION                        *
C   *                                                                  *
C   *    i.e. the spin-resolved   G- AND F-COEFFICIENTS                *
C   *                                                                  *
C   *   AMEOPC(LAM,LAM',MS) = CGC(K,MUE,MS)*CGC(K',MUE,MS)             *
C   *   AMEOPO(LAM,LAM',MS) = CGC(K,MUE,MS)*CGC(K',MUE,MS) * (MUE-MS)  *
C   *                                                                  *
C   * AMEOP*  allocated with dimension  NKMPMAX                        *
C   *         with AMEOP*(IKM) = 0 for IKM > NKM                       *
C   *         this allows NKM < IKM <= NKMPMAX for the small component *
C   *                                                                  *
C   *  27/11/07  HE                                                    *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:AMEOPC,AMEOPO,CGC,NKM,NLMAX,LTAB,KAPTAB,
     &    NMUETAB,NKMPMAX
      USE MOD_FILES,ONLY:IPRINT
      IMPLICIT NONE
C*--AMEOP23
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      INTEGER IA_ERR,IKM1,IKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,L2,M1,M2,
     &        MS,MUE1M05,MUE2M05,NK
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE (AMEOPC(NKMPMAX,NKMPMAX,2))
      ALLOCATE (AMEOPO(NKMPMAX,NKMPMAX,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: AMEOP -> AMEOPO'
C
      CALL RINIT(NKMPMAX**2*2,AMEOPC)
      CALL RINIT(NKMPMAX**2*2,AMEOPO)
C
C     JP05 = J +0.5     MUEP05 = MUE + 0.5
C     IKM  = L*2*(J+1/2) + J + MUE + 1
C
      NK = 2*NLMAX - 1
C
      IKM2 = 0
      DO K2 = 1,NK
         L2 = LTAB(K2)
         KAP2 = KAPTAB(K2)
         J2P05 = ABS(KAP2)
         MUE2M05 = -J2P05 - 1
C
         DO M2 = 1,NMUETAB(K2)
            MUE2M05 = MUE2M05 + 1
            IKM2 = IKM2 + 1
C
            IKM1 = 0
            DO K1 = 1,NK
               L1 = LTAB(K1)
               KAP1 = KAPTAB(K1)
               J1P05 = ABS(KAP1)
               MUE1M05 = -J1P05 - 1
C
               DO M1 = 1,NMUETAB(K1)
                  MUE1M05 = MUE1M05 + 1
                  IKM1 = IKM1 + 1
C
                  IF ( L1.EQ.L2 ) THEN
                     IF ( MUE1M05.EQ.MUE2M05 ) THEN
C
                        DO MS = 1,2
                           AMEOPC(IKM1,IKM2,MS) = CGC(IKM1,MS)
     &                        *CGC(IKM2,MS)
                           AMEOPO(IKM1,IKM2,MS) = CGC(IKM1,MS)
     &                        *CGC(IKM2,MS)*DBLE(MUE1M05-MS+2)
                        END DO
C
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO
C
      IF ( IPRINT.GT.1 ) THEN
         CALL RMATSTRUCT('AMEOPC  ms=1',AMEOPC(1,1,1),NKM,NKMPMAX,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AMEOPC  ms=2',AMEOPC(1,1,2),NKM,NKMPMAX,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AMEOPO  ms=1',AMEOPO(1,1,1),NKM,NKMPMAX,3,3,0,
     &                   1D-8,6)
         CALL RMATSTRUCT('AMEOPO  ms=2',AMEOPO(1,1,2),NKM,NKMPMAX,3,3,0,
     &                   1D-8,6)
      END IF
C
      END
