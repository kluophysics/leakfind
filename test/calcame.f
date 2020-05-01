C*==calcame.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE CALCAME
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements connected with           *
C   *                                                                  *
C   *   1: IDOS  < LAM |      1      | LAM' >    (l,s)-resolved DOS    *
C   *   2: ISMT  < LAM | sigma(ipol) | LAM' >    spin moment           *
C   *   3: IOMT  < LAM |     l(ipol) | LAM' >    orbital moment        *
C   *   4: IHFF  < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
C   *   5: ISDM  < LAM |     T(ipol) | LAM' >    spin dipole moment    *
C   *   6: IKDS              (dummy)             kappa-resolved DOS    *
C   *   7: IBND  < LAM |      1      | LAM' >    band energy           *
C   *   8: IODN  < LAM |P_dn l(ipol) | LAM' >    orb. mnt. spin down   *
C   *   9: IOUP  < LAM |P_up l(ipol) | LAM' >    orb. mnt. spin up     *
C   *  10: ITRQ  < LAM |     t(ipol) | LAM' >    torque operator       *
C   *  11: IWFD                                  Weiss field           *
C   *                                                                  *
C   *   IPOL= 1,2,3  ==  -1,0,+1  ==  (-),(z),(+)                      *
C   *                                                                  *
C   *   IKM = 2 * l * (j+1/2) + j + mj + 1                             *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   * NOTE: the matrix elements include the KAPPA corresponding to     *
C   *       KAPPA(f) = -KAPPA(g) = -( -l_max-1 ) = l_max+1             *
C   *                                                                  *
C   *   list revised with version 7.4.0  Oct. 2015                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:PI,SQRT_2,SQRT_4PI
      USE MOD_ANGMOM,ONLY:AME_G,AME_F,CGC,NMEMAX,TXT_OBS,NK_EXT,
     &    NKMP_EXT,IMKM_IKM,NKM,IDOS,ISMT,IOMT,IHFF,IBND,ISDM,IODN,IOUP,
     &    IKDS
      IMPLICIT NONE
C*--CALCAME39
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CALCAME')
      REAL*8 TOL
      PARAMETER (TOL=1D-12)
      INTEGER NMETAB
      PARAMETER (NMETAB=9)
C
C Local variables
C
      REAL*8 AIJ,AJI,AME1(:,:,:),AME2(:,:,:),MJ1,MJ2,MS,PRE,PREHF,RSUM,
     &       TDIA,TERM1,TOFF
      LOGICAL CHECK
      CHARACTER*1 CHPOL(3)
      REAL*8 GAUNT_CYLM
      INTEGER I,IA_ERR,IERR,IKM1,IKM2,IME,IMKM1,IMKM2,INTAB(NMETAB),
     &        IOS(-1:0),IPOL,IS,J,J1P05,J2P05,K1,K2,KAP1,KAP2,KK,L1,L2,
     &        LB1,LB2,LR,M1,M2,MPOL,MSM05,MUE1M05,MUE2M05,NDIFF
      CHARACTER*20 STR20
C
C*** End of declarations rewritten by SPAG
C
      DATA CHPOL/'-','z','+'/
      DATA INTAB/IDOS,ISMT,IOMT,IHFF,ISDM,IKDS,IBND,IODN,IOUP/
      DATA IOS/IODN,IOUP/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE AME1,AME2
C
      IF ( ALLOCATED(AME_G) ) DEALLOCATE (AME_G,AME_F)
      ALLOCATE (AME_G(NKMP_EXT,NKMP_EXT,3,NMEMAX))
      ALLOCATE (AME_F(NKMP_EXT,NKMP_EXT,3,NMEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: AME_G')
C
      AME_G(:,:,:,:) = 0D0
      AME_F(:,:,:,:) = 0D0
C
      ALLOCATE (AME1(NKMP_EXT,NKMP_EXT,3),AME2(NKMP_EXT,NKMP_EXT,3),
     &          STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: AME1')
C
C ----------------------------------------------------------------------
C                   check pointers for matrix elements
C ----------------------------------------------------------------------
C
      IERR = 0
      LR = LEN_TRIM(ROUTINE)
      DO I = 1,NMETAB
         IF ( INTAB(I).GT.NMETAB ) THEN
            WRITE (6,99001) ROUTINE(1:LR),I,INTAB(I),NMETAB
            IERR = 1
         END IF
         DO J = 1,NMETAB
            IF ( I.NE.J .AND. INTAB(I).EQ.INTAB(J) ) THEN
               WRITE (6,99002) ROUTINE(1:LR),I,INTAB(I),J,INTAB(J)
               IERR = 1
            END IF
         END DO
      END DO
      IF ( IERR.EQ.1 ) CALL STOP_MESSAGE(ROUTINE,'IERR.EQ.1 ')
C
      IF ( NMEMAX.GT.NMETAB )
     &      CALL STOP_MESSAGE(ROUTINE,'NMEMAX > NMETAB')
C
C ----------------------------------------------------------------------
C
      CHECK = .TRUE.
C
      PREHF = SQRT(8D0*PI/3D0)
C
      IF ( IREL.LE.1 ) RETURN
C
C=======================================================================
C                                 DOS                               IDOS
C=======================================================================
C
      DO IPOL = 1,3
         DO IKM1 = 1,NKMP_EXT
            AME_G(IKM1,IKM1,IPOL,IDOS) = 1D0
         END DO
      END DO
C
C=======================================================================
C               G(LAMBDA,LAMBDA') = <LAMBDA|sig_lambda|LAMBDA'>     ISMT
C=======================================================================
C
      IF ( NMEMAX.GE.ISMT ) THEN
C
         AME1(:,:,:) = 0D0
         AME2(:,:,:) = 0D0
C
         IKM1 = 0
         DO K1 = 1,NK_EXT + 1
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            LB1 = L1 - SIGN(1,KAP1)
            J1P05 = IABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               MJ1 = MUE1M05 + 0.5D0
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK_EXT + 1
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  J2P05 = IABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     MJ2 = MUE2M05 + 0.5D0
                     IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                     IF ( L1.EQ.L2 ) THEN
C
                        IF ( (MUE1M05-MUE2M05).EQ.-1 ) THEN
C
                           AME_G(IKM1,IKM2,1,ISMT) = +SQRT_2*CGC(IKM1,1)
     &                        *CGC(IKM2,2)
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.0 ) THEN
C
                           AME_G(IKM1,IKM2,2,ISMT) = CGC(IKM1,2)
     &                        *CGC(IKM2,2) - CGC(IKM1,1)*CGC(IKM2,1)
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.+1 ) THEN
C
                           AME_G(IKM1,IKM2,3,ISMT) = -SQRT_2*CGC(IKM1,2)
     &                        *CGC(IKM2,1)
C
                        END IF
                     END IF
C ----------------------------------------------------------------------
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     IF ( CHECK ) THEN
                        RSUM = 0D0
                        DO IS = 1,2
                           MS = -1.5D0 + IS
                           M1 = NINT(MJ1-MS)
                           M2 = NINT(MJ2-MS)
                           RSUM = RSUM + 2*MS*CGC(IKM1,IS)*CGC(IKM2,IS)
     &                            *GAUNT_CYLM(L1,M1,0,0,L2,M2)
                        END DO
                        AME1(IKM1,IKM2,2) = SQRT_4PI*RSUM
C
                        IF ( ABS(MJ1-MJ2).LT.1D-6 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME2(IKM1,IKM2,2) = -MJ1/(KAP1+0.5D0)
                           ELSE IF ( KAP1.EQ.-KAP2-1 ) THEN
                              AME2(IKM1,IKM2,2)
     &                           = -2*DSQRT(0.25D0-DBLE(MJ1/(KAP1-KAP2))
     &                           **2)
                           END IF
                        END IF
                     END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  END DO
               END DO
            END DO
         END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
            NDIFF = 0
            DO IKM2 = 1,NKMP_EXT
               DO IKM1 = 1,NKMP_EXT
                  IF ( ABS(AME_G(IKM1,IKM2,1,ISMT)+AME_G(IKM2,IKM1,3,
     &                 ISMT)).GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,ISMT)-AME1(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,ISMT)-AME2(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME1(IKM1,IKM2,2)-AME2(IKM1,IKM2,2)).GT.TOL )
     &                 NDIFF = NDIFF + 1
               END DO
            END DO
C
            IF ( NDIFF.GT.0 ) THEN
               WRITE (6,99004) 'A  '//TXT_OBS(ISMT)//'  ('//CHPOL(2)
     &                         //')'
               CALL RMATSTRUCT('<LAM| sig_z |LAM''>',AME_G(1,1,2,ISMT),
     &                         NKM,NKMP_EXT,3,3,0,1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (z)',AME1(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
               CALL RMATSTRUCT('ALTER 2 (z)',AME2(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
            END IF
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C-----------------------------------------------------------------------
C      check and enforce symmetry of the angular matrix elements
C-----------------------------------------------------------------------
C
      IPOL = 2
C
      DO J = 1,NKMP_EXT
         DO I = J + 1,NKMP_EXT
            AIJ = AME_G(I,J,IPOL,ISMT)
            AJI = AME_G(J,I,IPOL,ISMT)
            IF ( ABS(AIJ-AJI).GT.1D-14 ) THEN
               WRITE (6,99005) ISMT,IPOL,I,J,AIJ,AJI
               CALL STOP_MESSAGE(ROUTINE,'AME not SYMMATRIC')
            END IF
            AIJ = (AIJ+AJI)/2
            AME_G(I,J,IPOL,ISMT) = AIJ
            AME_G(J,I,IPOL,ISMT) = AIJ
         END DO
      END DO
C
C
C
C
C=======================================================================
C               F(LAMBDA,LAMBDA') = <LAMBDA|l_lambda|LAMBDA'>       IOMT
C=======================================================================
C
      IF ( NMEMAX.GE.IOMT ) THEN
C
         AME1(:,:,:) = 0D0
         AME2(:,:,:) = 0D0
C
         IKM1 = 0
         DO K1 = 1,NK_EXT + 1
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            LB1 = L1 - SIGN(1,KAP1)
            J1P05 = IABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               MJ1 = MUE1M05 + 0.5D0
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK_EXT + 1
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  J2P05 = IABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     MJ2 = MUE2M05 + 0.5D0
                     IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                     IF ( L1.EQ.L2 ) THEN
                        IF ( (MUE1M05-MUE2M05).EQ.-1 ) THEN
C
                           RSUM = 0D0
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = RSUM + CGC(IKM1,MSM05+2)
     &                               *CGC(IKM2,MSM05+2)
     &                               *SQRT(DBLE((L2+M2)*(L2-M2+1)))
                           END DO
                           AME_G(IKM1,IKM2,1,IOMT) = RSUM/SQRT_2
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.0 ) THEN
C
                           RSUM = 0D0
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = RSUM + CGC(IKM1,MSM05+2)
     &                               *CGC(IKM2,MSM05+2)*M2
                           END DO
                           AME_G(IKM1,IKM2,2,IOMT) = RSUM
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.+1 ) THEN
C
                           RSUM = 0D0
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = RSUM + CGC(IKM1,MSM05+2)
     &                               *CGC(IKM2,MSM05+2)
     &                               *SQRT(DBLE((L2-M2)*(L2+M2+1)))
                           END DO
                           AME_G(IKM1,IKM2,3,IOMT) = -RSUM/SQRT_2
C
                        END IF
                     END IF
C ----------------------------------------------------------------------
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     IF ( CHECK ) THEN
                        RSUM = 0D0
                        DO IS = 1,2
                           MS = -1.5D0 + IS
                           M1 = NINT(MJ1-MS)
                           M2 = NINT(MJ2-MS)
                           RSUM = RSUM + (MJ1-MS)*CGC(IKM1,IS)
     &                            *CGC(IKM2,IS)
     &                            *GAUNT_CYLM(L1,M1,0,0,L2,M2)
                        END DO
                        AME1(IKM1,IKM2,2) = SQRT_4PI*RSUM
C
                        IF ( ABS(MJ1-MJ2).LT.1D-6 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME2(IKM1,IKM2,2) = MJ1*(KAP1+1.0D0)
     &                           /(KAP1+0.5D0)
                           ELSE IF ( KAP1.EQ.-KAP2-1 ) THEN
                              AME2(IKM1,IKM2,2)
     &                           = DSQRT(0.25D0-(MJ1/(KAP1-KAP2))**2)
                           END IF
                        END IF
                     END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  END DO
               END DO
            END DO
         END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
            NDIFF = 0
            DO IKM2 = 1,NKMP_EXT
               DO IKM1 = 1,NKMP_EXT
                  IF ( ABS(AME_G(IKM1,IKM2,1,IOMT)+AME_G(IKM2,IKM1,3,
     &                 IOMT)).GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,IOMT)-AME1(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,IOMT)-AME2(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME1(IKM1,IKM2,2)-AME2(IKM1,IKM2,2)).GT.TOL )
     &                 NDIFF = NDIFF + 1
               END DO
            END DO
C
            IF ( NDIFF.GT.0 ) THEN
               WRITE (6,99004) 'A  '//TXT_OBS(IOMT)//'  ('//CHPOL(2)
     &                         //')'
               CALL RMATSTRUCT('<LAM| l_z |LAM''>',AME_G(1,1,2,IOMT),
     &                         NKM,NKMP_EXT,3,3,0,1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (z)',AME1(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
               CALL RMATSTRUCT('ALTER 2 (z)',AME2(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
            END IF
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C
C=======================================================================
C          A(LAMBDA,LAMBDA') = <LAMBDA| (sig x r)_lambda |LAMBDA'>  IHFF
C=======================================================================
C
      IF ( NMEMAX.GE.IHFF ) THEN
C
         AME1(:,:,:) = 0D0
         AME2(:,:,:) = 0D0
C
         IKM1 = 0
         DO K1 = 1,NK_EXT + 1
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            LB1 = L1 - SIGN(1,KAP1)
            J1P05 = IABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               MJ1 = MUE1M05 + 0.5D0
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK_EXT + 1
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  LB2 = L2 - SIGN(1,KAP2)
                  J2P05 = IABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     MJ2 = MUE2M05 + 0.5D0
                     IKM2 = IKM2 + 1
                     IMKM2 = 2*LB2*J2P05 + J2P05 + MUE2M05 + 1
C ----------------------------------------------------------------------
                     MPOL = -1
                     IS = NINT((MPOL+1D0)/2D0) + 1
                     MS = -1.5D0 + IS
                     M1 = NINT(MJ1-MS)
                     M2 = NINT(MJ2+MS)
                     TERM1 = (1D0/SQRT_2)*CGC(IKM1,IS)*CGC(IMKM2,3-IS)
     &                       *GAUNT_CYLM(L1,M1,1,0,LB2,M2)
                     RSUM = 0D0
                     DO IS = 1,2
                        MS = -1.5D0 + IS
                        M1 = NINT(MJ1-MS)
                        M2 = NINT(MJ2-MS)
                        RSUM = RSUM + MS*CGC(IKM1,IS)*CGC(IMKM2,IS)
     &                         *GAUNT_CYLM(L1,M1,1,MPOL,LB2,M2)
                     END DO
                     AME_G(IKM1,IKM2,1,IHFF)
     &                  = SQRT_2*PREHF*(TERM1+MPOL*RSUM)
C
                     RSUM = 0D0
                     DO IS = 1,2
                        MS = -1.5D0 + IS
                        M1 = NINT(MJ1-MS)
                        M2 = NINT(MJ2+MS)
                        RSUM = RSUM + CGC(IKM1,IS)*CGC(IMKM2,3-IS)
     &                         *GAUNT_CYLM(L1,M1,1,NINT(-2*MS),LB2,M2)
                     END DO
                     AME_G(IKM1,IKM2,2,IHFF) = PREHF*RSUM
C
                     MPOL = +1
                     IS = NINT((MPOL+1D0)/2D0) + 1
                     MS = -1.5D0 + IS
                     M1 = NINT(MJ1-MS)
                     M2 = NINT(MJ2+MS)
                     TERM1 = (1D0/SQRT_2)*CGC(IKM1,IS)*CGC(IMKM2,3-IS)
     &                       *GAUNT_CYLM(L1,M1,1,0,LB2,M2)
                     RSUM = 0D0
                     DO IS = 1,2
                        MS = -1.5D0 + IS
                        M1 = NINT(MJ1-MS)
                        M2 = NINT(MJ2-MS)
                        RSUM = RSUM + MS*CGC(IKM1,IS)*CGC(IMKM2,IS)
     &                         *GAUNT_CYLM(L1,M1,1,MPOL,LB2,M2)
                     END DO
                     AME_G(IKM1,IKM2,3,IHFF)
     &                  = SQRT_2*PREHF*(TERM1+MPOL*RSUM)
C
C ----------------------------------------------------------------------
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     IF ( CHECK ) THEN
                        IF ( NINT(MJ1-MJ2).EQ.-1 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME1(IKM1,IKM2,1)
     &                           = ((4D0*KAP1)/(4D0*KAP1**2-1D0))
     &                           *SQRT(KAP1**2-(MJ1+0.5D0)**2)/SQRT_2
                           ELSE IF ( KAP1.EQ.-KAP2-1 ) THEN
                              KK = -1
                              AME1(IKM1,IKM2,1) = -(1D0/(2D0*KAP1-KK))
     &                           *SQRT((KAP1-KK*(MJ1+0.5D0))
     &                           *(KAP1-KK*(MJ1+1.5D0)))/SQRT_2
                           ELSE IF ( KAP1.EQ.-KAP2+1 ) THEN
                              KK = +1
                              AME1(IKM1,IKM2,1) = -(1D0/(2D0*KAP1-KK))
     &                           *SQRT((KAP1-KK*(MJ1+0.5D0))
     &                           *(KAP1-KK*(MJ1+1.5D0)))/SQRT_2
                           END IF
                        END IF
C
                        IF ( ABS(MJ1-MJ2).LT.1D-6 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME1(IKM1,IKM2,2)
     &                           = 4D0*MJ1*KAP1/(4D0*KAP1**2-1D0)
                           ELSE IF ( KAP1.EQ.-KAP2-1 ) THEN
                              AME1(IKM1,IKM2,2)
     &                           = DSQRT(0.25D0-(MJ1/DBLE(KAP1-KAP2))
     &                           **2)
                           ELSE IF ( KAP1.EQ.-KAP2+1 ) THEN
                              AME1(IKM1,IKM2,2)
     &                           = -DSQRT(0.25D0-(MJ1/DBLE(KAP1-KAP2))
     &                           **2)
                           END IF
                        END IF
C
                        IF ( NINT(MJ1-MJ2).EQ.+1 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME1(IKM1,IKM2,3)
     &                           = -((4D0*KAP2)/(4D0*KAP2**2-1D0))
     &                           *SQRT(KAP2**2-(MJ2+0.5D0)**2)/SQRT_2
                           ELSE IF ( KAP1.EQ.-KAP2-1 ) THEN
                              KK = -1
                              AME1(IKM1,IKM2,3) = +(1D0/(2D0*KAP2-KK))
     &                           *SQRT((KAP2-KK*(MJ2+0.5D0))
     &                           *(KAP2-KK*(MJ2+1.5D0)))/SQRT_2
                           ELSE IF ( KAP1.EQ.-KAP2+1 ) THEN
                              KK = +1
                              AME1(IKM1,IKM2,3) = (1D0/(2D0*KAP2-KK))
     &                           *SQRT((KAP2-KK*(MJ2+0.5D0))
     &                           *(KAP2-KK*(MJ2+1.5D0)))/SQRT_2
                           END IF
                        END IF
C
                     END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  END DO
               END DO
            END DO
         END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
            NDIFF = 0
            DO IKM2 = 1,NKMP_EXT
               DO IKM1 = 1,NKMP_EXT
                  IF ( ABS(AME_G(IKM1,IKM2,1,IHFF)+AME_G(IKM2,IKM1,3,
     &                 IHFF)).GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,1,IHFF)-AME1(IKM1,IKM2,1))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,IHFF)-AME1(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,3,IHFF)-AME1(IKM1,IKM2,3))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
               END DO
            END DO
C
            IF ( NDIFF.GT.0 ) THEN
               WRITE (6,99004) 'A  '//TXT_OBS(IHFF)//'  ('//CHPOL(1)
     &                         //')'
               CALL RMATSTRUCT('<LAM|(sig x r)_(-)|LAM''>',
     &                         AME_G(1,1,1,IHFF),NKM,NKMP_EXT,3,3,0,
     &                         1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (-)',AME1(1,1,1),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               WRITE (6,99004) 'A  '//TXT_OBS(IHFF)//'  ('//CHPOL(2)
     &                         //')'
               CALL RMATSTRUCT('<LAM|(sig x r)_z|LAM''>',
     &                         AME_G(1,1,2,IHFF),NKM,NKMP_EXT,3,3,0,
     &                         1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (z)',AME1(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               WRITE (6,99004) 'A  '//TXT_OBS(IHFF)//'  ('//CHPOL(3)
     &                         //')'
               CALL RMATSTRUCT('<LAM|(sig x r)_(+)|LAM''>',
     &                         AME_G(1,1,3,IHFF),NKM,NKMP_EXT,3,3,0,
     &                         1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (+)',AME1(1,1,3),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
            END IF
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C=======================================================================
C               T(LAMBDA,LAMBDA') = <LAMBDA|T_lambda|LAMBDA'>       ISDM
C=======================================================================
C
      IF ( NMEMAX.GE.ISDM ) THEN
C
         AME1(:,:,:) = 0D0
         AME2(:,:,:) = 0D0
C
         PRE = -SQRT(4D0*PI/3D0)
         IKM1 = 0
         DO K1 = 1,NK_EXT + 1
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            LB1 = L1 - SIGN(1,KAP1)
            J1P05 = IABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               MJ1 = MUE1M05 + 0.5D0
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK_EXT + 1
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  LB2 = L2 - SIGN(1,KAP2)
                  J2P05 = IABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     MJ2 = MUE2M05 + 0.5D0
                     IKM2 = IKM2 + 1
                     IMKM2 = 2*LB2*J2P05 + J2P05 + MUE2M05 + 1
C ----------------------------------------------------------------------
C
                     MPOL = MUE1M05 - MUE2M05
C
                     IF ( ABS(MPOL).LE.1 ) THEN
C
                        IPOL = MPOL + 2
                        RSUM = 0D0
                        DO IS = 1,2
                           MS = -1.5D0 + IS
                           M1 = NINT(MJ1-MS)
                           M2 = NINT(MJ2-MS)
                           RSUM = RSUM + CGC(IKM1,IS)*CGC(IMKM2,IS)
     &                            *GAUNT_CYLM(L1,M1,1,MPOL,LB2,M2)
                        END DO
                        AME_G(IKM1,IKM2,IPOL,ISDM)
     &                     = AME_G(IKM1,IKM2,IPOL,ISMT) - 3*PRE*RSUM
C
                     END IF
C ----------------------------------------------------------------------
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     IF ( CHECK ) THEN
                        IF ( ABS(MJ1-MJ2).LT.1D-6 ) THEN
                           IF ( KAP1.EQ.KAP2 ) THEN
                              AME1(IKM1,IKM2,2) = -4*MJ1*(KAP1+1D0)
     &                           /(4D0*KAP1**2-1D0)
C
                              TDIA = 2*MJ1/DBLE((2*L1+1)*(2*LB1+1))
                              AME2(IKM1,IKM2,2)
     &                           = AME_G(IKM1,IKM2,2,ISMT) - 3.0D0*TDIA
                           ELSE IF ( KAP1.EQ.(-KAP2-1) ) THEN
C
                              AME1(IKM1,IKM2,2)
     &                           = SQRT(0.25D0-(MJ1/(KAP1-KAP2))**2)
C
                              TOFF = -SQRT((L1+0.5D0)**2-MJ1**2)
     &                               /DBLE(2*L1+1)
                              AME2(IKM1,IKM2,2)
     &                           = AME_G(IKM1,IKM2,2,ISMT) - 3.0D0*TOFF
                           ELSE IF ( KAP1.EQ.(-KAP2+1) ) THEN
C
                              AME1(IKM1,IKM2,2)
     &                           = 3*SQRT(0.25D0-(MJ1/(KAP1-KAP2))**2)
C
                              AME2(IKM1,IKM2,2) = AME1(IKM1,IKM2,2)
                           END IF
                        END IF
                     END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  END DO
               END DO
            END DO
         END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK .AND. NKMP_EXT.LT.0 ) THEN
            NDIFF = 0
            DO IKM2 = 1,NKMP_EXT
               DO IKM1 = 1,NKMP_EXT
                  IF ( ABS(AME_G(IKM1,IKM2,1,ISDM)+AME_G(IKM2,IKM1,3,
     &                 ISDM)).GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,ISDM)-AME1(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME_G(IKM1,IKM2,2,ISDM)-AME2(IKM1,IKM2,2))
     &                 .GT.TOL ) NDIFF = NDIFF + 1
                  IF ( ABS(AME1(IKM1,IKM2,2)-AME2(IKM1,IKM2,2)).GT.TOL )
     &                 NDIFF = NDIFF + 1
               END DO
            END DO
C
            IF ( NDIFF.GT.0 ) THEN
               WRITE (6,99004) 'A  '//TXT_OBS(ISDM)//'  ('//CHPOL(1)
     &                         //')'
               CALL RMATSTRUCT('<LAM| T_(-) |LAM''>',AME_G(1,1,1,ISDM),
     &                         NKM,NKMP_EXT,3,3,0,1D-8,6)
C           CALL  RMATSTR('ALTER 1 (-)',AME1(1,1,1),
C    &           NKM,NKMP_EXT,3,3,0,1D-8,6)
C
               WRITE (6,99004) 'A  '//TXT_OBS(ISDM)//'  ('//CHPOL(2)
     &                         //')'
               CALL RMATSTRUCT('<LAM| T_(z) |LAM''>',AME_G(1,1,2,ISDM),
     &                         NKM,NKMP_EXT,3,3,0,1D-8,6)
               CALL RMATSTRUCT('ALTER 1 (z)',AME1(1,1,2),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
C
               WRITE (6,99004) 'A  '//TXT_OBS(ISDM)//'  ('//CHPOL(3)
     &                         //')'
               CALL RMATSTRUCT('<LAM| T_(+) |LAM''>',AME_G(1,1,3,ISDM),
     &                         NKM,NKMP_EXT,3,3,0,1D-8,6)
C           CALL  RMATSTR('ALTER 1 (+)',AME1(1,1,3),
C    &          NKM,NKMP_EXT,3,3,0,1D-8,6)
C
               CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
            END IF
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C=======================================================================
C                       kappa-resolved  DOS                         IKDS
C=======================================================================
C
      IF ( NMEMAX.GE.IKDS ) THEN
C
         DO IPOL = 1,3
            DO IKM1 = 1,NKMP_EXT
               AME_G(IKM1,IKM1,IPOL,IKDS) = 1D0
            END DO
         END DO
C
      END IF
C
C=======================================================================
C                                 BND                               IBND
C=======================================================================
C
      IF ( NMEMAX.GE.IBND ) THEN
C
         AME_G(:,:,:,IBND) = AME_G(:,:,:,IDOS)
         AME_F(:,:,:,IBND) = AME_F(:,:,:,IDOS)
C
      END IF
C
C=======================================================================
C         F(LAMBDA,LAMBDA',SPIN) = <LAMBDA|P_spin l_lambda|LAMBDA'>  7,8
C=======================================================================
C
      IF ( NMEMAX.GE.IODN .AND. NMEMAX.GE.IOUP ) THEN
C
         IKM1 = 0
         DO K1 = 1,NK_EXT + 1
            L1 = K1/2
            IF ( MOD(K1,2).EQ.0 ) THEN
               KAP1 = L1
            ELSE
               KAP1 = -L1 - 1
            END IF
            LB1 = L1 - SIGN(1,KAP1)
            J1P05 = IABS(KAP1)
C
            DO MUE1M05 = -J1P05,J1P05 - 1
               MJ1 = MUE1M05 + 0.5D0
               IKM1 = IKM1 + 1
C
               IKM2 = 0
               DO K2 = 1,NK_EXT + 1
                  L2 = K2/2
                  IF ( MOD(K2,2).EQ.0 ) THEN
                     KAP2 = L2
                  ELSE
                     KAP2 = -L2 - 1
                  END IF
                  J2P05 = IABS(KAP2)
C
                  DO MUE2M05 = -J2P05,J2P05 - 1
                     MJ2 = MUE2M05 + 0.5D0
                     IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                     IF ( L1.EQ.L2 ) THEN
                        IF ( (MUE1M05-MUE2M05).EQ.-1 ) THEN
C
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                               *SQRT(DBLE((L2+M2)*(L2-M2+1)))
                              AME_G(IKM1,IKM2,1,IOS(MSM05))
     &                           = RSUM/SQRT_2
                           END DO
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.0 ) THEN
C
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                               *M2
                              AME_G(IKM1,IKM2,2,IOS(MSM05)) = RSUM
                           END DO
C
                        ELSE IF ( (MUE1M05-MUE2M05).EQ.+1 ) THEN
C
                           DO MSM05 = -1,0
                              M2 = MUE2M05 - MSM05
                              RSUM = CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                               *SQRT(DBLE((L2-M2)*(L2+M2+1)))
                              AME_G(IKM1,IKM2,3,IOS(MSM05))
     &                           = -RSUM/SQRT_2
                           END DO
C
                        END IF
                     END IF
C ----------------------------------------------------------------------
                  END DO
               END DO
            END DO
         END DO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( CHECK ) THEN
            NDIFF = 0
            DO IKM2 = 1,NKMP_EXT
               DO IKM1 = 1,NKMP_EXT
                  DO IPOL = 1,3
                     RSUM = 0D0
                     DO MSM05 = -1,0
                        RSUM = RSUM + AME_G(IKM1,IKM2,IPOL,IOS(MSM05))
                     END DO
C
                     IF ( ABS(AME_G(IKM1,IKM2,IPOL,IOMT)-RSUM).GT.TOL )
     &                    NDIFF = NDIFF + 1
                  END DO
               END DO
            END DO
C
            IF ( NDIFF.GT.0 ) THEN
               DO IPOL = 1,3
                  WRITE (6,99004) 'A  '//TXT_OBS(IOMT)
     &                            //'  ('//CHPOL(IPOL)//')'
                  CALL RMATSTRUCT('<LAM| l_z |LAM''>',
     &                            AME_G(1,1,IPOL,IOMT),NKM,NKMP_EXT,3,3,
     &                            0,1D-8,6)
                  CALL RMATSTRUCT('<LAM| l_z dn|LAM''>',
     &                            AME_G(1,1,IPOL,IODN),NKM,NKMP_EXT,3,3,
     &                            0,1D-8,6)
                  CALL RMATSTRUCT('<LAM| l_z up|LAM''>',
     &                            AME_G(1,1,IPOL,IOUP),NKM,NKMP_EXT,3,3,
     &                            0,1D-8,6)
               END DO
               CALL STOP_MESSAGE(ROUTINE,'inconsistencies occured')
            END IF
         END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C
C
C
C=======================================================================
      IF ( IPRINT.GE.1 ) THEN
C
         WRITE (6,99003) ROUTINE(1:LEN_TRIM(ROUTINE))
         DO IME = 1,5
            DO IPOL = 1,3
               STR20 = 'A  '//TXT_OBS(IME)//'  ('//CHPOL(IPOL)//')'
               CALL RMATSTRUCT(STR20,AME_G(1,1,IPOL,IME),NKM,NKMP_EXT,3,
     &                         3,0,1D-8,6)
            END DO
         END DO
C
      END IF
C=======================================================================
C
      DO IKM1 = 1,NKMP_EXT
         DO IKM2 = 1,NKMP_EXT
            IMKM1 = IMKM_IKM(IKM1)
            IMKM2 = IMKM_IKM(IKM2)
            IF ( IMKM1.LE.NKMP_EXT .AND. IMKM2.LE.NKMP_EXT ) THEN
C
               AME_F(IKM1,IKM2,1:3,1:NMEMAX)
     &            = AME_G(IMKM1,IMKM2,1:3,1:NMEMAX)
C
            ELSE
C
               AME_F(IKM1,IKM2,1:3,1:NMEMAX) = 0D0
C
            END IF
         END DO
      END DO
C
      DEALLOCATE (AME1,AME2)
C
99001 FORMAT (/,10X,'ERROR in <',A,'>',I3,'th  index =',I2,
     &        ' > NMEMAX =',I2)
99002 FORMAT (/,10X,'ERROR in <',A,'>',I3,'th  index =',I2,' = ',I3,
     &        'th  index =',I2)
99003 FORMAT (//,1X,79('*'),/,36X,'<',A,'>',/,1X,79('*'),//,10X,
     &        'calculation of angular matrix elements  AME ',/)
99004 FORMAT (//,10X,'PROBLEMS with angular matrix elements ',A,//)
99005 FORMAT (/,10X,'setting up angular matrix elements AME_RLM',/,10X,
     &        'for IMVEC =',I3,'  IPOL =',I3,'  I =',I3,'  J =',I3,/,
     &        10X,'AIJ = ',E25.15,/,:,10X,'AJI = ',E25.15,/)
      END
