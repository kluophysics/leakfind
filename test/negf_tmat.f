C*==negf_tmat.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NEGF_TMAT(IFIL_ZPLS,IFIL_ZMIN,IQBOT_LBAR,IQTOP_RBAR,
     &                     MEZZ,MEZZMM,MEZZMP,MEZZPM,STT,DNREL2)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the matrix elements                                   *
C   *                                                                  *
C   *   M(z,z') = < R(z) | R(z') >                                     *
C   *                                                                  *
C   *  for z unequal z':                                               *
C   *                                                                  *
C   *   MEZZMP = < R(-) | A | R(+) >  for A = 1, Sigma_z               *
C   *   MEZZPM = < R(+) | A | R(-) >  for A = 1, Sigma_z               *
C   *                                                                  *
C   *  coming from above (runssite for zpls and zmin):                 *
C   *                                                                  *
C   *   MEZZ   = < R(+) | A | R(+) >  for A = 1, Sigma_z               *
C   *   MEZZMM = < R(-) | A | R(-) >  for A = 1, Sigma_z               *
C   *                                                                  *
C   *  rotate matrix elements from local to global FOR if requested    *
C   *                                                                  *
C   *   KMROT = 0 (no rot.), 1 (individual site rot.), 2 (global rot.) *
C   *                                                                  *
C   *  calculate the matrix element of the exchange field along        *
C   *  the magnetisation axis in layer I expanded into spherical       *
C   *  harmonics, using the local FOR                                  *
C   *                                                                  *
C   *   DELTA_I = < R(+,r) | beta sigma_z B(r) | R(+,r) >              *
C   *             ---------------------------------------              *
C   *                       < R(+,r) | R(+,r) >                        *
C   *                                                                  *
C   *  NOTE: the routine assumes that the RH-convention is used        *
C   *  for the Green's function                                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:AME_G,IKM1LIN,IKM2LIN,NKM,NKMMAX,NKMQ,NLM,
     &    NLMMAX,NLQ,NMEMAX,WKM1,NCPLWF
      USE MOD_CALCMODE,ONLY:GF_CONV_RH,IREL,KMROT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_RMESH,ONLY:JRWS,NRMAX,R2DRDI
      USE MOD_SITES,ONLY:DROTQ,ITOQ,NOQ
      USE MOD_TYPES,ONLY:BT,IMT,NT,NTMAX,NCPLWFMAX,IKMCPLWF
      IMPLICIT NONE
C*--NEGF_TMAT45
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NEGF_TMAT')
      INTEGER ISMT,CHECK_INTEGRATION
      PARAMETER (ISMT=2,CHECK_INTEGRATION=0)
      LOGICAL CHECK_ME
      PARAMETER (CHECK_ME=.FALSE.)
C
C Dummy arguments
C
      INTEGER IFIL_ZMIN,IFIL_ZPLS,IQBOT_LBAR,IQTOP_RBAR
      LOGICAL STT
      COMPLEX*16 DNREL2(NLMMAX,NLMMAX,NTMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMM(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZMP(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZPM(NKMMAX,NKMMAX,NTMAX,NMEMAX)
C
C Local variables
C
      COMPLEX*16 BZZ(:,:),DELTAIT2(:,:),DELTAIT3(:,:),DZZMP(:,:),
     &           DZZPM(:,:),JF(:,:,:),JFPLS(:,:,:),JG(:,:,:),
     &           JGPLS(:,:,:),NORMZG1ZG2,SZZMP(:,:),SZZPM(:,:),ZBBZ(:,:)
     &           ,ZBSBZ(:,:),ZBZ(:,:),ZF(:,:,:),ZFPLS(:,:,:),
     &           ZFPLSZF(2,2),ZFPLSZFPLS(2,2),ZFZFPLS(2,2),ZG(:,:,:),
     &           ZGPLS(:,:,:),ZGPLSZG(2,2),ZGPLSZGPLS(2,2),ZGZGPLS(2,2),
     &           ZSBZ(:,:)
      REAL*8 CGF(2,2),CGG(2,2),DB_R2DRDI(:),MJ
      INTEGER G1,G2,I,I1,I2,IA_ERR,IFIL,IKM,IKM1,IKM2,IKM_SOL(2),IM,IME,
     &        IMKMI,IMKM_SOL(2),IO,IQ,IRTOP,IT,IZ,J,JKM,K1,K2,KAP1,KAP2,
     &        KK,KKA,KKB,L,L1,L2,LAM1,LAM2,LIN,M,MJM05,N,NSOL
      INTEGER IKAPMUE
      LOGICAL KDONE_T(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE BZZ,DB_R2DRDI,DELTAIT2,DELTAIT3
      ALLOCATABLE DZZMP,DZZPM,KDONE_T,SZZMP,SZZPM
      ALLOCATABLE ZBBZ,ZBSBZ,ZBZ,ZSBZ
      ALLOCATABLE JF,ZF,JG,ZG,JFPLS,JGPLS,ZFPLS,ZGPLS
C
      ALLOCATE (BZZ(NKMMAX,NKMMAX),DB_R2DRDI(NRMAX))
      ALLOCATE (DELTAIT2(NKMMAX,NKMMAX),DELTAIT3(NKMMAX,NKMMAX))
      ALLOCATE (DZZMP(NKMMAX,NKMMAX),DZZPM(NKMMAX,NKMMAX),KDONE_T(NT))
      ALLOCATE (SZZMP(NKMMAX,NKMMAX),SZZPM(NKMMAX,NKMMAX))
      ALLOCATE (ZBBZ(NKMMAX,NKMMAX),ZBSBZ(NKMMAX,NKMMAX))
      ALLOCATE (ZBZ(NKMMAX,NKMMAX),ZSBZ(NKMMAX,NKMMAX))
      ALLOCATE (JFPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JGPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGPLS(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      IF ( IA_ERR.NE.0 ) STOP 'alloc:ovlpint -> ZFP'
C
      WRITE (6,*) 'CHECK_INTEGRATION.EQ.',CHECK_INTEGRATION
C
      IF ( .NOT.GF_CONV_RH ) STOP 
     &                      '<NEGF_TMAT_MO>: GF convention has to be RH'
C
      KDONE_T(:) = .FALSE.
C
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      LOOP_IQ:DO IQ = IQBOT_LBAR,IQTOP_RBAR
         LOOP_IO:DO IO = 1,NOQ(IQ)
C
            IT = ITOQ(IO,IQ)
            IF ( KDONE_T(IT) ) EXIT LOOP_IO
            KDONE_T(IT) = .TRUE.
C
            IM = IMT(IT)
            IRTOP = JRWS(IM)
C
            IF ( CHECK_INTEGRATION.EQ.0 ) THEN
CSW            DB_R2DRDI(1:IRTOP) = -BT(1:IRTOP,IT)*R2DRDI(1:IRTOP,IM)
               DB_R2DRDI(1:IRTOP) = BT(1:IRTOP,IT)*R2DRDI(1:IRTOP,IM)
CSW ??????????????????????????? c scaling??? ???????????????????????????
            ELSE
               DB_R2DRDI(1:IRTOP) = R2DRDI(1:IRTOP,IM)
            END IF
C
            LOOP_IZ:DO IZ = 1,2
C
               IF ( IZ.EQ.1 ) THEN
                  IFIL = IFIL_ZPLS
               ELSE
                  IFIL = IFIL_ZMIN
               END IF
C
               CALL WAVFUN_READ_REL(IFIL,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                              IKMCPLWF)
C
               IF ( IZ.EQ.2 ) CALL WAVFUN_READ_REL(IFIL_ZPLS,IT,0,ZGPLS,
     &              ZFPLS,JGPLS,JFPLS,IRTOP,NCPLWF,IKMCPLWF)
C
               IF ( CHECK_INTEGRATION.EQ.1 ) THEN
                  ZG(:,:,:) = ZGPLS(:,:,:)
                  ZF(:,:,:) = ZFPLS(:,:,:)
               END IF
C
               DZZMP(:,:) = C0
               DZZPM(:,:) = C0
               SZZMP(:,:) = C0
               SZZPM(:,:) = C0
               BZZ(:,:) = C0
               ZBZ(:,:) = C0
               ZBBZ(:,:) = C0
               ZSBZ(:,:) = C0
               ZBSBZ(:,:) = C0
               DELTAIT2(:,:) = C0
               DELTAIT3(:,:) = C0
C
               LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
               DO L = 0,(NLQ(IQ)-1)
C
                  KAP1 = -L - 1
                  KAP2 = L
                  IF ( L.EQ.0 ) KAP2 = KAP1
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                  DO MJM05 = -L - 1, + L
                     MJ = DBLE(MJM05) + 0.5D0
C
                     IF ( ABS(MJ).GE.DBLE(L) ) THEN
                        NSOL = 1
                     ELSE
                        NSOL = 2
                     END IF
C-----------------------------------------------------------------------
                     IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
                     IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C-----------------------------------------------------------------------
                     IKM_SOL(1) = IKM1
                     IKM_SOL(2) = IKM2
                     IMKM_SOL(1) = IKAPMUE(-KAP1,NINT(MJ-0.5D0))
                     IMKM_SOL(2) = IKAPMUE(-KAP2,NINT(MJ-0.5D0))
C-----------------------------------------------------------------------
C
                     CALL RINIT(4,CGF)
C
                     DO I = 1,NSOL
                        IKM = IKM_SOL(I)
                        DO J = 1,NSOL
                           JKM = IKM_SOL(J)
                           CGG(I,J) = AME_G(IKM,JKM,2,ISMT)
                        END DO
C
                        IMKMI = IMKM_SOL(I)
                        CGF(I,I) = AME_G(IMKMI,IMKMI,2,ISMT)
                     END DO
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
                     DO K1 = 1,NSOL
                        LAM1 = IKM_SOL(K1)
C
                        DO K2 = 1,NSOL
                           LAM2 = IKM_SOL(K2)
                           LIN = LIN + 1
C
                           I1 = IKM1LIN(LIN)
                           I2 = IKM2LIN(LIN)
C
C***********************************************************************
C                      <Z(-)|Z(+)> and <Z(+)|Z(-)>
C***********************************************************************
C
                           IF ( IZ.EQ.2 ) THEN
C
                              CALL CINTABR(ZG(1,1,LAM1),ZGPLS(1,1,LAM2),
     &                           ZGZGPLS,ZF(1,1,LAM1),ZFPLS(1,1,LAM2),
     &                           ZFZFPLS,R2DRDI(1,IM),NSOL,NSOL,IRTOP,
     &                           NRMAX)
C
                              CALL CINTABR(ZGPLS(1,1,LAM1),ZG(1,1,LAM2),
     &                           ZGPLSZG,ZFPLS(1,1,LAM1),ZF(1,1,LAM2),
     &                           ZFPLSZF,R2DRDI(1,IM),NSOL,NSOL,IRTOP,
     &                           NRMAX)
C
                              DO KK = 1,NSOL
                                 DZZMP(I1,I2) = DZZMP(I1,I2)
     &                              + ZGZGPLS(KK,KK) + ZFZFPLS(KK,KK)
                                 DZZPM(I1,I2) = DZZPM(I1,I2)
     &                              + ZGPLSZG(KK,KK) + ZFPLSZF(KK,KK)
                              END DO
C
                              DO KKA = 1,NSOL
                                 DO KKB = 1,NSOL
                                    SZZMP(I1,I2) = SZZMP(I1,I2)
     &                                 + CGG(KKA,KKB)*ZGZGPLS(KKA,KKB)
     &                                 - CGF(KKA,KKB)*ZFZFPLS(KKA,KKB)
                                    SZZPM(I1,I2) = SZZPM(I1,I2)
     &                                 + CGG(KKA,KKB)*ZGPLSZG(KKA,KKB)
     &                                 - CGF(KKA,KKB)*ZFPLSZF(KKA,KKB)
                                 END DO
                              END DO
C
C***********************************************************************
C                      <Z(+,r)|beta sigma_z B(r)|Z(+,r)>
C***********************************************************************
                              IF ( STT ) THEN
C
                                 CALL CINTABR(ZGPLS(1,1,LAM1),
     &                              ZGPLS(1,1,LAM2),ZGPLSZGPLS,
     &                              ZFPLS(1,1,LAM1),ZFPLS(1,1,LAM2),
     &                              ZFPLSZFPLS,DB_R2DRDI,NSOL,NSOL,
     &                              IRTOP,NRMAX)
C
                                 DO KK = 1,NSOL
                                    ZBZ(I1,I2) = ZBZ(I1,I2)
     &                                 + ZGPLSZGPLS(KK,KK)
     &                                 + ZFPLSZFPLS(KK,KK)
                                    ZBBZ(I1,I2) = ZBBZ(I1,I2)
     &                                 + ZGPLSZGPLS(KK,KK)
     &                                 - ZFPLSZFPLS(KK,KK)
                                 END DO
C
                                 DO KKA = 1,NSOL
                                    DO KKB = 1,NSOL
                                       ZSBZ(I1,I2) = ZSBZ(I1,I2)
     &                                    + CGG(KKA,KKB)
     &                                    *ZGPLSZGPLS(KKA,KKB)
     &                                    - CGF(KKA,KKB)
     &                                    *ZFPLSZFPLS(KKA,KKB)
                                       ZBSBZ(I1,I2) = ZBSBZ(I1,I2)
     &                                    + CGG(KKA,KKB)
     &                                    *ZGPLSZGPLS(KKA,KKB)
     &                                    + CGF(KKA,KKB)
     &                                    *ZFPLSZFPLS(KKA,KKB)
                                    END DO
                                 END DO
C
                              END IF
C
                           END IF
C
                        END DO
                     END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
                  END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
               END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
               IF ( IZ.EQ.2 ) THEN
C
                  MEZZMP(:,:,IT,1) = DZZMP(:,:)
                  MEZZPM(:,:,IT,1) = DZZPM(:,:)
                  MEZZMP(:,:,IT,2) = SZZMP(:,:)
                  MEZZPM(:,:,IT,2) = SZZPM(:,:)
C
                  IF ( STT ) THEN
CSW                  BZZ(:,:) = ZBZ(:,:)
CSW                  BZZ(:,:) = ZBBZ(:,:)
CSW                  BZZ(:,:) = ZSBZ(:,:)
                     BZZ(:,:) = ZBSBZ(:,:)
                     CALL CMATSTRUCT('BZZ',BZZ(:,:),NKM,NKMMAX,3,3,0,
     &                               1D-8,6)
                  END IF
C
C rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
C   ********************************************************************
C   *                                                                  *
C   *  rotate the matrix elements MEZZMP and MEZZPM and BZZ            *
C   *  from the local to the global frame                              *
C   *                                                                  *
C   *  single site TB t-matrices and MEZZ are already rotated          *
C   *  MEZZMM as well, but currently not used (= L MEZZ^dagger L)      *
C   *                                                                  *
C   ********************************************************************
C
C
                  IF ( KMROT.NE.0 ) THEN
C
                     M = NKMMAX
C
                     IF ( IREL.NE.2 ) THEN
                        N = NKMQ(IQ)
                     ELSE
                        N = NKM
                     END IF
C
                     DO IME = 1,NMEMAX
C---------------------------------------------------------------- MEZZMP
                        CALL ROTATE(MEZZMP(1,1,IT,IME),'L->G',WKM1,N,
     &                              DROTQ(1,1,IQ),M)
                        MEZZMP(1:N,1:N,IT,IME) = WKM1(1:N,1:N)
C
C---------------------------------------------------------------- MEZZPM
                        CALL ROTATE(MEZZPM(1,1,IT,IME),'L->G',WKM1,N,
     &                              DROTQ(1,1,IQ),M)
                        MEZZPM(1:N,1:N,IT,IME) = WKM1(1:N,1:N)
                     END DO
C
C------------------------------------------------------------------- BZZ
                     IF ( STT ) THEN
                        CALL ROTATE(BZZ(1,1),'L->G',WKM1,N,DROTQ(1,1,IQ)
     &                              ,M)
                        BZZ(1:N,1:N) = WKM1(1:N,1:N)
                     END IF
C
                  END IF
C
C rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
C
                  IF ( STT ) THEN
C
C DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
C---------- normalize BZZ by MEZZs to get DELTAIT3, then go to l,m_l,m_s
C-------------------- and take the trace over spin to get DELTA^nrel_LL'
C
                     DO G1 = 1,NKM
                        DO G2 = 1,NKM
                           NORMZG1ZG2 = CDSQRT(MEZZ(G1,G1,IT,1)
     &                                  *MEZZ(G2,G2,IT,1))
                           DELTAIT3(G1,G2) = BZZ(G1,G2)/NORMZG1ZG2
                        END DO
                     END DO
C
                     CALL CHANGEREP(NKM,NKMMAX,DELTAIT3,'REL>RLM',
     &                              DELTAIT2)
                     WRITE (6,*) 'IT = ',IT
                     CALL CMATSTRUCT('DELTAIT',DELTAIT2(:,:),NKM,NKMMAX,
     &                               2,2,0,1D-8,6)
C
                     DO L1 = 1,NLM
                        DO L2 = 1,NLM
                           DNREL2(L1,L2,IT) = DELTAIT2(L1,L2)
     &                        + DELTAIT2(L1+NLM,L2+NLM)
                        END DO
                     END DO
C
C DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
                  END IF
C
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                     check consistency of data sets
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
                  IF ( CHECK_INTEGRATION.EQ.1 .OR. CHECK_ME ) THEN
C
                     WRITE (16,*) 'IT =',IT,':'
                     WRITE (26,*) 'IT =',IT,':'
                     WRITE (17,*) 'IT =',IT,':'
                     WRITE (18,*) 'IT =',IT,':'
                     WRITE (19,*) 'IT =',IT,':'
C
                     CALL CMATSTRUCT('MEZZ_1',MEZZ(:,:,IT,1),NKM,NKMMAX,
     &                               3,3,0,1D-8,16)
                     IF ( STT ) CALL CMATSTRUCT('BZZ',BZZ(:,:),NKM,
     &                    NKMMAX,3,3,0,1D-8,26)
                     CALL CMATSTRUCT('MEZZMM_1',MEZZMM(:,:,IT,1),NKM,
     &                               NKMMAX,3,3,0,1D-8,17)
                     CALL CMATSTRUCT('MEZZMP_1',MEZZMP(:,:,IT,1),NKM,
     &                               NKMMAX,3,3,0,1D-8,18)
                     CALL CMATSTRUCT('MEZZPM_1',MEZZPM(:,:,IT,1),NKM,
     &                               NKMMAX,3,3,0,1D-8,19)
                     CALL CMATSTRUCT('MEZZ_2',MEZZ(:,:,IT,2),NKM,NKMMAX,
     &                               3,3,0,1D-8,16)
                     CALL CMATSTRUCT('MEZZMM_2',MEZZMM(:,:,IT,2),NKM,
     &                               NKMMAX,3,3,0,1D-8,17)
                     CALL CMATSTRUCT('MEZZMP_2',MEZZMP(:,:,IT,2),NKM,
     &                               NKMMAX,3,3,0,1D-8,18)
                     CALL CMATSTRUCT('MEZZPM_2',MEZZPM(:,:,IT,2),NKM,
     &                               NKMMAX,3,3,0,1D-8,19)
C
                  END IF
C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
               END IF
C
            END DO LOOP_IZ
C
         END DO LOOP_IO
      END DO LOOP_IQ
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      IF ( STT ) THEN
         DO IQ = IQBOT_LBAR,IQTOP_RBAR
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               WRITE (6,*) 'IT = ',IT
               CALL CMATSTRUCT('DNREL2',DNREL2(:,:,IT),NLM,NLMMAX,1,1,0,
     &                         1D-8,6)
            END DO
         END DO
      END IF
C
      DEALLOCATE (ZF,ZG,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc <NEGF_TMAT>  -> ZG'
C
      END
