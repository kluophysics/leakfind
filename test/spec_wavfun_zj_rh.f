C*==wavfun_zj_rh.f    processed by SPAG 6.70Rc at 16:37 on 28 Feb 2017
      SUBROUTINE WAVFUN_ZJ_RH(IFILCBWF,IT,P,MSST,TSST,RG,RF,HG,HF,
     &                        ISTATE,RG1,RF1,NZERO)
C   ********************************************************************
C   *                                                                  *
C   *      read regular und irregular  wave functions  Z  and  J       *
C   *      and transfer to Juelich convention          R  and  H       *
C   *                                                                  *
C   *                ASA and FULL POTENTIAL version                    *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C0,CI
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_RMESH,ONLY:FULLPOT,JRWS,JRCRI,NRMAX
      USE MOD_TYPES,ONLY:NT,IMT,NCPLWFMAX,NTMAX
      USE MOD_SITES,ONLY:DROTQ,IQAT
      USE MOD_CALCMODE,ONLY:KMROT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WAVFUN_ZJ_RH')
C
C Dummy arguments
C
      INTEGER IFILCBWF,ISTATE,IT
      COMPLEX*16 P
      COMPLEX*16 HF(NRMAX,NKMMAX,NKMMAX),HG(NRMAX,NKMMAX,NKMMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),RF(NRMAX,NKMMAX,NKMMAX),
     &           RF1(NRMAX,NKMMAX,NKMMAX),RG(NRMAX,NKMMAX,NKMMAX),
     &           RG1(NRMAX,NKMMAX,NKMMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      REAL*8 NZERO(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER IKM,IKMCPLWF(:,:),IM,IQ,IR,IRTOP,IRTOPP,ISOL,ISOLP,ITP,
     &        JSOL,K,M,NKM_SOL(:),NOSWF
      COMPLEX*16 JF(:,:,:),JG(:,:,:),M_JI,T_JI,W1(:,:),W2(:,:),ZF(:,:,:)
     &           ,ZG(:,:,:)
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZG,ZF,JG,JF,W1,W2
      ALLOCATABLE IKMCPLWF,NKM_SOL
C
      M = NKMMAX
C
      ALLOCATE (NKM_SOL(NKMMAX))
      IF ( FULLPOT ) THEN
         ALLOCATE (IKMCPLWF(NCPLWFMAX,NKMMAX))
      ELSE
         ALLOCATE (IKMCPLWF(2,NKMMAX))
      END IF
C
      ALLOCATE (ZG(NRMAX,M,M),ZF(NRMAX,M,M))
      ALLOCATE (JG(NRMAX,M,M),JF(NRMAX,M,M))
C
      IM = IMT(IT)
C
      NOSWF = NT*NKM
C
C=======================================================================
C         READ    wave functions into  NKM x NKM matrix arrays
C=======================================================================
C
      IF ( FULLPOT ) THEN
         IRTOP = JRCRI(IM)
      ELSE
         IRTOP = JRWS(IM)
      END IF
C
      ZG(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      ZF(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      JG(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      JF(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
C
      DO ISOL = 1,NKM
C
         READ (IFILCBWF,REC=ISOL+(IT-1)*NKM) ITP,STRP,ISOLP,IRTOPP,
     &         NKM_SOL(ISOL),(IKMCPLWF(K,ISOL),K=1,NKM_SOL(ISOL)),
     &         ((ZG(IR,IKMCPLWF(K,ISOL),ISOL),IR=1,IRTOP),K=1,
     &         NKM_SOL(ISOL)),
     &         ((ZF(IR,IKMCPLWF(K,ISOL),ISOL),IR=1,IRTOP),K=1,
     &         NKM_SOL(ISOL))
C
         IF ( STRP.NE.'REG' .OR. IT.NE.ITP .OR. IRTOP.NE.IRTOPP ) THEN
            WRITE (6,*) 'IT STRP ISOL',IT,ITP,'REG',STRP,ISOL,ISOLP
            CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading REG WF')
         END IF
C
         READ (IFILCBWF,REC=ISOL+(IT-1)*NKM+NOSWF) ITP,STRP,ISOLP,
     &         IRTOPP,((JG(IR,IKMCPLWF(K,ISOL),ISOL),IR=1,IRTOP),K=1,
     &         NKM_SOL(ISOL)),
     &         ((JF(IR,IKMCPLWF(K,ISOL),ISOL),IR=1,IRTOP),K=1,
     &         NKM_SOL(ISOL))
C
         IF ( STRP.NE.'IRR' .OR. IT.NE.ITP .OR. IRTOP.NE.IRTOPP ) THEN
            WRITE (6,*) 'IT STRP ISOL',IT,ITP,'IRR',STRP,ISOL,ISOLP
            CALL STOP_MESSAGE(ROUTINE,'TROUBLE reading IRR WF')
         END IF
C
      END DO
C
C=======================================================================
C
C
C
C=======================================================================
C         convert  (Z,J)  to the Juelich convention  (R,H)
C=======================================================================
C
C-------------------------------------------------- R_i = SUM_j Z_j t_ji
C
      RG(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      RF(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      RG1(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      RF1(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
C
      NZERO = 0.0D0
      DO ISOL = 1,NKM
         DO JSOL = 1,NKM
            T_JI = TSST(JSOL,ISOL,IT)
            IF ( ABS(T_JI).GE.1D-12 ) NZERO(JSOL,ISOL) = 1.0D0
         END DO
      END DO
C
      DO ISOL = 1,NKM
         DO JSOL = 1,NKM
            T_JI = TSST(JSOL,ISOL,IT)
            IF ( ABS(T_JI).GE.1D-16 ) THEN
               DO IKM = 1,NKM
                  DO IR = 1,IRTOP
                     RG(IR,IKM,ISOL) = RG(IR,IKM,ISOL) + ZG(IR,IKM,JSOL)
     &                                 *T_JI
                     RF(IR,IKM,ISOL) = RF(IR,IKM,ISOL) + ZF(IR,IKM,JSOL)
     &                                 *T_JI
                  END DO
               END DO
            END IF
         END DO
      END DO
C
      IF ( ISTATE.EQ.1 ) THEN
         DO ISOL = 1,NKM
            DO IKM = 1,NKM
               DO IR = 1,IRTOP
                  RG1(IR,IKM,ISOL) = ZG(IR,IKM,ISOL)
                  RF1(IR,IKM,ISOL) = ZF(IR,IKM,ISOL)
               END DO
            END DO
         END DO
      END IF
C
C---------------------------------- H_i = SUM_j (R_j - J_j) m_ji / (-ip)
C
      HG(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
      HF(1:NRMAX,1:NKMMAX,1:NKMMAX) = C0
C
      DO ISOL = 1,NKM
         DO JSOL = 1,NKM
            M_JI = MSST(JSOL,ISOL,IT)/(-CI*P)
            IF ( ABS(M_JI).GE.1D-10 ) THEN
               DO IKM = 1,NKM
                  DO IR = 1,IRTOP
                     HG(IR,IKM,ISOL) = HG(IR,IKM,ISOL)
     &                                 + (RG(IR,IKM,JSOL)
     &                                 -JG(IR,IKM,JSOL))*M_JI
                     HF(IR,IKM,ISOL) = HF(IR,IKM,ISOL)
     &                                 + (RF(IR,IKM,JSOL)
     &                                 -JF(IR,IKM,JSOL))*M_JI
                  END DO
               END DO
            END IF
         END DO
      END DO
C
C=======================================================================
C        rotate the single site t-matrix and radial wave functions
C                 from local to global frame if necessary
C=======================================================================
C
      IF ( KMROT.GT.0 ) THEN
C
         ALLOCATE (W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX))
         IQ = IQAT(1,IT)
C
         W2(1:NKMMAX,1:NKMMAX) = TSST(1:NKMMAX,1:NKMMAX,IT)
         CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
C
         NZERO = 0.0D0
         DO ISOL = 1,NKM
            DO JSOL = 1,NKM
               T_JI = W1(JSOL,ISOL)
               IF ( ABS(T_JI).GE.1D-12 ) NZERO(JSOL,ISOL) = 1.0D0
            END DO
         END DO
C
         DO IR = 1,IRTOP
            W2(1:NKMMAX,1:NKMMAX) = RG(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            RG(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
C
C
            W2(1:NKMMAX,1:NKMMAX) = RF(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            RF(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
C
C
            W2(1:NKMMAX,1:NKMMAX) = RG1(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            RG1(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
C
C
            W2(1:NKMMAX,1:NKMMAX) = RF1(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            RF1(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
C
C
            W2(1:NKMMAX,1:NKMMAX) = HG(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            HG(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
C
C
            W2(1:NKMMAX,1:NKMMAX) = HF(IR,1:NKMMAX,1:NKMMAX)
            CALL ROTATE(W2,'L->G',W1,NKM,DROTQ(1,1,IQ),NKM)
            HF(IR,1:NKMMAX,1:NKMMAX) = W1(1:NKMMAX,1:NKMMAX)
         END DO
C
      END IF
C
      END
