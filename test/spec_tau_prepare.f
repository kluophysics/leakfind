C*==spec_tau_prepare.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_TAU_PREPARE(ISTATE,JE,NEW,LAYS,NATL,ESTAT,
     &                            UMAT_VT_PES,TAUMT,TSSTS,MSSTS,GAMMA1,
     &                            CPAPROJ,CPAATOM,P,MEZZ,MEZJ,MSSQ,SSST,
     &                            TAUQ,TAUT,PHASK,NOUT1,CALCDOS,ATA,IP,
     &                            PF,IBLOCH)
C
C   ********************************************************************
C   *                                                                  *
C   *  Interface to SPR-KKR                                            *
C   *  Main outpu of this routine are:                                 *
C   *  GAMMA1: Single site matrix                                      *
C   *  TAUMT: Tauq-t                                                   *
C   *  MSSTS: MSST for ISTATE                                          *
C   *  CPAPROJ: cpa projector                                          *
C   *  UMAT_VT_PES: U matrix for dislocations                          *
C   ********************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,NVFTPHOMAX
      USE MOD_THERMAL,ONLY:X_VFT,NVFO_Q,IVFT_VFOQ,NVFTMAX,UMAT_VT
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT
      USE MOD_CONSTANTS,ONLY:C0,C1
      USE MOD_CALCMODE,ONLY:THERMAL_VIBRA_FLUCT
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NQMAX
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_CPA,ONLY:NCPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATA,IBLOCH,IP,ISTATE,JE,LAYS,NEW,NOUT1
      LOGICAL CALCDOS
      COMPLEX*16 ESTAT,PF
      INTEGER CPAATOM(LAYSM,NATLM),NATL(LAYSM)
      COMPLEX*16 CPAPROJ(LAYSM,NATLM,MQD,MQD,2,NVFTPHOMAX),
     &           GAMMA1(LAYSM,NATLM,3,MQD,MQD),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSSTS(NKMMAX,NKMMAX,NTMAX,2),
     &           P(2),PHASK(NEMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUMT(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX+1),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TSSTS(NKMMAX,NKMMAX,NTMAX,2),
     &           UMAT_VT_PES(LAYSM,NATLM,MQD,MQD,NVFTPHOMAX)
C
C Local variables
C
      INTEGER ATOM,COUNTING,I,IQ,ITS,J,K,LAY
      COMPLEX*16 C2(:,:),C3(:,:),CM1(:,:),DMATTG(:,:,:),
     &           DMATTG_VFT(:,:,:),DOS(:,:,:),DTILTG(:,:,:),
     &           DTILTG_VFT(:,:,:),TAUT_GLO(:,:,:),TAUT_GLO_VFT(:,:,:),
     &           TSSQ(:,:,:),TSST_GLO(:,:,:),TSST_GLO_VFT(:,:,:)
      REAL*8 CONZ
C
C*** End of declarations rewritten by SPAG
C
C
C
      ALLOCATABLE C2,C3,CM1,DMATTG
      ALLOCATABLE     DMATTG_VFT,DOS,DTILTG
      ALLOCATABLE     DTILTG_VFT
      ALLOCATABLE     TAUT_GLO,TAUT_GLO_VFT
      ALLOCATABLE     TSSQ,TSST_GLO
      ALLOCATABLE     TSST_GLO_VFT
C
      ALLOCATE (CM1(MQD,MQD),C2(MQD,MQD),C3(MQD,MQD))
      ALLOCATE (TSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (DMATTG(NKMMAX,NKMMAX,NTMAX),DTILTG(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TSST_GLO(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TAUT_GLO(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (DOS(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TSST_GLO_VFT(NKMMAX,NKMMAX,NVFTMAX))
      ALLOCATE (TAUT_GLO_VFT(NKMMAX,NKMMAX,NVFTMAX))
      ALLOCATE (DMATTG_VFT(NKMMAX,NKMMAX,NVFTMAX))
      ALLOCATE (DTILTG_VFT(NKMMAX,NKMMAX,NVFTMAX))
C
      TSST_GLO = C0
      TAUT_GLO = C0
      TSST_GLO_VFT = C0
      TAUT_GLO_VFT = C0
      DMATTG_VFT = C0
      DTILTG_VFT = C0
      CM1 = C0
      C2 = C0
      C3 = C0
      DOS = C0
C
C ======================================================================
C     Initilise UMAT and CPAPROJ to diagonal C1 matrix
C     This is done in particular for ATAVIB=1
C ======================================================================
      UMAT_VT_PES = C0
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            IQ = IQ_SPR_LAYAT(ATOM,LAY)
            DO K = 1,NVFO_Q(IQ)
               FORALL(I=1:MQD)UMAT_VT_PES(LAY,ATOM,I,I,K) = C1
               FORALL(I=1:MQD)CPAPROJ(LAY,ATOM,I,I,ISTATE,K) = C1
            END DO
         END DO
      END DO
C
C ======================================================================
C     Driver to KKR single site and multiple scattering part
C ======================================================================
C
      IF ( THERMAL_VIBRA_FLUCT ) THEN
C
         STOP 'SPEC: THERMAL_VIBRA_FLUCT NOT YET IN THIS VERSION'
C         CALL SPEC_TAU_DRIVE_THERMAL(MEZJ,MEZZ,MSSQ,TSSQ,
C     &                               MSSTS(1,1,1,ISTATE),SSST,TAUQ,
C     &                               TSSTS(1,1,1,ISTATE),PHASK,
C     &                               DMATTG_VFT,DTILTG_VFT,2.0D0*ESTAT,
C     &                               ISTATE,JE,TSST_GLO_VFT,
C     &                               TAUT_GLO_VFT,DOS,CALCDOS,P(ISTATE),
C     &                               ATA)
CC
C         DO LAY = 1,LAYS
C            DO ATOM = 1,NATL(LAY)
C               IQ = IQ_SPR_LAYAT(ATOM,LAY)
C               DO K = 1,NVFO_Q(IQ)
C                  ITS = IVFT_VFOQ(K,IQ)
C                  CM1(:,:) = UMAT_VT(:,:,ITS)
C                  CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,0,PF)
C                  UMAT_VT_PES(LAY,ATOM,:,:,K) = C2(:,:)
C                  IF ( IP.GT.1 ) CALL  CMATSTR('UMAT_PES',C2,MQD,MQD,3,
C     &                 3,1,1D-8,NOUT1)
C               END DO
C            END DO
C         END DO
      ELSE
         CALL SPEC_TAU_DRIVE(MEZJ,MEZZ,MSSQ,TSSQ,MSSTS(1,1,1,ISTATE),
     &                       SSST,TAUQ,TAUT,TSSTS(1,1,1,ISTATE),PHASK,
     &                       DMATTG,DTILTG,2.0D0*ESTAT,ISTATE,JE,
     &                       TSST_GLO,TAUT_GLO,DOS,CALCDOS,NEW,P(ISTATE)
     &                       )
C
         TSST_GLO_VFT = TSST_GLO(:,:,1:NTMAX)
         TAUT_GLO_VFT = TAUT_GLO(:,:,1:NTMAX)
         DMATTG_VFT = DMATTG(:,:,1:NTMAX)
         DTILTG_VFT = DTILTG(:,:,1:NTMAX)
C
      END IF
C
C ======================================================================
C     Fine which atoms in layers are CPA atoms
C     In case of AI-PES force CPA
C ======================================================================
      IF ( NCPA.GT.0 ) THEN
         DO LAY = 1,LAYS
            DO ATOM = 1,NATL(LAY)
               IQ = IQ_SPR_LAYAT(ATOM,LAY)
               IF ( NVFO_Q(IQ).EQ.1 ) THEN
                  CPAATOM(LAY,ATOM) = 0
               ELSE IF ( NVFO_Q(IQ).GT.1 ) THEN
                  CPAATOM(LAY,ATOM) = 1
               END IF
            END DO
         END DO
      END IF
      IF ( IBLOCH.EQ.2 ) THEN
         DO LAY = 1,LAYS
            DO ATOM = 1,NATL(LAY)
               IQ = IQ_SPR_LAYAT(ATOM,LAY)
               CPAATOM(LAY,ATOM) = 1
            END DO
         END DO
      END IF
C
C ======================================================================
C
C
      IF ( IBLOCH.GT.1 ) THEN
         IF ( NCPA.GT.0 ) THEN
            DO LAY = 1,LAYS
               DO ATOM = 1,NATL(LAY)
                  IQ = IQ_SPR_LAYAT(ATOM,LAY)
                  DO K = 1,NVFO_Q(IQ)
                     ITS = IVFT_VFOQ(K,IQ)
                     CM1(:,:) = DMATTG_VFT(:,:,ITS)
                     CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,0,PF)
                     CPAPROJ(LAY,ATOM,:,:,ISTATE,K) = C2(:,:)
                  END DO
               END DO
            END DO
         END IF
         IF ( ISTATE.EQ.1 ) THEN
            DO LAY = 1,LAYS
               DO ATOM = 1,NATL(LAY)
                  IQ = IQ_SPR_LAYAT(ATOM,LAY)
C
                  IF ( CPAATOM(LAY,ATOM).EQ.1 ) THEN
C
                     DO K = 1,NVFO_Q(IQ)
                        ITS = IVFT_VFOQ(K,IQ)
                        CM1(:,:) = TAUT_GLO_VFT(:,:,ITS)
     &                             - TSST_GLO_VFT(:,:,ITS)
                        CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
                        TAUMT(LAY,ATOM,:,:,K) = C2(:,:)
                     END DO
C
                     CM1(:,:) = TSSQ(:,:,IQ)
                     CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
                     CALL INV(C2,C3,MQD)
C
                     CM1(:,:) = TAUQ(:,:,IQ) - TSSQ(:,:,IQ)
                     CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
                     CALL MULT(C3,C2,CM1,MQD)
                     CALL MULT(CM1,C3,C2,MQD)
                     TAUMT(LAY,ATOM,:,:,NVFO_Q(IQ)+1) = C2(:,:)
C
                  END IF
               END DO
            END DO
         END IF
      END IF
C
      DO LAY = 1,LAYS
         DO ATOM = 1,NATL(LAY)
            IQ = IQ_SPR_LAYAT(ATOM,LAY)
            CM1 = C0
            IF ( NCPA.GT.0 ) THEN
               IF ( ATA.EQ.1 ) THEN
                  IQ = IQ_SPR_LAYAT(ATOM,LAY)
                  DO K = 1,NVFO_Q(IQ)
                     ITS = IVFT_VFOQ(K,IQ)
                     CONZ = X_VFT(ITS)
                     CM1(:,:) = CM1(:,:) + CONZ*TSST_GLO_VFT(:,:,ITS)
                  END DO
                  CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
C
                  GAMMA1(LAY,ATOM,ISTATE,:,:) = C2(:,:)
               ELSE IF ( ATA.EQ.0 ) THEN
C
                  IF ( CPAATOM(LAY,ATOM).EQ.1 ) THEN
                     CM1(:,:) = TSSQ(:,:,IQ)
                     CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
                     GAMMA1(LAY,ATOM,ISTATE,:,:) = C2(:,:)
                  ELSE
                     DO K = 1,NVFO_Q(IQ)
                        ITS = IVFT_VFOQ(K,IQ)
                        CM1(:,:) = TSST_GLO_VFT(:,:,ITS)
                     END DO
C????
                     CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
C
                     GAMMA1(LAY,ATOM,ISTATE,:,:) = C2(:,:)
                  END IF
C
               END IF
C
            ELSE IF ( NCPA.EQ.0 ) THEN
               CM1(:,:) = TSSQ(:,:,IQ)
               CALL SPEC_TRANSPHO(CM1,C2,NKMMAX,1,PF)
               GAMMA1(LAY,ATOM,ISTATE,:,:) = C2(:,:)
            END IF
C
            IF ( IP.GE.2 ) THEN
               IF ( ISTATE.EQ.1 ) THEN
                  WRITE (NOUT1,99001)
               ELSE IF ( ISTATE.EQ.2 ) THEN
                  WRITE (NOUT1,99002)
               END IF
               WRITE (NOUT1,99003) ATOM,LAY,ESTAT
               COUNTING = 0
               DO I = 1,MQD
                  DO J = 1,MQD
                     IF ( CDABS(GAMMA1(LAY,ATOM,ISTATE,I,J)).GT.1.D-16 )
     &                    THEN
                        COUNTING = COUNTING + 1
                        WRITE (NOUT1,99004) COUNTING,I,J,
     &                         GAMMA1(LAY,ATOM,ISTATE,I,J),
     &                         GAMMA1(LAY,ATOM,ISTATE,J,I)
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
C
C
99001 FORMAT (2x,'****** initial state ******')
99002 FORMAT (2x,'****** final   state ******')
99003 FORMAT (2x,'scattering matrix for atom',2x,i3,2x,'layer',2x,i3,2x,
     &        'erel',2x,2(1x,e15.8))
99004 FORMAT (3(2x,i3),2x,2E14.7,2x,2E14.7)
C
C
      END
C
