C*==fullf.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE FULLF(LAY,ATOM,IO,ISTATE,EREL,CLIGHT,USEEULER,IBLOCH,
     &                 JE,VLM,MAXG,G1,G2,TSST,MSST,P)
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC
      USE MOD_SPEC_WAVE,ONLY:WFF,WFFX,WFFY,WFFM,WFFMATOM,WFFMZ,WFFXM,
     &    WFFYM,FFI_G1M,FFI_S1M,FF_G1,FF_G1M,FF_G1MZ,FF_S1,FF_S1M,
     &    FF_S1MZ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IBLOCH,IO,ISTATE,JE,LAY,MAXG,USEEULER
      REAL*8 CLIGHT
      COMPLEX*16 EREL,P
      INTEGER G1(MQK),G2(MQK)
      COMPLEX*16 MSST(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
C
C Local variables
C
      COMPLEX*16 CC(:,:,:),FF(:),IMM(:),IMP(:),IPM(:),IPP(:),SS(:,:,:)
      INTEGER CORESELRUL(:,:),FPSELRUL(:,:),G1_KKR(:),G1_NEW(:),
     &        G2_KKR(:),G2_NEW(:),I,IR,J,KAP(:),MAXG_KKR,MAXG_NEW,MAXPM,
     &        MAXPP,MUD(:),NR(:),NRL(:),NRM(:),PM1(:),PM2(:),PM3(:),
     &        PP1(:),PP2(:),PP3(:),REL(:,:),SELRUL(:,:)
      REAL*8 GRID(:),NRS(:),VORZ
C
C*** End of declarations rewritten by SPAG
C
C     ==================================================================
C     full-potential wave functions:
C     phase functions    :
C     regular solution   : phir <=> {\Phi}_{lm l'm'} (r)
C     irregular solution : phii <=> {\Phi}^{+}_{lm l'm'} (r)
C     derivatives of the phase functions    :
C     regular solution   : phirs <=> {\Phi}'_{lm l'm'} (r)
C     irregular solution : phiis <=> {\Phi}'^{+}_{lm l'm'} (r)
C
C     ==================================================================
C
      ALLOCATABLE FPSELRUL,CC,FF,NR,PM1,PM2,PM3,PP1,PP2,PP3
      ALLOCATABLE CORESELRUL,KAP,REL,MUD,NRL,NRM,NRS,GRID,G1_NEW
      ALLOCATABLE G2_NEW,SELRUL,G1_KKR,G2_KKR,SS
      ALLOCATABLE IMM,IMP,IPM,IPP
      ALLOCATE (FPSELRUL(MQD,MQD),CC(MQD,MQD,7),FF(RSTEP),NR(MQD))
      ALLOCATE (PM1(MAXIPP),PM2(MAXIPP),PM3(MAXIPP),PP1(MAXIPP))
      ALLOCATE (PP2(MAXIPP),PP3(MAXIPP),CORESELRUL(MQD,MQD))
      ALLOCATE (KAP(MQD),REL(MQD,2),MUD(MQD),NRL(MQD),NRM(MQD))
      ALLOCATE (NRS(MQD),GRID(RSTEPP5),G1_NEW(MQK),G2_NEW(MQK))
      ALLOCATE (SELRUL(MQD,MQD),SS(MQD,MQD,7))
      ALLOCATE (IMM(MAXIPP),IMP(MAXIPP),IPM(MAXIPP),IPP(MAXIPP))
C
      IF ( .NOT.ALLOCATED(FF_S1M) ) THEN
         ALLOCATE (FF_S1M(NRMAX,MQD,MQD),FF_G1M(NRMAX,MQD,MQD))
         ALLOCATE (FF_S1MZ(NRMAX,MQD,MQD))
         ALLOCATE (FF_G1MZ(NRMAX,MQD,MQD))
         ALLOCATE (FFI_S1M(NRMAX,MQD,MQD))
         ALLOCATE (FFI_G1M(NRMAX,MQD,MQD))
         ALLOCATE (FF_G1(RSTEPP5,MQD,MQD),FF_S1(RSTEPP5,MQD,MQD))
         FF_S1M = DCMPLX(0.0D0,0.0D0)
         FF_G1M = DCMPLX(0.0D0,0.0D0)
         FF_S1MZ = DCMPLX(0.0D0,0.0D0)
         FF_G1MZ = DCMPLX(0.0D0,0.0D0)
         FFI_S1M = DCMPLX(0.0D0,0.0D0)
         FFI_G1M = DCMPLX(0.0D0,0.0D0)
         FF_G1 = DCMPLX(0.0D0,0.0D0)
         FF_G1 = DCMPLX(0.0D0,0.0D0)
      END IF
C
C
      CALL RIND(ML,MQD,KAP,MUD,NRL,NRM,NRS,NR,REL)
C
      DO I = 1,MQK
         G1(I) = 0
         G2(I) = 0
      END DO
C
      CALL MAGFPSELRUL(LAY,ATOM,USEEULER,CORESELRUL,FPSELRUL,IO,VLM)
C
      CALL INDEXF(LAY,ATOM,IPP,IPM,IMP,IMM,PP1,PP2,PP3,PM1,PM2,PM3,
     &            MAXPP,MAXPM,SELRUL,IO,VLM)
C
      VORZ = 1.0D0
      CC = CZERO
      SS = CZERO
      GRID = 0.0D0
C
      CALL INITREF(LAY,ATOM,CC,SS,GRID,VORZ,RSTEP)
C
      CALL REGULARF(LAY,ATOM,CC,SS,G1,G2,MAXG,GRID,VORZ,PP1,PP2,PP3,PM1,
     &              PM2,PM3,MAXPP,MAXPM,KAP,NRL,IPP,IPM,IMP,IMM,EREL,
     &              CLIGHT,G1_NEW,G2_NEW,MAXG_NEW,FF_S1,FF_G1,IO,VLM)
C
      ALLOCATE (G1_KKR(MQK),G2_KKR(MQK))
      G1_KKR(:) = 0
      G2_KKR(:) = 0
C      MAXG_KKR = 0
C
      CALL SCATTERF(LAY,ATOM,IO,ISTATE,FF_S1M,FF_G1M,FFI_S1M,FFI_G1M,JE,
     &              FF_S1MZ,FF_G1MZ,G1_KKR,G2_KKR,MAXG_KKR,TSST,MSST,P)
C
C      WRITE(6,*) "-------------------------- Coupling scheme: KKR"
C      write (6,*) 'MAXG KKR',maxg_kkr
C      DO IR=1,MAXG_KKR
C         WRITE(6,*) IR,G1_KKR(IR),G2_KKR(IR)
C      END DO
C      WRITE(6,*) "-------------------------- Coupling scheme: SPEC"
C      DO IR=1,MAXG_NEW
C         WRITE(6,*) IR,G1_NEW(IR),G2_NEW(IR)
C      END DO
C      IF (MAXG_NEW .NE. MAXG_KKR) THEN
C         WRITE(6,*)
C     &        "WARNING: MAXG_NEW .NE. MAXG_KKR"
C      ELSE
C         DO IR=1,MAXG_KKR
C            IF (G1_NEW(IR).NE.G1_KKR(IR) .OR. G2_NEW(IR).NE.G2_KKR(IR))
C     &           THEN
C               WRITE(6,*) "WARNING: SPEC G1,G2:", IR,G1_NEW(IR),
C     &              G2_NEW(IR)
C               WRITE(6,*) "WARNING: KKR G1,G2:", IR,G1_KKR(IR),
C     &              G2_KKR(IR)
C            END IF
C         END DO
C      END IF
C
C     MAXG_NEW = MAXG_KKR
      MAXG_NEW = MAXG
C
      DO IR = 1,MAXG_NEW
C        G1_NEW(IR) = G1_KKR(IR)
C        G2_NEW(IR) = G2_KKR(IR)
         G1_NEW(IR) = G1(IR)
         G2_NEW(IR) = G2(IR)
      END DO
C     MAXG = MAXG_KKR
      MAXG_NEW = MAXG
      DO IR = 1,MAXG_NEW
C        G1(IR) = G1_KKR(IR)
C        G2(IR) = G2_KKR(IR)
         G1(IR) = G1_NEW(IR)
         G2(IR) = G2_NEW(IR)
      END DO
C
C
C      write (*,*) 'MAXG',maxg,maxg_new,MAXWF
C      write (*,*) 'G1G2', G1,G2
C      write (*,*) 'G1_NEW,G2_NEW',G1_NEW,G2_NEW
      CALL PRINT_WAVEM(FF_S1M,FF_G1M,G1_NEW,G2_NEW,MAXG_NEW,LAY,ATOM,
     &                 ISTATE,IO)
C
C
C
CC Setup MAXWF to actual value of MAXG
C      IF ( .NOT.ALLOCATED(WFFM) ) THEN
CC Setup MAXWF to actual value of MAXG
C         MAXWF = MAXG
CC
C         ALLOCATE (WFFM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
C         ALLOCATE (WFFMZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
C         ALLOCATE (WFFMATOM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
C         ALLOCATE (WFFXM(MQD,NSTATES,2,NTPHOMAX))
C         ALLOCATE (WFFYM(MAXWF,NSTATES,NTPHOMAX))
C         ALLOCATE (WFF(RSTEP,MAXWF,NSTATES,2))
C         ALLOCATE (WFFX(MQD,NSTATES,2),WFFY(MAXWF,NSTATES))
CC
C      END IF
      MAXWF = MAXG
      IF ( ISTATE.EQ.1 ) THEN
         IF ( .NOT.ALLOCATED(WFFM) ) THEN
C
            ALLOCATE (WFFM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFMZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFMATOM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFXM(MQD,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFYM(MAXWF,NSTATES,NTPHOMAX))
            ALLOCATE (WFF(RSTEP,MAXWF,NSTATES,2))
            ALLOCATE (WFFX(MQD,NSTATES,2),WFFY(MAXWF,NSTATES))
         ELSE
            DEALLOCATE (WFFM)
            DEALLOCATE (WFFMZ)
            DEALLOCATE (WFFMATOM)
            DEALLOCATE (WFFXM)
            DEALLOCATE (WFFYM)
            DEALLOCATE (WFF)
            DEALLOCATE (WFFX,WFFY)
C
            ALLOCATE (WFFM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFMZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFMATOM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFXM(MQD,NSTATES,2,NTPHOMAX))
            ALLOCATE (WFFYM(MAXWF,NSTATES,NTPHOMAX))
            ALLOCATE (WFF(RSTEP,MAXWF,NSTATES,2))
            ALLOCATE (WFFX(MQD,NSTATES,2),WFFY(MAXWF,NSTATES))
C
         END IF
      END IF
C
      CALL CPXPHIWFFM(WFFM,WFFXM,WFFYM,ISTATE,FF_S1M,FF_G1M,G1_NEW,
     &                G2_NEW,MAXG,ATOM,LAY,WFFMATOM,IO,FF_S1MZ,FF_G1MZ,
     &                WFFMZ,IBLOCH)
C
      IF ( ISTATE.EQ.1 .AND. IBLOCH.GT.1 ) THEN
C
         DO J = 1,RSTEP
            FF(J) = GRID(RSTEPP5+1-J)
         END DO
         DO J = 1,RSTEP
            GRID(J) = REAL(FF(J))
         END DO
C
         CALL PRINT_IRRWAVE(FFI_S1M,FFI_G1M,G1,G2,MAXG,LAY,ATOM,3,IO)
C
         CALL CPXPHIWFFM(WFFM,WFFXM,WFFYM,3,FFI_S1M,FFI_G1M,G1_NEW,
     &                   G2_NEW,MAXG,ATOM,LAY,WFFMATOM,IO,FF_S1MZ,
     &                   FF_G1MZ,WFFMZ,IBLOCH)
      END IF
C
      END
C*==cpxphiwffm.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE CPXPHIWFFM(WFF,WFFX,WFFY,STATE,FF_S1,FF_G1,G1,G2,MAXG,
     &                      ATOM,LAY,WFFATOM,IO,FF_S1Z,FF_G1Z,WFFZ,
     &                      IBLOCH)
C
C     # purpose      : transform matrix representation of munich
C                      wavefunction (ff_s1m,ff_g1m)
C                      to vector representation (wffm, wffxm, wffym)
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC
      USE MOD_SPEC_MESH,ONLY:RADMESHM
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_RINDC,ONLY:KAP,NRM,NRS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IBLOCH,IO,LAY,MAXG,STATE
      COMPLEX*16 FF_G1(NRMAX,MQD,MQD),FF_G1Z(NRMAX,MQD,MQD),
     &           FF_S1(NRMAX,MQD,MQD),FF_S1Z(NRMAX,MQD,MQD),
     &           WFF(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),
     &           WFFATOM(NRMAX,MAXWF,NSTATES,2,NTPHOMAX),
     &           WFFZ(NRMAX,MAXWF,NSTATES,2,NTPHOMAX)
      INTEGER G1(MQK),G2(MQK),WFFX(MQD,NSTATES,2,NTPHOMAX),
     &        WFFY(MAXWF,NSTATES,NTPHOMAX)
C
C Local variables
C
      INTEGER I,J,LAST
      REAL*8 MU,MUS
C
C*** End of declarations rewritten by SPAG
C
C      DO I = 1,MQD
C         WFFX(I,STATE,1,IO) = 0
C         WFFX(I,STATE,2,IO) = 0
C      END DO
      WFFX(:,STATE,:,IO) = 0
C
      WFFY(:,STATE,IO) = 0
      WFF(:,:,STATE,1,IO) = CZERO
      WFF(:,:,STATE,2,IO) = CZERO
      WFFZ(:,:,STATE,1,IO) = CZERO
      WFFZ(:,:,STATE,2,IO) = CZERO
      WFFATOM(:,:,STATE,1,IO) = CZERO
      WFFATOM(:,:,STATE,2,IO) = CZERO
C      DO I = 1,MAXWF
C         WFFY(I,STATE,IO) = 0
C         DO J = 1,NRMAX
C            WFF(J,I,STATE,1,IO) = CZERO
C            WFF(J,I,STATE,2,IO) = CZERO
C            WFFZ(J,I,STATE,1,IO) = CZERO
C            WFFZ(J,I,STATE,2,IO) = CZERO
C            WFFATOM(J,I,STATE,1,IO) = CZERO
C            WFFATOM(J,I,STATE,2,IO) = CZERO
C         END DO
C      END DO
C
      LAST = 0
      DO I = 1,MAXG
         IF ( G1(I).NE.LAST ) THEN
            WFFX(G1(I),STATE,1,IO) = I
            LAST = G1(I)
         END IF
         WFFX(G1(I),STATE,2,IO) = I
         WFFY(I,STATE,IO) = G2(I)
      END DO
C
      DO I = 1,MAXG
         DO J = 1,RSTEP
            WFFATOM(J,I,STATE,1,IO) = FF_S1(J,G2(I),G1(I))
            WFFATOM(J,I,STATE,2,IO) = FF_G1(J,G2(I),G1(I))
         END DO
      END DO
C
      IF ( STATE.EQ.1 .OR. STATE.EQ.3 ) THEN
C
         DO I = 1,MAXG
            DO J = 1,RSTEP
               IF ( IBLOCH.EQ.2 .AND. STATE.EQ.1 ) THEN
                  WFF(J,I,STATE,1,IO) = FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = FF_G1(J,G1(I),G2(I))
               ELSE
                  WFF(J,I,STATE,1,IO) = FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = FF_G1(J,G1(I),G2(I))
               END IF
C
               IF ( STATE.EQ.1 ) THEN
                  WFFZ(J,I,STATE,1,IO) = FF_S1Z(J,G1(I),G2(I))
                  WFFZ(J,I,STATE,2,IO) = FF_G1Z(J,G1(I),G2(I))
               END IF
            END DO
         END DO
C
         IF ( IP.GE.2 ) THEN
            WRITE (NOUT1,99001) ATOM,STATE,IO
            DO I = 1,MAXG
               MU = NRM(G1(I)) + NRS(G1(I))
               MUS = NRM(G2(I)) + NRS(G2(I))
               WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,STATE,
     &                             RADMESHM(RSTEP,ATOM,LAY),
     &                             RADMESHM(RSTEP,ATOM,LAY)
     &                             *WFF(RSTEP,I,STATE,1,IO)
            END DO
            WRITE (NOUT1,99002) ATOM,STATE,IO
            DO I = 1,MAXG
               MU = NRM(G1(I)) + NRS(G1(I))
               MUS = NRM(G2(I)) + NRS(G2(I))
               WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,STATE,
     &                             RADMESHM(RSTEP,ATOM,LAY),
     &                             RADMESHM(RSTEP,ATOM,LAY)
     &                             *WFF(RSTEP,I,STATE,2,IO)
            END DO
            IF ( STATE.EQ.1 ) THEN
               WRITE (NOUT1,99003) ATOM,STATE,IO
               DO I = 1,MAXG
                  MU = NRM(G1(I)) + NRS(G1(I))
                  MUS = NRM(G2(I)) + NRS(G2(I))
                  WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,
     &                                STATE,RADMESHM(RSTEP,ATOM,LAY),
     &                                RADMESHM(RSTEP,ATOM,LAY)
     &                                *WFFZ(RSTEP,I,STATE,1,IO)
               END DO
               WRITE (NOUT1,99004) ATOM,STATE,IO
               DO I = 1,MAXG
                  MU = NRM(G1(I)) + NRS(G1(I))
                  MUS = NRM(G2(I)) + NRS(G2(I))
                  WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,
     &                                STATE,RADMESHM(RSTEP,ATOM,LAY),
     &                                RADMESHM(RSTEP,ATOM,LAY)
     &                                *WFF(RSTEP,I,STATE,2,IO)
               END DO
            END IF
         END IF
C
      ELSE
C
C        ORBITAL DEPENDENT WEIGHTS FOR FINAL STATE WAVE FUNCTIONS
C
         DO I = 1,MAXG
            DO J = 1,RSTEP
               IF ( (KAP(G1(I)).LT.-4) .AND. (KAP(G2(I)).LT.-4) ) THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).LT.-4) .AND. (KAP(G2(I)).GT.3) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).GT.3) .AND. (KAP(G2(I)).LT.-4) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).GT.3) .AND. (KAP(G2(I)).GT.3) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).LT.-3) .AND. (KAP(G2(I)).LT.-3) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).LT.-3) .AND. (KAP(G2(I)).GT.2) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).GT.2) .AND. (KAP(G2(I)).LT.-3) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE IF ( (KAP(G1(I)).GT.2) .AND. (KAP(G2(I)).GT.2) )
     &                   THEN
                  WFF(J,I,STATE,1,IO) = 1.0D0*FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = 1.0D0*FF_G1(J,G1(I),G2(I))
               ELSE
                  WFF(J,I,STATE,1,IO) = FF_S1(J,G1(I),G2(I))
                  WFF(J,I,STATE,2,IO) = FF_G1(J,G1(I),G2(I))
               END IF
            END DO
         END DO
C
         IF ( IP.GE.2 ) THEN
            WRITE (NOUT1,99001) ATOM,STATE,IO
            DO I = 1,MAXG
               MU = NRM(G1(I)) + NRS(G1(I))
               MUS = NRM(G2(I)) + NRS(G2(I))
               WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,STATE,
     &                             RADMESHM(RSTEP,ATOM,LAY),
     &                             RADMESHM(RSTEP,ATOM,LAY)
     &                             *WFF(RSTEP,I,STATE,1,IO)
            END DO
            WRITE (NOUT1,99002) ATOM,STATE,IO
            DO I = 1,MAXG
               MU = NRM(G1(I)) + NRS(G1(I))
               MUS = NRM(G2(I)) + NRS(G2(I))
               WRITE (NOUT1,99005) I,KAP(G1(I)),KAP(G2(I)),MU,MUS,STATE,
     &                             RADMESHM(RSTEP,ATOM,LAY),
     &                             RADMESHM(RSTEP,ATOM,LAY)
     &                             *WFF(RSTEP,I,STATE,2,IO)
            END DO
         END IF
C
      END IF
C
99001 FORMAT (2x,'wff-munich solution at r=721 for atom',2x,i3,2x,
     &        'state',2x,i3,2x,'type',2x,i3)
99002 FORMAT (2x,'wfg-munich solution at r=721 for atom',2x,i3,2x,
     &        'state',2x,i3,2x,'type',2x,i3)
99003 FORMAT (2x,'wffz-munich solution at r=721 for atom',2x,i3,2x,
     &        'state',2x,i3,2x,'type',2x,i3)
99004 FORMAT (2x,'wfgz-munich solution at r=721 for atom',2x,i3,2x,
     &        'state',2x,i3,2x,'type',2x,i3)
99005 FORMAT (i3,2x,2I4,2x,2F5.2,2x,i3,2x,3E14.7)
C
      END
C*==besneu.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE BESNEU(X,DKAPPA,FJ,FN,FJP,FNP,MAXL)
C
C     /****************************************************************/
C      complex*16 spherical bessel functions from l=0 to l=maxl
C      for x in the upper half plane ( im(x) > -3)
C      tampered with to calculate neumann instead of hankel functions.
C     =================================================================
C       xj(l)   = j/l(x)          regular solution: xj(0)=sin(x)/x
C       xjp(l)  = d/dx j/l(x)
C       xh1(l)  = h(1)/l(x)       irregular hankel function:
C       xh1p(l) = d/dx h(1)/l(x)            xh1(0) = j0(x) + i. y0(x)
C                                                  =(sin(x)-i.cos(x))/x
C                                                  = -i.exp(i.x)/x
C       xn(l)   = n/l(x)
C       xnp(l)  = d/dx n/l(x)
C       using complex*16 cf1, and trigonometric forms for l=0 solutions.
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:ML,CZERO,CONE,CIMAG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER LIMIT
      PARAMETER (LIMIT=20000)
      REAL*8 TM30,ACCUR
      PARAMETER (TM30=1.0D-30,ACCUR=0.5D-11)
C
C Dummy arguments
C
      COMPLEX*16 DKAPPA,X
      INTEGER MAXL
      COMPLEX*16 FJ(0:ML),FJP(0:ML),FN(0:ML),FNP(0:ML)
C
C Local variables
C
      COMPLEX*16 B,C,D,DEL,F,PL,W,XH1(:),XH1P(:),XI,XJ0
      INTEGER IFAIL,L
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE XH1,XH1P
      ALLOCATE (XH1(0:ML),XH1P(0:ML))
C
      IFAIL = -1
      IF ( CDABS(X).GE.ACCUR .AND. DIMAG(X).GE.-3.D0 ) THEN
C
         XI = CONE/X
         W = XI + XI
         PL = MAXL*XI
         F = PL + XI
         B = F + F + XI
         D = CZERO
         C = F
C
         DO L = 1,LIMIT
            D = B - D
            C = B - CONE/C
            IF ( CDABS(D).LT.TM30 ) D = TM30
            IF ( CDABS(C).LT.TM30 ) C = TM30
            D = CONE/D
            DEL = D*C
            F = F*DEL
            B = B + W
            IF ( CDABS(DEL-CONE).LT.ACCUR ) GOTO 100
         END DO
C
         IFAIL = -2
      END IF
C
      WRITE (*,99001) IFAIL,CDABS(X),DIMAG(X)
      RETURN
C
 100  CONTINUE
      FJ(MAXL) = TM30
      FJP(MAXL) = F*FJ(MAXL)
C
C     downward recursion to l=0 (n.b.  coulomb functions)
      DO L = MAXL - 1,0, - 1
         FJ(L) = PL*FJ(L+1) + FJP(L+1)
         FJP(L) = PL*FJ(L) - FJ(L+1)
         PL = PL - XI
      END DO
C
C     calculate the l=0 bessel functions
      XJ0 = XI*CDSIN(X)
      XH1(0) = CDEXP(CIMAG*X)*XI*(-CIMAG)
      XH1P(0) = XH1(0)*(CIMAG-XI)
      FN(0) = -CDCOS(X)*XI
      FNP(0) = XJ0 - FN(0)*XI
C
C     rescale xj, xjp,  converting to spherical bessels.
C     recur   xh1,xh1p             as spherical bessels.
      W = CONE/FJ(0)
      PL = XI
C
      DO L = 0,MAXL
         FJ(L) = XJ0*(W*FJ(L))
         FJP(L) = XJ0*(W*FJP(L)) - XI*FJ(L)
         IF ( L.NE.0 ) THEN
            XH1(L) = (PL-XI)*XH1(L-1) - XH1P(L-1)
            PL = PL + XI
            XH1P(L) = -PL*XH1(L) + XH1(L-1)
            FN(L) = -(XH1(L)-FJ(L))*CIMAG
            FNP(L) = -(XH1P(L)-FJP(L))*CIMAG
         END IF
      END DO
C
      DO L = 0,MAXL
         FJP(L) = FJP(L)*DKAPPA
         FNP(L) = FNP(L)*DKAPPA
      END DO
C
      IFAIL = 0
      RETURN
C
99001 FORMAT (2x,'besneu : ifail = ',i3,2(1x,e14.7))
      END
C*==indexf.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE INDEXF(LAY,ATOM,IPP,IPM,IMP,IMM,PP1,PP2,PP3,PM1,PM2,
     &                  PM3,MAXPP,MAXPM,SELRUL,IO,VLM)
C
C     /****************************************************************/
C     # purpose      : calculate the number of non-zero                *
C                      gaunt coefficients symmetry adapted for         *
C                      the lattice structure                           *
C                      and                                             *
C                      the sums over the gaunt coefficients and        *
C                      the spin up and spin down potentials: vip,vim   *
C                                                                      *
C     # calls the following subroutines and functions:                 *
C       rind                                                           *
C       blm     clegor      gaunt                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,ML,MQD,RSTEP,MAXIPP,EPS12,CZERO,
     &    CONE,NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_POTLM,ONLY:EB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY,MAXPM,MAXPP
      COMPLEX*16 IMM(MAXIPP),IMP(MAXIPP),IPM(MAXIPP),IPP(MAXIPP),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
      INTEGER PM1(MAXIPP),PM2(MAXIPP),PM3(MAXIPP),PP1(MAXIPP),
     &        PP2(MAXIPP),PP3(MAXIPP),SELRUL(MQD,MQD)
C
C Local variables
C
      REAL*8 BFIELD,BX,BXY,BY,BZ,G,IM0,IM11,IM12,IM2,IM31,IM41,IP0,IP11,
     &       IP12,IP2,IP31,IP41,KA,KH,MH,MU,NRS(:),SD,SU,TB
      REAL*8 BLM,CLEGOR,GAUNT
      INTEGER H,H1(:),H2(:),H3(:),H4(:),I,IH,IS(:),ITB,J,K,KAP(:),LRI,
     &        LRJ,MAXH,MRI,MRJ,MUD(:),NBFIELD,NBXY,NBZ,NR(:),NRL(:),
     &        NRM(:),PM,PP,REL(:,:)
      COMPLEX*16 IBY,RBX,RBZ,TMP_IMM,TMP_IMP,TMP_IPM,TMP_IPP
      EXTERNAL BLM,CLEGOR,GAUNT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE H1,H2,H3,H4,IS,NR,KAP,REL,MUD,NRL,NRM,NRS
      ALLOCATE (H1(MAXIPP),H2(MAXIPP),H3(MAXIPP),H4(MAXIPP),IS(MQD))
      ALLOCATE (NR(MQD),KAP(MQD),REL(MQD,2),MUD(MQD),NRL(MQD))
      ALLOCATE (NRM(MQD),NRS(MQD))
C
C*** End of declarations rewritten by SPAG
C
      CALL RIND(ML,MQD,KAP,MUD,NRL,NRM,NRS,NR,REL)
C
      DO I = 1,MQD
         IS(I) = INT(NRS(I)+0.5D0)
      END DO
C
C
      BX = EB(1)
      BY = EB(2)
      BZ = EB(3)
C
      RBX = DCMPLX(BX,0.D0)
      IBY = DCMPLX(0.D0,BY)
      RBZ = DCMPLX(BZ,0.D0)
C
      DO I = 1,MQD
         DO J = 1,MQD
            SELRUL(I,J) = 0
         END DO
      END DO
C
      NBXY = 0
      NBZ = 0
      NBFIELD = 0
      BXY = SQRT(BX**2+BY**2)
      BFIELD = SQRT(BXY**2+BZ**2)
      IF ( BXY.GT.EPS12 ) NBXY = 1
      IF ( ABS(BZ).GT.EPS12 ) NBZ = 1
      IF ( BFIELD.GT.EPS12 ) NBFIELD = 1
C
      H = 0
C
C     prepare fullpot selection rules for the coupling integrals
C     ==========================================================
      DO I = 1,MQD
         DO K = 1,MQD
            DO J = 1,MQD
               TB = 0.D0
               IH = 0
C              SPAG: Pls keep this condition as it is.
               IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).NE.0.D0 ) THEN
C               IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).GT.1.0D-16 ) THEN
C
                  LRI = NRL(I)
                  MRI = NRM(I)
                  LRJ = NRL(J)
                  MRJ = NRM(J)
C
C                 elements with equal spin
                  TB = ABS(BLM(LRI,-MRI,NRL(K),NRM(K),LRJ,MRJ))
C                  IF ( NRS(I).EQ.NRS(J) .AND. TB.GT.EPS12 .AND.
C     &                 NBXY.EQ.0 ) IH = 1
                  IF ( ABS(NRS(I)-NRS(J)).LE.1.0D-16 .AND. (TB-EPS12)
     &                 .GT.1.0D-16 .AND. NBXY.EQ.0 ) IH = 1
C                 magnetisation couples states with same l
C                 here the small coupling (l+-2) may be respected
C                 if used below in imp and imm
                  IF ( NBFIELD.EQ.1 .AND. LRI.EQ.LRJ ) THEN
                     ITB = 0
C                     magnetisation in xy plane
                     IF ( NBXY.EQ.1 ) THEN
                        IF ( MUD(J).EQ.(MUD(I)+2) .AND. NRS(I).GE.NRS(J)
     &                       ) THEN
                           MRI = INT(MUD(I)/2.D0+0.5D0)
                           MRJ = INT(MUD(J)/2.D0-0.5D0)
                           ITB = 1
                        ELSE IF ( MUD(J).EQ.(MUD(I)-2) .AND. NRS(I)
     &                            .LE.NRS(J) ) THEN
                           MRI = INT(MUD(I)/2.D0-0.5D0)
                           MRJ = INT(MUD(J)/2.D0+0.5D0)
                        END IF
                        ITB = 1
                     END IF
C
C                     magnetisation along z
C                     IF ( NBZ.EQ.1 .AND. NRS(J).NE.NRS(I) .AND. MUD(I)
C     &                    .EQ.MUD(J) ) THEN
                     IF ( NBZ.EQ.1 .AND. ABS(NRS(J)-NRS(I))
     &                    .GT.1.0D-16 .AND. MUD(I).EQ.MUD(J) ) THEN
                        MRJ = NRM(J) + INT(2.D0*NRS(J))
                        MRI = NRM(I)
                        ITB = 1
                     END IF
C
C                     test valid m<=l etc.
                     IF ( ITB.EQ.1 .AND. ABS(MRJ).LE.ABS(LRJ) .AND. 
     &                    ABS(MRI).LE.ABS(LRI) )
     &                    TB = ABS(BLM(LRI,-MRI,NRL(K),NRM(K),LRJ,MRJ))
                  END IF
C
                  IF ( TB.GT.EPS12 ) THEN
                     H = H + 1
                     IF ( H.GT.MAXIPP ) THEN
                        WRITE (*,99007)
                        STOP
                     END IF
                     H1(H) = K
                     H2(H) = J
                     H3(H) = I
                     H4(H) = 0
                     IF ( IH.EQ.1 ) H4(H) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      MAXH = H
C
C     calculation of the coupling integrals
C     =====================================
      SU = 0.5D0
      SD = -0.5D0
      PP = 0
      PM = 0
C
      DO H = 1,MAXH
C
         IPP(H) = CZERO
         IPM(H) = CZERO
C
         KA = DBLE(KAP(H3(H)))
         MU = DBLE(MUD(H3(H)))/2.D0
         KH = DBLE(KAP(H2(H)))
         MH = DBLE(MUD(H2(H)))/2.D0
C
         G = GAUNT(H3(H),H2(H),H1(H),NBXY,NBZ,NBFIELD)
C
         IF ( (H2(H).NE.H3(H) .AND. H4(H).EQ.0) ) THEN
            IP11 = 0.D0
            IP12 = 0.D0
         ELSE
            IP11 = CLEGOR(KA,MU,SU)*CLEGOR(KH,MH,SU)*G
            IP12 = CLEGOR(KA,MU,SD)*CLEGOR(KH,MH,SD)*G
         END IF
         IP0 = IP11 + IP12
C
         IF ( NBFIELD.NE.0 ) THEN
C             magnetic case
C             =============
            IP2 = 0.D0
            IP31 = 0.D0
            IP41 = 0.D0
C            IF ( (KH.EQ.-KA-1.D0) .AND. (MH.EQ.MU) ) THEN
            IF ( ABS(KH-(-KA-1.D0)).LE.1.0D-16 .AND. ABS(MH-MU)
     &           .LE.1.0D-16 ) THEN
               IP2 = (CLEGOR(KA,MU,SU)*CLEGOR(KH,MH,SU)-CLEGOR(KA,MU,SD)
     &               *CLEGOR(KH,MH,SD))*G
            ELSE
               IP2 = 0.D0
            END IF
C
C            IF ( (KH.EQ.KA) .AND. (MH.EQ.MU-1.D0) ) THEN
            IF ( ABS(KH-KA).LE.1.0D-16 .AND. ABS(MH-(MU-1.D0))
     &           .LE.1.0D-16 ) THEN
               IP31 = CLEGOR(KA,MU,SU)*CLEGOR(KA,MH,SD)*G
C            ELSE IF ( (KH.EQ.-KA-1.D0) .AND. (MH.EQ.MU-1.D0) ) THEN
            ELSE IF ( ABS(KH-(-KA-1.D0)).LE.1.0D-16 .AND. 
     &                ABS(MH-(MU-1.D0)).LE.1.0D-16 ) THEN
               IP31 = CLEGOR(KA,MU,SU)*CLEGOR(KH,MH,SD)*G
            ELSE
               IP31 = 0.D0
            END IF
C
C            IF ( (KH.EQ.KA) .AND. (MH.EQ.MU+1.D0) ) THEN
            IF ( ABS(KH-KA).LE.1.0D-16 .AND. ABS(MH-(MU+1.D0))
     &           .LE.1.0D-16 ) THEN
               IP41 = CLEGOR(KA,MU,SD)*CLEGOR(KA,MH,SU)*G
C            ELSE IF ( (KH.EQ.-KA-1.D0) .AND. (MH.EQ.MU+1.D0) ) THEN
            ELSE IF ( ABS(KH-(-KA-1.D0)).LE.1.0D-16 .AND. 
     &                ABS(MH-(MU+1.D0)).LE.1.0D-16 ) THEN
               IP41 = CLEGOR(KA,MU,SD)*CLEGOR(KH,MH,SU)*G
            ELSE
               IP41 = 0.D0
            END IF
C
            TMP_IPM = ((CONE-RBZ)*IP11+(CONE+RBZ)*IP12-RBZ*IP2-(RBX-IBY)
     &                *IP31-(RBX+IBY)*IP41)*DBLE(IS(H1(H)))
C
            TMP_IPP = ((CONE+RBZ)*IP11+(CONE-RBZ)*IP12+RBZ*IP2+(RBX-IBY)
     &                *IP31+(RBX+IBY)*IP41)*DBLE(1-IS(H1(H)))
C
         ELSE
C             non-magnetic case
C             =================
            TMP_IPP = DCMPLX(IP0,0.D0)*DBLE(1-IS(H1(H)))
            TMP_IPM = DCMPLX(IP0,0.D0)*DBLE(IS(H1(H)))
         END IF
C
C         IF ( NRS(H1(H)).EQ.+0.5D0 ) TMP_IPP = CZERO
C         IF ( NRS(H1(H)).EQ.-0.5D0 ) TMP_IPM = CZERO
         IF ( ABS(NRS(H1(H))-0.5D0).LE.1.0D-16 ) TMP_IPP = CZERO
         IF ( ABS(NRS(H1(H))+0.5D0).LE.1.0D-16 ) TMP_IPM = CZERO
C
         IF ( (CDABS(TMP_IPM).GT.0.0D0) .OR. (CDABS(TMP_IPP)).GT.0.0D0 )
     &        THEN
            PP = PP + 1
            IF ( PP.GT.MAXIPP ) THEN
               WRITE (*,99007)
               STOP
            END IF
            IPP(PP) = TMP_IPP
            IPM(PP) = TMP_IPM
            PP1(PP) = H1(H)
            PP2(PP) = H2(H)
            PP3(PP) = H3(H)
            SELRUL(H2(H),H3(H)) = 1
         END IF
      END DO
C
      MAXPP = PP
C
      DO H = 1,MAXH
C
         IMP(H) = CZERO
         IMM(H) = CZERO
C
         KA = DBLE(KAP(H3(H)))
         MU = DBLE(MUD(H3(H)))/2.D0
         KH = DBLE(KAP(H2(H)))
         MH = DBLE(MUD(H2(H)))/2.D0
C
         G = GAUNT(H3(H),H2(H),H1(H),NBXY,NBZ,NBFIELD)
C
         IF ( (H2(H).NE.H3(H) .AND. H4(H).EQ.0) ) THEN
            IM11 = 0.D0
            IM12 = 0.D0
         ELSE
            IM11 = CLEGOR((-KA),MU,SU)*CLEGOR((-KH),MH,SU)*G
            IM12 = CLEGOR((-KA),MU,SD)*CLEGOR((-KH),MH,SD)*G
         END IF
         IM0 = IM11 + IM12
C
         IF ( NBFIELD.NE.0 ) THEN
C             magnetic case
C             =============
            IM2 = 0.D0
            IM31 = 0.D0
            IM41 = 0.D0
C            IF ( (KH.EQ.-KA-1.D0) .AND. (MH.EQ.MU) ) THEN
            IF ( ABS(KH-(-KA-1.D0)).LE.1.0D-16 .AND. ABS(MH-MU)
     &           .LE.1.0D-16 ) THEN
               IM2 = (CLEGOR((-KA),MU,SU)*CLEGOR((-KH),MH,SU)
     &               -CLEGOR((-KA),MU,SD)*CLEGOR((-KH),MH,SD))*G
            ELSE
               IM2 = 0.D0
            END IF
C            IF ( (KH.EQ.KA) .AND. (MH.EQ.MU-1.D0) )
C     &           IM31 = CLEGOR((-KA),MU,SU)*CLEGOR((-KH),MH,SD)*G
            IF ( ABS(KH-KA).LE.1.0D-16 .AND. ABS(MH-(MU-1.D0))
     &           .LE.1.0D-16 ) IM31 = CLEGOR((-KA),MU,SU)
     &                                *CLEGOR((-KH),MH,SD)*G
C            IF ( (KH.EQ.KA) .AND. (MH.EQ.MU+1.D0) )
            IF ( ABS(KH-KA).LE.1.0D-16 .AND. ABS(MH-(MU+1.D0))
     &           .LE.1.0D-16 ) IM41 = CLEGOR((-KA),MU,SD)
     &                                *CLEGOR((-KH),MH,SU)*G
C
            TMP_IMP = ((CONE+RBZ)*IM11+(CONE-RBZ)*IM12+RBZ*IM2+(RBX-IBY)
     &                *IM31+(RBX+IBY)*IM41)*DBLE(IS(H1(H)))
C
            TMP_IMM = ((CONE-RBZ)*IM11+(CONE+RBZ)*IM12-RBZ*IM2-(RBX-IBY)
     &                *IM31-(RBX+IBY)*IM41)*DBLE(1-IS(H1(H)))
         ELSE
C             non-magnetic case
C             =================
            TMP_IMP = DCMPLX(IM0,0.D0)*DBLE(IS(H1(H)))
            TMP_IMM = DCMPLX(IM0,0.D0)*DBLE(1-IS(H1(H)))
         END IF
C
C         IF ( NRS(H1(H)).EQ.-0.5D0 ) TMP_IMP = CZERO
C         IF ( NRS(H1(H)).EQ.+0.5D0 ) TMP_IMM = CZERO
         IF ( ABS(NRS(H1(H))+0.5D0).LE.1.0D-16 ) TMP_IMP = CZERO
         IF ( ABS(NRS(H1(H))-0.5D0).LE.1.0D-16 ) TMP_IMM = CZERO
C
         IF ( (CDABS(TMP_IMP).GT.0.0D0) .OR. (CDABS(TMP_IMM)).GT.0.0D0 )
     &        THEN
            PM = PM + 1
            IF ( PM.GT.MAXIPP ) THEN
               WRITE (*,99007)
               STOP
            END IF
            IMP(PM) = TMP_IMP
            IMM(PM) = TMP_IMM
            PM1(PM) = H1(H)
            PM2(PM) = H2(H)
            PM3(PM) = H3(H)
            SELRUL(H2(H),H3(H)) = 1
         END IF
      END DO
      MAXPM = PM
C
C     write (*,900)  maxh
C     write (*,1100) maxpp
C     write (*,1200) maxpm
C
      IF ( IP.GT.4 ) THEN
         IF ( BFIELD.LT.EPS12 ) THEN
            WRITE (NOUT1,*) ' paramagnetic'
         ELSE
            WRITE (NOUT1,*) ' ferromagnetic'
         END IF
         WRITE (NOUT1,99001) MAXH
C
         WRITE (NOUT1,99003) MAXPP
         IF ( MAXPP.GT.0 ) WRITE (NOUT1,99004)
         DO PP = 1,MAXPP
            WRITE (NOUT1,99002) PP1(PP),PP2(PP),PP3(PP),IPP(PP),IPM(PP)
         END DO
C
         WRITE (NOUT1,99005) MAXPM
         IF ( MAXPM.GT.0 ) WRITE (NOUT1,99006)
         DO PM = 1,MAXPM
            WRITE (NOUT1,99002) PM1(PM),PM2(PM),PM3(PM),IMP(PM),IMM(PM)
         END DO
      END IF
      RETURN
C
99001 FORMAT (1x,'            maxh  = ',i4)
99002 FORMAT (3I3,1x,2(1x,2E14.6))
99003 FORMAT (1x,'ipp,ipm for maxpp = ',i4)
99004 FORMAT (1x,'k1',2x,'k',1x,'ks',4x,'ipp',24x,'ipm')
99005 FORMAT (1x,'imp,imm for maxpm = ',i4)
99006 FORMAT (1x,'k1',2x,'k',1x,'ks',4x,'imp',24x,'imm')
99007 FORMAT (2x,'in indexf: maxipp not big enough...aborting')
      END
C*==gaunt.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      REAL*8 FUNCTION GAUNT(I,J,K,NBXY,NBZ,NBFIELD)
C
C     /****************************************************************/
C     # purpose      : calculate gaunt coefficients                    *
C                                                                      *
C     # calls subroutine and function:                                 *
C       rind       blm                                                 *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:EPS12
      USE MOD_SPEC_RINDC,ONLY:MUD,NRL,NRM,NRS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,J,K,NBFIELD,NBXY,NBZ
C
C Local variables
C
      REAL*8 BLM
      INTEGER ITB,LRI,LRJ,MRI,MRJ
      REAL*8 TBLM
      EXTERNAL BLM
C
C*** End of declarations rewritten by SPAG
C
      GAUNT = 0.D0
      TBLM = 0.D0
C
      LRI = NRL(I)
      MRI = NRM(I)
      LRJ = NRL(J)
      MRJ = NRM(J)
      ITB = 0
C
C      IF ( NRS(I).EQ.NRS(J) .AND. NBXY.EQ.0 )
      IF ( ABS(NRS(I)-NRS(J)).LE.1.0D-16 .AND. NBXY.EQ.0 )
     &     TBLM = DBLE(((-1)**NRM(I))
     &     *BLM(LRI,-MRI,NRL(K),NRM(K),LRJ,MRJ))
C
C     magnetisation couples states with equal l
      IF ( NBFIELD.NE.0 .AND. LRI.EQ.LRJ ) THEN
         ITB = 0
C         magnetisation in xy plane
         IF ( NBXY.NE.0 ) THEN
            IF ( MUD(J).EQ.(MUD(I)+2) .AND. NRS(I).GE.NRS(J) ) THEN
               MRI = INT(MUD(I)/2.D0+0.5D0)
               MRJ = INT(MUD(J)/2.D0-0.5D0)
            ELSE IF ( MUD(J).EQ.(MUD(I)-2) .AND. NRS(I).LE.NRS(J) ) THEN
               MRI = INT(MUD(I)/2.D0-0.5D0)
               MRJ = INT(MUD(J)/2.D0+0.5D0)
            END IF
            ITB = 1
         END IF
C         magnetisation along z
C         IF ( NBZ.NE.0 .AND. NRS(J).NE.NRS(I) .AND. MUD(I).EQ.MUD(J) )
C     &        THEN
         IF ( NBZ.NE.0 .AND. ABS(NRS(J)-NRS(I)).GT.1.0D-16 .AND. 
     &        ABS(MUD(I)-MUD(J)).LE.1.0D-16 ) THEN
            MRJ = NRM(J) + INT(2.D0*NRS(J))
            MRI = NRM(I)
            ITB = 1
         END IF
C         calculate new tblm if allowed by magnetic selection rules
         IF ( ITB.EQ.1 .AND. ABS(MRJ).LE.ABS(LRJ) .AND. ABS(MRI)
     &        .LE.ABS(LRI) ) TBLM = DBLE((-1)**MRI)
     &                              *BLM(LRI,-MRI,NRL(K),NRM(K),LRJ,MRJ)
      END IF
C
      IF ( ABS(TBLM).GT.EPS12 ) GAUNT = TBLM
C
      END
C*==regularf.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE REGULARF(LAY,ATOM,CC,SS,G1,G2,MAXG,GRID,VORZ,PP1,PP2,
     &                    PP3,PM1,PM2,PM3,MAXPP,MAXPM,KAP,NRL,IPP,IPM,
     &                    IMP,IMM,EREL,CLIGHT,G1_NEW,G2_NEW,MAXG_NEW,
     &                    FF_S1,FF_G1,IO,VLM)
C     /****************************************************************/
C     # purpose      : integrates the coupled set of differential      *
C                      equations using a predictor-corrector           *
C                      integration routine                             *
C                                                                      *
C     # calls the following subroutines:                               *
C       derivatf                                                       *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEPP5,MAXIPP,MQK,NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC_H3NARRAY,ONLY:H3N
      USE MOD_SPEC_POTLM,ONLY:ALPHA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY,MAXG,MAXG_NEW,MAXPM,MAXPP
      REAL*8 CLIGHT,VORZ
      COMPLEX*16 EREL
      COMPLEX*16 CC(MQD,MQD,7),FF_G1(RSTEPP5,MQD,MQD),
     &           FF_S1(RSTEPP5,MQD,MQD),IMM(MAXIPP),IMP(MAXIPP),
     &           IPM(MAXIPP),IPP(MAXIPP),SS(MQD,MQD,7),
     &           VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
      INTEGER G1(MQK),G1_NEW(MQK),G2(MQK),G2_NEW(MQK),KAP(MQD),NRL(MQD),
     &        PM1(MAXIPP),PM2(MAXIPP),PM3(MAXIPP),PP1(MAXIPP),
     &        PP2(MAXIPP),PP3(MAXIPP)
      REAL*8 GRID(RSTEPP5)
C
C Local variables
C
      REAL*8 ALPHA1,ALPHA2,CO(7),PR(6),RGRID
      COMPLEX*16 CINT(7),DC(:,:,:),DS(:,:,:),FLOW(:),HELP,SINT(7)
      INTEGER FIX,GP1(:),GP2(:),H3Q(:),I,II,ITEST,J,NMAXG,P,Q,RR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DC,DS,GP1,GP2,H3Q,FLOW
      ALLOCATE (DC(MQD,MQD,7),DS(MQD,MQD,7),GP1(MQK),GP2(MQK))
      ALLOCATE (H3Q(MQD),FLOW(MQD))
C
      IF ( .NOT.ALLOCATED(H3N) ) ALLOCATE (H3N(MQD))
      ALPHA1 = ALPHA(LAY,ATOM)/1440.D0
      ALPHA2 = ALPHA(LAY,ATOM)/60480.D0
      HELP = CLIGHT*CDSQRT(2.D0*EREL)/(EREL+2.D0*CLIGHT**2)
C
      IF ( VORZ.GT.0.D0 ) THEN
         ITEST = 1
C
         DO I = 1,MQD
            H3N(I) = NRL(I)
            H3Q(I) = H3N(I) - INT(KAP(I)/ABS(KAP(I)))
            FLOW(I) = HELP*DBLE(KAP(I)/ABS(KAP(I)))
         END DO
C
         MAXG = 0
         NMAXG = MQD*MQD
         DO RR = 1,6
            FIX = RR
            II = 0
            DO Q = 1,MQD
               DO P = 1,MQD
                  II = II + 1
                  GP1(II) = P
                  GP2(II) = Q
               END DO
            END DO
C
            RGRID = GRID(RR)
            CALL DERIVATF(LAY,ATOM,RR,FIX,RGRID,VORZ,PP1,PP2,PP3,PM1,
     &                    PM2,PM3,MAXPP,MAXPM,IPP,IPM,IMP,IMM,G1,G2,GP1,
     &                    GP2,MAXG,NMAXG,EREL,CLIGHT,CC,SS,DC,DS,IO,VLM,
     &                    ITEST,H3N,H3Q,FLOW,FF_S1,FF_G1)
         END DO
C
         FIX = 7
         DO RR = 7,20
            RGRID = GRID(RR)
C
C         *** predictor ***
            PR(6) = 4277.D0*GRID(RR-1)
            PR(5) = 7923.D0*GRID(RR-2)
            PR(4) = 9982.D0*GRID(RR-3)
            PR(3) = 7298.D0*GRID(RR-4)
            PR(2) = 2877.D0*GRID(RR-5)
            PR(1) = 475.D0*GRID(RR-6)
C
            DO I = 1,MAXG
               DO J = 1,6
                  CINT(J) = PR(J)*DC(G1(I),G2(I),J)
                  SINT(J) = PR(J)*DS(G1(I),G2(I),J)
               END DO
               CC(G1(I),G2(I),7) = (CINT(6)-CINT(5)+CINT(4)-CINT(3)+CINT
     &                             (2)-CINT(1))
     &                             *ALPHA1 + CC(G1(I),G2(I),6)
               SS(G1(I),G2(I),7) = (SINT(6)-SINT(5)+SINT(4)-SINT(3)+SINT
     &                             (2)-SINT(1))
     &                             *ALPHA1 + SS(G1(I),G2(I),6)
            END DO
C
            CALL DERIVATF(LAY,ATOM,RR,FIX,RGRID,VORZ,PP1,PP2,PP3,PM1,
     &                    PM2,PM3,MAXPP,MAXPM,IPP,IPM,IMP,IMM,G1,G2,GP1,
     &                    GP2,MAXG,NMAXG,EREL,CLIGHT,CC,SS,DC,DS,IO,VLM,
     &                    ITEST,H3N,H3Q,FLOW,FF_S1,FF_G1)
C
C         *** corrector ***
            CO(7) = 19087.D0*GRID(RR)
            CO(6) = 65112.D0*GRID(RR-1)
            CO(5) = 46461.D0*GRID(RR-2)
            CO(4) = 37504.D0*GRID(RR-3)
            CO(3) = 20211.D0*GRID(RR-4)
            CO(2) = 6312.D0*GRID(RR-5)
            CO(1) = 863.D0*GRID(RR-6)
C
            DO I = 1,MAXG
               DO J = 1,7
                  CINT(J) = CO(J)*DC(G1(I),G2(I),J)
                  SINT(J) = CO(J)*DS(G1(I),G2(I),J)
               END DO
               CC(G1(I),G2(I),7) = (CINT(7)+CINT(6)-CINT(5)+CINT(4)-CINT
     &                             (3)+CINT(2)-CINT(1))
     &                             *ALPHA2 + CC(G1(I),G2(I),6)
               SS(G1(I),G2(I),7) = (SINT(7)+SINT(6)-SINT(5)+SINT(4)-SINT
     &                             (3)+SINT(2)-SINT(1))
     &                             *ALPHA2 + SS(G1(I),G2(I),6)
            END DO
C
            CALL DERIVATF(LAY,ATOM,RR,FIX,RGRID,VORZ,PP1,PP2,PP3,PM1,
     &                    PM2,PM3,MAXPP,MAXPM,IPP,IPM,IMP,IMM,G1,G2,GP1,
     &                    GP2,MAXG,NMAXG,EREL,CLIGHT,CC,SS,DC,DS,IO,VLM,
     &                    ITEST,H3N,H3Q,FLOW,FF_S1,FF_G1)
            DO I = 1,6
               DO J = 1,MAXG
                  CC(G1(J),G2(J),I) = CC(G1(J),G2(J),I+1)
                  SS(G1(J),G2(J),I) = SS(G1(J),G2(J),I+1)
                  DC(G1(J),G2(J),I) = DC(G1(J),G2(J),I+1)
                  DS(G1(J),G2(J),I) = DS(G1(J),G2(J),I+1)
               END DO
            END DO
         END DO
      END IF
C
      DO I = 1,MQK
         G1_NEW(I) = G1(I)
         G2_NEW(I) = G2(I)
      END DO
      MAXG_NEW = MAXG
C
      END
C*==derivatf.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE DERIVATF(LAY,ATOM,STEP,FIX,RGRID,VORZ,PP1,PP2,PP3,PM1,
     &                    PM2,PM3,MAXPP,MAXPM,IPP,IPM,IMP,IMM,G1,G2,GP1,
     &                    GP2,MAXG,NMAXG,EREL,CLIGHT,CC,SS,DC,DS,IO,VLM,
     &                    ITEST,H3N,H3Q,FLOW,FF_S1,FF_G1)
C
C     /****************************************************************/
C     # purpose      : calculates the derivatives of the coefficient   *
C                      matrices c and s for all combinations of        *
C                      (l,m,l',m') at a given point of the mesh (step) *
C                                                                      *
C     # calls the following subroutines:                               *
C       besneu                                                         *
C     /****************************************************************/
C
      USE MOD_SPEC
      USE MOD_RMESH,ONLY:NRMAX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,FIX,IO,ITEST,LAY,MAXG,MAXPM,MAXPP,NMAXG,STEP
      REAL*8 CLIGHT,RGRID,VORZ
      COMPLEX*16 EREL
      COMPLEX*16 CC(MQD,MQD,7),DC(MQD,MQD,7),DS(MQD,MQD,7),
     &           FF_G1(RSTEPP5,MQD,MQD),FF_S1(RSTEPP5,MQD,MQD),FLOW(MQD)
     &           ,IMM(MAXIPP),IMP(MAXIPP),IPM(MAXIPP),IPP(MAXIPP),
     &           SS(MQD,MQD,7),VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
      INTEGER G1(MQK),G2(MQK),GP1(MQK),GP2(MQK),H3N(MQD),H3Q(MQD),
     &        PM1(MAXIPP),PM2(MAXIPP),PM3(MAXIPP),PP1(MAXIPP),
     &        PP2(MAXIPP),PP3(MAXIPP)
C
C Local variables
C
      COMPLEX*16 ARG,FACTOR,JJ(:),JJS(:),KAPPA,KM(:,:),KP(:,:),NN(:),
     &           NNS(:),PHIL(:,:),PHIU(:,:),VLMV
      INTEGER G,I,J,LLG,LLH,LQG,LQH,P,Q
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE JJ,KM,KP,NN,JJS,NNS,PHIL,PHIU
      ALLOCATE (JJ(0:ML),KM(MQD,MQD),KP(MQD,MQD),NN(0:ML),JJS(0:ML))
      ALLOCATE (NNS(0:ML),PHIL(MQD,MQD),PHIU(MQD,MQD))
C
      DO I = 1,MQD
         DO J = 1,MQD
            DC(I,J,FIX) = CZERO
            DS(I,J,FIX) = CZERO
            FF_S1(STEP,I,J) = CZERO
            FF_G1(STEP,I,J) = CZERO
         END DO
      END DO
C
      FACTOR = -VORZ*(EREL+2.D0*CLIGHT*CLIGHT)
     &         *RGRID*RGRID*CDSQRT(2.D0*EREL)/(CLIGHT*CLIGHT)
C
      ARG = CDSQRT(2.D0*EREL)*RGRID
      KAPPA = CDSQRT(2.D0*EREL)
      CALL BESNEU(ARG,KAPPA,JJ,NN,JJS,NNS,ML)
C
      IF ( MAXG.NE.0 .AND. ITEST.NE.1 ) THEN
         NMAXG = MAXG
         DO G = 1,NMAXG
            GP1(G) = G1(G)
            GP2(G) = G2(G)
         END DO
      END IF
      DO I = 1,MQD
         DO J = 1,MQD
            KP(I,J) = CZERO
            KM(I,J) = CZERO
            PHIU(I,J) = CZERO
            PHIL(I,J) = CZERO
         END DO
      END DO
C
      DO P = 1,NMAXG
         LLH = H3N(GP1(P))
         LQH = H3Q(GP1(P))
         PHIU(GP1(P),GP2(P)) = JJ(LLH)*CC(GP1(P),GP2(P),FIX) - NN(LLH)
     &                         *SS(GP1(P),GP2(P),FIX)
         PHIL(GP1(P),GP2(P)) = FLOW(GP1(P))*JJ(LQH)
     &                         *CC(GP1(P),GP2(P),FIX) - FLOW(GP1(P))
     &                         *NN(LQH)*SS(GP1(P),GP2(P),FIX)
      END DO
C
      IF ( STEP.GT.RSTEPP5 .OR. STEP.LT.0 ) THEN
         WRITE (*,*) 'step > rstep+5 in derivatf:',STEP,RSTEPP5
         STOP
      END IF
C
      DO I = 1,MQD
         DO P = 1,MAXPP
            IF ( ITEST.EQ.1 ) THEN
               VLMV = VLM(LAY,ATOM,PP1(P),RSTEP-20,IO)
            ELSE IF ( VORZ.GT.0.D0 ) THEN
               VLMV = VLM(LAY,ATOM,PP1(P),STEP,IO)
            ELSE IF ( STEP.LT.7 ) THEN
               VLMV = VLM(LAY,ATOM,PP1(P),RSTEP,IO)
            ELSE
               VLMV = VLM(LAY,ATOM,PP1(P),RSTEP-STEP+6,IO)
            END IF
C
            KP(PP3(P),I) = KP(PP3(P),I) + (IPM(P)+IPP(P))
     &                     *VLMV*PHIU(PP2(P),I)*CHALF
         END DO
      END DO
C
      DO I = 1,MQD
         DO P = 1,MAXPM
            IF ( ITEST.EQ.1 ) THEN
               VLMV = VLM(LAY,ATOM,PM1(P),RSTEP-20,IO)
            ELSE IF ( VORZ.GT.0.D0 ) THEN
               VLMV = VLM(LAY,ATOM,PM1(P),STEP,IO)
            ELSE IF ( STEP.LT.7 ) THEN
               VLMV = VLM(LAY,ATOM,PM1(P),RSTEP,IO)
            ELSE
               VLMV = VLM(LAY,ATOM,PM1(P),RSTEP-STEP+6,IO)
            END IF
C
            KM(PM3(P),I) = KM(PM3(P),I) + (IMP(P)+IMM(P))
     &                     *VLMV*PHIL(PM2(P),I)*CHALF
         END DO
      END DO
C
C     write (*,*) 'IPARA',ipara
C     *** test calculation ***
      IF ( ITEST.EQ.1 ) THEN
         G = 0
         DO P = 1,MQD
            DO Q = 1,MQD
C                if (fpselrul(p,q).ne.0) then
C                   if (cdabs(kp(p,q)) .gt. 0.0d0 .or.
C    1                 cdabs(km(p,q)) .gt. 0.0d0 .or.
C    2                 cdabs(cc(p,q,fix)) .ne. 0.d0 .or.
C    3                 cdabs(ss(p,q,fix)) .ne. 0.d0) then
C                      g=g+1
C                      g1(g)=p
C                      g2(g)=q
C                   end if
C                end if
C
C              SPAG: Pls keep this condition as it is.
               IF ( CDABS(KP(P,Q)).GT.0.0D0 .OR. CDABS(KM(P,Q))
     &              .GT.0.0D0 .OR. CDABS(CC(P,Q,FIX)).NE.0.D0 .OR. 
     &              CDABS(SS(P,Q,FIX)).NE.0.D0 ) THEN  !problem
C               IF ( CDABS(KP(P,Q)).GT.1.0D-16 .OR. CDABS(KM(P,Q))
C     &              .GT.1.0D-16 .OR. CDABS(CC(P,Q,FIX)).GE.1.0D-16 .OR.
C     &              CDABS(SS(P,Q,FIX)).GE.1.0D-16 ) THEN
                  G = G + 1
                  G1(G) = P
                  G2(G) = Q
               END IF
            END DO
         END DO
         MAXG = G
C         write (*,*) 'maxg',g, maxwf
C         IF ( MAXG.GT.MAXWF ) THEN
C            WRITE (*,99001)
C            WRITE (*,99002) G,MAXWF
C            STOP
C         END IF
      END IF
C
      DO G = 1,NMAXG
         LLG = H3N(GP1(G))
         LQG = H3Q(GP1(G))
         DC(GP1(G),GP2(G),FIX) = FACTOR*(NN(LLG)*KP(GP1(G),GP2(G))+FLOW(
     &                           GP1(G))*NN(LQG)*KM(GP1(G),GP2(G)))
         DS(GP1(G),GP2(G),FIX) = FACTOR*(JJ(LLG)*KP(GP1(G),GP2(G))+FLOW(
     &                           GP1(G))*JJ(LQG)*KM(GP1(G),GP2(G)))
      END DO
C
C     store phasefunctions:
C     =====================
      IF ( ITEST.NE.1 ) THEN
         DO I = 1,MAXG
            FF_S1(STEP,G1(I),G2(I)) = PHIU(G1(I),G2(I))
            FF_G1(STEP,G1(I),G2(I)) = PHIL(G1(I),G2(I))
         END DO
      END IF
C
      END
C*==scatterf.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE SCATTERF(LAY,ATOM,IO,ISTATE,FF_S1M,FF_G1M,FFI_S1M,
     &                    FFI_G1M,JE,FF_S1MZ,FF_G1MZ,G1_KKR,G2_KKR,
     &                    MAXG_KKR,TSST,MSST,P)
C
C
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,ISTATE,JE,LAY,MAXG_KKR
      COMPLEX*16 P
      COMPLEX*16 FFI_G1M(NRMAX,MQD,MQD),FFI_S1M(NRMAX,MQD,MQD),
     &           FF_G1M(NRMAX,MQD,MQD),FF_G1MZ(NRMAX,MQD,MQD),
     &           FF_S1M(NRMAX,MQD,MQD),FF_S1MZ(NRMAX,MQD,MQD),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      INTEGER G1_KKR(MQK),G2_KKR(MQK)
C
C Local variables
C
      COMPLEX*16 HF(:,:,:),HG(:,:,:),RF(:,:,:),RF1(:,:,:),RG(:,:,:),
     &           RG1(:,:,:)
      INTEGER I,IR,J,L
      REAL*8 NZERO(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE HF,HG,RF,RG,RF1,RG1,NZERO
      ALLOCATE (HF(NRMAX,MQD,MQD),HG(NRMAX,MQD,MQD))
      ALLOCATE (RF(NRMAX,MQD,MQD),RG(NRMAX,MQD,MQD))
      ALLOCATE (RF1(NRMAX,MQD,MQD),RG1(NRMAX,MQD,MQD))
      ALLOCATE (NZERO(MQD,MQD))
C
      DO I = 1,MQD
         DO J = 1,MQD
            DO L = 1,NRMAX
               RG(L,I,J) = CZERO
               RF(L,I,J) = CZERO
               RG1(L,I,J) = CZERO
               RF1(L,I,J) = CZERO
               HG(L,I,J) = CZERO
               HF(L,I,J) = CZERO
            END DO
         END DO
      END DO
C
      CALL SPEC_SSITE_DRIVE(P,TSST,RG,RF,HG,HF,ISTATE,LAY,ATOM,IO,JE,
     &                      RG1,RF1,MSST,NZERO)
C
      CALL TRANSPHOW(RG,RF,HG,HF,MQD,ML,FF_S1M,FF_G1M,FFI_S1M,FFI_G1M,
     &               RG1,RF1,FF_S1MZ,FF_G1MZ,NZERO)
C
C      NZERO = 0.0D0
C      DO I = 1,MQD
C         DO J = 1,MQD
C            DO IR = 30,100
C               NZERO(J,I) = NZERO(J,I) + ABS(FF_S1M(IR,J,I))
C            END DO
C         END DO
C      END DO
C
      IR = 0
      DO J = 1,MQD
         DO I = 1,MQD
            IF ( NZERO(J,I).GT.1D-16 ) THEN
               IR = IR + 1
               G1_KKR(IR) = J
               G2_KKR(IR) = I
            END IF
         END DO
      END DO
      MAXG_KKR = IR
      END
C*==initref.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE INITREF(LAY,ATOM,CC,SS,GRID,VORZ,RMAX)
C     /****************************************************************/
C     # purpose      : initializes the matrices c and s                *
C                      for the regular solution near the origin :      *
C                      c = unit matrix                                 *
C                      s = zero matrix                                 *
C                      initializes the integration direction           *
C                      for the regular solution:                       *
C                      0 -> rbs                                        *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,RSTEP,RSTEPP5,CZERO,CONE
      USE MOD_SPEC_POTLM,ONLY:MESH,ALPHA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,LAY,RMAX
      REAL*8 VORZ
      COMPLEX*16 CC(MQD,MQD,7),SS(MQD,MQD,7)
      REAL*8 GRID(RSTEPP5)
C
C Local variables
C
      INTEGER I,J,K
C
C*** End of declarations rewritten by SPAG
C
      VORZ = 1.D0
C
      DO I = 1,MQD
         DO J = 1,MQD
            DO K = 1,7
               CC(I,J,K) = CZERO
               SS(I,J,K) = CZERO
            END DO
         END DO
      END DO
C
      DO I = 1,MQD
         DO K = 1,6
            CC(I,I,K) = CONE
         END DO
      END DO
C
      DO J = 1,RMAX
         GRID(J) = MESH(LAY,ATOM,RSTEP)
     &             *DEXP(DBLE(J-RSTEP)*ALPHA(LAY,ATOM))
      END DO
C
      END
C*==magfpselrul.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE MAGFPSELRUL(LAY,ATOM,USEEULER,CORESELRUL,FPSELRUL,IO,
     &                       VLM)
C     /****************************************************************/
C     # purpose      : determine the non-zero coupling coefficients    *
C                                                                      *
C     # calls the following subroutines and functions:                 *
C       rind                                                           *
C       blm     clegor      gaunt                                      *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:LAYSM,NATLM,MQD,RSTEP,EPS12,NTPHOMAX
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      USE MOD_SPEC_POTLM,ONLY:EB,BB
      USE MOD_SPEC_RINDC,ONLY:MUD,NRL,NRM,NRS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY,USEEULER
      INTEGER CORESELRUL(MQD,MQD),FPSELRUL(MQD,MQD)
      COMPLEX*16 VLM(LAYSM,NATLM,MQD,NRMAX,NTPHOMAX)
C
C Local variables
C
      REAL*8 BFIELD,BX,BXY,BY,BZ,TB
      REAL*8 BLM
      INTEGER I,ICOUNT,ITB,J,K,LRI,LRJ,LRK,MRI,MRJ,MRK,NBFIELD,NBXY,NBZ
      EXTERNAL BLM
C
C*** End of declarations rewritten by SPAG
C
C     preset selection rules
C         diagonal elements are always present
C         clear off diagonal elements
      DO I = 1,MQD
         DO J = 1,MQD
            FPSELRUL(I,J) = 0
            CORESELRUL(I,J) = 0
         END DO
         FPSELRUL(I,I) = 1
      END DO
C
      NBXY = 0
      NBZ = 0
      NBFIELD = 0
      BX = EB(1)
      BY = EB(2)
      BZ = EB(3)
      BXY = SQRT(BX**2+BY**2)
      BFIELD = SQRT(BXY**2+BZ**2)
      IF ( BXY.GT.EPS12 ) NBXY = 1
      IF ( ABS(BZ).GT.EPS12 ) NBZ = 1
      IF ( BFIELD.GT.EPS12 ) NBFIELD = 1
      ICOUNT = 0
C
C     prepare fullpot selection rules for the coupling integrals
C     ==========================================================
      DO K = 1,MQD
C         loop on potential components (only y00 up/down if not fullpot)
C         IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).NE.0.D0 ) THEN
         IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).GT.1.0D-16 ) THEN
C
            DO I = 1,MQD
               DO J = 1,MQD
                  TB = 0.D0
C
                  LRK = NRL(K)
                  MRK = NRM(K)
                  LRI = NRL(I)
                  MRI = NRM(I)
                  LRJ = NRL(J)
                  MRJ = NRM(J)
C
C                 elements with equal spin
C                  IF ( NRS(I).EQ.NRS(J) .AND. NBXY.EQ.0 )
C     &                 TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
                  IF ( ABS(NRS(I)-NRS(J)).GT.1.0D-16 .AND. NBXY.EQ.0 )
     &                 TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
C                 magnetisation couples states with same l
C                 here the small coupling (l+-2) may be respected
C                 if used below in imp and imm
                  IF ( NBFIELD.EQ.1 .AND. LRI.EQ.LRJ ) THEN
                     ITB = 0
C                     magnetisation in xy plane
                     IF ( NBXY.EQ.1 ) THEN
                        IF ( MUD(J).EQ.(MUD(I)+2) .AND. NRS(I).GE.NRS(J)
     &                       ) THEN
                           MRI = INT(MUD(I)/2.D0+0.5D0)
                           MRJ = INT(MUD(J)/2.D0-0.5D0)
                           ITB = 1
                        ELSE IF ( MUD(J).EQ.(MUD(I)-2) .AND. NRS(I)
     &                            .LE.NRS(J) ) THEN
                           MRI = INT(MUD(I)/2.D0-0.5D0)
                           MRJ = INT(MUD(J)/2.D0+0.5D0)
                        END IF
                        ITB = 1
                     END IF
C
C                     magnetisation along z
C                     IF ( NBZ.EQ.1 .AND. NRS(J).NE.NRS(I) .AND. MUD(I)
C     &                    .EQ.MUD(J) ) THEN
                     IF ( NBZ.EQ.1 .AND. ABS(NRS(J)-NRS(I))
     &                    .GT.1.0D-16 .AND. ABS(MUD(I)-MUD(J))
     &                    .LE.1.0D-16 ) THEN
                        MRJ = NRM(J) + INT(2.D0*NRS(J))
                        MRI = NRM(I)
                        ITB = 1
                     END IF
C
C                     test valid m<=l etc.
                     IF ( ITB.EQ.1 .AND. ABS(MRJ).LE.ABS(LRJ) .AND. 
     &                    ABS(MRI).LE.ABS(LRI) )
     &                    TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
                  END IF
C
                  IF ( TB.GT.EPS12 ) THEN
                     FPSELRUL(I,J) = 1
                     ICOUNT = ICOUNT + 1
                  END IF
               END DO
            END DO
         END IF
      END DO
      IF ( IP.GT.3 ) WRITE (NOUT1,*) 'number of magnetic couplings',
     &                               ICOUNT/2
C
      IF ( USEEULER.EQ.0 ) THEN
         DO I = 1,MQD
            DO J = 1,MQD
               CORESELRUL(I,J) = FPSELRUL(I,J)
            END DO
         END DO
      ELSE
         NBXY = 0
         NBZ = 0
         NBFIELD = 0
         BX = BB(1)
         BY = BB(2)
         BZ = BB(3)
         BXY = SQRT(BX**2+BY**2)
         BFIELD = SQRT(BXY**2+BZ**2)
         IF ( BXY.GT.EPS12 ) NBXY = 1
         IF ( ABS(BZ).GT.EPS12 ) NBZ = 1
         IF ( BFIELD.GT.EPS12 ) NBFIELD = 1
C
C         prepare euler selection rules for the coupling integrals
C         ========================================================
         DO K = 1,MQD
C             loop on potential components (only y00 up/down if not fullpot)
C            IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).NE.0.D0 ) THEN
            IF ( CDABS(VLM(LAY,ATOM,K,RSTEP-20,IO)).GT.1.0D-16 ) THEN
C
               DO I = 1,MQD
                  DO J = 1,MQD
                     TB = 0.D0
C
                     LRK = NRL(K)
                     MRK = NRM(K)
                     LRI = NRL(I)
                     MRI = NRM(I)
                     LRJ = NRL(J)
                     MRJ = NRM(J)
C
C                     elements with equal spin
C                     IF ( NRS(I).EQ.NRS(J) .AND. NBXY.EQ.0 )
C     &                    TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
                     IF ( ABS(NRS(I)-NRS(J)).LE.1.0D-16 .AND. 
     &                    NBXY.EQ.0 )
     &                    TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
C                     magnetisation couples states with same l
C                     here the small coupling (l+-2) may be respected
C                     if used below in imp and imm
                     IF ( NBFIELD.EQ.1 .AND. LRI.EQ.LRJ ) THEN
                        ITB = 0
C                         magnetisation in xy plane
                        IF ( NBXY.EQ.1 ) THEN
                           IF ( MUD(J).EQ.(MUD(I)+2) .AND. NRS(I)
     &                          .GE.NRS(J) ) THEN
                              MRI = INT(MUD(I)/2.D0+0.5D0)
                              MRJ = INT(MUD(J)/2.D0-0.5D0)
                              ITB = 1
                           ELSE IF ( MUD(J).EQ.(MUD(I)-2) .AND. NRS(I)
     &                               .LE.NRS(J) ) THEN
                              MRI = INT(MUD(I)/2.D0-0.5D0)
                              MRJ = INT(MUD(J)/2.D0+0.5D0)
                           END IF
                           ITB = 1
                        END IF
C
C                         magnetisation along z
C                        IF ( NBZ.EQ.1 .AND. NRS(J).NE.NRS(I) .AND.
C     &                       MUD(I).EQ.MUD(J) ) THEN
                        IF ( NBZ.EQ.1 .AND. ABS(NRS(J)-NRS(I))
     &                       .GT.1.0D-16 .AND. ABS(MUD(I)-MUD(J))
     &                       .LE.1.0D-16 ) THEN
                           MRJ = NRM(J) + INT(2.D0*NRS(J))
                           MRI = NRM(I)
                           ITB = 1
                        END IF
C
C                         test valid m<=l etc.
                        IF ( ITB.EQ.1 .AND. ABS(MRJ).LE.ABS(LRJ) .AND. 
     &                       ABS(MRI).LE.ABS(LRI) )
     &                       TB = ABS(BLM(LRI,-MRI,LRK,MRK,LRJ,MRJ))
                     END IF
C
                     IF ( TB.GT.EPS12 ) THEN
                        CORESELRUL(I,J) = 1
                        CORESELRUL(J,I) = 1
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END IF
C
      END
C*==print_wavem.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE PRINT_WAVEM(FF_S1M,FF_G1M,G1,G2,MAXG,LAY,ATOM,STATE,IO)
C
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC
      USE MOD_SPEC_MESH,ONLY:RADMESHM
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY,MAXG,STATE
      COMPLEX*16 FF_G1M(NRMAX,MQD,MQD),FF_S1M(NRMAX,MQD,MQD)
      INTEGER G1(MQK),G2(MQK)
C
C Local variables
C
      INTEGER I,ICOUNT
C
C*** End of declarations rewritten by SPAG
C
      IF ( IP.GE.2 ) THEN
         WRITE (NOUT1,99001) ATOM,LAY,STATE,IO
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99003) ICOUNT,G1(I),G2(I),
     &                          RADMESHM(RSTEP,ATOM,LAY),
     &                          RADMESHM(RSTEP,ATOM,LAY)
     &                          *FF_S1M(RSTEP,G1(I),G2(I))
         END DO
         WRITE (NOUT1,99002) ATOM,LAY,STATE,IO
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99003) ICOUNT,G1(I),G2(I),
     &                          RADMESHM(RSTEP,ATOM,LAY),
     &                          RADMESHM(RSTEP,ATOM,LAY)
     &                          *FF_G1M(RSTEP,G1(I),G2(I))
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2x,'normalized munich-wave ff_s1(721) for atom',2x,i3,2x,
     &        'of layer',2x,i3,2x,'state',2x,i3,2x,'type',2x,i3)
99002 FORMAT (2x,'normalized munich-wave ff_g1(721) for atom',2x,i3,2x,
     &        'of layer',2x,i3,2x,'state',2x,i3,2x,'type',2x,i3)
99003 FORMAT (2x,i4,4x,i3,2x,i3,2x,e14.7,2x,2E14.7)
      END
C*==print_irrwave.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE PRINT_IRRWAVE(FF_S1M,FF_G1M,G1,G2,MAXG,LAY,ATOM,STATE,
     &                         IO)
C
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_SPEC,ONLY:MQD,RSTEP,MQK
      USE MOD_SPEC_MESH,ONLY:RADMESHM
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,IO,LAY,MAXG,STATE
      COMPLEX*16 FF_G1M(NRMAX,MQD,MQD),FF_S1M(NRMAX,MQD,MQD)
      INTEGER G1(MQK),G2(MQK)
C
C Local variables
C
      INTEGER I,ICOUNT
C
C*** End of declarations rewritten by SPAG
C
      IF ( IP.GE.2 ) THEN
         WRITE (NOUT1,99001) ATOM,LAY,STATE,IO
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99003) ICOUNT,G1(I),G2(I),
     &                          RADMESHM(RSTEP-100,ATOM,LAY),
     &                          RADMESHM(RSTEP-100,ATOM,LAY)
     &                          *FF_S1M(RSTEP-100,G1(I),G2(I))
         END DO
         WRITE (NOUT1,99002) ATOM,LAY,STATE,IO
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99003) ICOUNT,G1(I),G2(I),
     &                          RADMESHM(RSTEP-100,ATOM,LAY),
     &                          RADMESHM(RSTEP-100,ATOM,LAY)
     &                          *FF_G1M(RSTEP-100,G1(I),G2(I))
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2x,'normalized munich-irrwave ff_s1(721) for atom',2x,i3,
     &        2x,'of layer',2x,i3,2x,'state',2x,i3,2x,'type',2x,i3)
99002 FORMAT (2x,'normalized munich-irrwave ff_g1(721) for atom',2x,i3,
     &        2x,'of layer',2x,i3,2x,'state',2x,i3,2x,'type',2x,i3)
99003 FORMAT (2x,i4,4x,i3,2x,i3,2x,e14.7,2x,2E14.7)
      END
C*==rot_wave.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE ROT_WAVE(G1,G2,G1_NEW,G2_NEW,LAY,ATOM,H1,H2,H3,H4,ROT,
     &                    ROTM1,MAXG,MAXG_NEW,RMAX,FF_S1,FF_G1)
C     /****************************************************************/
C     # purpose      : rotate the wavefunctions according to b-field   *
C                      and normalize the wave functions                *
C                                                                      *
C     # calls the following subroutines:                               *
C       mult                                                           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,RSTEP,RSTEPP5,MQK,CZERO
      USE MOD_SPEC_OUTPUT,ONLY:IP,NOUT1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER ATOM,LAY,MAXG,MAXG_NEW,RMAX
      COMPLEX*16 FF_G1(RSTEPP5,MQD,MQD),FF_S1(RSTEPP5,MQD,MQD),
     &           H1(MQD,MQD),H2(MQD,MQD),H3(MQD,MQD),H4(MQD,MQD),
     &           ROT(MQD,MQD),ROTM1(MQD,MQD)
      INTEGER G1(MQK),G1_NEW(MQK),G2(MQK),G2_NEW(MQK)
C
C Local variables
C
      INTEGER I,ICOUNT,J,R
C
C*** End of declarations rewritten by SPAG
C
C     use rot and rotm1 to rotate wavefunctions according to
C     magnetic field direction; if x(r) is the matrix describing
C     the wavefunctions in z-direction, than the rotated matrix
C     is given by: x'(r)= rotm1 * x(r) * rot
C     ==========================================================
C
      IF ( IP.GT.4 ) THEN
         WRITE (NOUT1,99001) ATOM,LAY
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99005) ICOUNT,G1(I),G2(I),
     &                          FF_S1(RSTEP,G1(I),G2(I)),
     &                          FF_S1(RSTEP,G2(I),G1(I))
         END DO
         WRITE (NOUT1,99002) ATOM,LAY
         ICOUNT = 0
         DO I = 1,MAXG
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99005) ICOUNT,G1(I),G2(I),
     &                          FF_G1(RSTEP,G1(I),G2(I)),
     &                          FF_G1(RSTEP,G2(I),G1(I))
         END DO
      END IF
C
      DO R = 1,RMAX
         DO I = 1,MQD
            DO J = 1,MQD
               H1(I,J) = CZERO
               H2(I,J) = CZERO
            END DO
         END DO
         DO I = 1,MAXG
            H1(G1(I),G2(I)) = FF_S1(R,G1(I),G2(I))
            H2(G1(I),G2(I)) = FF_G1(R,G1(I),G2(I))
            FF_S1(R,G1(I),G2(I)) = CZERO
            FF_G1(R,G1(I),G2(I)) = CZERO
         END DO
C
C         build rotm1 * x * rot:
C         ======================
         CALL MULT(H1,ROT,H4,MQD)
         CALL MULT(ROTM1,H4,H3,MQD)
         CALL MULT(H2,ROT,H1,MQD)
         CALL MULT(ROTM1,H1,H4,MQD)
C
         DO I = 1,MAXG_NEW
            FF_S1(R,G1_NEW(I),G2_NEW(I)) = H3(G1_NEW(I),G2_NEW(I))
            FF_G1(R,G1_NEW(I),G2_NEW(I)) = H4(G1_NEW(I),G2_NEW(I))
         END DO
      END DO
C
      IF ( IP.GT.4 ) THEN
         WRITE (NOUT1,99003) ATOM,LAY
         ICOUNT = 0
         DO I = 1,MAXG_NEW
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99005) ICOUNT,G1_NEW(I),G2_NEW(I),
     &                          FF_S1(RSTEP,G1_NEW(I),G2_NEW(I)),
     &                          FF_S1(RSTEP,G2_NEW(I),G1_NEW(I))
         END DO
         WRITE (NOUT1,99004) ATOM,LAY
         ICOUNT = 0
         DO I = 1,MAXG_NEW
            ICOUNT = ICOUNT + 1
            WRITE (NOUT1,99005) ICOUNT,G1_NEW(I),G2_NEW(I),
     &                          FF_G1(RSTEP,G1_NEW(I),G2_NEW(I)),
     &                          FF_G1(RSTEP,G2_NEW(I),G1_NEW(I))
         END DO
      END IF
C
      RETURN
C
99001 FORMAT (2x,'original wave ff_s1(rstep) for atom:',i3,
     &        '   of layer:',i3)
99002 FORMAT (2x,'original wave ff_g1(rstep) for atom:',i3,
     &        '   of layer:',i3)
99003 FORMAT (2x,'rotated wave ff_s1(rstep) for atom:',i3,
     &        '   of layer:',i3)
99004 FORMAT (2x,'rotated wave ff_g1(rstep) for atom:',i3,
     &        '   of layer:',i3)
99005 FORMAT (2x,i4,4x,i3,2x,i3,2x,2E14.7,2x,2E14.7)
      END
C*==turn_cs.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE TURN_CS(MAXG_NEW,CC,SS,G1_NEW,G2_NEW,ROT,ROTM1,IPARA,
     &                   CORESELRUL)
C     /****************************************************************/
C     # purpose      : rotate c, s coefficients according to b-field   *
C                                                                      *
C     # calls the following subroutine:                                *
C       mult                                                           *
C     /****************************************************************/
C
      USE MOD_SPEC,ONLY:MQD,MQK,EPS12,CZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IPARA,MAXG_NEW
      COMPLEX*16 CC(MQD,MQD,7),ROT(MQD,MQD),ROTM1(MQD,MQD),SS(MQD,MQD,7)
      INTEGER CORESELRUL(MQD,MQD),G1_NEW(MQK),G2_NEW(MQK)
C
C Local variables
C
      INTEGER I,J
      COMPLEX*16 TMP_CC(:,:),TMP_CS(:,:),TMP_SS(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TMP_CC,TMP_CS,TMP_SS
      ALLOCATE (TMP_CC(MQD,MQD),TMP_CS(MQD,MQD),TMP_SS(MQD,MQD))
C
      DO I = 1,MQD
         DO J = 1,MQD
            TMP_CC(I,J) = CC(I,J,7)
            TMP_SS(I,J) = SS(I,J,7)
C
            TMP_CS(I,J) = CZERO
         END DO
      END DO
      CALL MULT(TMP_CC,ROT,TMP_CS,MQD)
      CALL MULT(ROTM1,TMP_CS,TMP_CC,MQD)
C
      DO I = 1,MQD
         DO J = 1,MQD
            TMP_CS(I,J) = CZERO
         END DO
      END DO
      CALL MULT(TMP_SS,ROT,TMP_CS,MQD)
      CALL MULT(ROTM1,TMP_CS,TMP_SS,MQD)
C
      DO I = 1,MQK
         G1_NEW(I) = 0
         G2_NEW(I) = 0
      END DO
      MAXG_NEW = 0
      DO I = 1,MQD
         DO J = 1,MQD
            CC(I,J,7) = CZERO
            SS(I,J,7) = CZERO
C
C             correct for numerical accuracy
            IF ( IPARA.NE.0 ) THEN
C                 ferromagnetic case for useuler=1
C                  if(coreselrul(i,j).eq.1) then
               IF ( CDABS(TMP_CC(I,J)).GT.EPS12 .OR. CDABS(TMP_SS(I,J))
     &              .GT.EPS12 ) THEN
                  MAXG_NEW = MAXG_NEW + 1
                  G1_NEW(MAXG_NEW) = I
                  G2_NEW(MAXG_NEW) = J
                  IF ( CORESELRUL(I,J).EQ.1 ) THEN
                     CC(I,J,7) = TMP_CC(I,J)
                     SS(I,J,7) = TMP_SS(I,J)
                  END IF
               END IF
C                  end if
C                 no selection rules applied
            ELSE IF ( CDABS(TMP_CC(I,J)).GT.EPS12 .OR. 
     &                CDABS(TMP_SS(I,J)).GT.EPS12 ) THEN
               MAXG_NEW = MAXG_NEW + 1
               G1_NEW(MAXG_NEW) = I
               G2_NEW(MAXG_NEW) = J
               CC(I,J,7) = TMP_CC(I,J)
               SS(I,J,7) = TMP_SS(I,J)
            END IF
         END DO
      END DO
C
      END
C*==tstreufull.f    processed by SPAG 6.70Rc at 08:28 on 19 Apr 2017
      SUBROUTINE TSTREUFULL(TEMP,DTEMP,E1,MASS,IREL,III,IAN,IT,IBLOCH,
     &                      TSEL,TSOL,TSEH,TSOH)
C**************************************************************
C subroutines and functions called from this routine          *
C                                                             *
C hank1     hank2     blm                                     *
C**************************************************************
C
      USE MOD_SPEC,ONLY:LAYSM,ML,MLQ,MLP,MLQNAT,PI,CZERO,CONE,CIMAG
      USE MOD_SPEC_POSI,ONLY:PTKE1,PTKO1
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 DTEMP,MASS,TEMP
      COMPLEX*16 E1
      INTEGER IAN,IBLOCH,III,IREL,IT
      COMPLEX*16 TSEH(MLQNAT,MLQNAT,LAYSM),TSEL(MLQNAT,MLQNAT,LAYSM),
     &           TSOH(MLQNAT,MLQNAT,LAYSM),TSOL(MLQNAT,MLQNAT,LAYSM)
C
C Local variables
C
      COMPLEX*16 AKQU,BLL1(2*ML,2*ML),EJOULE,FAKT1,FAKT2,TERM1(2*ML),
     &           TERM2(2*ML),TL(ML),TLM(ML),TLP(ML),TTLM(ML),TTLP(ML),
     &           WIG3J1,WIG3J2,XH1(MLP),XJ(MLP),XN(MLP),Z
      REAL*8 BLM
      INTEGER I,IAP,L,L1,L2,MAXL,MLQ1
      REAL*8 KBOLTZ,R,RINT,RSUM,U
      EXTERNAL BLM
C
C*** End of declarations rewritten by SPAG
C
      TL(:) = CZERO
      MAXL = ML - 1
      IF ( MAXL.GT.15 ) STOP 'ERROR TSTREU MAXL'
C
      MLQ1 = MLQ/2 + 2
C
      IF ( III.EQ.1 ) THEN
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TLP(PTKE1(I)) = TSEL(IAP,IAP,IT)
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TLM(PTKO1(I)) = TSOL(IAP,IAP,IT)
         END DO
         TLM(1) = CZERO
      ELSE IF ( III.EQ.2 ) THEN
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TLP(PTKE1(I)) = TSEH(IAP,IAP,IT)
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TLM(PTKO1(I)) = TSOH(IAP,IAP,IT)
         END DO
         TLM(1) = CZERO
      END IF
C
      EJOULE = E1*1.60219D+00
      KBOLTZ = 1.38062D+00
      R = TEMP/DTEMP
C
      IF ( R.LT.1.D-5 ) THEN
         AKQU = 0.3D+05*EJOULE*(1.64493D0*R*R+0.25D0)
     &          /(KBOLTZ*DTEMP*MASS)
         GOTO 100
      END IF
      IF ( 1.D0/R.GT.10.D0 ) THEN
         AKQU = 0.3D+05*EJOULE*(1.64493D0*R*R+0.25D0)
     &          /(KBOLTZ*DTEMP*MASS)
      ELSE
         U = DTEMP/TEMP
         RSUM = 1.0D0
         DO I = 1,100
            RSUM = RSUM + EXP((-REAL(I)*U))/(REAL(I+1)*REAL(I+1))
         END DO
         RINT = 1.64493D0 + U*LOG(1.D0-EXP((-U))) - EXP((-U))*RSUM
         AKQU = 0.3E+05*EJOULE*(RINT*R*R+0.25D0)/(KBOLTZ*DTEMP*MASS)
      END IF
C
C
 100  CONTINUE
      IF ( IBLOCH.EQ.1 ) AKQU = 0.D0
      AKQU = 0.D0
      Z = -2.D0*CIMAG*AKQU
      IF ( CDABS(Z).GT.0.1D0 ) THEN
         CALL HANK1(Z,XJ,XN,XH1)
      ELSE
         CALL HANK2(Z,XJ)
      END IF
      DO L = 0,MAXL
         DO L1 = 0,MAXL
            BLL1(L+1,L1+1) = CDEXP((-2.D0*AKQU))*REAL((2*L+1)*(2*L1+1))
     &                       *CIMAG**L1*XJ(L1+1)
         END DO
      END DO
      DO L = 0,MAXL
         TERM1(L+1) = CZERO
         TERM2(L+1) = CZERO
         DO L1 = 0,MAXL
            DO L2 = 0,MAXL
               WIG3J1 = DCMPLX(DSQRT(4.0*PI)*BLM(L1,0,L2,0,L,0)
     &                  /DSQRT(1.D0*((2*L1+1)*(2*L2+1)*(2*L+1))))
               FAKT1 = CONE
               IF ( L.NE.0 .AND. L2.NE.0 ) THEN
                  WIG3J2 = DCMPLX(DSQRT(4.0*PI)*BLM(L1,0,L2,1,L,-1)
     &                     /DSQRT(1.D0*((2*L1+1)*(2*L2+1)*(2*L+1))))
                  FAKT2 = DCMPLX(DSQRT(1.D0*(L2*(L2+1)))
     &                    /DSQRT(1.D0*(L*(L+1))))
               END IF
               IF ( IREL.EQ.2 ) THEN
                  TERM1(L+1) = TERM1(L+1) + BLL1(L+1,L1+1)
     &                         *WIG3J1*FAKT1*2.D0*(1.D0*(L2+1)*TLP(L2+1)
     &                         +1.D0*(L2)*TLM(L2+1))
               ELSE
                  TERM1(L+1) = TERM1(L+1) + BLL1(L+1,L1+1)
     &                         *WIG3J1*FAKT1*2.D0*(2*L2+1)*TL(L2+1)
               END IF
               IF ( L.NE.0 .AND. L2.NE.0 ) TERM2(L+1) = TERM2(L+1)
     &              + BLL1(L+1,L1+1)
     &              *WIG3J2*FAKT2*2.D0*(TLP(L2+1)-TLM(L2+1))
            END DO
         END DO
         IF ( IREL.EQ.2 ) THEN
            TTLM(L+1) = (TERM1(L+1)+1.D0*(L+1)*TERM2(L+1))/(2*(2*L+1))
            TTLP(L+1) = (TERM1(L+1)-L*TERM2(L+1))/(2*(2*L+1))
            IF ( L.EQ.0 ) TTLP(1) = TERM1(1)/2.D0
C            TTL(L+1) = TERM1(L+1)/1.D0*(2*(2*L+1))
         END IF
      END DO
C
      IF ( III.EQ.1 ) THEN
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TSEL(IAP,IAP,IT) = TTLP(PTKE1(I))
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TSEL(IAP,IAP,IT) = TTLM(PTKO1(I))
         END DO
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TSOL(IAP,IAP,IT) = TTLP(PTKE1(I))
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TSOL(IAP,IAP,IT) = TTLM(PTKO1(I))
         END DO
      ELSE IF ( III.EQ.2 ) THEN
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TSEH(IAP,IAP,IT) = TTLP(PTKE1(I))
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TSEH(IAP,IAP,IT) = TTLM(PTKO1(I))
         END DO
         DO I = 1,MLQ1
            IAP = I + (IAN-1)*MLQ
            TSOH(IAP,IAP,IT) = TTLP(PTKE1(I))
         END DO
         DO I = MLQ1 + 1,MLQ
            IAP = I + (IAN-1)*MLQ
            TSOH(IAP,IAP,IT) = TTLM(PTKO1(I))
         END DO
      END IF
C
      END
