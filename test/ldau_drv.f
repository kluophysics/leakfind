C*==ldau_drv.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LDAU_DRV(GFMAT,DMFTSIGMA,WETAB,EREFLDAU,ELDAU,OBS_LT,
     &                    SCLNOS,WRSIG)
C
C   ********************************************************************
C   *                                                                  *
C   *  Driver for LDAU                                                 *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:IDOS,NKMMAX,NKM,NLMAX,NOBSMAX
      USE MOD_ENERGY,ONLY:NEMAX,NETAB
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP,LOPT,TXT_T,LTXT_T
      USE MOD_CALCMODE,ONLY:IREL,LLOYD
      USE MOD_FILES,ONLY:IFILGFWF,IPRINT,IOTMP
      USE MOD_SCF,ONLY:ITRSCF
      USE MOD_DMFT_LDAU,ONLY:JEFF,UEFF,KSELF,DMFTDBLC,UMODE
      USE MOD_CONSTANTS,ONLY:PI,C0,CI,RY_EV
      IMPLICIT NONE
C*--LDAU_DRV20
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='LDAU_DRV')
      REAL*8 TOL
      PARAMETER (TOL=1.0D-10)
C
C Dummy arguments
C
      REAL*8 SCLNOS
      LOGICAL WRSIG
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),EREFLDAU(NTMAX),
     &           GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX),WETAB(NEMAX)
      REAL*8 ELDAU(NTMAX),OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 COCC(:,:),COCCR(:,:),CTMP,CWORK(:,:),CWORK1(:,:),
     &           GZ(:,:,:),NS,NS1,NS2,RWORK(:,:),RWORK2(:,:),
     &           SIGSTAT(:,:),TMP1,TMP2,UR(:,:,:,:),URREP(:,:,:,:)
      REAL*8 EDC,EGM,UJJ,ULDAU(:,:,:,:,:),UUU
      CHARACTER*80 FILNAM,TXTT
      INTEGER I1,I2,I3,I4,IE,IL,IR,IR1,IR2,IRP,IS1,IS2,IT,LFILNAM,LM,
     &        LMP,LMPP,LMPPP,ML,NEGF,NLM,NREP
      LOGICAL RESCALE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ULDAU,GZ,COCC,CWORK,UR,URREP,SIGSTAT
      ALLOCATABLE RWORK,RWORK2,CWORK1,COCCR
C
      RESCALE = .FALSE.
C
      IF ( 1.EQ.2 ) WRITE (6,*) ITRSCF,LLOYD,SCLNOS
C
C=======================================================================
C                  calculate U-matrix
C=======================================================================
C
      ALLOCATE (ULDAU(2*NLMAX,2*NLMAX,2*NLMAX,2*NLMAX,NTMAX))
C
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).GT.0 ) THEN
            I1 = ITBOT
            I2 = ITTOP
            ITBOT = IT
            ITTOP = IT
C
            CALL DMFT_U_MATRIX(EREFLDAU,IPRINT,JEFF,IFILGFWF,UEFF,ULDAU)
C
            ITBOT = I1
            ITTOP = I2
C
         END IF
      END DO
C
      ELDAU(:) = 0.0D0
      NEGF = NETAB(1)
      LOOP_IT:DO IT = ITBOT,ITTOP
         IF ( LOPT(IT).GE.0 .AND. KSELF(IT).EQ.1 ) THEN
            NLM = LOPT(IT)*2 + 1
            NREP = 2*NLM
            I1 = LOPT(IT)**2 + 1
            I2 = I1 + NLM - 1
            I3 = NKM/2 + I1
            I4 = I3 + NLM - 1
            UUU = UEFF(IT)/RY_EV
            UJJ = JEFF(IT)/RY_EV
C
C
C==============================================================
C Transform U-matrix:
C -- change representation from U12,34 into U1,3,2,4
C -- block up  U-matrix to include spin dependence
C==============================================================
C
            ALLOCATE (UR(NLM,NLM,NLM,NLM))
            UR(1:NLM,1:NLM,1:NLM,1:NLM) = C0
C
            ALLOCATE (URREP(NREP,NREP,NREP,NREP))
            URREP(1:NREP,1:NREP,1:NREP,1:NREP) = C0
C
            DO LM = 1,NLM
               DO LMP = 1,NLM
                  DO LMPP = 1,NLM
                     DO LMPPP = 1,NLM
                        UR(LMPPP,LMPP,LMP,LM)
     &                     = DCMPLX(ULDAU(LMPPP,LMP,LMPP,LM,IT),0.0D0)
                     END DO
                  END DO
               END DO
            END DO
C
            DO IS1 = 1,2
               DO IS2 = 1,2
                  URREP((IS1-1)*NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM,
     &               (IS1-1)*NLM+1:IS1*NLM,(IS2-1)*NLM+1:IS2*NLM)
     &               = UR(1:NLM,1:NLM,1:NLM,1:NLM)
               END DO
            END DO
C==============================================================
            ALLOCATE (COCC(NREP,NREP))
            ALLOCATE (GZ(NREP,NREP,NEGF))
            COCC = C0
C
            GZ(1:NLM,1:NLM,1:NEGF) = GFMAT(I1:I2,I1:I2,1:NEGF,IT)
            GZ(NLM+1:NREP,NLM+1:NREP,1:NEGF)
     &         = GFMAT(I3:I4,I3:I4,1:NEGF,IT)
            GZ(1:NLM,NLM+1:NREP,1:NEGF) = GFMAT(I1:I2,I3:I4,1:NEGF,IT)
            GZ(NLM+1:NREP,1:NLM,1:NEGF) = GFMAT(I3:I4,I1:I2,1:NEGF,IT)
C
C
C==============================================================
C Calculate occupation numbers:
C==============================================================
C
            DO IR = 1,NREP
               DO IRP = 1,NREP
                  DO IE = 1,NEGF
                     COCC(IR,IRP) = COCC(IR,IRP)
     &                              + (GZ(IR,IRP,IE)*WETAB(IE)
     &                              -DCONJG(GZ(IRP,IR,IE)*WETAB(IE)))
                  END DO
               END DO
            END DO
            COCC(:,:) = -1.0D0/(2.0D0*PI*CI)*COCC(:,:)
            CTMP = C0
            DO IR = 1,NREP
               CTMP = CTMP + COCC(IR,IR)
            END DO
C
            IF ( IPRINT.GT.0 ) THEN
               CALL CMATSTRUCT('Occupation numbers: in CMPLX Y_LM',
     &                         COCC(1,1),NREP,NREP,0,0,0,TOL,6)
               WRITE (*,'(/,10X,A,2F20.10)') 'Trace:',CTMP
C  Print occupation numbers in Real Y_LM
               IF ( .NOT.ALLOCATED(CWORK1) )
     &              ALLOCATE (CWORK1(NKMMAX,NKMMAX))
               IF ( .NOT.ALLOCATED(CWORK) )
     &              ALLOCATE (CWORK(NKMMAX,NKMMAX))
C
               DO IE = 1,NEGF
                  CWORK(1:NKMMAX,1:NKMMAX)
     &               = GFMAT(1:NKMMAX,1:NKMMAX,IE,IT)
                  CALL CHANGEREP(NKM,NKMMAX,CWORK(1,1),'CLM>RLM',CWORK1)
                  GZ(1:NLM,1:NLM,IE) = CWORK1(I1:I2,I1:I2)
                  GZ(NLM+1:NREP,NLM+1:NREP,IE) = CWORK1(I3:I4,I3:I4)
                  GZ(1:NLM,NLM+1:NREP,IE) = CWORK1(I1:I2,I3:I4)
                  GZ(NLM+1:NREP,1:NLM,IE) = CWORK1(I3:I4,I1:I2)
               END DO
               ALLOCATE (COCCR(NREP,NREP))
               COCCR = C0
               DO IR = 1,NREP
                  DO IRP = 1,NREP
                     DO IE = 1,NEGF
                        COCCR(IR,IRP) = COCCR(IR,IRP)
     &                                  + (GZ(IR,IRP,IE)*WETAB(IE)
     &                                  -DCONJG(GZ(IRP,IR,IE)*WETAB(IE))
     &                                  )
                     END DO
                  END DO
               END DO
               COCCR(:,:) = -1.0D0/(2.0D0*PI*CI)*COCCR(:,:)
               CTMP = C0
               DO IR = 1,NREP
                  CTMP = CTMP + COCCR(IR,IR)
               END DO
               CALL CMATSTRUCT('Occupation numbers: in REAL  Y_LM',
     &                         COCCR(1,1),NREP,NREP,0,0,0,TOL,6)
               WRITE (*,'(/,10X,A,2F20.10)') 'Trace:',CTMP
               DEALLOCATE (COCCR,CWORK,CWORK1)
            END IF
C
            IF ( RESCALE ) THEN
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
               IL = LOPT(IT) + 1
               CTMP = OBS_LT(0,IDOS,IL,IT)/CTMP
               DO IR = 1,NREP
                  COCC(IR,IR) = COCC(IR,IR)*CTMP
               END DO
               CTMP = C0
               DO IR = 1,NREP
                  CTMP = CTMP + COCC(IR,IR)
               END DO
               IF ( IPRINT.GT.0 ) WRITE (*,'(/,10X,A,2F20.10)')
     &               'Trace: Rescaled',CTMP
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            END IF
C
C
C
C=======================================================================
C                  calculate spin and orbital moments
C                  calculate occupations for spin 1 and 2
C                  calculate total occupation
C=======================================================================
C
            TMP1 = C0
            TMP2 = C0
            ML = -LOPT(IT)
            DO IR = 1,NREP/2
               TMP1 = TMP1 + ML*COCC(IR,IR)
               TMP2 = TMP2 + ML*COCC(IR+NREP/2,IR+NREP/2)
               ML = ML + 1
            END DO
C
            IF ( IPRINT.GT.0 ) WRITE (*,'(/,10X,A,2F20.10)')
     &                                 'Orb. moment',TMP1 - TMP2
C
            TMP1 = C0
            TMP2 = C0
            DO IR = 1,NREP/2
               TMP1 = TMP1 + COCC(IR,IR)
               TMP2 = TMP2 + COCC(IR+NREP/2,IR+NREP/2)
            END DO
C
            IF ( IPRINT.GT.0 ) WRITE (*,'(/,10X,A,2F20.10)')
     &                                 'Spin moment',TMP1 - TMP2
C
            NS1 = TMP1
            NS2 = TMP2
            NS = TMP1 + TMP2
            TMP1 = 2D0*TMP1/NREP
            TMP2 = 2D0*TMP2/NREP
C=============================================================
C
C==============================================================
C Calculate LDA+U static self energy:
C==============================================================
C
            ALLOCATE (SIGSTAT(NREP,NREP))
            ALLOCATE (RWORK(NREP,NREP))
            ALLOCATE (RWORK2(NREP,NREP))
            SIGSTAT(1:NREP,1:NREP) = C0
C
            RWORK = C0
            RWORK2 = C0
C
C
            IF ( UMODE.EQ.'DUDA' ) THEN
C=============================================================
C     Spherically avaraged form of
C     fully rotational invariant LDA+U
C     (See Petukhov et al. PRB 67, 153106 (2003) and
C      Dudarev et al. PRB 57, 1505 (1998) )
C=============================================================
               IF ( DMFTDBLC(1:3).EQ.'AAL' ) THEN
                  DO IR = 1,NREP/2
                     DO IRP = 1,NREP/2
                        IF ( IR.EQ.IRP ) THEN
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(0.5D0-COCC(IR,IRP))
                        ELSE
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(-COCC(IR,IRP))
                        END IF
                     END DO
                  END DO
C
                  DO IR = NREP/2 + 1,NREP
                     DO IRP = NREP/2 + 1,NREP
                        IF ( IR.EQ.IRP ) THEN
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(0.5D0-COCC(IR,IRP))
                        ELSE
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(-COCC(IR,IRP))
                        END IF
                     END DO
                  END DO
C
C
               ELSE IF ( DMFTDBLC(1:3).EQ.'AMF' ) THEN
C
                  DO IR = 1,NREP/2
                     DO IRP = 1,NREP/2
                        IF ( IR.EQ.IRP ) THEN
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(TMP1-COCC(IR,IRP))
                        ELSE
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(-COCC(IR,IRP))
                        END IF
                     END DO
                  END DO
C
                  DO IR = NREP/2 + 1,NREP
                     DO IRP = NREP/2 + 1,NREP
                        IF ( IR.EQ.IRP ) THEN
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(TMP2-COCC(IR,IRP))
                        ELSE
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP) + (UUU-UJJ)
     &                        *(-COCC(IR,IRP))
                        END IF
                     END DO
                  END DO
C
               ELSE
                  CALL STOP_MESSAGE(ROUTINE,
     &       'This doublecounting is not yet implemented for UMODE DUDA'
     &       )
               END IF
C
C
            ELSE IF ( UMODE.EQ.'ROTI' ) THEN
C=============================================================
C
C     full static electron interaction:
C     SIGMA_{12}=\sum_{34}(V_{1324}-V_{1342})*n_{34}
C subtract double-counting /Tr<SIGMA>/:
C        AMF case:
C        SIGMA_{11}=SIGMA_{11}-\sum_{2}(V_{1212}-V_{1221})*<n_{2}>
C        AAL case:
C        SIGMA_{11}=SIGMA_{11}-U(N-1/2)+J(N_{1}-1/2)
C
               DO IR = 1,NREP
                  DO IRP = 1,NREP
                     DO IR1 = 1,NREP
                        DO IR2 = 1,NREP
                           SIGSTAT(IR,IRP) = SIGSTAT(IR,IRP)
     &                        + (URREP(IR,IR1,IRP,IR2)
     &                        -URREP(IR,IR1,IR2,IRP))*COCC(IR1,IR2)
                        END DO
                     END DO
                     RWORK2(IR,IRP) = SIGSTAT(IR,IRP)
C
                     IF ( IR.EQ.IRP ) THEN
                        IF ( DMFTDBLC(1:3).EQ.'AMF' ) THEN
                           DO IR1 = 1,NREP/2
                              SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                           - (URREP(IR,IR1,IR,IR1)
     &                           -URREP(IR,IR1,IR1,IR))*TMP1
                           END DO
                           DO IR1 = NREP/2 + 1,NREP
                              SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                           - (URREP(IR,IR1,IR,IR1)
     &                           -URREP(IR,IR1,IR1,IR))*TMP2
                           END DO
                        ELSE IF ( DMFTDBLC(1:3).EQ.'AAL' ) THEN
                           SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                        - UUU*((TMP1+TMP2)*NREP/2-0.5D0)
                           RWORK(IR,IR) = RWORK(IR,IR)
     &                        - UUU*((TMP1+TMP2)*NREP/2-0.5D0)
                           IF ( IR.LE.NREP/2 ) THEN
                              SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                           + UJJ*(TMP1*NREP/2-0.5D0)
                              RWORK(IR,IR) = RWORK(IR,IR)
     &                           + UJJ*(TMP1*NREP/2-0.5D0)
                           ELSE
                              SIGSTAT(IR,IR) = SIGSTAT(IR,IR)
     &                           + UJJ*(TMP2*NREP/2-0.5D0)
                              RWORK(IR,IR) = RWORK(IR,IR)
     &                           + UJJ*(TMP2*NREP/2-0.5D0)
                           END IF
                        ELSE
                           STOP 'LDA_DRV: CHECK DBLC !'
                        END IF
                     END IF
                  END DO
               END DO
               IF ( IPRINT.GT.0 .AND. DMFTDBLC(1:3).EQ.'AAL' )
     &              CALL CMATSTRUCT('SIGMA DC PART',RWORK,NREP,NREP,0,0,
     &                              0,TOL,6)
C
C
C=====================================================================
C           TOTAL ENERY TERM
C=====================================================================
C
               DO IR = 1,NREP
                  DO IRP = 1,NREP
                     DO IR1 = 1,NREP
                        DO IR2 = 1,NREP
                           ELDAU(IT) = ELDAU(IT)
     &                                 + DREAL((URREP(IR,IR1,IRP,IR2)
     &                                 -URREP(IR,IR1,IR2,IRP))
     &                                 *COCC(IR1,IR2)*COCC(IR,IRP))
                        END DO
                     END DO
                  END DO
               END DO
C
               ELDAU(IT) = (1.0D0/2.0D0)*ELDAU(IT)
               EGM = 0.0D0
               DO IR1 = 1,NREP
                  DO IR2 = 1,NREP
                     EGM = EGM + 0.5D0*DREAL(RWORK2(IR1,IR2))
     &                     *DREAL(COCC(IR2,IR1))
                  END DO
               END DO
               IF ( IPRINT.GT.0 ) PRINT *,'ENERGY EGM WITHOUT DC ',IT,
     &                                  EGM
C
C
               IF ( IPRINT.GT.0 ) PRINT *,
     &                                'ENERGY LDAU WITHOUT DC FOR TYPE '
     &                                ,IT,ELDAU(IT)
               IF ( DMFTDBLC(1:3).EQ.'AAL' ) THEN
C
                  EDC = (1.0D0/2.0D0)*UUU*DREAL(NS)*(DREAL(NS)-1.0D0)
     &                  - (1.0D0/2.0D0)*UJJ*DREAL(NS1)
     &                  *(DREAL(NS1)-1.0D0) - (1.0D0/2.0D0)
     &                  *UJJ*DREAL(NS2)*(DREAL(NS2)-1.0D0)
                  ELDAU(IT) = ELDAU(IT) - EDC
                  IF ( IPRINT.GT.0 ) THEN
                     PRINT *,'ENERGY DC FOR TYPE ',IT,EDC
                     PRINT *,'ENERGY LDAU-DC TYPE ',IT,ELDAU(IT)
                  END IF
C     In KKR eldau=-1/2Tr(Sigma.G) (See PRB 79, 115111 (2009)))
C
                  ELDAU(IT) = -ELDAU(IT)
               ELSE
                  ELDAU(IT) = 0.0D0
                  IF ( IPRINT.GT.0 ) WRITE (6,*) 
     &                 'WARNING: not yet LDA+U total energy for this DC'
               END IF
C
C
            END IF
C
            IF ( IPRINT.GT.0 ) CALL CMATSTRUCT('Static self-energy:',
     &           SIGSTAT,NREP,NREP,0,0,0,TOL,6)
C
C
C============================================================
C     print-out of the self energy
C============================================================
            IF ( WRSIG ) THEN
               FILNAM = TXT_T(IT)(1:LTXT_T(IT))//'_SIGUNFORM'//'.sig'
               LFILNAM = LTXT_T(IT) + 10 + 4
C
               OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),STATUS='UNKNOWN',
     &               FORM='unformatted')
C
               REWIND IOTMP
               WRITE (IOTMP) NREP
               WRITE (IOTMP) (((SIGSTAT(LM,LMP)),LM=1,NREP),LMP=1,NREP)
               WRITE (IOTMP) EREFLDAU(IT)
               CLOSE (IOTMP)
            END IF
C
C=======================================================
C     Change the representation from lms to kappa-mue
C=======================================================
C
            IF ( .NOT.ALLOCATED(CWORK) ) ALLOCATE (CWORK(NKMMAX,NKMMAX))
            CWORK(1:NKMMAX,1:NKMMAX) = C0
C
            CWORK(I1:I2,I1:I2) = SIGSTAT(1:NLM,1:NLM)
            CWORK(I3:I4,I3:I4) = SIGSTAT(NLM+1:NREP,NLM+1:NREP)
            CWORK(I1:I2,I3:I4) = SIGSTAT(1:NLM,NLM+1:NREP)
            CWORK(I3:I4,I1:I2) = SIGSTAT(NLM+1:NREP,1:NLM)
C
            IF ( IREL.EQ.3 ) THEN
               TXTT(1:7) = 'CLM>REL'
            ELSE
               TXTT(1:7) = 'CLM>RLM'
            END IF
            IF ( .NOT.ALLOCATED(CWORK1) )
     &           ALLOCATE (CWORK1(NKMMAX,NKMMAX))
C
            CWORK1 = C0
            DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,1:NEGF) = C0
C
            CALL CHANGEREP(NKM,NKMMAX,CWORK(1,1),TXTT(1:7),CWORK1)
C
            DO IE = 1,NEGF
               DMFTSIGMA(1:NKMMAX,1:NKMMAX,IT,IE)
     &            = CWORK1(1:NKMMAX,1:NKMMAX)
            END DO
C
            IF ( IPRINT.GT.0 ) THEN
               WRITE (6,*) 'IT=',IT,'IE=',1
               CALL CMATSTRUCT('DMFTSIGMA_REL in repres of KKR',
     &                         DMFTSIGMA(1,1,IT,1),NKMMAX,NKMMAX,IREL,
     &                         IREL,0,1D-8,6)
            END IF
            DEALLOCATE (GZ,COCC,SIGSTAT,UR,URREP,RWORK,RWORK2)
         END IF
      END DO LOOP_IT
C
C=======================================================
C Calculate new reference energy
C=======================================================
C
      CALL DMFT_LDAU_EREF(OBS_LT,EREFLDAU)
C
      IF ( WRSIG .OR. IPRINT.GT.0 ) THEN
         DO IT = ITBOT,ITTOP
            IF ( KSELF(IT).EQ.1 ) WRITE (6,99001) IT,
     &                            DREAL(EREFLDAU(IT))
         END DO
         WRITE (6,'(/)')
      END IF
C
99001 FORMAT (10X,'LDAU reference energy for IT',I3,' = ',F20.8)
      END
C
C
