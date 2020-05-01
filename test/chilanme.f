C*==chilanme.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHILANME(NKMWF,IFILWF,IT,AMEA1,AMEA2,AMEB1,AMEB2,
     &                    AMEB1C,AMEB2C,MEA,MEB)
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements                 *
C   *                                                                  *
C   *               < LAM_f |     ->nabla * ->A_lam | LAM_i >          *
C   *                                                                  *
C   *               < LAM_f | y * ->nabla * ->A_lam | LAM_i >          *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_TYPES,ONLY:NT,IMT,NCPLWFMAX
      USE MOD_ANGMOM,ONLY:CGC,NKMMAX,NKM
      USE MOD_CONSTANTS,ONLY:CI,PI
      IMPLICIT NONE
C*--CHILANME19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='CHILANME')
      REAL*8 NON0_TOL
      PARAMETER (NON0_TOL=1D-8)
C
C Dummy arguments
C
      INTEGER IFILWF,IT,NKMWF
      REAL*8 AMEA1(NKMMAX,NKMMAX,3),AMEA2(NKMMAX,NKMMAX,3),
     &       AMEB1(NKMMAX,NKMMAX,3,3),AMEB1C(NKMMAX,NKMMAX,3,3),
     &       AMEB2(NKMMAX,NKMMAX,3,3),AMEB2C(NKMMAX,NKMMAX,3,3)
      COMPLEX*16 MEA(NKMMAX,NKMMAX,3,NT),MEB(NKMMAX,NKMMAX,3,3,NT)
C
C Local variables
C
      REAL*8 AF1,AF2,AF3,AF4,MJF,MJI,RARG
      COMPLEX*16 CCDN,CCUP,CINT1(:),CINT2(:),CINT3(:),CINT4(:),CINT5(:),
     &           CIW,MTMP(3,3),R3DRDI,RMEF1(2,2),RMEF2(2,2),RMEF3(2,2),
     &           RMEF4(2,2),RMEG1(2,2),RMEG2(2,2),RMEG3(2,2),RMEG4(2,2),
     &           RMEG5(2,2),WZ2,ZFF(:,:,:),ZFFP(:,:,:),ZFI(:,:,:),
     &           ZFIP(:,:,:),ZGF(:,:,:),ZGFP(:,:,:),ZGI(:,:,:),
     &           ZGIP(:,:,:)
      INTEGER I,IA_ERR,IKMF(2),IKMI(2),IKMTF,IKMTI,IM,IMKMF,IMKMI,IPOL,
     &        IPOL1,IPOL2,IPOLC,IRTOP,ITF,ITI,J,JF,JI,K,KAPF(2),KAPI(2),
     &        KF,KI,LBF,LBI,LF,LI,NPOL,NSOLF,NSOLI
      LOGICAL RNON0
      CHARACTER*3 STR3
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZFFP,ZGFP,ZGIP,ZFIP,CINT1,CINT2,CINT3,CINT4,CINT5
      ALLOCATABLE ZGF,ZFF,ZFI,ZGI
C
      RNON0(RARG) = ABS(RARG).GT.NON0_TOL
C
      ALLOCATE (CINT1(NRMAX),CINT2(NRMAX),CINT3(NRMAX))
      ALLOCATE (CINT4(NRMAX),CINT5(NRMAX))
      ALLOCATE (ZFFP(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGFP(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGIP(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFIP(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFI(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGI(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGI')
C
      NPOL = 3
      IM = IMT(IT)
      IRTOP = JRWS(IM)
C
      CALL CINIT(NKMMAX*NKMMAX*NPOL,MEA(1,1,1,IT))
      CALL CINIT(NKMMAX*NKMMAX*NPOL*NPOL,MEB(1,1,1,1,IT))
C
C=======================================================================
C
      DO IKMTF = 1,NKMWF
         READ (IFILWF,REC=IKMTF+(IT-1)*NKM) ITF,LF,MJF,NSOLF,STR3,
     &         (KAPF(K),IKMF(K),
     &         (ZGF(I,K,IKMTF),ZFF(I,K,IKMTF),I=1,IRTOP),K=1,NSOLF)
         IF ( IT.NE.ITF .OR. STR3.NE.'REG' ) THEN
            WRITE (6,*) 'it',IT,STR3
            WRITE (6,*) 'F: it l kap mj ',ITF,LF,MJF,STR3,KAPF
            CALL STOP_MESSAGE(ROUTINE,'WF(FIN) inconsistent')
         END IF
C
C-----------------------------------------------------------------------
C
         DO KF = 1,NSOLF
            CALL CDIFFER(IM,ZGF(1,KF,IKMTF),ZGFP(1,KF,IKMTF))
            CALL CDIFFER(IM,ZFF(1,KF,IKMTF),ZFFP(1,KF,IKMTF))
         END DO
C
C-----------------------------------------------------------------------
C
         DO IKMTI = 1,NKMWF
            READ (IFILWF,REC=IKMTI+(IT-1)*NKM) ITI,LI,MJI,NSOLI,STR3,
     &            (KAPI(K),IKMI(K),
     &            (ZGI(I,K,IKMTI),ZFI(I,K,IKMTI),I=1,IRTOP),K=1,NSOLI)
            IF ( IT.NE.ITI .OR. STR3.NE.'REG' ) THEN
               WRITE (6,*) 'it',IT,STR3
               WRITE (6,*) 'I: it l kap mj ',ITI,LI,MJI,STR3,KAPI
               CALL STOP_MESSAGE(ROUTINE,'WF(INI) inconsistent')
            END IF
C
            DO KI = 1,NSOLI
               CALL CDIFFER(IM,ZGI(1,KI,IKMTI),ZGIP(1,KI,IKMTI))
               CALL CDIFFER(IM,ZFI(1,KI,IKMTI),ZFIP(1,KI,IKMTI))
            END DO
C
C=======================================================================
C                     <LAM_f | ->nabla * ->A_lam | LAM_i>
C=======================================================================
C
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  DO IPOL = 1,NPOL
                     IF ( RNON0(AMEA1(JF,JI,IPOL)) ) GOTO 20
                     IF ( RNON0(AMEA2(JF,JI,IPOL)) ) GOTO 20
                  END DO
               END DO
            END DO
C -------------------------------------- all angular matrix elements = 0
            GOTO 40
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
 20         CONTINUE
            DO KI = 1,NSOLI
               LBI = LI - SIGN(1,KAPI(KI))
               DO KF = 1,NSOLF
                  LBF = LF - SIGN(1,KAPF(KF))
                  DO I = 1,IRTOP
                     CINT1(I) = DCONJG(ZGF(I,KF,IKMTF))*R2DRDI(I,IM)
     &                          *(ZGIP(I,KI,IKMTI)-ZGI(I,KI,IKMTI)
     &                          *LI/R(I,IM))
                     CINT2(I) = DCONJG(ZGF(I,KF,IKMTF))*R2DRDI(I,IM)
     &                          *(ZGIP(I,KI,IKMTI)+ZGI(I,KI,IKMTI)
     &                          *(LI+1)/R(I,IM))
                     CINT3(I) = (ZGI(I,KI,IKMTI)*R2DRDI(I,IM)*DCONJG(
     &                          ZGFP(I,KF,IKMTF)-ZGF(I,KF,IKMTF)
     &                          *LF/R(I,IM)))
                     CINT4(I) = (ZGI(I,KI,IKMTI)*R2DRDI(I,IM)*DCONJG(
     &                          ZGFP(I,KF,IKMTF)+ZGF(I,KF,IKMTF)*(LF+1)
     &                          /R(I,IM)))
                  END DO
                  CALL CRADINT(IM,CINT1,RMEG1(KF,KI))
                  CALL CRADINT(IM,CINT2,RMEG2(KF,KI))
                  CALL CRADINT(IM,CINT3,RMEG3(KF,KI))
                  CALL CRADINT(IM,CINT4,RMEG4(KF,KI))
C
                  DO I = 1,IRTOP
                     CINT1(I) = DCONJG(ZFF(I,KF,IKMTF))*R2DRDI(I,IM)
     &                          *(ZFIP(I,KI,IKMTI)-ZFI(I,KI,IKMTI)
     &                          *LBI/R(I,IM))
                     CINT2(I) = DCONJG(ZFF(I,KF,IKMTF))*R2DRDI(I,IM)
     &                          *(ZFIP(I,KI,IKMTI)+ZFI(I,KI,IKMTI)
     &                          *(LBI+1)/R(I,IM))
                     CINT3(I) = ZFI(I,KI,IKMTI)*R2DRDI(I,IM)
     &                          *DCONJG(ZFFP(I,KF,IKMTF)-ZFF(I,KF,IKMTF)
     &                          *LBF/R(I,IM))
                     CINT4(I) = ZFI(I,KI,IKMTI)*R2DRDI(I,IM)
     &                          *DCONJG(ZFFP(I,KF,IKMTF)+ZFF(I,KF,IKMTF)
     &                          *(LBF+1)/R(I,IM))
                  END DO
                  CALL CRADINT(IM,CINT1,RMEF1(KF,KI))
                  CALL CRADINT(IM,CINT2,RMEF2(KF,KI))
                  CALL CRADINT(IM,CINT3,RMEF3(KF,KI))
                  CALL CRADINT(IM,CINT4,RMEF4(KF,KI))
C
               END DO
            END DO
C
C -------------------------------------- calculate total matrix elements
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               LBI = LI - SIGN(1,KAPI(KI))
               IMKMI = LBI*2*ABS(KAPI(KI)) + ABS(KAPI(KI))
     &                 + NINT(-0.5D0+MJI) + 1
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  LBF = LF - SIGN(1,KAPF(KF))
                  IMKMF = LBF*2*ABS(KAPF(KF)) + ABS(KAPF(KF))
     &                    + NINT(-0.5D0+MJF) + 1
                  DO IPOL = 1,NPOL
                     IF ( IPOL.EQ.1 ) IPOLC = 2
                     IF ( IPOL.EQ.2 ) IPOLC = 1
                     IF ( (IMKMF.LE.NKMWF) .AND. (IMKMI.LE.NKMWF) ) THEN
                        AF1 = AMEA1(IMKMF,IMKMI,IPOL)
                        AF2 = AMEA2(IMKMF,IMKMI,IPOL)
                        AF3 = AMEA2(IMKMF,IMKMI,IPOL)
                        AF4 = AMEA1(IMKMF,IMKMI,IPOL)
                     ELSE
                        AF1 = 0.0D0
                        AF2 = 0.0D0
                        AF3 = 0.0D0
                        AF4 = 0.0D0
                     END IF
                     MEA(IKMTF,IKMTI,IPOL,IT) = MEA(IKMTF,IKMTI,IPOL,IT)
     &                  + (RMEG1(KF,KI)*AMEA1(JF,JI,IPOL)-RMEG2(KF,KI)
     &                  *AMEA2(JF,JI,IPOL)+RMEF1(KF,KI)*AF1-RMEF2(KF,KI)
     &                  *AF2+(RMEG3(KF,KI)*AMEA1(JI,JF,IPOLC))
     &                  -RMEG4(KF,KI)*AMEA2(JI,JF,IPOLC)+RMEF3(KF,KI)
     &                  *AF3-RMEF4(KF,KI)*AF4)*0.5D0
C
                  END DO
               END DO
            END DO
C
C=======================================================================
C                     <LAM_f | y ->nabla * ->A_lam | LAM_i>
C=======================================================================
C
 40         CONTINUE
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  DO IPOL1 = 1,NPOL
                     DO IPOL2 = 1,NPOL
                        IF ( RNON0(AMEB1(JF,JI,IPOL1,IPOL2)) ) GOTO 60
                        IF ( RNON0(AMEB2(JF,JI,IPOL1,IPOL2)) ) GOTO 60
                     END DO
                  END DO
               END DO
            END DO
C -------------------------------------- all angular matrix elements = 0
            CYCLE
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
 60         CONTINUE
            DO KI = 1,NSOLI
               LBI = LI - SIGN(1,KAPI(KI))
               DO KF = 1,NSOLF
                  LBF = LF - SIGN(1,KAPF(KF))
                  DO I = 1,IRTOP
                     R3DRDI = R(I,IM)*R2DRDI(I,IM)
                     CINT1(I) = DCONJG(ZGF(I,KF,IKMTF))
     &                          *R3DRDI*(ZGIP(I,KI,IKMTI)
     &                          -ZGI(I,KI,IKMTI)*LI/R(I,IM))
                     CINT2(I) = DCONJG(ZGF(I,KF,IKMTF))
     &                          *R3DRDI*(ZGIP(I,KI,IKMTI)
     &                          +ZGI(I,KI,IKMTI)*(LI+1)/R(I,IM))
                     CINT3(I) = (ZGI(I,KI,IKMTI)*R3DRDI*DCONJG(ZGFP(I,KF
     &                          ,IKMTF)-ZGF(I,KF,IKMTF)*LF/R(I,IM)))
                     CINT4(I) = (ZGI(I,KI,IKMTI)*R3DRDI*DCONJG(ZGFP(I,KF
     &                          ,IKMTF)+ZGF(I,KF,IKMTF)*(LF+1)/R(I,IM)))
                     CINT5(I) = ZGI(I,KI,IKMTI)*DCONJG(ZGF(I,KF,IKMTF))
     &                          *R2DRDI(I,IM)
                  END DO
                  CALL CRADINT(IM,CINT1,RMEG1(KF,KI))
                  CALL CRADINT(IM,CINT2,RMEG2(KF,KI))
                  CALL CRADINT(IM,CINT3,RMEG3(KF,KI))
                  CALL CRADINT(IM,CINT4,RMEG4(KF,KI))
                  CALL CRADINT(IM,CINT5,RMEG5(KF,KI))
C
                  DO I = 1,IRTOP
                     R3DRDI = R(I,IM)*R2DRDI(I,IM)
                     CINT1(I) = ZFF(I,KF,IKMTF)
     &                          *R3DRDI*(ZFIP(I,KI,IKMTI)-
     &                          ZFI(I,KI,IKMTI)*LBI/R(I,IM))
                     CINT2(I) = ZFF(I,KF,IKMTF)
     &                          *R3DRDI*(ZFIP(I,KI,IKMTI)+
     &                          ZFI(I,KI,IKMTI)*(LBI+1)/R(I,IM))
                     CINT3(I) = ZFI(I,KI,IKMTI)
     &                          *R3DRDI*(ZFFP(I,KF,IKMTF)-
     &                          ZFF(I,KF,IKMTF)*LBF/R(I,IM))
                     CINT4(I) = ZFI(I,KI,IKMTI)
     &                          *R3DRDI*(ZFFP(I,KF,IKMTF)+
     &                          ZFF(I,KF,IKMTF)*(LBF+1)/R(I,IM))
                  END DO
C
                  CALL CRADINT(IM,CINT1,RMEF1(KF,KI))
                  CALL CRADINT(IM,CINT2,RMEF2(KF,KI))
                  CALL CRADINT(IM,CINT3,RMEF3(KF,KI))
                  CALL CRADINT(IM,CINT4,RMEF4(KF,KI))
C
               END DO
            END DO
C
C -------------------------------------- calculate total matrix elements
            DO KI = 1,NSOLI
               JI = IKMI(KI)
               LBI = LI - SIGN(1,KAPI(KI))
               IMKMI = LBI*2*ABS(KAPI(KI)) + ABS(KAPI(KI))
     &                 + NINT(-0.5D0+MJI) + 1
               DO KF = 1,NSOLF
                  JF = IKMF(KF)
                  LBF = LF - SIGN(1,KAPF(KF))
                  IMKMF = LBF*2*ABS(KAPF(KF)) + ABS(KAPF(KF))
     &                    + NINT(-0.5D0+MJF) + 1
C ----------------------------------- Polarisation of radius-vector operator
                  DO IPOL1 = 1,NPOL
C ----------------------------------------- Polarisation of a field operator
                     DO IPOL2 = 1,NPOL
C
                        IF ( (IMKMF.LE.NKMWF) .AND. (IMKMI.LE.NKMWF) )
     &                       THEN
                           AF1 = AMEB1(IMKMF,IMKMI,IPOL1,IPOL2)
                           AF2 = AMEB2(IMKMF,IMKMI,IPOL1,IPOL2)
                           AF3 = AMEB1C(IMKMF,IMKMI,IPOL1,IPOL2)
                           AF4 = AMEB2C(IMKMF,IMKMI,IPOL1,IPOL2)
                        ELSE
                           AF1 = 0.0D0
                           AF2 = 0.0D0
                           AF3 = 0.0D0
                           AF4 = 0.0D0
                        END IF
                        MEB(IKMTF,IKMTI,IPOL1,IPOL2,IT)
     &                     = MEB(IKMTF,IKMTI,IPOL1,IPOL2,IT)
     &                     + (RMEG1(KF,KI)*AMEB1(JF,JI,IPOL1,IPOL2)
     &                     -RMEG2(KF,KI)*AMEB2(JF,JI,IPOL1,IPOL2)
     &                     +RMEF1(KF,KI)*AF1-RMEF2(KF,KI)
     &                     *AF2+RMEG3(KF,KI)*AMEB1C(JI,JF,IPOL1,IPOL2)
     &                     -RMEG4(KF,KI)*AMEB2C(JI,JF,IPOL1,IPOL2)
     &                     +RMEF3(KF,KI)*AF3-RMEF4(KF,KI)*AF4)*0.5D0
C
                        IF ( (JI.EQ.JF) .AND. (IPOL1+IPOL2.EQ.3) ) THEN
C
                           CCDN = CGC(IKMTF,1)*CGC(IKMTI,1)
                           CCUP = CGC(IKMTF,2)*CGC(IKMTI,2)
C
                           MEB(IKMTF,IKMTI,IPOL1,IPOL2,IT)
     &                        = MEB(IKMTF,IKMTI,IPOL1,IPOL2,IT)
     &                        + DSQRT(DBLE(1.)/DBLE(3.))*(-1.)
     &                        **IPOL2*DSQRT(DBLE(1.)/DBLE(4.)
     &                        /DBLE(3.1415))*(CCDN+CCUP)*RMEG5(KF,KI)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
C
C=======================================================================
C
         END DO
      END DO
C
C=======================================================================
C                   convert to cartesian coordinates
C=======================================================================
C
      WZ2 = DSQRT(2.0D0)
      CIW = CI*DSQRT(PI/3.0D0)
C
      DO J = 1,NKMMAX
         DO I = 1,NKMMAX
C
C                index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C                 M(X) =   [  M(+) + M(-) ] / SQRT(2)
C                 M(Y) = I*[ -M(+) + M(-) ] / SQRT(2)
C
            DO IPOL = 1,3
               MTMP(IPOL,1) = MEA(I,J,IPOL,IT)
            END DO
C
            MEA(I,J,1,IT) = (-MTMP(1,1)+MTMP(2,1))/WZ2
            MEA(I,J,2,IT) = CI*(MTMP(1,1)+MTMP(2,1))/WZ2
            MEA(I,J,3,IT) = MTMP(3,1)
C
            DO IPOL1 = 1,3
               DO IPOL2 = 1,3
                  MTMP(IPOL1,IPOL2) = MEB(I,J,IPOL1,IPOL2,IT)
               END DO
            END DO
C
C----------------------------------------- store only elements xy and yx
C
            DO IPOL1 = 1,3
               DO IPOL2 = 1,3
                  MEB(I,J,IPOL1,IPOL2,IT) = 0D0
               END DO
            END DO
C
            MEB(I,J,1,2,IT) = CIW*(MTMP(1,1)+MTMP(1,2)+MTMP(2,1)+MTMP(2,
     &                        2))
            MEB(I,J,2,1,IT) = CIW*(MTMP(1,1)-MTMP(1,2)-MTMP(2,1)+MTMP(2,
     &                        2))
C
         END DO
      END DO
C
C
      DEALLOCATE (ZFFP,ZGFP,ZGIP,ZFIP,CINT1,CINT2,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC: CINT2')
C
      DEALLOCATE (CINT3,CINT4,ZGF,ZFF,ZFI,ZGI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC: ZGI')
C
      END
