C*==mecoul.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE MECOUL(ICST,GCOR,FCOR,MM05COR,NKPCOR,IKMCOR,LCXRAY,
     &                  IKMB1,IKMB2,IFILB,IKMC1,IKMC2,IFILC,IKMD1,IKMD2,
     &                  IFILD,IT,ME,IEB,IEC,AMECIG,AMECIF,NEMEMAX,
     &                  NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   * read wave function and calculate matrix elements for             *
C   *                                                                  *
C   *                    COULOMB INTERACTION                           *
C   *                                                                  *
C   *  22/01/98  HE                                                    *
C   ********************************************************************
      USE MOD_TYPES,ONLY:NCPLWFMAX,IMT,NTMAX,NKM_T
      USE MOD_ANGMOM,ONLY:NKM,L_IKM,MUEM05_IKM,NKMMAX,NLMAX
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI
      USE MOD_CONSTANTS,ONLY:PI
      IMPLICIT NONE
C*--MECOUL19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='MECOUL')
      REAL*8 THIRD
      PARAMETER (THIRD=1D0/3D0)
C
C Dummy arguments
C
      INTEGER ICST,IEB,IEC,IFILB,IFILC,IFILD,IKMB1,IKMB2,IKMC1,IKMC2,
     &        IKMD1,IKMD2,IT,NCSTMAX,NEMEMAX
      REAL*8 AMECIF(NKMMAX,NKMMAX,2*NLMAX),AMECIG(NKMMAX,NKMMAX,2*NLMAX)
     &       ,FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),LCXRAY(NTMAX),MM05COR(NCSTMAX),
     &        NKPCOR(NCSTMAX)
      COMPLEX*16 ME(NKMMAX,NKMMAX,NKMMAX,NEMEMAX,NEMEMAX)
C
C Local variables
C
      REAL*8 AME(2),PRE,RPWL(:,:)
      COMPLEX*16 CINTF(:),CINTG(:),FDFC,GDGC,RAMEDC(:,:),RMEF,RMEG,
     &           SG(:,:),SL(:,:),TG(:,:),TL(:,:),WR,ZFB(:,:,:),
     &           ZFC(:,:,:),ZFD(:,:,:),ZGB(:,:,:),ZGC(:,:,:),ZGD(:,:,:)
      INTEGER IA_ERR,IKMB(:,:),IKMC(:,:),IKMD(:,:),IKMTB,IKMTC,IKMTD,
     &        IKMWA,IKMWB,IKMWC,IKMWD,ILR,IM,IR,IRTOP,KA,KB,KC,KD,L,L1,
     &        L2,L3,LA,LC,LD,LR,MR,N,NGF,NKM_T_SAV,NSOLB(:),NSOLC(:),
     &        NSOLD(:)
      LOGICAL TRIANGLE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE IKMB,IKMC,IKMD,RPWL,CINTF,CINTG
      ALLOCATABLE NSOLB,NSOLC,NSOLD,SG,RAMEDC,TG,SL,TL,ZFB
      ALLOCATABLE ZFC,ZFD,ZGB,ZGC,ZGD
C
      TRIANGLE(L1,L2,L3) = (L1.GE.ABS(L3-L2)) .AND. (L1.LE.(L3+L2))
     &                     .AND. (MOD((L1+L2+L3),2).EQ.0)
C
      ALLOCATE (IKMB(NCPLWFMAX,NKMMAX),NSOLB(NKMMAX))
      ALLOCATE (IKMC(NCPLWFMAX,NKMMAX),NSOLC(NKMMAX))
      ALLOCATE (IKMD(NCPLWFMAX,NKMMAX),NSOLD(NKMMAX))
      ALLOCATE (CINTF(NRMAX),CINTG(NRMAX),RPWL(NRMAX,0:2*NLMAX))
      ALLOCATE (SL(NRMAX,2),TL(NRMAX,2),TG(NRMAX,2))
      ALLOCATE (SG(NRMAX,2),RAMEDC(NRMAX,2*NLMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: MJC')
C
      ALLOCATE (ZFB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFC(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZFD(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZFC')
C
      ALLOCATE (ZGB(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGC(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZGD(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZGD')
C
      NGF = 2
C
      CALL CINIT(NRMAX*2*NLMAX,RAMEDC)
C
      IM = IMT(IT)
      IRTOP = JRWS(IM)
      DO IR = 1,NRMAX
         RPWL(IR,0) = 1.0D0
         RPWL(IR,1) = R(IR,IM)
         DO L = 2,2*NLMAX
            RPWL(IR,L) = RPWL(IR,L-1)*R(IR,IM)
         END DO
      END DO
C
C --------------------------- read in REGULAR wavefunctions for energy C
C
      CALL WAVFUN_READ_REL(IFILC,IT,0,ZGC,ZFC,ZGC,ZFC,IRTOP,NSOLC,IKMC)
C
C --------------------------- read in REGULAR wavefunctions for energy B
C
      CALL WAVFUN_READ_REL(IFILB,IT,0,ZGB,ZFB,ZGB,ZFB,IRTOP,NSOLB,IKMB)
C
C --------------------------- read in REGULAR wavefunctions for energy D
C
      NKM_T_SAV = NKM_T(IT)
      NKM_T(IT) = NKM
C
      CALL WAVFUN_READ_REL(IFILD,IT,0,ZGD,ZFD,ZGD,ZFD,IRTOP,NSOLD,IKMD)
C
      NKM_T(IT) = NKM_T_SAV
C
C DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      DO IKMTD = IKMD1,IKMD2
         DO KD = 1,NSOLD(IKMTD)
            IKMWD = IKMD(KD,IKMTD)
            LD = L_IKM(IKMTD)
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DO IKMTC = IKMC1,IKMC2
               DO KC = 1,NSOLC(IKMTC)
                  IKMWC = IKMC(KC,IKMTC)
C LRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLR
                  ILR = 0
                  CALL CINIT(NRMAX*2*NLMAX,RAMEDC)
                  LC = L_IKM(IKMTC)
                  DO LR = ABS(LD-LC),(LD+LC),2
                     ILR = LR + 1
                     AME(1) = AMECIG(IKMWD,IKMWC,ILR)
                     AME(2) = AMECIF(IKMWD,IKMWC,ILR)
C
C ---------------------------------------------------- set up integrands
                     DO IR = 1,IRTOP
                        GDGC = ZGD(IR,KD,IKMTD)*ZGC(IR,KC,IKMTC)
                        FDFC = ZFD(IR,KD,IKMTD)*ZFC(IR,KC,IKMTC)
                        WR = R2DRDI(IR,IM)*RPWL(IR,LR)
                        TL(IR,1) = GDGC*WR
                        TL(IR,2) = FDFC*WR
                        WR = R2DRDI(IR,IM)/RPWL(IR,LR+1)
                        TG(IR,1) = GDGC*WR
                        TG(IR,2) = FDFC*WR
                     END DO
C -------------------------------------------- evaluate radial integrals
                     DO N = 1,NGF
                        SL(1,N) = 0.0D0
                        SG(1,N) = 0.0D0
                        DO IR = 3,IRTOP,2
                           SL(IR,N) = SL(IR-2,N)
     &                                + THIRD*(TL(IR-2,N)+4.0D0*TL(IR-1,
     &                                N)+TL(IR,N))
                           SL(IR-1,N) = SL(IR,N)
     &                                  - 0.5D0*(TL(IR-1,N)+TL(IR,N))
C
                           SG(IR,N) = SG(IR-2,N)
     &                                + THIRD*(TG(IR-2,N)+4.0D0*TG(IR-1,
     &                                N)+TG(IR,N))
                           SG(IR-1,N) = SG(IR,N)
     &                                  - 0.5D0*(TG(IR-1,N)+TG(IR,N))
                        END DO
                     END DO
C ------------------------------------ add prefactor 1/r**(L+1) and r**L
                     DO N = 1,NGF
                        DO IR = 1,IRTOP
                           SL(IR,N) = SL(IR,N)/RPWL(IR,LR+1)
                           SG(IR,N) = (SG(IRTOP,N)-SG(IR,N))*RPWL(IR,LR)
                           RAMEDC(IR,ILR) = RAMEDC(IR,ILR)
     &                        + (SL(IR,N)+SG(IR,N))*AME(N)
                        END DO
                     END DO
C ----------------------------------------------------------------------
                  END DO
C LRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLR
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                  DO IKMTB = IKMB1,IKMB2
                     DO KB = 1,NSOLB(IKMTB)
                        IKMWB = IKMB(KB,IKMTB)
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
                        LA = LCXRAY(IT)
C
                        DO KA = 1,NKPCOR(ICST)
                           IKMWA = IKMCOR(ICST,KA)
C
C MRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMR
C
                           MR = MUEM05_IKM(IKMTD) - MUEM05_IKM(IKMTC)
C
                           IF ( -MR.EQ.(MM05COR(ICST)-MUEM05_IKM(IKMTB))
     &                          ) THEN
C
C LRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLR
                              ILR = 0
                              DO LR = ABS(LC-LD),(LC+LD),2
                                 ILR = LR + 1
                                 IF ( TRIANGLE(LR,LA,L_IKM(IKMTB)) )
     &                                THEN
C
                                    PRE = (4*PI/DBLE(2*LR+1))
                                    DO IR = 1,IRTOP
                                       CINTG(IR) = GCOR(IR,KA,ICST)
     &                                    *ZGB(IR,KB,IKMTB)
     &                                    *RAMEDC(IR,ILR)*R2DRDI(IR,IM)
                                       CINTF(IR) = FCOR(IR,KA,ICST)
     &                                    *ZFB(IR,KB,IKMTB)
     &                                    *RAMEDC(IR,ILR)*R2DRDI(IR,IM)
                                    END DO
C
                                    CALL CRADINT(IM,CINTG,RMEG)
                                    CALL CRADINT(IM,CINTF,RMEF)
C
                                    ME(IKMTB,IKMTC,IKMTD,IEB,IEC)
     &                                 = ME(IKMTB,IKMTC,IKMTD,IEB,IEC)
     &                                 + PRE*
     &                                 (RMEG*AMECIG(IKMWB,IKMWA,ILR)
     &                                 +RMEF*AMECIF(IKMWB,IKMWA,ILR))
                                 END IF
                              END DO
C LRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLRLR
C
                           END IF
C MRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMRMR
C
                        END DO
C     END DO
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
                     END DO
                  END DO
C BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               END DO
            END DO
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         END DO
      END DO
C DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
      DEALLOCATE (IKMB,IKMC,IKMD,RPWL,CINTF,CINTG)
      DEALLOCATE (NSOLB,NSOLC,NSOLD,SG,RAMEDC,TG,SL,TL,ZFB)
      DEALLOCATE (ZFC,ZFD,ZGB,ZGC,ZGD)
C
      END
