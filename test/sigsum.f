C*==sigsum.f    processed by SPAG 6.70Rc at 13:20 on 21 Apr 2017
      SUBROUTINE SIGSUM(SIG1Q_FSEA_NV,SIG1Q_FSEA_VC,SIG0Q_FSEA,
     &                  SIG1Q_SURF_NV,SIG1Q_SURF_VC,SIG0Q_SURF,
     &                  SIGMAAU_NV,SIGMAAU_VC,RHOAU_NV,RHOAU_VC,
     &                  SIGMA_NV,SIGMA_VC,RHO_NV,RHO_VC,NSPINPROJ,
     &                  SIGOFFQ,SIG_MODE,NEA,I_TEMP_LAT,CHANGE_FRAME,
     &                  MROT_FRAME)
C   ********************************************************************
C   *                                                                  *
C   *   sum up results                                                 *
C   *   calculate final resistivity and conductivity                   *
C   *   scalar  and  tensor  form                                      *
C   *                                                                  *
C   *   NOTE: JB gives HBAR etc. in SI units                           *
C   *                                                                  *
C   * 06/10/99  HE  based on JB's routine                              *
C   ********************************************************************
C   *                                                                  *
C   *  KMROT                                                           *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site IQ   *
C   *  2: global COMMON rotation of the magnetisation                  *
C   *  3: spin spiral    Theta =  90    not allowed                    *
C   *  4: spin spiral    Theta <> 90    not allowed                    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_RMESH,ONLY:VOL_AC
      USE MOD_LATTICE,ONLY:SYSTEM_TYPE
      USE MOD_THERMAL,ONLY:X_VFT
      USE MOD_FILES,ONLY:IPRINT,DATSET,LDATSET,IOTMP,IOTMP4,GIT_HASH,
     &    GIT_COMPILE_DATE,GIT_BRANCH,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CALCMODE,ONLY:KMROT,PUBLIC_VERSION
      USE MOD_SITES,ONLY:NQMAX,MROTQ,NQ,IQBOT_CHI,IQTOP_CHI,IMQ
      USE MOD_SIG,ONLY:STR_ISP_PROJ,LIST_ISPR,NSPR,CONSI,SIG_PREFAC_AU,
     &    SOTSI,SOT_PREFAC_AU,EESI,EE_PREFAC_AU,IRESPONSE_SOT,
     &    IRESPONSE_EDELSTEIN,IRESPONSE_ADA_ADA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGSUM')
C
C Dummy arguments
C
      LOGICAL CHANGE_FRAME
      INTEGER I_TEMP_LAT,NEA,NSPINPROJ
      REAL*8 RHOAU_NV,RHOAU_VC,SIGMAAU_NV,SIGMAAU_VC
      CHARACTER*10 SIG_MODE
      REAL*8 MROT_FRAME(3,3),RHO_NV(3,3),RHO_VC(3,3),SIGMA_NV(3,3),
     &       SIGMA_VC(3,3)
      COMPLEX*16 SIG0Q_FSEA(3,3,NSPINPROJ,NQMAX),
     &           SIG0Q_SURF(3,3,NSPINPROJ,NQMAX),
     &           SIG1Q_FSEA_NV(3,3,NSPINPROJ,NQMAX,NQMAX),
     &           SIG1Q_FSEA_VC(3,3,NSPINPROJ,NQMAX,NQMAX),
     &           SIG1Q_SURF_NV(3,3,NSPINPROJ,NQMAX,NQMAX),
     &           SIG1Q_SURF_VC(3,3,NSPINPROJ,NQMAX,NQMAX),
     &           SIGOFFQ(3,3,NSPINPROJ,NQMAX)
C
C Local variables
C
      COMPLEX*16 CSIGD(:,:),SIG0(3,3),SIG0EF(3,3),SIG1_NV(3,3),
     &           SIG1_NVEF(3,3),SIG1_VC(3,3),SIG1_VCEF(3,3),SIGOFF(3,3),
     &           SIGTOT_NV(3,3),SIGTOT_NVEF(3,3),SIGTOT_VC(3,3),
     &           SIGTOT_VCEF(3,3),SMA(2)
      REAL*8 DIFF,DUMMAT(3,3),MTMP(3,3),PREAU,PREAU_NQ,PRESI,PRESI_NQ,
     &       RESIST_NV,RESIST_VC,RHOL_NV(3,3),RHOL_VC(3,3),RHO_NVEF(3,3)
     &       ,RHO_NV_SI(3,3),RHO_VCEF(3,3),RHO_VC_SI(3,3),RSIGD(:,:),
     &       RSIGDL(:,:),SIG1AU_NV,SIG1AU_VC,SIGMAL_NV(3,3),
     &       SIGMAL_VC(3,3),SIGMA_NVEF(3,3),SIGMA_NVEFSAVE(3,3),
     &       SIGMA_NVSAVE(3,3),SIGMA_NVSAVE_SI(3,3),SIGMA_VCEF(3,3),
     &       SIGMA_VCEFSAVE(3,3),SIGMA_VCSAVE(3,3),SIGMA_VCSAVE_SI(3,3)
      CHARACTER*255 DS1,RESPSTR
      CHARACTER*80 FILNAM,LFILNAM
      INTEGER I,IM,IQ1,IQ2,IQB,IQT,ISPINPROJ,ISPR,IV,J,MUE,NQ_AVERAGE,
     &        NUE
      LOGICAL LPRSIG,L_SIG_SINGULAR
      CHARACTER*1 MSTR,NSTR
      REAL*8 RMAT3X3DET
      CHARACTER*2 STR_ISPR
      CHARACTER*5 STR_I_TEMP_LAT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE RSIGD,RSIGDL,CSIGD
C
      CALL TRACK_INFO(ROUTINE)
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      IF ( PUBLIC_VERSION ) THEN
         SIGOFFQ(:,:,:,:) = 0D0
         SIG0Q_FSEA(:,:,:,:) = 0D0
         SIG1Q_FSEA_NV(:,:,:,:,:) = 0D0
         SIG1Q_FSEA_VC(:,:,:,:,:) = 0D0
         SIG0Q_SURF(:,:,2:NSPINPROJ,:) = 0D0
         SIG1Q_SURF_NV(:,:,2:NSPINPROJ,:,:) = 0D0
         SIG1Q_SURF_VC(:,:,2:NSPINPROJ,:,:) = 0D0
C
         DO I = 1,3
            DO J = 1,3
               IF ( I.EQ.J ) CYCLE
               SIG0Q_SURF(I,J,1,:) = 0D0
               SIG1Q_SURF_NV(I,J,1,:,:) = 0D0
               SIG1Q_SURF_VC(I,J,1,:,:) = 0D0
            END DO
         END DO
      END IF
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      ALLOCATE (CSIGD(3,3),RSIGD(3,3),RSIGDL(3,3))
C
      WRITE (STR_I_TEMP_LAT,'(I5.5)') I_TEMP_LAT
      FILNAM = 'sigma'//STR_I_TEMP_LAT//'.dat'
      OPEN (IOTMP4,FILE=FILNAM,STATUS='REPLACE')
      WRITE (IOTMP4,'("# Git-Hash: ",A)') GIT_HASH(1:LEN_TRIM(GIT_HASH))
      WRITE (IOTMP4,'("# Git compile date: ",A)')
     &       GIT_COMPILE_DATE(1:LEN_TRIM(GIT_COMPILE_DATE))
      WRITE (IOTMP4,'("# Git branch: ",A)')
     &       GIT_BRANCH(1:LEN_TRIM(GIT_BRANCH))
      WRITE (IOTMP4,*)
C
      LPRSIG = .FALSE.
C
C-----------------------------------------------------------------------
C               devide by the atomic volume VOL_AC in a.u.
C-----------------------------------------------------------------------
C
      DO ISPR = 1,NSPR
C
         ISPINPROJ = LIST_ISPR(ISPR)
C------- no 1/Volume in torkance & Edelstein prefactor (I: no 1/VOL_AC)
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT .OR. 
     &        ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) CYCLE
C
         DO IQ1 = IQBOT_CHI,IQTOP_CHI
C
            IM = IMQ(IQ1)
C
C--------------------------------------------Contribution from Fermi-sea
            SIG0Q_FSEA(:,:,ISPINPROJ,IQ1)
     &         = SIG0Q_FSEA(:,:,ISPINPROJ,IQ1)/VOL_AC(IM)
C
            SIG1Q_FSEA_NV(:,:,ISPINPROJ,IQ1,:)
     &         = SIG1Q_FSEA_NV(:,:,ISPINPROJ,IQ1,:)/VOL_AC(IM)
C
            SIG1Q_FSEA_VC(:,:,ISPINPROJ,IQ1,:)
     &         = SIG1Q_FSEA_VC(:,:,ISPINPROJ,IQ1,:)/VOL_AC(IM)
C
C---------------------------------------------------Contribution from EF
            SIG0Q_SURF(:,:,ISPINPROJ,IQ1)
     &         = SIG0Q_SURF(:,:,ISPINPROJ,IQ1)/VOL_AC(IM)
C
            SIGOFFQ(:,:,ISPINPROJ,IQ1) = SIGOFFQ(:,:,ISPINPROJ,IQ1)
     &         /VOL_AC(IM)
C
            SIG1Q_SURF_NV(:,:,ISPINPROJ,IQ1,:)
     &         = SIG1Q_SURF_NV(:,:,ISPINPROJ,IQ1,:)/VOL_AC(IM)
C
            SIG1Q_SURF_VC(:,:,ISPINPROJ,IQ1,:)
     &         = SIG1Q_SURF_VC(:,:,ISPINPROJ,IQ1,:)/VOL_AC(IM)
C
         END DO
C
      END DO
C
C
C-----------------------------------------------------------------------
C          collect the perturbation from the various sites IQ2
C              take the average over the sites IQ1
C-----------------------------------------------------------------------
C
C
      DO ISPR = 1,NSPR
C
         ISPINPROJ = LIST_ISPR(ISPR)
C
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
C
C---------- torkance prefactor for SI output (C * m)
            PRESI = SOTSI
C---------- torkance prefactor for a.u. output
            PREAU = SOT_PREFAC_AU
C---------- no 1/Volume in torkance prefactor (II: no 1/NQ_AVERAGE)
            NQ_AVERAGE = 1
C
         ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
C
C---------- Edelstein prefactor for SI output (m/V)
            PRESI = EESI
C---------- Edelstein prefactor for a.u. output
            PREAU = EE_PREFAC_AU
C---------- no 1/Volume in Edelstein prefactor (II: no 1/NQ_AVERAGE)
            NQ_AVERAGE = 1
C
         ELSE
C
C---------- conductivity prefactor for S.I. output
C---------- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
            PRESI = 1D-8*CONSI
C---------- conductivity prefactor for a.u. output
            PREAU = SIG_PREFAC_AU
C---------- number of sites for 1/Volume-prefactor
            NQ_AVERAGE = IQTOP_CHI - IQBOT_CHI + 1
C
         END IF
C
         WRITE (6,99015)
         WRITE (6,'(/,A80)') STR_ISP_PROJ(ISPINPROJ)
         WRITE (6,99015)
C
C--------------------------------------------Contribution from Fermi-sea
         SIG0(:,:) = 0D0
         SIG1_NV(:,:) = 0D0
         SIG1_VC(:,:) = 0D0
C
         DO IQ1 = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
               DO NUE = 1,3
                  SIG0(MUE,NUE) = SIG0(MUE,NUE)
     &                            + SIG0Q_FSEA(MUE,NUE,ISPINPROJ,IQ1)
                  DO IQ2 = IQBOT_CHI,IQTOP_CHI
                     SIG1_NV(MUE,NUE) = SIG1_NV(MUE,NUE)
     &                                  + SIG1Q_FSEA_NV(MUE,NUE,
     &                                  ISPINPROJ,IQ1,IQ2)
                     SIG1_VC(MUE,NUE) = SIG1_VC(MUE,NUE)
     &                                  + SIG1Q_FSEA_VC(MUE,NUE,
     &                                  ISPINPROJ,IQ1,IQ2)
                  END DO
               END DO
            END DO
         END DO
C
         SIG0(:,:) = SIG0(:,:)/DBLE(NQ_AVERAGE)
         SIG1_NV(:,:) = SIG1_NV(:,:)/DBLE(NQ_AVERAGE)
         SIG1_VC(:,:) = SIG1_VC(:,:)/DBLE(NQ_AVERAGE)
C
C---------------------------------------------------Contribution from EF
         SIG0EF(:,:) = 0D0
         SIG1_NVEF(:,:) = 0D0
         SIG1_VCEF(:,:) = 0D0
         SIGOFF(:,:) = 0D0
C
         DO IQ1 = IQBOT_CHI,IQTOP_CHI
            DO MUE = 1,3
               DO NUE = 1,3
                  SIG0EF(MUE,NUE) = SIG0EF(MUE,NUE)
     &                              + SIG0Q_SURF(MUE,NUE,ISPINPROJ,IQ1)
                  SIGOFF(MUE,NUE) = SIGOFF(MUE,NUE)
     &                              + SIGOFFQ(MUE,NUE,ISPINPROJ,IQ1)
                  DO IQ2 = IQBOT_CHI,IQTOP_CHI
                     SIG1_NVEF(MUE,NUE) = SIG1_NVEF(MUE,NUE)
     &                  + SIG1Q_SURF_NV(MUE,NUE,ISPINPROJ,IQ1,IQ2)
                     SIG1_VCEF(MUE,NUE) = SIG1_VCEF(MUE,NUE)
     &                  + SIG1Q_SURF_VC(MUE,NUE,ISPINPROJ,IQ1,IQ2)
                  END DO
               END DO
            END DO
         END DO
C
         SIG0EF(:,:) = SIG0EF(:,:)/DBLE(NQ_AVERAGE)
         SIGOFF(:,:) = SIGOFF(:,:)/DBLE(NQ_AVERAGE)
         SIG1_NVEF(:,:) = SIG1_NVEF(:,:)/DBLE(NQ_AVERAGE)
         SIG1_VCEF(:,:) = SIG1_VCEF(:,:)/DBLE(NQ_AVERAGE)
C
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
C     calculate total conductivity tensor:
C     ------------------------------------
         DO MUE = 1,3
            DO NUE = 1,3
C
               SIGTOT_NV(MUE,NUE) = SIG1_NV(MUE,NUE) + SIG0(MUE,NUE)
     &                              + SIG1_NVEF(MUE,NUE)
     &                              + SIG0EF(MUE,NUE)
               SIGTOT_VC(MUE,NUE) = SIG1_VC(MUE,NUE) + SIG0(MUE,NUE)
     &                              + SIG1_VCEF(MUE,NUE)
     &                              + SIG0EF(MUE,NUE)
C
C           check if tensor elements are all real:
C           --------------------------------------
               IF ( ABS(DIMAG(SIGTOT_NV(MUE,NUE))).GT.1.D-5 )
     &              WRITE (6,99004) 'NVC',MUE,NUE,SIGTOT_NV(MUE,NUE)
               IF ( ABS(DIMAG(SIGTOT_VC(MUE,NUE))).GT.1.D-5 )
     &              WRITE (6,99004) ' VC',MUE,NUE,SIGTOT_VC(MUE,NUE)
C
               SIGMA_NV(MUE,NUE) = DREAL(SIGTOT_NV(MUE,NUE))
               SIGMA_VC(MUE,NUE) = DREAL(SIGTOT_VC(MUE,NUE))
C
            END DO
         END DO
C---------------------------------------------------Contribution from EF
         DO MUE = 1,3
            DO NUE = 1,3
C
               SIGTOT_NVEF(MUE,NUE) = SIG1_NVEF(MUE,NUE)
     &                                + SIG0EF(MUE,NUE)
               SIGTOT_VCEF(MUE,NUE) = SIG1_VCEF(MUE,NUE)
     &                                + SIG0EF(MUE,NUE)
C
C           check if tensor elements are all real:
C           --------------------------------------
               IF ( ABS(DIMAG(SIGTOT_NVEF(MUE,NUE))).GT.1.D-5 )
     &              WRITE (6,99008) 'NVC',MUE,NUE,SIGTOT_NVEF(MUE,NUE)
               IF ( ABS(DIMAG(SIGTOT_VCEF(MUE,NUE))).GT.1.D-5 )
     &              WRITE (6,99008) ' VC',MUE,NUE,SIGTOT_VCEF(MUE,NUE)
C
               IF ( ABS(DIMAG(SIGOFF(MUE,NUE))).GT.1.D-5 )
     &              WRITE (6,99008) 'bru',MUE,NUE,SIGOFF(MUE,NUE)
C
C
               SIGMA_NVEF(MUE,NUE) = DREAL(SIGTOT_NVEF(MUE,NUE))
               SIGMA_VCEF(MUE,NUE) = DREAL(SIGTOT_VCEF(MUE,NUE))
C
            END DO
         END DO
C
C
C
         SIGMA_NVEFSAVE(1:3,1:3) = SIGMA_NVEF(1:3,1:3)
         SIGMA_VCEFSAVE(1:3,1:3) = SIGMA_VCEF(1:3,1:3)
         SIGMA_NVSAVE(1:3,1:3) = SIGMA_NV(1:3,1:3)
         SIGMA_VCSAVE(1:3,1:3) = SIGMA_VC(1:3,1:3)
C
C
         L_SIG_SINGULAR = .FALSE.
         IF ( ABS(RMAT3X3DET(SIGMA_NVEF)).LT.1D-8 )
     &        L_SIG_SINGULAR = .TRUE.
C
         IF ( .NOT.L_SIG_SINGULAR ) THEN
            CALL RMATINV(3,3,SIGMA_NVEF,RHO_NVEF)
            CALL RMATINV(3,3,SIGMA_VCEF,RHO_VCEF)
            CALL RMATINV(3,3,SIGMA_NV,RHO_NV)
            CALL RMATINV(3,3,SIGMA_VC,RHO_VC)
         ELSE
            RHO_NVEF = 0D0
            RHO_VCEF = 0D0
            RHO_NV = 0D0
            RHO_VC = 0D0
         END IF
C
C     sum up diagonal elements:
C     -------------------------
         SIG1AU_NV = 0.D0
         SIG1AU_VC = 0.D0
         RHOAU_NV = 0.D0
         RHOAU_VC = 0.D0
         DO MUE = 1,3
            SIG1AU_NV = SIG1AU_NV + SIGMA_NVSAVE(MUE,MUE)/3D0
            SIG1AU_VC = SIG1AU_VC + SIGMA_VCSAVE(MUE,MUE)/3D0
            RHOAU_NV = RHOAU_NV + RHO_NV(MUE,MUE)/3D0
            RHOAU_VC = RHOAU_VC + RHO_VC(MUE,MUE)/3D0
         END DO
         SIGMAAU_NV = SIG1AU_NV
         SIGMAAU_VC = SIG1AU_VC
C
C        calculate resistivities in SI-units
C        -----------------------------------
         RESIST_NV = RHOAU_NV/PRESI
         RESIST_VC = RHOAU_VC/PRESI
         IF ( RESIST_NV.LT.0 ) WRITE (6,99006)
         IF ( RESIST_VC.LT.0 ) WRITE (6,99006)
C
         IF ( ISPR.EQ.1 ) DIFF = 100.D0*(SIGMAAU_VC-SIGMAAU_NV)
     &                           /SIGMAAU_NV
C
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
            RESPSTR = 'Torkance'
            WRITE (6,99009) TRIM(RESPSTR)
         ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
            RESPSTR = 'Edelstein response'
            WRITE (6,99009) TRIM(RESPSTR)
         ELSE
            RESPSTR = 'Sigma'
            WRITE (6,99009) 'Conductivity'
         END IF
         WRITE (6,99012)
         DO MUE = 1,3
            DO NUE = 1,3
               WRITE (6,99011) MUE,NUE,SIGMA_NVEFSAVE(MUE,NUE)*PREAU,
     &                         SIGMA_NVSAVE(MUE,NUE)*PREAU,
     &                         SIGMA_NVSAVE(MUE,NUE)*PRESI,
     &                         SIGMA_VCEFSAVE(MUE,NUE)*PREAU,
     &                         SIGMA_VCEFSAVE(MUE,NUE)*PRESI,
     &                         SIGMA_VCSAVE(MUE,NUE)*PREAU,
     &                         SIGMA_VCSAVE(MUE,NUE)*PRESI
            END DO
         END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
         IF ( WRBUILDBOT ) WRITE (IFILBUILDBOT,99021)
     &                            ROUTINE(1:LEN_TRIM(ROUTINE)),ISPR,
     &                            ISPINPROJ,STR_ISP_PROJ(ISPINPROJ)
     &                            (1:LEN_TRIM(STR_ISP_PROJ(ISPINPROJ))),
     &                            ((SIGMA_NVEFSAVE(MUE,NUE)*PREAU,
     &                            SIGMA_NVSAVE(MUE,NUE)*PREAU,
     &                            SIGMA_VCEFSAVE(MUE,NUE)*PREAU,
     &                            SIGMA_VCSAVE(MUE,NUE)*PREAU,NUE=1,3),
     &                            MUE=1,3)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C                  print sigma tensor w.r.t. to LOCAL frame
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
C------------------- if there is a global rotation of the magnetisation:
C----------------------------------- give the tensors in the local frame
C------------------------------------ use rotation matrix MROTQ for IQ=1
C
C       NOTE:      ->z_global = MROT_frame * ->z_special
C
         IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
C
            SIGMA_NVSAVE_SI(:,:) = 0D0
            SIGMA_VCSAVE_SI(:,:) = 0D0
C
            SIGMA_NVSAVE_SI = SIGMA_NVSAVE*PRESI
            SIGMA_VCSAVE_SI = SIGMA_VCSAVE*PRESI
C
            RHO_NV_SI(:,:) = 0D0
            RHO_VC_SI(:,:) = 0D0
C
            RHO_NV_SI = RHO_NV/PRESI
            RHO_VC_SI = RHO_VC/PRESI
C
            IF ( KMROT.EQ.2 .AND. .NOT.CHANGE_FRAME ) MROT_FRAME(:,:)
     &           = MROTQ(:,:,1)
C
            CALL DGEMM('N','T',3,3,3,1D0,SIGMA_NVSAVE_SI,3,MROT_FRAME,3,
     &                 0D0,MTMP,3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                 SIGMAL_NV,3)
C
            CALL DGEMM('N','T',3,3,3,1D0,SIGMA_VCSAVE_SI,3,MROT_FRAME,3,
     &                 0D0,MTMP,3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                 SIGMAL_VC,3)
C
            CALL DGEMM('N','T',3,3,3,1D0,RHO_NV_SI,3,MROT_FRAME,3,0D0,
     &                 MTMP,3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                 RHOL_NV,3)
C
            CALL DGEMM('N','T',3,3,3,1D0,RHO_VC_SI,3,MROT_FRAME,3,0D0,
     &                 MTMP,3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                 RHOL_VC,3)
C
            WRITE (6,99007) KMROT,CHANGE_FRAME,MROT_FRAME,SIGMAL_NV,
     &                      SIGMAL_VC,RHOL_NV,RHOL_VC
         END IF
Cfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
C
         IF ( ISPR.EQ.1 ) THEN
            SMA(1) = 100D0*(1.D8/SIGMA_VCSAVE(3,3)-1.D8/SIGMA_VCSAVE(1,1
     &               ))/(1D0/3D0*1.D8/SIGMA_VCSAVE(3,3)+2D0/3D0*1.D8/
     &               SIGMA_VCSAVE(1,1))
C
            WRITE (6,99013) DREAL(SMA(1))
         END IF
C
         IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
            WRITE (6,'(/,"a.u. (Rydberg):  [e*a_0]")')
            WRITE (6,'(/,"SI units      :  [10**-30 C*m]")')
         ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
            WRITE (6,'(/,"a.u. (Rydberg):  [a_0**2/e]")')
            WRITE (6,'(/,"SI units      :  [m/V]")')
         ELSE
            WRITE (6,'(/,"a.u. (Rydberg):  [(e**2)/(2*hbar*a_0)]")')
            WRITE (6,'(/,"SI units      :  [1/(muOhm*cm)]")')
         END IF
         WRITE (6,'(/,/,"In tensor form ( SI units )")')
C
         DUMMAT(:,:) = SIGMA_VCSAVE(:,:)*PRESI
         CALL SIGSTR(RESPSTR(1:3)//'_TOT_VC  ',12,DUMMAT,1D-8,6)
C
         IF ( NEA.GT.1 ) THEN
            DUMMAT(:,:) = (SIGMA_VCSAVE(:,:)-SIGMA_VCEFSAVE(:,:))*PRESI
            CALL SIGSTR(RESPSTR(1:3)//'_FSEA_VC ',12,DUMMAT,1D-8,6)
C
            DUMMAT(:,:) = SIGMA_VCEFSAVE(:,:)*PRESI
            CALL SIGSTR(RESPSTR(1:3)//'_FSURF_VC',12,DUMMAT,1D-8,6)
         END IF
C
         IQB = IQBOT_CHI
         IQT = IQTOP_CHI
         PRESI_NQ = PRESI/DBLE(NQ_AVERAGE)
         PREAU_NQ = PREAU/DBLE(NQ_AVERAGE)
C ----------------------------------------------------- print SIG0 FSURF
         CSIGD = SUM(SIG0Q_SURF(:,:,ISPINPROJ,IQB:IQT),DIM=3)
         RSIGD = DREAL(CSIGD)*PRESI_NQ
         DS1 = RESPSTR(1:3)//'0_FSURF'
         CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
         IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
            CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,MTMP,
     &                 3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,RSIGDL,
     &                 3)
            CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
         END IF
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',DS1)
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',DS1)
C -------------------------------------------------- print SIG1_NV FSURF
         CSIGD = SUM(SUM(SIG1Q_SURF_NV(:,:,ISPINPROJ,IQB:IQT,IQB:IQT),
     &           DIM=4),DIM=3)
         RSIGD = DREAL(CSIGD)*PRESI_NQ
         DS1 = RESPSTR(1:3)//'1_NV_FSURF'
         CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
         IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
            CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,MTMP,
     &                 3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,RSIGDL,
     &                 3)
            CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
         END IF
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',DS1)
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',DS1)
C -------------------------------------------------- print SIG1_VC FSURF
         CSIGD = SUM(SUM(SIG1Q_SURF_VC(:,:,ISPINPROJ,IQB:IQT,IQB:IQT),
     &           DIM=4),DIM=3)
         RSIGD = DREAL(CSIGD)*PRESI_NQ
         DS1 = RESPSTR(1:3)//'1_VC_FSURF'
         CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
         IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
            CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,MTMP,
     &                 3)
            CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,RSIGDL,
     &                 3)
            CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
         END IF
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',DS1)
         CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',DS1)
C
         IF ( NEA.GT.1 ) THEN
C ------------------------------------------------------ print SIG0 FSEA
            CSIGD = SUM(SIG0Q_FSEA(:,:,ISPINPROJ,IQB:IQT),DIM=3)
            RSIGD = DREAL(CSIGD)*PRESI_NQ
            DS1 = RESPSTR(1:3)//'0_FSEA'
            CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
            IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
               CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,
     &                    MTMP,3)
               CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                    RSIGDL,3)
               CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
            END IF
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',
     &                            DS1)
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',
     &                            DS1)
C --------------------------------------------------- print SIG1_NV FSEA
            CSIGD = SUM(SUM(SIG1Q_FSEA_NV(:,:,ISPINPROJ,IQB:IQT,IQB:IQT)
     &              ,DIM=4),DIM=3)
            RSIGD = DREAL(CSIGD)*PRESI_NQ
            DS1 = RESPSTR(1:3)//'1_NV_FSEA'
            CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
            IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
               CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,
     &                    MTMP,3)
               CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                    RSIGDL,3)
               CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
            END IF
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',
     &                            DS1)
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',
     &                            DS1)
C --------------------------------------------------- print SIG1_VC FSEA
            CSIGD = SUM(SUM(SIG1Q_FSEA_VC(:,:,ISPINPROJ,IQB:IQT,IQB:IQT)
     &              ,DIM=4),DIM=3)
            RSIGD = DREAL(CSIGD)*PRESI_NQ
            DS1 = RESPSTR(1:3)//'1_VC_FSEA'
            CALL SIGSTR(DS1,LEN_TRIM(DS1),RSIGD,1D-8,6)
            IF ( CHANGE_FRAME .OR. KMROT.EQ.2 ) THEN
               CALL DGEMM('N','T',3,3,3,1D0,RSIGD,3,MROT_FRAME,3,0D0,
     &                    MTMP,3)
               CALL DGEMM('N','N',3,3,3,1D0,MROT_FRAME,3,MTMP,3,0D0,
     &                    RSIGDL,3)
               CALL SIGSTR('in local FOR',12,RSIGDL,1D-8,6)
            END IF
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PREAU_NQ,ISPINPROJ,'au',
     &                            DS1)
            CALL DUMP_COND_TENSOR(IOTMP4,CSIGD,PRESI_NQ,ISPINPROJ,'SI',
     &                            DS1)
         END IF
C
C
         IF ( .NOT.L_SIG_SINGULAR .AND. ISPINPROJ.EQ.IRESPONSE_ADA_ADA )
     &        THEN
            WRITE (6,*) '   '
            WRITE (6,99010)
            WRITE (6,99012)
            DO MUE = 1,3
               DO NUE = 1,3
                  WRITE (6,99011) MUE,NUE,RHO_NVEF(MUE,NUE)/PREAU,
     &                            RHO_NV(MUE,NUE)/PREAU,RHO_NV(MUE,NUE)
     &                            /PRESI,RHO_VCEF(MUE,NUE)/PREAU,
     &                            RHO_VCEF(MUE,NUE)/PRESI,
     &                            RHO_VC(MUE,NUE)/PREAU,RHO_VC(MUE,NUE)
     &                            /PRESI
               END DO
            END DO
C
C
            IF ( ISPR.EQ.1 ) THEN
               SMA(2) = 100D0*(RHO_VC(3,3)-RHO_VC(1,1))
     &                  /(1D0/3D0*RHO_VC(3,3)+2D0/3D0*RHO_VC(1,1))
C
               WRITE (6,99014) DREAL(SMA(2))
            END IF
C
            WRITE (6,'(/,"a.u. (Rydberg):  [(hbar*a_0)*2/(e**2)]")')
            WRITE (6,'(/,"SI units      :  [muOhm*cm]")')
C
         END IF
C
C
         IF ( IPRINT.NE.0 .AND. ISPR.EQ.1 ) THEN
            WRITE (6,99001)
            WRITE (6,99005) SIGMAAU_NV*PREAU,SIGMAAU_NV*CONSI,
     &                      1.D8/(SIGMAAU_NV*CONSI),RESIST_NV,
     &                      ' (no VC)  '
            WRITE (6,99002)
            WRITE (6,99005) SIGMAAU_VC*PREAU,SIGMAAU_VC*CONSI,
     &                      1.D8/(SIGMAAU_VC*CONSI),RESIST_VC,
     &                      ' (with VC)'
            WRITE (6,99003) DIFF
            IF ( DIFF.LT.0D0 ) WRITE (6,*)
     &                                 'THIS IS STRANGE,PLEASE CHECK!!'
         END IF
C
C------------------------------------------ write conductivities to file
         IF ( ISPINPROJ.EQ.IRESPONSE_ADA_ADA .AND. 
     &        SIG_MODE.NE.'INTEGRAL  ' .AND. LPRSIG .AND. NQ.EQ.1 ) THEN
            DO MUE = 1,3
               DO NUE = 1,3
                  WRITE (MSTR(1:1),'(i1)') MUE
                  WRITE (NSTR(1:1),'(i1)') NUE
                  FILNAM = DATSET(1:LDATSET)//MSTR//NSTR//'_SIG1.dat'
                  OPEN (IOTMP,FILE=FILNAM,STATUS='replace')
                  WRITE (IOTMP,99016)
C
                  WRITE (IOTMP,99017) X_VFT(1),
     &                                PRESI*DREAL(SIG1_NVEF(MUE,NUE)
     &                                +SIG0EF(MUE,NUE)),
     &                                PRESI*DREAL(SIG1_VCEF(MUE,NUE)
     &                                +SIG0EF(MUE,NUE)),
     &                                PRESI*DREAL(SIG1_NVEF(MUE,NUE)),
     &                                PRESI*DREAL(SIG1_VCEF(MUE,NUE))
                  CLOSE (IOTMP)
               END DO ! NUE
            END DO ! MUE
         END IF
C-----------------------------------------------------------------------
      END DO ! ISPR
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C ======================================================================
C        write layer resolved conductivity in case of layered system
C ======================================================================
      IF ( SYSTEM_TYPE(1:3).EQ.'LIR' .OR. SYSTEM_TYPE(1:3).EQ.'LIV' .OR. 
     &     SYSTEM_TYPE(1:3).EQ.'VIV' ) THEN
C
         DO IV = 1,2
C
            DO ISPR = 1,NSPR
               ISPINPROJ = LIST_ISPR(ISPR)
C
               IF ( ISPINPROJ.EQ.IRESPONSE_SOT ) THEN
C---------------- torkance prefactor for SI output (C * m)
                  PRESI = SOTSI
               ELSE IF ( ISPINPROJ.EQ.IRESPONSE_EDELSTEIN ) THEN
C---------------- Edelstein prefactor for SI output (m/V)
                  PRESI = EESI
               ELSE
C----- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
                  PRESI = 1D-8*CONSI
               END IF
C
               WRITE (STR_ISPR,'(I2.2)') ISPINPROJ
C
               IF ( IV.EQ.1 ) THEN
                  LFILNAM = 'sigma'//STR_I_TEMP_LAT//'_layer_ISPR'//
     &                      STR_ISPR//'_NV.dat'
                  OPEN (IOTMP,FILE=LFILNAM)
               ELSE
                  LFILNAM = 'sigma'//STR_I_TEMP_LAT//'_layer_ISPR'//
     &                      STR_ISPR//'_VC.dat'
                  OPEN (IOTMP,FILE=LFILNAM)
               END IF
C
               WRITE (IOTMP,99019) IQBOT_CHI,IQTOP_CHI
               WRITE (IOTMP,99018) 0,(((PRESI*DREAL(SIG0Q_SURF(MUE,NUE,
     &                             ISPINPROJ,IQ1)),IQ1=IQBOT_CHI,
     &                             IQTOP_CHI),NUE=1,3),MUE=1,3)
C
               WRITE (IOTMP,99020) IQBOT_CHI,IQTOP_CHI
               DO IQ2 = IQBOT_CHI,IQTOP_CHI
C
                  IF ( IV.EQ.1 ) THEN
                     WRITE (IOTMP,99018) IQ2,
     &                      (((PRESI*DREAL(SIG1Q_SURF_NV(MUE,NUE,
     &                      ISPINPROJ,IQ1,IQ2)),IQ1=IQBOT_CHI,IQTOP_CHI)
     &                      ,NUE=1,3),MUE=1,3)
                  ELSE
                     WRITE (IOTMP,99018) IQ2,
     &                      (((PRESI*DREAL(SIG1Q_SURF_VC(MUE,NUE,
     &                      ISPINPROJ,IQ1,IQ2)),IQ1=IQBOT_CHI,IQTOP_CHI)
     &                      ,NUE=1,3),MUE=1,3)
                  END IF
               END DO
C
               CLOSE (IOTMP)
C
            END DO
C
         END DO
C
      END IF
C
      CLOSE (IOTMP4)
C
      IF ( PUBLIC_VERSION ) WRITE (8,99022)
C
C=======================================================================
99001 FORMAT (/,'Results without Vertex-corrections:',/,35('='))
99002 FORMAT (/,'Results including Vertex-corrections:',/,37('='))
99003 FORMAT ('Inclusion of vertex-correction increases',
     &        ' conductivity by ',f8.4,'%',/)
99004 FORMAT ('Sigma_',a3,' not real!!! mue/nue/sigma=',2I4,2E13.5)
99005 FORMAT (/,'Conductivity in a.u.:          ',f18.5,/,
     &        'Conductivity in [1/(Ohm*m)]    ',f18.5,/,
     &        'Inverse Conductivity [muOhm.cm]',f18.5,/,
     &        'Resistivity in [muOhm.cm]      ',f18.5,a10)
99006 FORMAT (/,45('!'),/,' Negative resistivity...bad luck... ',
     &        'try again',/,45('!'))
99007 FORMAT (/,5X,'galvano-magnetic tensors for the LOCAL frame ',
     &        'for  KMROT =',I2,' CHANGE_FRAME = ',L1,/,5X,
     &        'with rotation matrix:',//,
     &        3(25X,'(',F8.3,',',F8.3,',',F8.3,'  )',/),//,5X,
     &        'SIGMA  (NV)        ',3F14.8,/,5X,'[1/(muOhm*cm)]     ',
     &        3F14.8,/,24X,3F14.8,//,5X,'SIGMA  (VC)        ',3F14.8,/,
     &        5X,'[1/(muOhm*cm)]     ',3F14.8,/,24X,3F14.8,//,5X,
     &        'RHO    (NV)        ',3F14.8,/,5X,'[muOhm*cm]         ',
     &        3F14.8,/,24X,3F14.8,//,5X,'RHO    (VC)        ',3F14.8,/,
     &        5X,'[muOhm*cm]         ',3F14.8,/,24X,3F14.8,//)
99008 FORMAT ('Sigma_',a3,
     &        ' (Contribution from EF) not real!!!  mue/nue/sigma=',2I4,
     &        2E13.5)
99009 FORMAT (A,' (no Streda)',/,28('-'),/,25(' '),'without Vertex',
     &        43(' '),'with Vertex')
99010 FORMAT ('Resistivity (no Streda)',/,33('-'),/,25(' '),
     &        'without Vertex',43(' '),'with Vertex')
99011 FORMAT (I2,I4,f15.8,f15.8,f14.8,'     ',f15.8,'  ',f12.8,'  ',
     &        f15.8,'  ',f12.8)
99012 FORMAT ('MUE NUE        EF                  Complete',22(' '),
     &        '        EF                        Complete',/,14(' '),
     &        'a.u.           a.u.          SI    ',14(' '),
     &        'a.u.          SI               a.u.          SI     ')
99013 FORMAT ('SMA from inv. cond. (with VC)                    [%]:',
     &        3(' '),f10.6)
99014 FORMAT ('SMA from resistivity (with VC)                    [%]:',
     &        3(' '),f10.6)
99015 FORMAT (/,40('*'),/,40('*'),/)
99016 FORMAT ('# conc',12x,'signv',12x,'sigvc',12x,'sig1nv',12x,
     &        'sig1vc',12x,'sigcrep')
99017 FORMAT (99E17.8)
99018 FORMAT (I3,999E13.5)
99019 FORMAT ('#  sigma0_I   for layers I=',I3,'  to  ',I3,
     &        '     first columns xx, then xy, xz... and last zz',/,
     &        '#  dummy index 0 in first column indicates',/,
     &        '#  layer diagonal term')
99020 FORMAT ('#  sigma1_IJ  for layers I=',I3,'  to  ',I3,
     &        '     first columns xx, then xy, xz... and last zz',/,
     &        '#  1 line per neighboring layer J ',
     &        '(given in first column) ')
99021 FORMAT ('# BUILDBOT: ',A,':  SIG tensor (a.u.) for  ISPR =',I2,
     &        '  ISPINPROJ =',I2,2X,A,/,(1PE22.14))
99022 FORMAT (' ',10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &        '  PUBLIC VERSION does not supply:              *',/,10X,
     &        '               - off-diagonal condicutivity    *',/,10X,
     &        '               - spin conductivity etc.        *',/,10X,
     &        '*',60X,'*',/,10X,62('*'),//)
      END
C*==sigstr.f    processed by SPAG 6.70Rc at 13:20 on 21 Apr 2017
      SUBROUTINE SIGSTR(STR,LSTR,A,TOLP,K_FMT_FIL)
C   ********************************************************************
C   *                                                                  *
C   *   writes structure of   REAL    3x3   matrix   A                 *
C   *                                                                  *
C   *   TOL         tolerance for difference                           *
C   *   K_FMT_FIL   output channel                                     *
C   *               a negative sign suppresses table at the end        *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGSTR')
      INTEGER MCOL,MLIN,N,M
      PARAMETER (MCOL=0,MLIN=0,N=3,M=3)
C
C Dummy arguments
C
      INTEGER K_FMT_FIL,LSTR
      CHARACTER*(*) STR
      REAL*8 TOLP
      REAL*8 A(M,M)
C
C Local variables
C
      REAL*8 ARG,B(:,:),CA,CB,DTAB(:),TOL
      CHARACTER*1 CTAB(:),VZ(-1:+1)
      LOGICAL DEGENERATE,DIAGONAL,SYMMETRIC
      CHARACTER*250 FMT1,FMT2,FMT3,FMT4
      INTEGER I,I1,IA_ERR,IC0,ID,IFIL,IL,ILSEP(20),IPT(218),IQ,ISL,
     &        IW(:),J,J0,JP,JQ,K,L3,LF,MM,N1,N2,N3,NC,ND,NK,NM,NM1,NM2,
     &        NM3,NNON0,NSL
      LOGICAL SAME,SMALL
C
C*** End of declarations rewritten by SPAG
C
      DATA VZ/'-',' ',' '/
C
      ALLOCATABLE B,CTAB,DTAB,IW
C
      SMALL(ARG) = ABS(ARG*TOL).LT.1.0D0
C
      SAME(CA,CB) = SMALL(1.0D0-CA/CB)
C
      ALLOCATE (B(N,N),CTAB(0:N*N),DTAB(0:N*N),IW(M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'allocate B')
C
      DIAGONAL = .TRUE.
      SYMMETRIC = .TRUE.
      DEGENERATE = .TRUE.
C
      IFIL = ABS(K_FMT_FIL)
C
      TOL = 1.0D0/TOLP
C
      IQ = 1
      JQ = 1
C
C----------------------------------------------------- copy matrix block
C
      J0 = N*(JQ-1)
      DO J = 1,N
         I1 = N*(IQ-1) + 1
         JP = J0 + J
         CALL DCOPY(N,A(I1,JP),1,B(1,J),1)
      END DO
C
C------------------------------------------------ set up character table
C
      NC = 0
      DO I = 1,26
         NC = NC + 1
         IPT(NC) = 62 + I
      END DO
      DO I = 1,8
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 10,26
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 191,218
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 35,38
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 40,42
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 91,93
         NC = NC + 1
         IPT(NC) = I
      END DO
C
C---------------------------------------------------------------- header
      IC0 = ICHAR('0')
      N3 = N/100
      N2 = N/10 - N3*10
      N1 = N - N2*10 - N3*100
C
      IF ( N.LE.18 ) THEN
         FMT1 = '(8X,I3,''|'','
         FMT2 = '( 9X,''--+'','
         FMT3 = '( 9X,'' #|'','
         FMT4 = '( 9X,''  |'','
      ELSE
         FMT1 = '(   I4,''|'','
         FMT2 = '( 2X,''--|'','
         FMT3 = '( 2X,'' #|'','
         FMT4 = '( 2X,''  |'','
      END IF
C
      LF = 11
      L3 = 11
      IF ( MCOL.EQ.0 ) THEN
         FMT1 = FMT1(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'( 2A1),''|'',I3)'
         FMT2 = FMT2(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'(''--''),''+'',I3)'
         FMT3 = FMT3(1:LF)//'60(2X,I2))'
         FMT4 = FMT4(1:LF)//'60(I2,2X))'
         LF = 21
      ELSE
         IF ( MCOL.EQ.1 ) THEN
            NK = NINT(SQRT(DBLE(N)))
         ELSE IF ( MCOL.EQ.2 ) THEN
            NK = NINT(SQRT(DBLE(N/2)))
         ELSE IF ( MCOL.EQ.3 ) THEN
            NK = 2*NINT(SQRT(DBLE(N/2))) - 1
         END IF
         DO K = 1,NK
            IF ( MCOL.LE.2 ) THEN
               NM = 2*K - 1
            ELSE
               NM = 2*((K+1)/2)
            END IF
            NM2 = NM/10
            NM1 = NM - NM2*10
            NM3 = NM/2
            FMT1 = FMT1(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'( 2A1),''|'','
            FMT2 = FMT2(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'(''--''),''|'','
C
            IF ( MCOL.LE.2 ) THEN
               DO MM = 1,NM
                  IF ( MOD(MM,2).EQ.MOD(K,2) ) THEN
                     FMT3 = FMT3(1:L3)//'2X,'
                     FMT4 = FMT4(1:L3)//'I2,'
                  ELSE
                     FMT3 = FMT3(1:L3)//'I2,'
                     FMT4 = FMT4(1:L3)//'2X,'
                  END IF
                  L3 = L3 + 3
               END DO
               FMT3 = FMT3(1:L3)//'''|'','
               FMT4 = FMT4(1:L3)//'''|'','
               L3 = L3 + 4
            ELSE
               FMT3 = FMT3(1:LF)//CHAR(IC0+NM3)//'(2X,I2),''|'','
               FMT4 = FMT4(1:LF)//CHAR(IC0+NM3)//'(I2,2X),''|'','
               L3 = L3 + 13
            END IF
            LF = LF + 13
         END DO
         IF ( MCOL.EQ.2 ) THEN
            FMT1 = FMT1(1:LF)//FMT1(12:LF)
            FMT2 = FMT2(1:LF)//FMT2(12:LF)
C
            FMT3 = FMT3(1:L3)//FMT4(12:L3)
            FMT4 = FMT4(1:L3)//FMT3(12:L3)
            LF = 2*LF - 11
            L3 = 2*L3 - 11
         END IF
         FMT1 = FMT1(1:LF)//'I3)'
         FMT2 = FMT2(1:LF)//'I3)'
         FMT3 = FMT3(1:L3)//'I3)'
         FMT4 = FMT4(1:L3)//'I3)'
      END IF
      IF ( MLIN.EQ.0 ) THEN
         NSL = 1
         ILSEP(1) = N
      ELSE IF ( MLIN.EQ.1 ) THEN
         NSL = NINT(SQRT(DBLE(N)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
      ELSE IF ( MLIN.EQ.2 ) THEN
         NSL = NINT(SQRT(DBLE(N/2)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
         DO IL = 1,NSL
            ILSEP(NSL+IL) = ILSEP(NSL) + IL**2
         END DO
         NSL = 2*NSL
      ELSE IF ( MLIN.EQ.3 ) THEN
         NSL = 2*NINT(SQRT(DBLE(N/2))) - 1
         ILSEP(1) = 2
         DO K = 2,NSL
            ILSEP(K) = ILSEP(K-1) + 2*((K+1)/2)
         END DO
      END IF
C
      IF ( L3.GT.250 ) CALL STOP_MESSAGE(ROUTINE,'L3 > 250')
C
      WRITE (IFIL,99001) STR(1:LSTR)
      WRITE (IFIL,'("           | 1 2 3|    ")')
      WRITE (IFIL,FMT=FMT2)
C------------------------------------------------------------ header end
      NNON0 = 0
      ND = 0
      CTAB(0) = ' '
      DTAB(0) = 9999D0
C
      DO I = 1,N
         DO J = 1,N
            IF ( .NOT.SMALL(B(I,J)) ) THEN
C
               IF ( I.NE.J ) DIAGONAL = .FALSE.
C
               IF ( SYMMETRIC .AND. I.GT.J ) THEN
                  IF ( .NOT.SAME(B(J,I),B(I,J)) ) SYMMETRIC = .FALSE.
               END IF
C
               NNON0 = NNON0 + 1
               DO ID = 1,ND
                  IF ( SAME(B(I,J),+DTAB(ID)) ) THEN
                     IW(J) = +ID
                     GOTO 50
                  END IF
                  IF ( SAME(B(I,J),-DTAB(ID)) ) THEN
                     IW(J) = -ID
                     GOTO 50
                  END IF
               END DO
C----------------------------------------------------------- new element
               ND = ND + 1
               IW(J) = ND
               DTAB(ND) = B(I,J)
               IF ( ABS(DTAB(ND)-1.0D0)*TOL.LT.1.0D0 ) THEN
                  CTAB(ND) = '1'
               ELSE IF ( ABS(DTAB(ND)+1.0D0)*TOL.LT.1.0D0 ) THEN
                  DTAB(ND) = +1.0D0
                  CTAB(ND) = '1'
                  IW(J) = -ND
               ELSE
                  CTAB(ND) = CHAR(IPT(1+MOD((ND+1),NC)))
               END IF
            ELSE
               IW(J) = 0
            END IF
 50      END DO
C------------------------------------------------------------ write line
         WRITE (IFIL,FMT=FMT1) I,
     &                         (VZ(ISIGN(1,IW(J))),CTAB(ABS(IW(J))),J=1,
     &                         N)
C
         DO ISL = 1,NSL
            IF ( I.EQ.ILSEP(ISL) ) WRITE (IFIL,FMT=FMT2)
         END DO
      END DO
C
C------------------------------------------------------------------ foot
C
C      WRITE (IFIL,'("           | 1 2 3|    ")')
C
      IF ( K_FMT_FIL.GT.0 ) THEN
         WRITE (IFIL,99002) TOLP,(ID,CTAB(ID),DTAB(ID),ID=1,ND)
C         IF ( ND.NE.0 ) WRITE (IFIL,99003) NNON0,TOLP,N*N - NNON0,TOLP
         IF ( ND.NE.1 ) THEN
            DEGENERATE = .FALSE.
         ELSE
            DO J = 1,N
               IF ( ISIGN(1,IW(J)).NE.ISIGN(1,IW(1)) )
     &              DEGENERATE = .FALSE.
            END DO
         END IF
C
         IF ( .NOT.(DIAGONAL) ) THEN
C
            IF ( SYMMETRIC ) WRITE (IFIL,99003) 'SYMMETRIC '
            IF ( DEGENERATE ) WRITE (IFIL,99003) 'DEGENERATE'
C
         ELSE IF ( ND.EQ.0 ) THEN
C
            WRITE (IFIL,99004) '0-matrix  ',TOLP
C
         ELSE IF ( ND.EQ.1 .AND. ABS(DTAB(1)-1D0).LT.TOLP ) THEN
C
            WRITE (IFIL,99003) '1-matrix  '
C
         ELSE
C
            WRITE (IFIL,99003) 'DIAGONAL  '
            IF ( DEGENERATE ) WRITE (IFIL,99003) 'DEGENERATE'
C
         END IF
C
         WRITE (IFIL,*) ' '
      ELSE
         WRITE (IFIL,*) ' '
      END IF
C
      DEALLOCATE (B,CTAB,DTAB,IW,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC B')
C
99001 FORMAT (/,8X,A,/)
99002 FORMAT (/,8X,'Elements larger TOL=',E9.1,/,(8X,I3,3X,A1,2X,F20.12)
     &        )
99003 FORMAT (8X,'the matrix is  ',A,'  (within TOL)')
99004 FORMAT (8X,'the matrix is  ',A,'  (within TOL = ',E9.1,' )')
      END
C*==dump_cond_tensor.f    processed by SPAG 6.70Rc at 13:20 on 21 Apr 2017
      SUBROUTINE DUMP_COND_TENSOR(OUT,TENSOR,SCALEFAC,ISP,STRUNIT,DS1)
      USE MOD_SIG,ONLY:STR_ISP_PROJ
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) DS1
      INTEGER ISP,OUT
      REAL*8 SCALEFAC
      CHARACTER*2 STRUNIT
      COMPLEX*16 TENSOR(3,3)
C
C Local variables
C
      CHARACTER*255 DS2,DS3
      CHARACTER*100 FMT1
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DS2 = '# '//' '//ADJUSTL(STR_ISP_PROJ(ISP))
      DS2 = DS2(1:LEN_TRIM(DS2))//'_'//DS1(1:LEN_TRIM(DS1))
      DS3 = DS2(1:LEN_TRIM(DS2))//'_'//STRUNIT
      WRITE (OUT,'(A)') DS3(1:LEN_TRIM(DS3))
C
      FMT1 = '(2I2,2E20.12)'
C
      DO I = 1,3
         DO J = 1,3
            WRITE (OUT,FMT1) I,J,TENSOR(I,J)*SCALEFAC
         END DO
      END DO
C
      WRITE (OUT,*)
C
      END
