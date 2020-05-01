C*==sigsum_optics.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SIGSUM_OPTICS(SIG0Q_OPT,SIG1Q_OPT_NV,SIG1Q_OPT_VC,
     &                         NSPINPROJ,NOMEGA,OMEGATAB,I_TEMP_LAT)
C   ********************************************************************
C   *                                                                  *
C   *   sum up results and calculate final conductivity  tensor        *
C   *                                                                  *
C   *   NOTE: JB gives HBAR etc. in SI units                           *
C   *                                                                  *
C   * 06/10/99  HE  based on JB's routine                              *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:C0,RY_EV
      USE MOD_FILES,ONLY:IOTMP
      USE MOD_SITES,ONLY:NQMAX,IQBOT_CHI,IQTOP_CHI
      USE MOD_SIG,ONLY:LIST_ISPR,NSPR,CONSI
      IMPLICIT NONE
C*--SIGSUM_OPTICS18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIGSUM_OPTICS')
C
C Dummy arguments
C
      INTEGER I_TEMP_LAT,NOMEGA,NSPINPROJ
      REAL*8 OMEGATAB(NOMEGA)
      COMPLEX*16 SIG0Q_OPT(3,3,NSPINPROJ,NQMAX,NOMEGA),
     &           SIG1Q_OPT_NV(3,3,NSPINPROJ,NQMAX,NQMAX,NOMEGA),
     &           SIG1Q_OPT_VC(3,3,NSPINPROJ,NQMAX,NQMAX,NOMEGA)
C
C Local variables
C
      CHARACTER*80 FILNAM
      INTEGER I,IOBSE,IOM,IQ1,IQ2,ISPR
      REAL*8 PRESI
      COMPLEX*16 SIG0_OPT(:,:,:,:),SIG1_OPT_NV(:,:,:,:),
     &           SIG1_OPT_VC(:,:,:,:),SIGTOT_OPT_NV(:,:,:,:),
     &           SIGTOT_OPT_VC(:,:,:,:)
      CHARACTER*5 STR_I_TEMP_LAT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIG0_OPT,SIG1_OPT_NV,SIG1_OPT_VC
      ALLOCATABLE SIGTOT_OPT_NV,SIGTOT_OPT_VC
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (SIG0_OPT(3,3,NSPINPROJ,NOMEGA))
      ALLOCATE (SIG1_OPT_NV(3,3,NSPINPROJ,NOMEGA))
      ALLOCATE (SIG1_OPT_VC(3,3,NSPINPROJ,NOMEGA))
      ALLOCATE (SIGTOT_OPT_NV(3,3,NSPINPROJ,NOMEGA))
      ALLOCATE (SIGTOT_OPT_VC(3,3,NSPINPROJ,NOMEGA))
C
C---- multiply by 1d-8 to convert from 1/( Ohm * m ) to 1/(\mu Ohm * cm)
C
      PRESI = 1D-8*CONSI
C
      SIG0Q_OPT(:,:,:,:,:) = PRESI*SIG0Q_OPT(:,:,:,:,:)
      SIG1Q_OPT_NV(:,:,:,:,:,:) = PRESI*SIG1Q_OPT_NV(:,:,:,:,:,:)
      SIG1Q_OPT_VC(:,:,:,:,:,:) = PRESI*SIG1Q_OPT_VC(:,:,:,:,:,:)
C
      SIG0_OPT(:,:,:,:) = C0
      SIG1_OPT_NV(:,:,:,:) = C0
      SIG1_OPT_VC(:,:,:,:) = C0
C
      SIGTOT_OPT_NV(:,:,:,:) = C0
      SIGTOT_OPT_VC(:,:,:,:) = C0
C
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      LOOP_ISPR:DO ISPR = 1,NSPR
C
         IOBSE = LIST_ISPR(ISPR)
C
         LOOP_IOM:DO IOM = 1,NOMEGA
C
C***********************************************************************
            LOOP_IQ1:DO IQ1 = IQBOT_CHI,IQTOP_CHI
C
               SIG0_OPT(:,:,IOBSE,IOM) = SIG0_OPT(:,:,IOBSE,IOM)
     &            + SIG0Q_OPT(:,:,IOBSE,IQ1,IOM)
C
               LOOP_IQ2:DO IQ2 = IQBOT_CHI,IQTOP_CHI
C
                  SIG1_OPT_NV(:,:,IOBSE,IOM)
     &               = SIG1_OPT_NV(:,:,IOBSE,IOM)
     &               + SIG1Q_OPT_NV(:,:,IOBSE,IQ1,IQ2,IOM)
C
                  SIG1_OPT_VC(:,:,IOBSE,IOM)
     &               = SIG1_OPT_VC(:,:,IOBSE,IOM)
     &               + SIG1Q_OPT_VC(:,:,IOBSE,IQ1,IQ2,IOM)
C
               END DO LOOP_IQ2
            END DO LOOP_IQ1
C***********************************************************************
C
            SIGTOT_OPT_NV(:,:,IOBSE,IOM) = SIGTOT_OPT_NV(:,:,IOBSE,IOM)
     &         + SIG0_OPT(:,:,IOBSE,IOM) + SIG1_OPT_NV(:,:,IOBSE,IOM)
C
            SIGTOT_OPT_VC(:,:,IOBSE,IOM) = SIGTOT_OPT_VC(:,:,IOBSE,IOM)
     &         + SIG0_OPT(:,:,IOBSE,IOM) + SIG1_OPT_VC(:,:,IOBSE,IOM)
C
         END DO LOOP_IOM
      END DO LOOP_ISPR
COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      WRITE (STR_I_TEMP_LAT,'(I5.5)') I_TEMP_LAT
      FILNAM = 'sigma_NV'//STR_I_TEMP_LAT//'.dat'
      OPEN (IOTMP,FILE=FILNAM,STATUS='REPLACE')
C
      DO IOM = 1,NOMEGA
         WRITE (IOTMP,99001) OMEGATAB(IOM)*RY_EV,
     &                       (SIGTOT_OPT_NV(I,I,1,IOM),I=1,3),
     &                       SIGTOT_OPT_NV(1,2,1,IOM),
     &                       SIGTOT_OPT_NV(2,1,1,IOM)
      END DO
C
      CLOSE (IOTMP)
C
      FILNAM = 'sigma_VC'//STR_I_TEMP_LAT//'.dat'
      OPEN (IOTMP,FILE=FILNAM,STATUS='REPLACE')
C
      DO IOM = 1,NOMEGA
         WRITE (IOTMP,99001) OMEGATAB(IOM)*RY_EV,
     &                       (SIGTOT_OPT_VC(I,I,1,IOM),I=1,3),
     &                       SIGTOT_OPT_VC(1,2,1,IOM),
     &                       SIGTOT_OPT_VC(2,1,1,IOM)
      END DO
C
      CLOSE (IOTMP)
C
C=======================================================================
99001 FORMAT (f10.5,3(2X,2E16.4))
      END
