C*==dmft_ldau_eref.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_LDAU_EREF(OBS_LT,EREF)
C
C   ********************************************************************
C   *                                                                  *
C   *  Driver to calculate reference energy for DMFT and LDAU          *
C   *  IEREF=0:  EREF= EBAND_LOP/NOS_LOP, e.g. centre of mass for      *
C   *            LOP band (this is calculated from SCF run and does    *
C   *            not correspond to the localised basis set used in the *
C   *            DMFT/LDAU                                             *
C   *            This option does not work well for f-electrons        *
C   *                                                                  *
C   *  IEREF=1:  EREF is set to the resonance energy of the phase      *
C   *            shift. If EREF is not found, setting from             *
C   *            IEREF=0 is taken. Here additional calculations        *
C   *            are needed.                                           *
C   *                                                                  *
C   *  ELSE      EREF taken from input                                 *
C   *                                                                  *
C   *                                                                  *
C   * JM (+Laszlo Oroszlany)                             28.04.2015    *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NLMAX,NOBSMAX,IDOS,IBND
      USE MOD_TYPES,ONLY:BT,VNST,BNST,NTMAX,ITBOT,ITTOP,NT,NLMFPMAX,VT,
     &    CTL,CONC,NAT,LOPT
      USE MOD_DMFT_LDAU,ONLY:IEREF,KSELF,IBASIS,DMFTSIG,DMFT_FIX_BAS
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ENERGY,ONLY:EILOW,EFERMI
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRNSMIN,FLMSF
      USE MOD_CONSTANTS,ONLY:C0,PI
      USE MOD_CALCMODE,ONLY:IREL
      IMPLICIT NONE
C*--DMFT_LDAU_EREF34
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='DMFT_LDAU_EREF')
C
C Dummy arguments
C
      COMPLEX*16 EREF(NTMAX)
      REAL*8 OBS_LT(0:3,NOBSMAX,NLMAX,NTMAX)
C
C Local variables
C
      REAL*8 BNST0(:,:,:),BT0(:,:),COTDELTA,CTL0(:,:),D,E1,E2,EIMAG,
     &       ERES(:,:),PD(:,:),PD0(:,:),PDTAB(:,:,:,:),RAUX,RAUX1,RAUX2,
     &       SFTMP(:,:,:),VNST0(:,:,:),VT0(:,:)
      COMPLEX*16 CWORK(:,:),ETAB(:),MEZJ0(:,:,:,:),MEZZ0(:,:,:,:),
     &           MSST0(:,:,:,:),P,SSST(:,:,:),TSST(:,:,:),WETAB(:)
      LOGICAL FOUND,RESONANCE(:,:)
      INTEGER I,I1,I2,I3,I4,IE,IL,IT,ITBOTSAV,ITTOPSAV,J,K,NE,NLMOP
      CHARACTER*80 STR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MEZJ0,MEZZ0,MSST0,SSST,TSST
      ALLOCATABLE ETAB,WETAB,CWORK,PD,PD0,PDTAB
      ALLOCATABLE SFTMP,ERES,RESONANCE
C
      ALLOCATABLE BT0,VNST0,BNST0
      ALLOCATABLE VT0,CTL0
      ALLOCATE (BT0(NRMAX,NTMAX))
      ALLOCATE (VT0(NRMAX,NTMAX))
      IF ( FULLPOT ) THEN
         ALLOCATE (VNST0(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
         ALLOCATE (BNST0(JRNSMIN:NRMAX,NLMFPMAX,NTMAX))
      END IF
      ALLOCATE (CTL0(NTMAX,NLMAX))
      CALL TRACK_INFO(ROUTINE)
C
      IF ( DMFT_FIX_BAS ) RETURN
      IF ( IEREF.EQ.0 ) THEN
C
C=======================================================
C Calculate new reference energy as a center of mass
C=======================================================
C
         DO IT = ITBOT,ITTOP
            IF ( KSELF(IT).EQ.1 ) THEN
               RAUX = 0.0D0
               RAUX1 = 0.0D0
               RAUX2 = 0.0D0
               RAUX = RAUX + CONC(IT)*NAT(IT)
               IL = LOPT(IT) + 1
               RAUX1 = RAUX1 + OBS_LT(0,IBND,IL,IT)*CONC(IT)*NAT(IT)
               RAUX2 = RAUX2 + OBS_LT(0,IDOS,IL,IT)*CONC(IT)*NAT(IT)
               IF ( ABS(RAUX2).GT.1D-10 ) EREF(IT) = RAUX1/RAUX2
            END IF
         END DO
C
      ELSE IF ( IEREF.LT.0 ) THEN
C
         WRITE (6,*) 'EREF not updated: taken from input'
      ELSE
C
C=======================================================
C Calculate new reference energy from phase shift
C=======================================================
C
         NE = 1000
         EIMAG = 0.0D0
         E1 = 0.01D0
                !EMIN
         E2 = EFERMI + 0.5D0
         ALLOCATE (ETAB(NE),WETAB(NE))
         CALL EPATH(0,E1,E2,EIMAG,NE,ETAB,WETAB,EILOW,IPRINT,NE)
C
         ALLOCATE (MEZJ0(NKMMAX,NKMMAX,NTMAX,NMEMAX))
         ALLOCATE (MEZZ0(NKMMAX,NKMMAX,NTMAX,NMEMAX))
         ALLOCATE (TSST(NKMMAX,NKMMAX,NTMAX))
         ALLOCATE (MSST0(NKMMAX,NKMMAX,NTMAX,NE))
         ALLOCATE (SSST(NKMMAX,NKMMAX,NTMAX))
         ALLOCATE (CWORK(NKMMAX,NKMMAX))
         ALLOCATE (PDTAB(NKMMAX,NKMMAX,NTMAX,NE),PD(NKMMAX,NKMMAX),
     &             PD0(NKMMAX,NKMMAX))
         ALLOCATE (ERES(NKMMAX,NTMAX),RESONANCE(NKMMAX,NTMAX))
         ERES = 99999D0
         RESONANCE = .FALSE.
C
         I = SIZE(FLMSF,1)
         J = SIZE(FLMSF,2)
         K = SIZE(FLMSF,3)
         ALLOCATE (SFTMP(I,J,K))
         SFTMP = FLMSF
         IF ( IBASIS.EQ.1 ) FLMSF(:,1:J,:) = 0.0D0
         IF ( IBASIS.EQ.2 ) FLMSF(:,2:J,:) = 0.0D0
         IF ( IBASIS.EQ.3 ) FLMSF(:,1:J,:) = 0.0D0
         IF ( IBASIS.EQ.4 ) FLMSF(:,2:J,:) = 0.0D0
         CALL RVECCOP(NRMAX*NT,BT,BT0)
         CALL RINIT(NRMAX*NT,BT)
         CALL RVECCOP(NRMAX*NT,VT,VT0)
         CTL0 = CTL
         CTL = CTL/((1D-10)**0.5D0)
         IF ( FULLPOT ) THEN
            VNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
            BNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
C
            VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0.0D0
            BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX) = 0.0D0
         END IF
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         DMFTSIG(1:NKM,1:NKM,1:NT) = C0
         PD = 0.0D0
         PD0 = 0.0D0
         PDTAB = 0.0D0
C
         IF ( IREL.EQ.3 ) THEN
            STR(1:7) = 'REL>CLM'
         ELSE
            STR(1:7) = 'RLM>CLM'
         END IF
C
         DO IE = 1,NE
            DO IT = 1,NT
               IF ( KSELF(IT).EQ.1 ) THEN
                  ITBOT = IT
                  ITTOP = IT
                  CALL RUNSSITE(.FALSE.,0,0,888,.FALSE.,ETAB(IE),P,0,
     &                          TSST,MSST0(1,1,1,IE),SSST,MEZZ0,MEZJ0,
     &                          'NONE      ')
C
                  CALL CHANGEREP(NKM,NKMMAX,MSST0(1,1,IT,IE),STR(1:7),
     &                           CWORK)
                  DO I = 1,NKMMAX
                     COTDELTA = -DREAL(CWORK(I,I)/P)
                     D = ATAN(1D0/COTDELTA) - PD0(I,I)
                     PD(I,I) = D + PD0(I,I)
                     IF ( ABS(D).GT.ABS(D-PI) ) PD(I,I) = D - PI + 
     &                    PD0(I,I)
                     IF ( ABS(D).GT.ABS(D+PI) ) PD(I,I) = D + PI + 
     &                    PD0(I,I)
                     PD0(I,I) = PD(I,I)
                     PDTAB(I,I,IT,IE) = PD(I,I)
                  END DO
C               WRITE(888,99001) DREAL(ETAB(IE)),(PDTAB(I,I,IT,IE),I=1,
C     &              NKMMAX)
               END IF
            END DO
         END DO
C
         VT = VT0
         CTL = CTL0
C
         CALL RVECCOP(NRMAX*NT,BT0,BT)
         IF ( FULLPOT ) THEN
            VNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = VNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
            BNST(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
     &         = BNST0(JRNSMIN:NRMAX,1:NLMFPMAX,1:NTMAX)
         END IF
         FLMSF = SFTMP
C
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
C
         DO IT = ITBOT,ITTOP
            IF ( KSELF(IT).EQ.1 ) THEN
               DO I = 1,NKMMAX
                  DO IE = 2,NE
                     E1 = DREAL(ETAB(IE-1))
                     E2 = DREAL(ETAB(IE))
                     IF ( PDTAB(I,I,IT,IE).GT.PI/2D0 ) THEN
                        ERES(I,IT) = E1 + (E2-E1)
     &                               *(PI/2D0-PDTAB(I,I,IT,IE-1))
     &                               /(PDTAB(I,I,IT,IE-1)
     &                               -PDTAB(I,I,IT,IE))
                        RESONANCE(I,IT) = .TRUE.
                        EXIT
                     END IF
                  END DO
               END DO
               IF ( IPRINT.GT.0 ) THEN
                  WRITE (6,*) 'Resonance, IT,I,ERES'
                  DO I = 1,NKMMAX
                     WRITE (6,*) RESONANCE(I,IT),IT,I,ERES(I,IT)
                  END DO
               END IF
               NLMOP = LOPT(IT)*2 + 1
               I1 = LOPT(IT)**2 + 1
               I2 = I1 + NLMOP - 1
               I3 = NKM/2 + I1
               I4 = I3 + NLMOP - 1
               K = 0
               E1 = 0.0D0
               FOUND = .FALSE.
               DO I = I1,I2
                  IF ( RESONANCE(I,IT) ) THEN
                     K = K + 1
                     E1 = E1 + ERES(I,IT)
                     FOUND = .TRUE.
                  END IF
               END DO
               IF ( IPRINT.GT.0 ) WRITE (6,*) 'EREF SPIN 1:',E1/K
               E2 = 0.0D0
               J = 0
               DO I = I3,I4
                  IF ( RESONANCE(I,IT) ) THEN
                     J = J + 1
                     E2 = E2 + ERES(I,IT)
                     FOUND = .TRUE.
                  END IF
               END DO
               IF ( IPRINT.GT.0 ) WRITE (6,*) 'EREF SPIN 2:',E2/J
C
               IF ( FOUND ) THEN
                  EREF(IT) = (E1+E2)/(K+J)
               ELSE
                  WRITE (6,*) 
     &                       'Warning: EREF not found from phase shifts'
                  WRITE (6,*) ' IT=',IT,' LOP=',LOPT(IT)
                  WRITE (6,*) ' EREF taken from centre of mass'
                  RAUX = 0.0D0
                  RAUX1 = 0.0D0
                  RAUX2 = 0.0D0
                  RAUX = RAUX + CONC(IT)*NAT(IT)
                  IL = LOPT(IT) + 1
                  RAUX1 = RAUX1 + OBS_LT(0,IBND,IL,IT)*CONC(IT)*NAT(IT)
                  RAUX2 = RAUX2 + OBS_LT(0,IDOS,IL,IT)*CONC(IT)*NAT(IT)
                  IF ( ABS(RAUX2).GT.1D-10 ) EREF(IT) = RAUX1/RAUX2
               END IF
            END IF
         END DO
      END IF
C
      DO IT = ITBOT,ITTOP
         IF ( KSELF(IT).EQ.1 ) WRITE (6,99001) IT,DREAL(EREF(IT))
      END DO
      WRITE (6,'(/)')
C	  DMFT_FIX_BAS=.TRUE.
99001 FORMAT (10X,'DMFT and LDAU basis functions:',/,
     &        'reference energy for IT',I3,' = ',F20.8)
      END
