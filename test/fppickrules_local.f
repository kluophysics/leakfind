C*==fppickrules_local.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPPICKRULES_LOCAL(NFP_Q,KFP_LMQ,NFP_Q_LOC,KFP_LMQ_LOC,
     &                             NLMFPMAX,NSF_Q,KLM_SFQ,NLMSFMAX,
     &                             NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *  transfer the FP picking rules determined                        *
C   *  for the GLOBAL FRAME  to  the  LOCAL FRAME                      *
C   *                                                                  *
C   *  NOTE: the array sizes are chosen to deal with worst case        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_RMESH,ONLY:NLSF,NLMSF
      USE MOD_SITES,ONLY:IQBOT,IQTOP,QMPHI,QMTET,QMGAM
      USE MOD_TYPES,ONLY:NLFP,NLMFP
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      IMPLICIT NONE
C*--FPPICKRULES_LOCAL20
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NLMFPMAX,NLMSFMAX,NQMAX
      INTEGER KFP_LMQ(NLMFPMAX,NQMAX),KFP_LMQ_LOC(NLMFPMAX,NQMAX),
     &        KLM_SFQ(NLMSFMAX,NQMAX),NFP_Q(NQMAX),NFP_Q_LOC(NQMAX),
     &        NSF_Q(NQMAX)
C
C Local variables
C
      REAL*8 DROT_RYLM(:,:),D_21
      INTEGER IA_ERR,IFP,IQ,ISF,KLM_SFQ_GLO(:,:),LM,LM1,LM2,LM_FPQ(:,:),
     &        LM_FPQ_LOC(:,:),LM_SFQ(:,:),LM_SFQ_GLO(:,:),N,NFPMAX,
     &        NFPMAX_LOC,NSFMAX_GLO,NSFMAX_LOC,NSF_Q_GLO(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DROT_RYLM
      ALLOCATABLE LM_FPQ,LM_FPQ_LOC
      ALLOCATABLE NSF_Q_GLO,KLM_SFQ_GLO,LM_SFQ_GLO,LM_SFQ
C
      WRITE (6,99001)
C
C=======================================================================
C                     rotate the non-spherical potential
C=======================================================================
C
      ALLOCATE (LM_FPQ(NLMFPMAX,NQMAX),LM_FPQ_LOC(NLMFPMAX,NQMAX))
C
      NFP_Q_LOC(IQBOT:IQTOP) = 0
      KFP_LMQ_LOC(1:NLMFPMAX,IQBOT:IQTOP) = 0
C
C ----------------------------------------------------------------------
C
      ALLOCATE (DROT_RYLM(NLMFP,NLMFP),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: FPPICKRULES_LOCAL -> DROT_RYLM'
C
C-----------------------------------------------------------------------
C
      NFPMAX = 0
      DO IQ = IQBOT,IQTOP
         N = 0
         DO LM = 1,NLMFPMAX
            IF ( KFP_LMQ(LM,IQ).EQ.1 ) THEN
               N = N + 1
               LM_FPQ(N,IQ) = LM
            END IF
         END DO
         NFPMAX = MAX(NFPMAX,N)
         IF ( N.NE.NFP_Q(IQ) ) STOP '<FPPICKRULES_LOCAL>   >>>>   NFP_Q'
      END DO
C
      NFPMAX_LOC = 0
      DO IQ = IQBOT,IQTOP
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for REAL spherical harmonics
C ----------------------------------------------------------------------
C
         CALL ROTMAT_RYLM(NLFP,NLMFP,QMPHI(IQ),QMTET(IQ),QMGAM(IQ),
     &                    DROT_RYLM)
C
         DO LM1 = 1,NLMFP
C
            DO LM2 = 1,NLMFP
C
               D_21 = DROT_RYLM(LM2,LM1)
               IF ( ABS(D_21).GT.1D-8 .AND. KFP_LMQ(LM2,IQ).EQ.1 ) THEN
                  NFP_Q_LOC(IQ) = NFP_Q_LOC(IQ) + 1
                  LM_FPQ_LOC(NFP_Q_LOC(IQ),IQ) = LM1
                  KFP_LMQ_LOC(LM1,IQ) = 1
                  EXIT
               END IF
C
            END DO
C
         END DO
C
         NFPMAX_LOC = MAX(NFPMAX_LOC,NFP_Q_LOC(IQ))
C
      END DO
C-----------------------------------------------------------------------
C
      WRITE (6,99002) 'NON-SPHERICAL POTENTIAL ','NFPMAX',NFPMAX,
     &                NFPMAX_LOC
      DO IQ = IQBOT,IQTOP
         WRITE (6,99003) 'site',IQ,QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
         IF ( IPRINT.GE.-1 ) THEN
            WRITE (6,99007) 'global settings '
            WRITE (6,99005) 'NFP ',NFP_Q(IQ)
            WRITE (6,99006) (LM_FPQ(IFP,IQ),L_LM(LM_FPQ(IFP,IQ)),
     &                      M_LM(LM_FPQ(IFP,IQ)),IFP=1,NFP_Q(IQ))
            WRITE (6,99007) 'local settings '
         END IF
         WRITE (6,99005) 'NFP ',NFP_Q_LOC(IQ)
         WRITE (6,99006) (LM_FPQ_LOC(IFP,IQ),L_LM(LM_FPQ_LOC(IFP,IQ)),
     &                   M_LM(LM_FPQ_LOC(IFP,IQ)),IFP=1,NFP_Q_LOC(IQ))
      END DO
C
      WRITE (6,99004)
C
C=======================================================================
C                     rotate the shape functions
C=======================================================================
C
      ALLOCATE (NSF_Q_GLO(NQMAX),KLM_SFQ_GLO(NLMSFMAX,NQMAX))
      ALLOCATE (LM_SFQ_GLO(NLMSFMAX,NQMAX),LM_SFQ(NLMSFMAX,NQMAX))
C
      NSF_Q_GLO(IQBOT:IQTOP) = NSF_Q(IQBOT:IQTOP)
      KLM_SFQ_GLO(1:NLMSFMAX,IQBOT:IQTOP)
     &   = KLM_SFQ(1:NLMSFMAX,IQBOT:IQTOP)
      NSF_Q(IQBOT:IQTOP) = 0
      KLM_SFQ(1:NLMSFMAX,IQBOT:IQTOP) = 0
C
C ----------------------------------------------------------------------
C
      DEALLOCATE (DROT_RYLM)
      ALLOCATE (DROT_RYLM(NLMSF,NLMSF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: FPPICKRULES_LOCAL -> DROT_RYLM'
C
C-----------------------------------------------------------------------
C
      NSFMAX_GLO = 0
      DO IQ = IQBOT,IQTOP
         N = 0
         DO LM = 1,NLMSFMAX
            IF ( KLM_SFQ_GLO(LM,IQ).EQ.1 ) THEN
               N = N + 1
               LM_SFQ_GLO(N,IQ) = LM
            END IF
         END DO
         NSFMAX_GLO = MAX(NSFMAX_GLO,N)
         IF ( N.NE.NSF_Q_GLO(IQ) ) STOP 
     &                          '<FPPICKRULES_LOCAL>   >>>>   NSF_Q_GLO'
      END DO
C
      NSFMAX_LOC = 0
      DO IQ = IQBOT,IQTOP
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for REAL spherical harmonics
C ----------------------------------------------------------------------
C
         CALL ROTMAT_RYLM(NLSF,NLMSF,QMPHI(IQ),QMTET(IQ),QMGAM(IQ),
     &                    DROT_RYLM)
C
         DO LM1 = 1,NLMSF
C
            DO LM2 = 1,NLMSF
C
               D_21 = DROT_RYLM(LM2,LM1)
               IF ( ABS(D_21).GT.1D-8 .AND. KLM_SFQ_GLO(LM2,IQ).EQ.1 )
     &              THEN
C
                  NSF_Q(IQ) = NSF_Q(IQ) + 1
                  LM_SFQ(NSF_Q(IQ),IQ) = LM1
                  KLM_SFQ(LM1,IQ) = 1
                  EXIT
C
               END IF
C
            END DO
C
         END DO
C
         NSFMAX_LOC = MAX(NSFMAX_LOC,NSF_Q(IQ))
C
      END DO
C-----------------------------------------------------------------------
C
      WRITE (6,99002) 'SHAPE FUNCTIONS ','NSFMAX',NSFMAX_GLO,NSFMAX_LOC
      DO IQ = IQBOT,IQTOP
         WRITE (6,99003) 'site',IQ,QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
         IF ( IPRINT.GE.-1 ) THEN
            WRITE (6,99007) 'global settings '
            WRITE (6,99005) 'NSF ',NSF_Q_GLO(IQ)
            WRITE (6,99006) (LM_SFQ_GLO(ISF,IQ),L_LM(LM_SFQ_GLO(ISF,IQ))
     &                      ,M_LM(LM_SFQ_GLO(ISF,IQ)),ISF=1,
     &                      NSF_Q_GLO(IQ))
            WRITE (6,99007) 'local settings '
         END IF
         WRITE (6,99005) 'NSF ',NSF_Q(IQ)
         WRITE (6,99006) (LM_SFQ(ISF,IQ),L_LM(LM_SFQ(ISF,IQ)),
     &                   M_LM(LM_SFQ(ISF,IQ)),ISF=1,NSF_Q(IQ))
      END DO
C
      WRITE (6,99004)
      WRITE (6,99004) 'THESE ARE ONLY TEMPORARY SITE DEPENDENT SETTINGS'
C
C=======================================================================
99001 FORMAT (/,1X,79('*'),/,30X,'<FPPICKRULES_LOCAL>',/,1X,79('*'),//,
     &        10X,'transfering the site dependent picking rules ',/,10X,
     &        'from the global to the local frames of references ',/,
     &        10X,'according to the orientation of the magnetic moment',
     &        /)
99002 FORMAT (10X,A,//,10X,A,' =',I6,
     &        '  original setting for the global frame of reference',/,
     &        I24,
     &        '  new      setting for the local frames of references')
99003 FORMAT (/,:,10X,A,I4,5X,'PHI =',F7.2,5X,'TET =',F7.2,5X,'GAM =',
     &        F7.2,/)
99004 FORMAT (/,:,10X,A,I4,/)
99005 FORMAT (10X,A,I4,4X,A,I5)
99006 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99007 FORMAT (/,10X,A,/)
      END
