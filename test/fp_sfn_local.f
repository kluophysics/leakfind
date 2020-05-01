C*==fp_sfn_local.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FP_SFN_LOCAL
C   ********************************************************************
C   *                                                                  *
C   *  rotate the SHAPE FUNCTIONs to a local frame                     *
C   *  for which the exchange field B is spin diagonal                 *
C   *                                                                  *
C   *  NOTE: the shape functions a stored in a compressed way          *
C   *        i.e. the LM-index in mapped to ISF=1,..,NSF  via  ISFLM   *
C   *                                                                  *
C   *        NLSF = 4 * NL - 3                                         *
C   *                                                                  *
C   ********************************************************************
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_SITES,ONLY:IQBOT,IQTOP,QMPHI,QMTET,QMGAM,IMQ
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_RMESH,ONLY:NSF,NM,ISFLM,LMISF,KLMSF,NSFMAX,NMMAX,NRSFMAX,
     &    FLMSF,NRSFTOT,NLMSFMAX,NLSF,NLMSF,WINTLM
      USE MOD_SYMMETRY,ONLY:IQREPQ
      IMPLICIT NONE
C*--FP_SFN_LOCAL21
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 ALF,BET,DROT_SF(:,:),D_21,FLMSF_NEW(:),FLMSF_OLD(:,:,:),
     &       GAM,RSUM
      LOGICAL DONE_M(:),FLMSF_NE_0
      INTEGER IA_ERR,IM,IQ,IQM(:),IR,ISF,ISF_OLD,L,LM,LM1,LM2,
     &        LMISF_OLD(:,:),L_LAST,NRSF,NSFMAX_NEW,NSFMAX_OLD,NSF_NEW,
     &        NSF_OLD(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE FLMSF_NEW,NSF_OLD,FLMSF_OLD,LMISF_OLD,DROT_SF
      ALLOCATABLE DONE_M,IQM
C
C=======================================================================
C                             INITIALIZE
C=======================================================================
C
      IF ( NLMSFMAX.NE.NLMSF ) STOP '<FP_SFN_LOCAL>: NLMSFMAX <> NLMSF'
      IF ( NMMAX.NE.NM ) STOP '<FP_SFN_LOCAL> NMMAX <>  NM'
C
      ALLOCATE (DONE_M(NM),IQM(NM))
      DONE_M(1:NM) = .FALSE.
C
C ----------------------------------------------------------------------
C
      ALLOCATE (DROT_SF(NLMSF,NLMSF),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: FP_SFN_LOCAL -> DROT_SF'
C
C=======================================================================
C  all sites asociated with a mesh IM have to have same (TET,PHI,GAM)
C=======================================================================
C
      DO IM = 1,NM
         ALF = 999999D0
         DO IQ = IQBOT,IQTOP
            IF ( IM.EQ.IMQ(IQ) ) THEN
               ALF = QMTET(IQ)
               BET = QMPHI(IQ)
               GAM = QMGAM(IQ)
               EXIT
            END IF
         END DO
C
         DO IQ = IQBOT,IQTOP
            IF ( IM.EQ.IMQ(IQ) ) THEN
               IF ( ABS(ALF-QMTET(IQ)).GT.1D-6 .OR. ABS(BET-QMPHI(IQ))
     &              .GT.1D-6 .OR. ABS(GAM-QMGAM(IQ)).GT.1D-6 ) THEN
                  WRITE (6,99001)
                  STOP
               END IF
            END IF
         END DO
      END DO
C
C-----------------------------------------------------------------------
C      CHECK:  find maximum number of shape functions that may occur
C-----------------------------------------------------------------------
C
      NSFMAX_NEW = 0
      NSFMAX_OLD = 0
      DO IM = 1,NM
         NSFMAX_OLD = MAX(NSFMAX_OLD,NSF(IM))
         NSF_NEW = 0
         L_LAST = -10
         DO ISF = 1,NSF(IM)
            LM = LMISF(ISF,IM)
            L = L_LM(LM)
            IF ( L.NE.L_LAST ) NSF_NEW = NSF_NEW + (L+1)*(L+1)
            L_LAST = L
         END DO
         NSFMAX_NEW = MAX(NSFMAX_NEW,NSF_NEW)
      END DO
C
      IF ( NSFMAX_NEW.LT.NSFMAX_OLD )
     &      STOP 'in <FP_SFN_LOCAL> NSFMAX_NEW < NSFMAX_OLD'
      NSFMAX_OLD = NSFMAX
C
C-----------------------------------------------------------------------
C    swap the shape function information on temporary arrays  *_OLD
C-----------------------------------------------------------------------
C
      ALLOCATE (FLMSF_NEW(NRSFMAX))
C
      ALLOCATE (NSF_OLD(NM),FLMSF_OLD(NRSFMAX,NSFMAX,NM))
      ALLOCATE (LMISF_OLD(NLMSFMAX,NM))
C
      NSF_OLD(1:NM) = NSF(1:NM)
      FLMSF_OLD(1:NRSFMAX,1:NSFMAX,1:NM)
     &   = FLMSF(1:NRSFMAX,1:NSFMAX,1:NM)
      LMISF_OLD(1:NLMSFMAX,1:NM) = LMISF(1:NLMSFMAX,1:NM)
C
C=======================================================================
C                          find new  NSFMAX
C=======================================================================
C
      NSFMAX_NEW = 0
C
      DO IQ = IQBOT,IQTOP
         IF ( IQ.EQ.IQREPQ(IQ) ) THEN
            IM = IMQ(IQ)
            IQM(IM) = IQ
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for REAL spherical harmonics
C ----------------------------------------------------------------------
C
            CALL ROTMAT_RYLM(NLSF,NLMSF,QMPHI(IQ),QMTET(IQ),QMGAM(IQ),
     &                       DROT_SF)
C
C-----------------------------------------------------------------------
C
            NSF(IM) = 0
            DO LM1 = 1,NLMSF
C
               FLMSF_NE_0 = .FALSE.
C
               DO ISF_OLD = 1,NSF_OLD(IM)
                  LM2 = LMISF_OLD(ISF_OLD,IM)
                  D_21 = DROT_SF(LM2,LM1)
C
C
                  IF ( ABS(D_21).GT.1D-8 ) FLMSF_NE_0 = .TRUE.
               END DO
C
               IF ( FLMSF_NE_0 ) NSF(IM) = NSF(IM) + 1
C
            END DO
C
            NSFMAX_NEW = MAX(NSFMAX_NEW,NSF(IM))
C
C-----------------------------------------------------------------------
C
         END IF
      END DO
C
C=======================================================================
C       reallocate arrays
C=======================================================================
C
      NSFMAX = NSFMAX_NEW
C
      DEALLOCATE (WINTLM,LMISF,FLMSF)
      ALLOCATE (WINTLM(NRSFMAX,NSFMAX,NMMAX))
      ALLOCATE (LMISF(NLMSFMAX,NM),FLMSF(NRSFMAX,NSFMAX,1:NM))
C
      NSF(1:NM) = 0
      FLMSF(1:NRSFMAX,1:NSFMAX,1:NM) = 0D0
      ISFLM(1:NLMSFMAX,1:NM) = 0
      KLMSF(1:NLMSFMAX,1:NM) = 0
      LMISF(1:NLMSFMAX,1:NM) = 0
C
C=======================================================================
C                     rotate the shape functions
C=======================================================================
C
      NSFMAX_NEW = 0
C
      DO IQ = IQBOT,IQTOP
         IF ( IQ.EQ.IQREPQ(IQ) ) THEN
            IM = IMQ(IQ)
            IF ( .NOT.DONE_M(IM) ) THEN
C-----------------------------------------------------------------------
C
               DONE_M(IM) = .TRUE.
C
C ----------------------------------------------------------------------
C         rotation matrix  DROT for REAL spherical harmonics
C ----------------------------------------------------------------------
C
               CALL ROTMAT_RYLM(NLSF,NLMSF,QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
     &                          ,DROT_SF)
C
C-----------------------------------------------------------------------
C
               NSF(IM) = 0
               DO LM1 = 1,NLMSF
C
                  FLMSF_NEW(1:NRSFMAX) = 0D0
                  FLMSF_NE_0 = .FALSE.
                  NRSF = NRSFTOT(IM)
C
                  DO ISF_OLD = 1,NSF_OLD(IM)
                     LM2 = LMISF_OLD(ISF_OLD,IM)
                     D_21 = DROT_SF(LM2,LM1)
                     IF ( ABS(D_21).GT.1D-8 ) THEN
C
C
                        FLMSF_NEW(1:NRSF) = FLMSF_NEW(1:NRSF)
     &                     + FLMSF_OLD(1:NRSF,ISF_OLD,IM)*D_21
C
                        FLMSF_NE_0 = .TRUE.
                     END IF
                  END DO
C
C--------------------------------------------- check for vanishing terms
                  RSUM = 0D0
                  DO IR = 1,NRSF
                     RSUM = RSUM + ABS(FLMSF_NEW(IR))
                  END DO
                  IF ( RSUM.LT.1D-5 ) FLMSF_NE_0 = .FALSE.
C
                  IF ( FLMSF_NE_0 ) THEN
                     NSF(IM) = NSF(IM) + 1
                     ISF = NSF(IM)
                     LMISF(ISF,IM) = LM1
                     KLMSF(LM1,IM) = 1
                     ISFLM(LM1,IM) = ISF
                     FLMSF(1:NRSF,ISF,IM) = FLMSF_NEW(1:NRSF)
                  END IF
C
               END DO
C
               NSFMAX_NEW = MAX(NSFMAX_NEW,NSF(IM))
C
C-----------------------------------------------------------------------
C
            END IF
         END IF
      END DO
C
C=======================================================================
C
      WRITE (6,99002) NSFMAX,NSFMAX_OLD
      DO IM = 1,NM
C
         IQ = IQM(IM)
C
         WRITE (6,99003) IM,IQ,QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
C
         IF ( IPRINT.GE.-1 ) THEN
            WRITE (6,99006) 'old settings '
            WRITE (6,99004) 'NSF ',NSF_OLD(IM),'NRSFTOT ',NRSFTOT(IM)
            WRITE (6,99005) (LMISF_OLD(ISF,IM),L_LM(LMISF_OLD(ISF,IM)),
     &                      M_LM(LMISF_OLD(ISF,IM)),ISF=1,NSF_OLD(IM))
            WRITE (6,99006) 'new settings '
         END IF
C
         WRITE (6,99004) 'NSF ',NSF(IM),'NRSFTOT ',NRSFTOT(IM)
         WRITE (6,99005) (LMISF(ISF,IM),L_LM(LMISF(ISF,IM)),
     &                   M_LM(LMISF(ISF,IM)),ISF=1,NSF(IM))
      END DO
C
      WRITE (6,99003)
C
C-----------------------------------------------------------------------
C
99001 FORMAT (/,10X,'<FP_SFN_LOCAL>:  different rotation angles',
     &        ' found for equivalent sites ',/)
99002 FORMAT (/,1X,79('*'),/,33X,'<FP_SFN_LOCAL>',/,1X,79('*'),//,10X,
     &        'rotated shape functions generated',//,10X,'NSFMAX =',I6,
     &        '   new number of SF''s ',/,I24,
     &        '   setting prior to rotation',/)
99003 FORMAT (/,:,10X,'mesh type IM =',I4,/,:,10X,'for site  IQ =',I4,
     &        5X,'PHI =',F7.2,2X,'TET =',F7.2,2X,'GAM =',F7.2,/)
99004 FORMAT (10X,A,I4,4X,A,I5)
99005 FORMAT (10X,'LM  ',5(I4,' (',I2,',',I3,')'),:,/,
     &        (14X,5(I4,' (',I2,',',I3,')')))
99006 FORMAT (/,10X,A,/)
      END
