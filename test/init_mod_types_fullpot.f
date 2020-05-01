C*==init_mod_types_fullpot.f    processed by SPAG 6.70Rc at 16:58 on 19 Nov 2016
      SUBROUTINE INIT_MOD_TYPES_FULLPOT(FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT,NQCLU,KFP_LMQ,IMQ,IQBOT,IQTOP,NQ
      USE MOD_TYPES,ONLY:IKMCPLWF,IKMCPLWF_RA,IKMCPLWF_LA,IKMCPLWF_RB,
     &    IKMCPLWF_LB,VAMEF,VAMEG,TXT_T,LTXT_T,NTMAX,ISOLIKM,NLMFPT,
     &    NBLK,NSOLBLK,IKMSOLBLK,NFPT,LMIFP,KLMFP,NLMFPMAX,NCPLWFMAX,
     &    NAT,ITBOT,ITTOP,NT
      USE MOD_RMESH,ONLY:ISFLM,KLMSF,LMISF,NSF,NMMAX,NSFMAX,FULLPOT,
     &    NLMSFMAX,SPHERCELL
      USE MOD_ANGMOM,ONLY:NL,NKM,NKMMAX
      USE MOD_FILES,ONLY:IOTMP,IPRINT
      USE MOD_CALCMODE,ONLY:IREL,KMROT,LHS_SOL_EQ_RHS_SOL
      USE MOD_MPI,ONLY:MPI_ID,MPI
      USE MOD_LATTICE,ONLY:SYSTEM_DIMENSION
      USE MOD_SCF,ONLY:SCFSTATUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL FULLPOT5
C
C Local variables
C
      INTEGER I,IA_ERR,IBLK,IBLKTOP,IDIMS,IPOT,IQBOTSAV,IQTOPSAV,IT,
     &        ITBOTSAV,ITTOPSAV,J,J1,J1TOP,MC,NCPLWF
C
C*** End of declarations rewritten by SPAG
C
C-----------------------------------------------------------------------
C       ASA:   set NCPLWFMAX, dummy allocations  and  RETURN
C-----------------------------------------------------------------------
C
      IF ( .NOT.FULLPOT .AND. .NOT.FULLPOT5 ) THEN
C
         NCPLWFMAX = 2
C
         ALLOCATE (VAMEF(1,1,1,1,1,1),VAMEG(1,1,1,1,1,1))
         ALLOCATE (IKMCPLWF(2,NKMMAX))
         ALLOCATE (IKMCPLWF_RA(2,NKMMAX),IKMCPLWF_LA(2,NKMMAX))
         ALLOCATE (IKMCPLWF_RB(2,NKMMAX),IKMCPLWF_LB(2,NKMMAX))
         VAMEF = 999999D0
         VAMEG = 999999D0
         RETURN
      END IF
C
C-----------------------------------------------------------------------
C        run over all sites and types in case of  2D  SCF-start
C-----------------------------------------------------------------------
C
C     IF ( SYSTEM_DIMENSION
C    &    (1:2).EQ.'2D' .AND. SCFSTATUS(1:5).EQ.'START' ) THEN
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &     (SCFSTATUS(1:5).EQ.'START' .OR. SCFSTATUS(1:10)
     &     .EQ.'ITR-L-BULK') ) THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         IQBOTSAV = IQBOT
         IQTOPSAV = IQTOP
         ITBOT = 1
         ITTOP = NT
         IQBOT = 1
         IQTOP = NQ
      END IF
C
C-----------------------------------------------------------------------
C       rotate the shape functions and  non-spherical potentials
C         to the local frames of reference if necessary
C-----------------------------------------------------------------------
C
      IF ( KMROT.NE.0 .AND. IREL.EQ.3 .AND. FULLPOT ) THEN
C
         IF ( .NOT.SPHERCELL ) CALL FP_SFN_LOCAL
C
         CALL FP_POT_LOCAL
C
      END IF
C
C-----------------------------------------------------------------------
C          check or set up the FULL POTENTIAL PICKING RULES
C-----------------------------------------------------------------------
C
      IF ( FULLPOT ) THEN
C
         CALL FPPICKRULES('CHECK',NMMAX,TXT_T,LTXT_T,IMQ,IQAT,NAT,
     &                    KFP_LMQ,KLMFP,NLMFPT,NFPT,LMIFP,ISFLM,KLMSF,
     &                    NSF,LMISF,NLMFPMAX,NLMSFMAX)
C
      ELSE
C
C ------------------ SPHERCELL or start of FULLPOT calculations from ASA
C ------------------------------------ FULLPOT will be set TRUE in <SCF>
C
         CALL FPPICKRULES('SETUP',NMMAX,TXT_T,LTXT_T,IMQ,IQAT,NAT,
     &                    KFP_LMQ,KLMFP,NLMFPT,NFPT,LMIFP,ISFLM,KLMSF,
     &                    NSF,LMISF,NLMFPMAX,NLMSFMAX)
C
      END IF
C
C
      IF ( MPI_ID.EQ.0 ) THEN
C
         CALL FPCOUPL(IPRINT,NKM,NL)
C
         REWIND (IOTMP)
         READ (IOTMP) NCPLWF,IBLKTOP,J1TOP
C
         IF ( NQCLU.EQ.0 ) THEN
            NCPLWFMAX = NCPLWF
         ELSE
            NCPLWFMAX = NKMMAX
         END IF
C
         MC = NCPLWFMAX
         ALLOCATE (VAMEG(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEG'
         ALLOCATE (IKMCPLWF_RA(NCPLWFMAX,NKMMAX))
         ALLOCATE (IKMCPLWF_LA(NCPLWFMAX,NKMMAX))
         ALLOCATE (IKMCPLWF_RB(NCPLWFMAX,NKMMAX))
         ALLOCATE (IKMCPLWF_LB(NCPLWFMAX,NKMMAX))
         ALLOCATE (IKMCPLWF(NCPLWFMAX,NKMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> IKMCPLWF'
         IKMCPLWF = 0
C
         READ (IOTMP) ((((((VAMEG(I,J,IPOT,J1,IBLK,IT),I=1,NCPLWF),J=1,
     &                NCPLWF),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,
     &                IBLKTOP),IT=ITBOT,ITTOP)
         IF ( IREL.LE.2 ) THEN
            ALLOCATE (VAMEF(1,1,1,1,1,1),STAT=IA_ERR)
            VAMEF = 999999D0
         ELSE
            ALLOCATE (VAMEF(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
            IF ( IA_ERR.NE.0 ) STOP 'ALLOC: MAIN -> VAMEF'
            READ (IOTMP) ((((((VAMEF(I,J,IPOT,J1,IBLK,IT),I=1,NCPLWF),J=
     &                   1,NCPLWF),IPOT=1,NLMFPMAX),J1=1,J1TOP),IBLK=1,
     &                   IBLKTOP),IT=ITBOT,ITTOP)
         END IF
         CLOSE (IOTMP)
      END IF
C
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
      IF ( MPI ) THEN
C
         CALL DRV_MPI_BARRIER
         CALL DRV_MPI_BCAST_I(0,NCPLWFMAX,1)
         CALL DRV_MPI_BCAST_I(0,IBLKTOP,1)
         CALL DRV_MPI_BCAST_I(0,J1TOP,1)
         CALL DRV_MPI_BCAST_I(0,NBLK(1),NTMAX)
         CALL DRV_MPI_BCAST_I(0,NSOLBLK(1,1),NKMMAX*NTMAX)
         CALL DRV_MPI_BCAST_I(0,ISOLIKM(1,1),NKMMAX*NTMAX)
         CALL DRV_MPI_BCAST_I(0,IKMSOLBLK(1,1,1),NKMMAX*NKMMAX*NTMAX)
C
         MC = NCPLWFMAX
         IF ( MPI_ID.NE.0 ) THEN
            ALLOCATE (VAMEG(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),STAT=IA_ERR)
            ALLOCATE (IKMCPLWF(NCPLWFMAX,NKMMAX),STAT=IA_ERR)
            IF ( IREL.LE.2 ) THEN
               ALLOCATE (VAMEF(1,1,1,1,1,1),STAT=IA_ERR)
               VAMEF = 999999D0
            ELSE
               ALLOCATE (VAMEF(MC,MC,NLMFPMAX,2,NKMMAX,NTMAX),
     &                   STAT=IA_ERR)
            END IF
         END IF
C
         CALL DRV_MPI_BARRIER
C
         IDIMS = NCPLWFMAX*NCPLWFMAX*NLMFPMAX*2*NKMMAX*NTMAX
         CALL DRV_MPI_BCAST_C(0,VAMEG,IDIMS)
         IF ( IREL.GE.3 ) CALL DRV_MPI_BCAST_C(0,VAMEF,IDIMS)
         IF ( IREL.GE.3 ) CALL DRV_MPI_BCAST_L(0,LHS_SOL_EQ_RHS_SOL,1)
C
      END IF
C  MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPI
C
C     IF ( SYSTEM_DIMENSION
C     &    (1:2).EQ.'2D' .AND. SCFSTATUS(1:5).EQ.'START' ) THEN
C
      IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. 
     &     (SCFSTATUS(1:5).EQ.'START' .OR. SCFSTATUS(1:10)
     &     .EQ.'ITR-L-BULK') ) THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
         IQBOT = IQBOTSAV
         IQTOP = IQTOPSAV
      END IF
C
      END
