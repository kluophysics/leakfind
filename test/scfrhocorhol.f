C*==scfrhocorhol.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFRHOCORHOL(COREHOLE,ITRSCF,IKMCOR,NKPCOR,KAPCOR,
     &                        MM05COR,IZERO,ECOR,SZCOR,GCOR,FCOR,BCOR,
     &                        BCORS,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *   subroutine to support SCF-calculations with a CORE HOLE        *
C   *   for a slected atom type   ITHOLE                               *
C   *                                                                  *
C   *   First call: initialize calculation                             *
C   *               i.e. fix the atom type ITHOLE to deal with         *
C   *               and its involved core state ICSTHOLE               *
C   *               with NCXRAY  and  LCXRAY                           *
C   *                                                                  *
C   *   Subsequent calls:
C   *            calculate charge density for state  ITHOLE,ICSTHOLE   *
C   *            and subtract it from the core charge density RHOCHRC  *
C   *                                                                  *
C   *  OS updated May 2015                                             *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_SITES,ONLY:IQAT,ICPA
      USE MOD_RMESH,ONLY:NRMAX,FULLPOT,JRWS,R2DRDI
      USE MOD_FILES,ONLY:IPRINT,FOUND_SECTION,FOUND_INTEGER
      USE MOD_TYPES,ONLY:NVALTOT,NVALT,NTMAX,LCXRAY,NCXRAY,IMT,RHOSPNC,
     &    RHOCHRC
      USE MOD_MPI,ONLY:MPI_ID
      USE MOD_TYPES,ONLY:ITBOT,ITTOP
      IMPLICIT NONE
C*--SCFRHOCORHOL31
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL COREHOLE
      INTEGER ITRSCF,NCSTMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),SZCOR(NCSTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IZERO(NCSTMAX),KAPCOR(NCSTMAX),
     &        MM05COR(NCSTMAX),NKPCOR(NCSTMAX)
C
C Local variables
C
      INTEGER ICST,ICSTHOLE,IM,IR,IRTOP,IT,ITHOLE,K,N,NCSTHOLE
      LOGICAL INITIALIZE
      REAL*8 RAUX,RHOCHRC_SAV(:,:),RHOSPNC_SAV(:,:),RINT(:),SCLHOLE
      SAVE SCLHOLE
C
C*** End of declarations rewritten by SPAG
C
C
      DATA ICSTHOLE/0/,ITHOLE/0/,INITIALIZE/.TRUE./
C
      ALLOCATABLE RINT,RHOCHRC_SAV,RHOSPNC_SAV
C
C ======================================================================
C                 initialize treatment of core holes
C ======================================================================
C
      IF ( INITIALIZE ) THEN
C
         DO IT = 1,NTMAX
            NCXRAY(IT) = 0
            LCXRAY(IT) = 0
         END DO
C
         ITHOLE = 1
         NCXRAY(ITHOLE) = 2
         LCXRAY(ITHOLE) = 1
C
         CALL INPUT_FIND_SECTION('SCF',0)
C
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_SET_INTEGER('ITHOLE',ITHOLE,9999,0)
            COREHOLE = FOUND_INTEGER
C
            IF ( COREHOLE ) THEN
C
               CALL SECTION_SET_REAL('SCALEHOLE',SCLHOLE,1.0D0,0)
               CALL SECTION_SET_INTEGER('ICSTHOLE',ICSTHOLE,9999,0)
               CALL SECTION_SET_INTEGER('NQNHOLE',NCXRAY(ITHOLE),9999,0)
               CALL SECTION_SET_INTEGER('LQNHOLE',LCXRAY(ITHOLE),9999,0)
C
               WRITE (6,99003)
               WRITE (6,99002) 'ITHOLE    :  ',ITHOLE
               WRITE (6,99002) 'NQNHOLE   :  ',NCXRAY(ITHOLE)
               WRITE (6,99002) 'LQNHOLE   :  ',LCXRAY(ITHOLE)
               WRITE (6,99002) 'ICSTHOLE  :  ',ICSTHOLE
               WRITE (6,99005) 'SCALEHOLE :  ',SCLHOLE
            END IF
         END IF
C
         IF ( COREHOLE .AND. FULLPOT ) STOP 
     &       'in <SCFRHOCORHOL>  FULLPOT-mode not allowed for core hole'
C
         INITIALIZE = .FALSE.
         RETURN
      END IF
C=======================================================================
C
C
C ======================================================================
C                       account for CORE HOLE
C                      only for master process
C ======================================================================
C
      IF ( MPI_ID.EQ.0 .AND. COREHOLE ) THEN
C
         ALLOCATE (RINT(NRMAX))
         ALLOCATE (RHOCHRC_SAV(NRMAX,NTMAX),RHOSPNC_SAV(NRMAX,NTMAX))
C
C NOTE: for ITHOLE <> 0    <CORE> deals only with the core states
C       with quantum numbers: NCXRAY(ITHOLE) NCXRAY(ITHOLE)
C       and writes charge to  RHOCHR  and  RHOSPN
C
         IT = ITHOLE
         IM = IMT(IT)
         IRTOP = JRWS(IMT(IT))
         RHOCHRC_SAV(1:IRTOP,ITBOT:ITTOP) = RHOCHRC(1:IRTOP,ITBOT:ITTOP)
         RHOSPNC_SAV(1:IRTOP,ITBOT:ITTOP) = RHOSPNC(1:IRTOP,ITBOT:ITTOP)
C
         CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &             IKMCOR,IZERO,ITHOLE,BCOR,BCORS,NCSTMAX)
C
         RHOCHRC(1:IRTOP,ITBOT:ITTOP) = RHOCHRC_SAV(1:IRTOP,ITBOT:ITTOP)
         RHOSPNC(1:IRTOP,ITBOT:ITTOP) = RHOSPNC_SAV(1:IRTOP,ITBOT:ITTOP)
C
         NCSTHOLE = 4*LCXRAY(ITHOLE) + 2
         IF ( ICSTHOLE.NE.0 ) NCSTHOLE = 1
C
         RINT(1:NRMAX) = 0.0D0
C
         DO ICST = 1,NCSTHOLE
            IF ( ICST.EQ.ICSTHOLE .OR. ICSTHOLE.EQ.0 ) THEN
               DO K = 1,NKPCOR(ICST)
                  DO IR = 1,JRWS(IM)
                     RINT(IR) = RINT(IR)
     &                          + (GCOR(IR,K,ICST)**2+FCOR(IR,K,ICST)
     &                          **2)
                  END DO
               END DO
            END IF
         END DO
C
         DO N = 1,JRWS(IM)
            RHOCHRC(N,IT) = RHOCHRC(N,IT) - SCLHOLE*RINT(N)
     &                      /DBLE(NCSTHOLE)
            RINT(N) = R2DRDI(N,IM)*RHOCHRC(N,IT)
         END DO
C
         CALL RRADINT(IM,RINT,RAUX)
C
         WRITE (6,99001) IT,RAUX
C
C--------------------------- do not change NVALTOT in case of a CPA site
         IF ( ICPA(IQAT(1,IT)).EQ.1 .AND. ITRSCF.EQ.1 ) THEN
C
            NVALTOT = NVALTOT + 1.0*SCLHOLE
            NVALT(IT) = NVALT(IT) + 1
C
            WRITE (6,99004) IT,NVALTOT,NVALT(IT)
C
         END IF
C
      END IF
C=======================================================================
C
      RETURN
C=======================================================================
99001 FORMAT (/,10X,'CORE HOLE formation',/,
     &        ' integrated charge density for atom type ',I2,':',F12.8)
99002 FORMAT (10X,A,I10)
99003 FORMAT (//,1X,79('*'),/,28X,'inclusion of a core hole',/,1X,
     &        79('*'))
99004 FORMAT (/,10X,'CORE HOLE formation',/,
     &        ' number of valence electrons changed by  1: ',I2,':',
     &        F12.8,I4)
99005 FORMAT (10X,A,5F20.10)
      END
