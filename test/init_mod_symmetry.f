C*==init_mod_symmetry.f    processed by SPAG 6.70Rc at 08:14 on  8 Mar 2017
      SUBROUTINE INIT_MOD_SYMMETRY(TASK,INITELOOP,IREL,NONMAG,MOL,KMROT,
     &                             QMVEC,IPRINT,NT,NT0,NA_TAUX,IQ_ATAUX,
     &                             IT_OQAUX,IT0_TAUX,NL,BRAVAIS,NKMMAX,
     &                             NTMAX,NTLIM)
C   ********************************************************************
C   *                                                                  *
C   *  module to store all tables dependent ONLY on   SYMMETRY         *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_SYMMETRY,ONLY:SYMEULANG,SYMTVEC,DROT,IQPSYMQ,IQREPMSYM,
     &    IQREPQ,IQORGQP,ISYMGENQ,NSYM,NSYMACCEPTED,NSYMCRYSYS,NSFTSYMQ,
     &    SYMDET,SYMACCEPTED,SYMCRYSYS,SYMUNITARY,SYMSYMBL,IQEQ,NQEQ,
     &    IWEDGEROT,MROTK,MROTR,NWEDGE,NMSYM,NCL,NMB_CL,IQ_MBCL,NSYMMAX,
     &    MREP_Q
      USE MOD_ANGMOM,ONLY:NKM
      USE MOD_SITES,ONLY:NQ,NQMAX,QMPHI,QMTET,NOQ,QBAS
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='INIT_MOD_SYMMETRY')
C
C Dummy arguments
C
      INTEGER BRAVAIS,IPRINT,IREL,KMROT,NKMMAX,NL,NT,NT0,NTLIM,NTMAX
      LOGICAL INITELOOP,MOL,NONMAG
      CHARACTER*10 TASK
      INTEGER IQ_ATAUX(NQMAX,NTLIM),IT0_TAUX(NTLIM),
     &        IT_OQAUX(NTLIM,NQMAX),NA_TAUX(NTLIM)
      REAL*8 QMVEC(3)
C
C Local variables
C
      INTEGER I,IA_ERR,ICL,IMB,IQ,ISYM,J,NK,NMBMAX
      COMPLEX*16 W1(:,:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W1
C
      ALLOCATE (IQREPMSYM(NQMAX),IQREPQ(NQMAX))
      ALLOCATE (ISYMGENQ(NQMAX),IQPSYMQ(NSYMMAX,NQMAX))
      ALLOCATE (NSFTSYMQ(3,NSYMMAX,NQMAX),IQORGQP(NSYMMAX,NQMAX))
      ALLOCATE (IQEQ(NQMAX,NQMAX),NQEQ(NQMAX))
      ALLOCATE (DROT(NKMMAX,NKMMAX,NSYMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: INIT_MOD_SYMMETRY -> CREL'
C
      NSFTSYMQ = 999999
      IQREPQ = 999999
      IQPSYMQ = 999999
      IQREPMSYM = 999999
      IQEQ = 999999
      IQORGQP = 999999
C
      ALLOCATE (W1(NKMMAX,NKMMAX))
C
      CALL SYMINIT(IPRINT,BRAVAIS,NSYM,NSYMCRYSYS,SYMCRYSYS,SYMSYMBL,
     &             SYMEULANG)
C
C-----------------------------------------------------------------------
C
      CALL OPEN_IOTMP_SCRATCH(ROUTINE,IOTMP)
C
      CALL SYMLATTICE(IOTMP,MOL,IPRINT,NWEDGE,IWEDGEROT,NQ,QBAS,ABAS,
     &                NONMAG,IREL,KMROT,QMVEC,MROTR,MROTK,NSYM,
     &                NSYMCRYSYS,SYMACCEPTED,SYMUNITARY,SYMCRYSYS,
     &                SYMTVEC,SYMSYMBL,SYMEULANG,NOQ,IT_OQAUX,IT0_TAUX,
     &                NT0,NT,NA_TAUX,IQ_ATAUX,QMPHI,QMTET,IQORGQP,
     &                IQPSYMQ,NSFTSYMQ,NSYMACCEPTED,SYMDET,ISYMGENQ,
     &                IQREPQ,NQEQ,IQEQ,NMSYM,IQREPMSYM,NTMAX,NQMAX)
C
C--------------------- keep information on inequivalent sites or classes
C
      REWIND IOTMP
      READ (IOTMP) NCL,NMBMAX
      ALLOCATE (NMB_CL(NCL),IQ_MBCL(NMBMAX,NCL))
C
      DO ICL = 1,NCL
         READ (IOTMP) NMB_CL(ICL)
         READ (IOTMP) (IQ_MBCL(IMB,ICL),IMB=1,NMB_CL(ICL))
      END DO
C
      CLOSE (IOTMP)
C
C-----------------------------------------------------------------------
C   provide real space rotation matrices connecting equivanent sites
C-----------------------------------------------------------------------
C
      ALLOCATE (MREP_Q(3,3,NQMAX))
      MREP_Q(:,:,:) = 0D0
C
      DO IQ = 1,NQ
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            DO I = 1,3
               MREP_Q(I,I,IQ) = 1D0
            END DO
C
         ELSE
C
            ISYM = ISYMGENQ(IQ)
C
            CALL GETMROT(-SYMEULANG(3,ISYM),-SYMEULANG(2,ISYM),
     &                   -SYMEULANG(1,ISYM),MREP_Q(1,1,IQ))
C
         END IF
C
      END DO
C
C-----------------------------------------------------------------------
C
      IF ( INITELOOP .OR. TASK(1:5).EQ.'SIGMA' ) THEN
C
         IF ( IREL.LT.2 ) THEN
            NK = NL
            NKM = NL**2
         ELSE
            NK = 2*NL - 1
            NKM = 2*NL**2
         END IF
C
         CALL SYMROT(IPRINT,NL,MROTR,NSYM,SYMACCEPTED,SYMUNITARY,
     &               SYMEULANG,IREL,NKM,NK,DROT,NKMMAX)
C
C-------------------- for IREL <= 2 multiple scattering is treated using
C---------------------------- the non-relativistic (l,ml)-representation
C------------------------------------------ for REAL spherical harmonics
C
         IF ( IREL.LE.2 ) THEN
C
            DO ISYM = 1,NSYM
C
               CALL CHANGEREP(NKM,NKMMAX,DROT(1,1,ISYM),'CLM>RLM',W1)
C
               DROT(1:NKM,1:NKM,ISYM) = W1(1:NKM,1:NKM)
C
               DO I = 1,NKM
                  DO J = 1,NKM
                     IF ( ABS(DIMAG(DROT(I,J,ISYM))).GT.1D-8 )
     &                    CALL STOP_MESSAGE(ROUTINE,
     &                    '|Im(DROT(I,J,ISYM)| > 1D-8')
                  END DO
               END DO
C
            END DO
C
         END IF
C
      END IF
C
      DEALLOCATE (W1)
C
      END
