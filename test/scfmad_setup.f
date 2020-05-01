C*==scfmad_setup.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFMAD_SETUP(FULLPOT5)
C   ********************************************************************
C   *                                                                  *
C   *              set up the Madelung matrices                        *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ALAT,ABAS,SUB_SYSTEM,ABAS_L,ABAS_R,
     &    SYSTEM_DIMENSION,SYSTEM_TYPE
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_SITES,ONLY:NQMAX,QBAS,NQ_L,NQ_R,NQ_I,NLMAD,NLMMAD,NLQMAD,
     &    NLMQMAD,NLVMAD,NLMVMAD,AVMAD,IQ_LREF,IQ_RREF,NQ,KFP_LMQ,NQCLU,
     &    NQHOST,NOQ,ITOQ
      USE MOD_TYPES,ONLY:KLMFP,NTMAX,NLMFPMAX
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_ANGMOM,ONLY:NL,NLMMADMAX
      IMPLICIT NONE
C*--SCFMAD_SETUP19
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD_SETUP')
C
C Dummy arguments
C
      LOGICAL FULLPOT5
C
C Local variables
C
      REAL*8 AVMADTMP(:,:,:,:)
      INTEGER IA_ERR,IFLAG,IFLAG_Q,IO,IQ,IQ1,IQ1_R,IQ2,IQ2_R,IT,JQ,JQ1,
     &        JQ2,LMQ,LMV,NREPBAS_L,NREPBAS_R
      LOGICAL KDIPOLE_Q(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE KDIPOLE_Q,AVMADTMP
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( ALLOCATED(AVMAD) ) DEALLOCATE (AVMAD)
C
      ALLOCATE (KDIPOLE_Q(NQMAX))
      KDIPOLE_Q(1:NQMAX) = .FALSE.
C
      JQ1 = 999999
      JQ2 = 999999
C=======================================================================
C                  set up the Madelung matrices
C=======================================================================
C
      IF ( FULLPOT .OR. FULLPOT5 ) THEN
C
         NLMAD = 2*(NL-1) + 1
C
      ELSE
C
         NLMAD = 2
C
         CALL AMEY1M
C
      END IF
C
      NLQMAD = NLMAD
      NLVMAD = NLMAD
C
      NLMMAD = NLMAD**2
      NLMVMAD = NLVMAD**2
      NLMQMAD = NLQMAD**2
C
      NLMMADMAX = NLMMAD
C
C***********************************************************************
C                        HOST calculation
C***********************************************************************
C
      IF ( SYSTEM_TYPE(1:16).NE.'EMBEDDED-CLUSTER' ) THEN
C
C=======================================================================
         SELECT CASE (SYSTEM_DIMENSION(1:2))
C
C=======================================================================
         CASE ('3D')
C=======================================================================
C
            ALLOCATE (AVMAD(NQMAX,NQMAX,NLMMADMAX,NLMMADMAX),
     &                STAT=IA_ERR)
            IF ( IA_ERR.NE.0 )
     &            CALL STOP_MESSAGE(ROUTINE,'ALLOC: AVMAD A')
C
            CALL SCFMAD3D(IPRINT,ALAT,ABAS,QBAS,AVMAD,NQ,NQMAX)
C
            IQ1 = 1
            IQ2 = NQ
            JQ1 = IQ1
            JQ2 = IQ2
C
C=======================================================================
         CASE ('2D')
C=======================================================================
C
            SELECT CASE (SUB_SYSTEM(1:6))
C-----------------------------------------------------------------------
C              Left bulk system in case of  SYSTEM_TYPE=LIR
C-----------------------------------------------------------------------
            CASE ('L-BULK')
C
               ALLOCATE (AVMAD(NQMAX,NQMAX,NLMMADMAX,NLMMADMAX),
     &                   STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: AVMAD B')
C
               CALL SCFMAD3D(IPRINT,ALAT,ABAS_L,QBAS(1,1),AVMAD,NQ_L,
     &                       NQMAX)
C
               IQ1 = 1
               IQ2 = NQ
               JQ1 = IQ1
               JQ2 = IQ2
C
C-----------------------------------------------------------------------
C           Right bulk system in case of  SYSTEM_TYPE=LIR
C-----------------------------------------------------------------------
            CASE ('R-BULK')
C
               ALLOCATE (AVMAD(NQMAX,NQMAX,NLMMADMAX,NLMMADMAX),
     &                   STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: AVMAD C')
C
               CALL SCFMAD3D(IPRINT,ALAT,ABAS_R,QBAS(1,NQ_L+NQ_I+1),
     &                       AVMAD,NQ_R,NQMAX)
C
               IQ1_R = NQ_L + NQ_I + 1
               IQ2_R = NQ_L + NQ_I + NQ_R
C
               AVMAD(IQ1_R:IQ2_R,IQ1_R:IQ2_R,1:NLMMADMAX,1:NLMMADMAX)
     &            = AVMAD(1:NQ_R,1:NQ_R,1:NLMMADMAX,1:NLMMADMAX)
C
               IQ1 = IQ1_R
               IQ2 = IQ2_R
               JQ1 = IQ1
               JQ2 = IQ2
C
C-----------------------------------------------------------------------
C                SYSTEM_TYPE=LIR   or  slab (SYSTEM_TYPE=VIV)
C-----------------------------------------------------------------------
            CASE ('I-ZONE')
C
               ALLOCATE (AVMADTMP(NQMAX,NQMAX,NLMMADMAX,NLMMADMAX),
     &                   STAT=IA_ERR)
               IF ( IA_ERR.NE.0 )
     &               CALL STOP_MESSAGE(ROUTINE,'ALLOC: AVMAD D')
               ALLOCATE (AVMAD(NQMAX,NQMAX,NLMMADMAX,NLMMADMAX),
     &                   STAT=IA_ERR)
C
               AVMAD(1:NQMAX,1:NQMAX,1:NLMMADMAX,1:NLMMADMAX) = 0D0
C
C------------------------------------------ calculate R-R-block of AVMAD
C
               CALL SCFMAD3D(IPRINT,ALAT,ABAS_R,QBAS(1,NQ_L+NQ_I+1),
     &                       AVMADTMP,NQ_R,NQMAX)
C
               IQ1_R = NQ_L + NQ_I + 1
               IQ2_R = NQ_L + NQ_I + NQ_R
C
               AVMAD(IQ1_R:IQ2_R,IQ1_R:IQ2_R,1:NLMMADMAX,1:NLMMADMAX)
     &            = AVMADTMP(1:NQ_R,1:NQ_R,1:NLMMADMAX,1:NLMMADMAX)
C
C------------------------------------------ calculate L-L-block of AVMAD
C
               CALL SCFMAD3D(IPRINT,ALAT,ABAS_L,QBAS,AVMADTMP,NQ_L,
     &                       NQMAX)
C
               AVMAD(1:NQ_I,1:NQ_I,1:NLMMADMAX,1:NLMMADMAX)
     &            = AVMADTMP(1:NQ_I,1:NQ_I,1:NLMMADMAX,1:NLMMADMAX)
               NREPBAS_L = 10
               NREPBAS_R = 10
C
               IQ_LREF = NQ_L
               IQ_RREF = NQ_L + NQ_I + 1
C
C
               CALL SCFMAD2D(IQ_LREF,IQ_RREF,NREPBAS_L,NREPBAS_R)
C
C=======================================================================
C                         ASA calculation
C=======================================================================
C
               IF ( FULLPOT .OR. FULLPOT5 ) RETURN
C
C-------------------------- assume dipole moments to occur for all sites
C
C
               KLMFP(1:NLMFPMAX,1:NTMAX) = 0
               KLMFP(1,1:NTMAX) = 1
               KFP_LMQ(1:NLMFPMAX,1:NQMAX) = 0
               KFP_LMQ(1,1:NQMAX) = 1
C
               DO IQ = NQ_L + 1,NQ_L + NQ_I
                  DO LMV = 1,4
                     KFP_LMQ(LMV,IQ) = 1
                  END DO
C
                  DO IO = 1,NOQ(IQ)
                     IT = ITOQ(IO,IQ)
                     DO LMV = 2,4
                        KLMFP(LMV,IT) = KFP_LMQ(LMV,IQ)
                     END DO
                  END DO
               END DO
C
               RETURN
C
C-----------------------------------------------------------------------
            CASE DEFAULT
C
               CALL STOP_MESSAGE(ROUTINE,'CASE SELECT SUB_SYSTEM')
C
            END SELECT
C
C=======================================================================
         CASE ('0D')
C=======================================================================
C
C-----------------------------------------------------------------------
C                 free cluster
C-----------------------------------------------------------------------
C
            CALL SCFMAD0D
C
            IQ1 = 1
            IQ2 = NQ
            JQ1 = IQ1
            JQ2 = IQ2
C
C=======================================================================
         CASE DEFAULT
C
            CALL STOP_MESSAGE(ROUTINE,'CASE SELECT SYSTEM_DIMENSION')
C
         END SELECT
C
C***********************************************************************
C                 cluster embedded in host system
C***********************************************************************
      ELSE
C
         CALL SCFMAD0D_BACK
C
         CALL SCFMAD0D
C
         IQ1 = NQHOST + 1
         IQ2 = NQHOST + NQCLU
         JQ1 = IQ1
         JQ2 = IQ2
C
      END IF
C***********************************************************************
C
C
C
C***********************************************************************
C
C                         ASA calculation
C
C  check whether dipole moments (1,m) occur
C  indicated by   AVMAD(IQ,JQ,(1,m),LM2) <> 0 for m=-1,0,+1
C  information is stored by the flags  KLMFP  and  KFP_LMQ
C  set  NLQMAD =1 or 4 accordingly
C***********************************************************************
C
C=======================================================================
C                      3D or 0D (sub)systems treated
C=======================================================================
C
      IF ( FULLPOT .OR. FULLPOT5 ) RETURN
C
      KLMFP(1:NLMFPMAX,1:NTMAX) = 0
      KLMFP(1,1:NTMAX) = 1
      KFP_LMQ(1:NLMFPMAX,1:NQMAX) = 0
      KFP_LMQ(1,1:NQMAX) = 1
C
      IFLAG = 0
      DO IQ = IQ1,IQ2
         IFLAG_Q = 0
         DO JQ = JQ1,JQ2
            DO LMV = 2,4
               DO LMQ = 1,NLMMAD
                  IF ( ABS(AVMAD(IQ,JQ,LMV,LMQ)).GT.1D-10 ) THEN
                     KFP_LMQ(LMV,IQ) = 1
                     IFLAG_Q = 1
                     KDIPOLE_Q(IQ) = .TRUE.
                  END IF
               END DO
            END DO
         END DO
C
         DO IO = 1,NOQ(IQ)
            IT = ITOQ(IO,IQ)
            DO LMV = 2,4
               KLMFP(LMV,IT) = KFP_LMQ(LMV,IQ)
            END DO
         END DO
         IFLAG = IFLAG + IFLAG_Q
      END DO
C
C------------------------------- reduce NLMAD if no dipole moments occur
      IF ( IFLAG.EQ.0 ) THEN
         NLMAD = 1
         NLMMAD = 1
         NLQMAD = 1
         NLMQMAD = 1
         NLVMAD = 1
         NLMVMAD = 1
         WRITE (6,99001) NLVMAD,NLMVMAD,NLQMAD,NLMQMAD
      ELSE
         WRITE (6,99002) (IQ,KDIPOLE_Q(IQ),IQ=IQ1,IQ2)
      END IF
C
99001 FORMAT (//,1X,79('*'),/,33X,'<SCFMAD_SETUP>',/,1X,79('*'),//,10X,
     &        'for ASA  NO  dipole moments expected from AVMAD',//,10X,
     &        'angular momentum parameters reset to:',//,10X,'NLVMAD =',
     &        I4,5X,'NLMVMAD =',I4,/,10X,'NLQMAD =',I4,5X,'NLMQMAD =',
     &        I4,/)
99002 FORMAT (//,1X,79('*'),/,33X,'<SCFMAD_SETUP>',/,1X,79('*'),//,10X,
     &        'for ASA   FINITE  dipole moments expected from AVMAD ',
     &        'for sites IQ:',//,(10X,8(I3,L2,3X)))
      END
