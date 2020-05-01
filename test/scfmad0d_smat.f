C*==scfmad0d_smat.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD0D_SMAT(SMAT,L3MAX,NLM3MAX)
C   ********************************************************************
C   *                                                                  *
C   *       calculation of terms  SMAT  l <= 2*lpot :                  *
C   *                                                                  *
C   *                        ylm( q(i) - q(j) )                        *
C   *                     -------------------------                    *
C   *                      | q(i) - q(j) |**(l+1)                      *
C   *                                                                  *
C   *       ylm       : real spherical harmic to given l,m             *
C   *       q(i),q(j) : atomic sites                                   *
C   *                                                                  *
C   *                                                                  *
C   *       for a REAL SPACE CLUSTER                                   *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:A5BAS,ALAT
      USE MOD_SITES,ONLY:NQHOST,NQCLU,N5VEC_QCLU,IQ_QCLU,QBAS,NQMAX
      IMPLICIT NONE
C*--SCFMAD0D_SMAT22
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L3MAX,NLM3MAX
      REAL*8 SMAT(NQMAX,NQMAX,NLM3MAX)
C
C Local variables
C
      REAL*8 DNRM2
      REAL*8 DRVEC(3),RFAC,RIJ,RVEC_I(3),RVEC_J(3),YLM(:)
      INTEGER I,IA_ERR,IQCLU,IQ_HOST,IQ_IMP,JQCLU,JQ_HOST,JQ_IMP,L,LM,M
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE YLM
C
C-----------------------------------------------------------------------
C               initialize real spherical harmonics
C
      ALLOCATE (YLM(NLM3MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD0D_SMAT->G'
C
C-----------------------------------------------------------------------
C
      SMAT(1:NQMAX,1:NQMAX,1:NLM3MAX) = 0D0
C
C=======================================================================
C
C---> loop over atoms in the cluster
C
      DO IQCLU = 1,NQCLU
         IQ_IMP = NQHOST + IQCLU
         IQ_HOST = IQ_QCLU(IQCLU)
C
         RVEC_I(1:3) = 0D0
         DO I = 1,5
            RVEC_I(1:3) = RVEC_I(1:3) + A5BAS(1:3,I)*N5VEC_QCLU(I,IQCLU)
         END DO
C
         RVEC_I(1:3) = RVEC_I(1:3) + QBAS(1:3,IQ_HOST)
C
         DO JQCLU = 1,NQCLU
C
            JQ_IMP = NQHOST + JQCLU
            JQ_HOST = IQ_QCLU(JQCLU)
C
            RVEC_J(1:3) = 0D0
            DO I = 1,5
               RVEC_J(1:3) = RVEC_J(1:3) + A5BAS(1:3,I)
     &                       *N5VEC_QCLU(I,JQCLU)
            END DO
C
            RVEC_J(1:3) = RVEC_J(1:3) + QBAS(1:3,JQ_HOST)
C
            DRVEC(1:3) = (RVEC_I(1:3)-RVEC_J(1:3))*ALAT
C
            IF ( JQCLU.NE.IQCLU ) THEN
C
               RIJ = DNRM2(3,DRVEC,1)
C
               DRVEC(1:3) = DRVEC(1:3)/RIJ
C
               CALL CALC_RHPLM(DRVEC(1),DRVEC(2),DRVEC(3),YLM,L3MAX,
     &                         NLM3MAX)
C
               RFAC = 1D0
               DO L = 0,L3MAX
                  RFAC = RFAC/RIJ
C
                  DO M = -L,L
                     LM = L*(L+1) + M + 1
                     SMAT(IQ_IMP,JQ_IMP,LM) = YLM(LM)*RFAC
                  END DO
               END DO
C
            END IF
         END DO
      END DO
C
      DEALLOCATE (YLM,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD0D_SMAT->NSR'
C
      END
