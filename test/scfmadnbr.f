C*==scfmadnbr.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFMADNBR(AVMADNBR_Q,NQNBR_Q,RQNBR_Q,NQNBRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the structure dependent matrices A which are          *
C   *  used for the determination of the intercell potential.          *
C   *  the intercell-potential is expanded into spherical harmonics.   *
C   *  the lm-term of the intercell-potential V of the shell IQ        *
C   *  is given by                                                     *
C   *                                                                  *
C   *   V(r,lm,IQ) =  (-r)**l * SUM_{JQ,l'm'}                          *
C   *                             A(IQ,JQ,lm,l'm')*CMNTQ(JQ,l'm')      *
C   *                                                                  *
C   *   CMNTQ includes the nuclear contribution -Z/SQRT(4*pi)          *
C   *                                                                  *
C   *       l  <= NLVMAD - 1                                           *
C   *       l' <= NLQMAD - 1                                           *
C   *                                                                  *
C   *  summed over JQ (cluster sites) and l'm'                         *
C   *  using auxilary subroutine SCFMAD0D_SMAT                         *
C   *                                                                  *
C   *  NOTE: the A matrices contain the factor 2 = e^2                 *
C   *                                                                  *
C   *  in case of embedded cluster calculation only the (I-I) blocks   *
C   *  are calculated, i.e. IQ (JQ) = NQHOST+1, ..., NQHOST+NQCLU      *
C   *                                                                  *
C   ********************************************************************
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_SITES,ONLY:IQBOT,IQTOP
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      USE MOD_ANGMOM,ONLY:L_LM,NRGNT123TAB
      USE MOD_SITES,ONLY:NQMAX,NLVMAD,NLMVMAD,NLQMAD,NLMQMAD,NLMAD
      IMPLICIT NONE
C*--SCFMADNBR35
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NQNBRMAX
      REAL*8 AVMADNBR_Q(NQNBRMAX,NQMAX,NLMQMAD,NLMQMAD),
     &       RQNBR_Q(3,NQNBRMAX,NQMAX)
      INTEGER NQNBR_Q(NQMAX)
C
C Local variables
C
      REAL*8 DFAC(:,:),RGNTMAD(:),SMAT_Q(:,:,:)
      INTEGER I,IA_ERR,IQ,IRGNTMAD(:,:),JQNBR,L1,L2,L3MAX,LM3,LMQ,LMV,
     &        LQ,LV,NLM3MAX,NRGNTMAD,NRGNTMADMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SMAT_Q,DFAC,RGNTMAD,IRGNTMAD
C
      IF ( NLMQMAD.LT.NLMVMAD ) STOP '<SCFMAD0D>: NLMQMAD < NLMVMAD'
C
      L3MAX = 2*(NLMAD-1)
      NLM3MAX = (L3MAX+1)**2
      NRGNTMADMAX = NRGNT123TAB(NLMAD+1)
C
      WRITE (6,99003) NLVMAD,NLMVMAD,NLQMAD,NLMQMAD
C
      ALLOCATE (RGNTMAD(NRGNTMADMAX),IRGNTMAD(NRGNTMADMAX,3))
      ALLOCATE (DFAC(0:(NLMAD-1),0:(NLMAD-1)))
      ALLOCATE (SMAT_Q(NQNBRMAX,NQMAX,NLM3MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD0D->DFAC'
C
C-----------------------------------------------------------------------
C
      CALL SCFMADNBR_SMAT(SMAT_Q,NQNBR_Q,RQNBR_Q,NQNBRMAX,L3MAX,NLM3MAX)
C
C-----------------------------------------------------------------------
C
C                                            (2*(l+l')-1)!!
C                 dfac(l,l') = 4pi**2 *  ----------------------
C                                        (2*l+1)!! * (2*l'+1)!!
C-----------------------------------------------------------------------
C
      DFAC(0,0) = CONST_4PI*CONST_4PI
      DO L1 = 1,(NLMAD-1)
         DFAC(L1,0) = DFAC(L1-1,0)*DBLE(2*L1-1)/DBLE(2*L1+1)
         DFAC(0,L1) = DFAC(L1,0)
         DO L2 = 1,L1
            DFAC(L1,L2) = DFAC(L1,L2-1)*DBLE(2*(L1+L2)-1)/DBLE(2*L2+1)
            DFAC(L2,L1) = DFAC(L1,L2)
         END DO
      END DO
C
C-----------------------------------------------------------------------
C     set up of the Gaunt coefficients with an index field
C     recognize that they are needed here only for L3=LV+LQ
C-----------------------------------------------------------------------
C
      CALL SCFMAD_RGNT(NLVMAD,NLQMAD,RGNTMAD,IRGNTMAD,NRGNTMAD,
     &                 NRGNTMADMAX)
C
C-----------------------------------------------------------------------
C                   calculate A_Q(JQNBR,IQ,LMV,LMQ)
C-----------------------------------------------------------------------
C
      AVMADNBR_Q(:,:,:,:) = 0D0
C-----------------------------------------------------------------------
C
      DO IQ = IQBOT,IQTOP
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            DO JQNBR = 1,NQNBR_Q(IQ)
C
C---> this loop has to be calculated only for LV+LQ=L3
C
               DO I = 1,NRGNTMAD
                  LMV = IRGNTMAD(I,1)
                  LMQ = IRGNTMAD(I,2)
                  LM3 = IRGNTMAD(I,3)
                  LV = L_LM(LMV)
                  LQ = L_LM(LMQ)
C
                  AVMADNBR_Q(JQNBR,IQ,LMV,LMQ)
     &               = AVMADNBR_Q(JQNBR,IQ,LMV,LMQ) + 2.0D0*DFAC(LV,LQ)
     &               *SMAT_Q(JQNBR,IQ,LM3)*RGNTMAD(I)
               END DO
C
            END DO
C
         END IF
      END DO
C
      DEALLOCATE (SMAT_Q,DFAC,RGNTMAD,IRGNTMAD,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD0D'
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,99002)
         DO IQ = IQBOT,IQTOP
            IF ( IQREPQ(IQ).EQ.IQ ) WRITE (6,99001)
     &           (IQ,JQNBR,AVMADNBR_Q(JQNBR,IQ,1,1),JQNBR=1,NQNBR_Q(IQ))
         END DO
         WRITE (6,*) ' '
      END IF
C
99001 FORMAT ((5X,2I5,F10.6))
99002 FORMAT (//,10X,'Madelung-matrix  A  for cluster',/)
99003 FORMAT (//,1X,79('*'),/,35X,'<SCFMADNBR>',/,1X,79('*'),//,10X,
     &        'set up Madelung matrix for cluster calculations',//,10X,
     &        'NLVMAD =',I4,5X,'NLMVMAD =',I4,/,10X,'NLQMAD =',I4,5X,
     &        'NLMQMAD =',I4,/)
C
      END
C*==scfmadnbr_smat.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFMADNBR_SMAT(SMAT_Q,NQNBR_Q,RQNBR_Q,NQNBRMAX,L3MAX,
     &                          NLM3MAX)
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
      USE MOD_SYMMETRY,ONLY:IQREPQ
      USE MOD_LATTICE,ONLY:ALAT
      USE MOD_SITES,ONLY:IQBOT,IQTOP,NQMAX
      IMPLICIT NONE
C*--SCFMADNBR_SMAT183
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L3MAX,NLM3MAX,NQNBRMAX
      INTEGER NQNBR_Q(NQMAX)
      REAL*8 RQNBR_Q(3,NQNBRMAX,NQMAX),SMAT_Q(NQNBRMAX,NQMAX,NLM3MAX)
C
C Local variables
C
      REAL*8 DNRM2
      REAL*8 DRVEC(3),RFAC,RIJ,YLM(:)
      INTEGER IA_ERR,IQ,JQNBR,L,LM,M
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
      SMAT_Q(1:NQMAX,1:NQMAX,1:NLM3MAX) = 0D0
C
C=======================================================================
C
      DO IQ = IQBOT,IQTOP
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            DO JQNBR = 1,NQNBR_Q(IQ)
C
               DRVEC(1:3) = -RQNBR_Q(1:3,JQNBR,IQ)*ALAT
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
                     SMAT_Q(JQNBR,IQ,LM) = YLM(LM)*RFAC
                  END DO
               END DO
C
            END DO
C
         END IF
      END DO
C
      DEALLOCATE (YLM,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD0D_SMAT->NSR'
C
      END
