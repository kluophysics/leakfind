C*==scfmad0d.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFMAD0D
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
      USE MOD_FILES,ONLY:IPRINT
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      USE MOD_ANGMOM,ONLY:L_LM,NRGNT123TAB
      USE MOD_SITES,ONLY:AVMAD,NQHOST,NQCLU,NQMAX,NLVMAD,NLMVMAD,NLQMAD,
     &    NLMQMAD,NLMAD
      IMPLICIT NONE
C*--SCFMAD0D34
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL*8 DFAC(:,:),RGNTMAD(:),SMAT(:,:,:)
      INTEGER I,IA_ERR,IQ1,IQCLU,IQ_IMP,IRGNTMAD(:,:),JQ,JQCLU,JQ_IMP,
     &        L1,L2,L3MAX,LM3,LMQ,LMV,LQ,LV,M,NLM3MAX,NRGNTMAD,
     &        NRGNTMADMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SMAT,DFAC,RGNTMAD,IRGNTMAD
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
      ALLOCATE (SMAT(NQMAX,NQMAX,NLM3MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD0D->DFAC'
C
C-----------------------------------------------------------------------
C
      CALL SCFMAD0D_SMAT(SMAT,L3MAX,NLM3MAX)
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
C                   calculate A(IQ_IMP,JQ,LMV,LMQ)
C-----------------------------------------------------------------------
      IQ1 = NQHOST + 1
      M = NLMQMAD
      IF ( ALLOCATED(AVMAD) ) DEALLOCATE (AVMAD)
      ALLOCATE (AVMAD(IQ1:NQMAX,IQ1:NQMAX,M,M),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD0D->SMAT'
      AVMAD(IQ1:NQMAX,IQ1:NQMAX,1:M,1:M) = 0D0
C-----------------------------------------------------------------------
C
      DO IQCLU = 1,NQCLU
         IQ_IMP = NQHOST + IQCLU
C
         DO JQCLU = 1,NQCLU
C
            JQ_IMP = NQHOST + JQCLU
C
            IF ( JQCLU.NE.IQCLU ) THEN
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
                  AVMAD(IQ_IMP,JQ_IMP,LMV,LMQ)
     &               = AVMAD(IQ_IMP,JQ_IMP,LMV,LMQ) + 2.0D0*DFAC(LV,LQ)
     &               *SMAT(IQ_IMP,JQ_IMP,LM3)*RGNTMAD(I)
               END DO
C
            END IF
         END DO
      END DO
C
      DEALLOCATE (SMAT,DFAC,RGNTMAD,IRGNTMAD,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD0D'
C
      IF ( IPRINT.GT.0 ) THEN
         WRITE (6,99002)
         DO IQ_IMP = NQHOST + 1,NQHOST + NQCLU
            WRITE (6,99001) (IQ_IMP,JQ,AVMAD(IQ_IMP,JQ,1,1),JQ=NQHOST+1,
     &                      NQHOST+NQCLU)
         END DO
         WRITE (6,*) ' '
      END IF
C
99001 FORMAT ((5X,2I5,F10.6))
99002 FORMAT (//,10X,'Madelung-matrix  A  for cluster',/)
99003 FORMAT (//,1X,79('*'),/,33X,'<SCFMAD0D>',/,1X,79('*'),//,10X,
     &        'set up Madelung matrix for cluster calculations',//,10X,
     &        'NLVMAD =',I4,5X,'NLMVMAD =',I4,/,10X,'NLQMAD =',I4,5X,
     &        'NLMQMAD =',I4,/)
C
      END
