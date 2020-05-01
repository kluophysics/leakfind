C*==tauij_std.f    processed by SPAG 6.70Rc at 22:17 on 20 Dec 2016
      SUBROUTINE TAUIJ_STD(ERYD,P,TAUQ,MSSQ)
C   ********************************************************************
C   *                                                                  *
C   *  perform the Brillouin zone integration to get  TAU_ij           *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_LATTICE,ONLY:ABAS
      USE MOD_SITES,ONLY:NQMAX,NQ
      USE MOD_ANGMOM,ONLY:NKMMAX,NKKR,IND0Q,NKM,NLM,NXM_Q,NXM,NSPIN
      USE MOD_CALCMODE,ONLY:KMROT,IREL
      USE MOD_CONSTANTS,ONLY:C1,CI2PI
      USE MOD_TAUIJ,ONLY:NTAUIJ,NTAUIJ_QTAB2_TAUIJ,IQ_QTAB1_TAUIJ,
     &    NQTAB1_TAUIJ,N5VEC_TAUIJ,N123TAUIJMAX,TAUIJ,NKTABTAUIJ,
     &    IK_IKTABTAUIJ,ISYM_IKTABTAUIJ,KTABTAUIJ,JQ_QTAB2_TAUIJ,
     &    NQTAB2_TAUIJ
      USE MOD_MPI,ONLY:MPI_ID,MPI_KLOOP,MPI_ELOOP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 TOL
      PARAMETER (TOL=1D-8)
      LOGICAL CHECK
      PARAMETER (CHECK=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      COMPLEX*16 CWORK(:,:),EXIK1,EXIK2,EXIK3,EXIKRS,MAUX(:,:),
     &           PWEXIK1(-N123TAUIJMAX:N123TAUIJMAX),
     &           PWEXIK2(-N123TAUIJMAX:N123TAUIJMAX),
     &           PWEXIK3(-N123TAUIJMAX:N123TAUIJMAX),T1,T2,TAUK(:,:)
      REAL*8 DDOT
      INTEGER I,I1,IA_ERR,IFLAG,IK,INFO,IPIV(:),IPROCK(:),IPW,IQ,IQTAB,
     &        IS,ITAUIJ,ITAUIJ1,ITAUIJ2,J,J0,J1,JQ,JQTAB,LMSOFF,M,N1,N2,
     &        N3,NI,NJ,NSUM
      REAL*8 NORM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUK,IPIV,MAUX,IPROCK,CWORK
C
      IF ( KMROT.NE.0 ) STOP 'in <TAUIJ_STD>:  KMROT <> 0'
C
      ALLOCATE (IPIV(NKKR),CWORK(NXM,NXM))
      ALLOCATE (TAUK(NKKR,NKKR),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:TAUIJ_STD->MAUX'
C
C-----------------------------------------------------------------------
C                           initialize
C-----------------------------------------------------------------------
C
      IF ( .NOT.ALLOCATED(TAUIJ) ) THEN
         ALLOCATE (TAUIJ(NXM,NXM,NSPIN,NTAUIJ),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc: TAUIJ_STD->TAUIJ'
      END IF
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C      distribute k-points IK = 1, ..., NKTABTAUIJ
C      over processors;   IK=NKTABTAUIJ  dealt with by IPROC=0
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      ALLOCATE (IPROCK(NKTABTAUIJ))
C
      IPROCK(1:NKTABTAUIJ) = 0
C
      CALL MPI_DISTRIBUTE(IPROCK,NKTABTAUIJ,MPI_KLOOP,'K')
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
C=======================================================================
C  calculate spin-dependent   TAU[i,j]   via Brillouin zone integration
C=======================================================================
C
      CALL CINIT(NXM*NXM*NSPIN*NTAUIJ,TAUIJ)
C
      PWEXIK1(0) = 1.0D0
      PWEXIK2(0) = 1.0D0
      PWEXIK3(0) = 1.0D0
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO IS = 1,NSPIN
C
         IF ( IREL.EQ.2 ) THEN
            LMSOFF = NLM*(IS-1)
         ELSE
            LMSOFF = 0
         END IF
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
         DO IK = 1,NKTABTAUIJ
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
            IF ( MPI_ID.NE.IPROCK(IK) .AND. .NOT.MPI_ELOOP ) CYCLE
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
            CALL STRSET(IK,KTABTAUIJ(1,IK),MAUX,TAUK,P)
C
            DO IQ = 1,NQ
               NI = NXM_Q(IQ)
               I1 = IND0Q(IQ) + 1
               DO J = 1,NI
                  J1 = IND0Q(IQ) + J
                  CALL ZAXPY(NI,C1,MSSQ(LMSOFF+1,LMSOFF+J,IQ),1,
     &                       TAUK(I1,J1),1)
               END DO
            END DO
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
C -------------------------- table to calculate   exp( i ->k * ->R[ij] )
C
            EXIK1 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,1),1))
            EXIK2 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,2),1))
            EXIK3 = CDEXP(-CI2PI*DDOT(3,KTABTAUIJ(1,IK),1,ABAS(1,3),1))
C
            DO IPW = 1,N123TAUIJMAX
               PWEXIK1(IPW) = PWEXIK1(IPW-1)*EXIK1
               PWEXIK1(-IPW) = 1.0D0/PWEXIK1(IPW)
C
               PWEXIK2(IPW) = PWEXIK2(IPW-1)*EXIK2
               PWEXIK2(-IPW) = 1.0D0/PWEXIK2(IPW)
C
               PWEXIK3(IPW) = PWEXIK3(IPW-1)*EXIK3
               PWEXIK3(-IPW) = 1.0D0/PWEXIK3(IPW)
            END DO
C
C------------------------------------------------------------ store TAUQ
C
            ITAUIJ2 = 0
            DO IQTAB = 1,NQTAB1_TAUIJ
               IQ = IQ_QTAB1_TAUIJ(IQTAB)
               NI = NXM_Q(IQ)
               I1 = IND0Q(IQ) + 1
               DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
                  JQ = JQ_QTAB2_TAUIJ(IQTAB,JQTAB)
                  NJ = NXM_Q(JQ)
                  J0 = IND0Q(JQ)
                  ITAUIJ1 = ITAUIJ2 + 1
                  ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
                  DO ITAUIJ = ITAUIJ1,ITAUIJ2
                     N1 = N5VEC_TAUIJ(1,ITAUIJ)
                     N2 = N5VEC_TAUIJ(2,ITAUIJ)
                     N3 = N5VEC_TAUIJ(3,ITAUIJ)
                     EXIKRS = PWEXIK1(N1)*PWEXIK2(N2)*PWEXIK3(N3)
                     DO J = 1,NJ
                        J1 = J0 + J
                        CALL ZAXPY(NI,EXIKRS,TAUK(I1,J1),1,
     &                             TAUIJ(1,J,IS,ITAUIJ),1)
                     END DO
                  END DO
               END DO
            END DO
C-----------------------------------------------------------------------
C
         END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
      IF ( MPI_KLOOP ) THEN
C
         CALL DRV_MPI_BARRIER
C
         M = NXM*NXM
         DO IS = 1,NSPIN
            DO ITAUIJ = 1,NTAUIJ
               CALL DRV_MPI_REDUCE_C(TAUIJ(1,1,IS,ITAUIJ),CWORK(1,1),M)
            END DO
         END DO
      END IF
C MPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIMPIM
C
      NORM = 1D0/DBLE(NKTABTAUIJ)
C
      TAUIJ(1:NXM,1:NXM,1:NSPIN,1:NTAUIJ)
     &   = NORM*TAUIJ(1:NXM,1:NXM,1:NSPIN,1:NTAUIJ)
C
C=======================================================================
      IF ( .NOT.CHECK ) RETURN
C=======================================================================
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      DO IS = 1,NSPIN
C
         IF ( IREL.EQ.2 ) THEN
            LMSOFF = NLM*(IS-1)
         ELSE
            LMSOFF = 0
         END IF
C
         ITAUIJ2 = 0
         DO IQTAB = 1,NQTAB1_TAUIJ
            IQ = IQ_QTAB1_TAUIJ(IQTAB)
            NI = NXM_Q(IQ)
            DO JQTAB = 1,NQTAB2_TAUIJ(IQTAB)
               JQ = JQ_QTAB2_TAUIJ(IQTAB,JQTAB)
               NJ = NXM_Q(JQ)
               ITAUIJ1 = ITAUIJ2 + 1
               ITAUIJ2 = ITAUIJ2 + NTAUIJ_QTAB2_TAUIJ(IQTAB,JQTAB)
               DO ITAUIJ = ITAUIJ1,ITAUIJ2
                  N1 = N5VEC_TAUIJ(1,ITAUIJ)
                  N2 = N5VEC_TAUIJ(2,ITAUIJ)
                  N3 = N5VEC_TAUIJ(3,ITAUIJ)
                  NSUM = ABS(N1) + ABS(N2) + ABS(N3)
                  IF ( IQ.EQ.JQ .AND. NSUM.EQ.0 ) THEN
C
                     IFLAG = 0
                     WRITE (6,99002) IQ,IS
                     DO I = 1,NI
                        DO J = 1,NJ
                           T1 = TAUIJ(I,J,IS,ITAUIJ)
                           T2 = TAUQ(LMSOFF+I,LMSOFF+J,IQ)
                           IF ( ABS(T1-T2).GT.1D-8 ) THEN
                              WRITE (6,99001) IQ,IS,I,J,T1,T2
                              IFLAG = 1
                           END IF
                        END DO
                     END DO
C
                     IF ( IFLAG.NE.0 )
     &                    CALL CMATSTRUCT('TAU[q]',TAUQ(1,1,IQ),NKM,
     &                    NKMMAX,IREL,IREL,1,TOL,6)
C
                  END IF
               END DO
            END DO
         END DO
C
      END DO
C
      STOP
C-----------------------------------------------------------------------
99001 FORMAT (2I3,3X,2I3,2E18.10,/,15X,2E18.10,/)
99002 FORMAT (/,10X,'checking TAU for site IQ =',I3,'   IS =',I3,/)
      END
