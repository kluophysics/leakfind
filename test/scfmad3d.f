C*==scfmad3d.f    processed by SPAG 6.70Rc at 15:39 on 19 Dec 2016
      SUBROUTINE SCFMAD3D(IPRINT,ALAT,ABAS,QBAS,AVMAD,NQ_LOC,NQMAX)
C   ********************************************************************
C   *                                                                  *
C   *  calculate the structure dependent matrices A and B which are    *
C   *  used for the determination of the intercell potential.          *
C   *  the intercell-potential is expanded into spherical harmonics.   *
C   *  the lm-term of the intercell-potential V of the shell IQ        *
C   *  is given by                                                     *
C   *                                                                  *
C   *   V(r,lm,IQ) =  (-r)**l * SUM_{JQ,l'm'}                          *
C   *                             A(IQ,JQ,lm,l'm')*CMNTQ(JQ,l'm')      *
C   *                                                                  *
C   *  SCFMAD3D_SMAT calculates the lattice sum via Ewald technique    *
C   *                                                                  *
C   *  NOTE: the A matrices contain the factor 2 = e^2                 *
C   *        CMNTQ includes the nuclear contribution -Z/SQRT(4*pi)     *
C   *                                                                  *
C   *  based on   B. Drittler's  routines                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      USE MOD_ANGMOM,ONLY:L_LM,NRGNT123TAB
      USE MOD_SITES,ONLY:NLVMAD,NLMVMAD,NLQMAD,NLMQMAD,NLMAD,NLMMAD
      IMPLICIT NONE
C*--SCFMAD3D27
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFMAD3D')
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER IPRINT,NQMAX,NQ_LOC
      REAL*8 ABAS(3,3),AVMAD(NQMAX,NQMAX,NLMMAD,NLMMAD),QBAS(3,NQ_LOC)
C
C Local variables
C
      REAL*8 DFAC(:,:),RGAUNT,RGNTMAD(:),SMAT(:,:,:)
      REAL*8 GAUNT_RYLM
      INTEGER I,IA_ERR,IQ,IRGNTMAD(:,:),JQ,L1,L2,L3,L3MAX,LM3,LMQ,LMV,
     &        LQ,LV,M3,MQ,MV,NLM3MAX,NRGNTMAD,NRGNTMADMAX
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SMAT,DFAC,RGNTMAD,IRGNTMAD
C
      CALL TRACK_INFO(ROUTINE)
C
      IF ( NLMQMAD.LT.NLMVMAD ) STOP '<SCFMAD3D>: NLMQMAD < NLMVMAD'
C
      L3MAX = 2*(NLMAD-1)
      NLM3MAX = (L3MAX+1)**2
      NRGNTMADMAX = 2*NRGNT123TAB(NLMAD)
C
      WRITE (6,99003) NLVMAD,NLMVMAD,NLQMAD,NLMQMAD
C
      ALLOCATE (RGNTMAD(NRGNTMADMAX),IRGNTMAD(NRGNTMADMAX,3))
      ALLOCATE (DFAC(0:(NLMAD-1),0:(NLMAD-1)))
      ALLOCATE (SMAT(NQ_LOC,NQ_LOC,NLM3MAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:SCFMAD3D->DFAC'
C
C-----------------------------------------------------------------------
C
      CALL SCFMAD3D_SMAT(ALAT,ABAS,QBAS,NQ_LOC,IPRINT,NLMAD,SMAT,L3MAX,
     &                   NLM3MAX)
C
C-----------------------------------------------------------------------
C
C---> calculate:                             (2*(l+l')-1)!!
C                 DFAC(l,l') = 4pi**2 *  ----------------------
C                                        (2*l+1)!! * (2*l'+1)!!
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
C---> set up of the gaunt coefficients with an index field
C     recognize that they are needed here only for L3=LV+LQ
C
      I = 1
      DO LV = 0,(NLVMAD-1)
         DO LQ = 0,(NLQMAD-1)
            L3 = LV + LQ
            DO MV = -LV,LV
               DO MQ = -LQ,LQ
                  DO M3 = -L3,L3
C
                     RGAUNT = GAUNT_RYLM(LV,MV,LQ,MQ,L3,M3)
C
                     IF ( ABS(RGAUNT).GT.1.D-10 ) THEN
                        IF ( I.GT.NRGNTMADMAX ) THEN
                           WRITE (6,99002) I,NRGNTMADMAX
                           STOP 'in <SCFMAD3D> '
                        END IF
C
                        RGNTMAD(I) = RGAUNT
                        IRGNTMAD(I,1) = LV*(LV+1) + MV + 1
                        IRGNTMAD(I,2) = LQ*(LQ+1) + MQ + 1
                        IRGNTMAD(I,3) = L3*(L3+1) + M3 + 1
                        I = I + 1
                     END IF
C
                  END DO
               END DO
            END DO
         END DO
      END DO
      NRGNTMAD = I - 1
C
C-----------------------------------------------------------------------
C                  calculate A(IQ,JQ,LMV,LMQ)
C-----------------------------------------------------------------------
      AVMAD(1:NQMAX,1:NQMAX,1:NLMMAD,1:NLMMAD) = 0D0
C-----------------------------------------------------------------------
C
      DO JQ = 1,NQ_LOC
         DO IQ = 1,NQ_LOC
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
               AVMAD(IQ,JQ,LMV,LMQ) = AVMAD(IQ,JQ,LMV,LMQ)
     &                                + 2.0D0*DFAC(LV,LQ)
     &                                *SMAT(IQ,JQ,LM3)*RGNTMAD(I)
            END DO
         END DO
      END DO
C
      IF ( IPRINT.GE.1 ) THEN
         WRITE (6,*) ' Madelung-matrix (A(00) term)'
C
         DO IQ = 1,NQ_LOC
            WRITE (6,99001) (AVMAD(IQ,JQ,1,1),JQ=1,NQ_LOC)
         END DO
         WRITE (6,*) ' '
      END IF
C
      DEALLOCATE (SMAT,DFAC,RGNTMAD,IRGNTMAD,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:SCFMAD3D'
C
99001 FORMAT ((5X,7F10.6))
99002 FORMAT ('<SCFMAD3D>: I=',I5,' > NRGNTMADMAX =',I5)
99003 FORMAT (//,1X,79('*'),/,34X,'<SCFMAD3D>',/,1X,79('*'),//,10X,
     &        'set up Madelung matrix for 3D-bulk calculations',//,10X,
     &        'NLVMAD =',I4,5X,'NLMVMAD =',I4,/,10X,'NLQMAD =',I4,5X,
     &        'NLMQMAD =',I4,/)
C
      END
