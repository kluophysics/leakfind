C*==scfcmnt.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFCMNT(IPRINT,Q1M_OCCURS,IFLAG_NEUTRALITY,CMNTIST,
     &                   CMNTISQ,CMNTMTT,CMNTMTQ,CMNTQ,QNETT,DROT_QLM,
     &                   DSYM_RLM)
C   ********************************************************************
C   *                                                                  *
C   *  Calculate charge-moments  CMNTMTQ  w.r.t. muffin tin sphere     *
C   *                            CMNTISQ  w.r.t. interstitial regime   *
C   *                            CMNTQ    total                        *
C   *                            for each site IQ                      *
C   *  In case of  CPA  this is the concentration-average of the       *
C   *  atom-type-dependent charge-moments  CMNTMTT, CMNTIST            *
C   *  otherwise:         CMNTMTQ(IQ) = CMNTMTT(IT)                    *
C   *                 and CMNTISQ(IQ) = CMNTIST(IT)                    *
C   *                                                                  *
C   *  NOTE:  CMNTQ includes the nuclear contribution -Z/SQRT(4*pi)    *
C   *         used to set the up Madelung potential in <SCFMAD_POT>    *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_RMESH,ONLY:FULLPOT
      USE MOD_SYMMETRY,ONLY:NSYM,IQREPQ,ISYMGENQ
      USE MOD_SITES,ONLY:NQ,NQMAX,ITOQ,NOQ,IQBOT,IQTOP,MAGROT_Q,KFP_LMQ,
     &    NLMQMAD,QMPHI,QMTET,QMGAM,IQAT
      USE MOD_TYPES,ONLY:TXT_T,NTMAX,NAT,CONC,KLMFP,Z,LTXT_T,ITBOT,
     &    ITTOP,NLMFPMAX,QEL
      USE MOD_CONSTANTS,ONLY:SQRT_4PI
      USE MOD_CALCMODE,ONLY:BREAKPOINT,ROTATE_LATTICE,SOLVER_FP
      USE MOD_FILES,ONLY:IOTMP
      IMPLICIT NONE
C*--SCFCMNT32
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFCMNT')
C
C Dummy arguments
C
      INTEGER IFLAG_NEUTRALITY,IPRINT
      LOGICAL Q1M_OCCURS
      REAL*8 CMNTISQ(NLMFPMAX,NQMAX),CMNTIST(NLMFPMAX,NTMAX),
     &       CMNTMTQ(NLMFPMAX,NQMAX),CMNTMTT(NLMFPMAX,NTMAX),
     &       CMNTQ(NLMFPMAX,NQMAX),DROT_QLM(NLMQMAD,NLMQMAD,IQBOT:IQTOP)
     &       ,DSYM_RLM(NLMQMAD,NLMQMAD,NSYM),QNETT(NTMAX)
C
C Local variables
C
      REAL*8 D_LMLMP,QNETQ,RSUM,SUM_IS,SUM_MT,W_IS(:),W_MT(:),ZEFF_Q(:)
      CHARACTER*40 FMT1,FMT2,FMT_QLM_BREAK
      INTEGER IO,IQ,IQREP,ISYM,IT,LM,LMP
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE W_IS,W_MT,ZEFF_Q
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (W_IS(NLMFPMAX),W_MT(NLMFPMAX))
      ALLOCATE (ZEFF_Q(NQMAX))
C
      FMT1 = '(I7,2(F14.5,''  ('',F10.5,'')''),F14.5)'
      FMT2 = '(I7,2(F14.5,14X),F14.5)'
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT X
      IF ( FULLPOT .AND. BREAKPOINT.NE.0 .AND. BREAKPOINT.NE.5 ) THEN
         FMT_QLM_BREAK = '(A,I7,2I3,2X,2I3,2X,A,2F20.14)'
         FMT1 = '(I7,2I3,2F20.14)'
         FMT2 = '(I7,2(F20.14,14X),F20.14)'
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT X
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
      IF ( FULLPOT .AND. BREAKPOINT.EQ.2 ) THEN
         CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,'QLM_charge_moment.dat')
         WRITE (IOTMP,99011) ROTATE_LATTICE
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
C
      RSUM = 0.0D0
      DO IT = ITBOT,ITTOP
         QNETT(IT) = SQRT_4PI*(CMNTMTT(1,IT)+CMNTIST(1,IT)) - Z(IT)
         QEL(IT) = SQRT_4PI*(CMNTMTT(1,IT)+CMNTIST(1,IT))
         RSUM = RSUM + QNETT(IT)*NAT(IT)*CONC(IT)
      END DO
C
      WRITE (6,99003) RSUM
C
      IF ( ABS(RSUM).GT.1D-6 ) THEN
         DO IT = ITBOT,ITTOP
            WRITE (6,99004) IT,QNETT(IT)
         END DO
         IF ( ABS(RSUM).GT.1D-2 ) IFLAG_NEUTRALITY = 1
      END IF
C
      IF ( IPRINT.GE.0 ) THEN
         DO IT = ITBOT,ITTOP
            IQ = IQAT(1,IT)
            WRITE (6,99005) 'atom type  IT =',IT,TXT_T(IT)(1:LTXT_T(IT))
            DO LM = 1,NLMQMAD
               IF ( FULLPOT .AND. BREAKPOINT.NE.0 .AND. 
     &              BREAKPOINT.NE.5 ) THEN
                  IF ( ABS(CMNTMTT(LM,IT))+ABS(CMNTIST(LM,IT)).GT.1D-6 )
     &                 WRITE (6,FMT=FMT_QLM_BREAK) 'BREAK ',LM,L_LM(LM),
     &                        M_LM(LM),NINT(QMTET(IQ)),NINT(QMPHI(IQ)),
     &                        SOLVER_FP,SQRT_4PI*CMNTMTT(LM,IT),
     &                        SQRT_4PI*CMNTIST(LM,IT)
               ELSE
                  IF ( ABS(CMNTMTT(LM,IT))+ABS(CMNTIST(LM,IT)).GT.1D-6 )
     &                 WRITE (6,FMT=FMT1) LM,SQRT_4PI*CMNTMTT(LM,IT),
     &                        CMNTMTT(LM,IT),SQRT_4PI*CMNTIST(LM,IT),
     &                        CMNTIST(LM,IT),
     &                        SQRT_4PI*(CMNTMTT(LM,IT)+CMNTIST(LM,IT))
               END IF
            END DO
         END DO
      END IF
C
C-----------------------------------------------------------------------
C                    site IQ dependent charge moments
C-----------------------------------------------------------------------
C
      IF ( IPRINT.GE.0 .AND. NQ.GT.1 ) WRITE (6,*) ' '
C
      CMNTQ(1:NLMQMAD,IQBOT:IQTOP) = 0D0
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      LOOP_IQ:DO IQ = IQBOT,IQTOP
C
         ZEFF_Q(IQ) = 0.0D0
         CMNTMTQ(1:NLMQMAD,IQ) = 0.0D0
         CMNTISQ(1:NLMQMAD,IQ) = 0.0D0
         QNETQ = 0.0D0
C
C----------------------------------------------------------- IQ = IQREPQ
         IF ( IQREPQ(IQ).EQ.IQ ) THEN
C
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
               ZEFF_Q(IQ) = ZEFF_Q(IQ) + CONC(IT)*Z(IT)
               QNETQ = QNETQ + CONC(IT)*QNETT(IT)
               DO LM = 1,NLMQMAD
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     CMNTMTQ(LM,IQ) = CMNTMTQ(LM,IQ) + CONC(IT)
     &                                *CMNTMTT(LM,IT)
                     CMNTISQ(LM,IQ) = CMNTISQ(LM,IQ) + CONC(IT)
     &                                *CMNTIST(LM,IT)
                  END IF
               END DO
            END DO
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
            IF ( FULLPOT .AND. BREAKPOINT.EQ.2 ) THEN
               WRITE (IOTMP,99009) 'before',IQ
               WRITE (IOTMP,99010) (LM,L_LM(LM),M_LM(LM),CMNTMTQ(LM,IQ),
     &                             CMNTISQ(LM,IQ),LM=1,NLMQMAD)
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
C
C--------------------------- rotate Q-moments from local to global frame
C                 use inverse of DROT_QLM  with  1/DROT_QLM = DROT_QLM^T
C
            IF ( MAGROT_Q(IQ) ) THEN
C
               W_MT(1:NLMQMAD) = CMNTMTQ(1:NLMQMAD,IQ)
               W_IS(1:NLMQMAD) = CMNTISQ(1:NLMQMAD,IQ)
C
               DO LM = 1,NLMQMAD
                  SUM_MT = 0D0
                  SUM_IS = 0D0
                  DO LMP = 1,NLMQMAD
                     D_LMLMP = DROT_QLM(LM,LMP,IQ)
                     SUM_MT = SUM_MT + W_MT(LMP)*D_LMLMP
                     SUM_IS = SUM_IS + W_IS(LMP)*D_LMLMP
                  END DO
                  CMNTMTQ(LM,IQ) = SUM_MT
                  CMNTISQ(LM,IQ) = SUM_IS
               END DO
C
            END IF
            IF ( IPRINT.GE.0 .AND. NQ.GT.1 ) WRITE (6,99008) IQ,QNETQ
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
            IF ( FULLPOT .AND. BREAKPOINT.EQ.2 ) THEN
               WRITE (IOTMP,99009) 'after',IQ
               WRITE (IOTMP,99012) QMPHI(IQ),QMTET(IQ),QMGAM(IQ)
               WRITE (IOTMP,99010) (LM,L_LM(LM),M_LM(LM),CMNTMTQ(LM,IQ),
     &                             CMNTISQ(LM,IQ),LM=1,NLMQMAD)
            END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
C
         ELSE
C---------------------------------------------------------- IQ <> IQREPQ
C
            IQREP = IQREPQ(IQ)
            IF ( IQREP.GE.IQ ) CALL STOP_MESSAGE(ROUTINE,'IQREP >= IQ')
            ZEFF_Q(IQ) = ZEFF_Q(IQREP)
            ISYM = ISYMGENQ(IQ)
            DO LM = 1,NLMQMAD
               SUM_MT = 0D0
               SUM_IS = 0D0
               DO LMP = 1,NLMQMAD
                  D_LMLMP = DSYM_RLM(LM,LMP,ISYM)
                  SUM_MT = SUM_MT + CMNTMTQ(LMP,IQREP)*D_LMLMP
                  SUM_IS = SUM_IS + CMNTISQ(LMP,IQREP)*D_LMLMP
               END DO
               CMNTMTQ(LM,IQ) = SUM_MT
               CMNTISQ(LM,IQ) = SUM_IS
            END DO
C
         END IF
C---------------------------------------------------------------- IQREPQ
C
         CMNTQ(1:NLMQMAD,IQ) = CMNTMTQ(1:NLMQMAD,IQ)
     &                         + CMNTISQ(1:NLMQMAD,IQ)
         CMNTQ(1,IQ) = CMNTQ(1,IQ) - ZEFF_Q(IQ)/SQRT_4PI
C
      END DO LOOP_IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      IF ( Q1M_OCCURS ) THEN
         WRITE (6,99002)
         DO IQ = IQBOT,IQTOP
            WRITE (6,99001) IQ,CMNTMTQ(4,IQ),CMNTMTQ(2,IQ),CMNTMTQ(3,IQ)
     &                      ,CMNTQ(4,IQ),CMNTQ(2,IQ),CMNTQ(3,IQ)
         END DO
         WRITE (6,*) ' '
      END IF
C
C=======================================================================
      IF ( IPRINT.GE.0 ) THEN
         DO IQ = IQBOT,IQTOP
            WRITE (6,99005) 'atom site  IQ =',IQ,'(including  -Z)'
            DO LM = 1,NLMQMAD
               IF ( ABS(CMNTMTQ(LM,IQ)).GT.1D-6 .OR. ABS(CMNTISQ(LM,IQ))
     &              .GT.1D-6 ) THEN
                  IF ( FULLPOT .AND. BREAKPOINT.EQ.2 ) THEN
                     WRITE (6,FMT=FMT1) LM,L_LM(LM),M_LM(LM),
     &                                  SQRT_4PI*CMNTMTQ(LM,IQ),
     &                                  SQRT_4PI*CMNTISQ(LM,IQ)
                  ELSE
                     WRITE (6,FMT=FMT2) LM,SQRT_4PI*CMNTMTQ(LM,IQ),
     &                                  SQRT_4PI*CMNTISQ(LM,IQ),
     &                                  SQRT_4PI*CMNTQ(LM,IQ)
                  END IF
                  IF ( KFP_LMQ(LM,IQ).EQ.0 ) WRITE (6,99007)
               ELSE
                  IF ( KFP_LMQ(LM,IQ).EQ.1 .AND. IPRINT.GT.0 )
     &                 WRITE (6,99006) LM
               END IF
            END DO
         END DO
      END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
      IF ( FULLPOT .AND. BREAKPOINT.EQ.2 ) CALL STOP_BREAKPOINT(ROUTINE)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 2
C
      RETURN
C
99001 FORMAT (I8,2(1X,3F10.4,:,3X))
99002 FORMAT (/,5X,'site dependent dipole moments',//21X,'muffin-tin',
     &        26X,' total'//,6X,'IQ      q_x       q_y       q_z',
     &        '           q_x       q_y       q_z')
99003 FORMAT (/,5X,'charge neutrality    ',F16.12,' els.',/)
99004 FORMAT (5X,'WARNING: for atom type IT = ',I3,'  Q(netto) = ',
     &        F12.8,' els.')
99005 FORMAT (/,5X,'charge moments for ',A,I3,3X,A,//,5X,'LM',12X,
     &        'muffin-tin',17X,'interstitial',12X,'total')
99006 FORMAT (I7,3X,
     &        'WARNING: LM-contribution expected from picking rules')
99007 FORMAT (10X,'WARNING: LM-contribution violates picking rules')
99008 FORMAT (5X,'atom site IQ =',I3,'  Q(netto) = ',F12.8,' els.')
99009 FORMAT (/,10X,'Q(l,m) charge moments ',A,' rotation IQ=',I4,/)
99010 FORMAT (3I3,2F20.10)
99011 FORMAT (10X,'BREAKPOINT in <SCFCMNT>',//,10X,'ROTATE_LATTICE = ',
     &        L3,/)
99012 FORMAT ('PHI =',F7.2,5X,'TET =',F7.2,5X,'GAM =',F7.2,/)
      END
