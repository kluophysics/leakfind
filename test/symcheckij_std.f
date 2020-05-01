C*==symcheckij_std.f    processed by SPAG 6.70Rc at 22:17 on 20 Dec 2016
      SUBROUTINE SYMCHECKIJ_STD(ERYD,P,ICPAFLAG,CPACHNG,IPRINT,ITCPA,
     &                          ICPACONV,TSST,MSST,TSSQ,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *  PART 1:                                                         *
C   *                                                                  *
C   *  compare for each symmetry operation R:                          *
C   *          TAU(Rk) = R TAU(k) R+                                   *
C   *            G(Rk) = R   G(k) R+                                   *
C   *            m     = R   m    R+                                   *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  PART 4:                                                         *
C   *                                                                  *
C   *  <TAU_STD_KPOINTS> is run WITH and WITHOUT                       *
C   *                    making use of symmetry                        *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *  PART 5:                                                         *
C   *                                                                  *
C   *  check the relation   TAU(q',q') = R TAU(q,q) R+   for q' = R q  *
C   *                                                                  *
C   * 04/07/2000 HE - update 2004                                      *
C   ********************************************************************
C
      USE MOD_TYPES,ONLY:NTMAX,CONC
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_KSPACE,ONLY:NELMT
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NKKR,NKM
      USE MOD_SYMMETRY,ONLY:NWEDGE,SYMUNITARY,SYMDET,IQORGQP,DROT,
     &    NSYMACCEPTED,SYMACCEPTED,MROTK,MROTR,NSYM
      USE MOD_SITES,ONLY:QBAS,NQMAX,NQ,ITOQ,NOQ
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET
      USE MOD_CONSTANTS,ONLY:C0,C1,CI2PI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SYMCHECKIJ_STD')
      REAL*8 TOLR,TOLM
      PARAMETER (TOLR=1.0D-8,TOLM=1.0D-8)
      INTEGER NSTEP
      PARAMETER (NSTEP=1)
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD,P
      INTEGER ICPACONV,ICPAFLAG,IPRINT,ITCPA
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TSSQ(NKMMAX,NKMMAX,NQMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      CHARACTER*1 CNT
      REAL*8 DDOT
      REAL*8 DQVEC(3),KTAB(:,:),KVEC(3),KVECP(3),MVEC(3),MVECP(3),
     &       RVEC(3),RVECP(3),UKR,WKTAB(:)
      COMPLEX*16 EMIUKR,MAUX(:,:),PHASK(1),TAUK(:,:),TAUQA(:,:,:),
     &           TAUQB(:,:,:),TAUQQ(:,:,:,:),W1(:,:),W2(:,:)
      INTEGER I,I1,IA_ERR,IE,IK,INFO,IPIV(:),IQ,IQP,ISYM,IX,J,J1,JQ,JQP,
     &        M,NERR4,NERR5,NERRT,NERRTS,NI,NKPTS0,NKTAB,NKTABMAX,NSYM0,
     &        NSYMDUM
      CHARACTER*5 STRERR
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MAUX,TAUK,IPIV,TAUQQ,TAUQA,TAUQB,W1,W2,KTAB,WKTAB
C
      M = NKMMAX
C
      ALLOCATE (TAUQA(M,M,NQMAX),W1(M,M),W2(M,M),IPIV(NKKR))
      ALLOCATE (TAUQQ(M,M,NQMAX,NQMAX),TAUQB(M,M,NQMAX))
      ALLOCATE (MAUX(NKKR,NKKR))
      ALLOCATE (TAUK(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUQB')
C
      WRITE (6,99015)
C
      NSYM0 = NSYM
C
C=======================================================================
      RVEC(1) = 0.412D0
      RVEC(2) = 0.231D0
      RVEC(3) = 0.123D0
      KVEC(1) = 0.412D0
      KVEC(2) = 0.231D0
      KVEC(3) = 0.123D0
      MVEC(1) = 0.0D0
      MVEC(2) = 0.0D0
      MVEC(3) = 1.0D0
C=======================================================================
C                               PART  0
C=======================================================================
C
      WRITE (6,99001) 0,' check results of symmetry operations R'
      WRITE (6,99007)
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,RVEC,1,0D0,RVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
            CALL DGEMV('N',3,3,1D0,MROTR(1,1,ISYM),3,MVEC,1,0D0,MVECP,1)
            CALL DSCAL(3,DBLE(SYMDET(ISYM)),MVECP,1)
C
            WRITE (6,99003) ISYM,SYMUNITARY(ISYM),'M',MVECP,'R',RVECP,
     &                      'k',KVECP
C
         END IF
C
      END DO
C
C=======================================================================
C                               PART  1
C=======================================================================
C
      WRITE (6,99001) 1,' check results of symmetry operations R'
      WRITE (6,99008) KVEC
C
C=======================================================================
C ------------------ calculate energy - dependent terms of str.constants
C
      CALL STRCC(ERYD,.FALSE.)
C
      CALL CINIT(M*M*NQ,TAUQ)
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      NERRTS = 0
      WRITE (6,*) ' '
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
C-NOTE:  MROTK includes the INVERSION for ANTI - unitary transformations
C
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
C
C-----------------------------------------------------------------------
C                  calculate   TAU(Rk)
C-----------------------------------------------------------------------
C
            CALL STRSET(1,KVECP,TAUK,MAUX,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
            DO IQ = 1,NQ
               I1 = IND0Q(IQ) + 1
               NI = NKMQ(IQ)
               DO JQ = 1,NQ
                  DO J = 1,NKMQ(JQ)
                     J1 = IND0Q(JQ) + J
                     CALL ZCOPY(NI,TAUK(I1,J1),1,TAUQQ(1,J,IQ,JQ),1)
                  END DO
               END DO
            END DO
C
C-----------------------------------------------------------------------
C                  calculate  exp(-i Uk(R_q-R_q')) R TAU(k) R+
C-----------------------------------------------------------------------
C
            CALL DGEMV('N',3,3,1D0,MROTK(1,1,ISYM),3,KVEC,1,0D0,KVECP,1)
C
            CALL STRSET(1,KVECP,TAUK,MAUX,P)
C
            CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,M)
C
            CALL ZGETRF(NKKR,NKKR,TAUK,NKKR,IPIV,INFO)
            CALL ZGETRI(NKKR,TAUK,NKKR,IPIV,MAUX,NKKR*NKKR,INFO)
C
            NERRT = 0
            DO IQ = 1,NQ
               IQP = IQORGQP(ISYM,IQ)
C
               DO JQ = 1,NQ
                  JQP = IQORGQP(ISYM,JQ)
C
                  DQVEC(1:3) = QBAS(1:3,IQP) - QBAS(1:3,JQP)
                  UKR = DDOT(3,KVECP,1,DQVEC,1)
C
                  EMIUKR = EXP(-CI2PI*UKR)
                  EMIUKR = 1D0
C
                  WRITE (6,*) '********* EMIUKR ',EMIUKR
C
                  I1 = IND0Q(IQP) + 1
                  NI = NKMQ(IQP)
                  DO J = 1,NKMQ(JQP)
                     J1 = IND0Q(JQP) + J
                     CALL ZCOPY(NI,TAUK(I1,J1),1,W1(1,J),1)
                  END DO
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
                  IF ( SYMUNITARY(ISYM) ) THEN
                     CNT = 'N'
                  ELSE
                     CNT = 'T'
                  END IF
C-----------------------------------------------------------------------
                  CALL ZGEMM('N',CNT,NKM,NKM,NKM,EMIUKR,DROT(1,1,ISYM),
     &                       M,W1,M,C0,W2,M)
C
                  CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),
     &                       M,C0,W1,M)
C-----------------------------------------------------------------------
C
                  DO J = 1,NKM
                     DO I = 1,NKM
                        IF ( ABS(W1(I,J)-TAUQQ(I,J,IQ,JQ)).GT.TOLR )
     &                       THEN
                           NERRT = NERRT + 1
                           WRITE (6,99002) ISYM,IQ,JQ,I,J,'  TAU(R k) ',
     &                            TAUQQ(I,J,IQ,JQ),'R TAU(k) R+',W1(I,J)
     &                            ,W1(I,J) - TAUQQ(I,J,IQ,JQ)
                        END IF
                     END DO
                  END DO
C
               END DO
            END DO
C
C-----------------------------------------------------------------------
C
            NERRTS = NERRTS + NERRT
C
            IF ( NERRT.GT.0 ) THEN
               STRERR = 'ERROR'
            ELSE
               STRERR = 'OK   '
            END IF
            WRITE (6,99006) ISYM,SYMUNITARY(ISYM),'K',KVECP,NERRT,STRERR
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      WRITE (6,99010) 1,NERRTS
C
      IF ( NSTEP.LE.1 ) CALL STOP_MESSAGE(ROUTINE,'NSTEP.LE.1')
C
C
C=======================================================================
C                               PART  4
C=======================================================================
C
      WRITE (6,99001) 4,'   check integration of TAU - matrix'
      WRITE (6,99004)
C
      IE = 1
      NERR4 = 0
      NKPTS0 = 100
C
      CALL KMESHS(IOTMP,NSYM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB')
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
      CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                     ICPACONV,CONC,NOQ,ITOQ,PHASK,IE,NTMAX,TSST,
     &                     MSST,TSSQ,MSSQ,TAUQA)
C
      WRITE (6,99005)
      NSYMDUM = 1
C
      CALL KMESHS(IOTMP,NSYMDUM,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      DEALLOCATE (KTAB,WKTAB)
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: KTAB')
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
      NSYM = NSYMDUM
C
      CALL TAU_STD_KPOINTS(ICPAFLAG,CPACHNG,ERYD,P,IPRINT,ITCPA,
     &                     ICPACONV,CONC,NOQ,ITOQ,PHASK,IE,NTMAX,TSST,
     &                     MSST,TSSQ,MSSQ,TAUQB)
C
      NSYM = NSYM0
C
      DO IQ = 1,NQ
C
         WRITE (6,99009) 4,IQ
C
         DO J = 1,NKMQ(IQ)
            DO I = 1,NKMQ(IQ)
               TAUQB(I,J,IQ) = TAUQB(I,J,IQ)*NSYMACCEPTED - 
     &                         TAUQA(I,J,IQ)
               IF ( ABS(TAUQB(I,J,IQ)).GT.TOLM ) NERR4 = NERR4 + 1
            END DO
         END DO
C
         CALL CMATSTRUCT('TAU(IQ)  run A   ',TAUQA(1,1,IQ),NKMQ(IQ),M,
     &                   IREL,IREL,0,TOLM,6)
C
         CALL CMATSTRUCT('TAU(IQ) run B - A  SHOULD BE BLANK',
     &                   TAUQB(1,1,IQ),NKMQ(IQ),M,IREL,IREL,0,TOLM,6)
C
      END DO
C
      WRITE (6,99010) 4,NERR4
C
C=======================================================================
C                               PART  5
C=======================================================================
C
      WRITE (6,99001) 5,'check symmetry of site-diagonal TAU(q,q)'
      WRITE (6,99016)
C
      NERR5 = 0
C
      DO ISYM = 1,NSYM
         IF ( SYMACCEPTED(ISYM) ) THEN
C
            DO IQ = 1,NQ
C
               IQP = IQORGQP(ISYM,IQ)
C
C-----------------------------------------------------------------------
C                 unitary / ANTI - unitary transformation
C-----------------------------------------------------------------------
               IF ( SYMUNITARY(ISYM) ) THEN
                  CNT = 'N'
               ELSE
                  CNT = 'T'
               END IF
C-----------------------------------------------------------------------
               CALL ZGEMM('N',CNT,NKM,NKM,NKM,C1,DROT(1,1,ISYM),M,
     &                    TAUQA(1,1,IQ),M,C0,W2,M)
C
               CALL ZGEMM('N','C',NKM,NKM,NKM,C1,W2,M,DROT(1,1,ISYM),M,
     &                    C0,W1,M)
C-----------------------------------------------------------------------
C
               DO J = 1,NKM
                  DO I = 1,NKM
                     IF ( ABS(W1(I,J)-TAUQA(I,J,IQP)).GT.TOLR ) THEN
                        NERR5 = NERR5 + 1
                        WRITE (6,99002) ISYM,IQ,I,J,'  TAU(q'',q'') ',
     &                                  TAUQA(I,J,IQ),'R TAU(q,q) R+',
     &                                  W1(I,J),W1(I,J) - TAUQA(I,J,IQ)
                     END IF
                  END DO
               END DO
C
            END DO
C
C-----------------------------------------------------------------------
C
         END IF
C
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      WRITE (6,99010) 5,NERR5
C
C=======================================================================
C
      WRITE (6,99011) NQ,KMROT,NSYM,NSYMACCEPTED,NWEDGE
      IF ( KMROT.EQ.0 ) WRITE (6,99012) NELMT
      WRITE (6,99013) NERRTS
      WRITE (6,99014) NERR4,NERR5
C
C      STOP 'TEST 2 completed'
C
C=======================================================================
C
      DEALLOCATE (MAUX,TAUK,IPIV,TAUQQ,TAUQA,TAUQB,W1,W2,KTAB,WKTAB)
C
99001 FORMAT (//,1X,79('='),/,30X,'<SYMCHECKIJ_STD>   PART',I3,/,20X,A,
     &        /,1X,79('='),/)
99002 FORMAT (I3,' Q',2I3,2X,2I3,3X,A,2E13.5,/,22X,A,2E13.5,2F20.12,/)
99003 FORMAT (3X,'ISYM ',I3,' unitary ',L1,3X,A,1X,3F5.1,:,2X,A,1X,
     &        3F5.1,:,2X,A,1X,3F5.1)
99004 FORMAT (/,10X,'run A:  run <TAU_STD_KPOINTS> ',
     &        'making use of symmetry')
99005 FORMAT (/,10X,'run B:  run <TAU_STD_KPOINTS> ',
     &        'WITHOUT making use of symmetry  ---  NSYM = 0')
99006 FORMAT (3X,'ISYM ',I3,' unitary ',L1,3X,A,1X,3F5.1,:,4X,'T:',I4,
     &        3X,A)
99007 FORMAT (/,5X,'action of the symmetry operations  R  on:',/,7X,
     &        '- magnetisation vector ->M',/,7X,
     &        '- vector ->R in real space',/,7X,
     &        '- vector ->k in reciprocal space',/,5X,
     &        'antiunitary rotations of ->k include inversion I',/)
99008 FORMAT (/,5X,'action of the symmetry operations  R  on:',/,7X,
     &        '(T): scattering path operator TAU(->k)   ',5X,
     &        'TAU(Rk) = R TAU(k) R+',/,7X,
     &        '(m): single site m-matrix m_q            ',5X,
     &        '      m = R m R+',/,7X,
     &        '(G): KKR structure constant matrix G(->k)',5X,
     &        '  G(Rk) = R G(k) R+',/,5X,'for arbitrary vector ->k = (',
     &        2(F5.2,','),F5.2,' )',/)
99009 FORMAT (//,14x,53('='),/,14x,I2,
     &        ' check integration of TAU - matrix    for  IQ =',I3,/,
     &        14x,53('='),/)
99010 FORMAT (/,10X,'<SYMCHECKIJ_STD>  finished PART',I2,:,5X,
     &        'ERROR code ',3I5)
99011 FORMAT (//,14X,53('='),/,34X,'SUMMARY',/,14X,53('='),//,10X,
     &        'run parameters:',//,10X,'NQ               ',I5,/10X,
     &        'KMROT            ',I5,/,10X,'NSYM             ',I5,/,10X,
     &        'NSYMACCEPTED     ',I5,/,10X,'NWEDGE           ',I5)
99012 FORMAT (10X,'NELMT            ',I5)
99013 FORMAT (/,10X,'ERROR codes:',//,10X,
     &        'PART 1  TAU(Rk) = R TAU(k) R+                    :',3I5)
99014 FORMAT (10X,
     &     'PART 4  run <TAU_STD_KPOINTS> making use of symmetry      :'
     &     ,I5,/,10X,
     &     'PART 5  TAU(q,q) = R TAU(q,q) R+   for q'' = R q  :',I5,/)
99015 FORMAT (//,1X,79('*'),/,35X,'<SYMCHECKIJ_STD>',/,1X,79('*'),//,
     &        21X,'BZ-integral   Int d^3k  TAU(k)[L1,L2] ')
99016 FORMAT (/,10X,'check the relation     TAU(q,q) = R TAU(q,q) R+',
     &        '   for q'' = R q',/)
      END
