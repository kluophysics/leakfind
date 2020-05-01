C*==posaniprepare.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE POSANIPREPARE(ERYD,P,IWRREGWF,IWRIRRWF,IPRINT,MSSQ,
     &                         MSST,CALCINT,GETIRRSOL,TAUQ,TAUT,MEZZ,
     &                         MEZJ)
C    *******************************************************************
C    *                                                                 *
C    *  determine and WRITE out the information on the                 *
C    *  positron state at the bottom of the positron band              *
C    *                                                                 *
C    *******************************************************************
      USE MOD_KSPACE,ONLY:NKTAB
      USE MOD_SYMMETRY,ONLY:NSYM,MROTK,SYMACCEPTED,SYMUNITARY,NSYMMAX
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS
      USE MOD_CALCMODE,ONLY:IREL,KMROT
      USE MOD_ENERGY,ONLY:ETAB,NETAB
      USE MOD_LATTICE,ONLY:BBAS
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,IND0Q,NKMQ,NKKR,NCPLWF
      USE MOD_SITES,ONLY:NQ,NQMAX,DROTQ,IQAT
      USE MOD_TYPES,ONLY:NT,NTMAX,IMT,NCPLWFMAX,IKMCPLWF
      USE MOD_FILES,ONLY:DATSET,LSYSTEM,SYSTEM,IFILCBWF,LDATSET,IOTMP
      USE MOD_CONSTANTS,ONLY:C0,C1
      IMPLICIT NONE
C*--POSANIPREPARE23
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='POSANIPREPARE')
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      COMPLEX*16 ERYD,P
      INTEGER IPRINT,IWRIRRWF,IWRREGWF
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 DMAMC(:,:),DMATT(:,:,:),DROT(NKMMAX,NKMMAX,NSYMMAX),
     &           DTILT(:,:,:),JF(:,:,:),JG(:,:,:),MSSQL(:,:),TAUA(:,:,:)
     &           ,TAUAB(:,:,:,:),TAUQL(:,:),TAUQQ(:,:,:,:),TAURAT,
     &           ZF(:,:,:),ZF_ST(:,:),ZG(:,:,:),ZG_ST(:,:)
      CHARACTER*80 FILNAM
      INTEGER I,IA_ERR,IK,IM,IQ,IQA,IQB,IR,IRTOP,IT,ITA,ITB,IX,J,
     &        LFILNAM,M,N,NKPTS0,NKTABMAX,NSYMACCEPTEDTMP,NSYMTMP
      REAL*8 KTAB(:,:),WKTAB(:)
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE TAUA,DMAMC,TAUAB,DMATT,DTILT,MSSQL,TAUQL
      ALLOCATABLE KTAB,WKTAB,TAUQQ,JG,JF,ZG,ZF,ZF_ST,ZG_ST
C
      ALLOCATE (ZG_ST(NRMAX,NTMAX),ZF_ST(NRMAX,NTMAX))
      ALLOCATE (JF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (JG(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZF(NRMAX,NCPLWFMAX,NKMMAX))
      ALLOCATE (ZG(NRMAX,NCPLWFMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      ALLOCATE (MSSQL(NKMMAX,NKMMAX),TAUQL(NKMMAX,NKMMAX))
      ALLOCATE (DMATT(NKMMAX,NKMMAX,NTMAX),DTILT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TAUA(NKMMAX,NKMMAX,NTMAX),DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (TAUQQ(NKMMAX,NKMMAX,NQMAX,NQMAX))
      ALLOCATE (TAUAB(NKMMAX,NKMMAX,NTMAX,NTMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: TAUAB')
C
      IF ( IREL.LT.3 ) THEN
         WRITE (6,*) 'SORRY:  <POSANIPREPARE> does not work'
         WRITE (6,*) 'in the scalar relativistic mode'
         WRITE (6,*) 'use fully relativistic mode instead'
         STOP
      END IF
C
      IF ( KMROT.NE.0 ) THEN
         WRITE (6,*) ' #########################################'
         WRITE (6,*) ' <POSANIPREPAREARE> does not yet run for '
         WRITE (6,*) ' rotated magnetic moments    KMROT = ',KMROT
         STOP
      END IF
C
      IF ( IWRREGWF.NE.1 .OR. IWRIRRWF.NE.1 .OR. .NOT.CALCINT .OR. 
     &     .NOT.GETIRRSOL ) THEN
         WRITE (6,*) 'STOP in <POSANIPREPARE>'
         WRITE (6,*) 'positron data not available'
         WRITE (6,*) 'IWRREGWF,IWRIRRWF,CALCINT,GETIRRSOL: ',IWRREGWF,
     &               IWRIRRWF,CALCINT,GETIRRSOL
         STOP
      END IF
C
      IF ( NETAB(1).NE.1 ) STOP '<POSANIPREPARE>: check NETAB'
      IF ( ABS(ERYD-ETAB(1,1)).GT.1D-8 )
     &      STOP '<POSANIPREPARE>: check   ERYD  and P'
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C      read positronic wave functions -- keep only s_1/2 partial wave
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      DO IT = 1,NT
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         CALL WAVFUN_READ_REL(IFILCBWF,IT,0,ZG,ZF,JG,JF,IRTOP,NCPLWF,
     &                        IKMCPLWF)
C
         ZG_ST(:,IT) = ZG(:,1,1)
         ZF_ST(:,IT) = ZF(:,1,1)
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C=======================================================================
C                   calculate ALL elements of TAU-matrix
C=======================================================================
C
C----------------- ignore symmetry and use about the same k-mesh density
C--------------------------------- as used for calculation of TAU(IQ,IQ)
C
      NKPTS0 = NKTAB*NSYM
      NKTABMAX = NINT(NKPTS0*1.05)
C
      NSYMTMP = 1
      NSYMACCEPTEDTMP = 1
C
      CALL KMESHS(IOTMP,NSYMTMP,SYMACCEPTED,MROTK,IPRINT,NKPTS0,BBAS,
     &            DATSET,LDATSET)
C
      REWIND (IOTMP)
      READ (IOTMP) NKTAB
      NKTABMAX = NKTAB
      ALLOCATE (KTAB(3,NKTABMAX),WKTAB(NKTABMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:POSANIPREPARE -> KTAB'
      READ (IOTMP) ((KTAB(IX,IK),IX=1,3),WKTAB(IK),IK=1,NKTAB)
      CLOSE (IOTMP)
C
      LFILNAM = LDATSET + 12
      FILNAM = DATSET(1:LDATSET)//'_pos_dat.pan'
C
      CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:LFILNAM))
C
      WRITE (IOTMP,99005) SYSTEM(1:LSYSTEM),ERYD
C
C---------------------------------------------------- calculate TAU[i,j]
C
      CALL KLOOPSIJ(ERYD,P,NQ,NKKR,NKMQ,IND0Q,TAUQQ,DROT,SYMUNITARY,
     &              SYMACCEPTED,MSSQ,WKTAB,KTAB,NKTAB,NSYMTMP,
     &              NSYMACCEPTEDTMP,NQMAX,NKMMAX,NKTABMAX)
C
C-------------------------------------------------- CHECK BZ-integration
C
      WRITE (6,99001)
      DO IQ = 1,NQ
C
         N = NKMQ(IQ)
         DO J = 1,N
            DO I = 1,N
               IF ( ABS(TAUQ(I,J,IQ)).GT.1D-5 ) THEN
                  TAURAT = TAUQQ(I,J,IQ,IQ)/TAUQ(I,J,IQ)
                  IF ( ABS(1D0-TAURAT).GT.0.01D0 ) WRITE (6,99002) IQ,I,
     &                 J,TAUQ(I,J,IQ),TAUQQ(I,J,IQ,IQ)
               END IF
            END DO
         END DO
      END DO
C
C=======================================================================
C----------------------------------- set up the matrices TAUQL and MSSQL
C----------------------------------- that are referred to the local axis
C----- the resulting projection matrices DMATT and DTILT, as well as the
C------------------ matrices TAUA and TAUAB also refer to the local axis
C-------------------------  TAUQ(IQ) for equivalent sites are the same !
C
      M = NKMMAX
      N = NKMQ(1)
C
      DO ITA = 1,NT
         IQ = IQAT(1,ITA)
C
         IF ( NKMQ(IQ).NE.N ) THEN
            WRITE (6,*) '<POSANIPREPARE> assumes NKMQ the same',
     &                  'for all sites  IQ'
            WRITE (6,*) 'NKMQ(1) = N: ',N
            WRITE (6,*) 'NKMQ(IQ):    ',NKMQ(IQ)
            WRITE (6,*) 'IQ  IT:      ',IT,IQ
            STOP
         END IF
C
         IF ( KMROT.NE.0 ) THEN
C
            CALL ROTATE(MSSQ(1,1,IQ),'G->L',MSSQL,N,DROTQ(1,1,IQ),M)
C
            CALL ROTATE(TAUQ(1,1,IQ),'G->L',TAUQL,N,DROTQ(1,1,IQ),M)
C
         ELSE
C
            DO J = 1,N
               CALL ZCOPY(N,MSSQ(1,J,IQ),1,MSSQL(1,J),1)
               CALL ZCOPY(N,TAUQ(1,J,IQ),1,TAUQL(1,J),1)
            END DO
C
         END IF
C
         CALL GETDMAT(TAUQL,DMATT(1,1,ITA),DTILT(1,1,ITA),DMAMC,N,MSSQL,
     &                MSST(1,1,ITA),M)
C
         CALL ZGEMM('N','N',N,N,N,C1,DMATT(1,1,ITA),M,TAUQL,M,C0,
     &              TAUA(1,1,ITA),M)
C
      END DO
C
      DO ITA = 1,NT
         DO ITB = 1,NT
C
            CALL ZGEMM('N','N',N,N,N,C1,TAUA(1,1,ITA),M,DTILT(1,1,ITB),
     &                 M,C0,TAUAB(1,1,ITA,ITB),M)
C
         END DO
      END DO
C
C----------------------------------------------- CHECK TAUT against TAUA
C
      WRITE (6,99003)
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         N = NKMQ(IQ)
         DO J = 1,N
            DO I = 1,N
               IF ( ABS(TAUT(I,J,IT)).GT.1D-5 ) THEN
                  TAURAT = TAUA(I,J,IT)/TAUT(I,J,IT)
                  IF ( ABS(1D0-TAURAT).GT.0.01D0 ) WRITE (6,99004) IT,I,
     &                 J,TAUT(I,J,IT),TAUA(I,J,IT)
               END IF
            END DO
         END DO
      END DO
C
C=======================================================================
C                        write positron data
C=======================================================================
C
      DO IT = 1,NT
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         WRITE (IOTMP,99006) IT,IQ,IRTOP
         WRITE (IOTMP,99007) 'TAUT              ',(TAUA(I,I,IT),I=1,2)
         WRITE (IOTMP,99007) 'DMATT             ',(DMATT(I,I,IT),I=1,2)
         WRITE (IOTMP,99007) 'DTILT             ',(DTILT(I,I,IT),I=1,2)
         WRITE (IOTMP,99007) 'overlap matrix MZZ',(MEZZ(I,I,IT,1),I=1,2)
         WRITE (IOTMP,99007) 'overlap matrix MZJ',(MEZJ(I,I,IT,1),I=1,2)
         WRITE (IOTMP,99007) 's_1/2-like radial wave function g(r) f(r)'
         DO IR = 1,IRTOP
            WRITE (IOTMP,99009) R(IR,IM),ZG_ST(IR,IT),ZF_ST(IR,IT)
         END DO
C
      END DO
C
      WRITE (IOTMP,99007) 'TAUAB'
      DO ITA = 1,NT
         DO ITB = 1,NT
            WRITE (IOTMP,99008) ITA,ITB,(TAUAB(I,I,ITA,ITB),I=1,2)
         END DO
      END DO
C
      WRITE (IOTMP,99007) 'TAUQQ'
      DO IQA = 1,NT
         DO IQB = 1,NT
            WRITE (IOTMP,99008) IQA,IQB,(TAUQQ(I,I,IQA,IQB),I=1,2)
         END DO
      END DO
C
      WRITE (6,99010) FILNAM(1:LFILNAM)
C
C=======================================================================
C
      DEALLOCATE (ZG,ZF,TAUA,DMAMC,TAUAB,DMATT,DTILT,MSSQL,TAUQL)
      DEALLOCATE (KTAB,WKTAB,TAUQQ)
C
      STOP
99001 FORMAT (/,' comparing  TAUQ  and  TAUQQ ')
99002 FORMAT (' IQ=',I3,3X,'(I,J)=',2I3,3X,'TAUQ:',2(3X,2E12.5))
99003 FORMAT (/,' comparing  TAUT  and  TAUA ')
99004 FORMAT (' IT=',I3,3X,'(I,J)=',2I3,3X,'TAUT:',2(3X,2E12.5))
C
99005 FORMAT ('positron data for system:',A,/,'energy  E = ',2F10.6,
     &        ' Ry',/,'all variables restricted to  l=0')
99006 FORMAT (79('*'),/,'IT=',I3,3X,'IQ=',I3,3X,'IRTOP=',I3)
99007 FORMAT (A,:,(/,2(2E14.6,3x)))
99008 FORMAT (2I3,2(2E14.6,3x))
99009 FORMAT (5E14.6)
99010 FORMAT (/,1X,79('*'),/,33X,'<POSANIPREPARE>',/,1X,79('*'),//,10X,
     &        'positron data written to file    ',A,/)
      END
