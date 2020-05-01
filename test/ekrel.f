C*==ekrel.f    processed by SPAG 6.70Rc at 17:10 on 19 Apr 2017
      SUBROUTINE EKREL
C   ********************************************************************
C   *                                                                  *
C   *             determine the  E(k) - relation                       *
C   *                                                                  *
C   *   NKDIR   number of k-path segments created                      *
C   *   LBLKDIR name of each segment - contains name of first, last    *
C   *           and intermediate k-point (symmetry point in BZ)        *
C   *   INDKDIR index-table giving the index of the last k-point       *
C   *           of a segment in a continous table                      *
C   *                                                                  *
C   *                                                                  *
C   *   E-values occuring for more k-points as  NKTAB*0.9              *
C   *   are seen as spourious and are filtered out                     *
C   *                                                                  *
C   ********************************************************************
      USE MOD_CONSTANTS,ONLY:C1
      USE MOD_FILES,ONLY:IPRINT,IFILBUILDBOT,WRBUILDBOT
      USE MOD_CALCMODE,ONLY:IREL,ORBPOL,DMFT
      USE MOD_TYPES,ONLY:NT,NTMAX
      USE MOD_LATTICE,ONLY:BBAS,BRAVAIS
      USE MOD_ANGMOM,ONLY:NKMMAX,IND0Q,NKMQ,NLMQ,NLM,NKKR,MEZJ,MEZZ,
     &    TSSQ,MSSQ,MSST,SSST,TSST
      USE MOD_SITES,ONLY:NQMAX,NQ
      USE MOD_KSPACE,ONLY:KTAB,NKTABMAX
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      USE MOD_CPA,ONLY:NCPA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='EKREL')
      INTEGER NEKMAX,NKDIRMAX
      PARAMETER (NEKMAX=4000,NKDIRMAX=20)
      LOGICAL INTERPOLATE_TSS
      PARAMETER (INTERPOLATE_TSS=.TRUE.)
C
C Local variables
C
      LOGICAL CALCINT
      COMPLEX*16 CWORK(:),ERYD,MAUX(:,:),MQS(:,:,:,:),MSS_AUX(:,:,:,:),
     &           P,TAUK(:,:),TSS_AUX(:,:,:,:),VL(1),VR(1),WORK(:),
     &           ZEIG(:)
      REAL*8 DEL(3),EMAX,EMIN,ETOL,E_AUX(:),KA(3,NKDIRMAX),
     &       KE(3,NKDIRMAX),KP(:),RAT,RSUM,RWORK(:),SUMD
      REAL*8 DNRM2
      REAL E0,EJK(:,:,:),EJK0(:),EMIN4,ETOL4
      INTEGER I,I0,I1,IA_ERR,ID,IE,IE1,IE2,IFIL,IFLAG,IK,IK0,IKD,
     &        INDKDIR(NKDIRMAX),IQ,IS,IS0,J,J1,KPATH,N,NBK(:,:),NBK0,
     &        NBKMAX,NEK,NE_AUX,NKDIR,NKTAB,NKTAB0,NKTABD(NKDIRMAX),
     &        NOCCUR,NPEV(:,:,:),NPOSEV,NSPIN,NTHRESH
      CHARACTER*8 LBLKDIR(NKDIRMAX)
      CHARACTER*4 STR4
C
C*** End of declarations rewritten by SPAG
C
      DATA CALCINT/.FALSE./
C
      ALLOCATABLE ZEIG,TAUK,NPEV,MAUX,WORK,RWORK,KP,EJK,NBK
      ALLOCATABLE MQS,CWORK,EJK0,E_AUX, TSS_AUX,MSS_AUX
C
      IF ( NCPA.NE.0 ) CALL STOP_ERROR(ROUTINE,
     &      'no E(k) sensible for alloys >> use Bloch spectral function'
     &      )
C
      IF ( IREL.EQ.2 ) THEN
         ALLOCATE (MQS(NLM,NLM,NQ,2),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:ekrelspsrel -> MQS'
         NSPIN = 2
      ELSE
         NSPIN = 1
      END IF
C
      ALLOCATE (EJK(2*NKKR,NKTABMAX,NSPIN))
      ALLOCATE (TAUK(NKKR,NKKR),ZEIG(NKKR),NBK(NKTABMAX,NSPIN))
      ALLOCATE (WORK(50*NKKR),RWORK(2*NKKR),KP(NKTABMAX))
      ALLOCATE (NPEV(NKTABMAX,2,NSPIN),MAUX(NKKR,NKKR),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:ekrel -> NBK'
      WRITE (6,99001)
C
      KPATH = 1
      NKTAB = 1000
      EMIN = 0.01D0
      EMAX = 1.0D0
      NEK = 501
      NKDIR = 0
      DO ID = 1,NKDIRMAX
         LBLKDIR(ID) = '        '
      END DO
C
      IFIL = 10
C
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      CALL SECTION_SET_REAL('EMIN',EMIN,9999D0,0)
      CALL SECTION_SET_REAL('EMAX',EMAX,9999D0,0)
C
      CALL SECTION_SET_INTEGER('NE',NEK,9999,0)
      CALL SECTION_SET_INTEGER('NK',NKTAB,9999,0)
      CALL SECTION_SET_INTEGER('KPATH',KPATH,9999,0)
      CALL SECTION_SET_INTEGER('NKDIR',NKDIR,9999,0)
      NKDIR = MIN(NKDIR,9)
C
      IF ( NKDIR.NE.0 ) THEN
         DO ID = 1,NKDIR
            STR4 = 'KA'//CHAR(ICHAR('1')-1+ID)
            CALL SECTION_SET_REAL_ARRAY(STR4,KA(1,ID),3,3,1,9999D0,0)
            STR4 = 'KE'//CHAR(ICHAR('1')-1+ID)
            CALL SECTION_SET_REAL_ARRAY(STR4,KE(1,ID),3,3,1,9999D0,0)
            WRITE (6,99002) ID,(KA(I,ID),I=1,3),(KE(I,ID),I=1,3)
            LBLKDIR(ID) = '        '
         END DO
C
      ELSE
C
         CALL KDIRTAB(BRAVAIS,KPATH,NKDIR,LBLKDIR,KA,KE,BBAS)
C
      END IF
C
      IF ( NKTAB.GT.NKTABMAX ) THEN
         NKTAB = NKTABMAX
         WRITE (6,*) 'WARNING from <EKREL>:  NK reduced to array size',
     &               NKTABMAX
      END IF
      IF ( NEK.GT.NEKMAX ) THEN
         NEK = NEKMAX
         WRITE (6,*) 'WARNING from <EKREL>:  NEK reduced to array size',
     &               NEKMAX
      END IF
C
      DO IS = 1,NSPIN
         DO IK = 1,NKTAB
            NBK(IK,IS) = 0
         END DO
      END DO
C
      ETOL = (EMAX-EMIN)/DBLE(NEK-1)
C
C----------------------------------------------------------- avoid E = 0
      IF ( EMIN.LE.0D0 ) THEN
         IF ( ABS((EMIN/ETOL)-NINT(EMIN/ETOL)).LE.1D-10 ) THEN
            NEK = MAX(2,NEK-1)
            ETOL = (EMAX-EMIN)/DBLE(NEK-1)
         END IF
      END IF
C
      ETOL4 = REAL(ETOL)
      EMIN4 = REAL(EMIN)
C
C-----------------------------------------------------------------------
      IF ( NKDIR.EQ.0 ) THEN
         NKDIR = 1
         DO I = 1,3
            KA(I,1) = 0.0D0
            KE(I,1) = 0.0D0
         END DO
         KE(1,1) = 1.0D0
         NKTABD(1) = NKTAB
         LBLKDIR(1) = 'GG-x -1 '
      ELSE
         SUMD = 0.0D0
         DO ID = 1,NKDIR
            DO I = 1,3
               DEL(I) = KE(I,ID) - KA(I,ID)
            END DO
            RSUM = DNRM2(3,DEL,1)
            NKTABD(ID) = INT(1000000*RSUM)
            SUMD = SUMD + RSUM
         END DO
         NKTAB0 = 0
         DO ID = 1,NKDIR
            NKTABD(ID) = INT(NKTAB*NKTABD(ID)/(SUMD*1000000))
            NKTABD(ID) = MAX(2,NKTABD(ID))
            NKTAB0 = NKTAB0 + NKTABD(ID)
         END DO
         IF ( NKTAB0.GT.NKTABMAX ) THEN
            RAT = 0.95D0*DBLE(NKTABMAX)/DBLE(NKTAB0)
            NKTAB0 = 0
            DO ID = 1,NKDIR
               NKTABD(ID) = MAX(2,INT(RAT*NKTABD(ID)))
               NKTAB0 = NKTAB0 + NKTABD(ID)
            END DO
         END IF
C
         NKTAB = NKTAB0
      END IF
C
      INDKDIR(1) = NKTABD(1)
      DO ID = 2,NKDIR
         INDKDIR(ID) = INDKDIR(ID-1) + NKTABD(ID)
      END DO
C
      WRITE (6,99003) NKTAB
C
C-----------------------------------------------------------------------
C----------------------------------- set up k-points and linear array KP
C
      DEALLOCATE (KTAB)
      ALLOCATE (KTAB(3,NKTAB))
C
      IK = 0
      KP(1) = 0D0
      DO ID = 1,NKDIR
C
         DO I = 1,3
            DEL(I) = (KE(I,ID)-KA(I,ID))/DBLE(NKTABD(ID)-1)
         END DO
C
         DO IKD = 1,NKTABD(ID)
            IK = IK + 1
            DO I = 1,3
               KTAB(I,IK) = KA(I,ID) + DEL(I)*DBLE(IKD-1)
            END DO
            IF ( ABS(KTAB(1,IK))+ABS(KTAB(2,IK))+ABS(KTAB(3,IK))
     &           .LT.1D-8 ) KTAB(3,IK) = 1D-8
            IF ( IKD.NE.1 ) THEN
               KP(IK) = KP(IK-1) + DNRM2(3,DEL,1)
            ELSE IF ( IK.EQ.1 ) THEN
               KP(IK) = 0.0D0
            ELSE
               KP(IK) = KP(IK-1)
            END IF
         END DO
      END DO
C
      IF ( IPRINT.GE.1 ) THEN
         IK = 0
         DO ID = 1,NKDIR
            DO IKD = 1,NKTABD(ID)
               IK = IK + 1
               WRITE (6,'(A,4I5,3F8.3,F10.3)') ' ->K ',ID,NKTABD(ID),
     &                IKD,IK,(KTAB(I,IK),I=1,3),KP(IK)
            END DO
         END DO
      END IF
C
      IF ( DMFT ) THEN
         ALLOCATE (CWORK(NEK))
         IF ( ALLOCATED(DMFTSIGMA) ) THEN
            DEALLOCATE (DMFTSIGMA)
            ALLOCATE (DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEK))
         END IF
         DO IE = 1,NEK
            ERYD = EMIN + DBLE(IE-1)*ETOL
            CWORK(IE) = DCMPLX(DREAL(ERYD),0.0001D0)
         END DO
         CALL DMFT_READSIG(NEK,CWORK,IPRINT)
      END IF
C
C=======================================================================
      IF ( INTERPOLATE_TSS ) THEN
C
         NE_AUX = 100
         ALLOCATE (E_AUX(NE_AUX),TSS_AUX(NKMMAX,NKMMAX,NT,NE_AUX))
         ALLOCATE (MSS_AUX(NKMMAX,NKMMAX,NT,NE_AUX))
C
         DO IE = 1,NE_AUX
C
            ERYD = EMIN + DBLE(IE-1)*(EMAX-EMIN)/DBLE(NE_AUX-1)
C
            IF ( ABS(ERYD).LT.1D-6 ) ERYD = 1D-6
C
            E_AUX(IE) = DREAL(ERYD)
            IF ( IPRINT.GE.3 ) WRITE (6,*) 'E_AUX :',IE,E_AUX
C
            IF ( DMFT ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &           = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE))
     &           ,0.0D0)
C
            CALL RUNSSITE(.FALSE.,0,0,87,.FALSE.,ERYD,P,IPRINT,
     &                    TSS_AUX(1,1,1,IE),MSS_AUX(1,1,1,IE),SSST,MEZZ,
     &                    MEZJ,ORBPOL)
C
         END DO
C
      END IF
C=======================================================================
C
      WRITE (6,99004) EMIN,EMAX,ETOL,NEK
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      DO IE = 1,NEK
C
         IE2 = MOD(IE,2) + 1
         ERYD = EMIN + DBLE(IE-1)*ETOL
C
         IF ( IPRINT.GT.0 ) WRITE (6,'(A,I5,2F10.5)') 'E=',IE,
     &                             DREAL(ERYD)
C
         IF ( DMFT ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &        = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)),
     &        0.0D0)
C
         CALL RUNSSITE(CALCINT,0,0,87,.FALSE.,ERYD,P,IPRINT,TSST,MSST,
     &                 SSST,MEZZ,MEZJ,ORBPOL)
C
         CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
         IF ( IREL.EQ.2 ) THEN
            DO IQ = 1,NQ
               DO J = 1,NLMQ(IQ)
                  CALL ZCOPY(NLM,MSSQ(1,J,IQ),1,MQS(1,J,IQ,1),1)
                  CALL ZCOPY(NLM,MSSQ(NLMQ(IQ)+1,NLMQ(IQ)+J,IQ),1,
     &                       MQS(1,J,IQ,2),1)
               END DO
            END DO
         END IF
C
         CALL STRCC(ERYD,.FALSE.)
C
C ISISISIISIISISISISISISISIISISISISISISISISISISISISISIISISISISISISISISIS
         DO IS = 1,NSPIN
C
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
            DO IK = 1,NKTAB
C
               CALL STRSET(IK,KTAB(1,IK),MAUX,TAUK,P)
C
               IF ( IREL.NE.2 ) THEN
C
                  CALL SETKKR(NQ,NKMQ,IND0Q,TAUK,MSSQ,NQMAX,NKKR,NKMMAX)
C
               ELSE
C
                  DO IQ = 1,NQ
                     I1 = IND0Q(IQ) + 1
                     N = NLMQ(IQ)
                     DO J = 1,N
                        J1 = IND0Q(IQ) + J
                        CALL ZAXPY(N,C1,MQS(1,J,IQ,IS),1,TAUK(I1,J1),1)
                     END DO
                  END DO
C
               END IF
C
               CALL ZGEEV('N','N',NKKR,TAUK,NKKR,ZEIG,VL,1,VR,1,WORK,
     &                    50*NKKR,RWORK,IFLAG)
C
               NPOSEV = 0
               DO I = 1,NKKR
                  IF ( DREAL(ZEIG(I)).GT.0.0D0 ) NPOSEV = NPOSEV + 1
               END DO
C
               NPEV(IK,IE2,IS) = NPOSEV
C
            END DO
C KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
            IF ( IE.EQ.1 ) THEN
               NBK(IK,IS) = 0
            ELSE
               IE1 = 3 - IE2
               DO IK = 1,NKTAB
                  DO I = 1,(NPEV(IK,IE2,IS)-NPEV(IK,IE1,IS))
                     IF ( NBK(IK,IS).LE.2*NKKR ) THEN
                        NBK(IK,IS) = NBK(IK,IS) + 1
                        EJK(NBK(IK,IS),IK,IS) = EMIN4 + (REAL(IE-1)-0.5)
     &                     *ETOL4
                     ELSE
                        WRITE (6,99012) 2*NKKR,IK
                     END IF
                  END DO
               END DO
            END IF
C
         END DO
C
C ISISISIISIISISISISISISISIISISISISISISISISISISISISISIISISISISISISISISIS
C
      END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      WRITE (IFIL,99008)
      WRITE (IFIL,99005) 'NKDIR     ',NKDIR
      DO I = 1,NKDIR
         WRITE (IFIL,99006) LBLKDIR(I)
      END DO
      WRITE (IFIL,99008)
      DO I = 1,NKDIR
         WRITE (IFIL,99005) 'INDKDIR   ',INDKDIR(I)
      END DO
      WRITE (IFIL,99008)
      WRITE (IFIL,99005) 'NKTAB      ',NKTAB
      DO IK = 1,NKTAB
         WRITE (IFIL,99007) KP(IK)
      END DO
      WRITE (IFIL,99008)
      WRITE (IFIL,99005) 'NBAND     ',NKKR
C
      NBKMAX = 0
      DO IS = 1,NSPIN
         DO IK = 1,NKTAB
            NBKMAX = MAX(NBKMAX,NBK(IK,IS))
         END DO
      END DO
C
      WRITE (6,99009) NBKMAX
C
C************************************************************************
C              filter out spurious eigen values
C************************************************************************
C
      ALLOCATE (EJK0(NKTABMAX))
      ETOL = 1D-8
      NTHRESH = NINT(NKTAB*0.9)
C
      DO IS0 = 1,NSPIN
         IS = IS0
C
         IK0 = 1
         NBK0 = NBK(IK0,IS0)
         EJK0(1:NBK0) = EJK(1:NBK0,IK0,IS0)
C
         DO I0 = 1,NBK0
C
            E0 = EJK0(I0)
C
C=======================================================================
C          check how often  E0  occurs all k-vectors
C=======================================================================
C
            NOCCUR = 0
            DO IK = 1,NKTAB
C
               DO I = 1,NBK(IK,IS)
                  IF ( ABS(E0-EJK(I,IK,IS)).LT.ETOL ) THEN
                     NOCCUR = NOCCUR + 1
                     EXIT
                  END IF
               END DO
C
            END DO
C
C=======================================================================
C          remove spurious E-value E0 from tables
C=======================================================================
C
            IF ( NOCCUR.GT.NTHRESH ) THEN
C
               IF ( IPRINT.GE.-5 ) WRITE (6,*)
     &               'removing spurious E-value ',E0
               DO IK = 1,NKTAB
C
                  DO I = 1,NBK(IK,IS)
                     IF ( ABS(E0-EJK(I,IK,IS)).LT.ETOL ) THEN
                        DO J = I,NBK(IK,IS) - 1
                           EJK(J,IK,IS) = EJK(J+1,IK,IS)
                        END DO
                        NBK(IK,IS) = NBK(IK,IS) - 1
                        EXIT
                     END IF
                  END DO
C
               END DO
C
            END IF
C
C=======================================================================
C
         END DO
      END DO
C************************************************************************
C
      DO IS = 1,NSPIN
         DO IK = 1,NKTAB
C
            WRITE (IFIL,'(2I5)') IK,NBK(IK,IS)
            WRITE (IFIL,'(10F8.4)') (EJK(I,IK,IS),I=1,NBK(IK,IS))
C
            IF ( IPRINT.GT.0 ) THEN
               WRITE (6,99010) IK,(KTAB(I,IK),I=1,3),KP(IK)
               WRITE (6,99011) NBK(IK,IS),EMIN,EMAX,
     &                         (EJK(I,IK,IS),I=1,NBK(IK,IS))
            END IF
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
            IF ( IK.LT.3 .AND. WRBUILDBOT ) WRITE (IFILBUILDBOT,99013)
     &           ROUTINE(1:LEN_TRIM(ROUTINE)),IS,IK,
     &           (EJK(I,IK,IS),I=1,NBK(IK,IS))
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BUILDBOT
C
         END DO
      END DO
C
      DEALLOCATE (ZEIG,TAUK,NPEV,MAUX,WORK,RWORK,KP,EJK,NBK,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:ekrel -> NBK'
      WRITE (6,*) '          <EKREL> - finished'
      STOP
C
99001 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      ******  * *      *        *****   ******  *           *'
     &  ,/,10X,
     &  '*      *      *  *       *       *    *  *       *           *'
     &  ,/,10X,
     &  '*      *      *  *   *   *       *    *  *       *           *'
     &  ,/,10X,
     &  '*      ****   *  *  *    *  ***  *****   ****    *           *'
     &  ,/,10X,
     &  '*      *      *  * **    *       *  *    *       *           *'
     &  ,/,10X,
     &  '*      *      *  **  *   *       *   *   *       *           *'
     &  ,/,10X,
     &  '*      ******  * *    * *        *    *  ******  ******      *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99002 FORMAT (' DIR ',I3,' KA=',3F7.3,'  KE=',3F7.3)
99003 FORMAT (5X,'number of k-points created ',I6)
99004 FORMAT (/,5X,'calculation of energy-bands between ',F6.3,' and ',
     &        F6.3,' Ry',/,5X,'tolerance: ',F7.4,
     &        ' Ry    corresponding to',I5,'  E-points      ',I7,/)
99005 FORMAT (A10,I10)
99006 FORMAT (10X,A)
99007 FORMAT (8F10.4)
99008 FORMAT (80('#'))
99009 FORMAT (/,5X,'maximum number of eigenvalues ',I3,/)
99010 FORMAT (' -> K ',I3,3F7.3,F10.3)
99011 FORMAT (I3,' eigenvalues between ',F9.4,'  and ',F9.4,/,(8F10.4))
99012 FORMAT (5X,'number of bands >',I5,'   for IK=',I3,'   bands cut!')
99013 FORMAT ('# BUILDBOT: ',A,':  E(k) relation for  IS =',I3,'  IK =',
     &        I3,/,(1PE22.14))
      END
