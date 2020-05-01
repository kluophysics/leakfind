C*==spec_tau_drive.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_TAU_DRIVE(MEZJ,MEZZ,MSSQ,TSSQ,MSST,SSST,TAUQ,TAUT,
     &                          TSST,PHASK,DMATTG,DTILTG,ERYD,ISTATE,
     &                          IECURR,TSST_GLO,TAUT_GLO,DOS,CALCDOS,
     &                          NEW,P)
C
C
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NKMQ,NLMQ,NLM
      USE MOD_TYPES,ONLY:NTMAX,NT,TXT_T,CONC,ITBOT,ITTOP,LTXT_T,VT,BT,Z
      USE MOD_SITES,ONLY:NQMAX,IQAT,DROTQ,IQBOT,IQTOP,NQ
      USE MOD_ENERGY,ONLY:NEMAX,EFERMI
      USE MOD_RMESH,ONLY:NRMAX
      USE MOD_FILES,ONLY:IFILCBWF,IPRINT,LDATSET,DATSET,RDTAUMQ,
     &    IFILGFWF,FOUND_SECTION
      USE MOD_CALCMODE,ONLY:ORBPOL,DMFT,LDAU,KMROT,IREL
      USE MOD_TRANSPHO_LAYPOT,ONLY:IBLOCHTMP
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_TYPE
      USE MOD_CONSTANTS,ONLY:C0,C1,PI,RY_EV
      USE MOD_MPI,ONLY:MPI,MPI_ID
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C Dummy arguments
C
      LOGICAL CALCDOS
      COMPLEX*16 ERYD,P
      INTEGER IECURR,ISTATE,NEW
      COMPLEX*16 DMATTG(NKMMAX,NKMMAX,NTMAX),DOS(NKMMAX,NKMMAX,NTMAX),
     &           DTILTG(NKMMAX,NKMMAX,NTMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           PHASK(NEMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX),
     &           TAUT_GLO(NKMMAX,NKMMAX,NTMAX),TSSQ(NKMMAX,NKMMAX,NQMAX)
     &           ,TSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST_GLO(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      LOGICAL ATAFINAL,ATAINIT,CALCINT,FEGFINAL,GETIRRSOL
      REAL*8 BT0(:,:),CPACHNG,CPACHTAB(1),VT0(:,:)
      COMPLEX*16 DMAMC(:,:),MAUX(:,:),MQS(:,:,:),W1(:,:)
      CHARACTER*80 FILNAM
      INTEGER I,ICALL,ICPACONV,ICPAFLAG,IECPAFAIL(1),IFIL,INFO,IOTMP,
     &        IPIV(:),IQ,IQBOTSAV,IQTOPSAV,IS,IT,ITBOTSAV,ITCPA,
     &        ITTOPSAV,IWRIRRWF,IWRREGWF,J,LFN,M,N,NCPAFAIL,Z0(:)
      CHARACTER*10 IOUT
      CHARACTER*12 STR12
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
      ALLOCATABLE W1, DMAMC,IPIV,MAUX,MQS,VT0,BT0,Z0
C
      ALLOCATE (MAUX(NKMMAX,NKMMAX),IPIV(NKMMAX))
      IF ( IREL.EQ.2 ) ALLOCATE (MQS(NLM,NLM,2))
C
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      NCPAFAIL = 0
      ICPAFLAG = 0
      ATAFINAL = .FALSE.
      ATAINIT = .FALSE.
      IWRREGWF = 1
      IWRIRRWF = 1
      FEGFINAL = .FALSE.
C
      ALLOCATE (W1(NKMMAX,NKMMAX),DMAMC(NKMMAX,NKMMAX))
C
      TSSQ(1:NKMMAX,1:NKMMAX,1:NQMAX) = C0
C
      ICALL = ICALL + 1
C
C
C=======================================================================
C
C
      IF ( DMFT .OR. LDAU ) THEN
         IF ( ISTATE.EQ.2 ) THEN
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX) = C0
         ELSE IF ( IBLOCHTMP.EQ.0 ) THEN
            WRITE (*,*) 'Im(SIG)==0.0'
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &         = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,
     &         IECURR)),0.0D0)
         ELSE
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &         = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IECURR)
C            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
C     &         = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,
C     &         IECURR)),0.0D0)
         END IF
      END IF
C
C
C
C=======================================================================
      CALL INPUT_FIND_SECTION('SPEC ',0)
C
      IF ( ISTATE.EQ.2 ) THEN
         ATAFINAL = .FALSE.
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('ATAFINAL',ATAFINAL)
            IF ( ATAFINAL ) THEN
               WRITE (6,'(/)')
               WRITE (6,*) '  ATA FOR FINAL STATE  '
               WRITE (6,'(1X,79(''*''),/)')
            END IF
            CALL SECTION_FIND_KEYWORD('FEGFINAL',FEGFINAL)
            IF ( FEGFINAL ) THEN
               WRITE (6,'(/)')
               WRITE (6,*) '  FEG FOR FINAL STATE  '
               WRITE (6,'(1X,79(''*''),/)')
            END IF
C
         END IF
      ELSE
         ATAFINAL = .FALSE.
         ATAINIT = .FALSE.
         IF ( FOUND_SECTION ) THEN
            CALL SECTION_FIND_KEYWORD('ATAINIT',ATAINIT)
            IF ( ATAINIT ) THEN
               WRITE (6,'(/)')
               WRITE (6,*) '  ATA FOR INITIAL STATE  '
               WRITE (6,'(1X,79(''*''),/)')
            END IF
         END IF
C
      END IF
C
C=======================================================================
C
      IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOTSAV = ITBOT
         ITTOPSAV = ITTOP
         IQBOTSAV = IQBOT
         IQTOPSAV = IQTOP
         ITBOT = 1
         ITTOP = NT
         IQBOT = 1
         IQTOP = NQ
      END IF
C
      IF ( ISTATE.EQ.1 ) THEN
         IFIL = IFILCBWF
      ELSE
         IFIL = IFILGFWF
      END IF
C
      IF ( FEGFINAL .AND. ISTATE.EQ.2 ) THEN
         ALLOCATE (BT0(NRMAX,NTMAX),Z0(NTMAX))
         ALLOCATE (VT0(NRMAX,NTMAX))
         VT0 = VT
         BT0 = BT
         Z0 = Z
         VT = 0.0D0
         BT = 0.0D0
         Z = 0
      END IF
C
      CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL,GETIRRSOL,ERYD,P,
     &              IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
      IF ( FEGFINAL .AND. ISTATE.EQ.2 ) THEN
         VT = VT0
         BT = BT0
         Z = Z0
      END IF
C
      IF ( SYSTEM_TYPE(1:2).EQ.'LI' .AND. SUB_SYSTEM(1:6).EQ.'I-ZONE' )
     &     THEN
         ITBOT = ITBOTSAV
         ITTOP = ITTOPSAV
         IQBOT = IQBOTSAV
         IQTOP = IQTOPSAV
      END IF
C
      IF ( NCPA.GT.0 .AND. .NOT.ATAFINAL .OR. CALCDOS .AND. 
     &     .NOT.ATAINIT ) THEN
C
C=======================================================================
C
C
         IF ( RDTAUMQ ) THEN
            WRITE (6,*) 
     &               '<KKRSPEC> TAUQ and MSSQ will be read in from file'
            DO IQ = IQBOT,IQTOP
               CALL READTAU(9,ERYD,IECURR,NEW,W1,0,NT,.FALSE.,
     &                      TAUQ(1,1,IQ),MSSQ(1,1,IQ),IQ,NQ,2,NKMQ(IQ),
     &                      NKMMAX,IPRINT)
            END DO
         END IF
C
         IF ( .NOT.RDTAUMQ ) THEN
            CALL TAU_DRIVE(1,IPRINT,ERYD,P,TSSQ,MSSQ,TSST,MSST,TAUQ,
     &                     ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
C
            IF ( ICPAFLAG.NE.0 ) THEN
               NCPAFAIL = NCPAFAIL + 1
               CPACHTAB(NCPAFAIL) = CPACHNG
               IECPAFAIL(NCPAFAIL) = 1
            END IF
C
            CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
C
            IF ( NCPAFAIL.NE.0 ) THEN
               WRITE (6,99002) CPATOL,NCPAFAIL,
     &                         (IECPAFAIL(I),DREAL(ERYD),CPACHTAB(I),
     &                         I=1,NCPAFAIL)
               WRITE (6,'(1X,79(''*''),/)')
            ELSE IF ( NCPA.NE.0 ) THEN
               WRITE (6,99003)
            END IF
         END IF
      END IF
C
C
C=======================================================================
C     CALCULATE TSSQ FROM MSSQ IN THE GLOBAL COORDINATE SYSTEM
C
      DO IQ = IQBOT,IQTOP
C
         M = NKMMAX
         IF ( IREL.NE.2 ) THEN
            N = NKMQ(IQ)
         ELSE
            N = NKM
         END IF
C
         IF ( IREL.NE.2 ) THEN
C
            W1(1:N,1:N) = MSSQ(1:N,1:N,IQ)
            CALL ZGETRF(N,N,W1,NKMMAX,IPIV,INFO)
            CALL ZGETRI(N,W1,NKMMAX,IPIV,MAUX,NKMMAX*NKMMAX,INFO)
            TSSQ(1:N,1:N,IQ) = W1(1:N,1:N)
C
         ELSE
C
            DO J = 1,NLM
               CALL ZCOPY(NLM,MSSQ(1,J,IQ),1,MQS(1,J,1),1)
               CALL ZCOPY(NLM,MSSQ(NLM+1,NLM+J,IQ),1,MQS(1,J,2),1)
            END DO
C
            DO IS = 1,2
               CALL ZGETRF(NLMQ(IQ),NLMQ(IQ),MQS(1,1,IS),NLM,IPIV,INFO)
               CALL ZGETRI(NLMQ(IQ),MQS(1,1,IS),NLM,IPIV,MAUX,NLM*NLM,
     &                     INFO)
            END DO
C
            DO J = 1,NLMQ(IQ)
               CALL ZCOPY(NLMQ(IQ),MQS(1,J,1),1,TSSQ(1,J,IQ),1)
               CALL ZCOPY(NLMQ(IQ),MQS(1,J,2),1,TSSQ(NLM+1,NLM+J,IQ),1)
            END DO
C
         END IF
         IF ( ATAFINAL .OR. ATAINIT ) TAUQ(1:N,1:N,IQ)
     &        = TSSQ(1:N,1:N,IQ)
C
      END DO
C
C=======================================================================
C     CALCULATE PROJECTION MATRICES IN THE CASE OF CPA
C
C     calculate projection matrices   DMATT  and  DTILT
C     for preselected atom type  IT  on preselected site  IQ
C
C     DM = m(t)-m(c)
C
C     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
C     D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)
C
      IF ( NCPA.GT.0 .OR. CALCDOS ) THEN
C
         DO IT = ITBOT,ITTOP
            IQ = IQAT(1,IT)
            M = NKMMAX
            IF ( IREL.NE.2 ) THEN
               N = NKMQ(IQ)
            ELSE
               N = NKM
            END IF
C
            IF ( CONC(IT).LT.0.995D0 ) THEN
C
C ------------------------- rotate the single site m-matrix if necessary
               IF ( KMROT.NE.0 ) THEN
C
C                  CALL ROTATE(MSST(1,1,IT),'L->G',W1,N,DROTQ(1,1,IQ),M)
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMATTG(1,1,IT),
     &                         DTILTG(1,1,IT),DMAMC,N,MSSQ(1,1,IQ),
     &                         MSST(1,1,IT),M)
C
               ELSE
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMATTG(1,1,IT),
     &                         DTILTG(1,1,IT),DMAMC,N,MSSQ(1,1,IQ),
     &                         MSST(1,1,IT),M)
C
               END IF
            END IF
         END DO
         IF ( IPRINT.GT.1 ) THEN
            DO IT = ITBOT,ITTOP
C
               WRITE (6,99001) 'atom type   ',IT
C
               CALL CMATSTRUCT('DMATTG  ',DMATTG(1,1,IT),NKM,NKMMAX,3,3,
     &                         1,1.0D-9,6)
C
               CALL CMATSTRUCT('DTILTG  ',DTILTG(1,1,IT),NKM,NKMMAX,3,3,
     &                         1,1.0D-9,6)
            END DO
         END IF
C
      END IF
C=======================================================================
C
C
C=======================================================================
C     CALCULATE t and tau-matrix in the global coordinate system
C     in contrast to KKR we need MSST and TSST
C     in locall coordinate system
C
      DO IT = ITBOT,ITTOP
         IQ = IQAT(1,IT)
         IF ( KMROT.NE.0 ) THEN
C
            TSST_GLO(:,:,IT) = TSST(:,:,IT)
C            CALL ROTATE(TSST(1,1,IT),'L->G',TSST_GLO(1,1,IT),N,
C     &                  DROTQ(1,1,IQ),M)
            TAUT_GLO(:,:,IT) = TAUT(:,:,IT)
C            CALL ROTATE(TAUT(1,1,IT),'L->G',TAUT_GLO(1,1,IT),N,
C     &                  DROTQ(1,1,IQ),M)
C
C
            CALL ROTATE(MSST(1,1,IT),'G->L',W1,N,DROTQ(1,1,IQ),M)
            MSST(:,:,IT) = W1
            CALL ROTATE(TSST(1,1,IT),'G->L',W1,N,DROTQ(1,1,IQ),M)
            TSST(:,:,IT) = W1
C
         ELSE
            TSST_GLO(1:NKMMAX,1:NKMMAX,IT) = TSST(1:NKMMAX,1:NKMMAX,IT)
            TAUT_GLO(1:NKMMAX,1:NKMMAX,IT) = TAUT(1:NKMMAX,1:NKMMAX,IT)
         END IF
C
      END DO
C=======================================================================
C    CALCULATE \Lambda resolved DOS needed for XPS regime
C
      IF ( CALCDOS ) THEN
C
         DOS(1:NKMMAX,1:NKMMAX,1:NTMAX) = C0
C
         DO IT = ITBOT,ITTOP
C
            M = NKMMAX
            N = NKMQ(IQAT(1,IT))
C
            CALL ZGEMM('N','N',N,N,N,CPRE,MEZZ(1,1,IT,1),M,TAUT(1,1,IT),
     &                 M,C0,DOS(1,1,IT),M)
            DO J = 1,N
               DOS(J,J,IT) = DOS(J,J,IT) - CPRE*MEZJ(J,J,IT,1)
            END DO
C
         END DO
         IF ( IPRINT.GT.0 ) THEN
            IOTMP = 300
            STR12 = '_DOS_'//'.dat'
            IF ( MPI ) THEN
               IF ( MPI_ID.LT.10 ) THEN
                  WRITE (IOUT,'(I1)') MPI_ID
                  J = 1
               ELSE
                  WRITE (IOUT,'(I2)') MPI_ID
                  J = 2
               END IF
               STR12 = '_DOS_'//IOUT(1:J)//'.dat'
            END IF
C
            DO IT = ITBOT,ITTOP
               IOTMP = IOTMP + IT
               N = NKMQ(IQAT(1,IT))
               IF ( ICALL.EQ.1 ) CALL OPENFILT(DATSET,LDATSET,TXT_T,
     &              LTXT_T,IT,STR12,12,FILNAM,LFN,2,IOTMP,
     &              '(l,ml,ms)-DOS',13,NTMAX)
               WRITE (IOTMP,99004) DREAL(ERYD-EFERMI)*RY_EV,
     &                             (DIMAG(DOS(J,J,IT)),J=1,N)
            END DO
         END IF
      END IF
C-----------------------------------------------------------------------
C
      RETURN
C
99001 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),/)
99002 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99003 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
99004 FORMAT (60F18.10)
C
      END
