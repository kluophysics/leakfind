C*==projtau.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the component projected TAU - matrices               *
C   *                                                                  *
C   *      TAU(IT) =  TAU(IQ) * ( 1 + (m(t)-m(c))*TAU(IQ) )**(-1)      *
C   *                                                                  *
C   *   NOTE: it is assumed that all equivalent sites  IQ  have the    *
C   *   same TAU-matrix  TAUQ(IQ). To get  TAU(IT)  the first site IQ  *
C   *   occupied by type IT is taken to be representative for          *
C   *   all other (NAT(IT)-1) sites occupied by IT                     *
C   *                                                                  *
C   *   allows an atom type IT to have different orientation of        *
C   *   its moment on different but equivalent sites  IQ               *
C   *                                                                  *
C   * 01/11/00                                                         *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX,NKMQ
      USE MOD_SITES,ONLY:NQMAX,IQAT,ICPA,IQBOT,IQTOP
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP
      USE MOD_CALCMODE,ONLY:IREL,
c modified by XJQ: scf of vibrations
     &                      task
c end-mod-xjq
      USE MOD_FILES,ONLY:WRTAU,WRTAUMQ,IFILTAU
      USE MOD_CONSTANTS,ONLY:C0,C1
c modified by XJQ: scf of vibrations
      USE MOD_THERMAL,ONLY:UMAT_VT,NVIBRA
      use mod_scfvb_cpa_sigma,only:lscfvb
c end-mod-xjq
      IMPLICIT NONE
C*--PROJTAU26
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='PROJTAU')
      REAL*8 TOL
      PARAMETER (TOL=1.0D-10)
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD
      INTEGER ICPAFLAG
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           TAUQ(NKMMAX,NKMMAX,NQMAX),TAUT(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
c modified by XJQ: scf of vibrations
      integer ivt_extend
      complex*16 mss_vt(nkmmax,nkmmax)
c end-mod-xjq
      REAL*8 CPAC
      COMPLEX*16 DMAMC(:,:),DMATT(:,:),DTILT(:,:),RMSS,RTAU
      INTEGER I,IA_ERR,ICPAF,IQ,IT,J,M,N
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DMAMC,DMATT,DTILT
C
      ALLOCATE (DMAMC(NKMMAX,NKMMAX))
      ALLOCATE (DMATT(NKMMAX,NKMMAX),DTILT(NKMMAX,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: DTILT')
C
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      LOOP_IT:DO IT = ITBOT,ITTOP
C
C ---------- pick first site IQ occupied by type IT to be representative
C ----------- all other (NAT(IT)-1) occupied sites have to be equivalent
C
         IQ = IQAT(1,IT)
         M = NKMMAX
         IF ( IREL.NE.2 ) THEN
            N = NKMQ(IQ)
         ELSE
            N = NKM
         END IF
C
         IF ( ICPA(IQ).EQ.1 ) THEN
C
c modified by XJQ: scf of vibrations
           if(lscfvb .and. task(1:5) /= 'SIGMA') then
             if(mod(it,nvibra)==0) then
               ivt_extend = (it/nvibra-1) * nvibra*nvibra + nvibra
             else
               ivt_extend = (it/nvibra) * nvibra*nvibra +
     &                      mod(it,nvibra)
             endif
             call cmat_u_trans(umat_vt(1,1,ivt_extend),msst(1,1,it),
     &                         'UAUT',mss_vt)
             CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,
     &                    MSSQ(1,1,IQ),mss_vt,M)
           else
             CALL GETDMAT(TAUQ(1,1,IQ),DMATT,DTILT,DMAMC,N,MSSQ(1,1,IQ),
     &                    MSST(1,1,IT),M)
c end-mod-xjq
           endif
C
C     ----------------------------------------
C              TAU(t) = TAU * D~(t)
C     ----------------------------------------
            CALL ZGEMM('N','N',N,N,N,C1,TAUQ(1,1,IQ),M,DTILT,M,C0,
     &                 TAUT(1,1,IT),M)
c modified by XJQ: scf of vibrations
            if(lscfvb .and. task(1:5) /= 'SIGMA') then
              call cmat_u_trans(umat_vt(1,1,ivt_extend),taut(1,1,it),
     &                          'UTAU',mss_vt)
              taut(:,:,it) = mss_vt(:,:)
            endif
c end-mod-xjq
C
            ICPAF = ICPAFLAG
            CPAC = CPACHNG
C
         ELSE
C
C-------------------------------------------- no CPA:  copy TAUQ to TAUT
C
            TAUT(:,:,IT) = TAUQ(:,:,IQ)
C
            ICPAF = 0
            CPAC = 0D0
C
         END IF
C
CWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
         IF ( WRTAU ) THEN
            WRITE (IFILTAU,99001) ERYD,IT,IQ,ICPAF,CPAC
            DO I = 1,N
               DO J = 1,N
                  IF ( I.EQ.J ) THEN
                     WRITE (IFILTAU,99003) I,J,TAUT(I,J,IT)
                  ELSE
                     IF ( CDABS(TAUT(I,J,IT)/TAUT(I,I,IT)).GT.TOL )
     &                    WRITE (IFILTAU,99003) I,J,TAUT(I,J,IT)
                  END IF
               END DO
            END DO
         END IF
CWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
      END DO LOOP_IT
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
CWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      IF ( WRTAUMQ ) THEN
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
         LOOP_IQ:DO IQ = IQBOT,IQTOP
            WRITE (IFILTAU,99002) ERYD,IQ,ICPAFLAG,CPACHNG
            DO I = 1,N
               DO J = 1,N
                  IF ( I.EQ.J ) THEN
                     WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),MSSQ(I,J,IQ)
                  ELSE
                     IF ( CDABS(TAUQ(I,I,IQ)).LT.1.D-12 ) THEN
                        RTAU = C0
                     ELSE
                        RTAU = TAUQ(I,J,IQ)/TAUQ(I,I,IQ)
                     END IF
C
                     RMSS = MSSQ(I,J,IQ)/MSSQ(I,I,IQ)
                     IF ( (CDABS(RTAU).GT.TOL) .OR. (CDABS(RMSS).GT.TOL)
     &                    ) WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),
     &                             MSSQ(I,J,IQ)
                  END IF
C
               END DO
            END DO
C
         END DO LOOP_IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      END IF
CWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
      IF ( WRTAU .OR. WRTAUMQ ) CALL TAU_WRITE(ICPAFLAG,CPACHNG,ERYD,
     &     MSSQ,TAUQ)
C
      DEALLOCATE (DMAMC,DMATT,DTILT,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC: DTILT')
C--------------------------------------------------------- FORMAT IFMT=2
C99001 FORMAT (/,80('*')/,2F21.15,' RYD   TAU FOR IT=',I2,'  IQ=',I2,:,
C     &        '  CPA:',I2,F15.6)
C99002 FORMAT (/,80('*')/,2F21.15,' RYD   TAU-C M-C  FOR IQ=',I2,:,
C     &        '  CPA:',I2,F15.6)
C--------------------------------------------------------- FORMAT IFMT=3
99001 FORMAT (/,80('*')/,2F21.15,' RYD   TAU FOR IT=',I5,'  IQ=',I5,:,
     &        '  CPA:',I2,F15.8)
99002 FORMAT (/,80('*')/,2F21.15,' RYD   TAU-C M-C  FOR IQ=',I5,:,
     &        '  CPA:',I2,F15.8)
99003 FORMAT (2I5,1P,4E22.14)
      END
C*==tau_write.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE TAU_WRITE(ICPAFLAG,CPACHNG,ERYD,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *   write the current q-dependent matrices TAUQ and MSSQ           *
C   *   to disk for later use                                          *
C   *                                                                  *
C   *   these matrices refer always to the GLOBAL frame                *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SITES,ONLY:NQMAX,IQBOT,IQTOP
      USE MOD_FILES,ONLY:IOTMP,LDATSET,DATSET
      USE MOD_ENERGY,ONLY:NETAB
      IMPLICIT NONE
C*--TAU_WRITE194
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_WRITE')
C
C Dummy arguments
C
      REAL*8 CPACHNG
      COMPLEX*16 ERYD
      INTEGER ICPAFLAG
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      CHARACTER*80 FILNAM
      INTEGER I,IQ,J,LFILNAM
      LOGICAL INITIALIZE
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      FILNAM = DATSET(1:LDATSET)//'_TAUQ_MSSQ_GLO.tau'
      LFILNAM = LDATSET + 18
C
C=======================================================================
      IF ( INITIALIZE ) THEN
C
         OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),FORM='UNFORMATTED',
     &         STATUS='REPLACE')
C
         WRITE (IOTMP) IQBOT,IQTOP,NKM,NKMMAX,NETAB(1)
C
         CLOSE (IOTMP)
C
         IF ( NETAB(1).LE.0 ) CALL STOP_MESSAGE(ROUTINE,'NETAB < 0')
C
         INITIALIZE = .FALSE.
C
      END IF
C=======================================================================
C
      OPEN (IOTMP,FILE=FILNAM(1:LFILNAM),FORM='UNFORMATTED',
     &      POSITION='APPEND')
C
      WRITE (IOTMP) ERYD,ICPAFLAG,CPACHNG
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      LOOP_IQ:DO IQ = IQBOT,IQTOP
C
         WRITE (IOTMP) IQ,((TAUQ(I,J,IQ),I=1,NKM),J=1,NKM),
     &                 ((MSSQ(I,J,IQ),I=1,NKM),J=1,NKM)
C
      END DO LOOP_IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
      CLOSE (IOTMP)
C
      END
C*==tau_read.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE TAU_READ(NE,ERYD,MSSQ,TAUQ)
C   ********************************************************************
C   *                                                                  *
C   *   READ  the current q-dependent matrices TAUQ and MSSQ           *
C   *                                                                  *
C   *   these matrices refer always to the GLOBAL frame                *
C   *                                                                  *
C   *   NE = 0:  read the header information and return                *
C   *                                                                  *
C   ********************************************************************
      USE MOD_ANGMOM,ONLY:NKM,NKMMAX
      USE MOD_SITES,ONLY:NQMAX,IQBOT,IQTOP
      USE MOD_FILES,ONLY:LDATSET,DATSET,IFILTAU
      USE MOD_ENERGY,ONLY:NETAB
      IMPLICIT NONE
C*--TAU_READ287
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='TAU_READ')
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER NE
      COMPLEX*16 MSSQ(NKMMAX,NKMMAX,NQMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
C
C Local variables
C
      REAL*8 CPACHNG
      CHARACTER*80 FILNAM
      INTEGER I,ICPAFLAG,IQ,IQBOT_IN,IQTOP_IN,IQ_IN,J,LFILNAM,NKMMAX_IN,
     &        NKM_IN
C
C*** End of declarations rewritten by SPAG
C
C=======================================================================
      IF ( NE.EQ.0 ) THEN
C
         CLOSE (IFILTAU*33)
C
         FILNAM = DATSET(1:LDATSET)//'_TAUQ_MSSQ_GLO.tau'
         LFILNAM = LDATSET + 18
C
         OPEN (IFILTAU*33,FILE=FILNAM(1:LFILNAM),FORM='UNFORMATTED',
     &         STATUS='OLD',ERR=50)
C
         READ (IFILTAU*33,ERR=50,END=50) IQBOT_IN,IQTOP_IN,NKM_IN,
     &         NKMMAX_IN,NE
C
         IF ( IQBOT_IN.NE.IQBOT .OR. IQTOP_IN.NE.IQTOP .OR. 
     &        NKM_IN.NE.NKM .OR. NKMMAX_IN.NE.NKMMAX .OR. NE.NE.NETAB(1)
     &        ) THEN
C
            WRITE (6,99001) IQBOT_IN,IQBOT,IQTOP_IN,IQTOP,NKM_IN,NKM,
     &                      NKMMAX_IN,NKMMAX,NE,NETAB(1)
C
            CALL STOP_MESSAGE(ROUTINE,'INPUT NOT OK')
C
         END IF
C
         IF ( NE.LE.0 ) CALL STOP_MESSAGE(ROUTINE,'NE < 0')
C
         WRITE (6,99003) ROUTINE(1:LEN_TRIM(ROUTINE)),FILNAM(1:LFILNAM),
     &                   NE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 50      CONTINUE
         WRITE (6,99004) ROUTINE(1:LEN_TRIM(ROUTINE)),FILNAM(1:LFILNAM)
C
         NE = 0
C
         RETURN
C
      END IF
C=======================================================================
C
      READ (IFILTAU*33) ERYD,ICPAFLAG,CPACHNG
C
      IF ( ICPAFLAG.NE.0 ) WRITE (6,99002) ERYD,ICPAFLAG,CPACHNG
C
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
      LOOP_IQ:DO IQ = IQBOT,IQTOP
C
         READ (IFILTAU*33) IQ_IN,((TAUQ(I,J,IQ),I=1,NKM),J=1,NKM),
     &                     ((MSSQ(I,J,IQ),I=1,NKM),J=1,NKM)
C
         IF ( IQ_IN.NE.IQ ) CALL STOP_MESSAGE(ROUTINE,'IQ ?')
      END DO LOOP_IQ
CQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
99001 FORMAT (/,10X,'input versus internal settings:',/,10X,
     &        'IQBOT:    ',2I10,/,10X,'IQTOP:    ',2I10,/,10X,
     &        'NKM:      ',2I10,/,10X,'NKMMAX:   ',2I10,/,10X,
     &        'NETAB(1): ',2I10)
99002 FORMAT (/,10X,'WARNING for ERYD = ',2F10.5,5X,'ICPAFLAG =',I2,5X,
     &        'CPACHNG =',F12.6)
99003 FORMAT (//,10X,'INFO from <',A,'>: ',//,10X,'reading TAU-file ',A,
     &        //,10X,'with   NE =',I4,'   energies',/)
99004 FORMAT (//,10X,'INFO from <',A,'>: ',//,10X,'TAU-file ',A,
     &        ' not available or empty ',//,10X,
     &        'NE set to 0 --- try old format next')
      END
