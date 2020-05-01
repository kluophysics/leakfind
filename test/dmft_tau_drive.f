C*==dmft_tau_drive.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_TAU_DRIVE(ERYD,IECURR,LSSITE,GFMAT,DMFTSIGMA,
     &                          NEMAX)
C
      USE MOD_CPA,ONLY:CPATOL,NCPA
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM,NKMQ,NLMQ,NLM
      USE MOD_TYPES,ONLY:NTMAX
      USE MOD_SITES,ONLY:NQMAX,IQBOT,IQTOP
      USE MOD_FILES,ONLY:IFILCBWF,IPRINT
      USE MOD_CALCMODE,ONLY:ORBPOL,DMFT,IREL
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG
      USE MOD_CONSTANTS,ONLY:C0
      IMPLICIT NONE
C*--DMFT_TAU_DRIVE14
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMPLEX*16 ERYD
      INTEGER IECURR,NEMAX
      LOGICAL LSSITE
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),
     &           GFMAT(NKMMAX,NKMMAX,NEMAX,NTMAX)
C
C Local variables
C
      LOGICAL CALCINT,GETIRRSOL
      REAL*8 CPACHNG,CPACHTAB(1)
      INTEGER I,ICPACONV,ICPAFLAG,IECPAFAIL(1),INFO,IPIV(:),IQ,IS,ITCPA,
     &        IWRIRRWF,IWRREGWF,J,N,NCPAFAIL
      COMPLEX*16 MAUX(:,:),MEZJ(:,:,:,:),MEZZ(:,:,:,:),MQS(:,:,:),
     &           MSSQ(:,:,:),MSST(:,:,:),P,PHASK(:),SSST(:,:,:),
     &           TAUQ(:,:,:),TAUT(:,:,:),TSSQ(:,:,:),TSST(:,:,:),W1(:,:)
C
C*** End of declarations rewritten by SPAG
C
C
      ALLOCATABLE MEZJ,MEZZ,MSSQ,MSST,PHASK,SSST,TAUQ,TAUT,TSST,TSSQ,W1
      ALLOCATABLE IPIV,MAUX,MQS
C
      ALLOCATE (MAUX(NKMMAX,NKMMAX),IPIV(NKMMAX))
C
      IF ( IREL.EQ.2 ) ALLOCATE (MQS(NLM,NLM,2))
C
      CALCINT = .TRUE.
      GETIRRSOL = .TRUE.
      NCPAFAIL = 0
      ICPAFLAG = 0
      IWRREGWF = 1
      IWRIRRWF = 1
C
      ALLOCATE (MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX))
      ALLOCATE (MSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (MSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (SSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TAUQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (TAUT(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TSST(NKMMAX,NKMMAX,NTMAX))
      ALLOCATE (TSSQ(NKMMAX,NKMMAX,NQMAX))
      ALLOCATE (W1(NKMMAX,NKMMAX))
      ALLOCATE (PHASK(NEMAX))
C
      TSSQ(1:NKMMAX,1:NKMMAX,1:NQMAX) = C0
C=======================================================================
C
C
      IF ( DMFT ) DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &     = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IECURR)
C=======================================================================
C
      CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFILCBWF,GETIRRSOL,ERYD,P,
     &              IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
      CALL MSSINIT(TSST,MSST,TSSQ,MSSQ)
C
      IF ( LSSITE ) THEN
         DO IQ = IQBOT,IQTOP
C
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
                  CALL ZGETRF(NLMQ(IQ),NLMQ(IQ),MQS(1,1,IS),NLM,IPIV,
     &                        INFO)
                  CALL ZGETRI(NLMQ(IQ),MQS(1,1,IS),NLM,IPIV,MAUX,
     &                        NLM*NLM,INFO)
               END DO
C
               DO J = 1,NLMQ(IQ)
                  CALL ZCOPY(NLMQ(IQ),MQS(1,J,1),1,TSSQ(1,J,IQ),1)
                  CALL ZCOPY(NLMQ(IQ),MQS(1,J,2),1,TSSQ(NLM+1,NLM+J,IQ),
     &                       1)
               END DO
C
            END IF
            TAUQ(1:N,1:N,IQ) = TSSQ(1:N,1:N,IQ)
         END DO
      ELSE
C
C=======================================================================
C
C
         CALL TAU_DRIVE(1,IPRINT,ERYD,P,TSSQ,MSSQ,TSST,MSST,TAUQ,
     &                  ICPAFLAG,CPACHNG,ITCPA,ICPACONV,PHASK)
C
         IF ( ICPAFLAG.NE.0 ) THEN
            NCPAFAIL = NCPAFAIL + 1
            CPACHTAB(NCPAFAIL) = CPACHNG
            IECPAFAIL(NCPAFAIL) = 1
         END IF
C
         IF ( NCPAFAIL.NE.0 ) THEN
            WRITE (6,99001) CPATOL,NCPAFAIL,
     &                      (IECPAFAIL(I),DREAL(ERYD),CPACHTAB(I),I=1,
     &                      NCPAFAIL)
            WRITE (6,'(1X,79(''*''),/)')
         ELSE IF ( NCPA.NE.0 ) THEN
            WRITE (6,99002)
         END IF
C
      END IF
C
C
      CALL PROJTAU(ICPAFLAG,CPACHNG,ERYD,MSST,MSSQ,TAUQ,TAUT)
C
C=======================================================================
C                     Green's function matrix
C=======================================================================
C
C
      CALL GFUNMAT_DRIVE(ERYD,IECURR,TAUT,MSST,GFMAT)
C
      RETURN
C
99001 FORMAT (/,1X,79('*'),/,' tolerance for CPA-cycle:',F15.7,/,
     &        ' CPA not converged for',I3,' energies:',/,
     &        3(' E:',I3,F7.4,' C:',F8.5,:,2X))
99002 FORMAT (/,1X,79('*'),/,25X,'no problems with','  CPA-cycle !',/,
     &        1X,79('*'),/)
C
      END
