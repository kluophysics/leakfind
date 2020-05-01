C*==spec_ssite_drive.f    processed by SPAG 6.70Rc at 16:36 on 28 Feb 2017
      SUBROUTINE SPEC_SSITE_DRIVE(P,TSST,RG,RF,HG,HF,ISTATE,LAY,ATOM,IO,
     &                            IECURR,RG1,RF1,MSST,NZERO)
C   ********************************************************************
C   *                                                                  *
C   *  driver program to run the ssite routines                        *
C   *                                                                  *
C   *  used as interface of the  JB- and SPRKKR-packages               *
C   *                                                                  *
C   ********************************************************************
C
C
C
      USE MOD_RMESH,ONLY:NMMAX,FULLPOT,JRNSMIN,JRNS1,JRCRI,JRWS,NRMAX
      USE MOD_FILES,ONLY:IFILCBWF,IFILGFWF
      USE MOD_TYPES,ONLY:NTMAX,ITBOT,ITTOP
      USE MOD_ANGMOM,ONLY:NMEMAX,NKMMAX,NKM
      USE MOD_CALCMODE,ONLY:DMFT,ITEST,LDAU
      USE MOD_SITES,ONLY:ITOQ
      USE MOD_TRANSPHO_LAYPOT,ONLY:IQ_SPR_LAYAT,IBLOCHTMP
      USE MOD_DMFT_LDAU,ONLY:DMFTSIG,DMFTSIGMA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER ATOM,IECURR,IO,ISTATE,LAY
      COMPLEX*16 P
      COMPLEX*16 HF(NRMAX,NKMMAX,NKMMAX),HG(NRMAX,NKMMAX,NKMMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),RF(NRMAX,NKMMAX,NKMMAX),
     &           RF1(NRMAX,NKMMAX,NKMMAX),RG(NRMAX,NKMMAX,NKMMAX),
     &           RG1(NRMAX,NKMMAX,NKMMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      REAL*8 NZERO(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER IA_ERR,ICALL,IE,IFIL,IM,IT,ITBOTTMP,ITTOPTMP,M
      COMPLEX*16 MEZJ(:,:,:,:),MEZZ(:,:,:,:)
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE MEZJ,MEZZ
C
C------------------------------------------ variables depending on NLMAX
C
      M = NKMMAX
      ALLOCATE (MEZJ(M,M,NTMAX,NMEMAX))
      ALLOCATE (MEZZ(M,M,NTMAX,NMEMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'ALLOC: SSITE_DRIVE -> MSSQ'
      ICALL = ICALL + 1
C
C=======================================================================
C                                DMFT
C=======================================================================
C
      IF ( .NOT.FULLPOT ) THEN
         DO IM = 1,NMMAX
            JRNS1(IM) = JRNSMIN
            JRCRI(IM) = JRWS(IM)
         END DO
      END IF
C=======================================================================
C
      ITEST = 0
C
C-----------------------------------------------------------------------
C
      IE = IECURR
C
      IF ( DMFT .OR. LDAU ) THEN
         IF ( ISTATE.EQ.2 ) THEN
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX) = C0
         ELSE IF ( IBLOCHTMP.EQ.0 ) THEN
            WRITE (*,*) 'Im(Sig)==0: in SSITE_DRIVE'
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &         = DCMPLX(DREAL(DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)),
     &         0.0D0)
         ELSE
            DMFTSIG(1:NKMMAX,1:NKMMAX,1:NTMAX)
     &         = DMFTSIGMA(1:NKMMAX,1:NKMMAX,1:NTMAX,IE)
C
         END IF
      END IF
C
      ITBOTTMP = ITBOT
      ITTOPTMP = ITTOP
C
      ITBOT = ITOQ(IO,IQ_SPR_LAYAT(ATOM,LAY))
      ITTOP = ITOQ(IO,IQ_SPR_LAYAT(ATOM,LAY))
C
      IF ( ISTATE.EQ.1 ) THEN
         IFIL = IFILCBWF
      ELSE
         IFIL = IFILGFWF
      END IF
C
C
C      CALL RUNSSITE(CALCINT,IWRREGWF,IWRIRRWF,IFIL,GETIRRSOL,ERYD,P,
C     &     IPRINT,TSST,MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
C
      DO IT = ITBOT,ITTOP
         CALL WAVFUN_ZJ_RH(IFIL,IT,P,MSST,TSST,RG,RF,HG,HF,ISTATE,RG1,
     &                     RF1,NZERO)
      END DO
C
CCC      IF (ISTATE.EQ.1)THEN
CCC         write(6,*)
CCC     &        "<SPEC_SSITE_DRIVE> d-wave function for initial
CCC     &              states suppressed!!"
CCC         RG(1:NRMAX,9:18,9:18)=C0
CCC         RF(1:NRMAX,9:18,9:18)=C0
CCC         HG(1:NRMAX,9:18,9:18)=C0
CCC         HF(1:NRMAX,9:18,9:18)=C0
CCC         RG1(1:NRMAX,9:18,9:18)=C0
CCC         RF1(1:NRMAX,9:18,9:18)=C0
CCC      END IF
C
      IF ( ITEST.EQ.6 ) THEN
         DO IT = ITBOT,ITTOP
C
            WRITE (6,99001) 'atom type   ',IT
C
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,3,3,1,
     &                      1.0D-9,6)
C
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,6)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,3,3,
     &                      1,1.0D-9,6)
         END DO
      END IF
      ITBOT = ITBOTTMP
      ITTOP = ITTOPTMP
C
C         END DO
C      END DO
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C                      energy - loop    END
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
      RETURN
99001 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),/)
      END
