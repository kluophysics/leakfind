C*==lloydpt2.f    processed by SPAG 6.70Rc at 15:37 on 19 Dec 2016
      SUBROUTINE LLOYDPT2(TOTNOS,CTOTDOS,SCLNOS,EBANDLD,EFERMILD,SCFMIX,
     &                    NCPAFAIL,IECPAFAIL,PHASK,PHAST,PHASA)
C   ********************************************************************
C   *                                                                  *
C   *  complete evaluation of the Lloyd formula                        *
C   *                                                                  *
C   *             take care of jumps in PHAST and PHASA and sum up     *
C   *             extrapolate NOS for E=EMIN and E=EFERMI on real axis *
C   *             using ONLY the data from the NETAB E-countour points *
C   *             via a Pade approximation                             *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  set scaling factor for charge to   SCLNOS = TOTNOSLLOYD/TOTNOS  *
C   *  with  TOTNOS  from Greens function data                         *
C   *                                                                  *
C   *  use  Pade approximation  to find an update EFERMILD             *
C   *  for the Fermi energy                                            *
C   *                                                                  *
C   *------------------------------------------------------------------*
C   *                                                                  *
C   *  to avoid overwriting data for last E-mesh point                 *
C   *  several arrays are allocated locally                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_KSPACE,ONLY:NKTAB
      USE MOD_FILES,ONLY:IFILDOS
      USE MOD_SITES,ONLY:IQAT
      USE MOD_TYPES,ONLY:NAT,CONC,NTMAX,NT,NVALTOT
      USE MOD_ANGMOM,ONLY:NKMMAX,NKMQ,NKM
      USE MOD_CALCMODE,ONLY:ITEST,IREL
      USE MOD_ENERGY,ONLY:SEARCHEF,NETAB,WETAB,ETAB,NEMAX,EFERMI,EMIN,
     &    IGRID
      USE MOD_CONSTANTS,ONLY:C0,PI
      IMPLICIT NONE
C*--LLOYDPT238
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NEXPOL
      PARAMETER (NEXPOL=5)
      REAL*8 ESTEPTOL
      PARAMETER (ESTEPTOL=1D-8)
C
C Dummy arguments
C
      REAL*8 EBANDLD,EFERMILD,SCFMIX,SCLNOS,TOTNOS
      INTEGER NCPAFAIL
      COMPLEX*16 CTOTDOS(NEMAX),PHASA(NKMMAX,NTMAX,NEMAX),PHASK(NEMAX),
     &           PHAST(NKMMAX,NTMAX,NEMAX)
      INTEGER IECPAFAIL(NEMAX)
C
C Local variables
C
      LOGICAL CHECK
      COMPLEX*16 CSSNOSLD(:),CSUMA,CSUMSS,CTOTDOSLD(:),CTOTNOSLD(:),
     &           DELNOS1,DELNOSEF,DOSEF,EBOT,ETABLD(:,:),ETOP,ETRY,
     &           NOSBOT,NOSLD1,NOSTOP,NOSTRY,PCOEF(:,:),PHASALD(:,:,:),
     &           PHASKLD(:),PHASTLD(:,:,:),SSNOSBOT,WETABLD(:),WT,Y1,
     &           YEF1,YEF2
      REAL*8 DOS1,DOS2,DQ,DQ0,ESTEP,NOSLD,SCLMIX,SCLPREV,SUMDOS,
     &       TOTDOSEF,TOTNOSLLOYD,WSPIN,X
      INTEGER I,I1,I1TOP,IA_ERR,ICALL,IE,IELD,IEPATH,IT,ITRY,N,NELD,
     &        NELDMAX,NPAD
C
C*** End of declarations rewritten by SPAG
C
      DATA ICALL/0/,SCLPREV/1.005D0/
C
      ALLOCATABLE PHASALD,PHASTLD,PHASKLD,ETABLD,CTOTNOSLD,CTOTDOSLD
      ALLOCATABLE                     PCOEF,            CSSNOSLD,WETABLD
C
      CHECK = .FALSE.
C      CHECK = .TRUE.
      ICALL = ICALL + 1
      IEPATH = 1
      SCLMIX = MAX(0.1D0,SCFMIX)
      IF ( IREL.GT.1 ) THEN
         WSPIN = 1D0
      ELSE
         WSPIN = 2D0
      END IF
C
      NELD = NETAB(1)
C
      NELDMAX = NELD
C
      ALLOCATE (WETABLD(NELDMAX))
      ALLOCATE (CSSNOSLD(NELDMAX))
      ALLOCATE (CTOTNOSLD(NELDMAX),CTOTDOSLD(NELDMAX))
      ALLOCATE (PHASKLD(NELDMAX),ETABLD(NELDMAX,2))
      ALLOCATE (PHASALD(NKMMAX,NTMAX,NELDMAX))
      ALLOCATE (PHASTLD(NKMMAX,NTMAX,NELDMAX))
      ALLOCATE (PCOEF(NELDMAX,NELDMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc: PCOEF -> LLOYDPT2'
C
C-------------------------------------------------------- for lloyd test
      SUMDOS = 0D0
      DOS1 = 0D0
C
C--------------------- copy from standard E-mesh IE to Lloyd E-mesh IELD
C---------------------- accept only E-points for which the CPA converged
C
      N = NKMMAX*NT
C
      IELD = 0
      DO IE = 1,NETAB(1)
         DO I = 1,NCPAFAIL
            IF ( IECPAFAIL(I).EQ.IE ) GOTO 100
         END DO
         IELD = IELD + 1
         ETABLD(IELD,1) = ETAB(IE,1)
         WETABLD(IELD) = WETAB(IE,1)
         PHASKLD(IELD) = PHASK(IE)
         CTOTDOSLD(IELD) = CTOTDOS(IE)
         CALL ZCOPY(N,PHAST(1,1,IE),1,PHASTLD(1,1,IELD),1)
         CALL ZCOPY(N,PHASA(1,1,IE),1,PHASALD(1,1,IELD),1)
 100  END DO
      NELD = IELD
C
C------------------------------------- remove jumps from PHAST and PHASA
C
      CALL SBRANCH(PHASTLD,NKM,NT,NKMMAX,NTMAX,NELD)
      CALL SBRANCH(PHASALD,NKM,NT,NKMMAX,NTMAX,NELD)
C
C----------------------------------------------- sum up to get CTOTNOSLD
C
      DO IELD = 1,NELD
         CSUMA = C0
         CSUMSS = C0
         DO IT = 1,NT
            N = NKMQ(IQAT(1,IT))
            WT = NAT(IT)*CONC(IT)
            DO I = 1,N
               CSUMSS = CSUMSS - 
     &                  WT*(PHASTLD(I,IT,IELD)-PHASALD(I,IT,IELD))
               CSUMA = CSUMA + WT*PHASALD(I,IT,IELD)
            END DO
         END DO
C
         CTOTNOSLD(IELD) = (PHASKLD(IELD)+CSUMSS)/PI
         CSSNOSLD(IELD) = CSUMA/PI
C
         IF ( ITEST.EQ.5 ) THEN
            X = DREAL(ETABLD(IELD,IEPATH))
            WRITE (62,'(5F12.5)') X,DIMAG(-CSUMA)
            WRITE (63,'(2F12.5)') X,DIMAG(CSUMSS)
            IF ( IELD.EQ.1 ) THEN
               DOS1 = DIMAG(CTOTDOSLD(IELD))
               WRITE (64,'(6F12.5)') X,0D0,0D0,0D0,NVALTOT,DOS1
            ELSE
               ESTEP = DREAL(ETABLD(IELD,IEPATH)-ETABLD(IELD-1,IEPATH))
               DOS2 = DIMAG(CTOTDOSLD(IELD))
               SUMDOS = SUMDOS + (DOS1+DOS2)*ESTEP/2D0
               DOS1 = DOS2
               NOSLD = DIMAG(CTOTNOSLD(IELD)-CTOTNOSLD(1))
               WRITE (64,'(6F12.5)') X,NOSLD,SUMDOS,NOSLD - SUMDOS,
     &                               NVALTOT,DOS1
            END IF
         END IF
C
      END DO
C
      IF ( .NOT.SEARCHEF ) RETURN
      IF ( ITEST.EQ.5 ) STOP 'Lloyd-test completed'
C
C ----------------------------------------------------------------------
C        select energy points to be used for Pade approximation
C                  use ONLY data from standard E-loop
C ----------------------------------------------------------------------
C
      I1 = 1
C
      NPAD = MIN(NEXPOL,NETAB(1))
C
C ----------------------------------------------------------------------
C                  deal with single site term  CSSNOSLD
C ----------------------------------------------------------------------
C
      CALL PADECOEFF(CSSNOSLD(I1),ETABLD(I1,IEPATH),PCOEF,NELDMAX,NPAD)
C
      EBOT = EMIN
      CALL PADEAPPROX(SSNOSBOT,EBOT,ETABLD(I1,IEPATH),PCOEF,NELDMAX,
     &                NPAD)
C
C ----------------------------------------------------------------------
C                  deal with full NOS   CTOTNOSLD
C                   extrapolate NOS to real axis
C ----------------------------------------------------------------------
C
      EBOT = EMIN
C
      CALL PADECOEFF(CTOTNOSLD(I1),ETABLD(I1,IEPATH),PCOEF,NELDMAX,NPAD)
C
      CALL PADEAPPROX(NOSBOT,EBOT,ETABLD(I1,IEPATH),PCOEF,NELDMAX,NPAD)
C
      ETOP = EFERMI
      I1TOP = I1 + NETAB(1) - 1 - NPAD + 1
C
      CALL PADECOEFF(CTOTNOSLD(I1TOP),ETABLD(I1TOP,IEPATH),PCOEF,
     &               NELDMAX,NPAD)
C
      CALL PADEAPPROX(NOSTOP,ETOP,ETABLD(I1TOP,IEPATH),PCOEF,NELDMAX,
     &                NPAD)
C
C ----------------------------------------------------------------------
C              first and last E-mesh points coincide with
C              start and end  of E-contour in case of IGRID = 9
C ----------------------------------------------------------------------
C
      IF ( IGRID(IEPATH).EQ.9 ) THEN
C
         NOSBOT = CTOTNOSLD(1)
C
         NOSTOP = CTOTNOSLD(NETAB(IEPATH))
C
      END IF
C
C ----------------------------------------------------------------------
C                   find new Fermi energy EFERMILD
C ----------------------------------------------------------------------
C
      ETRY = EFERMI
C
      DQ = NVALTOT/WSPIN - DIMAG(NOSTOP-NOSBOT)
      DQ0 = DQ
      ESTEP = SIGN(0.01D0,DQ)
C
      ITRY = 0
 200  CONTINUE
      DO I = 1,200
         ETRY = ETRY + ESTEP
         ITRY = ITRY + 1
C
         CALL PADEAPPROX(NOSTRY,ETRY,ETABLD(I1TOP,IEPATH),PCOEF,NELDMAX,
     &                   NPAD)
C
         DQ = NVALTOT/WSPIN - DIMAG(NOSTRY-NOSBOT)
         IF ( CHECK ) WRITE (6,99001) ITRY,DREAL(ETRY),ESTEP,DQ,
     &                                DIMAG(NOSTRY-NOSBOT)
C
         IF ( DQ*DQ0.LT.0D0 ) THEN
            ESTEP = -ESTEP/2D0
            DQ0 = DQ
            IF ( ABS(ESTEP).GT.ESTEPTOL ) GOTO 200
            ETRY = ETRY + ESTEP
            GOTO 300
         END IF
      END DO
C
C-----------------------------------------------------------------------
      WRITE (6,*) 'WARNING from <LLOYDPT2>: no Fermi energy found'
C
      TOTDOSEF = DREAL(CTOTDOS(NETAB(1)))
      ETRY = EFERMI - (TOTNOS-NVALTOT/WSPIN)/TOTDOSEF
C
C     STOP 'in <LLOYDPT2> - no Fermi energy found'
C-----------------------------------------------------------------------
 300  CONTINUE
      EFERMILD = DREAL(ETRY)
C
      WRITE (6,99005) EFERMILD
C
C ------------------------------------------------ calculate band energy
C NOTE:  E_band = Imag( EBANDLD )
C
      EBANDLD = 0D0
      DO IELD = 1,NELD
         EBANDLD = EBANDLD - 
     &             DIMAG(WETABLD(IELD)*(CTOTNOSLD(IELD)-NOSBOT))
      END DO
      EBANDLD = NVALTOT*EFERMI + WSPIN*EBANDLD
C
C--------------------------------------------------- normalize CTOTNOSLD
C
      NOSLD1 = CTOTNOSLD(1)
      DO IELD = 1,NELD
         CTOTNOSLD(IELD) = CTOTNOSLD(IELD) - NOSLD1
      END DO
C
C--------------- delta NOS between real axis and 1st / last E-mesh point
C
      DELNOS1 = CTOTDOSLD(1)*(ETABLD(1,IEPATH)-EMIN)*0.5D0
      DELNOSEF = CTOTDOSLD(NELD)*(EFERMI-ETABLD(NELD,IEPATH))
C
C ---------------------------------------- total NOS from EMIN to EFERMI
C
      TOTNOSLLOYD = DIMAG(NOSTOP-NOSBOT)
C
C ---------------------- scale quantities derived from Green''s function
C
      SCLNOS = TOTNOSLLOYD/TOTNOS
      IF ( ICALL.NE.1 ) SCLNOS = SCLMIX*SCLNOS + (1D0-SCLMIX)*SCLPREV
C
      IF ( ABS(SCLNOS-1D0).GT.0.3D0 ) THEN
         SCLPREV = SCLNOS
         SCLNOS = 1D0 + SIGN(0.3D0,(SCLNOS-1D0))
         WRITE (6,99004) SCLPREV,SCLNOS
      END IF
C
      SCLPREV = SCLNOS
C
      WRITE (6,99003) DELNOS1,ETABLD(1,IEPATH),DELNOSEF,
     &                ETABLD(NELD,IEPATH),NOSTOP,NOSBOT,NOSTOP - NOSBOT,
     &                SSNOSBOT,CTOTNOSLD(NELD),TOTNOSLLOYD,TOTNOS,
     &                TOTNOSLLOYD - TOTNOS,SCLNOS
C
C ======================================================================
      IF ( CHECK ) THEN
C ----------------------------------------------------------------------
C                             extrapolate NOS
C ----------------------------------------------------------------------
C
         CALL PADEAPPROX(Y1,ETABLD(1,IEPATH),ETABLD(I1,IEPATH),PCOEF,
     &                   NELDMAX,NPAD)
C
         WRITE (6,99002) ' E(1)     ',ETABLD(1,IEPATH)
         WRITE (6,99002) ' N (Pade) ',Y1 - NOSLD1
         WRITE (6,*) '   '
C
         CALL PADEAPPROX(YEF1,ETABLD(NELD,IEPATH),ETABLD(I1,IEPATH),
     &                   PCOEF,NELDMAX,NPAD)
C
         WRITE (6,99002) ' E (TOP)  ',ETABLD(NELD,IEPATH)
         WRITE (6,99002) ' N (Pade) ',YEF1 - NOSLD1
         WRITE (6,99002) '   '
C
         ETRY = EFERMI
         CALL PADEAPPROX(YEF2,ETRY,ETABLD(I1,IEPATH),PCOEF,NELDMAX,NPAD)
C
         WRITE (6,99002) ' EFERMI   ',EFERMI
         WRITE (6,99002) ' N (Pade) ',YEF2 - NOSLD1
         WRITE (6,*) '   '
C
C ----------------------------------------------------------------------
C                    extrapolate DOS to Fermi energy
C ----------------------------------------------------------------------
C
         CALL PADECOEFF(CTOTDOSLD(I1),ETABLD(I1,IEPATH),PCOEF,NELDMAX,
     &                  NPAD)
C
         CALL PADEAPPROX(DOSEF,ETRY,ETABLD(I1,IEPATH),PCOEF,NELDMAX,
     &                   NPAD)
C
         WRITE (6,99002) ' DOS '
         WRITE (6,99002) ' n (Pade) ',DOSEF
         WRITE (6,99002) ' n (TOP)  ',CTOTDOSLD(NELD-1)
         WRITE (6,*) '   '
C
C ----------------------------------------------------------------------
C
         WRITE (6,'(''     XXLOOPXX'',i5,10f15.7)') NKTAB,
     &          DIMAG(CTOTNOSLD(NELD)),DIMAG(YEF1),DIMAG(YEF2),TOTNOS,
     &          DIMAG(CTOTDOSLD(NELD-1)),DIMAG(DOSEF),EFERMI,DREAL(ETRY)
         WRITE (6,*) '   '
C
C ----------------------------------------------------------------------
C
         OPEN (IFILDOS,FILE='AA_lloyd1.dat')
         DO IELD = 1,NELD
            WRITE (IFILDOS,'(5F12.5)') ETABLD(IELD,IEPATH),
     &                                 CTOTNOSLD(IELD)
         END DO
C
         OPEN (IFILDOS,FILE='AA_lloyd2.dat')
         X = 0D0
         WRITE (IFILDOS,'(5F12.5)') X,CTOTNOSLD(1)
         DO IELD = 2,NELD
            X = X + ABS(ETABLD(IELD,IEPATH)-ETABLD(IELD-1,IEPATH))
            WRITE (IFILDOS,'(5F12.5)') X,CTOTNOSLD(IELD)
         END DO
C
         CLOSE (IFILDOS)
C
      END IF
C ======================================================================
C
      DEALLOCATE (PHASALD,PHASKLD,PHASTLD,ETABLD,CTOTNOSLD,CTOTDOSLD)
      DEALLOCATE (PCOEF,CSSNOSLD,WETABLD)
C
C ----------------------------------------------------------------------
99001 FORMAT (2X,'<LLOYDPT2>: update E_F',I3,5F13.8)
99002 FORMAT (2X,'<LLOYDPT2>',A,5F15.10)
99003 FORMAT (/,1X,79('-'),/,10X,'evaluation of Lloyd''s formula ',/,
     &        10X,'DELNOS 1     ',2F12.6,5X,'E ',2F10.6,/,10X,
     &        'DELNOS EF    ',2F12.6,5X,'E ',2F10.6,/,10X,
     &        'TOTNOS TOP   ',2F12.6,/10X,'TOTNOS BOT   ',4F12.6,/,10X,
     &        'TOTNOS BOT SS',2F12.6,/,10X,'TOTNOS LLOYD ',2F12.6,4X,
     &        'N(NELD)   - N(1)    (Lloyd)',/,10X,'TOTNOS SUM   ',12X,
     &        F12.6,4X,'N(EFERMI) - N(EMIN) (Lloyd)',/,10X,
     &        'TOTNOS DOS   ',12X,F12.6,4X,'Green''s function  '/,10X,
     &        'TOTNOS DIFF  ',12X,F12.6,/,10X,'scale  LLOYD ',12X,F12.6,
     &        /,1X,79('-'),/)
99004 FORMAT (/,1X,79('!'),/,17X,'WARNING from  <LLOYDPT2>',/,5X,
     &        'Lloyd scaling factor questionable ',/,5X,'reduced from ',
     &        f8.3,'  to',f8.3,/,1X,79('!'),/)
99005 FORMAT (/' <LLOYDPT2>:  (2)      new Fermi energy',F11.5,/)
      END
