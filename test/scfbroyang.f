C*==scfbroyang.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBROYANG(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,BRYMIX,
     &                      ITOQ,DROTQ,MVPHI,MVTET,MVGAM,QMPHI,QMTET,
     &                      QMGAM,ANGFM,ANGSM,ANGFM1,ANGSM1,ANGINP,
     &                      ANGOUT,NQ,NK,ERRAVANG,NQMAX,NTMAX,NMVECMAX,
     &                      NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme                                       *
C   *                                                                  *
C   * applied to mix the angles specifying                             *
C   * the local frame of reference                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:IFILBROY_ANG
      IMPLICIT NONE
C*--SCFBROYANG18
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER ITDEPTMAX,ITDCMMAX
      PARAMETER (ITDEPTMAX=40,ITDCMMAX=3)
C
C Dummy arguments
C
      REAL*8 BRYMIX,ERRAVANG
      INTEGER IMIX,IPRINT,ISTBRY,ITDEPT,ITRSCF,NK,NKMMAX,NMVECMAX,NQ,
     &        NQMAX,NTMAX
      REAL*8 ANGFM(3*NQMAX),ANGFM1(3*NQMAX),ANGINP(3*NQMAX),
     &       ANGOUT(3*NQMAX),ANGSM(3*NQMAX),ANGSM1(3*NQMAX),
     &       MVGAM(NTMAX,NMVECMAX),MVPHI(NTMAX,NMVECMAX),
     &       MVTET(NTMAX,NMVECMAX),QMGAM(NQMAX),QMPHI(NQMAX),
     &       QMTET(NQMAX)
      COMPLEX*16 DROTQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER ITOQ(NTMAX,NQMAX)
C
C Local variables
C
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),DELPHI,DELTET,ERRANG,
     &       MIXING,STRMIX,UI(:,:),VI(:,:),WN,WO
      INTEGER IANG,IA_ERR,ICMBRY,IMV,IPF,IQ,IT,MIT,NMAP,NMAPMAX
      SAVE AM,BM,ERRANG,IANG,ICMBRY,IPF,MIT,MIXING,NMAP,NMAPMAX,STRMIX,
     &     WN,WO
      EXTERNAL SCFBROYPT2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE UI,VI
C
      ALLOCATE (UI(3*NQMAX,2:ITDCMMAX))
      ALLOCATE (VI(3*NQMAX,2:ITDCMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'alloc:scfbroyden -> VI'
C
      STRMIX = BRYMIX
      NMAPMAX = 3*NQMAX
C
      IF ( ITDEPT.LT.3 ) ITDEPT = 40
C   ====================================================================
C
C------------ array set up and definition of input parameter -----------
C
      IF ( ITRSCF.EQ.0 ) THEN
         MIT = 1
         IPF = 6
         ICMBRY = 0
C
         WRITE (6,FMT=99001)
         IF ( ITDCMMAX.NE.3 ) THEN
            WRITE (6,FMT=99002) ITDCMMAX
            STOP
C
         ELSE
C
            IF ( ITDEPT.GT.ITDEPTMAX ) THEN
               WRITE (6,FMT=99003) ITDEPTMAX
               ITDEPT = ITDEPTMAX
C
            END IF
C
            IF ( IMIX.EQ.3 .OR. IMIX.EQ.4 ) WRITE (6,FMT=99005) (IMIX-2)
     &           ,ITDEPT,ISTBRY
C
            IF ( IMIX.GE.5 ) WRITE (6,FMT=99006) (IMIX-4),ITDEPT - 1,
     &                              ISTBRY
C
            WRITE (6,FMT=99004) STRMIX,BRYMIX
C
            IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY ) THEN
               MIXING = BRYMIX
C
            ELSE
C
               MIXING = STRMIX
            END IF
C
C ----------------------------------------------------- store old angles
C
            NMAP = NQ*3
            IANG = 0
            DO IQ = 1,NQ
               ANGINP(IANG+1) = QMPHI(IQ)
               ANGINP(IANG+2) = QMTET(IQ)
               ANGINP(IANG+3) = QMGAM(IQ)
               IANG = IANG + 3
            END DO
            IF ( NMAP.NE.IANG ) STOP 'NMAP <> IANG'
            IF ( NMAP.GT.NMAPMAX ) STOP 'NMAP > NMAPMAX'
C
            IF ( ICMBRY.EQ.0 ) OPEN (IFILBROY_ANG,STATUS='SCRATCH',
     &                               FORM='UNFORMATTED')
C
         END IF
C
      ELSE
C   ====================================================================
C
C------------------------- set new local frame of reference according to
C-------------------------------------------- orientation of spin moment
C
         IMV = 1
C
         DO IQ = 1,NQ
C
            IT = ITOQ(1,IQ)
C
            QMPHI(IQ) = MVPHI(IT,IMV)
            QMTET(IQ) = MVTET(IT,IMV)
            QMGAM(IQ) = MVGAM(IT,IMV)
C
         END DO
C
C=======================================================================
C
C---> final construction of the charge density
C     first mixing scheme : straight mixing
C---> determination of the root mean square error
C
         WO = 1.0D0 - MIXING
         WN = MIXING
C
         IANG = 0
         DO IQ = 1,NQ
            ERRANG = 0D0
C
            IANG = IANG + 1
            ANGOUT(IANG) = WO*ANGINP(IANG) + WN*QMPHI(IQ)
C
            IANG = IANG + 1
            ANGOUT(IANG) = WO*ANGINP(IANG) + WN*QMTET(IQ)
C
            IANG = IANG + 1
            ANGOUT(IANG) = WO*ANGINP(IANG) + WN*QMGAM(IQ)
C
         END DO
C
C----> broyden updating schemes
C
         IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY )
     &        CALL SCFBROYPT2(IPRINT,ANGINP,ANGOUT,ANGFM,ANGFM1,ANGSM,
     &        ANGSM1,UI,VI,AM,BM,MIXING,ITDEPT,ITDEPTMAX,IMIX,
     &        IFILBROY_ANG,IPF,MIT,NMAP,NMAPMAX)
C
C----> reset to start new iteration     ang[old] = ang[new]
C
         WRITE (6,99008)
C
         ERRAVANG = 0.0D0
C
         IANG = 0
         DO IQ = 1,NQ
C
            IF ( ANGOUT(IANG+1).LT.0D0 ) ANGOUT(IANG+1) = ANGOUT(IANG+1)
     &           + 360D0
            IF ( ANGOUT(IANG+1).GT.360D0 ) ANGOUT(IANG+1)
     &           = ANGOUT(IANG+1) - 360D0
            ANGOUT(IANG+2) = MAX(ANGOUT(IANG+2),0.0D0)
            ANGOUT(IANG+2) = MIN(ANGOUT(IANG+2),180.0D0)
            DELPHI = ABS(ANGINP(IANG+1)-ANGOUT(IANG+1))
            DELTET = ABS(ANGINP(IANG+2)-ANGOUT(IANG+2))
            ERRAVANG = MAX(ERRANG,DELPHI,DELTET)
C
            WRITE (6,99009) IQ,ANGINP(IANG+1),ANGINP(IANG+2),
     &                      ANGOUT(IANG+1),ANGOUT(IANG+2),DELPHI,DELTET
C
            IANG = IANG + 1
            ANGINP(IANG) = ANGOUT(IANG)
            QMPHI(IQ) = ANGOUT(IANG)
C
            IANG = IANG + 1
            ANGINP(IANG) = ANGOUT(IANG)
            QMTET(IQ) = ANGOUT(IANG)
C
            IANG = IANG + 1
            ANGINP(IANG) = ANGOUT(IANG)
            QMGAM(IQ) = ANGOUT(IANG)
C
            CALL ROTMAT(NK,3,QMPHI(IQ),QMTET(IQ),0.0D0,DROTQ(1,1,IQ),
     &                  NKMMAX)
C
         END DO
C
         WRITE (IPF,FMT=99007) ITRSCF,ERRAVANG
         WRITE (6,'(i5,4f10.3,'' #  ANGLES'')') ITRSCF,
     &          (QMPHI(IQ),QMTET(IQ),IQ=1,2)
C
      END IF
C   ====================================================================
C
      DEALLOCATE (UI,VI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) STOP 'dealloc:scfbroyden -> UI,VI'
C
99001 FORMAT (/,1X,79('*'),/,28X,'BROYDEN - mixing scheme')
99002 FORMAT (/,' STOP IN  <BROYDEN>',/,' ARRAY LENGTHS OF LMTO-',
     &        ' AND BROYDEN-ROUTINES NOT COMPATIBLE !',/,
     &        '               EXPECTED     ACTUAL  ',/,
     &        ' ITDCMMAX           3       ',I5)
99003 FORMAT (/,10X,'ITDEPT to large  array-length ITDEPTMAX=',I5,/,
     &        ' ITDEPT reduced accordingly ',/)
99004 FORMAT (/,10X,'straight mixing factor used:    ',f12.7,/,10X,
     &        'parameter for broyden-update:   ',f12.7,/,1X,79('*'),/)
99005 FORMAT (/,10X,'broyden''s method # ',i1,' used',/,10X,
     &        'iterationdepth to accumulate the jacobian:  ',i3,/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99006 FORMAT (/,10X,'broyden''s method # ',i3,
     &        ' is used up to iteration depth: ',i3,/,10X,
     &        'then jacobian is fixed and angles ',
     &        'are updated using that jacobian',/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99007 FORMAT (5X,'iter.',i4,'     average DEL = ',1p,d11.4)
99008 FORMAT (//,5X,' setting new LOCAL frame of reference ',//,5X,
     &        'IQ  old   phi      tet    --> new   phi      tet',
     &        '     del   phi      tet')
99009 FORMAT (5X,I2,4X,2F9.4,8X,2F9.4,5X,2F9.4)
      END
