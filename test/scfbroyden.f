C*==scfbroyden.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBROYDEN(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,BRYMIX,R,
     &                      R2DRDI,JRWS,IMT,VT,BT,RMSAVV,RMSAVB,NMMAX,
     &                      NRMAX,NTMAX)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme                                       *
C   *                                                                  *
C   * coded by       s.bluegel   iff juelich                           *
C   * modified       b.drittler  iff juelich                           *
C   * -----------------------------------------------------------------*
C   * sprkkr-version h.ebert     lmu munich    july 1994               *
C   *                                                                  *
C   * - no metric function  GMET  used anymore                         *
C   * - no weighting of RMSERR with concentration                      *
C   * - XINP, XOUT  1-dimensional vectors > no mapping in <SCFBROYPT2> *
C   *                                                                  *
C   * -----------------------------------------------------------------*
C   *                                                                  *
C   * imix=(0,..,6) : key for using different mixing schemes           *
C   *                          0 means straight mixing                 *
C   *                          3 broyden's first method used           *
C   *                          4 broyden's second method used          *
C   *                          5 strmix with given jacobian            *
C   *                            of broyden's first method             *
C   *                          6 strmix with given jacobian            *
C   *                            of broyden's second method            *
C   *                          7 anderson's method used                *
C   * icmbry=(0,1): key to store broyden's update vectors in           *
C   *               core-memory  (only used if imix is greater than 2) *
C   *                          0 means vectors stored on buffer memory *
C   *                          1 means vectors stored in core memory   *
C   *                                                                  *
C   * icmbry is fixed to 0                                             *
C   * temporary data is stored on scratch file IFILBROY_RHOPOT         *
C   * modify record length   lngr8   according to machine              *
C   *                                                                  *
C   * _________________________________________________________________*
C   *                                                                  *
C   * on entry: itrscf=1 VT,BT = VOLD,BOLD  initialize tables          *
C   *           else     VT,BT = VNEW,BNEW  according to RHO[VOLD,BOLD]*
C   *                                                                  *
C   * on exit:           VT,BT = V[mix], B[mix]                        *
C   * _________________________________________________________________*
C   *                                                                  *
C   * recommended (b.d.):    imix         4                            *
C   *                        istbry       1                            *
C   *                        itdept      40                            *
C   *                        brymix    5 .. 20 %                       *
C   *                        strmix    == brymix                       *
C   ********************************************************************
      USE MOD_TYPES,ONLY:ITBOT,ITTOP
      USE MOD_FILES,ONLY:IFILBROY_RHOPOT
      IMPLICIT NONE
C*--SCFBROYDEN55
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBROYDEN')
      INTEGER NSMAX,ITDEPTMAX,ITDCMMAX
      PARAMETER (NSMAX=2,ITDEPTMAX=40,ITDCMMAX=3)
C
C Dummy arguments
C
      REAL*8 BRYMIX,RMSAVB,RMSAVV
      INTEGER IMIX,IPRINT,ISTBRY,ITDEPT,ITRSCF,NMMAX,NRMAX,NTMAX
      REAL*8 BT(NRMAX,NTMAX),R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX),
     &       VT(NRMAX,NTMAX)
      INTEGER IMT(NTMAX),JRWS(NMMAX)
C
C Local variables
C
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),D1(:),D2(:),DEL1,DEL2,
     &       FM(:),FM1(:),MIXING,RMSERB,RMSERV,RSQ,SM(:),SM1(:),STRMIX,
     &       UI(:,:),VI(:,:),WN,WO,XINP(:),XOUT(:)
      INTEGER I,IA_ERR,IBB,ICMBRY,IM,IPF,IT,IVV,MIT,NMAP,NMAPMAX
      LOGICAL INITIALIZE
      SAVE AM,BM,DEL1,DEL2,FM,FM1,I,IBB,ICMBRY,IM,IPF,IT,IVV,MIT,MIXING,
     &     NMAP,NMAPMAX,RMSERB,RMSERV,RSQ,SM,SM1,STRMIX,WN,WO,XINP,XOUT
      EXTERNAL SCFBROYPT2
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
C
      ALLOCATABLE FM,SM,FM1,SM1,XINP,XOUT
      ALLOCATABLE D1,D2,UI,VI
C
      CALL TRACK_INFO(ROUTINE)
C
      ALLOCATE (D1(NRMAX),D2(NRMAX))
      ALLOCATE (UI(NRMAX*2*NTMAX,2:ITDCMMAX))
      ALLOCATE (VI(NRMAX*2*NTMAX,2:ITDCMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: VI')
C
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
      IF ( INITIALIZE ) THEN
         ALLOCATE (SM(NRMAX*2*NTMAX),SM1(NRMAX*2*NTMAX))
         ALLOCATE (FM(NRMAX*2*NTMAX),FM1(NRMAX*2*NTMAX))
         ALLOCATE (XINP(NRMAX*2*NTMAX),XOUT(NRMAX*2*NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: FM')
         XINP(1:NRMAX*2*NTMAX) = 999999D0
         INITIALIZE = .FALSE.
      END IF
Ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
C
      STRMIX = BRYMIX
      NMAPMAX = NRMAX*2*NTMAX
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
         IF ( (NSMAX.NE.2) .OR. (ITDCMMAX.NE.3) ) THEN
            WRITE (6,FMT=99002) NSMAX,ITDCMMAX
            CALL STOP_MESSAGE(ROUTINE,'NSMAX.NE.2 .OR. ITDCMMAX.NE.3')
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
            IF ( IMIX.EQ.5 .OR. IMIX.EQ.6 ) WRITE (6,FMT=99006) (IMIX-4)
     &           ,ITDEPT - 1,ISTBRY
C
            IF ( IMIX.EQ.7 ) WRITE (6,FMT=99009) ITDEPT,ISTBRY
C
            IF ( IMIX.NE.7 ) THEN
               WRITE (6,FMT=99004) STRMIX,BRYMIX
            ELSE
               WRITE (6,FMT=99010) STRMIX
            END IF
C
            IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY ) THEN
               MIXING = BRYMIX
C
            ELSE
C
               MIXING = STRMIX
            END IF
C
C     store pot[old]
C
            NMAP = 0
            IVV = 0
            DO IT = ITBOT,ITTOP
               IM = IMT(IT)
               NMAP = NMAP + NSMAX*JRWS(IM)
               IBB = IVV + JRWS(IM)
               DO I = 1,JRWS(IM)
                  RSQ = R(I,IM)*R(I,IM)
                  IVV = IVV + 1
                  XINP(IVV) = VT(I,IT)*RSQ
                  IBB = IBB + 1
                  XINP(IBB) = BT(I,IT)*RSQ
               END DO
               IVV = IVV + JRWS(IM)
            END DO
            IF ( NMAP.NE.IBB ) CALL STOP_MESSAGE(ROUTINE,'NMAP.NE.IBB')
            IF ( NMAP.GT.NMAPMAX )
     &            CALL STOP_MESSAGE(ROUTINE,'NMAP.GT.NMAPMAX')
C
C --------------------------------------------------------------- APOLLO
C              lngr8 =  length of  real*8
C              apollo domain:  2 (double) * 4 (bytes/32-bit word)
C              LNGR8 = 2*4
C              OPEN (IFILBROY_RHOPOT,RECL=2*NMAP*LNGR8,STATUS='SCRATCH',
C    _              FORM='UNFORMATTED')
C ------------------------------------------------------------------ IBM
            IF ( ICMBRY.EQ.0 ) OPEN (IFILBROY_RHOPOT,STATUS='SCRATCH',
     &                               FORM='UNFORMATTED')
C
         END IF
C
      ELSE
C   ====================================================================
C
C---> final construction of the charge density
C     first mixing scheme : straight mixing
C---> determination of the root mean square error
C
         RMSAVV = 0.0D0
         RMSAVB = 0.0D0
C
         WO = 1.0D0 - MIXING
         WN = MIXING
C
         IVV = 0
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IBB = IVV + JRWS(IM)
C
            DO I = 1,JRWS(IM)
               IVV = IVV + 1
               IBB = IBB + 1
C pot[new]
               RSQ = R(I,IM)*R(I,IM)
               XOUT(IVV) = VT(I,IT)*RSQ
               XOUT(IBB) = BT(I,IT)*RSQ
               DEL1 = (XOUT(IVV)-XINP(IVV))
               DEL2 = (XOUT(IBB)-XINP(IBB))
               D1(I) = DEL1*DEL1*R2DRDI(I,IM)
               D2(I) = DEL2*DEL2*R2DRDI(I,IM)
C pot[mix]
               XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
               XOUT(IBB) = WO*XINP(IBB) + WN*XOUT(IBB)
            END DO
            IVV = IVV + JRWS(IM)
C
            CALL RRADINT(IM,D1,RMSERV)
            CALL RRADINT(IM,D2,RMSERB)
            RMSERV = SQRT(RMSERV)
            RMSERB = SQRT(RMSERB)
            RMSAVV = DMAX1(RMSERV,RMSAVV)
            RMSAVB = DMAX1(RMSERB,RMSAVB)
            WRITE (IPF,FMT=99007) IT,RMSERV,RMSERB
         END DO
C
         WRITE (IPF,FMT=99008) ITRSCF,RMSAVV,RMSAVB
C
C----> broyden updating schemes
C
         IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY )
     &        CALL SCFBROYPT2(IPRINT,XINP,XOUT,FM,FM1,SM,SM1,UI,VI,AM,
     &        BM,MIXING,ITDEPT,ITDEPTMAX,IMIX,IFILBROY_RHOPOT,IPF,MIT,
     &        NMAP,NMAPMAX)
C
C----> reset to start new iteration     pot[old] = pot[new]
C
         IVV = 0
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IBB = IVV + JRWS(IM)
C
            DO I = 1,JRWS(IM)
               IVV = IVV + 1
               IBB = IBB + 1
C
               XINP(IVV) = XOUT(IVV)
               XINP(IBB) = XOUT(IBB)
               RSQ = R(I,IM)*R(I,IM)
               VT(I,IT) = XOUT(IVV)/RSQ
               BT(I,IT) = XOUT(IBB)/RSQ
C
            END DO
            IVV = IVV + JRWS(IM)
C
         END DO
C
      END IF
C   ====================================================================
C
C
      DEALLOCATE (D1,D2,UI,VI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
99001 FORMAT (/,1X,79('*'),/,28X,'BROYDEN - mixing scheme')
99002 FORMAT (/,' STOP IN  <BROYDEN>',/,' ARRAY LENGTHS OF LMTO-',
     &        ' AND BROYDEN-ROUTINES NOT COMPATIBLE !',/,
     &        '               EXPECTED     ACTUAL  ',/,
     &        ' NSMAX (SPIN)       2       ',I5,/,
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
     &        'then jacobian is fixed and potential',
     &        ' is updated using that jacobian',/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99007 FORMAT (5X,'rms-error for type',i3,':  V = ',1p,d11.4,2x,
     &        '   B = ',1p,d11.4)
99008 FORMAT (5X,'iter.',i4,'     average:  V = ',1p,d11.4,2x,'   B = ',
     &        1p,d11.4)
99009 FORMAT (/,10X,'andersons''s method used',/,10X,
     &        'iterationdepth to accumulate the jacobian:  ',i3,/,10X,
     &        'anderson used after ',i3,' iteration different mix')
99010 FORMAT (/,10X,'straight mixing factor used:    ',f12.7,/,10X,
     &        'anderson''s mixing afterwards',/,1X,79('*'),/)
      END
