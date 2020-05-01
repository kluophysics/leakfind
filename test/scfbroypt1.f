C*==scfbroypt1.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCFBROYPT1(IPRINT,ITRSCF,IMIX,ISTBRY,ITDEPT,SCFMIX,VT,
     &                      BT,VNST,BNST,RMSAVV,RMSAVB,JRNS1,JRNSMIN)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme                                       *
C   *                                                                  *
C   * coded by       s.bluegel   iff juelich                           *
C   * modified       b.drittler  iff juelich                           *
C   * fp-spr-kkr-version h.ebert     lmu munich    july 2002           *
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
C   *                        scfmix    5 .. 20 %                       *
C   *                        strmix    == scfmix                       *
C   ********************************************************************
      USE MOD_RMESH,ONLY:NM,NMMAX,NRMAX,R,R2DRDI,JRWS,JRCRI,FULLPOT
      USE MOD_TYPES,ONLY:IMT,ITBOT,ITTOP,NLMFPT,KLMFP,NTMAX,NT,NLMFPMAX
      USE MOD_FILES,ONLY:IFILBROY_RHOPOT
      IMPLICIT NONE
C*--SCFBROYPT149
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCFBROYPT1')
      INTEGER NSMAX,ITDEPTMAX,ITDCMMAX
      PARAMETER (NSMAX=2,ITDEPTMAX=40,ITDCMMAX=3)
C
C Dummy arguments
C
      INTEGER IMIX,IPRINT,ISTBRY,ITDEPT,ITRSCF,JRNSMIN
      REAL*8 RMSAVB,RMSAVV,SCFMIX
      REAL*8 BNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),BT(NRMAX,NTMAX),
     &       VNST(JRNSMIN:NRMAX,NLMFPMAX,NTMAX),VT(NRMAX,NTMAX)
      INTEGER JRNS1(NMMAX)
C
C Local variables
C
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),D1(:),D2(:),DEL1,DEL2,
     &       FM(:),FM1(:),MIXING,RMSERB,RMSERBNS,RMSERV,RMSERVNS,
     &       RSQ(NRMAX,NMMAX),SM(:),SM1(:),STRMIX,UI(:,:),VI(:,:),WN,WO,
     &       XINP(:),XOUT(:)
      INTEGER IA_ERR,IBB,ICMBRY,IM,IPF,IR,IRTOP,IT,IVV,LM,MIT,NMAP,
     &        NMAPMAX,NRNS
      SAVE AM,BM,FM,FM1,ICMBRY,IPF,MIT,MIXING,NMAP,NMAPMAX,RMSERB,
     &     RMSERV,SM,SM1,STRMIX,XINP,XOUT
      EXTERNAL SCFBROYPT2
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE D1,D2,UI,VI,FM,SM,FM1,SM1,XINP,XOUT
C
      DO IM = 1,NM
         DO IR = 1,NRMAX
            RSQ(IR,IM) = R(IR,IM)*R(IR,IM)
         END DO
      END DO
C
      IF ( ITRSCF.EQ.0 ) THEN
C
         NMAP = 0
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
            NMAP = NMAP + NSMAX*IRTOP
C
            IF ( FULLPOT ) THEN
               NRNS = JRCRI(IM) - JRNS1(IM) + 1
               DO LM = 2,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) NMAP = NMAP + 2*NRNS
               END DO
            END IF
C
         END DO
C
         NMAPMAX = NMAP
C
C----------------------------- auxilary arrays to be ALLOCATED AND SAVED
C
         IF ( ALLOCATED(FM) ) DEALLOCATE (FM,SM,FM1,SM1,XINP,XOUT)
C
         ALLOCATE (FM(NMAPMAX),SM(NMAPMAX),FM1(NMAPMAX))
         ALLOCATE (SM1(NMAPMAX),XINP(NMAPMAX),XOUT(NMAPMAX))
C
      END IF
C
      ALLOCATE (D1(NRMAX),UI(NMAPMAX,2:ITDCMMAX))
      ALLOCATE (D2(NRMAX),VI(NMAPMAX,2:ITDCMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: VI')
C
      STRMIX = SCFMIX
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
         IF ( IMIX.LT.3 ) THEN
            WRITE (6,FMT=99001) 'STRAIGHT'
         ELSE
            WRITE (6,FMT=99001) 'BROYDEN'
         END IF
         IF ( ITDCMMAX.NE.3 ) CALL STOP_MESSAGE(ROUTINE,'ITDCMMAX != 3')
         IF ( ITDEPT.GT.ITDEPTMAX ) THEN
            WRITE (6,FMT=99002) ITDEPTMAX
            ITDEPT = ITDEPTMAX
         END IF
         IF ( IMIX.EQ.3 .OR. IMIX.EQ.4 ) WRITE (6,FMT=99005) (IMIX-2),
     &        ITDEPT,ISTBRY
         IF ( IMIX.GE.5 ) WRITE (6,FMT=99006) (IMIX-4),ITDEPT - 1,ISTBRY
         IF ( IMIX.LT.3 ) THEN
            WRITE (6,FMT=99003) STRMIX
         ELSE
            WRITE (6,FMT=99004) STRMIX,SCFMIX
         END IF
C
         IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY ) THEN
            MIXING = SCFMIX
         ELSE
            MIXING = STRMIX
         END IF
C
C     store pot[old]
C
         NMAP = 0
         IVV = 0
         DO IT = ITBOT,ITTOP
            IM = IMT(IT)
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
            NMAP = NMAP + NSMAX*IRTOP
            IBB = IVV + IRTOP
            DO IR = 1,IRTOP
               IVV = IVV + 1
               IBB = IBB + 1
               XINP(IVV) = VT(IR,IT)*RSQ(IR,IM)
               XINP(IBB) = BT(IR,IT)*RSQ(IR,IM)
            END DO
            IVV = IVV + IRTOP
C
            IF ( FULLPOT ) THEN
               NRNS = JRCRI(IM) - JRNS1(IM) + 1
               IBB = IVV + IRTOP - NRNS
               DO LM = 2,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     NMAP = NMAP + 2*NRNS
                     IBB = IVV + NRNS
                     DO IR = JRNS1(IM),JRCRI(IM)
                        IVV = IVV + 1
                        IBB = IBB + 1
                        XINP(IVV) = VNST(IR,LM,IT)*RSQ(IR,IM)
                        XINP(IBB) = BNST(IR,LM,IT)*RSQ(IR,IM)
                     END DO
                     IVV = IVV + NRNS
                  END IF
               END DO
            END IF
C
         END DO
C
         IF ( NMAP.NE.IBB ) CALL STOP_MESSAGE(ROUTINE,'NMAP.NE.IBB')
         IF ( NMAP.GT.NMAPMAX )
     &         CALL STOP_MESSAGE(ROUTINE,'NMAP.GT.NMAPMAX')
C
C --------------------------------------------------------------- APOLLO
C              lngr8 =  length of  real*8
C              apollo domain:  2 (double) * 4 (bytes/32-bit word)
C              LNGR8 = 2*4
C              OPEN (IFILBROY_RHOPOT,RECL=2*NMAP*LNGR8,STATUS='SCRATCH',
C    _              FORM='UNFORMATTED')
C ------------------------------------------------------------------ IBM
         IF ( ICMBRY.EQ.0 ) OPEN (IFILBROY_RHOPOT,STATUS='SCRATCH',
     &                            FORM='UNFORMATTED')
C
      ELSE
C   ====================================================================
C
C---> final construction of the potential
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
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
            IBB = IVV + IRTOP
C
            DO IR = 1,IRTOP
               IVV = IVV + 1
               IBB = IBB + 1
C pot[new]
               XOUT(IVV) = VT(IR,IT)*RSQ(IR,IM)
               XOUT(IBB) = BT(IR,IT)*RSQ(IR,IM)
               DEL1 = (XOUT(IVV)-XINP(IVV))
               DEL2 = (XOUT(IBB)-XINP(IBB))
               D1(IR) = DEL1*DEL1*R2DRDI(IR,IM)
               D2(IR) = DEL2*DEL2*R2DRDI(IR,IM)
C pot[mix]
               XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
               XOUT(IBB) = WO*XINP(IBB) + WN*XOUT(IBB)
            END DO
            DO IR = IRTOP + 1,NRMAX
               D1(IR) = 0D0
               D2(IR) = 0D0
            END DO
            IVV = IVV + IRTOP
C
C-------------------- NOTE: use of <RINTSIMP> not consistent for FULLPOT
C-------- should not cause problems because RMSERV should go to 0 anyhow
C
            CALL RRADINT(IM,D1,RMSERV)
            CALL RRADINT(IM,D2,RMSERB)
C
            RMSERV = SQRT(RMSERV)
            RMSERB = SQRT(RMSERB)
C
            IF ( FULLPOT ) THEN
               NRNS = JRCRI(IM) - JRNS1(IM) + 1
               IBB = IVV + IRTOP - NRNS
               DO LM = 2,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     IBB = IVV + NRNS
                     DO IR = JRNS1(IM),JRCRI(IM)
                        IVV = IVV + 1
                        IBB = IBB + 1
C
C pot[new]
                        XOUT(IVV) = VNST(IR,LM,IT)*RSQ(IR,IM)
                        XOUT(IBB) = BNST(IR,LM,IT)*RSQ(IR,IM)
                        DEL1 = (XOUT(IVV)-XINP(IVV))
                        DEL2 = (XOUT(IBB)-XINP(IBB))
                        D1(IR) = DEL1*DEL1*R2DRDI(IR,IM)
                        D2(IR) = DEL2*DEL2*R2DRDI(IR,IM)
C pot[mix]
                        XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
                        XOUT(IBB) = WO*XINP(IBB) + WN*XOUT(IBB)
                     END DO
                     IVV = IVV + NRNS
C
                     CALL RRADINT(IM,D1,RMSERVNS)
                     CALL RRADINT(IM,D2,RMSERBNS)
C
                     RMSERV = MAX(RMSERV,SQRT(RMSERVNS))
                     RMSERB = MAX(RMSERB,SQRT(RMSERBNS))
                  END IF
               END DO
            END IF
C
            IF ( NT.GT.1 ) WRITE (IPF,FMT=99007) IT,RMSERV,RMSERB
            RMSAVV = DMAX1(RMSERV,RMSAVV)
            RMSAVB = DMAX1(RMSERB,RMSAVB)
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
            IF ( FULLPOT ) THEN
               IRTOP = JRCRI(IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
            IBB = IVV + IRTOP
C
            DO IR = 1,IRTOP
               IVV = IVV + 1
               IBB = IBB + 1
C
               XINP(IVV) = XOUT(IVV)
               XINP(IBB) = XOUT(IBB)
               VT(IR,IT) = XOUT(IVV)/RSQ(IR,IM)
               BT(IR,IT) = XOUT(IBB)/RSQ(IR,IM)
C
            END DO
            IVV = IVV + IRTOP
C
            IF ( FULLPOT ) THEN
               NRNS = JRCRI(IM) - JRNS1(IM) + 1
               IBB = IVV + IRTOP - NRNS
               DO LM = 2,NLMFPT(IT)
                  IF ( KLMFP(LM,IT).NE.0 ) THEN
                     IBB = IVV + NRNS
                     DO IR = JRNS1(IM),JRCRI(IM)
                        IVV = IVV + 1
                        IBB = IBB + 1
C
                        XINP(IVV) = XOUT(IVV)
                        XINP(IBB) = XOUT(IBB)
                        VNST(IR,LM,IT) = XOUT(IVV)/RSQ(IR,IM)
                        BNST(IR,LM,IT) = XOUT(IBB)/RSQ(IR,IM)
                     END DO
                     IVV = IVV + NRNS
                  END IF
               END DO
            END IF
C
         END DO
C
      END IF
C   ====================================================================
C
      DEALLOCATE (D1,D2,UI,VI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
C   ====================================================================
99001 FORMAT (/,1X,79('*'),/,20X,A,' - mixing scheme for SCF - cycle')
99002 FORMAT (/,10X,'ITDEPT to large  array-length ITDEPTMAX=',I5,/,
     &        ' ITDEPT reduced accordingly ',/)
99003 FORMAT (/,10X,'straight mixing factor used:    ',f12.7,/1X,79('*')
     &        ,/)
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
     &        '   B = ',1P,D11.4)
99008 FORMAT (5X,'iter.',i4,'     average:  V = ',1p,d11.4,2x,'   B = ',
     &        1P,D11.4)
      END
