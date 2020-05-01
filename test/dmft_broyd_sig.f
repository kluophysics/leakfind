C*==dmft_broyd_sig.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE DMFT_BROYD_SIG(ITRSCF,IMIX,ISTBRY,ITDEPT,BRYMIX,
     &                          DMFTSIGMA,SCFDONE,WETAB,SHFTEF,RMSAVV)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme                                       *
C   *                                                                  *
C   * applied to mix the LDA+U and DMFT self energies                  *
C   *                                                                  *
C   ********************************************************************
C
C
      USE MOD_ENERGY,ONLY:NEMAX
      USE MOD_ANGMOM,ONLY:NKMMAX
      USE MOD_TYPES,ONLY:NTMAX,LOPT
      USE MOD_SCF,ONLY:SCFTOL
      USE MOD_DMFT_LDAU,ONLY:KSELF,SIGTOL,DMFT_FIX_DYN_SE
      USE MOD_CALCMODE,ONLY:LDAU,DMFT
      USE MOD_CONSTANTS,ONLY:C0
      USE MOD_LATTICE,ONLY:SUB_SYSTEM,SYSTEM_DIMENSION,SYSTEM_TYPE
      IMPLICIT NONE
C*--DMFT_BROYD_SIG22
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
      REAL*8 BRYMIX,RMSAVV,SHFTEF
      INTEGER IMIX,ISTBRY,ITDEPT,ITRSCF
      LOGICAL SCFDONE
      COMPLEX*16 DMFTSIGMA(NKMMAX,NKMMAX,NTMAX,NEMAX),WETAB(NEMAX)
C
C Local variables
C
      REAL*8 AM(:),BM(:),MIXING,RMSAVSIG,RMSERT,RMSERTREL,SIGFM(:),
     &       SIGFM1(:),SIGINP(:),SIGOUT(:),SIGSM(:),SIGSM1(:),STRMIX,
     &       UI(:,:),VI(:,:),WN,WO
      COMPLEX*16 CDEL,SIGI(:,:,:),SIGO(:,:,:)
      INTEGER IA_ERR,ICMBRY,IE,IK,IK1,IK2,IOBROY,IPF,ISIG,IT,MIT,NEBROY,
     &        NMAP,NMAPMAX
      SAVE AM,BM,IOBROY,MIT,MIXING,RMSAVSIG,SIGFM,SIGFM1,SIGI,SIGINP,
     &     SIGO,SIGOUT,SIGSM,SIGSM1,UI,VI
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SIGINP,SIGOUT,UI,VI,AM,BM
      ALLOCATABLE SIGFM,SIGFM1,SIGSM,SIGSM1,SIGI,SIGO
C
C
      IF ( LDAU ) THEN
         NEBROY = 1
      ELSE
         NEBROY = NEMAX
      END IF
C
      ISIG = 0
      DO IE = 1,NEBROY
         DO IT = 1,NTMAX
            IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
               DO IK1 = 1,NKMMAX
                  DO IK2 = 1,NKMMAX
                     ISIG = ISIG + 2
                  END DO
               END DO
            END IF
         END DO
      END DO
      NMAP = ISIG+1 !REC + one because of compiler problem
      NMAPMAX = NMAP
C
      STRMIX = BRYMIX
      IF ( .NOT.ALLOCATED(SIGINP) ) THEN
         ALLOCATE (SIGINP(NMAP),SIGOUT(NMAP))
         ALLOCATE (UI(NMAPMAX,2:ITDCMMAX))
         ALLOCATE (SIGFM(NMAPMAX),SIGFM1(NMAPMAX))
         ALLOCATE (SIGSM(NMAPMAX),SIGSM1(NMAPMAX))
         ALLOCATE (AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1))
         ALLOCATE (SIGO(NKMMAX*NKMMAX+1,NTMAX,NEMAX)) ! REC +1 because of compiler problem
         ALLOCATE (SIGI(NKMMAX*NKMMAX+1,NTMAX,NEMAX)) ! REC +1 because of compiler problem
         ALLOCATE (VI(NMAPMAX,2:ITDCMMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) STOP 'alloc:scfbroyden -> VI'
         SIGINP = 0.0D0
         SIGOUT = 0.0D0
         SIGFM = 0.0D0
         SIGFM1 = 0.0D0
         SIGSM = 0.0D0
         SIGSM1 = 0.0D0
         VI = 0.0D0
         UI = 0.0D0
         AM = 0.0D0
         BM = 0.0D0
      END IF
C
      IPF = 6
      IOBROY = 91
      ICMBRY = 0
      IF ( ITRSCF.EQ.1 ) THEN
         MIT = 1
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
         END IF
C
C
C ----------------------------------------------------- store old sigma
C
C
         ISIG = 0
         DO IE = 1,NEBROY
            DO IT = 1,NTMAX
               IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                  DO IK1 = 1,NKMMAX
                     DO IK2 = 1,NKMMAX
                        ISIG = ISIG + 1
                        SIGINP(ISIG) = DREAL(DMFTSIGMA(IK2,IK1,IT,IE))
                        ISIG = ISIG + 1
                        SIGINP(ISIG) = DIMAG(DMFTSIGMA(IK2,IK1,IT,IE))
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
         IF ( ISIG.NE.NMAP ) STOP '<DMFT_BROYD_SIG>: ISIG .NE. NMAP'
C
         IF ( ICMBRY.EQ.0 ) OPEN (IOBROY,STATUS='SCRATCH',FORM=
     &                            'UNFORMATTED')
C
      ELSE
C   ====================================================================
         WO = 1.0D0 - MIXING
         WN = MIXING
C
         ISIG = 0
         DO IE = 1,NEBROY
            DO IT = 1,NTMAX
               IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                  DO IK1 = 1,NKMMAX
                     DO IK2 = 1,NKMMAX
                        ISIG = ISIG + 1
                        SIGOUT(ISIG) = DREAL(DMFTSIGMA(IK2,IK1,IT,IE))
                        ISIG = ISIG + 1
                        SIGOUT(ISIG) = DIMAG(DMFTSIGMA(IK2,IK1,IT,IE))
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
         DO ISIG = 1,NMAP
            SIGOUT(ISIG) = WO*SIGINP(ISIG) + WN*SIGOUT(ISIG)
         END DO
C
C
C----> broyden updating schemes
         IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY )
     &        CALL SCFBROYPT2(1,SIGINP,SIGOUT,SIGFM,SIGFM1,SIGSM,SIGSM1,
     &        UI,VI,AM,BM,MIXING,ITDEPT,ITDEPTMAX,IMIX,IOBROY,IPF,MIT,
     &        NMAP,NMAPMAX)
C
C
         ISIG = 0
         DO IE = 1,NEBROY
            DO IT = 1,NTMAX
               IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                  DO IK1 = 1,NKMMAX
                     DO IK2 = 1,NKMMAX
                        ISIG = ISIG + 1
                        DMFTSIGMA(IK2,IK1,IT,IE)
     &                     = DCMPLX(SIGOUT(ISIG),SIGOUT(ISIG+1))
                        ISIG = ISIG + 1
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
C -------------------------------------------------------------
C        Calculate root mean quare error of self energy
C -------------------------------------------------------------
C
         SIGO = C0
         SIGI = C0
C
         ISIG = 0
         DO IE = 1,NEBROY
            DO IT = 1,NTMAX
               IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                  IK = 0
                  DO IK1 = 1,NKMMAX
                     DO IK2 = 1,NKMMAX
                        ISIG = ISIG + 1
                        IK = IK + 1
                        SIGO(IK,IT,IE)
     &                     = DCMPLX(SIGOUT(ISIG),SIGOUT(ISIG+1))
                        SIGI(IK,IT,IE)
     &                     = DCMPLX(SIGINP(ISIG),SIGINP(ISIG+1))
                        ISIG = ISIG + 1
                     END DO
                  END DO
               END IF
            END DO
         END DO
         RMSAVSIG = 0.0D0
         DO IT = 1,NTMAX
            IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
C
               IK = 0
               RMSERT = 0.0D0
               RMSERTREL = 0.0D0
               DO IK1 = 1,NKMMAX
                  DO IK2 = 1,NKMMAX
                     IK = IK + 1
                     DO IE = 1,NEBROY
                        CDEL = (SIGI(IK,IT,IE)-SIGO(IK,IT,IE))
                        RMSERT = MAX(RMSERT,ABS(CDEL))
                        IF ( ABS(SIGI(IK,IT,IE)).GT.1D-10 )
     &                       RMSERTREL = MAX(RMSERTREL,ABS(CDEL)
     &                       /ABS(SIGI(IK,IT,IE)))
                     END DO
                  END DO
               END DO
               WRITE (6,99009) IT,RMSERT,RMSERTREL
               RMSAVSIG = MAX(RMSAVSIG,RMSERT)
            END IF
         END DO
         WRITE (6,99007) ITRSCF,RMSAVSIG
C
         IF ( DMFT .AND. ITRSCF.GT.1 ) THEN
C
            IF ( SYSTEM_DIMENSION(1:2).EQ.'2D' .AND. SYSTEM_TYPE(1:3)
     &           .NE.'VIV' ) THEN
               IF ( SUB_SYSTEM(1:6).EQ.'I-ZONE' .OR. SUB_SYSTEM(1:6)
     &              .EQ.'L-BULK' ) THEN
                  IF ( RMSAVSIG.LE.SIGTOL ) DMFT_FIX_DYN_SE = .TRUE.
               END IF
            ELSE IF ( SYSTEM_TYPE(1:16).EQ.'EMBEDDED-CLUSTER' ) THEN
               IF ( RMSAVSIG.LE.SIGTOL ) DMFT_FIX_DYN_SE = .TRUE.
            ELSE
               IF ( RMSAVSIG.LE.SIGTOL .AND. ABS(SHFTEF).LE.SIGTOL )
     &              DMFT_FIX_DYN_SE = .TRUE.
            END IF
C
         END IF
C
C -------------------------------------------------------------
         IF ( SCFDONE .AND. RMSAVSIG.GT.SCFTOL ) THEN
            SCFDONE = .FALSE.
            WRITE (6,99008) 
     &             'DMFT or LDAU self energy not converged  -  continue'
         END IF
C Nur dummy now
         IF ( .FALSE. ) WRITE (6,*) RMSAVV,WETAB
C
C----> reset to start new iteration     pot[old] = pot[new]
C
C
         SIGINP = SIGOUT
         ISIG = 0
         DO IE = 1,NEBROY
            DO IT = 1,NTMAX
               IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                  DO IK1 = 1,NKMMAX
                     DO IK2 = 1,NKMMAX
                        ISIG = ISIG + 1
                        DMFTSIGMA(IK2,IK1,IT,IE)
     &                     = DCMPLX(SIGOUT(ISIG),SIGOUT(ISIG+1))
                        ISIG = ISIG + 1
                     END DO
                  END DO
               END IF
            END DO
         END DO
         IF ( LDAU ) THEN
            DO IE = 2,NEMAX
               DO IT = 1,NTMAX
                  IF ( LOPT(IT).GT.0 .AND. KSELF(IT).EQ.1 ) THEN
                     DO IK1 = 1,NKMMAX
                        DO IK2 = 1,NKMMAX
                           DMFTSIGMA(IK2,IK1,IT,IE)
     &                        = DMFTSIGMA(IK2,IK1,IT,1)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF
C
      END IF
C
99001 FORMAT (/,1X,79('*'),/,28X,
     &        'BROYDEN - mixing scheme - for Self energy')
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
     &        'then jacobian is fixed and sigma  ',
     &        'are updated using that jacobian',/,10X,
     &        'broyden used after ',i3,' iteration different mix')
99007 FORMAT (5X,'ERR iter.',i4,' Self energy average error = ',1p,
     &        d11.4)
99008 FORMAT (5X,A,/,1X,79('*'),/)
99009 FORMAT (5X,'rms-error for type',i3,':  SE = ',1p,d11.4,
     &        ' Relative Err:',1p,d11.4)
C
      END
C
C
