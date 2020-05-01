C*==scf_broyden_w_weiss.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE SCF_BROYDEN_W_WEISS(IPRINT,ITRSCF,SCFMIX,OBS_T,
     &                               W_WEISS_ERROR)
C   ********************************************************************
C   *                                                                  *
C   * BROYDEN's iteration scheme                                       *
C   *                                                                  *
C   * applied to mix the type dependent Weiss fields W_WEISS           *
C   *                                                                  *
C   * _________________________________________________________________*
C   *                                                                  *
C   * on entry: itrscf=1 W_T   = WOLD       initialize tables          *
C   *           else     W_T   = WNEW       according to WOLD          *
C   *                                                                  *
C   * on exit:           W_T   = W[mix]                                *
C   * _________________________________________________________________*
C   *                                                                  *
C   * recommended (b.d.):    imix         4                            *
C   *                        istbry       1                            *
C   *                        itdept      40                            *
C   *                        scfmix    5 .. 20 %                       *
C   *                        strmix    == scfmix                       *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NOBSMAX,IWFD
      USE MOD_TYPES,ONLY:ITBOT,ITTOP,NTMAX,NT,W_WEISS_T
      USE MOD_FILES,ONLY:IFILBROY_WEISS
      IMPLICIT NONE
C*--SCF_BROYDEN_W_WEISS30
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SCF_BROYDEN_W_WEISS')
      INTEGER ITDEPTMAX,ISTBRY,IMIX,ITDCMMAX
      PARAMETER (ITDEPTMAX=40,ISTBRY=1,IMIX=4,ITDCMMAX=3)
C
C Dummy arguments
C
      INTEGER IPRINT,ITRSCF
      REAL*8 SCFMIX,W_WEISS_ERROR
      REAL*8 OBS_T(0:3,NOBSMAX,NTMAX)
C
C Local variables
C
      REAL*8 AM(2:ITDEPTMAX-1),BM(2:ITDEPTMAX-1),FM(:),FM1(:),MIXING,
     &       SM(:),SM1(:),STRMIX,UI(:,:),VI(:,:),WN,WO,W_WEISS_ERROR_T,
     &       XINP(:),XOUT(:)
      INTEGER IA_ERR,ICMBRY,IPF,IT,ITDEPT,IVV,MIT,NMAP,NMAPMAX
      SAVE AM,BM,FM,FM1,ICMBRY,IPF,MIT,MIXING,NMAP,NMAPMAX,SM,SM1,
     &     STRMIX,W_WEISS_ERROR_T,XINP,XOUT
      EXTERNAL SCFBROYPT2
C
C*** End of declarations rewritten by SPAG
C
      DATA ITDEPT/40/
C
      ALLOCATABLE UI,VI,FM,SM,FM1,SM1,XINP,XOUT
C
      IF ( ITRSCF.EQ.0 ) THEN
C
         NMAP = ITTOP - ITBOT + 1
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
      ALLOCATE (UI(NMAPMAX,2:ITDCMMAX))
      ALLOCATE (VI(NMAPMAX,2:ITDCMMAX),STAT=IA_ERR)
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
C     store old Weiss field
C
         IVV = 0
         DO IT = ITBOT,ITTOP
            IVV = IVV + 1
            XINP(IVV) = OBS_T(0,IWFD,IT)
         END DO
C
         IF ( IVV.GT.NMAPMAX )
     &        CALL STOP_MESSAGE(ROUTINE,'IVV > NMAPMAX')
C
         IF ( ICMBRY.EQ.0 ) OPEN (IFILBROY_WEISS,STATUS='SCRATCH',
     &                            FORM='UNFORMATTED')
C
      ELSE
C   ====================================================================
C
C---> final construction of the Weiss field
C     first mixing scheme : straight mixing
C---> determination of the root mean square error
C
         W_WEISS_ERROR = 0.0D0
C
         WO = 1.0D0 - MIXING
         WN = MIXING
C
         IVV = 0
         DO IT = ITBOT,ITTOP
            IVV = IVV + 1
C W_WEISS[new]
            XOUT(IVV) = OBS_T(0,IWFD,IT)
            W_WEISS_ERROR_T = ABS(XOUT(IVV)-XINP(IVV))
C W_WEISS[mix]
            WRITE (6,'(A,F10.6)') ' OLD WEISS FIELD:',XINP(IVV)
            WRITE (6,'(A,F10.6)') ' NEW WEISS FIELD:',XOUT(IVV)
            XOUT(IVV) = WO*XINP(IVV) + WN*XOUT(IVV)
C
            IF ( NT.GT.1 ) WRITE (IPF,FMT=99007) IT,W_WEISS_ERROR_T
            W_WEISS_ERROR = DMAX1(W_WEISS_ERROR_T,W_WEISS_ERROR)
         END DO
C
         WRITE (IPF,FMT=99008) ITRSCF,W_WEISS_ERROR
C
C----> Broyden updating schemes
C
         IF ( IMIX.GE.3 .AND. ITRSCF.GE.ISTBRY )
     &        CALL SCFBROYPT2(IPRINT,XINP,XOUT,FM,FM1,SM,SM1,UI,VI,AM,
     &        BM,MIXING,ITDEPT,ITDEPTMAX,IMIX,IFILBROY_WEISS,IPF,MIT,
     &        NMAP,NMAPMAX)
C
C----> reset to start new iteration     W_WEISS[old] = W_WEISS[new]
C
         IVV = 0
         DO IT = ITBOT,ITTOP
            IVV = IVV + 1
C
            XINP(IVV) = XOUT(IVV)
            OBS_T(0,IWFD,IT) = XOUT(IVV)
            W_WEISS_T(IT) = OBS_T(0,IWFD,IT)
            WRITE (6,'(A,F10.6)') ' MIX WEISS FIELD:',W_WEISS_T(IT)
C
         END DO
C
      END IF
C   ====================================================================
C
      DEALLOCATE (UI,VI,STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'DEALLOC')
C
C   ====================================================================
99001 FORMAT (/,1X,79('*'),/,20X,A,' - mixing scheme for SCF - cycle')
99002 FORMAT (/,10X,'ITDEPT to large  array-length ITDEPTMAX=',I5,/,
     &        ' ITDEPT reduced accordingly ',/)
99003 FORMAT (/,10X,'straight mixing factor used:    ',f12.7,/1X,79('*')
     &        ,/)
99004 FORMAT (/,10X,'straight mixing factor used:    ',f12.7,/,10X,
     &        'parameter for Broyden-update:   ',f12.7,/,1X,79('*'),/)
99005 FORMAT (/,10X,'Broyden''s method # ',i1,' used',/,10X,
     &        'iterationdepth to accumulate the Jacobian:  ',i3,/,10X,
     &        'Broyden used after ',i3,' iteration different mix')
99006 FORMAT (/,10X,'Broyden''s method # ',i3,
     &        ' is used up to iteration depth: ',i3,/,10X,
     &        'then Jacobian is fixed and Weiss field',
     &        ' is updated using that Jacobian',/,10X,
     &        'Broyden used after ',i3,' iteration different mix')
99007 FORMAT (5X,'rms-error for type',i3,':  W = ',1p,d11.4)
99008 FORMAT (5X,'iter.',i4,'     average:  W = ',1p,d11.4)
      END
