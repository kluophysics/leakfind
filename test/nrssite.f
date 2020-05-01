C*==nrssite.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NRSSITE(IWRREGWF,IWRIRRWF,IFILSS,ERYD,P,IPRINT,TSST,
     &                   MSST,SSST,MEZZ,MEZJ)
C   ********************************************************************
C   *                                                                  *
C   *    call <RWFBS> to solve the dirac equation without              *
C   *    spin-orbit coupling and set up the radial integrals           *
C   *                                                                  *
C   *  IREL = 0 para-magnetic  non-relativistic     }                  *
C   *                                               } (l,m_l)-repres.  *
C   *  IREL = 1 para-magnetic  scalar-relativistic  }                  *
C   *                                                                  *
C   *  IREL = 2 spin-polarized  scalar-relativistic                    *
C   *           in (l,m_l,m_s)-representation                          *
C   *                                                                  *
C   *  NOTE: the minor component F is seen as an auxilary function     *
C   *        without meaning - it enters only the Wronskian            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_FILES,ONLY:FILNAM,LFILNAM,IOTMP
      USE MOD_RMESH,ONLY:R,NRMAX,JRWS,R2DRDI_W_RADINT,DRDI,NMMAX,NM
      USE MOD_ANGMOM,ONLY:NL,NMEMAX,NKMMAX,NLMAX,NKM
      USE MOD_TYPES,ONLY:CTL,NTMAX,IMT,BT,VT,Z,NT,NLT,LOPT,ITBOT,ITTOP
      USE MOD_CALCMODE,ONLY:IREL,SOLVER,GF_CONV_RH,U_CALCULATION,
     &    U_POT_SHIFT,BREAKPOINT,KMROT
      USE MOD_CONSTANTS,ONLY:C0,CI
      IMPLICIT NONE
C*--NRSSITE29
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NRSSITE')
      LOGICAL WRONSKI
      PARAMETER (WRONSKI=.FALSE.)
C
C Dummy arguments
C
      COMPLEX*16 ERYD,P
      INTEGER IFILSS,IPRINT,IWRIRRWF,IWRREGWF
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 BZJ,BZZ,CHL,CHLP1,CJL,CJL1,CJLP1,CRSQ,CV(:),DPRIRTOP,
     &           DZJ,DZZ,JFS(:),JGS(:,:,:),MSST0,NORM,PI(:),PLFAC,PR(:),
     &           QI(:),QR(:),SSST0,TLP,TSST0,W,WFAC,X,Y,ZFS(:),ZGJG,
     &           ZGS(:,:,:),ZGZG,ZTOP
      REAL*8 C,CSQR,DROVR,GAMMA(:,:),IGAMMA(:,:),RGAM(:,:,:),
     &       RIGAM(:,:,:),RWS_LAST(:),SPNWGT,VSPIN(:)
      COMPLEX*16 CJLZ,CNLZ
      LOGICAL DUMP_SSITE,INITIALIZE,RMESH_CHANGED
      INTEGER I,IA_ERR,IFLAG,IL,IM,IR,IRTOP,IS,IT,IW,J,JRCUT(0:1),JS,L,
     &        LM,LMS,LMSOFF,ML,NLM,NPAN,NSPIN,NSTEP,NT_LAST_CALL,NVIEW
      CHARACTER*10 SOLVER_LAST_CALL
      SAVE GAMMA,IGAMMA,RGAM,RIGAM,RWS_LAST
C
C*** End of declarations rewritten by SPAG
C
      DATA INITIALIZE/.TRUE./
      DATA NT_LAST_CALL/0/,SOLVER_LAST_CALL/'**********'/
      DATA DUMP_SSITE/.FALSE./
C
      ALLOCATABLE RWS_LAST,VSPIN
      ALLOCATABLE RIGAM,RGAM,GAMMA,IGAMMA
      ALLOCATABLE CV,JFS,JGS,PI,PR,QI,QR,ZFS,ZGS
C
      ALLOCATE (PI(NRMAX),PR(NRMAX),QI(NRMAX),QR(NRMAX))
      ALLOCATE (JFS(NRMAX),JGS(NRMAX,NLMAX,2),CV(NRMAX),VSPIN(NRMAX))
      ALLOCATE (ZFS(NRMAX),ZGS(NRMAX,NLMAX,2),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZFS')
C
      C = CTL(1,1)
      CSQR = C*C
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) THEN
         IPRINT = 5
         DUMP_SSITE = .TRUE.
      END IF
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C=======================================================================
C
      IF ( INITIALIZE ) THEN
         ALLOCATE (RWS_LAST(NMMAX))
         RWS_LAST(1:NMMAX) = 0D0
         INITIALIZE = .FALSE.
      END IF
C
      RMESH_CHANGED = .FALSE.
      DO IM = 1,NM
         IF ( ABS(RWS_LAST(IM)-R(JRWS(IM),IM)).GT.1D-12 )
     &        RMESH_CHANGED = .TRUE.
      END DO
C
      IF ( NT_LAST_CALL.LT.NT .OR. SOLVER_LAST_CALL.NE.SOLVER .OR. 
     &     RMESH_CHANGED ) THEN
C
         DO IM = 1,NM
            RWS_LAST(IM) = R(JRWS(IM),IM)
         END DO
C
         NT_LAST_CALL = NT
         SOLVER_LAST_CALL = SOLVER
C
         IF ( ALLOCATED(GAMMA) ) DEALLOCATE (GAMMA)
         IF ( ALLOCATED(IGAMMA) ) DEALLOCATE (IGAMMA)
         IF ( ALLOCATED(RGAM) ) DEALLOCATE (RGAM)
         IF ( ALLOCATED(RIGAM) ) DEALLOCATE (RIGAM)
         ALLOCATE (GAMMA(0:NLMAX,1:NTMAX))
         ALLOCATE (IGAMMA(0:NLMAX,1:NTMAX))
         ALLOCATE (RGAM(1:NRMAX,0:NLMAX,1:NTMAX))
         ALLOCATE (RIGAM(1:NRMAX,0:NLMAX,1:NTMAX),STAT=IA_ERR)
         IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: RIGAM')
C
         DO IT = 1,NT
            IM = IMT(IT)
            IRTOP = JRWS(IM)
            NPAN = 1
            JRCUT(0) = 0
            JRCUT(1) = IRTOP
C
            DO IL = 1,NL
               L = IL - 1
               IF ( IREL.EQ.0 ) THEN
                  GAMMA(L,IT) = L + 1.0D0
                  IGAMMA(L,IT) = -L
               ELSE IF ( Z(IT).EQ.0 ) THEN
                  GAMMA(L,IT) = L + 1.0D0
                  IGAMMA(L,IT) = -L
               ELSE
                  GAMMA(L,IT) = SQRT(DBLE(L*L+L+1)-(2.0D0*Z(IT)/C)**2)
                  IGAMMA(L,IT) = -GAMMA(L,IT)
               END IF
C
               IF ( SOLVER(1:2).EQ.'BS' ) THEN
                  DO IR = 1,IRTOP
                     RGAM(IR,L,IT) = R(IR,IM)**GAMMA(L,IT)
                     RIGAM(IR,L,IT) = R(IR,IM)**IGAMMA(L,IT)
                  END DO
               ELSE IF ( SOLVER.EQ.'RK' ) THEN
                  GAMMA(L,IT) = 0D0
                  IGAMMA(L,IT) = 0D0
                  IF ( ALLOCATED(RGAM) ) DEALLOCATE (RGAM,RIGAM)
               ELSE
C
                  WRITE (6,*) '#### SOLVER not found:  ',SOLVER
C
                  SOLVER = 'BS        '
                  DO IR = 1,IRTOP
                     RGAM(IR,L,IT) = R(IR,IM)**GAMMA(L,IT)
                     RIGAM(IR,L,IT) = R(IR,IM)**IGAMMA(L,IT)
                  END DO
C
                  CALL INFO_MESSAGE(ROUTINE,'SOLVER reset to BS')
C
               END IF
C
            END DO
         END DO
      END IF
C
C=======================================================================
C
      IF ( IREL.EQ.2 ) THEN
         NSPIN = 2
      ELSE
         NSPIN = 1
      END IF
C
C --------------------------------------------------- calculate momentum
C
      CALL GET_MOMENTUM(IREL,C,ERYD,P)
C
      IF ( IREL.EQ.0 ) THEN
         WFAC = C
      ELSE
         WFAC = C*(1.0D0+ERYD/CSQR)
      END IF
C
      NLM = NKM/2
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = ITBOT,ITTOP
C
         MEZZ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
         MEZJ(1:NKMMAX,1:NKMMAX,IT,1:NMEMAX) = C0
C
         IFLAG = 0
         DO IL = 2,NLT(IT)
            IF ( ABS(CTL(IT,IL)-CTL(IT,1)).GT.1D-6 ) IFLAG = 1
         END DO
         IF ( IFLAG.EQ.1 ) THEN
            WRITE (6,*) 'IT = ',IT,'  c =',CTL(IT,1)
            CALL STOP_MESSAGE(ROUTINE,'c inconsistent')
         END IF
C
         DO J = 1,NKMMAX
            DO I = 1,NKMMAX
               TSST(I,J,IT) = C0
               MSST(I,J,IT) = C0
               SSST(I,J,IT) = C0
            END DO
C----------- fill diagonal with 1 to allow matrix inversion for NLT < NL
C                 NOTE: this dummy setting will be overwritten up to NLT
            TSST(J,J,IT) = 1D0
            MSST(J,J,IT) = 1D0
            SSST(J,J,IT) = 1D0
         END DO
C
         IM = IMT(IT)
         IRTOP = JRWS(IM)
C
         NPAN = 1
         JRCUT(0) = 0
         JRCUT(1) = IRTOP
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
         DO IS = 1,NSPIN
C
            IF ( IREL.EQ.2 ) THEN
               SPNWGT = NINT((IS-1.5D0)*2D0)
            ELSE
               SPNWGT = 0D0
            END IF
            LMSOFF = NLM*(IS-1)
C
            DO IR = 1,IRTOP
               VSPIN(IR) = VT(IR,IT) + SPNWGT*BT(IR,IT)
            END DO
C
            LM = 0
C
            ZTOP = P*R(IRTOP,IM)
            PLFAC = P
            DROVR = DRDI(IRTOP,IM)/R(IRTOP,IM)
C
            DPRIRTOP = C0
            PR(IRTOP) = C0
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            DO IL = 1,NLT(IT)
C
               L = IL - 1
C
               PLFAC = PLFAC*DBLE(2*L+1)/P
C
               CV(1:IRTOP) = VSPIN(1:IRTOP)
C
CUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
               IF ( U_CALCULATION .AND. L.EQ.LOPT(IT) ) CV(1:IRTOP)
     &              = CV(1:IRTOP) + U_POT_SHIFT
CUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C
C --------------------------------------- calculate radial wave function
C
               IF ( SOLVER.EQ.'BS' ) THEN
C
                  CALL RWFBS(.TRUE.,C,ERYD,P,L,CV,Z(IT),R(1,IM),
     &                       DRDI(1,IM),IRTOP,JRCUT,NPAN,PR,QR,DPRIRTOP,
     &                       PI,QI,NRMAX,GAMMA(L,IT),IGAMMA(L,IT),IREL,
     &                       GF_CONV_RH)
C
               ELSE IF ( SOLVER.EQ.'RK' ) THEN
C
                  CALL RWFRK(.TRUE.,C,ERYD,P,L,CV,Z(IT),R(1,IM),
     &                       DRDI(1,IM),IRTOP,JRCUT,NPAN,PR,QR,DPRIRTOP,
     &                       PI,QI,NRMAX,GAMMA(L,IT),IGAMMA(L,IT),IREL,
     &                       GF_CONV_RH)
C
               ELSE
                  WRITE (6,*) 'SOLVER = ',SOLVER
                  CALL STOP_MESSAGE(ROUTINE,'SOLVER incorrect')
               END IF
C
C --------------------------------------------------- calculate t-matrix
C
               CJL = CJLZ(L,ZTOP)
               CHL = CJL + CI*CNLZ(L,ZTOP)
C
               CJLP1 = CJLZ(L+1,ZTOP)
               CHLP1 = CJLP1 + CI*CNLZ(L+1,ZTOP)
C
               X = (DBLE(L)/R(IRTOP,IM))*CHL - P*CHLP1
               Y = (DBLE(L)/R(IRTOP,IM))*CJL - P*CJLP1
C
               W = (DPRIRTOP/(PR(IRTOP)*DROVR)+GAMMA(L,IT)-1.0D0)
     &             /R(IRTOP,IM)
C
               TLP = -CI*((CJL*W-Y)/(CHL*W-X))
C
C -------------------------------------------- determine renormalisation
C
C              ALPHAL = (CJL-CI*CHL*TLP)*R(IRTOP,IM)
C    &                  /(PR(IRTOP)*R(IRTOP,IM)**GAMMA(L,IT))
C              SSST0 = ALPHAL*PLFAC/TSST0
C
               TSST0 = TLP/P
               MSST0 = 1.0D0/TSST0
C
C ----------------------------- calculate and store the radial integrals
C
               IF ( SOLVER(1:2).EQ.'BS' ) THEN
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
                  NORM = (CJL/TSST0-CI*P*CHL)
     &                   /(PR(IRTOP)*RGAM(IRTOP,L,IT)/R(IRTOP,IM))
C                                                basis functions R and H
                  IF ( GF_CONV_RH ) NORM = NORM*TSST0
C
                  DO I = 1,IRTOP
                     ZGS(I,IL,IS) = PR(I)*NORM*RGAM(I,L,IT)/R(I,IM)
                     JGS(I,IL,IS) = PI(I)*RIGAM(I,L,IT)/R(I,IM)
                  END DO
C
                  IF ( WRONSKI ) THEN
                     DO I = 1,IRTOP
                        ZFS(I) = QR(I)*NORM*RGAM(I,L,IT)/R(I,IM)/C
                        JFS(I) = QI(I)*RIGAM(I,L,IT)/R(I,IM)/C
                     END DO
                  END IF
C
               ELSE IF ( SOLVER(1:2).EQ.'RK' ) THEN
C
C----------------------------------------- convention for Green function
C                                                basis functions Z and J
                  NORM = (CJL/TSST0-CI*P*CHL)/(PR(IRTOP)/R(IRTOP,IM))
C                                                basis functions R and H
                  IF ( GF_CONV_RH ) NORM = NORM*TSST0
C
                  DO I = 1,IRTOP
                     ZGS(I,IL,IS) = PR(I)*NORM/R(I,IM)
                     JGS(I,IL,IS) = PI(I)/R(I,IM)
                  END DO
C
                  IF ( WRONSKI ) THEN
                     DO I = 1,IRTOP
                        ZFS(I) = QR(I)*NORM/R(I,IM)/C
                        JFS(I) = QI(I)/R(I,IM)/C
                     END DO
                  END IF
C
               END IF
C
C---------------------------------------------- check wronskian relation
C
               IF ( WRONSKI ) THEN
                  WRITE (6,99001) IT,L,IS,ERYD
                  NSTEP = 20
                  NVIEW = 10
                  I = 0
 5                CONTINUE
                  IF ( I.LT.NVIEW .OR. I.GE.(IRTOP-NVIEW) ) THEN
                     I = I + 1
                  ELSE IF ( I.LT.(IRTOP-NVIEW-NSTEP) ) THEN
                     I = I + NSTEP
                  ELSE
                     I = IRTOP - NVIEW
                  END IF
                  IF ( I.LE.IRTOP ) THEN
                     CRSQ = WFAC*R(I,IM)**2
                     WRITE (6,'(3I4,10F18.14)') IT,IL,I,
     &                      (JGS(I,IL,IS)*ZFS(I)-ZGS(I,IL,IS)*JFS(I))
     &                      *CRSQ
                     GOTO 5
                  END IF
               END IF
C
C-----------------------------------------------------------------------
C
C regular part    Z*Z    and    irregular part    Z*J
C
               DZZ = C0
               DZJ = C0
               DO I = 1,IRTOP
                  ZGZG = ZGS(I,IL,IS)*ZGS(I,IL,IS)
                  DZZ = DZZ + ZGZG*R2DRDI_W_RADINT(I,IM)
C
                  ZGJG = ZGS(I,IL,IS)*JGS(I,IL,IS)
                  DZJ = DZJ + ZGJG*R2DRDI_W_RADINT(I,IM)
               END DO
C
               CJL1 = CJLZ(L,P*R(1,IM))
               SSST0 = ZGS(1,IL,IS)/CJL1
C
               IF ( IL.EQ.1 ) THEN
                  CALL NRSSITE_HFF(IT,IL,IS,IM,IRTOP,C,ERYD,ZGS,JGS,BZZ,
     &                             BZJ)
               ELSE
                  BZZ = C0
                  BZJ = C0
               END IF
C
C -------------------------------------------------------- store results
C
               DO ML = -L, + L
                  LM = LM + 1
                  LMS = LMSOFF + LM
C
                  TSST(LMS,LMS,IT) = TSST0
                  MSST(LMS,LMS,IT) = MSST0
                  SSST(LMS,LMS,IT) = SSST0
C
                  MEZZ(LMS,LMS,IT,1) = DZZ
                  MEZJ(LMS,LMS,IT,1) = DZJ
                  MEZZ(LMS,LMS,IT,2) = DZZ*SPNWGT
                  MEZJ(LMS,LMS,IT,2) = DZJ*SPNWGT
                  MEZZ(LMS,LMS,IT,4) = BZZ
                  MEZJ(LMS,LMS,IT,4) = BZJ
               END DO
C
C ------------------------------------------------- write wave functions
C
               IF ( IS.EQ.NSPIN ) THEN
                  IF ( IWRREGWF.NE.0 ) WRITE (IFILSS,REC=IL+(IT-1)*NL)
     &                 IT,L,NL,NSPIN,'NRR',
     &                 ((ZGS(I,IL,JS),I=1,IRTOP),JS=1,NSPIN)
C
                  IF ( IWRIRRWF.NE.0 )
     &                 WRITE (IFILSS,REC=IL+(IT-1+NT)*NL) IT,L,NL,NSPIN,
     &                 'NRI',((JGS(I,IL,JS),I=1,IRTOP),JS=1,NSPIN)
               END IF
C
            END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
         IF ( IPRINT.GE.3 ) THEN
            IW = 6
            WRITE (IW,99003) 'atom type   ',IT,ERYD,CTL(IT,1)
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,IREL,
     &                      IREL,1,1.0D-9,IW)
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IW)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IW)
         END IF
         IF ( DUMP_SSITE ) THEN
            IF ( IT.EQ.1 ) THEN
               LFILNAM = LEN_TRIM(SOLVER)
               FILNAM = 'zzzzzz_ssite_'//SOLVER(1:LFILNAM)
               LFILNAM = LFILNAM + 13
               IF ( GF_CONV_RH ) THEN
                  FILNAM = FILNAM(1:LFILNAM)//'_RH.dat'
               ELSE
                  FILNAM = FILNAM(1:LFILNAM)//'_ZJ.dat'
               END IF
               CALL OPEN_IOTMP_FILE(ROUTINE,IOTMP,FILNAM(1:(LFILNAM+7)))
            END IF
            WRITE (IOTMP,99003) 'atom type   ',IT,ERYD,CTL(IT,1),SOLVER,
     &                          GF_CONV_RH
            CALL CMATSTRUCT('T-MATRIX',TSST(1,1,IT),NKM,NKMMAX,IREL,
     &                      IREL,1,1.0D-9,IOTMP)
            CALL CMATSTRUCT('MEZZ-MATRIX',MEZZ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IOTMP)
            CALL CMATSTRUCT('MEZJ-MATRIX',MEZJ(1,1,IT,1),NKM,NKMMAX,
     &                      IREL,IREL,1,1.0D-9,IOTMP)
            DO I = 1,NKM
               DO J = 1,NKM
                  WRITE (IOTMP,99002) I,J,TSST(I,J,IT),MEZZ(I,J,IT,1),
     &                                MEZJ(I,J,IT,1)
               END DO
            END DO
C
            IF ( IT.EQ.ITTOP ) THEN
               CLOSE (IOTMP)
               WRITE (6,99004) FILNAM(1:(LFILNAM+7))
            END IF
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C         in the case of spin-spiral calculations only
C         rotate ALL matrices from the LOCAL to the GLOBAL frame
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      IF ( KMROT.GE.3 ) CALL SSITE_ROTATE(TSST,MSST,SSST)
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
      IF ( WRONSKI ) CALL STOP_MESSAGE(ROUTINE,'check of Wronskian done'
     &                                 )
C
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
      IF ( BREAKPOINT.EQ.4 ) CALL STOP_BREAKPOINT(ROUTINE)
Cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb BREAKPOINT 4
C
99001 FORMAT (/,' wronski relation for IT =',I3,' L =',I3,' IS =',I3,
     &        ' E =',2F10.5)
99002 FORMAT (2I3,2X,'TSST ',2E25.15,/,8X,'MEZZ ',2E25.15,/,8X,'MEZJ ',
     &        2E25.15,/)
99003 FORMAT (//,1X,79('#'),/,10X,A,I4,/,1X,79('#'),//,10X,'energy',12X,
     &        2F12.6,/,10X,'speed of light',f16.6,/,:,10X,
     &        'solver            ',A,/,10X,'GF_CONV_RH        ',L10,/)
99004 FORMAT (//,5X,'>>>>>  DUMP written to file  ',A,//)
      END
C*==nrssite_hff.f    processed by SPAG 6.70Rc at 15:38 on 19 Dec 2016
      SUBROUTINE NRSSITE_HFF(IT,IL,IS,IM,IRTOP,C,ERYD,ZGS,JGS,BZZ,BZJ)
C   ********************************************************************
C   *                                                                  *
C   *    calculate non-relativistic hyperfine fields                   *
C   *                                                                  *
C   *   KB        1.3807D-16        erg/K                              *
C   *   HBAR      1.05459D-27       erg*s                              *
C   *   RY        13.6058D0         erg                                *
C   *   A0        0.529177D-08      cm                                 *
C   *   MB        9.2741D-21        erg/G                              *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_ANGMOM,ONLY:NLMAX
      USE MOD_RMESH,ONLY:NRMAX,DRDI,R
      USE MOD_CALCMODE,ONLY:IREL
      USE MOD_TYPES,ONLY:Z,VT
      USE MOD_CONSTANTS,ONLY:A0_CGS,MB_CGS,PI,CONST_4PI
      IMPLICIT NONE
C*--NRSSITE_HFF541
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='NRSSITE_HFF')
C
C Dummy arguments
C
      COMPLEX*16 BZJ,BZZ,ERYD
      REAL*8 C
      INTEGER IL,IM,IRTOP,IS,IT
      COMPLEX*16 JGS(NRMAX,NLMAX,2),ZGS(NRMAX,NLMAX,2)
C
C Local variables
C
      COMPLEX*16 CBEXTRA
      COMPLEX*16 DTR,HME,INTR0S_ZJ,INTR0S_ZZ,INTYZJ(:),INTYZZ(:),RPW3,
     &           YZJ(:),YZZ(:),ZJ,ZZ
      INTEGER IR,L
      REAL*8 OV_CSQ,RSQ,RTH
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE INTYZZ,YZZ,INTYZJ,YZJ
C
      ALLOCATE (INTYZZ(NRMAX),YZZ(NRMAX),INTYZJ(NRMAX),YZJ(NRMAX))
C
      IF ( IREL.EQ.0 ) THEN
         OV_CSQ = 0D0
      ELSE IF ( IREL.EQ.1 .OR. IREL.EQ.2 ) THEN
         OV_CSQ = 1D0/C**2
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'IREL > 2')
      END IF
C
C-----------------------------------------------------------------------
C          no hyperfine interaction in case of empty spheres (Z=0)
C-----------------------------------------------------------------------
      IF ( Z(IT).EQ.0 ) THEN
         BZZ = 0D0
         BZJ = 0D0
         RETURN
      END IF
C-----------------------------------------------------------------------
C
      L = IL - 1
C
      IF ( L.EQ.0 ) THEN
C
         IF ( IREL.EQ.0 ) THEN
C
            INTR0S_ZZ = 0.0D0
            INTR0S_ZJ = 0.0D0
            DO IR = 1,IRTOP
               ZZ = ZGS(IR,IL,IS)**2
               ZJ = ZGS(IR,IL,IS)*JGS(IR,IL,IS)
               INTYZZ(IR) = -ZZ/CONST_4PI
               INTYZJ(IR) = -ZJ/CONST_4PI
            END DO
C
         ELSE
C
            RTH = Z(IT)*2/(0.5D0*C**2)/2.0D0
C
            DO IR = 1,IRTOP
               RSQ = R(IR,IM)**2
               DTR = (1.D0/(4D0*PI*R(IR,IM)**2))
     &               *RTH/((1D0+ERYD/C**2)*R(IR,IM)+RTH)**2
               ZZ = ZGS(IR,IL,IS)**2
               ZJ = ZGS(IR,IL,IS)*JGS(IR,IL,IS)
               YZZ(IR) = ZZ*DTR*DRDI(IR,IM)*RSQ
               YZJ(IR) = ZJ*DTR*DRDI(IR,IM)*RSQ
            END DO
C
            CALL CRADINT_R(IM,YZZ,INTYZZ)
            INTR0S_ZZ = INTYZZ(IRTOP)
C
            CALL CRADINT_R(IM,YZJ,INTYZJ)
            INTR0S_ZJ = INTYZJ(IRTOP)
C
         END IF
C
         HME = CBEXTRA(R(1,IM),INTYZZ,INTR0S_ZZ)
         BZZ = (8*PI/3.0D0)*MB_CGS*HME/A0_CGS**3
C
         HME = CBEXTRA(R(1,IM),INTYZJ,INTR0S_ZJ)
         BZJ = (8*PI/3.0D0)*MB_CGS*HME/A0_CGS**3
C
      ELSE
C
         DO IR = 1,IRTOP
C           RPW3 = R(IR,IM)**3*(1D0+(ERYD-VT(IR,IT))*OV_CSQ)
            RPW3 = R(IR,IM)*(1D0+(ERYD-VT(IR,IT))*OV_CSQ)
            ZZ = ZGS(IR,IL,IS)**2
            ZJ = ZGS(IR,IL,IS)*JGS(IR,IL,IS)
            YZZ(IR) = ZZ*DRDI(IR,IM)/RPW3
            YZJ(IR) = ZJ*DRDI(IR,IM)/RPW3
         END DO
C
         CALL CRADINT_R(IM,YZZ,INTYZZ)
         INTR0S_ZZ = INTYZZ(IRTOP)
C
         HME = CBEXTRA(R(1,IM),INTYZZ,INTR0S_ZZ)
         BZZ = 2*MB_CGS*HME/A0_CGS**3
C
         CALL CRADINT_R(IM,YZJ,INTYZJ)
         INTR0S_ZJ = INTYZJ(IRTOP)
C
         HME = CBEXTRA(R(1,IM),INTYZJ,INTR0S_ZJ)
         BZJ = 2*MB_CGS*HME/A0_CGS**3
C
      END IF
C
      END
