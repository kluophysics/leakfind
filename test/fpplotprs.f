C*==fpplotprs.f    processed by SPAG 6.70Rc at 15:36 on 19 Dec 2016
      SUBROUTINE FPPLOTPRS
C   ********************************************************************
C   *                                                                  *
C   *   write potential, densities and shape functions                 *
C   *   in graphics files                                              *
C   *                                                                  *
C   *   PLOTPRS       ASA                     FULLPOT                  *
C   *    (1):      VT,    BT                + VNST,BNST                *
C   *    (2):      RHOCHR,RHOSPN            + RHO2NS                   *
C   *    (3):                                 SFN                      *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:NONMAG,TUTORIAL
      USE MOD_TYPES,ONLY:NT,NLMFPMAX,LTXT_T,TXT_T,BNST,VNST,BT,VT,
     &    RHO2NS,RHOSPN,RHOORB,RHOCHR,LMIFP,NFPT,IMT
      USE MOD_FILES,ONLY:IOTMP,LSYSTEM,SYSTEM,LDATSET0,DATSET0,PLOTPRS
      USE MOD_RMESH,ONLY:FULLPOT,NSFMAX,NRMAX,LMISF,NSF,FLMSF,NPAN,JRWS,
     &    JRCUT,JRNS1,R,NM
      USE MOD_ANGMOM,ONLY:L_LM,M_LM
      USE MOD_CONSTANTS,ONLY:CONST_4PI
      IMPLICIT NONE
C*--FPPLOTPRS24
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER*80 FILNAM,HEADER1,HEADER2,Y0TXT,Y1TXT
      CHARACTER*11 FMTLM
      INTEGER IFP,IM,IR,IRTOP,ISF,IT,J,JBOT,JOFF,JSF,JTOP,L,LFILNAM,
     &        LHEADER1,LHEADER2,LM,LS,LSTR20,LTXT_IT,LY0TXT,LY1TXT,M,
     &        NGRAPH,NX
      CHARACTER*20 LEG(:),STR20,STR5
      REAL*8 RHORSQ(:),SUMB,XMAX,XMIN,Y(:),Y0MAX,Y0MIN,Y1MAX,Y1MIN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE LEG,Y,RHORSQ
C
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF ( TUTORIAL ) THEN
         SUMB = 0D0
         DO IT = 1,NT
            IM = IMT(IT)
            IF ( FULLPOT ) THEN
               IRTOP = JRCUT(NPAN(IM),IM)
            ELSE
               IRTOP = JRWS(IM)
            END IF
            DO IR = 1,IRTOP
               SUMB = SUMB + ABS(BT(IR,IT))
            END DO
         END DO
         IF ( SUMB/DBLE(NT).LT.0.01D0 ) NONMAG = .TRUE.
      END IF
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
C=======================================================================
C                            POTENTIAL
C=======================================================================
C
      IF ( PLOTPRS(1) ) THEN
C
C ----------------------------------------------------------------------
C                     plot VT, BT
C ----------------------------------------------------------------------
C
         ALLOCATE (Y(NRMAX),RHORSQ(NRMAX))
C
         DO IT = 1,NT
C
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (a.u.)'
            LS = LTXT_T(IT) + 12
            Y0TXT = 'r V!s'//STR20(1:LS)
            LY0TXT = 5 + LS
            IF ( NONMAG ) THEN
               NGRAPH = 1
               Y1TXT = ' '
               LY1TXT = 1
            ELSE
               NGRAPH = 2
               STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
               LS = LTXT_T(IT) + 10
               Y1TXT = 'B!s'//STR20(1:LS)
               LY1TXT = 3 + LS
            END IF
C
            IM = IMT(IT)
            JBOT = 1
            JTOP = JRWS(IM)
C
            NX = JTOP - JBOT + 1
C
            XMIN = R(JBOT,IM)
            XMAX = R(JTOP,IM)
            Y0MAX = 0D0
            Y0MIN = 0D0
            Y1MIN = 0D0
            Y1MAX = 0D0
            DO J = JBOT,JTOP
               Y(J) = VT(J,IT)*R(J,IM)
               Y0MAX = MAX(Y0MAX,Y(J))
               Y0MIN = MIN(Y0MIN,Y(J))
               Y1MAX = MAX(Y1MAX,BT(J,IT))
               Y1MIN = MIN(Y1MIN,BT(J,IT))
            END DO
            IF ( ABS(Y1MAX-Y1MIN).LT.1D-8 ) THEN
               Y1MAX = Y0MAX
               Y1MIN = Y0MIN
            END IF
C
            HEADER2 = 'spherical potential of '//TXT_T(IT)(1:LTXT_T(IT))
            LHEADER2 = 23 + LTXT_T(IT)
            IF ( TUTORIAL ) THEN
               HEADER1 = 'SPR-KKR tutorial'
               LHEADER1 = 16
               LTXT_IT = 0
            ELSE
               HEADER1 = 'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
               LHEADER1 = 25 + LSYSTEM
               LTXT_IT = LTXT_T(IT)
               HEADER2 = HEADER2(1:LHEADER2)//' in '//SYSTEM(1:LSYSTEM)
               LHEADER2 = LHEADER2 + 4 + LSYSTEM
            END IF
C
            CALL XMGRHEAD(DATSET0,LDATSET0,'V',1,TXT_T(IT),LTXT_IT,
     &                    FILNAM,80,LFILNAM,IOTMP,NGRAPH,XMIN,1,XMAX,1,
     &                    Y0MIN,1,Y0MAX,1,Y1MIN,1,Y1MAX,1,
     &                    'radius r (a.u.)',15,Y0TXT,LY0TXT,Y1TXT,
     &                    LY1TXT,HEADER1,LHEADER1,HEADER2,LHEADER2,
     &                    .FALSE.)
C
            CALL XMGRTABLE(0,0,R(JBOT,IM),Y(JBOT),1D0,NX,IOTMP)
            IF ( .NOT.NONMAG ) CALL XMGRTABLE(1,0,R(JBOT,IM),BT(JBOT,IT)
     &           ,1D0,NX,IOTMP)
C
            WRITE (6,99002) 'spherical     potential V,   B    ',
     &                      FILNAM(1:LFILNAM)
            CLOSE (IOTMP)
C
         END DO
         DEALLOCATE (Y)
C
C ----------------------------------------------------------------------
C                     plot VNST, BNST
C ----------------------------------------------------------------------
C
         IF ( FULLPOT ) THEN
C
            ALLOCATE (LEG(NLMFPMAX))
C
            DO IT = 1,NT
C
               STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (Ry)'
C
               LS = LTXT_T(IT) + 10
               Y0TXT = 'V!s'//STR20(1:LS)
               LY0TXT = 3 + LS
               IF ( NONMAG ) THEN
                  NGRAPH = 1
                  Y1TXT = ' '
                  LY1TXT = 1
               ELSE
                  NGRAPH = 2
                  Y1TXT = 'B!s'//STR20(1:LS)
                  LY1TXT = 3 + LS
               END IF
C
               IM = IMT(IT)
               JBOT = JRNS1(IM)
               JTOP = JRCUT(NPAN(IM),IM)
C
               NX = JTOP - JBOT + 1
C
               XMIN = R(JBOT,IM)
               XMAX = R(JTOP,IM)
               Y0MAX = 0D0
               Y0MIN = 0D0
               Y1MIN = 0D0
               Y1MAX = 0D0
               DO J = JBOT,JTOP
                  DO IFP = 1,NFPT(IT)
                     LM = LMIFP(IFP,IT)
C
                     Y0MAX = MAX(Y0MAX,VNST(J,LM,IT))
                     Y0MIN = MIN(Y0MIN,VNST(J,LM,IT))
                     Y1MAX = MAX(Y1MAX,BNST(J,LM,IT))
                     Y1MIN = MIN(Y1MIN,BNST(J,LM,IT))
                  END DO
               END DO
               IF ( ABS(Y1MAX-Y1MIN).LT.1D-8 ) THEN
                  Y1MAX = Y0MAX
                  Y1MIN = Y0MIN
               END IF
C
C
               HEADER2 = 'non-spherical potential of '//TXT_T(IT)
     &                   (1:LTXT_T(IT))
               LHEADER2 = 27 + LTXT_T(IT)
               IF ( TUTORIAL ) THEN
                  HEADER1 = 'SPR-KKR tutorial'
                  LHEADER1 = 16
                  LTXT_IT = 0
               ELSE
                  HEADER1 = 'SPR-KKR calculations for '//
     &                      SYSTEM(1:LSYSTEM)
                  LHEADER1 = 25 + LSYSTEM
                  LTXT_IT = LTXT_T(IT)
                  HEADER2 = HEADER2(1:LHEADER2)
     &                      //' in '//SYSTEM(1:LSYSTEM)
                  LHEADER2 = LHEADER2 + 4 + LSYSTEM
               END IF
C
               CALL XMGRHEAD(DATSET0,LDATSET0,'VNS',3,TXT_T(IT),
     &                       LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,NGRAPH,
     &                       XMIN,1,XMAX,1,Y0MIN,1,Y0MAX,1,Y1MIN,1,
     &                       Y1MAX,1,'radius r (a.u.)',15,Y0TXT,LY0TXT,
     &                       Y1TXT,LY1TXT,HEADER1,LHEADER1,HEADER2,
     &                       LHEADER2,.FALSE.)
C
               DO IFP = 1,NFPT(IT)
                  LM = LMIFP(IFP,IT)
                  WRITE (LEG(IFP),99001) L_LM(LM),M_LM(LM)
               END DO
C
               CALL XMGRLEG1(IOTMP,0,NFPT(IT),LEG,0.18D0,0.48D0)
C
               DO IFP = 1,NFPT(IT)
                  LM = LMIFP(IFP,IT)
C
                  CALL XMGRTABLE(0,IFP-1,R(JBOT,IM),VNST(JBOT,LM,IT),
     &                           1D0,NX,IOTMP)
                  IF ( .NOT.NONMAG )
     &                 CALL XMGRTABLE(1,IFP-1,R(JBOT,IM),BNST(JBOT,LM,
     &                 IT),1D0,NX,IOTMP)
               END DO
C
               WRITE (6,99002) 'NON-spherical potential VNS, BNS  ',
     &                         FILNAM(1:LFILNAM)
               CLOSE (IOTMP)
C
            END DO
C
            DEALLOCATE (LEG)
         END IF
      END IF
C
C
C=======================================================================
C                            DENSITIES
C=======================================================================
C
      IF ( PLOTPRS(2) ) THEN
C
C ----------------------------------------------------------------------
C                     plot RHOCHR,RHOSPN,RHOORB
C ----------------------------------------------------------------------
C
         ALLOCATE (Y(NRMAX))
C
         DO IT = 1,NT
C
            STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (a.u.)'
C
            LS = LTXT_T(IT) + 12
            Y0TXT = 'r!S2!N !xr!0!s'//STR20(1:LS)
            LY0TXT = 14 + LS
            IF ( NONMAG ) THEN
               NGRAPH = 1
               Y1TXT = ' '
               LY1TXT = 1
            ELSE
               NGRAPH = 2
               Y1TXT = 'm!s'//STR20(1:LS)
               LY1TXT = 3 + LS
            END IF
C
            IM = IMT(IT)
            JBOT = 1
            JTOP = JRWS(IM)
C
            NX = JTOP - JBOT + 1
C
            XMIN = R(JBOT,IM)
            XMAX = R(JTOP,IM)
            Y0MAX = 0D0
            Y0MIN = 0D0
            Y1MIN = 0D0
            Y1MAX = 0D0
            DO J = JBOT,JTOP
               Y(J) = VT(J,IT)*R(J,IM)
               RHORSQ(J) = RHOCHR(J,IT)*R(J,IM)*R(J,IM)/CONST_4PI
               Y0MAX = MAX(Y0MAX,RHORSQ(J))
               Y0MIN = MIN(Y0MIN,RHORSQ(J))
               Y1MAX = MAX(Y1MAX,RHOSPN(J,IT))
               Y1MIN = MIN(Y1MIN,RHOSPN(J,IT))
            END DO
            IF ( ABS(Y1MAX-Y1MIN).LT.1D-8 ) THEN
               Y1MAX = Y0MAX
               Y1MIN = Y0MIN
            END IF
C
            IF ( ABS(Y0MAX)+ABS(Y0MIN).GE.1D-8 ) THEN
C
               HEADER2 = 'spherical charge density of '//TXT_T(IT)
     &                   (1:LTXT_T(IT))
               LHEADER2 = 28 + LTXT_T(IT)
               IF ( TUTORIAL ) THEN
                  HEADER1 = 'SPR-KKR tutorial'
                  LHEADER1 = 16
                  LTXT_IT = 0
               ELSE
                  HEADER1 = 'SPR-KKR calculations for '//
     &                      SYSTEM(1:LSYSTEM)
                  LHEADER1 = 25 + LSYSTEM
                  LTXT_IT = LTXT_T(IT)
                  HEADER2 = HEADER2(1:LHEADER2)
     &                      //' in '//SYSTEM(1:LSYSTEM)
                  LHEADER2 = LHEADER2 + 4 + LSYSTEM
               END IF
C
               CALL XMGRHEAD(DATSET0,LDATSET0,'RHO',3,TXT_T(IT),LTXT_IT,
     &                       FILNAM,80,LFILNAM,IOTMP,NGRAPH,XMIN,1,XMAX,
     &                       1,Y0MIN,1,Y0MAX,1,Y1MIN,1,Y1MAX,1,
     &                       'radius r (a.u.)',15,Y0TXT,LY0TXT,Y1TXT,
     &                       LY1TXT,HEADER1,LHEADER1,HEADER2,LHEADER2,
     &                       .FALSE.)
C
               CALL XMGRTABLE(0,0,R(JBOT,IM),RHORSQ,1D0,NX,IOTMP)
               IF ( .NOT.NONMAG ) THEN
                  CALL XMGRTABLE(1,0,R(JBOT,IM),RHOSPN(1,IT),1D0,NX,
     &                           IOTMP)
                  CALL XMGRTABLE(1,0,R(JBOT,IM),RHOORB(1,IT),1D0,NX,
     &                           IOTMP)
               END IF
C
               WRITE (6,99002) 'spherical densities RHOCHR, RHOSPN',
     &                         FILNAM(1:LFILNAM)
            END IF
            CLOSE (IOTMP)
C
         END DO
         DEALLOCATE (Y)
C
C ----------------------------------------------------------------------
C                     plot RHO2NS
C ----------------------------------------------------------------------
C
         IF ( FULLPOT ) THEN
C
            ALLOCATE (LEG(NLMFPMAX))
C
            DO IT = 1,NT
C
               STR20 = TXT_T(IT)(1:LTXT_T(IT))//'!N(r) (a.u.)'
C
               LS = LTXT_T(IT) + 12
               Y0TXT = '!xr!0!s'//STR20(1:LS)
               LY0TXT = 7 + LS
               IF ( NONMAG ) THEN
                  NGRAPH = 1
                  Y1TXT = ' '
                  LY1TXT = 1
               ELSE
                  NGRAPH = 2
                  Y1TXT = 'm!s'//STR20(1:LS)
                  LY1TXT = 3 + LS
               END IF
C
               IM = IMT(IT)
               JBOT = JRNS1(IM)
               JTOP = JRCUT(NPAN(IM),IM)
C
               NX = JTOP - JBOT + 1
C
               XMIN = R(JBOT,IM)
               XMAX = R(JTOP,IM)
               Y0MAX = 0D0
               Y0MIN = 0D0
               Y1MIN = 0D0
               Y1MAX = 0D0
               DO J = JBOT,JTOP
                  DO IFP = 1,NFPT(IT)
                     LM = LMIFP(IFP,IT)
C
                     Y0MAX = MAX(Y0MAX,RHO2NS(J,LM,IT,1))
                     Y0MIN = MIN(Y0MIN,RHO2NS(J,LM,IT,1))
                     Y1MAX = MAX(Y1MAX,RHO2NS(J,LM,IT,2))
                     Y1MIN = MIN(Y1MIN,RHO2NS(J,LM,IT,2))
                  END DO
               END DO
               IF ( ABS(Y1MAX-Y1MIN).LT.1D-8 ) THEN
                  Y1MAX = Y0MAX
                  Y1MIN = Y0MIN
               END IF
C
               IF ( ABS(Y0MAX)+ABS(Y0MIN).GE.1D-8 ) THEN
C
                  HEADER2 = 'non-spherical density of '//TXT_T(IT)
     &                      (1:LTXT_T(IT))
                  LHEADER2 = 25 + LTXT_T(IT)
                  IF ( TUTORIAL ) THEN
                     HEADER1 = 'SPR-KKR tutorial'
                     LHEADER1 = 16
                     LTXT_IT = 0
                  ELSE
                     HEADER1 = 'SPR-KKR calculations for '//
     &                         SYSTEM(1:LSYSTEM)
                     LHEADER1 = 25 + LSYSTEM
                     LTXT_IT = LTXT_T(IT)
                     HEADER2 = HEADER2(1:LHEADER2)
     &                         //' in '//SYSTEM(1:LSYSTEM)
                     LHEADER2 = LHEADER2 + 4 + LSYSTEM
                  END IF
C
                  CALL XMGRHEAD(DATSET0,LDATSET0,'RNS',3,TXT_T(IT),
     &                          LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,
     &                          NGRAPH,XMIN,1,XMAX,1,Y0MIN,1,Y0MAX,1,
     &                          Y1MIN,1,Y1MAX,1,'radius r (a.u.)',15,
     &                          Y0TXT,LY0TXT,Y1TXT,LY1TXT,HEADER1,
     &                          LHEADER1,HEADER2,LHEADER2,.FALSE.)
C
                  DO IFP = 1,NFPT(IT)
                     LM = LMIFP(IFP,IT)
                     WRITE (LEG(IFP),99001) L_LM(LM),M_LM(LM)
                  END DO
C
                  CALL XMGRLEG1(IOTMP,0,NFPT(IT),LEG,0.18D0,0.48D0)
C
                  DO IFP = 1,NFPT(IT)
                     LM = LMIFP(IFP,IT)
C
                     CALL XMGRTABLE(0,IFP-1,R(JBOT,IM),
     &                              RHO2NS(JBOT,LM,IT,1),1D0,NX,IOTMP)
                     IF ( .NOT.NONMAG )
     &                    CALL XMGRTABLE(1,IFP-1,R(JBOT,IM),
     &                    RHO2NS(JBOT,LM,IT,2),1D0,NX,IOTMP)
                  END DO
C
                  WRITE (6,99002) 'NON-spherical density   RHO2NS    ',
     &                            FILNAM(1:LFILNAM)
               END IF
               CLOSE (IOTMP)
C
            END DO
C
            DEALLOCATE (LEG)
         END IF
      END IF
C
C=======================================================================
C                    plot shape function  FLMSF
C=======================================================================
C
      IF ( PLOTPRS(3) .AND. FULLPOT ) THEN
C
         ALLOCATE (LEG(NSFMAX))
C
         DO IM = 1,NM
C
            WRITE (STR5,'(I2)') IM
            Y0TXT = 'f!sLM!N(r)'
            LY0TXT = 10
            NGRAPH = 1
            Y1TXT = ' '
            LY1TXT = 1
C
            JOFF = JRCUT(1,IM)
            JBOT = JRCUT(1,IM) + 1
            JTOP = JRCUT(NPAN(IM),IM)
C
            NX = JTOP - JBOT + 1
C
            XMIN = R(JBOT,IM)
            XMAX = R(JTOP,IM)
            Y0MAX = 0D0
            Y0MIN = 0D0
            Y1MIN = 0D0
            Y1MAX = 0D0
            DO ISF = 1,NSF(IM)
               DO JSF = JBOT - JOFF,JTOP - JOFF
                  Y0MAX = MAX(Y0MAX,FLMSF(JSF,ISF,IM))
                  Y0MIN = MIN(Y0MIN,FLMSF(JSF,ISF,IM))
               END DO
            END DO
C
            STR20 = 'SFN_M'
            CALL STRING_ADD_N(STR20,IM)
            LSTR20 = LEN_TRIM(STR20)
C
            CALL XMGRHEAD(DATSET0,LDATSET0,STR20,LSTR20,' ',0,FILNAM,80,
     &                    LFILNAM,IOTMP,NGRAPH,XMIN,1,XMAX,1,Y0MIN,1,
     &                    Y0MAX,1,Y1MIN,1,Y1MAX,1,'radius r (a.u.)',15,
     &                    Y0TXT,LY0TXT,Y1TXT,LY1TXT,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,'shape functions for mesh '//
     &                    STR5(1:2)//' of '//SYSTEM(1:LSYSTEM),
     &                    (25+2+4+LSYSTEM),.FALSE.)
C
            DO ISF = 1,NSF(IM)
               LM = LMISF(ISF,IM)
               L = L_LM(LM)
               M = M_LM(LM)
               IF ( L.LT.10 ) THEN
                  IF ( M.GE.0 ) THEN
                     FMTLM = '(I1,'','',I1)'
                  ELSE
                     FMTLM = '(I1,'','',I2)'
                  END IF
               ELSE IF ( M.GE.0 ) THEN
                  FMTLM = '(I2,'','',I2)'
               ELSE
                  FMTLM = '(I2,'','',I3)'
               END IF
               WRITE (LEG(ISF),FMT=FMTLM) L_LM(LM),M_LM(LM)
            END DO
C
            CALL XMGRLEG1(IOTMP,0,NSF(IM),LEG,0.65D0,0.80D0)
C
            DO ISF = 1,NSF(IM)
               CALL XMGRTABLE(0,ISF-1,R(JBOT,IM),FLMSF(JBOT-JOFF,ISF,IM)
     &                        ,1D0,NX,IOTMP)
            END DO
C
            WRITE (6,99002) 'shape functions         FLMSF     ',
     &                      FILNAM(1:LFILNAM)
            CLOSE (IOTMP)
C
         END DO
C
         DEALLOCATE (LEG)
C
      END IF
99001 FORMAT (I1,',',I2)
99002 FORMAT (/,10X,A,' written to file ',A)
      END
