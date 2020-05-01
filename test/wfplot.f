C*==wfplot.f    processed by SPAG 6.70Rc at 15:41 on 19 Dec 2016
      SUBROUTINE WFPLOT(TSST,MSST,IPRINT,SSST,MEZZ,MEZJ,NCSTMAX)
C   ********************************************************************
C   *                                                                  *
C   *  subroutine to plot    wave functions                            *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:ORBPOL,NONMAG
      USE MOD_SITES,ONLY:IQAT
      USE MOD_ENERGY,ONLY:ETAB,NETAB
      USE MOD_ANGMOM,ONLY:NKM,NMEMAX,NKMMAX,NLQ,NCPLWF,KAPPA_IKM,L_IKM,
     &    MUEM05_IKM
      USE MOD_TYPES,ONLY:NTMAX,NCPLWFMAX,IKMCPLWF,IMT,LTXT_T,TXT_T,
     &    LCXRAY,NCXRAY
      USE MOD_FILES,ONLY:IOTMP,LRECREAL8,LSYSTEM,SYSTEM,LDATSET,DATSET,
     &    FOUND_SECTION
      USE MOD_RMESH,ONLY:FULLPOT,NRMAX,R,JRWS,RWS,R2DRDI
      USE MOD_CONSTANTS,ONLY:CI,RY_EV
      IMPLICIT NONE
C*--WFPLOT21
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='WFPLOT')
      INTEGER FMATCH
      PARAMETER (FMATCH=4)
C
C Dummy arguments
C
      INTEGER IPRINT,NCSTMAX
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),SSST(NKMMAX,NKMMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX)
C
C Local variables
C
      COMPLEX*16 ARG,CHL,CJL,ERYD,P,ZF(:,:,:),ZG(:,:,:)
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),DR,ECOR(NCSTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),MJ,RINT(NRMAX),
     &       SZCOR(NCSTMAX),XMATCH(FMATCH*NRMAX),XNORM(2),
     &       YMATCH(FMATCH*NRMAX,2),YMAX,YMIN,YY(:,:,:)
      COMPLEX*16 CJLZ,CNLZ
      CHARACTER*2 CL
      CHARACTER*80 FILNAM,YTXT0,YTXT1
      CHARACTER*1 FUNTXTL
      CHARACTER*5 FUNTXTMJ
      LOGICAL GETIRRSOL
      INTEGER IA_ERR,ICST,IFIL,IG,IKM,IKMCOR(NCSTMAX,2),
     &        IKMCPLWFCOR(NCPLWFMAX),IKMWF,IM,IOL,IR,IRTOP,IS,ISEM,IT,
     &        ITSEL,ITXRAY,IZERO(NCSTMAX),J,K,KAP,KAPCOR(NCSTMAX),L,
     &        LFILNAM,LYTXT,MM05,MM05COR(NCSTMAX),N,NC,NCST,NCURVESMAX,
     &        NGRAPH,NKPCOR(NCSTMAX),NMATCH,NSG(0:1),NSOL
      CHARACTER*20 LEG(:,:),STR20
      CHARACTER*4 STATE
      CHARACTER*3 STR3
      CHARACTER*20 SYMBOLKM
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE ZF,ZG,YY,LEG
C
      IF ( FULLPOT ) THEN
         NC = NCPLWFMAX
         CALL STOP_MESSAGE(ROUTINE,'FULLPOT-option not yet available')
      ELSE
         NC = 2
      END IF
      IF ( NC.LE.0 ) CALL STOP_MESSAGE(ROUTINE,'array sizes ???')
C
      ALLOCATE (ZF(NRMAX,NC,NKMMAX),ZG(NRMAX,NC,NKMMAX),STAT=IA_ERR)
      IF ( IA_ERR.NE.0 ) CALL STOP_MESSAGE(ROUTINE,'ALLOC: ZG')
C
      WRITE (6,99004)
C
      GETIRRSOL = .TRUE.
C
C=======================================================================
C
      IFIL = 87
C
      IF ( FULLPOT ) THEN
         IOL = (3+NCPLWFMAX+1+NCPLWFMAX*2*NRMAX*2)*LRECREAL8
      ELSE
         IOL = (5+2*(2+2*NRMAX*2))*LRECREAL8
      END IF
C
      OPEN (UNIT=IFIL,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=IOL)
C
C-----------------------------------------------------------------------
C
      CALL INPUT_FIND_SECTION('TASK',0)
C
      IF ( FOUND_SECTION ) THEN
         ITSEL = 1
         CALL SECTION_SET_INTEGER('ITSEL',ITSEL,9999,1)
         IT = ITSEL
         IM = IMT(IT)
         IRTOP = JRWS(IM)
         WRITE (6,99008) ITSEL
C
         CALL SECTION_SET_INTEGER('PRINT',IPRINT,9999,0)
         CALL SECTION_SET_STRING('STATE',STATE,'BAND',0)
C
      ELSE
         CALL STOP_MESSAGE(ROUTINE,'TASK not found')
      END IF
C
      NSG(0) = 0
      NSG(1) = 0
C
C ======================================================================
C                             CORE state
C ======================================================================
      IF ( STATE.EQ.'CORE' ) THEN
C
         NCURVESMAX = 10
         ALLOCATE (YY(NRMAX,0:NCURVESMAX,0:1),LEG(NCURVESMAX,0:1))
C
         ITXRAY = IT
C
         CALL SECTION_GET_CORE_LEVEL_INFO(CL,NCXRAY(IT),LCXRAY(IT))
C
         CALL SECTION_GET_QUANTUM_NUMBER_MJ(MJ,LCXRAY(IT))
C
         MM05 = NINT(MJ-0.5D0)
C
         NCST = 4*LCXRAY(IT) + 2
C
         STR20 = FUNTXTMJ(MJ)
         WRITE (6,99005) NCXRAY(IT),LCXRAY(IT),STR20(1:5)
C
         IF ( FULLPOT ) THEN
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
            DO ICST = 1,NCST
               IF ( ABS(MM05COR(ICST)+0.5D0).GT.LCXRAY(IT) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
               DO J = 1,NCPLWFMAX
                  IKMCPLWFCOR(J) = 0
               END DO
               DO J = 1,NSOL
                  IKMCPLWFCOR(J) = IKMCOR(ICST,J)
               END DO
               WRITE (IFIL,REC=IKMCOR(ICST,1)+(IT-1)*NKM) IT,'COR',
     &                IKMCOR(ICST,1),JRWS(IM),
     &                (IKMCPLWFCOR(J),J=1,NCPLWFMAX),NSOL,
     &                ((DCMPLX(GCOR(IR,K,ICST),0.0D0),IR=1,JRWS(IM)),
     &                K=1,NSOL),
     &                ((DCMPLX(FCOR(IR,K,ICST),0.0D0),IR=1,JRWS(IM)),
     &                K=1,NSOL)
C
C
            END DO
C
         ELSE
C
            CALL CORE(IPRINT,GCOR,FCOR,ECOR,SZCOR,KAPCOR,MM05COR,NKPCOR,
     &                IKMCOR,IZERO,ITXRAY,BCOR,BCORS,NCSTMAX)
C
            WRITE (6,99001)
C
            DO ICST = 1,NCST
               KAP = KAPCOR(ICST)
               IF ( MM05.EQ.MM05COR(ICST) ) THEN
C
                  DO K = 1,NKPCOR(ICST)
                     DO N = 1,JRWS(IM)
                        RINT(N) = R2DRDI(N,IM)
     &                            *(GCOR(N,K,ICST)**2+FCOR(N,K,ICST)**2)
                     END DO
                     CALL RRADINT(IM,RINT,XNORM(K))
                  END DO
C
                  WRITE (6,99002) ICST,NCXRAY(IT),LCXRAY(IT),
     &                            KAPCOR(ICST),(2*MM05COR(ICST)+1),
     &                            IKMCOR(ICST,1),XNORM(1),ECOR(ICST),
     &                            ECOR(ICST)*RY_EV,SZCOR(ICST),
     &                            IZERO(ICST)
                  IF ( NKPCOR(ICST).EQ.2 ) WRITE (6,99003)
     &                 IKMCOR(ICST,2),XNORM(2)
C
                  IF ( KAPCOR(ICST).LT.0 ) THEN
                     IG = 0
                     IF ( NKPCOR(ICST).GT.1 ) THEN
                        IF ( IKMCOR(ICST,1).LT.IKMCOR(ICST,2) )
     &                       WRITE (6,*) '############ j=l+1/2 ????',
     &                       NKPCOR(ICST)
                        NSG(IG) = 4
                     ELSE
                        NSG(IG) = 2
                     END IF
                  ELSE
                     IF ( IKMCOR(ICST,1).GT.IKMCOR(ICST,2) ) WRITE (6,*)
     &                     '############ j=l-1/2 ????',NKPCOR(ICST)
                     IG = 1
                     NSG(IG) = 4
                  END IF
C
                  DO IR = 1,JRWS(IM)
                     YY(IR,0,IG) = GCOR(IR,1,ICST)
                     YY(IR,1,IG) = FCOR(IR,1,ICST)
                     YY(IR,2,IG) = GCOR(IR,2,ICST)
                     YY(IR,3,IG) = FCOR(IR,2,ICST)
                  END DO
                  STR20 = SYMBOLKM(KAP,MM05)
                  LEG(1,IG) = 'g(r) '//STR20(1:15)
                  LEG(2,IG) = 'f(r) '//STR20(1:15)
                  STR20 = SYMBOLKM(-KAP-1,MM05)
                  LEG(3,IG) = 'g(r) '//STR20(1:15)
                  LEG(4,IG) = 'f(r) '//STR20(1:15)
                  IF ( NONMAG ) THEN
                     ISEM = INDEX(LEG(1,IG),';')
                     IF ( ISEM.GT.1 ) LEG(1:2,IG)(ISEM:20) = ' '
                     LEG(3:4,IG)(1:20) = ' '
                  END IF
C
               END IF
C
            END DO
C
         END IF
C
         IF ( LCXRAY(IT).EQ.0 .OR. ABS(MJ).GT.LCXRAY(IT) ) THEN
            NGRAPH = 1
            IF ( ABS(MJ).GT.LCXRAY(IT) ) YTXT0 = 'g(r)    (j=l+1/2)'
            IF ( LCXRAY(IT).EQ.0 ) YTXT0 = 'g(r)  (l=0;j=1/2)'
            YTXT1 = ' '
         ELSE
            NGRAPH = 2
            YTXT0 = 'g(r)    (j=l+1/2)'
            YTXT1 = 'g(r)    (j=l-1/2)'
         END IF
         LYTXT = 17
C
         YMIN = 0D0
         YMAX = 0D0
         DO IG = 0,(NGRAPH-1)
            DO IS = 0,(NSG(IG)-1)
               DO IR = 1,IRTOP
                  YMIN = MIN(YMIN,YY(IR,IS,IG))
                  YMAX = MAX(YMAX,YY(IR,IS,IG))
               END DO
            END DO
         END DO
C
         CALL XMGRHEAD(DATSET,LDATSET,CL,2,TXT_T(IT),LTXT_T(IT),FILNAM,
     &                 80,LFILNAM,IOTMP,NGRAPH,0.0D0,1,RWS(IM),1,YMIN,0,
     &                 YMAX,1,YMIN,0,YMAX,0,'radius (a!s0!N)',15,YTXT0,
     &                 LYTXT,YTXT1,LYTXT,'SPR-KKR calculations for '//
     &                 SYSTEM(1:LSYSTEM),25+LSYSTEM,
     &                 CL//'-core level wave functions of '//TXT_T(IT)
     &                 (1:LTXT_T(IT)),(32+LTXT_T(IT)),.FALSE.)
C
         CALL XMGRLEGEND(IOTMP,NGRAPH,NSG(0),NSG(1),LEG(1,0),LEG(1,1))
C
         CALL XMGRCURVES(IOTMP,NGRAPH,NSG(0),NSG(1),2,1,0)
C
C
         DO IG = 0,(NGRAPH-1)
            DO IS = 0,(NSG(IG)-1)
               CALL XMGRTABLE(IG,IS,R(1,IM),YY(1,IS,IG),1.0D0,IRTOP,
     &                        IOTMP)
            END DO
         END DO
C
         WRITE (6,99006) 'core',FILNAM(1:LFILNAM)
         CLOSE (IOTMP)
C
C ======================================================================
C                             BAND state
C ======================================================================
C
      ELSE IF ( STATE.EQ.'BAND' ) THEN
C
         NCURVESMAX = MAX(NETAB(1)+1,4)
C
         ALLOCATE (YY(NRMAX,0:NCURVESMAX,0:1),LEG(NCURVESMAX,0:1))
C
         MJ = 0.5D0
         STR20 = FUNTXTMJ(MJ)
C
         ERYD = ETAB(1,1)
C
         CALL SECTION_SET_INTEGER('L',L,NLQ(IQAT(1,IT))-1,0)
C
         CALL SECTION_GET_QUANTUM_NUMBER_MJ(MJ,L)
C
         MM05 = NINT(MJ-0.5D0)
C
         IF ( ABS(ERYD).LT.1D-6 ) ERYD = (0.5D0,0.01D0)
C
         WRITE (6,99007) L,STR20(1:5),ERYD
C
         CALL RUNSSITE(.FALSE.,1,1,IFIL,GETIRRSOL,ERYD,P,IPRINT,TSST,
     &                 MSST,SSST,MEZZ,MEZJ,ORBPOL)
C
         IF ( FULLPOT ) THEN
C
            CALL WAVFUN_READ_REL(IFIL,IT,0,ZG,ZF,ZG,ZF,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C
         ELSE
C
            CALL WAVFUN_READ_REL(IFIL,IT,0,ZG,ZF,ZG,ZF,IRTOP,NCPLWF,
     &                           IKMCPLWF)
C
            DO IKM = 1,NKM
C
               IF ( L.EQ.L_IKM(IKM) .AND. MM05.EQ.MUEM05_IKM(IKM) ) THEN
C
                  IF ( IKM.NE.IKMCPLWF(1,IKM) ) WRITE (6,*)
     &                  '############  IKM ????',IKM,IKMCPLWF(1:2,IKM)
                  NGRAPH = NCPLWF(IKM)
                  IKMWF = IKM
C
                  KAP = KAPPA_IKM(IKMCPLWF(1,IKM))
                  IF ( KAP.LT.0 ) THEN
                     IG = 0
                     IF ( NCPLWF(IKM).GT.1 ) THEN
                        IF ( IKMCPLWF(1,IKM).LT.IKMCPLWF(2,IKM) )
     &                       WRITE (6,*) '############ j=l+1/2 ????',
     &                              IKM,IKMCPLWF(1:2,IKM)
                        NSG(IG) = 4
                     ELSE
                        NSG(IG) = 2
                     END IF
                  ELSE
                     IF ( IKMCPLWF(1,IKM).GT.IKMCPLWF(2,IKM) )
     &                    WRITE (6,*) '############ j=l-1/2 ????',IKM,
     &                                IKMCPLWF(1:2,IKM)
                     IG = 1
                     NSG(IG) = 4
                  END IF
C
                  DO IR = 1,IRTOP
                     YY(IR,0,IG) = DREAL(ZG(IR,1,IKM))
                     YY(IR,1,IG) = DREAL(ZF(IR,1,IKM))
                     YY(IR,2,IG) = DREAL(ZG(IR,2,IKM))
                     YY(IR,3,IG) = DREAL(ZF(IR,2,IKM))
                  END DO
                  STR20 = SYMBOLKM(KAP,MM05)
                  LEG(1,IG) = 'g(r) '//STR20(1:15)
                  LEG(2,IG) = 'f(r) '//STR20(1:15)
                  STR20 = SYMBOLKM(-KAP-1,MM05)
                  LEG(3,IG) = 'g(r) '//STR20(1:15)
                  LEG(4,IG) = 'f(r) '//STR20(1:15)
                  IF ( NONMAG ) THEN
                     ISEM = INDEX(LEG(1,IG),';')
                     IF ( ISEM.GT.1 ) LEG(1:2,IG)(ISEM:20) = ' '
                     LEG(3:4,IG)(1:20) = ' '
                  END IF
C
               END IF
C
            END DO
C
            IF ( NGRAPH.EQ.1 ) THEN
               IF ( ABS(MJ).GT.L ) YTXT0 = 'g(r)    (j=l+1/2)'
               IF ( L.EQ.0 ) YTXT0 = 'g(r)  (l=0;j=1/2)'
               YTXT1 = ' '
            ELSE
               YTXT0 = 'g(r)    (j=l+1/2)'
               YTXT1 = 'g(r)    (j=l-1/2)'
            END IF
            LYTXT = 17
C
            YMIN = 0D0
            YMAX = 0D0
            DO IG = 0,(NGRAPH-1)
               DO IS = 0,(NSG(IG)-1)
                  DO IR = 1,IRTOP
                     YMIN = MIN(YMIN,YY(IR,IS,IG))
                     YMAX = MAX(YMAX,YY(IR,IS,IG))
                  END DO
               END DO
            END DO
C99001 FORMAT (5x,30('#'),'  WARNING ',A,F5.1,A)
C            IF ( ABS(YMIN).GT.ABS(YMAX) ) THEN
C               WRITE (6,99001) ' sign of wave functions flipped'
C               SWAP = YMAX
C               YMAX = -YMIN
C               YMIN = -SWAP
C               DO IG = 0,(NGRAPH-1)
C                  DO IS = 0,(NSG(IG)-1)
C                     DO I = 1,IRTOP
C                        YY(I,IS,IG) = -YY(I,IS,IG)
C                     END DO
C                  END DO
C               END DO
C            END IF
C
            STR3 = FUNTXTL(L)
C
            CALL XMGRHEAD(DATSET,LDATSET,STR3(1:1),1,TXT_T(IT),
     &                    LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,NGRAPH,
     &                    0.0D0,1,RWS(IM),1,YMIN,0,YMAX,1,YMIN,0,YMAX,0,
     &                    'radius (a!s0!N)',15,YTXT0,LYTXT,YTXT1,LYTXT,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,
     &                    'valence band wave functions of '//TXT_T(IT)
     &                    (1:LTXT_T(IT)),(31+LTXT_T(IT)),.FALSE.)
C
            CALL XMGRLEGEND(IOTMP,NGRAPH,NSG(0),NSG(1),LEG(1,0),LEG(1,1)
     &                      )
C
            CALL XMGRCURVES(IOTMP,NGRAPH,NSG(0),NSG(1),2,1,0)
C
C
            DO IG = 0,(NGRAPH-1)
               DO IS = 0,(NSG(IG)-1)
                  CALL XMGRTABLE(IG,IS,R(1,IM),YY(1,IS,IG),1.0D0,IRTOP,
     &                           IOTMP)
               END DO
            END DO
C
            WRITE (6,99006) 'valence band',FILNAM(1:LFILNAM)
            CLOSE (IOTMP)
C
C-------------------------------------------------------------- MATCHING
C                            only for 1 wave function !!!!!!!!
C-------------------------------------------------------------- MATCHING
            NGRAPH = 1
C
            DR = R(IRTOP,IM) - R(IRTOP-1,IM)
            DO IR = 1,IRTOP
               XMATCH(IR) = R(IR,IM)
               YMATCH(IR,1) = YY(IR,0,0)
               ARG = P*XMATCH(IR)
               YMATCH(IR,2) = DREAL(CJLZ(L,ARG))
            END DO
C
            NMATCH = FMATCH*NRMAX
            DO IR = IRTOP + 1,NMATCH
               XMATCH(IR) = XMATCH(IR-1) + DR
               ARG = P*XMATCH(IR)
               CJL = CJLZ(L,ARG)
               CHL = CJL + CI*CNLZ(L,ARG)
               YMATCH(IR,1) = DREAL(CJL/TSST(IKMWF,IKMWF,IT)-CI*P*CHL)
               YMATCH(IR,2) = DREAL(CJL)
               IF ( XMATCH(IR).GT.15D0 ) THEN
                  NMATCH = IR
                  EXIT
               END IF
            END DO
C
            STR3 = FUNTXTL(L)
C
            CALL XMGRHEAD(DATSET,LDATSET,'WF_match_l'//FUNTXTL(L),11,
     &                    TXT_T(IT),LTXT_T(IT),FILNAM,80,LFILNAM,IOTMP,
     &                    NGRAPH,0.0D0,1,XMATCH(NMATCH),1,YMIN,0,YMAX,1,
     &                    YMIN,0,YMAX,0,'radius (a!s0!N)',15,YTXT0,
     &                    LYTXT,YTXT1,LYTXT,
     &                    'SPR-KKR calculations for '//SYSTEM(1:LSYSTEM)
     &                    ,25+LSYSTEM,
     &                    'matching of valence band wave function for '
     &                    //TXT_T(IT)(1:LTXT_T(IT)),(43+LTXT_T(IT)),
     &                    .FALSE.)
C
Cccc        CALL XMGRLEGEND(IOTMP,NGRAPH,NSG(0),NSG(1),LEG(1,0),LEG(1,1))
C
Cccc        CALL XMGRCURVES(IOTMP,NGRAPH,NSG(0),NSG(1),2,1,0)
C
C
            DO IG = 0,(NGRAPH-1)
               DO IS = 0,1
                  CALL XMGRTABLE(IG,IS,XMATCH,YMATCH(1,IS+1),1.0D0,
     &                           NMATCH,IOTMP)
               END DO
            END DO
C
            WRITE (6,99006) 'valence band matching ',FILNAM(1:LFILNAM)
            CLOSE (IOTMP)
C
C-------------------------------------------------------------- MATCHING
C
         END IF
C
      ELSE
         WRITE (6,*) 'STATE = ',STATE
         CALL STOP_MESSAGE(ROUTINE,'value should be  CORE  or  BAND ')
      END IF
C
C ======================================================================
C
      STOP
99001 FORMAT (/,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &        '     E(Ry)       E(eV)    <SIGMA_z>  I0')
99002 FORMAT (5I4,'/2',I4,F12.6,F12.4,F12.3,F12.4,I5)
99003 FORMAT (22X,I4,F12.6)
99004 FORMAT (' ',//,10X,62('*'),/,10X,'*',60X,'*',/,10X,
     &  '*      *     *  ******  *****   *        ****   *******      *'
     &  ,/,10X,
     &  '*      *     *  *       *    *  *       *    *     *         *'
     &  ,/,10X,
     &  '*      *     *  *       *    *  *       *    *     *         *'
     &  ,/,10X,
     &  '*      *  *  *  ***     *****   *       *    *     *         *'
     &  ,/,10X,
     &  '*      * * * *  *       *       *       *    *     *         *'
     &  ,/,10X,
     &  '*      **   **  *       *       *       *    *     *         *'
     &  ,/,10X,
     &  '*      *     *  *       *       ******   ****      *         *'
     &  ,/,10X,'*',60X,'*',/,10X,62('*'),//)
99005 FORMAT (/,10X,'core wave functions to be plotted for',
     &        ' quantum numbers:',//,22X,'n = ',I1,'  l = ',I1,
     &        '  mj = ',A5,/)
99006 FORMAT (/,10X,A,' wave functions written to file: ',A,/)
99007 FORMAT (/,10X,'valence band wave functions to be plotted for',
     &        ' quantum numbers:',//,15X,'l = ',I1,'  mj = ',A5,
     &        ' and energy  E =',F8.3,F6.3,/)
99008 FORMAT (//,10X,'selected atom type   IT =',I3,/)
      END
