C*==check_input.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHECK_INPUT(FULLPOT,POSANIPREP)
C   ********************************************************************
C   *                                                                  *
C   *  check input for compatibility                                   *
C   *                                                                  *
C   *  if the input is not compatible with the TASK                    *
C   *  the program will be stopped                                     *
C   *                                                                  *
C   ********************************************************************
C
      USE MOD_CALCMODE,ONLY:PROGNAME,PUBLIC_VERSION,IREL,DMFT,NONMAG,
     &    KMROT,TASK
      USE MOD_ENERGY,ONLY:IGRID
      IMPLICIT NONE
C*--CHECK_INPUT16
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER*40 ROUTINE
      PARAMETER (ROUTINE='SIG')
C
C Dummy arguments
C
      LOGICAL FULLPOT,POSANIPREP
C
C Local variables
C
      INTEGER IFLAG
C
C*** End of declarations rewritten by SPAG
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      IF ( PUBLIC_VERSION ) THEN
         IF ( DMFT ) CALL STOP_REGULAR(ROUTINE,
     &                        'PUBLIC VERSION does not run in DMFT mode'
     &                        )
         IF ( TASK.EQ.'PHONONS  ' ) CALL STOP_REGULAR(ROUTINE,
     &        'PUBLIC VERSION does not supply phonons')
      END IF
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      IFLAG = 0
C
C-----------------------------------------------------------------------
      IF ( TASK.EQ.'XAS       ' .OR. TASK.EQ.'XES       ' .OR. 
     &     TASK.EQ.'XMO       ' .OR. TASK.EQ.'XRS       ' ) THEN
C
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'AES' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'SCF' ) THEN
C
C         IF ( IGRID(1).NE.5 ) CALL CHKINPWR(TASK,'ONLY E-ARC',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'NRAES' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'APS  ' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'FSCAT' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:6).EQ.'SOCPAR' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:6).EQ.'PSHIFT' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'CLXPS' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'VBPES' .OR. TASK(1:5).EQ.'ARPES' .OR. 
     &          TASK(1:3).EQ.'BIS' ) THEN
C
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:7).EQ.'WFPLOT' ) THEN
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:7).EQ.'MECHECK' ) THEN
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:6).EQ.'MEPLOT' ) THEN
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:2).EQ.'T1' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'EKREL' ) THEN
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'BSF' ) THEN
C
Ccc         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:4).EQ.'FMAG' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
         IF ( IGRID(1).NE.5 ) CALL CHKINPWR(TASK,'ONLY E-ARC',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'JXC' ) THEN
C
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'SP-SREL   ',IFLAG)
         IF ( NONMAG ) CALL CHKINPWR(TASK,'NO PARAMAG',IFLAG)
         IF ( KMROT.NE.0 ) CALL CHKINPWR(TASK,'NO MAG-ROT',IFLAG)
         IF ( IGRID(1).NE.5 ) CALL CHKINPWR(TASK,'ONLY E-ARC',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'OPM' ) THEN
C
        CALL STOP_REGULAR(ROUTINE,'OPM scheme not available')
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:3).EQ.'CHI' ) THEN
C
C         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
Cc         IF ( IGRID(1).NE.5 ) CALL CHKINPWR(TASK,'ONLY E-ARC',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:5).EQ.'SIGMA' ) THEN
C
C         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
C         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
         IF ( IREL.EQ.2 ) CALL CHKINPWR(TASK,'NO IREL=2 ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( POSANIPREP ) THEN
         IF ( FULLPOT ) CALL CHKINPWR(TASK,'NO FULLPOT',IFLAG)
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
C
C-----------------------------------------------------------------------
      ELSE IF ( TASK(1:7).EQ.'COMPTON' .OR. TASK(1:6).EQ.'POSANI' ) THEN
C
         IF ( IREL.NE.3 ) CALL CHKINPWR(TASK,'ONLY SPR  ',IFLAG)
         IF ( IGRID(1).NE.5 ) CALL CHKINPWR(TASK,'ONLY E-ARC',IFLAG)
C
C-----------------------------------------------------------------------
      END IF
C
C-----------------------------------------------------------------------
      IF ( TASK(1:7).EQ.'TUTORIAL' ) THEN
C
         IF ( PROGNAME(1:6).NE.'KKRGEN' )
     &         CALL CHKINPWR(TASK,'USE KKRGEN',IFLAG)
C
      END IF
C-----------------------------------------------------------------------
C
      IF ( IFLAG.EQ.1 ) STOP ' in <CHKINP>'
      END
C*==chkinpwr.f    processed by SPAG 6.70Rc at 15:35 on 19 Dec 2016
      SUBROUTINE CHKINPWR(TASK,KW,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  write a short infor message according to the keyword   KW       *
C   *  informing on the inconsistency and what to do                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C*--CHKINPWR208
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER KMAX
      PARAMETER (KMAX=8)
C
C Dummy arguments
C
      INTEGER IFLAG
      CHARACTER*10 KW,TASK
C
C Local variables
C
      INTEGER K
      CHARACTER*10 L(KMAX)
      CHARACTER*58 T1(KMAX),T2(KMAX),TEXT1,TEXT2
C
C*** End of declarations rewritten by SPAG
C
      DATA L/'NO FULLPOT','ONLY SPR  ','SP-SREL   ','SREL      ',
     &     'NO PARAMAG','NO MAG-ROT','ONLY E-ARC','USE KKRGEN'/
      DATA T1/
     &     'NOT in  FULL POTENTIAL  mode                              ',
     &     'ONLY in  FULLY RELATIVISTIC (SPR)  mode                   ',
     &     'ONLY in  SPIN-POLARIZED SCALAR-RELATIVISTIC  mode         ',
     &     'ONLY in  NON-POLARIZED SCALAR-RELATIVISTIC  mode          ',
     &     'NOT in  PARAMAGNETIC CASE                                 ',
     &     'NOT for   NON-COLLINEAR   or   ROTATED   moments          ',
     &     'ONLY for E-path in complex E-plane                        ',
     &     'within the program  <<KKRGEN>>                            '/
      DATA T2/
     &     'use an  ASA - potential data set for this task            ',
     &     'remove NREL, SREL, or SP-SREL, resp., from  section MODE  ',
     &     'set the switch   SP-SREL   in section   MODE              ',
     &     'set the switch   SREL   in section   MODE                 ',
     &     'use spin-polarised potential or remove NONMAG in CONTROL  ',
     &     'remove all input that may rotate momenst to have  KMROT=0 ',
     &     'set   IGRID=5  and appropriate NE in the section  ENERGY  ',
     &     'call <<KKRGEN>> instead                                   '/
C
      DO K = 1,KMAX
         IF ( KW.EQ.L(K) ) THEN
            TEXT1 = T1(K)
            TEXT2 = T2(K)
            GOTO 100
         END IF
      END DO
      WRITE (6,99002) KW
      STOP ' in <CHKINPWR>'
C
C-----------------------------------------------------------------------
 100  CONTINUE
      IFLAG = 1
C
      WRITE (6,99001) TASK,TEXT1,TEXT2
C
99001 FORMAT (//,4(/,1X,79('#')),//,10X,'input check for task     ',A,
     &        '    failed',//,10X,
     &        'because the necessary subroutine works ',/,10X,A,//,10X,
     &        A,/,4(/,1X,79('#')),//)
99002 FORMAT (//,4(/,1X,79('#')),//,10X,
     &        'STOP in <CHKINPWR> when checking the input ',/,10X,
     &        'the keyword ',A,' is not known ',/,4(/,1X,79('#')),//)
      END
